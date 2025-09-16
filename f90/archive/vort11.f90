program vort11
! ======================================================================
! IVP for the disturbance
! Baseflow is frozen so that:
! Vb x W' + V' x Wb + V' x W'
! ======================================================================
   USE OMP_LIB
   USE MPI
   USE MOD_MISC
   USE MOD_BANDMAT
   USE MOD_EIG
   USE MOD_FD
   USE MOD_LIN_LEGENDRE
   USE MOD_SCALAR3
   USE MOD_FFT 
   USE MOD_LAYOUT
   USE MOD_LEGOPS
   USE MOD_MARCH
   USE MOD_DIAGNOSTICS
   USE MOD_INIT
   USE MOD_EVP, only: SAVE_PERTURB
! -------------------------
implicit none
! -------------------------
integer:: iii,it,mm,kk,pp,II,JJ,FILESTATUS
real(p8),dimension(1,1):: status
real(P8):: time_start,time_end,time0,tmp
type(scalar):: psi,chi,psio,chio,dpsi,dchi,psi_temp,chi_temp
complex(P8),dimension(:),allocatable:: vec_1, vec_2, vec_bar
complex(P8),dimension(:,:,:),allocatable:: psi_glb, chi_glb
logical:: file_save = .TRUE.

! REFINEMENT
character(LEN=200):: FileName_R
real(p8):: eig1, eig2, eig10, eig20
real(p8), allocatable, dimension(:,:):: EMK, EMK0
integer:: converg_flag1, converg_flag2

CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
   WRITE(*,*) 'The threading support is lesser than that demanded.'
   CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

CALL READCOM('NOECHO')
CALL READIN(5)
CALL SETUP_GRID()

status(1,1)=0
if (MPI_RANK.eq.0) then
    WRITE(*,*) 'PROGRAM STARTED'
    call msave(status, 'status.dat')
    CALL PRINT_REAL_TIME()
    WRITE(*,*) 'UZ: ', NADD%UZ
    WRITE(*,*) 'psi/chi(R) filename: '
    READ(*,10) FileName_R
endif
CALL MPI_BCAST(FileName_R,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)

! disturbance
call allocate(psi)
call allocate(chi)
! q-vortex + perturbation
call allocate(psio)
call allocate(chio)
! nonlin result
call allocate(dpsi)
call allocate(dchi)
! temporary file
call allocate(psi_temp)
call allocate(chi_temp)
! energy
allocate(EMK(NTCHOP,NXCHOPDIM))
allocate(EMK0(NTCHOP,NXCHOPDIM))

!> adjust MONITORDATA%MK if necessary
DO II = 1,SIZE(MONITORDATA%MK,1)
   IF (MONITORDATA%MK(II,1).LT.0) MONITORDATA%MK(II,:) = -MONITORDATA%MK(II,:)
   IF (MONITORDATA%MK(II,2).LT.0) MONITORDATA%MK(II,2) = 2*NXCHOP-1+MONITORDATA%MK(II,2)
ENDDO

!> debug
IF (MPI_RANK.EQ.0) THEN
   WRITE(*,*) 'TRACK MODES: (M,AK) - (#MM,#KK)'
   DO II = 1,SIZE(MONITORDATA%MK,1)
      WRITE(*,92) M(MONITORDATA%MK(II,1)+1),AK(MONITORDATA%MK(II,1)+1,MONITORDATA%MK(II,2)+1),MONITORDATA%MK(II,1)+1,MONITORDATA%MK(II,2)+1
   ENDDO
92 FORMAT('(',I2,',',F9.2,') - (#',I2,', #',I2,')')
ENDIF

!> load disturbance
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"1_psi.dat",psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"1_chi.dat",chi)

!> load base flow
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"0_psi.dat",psio)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"0_chi.dat",chio)

!> save the resonant pairs in PFF space
psi_temp%e = psi%e; psi_temp%ln = psi%ln; CALL RTRAN(psi_temp,1)
chi_temp%e = chi%e; chi_temp%ln = chi%ln; CALL RTRAN(chi_temp,1)
DO II = 1,SIZE(PSI%E,2)
   DO JJ = 1,SIZE(PSI%E,3)
      IF (((II+PSI%INTH .EQ. MONITORDATA%MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITORDATA%MK(2,2)+1)) &
         .OR. ((II+PSI%INTH .EQ. MONITORDATA%MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITORDATA%MK(3,2)+1))) THEN
         CALL save_mode(tim%t,psi_temp%E(:,II,JJ),chi_temp%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
      ENDIF
   ENDDO
ENDDO
CALL deallocate(psi_temp); CALL deallocate(chi_temp)

!> first diagnostic
call diagnost(psi,chi)
call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> richardson step
call RICH_LIN(psi,chi,dpsi,dchi,psio,chio)

!> 2nd diagnostic
call diagnost(psi,chi)
call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> save initial energy after Richardson step
EMK0 = ENERGY_SPEC_MODIFIED(psi,chi)

!> startup
!dpsi and dchi are initially empty, then they are assigned 
!the nonlinear part of the first step
time_start = mpi_wtime()
iii = tim%limit/tim%dt
files%n = 1
converg_flag1 = 0
converg_flag2 = 0
eig10 = 0.d0; eig20 = 0.d0

do it=1,iii

   !> admam-bashforth - degen renormalized
   !> psi(n+1),chi(n+1) = factor * |psi(n),chi(n)|_degen + [1.5*|psiN(n),chiN(n)|_degen - 0.5*|psiN(n-1),chiN(n-1)|_degen]
   call ADAMSB_LIN(psi,chi,dpsi,dchi,psio,chio)
   call diagnost(psi,chi)

   !> check convergence in PFF space 
   ! Use real growth rate as the indicator
   EMK = ENERGY_SPEC_MODIFIED(psi,chi)
   eig1 = LOG(EMK(MONITORDATA%MK(2,1)+1,MONITORDATA%MK(2,2)+1)/EMK0(MONITORDATA%MK(2,1)+1,MONITORDATA%MK(2,2)+1))/(2*TIM%DT)
   eig2 = LOG(EMK(MONITORDATA%MK(3,1)+1,MONITORDATA%MK(3,2)+1)/EMK0(MONITORDATA%MK(3,1)+1,MONITORDATA%MK(3,2)+1))/(2*TIM%DT)

   IF (MOD(it,100).EQ.0) THEN

      IF (MPI_RANK.EQ.0) THEN   

         WRITE(*,*) 'STEP: ',it,'/',iii
         WRITE(*,*) '{M1,K1}: ',eig1,' - ',EMK(MONITORDATA%MK(2,1)+1,MONITORDATA%MK(2,2)+1)
         WRITE(*,*) '{M2,K2}: ',eig2,' - ',EMK(MONITORDATA%MK(3,1)+1,MONITORDATA%MK(3,2)+1)

         ! save energy for the two degen modes
         open(UNIT=666,FILE='eig_eneg.dat',STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')
         WRITE(666,134) tim%t,eig1,eig2
         close(666)
         134 FORMAT(F10.3,2(',',SP,E26.15E3))

         ! save energy vs m
         OPEN(UNIT=667,FILE='especData_M.dat',STATUS='UNKNOWN',&
            & IOSTAT=FILESTATUS, POSITION='APPEND')
         WRITE(667, 788) tim%t, (SUM(EMK(MM,:)),MM=1,NTCHOP)          
         CLOSE(667)     

         ! save energy vs k
         OPEN(UNIT=668,FILE='especData_K.dat',STATUS='UNKNOWN',&
            & IOSTAT=FILESTATUS, POSITION='APPEND')
         WRITE(668, 788) tim%t, (SUM(EMK(:,KK)),KK=1,NXCHOPDIM)              
         CLOSE(668)  
         788 FORMAT(*(E23.15))

      ENDIF

   ENDIF
   ! if ((abs(eig1-eig10)/abs(eig10).lt.1e-8) .AND. (abs(eig2-eig20)/abs(eig20).lt.1e-8)) GOTO 999
   eig10 = eig1; eig20 = eig2

   !> save energy
   EMK0 = EMK

   !> save log
   IF (MOD(it,1000).EQ.0) THEN
      call allocate(psi_temp); psi_temp%e = psi%e; CALL RTRAN(psi_temp,1)
      call allocate(chi_temp); chi_temp%e = chi%e; CALL RTRAN(chi_temp,1)
      DO II = 1,SIZE(PSI%E,2)
         DO JJ = 1,SIZE(PSI%E,3)
            IF (((II+PSI%INTH .EQ. MONITORDATA%MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITORDATA%MK(2,2)+1)) &
               .OR. ((II+PSI%INTH .EQ. MONITORDATA%MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITORDATA%MK(3,2)+1))) THEN
               CALL save_mode(tim%t,psi_temp%E(:,II,JJ),chi_temp%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
            ENDIF
         ENDDO
      ENDDO
      CALL DEALLOCATE(psi_temp); CALL DEALLOCATE(chi_temp)
   ENDIF
   ! IF (it.eq.1001) CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)

   !> output
   if ((files%t(files%n).le.tim%t) .AND. (file_save)) then

      call msave(psi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi(files%n))
      call msave(chi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi(files%n))
      files%n = files%n + 1

      ! stop the simulation once the last file is saved
      ! if(files%n > files%ne) goto 999
      if(files%n > files%ne) file_save = .FALSE.
      
   endif

enddo

999 continue
time_end = mpi_wtime()

!> output
!> 1. FFF space: perturbation file
ALLOCATE(PSI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM),CHI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))
CALL MASSEMBLE(PSI%E,PSI_GLB,1); CALL MASSEMBLE(CHI%E,CHI_GLB,1)
IF (MPI_RANK.EQ.0) THEN
   ALLOCATE(vec_1(2*NRCHOPDIM),vec_2(2*NRCHOPDIM),vec_bar(2*NRCHOPDIM))
   vec_1(1:NRCHOPDIM) = PSI_GLB(:,MONITORDATA%MK(2,1)+1,MONITORDATA%MK(2,2)+1)
   vec_1(NRCHOPDIM+1:) = CHI_GLB(:,MONITORDATA%MK(2,1)+1,MONITORDATA%MK(2,2)+1)
   vec_2(1:NRCHOPDIM) = PSI_GLB(:,MONITORDATA%MK(3,1)+1,MONITORDATA%MK(3,2)+1)
   vec_2(NRCHOPDIM+1:) = CHI_GLB(:,MONITORDATA%MK(3,1)+1,MONITORDATA%MK(3,2)+1) 
   vec_bar(1:NRCHOPDIM) = PSI_GLB(:,MONITORDATA%MK(1,1)+1,MONITORDATA%MK(1,2)+1)
   vec_bar(NRCHOPDIM+1:) = CHI_GLB(:,MONITORDATA%MK(1,1)+1,MONITORDATA%MK(1,2)+1)
   CALL SAVE_PERTURB('./converg/new_perturb_final.input', &
      M(MONITORDATA%MK(2,1)), AK(MONITORDATA%MK(2,1),MONITORDATA%MK(2,2)),vec_1, CMPLX(eig1), &
      M(MONITORDATA%MK(3,1)), AK(MONITORDATA%MK(3,1),MONITORDATA%MK(3,2)),vec_2, CMPLX(eig2), &
      M(MONITORDATA%MK(1,1)), AK(MONITORDATA%MK(1,1),MONITORDATA%MK(1,2)),vec_bar, CMPLX(0.D0))
   DEALLOCATE(vec_1,vec_2,vec_bar)
ENDIF
DEALLOCATE(PSI_GLB,CHI_GLB)

!> 2. PFF space: PSI-CHI
CALL RTRAN(psi,1); CALL RTRAN(chi,1)
DO II = 1,SIZE(PSI%E,2)
   DO JJ = 1,SIZE(PSI%E,3)
      IF (((II+PSI%INTH .EQ. MONITORDATA%MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITORDATA%MK(2,2)+1)) &
         .OR. ((II+PSI%INTH .EQ. MONITORDATA%MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITORDATA%MK(3,2)+1))) THEN
         CALL save_mode(tim%t,psi%E(:,II,JJ),chi%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
      ENDIF
   ENDDO
ENDDO
CALL RTRAN(psi,-1); CALL RTRAN(chi,-1)

!> deallocate
call deallocate(psi)
call deallocate(chi)
call deallocate(psio)
call deallocate(chio)

!> final printout
IF (MPI_RANK.eq.0) THEN
    print *,tim%n,' steps'
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()
    WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

! ======================================================================
contains
! ======================================================================

subroutine save_mode(t,psi_new,chi_new,m_i,k_i)
! ======================================================================
   complex(p8),DIMENSION(:):: psi_new,chi_new
   real:: t
   integer:: m_i,k_i
   integer:: nn

   open(UNIT=777,FILE='mode_track_'//ITOA3(m_i)//'_'//ITOA3(k_i)//'.dat',&
   &STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')

   WRITE(777,320) t,psi_new(1:size(psi_new)),chi_new(1:size(psi_new))

   close(777)
320 FORMAT(F10.3,',',(S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))
end subroutine save_mode
! ======================================================================

SUBROUTINE NONLIN_LIN(PSI1,CHI1,PSIN,CHIN,PSIB,CHIB)
!=======================================================================
! [USAGE]: 
! GIVEN POLOIDAL-TOROIDAL TERMS PSI AND CHI, COMPUTE THE NON-LINEAR
! POLOIDAL-TOROIDAL TERMS PSIN AND CHIN FROM (U X OMEGA) (VPROD).
! SEE EQ. (101) IN MATSUSHIMA AND MARCUS (1997).
! [PARAMETERS]:
! PSI >> TOROIDAL TERM OF THE VELOCITY FIELD
! CHI >> POLOIDAL TERM OF THE VELOCITY FIELD
! PSIN >> TOROIDAL TERM OF THE (U X OMEGA) FIELD
! CHIN >> POLOIDAL TERM OF THE (U X OMEGA) FIELD
! [DEPENDENCIES]:
! 1. OPERATOR(.MUL.) @ MOD_EIG
! [UPDATES]:
! CODED BY JINGE WANG @ MAR 29 2023
!=======================================================================
   IMPLICIT NONE
   TYPE(SCALAR),INTENT(IN):: PSI1,CHI1,PSIB,CHIB
   TYPE(SCALAR),INTENT(INOUT)::PSIN,CHIN

   TYPE(SCALAR):: RUR1,RUP1,UZ1
   TYPE(SCALAR):: ROR1,ROP1,OZ1
   TYPE(SCALAR):: ROR ,ROP ,OZ
   TYPE(SCALAR):: W

   CALL ALLOCATE(ROR,PPP_SPACE); CALL ALLOCATE(ROP,PPP_SPACE); CALL ALLOCATE(OZ,PPP_SPACE)
   ROR%E = 0.D0; ROP%E = 0.D0; OZ%E = 0.D0

   CALL ALLOCATE(RUR1); CALL ALLOCATE(RUP1); CALL ALLOCATE(UZ1)
   CALL ALLOCATE(ROR1); CALL ALLOCATE(ROP1); CALL ALLOCATE(OZ1)

   CALL CHOPSET(3)

   !> Vb x W'
   CALL PC2VEL(PSIB,CHIB,RUR1,RUP1,UZ1)                                    ! NOTE: INPUT AND OUTPUTS ARE STILL IN FFF SPACE
   CALL PC2VOR(PSI1,CHI1,ROR1,ROP1,OZ1)                                    ! NOTE: INPUT AND OUTPUTS ARE STILL IN FFF SPACE
   CALL VPROD(RUR1,RUP1,UZ1,ROR1,ROP1,OZ1)                                 ! IN PPP SPACE
   ROR%E = ROR%E + ROR1%E
   ROP%E = ROP%E + ROP1%E
    OZ%E =  OZ%E +  OZ1%E
   
   !> V' x Wb
   CALL PC2VEL(PSI1,CHI1,RUR1,RUP1,UZ1)
   CALL PC2VOR(PSIB,CHIB,ROR1,ROP1,OZ1)
   CALL VPROD(RUR1,RUP1,UZ1,ROR1,ROP1,OZ1,IS_FREESTREAM=.FALSE.)
   ROR%E = ROR%E + ROR1%E
   ROP%E = ROP%E + ROP1%E
    OZ%E =  OZ%E +  OZ1%E

   ! !> V' x W'
   ! CALL PC2VEL(PSI1,CHI1,RUR1,RUP1,UZ1)
   ! CALL PC2VOR(PSI1,CHI1,ROR1,ROP1,OZ1)
   ! CALL VPROD(RUR1,RUP1,UZ1,ROR1,ROP1,OZ1,IS_FREESTREAM=.FALSE.)
   ! ROR%E = ROR%E + ROR1%E
   ! ROP%E = ROP%E + ROP1%E
   !  OZ%E =  OZ%E +  OZ1%E

   CALL CHOPSET(-3)

   !> CALCULATING NONLINEAR TERMS PSIN, CHIN
   CALL ALLOCATE(W)
   CALL PROJECT(ROR,ROP,OZ,PSIN,W)                                    ! NOW PSIN = PSI_NONLINEAR AND W = DELSQ(CHI_NONLINEAR)
   CALL IDEL2(W,CHIN)                                                 ! NOW CHIN = CHI_NONLINEAR
   CALL DEALLOCATE(W)

   CALL DEALLOCATE(RUR1); CALL DEALLOCATE(RUP1); CALL DEALLOCATE(UZ1)
   CALL DEALLOCATE(ROR1); CALL DEALLOCATE(ROP1); CALL DEALLOCATE(OZ1)
   CALL DEALLOCATE(ROR ); CALL DEALLOCATE(ROP ); CALL DEALLOCATE(OZ )

   RETURN
END SUBROUTINE NONLIN_LIN
!=======================================================================

SUBROUTINE ADAMSB_LIN(PSI1,CHI1,PSINO,CHINO,PSIB,CHIB)
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR),INTENT(IN):: PSIB,CHIB
   TYPE(SCALAR):: PSI1,CHI1,PSINO,CHINO
   TYPE(SCALAR):: PSIN,CHIN

   TYPE(SCALAR):: PSI0,CHI0
   REAL(P8):: DT
   INTEGER:: I1, J1

   DT = TIM%DT
   TIM%T = TIM%T + DT
   TIM%N = TIM%N + 1
   ADV%X = ADV%X + ADV%UX*DT
   ADV%Y = ADV%Y + ADV%UY*DT

   CALL ALLOCATE( PSIN )
   CALL ALLOCATE( CHIN )

   !> calculate nonlinear term
   CALL NONLIN_LIN(PSI1,CHI1,PSIN,CHIN,PSIB,CHIB)

!    CALL ALLOCATE( PSI0 )
!    CALL ALLOCATE( CHI0 )

!    PSI0 = PSI1
!    CHI0 = CHI1

   !> update psi & chi
   PSI1%E = PSI1%E +DT*(1.5D0*PSIN%E -0.5D0*PSINO%E)
   CHI1%E = CHI1%E +DT*(1.5D0*CHIN%E -0.5D0*CHINO%E)

   PSINO=PSIN
   CHINO=CHIN

   CALL DEALLOCATE( PSIN )
   CALL DEALLOCATE( CHIN )

   CALL CALC_INCOMP_HYPERV(PSI1,CHI1,DT)
   CALL HYPERV3(PSI1,CHI1)

!    CALL DEALLOCATE(PSI0)
!    CALL DEALLOCATE(CHI0)

   RETURN
END SUBROUTINE
! ======================================================================

SUBROUTINE RICH_LIN(PSI1,CHI1,PSIN,CHIN,PSIB,CHIB)
!=======================================================================
   IMPLICIT NONE
   TYPE(SCALAR),INTENT(IN):: PSIB,CHIB
   TYPE(SCALAR):: PSI1,CHI1,PSIN,CHIN

   TYPE(SCALAR):: PSI2,CHI2,PSI0,CHI0
   TYPE(SCALAR):: PN,CN
   REAL(P8):: DT, HDT
   INTEGER:: I1, J1

   DT = TIM%DT
   HDT = DT*0.5D0
   TIM%T = TIM%T + DT
   TIM%N = TIM%N + 1
   ADV%X = ADV%X + ADV%UX*DT
   ADV%Y = ADV%Y + ADV%UY*DT

   CALL ALLOCATE(PSI2)
   CALL ALLOCATE(CHI2)

   CALL ALLOCATE( PN )
   CALL ALLOCATE( CN )

   ! FIRST HALF-STEP
   CALL NONLIN_LIN(PSI1,CHI1,PSIN,CHIN,PSIB,CHIB)
   PSI2%E =PSI1%E + (HDT*PSIN%E)
   CHI2%E =CHI1%E + (HDT*CHIN%E)
   CALL CALC_INCOMP_HYPERV(PSI2,CHI2,HDT)
   CALL HYPERV3(PSI2,CHI2)

   ! SECOND HALF-STEP
   CALL NONLIN_LIN(PSI2,CHI2,PN,CN,PSIB,CHIB)

   PSI2%E =PSI2%E + (HDT*PN%E)
   CHI2%E =CHI2%E + (HDT*CN%E)
   CALL CALC_INCOMP_HYPERV(PSI2,CHI2,HDT)
   CALL HYPERV3(PSI2,CHI2)

   ! FULL STEP
   PSI1%E =PSI1%E +(DT)*PSIN%E
   CHI1%E =CHI1%E +(DT)*CHIN%E

   CALL CALC_INCOMP_HYPERV(PSI1,CHI1,DT)
   CALL HYPERV3(PSI1,CHI1)

   PSI1%E =(2.0D0*PSI2%E)-PSI1%E
   CHI1%E =(2.0D0*CHI2%E)-CHI1%E

   CALL DEALLOCATE( PSI2 )
   CALL DEALLOCATE( CHI2 )
   CALL DEALLOCATE( PN )
   CALL DEALLOCATE( CN )

   RETURN
END SUBROUTINE
! ======================================================================

end program vort11
