program vort9
! ======================================================================
! OPEN RE-FINED EIGENVECTORS AND RUN NONLINEAR SIMULATION 
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
type(scalar):: psi,chi,psio,chio,dpsi,dchi,dpsi2,dchi2
complex(P8),dimension(:),allocatable:: vec_1, vec_2, vec_bar
complex(P8),dimension(:,:,:),allocatable:: psi_glb, chi_glb
logical:: file_save = .TRUE.

! REFINEMENT
character(LEN=200):: FileName_R
real(p8):: eig0, eig1, eig2, eig00, eig10, eig20
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
CALL LEGINIT()

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

! degen triad + Q-vortex
call allocate(psi)
call allocate(chi)
! initial value
call allocate(psio)
call allocate(chio)
! nonlin result
call allocate(dpsi)
call allocate(dchi)
! energy
allocate(EMK(NTCHOP,NXCHOPDIM))
allocate(EMK0(NTCHOP,NXCHOPDIM))

!> load Q-vortex + resonant triad
! call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_psi.dat",psi)
! call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_chi.dat",chi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//"psi0.dat",psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//"chi0.dat",chi)

call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> adjust monitor_mk if necessary
DO II = 1,SIZE(MONITOR_MK,1)
   IF (MONITOR_MK(II,1).LT.0) MONITOR_MK(II,:) = -MONITOR_MK(II,:)
   IF (MONITOR_MK(II,2).LT.0) MONITOR_MK(II,2) = 2*NXCHOP-1+MONITOR_MK(II,2)
ENDDO

!> debug
IF (MPI_RANK.EQ.0) THEN
   WRITE(*,*) 'TRACK MODES: (M,AK) - (#MM,#KK)'
   DO II = 1,SIZE(MONITOR_MK,1)
      WRITE(*,92) M(MONITOR_MK(II,1)+1),AK(MONITOR_MK(II,1)+1,MONITOR_MK(II,2)+1),MONITOR_MK(II,1)+1,MONITOR_MK(II,2)+1
   ENDDO
92 FORMAT('(',I2,',',F9.2,') - (#',I2,', #',I2,')')
ENDIF

! !> save the resonant pairs in PFF space
! call allocate(psio); psio%e = psi%e; psio%ln = psi%ln; CALL RTRAN(psio,1)
! call allocate(chio); chio%e = chi%e; chio%ln = chi%ln; CALL RTRAN(chio,1)
! DO II = 1,SIZE(PSI%E,2)
!    DO JJ = 1,SIZE(PSI%E,3)
!       IF (((II+PSI%INTH .EQ. MONITOR_MK(1,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(1,2)+1)) &
!          .OR. ((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
!          .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
!          CALL save_mode(tim%t,psio%E(:,II,JJ),chio%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
!       ENDIF
!    ENDDO
! ENDDO
! CALL DEALLOCATE(psio); CALL DEALLOCATE(chio)

!> first diagnostic
call diagnost(psi,chi)
call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> richardson step
call rich(psi,chi,dpsi,dchi)

!> 2nd diagnostic
call diagnost(psi,chi)
call PRINT_ENERGY_SPECTRUM(psi,chi,1)

call DEALLOCATE(psi)
call DEALLOCATE(chi)

! CALL ALLOCATE(PSIO); CALL ALLOCATE(CHIO)

! ! FIRST HALF-STEP
! TIM%T = 5
! TIM%T = TIM%T + 1
! CALL NONLIN(psi,chi,dpsi,dchi)
! PSIO%E =PSI%E + (0.5D0*TIM%DT*dpsi%E)
! CHIO%E =CHI%E + (0.5D0*TIM%DT*dchi%E)
! CALL VISC1(PSIO,CHIO,0.5D0*TIM%DT)
! CALL HYPERV(PSIO,CHIO,0.5D0*TIM%DT)

! call diagnost(PSIO,CHIO)
! call PRINT_ENERGY_SPECTRUM(PSIO,CHIO,1)

! ! SECOND HALF-STEP
! TIM%T = 5
! TIM%T = TIM%T + 1
! CALL ALLOCATE(dpsi2); CALL ALLOCATE(dchi2)
! CALL NONLIN(PSIO,CHIO,dpsi2,dchi2)
! PSIO%E =PSIO%E + (0.5D0*TIM%DT*dpsi2%E)
! CHIO%E =CHIO%E + (0.5D0*TIM%DT*dchi2%E)
! CALL VISC1(PSIO,CHIO,0.5D0*TIM%DT)
! CALL HYPERV(PSIO,CHIO,0.5D0*TIM%DT)

! call diagnost(PSIO,CHIO)
! call PRINT_ENERGY_SPECTRUM(PSIO,CHIO,1)

! ! FULL STEP
! TIM%T = 5
! TIM%T = TIM%T + 1
! PSI%E =PSI%E +(TIM%DT)*dpsi%E
! CHI%E =CHI%E +(TIM%DT)*dchi%E
! CALL VISC1(PSI,CHI,TIM%DT)
! CALL HYPERV(PSI,CHI,TIM%DT)

! call diagnost(PSI,CHI)
! call PRINT_ENERGY_SPECTRUM(PSI,CHI,1)

! ! FINAL STEP
! TIM%T = 5
! TIM%T = TIM%T + 1
! PSI%E =(2.0D0*PSIO%E)-PSI%E
! CHI%E =(2.0D0*CHIO%E)-CHI%E

! call diagnost(PSI,CHI)
! call PRINT_ENERGY_SPECTRUM(PSI,CHI,1)

! !> 2nd diagnostic
! CALL DEALLOCATE(dpsi2)
! CALL DEALLOCATE(dchi2)
! CALL DEALLOCATE(PSIO)
! CALL DEALLOCATE(CHIO)
! CALL DEALLOCATE(PSI)
! CALL DEALLOCATE(CHI)


!> save initial energy after Richardson step
! EMK0 = ENERGY_SPEC_MODIFIED(psi,chi)

999 continue
time_end = mpi_wtime()

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

end program vort9
