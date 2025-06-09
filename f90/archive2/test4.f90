program test3
! ======================== TEST SUBROUTINE #3 ==========================
! [DESCRIPTION]
! OPEN SAVED PSI AND CHI, AND THEN PERFORM SOME OPERATIONS. EVENTUALLY,
! GATHER VECTORS/DATA ACROSS DIFFERENT MPI PROCS INTO ONE GLOBAL MATRIX
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
type(scalar):: psi,chi,psio,chio,dpsi,dchi
TYPE(SCALAR):: PSI2,CHI2
TYPE(SCALAR):: PN,CN
REAL(P8):: HDT
complex(P8),dimension(:),allocatable:: vec_1, vec_2, vec_bar
logical:: file_save = .TRUE.

! REFINEMENT
character(LEN=200):: FileName_R
real(p8):: eig0, eig1, eig2, eig00, eig10, eig20
real(p8), allocatable, dimension(:,:):: EMK, EMK0
integer:: converg_flag1, converg_flag2

! DEBUG
TYPE(SCALAR):: RUR,RUP,UZ
TYPE(SCALAR):: ROR,ROP,OZ
TYPE(SCALAR):: W
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: PSI_GLB,CHI_GLB,PSI_LOC, PSI_DIFF
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: RUR_GLB,RUP_GLB,UZ_GLB
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: ROR_GLB,ROP_GLB,OZ_GLB,W_GLB

! =============================== START ================================
CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
   WRITE(*,*) 'The threading support is lesser than that demanded.'
   CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

time_start = mpi_wtime()
CALL READCOM('NOECHO')
CALL READIN(5)
CALL LEGINIT()

status(1,1)=0
if (MPI_RANK.eq.0) then
    WRITE(*,*) 'PROGRAM STARTED'
    call msave(status, 'status.dat')
    CALL PRINT_REAL_TIME()
    ! WRITE(*,*) 'UZ: ', NADD%UZ
    ! WRITE(*,*) 'psi/chi(R) filename: '
    ! READ(*,10) FileName_R
endif
CALL MPI_BCAST(FileName_R,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)

! initilization ========================================================
call allocate(psi)
call allocate(chi)

! call allocate(dpsi)
! call allocate(dchi)
! ======================================================================

DO II = 1,SIZE(PSI%E,1)
   DO JJ = 1,SIZE(PSI%E,2)
      DO KK = 1,SIZE(PSI%E,3)

      PSI%E(II,JJ,KK) = (PSI%INR+II)*10**6 + (PSI%INTH+JJ)*10**3 + KK+PSI%INX

      ENDDO
   ENDDO
ENDDO

ALLOCATE(PSI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))
CALL MASSEMBLE(psi%E,PSI_GLB,1); CALL save_glb(PSI_GLB,"PSI_GLB_FFF")

ALLOCATE(PSI_LOC(SIZE(PSI%E,1),SIZE(PSI%E,2),SIZE(PSI%E,3)))
CALL MDISASSEMBLE(PSI_GLB,PSI_LOC,1);

ALLOCATE(PSI_DIFF(SIZE(PSI%E,1),SIZE(PSI%E,2),SIZE(PSI%E,3)))
PSI_DIFF = ABS(PSI%E - PSI_LOC)

IF (SUM(PSI_DIFF) .EQ. 0.D0) THEN
   WRITE(*,*) MPI_RANK,'PSI_DIFF = ',0

ELSE

   DO II = 1,SIZE(PSI%E,1)
      DO JJ = 1,SIZE(PSI%E,2)
         DO KK = 1,SIZE(PSI%E,3)

         IF (PSI_DIFF(II,JJ,KK).NE.0.D0) WRITE(*,*) MPI_RANK,II,JJ,KK,PSI_DIFF(II,JJ,KK)

         ENDDO
      ENDDO
   ENDDO

ENDIF

DEALLOCATE(PSI_GLB,PSI_LOC)

! ! ======================================================================
! ALLOCATE(PSI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM),CHI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))
! CALL MASSEMBLE(psi%E,PSI_GLB,1); CALL save_glb(PSI_GLB,"PSI_GLB_FFF")
! CALL MASSEMBLE(chi%E,CHI_GLB,1); CALL save_glb(CHI_GLB,"CHI_GLB_FFF")
! DEALLOCATE(PSI_GLB,CHI_GLB)
! ! ======================================================================

! ! ======================================================================
! ALLOCATE(ROR_GLB(NDIMR, NDIMTH, NDIMX),ROP_GLB(NDIMR, NDIMTH, NDIMX),OZ_GLB(NDIMR, NDIMTH, NDIMX))
! CALL MASSEMBLE(ROR%E,ROR_GLB,2); CALL MASSEMBLE(ROP%E,ROP_GLB,2); CALL MASSEMBLE( OZ%E, OZ_GLB,2);
! CALL save_glb(ROR_GLB,"ROR_GLB_PPP"); CALL save_glb(ROP_GLB,"ROP_GLB_PPP"); CALL save_glb( OZ_GLB, "OZ_GLB_PPP") 
! DEALLOCATE(ROR_GLB,ROP_GLB,OZ_GLB)
! ! ======================================================================

CALL DEALLOCATE(psi) 
CALL DEALLOCATE(chi)
! CALL DEALLOCATE(dpsi)
! CALL DEALLOCATE(dchi)

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
! ======================================================================

subroutine save_glb(GLB_DATA,FILENAME,ISERIAL)
! ======================================================================
   COMPLEX(P8),DIMENSION(:,:,:):: GLB_DATA
   character(*):: FILENAME
   INTEGER:: IM,IK

   LOGICAL,OPTIONAL:: ISERIAL

   if ((MPI_RANK.eq.0).OR.(PRESENT(ISERIAL).AND.(ISERIAL))) then
      open(UNIT=777,FILE=ADJUSTL(TRIM(FILENAME))//'.dat',&
      &STATUS='UNKNOWN',ACTION='WRITE')

      DO IM = 1,SIZE(GLB_DATA,2)
         DO IK = 1,SIZE(GLB_DATA,3)
            WRITE(777,256) GLB_DATA(:,IM,IK)
         ENDDO
      ENDDO
      256 FORMAT((S,E24.16E3,SP,E24.16E3,'i'),*(',',S,E24.16E3,SP,E24.16E3,'i'))
      close(777)
   endif

   IF (PRESENT(ISERIAL).AND.(ISERIAL)) RETURN
   call mpi_barrier(MPI_COMM_IVP,IERR)

end subroutine save_glb
! ======================================================================

SUBROUTINE DEL2_TEST(A,B)
!=======================================================================
! [USAGE]: 
! LAPLACIAN (DEL-SQUARE) OPERATOR.
! [PARAMETERS]:
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, DEL^2(A)
! [DEPENDENCIES]:
! 1. BAND_LEG_RAT_DEL2H(~) @ MOD_LIN_LEGENDRE
! 2. BANMUL(~) @ MOD_BANDMAT
! 3. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 4. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
TYPE(SCALAR),INTENT(IN):: A
TYPE(SCALAR),INTENT(INOUT):: B

INTEGER:: MMM,NN,N,KKK
REAL(P8),DIMENSION(:,:),ALLOCATABLE:: DEL2OP
REAL(P8),DIMENSION(:,:,:),ALLOCATABLE:: NORMS, LOGNORMS

IF(A%SPACE.NE. FFF_SPACE) THEN
  IF (MPI_RANK.EQ.0) WRITE(*,*) 'DEL2:INPUT NOT IN FFF_SPACE.'
  STOP
ENDIF

IF(B%SPACE.NE.FFF_SPACE) THEN
  CALL DEALLOCATE(B)
  CALL ALLOCATE(B)
ENDIF

IF ((A%INTH.LT.2).AND.(A%INTH+SIZE(A%E,2).GE.2)) THEN
  MMM = 2 - A%INTH
  IF(ABS(A%E(NRCHOPS(2),MMM,1)).GT.1.E-15) THEN
    WRITE(*,*) 'DEL2: OPERATION NOT EXACT'
    WRITE(*,*)  A%E(NRCHOPS(2),2,1)
    CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)
  ENDIF
ENDIF

ALLOCATE( DEL2OP(NRCHOP,5) )

B%LN = 0.0
CALL CHOPDO(B)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NN,KKK,DEL2OP)
DO MMM=1,SIZE(A%E,2) !NTCHOP
  NN = NRCHOPS(MMM+A%INTH)
  DEL2OP(:NN,:) = BAND_LEG_RAT_DEL2H(NN,M(MMM+A%INTH),ELL,TFM%NORM(:,MMM+A%INTH))

   IF (MMM+A%INTH .EQ. 86) THEN
      WRITE(*,*) "M(MMM+A%INTH): ", M(MMM+A%INTH), "; NN: ", NN, "; ELL: ", ELL
      WRITE(*,*) "TFM%NORM(:,MMM+A%INTH)"
      CALL MCAT(TFM%NORM(:,MMM+A%INTH))
      WRITE(*,*) "BAND_LEG_RAT_DEL2H"
      CALL MCAT(DEL2OP)

      ! SAVE NORM VECTOR AND LOGNORM VECTOR
      IF (A%INX.EQ.0) THEN
         ALLOCATE(NORMS(SIZE(TFM%NORM,1),SIZE(TFM%NORM,2),1), LOGNORMS(SIZE(TFM%LOGNORM,1),SIZE(TFM%LOGNORM,2),1))
         NORMS(:,:,1) = TFM%NORM(:,:)
         LOGNORMS(:,:,1) = TFM%LOGNORM(:,:)
         CALL save_glb(CMPLX(NORMS),'NORMS',.true.)
         CALL save_glb(CMPLX(LOGNORMS),'LOGNORMS',.true.)
         DEALLOCATE(NORMS,LOGNORMS)
      ENDIF

   ENDIF

  DO KKK=1,SIZE(A%E,3) !NXCHOPDIM
    B%E(:NN,MMM,KKK)=BANMUL(DEL2OP(:NN,:),3,A%E(:NN,MMM,KKK))
    B%E(:NN,MMM,KKK)=B%E(:NN,MMM,KKK)-(AK(MMM+A%INTH,KKK+A%INX)**2)*A%E(:NN,MMM,KKK)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

IF ((B%INTH.EQ.0).AND.(B%INX.EQ.0)) THEN
  B%E(1,1,1)=B%E(1,1,1) +4.0D0/3*A%LN/ELL2/TFM%NORM(1,1)
  B%E(2,1,1)=B%E(2,1,1) -2.0D0*A%LN/ELL2/TFM%NORM(2,1)
  B%E(3,1,1)=B%E(3,1,1) +2.0D0/3*A%LN/ELL2/TFM%NORM(3,1)
ENDIF

DEALLOCATE( DEL2OP )

RETURN
END SUBROUTINE DEL2_TEST
!=======================================================================

end program test3
