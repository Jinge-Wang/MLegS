program test_assemble
! ======================================================================
! TEST MASSABLE COMMAND
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
integer:: iii,it,mm,pp,II,JJ,KK,FILESTATUS
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
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: PSI_GLB,CHI_GLB
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: RUR_GLB,RUP_GLB,UZ_GLB
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: ROR_GLB,ROP_GLB,OZ_GLB,W_GLB

! =============================== START ================================
CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
   WRITE(*,*) 'The threading support is lesser than that demanded.'
   CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

CALL READCOM('NOECHO')
CALL READIN(5)
CALL LEGINIT()

if (MPI_RANK.eq.0) then
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()
endif

! ======================================================================
call allocate(psi,FFF_SPACE)
! call allocate(chi,FFF_SPACE)
ALLOCATE(PSI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))

DO II=1,SIZE(PSI%E,1)
  DO JJ=1,SIZE(PSI%E,2)
    DO KK=1,SIZE(PSI%E,3)
        PSI%E(II,JJ,KK) = (II+PSI%INR)*10**6+(JJ+PSI%INTH)*10**3+(KK+PSI%INX)
    ENDDO
  ENDDO
ENDDO

CALL MASSEMBLE(PSI%E,PSI_GLB,1); CALL save_glb(PSI_GLB,"PSI_GLB_FFF")

DEALLOCATE(PSI_GLB)
CALL DEALLOCATE(PSI)
! ======================================================================

! ======================================================================
CALL ALLOCATE(chi, PPP_SPACE)
ALLOCATE(CHI_GLB(NDIMR, NDIMTH, NDIMX))

DO II=1,SIZE(CHI%E,1)
  DO JJ=1,SIZE(CHI%E,2)
    DO KK=1,SIZE(CHI%E,3)
        CHI%E(II,JJ,KK) = (II+CHI%INR)*10**6+(JJ+CHI%INTH)*10**3+(KK+CHI%INX)
    ENDDO
  ENDDO
ENDDO

CALL MASSEMBLE(CHI%E,CHI_GLB,2); CALL save_glb(CHI_GLB,"CHI_GLB_PPP")

DEALLOCATE(CHI_GLB)
CALL DEALLOCATE(CHI)

! ======================================================================

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

subroutine save_glb(GLB_DATA,FILENAME)
! ======================================================================
   COMPLEX(P8),DIMENSION(:,:,:):: GLB_DATA
   character(*):: FILENAME
   INTEGER:: IM,IK

    IF (MPI_RANK.EQ.0) THEN
        open(UNIT=777,FILE=ADJUSTL(TRIM(FILENAME))//'.dat',&
        &STATUS='UNKNOWN',ACTION='WRITE')

        DO IM = 1,SIZE(GLB_DATA,2)
            DO IK = 1,SIZE(GLB_DATA,3)

            WRITE(777,256) REAL(GLB_DATA(:,IM,IK))

            ENDDO
        ENDDO

        close(777)
    ENDIF

    CALL MPI_BARRIER(MPI_COMM_IVP,IERR)
256 FORMAT((F10.0),*(',',F10.0))
end subroutine save_glb
! ======================================================================

end program test_assemble
