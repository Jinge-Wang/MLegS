program test_norm_lognorm
! ======================================================================
! With large M (e.g. M > 86), the normalization factors approach machine round-ff and eventually become identically zero. This creates problems in mod_lin_legendre, where subroutines use normlization factors to obtain normalized operator matrices or band matrices. To avoid division-by-zero issue, I created corresponding subroutines that use log(norm) rather than norm to obtain normalized operators. 

! This test program creates the (band) operator matrices using both the old and new subroutines and save them in folder ./test2X. 
! To run the program, use the following shell script:
!    srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/test3_exec,
! and each slurm task processor group will create the matrices for certain m and k. 
! To check the result, ./utility/test_legnorm.m can be used to convert the non-band matrices into band matrices, and the processed files can be compared to verify whether the new subroutine generate results that are consistent with the old subroutines. Additionally, ./utility/test_norm.m was also created to assist the creation and verification of the new routines.

! Custom subroutines save_glb_c and save_glb_r2 save complex 3d and real 2d matrices stored either in the current processor group or in the root (#0) processor group. They can be further expanded into an interface.
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
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: PSI_GLB,CHI_GLB
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: RUR_GLB,RUP_GLB,UZ_GLB
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: ROR_GLB,ROP_GLB,OZ_GLB,W_GLB

REAL(P8),DIMENSION(:,:),ALLOCATABLE:: GLB_ARRAY

! =============================== START ================================
CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
   WRITE(*,*) 'The threading support is lesser than that demanded.'
   CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
time_start = mpi_wtime()
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

CALL READCOM('NOECHO')
CALL READIN(5)
CALL LEGINIT()

status(1,1)=0

JJ = MPI_RANK+1
II = NRCHOPS(JJ)

! ======================================================================
ALLOCATE( GLB_ARRAY(NRCHOP,5) )
GLB_ARRAY = 0.D0
GLB_ARRAY(:II,:) = BAND_LEG_RAT_DEL2(II,M(JJ),AK(JJ,JJ),ELL,TFM%NORM(:,JJ))
CALL save_glb_r2(TRANSPOSE(GLB_ARRAY),'BAND_LEG_RAT_DEL2'//ITOA3(JJ),.true.)

GLB_ARRAY = 0.D0
GLB_ARRAY(:II,:) = BAND_LOGLEG_RAT_DEL2(II,M(JJ),AK(JJ,JJ),ELL,TFM%LOGNORM(:,JJ))
CALL save_glb_r2(TRANSPOSE(GLB_ARRAY),'BAND_LOGLEG_RAT_DEL2'//ITOA3(JJ),.true.)
DEALLOCATE(GLB_ARRAY)
! ======================================================================

ALLOCATE( GLB_ARRAY(II,II+1) )
GLB_ARRAY(:,:) = LEG_RAT_DEL2(II, II+1, M(JJ),AK(JJ,JJ),ELL, TFM%NORM(:,JJ))
CALL save_glb_r2(GLB_ARRAY,'MATX_LEG_RAT_DEL2'//ITOA3(JJ),.true.)

GLB_ARRAY(:,:) = LOGLEG_RAT_DEL2(II, II+1, M(JJ),AK(JJ,JJ),ELL, TFM%LOGNORM(:,JJ))
CALL save_glb_r2(GLB_ARRAY,'MATX_LOGLEG_RAT_DEL2'//ITOA3(JJ),.true.)
DEALLOCATE(GLB_ARRAY)
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

subroutine save_glb_c(GLB_DATA,FILENAME,ISERIAL)
! ======================================================================
   COMPLEX(P8),DIMENSION(:,:,:):: GLB_DATA
   character(*):: FILENAME
   INTEGER:: IM,IK

   LOGICAL,OPTIONAL:: ISERIAL

   if ((MPI_RANK.eq.0).OR.(PRESENT(ISERIAL).AND.(ISERIAL))) then
      open(UNIT=777,FILE='./test2X/'//ADJUSTL(TRIM(FILENAME))//'.dat',&
      &STATUS='UNKNOWN',ACTION='WRITE')

      DO IM = 1,SIZE(GLB_DATA,2)
         DO IK = 1,SIZE(GLB_DATA,3)
            WRITE(777,256) GLB_DATA(:,IM,IK)
         ENDDO
      ENDDO
    !   256 FORMAT((S,E24.16E3,SP,E24.16E3,'i'),*(',',S,E24.16E3,SP,E24.16E3,'i'))
      256 FORMAT((S,E24.14E3,SP,E24.14E3,'i'),*(',',S,E24.14E3,SP,E24.14E3,'i'))
      close(777)
   endif

   IF (PRESENT(ISERIAL).AND.(ISERIAL)) RETURN
   call mpi_barrier(MPI_COMM_IVP,IERR)

end subroutine save_glb_c
! ======================================================================

subroutine save_glb_r2(GLB_DATA,FILENAME,ISERIAL)
! ======================================================================
   REAL(P8),DIMENSION(:,:):: GLB_DATA
   character(*):: FILENAME
   INTEGER:: IM

   LOGICAL,OPTIONAL:: ISERIAL

   if ((MPI_RANK.eq.0).OR.(PRESENT(ISERIAL).AND.(ISERIAL))) then
      open(UNIT=777,FILE='./test2X/'//ADJUSTL(TRIM(FILENAME))//'.dat',&
      &STATUS='UNKNOWN',ACTION='WRITE')

      DO IM = 1,SIZE(GLB_DATA,2)
        WRITE(777,256) GLB_DATA(:,IM)
      ENDDO
      256 FORMAT((S,E24.14E3),*(',',S,E24.14E3))
      close(777)
   endif

   IF (PRESENT(ISERIAL).AND.(ISERIAL)) RETURN
   call mpi_barrier(MPI_COMM_IVP,IERR)

end subroutine save_glb_r2
! ======================================================================

end program test_norm_lognorm
