PROGRAM EVP_PRINT
!=======================================================================
! [USAGE]:
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
!=======================================================================
USE omp_lib
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
USE MOD_BOUSSINESQ
USE MOD_MARCH
USE MOD_INIT
USE MOD_EVP
USE MOD_DIAGNOSTICS

IMPLICIT NONE
INTEGER     :: II, JJ, KK
INTEGER     :: M_BAR, NR_MK
REAL(P8)    :: K_BAR

COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
! COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR, RUP_BAR, UZ_BAR, EIG_R_BAR
! COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_BAR, ROP_BAR, OZ_BAR, EIG_L_BAR

! SCANNING:
CHARACTER(len=6) :: K_VAL, H_VAL
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID
integer :: num_args
character(len=12), dimension(:), allocatable :: args

! MPI INITILIZATION
CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, MPI_THREAD_MODE, IERR)
IF (MPI_THREAD_MODE .LT. MPI_THREAD_SERIALIZED) THEN
    WRITE (*, *) 'The threading support is lesser than that demanded.'
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
END IF

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
END IF

! PARALLEL COMPUTATION
! - SPLIT THE FORMATION OF H INTO DIFFERENT PROCS
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS

! OBTAIN BASIC INFO OF THE FLOW FIELD
CALL READCOM('NOECHO')
CALL READIN(5)

! M AND K OF THE EVP & H OF THE Q-VORTEX
! - READ FROM COMMAND INPUT
num_args = command_argument_count()
IF (num_args .ge. 2) THEN
    allocate(args(num_args))
    call get_command_argument(1,args(1))
    read (args(1),'(I12)') M_BAR
    call get_command_argument(2,args(2))
    read (args(2),'(F12.6)') K_BAR
    IF (num_args .ge. 3) THEN
        call get_command_argument(3,args(3))
        read (args(3),'(F12.6)') QPAIR%H(1)
    ENDIF
! - DIRECT SETUP
ELSE
    M_BAR = 1
    K_BAR = 0.d0
ENDIF  

! I/O
WRITE(K_VAL,'(F06.2)') K_BAR
WRITE(H_VAL,'(F06.2)') QPAIR%H(1)
! IF (MPI_GLB_RANK.EQ.0) 
WRITE(*,*) MPI_GLB_RANK,': M = ',M_BAR,'; K = ',K_BAR

! BOUSSINESQ
BSNSQ%BV0 = 5.0
BSNSQ%OMEGA = -1.d0; !-1.2 

! OBTAIN EVP MATRIX, EIGENVALUE, AND EIGENVECTORS
! CALL EIG_MATRIX_SERIAL(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, print_switch=.true.)
CALL EIG_MATRIX_BSNSQ(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, serial_switch=.false.) 
NR_MK = NRCHOPS(2)

! OUTPUT
IF (MPI_GLB_RANK .NE. 0) THEN
    GOTO 129
ELSE
    ! SAVE EIGVALS
    ! open(FID,FILE='./converg/CriticalLayer_240618/qvortex_0.1/eig_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !         //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
    open(FID,FILE='./data/bsnsq_eig_MK_'//H_VAL//'_'//ITOA3(M_BAR)//'_'//&
    K_VAL//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',&
    ACTION='WRITE',IOSTAT=IERR)
    DO II = 1,SIZE(M_eig)
        ! WRITE(FID,*) 'II = ',II,':',REAL(M_eig(II)),AIMAG(M_eig(II))+M_BAR*BSNSQ%OMEGA,'-',EIGRES(EIG_R_mat(:,II),M_BAR)
        WRITE(FID,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_mat(:2*NR_MK,II),M_BAR)
    ENDDO
    close(FID)

    ! TEST SAVEPERTURB
    II = 1
    CALL SAVE_PERTURB_BSNSQ('./data/new_perturb.input', M_BAR, K_BAR, EIG_R_mat(:3*NR_MK,II), M_eig(II))

    ! ! SAVE EIGVECS
    ! ALLOCATE(EIG_R_BAR(SIZE(EIG_R_mat,1)))
    ! DO II = 1,SIZE(M_eig)

    !     EIG_R_BAR = EIG_R_mat(:,II)

    !     ! WRITE(*,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
    !     ! CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
    !     ! call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './qvortex/data/temp/temp_vel/vel_MK_'//H_VAL//'_'&
    !     !     //ITOA3(M_BAR)//'_'//K_VAL//'_IND_'//ITOA3(II) &
    !     !     //'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    !     ! ! call SAVE_VEC('./converg/CriticalLayer_240605/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !     !             ! //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
    !     DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
    ! ENDDO
    ! DEALLOCATE(EIG_R_BAR,M_eig)
ENDIF

129     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)


IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

END PROGRAM EVP_PRINT
!=======================================================================
