PROGRAM EVP_PRINT_Q
!=======================================================================
! [USAGE]:
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
!
! [UPDATES]:
! LAST UPDATE ON SEP 20, 2022
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
USE MOD_MARCH
USE MOD_INIT
USE MOD_EVP
USE MOD_DIAGNOSTICS

IMPLICIT NONE
INTEGER     :: II, JJ, KK, K_ind
INTEGER     :: M_BAR, EIG_BAR_IND, M_1
REAL(P8)    :: K_BAR
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_RANGE, H_RANGE

COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR, RUP_BAR, UZ_BAR, EIG_R_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_BAR, ROP_BAR, OZ_BAR, EIG_L_BAR

COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR2, RUP_BAR2, UZ_BAR2, EIG_R_BAR2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_BAR2, ROP_BAR2, OZ_BAR2, EIG_L_BAR2

! SCANNING:
CHARACTER(len=6) :: K_VAL, H_VAL
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID

! CALL MPI_INIT(IERR)
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
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS
CALL READCOM('NOECHO')
CALL READIN(5)

! K: -10:0.1:10
ALLOCATE(K_RANGE(201))
CALL LINSPACE(-10.d0,10.d0,K_RANGE)

! H: 0:0.1:2
ALLOCATE(H_RANGE(21))
CALL LINSPACE(0.d0,1.d0,H_RANGE)

! calculate operator matrix and its corresponding eigenvalues
! ======================================================================
DO JJ = 1,SIZE(H_RANGE)

    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    QPAIR%H(1) = H_RANGE(JJ)
    WRITE(H_VAL,'(F06.2)') H_RANGE(JJ)
    IF (ALLOCATED(RUR0)) THEN
        DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
    END IF

    DO M_BAR = 0,10
    DO KK = 1,SIZE(K_RANGE)

        K_BAR = K_RANGE(KK)
        WRITE(K_VAL,'(F06.2)') K_BAR

        ! IF (H_RANGE(JJ).NE.0.35) CYCLE
        ! IF (H_RANGE(JJ).EQ.0.35) THEN
        !     IF (M_BAR.LT.5) THEN
        !         CYCLE
        !     ELSEIF ((M_BAR.EQ.5).AND.(K_BAR.LT.8.90)) THEN
        !         CYCLE
        !     ENDIF
        ! ENDIF
        
        ! WRITE CURRENT 
        IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) RTOA6(QPAIR%H(1)),':',M_BAR,K_BAR

        ! SKIP M = 0, K = 0
        IF ((M_BAR.EQ.0).AND.(K_BAR.EQ.0.D0)) THEN
            IF (MPI_GLB_RANK.EQ.0) THEN
                open(FID,FILE='./data/eig_MK_'//H_VAL//'_'//ITOA3(M_BAR)//'_'//K_VAL &
                        //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
                WRITE(FID,*) 'II = ',1,':',0,'-',.FALSE.
                close(FID)
            ENDIF

            CYCLE
        ENDIF

        ! CALL EIG_MATRIX_SERIAL(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, print_switch=.true.)
        CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, serial_switch=.false.) 

        IF (MPI_GLB_RANK .NE. 0) THEN
            GOTO 129
        ELSE
            ! SAVE EIGVALS
            open(FID,FILE='./data/eig_MK_'//H_VAL//'_'//ITOA3(M_BAR)//'_'//K_VAL &
                    //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
            DO II = 1,SIZE(M_eig)
                ! WRITE(FID,*) REAL(M_eig(II)),',',AIMAG(M_eig(II)),',',EIGRES(EIG_R_mat(:,II),M_BAR)
                WRITE(FID,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
                ! IF (REAL(M_eig(II)).GT.0.10) WRITE(*,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
            ENDDO
            close(FID)

            ! ! SAVE EIGVECS
            ! ALLOCATE(EIG_R_BAR(SIZE(EIG_R_mat,1)))
            ! DO II = 1,SIZE(M_eig)
            !     ! II = 688
            !     EIG_R_BAR = EIG_R_mat(:,II)

            !     ! WRITE(*,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
            !     ! CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
            !     ! call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/CriticalLayer_240605/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
            !     !             //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
            !     ! ! call SAVE_VEC('./converg/CriticalLayer_240605/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
            !     !             ! //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
            !     ! DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)

            !     ! IF (EIGRES(EIG_R_BAR,M_BAR)) THEN
            !     ! ! IF ((EIGRES(EIG_R_BAR)).OR.(II.LE.4)) THEN
            !     !     ! WRITE(*,*) 'II = ',II,':',M_eig(II)
            !     !     CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
            !     !     call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/check5/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
            !     !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
            !     !     call SAVE_VEC('./converg/check5/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
            !     !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
            !     !     DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
            !     ! ENDIF
            ! ENDDO
            ! DEALLOCATE(EIG_R_BAR,M_eig)
        ENDIF

    129     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        ! ENDDO

        ! ENDDO
    ENDDO
    ENDDO
ENDDO

IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

END PROGRAM EVP_PRINT_Q
!=======================================================================
