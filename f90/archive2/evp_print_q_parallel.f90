PROGRAM EVP_PRINT_PARALLEL
!=======================================================================
! [USAGE]: 
! Parallel version of evp_print
!
! [UPDATES]:
! LAST UPDATE ON DEC 3, 2020
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

IMPLICIT NONE

! SCAN:
INTEGER     :: II, JJ, KK, K_ind
INTEGER     :: M_BAR, EIG_BAR_IND
REAL(P8)    :: K_BAR
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_RANGE, H_RANGE

COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR, RUP_BAR, UZ_BAR, EIG_R_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_BAR, ROP_BAR, OZ_BAR, EIG_L_BAR

! LOG:
INTEGER:: IS
CHARACTER(len=6) :: K_VAL, H_VAL
CHARACTER(LEN=72):: FILENAME
INTEGER:: LOGID, FID

! SCHEDULING:
INTEGER:: MPI_REC_ID, WORKER_STAT(MPI_STATUS_SIZE), REQUEST, WORKER_TASK(2)
LOGICAL:: WORKER

CALL MPI_INIT(IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
IF (MPI_GLB_RANK.EQ.0) THEN
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
ENDIF
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS
CALL READCOM('NOECHO')
CALL READIN(5)

! ============================= ALL PROCS ==============================
! K: -10:0.1:10
ALLOCATE(K_RANGE(201))
CALL LINSPACE(-10.d0,10.d0,K_RANGE)

! H: 0:0.1:2
ALLOCATE(H_RANGE(21))
CALL LINSPACE(0.d0,1.d0,H_RANGE)

! ======================== DYNAMIC SCHEDULING ==========================
! MANAGER:
IF (MPI_GLB_RANK.EQ.0) THEN

    WORKER_TASK = 0
    open(unit=LOGID,FILE='q_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)

    DO JJ = 1,SIZE(H_RANGE)
    DO II = 0,10

        WORKER_TASK(1) = JJ ! INDEX OF H
        WORKER_TASK(2) = II ! M

        ! Receive from ANY available worker
        CALL MPI_RECV(WORKER,1,MPI_LOGICAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,WORKER_STAT,IERR)
        MPI_REC_ID = WORKER_STAT(MPI_SOURCE)

        IF (WORKER) THEN
            ! Assign task to the worker
            CALL MPI_SEND(WORKER_TASK,2,MPI_INTEGER,MPI_REC_ID,1,MPI_COMM_WORLD,IERR)

            ! Update log
            WRITE(LOGID,*) 'RANK#',ITOA4(MPI_REC_ID),':',RTOA6(H_RANGE(JJ)),II

        ELSE
            CYCLE
        ENDIF

    ENDDO ! II
    ENDDO ! JJ

    ! Notify all worker job is done
    DO II = 1,MPI_GLB_PROCS-1
        WORKER_TASK(1) = -1
        ! CALL MPI_ISEND(WORKER_TASK,2,MPI_INTEGER,II,1,MPI_COMM_WORLD,REQUEST,IERR)
        CALL MPI_SEND(WORKER_TASK,2,MPI_INTEGER,II,1,MPI_COMM_WORLD,REQUEST,IERR) ! BLOCKING SEND
    ENDDO
    close(LOGID)
    WRITE(*,*) 'RANK#',ITOA4(MPI_GLB_RANK),' - FINISHED'

! WORKER:
ELSE

    ! Initilization
    JJ = 0

    DO WHILE (JJ.NE.-1)

        ! Tell manager the proc is available
        WORKER = .TRUE.
        CALL MPI_SEND(WORKER,1,MPI_LOGICAL,0,1,MPI_COMM_WORLD,IERR)

        ! Receive the task
        CALL MPI_RECV(WORKER_TASK,2,MPI_INTEGER,0,1,MPI_COMM_WORLD,WORKER_STAT,IERR)
        IF (WORKER_TASK(1).EQ.-1) EXIT ! ALL TASKS HAVE BEEN FINISHED
        WRITE(*,*) 'RECV BY RANK#',ITOA4(MPI_GLB_RANK),': H = ',H_RANGE(WORKER_TASK(1)),', M = ',WORKER_TASK(2)
        JJ = WORKER_TASK(1)
        
        ! SET H
        IF (QPAIR%H(1) .NE. H_RANGE(JJ)) THEN
            QPAIR%H(1) = H_RANGE(JJ)
            WRITE(H_VAL,'(F06.2)') H_RANGE(JJ)
            IF (ALLOCATED(RUR0)) THEN
                DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
            END IF
        ENDIF

        ! CALCULATE EVP
        ! DO M_BAR = 0,10
        M_BAR = WORKER_TASK(2)
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
    
            ! SKIP M = 0, K = 0
            IF ((M_BAR.EQ.0).AND.(K_BAR.EQ.0.D0)) THEN
                open(FID,FILE='./data/data_eig/eig_MK_'//H_VAL//'_'//ITOA3(M_BAR)//'_'//K_VAL &
                        //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
                WRITE(FID,*) 'II = ',1,':',0,'-',.FALSE.
                close(FID)

                CYCLE
            ENDIF

            ! CALL EIG_MATRIX_SERIAL(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, print_switch=.true.)
            CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, serial_switch=.true.) 

            ! SAVE EIGVALS
            open(FID,FILE='./data/data_eig/eig_MK_'//H_VAL//'_'//ITOA3(M_BAR)//'_'//K_VAL &
                    //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
            DO II = 1,SIZE(M_eig)
                ! WRITE(FID,*) REAL(M_eig(II)),',',AIMAG(M_eig(II)),',',EIGRES(EIG_R_mat(:,II),M_BAR)
                WRITE(FID,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_mat(:,II),M_BAR)
                ! IF (REAL(M_eig(II)).GT.0.10) WRITE(*,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
            ENDDO
            close(FID)

            ! SAVE EIGVECS
            ALLOCATE(EIG_R_BAR(SIZE(EIG_R_mat,1)))
            II = MAXLOC(real(M_eig),1)
            ! DO II = 1,SIZE(M_eig)
                ! II = 688
                EIG_R_BAR = EIG_R_mat(:,II)

                ! WRITE(*,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
                CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
                call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR,'./data/data_vel/vel_MK_'//H_VAL//'_'//ITOA3(M_BAR)//'_'//K_VAL &
                        //'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
                ! ! call SAVE_VEC('./converg/CriticalLayer_240605/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
                !             ! //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
                ! DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)

                ! IF (EIGRES(EIG_R_BAR,M_BAR)) THEN
                ! ! IF ((EIGRES(EIG_R_BAR)).OR.(II.LE.4)) THEN
                !     ! WRITE(*,*) 'II = ',II,':',M_eig(II)
                !     CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
                !     call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/check5/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
                !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
                !     call SAVE_VEC('./converg/check5/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
                !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
                DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
                ! ENDIF
            ! ENDDO
            DEALLOCATE(EIG_R_BAR)


        ENDDO
        ! ENDDO

    ENDDO

    WRITE(*,*) 'RANK#',ITOA4(MPI_GLB_RANK),' - FINISHED'

ENDIF

IF(ALLOCATED(RUR0)) DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

IF (MPI_GLB_RANK.EQ.0) THEN
    WRITE(*,*) ''
    WRITE(*,*) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
ENDIF
CALL MPI_FINALIZE(IERR)

! ======================================================================
END PROGRAM EVP_PRINT_PARALLEL
!=======================================================================
