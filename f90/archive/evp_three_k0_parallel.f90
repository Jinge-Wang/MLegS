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
INTEGER :: II, JJ, KK, R_IND, R_IND_END
! SCAN:
INTEGER :: M_0, M_1, M_2
REAL(P8):: K_0, K_1, K_2, SIG_0, SIG_1, SIG_2, R_CTL, OMEGA_CTL

COMPLEX(P8) :: J1, J2, J0
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig_0, M_eig_1, M_eig_2
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_0, EIG_R_mat_0, EIG_L_mat_0
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1, EIG_R_mat_1, EIG_L_mat_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_2, EIG_R_mat_2, EIG_L_mat_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_0, VECR_0, VECR_12, RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_0_COPY, ROP_0_COPY, OZ_0_COPY
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_1, VECR_1, VECR_01, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_1_COPY, ROP_1_COPY, OZ_1_COPY
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_2, VECR_2, VECR_02, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_2_COPY, ROP_2_COPY, OZ_2_COPY

! LOG:
INTEGER:: IS
CHARACTER(len=6) :: K_VAL, H_VAL
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
R_IND_END = 186
QPAIR%H(1) = 0.D0
K_0 = 0; K_1 = 0; K_2 = 0

! ======================== DYNAMIC SCHEDULING ==========================
! MANAGER:
IF (MPI_GLB_RANK.EQ.0) THEN

    WORKER_TASK = 0
    open(unit=LOGID,FILE='k0_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)

    DO M_0 = 1,9
    DO M_1 = 1,9

        IF (M_1 .LE. M_0) CYCLE
        IF (M_0+M_1 .GT. 10) CYCLE

        WORKER_TASK(1) = M_0
        WORKER_TASK(2) = M_1

        ! Receive from ANY available worker
        CALL MPI_RECV(WORKER,1,MPI_LOGICAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,WORKER_STAT,IERR)
        MPI_REC_ID = WORKER_STAT(MPI_SOURCE)

        IF (WORKER) THEN
            ! Assign task to the worker
            CALL MPI_SEND(WORKER_TASK,2,MPI_INTEGER,MPI_REC_ID,1,MPI_COMM_WORLD,IERR)

            ! Update log
            WRITE(LOGID,*) 'RANK#',ITOA4(MPI_REC_ID),':',ITOA4(M_0),',',ITOA4(M_1)

        ELSE
            CYCLE
        ENDIF

    ENDDO ! M_1
    ENDDO ! M_0

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
    M_0 = 0; M_1 = 0

    DO WHILE (M_0.NE.-1)

        ! Tell manager the proc is available
        WORKER = .TRUE.
        CALL MPI_SEND(WORKER,1,MPI_LOGICAL,0,1,MPI_COMM_WORLD,IERR)

        ! Receive the task
        CALL MPI_RECV(WORKER_TASK,2,MPI_INTEGER,0,1,MPI_COMM_WORLD,WORKER_STAT,IERR)
        IF (WORKER_TASK(1).EQ.-1) EXIT ! ALL TASKS HAVE BEEN FINISHED
        WRITE(*,*) 'RECV BY RANK#',ITOA4(MPI_GLB_RANK),':',WORKER_TASK(:)
        M_0 = WORKER_TASK(1); M_1 = WORKER_TASK(2); M_2 = M_0 + M_1

        ! Initilization
        IF (.NOT. ASSOCIATED(TFM%R)) THEN
            NTH = 2
            NX = 4
            NTCHOP = 2 ! M = 0, M = MREAD
            NXCHOP = 2 ! K = 0, K = AKREAD
            ZLEN = 1.D0; CALL LEGINIT(newcomm)
        END IF
       
        ! Triadic Resonance
        DO R_IND = 1, R_IND_END

            ! CRITICAL LAYER
            R_CTL =  TFM%R(R_IND) ! 0.21
            OMEGA_CTL = (1-EXP(-R_CTL**2))/(R_CTL**2)
            SIG_0 = -M_0*OMEGA_CTL; SIG_1 = -M_1*OMEGA_CTL; SIG_2 = -M_2*OMEGA_CTL

            ! CALCULATE EVP FOR EACH TRIAD
            CALL EIG_MATRIX(M_0, K_0, M_mat_0, M_eig_0, EIG_R_mat_0, EIG_L_mat_0, comm_grp=newcomm, serial_switch = .true.) 
            CALL EIG_MATRIX(M_1, K_1, M_mat_1, M_eig_1, EIG_R_mat_1, EIG_L_mat_1, comm_grp=newcomm, serial_switch = .true.)
            CALL EIG_MATRIX(M_2, K_2, M_mat_2, M_eig_2, EIG_R_mat_2, EIG_L_mat_2, comm_grp=newcomm, serial_switch = .true.)

            II = MINLOC(ABS(AIMAG(M_eig_0)-SIG_0),DIM=1)
            JJ = MINLOC(ABS(AIMAG(M_eig_1)-SIG_1),DIM=1)
            KK = MINLOC(ABS(AIMAG(M_eig_2)-SIG_2),DIM=1)

            ! {M0,K0}
            CALL EIG2VELVOR(M_0, K_0, II, EIG_R_mat_0, EIG_L_mat_0, RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0, VECR_0, VECL_0,comm_grp=newcomm)
            ALLOCATE(ROR_0_COPY(NDIMR), ROP_0_COPY(NDIMR), OZ_0_COPY(NDIMR))
            ROR_0_COPY = ROR_0; ROP_0_COPY = ROP_0; OZ_0_COPY = OZ_0

            ! {M1,K1}
            CALL EIG2VELVOR(M_1, K_1, JJ, EIG_R_mat_1, EIG_L_mat_1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1, VECR_1, VECL_1, comm_grp=newcomm)
            ALLOCATE(ROR_1_COPY(NDIMR),ROP_1_COPY(NDIMR),OZ_1_COPY(NDIMR))
            ROR_1_COPY = ROR_1; ROP_1_COPY = ROP_1; OZ_1_COPY = OZ_1

            ! {M2,K2}
            CALL EIG2VELVOR(M_2, K_2, KK, EIG_R_mat_2, EIG_L_mat_2, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, VECR_2, VECL_2, comm_grp=newcomm)
            ALLOCATE(ROR_2_COPY(NDIMR),ROP_2_COPY(NDIMR),OZ_2_COPY(NDIMR))
            ROR_2_COPY = ROR_2; ROP_2_COPY = ROP_2; OZ_2_COPY = OZ_2

            ! NONLINEAR TERMS:
            ! {M2,K2}: 0  X 1 -> 2
            CALL NONLIN_MK(M_2,K_2,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,RUR_1,RUP_1,UZ_1,ROR_1_COPY,ROP_1_COPY,OZ_1_COPY,VECR_01,comm_grp=newcomm)

            ! {M1,K1}: 0* X 2 -> 1
            CALL NONLIN_MK(M_1,K_1,CONJG(RUR_0),CONJG(RUP_0),CONJG(UZ_0),CONJG(ROR_0),CONJG(ROP_0),CONJG(OZ_0),RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_02,comm_grp=newcomm)

            ! {M0,K0}: 1* X 2 -> 0
            CALL NONLIN_MK(M_0,K_0,CONJG(RUR_1),CONJG(RUP_1),CONJG(UZ_1),CONJG(ROR_1),CONJG(ROP_1),CONJG(OZ_1),RUR_2,RUP_2,UZ_2,ROR_2_COPY,ROP_2_COPY,OZ_2_COPY,VECR_12,comm_grp=newcomm)
            
            ! ============================== RESULT ================================
            ! WRITE (*, *) "NEARLY DEGEN:"
            J0 = DOT_PRODUCT(VECL_0, VECR_12); ! WRITE (*,*) "J0: ", J0
            J1 = DOT_PRODUCT(VECL_1, VECR_02); ! WRITE (*,*) "J1: ", J1
            J2 = DOT_PRODUCT(VECL_2, VECR_01); ! WRITE (*,*) "J2: ", J2
            ! WRITE(*,*) "J1 + J0 - J2: ", J1 + J0 - J2

            IF (ABS(J0)+ABS(J1)+ABS(J2) .GT. 10**-16) THEN
                OPEN(unit=FID, FILE='./data/triad_'//ITOA3(M_0)//'_'//ITOA3(M_1) &
                        //'_NRCHOP_'//ITOA3(NRCHOP)//'.output', STATUS='unknown', ACTION='WRITE', ACCESS='APPEND', IOSTAT=IS)
                WRITE(FID,165) M_0, M_1, M_2, R_CTL, M_eig_0(II)/M_0, J0, J1, J2
                165  FORMAT(3(I3,","),1(F16.13),4(",",E20.13,SP,E20.13,"i"))
                CLOSE(FID)
            ENDIF

            ! DEALLOCATE
            !=======================================================================
        1   DEALLOCATE(RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0)
            DEALLOCATE(RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1)
            DEALLOCATE(RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2)
            DEALLOCATE(ROR_0_COPY, ROP_0_COPY, OZ_0_COPY)
            DEALLOCATE(ROR_1_COPY, ROP_1_COPY, OZ_1_COPY)
            DEALLOCATE(ROR_2_COPY, ROP_2_COPY, OZ_2_COPY)
            DEALLOCATE(VECR_0, VECL_0, VECR_1, VECR_2, VECL_1, VECL_2, VECR_01, VECR_02, VECR_12)
            DEALLOCATE(EIG_R_mat_0, EIG_L_mat_0, EIG_R_mat_1, EIG_L_mat_1, EIG_R_mat_2, EIG_L_mat_2)
            DEALLOCATE(M_eig_0, M_eig_1, M_eig_2)
            IF (ALLOCATED(M_mat_1)) THEN
                DEALLOCATE(M_mat_0,M_mat_1,M_mat_2)
            ENDIF

        ENDDO ! R_IND

    ENDDO

    WRITE(*,*) 'RANK#',ITOA4(MPI_GLB_RANK),' - FINISHED'

ENDIF

IF(ALLOCATED(RUR0)) DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
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
