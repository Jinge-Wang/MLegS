PROGRAM EVP_ALMOST_DISCRETE_SCAN5
!=======================================================================
! [USAGE]: 
! Faster Version: the bar term is calculated once for each wavenumber
! Modified to include M=0 mode
! Check resolved or not using EIGRES
! 
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
! USE VPROD_PFF
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
    INTEGER     :: N_mat1,N_mat2,II,JJ,KK
    INTEGER     :: M_IND1,M_IND2
    INTEGER     :: M_BAR,M_1,M_2,EIG_BAR_IND
    REAL(P8)    :: K_BAR,K_VAL_1_0,K_VAL_1_1,K_VAL_2_0,K_VAL_2_1
    REAL(P8)    :: K_shift,tol,tol_diff,K_delta
    REAL(P8)    :: EIG_LMR_1, ZLEN_BAR
    COMPLEX(P8) :: EIG_BAR, EIG_PRIME, EIG_PRIME_R, ALPHA, BETA
    REAL(P8),DIMENSION(:),ALLOCATABLE:: K_tune, K_tune_JJ, tol_diff_JJ
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_1_0, M_eig_1_1
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_2_0, M_eig_2_1
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_diff_1,M_mat_diff_2
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_2_0, EIG_R_mat_2_0, EIG_L_mat_2_0
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_BAR,VECR_BAR,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_1, VECR_1, VECR_B1, RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_2, VECR_2, VECR_B2, RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2

    ! ADD PERTURB:
    INTEGER:: IS,NDIM
    REAL(P8):: rp1,ip1
    ! SCANNING:
    CHARACTER(LEN=72):: FILENAME
    INTEGER:: FID, LOGID, K_BAR_SCAN, K_1_SCAN, JUMP_FLAG


    CALL MPI_INIT(IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
    IF (MPI_GLB_RANK.EQ.0) THEN
        WRITE(*,*) 'PROGRAM STARTED'
        CALL PRINT_REAL_TIME()   ! @ MOD_MISC
        ! open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
    ENDIF
    CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS
    CALL READCOM('NOECHO')
    CALL READIN(5)

! ======================================================================
! ============================== DO LOOP ===============================
! {M_bar, AK_bar}: PERTURBATION
!=======================================================================    
DO M_BAR = 0,10; DO K_BAR_SCAN = -10.D0,10.D0,1.D0
    
    IF ((M_BAR.EQ.0).AND.(K_BAR_SCAN.EQ.0.D0)) CYCLE
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    K_BAR = REAL(K_BAR_SCAN)
    IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) 'M_BAR,KBAR: ',ITOA3(M_BAR),',',ITOA3(K_BAR_SCAN)

    IF (ALLOCATED(M_mat)) DEALLOCATE(M_mat)
    IF (ALLOCATED(M_eig)) DEALLOCATE(M_eig)
    IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
    IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)    

    ! calculate operator matrix and its corresponding eigenvalues
    CALL EIG_MATRIX(M_BAR,K_BAR,M_mat,M_eig,EIG_R_mat,EIG_L_mat,comm_grp=newcomm)

! ======================================================================
loop_ind: DO EIG_BAR_IND = 1,size(M_mat,1)

    JUMP_FLAG = 0

    ! Pick a target Sigma_bar
    IF (MPI_GLB_RANK.EQ.0) THEN

        EIG_BAR = M_eig(EIG_BAR_IND)
        IF (.NOT.((EIGRES(EIG_R_mat(:,EIG_BAR_IND),M_BAR)))) THEN
            JUMP_FLAG = 39
        ELSE
            WRITE(*,*) 'TARGET EIG_BAR:', EIG_BAR, '(# ', EIG_BAR_IND, ')'
        ENDIF
        
    ENDIF
    CALL MPI_BCAST(JUMP_FLAG,1,MPI_INT,0,MPI_COMM_WORLD,IERR)

    IF (JUMP_FLAG.EQ.39) THEN
        ! IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) 'CYCLE IND'
        CYCLE loop_ind
    ENDIF
    ! ======================================================================

    IF (MPI_GLB_RANK.EQ.0) THEN
        ! Save rur_bar, rup_bar, uz_bar, ror_bar, rop_bar, oz_bar for DEGENERATE PERTURBATION THEORY
        CALL EIG2VELVOR(M_BAR,K_BAR,EIG_BAR_IND,EIG_R_mat,EIG_L_mat,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECR_BAR,VECL_BAR,comm_grp=newcomm)
    ENDIF

! ! {M_1, AK_1} & {M_2}
! !=======================================================================
DO M_1 = -M_BAR/2,(10-M_BAR); DO K_1_SCAN = -10.D0,10.D0,1.D0

    ! IF ((M_BAR.EQ.M_1).AND.(K_BAR_SCAN.EQ.K_1_SCAN)) CYCLE
    IF ((M_1.EQ.0).AND.(K_1_SCAN.EQ.0.D0)) CYCLE                    ! Not interested in m = 0, k = 0 mode
    IF ((M_1+M_BAR.EQ.0).AND.(K_1_SCAN+K_BAR_SCAN.EQ.0.D0)) CYCLE   ! Not interested in m = 0, k = 0 mode
    IF ((M_1.EQ.0).AND.(M_BAR.EQ.0.D0).AND.(K_BAR.LE.0.D0)) CYCLE   ! for m_b = m1 = m2 = 0, only need to check positive k_1

    IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) 'M_1,K_1: ', ITOA3(M_1),',',ITOA3(K_1_SCAN)
    K_VAL_1_0 = REAL(K_1_SCAN)
    JUMP_FLAG = 0

    CALL MPI_BCAST(EIG_BAR,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(M_1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    M_2 = M_1 + M_BAR

    ! K_VAL_1_0 and K_VAL_1_1 are both for M1,AK1
    K_delta = 0.001
    tol = 1.0E-8 !13  ! tolerance

    DO II = 1,6
    ! For M1,K1:
        K_VAL_1_1 = K_VAL_1_0 + K_delta

        CALL EIG_MATRIX(M_1,K_VAL_1_0,M_mat_1_0,M_eig_1_0,EIG_R_mat_1_0,EIG_L_mat_1_0,comm_grp=newcomm)
        CALL EIG_MATRIX(M_1,K_VAL_1_1,M_mat_1_1,comm_grp=newcomm)

    ! For M2,K2:
        K_VAL_2_0 = K_VAL_1_0 + K_BAR
        K_VAL_2_1 = K_VAL_2_0 + K_delta

        CALL EIG_MATRIX(M_2,K_VAL_2_0,M_mat_2_0,M_eig_2_0,EIG_R_mat_2_0,EIG_L_mat_2_0,comm_grp=newcomm)
        CALL EIG_MATRIX(M_2,K_VAL_2_1,M_mat_2_1,comm_grp=newcomm)

    ! MASTER RANK:
    IF (MPI_GLB_RANK.EQ.0) THEN

        IF(.NOT.ALLOCATED(M_mat_diff_1)) THEN
            N_mat1 = SIZE(EIG_R_mat_1_0,1)
            ALLOCATE(M_mat_diff_1(N_mat1,N_mat1))
        ENDIF
        M_mat_diff_1 = 0.D0
        M_mat_diff_1 = MATMUL((M_mat_1_1 - M_mat_1_0), EIG_R_mat_1_0)

        IF(.NOT.ALLOCATED(M_mat_diff_2)) THEN
            N_mat2 = SIZE(EIG_R_mat_2_0,1)
            ALLOCATE(M_mat_diff_2(N_mat2,N_mat2))
        ENDIF
        M_mat_diff_2 = 0.D0
        M_mat_diff_2 = MATMUL((M_mat_2_1 - M_mat_2_0), EIG_R_mat_2_0)        

    ! Find possible K_shift
        IF(.NOT.ALLOCATED(K_tune)) ALLOCATE(K_tune(N_mat1))
        IF(.NOT.ALLOCATED(K_tune_JJ)) ALLOCATE(K_tune_JJ(N_mat2))

        !$OMP parallel do private(EIG_LMR_1,KK,K_tune_JJ)
        DO JJ = 1,N_mat1
            ! eigc_b = -(m_bar*Q+k_bar*H);
            IF (.NOT.((EIGRES(EIG_R_mat_1_0(:,JJ),M_1)))) THEN          ! UNRESOLVED MODES
                K_tune(JJ) = 100.D0
            ! ELSEIF ((M_BAR.EQ.M_1).AND.(ABS(K_BAR_SCAN-K_1_SCAN).LE.0.1)) THEN   
            !     IF (ABS(EIG_BAR-M_eig_1_0(JJ)).LE.0.01) K_tune(JJ) = 100.D0
            ELSEIF ((M_BAR.EQ.M_1).AND.(ABS(K_BAR_SCAN-K_VAL_1_0).LE.0.1)) THEN   
                IF (ABS(EIG_BAR-M_eig_1_0(JJ)).LE.0.01) K_tune(JJ) = 100.D0
            ELSE
                EIG_LMR_1 = 0.D0
                EIG_LMR_1 = &
                & AIMAG(DOT_PRODUCT(EIG_L_mat_1_0(:,JJ),M_mat_diff_1(:,JJ))/DOT_PRODUCT(EIG_L_mat_1_0(:,JJ),EIG_R_mat_1_0(:,JJ)))
    
                K_tune_JJ = 0.D0
                DO KK = 1,N_mat2
                    IF (.NOT.((EIGRES(EIG_R_mat_2_0(:,KK),M_2)))) THEN
                        K_tune_JJ(KK) = 100.D0
                    ELSE
                        K_tune_JJ(KK) = K_delta*AIMAG(M_eig_2_0(KK) - M_eig_1_0(JJ) - EIG_BAR) / ( &
                                    & EIG_LMR_1 - &
                                    & AIMAG(DOT_PRODUCT(EIG_L_mat_2_0(:,KK),M_mat_diff_2(:,KK))/DOT_PRODUCT(EIG_L_mat_2_0(:,KK),EIG_R_mat_2_0(:,KK))))
                    ENDIF
                ENDDO
                K_tune(JJ) = K_tune_JJ(MINLOC(ABS(K_tune_JJ),1))
            ENDIF
        ENDDO
        !$OMP end parallel do

        K_shift = K_tune(MINLOC(ABS(K_tune),1))

        WRITE(*,*) 'step#',ITOA3(II),'- K_shift:', K_shift
        IF (ABS(K_shift).GT.1) JUMP_FLAG = 39 !CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    ENDIF
    CALL MPI_BCAST(JUMP_FLAG,1,MPI_INT,0,MPI_COMM_WORLD,IERR)
    IF (JUMP_FLAG.EQ.39) GOTO 10

    ! Calculate the actual 
        K_VAL_1_1 = K_VAL_1_0 + K_shift
        CALL EIG_MATRIX(M_1,K_VAL_1_1,M_mat_1_1,M_eig_1_1,EIG_R_mat_1_1,comm_grp=newcomm)

        K_VAL_2_1 = K_VAL_1_1 + K_BAR
        CALL EIG_MATRIX(M_2,K_VAL_2_1,M_mat_2_1,M_eig_2_1,EIG_R_mat_2_1,comm_grp=newcomm)

    ! Check if tolerance is reached
    ! Match the IMAG part
    IF (MPI_GLB_RANK.EQ.0) THEN

        IF(.NOT.ALLOCATED(tol_diff_JJ)) ALLOCATE(tol_diff_JJ(N_mat1))
        DO KK = 1,N_mat2
            IF (.NOT.((EIGRES(EIG_R_mat_2_1(:,KK),M_2)))) M_eig_2_1(KK) = 999.D0*IU
        ENDDO
        DO JJ = 1,N_mat1
            IF (.NOT.((EIGRES(EIG_R_mat_1_1(:,JJ),M_1)))) THEN
                tol_diff_JJ(JJ) = 100.D0
                CYCLE
            ENDIF
            tol_diff_JJ(JJ) = MINVAL(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)))
        ENDDO

        tol_diff = MINVAL(ABS(tol_diff_JJ))
        tol_diff = tol_diff/ABS(AIMAG(EIG_BAR))
        WRITE(*,*) tol_diff
        IF ((tol_diff.LT.tol).OR.(K_shift.LE.1.0E-4)) THEN
            WRITE(*,*) 'Desired tol reached'
            WRITE(*,*) 'tol_diff:', tol_diff, ', Number of steps:', II
            tol = 999.0
        ENDIF
        K_VAL_1_0 = K_VAL_1_1
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    IF (tol.EQ.999.0) goto 20

    ENDDO

! Maximum number of steps allowed
!=======================================================================
    IF (II.EQ.6) THEN ! Exceeds maximum number of steps allowed
        IF (MPI_GLB_RANK.EQ.0) THEN
            WRITE(*,*) 'Seem trapped'
            WRITE(*,*) 'tol_diff:', tol_diff
        ENDIF
        GOTO 10!CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    ENDIF

! Saving and outputing the results in FFF Space
!=======================================================================   
    ! Find JJth Eigvector of {M1,K1} & KKth Eigvector of {M2,K2} 
20  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    CALL EIG_MATRIX(M_1,K_VAL_1_1,M_mat_1_1,M_eig_1_1,EIG_R_mat_1_1,EIG_L_mat_1_1,comm_grp=newcomm)
    CALL EIG_MATRIX(M_2,K_VAL_2_1,M_mat_2_1,M_eig_2_1,EIG_R_mat_2_1,EIG_L_mat_2_1,comm_grp=newcomm)

    IF (MPI_GLB_RANK.EQ.0) THEN

        DO KK = 1,N_mat2
            IF (.NOT.((EIGRES(EIG_R_mat_2_1(:,KK),M_2)))) M_eig_2_1(KK) = 999.D0*IU
        ENDDO
        DO JJ = 1,N_mat1
            IF (.NOT.((EIGRES(EIG_R_mat_1_1(:,JJ),M_1)))) THEN
                tol_diff_JJ(JJ) = 100.D0
                K_tune(JJ) = 100
                CYCLE
            ENDIF
            tol_diff_JJ(JJ) = MINVAL(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)))
            K_tune(JJ) = MINLOC(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)),1)
        ENDDO
        ! DO JJ = 1,N_mat1

        !     tol_diff_JJ(JJ) = MINVAL(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)))
        !     K_tune(JJ) = MINLOC(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)),1)
        ! ENDDO
        
        JJ = MINLOC(ABS(tol_diff_JJ),1)
        KK = K_tune(JJ)
        !CALL SORT(0,M_eig_1_1)
        ! CALL MSAVE(M_eig_1_1,'M_eig_1_1.dat')
        !CALL SORT(0,M_eig_2_1)
        ! CALL MSAVE(M_eig_2_1,'M_eig_2_1.dat')

        WRITE(*,*) 'DEGEN:'
        WRITE(*,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
        WRITE(*,201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
        WRITE(*,201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)

201  FORMAT('Eig: M# = ',I3,'; AK = ',F20.13,'; Sig (#',I3') =',1P8E20.13)
    !=======================================================================   

    ! Calculate the FINAL RESULT:
    ! M_bar, K_bar: RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECL_BAR

    ! {M1,K1}
    CALL EIG2VELVOR(M_1,K_VAL_1_1,JJ,EIG_R_mat_1_1,EIG_L_mat_1_1,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_1,VECL_1,comm_grp=newcomm)

    ! {M2,K2}
    CALL EIG2VELVOR(M_2,K_VAL_2_1,KK,EIG_R_mat_2_1,EIG_L_mat_2_1,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_2,VECL_2,comm_grp=newcomm)
    CALL NONLIN_MK(M_2,K_VAL_2_1,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_B1,comm_grp=newcomm)
    ! CALL MCAT(UZ_2)

    ! {M1,K1}
    ALLOCATE(RUR_BARC(NDIMR),RUP_BARC(NDIMR),UZ_BARC(NDIMR),ROR_BARC(NDIMR),ROP_BARC(NDIMR),OZ_BARC(NDIMR))
    RUR_BARC = CONJG(RUR_BAR); RUP_BARC = CONJG(RUP_BAR); UZ_BARC = CONJG( UZ_BAR)
    ROR_BARC = CONJG(ROR_BAR); ROP_BARC = CONJG(ROP_BAR); OZ_BARC = CONJG( OZ_BAR)
    CALL NONLIN_MK(M_1,K_VAL_1_1,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_B2,comm_grp=newcomm)
    ! call mcat(VECR_B2)

    ! Sigma(1):
    EIG_PRIME_R = DOT_PRODUCT(VECL_1,VECR_B2)*DOT_PRODUCT(VECL_2,VECR_B1)
    WRITE(*,*) ''
    WRITE(*,*) " SIGMA(degen)'^2= ",EIG_PRIME_R
    WRITE(*,*) " SIGMA(degen)'^2= ",EIG_PRIME_R/(DOT_PRODUCT(VECL_1,VECR_1)*DOT_PRODUCT(VECL_2,VECR_2))
    WRITE(*,*) " SIGMA(degen)'  = ",EIG_PRIME_R**0.5
    ! Beta/Alpha: ALPHA * {M1,K1} + BETA * {M2,K2}
    ALPHA = 1.D0
    ! BETA = (DOT_PRODUCT(VECL_2,VECR_B1)/DOT_PRODUCT(VECL_1,VECR_B2))**0.5
    BETA = EIG_PRIME_R**0.5/DOT_PRODUCT(VECL_1, VECR_B2)
    EIG_R_mat_1_1(:,JJ) = EIG_R_mat_1_1(:,JJ)*ALPHA
    EIG_R_mat_2_1(:,KK) = EIG_R_mat_2_1(:,KK)*BETA
    WRITE(*,*) " BETA(degen)    = ",BETA

    ! ======================================================================
        IF (abs(REAL(EIG_PRIME_R**0.5)).GT.0.01) THEN
            open(unit=FID,FILE='perturb_scan.GROW',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
            WRITE(FID,*) ''
            WRITE(FID,*) 'DEGEN: *GROWTH*'
        ELSE
            open(unit=FID,FILE='perturb_scan.NEGLIGIBLE',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
            WRITE(FID,*) ''
            WRITE(FID,*) 'DEGEN: *NEGLIGIBLE*'
        ENDIF
        WRITE(FID,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
        WRITE(FID,201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
        WRITE(FID,201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)
        WRITE(FID,*) "SIGMA(degen)'= RE:",abs(real(EIG_PRIME_R**0.5)),',IM:',imag(EIG_PRIME_R**0.5)
        WRITE(FID,*) "BETA/ALPHA   = RE:",abs(real(BETA)),',IM:',imag(BETA)
        WRITE(FID,*) ''
        close(FID)

    ENDIF

! DEALLOCATE 
!=======================================================================     
10  IF(ALLOCATED(M_mat_1_0)) DEALLOCATE(M_mat_1_0)
    IF(ALLOCATED(M_mat_1_1)) DEALLOCATE(M_mat_1_1)
    IF(ALLOCATED(M_mat_2_0)) DEALLOCATE(M_mat_2_0)
    IF(ALLOCATED(M_mat_2_1)) DEALLOCATE(M_mat_2_1)
    ! ======================================================================
    IF(ALLOCATED(RUR_BARC)) DEALLOCATE(RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC)
    IF(ALLOCATED(RUR_1)) DEALLOCATE(RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1)
    IF(ALLOCATED(RUR_2)) DEALLOCATE(RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2)
    IF(ALLOCATED(VECR_B1)) DEALLOCATE(VECR_B1,VECR_B2,VECL_1,VECR_1,VECL_2,VECR_2)
    IF(ALLOCATED(EIG_R_mat_1_1)) DEALLOCATE(EIG_R_mat_1_1)
    IF(ALLOCATED(EIG_L_mat_1_1)) DEALLOCATE(EIG_L_mat_1_1)
    IF(ALLOCATED(EIG_R_mat_2_1)) DEALLOCATE(EIG_R_mat_2_1)
    IF(ALLOCATED(EIG_L_mat_2_1)) DEALLOCATE(EIG_L_mat_2_1)
    ! ======================================================================
    IF(ALLOCATED(K_tune)) DEALLOCATE(K_tune)
    IF(ALLOCATED(K_tune_JJ)) DEALLOCATE(K_tune_JJ)
    IF(ALLOCATED(tol_diff_JJ)) DEALLOCATE(tol_diff_JJ)
    IF(ALLOCATED(M_eig_1_0)) DEALLOCATE(M_eig_1_0)
    IF(ALLOCATED(M_eig_1_1)) DEALLOCATE(M_eig_1_1)
    IF(ALLOCATED(M_eig_2_0)) DEALLOCATE(M_eig_2_0)
    IF(ALLOCATED(M_eig_2_1)) DEALLOCATE(M_eig_2_1)
    IF(ALLOCATED(EIG_R_mat_1_0)) DEALLOCATE(EIG_R_mat_1_0)
    IF(ALLOCATED(EIG_L_mat_1_0)) DEALLOCATE(EIG_L_mat_1_0)
    IF(ALLOCATED(EIG_R_mat_2_0)) DEALLOCATE(EIG_R_mat_2_0)
    IF(ALLOCATED(EIG_L_mat_2_0)) DEALLOCATE(EIG_L_mat_2_0)
    IF(ALLOCATED(M_mat_diff_1)) DEALLOCATE(M_mat_diff_1)  
    IF(ALLOCATED(M_mat_diff_2)) DEALLOCATE(M_mat_diff_2)
        
ENDDO ! AK_1
ENDDO ! M_1
    IF(ALLOCATED(RUR_BAR)) DEALLOCATE(RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR)
    IF(ALLOCATED(VECL_BAR)) DEALLOCATE(VECL_BAR,VECR_BAR)
ENDDO loop_ind! IND
    IF(ALLOCATED(M_mat)) DEALLOCATE(M_mat)
    IF(ALLOCATED(M_mat)) DEALLOCATE(M_eig)
    IF(ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
    IF(ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)
ENDDO ! AK_BAR
ENDDO ! M_BAR
IF(ALLOCATED(RUR0)) DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)  

IF (MPI_GLB_RANK.EQ.0) THEN    
    WRITE(*,*) ''
    WRITE(*,*) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
ENDIF
CALL MPI_FINALIZE(IERR)

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================
FUNCTION TARGET_DIST(KVAL,MIND1,MIND2,KBAR,SIG)
    ! Calcualte Eig1 + Sig_bar - Eig2
    IMPLICIT NONE
    INTEGER:: MIND1,MIND2,III
    REAL(P8):: KVAL,KVAL2,KBAR,TARGET_DIST
    REAL(P8),DIMENSION(:),ALLOCATABLE:: EIG_DIST_LIST, EIG_DIST_LIST_I
    COMPLEX(P8):: SIG
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: Meig1,Meig2
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: Mmat1,Mmat2

    CALL EIG_MATRIX(MIND1,KVAL,Mmat1,Meig1,comm_grp = newcomm)
    IF (MPI_GLB_RANK.EQ.0) THEN
        WRITE(*,*) ''
        WRITE(*,*) KVAL
    ENDIF

    KVAL2 = KVAL+KBAR
    CALL EIG_MATRIX(MIND2,KVAL2,Mmat2,Meig2,comm_grp = newcomm)

    IF (MPI_GLB_RANK.EQ.0) THEN
        ALLOCATE(EIG_DIST_LIST(SIZE(Meig1)))
        ALLOCATE(EIG_DIST_LIST_I(SIZE(Meig2)))

        DO III = 1,SIZE(Meig1)
            EIG_DIST_LIST_I = AIMAG(Meig1(III)+SIG-Meig2)
            EIG_DIST_LIST(III) = EIG_DIST_LIST_I(MINLOC(ABS(EIG_DIST_LIST_I),1)) ! CLOSEST DISTANCE for III
        ENDDO
        TARGET_DIST = EIG_DIST_LIST(MINLOC(ABS(EIG_DIST_LIST),1)) ! CLOSEST DISTANCE overall
        WRITE(*,*) TARGET_DIST

        DEALLOCATE(Meig1,Meig2,Mmat1,Mmat2,EIG_DIST_LIST,EIG_DIST_LIST_I)
    ENDIF
    CALL MPI_BCAST(TARGET_DIST,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

END FUNCTION TARGET_DIST
!=======================================================================
SUBROUTINE ROOT_RANGE(LB,UB,UMAX,STEP,MIND1,MIND2,KBAR,SIG,LMIN)
!CALL ROOT_RANGE(K_VAL_1_1,K_VAL_1_0,20.D0,2.D0,M_IND1,M_IND2,K_BAR,EIG_BAR)
    IMPLICIT NONE
    REAL(P8),INTENT(INOUT):: LB,UB
    REAL(P8),OPTIONAL:: LMIN
    
    INTEGER:: III,MIND1,MIND2
    COMPLEX(P8):: SIG
    REAL(P8):: STEP,UMAX,LB_DIST,UB_DIST,KBAR

    IF (.NOT.(PRESENT(LMIN))) LMIN = 2.0 ! Starting lower-bound for K
    LB_DIST = TARGET_DIST(LMIN,MIND1,MIND2,KBAR,SIG)

    DO III = 1, CEILING((UMAX-LMIN)/STEP)
        UB = LMIN + STEP*III
        LB = LMIN + STEP*(III-1)
        UB_DIST = TARGET_DIST(UB,MIND1,MIND2,KBAR,SIG)
        IF (LB_DIST*UB_DIST.LT.0.D0) THEN
            IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) 'ROOT_RANGE completed'
            GOTO 211
        ELSEIF (UB.GE.UMAX) THEN
            IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) 'ROOT_RANGE: CANT FIND ZERO'
            CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        ENDIF
    ENDDO

211 RETURN
END SUBROUTINE ROOT_RANGE   
!=======================================================================

END PROGRAM EVP_ALMOST_DISCRETE_SCAN5
!=======================================================================
