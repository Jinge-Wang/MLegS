PROGRAM EVP_ALMOST_DISCRETE
!=======================================================================
! [USAGE]: 
! USE VPROD_PFF RATHER THAN VPROD
! 
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
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
    INTEGER     :: N_mat1,N_mat2,II,JJ,KK
    INTEGER     :: M_IND1,M_IND2
    INTEGER     :: M_BAR,M_1,M_2,EIG_BAR_IND
    REAL(P8)    :: K_BAR,K_VAL_1_0,K_VAL_1_1,K_VAL_2_0,K_VAL_2_1
    REAL(P8)    :: K_shift,tol,tol_diff,K_delta, LAMB_CUST
    REAL(P8)    :: EIG_LMR_1, ZLEN_BAR, EIG_PRIME_R_L, EIG_PRIME_R_U, EIG_PRIME_R_custom
    REAL(P8)    :: LAMB_LOW, LAMB_UP
    REAL(P8),DIMENSION(:),ALLOCATABLE:: LAMB_1, LAMB_2
    COMPLEX(P8) :: EIG_BAR, EIG_PRIME, EIG_PRIME_R, ALPHA, BETA
    REAL(P8),DIMENSION(:),ALLOCATABLE:: K_tune, K_tune_JJ, tol_diff_JJ
    REAL(P8),DIMENSION(:),ALLOCATABLE:: RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_1_0, M_eig_1_1
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_2_0, M_eig_2_1
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_diff_1,M_mat_diff_2
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_2_0, EIG_R_mat_2_0, EIG_L_mat_2_0
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_BAR,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_1, VECR_B1, RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_2, VECR_B2, RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2
    ! MPI:
    INTEGER:: MPI_GLB_PROCS, MPI_GLB_RANK, newcomm
    ! ADD PERTURB:
    INTEGER:: IS,NDIM
    REAL(P8):: rp1,ip1
    ! ZERO-CROSSING LIMIT:
    INTEGER,PARAMETER:: CROSS_LIMIT = 4
    ! SCANNING:
    CHARACTER(LEN=72):: FILENAME
    INTEGER:: FID


    ! CALL MPI_INIT(IERR)
    CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
    IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
        WRITE(*,*) 'The threading support is lesser than that demanded.'
        CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    ENDIF

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
    IF (MPI_GLB_RANK.EQ.0) THEN
        WRITE(*,*) 'PROGRAM STARTED'
        CALL PRINT_REAL_TIME()   ! @ MOD_MISC
    ENDIF
    CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS
    CALL READCOM('NOECHO')
    CALL READIN(5)


! {M_bar, AK_bar}: BASE FLOW
!=======================================================================
    IF (MPI_GLB_RANK.EQ.0) THEN
    !     WRITE(*,*) ''        
    !     WRITE(*,*) 'STEP 1: SELECT {M_bar,AK_bar}'
    !     WRITE(*,*) ''
    !     WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M >= 0) ; M ='
    !     READ(5,'(I72)') M_BAR

    !     WRITE(*,*) ''
    !     WRITE(*,103) ZLEN
    ! 103  FORMAT(' AK = 2*PI/', F0.3)    
    !     WRITE(*,*) 'CHOOSE AN AXIAL WAVENUMBER (AK) INDEX; AK ='
    !     READ(5,'(F5.3)') K_BAR

    !     WRITE(*,*) ''
    !     WRITE(*,*) 'CALCULATING EIGENVALUES FOR {M_bar,AK_bar}. MIGHT TAKE MINUTES' 
        
        FILENAME = 'almostdegen_read.input'
        CALL READIN(FILENAME,M_BAR,K_BAR,EIG_BAR_IND,M_1,K_VAL_1_0)
    ENDIF

    CALL EIG_MATRIX(M_BAR,K_BAR,M_mat,M_eig,EIG_R_mat,EIG_L_mat,comm_grp=newcomm) ! calculate operator matrix and its corresponding eigenvalues

    ! Pick a target Sigma_bar
    IF (MPI_GLB_RANK.EQ.0) THEN
        ! WRITE(*,*) ''
        ! WRITE(*,*) 'EIGVALS FOR {M_bar,K_bar}:'
        ! CALL MCAT(M_eig) ! print all eigenvalues
        ! WRITE(*,*) 'CHOOSE THE TARGET EIGVAL #:'
        ! READ(5,'(I72)') EIG_BAR_IND
        EIG_BAR = M_eig(EIG_BAR_IND)
        WRITE(*,*) 'TARGET EIG_BAR:', EIG_BAR
        DEALLOCATE(M_eig)
        ! =================== CHECK IF CRITICAL LAYER MODE =====================
        IF (M_BAR.LE.3) THEN
            IF (IMAG(EIG_BAR).GT.-(M_BAR*QPAIR%Q(1) + K_BAR*QPAIR%H(1))) CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ELSE
            IF (IMAG(EIG_BAR).GT.-(M_BAR*QPAIR%Q(1) + K_BAR*QPAIR%H(1))*1.005) CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ENDIF
        ! ======================================================================

        ! Save rur_bar, rup_bar, uz_bar, ror_bar, rop_bar, oz_bar for DEGENERATE PERTURBATION THEORY
        CALL EIG2VELVOR(M_BAR,K_BAR,EIG_BAR_IND,EIG_R_mat,EIG_L_mat,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECL_BAR,comm_grp=newcomm)
        ! WRITE(*,*) size(RUR_BAR)
        ! WRITE(*,*) size(TFM%R)
        ! WRITE(*,*) size(TFM%PF,1),size(TFM%PF,2)
        ! ! WRITE(*,*) RUR_BAR
        ! ! WRITE(*,*) ''
        ! ! DO II = 1,NR
        ! !     RUR_BAR(II) = RUR_BAR(II)/TFM%R(II)
        ! ! ENDDO
        ! ! WRITE(*,*) RUR_BAR
        ! CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)


! ! {M_1, AK_1} & {M_2}
! !=======================================================================
!         ! {M',AK'}: M' = M + M_bar
!         WRITE(*,*) ''
!         WRITE(*,*) 'STEP 2: SELECT M1'
!         WRITE(*,*) ''
!         WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M >= 0) ; M ='
!         READ(5,'(I72)') M_1
        
    ENDIF
    DEALLOCATE(M_mat)
    CALL MPI_BCAST(EIG_BAR,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(M_1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    M_2 = M_1 + M_BAR

! ! TUNE K
! !=======================================================================
!     ! CALL ROOT_RANGE(K_VAL_1_1,K_VAL_1_0,20.D0,1.0,M_1,M_2,K_BAR,EIG_BAR,1.0)
!     IF (MPI_GLB_RANK.EQ.0) THEN
!         WRITE(*,*) 'CHOOSE K_1 GUESS:'
!         READ(5,'(F5.3)') K_VAL_1_0
!     ENDIF
    CALL MPI_BCAST(K_VAL_1_0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

    ! K_VAL_1_0 and K_VAL_1_1 are both for M1,AK1
    ! K_VAL_1_0 = 8.43 ! lower bound
    K_delta = 0.001
    tol = 1.0E-13  ! tolerance

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
        IF(.NOT.ALLOCATED(K_tune)) ALLOCATE(K_tune(CROSS_LIMIT)) !N_mat1)) !*2))
        IF(.NOT.ALLOCATED(K_tune_JJ)) ALLOCATE(K_tune_JJ(CROSS_LIMIT)) !N_mat2)) !*2))

        !$OMP parallel do private(EIG_LMR_1,KK,K_tune_JJ)
        DO JJ = 1,CROSS_LIMIT !N_mat1!2*N_mat1
            ! eigc_b = -(m_bar*Q+k_bar*H);
            IF ((M_1.LE.3).AND.(IMAG(M_eig_1_0(JJ)).GT.-(M_1*QPAIR%Q(1) + K_VAL_1_0*QPAIR%H(1)))) THEN
                K_tune(JJ) = 100.D0
            ELSEIF ((M_1.GT.3).AND.(IMAG(M_eig_1_0(JJ)).GT.-(M_1*QPAIR%Q(1) + K_VAL_1_0*QPAIR%H(1))*1.005)) THEN
                K_tune(JJ) = 100.D0
            ELSE
                EIG_LMR_1 = 0.D0
                EIG_LMR_1 = &
                & AIMAG(DOT_PRODUCT(EIG_L_mat_1_0(:,JJ),M_mat_diff_1(:,JJ))/DOT_PRODUCT(EIG_L_mat_1_0(:,JJ),EIG_R_mat_1_0(:,JJ)))
    
                K_tune_JJ = 0.D0
                DO KK = 1,CROSS_LIMIT !N_mat2!2*N_mat2
                    
                    IF (ABS(REAL(M_eig_2_0(KK) - M_eig_1_0(JJ) - EIG_BAR)) .GT. 1.0E-8) THEN ! Arbitrary threshold
                        K_tune_JJ(KK) = 100.D0
                    ELSEIF ((M_2.LE.3).AND.(IMAG(M_eig_2_0(KK)).GT.-(M_2*QPAIR%Q(1) + K_VAL_2_0*QPAIR%H(1)))) THEN
                        K_tune_JJ(KK) = 100.D0
                    ELSEIF ((M_2.GT.3).AND.(IMAG(M_eig_2_0(KK)).GT.-(M_2*QPAIR%Q(1) + K_VAL_2_0*QPAIR%H(1))*1.005)) THEN
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

        WRITE(*,*) 'K_shift:'
        WRITE(*,*) K_shift
        ! WRITE(*,*) 'Original diff:'
        ! WRITE(*,*) MINVAL(ABS(AIMAG(M_eig_1_0-EIG_BAR))/ABS(AIMAG(EIG_BAR)))
        IF (ABS(K_shift).GT.1) CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    ENDIF

    ! Calculate the actual 
        K_VAL_1_1 = K_VAL_1_0 + K_shift
        CALL EIG_MATRIX(M_1,K_VAL_1_1,M_mat_1_1,M_eig_1_1,comm_grp=newcomm)

        K_VAL_2_1 = K_VAL_1_1 + K_BAR
        CALL EIG_MATRIX(M_2,K_VAL_2_1,M_mat_2_1,M_eig_2_1,comm_grp=newcomm)

    ! Check if tolerance is reached
    ! Match the IMAG part
    IF (MPI_GLB_RANK.EQ.0) THEN
        IF(.NOT.ALLOCATED(tol_diff_JJ)) ALLOCATE(tol_diff_JJ(CROSS_LIMIT)) !N_mat1))!*2))

        DO JJ = 1,CROSS_LIMIT !N_mat1 !2*N_mat1
            tol_diff_JJ(JJ) = MINVAL(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)))
        ENDDO

        tol_diff = MINVAL(ABS(tol_diff_JJ))
        tol_diff = tol_diff/ABS(AIMAG(EIG_BAR))
        WRITE(*,*) II
        WRITE(*,*) tol_diff
        IF (tol_diff.LT.tol) THEN
            WRITE(*,*) 'Desired tol reached'
            WRITE(*,*) 'tol_diff:'
            WRITE(*,*) tol_diff
            WRITE(*,*) 'Number of steps:'
            WRITE(*,*) II
            tol = 999.0
            ! GOTO 20
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
            WRITE(*,*) 'tol_diff:'
            WRITE(*,*) tol_diff
        ENDIF
        CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    ENDIF

! Saving and outputing the results in FFF Space
!=======================================================================   
    ! Find JJth Eigvector of {M1,K1} & KKth Eigvector of {M2,K2} 
20  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

CALL EIG_MATRIX(M_1,K_VAL_1_1,M_mat_1_1,M_eig_1_1,EIG_R_mat_1_1,EIG_L_mat_1_1,comm_grp=newcomm)
CALL EIG_MATRIX(M_2,K_VAL_2_1,M_mat_2_1,M_eig_2_1,EIG_R_mat_2_1,EIG_L_mat_2_1,comm_grp=newcomm)

    IF (MPI_GLB_RANK.EQ.0) THEN
        DO JJ = 1,CROSS_LIMIT !N_mat1 !2*N_mat1
            tol_diff_JJ(JJ) = MINVAL(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)))
            K_tune(JJ) = MINLOC(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)),1)
        ENDDO
        
        JJ = MINLOC(ABS(tol_diff_JJ),1)
        KK = K_tune(JJ)
        !CALL SORT(0,M_eig_1_1)
        ! CALL MSAVE(M_eig_1_1,'M_eig_1_1.dat')
        !CALL SORT(0,M_eig_2_1)
        ! CALL MSAVE(M_eig_2_1,'M_eig_2_1.dat')
        WRITE(*,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
        WRITE(*,201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
        WRITE(*,201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)
201  FORMAT('Eig: M# = ',I3,'; AK = ',F20.13,'; Sig (#',I3') =',1P8E20.13)
    ENDIF

! NEARLY DEGEN
K_VAL_1_1 = real(NINT(K_VAL_1_1*4))/4
K_VAL_2_1 = K_VAL_1_1 + K_BAR
CALL EIG_MATRIX(M_1,K_VAL_1_1,M_mat_1_1,M_eig_1_1,EIG_R_mat_1_1,EIG_L_mat_1_1,comm_grp=newcomm)
CALL EIG_MATRIX(M_2,K_VAL_2_1,M_mat_2_1,M_eig_2_1,EIG_R_mat_2_1,EIG_L_mat_2_1,comm_grp=newcomm)

IF (MPI_GLB_RANK.EQ.0) THEN
! ==================== NEARLY DEGEN PERTURB THEORY =====================
    WRITE(*,*) 'NEARLY DEGEN:'
    WRITE(*,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
    WRITE(*,201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
    WRITE(*,201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)
!=======================================================================   

! Calculate the FINAL RESULT:
! M_bar, K_bar: RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECL_BAR

! {M1,K1}
CALL EIG2VELVOR(M_1,K_VAL_1_1,JJ,EIG_R_mat_1_1,EIG_L_mat_1_1,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECL_1,comm_grp=newcomm)

! {M2,K2}
CALL EIG2VELVOR(M_2,K_VAL_2_1,KK,EIG_R_mat_2_1,EIG_L_mat_2_1,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECL_2,comm_grp=newcomm)
CALL NONLIN_MK(M_2,K_VAL_2_1,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_B1,comm_grp=newcomm)
! CALL MCAT(UZ_2)

! {M1,K1}
ALLOCATE(RUR_BARC(NDIMR),RUP_BARC(NDIMR),UZ_BARC(NDIMR),ROR_BARC(NDIMR),ROP_BARC(NDIMR),OZ_BARC(NDIMR))
RUR_BARC = CONJG(RUR_BAR)
RUP_BARC = CONJG(RUP_BAR)
 UZ_BARC = CONJG( UZ_BAR)
ROR_BARC = CONJG(ROR_BAR)
ROP_BARC = CONJG(ROP_BAR)
 OZ_BARC = CONJG( OZ_BAR)
CALL NONLIN_MK(M_1,K_VAL_1_1,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_B2,comm_grp=newcomm)
! call mcat(VECR_B2)

! Lambda range
LAMB_LOW =  MAX(abs(M_eig_2_1(KK) - M_eig_1_1(JJ) - EIG_BAR)/abs(DOT_PRODUCT(VECL_2,VECR_B1)), &
                abs(M_eig_2_1(KK) - M_eig_1_1(JJ) - EIG_BAR)/abs(DOT_PRODUCT(VECL_1,VECR_B2)) )

ALLOCATE(LAMB_2(SIZE(M_eig_2_1)))
FORALL(II=1:SIZE(M_eig_2_1),II/=KK) LAMB_2(II)=abs(M_eig_2_1(II)-M_eig_1_1(JJ)-EIG_BAR)/abs(DOT_PRODUCT(EIG_L_mat_2_1(:,II),VECR_B1))
LAMB_2(KK) = 999
LAMB_UP = MINVAL(LAMB_2)
!CALL MCAT(LAMB_2)
DEALLOCATE(LAMB_2)

ALLOCATE(LAMB_1(SIZE(M_eig_1_1)))
FORALL(II=1:SIZE(M_eig_1_1),II/=JJ) LAMB_1(II)=abs(M_eig_2_1(KK)-M_eig_1_1(II)-EIG_BAR)/abs(DOT_PRODUCT(EIG_L_mat_1_1(:,II),VECR_B2))
LAMB_1(JJ) = 999
LAMB_UP = MIN(LAMB_UP,MINVAL(LAMB_1))
!CALL MCAT(LAMB_1)
DEALLOCATE(LAMB_1)

WRITE(*,*) 'LAMBDA RANGE: '
WRITE(*,*) 'MIN: ', LAMB_LOW
WRITE(*,*) 'MAX: ', LAMB_UP

! Sigma(1):
EIG_PRIME_R = DOT_PRODUCT(VECL_1,VECR_B2)*DOT_PRODUCT(VECL_2,VECR_B1)
! WRITE(*,*) "L1B*R2-L2BR1 = ", DOT_PRODUCT(VECL_1,VECR_B2)-DOT_PRODUCT(VECL_2,VECR_B1)
EIG_PRIME_R_L = (0.25*(EIG_BAR+M_eig_1_1(JJ)-M_eig_2_1(KK))**2+LAMB_LOW**2*EIG_PRIME_R)**0.5
EIG_PRIME_R_U = (0.25*(EIG_BAR+M_eig_1_1(JJ)-M_eig_2_1(KK))**2+ LAMB_UP**2*EIG_PRIME_R)**0.5
! WRITE(*,*) 'Customized Lambda value: '
! READ(5,'(F8.6)') LAMB_CUST
LAMB_CUST = 0.01
EIG_PRIME_R_custom = (0.25*(EIG_BAR+M_eig_1_1(JJ)-M_eig_2_1(KK))**2+LAMB_CUST**2*EIG_PRIME_R)**0.5
WRITE(*,*) ''
WRITE(*,*) "Lambda*SIGMA'(nearly degen)= +/- ",EIG_PRIME_R_L,'~',EIG_PRIME_R_U
WRITE(*,*) "Lambda*SIGMA'(nearly degen)= +/- ",EIG_PRIME_R_custom,' (',LAMB_CUST,')'
WRITE(*,*) " SIGMA(degen)'^2= ",EIG_PRIME_R
WRITE(*,*) " SIGMA(degen)'  = ",EIG_PRIME_R**0.5
! Beta/Alpha: ALPHA * {M1,K1} + BETA * {M2,K2}
ALPHA = 1.D0
BETA = EIG_PRIME_R**0.5/DOT_PRODUCT(VECL_1,VECR_B2)
WRITE(*,*) " BETA(degen)    = ",BETA
! BETA = (DOT_PRODUCT(VECL_2,VECR_B1)/DOT_PRODUCT(VECL_1,VECR_B2))**0.5
EIG_R_mat_1_1(:,JJ) = EIG_R_mat_1_1(:,JJ)*ALPHA
EIG_R_mat_2_1(:,KK) = EIG_R_mat_2_1(:,KK)*BETA
! WRITE(*,*) " BETA(degen)    = ",BETA

! ======================================================================
    IF (abs(REAL(EIG_PRIME_R**0.5)).GT.0.01) THEN
        open(unit=FID,FILE='./converg/perturb_scan.GROW',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
        WRITE(FID,*) ''
        WRITE(FID,*) 'NEARLY DEGEN: *GROWTH*'
    ELSE
        open(unit=FID,FILE='./converg/perturb_scan.NEGLIGIBLE',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
        WRITE(FID,*) ''
        WRITE(FID,*) 'NEARLY DEGEN: *NEGLIGIBLE*'
    ENDIF
    WRITE(FID,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
    WRITE(FID,201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
    WRITE(FID,201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)
    WRITE(FID,*) 'LAMBDA RANGE: '
    WRITE(FID,*) 'MIN: ', LAMB_LOW,'; MAX: ', LAMB_UP
    WRITE(FID,*) "Lambda*SIGMA'(nearly degen)= +/- ",EIG_PRIME_R_L,'~',EIG_PRIME_R_U
    WRITE(FID,*) "Lambda*SIGMA'(nearly degen)= +/- ",EIG_PRIME_R_custom,' (',LAMB_CUST,')'
    WRITE(FID,*) "SIGMA(degen)'= RE:",abs(real(EIG_PRIME_R**0.5)),',IM:',imag(EIG_PRIME_R**0.5)
    WRITE(FID,*) "BETA/ALPHA   = RE:",abs(real(BETA)),',IM:',imag(BETA)
    WRITE(FID,*) ''
    close(FID)
! ======================================================================


! CREATE PERTURBATION FILE
! ======================================================================
open(10,FILE='./converg/new_perturb.input',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)
if (IS.ne.0) then
    print *, 'ERROR: createperturb -- Could not creat new file new_perturb.dat'
    GOTO 1
end if

! Number of perturbation 
II = 3
WRITE(10,'(I2)') II

! Scale
WRITE(6,*) ''
WRITE(6,*) 'SELECT SCALE(EPSILON)'
WRITE(6,*) 'REAL(SCALE) ='
READ(5,*) rp1
WRITE(6,*) 'IMAG(SCALE) ='
READ(5,*) ip1
WRITE(10,'(F9.7,1X,F9.7)') rp1, ip1
WRITE(10,'(F9.7,1X,F9.7)') rp1, ip1
WRITE(10,'(F9.7,1X,F9.7)') rp1, ip1
            
! WRITE(*,*) 'NEARLY DEGEN:'
! WRITE(*,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
! WRITE(*,201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
! WRITE(*,201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)

! EIB_BAR
WRITE(10,*) 'right eigenfunction'
WRITE(10,108) 'ndim','m','k'
NDIM = size(EIG_R_mat,1)/2
WRITE(10,109) NDIM, M_BAR, K_BAR
    ! Below are not actually used:
    WRITE(10,108) 'q','h','b'
    WRITE(10,110) QPAIR%Q(1),QPAIR%H(1),QPAIR%B(1) ! MIGHT NEED TO CHANGE!!!!!! 
    WRITE(10,108) 'nu_pow','nu','lmap'
    WRITE(10,110) VISC%NUP,VISC%NU,ELL
    WRITE(10,'(9X,A9,9X,A9)') 're(sigma)','im(sigma)'
    WRITE(10,113) REAL(EIG_BAR),AIMAG(EIG_BAR)
    WRITE(10,*) 'radial spectral coefficients'
    WRITE(10,111) 're(psi)','im(psi)','re(chi)','im(chi)'
DO II = 1,NDIM
    WRITE(10,112) REAL(EIG_R_mat(II,EIG_BAR_IND)),AIMAG(EIG_R_mat(II,EIG_BAR_IND)),REAL(EIG_R_mat(NDIM+II,EIG_BAR_IND)),AIMAG(EIG_R_mat(NDIM+II,EIG_BAR_IND))
ENDDO
! EIG_1
WRITE(10,*) 'right eigenfunction'
WRITE(10,108) 'ndim','m','k'
NDIM = size(EIG_R_mat_1_1,1)/2
WRITE(10,109) NDIM, M_1, K_VAL_1_1
    ! Below are not actually used:
    WRITE(10,108) 'q','h','b'
    WRITE(10,110) QPAIR%Q(1),QPAIR%H(1),QPAIR%B(1) ! MIGHT NEED TO CHANGE!!!!!! 
    WRITE(10,108) 'nu_pow','nu','lmap'
    WRITE(10,110) VISC%NUP,VISC%NU,ELL
    WRITE(10,'(9X,A9,9X,A9)') 're(sigma)','im(sigma)'
    WRITE(10,113) REAL(M_eig_1_1(JJ)),AIMAG(M_eig_1_1(JJ))
    WRITE(10,*) 'radial spectral coefficients'
    WRITE(10,111) 're(psi)','im(psi)','re(chi)','im(chi)'
DO II = 1,NDIM
    WRITE(10,112) REAL(EIG_R_mat_1_1(II,JJ)),AIMAG(EIG_R_mat_1_1(II,JJ)),REAL(EIG_R_mat_1_1(NDIM+II,JJ)),AIMAG(EIG_R_mat_1_1(NDIM+II,JJ))
ENDDO
! EIG_2
WRITE(10,*) 'right eigenfunction'
WRITE(10,108) 'ndim','m','k'
NDIM = size(EIG_R_mat_2_1,1)/2
WRITE(10,109) NDIM, M_2, K_VAL_2_1
    ! Below are not actually used:
    WRITE(10,108) 'q','h','b'
    WRITE(10,110) QPAIR%Q(1),QPAIR%H(1),QPAIR%B(1) ! MIGHT NEED TO CHANGE!!!!!! 
    WRITE(10,108) 'nu_pow','nu','lmap'
    WRITE(10,110) VISC%NUP,VISC%NU,ELL
    WRITE(10,'(9X,A9,9X,A9)') 're(sigma)','im(sigma)'
    WRITE(10,113) REAL(M_eig_2_1(KK)),AIMAG(M_eig_2_1(KK))
    WRITE(10,*) 'radial spectral coefficients'
    WRITE(10,111) 're(psi)','im(psi)','re(chi)','im(chi)'
DO II = 1,NDIM
    WRITE(10,112) REAL(EIG_R_mat_2_1(II,KK)),AIMAG(EIG_R_mat_2_1(II,KK)),REAL(EIG_R_mat_2_1(NDIM+II,KK)),AIMAG(EIG_R_mat_2_1(NDIM+II,KK))
ENDDO

108  FORMAT(3(9X,A9))
109  FORMAT(2(9X,I9),9X,F9.5)
110  FORMAT(3(9X,F9.5))
111  FORMAT(4(9X,A9))
112  FORMAT(4(G20.12,:,''))
113  FORMAT(9X,F9.5,9X,F9.5)

close(10)
! ======================================================================

! ANALYSIS
! 1. SAVE UR, UP, UZ
call SAVE_VEL(RUR_BAR,RUP_BAR,UZ_BAR,'./converg/vel_MK_'//ITOA3(M_BAR)//'_'//ITOA4(INT(K_BAR*100))&
//'_IND_'//ITOA3(EIG_BAR_IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
call SAVE_VEL(RUR_1,RUP_1,UZ_1,'./converg/vel_MK_'//ITOA3(M_1)//'_'//ITOA4(INT(K_VAL_1_1*100))&
//'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
call SAVE_VEL(RUR_2,RUP_2,UZ_2,'./converg/vel_MK_'//ITOA3(M_2)//'_'//ITOA4(INT(K_VAL_2_1*100))&
//'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')


! DEALLOCATE 
!======================================================================= 
1    DEALLOCATE(RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR)
DEALLOCATE(RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC)
DEALLOCATE(RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1)
DEALLOCATE(RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2)
DEALLOCATE(VECR_B1,VECR_B2,VECL_1,VECL_2,VECL_BAR)
DEALLOCATE(EIG_R_mat_1_1,EIG_L_mat_1_1,EIG_R_mat_2_1,EIG_L_mat_2_1)

    ENDIF

! DEALLOCATE 
!=======================================================================     
10  DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
    DEALLOCATE(M_mat_1_0, M_mat_1_1, M_mat_2_0, M_mat_2_1)

    IF (MPI_GLB_RANK.EQ.0) THEN
        DEALLOCATE(EIG_R_mat,EIG_L_mat)
        DEALLOCATE(K_tune, K_tune_JJ, tol_diff_JJ)
        DEALLOCATE(M_eig_1_0, M_eig_1_1)
        DEALLOCATE(M_eig_2_0, M_eig_2_1)
        DEALLOCATE(EIG_R_mat_1_0, EIG_L_mat_1_0)
        DEALLOCATE(EIG_R_mat_2_0, EIG_L_mat_2_0)
        DEALLOCATE(M_mat_diff_1)  
        DEALLOCATE(M_mat_diff_2)
        
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
SUBROUTINE EIG_MATRIX(MREAD,AKREAD,H,EIG_VAL,EIG_VEC_R,EIG_VEC_L,comm_grp)
    IMPLICIT NONE
    INTEGER,INTENT(IN)    :: MREAD, comm_grp
    REAL(P8),INTENT(INOUT)    :: AKREAD
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: H
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT),OPTIONAL:: EIG_VAL
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT),OPTIONAL:: EIG_VEC_R, EIG_VEC_L

    INTEGER     :: NR_MK, I
    REAL(P8)    :: REI
    INTEGER     :: MPI_NR_SIZE,MPI_NR_INDEX
    ! REAL(P8),DIMENSION(:),ALLOCATABLE:: RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: RURU,RUPU,UZU,RORU,ROPU,OZU
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSIU,CHIU,PSI1,CHI1,PSI2,CHI2,PSI3,CHI3

    ZLEN = 2*PI/AKREAD
    CALL MPI_BCAST(MREAD,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(ZLEN,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    
    ! INIT 
    NTH = 2
    NX = 4
    NTCHOP = 2 ! M = 0, M = MREAD
    NXCHOP = 2 ! K = 0, K = AKREAD
    CALL LEGINIT(comm_grp,MREAD)
    IF (.NOT.ALLOCATED(RUR0)) THEN
        ALLOCATE(RUR0(NR),RUP0(NR),UZ0(NR),ROR0(NR),ROP0(NR),OZ0(NR))
        CALL INIT_LOOP(RUR0,RUP0,UZ0,ROR0,ROP0,OZ0) ! ALL IN PFF space
    ENDIF

    ! SPLIT NRCHOPDIM
    NR_MK = NRCHOPS(2)
    IF (ALLOCATED(H).AND.(SIZE(H,1).NE.NR_MK)) THEN
        DEALLOCATE(H)
        ALLOCATE(H(2*NR_MK,2*NR_MK))
    ELSEIF (.NOT.(ALLOCATED(H))) THEN
        ALLOCATE(H(2*NR_MK,2*NR_MK))
    ENDIF 
    H = CMPLX(0.D0,0.D0,P8)
    CALL DECOMPOSE(NR_MK, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)

    ! MAIN JOB
    ALLOCATE(PSIU(NRCHOPDIM),CHIU(NRCHOPDIM))
    ALLOCATE(PSI1(NRCHOPDIM),CHI1(NRCHOPDIM))
    ALLOCATE(PSI2(NRCHOPDIM),CHI2(NRCHOPDIM))
    ALLOCATE(PSI3(NRCHOPDIM),CHI3(NRCHOPDIM))
    ALLOCATE(RURU(NRCHOPDIM),RUPU(NRCHOPDIM),UZU(NRCHOPDIM))
    ALLOCATE(RORU(NRCHOPDIM),ROPU(NRCHOPDIM),OZU(NRCHOPDIM))

    DO I = MPI_NR_INDEX+1,MPI_NR_INDEX+MPI_NR_SIZE
    ! ================================ PSI =================================
        PSIU = 0.D0
        CHIU = 0.D0
        PSIU(I) = 1.D0

        ! VISCOSITY/HYPERV
        CALL CHOPSET(2)
        IF (VISC%SW .EQ. 1) THEN
        CALL DEL2_MK(PSIU,PSI3)
        PSI3 = PSI3 * VISC%NU
        ELSE IF (VISC%SW .EQ. 2) THEN
        CALL CHOPSET(-2+VISC%P)
        CALL HELMP_MK(VISC%P,PSIU,PSI3,0.D0)
        IF (MOD(VISC%P/2,2).EQ.0) PSI3 = -PSI3
        PSI3 = PSI3 * VISC%NUP
        CALL CHOPSET( 2-VISC%P)
        ELSE
        PSI3 = 0.D0
        ENDIF
        CHI3 = 0.D0
        CALL CHOPSET(-2)

        CALL CHOPSET(3)
        ! NONLINEAR TERM: W0 X U'
        CALL PC2VEL_MK(PSIU,CHIU,RURU,RUPU,UZU)
        CALL RTRAN_MK(RURU,1)
        CALL RTRAN_MK(RUPU,1)
        CALL RTRAN_MK(UZU,1)
        CALL VPROD_MK(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT_MK(RURU,RUPU,UZU,PSI1,CHI1)

        ! NONLINEAR TERM: W' X U0
        CALL PC2VOR_MK(PSIU,CHIU,RORU,ROPU,OZU)
        CALL RTRAN_MK(RORU,1)
        CALL RTRAN_MK(ROPU,1)
        CALL RTRAN_MK(OZU,1)
        CALL VPROD_MK(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT_MK(RORU,ROPU,OZU,PSI2,CHI2)
        CALL CHOPSET(-3)
        ! IF (I.EQ.2) THEN
        !     CALL XXDX_MK(PSIU,RURU)
        !     CALL MCAT(ROP0)
        ! ENDIF

        PSI1 = -PSI1 + PSI2 + PSI3
        CHI1 = -CHI1 + CHI2 + CHI3
        PSI1(NRCHOPS(2)+1:) = 0.D0
        CHI1(NRCHOPS(2)+1:) = 0.D0

        H(:NR_MK,  I) = PSI1(:NR_MK)
        H(NR_MK+1:,I) = CHI1(:NR_MK)

    ! ================================ CHI =================================
        PSIU = 0.D0
        CHIU = 0.D0
        CHIU(I) = 1.D0

        CALL IDEL2_MK(CHIU,CHI1)
        CHIU = CHI1
        CHI1 = 0.D0

        IF (VISC%SW .NE. 0) THEN
        CHI3 = PSI3
        ELSE
        CHI3 = 0.D0
        ENDIF
        PSI3 = 0.D0

        CALL CHOPSET(3)
        ! NONLINEAR TERM: W0 X U'
        CALL PC2VEL_MK(PSIU,CHIU,RURU,RUPU,UZU)
        CALL RTRAN_MK(RURU,1)
        CALL RTRAN_MK(RUPU,1)
        CALL RTRAN_MK(UZU,1)
        CALL VPROD_MK(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT_MK(RURU,RUPU,UZU,PSI1,CHI1)

        ! NONLINEAR TERM: W' X U0
        CALL PC2VOR_MK(PSIU,CHIU,RORU,ROPU,OZU)
        CALL RTRAN_MK(RORU,1)
        CALL RTRAN_MK(ROPU,1)
        CALL RTRAN_MK(OZU,1)
        CALL VPROD_MK(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT_MK(RORU,ROPU,OZU,PSI2,CHI2)
        CALL CHOPSET(-3)

        PSI1 = -PSI1 + PSI2 + PSI3
        CHI1 = -CHI1 + CHI2 + CHI3
        PSI1(NRCHOPS(2)+1:) = 0.D0
        CHI1(NRCHOPS(2)+1:) = 0.D0

        H(:NR_MK,  I+NR_MK) = PSI1(:NR_MK)
        H(NR_MK+1:,I+NR_MK) = CHI1(:NR_MK)
    ENDDO
    ! DEALLOCATE( RUR0,RUP0,UZ0,ROR0,ROP0,OZ0 )
    DEALLOCATE( RUPU,UZU,RORU,ROPU,OZU )
    DEALLOCATE( PSIU,CHIU,PSI1,CHI1,PSI2,CHI2,PSI3,CHI3 )

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,H,SIZE(H),MPI_DOUBLE_COMPLEX,MPI_SUM,&
    MPI_COMM_WORLD,IERR)

! OUTPUT:
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    IF (MPI_GLB_RANK.EQ.0) THEN  ! START OF THE SERIAL PART
        !WRITE OUT THE 3D MATRIX IN THE SCALAR USING UNFORMATTED MODE
        
        IF (PRESENT(EIG_VAL)) THEN
            IF (ALLOCATED(EIG_VAL)) DEALLOCATE(EIG_VAL)
            ALLOCATE( EIG_VAL(2*NR_MK) )

            IF (PRESENT(EIG_VEC_L)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE(EIG_VEC_R) 
                IF (ALLOCATED(EIG_VEC_L)) DEALLOCATE(EIG_VEC_L) 
                ALLOCATE( EIG_VEC_R(2*NR_MK,2*NR_MK),EIG_VEC_L(2*NR_MK,2*NR_MK) )
                CALL EIGENDECOMPOSE(H,EIG_VAL,ER = EIG_VEC_R,EL = EIG_VEC_L)
            ELSEIF (PRESENT(EIG_VEC_R)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE(EIG_VEC_R) 
                ALLOCATE( EIG_VEC_R(2*NR_MK,2*NR_MK) )
                CALL EIGENDECOMPOSE(H,EIG_VAL,ER = EIG_VEC_R)
            ELSE
                CALL EIGENDECOMPOSE(H,EIG_VAL)
            ENDIF
        ENDIF

        ! WRITE(*,*) 'NON-HERMITIAN MATRIX EIGENVALUE SOLVER'
        ! WRITE(*,*) 'CURRENTLY RUNNING IN SERIAL...'
        ! WRITE(*,*) 'EIGVALS:'
        ! CALL MCAT(EIG_VAL)
    ENDIF

    RETURN

END SUBROUTINE EIG_MATRIX
!=======================================================================
! SUBROUTINE EIG2VELVOR(MREAD,AKREAD,EIG_IND,EIG_VEC_R,EIG_VEC_L,RUR,RUP,UZ,ROR,ROP,OZ,VEC_L,comm_grp)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)    :: MREAD, EIG_IND, comm_grp
!     REAL(P8),INTENT(INOUT)    :: AKREAD
!     COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: EIG_VEC_R, EIG_VEC_L
!     COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: RUR,RUP,UZ,ROR,ROP,OZ,VEC_L

!     INTEGER:: NR_MK
!     COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSI,CHI,CHI_TEMP
!     REAL(P8):: ENE_MK

!     REAL(P8):: ANG
!     COMPLEX(P8):: PHASE_FACTOR

!     IF ((M(2).NE.MREAD).OR.(ABS(ZLEN-2*PI/AKREAD).GE.1.0E-13)) THEN
!         WRITE(*,*) 'EIG2VELVOR: SET NEW M AND AK'
!         ZLEN = 2*PI/AKREAD
!         NTH = 2
!         NX = 4
!         NTCHOP = 2 ! M = 0, M = MREAD
!         NXCHOP = 2 ! K = 0, K = AKREAD
!         CALL LEGINIT(comm_grp,MREAD)
!     ENDIF

!     NR_MK = SIZE(EIG_VEC_R,1)/2
!     ALLOCATE(VEC_L(SIZE(EIG_VEC_R,1)))
!     VEC_L(:) = EIG_VEC_L(:,EIG_IND)

!     ALLOCATE(PSI(NRCHOPDIM),CHI(NRCHOPDIM),CHI_TEMP(NRCHOPDIM))  
!     ALLOCATE(RUR(NRCHOPDIM),RUP(NRCHOPDIM),UZ(NRCHOPDIM))
!     ALLOCATE(ROR(NRCHOPDIM),ROP(NRCHOPDIM),OZ(NRCHOPDIM))  

!     PSI = 0.D0
!     CHI = 0.D0
!     PSI(:NR_MK) = EIG_VEC_R(:NR_MK,EIG_IND) !PSI
!     CHI(:NR_MK) = EIG_VEC_R(NR_MK+1:,EIG_IND) !DEL2CHI
!     CALL IDEL2_MK(CHI,CHI_TEMP)

!     ! =========================== NORMALIZATION ============================
!     ! 1). NORMALIZE THE ENERGY =============================================
!     ! ENE_MK = ENERGY_MK(PSI,CHI_TEMP)

!     ! PSI = PSI/ENE_MK**0.5
!     ! CHI_TEMP = CHI_TEMP/ENE_MK**0.5
!     ! VEC_L(:) = VEC_L(:)*ENE_MK**0.5
!     ! CALL CHOPSET(3)
!     ! ======================================================================

!     ! 2). NORMALIZE THE MAX UP =============================================
!     CALL CHOPSET(3)
!     CALL PC2VEL_MK(PSI,CHI_TEMP,RUR,RUP,UZ)
!     ! CALCULATE THE PHASE_FACTOR:
!     CALL RTRAN_MK(RUP,1)
!     ANG = ATAN2(IMAG(RUP(1)),REAL(RUP(1)))
!     PHASE_FACTOR = EXP(-IU*ANG)
!     RUP = RUP*PHASE_FACTOR ! RUP NOW PURELY REAL
!     RUP = RUP/TFM%R ! RUP -> UP
!     PHASE_FACTOR = PHASE_FACTOR/MAXVAL(ABS(RUP(:)))
!     ! APPLY THE PHASE_FACTOR:
!     PSI = PSI*PHASE_FACTOR
!     CHI_TEMP = CHI_TEMP*PHASE_FACTOR
!     VEC_L(:) = VEC_L(:)/PHASE_FACTOR
!     RUR = 0.D0
!     RUP = 0.D0
!      UZ = 0.D0
!     ! ======================================================================

!     CALL PC2VEL_MK(PSI,CHI_TEMP,RUR,RUP,UZ)    
!     CALL PC2VOR_MK(PSI,CHI_TEMP,ROR,ROP,OZ)
!     CALL RTRAN_MK(RUR,1)
!     CALL RTRAN_MK(RUP,1)
!     CALL RTRAN_MK( UZ,1)
!     CALL RTRAN_MK(ROR,1)
!     CALL RTRAN_MK(ROP,1)
!     CALL RTRAN_MK( OZ,1)
!     CALL CHOPSET(-3)

!     ! ========================= ADJUST EIG_VEC_R ===========================
!     ! ================= ONLY CHANGE EIG_VEC_R(:,EIG_IND) ===================
!     EIG_VEC_R(:NR_MK,EIG_IND) = PSI(:NR_MK) !PSI
!     EIG_VEC_R(NR_MK+1:,EIG_IND) = CHI_TEMP(:NR_MK) !CHI
!     ! WRITE(*,*) 'NORMED ENE:', ENERGY_MK(PSI,CHI_TEMP)      

!     ! DEBUG:
!    ! WRITE(*,*) 'MAX:'
!    ! WRITE(*,*) MAXVAL(REAL(RUP/TFM%R))
!    ! WRITE(*,*) MAXVAL(-REAL(RUP/TFM%R))

!     DEALLOCATE(PSI,CHI,CHI_TEMP)

! END SUBROUTINE EIG2VELVOR
! ! ======================================================================
! SUBROUTINE NONLIN_MK(MREAD,AKREAD,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,RUR,RUP,UZ,ROR,ROP,OZ,VEC_R,comm_grp)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)    :: MREAD, comm_grp
!     REAL(P8),INTENT(INOUT):: AKREAD
!     COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(IN):: RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0
!     COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(IN):: RUR,RUP,UZ
!     COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: ROR,ROP,OZ
!     COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: VEC_R

!     INTEGER:: NN,NR_MK
!     REAL(P8):: A1,A2,A3,C1,C2,C3,B1,B2,B3,D1,D2,D3
!     COMPLEX(P8),DIMENSION(NRCHOPDIM):: PSI,CHI

!     IF ((M(2).NE.MREAD).OR.(ABS(ZLEN-2*PI/AKREAD).GE.1.0E-13)) THEN
!         WRITE(*,*) 'NONLIN_MK: SET NEW M AND AK'
!         ZLEN = 2*PI/AKREAD
!         NTH = 2
!         NX = 4
!         NTCHOP = 2 ! M = 0, M = MREAD
!         NXCHOP = 2 ! K = 0, K = AKREAD
!         CALL LEGINIT(comm_grp,MREAD)
!     ENDIF

!     NR_MK = NRCHOPS(2)
!     !$OMP PARALLEL DO DEFAULT(SHARED) &
!     !$OMP& PRIVATE(A1,A2,A3,C1,C2,C3,B1,B2,B3,D1,D2,D3)
!     DO NN=1,NR
!         ! U_BAR X W'
!         A1=REAL(RUR_0(NN))
!         A2=REAL(RUP_0(NN))
!         A3=REAL(UZ_0 (NN))
!         C1=AIMAG(RUR_0(NN))
!         C2=AIMAG(RUP_0(NN))
!         C3=AIMAG(UZ_0 (NN))
!         B1=REAL  (ROR(NN))
!         B2=REAL  (ROP(NN))
!         B3=REAL  (OZ (NN))
!         D1=AIMAG (ROR(NN))
!         D2=AIMAG (ROP(NN))
!         D3=AIMAG (OZ (NN))
!         ROR(NN)=CMPLX(A2*B3-A3*B2+C3*D2-C2*D3,B3*C2-B2*C3+A2*D3-A3*D2,P8)
!         ROP(NN)=CMPLX(A3*B1-A1*B3+C1*D3-C3*D1,B1*C3-B3*C1+A3*D1-A1*D3,P8)
!             OZ(NN)=CMPLX(A1*B2-A2*B1+C2*D1-C1*D2,B2*C1-B1*C2+A1*D2-A2*D1,P8)&
!                     /TFM%R(NN)**2.D0
!         ! ROR(NN)= RUP_0(NN)* OZ(NN)- UZ_0(NN)*ROP(NN)
!         ! ROP(NN)=  UZ_0(NN)*ROR(NN)-RUR_0(NN)* OZ(NN)
!         !  OZ(NN)=(RUR_0(NN)*ROP(NN)-RUP_0(NN)*ROR(NN))/TFM%R(NN)**2.D0
!     ! enddo
!     ! !$OMP END PARALLEL DO
!     ! call mcat(oz)
!     ! !$OMP PARALLEL DO DEFAULT(SHARED) &
!     ! !$OMP& PRIVATE(A1,A2,A3,C1,C2,C3,B1,B2,B3,D1,D2,D3)
!     ! do nn = 1,nr

!         ! W_BAR X U'
!         A1=REAL(ROR_0(NN))
!         A2=REAL(ROP_0(NN))
!         A3=REAL(OZ_0 (NN))
!         C1=AIMAG(ROR_0(NN))
!         C2=AIMAG(ROP_0(NN))
!         C3=AIMAG(OZ_0 (NN))
!         B1=REAL  (RUR(NN))
!         B2=REAL  (RUP(NN))
!         B3=REAL  (UZ (NN))
!         D1=AIMAG (RUR(NN))
!         D2=AIMAG (RUP(NN))
!         D3=AIMAG (UZ (NN))
!         ROR(NN)=ROR(NN)-CMPLX(A2*B3-A3*B2+C3*D2-C2*D3,B3*C2-B2*C3+A2*D3-A3*D2,P8)
!         ROP(NN)=ROP(NN)-CMPLX(A3*B1-A1*B3+C1*D3-C3*D1,B1*C3-B3*C1+A3*D1-A1*D3,P8)
!             OZ(NN)= OZ(NN)-(CMPLX(A1*B2-A2*B1+C2*D1-C1*D2,B2*C1-B1*C2+A1*D2-A2*D1,P8)&
!                     /TFM%R(NN)**2.D0)
!         ! ROR(NN)=ROR(NN) - (ROP_0(NN)* UZ(NN)- OZ_0(NN)*RUP(NN))
!         ! ROP(NN)=ROP(NN) - ( OZ_0(NN)*RUR(NN)-ROR_0(NN)* UZ(NN))
!         !  OZ(NN)= OZ(NN) -((ROR_0(NN)*RUP(NN)-ROP_0(NN)*RUR(NN))/TFM%R(NN)**2.D0)
!     ENDDO
!     !$OMP END PARALLEL DO

!     CALL PROJECT_MK(ROR,ROP,OZ,PSI,CHI)
!     ALLOCATE(VEC_R(2*NR_MK))
!     VEC_R(:NR_MK) = PSI(:NR_MK)
!     VEC_R(NR_MK+1:) = CHI(:NR_MK)

! END SUBROUTINE NONLIN_MK
! ! ======================================================================
! SUBROUTINE SAVE_VEL(RUR,RUP,UZ,VEL_FILENAME)
!     IMPLICIT NONE

!     INTEGER:: VEL_FID=17
!     CHARACTER(LEN=*):: VEL_FILENAME
!     COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(IN):: RUR,RUP,UZ

!     open(VEL_FID,FILE=VEL_FILENAME,STATUS='unknown',ACTION='WRITE',IOSTAT=IS)

!     DO II = 1,NR
!         WRITE(VEL_FID,117) TFM%R(II),REAL(RUR(II))/TFM%R(II),AIMAG(RUR(II))/TFM%R(II) &
!                                     ,REAL(RUP(II))/TFM%R(II),AIMAG(RUP(II))/TFM%R(II) &
!                                     ,REAL(UZ(II)),AIMAG(UZ(II))
!     ENDDO
    
! 117  FORMAT(7(G20.12,:,','))
!     close(VEL_FID)

! END SUBROUTINE SAVE_VEL

END PROGRAM EVP_ALMOST_DISCRETE
!=======================================================================
