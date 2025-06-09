PROGRAM EVP_SCAN_0915
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
    INTEGER     :: M_BAR,M_1,M_2,EIG_BAR_IND,tol_ind(2)
    REAL(P8)    :: K_BAR,K_VAL_1_0,K_VAL_2_0
    REAL(P8)    :: K_shift,tol,tol_diff,K_delta
    REAL(P8)    :: EIG_LMR_1, ZLEN_BAR
    COMPLEX(P8) :: EIG_BAR, EIG_PRIME, EIG_PRIME_R, ALPHA, BETA
    REAL(P8),DIMENSION(:),ALLOCATABLE:: K_tune, K_tune_JJ, tol_diff_JJ
    LOGICAL, DIMENSION(:), ALLOCATABLE:: MASK_1, MASK_2
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_1_0, M_eig_1_1
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_2_0, M_eig_2_1
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_diff_1,M_mat_diff_2
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_2_0, EIG_R_mat_2_0, EIG_L_mat_2_0
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
DO M_BAR = 0,10; K_BAR = -10.10; DO WHILE (K_BAR.LT.10.10)
    
    K_BAR = NINT(10*(K_BAR + 0.10))/10.D0
    IF ((M_BAR.EQ.0).AND.(K_BAR.EQ.0.D0)) CYCLE
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    IF (ALLOCATED(M_mat)) DEALLOCATE(M_mat)
    IF (ALLOCATED(M_eig)) DEALLOCATE(M_eig)
    IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
    IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)    

    ! calculate operator matrix and its corresponding eigenvalues
    IF (.NOT. ALLOCATED(RUR0)) THEN 
        CALL EIG_MATRIX(M_BAR,K_BAR,M_mat,M_eig,EIG_R_mat,EIG_L_mat,comm_grp=newcomm)
    ELSE
        CALL READ_MODE( M_BAR,K_BAR,M_mat,M_eig,EIG_R_mat,EIG_L_mat )
    ENDIF

! ======================================================================
loop_ind: DO EIG_BAR_IND = 1,size(M_mat,1)

    JUMP_FLAG = 0

    ! Pick a target Sigma_bar
    IF (MPI_GLB_RANK.EQ.0) THEN

        EIG_BAR = M_eig(EIG_BAR_IND)
        IF (.NOT.((EIGRES(EIG_R_mat(:,EIG_BAR_IND),M_BAR)))) THEN
            JUMP_FLAG = 39
        ELSE
            open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
            WRITE(LOGID,*) 'TARGET EIG_BAR:', EIG_BAR, '(# ', EIG_BAR_IND, ')'
            close(LOGID)
        ENDIF
        
    ENDIF
    CALL MPI_BCAST(JUMP_FLAG,1,MPI_INT,0,MPI_COMM_WORLD,IERR)

    IF (JUMP_FLAG.EQ.39) THEN
        CYCLE loop_ind
    ENDIF
    ! ======================================================================

    IF (MPI_GLB_RANK.EQ.0) THEN
        ! Save rur_bar, rup_bar, uz_bar, ror_bar, rop_bar, oz_bar for DEGENERATE PERTURBATION THEORY
        CALL EIG2VELVOR(M_BAR,K_BAR,EIG_BAR_IND,EIG_R_mat,EIG_L_mat,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECR_BAR,VECL_BAR,comm_grp=newcomm)
    ENDIF

! ! {M_1, AK_1} & {M_2}
! !=======================================================================
DO M_1 = -M_BAR/2,(10-M_BAR); K_VAL_1_0 = -10.10; DO WHILE (K_VAL_1_0 .LT. 10.10)

    K_VAL_1_0 = NINT(10*(K_VAL_1_0 + 0.10))/10.D0
    ! IF ((M_BAR.EQ.M_1).AND.(K_BAR_SCAN.EQ.K_1_SCAN)) CYCLE
    IF ((M_1.EQ.0).AND.(K_VAL_1_0.EQ.0.D0)) CYCLE                   ! Not interested in m = 0, k = 0 mode
    IF ((M_1+M_BAR.EQ.0).AND.(K_VAL_1_0+K_BAR.EQ.0.D0)) CYCLE       ! Not interested in m = 0, k = 0 mode
    IF ((M_1.EQ.0).AND.(M_BAR.EQ.0.D0).AND.(K_BAR.LE.0.D0)) CYCLE   ! for m_b = m1 = m2 = 0, only need to check positive k_1

    ! ======================================================================
    CALL MPI_BCAST(EIG_BAR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)

    ! ! TUNE K
    ! !=======================================================================
    ! 0. Parameters:
    JUMP_FLAG = 0    
    tol = 1.0E-3  ! tolerance

    ! 1. Obtain initial eigvals and eigvecs:
    ! For M1,K1:
    IF (MPI_GLB_RANK.EQ.0) THEN
        open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
        WRITE(LOGID,*) 'M_BAR,KBAR: ',ITOA3(M_BAR),',',K_BAR,'; M_1,K_1: ', ITOA3(M_1),',',K_VAL_1_0
        close(LOGID)
    ENDIF
    ! CALL EIG_MATRIX(M_1,K_VAL_1_0,M_mat_1_0,M_eig_1_0,EIG_R_mat_1_0,EIG_L_mat_1_0,comm_grp=newcomm)
    CALL READ_MODE( M_1,K_VAL_1_0,M_mat_1_0,M_eig_1_0,EIG_R_mat_1_0,EIG_L_mat_1_0)
    ! For M2,K2:
    M_2 = M_1 + M_BAR
    K_VAL_2_0 = K_VAL_1_0 + K_BAR
    ! CALL EIG_MATRIX(M_2,K_VAL_2_0,M_mat_2_0,M_eig_2_0,EIG_R_mat_2_0,EIG_L_mat_2_0,comm_grp=newcomm)
    CALL READ_MODE( M_2,K_VAL_2_0,M_mat_2_0,M_eig_2_0,EIG_R_mat_2_0,EIG_L_mat_2_0)

    ! 2. Check if tolerence is reached:
    ! Match the IMAG part
    IF (MPI_GLB_RANK.EQ.0) THEN
        open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)

        ! Create resolving mask
        IF(.NOT.ALLOCATED(MASK_1)) THEN
            ALLOCATE(MASK_1(SIZE(EIG_R_mat_1_0,2)),MASK_2(SIZE(EIG_R_mat_2_0,2)))
        ENDIF
        MASK_1 = EIGRES(EIG_R_mat_1_0,M_1)
        MASK_2 = EIGRES(EIG_R_mat_2_0,M_2)

        ! Avoid same eigenmode for {mbar,kbar} and {m1,k1}
        IF ((M_BAR.EQ.M_1).AND.(ABS(K_BAR-K_VAL_1_0).LE.0.1)) THEN
            MASK_1 = MASK_1 .AND. (ABS(EIG_BAR-M_eig_1_0).GT.0.01)
        ENDIF
        
        ! Find the minimum distance
        tol_diff = MINVAL(ABS(AIMAG(MESHGRID(M_eig_1_0, M_eig_2_0) - EIG_BAR)),MESHGRID(MASK_1,MASK_2))/ABS(AIMAG(EIG_BAR))
        WRITE(LOGID,*) tol_diff
        IF (tol_diff.LT.tol) THEN

            tol_ind = MINLOC(ABS(AIMAG(MESHGRID(M_eig_1_0, M_eig_2_0) - EIG_BAR)),MESHGRID(MASK_1,MASK_2))
            JJ = tol_ind(1); KK = tol_ind(2)

            WRITE(LOGID,*) 'Desired tol reached: tol_diff = ', tol_diff
            WRITE(LOGID,*) 'DEGEN:'
            WRITE(LOGID,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
            WRITE(LOGID,201) M_1, K_VAL_1_0, JJ, M_eig_1_0(JJ)
            WRITE(LOGID,201) M_2, K_VAL_2_0, kk, M_eig_2_0(KK)
            201  FORMAT('Eig: M# = ',I3,'; AK = ',F20.13,'; Sig (#',I3') =',1P8E20.13)

            DEALLOCATE(MASK_1,MASK_2)
            tol = 999.0
        ENDIF

        close(LOGID)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    IF (tol.EQ.999.0) THEN
        goto 20
    ELSE
        goto 10
    ENDIF

! Saving and outputing the results in FFF Space
!=======================================================================   
    ! Find JJth Eigvector of {M1,K1} & KKth Eigvector of {M2,K2} 
20  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    IF (MPI_GLB_RANK.EQ.0) THEN
    open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)

    ! Calculate the FINAL RESULT:
    ! M_bar, K_bar: RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECL_BAR

    ! {M1,K1}
    CALL EIG2VELVOR(M_1,K_VAL_1_0,JJ,EIG_R_mat_1_0,EIG_L_mat_1_0,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_1,VECL_1,comm_grp=newcomm)

    ! {M2,K2}
    CALL EIG2VELVOR(M_2,K_VAL_2_0,KK,EIG_R_mat_2_0,EIG_L_mat_2_0,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_2,VECL_2,comm_grp=newcomm)
    CALL NONLIN_MK(M_2,K_VAL_2_0,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_B1,comm_grp=newcomm)
    ! CALL MCAT(UZ_2)

    ! {M1,K1}
    ALLOCATE(RUR_BARC(NDIMR),RUP_BARC(NDIMR),UZ_BARC(NDIMR),ROR_BARC(NDIMR),ROP_BARC(NDIMR),OZ_BARC(NDIMR))
    RUR_BARC = CONJG(RUR_BAR); RUP_BARC = CONJG(RUP_BAR); UZ_BARC = CONJG( UZ_BAR)
    ROR_BARC = CONJG(ROR_BAR); ROP_BARC = CONJG(ROP_BAR); OZ_BARC = CONJG( OZ_BAR)
    CALL NONLIN_MK(M_1,K_VAL_1_0,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_B2,comm_grp=newcomm)
    ! call mcat(VECR_B2)

    ! Sigma(1):
    EIG_PRIME_R = DOT_PRODUCT(VECL_1,VECR_B2)*DOT_PRODUCT(VECL_2,VECR_B1)
    WRITE(LOGID,*) ''
    WRITE(LOGID,*) " SIGMA(degen)'^2= ",EIG_PRIME_R
    WRITE(LOGID,*) " SIGMA(degen)'^2= ",EIG_PRIME_R/(DOT_PRODUCT(VECL_1,VECR_1)*DOT_PRODUCT(VECL_2,VECR_2))
    WRITE(LOGID,*) " SIGMA(degen)'  = ",EIG_PRIME_R**0.5
    ! Beta/Alpha: ALPHA * {M1,K1} + BETA * {M2,K2}
    ALPHA = 1.D0
    ! BETA = (DOT_PRODUCT(VECL_2,VECR_B1)/DOT_PRODUCT(VECL_1,VECR_B2))**0.5
    BETA = EIG_PRIME_R**0.5/DOT_PRODUCT(VECL_1, VECR_B2)
    ! EIG_R_mat_1_1(:,JJ) = EIG_R_mat_1_1(:,JJ)*ALPHA
    ! EIG_R_mat_2_1(:,KK) = EIG_R_mat_2_1(:,KK)*BETA
    WRITE(LOGID,*) " BETA(degen)    = ",BETA

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
        WRITE(FID,201) M_1, K_VAL_1_0, JJ, M_eig_1_0(JJ)
        WRITE(FID,201) M_2, K_VAL_2_0, kk, M_eig_2_0(KK)
        WRITE(FID,*) "SIGMA(degen)'= RE:",abs(real(EIG_PRIME_R**0.5)),',IM:',imag(EIG_PRIME_R**0.5)
        WRITE(FID,*) "BETA/ALPHA   = RE:",abs(real(BETA)),',IM:',imag(BETA)
        WRITE(FID,*) ''
        close(FID)

    close(LOGID)
    ENDIF

! DEALLOCATE 
!=======================================================================     
10  IF(ALLOCATED(M_mat_1_0)) DEALLOCATE(M_mat_1_0)
    IF(ALLOCATED(M_mat_2_0)) DEALLOCATE(M_mat_2_0)
    ! ======================================================================
    IF(ALLOCATED(RUR_BARC)) DEALLOCATE(RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC)
    IF(ALLOCATED(RUR_1)) DEALLOCATE(RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1)
    IF(ALLOCATED(RUR_2)) DEALLOCATE(RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2)
    IF(ALLOCATED(VECR_B1)) DEALLOCATE(VECR_B1,VECR_B2,VECL_1,VECR_1,VECL_2,VECR_2)
    IF(ALLOCATED(EIG_R_mat_1_0)) DEALLOCATE(EIG_R_mat_1_0)
    IF(ALLOCATED(EIG_L_mat_1_0)) DEALLOCATE(EIG_L_mat_1_0)    
    IF(ALLOCATED(EIG_R_mat_2_0)) DEALLOCATE(EIG_R_mat_2_0)
    IF(ALLOCATED(EIG_L_mat_2_0)) DEALLOCATE(EIG_L_mat_2_0)    
    ! ======================================================================
    IF(ALLOCATED(K_tune)) DEALLOCATE(K_tune)
    IF(ALLOCATED(K_tune_JJ)) DEALLOCATE(K_tune_JJ)
    IF(ALLOCATED(MASK_1)) DEALLOCATE(MASK_1,MASK_2)
    IF(ALLOCATED(M_eig_1_0)) DEALLOCATE(M_eig_1_0)
    IF(ALLOCATED(M_eig_2_0)) DEALLOCATE(M_eig_2_0)
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

END PROGRAM EVP_SCAN_0915
!=======================================================================
