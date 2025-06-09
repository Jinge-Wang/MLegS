PROGRAM EVP_K_TUNE
!=======================================================================
! [USAGE]: 
! USE VPROD_PFF RATHER THAN VPROD
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
    REAL(P8)    :: EIG_LMR_1, ZLEN_BAR, EIG_PRIME
    COMPLEX(P8) :: EIG_BAR
    REAL(P8),DIMENSION(:),ALLOCATABLE:: K_tune, K_tune_JJ, tol_diff_JJ
    REAL(P8),DIMENSION(:),ALLOCATABLE:: RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_1_0, M_eig_1_1, M_eig_1_1_ND
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_2_0, M_eig_2_1, M_eig_2_1_ND
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: Eig_01, Eig_02
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat, EIG_R_mat
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_diff_1,M_mat_diff_2
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_2_0, EIG_R_mat_2_0, EIG_L_mat_2_0
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1
    TYPE(SCALAR) :: PSI_BAR,CHI_BAR
    TYPE(SCALAR) :: PSI_01, PSI_11, PSI_02, PSI_21
    TYPE(SCALAR) :: CHI_01, CHI_11, CHI_02, CHI_21
    
    INTEGER:: MPI_GLB_PROCS, MPI_GLB_RANK, newcomm


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


! {M_bar, AK_bar}: BASE FLOW
!=======================================================================
    IF (MPI_GLB_RANK.EQ.0) THEN
        WRITE(*,*) ''        
        WRITE(*,*) 'STEP 1: SELECT {M_bar,AK_bar}'
        WRITE(*,*) ''
        WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M >= 0) ; M ='
        READ(5,'(I72)') M_BAR

        WRITE(*,*) ''
        WRITE(*,103) ZLEN
    103  FORMAT(' AK = 2*PI/', F0.3)    
        WRITE(*,*) 'CHOOSE AN AXIAL WAVENUMBER (AK) INDEX; AK ='
        READ(5,'(F5.3)') K_BAR

        WRITE(*,*) ''
        WRITE(*,*) 'CALCULATING EIGENVALUES FOR {M_bar,AK_bar}. MIGHT TAKE MINUTES'    
    ENDIF

    CALL EIG_MATRIX(M_BAR,K_BAR,M_mat,M_eig,EIG_R_mat,comm_grp=newcomm) ! calculate operator matrix and its corresponding eigenvalues

    ! Pick a target Sigma_bar
    IF (MPI_GLB_RANK.EQ.0) THEN
        WRITE(*,*) ''
        WRITE(*,*) 'EIGVALS FOR {M_bar,K_bar}:'
        CALL MCAT(M_eig) ! print all eigenvalues
        WRITE(*,*) 'CHOOSE THE TARGET EIGVAL #:'
        READ(5,'(I72)') EIG_BAR_IND
        EIG_BAR = M_eig(EIG_BAR_IND)
        WRITE(*,*) EIG_BAR
        DEALLOCATE(M_eig)


! {M_1, AK_1} & {M_2}
!=======================================================================
        ! {M',AK'}: M' = M + M_bar
        WRITE(*,*) ''
        WRITE(*,*) 'STEP 2: SELECT M1'
        WRITE(*,*) ''
        WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M >= 0) ; M ='
        READ(5,'(I72)') M_1
        
    ENDIF
    DEALLOCATE(M_mat)
    CALL MPI_BCAST(EIG_BAR,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(M_1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    M_2 = M_1 + M_BAR

! TUNE K
!=======================================================================
    CALL ROOT_RANGE(K_VAL_1_1,K_VAL_1_0,20.D0,2.D0,M_1,M_2,K_BAR,EIG_BAR)
    ! K_VAL_1_0 and K_VAL_1_1 are both for M1,AK1
    ! K_VAL_1_0 = 8.43 ! lower bound
    K_delta = 0.001
    tol = 1.0E-13  ! tolerance

    DO II = 1,100
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
            N_mat1 = SIZE(M_eig_1_0)
            ALLOCATE(M_mat_diff_1(2*N_mat1,2*N_mat1))
        ENDIF
        M_mat_diff_1 = 0.D0
        M_mat_diff_1 = MATMUL((M_mat_1_1 - M_mat_1_0), EIG_R_mat_1_0)

        IF(.NOT.ALLOCATED(M_mat_diff_2)) THEN
            N_mat2 = SIZE(M_eig_2_0)
            ALLOCATE(M_mat_diff_2(2*N_mat2,2*N_mat2))
        ENDIF
        M_mat_diff_2 = 0.D0
        M_mat_diff_2 = MATMUL((M_mat_2_1 - M_mat_2_0), EIG_R_mat_2_0)  
        

    ! Find possible K_shift
        IF(.NOT.ALLOCATED(K_tune)) ALLOCATE(K_tune(N_mat1*2))
        IF(.NOT.ALLOCATED(K_tune_JJ)) ALLOCATE(K_tune_JJ(N_mat2*2))

        !$OMP parallel do private(EIG_LMR_1,KK,K_tune_JJ)
        DO JJ = 1,2*N_mat1
            
            EIG_LMR_1 = 0.D0
            EIG_LMR_1 = &
            & AIMAG(DOT_PRODUCT(EIG_L_mat_1_0(:,JJ),M_mat_diff_1(:,JJ))/DOT_PRODUCT(EIG_L_mat_1_0(:,JJ),EIG_R_mat_1_0(:,JJ)))

            K_tune_JJ = 0.D0
            DO KK = 1,2*N_mat2
                
                IF (ABS(REAL(M_eig_2_0(KK) - M_eig_1_0(JJ) - EIG_BAR)) .GT. 1.0E-8) THEN ! Arbitrary threshold
                    K_tune_JJ(KK) = 100.D0
                ELSE
                    K_tune_JJ(KK) = K_delta*AIMAG(M_eig_2_0(KK) - M_eig_1_0(JJ) - EIG_BAR) / ( &
                                & EIG_LMR_1 - &
                                & AIMAG(DOT_PRODUCT(EIG_L_mat_2_0(:,KK),M_mat_diff_2(:,KK))/DOT_PRODUCT(EIG_L_mat_2_0(:,KK),EIG_R_mat_2_0(:,KK))))
                ENDIF

            ENDDO

            K_tune(JJ) = K_tune_JJ(MINLOC(ABS(K_tune_JJ),1))

        ENDDO
        !$OMP end parallel do

        K_shift = K_tune(MINLOC(ABS(K_tune),1))

        WRITE(*,*) 'K_shift:'
        WRITE(*,*) K_shift
        ! WRITE(*,*) 'Original diff:'
        ! WRITE(*,*) MINVAL(ABS(AIMAG(M_eig_1_0-EIG_BAR))/ABS(AIMAG(EIG_BAR)))
    ENDIF

    ! Calculate the actual 
        K_VAL_1_1 = K_VAL_1_0 + K_shift
        CALL EIG_MATRIX(M_1,K_VAL_1_1,M_mat_1_1,M_eig_1_1,comm_grp=newcomm)

        K_VAL_2_1 = K_VAL_1_1 + K_BAR
        CALL EIG_MATRIX(M_2,K_VAL_2_1,M_mat_2_1,M_eig_2_1,comm_grp=newcomm)

    ! Check if tolerance is reached
    ! Match the IMAG part
    IF (MPI_GLB_RANK.EQ.0) THEN
        IF(.NOT.ALLOCATED(tol_diff_JJ)) ALLOCATE(tol_diff_JJ(N_mat1*2))

        DO JJ = 1,2*N_mat1
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
    IF (II.EQ.100) THEN ! Exceeds maximum number of steps allowed
        IF (MPI_GLB_RANK.EQ.0) THEN
            WRITE(*,*) 'Seem trapped'
            WRITE(*,*) 'tol_diff:'
            WRITE(*,*) tol_diff
        ENDIF
    ENDIF

! Saving and outputing the results in FFF Space
!=======================================================================   
    ! Find JJth Eigvector of {M1,K1} & KKth Eigvector of {M2,K2} 
20  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    IF (MPI_GLB_RANK.EQ.0) THEN
        DO JJ = 1,2*N_mat1
            tol_diff_JJ(JJ) = MINVAL(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)))
            K_tune(JJ) = MINLOC(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)),1)
        ENDDO
        
        JJ = MINLOC(ABS(tol_diff_JJ),1)
        KK = K_tune(JJ)
        !CALL SORT(0,M_eig_1_1)
        CALL MSAVE(M_eig_1_1,'M_eig_1_1.dat')
        !CALL SORT(0,M_eig_2_1)
        CALL MSAVE(M_eig_2_1,'M_eig_2_1.dat')
        WRITE(*,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
        WRITE(*,201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
        WRITE(*,201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)
201  FORMAT('Eig: M# = ',I3,'; AK = ',F8.3,'; Sig (#',I3') =',1P8E13.6)
    ENDIF
    
! ! Obtain Sigma(1) in PFF space
! !=======================================================================    
! ! Calculate the FINAL RESULT:

!     ! M_bar, K_bar -> Obtain PSI_BAR, CHI_BAR
!     ZLEN = 2*PI/K_BAR
!     CALL INIT_LOOP() ! Reset Psi0 and Chi0 to read.input; Adjust operators 
!     CALL ALLOCATE(PSI_BAR)
!     CALL ALLOCATE(CHI_BAR)  
!     PSI_BAR%LN = 0.D0
!     CHI_BAR%LN = 0.D0
!     PSI_BAR%E = 0.D0
!     CHI_BAR%E = 0.D0
!     PSI_BAR%E(:NRCHOPS(M_IND_bar),M_IND_bar,2) = EIG_R_mat(:NRCHOPS(M_IND_bar),EIG_BAR_IND)
!     CHI_BAR%E(:NRCHOPS(M_IND_bar),M_IND_bar,2) = EIG_R_mat(NRCHOPS(M_IND_bar)+1:,EIG_BAR_IND) !DEL2(CHI)
!     CALL IDEL2(CHI_BAR,CHI_BAR)
!     CALL RTRAN(PSI_BAR,1)
!     CALL RTRAN(CHI_BAR,1) ! PFF
!     PSI_BAR%E(:,1,1) = PSI_BAR%E(:,M_IND_bar,2)
!     CHI_BAR%E(:,1,1) = CHI_BAR%E(:,M_IND_bar,2)
!     PSI_BAR%E(:,M_IND_bar,2) = 0.D0
!     CHI_BAR%E(:,M_IND_bar,2) = 0.D0
!     CALL TOFF(PSI_BAR)
!     CALL TOFF(CHI_BAR) ! FFF; New base flow (NEED DOUBLE CHECK)

    
!     ! M1,K1
!     ZLEN = 2*PI/K_VAL_1_1
!     CALL INIT_LOOP()
!     CALL ALLOCATE(PSI_11)
!     CALL ALLOCATE(CHI_11) 
!     CALL EIG_MATRIX(M_IND1,2,K_VAL_1_1,M_mat_1_1,M_eig_1_1,EIG_R_mat_1_1,EIG_L_mat_1_1)
!     PSI_11%E(:NRCHOPS(M_IND1),M_IND1,2) = EIG_L_mat_1_1(:NRCHOPS(M_IND1),JJ)
!     CHI_11%E(:NRCHOPS(M_IND1),M_IND1,2) = EIG_L_mat_1_1(NRCHOPS(M_IND1)+1:,JJ)
!     CALL RTRAN(PSI_11,1)
!     CALL RTRAN(CHI_11,1)        

!     CALL MSAVE(PSI_BAR, TRIM(ADJUSTL(FILES%SAVEDIR))// FILES%PSI0)
!     CALL MSAVE(CHI_BAR, TRIM(ADJUSTL(FILES%SAVEDIR))// FILES%CHI0)
!     CALL EIG_MATRIX(M_IND1,2,K_VAL_1_1,M_mat_1_1)

!     CALL ALLOCATE(PSI_01)
!     CALL ALLOCATE(CHI_01) 
!     ALLOCATE(EIG_01(2*N_mat1))
!     EIG_01 = MATMUL(M_mat_1_1(:,:),EIG_R_mat_1_1(:,JJ))
!     PSI_01%E(:NRCHOPS(M_IND1),M_IND1,2) = EIG_01(:NRCHOPS(M_IND1))
!     CHI_01%E(:NRCHOPS(M_IND1),M_IND1,2) = EIG_01(NRCHOPS(M_IND1)+1:)
!     CALL RTRAN(PSI_01,1)
!     CALL RTRAN(CHI_01,1)


!     ! M2,K2
!     ZLEN = 2*PI/K_VAL_2_1
!     CALL INIT_LOOP()
!     CALL ALLOCATE(PSI_21)
!     CALL ALLOCATE(CHI_21) 
!     CALL EIG_MATRIX(M_IND2,2,K_VAL_2_1,M_mat_2_1,M_eig_2_1,EIG_R_mat_2_1,EIG_L_mat_2_1)
!     PSI_21%E(:NRCHOPS(M_IND2),M_IND2,2) = EIG_L_mat_2_1(:NRCHOPS(M_IND2),KK)
!     CHI_21%E(:NRCHOPS(M_IND2),M_IND2,2) = EIG_L_mat_2_1(NRCHOPS(M_IND2)+1:,KK)
!     CALL RTRAN(PSI_21,1)
!     CALL RTRAN(CHI_21,1)

!     PSI_BAR%E(:,1,1) = CONJG(PSI_BAR%E(:,1,1))
!     CHI_BAR%E(:,1,1) = CONJG(CHI_BAR%E(:,1,1))
!     CALL MSAVE(PSI_BAR, TRIM(ADJUSTL(FILES%SAVEDIR))// FILES%PSI0)
!     CALL MSAVE(CHI_BAR, TRIM(ADJUSTL(FILES%SAVEDIR))// FILES%CHI0)
!     CALL EIG_MATRIX(M_IND2,2,K_VAL_2_1,M_mat_2_1)

!     CALL ALLOCATE(PSI_02)
!     CALL ALLOCATE(CHI_02) 
!     ALLOCATE(EIG_02(2*N_mat2))
!     EIG_02 = MATMUL(M_mat_2_1(:,:),EIG_R_mat_2_1(:,KK))
!     PSI_02%E(:NRCHOPS(M_IND2),M_IND2,2) = EIG_02(:NRCHOPS(M_IND2))
!     CHI_02%E(:NRCHOPS(M_IND2),M_IND2,2) = EIG_02(NRCHOPS(M_IND2)+1:)
!     CALL RTRAN(PSI_02,1)
!     CALL RTRAN(CHI_02,1)

!     ! ! Find JJth Eigvector of {M1,K1} & KKth Eigvector of {M2,K2}
!     ! DO JJ = 1,2*N_mat1
!     !     tol_diff_JJ(JJ) = MINVAL(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)))
!     !     K_tune(JJ) = MINLOC(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)),1)
!     ! ENDDO
!     ! JJ = MINLOC(ABS(tol_diff_JJ),1)
!     ! KK = K_tune(JJ)

!     EIG_PRIME = &
!     ( ( DOT_PRODUCT(PSI_01%E(:,M_IND1,2),PSI_21%E(:,M_IND2,2)) + DOT_PRODUCT(CHI_01%E(:,M_IND1,2),CHI_21%E(:,M_IND2,2)) ) * &
!     & ( DOT_PRODUCT(PSI_02%E(:,M_IND2,2),PSI_11%E(:,M_IND1,2)) + DOT_PRODUCT(CHI_02%E(:,M_IND2,2),CHI_11%E(:,M_IND1,2)) ) / & 
!     & ( DOT_PRODUCT(EIG_L_mat_1_1(:,JJ),EIG_R_mat_1_1(:,JJ))   * DOT_PRODUCT(EIG_L_mat_2_1(:,KK),EIG_R_mat_2_1(:,KK))   ) )**0.5

!     WRITE(*,*) 'EIG_PRIME = '
!     WRITE(*,*) EIG_PRIME
!     WRITE(*,*) 'M1, K1: '
!     WRITE(*,*) M_eig_1_1(JJ)
!     WRITE(*,*) 'M2, K2: '
!     WRITE(*,*) M_eig_2_1(KK)
!     WRITE(*,*) 'M_bar, K_bar: '
!     WRITE(*,*) EIG_BAR


!     ! (  DOT_PRODUCT(EIG_L_mat_1_1(:,JJ),MATMUL(M_mat_2_1(:,:),EIG_R_mat_2_1(:,KK))) * &
!     !  & DOT_PRODUCT(EIG_L_mat_2_1(:,KK),MATMUL(M_mat_1_1(:,:),EIG_R_mat_1_1(:,JJ))) / &
!     !  & DOT_PRODUCT(EIG_L_mat_1_1(:,JJ),EIG_R_mat_1_1(:,JJ)) / &
!     !  & DOT_PRODUCT(EIG_L_mat_2_1(:,KK),EIG_R_mat_2_1(:,KK)) )**0.5 ! PFF





! DEALLOCATE 
!======================================================================= 
    ! CALL DEALLOCATE(PSI_BAR)
    ! CALL DEALLOCATE(CHI_BAR) 
    ! CALL DEALLOCATE(PSI_01)
    ! CALL DEALLOCATE(CHI_01)
    ! CALL DEALLOCATE(PSI_02)
    ! CALL DEALLOCATE(CHI_02)
    ! CALL DEALLOCATE(PSI_11)
    ! CALL DEALLOCATE(CHI_11)
    ! CALL DEALLOCATE(PSI_21)
    ! CALL DEALLOCATE(CHI_21)

    ! DEALLOCATE(Eig_01, Eig_02)
    
10  DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
    DEALLOCATE(M_mat_1_0, M_mat_1_1, M_mat_2_0, M_mat_2_1)

   IF (MPI_GLB_RANK.EQ.0) THEN
        DEALLOCATE(EIG_R_mat)
        DEALLOCATE(K_tune, K_tune_JJ, tol_diff_JJ)
        DEALLOCATE(M_eig_1_0, M_eig_1_1)
        DEALLOCATE(M_eig_2_0, M_eig_2_1)
        DEALLOCATE(EIG_R_mat_1_0, EIG_L_mat_1_0)
        DEALLOCATE(EIG_R_mat_2_0, EIG_L_mat_2_0)
        DEALLOCATE(M_mat_diff_1) !, EIG_R_mat_1_1, EIG_L_mat_1_1)     
        DEALLOCATE(M_mat_diff_2) !, EIG_R_mat_2_1, EIG_L_mat_2_1)
        
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
SUBROUTINE ROOT_RANGE(LB,UB,UMAX,STEP,MIND1,MIND2,KBAR,SIG)
!CALL ROOT_RANGE(K_VAL_1_1,K_VAL_1_0,20.D0,2.D0,M_IND1,M_IND2,K_BAR,EIG_BAR)
    IMPLICIT NONE
    REAL(P8),INTENT(INOUT):: LB,UB
    
    INTEGER:: III,MIND1,MIND2
    COMPLEX(P8):: SIG
    REAL(P8):: STEP,LMIN,UMAX,LB_DIST,UB_DIST,KBAR

    LMIN = 4.0 ! Starting lower-bound for K
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
                CALL EIGENDECOMPOSE(H,EIG_VAL,EIG_VEC_R,EIG_VEC_L)
            ELSEIF (PRESENT(EIG_VEC_R)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE(EIG_VEC_R) 
                ALLOCATE( EIG_VEC_R(2*NR_MK,2*NR_MK) )
                CALL EIGENDECOMPOSE(H,EIG_VAL,EIG_VEC_R)
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
END PROGRAM EVP_K_TUNE
!=======================================================================
