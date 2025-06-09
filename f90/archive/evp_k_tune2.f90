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
    USE MOD_MISC
    USE MOD_BANDMAT
    USE MOD_EIG
    USE MOD_FD
    USE MOD_LIN_LEGENDRE
    USE MOD_SCALAR3
    USE MOD_LAYOUT
    USE MOD_LEGOPS
    USE MOD_MARCH
    USE MOD_DIAGNOSTICS
    USE MOD_INIT
    USE MOD_INIT_LOOP

    IMPLICIT NONE
    INTEGER     :: N_mat1,N_mat2,II,JJ,KK
    INTEGER     :: M_IND_bar,K_IND_bar,M_IND1,K_IND1,M_IND2,K_IND2,EIG_BAR_IND
    REAL(P8)    :: K_BAR,K_VAL_1_0,K_VAL_1_1,K_VAL_2_0,K_VAL_2_1
    REAL(P8)    :: K_shift,tol,tol_diff,K_delta
    REAL(P8)    :: EIG_LMR_1, ZLEN_BAR, EIG_PRIME
    COMPLEX(P8) :: EIG_BAR
    REAL(P8),DIMENSION(:),ALLOCATABLE:: K_tune, K_tune_JJ, tol_diff_JJ
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


    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC

    CALL execute_command_line('echo "%read.input" | ./bin/init_exec')
    CALL READIN(5)           ! @ MOD_INIT


! {M_bar, AK_bar}: BASE FLOW
!=======================================================================
    WRITE(*,*) ''        
    WRITE(*,*) 'STEP 1: SELECT {M_bar,AK_bar}'
    WRITE(*,*) ''
    WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M) INDEX; M_IND ='
    WRITE(*,102) MINC
    WRITE(*,101) NTCHOPDIM-1
    READ(5,'(I72)') M_IND_bar
    M_IND_bar = M_IND_bar + 1

    WRITE(*,*) ''
    WRITE(*,*) 'CHOOSE AN AXIAL WAVENUMBER (AK) INDEX; AK_IND ='
    WRITE(*,103) ZLEN
    WRITE(*,104) NXCHOP-1
    WRITE(*,105) NXCHOPDIM, NXCHOP-1
    WRITE(*,101) NXCHOPDIM-1
    READ(5,'(I72)') K_IND_bar
    WRITE(*,*) 'AK = '
    WRITE(*,'(F8.3)') K_IND_bar*2*PI/ZLEN
    K_IND_bar = K_IND_bar + 1
    ZLEN_BAR = ZLEN

! FORMAT LIST    
101  FORMAT(' *** POSSIBLE INDICES: 0 ~ ', I3)
102  FORMAT(' *** M  = ', I3, ' *M_IND')
103  FORMAT(' *** AK = 2*PI*K_IND/', F0.3)
104  FORMAT(' *** K_IND = INPUT          (INPUT <=', I3, ')')
105  FORMAT(' ***         -(', I3,' - INPUT) (INPUT  >', I3, ')')

    ! Pick a target Sigma_bar
    WRITE(*,*) ''
    WRITE(*,*) 'CALCULATING EIGENVALUES FOR {M_bar,AK_bar}. MIGHT TAKE MINUTES'
    CALL EIG_MATRIX(M_IND_bar,K_IND_bar,K_BAR,M_mat,M_eig,EIG_R_mat) ! calculate operator matrix and its corresponding eigenvalues
    WRITE(*,*) ''
    WRITE(*,*) 'EIGVALS FOR {M_bar,K_bar}:'
    CALL MCAT(M_eig) ! print all eigenvalues
    WRITE(*,*) 'CHOOSE THE TARGET EIGVAL #:'
    READ(5,'(I72)') EIG_BAR_IND
    EIG_BAR = M_eig(EIG_BAR_IND)
    WRITE(*,*) EIG_BAR
    DEALLOCATE(M_mat,M_eig)


! {M_1, AK_1} & {M_2}
!=======================================================================
    ! {M',AK'}: M' = M + M_bar
    WRITE(*,*) ''
    WRITE(*,*) 'STEP 2: SELECT M1'
    WRITE(*,*) ''
    WRITE(*,102) MINC
    WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M) INDEX; M_IND ='
    WRITE(*,101) NTCHOPDIM-1
    READ(5,'(I72)') M_IND1

    M_IND1 = M_IND1 + 1
    N_mat1 = NRCHOPS(M_IND1) ! determine the operator size for m1

    M_IND2 = M_IND1 + M_IND_bar - 1 ! index; MODE #(M_IND2 - 1)
    N_mat2 = NRCHOPS(M_IND2) ! determine the operator size for m2
        

! TUNE K
!=======================================================================
    CALL ROOT_RANGE(K_VAL_1_1,K_VAL_1_0,20.D0,2.D0,M_IND1,M_IND2,K_BAR,EIG_BAR)
    ! K_VAL_1_0 and K_VAL_1_1 are both for M1,AK1
    ! K_VAL_1_0 = 8.43 ! lower bound
    K_delta = 0.001
    tol = 1.0E-13  ! tolerance

    DO II = 1,100
    ! For M1,K1:
        K_VAL_1_1 = K_VAL_1_0 + K_delta

        ZLEN = 2*PI/K_VAL_1_0
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND1,2,K_VAL_1_0,M_mat_1_0,M_eig_1_0,EIG_R_mat_1_0,EIG_L_mat_1_0)

        ZLEN = 2*PI/K_VAL_1_1
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND1,2,K_VAL_1_1,M_mat_1_1)

        IF(.NOT.ALLOCATED(M_mat_diff_1)) ALLOCATE(M_mat_diff_1(2*N_mat1,2*N_mat1))
        M_mat_diff_1 = 0.D0
        M_mat_diff_1 = MATMUL((M_mat_1_1 - M_mat_1_0), EIG_R_mat_1_0)
    
    ! For M2,K2:
        K_VAL_2_0 = K_VAL_1_0 + K_BAR
        K_VAL_2_1 = K_VAL_2_0 + K_delta

        ZLEN = 2*PI/K_VAL_2_0
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND2,2,K_VAL_2_0,M_mat_2_0,M_eig_2_0,EIG_R_mat_2_0,EIG_L_mat_2_0)

        ZLEN = 2*PI/K_VAL_2_1
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND2,2,K_VAL_2_1,M_mat_2_1)

        IF(.NOT.ALLOCATED(M_mat_diff_2)) ALLOCATE(M_mat_diff_2(2*N_mat2,2*N_mat2))
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

    ! Calculate the actual 
        K_VAL_1_1 = K_VAL_1_0 + K_shift
        ZLEN = 2*PI/K_VAL_1_1
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND1,2,K_VAL_1_1,M_mat_1_1,M_eig_1_1)

        K_VAL_2_1 = K_VAL_1_1 + K_BAR
        ZLEN = 2*PI/K_VAL_2_1
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND2,2,K_VAL_2_1,M_mat_2_1,M_eig_2_1)

    ! Check if tolerance is reached
    ! Match the IMAG part
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
            GOTO 20
        ENDIF
        K_VAL_1_0 = K_VAL_1_1

    ENDDO

! Maximum number of steps allowed
!=======================================================================
    IF (II.EQ.100) THEN ! Exceeds maximum number of steps allowed
        WRITE(*,*) 'Seem trapped'
        WRITE(*,*) 'tol_diff:'
        WRITE(*,*) tol_diff
    ENDIF

! Saving and outputing the results in FFF Space
!=======================================================================   
    ! Find JJth Eigvector of {M1,K1} & KKth Eigvector of {M2,K2} 
20  DO JJ = 1,2*N_mat1
        tol_diff_JJ(JJ) = MINVAL(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)))
        K_tune(JJ) = MINLOC(ABS(AIMAG(M_eig_2_1 - M_eig_1_1(JJ) - EIG_BAR)),1)
    ENDDO
    
    JJ = MINLOC(ABS(tol_diff_JJ),1)
    KK = K_tune(JJ)
    !CALL SORT(0,M_eig_1_1)
    CALL MSAVE(M_eig_1_1,'M_eig_1_1.dat')
    !CALL SORT(0,M_eig_2_1)
    CALL MSAVE(M_eig_2_1,'M_eig_2_1.dat')
    WRITE(*,201) M_IND_bar-1, K_BAR, EIG_BAR_IND, EIG_BAR
    WRITE(*,201) M_IND1-1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
    WRITE(*,201) M_IND2-1, K_VAL_2_1, kk, M_eig_2_1(KK)
201  FORMAT('Eig: M# = ',I3,'; AK = ',F8.3,'; Sig (#',I3') =',1P8E13.6)
    
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

    DEALLOCATE(EIG_R_mat, Eig_01, Eig_02)
    DEALLOCATE(K_tune, K_tune_JJ, tol_diff_JJ)
    DEALLOCATE(M_eig_1_0, M_eig_1_1)
    DEALLOCATE(M_eig_2_0, M_eig_2_1)
    DEALLOCATE(M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0)
    DEALLOCATE(M_mat_2_0, EIG_R_mat_2_0, EIG_L_mat_2_0)
    DEALLOCATE(M_mat_1_1, M_mat_diff_1, EIG_R_mat_1_1, EIG_L_mat_1_1)     
    DEALLOCATE(M_mat_2_1, M_mat_diff_2, EIG_R_mat_2_1, EIG_L_mat_2_1)

10    WRITE(*,*) ''
    WRITE(*,*) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC

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

    ZLEN = 2*PI/KVAL ! K1
    CALL INIT_LOOP()
    CALL EIG_MATRIX(MIND1,2,KVAL,Mmat1,Meig1)
    WRITE(*,*) ''
    WRITE(*,*) KVAL

    ZLEN = 2*PI/(KVAL+KBAR) ! K2
    CALL INIT_LOOP()
    CALL EIG_MATRIX(MIND2,2,KVAL2,Mmat2,Meig2)

    ALLOCATE(EIG_DIST_LIST(SIZE(Meig1)))
    ALLOCATE(EIG_DIST_LIST_I(SIZE(Meig2)))

    DO III = 1,SIZE(Meig1)
        EIG_DIST_LIST_I = AIMAG(Meig1(III)+SIG-Meig2)
        EIG_DIST_LIST(III) = EIG_DIST_LIST_I(MINLOC(ABS(EIG_DIST_LIST_I),1)) ! CLOSEST DISTANCE for III
    ENDDO
    TARGET_DIST = EIG_DIST_LIST(MINLOC(ABS(EIG_DIST_LIST),1)) ! CLOSEST DISTANCE overall
    WRITE(*,*) TARGET_DIST

    DEALLOCATE(Meig1,Meig2,Mmat1,Mmat2,EIG_DIST_LIST)
END FUNCTION TARGET_DIST
!=======================================================================
SUBROUTINE ROOT_RANGE(LB,UB,UMAX,STEP,MIND1,MIND2,KBAR,SIG)
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
            WRITE(*,*) 'ROOT_RANGE completed'
            GOTO 211
        ELSEIF (UB.GE.UMAX) THEN
            WRITE(*,*) 'ROOT_RANGE: CANT FIND ZERO'
            RETURN
        ENDIF
    ENDDO

211 RETURN
END SUBROUTINE ROOT_RANGE   
!=======================================================================
SUBROUTINE EIG_MATRIX(M_INDEX,K_INDEX,AKREAD,H,EIG_VAL,EIG_VEC_R,EIG_VEC_L)
    IMPLICIT NONE
    INTEGER,INTENT(IN)    :: M_INDEX,K_INDEX
    REAL(P8),INTENT(INOUT)    :: AKREAD
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: H
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT),OPTIONAL:: EIG_VAL
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT),OPTIONAL:: EIG_VEC_L,EIG_VEC_R

    INTEGER     :: I,MREAD,CHOP_N,NN
    REAL(P8)    :: REI
    CHARACTER*72:: STR1, STR2, STR3, STR4    
    TYPE(SCALAR):: PSI0,CHI0,RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
    TYPE(SCALAR):: PSIU,CHIU,RURU,RUPU,UZU,RORU,ROPU,OZU
    TYPE(SCALAR):: PSI1,CHI1,PSI2,CHI2,PSI3,CHI3


    ! ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
    ! CALL READIN(5)           ! @ MOD_INIT
    ! CALL execute_command_line('echo "%read.input"')

    CALL ALLOCATE(PSI0)
    CALL ALLOCATE(CHI0)
    CALL ALLOCATE(RUR0)
    CALL ALLOCATE(RUP0)
    CALL ALLOCATE( UZ0)
    CALL ALLOCATE(ROR0)
    CALL ALLOCATE(ROP0)
    CALL ALLOCATE( OZ0)

    ! LOAD INITIAL PSI (PSI0) AND CHI (CHI0) DATA
    ! PERFORM 'INIT' AS A PREREQUISITE TO MAKREADE PSI0 AND CHI0
    CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI0,PSI0)
    CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI0,CHI0)

    CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)
    CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0)

    MREAD = M(M_INDEX)
    AKREAD = AK(M_INDEX,K_INDEX)
    CHOP_N = 3

    NN = NRCHOPS(M_INDEX)
    IF (.NOT. ALLOCATED(H))         ALLOCATE(H(2*NN,2*NN))
    H = 0.D0

    CALL ALLOCATE(PSIU)
    CALL ALLOCATE(CHIU)
    CALL ALLOCATE(RURU)
    CALL ALLOCATE(RUPU)
    CALL ALLOCATE( UZU)
    CALL ALLOCATE(RORU)
    CALL ALLOCATE(ROPU)
    CALL ALLOCATE( OZU)

    CALL ALLOCATE(PSI1)
    CALL ALLOCATE(CHI1)

    CALL ALLOCATE(PSI2)
    CALL ALLOCATE(CHI2)

    CALL ALLOCATE(PSI3)
    CALL ALLOCATE(CHI3)

    ! Use test vector to determine operator H
    DO I = 1,NN
        PSIU%LN = 0.D0 
        CHIU%LN = 0.D0 ! DEL2CHI
        PSIU%E = 0.D0
        CHIU%E = 0.D0 ! DEL2CHI
        PSIU%E(I,M_INDEX,K_INDEX) = 1.D0

        CALL CHOPSET(CHOP_N)

        ! Del2[Psi,Del2Chi]/Re term: PSI
        IF (VISC%SW.EQ.1) THEN
        IF(VISC%NU.EQ.0.0) THEN
            WRITE(*,*) 'VISC: VISC%NU CANNOT BE ZERO'
            STOP
        ENDIF
        REI = VISC%NU
        CALL DEL2(PSIU,PSI3)
        CALL DEL2(CHIU,CHI3) ! DEL2[DEL2CHI]
        PSI3%E = PSI3%E * VISC%NU ! NU = 1/RE
        CHI3%E = CHI3%E * VISC%NU ! DEL2[DEL2CHI]/RE
        ELSEIF (VISC%SW.EQ.2) THEN
        REI = VISC%NUP
            CALL HELMP(VISC%P,PSIU,PSI3,0.D0)
            CALL HELMP(VISC%P,CHIU,CHI3,0.D0)
            PSI3%E = PSI3%E * VISC%NUP ! NU = 1/RE
            CHI3%E = CHI3%E * VISC%NUP ! DEL2[DEL2CHI]/RE
        ELSE
        REI = 0.D0
        PSI3%E = 0.D0
        CHI3%E = 0.D0          
        PSI3%LN = 0.D0
        CHI3%LN = 0.D0
        ENDIF

        ! -W x U' + U X W' term: PSI
        ! NOTE: For this test vecotr, DEL2CHI = 0, so CHI = 0.
        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU) ! CHI (= 0 = DEL2CHI)
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU) ! CHI (= 0 = DEL2CHI)
        
        ! CALL VEL2VOR(RURU,RUPU,UZU,RORU,ROPU,OZU)
        
        CALL VPROD_PFF(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        ! CALL TOFP(RURU)
        ! CALL TOFP(RUPU)
        ! CALL TOFP(UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)

        CALL VPROD_PFF(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        ! CALL TOFP(RORU)
        ! CALL TOFP(ROPU)
        ! CALL TOFP(OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-CHOP_N)

        PSI1%E = -PSI1%E + PSI2%E + PSI3%E
        CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NN,  I) = PSI1%E(:NN,M_INDEX,K_INDEX)
        H(NN+1:,I) = CHI1%E(:NN,M_INDEX,K_INDEX) ! DEL2CHI
        
        PSIU%LN = 0.D0
        CHIU%LN = 0.D0
        PSIU%E = 0.D0
        CHIU%E = 0.D0
        CHIU%E(I,M_INDEX,K_INDEX) = 1.D0 ! DEL2CHI

        ! Del2[Psi,Del2Chi]/Re term: CHI
        CALL CHOPSET(CHOP_N)
        IF (VISC%SW.EQ.1) THEN
        CALL DEL2(PSIU,PSI3)
        CALL DEL2(CHIU,CHI3) ! DEL2[DEL2CHI]
        PSI3%E = PSI3%E * VISC%NU ! NU = 1/RE
        CHI3%E = CHI3%E * VISC%NU ! DEL2[DEL2CHI]/RE
        ELSEIF (VISC%SW.EQ.2) THEN
            CALL HELMP(VISC%P,PSIU,PSI3,0.D0)
            CALL HELMP(VISC%P,CHIU,CHI3,0.D0)
            PSI3%E = PSI3%E * VISC%NUP ! NU = 1/RE
            CHI3%E = CHI3%E * VISC%NUP ! DEL2[DEL2CHI]/RE
        ELSE
        PSI3%E = 0.D0
        CHI3%E = 0.D0          
        PSI3%LN = 0.D0
        CHI3%LN = 0.D0
        ENDIF
        
        CALL CHOPSET(-CHOP_N)
        CALL IDEL2(CHIU,CHIU) ! DEL2CHI -> CHI so that PC2VOR and PC2VEL calcualte the correct quantities  
        CALL CHOPSET(CHOP_N)

        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU) ! CHI
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU) ! CHI

        ! CALL VEL2VOR(RURU,RUPU,UZU,RORU,ROPU,OZU)

        CALL VPROD_PFF(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        ! CALL TOFP(RURU)
        ! CALL TOFP(RUPU)
        ! CALL TOFP(UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)

        CALL VPROD_PFF(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        ! CALL TOFP(RORU)
        ! CALL TOFP(ROPU)
        ! CALL TOFP(OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-CHOP_N)

        PSI1%E = -PSI1%E + PSI2%E + PSI3%E
        CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NN,  I+NN) &
                            = PSI1%E(:NN,M_INDEX,K_INDEX)
        H(NN+1:,I+NN) &
                            = CHI1%E(:NN,M_INDEX,K_INDEX) ! DEL2CHI

    ENDDO

! OUTPUT:
WRITE( STR1, '(f10.2)' )  AKREAD  
WRITE( STR2, '(f10.2)' )  ELL
WRITE( STR3, '(f10.6)' )  REI
WRITE( STR4, '(f10.3)' )  QPAIR%H(1)

! CALL MSAVE(H, TRIM(ADJUSTL(FILES%SAVEDIR))// 'mtrx_r_m_' // ITOA4(MREAD) //&
!         '_k_'       // TRIM(ADJUSTL(STR1))  //&
!         '_L_'       // TRIM(ADJUSTL(STR2))  //&
!         '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
!         '_QI_'      // TRIM(ADJUSTL(STR4))  //&
!         '_N_'       // ITOA4(NRCHOPS(M_INDEX)) // '.dat')

IF (PRESENT(EIG_VAL)) THEN    
    IF (.NOT. ALLOCATED(EIG_VAL))   ALLOCATE(EIG_VAL(2*NN))
    EIG_VAL = 0.D0
    EIG_VAL = GENEIG(H)
    ! CALL MSAVE(EIG_VAL, TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigval_m_'// ITOA4(MREAD) //&
    !         '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !         '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !         '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !         '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !         '_N_'       // ITOA4(NRCHOPS(M_INDEX)) // '.dat')
ENDIF

IF (PRESENT(EIG_VEC_R)) THEN
    IF (.NOT. ALLOCATED(EIG_VEC_R)) ALLOCATE(EIG_VEC_R(2*NN,2*NN))
    EIG_VEC_R = 0.D0
    EIG_VEC_R = EIGVEC('R',H)
    ! CALL MSAVE(EIG_VEC_R, TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_FFF_r_m_'// ITOA4(MREAD) //&
    !     '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !     '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !     '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !     '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !     '_N_'       // ITOA4(NRCHOPS(M_INDEX)) // '.dat')
ENDIF

IF (PRESENT(EIG_VEC_L)) THEN
    IF (.NOT. ALLOCATED(EIG_VEC_L)) ALLOCATE(EIG_VEC_L(2*NN,2*NN))
    EIG_VEC_L = 0.D0
    EIG_VEC_L = EIGVEC('L',H)
    ! CALL MSAVE(EIG_VEC_L, TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_FFF_l_m_'// ITOA4(MREAD) //&
    !     '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !     '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !     '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !     '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !     '_N_'       // ITOA4(NRCHOPS(M_INDEX)) // '.dat')
ENDIF

! DEALLOCATE:
    CALL DEALLOCATE(PSI1)
    CALL DEALLOCATE(CHI1)
    CALL DEALLOCATE(PSI2)
    CALL DEALLOCATE(CHI2)
    CALL DEALLOCATE(PSI3)
    CALL DEALLOCATE(CHI3)

    CALL DEALLOCATE(PSIU)
    CALL DEALLOCATE(CHIU)
    CALL DEALLOCATE(RURU)
    CALL DEALLOCATE(RUPU)
    CALL DEALLOCATE( UZU)
    CALL DEALLOCATE(RORU)
    CALL DEALLOCATE(ROPU)
    CALL DEALLOCATE( OZU)

    CALL DEALLOCATE(PSI0)
    CALL DEALLOCATE(CHI0)
    CALL DEALLOCATE(RUR0)
    CALL DEALLOCATE(RUP0)
    CALL DEALLOCATE( UZ0)
    CALL DEALLOCATE(ROR0)
    CALL DEALLOCATE(ROP0)
    CALL DEALLOCATE( OZ0)

    RETURN

END SUBROUTINE EIG_MATRIX
!=======================================================================
END PROGRAM EVP_K_TUNE
!=======================================================================
