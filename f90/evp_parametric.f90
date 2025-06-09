PROGRAM EVP_PARAMETRIC
!=======================================================================
! [USAGE]:
! Calculate interaction coefficients of triadic resonance between
! three prescribed modes.
! 
! Treat 0-mode as constant and form parametric instability between the 
! other two. Calculate the coupling factors between the two modes and
! create perturbation files for the triad and their correction terms
!
! [UPDATES]:
! LAST UPDATE ON FEB 28, 2023
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
INTEGER     :: II, JJ, KK
INTEGER     :: M_0, M_1, M_2
REAL(P8)    :: K_0, K_1, K_2, SIG_0, SIG_1, SIG_2, R_CTL, OMEGA_CTL, VZ_CTL
COMPLEX(P8) :: EIG_BAR, J1, J2, J0
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

! PARAMETRIC INSTABILITY
COMPLEX(P8) :: EIG_PRIME_R, ALPHA, BETA, ALPHA_L, BETA_L
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_1_CORRECT, VECR_2_CORRECT

! I/O:
INTEGER:: IS, FID
CHARACTER(LEN=72):: FILENAME
CHARACTER(len=6) :: K_VAL
LOGICAL:: DIAGNOST_SWITCH = .FALSE.

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

! SPATIAL SYMMETRY
QPAIR%H(1) = 0.40
M_0 = 1; M_1 = 2; M_2 = M_0 + M_1
K_0 = -6.0; K_1 = -3.0; K_2 = K_0 + K_1 ! ALL K = 0

! ! CRITICAL LAYER
! R_CTL = 0.5
! OMEGA_CTL = (1-EXP(-R_CTL**2))/(R_CTL**2)
! VZ_CTL = QPAIR%H(1) * EXP(-R_CTL**2)
! SIG_0 = -(M_0*OMEGA_CTL+K_0*VZ_CTL)
! SIG_1 = -(M_1*OMEGA_CTL+K_1*VZ_CTL)
! SIG_2 = -(M_2*OMEGA_CTL+K_2*VZ_CTL)
SIG_0 = 2.834
SIG_1 = -1.246
SIG_2 = 1.589

! CALCULATE EVP FOR EACH TRIAD
CALL EIG_MATRIX(M_0, K_0, M_mat_0, M_eig_0, EIG_R_mat_0, EIG_L_mat_0, comm_grp=newcomm) 
CALL EIG_MATRIX(M_1, K_1, M_mat_1, M_eig_1, EIG_R_mat_1, EIG_L_mat_1, comm_grp=newcomm)
CALL EIG_MATRIX(M_2, K_2, M_mat_2, M_eig_2, EIG_R_mat_2, EIG_L_mat_2, comm_grp=newcomm)

IF (MPI_GLB_RANK.EQ.0) THEN

    II = MINLOC(ABS(AIMAG(M_eig_0)-SIG_0),DIM=1)
    JJ = MINLOC(ABS(AIMAG(M_eig_1)-SIG_1),DIM=1)
    KK = MINLOC(ABS(AIMAG(M_eig_2)-SIG_2),DIM=1)
    ! II = MAXLOC(REAL(M_eig_0),DIM=1)
    ! JJ = MAXLOC(REAL(M_eig_1),DIM=1)
    ! KK = MAXLOC(REAL(M_eig_2),DIM=1)

    ! PRINT TRIAD INFO
    WRITE (*, *) 'RESONANT TRIAD:'
    WRITE (*, 201) M_0, K_0, II, M_eig_0(II)
    WRITE (*, 201) M_1, K_1, JJ, M_eig_1(JJ)
    WRITE (*, 201) M_2, K_2, kk, M_eig_2(KK)
201  FORMAT('Eig: M# = ',I3,'; AK = ',F20.13,'; Sig (#',I3') =',1P8E20.13)

    ! FINAL RESULT: 
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

    ! ! SAVE PERTURBATION
    ! ! AFTER EIG2VELVOR TO NORMALIZE THE PERTURBATION FIELDS
    ! CALL SAVE_PERTURB('./data/new_perturb.input', &
    !                     M_0, K_0, EIG_R_mat_0(:, II), M_eig_0(II), &
    !                     M_1, K_1, EIG_R_mat_1(:, JJ), M_eig_1(JJ), &
    !                     M_2, K_2, EIG_R_mat_2(:, KK), M_eig_2(kk))

    ! SAVE VELOCITTY
    WRITE(K_VAL,'(F06.2)') K_0; 
    CALL SAVE_VEL(RUR_0, RUP_0, UZ_0, './data/vel_MK_'//ITOA3(M_0)//'_'//K_VAL//'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    WRITE(K_VAL,'(F06.2)') K_1; 
    CALL SAVE_VEL(RUR_1, RUP_1, UZ_1, './data/vel_MK_'//ITOA3(M_1)//'_'//K_VAL//'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    WRITE(K_VAL,'(F06.2)') K_2; 
    CALL SAVE_VEL(RUR_2, RUP_2, UZ_2, './data/vel_MK_'//ITOA3(M_2)//'_'//K_VAL//'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! NONLINEAR TERMS:
    ! {M2,K2}: 0  X 1 -> 2
    CALL NONLIN_MK(M_2,K_2,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,RUR_1,RUP_1,UZ_1,ROR_1_COPY,ROP_1_COPY,OZ_1_COPY,VECR_01,comm_grp=newcomm)

    ! {M1,K1}: 0* X 2 -> 1
    CALL NONLIN_MK(M_1,K_1,CONJG(RUR_0),CONJG(RUP_0),CONJG(UZ_0),CONJG(ROR_0),CONJG(ROP_0),CONJG(OZ_0),RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_02,comm_grp=newcomm)

    ! {M0,K0}: 1* X 2 -> 0
    CALL NONLIN_MK(M_0,K_0,CONJG(RUR_1),CONJG(RUP_1),CONJG(UZ_1),CONJG(ROR_1),CONJG(ROP_1),CONJG(OZ_1),RUR_2,RUP_2,UZ_2,ROR_2_COPY,ROP_2_COPY,OZ_2_COPY,VECR_12,comm_grp=newcomm)
    
ENDIF ! MPI_GLB_RANK.EQ.0

RANK0: IF (MPI_GLB_RANK .EQ. 0) THEN

! ============================== RESULT ================================
WRITE (*, *) "NEARLY DEGEN:"
J0 = DOT_PRODUCT(VECL_0, VECR_12); WRITE (*,*) "J0: ", J0
J1 = DOT_PRODUCT(VECL_1, VECR_02); WRITE (*,*) "J1: ", J1
J2 = DOT_PRODUCT(VECL_2, VECR_01); WRITE (*,*) "J2: ", J2
! WRITE(*,*) "J1 + JB - J2: ", J1 + JB - J2

EIG_PRIME_R = (J1*J2)**0.5
ALPHA = 1.D0; BETA = 1.D0; !EIG_PRIME_R/J1
WRITE (*, *) " BETA(degen) = ", BETA

! ========================= CORRECTION TERMS ===========================
VECR_1_CORRECT = &
EIG_CORRECT_3(JJ, M_eig_1, EIG_R_mat_1, EIG_L_mat_1, VECR_02, M_1, K_1, ALPHA, BETA, newcomm, .FALSE.)
VECR_2_CORRECT = &
EIG_CORRECT_3(KK, M_eig_2, EIG_R_mat_2, EIG_L_mat_2, VECR_01, M_2, K_2, BETA, ALPHA, newcomm, .FALSE.)
VECR_1_CORRECT = VECR_1_CORRECT*BETA
VECR_2_CORRECT = VECR_2_CORRECT*ALPHA
EIG_R_mat_1(:, JJ) = EIG_R_mat_1(:, JJ)*ALPHA
EIG_R_mat_2(:, KK) = EIG_R_mat_2(:, KK)*BETA

! ===================== CREATE PERTURBATION FILES ======================
CALL SAVE_PERTURB('./data/new_perturb.input', M_1, K_1, EIG_R_mat_1(:, JJ), M_eig_1(JJ), &
    M_2, K_2, EIG_R_mat_2(:, KK), M_eig_2(kk), M_0, K_0, EIG_R_mat_0(:, II), M_eig_0(II))
CALL SAVE_PERTURB('./data/cor_perturb.input', M_1, K_1, VECR_1_CORRECT, M_eig_1(JJ), &
    M_2, K_2, VECR_2_CORRECT, M_eig_2(kk)) 

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
    DEALLOCATE(VECR_1_CORRECT, VECR_2_CORRECT)

END IF RANK0

! DEALLOCATE
!=======================================================================
10  CONTINUE
DEALLOCATE (RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
IF (ALLOCATED(M_mat_1)) THEN
    DEALLOCATE(M_mat_0,M_mat_1,M_mat_2)
ENDIF

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================

FUNCTION EIG_CORRECT_3(EIG_IND_DEGEN,EIG_VAL,EIG_VEC_R,EIG_VEC_L,VECR_B_DEGEN,MREAD,AKREAD,factor_self,factor_other,comm_grp,save_switch_in)
! ======================================================================
! NOTE: VECR_B_DEGEN AND EIG_VEC_R CORRESPOND TO THE TWO WAVENUMBERS
!       e.g. VECR_B_DEGEN = N(V_BAR)V2, EIG_VEC_R/L & EIG_IND_DEGEN & EIG_VAL ~ V1
! ======================================================================
    IMPLICIT NONE
    INTEGER:: EIG_IND_DEGEN, MREAD, comm_grp
    REAL(P8):: AKREAD,AK_ACTUAL
    COMPLEX(P8), DIMENSION(:):: EIG_VAL, VECR_B_DEGEN
    COMPLEX(P8), DIMENSION(:, :):: EIG_VEC_R, EIG_VEC_L
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_CORRECT_3
    LOGICAL, OPTIONAL:: save_switch_in
    LOGICAL:: save_switch, flip_switch

    INTEGER:: IND, EIG_NUM, NR_MK, M_ACTUAL
    COMPLEX(P8):: EIG_VAL_DEGEN, COEFF, factor_self, factor_other
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: CORRECTION, EIG_LEFTOVER
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_new, RUP_new, UZ_new
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_old, RUP_old, UZ_old
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSI, DEL2CHI, CHI

    IF ((PRESENT(save_switch_in)) .AND. (save_switch_in)) THEN
        save_switch = .TRUE.
    ELSE
        save_switch = .FALSE.
    END IF

    IF (MREAD.LT.0) THEN
        flip_switch = .TRUE.
        M_ACTUAL = -MREAD; AK_ACTUAL = -AKREAD
    ELSE
        flip_switch = .FALSE.
        M_ACTUAL = MREAD; AK_ACTUAL = AKREAD
    ENDIF
    IF ((M(2) .NE. M_ACTUAL) .OR. (ABS(ZLEN - 2*PI/AK_ACTUAL) .GE. 1.0E-13)) THEN
        WRITE (*, *) 'EIG_CORRECT_3: SET NEW M AND AK'
        ZLEN = 2*PI/AK_ACTUAL
        NTH = 2
        NX = 4
        NTCHOP = 2 ! M = 0, M = MREAD
        NXCHOP = 2 ! K = 0, K = AKREAD
        CALL LEGINIT(comm_grp, M_ACTUAL)
    ENDIF

    EIG_VAL_DEGEN = EIG_VAL(EIG_IND_DEGEN)
    ALLOCATE (EIG_CORRECT_3(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (CORRECTION(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (EIG_LEFTOVER(SIZE(EIG_VEC_R, 1)))

    EIG_CORRECT_3 = 0.D0
    EIG_LEFTOVER = 0.D0
    EIG_NUM = SIZE(EIG_VEC_R, 1)
    NR_MK = SIZE(EIG_VEC_R, 2)/2

    DO IND = 1, SIZE(EIG_VEC_R, 2)
        IF (IND .EQ. EIG_IND_DEGEN) CYCLE
        ! IF (SUM(ABS(EIG_VEC_L(:, IND))) .GT. 10E50) CYCLE ! TEMPORARY - FOR EIGVAL = 0 MODES
        IF (EIG_VAL(IND).EQ.0.D0) THEN
            IF (ISZERO(EIG_VEC_R,IND,NR_MK,flip_switch)) THEN
                WRITE(*,*) 'ISZERO:',IND; CYCLE
            ENDIF
        ENDIF

        ! calculate coefficient
        COEFF = DOT_PRODUCT(EIG_VEC_L(:, IND), VECR_B_DEGEN(:))&
        &/DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND))&
        &/(EIG_VAL_DEGEN - EIG_VAL(IND))

        ! debug
        ! write(*,*) IND,DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND)) !,EIG_VAL_DEGEN - EIG_VAL(IND)
        ! write(*,*) IND,SUM(ABS(EIG_VEC_R(2:, IND))),SUM(ABS(EIG_VEC_L(2:, IND)))

        ! obtain correction term
        CORRECTION = 0.D0
        CORRECTION = EIG2PTVEL(EIG_VEC_R,IND,NR_MK,COEFF,factor_other,flip_switch,save_switch,MREAD,AKREAD,EIG_VAL(IND))
        EIG_CORRECT_3 = EIG_CORRECT_3 + CORRECTION

        ! save left-over
        IF (IND .GT. 61) THEN
            EIG_LEFTOVER = EIG_LEFTOVER + factor_other*CORRECTION
        END IF     

    END DO

    ! save high-order terms
    IF (save_switch) THEN
        ALLOCATE (RUR_new(NRCHOPDIM), RUP_new(NRCHOPDIM), UZ_new(NRCHOPDIM))
        ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))
        RUR_new = 0.D0; RUP_new = 0.D0; UZ_new = 0.D0
        PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
        PSI(:NR_MK) = EIG_LEFTOVER(:NR_MK) !PSI
        DEL2CHI(:NR_MK) = EIG_LEFTOVER(NR_MK + 1:) !DEL2CHI
        IF (flip_switch) THEN
            PSI = CONJG(PSI); DEL2CHI = CONJG(DEL2CHI)
        ENDIF
        CALL IDEL2_MK(DEL2CHI, CHI)
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_new, RUP_new, UZ_new)
        CALL RTRAN_MK(RUR_new, 1); CALL RTRAN_MK(RUP_new, 1); CALL RTRAN_MK(UZ_new, 1)
        CALL CHOPSET(-3)
        IF (flip_switch) THEN
            RUR_new = CONJG(RUR_new); RUP_new = CONJG(RUP_new); UZ_new = CONJG(UZ_new)
        ENDIF
        CALL SAVE_VEL(RUR_new, RUP_new, UZ_new, './converg/check/lft_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
                        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
        DEALLOCATE(PSI, DEL2CHI, CHI, RUR_new, RUP_new, UZ_new)
    END IF

    DEALLOCATE(CORRECTION, EIG_LEFTOVER)

END FUNCTION EIG_CORRECT_3
! ======================================================================

FUNCTION EIG2PTVEL(EIG_VEC_R,IND,NR_MK,coeff,factor,flip_switch,save_switch,MREAD,AKREAD,eig_val)
! ======================================================================
! SUBROUTINE FOR EIG_CORRECT_3
! CONVERT EIG_VEC_R(:,IND) TO PSI,CHI
! IF SAVE_SWITH == TRUE:
!   SAVE THE CORRESPONDING VELOCITY FIELDS
!   IF FACTOR IS PRESENTED, SAVE THE FACTORED VELOCITY FIELDS AS WELL
! ======================================================================
    IMPLICIT NONE
    COMPLEX(P8):: coeff, factor, eig_val
    INTEGER:: IND, NR_MK, MREAD
    REAL(P8):: AKREAD
    LOGICAL:: flip_switch, save_switch
    COMPLEX(P8), DIMENSION(:,:):: EIG_VEC_R
    COMPLEX(P8):: EIG2PTVEL(2*NR_MK)

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSI,CHI,DEL2CHI
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_new, RUP_new, UZ_new
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_old, RUP_old, UZ_old
      
    ! convert to PSI,CHI
    ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))  
    PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
    PSI(:NR_MK) = EIG_VEC_R(:NR_MK, IND) !PSI
    DEL2CHI(:NR_MK) = EIG_VEC_R(NR_MK + 1:, IND) !DEL2CHI
    IF (flip_switch) THEN
        PSI = CONJG(PSI); DEL2CHI = CONJG(DEL2CHI)
    ENDIF
    CALL IDEL2_MK(DEL2CHI, CHI)

    ! save original
    IF (save_switch) THEN
        ALLOCATE (RUR_old(NRCHOPDIM), RUP_old(NRCHOPDIM), UZ_old(NRCHOPDIM))
        RUR_old = 0.D0; RUP_old = 0.D0; UZ_old = 0.D0
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_old, RUP_old, UZ_old)
        CALL RTRAN_MK(RUR_old, 1); CALL RTRAN_MK(RUP_old, 1); CALL RTRAN_MK(UZ_old, 1)
        CALL CHOPSET(-3)
        IF (flip_switch) THEN
            RUR_old = CONJG(RUR_old); RUP_old = CONJG(RUP_old); UZ_old = CONJG(UZ_old)
        ENDIF
        IF (IND .LE. 61) CALL SAVE_VEL(RUR_old, RUP_old, UZ_old, './converg/check/before_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
            //'_IND_'//ITOA3(IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ENDIF

    ! times coefficient
    EIG2PTVEL(:NR_MK) = PSI(:NR_MK)
    EIG2PTVEL(NR_MK + 1:) = CHI(:NR_MK)
    IF (flip_switch) THEN
        EIG2PTVEL = coeff*CONJG(EIG2PTVEL)
    ELSE
        EIG2PTVEL = coeff*EIG2PTVEL
    ENDIF

    ! save new
    IF (save_switch) THEN
        ALLOCATE (RUR_new(NRCHOPDIM), RUP_new(NRCHOPDIM), UZ_new(NRCHOPDIM))
        RUR_new = 0.D0; RUP_new = 0.D0; UZ_new = 0.D0
        IF (flip_switch) THEN
            PSI(:) = CONJG(factor*coeff)*PSI(:)
            CHI(:) = CONJG(factor*coeff)*CHI(:)
        ELSE
            PSI(:) = factor*coeff*PSI(:)
            CHI(:) = factor*coeff*CHI(:)
        ENDIF
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_new, RUP_new, UZ_new)
        CALL RTRAN_MK(RUR_new, 1); CALL RTRAN_MK(RUP_new, 1); CALL RTRAN_MK(UZ_new, 1)
        CALL CHOPSET(-3)
        IF (flip_switch) THEN
            RUR_new = CONJG(RUR_new); RUP_new = CONJG(RUP_new); UZ_new = CONJG(UZ_new)
        ENDIF
        IF (IND .LE. 61) CALL SAVE_VEL(RUR_new, RUP_new, UZ_new, './converg/check/after_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
            //'_IND_'//ITOA3(IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    END IF

    ! print data
    IF (save_switch) THEN
        open (10, FILE='./converg/check/contribution_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
            //'.output', STATUS='unknown', ACTION='WRITE', IOSTAT=IS)
        if (IS .ne. 0) then
            print *, 'ERROR: createperturb -- Could not creat new file contribution.output'
            RETURN
        end if
        WRITE (10, 1345) IND, ABS(EIG_VAL), ABS(factor*coeff),&
        &MAXVAL(ABS(RUP_old/TFM%R(1:NDIMR))),&
        &MAXVAL(ABS(RUP_new/TFM%R(1:NDIMR))),&
        &MAXVAL(ABS(EIG_VEC_R(:, IND))),&
        &MAXVAL(ABS(factor*COEFF*EIG_VEC_R(:, IND))),&
        &MAXVAL(ABS(RUP_old)),&
        &MAXVAL(ABS(RUP_new))
        close(10)
        DEALLOCATE(RUR_old, RUP_old, UZ_old)
        DEALLOCATE(RUR_new,RUP_new,UZ_new)
    END IF
    1345    FORMAT(I4, 8(',', G20.12))

    DEALLOCATE(PSI,CHI,DEL2CHI)

END FUNCTION EIG2PTVEL
! ======================================================================

FUNCTION ISZERO(EIG_VEC_R,IND,NR_MK,flip_switch)
! ======================================================================
! CHECK IF AN ERIGENVECTOR IS PRURELY ZERO
! ======================================================================
    IMPLICIT NONE
    INTEGER:: IND, NR_MK
    LOGICAL:: flip_switch, ISZERO
    COMPLEX(P8), DIMENSION(:,:):: EIG_VEC_R

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSI,CHI,DEL2CHI
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_new, RUP_new, UZ_new
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_old, RUP_old, UZ_old

    ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))  
    PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
    PSI(:NR_MK) = EIG_VEC_R(:NR_MK, IND) !PSI
    DEL2CHI(:NR_MK) = EIG_VEC_R(NR_MK + 1:, IND) !DEL2CHI
    IF (flip_switch) THEN
        PSI = CONJG(PSI); DEL2CHI = CONJG(DEL2CHI)
    ENDIF
    CALL IDEL2_MK(DEL2CHI, CHI)

    ! CALL RTRAN_MK(PSI,1); CALL RTRAN_MK(CHI,1);
    ! IF (SUM(ABS(PSI(1:NR)))+SUM(ABS(CHI(1:NR))) .EQ. 0.D0) THEN
    !     ISZERO = .TRUE.
    ! ELSE
    !     ISZERO = .FALSE.
    !     WRITE(*,*) 'PSI:',IND,PSI(1:NR)
    !     WRITE(*,*) 'CHI:',CHI(1:NR)
    ! ENDIF

    ALLOCATE (RUR_old(NRCHOPDIM), RUP_old(NRCHOPDIM), UZ_old(NRCHOPDIM))
    RUR_old = 0.D0; RUP_old = 0.D0; UZ_old = 0.D0
    CALL CHOPSET(3)
    CALL PC2VEL_MK(PSI, CHI, RUR_old, RUP_old, UZ_old)
    CALL RTRAN_MK(RUR_old, 1); CALL RTRAN_MK(RUP_old, 1); CALL RTRAN_MK(UZ_old, 1)
    CALL CHOPSET(-3)
    IF (flip_switch) THEN
        RUR_old = CONJG(RUR_old); RUP_old = CONJG(RUP_old); UZ_old = CONJG(UZ_old)
    ENDIF
    IF (SUM(ABS(RUR_old(1:NR)))+SUM(ABS(RUP_old(1:NR)))+SUM(ABS(UZ_old(1:NR))) .EQ. 0.D0) THEN
        ISZERO = .TRUE.
    ELSE
        ISZERO = .FALSE.
        WRITE(*,*) IND,SUM(ABS(RUR_old(1:NR)))+SUM(ABS(RUP_old(1:NR))),SUM(ABS(UZ_old(1:NR)))
    ENDIF
    DEALLOCATE(RUR_old,RUP_old,UZ_old)


END FUNCTION ISZERO
! ======================================================================


END PROGRAM EVP_PARAMETRIC
!=======================================================================
