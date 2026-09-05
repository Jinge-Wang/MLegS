PROGRAM EVP_TUNING
!=======================================================================
! [USAGE]:
! TUNE RESONANT TRIADS BY FINDING THE ROOT OF DELTA_OMEGA = 0
! WHERE DELTA_OMEGA = OMEGA_0 + OMEGA_1 - OMEGA_2
!
! USES NEWTON-RAPHSON METHOD WITH EIGENVALUE PERTURBATIONS:
!   delta_k = -delta_omega / (dsigma_1/dk - dsigma_2/dk)
!
! RESONANCE CONDITIONS:
!   m_2 = m_0 + m_1  (azimuthal resonance)
!   k_2 = k_0 + k_1  (axial resonance)
!   omega_2 = omega_0 + omega_1  (frequency resonance)
!
! PARAMETERS:
!   m0=<int>        : Azimuthal wavenumber of mode 0
!   k0=<float>      : Axial wavenumber of mode 0
!   omega0=<float>  : Target imaginary eigenvalue of mode 0
!   m1=<int>        : Azimuthal wavenumber of mode 1
!   k1=<float>      : Initial guess for k1 (default: k0)
!   omega1=<float>  : Initial guess for omega1 (optional)
!   nr=<int>        : Number of radial modes
!   tol=<float>     : Convergence tolerance (default: 1e-14)
!   maxiter=<int>   : Maximum iterations (default: 50)
!
! OUTPUT: ./data/tune_m0_<m0>_m1_<m1>_k0_<k0>.output
!   Contains: converged k1, eigenvalue indices, and final delta_omega
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
INTEGER     :: II, JJ, KK, ITER
INTEGER     :: M0, M1, M2, NR_MK, IDX0, IDX1, IDX2
INTEGER     :: MAX_ITER
REAL(P8)    :: K0, K1, K2, K1_PERT, K2_PERT, DELTA_K
REAL(P8)    :: OMEGA0_TARGET, OMEGA1_GUESS, TOLERANCE
REAL(P8)    :: DELTA_OMEGA, DSIGMA1_DK, DSIGMA2_DK, DELTA_K_UPDATE
LOGICAL     :: OMEGA1_PROVIDED, CONVERGED
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG0, EIG1, EIG2, EIG1_PERT, EIG2_PERT
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: MAT0, MAT1, MAT2, MAT1_PERT, MAT2_PERT
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: VEC_R0, VEC_L0, VEC_R1, VEC_L1, VEC_R2, VEC_L2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VEC_R1_SAVE, VEC_L1_SAVE, VEC_R2_SAVE, VEC_L2_SAVE
COMPLEX(P8) :: SIGMA0, SIGMA1, SIGMA2, PERTURB1, PERTURB2, NORM_LR

! SCANNING:
CHARACTER(len=6) :: K_VAL, H_VAL
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID, num_args, eq_pos, arg_len, stat
character(len=12), dimension(:), allocatable :: args
CHARACTER(LEN=32) :: arg_buffer
CHARACTER(LEN=16) :: param_name
CHARACTER(LEN=16) :: param_value
CHARACTER(LEN=256) :: dir_path
REAL(P8) :: temp_real
LOGICAL :: file_exists

! MPI INITILIZATION
CALL SETUP_ENVIRONMENT('NOECHO')
! PARALLEL COMPUTATION
! - SPLIT THE FORMATION OF H INTO DIFFERENT PROCS
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS

! RESONANT TRIAD TUNING PARAMETERS
! - Default values
M0 = 1
M1 = 1
K0 = 1.0d0
K1 = 0.0d0  ! Will default to K0 if not specified
OMEGA0_TARGET = 0.0d0
OMEGA1_GUESS = 0.0d0
OMEGA1_PROVIDED = .FALSE.
TOLERANCE = 1.0d-14
MAX_ITER = 10
DELTA_K = 0.0d0  ! Will be set to 0.01% of k if not specified
BSNSQ%BV0 = 0.D0
BSNSQ%OMEGA = 0.D0

! - READ FROM COMMAND INPUT
num_args = command_argument_count()
IF (num_args .ge. 1) THEN
    ! Process named arguments in param=value format
    DO II = 1, num_args
        CALL get_command_argument(II, arg_buffer)
        arg_len = LEN_TRIM(arg_buffer)        
        eq_pos = INDEX(arg_buffer, '=')
        
        ! Skip if no equals sign or it's at beginning/end
        IF (eq_pos <= 1 .OR. eq_pos >= arg_len) THEN
            IF (MPI_GLB_RANK .EQ. 0) THEN
                WRITE(*,*) 'Warning: Skipping invalid argument format:', TRIM(arg_buffer)
                WRITE(*,*) 'Expected format: param=value'
            END IF
            CYCLE
        END IF
        
        ! Extract parameter name and value
        param_name = arg_buffer(1:eq_pos-1)
        param_value = arg_buffer(eq_pos+1:arg_len)        
        ! Convert parameter name to lowercase for case-insensitive comparison
        CALL lowercase(param_name)
        SELECT CASE(TRIM(param_name))
            CASE('m0')
                READ(param_value, *, IOSTAT=stat) temp_real
                IF (stat == 0) THEN
                    M0 = NINT(temp_real)
                ELSE
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading m0 parameter:', TRIM(param_value)
                END IF
                
            CASE('m1')
                READ(param_value, *, IOSTAT=stat) temp_real
                IF (stat == 0) THEN
                    M1 = NINT(temp_real)
                ELSE
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading m1 parameter:', TRIM(param_value)
                END IF
                
            CASE('k0')
                READ(param_value, *, IOSTAT=stat) K0
                IF (stat /= 0 .AND. MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading k0 parameter:', TRIM(param_value)
                
            CASE('k1')
                READ(param_value, *, IOSTAT=stat) K1
                IF (stat /= 0 .AND. MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading k1 parameter:', TRIM(param_value)
                
            CASE('omega0')
                READ(param_value, *, IOSTAT=stat) OMEGA0_TARGET
                IF (stat /= 0 .AND. MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading omega0 parameter:', TRIM(param_value)
                
            CASE('omega1')
                READ(param_value, *, IOSTAT=stat) OMEGA1_GUESS
                IF (stat == 0) THEN
                    OMEGA1_PROVIDED = .TRUE.
                ELSE
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading omega1 parameter:', TRIM(param_value)
                END IF
                
            CASE('nr')
                READ(param_value, *, IOSTAT=stat) temp_real
                IF (stat == 0) THEN
                    NRCHOP = NINT(temp_real)
                ELSE
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading nr parameter:', TRIM(param_value)
                END IF
                
            CASE('tol', 'tolerance')
                READ(param_value, *, IOSTAT=stat) TOLERANCE
                IF (stat /= 0 .AND. MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading tol parameter:', TRIM(param_value)
                
            CASE('maxiter')
                READ(param_value, *, IOSTAT=stat) temp_real
                IF (stat == 0) THEN
                    MAX_ITER = NINT(temp_real)
                ELSE
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading maxiter parameter:', TRIM(param_value)
                END IF
                
            CASE('dk', 'deltak')
                READ(param_value, *, IOSTAT=stat) DELTA_K
                IF (stat /= 0 .AND. MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Error reading dk parameter:', TRIM(param_value)
                
            CASE DEFAULT
                IF (MPI_GLB_RANK .EQ. 0) THEN
                    WRITE(*,*) 'Unknown parameter:', TRIM(param_name)
                    WRITE(*,*) 'Supported parameters: m0, m1, k0, k1, omega0, omega1, nr, tol, maxiter, dk'
                END IF
        END SELECT
    END DO
    
    ! Set defaults for unspecified parameters
    IF (K1 .EQ. 0.d0) K1 = K0
    IF (DELTA_K .EQ. 0.d0) DELTA_K = 0.0001d0 * MAX(ABS(K0), ABS(K1))
    
    ! Calculate derived values
    M2 = M0 + M1
    K2 = K0 + K1
    
    ! Display parameters
    IF (MPI_GLB_RANK .EQ. 0) THEN
        WRITE(*,*) '======================================================'
        WRITE(*,*) 'RESONANT TRIAD TUNING'
        WRITE(*,*) '======================================================'
        WRITE(*,*) 'Mode 0: m0 =', M0, ', k0 =', K0, ', omega0_target =', OMEGA0_TARGET
        WRITE(*,*) 'Mode 1: m1 =', M1, ', k1_init =', K1
        IF (OMEGA1_PROVIDED) WRITE(*,*) '        omega1_guess =', OMEGA1_GUESS
        WRITE(*,*) 'Mode 2: m2 =', M2, ', k2_init =', K2
        WRITE(*,*) 'Numerical parameters:'
        WRITE(*,*) '  nr =', NRCHOP
        WRITE(*,*) '  dk =', DELTA_K
        WRITE(*,*) '  tolerance =', TOLERANCE
        WRITE(*,*) '  max_iterations =', MAX_ITER
        WRITE(*,*) '======================================================'
    END IF
ENDIF  

! ======================================================================
! STEP 1: FIND MODE 0 EIGENVALUE CLOSEST TO OMEGA0_TARGET
! ======================================================================
IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Step 1: Finding mode 0 eigenvalue...'
CALL EIG_MATRIX(M0, K0, MAT0, EIG0, EIG_VEC_R = VEC_R0, EIG_VEC_L = VEC_L0, comm_grp=newcomm, serial_switch=.false.)
NR_MK = NRCHOPS(2)

! Find closest resolved eigenvalue to OMEGA0_TARGET
IF (MPI_GLB_RANK .EQ. 0) THEN
    IDX0 = -1
    temp_real = HUGE(0.D0)
    DO II = 1, SIZE(EIG0)
        IF (EIGRES(VEC_R0(:,II), M0)) THEN  ! Check if resolved
            IF (ABS(AIMAG(EIG0(II)) - OMEGA0_TARGET) < temp_real) THEN
                temp_real = ABS(AIMAG(EIG0(II)) - OMEGA0_TARGET)
                IDX0 = II
            END IF
        END IF
    END DO
    
    IF (IDX0 .EQ. -1) THEN
        WRITE(*,*) 'ERROR: No resolved eigenvalue found for mode 0 near omega0 =', OMEGA0_TARGET
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    END IF
    
    SIGMA0 = EIG0(IDX0)
    WRITE(*,*) 'Mode 0: idx =', IDX0, ', sigma =', SIGMA0, ', resolved =', EIGRES(VEC_R0(:,IDX0), M0)
END IF

CALL MPI_BCAST(IDX0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
CALL MPI_BCAST(SIGMA0, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, IERR)

! ======================================================================
! STEP 2: INITIALIZE OUTPUT FILE
! ======================================================================
IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE(FILENAME,'(A)') './data/tune.log'
    
    ! Check if file exists to determine if we need to write header
    INQUIRE(FILE=TRIM(FILENAME), EXIST=file_exists)
    
    IF (file_exists) THEN
        OPEN(FID, FILE=TRIM(FILENAME), STATUS='old', POSITION='append', ACTION='WRITE', IOSTAT=IERR)
        WRITE(FID,*) '# ===== NEW RUN ====='
    ELSE
        OPEN(FID, FILE=TRIM(FILENAME), STATUS='new', ACTION='WRITE', IOSTAT=IERR)
        WRITE(FID,*) '# Resonant triad tuning convergence history'
        WRITE(FID,*) '# Columns: iter, m0, k0, idx0, m1, k1, idx1, m2, k2, idx2, sigma0_real, sigma0_imag, sigma1_real, sigma1_imag, sigma2_real, sigma2_imag, delta_omega, converged'
    END IF
    CLOSE(FID)
END IF

! ======================================================================
! STEP 3: ITERATE TO FIND RESONANT TRIAD
! ======================================================================
CONVERGED = .FALSE.
ITER = 0

DO WHILE ((.NOT. CONVERGED) .AND. (ITER < MAX_ITER))
    ITER = ITER + 1
    K2 = K0 + K1
    
    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) ''
    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Iteration', ITER, ': k1 =', K1, ', k2 =', K2
    
    ! Compute eigenvalues for modes 1 and 2
    CALL EIG_MATRIX(M1, K1, MAT1, EIG1, EIG_VEC_R = VEC_R1, EIG_VEC_L = VEC_L1, comm_grp=newcomm, serial_switch=.false.)
    CALL EIG_MATRIX(M2, K2, MAT2, EIG2, EIG_VEC_R = VEC_R2, EIG_VEC_L = VEC_L2, comm_grp=newcomm, serial_switch=.false.)
    
    ! Find best eigenvalue pair
    IF (MPI_GLB_RANK .EQ. 0) THEN
        temp_real = HUGE(0.D0)
        
        ! For iterations after the first, use predicted omega values to track branches
        ! BUT STILL ONLY USE RESOLVED MODES
        IF (OMEGA1_PROVIDED .AND. ITER .GT. 1) THEN
            ! Find RESOLVED mode 1 eigenvalue closest to predicted OMEGA1_GUESS
            IDX1 = -1
            DO II = 1, SIZE(EIG1)
                IF (EIGRES(VEC_R1(:,II), M1)) THEN
                    IF (ABS(AIMAG(EIG1(II)) - OMEGA1_GUESS) < temp_real) THEN
                        temp_real = ABS(AIMAG(EIG1(II)) - OMEGA1_GUESS)
                        IDX1 = II
                    END IF
                END IF
            END DO
            
            IF (IDX1 .EQ. -1) THEN
                WRITE(*,*) 'WARNING: No resolved eigenvalue found for mode 1 near predicted omega1 =', OMEGA1_GUESS
                WRITE(*,*) 'Falling back to best pair search...'
            ELSE
                SIGMA1 = EIG1(IDX1)
                ! Find RESOLVED mode 2 eigenvalue closest to predicted OMEGA2_GUESS
                temp_real = HUGE(0.D0)
                IDX2 = -1
                DO JJ = 1, SIZE(EIG2)
                    IF (EIGRES(VEC_R2(:,JJ), M2)) THEN
                        IF (ABS(AIMAG(EIG2(JJ)) - OMEGA1_GUESS - AIMAG(SIGMA0)) < temp_real) THEN
                            temp_real = ABS(AIMAG(EIG2(JJ)) - OMEGA1_GUESS - AIMAG(SIGMA0))
                            IDX2 = JJ
                        END IF
                    END IF
                END DO
                
                IF (IDX2 .EQ. -1) THEN
                    WRITE(*,*) 'WARNING: No resolved eigenvalue found for mode 2 near predicted omega2'
                    WRITE(*,*) 'Falling back to best pair search...'
                    IDX1 = -1  ! Reset so we do full search
                END IF
            END IF
        END IF
        
        ! For first iteration or if prediction failed, find best pair minimizing |delta_omega|
        ! ONLY USE RESOLVED MODES
        IF (.NOT. OMEGA1_PROVIDED .OR. IDX1 .EQ. -1 .OR. ITER .EQ. 1) THEN
            temp_real = HUGE(0.D0)
            IDX1 = -1
            IDX2 = -1
            DO II = 1, SIZE(EIG1)
                IF (EIGRES(VEC_R1(:,II), M1)) THEN
                    DO JJ = 1, SIZE(EIG2)
                        IF (EIGRES(VEC_R2(:,JJ), M2)) THEN
                            DELTA_OMEGA = ABS(AIMAG(SIGMA0) + AIMAG(EIG1(II)) - AIMAG(EIG2(JJ)))
                            IF (DELTA_OMEGA < temp_real) THEN
                                temp_real = DELTA_OMEGA
                                IDX1 = II
                                IDX2 = JJ
                            END IF
                        END IF
                    END DO
                END IF
            END DO
            
            ! After first iteration, set OMEGA1_PROVIDED flag for branch tracking
            IF (ITER .EQ. 1) OMEGA1_PROVIDED = .TRUE.
        END IF
        
        IF (IDX1 .EQ. -1 .OR. IDX2 .EQ. -1) THEN
            WRITE(*,*) 'ERROR: No resolved eigenvalue pair found for modes 1 and 2'
            CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        END IF
        
        SIGMA1 = EIG1(IDX1)
        SIGMA2 = EIG2(IDX2)
        DELTA_OMEGA = AIMAG(SIGMA0) + AIMAG(SIGMA1) - AIMAG(SIGMA2)
        
        WRITE(*,*) 'Mode 1: idx =', IDX1, ', sigma =', SIGMA1, ', resolved =', EIGRES(VEC_R1(:,IDX1), M1)
        WRITE(*,*) 'Mode 2: idx =', IDX2, ', sigma =', SIGMA2, ', resolved =', EIGRES(VEC_R2(:,IDX2), M2)
        WRITE(*,*) 'delta_omega =', DELTA_OMEGA
        
        ! Write iteration data to file
        OPEN(FID, FILE=TRIM(FILENAME), STATUS='old', POSITION='append', ACTION='WRITE', IOSTAT=IERR)
        WRITE(FID,'(I0,1X,I0,1X,ES23.16,1X,I0,1X,I0,1X,ES23.16,1X,I0,1X,I0,1X,ES23.16,1X,I0,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,L1)') &
            ITER, M0, K0, IDX0, M1, K1, IDX1, M2, K2, IDX2, &
            REAL(SIGMA0), AIMAG(SIGMA0), REAL(SIGMA1), AIMAG(SIGMA1), &
            REAL(SIGMA2), AIMAG(SIGMA2), DELTA_OMEGA, .FALSE.
        CLOSE(FID)
        
        ! Check convergence
        IF (ABS(DELTA_OMEGA) < TOLERANCE) THEN
            CONVERGED = .TRUE.
            WRITE(*,*) 'CONVERGED!'
            
            ! Write final converged iteration data to file
            OPEN(FID, FILE=TRIM(FILENAME), STATUS='old', POSITION='append', ACTION='WRITE', IOSTAT=IERR)
            WRITE(FID,'(I0,1X,I0,1X,ES23.16,1X,I0,1X,I0,1X,ES23.16,1X,I0,1X,I0,1X,ES23.16,1X,I0,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,ES23.16,1X,L1)') &
                ITER, M0, K0, IDX0, M1, K1, IDX1, M2, K2, IDX2, &
                REAL(SIGMA0), AIMAG(SIGMA0), REAL(SIGMA1), AIMAG(SIGMA1), &
                REAL(SIGMA2), AIMAG(SIGMA2), DELTA_OMEGA, .TRUE.
            WRITE(FID,'(A,I0,A)') '# Converged after ', ITER, ' iterations'
            CLOSE(FID)
        END IF
    END IF
    
    CALL MPI_BCAST(CONVERGED, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)
    IF (CONVERGED) EXIT
    
    ! ======================================================================
    ! STEP 3: COMPUTE PERTURBATIONS AND UPDATE K1
    ! ======================================================================
    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Computing perturbations...'
    
    ! Compute perturbed eigenmatrices
    K1_PERT = K1 + DELTA_K
    K2_PERT = K2 + DELTA_K
    CALL EIG_MATRIX(M1, K1_PERT, MAT1_PERT, EIG1_PERT, comm_grp=newcomm, serial_switch=.false.)
    CALL EIG_MATRIX(M2, K2_PERT, MAT2_PERT, EIG2_PERT, comm_grp=newcomm, serial_switch=.false.)

    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Computing dsigma/dk...'
    
    ! Compute dsigma/dk for modes 1 and 2
    IF (MPI_GLB_RANK .EQ. 0) THEN
        ! dsigma1/dk using non-degenerate perturbation theory
        NORM_LR = DOT_PRODUCT(VEC_L1(:, IDX1), VEC_R1(:, IDX1))
        PERTURB1 = DOT_PRODUCT(VEC_L1(:, IDX1), &
                               MATMUL(MAT1_PERT(:, :), VEC_R1(:, IDX1)))
        DSIGMA1_DK = AIMAG((PERTURB1 / NORM_LR - SIGMA1) / DELTA_K)
        
        ! dsigma2/dk using non-degenerate perturbation theory
        NORM_LR = DOT_PRODUCT(VEC_L2(:, IDX2), VEC_R2(:, IDX2))
        PERTURB2 = DOT_PRODUCT(VEC_L2(:, IDX2), &
                               MATMUL(MAT2_PERT(:, :), VEC_R2(:, IDX2)))
        DSIGMA2_DK = AIMAG((PERTURB2 / NORM_LR - SIGMA2) / DELTA_K)
        
        WRITE(*,*) 'dsigma1/dk =', DSIGMA1_DK
        WRITE(*,*) 'dsigma2/dk =', DSIGMA2_DK
        
        ! Newton-Raphson update: delta_k = -delta_omega / (dsigma1/dk - dsigma2/dk)
        ! Note: k2 = k0 + k1, so dk2 = dk1
        DELTA_K_UPDATE = -DELTA_OMEGA / (DSIGMA1_DK - DSIGMA2_DK)
        
        ! ! Apply damping for stability
        ! IF (ABS(DELTA_K_UPDATE) > 0.5 * ABS(K1)) THEN
        !     DELTA_K_UPDATE = SIGN(0.5 * ABS(K1), DELTA_K_UPDATE)
        !     WRITE(*,*) 'Damping applied: delta_k limited to 50% of k1'
        ! END IF
        
        K1 = K1 + DELTA_K_UPDATE
        WRITE(*,*) 'delta_k =', DELTA_K_UPDATE, ', new k1 =', K1
        
        ! Predict omega values for next iteration (for branch tracking)
        OMEGA1_GUESS = AIMAG(SIGMA1) + DSIGMA1_DK * DELTA_K_UPDATE
        WRITE(*,*) 'Predicted omega1 for next iteration =', OMEGA1_GUESS
    END IF
    
    CALL MPI_BCAST(K1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(OMEGA1_GUESS, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(OMEGA1_PROVIDED, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)
    
    ! Clean up iteration matrices
    IF (ALLOCATED(MAT1)) DEALLOCATE(MAT1)
    IF (ALLOCATED(MAT2)) DEALLOCATE(MAT2)
    IF (ALLOCATED(MAT1_PERT)) DEALLOCATE(MAT1_PERT)
    IF (ALLOCATED(MAT2_PERT)) DEALLOCATE(MAT2_PERT)
    IF (ALLOCATED(EIG1)) DEALLOCATE(EIG1)
    IF (ALLOCATED(EIG2)) DEALLOCATE(EIG2)
    IF (ALLOCATED(EIG1_PERT)) DEALLOCATE(EIG1_PERT)
    IF (ALLOCATED(EIG2_PERT)) DEALLOCATE(EIG2_PERT)
    IF (ALLOCATED(VEC_R1)) DEALLOCATE(VEC_R1)
    IF (ALLOCATED(VEC_L1)) DEALLOCATE(VEC_L1)
    IF (ALLOCATED(VEC_R2)) DEALLOCATE(VEC_R2)
    IF (ALLOCATED(VEC_L2)) DEALLOCATE(VEC_L2)
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

! ======================================================================
! OUTPUT RESULTS
! ======================================================================
IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE(*,*) ''
    WRITE(*,*) '======================================================'
    IF (CONVERGED) THEN
        WRITE(*,*) 'TUNING CONVERGED IN', ITER, 'ITERATIONS'
    ELSE
        WRITE(*,*) 'WARNING: MAXIMUM ITERATIONS REACHED WITHOUT CONVERGENCE'
    END IF
    WRITE(*,*) '======================================================'
    WRITE(*,*) 'Final results:'
    WRITE(*,*) '  k1 =', K1
    WRITE(*,*) '  k2 =', K2
    WRITE(*,*) '  Mode 0: idx =', IDX0, ', sigma =', SIGMA0
    WRITE(*,*) '  Mode 1: idx =', IDX1, ', sigma =', SIGMA1
    WRITE(*,*) '  Mode 2: idx =', IDX2, ', sigma =', SIGMA2
    WRITE(*,*) '  delta_omega =', DELTA_OMEGA
    WRITE(*,*) '======================================================'
    
    ! Write final summary to file
    OPEN(FID, FILE=TRIM(FILENAME), STATUS='old', POSITION='append', ACTION='WRITE', IOSTAT=IERR)
    IF (.NOT. CONVERGED) THEN
        WRITE(FID,'(A,I0,A)') '# WARNING: Did not converge after ', MAX_ITER, ' iterations'
    END IF
    WRITE(FID,*) ''
    CLOSE(FID)
    
    WRITE(*,*) 'Convergence history saved to:', TRIM(FILENAME)
END IF

129 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

CONTAINS
  SUBROUTINE lowercase(str)
    CHARACTER(LEN=*), INTENT(INOUT) :: str
    INTEGER :: i, diff
    
    diff = IACHAR('a') - IACHAR('A')
    DO i = 1, LEN_TRIM(str)
      IF (str(i:i) >= 'A' .AND. str(i:i) <= 'Z') THEN
        str(i:i) = ACHAR(IACHAR(str(i:i)) + diff)
      END IF
    END DO
  END SUBROUTINE lowercase

END PROGRAM EVP_TUNING
!=======================================================================
