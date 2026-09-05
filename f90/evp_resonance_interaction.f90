PROGRAM EVP_RESONANCE_INTERACTION
!=======================================================================
! [USAGE]:
! Calculate interaction coefficients J0, J1, J2 for resonant triads
! from resonance_triads.dat file using MPI parallelization
!
! [DESCRIPTION]:
! - Root process reads resonance_triads.dat
! - For each unique {m0,k0}, assigns task to available MPI process
! - Worker processes:
!   * Solve EVP for {m0,k0}, {m1,k1}, {m2,k2}
!   * Find modes closest to listed eigenvalues (within tolerance)
!   * If no mode within tolerance, use closest mode with |Re(omega)| < 1e-14
!   * Calculate interaction coefficients J0, J1, J2
!   * Save results to process-specific file
! - Root process assembles all files at the end
!
! [COMMAND LINE ARGUMENTS]:
! m0=<value>  : Process only this azimuthal wavenumber
! k0=<value>  : Process only this axial wavenumber
! nr=<value>  : Number of radial grid points (default from input file)
!
! [BASED ON]:
! evp_print_org.f90 (current version structure)
! evp_three_k0_parallel.f90 (MPI task scheduling pattern)
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

! TRIAD DATA STRUCTURE
TYPE :: TRIAD_DATA
    INTEGER :: m0, m1, m2
    INTEGER :: n0, n1, n2
    REAL(P8) :: k0, k1, k2
    REAL(P8) :: omega0, omega1, omega2
END TYPE TRIAD_DATA

! EVP VARIABLES
INTEGER :: M_BAR, NR_MK
REAL(P8) :: K_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: M_eig_0, M_eig_1, M_eig_2
COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: M_mat_0, EIG_R_mat_0, EIG_L_mat_0
COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: M_mat_1, EIG_R_mat_1, EIG_L_mat_1
COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: M_mat_2, EIG_R_mat_2, EIG_L_mat_2

! VELOCITY AND VORTICITY FIELDS
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: VECL_0, VECR_0, VECR_01, VECR_12
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: ROR_0_COPY, ROP_0_COPY, OZ_0_COPY
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: VECL_1, VECR_1, VECR_02
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: ROR_1_COPY, ROP_1_COPY, OZ_1_COPY
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: VECL_2, VECR_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: ROR_2_COPY, ROP_2_COPY, OZ_2_COPY

! INTERACTION COEFFICIENTS
COMPLEX(P8) :: J0, J1, J2

! DATA MANAGEMENT
TYPE(TRIAD_DATA), DIMENSION(:), ALLOCATABLE :: TRIADS
TYPE(TRIAD_DATA) :: CURRENT_TRIAD
INTEGER :: NTRIADS, TRIAD_IDX
INTEGER, DIMENSION(:,:), ALLOCATABLE :: UNIQUE_M0K0  ! (m0, k0_idx) pairs
REAL(P8), DIMENSION(:), ALLOCATABLE :: UNIQUE_K0_VALS
INTEGER :: NUNIQUE, NUNIQUE_K0, TASK_IDX
INTEGER, DIMENSION(2) :: TASK_INFO  ! (m0, k0_idx)

! FILE I/O
CHARACTER(LEN=256) :: INPUT_FILE, OUTPUT_FILE, FINAL_OUTPUT
INTEGER :: FID, STAT, LINE_COUNT, II, JJ, KK
CHARACTER(LEN=512) :: LINE
LOGICAL :: FILE_EXISTS

! MPI SCHEDULING
INTEGER :: MPI_REC_ID, WORKER_STAT(MPI_STATUS_SIZE)
INTEGER :: REQUEST
LOGICAL :: WORKER

! COMMAND LINE ARGUMENTS
INTEGER :: num_args, eq_pos, arg_len
CHARACTER(LEN=256) :: arg_buffer
CHARACTER(LEN=32) :: param_name, param_value
REAL(P8) :: temp_real
INTEGER :: M0_FILTER, K0_IDX_FILTER
REAL(P8) :: K0_FILTER
LOGICAL :: USE_M0_FILTER, USE_K0_FILTER

! TIME LIMIT VARIABLES
REAL(P8) :: TIME_LIMIT_MINUTES, ELAPSED_TIME_MINUTES
LOGICAL :: USE_TIME_LIMIT
INTEGER :: TIME_START(8), TIME_NOW(8)
REAL(P8) :: ELAPSED_SECONDS

! EIGENVALUE SEARCH
REAL(P8), PARAMETER :: OMEGA_TOL_MAX = 1.0D-2  ! Max tolerance for frequency match
REAL(P8), PARAMETER :: REAL_TOL = 1.0D-14      ! Max |Re(omega)| for acceptable mode
INTEGER :: EIG_IDX_0, EIG_IDX_1, EIG_IDX_2
REAL(P8) :: OMEGA_ERR_0, OMEGA_ERR_1, OMEGA_ERR_2
LOGICAL :: FOUND_MODE_0, FOUND_MODE_1, FOUND_MODE_2
LOGICAL :: INTERACTION_FLAG

!=======================================================================
! MPI INITIALIZATION
!=======================================================================
CALL SETUP_ENVIRONMENT('NOECHO')
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR)

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE(*,*) '======================================================='
    WRITE(*,*) 'RESONANCE TRIAD INTERACTION COEFFICIENT CALCULATOR'
    WRITE(*,*) '======================================================='
    CALL PRINT_REAL_TIME()
ENDIF

!=======================================================================
! PARSE COMMAND LINE ARGUMENTS
!=======================================================================
USE_M0_FILTER = .FALSE.
USE_K0_FILTER = .FALSE.
USE_TIME_LIMIT = .FALSE.

num_args = command_argument_count()
IF (num_args .GE. 1) THEN
    DO II = 1, num_args
        CALL get_command_argument(II, arg_buffer)
        arg_len = LEN_TRIM(arg_buffer)
        eq_pos = INDEX(arg_buffer, '=')
        
        IF (eq_pos <= 1 .OR. eq_pos >= arg_len) CYCLE
        
        param_name = arg_buffer(1:eq_pos-1)
        param_value = arg_buffer(eq_pos+1:arg_len)
        CALL lowercase(param_name)
        
        SELECT CASE(TRIM(param_name))
            CASE('m0')
                READ(param_value, *, IOSTAT=STAT) temp_real
                IF (STAT == 0) THEN
                    M0_FILTER = NINT(temp_real)
                    USE_M0_FILTER = .TRUE.
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Filtering for m0 =', M0_FILTER
                ENDIF
                
            CASE('k0')
                READ(param_value, *, IOSTAT=STAT) K0_FILTER
                IF (STAT == 0) THEN
                    USE_K0_FILTER = .TRUE.
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Filtering for k0 =', K0_FILTER
                ENDIF
                
            CASE('nr')
                READ(param_value, *, IOSTAT=STAT) temp_real
                IF (STAT == 0) THEN
                    NRCHOP = NINT(temp_real)
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Using nr =', NRCHOP
                ENDIF
                
            CASE('time_limit')
                READ(param_value, *, IOSTAT=STAT) TIME_LIMIT_MINUTES
                IF (STAT == 0) THEN
                    USE_TIME_LIMIT = .TRUE.
                    IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Time limit =', TIME_LIMIT_MINUTES, 'minutes'
                ENDIF
        END SELECT
    END DO
ENDIF

!=======================================================================
! ROOT PROCESS: READ RESONANCE TRIAD DATA
!=======================================================================
IF (MPI_GLB_RANK .EQ. 0) THEN
    INPUT_FILE = './resonance_triads.dat'
    
    ! Count non-comment lines
    OPEN(UNIT=FID, FILE=TRIM(INPUT_FILE), STATUS='OLD', ACTION='READ', IOSTAT=STAT)
    IF (STAT /= 0) THEN
        WRITE(*,*) 'ERROR: Cannot open file ', TRIM(INPUT_FILE)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    ENDIF
    
    LINE_COUNT = 0
    DO
        READ(FID, '(A)', IOSTAT=STAT) LINE
        IF (STAT /= 0) EXIT
        IF (LEN_TRIM(LINE) == 0) CYCLE
        IF (LINE(1:1) == '#') CYCLE
        LINE_COUNT = LINE_COUNT + 1
    END DO
    REWIND(FID)
    
    WRITE(*,*) 'Found', LINE_COUNT, 'triad entries in', TRIM(INPUT_FILE)
    
    ! Read triad data
    ALLOCATE(TRIADS(LINE_COUNT))
    NTRIADS = 0
    
    DO
        READ(FID, '(A)', IOSTAT=STAT) LINE
        IF (STAT /= 0) EXIT
        IF (LEN_TRIM(LINE) == 0) CYCLE
        IF (LINE(1:1) == '#') CYCLE
        
        NTRIADS = NTRIADS + 1
        READ(LINE, *) TRIADS(NTRIADS)%m0, TRIADS(NTRIADS)%k0, TRIADS(NTRIADS)%n0, TRIADS(NTRIADS)%omega0, &
                      TRIADS(NTRIADS)%m1, TRIADS(NTRIADS)%k1, TRIADS(NTRIADS)%n1, TRIADS(NTRIADS)%omega1, &
                      TRIADS(NTRIADS)%m2, TRIADS(NTRIADS)%k2, TRIADS(NTRIADS)%n2, TRIADS(NTRIADS)%omega2
    END DO
    CLOSE(FID)
    
    WRITE(*,*) 'Loaded', NTRIADS, 'triads successfully'
    
    ! Extract unique (m0, k0) pairs
    CALL EXTRACT_UNIQUE_M0K0(TRIADS, NTRIADS, UNIQUE_M0K0, UNIQUE_K0_VALS, NUNIQUE, &
                              USE_M0_FILTER, M0_FILTER, USE_K0_FILTER, K0_FILTER, K0_IDX_FILTER)
    
    WRITE(*,*) 'Found', NUNIQUE, 'unique (m0, k0) pairs to process'
    IF (USE_M0_FILTER .OR. USE_K0_FILTER) THEN
        WRITE(*,*) 'After filtering'
    ENDIF
    
    NUNIQUE_K0 = SIZE(UNIQUE_K0_VALS)
    
    ! Open log file
    OPEN(UNIT=99, FILE='resonance_interaction.log', STATUS='unknown', ACTION='WRITE', IOSTAT=STAT)
    WRITE(99,*) 'Resonance Interaction Calculation Log'
    WRITE(99,*) 'Started at:'
    CALL PRINT_REAL_TIME_TO_UNIT(99)
    WRITE(99,*) 'Total unique (m0,k0) pairs:', NUNIQUE
    
    ! Broadcast size of unique_k0_vals array to all processes
    CALL MPI_BCAST(NUNIQUE_K0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    
    ! Broadcast unique k0 values to all processes
    CALL MPI_BCAST(UNIQUE_K0_VALS, NUNIQUE_K0, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    
ELSE
    ! Workers receive broadcast data
    CALL MPI_BCAST(NUNIQUE_K0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    
    ! Allocate on worker processes
    ALLOCATE(UNIQUE_K0_VALS(NUNIQUE_K0))
    
    ! Receive broadcast of unique k0 values
    CALL MPI_BCAST(UNIQUE_K0_VALS, NUNIQUE_K0, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
ENDIF

!=======================================================================
! MPI TASK DISTRIBUTION
!=======================================================================
IF (MPI_GLB_RANK .EQ. 0) THEN
    ! ==================== ROOT: TASK MANAGER ====================
    TASK_INFO = 0
    
    ! Record start time if using time limit
    IF (USE_TIME_LIMIT) THEN
        CALL DATE_AND_TIME(VALUES=TIME_START)
        WRITE(*,*) 'Starting timer - will terminate after', TIME_LIMIT_MINUTES, 'minutes'
    ENDIF
    
    DO TASK_IDX = 1, NUNIQUE
        ! Check time limit before assigning new task
        IF (USE_TIME_LIMIT) THEN
            CALL DATE_AND_TIME(VALUES=TIME_NOW)
            ! Calculate elapsed time in seconds
            ELAPSED_SECONDS = (TIME_NOW(5) - TIME_START(5)) * 3600.0D0 + &
                             (TIME_NOW(6) - TIME_START(6)) * 60.0D0 + &
                             (TIME_NOW(7) - TIME_START(7)) + &
                             (TIME_NOW(8) - TIME_START(8)) * 0.001D0
            ! Handle day rollover
            IF (TIME_NOW(3) /= TIME_START(3)) THEN
                ELAPSED_SECONDS = ELAPSED_SECONDS + 86400.0D0
            ENDIF
            ELAPSED_TIME_MINUTES = ELAPSED_SECONDS / 60.0D0
            
            IF (ELAPSED_TIME_MINUTES >= TIME_LIMIT_MINUTES) THEN
                WRITE(*,*) ''
                WRITE(*,*) 'Time limit reached (', ELAPSED_TIME_MINUTES, 'minutes)'
                WRITE(*,*) 'Completed', TASK_IDX-1, 'of', NUNIQUE, 'tasks'
                WRITE(*,*) 'Sending termination signal to all workers...'
                WRITE(99,*) 'Time limit reached after', ELAPSED_TIME_MINUTES, 'minutes'
                WRITE(99,*) 'Completed', TASK_IDX-1, 'of', NUNIQUE, 'tasks'
                EXIT  ! Exit task assignment loop
            ENDIF
        ENDIF
        
        ! Wait for available worker
        CALL MPI_RECV(WORKER, 1, MPI_LOGICAL, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                      MPI_COMM_WORLD, WORKER_STAT, IERR)
        MPI_REC_ID = WORKER_STAT(MPI_SOURCE)
        
        IF (WORKER) THEN
            ! Assign task to worker
            TASK_INFO(1) = UNIQUE_M0K0(TASK_IDX, 1)  ! m0
            TASK_INFO(2) = UNIQUE_M0K0(TASK_IDX, 2)  ! k0_idx
            
            CALL MPI_SEND(TASK_INFO, 2, MPI_INTEGER, MPI_REC_ID, 1, MPI_COMM_WORLD, IERR)
            
            ! Log assignment
            K0_FILTER = UNIQUE_K0_VALS(TASK_INFO(2))
            WRITE(99,'(A,I4,A,I3,A,F10.6,A,I4,A,I4,A)') &
                'Rank ', MPI_REC_ID, ': m0=', TASK_INFO(1), ', k0=', K0_FILTER, &
                ' (Task ', TASK_IDX, '/', NUNIQUE, ')'
            WRITE(*,'(A,I4,A,I3,A,F10.6,A,I4,A,I4,A)') &
                'Assigned to Rank ', MPI_REC_ID, ': m0=', TASK_INFO(1), ', k0=', K0_FILTER, &
                ' (', TASK_IDX, '/', NUNIQUE, ')'
        ENDIF
    END DO
    
    ! Notify all workers: job done
    IF (.NOT. USE_TIME_LIMIT .OR. TASK_IDX > NUNIQUE) THEN
        WRITE(*,*) 'All tasks assigned. Sending termination signal to workers...'
    ENDIF
    DO II = 1, MPI_GLB_PROCS - 1
        TASK_INFO(1) = -1
        CALL MPI_RECV(WORKER, 1, MPI_LOGICAL, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                      MPI_COMM_WORLD, WORKER_STAT, IERR)
        MPI_REC_ID = WORKER_STAT(MPI_SOURCE)
        CALL MPI_SEND(TASK_INFO, 2, MPI_INTEGER, MPI_REC_ID, 1, MPI_COMM_WORLD, IERR)
    END DO
    
    CLOSE(99)
    IF (USE_TIME_LIMIT .AND. TASK_IDX <= NUNIQUE) THEN
        WRITE(*,*) 'Root process: Terminated early due to time limit'
    ELSE
        WRITE(*,*) 'Root process: All tasks distributed'
    ENDIF
    
ELSE
    ! ==================== WORKER PROCESSES ====================
    
    ! Read all triads from file (each worker reads independently)
    INPUT_FILE = './resonance_triads.dat'
    OPEN(UNIT=FID, FILE=TRIM(INPUT_FILE), STATUS='OLD', ACTION='READ', IOSTAT=STAT)
    
    LINE_COUNT = 0
    DO
        READ(FID, '(A)', IOSTAT=STAT) LINE
        IF (STAT /= 0) EXIT
        IF (LEN_TRIM(LINE) == 0) CYCLE
        IF (LINE(1:1) == '#') CYCLE
        LINE_COUNT = LINE_COUNT + 1
    END DO
    REWIND(FID)
    
    ALLOCATE(TRIADS(LINE_COUNT))
    NTRIADS = 0
    
    DO
        READ(FID, '(A)', IOSTAT=STAT) LINE
        IF (STAT /= 0) EXIT
        IF (LEN_TRIM(LINE) == 0) CYCLE
        IF (LINE(1:1) == '#') CYCLE
        
        NTRIADS = NTRIADS + 1
        READ(LINE, *) TRIADS(NTRIADS)%m0, TRIADS(NTRIADS)%k0, TRIADS(NTRIADS)%n0, TRIADS(NTRIADS)%omega0, &
                      TRIADS(NTRIADS)%m1, TRIADS(NTRIADS)%k1, TRIADS(NTRIADS)%n1, TRIADS(NTRIADS)%omega1, &
                      TRIADS(NTRIADS)%m2, TRIADS(NTRIADS)%k2, TRIADS(NTRIADS)%n2, TRIADS(NTRIADS)%omega2
    END DO
    CLOSE(FID)
    
    ! Initialize task tracking
    TASK_INFO(1) = 0
    
    DO WHILE (TASK_INFO(1) /= -1)
        ! Signal availability to root
        WORKER = .TRUE.
        CALL MPI_SEND(WORKER, 1, MPI_LOGICAL, 0, 1, MPI_COMM_WORLD, IERR)
        
        ! Receive task
        CALL MPI_RECV(TASK_INFO, 2, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, WORKER_STAT, IERR)
        
        IF (TASK_INFO(1) == -1) EXIT  ! Termination signal
        
        ! Extract task parameters: m0 value and k0 index
        M0_FILTER = TASK_INFO(1)
        K0_IDX_FILTER = TASK_INFO(2)
        K0_FILTER = UNIQUE_K0_VALS(K0_IDX_FILTER)
        
        ! Initialize if first task
        IF (.NOT. ASSOCIATED(TFM%R)) THEN
            NTH = 2
            NX = 4
            NTCHOP = 2
            NXCHOP = 2
            ZLEN = 1.D0
            CALL LEGINIT(newcomm)
        ENDIF
        
        ! Open output file for this worker
        WRITE(OUTPUT_FILE, '(A,I0.4,A)') './data/interaction_rank_', MPI_GLB_RANK, '.dat'
        OPEN(UNIT=FID, FILE=TRIM(OUTPUT_FILE), STATUS='unknown', ACTION='WRITE', ACCESS='APPEND', IOSTAT=STAT)
        
        ! Solve EVP for (m0, k0) - only once per task
        CALL EIG_MATRIX(M0_FILTER, K0_FILTER, M_mat_0, M_eig_0, &
                        EIG_VEC_R=EIG_R_mat_0, EIG_VEC_L=EIG_L_mat_0, &
                        comm_grp=newcomm, serial_switch=.TRUE.)
        NR_MK = NRCHOPS(2)
        
        ! Count and process triads with this (m0, k0)
        TRIAD_IDX = 0
        DO II = 1, NTRIADS
            ! Check if this triad matches current (m0, k0)
            IF (TRIADS(II)%m0 /= M0_FILTER) CYCLE
            IF (ABS(TRIADS(II)%k0 - K0_FILTER) > 1.0D-10) CYCLE
            
            TRIAD_IDX = TRIAD_IDX + 1
            CURRENT_TRIAD = TRIADS(II)
            
            ! Solve EVP for (m1, k1) and (m2, k2)
            CALL EIG_MATRIX(CURRENT_TRIAD%m1, CURRENT_TRIAD%k1, M_mat_1, M_eig_1, &
                            EIG_VEC_R=EIG_R_mat_1, EIG_VEC_L=EIG_L_mat_1, &
                            comm_grp=newcomm, serial_switch=.TRUE.)
            
            CALL EIG_MATRIX(CURRENT_TRIAD%m2, CURRENT_TRIAD%k2, M_mat_2, M_eig_2, &
                            EIG_VEC_R=EIG_R_mat_2, EIG_VEC_L=EIG_L_mat_2, &
                            comm_grp=newcomm, serial_switch=.TRUE.)
            
            ! Find modes closest to target frequencies
            CALL FIND_MODE(M_eig_0, CURRENT_TRIAD%omega0, OMEGA_TOL_MAX, REAL_TOL, &
                          EIG_IDX_0, OMEGA_ERR_0, FOUND_MODE_0)
            CALL FIND_MODE(M_eig_1, CURRENT_TRIAD%omega1, OMEGA_TOL_MAX, REAL_TOL, &
                          EIG_IDX_1, OMEGA_ERR_1, FOUND_MODE_1)
            CALL FIND_MODE(M_eig_2, CURRENT_TRIAD%omega2, OMEGA_TOL_MAX, REAL_TOL, &
                          EIG_IDX_2, OMEGA_ERR_2, FOUND_MODE_2)
            
            IF (FOUND_MODE_0 .AND. FOUND_MODE_1 .AND. FOUND_MODE_2) THEN
                ! Calculate velocity and vorticity fields
                CALL EIG2VELVOR(CURRENT_TRIAD%m0, CURRENT_TRIAD%k0, EIG_IDX_0, &
                               EIG_R_mat_0, EIG_L_mat_0, &
                               RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0, &
                               VECR_0, VECL_0, comm_grp=newcomm)
                ALLOCATE(ROR_0_COPY(NDIMR), ROP_0_COPY(NDIMR), OZ_0_COPY(NDIMR))
                ROR_0_COPY = ROR_0; ROP_0_COPY = ROP_0; OZ_0_COPY = OZ_0
                
                CALL EIG2VELVOR(CURRENT_TRIAD%m1, CURRENT_TRIAD%k1, EIG_IDX_1, &
                               EIG_R_mat_1, EIG_L_mat_1, &
                               RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1, &
                               VECR_1, VECL_1, comm_grp=newcomm)
                ALLOCATE(ROR_1_COPY(NDIMR), ROP_1_COPY(NDIMR), OZ_1_COPY(NDIMR))
                ROR_1_COPY = ROR_1; ROP_1_COPY = ROP_1; OZ_1_COPY = OZ_1
                
                CALL EIG2VELVOR(CURRENT_TRIAD%m2, CURRENT_TRIAD%k2, EIG_IDX_2, &
                               EIG_R_mat_2, EIG_L_mat_2, &
                               RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, &
                               VECR_2, VECL_2, comm_grp=newcomm)
                ALLOCATE(ROR_2_COPY(NDIMR), ROP_2_COPY(NDIMR), OZ_2_COPY(NDIMR))
                ROR_2_COPY = ROR_2; ROP_2_COPY = ROP_2; OZ_2_COPY = OZ_2
                
                ! Calculate nonlinear interactions
                ! {M2,K2}: 0 X 1 -> 2
                CALL NONLIN_MK(CURRENT_TRIAD%m2, CURRENT_TRIAD%k2, &
                              RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0, &
                              RUR_1, RUP_1, UZ_1, ROR_1_COPY, ROP_1_COPY, OZ_1_COPY, &
                              VECR_01, comm_grp=newcomm)
                
                ! {M1,K1}: 0* X 2 -> 1
                CALL NONLIN_MK(CURRENT_TRIAD%m1, CURRENT_TRIAD%k1, &
                              CONJG(RUR_0), CONJG(RUP_0), CONJG(UZ_0), &
                              CONJG(ROR_0), CONJG(ROP_0), CONJG(OZ_0), &
                              RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, &
                              VECR_02, comm_grp=newcomm)
                
                ! {M0,K0}: 1* X 2 -> 0
                CALL NONLIN_MK(CURRENT_TRIAD%m0, CURRENT_TRIAD%k0, &
                              CONJG(RUR_1), CONJG(RUP_1), CONJG(UZ_1), &
                              CONJG(ROR_1), CONJG(ROP_1), CONJG(OZ_1), &
                              RUR_2, RUP_2, UZ_2, ROR_2_COPY, ROP_2_COPY, OZ_2_COPY, &
                              VECR_12, comm_grp=newcomm)
                
                ! Calculate interaction coefficients
                J0 = DOT_PRODUCT(VECL_0, VECR_12)
                J1 = DOT_PRODUCT(VECL_1, VECR_02)
                J2 = DOT_PRODUCT(VECL_2, VECR_01)
                
                ! Check if Re(J0*J2) > 0 or Re(J1*J2) > 0
                ! For purely imaginary J, this means opposite signs in imaginary parts
                INTERACTION_FLAG = (REAL(J0*J2) > 0.0D0) .OR. (REAL(J1*J2) > 0.0D0)
                
                ! Write results (using mode indices instead of branch indices)
                WRITE(FID, 100) CURRENT_TRIAD%m0, CURRENT_TRIAD%k0, EIG_IDX_0, &
                               CURRENT_TRIAD%m1, CURRENT_TRIAD%k1, EIG_IDX_1, &
                               CURRENT_TRIAD%m2, CURRENT_TRIAD%k2, EIG_IDX_2, &
                               M_eig_0(EIG_IDX_0), M_eig_1(EIG_IDX_1), M_eig_2(EIG_IDX_2), &
                               J0, J1, J2, INTERACTION_FLAG
100             FORMAT(3(I4,1X,F12.6,1X,I4,1X),3(SP,E20.12,1X,E20.12,'i',1X),3(SP,E20.12,1X,E20.12,'i',1X),L1)
                
                ! Deallocate
                DEALLOCATE(RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0)
                DEALLOCATE(RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1)
                DEALLOCATE(RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2)
                DEALLOCATE(ROR_0_COPY, ROP_0_COPY, OZ_0_COPY)
                DEALLOCATE(ROR_1_COPY, ROP_1_COPY, OZ_1_COPY)
                DEALLOCATE(ROR_2_COPY, ROP_2_COPY, OZ_2_COPY)
                DEALLOCATE(VECR_0, VECL_0, VECR_1, VECL_1, VECR_2, VECL_2)
                DEALLOCATE(VECR_01, VECR_02, VECR_12)
            ENDIF
            
            ! Clean up EVP for m1, m2
            IF (ALLOCATED(M_eig_1)) DEALLOCATE(M_eig_1)
            IF (ALLOCATED(M_eig_2)) DEALLOCATE(M_eig_2)
            IF (ALLOCATED(EIG_R_mat_1)) DEALLOCATE(EIG_R_mat_1)
            IF (ALLOCATED(EIG_L_mat_1)) DEALLOCATE(EIG_L_mat_1)
            IF (ALLOCATED(EIG_R_mat_2)) DEALLOCATE(EIG_R_mat_2)
            IF (ALLOCATED(EIG_L_mat_2)) DEALLOCATE(EIG_L_mat_2)
            IF (ALLOCATED(M_mat_1)) DEALLOCATE(M_mat_1)
            IF (ALLOCATED(M_mat_2)) DEALLOCATE(M_mat_2)
            
        END DO  ! II (loop over all triads)
        
        CLOSE(FID)
        
        ! Clean up EVP for m0
        IF (ALLOCATED(M_eig_0)) DEALLOCATE(M_eig_0)
        IF (ALLOCATED(EIG_R_mat_0)) DEALLOCATE(EIG_R_mat_0)
        IF (ALLOCATED(EIG_L_mat_0)) DEALLOCATE(EIG_L_mat_0)
        IF (ALLOCATED(M_mat_0)) DEALLOCATE(M_mat_0)
        
    END DO  ! WHILE
    
    WRITE(*,*) 'Rank', MPI_GLB_RANK, '- FINISHED'
    
ENDIF

!=======================================================================
! SYNCHRONIZE AND ASSEMBLE RESULTS
!=======================================================================
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE(*,*) ''
    WRITE(*,*) 'Assembling results from all processes...'
    
    FINAL_OUTPUT = './data/resonance_interactions_all.dat'
    OPEN(UNIT=FID, FILE=TRIM(FINAL_OUTPUT), STATUS='unknown', ACTION='WRITE', IOSTAT=STAT)
    
    ! Write header
    WRITE(FID, '(A)') '# Resonance Triad Interaction Coefficients'
    WRITE(FID, '(A)') '# Format: m0 k0 idx0 m1 k1 idx1 m2 k2 idx2 omega0 omega1 omega2 J0 J1 J2 Flag'
    WRITE(FID, '(A)') '# idx = eigenvalue mode index (not branch index)'
    WRITE(FID, '(A)') '# Flag: T if Re(J0*J2)>0 or Re(J1*J2)>0, F otherwise'
    WRITE(FID, '(A)') '#'
    
    ! Concatenate all worker files
    DO II = 1, MPI_GLB_PROCS - 1
        WRITE(OUTPUT_FILE, '(A,I0.4,A)') './data/interaction_rank_', II, '.dat'
        INQUIRE(FILE=TRIM(OUTPUT_FILE), EXIST=FILE_EXISTS)
        
        IF (FILE_EXISTS) THEN
            OPEN(UNIT=98, FILE=TRIM(OUTPUT_FILE), STATUS='OLD', ACTION='READ', IOSTAT=STAT)
            IF (STAT == 0) THEN
                DO
                    READ(98, '(A)', IOSTAT=STAT) LINE
                    IF (STAT /= 0) EXIT
                    WRITE(FID, '(A)') TRIM(LINE)
                END DO
                CLOSE(98)
            ENDIF
        ENDIF
    END DO
    
    CLOSE(FID)
    WRITE(*,*) 'Results assembled in:', TRIM(FINAL_OUTPUT)
ENDIF

!=======================================================================
! FINALIZE
!=======================================================================
IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE(*,*) ''
    WRITE(*,*) '======================================================='
    WRITE(*,*) 'PROGRAM FINISHED'
    WRITE(*,*) '======================================================='
    CALL PRINT_REAL_TIME()
ENDIF

CALL MPI_FINALIZE(IERR)

!=======================================================================
CONTAINS
!=======================================================================

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

!-----------------------------------------------------------------------
SUBROUTINE EXTRACT_UNIQUE_M0K0(TRIADS, NTRIADS, UNIQUE_M0K0, UNIQUE_K0_VALS, NUNIQUE, &
                                 USE_M0_FILTER, M0_FILTER, USE_K0_FILTER, K0_FILTER, K0_IDX_FILTER)
    TYPE(TRIAD_DATA), DIMENSION(:), INTENT(IN) :: TRIADS
    INTEGER, INTENT(IN) :: NTRIADS
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: UNIQUE_M0K0
    REAL(P8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: UNIQUE_K0_VALS
    INTEGER, INTENT(OUT) :: NUNIQUE
    LOGICAL, INTENT(IN) :: USE_M0_FILTER, USE_K0_FILTER
    INTEGER, INTENT(IN) :: M0_FILTER
    REAL(P8), INTENT(IN) :: K0_FILTER
    INTEGER, INTENT(OUT) :: K0_IDX_FILTER
    
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: TEMP_M0K0
    REAL(P8), DIMENSION(:), ALLOCATABLE :: TEMP_K0_VALS
    INTEGER :: I, J, K0_IDX, NUM_K0
    LOGICAL :: FOUND
    REAL(P8), PARAMETER :: K_TOL = 1.0D-10
    
    ALLOCATE(TEMP_M0K0(NTRIADS, 2))
    ALLOCATE(TEMP_K0_VALS(NTRIADS))
    NUNIQUE = 0
    NUM_K0 = 0
    K0_IDX_FILTER = 0
    
    DO I = 1, NTRIADS
        ! Apply filters if requested
        IF (USE_M0_FILTER .AND. TRIADS(I)%m0 /= M0_FILTER) CYCLE
        IF (USE_K0_FILTER .AND. ABS(TRIADS(I)%k0 - K0_FILTER) > K_TOL) CYCLE
        
        ! Check if (m0, k0) already in list
        FOUND = .FALSE.
        DO J = 1, NUNIQUE
            IF (TEMP_M0K0(J,1) == TRIADS(I)%m0 .AND. &
                ABS(TEMP_K0_VALS(TEMP_M0K0(J,2)) - TRIADS(I)%k0) < K_TOL) THEN
                FOUND = .TRUE.
                EXIT
            ENDIF
        END DO
        
        IF (.NOT. FOUND) THEN
            ! Find or add k0 value to unique list
            K0_IDX = 0
            DO J = 1, NUM_K0
                IF (ABS(TEMP_K0_VALS(J) - TRIADS(I)%k0) < K_TOL) THEN
                    K0_IDX = J
                    EXIT
                ENDIF
            END DO
            
            IF (K0_IDX == 0) THEN
                ! New k0 value
                NUM_K0 = NUM_K0 + 1
                K0_IDX = NUM_K0
                TEMP_K0_VALS(K0_IDX) = TRIADS(I)%k0
            ENDIF
            
            ! Add new (m0, k0_idx) pair
            NUNIQUE = NUNIQUE + 1
            TEMP_M0K0(NUNIQUE, 1) = TRIADS(I)%m0
            TEMP_M0K0(NUNIQUE, 2) = K0_IDX
            
            IF (USE_K0_FILTER .AND. ABS(TRIADS(I)%k0 - K0_FILTER) < K_TOL) THEN
                K0_IDX_FILTER = K0_IDX
            ENDIF
        ENDIF
    END DO
    
    ! Copy to output arrays
    ALLOCATE(UNIQUE_M0K0(NUNIQUE, 2))
    ALLOCATE(UNIQUE_K0_VALS(NUM_K0))
    UNIQUE_M0K0 = TEMP_M0K0(1:NUNIQUE, :)
    UNIQUE_K0_VALS = TEMP_K0_VALS(1:NUM_K0)
    
    DEALLOCATE(TEMP_M0K0, TEMP_K0_VALS)
END SUBROUTINE EXTRACT_UNIQUE_M0K0

!-----------------------------------------------------------------------
SUBROUTINE FIND_MODE(M_eig, omega_target, omega_tol_max, real_tol, &
                      eig_idx, omega_err, found)
    COMPLEX(P8), DIMENSION(:), INTENT(IN) :: M_eig
    REAL(P8), INTENT(IN) :: omega_target, omega_tol_max, real_tol
    INTEGER, INTENT(OUT) :: eig_idx
    REAL(P8), INTENT(OUT) :: omega_err
    LOGICAL, INTENT(OUT) :: found
    
    INTEGER :: I, best_idx
    REAL(P8) :: min_err, curr_err
    LOGICAL :: has_acceptable
    
    found = .FALSE.
    eig_idx = 0
    omega_err = 1.0D10
    min_err = 1.0D10
    best_idx = 0
    has_acceptable = .FALSE.
    
    ! First pass: look for modes with |Re(omega)| < real_tol
    DO I = 1, SIZE(M_eig)
        IF (ABS(REAL(M_eig(I))) < real_tol) THEN
            curr_err = ABS(AIMAG(M_eig(I)) - omega_target)
            IF (curr_err < min_err) THEN
                min_err = curr_err
                best_idx = I
                has_acceptable = .TRUE.
            ENDIF
        ENDIF
    END DO
    
    ! Second pass: if no acceptable mode found, use closest regardless
    IF (.NOT. has_acceptable) THEN
        DO I = 1, SIZE(M_eig)
            curr_err = ABS(AIMAG(M_eig(I)) - omega_target)
            IF (curr_err < min_err) THEN
                min_err = curr_err
                best_idx = I
            ENDIF
        END DO
    ENDIF
    
    IF (best_idx > 0) THEN
        eig_idx = best_idx
        omega_err = min_err
        ! Accept if within tolerance OR if Re(omega) small enough
        IF (omega_err < omega_tol_max .OR. ABS(REAL(M_eig(best_idx))) < real_tol) THEN
            found = .TRUE.
        ENDIF
    ENDIF
    
END SUBROUTINE FIND_MODE

!-----------------------------------------------------------------------
SUBROUTINE PRINT_REAL_TIME_TO_UNIT(UNIT_NUM)
    INTEGER, INTENT(IN) :: UNIT_NUM
    INTEGER :: TIME_ARRAY(8)
    CALL DATE_AND_TIME(VALUES=TIME_ARRAY)
    WRITE(UNIT_NUM, '(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
        TIME_ARRAY(1), '/', TIME_ARRAY(2), '/', TIME_ARRAY(3), ' ', &
        TIME_ARRAY(5), ':', TIME_ARRAY(6), ':', TIME_ARRAY(7)
END SUBROUTINE PRINT_REAL_TIME_TO_UNIT

END PROGRAM EVP_RESONANCE_INTERACTION
!=======================================================================
