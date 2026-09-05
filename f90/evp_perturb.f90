PROGRAM EVP_PERTURB
!=======================================================================
! [USAGE]:
! COMPUTE EIGENVALUE PERTURBATIONS WITH RESPECT TO AXIAL WAVENUMBER K
! FOR THE LINEARIZED N-S EQUATIONS IN POLOIDAL-TOROIDAL DECOMPOSITION
!
! CALCULATES: (L^H_i * M(k+dk) * R_i) / (L^H_i * R_i) - sigma_i
! WHERE:
!   - M(k) IS THE EIGENMATRIX AT WAVENUMBER k
!   - L_i, R_i ARE LEFT AND RIGHT EIGENVECTORS FOR EIGENVALUE sigma_i
!   - dk IS THE PERTURBATION IN k (default: 0.01% of k)
!
! OUTPUT: ./data/per_m_{m}_k_{k}_dk_{dk}_nr_{nr}.output
!   FORMAT: II = {index} : {eigenvalue} - {residual} - {d(sigma)/dk}
!
! PARAMETERS:
!   m=<int>       : Azimuthal wavenumber
!   k=<float>     : Axial wavenumber
!   dk=<float>    : Perturbation in k (optional, default 0.01% of k)
!   nr=<int>      : Number of radial modes
!   base=<bool>   : Save baseflow (optional, default false)
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
INTEGER     :: II, JJ, KK
INTEGER     :: M_BAR, NR_MK
REAL(P8)    :: K_BAR, DELTA_K, K_PERT
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_pert
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: VEC_L, VEC_R, PERTURBATION
COMPLEX(P8) :: NORM_LR, PERTURB_VAL

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
LOGICAL :: file_exists, SAVE_BASE

! MPI INITILIZATION
CALL SETUP_ENVIRONMENT('NOECHO')
! PARALLEL COMPUTATION
! - SPLIT THE FORMATION OF H INTO DIFFERENT PROCS
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS

! M AND K OF THE EVP & H OF THE Q-VORTEX
! - Default values
M_BAR = 1
K_BAR = 0.d0
DELTA_K = 0.d0  ! Will be set to 1% of K_BAR if not specified
SAVE_BASE = .FALSE.  ! Only save baseflow if explicitly requested
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
            CASE('m')
                ! For m, read as real and convert to nearest integer
                READ(param_value, *, IOSTAT=stat) temp_real
                IF (stat == 0) THEN
                    M_BAR = NINT(temp_real)  ! Convert to nearest integer
                ELSE
                    IF (MPI_GLB_RANK .EQ. 0) THEN
                        WRITE(*,*) 'Error reading m parameter:', TRIM(param_value)
                    END IF
                END IF
                
            CASE('k')
                ! Read k as a real number
                READ(param_value, *, IOSTAT=stat) K_BAR
                IF (stat /= 0 .AND. MPI_GLB_RANK .EQ. 0) THEN
                    WRITE(*,*) 'Error reading k parameter:', TRIM(param_value)
                END IF
                
            CASE('nr')
                ! Read nr parameter as integer
                READ(param_value, *, IOSTAT=stat) temp_real
                IF (stat == 0) THEN
                    NRCHOP = NINT(temp_real)  ! Convert to nearest integer
                ELSE
                    IF (MPI_GLB_RANK .EQ. 0) THEN
                        WRITE(*,*) 'Error reading nr parameter:', TRIM(param_value)
                    END IF
                END IF
                
            CASE('dk', 'deltak')
                ! Read delta_k parameter as real number
                READ(param_value, *, IOSTAT=stat) DELTA_K
                IF (stat /= 0 .AND. MPI_GLB_RANK .EQ. 0) THEN
                    WRITE(*,*) 'Error reading dk parameter:', TRIM(param_value)
                END IF
                
            CASE('base')
                ! Read base parameter as logical
                CALL lowercase(param_value)
                IF (TRIM(param_value) == 'true' .OR. TRIM(param_value) == 't' .OR. TRIM(param_value) == '1') THEN
                    SAVE_BASE = .TRUE.
                ELSE IF (TRIM(param_value) == 'false' .OR. TRIM(param_value) == 'f' .OR. TRIM(param_value) == '0') THEN
                    SAVE_BASE = .FALSE.
                ELSE
                    IF (MPI_GLB_RANK .EQ. 0) THEN
                        WRITE(*,*) 'Error reading base parameter:', TRIM(param_value)
                        WRITE(*,*) 'Expected: true/false, t/f, or 1/0'
                    END IF
                END IF
                
            CASE DEFAULT
                ! Unknown parameter
                IF (MPI_GLB_RANK .EQ. 0) THEN
                    WRITE(*,*) 'Unknown parameter:', TRIM(param_name)
                    WRITE(*,*) 'Supported parameters: m, k, h, nr, dk, base'
                END IF
        END SELECT
    END DO
    
    ! Set default delta_k if not specified (1% of k)
    IF (DELTA_K .EQ. 0.d0) THEN
        DELTA_K = 0.0001d0 * K_BAR
    END IF
    
    ! Calculate perturbed k
    K_PERT = K_BAR + DELTA_K
    
    ! Display the parameters on rank 0
    IF (MPI_GLB_RANK .EQ. 0) THEN
        WRITE(*,*) 'Using parameters:'
        WRITE(*,*) '  m =', M_BAR
        WRITE(*,*) '  k =', K_BAR
        WRITE(*,*) '  dk =', DELTA_K
        WRITE(*,*) '  k_pert =', K_PERT
        WRITE(*,*) '  save_base =', SAVE_BASE
        WRITE(*,*) '  bv =', BSNSQ%BV0
        WRITE(*,*) '  w =', BSNSQ%OMEGA
        WRITE(*,*) '  nr =', NRCHOP
        WRITE(*,*) '  Q =', QPAIR%Q(1)
        WRITE(*,*) '  H =', QPAIR%H(1)
    END IF
ENDIF  

! I/O
WRITE(K_VAL,'(F06.3)') K_BAR
WRITE(H_VAL,'(F06.3)') QPAIR%H(1)
! WRITE(*,*) MPI_GLB_RANK,': M = ',M_BAR,'; K = ',K_BAR

! OBTAIN EVP MATRIX, EIGENVALUE, AND EIGENVECTORS FOR UNPERTURBED CASE
IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Computing eigenmatrix for k =', K_BAR
CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, serial_switch=.false.) 
NR_MK = NRCHOPS(2)

! OBTAIN EVP MATRIX FOR PERTURBED CASE (k + delta_k)
IF (MPI_GLB_RANK .EQ. 0) WRITE(*,*) 'Computing eigenmatrix for k_pert =', K_PERT
CALL EIG_MATRIX(M_BAR, K_PERT, M_mat_pert, comm_grp=newcomm, serial_switch=.false.)

! CALCULATE PERTURBATIONS: L^H_i * M(k+dk) * R_i / (L^H_i * R_i) - sigma_i
IF (MPI_GLB_RANK .EQ. 0) THEN
    ALLOCATE(PERTURBATION(SIZE(M_eig)))
    ALLOCATE(VEC_L(2*NR_MK))
    ALLOCATE(VEC_R(2*NR_MK))
    
    WRITE(*,*) 'Computing perturbations...'
    DO II = 1, SIZE(M_eig)
        ! Extract left and right eigenvectors for this eigenvalue
        VEC_L(:) = EIG_L_mat(:2*NR_MK, II)
        VEC_R(:) = EIG_R_mat(:2*NR_MK, II)
        
        ! Compute L^H_i * R_i (normalization factor)
        NORM_LR = DOT_PRODUCT(VEC_L, VEC_R)
        
        ! Compute L^H_i * M(k+dk) * R_i
        PERTURB_VAL = DOT_PRODUCT(VEC_L, MATMUL(M_mat_pert(:2*NR_MK, :2*NR_MK), VEC_R))
        
        ! Calculate perturbation: (L^H * M_pert * R) / (L^H * R) - sigma
        PERTURBATION(II) = (PERTURB_VAL / NORM_LR - M_eig(II)) / DELTA_K
    END DO
    
    DEALLOCATE(VEC_L, VEC_R)
ENDIF

! OUTPUT
IF (MPI_GLB_RANK .NE. 0) THEN
    GOTO 129
ELSE
    ! SAVE BASEFLOW (only if requested)
    IF (SAVE_BASE) THEN
        WRITE(FILENAME,'(A)') './data/base.output'
        INQUIRE(FILE=TRIM(FILENAME), EXIST=file_exists)
        IF (.NOT. file_exists) THEN
            open(FID, FILE=TRIM(FILENAME), STATUS='unknown', ACTION='WRITE', IOSTAT=IERR)
            DO II = 1, SIZE(TFM%R)
                WRITE(FID,*) TFM%R(II), RUR0(II)/TFM%R(II), RUP0(II)/TFM%R(II), UZ0(II)
            ENDDO
            close(FID)
            WRITE(*,*) 'Saved baseflow to: ', TRIM(FILENAME)
        ENDIF
    ENDIF

    ! SAVE PERTURBATION RESULTS (with eigenvalues, residuals, and perturbations)
    WRITE(FILENAME,'(A,I0,A,F0.3,A,F0.6,A,I0,A)') './data/per_m_', M_BAR, '_k_', K_BAR, '_dk_', DELTA_K, '_nr_', NRCHOP, '.output'
    open(FID, FILE=TRIM(FILENAME), STATUS='unknown', ACTION='WRITE', IOSTAT=IERR)
    DO II = 1, SIZE(M_eig)
        WRITE(FID,*) 'II = ', II, ':', M_eig(II), '-', EIGRES(EIG_R_mat(:2*NR_MK,II), M_BAR), '-', PERTURBATION(II)
    ENDDO
    close(FID)
    WRITE(*,*) 'Saved perturbation results to: ', TRIM(FILENAME)
    
    DEALLOCATE(PERTURBATION, M_eig)
ENDIF

129     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)


IF (ALLOCATED(M_mat)) DEALLOCATE(M_mat)
IF (ALLOCATED(M_mat_pert)) DEALLOCATE(M_mat_pert)
IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)

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

END PROGRAM EVP_PERTURB
!=======================================================================
