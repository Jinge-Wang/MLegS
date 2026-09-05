PROGRAM EVP_PRINT
!=======================================================================
! [USAGE]:
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
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
REAL(P8)    :: K_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR

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

! M AND K OF THE EVP & H OF THE Q-VORTEX
! - Default values
M_BAR = 1
K_BAR = 0.d0
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
                
            CASE DEFAULT
                ! Unknown parameter
                IF (MPI_GLB_RANK .EQ. 0) THEN
                    WRITE(*,*) 'Unknown parameter:', TRIM(param_name)
                    WRITE(*,*) 'Supported parameters: m, k, h, nr'
                END IF
        END SELECT
    END DO
    
    ! Display the parameters on rank 0
    IF (MPI_GLB_RANK .EQ. 0) THEN
        WRITE(*,*) 'Using parameters:'
        WRITE(*,*) '  m =', M_BAR
        WRITE(*,*) '  k =', K_BAR
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

! OBTAIN EVP MATRIX, EIGENVALUE, AND EIGENVECTORS
! CALL EIG_MATRIX_SERIAL(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, print_switch=.true.)
CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, serial_switch=.false.) 
NR_MK = NRCHOPS(2)

! OUTPUT
IF (MPI_GLB_RANK .NE. 0) THEN
    GOTO 129
ELSE
    ! SAVE BASEFLOW
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

    ! SAVE EIGVALS
    WRITE(FILENAME,'(A,I0,A,F0.3,A,I0,A)') './data/eig_m_', M_BAR, '_k_', K_BAR, '_nr_', NRCHOP, '.output'
    open(FID, FILE=TRIM(FILENAME), STATUS='unknown', ACTION='WRITE', IOSTAT=IERR)
    DO II = 1, SIZE(M_eig)
        WRITE(FID,*) 'II = ', II, ':', M_eig(II), '-', EIGRES(EIG_R_mat(:2*NR_MK,II), M_BAR)
    ENDDO
    close(FID)

    ! Create subfolder for velocity files
    WRITE(dir_path,'(A,I0,A,F0.3)') './data/vel_m_', M_BAR, '_k_', K_BAR
    CALL SYSTEM('mkdir -p ' // TRIM(dir_path))
    WRITE(*,*) 'Created directory: ', TRIM(dir_path)
    
    ALLOCATE(EIG_R_BAR(2*NR_MK))
    DO II = 1,SIZE(M_eig)
        EIG_R_BAR = EIG_R_mat(:2*NR_MK,II)
        CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
        WRITE(FILENAME,'(A,A,I0,A)') TRIM(dir_path), '/ind_', II, '.output'
        CALL SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, TRIM(FILENAME))
        DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
    ENDDO
    DEALLOCATE(EIG_R_BAR,M_eig)
ENDIF

129     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)


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

END PROGRAM EVP_PRINT
!=======================================================================
