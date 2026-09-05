PROGRAM EVP_RESONANCE_COMBINE
!=======================================================================
! [USAGE]:
! Combine individual interaction data files from parallel processes
! into a single unified output file
!
! [DESCRIPTION]:
! This program collects all interaction_rank_*.dat files from the
! data/ directory and combines them into resonance_interactions_all.dat
! Useful for recovering data from interrupted runs or manually
! combining results from multiple runs.
!
! The original individual files are preserved.
!
! [INPUT FILES]:
! ./data/interaction_rank_0001.dat
! ./data/interaction_rank_0002.dat
! ... (all rank files found in data/ directory)
!
! [OUTPUT FILE]:
! ./data/resonance_interactions_all.dat
!
! [COMMAND LINE ARGUMENTS]:
! output=<filename>  : Specify custom output filename (default: resonance_interactions_all.dat)
!
! [BASED ON]:
! evp_resonance_interaction.f90 (assembly section)
!=======================================================================
IMPLICIT NONE

! FILE I/O
CHARACTER(LEN=256) :: INPUT_FILE, OUTPUT_FILE, DATA_DIR
CHARACTER(LEN=512) :: LINE, SYSTEM_CMD
INTEGER :: FID_IN, FID_OUT, STAT, RANK_NUM
INTEGER :: NUM_FILES, TOTAL_LINES
LOGICAL :: FILE_EXISTS
INTEGER :: II, JJ

! COMMAND LINE ARGUMENTS
INTEGER :: num_args, eq_pos, arg_len
CHARACTER(LEN=256) :: arg_buffer
CHARACTER(LEN=32) :: param_name, param_value

! STATISTICS
INTEGER :: FILES_FOUND, FILES_PROCESSED, LINES_WRITTEN

!=======================================================================
! PARSE COMMAND LINE ARGUMENTS
!=======================================================================
OUTPUT_FILE = './data/resonance_interactions_all.dat'

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
            CASE('output')
                OUTPUT_FILE = TRIM(param_value)
                WRITE(*,*) 'Output file:', TRIM(OUTPUT_FILE)
        END SELECT
    END DO
ENDIF

!=======================================================================
! INITIALIZE
!=======================================================================
WRITE(*,*) '======================================================='
WRITE(*,*) 'RESONANCE INTERACTION DATA COMBINER'
WRITE(*,*) '======================================================='
CALL PRINT_REAL_TIME()
WRITE(*,*) ''

DATA_DIR = './data/'
FILES_FOUND = 0
FILES_PROCESSED = 0
LINES_WRITTEN = 0

! Assign explicit file unit numbers
FID_OUT = 20
FID_IN = 21

!=======================================================================
! SCAN FOR INPUT FILES
!=======================================================================
WRITE(*,*) 'Scanning for interaction_rank_*.dat files in', TRIM(DATA_DIR)

! Count available files by checking rank numbers from 0001 to 9999
DO RANK_NUM = 1, 9999
    WRITE(INPUT_FILE, '(A,I0.4,A)') TRIM(DATA_DIR)//'interaction_rank_', RANK_NUM, '.dat'
    INQUIRE(FILE=TRIM(INPUT_FILE), EXIST=FILE_EXISTS)
    IF (FILE_EXISTS) THEN
        FILES_FOUND = FILES_FOUND + 1
    ENDIF
END DO

WRITE(*,*) 'Found', FILES_FOUND, 'data files'

IF (FILES_FOUND == 0) THEN
    WRITE(*,*) 'ERROR: No interaction_rank_*.dat files found in', TRIM(DATA_DIR)
    WRITE(*,*) 'Please check that the directory exists and contains data files.'
    STOP
ENDIF

!=======================================================================
! CREATE OUTPUT FILE WITH HEADER
!=======================================================================
WRITE(*,*) ''
WRITE(*,*) 'Creating combined output file:', TRIM(OUTPUT_FILE)

OPEN(UNIT=FID_OUT, FILE=TRIM(OUTPUT_FILE), STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=STAT)
IF (STAT /= 0) THEN
    WRITE(*,*) 'ERROR: Cannot create output file:', TRIM(OUTPUT_FILE)
    STOP
ENDIF

! Write header
WRITE(FID_OUT, '(A)') '# Resonance Triad Interaction Coefficients'
WRITE(FID_OUT, '(A)') '# Combined from individual process files'
WRITE(FID_OUT, '(A)') '# Format: m0 k0 idx0 m1 k1 idx1 m2 k2 idx2 omega0 omega1 omega2 J0 J1 J2 Flag'
WRITE(FID_OUT, '(A)') '# idx = eigenvalue mode index (not branch index)'
WRITE(FID_OUT, '(A)') '# Flag: T if Re(J0*J2)>0 or Re(J1*J2)>0, F otherwise'
WRITE(FID_OUT, '(A)') '#'

!=======================================================================
! COMBINE ALL FILES
!=======================================================================
WRITE(*,*) ''
WRITE(*,*) 'Combining files...'

DO RANK_NUM = 1, 9999
    WRITE(INPUT_FILE, '(A,I0.4,A)') TRIM(DATA_DIR)//'interaction_rank_', RANK_NUM, '.dat'
    INQUIRE(FILE=TRIM(INPUT_FILE), EXIST=FILE_EXISTS)
    
    IF (FILE_EXISTS) THEN
        FILES_PROCESSED = FILES_PROCESSED + 1
        
        ! Open input file
        OPEN(UNIT=FID_IN, FILE=TRIM(INPUT_FILE), STATUS='OLD', ACTION='READ', IOSTAT=STAT)
        
        IF (STAT == 0) THEN
            ! Count lines in this file
            NUM_FILES = 0
            DO
                READ(FID_IN, '(A)', IOSTAT=STAT) LINE
                IF (STAT /= 0) EXIT
                IF (LEN_TRIM(LINE) > 0) THEN
                    NUM_FILES = NUM_FILES + 1
                ENDIF
            END DO
            REWIND(FID_IN)
            
            ! Copy all lines to output
            DO
                READ(FID_IN, '(A)', IOSTAT=STAT) LINE
                IF (STAT /= 0) EXIT
                WRITE(FID_OUT, '(A)') TRIM(LINE)
                LINES_WRITTEN = LINES_WRITTEN + 1
            END DO
            
            CLOSE(FID_IN)
            
            WRITE(*,'(A,I0.4,A,I8,A)') '  Processed rank ', RANK_NUM, ': ', NUM_FILES, ' lines'
        ELSE
            WRITE(*,*) '  WARNING: Could not open', TRIM(INPUT_FILE)
        ENDIF
    ENDIF
END DO

CLOSE(FID_OUT)

!=======================================================================
! SUMMARY
!=======================================================================
WRITE(*,*) ''
WRITE(*,*) '======================================================='
WRITE(*,*) 'COMBINATION COMPLETE'
WRITE(*,*) '======================================================='
WRITE(*,*) 'Files processed:', FILES_PROCESSED
WRITE(*,*) 'Total lines written:', LINES_WRITTEN
WRITE(*,*) 'Output file:', TRIM(OUTPUT_FILE)
WRITE(*,*) ''
CALL PRINT_REAL_TIME()
WRITE(*,*) '======================================================='

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
SUBROUTINE PRINT_REAL_TIME()
    INTEGER :: TIME_ARRAY(8)
    CALL DATE_AND_TIME(VALUES=TIME_ARRAY)
    WRITE(*, '(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
        TIME_ARRAY(1), '/', TIME_ARRAY(2), '/', TIME_ARRAY(3), ' ', &
        TIME_ARRAY(5), ':', TIME_ARRAY(6), ':', TIME_ARRAY(7)
END SUBROUTINE PRINT_REAL_TIME

END PROGRAM EVP_RESONANCE_COMBINE
!=======================================================================
