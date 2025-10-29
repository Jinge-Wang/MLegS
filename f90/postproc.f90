!=======================================================================
!
!     WRITTEN BY JINGE WANG
!     DEPT. OF MECHANICAL ENGINEERING
!     UNIV. OF CALIFORNIA AT BERKELEY
!     EMAIL: JINGE@BERKELEY.EDU
!     
!     UC BERKELEY CFD LAB
!     HTTPS://CFD.ME.BERKELEY.EDU/
!
!     NONCOMMERCIAL USE WITH COPYRIGHTED (C) MARK
!     UNDER DEVELOPMENT FOR RESEARCH PURPOSE 
!
!=======================================================================
PROGRAM POSTPROC
!=======================================================================
! [USAGE]:
! This program post-processes simulation data snapshots. It reads the 
! spectral coefficients for the poloidal (CHI) and toroidal (PSI) 
! velocity potentials, and buoyancy (B), then computes and saves 
! physical fields like velocity and vorticity. 
!
! Postproc is defined by the 'POSTPROCESS' section in 'read.input'.
!
! [INPUT PARAMETERS (from read.input)]:
!
! POSTPROCESS%JOB: A 3-digit integer XYZ that specifies the job.
!   - X (Field Type): 
!     - 0: Process both Velocity and Vorticity fields. 
!     - 1: Process Velocity fields only. 
!     - 2: Process Vorticity fields only. 
!   - Y (Component): 
!     - 0: Process all components (R, Theta, Z) and Buoyancy. 
!     - 1: Process R-component only. 
!     - 2: Process Theta-component only. 
!     - 3: Process Z-component only. 
!     - 4: Process Buoyancy field only. 
!   - Z (Slice Orientation): 
!     - 1: Save R-Theta plane(s). 
!     - 2: Save R-Z plane(s). 
!
! POSTPROCESS%SLICEINT: Specifies the index of the plane to be saved.
!   - For R-Theta slices (Z=1), this is the Z-index. 
!   - For R-Z slices (Z=2), this is the Theta-index. 
!   - Set to 999 to save the full 3D field. 
!
! POSTPROCESS%START, POSTPROCESS%FINISH: The range of snapshot indices 
! to process.
!   - Set START = -1 to save the postprocessed data only.
!=======================================================================
USE OMP_LIB
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
USE MOD_DIAGNOSTICS
USE MOD_INIT

IMPLICIT NONE
TYPE(SCALAR):: PSI,CHI,B
INTEGER :: I, JOB_X, JOB_Y, JOB_Z, JOB_NUM, JOB_IND
LOGICAL:: SAVE_HIGHFREQ
! I/O
real(P8):: time_start,time_end

CALL SETUP_ENVIRONMENT('noecho')
CALL SETUP_GRID(FILES%SAVEDIR)
IF (MPI_RANK.EQ.0) time_start = mpi_wtime()

IF (.NOT. ALLOCATED(POSTPROCESS%START)) THEN
    IF (MPI_RANK.EQ.0) WRITE(*,*) 'POSTPROCESS: No postprocess jobs defined.'
    CALL MPI_FINALIZE(IERR)
    STOP
ENDIF

JOB_NUM = SIZE(POSTPROCESS%START)
DO JOB_IND = 1, JOB_NUM

    SAVE_HIGHFREQ = POSTPROCESS%JOB(JOB_IND) .LT. 0
    JOB_X = ABS(POSTPROCESS%JOB(JOB_IND)) / 100
    JOB_Y = MOD(ABS(POSTPROCESS%JOB(JOB_IND)), 100) / 10
    JOB_Z = MOD(ABS(POSTPROCESS%JOB(JOB_IND)), 10)

    IF (MPI_RANK.EQ.0) THEN
        WRITE(*,*) 'Postprocessing Job:', POSTPROCESS%JOB(JOB_IND)
        WRITE(*,*) 'Field (X):', JOB_X, ' Component (Y):', JOB_Y, ' Slice (Z):', JOB_Z
    ENDIF

    DO I=POSTPROCESS%START(JOB_IND),POSTPROCESS%FINISH(JOB_IND)

        CALL MPI_BARRIER(MPI_COMM_IVP,IERR)

        ! Check if required PSI file exists before attempting to load
        IF (POSTPROCESS%START(JOB_IND) .EQ. -1) THEN
            IF (.NOT. FILE_EXISTS(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSII)) THEN
            IF (MPI_RANK.EQ.0) THEN
                WRITE(*,*) 'POSTPROCESS: Initial condition PSI file missing, skipping...'
                WRITE(*,*) '  PSI:', TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSII
            ENDIF
            CYCLE
            ENDIF
        ELSE
            IF (.NOT. FILE_EXISTS(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI(I))) THEN
            IF (MPI_RANK.EQ.0) THEN
                WRITE(*,*) 'POSTPROCESS: Snapshot', I, 'PSI file missing, skipping...'
                WRITE(*,*) '  PSI:', TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI(I)
            ENDIF
            CYCLE
            ENDIF
        ENDIF

        ! Load data for the current snapshot
        CALL ALLOCATE(PSI); CALL ALLOCATE(CHI); CALL ALLOCATE(B)
        IF (POSTPROCESS%START(JOB_IND) .EQ. -1) THEN
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSII,PSI)
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHII,CHI)
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%BI,B)
        ELSE
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI(I),PSI)
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI(I),CHI)
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%B(I),B)
        ENDIF

        ! Optional high-pass filtering
        ! This step keeps only the high-frequency components for debug visualization.
        IF (SAVE_HIGHFREQ) THEN
            ! CALL HIGH_PASS_FILTER(PSI)
            PSI%E = (0.0D0, 0.0D0); PSI%LN = 0.D0
            CALL HIGH_PASS_FILTER(CHI)
            CALL HIGH_PASS_FILTER(B)
        ENDIF
        
        ! --- 1. Process Velocity Fields ---
        IF (JOB_X == 0 .OR. JOB_X == 1) THEN 
            BLOCK
                TYPE(SCALAR) :: RUR, RUP, UZ
                IF (MPI_RANK.EQ.0) WRITE(*,*) 'Computing velocity fields for snapshot', I, '...'
                CALL ALLOCATE(RUR); CALL ALLOCATE(RUP); CALL ALLOCATE(UZ)
                CALL CHOPSET(3)
                CALL PC2VEL(PSI, CHI, RUR, RUP, UZ)

                ! Save Radial Velocity (velR)
                IF (JOB_Y == 0 .OR. JOB_Y == 1) THEN
                    CALL SAVE_FIELD('velR', I, RUR, .TRUE.)
                ENDIF
                ! Save Azimuthal Velocity (velTh)
                IF (JOB_Y == 0 .OR. JOB_Y == 2) THEN
                    CALL SAVE_FIELD('velTh', I, RUP, .TRUE.)
                ENDIF
                ! Save Axial Velocity (velZ)
                IF (JOB_Y == 0 .OR. JOB_Y == 3) THEN
                     CALL SAVE_FIELD('velZ', I, UZ, .FALSE.)
                ENDIF
                
                CALL CHOPSET(-3)
                CALL DEALLOCATE(RUR); CALL DEALLOCATE(RUP); CALL DEALLOCATE(UZ)
            END BLOCK
        ENDIF

        ! --- 2. Process Vorticity Fields ---
        IF (JOB_X == 0 .OR. JOB_X == 2) THEN 
            BLOCK
                TYPE(SCALAR) :: ROR, ROP, OZ
                IF (MPI_RANK.EQ.0) WRITE(*,*) 'Computing vorticity fields for snapshot', I, '...'
                CALL ALLOCATE(ROR); CALL ALLOCATE(ROP); CALL ALLOCATE(OZ)
                CALL CHOPSET(3)
                CALL PC2VOR(PSI, CHI, ROR, ROP, OZ)

                ! Save Radial Vorticity (vorR)
                IF (JOB_Y == 0 .OR. JOB_Y == 1) THEN
                    CALL SAVE_FIELD('vorR', I, ROR, .TRUE.)
                ENDIF
                ! Save Azimuthal Vorticity (vorTh)
                IF (JOB_Y == 0 .OR. JOB_Y == 2) THEN
                    CALL SAVE_FIELD('vorTh', I, ROP, .TRUE.)
                ENDIF
                ! Save Axial Vorticity (vorZ)
                IF (JOB_Y == 0 .OR. JOB_Y == 3) THEN
                    CALL SAVE_FIELD('vorZ', I, OZ, .FALSE.)
                ENDIF
                
                CALL CHOPSET(-3)
                CALL DEALLOCATE(ROR); CALL DEALLOCATE(ROP); CALL DEALLOCATE(OZ)
            END BLOCK
        ENDIF

        ! --- 3. Process Buoyancy Field ---
        ! This logic saves buoyancy if requested (JOB_Y 0 or 4)
        ! and any valid field type (JOB_X 0, 1, or 2) was selected.
        IF (JOB_Y == 0 .OR. JOB_Y == 4) THEN
            IF (JOB_X >= 0 .AND. JOB_X <= 2) THEN
                IF (MPI_RANK.EQ.0) WRITE(*,*) 'Computing buoyancy field for snapshot', I, '...'
                ! Buoyancy is loaded as 'B' and is already in spectral space.
                ! We pass .FALSE. for DO_DIVR as it doesn't need 1/r scaling.
                ! NOTE: This modifies 'B' in-place to save memory.
                CALL SAVE_FIELD('buoyancy', I, B, .FALSE.)
            ENDIF
        ENDIF

        ! --- 4. Handle Invalid JOB_X ---
        IF ((JOB_X > 2) .OR. (JOB_X < 0)) THEN
            IF (MPI_RANK.EQ.0) WRITE(*,*) 'POSTPROCESS: Invalid JOB field type (X digit).'
        ENDIF

        CALL DEALLOCATE(PSI); CALL DEALLOCATE(CHI); CALL DEALLOCATE(B)
        IF (POSTPROCESS%START(JOB_IND) .EQ. -1) EXIT

    ENDDO
ENDDO

!> final printout
IF (MPI_RANK.eq.0) THEN
    WRITE(*,*) 'PROGRAM ENDED'
    CALL PRINT_REAL_TIME()
    time_end = mpi_wtime()
    WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================
LOGICAL FUNCTION FILE_EXISTS(FILENAME)
!=======================================================================
!> @brief checks if a file exists on the filesystem.
!> @param FILENAME (in) full path to the file to check
!> @retval FILE_EXISTS .true. if file exists, .false. otherwise
!=======================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    
    INQUIRE(FILE=FILENAME, EXIST=FILE_EXISTS)
    
END FUNCTION FILE_EXISTS
!=======================================================================
SUBROUTINE HIGH_PASS_FILTER(A, CUTOFF_RATIO)
!=======================================================================
!> @brief zeroes out low-frequency spectral coefficients.
!> @param A (inout) scalar field to be filtered.
!> @param CUTOFF_RATIO (in, optional) cutoff ratio, default 2/3.
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR), INTENT(INOUT) :: A
    REAL(P8), INTENT(IN), OPTIONAL :: CUTOFF_RATIO
    
    REAL(P8) :: CUTOFF
    INTEGER :: NN, MM, KK, NR_CUTOFF, NTH_CUTOFF, NX_CUTOFF_MIN, NX_CUTOFF_MAX
    INTEGER :: GLOBAL_MM, GLOBAL_KK
    
    IF (.NOT. PRESENT(CUTOFF_RATIO)) THEN
        CUTOFF = 2.0D0/3.0D0
    ELSE
        CUTOFF = CUTOFF_RATIO
    ENDIF
    
    NR_CUTOFF = INT(NRCHOP / 2.0d0)
    NTH_CUTOFF = INT(NTCHOP * CUTOFF)
    NX_CUTOFF_MIN = INT(NXCHOP * CUTOFF)
    NX_CUTOFF_MAX = 2*NXCHOP+1-NX_CUTOFF_MIN
    
    CALL TOFF(A)
    
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(KK, MM, NN, GLOBAL_MM, GLOBAL_KK) SCHEDULE(STATIC)
    DO KK = 1, SIZE(A%E, 3)
        DO MM = 1, SIZE(A%E, 2)
            DO NN = 1, SIZE(A%E, 1)
                
                GLOBAL_MM = MM + A%INTH
                GLOBAL_KK = KK + A%INX
                
                ! IF (NN < NR_CUTOFF .OR. GLOBAL_MM < NTH_CUTOFF .OR. GLOBAL_KK <= NX_CUTOFF_MIN .OR. GLOBAL_KK >= NX_CUTOFF_MAX) THEN
                !     A%E(NN, MM, KK) = (0.0D0, 0.0D0)
                ! ENDIF

                IF ((GLOBAL_MM .NE. 1) .OR. (GLOBAL_KK .NE. 1)) THEN
                    A%E(NN, MM, KK) = (0.0D0, 0.0D0)
                ENDIF
                
            ENDDO
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    
END SUBROUTINE HIGH_PASS_FILTER
!=======================================================================
SUBROUTINE SAVE_FIELD(FIELD_NAME, INDEX, FIELD, DO_DIVR)
!=======================================================================
!> @brief takes a spectral field, transforms it to physical space,
!>        optionally scales it, and saves the requested slice.
!>
!> @note this subroutine modifies FIELD in-place.
!>       at exit, FIELD will be in physical space and,
!>       if DO_DIVR is .true., will have been divided by 'r'.
!>
!> @param FIELD_NAME (in) name of the field (e.g., 'velR', 'vorTh')
!> @param INDEX      (in) snapshot index (or -1 for initial)
!> @param FIELD      (inout) the scalar field (in spectral space)
!> @param DO_DIVR    (in) logical flag. if .true., call DIVR (divide by 'r').
!=======================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FIELD_NAME
    INTEGER, INTENT(IN) :: INDEX
    TYPE(SCALAR), INTENT(INOUT) :: FIELD
    LOGICAL, INTENT(IN) :: DO_DIVR

    CHARACTER(LEN=72) :: OUT_FILENAME
    CHARACTER(LEN=3) :: NUM_STR
    CHARACTER(LEN=15) :: SLICE_TYPE_STR

    IF (INDEX .EQ. -1) THEN
        NUM_STR = 'ini'
    ELSE
        WRITE(NUM_STR, '(I3.3)') INDEX
    ENDIF

    CALL TOFP(FIELD)
    IF (DO_DIVR) THEN
        CALL DIVR(FIELD)
    ENDIF

    ! Determine filename and save
    SELECT CASE(JOB_Z)
    CASE(1) ! r-theta slice
        IF (POSTPROCESS%SLICEINT(JOB_IND) == 999) THEN
            SLICE_TYPE_STR = '_3D_'
        ELSE
            SLICE_TYPE_STR = '_RTplane_'
        ENDIF
        OUT_FILENAME = TRIM(FIELD_NAME)//TRIM(SLICE_TYPE_STR)//NUM_STR//'.dat'
        CALL MSAVE_SLICES_IN_RTHETA_PLANE(FIELD, TRIM(ADJUSTL(FILES%SAVEDIR))//OUT_FILENAME, POSTPROCESS%SLICEINT(JOB_IND))
    CASE(2) ! r-z slice
        IF (POSTPROCESS%SLICEINT(JOB_IND) == 999) THEN
            SLICE_TYPE_STR = '_3D_'
        ELSE
            SLICE_TYPE_STR = '_RZplane_'
        ENDIF
        OUT_FILENAME = TRIM(FIELD_NAME)//TRIM(SLICE_TYPE_STR)//NUM_STR//'.dat'
        CALL MSAVE_SLICES_IN_RZ_PLANE(FIELD, TRIM(ADJUSTL(FILES%SAVEDIR))//OUT_FILENAME, POSTPROCESS%SLICEINT(JOB_IND))
    CASE DEFAULT
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'SAVE_FIELD: Invalid JOB_Z:', JOB_Z
    END SELECT

    CALL MPI_BARRIER(MPI_COMM_IVP, IERR)
    
END SUBROUTINE SAVE_FIELD
!=======================================================================
SUBROUTINE MSAVE_SLICES_IN_RTHETA_PLANE(A,FILENAME,ZPLANE)
!=======================================================================
!> @brief assembles and saves r-theta slice(s) of a 3d field.
!> @param A (in) scalar field in ppp_space (physical space).
!> @param FILENAME (in) output file name.
!> @param ZPLANE (in) z-index of the plane to save (999 for 3d).
!=======================================================================
TYPE(SCALAR):: A
CHARACTER(LEN=*):: FILENAME
INTEGER:: ZPLANE

INTEGER:: MM,KK
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: A_GLB

IF(A%SPACE.NE.PPP_SPACE) THEN
    IF (MPI_RANK.EQ.0) WRITE(*,*) 'MSAVE_SLICES_IN_RTHETA_PLANE: INPUT NOT IN PPP SPACE'
    CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)
ENDIF

ALLOCATE(A_GLB(NDIMR,NDIMTH,NDIMX))
CALL MASSEMBLE(A%E,A_GLB,2)

IF (MPI_RANK.EQ.0) THEN
    IF (ZPLANE.EQ.999) THEN
        CALL MSAVE(A_GLB(:NR,:NTH,:NX), FILENAME)
    ELSE IF (ZPLANE > 0 .AND. ZPLANE <= NX) THEN
        CALL MSAVE(A_GLB(:NR,:NTH,ZPLANE), FILENAME)
    ELSE
        WRITE(*,*) 'MSAVE_SLICES_IN_RTHETA_PLANE: Invalid ZPLANE index:', ZPLANE
        WRITE(*,*) 'Valid range is 1 to', NX, 'or 999 for 3D.'
    ENDIF
ENDIF
CALL MPI_BARRIER(MPI_COMM_IVP,IERR)
DEALLOCATE(A_GLB)

RETURN
END SUBROUTINE MSAVE_SLICES_IN_RTHETA_PLANE
!=======================================================================
SUBROUTINE MSAVE_SLICES_IN_RZ_PLANE(A,FILENAME,THETAPLANE)
!=======================================================================
!> @brief assembles and saves r-z slice(s) of a 3d field.
!> @param A (in) scalar field in ppp_space (physical space).
!> @param FILENAME (in) output file name.
!> @param THETAPLANE (in) theta-index of the plane to save (999 for 3d).
!=======================================================================
TYPE(SCALAR):: A
CHARACTER(LEN=*):: FILENAME
INTEGER:: THETAPLANE

INTEGER:: MM,KK
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: A_GLB

IF(A%SPACE.NE.PPP_SPACE) THEN
    IF (MPI_RANK.EQ.0) WRITE(*,*) 'MSAVE_SLICES_IN_RZ_PLANE: INPUT NOT IN PPP SPACE'
    CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)
ENDIF

ALLOCATE(A_GLB(NDIMR,NDIMTH,NDIMX))
CALL MASSEMBLE(A%E,A_GLB,2)

IF (MPI_RANK.EQ.0) THEN
    IF (THETAPLANE.EQ.999) THEN
        CALL MSAVE(A_GLB(:NR,:NTH,:NX), FILENAME)
    ELSE IF (THETAPLANE > 0 .AND. THETAPLANE <= NTH) THEN
        CALL MSAVE(A_GLB(:NR,THETAPLANE,:NX), FILENAME)
    ELSE
        WRITE(*,*) 'MSAVE_SLICES_IN_RZ_PLANE: Invalid THETAPLANE index:', THETAPLANE
        WRITE(*,*) 'Valid range is 1 to', NTH, 'or 999 for 3D.'
    ENDIF
ENDIF
CALL MPI_BARRIER(MPI_COMM_IVP,IERR)
DEALLOCATE(A_GLB)

RETURN

END SUBROUTINE MSAVE_SLICES_IN_RZ_PLANE
!=======================================================================
SUBROUTINE DIVR(A)
!=======================================================================
IMPLICIT NONE
TYPE(SCALAR), INTENT(INOUT) :: A
INTEGER :: NN, MM, KK, INR, NN_BOUND

IF (A%SPACE.NE.PPP_SPACE) THEN
    CALL MPRINT('DIVR: INPUT NOT IN PPP SPACE')
    CALL MPI_ABORT(MPI_COMM_IVP,ERR_FLAGS%SCALAR3,IERR)
ENDIF

INR = A%INR
NN_BOUND = MIN(SIZE(A%E,1), NR - INR)

!$OMP PARALLEL DO COLLAPSE(3) 
DO KK = 1, SIZE(A%E,3)
    DO MM = 1, SIZE(A%E,2)
        DO NN = 1, NN_BOUND
            A%E(NN,MM,KK) = A%E(NN,MM,KK) / TFM%R(NN+INR)
        ENDDO
    ENDDO
ENDDO
!$OMP END PARALLEL DO

RETURN

END SUBROUTINE DIVR
!=======================================================================
END PROGRAM POSTPROC
!=======================================================================