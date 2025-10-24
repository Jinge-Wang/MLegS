MODULE MOD_FFT ! LEVEL 2.5 MODULE
! FFT ROUTINES & MPI OPERATIONS
    USE OMP_LIB
    USE MPI
    USE MOD_MISC, ONLY : P4,P8,PI,IU,ITOA3,CP8_SIZE,MSAVE              ! LEVEL 0
    USE MOD_EIG                                                        ! LEVEL 1
    USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
    USE MOD_SCALAR3                                                    ! LEVEL 2
    IMPLICIT NONE
    PRIVATE
!==============================================================================
! TRANSFORM STRATEGY SUMMARY:
!==============================================================================
!
! DOMAIN DECOMPOSITION:
! ---------------------
! - N1: Number of processors in SUBCOMM_1 (typically larger dimension)
! - N2: Number of processors in SUBCOMM_2 (typically smaller dimension)
!
! SPACE DEFINITIONS:
! -----------------
! - PPP: Physical space in (R,th,Z)      [NDIMR/N1, NDIMTH, NDIMX/N2]
! - PFP: Fourier in th only              [NDIMR/N1, NDIMTH, NDIMX/N2]
! - PFF: Fourier in th and Z             [NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1]
! - FFF: Spectral in all directions     [NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1]
!
! TRANSFORMATION SEQUENCES:
! -----------------------
! 1. PHYSICAL TO SPECTRAL: PPP -> PFP -> PFF -> FFF
!    a) HORFFT(A,-1): th-direction FFT transforms PPP -> PFP
!    b) VERFFT(A,-1): Z-direction FFT transforms PFP -> PFF
!    c) RTRAN(A,-1):  R-direction Legendre transform PFF -> FFF
!    d) Alternative: TOFF(A) performs all steps in sequence
!
! 2. SPECTRAL TO PHYSICAL: FFF -> PFF -> PFP -> PPP
!    a) RTRAN(A,1):  R-direction inverse Legendre transform FFF -> PFF
!    b) VERFFT(A,1): Z-direction inverse FFT transforms PFF -> PFP
!    c) HORFFT(A,1): th-direction inverse FFT transforms PFP -> PPP
!    d) Alternative: TOFP(A) performs all steps in sequence
!
! INTERMITTENT MPI DATATYPES:
! --------------------------
! During VERFFT, several specialized MPI datatypes facilitate efficient 
! parallel data exchange:
!
! 1. TYPE_PFP0: [NDIMR/N1, NTCHOPDIM, NDIMX/N2]
!    - Used in VERFFT for initial data layout where R is distributed across SUBCOMM_1 
!      and Z is distributed across SUBCOMM_2
!    - Used as source type in EXCHANGE_3DCOMPLEX_FAST to redistribute data
!
! 2. TYPE_PFP1: [NDIMR/N1, NTCHOPDIM/N2, NDIMX]
!    - Target format after first redistribution in VERFFT
!    - Consolidates Z dimension locally while distributing theta across processors
!    - Enables efficient Z-direction FFT with contiguous memory access
!
! 3. TYPE_PFF0: [NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM]
!    - Format after Z-direction FFT and chopping in VERFFT
!    - Used as source type for second redistribution
!    - Prepares data for radial dimension consolidation
!
! 4. TYPE_PFF1: [NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1]
!    - Final distributed format before Legendre transform
!    - Allows RTRAN to operate on complete radial dimension locally
!
! VERFFT TRANSFORMATION PROCESS IN DETAIL:
! ---------------------------------------
! 1. Forward transform (PPP->PFF, IS=-1):
!    a) Input: A%E in PFP_SPACE format [NDIMR/N1, NDIMTH, NDIMX/N2]
!    b) Extract chopped theta data into A_PFP0 [NDIMR/N1, NTCHOPDIM, NDIMX/N2]
!    c) Call EXCHANGE_3DCOMPLEX_FAST with TYPE_PFP0/TYPE_PFP1 to redistribute
!       producing A_PFP1 [NDIMR/N1, NTCHOPDIM/N2, NDIMX]
!    d) Perform Z-direction FFT using DFFTW_EXECUTE_DFT on each processor
!    e) Chop to spectral resolution, forming A_PFF0 [NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM]
!    f) Call EXCHANGE_3DCOMPLEX_FAST with TYPE_PFF0/TYPE_PFF1 for final redistribution
!       producing A_PFF1 [NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1]
!    g) Reallocate A in PFF_SPACE and copy A_PFF1 data
!
! 2. Backward transform (PFF->PFP, IS=1):
!    a) Input: A%E in PFF_SPACE format [NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1]
!    b) Store data in A_PFF1 [NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1]
!    c) Call EXCHANGE_3DCOMPLEX_FAST with TYPE_PFF1/TYPE_PFF0 to redistribute
!       producing A_PFF0 [NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM]
!    d) Expand to full Z dimension in A_PFP1 [NDIMR/N1, NTCHOPDIM/N2, NDIMX]
!    e) Perform inverse Z-direction FFT using DFFTW_EXECUTE_DFT
!    f) Call EXCHANGE_3DCOMPLEX_FAST with TYPE_PFP1/TYPE_PFP0 to redistribute
!       producing A_PFP0 [NDIMR/N1, NTCHOPDIM, NDIMX/N2]
!    g) Reallocate A in PFP_SPACE and copy expanded data
!
! The EXCHANGE_3DCOMPLEX_FAST function utilizes MPI_ALLTOALLW with pre-created
! datatypes to redistribute data with minimal memory copies.
!
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================

! =========================== FFT UTILITIES ============================
! INITIALIZE ALL SPECTRAL METHOD PARAMETERS AND SET THE TFM KIT:
PUBLIC:: SETUP_GRID, LEGINIT, PRINT_COLLOC_INFO
! ALLOCATE/DEALLOCATE MEMORY TO VARIABLES OF "TYPE(SCALAR)":
PUBLIC:: ALLOCATE, DEALLOCATE
! WRAPPER FOR "EQUAL" OPERATOR (COPY0) AMONG TYPE(SCALAR) VARS.:
PUBLIC:: ASSIGNMENT(=)
! PERFORM THE FFT IN THETA (=HORFFT) AND Z (=VERFFT) DIRECTIONS:
PUBLIC:: HORFFT
PUBLIC:: VERFFT
PUBLIC:: PLNFFT
! PERFORM MLEGS TRANSFORMATION
PUBLIC:: RTRAN
PUBLIC:: TOFF,TOFP
    
! =========================== MPI UTILITIES ============================
! CREATE CARTESIAN SUB-GRID COMMUNICATORS:
PUBLIC:: SUBCOMM_CART, DECOMPOSE
! ASSEMBLE LOCAL ARRAY INTO GLOBAL ARRAY:
PUBLIC:: MASSEMBLE, MDISASSEMBLE
! LOAD/SAVE A MATRIX WITH TYPE(SCALAR) FROM/INTO A SPECIFIED FILE
PUBLIC:: MLOAD, MSAVE
PUBLIC:: EXCHANGE_3DCOMPLEX_FAST
PUBLIC:: PRINT_MPI_STRATEGY, MPRINT

! ========================== MOD_FFT PRIVATE ===========================
! ! DECOMPOSE A 1D DOMAIN:
! PUBLIC:: DECOMPOSE
! ! EXCHANGE ARRAY DIMENSIONS:
! PUBLIC:: EXCHANGE
! ! CREATE SUBARRAY DATATYPE FOR EACH PROC THAT SPLITS THE ORIGINAL AR
! ! RAY INTO NPROC PRATS ALONG SPECIFIED DIM:
! PUBLIC:: SUBARRAY

! =========================== UTILITY FUNCS ============================
PUBLIC:: local_size, local_index, local_proc, count_proc

!=======================================================================
!============================ INTERFACES ===============================
!=======================================================================
    INTERFACE ALLOCATE
        MODULE PROCEDURE SALLOC
    END INTERFACE

    INTERFACE DEALLOCATE
        MODULE PROCEDURE SFREE
    END INTERFACE

    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE COPY0
    END INTERFACE

    INTERFACE EXCHANGE
        MODULE PROCEDURE EXCHANGE_3DCOMPLEX
        MODULE PROCEDURE EXCHANGE_3DCOMPLEX_FAST
    END INTERFACE

    INTERFACE MLOAD
    MODULE PROCEDURE MLOAD0
    END INTERFACE

    INTERFACE MSAVE
        MODULE PROCEDURE MSAVE0
    END INTERFACE

CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
SUBROUTINE SETUP_GRID(SAVEDIR)
!=======================================================================
! [USAGE]:
! INITIALIZE PHYSICAL AND SPECTRAL DOMAINS AND TRANSFORMATION KITS
! [PARAMETERS]:
! SAVEDIR >> (OPTIONAL) CHARACTER STRING FOR OUTPUT DIRECTORY
!=======================================================================
IMPLICIT NONE
CHARACTER(LEN=*), OPTIONAL :: SAVEDIR

CALL LEGINIT()
IF (PRESENT(SAVEDIR)) CALL PRINT_MPI_STRATEGY(SAVEDIR)

END SUBROUTINE SETUP_GRID
!=======================================================================
SUBROUTINE LEGINIT(MPI_COMM_INPUT, M_INPUT)
!=======================================================================
! [USAGE]: 
! INITIALIZE ALL SPECTRAL METHOD PARAMETERS AND 
! SET THE TRANSFORM KIT VIA TFM VARIABLE
! [PARAMETERS]:
! NRIN >> INTEGER VALUE USED FOR NR (radial)
! NTHIN >> INTEGER VALUE USED FOR NTH
! NXIN >> INTEGER VALUE USED FOR NX (axial)
! NRCHOPIN >> INTEGER VALUE USED FOR NRCHOP
! NTCHOPIN >> INTEGER VALUE USED FOR NTCHOP
! NXCHOPIN >> INTEGER VALUE USED FOR NXCHOP
! ZLENIN >> REAL VALUE USED FOR ZLEN
! ELLIN >> REAL VALUE USED FOR ELL
! MKLINKIN >> (OPTIONAL) INTEGER VALUE USED FOR MKLINK (DEFAULT: 0)
! MINCIN >> (OPTIONAL) INTEGER VALUE USED FOR MKLINK (DEFAULT: 1)
! [NOTE]:
! SEE ANNOTATIONS IN THE 'PARAMETERS' SECTION OF THIS MODULE TO
! COMPREHEND WHICH INPUT ASSIGNED TO WHICH PUBLIC PARAMETER AND WHAT
! THE MEANING OF EACH PARAMETER IS
! [DEPENDENCIES]:
! 1. LEG_ZERO(~)  @ MOD_LIN_LEGENDRE
! 2. LEG_NORM(~)  @ MOD_LIN_LEGENDRE
! 3. LEG_TBL(~)   @ MOD_LIN_LEGENDRE
! 4. EXSET(~)     @ MOD_SCALAR3
! [NOTES]:
! SUBCOMM_CART WILL REODER THE PROC INDEX
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    INTEGER :: MM,KK,KV,I
    INTEGER, DIMENSION(:), ALLOCATABLE:: SUB_GROUPS

    INTEGER, OPTIONAL:: M_INPUT,MPI_COMM_INPUT

    IF (PRESENT(MPI_COMM_INPUT)) THEN
        MPI_COMM_IVP = MPI_COMM_INPUT
    ELSE
        MPI_COMM_IVP = MPI_COMM_WORLD
    ENDIF
    CALL MPI_COMM_RANK(MPI_COMM_IVP,MPI_RANK,IERR)

    NRH = NR/2
    ELL2 = ELL**2.D0

    NDIMR  = NR+3
    NDIMTH = NTH+1
    NDIMX  = NX
    IF(NDIMTH.EQ.2) NDIMTH=1

    ! NRCHOPDIM  = MIN(NRCHOP+4,NR)       ! For Init_ene
    NRCHOPDIM  = MIN(NRCHOP+6,NR)         ! Original: +3
    NTCHOPDIM  = MIN(NTCHOP,NTH)
    NXCHOPDIM  = NXCHOP*2-1
    NXCHOPH    = NX-NXCHOP+2                                           ! CHOP LOCATION FOR CONJUGATE SIDE
    IF (NXCHOP==1) NXCHOPH=1

    IF ((2*NXCHOP>NX).AND.((NX>1).OR.(NXCHOP>1))) THEN
        IF (MPI_RANK.EQ.0) THEN
            WRITE(*,*) 'ERROR: SCALAR3() -- NXCHOP MUST BE < NX/2,'
            WRITE(*,*) '                    UNLESS NX = NXCHOP = 1'
        ENDIF
      STOP
    ENDIF

    IF (.NOT.(ALLOCATED(M))) ALLOCATE(M(NTCHOPDIM))
    M = (/ (MINC*(MM-1),MM=1,NTCHOPDIM) /)
    IF (PRESENT(M_INPUT).AND.(NTCHOPDIM.EQ.2)) M(2) = M_INPUT

    IF (.NOT.(ALLOCATED(NRCHOPS))) THEN
        ALLOCATE(NRCHOPS(NTH))
        ALLOCATE(NTCHOPS(NR))
    ENDIF
    CALL CHOPSET(0)

    IF (.NOT.(ALLOCATED(AK))) ALLOCATE( AK(NTCHOPDIM,NXCHOPDIM) )
    DO MM=1,NTCHOPDIM
      DO KK=1,NXCHOPDIM
        KV=KK-1 + MKLINK*(MM-1)
        IF(KK.GT.NXCHOP) THEN
          KV=-(NXCHOPDIM-KK+1) ! This needs to be fixed for mklink != 0 case
        ENDIF
        AK(MM,KK) = 2*PI/ZLEN*KV
      ENDDO
    ENDDO

    IF (.NOT.(ASSOCIATED(TFM%R))) THEN
        ALLOCATE( TFM%R(NR) )                                              ! START SETTING UP TRANSFORM PARAMETERS INTO THE TFM VARIABLE
        ALLOCATE( TFM%dR(NR))
        ALLOCATE( TFM%W(NR) )
        ALLOCATE( TFM%X(NR) )
        ALLOCATE( TFM%LN(NR) )
        ALLOCATE( TFM%NORM(NRCHOPDIM+14,NTCHOPDIM) )
        ALLOCATE( TFM%LOGNORM(NRCHOPDIM+14,NTCHOPDIM) )
        ALLOCATE( TFM%PF(NRH,NRCHOPDIM+1,NTCHOPDIM) )
        ALLOCATE( TFM%AT0(NRCHOPDIM) )
        ALLOCATE( TFM%AT1(NRCHOPDIM) )
        ALLOCATE( TFM%THR(MINC*NTH+1) )
        ALLOCATE( TFM%THI(MINC*NTH+1) )
        ALLOCATE( TFM%TH(2*MINC*NTH+1) )
    ENDIF

    CALL LEG_ZERO(TFM%X,TFM%R,ELL,TFM%W)
    TFM%dR(1:NR-1) = TFM%R(2:NR)-TFM%R(1:NR-1)
    TFM%dR(NR) = TFM%dR(NR-1)

    !ALLOCATE( TFM%LN(NR) )
    TFM%LN = -LOG(1-TFM%X) ! -LOG(2L^2/(r^2+L^2)) = LOG((r^2+L^2)/(2L^2)) IS P_L(r) AS IN EQN (54) OF MATSUSHIMA & MARCUS (1997)

    !ALLOCATE( TFM%NORM(NRCHOPDIM+14,NTCHOPDIM) )
    TFM%NORM = LEG_NORM(NRCHOPDIM+14, M )   ! BECOME INACCURATE AT HIGHER SPECTRAL MODES AND SHOULD BE REPLACED BY TFM%LOGNORM INSTEAD

    !ALLOCATE( TFM%LOGNORM(NRCHOPDIM+14,NTCHOPDIM) )
    TFM%LOGNORM = LEG_LOG_NORM(NRCHOPDIM+14, M )

    !ALLOCATE( TFM%PF(NRH,NRCHOPDIM+1,NTCHOPDIM) )
    ! TFM%PF = LEG_TBL(TFM%X(:NRH),NRCHOPDIM+1, M ,TFM%NORM)
    TFM%PF = LEG_TBL(TFM%X(:NRH), NRCHOPDIM+1, M ,TFM%LOGNORM)

    !ALLOCATE( TFM%AT0(NRCHOPDIM) )
    !ALLOCATE( TFM%AT1(NRCHOPDIM) )
    ! TFM%AT0 = LEG_TBL( -1.D0,NRCHOPDIM, 0 , TFM%NORM(:,1))
    ! TFM%AT1 = LEG_TBL(  1.D0,NRCHOPDIM, 0 , TFM%NORM(:,1))
    TFM%AT0 = LEG_TBL( -1.D0, NRCHOPDIM, 0 , TFM%LOGNORM(:,1))
    TFM%AT1 = LEG_TBL(  1.D0, NRCHOPDIM, 0 , TFM%LOGNORM(:,1))

    ! TFM%NMAX = MAX(NTH,NX)*2
    ! ALLOCATE( TFM%EX(TFM%NMAX) )
    ! CALL EXSET( TFM%EX, TFM%NMAX )

    IF(MKLINK.EQ.0) THEN
        IF (.NOT.(ASSOCIATED(TFM%Z))) ALLOCATE( TFM%Z(NX+1) )
        TFM%Z = (/ (ZLEN*(KK-1)/NX, KK=1,NX+1) /)
    ELSE IF(MKLINK.EQ.1 .OR. MKLINK.EQ.-1) THEN
        IF (.NOT.(ASSOCIATED(TFM%Z))) ALLOCATE( TFM%Z(2*NTH+1) )
        TFM%Z = (/ (ZLEN*(KK-1)/(2*NTH), KK=1,2*NTH+1) /)
    ELSE
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'LEGINIT: INVALID MKLINK. MKLINK=',MKLINK
      STOP
    ENDIF

    ! ALLOCATE( TFM%THR(MINC*NTH+1) )
    ! ALLOCATE( TFM%THI(MINC*NTH+1) )
    TFM%THR = (/ (2*PI/MINC/NTH*(I-1), I=1,MINC*NTH+1) /)
    TFM%THI = (/ (2*PI/MINC/NTH*(I-0.5D0), I=1,MINC*NTH+1) /)

    ! ALLOCATE( TFM%TH(2*MINC*NTH+1) )
    TFM%TH = (/ (PI/MINC/NTH*(I-1), I=1,2*MINC*NTH+1) /)

    ! FORM MPI COMMUNICATOR GROUPS:
    IF (PRESENT(M_INPUT)) THEN ! FREE PREVIOUSLY USED COMMUNICATORS
        IF (SUBCOMM_1.NE.0) CALL MPI_COMM_FREE(SUBCOMM_1,IERR)
        IF (SUBCOMM_2.NE.0) CALL MPI_COMM_FREE(SUBCOMM_2,IERR)
    ENDIF
    CALL SUBCOMM_CART(MPI_COMM_IVP, 2, SUB_GROUPS)
! ============= FROM THIS POINT, PROC INDEX IS REORDERED ===============
    IF (NXCHOPDIM.GE.NTCHOPDIM) THEN
        SUBCOMM_1 = SUB_GROUPS(1) ! FFF: NX/N1
        SUBCOMM_2 = SUB_GROUPS(2) ! FFF: NT/N2
    ELSE
        SUBCOMM_1 = SUB_GROUPS(2)
        SUBCOMM_2 = SUB_GROUPS(1)       
    ENDIF
    DEALLOCATE(SUB_GROUPS)

    ! FORM MPI SUBARRAY DATA TYPES:
    SIZE_PFP0 = (/local_size(NDIMR,SUBCOMM_1), &
                                    NTCHOPDIM, &
                  local_size(NDIMX,SUBCOMM_2) /)
    TYPE_PFP0 = create_new_type3D(SUBCOMM_2, SIZE_PFP0, 2, MPI_DOUBLE_COMPLEX)

    SIZE_PFP1 = (/local_size(NDIMR,SUBCOMM_1), &
              local_size(NTCHOPDIM,SUBCOMM_2), &
                                        NDIMX /)
    TYPE_PFP1 = create_new_type3D(SUBCOMM_2, SIZE_PFP1, 3, MPI_DOUBLE_COMPLEX)

    SIZE_PFF0 = (/local_size(NDIMR,SUBCOMM_1), &
              local_size(NTCHOPDIM,SUBCOMM_2), &
                                    NXCHOPDIM /)
    TYPE_PFF0 = create_new_type3D(SUBCOMM_1, SIZE_PFF0, 3, MPI_DOUBLE_COMPLEX)

    SIZE_PFF1 = (/                      NDIMR, &
              local_size(NTCHOPDIM,SUBCOMM_2), &
              local_size(NXCHOPDIM,SUBCOMM_1) /)
    TYPE_PFF1 = create_new_type3D(SUBCOMM_1, SIZE_PFF1, 1, MPI_DOUBLE_COMPLEX)

    RETURN
END SUBROUTINE LEGINIT
!=======================================================================
SUBROUTINE PRINT_COLLOC_INFO(SAVEDIR)
!=======================================================================
! [USAGE]:
! SAVE PHYSICAL COLLOCATION POINTS
!=======================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*) :: SAVEDIR

    IF (MPI_RANK.EQ.0) THEN
        ! SAVE THE RADIAL COLLOCATION POINTS (R)
        CALL MSAVE(TFM%R, TRIM(ADJUSTL(SAVEDIR))//&
            'r_colloc_pts.info')

        ! SAVE THE X COLLOCATION POINTS (X)
        CALL MSAVE(TFM%X, TRIM(ADJUSTL(SAVEDIR))//&
            'x_colloc_pts.info')

        ! SAVE THE GAUSS_LEGENDRE WEIGHTS (W)
        CALL MSAVE(TFM%W, TRIM(ADJUSTL(SAVEDIR))//&
            'gau_leg_weights.info')

        ! SAVE THE AZIMUTHAL COLLOCATION POINTS (THETA)
        CALL MSAVE(TFM%TH, TRIM(ADJUSTL(SAVEDIR))//&
            't_colloc_pts.info')

        ! SAVE THE AXIAL COLLOCATION POINTS (Z)
        CALL MSAVE(TFM%Z, TRIM(ADJUSTL(SAVEDIR))//&
            'z_colloc_pts.info')
    ENDIF

    RETURN
END SUBROUTINE PRINT_COLLOC_INFO
!=======================================================================
SUBROUTINE PRINT_MPI_STRATEGY(SAVEDIR)
!=======================================================================
! [USAGE]: 
! PRINT DEBUG INFORMATION SHOWING HOW SCALAR3 VARIABLES ARE DISTRIBUTED
! AMONG MPI PROCESSORS FOR EACH SPACE TYPE (FFF, PFF, PFP, PPP)
! [OUTPUTS]:
! FOR EACH SPACE TYPE:
! - LOCAL ARRAY DIMENSIONS FOR EACH PROCESSOR
! - STARTING INDICES FOR EACH PROCESSOR  
! - GLOBAL ARRAY DIMENSIONS
! - MPI COMMUNICATOR INFORMATION
! [NOTE]:
! OUTPUT IS WRITTEN TO 'savedir/mpi.info' FILE
!=======================================================================
    IMPLICIT NONE

    ! MPI process information
    INTEGER :: TOTAL_PROCS, N1, N2
    
    ! Variables for testing and looping
    CHARACTER(LEN=20), DIMENSION(4) :: SPACE_NAME_LIST = ['PPP_SPACE', 'PFP_SPACE', 'PFF_SPACE', 'FFF_SPACE']
    INTEGER :: SPACE_ID, II, RANK
    CHARACTER(LEN=20) :: SPACE_NAME
    TYPE(SCALAR) :: TEST_SCALAR
    
    ! File I/O variables
    INTEGER :: MPI_INFO_UNIT = 99, STATUS
    CHARACTER(LEN=200) :: MPI_INFO_FILENAME
    CHARACTER(LEN=*), OPTIONAL :: SAVEDIR

    IF (.NOT. PRESENT(SAVEDIR)) THEN
        MPI_INFO_FILENAME = 'mpi.info'
    ELSE
        MPI_INFO_FILENAME = TRIM(ADJUSTL(SAVEDIR)) // 'mpi.info'
    ENDIF

    CALL MPI_COMM_SIZE(MPI_COMM_IVP, TOTAL_PROCS, IERR)
    N1 = count_proc(SUBCOMM_1)
    N2 = count_proc(SUBCOMM_2)
    IF (MPI_RANK == 0) THEN
        OPEN(UNIT=MPI_INFO_UNIT, FILE=TRIM(MPI_INFO_FILENAME), STATUS='UNKNOWN', &
            FORM='FORMATTED', IOSTAT=STATUS)
        
        IF (STATUS /= 0) THEN
            WRITE(*,*) 'PRINT_MPI_STRATEGY: Failed to open file: ', TRIM(MPI_INFO_FILENAME)
            WRITE(*,*) 'Writing to stdout instead...'
            MPI_INFO_UNIT = 6  ! Use stdout as fallback
        ELSE
            WRITE(*,*) 'MPI strategy written to: ', TRIM(MPI_INFO_FILENAME)
        ENDIF

        WRITE(MPI_INFO_UNIT,*) '=========================================================='
        WRITE(MPI_INFO_UNIT,*) '              MPI SUBARRAY DISTRIBUTION DEBUG             '
        WRITE(MPI_INFO_UNIT,*) '=========================================================='
        WRITE(MPI_INFO_UNIT,*) 'Total MPI Processes:', TOTAL_PROCS
        WRITE(MPI_INFO_UNIT,*) 'Processor Grid: N1 x N2 =', N1, 'x', N2
        WRITE(MPI_INFO_UNIT,*) 'SUBCOMM_1 size:', N1, ' (typically for X/Z dimension)'
        WRITE(MPI_INFO_UNIT,*) 'SUBCOMM_2 size:', N2, ' (typically for Theta dimension)'
        WRITE(MPI_INFO_UNIT,*) ''
        WRITE(MPI_INFO_UNIT,*) 'Global Dimensions:'
        WRITE(MPI_INFO_UNIT,*) '  NDIMR    =', NDIMR,    '  NDIMTH   =', NDIMTH,   '  NDIMX     =', NDIMX
        WRITE(MPI_INFO_UNIT,*) '  NRCHOPDIM=', NRCHOPDIM,'  NTCHOPDIM=', NTCHOPDIM,'  NXCHOPDIM =', NXCHOPDIM
        WRITE(MPI_INFO_UNIT,*) '=========================================================='
        
        IF (STATUS == 0) CLOSE(MPI_INFO_UNIT)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_IVP, IERR)

    DO II = 1, SIZE(SPACE_NAME_LIST)
        
        SPACE_NAME = SPACE_NAME_LIST(II)
        SPACE_ID = SPACE_MAPPER%GET_ID(SPACE_NAME)        
        CALL ALLOCATE(TEST_SCALAR, SPACE_ID)
        
        DO RANK = 0, TOTAL_PROCS - 1
            IF (MPI_RANK == RANK) THEN
                IF (RANK == 0) THEN
                    OPEN(UNIT=MPI_INFO_UNIT, FILE=TRIM(MPI_INFO_FILENAME), STATUS='OLD', &
                    POSITION='APPEND', FORM='FORMATTED', IOSTAT=STATUS)
                    
                    IF (STATUS == 0) THEN
                    WRITE(MPI_INFO_UNIT,*) ''
                    WRITE(MPI_INFO_UNIT,*) '----------------------------------------------------------'
                    WRITE(MPI_INFO_UNIT,*) TRIM(SPACE_NAME), ' DISTRIBUTION:'
                    WRITE(MPI_INFO_UNIT,*) '----------------------------------------------------------'
                    
                    ! Use the actual space ID for comparison
                    IF ((SPACE_ID == PPP_SPACE).OR.(SPACE_ID == PFP_SPACE)) THEN
                        WRITE(MPI_INFO_UNIT,*) 'Layout: [NDIMR/N1, NDIMTH, NDIMX/N2]'
                        WRITE(MPI_INFO_UNIT,*) 'Global: [', NDIMR, ',', NDIMTH, ',', NDIMX, ']'
                    ELSEIF (SPACE_ID == PFF_SPACE) THEN
                        WRITE(MPI_INFO_UNIT,*) 'Layout: [NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1]'
                        WRITE(MPI_INFO_UNIT,*) 'Global: [', NDIMR, ',', NTCHOPDIM, ',', NXCHOPDIM, ']'
                    ELSEIF (SPACE_ID == FFF_SPACE) THEN
                        WRITE(MPI_INFO_UNIT,*) 'Layout: [NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1]'
                        WRITE(MPI_INFO_UNIT,*) 'Global: [', NRCHOPDIM, ',', NTCHOPDIM, ',', NXCHOPDIM, ']'
                    ENDIF

                    WRITE(MPI_INFO_UNIT,*) ''
                    WRITE(MPI_INFO_UNIT,'(A4,A3,A19,A3,A19,A3,A18)') &
                        'Rank', ' | ', 'Local_Size(1,2,3)', ' | ', 'Start_Index(1,2,3)', ' | ', 'SUBCOMM_Ranks(1,2)'
                    WRITE(MPI_INFO_UNIT,'(A4,A3,A19,A3,A19,A3,A18)') &
                        '----', '-|-', '-------------------', '-|-', '-------------------', '-|-', '------------------'
                    CLOSE(MPI_INFO_UNIT)
                    ENDIF
                ENDIF
                
                OPEN(UNIT=MPI_INFO_UNIT, FILE=TRIM(MPI_INFO_FILENAME), STATUS='OLD', &
                    POSITION='APPEND', FORM='FORMATTED', IOSTAT=STATUS)
                
                IF (STATUS == 0) THEN
                    IF ((SPACE_ID == PPP_SPACE).OR.(SPACE_ID == PFP_SPACE)) THEN
                        WRITE(MPI_INFO_UNIT,'(I4,A3,A19,A3,A19,A3,A18)') &
                            MPI_RANK, ' | ', &
                            ADJUSTR('(' // TRIM(ITOA3(SIZE(TEST_SCALAR%E,1))) // ',' // &
                                         TRIM(ITOA3(SIZE(TEST_SCALAR%E,2))) // ',' // &
                                         TRIM(ITOA3(SIZE(TEST_SCALAR%E,3))) // ')'), ' | ', &
                            ADJUSTR('(' // TRIM(ITOA3(TEST_SCALAR%INR)) // ',' // &
                                         TRIM(ITOA3(0)) // ',' // &
                                         TRIM(ITOA3(TEST_SCALAR%INX)) // ')'), ' | ', &
                            ADJUSTR('(' // TRIM(ITOA3(local_proc(SUBCOMM_1))) // ',' // &
                                         TRIM(ITOA3(local_proc(SUBCOMM_2))) // ')')
                    ELSEIF ((SPACE_ID == PFF_SPACE).OR.(SPACE_ID == FFF_SPACE)) THEN
                        WRITE(MPI_INFO_UNIT,'(I4,A3,A19,A3,A19,A3,A18)') &
                            MPI_RANK, ' | ', &
                            ADJUSTR('(' // TRIM(ITOA3(SIZE(TEST_SCALAR%E,1))) // ',' // &
                                         TRIM(ITOA3(SIZE(TEST_SCALAR%E,2))) // ',' // &
                                         TRIM(ITOA3(SIZE(TEST_SCALAR%E,3))) // ')'), ' | ', &
                            ADJUSTR('(' // TRIM(ITOA3(0)) // ',' // &
                                         TRIM(ITOA3(TEST_SCALAR%INTH)) // ',' // &
                                         TRIM(ITOA3(TEST_SCALAR%INX)) // ')'), ' | ', &
                            ADJUSTR('(' // TRIM(ITOA3(local_proc(SUBCOMM_1))) // ',' // &
                                         TRIM(ITOA3(local_proc(SUBCOMM_2))) // ')')
                    ENDIF
                    CLOSE(MPI_INFO_UNIT)
                ENDIF
            
            ENDIF
            
            CALL MPI_BARRIER(MPI_COMM_IVP, IERR)

        ENDDO
        
        IF (MPI_RANK == 0) THEN
            OPEN(UNIT=MPI_INFO_UNIT, FILE=TRIM(MPI_INFO_FILENAME), STATUS='OLD', &
                POSITION='APPEND', FORM='FORMATTED', IOSTAT=STATUS)
            
            IF (STATUS == 0) THEN
                WRITE(MPI_INFO_UNIT,*) '----------------------------------------------------------'
                
                IF ((SPACE_ID == PPP_SPACE).OR.(SPACE_ID == PFP_SPACE)) THEN
                    WRITE(MPI_INFO_UNIT,*) 'Summary: R-dimension distributed across', N1, 'processors'
                    WRITE(MPI_INFO_UNIT,*) '         Theta-dimension NOT distributed (all procs have full)'
                    WRITE(MPI_INFO_UNIT,*) '         X-dimension distributed across', N2, 'processors'
                ELSEIF ((SPACE_ID == PFF_SPACE).OR.(SPACE_ID == FFF_SPACE)) THEN
                    WRITE(MPI_INFO_UNIT,*) 'Summary: R-dimension NOT distributed (all procs have full)'
                    WRITE(MPI_INFO_UNIT,*) '         Theta-dimension distributed across', N2, 'processors'
                    WRITE(MPI_INFO_UNIT,*) '         X-dimension distributed across', N1, 'processors'
                ENDIF
                CLOSE(MPI_INFO_UNIT)
            ENDIF
        ENDIF
        
        CALL DEALLOCATE(TEST_SCALAR)        
        CALL MPI_BARRIER(MPI_COMM_IVP, IERR)
    ENDDO

    ! Print transformation flow (only from rank 0)
    ! And save collocation info
    IF (MPI_RANK == 0) THEN

        OPEN(UNIT=MPI_INFO_UNIT, FILE=TRIM(MPI_INFO_FILENAME), STATUS='OLD', &
            POSITION='APPEND', FORM='FORMATTED', IOSTAT=STATUS)
        
        IF (STATUS == 0) THEN
            WRITE(MPI_INFO_UNIT,*) ''
            WRITE(MPI_INFO_UNIT,*) '=========================================================='
            WRITE(MPI_INFO_UNIT,*) '                   TRANSFORMATION FLOW                   '
            WRITE(MPI_INFO_UNIT,*) '=========================================================='
            WRITE(MPI_INFO_UNIT,*) 'PHYSICAL -> SPECTRAL:'
            WRITE(MPI_INFO_UNIT,*) '  PPP -> PFP -> PFF -> FFF'
            WRITE(MPI_INFO_UNIT,*) '  (HORFFT) (VERFFT) (RTRAN)'
            WRITE(MPI_INFO_UNIT,*) ''
            WRITE(MPI_INFO_UNIT,*) 'SPECTRAL -> PHYSICAL:'
            WRITE(MPI_INFO_UNIT,*) '  FFF -> PFF -> PFP -> PPP'
            WRITE(MPI_INFO_UNIT,*) '  (RTRAN) (VERFFT) (HORFFT)'
            WRITE(MPI_INFO_UNIT,*) ''
            WRITE(MPI_INFO_UNIT,*) 'MPI DATA EXCHANGES occur during:'
            WRITE(MPI_INFO_UNIT,*) '  - VERFFT: Redistribution between PFP<->PFF spaces'
            WRITE(MPI_INFO_UNIT,*) '  - Uses TYPE_PFP0, TYPE_PFP1, TYPE_PFF0, TYPE_PFF1 datatypes'
            WRITE(MPI_INFO_UNIT,*) '=========================================================='
            CLOSE(MPI_INFO_UNIT)
        ENDIF

        CALL PRINT_COLLOC_INFO(SAVEDIR)

    ENDIF

    RETURN
END SUBROUTINE PRINT_MPI_STRATEGY
!=======================================================================
SUBROUTINE SALLOC(A,SP)
!=======================================================================
! [USAGE]: 
! ALLOCATE AN ARRAY AT A%E WHERE A IS SCALAR-TYPE VARIABLE WITH SP TAG
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! SP >> (OPTIONAL) SPACE TAG TO BE ATTACHED. IF NONE, FFF IS ASSIGNED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! [NOTE]:
! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
! PFP: NDIMR/N1, NDIMTH, NDIMX/N2
! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1
! FFF: NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    INTEGER,OPTIONAL:: SP

    INTEGER:: MM

    IF (ASSOCIATED(A%E)) THEN
        CALL SFREE(A)
    ENDIF

    IF (PRESENT(SP)) THEN
        IF (SP < 0 .OR. SP > 3) THEN
            WRITE(*,*) "SALLOC: INVALID SPACE TAG - ", SP
            STOP
        ENDIF
        A%SPACE = SP
    ELSE
        A%SPACE = FFF_SPACE
    ENDIF

    ! Initialize all fields
    A%INR = 0
    A%INTH = 0
    A%INX = 0
    A%LN = 0.D0

    SELECT CASE(A%SPACE)
    CASE (PFF_SPACE)
        ALLOCATE( A%E(NDIMR    , local_size(NTCHOPDIM,SUBCOMM_2), local_size(NXCHOPDIM,SUBCOMM_1)),STAT=IERR )
        A%INTH = local_index(NTCHOPDIM,SUBCOMM_2)
        A%INX  = local_index(NXCHOPDIM,SUBCOMM_1)

    CASE (FFF_SPACE)
        ALLOCATE( A%E(NRCHOPDIM, local_size(NTCHOPDIM,SUBCOMM_2), local_size(NXCHOPDIM,SUBCOMM_1)),STAT=IERR )
        A%INTH = local_index(NTCHOPDIM,SUBCOMM_2)
        A%INX  = local_index(NXCHOPDIM,SUBCOMM_1)

    ! PPP_SPACE OR PFP_SPACE:
    CASE DEFAULT 
        ALLOCATE( A%E(local_size(NDIMR,SUBCOMM_1), NDIMTH, local_size(NDIMX,SUBCOMM_2)),STAT=IERR )
        A%INR  = local_index(NDIMR,SUBCOMM_1)
        A%INX  = local_index(NDIMX,SUBCOMM_2)
    END SELECT

    IF (IERR.NE.0) THEN
        WRITE(*,*) "SALLOC: ERROR OCCURED - COULD NOT ALLOCATE ARRAY."
        STOP
    ENDIF
    A%IS_ALLOCATED = .TRUE.

    RETURN
    END SUBROUTINE SALLOC
!=======================================================================
    SUBROUTINE SFREE(A)
!=======================================================================
! [USAGE]: 
! FREE AN ARRAY WHERE A SCLAR-TYPE VARIABLE A INDICATES VIA A%E
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! [NOTES]:
! DEALLOCATE TARGET OF A%E
! NULLIFY POINTER A%E
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR),INTENT(INOUT):: A
    INTEGER:: DEALLOC_STATUS

    IF (.NOT. ASSOCIATED( A%E )) THEN
    WRITE(*,*) 'SFREE: TRIED TO DEALLOCATE THOUGH NOT ASSOCIATED.'
    STOP
    ENDIF

    IF (A%IS_ALLOCATED) DEALLOCATE( A%E , STAT=DEALLOC_STATUS )
    NULLIFY ( A%E )

    ! Reset all fileds to safe default values
    A%SPACE = -1
    A%INR = 0
    A%INTH = 0
    A%INX = 0
    A%LN = 0.D0
    A%IS_ALLOCATED = .FALSE.

    RETURN
END SUBROUTINE SFREE
!=======================================================================
SUBROUTINE COPY0(B,A)
!=======================================================================
! [USAGE]: 
! COPY THE SCALAR-TYPE VARIABLE A TO B
! WRAPPED BY 'EQUAL(=)' SYMBOL, THUS ABLE TO BE CALLED VIA 'B = A'
! [PARAMETERS]:
! A >> ORIGINAL VARIABLE 
! B >> SCALAR-TYPE VARIABLE TO BE COPIED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR),INTENT(INOUT):: B
    TYPE(SCALAR),INTENT(IN):: A
    INTEGER:: MM

    IF(.NOT.ASSOCIATED(B%E)) THEN
        WRITE(*,*) 'COPY0: NO MEMORY FOR LHS.'
        STOP
    ENDIF
    IF (.NOT. A%IS_ALLOCATED) THEN
        WRITE(*,*) 'COPY0: RHS IS NOT ALLOCATED.'
        STOP
    ENDIF

    ! IF(SIZE(B%E)<SIZE(A%E)) THEN
    ! CALL DEALLOCATE(B)
    ! CALL ALLOCATE(B,A%SPACE)
    ! ENDIF
    IF(B%SPACE.NE.A%SPACE) THEN
        CALL DEALLOCATE(B)
        CALL ALLOCATE(B,A%SPACE)
    ENDIF

    ! B%SPACE = A%SPACE
    B%LN = A%LN
    B%INR = A%INR
    B%INTH = A%INTH
    B%INX = A%INX

    ! IF(A%SPACE .EQ. FFF_SPACE) THEN
    !     DO MM=1,NTCHOP
    !         B%E(:NRCHOPS(MM),MM,:) = A%E(:NRCHOPS(MM),MM,:)
    !     ENDDO
    ! CALL CHOPDO(B)
    ! ELSE
    !     B%E(:SIZE(A%E,1),:SIZE(A%E,2),:SIZE(A%E,3)) = A%E
    !     B%E(:,:,SIZE(A%E,3)+1:)=0
    !     B%E(:,SIZE(A%E,2)+1:,:SIZE(A%E,3))=0
    !     B%E(SIZE(A%E,1)+1:,:SIZE(A%E,2),:SIZE(A%E,3))=0
    ! ENDIF
    IF(A%SPACE .EQ. FFF_SPACE) THEN
        B%E = A%E
    CALL CHOPDO(B)
    ELSE
        B%E(:SIZE(A%E,1),:SIZE(A%E,2),:SIZE(A%E,3)) = A%E
        B%E(:,:,SIZE(A%E,3)+1:)=0
        B%E(:,SIZE(A%E,2)+1:,:SIZE(A%E,3))=0
        B%E(SIZE(A%E,1)+1:,:SIZE(A%E,2),:SIZE(A%E,3))=0
    ENDIF

    RETURN
END SUBROUTINE COPY0
!=======================================================================
SUBROUTINE HORFFT(A,IS)
!=======================================================================
! [USAGE]: 
! PERFORM HORIZONTAL FFT WITH RESPECT TO AZIMUTHAL (THETA) DIRECTION
! EXPECTED TO BE CALLED WHEN DATA IS IN PPP SPACE
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (R,THETA,Z) VALUES
! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL RTRAN(A,1) (FFF -> PFF)
!    CALL VERFFT(A,1) (PFF -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,-1) (PPP -> PFP)
!    CALL VERFFT(A,-1) (PFP -> PFF)
!    CALL RTRAN(A,-1) (PFF -> FFF)
! [DEPENDENCIES]
! 1. DFFTW_PLAN_DFT                           @ FFTW3 
! 2. DFFTW_EXECUTE_DFT                        @ FFTW3
! 3. DFFTW_DESTROY_PLAN                       @ FFTW3
! 4. FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE @ FFTW3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    TYPE(SCALAR):: A
    INTEGER,INTENT(IN):: IS

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: B
    REAL(P8),DIMENSION(:),ALLOCATABLE:: C
    INTEGER:: II,KK,MM,RSIZE,XSIZE
    REAL(P8):: FAC
    INTEGER(P8):: PLAN

    ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
    ! PFP: NDIMR/N1, NDIMTH, NDIMX/N2

    RSIZE = SIZE(A%E,1)
    XSIZE = SIZE(A%E,3)
    CALL CHOPDO(A) 

    ! PHYSICAL TO FOURIER SPACE:
    IF(IS .EQ.-1) THEN
        IF(A%SPACE.NE.PPP_SPACE)THEN
            IF (MPI_RANK.EQ.0) WRITE(*,*) 'HORFFT: NOT IN PHYSICAL SPACE'
            STOP
        ENDIF
        
        IF(NTH.EQ.1) THEN
            A%SPACE=PFP_SPACE
            RETURN
        ENDIF

        ! CREATE PLAN OUTSIDE OF PARALLEL REGION (THREAD-SAFETY)
        ALLOCATE(B(NDIMTH))
        B = CMPLX(0.D0)
        ALLOCATE(C(2*NTH))
        C = 0.D0
        CALL DFFTW_PLAN_DFT_R2C_1D(PLAN,2*NTH,C,B,FFTW_ESTIMATE)
        DEALLOCATE(B)
        DEALLOCATE(C)

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(B,C,KK,II,MM)
        ALLOCATE(B(NDIMTH))
        ALLOCATE(C(2*NTH))
        
        !$OMP DO COLLAPSE(2)
        DO KK=1,XSIZE !NX
            DO II=1,RSIZE !NR
                B = A%E(II,:NDIMTH,KK)
                DO MM = 1,NTH
                    C(2*MM-1) = REAL(B(MM))
                    C(2*MM  ) = AIMAG(B(MM))
                ENDDO
                
                CALL DFFTW_EXECUTE_DFT_R2C(PLAN,C,B)
                A%E(II,:NDIMTH,KK) = B/(2*NTH)
            ENDDO
        ENDDO
        !$OMP END DO

        DEALLOCATE(B)
        DEALLOCATE(C)
        !$OMP END PARALLEL

        CALL DFFTW_DESTROY_PLAN(PLAN)

        ! CHOPPING IN THETA DIRECTION:
        A%E(:,NTCHOP+1:,:) = CMPLX(0.D0,0.D0)

        ! ELIMINATE COMPLEX RESIDUAL DUE TO MACHINE ERROR
        A%E(:,1,:) = CMPLX(REAL(A%E(:,1,:)),0.D0)

        A%SPACE = PFP_SPACE

    ! FOURIER TO PHYSICAL SPACE:
    ELSE
        IF(A%SPACE.NE.PFP_SPACE)THEN
            IF (MPI_RANK.EQ.0) WRITE(*,*) 'HORFFT: NOT IN PFP SPACE'
            STOP
        ENDIF
        
        IF(NTH.EQ.1) THEN
            A%SPACE=PPP_SPACE
            RETURN
        ENDIF

        ! ELIMINATE COMPLEX RESIDUAL DUE TO MACHINE ERROR
        A%E(:,1,:) = CMPLX(REAL(A%E(:,1,:)),0.D0)

        ! CHOPPING IN THETA DIRECTION:
        A%E(:,NTCHOP+1:,:) = CMPLX(0.D0,0.D0)

        ! CREATE PLAN OUTSIDE OF PARALLEL REGION (THREAD-SAFETY)
        ALLOCATE(B(NDIMTH))
        B = CMPLX(0.D0)
        ALLOCATE(C(2*NTH))
        C = 0.D0
        CALL DFFTW_PLAN_DFT_C2R_1D(PLAN,2*NTH,B,C,FFTW_ESTIMATE)
        DEALLOCATE(B)
        DEALLOCATE(C)

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(B,C,KK,II,MM)
        ALLOCATE(B(NDIMTH))
        ALLOCATE(C(2*NTH))
        
        !$OMP DO COLLAPSE(2)
        DO KK=1,XSIZE !NX
            DO II=1,RSIZE !NR
                B = A%E(II,:NDIMTH,KK)
                C = 0.D0
                CALL DFFTW_EXECUTE_DFT_C2R(PLAN,B,C)
                
                DO MM = 1,NTH
                    A%E(II,MM,KK) = CMPLX(C(2*MM-1),C(2*MM),P8)
                ENDDO
                A%E(II,NDIMTH,KK) = CMPLX(0.D0, 0.D0)
            ENDDO
        ENDDO
        !$OMP END DO

        DEALLOCATE(B)
        DEALLOCATE(C)
        !$OMP END PARALLEL

        CALL DFFTW_DESTROY_PLAN(PLAN)
        A%SPACE = PPP_SPACE
    ENDIF

    CALL CHOPDO(A)
    RETURN
END SUBROUTINE HORFFT
!=======================================================================
SUBROUTINE VERFFT(A,IS)
!=======================================================================
! [USAGE]: 
! PERFORM VERTICAL FFT WITH RESPECT TO AXIAL (Z) DIRECTION
! EXPECTED TO BE CALLED WHEN DATA IS IN PFP SPACE
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (R,THETA,Z) VALUES
! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL RTRAN(A,1) (FFF -> PFF)
!    CALL VERFFT(A,1) (PFF -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,-1) (PPP -> PFP)
!    CALL VERFFT(A,-1) (PFP -> PFF)
!    CALL RTRAN(A,-1) (PFF -> FFF)
! [DEPENDENCIES]
! 1. DFFTW_PLAN_DFT                           @ FFTW3 
! 2. DFFTW_EXECUTE_DFT                        @ FFTW3
! 3. DFFTW_DESTROY_PLAN                       @ FFTW3
! 4. FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE @ FFTW3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    TYPE(SCALAR),INTENT(INOUT):: A
    INTEGER,INTENT(IN):: IS

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: B
    INTEGER:: JJ, KK
    REAL(P8):: FAC
    INTEGER(P8):: PLAN

    COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: A_PFP0, A_PFP1, A_PFF0, A_PFF1
    REAL(P8):: LN

    ! PFP: NDIMR/N1, NDIMTH, NDIMX/N2
    ! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1

    ! MPI PREP =============================================================
    LN = A%LN
    CALL CHOPDO(A)
    CALL MPI_BARRIER(MPI_COMM_IVP, IERR)

    ! ALLOCATE: TWO UTILITY ARRAYS IN PFP0 AND PFP1
    ALLOCATE(A_PFP0(SIZE_PFP0(1), SIZE_PFP0(2), SIZE_PFP0(3)))
    ALLOCATE(A_PFP1(SIZE_PFP1(1), SIZE_PFP1(2), SIZE_PFP1(3)))
    ! ALLOCATE: TWO UTILITY ARRAYS IN PFF0 AND PFF1
    ALLOCATE(A_PFF0(SIZE_PFF0(1), SIZE_PFF0(2), SIZE_PFF0(3)))
    ALLOCATE(A_PFF1(SIZE_PFF1(1), SIZE_PFF1(2), SIZE_PFF1(3)))  

    ! PHYSICAL TO FOURIER SPACE: ===========================================
    IF(IS .EQ. -1) THEN
        IF(A%SPACE .NE. PFP_SPACE) THEN
            IF (MPI_RANK .EQ. 0) WRITE(*,*) 'VERFFT: NOT IN PFP SPACE'
            STOP
        ENDIF
        
        IF(NX .EQ. 1) THEN
            A%SPACE = PFF_SPACE
            RETURN
        ENDIF
     
        ! CHOP: NDIMTH -> NTCHOPDIM (NDIMR/N1, NTCHOPDOM, NDIMX/N2)
        A_PFP0 = A%E(:, :NTCHOPDIM, :)

        ! EXCHANGE: A_PFP0 -> A_PFP1(NDIMR/N1, NTCHOPDIM/N2, NDIMX)
        CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_2, A_PFP0, TYPE_PFP0, A_PFP1, TYPE_PFP1)

        DEALLOCATE(A_PFP0)
        
        ! Create FFT plan outside of the parallel region
        ALLOCATE(B(NX))
        B = CMPLX(0.D0)
        CALL DFFTW_PLAN_DFT_1D(PLAN, NX, B, B, FFTW_FORWARD, FFTW_ESTIMATE)
        DEALLOCATE(B)

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(B, JJ, KK) 
        ALLOCATE(B(NX))
        
        !$OMP DO COLLAPSE(2)
        DO JJ = 1, SIZE(A_PFP1, 2) !NTCHOP
            DO KK = 1, SIZE(A_PFP1, 1) !NR
                B = A_PFP1(KK, JJ, :NX)
                CALL DFFTW_EXECUTE_DFT(PLAN, B, B)
                A_PFP1(KK, JJ, :NX) = B / NX
            ENDDO
        ENDDO
        !$OMP END DO

        DEALLOCATE(B)
        !$OMP END PARALLEL

        CALL DFFTW_DESTROY_PLAN(PLAN)

        ! CHOP: NDIMX -> NXCHOPDIM (NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM)
        A_PFF0(:, :, :NXCHOP) = A_PFP1(:, :, :NXCHOP)
        A_PFF0(:, :, NXCHOP+1:) = A_PFP1(:, :, NXCHOPH:)

        ! EXCHANGE: A_PFF0 -> A_PFF1(NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1)
        CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_1, A_PFF0, TYPE_PFF0, A_PFF1, TYPE_PFF1)

        DEALLOCATE(A_PFP1, A_PFF0)

        ! OUTPUT A BACK (NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1)
        CALL DEALLOCATE(A)
        CALL ALLOCATE(A, PFF_SPACE) ! RESET A%INR, INTH, INX, SPACE
        A%E = A_PFF1
        A%LN = LN
        DEALLOCATE(A_PFF1)

    ! FOURIER TO PHYSICAL ==================================================
    ELSE                                                               
        IF(A%SPACE .NE. PFF_SPACE) THEN
            IF (MPI_RANK .EQ. 0) WRITE(*,*) 'VERFFT: NOT IN PFF SPACE'
            STOP
        ENDIF
        
        IF(NX .EQ. 1) THEN
            A%SPACE = PFP_SPACE
            RETURN
        ENDIF
        
        ! A%E -> A_PFF1(NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1)
        A_PFF1 = A%E

        ! EXCHANGE: A_PFF1 -> A_PFF0(NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM)
        CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_1, A_PFF1, TYPE_PFF1, A_PFF0, TYPE_PFF0)

        ! UNCHOP: NXCHOPDIM -> NDIMX(NDIMR/N1, NTCHOPDIM/N2, NDIMX)
        A_PFP1 = CMPLX(0.D0, 0.D0)
        A_PFP1(:, :, :NXCHOP) = A_PFF0(:, :, :NXCHOP)
        A_PFP1(:, :, NXCHOPH:) = A_PFF0(:, :, NXCHOP+1:)
        A_PFP1(:, :, NXCHOP+1:NXCHOPH-1) = CMPLX(0.D0, 0.D0)
        
        DEALLOCATE(A_PFF0, A_PFF1)

        ! Create FFT plan outside of the parallel region
        ALLOCATE(B(NX))
        B = CMPLX(0.D0)
        CALL DFFTW_PLAN_DFT_1D(PLAN, NX, B, B, FFTW_BACKWARD, FFTW_ESTIMATE)
        DEALLOCATE(B)

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(B, JJ, KK) 
        ALLOCATE(B(NX))
        
        !$OMP DO COLLAPSE(2)
        DO JJ = 1, SIZE(A_PFP1, 2) !NTCHOP
            DO KK = 1, SIZE(A_PFP1, 1) !NR
                B = A_PFP1(KK, JJ, :NX)
                CALL DFFTW_EXECUTE_DFT(PLAN, B, B)
                A_PFP1(KK, JJ, :NX) = B
            ENDDO
        ENDDO
        !$OMP END DO

        DEALLOCATE(B)
        !$OMP END PARALLEL

        CALL DFFTW_DESTROY_PLAN(PLAN)
        
        ! EXCHANGE: A_PFP1 -> A_PFP0(NDIMR/N1, NTCHOPDIM, NDIMX/N2)
        CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_2, A_PFP1, TYPE_PFP1, A_PFP0, TYPE_PFP0)

        DEALLOCATE(A_PFP1)
        
        ! OUTPUT A BACK (NDIMR/N1, NDIMTH, NDIMX/N2)
        CALL DEALLOCATE(A)
        CALL ALLOCATE(A, PFP_SPACE) ! RESET A%INR, INTH, INX
        ! UNCHOP: NTCHOPDIM -> NDIMTH (NDIMR/N1, NDIMTH, NDIMX/N2)
        A%E(:, :NTCHOPDIM, :) = A_PFP0
        A%E(:, NTCHOPDIM+1:, :) = CMPLX(0.D0, 0.D0)
        A%LN = LN
        DEALLOCATE(A_PFP0)
        A%SPACE = PFP_SPACE
    ENDIF

    CALL CHOPDO(A)
    CALL MPI_BARRIER(MPI_COMM_IVP, IERR)

    RETURN
END SUBROUTINE VERFFT
!=======================================================================
SUBROUTINE PLNFFT(A,IS)
!=======================================================================
! [USAGE]:
! PERFORM PLANAR FFT WITH RESPECT TO AZIMUTHAL (THETA) AND AXIAL (Z)
! DIRECTIONS
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (R,THETA,Z) VALUES
! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR),INTENT(INOUT):: A
    INTEGER,INTENT(IN):: IS

    ! PHYSICAL TO FOURIER SPACE
    IF (IS.EQ.-1) THEN
        IF (A%SPACE.NE.PPP_SPACE) THEN
            CALL MPRINT('PLNFFT: NOT IN PPP SPACE.')
            CALL MPI_ABORT(MPI_COMM_WORLD, ERR_FLAGS%TRANSFORM, IERR)
        ENDIF
        CALL HORFFT(A, -1) ! PPP -> PFP
        CALL VERFFT(A, -1) ! PFP -> PFF

    ! FOURIER TO PHYSICAL SPACE
    ELSE
        IF (A%SPACE.NE.PFF_SPACE) THEN
            CALL MPRINT('PLNFFT: NOT IN PFF SPACE.')
            CALL MPI_ABORT(MPI_COMM_WORLD, ERR_FLAGS%TRANSFORM, IERR)
        ENDIF
        CALL VERFFT(A, 1) ! PFF -> PFP
        CALL HORFFT(A, 1) ! PFP -> PPP
    ENDIF

    RETURN
END SUBROUTINE PLNFFT
!=======================================================================
SUBROUTINE RTRAN(A,IS)
!=======================================================================
! [USAGE]: 
! PERFORM MAPPED LEGENDRE TRANSFORM WITH RESPECT TO X (OR R) DIRECTION
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL RTRAN(A,1) (FFF -> PFF)
!    CALL VERFFT(A,1) (PFF -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,-1) (PPP -> PFP)
!    CALL VERFFT(A,-1) (PFP -> PFF)
!    CALL RTRAN(A,-1) (PFF -> FFF)
! [DEPENDENCIES]
! 1. ALLOCATE(~) @ MOD_SCALAR3
! 2. CHOPDO(~) @ MOD_SCALAR3
! 3. DEALLOCATE(~)  @ MOD_SCALAR3
! 4. OPERATOR(.MUL.) @ MOD_EIG
! [UPDATES]:
! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR),INTENT(INOUT):: A
    INTEGER,INTENT(IN):: IS

    TYPE(SCALAR):: B
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: BE,BO
    INTEGER:: NN,I,MM
    INTEGER:: XSIZE, THSIZE

    ! MPI PREP
    XSIZE = SIZE(A%E,3)
    THSIZE = SIZE(A%E,2)

    ! PHYSICAL TO MAPPED LEGENDRE SPACE:
    IF(IS.EQ.-1) THEN                                                  
    IF(A%SPACE .NE. PFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'RTRAN: NOT IN PFF_SPACE.'
        STOP
    ENDIF

    CALL ALLOCATE(B,FFF_SPACE)
    CALL CHOPDO(B)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(NN,BE,BO,MM,I)

    ! ALLOCATE( BE(NRH,NXCHOPDIM) )
    ! ALLOCATE( BO(NRH,NXCHOPDIM) )
    ! ALLOCATE( BE(NRH,XSIZE) )
    ! ALLOCATE( BO(NRH,XSIZE) )
    ALLOCATE( BE(XSIZE,NRH) )
    ALLOCATE( BO(XSIZE,NRH) )

    !$OMP DO
    DO MM=1, THSIZE !NTCHOP

        ! NN = NRCHOPS(MM+A%INTH)
        NN = NRCHOPS(MM+B%INTH)
        IF (NN.LT.1) CYCLE

        DO I=1,NRH
            ! BE(I,:NXCHOP) = (A%E(I,MM,:NXCHOP)+A%E(NR-I+1,MM,:NXCHOP))&
            !                 *TFM%W(I)
            ! BO(I,:NXCHOP) = (A%E(I,MM,:NXCHOP)-A%E(NR-I+1,MM,:NXCHOP))&
            !                 *TFM%W(I)

            ! IF(NXCHOP .NE. 1) THEN
            !     BE(I,NXCHOP+1:) = (A%E(I,MM,NXCHOPH:NX)+A%E(NR-I+1,MM,&
            !                         NXCHOPH:NX))*TFM%W(I)  
            !     BO(I,NXCHOP+1:) = (A%E(I,MM,NXCHOPH:NX)-A%E(NR-I+1,MM,&
            !                         NXCHOPH:NX))*TFM%W(I)
            ! ENDIF

            ! A%E: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1
            ! BE(I,:) = (A%E(I,MM,:)+A%E(NR-I+1,MM,:))*TFM%W(I)
            ! BO(I,:) = (A%E(I,MM,:)-A%E(NR-I+1,MM,:))*TFM%W(I)
            BE(:,I) = (A%E(I,MM,:)+A%E(NR-I+1,MM,:))*TFM%W(I)
            BO(:,I) = (A%E(I,MM,:)-A%E(NR-I+1,MM,:))*TFM%W(I)

        ENDDO

        ! TFM%PF(RadialCollocPts, RadialModes, AzimuthalModes)
        ! IF(NN.GE.1) THEN
        ! B%E(1:NN:2,MM,:)=TRANSPOSE(TFM%PF(:NRH,1:NN:2,MM+A%INTH)) .MUL. BE
        ! ENDIF
        B%E(1:NN:2,MM,:)=TRANSPOSE(BE .MUL. TFM%PF(:NRH,1:NN:2,MM+B%INTH))
        
        IF (NN.GE.2) THEN
        ! B%E(2:NN:2,MM,:)=TRANSPOSE(TFM%PF(:NRH,2:NN:2,MM+A%INTH)) .MUL. BO
        B%E(2:NN:2,MM,:)=TRANSPOSE(BO .MUL. TFM%PF(:NRH,2:NN:2,MM+B%INTH))
        ENDIF

    ENDDO
    !$OMP END DO

    DEALLOCATE( BE,BO )

!$OMP END PARALLEL

    CALL DEALLOCATE(A)

    A%E => B%E
    A%INR = B%INR
    A%INTH = B%INTH
    A%INX = B%INX
    NULLIFY(B%E)
    A%SPACE = FFF_SPACE
    A%IS_ALLOCATED = .TRUE.

    ! MAPPED LEGENDRE SPACE TO PHYSICAL:
    ELSE                                                               
    IF(A%SPACE .NE. FFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'RTRAN: NOT IN FFF_SPACE'
        STOP
    ENDIF

    CALL ALLOCATE(B,PFF_SPACE)
    CALL CHOPDO(B)
    B%LN=0.D0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(NN,BE,BO,MM)

    ALLOCATE( BE(NRH,XSIZE) ) !NXCHOPDIM) )
    ALLOCATE( BO(NRH,XSIZE) ) !NXCHOPDIM) )

    !$OMP DO
    DO MM=1,THSIZE !NTCHOP
        NN = NRCHOPS(MM+A%INTH)

        IF(NN.GE.1) THEN
        BE= TFM%PF(:NRH,1:NN:2,MM+A%INTH) .MUL. A%E(1:NN:2,MM,:)!NXCHOPDIM)
        ELSE
        BE=0
        ENDIF

        IF(NN.GE.2) THEN
        BO= TFM%PF(:NRH,2:NN:2,MM+A%INTH) .MUL. A%E(2:NN:2,MM,:)!NXCHOPDIM)
        ELSE
        BO=0
        ENDIF

        ! B%E(1:NRH,MM,1:NXCHOP)= BE(:,:NXCHOP)+BO(:,:NXCHOP)

        ! IF(NXCHOP .NE. 1) THEN
        ! B%E(1:NRH,MM,NXCHOPH:NX)= BE(:,NXCHOP+1:)+BO(:,NXCHOP+1:)
        ! B%E(NR:NRH+1:-1,MM,NXCHOPH:NX)= BE(:,NXCHOP+1:)-BO(:,NXCHOP+1:)
        ! ENDIF

        ! B%E(NR:NRH+1:-1,MM,:NXCHOP)= BE(:,:NXCHOP)-BO(:,:NXCHOP)

        B%E(1:NRH,MM,:)= BE(:,:)+BO(:,:)
        B%E(NR:NRH+1:-1,MM,:)= BE(:,:)-BO(:,:)

    ENDDO
    !$OMP END DO

    DEALLOCATE( BE,BO )

!$OMP END PARALLEL
    
    ! B%E(:,NTCHOP+1:,:)=0
    ! B%E(:,:NTCHOP,NXCHOP+1:NXCHOPH-1)=0
    B%E(:,MAX(1,NTCHOP-A%INTH+1):,:) = CMPLX(0.D0,0.D0)

    IF(A%INTH.EQ.0 .AND. A%INX.EQ.0) THEN
        IF(A%LN.NE.0.0) THEN
            WRITE(*,*) 'RTRAN: TO PHYSICAL SPACE THOUGH LOGTERM IS NONZERO.'
            WRITE(*,*) 'LOGTERM=',A%LN
            B%E(1:NR,1,1)=B%E(1:NR,1,1)+A%LN*TFM%LN(1:NR) ! ADD LOG TERM BACK TO THE FUNCTION IN PHYSICAL SPACE
        ENDIF
    ENDIF

    CALL DEALLOCATE(A)

    A%E => B%E
    NULLIFY(B%E)
    A%SPACE = PFF_SPACE
    A%INR = B%INR
    A%INTH = B%INTH
    A%INX = B%INX
    A%IS_ALLOCATED = .TRUE.

    ENDIF
    
    CALL CHOPDO(A)
    RETURN
END SUBROUTINE RTRAN
!=======================================================================
SUBROUTINE TOFF(A)
!=======================================================================
! [USAGE]: 
! TRANSFORM A FROM PPP SPACE TO FFF SPACE
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN PPP SPACE
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    INTEGER:: N

    N = A%SPACE
    IF (N.EQ.PPP_SPACE) THEN
    CALL HORFFT(A,-1)
    CALL VERFFT(A,-1)
    CALL RTRAN(A,-1)
    ELSEIF (N.EQ.PFP_SPACE) THEN
    CALL VERFFT(A,-1)
    CALL RTRAN(A,-1)
    ELSEIF (N.EQ.PFF_SPACE) THEN
    CALL RTRAN(A,-1)
    ENDIF

    RETURN
END SUBROUTINE TOFF
!=======================================================================
SUBROUTINE TOFP(A)
!=======================================================================
! [USAGE]: 
! TRANSFORM A FROM FFF SPACE TO PPP SPACE
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN FFF SPACE
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    INTEGER:: N

    N = A%SPACE
    IF (N.EQ.FFF_SPACE) THEN
    CALL RTRAN(A,1)
    CALL VERFFT(A,1)
    CALL HORFFT(A,1)
    ELSEIF (N.EQ.PFF_SPACE) THEN
    CALL VERFFT(A,1)
    CALL HORFFT(A,1)
    ELSEIF (N.EQ.PFP_SPACE) THEN
    CALL HORFFT(A,1)
    ENDIF
    
    RETURN
END SUBROUTINE TOFP
! ======================================================================
SUBROUTINE DECOMPOSE(NSIZE,NPROCS,PROC_NUM,NSIZE_PROC,INDEX_PROC)
! ======================================================================
! [USAGE]:
! DECOMPOSE 1D DOMAIN INTO ALL PROCESSORS
! [PARAMETERS]:
! NSIZE >> TOT NUM OF ELEMENTS ALONG THAT DIMENSION 
! NPROCS >> TOT NUM OF PROCS ALONG THAT DIMENSION
! PROC_NUM >> RANK OF THAT PROCESSOR
! NSIZE_PROC >> LOC NUM OF ELEMENTS ALONG THAT DIMENSION
! INDEX_PROC >> START INDEX OF THAT PROCESSOR ALONG THAT DIMENSION
! [NOTE]:
! INDEX STARTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NSIZE, NPROCS, PROC_NUM, NSIZE_PROC, INDEX_PROC
    INTEGER:: Q,R

    Q = NSIZE/NPROCS
    R = MOD(NSIZE,NPROCS)
    IF (R.GT.PROC_NUM) THEN
        NSIZE_PROC = Q+1
        INDEX_PROC = NSIZE_PROC*PROC_NUM
    ELSE
        NSIZE_PROC = Q
        INDEX_PROC = NSIZE_PROC*PROC_NUM + R
    ENDIF

    ! NSIZE_PROC = Q + (R.GT.PROC_NUM) ! DOESN'T WORK
    ! INDEX_PROC = Q * PROC_NUM + MIN(R,PROC_NUM)

    RETURN
end SUBROUTINE DECOMPOSE
! ======================================================================
SUBROUTINE SUBARRAY(ELEMENT_DATA_TYPE,NDIM,DATASIZE_PROC,DIM,NPROCS,NEW_DATA_TYPE)
! ======================================================================
! [USAGE]: CREATE SUBARRAY DATATYPES FOR EACH PROC. ORIGINAL ARRAY IS DE
! COMPOSED ALONG DIM INTO SUBARRAYS.
! [PARAMETERS]: 
! OLD_DATA_TYPE >> ORIGINAL MPI ELEMENTARY DATATYPE OF THE ARRAY
! NDIM >> NUMBER OF DIMENSIONS OF THE LOCAL ARRAY
! DATASIZE_PROC >> ORIGINAL DATASIZE OF THE LOCAL ARRAY
! DIM >> DIMENSION ALONG WHICH THE LOCAL ARRAY IS FURTHER DECOMPOSED
! NPROCS >> NUMBER OF PROCS ALONG THAT DIM
! NEW_DATA_TYPE >> MPI DERIVED DATATYPE FOR ALL PROCS (DIM = NPROCS)
! [NOTE]:
! 1. PROC INDEX STARTS FROM 0
! 2. DIMN INDEX STARTS FROM 1
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: ELEMENT_DATA_TYPE, NDIM, DIM, NPROCS
    INTEGER,DIMENSION(NDIM):: DATASIZE_PROC, SUBDSIZE_PROC, SUBSTARTS_PROC
    INTEGER,DIMENSION(0:NPROCS-1):: NEW_DATA_TYPE

    INTEGER:: I 

    DO I = 1,NDIM ! MEMO >>> MAY NEED TO CHANGE THIS TO 0~NDIM-1
        SUBDSIZE_PROC(I) = DATASIZE_PROC(I)
        SUBSTARTS_PROC(I) = 0
    ENDDO

    ! ALONG THAT DIM
    DO I = 0,NPROCS-1
        ! DECOMPOSE A DOMAIN OF DATASIZE_PROC(DIM) INTO NPROCS
        ! -> DECOMPOSED DOMAIN SIZE ALONG THAT DIM: NSIZE_PROC(DIM)
        ! -> OTHER DIMS KEEP THE ORIGINAL DOMAIN SIZE
        ! -> DECOMPOSED DOMAIN START INDEX ALONG THAT DIM: SUBSTARTS_PROC(DIM)
        ! -> OTHER DIMS KEEP THE ORIGINAL START INDEX
        CALL DECOMPOSE(DATASIZE_PROC(DIM),NPROCS,I,SUBDSIZE_PROC(DIM),SUBSTARTS_PROC(DIM))
        ! CREATE A NDIM-DIMENSION SUBARRAY DATATYPE FOR EACH PROC
        CALL MPI_TYPE_CREATE_SUBARRAY(NDIM,DATASIZE_PROC,SUBDSIZE_PROC,SUBSTARTS_PROC, & 
            MPI_ORDER_FORTRAN,ELEMENT_DATA_TYPE,NEW_DATA_TYPE(I),IERR)
        CALL MPI_TYPE_COMMIT(NEW_DATA_TYPE(I),IERR)
    ENDDO

    RETURN
END SUBROUTINE SUBARRAY
! ======================================================================
SUBROUTINE EXCHANGE_3DCOMPLEX(COMM,DATA_SIZE_PROC_OLD,ARRAY_PROC_OLD,DIM_OLD &
                                    ,DATA_SIZE_PROC_NEW,ARRAY_PROC_NEW,DIM_NEW)
! ======================================================================
! [USAGE]:
! FOR PROCESSOR GROUP 'COMM', SWAP DIM_OLD AND DIM_NEW.
! [PARAMETERS]:
! COMM >> PROCESSOR GROUP
! DATA_SIZE_PROC >> SIZE OF THE LOCAL ARRAY
! ARRAY_PROC >> LOCAL ARRAY TO BE SWAPPED
! DIM >> DIMENSIONS TO BE SWAPPED
! [NOTE]
! 1. USE THIS ONLY IF THE DATATYPE ONLY NEEDS TO BE USED ONCE
! 2. PROC INDEX STARTS FROM 0
! 3. DIMN INDEX STARTS FROM 1
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NDIM = 3
    INTEGER:: COMM, DIM_OLD, DIM_NEW
    INTEGER,DIMENSION(3),INTENT(IN):: DATA_SIZE_PROC_OLD, DATA_SIZE_PROC_NEW
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: ARRAY_PROC_OLD
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(INOUT):: ARRAY_PROC_NEW

    INTEGER:: I , NPROC_COMM
    INTEGER,DIMENSION(:),ALLOCATABLE:: DATA_TYPE_OLD, DATA_TYPE_NEW, counts, displs

    ! DETERMINE THE NUMBER OF PROCS IN PROCESSOR GROUP: COMM
    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    ALLOCATE(DATA_TYPE_OLD(0:NPROC_COMM-1))
    ALLOCATE(DATA_TYPE_NEW(0:NPROC_COMM-1))
    ALLOCATE(counts(0:NPROC_COMM-1))
    ALLOCATE(displs(0:NPROC_COMM-1))

    ! CREATE SUBARRAY DATATYPES
    ! OLD LOCAL ARRAY: ARRAY_OLD OF SIZE_OLD IS DECOMPOSED INTO SUBARRAY ALONG DIM_OLD
    CALL SUBARRAY(MPI_DOUBLE_COMPLEX,NDIM,DATA_SIZE_PROC_OLD,DIM_OLD,NPROC_COMM,DATA_TYPE_OLD)
    ! NEW LOCAL ARRAY: ARRAY_NEW OF SIZE_NEW IS DECOMPOSED INTO SUBARRAY ALONG DIM_NEW
    CALL SUBARRAY(MPI_DOUBLE_COMPLEX,NDIM,DATA_SIZE_PROC_NEW,DIM_NEW,NPROC_COMM,DATA_TYPE_NEW)

    ! SWAP SUBARRAYS
    ! EACH SUBARRAY IS COUNTED AS ONE UNIT
    DO I = 0,NPROC_COMM-1
        counts(I) = 1; ! SWAP ONE SUBARRAY A TIME
        displs(I) = 0; ! DIRECTLY START FROM THE FIRST MEMORY LOC OF THE ARRAY
    enddo

    CALL MPI_ALLTOALLW(ARRAY_PROC_OLD, counts, displs, DATA_TYPE_OLD, &
                    ARRAY_PROC_NEW, counts, displs, DATA_TYPE_NEW, COMM, IERR)

    ! FREE THE DATATYPE
    DO I = 0,NPROC_COMM-1
        CALL MPI_TYPE_FREE(DATA_TYPE_OLD(I),IERR)
        CALL MPI_TYPE_FREE(DATA_TYPE_NEW(I),IERR)
    ENDDO

    DEALLOCATE(DATA_TYPE_OLD, DATA_TYPE_NEW, counts, displs)

end SUBROUTINE EXCHANGE_3DCOMPLEX
! ======================================================================
SUBROUTINE EXCHANGE_3DCOMPLEX_FAST(COMM,ARRAY_PROC_OLD,DATA_TYPE_OLD &
                                        ,ARRAY_PROC_NEW,DATA_TYPE_NEW)
! ======================================================================
! [USAGE]:
! FOR PROCESSOR GROUP 'COMM', SWAP DIM_OLD AND DIM_NEW.
! [PARAMETERS]:
! COMM >> PROCESSOR GROUP
! DATA_SIZE_PROC >> SIZE OF THE LOCAL ARRAY
! ARRAY_PROC >> LOCAL ARRAY TO BE SWAPPED
! DIM >> DIMENSIONS TO BE SWAPPED
! [NOTE]
! 1. USE THIS ONLY IF THE DATATYPE ONLY NEEDS TO BE USED ONCE
! 2. PROC INDEX STARTS FROM 0
! 3. DIMN INDEX STARTS FROM 1
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NDIM = 3
    INTEGER:: COMM
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: ARRAY_PROC_OLD
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(INOUT):: ARRAY_PROC_NEW
    INTEGER,DIMENSION(:),INTENT(IN):: DATA_TYPE_OLD, DATA_TYPE_NEW

    INTEGER:: I , NPROC_COMM
    INTEGER,DIMENSION(:),ALLOCATABLE:: counts, displs, DATA_TYPE_OLD_COPY

    ! CHECK
    IF (RANK(ARRAY_PROC_OLD).NE.NDIM) THEN
    ! ARRAY_PROC_OLD must have NDIM ranks
    ! i.e NDIM = 3, ARRAY_PROC_OLD(:,:,:)
        CALL MPI_COMM_RANK(MPI_COMM_IVP, MPI_RANK, IERR)
        WRITE(*,*) 'ERROR in GLOBAL PROC#'//ITOA3(MPI_RANK)//': EXCHANGE - INCOMPATIBLE INPUT ARRAY DIMS'
    ENDIF
    ! IF (LBOUND(DATA_TYPE_OLD,1).NE.0) THEN
    ! ! DATA_TYPE_OLD(NPROC_COMM) should count from 0
    ! ! i.e DATA_TYPE_OLD(0), ... , DATA_TYPE_OLD(NPROC_COMM-1)
    !     CALL MPI_COMM_RANK(MPI_COMM_IVP,PROC_NUM,IERR)
    !     IF (PROC_NUM.EQ.0) THEN
    !         WRITE(*,*) 'ERROR in GLOBAL PROC#'//ITOA3(PROC_NUM)//': EXCHANGE - SUBARRAY DATATYPE MUST START FROM 0'
    !     ENDIF
    ! ENDIF

    ! SWAP SUBARRAYS
    ! EACH SUBARRAY IS COUNTED AS ONE UNIT
    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    ALLOCATE(counts(0:NPROC_COMM-1))
    ALLOCATE(displs(0:NPROC_COMM-1))
    ALLOCATE(DATA_TYPE_OLD_COPY(0:NPROC_COMM-1))
    DATA_TYPE_OLD_COPY(0:NPROC_COMM-1) = DATA_TYPE_OLD(:)

    DO I = 0,NPROC_COMM-1
    counts(I) = 1; ! SWAP ONE SUBARRAY A TIME
    displs(I) = 0; ! DIRECTLY START FROM THE FIRST MEMORY LOC OF THE ARRAY
    enddo

    CALL MPI_ALLTOALLW(ARRAY_PROC_OLD, counts, displs, DATA_TYPE_OLD, &
    ARRAY_PROC_NEW, counts, displs, DATA_TYPE_NEW, COMM, IERR)

    DEALLOCATE(counts, displs)

end SUBROUTINE EXCHANGE_3DCOMPLEX_FAST
! ======================================================================
SUBROUTINE SUBCOMM_CART(COMM,NDIM,SUBCOMMS)
! ======================================================================
! [USAGE]:
! CREATE COMMUNICATORS (SUBCOMMS(I)) FOR EACH DIMENSION(I) THAT HAS CART
! ESIAN TOPOLOGY. 
! [EXAMPLE]: 
! FOR N1xN2xN3 = NPROC_COMM (NDIM = 3) PROCESSORS,
! SUBCOMMS(1) >> N2XN3 GROUPS EACH HAS N1 PROCS
! SUBCOMMS(2) >> N1XN3 GROUPS EACH HAS N2 PROCS
! [NOTES]:
! EA PROC CALCULATES ITS CORRESPONDING SUBCOMMS. PROCS ON THE SAME LINE
! IN I-TH DIM WILL SHARE THE SAME SUBCOMMS(I)
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM
    INTEGER:: NDIM
    INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT)::SUBCOMMS

    INTEGER:: COMM_CART, NPROC_COMM, I 
    INTEGER,DIMENSION(1:NDIM):: DIMS
    LOGICAL,DIMENSION(1:NDIM):: PERIODS, REMDIMS

    DIMS = 0
    PERIODS = .FALSE.
    REMDIMS = .FALSE.
    ALLOCATE(SUBCOMMS(NDIM))

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)

    ! CREATES A DIVISION OF PROCESSORS IN A CARTESIAN GRID
    ! DIMS: NUMBER OF PROCESSORS IN EACH DIM
    CALL MPI_DIMS_CREATE(NPROC_COMM,NDIM,DIMS,IERR)
    ! write(*,*) MPI_RANK,'-',dims(1),'x',dims(2)

    ! MAKES PROCESSOR GROUP THAT ATTACHES THE CARTESIAN TOPOLOGY
    ! COMM: ORIGINAL PROCESSOR GROUP
    ! NDIM: NUMBER OF DIMS IN CARTESIAN GRID
    ! DIMS: NUMBER OF PROCESSORS IN EACH DIM
    ! PERI: WHETHER EA DIM IS PERIODIC
    ! REOR: REORDER OR NOT
    ! COMM_CART: NEW PROCESSOR GROUP
    CALL MPI_CART_CREATE(COMM,NDIM,DIMS,PERIODS,.TRUE.,COMM_CART,IERR)

    ! CREATES SUB-PROCESSOR GROUP
    ! REMDIMS: TELLS WHICH DIM(S) IS KEPT IN THE SUBGRID
    ! E.X: FOR 3D, REMDIMS = .FALSE.,.TRUE.,.FALSE.
    !      THEN CREATES DIM1xDIM3 SUBGROUPS, EA WITH DIM2 PROCESSORS
    DO I = 1,NDIM
        REMDIMS(I) = .TRUE.
        CALL MPI_CART_SUB(COMM_CART,REMDIMS,SUBCOMMS(I),IERR)
        REMDIMS(I) = .FALSE.
    ENDDO

    ! FOR 2D GRID ONLY:
    ! MAKE SURE SUBCOMMS(1) IS ALWAYS BIGGER THAN SUBCOMMS(2)
    IF (NDIM.EQ.2) THEN
        IF (count_proc(SUBCOMMS(1)).LT.count_proc(SUBCOMMS(2))) THEN
            IERR = SUBCOMMS(1)
            SUBCOMMS(1) = SUBCOMMS(2)
            SUBCOMMS(2) = IERR
        ENDIF
    ENDIF
    ! write(*,*) MPI_RANK,'-',count_proc(SUBCOMMS(1)),'x',count_proc(SUBCOMMS(2)), &
    !     '-',local_proc(SUBCOMMS(1)), &
    !     '-',local_proc(SUBCOMMS(2))

    ! FUTURE EDIT NOTE:
    ! FOR NTH = 1 OR NX = 1, MAKE SUBCOMMS(I) HAVE ONE PROC IN EA GROUP
    ! FOR SLAB DECOMP, MAKE SUBCOMMS(I) HAVE ONE PROC IN EA GROUP
    ! (USE MPI_COMM_SPLIT)

    MPI_COMM_IVP = COMM_CART
    CALL MPI_COMM_RANK(MPI_COMM_IVP,MPI_RANK,IERR)

    ! CALL MPI_COMM_FREE(COMM_CART,IERR)      

end SUBROUTINE SUBCOMM_CART
! ======================================================================
SUBROUTINE MASSEMBLE(LOCAL_ARRAY, GLOBAL_ARRAY, axis)
! ======================================================================
! [USAGE]: 
! ASSEMBLE COMPLEX LOCAL ARRAYS INTO GLOBAL ARRAY IN PROC#0
! [PARAMETERS]:
! axis >> axis along which the domain is NOT chopped    
! MPI_Datatype >> element type of local array (e.g. MPI_INTEGER)
! [UPDATES]:
! CODED BY JINGE WANG @ SEP 19 2021
! [NOTE1]:
! SUBCOMMS_L,axis,SUBCOMMS_R
! or axis,SUBCOMMS_L,SUBCOMMS_R
! DOES NOT work with SUBCOMMS_L,SUBCOMMS_R,axis (not needed in IVP)
! [NOTE2]:
! USE public SUBCOMM_1 and SUBCOMM_2. Need to change their definition to
! match the actual usage.
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: LOCAL_ARRAY
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(INOUT):: GLOBAL_ARRAY
    integer,INTENT(IN):: axis 

    INTEGER:: II, JJ, KK, N1_glb, N2_glb, N3_glb, N1_loc, N2_loc, N3_loc, ELEMENT_SIZE 
    INTEGER(KIND=MPI_ADDRESS_KIND):: EXTEND_SIZE
    INTEGER:: SUBCOMM_L, SUBCOMM_R, MPI_Datatype
    INTEGER:: SUBARRAY_TYPE, SUBARRAY_TYPE_resized, DISPLACEMENT_loc, RECVCOUNT_loc
    integer,DIMENSION(:),ALLOCATABLE:: DISPLACEMENT, RECVCOUNT
    COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: LOCAL_ARRAY2, SUBGLOBAL_ARRAY, GLOBAL_ARRAY_COPY

    N1_glb = SIZE(GLOBAL_ARRAY,1)
    N2_glb = SIZE(GLOBAL_ARRAY,2)
    N3_glb = SIZE(GLOBAL_ARRAY,3)

    N1_loc = SIZE(LOCAL_ARRAY,1)
    N2_loc = SIZE(LOCAL_ARRAY,2)
    N3_loc = SIZE(LOCAL_ARRAY,3)

    MPI_Datatype = MPI_DOUBLE_COMPLEX
    IF (CP8_SIZE.EQ.0) THEN
    CALL MPI_TYPE_SIZE(MPI_Double_Complex,CP8_SIZE,IERR)
    ENDIF
    ELEMENT_SIZE = CP8_SIZE

    
    IF (axis.EQ.1) THEN ! axis,SUBCOMMS_L,SUBCOMMS_R => FFF,PFF
    ! PFF: NDIMR    , NTCHOPDIM/N2, NXCHOPDIM/N1 << PFF
    ! FFF: NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1 << FFF

        SUBCOMM_L = SUBCOMM_2
        SUBCOMM_R = SUBCOMM_1

        ALLOCATE(SUBGLOBAL_ARRAY(N1_glb,N2_glb,N3_loc))
        ALLOCATE(LOCAL_ARRAY2(N1_glb,N3_loc,N2_loc))

        LOCAL_ARRAY2 = RESHAPE(LOCAL_ARRAY,SHAPE(LOCAL_ARRAY2),ORDER = [1,3,2])

        CALL MPI_TYPE_VECTOR(N3_loc, N1_glb, N2_glb*N1_glb, &
                            MPI_Datatype, SUBARRAY_TYPE, IERR)
        EXTEND_SIZE = ELEMENT_SIZE*N1_glb
        RECVCOUNT_loc = N2_loc
        DISPLACEMENT_loc = local_index(N2_glb,SUBCOMM_L)

    ELSEIF (axis.EQ.2) THEN ! SUBCOMMS_L,axis,SUBCOMMS_R => PPP,PFP
    ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2 << PPP
    ! PFP: NDIMR/N1, NDIMTH, NDIMX/N2 << PFP

        SUBCOMM_L = SUBCOMM_1
        SUBCOMM_R = SUBCOMM_2

        ALLOCATE(SUBGLOBAL_ARRAY(N2_glb,N1_glb,N3_loc))
        ALLOCATE(LOCAL_ARRAY2(N2_glb,N3_loc,N1_loc))

        DO KK = 1,N1_loc
            DO JJ = 1,N3_loc
                DO II = 1,N2_glb
                    LOCAL_ARRAY2(II,JJ,KK) = LOCAL_ARRAY(KK,II,JJ)
                ENDDO
            ENDDO
        ENDDO

        CALL MPI_TYPE_VECTOR(N3_loc, N2_glb, N2_glb*N1_glb, &
                            MPI_Datatype, SUBARRAY_TYPE, IERR)
        EXTEND_SIZE = ELEMENT_SIZE*N2_glb
        RECVCOUNT_loc = N1_loc
        DISPLACEMENT_loc = local_index(N1_glb,SUBCOMM_L)

    ELSE
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'MASSEMBLE: undefined for the current axis'
        RETURN
    ENDIF

    ! CREATE ELEMENTAL SUBARRAY OF LOCAL_ARRAY
    CALL MPI_TYPE_CREATE_RESIZED(SUBARRAY_TYPE, INT(0, MPI_ADDRESS_KIND), EXTEND_SIZE, SUBARRAY_TYPE_resized, IERR)
    CALL MPI_TYPE_COMMIT(SUBARRAY_TYPE_resized,IERR)

    ! IN EA SUBCOMM_L, GATHER LOCAL_ARRAY TO FORM SUBGLOBAL_ARRAY in SUBCOMM_L's #0 PROC
    ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMM_L)-1))
    ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMM_L)-1))
    CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMM_L,IERR)
    CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMM_L,IERR)
    CALL MPI_GATHERV(LOCAL_ARRAY2, N1_loc*N2_loc*N3_loc, MPI_Datatype, &
                    SUBGLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, SUBARRAY_TYPE_resized, &
                    0, SUBCOMM_L, IERR)

    ! IN THE SUBCOMM_R THAT CONTAINS SUBCOMM_Ls' #0 PROCs: 
    ! GATHER SUBGLOBAL_ARRAY TO FORM GLOBAL_ARRAY in SUBCOMM_R's #0 proc (=> GLOBAL #0 PROC)
    IF (local_proc(SUBCOMM_L).EQ.0) THEN

        DEALLOCATE(RECVCOUNT,DISPLACEMENT)
        ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMM_R)-1))
        ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMM_R)-1))
        RECVCOUNT_loc = N1_glb*N2_glb*N3_loc
        DISPLACEMENT_loc = N1_glb*N2_glb*local_index(N3_glb,SUBCOMM_R)
        CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMM_R,IERR)
        CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMM_R,IERR)

        IF (axis.EQ.1) THEN
        
            CALL MPI_GATHERV(SUBGLOBAL_ARRAY, N1_glb*N2_glb*N3_loc, MPI_Datatype, &
                            GLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, MPI_Datatype, &
                            0, SUBCOMM_R, IERR)

        ELSEIF (axis.EQ.2) THEN

            ALLOCATE(GLOBAL_ARRAY_COPY(N2_glb,N1_glb,N3_glb))
            CALL MPI_GATHERV(SUBGLOBAL_ARRAY, N1_glb*N2_glb*N3_loc, MPI_Datatype, &
                            GLOBAL_ARRAY_COPY, RECVCOUNT, DISPLACEMENT, MPI_Datatype, &
                            0, SUBCOMM_R, IERR)
            GLOBAL_ARRAY = RESHAPE(GLOBAL_ARRAY_COPY,SHAPE(GLOBAL_ARRAY),ORDER = [2,1,3]) ! NEED TO REORDER
            DEALLOCATE(GLOBAL_ARRAY_COPY)

        ENDIF    
    ENDIF                 

    CALL MPI_TYPE_FREE(SUBARRAY_TYPE_resized,IERR)
    DEALLOCATE(LOCAL_ARRAY2,SUBGLOBAL_ARRAY,RECVCOUNT,DISPLACEMENT)

END SUBROUTINE MASSEMBLE
! ======================================================================
SUBROUTINE MDISASSEMBLE(GLOBAL_ARRAY, LOCAL_ARRAY, axis)
! ======================================================================
! [USAGE]:
! DISASSEMBLE COMPLEX GLOBAL ARRAY IN PROC#0 INTO LOCAL ARRAYS
! [PARAMETERS]:
! axis >> axis along which the domain is NOT chopped
! MPI_Datatype >> element type of local array (e.g. MPI_INTEGER)
! [UPDATES]:
! CODED BY JINGE WANG @ SEP 19 2021
! [NOTE1]:
! SUBCOMMS_L,axis,SUBCOMMS_R
! or axis,SUBCOMMS_L,SUBCOMMS_R
! DOES NOT work with SUBCOMMS_L,SUBCOMMS_R,axis (not needed in IVP)
! [NOTE2]:
! USE public SUBCOMM_1 and SUBCOMM_2. Need to change their definition to
! match the actual usage.
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(INOUT):: LOCAL_ARRAY
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: GLOBAL_ARRAY
    integer,INTENT(IN):: axis

    INTEGER:: II, JJ, KK, N1_glb, N2_glb, N3_glb, N1_loc, N2_loc, N3_loc, ELEMENT_SIZE
    INTEGER(KIND=MPI_ADDRESS_KIND):: EXTEND_SIZE
    INTEGER:: SUBCOMM_L, SUBCOMM_R, MPI_Datatype
    INTEGER:: SUBARRAY_TYPE, SUBARRAY_TYPE_resized, DISPLACEMENT_loc, RECVCOUNT_loc
    integer,DIMENSION(:),ALLOCATABLE:: DISPLACEMENT, RECVCOUNT
    COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: LOCAL_ARRAY2, SUBGLOBAL_ARRAY, GLOBAL_ARRAY_COPY

    N1_glb = SIZE(GLOBAL_ARRAY,1)
    N2_glb = SIZE(GLOBAL_ARRAY,2)
    N3_glb = SIZE(GLOBAL_ARRAY,3)

    N1_loc = SIZE(LOCAL_ARRAY,1)
    N2_loc = SIZE(LOCAL_ARRAY,2)
    N3_loc = SIZE(LOCAL_ARRAY,3)

    MPI_Datatype = MPI_DOUBLE_COMPLEX
    IF (CP8_SIZE.EQ.0) THEN
    CALL MPI_TYPE_SIZE(MPI_Double_Complex,CP8_SIZE,IERR)
    ENDIF
    ELEMENT_SIZE = CP8_SIZE

    ! SUB-COMMUNICATOR GROUP DEFINEMENT
    IF (axis.EQ.1) THEN ! axis,SUBCOMMS_L,SUBCOMMS_R => FFF,PFF
    ! PFF: NDIMR    , NTCHOPDIM/N2, NXCHOPDIM/N1 << PFF
    ! FFF: NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1 << FFF

        SUBCOMM_L = SUBCOMM_2
        SUBCOMM_R = SUBCOMM_1

    ELSEIF (axis.EQ.2) THEN ! SUBCOMMS_L,axis,SUBCOMMS_R => PPP,PFP
    ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2 << PPP
    ! PFP: NDIMR/N1, NDIMTH, NDIMX/N2 << PFP

        SUBCOMM_L = SUBCOMM_1 ! >>> SHOULD BE SUBCOMM_2
        SUBCOMM_R = SUBCOMM_2 ! >>> SHOULD BE SUBCOMM_1

    ELSE
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'MASSEMBLE: undefined for the current axis'
        RETURN
    ENDIF

    ! UNPACK GLOBAL ARRAY IN GLOBAL PROC #0 TO 
    ! SUBGLOBAL ARRAYS (SLABS DECOMPOSED IN THE 3RD DIM) STORED IN PROC #0 OF EACH SUBCOMM_L
    IF (local_proc(SUBCOMM_L).EQ.0) THEN

        ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMM_R)-1))
        ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMM_R)-1))
        RECVCOUNT_loc = N1_glb*N2_glb*N3_loc
        DISPLACEMENT_loc = N1_glb*N2_glb*local_index(N3_glb,SUBCOMM_R)
        CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMM_R,IERR)
        CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMM_R,IERR)

        IF (axis.EQ.1) THEN ! axis,SUBCOMMS_L,SUBCOMMS_R

            ! UNPACK GLOBAL ARRAY:
            ALLOCATE(SUBGLOBAL_ARRAY(N1_glb,N2_glb,N3_loc))
            CALL MPI_SCATTERV(GLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, MPI_Datatype, &
                            SUBGLOBAL_ARRAY, N1_glb*N2_glb*N3_loc, MPI_Datatype, &
                            0, SUBCOMM_R, IERR)

        ELSEIF (axis.EQ.2) THEN ! SUBCOMMS_L,axis,SUBCOMMS_R

            ! REORDER:
            ! SUBCOMMS_L,axis,SUBCOMMS_R -> axis,SUBCOMMS_L,SUBCOMMS_R
            ALLOCATE(GLOBAL_ARRAY_COPY(N2_glb,N1_glb,N3_glb))
            GLOBAL_ARRAY_COPY = RESHAPE(GLOBAL_ARRAY,SHAPE(GLOBAL_ARRAY_COPY),ORDER = [2,1,3])

            ! UNPACK GLOBAL ARRAY:
            ALLOCATE(SUBGLOBAL_ARRAY(N2_glb,N1_glb,N3_loc))
            CALL MPI_SCATTERV(GLOBAL_ARRAY_COPY, RECVCOUNT, DISPLACEMENT, MPI_Datatype, &
                            SUBGLOBAL_ARRAY, N1_glb*N2_glb*N3_loc, MPI_Datatype, &
                            0, SUBCOMM_R, IERR)
            DEALLOCATE(GLOBAL_ARRAY_COPY)

        ENDIF

         DEALLOCATE(RECVCOUNT,DISPLACEMENT)

    ENDIF

    IF (axis.EQ.1) THEN ! axis,SUBCOMMS_L,SUBCOMMS_R

        ! CUSTOM VECTOR TYPE:
        ! N2_loc x Vector {1:N1_glb,n2,1:N3_loc} per PROC
        CALL MPI_TYPE_VECTOR(N3_loc, N1_glb, N2_glb*N1_glb, &
                            MPI_Datatype, SUBARRAY_TYPE, IERR)
        EXTEND_SIZE = ELEMENT_SIZE*N1_glb
        RECVCOUNT_loc = N2_loc
        DISPLACEMENT_loc = local_index(N2_glb,SUBCOMM_L)

        ! RECEIVED DATA ORDER:
        ! {1:N1_glb,1,1:N3_loc}, ... , {1:N1_glb,N2_loc,1:N3_loc}
        ALLOCATE(LOCAL_ARRAY2(N1_glb,N3_loc,N2_loc))

    ELSE

        ! CUSTOM VECTOR TYPE:
        ! N3_loc x Vector {1:N2_glb,n1,1:N3_loc} per PROC
        CALL MPI_TYPE_VECTOR(N3_loc, N2_glb, N2_glb*N1_glb, &
                            MPI_Datatype, SUBARRAY_TYPE, IERR)
        EXTEND_SIZE = ELEMENT_SIZE*N2_glb
        RECVCOUNT_loc = N1_loc
        DISPLACEMENT_loc = local_index(N1_glb,SUBCOMM_L)

        ! RECEIVED DATA ORDER:
        ! {1:N2_glb,1,1:N3_loc}, ... , {1:N2_glb,N1_loc,1:N3_loc}
        ALLOCATE(LOCAL_ARRAY2(N2_glb,N3_loc,N1_loc))

    ENDIF

    ! COMMIT VECTOR TYPE
    CALL MPI_TYPE_CREATE_RESIZED(SUBARRAY_TYPE, INT(0, MPI_ADDRESS_KIND), EXTEND_SIZE, SUBARRAY_TYPE_resized, IERR)
    CALL MPI_TYPE_COMMIT(SUBARRAY_TYPE_resized,IERR)

    ! UNPACK SUBGLOBAL ARRAY (SLABS DECOMPOSED IN THE 3RD DIM) IN EA SUBCOMM_L's PROC #0 TO 
    ! LOCAL ARRAYS (PENCILS) STORED IN INDIVIDUAL PROCs IN THAT SUBCOMM_L
    ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMM_L)-1))
    ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMM_L)-1))
    CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMM_L,IERR)
    CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMM_L,IERR)
    CALL MPI_SCATTERV(SUBGLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, SUBARRAY_TYPE_resized, &
                    LOCAL_ARRAY2, N1_loc*N2_loc*N3_loc, MPI_Datatype, &
                    0, SUBCOMM_L, IERR)

    IF (axis.EQ.1) THEN
        !               #1   #2   #3
        ! LOCAL_ARRAY2: N1 X N3 X N2
        ! LOCAL_ARRAY : N1 X N2 X N3
        LOCAL_ARRAY = RESHAPE(LOCAL_ARRAY2,SHAPE(LOCAL_ARRAY),ORDER = [1,3,2])

    ELSEIF (axis.EQ.2) THEN
        !               #1   #2   #3
        ! LOCAL_ARRAY2: N2 X N3 X N1
        ! LOCAL_ARRAY : N1 X N2 X N3
        LOCAL_ARRAY = RESHAPE(LOCAL_ARRAY2,SHAPE(LOCAL_ARRAY),ORDER = [2,3,1])

        ! DO KK = 1,N1_loc
        !     DO JJ = 1,N3_loc
        !         DO II = 1,N2_glb
        !             LOCAL_ARRAY(KK,II,JJ) = LOCAL_ARRAY2(II,JJ,KK)
        !         ENDDO
        !     ENDDO
        ! ENDDO

    ENDIF

    CALL MPI_TYPE_FREE(SUBARRAY_TYPE_resized,IERR)
    IF (ALLOCATED(SUBGLOBAL_ARRAY)) DEALLOCATE(SUBGLOBAL_ARRAY)
    DEALLOCATE(LOCAL_ARRAY2,RECVCOUNT,DISPLACEMENT)

END SUBROUTINE MDISASSEMBLE
! ======================================================================
SUBROUTINE MSAVE0(A,FN,GLB)
!=======================================================================
! [USAGE]: 
! WRAPPER OF MSAVEX. PASSING ARGUMENTS A AND FN INTO MSAVEX
! [PARAMETERS]:
! A >> SCALAR-TYPE VARIABLE TO SAVE
! FN >> FILENAME
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    CHARACTER(LEN=*),INTENT(IN):: FN
    LOGICAL, OPTIONAL:: GLB

    IF (PRESENT(GLB)) THEN
        CALL MSAVEX(A,FN)
        WRITE(*,*) 'WRITTEN TO ',TRIM(FN)
    ELSE
        CALL MPISAVEX(A,FN)
    ENDIF

    RETURN
END SUBROUTINE MSAVE0
!=======================================================================
SUBROUTINE MSAVEX(A,FN)
!=======================================================================
! [USAGE]: 
! SAVE THE GLOBAL SCALAR-TYPE VARIABLE A INTO FN
! [PARAMETERS]:
! A >> SCALAR-TYPE VARIABLE TO SAVE
! FN >> FILENAME
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
! MUST BE CALLED IN ROOT PROC
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    CHARACTER(LEN=*):: FN

    INTEGER:: STATUS,MM,KK,I,PROC_NUM

    CALL MPI_COMM_RANK(MPI_COMM_IVP, PROC_NUM, IERR)
    IF (PROC_NUM.NE.0) THEN
        WRITE(*,*) 'MSAVEX_IN_SCALAR3: NOT IN ROOT PROC'
        STOP
    ENDIF

    OPEN(UNIT=7,FILE=FN,STATUS='UNKNOWN',&
        FORM='UNFORMATTED',IOSTAT=STATUS)

    IF(STATUS.NE.0) THEN
        WRITE(*,*) 'MSAVEX_IN_SCALAR3: FAILED TO OPEN ',TRIM(FN)
        STOP
    ENDIF

    WRITE(7) NR,NTH,NX
    WRITE(7) NRCHOP,NTCHOP,NXCHOP
    WRITE(7) A%SPACE,A%LN
    WRITE(7) ZLEN,ELL
    WRITE(7) MINC,MKLINK
    WRITE(7) SIZE(A%E,1),SIZE(A%E,2),SIZE(A%E,3)

    ! IF(A%SPACE.EQ.FFF_SPACE) THEN
    ! DO MM=1,NTCHOP
    !     WRITE(7) A%E(:NRCHOPS(MM),MM,:NXCHOPDIM)
    ! ENDDO
    ! ELSE
    ! WRITE(7) A%E(:NR,:NTH,:NX)
    ! ENDIF

    WRITE(7) A%E(:,:,:)

    CLOSE(7)

    RETURN
END SUBROUTINE MSAVEX
!=======================================================================
SUBROUTINE MPISAVEX(local_scalar,FN)
! ======================================================================
! [USAGE]: 
! SAVE THE LOCAL SCALAR-TYPE VARIABLE local_scalar INTO FN
! [PARAMETERS]:
! local_scalar >> SCALAR-TYPE VARIABLE TO SAVE
! FN >> FILENAME
! [UPDATES]:
! WRITTEN BY JINGE WANG @ SEP 22 2021
! SAVE TWO FILES:
! >>> FN.info for scalar info
! >>> FN for local_scalar%E
!=======================================================================
    TYPE(SCALAR):: local_scalar
    CHARACTER(LEN=*):: FN

    INTEGER:: PROC_NUM 
    INTEGER:: ARRAYSIZE, FILEBLK, LOCAL_ARRAYSIZE
    INTEGER:: STATUS, MPIFILE, INDEX, MPISTATUS(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND):: DISPLACEMENT

    ! 1. save basic info:
    CALL MPI_COMM_RANK(MPI_COMM_IVP, PROC_NUM, IERR)
    IF (PROC_NUM.EQ.0) THEN
    
        OPEN(UNIT=7,FILE=TRIM(FN)//'.info',STATUS='UNKNOWN',&
        FORM='UNFORMATTED',IOSTAT=STATUS)
    
        IF(STATUS.NE.0) THEN
            WRITE(*,*) 'MPISAVEX_IN_SCALAR3: ROOT PROC FAILED TO OPEN ',TRIM(FN),'.info'
            STOP
        ENDIF
    
        WRITE(7) NDIMR,NDIMTH,NDIMX
        WRITE(7) NRCHOPDIM,NTCHOPDIM,NXCHOPDIM
        WRITE(7) local_scalar%SPACE,local_scalar%LN
        WRITE(7) ZLEN,ELL
        WRITE(7) MINC,MKLINK
    
        CLOSE(7)
        WRITE(*,*) 'WRITTEN TO ',TRIM(FN)
    ENDIF

    ! 2. determine element size:
    ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
    ! PFP: NDIMR/N1, NDIMTH, NDIMX/N2
    ! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1
    ! FFF: NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1
    SELECT CASE(local_scalar%SPACE)
    CASE (PFF_SPACE)
        ARRAYSIZE = NDIMR * CEILING(REAL(NTCHOPDIM)/count_proc(SUBCOMM_2)) * CEILING(REAL(NXCHOPDIM)/count_proc(SUBCOMM_1))
    CASE (FFF_SPACE)
        ARRAYSIZE = NRCHOPDIM * CEILING(REAL(NTCHOPDIM)/count_proc(SUBCOMM_2)) * CEILING(REAL(NXCHOPDIM)/count_proc(SUBCOMM_1))
    ! PPP_SPACE OR PFP_SPACE:
    CASE DEFAULT
        ARRAYSIZE = CEILING(REAL(NDIMR)/count_proc(SUBCOMM_1)) * NDIMTH * CEILING(REAL(NDIMX)/count_proc(SUBCOMM_2))
    END SELECT
    LOCAL_ARRAYSIZE = SIZE(local_scalar%E,1)*SIZE(local_scalar%E,2)*SIZE(local_scalar%E,3)
    
    ! 3. determine Processor index
    IF (NXCHOPDIM.LT.NTCHOPDIM) THEN
        INDEX = local_proc(SUBCOMM_1) + local_proc(SUBCOMM_2) * count_proc(SUBCOMM_1)
    ELSE
        INDEX = local_proc(SUBCOMM_2) + local_proc(SUBCOMM_1) * count_proc(SUBCOMM_2)
    ENDIF

    ! 4. determine element size of MPI_DOUBLE_COMPLEX
    IF (CP8_SIZE.EQ.0) THEN
        CALL MPI_TYPE_SIZE(MPI_Double_Complex,CP8_SIZE,IERR)
    ENDIF
    
    ! 5. create FILEBLK (derived datatype)
    call MPI_TYPE_CONTIGUOUS(ARRAYSIZE, MPI_DOUBLE_COMPLEX, FILEBLK, IERR)
    call MPI_TYPE_COMMIT(FILEBLK, IERR)
    
    ! 4. save A%E in each proc
    DISPLACEMENT = INDEX*ARRAYSIZE*CP8_SIZE
    CALL MPI_FILE_OPEN(MPI_COMM_IVP, FN, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, MPIFILE, IERR)
    CALL MPI_FILE_SET_VIEW(MPIFILE, DISPLACEMENT, MPI_DOUBLE_COMPLEX, FILEBLK, &
                        'NATIVE', MPI_INFO_NULL, IERR)
    CALL MPI_FILE_WRITE_ALL(MPIFILE, local_scalar%E, LOCAL_ARRAYSIZE, MPI_DOUBLE_COMPLEX, MPISTATUS, IERR)
    CALL MPI_FILE_CLOSE(MPIFILE, IERR)
    CALL MPI_TYPE_FREE(FILEBLK, IERR)

END SUBROUTINE MPISAVEX
! ======================================================================
SUBROUTINE MLOAD0(FN,A,GLB)
!=======================================================================
! [USAGE]: 
! WRAPPER OF MLOADX. PASSING ARGUMENTS FN AND A INTO MLOADX
! [PARAMETERS]:
! FN >> FILENAME
! A >> SCALAR-TYPE VARIABLE WHERE THE LOADED VALUES FROM FN ARE STROED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*):: FN
    TYPE(SCALAR):: A
    LOGICAL, OPTIONAL:: GLB

    IF (PRESENT(GLB)) THEN
        CALL MLOADX(FN,A)
        WRITE(*,*) 'READ FROM ',TRIM(FN)
    ELSE
        CALL MPILOADX(FN,A)
    ENDIF

    RETURN
END SUBROUTINE MLOAD0
!=======================================================================
SUBROUTINE MLOADX(FN,A)
!=======================================================================
! [USAGE]: 
! READ VALUES FROM FN AND THEN STORE THEM INTO A SCALR-TYPE VARIABLE A
! [PARAMETERS]:
! FN >> FILENAME
! A >> SCALAR-TYPE VARIABLE WHERE THE LOADED VALUES FROM FN ARE STROED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*):: FN
    TYPE(SCALAR):: A

    INTEGER:: STATUS,MM
    INTEGER:: IZ,ITH,IX
    INTEGER:: IZCHOP,ITCHOP,IXCHOP
    REAL(P8):: IZLEN,IELL
    INTEGER:: IINC,IKLINK,IOF
    INTEGER:: SIZE1,SIZE2,SIZE3

    OPEN(UNIT=7,FILE=FN,STATUS='OLD',&
        FORM='UNFORMATTED',IOSTAT=STATUS)

    IF(STATUS.NE.0) THEN
    WRITE(*,*) 'MLOADX_IN_SCALAR3: FAILED TO OPEN ',TRIM(FN)
    STOP
    ENDIF

    READ(7) IZ,ITH,IX
    READ(7) IZCHOP,ITCHOP,IXCHOP
    READ(7) A%SPACE,A%LN
    READ(7) IZLEN,IELL
    READ(7) IINC,IKLINK
    READ(7) SIZE1,SIZE2,SIZE3

    ! IF(A%SPACE.EQ.FFF_SPACE) THEN
    ! IOF = IZCHOP-NRCHOP
    ! CALL CHOPSET(IOF)
    ! A%E=0
    ! DO MM=1,ITCHOP
    ! READ(7) A%E(:NRCHOPS(MM),MM,:NXCHOPDIM)
    ! ENDDO
    ! CALL CHOPSET(-IOF)
    ! CALL CHOPDO(A)
    ! ELSE
    ! READ(7) A%E(:IZ,:ITH,:IX)
    ! ENDIF

    NULLIFY(A%E)
    ALLOCATE(A%E(SIZE1,SIZE2,SIZE3))
    A%E = 0.D0
    READ(7) A%E(:,:,:)
    A%INR = 0
    A%INTH = 0
    A%INX = 0
    A%IS_ALLOCATED = .TRUE.

    CLOSE(7)

    IF(IZLEN.NE.ZLEN) THEN
    WRITE(*,*) 'MLOADX: ZLEN INCONSISTENT.'
    WRITE(*,*) 'PROG =',ZLEN
    WRITE(*,*) 'FILE =',IZLEN
    ENDIF

    IF(IELL.NE.ELL) THEN
    WRITE(*,*) 'MLOADX: ELL INCONSISTENT.'
    WRITE(*,*) 'PROG =',ELL
    WRITE(*,*) 'FILE =',IELL
    ENDIF

    IF(MINC.NE.IINC) THEN
    WRITE(*,*) 'MLOADX: MINC INCONSISTENT.'
    WRITE(*,*) 'PROG =',MINC
    WRITE(*,*) 'FILE =',IINC
    ENDIF

    IF(MKLINK.NE.IKLINK) THEN
    WRITE(*,*) 'MLOADX: MKLINK INCONSISTENT.'
    WRITE(*,*) 'PROG =',MKLINK
    WRITE(*,*) 'FILE =',IKLINK
    ENDIF

    IF (NRCHOP.NE.IZCHOP .OR. NTCHOP.NE.ITCHOP .OR. NXCHOP.NE.IXCHOP)&
    THEN
    WRITE(*,*) 'MLOADX: CHOPPING DATA INCONSISTENT.'
    WRITE(*,*) 'PROG NR,NT,NXCHOP=',NRCHOP,NTCHOP,NXCHOP
    WRITE(*,*) 'FILE NR,NT,NXCHOP=',IZCHOP,ITCHOP,IXCHOP
    ENDIF

    RETURN
END SUBROUTINE MLOADX
! ======================================================================
SUBROUTINE MPILOADX(FN,local_scalar)
! ======================================================================
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    TYPE(SCALAR):: local_scalar
    CHARACTER(LEN=*):: FN

    INTEGER:: PROC_NUM 
    INTEGER:: ARRAYSIZE, FILEBLK, LOCAL_ARRAYSIZE
    INTEGER:: STATUS, MPIFILE, INDEX, MPISTATUS(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND):: DISPLACEMENT

    INTEGER:: INDIMR, INDIMTH, INDIMX, INRCHOPDIM, INTCHOPDIM, INXCHOPDIM
    INTEGER:: IMINC, IMKLINK
    REAL(P8):: IZLEN, IELL
    
    ! 1. root proc collects all info and bcast to all procs
    CALL MPI_COMM_RANK(MPI_COMM_IVP, PROC_NUM, IERR)
    IF (PROC_NUM.EQ.0) THEN
    
        OPEN(UNIT=7,FILE=TRIM(FN)//'.info',STATUS='OLD',&
        FORM='UNFORMATTED',IOSTAT=STATUS)
    
        IF(STATUS.NE.0) THEN
            WRITE(*,*) 'MSAVEX_IN_SCALAR3: ROOT PROC FAILED TO READ ',TRIM(FN),'.info'
            CALL MPI_ABORT(MPI_COMM_IVP,ERR_FLAGS%MPIERROR,IERR)
        ENDIF
    
        READ(7) INDIMR,INDIMTH,INDIMX
        READ(7) INRCHOPDIM,INTCHOPDIM,INXCHOPDIM
        READ(7) local_scalar%SPACE,local_scalar%LN
        READ(7) IZLEN,IELL
        READ(7) IMINC,IMKLINK
        
        CLOSE(7)
        WRITE(*,*) 'MPIREAD FROM ',TRIM(FN)

        ! Check consistency
        IF(IZLEN.NE.ZLEN) THEN
            WRITE(*,*) 'MLOADX: ZLEN INCONSISTENT.'
            WRITE(*,*) 'PROG =',ZLEN
            WRITE(*,*) 'FILE =',IZLEN
        ENDIF

        IF(IELL.NE.ELL) THEN
            WRITE(*,*) 'MLOADX: ELL INCONSISTENT.'
            WRITE(*,*) 'PROG =',ELL
            WRITE(*,*) 'FILE =',IELL
        ENDIF

        IF(MINC.NE.IMINC) THEN
            WRITE(*,*) 'MLOADX: MINC INCONSISTENT.'
            WRITE(*,*) 'PROG =',MINC
            WRITE(*,*) 'FILE =',IMINC
        ENDIF

        IF(MKLINK.NE.IMKLINK) THEN
            WRITE(*,*) 'MLOADX: MKLINK INCONSISTENT.'
            WRITE(*,*) 'PROG =',MKLINK
            WRITE(*,*) 'FILE =',IMKLINK
        ENDIF
        ! ! DEBUG:
        ! WRITE(*,*) IZLEN
        ! WRITE(*,*) IELL
        ! WRITE(*,*) INDIMR,INDIMTH,INDIMX
        ! WRITE(*,*) local_scalar%SPACE,local_scalar%LN
    ENDIF
    CALL MPI_BCAST(INDIMR,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
    CALL MPI_BCAST(INDIMTH,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
    CALL MPI_BCAST(INDIMX,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
    CALL MPI_BCAST(INRCHOPDIM,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
    CALL MPI_BCAST(INTCHOPDIM,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
    CALL MPI_BCAST(INXCHOPDIM,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
    CALL MPI_BCAST(local_scalar%SPACE,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)

    ! 2. determine element size
    IF (ASSOCIATED(local_scalar%E)) NULLIFY(local_scalar%E)
    local_scalar%INR = 0
    local_scalar%INTH = 0
    local_scalar%INX = 0

    SELECT CASE(local_scalar%SPACE)
    CASE (PFF_SPACE)
        ARRAYSIZE = INDIMR * CEILING(REAL(INTCHOPDIM)/count_proc(SUBCOMM_2)) * CEILING(REAL(INXCHOPDIM)/count_proc(SUBCOMM_1))
        ALLOCATE(local_scalar%E(INDIMR,local_size(INTCHOPDIM,SUBCOMM_2),local_size(INXCHOPDIM,SUBCOMM_1)))
        local_scalar%INTH = local_index(INTCHOPDIM,SUBCOMM_2)
        local_scalar%INX  = local_index(INXCHOPDIM,SUBCOMM_1)

    CASE (FFF_SPACE)
        ARRAYSIZE = INRCHOPDIM * CEILING(REAL(INTCHOPDIM)/count_proc(SUBCOMM_2)) * CEILING(REAL(INXCHOPDIM)/count_proc(SUBCOMM_1))
        ALLOCATE(local_scalar%E(INRCHOPDIM,local_size(INTCHOPDIM,SUBCOMM_2),local_size(INXCHOPDIM,SUBCOMM_1)))
        local_scalar%INTH = local_index(INTCHOPDIM,SUBCOMM_2)
        local_scalar%INX  = local_index(INXCHOPDIM,SUBCOMM_1)

    ! PPP_SPACE OR PFP_SPACE:
    CASE DEFAULT
        ARRAYSIZE = CEILING(REAL(INDIMR)/count_proc(SUBCOMM_1)) * INDIMTH * CEILING(REAL(INDIMX)/count_proc(SUBCOMM_2))
        ALLOCATE(local_scalar%E(local_size(INDIMR,SUBCOMM_1),INDIMTH,local_size(INDIMX,SUBCOMM_2)))
        local_scalar%INR  = local_index(INDIMR,SUBCOMM_1)
        local_scalar%INX  = local_index(INDIMX,SUBCOMM_2)

    END SELECT
    LOCAL_ARRAYSIZE = SIZE(local_scalar%E,1)*SIZE(local_scalar%E,2)*SIZE(local_scalar%E,3)
    
    ! 3. determine Processor index
    IF (NXCHOPDIM.LT.NTCHOPDIM) THEN
        INDEX = local_proc(SUBCOMM_1) + local_proc(SUBCOMM_2) * count_proc(SUBCOMM_1)
    ELSE
        INDEX = local_proc(SUBCOMM_2) + local_proc(SUBCOMM_1) * count_proc(SUBCOMM_2)
    ENDIF

    ! 4. determine element size of MPI_DOUBLE_COMPLEX
    IF (CP8_SIZE.EQ.0) THEN
        CALL MPI_TYPE_SIZE(MPI_Double_Complex,CP8_SIZE,IERR)
    ENDIF
    
    ! 5. create FILEBLK (derived datatype)
    call MPI_TYPE_CONTIGUOUS(ARRAYSIZE, MPI_DOUBLE_COMPLEX, FILEBLK, IERR)
    call MPI_TYPE_COMMIT(FILEBLK, IERR)
    
    ! 6. save A%E in each proc
    CALL MPI_FILE_OPEN(MPI_COMM_IVP, FN, MPI_MODE_RDONLY, MPI_INFO_NULL, MPIFILE, IERR)
    DISPLACEMENT = INDEX*ARRAYSIZE*CP8_SIZE
    CALL MPI_FILE_SET_VIEW(MPIFILE, DISPLACEMENT, MPI_DOUBLE_COMPLEX, FILEBLK, &
                        'NATIVE', MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ_ALL(MPIFILE, local_scalar%E, LOCAL_ARRAYSIZE, MPI_DOUBLE_COMPLEX, MPISTATUS, IERR)
    CALL MPI_FILE_CLOSE(MPIFILE, IERR)
    CALL MPI_TYPE_FREE(FILEBLK, IERR)
    
    local_scalar%IS_ALLOCATED = .TRUE.

END SUBROUTINE MPILOADX
! ======================================================================
SUBROUTINE MPRINT(MESSAGE)
!=======================================================================
! [USAGE]:
! PAUSE ALL MPI RANK AND PRINT MESSAGE. MUST BE CALLED BY ALL RANKS.
! [PARAMETERS]:
! MESSAGE >> A CHARACTER STRING TO BE PRINTED
!=======================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    CALL MPI_BARRIER(MPI_COMM_IVP,IERR)
    IF (MPI_RANK == 0) THEN
    WRITE(*,*) TRIM(MESSAGE)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_IVP,IERR)

    RETURN
END SUBROUTINE MPRINT

! ======================================================================
!                           UTILITY FUNCTIONS                           
! ======================================================================
    function local_size(NSIZE,COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL SIZE
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, NSIZE
    INTEGER:: local_size

    INTEGER:: NPROC_COMM, RANK_COMM, INDEX_COMM 

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    CALL MPI_COMM_RANK(COMM,RANK_COMM,IERR)
    CALL DECOMPOSE(NSIZE,NPROC_COMM,RANK_COMM,local_size,INDEX_COMM)

    end function local_size
! ======================================================================
function local_index(NSIZE,COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL INDEX
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! 2. INDEX COUNTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, NSIZE
    INTEGER:: local_index

    INTEGER:: NPROC_COMM, RANK_COMM, SIZE_COMM 

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    CALL MPI_COMM_RANK(COMM,RANK_COMM,IERR)
    CALL DECOMPOSE(NSIZE,NPROC_COMM,RANK_COMM,SIZE_COMM,local_index)

    end function local_index
! ======================================================================
function local_proc(COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL MPI RANK IN THAT COMM GROUP
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! 2. RANK COUNTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM
    INTEGER:: local_proc

    CALL MPI_COMM_RANK(COMM,local_proc,IERR)

    end function local_proc
! ======================================================================
    function count_proc(COMM)
! ======================================================================
! [USAGE]:
! COUNT THE NUMBER OF PROCS IN THAT PROC_GROUP
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, count_proc 

    CALL MPI_COMM_SIZE(COMM, count_proc, IERR)

    end function count_proc
! ======================================================================
    function create_new_type3D(COMM,DATA_SIZE_PROC_OLD,DIM_OLD,DATA_TYPE)
! ======================================================================
! [USAGE]:
! CREATE SUBARRAY DATATYPE FOR 3D ARRAY OF DATA_SIZE_PROC_OLD STORED IN
! COMMUNICATOR GROUP: COMM
! [NOTE]:
! COMPANION FUNC FOR EXCHANGE_3DCOMPLEX_FAST
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NDIM = 3
    INTEGER:: COMM, DIM_OLD, DATA_TYPE
    INTEGER,DIMENSION(3),INTENT(IN):: DATA_SIZE_PROC_OLD
    INTEGER,DIMENSION(:),ALLOCATABLE:: create_new_type3D

    INTEGER:: I , NPROC_COMM

    ! DETERMINE THE NUMBER OF PROCS IN PROCESSOR GROUP: COMM
    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    ALLOCATE(create_new_type3D(0:NPROC_COMM-1))

    ! CREATE SUBARRAY DATATYPES
    ! OLD LOCAL ARRAY: ARRAY_OLD OF SIZE_OLD IS DECOMPOSED INTO SUBARRAY ALONG DIM_OLD
    CALL SUBARRAY(DATA_TYPE,NDIM,DATA_SIZE_PROC_OLD,DIM_OLD,NPROC_COMM,create_new_type3D)

    end function create_new_type3D
! ======================================================================
END MODULE MOD_FFT