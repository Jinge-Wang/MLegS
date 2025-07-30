MODULE MOD_SCALAR3 ! LEVEL 2 MODULE
! TYPE SETUP: SCALAR & TRANSFORM
USE omp_lib
USE MPI
USE MOD_MISC, ONLY : P4,P8,PI,IU, CP8_SIZE                         ! LEVEL 0
!XUSE USE MOD_FD                                                         ! LEVEL 1
USE MOD_EIG                                                        ! LEVEL 1
USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
!XUSE USE MOD_BANDMAT                                                    ! LEVEL 1
IMPLICIT NONE
PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
! 1.1> DIMENSIONS
INTEGER, PUBLIC:: NR,NRH                                           ! # OF COLLOCATION POINTS (= SPECTRAL COEFFICIENTS IF IN FUNC. SPACE) IN THE RADIAL DIRECTION R (OR X WHEN UNMAPPED); NRH = NR/2
INTEGER, PUBLIC:: NX                                               ! # OF COLLOCATION POINTS (= SPECTRAL COEFFICIENTS IF IN FUNC. SPACE) IN THE AXIAL DIRECTION Z
INTEGER, PUBLIC:: NTH                                              ! # OF COLLOCATION POINTS (= SPECTRAL COEFFICIENTS IF IN FUNC. SPACE) IN THE AZIMUTHAL DIRECTION THETA
INTEGER, PUBLIC:: NRCHOP,NTCHOP                                    ! CHOP INDICES IN RADIAL AND TH FOR TRUCATION IN ORDER TO LET THE MAIMUM DEGREE OF THE BASIS BE LIMITED.
INTEGER, PUBLIC,DIMENSION(:),ALLOCATABLE:: NRCHOPS,NTCHOPS         ! CHOP INDICES IN RADIAL AND TH FOR EACH N OF P^M_N
INTEGER, PUBLIC:: NXCHOP,NXCHOPH                                   ! CHOP INDEX IN Z DETERMINING THE MAXIMUM DEGREE USED IN THIS SPECTRAL ANALYSIS) ; NXCHOPH = CHOPPING LOCATION FOR THE CONJUGATE SIDE
REAL(P8),PUBLIC:: ZLEN,ZLEN0                                       ! PERIOD IN THE AXIAL DIRECTION (NOTE: Z-PERIODICITY IS ASSUMED SO THAT WE USE FOURIER SPECTRAL METHOD IN Z. THE EVP CODE CHANGES ZLEN FOR FAST CALCULATION, SO THE ORIGINAL ZLEN IS SAVED AS ZLEN0 FOR ENERGY CALCULATION.)
REAL(P8),PUBLIC:: ELL,ELL2                                         ! THE MAP PARAMTER ELL (MATSUSHIMA AND MARCUS, 1997); ELL2 = ELL**2
INTEGER, PUBLIC:: MKLINK                                           ! MKLINK = 1 OR -1 IF M AND K ARE LINKED TOGETHER. OTHERWISE, SET 0 AS DEFAULT.
INTEGER, PUBLIC:: MINC                                             ! THE NUMBER OF SYMMETRY OCCURRING IN THE AZIMUTHAL DIRECTION (E.G., SEEING 4 IDENTICAL PATTERNS IN THETA, THEN MINC = 4). DEFAULT IS 1 (NO INNER SYMMETRY).
INTEGER, PUBLIC:: NDIMR,NDIMTH,NDIMX                               ! # OF DIMENSIONS IN EACH DIRECTION
INTEGER, PUBLIC:: NRCHOPDIM,NTCHOPDIM,NXCHOPDIM                    ! # OF DIMENSIONS IN EACH DIRECTION AFTER CHOPPED
INTEGER, PUBLIC,DIMENSION(:),ALLOCATABLE:: M                       ! ACCESSIBLE AZIMUTHAL WAVENUMBERS
REAL(P8),PUBLIC,DIMENSION(:,:),ALLOCATABLE:: AK                    ! ACCESSIBLE AXIAL WAVENUMBERS GIVEN LEGENDRE ORDER N AND AZIMUTHAL WAVENUMBER M
! 1.2> MPI-SPECIFIC
INTEGER, PUBLIC:: MPI_RANK, MPI_COMM_IVP, MPI_THREAD_MODE
! INTEGER, PUBLIC:: SUBCOMM_1, SUBCOMM_2                             ! SUBCOMMUNICATOR GROUPS FOR MPI EXCHANGE
INTEGER, PUBLIC:: SUBCOMM_1 = 0, SUBCOMM_2 = 0                     ! SUBCOMMUNICATOR GROUPS FOR MPI EXCHANGE
INTEGER, PUBLIC,DIMENSION(:),ALLOCATABLE:: &
                  TYPE_PFP0,TYPE_PFP1,TYPE_PFF0,TYPE_PFF1          ! SUBARRAY DATA TYPES FOR MPI EXCHANGE
INTEGER, PUBLIC,DIMENSION(3):: &
                  SIZE_PFP0,SIZE_PFP1,SIZE_PFF0,SIZE_PFF1          ! SUBARRAY DATA TYPES FOR MPI EXCHANGE
INTEGER, PUBLIC:: IERR                                             ! MPI ERROR STATUS

! 2> TYPES
PUBLIC:: SCALAR                                                    ! ANY PHYSICAL QUANTITY IN INTEREST WILL BE DECLARED THORUGH THIS TYPE (U,V,W,P,ETC.)
TYPE SCALAR
COMPLEX(P8),DIMENSION(:,:,:),POINTER:: E                           ! PHYSICAL QUANTITY (OR FUNCTION COEFFICIENT IF IN FUNCTION SPACE) IN EACH NODE; ORDER CONVENTION = (Z, THETA, X)
REAL(P8)                            :: LN                          ! LOCARITHMIC TERM UTILIZED FOR THE TOROIDAL-POLOIDAL DECOMPOSITION (EQ. 54 OF MATSUSHIMA AND MARCUS (1997))
INTEGER                             :: SPACE                       ! TAG SHOWING IN WHICH SPACE THIS SCALAR-TYPE VALUE IS. SEE 4> SPACE TAG 
INTEGER                             :: INR,INTH,INX                ! LOCAL ARRAY STARTING INDEX
END TYPE

! 3> TRANSFORM KIT
PUBLIC:: TRANSFORM                                                 ! PHY. -> FUNC. OR FUNC -> PHY. TRANSFORMATION
TYPE TRANSFORM
REAL(P8),   DIMENSION(:),    POINTER:: R,dR                        ! COLLOCATION POINTS IN R (EVALUATED VIA GAUSS-LEGENDRE QUADRATURE)
REAL(P8),   DIMENSION(:),    POINTER:: X                           ! UNMAPPED COLLOCATION POINTS IN R (X=(R^2-ELL^2)/(R^2+ELL^2))
REAL(P8),   DIMENSION(:),    POINTER:: W                           ! GAUSS-LEGENDRE WEIGHT. W(I) IS FOR X(I)
REAL(P8),   DIMENSION(:,:),  POINTER:: NORM,LOGNORM                ! NORMALIZATION FACTORS OF P^M_N
REAL(P8),   DIMENSION(:,:,:),POINTER:: PF                          ! LEGENDRE POLYNOMIAL VALUE TABLE FOR P^M_N(X_I) GIVEN M, N AND I
REAL(P8),   DIMENSION(:),    POINTER:: AT0,AT1                     ! LEGENDRE POLYNOMIAL EVALUATION AT X=-1(R=0 -> AT0) AND X=1(R=INF -> AT1)
REAL(P8),   DIMENSION(:),    POINTER:: LN                          ! MINUS LOG VALUE OF X 
REAL(P8),   DIMENSION(:),    POINTER:: Z                           ! EQUISPACED COLLOCATION POINTS IN Z (FOR FOURIER SPECTRAL ANALYSIS)
REAL(P8),   DIMENSION(:),    POINTER:: THR,THI,TH                  ! EQUISPACED COLLOCATION POINTS IN THETA (FOR FOURIER SPECTRAL ANALYSIS)
! COMPLEX(P8),DIMENSION(:),    POINTER:: EX                          ! STORES EXP((I*2*PI/N)*J) FOR FFT
! INTEGER(P8):: NMAX                                                 ! SET TO BE MAX(NTH,NX)*2
END TYPE
TYPE(TRANSFORM),PUBLIC:: TFM                                       ! UNIVERSAL TRANSFORM KIT DECLARED IN THIS CLASS ACROSS THE CODE

! 4> SPACE TAG (KIT)
TYPE :: SPACE_MAP_ENTRY
    CHARACTER(LEN=20) :: NAME
    INTEGER :: ID
END TYPE
TYPE :: SPACE_MAP
    TYPE(SPACE_MAP_ENTRY), DIMENSION(4) :: ENTRIES
    LOGICAL :: INITIALIZED = .FALSE.
CONTAINS
    PROCEDURE :: INIT => SPACE_MAP_INIT
    PROCEDURE :: GET_ID => SPACE_MAP_GET_ID
    PROCEDURE :: GET_NAME => SPACE_MAP_GET_NAME
    PROCEDURE :: HAS_KEY => SPACE_MAP_HAS_KEY
    PROCEDURE :: SIZE => SPACE_MAP_SIZE
END TYPE
TYPE(SPACE_MAP), PUBLIC :: SPACE_MAPPER
INTEGER,PUBLIC,PARAMETER:: PPP_SPACE=0                             ! (X (OR R), THETA, Z) = (PHYS.          , PHYS.          , PHYS.          )
INTEGER,PUBLIC,PARAMETER:: FFF_SPACE=1                             ! (X (OR R), THETA, Z) = (FUNC. (MAP_LEG), FUNC. (FOURIER), FUNC. (FOURIER))
INTEGER,PUBLIC,PARAMETER:: PFP_SPACE=2                             ! (X (OR R), THETA, Z) = (PHYS.          , FUNC. (FOURIER), PHYS.          )
INTEGER,PUBLIC,PARAMETER:: PFF_SPACE=3                             ! (X (OR R), THETA, Z) = (PHYS.          , FUNC. (FOURIER), FUNC. (FOURIER))

!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
! CHOP BASED ON THE INITIALZED CHOP VALUES IN EACH DIRECTION
PUBLIC:: CHOPSET
PUBLIC:: CHOPDO
! RETURNS THE FUNCTION VALUE AT INFINITY(CALCAT1) OR ORIGIN(0) (ONLY BE CALLED BY ROOT PROC)
PUBLIC:: CALCAT1,CALCAT0                                           ! A VALUE IS RETURNED FOR EACH K.
! INTEGRATION OF A FUNCTION OVER A FULL-(~) OR HALF-(~H) DOMAIN (ONLY BE CALLED BY ROOT PROC)
PUBLIC:: INTEG                                                     ! AVAILABLE ONLY WHEN SPACE = FFF
PUBLIC:: INTEGH                                                    ! AVAILABLE ONLY WHEN SPACE = PFF
! PRODUCT & INTEGRATE F=A*B*(1-X)**2. OVER THE DOMAIN
PUBLIC:: PRODCT                                                    ! AVAILABLE ONLY WHEN SPACE = PFF
! TEST IF A SCALAR CONTAINS NAN
PUBLIC:: TEST_NAN

CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
SUBROUTINE EXSET(EX, NPTS)
!=======================================================================
! [USAGE]: 
! CALCULATE EXP((i*2*PI/NPTS)*J) INTENSIVELY USED FOR FOURIER SPACE
! WHERE J RANGES FROM -NPTS/2 TO NPTS/2
! [PARAMETERS]:
! EX >> AN ARRAY WHERE THE EXPONENTS ARE STORED
! NPTS >> THE SIZE OF THE ARRAY EX
! [UPDATES]:
! RE-CODED IN MODERN FORTRAN FORMAT BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
IMPLICIT NONE
INTEGER(P8):: NPTS
COMPLEX(P8),DIMENSION(1:NPTS):: EX

INTEGER(P8):: NH,J
REAL(P8)   :: W

W = 2*PI/FLOAT(NPTS)
NH = NPTS/2

DO J = 1,NH
  EX(J) = EXP(CMPLX(0.D0,DFLOAT(J-1)*W))
  EX(J+NH) = EXP(CMPLX(0.D0,DFLOAT(1-J)*W))
ENDDO

RETURN
END SUBROUTINE EXSET
!=======================================================================
SUBROUTINE CHOPSET(IOF)
!=======================================================================
! [USAGE]: 
! SETTING UP THE CHOP LOCATIONS FOR AZIMUTHAL AND AXIAL FUNCTION COEFFS.
! [PARAMETERS]:
! IOF >> OFFSET OF NRCHOP
! [NOTES]:
! SEE THE BELOW DIAGRAM TO UNDERSTAND HOW CHOPPING IS DONE.
!     NR-> +----------------------------------------+---------+
! NRCHOP-> +--+                                     +         +
! (NRDEG)  +  +-----+              INACTIVE         +    I    +
!          +        +-----+                         +    N    +
!          +              +-----+                   +    A    +
!          +                    +-----+ GRAD=-1     +    C    +
!          +                          +-----+       +    T    +
!          +                                +-----+ +    I    +
!          +         ACTIVE COEFFICIENTS          +-+    V    +
!          +                                        +    E    +
!          +                                        +         +
!          +----------------------------------------+---------+
!                                                   ↑         ↑ 
! [UPDATES]:                                   NTCHOP       NTH
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
INTEGER:: IOF
INTEGER:: MM,NN,NRDEG
NRCHOP = NRCHOP + IOF
NRDEG  = NRCHOP 
IF(NRCHOP.GT.NRCHOPDIM) THEN
  IF (MPI_RANK.EQ.0) THEN
  WRITE(*,*) 'CHOPSET: NRCHOP TOO LARGE.'
  WRITE(*,*) '  NRCHOP=',NRCHOP,'  NRCHOPDIM=',NRCHOPDIM
  ENDIF
STOP
ENDIF

DO MM=1,NTCHOPDIM
  NRCHOPS(MM) = MAX(MIN(NRCHOP,NRDEG-M(MM)),0)
ENDDO
NRCHOPS(NTCHOPDIM+1:)=0
DO NN=1,NR
  NTCHOPS(NN) = MAX(MIN(NTCHOP,(NRDEG-NN+MINC)/MINC),0)
ENDDO

RETURN
END SUBROUTINE CHOPSET
!=======================================================================
SUBROUTINE CHOPDO(A)
!=======================================================================
! [USAGE]: 
! PERFORM CHOPPING FUNCTION SPACE COEFFFICIENTS.
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN FUNCTION SPACE INFORMATION
!      HENCE, EXPECTED TO BE A%SPACE = FFF_SPACE
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
TYPE(SCALAR):: A
INTEGER::MM,NTCHOP_LOC,NR_LOC

! FFF_SPACE:
IF (A%SPACE.EQ.FFF_SPACE) THEN
  NTCHOP_LOC = NTCHOP-A%INTH
  ! IF (NTCHOP_LOC.GT.0) THEN
  DO MM = 1,MIN(SIZE(A%E,2),NTCHOP_LOC)
    A%E(NRCHOPS(A%INTH+MM)+1:,MM,:) = CMPLX(0.D0,0.D0)
  ENDDO
  ! ENDIF
  A%E(:,MAX(1,NTCHOP_LOC+1):,:) = CMPLX(0.D0,0.D0)
  ! DO MM=1,NTCHOP
  !   A%E(NRCHOPS(MM)+1:,MM,:)= CMPLX(0.D0,0.D0)
  ! ENDDO
  ! A%E(:,NTCHOP+1:,:) = CMPLX(0.D0,0.D0)

! PFF_SPACE:
ELSEIF (A%SPACE.EQ.PFF_SPACE) THEN
  A%E(NR+1,:,:)  = 0.D0

! PPP_SPACE OR PFP_SPACE:
ELSE
  IF (A%SPACE.EQ.PPP_SPACE) A%E(:,NTH+1,:) = 0.D0
  NR_LOC = MAX(NR - A%INR,0)
  IF (NR_LOC .LT. SIZE(A%E,1)) THEN
    A%E(NR_LOC+1:,:,:) = 0.D0
  ENDIF
  ! A%E(:,NTH+1,:) = 0.D0
  ! A%E(NR+1,:,:)  = 0.D0
ENDIF

RETURN
END SUBROUTINE CHOPDO

!=======================================================================
!============================ FUNCTIONS ================================
!=======================================================================
FUNCTION CALCAT0(A)
!=======================================================================
! [USAGE]:
! RETURNS THE FUNCTION VALUE AT ORIGIN (R = 0)
! [INPUTS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN FFF SPACE
! [OUTPUTS]:
! CALCAT0 >> FUNCTION VALUE AT ORIGIN (R = 0)
! [NOTE]:
! SHOULD ONLY CALLED BY PROCS THAT HAVE M=0 MODE
!=======================================================================
TYPE(SCALAR),INTENT(IN):: A

COMPLEX(P8),DIMENSION(SIZE(A%E,3)):: CALCAT0

IF (A%INTH.NE.0) THEN
  WRITE(*,*) 'CALCAT0: INPUT MUST CONTAIN M=0 MODE'
ENDIF

CALCAT0 = TRANSPOSE(A%E(:NRCHOP,1,:)).MUL. TFM%AT0(:NRCHOP)

RETURN
END FUNCTION CALCAT0
!=======================================================================
FUNCTION CALCAT1(A)
!=======================================================================
! [USAGE]:
! RETURNS THE FUNCTION VALUE AT INFINITY (R = INF)
! [INPUTS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN FFF SPACE
! [OUTPUTS]:
! CALCAT0 >> FUNCTION VALUE AT ORIGIN (R = INF)
! [NOTE]:
! SHOULD ONLY CALLED BY PROCS THAT HAVE M=0 MODE
! SAFE TO CALL IN EA LOCALARRAY WITH INTH = 0
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
TYPE(SCALAR),INTENT(IN):: A

COMPLEX(P8),DIMENSION(SIZE(A%E,3)):: CALCAT1

IF (A%INTH.NE.0) THEN
  WRITE(*,*) 'CALCAT1: INPUT MUST CONTAIN M=0 MODE'
ENDIF

CALCAT1 = TRANSPOSE(A%E(:NRCHOP,1,:)) .MUL. TFM%AT1(:NRCHOP)

RETURN
END FUNCTION CALCAT1
!=======================================================================
FUNCTION INTEG(F)
!=======================================================================
! [USAGE]:
! INTEGRATION OF A FUNCTION OVER A DOMAIN 
!     0  <  R   < INFTY
!     0  <  PHI < 2*PI
!     0  <  Z   < ZLEN
! CALL IN F/F SPACE
! INTEGRATED FUNCTION G IS IN THE FORM G = F*(1-MU)^2
! [INPUTS]:
! F >> A FUNTION INPUT FOR INTEGRATION (ACTUAL INTEGRATION IS DONE ON G)
! [OUTPUTS]:
! INTEG >> INTEGRATION EQUAL TO INT_[FULL_RANGE](G(R,PHI,Z) DR DPHI DZ)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 18 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
IMPLICIT NONE
TYPE(SCALAR),INTENT(IN):: F
REAL(P8)::INTEG

IF(F%SPACE.NE.FFF_SPACE) THEN
  IF (MPI_RANK.EQ.0) WRITE(*,*) 'INTEG: NOT IN FFF_SPACE'
  STOP
ENDIF      

IF ((F%INTH.EQ.0).AND.(F%INX.EQ.0)) THEN
  IF(F%LN.NE.0.0) THEN
    WRITE(*,*) 'INTEG: LOGTERM NOT ZERO'
    WRITE(*,*) 'LOGTERM=',F%LN
  ENDIF
  INTEG = 4*PI*ZLEN*ELL2*REAL(F%E(1,1,1))*TFM%NORM(1,1)
ELSE
  INTEG = 0
ENDIF
CALL MPI_ALLREDUCE(MPI_IN_PLACE, INTEG, 1, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_IVP, IERR)

RETURN
END FUNCTION INTEG
!=======================================================================
FUNCTION PRODCT(A,B)
!=======================================================================
! [USAGE]:
! CALCULATE THE PRODUCT AND INTEGRATE OVER THE DOMAIN
! CALL IN R-PHYSICAL / PHI,Z-FOURIER SPACE
! WHAT IS ACTUALLY INTEGRATED IS G = A*B*(1-MU)^2
! [INPUTS]:
! A >> FIRST SCALAR-TYPE VARIABLE FOR INTEGRATION
! B >> SECOND SCALAR-TYPE VARIABLE FOR INTEGRATION
! [OUTPUTS]:
! PRODCT >> INTEGRATION EQUAL TO INT_[FULL_RANGE](G(R,PHI,Z) DR DPHI DZ)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 18 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
IMPLICIT NONE
TYPE(SCALAR),INTENT(IN):: A,B
REAL(P8):: PRODCT

COMPLEX(P8):: PROD(NR)
INTEGER:: MM,KK

IF(A%SPACE.NE.PFF_SPACE .OR. B%SPACE.NE.PFF_SPACE) THEN
  IF (MPI_RANK.EQ.0) THEN
    WRITE(*,*) 'PRODCT:NOT IN PFF_SPACE'
    WRITE(*,*) 'A%SPACE,B%SPACE=',A%SPACE,B%SPACE
  ENDIF
  STOP
ENDIF

IF(A%LN.NE.0.0 .OR. B%LN.NE.0.0) THEN
  WRITE(*,*) 'PRODCT:LOGTERM NOT ZERO'
ENDIF

! PROD=0
! DO KK = 1,NX
!   IF(KK.GT.NXCHOP .AND. KK.LT.NXCHOPH) CYCLE
!   PROD = PROD+A%E(:NR,1,KK)*CONJG(B%E(:NR,1,KK))
! ENDDO

! DO KK = 1,NX
!   IF(KK.GT.NXCHOP .AND. KK.LT.NXCHOPH) CYCLE
!   DO MM = 2,NTCHOP
!     PROD = PROD + 2*(REAL(A%E(:NR,MM,KK))*REAL(B%E(:NR,MM,KK)) &
!                 + AIMAG(A%E(:NR,MM,KK))*AIMAG(B%E(:NR,MM,KK)))
!   ENDDO
! ENDDO

PROD=0
DO KK = 1,SIZE(A%E,3)
  DO MM = 1,SIZE(A%E,2) !NTCHOP
    IF (MM+A%INTH.EQ.1) THEN ! M(1) = 0
    PROD = PROD+A%E(:NR,1,KK)*CONJG(B%E(:NR,1,KK))
    ELSE
    PROD = PROD + 2*(REAL(A%E(:NR,MM,KK))*REAL(B%E(:NR,MM,KK)) &
                + AIMAG(A%E(:NR,MM,KK))*AIMAG(B%E(:NR,MM,KK)))
    ENDIF
  ENDDO
ENDDO
PRODCT = SUM((PROD*TFM%W)*TFM%PF(1,1,1))

CALL MPI_ALLREDUCE(MPI_IN_PLACE, PRODCT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_IVP, IERR)
PRODCT = 4*PI*ZLEN*ELL2*PRODCT*TFM%NORM(1,1)

RETURN
END FUNCTION PRODCT
!=======================================================================
FUNCTION INTEGH(F,AIN,BIN)
!=======================================================================
! [USAGE]:
! INTEGRATION OF A FUNCTION OVER A PART OF THE FULL DOMAIN 
!     0   <  R   < INFTY
!     AIN <  PHI < BIN
!     0   <  Z   < ZLEN
! CALL IN P/F SPACE, MEANING THAT F IS IN HAT FORM
! INTEGRATED FUNCTION G IS IN THE FORM G = F*(1-MU)^2
! [INPUTS]:
! F >> A FUNTION INPUT FOR INTEGRATION (ACTUAL INTEGRATION IS DONE ON G)
! AIN >> (OPTIONAL) INFIMUM OF PHI. DEFAULT IS -PI/2.
! BIN >> (OPTIONAL) SUPREMUM OF PHI. DEFAULT IS PI/2.
! [OUTPUTS]:
! INTEGH >> INTEGRATION EQUAL TO INT_[PART_RANGE](G(R,PHI,Z) DR DPHI DZ)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 18 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
IMPLICIT NONE
TYPE(SCALAR),INTENT(IN):: F
REAL(P8),INTENT(IN),OPTIONAL:: AIN, BIN
REAL(P8):: INTEGH

COMPLEX(P8):: FAB
! COMPLEX(P8),PARAMETER:: IU = (0.0D0,1.0D0) !Declared in mod_misc.f90
REAL(P8):: A,B
INTEGER:: MM,M0

IF(F%SPACE.NE.PFF_SPACE) THEN
  IF (MPI_RANK.EQ.0) WRITE(*,*) 'INTEGH: NOT IN PFF_SPACE'
  STOP
ENDIF

IF(PRESENT(AIN)) THEN
  A=AIN
ELSE
  A=-PI/2                                                          ! DEFAULT: RIGHT HAND SIDE HALF PLANE
ENDIF

IF(PRESENT(BIN)) THEN
  B=BIN
ELSE
  B=PI/2
ENDIF

! INTEGH = (B-A)*SUM(TFM%W*F%E(:NR,1,1))

! DO MM=2,NTCHOP
!   M0=M(MM)
!   FAB = (EXP(IU*M0*B)-EXP(IU*M0*A))/(IU*M0)
!   INTEGH = INTEGH+2*REAL(FAB*SUM(TFM%W*F%E(:NR,MM,1)))
! ENDDO

! INTEGH = INTEGH*ZLEN*ELL2

INTEGH = 0.D0
IF (F%INX.EQ.0) THEN
  !INTEGH = (B-A)*SUM(TFM%W*F%E(:NR,1,1))
  DO MM=1,SIZE(F%E,2) !NTCHOP
    M0=M(MM+F%INTH) ! M(1) = 0
    IF (M0.EQ.0) THEN
      INTEGH = INTEGH + (B-A)*SUM(TFM%W*F%E(:NR,1,1))
    ELSE
      FAB = (EXP(IU*M0*B)-EXP(IU*M0*A))/(IU*M0)
      INTEGH = INTEGH + 2*REAL(FAB*SUM(TFM%W*F%E(:NR,MM,1)))
    ENDIF
  ENDDO
ENDIF

CALL MPI_ALLREDUCE(MPI_IN_PLACE, INTEGH, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_IVP, IERR)
INTEGH = INTEGH*ZLEN*ELL2

RETURN
END FUNCTION INTEGH
!=======================================================================
SUBROUTINE SPACE_MAP_INIT(THIS)
!=======================================================================
! [USAGE]: Initialize the space map with name-ID pairs
!=======================================================================
CLASS(SPACE_MAP), INTENT(INOUT) :: THIS

THIS%ENTRIES(1) = SPACE_MAP_ENTRY('PPP_SPACE', PPP_SPACE)
THIS%ENTRIES(2) = SPACE_MAP_ENTRY('FFF_SPACE', FFF_SPACE)
THIS%ENTRIES(3) = SPACE_MAP_ENTRY('PFP_SPACE', PFP_SPACE)
THIS%ENTRIES(4) = SPACE_MAP_ENTRY('PFF_SPACE', PFF_SPACE)
THIS%INITIALIZED = .TRUE.

END SUBROUTINE SPACE_MAP_INIT
!=======================================================================
FUNCTION SPACE_MAP_GET_ID(THIS, NAME) RESULT(ID)
!=======================================================================
! [USAGE]: Get space ID from name (like map[key] in C++)
!=======================================================================
CLASS(SPACE_MAP), INTENT(INOUT) :: THIS
CHARACTER(LEN=*), INTENT(IN) :: NAME
INTEGER :: ID
INTEGER :: I

IF (.NOT. THIS%INITIALIZED) CALL THIS%INIT()

ID = -1  ! Default: not found
DO I = 1, SIZE(THIS%ENTRIES)
IF (TRIM(ADJUSTL(THIS%ENTRIES(I)%NAME)) == TRIM(ADJUSTL(NAME))) THEN
  ID = THIS%ENTRIES(I)%ID
  EXIT
ENDIF
ENDDO

END FUNCTION SPACE_MAP_GET_ID
!=======================================================================
FUNCTION SPACE_MAP_GET_NAME(THIS, ID) RESULT(NAME)
!=======================================================================
! [USAGE]: Get space name from ID
!=======================================================================
CLASS(SPACE_MAP), INTENT(INOUT) :: THIS
INTEGER, INTENT(IN) :: ID
CHARACTER(LEN=20) :: NAME
INTEGER :: I

IF (.NOT. THIS%INITIALIZED) CALL THIS%INIT()

NAME = 'UNKNOWN'  ! Default: not found
DO I = 1, SIZE(THIS%ENTRIES)
IF (THIS%ENTRIES(I)%ID == ID) THEN
  NAME = THIS%ENTRIES(I)%NAME
  EXIT
ENDIF
ENDDO

END FUNCTION SPACE_MAP_GET_NAME
!=======================================================================
FUNCTION SPACE_MAP_HAS_KEY(THIS, NAME) RESULT(HAS_KEY)
!=======================================================================
! [USAGE]: Check if space name exists in map
!=======================================================================
CLASS(SPACE_MAP), INTENT(INOUT) :: THIS
CHARACTER(LEN=*), INTENT(IN) :: NAME
LOGICAL :: HAS_KEY

HAS_KEY = (THIS%GET_ID(NAME) /= -1)

END FUNCTION SPACE_MAP_HAS_KEY
!=======================================================================
FUNCTION SPACE_MAP_SIZE(THIS) RESULT(MAP_SIZE)
!=======================================================================
! [USAGE]: Return number of entries in map
!=======================================================================
CLASS(SPACE_MAP), INTENT(IN) :: THIS
INTEGER :: MAP_SIZE

MAP_SIZE = SIZE(THIS%ENTRIES)

END FUNCTION SPACE_MAP_SIZE
!=======================================================================
SUBROUTINE TEST_NAN(A, VAR_NAME)
! ======================================================================
! [USAGE]:
! Check if a scalar variable contains NaN values and report their locations
! [PARAMETERS]:
! A        >> SCALAR-TYPE VARIABLE TO CHECK
! VAR_NAME >> NAME OF VARIABLE FOR REPORTING
! ======================================================================
IMPLICIT NONE
TYPE(SCALAR), INTENT(IN) :: A
CHARACTER(LEN=*), INTENT(IN) :: VAR_NAME
LOGICAL :: HAS_NAN
INTEGER :: NN, MM, KK, NAN_COUNT

! Safety check - return if array not allocated
IF (.NOT. ASSOCIATED(A%E)) THEN
  WRITE(*,'(A,A,A)') 'WARNING: Variable ', TRIM(VAR_NAME), ' not allocated!'
  RETURN
END IF

! Check for NaNs efficiently
HAS_NAN = .FALSE.
NAN_COUNT = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NN,MM,KK) REDUCTION(+:NAN_COUNT) REDUCTION(.OR.:HAS_NAN)
DO KK = 1, SIZE(A%E, 3)
  DO MM = 1, SIZE(A%E, 2)
    DO NN = 1, SIZE(A%E, 1)
      IF (ISNAN(REAL(A%E(NN,MM,KK))) .OR. ISNAN(AIMAG(A%E(NN,MM,KK)))) THEN
        HAS_NAN = .TRUE.
        NAN_COUNT = NAN_COUNT + 1
        ! Only print first few occurrences to avoid flooding output
        IF (NAN_COUNT <= 10) THEN
          WRITE(*,'(A,I4,A,A,A,3I5,A,2ES15.6)') 'Rank ', MPI_RANK, ': NaN in ', TRIM(VAR_NAME), &
            ' at (NN,MM,KK): ', NN, MM, KK, ' Value: ', A%E(NN,MM,KK)
        END IF
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Summary report
IF (HAS_NAN) THEN
  WRITE(*,'(A,I4,A,A,A,I8,A)') 'Rank ', MPI_RANK, ': ERROR: Variable ', TRIM(VAR_NAME), &
    ' contains ', NAN_COUNT, ' NaN values!'
  IF (NAN_COUNT > 10) THEN
    WRITE(*,'(A,I4,A,I8,A)') 'Rank ', MPI_RANK, ': (Only first 10 of ', NAN_COUNT, ' NaNs shown)'
  END IF
  CALL MPI_ABORT(MPI_COMM_IVP, 1, IERR)  ! Abort the program if NaNs found
END IF

! ! If no NaNs found, print confirmation
! IF (.NOT. HAS_NAN) THEN
!   WRITE(*,'(A,I4,A,A)') 'Rank ', MPI_RANK, ': Variable ', TRIM(VAR_NAME), ' has no NaN values.'
! END IF

RETURN
END SUBROUTINE TEST_NAN
!=======================================================================

END MODULE MOD_SCALAR3