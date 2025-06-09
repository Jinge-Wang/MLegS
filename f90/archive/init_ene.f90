!=======================================================================
!
!     WRITTEN BY SANGJOON (JOON) LEE
!     DEPT. OF MECHANICAL ENGINEERING
!     UNIV. OF CALIFORNIA AT BERKELEY
!     EMAIL: SANGJOONLEE@BERKELEY.EDU
!     
!     UC BERKELEY CFD LAB
!     HTTPS://CFD.ME.BERKELEY.EDU/
!
!     NONCOMMERCIAL USE WITH COPYRIGHTED (C) MARK
!     UNDER DEVELOPMENT FOR RESEARCH PURPOSE 
!
!=======================================================================
PROGRAM INIT_ENE
!=======================================================================
! [USAGE]: 
! MAKE AN INITIAL FIELD IN A POLOIDAL-TOROIDAL FORM 
! BASED ON THE INITIAL Q-VORTEX PROFILE, CREATE PSI0 AND CHI0
! [UPDATES]:
! BASED ON TATSU'S ORIGINAL CODE PAIRMAKE
! LAST UPDATE ON NOV 23, 2020
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

      IMPLICIT NONE
      TYPE(SCALAR):: A, B, C   ! EXTRA SCALAR-TYPE VARIABLES
      INTEGER     :: I, J, K   ! EXTRA INTEGER VARIABLES (FOR ITER.)
      REAL(P8)    :: X, Y, Z   ! EXTRA REAL(P8) VARIABLES

      WRITE(*,*) 'PROGRAM STARTED'
      CALL PRINT_REAL_TIME()   ! @ MOD_MISC

      ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
      CALL READIN(5)           ! @ MOD_INIT

      CALL COLLOC_INFO()

      ! NR: # OF RADIAL COLLOC. PTS
      ! NTH: # OF AZIMUTHAL COLLOC. PTS
      ! NX: # OF AXIAL COLLOC. PTS
      WRITE(*,*) 'NR,NTH,NX=',NR,NTH,NX

      DO I=1,QPAIR%N
        ! INDICATE THE INITIAL VORTEX POSITIONS AND STRENGTHS
        WRITE(6,61) 'VORTEX:',I
        WRITE(6,60) 'X=',QPAIR%X(I),'Y=',QPAIR%Y(I)
        WRITE(6,60) 'Q=',QPAIR%Q(I),'H=',QPAIR%H(I),'B=',QPAIR%B(I)
        
        ! THESE QUANTITIES ARE FOR PERTURBING THE VORTEX 
        ! K=1, R=1, REST=0 (DEFAULT)
        WRITE(6,60) 'DX=',QPAIR%DX(I),'DY=',QPAIR%DY(I)
        WRITE(6,60) 'K=' ,QPAIR%K(I), 'PHAS=',QPAIR%DP(I)
        WRITE(6,60) 'R=' ,QPAIR%R(I), 'DR=',QPAIR%DR(I)
        WRITE(6,60) 'RK=',QPAIR%RK(I),'PHAS=',QPAIR%RP(I)
      ENDDO
      WRITE(6,61) 'NOISE=',QPAIR%NOISE ! NOISE=0 (DEFAULT)

  61   FORMAT(A,I2)
  60   FORMAT(A5,1PE13.6,1X,A5,1PE13.6)
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!! STAGE 1. GENERATE PSI0 !!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL ALLOCATE(A,PPP_SPACE)

      ! INITIALIZE
      A%E = 0.D0

      DO I=1,QPAIR%N
        ! INPUT OMEGA_Z (=2Q*EXP(-R**2)) INTO GAUSDATA
        ! IN ORDER TO OBTAIN PSI0, WE USE THE FORMULA
        ! OMEGA_Z = - DELSQH PSI0      <-->
        ! PSI0    = DELSQH^(-1)(-OMEGA_Z).
        GAUSDATA%XC = QPAIR%X(I)
        GAUSDATA%YC = QPAIR%Y(I)
        GAUSDATA%Q  = QPAIR%Q(I)*2.
        GAUSDATA%B  = 1.D0
        GAUSDATA%DX = QPAIR%DX(I)
        GAUSDATA%DY = QPAIR%DY(I)
        GAUSDATA%K  = QPAIR%K(I)
        GAUSDATA%DP = QPAIR%DP(I)
        GAUSDATA%R  = QPAIR%R(I)
        GAUSDATA%DR = QPAIR%DR(I)
        GAUSDATA%RK = QPAIR%RK(I)
        GAUSDATA%RP = QPAIR%RP(I)
        
        ! CREATE A GAUSSIAN VORTEX MULTIPLIED BY (1-X)^(-2)
        ! NOW A%E HAS THE FIELD F(R,T,Z) = 2Q*EXP(-R**2)*(1-X)^(-2)
        ! WHERE X = (R^2-ELL^2)/(R^2+ELL^2)
        ! CALL LAY(A,GAUS,2)
        A%E = 0.D0
        
        ! CHECK WHETHER THE GAUSSIAN IS ACTUALLY GAUSSIAN
        ! FOR DEBUGGIN PURPOSES. TURN IT OFF IN PRACTICE
        ! CALL CHECKGAUSSIAN(A,GAUSDATA%XC,GAUSDATA%YC,2)
      ENDDO

      ! (1-X)^(-2)*OMEGA_Z (PPP) --> (1-X)^(-2)*OMEGA_Z (FFF)
      ! THE FUNCTION IS TRUNCATED INTO THE FINITE BASIS FUNCTION SET,
      ! CAUSING THE INFORMATION LOSS DUE TO APPROXIMATION, WHICH IS
      ! INEVITABLE UNLESS WE USE THE INFINITE-SIZE FULL BASIS SET
      ! STILL, THE TRANSFORMED COEFFICIENT MATRIX IN FFF SPACE
      ! IS THE BEST APPROXIMATION OF THE INPUT FIELD IN PPP SPACE.
      CALL TOFF(A)

      ! (1-X)^(-2)*OMEGA_Z (FFF) --> -(1-X)^(-2)*OMEGA_Z (FFF)
      A%E=-A%E
      
      ! MUTE A LOGTERM
      A%LN=0

      ! ADD RANDOM NOISE IF THE NOISE OPTION IS ON
      IF(QPAIR%NOISE .EQ. 1) THEN
        CALL ALLOCATE(B)
        B=A
        CALL NOISE(A,0.02_P8)
        CALL TOFP(A)
        CALL TOFP(B)
        CALL FLATTEN(A,B)
        CALL TOFF(A)
        CALL DEALLOCATE(B)
      ENDIF

      ! -(1-X)^(-2)*OMEGA_Z (FFF) -> DELSQH^(-1)(-OMEGA_Z) = PSI0 (FFF)
      ! NOTE: IDELSQH = ((1-X)^(-2) DELSQH)^(-1) = DELSQH^(-1)*((1-X)^(-2))^(-1)
      ! CREATES INVERSE HORIZONTAL DEL-SQUARE OPERATOR WITH EXTRA (1-X)^2
      ! B NOW HAS A LOG TERM ASSOCIATED WITH THE R-DERIVATIVE 
      ! IN THE INVERSE DEL SQ.
      CALL ALLOCATE(B)
      CALL IDELSQH(A,B)

      !WRITE OUT THE 3D MATRIX IN THE SCALAR USING UNFORMATTED MODE
      CALL MSAVE(B,FILES%PSI0)

      CALL DEALLOCATE(A)
      CALL DEALLOCATE(B)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!! STAGE 1. GENERATE CHI0 !!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL ALLOCATE(A,PPP_SPACE)

      ! INITIALIZE
      A%E=0

      DO I=1,QPAIR%N
        ! INPUT V_Z (=H*EXP(-B*R**2)) INTO GAUSDATA
        ! IN ORDER TO OBTAIN CHI0, WE USE THE FORMULA
        ! V_Z = - DELSQH CHI0      <-->
        ! CHI0    = DELSQH^(-1)(-V_Z).
        GAUSDATA%XC = QPAIR%X(I)
        GAUSDATA%YC = QPAIR%Y(I)
        ! GAUSDATA%Q  = QPAIR%H(I)
        GAUSDATA%Q  = QPAIR%Q(I)
        ! GAUSDATA%B  = QPAIR%B(I)
        GAUSDATA%B  = 1.D0
        GAUSDATA%DX = QPAIR%DX(I)
        GAUSDATA%DY = QPAIR%DY(I)
        GAUSDATA%K  = QPAIR%K(I)
        GAUSDATA%DP = QPAIR%DP(I)
        GAUSDATA%R  = QPAIR%R(I)
        GAUSDATA%DR = QPAIR%DR(I)
        GAUSDATA%RK = QPAIR%RK(I)
        GAUSDATA%RP = QPAIR%RP(I)

        ! CREATE A GAUSSIAN VORTEX MULTIPLIED BY (1-X)^(-2)
        ! NOW A%E HAS THE FIELD F(R,T,Z) = H*EXP(-B*R**2)*(1-X)^(-2)
        ! WHERE X = (R^2-ELL^2)/(R^2+ELL^2)
        CALL LAY(A,GAUS,2)

        ! UZ = (COS(M*TH)+2)*(1-X)
        CALL LAYP(A,COSM)

        ! CHECK WHETHER THE GAUSSIAN IS ACTUALLY GAUSSIAN
        ! FOR DEBUGGIN PURPOSES. TURN IT OFF IN PRACTICE
        ! CALL CHECKGAUSSIAN(A,GAUSDATA%XC,GAUSDATA%YC,2)
      ENDDO

      ! (1-X)^(-2)*V_Z (PPP) --> (1-X)^(-2)*V_Z (FFF)
      ! THE FUNCTION IS TRUNCATED INTO THE FINITE BASIS FUNCTION SET,
      ! CAUSING THE INFORMATION LOSS DUE TO APPROXIMATION, WHICH IS
      ! INEVITABLE UNLESS WE USE THE INFINITE-SIZE FULL BASIS SET
      ! STILL, THE TRANSFORMED COEFFICIENT MATRIX IN FFF SPACE
      ! IS THE BEST APPROXIMATION OF THE INPUT FIELD IN PPP SPACE.
      CALL TOFF(A)

      ! (1-X)^(-2)*V_Z (FFF) --> -(1-X)^(-2)*V_Z (FFF)
      A%E=-A%E

      ! MUTE A LOGTERM
      A%LN=0

      ! ADD RANDOM NOISE IF THE NOISE OPTION IS ON
      IF(QPAIR%NOISE.EQ.1) THEN
        CALL ALLOCATE(B)
        B=A
        CALL NOISE(A,0.02_P8)
        CALL TOFP(A)
        CALL TOFP(B)
        CALL FLATTEN(A,B)
        CALL TOFF(A)
        CALL DEALLOCATE(B)
      ENDIF

      ! -(1-X)^(-2)*V_Z (FFF) -> DELSQH^(-1)(-V_Z) = CHI0 (FFF)
      ! NOTE: IDELSQH = ((1-X)^(-2) DELSQH)^(-1) = DELSQH^(-1)*((1-X)^(-2))^(-1)
      ! CREATES INVERSE HORIZONTAL DEL-SQUARE OPERATOR WITH EXTRA (1-X)^2
      ! B NOW HAS A LOG TERM ASSOCIATED WITH THE R-DERIVATIVE 
      ! IN THE INVERSE DEL SQ.
      CALL ALLOCATE(B)
      CALL IDELSQH(A,B)

      !WRITE OUT THE 3D MATRIX IN THE SCALAR USING UNFORMATTED MODE
      CALL MSAVE(B,FILES%CHI0)


      CALL DEALLOCATE(A)
      CALL DEALLOCATE(B)


      WRITE(*,*) ''
      WRITE(*,*) 'PROGRAM FINISHED'
      CALL PRINT_REAL_TIME()  ! @ MOD_MISC

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================
! N/A
!=======================================================================
END PROGRAM INIT_ENE
!=======================================================================
    