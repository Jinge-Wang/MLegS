MODULE MOD_LEGOPS ! LEVEL 3 MODULE
      USE MOD_MISC, ONLY : P4,P8,PI,IU,MCAT                                   ! LEVEL 0
!XUSE USE MOD_FD                                                         ! LEVEL 1
      USE MOD_EIG                                                        ! LEVEL 1
      USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
      USE MOD_BANDMAT                                                    ! LEVEL 1
      USE MOD_SCALAR3                                                    ! LEVEL 2
      IMPLICIT NONE
      PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
      ! 1> NONLINEAR TERMS
      PUBLIC:: NONLINADD
      TYPE NONLINADD
        REAL(P8):: STRAIN                                                ! STRAIN
        REAL(P8):: ROT                                                   ! ROTATION
        REAL(P8):: U                                                     ! ADVECTIVE VELOCITY (SQRT(UX**2+UY**2))
        REAL(P8):: ANG                                                   ! ADVECTIVE VELOCITY ANGLE (ACOS(UX/U))
        REAL(P8):: UZ                                                    ! ADVECTIVE VELOCITY IN Z DIRECTION
      END TYPE
      TYPE(NONLINADD),PUBLIC:: NADD

      ! 2> VELOCITY MONITOR
      PUBLIC:: MONITOR
      TYPE MONITOR
        INTEGER:: N,INT
        REAL(P8):: URMAX,UPMAX,UZMAX
      END TYPE
      TYPE(MONITOR),PUBLIC:: VELMON
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
      ! HORIZONTAL LAPLACIAN (DEL-SQUARE-H) OPERATORS (FORWARD, INVERSE)
      PUBLIC:: DELSQH, IDELSQH
      ! (1-X)*, (1+X)*, (1+X)^(-1)* AND (1-X)^(-1)
      PUBLIC:: MULXM, MULXP, DIVXP, DIVXM
      ! [(1-X)/(1+X)]* OPERATOR
      PUBLIC:: MULXMDIVXP
      ! (1-X^2)*D/DX OPERATOR
      PUBLIC:: XXDX
      ! LAPLACIAN (DEL-SQUARE) OPERATORS (FORWARD, INVERSE, LOG-TERM)
      PUBLIC:: DEL2, IDEL2, IDEL2LN                                      ! LOG-TERM ENABLES THE SOLUTION TO SATISFY THE BC OF F->0 AS R->INFTY (EQ.54 IN MATSUSHIMA & MARCUS, 1997)
      ! HELMHOLTZ OPERATORS (FORWARD, INVERSE)
      PUBLIC:: HELM, IHELM
      ! DEL^P OPERATORS WHERE P IS EVEN AND >=4 (FORWARD, INV, INV_VER2)
      PUBLIC:: HELMP, IHELMP, IHELMPO
      ! PT-DECOMPOSED VECTOR (PSI,CHI) TO VELOCITY OR VORTICITY
      PUBLIC:: PC2VEL                                                    !  PC2VEL: (PSI,CHI) TO VELOCITY (R*U_R,R*U_PHI,U_Z)
      PUBLIC:: PC2VOR                                                    !  PC2VOR: (PSI,CHI) TO VORTICITY (R*O_R,R*O_PHI,O_Z) 
      ! COMPUTE VECTOR PRODUCT
      PUBLIC:: VPROD
      ! COMPUTE VECTOR PRODUCT with Q-vortex
      PUBLIC:: VPROD_PFF      
      ! PROJECTION OPERATOR INTO TOROIDAL AND POLOIDAL SPACE
      PUBLIC:: PROJECT
      ! NONLINEAR TERMS
      PUBLIC:: NONLIN
      ! VELOCITY(FFF) TO VORTICITY(PFF)
      PUBLIC:: VEL2VOR
CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
      SUBROUTINE DELSQH(A,B)
!=======================================================================
! [USAGE]: 
! MODIFIED HORIZONTAL LAPLACIAN (DEL-SQUARE-H) OPERATOR
! [PARAMETERS]:
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, (1-X)**(-2)*DEL^2_H(A)
! [NOTE]:
! RATHER THAN DIRECTLY CALCULATING PARTIAL DERIVATIVES TO GET THIS
! OPERATOR, THIS SUBROUTINE USES THE STURM-LIOUVILLE EQUATION OF P^M_L_N
! THAT IS DEL^2_H(P^M_L_N) = - 4N(N+1)ELL^2/(ELL^2+R^2)^2 P^M_L_N)
! OR EQUIVALENTLY 
! (ELL^2+R^2)^2/(2*ELL)^2*DEL^2_H(P^M_L_N) = -N(N+1)/ELL^2*P^M_L_N
! <-> (1-X)**(-2)*DEL^2_H(P^M_L_N) = -N(N+1)/ELL^2 * P^M_L_N
! <-> (1-X)**(-2)*DEL^2_H = - N(N+1)/ELL^2.
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR),INTENT(IN):: A
      TYPE(SCALAR),INTENT(INOUT):: B

      INTEGER:: MM,NN,N

      IF(A%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'DELSQH: INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF
     
      IF(B%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(B)
        CALL ALLOCATE(B)
      ENDIF
      
      IF(M(1).NE.0) THEN
        WRITE(*,*) 'DELSQH: M(1) NOT EQUAL 0.'
        STOP
      ENDIF

      B%LN = 0.D0
      CALL CHOPDO(B)

      DO MM=1,NTCHOP
        DO NN=1,NZCHOPS(MM)
          N = M(MM) + (NN-1)
          B%E(NN,MM,:) = -A%E(NN,MM,:)*(N*(N+1)/ELL2)
        ENDDO
      ENDDO

      B%E(1,1,1) = A%LN/ELL2/TFM%NORM(1,1)

      RETURN
      END SUBROUTINE DELSQH
!=======================================================================
      SUBROUTINE IDELSQH(B,A)
!=======================================================================
! [USAGE]: 
! MODIFIED INVERSE HORIZONTAL LAPLACIAN (DEL-SQUARE-H) OPERATOR.
! LOG-TERM RETURNS IN A%LN.
! [PARAMETERS]:
! B >> INPUT IN FFF SPACE
! A >> ON EXIT, ((1-X)**(-2)*DEL^2_H)^(-1)(B)
! [NOTE]:
! RATHER THAN DIRECTLY CALCULATING PARTIAL DERIVATIVES TO GET THIS
! OPERATOR, THIS SUBROUTINE USES THE STURM-LIOUVILLE EQUATION OF P^M_L_N
! THAT IS DEL^2_H(P^M_L_N) = - 4N(N+1)ELL^2/(ELL^2+R^2)^2 P^M_L_N)
! OR EQUIVALENTLY 
! (ELL^2+R^2)^2/(2*ELL)^2*DEL^2_H(P^M_L_N) = -N(N+1)/ELL^2*P^M_L_N
! <-> (1-X)**(-2)*DEL^2_H(P^M_L_N) = -N(N+1)/ELL^2 * P^M_L_N
! <-> (1-X)**(-2)*DEL^2_H = - N(N+1)/ELL^2.
! [DEPENDENCIES]:
! 1. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 2. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR),INTENT(IN):: B
      TYPE(SCALAR),INTENT(INOUT):: A

      INTEGER:: MM,N,NN

      IF(B%LN.NE.0.D0) THEN
        WRITE(*,*) 'IDELSQH: INPUT WITH LOGTERM'
      ENDIF

      IF(M(1).NE.0) THEN
        WRITE(*,*) 'M(1) NEQ 0'
        STOP
      ENDIF

      IF(B%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'IDELSQH:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(A%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(A)
        CALL ALLOCATE(A)
      ENDIF

      CALL CHOPDO(A)
      A%LN = B%E(1,1,1)*TFM%NORM(1,1)*ELL2 ! Find log-term
      A%E(1,1,:)=0.D0

      DO NN=2,NZCHOPS(1)
        N = NN - 1 
        A%E(NN,1,:) = -B%E(NN,1,:)*(ELL2/(N*(N+1)))
      ENDDO

      DO MM=2,NTCHOP
        DO NN=1,NZCHOPS(MM)
          N = M(MM) + (NN-1)
          A%E(NN,MM,:) = -B%E(NN,MM,:)*(ELL2/(N*(N+1)))
        ENDDO
      ENDDO

      A%E(1,1,:) = A%E(1,1,:) - CALCAT1(A)/TFM%AT1(1) ! Fix A at r=inft

      RETURN
      END SUBROUTINE IDELSQH
!=======================================================================
      SUBROUTINE MULXP(A,B)
!=======================================================================
! [USAGE]: 
! (1+X)* (OR 2R^2/(R^2+ELL^2)*) OPERATOR.
! [PARAMETERS]:
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, (1+X)*A
! [DEPENDENCIES]:
! 1. BAND_LEG_XP(~) @ MOD_LIN_LEGENDRE
! 2. BANMUL(~) @ MOD_BANDMAT
! 3. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 4. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: A
      TYPE(SCALAR),INTENT(INOUT):: B

      INTEGER:: MM,NN,N,KK
      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XP

      IF(A%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'MULXP:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(B%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(B)
        CALL ALLOCATE(B)
      ENDIF

      IF(A%LN.NE.0.D0) WRITE(*,*) 'MULXP: A%LN IS NONZERO.'

      ALLOCATE( XP(NZCHOP,3) )

      B%LN = 0.0
      CALL CHOPDO(B)

      DO MM=1,NTCHOP
        NN = NZCHOPS(MM)
        XP(:NN,:) = BAND_LEG_XP(NN,M(MM),TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          B%E(:NN,MM,KK)=BANMUL(XP(:NN,:),2,A%E(:NN,MM,KK))
        ENDDO
      ENDDO

      DEALLOCATE( XP )

      RETURN
      END SUBROUTINE MULXP
!=======================================================================
      SUBROUTINE DIVXP(B,A)
!=======================================================================
! [USAGE]: 
! (1+X)^(-1)* (OR (R^2+ELL^2)/2R^2*) OPERATOR.
! [PARAMETERS]:
! B >> INPUT IN FFF SPACE
! A >> ON EXIT, (1+X)^(-1)*B
! [DEPENDENCIES]:
! 1. BAND_LEG_XP(~) @ MOD_LIN_LEGENDRE
! 2. BANMUL(~) @ MOD_BANDMAT
! 3. LUB(~) @ MOD_BANDMAT
! 4. SOLVEB(~) @ MOD_BANDMAT
! 5. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 6. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! WRITTEN BY SANGJOON LEE @ NOV 28 2020
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: B
      TYPE(SCALAR),INTENT(INOUT):: A

      INTEGER:: MM,NN,N,KK
      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XP

      IF(B%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'DIVXP:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(A%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(B)
        CALL ALLOCATE(B)
      ENDIF

      IF(B%LN.NE.0.D0) WRITE(*,*) 'DIVXP: A%LN IS NONZERO.'
      A = B

      ALLOCATE( XP(NZCHOP,3) )

      A%LN = 0.0
      CALL CHOPDO(A)

      DO MM=1,NTCHOP
        NN = NZCHOPS(MM)
        XP(:NN,:) = BAND_LEG_XP(NN,M(MM),TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          CALL LUB(XP(:NN,:),2)
          CALL SOLVEB(XP(:NN,:),2,A%E(:NN,MM,KK))
        ENDDO
      ENDDO

      DEALLOCATE( XP )

      RETURN
      END SUBROUTINE DIVXP
!=======================================================================
      SUBROUTINE MULXM(A,B)
!=======================================================================
! [USAGE]: 
! (1-X)* (OR 2ELL^2/(R^2+ELL^2)*) OPERATOR.
! [PARAMETERS]:
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, (1-X)*A
! [DEPENDENCIES]:
! 1. BAND_LEG_XM(~) @ MOD_LIN_LEGENDRE
! 2. BANMUL(~) @ MOD_BANDMAT
! 3. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 4. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: A
      TYPE(SCALAR),INTENT(INOUT):: B

      INTEGER:: MM,NN,N,KK
      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XM

      IF(A%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'MULXM:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(B%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(B)
        CALL ALLOCATE(B)
      ENDIF

      IF(A%LN.NE.0.D0) WRITE(*,*) 'MULXM: A%LN IS NONZERO.'

      ALLOCATE( XM(NZCHOP,3) )

      B%LN = 0.0
      CALL CHOPDO(B)

      DO MM=1,NTCHOP
        NN = NZCHOPS(MM)
        XM(:NN,:) = BAND_LEG_XM(NN,M(MM),TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          B%E(:NN,MM,KK)=BANMUL(XM(:NN,:),2,A%E(:NN,MM,KK))
        ENDDO
      ENDDO

      DEALLOCATE( XM )

      RETURN
      END SUBROUTINE MULXM
!=======================================================================
      SUBROUTINE DIVXM(B,A)
!=======================================================================
! [USAGE]: 
! (1-X)^(-1)* (OR (R^2+ELL^2)/2ELL^2*) OPERATOR.
! [PARAMETERS]:
! B >> INPUT IN FFF SPACE
! A >> ON EXIT, (1-X)^(-1)*B
! [DEPENDENCIES]:
! 1. BAND_LEG_XM(~) @ MOD_LIN_LEGENDRE
! 2. BANMUL(~) @ MOD_BANDMAT
! 3. LUB(~) @ MOD_BANDMAT
! 4. SOLVEB(~) @ MOD_BANDMAT
! 5. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 6. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! WRITTEN BY SANGJOON LEE @ NOV 28 2020
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: B
      TYPE(SCALAR),INTENT(INOUT):: A

      INTEGER:: MM,NN,N,KK
      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XM

      IF(B%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'DIVXM:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(A%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(B)
        CALL ALLOCATE(B)
      ENDIF

      IF(B%LN.NE.0.D0) WRITE(*,*) 'DIVXM: A%LN IS NONZERO.'
      A = B

      ALLOCATE( XM(NZCHOP,3) )

      A%LN = 0.0
      CALL CHOPDO(A)

      DO MM=1,NTCHOP
        NN = NZCHOPS(MM)
        XM(:NN,:) = BAND_LEG_XM(NN,M(MM),TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          CALL LUB(XM(:NN,:),2)
          CALL SOLVEB(XM(:NN,:),2,A%E(:NN,MM,KK))
        ENDDO
      ENDDO

      DEALLOCATE( XM )

      RETURN
      END SUBROUTINE DIVXM
!=======================================================================
      SUBROUTINE MULXMDIVXP(A,B,SP)
!=======================================================================
! [USAGE]: 
! [(1-X)/(1+X)]* (OR [L^2/R^2]*) OPERATOR.
! [PARAMETERS]:
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, [(1-X)/(1+X)]*A
! SP >> (OPTIONAL) IF EXIST, DIVIDED BY ELL^2
! [DEPENDENCIES]:
! 1. MULXM(~) @ MOD_LEGOPS
! 2. DIVXP(~) @ MOD_LEGOPS
! 5. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! [UPDATES]:
! WRITTEN BY SANGJOON LEE @ NOV 28 2020
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: A
      TYPE(SCALAR),INTENT(INOUT):: B
      LOGICAL,OPTIONAL,INTENT(IN):: SP

      TYPE(SCALAR):: C

      CALL ALLOCATE( C )

      CALL MULXM(A,C)
      CALL DIVXP(C,B)

      IF(PRESENT(SP) .AND. SP .EQ. .TRUE.) B%E = B%E / ELL2

      CALL DEALLOCATE( C )

      RETURN
      END SUBROUTINE MULXMDIVXP
!=======================================================================
      SUBROUTINE XXDX(A,B)
!=======================================================================
! [USAGE]: 
! (1-X**2)*D/DX (OR R*D/DR) OPERATOR.
! [PARAMETERS]:
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, (1-X)**2*D/DX(A)
! [DEPENDENCIES]:
! 1. BAND_LEG_XXDX(~) @ MOD_LIN_LEGENDRE
! 2. BANMUL(~) @ MOD_BANDMAT
! 3. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 4. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: A
      TYPE(SCALAR),INTENT(INOUT):: B

      INTEGER:: MM,NN,N,KK
      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XXDXOP

      IF(A%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'XXDX:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(B%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(B)
        CALL ALLOCATE(B)
      ENDIF

      ! NOTE: temporarily muted
      IF(ABS(A%E(NZCHOPS(2),2,1)).GT.1.E-15) THEN
       WRITE(*,*) 'XXDX: OPERATION NOT EXACT'
       WRITE(*,*)  A%E(NZCHOPS(2),2,1)
      !  STOP
      ENDIF

      ALLOCATE( XXDXOP(NZCHOP,3) )

      B%LN = 0.D0
      CALL CHOPDO(B)

      DO MM=1,NTCHOP
        NN = NZCHOPS(MM)
        XXDXOP(:NN,:) = BAND_LEG_XXDX(NN,M(MM),TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          B%E(:NN,MM,KK)=BANMUL(XXDXOP(:NN,:),2,A%E(:NN,MM,KK))
        ENDDO
      ENDDO

      B%E(1,1,1) = B%E(1,1,1) + A%LN/TFM%NORM(1,1)
      B%E(2,1,1) = B%E(2,1,1) + A%LN/TFM%NORM(2,1)

      DEALLOCATE( XXDXOP )

      RETURN
      END SUBROUTINE XXDX
!=======================================================================
      SUBROUTINE DEL2(A,B)
!=======================================================================
! [USAGE]: 
! LAPLACIAN (DEL-SQUARE) OPERATOR.
! [PARAMETERS]:
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, DEL^2(A)
! [DEPENDENCIES]:
! 1. BAND_LEG_RAT_DEL2H(~) @ MOD_LIN_LEGENDRE
! 2. BANMUL(~) @ MOD_BANDMAT
! 3. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 4. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: A
      TYPE(SCALAR),INTENT(INOUT):: B

      INTEGER:: MM,NN,N,KK
      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: DEL2OP

      IF(A%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'DEL2:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(B%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(B)
        CALL ALLOCATE(B)
      ENDIF

      IF(ABS(A%E(NZCHOPS(2),2,1)).GT.1.E-15) THEN
        WRITE(*,*) 'DEL2: OPERATION NOT EXACT'
        WRITE(*,*)  A%E(NZCHOPS(2),2,1)
        ! STOP
      ENDIF

      ALLOCATE( DEL2OP(NZCHOP,5) )

      B%LN = 0.0
      CALL CHOPDO(B)

      DO MM=1,NTCHOP
        NN = NZCHOPS(MM)
        DEL2OP(:NN,:) = BAND_LEG_RAT_DEL2H(NN,M(MM),ELL,TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          B%E(:NN,MM,KK)=BANMUL(DEL2OP(:NN,:),3,A%E(:NN,MM,KK))
          B%E(:NN,MM,KK)=B%E(:NN,MM,KK)-(AK(MM,KK)**2)*A%E(:NN,MM,KK)
        ENDDO
      ENDDO

      B%E(1,1,1)=B%E(1,1,1) +4.0D0/3*A%LN/ELL2/TFM%NORM(1,1)
      B%E(2,1,1)=B%E(2,1,1) -2.0D0*A%LN/ELL2/TFM%NORM(2,1)
      B%E(3,1,1)=B%E(3,1,1) +2.0D0/3*A%LN/ELL2/TFM%NORM(3,1)

      DEALLOCATE( DEL2OP )

      RETURN
      END SUBROUTINE DEL2
!=======================================================================
      SUBROUTINE IDEL2(B,A,LN)
!=======================================================================
! [USAGE]: 
! INVERSE LAPLACIAN (DEL-SQUARE^(-1)) OPERATOR.
! LOG-TERM IS ASSUMED AS BOUNDARY CONDITION (DEFAULT SETTING IS LN=0).
! [PARAMETERS]:
! B >> INPUT IN FFF SPACE
! A >> ON EXIT, (DEL^2)^(-1)(B)
! LN >> (OPTOINAL) ASSUMED AS BOUNDARY CONDITION. (DEFAULT IS 0)
! [DEPENDENCIES]:
! 1. BAND_LEG_RAT_DEL2H(~) @ MOD_LIN_LEGENDRE
! 2. LUB(~) @ MOD_BANDMAT
! 3. SOLVEB(~) @ MOD_BANDMAT
! 4. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 5. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: A, B
      REAL(P8),DIMENSION(:,:,:),ALLOCATABLE :: DEL2OP
      REAL(P8),OPTIONAL:: LN

      INTEGER::KK,NN,MM,KV,NP,NM

      IF(B%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'IDEL2:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(A%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(A)
        CALL ALLOCATE(A)
      ENDIF

      IF(B%LN.NE.0.0D0) THEN
        WRITE(*,*) 'IDEL2:LOGTERM FOR INPUT OF IDEL2 NONZERO.'
        WRITE(*,*) 'LOGTERM=',B%LN
      ENDIF
      A = B

      !> M .EQ. 0
      NN=NZCHOPS(1)
      NP=NN+1
      NM=NN-1

      ALLOCATE( DEL2OP(NP,7,NXCHOPDIM) )

      DEL2OP(:NP,2:6,1)  = BAND_LEG_RAT_DEL2H(NP,M(1),ELL,TFM%NORM(:,1))

      DO KK=2,NXCHOPDIM
        DEL2OP(:NP,2:3,KK)=DEL2OP(:NP,2:3,1)
        DEL2OP(:NP,4  ,KK)=DEL2OP(:NP,4  ,1)-AK(1,KK)**2.
        DEL2OP(:NP,5:6,KK)=DEL2OP(:NP,5:6,1)
        DEL2OP(:NP,1  ,KK)=0
        DEL2OP(:NP,7  ,KK)=0
      ENDDO

      DEL2OP(1,4,1)=DEL2OP(1,4,1)+4.0D0/3/ELL2
      DEL2OP(2,3,1)=DEL2OP(2,3,1)-2*TFM%NORM(1,1)/TFM%NORM(2,1)/ELL2
      DEL2OP(3,2,1)=DEL2OP(3,2,1)+2.0D0/3*TFM%NORM(1,1)/TFM%NORM(3,1)/ELL2
      DEL2OP(2:NP,1:5,1)=DEL2OP(1:NP-1,2:6,1)
      DEL2OP(:,6:7,1)=0
      DEL2OP(1,:,1)=0
      DEL2OP(1,4,1)=1

      A%E(2:NN,1,1)=A%E(1:NM,1,1)

      IF(PRESENT(LN)) THEN
        A%E(1,1,1)=LN/TFM%NORM(1,1)
      ELSE
        A%E(1,1,1)=0
      ENDIF

      CALL LUB(DEL2OP(:NN,:,:),4)
      CALL SOLVEB(DEL2OP(:NN,:,:),4,A%E(:NN,1,:NXCHOPDIM))

      DEALLOCATE( DEL2OP )

      !> M .NE. 0
      ALLOCATE( DEL2OP(NZCHOP,5,NXCHOPDIM) )

      DO MM=2,NTCHOP
        NN=NZCHOPS(MM)
        DEL2OP(:NN,1:5,1) = BAND_LEG_RAT_DEL2H(NN,M(MM),ELL,TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          DEL2OP(:NN,1:2,KK)=DEL2OP(:NN,1:2,1)
          DEL2OP(:NN,3  ,KK)=DEL2OP(:NN,3  ,1)-AK(MM,KK)**2
          DEL2OP(:NN,4:5,KK)=DEL2OP(:NN,4:5,1)
        ENDDO
        CALL LUB(DEL2OP(:NN,1:5,1:NXCHOPDIM),3)
        CALL SOLVEB(DEL2OP(:NN,1:5,1:NXCHOPDIM),3,A%E(:NN,MM,1:NXCHOPDIM))
      ENDDO

      DEALLOCATE( DEL2OP )

      IF(PRESENT(LN)) THEN
        A%LN = A%E(1,1,1)*TFM%NORM(1,1)
      ELSE
        A%LN=0.0
      ENDIF

      CALL CHOPDO(A)

      DO KK=1,NXCHOPDIM
        A%E(1,1,KK)=A%E(1,1,KK)  &
        -SUM(TFM%AT1(:NZCHOP)*A%E(:NZCHOP,1,KK))/TFM%AT1(1)
      ENDDO

      RETURN
      END SUBROUTINE IDEL2
!=======================================================================
      SUBROUTINE IDEL2LN(B,A)
!=======================================================================
! [USAGE]: 
! INVERSE LAPLACIAN (DEL-SQUARE^(-1)) OPERATOR.
! LOG-TERM IS COMPUTED INTO A%LN TO SATISFY FIELD BEHAVIOR AT INFINITY.
! [PARAMETERS]:
! B >> INPUT IN FFF SPACE
! A >> ON EXIT, (DEL^2)^(-1)(B)
! [DEPENDENCIES]:
! 1. BAND_LEG_RAT_DEL2H(~) @ MOD_LIN_LEGENDRE
! 2. LUB(~) @ MOD_BANDMAT
! 3. SOLVEB(~) @ MOD_BANDMAT
! 4. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 5. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: A, B
      REAL(P8),DIMENSION(:,:,:),ALLOCATABLE :: DEL2OP
      INTEGER::KK,NN,MM,KV,NP,NM

      IF(B%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'IDEL2:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(A%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(A)
        CALL ALLOCATE(A)
      ENDIF

      IF(B%LN.NE.0.0D0) THEN
        WRITE(*,*) 'IDEL2:LOGTERM NONZERO.'
        WRITE(*,*) 'LOGTERM=',B%LN
      ENDIF
      A = B

      !> M .EQ. 0
      NN=NZCHOPS(1)
      NP=NN+1
      NM=NN-1

      ALLOCATE( DEL2OP(NP,5,NXCHOPDIM) )

      DEL2OP(:NP,1:5,1)  = BAND_LEG_RAT_DEL2H(NP,M(1),ELL,TFM%NORM(:,1))

      DO KK=2,NXCHOPDIM
        DEL2OP(:NP,1:2,KK)=DEL2OP(:NP,1:2,1)
        DEL2OP(:NP,3  ,KK)=DEL2OP(:NP,3  ,1)-AK(1,KK)**2
        DEL2OP(:NP,4:5,KK)=DEL2OP(:NP,4:5,1)
      ENDDO

      DEL2OP(1,3,1)=DEL2OP(1,3,1)+4.0D0/3/ELL2
      DEL2OP(2,2,1)=DEL2OP(2,2,1)-2*TFM%NORM(1,1)/TFM%NORM(2,1)/ELL2
      DEL2OP(3,1,1)=DEL2OP(3,1,1)+2.0D0/3*TFM%NORM(1,1)/TFM%NORM(3,1)/ELL2

      CALL LUB(DEL2OP(:NN,:,:),3)
      CALL SOLVEB(DEL2OP(:NN,:,:),3,A%E(:NN,1,:NXCHOPDIM))

      A%LN = A%E(1,1,1)*TFM%NORM(1,1)

      DEALLOCATE( DEL2OP )

      !> M .NE. 0
      ALLOCATE( DEL2OP(NZCHOP,5,NXCHOPDIM) )

      DO MM=2,NTCHOP
        NN=NZCHOPS(MM)
        DEL2OP(:NN,1:5,1) = BAND_LEG_RAT_DEL2H(NN,M(MM),ELL,TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          DEL2OP(:NN,1:2,KK)=DEL2OP(:NN,1:2,1)
          DEL2OP(:NN,3  ,KK)=DEL2OP(:NN,3  ,1)-AK(MM,KK)**2
          DEL2OP(:NN,4:5,KK)=DEL2OP(:NN,4:5,1)
        ENDDO
        CALL LUB(DEL2OP(:NN,1:5,1:NXCHOPDIM),3)
        CALL SOLVEB(DEL2OP(:NN,1:5,1:NXCHOPDIM),3,A%E(:NN,MM,1:NXCHOPDIM))
      ENDDO

      DEALLOCATE( DEL2OP )

      CALL CHOPDO(A)

      DO KK=1,NXCHOPDIM
        A%E(1,1,KK)=A%E(1,1,KK) &
        -SUM(TFM%AT1(:NZCHOP)*A%E(:NZCHOP,1,KK))/TFM%AT1(1)
      ENDDO

      RETURN
      END SUBROUTINE IDEL2LN
!=======================================================================
      SUBROUTINE HELM(A,B,ALP)
!=======================================================================
! [USAGE]: 
! HELMHOLTZ OPERATOR (DEL^2 - ALPHA)*.
! [PARAMETERS]:
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, (DEL^2 - ALPHA)(A)
! ALP >> ALPHA IN THE HELMHOLTZ OPERATOR
! [DEPENDENCIES]:
! 1. DEL2(~) @ MOD_LEGOPS
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: A, B
      REAL(P8):: ALP

      CALL DEL2(A,B)

      B%E=B%E-ALP*A%E
      B%LN = -ALP*A%LN

      RETURN
      END SUBROUTINE HELM
!=======================================================================
      SUBROUTINE IHELM(B,A,ALP)
!=======================================================================
! [USAGE]: 
! INVERSE HELMHOLTZ OPERATOR (DEL^2 - ALPHA)^(-1)*.
! [PARAMETERS]:
! B >> INPUT IN FFF SPACE
! A >> ON EXIT, (DEL^2 - ALPHA)^(-1)(A)
! ALP >> ALPHA IN THE HELMHOLTZ OPERATOR
! [DEPENDENCIES]:
! 1. BAND_LEG_RAT_DEL2H(~) @ MOD_LIN_LEGENDRE
! 2. LUB(~) @ MOD_BANDMAT
! 3. SOLVEB(~) @ MOD_BANDMAT
! 4. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 5. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: A, B
      REAL(P8),INTENT(IN):: ALP
      REAL(P8),DIMENSION(:,:,:),ALLOCATABLE :: HELMOP
      REAL(P8),DIMENSION(:,:),ALLOCATABLE :: XM
      INTEGER::KK,NN,MM,KV,NP,NM

      IF(B%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'IDEL2:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(A%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(A)
        CALL ALLOCATE(A)
      ENDIF

      IF(ALP.EQ.0.0) THEN
        WRITE(*,*) 'IHELM: ALP CAN''T BE ZERO.'
        STOP
      ENDIF

      A=B
      A%LN= -B%LN/ALP
      A%E(1,1,1)=A%E(1,1,1)-A%LN*(4.0D0/3/ELL2)/TFM%NORM(1,1)
      A%E(2,1,1)=A%E(2,1,1)-A%LN*(-2/TFM%NORM(2,1)/ELL2)
      A%E(3,1,1)=A%E(3,1,1)-A%LN*(+2.0D0/3/TFM%NORM(3,1)/ELL2)

      !> M .EQ. 0
      NN=NZCHOPS(1)
      NP=NN+1
      NM=NN-1

      ALLOCATE( XM(NP,3) )

      XM = BAND_LEG_XM(NP,M(1),TFM%NORM(:,1))

      ALLOCATE( HELMOP(NP,7,NXCHOPDIM) )

      HELMOP(:NP,2:6,1) = BAND_LEG_RAT_DEL2H(NP,M(1),ELL,TFM%NORM(:,1))
      HELMOP(:NP, 4 ,1) = HELMOP(:NP, 4 ,1)-ALP

      DO KK=2,NXCHOPDIM
        HELMOP(:NP,2:3,KK)= HELMOP(:NP,2:3,1)
        HELMOP(:NP,4  ,KK)= HELMOP(:NP,4  ,1)-AK(1,KK)**2
        HELMOP(:NP,5:6,KK)= HELMOP(:NP,5:6,1)
        HELMOP(:,:,KK)=BANMUL(HELMOP(:,2:6,KK),3,XM,2)
      ENDDO

      HELMOP(:,:,1)=BANMUL(HELMOP(:,2:6,1),3,XM,2)

      CALL LUB(HELMOP(:NM,:,:),4)
      CALL SOLVEB(HELMOP(:NM,:,:),4,A%E(:NM,1,:NXCHOPDIM))

      A%E(NN,1,:NXCHOPDIM)=0

      DO KK=1,NXCHOPDIM
        A%E(:NN,1,KK)=BANMUL(XM(:NN,:),2,A%E(1:NN,1,KK))
      ENDDO

      DEALLOCATE( XM )
      DEALLOCATE( HELMOP )

      !> M .NE. 0
      ALLOCATE( HELMOP(NZCHOP,5,NXCHOPDIM) )

      DO MM=2,NTCHOP
        NN=NZCHOPS(MM)
        HELMOP(:NN,1:5,1) = BAND_LEG_RAT_DEL2H(NN,M(MM),ELL,TFM%NORM(:,MM))
        HELMOP(:NN, 3 ,1) = HELMOP(:NN, 3 ,1)-ALP
        DO KK=1,NXCHOPDIM
          HELMOP(:NN,1:2,KK)=HELMOP(:NN,1:2,1)
          HELMOP(:NN,3  ,KK)=HELMOP(:NN,3  ,1)-AK(MM,KK)**2
          HELMOP(:NN,4:5,KK)=HELMOP(:NN,4:5,1)
        ENDDO
        CALL LUB(HELMOP(:NN,1:5,1:NXCHOPDIM),3)
        CALL SOLVEB(HELMOP(:NN,1:5,1:NXCHOPDIM),3,A%E(:NN,MM,1:NXCHOPDIM))
      ENDDO

      DEALLOCATE( HELMOP )

      CALL CHOPDO(A)

      RETURN
      END SUBROUTINE IHELM
!=======================================================================
      SUBROUTINE HELMP(P,A,B,ALP,NU)
!=======================================================================
! [USAGE]: 
! REPEATED HELMHOLTZ OPERATOR (DEL^P+NU*DEL^2+ALPHA)* 
! WHERE P >= 4 AND EVEN.
! [PARAMETERS]:
! P >> EVEN NUMBER AND >=4
! A >> INPUT IN FFF SPACE
! B >> ON EXIT, (DEL^P+NU*DEL^2+ALPHA)(A)
! ALP >> ALPHA IN THE HELMHOLTZ OPERATOR
! NU >> (OPTIONAL) DEFAULT IS 0.
! [DEPENDENCIES]:
! 1. DEL2(~) @ MOD_LEGOPS
! 2. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER:: P
      TYPE(SCALAR):: A,B
      REAL(P8),INTENT(IN):: ALP
      REAL(P8),OPTIONAL:: NU

      TYPE(SCALAR):: W,D2
      INTEGER:: I

      CALL ALLOCATE(W)
      CALL ALLOCATE(D2)
      CALL DEL2(A,D2)

      B=D2

      DO I=4,P,2
        CALL DEL2(B,W)
        B=W
      ENDDO

      B%E = B%E -ALP*A%E

      IF(PRESENT(NU)) THEN
        B%E=B%E +NU*D2%E
      ENDIF

      B%LN = -ALP*A%LN

      CALL DEALLOCATE(W)
      CALL DEALLOCATE(D2)

      RETURN
      END SUBROUTINE HELMP
!=======================================================================

      SUBROUTINE IHELMP(P,B,A,ALP,NUIN)
!=======================================================================
! [USAGE]: 
! INVERSE REPEATED HELMHOLTZ OPERATOR (DEL^P+NUIN*DEL^2-ALPHA)^(-1)* 
! WHERE P >= 4 AND EVEN.
! [PARAMETERS]:
! P >> EVEN NUMBER AND >=4
! B >> INPUT IN FFF SPACE
! A >> ON EXIT, (DEL^P+NU*DEL^2+ALPHA)^(-1)(B)
! ALP >> ALPHA IN THE HELMHOLTZ OPERATOR
! NUIN >> (OPTIONAL) DEFAULT IS 0.
! [DEPENDENCIES]:
! 1. BAND_LEG_RAT_DEL2H(~) @ MOD_LIN_LEGENDRE
! 2. LUB(~) @ MOD_BANDMAT
! 3. SOLVEB(~) @ MOD_BANDMAT
! 4. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 5. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
        IMPLICIT NONE
        INTEGER,INTENT(IN):: P
        TYPE(SCALAR):: A, B
        REAL(P8),INTENT(IN):: ALP
        REAL(P8),INTENT(IN),OPTIONAL:: NUIN
  
        REAL(P8),DIMENSION(:,:,:),ALLOCATABLE :: HELMP_VAL
        REAL(P8),DIMENSION(:,:),ALLOCATABLE :: DEL2H,DEL2_VAL
        REAL(P8),DIMENSION(:,:),ALLOCATABLE :: DEL4H,DEL6H
        REAL(P8),DIMENSION(:,:),ALLOCATABLE :: DEL8H,DEL10H
        REAL(P8),DIMENSION(:,:),ALLOCATABLE :: DEL12H,DEL14H
        REAL(P8),DIMENSION(NZCHOP):: AL2,AL
  
        INTEGER::KK,NN,MM,KV,NP,NM,I
        REAL(P8):: NU,ELL4,K2
  
        IF(B%SPACE.NE. FFF_SPACE) THEN
          WRITE(*,*) 'IDEL2:INPUT NOT IN FFF_SPACE.'
          STOP
        ENDIF
  
        IF(A%SPACE.NE.FFF_SPACE) THEN
          CALL DEALLOCATE(A)
          CALL ALLOCATE(A)
        ENDIF
  
        IF(M(1).NE.0) THEN
          WRITE(*,*) 'IHELMP: M(1) MUST BE ZERO.'
          STOP
        ENDIF
  
        IF(ALP.EQ.0.0) THEN
          WRITE(*,*) 'IHELMP: ALP CANNOT BE ZERO.'
          STOP
        ENDIF
  
        IF(PRESENT(NUIN)) THEN
          NU=NUIN
        ELSE
          NU=0.0D0
        ENDIF
        A=B
  
        !> SUBTRACT LOGARITHMIC TERM CONTRIBUTION
        A%LN= -B%LN/ALP
        ELL4=ELL2**2
        AL2=0
        AL2(1)=A%LN*(4.0D0/3/ELL2)/TFM%NORM(1,1)
        AL2(2)=A%LN*(-2/TFM%NORM(2,1)/ELL2)
        AL2(3)=A%LN*(+2.0D0/3/TFM%NORM(3,1)/ELL2)
        AL=AL2
  
        ALLOCATE( DEL2H(NZCHOP,5) )
  
        DEL2H = BAND_LEG_RAT_DEL2H(NZCHOP,M(1),ELL,TFM%NORM(:,1))
  
        DO I=4,P,2
          AL=BANMUL(DEL2H,3,AL)
        ENDDO
  
        A%E(:NZCHOP,1,1)=A%E(:NZCHOP,1,1)-AL-NU*AL2
  
        DEALLOCATE( DEL2H )
  
        !> ALL M
        ALLOCATE( DEL2_VAL(NZCHOP+P,5) )
        ALLOCATE( HELMP_VAL(NZCHOP+P,2*P+1,NXCHOPDIM) )
        ALLOCATE( DEL2H(NZCHOP+P,5) )
        ALLOCATE( DEL4H(NZCHOP+P,9) )
        IF(P >=6 ) ALLOCATE( DEL6H(NZCHOP+P,13) )
        IF(P >=8 ) ALLOCATE( DEL8H(NZCHOP+P,17) )
        IF(P >=10) ALLOCATE( DEL10H(NZCHOP+P,21) )
        IF(P >=12) ALLOCATE( DEL12H(NZCHOP+P,25) )
        IF(P >=14) ALLOCATE( DEL14H(NZCHOP+P,29) )
  
        DO MM=1,NTCHOP
          NN=NZCHOPS(MM)
          NP=NN+P
          DEL2H(:NP,:) = BAND_LEG_RAT_DEL2H(NP,M(MM),ELL,TFM%NORM(:,MM))
          DEL4H(:NP,:) = BANMUL(DEL2H(:NP,:),3,DEL2H(:NP,:),3)
          IF(P >=6 ) DEL6H(:NP,:) =BANMUL(DEL4H(:NP,:),  5,DEL2H(:NP,:),3)
          IF(P >=8 ) DEL8H(:NP,:) =BANMUL(DEL6H(:NP,:),  7,DEL2H(:NP,:),3)
          IF(P >=10) DEL10H(:NP,:)=BANMUL(DEL8H(:NP,:),  9,DEL2H(:NP,:),3)
          IF(P >=12) DEL12H(:NP,:)=BANMUL(DEL10H(:NP,:),11,DEL2H(:NP,:),3)
          IF(P >=14) DEL14H(:NP,:)=BANMUL(DEL12H(:NP,:),13,DEL2H(:NP,:),3)
          DO KK=1,NXCHOPDIM
            K2 = AK(MM,KK)**2
            DEL2_VAL(:NP,:)=DEL2H(:NP,:)
            DEL2_VAL(:NP,3)=DEL2_VAL(:NP,3)-K2
            IF(P.EQ.4) THEN
              HELMP_VAL(:NP, : ,KK)=                          DEL4H(:NP,:)
              HELMP_VAL(:NP,3:7,KK)=HELMP_VAL(:NP,3:7,KK)-2*K2   *DEL2H(:NP,:)
              HELMP_VAL(:NP,  5,KK)=HELMP_VAL(:NP,5,KK)  +  K2**2
            ELSE IF(P.EQ.6) THEN
              HELMP_VAL(:NP, :  ,KK)=                           DEL6H(:NP,:)
              HELMP_VAL(:NP,3:11,KK)=HELMP_VAL(:NP,3:11,KK)-3*K2   *DEL4H(:NP,:)
              HELMP_VAL(:NP,5: 9,KK)=HELMP_VAL(:NP,5: 9,KK)+3*K2**2*DEL2H(:NP,:)
              HELMP_VAL(:NP,   7,KK)=HELMP_VAL(:NP,   7,KK)-  K2**3
            ELSE IF(P.EQ.8) THEN
              HELMP_VAL(:NP, : ,KK)   =                         DEL8H(:NP,:)
              HELMP_VAL(:NP,3:15,KK)=HELMP_VAL(:NP,3:15,KK)-4*K2   *DEL6H(:NP,:)
              HELMP_VAL(:NP,5:13,KK)=HELMP_VAL(:NP,5:13,KK)+6*K2**2*DEL4H(:NP,:)
              HELMP_VAL(:NP,7:11,KK)=HELMP_VAL(:NP,7:11,KK)-4*K2**3*DEL2H(:NP,:)
              HELMP_VAL(:NP,   9,KK)=HELMP_VAL(:NP,   9,KK)+  K2**4
            ELSE IF(P.EQ.10) THEN
              HELMP_VAL(:NP,:,KK)   =                           DEL10H(:NP,:)
              HELMP_VAL(:NP,3:19,KK)=HELMP_VAL(:NP,3:19,KK)- 5*K2   *DEL8H(:NP,:)
              HELMP_VAL(:NP,5:17,KK)=HELMP_VAL(:NP,5:17,KK)+10*K2**2*DEL6H(:NP,:)
              HELMP_VAL(:NP,7:15,KK)=HELMP_VAL(:NP,7:15,KK)-10*K2**3*DEL4H(:NP,:)
              HELMP_VAL(:NP,9:13,KK)=HELMP_VAL(:NP,9:13,KK)+ 5*K2**4*DEL2H(:NP,:)
              HELMP_VAL(:NP,  11,KK)=HELMP_VAL(:NP,  11,KK)-   K2**5
            ELSE IF(P.EQ.12) THEN
              HELMP_VAL(:NP,:,KK)   =                            DEL12H(:NP,:)
              HELMP_VAL(:NP,3:23,KK)=HELMP_VAL(:NP, 3:23,KK)- 6*K2  *DEL10H(:NP,:)
              HELMP_VAL(:NP,5:21,KK)=HELMP_VAL(:NP, 5:21,KK)+15*K2**2*DEL8H(:NP,:)
              HELMP_VAL(:NP,7:19,KK)=HELMP_VAL(:NP, 7:19,KK)-20*K2**3*DEL6H(:NP,:)
              HELMP_VAL(:NP,9:17,KK)=HELMP_VAL(:NP, 9:17,KK)+15*K2**4*DEL4H(:NP,:)
              HELMP_VAL(:NP,11:15,KK)=HELMP_VAL(:NP,11:15,KK)-6*K2**5*DEL2H(:NP,:)
              HELMP_VAL(:NP,   13,KK)=HELMP_VAL(:NP,  13,KK)+    K2**6
            ELSE IF(P.EQ.14) THEN
              HELMP_VAL(:NP,:,KK)   =                            DEL14H(:NP,:)
              HELMP_VAL(:NP,3:27,KK)=HELMP_VAL(:NP,3:27,KK)- 7*K2   *DEL12H(:NP,:)
              HELMP_VAL(:NP,5:25,KK)=HELMP_VAL(:NP,5:25,KK)+21*K2**2*DEL10H(:NP,:)
              HELMP_VAL(:NP, 7:23,KK)=HELMP_VAL(:NP,7:23,KK)-35*K2**3*DEL8H(:NP,:)
              HELMP_VAL(:NP, 9:21,KK)=HELMP_VAL(:NP,9:21,KK)+35*K2**4*DEL6H(:NP,:)
              HELMP_VAL(:NP,11:19,KK)=HELMP_VAL(:NP,11:19,KK)-21*K2**5&
                                                             *DEL4H(:NP,:)
              HELMP_VAL(:NP,13:17,KK)=HELMP_VAL(:NP,13:17,KK)+ 7*K2**6&
                                                             *DEL2H(:NP,:)
              HELMP_VAL(:NP,   15,KK)=HELMP_VAL(:NP,   15,KK)-   K2**7
            ELSE 
              WRITE(*,*) 'P OUT OF RANGE. P=',P
              STOP
            ENDIF
            HELMP_VAL(:NN,P-1:P+3,KK)=HELMP_VAL(:NN,P-1:P+3,KK)+NU*DEL2_VAL(:NN,:)
            HELMP_VAL(:NN,P+1,KK)=HELMP_VAL(:NN,P+1,KK)-ALP
          ENDDO
          CALL LUB(HELMP_VAL(:NN,:,:),P+1)
          CALL SOLVEB(HELMP_VAL(:NN,:,:),P+1,A%E(:NN,MM,1:NXCHOPDIM))
        ENDDO
  
        DEALLOCATE( DEL2H, DEL4H, DEL2_VAL, HELMP_VAL )
        IF(P >=6 ) DEALLOCATE( DEL6H )
        IF(P >=8 ) DEALLOCATE( DEL8H )
        IF(P >=10) DEALLOCATE( DEL10H )
        IF(P >=12) DEALLOCATE( DEL12H )
        IF(P >=14) DEALLOCATE( DEL14H )
  
        !> ZERO AT INFINITY
        A%E(1,1,:)=A%E(1,1,:)-CALCAT1(A)/TFM%PF(1,1,1)
  
        CALL CHOPDO(A)
  
        RETURN
        END SUBROUTINE IHELMP

!=======================================================================

      SUBROUTINE IHELMPO(P,B,A,ALP,NUIN)
!=======================================================================
! [USAGE]: 
! INVERSE REPEATED HELMHOLTZ OPERATOR (DEL^P+NUIN*DEL^2+ALPHA)^(-1)* 
! WHERE P >= 4 AND EVEN. SHORTER VERSION OF IHEMLP BUT INCLUDING THE
! TRIPLE NESTED LOOP (WILL BE SLOWER)
! [PARAMETERS]:
! P >> EVEN NUMBER AND >=4
! B >> INPUT IN FFF SPACE
! A >> ON EXIT, (DEL^P+NU*DEL^2+ALPHA)^(-1)(B)
! ALP >> ALPHA IN THE HELMHOLTZ OPERATOR
! NUIN >> (OPTIONAL) DEFAULT IS 0.
! [DEPENDENCIES]:
! 1. BAND_LEG_RAT_DEL2H(~) @ MOD_LIN_LEGENDRE
! 2. LUB(~) @ MOD_BANDMAT
! 3. SOLVEB(~) @ MOD_BANDMAT
! 4. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 5. CHOPDO(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER,INTENT(IN):: P
      TYPE(SCALAR):: A, B
      REAL(P8),INTENT(IN):: ALP
      REAL(P8),INTENT(IN),OPTIONAL:: NUIN

      REAL(P8),DIMENSION(:,:,:),ALLOCATABLE :: HELMP_VAL
      REAL(P8),DIMENSION(:,:),ALLOCATABLE :: DEL2H,DEL2_VAL
      REAL(P8),DIMENSION(NZCHOP):: AL2,AL


      INTEGER::KK,NN,MM,KV,NP,NM,I
      REAL(P8):: NU,ELL4

      IF(B%SPACE.NE. FFF_SPACE) THEN
        WRITE(*,*) 'IDEL2:INPUT NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(A%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(A)
        CALL ALLOCATE(A)
      ENDIF

      IF(M(1).NE.0) THEN
        WRITE(*,*) 'IHELMP: M(1) MUST BE ZERO.'
        STOP
      ENDIF

      IF(ALP.EQ.0.0) THEN
        WRITE(*,*) 'IHELMP: ALP CAN''T BE ZERO.'
        STOP
      ENDIF

      IF(PRESENT(NUIN)) THEN
        NU=NUIN
      ELSE
        NU=0.0D0
      ENDIF
      A=B

      !> SUBTRACT LOGARITHMIC TERM CONTRIBUTION
      A%LN= -B%LN/ALP
      ELL4=ELL2**2
      AL2=0
      AL2(1)=A%LN*(4.0D0/3/ELL2)/TFM%NORM(1,1)
      AL2(2)=A%LN*(-2/TFM%NORM(2,1)/ELL2)
      AL2(3)=A%LN*(+2.0D0/3/TFM%NORM(3,1)/ELL2)
      AL=AL2

      ALLOCATE( DEL2H(NZCHOP,5) )

      DEL2H = BAND_LEG_RAT_DEL2H(NZCHOP,M(1),ELL,TFM%NORM(:,1))

      DO I=4,P,2
        AL=BANMUL(DEL2H,3,AL)
      ENDDO

      A%E(:NZCHOP,1,1)=A%E(:NZCHOP,1,1)-AL-NU*AL2

      DEALLOCATE( DEL2H )

      !> ALL M
      ALLOCATE( DEL2H(NZCHOP+P,5) )
      ALLOCATE( DEL2_VAL(NZCHOP+P,5) )
      ALLOCATE( HELMP_VAL(NZCHOP+P,2*P+1,NXCHOPDIM) )

      DO MM=1,NTCHOP
        NN=NZCHOPS(MM)
        NP=NN+P
        DEL2H(:NP,:) = BAND_LEG_RAT_DEL2H(NP,M(MM),ELL,TFM%NORM(:,MM))
        DO KK=1,NXCHOPDIM
          DEL2_VAL(:NP,:)=DEL2H(:NP,:)
          DEL2_VAL(:NP,3)=DEL2_VAL(:NP,3)-AK(MM,KK)**2
          HELMP_VAL(:NP,P-1:P+3,KK)=DEL2_VAL(:NP,:)
          DO I=4,P,2
            HELMP_VAL(:NP,P+1-I:P+1+I,KK)= &
            BANMUL(HELMP_VAL(:NP,P+3-I:P-1+I,KK),I-1, DEL2_VAL(:NP,:),3 )
          ENDDO
          HELMP_VAL(:NN,P-1:P+3,KK)=HELMP_VAL(:NN,P-1:P+3,KK)+NU*DEL2_VAL(:NN,:)
          HELMP_VAL(:NN,P+1,KK)=HELMP_VAL(:NN,P+1,KK)-ALP
        ENDDO
        CALL LUB(HELMP_VAL(:NN,:,:),P+1)
        CALL SOLVEB(HELMP_VAL(:NN,:,:),P+1,A%E(:NN,MM,1:NXCHOPDIM))
      ENDDO

      DEALLOCATE( DEL2H, DEL2_VAL, HELMP_VAL )

      !> ZERO AT INFINITY
      A%E(1,1,:)=A%E(1,1,:)-CALCAT1(A)/TFM%PF(1,1,1)

      CALL CHOPDO(A)

      RETURN
      END SUBROUTINE IHELMPO

!=======================================================================
      SUBROUTINE PC2VEL(PSI,CHI,RUR,RUP,UZ,C)
!=======================================================================
! [USAGE]: 
! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VELOCITY VARIALBES
! RETURNS (R*UR,R*UP,UZ) IF C IS NOT PRESENT
! RETURNS (R*UR,R*UP,UZ/(1-MU)^2 IF C IS PRESENT
! [PARAMETERS]:
! PSI >> TOROIDAL TERM IN A SCALAR-TYPE VARIABLE
! CHI >> POLOIDAL TERM IN A SCALAR-TYPE VARIABLE
! RUR >> ON EXIT, R*U_R (RADIAL VEL.) TERM IN A SCALAR-TYPE VARIABLE
! RUP >> ON EXIT, R*U_P (AZIMUTHAL VEL.) TERM IN A SCALAR-TYPE VARIABLE
! UZ >> ON EXIT, U_Z (AXIAL VEL.) TERM IN A SCALAR-TYPE VARIABLE
! C >> (OPTIONAL) IF NOT EXIST, RETURNS (RUR,RUP,UZ)
!                 IF EXIST, RETURNS (RUR,RUP,UZ/(1-X)**2)
! [NOTE]:
! FOR A SOLENOIDAL VECTOR FIELD (DEL.V=0), THE FIELD CAN BE DECOMPOSED
! BY V = DEL X (PSI E_Z) + DEL X (DEL X CHI E_Z)
! WHERE PSI AND CHI ARE CALLED TOROIDAL AND POLOIDAL TERMS, RSPCTVLY.
! THE FOLLOWING RELATIONS HOLD TRUE.
! UR = 1/R*D/DPHI(PSI) + D/DR(D/DZ(CHI))
! UPHI = -D/DR(PSI) = 1/R*D/DPHI(D/DZ(CHI))
! UZ = -DEL^2_H(CHI)
! [DEPENDENCIES]:
! 1. XXDX(~) @ MOD_LEGOPS
! 2. DELSQH(~) @ MOD_LEGOPS
! 3. MULXM(~) @ MOD_BANDMAT
! 4. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: PSI,CHI,RUR,RUP,UZ
      INTEGER:: MM,KK,NN,MV
      REAL(P8):: KV
      TYPE(SCALAR):: W
      CHARACTER(LEN=1),OPTIONAL:: C

      CALL XXDX(CHI,RUR)                                                 ! RUR NOW HAS (1-X^2)D/DX(CHI) = R*D/DR(CHI)
      CALL XXDX(PSI,RUP)                                                 ! RUP NOW HAS (1-X^2)D/DX(PSI) = R*D/DR(PSI)
      
      DO MM=1,NTCHOP
        NN=NZCHOPS(MM)
        MV=M(MM)
        DO KK=1,NXCHOPDIM
          KV=AK(MM,KK)
          RUR%E(:NN,MM,KK)=IU*MV*PSI%E(:NN,MM,KK)&                       ! RUR = D/DPHI(PSI)+ R*D/DR(D/DZ(CHI))
                           +(IU*KV)*RUR%E(:NN,MM,KK)
          RUP%E(:NN,MM,KK)=     -RUP%E(:NN,MM,KK)&                       ! RUP =-R*D/DR(PSI)+ D/DPHI(D/DZ(CHI))
                           -(MV*KV)*CHI%E(:NN,MM,KK)                     ! (IU*MV)*(IU*KV) = -(MV*KV)
        ENDDO
      ENDDO

      CALL DELSQH(CHI,UZ)                                                ! UZ/(1-X)**2 = - DELSQH(CHI) (MINUS SIGN ASSIGNED AT THE END)
      IF(.NOT.PRESENT(C)) THEN
        CALL ALLOCATE(W)
        CALL MULXM(UZ,W)                                                 ! W <- (1-X)*[UZ/(1-X)**2] = UZ/(1-X)
        CALL MULXM(W,UZ)                                                 ! UZ <- (1-X)*[W] = (1-X)*[UZ/(1-X)]= UZ
        CALL DEALLOCATE(W)
      ENDIF
      UZ%E = -UZ%E

      RETURN
      END SUBROUTINE PC2VEL
!=======================================================================
      SUBROUTINE PC2VOR(PSI,CHI,ROR,ROP,OZ,C)
!=======================================================================
! [USAGE]: 
! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VELOCITY VARIALBES
! RETURNS (R*OR,R*OP,OZ) IF C IS NOT PRESENT
! RETURNS (R*OR,R*OP,OZ/(1-MU)^2 IF C IS PRESENT
! [PARAMETERS]:
! PSI >> TOROIDAL TERM IN A SCALAR-TYPE VARIABLE
! CHI >> POLOIDAL TERM IN A SCALAR-TYPE VARIABLE
! ROR >> ON EXIT, R*OMEGA_R TERM IN A SCALAR-TYPE VARIABLE
! ROP >> ON EXIT, R*OMEGA_P TERM IN A SCALAR-TYPE VARIABLE
! OZ >> ON EXIT, OMEGA_Z TERM IN A SCALAR-TYPE VARIABLE
! C >> (OPTIONAL) IF NOT EXIST, RETURNS (ROR,ROP,OZ)
!                 IF EXIST, RETURNS (ROR,ROP,OZ/(1-X)**2)
! [NOTE]:
! SUBROUTINE STRUCTURE IS ALMOST THE SAME AS PC2VEL EXCEPT THAT
! (-DEL^2(CHI)) AND PSI IS SUBSTITUTED FOR PSI AND CHI @ PC2VEL
! BECAUSE VELOCITY  = DEL X (PSI         E_Z) + DEL X (DEL X CHI E_Z)
! WHEREAS VORTICITY = DEL X (-DEL^2(CHI) E_Z) + DEL X (DEL X PSI E_Z). 
! [DEPENDENCIES]:
! 1. XXDX(~) @ MOD_LEGOPS
! 2. DELSQH(~) @ MOD_LEGOPS
! 3. MULXM(~) @ MOD_BANDMAT
! 4. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: PSI,CHI,ROR,ROP,OZ
      INTEGER:: MM,KK,NN,MV
      REAL(P8):: KV
      TYPE(SCALAR):: W
      CHARACTER(LEN=1),OPTIONAL:: C

      CALL DEL2(CHI,OZ)                                                  ! OZ NOW HAS (DEL^2)(CHI)
      CALL XXDX(PSI,ROR)                                                 ! ROR NOW HAS (1-X^2)D/DX(PSI) = R*D/DR(PSI)
      CALL XXDX(OZ,ROP)                                                  ! ROP NOW HAS (1-X^2)D/DX((DEL^2)(CHI)) = R*D/DR((DEL^2)(CHI))

      ! if (psi%e(2,2,2).eq.1) then
      !   ror%e(:,2,2) = 99
      !   call mcat(ROR%e(:,2,2))
      ! endif

      DO MM=1,NTCHOP
        NN=NZCHOPS(MM)
        MV=M(MM)
        DO KK=1,NXCHOPDIM
          KV=AK(MM,KK)
          ! if (psi%e(2,2,2).eq.1) then
          !   if (mm.le.1) then
          !   WRITE(*,*) ''
          !   WRITE(*,*) mm,kk
          !   WRITE(*,*) NN
          !   ! if ((mm.eq.2).and.(kk.eq.2)) then
          !     ! write(*,*) iu
          !     ! write(*,*) kv
          !     ! WRITE(*,*) ''
          !     ! WRITE(*,*) 'INSDIE DO LOOP'
          !     call mcat(ROR%E(:,2,2))
          !   ! endif
          !   endif
          ! endif          
          ROR%E(:NN,MM,KK)=-IU*MV*OZ%E(:NN,MM,KK)&                       ! ROR = -D/DPHI(DEL^2(CHI))+ R*D/DR(D/DZ(PSI))
                           +IU*KV*ROR%E(:NN,MM,KK)

          ROP%E(:NN,MM,KK)=       ROP%E(:NN,MM,KK)&                      ! ROP = R*D/DR(DEL^2(CHI))+ D/DPHI(D/DZ(PSI))
                           -(KV*MV)*PSI%E(:NN,MM,KK)                     ! (IU*MV)*(IU*KV) = -(MV*KV)
        ENDDO
      ENDDO
      !if (psi%e(2,2,2).eq.1) call mcat(ror%e(:,2,2)) 

      CALL DELSQH(PSI,OZ)                                                ! OZ/(1-X)**2 = - DELSQH(PSI) (MINUS SIGN ASSIGNED AT THE END)

      IF(.NOT.PRESENT(C)) THEN
        CALL ALLOCATE(W)
        CALL MULXM(OZ,W)                                                 ! W <- (1-X)*[OZ/(1-X)**2] = OZ/(1-X)
        CALL MULXM(W,OZ)                                                 ! OZ <- (1-X)*[W] = (1-X)*[OZ/(1-X)]= OZ
        CALL DEALLOCATE(W)
      ENDIF
      OZ%E = -OZ%E

      RETURN
      END SUBROUTINE PC2VOR
!=======================================================================
      SUBROUTINE VPROD(RUR,RUP,UZ,ROR,ROP,OZ)
!=======================================================================
! [USAGE]: 
! WRAPPER OF VPRODSUB(~) COMPUTING A VECTOR PRODUCT (PPP) 
! (PR,PP,PZ) = (UR,UP,UZ) X (OR,OP,OZ)
! [PARAMETERS]:
! RUR >> ON ENTRY, R*U_R TERM FROM PC2VEL
! RUP >> ON ENTRY, R*U_P TERM FROM PC2VEL
! UZ >> ON ENTRY, U_Z TERM FROM PC2VEL
! ROR >> ON ENTRY, R*OMEGA_R TERM FROM PC2VOR
!        ON EXIT, R*[THE FIRST COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)]
!                 IN PHYSICAL-PHYSICAL-PHYSICAL SPACE
! ROP >> ON ENTRY, R*OMEGA_P TERM FROM PC2VOR
!        ON EXIT, R*[THE SECOND COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)]
!                 IN PHYSICAL-PHYSICAL-PHYSICAL SPACE
! OZ >> ON ENTRY, OMEGA_Z TERM FROM PC2OZ
!        ON EXIT, THE THIRD COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)
!                 IN PHYSICAL-PHYSICAL-PHYSICAL SPACE
! [DEPENDENCIES]:
! 1. TOFP(~) @ MOD_SCALAR3
! 2. VPRODSUB(~) @ MOD_LEGOPS
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      TYPE(SCALAR):: RUR,RUP,UZ
      TYPE(SCALAR):: ROR,ROP,OZ
      INTEGER:: NI,NJ,NK

      CALL TOFP(RUR)                                                     ! FFF -> PPP SPACE
      CALL TOFP(RUP)                                                     ! FFF -> PPP SPACE
      CALL TOFP(UZ)                                                      ! FFF -> PPP SPACE
      CALL TOFP(ROR)                                                     ! FFF -> PPP SPACE
      CALL TOFP(ROP)                                                     ! FFF -> PPP SPACE
      CALL TOFP(OZ)                                                      ! FFF -> PPP SPACE

      NI=SIZE(RUR%E,1)
      NJ=SIZE(RUR%E,2)
      NK=SIZE(RUR%E,3)

      CALL VPRODSUB(RUR%E,RUP%E,UZ%E,ROR%E,ROP%E,OZ%E,NI,NJ,NK)
      ! CALL MCAT(UZ%E(:,1,:))
      ! WRITE(*,*) SHAPE(UZ%e)

      RETURN
      END SUBROUTINE VPROD
!=======================================================================
      SUBROUTINE VPRODSUB(RUR,RUP,UZ,ROR,ROP,OZ,NI,NJ,NK)
!=======================================================================
! [USAGE]: 
! COMPUTE A VECTOR PRODUCT (PR,PP,PZ) = (UR,UP,UZ) X (OR,OP,OZ)
! THE RESULTS ARE THEN STORED INTO ROR, ROP AND OZ.
! [PARAMETERS]:
! RUR >> ON ENTRY, R*U_R IN PPP SPACE
! RUP >> ON ENTRY, R*U_P IN PPP SPACE
! UZ >> ON ENTRY, U_Z IN PPP SPACE
! ROR >> ON ENTRY, R*OMEGA_R IN PPP SPACE
!        ON EXIT, R*[THE FIRST COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)]
!                 IN PHYSICAL-PHYSICAL-PHYSICAL SPACE
! ROP >> ON ENTRY, R*OMEGA_P IN PPP SPACE
!        ON EXIT, R*[THE SECOND COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)]
!                 IN PHYSICAL-PHYSICAL-PHYSICAL SPACE
! OZ >> ON ENTRY, OMEGA_Z IN PPP SPACE
!        ON EXIT, THE THIRD COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)
!                 IN PHYSICAL-PHYSICAL-PHYSICAL SPACE
! NI >> 1ST DIMENSION SIZE (IN THE X(OR R) DIRECTION)
! NJ >> 2ND DIMENSION SIZE (IN THE AZIMUTHAL DIRECTION)
! NK >> 3RD DIMENSION SIZE (IN THE AXIAL DIRECTION)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER:: NI,NJ,NK
      COMPLEX(P8),DIMENSION(NI,NJ,NK):: RUR,RUP,UZ,ROR,ROP,OZ

      REAL(P8):: A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3
      INTEGER:: KK,MM,NN
      REAL(P8):: TH1,TH2,FUR1,FUP1,FUR2,FUP2,FSR1,FSP1,FSR2,FSP2,R
      COMPLEX(P8):: CM1,CM2
      REAL(P8):: AAA,BBB
      
      !> FREESTREAM & STRAIN
      UZ(:NZ,:NTH,:NX) = UZ(:NZ,:NTH,:NX)+CMPLX(NADD%UZ,NADD%UZ)
      DO MM=1,NTH
        TH1=TFM%THR(MM)
        TH2=TFM%THI(MM)
        FUR1 =  NADD%U*COS(TH1-NADD%ANG)
        FUP1 = -NADD%U*SIN(TH1-NADD%ANG)
        FUR2 =  NADD%U*COS(TH2-NADD%ANG)
        FUP2 = -NADD%U*SIN(TH2-NADD%ANG)
        FSR1 =  NADD%STRAIN*SIN(2*TH1)
        FSP1 =  NADD%STRAIN*COS(2*TH1)
        FSR2 =  NADD%STRAIN*SIN(2*TH2)
        FSP2 =  NADD%STRAIN*COS(2*TH2)
        DO KK=1,NX
          DO NN=1,NZ
            R = TFM%R(NN)
            CM1 = CMPLX(FUR1,FUR2)*R+CMPLX(FSR1,FSR2)*R**2
            CM2 = CMPLX(FUP1,FUP2)*R+CMPLX(FSP1,FSP2)*R**2
            RUR(NN,MM,KK)= RUR(NN,MM,KK)+CM1
            RUP(NN,MM,KK)= RUP(NN,MM,KK)+CM2
          ENDDO
        ENDDO
      ENDDO

      !> MAXIMUM VELOCITY MONITOR
      IF(VELMON%N.EQ.0) THEN
        VELMON%URMAX=0
        VELMON%UPMAX=0
        VELMON%UZMAX=0
        DO NN=1,NZ
          AAA = MAXVAL(REAL(RUR(NN,:NTH,:NX)))/TFM%R(NN)
          BBB = MINVAL(REAL(RUR(NN,:NTH,:NX)))/TFM%R(NN)
          VELMON%URMAX=MAX(ABS(AAA),ABS(BBB),VELMON%URMAX)
          AAA = MAXVAL(REAL(RUP(NN,:NTH,:NX)))/TFM%R(NN)
          BBB = MINVAL(REAL(RUP(NN,:NTH,:NX)))/TFM%R(NN)
          VELMON%UPMAX=MAX(ABS(AAA),ABS(BBB),VELMON%UPMAX)
          AAA = MAXVAL(REAL(UZ(NN,:NTH,:NX)))
          BBB = MINVAL(REAL(UZ(NN,:NTH,:NX)))
          VELMON%UZMAX=MAX(ABS(AAA),ABS(BBB),VELMON%UZMAX)
        ENDDO
        WRITE(6,60) VELMON%URMAX, VELMON%UPMAX, VELMON%UZMAX
 60     FORMAT('NONLIN: MAXVEL=',3E15.7)
      ENDIF
      VELMON%N = MOD(VELMON%N+1,VELMON%INT)
      
      !> CROSS PRODUCT
      DO KK=1,NX
         DO MM=1,NTH
            DO NN=1,NZ
               A1=REAL(RUR(NN,MM,KK))
               A2=REAL(RUP(NN,MM,KK))
               A3=REAL(UZ(NN,MM,KK))
               C1=AIMAG(RUR(NN,MM,KK))
               C2=AIMAG(RUP(NN,MM,KK))
               C3=AIMAG(UZ(NN,MM,KK))
               B1=REAL(ROR(NN,MM,KK))
               B2=REAL(ROP(NN,MM,KK))
               B3=REAL(OZ(NN,MM,KK))
               D1=AIMAG(ROR(NN,MM,KK))
               D2=AIMAG(ROP(NN,MM,KK))
               D3=AIMAG(OZ(NN,MM,KK))
               ROR(NN,MM,KK)=CMPLX(A2*B3-A3*B2, C2*D3-C3*D2,P8)
               ROP(NN,MM,KK)=CMPLX(A3*B1-A1*B3, C3*D1-C1*D3,P8)
               OZ(NN,MM,KK)=CMPLX(A1*B2-A2*B1, C1*D2-C2*D1,P8)&
                            /TFM%R(NN)**2.
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE VPRODSUB
!=======================================================================
      SUBROUTINE PROJECT(RUR,RUP,UZ,PSI,CHI)
!=======================================================================
! [USAGE]: 
! PROJECT A VECTOR (RUR,RUP,UZ) IN PPP SPACE
! TO A POLOIDAL-TOROIDAL FIELD (PSI, DEL^2(CHI))
! REFER TO EQS. (64,65,68,69 AND 70) IN MATSUSHIMA AND MARCUS (1997).
! THE RESULTS ARE THEN STORED INTO ROR, ROP AND OZ.
! [PARAMETERS]:
! RUR >> R*U_R IN PPP SPACE
! RUP >> R*U_P IN PPP SPACE
! UZ >> U_Z IN PPP SPACE
! PSI >> TOROIDAL TERM OF THE VELOCITY FIELD IN FFF SPACE
! CHI >> POLOIDAL TERM OF THE VELOCITY FIELD IN FFF SPACE
! NI >> 1ST DIMENSION SIZE (IN THE X(OR R) DIRECTION)
! NJ >> 2ND DIMENSION SIZE (IN THE AZIMUTHAL DIRECTION)
! NK >> 3RD DIMENSION SIZE (IN THE AXIAL DIRECTION)
! [DEPENDENCIES]:
! 1. HORFFT(~) @ MOD_SCALAR3
! 2. VERFFT(~) @ MOD_SCALAR3
! 3. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! 4. CHOPDO(~) @ MOD_SCALAR3
! 5. EOMUL(~) @ MOD_LEGOPS
! 6. OEMUL(~) @ MOD_LEGOPS
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 19 2020
!=======================================================================
      TYPE(SCALAR):: RUR,RUP,UZ,PSI,CHI

      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: V,D,T,PFD
      COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: W1,W2

      INTEGER:: I,J,KK,MM,NN,MV,N
      REAL(P8)::KV

      REAL(P8),DIMENSION(NZH):: FFF

      IF(RUR%SPACE.EQ.PPP_SPACE .OR. UZ%SPACE.EQ.PPP_SPACE) THEN
        CALL HORFFT(RUR,-1)                                                ! FORWARD TRANSFORM. PPP -> PFP
        CALL VERFFT(RUR,-1)                                                ! FORWARD TRANSFORM. PFP -> PFF
        CALL HORFFT(RUP,-1)                                                ! FORWARD TRANSFORM. PPP -> PFP
        CALL VERFFT(RUP,-1)                                                ! FORWARD TRANSFORM. PFP -> PFF
        CALL HORFFT(UZ,-1)                                                 ! FORWARD TRANSFORM. PPP -> PFP
        CALL VERFFT(UZ,-1)                                                 ! FORWARD TRANSFORM. PFP -> PFF
      ELSEIF (RUR%SPACE.NE.PFF_SPACE .OR. UZ%SPACE.NE.PFF_SPACE) THEN
        WRITE(*,*) 'PROJECT: RUR,RUP,UZ NOT IN PPP_SPACE OR PFF_SPACE.'
        STOP
      ENDIF

      IF(.NOT.ASSOCIATED(PSI%E).OR. .NOT.ASSOCIATED(CHI%E)) THEN
        WRITE(*,*) 'PROJECT: PSI/CHI NOT ASSOCIATED'
        STOP
      ENDIF

      IF(PSI%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(PSI)
        CALL ALLOCATE(PSI)
      ENDIF

      IF(CHI%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(CHI)
        CALL ALLOCATE(CHI)
      ENDIF

      IF(SIZE(PSI%E,3).NE.NXCHOPDIM.OR.SIZE(CHI%E,3).NE.NXCHOPDIM) THEN
        WRITE(*,*) 'PROJECT:SIZE OF PSI,CHI WRONG.'
        WRITE(*,*) 'NXCHOPDIM=',NXCHOPDIM
        WRITE(*,*) 'SIZEOF (PSI)=[',SIZE(PSI%E,1),SIZE(PSI%E,2),&
                                    SIZE(PSI%E,3),']'
        WRITE(*,*) 'SIZEOF (CHI)=[',SIZE(CHI%E,1),SIZE(CHI%E,2),&
                                    SIZE(CHI%E,3),']'
        STOP
      ENDIF

      RUR%E(:NZ,:NTCHOP,NXCHOP+1:NXCHOPDIM)=&
                                           RUR%E(:NZ,:NTCHOP,NXCHOPU:NX) ! TAKING CARE OF THE CONJUGATE COEFFICIENTS.
      RUP%E(:NZ,:NTCHOP,NXCHOP+1:NXCHOPDIM)=&
                                           RUP%E(:NZ,:NTCHOP,NXCHOPU:NX) ! COPYING THEM OVER SO THAT : K:[0,1,2,3,4,-3,-2,-1]
       UZ%E(:NZ,:NTCHOP,NXCHOP+1:NXCHOPDIM)=&
                                            UZ%E(:NZ,:NTCHOP,NXCHOPU:NX)

      ALLOCATE( V(NZCHOP,NZH) )
      ALLOCATE( D(NZCHOP,NZH) )
      ALLOCATE( T(NZCHOP,NZH) )
      ALLOCATE( PFD(NZH,NZCHOP+1) )
      ALLOCATE( W1(NZCHOP,NXCHOPDIM) )
      ALLOCATE( W2(NZCHOP,NXCHOPDIM) )
     
      W1=0
      FFF=(1-TFM%X(:NZH)**2)/TFM%W(:NZH)

      DO MM=1,NTCHOP
        MV=M(MM)
        NN=NZCHOPS(MM)
        PFD(:NZH,:NN+1)=BANMUL( TFM%PF(:NZH,:NN+1,MM), &
             BAND_LEG_XXDX(NN+1,MV,TFM%NORM(:,MM)),2 )
        DO I=1,NN
          N=MAX(1,MV+I-1)
          V(I,:NZH)= TFM%PF(:NZH,I,MM)/(N*(N+1)*FFF)
          D(I,:NZH)= PFD(:NZH,I)/(N*(N+1)*FFF)
          T(I,:NZH)= TFM%PF(:NZH,I,MM)*TFM%W(:NZH)   
        ENDDO

        ! Debug
!        WRITE(6,*) 'D,1:'
!        WRITE(6,*) SIZE(D(:NN,:),1)
!        WRITE(6,*) 'D,2:'
!        write(6,*) SIZE(D(:NN,:),2)

        IF(MV.NE.0) THEN
          CALL EOMUL(V(:NN,:),RUR%E(:NZ,MM,:NXCHOPDIM),W1(:NN,:))
        ENDIF
        CALL OEMUL(D(:NN,:),RUP%E(:NZ,MM,:NXCHOPDIM),W2(:NN,:))

        PSI%E(:NN,MM,:) = -(IU*MV)*W1(:NN,:) -W2(:NN,:)
        CALL OEMUL(D(:NN,:),RUR%E(:NZ,MM,:NXCHOPDIM),W1(:NN,:))

        IF(MV.NE.0) THEN
          CALL EOMUL(V(:NN,:),RUP%E(:NZ,MM,:NXCHOPDIM),W2(:NN,:))
        ENDIF
        CALL EOMUL(T(:NN,:),UZ%E(:NZ,MM,:NXCHOPDIM),CHI%E(:NN,MM,:))

        DO KK=1,NXCHOPDIM
           KV=AK(MM,KK)
           CHI%E(:NN,MM,KK)=(IU*KV)*W1(:NN,KK)&
                            +MV*KV*W2(:NN,KK)-CHI%E(:NN,MM,KK)
        ENDDO
      ENDDO

      DEALLOCATE( V,D,T,PFD )
      DEALLOCATE( W1,W2 )
      
      CALL CHOPDO(PSI)
      CALL CHOPDO(CHI)

      PSI%E(1,1,:)=PSI%E(1,1,:)-CALCAT1(PSI)/TFM%PF(1,1,1)

      RETURN
      END SUBROUTINE PROJECT
!=======================================================================
      SUBROUTINE EOMUL(A,B,C)
!=======================================================================
! [USAGE]: 
! COMPUTE C(:NI,:NK)=A(:NI,:NJ)*B(:NJ,:NK) WHERE A HAS PATTERN, I.E.,
!    A(2*I  :J)= A(2*I  :NJ-J+1) AND 
!    A(2*I-1:J)=-A(2*I-1:NI-J+1).
! ONLY HALF OF A SHOULD BE GIVEN ON INPUT.
! [PARAMETERS]:
! A >> A SPECIFIC-PATTERN (SEE USAGE) MATRIX OF THE DIMENSION OF NI X NJ
! B >> A MATRIX OF THE DIMENSION OF NJ X NK
! C >> ON EXIT, RETURNS MATRIX MULTIPLICATION OF A*B
! [DEPENDENCIES]:
! 1. OPERATOR(.MUL.) @ MOD_EIG
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 19 2020
!=======================================================================
      REAL(P8),DIMENSION(:,:):: A
      COMPLEX(P8),DIMENSION(:,:):: B,C

      COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: BE,BO
      INTEGER:: NI,NJ,NK,NJH

      NI  = SIZE(A,1)
      NJH = SIZE(A,2)
      NJ  = SIZE(B,1)
      NK  = SIZE(B,2)

      IF(NJH*2 .NE. NJ) THEN
        WRITE(*,*) 'EOMUL: SIZE MISMATCH.'
        WRITE(*,*) 'NJ=',NJ,' NJH=',NJH
        STOP
      ENDIF

      ALLOCATE( BE(NJH,NK) )
      ALLOCATE( BO(NJH,NK) )

      BE = ( B(1:NJH,:) + B(NJ:NJH+1:-1,:) )
      BO = ( B(1:NJH,:) - B(NJ:NJH+1:-1,:) )

      C(1::2,:)= A(1::2,:) .MUL. BE

      IF (NI .GT. 1) THEN
        C(2::2,:)= A(2::2,:) .MUL. BO
      ENDIF

      DEALLOCATE( BE,BO )

      RETURN
      END SUBROUTINE EOMUL
!=======================================================================
      SUBROUTINE OEMUL(A,B,C)
!=======================================================================
! [USAGE]: 
! COMPUTE C(:NI,:NK)=A(:NI,:NJ)*B(:NJ,:NK) WHERE A HAS PATTERN, I.E.,
!    A(2*I  :J)=-A(2*I  :NJ-J+1) AND 
!    A(2*I-1:J)= A(2*I-1:NI-J+1).
! ONLY HALF OF A SHOULD BE GIVEN ON INPUT.
! [PARAMETERS]:
! A >> A SPECIFIC-PATTERN (SEE USAGE) MATRIX OF THE DIMENSION OF NI X NJ
! B >> A MATRIX OF THE DIMENSION OF NJ X NK
! C >> ON EXIT, RETURNS MATRIX MULTIPLICATION OF A*B
! [DEPENDENCIES]:
! 1. OPERATOR(.MUL.) @ MOD_EIG
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 19 2020
!=======================================================================
      REAL(P8),DIMENSION(:,:):: A
      COMPLEX(P8),DIMENSION(:,:):: B,C

      COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: BE,BO
      INTEGER:: NI,NJ,NK,NJH

      NI  = SIZE(A,1)
      NJH = SIZE(A,2)
      NJ  = SIZE(B,1)
      NK  = SIZE(B,2)

      IF(NJH*2 .NE. NJ) THEN
        WRITE(*,*) 'OEMUL: SIZE MISMATCH.'
        WRITE(*,*) 'NJ=',NJ,' NJH=',NJH
        STOP
      ENDIF

      ALLOCATE( BE(NJH,NK) )
      ALLOCATE( BO(NJH,NK) )

      BE = ( B(1:NJH,:) + B(NJ:NJH+1:-1,:) )
      BO = ( B(1:NJH,:) - B(NJ:NJH+1:-1,:) )

      C(1::2,:)= A(1::2,:) .MUL. BO

      C(2::2,:)= A(2::2,:) .MUL. BE

      DEALLOCATE( BE,BO )

      RETURN
      END SUBROUTINE OEMUL
!=======================================================================
      SUBROUTINE NONLIN(PSI,CHI,PSIN,CHIN)
!=======================================================================
! [USAGE]: 
! GIVEN POLOIDAL-TOROIDAL TERMS PSI AND CHI, COMPUTE THE NON-LINEAR
! POLOIDAL-TOROIDAL TERMS PSIN AND CHIN FROM (U X OMEGA) (VPROD).
! SEE EQ. (101) IN MATSUSHIMA AND MARCUS (1997).
! [PARAMETERS]:
! PSI >> TOROIDAL TERM OF THE VELOCITY FIELD
! CHI >> POLOIDAL TERM OF THE VELOCITY FIELD
! PSIN >> TOROIDAL TERM OF THE (U X OMEGA) FIELD
! CHIN >> POLOIDAL TERM OF THE (U X OMEGA) FIELD
! [DEPENDENCIES]:
! 1. OPERATOR(.MUL.) @ MOD_EIG
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 19 2020
!=======================================================================
!   COMPUTE NONLINEAR TERM
!   OUTPUT RETURNS IN PSIN, CHIN
! ----------------------------------------------------
      TYPE(SCALAR),INTENT(IN):: PSI,CHI
      TYPE(SCALAR),INTENT(INOUT)::PSIN,CHIN

      TYPE(SCALAR):: RUR,RUP,UZ
      TYPE(SCALAR):: ROR,ROP,OZ
      TYPE(SCALAR):: W
    real(P8):: count_rate_vprod,count_rate_proj,count_rate_pc2
    integer:: scount_vprod,ecount_vprod,scount_proj,ecount_proj,scount_pc2,ecount_pc2

      CALL ALLOCATE( RUR )
      CALL ALLOCATE( RUP )
      CALL ALLOCATE(  UZ )
      CALL ALLOCATE( ROR )
      CALL ALLOCATE( ROP )
      CALL ALLOCATE(  OZ )
      
      !> DEALIASING
      CALL CHOPSET(3)

    ! open(27,FILE='time_pc2.dat',STATUS='unknown',ACTION='WRITE',POSITION='APPEND')
    ! call system_clock(scount_pc2,count_rate_pc2)
      CALL PC2VEL(PSI,CHI,RUR,RUP,UZ)                                    ! NOTE: INPUT AND OUTPUTS ARE STILL IN FFF SPACE
      CALL PC2VOR(PSI,CHI,ROR,ROP,OZ)                                    ! NOTE: INPUT AND OUTPUTS ARE STILL IN FFF SPACE
    ! call system_clock(ecount_pc2)
    ! WRITE(27,*) "vprod",real(ecount_pc2-scount_pc2,p8)/real(count_rate_pc2,p8)
    ! close(27)

    ! open(28,FILE='time_vprod.dat',STATUS='unknown',ACTION='WRITE',POSITION='APPEND')
    ! call system_clock(scount_vprod,count_rate_vprod)
      CALL VPROD(RUR,RUP,UZ,ROR,ROP,OZ)                                  ! NOTE: VPROD BRINGS EVERY SCALAR-TYPE VARIABLES TO PPP SPACE AND REMEMBER THAT ROR, ROP AND OZ CONTAIN THE RESULT OF THE PRODUCT OF (UR,UP,UZ) X (OR,OP,OZ) IN PPP SPACE
    ! call system_clock(ecount_vprod)
    ! WRITE(28,*) "vprod",real(ecount_vprod-scount_vprod,p8)/real(count_rate_vprod,p8)
    ! close(28)

      CALL CHOPSET(-3)

      CALL DEALLOCATE( RUR )
      CALL DEALLOCATE( RUP )
      CALL DEALLOCATE( UZ  )

      !> CALCULATING NONLINEAR TERMS PSIN, CHIN
      CALL ALLOCATE(W)
      W%LN=0

    ! open(29,FILE='time_proj.dat',STATUS='unknown',ACTION='WRITE',POSITION='APPEND')
    ! call system_clock(scount_proj,count_rate_proj)
      CALL PROJECT(ROR,ROP,OZ,PSIN,W)                                    ! NOW PSIN = PSI_NONLINEAR AND W = DELSQ(CHI_NONLINEAR)
    ! call system_clock(ecount_proj)
    ! WRITE(29,*) "proj",real(ecount_proj-scount_proj,p8)/real(count_rate_proj,p8)
    ! close(29)
      CALL IDEL2(W,CHIN)                                                 ! NOW CHIN = CHI_NONLINEAR

      CALL DEALLOCATE( ROR )
      CALL DEALLOCATE( ROP )
      CALL DEALLOCATE( OZ  )
      CALL DEALLOCATE( W   )

      RETURN
      END SUBROUTINE NONLIN
!=======================================================================
      SUBROUTINE VEL2VOR(RUR,RUP,UZ,ROR,ROP,OZ)
!=======================================================================
! [USAGE]: 
! CONVERT (R*UR,R*UP,UZ) IN FFF SPACE
! TO      (R*OR,R*OP,OZ) IN FFF SPACE
! [PARAMETERS]:
! RUR >> R*U_R TERM IN A SCALAR-TYPE VARIABLE (FFF SPACE)
! RUP >> R*U_P TERM IN A SCALAR-TYPE VARIABLE (FFF SPACE)
! UZ  >>   U_Z TERM IN A SCALAR-TYPE VARIABLE (FFF SPACE)
! ROR >> ON EXIT, R*OMEGA_R TERM IN A SCALAR-TYPE VARIABLE (FFF SPACE)
! ROP >> ON EXIT, R*OMEGA_P TERM IN A SCALAR-TYPE VARIABLE (FFF SPACE)
! OZ  >> ON EXIT,   OMEGA_Z TERM IN A SCALAR-TYPE VARIABLE (FFF SPACE)
! [NOTE]:
! ROR = D/DPHI(UZ) - D/DZ(RUP)
! ROP = D/DZ(RUR) - R*D/DR(UZ)
!  OZ = R^-2*(R*D/DR(RUP) - D/DPHI(RUR))
! [DEPENDENCIES]:
! 1. XXDX(~) @ MOD_LEGOPS
! 2. MULXMDIVXP(~) @ MOD_LEGOPS
! 3. (DE)ALLOCATE(SCALAR) @ MOD_SCALAR3
! [UPDATES]:
! CODED BY JINGE WANG @ APR 16 2021
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: RUR,RUP,UZ,ROR,ROP,OZ,OZ2
      INTEGER:: MM,KK,NN,MV
      REAL(P8):: KV
      TYPE(SCALAR):: W
      
      IF(RUP%SPACE.NE.FFF_SPACE .OR. UZ%SPACE.NE.FFF_SPACE) THEN
        WRITE(*,*) 'VEL2VOR: RUR,RUP,UZ NOT IN FFF_SPACE.'
        STOP
      ENDIF

      IF(ROR%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(ROR)
        CALL ALLOCATE(ROR)
      ENDIF

      IF(ROP%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(ROP)
        CALL ALLOCATE(ROP)
      ENDIF

      IF(OZ%SPACE.NE.FFF_SPACE) THEN
        CALL DEALLOCATE(OZ)
        CALL ALLOCATE(OZ)
      ENDIF
      CALL ALLOCATE(OZ2)

      ! R*D/DR:
      CALL XXDX(RUP,OZ2)                                                !  OZ2 = R*D/DR(RUP)
      CALL XXDX(UZ,ROP)                                                 ! ROP = R*D/DR(UZ)
      
      ROR%LN = 0.D0
      ! FFF SPACE:
      DO MM=1,NTCHOP
        NN=NZCHOPS(MM)
        MV=M(MM)
        DO KK=1,NXCHOPDIM
          KV=AK(MM,KK)
          ROR%E(:NN,MM,KK)=IU*MV*UZ%E(:NN,MM,KK)&                       ! ROR = D/DPHI(UZ) - D/DZ(RUP)
                            -(IU*KV)*RUP%E(:NN,MM,KK)
          ROP%E(:NN,MM,KK)=    -ROP%E(:NN,MM,KK)&                       ! RUP =-R*D/DR(UZ) + D/DZ(RUR)
                            +(IU*KV)*RUR%E(:NN,MM,KK) 
          OZ2%E(:NN,MM,KK)=     OZ2%E(:NN,MM,KK)&                       ! OZ2*R^2 = R*D/DR(RUP) - D/DPHI(RUR)
                            -(IU*MV)*RUR%E(:NN,MM,KK)                  
        ENDDO
      ENDDO

      ! CALL MULXMDIVXP(OZ2, OZ,.TRUE.)
      CALL RTRAN(OZ2,1)
      DO NN = 1,SIZE(OZ2%E,1)
        OZ2%E(NN,:,:) = OZ2%E(NN,:,:)/TFM%R(NN)**2.D0                   ! OZ2*R^2 -> OZ2
      ENDDO
      CALL RTRAN(OZ2,-1)
      OZ = OZ2

      CALL DEALLOCATE(OZ2)


      RETURN
      END SUBROUTINE VEL2VOR      
!=======================================================================
subroutine VPROD_PFF(RUR,RUP,UZ,ROR,ROP,OZ)
!=======================================================================
! [USAGE]: 
! COMPUTE A VECTOR PRODUCT IN PFF space
! (PR,PP,PZ) = (UR,UP,UZ) X (OR,OP,OZ)
! [PARAMETERS]:
! RUR >> ON ENTRY, R*U_R TERM FROM PC2VEL (FFF)
! RUP >> ON ENTRY, R*U_P TERM FROM PC2VEL (FFF)
! UZ >> ON ENTRY, U_Z TERM FROM PC2VEL (FFF)
! ROR >> ON ENTRY, R*OMEGA_R TERM FROM VEL2VOR
!        ON EXIT, R*[THE FIRST COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)]
!                 IN PHYSICAL-FUNCTION-FUNCTION SPACE
! ROP >> ON ENTRY, R*OMEGA_P TERM FROM VEL2VOR
!        ON EXIT, R*[THE SECOND COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)]
!                 IN PHYSICAL-FUNCTION-FUNCTION SPACE
! OZ >> ON ENTRY, OMEGA_Z TERM FROM VEL2VOR
!        ON EXIT, THE THIRD COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)
!                 IN PHYSICAL-FUNCTION-FUNCTION SPACE
! [DEPENDENCIES]:
! 1. RTRAN(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY JIGNE WANG @ APR 18 2021
!=======================================================================
      TYPE(SCALAR):: RUR,RUP,UZ
      TYPE(SCALAR):: ROR,ROP,OZ
      TYPE(SCALAR):: W1,W2,W3
      INTEGER:: NI,NJ,NK
      
      IF(RUP%SPACE.NE.FFF_SPACE .OR. UZ%SPACE.NE.FFF_SPACE) THEN
          WRITE(*,*) 'VPROD_PFF: RUR,RUP,UZ NOT IN FFF_SPACE.'
          ! WRITE(*,*) FFF_SPACE
          ! WRITE(*,*) RUP%SPACE,UZ%SPACE
          STOP
      ENDIF

      IF(ROP%SPACE.NE.FFF_SPACE .OR. OZ%SPACE.NE.FFF_SPACE) THEN
          WRITE(*,*) 'VPROD_PFF: ROR,ROP,OZ NOT IN FFF_SPACE.'
          STOP
      ENDIF

      CALL ALLOCATE(W1)
      CALL ALLOCATE(W2)
      CALL ALLOCATE(W3)
      W1 = RUR
      W2 = RUP
      W3 = UZ
      CALL RTRAN(W1,1)                                                     ! FFF -> PFF SPACE
      CALL RTRAN(W2,1)                                                     ! FFF -> PFF SPACE
      CALL RTRAN(W3,1)                                                      ! FFF -> PFF SPACE
      CALL RTRAN(ROR,1)                                                     ! FFF -> PFF SPACE
      CALL RTRAN(ROP,1)                                                     ! FFF -> PFF SPACE
      CALL RTRAN(OZ,1)                                                      ! FFF -> PFF SPACE
    
      NI=SIZE(W1%E,1)
      NJ=SIZE(W1%E,2)
      NK=SIZE(W1%E,3)

      IF(ABS(W2%E(2,2,2)).GT.1.D-15 .OR.  ABS(W3%E(2,2,2)).GT.1.D-15) THEN
          WRITE(*,*) 'VPROD_PFF: First 3 entries NOT axisymmetric and axial-invariant.'
          !STOP
      ENDIF
      
      CALL VPRODSUB_PFF(W1%E(:NI,1,1),W2%E(:NI,1,1),W3%E(:NI,1,1),ROR%E,ROP%E,OZ%E,NI,NJ,NK)

      CALL DEALLOCATE(W1)
      CALL DEALLOCATE(W2)
      CALL DEALLOCATE(W3)

      RETURN
      END subroutine VPROD_PFF
!=======================================================================
SUBROUTINE VPRODSUB_PFF(RUR,RUP,UZ,ROR,ROP,OZ,NI,NJ,NK)
!=======================================================================
! [USAGE]: 
! COMPUTE A VECTOR PRODUCT (PR,PP,PZ) = (UR,UP,UZ) X (OR,OP,OZ)
! THE RESULTS ARE THEN STORED INTO ROR, ROP AND OZ.
! [PARAMETERS]:
! RUR >> ON ENTRY, R*U_R IN PFF SPACE
! RUP >> ON ENTRY, R*U_P IN PFF SPACE
! UZ >> ON ENTRY, U_Z IN PFF SPACE
! ROR >> ON ENTRY, R*OMEGA_R IN PFF SPACE
!        ON EXIT, R*[THE FIRST COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)]
!                 IN PHYSICAL-FUNCTION-FUNCTION SPACE
! ROP >> ON ENTRY, R*OMEGA_P IN PFF SPACE
!        ON EXIT, R*[THE SECOND COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)]
!                 IN PHYSICAL-FUNCTION-FUNCTION SPACE
! OZ >> ON ENTRY, OMEGA_Z IN PFF SPACE
!        ON EXIT, THE THIRD COMPONENT OF (UR,UP,UZ) X (OR,OP,OZ)
!                 IN PHYSICAL-FUNCTION-FUNCTION SPACE
! NI >> 1ST DIMENSION SIZE (IN THE X(OR R) DIRECTION)
! NJ >> 2ND DIMENSION SIZE (IN THE AZIMUTHAL DIRECTION)
! NK >> 3RD DIMENSION SIZE (IN THE AXIAL DIRECTION)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER:: NI,NJ,NK
      COMPLEX(P8),DIMENSION(NI):: RUR,RUP,UZ
      COMPLEX(P8),DIMENSION(NI,NJ,NK):: ROR,ROP,OZ

      REAL(P8):: A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3
      INTEGER:: KK,MM,NN
      
      !> CROSS PRODUCT
      DO KK=1,NK
          DO MM=1,NJ
              DO NN=1,NI
                  A1=REAL(RUR(NN))
                  A2=REAL(RUP(NN))
                  A3=REAL(UZ(NN))
                  C1=AIMAG(RUR(NN))
                  C2=AIMAG(RUP(NN))
                  C3=AIMAG(UZ(NN))
                  B1=REAL(ROR(NN,MM,KK))
                  B2=REAL(ROP(NN,MM,KK))
                  B3=REAL(OZ(NN,MM,KK))
                  D1=AIMAG(ROR(NN,MM,KK))
                  D2=AIMAG(ROP(NN,MM,KK))
                  D3=AIMAG(OZ(NN,MM,KK))
                  ROR(NN,MM,KK)=CMPLX(A2*B3-A3*B2-C2*D3+C3*D2, A2*D3+C2*B3-A3*D2-C3*B2,P8)
                  ROP(NN,MM,KK)=CMPLX(A3*B1-A1*B3-C3*D1+C1*D3, A3*D1+C3*B1-A1*D3-C1*B3,P8)
                   OZ(NN,MM,KK)=CMPLX(A1*B2-A2*B1-C1*D2+C2*D1, A1*D2+C1*B2-A2*D1-C2*B1,P8)&
                              /TFM%R(NN)**2.
              ENDDO
          ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE VPRODSUB_PFF
!=======================================================================
END MODULE MOD_LEGOPS