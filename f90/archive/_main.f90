!=======================================================================
!
!     ORIGINAL CODE BY TATSUHITO MATSUSHIMA
!     RE-CODED BY SANGJOON LEE
!     DEPARTMENT OF MECHANICAL ENGINEERING
!     UNIVERSITY OF CALIFORNIA AT BERKELEY
!     E-MAIL: SANGJOONLEE@BERKELEY.EDU
!     
!     UC BERKELEY CFD LAB
!     URL: HTTPS://CFD.ME.BERKELEY.EDU/
!
!     NONCOMMERCIAL USE WITH COPYRIGHTED (C) MARK
!     UNDER DEVELOPMENT FOR RESEARCH PURPOSE 
!     LAST UPDATE ON NOVEMBER 22, 2020
!
!=======================================================================
      PROGRAM MAIN
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
      INTEGER:: I,J,K,N

      WRITE(*,*) 'PROGRAM STARTED'
      CALL PRINT_REAL_TIME()  ! @ MOD_MISC

      N = 8           ! # OF DEGREES 8 - TRUNCATION OCCURS
      ELL = 1.D0      ! MAP PARAMETER SET TO BE 1.0
      M = 3           ! USE P^3_3 THROUGH P^3_10

      ALLOCATE(X(1:N))     ! COLLOC. POINTS RANGING FROM -1 TO 1
      ALLOCATE(R(1:N))     ! MAPPED COLLOC. POINTS RANING FROM 0 TO INFTY 
      ALLOCATE(W(1:N))     ! QUADRATURE WEIGHT
      ALLOCATE(F(1:N))     ! TEST FUNCTION. HERE F(R) = (1-EXP(-R**2))/R.
      ALLOCATE(FF(1:N))    ! DIRECT CALCULATION OF L(F) WHERE L = (R^2-ELL^2)/(R^2+ELL^2)*
      ALLOCATE(F_FUN(1:N)) ! FUNCTION IN FUNCTION SPACE (TRANSFORM FROM F)
      ALLOCATE(NORM(1:N))  ! NORMALIZATION FACTORS TO USE NORMALIZED POLYNOMIALS
      ALLOCATE(P_TAB(1:SIZE(X),1:N)) ! NORMALIZED LEGENDRE POLYNOMIAL TABLE

! ! ! EXAMPLE 0 (GAUSSIAN QUADRATURE DEBUG) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       CALL LEG_ZERO(X,R,ELL,W) ! COLLOCATION POINTS IN X [-1, 1]

!       WRITE(*,*) 'N = ', N
!       WRITE(*,*) 'M = ', M

!       WRITE(*,*) "CALCULATE INT_{-1}^{+1}X^I DX: LEG-GAU QUAD VS DIRECT"
!       WRITE(*,*) ""
!       WRITE(*,*) "  I      LEG-GAU             DIRECT_INT           ERR"

!       DO J = 1, 4*N-1
!         DO I = 1, N
!             F(I) = X(I)**J ! F(X) = X**J. J = 0 ~ N-1
!         ENDDO

!                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   ! GAUSS LEGENDRE QUADRATURE INT !
!         F_FUN = 0.                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DO I = 1, N
!         F_FUN(1) = F_FUN(1) + F(I) * W(I)
!       ENDDO

!       WRITE(*,"(I4,3F20.12)") J, F_FUN(1), (1.D0 - (-1.D0)**(J+1))/(J+1.D0),&
!                              (F_FUN(1)-(1.D0 - (-1.D0)**(J+1))/(J+1.D0))   
!       ENDDO
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! EXAMPLE 1 (PROFESSOR'S ITEM 2) (UNINDENT LINE 43 - 100)
! !!!PHYSICAL -> FUNCTION -> OPERATOR ON FUNCTION SP. -> PHYSICAL TEST!!!

!       OPEN(101,FILE='./leg_x_test.out', ACTION="WRITE")

!       CALL LEG_ZERO(X,R,ELL,W) ! COLLOCATION POINTS IN X [-1, 1]

!       WRITE(*,*) 'N = ', N
!       WRITE(*,*) 'M = ', M

!       DO I = 1, N
!         F(I) = (1-EXP(-R(I)**2.))/R(I)
!         ! EXAMPLE F(R) = (1-EXP(R)^2)/R. 
!         ! F(R) ~ R AT 0 AND F(R) ~ 0 AS R -> INF. GOOD FOR THE CURRENT BASIS.
!       ENDDO

!       WRITE(*,*) '                  X                   R'//&
!                  '                   W                   F'
!       DO I= 1, N
!         WRITE(*,1001) X(I), R(I), W(I), F(I)
!       ENDDO

!       NORM = LEG_NORM(N,M)        ! GET NORMALIZATION FACTORS
!       P_TAB = LEG_TBL(X,N,M,NORM) ! FIND NORMALIZED P^M_N VALS

!                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DO J = 1, N                 ! TRANSFORM OF F (PHY -> FUNC) !
!         F_FUN(J) = 0              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DO I = 1, N
!         F_FUN(J) = F_FUN(J) + F(I) * W(I) * P_TAB(I,J)
!       ENDDO
!       ENDDO

!       F_FUN = MATMUL(LEG_X(N, N, M, NORM), F_FUN) ! L = X* OPERATOR IN FUNCTION SPACE.  
!                                                   ! NOTE: L = X* = (R^2-L^2)/(R^2+L^2)*

!                                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DO J = 1, N                ! INVERSE TRANSFORM OF F (FUNC -> PHY) !
!         F(J) = 0                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DO I = 1, N
!         F(J) = F(J) + F_FUN(I)*P_TAB(J,I)
!       ENDDO                      ! EXPECTED TO BE F(R) = (R^2-1^2)/(R^2+1^2)*(1-EXP(R)^2)/R
!       ENDDO

!       DO I = 1, N                ! ... LET'S VALIDATE
!         FF(I) = (R(I)**2. - ELL**2.)/(R(I)**2. + ELL**2.)*(1-EXP(-R(I)**2.))/R(I)
!       ENDDO

!       WRITE(*,*) ''
!       WRITE(101,*) '                R  P>F&OPER>P__L(F) F'//&
!                  '         DIRECT_L(F)               ERROR'
!       WRITE(*,*) '                  R  P>F&OPER>P__L(F) F'//&
!                  '         DIRECT_L(F)               ERROR'
!       DO I= 1, N
!         WRITE(101,1001) R(I), F(I), FF(I), FF(I)-F(I) 
!         WRITE(*,1001) R(I),F(I), FF(I), FF(I)-F(I)
!       ENDDO
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! EXAMPLE 2 (PROFESSOR'S ITEM 6) (UNINDENT LINE 103 - 169)
!!!    HORIZONTAL LAPLACIAN OPERATOR TEST USING THE S-L EQATION     !!!

      ! OPEN(101,FILE='./leg_rat_del2h_test.out', ACTION="WRITE")

      ! CALL LEG_ZERO(X,R,ELL,W) ! COLLOCATION POINTS IN X [-1, 1]

      ! WRITE(*,*) 'N = ', N
      ! WRITE(*,*) 'M = ', M

      ! NORM = LEG_NORM(N,M)        ! GET NORMALIZATION FACTORS
      ! P_TAB = LEG_TBL(X,N,M,NORM) ! FIND NORMALIZED P^M_N VALS

      ! K = 5   ! FIX N AT |M| + 10 -1 (YOU MAY CHANGE THIS)

      ! DO I = 1, N
      ! F(I) = P_TAB(I,K)
      ! ! EXAMPLE F(R) = P^M_L_N(R).
      ! ! FROM S-L, WE KNOW DEL^2_PERP(P^M_L_N) = -(4*N*(N+1)*ELL^2)/(ELL^2+R^2)^2 * P^M_L_N
      ! ! LET'S CHECK THIS EQUALITY
      ! ENDDO


      ! WRITE(*,*) '                  X                   R'//&
      !            '                   W                   F'
      ! DO I= 1, N
      !   WRITE(*,1001) X(I), R(I), W(I), F(I)
      ! ENDDO



      !                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO J = 1, N                 ! TRANSFORM OF F (PHY -> FUNC) !
      !   F_FUN(J) = 0              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO I = 1, N
      !   F_FUN(J) = F_FUN(J) + F(I) * W(I) * P_TAB(I,J)
      ! ENDDO
      ! ENDDO

      ! WRITE(*,*) ''
      ! WRITE(*,*) 'CHECK 1 APPEARS ONLY AT 5TH COEFF. OTHERWISE 0'
      ! WRITE(*,1002) F_FUN

      ! F_FUN = MATMUL(LEG_RAT_DEL2H(N, N, M, ELL, NORM), F_FUN) ! L = DEL^2_PERP OPERATOR IN FUNCTION SPACE.  
                                                

      !                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO J = 1, N                ! INVERSE TRANSFORM OF F (FUNC -> PHY) !
      !   F(J) = 0                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO I = 1, N
      !   F(J) = F(J) + F_FUN(I)*P_TAB(J,I)
      ! ENDDO                      ! EXPECTED TO BE F(R) = -(4*N*(N+1)*ELL^2)/(ELL^2+R^2)^2 * P^M_L_N
      ! ENDDO

      ! DO I = 1, N                ! ... LET'S VALIDATE
      !   FF(I) = -(4*(M+K-1)*(M+K)*ELL**2.)/(ELL**2.+R(I)**2.)**2. * P_TAB(I,K)
      ! ENDDO

      ! WRITE(*,*) ''
      ! WRITE(101,*) '                R  P>F>DEL2P>P_L(F) F'//&
      !            '         EQ(10) RHS          RATIO_ERROR'
      ! WRITE(*,*) '                  R  P>F>DEL2P>P_L(F) F'//&
      !            '         EQ(10) RHS          RATIO_ERROR'
      ! DO I= 1, N
      !   WRITE(101,1001) R(I), F(I), FF(I), FF(I)/F(I) 
      !   WRITE(*,1001) R(I),F(I), FF(I), FF(I)/F(I)
      ! ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! EXAMPLE 3 (CONSTRUCT TRANS. & INV. TRANS. MATRICES) 
!!PHYSICAL -> FUNCTION -> OPERATOR ON FUNCTION SP. -> PHYSICAL TEST!!!
      ALLOCATE(TRANS(1:SIZE(X),1:N)) ! TRANSFORMATION MATRIX (P->F)
      ALLOCATE(INVTR(1:SIZE(X),1:N)) ! INVERSE TRANSFORMATION MATRIX (F->P)
      ALLOCATE(C(1:SIZE(X),1:N))
      OPEN(101,FILE='./leg_x_test.out', ACTION="WRITE")

      CALL LEG_ZERO(X,R,ELL,W) ! COLLOCATION POINTS IN X [-1, 1]

      NORM = LEG_NORM(N,M)        ! GET NORMALIZATION FACTORS
      P_TAB = LEG_TBL(X,N,M,NORM) ! FIND NORMALIZED P^M_N VALS

      WRITE(*,1005) NORM

      WRITE(*,*) 'N = ', N
      WRITE(*,*) 'M = ', M

      DO I = 1, N
        F(I) = P_TAB(I,4)+ 3*P_TAB(I,2)
        ! F(I) = (1-EXP(-R(I)**2.))/R(I)
        ! EXAMPLE F(R) = (1-EXP(R)^2)/R. 
        ! F(R) ~ R AT 0 AND F(R) ~ 0 AS R -> INF. GOOD FOR THE CURRENT BASIS.
      ENDDO
      WRITE(*,*) F
      ! WRITE(*,*) '                  X                   R'//&
      !            '                   W                   F'
      ! DO I= 1, N
      !   WRITE(*,1001) X(I), R(I), W(I), F(I)
      ! ENDDO

      DO I = 1 , SIZE(X)
        DO J = 1 , N
          TRANS(I,J) = W(I)*P_TAB(I,J)
          INVTR(I,J) = P_TAB(J,I)
        ENDDO
      ENDDO
      TRANS(:,6:) = 0.
      ! INVTR(6:,:) = 0.

      WRITE(*, *) 'P2F TRANS'
      WRITE(*, 1005) TRANS
      WRITE(*, *) ''
      WRITE(*, *) 'F2P TRANS'
      WRITE(*, 1005) INVTR
      WRITE(*, *) ''
      WRITE(*, *) 'P2F2P (F2P X P2F)'
      WRITE(*, 1005) TRANSPOSE(MATMUL(TRANSPOSE(INVTR),TRANSPOSE(TRANS)))
      WRITE(*, *) ''
      WRITE(*, *) 'F2P2F (P2F X F2P)'
      WRITE(*, 1005) TRANSPOSE(MATMUL(TRANSPOSE(TRANS),TRANSPOSE(INVTR)))
1005  FORMAT(8F7.3)


      !                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO J = 1, N-M+1             ! TRANSFORM OF F (PHY -> FUNC) !
      !   F_FUN(J) = 0              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO I = 1, N
      !   F_FUN(J) = F_FUN(J) + F(I) * W(I) * P_TAB(I,J)
      ! ENDDO
      ! ENDDO

      F_FUN = MATMUL(TRANSPOSE(TRANS),F)

      WRITE(*,*) F_FUN
      ! F_FUN = MATMUL(LEG_X(N, N, M, NORM), F_FUN) ! L = X* OPERATOR IN FUNCTION SPACE.  
                                                  ! NOTE: L = X* = (R^2-L^2)/(R^2+L^2)*

      !                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO J = 1, N                ! INVERSE TRANSFORM OF F (FUNC -> PHY) !
      !   F(J) = 0                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO I = 1, N
      !   F(J) = F(J) + F_FUN(I)*P_TAB(J,I)
      ! ENDDO                      ! EXPECTED TO BE F(R) = (R^2-1^2)/(R^2+1^2)*(1-EXP(R)^2)/R
      ! ENDDO

      ! F = MATMUL(TRANSPOSE(INVTR),MATMUL(TRANSPOSE(TRANS),F))
      C = MATMUL(TRANSPOSE(INVTR),TRANSPOSE(TRANS))

      WRITE(*, 1005) C
      WRITE(*, *) F
      WRITE(*, *) MATMUL(C,F)
      WRITE(*, *) MATMUL(C,MATMUL(C,F))
      WRITE(*, *) MATMUL(C,MATMUL(C,MATMUL(C,F)))

      DO I = 1, N                ! ... LET'S VALIDATE
        FF(I) = P_TAB(I,4) + 3*P_TAB(I,2) +.1
        ! FF(I) = (R(I)**2. - ELL**2.)/(R(I)**2. + ELL**2.)*(1-EXP(-R(I)**2.))/R(I)
      ENDDO

      WRITE(*,*) ''
      WRITE(101,*) '                R  P>F&OPER>P__L(F) F'//&
                 '         DIRECT_L(F)               ERROR'
      WRITE(*,*) '                  R  P>F&OPER>P__L(F) F'//&
                 '         DIRECT_L(F)               ERROR'
      DO I= 1, N
        WRITE(101,1001) R(I), F(I), FF(I), FF(I)-F(I) 
        WRITE(*,1001) R(I),F(I), FF(I), FF(I)-F(I)
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1001  FORMAT(4F20.12)
1002  FORMAT( F20.12)


      CLOSE(101)
      DEALLOCATE(X,R,W,F,FF,F_FUN,NORM,P_TAB)

      WRITE(*,*) ''
      WRITE(*,*) 'PROGRAM FINISHED'
      CALL PRINT_REAL_TIME()  ! @ MOD_MISC

END PROGRAM MAIN
!=======================================================================
!=======================================================================