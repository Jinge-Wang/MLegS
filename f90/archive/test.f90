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
!     LAST UPDATE ON OCTOBER 27, 2020
!
!=======================================================================
PROGRAM TEST
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
      TYPE(SCALAR):: RUR,RUP,UZ,RUR2,RUP2,UZ2,UP
      INTEGER     :: I, J, K   ! EXTRA INTEGER VARIABLES (FOR ITER.)
      REAL(P8)    :: X, Y, Z   ! EXTRA REAL(P8) VARIABLES
      ! REAL(P8),DIMENSION(:),ALLOCATABLE    :: ELL_RANGE
      REAL(P8),DIMENSION(:),ALLOCATABLE:: EM,EK
      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: ENE
      INTEGER(P4) :: ANUMBER

      WRITE(*,*) 'PROGRAM STARTED'
      WRITE(*,*) HUGE(ANUMBER)

      ! allocate(ELL_RANGE(10))
      ! CALL LINSPACE(0.D0,1.D0,10,ELL_RANGE)

      CALL PRINT_REAL_TIME()  ! @ MOD_MISC

      CALL READIN(5)

      CALL ALLOCATE(A)
      CALL ALLOCATE(B)

      ! !Test: Msave
      ! CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI0, B)
      ! CALL TOFP(B)
      ! WRITE(*,*) SIZE(B%E,1),SIZE(B%E,2),SIZE(B%E,3)
      ! WRITE(*,*) B%E(NR,NTH,NX)
      ! call MSAVE(B%E(:NR,:NTH,:NX), "TEST.dat")

      ! Test Energy: Case 1 --- Gauss field w/ q = 1
      CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI0, A)
      CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI0, B)
      

      ALLOCATE(EM(NTCHOP))
      ALLOCATE(EK(NXCHOP))

      CALL ENEMON(A,B,EM,EK)
      WRITE(*,*) SUM(EM)
      WRITE(*,*) SUM(EK)
      ! WRITE(*,*) TFM%W(1)
      ! WRITE(*,*) ZLEN
      

      CALL PRINT_ENERGY_SPECTRUM(A,B,0)

      CALL CHOPSET(1) ! in LEGINIT, I manually added 1 to NRchopdim
      CALL PRINT_ENERGY_SPECTRUM(A,B,1)

      ! WRITE(*,*) (1-EXP(-TFM%R**2))/TFM%R

      CALL ALLOCATE(RUR)
      CALL ALLOCATE(RUP)
      CALL ALLOCATE(UZ)
      CALL ALLOCATE(UP)

      ! goto 101
      ! CALL CHOPSET(3)
      CALL PC2VEL(A,B,RUR,RUP,UZ)
      CALL MULXMDIVXP(RUP,UP,.TRUE.)

      CALL PC2VEL(A,B,RUR,RUP,UZ)
      CALL DIVXM(RUP,UZ)
      CALL DIVXM(UZ,RUP)
      CALL CHOPSET(-1)

      ! CALL RTRAN(RUR,1)
      ! CALL RTRAN(RUP,1)
      ! CALL RTRAN(UZ,1)
      ! WRITE(*,*) RUR%LN 
      ! WRITE(*,*) RUP%LN 
      ! WRITE(*,*) UZ%LN 

      CALL TOFP(RUP)
      CALL TOFP(RUR)
      CALL TOFP(UZ)
      CALL TOFP(UP)

      ! CALL ALLOCATE(RUR2)
      ! CALL ALLOCATE(RUP2)
      ! CALL ALLOCATE(UZ2)
      ! CALL PC2VEL(A,B,RUR2,RUP2,UZ2)
      ! CALL TOFP(RUP2)
      ! CALL TOFP(RUR2)
      ! CALL TOFP(UZ2)

      ! WRITE(*,*) SIZE(TFM%R,1)
      ! WRITE(*,*) SIZE(RUP%E,1),SIZE(RUP%E,2),SIZE(RUP%E,3)

      WRITE(*,*) "U_theta/R:"
      ! WRITE(*,*) RUP%E(:NR,2,2)/TFM%R
      WRITE(*,*) REAL(UP%E(:NR,1,1))

      WRITE(*,*) "rU_theta/(1-x)^2:"
      WRITE(*,*) REAL(RUP%E(:NR,1,1))
      ! WRITE(*,*) DOT_PRODUCT(REAL((RUP%E(:NR,1,1)/TFM%R))**2,TFM%W)*2*PI*ZLEN*ELL2/2
      ! WRITE(*,*) ZLEN
      ! WRITE(*,*) ELL2

      ! WRITE(*,*) SUM(RUR%E(:NR,1,1)/TFM%R)
      ! WRITE(*,*) SUM(UZ%E(:NR,1,1)/TFM%R)


101   DEALLOCATE(EM,EK)
      
      ! A%E = CMPLX(1.D0, 1.D0)

      ! CALL TOFF(A)
      
      ! CALL MULXMDIVXP(A,A)
      ! CALL TOFP(A)

      ! CALL MCAT(A%E(:,:,1))


      CALL DEALLOCATE( A )
      CALL DEALLOCATE( B )
      CALL DEALLOCATE(RUR)
      CALL DEALLOCATE(RUP)
      CALL DEALLOCATE(UZ)
      CALL DEALLOCATE(UP)

      WRITE(*,*) ''
      WRITE(*,*) 'PROGRAM FINISHED'
      CALL PRINT_REAL_TIME()  ! @ MOD_MISC

END PROGRAM TEST
!=======================================================================