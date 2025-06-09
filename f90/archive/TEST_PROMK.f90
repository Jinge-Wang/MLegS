PROGRAM TEST_PROMK
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
          REAL(P8),DIMENSION(:,:),ALLOCATABLE:: ENE_MK
    
          WRITE(*,*) 'PROGRAM STARTED'
          CALL PRINT_REAL_TIME()  ! @ MOD_MISC


          CALL READIN(5)
          ! CREATE PPP SCALAR TYPE
          CALL ALLOCATE(A,PPP_SPACE)
          CALL ALLOCATE(B,PPP_SPACE) 


          CALL LAY(B,COSM)

          CALL LAY(A,COSM)

      !     CALL LAYP(A,COSM)
      !     CALL LAYP(B,COSM)
          
          WRITE(*,*) 'A:'
          WRITE(*,*) A%E(2,1,1)

          WRITE(*,*) 'B:'
          WRITE(*,*) B%E(1,1,1)

          CALL HORFFT(A,-1)
          CALL VERFFT(A,-1)

          CALL HORFFT(B,-1)
          CALL VERFFT(B,-1)

          ALLOCATE(ENE_MK(NTCHOP,NXCHOP))
          ENE_MK = PRODUCT_MK(A,B) ! integrate [F*G/(1-X)^2] * (1-x)^2 over the whole domain
          PRINT *, 'ENERGY SPECTRUM  E(M,K=0,K=1,K=2...):'
          WRITE(*,FMT=111) ((ENE_MK(I,K),K=1,7),I=1,7)
          WRITE(*,*) SUM(ENE_MK)
          WRITE(*,*) ENE_MK(:NTH,1)

          WRITE(*,*) B%E(1,:,1)

      !     ENE_MK = PRODUCT_MK_MODIFIED(A,B) ! integrate [F*G/(1-X)^2] * (1-x)^2 over the whole domain
      !     PRINT *, 'ENERGY SPECTRUM MODIFIED  E(M,K=0,K=1,K=2...):'
      !     WRITE(*,FMT=111) ((ENE_MK(I,K),K=1,7),I=1,7)
      !     PRINT *,'---'
      !     PRINT *, 'ENERGY SPECTRUM  E(M,K=-1,K=-2,K=-3...):'
      !     WRITE(*,FMT=111) ((ENE_MK(I,NX+1-K),K=1,7),I=1,7)

      !     WRITE(*,*) 'A:'
      !     WRITE(*,*) A%E(1,1,1)
      !     WRITE(*,*) A%E(1,2,1)
      !     WRITE(*,*) A%E(1,1,2)
      !     WRITE(*,*) A%E(1,2,2)

      !     WRITE(*,*) 'B:'
      !     WRITE(*,*) B%E(1,1,1)
      !     WRITE(*,*) B%E(1,3,1)
      !     WRITE(*,*) B%E(1,1,2)
      !     WRITE(*,*) B%E(1,2,2)

    111 FORMAT(7E14.6)
    
          WRITE(*,*) ''
          WRITE(*,*) 'PROGRAM FINISHED'
          CALL PRINT_REAL_TIME()  ! @ MOD_MISC
    
    END PROGRAM TEST_PROMK
    !=======================================================================