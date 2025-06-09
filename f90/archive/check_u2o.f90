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
PROGRAM check_u2o
    !=======================================================================
    ! [USAGE]: 
    ! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
    ! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
    ! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
    ! [UPDATES]:
    ! LAST UPDATE ON DEC 3, 2020
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
          TYPE(SCALAR):: PSI0,CHI0,RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
          TYPE(SCALAR):: RORU,ROPU,OZU
          INTEGER     :: MIND,AKIND,MREAD
          REAL(P8)    :: AKREAD    
    
          INTEGER     :: I,J,K
    
          WRITE(*,*) 'PROGRAM STARTED'
          CALL PRINT_REAL_TIME()   ! @ MOD_MISC
    
          ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
          CALL READIN(5)           ! @ MOD_INIT
    
          CALL ALLOCATE(PSI0)
          CALL ALLOCATE(CHI0)
          CALL ALLOCATE(RUR0)
          CALL ALLOCATE(RUP0)
          CALL ALLOCATE( UZ0)
          CALL ALLOCATE(ROR0)
          CALL ALLOCATE(ROP0)
          CALL ALLOCATE( OZ0)
          CALL ALLOCATE(RORU)
          CALL ALLOCATE(ROPU)
          CALL ALLOCATE( OZU)
    
          ! LOAD INITIAL PSI (PSI0) AND CHI (CHI0) DATA
          ! PERFORM 'INIT' AS A PREREQUISITE TO MAKE PSI0 AND CHI0
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI0,PSI0)
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI0,CHI0)
    
          CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)
          CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0)
          CALL VEL2VOR(RUR0,RUP0,UZ0,RORU,ROPU,OZU)

          CALL RTRAN(ROR0,1)
          CALL RTRAN(ROP0,1)
          CALL RTRAN( OZ0,1)
          CALL RTRAN(RORU,1)
          CALL RTRAN(ROPU,1)
          CALL RTRAN( OZU,1)
          CALL RTRAN(RUR0,1)
          CALL RTRAN(RUP0,1)
          CALL RTRAN( UZ0,1)

          CALL CHOPDO(ROR0)
          CALL MSAVE(ROR0%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'ROR0.dat')
          CALL CHOPDO(ROP0)
          CALL MSAVE(ROP0%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'ROP0.dat')
          CALL CHOPDO(OZ0)
          CALL MSAVE(OZ0%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'OZ0.dat')
          CALL CHOPDO(RORU)
          CALL MSAVE(RORU%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'RORU.dat')
          CALL CHOPDO(ROPU)
          CALL MSAVE(ROPU%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'ROPU.dat')
          CALL CHOPDO(OZU)
          CALL MSAVE(OZU%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'OZU.dat')
          CALL CHOPDO(RUR0)
          CALL MSAVE(RUR0%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'RUR0.dat')
          CALL CHOPDO(RUP0)
          CALL MSAVE(RUP0%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'RUP0.dat')
          CALL CHOPDO(UZ0)
          CALL MSAVE(UZ0%E(:NR,:NTH,:NX), TRIM(ADJUSTL(FILES%SAVEDIR))//'UZ0.dat')

          CALL DEALLOCATE(RORU)
          CALL DEALLOCATE(ROPU)
          CALL DEALLOCATE( OZU)
    
          CALL DEALLOCATE(PSI0)
          CALL DEALLOCATE(CHI0)
          CALL DEALLOCATE(RUR0)
          CALL DEALLOCATE(RUP0)
          CALL DEALLOCATE( UZ0)
          CALL DEALLOCATE(ROR0)
          CALL DEALLOCATE(ROP0)
          CALL DEALLOCATE( OZ0)
    
          WRITE(*,*) ''
          WRITE(*,*) 'PROGRAM FINISHED'
          CALL PRINT_REAL_TIME()  ! @ MOD_MISC
    
    CONTAINS
    !=======================================================================
    !=================== PROGRAM-DEPENDENT SUBROUTINES =====================
    !=======================================================================
    ! N/A
    !=======================================================================
    END PROGRAM check_u2o
    !=======================================================================