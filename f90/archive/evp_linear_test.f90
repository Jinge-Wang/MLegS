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
PROGRAM EVP_LINEAR
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
          TYPE(SCALAR):: PSIU,CHIU,RURU,RUPU,UZU,RORU,ROPU,OZU
          TYPE(SCALAR):: PSI1,CHI1,PSI2,CHI2,PSI3,CHI3
          INTEGER     :: MIND,AKIND,MREAD
          REAL(P8)    :: AKREAD
          COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: H ! Assembled NS operator
          COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: EIG_VEC_3D ! Assembled EigVec matrix of (rur,rup,uz) in PPP space 
          LOGICAL:: SP
    
    
          INTEGER     :: I,J,K
    
          WRITE(*,*) 'PROGRAM STARTED'
          CALL PRINT_REAL_TIME()   ! @ MOD_MISC

          call execute_command_line('echo "%read.input" | ./bin/init_exec')
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
    
          ! LOAD INITIAL PSI (PSI0) AND CHI (CHI0) DATA
          ! PERFORM 'INIT' AS A PREREQUISITE TO MAKE PSI0 AND CHI0
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI0,PSI0)
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI0,CHI0)
    
          CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)
          CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0) 
          CALL RTRAN(ROR0,1) 
          CALL RTRAN(ROP0,1) 
          CALL RTRAN( OZ0,1) 
          CALL RTRAN(RUR0,1) 
          CALL RTRAN(RUP0,1) 
          CALL RTRAN( UZ0,1) 
          CALL MSAVE(ROR0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'ROR0.dat')
          CALL MSAVE(ROP0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'ROP0.dat')
          CALL MSAVE( OZ0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'OZ0.dat')
          CALL MSAVE(RUR0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'RUR0.dat')
          CALL MSAVE(RUP0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'RUP0.dat')
          CALL MSAVE( UZ0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'UZ0.dat')

          CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)
          CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0)

          ! Test: MULXMDIVXP for test function: test failed. MULXMDIVXP does not give correct result

          ! SP = .TRUE.
          ! CALL ALLOCATE(OZU)
          ! CALL MULXMDIVXP(OZ0,OZU, SP)
          ! CALL RTRAN(OZU,1)
          ! CALL MSAVE( OZU%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'OZU_MULX.dat')

          ! CALL RTRAN(OZ0,1)
          ! DO I = 1,SIZE(OZ0%E,1)
          !   OZ0%E(I,:,:) = OZ0%E(I,:,:)/TFM%R(I)**2.D0
          ! ENDDO
          ! CALL MSAVE( OZ0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'OZ0_DRSQ.dat' )
          ! CALL MSAVE( OZU%E - OZ0%E, TRIM(ADJUSTL(FILES%SAVEDIR))//'OZU_DIFF.dat')
          ! CALL DEALLOCATE(OZU)

          SP = .TRUE.
          CALL ALLOCATE(OZU)
          CALL MULXMDIVXP(UZ0,OZU, SP)
          CALL RTRAN(OZU,1)
          CALL MSAVE( OZU%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'OZU_MULX.dat')

          CALL RTRAN(UZ0,1)
          DO I = 1,SIZE(UZ0%E,1)
            UZ0%E(I,:,:) = UZ0%E(I,:,:)/TFM%R(I)**2.D0
          ENDDO
          CALL MSAVE( UZ0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'OZ0_DRSQ.dat' )
          CALL MSAVE( OZU%E - UZ0%E, TRIM(ADJUSTL(FILES%SAVEDIR))//'OZU_DIFF.dat')
          CALL DEALLOCATE(OZU)



        ! ! Test VPROD_PFF: test passed
        !   RUR0%E = cmplx(1.D0,2.D0)
        !   RUP0%E = 1.D0
        !   UZ0%E = 1.D0 
        !   CALL VPROD_PFF(ROR0,ROP0,OZ0,RUR0,RUP0,UZ0)
        ! !   CALL TOFF(RUR0)
        ! !   CALL TOFF(RUP0)
        ! !   CALL TOFF(UZ0)
        ! !   CALL TOFF(RUR0)
        ! !   CALL TOFF(RUP0)
        ! !   CALL TOFF(UZ0)
        ! !   CALL RTRAN(RUR0,1)
        ! !   CALL RTRAN(RUP0,1)
        ! !   CALL RTRAN(UZ0,1)
        !   CALL MSAVE(RUP0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'VPROD_PFF_RESULT.dat')  

        !   CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)
        !   CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0)
        ! !   RUR0%E = cmplx(1.D0,2.D0)
        ! !   RUP0%E = 1.D0
        ! !   UZ0%E = 1.D0

        !   CALL VPROD(ROR0,ROP0,OZ0,RUR0,RUP0,UZ0)
        !   CALL TOFF(RUR0)
        !   CALL TOFF(RUP0)
        !   CALL TOFF(UZ0)
        !   CALL RTRAN(RUR0,1)
        !   CALL RTRAN(RUP0,1)
        !   CALL RTRAN(UZ0,1)
        !   CALL MSAVE(RUP0%E,TRIM(ADJUSTL(FILES%SAVEDIR))//'VPRODRESULT.dat')  

    
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
    END PROGRAM EVP_LINEAR
    !=======================================================================