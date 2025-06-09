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
PROGRAM EVP_LINEAR_NS
!=======================================================================
! [USAGE]: 
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! [UPDATES]:
! LAST UPDATE ON DEC 21, 2020
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
      TYPE(SCALAR):: RDURU,RDUPU,DUZU
      INTEGER     :: MIND,AKIND,MREAD
      REAL(P8)    :: AKREAD,RE
      COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: H
      COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: HH
      REAL(P8),DIMENSION(:,:),ALLOCATABLE:: EF, EFG

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

      ! LOAD INITIAL PSI (PSI0) AND CHI (CHI0) DATA
      ! PERFORM 'INIT' AS A PREREQUISITE TO USE Q-VORTEX AS PSI0/CHI0
      CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI0,PSI0)
      CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI0,CHI0)

      ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
      ! RETURNS (R*OMEGA0_R,R*OMEGA0_THETA,OMEGA0_Z) IN FFF SPACE
      CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)

      ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
      ! RETURNS (R*U0_R,R*U0_THETA,U0_Z) IN FFF SPACE
      CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0)

      WRITE(*,*) ''
      WRITE(*,102) MINC
      WRITE(*,*) 'M_IND = INPUT - 1'
      WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M) INDEX; M_IND ='
      WRITE(*,101) NTCHOPDIM
      READ(5,'(I72)') MIND
      WRITE(*,*) ''
      WRITE(*,103) ZLEN
      WRITE(*,104) NXCHOP
      WRITE(*,105) NXCHOPDIM, NXCHOP
      WRITE(*,*) 'CHOOSE AN AXIAL WAVENUMBER (AK) INDEX; AK_IND ='
      WRITE(*,101) NXCHOPDIM
      READ(5,'(I72)') AKIND
      WRITE(*,*) ''
      WRITE(*,*) 'DEFINE THE REYNOLDS #; RE = '
      WRITE(*,*) 'IF INPUT <= 0, SET RE = +INF (INVISCID)'
      READ(5,'(F72.0)') RE
 101  FORMAT(' *** POSSIBLE INDICES: 1 ~ ', I3)
 102  FORMAT(' M  = ', I3, ' *(M_IND - 1)')
 103  FORMAT(' AK = 2*PI*AK_IND/', F0.3)
 104  FORMAT(' AK_IND = INPUT - 1          (INPUT <=', I3, ')')
 105  FORMAT('          -(', I3,' - INPUT + 1) (INPUT >', I3, ')')

      MREAD = M(MIND)
      AKREAD = AK(MIND,AKIND)

      WRITE(*,*) 'AZIMUTHAL AND AXIAL WAVENUMBERS:'
      WRITE(*,106) MREAD, AKREAD
 106  FORMAT('M = ',I8,' & K = ',F8.3)

      ALLOCATE( H(2*NRCHOPS(MIND),2*NRCHOPS(MIND)) )
      ALLOCATE( HH(NR,2*NRCHOPS(MIND),3) )
      ALLOCATE( EF(NRCHOPS(MIND),2*NRCHOPS(MIND)) )
      ALLOCATE( EFG(NR,2*NRCHOPS(MIND)) )

      CALL ALLOCATE(PSIU)
      CALL ALLOCATE(CHIU)
      CALL ALLOCATE(RURU)
      CALL ALLOCATE(RUPU)
      CALL ALLOCATE( UZU)
      CALL ALLOCATE(RORU)
      CALL ALLOCATE(ROPU)
      CALL ALLOCATE( OZU)

      CALL ALLOCATE(RDURU)
      CALL ALLOCATE(RDUPU)
      CALL ALLOCATE(DUZU )

      CALL ALLOCATE(PSI3)
      CALL ALLOCATE(CHI3)

      CALL ALLOCATE(PSI1)
      CALL ALLOCATE(CHI1)

      CALL ALLOCATE(PSI2)
      CALL ALLOCATE(CHI2)

      DO I = 1,NRCHOPS(MIND)
        PSIU%LN = 0.D0
        CHIU%LN = 0.D0
        PSIU%E = 0.D0
        CHIU%E = 0.D0
        PSIU%E(I,MIND,AKIND) = 1.D0

        CALL CHOPSET(3)

        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU)
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU)

        IF (RE .GT. 0.D0) THEN
          CALL DEL2(PSIU,PSI3)
          PSI3%E = PSI3%E * RE**(-1.D0)
        ELSE
          PSI3%E = 0.D0
        ENDIF
        CHI3%E = 0.D0

        CALL VPROD(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)
        
        CALL VPROD(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-3)

        PSI1%E = -PSI1%E + PSI2%E + PSI3%E
        CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NRCHOPS(MIND),  I) = PSI1%E(:NRCHOPS(MIND),MIND,AKIND)
        H(NRCHOPS(MIND)+1:,I) = CHI1%E(:NRCHOPS(MIND),MIND,AKIND)
        
        PSIU%LN = 0.D0
        CHIU%LN = 0.D0
        PSIU%E = 0.D0
        CHIU%E = 0.D0
        CHIU%E(I,MIND,AKIND) = 1.D0

        CALL IDEL2(CHIU,CHIU)

        CALL CHOPSET(3)

        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU)
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU)

        IF (RE .GT. 0.D0) THEN
          CHI3%E = PSI3%E
        ELSE
          CHI3%E = 0.D0
        ENDIF
        PSI3%E = 0.D0

        CALL VPROD(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)
        
        CALL VPROD(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-3)

        PSI1%E = -PSI1%E + PSI2%E + PSI3%E
        CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NRCHOPS(MIND),  I+NRCHOPS(MIND)) &
                            = PSI1%E(:NRCHOPS(MIND),MIND,AKIND)
        H(NRCHOPS(MIND)+1:,I+NRCHOPS(MIND)) &
                            = CHI1%E(:NRCHOPS(MIND),MIND,AKIND)

        WRITE(*,107) I, NRCHOPS(MIND)
      ENDDO
 107  FORMAT('PROGRESS: ',I3.3,'/',I3.3)

      CALL MSAVE(H, TRIM(ADJUSTL(FILES%SAVEDIR))// 'matrix_N=' // ITOA3(NRCHOP-1) // '.dat')

      CALL MSAVE(EIGVEC('R',H), TRIM(ADJUSTL(FILES%SAVEDIR))// 'evec_r_N=' // ITOA3(NRCHOP-1) // '.dat')
      CALL MSAVE(EIGVEC('L',H), TRIM(ADJUSTL(FILES%SAVEDIR))// 'evec_l_N=' // ITOA3(NRCHOP-1) // '.dat')
      
      CALL MSAVE(GENEIG(H), TRIM(ADJUSTL(FILES%SAVEDIR))// 'eval_N=' // ITOA3(NRCHOP-1) // '.dat')

      H = EIGVEC('R',H)
      HH = 0.D0

      DO I = 1,2*NRCHOPS(MIND)
        PSIU%LN = 0.D0
        CHIU%LN = 0.D0
        PSIU%E = 0.D0
        PSIU%E(:NRCHOPS(MIND),MIND,AKIND) = H(:NRCHOPS(MIND),I)
        CHIU%E = 0.D0
        CHIU%E(:NRCHOPS(MIND),MIND,AKIND) = H(NRCHOPS(MIND)+1:,I)
        CALL IDEL2(CHIU,CHIU)

        CALL CHOPSET(3)
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU)
        CALL CHOPSET(-3)

        CALL RTRAN(RURU,1)
        HH(:NR,I,1) = RURU%E(:NR,MIND,AKIND)
        CALL RTRAN(RUPU,1)
        HH(:NR,I,2) = RUPU%E(:NR,MIND,AKIND)
        CALL RTRAN(UZU ,1)
        HH(:NR,I,3) = UZU%E (:NR,MIND,AKIND)
      ENDDO

      CALL MSAVE(HH(:,:,1), TRIM(ADJUSTL(FILES%SAVEDIR))//'efun_phys_rur_N=' // ITOA3(NRCHOP-1) // '.dat')
      CALL MSAVE(HH(:,:,2), TRIM(ADJUSTL(FILES%SAVEDIR))//'efun_phys_rup_N=' // ITOA3(NRCHOP-1) // '.dat')
      CALL MSAVE(HH(:,:,3), TRIM(ADJUSTL(FILES%SAVEDIR))//'efun_phys__uz_N=' // ITOA3(NRCHOP-1) // '.dat')

      DEALLOCATE( EF,EFG )
      DEALLOCATE( H,HH )

      CALL DEALLOCATE(PSI1)
      CALL DEALLOCATE(CHI1)
      CALL DEALLOCATE(PSI2)
      CALL DEALLOCATE(CHI2)

      CALL DEALLOCATE(RDURU)
      CALL DEALLOCATE(RDUPU)
      CALL DEALLOCATE(DUZU )

      CALL DEALLOCATE(PSI3)
      CALL DEALLOCATE(CHI3)

      CALL DEALLOCATE(PSIU)
      CALL DEALLOCATE(CHIU)
      CALL DEALLOCATE(RURU)
      CALL DEALLOCATE(RUPU)
      CALL DEALLOCATE( UZU)
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
END PROGRAM EVP_LINEAR_NS
!=======================================================================