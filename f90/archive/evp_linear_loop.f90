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
PROGRAM EVP_LINEAR_LOOP
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
      USE MOD_INIT_LOOP

      IMPLICIT NONE
      TYPE(SCALAR):: PSI0,CHI0,RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
      TYPE(SCALAR):: PSIU,CHIU,RURU,RUPU,UZU,RORU,ROPU,OZU
      TYPE(SCALAR):: PSI1,CHI1,PSI2,CHI2,PSI3,CHI3
      INTEGER     :: MIND,AKIND,MREAD
      REAL(P8)    :: AKREAD
      COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: H ! Assembled NS operator
    !   COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: EIG_VEC_3D ! Assembled EigVec matrix of (rur,rup,uz) in PPP space 


      INTEGER     :: I,J,K
      INTEGER     :: N_ELL, I_ELL
      REAL(P8),DIMENSION(:),ALLOCATABLE    :: ELL_RANGE
      COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: EIG_VAL_MAT


      WRITE(*,*) 'PROGRAM STARTED'
      CALL PRINT_REAL_TIME()   ! @ MOD_MISC

      ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
      CALL READIN(5)           ! @ MOD_INIT

      
      ! START LOOP:

      ! STEP 1: ELL_RANGE
      WRITE(*,*) 'Number of increments for ELL_RANGE:'
      READ(5,'(I72)') N_ELL
      N_ELL = N_ELL + 1
      ALLOCATE(ELL_RANGE(N_ELL))
      CALL  LINSPACE(1.D0,15.D0,ELL_RANGE)
      WRITE(*,*) 'STEP 1 finished'

      ! STEP 2: LOOP
      DO I_ELL = 1,N_ELL
        ELL = ELL_RANGE(I_ELL)
        ! DEALLOCATE(NRCHOPS)
        ! DEALLOCATE(NTCHOPS)
        ! DEALLOCATE(AK)
        ! DEALLOCATE(M)
        ! DEALLOCATE(TFM%R,TFM%W,TFM%X, &
        !            TFM%LN,TFM%NORM,TFM%PF,TFM%AT0,TFM%AT1,TFM%EX, &
        !            TFM%Z,TFM%THR,TFM%THI,TFM%TH)
        ! CALL LEGINIT(NR,NTH,NX,NRCHOP,NTCHOP,NXCHOP,&
        ! ZLEN,ELL,MKLINK,MINC)
        ! WRITE(*,*) 'ELL updated'
        CALL INIT_LOOP()

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

      ! Diagnostics
      !CALL DIAGNOST(PSI0,CHI0)
      !call PRINT_ENERGY_SPECTRUM(PSI0,CHI0,1)

      ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
      ! RETURNS (R*OMEGA0_R,R*OMEGA0_THETA,OMEGA0_Z) IN FFF SPACE
      CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)

      ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
      ! RETURNS (R*U0_R,R*U0_THETA,U0_Z) IN FFF SPACE
      CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0)

      IF (I_ELL == 1) THEN
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
      WRITE(*,*) 'AK = '
      WRITE(*,'(F8.3)') (AKIND-1)*2*PI/ZLEN
  101  FORMAT(' *** POSSIBLE INDICES: 1 ~ ', I3)
  102  FORMAT(' M  = ', I3, ' *(M_IND - 1)')
  103  FORMAT(' AK = 2*PI*AK_IND/', F0.3)
  104  FORMAT(' AK_IND = INPUT - 1          (INPUT <=', I3, ')')
  105  FORMAT('          -(', I3,' - INPUT + 1) (INPUT >', I3, ')')
      ALLOCATE(EIG_VAL_MAT(2*NRCHOPS(MIND),N_ELL))
      ENDIF

      MREAD = M(MIND)
      AKREAD = AK(MIND,AKIND)

      WRITE(*,*) 'AZIMUTHAL AND AXIAL WAVENUMBERS:'
      WRITE(*,106) MREAD, AKREAD
  106  FORMAT('M = ',I8,' & K = ',F8.3)

      ALLOCATE(H(2*NRCHOPS(MIND),2*NRCHOPS(MIND)))
    !   ALLOCATE(EIG_VEC_3D(NR,2*NRCHOPS(MIND),3))

      CALL ALLOCATE(PSIU)
      CALL ALLOCATE(CHIU)
      CALL ALLOCATE(RURU)
      CALL ALLOCATE(RUPU)
      CALL ALLOCATE( UZU)
      CALL ALLOCATE(RORU)
      CALL ALLOCATE(ROPU)
      CALL ALLOCATE( OZU)

      CALL ALLOCATE(PSI1)
      CALL ALLOCATE(CHI1)

      CALL ALLOCATE(PSI2)
      CALL ALLOCATE(CHI2)

      CALL ALLOCATE(PSI3)
      CALL ALLOCATE(CHI3)

      ! Use test vector to determine operator H
      DO I = 1,NRCHOPS(MIND)
        PSIU%LN = 0.D0 
        CHIU%LN = 0.D0 ! DEL2CHI
        PSIU%E = 0.D0
        CHIU%E = 0.D0 ! DEL2CHI
        PSIU%E(I,MIND,AKIND) = 1.D0

        CALL CHOPSET(3)

        ! Del2[Psi,Del2Chi]/Re term: PSI
        IF (VISC%SW.EQ.1) THEN
          IF(VISC%NU.EQ.0.0) THEN
            WRITE(*,*) 'VISC: VISC%NU CANNOT BE ZERO'
            STOP
          ENDIF
          CALL DEL2(PSIU,PSI3)
          CALL DEL2(CHIU,CHI3) ! DEL2[DEL2CHI]
          PSI3%E = PSI3%E * VISC%NU ! NU = 1/RE
          CHI3%E = CHI3%E * VISC%NU ! DEL2[DEL2CHI]/RE
        ELSE
          PSI3%E = 0.D0
          CHI3%E = 0.D0          
          PSI3%LN = 0.D0
          CHI3%LN = 0.D0
        ENDIF

        ! W x U' + U X W' term: PSI
        ! NOTE: For this test vecotr, DEL2CHI = 0, so CHI = 0.
        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU) ! CHI (= 0 = DEL2CHI)
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU) ! CHI (= 0 = DEL2CHI)
        
        CALL VPROD(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)
        
        CALL VPROD(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-3)

        PSI1%E = -PSI1%E + PSI2%E + PSI3%E
        CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NRCHOPS(MIND),  I) = PSI1%E(:NRCHOPS(MIND),MIND,AKIND)
        H(NRCHOPS(MIND)+1:,I) = CHI1%E(:NRCHOPS(MIND),MIND,AKIND) ! DEL2CHI
        
        PSIU%LN = 0.D0
        CHIU%LN = 0.D0
        PSIU%E = 0.D0
        CHIU%E = 0.D0
        CHIU%E(I,MIND,AKIND) = 1.D0 ! DEL2CHI

        ! Del2[Psi,Del2Chi]/Re term: CHI
        CALL CHOPSET(3)
        IF (VISC%SW.EQ.1) THEN
          CALL DEL2(PSIU,PSI3)
          CALL DEL2(CHIU,CHI3) ! DEL2[DEL2CHI]
          PSI3%E = PSI3%E * VISC%NU ! NU = 1/RE
          CHI3%E = CHI3%E * VISC%NU ! DEL2[DEL2CHI]/RE
        ELSE
          PSI3%E = 0.D0
          CHI3%E = 0.D0          
          PSI3%LN = 0.D0
          CHI3%LN = 0.D0
        ENDIF
        CALL CHOPSET(-3)  !IDEL2 must have the original chopset.

        CALL IDEL2(CHIU,CHIU) ! DEL2CHI -> CHI so that PC2VOR and PC2VEL calcualte the correct quantities

        CALL CHOPSET(3)

        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU) ! CHI
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU) ! CHI

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
                            = CHI1%E(:NRCHOPS(MIND),MIND,AKIND) ! DEL2CHI

        WRITE(*,107) I, NRCHOPS(MIND),I_ELL
      ENDDO
  107  FORMAT('PROGRESS: ',I3.3,'/',I3.3,'OF #',I3.3)

      ! CALL MCAT(EIGVEC('R',H))
      CALL MCAT(GENEIG(H))
      EIG_VAL_MAT(:,I_ELL) = GENEIG(H);
      CALL MCAT(EIG_VAL_MAT)

      !CALL MSAVE(GENEIG(H), 'eigval_N='//ITOA3(NRCHOP-1)//'M='//ITOA3(MIND)//'K='//ITOA3(AKIND)//'.dat')
      !CALL MSAVE(GENEIG(H), 'eigval_N='//ITOA3(NRCHOP-1)//'M='//ITOA3(MIND)//'K='//ITOA3(AKIND)//'ELL='//ITOA3(INT(ELL))//'.dat')

    !   ! CREATE EIG_VEC_3D Matrix in PPP space:
    !   H = EIGVEC('R',H)
    !   EIG_VEC_3D = 0.D0
    !   DO I = 1,2*NRCHOPS(MIND)
    !     PSIU%LN = 0.D0
    !     CHIU%LN = 0.D0
    !     PSIU%E = 0.D0
    !     CHIU%E = 0.D0     
    !     PSIU%E(:NRCHOPS(MIND),MIND,AKIND) = H(:NRCHOPS(MIND),I)
    !     CHIU%E(:NRCHOPS(MIND),MIND,AKIND) = H(NRCHOPS(MIND)+1:,I)
    !     CALL IDEL2(CHIU,CHIU)

    !     CALL CHOPSET(3)
    !     CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU) ! FFF
    !     CALL CHOPSET(-3)

    !     ! Convert to PFF
    !     CALL RTRAN(RURU,1)
    !     EIG_VEC_3D(:NR,I,1) = RURU%E(:NR,MIND,AKIND)
    !     CALL RTRAN(RUPU,1)
    !     EIG_VEC_3D(:NR,I,2) = RUPU%E(:NR,MIND,AKIND)
    !     CALL RTRAN(UZU ,1)
    !     EIG_VEC_3D(:NR,I,3) = UZU%E (:NR,MIND,AKIND)
    !   ENDDO
    !   CALL MSAVE(EIG_VEC_3D(:NR,:,1),'eigvec_ppp_rur_N=' // ITOA3(NRCHOP-1) // '.dat')
    !   CALL MSAVE(EIG_VEC_3D(:NR,:,2),'eigvec_ppp_rup_N=' // ITOA3(NRCHOP-1) // '.dat')
    !   CALL MSAVE(EIG_VEC_3D(:NR,:,3),'eigvec_ppp__uz_N=' // ITOA3(NRCHOP-1) // '.dat')

      ! DEALLOCATE:
      DEALLOCATE( H )
    !   DEALLOCATE(EIG_VEC_3D)

      CALL DEALLOCATE(PSI1)
      CALL DEALLOCATE(CHI1)
      CALL DEALLOCATE(PSI2)
      CALL DEALLOCATE(CHI2)
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

    END DO

      CALL MSAVE(EIG_VAL_MAT, TRIM(ADJUSTL(FILES%SAVEDIR))//'eigval_matrix.dat')
      WRITE(*,*) ''
      WRITE(*,*) 'PROGRAM FINISHED'
      CALL PRINT_REAL_TIME()  ! @ MOD_MISC

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================
! N/A
!=======================================================================
END PROGRAM EVP_LINEAR_LOOP
!=======================================================================