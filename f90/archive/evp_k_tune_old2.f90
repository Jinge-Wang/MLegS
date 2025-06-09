PROGRAM EVP_K_TUNE
!=======================================================================
! [USAGE]: 
! USE VPROD_PFF RATHER THAN VPROD
! TUNE K SUCH THAT SIGMA_1 = SIGMA_2
! 
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
! USE VPROD_PFF
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
    INTEGER     :: N_mat,II,JJ,KK
    INTEGER     :: M_IND,K_IND,M_IND2,K_IND2,EIG_TARGET_IND
    REAL(P8)    :: K_VAL,K_VAL0,K_VAL1,K_shift,tol,tol_diff,K_delta,K_step
    COMPLEX(P8) :: EIG_TARGET
    REAL(P8),DIMENSION(:),ALLOCATABLE:: K_tune
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig0, M_eig1, M_eig1_ND
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat      !, EIG_R_mat, EIG_L_mat
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat0, EIG_R_mat0, EIG_L_mat0
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat1, EIG_R_mat1, EIG_L_mat1

    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC

    CALL execute_command_line('echo "%read.input" | ./bin/init_exec')
    CALL READIN(5)           ! @ MOD_INIT

    ! {M,AK}
    WRITE(*,*) 'STEP 1: SELECT {M,AK}'
    WRITE(*,*) ''
    WRITE(*,102) MINC
    WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M) INDEX; M_IND ='
    WRITE(*,101) NTCHOPDIM-1
    READ(5,'(I72)') M_IND
    M_IND = M_IND + 1

    WRITE(*,*) ''
    WRITE(*,103) ZLEN
    WRITE(*,104) NXCHOP-1
    WRITE(*,105) NXCHOPDIM, NXCHOP-1
    WRITE(*,*) 'CHOOSE AN AXIAL WAVENUMBER (AK) INDEX; AK_IND ='
    WRITE(*,101) NXCHOPDIM-1
    READ(5,'(I72)') K_IND
    WRITE(*,*) 'AK = '
    WRITE(*,'(F8.3)') K_IND*2*PI/ZLEN
    K_IND = K_IND + 1

101  FORMAT(' *** POSSIBLE INDICES: 0 ~ ', I3)
102  FORMAT(' M  = ', I3, ' *M_IND')
103  FORMAT(' AK = 2*PI*K_IND/', F0.3)
104  FORMAT(' K_IND = INPUT          (INPUT <=', I3, ')')
105  FORMAT('         -(', I3,' - INPUT) (INPUT  >', I3, ')')

    ! Pick a target Sigma
    WRITE(*,*) ''
    WRITE(*,*) 'CALCULATING EIGENVALUES FOR {M,AK}. MIGHT TAKE MINUTES'
    CALL EIG_MATRIX(M_IND,K_IND,K_VAL,M_mat,M_eig)
    WRITE(*,*) ''
    WRITE(*,*) 'EIGVALS FOR {M,K}:'
    CALL MCAT(M_eig)
    WRITE(*,*) 'CHOOSE THE TARGET EIGVAL #:'
    READ(5,'(I72)') EIG_TARGET_IND
    EIG_TARGET = M_eig(EIG_TARGET_IND)
    WRITE(*,*) EIG_TARGET
    DEALLOCATE(M_mat,M_eig)

    ! {M',AK'}: M' = M + M_bar
    WRITE(*,*) ''
    WRITE(*,*) 'STEP 2: SELECT {M + M_bar}'
    WRITE(*,*) ''
    WRITE(*,102) MINC
    WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M) INDEX; M_IND2 ='
    WRITE(*,101) NTCHOPDIM-1
    READ(5,'(I72)') M_IND2
    M_IND2 = M_IND2 + 1
    N_mat = NRCHOPS(M_IND2)


    K_VAL0 = 8.47 ! A random starting guess
    tol = 1.0E-6 ! 0.1%
    K_delta = 0.001
    K_step = 100

    DO JJ = 1,100
    ! Find how far away
        K_VAL1 = K_VAL0 + K_delta

        ZLEN = 2*PI/K_VAL0
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND2,2,K_VAL0,M_mat0,M_eig0,EIG_R_mat0,EIG_L_mat0)
        WRITE(*,*) K_VAL0

        ZLEN = 2*PI/K_VAL1
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND2,2,K_VAL1,M_mat1,M_eig1,EIG_R_mat1,EIG_L_mat1)
        WRITE(*,*) K_VAL1

        ! DEALLOCATE(EIG_R_mat1)
        ! ALLOCATE(EIG_R_mat1(size(M_mat0,1),size(EIG_L_mat0,2)))
        EIG_R_mat1 = 0.D0
        EIG_R_mat1 = MATMUL((M_mat1 - M_mat0), EIG_R_mat0)
    
    ! Find possible K_shift
        IF(.NOT.ALLOCATED(K_tune)) ALLOCATE(K_tune(N_mat*2))
        DO II = 1,2*N_mat
            ! Gradient Descent: stuck at local minimum
            K_tune(II) = ABS( AIMAG(M_eig0(II)) + AIMAG(DOT_PRODUCT(EIG_L_mat0(:,II),EIG_R_mat1(:,II))) - AIMAG(EIG_TARGET) )
            K_tune(II) = K_tune(II) - ABS( AIMAG(M_eig0(II)) - AIMAG(EIG_TARGET) )
            K_tune(II) = K_tune(II) / ABS(AIMAG(EIG_TARGET))
            K_tune(II) = K_tune(II) / K_delta
            K_tune(II) = K_tune(II) * K_step

            ! ! Newton's method
            ! K_tune(II) = ABS( AIMAG(M_eig0(II)) + AIMAG(DOT_PRODUCT(EIG_L_mat0(II,:),EIG_R_mat1(II,:))) - AIMAG(EIG_TARGET) )
            ! K_tune(II) = K_tune(II) - ABS( AIMAG(M_eig0(II)) - AIMAG(EIG_TARGET) )
            ! K_tune(II) = K_tune(II) / K_delta
            ! K_tune(II) = ABS( AIMAG(M_eig0(II)) - AIMAG(EIG_TARGET) ) / K_tune(II) ! F/F'

            IF (K_tune(II).GT.1.D0) THEN 
                K_tune(II) = 1.D0
            ELSEIF (K_tune(II).LT.-1.D0) THEN
                K_tune(II) = -1.D0
            ENDIF
            IF ((K_VAL0-K_tune(II)).LT.0.D0) THEN
                K_tune(II) = 999 ! GET RID OF NEGATIVE K_VAL
            ENDIF
            ! WRITE(*,*) ''
            ! WRITE(*,*) -K_tune(II)
        ENDDO
        
    ! Find the new K_shift: currently to be the closest
        K_shift = K_tune(MINLOC( ABS(K_tune) ,1))
        WRITE(*,*) 'K_shift:'
        WRITE(*,*) K_shift

    ! ! Calculate the Non-Degen Predicted
    !     ALLOCATE(M_eig1_ND(N_mat*2))
    !     DO II = 1,2*N_mat
    !         M_eig1_ND(II) = M_eig0(II) + K_shift*DOT_PRODUCT(EIG_L_mat0(II,:),EIG_R_mat1(II,:)) ! DOT_PRODUCT(A,B) = SUM(CONJG(A)*B)
    !     ENDDO

    ! Calculate the actual 
        K_VAL1 = K_VAL0 - K_shift
        ZLEN = 2*PI/K_VAL1
        CALL INIT_LOOP()
        CALL EIG_MATRIX(M_IND2,2,K_VAL1,M_mat1,M_eig1,EIG_R_mat1,EIG_L_mat1)

    ! Check if tolerance is reached
    ! Match the IMAG part
        tol_diff = MINVAL(ABS(AIMAG(M_eig1)-AIMAG(EIG_TARGET)))
        tol_diff = tol_diff/ABS(AIMAG(EIG_TARGET))
        WRITE(*,*) 'Difference:'
        WRITE(*,*) tol_diff
        IF (tol_diff.LE.tol) THEN
            WRITE(*,*) 'Target Reached!'
            GOTO 20
        ELSE
            K_VAL0 = K_VAL1

            ! To be safe, set everything to 0.D0
            M_mat0 = 0.D0
            M_eig0 = 0.D0
            EIG_R_mat0 = 0.D0
            EIG_L_mat0 = 0.D0
            M_mat1 = 0.D0
            M_eig1 = 0.D0
            EIG_R_mat1 = 0.D0
            EIG_L_mat1 = 0.D0
            K_tune = 0.D0
        ENDIF

    ENDDO

    ! Sort and Print eigvals
    ! CALL SORT(0,M_eig1_ND) ! Predicted value
    ! CALL MSAVE(M_eig1_ND,'M_eig1_ND.dat')
    CALL EIG_MATRIX(M_IND2,2,K_VAL1,M_mat1,M_eig1,EIG_R_mat1,EIG_L_mat1)
20  CALL SORT(0,M_eig1)
    CALL MSAVE(M_eig1,'M_eig1.dat')
    WRITE(*,*) 'Final AK value:'
    WRITE(*,*) K_VAL1


    ! DEALLOCATE(K_tune)
    ! DEALLOCATE(M_eig)
    ! DEALLOCATE(M_eig0, M_eig1, M_eig1_ND)
    ! DEALLOCATE(M_mat)      !, EIG_R_mat, EIG_L_mat
    ! DEALLOCATE(M_mat0, EIG_R_mat0, EIG_L_mat0)
    ! DEALLOCATE(M_mat1, EIG_R_mat1, EIG_L_mat1)

10    WRITE(*,*) ''
    WRITE(*,*) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================
SUBROUTINE EIG_MATRIX(M_INDEX,K_INDEX,AKREAD,H,EIG_VAL,EIG_VEC_R,EIG_VEC_L)

    IMPLICIT NONE
    INTEGER,INTENT(IN)    :: M_INDEX,K_INDEX
    REAL(P8),INTENT(INOUT)    :: AKREAD
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: H
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT),OPTIONAL:: EIG_VAL
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT),OPTIONAL:: EIG_VEC_L,EIG_VEC_R

    INTEGER     :: I,MREAD,CHOP_N,NN
    REAL(P8)    :: REI
    CHARACTER*72:: STR1, STR2, STR3, STR4    
    TYPE(SCALAR):: PSI0,CHI0,RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
    TYPE(SCALAR):: PSIU,CHIU,RURU,RUPU,UZU,RORU,ROPU,OZU
    TYPE(SCALAR):: PSI1,CHI1,PSI2,CHI2,PSI3,CHI3


    ! ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
    ! CALL READIN(5)           ! @ MOD_INIT
    ! CALL execute_command_line('echo "%read.input"')

    CALL ALLOCATE(PSI0)
    CALL ALLOCATE(CHI0)
    CALL ALLOCATE(RUR0)
    CALL ALLOCATE(RUP0)
    CALL ALLOCATE( UZ0)
    CALL ALLOCATE(ROR0)
    CALL ALLOCATE(ROP0)
    CALL ALLOCATE( OZ0)

    ! LOAD INITIAL PSI (PSI0) AND CHI (CHI0) DATA
    ! PERFORM 'INIT' AS A PREREQUISITE TO MAKREADE PSI0 AND CHI0
    CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI0,PSI0)
    CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI0,CHI0)

    CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)
    CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0)

    MREAD = M(M_INDEX)
    AKREAD = AK(M_INDEX,K_INDEX)
    CHOP_N = 3

    NN = NRCHOPS(M_INDEX)
    IF (.NOT. ALLOCATED(H))         ALLOCATE(H(2*NN,2*NN))
    H = 0.D0

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
    DO I = 1,NN
        PSIU%LN = 0.D0 
        CHIU%LN = 0.D0 ! DEL2CHI
        PSIU%E = 0.D0
        CHIU%E = 0.D0 ! DEL2CHI
        PSIU%E(I,M_INDEX,K_INDEX) = 1.D0

        CALL CHOPSET(CHOP_N)

        ! Del2[Psi,Del2Chi]/Re term: PSI
        IF (VISC%SW.EQ.1) THEN
        IF(VISC%NU.EQ.0.0) THEN
            WRITE(*,*) 'VISC: VISC%NU CANNOT BE ZERO'
            STOP
        ENDIF
        REI = VISC%NU
        CALL DEL2(PSIU,PSI3)
        CALL DEL2(CHIU,CHI3) ! DEL2[DEL2CHI]
        PSI3%E = PSI3%E * VISC%NU ! NU = 1/RE
        CHI3%E = CHI3%E * VISC%NU ! DEL2[DEL2CHI]/RE
        ELSEIF (VISC%SW.EQ.2) THEN
        REI = VISC%NUP
            CALL HELMP(VISC%P,PSIU,PSI3,0.D0)
            CALL HELMP(VISC%P,CHIU,CHI3,0.D0)
            PSI3%E = PSI3%E * VISC%NUP ! NU = 1/RE
            CHI3%E = CHI3%E * VISC%NUP ! DEL2[DEL2CHI]/RE
        ELSE
        REI = 0.D0
        PSI3%E = 0.D0
        CHI3%E = 0.D0          
        PSI3%LN = 0.D0
        CHI3%LN = 0.D0
        ENDIF

        ! -W x U' + U X W' term: PSI
        ! NOTE: For this test vecotr, DEL2CHI = 0, so CHI = 0.
        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU) ! CHI (= 0 = DEL2CHI)
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU) ! CHI (= 0 = DEL2CHI)
        
        ! CALL VEL2VOR(RURU,RUPU,UZU,RORU,ROPU,OZU)
        
        CALL VPROD_PFF(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        ! CALL TOFP(RURU)
        ! CALL TOFP(RUPU)
        ! CALL TOFP(UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)

        CALL VPROD_PFF(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        ! CALL TOFP(RORU)
        ! CALL TOFP(ROPU)
        ! CALL TOFP(OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-CHOP_N)

        PSI1%E = -PSI1%E + PSI2%E + PSI3%E
        CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NN,  I) = PSI1%E(:NN,M_INDEX,K_INDEX)
        H(NN+1:,I) = CHI1%E(:NN,M_INDEX,K_INDEX) ! DEL2CHI
        
        PSIU%LN = 0.D0
        CHIU%LN = 0.D0
        PSIU%E = 0.D0
        CHIU%E = 0.D0
        CHIU%E(I,M_INDEX,K_INDEX) = 1.D0 ! DEL2CHI

        ! Del2[Psi,Del2Chi]/Re term: CHI
        CALL CHOPSET(CHOP_N)
        IF (VISC%SW.EQ.1) THEN
        CALL DEL2(PSIU,PSI3)
        CALL DEL2(CHIU,CHI3) ! DEL2[DEL2CHI]
        PSI3%E = PSI3%E * VISC%NU ! NU = 1/RE
        CHI3%E = CHI3%E * VISC%NU ! DEL2[DEL2CHI]/RE
        ELSEIF (VISC%SW.EQ.2) THEN
            CALL HELMP(VISC%P,PSIU,PSI3,0.D0)
            CALL HELMP(VISC%P,CHIU,CHI3,0.D0)
            PSI3%E = PSI3%E * VISC%NUP ! NU = 1/RE
            CHI3%E = CHI3%E * VISC%NUP ! DEL2[DEL2CHI]/RE
        ELSE
        PSI3%E = 0.D0
        CHI3%E = 0.D0          
        PSI3%LN = 0.D0
        CHI3%LN = 0.D0
        ENDIF
        
        CALL CHOPSET(-CHOP_N)
        CALL IDEL2(CHIU,CHIU) ! DEL2CHI -> CHI so that PC2VOR and PC2VEL calcualte the correct quantities  
        CALL CHOPSET(CHOP_N)

        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU) ! CHI
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU) ! CHI

        ! CALL VEL2VOR(RURU,RUPU,UZU,RORU,ROPU,OZU)

        CALL VPROD_PFF(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        ! CALL TOFP(RURU)
        ! CALL TOFP(RUPU)
        ! CALL TOFP(UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)

        CALL VPROD_PFF(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        ! CALL TOFP(RORU)
        ! CALL TOFP(ROPU)
        ! CALL TOFP(OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-CHOP_N)

        PSI1%E = -PSI1%E + PSI2%E + PSI3%E
        CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NN,  I+NN) &
                            = PSI1%E(:NN,M_INDEX,K_INDEX)
        H(NN+1:,I+NN) &
                            = CHI1%E(:NN,M_INDEX,K_INDEX) ! DEL2CHI

    ENDDO

! OUTPUT:
WRITE( STR1, '(f10.2)' )  AKREAD  
WRITE( STR2, '(f10.2)' )  ELL
WRITE( STR3, '(f10.6)' )  REI
WRITE( STR4, '(f10.3)' )  QPAIR%H(1)

! CALL MSAVE(H, TRIM(ADJUSTL(FILES%SAVEDIR))// 'mtrx_r_m_' // ITOA4(MREAD) //&
!         '_k_'       // TRIM(ADJUSTL(STR1))  //&
!         '_L_'       // TRIM(ADJUSTL(STR2))  //&
!         '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
!         '_QI_'      // TRIM(ADJUSTL(STR4))  //&
!         '_N_'       // ITOA4(NRCHOPS(M_INDEX)) // '.dat')

IF (PRESENT(EIG_VAL)) THEN    
    IF (.NOT. ALLOCATED(EIG_VAL))   ALLOCATE(EIG_VAL(2*NN))
    EIG_VAL = 0.D0
    EIG_VAL = GENEIG(H)
    ! CALL MSAVE(EIG_VAL, TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigval_m_'// ITOA4(MREAD) //&
    !         '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !         '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !         '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !         '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !         '_N_'       // ITOA4(NRCHOPS(M_INDEX)) // '.dat')
ENDIF

IF (PRESENT(EIG_VEC_R)) THEN
    IF (.NOT. ALLOCATED(EIG_VEC_R)) ALLOCATE(EIG_VEC_R(2*NN,2*NN))
    EIG_VEC_R = 0.D0
    EIG_VEC_R = EIGVEC('R',H)
    ! CALL MSAVE(EIG_VEC_R, TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_FFF_r_m_'// ITOA4(MREAD) //&
    !     '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !     '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !     '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !     '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !     '_N_'       // ITOA4(NRCHOPS(M_INDEX)) // '.dat')
ENDIF

IF (PRESENT(EIG_VEC_L)) THEN
    IF (.NOT. ALLOCATED(EIG_VEC_L)) ALLOCATE(EIG_VEC_L(2*NN,2*NN))
    EIG_VEC_L = 0.D0
    EIG_VEC_L = EIGVEC('L',H)
    ! CALL MSAVE(EIG_VEC_L, TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_FFF_l_m_'// ITOA4(MREAD) //&
    !     '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !     '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !     '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !     '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !     '_N_'       // ITOA4(NRCHOPS(M_INDEX)) // '.dat')
ENDIF

! DEALLOCATE:
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

    RETURN

END SUBROUTINE EIG_MATRIX
!=======================================================================
END PROGRAM EVP_K_TUNE
!=======================================================================

! ! NOTE1: Adjusting Kmode#
! K_IND2 = 1 + (K_shift-1) ! actual K modenumber (not index)
! WRITE(*,*) ''
! WRITE(*,*) 'Nearest K#: '
! WRITE(*,*) K_IND2
! ! Check if K modenumber is allowed
! IF (ABS(K_IND2).GT.(NXCHOP-1)) THEN
!     WRITE(*,*) ''
!     WRITE(*,*) 'K out of range'
!     GOTO 10
! ENDIF
! IF (K_IND2 .LT. 0) THEN
!     K_IND2 = NXCHOPDIM - K_IND2 + 1 ! find index for negative K modenumber
! ELSE
!     K_IND2 = K_IND2 + 1
! ENDIF