!=======================================================================
!
!     WRITTEN BY SANGJOON (JOON) LEE
!     MODIFIED BY JINGE WANG
!     DEPT. OF MECHANICAL ENGINEERING
!     UNIV. OF CALIFORNIA AT BERKELEY
!     EMAIL: SANGJOONLEE@BERKELEY.EDU
!     EMAIL: JINGE@BERKELEY.EDU
!     
!     UC BERKELEY CFD LAB
!     HTTPS://CFD.ME.BERKELEY.EDU/
!
!     NONCOMMERCIAL USE WITH COPYRIGHTED (C) MARK
!     UNDER DEVELOPMENT FOR RESEARCH PURPOSE 
!
!=======================================================================
PROGRAM EVP_LINEAR_PRIMITIVE
!=======================================================================
! [USAGE]: 
! CURRENTLY INVISCID
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN TERMS OF PRIMITIVE VARIABLES (U & P)
! [UPDATES]:
! LAST UPDATE ON MAR 28, 2021
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
    TYPE(SCALAR):: PSI0,CHI0,PSI1,CHI1
    TYPE(SCALAR):: RUR0,RUP0,UZ0,ROR0,ROP0,OZ0 ! Q-vortex
    TYPE(SCALAR):: RURU,RUPU,UZU,RORU,ROPU,OZU ! Test-function
    TYPE(SCALAR):: DEL2R,DEL2P,DEL2Z           ! DEL2[U] terms
    INTEGER     :: MIND,AKIND,MREAD
    REAL(P8)    :: AKREAD, REI
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: H ! Assembled NS operator
    COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: EIG_VEC_3D ! Assembled EigVec matrix of (rur,rup,uz) in PPP space 

    INTEGER     :: I, J, K, IS
    INTEGER     :: N ! H: 3N x 3N
    CHARACTER*72:: STR1, STR2, STR3, STR4

    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC

    ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
    CALL READIN(5)           ! @ MOD_INIT

    CALL ALLOCATE(PSI0)
    CALL ALLOCATE(CHI0)
    CALL ALLOCATE(PSI1)
    CALL ALLOCATE(CHI1)
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

    ! DEALLOCATE T-P FIELDS SINCE ONLY PRIMITIVE VARIABLES ARE USED
    CALL DEALLOCATE(PSI0)
    CALL DEALLOCATE(CHI0)

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

    MREAD = M(MIND)
    AKREAD = AK(MIND,AKIND)
    N = NRCHOPS(MIND)

    WRITE(*,*) 'AZIMUTHAL AND AXIAL WAVENUMBERS:'
    WRITE(*,106) MREAD, AKREAD
106  FORMAT('M = ',I8,' & K = ',F8.3) 

    ALLOCATE( H(3*N,3*N))
    ALLOCATE(EIG_VEC_3D(NR,3*N,3))
    H = 0.D0

    CALL ALLOCATE(RURU)
    CALL ALLOCATE(RUPU)
    CALL ALLOCATE( UZU)
    CALL ALLOCATE(RORU)
    CALL ALLOCATE(ROPU)
    CALL ALLOCATE( OZU)

    CALL ALLOCATE(DEL2R)
    CALL ALLOCATE(DEL2P)
    CALL ALLOCATE(DEL2Z)

    DO I = 1,N
        RURU%LN = 0.D0
        RUPU%LN = 0.D0
        UZU%LN  = 0.D0
        RURU%E  = 0.D0
        RUPU%E  = 0.D0
        UZU%E   = 0.D0

        RORU%LN = 0.D0
        ROPU%LN = 0.D0
        OZU%LN  = 0.D0
        RORU%E  = 0.D0
        ROPU%E  = 0.D0
        OZU%E   = 0.D0

        RURU%E(I,MIND,AKIND) = 1.D0

        CALL CHOPSET(3)
        CALL VEL2VOR(RURU,RUPU,UZU,RORU,ROPU,OZU)
        ! RORU = RURU
        ! ROPU = RUPU
        !  OZU =  UZU
        ! CALL RTRAN(RORU,1)
        ! CALL RTRAN(ROPU,1)
        ! CALL RTRAN( OZU,1)
        ! CALL PROJECT(RORU,ROPU,OZU,PSI1,CHI1)
        ! CALL IDEL2(CHI1,CHI1)
        ! CALL PC2VOR(PSI1,CHI1,RORU,ROPU,OZU)

        ! NU*DEL2[U]
        IF (VISC%SW.EQ.1) THEN
            CALL DEL2_PRIM(RURU,RUPU,UZU,DEL2R,DEL2P,DEL2Z)
            DEL2R%E = DEL2R%E*VISC%NU
            DEL2P%E = DEL2P%E*VISC%NU
            DEL2Z%E = DEL2Z%E*VISC%NU
        ELSE
            DEL2R%LN = 0.D0
            DEL2R%E  = 0.D0
            DEL2P%LN = 0.D0
            DEL2P%E  = 0.D0
            DEL2Z%LN = 0.D0
            DEL2Z%E  = 0.D0
        ENDIF

        ! U X W'
        CALL VPROD_PFF(RUR0,RUP0,UZ0,RORU,ROPU,OZU) ! (RORU,ROPU,OZ) = (RUR0,RUP0,UZ0) x (RORU,ROPU,OZ)
        CALL TOFF(RORU)
        CALL TOFF(ROPU)
        CALL TOFF(OZU)
        
        ! W X U'
        CALL VPROD_PFF(ROR0,ROP0,OZ0,RURU,RUPU,UZU) ! (RURU,RUPU,UZ) = (ROR0,ROP0,OZ0) x (RURU,RUPU,UZ)
        CALL TOFF(RURU)
        CALL TOFF(RUPU)
        CALL TOFF(UZU)

        CALL CHOPSET(-3)

        ! U X W' - W X U'
        H(    1:  N, I) = RORU%E(1:N,MIND,AKIND) - RURU%E(1:N,MIND,AKIND) + DEL2R%E(1:N,MIND,AKIND)
        H(  N+1:2*N, I) = ROPU%E(1:N,MIND,AKIND) - RUPU%E(1:N,MIND,AKIND) + DEL2P%E(1:N,MIND,AKIND)
        H(2*N+1:3*N, I) =  OZU%E(1:N,MIND,AKIND) -  UZU%E(1:N,MIND,AKIND) + DEL2Z%E(1:N,MIND,AKIND)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        RURU%LN = 0.D0
        RUPU%LN = 0.D0
        UZU%LN  = 0.D0
        RURU%E  = 0.D0
        RUPU%E  = 0.D0
        UZU%E   = 0.D0

        RORU%LN = 0.D0
        ROPU%LN = 0.D0
        OZU%LN  = 0.D0
        RORU%E  = 0.D0
        ROPU%E  = 0.D0
        OZU%E   = 0.D0

        RUPU%E(I,MIND,AKIND) = 1.D0

        CALL CHOPSET(3)
        CALL VEL2VOR(RURU,RUPU,UZU,RORU,ROPU,OZU)
        ! RORU = RURU
        ! ROPU = RUPU
        !  OZU =  UZU
        ! CALL RTRAN(RORU,1)
        ! CALL RTRAN(ROPU,1)
        ! CALL RTRAN( OZU,1)
        ! CALL PROJECT(RORU,ROPU,OZU,PSI1,CHI1)
        ! CALL IDEL2(CHI1,CHI1)
        ! CALL PC2VOR(PSI1,CHI1,RORU,ROPU,OZU)

        ! NU*DEL2[U]
        IF (VISC%SW.EQ.1) THEN
            CALL DEL2_PRIM(RURU,RUPU,UZU,DEL2R,DEL2P,DEL2Z)
            DEL2R%E = DEL2R%E*VISC%NU
            DEL2P%E = DEL2P%E*VISC%NU
            DEL2Z%E = DEL2Z%E*VISC%NU
        ELSE
            DEL2R%LN = 0.D0
            DEL2R%E  = 0.D0
            DEL2P%LN = 0.D0
            DEL2P%E  = 0.D0
            DEL2Z%LN = 0.D0
            DEL2Z%E  = 0.D0
        ENDIF

        ! U X W'
        CALL VPROD_PFF(RUR0,RUP0,UZ0,RORU,ROPU,OZU) ! (RORU,ROPU,OZ) = (RUR0,RUP0,UZ0) x (RORU,ROPU,OZ)
        CALL TOFF(RORU)
        CALL TOFF(ROPU)
        CALL TOFF(OZU)
        
        ! W X U'
        CALL VPROD_PFF(ROR0,ROP0,OZ0,RURU,RUPU,UZU) ! (RURU,RUPU,UZ) = (ROR0,ROP0,OZ0) x (RURU,RUPU,UZ)
        CALL TOFF(RURU)
        CALL TOFF(RUPU)
        CALL TOFF(UZU)

        CALL CHOPSET(-3)

        ! U X W' - W X U'
        H(    1:  N, I+N) = RORU%E(1:N,MIND,AKIND) - RURU%E(1:N,MIND,AKIND) + DEL2R%E(1:N,MIND,AKIND)
        H(  N+1:2*N, I+N) = ROPU%E(1:N,MIND,AKIND) - RUPU%E(1:N,MIND,AKIND) + DEL2P%E(1:N,MIND,AKIND)
        H(2*N+1:3*N, I+N) =  OZU%E(1:N,MIND,AKIND) -  UZU%E(1:N,MIND,AKIND) + DEL2Z%E(1:N,MIND,AKIND)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        RURU%LN = 0.D0
        RUPU%LN = 0.D0
        UZU%LN  = 0.D0
        RURU%E  = 0.D0
        RUPU%E  = 0.D0
        UZU%E   = 0.D0

        RORU%LN = 0.D0
        ROPU%LN = 0.D0
        OZU%LN  = 0.D0
        RORU%E  = 0.D0
        ROPU%E  = 0.D0
        OZU%E   = 0.D0
        
        UZU%E(I,MIND,AKIND) = 1.D0

        CALL CHOPSET(3)
        CALL VEL2VOR(RURU,RUPU,UZU,RORU,ROPU,OZU)
        ! RORU = RURU
        ! ROPU = RUPU
        !  OZU =  UZU
        ! CALL RTRAN(RORU,1)
        ! CALL RTRAN(ROPU,1)
        ! CALL RTRAN( OZU,1)
        ! CALL PROJECT(RORU,ROPU,OZU,PSI1,CHI1)
        ! CALL IDEL2(CHI1,CHI1)
        ! CALL PC2VOR(PSI1,CHI1,RORU,ROPU,OZU)

        ! NU*DEL2[U]
        IF (VISC%SW.EQ.1) THEN
            CALL DEL2_PRIM(RURU,RUPU,UZU,DEL2R,DEL2P,DEL2Z)
            DEL2R%E = DEL2R%E*VISC%NU
            DEL2P%E = DEL2P%E*VISC%NU
            DEL2Z%E = DEL2Z%E*VISC%NU
        ELSE
            DEL2R%LN = 0.D0
            DEL2R%E  = 0.D0
            DEL2P%LN = 0.D0
            DEL2P%E  = 0.D0
            DEL2Z%LN = 0.D0
            DEL2Z%E  = 0.D0
        ENDIF

        ! U X W'
        CALL VPROD_PFF(RUR0,RUP0,UZ0,RORU,ROPU,OZU) ! (RORU,ROPU,OZ) = (RUR0,RUP0,UZ0) x (RORU,ROPU,OZ)
        CALL TOFF(RORU)
        CALL TOFF(ROPU)
        CALL TOFF(OZU)
        
        ! W X U'
        CALL VPROD_PFF(ROR0,ROP0,OZ0,RURU,RUPU,UZU) ! (RURU,RUPU,UZ) = (ROR0,ROP0,OZ0) x (RURU,RUPU,UZ)
        CALL TOFF(RURU)
        CALL TOFF(RUPU)
        CALL TOFF(UZU)

        CALL CHOPSET(-3)

        ! U X W' - W X U'
        H(    1:  N, I+2*N) = RORU%E(1:N,MIND,AKIND) - RURU%E(1:N,MIND,AKIND) + DEL2R%E(1:N,MIND,AKIND)
        H(  N+1:2*N, I+2*N) = ROPU%E(1:N,MIND,AKIND) - RUPU%E(1:N,MIND,AKIND) + DEL2P%E(1:N,MIND,AKIND)
        H(2*N+1:3*N, I+2*N) =  OZU%E(1:N,MIND,AKIND) -  UZU%E(1:N,MIND,AKIND) + DEL2Z%E(1:N,MIND,AKIND)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        WRITE(*,107) I, N
    ENDDO
107  FORMAT('PROGRESS: ',I3.3,'/',I3.3)

    ! OUTPUT:
    REI = VISC%NU
    WRITE( STR1, '(f10.2)' )  AKREAD  
    WRITE( STR2, '(f10.2)' )  ELL
    WRITE( STR3, '(f10.6)' )  REI
    WRITE( STR4, '(f10.3)' )  QPAIR%H(1)

    ! CALL MCAT(GENEIG(H))

    WRITE(*,*) "Print output:"
    open(10,FILE= TRIM(ADJUSTL(FILES%SAVEDIR))// 'output_name.dat',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)
    if (IS.ne.0) then
        print *, 'ERROR: evp_linear -- Could not creat new file output_name.dat'
        GOTO 1
    end if
    WRITE(10,*) 'mtrx_r_m_' // ITOA4(MREAD) //&
                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    WRITE(10,*) 'eigval_m_'// ITOA4(MREAD) //&
                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    WRITE(10,*) 'eigvec_FFF_r_m_'// ITOA4(MREAD) //&
                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    WRITE(10,*) 'eigvec_FFF_l_m_'// ITOA4(MREAD) //&
                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    WRITE(10,*) 'eigvec_PFF_rur_m_'// ITOA4(MREAD) //&
                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    WRITE(10,*) 'eigvec_PFF_rup_m_'// ITOA4(MREAD) //&
                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    WRITE(10,*) 'eigvec_PFF__uz_m_'// ITOA4(MREAD) //&
                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    CLOSE(10)

    ! SAVE M_MATRIX, EIGVAL, AND EIGVEC
    CALL MSAVE(H, TRIM(ADJUSTL(FILES%SAVEDIR))// 'mtrx_r_m_' // ITOA4(MREAD) //&
                                                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                                                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                                                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                                                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                                                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat')

    CALL MSAVE(GENEIG(H), TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigval_m_'// ITOA4(MREAD) //&
                                                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                                                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                                                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                                                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                                                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat')

    CALL MSAVE(EIGVEC('R',H), TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_FFF_r_m_'// ITOA4(MREAD) //&
                                                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                                                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                                                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                                                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                                                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat')

    CALL MSAVE(EIGVEC('L',H), TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_FFF_l_m_'// ITOA4(MREAD) //&
                                                '_k_'       // TRIM(ADJUSTL(STR1))  //&
                                                '_L_'       // TRIM(ADJUSTL(STR2))  //&
                                                '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                                                '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                                                '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat')

    ! CREATE AND SAVE EIGVEC IN PFF space:
    H = EIGVEC('R', H)
    EIG_VEC_3D = 0.D0
    DO I = 1,3*N
        CALL TOFF(RURU)
        CALL TOFF(RUPU)
        CALL TOFF( UZU)
        RURU%LN = 0.D0
        RUPU%LN = 0.D0
        UZU%LN  = 0.D0
        RURU%E  = 0.D0
        RUPU%E  = 0.D0
        UZU%E   = 0.D0

        RURU%E(:N,MIND,AKIND) =        H(:N,I)
        RUPU%E(:N,MIND,AKIND) =   H(N+1:2*N,I)
         UZU%E(:N,MIND,AKIND) = H(2*N+1:3*N,I)

        CALL RTRAN(RURU,1)
        CALL RTRAN(RUPU,1)
        CALL RTRAN( UZU,1)

        EIG_VEC_3D(:NR,I,1) = RURU%E(:N,MIND,AKIND)
        EIG_VEC_3D(:NR,I,2) = RUPU%E(:NR,MIND,AKIND)
        EIG_VEC_3D(:NR,I,3) = UZU%E (:NR,MIND,AKIND)
    ENDDO
    CALL MSAVE(EIG_VEC_3D(:NR,:,1), TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_PFF_rur_m_'// ITOA4(MREAD) //&
                                                    '_k_'       // TRIM(ADJUSTL(STR1))  //&
                                                    '_L_'       // TRIM(ADJUSTL(STR2))  //&
                                                    '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                                                    '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                                                    '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat')
    CALL MSAVE(EIG_VEC_3D(:NR,:,2), TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_PFF_rup_m_'// ITOA4(MREAD) //&
                                                    '_k_'       // TRIM(ADJUSTL(STR1))  //&
                                                    '_L_'       // TRIM(ADJUSTL(STR2))  //&
                                                    '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                                                    '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                                                    '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat')
    CALL MSAVE(EIG_VEC_3D(:NR,:,3), TRIM(ADJUSTL(FILES%SAVEDIR))// 'eigvec_PFF__uz_m_'// ITOA4(MREAD) //&
                                                    '_k_'       // TRIM(ADJUSTL(STR1))  //&
                                                    '_L_'       // TRIM(ADJUSTL(STR2))  //&
                                                    '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
                                                    '_QI_'      // TRIM(ADJUSTL(STR4))  //&
                                                    '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat')


    ! open(10,FILE= TRIM(ADJUSTL(FILES%SAVEDIR))// 'output_name.dat',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)
    ! if (IS.ne.0) then
    !     print *, 'ERROR: evp_linear -- Could not creat new file output_name.dat'
    !     GOTO 1
    ! end if
    ! WRITE(10,*) 'mtrx_r_m_' // ITOA4(MREAD) //&
    !             '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !             '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !             '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !             '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !             '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    ! WRITE(10,*) 'eigval_m_'// ITOA4(MREAD) //&
    !             '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !             '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !             '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !             '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !             '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    ! WRITE(10,*) 'eigvec_FFF_r_m_'// ITOA4(MREAD) //&
    !             '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !             '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !             '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !             '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !             '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    ! WRITE(10,*) 'eigvec_FFF_l_m_'// ITOA4(MREAD) //&
    !             '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !             '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !             '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !             '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !             '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    ! WRITE(10,*) 'eigvec_PFF_rur_m_'// ITOA4(MREAD) //&
    !             '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !             '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !             '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !             '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !             '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    ! WRITE(10,*) 'eigvec_PFF_rup_m_'// ITOA4(MREAD) //&
    !             '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !             '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !             '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !             '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !             '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    ! WRITE(10,*) 'eigvec_PFF__uz_m_'// ITOA4(MREAD) //&
    !             '_k_'       // TRIM(ADJUSTL(STR1))  //&
    !             '_L_'       // TRIM(ADJUSTL(STR2))  //&
    !             '_ReI_'     // TRIM(ADJUSTL(STR3))  //&
    !             '_QI_'      // TRIM(ADJUSTL(STR4))  //&
    !             '_N_'       // ITOA4(NRCHOPS(MIND)) // '.dat'
    ! CLOSE(10)

1    DEALLOCATE( H )
    DEALLOCATE(EIG_VEC_3D)
    
    CALL DEALLOCATE(DEL2R)
    CALL DEALLOCATE(DEL2P)
    CALL DEALLOCATE(DEL2Z)

    CALL DEALLOCATE(RURU)
    CALL DEALLOCATE(RUPU)
    CALL DEALLOCATE( UZU)
    CALL DEALLOCATE(RORU)
    CALL DEALLOCATE(ROPU)
    CALL DEALLOCATE( OZU)

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
SUBROUTINE DEL2_PRIM(RUR,RUP,UZ,D2R,D2P,D2Z)
    IMPLICIT NONE
    TYPE(SCALAR),INTENT(IN)::RUR,RUP,UZ
    TYPE(SCALAR),INTENT(INOUT):: D2R,D2P,D2Z
    TYPE(SCALAR)::W1,W2,W3,W4,W5

    CALL ALLOCATE(W1)
    CALL ALLOCATE(W2)
    CALL ALLOCATE(W3)
    CALL ALLOCATE(W4)
    CALL ALLOCATE(W5)

    W1 = RUR
    W2 = RUP
    W3 =  UZ
    CALL RTRAN(W1,1)
    CALL RTRAN(W2,1)
    CALL RTRAN(W3,1)

    CALL PROJECT(W1,W2,W3,W4,W5)
    CALL DEL2(W4,W1)
    ! CALL DEL2(W5,W2) W5 from PROJECT is already DEL2CHI 

    CALL PC2VEL(W1,W5,D2R,D2P,D2Z)

    RETURN
END SUBROUTINE DEL2_PRIM

!=======================================================================
END PROGRAM EVP_LINEAR_PRIMITIVE
!=======================================================================
