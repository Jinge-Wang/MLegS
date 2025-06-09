PROGRAM EVP_STRAIN_MATRIX
!=======================================================================
! [USAGE]:
! USE VPROD_PFF RATHER THAN VPROD
!
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
!
! [UPDATES]:
! LAST UPDATE ON DEC 3, 2020
!=======================================================================
USE omp_lib
USE MPI
USE MOD_MISC
USE MOD_BANDMAT
USE MOD_EIG
USE MOD_FD
USE MOD_LIN_LEGENDRE
USE MOD_SCALAR3
USE MOD_FFT
USE MOD_LAYOUT
USE MOD_LEGOPS
USE MOD_MARCH
USE MOD_INIT
USE MOD_EVP
USE MOD_DIAGNOSTICS

IMPLICIT NONE
INTEGER     :: II, JJ, KK, K_ind
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_R
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_L

! MPI:
INTEGER:: MPI_GLB_PROCS, MPI_GLB_RANK, newcomm
! ADD PERTURB:
INTEGER:: IS, NDIM
! SCANNING:
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID

! CALL MPI_INIT(IERR)
CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, MPI_THREAD_MODE, IERR)
IF (MPI_THREAD_MODE .LT. MPI_THREAD_SERIALIZED) THEN
    WRITE (*, *) 'The threading support is lesser than that demanded.'
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
END IF

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
END IF
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS
CALL READCOM('NOECHO')
CALL READIN(5)

! CALL STRAIN_MATRIX_F(M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm) 

! IF (MPI_GLB_RANK .EQ. 0) THEN

!     ALLOCATE(EIG_R(NRCHOPDIM))

!     open(FID,FILE='./converg/check4/r_info.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)
!     DO II = 1,NR
!         WRITE(FID,'(G20.12)') TFM%R(II)
!     ENDDO

!     II = MINLOC(ABS(M_eig),1)
!     EIG_R(1:SIZE(M_eig)) = EIG_R_mat(:,II)
!     WRITE(*,*) 'II = ',II,':',M_eig(II)
!     WRITE(*,*) EIG_R
!     CALL RTRAN_MK(EIG_R,1)

!     WRITE(*,*) 'PHYSICAL:'
!     WRITE(*,*) EIG_R
!     call SAVE_VEC('./converg/check4/Strain_Vec_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R)
!     DEALLOCATE(EIG_R,M_eig)

! ENDIF

! CALL STRAIN_FIELD1(newcomm)
WRITE(*,*) STRAIN_FIELD2(newcomm)

IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================

SUBROUTINE STRAIN_MATRIX_P(H, EIG_VAL, EIG_VEC_R, EIG_VEC_L, comm_grp)
! ======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: comm_grp
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT)            :: H
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT), OPTIONAL     :: EIG_VAL
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT), OPTIONAL  :: EIG_VEC_R, EIG_VEC_L

    INTEGER     :: NR_MK, I, M_ACTUAL
    REAL(P8)    :: REI, AK_ACTUAL
    INTEGER     :: MPI_NR_SIZE, MPI_NR_INDEX
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: TU, T1, T2

    REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XXDXOP

    ! INIT
    NTH = 2
    NX = 4
    NTCHOP = 2 ! M = 0, whatever
    NXCHOP = 2 ! K = 0, whatever
    CALL LEGINIT(comm_grp,0)

    ! FORM RD()/DR OPERATOR FOR M = 0:
    ALLOCATE( XXDXOP(NRCHOP,3) )
    NR_MK = NRCHOPS(1)
    XXDXOP(:NR_MK,:) = BAND_LEG_XXDX(NR_MK,M(1),TFM%NORM(:,1))
    ! B(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,A(:NR_MK))

    ! IF (MPI_GLB_RANK.EQ.0) THEN
    !     ! CALL MCAT(TFM%PF(:NRH,1:NR_MK:2,1)-TFM%PF(:NRH,1:NR_MK:2,2))
    !     write(*,*) sum(abs(TFM%PF(:NRH,1:NR_MK:2,1)-TFM%PF(:NRH,1:NR_MK:2,2)))
    !     CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    ! ENDIF

    ! SPLIT NR
    IF (ALLOCATED(H) .AND. (SIZE(H, 1) .NE. NR)) THEN
        DEALLOCATE (H)
        ALLOCATE (H(NR,NR))
    ELSEIF (.NOT. (ALLOCATED(H))) THEN
        ALLOCATE (H(NR,NR))
    END IF
    H = CMPLX(0.D0, 0.D0, P8)
    CALL DECOMPOSE(NR, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)

    ! MAIN JOB
    DO I = MPI_NR_INDEX + 1, MPI_NR_INDEX + MPI_NR_SIZE

        ALLOCATE (TU(NR),T1(NRCHOPDIM),T2(NRCHOPDIM))
        TU = 0.D0; T1 = 0.D0; T2 = 0.D0

        ! TEST FUNCTION:
        TU(I) = 1.D0
        CALL RTRAN_MK(TU,-1)

        ! OPERATOR:
        T1(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,TU(:NR_MK)) ! T1 = R*D[TU]/DR
        T2(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,T1(:NR_MK)) ! T2 = R*D[R*D[TU]/DR]/DR

        ! BACK TO PHYSICAL
        CALL RTRAN_MK(T1, 1)
        CALL RTRAN_MK(T2, 1)

        ! FILL IN MATRIX
        H(:, I) = T2+4*T1
        H(I, I) = H(I, I) + 4*(TFM%R(I))**4/(EXP((TFM%R(I))**2)-1)

        ! DEALLOCATE
        DEALLOCATE (TU,T1,T2)

    END DO
    DEALLOCATE(XXDXOP)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, H, SIZE(H), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                        MPI_COMM_WORLD, IERR)

    ! OUTPUT:
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    IF (MPI_GLB_RANK .EQ. 0) THEN  ! START OF THE SERIAL PART
    !WRITE OUT THE 3D MATRIX IN THE SCALAR USING UNFORMATTED MODE

        IF (PRESENT(EIG_VAL)) THEN
            IF (ALLOCATED(EIG_VAL)) DEALLOCATE (EIG_VAL)
            ALLOCATE (EIG_VAL(NR))

            IF (PRESENT(EIG_VEC_L)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE (EIG_VEC_R)
                IF (ALLOCATED(EIG_VEC_L)) DEALLOCATE (EIG_VEC_L)
                ALLOCATE (EIG_VEC_R(NR, NR), EIG_VEC_L(NR, NR))
                CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R, EL=EIG_VEC_L)
            ELSEIF (PRESENT(EIG_VEC_R)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE (EIG_VEC_R)
                ALLOCATE (EIG_VEC_R(NR, NR))
                CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R)
            ELSE
                CALL EIGENDECOMPOSE(H, EIG_VAL)
            END IF
        END IF
    END IF

    RETURN

END SUBROUTINE STRAIN_MATRIX_P
!=======================================================================

SUBROUTINE STRAIN_MATRIX_F(H, EIG_VAL, EIG_VEC_R, EIG_VEC_L, comm_grp)
! ======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: comm_grp
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT)            :: H
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT), OPTIONAL     :: EIG_VAL
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT), OPTIONAL  :: EIG_VEC_R, EIG_VEC_L

    INTEGER     :: NR_MK, I, M_ACTUAL
    REAL(P8)    :: REI, AK_ACTUAL
    INTEGER     :: MPI_NR_SIZE, MPI_NR_INDEX
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: TU, T1, T2, T3

    REAL(P8),DIMENSION(:),ALLOCATABLE:: func
    REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XXDXOP

    ! INIT
    NTH = 2
    NX = 2
    NTCHOP = 2 ! M = 0, whatever
    NXCHOP = 1 ! K = 0, whatever
    CALL LEGINIT(comm_grp,0)

    ! FORM Wz'/U_theta:
    ALLOCATE(func(NR))
    func = 4*(TFM%R)**4/(EXP((TFM%R)**2)-1)

    ! FORM RD()/DR OPERATOR FOR M = 0:
    ALLOCATE( XXDXOP(NRCHOP,3) )
    NR_MK = NRCHOPS(1)
    XXDXOP(:NR_MK,:) = BAND_LEG_XXDX(NR_MK,M(1),TFM%NORM(:,1))
    ! B(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,A(:NR_MK))

    IF (MPI_GLB_RANK.EQ.0) THEN
        ! ! CALL MCAT(TFM%PF(:NRH,1:NR_MK:2,1)-TFM%PF(:NRH,1:NR_MK:2,2))
        ! write(*,*) sum(abs(TFM%PF(:NRH,1:NR_MK:2,1)-TFM%PF(:NRH,1:NR_MK:2,2)))
        ! CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)

        write(*,*) NR_MK
    ENDIF

    ! SPLIT NR_MK
    IF (ALLOCATED(H) .AND. (SIZE(H, 1) .NE. NR_MK)) THEN
        DEALLOCATE (H)
        ALLOCATE (H(NR_MK,NR_MK))
    ELSEIF (.NOT. (ALLOCATED(H))) THEN
        ALLOCATE (H(NR_MK,NR_MK))
    END IF
    H = CMPLX(0.D0, 0.D0, P8)
    CALL DECOMPOSE(NR_MK, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)

    ! MAIN JOB
    DO I = MPI_NR_INDEX + 1, MPI_NR_INDEX + MPI_NR_SIZE

        ALLOCATE (TU(NRCHOPDIM),T1(NRCHOPDIM),T2(NRCHOPDIM),T3(NDIMR))
        TU = 0.D0; T1 = 0.D0; T2 = 0.D0; T3 = 0.D0

        ! TEST FUNCTION:
        TU(I) = 1.D0

        ! OPERATOR:
        T1(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,TU(:NR_MK)) ! T1 = R*D[TU]/DR
        T2(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,T1(:NR_MK)) ! T2 = R*D[R*D[TU]/DR]/DR

        ! BACK TO PHYSICAL
        CALL RTRAN_MK(TU, 1)
        T3(:NR) = FUNC(:NR)*TU(:NR)
        CALL RTRAN_MK(T3,-1)

        ! FILL IN MATRIX
        H(:, I) = T2+4*T1+T3

        ! DEALLOCATE
        DEALLOCATE (TU,T1,T2,T3)

    END DO
    DEALLOCATE(XXDXOP,FUNC)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, H, SIZE(H), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                        MPI_COMM_WORLD, IERR)

    ! OUTPUT:
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    IF (MPI_GLB_RANK .EQ. 0) THEN  ! START OF THE SERIAL PART
    !WRITE OUT THE 3D MATRIX IN THE SCALAR USING UNFORMATTED MODE

        IF (PRESENT(EIG_VAL)) THEN
            IF (ALLOCATED(EIG_VAL)) DEALLOCATE (EIG_VAL)
            ALLOCATE (EIG_VAL(NR_MK))

            IF (PRESENT(EIG_VEC_L)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE (EIG_VEC_R)
                IF (ALLOCATED(EIG_VEC_L)) DEALLOCATE (EIG_VEC_L)
                ALLOCATE (EIG_VEC_R(NR_MK, NR_MK), EIG_VEC_L(NR_MK, NR_MK))
                CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R, EL=EIG_VEC_L)
            ELSEIF (PRESENT(EIG_VEC_R)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE (EIG_VEC_R)
                ALLOCATE (EIG_VEC_R(NR_MK, NR_MK))
                CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R)
            ELSE
                CALL EIGENDECOMPOSE(H, EIG_VAL)
            END IF
        END IF
    END IF

    RETURN

END SUBROUTINE STRAIN_MATRIX_F
!=======================================================================

SUBROUTINE STRAIN_FIELD1(comm_grp)
! ======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN):: comm_grp

    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: EIG_VAL
    COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: H, EIG_VEC_R
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: RUR,RUP,UZ
    
    INTEGER     :: NR_MK, I, M_ACTUAL
    REAL(P8)    :: REI, AK_ACTUAL
    COMPLEX(P8) :: COEFF
    INTEGER     :: MPI_NR_SIZE, MPI_NR_INDEX
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: TU, T1, T2, T3, G

    REAL(P8),DIMENSION(:),ALLOCATABLE:: func
    REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XXDXOP

    ! INIT
    NTH = 2; NX = 2
    NTCHOP = 2 ! M = 0, whatever
    NXCHOP = 1 ! K = 0, whatever
    CALL LEGINIT(comm_grp,0)

    ! FORM Wz'/U_theta:
    ALLOCATE(func(NR))
    func = 4*(TFM%R)**4/(EXP((TFM%R)**2)-1) ! Q-VORTEX

    ! FORM RD()/DR OPERATOR FOR M = 0:
    ALLOCATE( XXDXOP(NRCHOP,3) )
    NR_MK = NRCHOPS(1)
    XXDXOP(:NR_MK,:) = BAND_LEG_XXDX(NR_MK,M(1),TFM%NORM(:,1))
    ! B(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,A(:NR_MK))

    IF (MPI_GLB_RANK.EQ.0) THEN
        write(*,*) NR_MK
    ENDIF

    ALLOCATE(H(NR_MK,NR_MK))
    H = CMPLX(0.D0, 0.D0, P8)

    ! MAIN JOB
    CALL DECOMPOSE(NR_MK, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)
    DO I = MPI_NR_INDEX + 1, MPI_NR_INDEX + MPI_NR_SIZE

        ALLOCATE (TU(NRCHOPDIM),T1(NRCHOPDIM),T2(NRCHOPDIM),T3(NDIMR))
        TU = 0.D0; T1 = 0.D0; T2 = 0.D0; T3 = 0.D0

        ! TEST FUNCTION:
        TU(I) = 1.D0

        ! OPERATOR:
        T1(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,TU(:NR_MK)) ! T1 = R*D[TU]/DR
        T2(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,T1(:NR_MK)) ! T2 = R*D[R*D[TU]/DR]/DR

        ! BACK TO PHYSICAL
        CALL RTRAN_MK(TU, 1)
        T3(:NR) = FUNC(:NR)*TU(:NR)
        CALL RTRAN_MK(T3,-1)

        ! FILL IN MATRIX
        H(:, I) = T2+4*T1+T3

        ! DEALLOCATE
        DEALLOCATE (TU,T1,T2,T3)

    END DO

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, H, SIZE(H), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                        MPI_COMM_WORLD, IERR)
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)                    

    ! SOLVE EVP:
    IF (MPI_GLB_RANK .EQ. 0) THEN

        ALLOCATE(EIG_VAL(NR_MK), EIG_VEC_R(NR_MK, NR_MK))
        CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R)

        ! LOCATE EIG=0
        ALLOCATE(G(NRCHOPDIM),RUR(NRCHOPDIM),RUP(NRCHOPDIM),UZ(NDIMR))
        I = MINLOC(ABS(EIG_VAL),1)
        G(1:SIZE(EIG_VAL)) = EIG_VEC_R(:,I)
        DEALLOCATE(EIG_VAL,EIG_VEC_R,H)

        ! FORM RUR
        RUR = G
        CALL RTRAN_MK(RUR,1)
        COEFF = 1/RUR(NR)
        RUR(1:NR) = IU*(-COEFF*(RUR(1:NR)*(TFM%R(1:NR))**2)/2)

        ! FORM RUP
        RUP = G + 0.5*BANMUL(XXDXOP(:NR_MK,:),2,G(:NR_MK))
        CALL RTRAN_MK(RUP,1)
        RUP(1:NR) = COEFF*(RUP(1:NR)*(TFM%R(1:NR))**2)/4

        ! FORM UZ
        UZ = 0.D0

        ! WRITE(*,*) 'RUR:'
        ! WRITE(*,*) RUR

        ! WRITE(*,*) 'RUP:'
        ! WRITE(*,*) RUP

        ! WRITE(*,*) 'COEFF:'
        ! write(*,*) COEFF
    END IF

    DEALLOCATE(XXDXOP,FUNC)
    RETURN

END SUBROUTINE STRAIN_FIELD1
!=======================================================================

FUNCTION STRAIN_FIELD2(comm_grp)
! ======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN):: comm_grp
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: STRAIN_FIELD2

    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: EIG_VAL
    COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: H, EIG_VEC_R
    
    INTEGER     :: NR_MK, I, M_ACTUAL
    REAL(P8)    :: REI, AK_ACTUAL
    COMPLEX(P8) :: COEFF
    INTEGER     :: MPI_NR_SIZE, MPI_NR_INDEX
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: TU, T1, T2, T3, G

    ! OPERATOR:
    REAL(P8),DIMENSION(:),ALLOCATABLE:: func
    REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XXDXOP

    ! DAMPING:
    REAL(P8):: RC = 4, E = 1.D0, dR = 2
    REAL(P8),DIMENSION(:),ALLOCATABLE:: R

    ! INIT
    NTH = 2; NX = 2
    NTCHOP = 2 ! M = 0, whatever
    NXCHOP = 1 ! K = 0, whatever
    CALL LEGINIT(comm_grp,0)

    ! FORM Wz'/U_theta:
    ALLOCATE(func(NR))
    func = 4*(TFM%R)**4/(EXP((TFM%R)**2)-1) ! Q-VORTEX

    ! FORM RD()/DR OPERATOR FOR M = 0:
    ALLOCATE( XXDXOP(NRCHOP,3) )
    NR_MK = NRCHOPS(1)
    XXDXOP(:NR_MK,:) = BAND_LEG_XXDX(NR_MK,M(1),TFM%NORM(:,1))
    ! B(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,A(:NR_MK))

    IF (MPI_GLB_RANK.EQ.0) THEN
        write(*,*) NR_MK
    ENDIF

    ALLOCATE(H(NR_MK,NR_MK))
    H = CMPLX(0.D0, 0.D0, P8)

    ! MAIN JOB
    CALL DECOMPOSE(NR_MK, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)
    DO I = MPI_NR_INDEX + 1, MPI_NR_INDEX + MPI_NR_SIZE

        ALLOCATE (TU(NRCHOPDIM),T1(NRCHOPDIM),T2(NRCHOPDIM),T3(NDIMR))
        TU = 0.D0; T1 = 0.D0; T2 = 0.D0; T3 = 0.D0

        ! TEST FUNCTION:
        TU(I) = 1.D0

        ! OPERATOR:
        T1(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,TU(:NR_MK)) ! T1 = R*D[TU]/DR
        T2(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,T1(:NR_MK)) ! T2 = R*D[R*D[TU]/DR]/DR

        ! BACK TO PHYSICAL
        CALL RTRAN_MK(TU, 1)
        T3(:NR) = FUNC(:NR)*TU(:NR)
        CALL RTRAN_MK(T3,-1)

        ! FILL IN MATRIX
        H(:, I) = T2+4*T1+T3

        ! DEALLOCATE
        DEALLOCATE (TU,T1,T2,T3)

    END DO

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, H, SIZE(H), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                        MPI_COMM_WORLD, IERR)
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)                    

    ! SOLVE EVP:
    IF (MPI_GLB_RANK .EQ. 0) THEN

        ALLOCATE(EIG_VAL(NR_MK), EIG_VEC_R(NR_MK, NR_MK))
        CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R)

        ! LOCATE EIG=0
        I = MINLOC(ABS(EIG_VAL),1)
        ALLOCATE(G(NRCHOPDIM)); G(1:SIZE(EIG_VAL)) = EIG_VEC_R(:,I)
        CALL RTRAN_MK(G,1); G = G/G(NR)
        WRITE(*,*) G
        DEALLOCATE(EIG_VAL,EIG_VEC_R,H)

        ! FORM STREAM FUNCTION
        ALLOCATE(R(SIZE(TFM%R)),STRAIN_FIELD2(SIZE(TFM%R)))
        R = TFM%R
        STRAIN_FIELD2 = -0.5*E*R**2*G(1:SIZE(R))
        STRAIN_FIELD2 = STRAIN_FIELD2*0.5*(1-TANH((R-RC)/dR))
        DEALLOCATE(G)

    END IF

    DEALLOCATE(XXDXOP,FUNC)
    RETURN

END FUNCTION STRAIN_FIELD2
!=======================================================================

END PROGRAM EVP_STRAIN_MATRIX
!=======================================================================
