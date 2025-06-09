PROGRAM EVP_SCAN_0920
!=======================================================================
! [USAGE]: 
! Faster Version: the bar term is calculated once for each wavenumber
! Modified to include M=0 mode
! Check resolved or not using EIGRES
! 
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
! USE VPROD_PFF
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

IMPLICIT NONE
INTEGER :: II, JJ, KK
INTEGER :: M_1, M_2
INTEGER :: M_BAR, EIG_BAR_IND, tol_ind(2)
INTEGER :: K_BAR_II, K_VAL_1_II, K_VAL_2_II, K_VAL_INDEX, K_VAL_COUNT
REAL(P8):: tol, tol_diff
REAL(P8):: K_BAR, K_BAR_START, K_BAR_END, K_VAL_1_0, K_VAL_2_0
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_BAR_RANGE, K_VAL_1_RANGE, K_VAL_2_RANGE
LOGICAL, DIMENSION(:),ALLOCATABLE:: MASK_1, MASK_2

COMPLEX(P8):: EIG_BAR, EIG_PRIME, EIG_PRIME_R, ALPHA, BETA
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_1_0
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig_2_0
LOGICAL,    DIMENSION(:,:),ALLOCATABLE:: M_res_1_mat, M_res_2_mat
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_eig_1_mat, M_eig_2_mat
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat_2_0, EIG_R_mat_2_0, EIG_L_mat_2_0

COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_BAR,VECR_BAR,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_1, VECR_1, VECR_B1, RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: VECL_2, VECR_2, VECR_B2, RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2

INTEGER:: IS
CHARACTER(LEN=72):: FILENAME
CHARACTER(LEN=300):: PRTSTR
CHARACTER(:),ALLOCATABLE:: LONGSTR
INTEGER:: LOGID, K_BAR_SCAN, K_1_SCAN, JUMP_FLAG
! TYPE(MPI_File):: FID, FID_GROW, FID_NEGL
INTEGER:: FID, FID_GROW, FID_NEGL, ST

CALL MPI_INIT(IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
IF (MPI_GLB_RANK.EQ.0) THEN
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
    ! open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
ENDIF
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS
CALL READCOM('NOECHO')
CALL READIN(5)

open(unit=LOGID,FILE='discrete_scan.'//TRIM(ITOA3(MPI_GLB_RANK))//'.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)

! open(unit=FID_GROW,FILE='perturb_scan.GROW',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'perturb_scan.GROW', MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, FID_GROW, IS)
IF (IS.NE.0) WRITE(*,*) 'RANK#',MPI_GLB_RANK,'ERROR open perturb_scan.GROW'
CALL MPI_FILE_SET_ATOMICITY(FID_GROW, .TRUE., IERR)
write(*,*) IERR,MPI_GLB_RANK

! open(unit=FID_NEGL,FILE='perturb_scan.NEGLIGIBLE',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'perturb_scan.NEGLIGIBLE', MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, FID_NEGL, IS)
IF (IS.NE.0) WRITE(*,*) 'RANK#',MPI_GLB_RANK,'ERROR open perturb_scan.NEGLIGIBLE'
CALL MPI_FILE_SET_ATOMICITY(FID_NEGL, .TRUE., IERR)
write(*,*) IERR,MPI_GLB_RANK

! ======================================================================
! ============================== DO LOOP ===============================
! {M_bar, AK_bar}: PERTURBATION
!=======================================================================
! create K_BAR_RANGE
K_BAR_START = -10.D0; K_BAR_END = 10.D0
ALLOCATE(K_BAR_RANGE(NINT((K_BAR_END-K_BAR_START)/1.D0) + 1))
CALL LINSPACE(K_BAR_START,K_BAR_END,K_BAR_RANGE)
IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) K_BAR_RANGE

! loop
DO M_BAR = 0,10; DO K_BAR_II = 1,SIZE(K_BAR_RANGE)

    ! obtain K_BAR
    K_BAR = K_BAR_RANGE(K_BAR_II)

    ! skip {0,0} mode
    IF ((M_BAR.EQ.0).AND.(K_BAR.EQ.0.D0)) CYCLE
    ! CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    IF (ALLOCATED(M_mat)) DEALLOCATE(M_mat)
    IF (ALLOCATED(M_eig)) DEALLOCATE(M_eig)
    IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
    IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)    

    ! calculate operator matrix and its corresponding eigenvalues for ea MPI proc
    CALL EIG_MATRIX(M_BAR,K_BAR,M_mat,M_eig,EIG_R_mat,EIG_L_mat,comm_grp=newcomm,serial_switch=.true.)

! ======================================================================
loop_ind: DO EIG_BAR_IND = 1,size(M_mat,1)

    ! set a target Sigma_bar
    EIG_BAR = M_eig(EIG_BAR_IND)
    ! CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    IF (.NOT.((EIGRES(EIG_R_mat(:,EIG_BAR_IND),M_BAR)))) THEN
        CYCLE loop_ind
    ELSE
        ! open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
        ! WRITE(LOGID,*) 'TARGET EIG_BAR:', EIG_BAR, '(# ', ITOA3(EIG_BAR_IND), ') - RANK#',ITOA3(MPI_GLB_RANK)
        ! close(LOGID)
        WRITE(LOGID,*) '(RANK#',ITOA3(MPI_GLB_RANK),') TARGET EIG_BAR: (#',ITOA3(EIG_BAR_IND),')',EIG_BAR
    ENDIF

    ! Save rur_bar, rup_bar, uz_bar, ror_bar, rop_bar, oz_bar for DEGENERATE PERTURBATION THEORY
    CALL EIG2VELVOR(M_BAR,K_BAR,EIG_BAR_IND,EIG_R_mat,EIG_L_mat,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECR_BAR,VECL_BAR,comm_grp=newcomm)
    ALLOCATE(RUR_BARC(NDIMR),RUP_BARC(NDIMR),UZ_BARC(NDIMR),ROR_BARC(NDIMR),ROP_BARC(NDIMR),OZ_BARC(NDIMR))
    RUR_BARC = CONJG(RUR_BAR); RUP_BARC = CONJG(RUP_BAR); UZ_BARC = CONJG( UZ_BAR)
    ROR_BARC = CONJG(ROR_BAR); ROP_BARC = CONJG(ROP_BAR); OZ_BARC = CONJG( OZ_BAR)  

    ! {M_1, AK_1} & {M_2}
    !=======================================================================
    DO M_1 = -M_BAR/2,(10-M_BAR); 
    
        ! obtain M_1 &  M_2 info:
        M_2 = M_1 + M_BAR
        CALL READ_MODE_M(M_1, K_VAL_1_RANGE, M_eig_1_mat, M_res_1_mat)
        CALL READ_MODE_M(M_2, K_VAL_2_RANGE, M_eig_2_mat, M_res_2_mat)

        ! split K_VAL_1_RANGE
        CALL DECOMPOSE(SIZE(K_VAL_1_RANGE), MPI_GLB_PROCS, MPI_GLB_RANK, K_VAL_COUNT, K_VAL_INDEX)

        ! IF (MPI_GLB_RANK.NE.3) GOTO 294

        DO K_VAL_1_II = K_VAL_INDEX+1, K_VAL_INDEX+K_VAL_COUNT

            ! obtain K_VAL_1_0 & K_VAL_2_0
            K_VAL_1_0 = K_VAL_1_RANGE(K_VAL_1_II)
            K_VAL_2_0 = K_VAL_1_0 + K_BAR

            ! skip unnecessary modes:
            ! IF ((M_BAR.EQ.M_1).AND.(K_BAR_SCAN.EQ.K_1_SCAN)) CYCLE
            IF ((M_1.EQ.0).AND.(K_VAL_1_0.EQ.0.D0)) CYCLE                   ! Not interested in m = 0, k = 0 mode
            IF ((M_1+M_BAR.EQ.0).AND.(K_VAL_1_0+K_BAR.EQ.0.D0)) CYCLE       ! Not interested in m = 0, k = 0 mode
            IF ((M_1.EQ.0).AND.(M_BAR.EQ.0.D0).AND.(K_BAR.LE.0.D0)) CYCLE   ! for m_b = m1 = m2 = 0, only need to check positive k_1
            IF (ABS(K_VAL_2_0).GT.MAXVAL(ABS(K_VAL_2_RANGE))) CYCLE         ! k_2 out of range

            ! ! TUNE K
            ! !=======================================================================
            ! 0. Parameters:
            tol = 1.0E-3  ! tolerance

            ! 1. Obtain initial eigvals and eigvecs
            ! K_VAL_2_II = FINDLOC(K_VAL_2_RANGE,K_VAL_2_0,DIM=1)
            K_VAL_2_II = MINLOC(ABS(K_VAL_2_RANGE-K_VAL_2_0),DIM=1)
            IF (ABS(K_VAL_2_RANGE(K_VAL_2_II)-K_VAL_2_0).GE.1E-2) THEN
                WRITE(*,*) 'RANK#',MPI_GLB_RANK,': K_VAL_2 mismatch!'
                CYCLE
            ENDIF
            ALLOCATE(M_eig_1_0(SIZE(M_eig_1_mat,1)), M_eig_2_0(SIZE(M_eig_2_mat,1)))
            M_eig_1_0 = M_eig_1_mat(:,K_VAL_1_II); M_eig_2_0 = M_eig_2_mat(:,K_VAL_2_II)
            ALLOCATE(MASK_1(SIZE(M_eig_1_mat,1)), MASK_2(SIZE(M_eig_2_mat,1)))
            MASK_1 = M_res_1_mat(:,K_VAL_1_II); MASK_2 = M_res_2_mat(:,K_VAL_2_II)

            ! 2. Open log file
            ! open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
            ! WRITE(LOGID,*) '(RANK#',ITOA3(MPI_GLB_RANK),') M_BAR,KBAR: (#',ITOA3(EIG_BAR_IND),')',ITOA3(M_BAR),',',RTOA6(K_BAR),'; M_1,K_1: ', ITOA3(M_1),',',RTOA6(K_VAL_1_0)
            ! close(LOGID)
            WRITE(LOGID,*) '(RANK#',ITOA3(MPI_GLB_RANK),') M_BAR,KBAR: (#',ITOA3(EIG_BAR_IND),')',ITOA3(M_BAR),',',RTOA6(K_BAR),'; M_1,K_1: ', ITOA3(M_1),',',RTOA6(K_VAL_1_0)

            ! 3. Avoid same eigenmode for {mbar,kbar} and {m1,k1}
            IF ((M_BAR.EQ.M_1).AND.(ABS(K_BAR-K_VAL_1_0).LE.0.1)) THEN
                MASK_1 = MASK_1 .AND. (ABS(EIG_BAR-M_eig_1_0).GT.0.01)
            ENDIF
            
            ! 4. Find the minimum distance
            ! tol_diff = MINVAL(ABS(AIMAG(MESHGRID(M_eig_1_0, M_eig_2_0) - EIG_BAR)),MESHGRID(MASK_1,MASK_2))/ABS(AIMAG(EIG_BAR))
            tol_diff = MINVAL(ABS(AIMAG(MESHGRID(M_eig_1_0, M_eig_2_0) - EIG_BAR)),MESHGRID(MASK_1,MASK_2))

            ! 5. Save the result
            IF (tol_diff.LT.tol) THEN

                tol_ind = MINLOC(ABS(AIMAG(MESHGRID(M_eig_1_0, M_eig_2_0) - EIG_BAR)),MESHGRID(MASK_1,MASK_2))
                JJ = tol_ind(1); KK = tol_ind(2)
            
                ! Calculate the FINAL RESULT:
                ! M_bar, K_bar: RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECL_BAR
                CALL EIG_MATRIX(M_1,K_VAL_1_0,M_mat_1_0,M_eig_1_0,EIG_R_mat_1_0,EIG_L_mat_1_0,comm_grp=newcomm,serial_switch=.true.)
                CALL EIG2VELVOR(M_1,K_VAL_1_0,JJ,EIG_R_mat_1_0,EIG_L_mat_1_0,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_1,VECL_1,comm_grp=newcomm)
                CALL EIG_MATRIX(M_2,K_VAL_2_0,M_mat_2_0,M_eig_2_0,EIG_R_mat_2_0,EIG_L_mat_2_0,comm_grp=newcomm,serial_switch=.true.)
                CALL EIG2VELVOR(M_2,K_VAL_2_0,KK,EIG_R_mat_2_0,EIG_L_mat_2_0,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_2,VECL_2,comm_grp=newcomm)
                CALL NONLIN_MK(M_2,K_VAL_2_0,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_B1,comm_grp=newcomm)
                CALL NONLIN_MK(M_1,K_VAL_1_0,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_B2,comm_grp=newcomm)

                ! Sigma(1):
                EIG_PRIME_R = DOT_PRODUCT(VECL_1,VECR_B2)*DOT_PRODUCT(VECL_2,VECR_B1)
                ! Beta/Alpha: ALPHA * {M1,K1} + BETA * {M2,K2}
                ALPHA = 1.D0
                ! BETA = (DOT_PRODUCT(VECL_2,VECR_B1)/DOT_PRODUCT(VECL_1,VECR_B2))**0.5
                BETA = EIG_PRIME_R**0.5/DOT_PRODUCT(VECL_1, VECR_B2)

                ! ! Save the result
                ! DO II = 0,MPI_GLB_PROCS-1; IF (MPI_GLB_RANK.EQ.II) THEN
                ! ! open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
                ! WRITE(LOGID,*) ''
                ! WRITE(LOGID,*) 'Desired tol reached: tol_diff = ', tol_diff
                ! WRITE(LOGID,*) 'DEGEN:'
                ! WRITE(LOGID,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
                ! WRITE(LOGID,201) M_1, K_VAL_1_0, JJ, M_eig_1_0(JJ)
                ! WRITE(LOGID,201) M_2, K_VAL_2_0, kk, M_eig_2_0(KK)
                ! WRITE(LOGID,*) " SIGMA(degen)'^2= ",EIG_PRIME_R
                ! WRITE(LOGID,*) " SIGMA(degen)'^2= ",EIG_PRIME_R/(DOT_PRODUCT(VECL_1,VECR_1)*DOT_PRODUCT(VECL_2,VECR_2))
                ! WRITE(LOGID,*) " SIGMA(degen)'  = ",EIG_PRIME_R**0.5
                ! WRITE(LOGID,*) " BETA(degen)    = ",BETA
                ! 201  FORMAT('Eig: M# = ',I3,'; AK = ',F20.13,'; Sig (#',I3') =',1P8E20.13)
                ! ! close(LOGID)
                ! ! ======================================================================
                ! IF (abs(REAL(EIG_PRIME_R**0.5)).GT.0.01) THEN
                !     open(unit=FID,FILE='perturb_scan.GROW',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
                !     WRITE(FID,*) ''
                !     WRITE(FID,*) 'DEGEN: *GROWTH*'
                ! ELSE
                !     open(unit=FID,FILE='perturb_scan.NEGLIGIBLE',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
                !     WRITE(FID,*) ''
                !     WRITE(FID,*) 'DEGEN: *NEGLIGIBLE*'
                ! ENDIF
                ! WRITE(FID,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
                ! WRITE(FID,201) M_1, K_VAL_1_0, JJ, M_eig_1_0(JJ)
                ! WRITE(FID,201) M_2, K_VAL_2_0, kk, M_eig_2_0(KK)
                ! WRITE(FID,*) "SIGMA(degen)'= RE:",abs(real(EIG_PRIME_R**0.5)),',IM:',imag(EIG_PRIME_R**0.5)
                ! WRITE(FID,*) "BETA/ALPHA   = RE:",abs(real(BETA)),',IM:',imag(BETA)
                ! close(FID)
                ! ENDIF; ENDDO

                ! Save the result
                ! open(unit=LOGID,FILE='discrete_scan.log',STATUS='unknown',ACTION='WRITE',ACCESS='APPEND',IOSTAT=IS)
                WRITE(LOGID,201) M_BAR, K_BAR, EIG_BAR_IND, REAL(EIG_BAR), AIMAG(EIG_BAR), &
                                M_1, K_VAL_1_0, JJ, REAL(M_eig_1_0(JJ)), AIMAG(M_eig_1_0(JJ)), &
                                M_2, K_VAL_2_0, kk, REAL(M_eig_2_0(KK)), AIMAG(M_eig_2_0(KK)), &
                                ABS(REAL(EIG_PRIME_R**0.5)),AIMAG(EIG_PRIME_R**0.5), REAL(BETA), AIMAG(BETA)
                201 FORMAT(3(I3,',',F06.2,',',I3,',',2(G20.12,',')),3(G20.12,','),G20.12)

                WRITE(PRTSTR,201) M_BAR, K_BAR, EIG_BAR_IND, REAL(EIG_BAR), AIMAG(EIG_BAR), &
                                M_1, K_VAL_1_0, JJ, REAL(M_eig_1_0(JJ)), AIMAG(M_eig_1_0(JJ)), &
                                M_2, K_VAL_2_0, kk, REAL(M_eig_2_0(KK)), AIMAG(M_eig_2_0(KK)), &
                                ABS(REAL(EIG_PRIME_R**0.5)),AIMAG(EIG_PRIME_R**0.5), REAL(BETA), AIMAG(BETA)
                
                LONGSTR = TRIM(PRTSTR) // NEW_LINE('A')
                IF (abs(REAL(EIG_PRIME_R**0.5)).GT.0.01) THEN
                    CALL MPI_FILE_WRITE_SHARED(FID_GROW, TRIM(LONGSTR), LEN(TRIM(LONGSTR)) , MPI_CHAR, ST, IERR)
                ELSE
                    CALL MPI_FILE_WRITE_SHARED(FID_NEGL, TRIM(LONGSTR), LEN(TRIM(LONGSTR)) , MPI_CHAR, ST, IERR)
                ENDIF

                ! Deallocate:
                ! ======================================================================
                IF(ALLOCATED(M_mat_1_0)) DEALLOCATE(M_mat_1_0)
                IF(ALLOCATED(M_mat_2_0)) DEALLOCATE(M_mat_2_0)
                ! ======================================================================
                IF(ALLOCATED(EIG_R_mat_1_0)) DEALLOCATE(EIG_R_mat_1_0)
                IF(ALLOCATED(EIG_L_mat_1_0)) DEALLOCATE(EIG_L_mat_1_0)    
                IF(ALLOCATED(EIG_R_mat_2_0)) DEALLOCATE(EIG_R_mat_2_0)
                IF(ALLOCATED(EIG_L_mat_2_0)) DEALLOCATE(EIG_L_mat_2_0)               
                IF(ALLOCATED(RUR_1)) DEALLOCATE(RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1)
                IF(ALLOCATED(RUR_2)) DEALLOCATE(RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2)
                IF(ALLOCATED(VECR_B1)) DEALLOCATE(VECR_B1,VECR_B2,VECL_1,VECR_1,VECL_2,VECR_2)
                ! ======================================================================

                ! ! debug:
                ! close(LOGID)
                ! CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
                ! CALL MPI_FILE_CLOSE(FID_GROW,IERR)
                ! CALL MPI_FILE_CLOSE(FID_NEGL,IERR)
                ! CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)            

            ENDIF

            ! DEALLOCATE 
            !======================================================================= 
            IF(ALLOCATED(MASK_1)) DEALLOCATE(MASK_1,MASK_2)
            IF(ALLOCATED(M_eig_1_0)) DEALLOCATE(M_eig_1_0)
            IF(ALLOCATED(M_eig_2_0)) DEALLOCATE(M_eig_2_0)

        ENDDO ! AK_1
    ENDDO ! M_1

! 294     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!         WRITE(*,*) MPI_GLB_RANK,K_VAL_INDEX+1, K_VAL_INDEX+K_VAL_COUNT
!         CALL MPI_FILE_CLOSE(FID_GROW,IERR)
!         CALL MPI_FILE_CLOSE(FID_NEGL,IERR)
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR) 

    IF(ALLOCATED(RUR_BAR)) DEALLOCATE(RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR)
    IF(ALLOCATED(RUR_BARC)) DEALLOCATE(RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC)
    IF(ALLOCATED(VECL_BAR)) DEALLOCATE(VECL_BAR,VECR_BAR)
ENDDO loop_ind! IND
    IF(ALLOCATED(M_mat)) DEALLOCATE(M_mat)
    IF(ALLOCATED(M_eig)) DEALLOCATE(M_eig)
    IF(ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
    IF(ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)
ENDDO ! AK_BAR
ENDDO ! M_BAR

IF(ALLOCATED(RUR0)) DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)

close(LOGID)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL MPI_FILE_CLOSE(FID_GROW,IERR)
CALL MPI_FILE_CLOSE(FID_NEGL,IERR)

IF (MPI_GLB_RANK.EQ.0) THEN    
    WRITE(*,*) ''
    WRITE(*,*) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
ENDIF
CALL MPI_FINALIZE(IERR)

! ======================================================================
CONTAINS
! ======================================================================
SUBROUTINE READ_MODE_M(M_READ, K_RANGE, EIG_MAT, RES_MAT)
! ======================================================================

! ======================================================================
IMPLICIT NONE
INTEGER:: M_READ, M_FILE
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: EIG_MAT
LOGICAL,DIMENSION(:,:),ALLOCATABLE:: RES_MAT

INTEGER:: K_COUNT, EIG_SIZE
REAL(P8):: K_START,K_END
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_RANGE

CHARACTER(LEN=72):: FN
INTEGER:: STATUS

! open file
FN = 'EIG_'//ITOA4(M_READ)//'.dat'
OPEN(UNIT=7,FILE='./data/'//TRIM(FN),STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',IOSTAT=STATUS)
IF (STATUS.NE.0) THEN
    WRITE(*,*) 'READ_MODE_M: FAILED TO OPEN ',TRIM(FN),' - RANK#',MPI_GLB_RANK
    STOP
ENDIF

READ(7) M_FILE,K_START,K_END
READ(7) EIG_SIZE,K_COUNT

IF (M_FILE.NE.M_READ) THEN
    WRITE(*,*) 'READ_MODE_M: M_READ = ', M_READ,' .NE. M_FILE = ', M_FILE,' - RANK#',MPI_GLB_RANK
    STOP
ENDIF

! obtain K_RANGE
IF (ALLOCATED(K_RANGE)) DEALLOCATE(K_RANGE)
ALLOCATE(K_RANGE(K_COUNT))
CALL LINSPACE(K_START,K_END,K_RANGE)

! read EIG_MAT and 
IF (ALLOCATED(EIG_MAT)) DEALLOCATE(EIG_MAT)
IF (ALLOCATED(RES_MAT)) DEALLOCATE(RES_MAT)
ALLOCATE(EIG_MAT(EIG_SIZE,K_COUNT))
ALLOCATE(RES_MAT(EIG_SIZE,K_COUNT))

READ(7) EIG_MAT(:,:)
READ(7) RES_MAT(:,:)
CLOSE(7)

END SUBROUTINE
! ======================================================================

END PROGRAM EVP_SCAN_0920
!=======================================================================
