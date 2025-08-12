program TEST_EXCHANGE
! ======================================================================
! test 3d exchange subroutines which is essential for the pencil decomp
! via MPI.
! ======================================================================
   USE OMP_LIB
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
   USE MOD_DIAGNOSTICS
   USE MOD_INIT
   USE MOD_EVP, only: SAVE_PERTURB
! -------------------------
implicit none
! -------------------------
integer:: ii, jj, kk, ni, nj, nk, si, sj, sk, fileid, cart_coor(2), total_procs
complex(P8),dimension(:,:,:),allocatable:: A, B, C, D
complex(P8):: DIFF

CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
   WRITE(*,*) 'The threading support is lesser than that demanded.'
   CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

CALL READCOM('NOECHO')
CALL READIN(5)
CALL SETUP_GRID()

CALL MPI_COMM_SIZE(MPI_COMM_IVP, total_procs, IERR)
CALL MPI_CART_COORDS(MPI_COMM_IVP,MPI_RANK,2,cart_coor,IERR);

fileid = 777

open(UNIT=fileid,FILE='./test/'//ITOA3(total_procs)//'_check_'//ITOA3(cart_coor(1))//'_'//ITOA3(cart_coor(2))//'.dat',&
   &STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')

CALL DISPLAY_COMM_INFO(fileid)

! ======================================================================
ALLOCATE(A(SIZE_PFF1(1),SIZE_PFF1(2),SIZE_PFF1(3)),D(SIZE_PFF1(1),SIZE_PFF1(2),SIZE_PFF1(3)))
ALLOCATE(B(SIZE_PFF0(1),SIZE_PFF0(2),SIZE_PFF0(3)),C(SIZE_PFF0(1),SIZE_PFF0(2),SIZE_PFF0(3)))

! PFF1: A & D
! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1 -> TYPE_PFF1
SI = 0
SJ = local_index(NTCHOPDIM,SUBCOMM_2)
SK = local_index(NXCHOPDIM,SUBCOMM_1)
DO II = 1,SIZE(A,1)
    DO JJ = 1,SIZE(A,2)
        DO KK = 1,SIZE(A,3)
            NI = II + SI
            NJ = JJ + SJ
            NK = KK + SK
            A(II,JJ,KK) = 10**6 * NI + 10**3 * NJ + NK! + MPI_RANK * IU
            A(II,JJ,KK) = A(II,JJ,KK) + IU * (10**6 * NI + 10**3 * NJ + NK)
        ENDDO
    ENDDO
ENDDO

CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_1, A, TYPE_PFF1, B, TYPE_PFF0)
call MPI_BARRIER(MPI_COMM_IVP,IERR)

! PFF0: B & C
! PFF: NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM -> TYPE_PFF0
WRITE(fileid,*) "TYPE_PFF1 EXCHANGED TO TYPE_PFF0"
SI = local_index(NDIMR,SUBCOMM_1)
SJ = local_index(NTCHOPDIM,SUBCOMM_2)
SK = 0
DO II = 1,SIZE(C,1)
    DO JJ = 1,SIZE(C,2)
        DO KK = 1,SIZE(C,3)
            NI = II + SI
            NJ = JJ + SJ
            NK = KK + SK
            C(II,JJ,KK) = 10**6 * NI + 10**3 * NJ + NK! + MPI_RANK * IU
            C(II,JJ,KK) = C(II,JJ,KK) + IU * (10**6 * NI + 10**3 * NJ + NK)

            DIFF = C(II,JJ,KK) - B(II,JJ,KK)
            IF (ABS(DIFF).NE.0.D0) THEN
                WRITE(fileid,*) "difference at i = ",ITOA3(NI), &
                                            ", j = ",ITOA3(NJ), &
                                            ", k = ",ITOA3(NK)
            ENDIF

        ENDDO
    ENDDO
ENDDO

CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_1, C, TYPE_PFF0, D, TYPE_PFF1)
call MPI_BARRIER(MPI_COMM_IVP,IERR)

WRITE(fileid,*) "TYPE_PFF0 EXCHANGED TO TYPE_PFF1"
DO II = 1,SIZE(D,1)
    DO JJ = 1,SIZE(D,2)
        DO KK = 1,SIZE(D,3)
            NI = II + SI
            NJ = JJ + SJ
            NK = KK + SK

            DIFF = D(II,JJ,KK) - A(II,JJ,KK)
            IF (ABS(DIFF).NE.0.D0) THEN
                WRITE(fileid,*) "difference at i = ",ITOA3(NI), &
                                            ", j = ",ITOA3(NJ), &
                                            ", k = ",ITOA3(NK)
            ENDIF

        ENDDO
    ENDDO
ENDDO

DEALLOCATE(A,B,C,D)

! ======================================================================

ALLOCATE(A(SIZE_PFP0(1),SIZE_PFP0(2),SIZE_PFP0(3)),D(SIZE_PFP0(1),SIZE_PFP0(2),SIZE_PFP0(3)))
ALLOCATE(B(SIZE_PFP1(1),SIZE_PFP1(2),SIZE_PFP1(3)),C(SIZE_PFP1(1),SIZE_PFP1(2),SIZE_PFP1(3)))

! PFF1: A & D
! PFP: NDIMR/N1, NTCHOPDIM, NDIMX/N2 -> TYPE_PFP0
SI = local_index(NDIMR,SUBCOMM_1)
SJ = 0
SK = local_index(NDIMX,SUBCOMM_2)
DO II = 1,SIZE(A,1)
    DO JJ = 1,SIZE(A,2)
        DO KK = 1,SIZE(A,3)
            NI = II + SI
            NJ = JJ + SJ
            NK = KK + SK
            A(II,JJ,KK) = 10**6 * NI + 10**3 * NJ + NK! + MPI_RANK * IU
            A(II,JJ,KK) = A(II,JJ,KK) + IU * (10**6 * NI + 10**3 * NJ + NK)
        ENDDO
    ENDDO
ENDDO

CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_2, A, TYPE_PFP0, B, TYPE_PFP1)
call MPI_BARRIER(MPI_COMM_IVP,IERR)

! PFF0: B & C
! PFP: NDIMR/N1, NTCHOPDIM/N2, NDIMX -> TYPE_PFP1
WRITE(fileid,*) "TYPE_PFP0 EXCHANGED TO TYPE_PFP1"
SI = local_index(NDIMR,SUBCOMM_1)
SJ = local_index(NTCHOPDIM,SUBCOMM_2)
SK = 0
DO II = 1,SIZE(C,1)
    DO JJ = 1,SIZE(C,2)
        DO KK = 1,SIZE(C,3)
            NI = II + SI
            NJ = JJ + SJ
            NK = KK + SK
            C(II,JJ,KK) = 10**6 * NI + 10**3 * NJ + NK! + MPI_RANK * IU
            C(II,JJ,KK) = C(II,JJ,KK) + IU * (10**6 * NI + 10**3 * NJ + NK)

            DIFF = B(II,JJ,KK) - C(II,JJ,KK)
            IF (ABS(DIFF).NE.0.D0) THEN
                WRITE(fileid,*) "difference at i = ",ITOA3(NI), &
                                            ", j = ",ITOA3(NJ), &
                                            ", k = ",ITOA3(NK)
            ENDIF

        ENDDO
    ENDDO
ENDDO

CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_2, C, TYPE_PFP1, D, TYPE_PFP0)
call MPI_BARRIER(MPI_COMM_IVP,IERR)

WRITE(fileid,*) "TYPE_PFP1 EXCHANGED TO TYPE_PFP1"
DO II = 1,SIZE(D,1)
    DO JJ = 1,SIZE(D,2)
        DO KK = 1,SIZE(D,3)
            NI = II + SI
            NJ = JJ + SJ
            NK = KK + SK

            DIFF = D(II,JJ,KK) - A(II,JJ,KK)
            IF (ABS(DIFF).NE.0.D0) THEN
                WRITE(fileid,*) "difference at i = ",ITOA3(NI), &
                                            ", j = ",ITOA3(NJ), &
                                            ", k = ",ITOA3(NK)
            ENDIF

        ENDDO
    ENDDO
ENDDO

DEALLOCATE(A,B,C,D)

close(fileid)

!> final printout
IF (MPI_RANK.eq.0) THEN
    WRITE(*,*) 'PROGRAM STARTED'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

! ======================================================================
CONTAINS
! ======================================================================
SUBROUTINE DISPLAY_COMM_INFO(FID)
! ======================================================================
! WRITE COMMUNICATOR INFO TO FILE
! ======================================================================
INTEGER:: FID
INTEGER:: nproc_x, nproc_y, iproc_x, iproc_y

write(FID,*) "COMM" 
write(FID,*) "rank = ", ITOA3(MPI_RANK)

CALL MPI_COMM_SIZE(SUBCOMM_1, nproc_x, IERR);
CALL MPI_COMM_RANK(SUBCOMM_1, iproc_x, IERR);
write(FID,*) "rank_x = ",ITOA3(iproc_x),' out of ',ITOA3(nproc_x)

CALL MPI_COMM_SIZE(SUBCOMM_2, nproc_y, IERR);
CALL MPI_COMM_RANK(SUBCOMM_2, iproc_y, IERR);
write(FID,*) "rank_y = ",ITOA3(iproc_y),' out of ',ITOA3(nproc_y)

! PFP: NDIMR/N1, NTCHOPDIM, NDIMX/N2 -> TYPE_PFP0
WRITE(FID,*) "PFP0: ",ITOA3(NDIMR),' X ',ITOA3(NTCHOPDIM),' X ',ITOA3(NDIMX)
WRITE(FID,*) "------",ITOA3(local_index(NDIMR,SUBCOMM_1)), &
                ' X ',ITOA3(0), &
                ' X ',ITOA3(local_index(NDIMX,SUBCOMM_2))
WRITE(FID,*) "------",ITOA3(SIZE_PFP0(1)),' X ',ITOA3(SIZE_PFP0(2)),' X ',ITOA3(SIZE_PFP0(3))

! PFP: NDIMR/N1, NTCHOPDIM/N2, NDIMX -> TYPE_PFP1
WRITE(FID,*) "PFP1: ",ITOA3(NDIMR),' X ',ITOA3(NTCHOPDIM),' X ',ITOA3(NDIMX)
WRITE(FID,*) "------",ITOA3(local_index(NDIMR,SUBCOMM_1)), &
                ' X ',ITOA3(local_index(NTCHOPDIM,SUBCOMM_2)), &
                ' X ',ITOA3(0)
WRITE(FID,*) "------",ITOA3(SIZE_PFP1(1)),' X ',ITOA3(SIZE_PFP1(2)),' X ',ITOA3(SIZE_PFP1(3))

! PFF: NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM -> TYPE_PFF0
WRITE(FID,*) "PFF0: ",ITOA3(NDIMR),' X ',ITOA3(NTCHOPDIM),' X ',ITOA3(NXCHOPDIM)
WRITE(FID,*) "------",ITOA3(local_index(NDIMR,SUBCOMM_1)), &
                ' X ',ITOA3(local_index(NTCHOPDIM,SUBCOMM_2)), &
                ' X ',ITOA3(0)
WRITE(FID,*) "------",ITOA3(SIZE_PFF0(1)),' X ',ITOA3(SIZE_PFF0(2)),' X ',ITOA3(SIZE_PFF0(3))

! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1 -> TYPE_PFF1
WRITE(FID,*) "PFF1: ",ITOA3(NDIMR),' X ',ITOA3(NTCHOPDIM),' X ',ITOA3(NXCHOPDIM)
WRITE(FID,*) "------",ITOA3(0), &
                ' X ',ITOA3(local_index(NTCHOPDIM,SUBCOMM_2)), &
                ' X ',ITOA3(local_index(NXCHOPDIM,SUBCOMM_1))
WRITE(FID,*) "------",ITOA3(SIZE_PFF1(1)),' X ',ITOA3(SIZE_PFF1(2)),' X ',ITOA3(SIZE_PFF1(3))

END SUBROUTINE

END PROGRAM TEST_EXCHANGE