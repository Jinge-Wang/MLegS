program sample
    ! PROGRAM TO ASSEMBLE LOCAL ARRAYS INTO GLOBAL ARRAYS
        USE OMP_LIB
        USE MPI
        USE MOD_MISC, ONLY : P4,P8,PI,IU, ITOA3, MCAT                      ! LEVEL 0
        USE MOD_EIG                                                        ! LEVEL 1
        USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
        USE MOD_SCALAR3                                                    ! LEVEL 2
        USE MOD_FFT
    implicit none
    
    ! PROC INFO:
    integer:: PROC_NUM, IERR
    integer,DIMENSION(:),ALLOCATABLE:: SUBCOMMS
    
    
    ! ============================ LOCAL ARRAY =============================
    integer,DIMENSION(:,:,:),ALLOCATABLE:: LOCAL_ARRAY2, LOCAL_ARRAY
    ! INDEX inside LOCAL ARRAY:
    INTEGER:: I, J, K
    ! SUBARRAY inside LOCAL ARRAAY:
    INTEGER:: INT_SIZE
    INTEGER(KIND=MPI_ADDRESS_KIND):: EXTEND_SIZE
    INTEGER:: SUBARRAY_TYPE, SUBARRAY_TYPE_resized
    ! ========================== SUBGLOBAL ARRAY ===========================
    integer,DIMENSION(:,:,:),ALLOCATABLE:: SUBGLOBAL_ARRAY
    ! INDEX inside SUBGLOBAL ARRAY:
    INTEGER:: DISPLACEMENT_loc, RECVCOUNT_loc
    integer,DIMENSION(:),ALLOCATABLE:: DISPLACEMENT, RECVCOUNT
    ! =========================== GLOBAL ARRAY =============================
    INTEGER:: N1, N2, N3
    integer,DIMENSION(:,:,:),ALLOCATABLE:: GLOBAL_ARRAY, GLOBAL_ARRAY_COPY
    
    
    
    CALL MPI_INIT(IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_NUM, IERR)
    
    ! CREATE CARTESIAN GRID OF PROCS
    CALL SUBCOMM_CART(MPI_COMM_WORLD, 2, SUBCOMMS)
    
    ! TEST ASSEMBLE SUBARRAY
    N1 = 12
    N2 = 4
    N3 = 17
    
    ! CREATE LOCAL ARRAY
    ! ALLOCATE(LOCAL_ARRAY(1:N1,1:local_size(N2,SUBCOMMS(1)),1:local_size(N3,SUBCOMMS(2))))
    ! ALLOCATE(LOCAL_ARRAY2(1:N1,1:local_size(N3,SUBCOMMS(2)),1:local_size(N2,SUBCOMMS(1))))
    ALLOCATE(LOCAL_ARRAY(1:local_size(N1,SUBCOMMS(1)),1:N2,1:local_size(N3,SUBCOMMS(2))))
    ALLOCATE(LOCAL_ARRAY2(1:N2,1:local_size(N3,SUBCOMMS(2)),1:local_size(N1,SUBCOMMS(1))))
    DO K = 1,SIZE(LOCAL_ARRAY,3)
        DO J = 1,SIZE(LOCAL_ARRAY,2)
            DO I = 1,SIZE(LOCAL_ARRAY,1)
                LOCAL_ARRAY(I,J,K) = (I+local_index(N1,SUBCOMMS(1)))*10000 & 
                                    + J*100 &
                                    + (K+local_index(N3,SUBCOMMS(2)))
            ENDDO
        ENDDO
    ENDDO
    
    
    ! ===================== GATHER INSIDE SUBCOMMS(2) ======================
    ! LOCAL_ARRAY2 = RESHAPE(LOCAL_ARRAY,SHAPE(LOCAL_ARRAY2),ORDER = [2,3,1])
    DO K = 1,SIZE(LOCAL_ARRAY2,3)
        DO J = 1,SIZE(LOCAL_ARRAY2,2)
            DO I = 1,SIZE(LOCAL_ARRAY2,1)
                LOCAL_ARRAY2(I,J,K) = LOCAL_ARRAY(K,I,J)
            ENDDO
        ENDDO
    ENDDO
    
    ! VECTOR TYPE
    CALL MPI_TYPE_VECTOR(local_size(N3,SUBCOMMS(2)), N2, N2*N1, &
                        MPI_INTEGER, SUBARRAY_TYPE, IERR)
    CALL MPI_TYPE_SIZE(MPI_INTEGER,INT_SIZE,IERR)
    EXTEND_SIZE = INT_SIZE*N2
    CALL MPI_TYPE_CREATE_RESIZED(SUBARRAY_TYPE, 0, EXTEND_SIZE, SUBARRAY_TYPE_resized, IERR)
    CALL MPI_TYPE_COMMIT(SUBARRAY_TYPE_resized,IERR)
    
    ! PREP GATHERV
    ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMMS(1))-1))
    ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMMS(1))-1))
    RECVCOUNT_loc = local_size(N1,SUBCOMMS(1))
    CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMMS(1),IERR)
    DISPLACEMENT_loc = local_index(N1,SUBCOMMS(1))
    CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMMS(1),IERR)
    
    ! GATHERV
    ALLOCATE(SUBGLOBAL_ARRAY(N2,N1,local_size(N3,SUBCOMMS(2))))
    CALL MPI_GATHERV(LOCAL_ARRAY2, SIZE(LOCAL_ARRAY,1)*SIZE(LOCAL_ARRAY,2)*SIZE(LOCAL_ARRAY,3), &
        MPI_INTEGER, SUBGLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, SUBARRAY_TYPE_resized,0, &
        SUBCOMMS(1), IERR)
    
    ! ======================== GATHER GLOBAL ARRAY =========================
    ALLOCATE(GLOBAL_ARRAY(N2,N1,N3))
    IF (local_proc(SUBCOMMS(1)).EQ.0) THEN
        CALL MPI_GATHER(SUBGLOBAL_ARRAY,N1*N2*SIZE(LOCAL_ARRAY,3),MPI_INTEGER, &
                        GLOBAL_ARRAY,N1*N2*SIZE(LOCAL_ARRAY,3),MPI_INTEGER,0, &
                        SUBCOMMS(2), IERR)
    ENDIF
    ALLOCATE(GLOBAL_ARRAY_COPY(N1,N2,N3))
    GLOBAL_ARRAY_COPY = RESHAPE(GLOBAL_ARRAY,SHAPE(GLOBAL_ARRAY_COPY),ORDER = [2,1,3])
    
    ! ============================== OUTPUT ================================
    IF (PROC_NUM.EQ.0) THEN
        ! WRITE(*,*) count_proc(SUBCOMMS(1))
        ! WRITE(*,*) count_proc(SUBCOMMS(2))
        ! WRITE(*,*) GLOBAL_ARRAY(1,:,:)
        CALL MCAT(GLOBAL_ARRAY_COPY(3,:,:))
        WRITE(*,*) ''
        ! WRITE(*,*) LOCAL_ARRAY(1,:,:)
        CALL MCAT(SUBGLOBAL_ARRAY(:,1,:))
        WRITE(*,*) ''
        CALL MCAT(LOCAL_ARRAY(2,:,:))
        WRITE(*,*) ''
        CALL MCAT(LOCAL_ARRAY2(:,:,2))
        WRITE(*,*) ''
    ENDIF
        ! WRITE(*,*) DISPLACEMENT
        ! WRITE(*,*) 'RECV:'
        ! WRITE(*,*) RECVCOUNT
    
        ! WRITE(*,*) ITOA3(PROC_NUM)//':'//ITOA3(local_size(N2,SUBCOMMS(1)))
        ! WRITE(*,*) ITOA3(PROC_NUM)//':'//ITOA3(local_size(N3,SUBCOMMS(2)))
    
    
    ! IF (PROC_NUM.EQ.2) THEN
    !     WRITE(*,*) LOCAL_ARRAY(1,:,:)
    ! ENDIF
    CALL MPI_FINALIZE(IERR)
    
    CONTAINS
    
    
    end program sample