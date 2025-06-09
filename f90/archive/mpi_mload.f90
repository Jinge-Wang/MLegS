program sample2
    ! PROGRAM TO ASSEMBLE LOCAL ARRAYS INTO GLOBAL ARRAYS
        USE OMP_LIB
        USE MPI
        USE MOD_MISC, ONLY : P4,P8,PI,IU, ITOA3, MCAT                      ! LEVEL 0
        USE MOD_EIG                                                        ! LEVEL 1
        USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
        USE MOD_SCALAR3                                                    ! LEVEL 2
        USE MOD_FFT
    implicit none
    
    ! Local Scalar:
    TYPE(SCALAR):: local_scalar
    CHARACTER(LEN=72):: FN
    INTEGER:: STATUS, MPIFILE, INDEX, DISP
    INTEGER:: SIZE1, SIZE2, SIZE3, ARRAYSIZE, ELEMENT_SIZE, FILEBLK, FILESTATUS, LOCAL_ARRAYSIZE
    REAL(P8):: IZLEN, IELL
    INTEGER:: IMINC, IMKLINK
    INTEGER:: IN1, IN2, IN3, IN1_2, IN2_2, IN3_3
    INTEGER(KIND=MPI_ADDRESS_KIND):: DISPLACEMENT

    
    ! PUBLIC:
    integer,DIMENSION(:),ALLOCATABLE:: SUBCOMMS
    
    ! TEMPORARY:
    INTEGER:: N1, N2, N3, PROC_NUM, ERR, I, J, K
    
    CALL MPI_INIT(ERR)
    
    ! CREATE CARTESIAN GRID OF PROCS
    CALL SUBCOMM_CART(MPI_COMM_WORLD, 2, SUBCOMMS)

    ZLEN = 1.D0
    ELL = 1.D0
    MINC = 1
    MKLINK = 0
    
    SUBCOMM_1 = SUBCOMMS(1)
    SUBCOMM_2 = SUBCOMMS(2)
    
    ! ! CREATE LOCAL SCALAR
    ! local_scalar%SPACE = PPP_SPACE
    ! local_scalar%LN = 0.D0
    ! local_scalar%INR = local_index(N1,SUBCOMM_1)
    ! local_scalar%INTH = 0
    ! local_scalar%INX = local_index(N3,SUBCOMM_2)
    ! ALLOCATE(local_scalar%E(local_size(N1,SUBCOMM_1),N2,local_size(N3,SUBCOMM_2)))
    
    ! MLOAD
    FN = 'local_scalar.dat'
    CALL MLOAD(FN,local_scalar)
    ! ! 1. save basic info:
    ! CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_NUM, ERR)
    ! IF (PROC_NUM.EQ.0) THEN
    
    !     OPEN(UNIT=7,FILE=TRIM(FN)//'.info',STATUS='OLD',&
    !     FORM='UNFORMATTED',IOSTAT=STATUS)
    
    !     IF(STATUS.NE.0) THEN
    !         WRITE(*,*) 'MSAVEX_IN_SCALAR3: ROOT PROC FAILED TO READ ',TRIM(FN),'.info'
    !         STOP
    !     ENDIF
    
    !     READ(7) IN1,IN2,IN3
    !     READ(7) IN1_2,IN2_2,IN3_3
    !     READ(7) local_scalar%SPACE,local_scalar%LN
    !     READ(7) IZLEN,IELL
    !     READ(7) IMINC,IMKLINK   
    !     CLOSE(7)

    !     ! Check consistency
    !     IF(IZLEN.NE.ZLEN) THEN
    !         WRITE(*,*) 'MLOADX: ZLEN INCONSISTENT.'
    !         WRITE(*,*) 'PROG =',ZLEN
    !         WRITE(*,*) 'FILE =',IZLEN
    !     ENDIF

    !     IF(IELL.NE.ELL) THEN
    !         WRITE(*,*) 'MLOADX: ELL INCONSISTENT.'
    !         WRITE(*,*) 'PROG =',ELL
    !         WRITE(*,*) 'FILE =',IELL
    !     ENDIF

    !     IF(MINC.NE.IMINC) THEN
    !         WRITE(*,*) 'MLOADX: MINC INCONSISTENT.'
    !         WRITE(*,*) 'PROG =',MINC
    !         WRITE(*,*) 'FILE =',IMINC
    !     ENDIF

    !     IF(MKLINK.NE.IMKLINK) THEN
    !         WRITE(*,*) 'MLOADX: MKLINK INCONSISTENT.'
    !         WRITE(*,*) 'PROG =',MKLINK
    !         WRITE(*,*) 'FILE =',IMKLINK
    !     ENDIF

    !     ! DEBUG:
    !     ! WRITE(*,*) IZLEN
    !     ! WRITE(*,*) IELL
    !     ! WRITE(*,*) IN1,IN2,IN3
    !     ! WRITE(*,*) local_scalar%SPACE,local_scalar%LN
    ! ENDIF
    ! CALL MPI_BCAST(IN1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ERR)
    ! CALL MPI_BCAST(IN2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ERR)
    ! CALL MPI_BCAST(IN3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ERR)
    ! CALL MPI_BCAST(local_scalar%SPACE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ERR)

    ! ! 2. determine element size
    ! SIZE1 = CEILING(REAL(IN1)/count_proc(SUBCOMM_1))
    ! SIZE2 = IN2
    ! SIZE3 = CEILING(REAL(IN3)/count_proc(SUBCOMM_2))
    ! ARRAYSIZE = SIZE1*SIZE2*SIZE3
    
    ! ! 3. determine Processor index
    ! INDEX = local_proc(SUBCOMM_2) + local_proc(SUBCOMM_1) * count_proc(SUBCOMM_2)
    
    ! ! 4. determine element size of MPI_DOUBLE_COMPLEX
    ! CALL MPI_TYPE_SIZE(MPI_Double_Complex,ELEMENT_SIZE,ERR)
    
    ! ! 5. create FILEBLK (derived datatype)
    ! call MPI_TYPE_CONTIGUOUS(ARRAYSIZE, MPI_DOUBLE_COMPLEX, FILEBLK, ERR)
    ! call MPI_TYPE_COMMIT(FILEBLK, ERR)
    
    ! ! 6. save A%E in each proc
    ! ALLOCATE(local_scalar%E(local_size(IN1,SUBCOMM_1),IN2,local_size(IN3,SUBCOMM_2)))
    ! LOCAL_ARRAYSIZE = SIZE(local_scalar%E,1)*SIZE(local_scalar%E,2)*SIZE(local_scalar%E,3)
    ! CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FN, MPI_MODE_RDONLY, MPI_INFO_NULL, MPIFILE, ERR)
    ! DISPLACEMENT = INDEX*ARRAYSIZE*ELEMENT_SIZE
    ! CALL MPI_FILE_SET_VIEW(MPIFILE, DISPLACEMENT, MPI_DOUBLE_COMPLEX, FILEBLK, &
    !                     'NATIVE', MPI_INFO_NULL, ERR)
    ! CALL MPI_FILE_READ_ALL(MPIFILE, local_scalar%E, LOCAL_ARRAYSIZE, MPI_DOUBLE_COMPLEX, FILESTATUS, ERR)
    ! CALL MPI_FILE_CLOSE(MPIFILE, ERR)
    
    ! check result:
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_NUM, ERR)
    IF (PROC_NUM.EQ.2) THEN
        CALL MCAT(INT(REAL(local_scalar%E(:,2,:))))
    ENDIF
    ! IF (PROC_NUM.EQ.1) THEN
    !     CALL MCAT(REAL(local_scalar%E(2,:,:)))
    ! ENDIF

    CALL MPI_FINALIZE(ERR)
    
    
    ! CONTAINS
    
    ! SUBROUTINE ASSEMBLE(LOCAL_ARRAY, GLOBAL_ARRAY, axis, MPI_Datatype)
    ! ! ======================================================================
    ! ! [USAGE]: 
    ! ! ASSEMBLE LOCAL ARRAYS INTO GLOBAL ARRAY IN PROC#0
    ! ! [PARAMETERS]:
    ! ! axis >> axis along which the domain is NOT chopped    
    ! ! MPI_Datatype >> element type of local array (e.g. MPI_INTEGER)
    ! ! [UPDATES]:
    ! ! CODED BY JINGE WANG @ SEP 19 2021
    ! ! [NOTE1]:
    ! ! SUBCOMMS_L,axis,SUBCOMMS_R
    ! ! or axis,SUBCOMMS_L,SUBCOMMS_R
    ! ! DOES NOT work with SUBCOMMS_L,SUBCOMMS_R,axis (not needed in IVP)
    ! ! [NOTE2]:
    ! ! USE public SUBCOMM_1 and SUBCOMM_2. Need to change their definition to
    ! ! match the actual usage.
    ! ! ======================================================================
    ! integer,DIMENSION(:,:,:),INTENT(IN):: LOCAL_ARRAY
    ! integer,DIMENSION(:,:,:),INTENT(INOUT):: GLOBAL_ARRAY
    ! integer,INTENT(IN):: axis, MPI_Datatype
    
    ! INTEGER:: II, JJ, KK, N1_glb, N2_glb, N3_glb, N1_loc, N2_loc, N3_loc, ELEMENT_SIZE, IERR
    ! INTEGER(KIND=MPI_ADDRESS_KIND):: EXTEND_SIZE
    ! INTEGER:: SUBCOMM_L, SUBCOMM_R
    ! INTEGER:: SUBARRAY_TYPE, SUBARRAY_TYPE_resized, DISPLACEMENT_loc, RECVCOUNT_loc
    ! integer,DIMENSION(:),ALLOCATABLE:: DISPLACEMENT, RECVCOUNT
    ! integer,DIMENSION(:,:,:),ALLOCATABLE:: LOCAL_ARRAY2, SUBGLOBAL_ARRAY, GLOBAL_ARRAY_COPY
    
    ! N1_glb = SIZE(GLOBAL_ARRAY,1)
    ! N2_glb = SIZE(GLOBAL_ARRAY,2)
    ! N3_glb = SIZE(GLOBAL_ARRAY,3)
    
    ! N1_loc = SIZE(LOCAL_ARRAY,1)
    ! N2_loc = SIZE(LOCAL_ARRAY,2)
    ! N3_loc = SIZE(LOCAL_ARRAY,3)
    
    ! CALL MPI_TYPE_SIZE(MPI_Datatype,ELEMENT_SIZE,IERR)
    
    ! ! SUBCOMMS_L,axis,SUBCOMMS_R
    ! ! or axis,SUBCOMMS_L,SUBCOMMS_R
    ! IF (axis.EQ.1) THEN
    !     SUBCOMM_L = SUBCOMM_1
    !     SUBCOMM_R = SUBCOMM_2
    
    !     ALLOCATE(SUBGLOBAL_ARRAY(N1_glb,N2_glb,N3_loc))
    !     ALLOCATE(LOCAL_ARRAY2(N1_glb,N3_loc,N2_loc))
    
    !     LOCAL_ARRAY2 = RESHAPE(LOCAL_ARRAY,SHAPE(LOCAL_ARRAY2),ORDER = [1,3,2])
    
    !     CALL MPI_TYPE_VECTOR(N3_loc, N1_glb, N2_glb*N1_glb, &
    !                          MPI_Datatype, SUBARRAY_TYPE, IERR)
    !     EXTEND_SIZE = ELEMENT_SIZE*N1_glb
    !     RECVCOUNT_loc = N2_loc
    !     DISPLACEMENT_loc = local_index(N2_glb,SUBCOMM_L)
    
    ! ELSEIF (axis.EQ.2) THEN
    !     SUBCOMM_L = SUBCOMM_1
    !     SUBCOMM_R = SUBCOMM_2
    
    !     ALLOCATE(SUBGLOBAL_ARRAY(N2_glb,N1_glb,N3_loc))
    !     ALLOCATE(LOCAL_ARRAY2(N2_glb,N3_loc,N1_loc))
    
    !     DO KK = 1,N1_loc
    !         DO JJ = 1,N3_loc
    !             DO II = 1,N2_glb
    !                 LOCAL_ARRAY2(II,JJ,KK) = LOCAL_ARRAY(KK,II,JJ)
    !             ENDDO
    !         ENDDO
    !     ENDDO
    
    !     CALL MPI_TYPE_VECTOR(N3_loc, N2_glb, N2_glb*N1_glb, &
    !                          MPI_Datatype, SUBARRAY_TYPE, IERR)
    !     EXTEND_SIZE = ELEMENT_SIZE*N2_glb
    !     RECVCOUNT_loc = N1_loc
    !     DISPLACEMENT_loc = local_index(N1_glb,SUBCOMM_L)
    
    ! ELSE
    !     WRITE(*,*) 'GATHER: undefined for the current axis'
    !     RETURN
    ! ENDIF
    
    ! ! CREATE ELEMENTAL SUBARRAY OF LOCAL_ARRAY
    ! CALL MPI_TYPE_CREATE_RESIZED(SUBARRAY_TYPE, 0, EXTEND_SIZE, SUBARRAY_TYPE_resized, IERR)
    ! CALL MPI_TYPE_COMMIT(SUBARRAY_TYPE_resized,IERR)
    
    ! ! IN EA SUBCOMM_L, GATHER LOCAL_ARRAY TO FORM SUBGLOBAL_ARRAY in SUBCOMM_L's #0 PROC
    ! ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMM_L)-1))
    ! ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMM_L)-1))
    ! CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMM_L,IERR)
    ! CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMM_L,IERR)
    ! CALL MPI_GATHERV(LOCAL_ARRAY2, N1_loc*N2_loc*N3_loc, MPI_Datatype, &
    !                  SUBGLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, SUBARRAY_TYPE_resized, &
    !                  0, SUBCOMM_L, IERR)
    
    ! ! IN THE SUBCOMM_R THAT CONTAINS SUBCOMM_Ls' #0 PROCs: 
    ! ! GATHER SUBGLOBAL_ARRAY TO FORM GLOBAL_ARRAY in SUBCOMM_R's #0 proc (=> GLOBAL #0 PROC)
    ! IF (local_proc(SUBCOMM_L).EQ.0) THEN
    
    !     DEALLOCATE(RECVCOUNT,DISPLACEMENT)
    !     ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMM_R)-1))
    !     ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMM_R)-1))
    !     RECVCOUNT_loc = N1_glb*N2_glb*N3_loc
    !     DISPLACEMENT_loc = N1_glb*N2_glb*local_index(N3,SUBCOMM_R)
    !     CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMM_R,IERR)
    !     CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMM_R,IERR)
    
    !     IF (axis.EQ.1) THEN
       
    !         ! CALL MPI_GATHER(SUBGLOBAL_ARRAY,N1_glb*N2_glb*N3_loc,MPI_Datatype, &
    !         !                 GLOBAL_ARRAY,N1_glb*N2_glb*N3_loc,MPI_Datatype, &
    !         !                 0, SUBCOMM_R, IERR)
    
    !         CALL MPI_GATHERV(SUBGLOBAL_ARRAY, N1_glb*N2_glb*N3_loc, MPI_Datatype, &
    !                          GLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, MPI_Datatype, &
    !                          0, SUBCOMM_R, IERR)
    
    !         ! DEBUG:
    !         ! CALL MCAT(SUBGLOBAL_ARRAY(1,:,:))
    
    !     ELSEIF (axis.EQ.2) THEN
    
    !         ALLOCATE(GLOBAL_ARRAY_COPY(N2_glb,N1_glb,N3_glb))
    
    !         ! CALL MPI_GATHER(SUBGLOBAL_ARRAY,N1_glb*N2_glb*N3_loc,MPI_Datatype, &
    !         !                 GLOBAL_ARRAY_COPY,N1_glb*N2_glb*N3_loc,MPI_Datatype, &
    !         !                 0, SUBCOMM_R, IERR)
    
    !         CALL MPI_GATHERV(SUBGLOBAL_ARRAY, N1_glb*N2_glb*N3_loc, MPI_Datatype, &
    !                         GLOBAL_ARRAY_COPY, RECVCOUNT, DISPLACEMENT, MPI_Datatype, &
    !                         0, SUBCOMM_R, IERR)
    !         GLOBAL_ARRAY = RESHAPE(GLOBAL_ARRAY_COPY,SHAPE(GLOBAL_ARRAY),ORDER = [2,1,3]) ! NEED TO REORDER
    
    !         ! DEBUG:
    !         ! CALL MCAT(SUBGLOBAL_ARRAY(:,3,:))        
    
    !     ENDIF    
    ! ENDIF                 
    
    ! CALL MPI_TYPE_FREE(SUBARRAY_TYPE_resized,IERR)
    ! DEALLOCATE(LOCAL_ARRAY2,SUBGLOBAL_ARRAY,RECVCOUNT,DISPLACEMENT)
    ! END SUBROUTINE ASSEMBLE
    
    end program sample2