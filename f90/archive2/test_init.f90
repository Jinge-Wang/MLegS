program test_init
    ! PROGRAM TO ASSEMBLE LOCAL ARRAYS INTO GLOBAL ARRAYS
        USE OMP_LIB
        USE MPI
        USE MOD_MISC                                                       ! LEVEL 0
        USE MOD_EIG                                                        ! LEVEL 1
        !USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
        USE MOD_SCALAR3                                                    ! LEVEL 2
        USE MOD_FFT
        USE MOD_LEGOPS
        USE MOD_INIT
        !USE MOD_LAYOUT
    implicit none
    
    ! INPUT:
    TYPE(SCALAR):: A,B,C
    ! OUTPUT:
    COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: GLOBAL_ARRAY
    CHARACTER(LEN=72) :: A_GLB_FILENAME
    
    ! TEMPORARY:
    INTEGER:: N1, N2, N3, I, J, K
    REAL:: RANDNUM_RE, RANDNUM_IM
    
    CALL MPI_INIT(IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

    CALL READCOM('NOECHO')
    CALL READIN(5)
    CALL LEGINIT()
    
    ALLOCATE(GLOBAL_ARRAY(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))

    CALL ALLOCATE(A,FFF_SPACE)
    !CALL RANDOM_SEED
    DO K = 1,SIZE(A%E,3)
        DO J = 1,SIZE(A%E,2)
            DO I = 1,SIZE(A%E,1)

                CALL RANDOM_NUMBER(RANDNUM_RE)
                CALL RANDOM_NUMBER(RANDNUM_IM)
                A%E(I,J,K) = RANDNUM_RE+RANDNUM_IM*IU

            ENDDO
        ENDDO
    ENDDO
    IF ((A%INTH.EQ.0).AND.(A%INX.EQ.0)) A%LN = 2.D0
    

    CALL MSAVE(A, 'A_loc.dat')
    CALL MASSEMBLE(A%E,GLOBAL_ARRAY,1)
    
    ! ============================== OUTPUT ================================
    IF (MPI_RANK.EQ.0) THEN
        A_GLB_FILENAME = 'A_GLB_E.dat'
        CALL MCAT(REAL(A%E(:,1,2)))
        WRITE(*,*) ''

        CALL ALLOCATE(B, FFF_SPACE)
        NULLIFY(B%E)
        ALLOCATE(B%E(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))
        B%E = GLOBAL_ARRAY

        CALL MSAVE(B, A_GLB_FILENAME, .TRUE.)
        CALL MCAT(REAL(B%E(:,1,2)))
        WRITE(*,*) ''
        WRITE(*,*) ''
        ! CALL RANDOM_NUMBER(RANDNUM)
        ! WRITE(*,*) RANDNUM
        ! CALL RANDOM_NUMBER(RANDNUM)
        ! WRITE(*,*) RANDNUM
        CALL ALLOCATE(C, FFF_SPACE)
        CALL MLOAD(A_GLB_FILENAME, C, .TRUE.)
        CALL MCAT(REAL(C%E(:,1,2)))
    ENDIF
    
    CALL DEALLOCATE(A)
    DEALLOCATE(GLOBAL_ARRAY)
    CALL MPI_FINALIZE(IERR)
    
    end program test_init