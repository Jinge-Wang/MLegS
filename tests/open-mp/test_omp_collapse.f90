! mpifort -cpp -fopenmp test_omp_collapse.f90 -o test_omp_collapse \
!   -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib
program test_omp_collapse
    use omp_lib
    use mpi
    implicit none
    
    integer, parameter :: n = 100, m = 100
    integer :: i, j, count, ierr, rank, size
    integer :: array(n, m), test_array(n, m)
    character(len=100) :: version
    
    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    
    ! Initialize arrays
    array = 0
    test_array = 0
    
    if (rank == 0) then
        write(*,*) "Testing OpenMP collapse clause with non-rectangular loops..."
        
        ! Get compiler version (without preprocessor directives)
#ifdef __GFORTRAN__
        version = "GFortran"
#else
        version = "Unknown compiler"
#endif
        write(*,*) "Compiler:", trim(version)
        
        ! Check OpenMP version
#ifdef _OPENMP
        write(*,*) "OpenMP version:", _OPENMP
#else
        write(*,*) "OpenMP version: Not available"
#endif
    end if
    
    ! Test 1: Standard rectangular loop with collapse
    !$omp parallel do collapse(2) private(i,j)
    do i = 1, n
        do j = 1, m
            array(i,j) = i + j
        end do
    end do
    !$omp end parallel do
    
    ! Serial calculation for verification
    do i = 1, n
        do j = 1, m
            test_array(i,j) = i + j
        end do
    end do
    
    ! Check correctness
    if (all(array == test_array)) then
        if (rank == 0) write(*,*) "Test 1 (rectangular loop): PASSED"
    else
        if (rank == 0) write(*,*) "Test 1 (rectangular loop): FAILED"
    end if
    
    ! Reset arrays
    array = 0
    test_array = 0
    
    ! Test 2: Non-rectangular loop with collapse
    count = 0
    !$omp parallel do collapse(2) private(i,j) reduction(+:count)
    do i = 1, n
        do j = 1, i  ! Non-rectangular loop
            array(i,j) = i + j
            count = count + 1
        end do
    end do
    !$omp end parallel do
    
    ! Serial calculation for verification
    do i = 1, n
        do j = 1, i  ! Non-rectangular loop
            test_array(i,j) = i + j
        end do
    end do
    
    ! Check correctness
    if (all(array(:n,:n) == test_array(:n,:n))) then
        if (rank == 0) write(*,*) "Test 2 (non-rectangular loop): PASSED"
        if (rank == 0) write(*,*) "Total iterations executed:", count
    else
        if (rank == 0) write(*,*) "Test 2 (non-rectangular loop): FAILED"
    end if
    
    call MPI_Finalize(ierr)
    
end program test_omp_collapse