!-------------------------------------------------------------------------------
!> @brief Test program for product and integral operations on scalar objects.
!>
!> This program demonstrates computation of products and integrals using custom
!> scalar types in a mapped Legendre polynomial framework, with parallelization
!> via MPI.
!>
!> Main steps:
!>   - Initializes MPI with serialized threading support.
!>   - Reads input parameters and sets up mapped Legendre polynomials.
!>   - Allocates and initializes scalar objects `a` and `b` to unity.
!>   - Computes:
!>       - Scalar product (`prod`)
!>       - Product over m (`prod_m`)
!>       - Product over k (`prod_k`)
!>       - Product over m and k (`prod_mk`)
!>       - Modified product over m and k (`prod_mk_mod`)
!>       - Integral of scalar object (`inte`)
!>   - Outputs results to standard output if rank is 0.
!>   - Cleans up allocated memory and finalizes MPI.
!>
!> @note
!>   All products should be consistent. The original `prodct_m` and `prodct_k`
!>   include an extra multiplication by `2*TFM%PF(1,1,1)*TFM%NORM(1,1)`, which
!>   is unnecessary for dot products in PFF space. The correct approach uses
!>   Gauss-Legendre quadrature and is independent of radial spectral projection
!>   parameters.
!>
!>   This extra factor compensates for normalization differences between the
!>   original M&M 1997 definition of mapped Legendre polynomials and the
!>   re-normalized version used in this code. In physical space, a spectral
!>   coefficient of `A%E(1,1,1) = 1` corresponds to a function `f(r) = 1/N^0_0`,
!>   where `N^0_0` follows the M&M 1997 definition. M&M 1997 uses a non-normalized
!>   polynomial, yielding `f(r) = 1`, so if a user assigns spectral coefficients
!>   under the impression that he uses the original M&M 1997 un-normalized
!>   polynomials, the product would require a factor of `TFM%NORM(1,1)**2` for
!>   consistency. However, a correct TFM table should always have the correct
!>   function values in the PFF space, so it is still hard to justify the need for
!>   this factor.
!>
!>   Back to the discrepancy in normalization, the reason we use the normalized
!>   version is to allow the use of the same TFM%PF for both forward and backward
!>   spectral projection without additional factors in either step. Without 
!>   re-normalization, an additional division by `(N^m_n)^2` would be needed when
!>   projecting to spectral space.
!-------------------------------------------------------------------------------
program test_prodct
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

implicit none
type(scalar):: a, b
integer:: mm,nn,kk
real(p8):: prod, inte
real(p8), allocatable:: prod_mk(:,:), prod_mk_mod(:,:), prod_m(:), prod_k(:)

! INITITALIZE MPI ENVIRONMENT
CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
   WRITE(*,*) 'The threading support is lesser than that demanded.'
   CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

! READ INPUTS AND INITIALIZE MAPPED LEGENDRE POLYNOMIALS
CALL READCOM('NOECHO')
CALL READIN(5)
CALL LEGINIT()
CALL PRINT_MPI_STRATEGY(FILES%SAVEDIR)
IF (MPI_RANK.EQ.0) CALL SAVEDIN()

! MAIN
call allocate(a, PFF_SPACE)
call allocate(b, PFF_SPACE)

a%E = 0.d0; b%E = 0.d0
if ((a%INTH .eq. 0).and.(a%inx .eq. 0)) then
    a%e(:,1,1) = 1.d0
    b%e(:,1,1) = 1.d0
endif

prod = prodct(a,b) ! multiply tfm%pf(1,1,1)
if (MPI_RANK.EQ.0) then
    WRITE(*,*) 'Product computed successfully.'
    WRITE(*,*) 'Product = ', prod
endif

allocate(prod_m(ntchop), prod_k(nxchop))
prod_m = PRODCTM(a,b)
prod_k = PRODCTK(a,b)
if (MPI_RANK.EQ.0) then
    WRITE(*,*) 'Product (m) computed successfully.'
    WRITE(*,*) 'Product (m) = ', prod_m(:5)
    WRITE(*,*) 'Product (k) computed successfully.'
    WRITE(*,*) 'Product (k) = ', prod_k(:5)
endif

allocate(prod_mk(ntchop, nxchop))
prod_mk = product_mk(a,b) ! < correct
if (MPI_RANK.EQ.0) then
    WRITE(*,*) 'PF(1,1,1) = ', TFM%PF(1,1,1)
    WRITE(*,*) 'NORM(1,1) = ', TFM%NORM(1,1)
    WRITE(*,*) 'Product (mk) computed successfully.'
    WRITE(*,*) 'Product (mk) = ', prod_mk(1,:5)
endif

allocate(prod_mk_mod(ntchop, nxchopdim))
prod_mk_mod = prodct_mk(a,b) ! < correct
if (MPI_RANK.EQ.0) then
    WRITE(*,*) 'Product (mk_mod) computed successfully.'
    WRITE(*,*) 'Product (mk_mod) = ', prod_mk_mod(1,:5)
endif

call rtran(a,-1) ! a to FFF_SPACE
inte = integ(a) ! < correct
if (MPI_RANK.EQ.0) then
    WRITE(*,*) 'Integral computed successfully.'
    WRITE(*,*) 'Integral = ', inte
endif

CALL DEALLOCATE(A)
CALL DEALLOCATE(B)
DEALLOCATE(prod_mk, prod_mk_mod, prod_m, prod_k)
call mpi_barrier(MPI_COMM_IVP,IERR)
call mpi_finalize(IERR)

end program test_prodct