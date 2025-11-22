program test_shear
! ======================================================================
! TEST PROGRAM FOR SHEAR PROFILE INJECTION VIA CALC_BOUSSI_BASEFLOW
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
   USE MOD_BOUSSINESQ
   USE MOD_MARCH
   USE MOD_DIAGNOSTICS
   USE MOD_INIT
! -------------------------
implicit none
! -------------------------
integer:: ir
real(p8), dimension(:), allocatable:: shear_r
real(p8):: s0, r1, r2, w

call setup_environment('noecho')
call setup_grid(files%savedir)

if (mpi_rank .eq. 0) then
   write(*,*)
   write(*,*) '========================================================'
   write(*,*) 'TEST SHEAR PROFILE INJECTION'
   write(*,*) '========================================================'
   write(*,*)
endif

! ALLOCATE ARRAY
allocate(shear_r(NR))

! SET PARAMETERS FOR SMOOTHENED TOP HAT SHEAR FUNCTION
s0 = 1.0d0
r1 = 0.5d0 * ELL
r2 = 1.5d0 * ELL
w = 0.1d0 * ELL

if (mpi_rank .eq. 0) then
   write(*,*) 'SHEAR FUNCTION PARAMETERS:'
   write(*,*) '  s0 =', s0
   write(*,*) '  r1 =', r1
   write(*,*) '  r2 =', r2
   write(*,*) '  w  =', w
   write(*,*)
endif

! CREATE TEST SHEAR FUNCTION - SMOOTHENED TOP HAT
do ir = 1, NR
   shear_r(ir) = s0 * (tanh((TFM%R(ir) - r1)/w) - tanh((TFM%R(ir) - r2)/w))
enddo

! SAVE ORIGINAL SHEAR PROFILE
if (mpi_rank .eq. 0) then
   open(unit=100, file=TRIM(FILES%SAVEDIR)//'shear_input.info', status='replace')
   write(100, '(A)') '# INPUT SHEAR PROFILE'
   write(100, '(A)') '# Column 1: r (radial collocation points)'
   write(100, '(A)') '# Column 2: shear(r) (prescribed shear profile)'
   
   do ir = 1, NR
      write(100, '(2E20.12)') TFM%R(ir), shear_r(ir)
   enddo
   
   close(100)
   
   write(*,*) 'INPUT SHEAR SAVED TO: ', TRIM(FILES%SAVEDIR)//'shear_input.info'
   write(*,*)
endif

! USE CALC_BOUSSI_BASEFLOW TO CONSTRUCT ALL BASEFLOW QUANTITIES
! THIS WILL:
! 1. INTEGRATE SHEAR TO GET ANGULAR VELOCITY
! 2. CHECK BOUNDEDNESS
! 3. COMPUTE OMEGA0, OZ0, SIGMA0
! 4. SAVE TO baseflow.dat
call CALC_BOUSSI_BASEFLOW(FILES%SAVEDIR, SHEAR_PROFILE=shear_r)

if (mpi_rank .eq. 0) then
   write(*,*)
   write(*,*) '========================================================'
   write(*,*) 'TEST COMPLETED'
   write(*,*) '========================================================'
   write(*,*)
   write(*,*) 'OUTPUT FILES:'
   write(*,*) '  ', TRIM(FILES%SAVEDIR)//'shear_input.info'
   write(*,*) '  ', TRIM(FILES%SAVEDIR)//'baseflow.dat'
   write(*,*)
endif

! CLEANUP
deallocate(shear_r)

call MPI_BARRIER(MPI_COMM_IVP, IERR)
call MPI_FINALIZE(IERR)

end program test_shear
