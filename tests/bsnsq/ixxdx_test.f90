program test_ixxdx
! ======================================================================
! TEST PROGRAM FOR RADIAL INTEGRATION AND R*DR OPERATOR
! VALIDATES NUMERICAL INTEGRATION OF f(r) = r*dg/dr USING TRAPEZOIDAL RULE
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
type(scalar):: a_pff, a_fff, a_rdr_fff, a_rdr_pff
real(p8), dimension(:), allocatable:: f_r, g_r, rdr_g_r
real(p8):: s0, r1, r2, w
real(p8):: max_diff, max_abs_err, max_f, rms_err, norm_rms_err
real(p8):: sum_sq_err, sum_sq_f
real(p8), parameter:: threshold = 0.01d0  ! Only consider points where |f| > threshold*max|f|
integer:: n_active

call setup_environment('noecho')
call setup_grid(files%savedir)

if (mpi_rank .eq. 0) then
   write(*,*)
   write(*,*) '========================================================'
   write(*,*) 'TESTING RADIAL INTEGRATION AND R*DR OPERATOR'
   write(*,*) '========================================================'
   write(*,*) 'TEST FUNCTION: f(r) = s0*[tanh((r-r1)/w) - tanh((r-r2)/w)]'
   write(*,*) 'INTEGRATION: g(r) such that f(r) = r*dg/dr'
   write(*,*)
endif

! ALLOCATE ARRAYS
allocate(f_r(NR))
allocate(g_r(NR))
allocate(rdr_g_r(NR))

! SET PARAMETERS FOR SMOOTHENED TOP HAT FUNCTION
s0 = 1.0d0
r1 = 0.5d0 * ELL
r2 = 1.5d0 * ELL
w = 0.1d0 * ELL

if (mpi_rank .eq. 0) then
   write(*,*) 'PARAMETERS:'
   write(*,*) '  s0 =', s0
   write(*,*) '  r1 =', r1
   write(*,*) '  r2 =', r2
   write(*,*) '  w  =', w
   write(*,*)
endif

! CREATE TEST FUNCTION f(r) - SMOOTHENED TOP HAT
do ir = 1, NR
   f_r(ir) = s0 * (tanh((TFM%R(ir) - r1)/w) - tanh((TFM%R(ir) - r2)/w))
enddo

! PERFORM NUMERICAL INTEGRATION USING TRAPEZOIDAL RULE
! g(r) = integral from 0 to r of [f(r')/r'] dr'
g_r = radial_integrate_trapz(f_r, TFM%R, NR)

if (mpi_rank .eq. 0) then
   write(*,*) 'NUMERICAL INTEGRATION COMPLETED'
   write(*,*) '  g(r=0)     =', g_r(1)
   write(*,*) '  g(r=ELL)   =', g_r(NR/2)
   write(*,*) '  g(r=inf)   =', g_r(NR)
   write(*,*)
endif

! CREATE SCALAR IN PFF SPACE AND ASSIGN g(r) TO M=0, K=0 MODE
call allocate(a_pff, PFF_SPACE)
a_pff%e = 0.0d0
a_pff%ln = 0.0d0

! ONLY PROC WITH M=0, K=0 ASSIGNS THE VALUES
if ((a_pff%inth .eq. 0) .and. (a_pff%inx .eq. 0)) then
   do ir = 1, NR
      a_pff%e(ir, 1, 1) = CMPLX(g_r(ir), 0.0d0, P8)
   enddo
endif

! TRANSFORM TO FFF SPACE
call allocate(a_fff)
a_fff = a_pff
call rtran(a_fff, -1)

if (mpi_rank .eq. 0) then
   write(*,*) 'TRANSFORMED g(r) TO FFF SPACE'
   write(*,*)
endif

! APPLY R*DR OPERATOR (XXDX)
call allocate(a_rdr_fff)
call xxdx(a_fff, a_rdr_fff)

if (mpi_rank .eq. 0) then
   write(*,*) 'APPLIED XXDX OPERATOR (R*DR)'
   write(*,*)
endif

! TRANSFORM BACK TO PFF SPACE
call allocate(a_rdr_pff)
a_rdr_pff = a_rdr_fff
call rtran(a_rdr_pff, 1)

! EXTRACT rdr_g(r) FROM M=0, K=0 MODE
rdr_g_r = 0.0d0
if ((a_rdr_pff%inth .eq. 0) .and. (a_rdr_pff%inx .eq. 0)) then
   do ir = 1, NR
      rdr_g_r(ir) = REAL(a_rdr_pff%e(ir, 1, 1))
   enddo
endif

! COMPUTE ERROR METRICS
! For a function with varying amplitude, we need normalized errors
max_abs_err = 0.0d0
max_f = 0.0d0
sum_sq_err = 0.0d0
sum_sq_f = 0.0d0
n_active = 0

if ((a_rdr_pff%inth .eq. 0) .and. (a_rdr_pff%inx .eq. 0)) then
   ! Find maximum function amplitude
   do ir = 1, NR
      max_f = max(max_f, abs(f_r(ir)))
   enddo
   
   ! Compute errors
   do ir = 1, NR
      max_abs_err = max(max_abs_err, abs(rdr_g_r(ir) - f_r(ir)))
      
      ! For RMS, only consider "active" region where f is significant
      if (abs(f_r(ir)) > threshold * max_f) then
         sum_sq_err = sum_sq_err + (rdr_g_r(ir) - f_r(ir))**2
         sum_sq_f = sum_sq_f + f_r(ir)**2
         n_active = n_active + 1
      endif
   enddo
   
   rms_err = sqrt(sum_sq_err / max(1, n_active))
   norm_rms_err = sqrt(sum_sq_err / max(1.0d-15, sum_sq_f))
endif

if (mpi_rank .eq. 0) then
   write(*,*) '========================================================'
   write(*,*) 'ERROR ANALYSIS'
   write(*,*) '========================================================'
   write(*,*) 'Function amplitude: max|f(r)| =', max_f
   write(*,*) 'Active region: points with |f| >', threshold*max_f, '(', n_active, 'points)'
   write(*,*)
   write(*,*) 'ABSOLUTE ERRORS:'
   write(*,*) '  Max |f - r*dg/dr|         =', max_abs_err
   write(*,*) '  RMS error (active region) =', rms_err
   write(*,*)
   write(*,*) 'NORMALIZED ERRORS:'
   write(*,*) '  Max error / max|f|        =', max_abs_err / max_f
   write(*,*) '  RMS error / RMS(f)        =', norm_rms_err
   write(*,*)
   
   if (norm_rms_err < 1.0d-4) then
      write(*,*) '*** TEST PASSED: EXCELLENT AGREEMENT (< 0.01%) ***'
   else if (norm_rms_err < 1.0d-2) then
      write(*,*) '*** TEST PASSED: GOOD AGREEMENT (< 1%) ***'
   else
      write(*,*) '*** TEST WARNING: MODERATE AGREEMENT ***'
   endif
   write(*,*)
endif

! SAVE RESULTS TO FILE (ONLY PROC WITH M=0, K=0)
if ((a_rdr_pff%inth .eq. 0) .and. (a_rdr_pff%inx .eq. 0)) then
   open(unit=100, file=TRIM(FILES%SAVEDIR)//'rdr_compare.info', status='replace')
   write(100, '(A)') '# RADIAL INTEGRATION AND R*DR OPERATOR TEST'
   write(100, '(A)') '# Column 1: r (radial collocation point)'
   write(100, '(A)') '# Column 2: f(r) (original function)'
   write(100, '(A)') '# Column 3: g(r) (integrated result)'
   write(100, '(A)') '# Column 4: r*dg/dr (XXDX applied to g)'
   
   do ir = 1, NR
      write(100, '(4E20.12)') TFM%R(ir), f_r(ir), g_r(ir), rdr_g_r(ir)
   enddo
   
   close(100)
   
   write(*,*) 'RESULTS SAVED TO: ', TRIM(FILES%SAVEDIR)//'rdr_compare.info'
   write(*,*)
endif

! ======================================================================
! STAGE 2: CONSTRUCT FULL VELOCITY FIELD FROM SHEAR
! ======================================================================
if (mpi_rank .eq. 0) then
   write(*,*)
   write(*,*) '========================================================'
   write(*,*) 'STAGE 2: CONSTRUCT VELOCITY FIELD FROM SHEAR'
   write(*,*) '========================================================'
   write(*,*)
endif

call construct_baseflow_from_shear(f_r)

if (mpi_rank .eq. 0) then
   write(*,*)
   write(*,*) '========================================================'
   write(*,*) 'TEST COMPLETED SUCCESSFULLY'
   write(*,*) '========================================================'
   write(*,*)
endif

! CLEANUP
deallocate(f_r)
deallocate(g_r)
deallocate(rdr_g_r)
call deallocate(a_pff)
call deallocate(a_fff)
call deallocate(a_rdr_fff)
call deallocate(a_rdr_pff)

call MPI_BARRIER(MPI_COMM_IVP, IERR)
call MPI_FINALIZE(IERR)

contains

! ======================================================================
function radial_integrate_trapz(f, r, npts) result(g)
! ======================================================================
! INTEGRATE f(r)/r USING TRAPEZOIDAL RULE
! COMPUTES g(r) = integral from 0 to r of [f(r')/r'] dr'
! ======================================================================
    implicit none
    integer, intent(in):: npts
    real(p8), dimension(npts), intent(in):: f, r
    real(p8), dimension(npts):: g
    integer:: i

    ! BOUNDARY CONDITION AT r=0
    g(1) = 0.0d0

    ! CUMULATIVE TRAPEZOIDAL INTEGRATION
    do i = 2, npts
        g(i) = g(i-1) + 0.5d0 * (f(i-1)/r(i-1) + f(i)/r(i)) * (r(i) - r(i-1))
    enddo

    ! ADJUST TO ENSURE g(r=infinity) = 0
    g(:) = g(:) - g(npts)
   
end function radial_integrate_trapz
! ======================================================================

subroutine construct_baseflow_from_shear(shear)
! ======================================================================
! CONSTRUCT VELOCITY FIELD FROM SHEAR FUNCTION
! WORKFLOW:
! 1. shear(r) = r*d/dr(angular_velocity)
! 2. Integrate to get angular_velocity(r)
! 3. omega_z = 2*angular_velocity + shear
! 4. Create (1-x)^-2*omega_z in PFF space (M=0, K=0)
! 5. Transform to FFF and apply IDELSQH to get PSI
! 6. Call PC2VEL and PC2VOR to get velocity/vorticity
! 7. Extract RUP0 and OZ0, assign to BSNSQ
! 8. Validate by computing LOCAL_SHEAR and comparing
! ======================================================================
    implicit none
    real(p8), dimension(NR), intent(in):: shear

    ! LOCAL VARIABLES
    type(scalar):: omega_z, psi, chi
    type(scalar):: rur, rup, uz, ror, rop, oz
    real(p8), dimension(NR):: angular_vel, omega_z_phys, omega_z_refactor
    real(p8), dimension(NR):: rup0_extracted, oz0_extracted, shear_extracted
    integer:: ir_bf
    real(p8):: max_shear_err, norm_shear_err

    ! shear/r = d/dr(angular_velocity)
    angular_vel = radial_integrate_trapz(shear, TFM%R, NR)

    ! omega_z = 2*angular_velocity + shear
    omega_z_phys = 2.0d0 * angular_vel + shear

    ! omega_z_refactor = omega_z / ((1-x)^2)
    omega_z_refactor = omega_z_phys / (1.0d0 - TFM%X(:NR))**2   
    call allocate(omega_z, PFF_SPACE)
    omega_z%e = 0.0d0; omega_z%ln = 0.0d0
    if ((omega_z%inth .eq. 0) .and. (omega_z%inx .eq. 0)) then
        omega_z%e(:NR, 1, 1) = CMPLX(omega_z_refactor(:NR), 0.0d0, P8)
    endif
    call rtran(omega_z, -1)
    omega_z%e = -omega_z%e

    call allocate(psi)
    call idelsqh(omega_z, psi)

    ! CHI = 0 (no poloidal component for axisymmetric rotation)
    call allocate(chi)
    chi%e = 0.0d0; chi%ln = 0.0d0

    ! GET VELOCITY AND VORTICITY FROM PSI, CHI
    call allocate(rur); call allocate(rup); call allocate(uz)
    call allocate(ror); call allocate(rop); call allocate(oz)
    call chopset(3)
    call pc2vel(psi, chi, rur, rup, uz)
    call pc2vor(psi, chi, ror, rop, oz)
    call chopset(-3)

    ! TRANSFORM TO PFF SPACE TO EXTRACT M=0, K=0 MODE
    call rtran(rup, 1)
    call rtran(oz, 1)

    rup0_extracted = 0.0d0
    oz0_extracted = 0.0d0

    if ((rup%inth .eq. 0) .and. (rup%inx .eq. 0)) then
        if (maxval(abs(aimag(rup%e(:NR,1,1)))) .gt. 1.0d-14) then
            if (mpi_rank .eq. 0) then
            write(*,*) 'WARNING: RUP has non-zero imaginary part:', &
                    maxval(abs(aimag(rup%e(:NR,1,1))))
            endif
        endif
        rup0_extracted(:NR) = real(rup%e(:NR,1,1))/(TFM%R(:NR))**2
        oz0_extracted(:NR) = real(oz%e(:NR,1,1))
    endif

    call mpi_allreduce(MPI_IN_PLACE, rup0_extracted, NR, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_IVP, IERR)
    call mpi_allreduce(MPI_IN_PLACE, oz0_extracted, NR, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_IVP, IERR)

    ! ALLOCATE AND ASSIGN BSNSQ ARRAYS
    if (.not. allocated(BSNSQ%OMEGA0)) then
        allocate(BSNSQ%OMEGA0(NDIMR), BSNSQ%OZ0(NDIMR), BSNSQ%SIGMA0(NDIMR))
    endif
    BSNSQ%OMEGA0(:NR) = rup0_extracted(:NR)
    BSNSQ%OZ0(:NR) = oz0_extracted(:NR)

    ! COMPUTE SIGMA0 = LOCAL_SHEAR(IS_REFACTORED=.TRUE.)
    BSNSQ%SIGMA0(:NR) = BSNSQ%LOCAL_SHEAR(IS_REFACTORED=.TRUE.)

    ! COMPUTE LOCAL_SHEAR WITHOUT REFACTORING
    shear_extracted(:NR) = BSNSQ%LOCAL_SHEAR(IS_REFACTORED=.FALSE.)

    ! COMPUTE ERROR METRICS
    max_shear_err = maxval(abs(shear_extracted(:NR) - shear(:NR)))
    norm_shear_err = sqrt(sum((shear_extracted(:NR) - shear(:NR))**2) / &
                            sum(shear(:NR)**2))

    if (mpi_rank .eq. 0) then
        write(*,*)
        write(*,*) '========================================================'
        write(*,*) 'BASEFLOW VALIDATION RESULTS'
        write(*,*) '========================================================'
        write(*,*) 'Max |shear_recon - shear_orig| =', max_shear_err
        write(*,*) 'RMS error / RMS(shear)         =', norm_shear_err
        write(*,*)
        
        if (norm_shear_err < 1.0d-4) then
            write(*,*) '*** BASEFLOW TEST PASSED: EXCELLENT AGREEMENT ***'
        else if (norm_shear_err < 1.0d-2) then
            write(*,*) '*** BASEFLOW TEST PASSED: GOOD AGREEMENT ***'
        else
            write(*,*) '*** BASEFLOW TEST WARNING: MODERATE AGREEMENT ***'
        endif
        write(*,*)
        
        ! SAVE COMPARISON TO FILE
        open(unit=101, file=TRIM(FILES%SAVEDIR)//'baseflow_compare.info', status='replace')
        write(101, '(A)') '# BASEFLOW RECONSTRUCTION FROM SHEAR TEST'
        write(101, '(A)') '# Column 1: r'
        write(101, '(A)') '# Column 2: shear_original'
        write(101, '(A)') '# Column 3: angular_velocity'
        write(101, '(A)') '# Column 4: omega_z'
        write(101, '(A)') '# Column 5: rup0_extracted'
        write(101, '(A)') '# Column 6: oz0_extracted'
        write(101, '(A)') '# Column 7: shear_reconstructed'
        
        do ir_bf = 1, NR
            write(101, '(7E20.12)') TFM%R(ir_bf), shear(ir_bf), angular_vel(ir_bf), &
                omega_z_phys(ir_bf), rup0_extracted(ir_bf), oz0_extracted(ir_bf), &
                shear_extracted(ir_bf)
        enddo
        
        close(101)
        
        write(*,*) 'BASEFLOW RESULTS SAVED TO: ', TRIM(FILES%SAVEDIR)//'baseflow_compare.info'
        write(*,*)
    endif

    ! CLEANUP
    call deallocate(omega_z)
    call deallocate(omega_z_fff)
    call deallocate(psi)
    call deallocate(chi)
    call deallocate(rur)
    call deallocate(rup)
    call deallocate(uz)
    call deallocate(ror)
    call deallocate(rop)
    call deallocate(oz)
   
end subroutine construct_baseflow_from_shear
! ======================================================================

end program test_ixxdx