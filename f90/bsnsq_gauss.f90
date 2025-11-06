program bsnsq_gauss
! ======================================================================
! NONLINEAR SIMULATION W/ FROZEN BACKGROUND Q-VORTEX + GAUSSIAN VORTEX
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
integer:: iii,it
type(scalar):: psi_per,chi_per,b_per

call setup_environment('noecho')
call setup_grid(files%savedir)

! set background flow as q-vortex defined in read.input
! set initial condition as a gaussian vortex perturbation
call allocate(psi_per)
call allocate(chi_per) 
call allocate(b_per)
call initialize_gaussian_vortex_example(psi_per, chi_per, b_per)

! set up monitoring modes
do iii = 1,size(MONITORDATA%MK,1)
   if (MONITORDATA%MK(iii,1).lt.0) MONITORDATA%MK(iii,:) = -MONITORDATA%MK(iii,:)
   if (MONITORDATA%MK(iii,2).lt.0) MONITORDATA%MK(iii,2) = 2*nxchop-1+MONITORDATA%MK(iii,2)
enddo

!> first diagnostic
call diagnost(psi_per,chi_per)
call CALC_BOUSSI_ENERGY(psi_per,chi_per,b_per,tim%t,files%savedir)
call inspect(psi_per,1)
call inspect(chi_per,1)
call inspect(b_per,1)

!> initialize solver
CALL PT_SOLVER%INITIALIZE(psi_per, chi_per, b_per)
call diagnost(psi_per,chi_per)
CALL CALC_BOUSSI_ENERGY(psi_per,chi_per,b_per,TIM%T,FILES%SAVEDIR)

!> startup
iii = tim%limit/tim%dt
files%n = 1
do it=1,iii

   !> time-stepping
   CALL PT_SOLVER%TIME_STEPPING(psi_per, chi_per, b_per)
   
enddo

!> final printout
CALL PT_SOLVER%FINALIZE(psi_per, chi_per, b_per)

! ======================================================================
contains
! ======================================================================
!> @brief Saves mode data to a file for tracking purposes
!>
!> This subroutine writes complex field data for specific azimuthal and
!> axial mode numbers to a tracking file. The data is appended to an
!> existing file or creates a new one if it doesn't exist.
!>
!> @details The output file is named 'mode_track_mmm_kkk.dat' where mmm
!> and kkk are the azimuthal and axial mode numbers respectively. Each
!> line contains the time followed by the complex values of psi_new and
!> chi_new arrays.
!>
!> @param[in] t        Current time value
!> @param[in] psi_new  Complex array containing psi field data
!> @param[in] chi_new  Complex array containing chi field data  
!> @param[in] m_i      Azimuthal mode number for file naming
!> @param[in] k_i      Axial mode number for file naming
!>
!> @note The format uses scientific notation with 6 digits precision
!> @note File is opened in append mode to preserve existing data
!>
!> @author Jinge Wang
!> @date AUG 2025
subroutine save_mode(t,psi_new,chi_new,m_i,k_i)
! ======================================================================
   complex(p8),DIMENSION(:):: psi_new,chi_new
   real:: t
   integer:: m_i,k_i
   integer:: nn

   open(UNIT=777,FILE=TRIM(ADJUSTL(FILES%SAVEDIR))//'mode_track_'//ITOA3(m_i)//'_'//ITOA3(k_i)//'.dat',&
   &STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')

   WRITE(777,320) t,psi_new(1:size(psi_new)),chi_new(1:size(psi_new))

   close(777)
320 FORMAT(F10.3,',',(S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))
end subroutine save_mode

! ======================================================================
!> @brief Generates random noise with Kolmogorov-like energy spectrum
!>
!> This subroutine fills a scalar field with random noise that follows
!> a Kolmogorov-like energy spectrum. The routine handles MPI distribution
!> and ensures proper conjugate symmetry for m=0 modes across different
!> processor ranks.
!>
!> @details The noise generation uses an exponential decay based on total
!> wavenumber. For m=0 modes, deterministic seeding ensures conjugate
!> symmetry between positive and negative axial modes across MPI ranks.
!> The decay rate and amplitude can be controlled through parameters.
!>
!> @param[inout] a            Scalar field object to be filled with noise.
!>                           Contains spectral data array a%e and local
!>                           index offsets a%inth and a%inx
!> @param[in]    noise_level  Optional noise amplitude (default: 1.0e-6)
!> @param[in]    is_save      Optional flag to save debugging output files
!>                           for m=0 modes if present and .true.
!>
!> @note Uses MPI barriers for synchronization
!> @note Conjugate symmetry is manually enforced for m=0 modes
!> @note Random seed is based on system clock and MPI rank
!> @note Debugging files are saved as 'm0_mode_rank_XXX.output'
!>
!> @author Jinge WANG  
!> @date AUG 2025
subroutine random_noise(a, noise_level, is_save)
! ======================================================================
   implicit none

   type(scalar), intent(inout)   :: a
   real(p8), optional :: noise_level
   logical, optional, intent(in) :: is_save

   real(p8)    :: rand_real, rand_imag, magnitude, k2_total
   real(p8)    :: k_radial_norm, k_azimuthal_norm, k_axial_norm
   real(p8)    :: noise_amplitude, decay_rate
   ! real(p8), allocatable, dimension(:,:,:):: kolm_data
   integer     :: nnn, mmm, kkk, k_ind
   real(p8)    :: global_k, global_m, kolm_factor
   integer     :: seed_size, clock, base_seed
   integer, allocatable :: seed(:), saved_seed(:), deterministic_seed(:)

   if (present(noise_level)) then
      noise_amplitude = noise_level
   else
      noise_amplitude = 1.0e-6_p8
   endif

   base_seed = 987654321
   decay_rate = 12.0_p8

   call mpi_barrier(MPI_COMM_IVP,IERR)
   ! for all m>0 modes
   call random_seed(size=seed_size)
   allocate(seed(seed_size), saved_seed(seed_size), deterministic_seed(seed_size))
   call system_clock(count=clock)
   seed = clock + 37 * MPI_RANK ! Unique seed per rank
   call random_seed(put=seed)

   ! allocate(kolm_data(size(a%e,1), size(a%e,2), size(a%e,3)))
   ! kolm_data = 0.d0
   do mmm = 1, size(a%e, 2)
      global_m = real(m(mmm + a%inth))
      do kkk = 1, size(a%e, 3)
         global_k = ak(mmm + a%inth, kkk + a%inx)

         if (global_m == 0.0) then ! manual conjugate symmetry is required
               
               call random_seed(get=saved_seed)

               ! Create a deterministic seed from the global axial mode number ONLY.
               ! This seed is identical for the rank owning (+k) and the rank owning (-k).
               deterministic_seed = 0
               deterministic_seed(1) = base_seed + 1000 * abs(nint(global_k * 100.0_p8))
               call random_seed(put=deterministic_seed)
               do nnn = 1, size(a%e, 1)
                  call random_number(rand_real)
                  call random_number(rand_imag)
                  rand_real = (rand_real - 0.5_p8) * 2.0_p8
                  rand_imag = (rand_imag - 0.5_p8) * 2.0_p8

                  k_azimuthal_norm = 0.0_p8
                  k_axial_norm     = abs(global_k)/maxval(abs(ak))
                  k_radial_norm    = real(nnn) / real(nrchop)
                  k2_total    = k_radial_norm**2 + k_azimuthal_norm**2 + k_axial_norm**2

                  if (k2_total > 1.0e-12_p8) then
                     kolm_factor = 0.1*exp(-decay_rate*k2_total)
                     magnitude = noise_amplitude * kolm_factor
                     ! kolm_data(nnn, mmm, kkk) = magnitude
                  else
                     magnitude = 0.0_p8
                  endif
                  
                  ! Apply conjugate logic based on the sign of k
                  if (abs(global_k) < 1.0e-12_p8) then ! n=0 mode
                     a%e(nnn, mmm, kkk) = cmplx(rand_real * magnitude, 0.0_p8, kind=p8)
                  else if (global_k > 0.0_p8) then ! n>0 modes
                     a%e(nnn, mmm, kkk) = cmplx(rand_real * magnitude,  rand_imag * magnitude, kind=p8)
                  else ! n<0 modes
                     a%e(nnn, mmm, kkk) = cmplx(rand_real * magnitude, -rand_imag * magnitude, kind=p8)
                  endif
               end do

               call random_seed(put=saved_seed)

         else ! m > 0
               do nnn = 1, size(a%e, 1)
                  k_azimuthal_norm = global_m/maxval(m)
                  k_axial_norm     = abs(global_k)/maxval(abs(ak))
                  k_radial_norm    = real(nnn) / real(nrchop)
                  k2_total    = k_radial_norm**2 + k_azimuthal_norm**2 + k_axial_norm**2

                  if (k2_total > 1.0e-12_p8) then
                     kolm_factor = 0.1*exp(-decay_rate*k2_total)
                     magnitude = noise_amplitude * kolm_factor
                     ! kolm_data(nnn, mmm, kkk) = magnitude
                  else
                     magnitude = 0.0_p8
                  endif
                  
                  call random_number(rand_real)
                  call random_number(rand_imag)
                  rand_real = (rand_real - 0.5_p8) * 2.0_p8
                  rand_imag = (rand_imag - 0.5_p8) * 2.0_p8
                  a%e(nnn, mmm, kkk) = cmplx(rand_real * magnitude, rand_imag * magnitude, kind=p8)
               end do
         endif
      end do ! k loop

      if ((global_m == 0.0).and.(present(is_save))) then
         open(unit=777, file='./output/m0_mode_rank_'//trim(itoa3(MPI_RANK))//'.output', status='unknown', action='write')
         ! Write a placeholder for alignment, then the global k indices
         write(777, '(I4,1x,*(I8,1x))') 0, (k_ind+a%INX, k_ind=1,size(a%e,3))
         do nnn = 1, size(a%e, 1)
            write(777, '(I4,",",*(ES14.6,SP,ES14.6,"i",","))') nnn, &
               (real(a%e(nnn,mmm,k_ind)), imag(a%e(nnn,mmm,k_ind)), k_ind=1,size(a%e,3))
         end do
         close(777)
      endif

   end do ! m loop

   deallocate(seed, saved_seed, deterministic_seed)
    
   ! ! save kolmogorov factor for each rank
   ! if (present(is_save)) then
   !    open(unit=778, file='./output/kolmo_factor_rank_'//trim(itoa3(MPI_RANK))//'.output', status='unknown', action='write')
   !    do mmm = 1, size(kolm_data, 2)
   !       do kkk = 1, size(kolm_data, 3)
   !          do nnn = 1, size(kolm_data, 1)
   !             write(778, '(I4,1x,I4,1x,F12.6,1x,ES14.6)') nnn, m(mmm+a%INTH), ak(mmm+a%inth,kkk+a%INX), kolm_data(nnn,mmm,kkk)
   !          end do
   !       end do
   !    end do
   !    close(778)
   ! endif

   ! deallocate(kolm_data)
   call mpi_barrier(MPI_COMM_IVP,ierr)
   call chopdo(a)
    
   if (MPI_RANK == 0) then
      a%ln = 0.d0
      write(*,*) 'Random noise with Kolmogorov spectrum generated.'
   endif

end subroutine random_noise

! ======================================================================
!> @brief Extract local angular velocity from physical r*u_theta field
!>
!> @details Computes the angular velocity Ω = (r*u_θ)/r² at a specified
!> radial location r0, averaged over all azimuthal and axial positions.
!> Handles MPI-distributed data using collective operations:
!> - Each rank checks if it owns the target radial coordinate
!> - Local averaging over theta (mm) and z (kk) indices
!> - MPI_ALLREDUCE to gather contributions from all ranks
!>
!> @param[in]  r0        Target radial coordinate for extraction
!> @param[in]  rup_phys  r*u_theta field in PHYSICAL space (PPP_SPACE)
!> @param[out] omega0    Extracted angular velocity (averaged over θ, z)
!>
!> @note rup_phys must be in physical (PPP) space with MPI pencil decomposition
!> @note Averages real and imaginary parts: 0.5*(THR + THI) for each point
!> @note Aborts if target radius not found in any rank (should never happen)
!>
!> @author Jinge Wang  
!> @date Oct 2025
subroutine get_omega_from_rup(r0, rup_phys, omega0)
  implicit none
  real(p8), intent(in) :: r0
  type(scalar), intent(in) :: rup_phys
  real(p8), intent(out) :: omega0
  
  integer :: nn, mm, kk, nn_local, nn_global
  real(p8) :: rg, rup_val, omega_sum, cnt
  
  if (rup_phys%space .ne. PPP_SPACE) then
    if (MPI_RANK == 0) then
      write(*,*) 'Error: rup_phys must be in physical space (PPP_SPACE).'
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
    endif
  endif

  ! Average over theta and z at the radial location (MPI-distributed)
  omega_sum = 0.0_p8
  cnt = 0.0_p8

  nn_global = minloc(abs(tfm%r - r0), 1)
  rg = tfm%r(nn_global)
  if (abs(rg) .le. 1.0e-12_p8) then
    if (MPI_RANK == 0) then
      write(*,*) 'Error: Target radius r0 is too close to zero.'
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
    endif
  endif

  nn_local = nn_global - rup_phys%inr
  if (nn_local .ge. 1 .and. nn_local .le. size(rup_phys%e, 1)) then
    do kk = 1, size(rup_phys%e,3)
      do mm = 1, size(rup_phys%e,2)
        ! Average real and imaginary parts (THR and THI values)
        rup_val = 0.5_p8 * (real(rup_phys%e(nn_local,mm,kk), p8) + &
                            aimag(rup_phys%e(nn_local,mm,kk)))
        omega_sum = omega_sum + rup_val / (rg**2)
        cnt = cnt + 1.0_p8
      enddo
    enddo
  endif
  
  ! Gather contributions from all MPI ranks
  call mpi_allreduce(MPI_IN_PLACE, omega_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_IVP, ierr)
  call mpi_allreduce(MPI_IN_PLACE, cnt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_IVP, ierr)
  
  ! Compute average
  if (cnt > 0.0_p8) then
    omega0 = omega_sum / cnt
  else
    call mpi_abort(MPI_COMM_WORLD, 1, ierr)
  endif
  
end subroutine get_omega_from_rup

! ======================================================================
!> @brief Initialize quasi-equilibrium Gaussian vortex perturbation
!>
!> @details Creates a balanced Gaussian vortex in physical space using
!> cyclogeostrophic balance. The vortex is characterized by:
!> - Azimuthal velocity: u_theta ~ r_local * exp(-r_local^2/L^2 - z_local^2/H^2)
!> - Radial & Axial velocity: u_r, u_z = 0 (quasi-balanced state)
!> - Buoyancy: b ~ z_local * exp(-r_local^2/L^2 - z_local^2/H^2)
!> 
!> The vertical scale H is determined from cyclostrophic balance and
!> stratification constraints. The routine extracts background angular
!> velocity omega0 from the input rup_bg field at the vortex center.
!>
!> @param[in]  center    Vortex center [r0, theta0, z0] in cylindrical coords
!> @param[in]  Ro        Rossby number (controls vortex strength)
!> @param[in]  L         Horizontal length scale (pressure anomaly radius)
!> @param[in]  Nc2_in    Squared Brunt-Väisälä frequency of the vortex core
!> @param[in]  rup_bg    Background r*u_theta field in PHYSICAL space (PPP_SPACE)
!> @param[out] rur       r*u_r perturbation in physical space
!> @param[out] rup       r*u_theta perturbation in physical space  
!> @param[out] uz        u_z perturbation in physical space (zero)
!> @param[out] b         Buoyancy perturbation in physical space
!>
!> @note rup_bg must be in physical (PPP) space, not spectral (FFF) space.
!> @note Output fields are in physical space with MPI pencil decomposition.
!> @note Physical space convention: complex value = cmplx(THR_value, THI_value).
!> @note Uses BSNSQ%OMEGA and BSNSQ%BV0 for background rotation and stratification.
!>
!> @author Jinge Wang
!> @date Oct 2025
subroutine initialize_gaussian_vortex( &
    center, Ro, L, Nc2, &
    rur, rup, uz, b)

  implicit none
  
  ! Input parameters
  real(p8), intent(in) :: center(3)  ! (r0, theta0, z0)
  real(p8), intent(in) :: Ro, L
  real(p8), intent(in) :: Nc2
  
  ! Output perturbation fields
  type(scalar), intent(inout) :: rur, rup, uz, b
  
  ! Local variables
  real(p8) :: r0, th0, z0, omega0, H, f_eff, f, N, x0, y0
  real(p8) :: rg, zg, thr_g, thi_g, xg, yg, rh, zh, gauss_exp, v_th_local
  real(p8) :: v_r_glob, v_th_glob, cos_phi, sin_phi
  real(p8) :: rur_r, rup_r, b_r, rur_i, rup_i, b_i
  integer :: nn, mm, kk
  
  ! Extract center and get background state
  r0 = center(1); th0 = center(2); z0 = center(3)
  f = 2.0_p8 * BSNSQ%OMEGA
  N = BSNSQ%BV0
  
  ! old: extract omega0 from rup(scalar) using theta-z averaging
  ! call get_omega_from_rup(r0, rup_bg_phys, omega0)

  ! Directly use BSNSQ%RUP0 to extract omega0
  nn = minloc(abs(tfm%r - r0), 1)
  if (abs(tfm%r(nn) - r0) .le. 1.0e-12_p8) then
    call mprint('Error: Vortex center r0 is too close to zero.')
    call mpi_abort(MPI_COMM_WORLD, 1, ierr)
  endif
  omega0 = BSNSQ%RUP0(nn)/(tfm%r(nn)**2)

  if (mpi_rank == 0) then
    write(*,'(A)') '=== Gaussian Vortex ==='
    write(*,'(A,3F8.3)') 'Center (r,θ,z): ', r0, th0, z0
    write(*,'(A,F8.3)') 'Background Ω: ', omega0
    write(*,'(A,F8.3)') 'Coriolis f: ', f
    write(*,'(A,F8.3)') 'Strat N: ', N
    write(*,'(A,F8.3)') 'Rossby number: ', Ro
    write(*,'(A,F8.3)') 'Input Nc^2:    ', Nc2
    write(*,'(A,2F8.3)') 'Scales (L,H): ', L, H
  endif

  call calc_vortex_params(Ro, L, omega0, f, N, Nc2, H)
  
  if (mpi_rank == 0) then
    write(*,'(A,F8.3)') 'Aspect H/L: ', H/L
  endif
  
  ! Physical parameters
  f_eff = f + 2.0_p8 * omega0
  x0 = r0 * cos(th0)
  y0 = r0 * sin(th0)
  
  ! Allocate and initialize fields
  if (rur%space .ne. PPP_SPACE) then
    call mprint('Error: Output fields must be in physical space (PPP_SPACE).')
    call mpi_abort(MPI_COMM_IVP, ERR_FLAGS%SOLVER, ierr)
  endif
  
  ! Loop over MPI-local physical space: e(nn,mm,kk) with offsets inr,inx (inth=0 in PPP)
  !$omp parallel do collapse(2) private(nn,mm,kk,rg,zg,thr_g,thi_g, &
  !$omp xg,yg,rh,zh,gauss_exp,v_th_local,v_r_glob,v_th_glob, &
  !$omp cos_phi,sin_phi,rur_r,rup_r,b_r,rur_i,rup_i,b_i)
  do kk = 1, size(rup%e,3)
    do mm = 1, size(rup%e,2)
      do nn = 1, size(rup%e,1)
        
        if (kk + rup%inx .gt. ndimx) cycle
        if (nn + rup%inr .gt. ndimr) cycle

        ! Global grid coordinates (MPI pencil: full theta, distributed r and z)
        rg = tfm%r(nn + rup%inr)
        zg = tfm%z(kk + rup%inx)
        thr_g = tfm%thr(mm)
        thi_g = tfm%thi(mm)
        
        ! === THR component (real part) ===
        xg = rg * cos(thr_g); yg = rg * sin(thr_g)
        rh = sqrt((xg - x0)**2 + (yg - y0)**2); zh = zg - z0
        gauss_exp = exp(-(rh/L)**2 - (zh/H)**2)
        
        ! Set velocity to zero if too close to vortex center or axis
        if (rh < 0.01_p8 * L .or. rg < 1.0e-12_p8) then
          rur_r = 0.0_p8; rup_r = 0.0_p8
        else
          v_th_local = (H**2 * (N**2 - Nc2) / f_eff) * (-2.0_p8 * rh / L**2) * gauss_exp
          v_r_glob = v_th_local * (-(yg - y0) / rh) * cos(thr_g) + v_th_local * ((xg - x0) / rh) * sin(thr_g)
          v_th_glob = -v_th_local * (-(yg - y0) / rh) * sin(thr_g) + v_th_local * ((xg - x0) / rh) * cos(thr_g)
          rur_r = rg * v_r_glob; rup_r = rg * v_th_glob
        endif
        b_r = 2.0_p8 * zh * (N**2 - Nc2) * gauss_exp
        
        ! === THI component (imaginary part) ===
        xg = rg * cos(thi_g); yg = rg * sin(thi_g)
        rh = sqrt((xg - x0)**2 + (yg - y0)**2)
        gauss_exp = exp(-(rh/L)**2 - (zh/H)**2)
        
        ! Set velocity to zero if too close to vortex center or axis
        if (rh < 0.01_p8 * L .or. rg < 1.0e-12_p8) then
          rur_i = 0.0_p8; rup_i = 0.0_p8
        else
          v_th_local = (H**2 * (N**2 - Nc2) / f_eff) * (-2.0_p8 * rh / L**2) * gauss_exp
          v_r_glob = v_th_local * (-(yg - y0) / rh) * cos(thi_g) + v_th_local * ((xg - x0) / rh) * sin(thi_g)
          v_th_glob = -v_th_local * (-(yg - y0) / rh) * sin(thi_g) + v_th_local * ((xg - x0) / rh) * cos(thi_g)
          rur_i = rg * v_r_glob; rup_i = rg * v_th_glob
        endif
        b_i = 2.0_p8 * zh * (N**2 - Nc2) * gauss_exp
        
        ! Store in physical space convention
        rur%e(nn,mm,kk) = rur%e(nn,mm,kk) + cmplx(rur_r, rur_i, p8)
        rup%e(nn,mm,kk) = rup%e(nn,mm,kk) + cmplx(rup_r, rup_i, p8)
        uz%e(nn,mm,kk) = uz%e(nn,mm,kk) + cmplx(0.0_p8, 0.0_p8, p8)
        b%e(nn,mm,kk) = b%e(nn,mm,kk) + cmplx(b_r, b_i, p8)
        
      enddo
    enddo
  enddo
  !$omp end parallel do
  
  call chopdo(rur); call chopdo(rup); call chopdo(uz); call chopdo(b)
  if (mpi_rank == 0) write(*,'(A,3F8.3)') 'Vortex initialized at ', r0, th0, z0
  
end subroutine initialize_gaussian_vortex

! ======================================================================
!> @brief Calculate vortex physical parameters from balance equations
!>
!> @details Computes the vertical scale H for a Gaussian vortex using
!> the cyclostrophic and hydrostatic balance derived from Hassanzadeh et al. (2012),
!> adapted for a non-static background flow.
!>
!> @param[in]  Ro      Rossby number (controls vortex strength and sign)
!> @param[in]  L       Horizontal length scale (pressure anomaly radius)
!> @param[in]  omega0  Background angular velocity at vortex center
!> @param[in]  f       System Coriolis parameter (= 2*BSNSQ%OMEGA)
!> @param[in]  N       Background Brunt-Väisälä frequency
!> @param[in]  Nc2     Squared Brunt-Väisälä frequency of the vortex core
!> @param[out] H       Vertical length scale (computed from balance)
!>
!> @note The aspect ratio formula accounts for Rv != L for a Gaussian model.
!> @note A positive alpha^2 is required for a physically realizable vortex.
!>
!> @author Jinge Wang
!> @date Oct 2025
subroutine calc_vortex_params(Ro, L, omega0, f, N, Nc2, H)
  implicit none
  real(p8), intent(in) :: Ro, L, omega0, f, N, Nc2
  real(p8), intent(out) :: H
  
  real(p8) :: f_eff, f_eff_over_f, L_over_Rv, alpha2
  
  f_eff = f + 2.0_p8 * omega0
  
  ! For a vortex derived from a Gaussian pressure anomaly, the radius
  ! of maximum velocity (Rv) is related to the pressure scale (L) by Rv = L/sqrt(2).
  ! f_eff_over_f = f_eff / f
  L_over_Rv = sqrt(2.0_p8)
  
  ! Denominator must be non-zero for a baroclinic vortex
  if (abs(Nc2 - N**2) < 1.0e-12_p8) then
     if (mpi_rank == 0) write(*,'(A)') 'Error: N^2 and Nc^2 are equal. Cannot form baroclinic vortex.'
     H = 0.0_p8
     return
  endif

  ! Aspect ratio from the corrected balance equation
  ! alpha2 = (Ro * (f_eff_over_f + Ro * L_over_Rv) * f**2) / (Nc2 - N**2)
  alpha2 = (Ro * (1 + Ro * L_over_Rv) * f_eff**2) / (Nc2 - N**2)
  
  if (alpha2 <= 0.0_p8) then
    if (mpi_rank == 0) write(*,'(A,E12.4)') 'Warning: alpha^2 is negative or zero: ', alpha2
    if (mpi_rank == 0) write(*,'(A)') '         Vortex cannot exist in this configuration. Setting H=0.'
    H = 0.0_p8
    call mpi_abort(MPI_COMM_IVP, ERR_FLAGS%SOLVER, ierr)
  endif
  
  H = sqrt(alpha2) * L
  
end subroutine calc_vortex_params

! ======================================================================
!> @brief Example wrapper for Gaussian vortex initialization
!> 
!> @details Creates a quasi-equilibrium Gaussian vortex superposed on a
!> background q-vortex flow. The workflow is:
!> 1. Generate background q-vortex (psi_bg, chi_bg) using QVORTEX
!> 2. Extract local omega at vortex center
!> 3. Initialize Gaussian vortex perturbation in physical space
!> 4. Add perturbation to background velocity fields
!> 5. Project combined velocity back to (psi, chi) spectral representation
!> @param[out] psi_field  Toroidal streamfunction (background + perturbation)
!> @param[out] chi_field  Poloidal streamfunction (background + perturbation)  
!> @param[out] b_field    Buoyancy field (perturbation only)
!>
!> @note Vortex center is set off-axis to trigger Zombie Vortex Instability
!> @note Background field from QVORTEX must be compatible with domain
!>
!> @author Jinge Wang
!> @date Oct 2025
subroutine initialize_gaussian_vortex_example(psi_field, chi_field, b_field)
  implicit none
  type(scalar), intent(inout) :: psi_field, chi_field, b_field
  
  real(p8) :: center(3)
  real(p8), parameter :: Ro = -0.3_p8  ! Rossby number
  real(p8), parameter :: L = 0.3_p8  ! Horizontal length scale
  real(p8), parameter :: Nc2 = 0.0_p8  ! Core stratification
  type(scalar) :: rur, rup, uz, b_pert, w
  
  if (mpi_rank == 0) write(*,'(A)') 'Setting up Gaussian vortex...'
  
  center = [5.0_p8, PI/4.0_p8, zlen/2.0_p8]

  ! Step 1: Initialize background q-vortex
  ! NOTE: we only use the m=0, k=0 component for background flow
  call calc_boussi_baseflow(files%savedir)
  
  ! Step 2: Initialize Gaussian vortex perturbation in physical space
  ! This extracts angular velocity of background flow at center, computes 
  ! vortex parameters, and creates perturbation velocity and buoyancy fields
  call allocate(rur, PPP_SPACE); rur%e = cmplx(0.0_p8, 0.0_p8, kind=p8)
  call allocate(rup, PPP_SPACE); rup%e = cmplx(0.0_p8, 0.0_p8, kind=p8)
  call allocate(uz, PPP_SPACE); uz%e = cmplx(0.0_p8, 0.0_p8, kind=p8)
  call allocate(b_pert, PPP_SPACE); b_pert%e = cmplx(0.0_p8, 0.0_p8, kind=p8)
  call initialize_gaussian_vortex(center, Ro, L, Nc2, rur, rup, uz, b_pert)

  ! ! Step 2.5: Add mirroring vortex
  ! center = [5.0_p8, PI/4.0_p8 + PI, zlen/2.0_p8]
  ! call initialize_gaussian_vortex(center, Ro, L, Nc2, rur, rup, uz, b_pert)
  
  ! ! DEBUG: Save perturbation velocity fields for inspection
  ! call save_perturbation_velocity(rur, rup, uz, b_pert)

  ! Step 3: Project velocity perturbations back to (psi, chi) spectral space
  call allocate(w)
  call project(rur, rup, uz, psi_field, w, ln=0.0_p8)
  call deallocate(rur); call deallocate(rup); call deallocate(uz)
  call idel2ln(w, chi_field)
  call deallocate(w)

  ! Step 4: Assign perturbation density
  call toff(b_pert)
  call chopdo(b_pert)
  b_field = b_pert
  call deallocate(b_pert)

  ! Step 5: Smooth perturbation fields
  call hyp3init(.true.)
  call hyperv3(psi_field, chi_field, b_field, DT_SCALE = 2.0_p8) ! two exponential filters at highest
  call hyp3free()

  if (mpi_rank == 0) then
    write(*,'(A)') 'Gaussian vortex setup complete.'
  endif  
  
end subroutine initialize_gaussian_vortex_example

end program bsnsq_gauss
