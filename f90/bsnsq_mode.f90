program bsnsq_mode
! ======================================================================
! NONLINEAR SIMULATION W/ FROZEN BACKGROUND Q-VORTEX +
! LINEAR EIGENMODE PERTURBATION
!
! The initial perturbation is loaded from a perturb.input file
! produced by bsnsq_evp (using the --save-mode flag).
!
! Usage:
!   mpirun -n <N> ./bin/bsnsq_mode_exec [perturb_input_path]
!
!   perturb_input_path  (optional) path to the perturb.input file.
!                       Defaults to './data/perturb.input'.
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
integer:: iii, it
type(scalar):: psi_per, chi_per, b_per
character(len=256) :: perturb_path
integer :: nargs_cmd

call setup_environment('noecho')
call setup_grid(files%savedir)

! Determine path to perturb.input: use command-line argument if
! provided, otherwise fall back to the default location.
nargs_cmd = command_argument_count()
perturb_path = './data/perturb.input'   ! default
if (nargs_cmd > 0) then
   call get_command_argument(1, perturb_path)
endif
if (MPI_RANK == 0) then
   write(*,'(2A)') '  Eigenmode input file: ', trim(perturb_path)
endif

! intialize background shear
call calc_boussi_baseflow(files%savedir)

! set initial perturbation from linear eigenmode
call allocate(psi_per); psi_per%e = 0.d0; psi_per%ln = 0.d0
call allocate(chi_per); chi_per%e = 0.d0; chi_per%ln = 0.d0
call allocate(b_per);     b_per%e = 0.d0;   b_per%ln = 0.d0
call load_eigenmode(perturb_path, psi_per, chi_per, b_per)

! set up monitoring modes
do iii = 1,size(MONITORDATA%MK,1)
   if (MONITORDATA%MK(iii,1).lt.0) MONITORDATA%MK(iii,:) = -MONITORDATA%MK(iii,:)
   if (MONITORDATA%MK(iii,2).lt.0) MONITORDATA%MK(iii,2) = 2*nxchop-1+MONITORDATA%MK(iii,2)
enddo

!> first diagnostic
call diagnost(psi_per,chi_per)
call CALC_BOUSSI_ENERGY(psi_per,chi_per,b_per,tim%t,files%savedir) ! <- requires the baseflow to be initialized
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
!> @brief Load a linear eigenmode perturbation from a perturb.input file
!>
!> Reads the mode profile written by SAVE_PERTURB_BSNSQ (MOD_EVP) and
!> sets the corresponding spectral coefficients in psi_per, chi_per, and
!> b_per. The file format and loading logic follow addperturb_non.f90.
!>
!> The perturb.input contains:
!>   - Number of perturbation modes
!>   - Scale factor (complex) per mode
!>   - For each mode: NDIM, azimuthal m, physical axial wavenumber k,
!>     followed by NDIM lines of spectral coefficients
!>     re(psi), im(psi), re(del2chi), im(del2chi), re(b), im(b)
!>
!> @param[in]    perturb_path  Path to the perturb.input file
!> @param[inout] psi_per       Scalar field for poloidal stream function
!> @param[inout] chi_per       Scalar field for del^2(chi) (toroidal)
!> @param[inout] b_per         Scalar field for buoyancy perturbation
!>
!> @note All MPI ranks read the file independently (same as addperturb_non).
!> @note Conjugate symmetry is enforced for m=0, k>0 modes.
!> @note Coefficients are already scaled by the complex scale factor.
!>
!> @author Jinge Wang
!> @date MAR 2026
subroutine load_eigenmode(perturb_path, psi_per, chi_per, b_per)
! ======================================================================
   implicit none
   character(len=*), intent(in)    :: perturb_path
   type(scalar),     intent(inout) :: psi_per, chi_per, b_per

   integer  :: num, j, jj, ii, nn, err, ios, flag
   real(p8) :: kk, rp1, ip1, rp2, ip2, rp3, ip3
   integer,     allocatable :: mp(:), kp(:)
   complex(p8), allocatable :: perturb_psi(:,:), perturb_chi(:,:), perturb_b(:,:)
   complex(p8), allocatable :: scale(:)
   integer :: mm_local, kk_local

   ! --- Open perturb.input file ---
   open(11, FILE=TRIM(ADJUSTL(perturb_path)), STATUS='OLD', ACTION='READ', IOSTAT=err)
   if (err .ne. 0) then
      if (MPI_RANK .eq. 0) &
         write(*,'(2A)') 'ERROR: load_eigenmode -- cannot open: ', trim(perturb_path)
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
   end if

   ! --- Read number of perturbation modes ---
   read(11,*) num
   allocate(perturb_psi(NRchop, num), perturb_chi(NRchop, num), &
            perturb_b(NRchop, num), mp(num), kp(num), scale(num), STAT=err)
   perturb_psi = 0.0_p8; perturb_chi = 0.0_p8; perturb_b = 0.0_p8

   ! --- Read scale factors ---
   do j = 1, num
      read(11,*) rp1, ip1
      scale(j) = cmplx(rp1, ip1, p8)
   end do

   ! --- Read eigenfunction data for each mode ---
   do j = 1, num
      read(11,*)           ! skip: 'eigenfunction'
      read(11,*)           ! skip: column header (ndim / m / k labels)
      read(11,*) nn, mp(j), kk   ! ndim, azimuthal m, physical axial wavenumber

      ! Skip 8 metadata/header lines (q/h/b, nu params, sigma, coeff header…)
      do jj = 1, 8
         read(11,*)
      end do

      ! Convert physical wavenumber k to 1-based axial index
      kk = zlen / (2.0_p8 * pi) * kk
      kp(j) = nint(kk)

      if (MPI_RANK .eq. 0) then
         write(*,'(A,I0,A,I0,A,2ES13.5)') &
            '  m=', mp(j), '  k=', kp(j), '  scale=', real(scale(j)), aimag(scale(j))
      end if

      ! Sanity checks
      if (abs(real(kp(j), p8) - kk) > 1.0e-6_p8) then
         if (MPI_RANK .eq. 0) then
            write(*,'(A)') 'WARNING: Axial wavenumber is not an integer.'
            write(*,'(A,F14.8)') '               zlen = ', zlen
            write(*,'(A,F14.8)') '         Difference = ', abs(real(kp(j), p8) - kk)
            write(*,'(A)') '         Check that zlen in read.input matches the EVP run.'
         end if
         call mpi_abort(MPI_COMM_WORLD, 1, ierr)
      end if
      if (abs(kp(j)) >= NXCHOP) then
         if (MPI_RANK .eq. 0) &
            write(*,'(A,I0,A)') 'WARNING: |k|=', abs(kp(j)), ' >= NXCHOP. Check NXCHOP in read.input.'
         call mpi_abort(MPI_COMM_WORLD, 1, ierr)
      end if

      ! Handle negative m: use conjugate symmetry
      flag = 0
      if (mp(j) < 0) then
         mp(j) = -mp(j)
         kp(j) = -kp(j)
         flag  = 1
      end if

      ! Shift to 1-based array indices
      mp(j) = mp(j) + 1
      kp(j) = kp(j) + 1
      if (kp(j) <= 0) kp(j) = kp(j) + 2*NXCHOP - 1

      ! Read spectral coefficients: re(psi), im(psi), re(del2chi), im(del2chi), re(b), im(b)
      do ii = 1, nn
         if (ii <= NRchop) then
            read(11,*) rp1, ip1, rp2, ip2, rp3, ip3
            perturb_psi(ii, j) = cmplx(rp1, ip1, p8)
            perturb_chi(ii, j) = cmplx(rp2, ip2, p8)
            perturb_b(ii, j)   = cmplx(rp3, ip3, p8)
         else
            read(11,*)
         end if
      end do

      ! Apply scale factor
      perturb_psi(:, j) = perturb_psi(:, j) * scale(j)
      perturb_chi(:, j) = perturb_chi(:, j) * scale(j)
      perturb_b(:, j)   = perturb_b(:, j)   * scale(j)

      ! For negative-m modes, take complex conjugate
      if (flag == 1) then
         perturb_psi(:, j) = conjg(perturb_psi(:, j))
         perturb_chi(:, j) = conjg(perturb_chi(:, j))
         perturb_b(:, j)   = conjg(perturb_b(:, j))
      end if
   end do
   close(11)

   ! --- Set MPI-local spectral coefficients in the scalar fields ---
   do j = 1, num
      mm_local = mp(j) - psi_per%inth
      kk_local = kp(j) - psi_per%inx

      if (0 < mm_local .and. size(psi_per%e, 2) >= mm_local) then
         if (0 < kk_local .and. size(psi_per%e, 3) >= kk_local) then
            psi_per%e(:NRchop, mm_local, kk_local) = &
               psi_per%e(:NRchop, mm_local, kk_local) + perturb_psi(:NRchop, j)
            chi_per%e(:NRchop, mm_local, kk_local) = &
               chi_per%e(:NRchop, mm_local, kk_local) + perturb_chi(:NRchop, j)
            b_per%e(:NRchop, mm_local, kk_local) = &
               b_per%e(:NRchop, mm_local, kk_local) + perturb_b(:NRchop, j)
         end if

         ! For m=0, k>0: enforce conjugate symmetry at the -k counterpart
         if (mp(j) .eq. 1 .and. kp(j) .gt. 1) then
            kk_local = (2*NXCHOP + 1 - kp(j)) - psi_per%inx
            if (0 < kk_local .and. size(psi_per%e, 3) >= kk_local) then
               psi_per%e(:NRchop, mm_local, kk_local) = &
                  psi_per%e(:NRchop, mm_local, kk_local) + conjg(perturb_psi(:NRchop, j))
               chi_per%e(:NRchop, mm_local, kk_local) = &
                  chi_per%e(:NRchop, mm_local, kk_local) + conjg(perturb_chi(:NRchop, j))
               b_per%e(:NRchop, mm_local, kk_local) = &
                  b_per%e(:NRchop, mm_local, kk_local) + conjg(perturb_b(:NRchop, j))
            end if
         end if
      end if
   end do

   ! Zero out numerical noise below machine-epsilon level
   where (abs(psi_per%e) < 1.0e-24_p8) psi_per%e = 0.0_p8
   where (abs(chi_per%e) < 1.0e-24_p8) chi_per%e = 0.0_p8
   where (abs(b_per%e)   < 1.0e-24_p8)   b_per%e = 0.0_p8

   call chopdo(psi_per)
   call chopdo(chi_per)
   call chopdo(b_per)

   deallocate(perturb_psi, perturb_chi, perturb_b, scale, mp, kp)

   if (MPI_RANK .eq. 0) then
      write(*,'(2A)') '  Eigenmode loaded from: ', trim(perturb_path)
   end if

end subroutine load_eigenmode

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

end program bsnsq_mode
