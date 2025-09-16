program uz_testing
! ======================================================================
! TEST B-UZ CONVERSION
! by JINGE WANG, Aug 2025
! ======================================================================
! [PURPOSE]:
! To test the numerical conversion between a buoyancy perturbation `b`
! and the resulting vertical velocity `uz`. This conversion is a key
! part of the Boussinesq dynamics, where the buoyancy term `-b*z_hat`
! in the momentum equation acts as a source for the poloidal potential
! `chi` via an inverse Laplacian. A crucial aspect of this inversion is
! the handling of a potential logarithmic term in `chi`, which is
! physically required for axisymmetric (m=0), axially-uniform (k=0)
! buoyancy distributions. This test compares two inverse Laplacian
! solvers: `IDEL2`, which requires a pre-calculated log term, and
! `IDEL2LN`, which solves for the log term internally.
! [FINDINGS]:
! 1. The standard solver `IDEL2`, when used without a pre-supplied
!    logarithmic term, fails to correctly invert the buoyancy field.
!    This failure manifests as a significant, unphysical magnitude
!    increase in the highest radial mode of the m=0, k=0 component of
!    the resulting `uz` field. This is the expected numerical artifact.
! 2. The `IDEL2LN` solver, designed to handle this specific case,
!    correctly computes the logarithmic term and produces a physically
!    sound `uz` field without numerical artifacts.
! 3. A manual calculation of the required log term from the initial `uz`
!    field confirms that its value is consistent with the one computed
!    internally by `IDEL2LN`, cross-verifying the implementation. The
!    log term is -ELL2 * TFM%NORM(1,1) * F_000, where f = uz/(1-x)^2.
!    Manually assigning `IDEL2` the correct log term allows it to
!    produce a valid `uz` field.
! [CONCLUSION]:
! Using either the `IDEL2LN` or 'IDEL2' with manually calculated LN term
! leads to the correct result. However, the `IDEL2LN` can be sensitive
! to the last few modes, while manually calculating the LN term involves
! division by (1-TFM%X)^2 in physical space, which may not be as accurate.
! We still need to decide which method to use.s
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
type(scalar):: psi_tot,chi_tot,b_per,b_divxm2_per
type(scalar):: uz_test, rur_test, rup_test, w
type(scalar):: oz_test, ror_test, rop_test
complex(p8), dimension(4,4):: ETD_E, ETD_NL
real(p8):: chiln

call setup_environment('noecho')
call setup_grid(files%savedir)

! initialize fields
call allocate(psi_tot); psi_tot%e = 0.d0
call allocate(chi_tot); chi_tot%e = 0.d0
call allocate(b_per); call gaussian_blob(b_per, 1.0d-3, 0.5d0, 0.0d0, 0.5*zlen, 1.d0, 0.1*zlen)
call allocate(b_divxm2_per); call gaussian_blob_divxm2(b_divxm2_per, 1.0d-3, 0.5d0, 0.0d0, 0.5*zlen, 1.d0, 0.1*zlen)

! set up monitoring modes
do iii = 1,size(MONITORDATA%MK,1)
   if (MONITORDATA%MK(iii,1).lt.0) MONITORDATA%MK(iii,:) = -MONITORDATA%MK(iii,:)
   if (MONITORDATA%MK(iii,2).lt.0) MONITORDATA%MK(iii,2) = 2*nxchop-1+MONITORDATA%MK(iii,2)
enddo

! initial conditions
if (mpi_rank.eq.0) then
   write(*,*) 'program started'
   call print_real_time()
endif

call allocate(uz_test)
uz_test = b_per
call rtran(uz_test,1)
! call inspect_b(uz_test,1) ! PFF
chiln = z2ln(uz_test)
if ((uz_test%inth.eq.0).and.(uz_test%inx.eq.0)) then
   write(*,*) 'chiln = ', chiln
endif

!> b:
tim%t = tim%t + tim%dt
call inspect_b(b_per,1) ! FFF
call allocate(rur_test,PFF_SPACE)
call allocate(rup_test,PFF_SPACE)
rur_test%e = 0.d0
rup_test%e = 0.d0

call allocate(w)
call project(rur_test,rup_test,uz_test,psi_tot,w,0.D0)
call idel2(w,chi_tot,chiln)
if (mpi_rank .eq. 0) write(*,*) 'idel2: chi_tot%ln = ', chi_tot%ln
call chopset(3)
call pc2vel(psi_tot,chi_tot,rur_test,rup_test,uz_test)
call chopset(-3)
call chopdo(uz_test)
call inspect_b(uz_test,1) ! FFF
! call rtran(uz_test,1)
! call inspect_b(uz_test,1) ! PFF

call idel2ln(w,chi_tot)
if (mpi_rank .eq. 0) write(*,*) 'idel2ln: chi_tot%ln = ', chi_tot%ln
call chopset(3)
call pc2vel(psi_tot,chi_tot,rur_test,rup_test,uz_test)
call chopset(-3)
call chopdo(uz_test)
call inspect_b(uz_test,1)

call deallocate(w)
call deallocate(rur_test)
call deallocate(rup_test)
call deallocate(uz_test)

!> b_divxm2
tim%t = tim%t + tim%dt
psi_tot%e = 0.d0
chi_tot%e = 0.d0
if ((b_divxm2_per%INTH.eq.0).and.(b_divxm2_per%INX.eq.0)) then
   write(*,*) '-ell2*b_divxm2_per_(0,0,0)*norm(1,1) = ', -ell2*b_divxm2_per%e(1,1,1)*tfm%norm(1,1)
end if
b_divxm2_per%e = -b_divxm2_per%e
! call inspect_b(b_divxm2_per,1) ! FFF
call idelsqh(b_divxm2_per,chi_tot)
if (mpi_rank .eq. 0) write(*,*) 'chi_tot%ln = ', chi_tot%ln

call allocate(rur_test)
call allocate(rup_test)
call allocate(uz_test)
call chopset(3)
call pc2vel(psi_tot,chi_tot,rur_test,rup_test,uz_test)
call chopset(-3)
call chopdo(uz_test)
call inspect_b(uz_test,1) ! FFF
! call rtran(uz_test,1)
! call inspect_b(uz_test,1) ! PFF

call deallocate(rur_test)
call deallocate(rup_test)
call deallocate(uz_test)

call deallocate(psi_tot)
call deallocate(chi_tot)
call deallocate(b_per)
call deallocate(b_divxm2_per)

call mpi_finalize(ierr)

! ======================================================================
contains
! ======================================================================
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
subroutine random_noise(a, noise_level, is_save)
! ======================================================================
! [USAGE]:
! GENERATES RANDOM NOISE WITH A KOLMOGOROV-LIKE ENERGY SPECTRUM.
! [PARAMETERS]:
! A >> (INOUT) TYPE(SCALAR), THE SCALAR FIELD OBJECT TO BE FILLED WITH NOISE.
!      CONTAINS THE SPECTRAL DATA ARRAY `A%E` AND THE LOCAL
!      INDEX OFFSETS `A%INTH` AND `A%INX`.
! IS_SAVE >> (IN) LOGICAL, OPTIONAL. IF PRESENT AND .TRUE., SAVES THE
!            KOLMOGOROV SCALING FACTOR FOR DEBUGGING.
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
   integer, allocatable :: seed(:), saved_seed(:)
   integer     :: deterministic_seed(1)

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
   allocate(seed(seed_size), saved_seed(seed_size))
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

   deallocate(seed, saved_seed)
    
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

subroutine combine_mode_files(mode_num, k_ind_tot, output_dir)
! ======================================================================
! [USAGE]:
! Combines mode data from individual MPI rank output files into a single
! combined file for a specific mode number. This is primarily used for
! m=0 modes where manual conjugate symmetry handling requires each rank
! to write its portion of the spectrum.
!
! [PARAMETERS]:
! mode_num   (in) : integer, The azimuthal mode number (e.g., 0 for m=0).
! k_ind_tot  (in) : integer, The total number of axial (k) indices in the
!                   global spectrum.
! output_dir (in) : character(len=*), The directory where individual rank
!                   output files are located
! ======================================================================

   implicit none

   ! Input parameters
   integer, intent(in) :: mode_num        ! Mode number (e.g., 0 for m0)
   integer, intent(in) :: k_ind_tot       ! Global ending k index
   character(len=*), intent(in) :: output_dir  ! Directory containing the files
   
   ! Local variables
   integer :: num_procs, temp_k_indices(nxchopdim)
   integer :: rank, row_idx, i, j, k_parse_loop, k_global, ios
   character(len=2048) :: input_filename, output_filename, line
   character(len=32) :: rank_str, mode_str
   character(len=100) :: complex_str
   integer :: input_unit, output_unit
   logical :: file_exists
   character(len=2048) :: header_line
   
   ! Arrays to store data
   complex(kind=8), allocatable :: mode_data(:,:)  ! (N_radial, K_global)
   integer, allocatable :: k_indices(:)            ! k indices for current file
   integer :: num_k_local, num_radial
   real(kind=8) :: real_part, imag_part
   integer :: comma_pos, i_pos, sign_pos
   
   if (MPI_RANK /= 0) return
   call mpi_comm_size(MPI_COMM_IVP, num_procs, ierr)
   
   num_radial = nrchop
   allocate(mode_data(num_radial, k_ind_tot))
   mode_data = (0.0d0, 0.0d0)
   
   ! Create output filename
   write(mode_str, '(I0)') mode_num
   output_filename = trim(output_dir) // '/m' // trim(mode_str) // '_mode_combined.output'
   
   print *, 'Combining mode files for m', mode_num
   print *, 'Output file: ', trim(output_filename)
   
   ! Read data from each rank's output file
   do rank = 0, num_procs - 1
      write(rank_str, '(I3.3)') rank  ! 3-digit rank with leading zeros
      write(mode_str, '(I0)') mode_num
      input_filename = trim(output_dir) // '/m' // trim(mode_str) // '_mode_rank_' // trim(rank_str) // '.output'
      
      ! Check if file exists
      inquire(file=trim(input_filename), exist=file_exists)
      if (.not. file_exists) then
         print *, 'Warning: File does not exist: ', trim(input_filename)
         cycle
      end if
      
      print *, 'Reading file: ', trim(input_filename)
      
      ! Open input file
      input_unit = 100 + rank
      open(unit=input_unit, file=trim(input_filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print *, 'Error opening file: ', trim(input_filename)
         cycle
      end if
      
      ! Read the header line to get k indices for this rank
      read(input_unit, '(A)', iostat=ios) header_line
      if (ios /= 0) then
         print *, 'Error reading header from: ', trim(input_filename)
         close(input_unit)
         cycle
      end if
      
      ! Count how many k indices are in this file
      num_k_local = 0
      do j = 1, nxchopdim
         read(header_line, *, iostat=ios) (temp_k_indices(i), i=1,j)
         if (ios == 0) num_k_local = j
      enddo

      ! First integer is placeholder, rest are k indices
      num_k_local = num_k_local - 1
      if (num_k_local <= 0) then
         close(input_unit)
         cycle
      endif

      allocate(k_indices(num_k_local))
      k_indices = temp_k_indices(2:num_k_local+1)  ! Skip first placeholder

      ! Read num_radial data lines
      do i = 1, num_radial
         read(input_unit, '(A)', iostat=ios) line
         if (ios /= 0) then
            if (ios > 0) then
               print *, 'Error reading data line', i, 'from file:', trim(input_filename)
            end if
            exit  ! End of file or error
         end if
         
         ! Parse the line: "   i real+imagi, real+imagi, ..."
         ! First extract the row index
         read(line, *, iostat=ios) row_idx
         if (ios /= 0) then
            print *, 'Error parsing row index from line:', trim(line)
            cycle
         end if
         
         ! Remove the row index from the line
         comma_pos = index(line, ',')
         line = adjustl(line(comma_pos+1:))
         
         ! Parse each complex number and store in correct global k position
         do j = 1, num_k_local
            comma_pos = index(line, ',')
            if (comma_pos == 0) then
               ! Last complex number (no comma after it)
               complex_str = trim(adjustl(line))
               line = ''
            else
               complex_str = line(1:comma_pos-1)
               line = adjustl(line(comma_pos+1:))
            endif
            
            if (len_trim(complex_str) == 0) exit
            
            ! Parse complex number in format: real+imagi or real-imagi
            i_pos = index(complex_str, 'i')
            if (i_pos == 0) then
               print *, 'Error: no i found in complex number:', trim(complex_str)
               real_part = 0.0d0
               imag_part = 0.0d0
            else
               sign_pos = 0
               do k_parse_loop = i_pos - 1, 1, -1
                  if (complex_str(k_parse_loop:k_parse_loop) == '+' .or. &
                      complex_str(k_parse_loop:k_parse_loop) == '-') then
                     ! Check if it's part of an exponent
                     if (k_parse_loop > 1 .and. (complex_str(k_parse_loop-1:k_parse_loop-1) == 'E')) then
                        continue ! It's an exponent sign, skip it
                     else
                        sign_pos = k_parse_loop
                        exit
                     endif
                  endif
               enddo

               if (sign_pos > 1) then
                  read(complex_str(1:sign_pos-1), *, iostat=ios) real_part
                  read(complex_str(sign_pos:i_pos-1), *, iostat=ios) imag_part
               elseif (sign_pos == 1) then ! Starts with a sign, so it's all imaginary
                  real_part = 0.0_p8
                  read(complex_str(1:i_pos-1), *, iostat=ios) imag_part
               else
                  print *, 'Error parsing complex number:', trim(complex_str)
                  real_part = 0.0d0
                  imag_part = 0.0d0
               endif
               if (ios /= 0) print *, 'Warning: I/O error parsing complex substring: ', trim(complex_str)
            endif
            
            ! Store in the correct global position
            k_global = k_indices(j)
            if (k_global >= 1 .and. k_global <= k_ind_tot) then
               mode_data(i, k_global) = cmplx(real_part, imag_part, kind=8)
            endif
            
            if (len_trim(line) == 0) exit  ! No more data
         enddo
      enddo
      
      deallocate(k_indices)
      close(input_unit)
   enddo
   
   ! Write combined output file
   output_unit = 200
   open(unit=output_unit, file=trim(output_filename), status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      print *, 'Error creating output file: ', trim(output_filename)
      deallocate(mode_data)
      return
   end if
   
   ! Write header line with all k indices (1 to k_ind_tot)
   write(output_unit, '(I4,1x,*(I8,1x))') 0, (j, j=1,k_ind_tot)
   
   ! Write data lines in the same simple format
   do i = 1, num_radial
      write(output_unit, '(I4,1x,*(ES14.6,SP,ES14.6,"i",","))') i, &
         (real(mode_data(i,j)), aimag(mode_data(i,j)), j=1,k_ind_tot)
   enddo
   
   close(output_unit)
   deallocate(mode_data)
   
   print *, 'Successfully created combined file: ', trim(output_filename)

   ! Remove all files in output_dir matching m<mode_num>_mode_rank_*.output
   call system('rm -f '//trim(output_dir)//'/m'//trim(mode_str)//'_mode_rank_*.output')

   return
   
end subroutine combine_mode_files

subroutine gaussian_blob(field, amplitude, r0, phi0, z0, wr, wz)
!=======================================================================
! [USAGE]:
! Initializes a scalar field with a single Gaussian blob. This version
! is corrected for a 2D (r-z) MPI decomposition and populates the
! staggered azimuthal grid required by the pseudo-spectral method.
!
! [PARAMETERS]:
! field    (OUT): The type(SCALAR) object to be initialized.
! amplitude(IN) : The amplitude of the Gaussian.
! r0       (IN) : The radial distance of the blob's center.
! phi0     (IN) : The azimuthal angle (in radians) of the blob's center.
! z0       (IN) : The vertical position of the blob's center.
! wr       (IN) : The characteristic width of the blob in the xy-plane.
! wz       (IN) : The characteristic width of the blob in the z-direction.
!=======================================================================
   implicit none
   type(scalar), intent(out)    :: field
   real(p8), intent(in)        :: amplitude, r0, phi0, z0, wr, wz

   type(scalar) :: b_phys  ! Temporary field in physical space
   integer      :: nn, mm, kk
   real(p8)     :: r_val, z_val, x_val, y_val
   real(p8)     :: phi_val_real, phi_val_imag
   real(p8)     :: x_val_real, y_val_real, x_val_imag, y_val_imag
   real(p8)     :: x0, y0
   real(p8)     :: dist_sq_xy_real, dist_sq_xy_imag, dist_sq_z
   real(p8)     :: real_part, imag_part
   integer      :: local_nr, local_nx

   ! blob's center
   x0 = r0 * cos(phi0)
   y0 = r0 * sin(phi0)

   call allocate(b_phys, ppp_space)
   b_phys%e = 0.d0

   local_nr = size(b_phys%e, 1)
   local_nx = size(b_phys%e, 3)

   !$omp parallel do default(shared) &
   !$omp private(nn, mm, kk, r_val, z_val, phi_val_real, phi_val_imag, &
   !$omp x_val_real, y_val_real, x_val_imag, y_val_imag, &
   !$omp dist_sq_xy_real, dist_sq_xy_imag, dist_sq_z, &
   !$omp real_part, imag_part) collapse(3)
   do kk = 1, local_nx
      do mm = 1, ndimth
         do nn = 1, local_nr
            r_val   = tfm%r(nn + b_phys%inr)
            z_val   = tfm%z(kk + b_phys%inx)
            
            ! --- Calculate value on the primary grid for the REAL part ---
            phi_val_real = tfm%thr(mm)
            x_val_real = r_val * cos(phi_val_real)
            y_val_real = r_val * sin(phi_val_real)
            dist_sq_xy_real = (x_val_real - x0)**2 + (y_val_real - y0)**2
            dist_sq_z = (z_val - z0)**2
            ! real_part = amplitude * exp(-(dist_sq_xy_real / wr**2) - (dist_sq_z / wz**2))
            real_part = amplitude * exp(-(dist_sq_xy_real / wr**2))

            ! --- Calculate value on the staggered grid for the IMAGINARY part ---
            phi_val_imag = tfm%thi(mm)
            x_val_imag = r_val * cos(phi_val_imag)
            y_val_imag = r_val * sin(phi_val_imag)
            dist_sq_xy_imag = (x_val_imag - x0)**2 + (y_val_imag - y0)**2
            ! imag_part = amplitude * exp(-(dist_sq_xy_imag / wr**2) - (dist_sq_z / wz**2))
            imag_part = amplitude * exp(-(dist_sq_xy_imag / wr**2))
            
            b_phys%e(nn, mm, kk) = cmplx(real_part, imag_part, p8)
         end do
      end do
   end do
   !$omp end parallel do

   call chopdo(b_phys)
   call toff(b_phys)
   call allocate(field)
   field = b_phys
   call deallocate(b_phys)

   field%ln = 0.0_p8
   call chopdo(field)

end subroutine gaussian_blob

subroutine gaussian_blob_divxm2(field, amplitude, r0, phi0, z0, wr, wz)
!=======================================================================
! [USAGE]:
! Initializes a scalar field with a single Gaussian blob. This version
! is corrected for a 2D (r-z) MPI decomposition and populates the
! staggered azimuthal grid required by the pseudo-spectral method.
!
! [PARAMETERS]:
! field    (OUT): The type(SCALAR) object to be initialized.
! amplitude(IN) : The amplitude of the Gaussian.
! r0       (IN) : The radial distance of the blob's center.
! phi0     (IN) : The azimuthal angle (in radians) of the blob's center.
! z0       (IN) : The vertical position of the blob's center.
! wr       (IN) : The characteristic width of the blob in the xy-plane.
! wz       (IN) : The characteristic width of the blob in the z-direction.
!=======================================================================
   implicit none
   type(scalar), intent(out)    :: field
   real(p8), intent(in)        :: amplitude, r0, phi0, z0, wr, wz

   type(scalar) :: b_phys  ! Temporary field in physical space
   integer      :: nn, mm, kk
   real(p8)     :: r_val, z_val, x_val, y_val, xm2_val
   real(p8)     :: phi_val_real, phi_val_imag
   real(p8)     :: x_val_real, y_val_real, x_val_imag, y_val_imag
   real(p8)     :: x0, y0
   real(p8)     :: dist_sq_xy_real, dist_sq_xy_imag, dist_sq_z
   real(p8)     :: real_part, imag_part
   integer      :: local_nr, local_nx

   ! blob's center
   x0 = r0 * cos(phi0)
   y0 = r0 * sin(phi0)

   call allocate(b_phys, ppp_space)
   b_phys%e = 0.d0

   local_nr = size(b_phys%e, 1)
   local_nx = size(b_phys%e, 3)

   !$omp parallel do default(shared) &
   !$omp private(nn, mm, kk, r_val, z_val, phi_val_real, phi_val_imag, &
   !$omp x_val_real, y_val_real, x_val_imag, y_val_imag, &
   !$omp dist_sq_xy_real, dist_sq_xy_imag, dist_sq_z, &
   !$omp real_part, imag_part) collapse(3)
   do kk = 1, local_nx
      do mm = 1, ndimth
         do nn = 1, local_nr
            r_val   = tfm%r(nn + b_phys%inr)
            z_val   = tfm%z(kk + b_phys%inx)
            xm2_val = (1-tfm%x(nn + b_phys%inr))**2

            ! --- Calculate value on the primary grid for the REAL part ---
            phi_val_real = tfm%thr(mm)
            x_val_real = r_val * cos(phi_val_real)
            y_val_real = r_val * sin(phi_val_real)
            dist_sq_xy_real = (x_val_real - x0)**2 + (y_val_real - y0)**2
            dist_sq_z = (z_val - z0)**2
            ! real_part = amplitude * exp(-(dist_sq_xy_real / wr**2) - (dist_sq_z / wz**2))
            real_part = amplitude * exp(-(dist_sq_xy_real / wr**2)) / xm2_val

            ! --- Calculate value on the staggered grid for the IMAGINARY part ---
            phi_val_imag = tfm%thi(mm)
            x_val_imag = r_val * cos(phi_val_imag)
            y_val_imag = r_val * sin(phi_val_imag)
            dist_sq_xy_imag = (x_val_imag - x0)**2 + (y_val_imag - y0)**2
            ! imag_part = amplitude * exp(-(dist_sq_xy_imag / wr**2) - (dist_sq_z / wz**2))
            imag_part = amplitude * exp(-(dist_sq_xy_imag / wr**2)) / xm2_val
            
            b_phys%e(nn, mm, kk) = cmplx(real_part, imag_part, p8)
         end do
      end do
   end do
   !$omp end parallel do

   call chopdo(b_phys)
   call toff(b_phys)
   call allocate(field)
   field = b_phys
   call deallocate(b_phys)

   field%ln = 0.0_p8
   call chopdo(field)

end subroutine gaussian_blob_divxm2

subroutine inspect_b(field, nk)
!=======================================================================
! [USAGE]:
! Save the field data of all radial modes for m = 0, and k = 1, nk
!=======================================================================
    implicit none
    type(scalar), intent(in) :: field
    integer, intent(in), optional :: nk
    integer :: iunit, nrad, nsave, i, nn
    character(len=256) :: fname
    if ((field%INTH.eq.0).and.(field%INX.eq.0)) then

    nrad = size(field%e,1)
    if (present(nk)) then
        nsave = min(nk, size(field%e,3))
    else
        nsave = min(10, size(field%e,3))
    end if
    fname = trim(files%savedir)//'inspect_b_output.dat'
    open(newunit=iunit, file=trim(fname), status='unknown', action='write', position='append')
    write(iunit,*) TIM%T
    do i = 1, nsave
      !   write(iunit,'(I4,1x,*(ES14.6,","))') i, (abs(field%e(nn,1,i))**2, nn=1,nrad)
        write(iunit,'(I4,1x,*(ES14.6,SP,ES14.6,"i",","))') i, (field%e(nn,1,i), nn=1,nrad)
    end do
    close(iunit)

    end if
    call mpi_barrier(MPI_COMM_IVP,IERR)
    return
end subroutine inspect_b

end program uz_testing
