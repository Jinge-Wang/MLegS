program bsnsq_test
! ======================================================================
! OPEN RE-FINED EIGENVECTORS AND RUN NONLINEAR SIMULATION 
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
integer:: iii,it,mm,kk,pp,II,JJ,FILESTATUS
real(p8),dimension(1,1):: status
real(P8):: time_start,time_end,time0,tmp
type(scalar):: psi_tot,chi_tot,b_per
type(scalar):: dpsi,dchi,db
logical:: file_save = .TRUE.

CALL SETUP_ENVIRONMENT('NOECHO')
CALL SETUP_GRID(FILES%SAVEDIR)

! INITIALIZE FIELDS
CALL allocate(psi_tot); call random_noise(psi_tot)
CALL allocate(chi_tot); call random_noise(chi_tot)
CALL allocate(b_per); call random_noise(b_per, is_save = .true.)
CALL allocate(dpsi)
CALL allocate(dchi)
CALL allocate(db)

! SET UP MONITORING MODES
DO II = 1,SIZE(MONITOR_MK,1)
   IF (MONITOR_MK(II,1).LT.0) MONITOR_MK(II,:) = -MONITOR_MK(II,:)
   IF (MONITOR_MK(II,2).LT.0) MONITOR_MK(II,2) = 2*NXCHOP-1+MONITOR_MK(II,2)
ENDDO
IF (MPI_RANK.EQ.0) THEN
   WRITE(*,*) 'TRACK MODES: (M,AK) - (#MM,#KK)'
   DO II = 1,SIZE(MONITOR_MK,1)
      WRITE(*,92) M(MONITOR_MK(II,1)+1),AK(MONITOR_MK(II,1)+1,MONITOR_MK(II,2)+1),MONITOR_MK(II,1)+1,MONITOR_MK(II,2)+1
   ENDDO
92 FORMAT('(',I3,',',F9.2,') - (#',I3,', #',I3,')')
ENDIF

! INITIAL CONDITIONS
status(1,1)=0
if (MPI_RANK.eq.0) then
   WRITE(*,*) 'PROGRAM STARTED'
   call msave(status, 'status.dat')
   CALL PRINT_REAL_TIME()
endif

!> print psi_tot%e
IF (MPI_RANK.EQ.0) THEN
   WRITE(*,*) 'INITIAL PSI_TOT: '
   DO III = 1,10
      WRITE(*,*) '( MM = ',III,') - '
      CALL MCAT(psi_tot%E(1:10,III,1:10))
   ENDDO
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

!> first diagnostic
call diagnost(psi_tot,chi_tot)
call CALC_ENERGY(psi_tot,chi_tot,b_per,TIM%T,FILES%SAVEDIR)

!> richardson step
call rich(psi_tot,chi_tot,b_per,dpsi,dchi,db)

!> 2nd diagnostic
call diagnost(psi_tot,chi_tot)

!> save initial energy after Richardson step
CALL CALC_ENERGY(psi_tot,chi_tot,b_per,TIM%T,FILES%SAVEDIR)

!> startup
!dpsi and dchi are initially empty, then they are assigned 
!the nonlinear part of the first step
time_start = mpi_wtime()
iii = tim%limit/tim%dt
files%n = 1

do it=1,2

   !> admam-bashforth
   call ADAMSB(psi_tot,chi_tot,b_per,dpsi,dchi,db)

   !> hyperviscosity
   call HYPERV3(psi_tot,chi_tot,b_per)
   call DIAGNOST(psi_tot,chi_tot)

   !> save energy spectrum
   CALL CALC_ENERGY(psi_tot,chi_tot,b_per,TIM%T,FILES%SAVEDIR)

   !> output
   if ((files%t(files%n).le.tim%t) .AND. (file_save)) then

      call msave(psi_tot, TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi(files%n))
      call msave(chi_tot, TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi(files%n))
      call msave(b_per, TRIM(ADJUSTL(FILES%SAVEDIR))//files%b(files%n))
      files%n = files%n + 1

      ! stop the simulation once the last file is saved
      ! if(files%n > files%ne) goto 999
      if(files%n > files%ne) file_save = .FALSE.
      
   endif

enddo

999 continue
time_end = mpi_wtime()

call combine_mode_files(0, NXCHOPDIM, './output')

!> final printout
IF (MPI_RANK.eq.0) THEN
    print *,tim%n,' steps'
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()
    WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

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

end program bsnsq_test
