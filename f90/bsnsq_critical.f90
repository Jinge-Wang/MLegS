program bsnsq_critical
! ======================================================================
! BOUSSINESQ SIMULATION WITH QVORTEX AS INITIAL CONDITION
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
type(scalar):: psi_tot,chi_tot,b_per

call setup_environment('noecho')
call setup_grid(files%savedir)

! initialize fields
! call allocate(psi_tot); call allocate(chi_tot);
call calc_qvortex(psi_tot, chi_tot)
call allocate(b_per); b_per%e = 0.d0

! set up monitoring modes
do iii = 1,size(monitor_mk,1)
   if (monitor_mk(iii,1).lt.0) monitor_mk(iii,:) = -monitor_mk(iii,:)
   if (monitor_mk(iii,2).lt.0) monitor_mk(iii,2) = 2*nxchop-1+monitor_mk(iii,2)
enddo

! initial conditions
if (mpi_rank.eq.0) then
   write(*,*) 'program started'
   call print_real_time()
endif

!> first diagnostic
call diagnost(psi_tot,chi_tot)
call CALC_BOUSSI_ENERGY(psi_tot,chi_tot,b_per,tim%t,files%savedir)
call inspect(psi_tot,1)
call inspect(chi_tot,1)

!> initialize solver
CALL PT_SOLVER%INITIALIZE(psi_tot, chi_tot, b_per)
call diagnost(psi_tot,chi_tot)
CALL CALC_BOUSSI_ENERGY(psi_tot,chi_tot,b_per,TIM%T,FILES%SAVEDIR)
call inspect(psi_tot,1)
call inspect(chi_tot,1)
call inspect(b_per,1)

!> startup
iii = tim%limit/tim%dt
files%n = 1
do it=1,iii

   !> time-stepping
   CALL PT_SOLVER%TIME_STEPPING(psi_tot, chi_tot, b_per)
   call inspect(psi_tot,1)
   call inspect(chi_tot,1)
   call inspect(b_per,1)

enddo

999 continue

!> final printout
CALL PT_SOLVER%FINALIZE(psi_tot, chi_tot, b_per)

! ======================================================================
contains
! ======================================================================
!> @brief Inspects and outputs field data to a file for analysis
!>
!> This subroutine examines a scalar field and writes its magnitude
!> squared values to an output file for inspection purposes. The
!> routine only processes distributed fields where both azimuthal 
!> and axial wavenumber starting indices (INTH and INX) are zero.
!>
!> @details The subroutine writes data in a structured format where
!> each line contains the current time followed by radial data for
!> selected vertical modes. The output is appended to 'inspect.output'
!> in the save directory. A MPI barrier ensures synchronization across
!> all processes after execution. Note that field%e is distributed
!> across processors.
!>
!> @param[in] field    Scalar field of type scalar containing the
!>                     distributed data to inspect
!> @param[in] nk       Optional integer specifying number of vertical
!>                     modes to save. If not provided, defaults to
!>                     minimum of 10 or total available modes
!>
!> @note Only processes fields with field%INTH = 0 and field%INX = 0
!> @note Output format: time on first line, then mode index followed
!>       by |field|² values
!> @note Uses MPI barrier for process synchronization
!> @note field%e is distributed across processors
!>
!> @author Jinge WANG
!> @date AUG 2025
subroutine inspect(field, nk)
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
   fname = trim(files%savedir)//'inspect.output'
   open(newunit=iunit, file=trim(fname), status='unknown', action='write', position='append')
   write(iunit,*) TIM%T
   do i = 1, nsave
      write(iunit,'(I4,1x,*(ES14.6,","))') i, (abs(field%e(nn,1,i))**2, nn=1,nrad)
   end do
   close(iunit)

end if
call mpi_barrier(MPI_COMM_IVP,IERR)
return
end subroutine inspect

end program bsnsq_critical
