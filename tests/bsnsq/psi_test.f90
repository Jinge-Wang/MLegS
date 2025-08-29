program psi_test
! ======================================================================
! TEST VEL-PSI PROJECTION
! by JINGE WANG, Aug 2025
! ======================================================================
! [PURPOSE]:
! I noticed that the code was having the following issue:
! For a psi-chi field, constructing a velocity field using PC2VEL, and 
! then immediately projecting them back to psi and chi using PROJECT and
! IDEL2(or IDEL2LN), the psi-field is totally off for all the spectral
! coefficients of the (m = 0, k = 0) mode. Chi calculation is fine.
! The test below shows the cause and the solution.
!
! [FINDINGS]:
! 1. As said in Matushima paper starting Eqn (60), to obtain the poloidal
!     toroidal components of a vector field, one first finds the log
!     term of psi using the relation that psi%ln = - rup_00(r->infty)/2.
!     The contribution from this log term is then removed from the vector
!     field, and one will use Eqn (62~65) to obtain the main psi and 
!     del2chi. Project operator (before the fix) did not consider the 
!     log term.
!
! 2. The reason we never considered the log term for psi is because, in
!     the original AB-CN scheme with incompressible fluid, project is 
!     only used for projection (U X W) into PSIN and CHIN, and one can
!     show that PSIN will always have a zero log term. This is not true
!     for the Boussinesq approximations. Nor it holds true is we add
!     background rotation.
!
! 3. One might be tempted to use the value of rup_00(NR) in PFF space
!     to infer the log term contribution. However, this approach is
!     not accurate since r(NR) is not really at infinity. With our
!     mapping, r->infty corresponds to x = 1, so we must evaluate rup
!     at x = 1. This can be done by recognizing the fact that the
!     un-normalized Legendre Polynomials (m = 0, k = 0) all equal 1
!     at x = 1, so we simply need to sum the un-normalized spectral
!     coefficients of rup_00 to obtain its limit at r -> infty. The
!     only minor complication is that we need to ensure the correct
!     normalization is applied to the spectral coefficients before
!     summing them.
!
! [CONCLUSION]:
! We create VEL2PLN subroutine to exactly calculate the psi log term,
! and remove its contribution from rup to allow further projection. 
! The test below demonstrates that it correctly reconstructs the psi 
! field and the subsequent velocity fields with machine precision. The 
! test also demonstrates that using rup_00(NR) is not as accurate.
!
! VEL2PLN and the modified PROJECT are now part of MOD_LEGOPS.
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
type(scalar):: uz_test, rur_test, rup_test, w
complex(p8), dimension(4,4):: ETD_E, ETD_NL

call setup_environment('noecho')
call setup_grid(files%savedir)

! initialize fields
! a. gaussian vortex
! call calc_qvortex(psi_tot, chi_tot)
! call allocate(b_per); b_per%e = 0.d0
! b. psi with non-zero log
call allocate(psi_tot); psi_tot%e = 0.d0; psi_tot%ln = 1.d0
call allocate(chi_tot); chi_tot%e = 0.d0
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
call inspect(b_per,1)

! ======================================================================
!> psi log test
tim%t = tim%t + tim%dt
if ((psi_tot%inth.eq.0).and.(psi_tot%inx.eq.0)) then
   write(*,*) '(original) psi_ln = ', psi_tot%ln
endif

call allocate(rur_test); call allocate(rup_test); call allocate(uz_test)
call pc2vel(psi_tot,chi_tot,rur_test,rup_test,uz_test)

call tofp(rur_test)
call tofp(rup_test)
call tofp(uz_test)
call velmax(rur_test,rup_test,uz_test)

!> use original project
call allocate(w)
call project(rur_test,rup_test,uz_test,psi_tot,w)
call idel2ln(w,chi_tot)
call deallocate(w)

call inspect(psi_tot,1)
call inspect(chi_tot,1)
call inspect(b_per,1)

!> check if -.5*rup_test(rmax,00) = psi_l
call mpi_barrier(mpi_comm_ivp,ierr)
if ((rup_test%inth.eq.0).and.(rup_test%inx.eq.0)) then ! rup in PFF space
   write(*,*) '(project) psi_ln = ', psi_tot%ln ! psi%ln is not modified in project
   write(*,*) 'rup_test to psi_ln = ', -rup_test%e(nr,1,1)/2.d0
endif

!> calculate psi%ln and remove its contribution for rup
tim%t = tim%t + tim%dt
call allocate(w)
call project(rur_test,rup_test,uz_test,psi_tot,w,.true.)
call idel2ln(w,chi_tot)
call deallocate(w)

if ((psi_tot%inth.eq.0).and.(psi_tot%inx.eq.0)) then
   write(*,*) '(new) psi_ln = ', psi_tot%ln
endif

call inspect(psi_tot,1)
call inspect(chi_tot,1)
call inspect(b_per,1)

call pc2vel(psi_tot,chi_tot,rur_test,rup_test,uz_test)
call tofp(rur_test)
call tofp(rup_test)
call tofp(uz_test)
call velmax(rur_test,rup_test,uz_test)

call deallocate(rur_test)
call deallocate(rup_test)
call deallocate(uz_test)

! ======================================================================
! !> project test
! call allocate(rur_test); call allocate(rup_test); call allocate(uz_test)
! call pc2vel(psi_tot,chi_tot,rur_test,rup_test,uz_test)

! call tofp(rur_test)
! call tofp(rup_test)
! call tofp(uz_test)
! call velmax(rur_test,rup_test,uz_test)

! call allocate(w)
! call project(rur_test,rup_test,uz_test,psi_tot,w)
! call idel2ln(w,chi_tot)
! call deallocate(w)

! !> check if -.5*rup_test(rmax,00) = psi_l
! call rtran(oz_test,1)
! call mpi_barrier(mpi_comm_ivp,ierr)
! if ((oz_test%inth.eq.0).and.(oz_test%inx.eq.0)) then ! both in PFF space
!    write(*,*) 'psi_ln = ', psi_tot%ln
!    write(*,*) 'rup_test to psi_ln = ', -rup_test%e(nr,1,1)/2.d0
! endif

! call inspect(psi_tot,1)
! call inspect(chi_tot,1)
! call inspect(b_per,1)

! call pc2vel(psi_tot,chi_tot,rur_test,rup_test,uz_test)
! call tofp(rur_test)
! call tofp(rup_test)
! call tofp(uz_test)
! call velmax(rur_test,rup_test,uz_test)

! call deallocate(rur_test)
! call deallocate(rup_test)
! call deallocate(uz_test)

! ======================================================================
! !> half step test
! !> initialize etd operators
! call CALC_BOUSSI_ETD_OP(etd_e, etd_nl, tim%dt)
! etd_nl = 0.d0
! etd_e = 0.d0
! do iii=1,4
!    etd_e(iii,iii) = 1.d0
! end do

! !> first half step
! tim%t = tim%t + tim%dt
! call STEP_BOUSSI_ETDFE_BE(psi_tot, chi_tot, b_per, etd_e, etd_nl)
! call CALC_BOUSSI_ENERGY(psi_tot,chi_tot,b_per,tim%t,files%savedir)
! call inspect(psi_tot,1)
! call inspect(chi_tot,1)
! call inspect(b_per,1)
! if (mpi_rank.eq.0) then
!    write(*,*) 'b_per%ln = ', b_per%ln
!    write(*,*) 'psi_tot%ln = ', psi_tot%ln
!    write(*,*) 'chi_tot%ln = ', chi_tot%ln
! endif
! call anynan(b_per,'b_per')
! call anynan(psi_tot,'psi_tot')
! call anynan(chi_tot,'chi_tot')

! !> second half step
! tim%t = tim%t + tim%dt
! call STEP_BOUSSI_ETDFE_BE(psi_tot, chi_tot, b_per, etd_e, etd_nl)
! call CALC_BOUSSI_ENERGY(psi_tot,chi_tot,b_per,tim%t,files%savedir)
! call inspect(psi_tot,1)
! call inspect(chi_tot,1)
! call inspect(b_per,1)
! if (mpi_rank.eq.0) then
!    write(*,*) 'b_per%ln = ', b_per%ln
!    write(*,*) 'psi_tot%ln = ', psi_tot%ln
!    write(*,*) 'chi_tot%ln = ', chi_tot%ln
! endif
! call anynan(b_per,'b_per')
! call anynan(psi_tot,'psi_tot')
! call anynan(chi_tot,'chi_tot')

!> finalization
call mprint('half step completed')
call deallocate(psi_tot)
call deallocate(chi_tot)
call deallocate(b_per)
call mpi_finalize(ierr)

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

end program psi_test
