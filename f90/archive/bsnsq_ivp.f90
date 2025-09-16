program bsnsq_ivp
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
   USE MOD_EVP, only: SAVE_PERTURB
! -------------------------
implicit none
! -------------------------
integer:: iii,it,mm,kk,pp,II,JJ,FILESTATUS
real(p8),dimension(1,1):: status
real(P8):: time_start,time_end,time0,tmp
type(scalar):: psi_tot,chi_tot,b_per
type(scalar):: dpsi,dchi,db
logical:: file_save = .TRUE.

! REFINEMENT
character(LEN=200):: FILENAME

! INITIALIZE MPI
CALL SETUP_ENVIRONMENT('NOECHO')
CALL SETUP_GRID(FILES%SAVEDIR)

! INITIALIZE FIELDS
CALL allocate(psi_tot)
CALL allocate(chi_tot)
CALL allocate(b_per)
CALL allocate(dpsi)
CALL allocate(dchi)
CALL allocate(db)

! SET UP MONITORING MODES
DO II = 1,SIZE(MONITORDATA%MK,1)
   IF (MONITORDATA%MK(II,1).LT.0) MONITORDATA%MK(II,:) = -MONITORDATA%MK(II,:)
   IF (MONITORDATA%MK(II,2).LT.0) MONITORDATA%MK(II,2) = 2*NXCHOP-1+MONITORDATA%MK(II,2)
ENDDO
IF (MPI_RANK.EQ.0) THEN
   WRITE(*,*) 'TRACK MODES: (M,AK) - (#MM,#KK)'
   DO II = 1,SIZE(MONITORDATA%MK,1)
      WRITE(*,92) M(MONITORDATA%MK(II,1)+1),AK(MONITORDATA%MK(II,1)+1,MONITORDATA%MK(II,2)+1),MONITORDATA%MK(II,1)+1,MONITORDATA%MK(II,2)+1
   ENDDO
92 FORMAT('(',I3,',',F9.2,') - (#',I3,', #',I3,')')
ENDIF

! INITIAL CONDITIONS
status(1,1)=0
if (MPI_RANK.eq.0) then
   WRITE(*,*) 'PROGRAM STARTED'
   call msave(status, 'status.dat')
   CALL PRINT_REAL_TIME()
   WRITE(*,*) 'psi/chi/b filename (leave empty for default): '
   READ(*,10) FILENAME
   if (LEN_TRIM(FILENAME) == 0) then
      WRITE(*,*) 'Using default filename. '
   else
      WRITE(*,*) 'Using provided filename: ', TRIM(FILENAME)
   endif
endif
CALL MPI_BCAST(FILENAME,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)
! Check if filename is empty and set default if needed
if (LEN_TRIM(FILENAME) == 0) then
   call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//files%PSI0,psi_tot)
   call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//files%CHI0,chi_tot)
   call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//files%B0,b_per)
else
   call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FILENAME))//"_psi.dat",psi_tot)
   call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FILENAME))//"_chi.dat",chi_tot)
   call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FILENAME))//"_b.dat",b_per)
endif

!> first diagnostic
call diagnost(psi_tot,chi_tot)
call PRINT_ENERGY_SPECTRUM(psi_tot,chi_tot,1)

!> richardson step
call RICH_BOUSSI_FE_BE(psi_tot,chi_tot,b_per,dpsi,dchi,db)

!> 2nd diagnostic
call diagnost(psi_tot,chi_tot)
! call PRINT_ENERGY_SPECTRUM(psi_tot,chi_tot,1)

!> save initial energy after Richardson step
CALL CALC_BOUSSI_ENERGY(psi_tot,chi_tot,b_per,TIM%T,FILES%SAVEDIR)

!> startup
!dpsi and dchi are initially empty, then they are assigned 
!the nonlinear part of the first step
time_start = mpi_wtime()
iii = tim%limit/tim%dt
files%n = 1

do it=1,iii

   !> adams-bashforth
   call STEP_BOUSSI_AB_CN(psi_tot,chi_tot,b_per,dpsi,dchi,db)

   !> hyperviscosity
   call HYPERV3(psi_tot,chi_tot,b_per)
   call DIAGNOST(psi_tot,chi_tot)

   !> save energy spectrum
   CALL CALC_BOUSSI_ENERGY(psi_tot,chi_tot,b_per,TIM%T,FILES%SAVEDIR)

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

end program bsnsq_ivp
