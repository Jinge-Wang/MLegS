program vort2
! WIPE OUT OTHER MODES: TESTING GROWTH RATE
! -------------------------
   USE MOD_MISC
   USE MOD_BANDMAT
   USE MOD_EIG
   USE MOD_FD
   USE MOD_LIN_LEGENDRE
   USE MOD_SCALAR3
   USE MOD_LAYOUT
   USE MOD_LEGOPS
   USE MOD_MARCH
   USE MOD_DIAGNOSTICS
   USE MOD_INIT
! -------------------------
implicit none
! -------------------------
type(scalar):: psi,chi,psio,chio,dpsi,dchi,psi_o,chi_o
integer:: iii,it,mm,kk,pp, MIND,KIND
real(p8),dimension(1,1):: status
real(p8), allocatable, dimension(:,:):: tempMonitor_mk
real(p8):: time0, tmp, rp1,ip1
REAL(P8), ALLOCATABLE,DIMENSION(:):: ENE_M
character(len=72) :: psiFileName, chiFileName, num


!call readcom('echo')
call readin(5)
!call readcom('noecho')
status(1,1)=0
call msave(status, TRIM(ADJUSTL(FILES%SAVEDIR))//'status.dat')

call allocate(psi)
call allocate(chi)
call allocate(psio)
call allocate(chio)
call allocate(psi_o)
call allocate(chi_o)
call allocate(dpsi)
call allocate(dchi)
call allocate(psio)
call allocate(chio)
ALLOCATE(ENE_M(NTCHOP))

!> initial values
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%psii,psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%chii,chi)

call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> startup : very first time step
call diagnost(psi,chi)

!dpsi and dchi are initially empty, then they are assigned 
!the nonlinear part of the first step 
call rich(psi,chi,dpsi,dchi)
call diagnost(psi,chi)

if(monitorData%start.eq.0) then 
   psi_o = psi
   chi_o = chi
   time0 = tim%t
end if

iii = tim%limit/tim%dt
files%n = 1

do it=1,iii

   if(monitorData%start.eq.tim%n) then 
      ! use this to calculate eigenvalues. this is the psi/chi at t1
      psi_o = psi
      chi_o = chi
      time0 = tim%t
   end if
   
   
   !> integrate
   call adamsb(psi,chi,dpsi,dchi)
   call diagnost(psi,chi)
   
   !> monitor eigenmodes
   call WRITE_ENERGY_DATA_TO_FILE(psi,chi)
!     call WRITE_MODE_DATA_TO_FILE(psi,chi)

   tmp = (real(tim%n)-real(monitorData%start))/real(monitorData%interval)
   if (tmp.ge.1.0_p8 .and. abs(tmp-int(tmp)).lt.1e-08) then
      call egrowth_modified(psi_o,chi_o,psi,chi,time0,tim%t)
      call monitor_eig(psi_o,chi_o,psi,chi,time0,tim%t)
   !   call PRINT_ENERGY_SPECTRUM(psi,chi,1)
      call PRINT_ENERGY_SPECTRUM(psi,chi,0)
   end if
   
   !> freestream adjust
   if(adv%sw.eq.1 .and. mod(it,adv%int).eq.0) then
      call freeadj
      call rich(psi,chi,dpsi,dchi)
      call diagnost(psi,chi)
   endif
   
   !> vorticity removal
   if(rmv%sw.ne.0 .and. mod(it,rmv%int).eq.0) then
      call remove(psi,chi)
      call rich(psi,chi,dpsi,dchi)
      call diagnost(psi,chi)
   endif
   
   !> hyperviscosity adjust
   if(visc%adjsw.ne.0 .and. mod(it,visc%adjint).eq.0) then
      call hypadj(psi,chi,visc%nup)
      call diagnost(psi,chi)
   endif
   
   !> output
   if(files%t(files%n).le.tim%t) then
      call msave(psi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi(files%n))
      call msave(chi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi(files%n))
      files%n = files%n + 1
      if(files%n > files%ne) goto 999
      call mload('status.dat',status)
      if(status(1,1).eq.1.0) goto 999
   endif


   !> Clean modes
   psio%ln = 0.D0 ! Wipe-out other modes
   chio%ln = 0.D0
   psio%e = 0.D0
   chio%e = 0.D0
   psio%e(:,MIND,KIND) = psi%e(:,MIND,KIND)
   chio%e(:,MIND,KIND) = chi%e(:,MIND,KIND)
   ! WRITE(*,*) "ENEMON"
   CALL ENEMON(psio,chio,ENE_M)
   psio%e(:,MIND,KIND) = cmplx(rp1,ip1,p8)*psio%e(:,MIND,KIND)/(sum(ENE_M)**0.5) ! Renormalize 1,1 mode
   chio%e(:,MIND,KIND) = cmplx(rp1,ip1,p8)*chio%e(:,MIND,KIND)/(sum(ENE_M)**0.5)

   psio%e(:,1,1) = psi_o%e(:,1,1) ! Refix 0,0 mode
   chio%e(:,1,1) = chi_o%e(:,1,1)

   psi%e = psio%e
   chi%e = chio%e

enddo

999 continue

print *,tim%n,' steps'
call DEALLOCATE(psio)
call DEALLOCATE(chio)
stop

CONTAINS


end program vort2


