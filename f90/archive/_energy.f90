program viewEnergy
  !     read psi, chi
  !     write velocity, vorticity
  !    ----------------------------------------------------
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
  !    ----------------------------------------------------
  implicit none
  !    ----------------------------------------------------
  type(scalar):: psi, chi 
  type(scalar):: ror,rop,oz
  integer :: i,jj,kk
  real(p8) :: val
  real(p8),allocatable,dimension(:,:):: data

  !------------------------------------------------------------------	  
  call readcom('noecho')
  call readin(5)

  !-------------------------------
  call allocate(psi)
  call allocate(chi)
  call allocate(ror)
  call allocate(rop)
  call allocate(oz)

  !-------------------------------
  allocate(data(ntchop,nxchop))

  call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//'psi0.dat',psi)
  call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//'chi0.dat',chi)

!  call pc2vor(psi,chi,ror,rop,oz)
!  call tofp(oz)
!  call checkGaussian(oz,2.0_p8,0.0_p8,0)
!  stop

  data = energy_spec(psi,chi)
  
  call PRINT_ENERGY_SPECTRUM(psi,chi,0)
  
  val = 0.0_p8
  do jj=1,ntchop
     do kk=1,nxchop
        val = val + data(jj,kk)
     end do
  end do
  print *,'Total energy =', val
  print *,'log-term     =',psi%ln

!  psi%ln = 0.0_p8
!  call tofp(psi)
!  call msave(psi%e(:,1,:), TRIM(ADJUSTL(FILES%SAVEDIR))//'psiZerodeg3')
!  call msave(psi%e(:,nth/2+1,:), TRIM(ADJUSTL(FILES%SAVEDIR))//'psiPideg3')
  call msave(tfm%r, TRIM(ADJUSTL(FILES%SAVEDIR))//'r.dat')
  call msave(tfm%th, TRIM(ADJUSTL(FILES%SAVEDIR))//'th.dat')
  
  call deallocate(psi)
  call deallocate(chi)
  
end program viewEnergy
