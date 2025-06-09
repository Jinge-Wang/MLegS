
  program energyspec

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
  implicit none

  type(scalar):: psi,chi
  real(p8), allocatable, dimension(:,:):: espec

  call readin(5)

  call allocate(psi)
  call allocate(chi)

  call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//'psi0.dat',psi)
  call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//'chi0.dat',chi)

  allocate(espec(ntchop,nxchop))

  espec = energy_spec(psi,chi)

  call msave(espec, TRIM(ADJUSTL(FILES%SAVEDIR))//'espec.dat')

  end
