program mirror
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
  type(scalar)::b
  real(p8) :: d,scale 
  integer::i,j,k

  !Distance of separation between 2 vortices, scale of translation
  ! -------------------------
  
  type(scalar):: psi,chi,ror,rop,oz
  call readin(5)
 
  !> read in initial values
  call allocate(psi)
  call allocate(chi)
  call allocate(ror)
  call allocate(rop)
  call allocate(oz)

  call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//'psi0.dat',psi)
  call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//'chi0.dat',chi)

  !Distance between 2 vortex centers
  d = 0.0_p8
  print *, '-----------------------'
  print *, 'vortex x-location:',d/2.0_p8

  call computeMirror(psi,chi,d)
  ! This is an exact way of doing it except the one aliasing term
  call pc2vor(psi,chi,ror,rop,oz)
  call tofp(oz)
  call checkGaussian(oz,d/2.0_p8,0.0_p8,0)
  
  call msave(psi, TRIM(ADJUSTL(FILES%SAVEDIR))//'psi0.dat')
  call msave(chi, TRIM(ADJUSTL(FILES%SAVEDIR))//'chi0.dat')
  
  call deallocate(psi)
  call deallocate(chi)
  
end program mirror

!------------- test 1 (passed)
!------------- should get a cosine
!  psi%e(:,:,:) = 0.0_p8
!  psi%e(1,1,2) = 0.5_p8*dsqrt(2.0_p8)
!  psi%e(1,1,nxchopdim) = 0.5_p8*dsqrt(2.0_p8)
!-------------
!------------- test 2 (passed)
!------------- should get (1-x^2)^0.5/sqrt(4/3) in (x=-0.9894,phi=0,z=0)
!  psi%e(:,:,:) = 0.0_p8
!  psi%e(1,2,1) = 1.0_p8
!  print *, '1:', 2.0_p8*dsqrt(1-tfm%x(1)**2)/dsqrt(4.0_p8/3.0_p8)
!-------------

!------------- test 3
!------------- should get 3*(1-x^2)/sqrt(9.6) in (x=-0.9894,phi=0,z=0)
!  psi%e(:,:,:) = 0.0_p8
!  psi%e(1,3,1) = 1.0_p8
!  print *, '1:', 3.0_p8*(1-tfm%x(1)**2)/dsqrt(9.6_p8)
!-------------

!  psi%ln=0.0_p8
!  chi%ln=0.0_p8
!  call tofp(psi)
!  call tofp(chi)
!  call msave(psi%e,'psiPPP')
!  print *,'2:', psi%e(1,1,1)
!  call msave(chi%e,'chiPPP')
!  call toff(psi)
