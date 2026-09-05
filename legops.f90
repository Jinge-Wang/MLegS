!----------------
 module legops
!----------------
 use eig
 use scalar3
 use bandmat
 use lin_legendre
 implicit none
 integer, parameter :: p8=selected_real_kind(p=11)
 real(p8),parameter:: pi=3.141592653589793238462643_p8
!----------------------------------------------------------
 private
 public:: scalar
 public:: delsqh, idelsqh

 public:: mulxm, mulxp

 public:: xxdx

 public:: del2, idel2, idel2ln
    !  del2(a,b) computes b=del^2(a).
    !  idel2(b,a) inverts Laplacian. boundary condition is
    !     that a at r=\infty is zero (including log-term).
    !  idel2ln(b,a) computes logarithmic term.

 public:: helm, ihelm

 ! public:: helm4, ihelm4

 public:: helmp, ihelmp, ihelmpo

 public:: pc2vel, pc2vor
    !  pc2vel: (psi,chi) to velocity (r*u_r,r*u_phi,u_z)
    !  pc2vor: (psi,chi) to vorticity (r*o_r,r*o_phi,o_z)

 public:: vprod, project, nonlin
 public:: nonlinadd,nadd
    !  compute vector product  
    !  project into toroidal and poloidal space

 public:: monitor, velmon

!----------------------------------------------------------
 type nonlinadd
   real(p8):: strain
   real(p8):: rot	
   real(p8):: u,ang,uz
 endtype
 type(nonlinadd):: nadd
 type monitor
   integer:: n,int
   real(p8):: urmax,upmax,uzmax
 endtype
 type(monitor):: velmon

!----------------------------------------------------------
 contains
!---------
 
!
!
!
      subroutine delsqh(a,b)
!---------------------------------------------------------------
!   modified horizontal delta square operator  (delsq)/(1-x)^2
!   call in polynomial/fourier space
!---------------------------------------------------------------
      type(scalar),intent(in):: a
      type(scalar),intent(inout):: b
      integer:: mm,nn,n

      if(a%space.ne. fff_space) then
        print *,'delsqh:input not in fff_space.'
	stop
      endif
      if(b%space.ne.fff_space) then
        call deallocate(b)
        call allocate(b)
      endif
      if(m(1).ne.0) then
	print *,'delsqh: m(1) neq 0. trouble.'
	stop
      endif

      b%ln = 0.0
      call chopdo(b)
      do mm=1,ntchop
      do nn=1,nzchops(mm)
	n = m(mm) + (nn-1)
        b%e(nn,mm,:) = -a%e(nn,mm,:)*(n*(n+1)/eel2)
      enddo
      enddo
      b%e(1,1,1) = a%ln/eel2/tfm%norm(1,1)

      return
      end subroutine
!
!
!
      subroutine idelsqh(b,a)
!    ----------------------------------------------------
!     inverse horizontal delsq operator ((delsq)/(1-x)^2)^(-1)
!     logarithmic term returns in a%ln.
!     call in polynomial/fourier space. 
!    ----------------------------------------------------
      type(scalar),intent(in):: b
      type(scalar),intent(inout):: a
      integer:: mm,n,nn
 
      if(b%ln.ne.0.0_p8) then
        print *,'idelsqh: input with logterm'
      endif
      if(m(1).ne.0) then
	print *,'m(1) neq 0'
	stop
      endif
      if(b%space.ne. fff_space) then
        print *,'idelsqh:input not in fff_space.'
	stop
      endif
      if(a%space.ne.fff_space) then
        call deallocate(a)
        call allocate(a)
      endif
 
      call chopdo(a)
      a%ln = b%e(1,1,1)*tfm%norm(1,1)*eel2
 
      a%e(1,1,:)=0.0
      do nn=2,nzchops(1)
        n = nn-1
        a%e(nn,1,:) = -b%e(nn,1,:)*(eel2/(n*(n+1)))
      enddo

      do mm=2,ntchop
      do nn=1,nzchops(mm)
        n = m(mm) + (nn-1)
        a%e(nn,mm,:) = -b%e(nn,mm,:)*(eel2/(n*(n+1)))
      enddo
      enddo
 
      a%e(1,1,:)=a%e(1,1,:)-calcat1(a)/tfm%at1(1)
 
      return
      end subroutine

!
!
!
      subroutine mulxp(a,b)
!---------------------------------------------------------------
!   operator  (1+x)
!   call in function/fourier space
!---------------------------------------------------------------
      type(scalar),intent(in):: a
      type(scalar),intent(inout):: b
      integer:: mm,nn,n,kk
      real(p8),dimension(:,:),allocatable:: xp

      if(a%space.ne. fff_space) then
        print *,'mulxp:input not in fff_space.'
	stop
      endif
      if(b%space.ne.fff_space) then
        call deallocate(b)
        call allocate(b)
      endif
      if(a%ln.ne.0.0) print *,'mulxp: a%ln nonzero'

      allocate( xp(nzchop,3) )

      b%ln = 0.0
      call chopdo(b)
      do mm=1,ntchop
	nn = nzchops(mm)
	xp(:nn,:) = band_leg_xp(nn,m(mm),tfm%norm(:,mm))
        do kk=1,nxchopdim
	   b%e(:nn,mm,kk)=banmul(xp(:nn,:),2,a%e(:nn,mm,kk))
        enddo
      enddo

      deallocate( xp )

      return
      end subroutine

!
!
!
      subroutine mulxm(a,b)
!---------------------------------------------------------------
!   operator  (1-x)
!   call in function/fourier space
!---------------------------------------------------------------
      type(scalar),intent(in):: a
      type(scalar),intent(inout):: b
      integer:: mm,nn,n,kk
      real(p8),dimension(:,:),allocatable:: xm

      if(a%space.ne. fff_space) then
        print *,'mulxm:input not in fff_space.'
	stop
      endif
      if(b%space.ne.fff_space) then
        call deallocate(b)
        call allocate(b)
      endif
      if(a%ln.ne.0.0)         print *,'mulxm: a%ln nonzero'

      allocate( xm(nzchop,3) )

      b%ln = 0.0
      call chopdo(b)
      do mm=1,ntchop
	nn = nzchops(mm)
	xm(:nn,:) = band_leg_xm(nn,m(mm),tfm%norm(:,mm))
        do kk=1,nxchopdim
	   b%e(:nn,mm,kk)=banmul(xm(:nn,:),2,a%e(:nn,mm,kk))
        enddo
      enddo

      deallocate( xm )

      return
      end subroutine
!
!
!
      subroutine xxdx(a,b)
!---------------------------------------------------------------
!   operator  (1-x^2)d/dx
!   call in function/fourier space
!---------------------------------------------------------------
      type(scalar),intent(in):: a
      type(scalar),intent(inout):: b
      integer:: mm,nn,n,kk
      real(p8),dimension(:,:),allocatable:: xxdxop

      if(a%space.ne. fff_space) then
        print *,'xxdx:input not in fff_space.'
	stop
      endif
      if(b%space.ne.fff_space) then
        call deallocate(b)
        call allocate(b)
      endif

      if(a%e(nzchops(2),2,1).ne.0.0) then
	!print *,'xxdx: operation not exact'
        !print *, a%e(nzchops(2),2,1)
        !stop
      endif

      allocate( xxdxop(nzchop,3) )

      b%ln = 0.0
      call chopdo(b)
      do mm=1,ntchop
	nn = nzchops(mm)
	xxdxop(:nn,:) = band_leg_xxdx(nn,m(mm),tfm%norm(:,mm))
        do kk=1,nxchopdim
	   b%e(:nn,mm,kk)=banmul(xxdxop(:nn,:),2,a%e(:nn,mm,kk))
        enddo
      enddo

      b%e(1,1,1)=b%e(1,1,1) +a%ln/tfm%norm(1,1)
      b%e(2,1,1)=b%e(2,1,1) +a%ln/tfm%norm(2,1)

      deallocate( xxdxop )

      return
      end subroutine
!
!
!
      subroutine del2(a,b)
!---------------------------------------------------------------
!   laplacian operator
!   call in function/fourier space
!---------------------------------------------------------------
      type(scalar),intent(in):: a
      type(scalar),intent(inout):: b
      integer:: mm,nn,n,kk
      real(p8),dimension(:,:),allocatable:: del2op

      if(a%space.ne. fff_space) then
        print *,'del2:input not in fff_space.'
	stop
      endif
      if(b%space.ne.fff_space) then
        call deallocate(b)
        call allocate(b)
      endif
      if(a%e(nzchops(2),2,1).ne.0.0) then
	print *,'del2: operation not exact'
      endif

      allocate( del2op(nzchop,5) )

      b%ln = 0.0
      call chopdo(b)
      do mm=1,ntchop
	nn = nzchops(mm)
	del2op(:nn,:) = band_leg_rat_del2h(nn,m(mm),eel,tfm%norm(:,mm))
        do kk=1,nxchopdim
	   b%e(:nn,mm,kk)=banmul(del2op(:nn,:),3,a%e(:nn,mm,kk))
	   b%e(:nn,mm,kk)=b%e(:nn,mm,kk)-(ak(mm,kk)**2)*a%e(:nn,mm,kk)
        enddo
      enddo

      b%e(1,1,1)=b%e(1,1,1) +4.0_p8/3*a%ln/eel2/tfm%norm(1,1)
      b%e(2,1,1)=b%e(2,1,1) -2.0_p8*a%ln/eel2/tfm%norm(2,1)
      b%e(3,1,1)=b%e(3,1,1) +2.0_p8/3*a%ln/eel2/tfm%norm(3,1)

      deallocate( del2op )

      return
      end subroutine
!
!
!
      subroutine idel2(b,a,ln)
! ----------------------------------------------------
!     inverse laplacian operator
!     call in polynomial/fourier/fourier space. 
!     logterm=ln is assumed as boundary condition (default ln=0).
! ----------------------------------------------------
      type(scalar):: a, b
      real(p8),dimension(:,:,:),allocatable :: del2op
      real(p8),optional:: ln
      integer::kk,nn,mm,kv,np,nm

      if(b%space.ne. fff_space) then
        print *,'idel2:input not in fff_space.'
	stop
      endif
      if(a%space.ne.fff_space) then
        call deallocate(a)
        call allocate(a)
      endif
      if(b%ln.ne.0.0_p8) then
 	print *,'idel2:logterm for input of idel2 nonzero.'
	print *,'logterm=',b%ln
      endif
      a = b

      !> m=0
      nn=nzchops(1)
      np=nn+1
      nm=nn-1
      allocate( del2op(np,7,nxchopdim) )
      del2op(:np,2:6,1)  = band_leg_rat_del2h(np,m(1),eel,tfm%norm(:,1))
      do kk=2,nxchopdim
	del2op(:np,2:3,kk)=del2op(:np,2:3,1)
	del2op(:np,4  ,kk)=del2op(:np,4  ,1)-ak(1,kk)**2
	del2op(:np,5:6,kk)=del2op(:np,5:6,1)
	del2op(:np,1  ,kk)=0
	del2op(:np,7  ,kk)=0
      enddo
      del2op(1,4,1)=del2op(1,4,1)+4.0_p8/3/eel2
      del2op(2,3,1)=del2op(2,3,1)-2*tfm%norm(1,1)/tfm%norm(2,1)/eel2
      del2op(3,2,1)=del2op(3,2,1)+2.0_p8/3*tfm%norm(1,1)/tfm%norm(3,1)/eel2
      del2op(2:np,1:5,1)=del2op(1:np-1,2:6,1)
      del2op(:,6:7,1)=0
      del2op(1,:,1)=0
      del2op(1,4,1)=1
      a%e(2:nn,1,1)=a%e(1:nm,1,1)
      if(present(ln)) then
        a%e(1,1,1)=ln/tfm%norm(1,1)
      else
        a%e(1,1,1)=0
      endif

      call lub(del2op(:nn,:,:),4)
      call solveb(del2op(:nn,:,:),4,a%e(:nn,1,:nxchopdim))
      deallocate( del2op )

      !> m neq 0
      allocate( del2op(nzchop,5,nxchopdim) )
      do mm=2,ntchop
        nn=nzchops(mm)
	del2op(:nn,1:5,1) = band_leg_rat_del2h(nn,m(mm),eel,tfm%norm(:,mm))
        do kk=1,nxchopdim
	   del2op(:nn,1:2,kk)=del2op(:nn,1:2,1)
	   del2op(:nn,3  ,kk)=del2op(:nn,3  ,1)-ak(mm,kk)**2
	   del2op(:nn,4:5,kk)=del2op(:nn,4:5,1)
	enddo
        call lub(del2op(:nn,1:5,1:nxchopdim),3)
        call solveb(del2op(:nn,1:5,1:nxchopdim),3,a%e(:nn,mm,1:nxchopdim))
      enddo
      deallocate( del2op )

      if(present(ln)) then
        a%ln = a%e(1,1,1)*tfm%norm(1,1)
      else
        a%ln=0.0
      endif
      call chopdo(a)
      do kk=1,nxchopdim
        a%e(1,1,kk)=a%e(1,1,kk)  &
	    -sum(tfm%at1(:nzchop)*a%e(:nzchop,1,kk))/tfm%at1(1)
      enddo

      return
      end subroutine
!
!
!
      subroutine idel2ln(b,a)
! ----------------------------------------------------
!     inverse laplacian operator
!     call in polynomial/fourier/fourier space. 
!     logterm is computed.
! ----------------------------------------------------
      type(scalar):: a, b
      real(p8),dimension(:,:,:),allocatable :: del2op
      integer::kk,nn,mm,kv,np,nm

      if(b%space.ne. fff_space) then
        print *,'idel2:input not in fff_space.'
	stop
      endif
      if(a%space.ne.fff_space) then
        call deallocate(a)
        call allocate(a)
      endif
      if(b%ln.ne.0.0_p8) then
 	print *,'idel2:logterm nonzero.'
	print *,'logterm=',b%ln
      endif
      a = b

      !> m=0
      nn=nzchops(1)
      np=nn+1
      nm=nn-1
      allocate( del2op(np,5,nxchopdim) )
      del2op(:np,1:5,1)  = band_leg_rat_del2h(np,m(1),eel,tfm%norm(:,1))
      do kk=2,nxchopdim
	del2op(:np,1:2,kk)=del2op(:np,1:2,1)
	del2op(:np,3  ,kk)=del2op(:np,3  ,1)-ak(1,kk)**2
	del2op(:np,4:5,kk)=del2op(:np,4:5,1)
      enddo
      del2op(1,3,1)=del2op(1,3,1)+4.0_p8/3/eel2
      del2op(2,2,1)=del2op(2,2,1)-2*tfm%norm(1,1)/tfm%norm(2,1)/eel2
      del2op(3,1,1)=del2op(3,1,1)+2.0_p8/3*tfm%norm(1,1)/tfm%norm(3,1)/eel2

      call lub(del2op(:nn,:,:),3)
      call solveb(del2op(:nn,:,:),3,a%e(:nn,1,:nxchopdim))
      a%ln = a%e(1,1,1)*tfm%norm(1,1)
      deallocate( del2op )

      !> m neq 0
      allocate( del2op(nzchop,5,nxchopdim) )
      do mm=2,ntchop
        nn=nzchops(mm)
	del2op(:nn,1:5,1) = band_leg_rat_del2h(nn,m(mm),eel,tfm%norm(:,mm))
        do kk=1,nxchopdim
	   del2op(:nn,1:2,kk)=del2op(:nn,1:2,1)
	   del2op(:nn,3  ,kk)=del2op(:nn,3  ,1)-ak(mm,kk)**2
	   del2op(:nn,4:5,kk)=del2op(:nn,4:5,1)
	enddo
        call lub(del2op(:nn,1:5,1:nxchopdim),3)
        call solveb(del2op(:nn,1:5,1:nxchopdim),3,a%e(:nn,mm,1:nxchopdim))
      enddo
      deallocate( del2op )

      call chopdo(a)
      do kk=1,nxchopdim
        a%e(1,1,kk)=a%e(1,1,kk) &
		  -sum(tfm%at1(:nzchop)*a%e(:nzchop,1,kk))/tfm%at1(1)
      enddo

      return
      end subroutine
!
!
!
    subroutine helm(a,b,alp)
!----------------------------------------------------
!   helmholtz operator  b=(delsq-alp)a
!   call in polynomial/fourier space
!----------------------------------------------------
    type(scalar):: a, b
    real(p8):: alp
    call del2(a,b)
    b%e=b%e-alp*a%e
    b%ln = -alp*a%ln
    return
    end subroutine
!
!
!
      subroutine ihelm(b,a,alp)
! ----------------------------------------------------
!     inverse helmholtz operator inv(del^2-alp)
!     call in polynomial/fourier/fourier space. 
! ----------------------------------------------------
      type(scalar):: a, b
      real(p8),intent(in):: alp
      real(p8),dimension(:,:,:),allocatable :: helmop
      real(p8),dimension(:,:),allocatable :: xm
      integer::kk,nn,mm,kv,np,nm

      if(b%space.ne. fff_space) then
        print *,'idel2:input not in fff_space.'
	stop
      endif
      if(a%space.ne.fff_space) then
        call deallocate(a)
        call allocate(a)
      endif
      if(alp.eq.0.0) then
         print *,'ihelm: alp can''t be zero.'
         stop
      endif

      a=b
      a%ln= -b%ln/alp
      a%e(1,1,1)=a%e(1,1,1)-a%ln*(4.0_p8/3/eel2)/tfm%norm(1,1)
      a%e(2,1,1)=a%e(2,1,1)-a%ln*(-2/tfm%norm(2,1)/eel2)
      a%e(3,1,1)=a%e(3,1,1)-a%ln*(+2.0_p8/3/tfm%norm(3,1)/eel2)

      !> m=0
      nn=nzchops(1)
      np=nn+1
      nm=nn-1
      allocate( xm(np,3) )
      xm = band_leg_xm(np,m(1),tfm%norm(:,1))
      allocate( helmop(np,7,nxchopdim) )
      helmop(:np,2:6,1) = band_leg_rat_del2h(np,m(1),eel,tfm%norm(:,1))
      helmop(:np, 4 ,1) = helmop(:np, 4 ,1)-alp
      do kk=2,nxchopdim
	helmop(:np,2:3,kk)=helmop(:np,2:3,1)
	helmop(:np,4  ,kk)=helmop(:np,4  ,1)-ak(1,kk)**2
	helmop(:np,5:6,kk)=helmop(:np,5:6,1)
	helmop(:,:,kk)=banmul(helmop(:,2:6,kk),3,xm,2)
      enddo
      helmop(:,:,1)=banmul(helmop(:,2:6,1),3,xm,2)

      call lub(helmop(:nm,:,:),4)
      call solveb(helmop(:nm,:,:),4,a%e(:nm,1,:nxchopdim))
      a%e(nn,1,:nxchopdim)=0
      do kk=1,nxchopdim
	a%e(:nn,1,kk)=banmul(xm(:nn,:),2,a%e(1:nn,1,kk))
      enddo
      deallocate( xm )
      deallocate( helmop )

      !> m neq 0
      allocate( helmop(nzchop,5,nxchopdim) )
      do mm=2,ntchop
        nn=nzchops(mm)
	helmop(:nn,1:5,1) = band_leg_rat_del2h(nn,m(mm),eel,tfm%norm(:,mm))
	helmop(:nn, 3 ,1) = helmop(:nn, 3 ,1)-alp
        do kk=1,nxchopdim
	   helmop(:nn,1:2,kk)=helmop(:nn,1:2,1)
	   helmop(:nn,3  ,kk)=helmop(:nn,3  ,1)-ak(mm,kk)**2
	   helmop(:nn,4:5,kk)=helmop(:nn,4:5,1)
	enddo
        call lub(helmop(:nn,1:5,1:nxchopdim),3)
        call solveb(helmop(:nn,1:5,1:nxchopdim),3,a%e(:nn,mm,1:nxchopdim))
      enddo

      deallocate( helmop )

      call chopdo(a)

      return
      end subroutine
!
!
!
    subroutine helmp(p,a,b,alp,nu)
! ----------------------------------------------------
!   del^p helmholtz operator  b=(del^p +nu*del^2 -alp)a
!   p even, p>=4.
!   call in polynomial/fourier space
! ----------------------------------------------------
    type(scalar):: a,b
    real(p8):: alp
    real(p8),optional:: nu
    type(scalar):: w,d2
    integer:: i,p

    call allocate(w)
    call allocate(d2)
    call del2(a,d2)
    b=d2
    do i=4,p,2
      call del2(b,w)
      b=w
    enddo

    b%e=b%e -alp*a%e
    if(present(nu)) then
      b%e=b%e +nu*d2%e
    endif
    b%ln=-alp*a%ln

    call deallocate(w)
    call deallocate(d2)

    return
    end subroutine
!
!
!
      subroutine ihelmp(p,b,a,alp,nuin)
! ----------------------------------------------------
!     inverse (del^p +nu*del^2 -alp) operator
!     p is even and p>=4.
!     call in polynomial/fourier/fourier space. 
! ----------------------------------------------------
      type(scalar):: a, b
      real(p8),intent(in):: alp
      real(p8),intent(in),optional:: nuin
      real(p8),dimension(:,:,:),allocatable :: helmp
      real(p8),dimension(:,:),allocatable :: del2h,del2
      real(p8),dimension(:,:),allocatable :: del4h,del6h
      real(p8),dimension(:,:),allocatable :: del8h,del10h
      real(p8),dimension(:,:),allocatable :: del12h,del14h
      real(p8),dimension(nzchop):: al2,al
      integer,intent(in):: p
      integer::kk,nn,mm,kv,np,nm,i
      real(p8):: nu,eel4,k2

      if(b%space.ne. fff_space) then
        print *,'idel2:input not in fff_space.'
	stop
      endif
      if(a%space.ne.fff_space) then
        call deallocate(a)
        call allocate(a)
      endif
      if(m(1).ne.0) then
         print *,'ihelmp: m(1) must be zero.'
         stop
      endif
      if(alp.eq.0.0) then
         print *,'ihelmp: alp can''t be zero.'
         stop
      endif
      if(present(nuin)) then
	nu=nuin
      else
	nu=0.0_p8
      endif
      a=b

      !> subtract logarithmic term contribution
      a%ln= -b%ln/alp
      eel4=eel2**2
      al2=0
      al2(1)=a%ln*(4.0_p8/3/eel2)/tfm%norm(1,1)
      al2(2)=a%ln*(-2/tfm%norm(2,1)/eel2)
      al2(3)=a%ln*(+2.0_p8/3/tfm%norm(3,1)/eel2)
      al=al2
      allocate( del2h(nzchop,5) )
      del2h = band_leg_rat_del2h(nzchop,m(1),eel,tfm%norm(:,1))
      do i=4,p,2
	al=banmul(del2h,3,al)
      enddo
      a%e(:nzchop,1,1)=a%e(:nzchop,1,1)-al-nu*al2
      deallocate( del2h )

      !> all m
      allocate( del2(nzchop+p,5) )
      allocate( helmp(nzchop+p,2*p+1,nxchopdim) )
      allocate( del2h(nzchop+p,5) )
      allocate( del4h(nzchop+p,9) )
      if(p >=6 ) allocate( del6h(nzchop+p,13) )
      if(p >=8 ) allocate( del8h(nzchop+p,17) )
      if(p >=10) allocate( del10h(nzchop+p,21) )
      if(p >=12) allocate( del12h(nzchop+p,25) )
      if(p >=14) allocate( del14h(nzchop+p,29) )
      do mm=1,ntchop
        nn=nzchops(mm)
	np=nn+p
	del2h(:np,:) = band_leg_rat_del2h(np,m(mm),eel,tfm%norm(:,mm))
	del4h(:np,:) = banmul(del2h(:np,:),3,del2h(:np,:),3)
        if(p >=6 ) del6h(:np,:)=banmul(del4h(:np,:),5,del2h(:np,:),3)
        if(p >=8 ) del8h(:np,:)=banmul(del6h(:np,:),7,del2h(:np,:),3)
        if(p >=10) del10h(:np,:)=banmul(del8h(:np,:),9,del2h(:np,:),3)
        if(p >=12) del12h(:np,:)=banmul(del10h(:np,:),11,del2h(:np,:),3)
        if(p >=14) del14h(:np,:)=banmul(del12h(:np,:),13,del2h(:np,:),3)
        do kk=1,nxchopdim
	   k2 = ak(mm,kk)**2
	   del2(:np,:)=del2h(:np,:)
	   del2(:np,3)=del2(:np,3)-k2
	   if(p.eq.4) then
	     helmp(:np, : ,kk)=                          del4h(:np,:)
	     helmp(:np,3:7,kk)=helmp(:np,3:7,kk)-2*k2   *del2h(:np,:)
	     helmp(:np,  5,kk)=helmp(:np,5,kk)  +  k2**2
	   else if(p.eq.6) then
	     helmp(:np, :  ,kk)=                           del6h(:np,:)
	     helmp(:np,3:11,kk)=helmp(:np,3:11,kk)-3*k2   *del4h(:np,:)
	     helmp(:np,5: 9,kk)=helmp(:np,5: 9,kk)+3*k2**2*del2h(:np,:)
	     helmp(:np,   7,kk)=helmp(:np,   7,kk)-  k2**3
	   else if(p.eq.8) then
	     helmp(:np, : ,kk)   =                         del8h(:np,:)
	     helmp(:np,3:15,kk)=helmp(:np,3:15,kk)-4*k2   *del6h(:np,:)
	     helmp(:np,5:13,kk)=helmp(:np,5:13,kk)+6*k2**2*del4h(:np,:)
	     helmp(:np,7:11,kk)=helmp(:np,7:11,kk)-4*k2**3*del2h(:np,:)
	     helmp(:np,   9,kk)=helmp(:np,   9,kk)+  k2**4
	   else if(p.eq.10) then
	     helmp(:np,:,kk)   =                            del10h(:np,:)
	     helmp(:np,3:19,kk)=helmp(:np,3:19,kk)- 5*k2   *del8h(:np,:)
	     helmp(:np,5:17,kk)=helmp(:np,5:17,kk)+10*k2**2*del6h(:np,:)
	     helmp(:np,7:15,kk)=helmp(:np,7:15,kk)-10*k2**3*del4h(:np,:)
	     helmp(:np,9:13,kk)=helmp(:np,9:13,kk)+ 5*k2**4*del2h(:np,:)
	     helmp(:np,  11,kk)=helmp(:np,  11,kk)-   k2**5
	   else if(p.eq.12) then
	     helmp(:np,:,kk)   =                            del12h(:np,:)
	     helmp(:np, 3:23,kk)=helmp(:np, 3:23,kk)- 6*k2   *del10h(:np,:)
	     helmp(:np, 5:21,kk)=helmp(:np, 5:21,kk)+15*k2**2*del8h(:np,:)
	     helmp(:np, 7:19,kk)=helmp(:np, 7:19,kk)-20*k2**3*del6h(:np,:)
	     helmp(:np, 9:17,kk)=helmp(:np, 9:17,kk)+15*k2**4*del4h(:np,:)
	     helmp(:np,11:15,kk)=helmp(:np,11:15,kk)- 6*k2**5*del2h(:np,:)
	     helmp(:np,   13,kk)=helmp(:np,  13,kk)+    k2**6
	   else if(p.eq.14) then
	     helmp(:np,:,kk)   =                              del14h(:np,:)
	     helmp(:np, 3:27,kk)=helmp(:np, 3:27,kk)- 7*k2   *del12h(:np,:)
	     helmp(:np, 5:25,kk)=helmp(:np, 5:25,kk)+21*k2**2*del10h(:np,:)
	     helmp(:np, 7:23,kk)=helmp(:np, 7:23,kk)-35*k2**3*del8h(:np,:)
	     helmp(:np, 9:21,kk)=helmp(:np, 9:21,kk)+35*k2**4*del6h(:np,:)
	     helmp(:np,11:19,kk)=helmp(:np,11:19,kk)-21*k2**5*del4h(:np,:)
	     helmp(:np,13:17,kk)=helmp(:np,13:17,kk)+ 7*k2**6*del2h(:np,:)
	     helmp(:np,   15,kk)=helmp(:np,   15,kk)-   k2**7
	   else 
	     print *,'p out of range. p=',p
	     stop
	   endif
	   helmp(:nn,p-1:p+3,kk)=helmp(:nn,p-1:p+3,kk)+nu*del2(:nn,:)
	   helmp(:nn,p+1,kk)=helmp(:nn,p+1,kk)-alp
	enddo
        call lub(helmp(:nn,:,:),p+1)
        call solveb(helmp(:nn,:,:),p+1,a%e(:nn,mm,1:nxchopdim))
      enddo

      deallocate( del2h, del4h, del2, helmp )
      if(p >=6 ) deallocate( del6h )
      if(p >=8 ) deallocate( del8h )
      if(p >=10) deallocate( del10h )
      if(p >=12) deallocate( del12h )
      if(p >=14) deallocate( del14h )
      !> zero at infinity
      a%e(1,1,:)=a%e(1,1,:)-calcat1(a)/tfm%pf(1,1,1)

      call chopdo(a)

      return
      end subroutine
!
!
!
      subroutine ihelmpo(p,b,a,alp,nuin)
! ----------------------------------------------------
!     inverse (del^p +nu*del^2 -alp) operator
!     p is even and p>=4.
!     call in polynomial/fourier/fourier space. 
! ----------------------------------------------------
      type(scalar):: a, b
      real(p8),intent(in):: alp
      real(p8),intent(in),optional:: nuin
      real(p8),dimension(:,:,:),allocatable :: helmp
      real(p8),dimension(:,:),allocatable :: del2h,del2
      real(p8),dimension(nzchop):: al2,al
      integer,intent(in):: p
      integer::kk,nn,mm,kv,np,nm,i
      real(p8):: nu,eel4

      if(b%space.ne. fff_space) then
        print *,'idel2:input not in fff_space.'
	stop
      endif
      if(a%space.ne.fff_space) then
        call deallocate(a)
        call allocate(a)
      endif
      if(m(1).ne.0) then
         print *,'ihelmp: m(1) must be zero.'
         stop
      endif
      if(alp.eq.0.0) then
         print *,'ihelmp: alp can''t be zero.'
         stop
      endif
      if(present(nuin)) then
	nu=nuin
      else
	nu=0.0_p8
      endif
      a=b

      !> subtract logarithmic term contribution
      a%ln= -b%ln/alp
      eel4=eel2**2
      al2=0
      al2(1)=a%ln*(4.0_p8/3/eel2)/tfm%norm(1,1)
      al2(2)=a%ln*(-2/tfm%norm(2,1)/eel2)
      al2(3)=a%ln*(+2.0_p8/3/tfm%norm(3,1)/eel2)
      al=al2
      allocate( del2h(nzchop,5) )
      del2h = band_leg_rat_del2h(nzchop,m(1),eel,tfm%norm(:,1))
      do i=4,p,2
	al=banmul(del2h,3,al)
      enddo
      a%e(:nzchop,1,1)=a%e(:nzchop,1,1)-al-nu*al2
      deallocate( del2h )

      !> all m
      allocate( del2h(nzchop+p,5) )
      allocate( del2(nzchop+p,5) )
      allocate( helmp(nzchop+p,2*p+1,nxchopdim) )
      do mm=1,ntchop
        nn=nzchops(mm)
	np=nn+p
	del2h(:np,:) = band_leg_rat_del2h(np,m(mm),eel,tfm%norm(:,mm))
        do kk=1,nxchopdim
	   del2(:np,:)=del2h(:np,:)
	   del2(:np,3)=del2(:np,3)-ak(mm,kk)**2
	   helmp(:np,p-1:p+3,kk)=del2(:np,:)
	   do i=4,p,2
	     helmp(:np,p+1-i:p+1+i,kk)= &
	       banmul(helmp(:np,p+3-i:p-1+i,kk),i-1, del2(:np,:),3 )
	   enddo
	   helmp(:nn,p-1:p+3,kk)=helmp(:nn,p-1:p+3,kk)+nu*del2(:nn,:)
	   helmp(:nn,p+1,kk)=helmp(:nn,p+1,kk)-alp
	enddo
        call lub(helmp(:nn,:,:),p+1)
        call solveb(helmp(:nn,:,:),p+1,a%e(:nn,mm,1:nxchopdim))
      enddo

      deallocate( del2h, del2, helmp )
      !> zero at infinity
      a%e(1,1,:)=a%e(1,1,:)-calcat1(a)/tfm%pf(1,1,1)

      call chopdo(a)

      return
      end subroutine
!! !
!! !
!! !
!!     subroutine helm4(a,b,alp,nu)
!! ! ----------------------------------------------------
!! !   del4 helmholtz operator  b=(del^4 +nu*del^2 -alp)a
!! !   call in polynomial/fourier space
!! ! ----------------------------------------------------
!!     type(scalar):: a,b
!!     real(p8):: alp
!!     real(p8),optional:: nu
!!     type(scalar):: w
!! 
!!     call allocate(w)
!!     call del2(a,w)
!!     call del2(w,b)
!! 
!!     b%e=b%e -alp*a%e
!!     if(present(nu)) then
!!       b%e=b%e +nu*w%e
!!     endif
!!     b%ln=-alp*a%ln
!! 
!!     call deallocate(w)
!! 
!!     return
!!     end subroutine
!! !
!! !
!! !
!!       subroutine ihelm4(b,a,alp,nuin)
!! ! ----------------------------------------------------
!! !     inverse (del^4 +nu*del^2 -alp) operator
!! !     call in polynomial/fourier/fourier space. 
!! ! ----------------------------------------------------
!!       type(scalar):: a, b
!!       real(p8),intent(in):: alp
!!       real(p8),intent(in),optional:: nuin
!!       real(p8),dimension(:,:,:),allocatable :: helm4
!!       real(p8),dimension(:,:),allocatable :: del2h,del4h,xm
!!       integer::kk,nn,mm,kv,np,nm
!!       real(p8):: nu,eel4
!! 
!!       if(b%space.ne. fff_space) then
!!         print *,'idel2:input not in fff_space.'
!! 	stop
!!       endif
!!       if(a%space.ne.fff_space) then
!!         call deallocate(a)
!!         call allocate(a)
!!       endif
!!       if(m(1).ne.0) then
!!          print *,'ihelm4: m(1) must be zero.'
!!          stop
!!       endif
!!       if(alp.eq.0.0) then
!!          print *,'ihelm4: alp can''t be zero.'
!!          stop
!!       endif
!!       if(present(nuin)) then
!! 	nu=nuin
!!       else
!! 	nu=0.0_p8
!!       endif
!!       a=b
!! 
!!       !> subtract logarithmic term contribution
!!       a%ln= -b%ln/alp
!!       eel4=eel2**2
!!       a%e(1,1,1)=a%e(1,1,1)-a%ln*(-16.0_p8/5/eel4)/tfm%norm(1,1)
!!       a%e(2,1,1)=a%e(2,1,1)-a%ln*( 48.0_p8/5/eel4)/tfm%norm(2,1)
!!       a%e(3,1,1)=a%e(3,1,1)-a%ln*(-80.0_p8/7/eel4)/tfm%norm(3,1)
!!       a%e(4,1,1)=a%e(4,1,1)-a%ln*( 32.0_p8/5/eel4)/tfm%norm(4,1)
!!       a%e(5,1,1)=a%e(5,1,1)-a%ln*(-48.0_p8/35/eel4)/tfm%norm(5,1)
!!       !
!!       a%e(1,1,1)=a%e(1,1,1)-nu*a%ln*(4.0_p8/3/eel2)/tfm%norm(1,1)
!!       a%e(2,1,1)=a%e(2,1,1)-nu*a%ln*(-2/tfm%norm(2,1)/eel2)
!!       a%e(3,1,1)=a%e(3,1,1)-nu*a%ln*(+2.0_p8/3/tfm%norm(3,1)/eel2)
!! 
!!       !> m=0
!!       nn=nzchops(1)
!!       np=nn+2
!!       nm=nn-1
!!       allocate( helm4(nn,11,nxchopdim) )
!!       allocate( xm(np,3) )
!!       xm = band_leg_xm(np,m(1),tfm%norm(:,1))
!!       allocate( del2h(np,5) )
!!       allocate( del4h(np,9) )
!!       del2h = band_leg_rat_del2h(np,m(1),eel,tfm%norm(:,1))
!!       del4h = banmul(del2h(:np,:),3,del2h(:np,:),3)
!!       do kk=1,nxchopdim
!! 	helm4(:nn,2:10,kk)=del4h(:nn,:)
!! 	helm4(:nn,4:8,kk)=helm4(:nn,4:8,kk)  &
!! 		      -(2*ak(1,kk)**2-nu)*del2h(:nn,:)
!! 	helm4(:nn,6  ,kk)=helm4(:nn,6  ,kk) &
!! 		      +(ak(1,kk)**4-nu*ak(1,kk)**2-alp)
!! 	helm4(:nn,:,kk)=banmul(helm4(:nn,2:10,kk),5,xm(:nn,:),2)
!!       enddo
!!       call lub(helm4(:nm,:,:),6)
!!       call solveb(helm4(:nm,:,:),6,a%e(:nm,1,:nxchopdim))
!!       a%e(nn,1,:nxchopdim)=0
!!       do kk=1,nxchopdim
!! 	a%e(1:nn,1,kk)=banmul(xm(:nn,:),2,a%e(1:nn,1,kk))
!!       enddo
!!       deallocate( xm )
!!       deallocate( del2h, del4h, helm4 )
!! 
!!       !> m neq 0
!!       allocate( del2h(nzchop+2,5) )
!!       allocate( del4h(nzchop+2,9) )
!!       allocate( helm4(nzchop,9,nxchopdim) )
!!       do mm=2,ntchop
!!         nn=nzchops(mm)
!! 	np=nn+2
!! 	del2h(:np,:) = band_leg_rat_del2h(np,m(mm),eel,tfm%norm(:,mm))
!! 	del4h(:np,:) = banmul(del2h(:np,:),3,del2h(:np,:),3)
!!         do kk=1,nxchopdim
!! 	   helm4(:nn,1:9,kk)=del4h(:nn,:)
!! 	   helm4(:nn,3:7,kk)=helm4(:nn,3:7,kk)  &
!! 			      -(2*ak(mm,kk)**2-nu)*del2h(:nn,:)
!! 	   helm4(:nn,5  ,kk)=helm4(:nn,5  ,kk) &
!! 			      +(ak(mm,kk)**4-nu*ak(mm,kk)**2-alp)
!! 	enddo
!!         call lub(helm4(:nn,:,:),5)
!!         call solveb(helm4(:nn,:,:),5,a%e(:nn,mm,1:nxchopdim))
!!       enddo
!! 
!!       deallocate( del2h, del4h, helm4 )
!! 
!!       call chopdo(a)
!! 
!!       return
!!       end subroutine
!! !
!! !
!! !
!!       subroutine ihelm4(b,a,alp,nuin)
!! ! ----------------------------------------------------
!! !     inverse (del^4 +nu*del^2 -alp) operator
!! !     call in polynomial/fourier/fourier space. 
!! ! ----------------------------------------------------
!!       type(scalar):: a, b
!!       real(p8),intent(in):: alp
!!       real(p8),intent(in),optional:: nuin
!!       real(p8),dimension(:,:,:),allocatable :: helm4
!!       real(p8),dimension(:,:),allocatable :: del2h,del4h
!!       integer::kk,nn,mm,kv,np,nm
!!       real(p8):: nu,eel4
!! 
!!       if(b%space.ne. fff_space) then
!!         print *,'idel2:input not in fff_space.'
!! 	stop
!!       endif
!!       if(a%space.ne.fff_space) then
!!         call deallocate(a)
!!         call allocate(a)
!!       endif
!!       if(m(1).ne.0) then
!!          print *,'ihelm4: m(1) must be zero.'
!!          stop
!!       endif
!!       if(alp.eq.0.0) then
!!          print *,'ihelm4: alp can''t be zero.'
!!          stop
!!       endif
!!       if(present(nuin)) then
!! 	nu=nuin
!!       else
!! 	nu=0.0_p8
!!       endif
!!       a=b
!! 
!!       !> subtract logarithmic term contribution
!!       a%ln= -b%ln/alp
!!       eel4=eel2**2
!!       a%e(1,1,1)=a%e(1,1,1)-a%ln*(-16.0_p8/5/eel4)/tfm%norm(1,1)
!!       a%e(2,1,1)=a%e(2,1,1)-a%ln*( 48.0_p8/5/eel4)/tfm%norm(2,1)
!!       a%e(3,1,1)=a%e(3,1,1)-a%ln*(-80.0_p8/7/eel4)/tfm%norm(3,1)
!!       a%e(4,1,1)=a%e(4,1,1)-a%ln*( 32.0_p8/5/eel4)/tfm%norm(4,1)
!!       a%e(5,1,1)=a%e(5,1,1)-a%ln*(-48.0_p8/35/eel4)/tfm%norm(5,1)
!!       !
!!       a%e(1,1,1)=a%e(1,1,1)-nu*a%ln*(4.0_p8/3/eel2)/tfm%norm(1,1)
!!       a%e(2,1,1)=a%e(2,1,1)-nu*a%ln*(-2/tfm%norm(2,1)/eel2)
!!       a%e(3,1,1)=a%e(3,1,1)-nu*a%ln*(+2.0_p8/3/tfm%norm(3,1)/eel2)
!! 
!!       !> all m
!!       allocate( del2h(nzchop+2,5) )
!!       allocate( del4h(nzchop+2,9) )
!!       allocate( helm4(nzchop,9,nxchopdim) )
!!       do mm=1,ntchop
!!         nn=nzchops(mm)
!! 	np=nn+2
!! 	del2h(:np,:) = band_leg_rat_del2h(np,m(mm),eel,tfm%norm(:,mm))
!! 	del4h(:np,:) = banmul(del2h(:np,:),3,del2h(:np,:),3)
!!         do kk=1,nxchopdim
!! 	   helm4(:nn,1:9,kk)=del4h(:nn,:)
!! 	   helm4(:nn,3:7,kk)=helm4(:nn,3:7,kk)  &
!! 			      -(2*ak(mm,kk)**2-nu)*del2h(:nn,:)
!! 	   helm4(:nn,5  ,kk)=helm4(:nn,5  ,kk) &
!! 			      +(ak(mm,kk)**4-nu*ak(mm,kk)**2-alp)
!! 	enddo
!!         call lub(helm4(:nn,:,:),5)
!!         call solveb(helm4(:nn,:,:),5,a%e(:nn,mm,1:nxchopdim))
!!       enddo
!! 
!!       deallocate( del2h, del4h, helm4 )
!! 
!!       a%e(1,1,:)=a%e(1,1,:)-calcat1(a)/tfm%pf(1,1,1)
!! 
!!       call chopdo(a)
!! 
!!       return
!!       end subroutine
!
!
!
    subroutine pc2vel(psi,chi,rur,rup,uz,c)
!   psi,chi to velocity
!   returns (r*ur,r*up,uz) if c is not present
!   returns (r*ur,r*up,uz/(1-mu)^2 if c is present
!---------------------------------------------
    type(scalar):: psi,chi,rur,rup,uz
    integer:: mm,kk,nn,mv
    complex(p8),parameter:: iu=(0.0_p8,1.0_p8)
    real(p8):: kv
    type(scalar):: w
    character(len=1),optional:: c

    call xxdx(chi,rur)
    !rur now has  (1-x^2)d/dx chi
    call xxdx(psi,rup)
    !rup now has  (1-x^2)d/dx psi
    do mm=1,ntchop
       nn=nzchops(mm)
       mv=m(mm)
       do kk=1,nxchopdim
          kv=ak(mm,kk)
          rur%e(:nn,mm,kk)=iu*mv*psi%e(:nn,mm,kk)+(iu*kv)*rur%e(:nn,mm,kk)
          rup%e(:nn,mm,kk)=     -rup%e(:nn,mm,kk)-(mv*kv)*chi%e(:nn,mm,kk)
       enddo
    enddo
    call delsqh(chi,uz)
    if(.not.present(c)) then
      call allocate(w)
      call mulxm(uz,w)
      call mulxm(w,uz)
      call deallocate(w)
    endif
    uz%e = -uz%e

    end subroutine
!
!
!
    subroutine pc2vor(psi,chi,ror,rop,oz,c)
!   psi,chi to vorticity
!   returns (r*or,r*op,oz) if c is not present
!   returns (r*or,r*op,oz/(1-mu)^2 if c is present
!---------------------------------------------
    type(scalar):: psi,chi,ror,rop,oz
    integer:: mm,kk,nn,mv
    complex(p8),parameter:: iu=(0.0_p8,1.0_p8)
    real(p8):: kv
    type(scalar):: w
    character(len=1),optional:: c

    call del2(chi,oz)
    call xxdx(psi,ror)
    call xxdx(oz,rop)
    do mm=1,ntchop
    nn=nzchops(mm)
    mv=m(mm)
    do kk=1,nxchopdim
      kv=ak(mm,kk)
      ror%e(:nn,mm,kk)=-iu*mv*oz%e(:nn,mm,kk)+iu*kv*ror%e(:nn,mm,kk)
      rop%e(:nn,mm,kk)=       rop%e(:nn,mm,kk)-(kv*mv)*psi%e(:nn,mm,kk)
    enddo
    enddo
    call delsqh(psi,oz)
    if(.not.present(c)) then
      call allocate(w)
      call mulxm(oz,w)
      call mulxm(w,oz)
      call deallocate(w)
    endif
    oz%e = -oz%e

    end subroutine
!
!
!
    subroutine vprod(rur,rup,uz,ror,rop,oz)
!----------------------------------------------------
!   compute a vector product (wrapper)
!     (pr,pp,pz) = (ur,up,uz) x (or,op,oz)
!
!     input:
!       def h:=(1-mu)^2
!       rur,rup,uz: (r*ur,r*up,uz) 
!       ror,rop,oz: (r*or,r*op,oz)
!
!     output: 
!       (rpr,rpp,pz)=(r*pr,r*pp,pz)
!
!   input: any space
!   output: PHYSICAL space and the result returns in 
!       rpr => ror
!       rpp => rop
!       pz  => oz
!   contents of rur,rup,uz,ror,rop,oz destroyed
!----------------------------------------------------
    type(scalar):: rur,rup,uz
    type(scalar):: ror,rop,oz
    integer:: ni,nj,nk

    call tofp(rur)
    call tofp(rup)
    call tofp(uz)
    call tofp(ror)
    call tofp(rop)
    call tofp(oz)

    ni=size(rur%e,1)
    nj=size(rur%e,2)
    nk=size(rur%e,3)
    call vprodsub(rur%e(1,1,1),rup%e(1,1,1),uz%e(1,1,1), &
		  ror%e(1,1,1),rop%e(1,1,1),oz%e(1,1,1),ni,nj,nk)

    return
    end subroutine
!
!
!
    subroutine project(rur,rup,uz,psi,chi)
! -----------------------------------------------
!   project a vector in physical space (r*ur,r*up,uz) to
!   toroidal and poloidal field (psi,del2(chi)).
! -----------------------------------------------
    type(scalar):: rur,rup,uz,psi,chi
    real(p8),dimension(:,:),allocatable:: v,d,t,pfd
    complex(p8),dimension(:,:),allocatable:: w1,w2
    integer:: i,j,kk,mm,nn,mv,n
    real(p8)::kv
    complex(p8),parameter:: iu=(0.0_p8,1.0_p8)
    real(p8),dimension(nzh):: fff

    if(rur%space.ne.ppp_space .or. uz%space.ne.ppp_space) then
      print *,'project: rur,rup,uz not in ppp_space.'
      stop
    endif
    if(.not.associated(psi%e).or. .not.associated(chi%e)) then
      print *,'project: psi/chi not associated'
      stop
    endif
    if(psi%space.ne.fff_space) then
      call deallocate(psi)
      call allocate(psi)
    endif
    if(chi%space.ne.fff_space) then
      call deallocate(chi)
      call allocate(chi)
    endif
    if(size(psi%e,3).ne.nxchopdim.or.size(psi%e,3).ne.nxchopdim) then
      print *,'project:size of psi,chi wrong.'
      print *,'nxchopdim=',nxchopdim
      print *,'sizeof(psi)=[',size(psi%e,1),size(psi%e,2),size(psi%e,3),']'
      print *,'sizeof(psi)=[',size(psi%e,1),size(psi%e,2),size(psi%e,3),']'
      stop
    endif

    call horfft(rur,-1)
    call verfft(rur,-1)
    call horfft(rup,-1)
    call verfft(rup,-1)
    call horfft(uz,-1)
    call verfft(uz,-1)
    !taking care of the conjugate coefficients.
    !copying them over so that : k:[0,1,2,3,4,-3,-2,-1]
    rur%e(:nz,:ntchop,nxchop+1:nxchopdim)=rur%e(:nz,:ntchop,nxchopu:nx)
    rup%e(:nz,:ntchop,nxchop+1:nxchopdim)=rup%e(:nz,:ntchop,nxchopu:nx)
     uz%e(:nz,:ntchop,nxchop+1:nxchopdim)= uz%e(:nz,:ntchop,nxchopu:nx)

    allocate( v(nzchop,nzh) )
    allocate( d(nzchop,nzh) )
    allocate( t(nzchop,nzh) )
    allocate( pfd(nzh,nzchop+1) )
    allocate( w1(nzchop,nxchopdim) )
    allocate( w2(nzchop,nxchopdim) )
   
    w1=0
    fff=(1-tfm%x(:nzh)**2)/tfm%w(:nzh)
    do mm=1,ntchop
       mv=m(mm)
       nn=nzchops(mm)
       pfd(:nzh,:nn+1)=banmul( tfm%pf(:nzh,:nn+1,mm), &
            band_leg_xxdx(nn+1,mv,tfm%norm(:,mm)),2 )
       do i=1,nn
          n=max(1,mv+i-1)
          v(i,:nzh)= tfm%pf(:nzh,i,mm)/(n*(n+1)*fff)
          d(i,:nzh)= pfd(:nzh,i)/(n*(n+1)*fff)
          t(i,:nzh)= tfm%pf(:nzh,i,mm)*tfm%w(:nzh)
       enddo
       if(mv.ne.0) call eomul(v(:nn,:),rur%e(:nz,mm,:nxchopdim),w1(:nn,:))
       call oemul(d(:nn,:),rup%e(:nz,mm,:nxchopdim),w2(:nn,:))
       psi%e(:nn,mm,:) = -(iu*mv)*w1(:nn,:) -w2(:nn,:)
       !
       call oemul(d(:nn,:),rur%e(:nz,mm,:nxchopdim),w1(:nn,:))
       if(mv.ne.0) call eomul(v(:nn,:),rup%e(:nz,mm,:nxchopdim),w2(:nn,:))
       call eomul(t(:nn,:),uz%e(:nz,mm,:nxchopdim),chi%e(:nn,mm,:))
       do kk=1,nxchopdim
          kv=ak(mm,kk)
          chi%e(:nn,mm,kk)=(iu*kv)*w1(:nn,kk)+mv*kv*w2(:nn,kk)-chi%e(:nn,mm,kk)
       enddo
    enddo

    deallocate( v,d,t,pfd )
    deallocate( w1,w2 )
    call chopdo(psi)
    call chopdo(chi)
    
    psi%e(1,1,:)=psi%e(1,1,:)-calcat1(psi)/tfm%pf(1,1,1)
    
  end subroutine project
!
!
!
    subroutine nonlin(psi,chi,psin,chin)
!   compute nonlinear term
!   output returns in psin, chin
! ----------------------------------------------------
    type(scalar):: psi,chi,psin,chin
    type(scalar):: rur,rup,uz
    type(scalar):: ror,rop,oz
    type(scalar):: w

    call allocate(rur)
    call allocate(rup)
    call allocate( uz)
    call allocate(ror)
    call allocate(rop)
    call allocate( oz)
    !
    !dealiasing
    call chopset(3)
    call pc2vel(psi,chi,rur,rup,uz)
    !input and outputs are still FFF
    call pc2vor(psi,chi,ror,rop,oz)
    !input and outputs are still FFF
    call vprod(rur,rup,uz,ror,rop,oz)
    !vprod takes all to PPP and ... 
    !Remember that ror,rop,oz contain the result of the product, in PPP
    call chopset(-3)
    !
    call deallocate(rur)
    call deallocate(rup)
    call deallocate( uz)
    !
    call allocate(w)
    w%ln=0
    call project(ror,rop,oz,psin,w)
    !now psin = psi_nonlinear and w = delsq(chi_nonlinear)
    call idel2(w,chin)
    !now chin = chi_nonlinear

    call deallocate(ror)
    call deallocate(rop)
    call deallocate( oz)
    call deallocate( w )

    return
    end subroutine
!
!
!
  subroutine eomul(a,b,c)
! compute c(:ni,:nk)=a(:ni,:nj)*b(:nj,:nk) where a has symmetry, i.e.
!    a(2*i  :j)= a(2*i  :nj-j+1) and 
!    a(2*i-1:j)=-a(2*i-1:ni-j+1).
! only half of a should be given on input.
! ---------------------------------------------------------------
  real(p8),dimension(:,:):: a
  complex(p8),dimension(:,:):: b,c
  complex(p8),dimension(:,:),allocatable:: be,bo
  integer:: ni,nj,nk,njh

  ni  = size(a,1)
  njh = size(a,2)
  nj  = size(b,1)
  nk  = size(b,2)
  if(njh*2 .ne. nj) then
    print *,'eomul: size mismatch.'
    print *,'nj=',nj,' njh=',njh
    stop
  endif
  allocate( be(njh,nk) )
  allocate( bo(njh,nk) )
  be = ( b(1:njh,:) + b(nj:njh+1:-1,:) )
  bo = ( b(1:njh,:) - b(nj:njh+1:-1,:) )
  c(1::2,:)= a(1::2,:) .mul. be
  if (ni .gt. 1) then
	c(2::2,:)= a(2::2,:) .mul. bo
  endif
  deallocate( be,bo )

  return
  end subroutine
!
!
!
  subroutine oemul(a,b,c)
! compute c(:ni,:nk)=a(:ni,:nj)*b(:nj,:nk) where a has symmetry, i.e.
!    a(2*i  :j)=-a(2*i  :nj-j+1) and 
!    a(2*i-1:j)= a(2*i-1:ni-j+1).
! only half of a should be given on input.
! ---------------------------------------------------------------
  real(p8),dimension(:,:):: a
  complex(p8),dimension(:,:):: b,c
  complex(p8),dimension(:,:),allocatable:: be,bo
  integer:: ni,nj,nk,njh

  ni  = size(a,1)
  njh = size(a,2)
  nj  = size(b,1)
  nk  = size(b,2)
  if(njh*2 .ne. nj) then
    print *,'oemul: size mismatch.'
    print *,'nj=',nj,' njh=',njh
    stop
  endif
  allocate( be(njh,nk) )
  allocate( bo(njh,nk) )
  be = ( b(1:njh,:) + b(nj:njh+1:-1,:) )
  bo = ( b(1:njh,:) - b(nj:njh+1:-1,:) )
  c(1::2,:)= a(1::2,:) .mul. bo
  c(2::2,:)= a(2::2,:) .mul. be
  deallocate( be,bo )

  return
  end subroutine

!-----------
 end module 
!-----------
