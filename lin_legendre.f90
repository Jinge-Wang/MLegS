!
!---------------------
  module lin_legendre
!---------------------
  implicit none
  integer,parameter:: p8=selected_real_kind(p=11)
!======================================================================
  private
!
! Functions
! ---------
  public:: leg  
     !  evaluate legendre polynomial.
     !
     !  leg(m,n,x) evaluates the (unnormalized) legendre polynomial 
     !          of degree n and order m.  x can be either real scalar 
     !          or real one dimensional vector. The size of the result 
     !          is the same as the size of x.
         
  public:: nleg  
     !  nleg(m,n,x) evaluates normalized legendre polynomial.
     !          this routine is the normalized version of leg.

  public:: leg_norm     
     !  normalization factor. 
     !  leg_norm(ndim,m) with integer ndim and m retuns an one dimensional 
     !          real array containing the normalization factor of 
     !          legendre polynomials of order m.  
     !          if m(:) is one dimensional integer array, a two dimensional
     !          real array whose i-th column contains normalization factor for 
     !          order m(i).
     !          normalized function \Phi^m_n is obtained by multiplying
     !          the factor to Legendre functions P^m_n, i.e.,
     !          \Phi^m_n = norm(#,#)*P^m_n.
      
  public:: legtbl 
     !  table of legendre polynomial 
     !
     !  legtbl(x,ne,mdim) computes table of legendre polynomials
     !          and returns a three dimensional real array of size
     !          size(x)*ne*size(mdim).  x is real one dimensional array
     !          or a scalar containing the collocation points. 
     !          (ne-1) specifies the
     !          highest degree to be computed. mdim is one dimensional 
     !          integer array containing the orders of legendre polynomials.
     !
     !  legtbl(x,ne,mdim,norm) computes the table of normalized legendre
     !          polynomials. norm is the table of normalization factor and
     !          need to be set beforehand by leg_norm.
     !        
     !          if mdim is a integer, then legtbl returns a two dimensional
     !          table for the polynomials of order mdim.
     ! 
     !          example to obtain a normalized table: 
     !              norm = leg_norm(ni,(/0,mmax-1/))
     !              tbl  = legtbl(x,ni,(/0,mmax-1/),norm)

  public:: leg_xm  	! operator to multiply (1-x)
  public:: leg_xp  	! operator to multiply (1+x)
  public:: leg_x   	! operator to multiply x
     !  leg_xm(ni,nj,m) returns a real two dimensional array of size ni*nj
     !          representing (1-x) in legendre space. the operator is 
     !          tridiagonal.
     !
     !  leg_xm(ni,ni,m,norm) returns a normalized version of the operator,
     !          where norm is the normalization table.
     !
     !  leg_xp and leg_x with similar calling procedures return 
     !          the operators for (1+x) and x.

  public:: leg_xxdx	
     !  operator to apply (1-x^2)d/dx
     !
     !  leg_xxdx(ni,nj,m) returns a real two dimensional array of size ni*nj
     !          representing (1-x^2)d/dx in legendre function space. 
     !          the operator is tridiagonal.
     !
     !  leg_xxdx(ni,ni,m,norm) returns a normalized version of the operator.
     !          norm is the normalization table.

  public:: leg_rat_del2h
     !  horizontal Laplacian operator for mapped rational functions.
     !  
     !  leg_rat_del2h(ni,nj,m,eel) returns ni*nj unnormalized operator.
     !  where m is the azimuthal wave number, eel is map parameter.
     !  
     !  leg_rat_del2h(ni,nj,m,eel,norm) returns normalized operator.
     ! 
     !  (note: leg_rat_del2 is not supported. it can be obtained by
     !   mtrx(band_leg_rat_del2(...),3) where mtrx() is defined in
     !   bandmat.f90).

  public:: band_leg_xm     ! tridiagonal matrix to multiply (1-x)
  public:: band_leg_xp     ! tridiagonal matrix to multiply (1+x)
  public:: band_leg_x      ! tridiagonal matrix to multiply x
     !  band_leg_xm(ni,m) returns a ni*3 array containing the nonzero part of
     !          (1-x) operator in legendre space. the operator is 
     !          tridiagonal and the diagonal is band_leg_xm(:,2). 
     !          band_leg_xm(i,:) contains the i-th row of the full matrix.
     !
     !  band_leg_xm(ni,m,norm) returns a normalized version of the matrix,
     !          where norm is the normalization table.
     !
     !  band_leg_xp and leg_x with similar calling procedures return 
     !          the tridiagonal band matrices for (1+x) and x.

  public:: band_leg_xxdx	
     !  tridiagonal matrix to apply (1-x^2)d/dx
     !
     !  band_leg_xxdx(ni,m) returns a real tridiagonal array of size ni*3
     !          representing (1-x^2)d/dx in legendre function space. 
     !          band_leg_xxdx(:,2) corresponds to diagonal elements.
     !
     !  band_leg_xxdx(ni,m,norm) returns a normalized version of the matrix.
     !          norm is the normalization table.

  public:: band_leg_rat_del2h
     !  horizontal Laplacian operator for mapped rational functions.
     !  
     !  band_leg_rat_del2h(ni,m,eel) returns ni*5 unnormalized operator.
     !  band width is 5. 
     !  m is the azimuthal wave number.
     !  eel is the map parameter.
     !  
     !  band_leg_rat_del2h(ni,m,eel,norm) returns normalized operator.

  public:: band_leg_rat_del2
     !  Laplacian operator for mapped rational functions.
     !  band_leg_rat_del2(ni,m,ak,eel) returns ni*5 unnormalized operator.
     !  band width is 5.
     !  band_leg_rat_del2(ni,m,ak,eel,norm) returns normalized operator.
     !  ak is the axial wave number.
     !  if m and ak are both zero, logarithmic term is assumed at the 
     !  lowest coefficient (this is not the case for band_leg_rat_del2h).
     !  

  public:: leg_rat_proj
     !  projection operator for rational Legendre functions.
     !
     !  leg_rat_proj(x,w,pf,pfd,m,ak) returns a projection operator.
     !  Let ni=size(pf,2), nj=size(x).
     !  The size of returning matrix is complex(p8), (2*ni,3*nj).
     !  x(:nj) contains the collocation points. 
     !  w(:nj) contains the weight factor for quadrature. 
     !  pf(i,n) is the table of  P^m_n(mu_i)
     !  pfd(i,n) is the table of (1-mu_i^2)(d/dmu)P^m_n(mu_i).
     !  m,ak are the azimuthal and axial wave number.
     
! ------------------------------------------------------------------

!
! Subroutines:
! ------------
  public:: leg_zeros   
	! zeros and weights for gauss-legendre integration
	!
        ! leg_zeros(x) computes the zeros of legendre polynomials
	!      of degree size(x). 
	!
	! leg_zeros(x,w) returns the zeros
	!      and weights of gauss_legendre quadrature.
	!
	! leg_zeros(x,w,r,eel) additionally computes the 
	!      collocation points of rational legendre functions
	!      of with map parameter eel. 
	!      x,w,r should all be real(p8) one dimensional array 
	!      and have the same size. 
	!      eel should be real(p8).





!======================================================================

  interface leg
    module procedure epol
    module procedure epolv
  end interface

  interface nleg
    module procedure enpol
    module procedure enpolv
  end interface

  interface leg_norm
    module procedure leg_norm1
    module procedure leg_normt
  end interface

  interface legtbl
    module procedure legtbl3
    module procedure legtbl2
    module procedure legtbl1
  end interface


 contains
!--------
!
!
!
      subroutine gauleg(x1,x2,x,w,n)
!     ------------------------------------
!     a routine from Numerical Recipes:
!     calculate the collocation pts (x) and weights(w)
!     of Gauss-Legendre n-point quadratures
!     ------------------------------------
      real(p8),intent(in):: x1,x2 ! x1 and x2 are lower and upper bounds
      real(p8),dimension(:) :: x,w
      integer:: n

      real(p8),parameter ::eps=1.e-13
      real(p8):: xm,xl,p1,p2,p3,z,z1,pp
      integer:: m,i,j

      m=(n+1)/2
      xm=0.5_p8*(x2+x1)
      xl=0.5_p8*(x2-x1)
      do i=1,m ! refinement by Newton's method
        z=cos(3.141592653589793238462643_p8*(i-.25_p8)/(n+.5_p8))
1       continue
          p1=1
          p2=0
          do j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1)*z*p2-(j-1)*p3)/j
          enddo
          pp=n*(z*p1-p2)/(z*z-1)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps)go to 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2*xl/((1-z*z)*pp*pp)
        w(n+1-i)=w(i)
      enddo
      return
      end subroutine
!
!
!
      subroutine leg_zeros(zero,wei,rzero,eel)
!    -------------------------------------------------------------------------
!     finds the zeros and weights of the legendre polynomial of degree nn, 
!     where n=size(zero) is the degree of legendre polynomial.
!     number of zeros = n.
!     gauleg uses a precision-increased version 
!     of the routine found in Numerical Recipes 
!    -------------------------------------------------------------------------
      real(p8),dimension(:):: zero	      !zeros of legendre polynomials
      real(p8),dimension(:),optional:: wei    !weights for quadrature
      real(p8),dimension(:),optional:: rzero  !colloc.pts for rational mapping
      real(p8),optional:: eel                 !map parameter

      real(p8),dimension(size(zero)):: wk
!
      if(present(rzero).neqv.present(eel)) then
	print *,'leg_zeros: both rzeros and eel have to be specified.'
	print *,'to obtain the mapped collocation points.'
	stop
      endif
!
      if(present(wei)) then
        call gauleg(-1.0_p8,1.0_p8,zero,wei,size(zero))
      else
        call gauleg(-1.0_p8,1.0_p8,zero,wk,size(zero))
      endif
!
      if(present(rzero)) then
        rzero = eel*sqrt((1+zero)/(1-zero))
      endif
!
      end subroutine
!
!
!
      function leg_norm1(ndim,m)
!     --------------------------------------------
!     table of normalization factor
!     ex: call leg_norm(64,m)
!     --------------------------------------------
      integer,intent(in) &
	:: ndim           	   ! number of coefficients
      integer,intent(in) & 	   ! m
	:: m   			   ! ex: (/0:mmax/)
      real(p8), dimension(ndim) &  ! return array
	:: leg_norm1

      leg_norm1=(/ leg_normt(ndim, (/abs(m)/) ) /)
      end function
!
!
!
      function leg_normt(ndim,mdim)
!     ----------------------------------------------------
!     table of normalization factor
!     ex: call leg_norm(64,(/0,mdim/))
!
!     if n=nn-1, m=mdim(mm), leg_normt(nn,mm)*P^m_n(x) is 
!     an normalized function.
!     ----------------------------------------------------
      integer,intent(in) &
	:: ndim           	   ! number of coefficients
      integer,intent(in), &
      dimension(:) :: mdim 	   ! specify m's in increasing order
      real(p8), &		   ! return array
      dimension(ndim,size(mdim)) &
	:: leg_normt

      integer ::nn,mm,n,m,ms,me
      real(p8),dimension(:),allocatable:: wk

      !me=mdim(size(mdim))
      me=maxval(abs(mdim))
      allocate(wk(0:me))

      wk(0) = 0.5_p8
      do m=1,me
	wk(m) = wk(m-1)*(2*m+1.0_p8)/(2*m*(2*m-1.0_p8)**2)
      enddo
!
!     > avold underflow
      do mm=1,size(mdim)
	m=abs(mdim(mm))
	leg_normt(1,mm) = wk(m)*100.0_p8**(m)
      enddo
!
      do mm=1,size(mdim)
      m = abs(mdim(mm))
      do nn=2,ndim
	n = m+(nn-1)
	leg_normt(nn,mm) = leg_normt(nn-1,mm)*(2*n+1.0_p8)/(2*n-1.0_p8)*   &
                       (n-m)/(n+m)
      enddo
      enddo
!
      do mm=1,size(mdim)
      m=abs(mdim(mm))
      leg_normt(:,mm) = sqrt(leg_normt(:,mm))*10.0_p8**(-m)
      enddo

      deallocate(wk)

      end function

!
!
!
      function legtbl1(x,ne,m,norm)
!-------------------------------------------------------------------------
!  legendre polynomial table for only one m and x
!-------------------------------------------------------------------------
      real(p8),intent(in):: x                   ! pts to evaluate
      integer,intent(in) :: ne			! max degree
      integer,intent(in) :: m                ! m
      real(p8),dimension(ne) :: legtbl1         ! table
      real(p8),dimension(ne,1) :: norm1         ! table
      real(p8),dimension(:), &
		intent(in),optional   :: norm   ! normalization table
      integer:: n

      if(present(norm)) then
	norm1(:,1)=norm(:ne)
        legtbl1 = (/legtbl3( (/x/), ne, (/abs(m)/), norm1) /)
      else
        legtbl1 = (/ legtbl3( (/x/), ne, (/abs(m)/) ) /)
      endif

      end function
!
!
!
      function legtbl2(x,ne,m,norm)
!-------------------------------------------------------------------------
!  legendre polynomial table for only one m 
!-------------------------------------------------------------------------
      real(p8),dimension(:),intent(in):: x      ! pts to evaluate
      integer,intent(in) :: ne			! max degree
      integer,intent(in) :: m                ! m
      real(p8),dimension(size(x),ne) :: legtbl2 ! table
      real(p8),dimension(:), &
		intent(in),optional   :: norm   ! normalization table

      if(present(norm)) then
        legtbl2 = reshape(legtbl3(x, ne, (/abs(m)/),  &
			reshape(norm,(/size(norm),1/)) ) &
		        , (/size(x),ne/))
      else
        legtbl2 = reshape(legtbl3(x, ne, (/abs(m)/)), (/size(x),ne/))
      endif

      end function
!
!
!
      function legtbl3(x,ne,mdim,norm)
!-------------------------------------------------------------------------
!  Make a table of associate legendre polynomials
!  at location x(:). Normalization is not performed if norm is not present.
!  Evaluation is performed for order mdim(:).
!  result is a 3-dim array and legtbl3(:,nn,mm) are the values
!  of legendre polynomials of order mdim(mm) and degree (nn-1)+mdim(mm).
!  ex: legtbl3(nend,(/0:mend/),x)
!     
!     ex: normalize=leg_norm(ni,(/0,mmax-1/))
!         tbl =legtbl3(x,ni,(/0,mmax-1/),normalize)
!-------------------------------------------------------------------------
      real(p8),dimension(:),intent(in):: x      ! pts to evaluate
      integer,intent(in) :: ne			! max degree
      integer,intent(in),dimension(:) :: mdim   ! array containing m's
      real(p8),dimension(:,:), &
		intent(in),optional   :: norm   ! normalization table
      real(p8),dimension(size(x),ne,size(mdim))&
	   			      :: legtbl3 ! table
!
      integer :: m,me,mm,nn,n
!
      me = size(mdim)
!
!     > orthogonal-unnormalized polynomial table
!     > starting values of legtbl3
      do mm=1,me
        m = abs(mdim(mm))
	if(m.eq.0) then
	  legtbl3(:,1,mm) = 1.0_p8
	else
	  legtbl3(:,1,mm) = fact1(m)*(1-x**2)**(0.5_p8*m)
	endif
      enddo

      do mm=1,me
        m = abs(mdim(mm))
	legtbl3(:,2,mm) = ((2*m+1)*x)*legtbl3(:,1,mm)
      enddo
!
!     > recurrence legtbl3
      do mm=1,me
      m = abs(mdim(mm))
      do nn=3,ne
        n = m+(nn-1)
	  legtbl3(:,nn,mm) = 1.0_p8/(n-m)*(   &
      	     legtbl3(:,nn-1,mm)*((2*n-1)*x)  &
      	    -legtbl3(:,nn-2,mm)*(n+m-1))
      enddo
      enddo
!
!     > normalize
      if(present(norm)) then
	if(size(norm,1).lt.ne .or. size(norm,2).lt.me) then
	  print *,'size of norm too small'
	  print *,'norm=',size(norm,1),'x',size(norm,2)
	  stop
	endif
	do mm=1,me
	do nn=1,ne
	  legtbl3(:,nn,mm)=legtbl3(:,nn,mm)*norm(nn,mm)
	enddo
	enddo
      endif
      end function
!
!
!
  function leg_xm(ni,nj,m,norm)
! -----------------------------------------------------
! operator to multiply (1-x) in function space
! m is the order of legendre polynomials to operate on.
! if norm is present in the argument, mormalized 
! operator returns.
! -----------------------------------------------------
  integer,intent(in):: ni,nj,m
  real(p8),dimension(ni,nj):: leg_xm
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  am=abs(m)
  leg_xm=0
  do nn=1,ni
    n = am + (nn-1)
    if(nn-1.ge.1)  leg_xm(nn,nn-1) = -(n-am)/(2*n-1.0_p8)
    if(nn.le.nj)   leg_xm(nn,nn)   = 1.0
    if(nn+1.le.nj) leg_xm(nn,nn+1) = -(n+am+1.0_p8)/(2*n+3)
  enddo

  if(present(norm)) then
    if(size(norm) < max(ni,nj)) then
      print *,'leg_xm: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni,'  nj=',nj
      stop
    endif
    do i=1,ni
    do j=1,nj
      leg_xm(i,j)=leg_xm(i,j)*(norm(j)/norm(i))
    enddo
    enddo
  endif
  end function
!
!
!
  function leg_xp(ni,nj,m,norm)
! -----------------------------------------------------
! operator to multiply (1+x) in function space
! m is the order of legendre polynomials to operate on.
! if norm is present in the argument, mormalized 
! operator returns.
! -----------------------------------------------------
  integer,intent(in):: ni,nj,m
  real(p8),dimension(ni,nj):: leg_xp
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  am=abs(m)
  leg_xp=0
  do nn=1,ni
    n = am + (nn-1)
    if(nn-1.ge.1)  leg_xp(nn,nn-1) = (n-am)/(2*n-1.0_p8)
    if(nn.le.nj)   leg_xp(nn,nn)   = 1.0_p8
    if(nn+1.le.nj) leg_xp(nn,nn+1) = (n+am+1.0_p8)/(2*n+3)
  enddo

  if(present(norm)) then
    if(size(norm) < max(ni,nj)) then
      print *,'leg_xp: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni,'  nj=',nj
      stop
    endif
    do i=1,ni
    do j=1,nj
      leg_xp(i,j)=leg_xp(i,j)*(norm(j)/norm(i))
    enddo
    enddo
  endif
  end function
!
!
!
  function leg_x(ni,nj,m,norm)
! -----------------------------------------------------
! operator to multiply x in function space
! m is the order of legendre polynomials to operate on.
! if norm is present in the argument, mormalized 
! operator returns.
! -----------------------------------------------------
  integer,intent(in):: ni,nj,m
  real(p8),dimension(ni,nj):: leg_x
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  am=abs(m)
  leg_x=0
  do nn=1,ni
    n = am + (nn-1)
    if(nn-1.ge.1)  leg_x(nn,nn-1) = (n-am)/(2*n-1.0_p8)
    if(nn+1.le.nj) leg_x(nn,nn+1) = (n+am+1.0_p8)/(2*n+3)
  enddo

  if(present(norm)) then
    if(size(norm) < max(ni,nj)) then
      print *,'leg_x: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni,'  nj=',nj
      stop
    endif
    do i=1,ni
    do j=1,nj
      leg_x(i,j)=leg_x(i,j)*(norm(j)/norm(i))
    enddo
    enddo
  endif
  end function
!
!
!
  function leg_xxdx(ni,nj,m,norm)
! -----------------------------------------------------
! operator to apply (1-x^2)*d/dx in function space
! m is the order of legendre polynomials to operate on.
! if norm is present in the argument, mormalized 
! operator returns.
! -----------------------------------------------------
  integer,intent(in):: ni,nj,m
  real(p8),dimension(ni,nj):: leg_xxdx
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  am=abs(m)
  leg_xxdx=0
  do nn=1,ni
    n = am + (nn-1)
    if(nn-1.ge.1)  leg_xxdx(nn,nn-1) = -(n-1.0_p8)*(n-am)/(2*n-1.0_p8)
    if(nn+1.le.nj) leg_xxdx(nn,nn+1) = (n+2)*(n+am+1.0_p8)/(2*n+3.0_p8)
  enddo

  if(present(norm)) then
    if(size(norm) < max(ni,nj)) then
      print *,'leg_xxdx: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni,'  nj=',nj
      stop
    endif
    do i=1,ni
    do j=1,nj
      leg_xxdx(i,j)=leg_xxdx(i,j)*(norm(j)/norm(i))
    enddo
    enddo
  endif
  end function
!
!
!
  function leg_rat_del2h(ni,nj,m,eel,norm)
! ------------------------------------------------------------------
! horizontal Laplacian operator for rational legendre equation.
! m is the azimuthal wave number.
! eel is the map parameter.
! if norm is present in the argument, mormalized operator returns.
! ------------------------------------------------------------------
  integer,intent(in):: ni,nj,m
  real(p8),dimension(ni,nj):: leg_rat_del2h
  real(p8) :: eel
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  real(p8) :: eel2
  real(p8),parameter:: s=0.0_p8
  am=abs(m)
!
  leg_rat_del2h = 0
  eel2=eel**2
  do nn=1,ni
    n = am + (nn-1)
    if(nn-2.ge.1) leg_rat_del2h(nn,nn-2) = -(n-am-1.0_p8)*(n-am)  &
                  *(n-2+s)*(n-1+s)/(2*n-3.)/(2*n-1.0_p8)/eel2
    if(nn-1.ge.1) leg_rat_del2h(nn,nn-1) = 2.0_p8*n*(n-am)*(n-1+s) &
		  /(2*n-1.0_p8) /eel2
    leg_rat_del2h(nn,nn) = (-2.0_p8*n*(n+1._p8)*(3.0_p8*n*n+3*n-am**2-2)  &
                  +2*s*(s-2)*(n*n+n+m*m-1._p8))/(2*n-1.0_p8)/(2*n+3.0_p8)/eel2  
    if(nn+1.le.nj) leg_rat_del2h(nn,nn+1) =  &
                  +2.0_p8*(n+1)*(n+am+1.0_p8)*(n+2-s)/(2*n+3.0_p8)/eel2
    if(nn+2.le.nj) leg_rat_del2h(nn,nn+2) =  &
                   -(n+am+1.0_p8)*(n+am+2.0_p8)*(n+3-s)*(n+2-s) &
		   /(2*n+3.0_p8)/eel2/(2*n+5.0_p8)
  enddo

  if(present(norm)) then
    if(size(norm) < max(ni,nj)) then
      print *,'leg_rat_del2h: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni,'  nj=',nj
      stop
    endif
    do i=1,ni
    do j=1,nj
      leg_rat_del2h(i,j)=leg_rat_del2h(i,j)*(norm(j)/norm(i))
    enddo
    enddo
  endif
  end function
!
!
!
  function band_leg_xm(ni,m,norm)
! -----------------------------------------------------
! band matrix to multiply (1-x) in function space.
! m is the order of legendre polynomials to operate on.
! if norm is present in the argument, mormalized 
! operator returns. diagonal corresponds to band_leg_xm(:,2).
! -----------------------------------------------------
  integer,intent(in):: ni,m
  real(p8),dimension(ni,3):: band_leg_xm
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  am=abs(m)
  band_leg_xm=0
  do nn=1,ni
    n = am + (nn-1)
    if(nn-1.ge.1)  band_leg_xm(nn,1) = -(n-am)/(2*n-1.0_p8)
    if(nn.le.ni)   band_leg_xm(nn,2) = 1.0
    if(nn+1.le.ni) band_leg_xm(nn,3) = -(n+am+1.0_p8)/(2*n+3)
  enddo

  if(present(norm)) then
    if(size(norm) < ni) then
      print *,'band_leg_xm: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni
      stop
    endif
    do i=1,ni
      if(i.gt.1) band_leg_xm(i,1)=band_leg_xm(i,1)*(norm(i-1)/norm(i))
      if(i.lt.ni) band_leg_xm(i,3)=band_leg_xm(i,3)*(norm(i+1)/norm(i))
    enddo
  endif
  end function
!
!
!
  function band_leg_xp(ni,m,norm)
! -----------------------------------------------------
! band matrix to multiply (1-x) in function space.
! m is the order of legendre polynomials to operate on.
! if norm is present in the argument, mormalized 
! operator returns. diagonal corresponds to band_leg_xp(:,2).
! -----------------------------------------------------
  integer,intent(in):: ni,m
  real(p8),dimension(ni,3):: band_leg_xp
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  am=abs(m)
  band_leg_xp=0
  do nn=1,ni
    n = am + (nn-1)
    if(nn-1.ge.1)  band_leg_xp(nn,1) = (n-am)/(2*n-1.0_p8)
    if(nn.le.ni)   band_leg_xp(nn,2) = 1.0
    if(nn+1.le.ni) band_leg_xp(nn,3) = (n+am+1.0_p8)/(2*n+3)
  enddo

  if(present(norm)) then
    if(size(norm) < ni) then
      print *,'band_leg_xp: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni
      stop
    endif
    do i=1,ni
      if(i.gt.1)  band_leg_xp(i,1)=band_leg_xp(i,1)*(norm(i-1)/norm(i))
      if(i.lt.ni) band_leg_xp(i,3)=band_leg_xp(i,3)*(norm(i+1)/norm(i))
    enddo
  endif
  end function
!
!
!
  function band_leg_x(ni,m,norm)
! -----------------------------------------------------
! band matrix to multiply (1-x) in function space.
! m is the order of legendre polynomials to operate on.
! if norm is present in the argument, mormalized 
! operator returns. diagonal corresponds to band_leg_x(:,2).
! -----------------------------------------------------
  integer,intent(in):: ni,m
  real(p8),dimension(ni,3):: band_leg_x
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  am=abs(m)
  band_leg_x=0
  do nn=1,ni
    n = am + (nn-1)
    if(nn-1.ge.1)  band_leg_x(nn,1) =  (n-am)/(2*n-1.0_p8)
    if(nn+1.le.ni) band_leg_x(nn,3) =  (n+am+1.0_p8)/(2*n+3)
  enddo

  if(present(norm)) then
    if(size(norm) < ni) then
      print *,'band_leg_x: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni
      stop
    endif
    do i=1,ni
      if(i.gt.1)  band_leg_x(i,1)=band_leg_x(i,1)*(norm(i-1)/norm(i))
      if(i.lt.ni) band_leg_x(i,3)=band_leg_x(i,3)*(norm(i+1)/norm(i))
    enddo
  endif
  end function
!
!
!
  function band_leg_xxdx(ni,m,norm)
! -----------------------------------------------------
! matrix to apply (1-x^2)*d/dx in function space
! m is the order of legendre polynomials to operate on.
! if norm is present in the argument, mormalized 
! operator returns. diagonal corresponds to band_leg_x(:,2).
! -----------------------------------------------------
  integer,intent(in):: ni,m
  real(p8),dimension(ni,3):: band_leg_xxdx
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  am=abs(m)
  band_leg_xxdx=0
  do nn=1,ni
    n = am + (nn-1)
    if(nn-1.ge.1)  band_leg_xxdx(nn,1) = -(n-1.0_p8)*(n-am)/(2*n-1.0_p8)
    if(nn+1.le.ni) band_leg_xxdx(nn,3) = (n+2)*(n+am+1.0_p8)/(2*n+3.0_p8)
  enddo

  if(present(norm)) then
    if(size(norm) < ni) then
      print *,'band_leg_xxdx: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni
      stop
    endif
    do i=1,ni
      if(i.gt.1)  band_leg_xxdx(i,1)=band_leg_xxdx(i,1)*(norm(i-1)/norm(i))
      if(i.lt.ni) band_leg_xxdx(i,3)=band_leg_xxdx(i,3)*(norm(i+1)/norm(i))
    enddo
  endif
  end function
!
!
!
  function band_leg_rat_del2(ni,m,ak,eel,norm)
! ------------------------------------------------------------------
! Laplacian operator for rational legendre equation.
! m is the azimuthal wave number.
! ak is the axial Fourier wave number.
! eel is the map parameter.
! if m and ak are both zero, logarithmic term is assumed at the 
! lowest coefficient.
! if norm is present in the argument, mormalized operator returns.
! diagonal corresponds to band_leg_rat_del2(:,3).
! ------------------------------------------------------------------
  integer,intent(in):: ni,m
  real(p8),dimension(ni,5):: band_leg_rat_del2
  real(p8) :: eel,ak
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  real(p8) :: eel2
  real(p8),parameter:: s=0.0_p8
  am=abs(m)
  band_leg_rat_del2 = 0
  eel2=eel**2
  do nn=1,ni
    n = am + (nn-1)
    if(nn-2.ge.1) band_leg_rat_del2(nn,1) = -(n-am-1.0_p8)*(n-am)  &
                  *(n-2+s)*(n-1+s)/(2*n-3.)/(2*n-1.0_p8)/eel2
    if(nn-1.ge.1) band_leg_rat_del2(nn,2) = 2.0_p8*n*(n-am)*(n-1+s) &
		  /(2*n-1.0_p8) /eel2
    band_leg_rat_del2(nn,3) = (-2.0_p8*n*(n+1._p8)*(3.0_p8*n*n+3*n-am**2-2)  &
                  +2*s*(s-2)*(n*n+n+m*m-1._p8))/(2*n-1.0_p8)/(2*n+3.0_p8)/eel2  
    if(nn+1.le.ni) band_leg_rat_del2(nn,4) =  &
                  +2.0_p8*(n+1)*(n+am+1.0_p8)*(n+2-s)/(2*n+3.0_p8)/eel2
    if(nn+2.le.ni) band_leg_rat_del2(nn,5) =  &
                   -(n+am+1.0_p8)*(n+am+2.0_p8)*(n+3-s)*(n+2-s) &
		   /(2*n+3.0_p8)/eel2/(2*n+5.0_p8)
  enddo
  band_leg_rat_del2(:,3)= band_leg_rat_del2(:,3)-ak**2

  if(present(norm)) then
    if(size(norm) < ni) then
      print *,'band_leg_rat_del2: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni
      stop
    endif
    do i=1,ni
      if(i.gt.2) band_leg_rat_del2(i,1)= &
		 band_leg_rat_del2(i,1)*(norm(i-2)/norm(i))
      if(i.gt.1) band_leg_rat_del2(i,2)= &
		 band_leg_rat_del2(i,2)*(norm(i-1)/norm(i))
      if(i.lt.ni) band_leg_rat_del2(i,4)= &
		  band_leg_rat_del2(i,4)*(norm(i+1)/norm(i))
      if(i.lt.ni-1) band_leg_rat_del2(i,5)= &
		    band_leg_rat_del2(i,5)*(norm(i+2)/norm(i))
    enddo
    if(ak.eq.0 .and. m.eq.0) then
      band_leg_rat_del2(1,3)=  4.0_p8 /3/eel2/norm(1)
      band_leg_rat_del2(2,2)= -2.0_p8   /eel2/norm(2)
      band_leg_rat_del2(3,1)=  2.0_p8 /3/eel2/norm(3)
    endif
  else
    if(ak.eq.0 .and. m.eq.0) then
      band_leg_rat_del2(1,3)=  4.0_p8 /3/eel2
      band_leg_rat_del2(2,2)= -2.0_p8   /eel2
      band_leg_rat_del2(3,1)=  2.0_p8 /3/eel2
    endif
  endif
  end function
!
!
!
  function band_leg_rat_del2h(ni,m,eel,norm)
! ------------------------------------------------------------------
! horizontal Laplacian operator for rational legendre equation.
! m is the azimuthal wave number.
! eel is the map parameter.
! if norm is present in the argument, mormalized operator returns.
! diagonal corresponds to band_leg_rat_del2h(:,3).
! ------------------------------------------------------------------
  integer,intent(in):: ni,m
  real(p8),dimension(ni,5):: band_leg_rat_del2h
  real(p8) :: eel
  real(p8),dimension(:),optional:: norm
  integer :: nn,n,i,j,am
  real(p8) :: eel2
  real(p8),parameter:: s=0.0_p8
  am=abs(m)
!
  band_leg_rat_del2h = 0
  eel2=eel**2
  do nn=1,ni
    n = am + (nn-1)
    if(nn-2.ge.1) band_leg_rat_del2h(nn,1) = -(n-am-1.0_p8)*(n-am)  &
                  *(n-2+s)*(n-1+s)/(2*n-3.)/(2*n-1.0_p8)/eel2
    if(nn-1.ge.1) band_leg_rat_del2h(nn,2) = 2.0_p8*n*(n-am)*(n-1+s) &
		  /(2*n-1.0_p8) /eel2
    band_leg_rat_del2h(nn,3) = (-2.0_p8*n*(n+1._p8)*(3.0_p8*n*n+3*n-m**2-2)  &
                  +2*s*(s-2)*(n*n+n+m*m-1._p8))/(2*n-1.0_p8)/(2*n+3.0_p8)/eel2  
    if(nn+1.le.ni) band_leg_rat_del2h(nn,4) =  &
                  +2.0_p8*(n+1)*(n+am+1.0_p8)*(n+2-s)/(2*n+3.0_p8)/eel2
    if(nn+2.le.ni) band_leg_rat_del2h(nn,5) =  &
                   -(n+am+1.0_p8)*(n+am+2.0_p8)*(n+3-s)*(n+2-s) &
		   /(2*n+3.0_p8)/eel2/(2*n+5.0_p8)
  enddo

  if(present(norm)) then
    if(size(norm) < ni) then
      print *,'band_leg_rat_del2h: size(norm) too small'
      print *,'size(norm)=',size(norm)
      print *,'ni=',ni
      stop
    endif
    do i=1,ni
      if(i.gt.2) band_leg_rat_del2h(i,1)= &
		 band_leg_rat_del2h(i,1)*(norm(i-2)/norm(i))
      if(i.gt.1) band_leg_rat_del2h(i,2)= &
		 band_leg_rat_del2h(i,2)*(norm(i-1)/norm(i))
      if(i.lt.ni) band_leg_rat_del2h(i,4)= &
		  band_leg_rat_del2h(i,4)*(norm(i+1)/norm(i))
      if(i.lt.ni-1) band_leg_rat_del2h(i,5)= &
		    band_leg_rat_del2h(i,5)*(norm(i+2)/norm(i))
    enddo
  endif
  end function
!
!
!
  function leg_rat_proj(x,w,pf,pfd,m,ak)
! ------------------------------------------------------------------
! projection operator for rational legendre equation.
! m,ak is the azimuthal and axial wave number.
! Let ni=size(pf,2), nj=size(x).
! pf(i,n) is the table of  P^m_n(mu_i)
! pfd(i,n) is the table of (1-mu_i^2)(d/dmu)P^m_n(mu_i).
! x(:nj) is the collocation points. 
! the size of returning matrix is complex(p8), (2*ni,3*nj).
! ------------------------------------------------------------------
  real(p8),dimension(:,:),intent(in):: pf,pfd
  real(p8),dimension(:),intent(in):: x,w
  complex(p8),dimension(2*size(pf,2),3*size(pf,1)):: leg_rat_proj
  integer,intent(in):: m
  real(p8),intent(in):: ak
  integer :: ni,nj,njh
  integer :: i,j,n
  complex(p8),parameter:: iu=(0.0_p8,1.0_p8)
  real(p8),dimension(:,:),allocatable:: v,d,t

  ni=size(pf,2)
  nj=size(x,1)
  if(size(pfd,1).ne.nj) then
    print *,'leg_rat_proj: pfd size not equal to sizeof(pf).'
    stop
  endif

  allocate( v(ni,nj) )
  allocate( d(ni,nj) )
  allocate( t(ni,nj) )
  do j=1,nj
  do i=1,ni
    n=max(1,abs(m)+i-1)
    v(i,j)= pf(j,i)/(n*(n+1)*(1-x(j)**2))*w(j)
    d(i,j)= pfd(j,i)/(n*(n+1)*(1-x(j)**2))*w(j)
    t(i,j)= pf(j,i)*w(j)
  enddo
  enddo
  if(m.eq.0) then
    v(1,:)=0
  endif
  leg_rat_proj(1:ni,     1:nj)        = -iu*m*v
  leg_rat_proj(1:ni,     nj+1:2*nj)   = -d
  leg_rat_proj(1:ni,     2*nj+1:3*nj) = 0
  leg_rat_proj(ni+1:2*ni,1:nj)        = iu*ak*d
  leg_rat_proj(ni+1:2*ni,nj+1:2*nj)   = m*ak*v
  leg_rat_proj(ni+1:2*ni,2*nj+1:3*nj) = -t

  deallocate( v,d,t )

  end function
!
!
!
      real(p8) function fact1(m)
!    ----------------------------------------------------
!     compute (2m)!/2^m/m! or (2m-1)!!
!    ----------------------------------------------------
      integer :: m
      real(p8),dimension(129),save:: val=0
      integer :: i

      if(m.gt.129) then
        print *,'fact1: argument out of range'
        stop
      endif
      
      if(val(2).ne.3.0) then
        val(1) = 1.0
        do i=2,129
  	  val(i) = val(i-1)*(2*i-1)
        enddo
      endif
!
      fact1 = val(abs(m))
      return
      end function
!
!
!
      real(p8) function cfact(m,n)
!    ----------------------------------------------------
!     compute (2*n+1)/2 *(n-m)! / (n+m)!
!    ----------------------------------------------------
      integer :: m,n
      integer :: i,am
      real(p8) :: aaa
      am = abs(m)

      aaa = 1.0
      do i=(n-am)+1,am+n
	aaa = aaa/i
      enddo    
      cfact = aaa*(2*n+1.0_p8)/2.0_p8
      return
      end function
!
!
!
      function enpolv(m,n,x)
!----------------------------------------------------------
!  evaluate the orthonormalized polynomial
!  m,n mathematical numbers.(not fortran coefficients)
!  this routine is very slow
!----------------------------------------------------------
      integer,intent(in) :: m,n
      real(p8),intent(in),dimension(:):: x
      real(p8),dimension(size(x)):: enpolv
      real(p8):: fac
      integer:: am
      am = abs(m)

      if(n.lt.am) then
	enpolv=0
	return
      endif
      enpolv = epolv(am,n,x)
      fac = cfact(am,n)
      enpolv = enpolv * sqrt(fac) 
      end function
!
!
!
      real(p8) function enpol(m,n,x)
!---------------------------------------------------------
!  evaluate the orthonormalized polynomial
!  m,n mathematical numbers.(not fortran coefficients)
!  this routine is very slow
!---------------------------------------------------------
      integer :: m,n,am
      real(p8):: x
      real(p8):: fac
      am = abs(m)

      if(n.lt.am) then
	enpol=0
	return
      endif
      enpol = epol(am,n,x)
      fac = cfact(am,n)
      enpol = enpol * sqrt(fac) 
      return
      end function
!
!
!
      function epolv(m,n,x)
!----------------------------------------------------
!  evaluate the associate legendre polynomial 
!  (unnormalized) on the condition 
!        -1 <= x <= 1
!        n >= m >= 0
!  this routine is slow
!----------------------------------------------------
      integer,intent(in) :: m,n
      real(p8),intent(in),dimension(:) :: x
      real(p8),dimension(size(x)):: epolv
      real(p8),dimension(size(x)) :: y1,y2,y3
      integer:: nn,am
      am=abs(m)
      if(n.lt.am) then
	epolv=0
	return
      endif
      if(m.eq.0) then
	y1 = 1.0_p8
        y2 = x
      else
        y1 = fact1(am)*(1-x*x)**(0.5_p8*am)
        y2 = ((2*am+1)*x)*y1
      endif
      do nn=am+2,n
        y3 = 1.0_p8/(nn-am)*(((2*nn-1)*x)*y2 -(nn+am-1)*y1)
	y1 = y2
	y2 = y3
      enddo
      if(am.eq.n) then
	 epolv = y1
      else if(am.eq.n-1) then
	 epolv = y2
      else
         epolv = y3
      endif
      return
      end function
!
!
!
      real(p8) function epol(m,n,x)
!---------------------------------------------
!  evaluate the associate legendre polynomial 
!  (unnormalized) on the condition 
!        -1 <= x <= 1
!        n >= m >= 0
!  this routine is very slow
!---------------------------------------------
      integer :: m,n,am
      real(p8) :: x
      real(p8) :: y1,y2,y3
      integer:: nn
      am=abs(m)

      if(n.lt.am) then
	epol=0
	return
      endif

      if(am.eq.0) then
	y1 = 1.0_p8
        y2 = x
      else
        y1 = fact1(am)*(1-(x*x))**(0.5_p8*am)
        y2 = ((2*am+1)*x)*y1
      endif
      do 10 nn=am+2,n
        y3 = 1.0_p8/(nn-am)*(((2*nn-1)*x)*y2 - (nn+am-1)*y1)
	y1 = y2
	y2 = y3
   10 continue
      if(am.eq.n) then
	 epol = y1
      else if(am.eq.n-1) then
	 epol = y2
      else
         epol = y3
      endif
      return
      end function

!------------
  end module
!------------
