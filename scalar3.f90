!---------------
 module scalar3
!---------------
   use eig
   use lin_legendre
   implicit none
   integer, parameter :: p8=selected_real_kind(p=11)
   real(p8),parameter:: pi=3.141592653589793238462643_p8

!-----------------------------------------------------
!  routines
!-----------------------------------------------------
   private

   public:: leginit
      !  nzin,nthin,nxin
      !         the problem size. these 3 parameters are powers of 2.
      !  nzchopin,ntchopin,nxchopin
      !         chopping index. 
      !  zlenin
      !         period in the axial direction
      !  eelin
      !         map parameter
      !  mklinkin,mincin  (optional)
      !         mklink=1 or -1 if m and k is linked. otherwise 0
      !         plus or minus denotes spiral direction.
      !         minc is the number of symmetry in the azimuthal direction.
      !         if no symmetry, set minc=1.
      !         default is mklink=0, minc=1.

   public:: chopset
   public:: chopdo

   public:: horfft
   public:: verfft

   public:: rtran
      !  transforms. for example, rtran(a,1) transforms from 
      !  fff_space(function) to pff_space(physical).

   public:: toff,tofp

   public:: allocate, deallocate
      ! allocate/deallocate memory to type(scalar)
      ! allocate(a,space) allocates memory for the ???_space
      ! where ??? is either ppp,pfp,pff,fff. p,f denote physical, function,
      ! respectively.
      ! allocate(a) allocates memory for function space.

   public:: assignment(=)
      ! copy type(scalar)

   public:: calcat1,calcat0
      ! calcat1(a) returns the function value at infinity. for each k.
      ! calcat0(a) returns the function value at origin. for each k.

   public:: integ
   public:: integh
   public:: prodct

   public:: mload,msave

!----------------------------------------------------------
!  public values
!----------------------------------------------------------
   ! 1> dimensions
   integer,public:: nz,nth,nx,nzh
   integer,public:: nzchop,ntchop,nxchop,nxchopu
   real(p8),public:: zlen,eel,eel2
   integer,public:: mklink,minc
   integer,public:: ndimz,ndimth,ndimx
   integer,public:: nzchopdim,ntchopdim,nxchopdim
   integer,public,dimension(:),allocatable:: nzchops,ntchops
   integer,public,dimension(:),allocatable:: m
   real(p8),public,dimension(:,:),allocatable:: ak

   ! 2> types
   public:: scalar
   type scalar
     complex(p8),dimension(:,:,:),pointer:: e
     real(p8):: ln
     integer:: space
   end type

   ! 3> transform kit
   public:: transform,tfm
   type transform
     real(p8),dimension(:),pointer:: r
     !collocation points in r

     real(p8),dimension(:),pointer:: x
     !mapped collocation points in r

     real(p8),dimension(:),pointer:: w
     real(p8),dimension(:,:),pointer:: norm
     real(p8),dimension(:,:,:),pointer:: pf
     real(p8),dimension(:),pointer:: at0,at1
     real(p8),dimension(:),pointer:: ln
     real(p8),dimension(:),pointer:: z
     !collocation points in z

     real(p8),dimension(:),pointer:: thr,thi,th
     !collocation points in theta

     complex(p8),dimension(:),pointer:: ex
     integer:: nmax
   end type
   type(transform):: tfm

   ! 4> space tag
   integer,public,parameter:: ppp_space=0   ! full physical
   integer,public,parameter:: fff_space=1   ! full function
   integer,public,parameter:: pfp_space=2   
   integer,public,parameter:: pff_space=3

!============================================================
!  interfaces
!============================================================
   interface allocate
     module procedure salloc
   end interface

   interface deallocate
     module procedure sfree
   end interface

   interface mload
     module procedure mload0
   end interface

   interface msave
     module procedure msave0
   end interface

   interface assignment(=)
     module procedure copy0
   end interface

!============================================================

 contains
!--------
!
!
!
 subroutine leginit(nzin,nthin,nxin,nzchopin,ntchopin,nxchopin,&
		    zlenin,eelin,mklinkin,mincin)
!--------------------------------------------------------------------
   integer,intent(in):: nzin,nthin,nxin,nzchopin,ntchopin,nxchopin
   integer,intent(in),optional:: mklinkin,mincin
   integer:: mm,kk,kv,i
   real(p8):: zlenin,eelin

   nz = nzin
   nzh = nz/2
   nth= nthin
   nx = nxin
   nzchop = nzchopin
   ntchop = ntchopin
   nxchop = nxchopin
   zlen = zlenin
   eel  = eelin
   eel2 = eel**2
   if(present(mklinkin)) then
     mklink = mklinkin       ! spiral computation
   else
     mklink = 0
   endif
   if(present(mincin)) then
     minc = mincin           ! impose symmetry
     ! This is some sort of m increment
   else
     minc = 1
   endif


   ndimz  = nz+1
   ndimth = nth+1
   ndimx  = nx
   if(ndimth.eq.2) ndimth=1

   nzchopdim  = min(nzchop+3,nz)
   ntchopdim  = min(ntchop,nth)
   nxchopdim  = nxchop*2-1
   nxchopu    = nx-nxchop+2  ! chop location for conjugate side
   if (nxchop==1) nxchopu=1

   if ((2*nxchop>nx).and.((nx>1).or.(nxchop>1))) then
      print *, 'ERROR: scalar3() -- nxchop must be < nx/2,'
      print *, '                    unless nx=nxchop=1.'
      STOP
   end if

   allocate( m(ntchopdim) )
   m = (/ (minc*(mm-1),mm=1,ntchopdim) /)

   allocate(nzchops(nth))
   allocate(ntchops(nz))
   call chopset(0)

   allocate( ak(ntchopdim,nxchopdim) )
   do mm=1,ntchopdim
   do kk=1,nxchopdim
     kv=kk-1 + mklink*(mm-1)
     if(kk.gt.nxchop) then
       kv=-(nxchopdim-kk+1)
     endif
     ak(mm,kk) = 2*pi/zlen*kv
   enddo
   enddo

   !> tfm 
   allocate( tfm%r(nz) )
   allocate( tfm%w(nz) )
   allocate( tfm%x(nz) )
   call leg_zeros(tfm%x,tfm%w,tfm%r,eel)

   allocate( tfm%ln(nz) )
   tfm%ln = -log(1-tfm%x)

   allocate( tfm%norm(nzchopdim+14,ntchopdim) )
   tfm%norm = leg_norm(nzchopdim+14, m )

   allocate( tfm%pf(nzh,nzchopdim+1,ntchopdim) )
   tfm%pf = legtbl(tfm%x(:nzh),nzchopdim+1, m ,tfm%norm)

   allocate( tfm%at0(nzchopdim) )
   allocate( tfm%at1(nzchopdim) )
   tfm%at0 = legtbl( -1.0_p8,nzchopdim, 0 , tfm%norm(:,1))
   tfm%at1 = legtbl(  1.0_p8,nzchopdim, 0 , tfm%norm(:,1))

   tfm%nmax = max(nth,nx)*2
   allocate( tfm%ex(tfm%nmax) )
   call exset( tfm%ex, tfm%nmax )

   if(mklink.eq.0) then
     allocate( tfm%z(nx+1) )
     tfm%z = (/ (zlen*(kk-1)/nx, kk=1,nx+1) /)
   else if(mklink.eq.1 .or. mklink.eq.-1) then
     allocate( tfm%z(2*nth+1) )
     tfm%z = (/ (zlen*(kk-1)/(2*nth), kk=1,2*nth+1) /)
   else
     print *,'leginit: invalid mklink. mklink=',mklink
     stop
   endif

   allocate( tfm%thr(minc*nth+1) )
   allocate( tfm%thi(minc*nth+1) )
   tfm%thr = (/ (2*pi/minc/nth*(i-1), i=1,minc*nth+1) /)
   tfm%thi = (/ (2*pi/minc/nth*(i-0.5_p8), i=1,minc*nth+1) /)

   allocate( tfm%th(2*minc*nth+1) )
   tfm%th = (/ (pi/minc/nth*(i-1), i=1,2*minc*nth+1) /)

   end subroutine



   subroutine chopset(iof)
!-----------------------------------------------------------------
!  iof : offset of nzchop
!          +----------------------------------------+
!  nzchop->+--+                                     +
!  (nzdeg) +  +-----+              inactive         +
!          +        +-----+                         +
!          +              +-----+                   +
!          +                    +-----+ grad=1      +
!          +                          +-----+       +
!          +                                +-----+ +
!          +         active coefficients          +-+
!          +                                        +
!          +                                        +
!          +----------------------------------------+
!                                                   ^ ntchop
!-----------------------------------------------------------------
   integer:: iof
   integer:: mm,nn,nzdeg
   nzchop = nzchop + iof
   nzdeg  = nzchop 
   if(nzchop.gt.nzchopdim) then
     print *,'chopset: nzchop too large.'
     print *,'  nzchop=',nzchop,'  nzchopdim=',nzchopdim
     stop
   endif

   do mm=1,ntchopdim
     nzchops(mm) = max(min(nzchop,nzdeg-m(mm)),0)
   enddo
   nzchops(ntchopdim+1:)=0
   do nn=1,nz
     ntchops(nn) = max(min(ntchop,(nzdeg-nn+minc)/minc),0)
   enddo
   return
   end subroutine



   subroutine chopdo(a)
!-----------------------------------------------------------------
   type(scalar):: a
   integer::mm

   if(a%space .eq. fff_space) then
     do mm=1,ntchop
       a%e(nzchops(mm)+1:,mm,:)=0
     enddo
     a%e(:,ntchop+1:,:) = 0
   else 
!     a%e(:,:,nx+1:)=0
!     a%e(:,nth+1:,:nx)=0
!     a%e(nz+1:,:nth,:nx)=0
     a%e(:,nth+1,:) = 0._p8
     a%e(nz+1,:,:)  = 0._p8
   endif

   end subroutine



   subroutine salloc(a,sp)
!---------------------------------------------------------
! allocate array. sp: ppp_space or fff_space
!---------------------------------------------------------
   type(scalar):: a
   integer,optional:: sp
   integer:: mm, err

   if(present(sp)) then
     a%space = sp
   else
     a%space = fff_space
   endif

   if(associated(a%e)) then 
     print *,'salloc: already allocated.'
     stop
   endif
   if(a%space .eq. fff_space) then
     allocate( a%e(nzchopdim,ntchopdim,nxchopdim),STAT=err )
   else
     allocate( a%e(ndimz,ndimth,ndimx),STAT=err )
   endif
   if (err.ne.0) then
	print *, "ERROR: salloc() -- Could not allocate array."
	STOP
   end if

   return
   end subroutine




   subroutine sfree(a)
!---------------------------------------------------------
! free array.
!---------------------------------------------------------
   type(scalar):: a

   if(.not.associated( a%e )) then
     print *,'sfree:tried to deallocate though not associated.'
     stop
   endif
  
   deallocate( a%e)

   return
   end subroutine



 subroutine msave0(a,fn)
!---------------------------------------------------------
 type(scalar):: a
 character(len=*),intent(in):: fn
 call msavex(a,fn)
 print *,'written to ',trim(fn)
 end subroutine



 subroutine msavex(a,fn)
!---------------------------------------------------------
 type(scalar):: a
 character(len=*):: fn
 integer:: status,mm,kk,i

 open(unit=7,file=fn,status='unknown',form='unformatted',iostat=status)
 if(status.ne.0) then
   print *,'msave0_in_scalar3: failed to open ',trim(fn)
   stop
 endif

 write(7) nz,nth,nx
 write(7) nzchop,ntchop,nxchop
 write(7) a%space,a%ln
 write(7) zlen,eel
 write(7) minc,mklink
 if(a%space.eq.fff_space) then
   do mm=1,ntchop
     write(7) a%e(:nzchops(mm),mm,:nxchopdim)
   enddo
 else
   write(7) a%e(:nz,:nth,:nx)
 endif
 close(7)

 return
 end subroutine


 subroutine mload0(fn,a)
!---------------------------------------------------------
 character(len=*):: fn
 type(scalar):: a
 call mloadx(fn,a)
 print *,'read from ',trim(fn)
 end subroutine



 subroutine mloadx(fn,a)
!---------------------------------------------------------
 character(len=*):: fn
 type(scalar):: a
 integer:: status,mm
 integer:: iz,ith,ix
 integer:: izchop,itchop,ixchop
 real(p8):: izlen,ieel
 integer:: iinc,iklink,iof

 open(unit=7,file=fn,status='old',form='unformatted',iostat=status)
 if(status.ne.0) then
   print *,'msave0_in_scalar3: failed to open ',trim(fn)
   stop
 endif

 read(7) iz,ith,ix
 read(7) izchop,itchop,ixchop
 read(7) a%space,a%ln
 read(7) izlen,ieel
 read(7) iinc,iklink
 if(a%space.eq.fff_space) then
   iof = izchop-nzchop
   call chopset(iof)
   a%e=0
   do mm=1,itchop
     read(7) a%e(:nzchops(mm),mm,:nxchopdim)
   enddo
   call chopset(-iof)
   call chopdo(a)
 else
   read(7) a%e(:iz,:ith,:ix)
 endif
 close(7)

 if(izlen.ne.zlen) then
   print *,'mloadx: zlen inconsistent.'
   print *,'prog =',zlen
   print *,'file =',izlen
 endif

 if(ieel.ne.eel) then
   print *,'mloadx: eel inconsistent.'
   print *,'prog =',eel
   print *,'file =',ieel
 endif

 if(minc.ne.iinc) then
   print *,'mloadx: minc inconsistent.'
   print *,'prog =',minc
   print *,'file =',iinc
 endif

 if(mklink.ne.iklink) then
   print *,'mloadx: mklink inconsistent.'
   print *,'prog =',mklink
   print *,'file =',iklink
 endif

 if(nzchop.ne.izchop .or. ntchop.ne.itchop .or. nxchop.ne.ixchop) then
   print *,'mloadx: chopping data inconsistent.'
   print *,'prog nz,nt,nxchop=',nzchop,ntchop,nxchop
   print *,'file nz,nt,nxchop=',izchop,itchop,ixchop
 endif

 return
 end subroutine

  subroutine copy0(b,a)
  type(scalar),intent(inout):: b
  type(scalar),intent(in):: a
  integer:: mm

  if(.not.associated(b%e)) then
    print *,'copy0: no memory for lhs.'
    stop
  endif
  if(size(b%e)<size(a%e)) then
    call deallocate(b)
    call allocate(b,a%space)
  endif
 
  b%space = a%space
  b%ln = a%ln
  if(a%space .eq. fff_space) then
    do mm=1,ntchop
      b%e(:nzchops(mm),mm,:) = a%e(:nzchops(mm),mm,:)
    enddo
    call chopdo(b)
  else
    b%e(:size(a%e,1),:size(a%e,2),:size(a%e,3)) = a%e
    b%e(:,:,size(a%e,3)+1:)=0
    b%e(:,size(a%e,2)+1:,:size(a%e,3))=0
    b%e(size(a%e,1)+1:,:size(a%e,2),:size(a%e,3))=0
  endif
  end subroutine


! THE FOLLOWING IS AN OLDER VERSION OF HORFFT()
! rcs(), csr(), and log2() are assembly-language routines on the Cray C90
! fortran versions are available, but they are relatively slow
   subroutine horfft(a,is)
!----------------------------------------------------
!  call in r physicl space
!  expected calling sequence:
!  function space  ->  physical space
!    call rtran(a,1)
!    call verfft(a,1)
!    call horfft(a,1)
!  physical space  ->  function space
!    call horfft(a,-1)
!    call verfft(a,-1)
!    call rtran(a,-1)
!-------------------------------------------------------
   common /ft/ npts,nskip,mtrn,mskip,isign,log,iex
   integer::   npts,nskip,mtrn,mskip,isign,log,iex
!-------------------------------------------------------
   type(scalar):: a
   integer,intent(in):: is
   integer,parameter:: p=1
   complex(p8),dimension(:,:),allocatable:: b
   integer:: mm,kk,k1,k2
   real(p8):: fac
   integer,external:: log2

   allocate( b(nth+p,nz*nx) )
 
   if(is .eq.-1) then
      !> physical to foruier space
      if(a%space.ne.ppp_space)then
	print *,'horfft: not in physical space'
      endif
      if(nth.eq.1) then
	a%space=pfp_space
	return
      endif
      npts = 2*nth
      nskip = 1
      mtrn = nz*nx
      mskip = nth+p
      isign = is
      log = log2(npts)
      iex = tfm%nmax/npts
      do mm=1,nth
      do kk=1,nx
	k1=nz*(kk-1)+1
	k2=nz*kk
 	b(mm,k1:k2) = a%e(:nz,mm,kk)
      enddo
      enddo

      call rcs(b,tfm%ex)
      !Real to complex symmetric transform
      !Keeps only half of the symmetric array
      !Packs b(0) and b(n/2) into b(0) since they are both real.
      
      fac = 1.0_p8/(4*nth)
      do mm=1,ntchop
      do kk=1,nx
	k1=nz*(kk-1)+1
	k2=nz*kk
        a%e(:nz,mm,kk) = b(mm,k1:k2)*fac
      enddo
      enddo
      a%e(:,ntchop+1:,:) = 0.0

      call norm(a)

      a%space = pfp_space
 
   else
      !> fourier to physical
      if(a%space.ne.pfp_space)then
	print *,'horfft: not in pfp space'
      endif
      if(nth.eq.1) then
	a%space=ppp_space
	return
      endif
      call norm(a)
      npts = nth*2
      nskip = 1
      mtrn = nz*nx
      mskip = nth+p
      isign = is
      log = log2(npts)
      iex = tfm%nmax/npts

      do mm=1,ntchop
      do kk=1,nx
	k1=nz*(kk-1)+1
	k2=nz*kk
 	b(mm,k1:k2) = a%e(:nz,mm,kk)
      enddo
      enddo
      b(ntchop+1:,:) = 0.0
 
      call csr(b,tfm%ex)
 
      do mm=1,nth
      do kk=1,nx
	k1=nz*(kk-1)+1
	k2=nz*kk
        a%e(:nz,mm,kk) = b(mm,k1:k2)
      enddo
      enddo

      a%space = ppp_space
 
   endif
   deallocate( b )
   return
   end subroutine
 
      subroutine norm(q)
!    ----------------------------------------------------
      type(scalar):: q
      q%e(:,1,:) = real(q%e(:,1,:))
      return
      end subroutine

! THE FOLLOWING IS AN OLDER VERSION OF VERFFT()
! ftrvmt() and log2() are assembly-language routines on the Cray C90
! fortran versions are available, but are relatively slow
      subroutine verfft(a,is)
!----------------------------------------------------
!  call in r physicl space
!  expected calling sequence:
!  function space  ->  physical space
!    call rtran(a,1)
!    call verfft(a,1)
!    call horfft(a,1)
!  physical space  ->  function space
!    call horfft(a,-1)
!    call verfft(a,-1)
!    call rtran(a,-1)
!----------------------------------------------------
   common /ft/ npts,nskip,mtrn,mskip,isign,log,iex
   integer::   npts,nskip,mtrn,mskip,isign,log,iex
!----------------------------------------------------
   type(scalar):: a
   integer,intent(in):: is
   complex(p8),dimension(:,:),allocatable:: b
   integer,external:: log2
   integer:: kk,mm,m1,m2
   integer,parameter:: p=1
   real(p8):: fac

   allocate( b(nz*ntchop+p,nx) )
 
   if(is .eq.-1) then
      !> physical to foruier space
      if(a%space.ne.pfp_space) then
	print *,'verfft: not in pfp space'
      endif
      if(nx.eq.1) then
	a%space=pff_space
	return
      endif
      npts = nx
      nskip = nz*ntchop+p
      mtrn = nz*ntchop
      mskip = 1
      isign = is
      log = log2(npts)
      iex = tfm%nmax/npts

!     do kk=1,nx
!       b(:nz*ntchop,kk) = (/ a%e(:nz,:ntchop,kk) /)
!     enddo
      do kk=1,nx
      do mm=1,ntchop
	m1=nz*(mm-1)+1
	m2=nz*mm
        b(m1:m2,kk) = a%e(:nz,mm,kk)
      enddo
      enddo

      call ftrvmt(b,tfm%ex)
      !complex to complex FFT
      
      fac = 1./float(nx)
!     do kk=1,nx
!       a%e(:nz,:ntchop,kk) = reshape(b(:nz*ntchop,kk) *fac, (/nz,ntchop/))
!     enddo
      do kk=1,nx
      do mm=1,ntchop
	m1=nz*(mm-1)+1
	m2=nz*mm
        a%e(:nz,mm,kk) = b(m1:m2,kk)*fac
      enddo
      enddo

      do kk = nxchop+1,nxchopu-1
        a%e(:nz,:ntchop,nxchop+1:nxchopu-1) = 0.0
      enddo

      a%space = pff_space
 
 
  else          
      !> fourier to physical
      if(a%space.ne.pff_space) then
	print *,'verfft: not in pff space'
      endif
      if(nx.eq.1) then
	a%space=pfp_space
	return
      endif
      npts = nx
      nskip = nz*ntchop+p
      mtrn = nz*ntchop
      mskip = 1
      isign = is
      log = log2(npts)
      iex = tfm%nmax/npts
 
      a%e(:nz,:ntchop,nxchop+1:nxchopu-1) = 0
!     do kk=1,nx
!       b(:nz*ntchop,kk) = (/ a%e(:nz,:ntchop,kk) /)
!     enddo
      do kk=1,nx
      do mm=1,ntchop
	m1=nz*(mm-1)+1
	m2=nz*mm
        b(m1:m2,kk) = a%e(:nz,mm,kk)
      enddo
      enddo
 
      call ftrvmt(b,tfm%ex)
 
!     do kk=1,nx
!       a%e(:nz,:ntchop,kk) = reshape(b(:nz*ntchop,kk) ,(/nz,ntchop/) )
!     enddo
      do kk=1,nx
      do mm=1,ntchop
	m1=nz*(mm-1)+1
	m2=nz*mm
        a%e(:nz,mm,kk) = b(m1:m2,kk)
      enddo
      enddo

      a%space = pfp_space
   endif
 
   deallocate( b )
 
   return
   end subroutine

   subroutine rtran(a,is)
!----------------------------------------------------
   type(scalar):: a
   type(scalar):: b
   integer,intent(in):: is
   complex(p8),dimension(:,:),allocatable:: be,bo
   integer:: nn,i,mm
 
   if(is.eq.1) then
     !> polynomial to physical
     if(a%space .ne. fff_space) then
       print *,'rtran: not in fff_space'
     endif
     call salloc(b,ppp_space)
     call chopdo(b)
     b%ln=0

     allocate( be(nzh,nxchopdim) )
     allocate( bo(nzh,nxchopdim) )

     do mm=1,ntchop
	nn = nzchops(mm)
	if(nn.ge.1) then
	  be= tfm%pf(:nzh,1:nn:2,mm) .mul. a%e(1:nn:2,mm,:nxchopdim)
	else
	  be=0
	endif
	if(nn.ge.2) then
	  bo= tfm%pf(:nzh,2:nn:2,mm) .mul. a%e(2:nn:2,mm,:nxchopdim)
	else
	  bo=0
	endif
	b%e(1:nzh,mm,1:nxchop)= be(:,:nxchop)+bo(:,:nxchop)
	
	!-----------change below: was hanging on 2-d, where nxchop=1  ------------
	!------------ code remains same if 3-d, i.e. if nxchop > 1 ---------
	!-----change:  if statement added:  code inside "if" unchanged------
	if(nxchop .ne. 1) then
           b%e(1:nzh,mm,nxchopu:nx)= be(:,nxchop+1:)+bo(:,nxchop+1:)
           b%e(nz:nzh+1:-1,mm,nxchopu:nx)= be(:,nxchop+1:)-bo(:,nxchop+1:)
    endif

	b%e(nz:nzh+1:-1,mm,:nxchop)= be(:,:nxchop)-bo(:,:nxchop)
	
     enddo
     b%e(:,ntchop+1:,:)=0
     b%e(:,:ntchop,nxchop+1:nxchopu-1)=0
 
     if(a%ln.ne.0.0) then
        print *,'rtran:to physical space though logterm is nonzero.'
        print *,'logterm=',a%ln
        b%e(1:nz,1,1)=b%e(1:nz,1,1)+a%ln*tfm%ln(1:nz)
     endif

     call sfree(a)
     a%e => b%e
     nullify(b%e)
     a%space = pff_space

     deallocate( be,bo )
 
  else
     !> physical to polynomial
     if(a%space .ne. pff_space) then
       print *,'rtran: not in pff_space.'
     endif
     call salloc(b,fff_space)
     call chopdo(b)

     allocate( be(nzh,nxchopdim) )
     allocate( bo(nzh,nxchopdim) )

     do mm=1,ntchop
       nn = nzchops(mm)
       do i=1,nzh
	 be(i,:nxchop) = (a%e(i,mm,:nxchop)+a%e(nz-i+1,mm,:nxchop))*tfm%w(i)
	 bo(i,:nxchop) = (a%e(i,mm,:nxchop)-a%e(nz-i+1,mm,:nxchop))*tfm%w(i)

!-----------change below: was hanging on 2-d, where nxchop=1  ------------
!------------ code remains same if 3-d, i.e. if nxchop > 1 ---------
!-----change:  if statement added:  code inside "if" unchanged------
         if(nxchop .ne. 1) then
            be(i,nxchop+1:) = (a%e(i,mm,nxchopu:nx)+a%e(nz-i+1,mm,nxchopu:nx)) &
                 *tfm%w(i)	
            bo(i,nxchop+1:) = (a%e(i,mm,nxchopu:nx)-a%e(nz-i+1,mm,nxchopu:nx)) &
                 *tfm%w(i)
         endif
         
      enddo
      if(nn.ge.1) b%e(1:nn:2,mm,:)=transpose(tfm%pf(:nzh,1:nn:2,mm)).mul. be
      if(nn.ge.2) b%e(2:nn:2,mm,:)=transpose(tfm%pf(:nzh,2:nn:2,mm)).mul. bo
     enddo

     call sfree(a)
     a%e => b%e
     nullify(b%e)
     a%space=fff_space

     deallocate( be, bo )

  endif
!
  return
end subroutine rtran


   subroutine toff(a)
   type(scalar):: a
   if (a%space.eq.ppp_space) call horfft(a,-1)
   if (a%space.eq.pfp_space) call verfft(a,-1)
   if (a%space.eq.pff_space) call rtran(a,-1)
   end subroutine
  
  
   subroutine tofp(a)
   type(scalar):: a
   if (a%space.eq.fff_space) call rtran(a,1)
   if (a%space.eq.pff_space) call verfft(a,1)
   if (a%space.eq.pfp_space) call horfft(a,1)
   end subroutine


   function calcat0(a)
!  -------------------
   type(scalar),intent(in):: a
   complex(p8),dimension(size(a%e,3)):: calcat0
   calcat0 = transpose(a%e(:nzchop,1,:)).mul. tfm%at0(:nzchop)
   end function



   function calcat1(a)
!  -------------------
   type(scalar),intent(in):: a
   complex(p8),dimension(size(a%e,3)):: calcat1
   calcat1 = transpose(a%e(:nzchop,1,:)).mul. tfm%at1(:nzchop)
   end function
!
!
!
   real(p8) function integ(f)
!----------------------------------------------------
!  integration of a function over a domain 
!      0  <  r  < \infty
!      0  <  phi < 2*pi
!      0  <  z   < zlen
!  call in f/f space
!  integrated function g is in the form f=g/(1-mu)^2
!----------------------------------------------------
   type(scalar):: f

   if(f%space.ne.fff_space) then
      print *,'integ: not in fff_space'
      stop
   endif
   if(f%ln.ne.0.0) then
      print *,'integ: logterm not zero'
      print *,'logterm=',f%ln
   endif
 
   integ = 4*pi*zlen*eel2*real(f%e(1,1,1))*tfm%norm(1,1)

   return
   end function
!
!
!
   real(p8) function prodct(a,b)
!----------------------------------------------------
!  calculate the product and integrate over the domain
!  call in r-physical / phi,z-fourier space
!  what's integrated is f=a*b*(1-mu)^2
!----------------------------------------------------
   type(scalar):: a,b
   complex(p8):: prod(nz)
   integer:: mm,kk
       
   if(a%space.ne.pff_space .or. b%space.ne.pff_space) then
     print *,'prodct:not in pff_space'
     print *,'a%space,b%space=',a%space,b%space
     stop
   endif
   if(a%ln.ne.0.0 .or. b%ln.ne.0.0) then
     print *,'prodct:logterm not zero'
   endif

   prod=0
   do kk = 1,nx
     if(kk.gt.nxchop .and. kk.lt.nxchopu) cycle
     prod = prod+a%e(:nz,1,kk)*conjg(b%e(:nz,1,kk))
   enddo

   do kk = 1,nx
     if(kk.gt.nxchop .and. kk.lt.nxchopu) cycle
     do mm = 2,ntchop
        prod = prod + 2*(real(a%e(:nz,mm,kk))*real(b%e(:nz,mm,kk)) &
                      + aimag(a%e(:nz,mm,kk))*aimag(b%e(:nz,mm,kk)))
     enddo
   enddo

   prodct = sum((prod*tfm%w)*tfm%pf(1,1,1))
   prodct = 4*pi*zlen*eel2*prodct*tfm%norm(1,1)

   return
   end function
!
!
!
   real(p8) function integh(f,ain,bin)
!----------------------------------------------------
!  integration of a function over a half domain 
!  f is in hat form (i.e. g=(1-mu)^2*f where g is the function
!  you really want to integrate.
!       0  <  r  < \infty
!       a  <  phi < b
!       0  <  z   < zlen
!  call in p/f space
!----------------------------------------------------
   type(scalar):: f
   real(p8),intent(in),optional:: ain, bin
   complex(p8):: fab
   complex(p8),parameter:: iu = (0.0_p8,1.0_p8)
   real(p8):: a,b
   integer:: mm,m0
   if(f%space.ne.pff_space) then
     print *,'integh: not in pff_space'
     stop
   endif
   if(present(ain)) then
     a=ain
   else
     a=-pi/2        ! default: right hand side half plane
   endif
   if(present(bin)) then
     b=bin
   else
     b=pi/2
   endif

   integh = (b-a)*sum(tfm%w*f%e(:nz,1,1))
   do mm=2,ntchop
     m0=m(mm)
     fab = (exp(iu*m0*b)-exp(iu*m0*a))/(iu*m0)
     integh = integh+2*real(fab*sum(tfm%w*f%e(:nz,mm,1)))
   enddo
   integh = integh*zlen*eel2

   return
   end function

end module
