!-------------------------------------------------------------
 module bandmat
!-------------------------------------------------------------
!  SUBROUTINES:
!   lub(amat,nc) LU-decomposes amat(1:ni,1:nb) where amat is real
!                two dimensional banded array. no pivoting is performed.
!                ni is the size of the matrix. the size is assumed to be
!                ni*ni. nb is the bandwidth. amat(i,:) contains the i-th 
!                row of the matrix.  nc is the location of 
!                the center of the matrix so that amat(:,1:nc-1) is the 
!                subdiagonal part and amat(:,nc+1:nb) is the superdiagonal
!		 part of amat. the result returns in amat. 
!
!   solveb(amat,nc,b) solves linear equation (amat)x = b.  amat is the 
!                result of lub. b(:) is real or complex.
!
!   lub(amat,nc) if amat is three dimensional: amat(:,:,1:k), the third 
!                argument is interpreted as matrix number so that k matrices 
!                are LU-decomposed at once. nc must be the same for all k
!                matrices.
!
!   solveb(amat,nc,b) where amat is three dimensional and b is two 
!                dimensional: amat(:,:,1:k) and b(:,1:k) solves 
!                k linear equations (amat)x = b.
!
!   lub(amat,nc,msize) LU-decomposes multiple matrices of different sizes
!                at once. msize(i) specifies the size of amat(:,:,i).
!                msize(i) cannot increase as i becomes larger.
!                this multiple treatment of lub and solveb is efficient in
!                vector processors.
!
!   solveb(amat,nc,b,msize) solves k linear equations of different sizes,
!		 where amat(:,:,1:k), b(:,1:k).
!                msize(k) specifies the size of k-th matrix b(:,k).
!
!  FUNCTION:
!   banmul(amat,nc,b) multiplies a band matrix [amat,nc] to a matrix b.
!                amat(:,:), b(:,:). amat is real. b can be b(:).
!                b can be real or complex. the size and type of the resulting
!                matrix is the same as those of b.
!
!   banmul(b,amat,nc) multiplies a matrix b to a band matrix [amat,nc].
!                amat(:,:), b(:,:). amat is real.
!                b can be real or complex. the size and type of the resulting
!                matrix is the same as those of b.
!
!   banmul(amat,na,bmat,nb) 
!   banmul(amat,na,bmat,nb,nc) 
!                returns a banded matrix [amat,na]*[bmat,nb]. the band 
!                width of the resulting matrix is size(a,2)+size(b,2)-1.
!                optional output nc returns the index for the diagonal of the 
!                resulting matrix: nc=na+nb-1.
!  
!   mtrx(amat,nc) converts the [amat,nc] to a square matrix.
!
!-------------------------------------------------------------
 implicit none
 integer,parameter:: p8=selected_real_kind(p=11)
!----------------------------------------------------
 private
 public :: lub
 public :: solveb
 public :: banmul
 public :: mtrx
!----------------------------------------------------
 interface lub
   module procedure ludcmpb
   module procedure ludcmp1
   module procedure cludcmpb
   module procedure cludcmp1
 end interface

 interface solveb
   module procedure lubksbb
   module procedure lubksbbr
   module procedure lubksb1
   module procedure lubksb1r
   module procedure clubksbb
   module procedure clubksbbr
   module procedure clubksb1
   module procedure clubksb1r
 end interface

 interface banmul
   module procedure banmulr     ! real(p8)_band * real(p8)(:)
   module procedure banmulc     ! real(p8)_band * complex(p8)(:)
   module procedure banmulrm    ! real(p8)_band * real(p8)(:,:)
   module procedure banmulcm    ! real(p8)_band * complex(p8)(:,:)
   module procedure cbanmulr    ! complex(p8)_band * real(p8)(:)
   module procedure cbanmulc    ! complex(p8)_band * complex(p8)(:)
   module procedure cbanmulrm   ! complex(p8)_band * real(p8)(:,:)
   module procedure cbanmulcm   ! complex(p8)_band * complex(p8)(:,:)
   module procedure mr_ban_mul  ! real(p8)(:,:) * real(p8)_band 
   module procedure vr_ban_mul  ! real(p8)(:) * real(p8)_band 
   module procedure mc_ban_mul  ! complex(p8)(:,:) * real(p8)_band 
   module procedure vc_ban_mul  ! complex(p8)(:) * real(p8)_band 
   module procedure banbanmul   ! real(p8)_band * real(p8)_band
 end interface

 interface mtrx
   module procedure mtrxr
   module procedure mtrxc
 end interface
!
!
!
      contains
!     --------
!
      subroutine ludcmpb(amat,nc,msize,auxsize)
!    ----------------------------------------------------
!     LU decomposition of banded matrix by Gaussian elimination
!     number of superdiagonal and subdiagonal elements can be different.
!
!     amat(1:nd,1:ne,1:k) : input banded matrix
!             nd is the dimension of the first slot. 
!             ne is the band width of the matrix. 
!             k is the number of the matrices(this routine LU decomposes
!             multiple matrices at once)
!
!     nc : element number which corresponds to the diagonal element.
!
!     msize(1:k): (optional)the sizes of the matrices. the size must not be
!             increasing : msize(i) >= msize(i+1).
!             if msize is not present, then the matrix size is taken as nd.
!
!     auxsize(1:nd): (optional) auxiliary array which can be derived from
!             msize(see bandaux). auxsize can be specified to save
!             time to compute auxsize.
!
      implicit real(p8)(a-h,o-z)
      implicit integer(i-n)
      real(p8),dimension(:,:,:),intent(inout) :: amat
      integer,intent(in) :: nc
      integer,dimension(size(amat,3)),optional,intent(in) :: msize
      integer,dimension(size(amat,1)),optional :: auxsize
      integer,dimension(size(amat,3)):: nzchops
      integer,dimension(size(amat,1)+size(amat,2)):: ntchops
      integer :: ne,nd,ntchop,nzchop

      nd = size(amat,1)
      ne = size(amat,2)
      ntchop = size(amat,3)

      if(present(msize)) then
	nzchops = msize
	if(minval(msize) <1 ) then
	  print *,'ludcmpb: size of some matrices is zero.'
	  print *,'msize=',msize
	  stop
	endif
      else
	nzchops = size(amat,1)
      endif
      if(present(auxsize)) then
	ntchops = auxsize
      else
	call bandaux(nzchops,ntchops)
      endif

      nzchop = maxval(nzchops)

      do 10 nn=1,nzchop
      do 20 i=nc+1,ne
      do 20 mm=1,ntchops(nn+i-nc)
	amat(nn,i,mm) = amat(nn,i,mm)/amat(nn,nc,mm)
   20 continue
!
      do 30 i=nc+1,ne
      do 30 j=1,nc-1
      k  = max(i-nc,j)
      do 30 mm=1,ntchops(nn+k)
	amat(nn+j,i-j,mm)=amat(nn+j,i-j,mm) &
                       -amat(nn,i,mm)*amat(nn+j,nc-j,mm)
   30 continue
!
      do 40 j=1,nc-1
      do 40 mm=1,ntchops(nn+j)
	amat(nn+j,nc-j,mm) = amat(nn+j,nc-j,mm)/amat(nn,nc,mm)
   40 continue
   10 continue
!
      do 50 nn=1,nzchop
      do 50 mm=1,ntchops(nn)
	amat(nn,nc,mm) = 1.0/amat(nn,nc,mm)
   50 continue
!
      return
      end subroutine
!
!
!
      subroutine cludcmpb(amat,nc,msize,auxsize)
!    ----------------------------------------------------
!     LU decomposition of banded matrix by Gaussian elimination
!     number of superdiagonal and subdiagonal elements can be different.
!
!     amat(1:nd,1:ne,1:k) : input banded matrix
!             nd is the dimension of the first slot. 
!             ne is the band width of the matrix. 
!             k is the number of the matrices(this routine LU decomposes
!             multiple matrices at once)
!
!     nc : element number which corresponds to the diagonal element.
!
!     msize(1:k): (optional)the sizes of the matrices. the size must not be
!             increasing : msize(i) >= msize(i+1).
!             if msize is not present, then the matrix size is taken as nd.
!
!     auxsize(1:nd): (optional) auxiliary array which can be derived from
!             msize(see bandaux). auxsize can be specified to save
!             time to compute auxsize.
!
      implicit real(p8)(a-h,o-z)
      implicit integer(i-n)
      complex(p8),dimension(:,:,:),intent(inout) :: amat
      integer,intent(in) :: nc
      integer,dimension(size(amat,3)),optional,intent(in) :: msize
      integer,dimension(size(amat,1)),optional :: auxsize
      integer,dimension(size(amat,3)):: nzchops
      integer,dimension(size(amat,1)+size(amat,2)):: ntchops
      integer :: ne,nd,ntchop,nzchop

      nd = size(amat,1)
      ne = size(amat,2)
      ntchop = size(amat,3)

      if(present(msize)) then
	nzchops = msize
	if(minval(msize) <1 ) then
	  print *,'cludcmpb: size of some matrices is zero.'
	  print *,'msize=',msize
	  stop
	endif
      else
	nzchops = size(amat,1)
      endif
      if(present(auxsize)) then
	ntchops = auxsize
      else
	call bandaux(nzchops,ntchops)
      endif

      nzchop = maxval(nzchops)

      do 10 nn=1,nzchop
      do 20 i=nc+1,ne
      do 20 mm=1,ntchops(nn+i-nc)
	amat(nn,i,mm) = amat(nn,i,mm)/amat(nn,nc,mm)
   20 continue
!
      do 30 i=nc+1,ne
      do 30 j=1,nc-1
      k  = max(i-nc,j)
      do 30 mm=1,ntchops(nn+k)
	amat(nn+j,i-j,mm)=amat(nn+j,i-j,mm) &
                       -amat(nn,i,mm)*amat(nn+j,nc-j,mm)
   30 continue
!
      do 40 j=1,nc-1
      do 40 mm=1,ntchops(nn+j)
	amat(nn+j,nc-j,mm) = amat(nn+j,nc-j,mm)/amat(nn,nc,mm)
   40 continue
   10 continue
!
      do 50 nn=1,nzchop
      do 50 mm=1,ntchops(nn)
	amat(nn,nc,mm) = 1.0/amat(nn,nc,mm)
   50 continue
!
      return
      end subroutine

!
      subroutine lubksbb(amat,nc,a,msize,auxsize)
!    ----------------------------------------------------
!     LU back substitution(for complex array)
!
      implicit real(p8)(a-h,o-z)
      implicit integer(i-n)
      real(p8),dimension(:,:,:) :: amat
      integer :: nc
      complex(p8),dimension(:,:),intent(inout) :: a
      integer,dimension(size(amat,3)),optional :: msize
      integer,dimension(size(amat,1)),optional :: auxsize
      integer,dimension(size(amat,3)) :: nzchops
      integer,dimension(size(amat,1)+size(amat,2)) :: ntchops
      integer :: ne,nd,ntchop,nzchop
!
      if(present(msize)) then
	nzchops = msize
      else
	nzchops = size(a,1)
      endif
      if(present(auxsize)) then
	ntchops = auxsize
      else
	call bandaux(nzchops,ntchops)
      endif
!
      nd = size(amat,1)
      ne = size(amat,2)
      ntchop = size(amat,3)
      nzchop = maxval(nzchops)
!
      do 10 nn=2,nzchop
      nne = min(nn-1,nc-1)
      do 10 j=1,nne
      do 10 mm=1,ntchops(nn)
	a(nn,mm) = a(nn,mm)-a(nn-j,mm)*amat(nn,nc-j,mm)
   10 continue
!
      do 20 nn=1,nzchop
      do 20 mm=1,ntchops(nn)
	a(nn,mm) = a(nn,mm)*amat(nn,nc,mm)
   20 continue
!
      do 40 nn=nzchop-1,1,-1
      do 40 i=1,ne-nc
      do 40 mm=1,ntchops(nn+i)
	a(nn,mm) = a(nn,mm)-a(nn+i,mm)*amat(nn,nc+i,mm)
   40 continue
!
      return
      end subroutine
!
!
!
      subroutine clubksbb(amat,nc,a,msize,auxsize)
!    ----------------------------------------------------
!     LU back substitution(for complex array)
!
      implicit real(p8)(a-h,o-z)
      implicit integer(i-n)
      complex(p8),dimension(:,:,:) :: amat
      integer :: nc
      complex(p8),dimension(:,:),intent(inout) :: a
      integer,dimension(size(amat,3)),optional :: msize
      integer,dimension(size(amat,1)),optional :: auxsize
      integer,dimension(size(amat,3)) :: nzchops
      integer,dimension(size(amat,1)+size(amat,2)) :: ntchops
      integer :: ne,nd,ntchop,nzchop
!
      if(present(msize)) then
	nzchops = msize
      else
	nzchops = size(a,1)
      endif
      if(present(auxsize)) then
	ntchops = auxsize
      else
	call bandaux(nzchops,ntchops)
      endif
!
      nd = size(amat,1)
      ne = size(amat,2)
      ntchop = size(amat,3)
      nzchop = maxval(nzchops)
!
      do 10 nn=2,nzchop
      nne = min(nn-1,nc-1)
      do 10 j=1,nne
      do 10 mm=1,ntchops(nn)
	a(nn,mm) = a(nn,mm)-a(nn-j,mm)*amat(nn,nc-j,mm)
   10 continue
!
      do 20 nn=1,nzchop
      do 20 mm=1,ntchops(nn)
	a(nn,mm) = a(nn,mm)*amat(nn,nc,mm)
   20 continue
!
      do 40 nn=nzchop-1,1,-1
      do 40 i=1,ne-nc
      do 40 mm=1,ntchops(nn+i)
	a(nn,mm) = a(nn,mm)-a(nn+i,mm)*amat(nn,nc+i,mm)
   40 continue
!
      return
      end subroutine


      subroutine lubksbbr(amat,nc,a,msize,auxsize)
!    ----------------------------------------------------
!     LU back substitution(for real array)
!
      implicit real(p8)(a-h,o-z)
      implicit integer(i-n)
      real(p8),dimension(:,:,:) :: amat
      integer :: nc
      real(p8),dimension(:,:),intent(inout) :: a
      integer,dimension(size(amat,3)),optional :: msize
      integer,dimension(size(amat,1)),optional :: auxsize
      integer,dimension(size(amat,3)) :: nzchops
      integer,dimension(size(amat,1)+size(amat,2)) :: ntchops
      integer :: ne,nd,ntchop,nzchop
!
      if(present(msize)) then
	nzchops = msize
      else
	nzchops = size(a,1)
      endif
      if(present(auxsize)) then
	ntchops = auxsize
      else
	call bandaux(nzchops,ntchops)
      endif
!
      nd = size(amat,1)
      ne = size(amat,2)
      ntchop = size(amat,3)
      nzchop = maxval(nzchops)
!
      do 10 nn=2,nzchop
      nne = min(nn-1,nc-1)
      do 10 j=1,nne
      do 10 mm=1,ntchops(nn)
	a(nn,mm) = a(nn,mm)-a(nn-j,mm)*amat(nn,nc-j,mm)
   10 continue
!
      do 20 nn=1,nzchop
      do 20 mm=1,ntchops(nn)
	a(nn,mm) = a(nn,mm)*amat(nn,nc,mm)
   20 continue
!
      do 40 nn=nzchop-1,1,-1
      do 40 i=1,ne-nc
      do 40 mm=1,ntchops(nn+i)
	a(nn,mm) = a(nn,mm)-a(nn+i,mm)*amat(nn,nc+i,mm)
   40 continue
!
      return
      end subroutine
!
!
!
      subroutine clubksbbr(amat,nc,a,msize,auxsize)
!    ----------------------------------------------------
!     LU back substitution(for real array)
!
      implicit real(p8)(a-h,o-z)
      implicit integer(i-n)
      complex(p8),dimension(:,:,:) :: amat
      integer :: nc
      real(p8),dimension(:,:),intent(inout) :: a
      integer,dimension(size(amat,3)),optional :: msize
      integer,dimension(size(amat,1)),optional :: auxsize
      integer,dimension(size(amat,3)) :: nzchops
      integer,dimension(size(amat,1)+size(amat,2)) :: ntchops
      integer :: ne,nd,ntchop,nzchop
!
      if(present(msize)) then
	nzchops = msize
      else
	nzchops = size(a,1)
      endif
      if(present(auxsize)) then
	ntchops = auxsize
      else
	call bandaux(nzchops,ntchops)
      endif
!
      nd = size(amat,1)
      ne = size(amat,2)
      ntchop = size(amat,3)
      nzchop = maxval(nzchops)
!
      do 10 nn=2,nzchop
      nne = min(nn-1,nc-1)
      do 10 j=1,nne
      do 10 mm=1,ntchops(nn)
	a(nn,mm) = a(nn,mm)-a(nn-j,mm)*amat(nn,nc-j,mm)
   10 continue
!
      do 20 nn=1,nzchop
      do 20 mm=1,ntchops(nn)
	a(nn,mm) = a(nn,mm)*amat(nn,nc,mm)
   20 continue
!
      do 40 nn=nzchop-1,1,-1
      do 40 i=1,ne-nc
      do 40 mm=1,ntchops(nn+i)
	a(nn,mm) = a(nn,mm)-a(nn+i,mm)*amat(nn,nc+i,mm)
   40 continue
!
      return
      end subroutine
!
!

!
      subroutine ludcmp1(amat,nc)
!     ---------------------------
      real(p8),dimension(:,:),intent(inout) :: amat
      integer,intent(in) :: nc
      real(p8),dimension(size(amat,1),size(amat,2),1) :: b
      b(:,:,1) = amat
      call ludcmpb(b,nc)
      amat = b(:,:,1)
      return
      end subroutine
!
!
!
      subroutine cludcmp1(amat,nc)
!     ---------------------------
      complex(p8),dimension(:,:),intent(inout) :: amat
      integer,intent(in) :: nc
      complex(p8),dimension(size(amat,1),size(amat,2),1) :: b
      b(:,:,1) = amat
      call cludcmpb(b,nc)
      amat = b(:,:,1)
      return
      end subroutine
!
!
!
!
!
!
      subroutine lubksb1(amat,nc,a)
!     -----------------------------
      real(p8),dimension(:,:),intent(inout) :: amat
      integer,intent(in) :: nc
      complex(p8),dimension(:),intent(inout) :: a
      complex(p8),dimension(size(a,1),1):: b
      b(:,1)=a
      call lubksbb(reshape(amat,(/size(amat,1),size(amat,2),1/)),nc, b)
      a=b(:,1)
      end subroutine
!
!
!
      subroutine clubksb1(amat,nc,a)
!     -----------------------------
      complex(p8),dimension(:,:),intent(inout) :: amat
      integer,intent(in) :: nc
      complex(p8),dimension(:),intent(inout) :: a
      complex(p8),dimension(size(a,1),1):: b
      b(:,1)=a
      call clubksbb(reshape(amat,(/size(amat,1),size(amat,2),1/)),nc, b)
      a=b(:,1)
      end subroutine
!
!
!
      subroutine lubksb1r(amat,nc,a)
!     -----------------------------
      real(p8),dimension(:,:),intent(inout) :: amat
      integer,intent(in) :: nc
      real(p8),dimension(:),intent(inout) :: a
      real(p8),dimension(size(a,1),1):: b
      b(:,1)=a
      call lubksbbr(reshape(amat,(/size(amat,1),size(amat,2),1/)),nc, b)
      a=b(:,1)
      end subroutine
!
!
!
      subroutine clubksb1r(amat,nc,a)
!     -----------------------------
      complex(p8),dimension(:,:),intent(inout) :: amat
      integer,intent(in) :: nc
      real(p8),dimension(:),intent(inout) :: a
      real(p8),dimension(size(a,1),1):: b
      b(:,1)=a
      call clubksbbr(reshape(amat,(/size(amat,1),size(amat,2),1/)),nc, b)
      a=b(:,1)
      end subroutine
!
!
!
!
!
!
!
 subroutine bandaux(nzchop,ntchop)
 integer,dimension(:) :: nzchop
 integer,dimension(:) :: ntchop
 integer :: nz,nt,i

 nt = size(nzchop)
 nz = size(ntchop)
 if(maxval(nzchop(2:nt)-nzchop(1:nt-1)) > 0) then
   print *,'ntset: nzchop cannot increase.'
   print *,'nzchop:'
   print *, nzchop
   stop
 endif
 ntchop(nzchop(1):nz) = 0
 do i=1,nt-1
   ntchop(nzchop(i+1):nzchop(i)) = i
 enddo
 ntchop(1:nzchop(nt)) = nt

 end subroutine
!
!
!
      function mr_ban_mul(mat,band,nc)
      real(p8),dimension(:,:):: mat
      real(p8),dimension(:,:):: band
      integer,intent(in):: nc
      real(p8),dimension(size(mat,1),size(mat,2)):: mr_ban_mul
      integer :: i
      do i=1,size(mat,1)
	mr_ban_mul(i,:) = vr_ban_mul(mat(i,:),band,nc)
      enddo
      end function
!
!
!
      function mc_ban_mul(mat,band,nc)
      complex(p8),dimension(:,:):: mat
      real(p8),dimension(:,:):: band
      integer,intent(in):: nc
      complex(p8),dimension(size(mat,1),size(mat,2)):: mc_ban_mul
      integer :: i
      do i=1,size(mat,1)
	mc_ban_mul(i,:) = vc_ban_mul(mat(i,:),band,nc)
      enddo
      end function
!
!
!
      function vr_ban_mul(vec,a,nc)
!     --------------------------------
!     multiply a row vector and a banded matrix  y=vec.a
      real(p8),dimension(:),intent(in) :: vec
      real(p8),dimension(:,:),intent(in) :: a
      integer,intent(in) :: nc
      real(p8),dimension(size(vec)) :: vr_ban_mul
      integer :: i,is,c0,c1,b0,b1,n

      n=size(vec)
      if(n.gt.size(a,1)) then
	print *,'vr_ban_mul:sizeof bandmat small compared to vec'
        stop
      endif

      vr_ban_mul=0
      do i=1,size(a,2) 
	is = i-nc
	c0 = max(1,1-is)
	c1 = min(n,n-is)
	b0 = max(1,1+is)
	b1 = min(n,n+is)
	vr_ban_mul(b0:b1)=vr_ban_mul(b0:b1)+a(c0:c1,i)*vec(c0:c1)
      enddo

      end function
!
!
!
      function vc_ban_mul(vec,a,nc)
!     --------------------------------
!     multiply a row vector and a banded matrix  y=vec.a
      complex(p8),dimension(:),intent(in) :: vec
      real(p8),dimension(:,:),intent(in) :: a
      integer,intent(in) :: nc
      complex(p8),dimension(size(vec)) :: vc_ban_mul
      integer :: i,is,c0,c1,b0,b1,n

      n=size(vec)
      if(n.gt.size(a,1)) then
	print *,'vc_ban_mul:sizeof bandmat small compared to vec'
        stop
      endif

      vc_ban_mul=0
      do i=1,size(a,2) 
	is = i-nc
	c0 = max(1,1-is)
	c1 = min(n,n-is)
	b0 = max(1,1+is)
	b1 = min(n,n+is)
	vc_ban_mul(b0:b1)=vc_ban_mul(b0:b1)+a(c0:c1,i)*vec(c0:c1)
      enddo

      end function
!
!
!
 function banmulcm(a,nc,x)
!--------------------
!multiply a banded matrix and a matrix  y=a.x
 real(p8),dimension(:,:) :: a
 integer :: nc
 complex(p8),dimension(:,:) :: x
 complex(p8),dimension(size(x,1),size(x,2)) :: banmulcm
 integer :: k
 do k=1,size(x,2)
   banmulcm(:,k) = banmulc(a,nc,x(:,k))
 enddo
 end function
!
!
!
 function cbanmulcm(a,nc,x)
!--------------------
 complex(p8),dimension(:,:) :: a
 integer :: nc
 complex(p8),dimension(:,:) :: x
 complex(p8),dimension(size(x,1),size(x,2)) :: cbanmulcm
 integer :: k
 do k=1,size(x,2)
   cbanmulcm(:,k) = cbanmulc(a,nc,x(:,k))
 enddo
 end function
!
!
!
 function banmulrm(a,nc,x)
!--------------------
!multiply a banded matrix and a matrix  y=a.x
 real(p8),dimension(:,:) :: a
 integer :: nc
 real(p8),dimension(:,:) :: x
 real(p8),dimension(size(x,1),size(x,2)) :: banmulrm
 integer :: k
 do k=1,size(x,2)
   banmulrm(:,k) = banmulr(a,nc,x(:,k))
 enddo
 end function
!
!
!
 function cbanmulrm(a,nc,x)
!--------------------
 complex(p8),dimension(:,:) :: a
 integer :: nc
 real(p8),dimension(:,:) :: x
 real(p8),dimension(size(x,1),size(x,2)) :: cbanmulrm
 integer :: k
 do k=1,size(x,2)
   cbanmulrm(:,k) = cbanmulr(a,nc,x(:,k))
 enddo
 end function
!
!
!
      function banmulr(band,nc,x)
!     --------------------
!     multiply a banded matrix and a vector  y=a.x
      real(p8),dimension(:,:) :: band
      integer :: nc
      real(p8),dimension(:) :: x
      real(p8),dimension(size(x)) :: banmulr
      integer :: i,is,c0,c1,b0,b1,n

      n=size(x)
      if(n.gt.size(band,1)) then
	print *,'banmulr:sizeof bandmat small compared to x'
        stop
      endif

      banmulr=0
      do i=1,size(band,2) 
	is = i-nc
	c0 = max(1,1-is)
	c1 = min(n,n-is)
	b0 = max(1,1+is)
	b1 = min(n,n+is)
	banmulr(c0:c1)=banmulr(c0:c1)+band(c0:c1,i)*x(b0:b1)
      enddo

      end function
!
!
!
      function cbanmulr(band,nc,x)
!     --------------------
      complex(p8),dimension(:,:),intent(in) :: band
      integer :: nc
      real(p8),dimension(:) :: x
      real(p8),dimension(size(x)) :: cbanmulr
      integer :: i,is,c0,c1,b0,b1,n

      n=size(x)
      if(n.gt.size(band,1)) then
	print *,'cbanmulr:sizeof bandmat small compared to x'
        stop
      endif

      cbanmulr=0
      do i=1,size(band,2) 
	is = i-nc
	c0 = max(1,1-is)
	c1 = min(n,n-is)
	b0 = max(1,1+is)
	b1 = min(n,n+is)
	cbanmulr(c0:c1)=cbanmulr(c0:c1)+band(c0:c1,i)*x(b0:b1)
      enddo

      end function
!
!
!
      function banmulc(a,nc,x)
!     --------------------
!     multiply a banded matrix and a vector  y=a.x
      real(p8),dimension(:,:) :: a
      integer :: nc
      complex(p8),dimension(:) :: x
      complex(p8),dimension(size(x)) :: banmulc
      integer :: i,is,c0,c1,b0,b1,n
      real(p8):: dum1,dum2

      n=size(x)
      if(n.gt.size(a,1)) then
	print *,'banmulc:sizeof bandmat small compared to x'
        stop
      endif

      banmulc=a(1:n,nc)*x
      do i=1,size(a,2) 
	if(i.eq.nc) cycle
	is = i-nc
	c0 = max(1,1-is)
	c1 = min(n,n-is)
	b0 = max(1,1+is)
	b1 = min(n,n+is)
	banmulc(c0:c1)=banmulc(c0:c1)+a(c0:c1,i)*x(b0:b1)
      enddo

      end function
!
!
!
      function cbanmulc(a,nc,x)
!     --------------------
      complex(p8),dimension(:,:) :: a
      integer :: nc
      complex(p8),dimension(:) :: x
      complex(p8),dimension(size(x)) :: cbanmulc
      integer :: i,is,c0,c1,b0,b1,n

      n=size(x)
      if(n.gt.size(a,1)) then
	print *,'cbanmulc:sizeof bandmat small compared to x'
        stop
      endif

      cbanmulc=0
      do i=1,size(a,2) 
	is = i-nc
	c0 = max(1,1-is)
	c1 = min(n,n-is)
	b0 = max(1,1+is)
	b1 = min(n,n+is)
	cbanmulc(c0:c1)=cbanmulc(c0:c1)+a(c0:c1,i)*x(b0:b1)
      enddo

      end function
!
!
!
      function mtrxr(a,nc)
!     --------------------
      real(p8),dimension(:,:) :: a
      integer :: nc
      real(p8),dimension(size(a,1),size(a,1)) :: mtrxr
      integer :: i,j,n,jj

      mtrxr = 0
      n = size(a,1)

      do j=1,size(a,2)
       do i=1,n
	jj = j+i-nc
	if(jj.ge.1 .and. jj.le.n) mtrxr(i,jj) = a(i,j)
       enddo
      enddo
      end function
!
!
!
      function mtrxc(a,nc)
!     --------------------
      complex(p8),dimension(:,:) :: a
      integer :: nc
      complex(p8),dimension(size(a,1),size(a,1)) :: mtrxc
      integer :: i,j,n,jj

      mtrxc = 0
      n = size(a,1)

      do j=1,size(a,2)
       do i=1,n
	jj = j+i-nc
	if(jj.ge.1 .and. jj.le.n) mtrxc(i,jj) = a(i,j)
       enddo
      enddo
      end function
!
!
!
      function banbanmul(a,ac,b,bc,cc)
!     --------------------------------
!     multiply banded-banded matrix
      real(p8),dimension(:,:),intent(in) :: a
      real(p8),dimension(:,:),intent(in) :: b
      integer,intent(in) :: ac,bc
      integer,intent(out),optional :: cc
      real(p8),dimension(size(a,1),size(a,2)+size(b,2)-1) :: banbanmul
      integer :: i,j,k,aw,bw,cw,a0,a1,b0,b1,c0,c1,as,bs,nn,ii

      aw=size(a,2)
      bw=size(b,2)
      nn=size(a,1)
      cw=aw+bw-1
      if(present(cc)) then
        cc=ac+bc-1
      endif

      if(nn.ne.size(b,1)) then
	print *,'banbanmul: size mismatch.'
	print *,'sizeof a=',size(a,1)
	print *,'sizeof b=',size(b,1)
	stop
      endif

      banbanmul = 0
 
      do i=1,aw
        as = i-ac
        a0 = max(1,1-as)
        a1 = min(nn,nn-as)
        do j=1,bw
	  k = i+j-1
	  bs = j-bc
	  b0 = max(1,1-bs)
	  b1 = min(nn,nn-bs)
	  c0 = max(a0,b0-as)
	  c1 = min(a1,b1-as)
	  banbanmul(c0:c1,k)= banbanmul(c0:c1,k)+a(c0:c1,i)*b(c0+as:c1+as,j)
        enddo
      enddo
!
      end function

subroutine clean(a,ac)
real(p8),dimension(:,:):: a
integer:: ac,aw,nn,i
aw=size(a,2)
nn=size(a,1)
do i=1,ac-1
  a(i,:ac-i)=0
enddo
do i=1,aw-ac
  a(nn-i+1,ac+i)=0
enddo
end subroutine


end module

! 
!  program main
!  use bandmat
!  use lin_legendre
!  use eig
!  integer,parameter:: p8=selected_real_kind(p=11)
!  real(p8),dimension(40):: x
!  real(p8),dimension(1:40,3):: a
!  real(p8),dimension(0:41,3):: b
!  real(p8),dimension(40,5):: c1,c2
!  
!  a=1
!  b=1
!  a(2:11,2)=1.1_p8
!  b(2:13,2)=1.3_p8
!  b(:,2)=1.5_p8
!  b(:,3)=1.7_p8
!  c1=banmul(a,2,b(1:40,:),2)
!  c2=mul5_5(a,2,b,2)
! 
!  call spy(mtrx(c1-c2,2))
! 
!  
!  end program
