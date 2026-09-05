!========================================================
!     USE EIG  makes the interface explicit and allows
!     the use of the generic names
!========================================================
!   Operators:
!     .mul.   matrix/vector multiplication
!     .lr.    left-right matrix paste
!     .ud.    up-down matrix paste
!
!   Functions:
!     geneig(a) returns an one-dimensional complex array 
!             containing the eigenvalues of matrix a(:,:).
!             matrix a must be square and complex.
!     geneig(a,b) 
!             solves generalized eigenvalue problem a.x=(lamb)b.x
!
!     eigvec(a,eig) returns a matrix whose columns are the 
!             eigenvectors corresponding to the eigenvalues eig.
!             a(:,:) and eig(:) are complex arrays. number of 
!             columns in the result is the same as the number of 
!             elements in eig. eigenvectors are normalized such that
!             its norm is one sum(conjg(eivec)*eivec)=1.
!     eigvec(a,b,eig) is for generalized eigenvalue problem.
!     eigvec(a,eig,adj), eigvec(a,b,eig,adj) in addition returns the 
!             adjoint eigenvectors in columns of complex(p8) adj(:,:).
!             adjoint is normalized such that sum(adj*eivec) = 1
!             (not sum(conjg(adj)*eivec), nor sum(conjg(adj)*a*adj).
!
!     inv(a)  returns the inverse matrix of a. 
!             matrix a must be square and can be either 
!             complex or real.  
!
!     zeros(ni,nj)  returns a real ni*nj matrix whose elements are zero.
!             ni,nj are integers.
!
!     eye(ni) returns a real identity matrix of order ni. ni is integer.
!
!     diag(X) if X is one dimensional, diag(X) is a two dimensional
!             array whose diagonals are the elements of X. if X is two
!             dimensional, diag(X) is one dimensional whose elements are
!             the diagonals of X.
!
!     fliplr(X) flips matrix X in the left/right direction.
!
!     flipud(X) flips matrix X in the up/down direction.
!
!  Subroutines
!     lu(a,w) decomposes real(p8)/complex(p8) a(:,:) into lu. 
!             integer w(n,2) is a work area with n=size(a,1)
!
!     solve(a,w,b) solves a.x=b where a(:,:), w(:,:) is the result of lu 
!             above, real(p8)/complex(p8) b(:,:). On output, b is replaced
!             by x.
!
!     spy(a)  shows the structure of the matrix a(:,:) based on the values
!             of abs(a). matrix a can be complex or real. 
!
!     mcat(a) prints the values of matrix a(:,:). if a is one dimensional, 
!             values are printed out in a greater(~16 digits) accuracy.
!             matrix a can be complex or real.
!     
!     msave(a,'aa') saves matrix a to disk under the name 'aa'.
!             matrix a can be complex or real and can be one/two dimensional.
!             if a is complex valued, the number of rows are doubled. the 
!             format is
!                real(a(1,1)) imag(a(1,1)) real(a(1,2)) ....
!                real(a(2,1)) imag(a(2,1)) real(a(2,2)) ....
!                                   .
!                                   .
!
!             The file can be read by mload below. 
!             Also, MATLAB can read the file as follows. 
!
!             >> load aa
!
!             If the matrix is complex, a m-file
!
!             function [x]=mloadc(a)
!             %  MLOADC   convert ascii file to complex matrix
!             %         example:   
!             %         >> load abc
!             %         >> abc=mloadc(abc);
!             %
!             %         see also MSAVEC, fortran routines msave, mload
!             %
!             [ni,nj]=size(a);
!             x=a(:,1:2:nj-1)+1i*a(:,2:2:nj);
!          
!             can recover the matrix.
!
!     mload('aa',a) loads the matrix a(:,:). matrix a can be complex or real.
!             the size of a must agree with the size of the file
!             to be read. mload can load a matrix saved by MATLAB by
!
!             >> save aa aa -ascii -double
!
!             To save a complex matrix, following preprocessing by a m-file
!             is useful.
!
!             function [x]=msavec(a)
!             %  MSAVEC   convert a complex matrix to a real matrix 
!             %         for ascii save.
!             %         example:   
!             %         >> a=msavec(abc);
!             %         >> save a a -ascii -double
!             %
!             %         see also MLOADC, fortran routines msave, mload
!             %
!             [ni,nj]=size(a);
!             x=zeros(ni,nj*2);
!             for j=1:nj
!               x(:,j*2-1)=real(a(:,j));
!               x(:,j*2  )=imag(a(:,j));
!             end
!
!     mload('aa',a,ni,nj) loads a real or complex matrix of
!             unknown size. outputs ni,nj are the size of the loaded matrix.
!             array a must be large enough to store the ni*nj matrix.
!
!     --------------------------------------------------
      module eig
      implicit none
!     --------------------------------------------------
      private
      public:: operator(.mul.)
      public:: operator(.lr.)
      public:: operator(.ud.)
      public:: inv
      public:: geneig
      public:: eigvec
      public:: zeros,eye
      public:: diag
      public:: flipud
      public:: fliplr
      public:: spy, mcat
      public:: mload, msave
      public:: lu, solve
!     --------------------------------------------------
      integer, parameter :: p8=selected_real_kind(p=11)

!     -------------------------
      interface operator(.mul.)
!     -------------------------
      module procedure mulrr
      module procedure mulcc
      module procedure mulcr
      module procedure mulrc
      module procedure mulrr1
      module procedure mulrc1
      module procedure mulcr1
      module procedure mulcc1
      module procedure mulc1r
      end interface

!     -------------------------
      interface operator(.lr.)
!     -------------------------
      module procedure paste_lr_rr
      module procedure paste_lr_cc
      end interface

!     -------------------------
      interface operator(.ud.)
!     -------------------------
      module procedure paste_ud_rr
      module procedure paste_ud_cc
      end interface

!     ----------------
      interface inv
!     ----------------
      module procedure invdr
      module procedure invdc
      end interface

!     ----------------
      interface geneig
!     ----------------
      module procedure geneig2
      end interface

!     ----------------
      interface zeros
!     ----------------
      module procedure zeros2
      end interface

!     ----------------
      interface eye
!     ----------------
      module procedure eye1
      end interface

!     ----------------
      interface flipud
!     ----------------
      module procedure flipudr
      module procedure flipudc
      end interface

!     ----------------
      interface fliplr
!     ----------------
      module procedure fliplrr
      module procedure fliplrc
      end interface

!     -------------
      interface spy
!     -------------
      subroutine spyr(a,ni,nj,ndim)
        implicit none
        integer, parameter :: p8=selected_real_kind(p=11)
        integer :: ni,nj,ndim
        real(p8),dimension(ndim,1) :: a
      end subroutine

      subroutine spyc(a,ni,nj,ndim)
        implicit none
        integer, parameter :: p8=selected_real_kind(p=11)
        integer :: ni,nj,ndim
        complex(p8),dimension(ndim,1):: a
      end subroutine

      module procedure spydr
      module procedure spydc
      end interface

!     --------------
      interface mcat
!     --------------
      subroutine mcatr(a,ni,nj,ndim)
        implicit none
        integer, parameter :: p8=selected_real_kind(p=11)
        integer :: ni,nj,ndim
        real(p8),dimension(ndim,1):: a
      end subroutine

      subroutine mcatc(a,ni,nj,ndim)
        implicit none
        integer, parameter :: p8=selected_real_kind(p=11)
        integer :: ni,nj,ndim
        complex(p8),dimension(ndim,1):: a
      end subroutine

      module procedure mcatdr
      module procedure mcatdc
      module procedure mcat1dr
      module procedure mcat1dc
      end interface

!     ---------------
      interface mload
!     ---------------
!     module procedure mloadr
!     module procedure mloadc
      module procedure mload1r
      module procedure mload1c
      module procedure mloadr0
      module procedure mloadc0
      module procedure mloaddr
      module procedure mloaddc
      end interface

!     ---------------
      interface msave
!     ---------------
!     module procedure msaver
!     module procedure msavec
      module procedure msavedr
      module procedure msavedc
      module procedure msave1r
      module procedure msave1c
      module procedure msave3DimensionalDataComplex
      end interface

!     ----------------
      interface eigvec
!     ----------------
      module procedure eigvec10
      module procedure eigvec20
      module procedure eigvec1
      module procedure eigvec2
      end interface

!     --------------
      interface diag
!     --------------
      module procedure diag1r
      module procedure diag2r
      module procedure diag1c
      module procedure diag2c
      end interface

!     --------------
      interface lu
!     --------------
       module procedure lu_real0
       module procedure lu_complex0
      end interface

!     ---------------
      interface solve
!     ---------------
       module procedure solve_real1
       module procedure solve_real2
       module procedure solve_real1c
       module procedure solve_real2c
       module procedure solve_complex1
       module procedure solve_complex2
      end interface

!     ==================================================
      contains
!     --------


      subroutine lu_real0(a,w)
!     -------------------------------------------
      real(p8),dimension(:,:),intent(inout):: a
      integer,dimension(:,:),intent(inout):: w
      call lur(a,size(a,1),size(a,1),w(:,1),w(:,2))
      end subroutine


      subroutine lu_complex0(a,w)
!     -------------------------------------------
      complex(p8),dimension(:,:),intent(inout):: a
      integer,dimension(:,:),intent(inout):: w
      call luc(a,size(a,1),size(a,1),w(:,1),w(:,2))
      end subroutine



      subroutine solve_real1(a,w,b)
!     -------------------------------------------
      real(p8),dimension(:,:),intent(in):: a
      real(p8),dimension(:),intent(inout):: b
      integer,dimension(:,:),intent(in):: w
      integer:: ni,nj
      ni = size(a,1)
      nj = 1
      call solver(b,a,nj,ni,ni,w(:,1),w(:,2))
      end subroutine



      subroutine solve_real1c(a,w,b)
!     -------------------------------------------
      real(p8),dimension(:,:),intent(in):: a
      complex(p8),dimension(:),intent(inout):: b
      integer,dimension(:,:),intent(in):: w
      real(p8),dimension(size(b)):: wr,wi
      integer:: ni,nj
      ni = size(a,1)
      nj = 1
      wr = real(b)
      wi = aimag(b)
      call solver(wr,a,nj,ni,ni,w(:,1),w(:,2))
      call solver(wi,a,nj,ni,ni,w(:,1),w(:,2))
      b = cmplx(wr,wi,p8)
      end subroutine




      subroutine solve_real2c(a,w,b)
!     -------------------------------------------
      real(p8),dimension(:,:),intent(in):: a
      complex(p8),dimension(:,:),intent(inout):: b
      integer,dimension(:,:),intent(in):: w
      integer:: ni,nj,i
      ni = size(a,1)
      nj = size(b,2)
      do i=1,nj
        call solve_real1c(a,w,b(:,i))
      enddo
      end subroutine



      subroutine solve_real2(a,w,b)
!     -------------------------------------------
      real(p8),dimension(:,:),intent(in):: a
      real(p8),dimension(:,:),intent(inout):: b
      integer,dimension(:,:),intent(in):: w
      integer:: ni,nj
      ni = size(a,1)
      nj = size(b,2)
      call solver(b,a,nj,ni,ni,w(:,1),w(:,2))
      end subroutine



      subroutine solve_complex1(a,w,b)
!     -------------------------------------------
      complex(p8),dimension(:,:),intent(in):: a
      complex(p8),dimension(:),intent(inout):: b
      integer,dimension(:,:),intent(in):: w
      integer:: ni,nj
      ni = size(a,1)
      nj = 1
      call solvec(b,a,nj,ni,ni,w(:,1),w(:,2))
      end subroutine



      subroutine solve_complex2(a,w,b)
!     -------------------------------------------
      complex(p8),dimension(:,:),intent(in):: a
      complex(p8),dimension(:,:),intent(inout):: b
      integer,dimension(:,:),intent(in):: w
      integer:: ni,nj
      ni = size(a,1)
      nj = size(b,2)
      call solvec(b,a,nj,ni,ni,w(:,1),w(:,2))
      end subroutine
!
!
!
      function mulcr1(a,b)
!     -------------------------------------------
!     complex matrix - real vector multiply
!     -------------------------------------------
      complex(p8),dimension(:,:),intent(in):: a
      real(p8),dimension(:),intent(in):: b
      complex(p8),dimension(size(a,1)):: mulcr1
      integer :: nv
      nv = size(b)
      if(nv.ne.size(a,2)) then
	print *,'mulcr1: size mismatch.'
	print *,'[',size(a,1),',',size(a,2),']*[',nv,']'
	stop
      endif
      mulcr1 = matmul(a(:,1:nv),b)
      end function
!
!
!
      function mulrr1(a,b)
!     -------------------------------------------
!     real matrix - real vector multiply
!     -------------------------------------------
      real(p8),dimension(:,:),intent(in):: a
      real(p8),dimension(:),intent(in):: b
      real(p8),dimension(size(a,1)):: mulrr1
      integer :: nv
      nv = size(b)
      if(nv.ne.size(a,2)) then
	print *,'mulrr1: size mismatch.'
	print *,'[',size(a,1),',',size(a,2),']*[',nv,']'
	stop
      endif
      mulrr1 = matmul(a(:,1:nv),b)
      end function
!
!
!
      function mulrc1(a,b)
!     -------------------------------------------
!     real matrix - complex vector multiply
!     -------------------------------------------
      real(p8),dimension(:,:),intent(in):: a
      complex(p8),dimension(:),intent(in):: b
      complex(p8),dimension(size(a,1)):: mulrc1
      integer :: nv
      nv = size(b)
      if(nv.ne.size(a,2)) then
	print *,'mulrc1: size mismatch.'
	print *,'[',size(a,1),',',size(a,2),']*[',nv,']'
	stop
      endif
      mulrc1 = matmul(a(:,1:nv),b)
      end function
!
!
!
      function mulc1r(a,b)
!     -------------------------------------------
!     complex vector - real matrix multiply
!     -------------------------------------------
      complex(p8),dimension(:),intent(in):: a
      real(p8),dimension(:,:),intent(in):: b
      complex(p8),dimension(size(a,1)):: mulc1r
      if(size(a).ne.size(b,1)) then
	print *,'mulc1r: size mismatch.'
	print *,'[',size(a),']*[',size(b,1),',',size(b,2),']'
	stop
      endif
      mulc1r = (/ matmul(a,b) /)
      end function
!
!
!
      function mulcc1(a,b)
!     -------------------------------------------
!     real matrix - real vector multiply
!     -------------------------------------------
      complex(p8),dimension(:,:),intent(in):: a
      complex(p8),dimension(:),intent(in):: b
      complex(p8),dimension(size(a,1)):: mulcc1
      integer :: nv
      nv = size(b)
      if(nv.gt.size(a,2)) then
	print *,'mulcc1: size mismatch.'
	print *,'[',size(a,1),',',size(a,2),']*[',nv,']'
	stop
      endif
      mulcc1 = matmul(a(:,1:nv),b)
      end function
!
!
!
      function mulrr(a,b)
!     -------------------------------------------
!     real-real matrix multiply
!     -------------------------------------------
      real(p8),dimension(:,:),intent(in):: a,b
      real(p8),dimension(size(a,1),size(b,2)):: mulrr
      integer:: i,j,k
      real(p8),parameter:: one=1.0_p8
      if (size(a,2).ne.size(b,1)) then
	print *,'mul: invalid shape.'
	print *,'argument1 = ',size(a,1),'x',size(a,2)
	print *,'argument2 = ',size(b,1),'x',size(b,2)
	stop
      endif
!#BEGIN_SCALAR
      mulrr = matmul(a,b)
!#END_SCALAR
!#CRAY      i=size(a,1)
!#CRAY      j=size(b,2)
!#CRAY      k=size(a,2)
!#CRAY      call sgemm('n','n',i,j,k, one,a,i,b,k,0.0_p8, mulrr,i)
      
      end function
!
!
!
      function mulcc(a,b)
!     -------------------------------------------
!     complex-complex matrix multiply
!     -------------------------------------------
      complex(p8),dimension(:,:),intent(in):: a,b
      complex(p8),dimension(size(a,1),size(b,2)):: mulcc
      integer:: i,j,k
      complex(p8),parameter:: one =(1.0_p8,0.0_p8)
      complex(p8),parameter:: zero=(0.0_p8,0.0_p8)
      if (size(a,2).ne.size(b,1)) then
	print *,'mul: invalid shape.'
	print *,'argument1 = ',size(a,1),'x',size(a,2)
	print *,'argument2 = ',size(b,1),'x',size(b,2)
	stop
      endif
!#BEGIN_SCALAR
      !mulcc = matmul(a,b)
      mulcc = 0
      do i=1,size(a,1)
      do j=1,size(b,2)
      do k=1,size(a,2)
	mulcc(i,j) = mulcc(i,j)+a(i,k)*b(k,j)
      enddo
      enddo
      enddo
!#END_SCALAR
!#CRAY      i=size(a,1)
!#CRAY      j=size(b,2)
!#CRAY      k=size(a,2)
!#CRAY      call cgemm('n','n',i,j,k, one,a,i,b,k,zero, mulcc,i)
      end function
!
!
!
      function mulrc(a,b)
!     -------------------------------------------
!     real-complex matrix multiply
!     -------------------------------------------
      real(p8),dimension(:,:),intent(in):: a
      complex(p8),dimension(:,:),intent(in):: b
      complex(p8),dimension(size(a,1),size(b,2)):: mulrc
      real(p8),dimension(size(a,1),size(b,2)):: w1,w2
      if (size(a,2).ne.size(b,1)) then
	print *,'mul: invalid shape.'
	print *,'argument1 = ',size(a,1),'x',size(a,2)
	print *,'argument2 = ',size(b,1),'x',size(b,2)
	stop
      endif
      w1 = mulrr(a,real(b))
      w2 = mulrr(a,aimag(b))
      mulrc = cmplx(w1,w2,p8)
      end function
!
!
!
      function mulcr(a,b)
!     -------------------------------------------
!     complex-real matrix multiply
!     -------------------------------------------
      complex(p8),dimension(:,:),intent(in):: a
      real(p8),dimension(:,:),intent(in):: b
      complex(p8),dimension(size(a,1),size(b,2)):: mulcr
      real(p8),dimension(size(a,1),size(b,2)):: w1,w2
      if (size(a,2).ne.size(b,1)) then
	print *,'mul: invalid shape.'
	print *,'argument1 = ',size(a,1),'x',size(a,2)
	print *,'argument2 = ',size(b,1),'x',size(b,2)
	stop
      endif
      w1 = mulrr(real(a),b)
      w2 = mulrr(aimag(a),b)
      mulcr = cmplx(w1,w2,p8)
      end function

!
!
!
      function geneig2(a,b)
!     -------------------------------------------
!     solve  [a]=lamb[b]
!     for standard eigenvalue problem, set b=eye
!     -------------------------------------------
      complex(p8),dimension(:,:):: a
      complex(p8),dimension(:,:),optional:: b
      complex(p8),dimension(size(a,1),size(a,2)):: w1,w2
      complex(p8),dimension(size(a,1)):: geneig2
      real(p8),dimension(size(a,1)):: reig,aeig
      integer,dimension(size(a,1)):: indic
      integer :: ni,i

      if(size(a,1).ne.size(a,2)) then
	print *,'geneig2: matrix not square'
	print *,size(a,1),'x',size(a,2)
	return
      endif
      if(present(b)) then
      if(size(a).ne.size(b)) then
	print *,'geneig2: matrix size inconsistent'
	return
      endif
      endif

      ni=size(a,1)
      if(present(b)) then
        w1=a
        w2=b
      else
	w1=a
	w2=eye1(ni)
      endif
      call geneigc(w1,w2,ni,ni,reig,aeig,indic)
      geneig2 = cmplx(reig,aeig,p8)
      end function
!
!
!
      function eigvec10(a,eiv,adjm)
      complex(p8),dimension(:,:):: a
      complex(p8):: eiv
      complex(p8),dimension(:),optional:: adjm
      complex(p8),dimension(1):: eivwk
      complex(p8),dimension(size(a,1)):: eigvec10
      complex(p8),dimension(size(a,1),1):: adjwk
      eivwk(1)=eiv
      if(present(adjm)) then
        eigvec10=(/ eigvec1(a,eivwk,adjwk) /)
	adjm = adjwk(:,1)
      else
        eigvec10=(/ eigvec1(a,eivwk) /)
      endif
      end function
!
!
!
      function eigvec1(a,eiv,adjm)
!     -------------------------------------------
!     find eigenvectors for eiv
!     for standard eigenvalue problem, set b=eye
!     -------------------------------------------
      complex(p8),dimension(:,:):: a
      complex(p8),dimension(:):: eiv
      complex(p8),dimension(:,:),optional:: adjm
      complex(p8),dimension(size(a,1),size(eiv)):: eigvec1
      complex(p8),dimension(size(a,1),size(a,2)):: wk,b
      complex(p8),dimension(size(a,1)):: adj,vec
      complex(p8):: s=(1.0,1.0)
      integer :: ni,i,j

      ni=size(a,1)
      b=eye1(ni)
      do i=1,size(eiv)
        vec=s
        adj=s
        do j=1,3
          call eigrevc(a,b,ni,ni,wk,eiv(i),3,vec,adj)
        enddo
	eigvec1(:,i)=vec
	if(present(adjm)) then
	  adjm(:,i)=adj
	endif
      enddo
      end function
!
!
!
      function eigvec20(a,b,eiv,adjm)
      complex(p8),dimension(:,:):: a,b
      complex(p8):: eiv
      complex(p8),dimension(:),optional:: adjm
      complex(p8),dimension(1):: eivwk
      complex(p8),dimension(size(a,1)):: eigvec20
      complex(p8),dimension(size(a,1),1):: adjwk
      eivwk(1)=eiv
      if(present(adjm)) then
        eigvec20=(/ eigvec2(a,b,eivwk,adjwk) /)
	adjm = adjwk(:,1)
      else
        eigvec20=(/ eigvec2(a,b,eivwk) /)
      endif
      end function
!
!
!
      function eigvec2(a,b,eiv,adjm)
!     -------------------------------------------
!     find eigenvectors for eiv
!     for standard eigenvalue problem, set b=eye
!     -------------------------------------------
      complex(p8),dimension(:,:):: a,b
      complex(p8),dimension(:):: eiv
      complex(p8),dimension(:,:),optional:: adjm
      complex(p8),dimension(size(a,1),size(eiv,1)):: eigvec2
      complex(p8),dimension(size(a,1),size(a,2)):: wk
      complex(p8),dimension(size(a,1)):: adj,vec
      complex(p8):: s=(1.0,1.0)
      integer :: ni,i,j

      if(size(a).ne.size(b)) then
	print *,'eigvec2:size of matrices different'
	return
      endif

      ni=size(a,1)
      do i=1,size(eiv)
        vec=s
        adj=s
        do j=1,3
          call eigrevc(a,b,ni,ni,wk,eiv(i),3,vec,adj)
        enddo
	eigvec2(:,i)=vec
	if(present(adjm)) then
	  adjm(:,i)=adj
	endif
      enddo
      end function
!
!
!
      function eye1(ni)
!     -------------------------------
!     identity matrix of order ni
!     -------------------------------
      integer :: ni
      real(p8),dimension(ni,ni):: eye1
      integer :: i
      eye1=0
      do i=1,ni
	eye1(i,i)=1.0
      enddo
      end function
!
!
!
      function zeros2(ni,nj)
      integer :: ni,nj
      real(p8),dimension(ni,nj):: zeros2
      zeros2=0
      end function
!
!
!
      function diag1r(a)
!     -----------------------------------------
!     two dim array to one dim diagonal matrix.
!     -----------------------------------------
      real(p8),dimension(:,:) :: a
      real(p8),dimension(min(size(a,1),size(a,2))):: diag1r
      integer :: i

      do i=1,min(size(a,1),size(a,2))
	diag1r(i)=a(i,i)
      enddo
      end function
!
!
!
      function diag1c(a)
!     -----------------------------------------
!     two dim array to one dim diagonal matrix.
!     -----------------------------------------
      complex(p8),dimension(:,:) :: a
      complex(p8),dimension(min(size(a,1),size(a,2))):: diag1c
      integer :: i

      do i=1,min(size(a,1),size(a,2))
	diag1c(i)=a(i,i)
      enddo
      end function
!
!
!
      function diag2r(a)
!     -----------------------------------------
!     one dim array to two dim diagonal matrix.
!     -----------------------------------------
      real(p8),dimension(:) :: a
      real(p8),dimension(size(a),size(a)):: diag2r
      integer :: i
      diag2r=0
      do i=1,size(a)
	diag2r(i,i)=a(i)
      enddo
      end function
!
!
!
      function diag2c(a)
!     -----------------------------------------
!     one dim array to two dim diagonal matrix.
!     -----------------------------------------
      complex(p8),dimension(:) :: a
      complex(p8),dimension(size(a),size(a)):: diag2c
      integer :: i
      diag2c=0
      do i=1,size(a)
	diag2c(i,i)=a(i)
      enddo
      end function
!
!
!
      function flipudr(a)
!     -----------------------------------------
      real(p8),dimension(:,:),intent(in) :: a
      real(p8),dimension(size(a,1),size(a,2)):: flipudr
      integer :: i,n

      n = size(a,1)
      do i=1,n
	flipudr(i,:) = a(n-i+1,:)
      enddo
      end function
!
!
!
      function flipudc(a)
!     -----------------------------------------
      complex(p8),dimension(:,:),intent(in) :: a
      complex(p8),dimension(size(a,1),size(a,2)):: flipudc
      integer :: i,n

      n = size(a,1)
      do i=1,n
	flipudc(i,:) = a(n-i+1,:)
      enddo
      end function
!
!
!
      function fliplrr(a)
!     -----------------------------------------
      real(p8),dimension(:,:),intent(in) :: a
      real(p8),dimension(size(a,1),size(a,2)):: fliplrr
      integer :: j,n

      n = size(a,2)
      do j=1,n
	fliplrr(:,j) = a(:,n-j+1)
      enddo
      end function
!
!
!
      function fliplrc(a)
!     -----------------------------------------
      complex(p8),dimension(:,:),intent(in) :: a
      complex(p8),dimension(size(a,1),size(a,2)):: fliplrc
      integer :: j,n

      n = size(a,2)
      do j=1,n
	fliplrc(:,j) = a(:,n-j+1)
      enddo
      end function
!
!
!
      function invdr(a)
!     ----------------------------
!     invert a real matrix
!     ----------------------------
      real(p8),dimension(:,:):: a
      real(p8),dimension(size(a,1),size(a,2)):: invdr
      real(p8),dimension(size(a,1),size(a,2)):: c
      integer,dimension(size(a,1)):: ir,ic
      integer:: ni,nj

      ni=size(a,1)
      nj=size(a,2)
      if(ni.ne.nj) then
	print *,'invdr: matrix not square.'
	return
      endif
      c=a
      call lur(c,ni,ni,ir,ic)
      call eyer(invdr,ni,ni,size(invdr,1))
      call solver(invdr,c,ni,ni,size(invdr,1),ir,ic)
      end function
!
!
!
      function invdc(a)
!     invert a complex matrix
      complex(p8),dimension(:,:):: a
      complex(p8),dimension(size(a,1),size(a,2)):: invdc
      complex(p8),dimension(size(a,1),size(a,2)):: c
      integer,dimension(size(a,1)):: ir,ic
      integer:: ni,nj

      ni=size(a,1)
      nj=size(a,2)
      if(ni.ne.nj) then
	print *,'invdc: matrix not square.'
	return
      endif
      c=a
      call luc(c,ni,ni,ir,ic)
      call eyec(invdc,ni,ni,size(invdc,1))
      call solvec(invdc,c,ni,ni,size(invdc,1),ir,ic)
      end function
!
!
!
      subroutine spydr(a)
      real(p8),dimension(:,:),intent(in):: a
      call spy(a,size(a,1),size(a,2),size(a,1))
      end subroutine
!
!
!
      subroutine spydc(a)
      complex(p8),dimension(:,:),intent(in):: a
      call spyc(a,size(a,1),size(a,2),size(a,1))
      end subroutine
!
!
!
      subroutine mcatdr(a)
      real(p8),dimension(:,:),intent(in):: a
      call mcatr(a,size(a,1),size(a,2),size(a,1))
      end subroutine
!
!
!
      subroutine mcatdc(a)
      complex(p8),dimension(:,:),intent(in):: a
      call mcatc(a,size(a,1),size(a,2),size(a,1))
      end subroutine
!
!
!
      subroutine mcat1dr(a)
      real(p8),dimension(:),intent(in):: a
      integer:: i
      write(6,60) (i,a(i),i=1,size(a))
   60 format(i4,' ',1pe24.15e3)
      end subroutine
!
!
!
      subroutine mcat1dc(a)
      complex(p8),dimension(:),intent(in):: a
      integer:: i
      write(6,60) (i,a(i),i=1,size(a))
   60 format(i4,' (',1pe24.15e3,',',e24.15e3,')')
      end subroutine
!
!
!
      subroutine msavedr(a,fn)
      real(p8),dimension(:,:),intent(in) :: a
      character(len=*),intent(in):: fn
      call msaver(fn,a,size(a,1),size(a,2),size(a,1))
      end subroutine
!
!
!
      subroutine msavedc(a,fn)
      complex(p8),dimension(:,:),intent(in) :: a
      character(len=*),intent(in):: fn
      call msavec(fn,a,size(a,1),size(a,2),size(a,1))
      end subroutine
!
!
!
      subroutine msave1r(a,fn)
      real(p8),dimension(:),intent(in) :: a
      character(len=*),intent(in):: fn
      call msaver(fn,a,size(a),1,size(a))
      end subroutine
!
!
!
      subroutine msave1c(a,fn)
      complex(p8),dimension(:),intent(in) :: a
      character(len=*),intent(in):: fn
      call msavec(fn,a,size(a),1,size(a))
      end subroutine
!
!
!
      subroutine msaver(fn,a,ni,nj,ndim)
!     output a matrix for matlab.
!     use mloadc.m and msavec.m for complex matrices
!    ----------------------------------------------------
      implicit none
      integer :: ndim,ni,nj
      real(p8),dimension(ndim,nj):: a
      character(len=*),intent(in):: fn
      integer :: i,j
!
      call labopen(fn,nj)
      do 20 i=1,ni
      do 10 j=1,nj
	call labpr(a(i,j))
   10 continue
      call labcr
   20 continue
      call labclose
      write(6,60) ni,nj,trim(fn)
   60 format(i3,'x',i3,' written to ',a)
!
      return
      end subroutine
!
!
!
      subroutine msavec(fn,a,ni,nj,ndim)
!     output a matrix for matlab.
!     use mload.m to convert to complex matrix in matlab
!    ----------------------------------------------------
      integer :: ni,nj,ndim
      complex(p8),dimension(ndim,nj):: a
      character(len=*),intent(in):: fn
      integer i,j
!
      call labopen(fn,nj*2)
      do 20 i=1,ni
      do 10 j=1,nj
	call labpr(real(a(i,j)))
	call labpr(aimag(a(i,j)))
   10 continue
      call labcr
   20 continue
      call labclose
      write(6,60) ni,nj,trim(fn)
   60 format(i3,'x',i3,'c written to ',a)
!
      return
      end subroutine
!

      subroutine msave3DimensionalDataComplex(a, fn)
        !this is a wrapper
        complex(p8),dimension(:,:,:),intent(in) ::a
        character(len=*),intent(in):: fn
        call msaveDo3D(a,fn,size(a,1),size(a,2),size(a,3))        
      end subroutine msave3DimensionalDataComplex


      subroutine msaveDo3D(input3D, fileName, dim1, dim2, dim3)
!     output a 3D matrix for matlab to read.
!     Save one (2D) plane after another in one file 
!    ----------------------------------------------------
        complex(p8),dimension(dim1,dim2,dim3),intent(in):: input3D
        character(len=*),intent(in):: fileName
        integer, intent(in):: dim1, dim2, dim3
        integer i,j,k
        !
        !print *, dim1, dim2, dim3
        !stop
        call labopen(fileName,dim2*2)
        
        do k=1,dim3
           
           do i=1,dim1
              do j=1,dim2
                 call labpr(real(input3D(i,j,k)))
                 call labpr(aimag(input3D(i,j,k)))
              end do
              call labcr
           end do
        end do
        
        call labclose
        print *, dim1,'x',dim2,'x',dim3,'c written to ',trim(fileName)
        !write(6,60) dim1,dim2,dim3,trim(fileName)
        !60      format(i3,'x',i3,'x',i3,'c written to ',input3D)
        !
        return
      end subroutine msaveDo3D
      !
!
      subroutine mload1r(fn,a)
      character(len=*),intent(in):: fn
      real(p8),dimension(:) :: a
      real(p8),dimension(size(a),1) :: aa
      integer :: ni,nj
      ni=size(a,1)
      nj=1
      call mloadr(fn,aa,ni,nj,size(a))
      a = aa(:,1)
      end subroutine
!
!
!
      subroutine mload1c(fn,a)
      character(len=*),intent(in):: fn
      complex(p8),dimension(:) :: a
      complex(p8),dimension(size(a),1) :: aa
      integer :: ni,nj
      ni=size(a,1)
      nj=1
      call mloadc(fn,aa,ni,nj,size(a))
      a = aa(:,1)
      end subroutine
!
!
!
      subroutine mloaddr(fn,a)
      character(len=*),intent(in):: fn
      real(p8),dimension(:,:) :: a
      integer :: ni,nj
      ni=size(a,1)
      nj=size(a,2)
      call mloadr(fn,a,ni,nj,ni)
      end subroutine
!
!
!
      subroutine mloaddc(fn,a)
      character(len=*),intent(in):: fn
      complex(p8),dimension(:,:) :: a
      integer :: ni,nj
      ni=size(a,1)
      nj=size(a,2)
      call mloadc(fn,a,ni,nj,ni)
      end subroutine
!
!
!
      subroutine mloadr0(fn,a,ni,nj)
      character(len=*),intent(in):: fn
      real(p8),dimension(:,:) :: a
      integer,intent(out) :: ni,nj
      call mloadr(fn,a,ni,nj,size(a,1))
      end subroutine
!
!
!
      subroutine mloadc0(fn,a,ni,nj)
      character(len=*),intent(in):: fn
      complex(p8),dimension(:,:) :: a
      integer,intent(out) :: ni,nj
      call mloadc(fn,a,ni,nj,size(a,1))
      end subroutine
!
!
!
      subroutine mloadr(fn,a,ni,nj,ndim)
!     input a matrix from matlab ascii file
!     ni, nj is output
!    ----------------------------------------------------
      integer,intent(in) :: ndim
      integer,intent(out) :: ni,nj
!     real(p8),dimension(ndim,1):: a
      real(p8),dimension(:,:):: a
      character(len=*),intent(in):: fn
      integer :: i,j,is
      real(p8) :: val
!
      call labopen(fn,nj)
      i=1
      j=0
      nj=-1

      loop:do
	call labrd(val,is)
!       > value read normally
	if(is.eq.0) then
	  j=j+1
	  if(j > size(a,2).or. i>size(a,1)) then
	    print *,'mloadr:file size too large.'
	    print *,'mem:',size(a,1),'x',size(a,2)
	    print *,'fil:',i,'x',j
	    stop
	  endif
	  a(i,j)=val
	  cycle loop
!       > value read but EOR occured
	else if(is.eq.-1) then
!         >> first line ?
	  if(nj.eq.-1) then
	    j=j+1
	    if(j > size(a,2).or. i>size(a,1)) then
	      print *,'mloadr:file size too large.'
	      print *,'mem:',size(a,1),'x',size(a,2)
	      print *,'fil:',i,'x',j
	      stop
	    endif
	    a(i,j)=val
	    i=i+1
	    nj=j
	    j=0
	    cycle loop
	  else if(j.eq.nj-1) then
	    j=j+1
	    if(j > size(a,2).or. i>size(a,1)) then
	      print *,'mloadr:file size too large.'
	      print *,'mem:',size(a,1),'x',size(a,2)
	      print *,'fil:',i,'x',j
	      stop
	    endif
	    a(i,j)=val
	    i=i+1
	    j=0
	    cycle loop
	  else
	    print *,'mloadr: EOF at strange place.'
	    ni=i-1
	    exit loop
	  endif
!       > no value read and EOR occured
	else if(is.eq.-2) then
	  if(j.ne.0.and.nj.eq.-1) then
	    nj=j
	    j=0
	    i=i+1
	    cycle loop
	  else if(j.eq.nj) then
	    i=i+1
	    j=0
	    cycle loop
	  else if(j.eq.0) then
	    ni=i-1
	    exit loop
	  else
	    print *,'mload: EOF at strange place.'
	    ni=i-1
	    exit loop
	  endif
!       > unexpected character
        else
	  print *,'mload: strange character ?'
	  ni=0
	  nj=0
	  exit loop
	endif
      enddo loop

      call labclose
      write(6,60) ni,nj,trim(fn)
   60 format(i3,'x',i3,' read from ',a)
      return
      end subroutine
!
!
!
      subroutine mloadc(fn,a,ni,nj,ndim)
!     input a complex matrix from matlab ascii file
!     ni, nj is output
!    ----------------------------------------------------
      integer,intent(in) :: ndim
      integer,intent(out) :: ni,nj
!     complex(p8),dimension(ndim,1):: a
      complex(p8),dimension(:,:):: a
      character(len=*),intent(in):: fn
      integer :: i,j,is
      real(p8) :: val1,val2
!
      call labopen(fn,nj)
      i=1
      j=0
      nj=-1

      loop:do
	call labrd(val1,is)
	if(is.eq.0) then
	  call labrd(val2,is)
	endif
!       > value read normally
	if(is.eq.0) then
	  j=j+1
	  if(j>size(a,2).or. i>size(a,1)) then
	    print *,'mloadc:file size too large.'
	    print *,'mem:',size(a,1),'x',size(a,2)
	    print *,'fil:',i,'x',j
	    stop
	  endif
	  a(i,j)=cmplx(val1,val2,p8)
	  cycle loop
!       > value read but EOR occured
	else if(is.eq.-1) then
!         >> first line ?
	  if(nj.eq.-1) then
	    j=j+1
	    if(j>size(a,2).or. i>size(a,1)) then
	      print *,'mloadc:file size too large.'
	      print *,'mem:',size(a,1),'x',size(a,2)
	      print *,'fil:',i,'x',j
	      stop
	    endif
	    a(i,j)=cmplx(val1,val2,p8)
	    i=i+1
	    nj=j
	    j=0
	    cycle loop
	  else if(j.eq.nj-1) then
	    j=j+1
	    if(j>size(a,2).or. i>size(a,1)) then
	      print *,'mloadc:file size too large.'
	      print *,'mem:',size(a,1),'x',size(a,2)
	      print *,'fil:',i,'x',j
	      stop
	    endif
	    a(i,j)=cmplx(val1,val2,p8)
	    i=i+1
	    j=0
	    cycle loop
	  else
	    print *,'mload: EOF at strange place.'
	    ni=i-1
	    exit loop
	  endif
!       > no value read and EOR occured
	else if(is.eq.-2) then
	  if(j.ne.0.and.nj.eq.-1) then
	    nj=j
	    j=0
	    i=i+1
	    cycle loop
	  else if(j.eq.nj) then
	    i=i+1
	    j=0
	    cycle loop
	  else if(j.eq.0) then
	    ni=i-1
	    exit loop
	  else
	    print *,'mload: EOF at strange place.'
	    ni=i-1
	    exit loop
	  endif
!       > unexpected character
        else
	  print *,'mload: strange character ?'
	  ni=0
	  nj=0
	  exit loop
	endif
      enddo loop

      call labclose
      write(6,60) ni,nj,trim(fn)
   60 format(i3,'x',i3,'c read from ',a)
      return
      end subroutine
!
!
!
      subroutine labopen(fn,nc)
!     ----------------------------------------------
!     open a file
!     ----------------------------------------------
      character(len=*) :: fn
      integer :: nc
      integer :: idmy,reclen
      character(len=1):: dmy
      integer:: is
      reclen = 25*(nc+1)
      open(unit=11,file=fn,status='unknown',iostat=is,recl=reclen)
      if(is.ne.0) then
	print *,'labopen: can''t open',fn
      endif
!     > initialize
      call readchar(dmy,idmy,1)
      return
    end subroutine labopen
!
!
!
      subroutine labpr(a)
      real(p8),intent(in) :: a
      write(unit=11,fmt=60,advance="no") a
!     !!! if you modify the format(25 column), modify reclen in labopen.
   60 format(1pe25.16e3)
      return
      end subroutine
!
!
!
      subroutine labrd(a,is)
!     a : read value
!     is= 0 : number read, normal
!     is=-1 : number read, end of record/file
!     is=-2 : number not read, end of record/file
!     is=-3 : unexpected character encounter
!
      real(p8),intent(out) :: a
      integer :: is

      character(len=1) :: c
      character(len=72) :: str
!     real(p8) :: atof
      integer :: count

      count=0
      str=" "
      do
!       read(unit=11,fmt=60,advance="no",iostat=is) c
        call readchar(c,is)
   60   format(a1)
	if(is.ne.0) then
	  a=0.0
	  is=-2
          return
	else if(index(" ,	",c).ne.0) then
          cycle
	else
	  exit
	endif
      enddo

      count=1
      str(count:count)=c
      do
!       read(unit=11,fmt='(a1)',advance="no",iostat=is) c
        call readchar(c,is)
	if(is.ne.0) then
          a=atof(str)
	  is=-1
	  return
	else if(index("0123456789.eEdD+-",c).ne.0) then
	  count=count+1
          str(count:count)=c
	else if(index(" ,	",c).ne.0) then
          a=atof(str)
	  is=0
	  return
	else
	  print *,'labrd: unexpected character : ',c
	  is=-3
	endif
      enddo
!
      end subroutine
!
!
!
      subroutine readchar(c,is,init)
      character(len=1) :: c
      integer :: is
      integer,intent(in),optional:: init
      integer,save:: ii,ie,issav
      character(len=256),save:: buf=" "

      if(present(init)) then
	ii=0
	ie=0
	issav=0
	buf=" "
	return
      endif
!     > output buffered character if not EOR
      if(ii.lt.ie) then
        ii=ii+1
        c=buf(ii:ii)
        is=0
	return
      endif
!     > end of record ?
      if(issav.ne.0)then
	is=issav
	c=' '
	ii=0
	ie=0
	issav=0
	return
      endif
!     > read buffer
      ii=0
      buf=" "
      read(unit=11,fmt='(a256)',advance="no",iostat=issav) buf
      if (issav.eq.0) then
	ie=256
	ii=1
	c=buf(ii:ii)
	is=0
	return
      else
	ie=len(trim(buf))
	if(ie.eq.0) then
	  is=issav
	  c=' '
	  ii=0
	  issav=0
	  return
	endif
	ii=1
	c=buf(ii:ii)
	is=0
	return
      endif
      end subroutine
!
!
!
      subroutine labcr
      write(11,60) ' '
   60 format(a1)
      return
      end subroutine
!
!
!
      subroutine labclose
      close(11)
      return
      end subroutine
!
!
!
      real(p8) function atof(comm)
!    ----------------------------------------------------
!     ascii to floating number
      character(len=*),intent(in):: comm
      character(len=72) :: str
      real(p8) :: value
      str=comm(1:len(comm))//repeat(" ",72)
      read(unit=str,fmt='(d72.0)',err=911) value
!   8 format(d72.0)
      atof = value
      return
  911 continue
      write(0,*) 'atof: error in conversion'
      print *,'comm=',comm
      atof = 0.0
      return
      end function



      function paste_ud_rr(a,b)
      real(p8),dimension(:,:),intent(in):: a,b
      real(p8),dimension(size(a,1)+size(b,1),size(a,2)):: paste_ud_rr
      if(size(a,2).ne.size(b,2)) then
	print *,'paste_ud_rr: size mismatch.'
	print *,'sizeof(a)=[',size(a,1),',',size(a,2),']'
	print *,'sizeof(b)=[',size(b,1),',',size(b,2),']'
      endif
      paste_ud_rr(:size(a,1),:) = a
      paste_ud_rr(size(a,1)+1:,:) = b
      end function



      function paste_ud_cc(a,b)
      complex(p8),dimension(:,:),intent(in):: a,b
      complex(p8),dimension(size(a,1)+size(b,1),size(a,2)):: paste_ud_cc
      if(size(a,2).ne.size(b,2)) then
	print *,'paste_ud_cc: size mismatch.'
	print *,'sizeof(a)=[',size(a,1),',',size(a,2),']'
	print *,'sizeof(b)=[',size(b,1),',',size(b,2),']'
      endif
      paste_ud_cc(:size(a,1),:) = a
      paste_ud_cc(size(a,1)+1:,:) = b
      end function



      function paste_lr_rr(a,b)
      real(p8),dimension(:,:),intent(in):: a,b
      real(p8),dimension(size(a,1),size(a,2)+size(b,2)):: paste_lr_rr
      if(size(a,1).ne.size(b,1)) then
	print *,'paste_lr_rr: size mismatch.'
	print *,'sizeof(a)=[',size(a,1),',',size(a,2),']'
	print *,'sizeof(b)=[',size(b,1),',',size(b,2),']'
      endif
      paste_lr_rr(:,:size(a,2)) = a
      paste_lr_rr(:,size(a,2)+1:) = b
      end function



      function paste_lr_cc(a,b)
      complex(p8),dimension(:,:),intent(in):: a,b
      complex(p8),dimension(size(a,1),size(a,2)+size(b,2)):: paste_lr_cc
      if(size(a,1).ne.size(b,1)) then
	print *,'paste_lr_cc: size mismatch.'
	print *,'sizeof(a)=[',size(a,1),',',size(a,2),']'
	print *,'sizeof(b)=[',size(b,1),',',size(b,2),']'
      endif
      paste_lr_cc(:,:size(a,2)) = a
      paste_lr_cc(:,size(a,2)+1:) = b
      end function


      end module
!     ----------
!========================================================
!========================================================
!========================================================
!========================================================
!
!
!  UTILITIES: t.matsushima 8/3/95
!  ------------------------------
!
!
      subroutine invr(a,b,ni,ndim)
!     -------------------------------------
!     invert a matrix : b=inv(a)
!     a is destroyed
!     -------------------------------------
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: ni,ndim
      real(p8),dimension(ndim,ni):: a,b
      integer,parameter:: id=1024
      integer,dimension(id):: ir,ic
      if(ni.gt.id) then
	print *,'inv: matrix too large. increase id.'
	stop
      endif
      call lur(a,ni,ndim,ir,ic)
      call eyer(b,ni,ni,ndim)
      call solver(b,a,ni,ni,ndim,ir,ic)
      return
      end
!
!
!
      subroutine invc(a,b,ni,ndim)
!     -------------------------------------
!     invert complex matrix : b=inv(a)
!     a is destroyed
!     -------------------------------------
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer,parameter:: id=512
      integer :: ni,ndim
      complex(p8),dimension(ndim,ni):: a,b
      integer,dimension(id):: ir,ic
      if(ni.gt.id) then
	print *,'inv: matrix too large. increase id.'
	stop
      endif
      call luc(a,ni,ndim,ir,ic)
      call eyec(b,ni,ni,ndim)
      call solvec(b,a,ni,ni,ndim,ir,ic)
      return
      end
!
!
!
      subroutine eyer(a,ni,nj,ndim)
!     -------------------------------------
!     identity matrix
!     -------------------------------------
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer ni,nj,ndim
      real(p8),dimension(ndim,1) :: a
      integer i,j

      do 10 i=1,ni
      do 10 j=1,nj
        a(i,j)=0.0
   10 continue
      do 20 i=1,min(ni,nj)
        a(i,i)=1.0
   20 continue
      return
      end
!
!
!
      subroutine eyec(a,ni,nj,ndim)
!     -------------------------------------
!     (complex) identity matrix
!     -------------------------------------
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer ni,nj,ndim
      complex(p8),dimension(ndim,1) :: a
      integer i,j

      do 10 i=1,ni
      do 10 j=1,nj
        a(i,j)=0.0
   10 continue
      do 20 i=1,min(ni,nj)
        a(i,i)=1.0
   20 continue
      return
      end

!
!
!
      subroutine mcatr(a,ni,nj,ndim)
!     -------------------------------------
!     print real matrix
!     -------------------------------------
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: ni,nj,ndim
      real(p8),dimension(ndim,1):: a
      integer :: i,j,jb,je,jj

      do j=1,(nj+7)/8
        jb=8*j-7
        je=min(8*j,nj)
        write(6,*) ' '
        write(6,61) (jj,jj=jb,je)
        61 format(8i9)
        do i=1,ni
	  write(6,60) (a(i,jj),jj=jb,je)
          60 format(1p8e9.2)
        enddo
      enddo
      return
      end
!
!
!
      subroutine mcatc(a,ni,nj,ndim)
!     -------------------------------------
!     print complex matrix
!     -------------------------------------
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer ni,nj,ndim
      complex(p8),dimension(ndim,1) :: a
      integer :: i,j,jb,je,jj

      do 10 j=1,(nj+3)/4
        jb=4*j-3
        je=min(4*j,nj)
        write(6,*) ' '
        write(6,61) (jj,jj=jb,je)
   61   format(8i18)
        do 10 i=1,ni
	  write(6,60) (a(i,jj),jj=jb,je)
   60     format(1p8e9.2)
   10 continue
      return
      end
!
!
!
      subroutine spyr(a,ni,nj,ndim)
!     -------------------------------------
!     print structure of a real matrix
!     printed according to abs(a)
!     -------------------------------------
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: ni,nj,ndim
      real(p8),dimension(ndim,1) :: a
      integer,parameter:: nc=70,ncm=nc-1
      real(p8),parameter:: eps=1e-11
      character(len=80) :: b,c
      character(len=60) :: csr,tbl
      integer :: i,j,jb,je,jc,jj,ii,l,at
      real(p8) :: x

      b='----+----1----+----2----+----3----+----4'// &
        '----+----5----+----6----+----7----+----8'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      csr='log10:-2----+----1----+----0----+----1----+----2----'
      tbl='   (-) ________...ihgfedcba@ABCDEFGHI!!!!~~~~~~~~(+)'
      write(6,600) csr
      write(6,600) tbl
  600 format(a60)
      at=index(tbl,'@')
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do 10 j=1,(nj+ncm)/nc
        jb=nc*j-ncm
        je=min(nc*j,nj)
	write(6,*) ' '
        write(6,62) b(1:mod(je-1,nc)+1)
   62   format(3x,':',a)
	do 10 i=1,ni
	ii=1
	do 20 jj=jb,je
	  x=abs(a(i,jj))
	  jc=mod(jj-1,nc)+1
	  if(a(i,jj).gt.1.0_p8-eps .and. a(i,jj).lt.1.0_p8+eps) then
	    c(jc:jc)='1'
	  else if(x.eq.0.0) then
	    c(jc:jc)=' '
	  else 
	    l=log10(x)
	    if(l.lt.-20) then
	      c(jc:jc)='_'
	    else if(l.gt.20) then
	      c(jc:jc)='~'
	    else
	      l=l+at
	      c(jc:jc)=tbl(l:l)
	    endif
	  endif
   20   continue
	write(6,61) i,c(1:mod(je-1,nc)+1)
   61   format(i3,':',a)
   10 continue
      return
      end
!
!
!
      subroutine spyc(a,ni,nj,ndim)
!     -------------------------------------
!     print structure of a complex matrix
!     printed according to abs(a)
!     -------------------------------------
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: ni,nj,ndim
      complex(p8),dimension(ndim,1) :: a
      integer,parameter:: nc=70,ncm=nc-1
      real(p8),parameter:: eps=1e-11
      character(len=80) :: b,c
      character(len=60) :: csr,tbl
      integer :: i,j,jb,je,jc,jj,ii,l,at
      real(p8) :: x

      b='----+----1----+----2----+----3----+----4'// &
        '----+----5----+----6----+----7----+----8'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      csr='log10:-2----+----1----+----0----+----1----+----2----'
      tbl='   (-) ________...ihgfedcba@ABCDEFGHI!!!!~~~~~~~~(+)'
      write(6,600) csr
      write(6,600) tbl
  600 format(a60)
      at = index(tbl,'@')
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do 10 j=1,(nj+ncm)/nc
        jb=nc*j-ncm
        je=min(nc*j,nj)
	write(6,*) ' '
        write(6,62) b(1:mod(je-1,nc)+1)
   62   format(3x,':',a)
	do 10 i=1,ni
	ii=1
	do 20 jj=jb,je
	  x=abs(a(i,jj))
	  jc=mod(jj-1,nc)+1
	  if(x.gt.1.0_p8-eps .and. x.lt.1.0_p8+eps) then
	    if(real(a(i,jj)).gt.1.0_p8-eps) then
	      c(jc:jc)='1'
	    else
	      c(jc:jc)=']'
	    endif
	  else if(x.eq.0.0) then
	    c(jc:jc)=' '
	  else 
	    l=log10(x)
	    if(l.lt.-20) then
	      c(jc:jc)='_'
	    else if(l.gt.20) then
	      c(jc:jc)='~'
	    else
	      l=l+at
	      c(jc:jc)=tbl(l:l)
	    endif
	  endif
   20   continue
	write(6,61) i,c(1:mod(je-1,nc)+1)
   61   format(i3,':',a)
   10 continue
      return
      end
!========================================================
!========================================================
!  eig library:  to reformulate for f90, intrinsic function
!  cmplx(a,b) needed to be modified to cmplx(a,b,p8).
!  also all the numerical constants needed to be put _p8 at
!  their tails.
!========================================================
!========================================================
!
!
!
      subroutine atx(a,x,y,n,ndim)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: ndim,n
      real(kind=p8),dimension(ndim,1) :: a
      real(kind=p8),dimension(1):: x, y
      integer :: i,j
!
      do 1 i=1,n
      y(i)=0._p8
      do 2 j=1,n
    2 y(i) = y(i) + a(j,i) * x(j)
    1 continue
      return
      end
!
!
!
      subroutine atxc(a,x,y,n,ndim)
!    program to compute a transpose * x
!    the results are returned in y
      implicit none
      integer :: n,ndim
      integer, parameter :: p8=selected_real_kind(p=11)
      complex(kind=p8),dimension(ndim,1):: a
      complex(kind=p8),dimension(1) :: x, y
      integer :: i,j

      do 1 i=1,n
      y(i)=0.0_p8
      do 2 j=1,n
    2 y(i) = y(i) + a(j,i) * x(j)
    1 continue
      return
      end
!
!
!
      subroutine ax(a,x,y,n,ndim)
      implicit none
      integer :: n,ndim
      integer, parameter :: p8=selected_real_kind(p=11)
      real(kind=p8),dimension(ndim,1) :: a
      real(kind=p8),dimension(1) :: x, y
      integer :: i,j
      do 1 i=1,n
      y(i)=0.0_p8
      do 2 j=1,n
    2 y(i) = y(i) + a(i,j) * x(j)
    1 continue
      return
      end
!
!
!
      subroutine axc(a,x,y,n,ndim)
!    program to compute a*x where a is an n x n matrix and x is a vector
!    the results are returned in y
      implicit none
      integer :: n,ndim
      integer, parameter :: p8=selected_real_kind(p=11)
      complex(kind=p8),dimension(ndim,1) :: a
      complex(kind=p8),dimension(1) :: x,y
      integer :: i,j
      do 1 i=1,n
      y(i)=0.0_p8
      do 2 j=1,n
    2 y(i) = y(i) + a(i,j) * x(j)
    1 continue
      return
      end
!c
!c
!c
!c
!      function del(k,j,n)
!      implicit none
!      integer :: k,j,n,kr,jr,nr
!      integer, parameter :: p8=selected_real_kind(p=11)
!      real(kind=p8) :: fac,del,xnum,den
!      kr = k - 1
!      jr = j - 1
!      nr = n - 1
!      fac = 1._p8
!      if( kr .ne. nr ) go to 10
!      kr = 0
!      jr = nr - jr
!      fac = -1._p8
!10    del = fac*xnum(kr,jr,nr)/den(kr,jr,nr)
!      return
!      end
!c
!c
!c
!      function den(kr,jr,nr)
!      implicit none
!      integer, parameter :: p8=selected_real_kind(p=11)
!      integer :: kr,jr,nr
!      real(kind=p8) :: pi,den
!
!      pi = 4.0_p8*atan(1.0_p8)
!      if(kr .eq. 0) go to 20
!      den = real(nr,p8)*sin(pi*real(kr,p8)/real(nr,p8))
!      if(jr .eq. 0 .or. jr .eq. nr) den = den*2._p8
!      return
!20    den = .5_p8*real(nr,p8)
!      if(mod(jr,2) .eq. 1) den = - den
!      if(jr .eq. 0 .or. jr .eq. nr) den = 1._p8
!      return
!      end
!
!
!
!
      SUBROUTINE EIGCG1(NM,N,AR,AI,IFLAG,WR,WI,ZR,ZI,WORK1,WORK2,IERR)
!
!
! DIMENSION OF           AR(NM,N),AI(NM,N),WR(N),WI(N),ZR(NM,N),
! ARGUMENTS              ZI(NM,N),WORK1(NM,2*N),WORK2(2*N*N+4*N+2)
!
! LATEST REVISION        JANUARY 1985
!
! PURPOSE                EIGCG1 COMPUTES ALL EIGENVALUES OF A GENERAL
!                        COMPLEX MATRIX BY THE MODIFIED LR METHOD.  IT
!                        ALSO CAN COMPUTE ALL THE EIGENVECTORS OR IT
!                        CAN COMPUTE THE INFORMATION NECESSARY FOR
!                        SUBROUTINE EIGCG2 TO DETERMINE SELECTED
!                        EIGENVECTORS.
!
! USAGE                  CALL EIGCG1(NM,N,AR,AI,IFLAG,WR,WI,ZR,ZI,WORK1,
!                                 WORK2,IERR)
!
! ARGUMENTS
!
! ON INPUT               NM
!                          THE ROW (FIRST) DIMENSION IN THE CALLING
!                          PROGRAM OF THE TWO-DIMENSIONAL ARRAYS AR, AI,
!                          ZR, ZI, AND WORK1.
!
!                        N
!                          THE ORDER OF THE COMPLEX MATRIX REPRESENTED
!                          BY ARRAYS AR AND AI.
!
!                        AR, AI
!                          ARRAYS CONTAINING THE REAL AND IMAGINARY
!                          PARTS, RESPECTIVELY, OF THE MATRIX A.  THE
!                          ROW (FIRST) DIMENSION OF THESE ARRAYS MUST BE
!                          EQUAL TO NM.
!
!                        IFLAG
!                          AN INTEGER VARIABLE INDICATING THE NUMBER OF
!                          EIGENVECTORS TO BE COMPUTED.
!
!                          IFLAG = 0       NO EIGENVECTORS TO BE
!                                          DETERMINED BY EIGCG1
!                          0 < IFLAG < N   IFLAG SELECTED EIGENVECTORS
!                                          ARE TO BE FOUND USING
!                                          SUBROUTINE EIGCG2.
!                          IFLAG = N       ALL EIGENVECTORS TO BE
!                                          DETERMINED BY EIGCG1
!
!                        WORK1
!                          A DOUBLY-DIMENSIONED WORK ARRAY WITH
!                          ROW (FIRST) DIMENSION NM AND COLUMN DIMENSION
!                          AT LEAST 2*N.  EIGCG1 COPIES  AR  AND  AI
!                          TO WORK1 SO THAT THE ORIGINAL MATRIX IS
!                          NOT DESTROYED.  AFTER THE COPY,
!                          AR  AND  AI  ARE NOT USED AGAIN.  IF THE
!                          USER DOES NOT MIND IF HIS ORIGINAL
!                          MATRIX IS DESTROYED, HE CAN USE THE SAME
!                          CORE LOCATIONS FOR  AR, AI, AND WORK1
!                          WITH SOME ATTENTION GIVEN TO ENSURING
!                          THAT  AR AND AI  ARE STORED CONSECUTIVELY.
!                          IF EIGCG2 IS CALLED, WORK1 MUST BE INPUT TO
!                          EIGCG2 EXACTLY AS IT WAS OUTPUT FROM
!                          EIGCG1.
!
!                        WORK2
!                          A WORK ARRAY.  IF IFLAG = 0 (NO EIGENVECTOR
!                          DESIRED) OR IFLAG = N (ALL EIGENVECTORS
!                          DESIRED), THEN WORK2 MUST BE DIMENSIONED AT
!                          LEAST 2*N+2.  IF 0 .LT. IFLAG .LT. N (SOME
!                          EIGENVECTORS DESIRED), THEN WORK2 MUST BE
!                          DIMENSIONED AT LEAST 2*N*N + 4*N + 2.  IN
!                          THIS CASE, WORK2 MUST BE INPUT TO EIGCG2
!                          EXACTLY AS IT WAS OUTPUT FROM EIGCG1.
!
! ON OUTPUT              WR, WI
!                          ARRAYS OF DIMENSION N CONTAINING THE REAL AND
!                          IMAGINARY PARTS, RESPECTIVELY, OF THE
!                          EIGENVALUES.
!
!                        ZR, ZI
!                          IF IFLAG = 0 (NO EIGENVECTORS DESIRED), THEN
!                          ZR AND ZI ARE NOT USED, SO THEY MAY BE INPUT
!                          TO EIGCG1 AS SINGLE VARIABLES.  IF IFLAG > 0
!                          (SOME OR ALL EIGENVECTORS DESIRED), THEN ZR
!                          AND ZI ARE DOUBLY DIMENSIONED ARRAYS.  THE
!                          ROW (FIRST) DIMENSION IN THE CALLING PROGRAM
!                          MUST BE NM AND THE COLUMN DIMENSION AT LEAST
!                          N.  IF 0 < IFLAG < N (SOME EIGENVECTORS
!                          DESIRED), ZR AND ZI ARE USED AS WORK ARRAYS.
!                          IF IFLAG = N (ALL EIGENVECTORS DESIRED),
!                          THE COLUMNS OF ZR AND ZI CONTAIN THE REAL AND
!                          IMAGINARY PARTS, RESPECTIVELY, OF THE
!                          COMPUTED EIGENVECTORS.  THE EIGENVECTORS ARE
!                          UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE
!                          OF THE EIGENVECTORS HAS BEEN FOUND.
!
!                        IERR
!                          IS SET TO
!
!                          ZERO    FOR NORMAL RETURN.
!                          32+J    IF THE J-TH EIGENVALUE HAS NOT BEEN
!                                  DETERMINED AFTER 30 ITERATIONS.
!
! SPECIAL CONDITIONS     ALL EIGENVALUES ARE FOUND UNLESS AN ERROR
!                        MESSAGE INDICATES OTHERWISE.  IF THIS OCCURS,
!                        NO EIGENVECTORS WILL BE FOUND, AND SUBROUTINE
!                        EIGCG2 SHOULD NOT BE CALLED TO DETERMINE
!                        EIGENVECTORS.
!
! I/O                    IF AN ERROR OCCURS, AN ERROR MESSAGE WILL BE
!                        PRINTED.
!
! PRECISION              SINGLE
!
! REQUIRED LIBRARY       ULIBER AND Q8QST4, WHICH ARE LOADED BY
! FILES                  DEFAULT ON NCAR'S CRAY MACHINES.
!
! LANGUAGE               FORTRAN
!
! HISTORY                THIS PACKAGE WAS WRITTEN AT NCAR BY
!                        MEMBERS OF THE SCIENTIFIC COMPUTING
!                        DIVISION IN THE EARLY 1970'S.
!
! ALGORITHM              THIS SUBROUTINE IS PRINCIPALLY COMPOSED OF
!                        SUBROUTINES OBTAINED FROM EISPACK.
!
!                        THE MATRIX IS BALANCED, AND EIGENVALUES ARE
!                        ISOLATED WHENEVER POSSIBLE.  THE MATRIX IS
!                        THEN REDUCED TO COMPLEX UPPER HESSENBERG
!                        FORM USING STABILIZED ELEMENTARY SIMILARITY
!                        TRANSFORMATIONS.  THE EIGENVALUES ARE FOUND
!                        USING THE MODIFIED  LR  METHOD.  IF ALL
!                        EIGENVECTORS ARE DESIRED, THE  LR
!                        TRANSFORMATIONS ARE ACCUMULATED AND USED
!                        TO COMPUTE THE EIGENVECTORS OF THE
!                        BALANCED MATRIX.  THE EIGENVECTORS OF THE
!                        ORIGINAL MATRIX ARE FOUND FROM THOSE OF THE
!                        BALANCED MATRIX BY APPLYING THE BACK
!                        TRANSFORMATIONS USED IN BALANCING THE MATRIX.
!                        IF SOME EIGENVECTORS ARE DESIRED, INFORMATION
!                        ON THE TRANSFORMATIONS IS SAVED FOR USE BY
!                        EIGCG2  IN FINDING THE EIGENVECTORS.
!
! PORTABILITY            FORTRAN 66.  THERE ARE TWO MACHINE-DEPENDENT
!                        CONSTANTS DEFINED IN SUBROUTINES  CBAL, COMLR,
!                        AND  COMLR2: RADIX = 2. AND MACHEP = 2.**(-47).
!
! ACCURACY               DEPENDS ON THE TYPE OF MATRIX
!
! TIMING                 NO TIMING TESTS HAVE BEEN DONE.
!***********************************************************************
      IMPLICIT NONE
      INTEGER, PARAMETER :: P8=SELECTED_REAL_KIND(P=11)
      INTEGER :: NM,N
      REAL(kind=p8):: AR(NM,N),AI(NM,N),WR(N),WI(N),ZR(NM,1),ZI(NM,1)
      REAL(kind=p8)::  WORK1(NM,1),WORK2(1)
      INTEGER :: I,J,JPN,IFLAG,IERR

      DO 50 I = 1,N
      DO 50 J = 1,N
      WORK1(I,J) = AR(I,J)
      JPN = J+N
      WORK1(I,JPN) = AI(I,J)
   50 CONTINUE
      CALL CBAL(NM,N,WORK1,WORK1(1,N+1),WORK2(N+1),WORK2(N+2),WORK2)
      CALL COMHES(NM,N,WORK2(N+1),WORK2(N+2),WORK1,WORK1(1,N+1), &
      WORK2(N+3))
!
!     CHECK FOR NUMBER OF EIGENVECTORS DESIRED.
!
      IF (IFLAG .NE. 0) GO TO 100
!
!     NO EIGENVECTORS DESIRED.
!
      CALL COMLR(NM,N,WORK2(N+1),WORK2(N+2),WORK1,WORK1(1,N+1),&
      WR,WI,IERR)
      IF (IERR .NE. 0) GO TO 400
      RETURN
!
!     CHECK AGAIN FOR NUMBER OF EIGENVECTORS DESIRED.
!
  100 IF (IFLAG .NE. N) GO TO 200
!
!     ALL EIGENVECTORS DESIRED.
!
      CALL COMLR2(NM,N,WORK2(N+1),WORK2(N+2),WORK2(N+3),WORK1,&
      WORK1(1,N+1),WR,WI,ZR,ZI,IERR)
      IF (IERR .NE. 0) GO TO 400
      CALL CBABK2(NM,N,WORK2(N+1),WORK2(N+2),WORK2,N,ZR,ZI)
      RETURN
!
!     SOME EIGENVECTORS DESIRED.
!
  200 DO 300 I = 1,N
      DO 300 J = 1,N
      ZR(I,J) = WORK1(I,J)
      JPN = J+N
      ZI(I,J) = WORK1(I,JPN)
  300 CONTINUE
      CALL COMLR(NM,N,WORK2(N+1),WORK2(N+2),ZR,ZI,WR,WI,IERR)
      IF (IERR .EQ. 0) RETURN
  400 IERR = IABS(IERR) + 32
      WRITE(6,1212) IERR
1212  FORMAT(2X,I6,42H EIGCG1 FAILED TO COMPUTE ALL EIGENVALUES.  ,/)
!     CALL ULIBER(IERR,42H EIGCG1 FAILED TO COMPUTE ALL EIGENVALUES. ,
!     1            42)
      RETURN
!
!
!
      END
!
!
!
      SUBROUTINE CBAL(NM,N,AR,AI,LOW,IGH,SCALE)
      IMPLICIT NONE
      INTEGER, PARAMETER :: P8=SELECTED_REAL_KIND(P=11)
!
      INTEGER :: I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL(KIND=P8) AR(NM,N),AI(NM,N),SCALE(N)
      REAL(KIND=P8) C,F,G,R,S,B2,RADIX
!     REAL(KIND=P8) ABS
      LOGICAL NOCONV
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
!     CBALANCE, WHICH IS A COMPLEX VERSION OF BALANCE,
!     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
!
!     THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES
!     EIGENVALUES WHENEVER POSSIBLE.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE COMPLEX MATRIX TO BE BALANCED.
!
!     ON OUTPUT-
!
!        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE BALANCED MATRIX,
!
!        LOW AND IGH ARE TWO INTEGERS SUCH THAT AR(I,J) AND AI(I,J)
!          ARE EQUAL TO ZERO IF
!           (1) I IS GREATER THAN J AND
!           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,
!
!        SCALE CONTAINS INFORMATION DETERMINING THE
!           PERMUTATIONS AND SCALING FACTORS USED.
!
!     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
!     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
!     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
!     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
!        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
!                 = D(J,J)       J = LOW,...,IGH
!                 = P(J)         J = IGH+1,...,N.
!     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
!     THEN 1 TO LOW-1.
!
!     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
!
!     THE ALGOL PROCEDURE EXC CONTAINED IN CBALANCE APPEARS IN
!     CBAL  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
!     K,L HAVE BEEN REVERSED.)
!
!     ARITHMETIC IS REAL THROUGHOUT.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
!
!                **********
      RADIX = 2._p8
!
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
!     ********** IN-LINE PROCEDURE FOR ROW AND
!                COLUMN EXCHANGE **********
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
!
      DO 30 I = 1, L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   30 CONTINUE
!
      DO 40 I = K, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
!
   50 GO TO (80,130), IEXC
!     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                AND PUSH THEM DOWN **********
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
!     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
!
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (AR(J,I) .NE. 0.0 .OR. AI(J,I) .NE. 0.0) GO TO 120
  110    CONTINUE
!
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
!
      GO TO 140
!     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
!                AND PUSH THEM LEFT **********
  130 K = K + 1
!
  140 DO 170 J = K, L
!
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (AR(I,J) .NE.0.0_p8.OR. AI(I,J) .NE.0.0_p8) GO TO 170
  150    CONTINUE
!
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
!     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********
      DO 180 I = K, L
  180 SCALE(I) = 1.0_p8
!     ********** ITERATIVE LOOP FOR NORM REDUCTION **********
  190 NOCONV = .FALSE.
!
      DO 270 I = K, L
         C = 0.0_p8
         R = 0.0_p8
!
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(AR(J,I)) + ABS(AI(J,I))
            R = R + ABS(AR(I,J)) + ABS(AI(I,J))
  200    CONTINUE
!
         G = R / RADIX
         F = 1.0_p8
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
!     ********** NOW BALANCE **********
  240    IF ((C + R) / F .GE. 0.95_p8 * S) GO TO 270
         G = 1.0_p8 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
!
         DO 250 J = K, N
            AR(I,J) = AR(I,J) * G
            AI(I,J) = AI(I,J) * G
  250    CONTINUE
!
         DO 260 J = 1, L
            AR(J,I) = AR(J,I) * F
            AI(J,I) = AI(J,I) * F
  260    CONTINUE
!
  270 CONTINUE
!
      IF (NOCONV) GO TO 190
!
  280 LOW = K
      IGH = L
      RETURN
!     ********** LAST CARD OF CBAL **********
!
!
!
      END
!
!
!
      SUBROUTINE COMHES(NM,N,LOW,IGH,AR,AI,INT)
      IMPLICIT NONE
      INTEGER, PARAMETER :: P8=SELECTED_REAL_KIND(P=11)

!
      INTEGER :: I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL(KIND=p8) :: AR(NM,N),AI(NM,N)
      REAL(KIND=p8) :: XR,XI,YR,YI
!     REAL(KIND=p8) :: ABS
      INTEGER :: INT(IGH)
      COMPLEX(KIND=p8) :: Z3
!     COMPLEX(KIND=p8) :: CMPLX
!     REAL(KIND=p8) ::  REAL,AIMAG
!
!
!
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMHES,
!     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
!     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
!     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
!     STABILIZED ELEMENTARY SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
!          SET LOW=1, IGH=N,
!
!        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
!
!     ON OUTPUT-
!
!        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE HESSENBERG MATRIX.  THE
!          MULTIPLIERS WHICH WERE USED IN THE REDUCTION
!          ARE STORED IN THE REMAINING TRIANGLES UNDER THE
!          HESSENBERG MATRIX,
!
!        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
!          INTERCHANGED IN THE REDUCTION.
!          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
!
!     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL
!     PROCEDURE CDIV BY COMPLEX DIVISION USING SUBROUTINE CMPLX.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
!
      DO 180 M = KP1, LA
         MM1 = M - 1
         XR = 0.0_p8
         XI = 0.0_p8
         I = M
!
         DO 100 J = M, IGH
            IF (ABS(AR(J,MM1)) + ABS(AI(J,MM1)) &
               .LE. ABS(XR) + ABS(XI)) GO TO 100
            XR = AR(J,MM1)
            XI = AI(J,MM1)
            I = J
  100    CONTINUE
!
         INT(M) = I
         IF (I .EQ. M) GO TO 130
!     ********** INTERCHANGE ROWS AND COLUMNS OF AR AND AI **********
         DO 110 J = MM1, N
            YR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = YR
            YI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = YI
  110    CONTINUE
!
         DO 120 J = 1, IGH
            YR = AR(J,I)
            AR(J,I) = AR(J,M)
            AR(J,M) = YR
            YI = AI(J,I)
            AI(J,I) = AI(J,M)
            AI(J,M) = YI
  120    CONTINUE
!     ********** END INTERCHANGE **********
  130    IF (XR .EQ. 0.0 .AND. XI .EQ. 0.0) GO TO 180
         MP1 = M + 1
!
         DO 160 I = MP1, IGH
            YR = AR(I,MM1)
            YI = AI(I,MM1)
            IF (YR .EQ. 0.0 .AND. YI .EQ. 0.0) GO TO 160
            Z3 = CMPLX(YR,YI,p8) / CMPLX(XR,XI,p8)
            YR = REAL(Z3)
            YI = AIMAG(Z3)
            AR(I,MM1) = YR
            AI(I,MM1) = YI
!
            DO 140 J = M, N
               AR(I,J) = AR(I,J) - YR * AR(M,J) + YI * AI(M,J)
               AI(I,J) = AI(I,J) - YR * AI(M,J) - YI * AR(M,J)
  140       CONTINUE
!
            DO 150 J = 1, IGH
               AR(J,M) = AR(J,M) + YR * AR(J,I) - YI * AI(J,I)
               AI(J,M) = AI(J,M) + YR * AI(J,I) + YI * AR(J,I)
  150       CONTINUE
!
  160    CONTINUE
!
  180 CONTINUE
!
  200 RETURN
!     ********** LAST CARD OF COMHES **********
!
!
!
      END
!
!
!
      SUBROUTINE COMLR(NM,N,LOW,IGH,HR,HI,WR,WI,IERR)
      IMPLICIT NONE
      INTEGER, PARAMETER :: P8=SELECTED_REAL_KIND(P=11)
!
      INTEGER ::I,J,L,M,N,EN,LL,MM,NM,IGH,IM1,ITS,LOW,MP1,ENM1,IERR
      REAL(KIND=p8):: HR(NM,N),HI(NM,N),WR(N),WI(N)
      REAL(KIND=p8):: SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,MACHEP
!     REAL(KIND=p8):: ABS
      COMPLEX(KIND=p8):: Z3
!     COMPLEX(KIND=p8):: CSQRT,CMPLX
!     REAL(KIND=p8):: REAL,AIMAG
!
!
!
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMLR,
!     NUM. MATH. 12, 369-376(1968) BY MARTIN AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES OF A COMPLEX
!     UPPER HESSENBERG MATRIX BY THE MODIFIED LR METHOD.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
!          SET LOW=1, IGH=N,
!
!        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
!          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN THE
!          MULTIPLIERS WHICH WERE USED IN THE REDUCTION BY  COMHES,
!          IF PERFORMED.
!
!     ON OUTPUT-
!
!        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
!          DESTROYED.  THEREFORE, THEY MUST BE SAVED BEFORE
!          CALLING  COMLR  IF SUBSEQUENT CALCULATION OF
!          EIGENVECTORS IS TO BE PERFORMED,
!
!        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
!          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
!          FOR INDICES IERR+1,...,N,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL
!     PROCEDURE CDIV BY COMPLEX DIVISION AND USE OF THE SUBROUTINES
!     CSQRT AND CMPLX IN COMPUTING COMPLEX SQUARE ROOTS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!
!                **********
      MACHEP = 2._p8**(-47)
!
      IERR = 0
!     ********** STORE ROOTS ISOLATED BY CBAL **********
      DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
!
      EN = IGH
      TR = 0.0_p8
      TI = 0.0_p8
!     ********** SEARCH FOR NEXT EIGENVALUE **********
  220 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      ENM1 = EN - 1
!     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!                FOR L=EN STEP -1 UNTIL LOW  -- **********
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         IF (ABS(HR(L,L-1)) + ABS(HI(L,L-1)) .LE. &
            MACHEP * (ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1)) &
                   + ABS(HR(L,L)) + ABS(HI(L,L)))) GO TO 300
  260 CONTINUE
!     ********** FORM SHIFT **********
  300 IF (L .EQ. EN) GO TO 660
      IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)
      XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0_p8 .AND. XI .EQ. 0.0_p8) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0_p8
      YI = (HI(ENM1,ENM1) - SI) / 2.0_p8
      Z3 = SQRT(CMPLX(YR**2-YI**2+XR,2.0_p8*YR*YI+XI,p8))
      ZZR = REAL(Z3)
      ZZI = AIMAG(Z3)
      IF (YR * ZZR + YI * ZZI .GE. 0.0_p8) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 Z3 = CMPLX(XR,XI,p8) / CMPLX(YR+ZZR,YI+ZZI,p8)
      SR = SR - REAL(Z3)
      SI = SI - AIMAG(Z3)
      GO TO 340
!     ********** FORM EXCEPTIONAL SHIFT **********
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))
!
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
!
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
!     ********** LOOK FOR TWO CONSECUTIVE SMALL
!                SUB-DIAGONAL ELEMENTS **********
      XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))
      YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))
      ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))
!     ********** FOR M=EN-1 STEP -1 UNTIL L DO -- **********
      DO 380 MM = L, ENM1
         M = ENM1 + L - MM
         IF (M .EQ. L) GO TO 420
         YI = YR
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
         XI = ZZR
         ZZR = XR
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
         IF (YR .LE. MACHEP * ZZR / YI * (ZZR + XR + XI)) GO TO 420
  380 CONTINUE
!     ********** TRIANGULAR DECOMPOSITION H=L*R **********
  420 MP1 = M + 1
!
      DO 520 I = MP1, EN
         IM1 = I - 1
         XR = HR(IM1,IM1)
         XI = HI(IM1,IM1)
         YR = HR(I,IM1)
         YI = HI(I,IM1)
         IF (ABS(XR) + ABS(XI) .GE. ABS(YR) + ABS(YI)) GO TO 460
!     ********** INTERCHANGE ROWS OF HR AND HI **********
         DO 440 J = IM1, EN
            ZZR = HR(IM1,J)
            HR(IM1,J) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(IM1,J)
            HI(IM1,J) = HI(I,J)
            HI(I,J) = ZZI
  440    CONTINUE
!
         Z3 = CMPLX(XR,XI,p8) / CMPLX(YR,YI,p8)
         WR(I) = 1.0_p8
         GO TO 480
  460    Z3 = CMPLX(YR,YI,p8) / CMPLX(XR,XI,p8)
         WR(I) = -1.0_p8
  480    ZZR = REAL(Z3)
         ZZI = AIMAG(Z3)
         HR(I,IM1) = ZZR
         HI(I,IM1) = ZZI
!
         DO 500 J = I, EN
            HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)
            HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)
  500    CONTINUE
!
  520 CONTINUE
!     ********** COMPOSITION R*L=H **********
      DO 640 J = MP1, EN
         XR = HR(J,J-1)
         XI = HI(J,J-1)
         HR(J,J-1) = 0.0_p8
         HI(J,J-1) = 0.0_p8
!     ********** INTERCHANGE COLUMNS OF HR AND HI,
!                IF NECESSARY **********
         IF (WR(J) .LE. 0.0_p8) GO TO 580
!
         DO 540 I = L, J
            ZZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZZI
  540    CONTINUE
!
  580    DO 600 I = L, J
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
!
  640 CONTINUE
!
      GO TO 240
!     ********** A ROOT FOUND **********
  660 WR(EN) = HR(EN,EN) + TR
      WI(EN) = HI(EN,EN) + TI
      EN = ENM1
      GO TO 220
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = EN
 1001 RETURN
!     ********** LAST CARD OF COMLR **********
!
!
!
      END
!
!
!
      SUBROUTINE COMLR2(NM,N,LOW,IGH,INT,HR,HI,WR,WI,ZR,ZI,IERR)
      IMPLICIT NONE
      INTEGER, PARAMETER :: P8=SELECTED_REAL_KIND(P=11)
!
      INTEGER ::I,J,K,L,M,N,EN,II,JJ,LL,MM,NM,NN,IGH,IM1,IP1, &
              ITS,LOW,MP1,ENM1,IEND,IERR
      REAL(KIND=p8):: HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N)
      REAL(KIND=p8):: SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,MACHEP
!     REAL(KIND=p8):: ABS
      INTEGER:: INT(IGH)
!     INTEGER:: MIN0
      COMPLEX(KIND=p8):: Z3
!     COMPLEX(KIND=p8):: CSQRT,CMPLX
!     REAL(KIND=p8):: REAL,AIMAG
!
!
!
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMLR2,
!     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE MODIFIED LR
!     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
!     CAN ALSO BE FOUND IF  COMHES  HAS BEEN USED TO REDUCE
!     THIS GENERAL MATRIX TO HESSENBERG FORM.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
!          SET LOW=1, IGH=N,
!
!        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS INTERCHANGED
!          IN THE REDUCTION BY  COMHES, IF PERFORMED.  ONLY ELEMENTS
!          LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS OF THE HESSEN-
!          BERG MATRIX ARE DESIRED, SET INT(J)=J FOR THESE ELEMENTS,
!
!        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
!          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN THE
!          MULTIPLIERS WHICH WERE USED IN THE REDUCTION BY  COMHES,
!          IF PERFORMED.  IF THE EIGENVECTORS OF THE HESSENBERG
!          MATRIX ARE DESIRED, THESE ELEMENTS MUST BE SET TO ZERO.
!
!     ON OUTPUT-
!
!        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
!          DESTROYED, BUT THE LOCATION HR(1,1) CONTAINS THE NORM
!          OF THE TRIANGULARIZED MATRIX,
!
!        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
!          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
!          FOR INDICES IERR+1,...,N,
!
!        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
!          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
!          THE EIGENVECTORS HAS BEEN FOUND,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL
!     PROCEDURE CDIV BY COMPLEX DIVISION AND USE OF THE SUBROUTINES
!     CSQRT AND CMPLX IN COMPUTING COMPLEX SQUARE ROOTS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!
!                **********
      MACHEP = 2._p8**(-47)
!
      IERR = 0
!     ********** INITIALIZE EIGENVECTOR MATRIX **********
      DO 100 I = 1, N
!
         DO 100 J = 1, N
            ZR(I,J) = 0.0_p8
            ZI(I,J) = 0.0_p8
            IF (I .EQ. J) ZR(I,J) = 1.0_p8
  100 CONTINUE
!     ********** FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
!                FROM THE INFORMATION LEFT BY COMHES **********
      IEND = IGH - LOW - 1
      IF (IEND .LE. 0) GO TO 180
!     ********** FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********
      DO 160 II = 1, IEND
         I = IGH - II
         IP1 = I + 1
!
         DO 120 K = IP1, IGH
            ZR(K,I) = HR(K,I-1)
            ZI(K,I) = HI(K,I-1)
  120    CONTINUE
!
         J = INT(I)
         IF (I .EQ. J) GO TO 160
!
         DO 140 K = I, IGH
            ZR(I,K) = ZR(J,K)
            ZI(I,K) = ZI(J,K)
            ZR(J,K) = 0.0_p8
            ZI(J,K) = 0.0_p8
  140    CONTINUE
!
         ZR(J,I) = 1.0_p8
  160 CONTINUE
!     ********** STORE ROOTS ISOLATED BY CBAL **********
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
!
      EN = IGH
      TR = 0.0_p8
      TI = 0.0_p8
!     ********** SEARCH FOR NEXT EIGENVALUE **********
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
!     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!                FOR L=EN STEP -1 UNTIL LOW DO -- **********
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         IF (ABS(HR(L,L-1)) + ABS(HI(L,L-1)) .LE. &
            MACHEP * (ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1)) &
                   + ABS(HR(L,L)) + ABS(HI(L,L)))) GO TO 300
  260 CONTINUE
!     ********** FORM SHIFT **********
  300 IF (L .EQ. EN) GO TO 660
      IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)
      XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0_p8 .AND. XI .EQ. 0.0_p8) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0_p8
      YI = (HI(ENM1,ENM1) - SI) / 2.0_p8
      Z3 = SQRT(CMPLX(YR**2-YI**2+XR,2.0_p8*YR*YI+XI,p8))
      ZZR = REAL(Z3)
      ZZI = AIMAG(Z3)
      IF (YR * ZZR + YI * ZZI .GE. 0.0_p8) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 Z3 = CMPLX(XR,XI,p8) / CMPLX(YR+ZZR,YI+ZZI,p8)
      SR = SR - REAL(Z3)
      SI = SI - AIMAG(Z3)
      GO TO 340
!     ********** FORM EXCEPTIONAL SHIFT **********
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))
!
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
!
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
!     ********** LOOK FOR TWO CONSECUTIVE SMALL
!                SUB-DIAGONAL ELEMENTS **********
      XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))
      YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))
      ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))
!     ********** FOR M=EN-1 STEP -1 UNTIL L DO -- **********
      DO 380 MM = L, ENM1
         M = ENM1 + L - MM
         IF (M .EQ. L) GO TO 420
         YI = YR
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
         XI = ZZR
         ZZR = XR
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
         IF (YR .LE. MACHEP * ZZR / YI * (ZZR + XR + XI)) GO TO 420
  380 CONTINUE
!     ********** TRIANGULAR DECOMPOSITION H=L*R **********
  420 MP1 = M + 1
!
      DO 520 I = MP1, EN
         IM1 = I - 1
         XR = HR(IM1,IM1)
         XI = HI(IM1,IM1)
         YR = HR(I,IM1)
         YI = HI(I,IM1)
         IF (ABS(XR) + ABS(XI) .GE. ABS(YR) + ABS(YI)) GO TO 460
!     ********** INTERCHANGE ROWS OF HR AND HI **********
         DO 440 J = IM1, N
            ZZR = HR(IM1,J)
            HR(IM1,J) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(IM1,J)
            HI(IM1,J) = HI(I,J)
            HI(I,J) = ZZI
  440    CONTINUE
!
         Z3 = CMPLX(XR,XI,p8) / CMPLX(YR,YI,p8)
         WR(I) = 1.0_p8
         GO TO 480
  460    Z3 = CMPLX(YR,YI,p8) / CMPLX(XR,XI,p8)
         WR(I) = -1.0_p8
  480    ZZR = REAL(Z3)
         ZZI = AIMAG(Z3)
         HR(I,IM1) = ZZR
         HI(I,IM1) = ZZI
!
         DO 500 J = I, N
            HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)
            HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)
  500    CONTINUE
!
  520 CONTINUE
!     ********** COMPOSITION R*L=H **********
      DO 640 J = MP1, EN
         XR = HR(J,J-1)
         XI = HI(J,J-1)
         HR(J,J-1) = 0.0_p8
         HI(J,J-1) = 0.0_p8
!     ********** INTERCHANGE COLUMNS OF HR, HI, ZR, AND ZI,
!                IF NECESSARY **********
         IF (WR(J) .LE. 0.0_p8) GO TO 580
!
         DO 540 I = 1, J
            ZZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZZI
  540    CONTINUE
!
         DO 560 I = LOW, IGH
            ZZR = ZR(I,J-1)
            ZR(I,J-1) = ZR(I,J)
            ZR(I,J) = ZZR
            ZZI = ZI(I,J-1)
            ZI(I,J-1) = ZI(I,J)
            ZI(I,J) = ZZI
  560    CONTINUE
!
  580    DO 600 I = 1, J
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
!     ********** ACCUMULATE TRANSFORMATIONS **********
         DO 620 I = LOW, IGH
            ZR(I,J-1) = ZR(I,J-1) + XR * ZR(I,J) - XI * ZI(I,J)
            ZI(I,J-1) = ZI(I,J-1) + XR * ZI(I,J) + XI * ZR(I,J)
  620    CONTINUE
!
  640 CONTINUE
!
      GO TO 240
!     ********** A ROOT FOUND **********
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
!     ********** ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
!                VECTORS OF UPPER TRIANGULAR FORM **********
  680 NORM = 0.0_p8
!
      DO 720 I = 1, N
!
         DO 720 J = I, N
            NORM = NORM + ABS(HR(I,J)) + ABS(HI(I,J))
  720 CONTINUE
!
      HR(1,1) = NORM
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0_p8) GO TO 1001
!     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         ENM1 = EN - 1
!     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = HR(I,EN)
            ZZI = HI(I,EN)
            IF (I .EQ. ENM1) GO TO 760
            IP1 = I + 1
!
            DO 740 J = IP1, ENM1
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
!
  760       YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .EQ. 0.0_p8 .AND. YI .EQ. 0.0_p8) YR = MACHEP * NORM
            Z3 = CMPLX(ZZR,ZZI,p8) / CMPLX(YR,YI,p8)
            HR(I,EN) = REAL(Z3)
            HI(I,EN) = AIMAG(Z3)
  780    CONTINUE
!
  800 CONTINUE
!     ********** END BACKSUBSTITUTION **********
      ENM1 = N - 1
!     ********** VECTORS OF ISOLATED ROOTS **********
      DO 840 I = 1, ENM1
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
         IP1 = I + 1
!
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
!
  840 CONTINUE
!     ********** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
!                VECTORS OF ORIGINAL FULL MATRIX.
!                FOR J=N STEP -1 UNTIL LOW+1 DO -- **********
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN0(J-1,IGH)
!
         DO 880 I = LOW, IGH
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
!
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
!
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
!
      GO TO 1001
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = EN
 1001 RETURN
!     ********** LAST CARD OF COMLR2 **********
!
!
!
      END
!
!
!
      SUBROUTINE CBABK2(NM,N,LOW,IGH,SCALE,M,ZR,ZI)
      IMPLICIT NONE
      INTEGER, PARAMETER :: P8=SELECTED_REAL_KIND(P=11)
!
      INTEGER :: I,J,K,M,N,II,NM,IGH,LOW
      REAL(KIND=P8):: SCALE(N),ZR(NM,M),ZI(NM,M)
      REAL(KIND=P8):: S
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
!     CBABK2, WHICH IS A COMPLEX VERSION OF BALBAK,
!     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
!
!     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL
!     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!     BALANCED MATRIX DETERMINED BY  CBAL.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        LOW AND IGH ARE INTEGERS DETERMINED BY  CBAL,
!
!        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
!          AND SCALING FACTORS USED BY  CBAL,
!
!        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
!
!        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE EIGENVECTORS TO BE
!          BACK TRANSFORMED IN THEIR FIRST M COLUMNS.
!
!     ON OUTPUT-
!
!        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
!          IN THEIR FIRST M COLUMNS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
!
      DO 110 I = LOW, IGH
         S = SCALE(I)
!     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
!                IF THE FOREGOING STATEMENT IS REPLACED BY
!                S=1.0/SCALE(I). **********
         DO 100 J = 1, M
            ZR(I,J) = ZR(I,J) * S
            ZI(I,J) = ZI(I,J) * S
  100    CONTINUE
!
  110 CONTINUE
!     ********** FOR I=LOW-1 STEP -1 UNTIL 1,
!                IGH+1 STEP 1 UNTIL N DO -- **********
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
!
         DO 130 J = 1, M
            S = ZR(I,J)
            ZR(I,J) = ZR(K,J)
            ZR(K,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(K,J)
            ZI(K,J) = S
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
!
! REVISION HISTORY---
!
!
! JANUARY 1978     DELETED REFERENCES TO THE  *COSY  CARDS AND
!                  ADDED REVISION HISTORY AND REMOVED COMMENT CARDS
!                  BETWEEN SUBROUTINES
!
! APRIL   1985     FIXED THE DIMENSION OF ARGUMENTS IN DOCUMENTATION.
!-----------------------------------------------------------------------
!     ********** LAST CARD OF CBABK2 **********
      END
!
!
!
      subroutine eigrev(a,b,n,ndim,c,lamda,lll,eig,adj)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
!        this code assumes ndim is 33 but computes
!       fac1 sensibly
      integer :: ndim,n,lll
      real(kind=p8),dimension(ndim,1):: a, b, c
!     dimension eig(1),adj(1),t1(ndim),t2(ndim)
      real(kind=p8),dimension(1):: eig(1),adj(1)
      real(kind=p8),dimension(1024):: t1(1024),t2(1024)
      real(kind=p8):: lamda
      real(kind=p8),dimension(1024):: ir,ic
      real(kind=p8):: g,h,sum1,sum2,fac1,fac2
      integer :: l,i,j
!
!      if n (which is equal to the actual length of the eigenvector)
!      is greater than 100 than you must modify the above card
!      by change #100#=#some integer >= n#
!      n.b.  that the dimension of the matrices is ndim >= n
!
!      the above mod is for eigrev which is called by eigset
!      so the length of the eigenvector, n, is the number of collocation
!      points in the radial direction
!
!
      call ax(b,eig,t1,n,ndim)
      l=0
      do 100 i=1,n
      do 100 j=1,n
100   c(i,j)=a(i,j)-lamda*b(i,j)
      call lur(c,n,ndim,ir,ic)
  500 l=l+1
      call solver(t1,c,1,n,ndim,ir,ic)
      call atx(b,adj,t2,n,ndim)
      call solvtrn(t2,c,n,ndim,ir,ic)
       sum1=0._p8
       sum2=0._p8
       do 700 i=1,n
       sum1=sum1 + t1(i)*t1(i)
700    sum2=sum2 + t2(i)*t2(i)
       fac1=1/sqrt(sum1)
       fac2=1/sqrt(sum2)
      do 200 i=1,n
      eig(i)=t1(i)*fac1
  200 adj(i)=t2(i)*fac2
      call ax(a,eig,t2,n,ndim)
      call ax(b,eig,t1,n,ndim)
      if(l.lt.lll)go to 500
      g=(0.0_p8,0.0_p8)
      h=(0.0_p8,0.0_p8)
      do 10 i=1,n
      g = g + adj(i) * t2(i)
   10 h = h + adj(i) * t1(i)
      lamda = g / h
!20    format(* lambda=*,2e17.10)
!20    format(' lambda=',2e17.10)
       sum1=0.0_p8
       sum2=0.0_p8
       do 701 i=1,n
       sum1=sum1 + eig(i)*eig(i)
701    sum2=sum2 + adj(i)*adj(i)
       fac1=1.0_p8/sqrt(sum1)
       fac2=1.0_p8/sqrt(sum2)
      do 702 i=1,n
      eig(i)=eig(i)*fac1
702   adj(i)=adj(i)*fac2
      return
      end
!
!
!
      subroutine eigrevc(a,b,n,ndim,c,lamda,lll,eig,adj)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer:: ndim,n
      complex(kind=p8),dimension(ndim,ndim):: a,b,c
      complex(kind=p8):: fac1,fac2
      complex(kind=p8),dimension(ndim):: eig,adj
      complex(kind=p8),dimension(1000):: t1(1000),t2(1000)
      complex(kind=p8) :: g,h,lamda
      integer,dimension(1000):: ir,ic
      integer :: i,lll,j,l
      real(kind=p8):: sum0,sum1,sum2
      complex(p8):: fac
!
!      this routine assumes ndim <= 256
!
!      if the length of the complex eigenvector, n, is greater
!      than 1000 than you must modify the above.
!
!     note that this modifies eigrevc which is called by linear
!     so n is the length of the eigenvector you are finding in
!     subroutine linear  so in general n = 3 * the number of radial
!     collocations points.  in the convection code it might be
!       4* the number of vertical collocations points.
!
!     (matrix a)*(eigvec) = (lamda)*(matrix b)*(eigvec)
!
!     matrix a nonsingular
!
!     adj : adjoint of eig such that (eig,adj)=1, (eig,other adj)=0
!
      call axc(b,eig,t1,n,ndim)
      l=0
      do 100 i=1,n
      do 100 j=1,n
100   c(i,j)=a(i,j)-lamda*b(i,j)
      call luc(c,n,ndim,ir,ic)
  500 l=l+1
      call solvec(t1,c,1,n,ndim,ir,ic)
      call atxc(b,adj,t2,n,ndim)
      call solvtrc(t2,c,n,ndim,ir,ic)
      sum1 = 0.0_p8
      sum2 = 0.0_p8
      do 199 i = 1,n
      sum1 = sum1 + abs(t1(i))
      sum2 = sum2 + abs(t2(i))
199   continue
      do 200 i = 1,n
      eig(i) = t1(i)/sum1
      adj(i) = t2(i)/sum2
200   continue
      call axc(a,eig,t2,n,ndim)
      call axc(b,eig,t1,n,ndim)
      if(l.lt.lll)go to 500
      g=(0.0_p8,0.0_p8)
      h=(0.0_p8,0.0_p8)
      do 10 i=1,n
      g = g + adj(i) * t2(i)
   10 h = h + adj(i) * t1(i)
      lamda = g / h
!     print 20,  lamda
!     write(6,20) lamda
!20    format(* lambda=*,2e17.10)
 20    format(' lambda=',2e25.18)
      sum0 = 0.0_p8
! ==================================
      eig(1:n)=eig(1:n)/eig(1)
      sum0 = sqrt(sum(real(eig(1:n))**2+aimag(eig(1:n))**2))
      eig(1:n)=eig(1:n)/sum0
      fac = sum(adj(1:n)*eig(1:n))
      adj(1:n) = adj(1:n)/fac 
! ==================================
! ??????????????????????????????????
!       do 80 i = 1,n
! 80    sum = sum + abs(eig(i))
!       do 90 i = 1,n
!       eig(i) = eig(i)/sum
! 90    adj(i) = adj(i)/sum
! ??????????????????????????????????
      return
      end
!Cc
!Cc
!Cc
!C       subroutine eigset(e,ei,ecol,a,nrk,iflg,b,c,wk,s1,s2,
!C     1 eigr,eigi,ir,ic)
!C      implicit none
!C      integer, parameter :: p8=selected_real_kind(p=11)
!C      integer :: nrk,iflg
!C
!Cc         this code assumes that ndim is 65 but
!Cc          calls eigrev 3 times instead of once
!C       real(kind=p8),dimension(65,1):: e,ei,a,b,c
!C       real(kind=p8),dimension(1):: ecol(1)
!Cc
!Cc      the above mod is in eigset which is called by setsuff
!Cc      if the first dimension (not the actual length which is = nrk)
!Cc      of the matrix you are decomposing. if not identically
!Cc      equal to 65 then you must change 65 to ndim=the first dimension
!Cc
!Cc      n.b. the first dimension of the matrix you are decomposing
!Cc      nust be equal to the first dimension of the matrices
!Cc      e and ei into which you are storing the factors
!Cc
!C
!C      real(kind=p8),dimension(1):: s1,s2,wk,eigr,eigi
!C      integer,dimension(1):: ir,ic
!C      integer :: i,j,k,ierr,inc
!C      real(kind=p8) :: zdum
!Cc
!C      do 25 i = 1,nrk
!C      do 20 j = 1,nrk
!C20    b(i,j) = 0.
!C25    b(i,i) = 1.
!Cc      note this next line is diffrent  normally #a#=#a(2,2)#
!C       call eigrg1(65,nrk,a(2,2),0,eigr,eigi,zdum,c,wk,ierr)
!Cc
!Cc      the above mod is in eigset which is called by linear
!Cc      eigrg1 is an ncar subroutine
!Cc
!Cc      if the first dimension of the matrix you are factoring is not
!Cc      equal to 65 than you must change 65 in the above card
!Cc      to be identically euall to the first dimension=ndim
!Cc
!C      do 30 i = 1,nrk
!C      s1(i) = 1.0
!C      s2(i) = 1.0
!C30    continue
!C      inc = 0
!C      do 40 i = 1,nrk
!C      if(iflg .eq. 0) go to 32
!C      if(abs(eigr(i)) .le. 1.0e-08) go to 40
!C32    inc = inc + 1
!C      eigr(i)=eigr(i)+1.0e-9
!Cc      note this next line is diffrent  normally #a#=#a(2,2)#
!C       call eigrev(a(2,2),b,nrk,65,c,eigr(i),5,s1,s2)
!Cc
!Cc      the above mod is in eigset which is called by  setsuff
!Cc
!Cc      if the first dimension of the matrix you are factoring
!Cc      is not equal to 65 change the above card so that
!Cc       65 is replaced by ndim
!Cc
!Cc      note that the 5 that appears in the argument list before
!Cc      s1 and s2 tells how many iterations are done in eigrev
!Cc      one could increase the accuracy by replacing the 5 with
!Cc      6.  one could also call eigrev several times (e.g. 3 times)
!Cc      with updated values of the eigenvalue for better accuracy
!Cc      c.f.  how subroutine linear calls eigrevc 3 times
!Cc
!C      ecol(inc) = eigr(i)
!C      do 35 k = 1,nrk
!C      e(k,inc) = s1(k)
!C      ei(k,inc) = s1(k)
!C35    continue
!C40    continue
!C       call lur(ei,inc,65,ir,ic)
!C       call solver(b,ei,inc,inc,65,ir,ic)
!Cc
!Cc     the above mod is in eigset which is called by setsuff
!Cc
!Cc     if the first dimension of the matrix you are factoring
!Cc     is not identically equal to 65 then replace 65 with ndim
!Cc
!C      do 50 i = 1,inc
!C      do 50 j = 1,inc
!C50    ei(i,j) = b(i,j)
!C       write(6,7676) inc
!Cc7676     format(/,*  the number of rows correctly*,
!Cc     1     * processed in eigset is *,i4,/)
!C 7676     format(/,'  the number of rows correctly',
!C     1     ' processed in eigset is ',i4,/)
!C      return
!C1000  format(1x,e15.7,4x,e15.7)
!C      end
!Cc
!Cc
!Cc
!C      subroutine eigsol(a,linsk,numb,e,ei,ecol,nrk,dfac)
!C      implicit none
!C      integer, parameter :: p8=selected_real_kind(p=11)
!C      integer :: numb,nrk,linsk
!C       real(kind=p8),dimension(1):: a
!C       real(kind=p8),dimension(130):: t
!C       real(kind=p8),dimension(65,1):: e,ei
!C       real(kind=p8),dimension(1)::ecol,dfac
!C       real(kind=p8),dimension(65):: dtfac
!C       integer :: i,j,ii,lbi
!C       real(kind=p8):: dtp
!Cc
!Cc      the above is for eigsol which is called by press and visc
!Cc      t(130) should have dimension >= 2* the number of radial
!Cc      collocation points. dtfac should have dimension >= the
!Cc      number of radial colloation points.
!Cc
!Cc      e and ei must have first dimension the same as in the calling
!Cc      program. i.e.  the same as the first dimension of the
!Cc      matrix you are inverting.
!Cc
!C      do 100 ii = 1,numb
!C      lbi = 2*(ii-1)*linsk
!C      dtp = dfac(ii)
!C      do 10 i = 1,nrk
!C      dtfac(i) = 0.
!C      if(abs(ecol(i)+dtp) .ge. 1.0e-04) dtfac(i) = 1./(ecol(i)+dtp)
!C10    continue
!C      call matmult(t,ei,a(lbi+1),nrk,1,nrk,1,65,1,2,2)
!C      call matmult(t(2),ei,a(lbi+2),nrk,1,nrk,1,65,1,2,2)
!C      do 40 i = 1,nrk
!C40    a(lbi+2*i-1) = t(2*i-1)*dtfac(i)
!C      do 41 i = 1,nrk
!C41    a(lbi+2*i) = t(2*i)*dtfac(i)
!C      call matmult(t,e,a(lbi+1),nrk,1,nrk,1,65,1,2,2)
!C      call matmult(t(2),e,a(lbi+2),nrk,1,nrk,1,65,1,2,2)
!Cc
!Cc      the above 4 mods are in eigsol.
!Cc      the calls to matmult contain the parameter 65 which
!Cc      must be identically equal to the first dimension of
!Cc      the matrix that has been factored and which are inverting.
!Cc
!C      do 70 j = 1,nrk*2
!C70    a(lbi+j) = t(j)
!C100   continue
!C      return
!C      end
!c
!c
!c
      function ff(i,nr)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: i,nr
      real(kind=p8) :: pi,ff

      pi = 4.0_p8*atan(1.0_p8)
      if(i .eq. 0) go to 20
      ff = nr*0.5_p8/tan(pi*i/(2*nr))
      if(mod(i,2) .eq. 0) ff = -ff
      return
20    ff = 0.0_p8
      return
      end
!
!
!
!
      subroutine geneigc(a,b,nn,ndim,reig,aeig,indic)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: nn,ndim
      complex(kind=p8),dimension(ndim,1):: a,b
      real(kind=p8),dimension(1):: reig,aeig
      real(kind=p8),dimension(2000):: wk2
      integer :: i,j,k,k1,kp,km1,jp,ip,ierr,kk,m,nm1,l,n
!
!      this routine assumes ndim <= 256
!      the above mod is in geneigc which is called by linear
!      wk2 is workspace passed on to eigcg1 via subroutine shuffle
!      eigcg1 is an ncar subroutine that requires that wk2
!      be of dimension (2*n + 2)
!      so if the dimension of the eigenvector that you are finding
!      in subroutine linear is greater than n=999 than you must
!      change 2000 to an integer >= (2*n +2)
!
!   a and b are matrices
!   this program solves the generalized eiginvalue problem
!   a = (eig)b
!   indic is an array which contains informations on the success
!
!   reig : real part of the eigenvalues
!   aeig : imaginary part of the eigenvalues
!
!   of finding the kth eigenvalue
      complex(kind=p8),dimension(1):: eig
      complex(kind=p8):: c,d1,det
      integer,dimension(1):: indic
      real(kind=p8):: xmax,y,eps
      real(kind=p8),parameter :: epss=1.0e-12_p8
      n=nn
      d1=cmplx(1.0_p8,0.0_p8,p8)
1     nm1=n-1
!   reduce b to upper triangular form using gaussian elimination
!   with complete pivoting
      do 101 k=1,nm1
      l=k
      m=k
      xmax=abs(b(l,m))
      do 102 i=k,n
      do 102 j=k,n
      y=abs(b(i,j))
      if(xmax.ge.y) go to 102
      l=i
      m=j
      xmax=y
102   continue
      if(k.eq.1)eps=epss*xmax
      if(xmax.lt.eps)go to 200
      if(l.eq.k)go to 105
!   row operations
      do 103 i=1,n
      c=a(k,i)
      a(k,i)=a(l,i)
103   a(l,i)=c
      do 104 i=k,n
      c=b(k,i)
      b(k,i)=b(l,i)
104   b(l,i)=c
105   if(m.eq.k)go to 107
!   column operations
      do 106 i=1,n
      c=a(i,k)
      a(i,k)=a(i,m)
      a(i,m)=c
      c=b(i,k)
      b(i,k)=b(i,m)
106   b(i,m)=c
107   det=d1/b(k,k)
      kp=k+1
      do 108 i=kp,n
      c=b(i,k)*det
      do 109 j=1,n
109   a(i,j)=a(i,j)-c*a(k,j)
      do 110 j=k,n
110   b(i,j)=b(i,j)-c*b(k,j)
108   continue
101   continue
!   b is now triangularized
      if(abs(b(n,n)).gt.eps)go to 300
      k=n
!   triangularize the lower part of a
200   km1=k-1
      do 201 kk=k,n
      k1=n+k-kk
      l=k1
      m=k1
      xmax=abs(a(l,m))
      do 202 i=k,k1
      do 202 j=1,k1
      y=abs(a(i,j))
      if(xmax.ge.y)go to 202
      l=i
      m=j
      xmax=y
202   continue
      if(xmax.lt.eps)go to 250
      if(m.eq.k1)go to 205
      do 203 i=1,k1
      c=a(i,k1)
      a(i,k1)=a(i,m)
203   a(i,m)=c
      do 204 i=1,km1
      c=b(i,k1)
      b(i,k1)=b(i,m)
204   b(i,m)=c
205   if(l.eq.k1)go to 207
      do 206 i=1,n
      c=a(k1,i)
      a(k1,i)=a(l,i)
206   a(l,i)=c
207   det=d1/a(k1,k1)
      do 208 j=1,k1
      c=a(k1,j)*det
      do 209 i=1,k1
      a(i,j)=a(i,j)-c*a(i,k1)
209   b(i,j)=b(i,j)-c*b(i,k1)
208   continue
201   continue
210   n=k-1
      go to 1
300   c=d1/b(n,n)
      do 301 k=1,n
301   a(n,k)=a(n,k)*c
      do 302 i=1,nm1
      j=n-i
      c=d1/b(j,j)
      jp=j+1
      do 303 k=1,n
      det=a(j,k)
      do 304 ip=jp,n
304   det=det-b(j,ip)*a(ip,k)
303   a(j,k)=c*det
302   continue
      call shuffle(a,b,n,ndim,reig,aeig,wk2,ierr)
      return
250   print 251,xmax
!251   format(* error--det(a-(eig)b)=0 for all values of eig--xmax=*
!     1e20.8)
251   format(' error--det(a-(eig)b)=0 for all values of eig--xmax=',e20.8)
      return
      end
!
!
!
      subroutine lur(a,n,ndim,ir,ic)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: n,ndim
      real(kind=p8),dimension(ndim,1):: a
      integer,dimension(1):: ir, ic
      integer :: i,j,k,l,m,irl,icm,k1
      real(kind=p8) :: xmax,y,b,c

      do 10 i = 1,n
      ir(i) = i
 10   ic(i) = i
      k = 1
      l = k
      m = k
      xmax = abs(a(k,k))
      do 100 i = k,n
      do 100 j = k,n
      y = abs(a(i,j))
      if( xmax .ge. y)  go to 100
      xmax = y
      l = i
      m = j
 100  continue
      do 1000 k = 1,n
      irl = ir(l)
      ir(l) = ir(k)
      ir(k) = irl
      icm = ic(m)
      ic(m) = ic(k)
      ic(k) = icm
      if( l .eq. k )  go to 300
      do 200 j = 1,n
      b = a(k,j)
      a(k,j) = a(l,j)
 200  a(l,j) = b
 300  if( m .eq. k )  go to 500
      do 400 i = 1,n
      b = a(i,k)
      a(i,k) = a(i,m)
 400  a(i,m) = b
 500  c = 1._p8/a(k,k)
      a(k,k) = c
      if( k .eq. n)  go to 1000
      k1 = k + 1
      xmax = abs(a(k1,k1))
      l = k1
      m = k1
      do 600 i = k1,n
 600  a(i,k) = c*a(i,k)
      do 800 i = k1,n
      b = a(i,k)
      do 800 j = k1,n
      a(i,j) = a(i,j) - b*a(k,j)
      y = abs(a(i,j))
      if( xmax .ge. y) go to 800
      xmax = y
      l = i
      m = j
 800  continue
 1000 continue
      return
      end
!
!
!
      subroutine luc(a,n,ndim,ir,ic)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
!
!   program to perform fully pivoted lu decomposition of general
!   complex array a
!
      integer :: n,ndim
      complex(kind=p8),dimension(ndim,1):: a
      complex(kind=p8):: b, c
      integer,dimension(1):: ir, ic
      integer i,j,k,l,m,irl,icm,k1
      real(kind=p8) :: xmax,y

      do 10 i = 1,n
      ir(i) = i
 10   ic(i) = i
      k = 1
      l = k
      m = k
      xmax = abs(real(a(k,k))) + abs(aimag(a(k,k)))
      do 100 i = k,n
      do 100 j = k,n
      y = abs(real(a(i,j))) + abs(aimag(a(i,j)))
      if( xmax .ge. y)  go to 100
      xmax = y
      l = i
      m = j
 100  continue
      do 1000 k = 1,n
      irl = ir(l)
      ir(l) = ir(k)
      ir(k) = irl
      icm = ic(m)
      ic(m) = ic(k)
      ic(k) = icm
      if( l .eq. k )  go to 300
      do 200 j = 1,n
      b = a(k,j)
      a(k,j) = a(l,j)
 200  a(l,j) = b
 300  if( m .eq. k )  go to 500
      do 400 i = 1,n
      b = a(i,k)
      a(i,k) = a(i,m)
 400  a(i,m) = b
 500  c = 1._p8/a(k,k)
      a(k,k) = c
      if( k .eq. n)  go to 1000
      k1 = k + 1
      xmax = abs(real(a(k1,k1))) + abs(aimag(a(k1,k1)))
      l = k1
      m = k1
      do 600 i = k1,n
 600  a(i,k) = c*a(i,k)
      do 800 i = k1,n
      b = a(i,k)
      do 800 j = k1,n
      a(i,j) = a(i,j) - b*a(k,j)
      y = abs(real(a(i,j))) + abs(aimag(a(i,j)))
      if( xmax .ge. y) go to 800
      xmax = y
      l = i
      m = j
 800  continue
 1000 continue
      return
      end

!Cc
!Cc
!Cc
!C      subroutine neigsol(a,linsk,numb,e,ei,ecol,nrk,dfac,ndim)
!C      implicit none
!C      integer, parameter :: p8=selected_real_kind(p=11)
!C      integer :: linsk,numb,nrk,ndim
!C
!Cc     this routine assumes ndim <= 256
!Cc     1000 is greater than 2*ndim
!Cc
!C      real(kind=p8),dimension(1):: a
!C      real(kind=p8),dimension(1000):: t
!Cc     dimension e(ndim,1),ei(ndim,1),ecol(1),dfac(1),dtfac(ndim)
!C      real(kind=p8),dimension(ndim,1):: e,ei
!C      real(kind=p8),dimension(1):: ecol,dfac
!C      real(kind=p8),dimension(1024):: dtfac
!C      integer :: ii,i,j,lbi
!C      real(kind=p8):: dtp
!C
!C      do 100 ii = 1,numb
!C         lbi = 2*(ii-1)*linsk
!C         dtp = dfac(ii)
!C
!Cc        if(nrk.gt.100) then
!Cc           write(6,*) 'error in neigsol - exiting'
!Cc           stop
!Cc        endif
!C         do 10 i = 1,nrk
!C            dtfac(i) = 0.
!C10          if(abs(ecol(i)+dtp).ge.1.0e-04) dtfac(i)=1./(ecol(i)+dtp)
!C
!C         call mxva(ei,1,ndim,a(lbi+1),2,t(1),2,nrk,nrk)
!C         call mxva(ei,1,ndim,a(lbi+2),2,t(2),2,nrk,nrk)
!C
!C         do 40 i = 1,nrk
!C40          a(lbi+2*i-1) = t(2*i-1)*dtfac(i)
!C
!C         do 41 i = 1,nrk
!C41          a(lbi+2*i) = t(2*i)*dtfac(i)
!C
!C         call mxva(e,1,ndim,a(lbi+1),2,t(1),2,nrk,nrk)
!C         call mxva(e,1,ndim,a(lbi+2),2,t(2),2,nrk,nrk)
!C
!C         do 70 j = 1,nrk*2
!C70          a(lbi+j) = t(j)
!C
!C100   continue
!C
!C      return
!C      end
!
!
!
      subroutine shuffle(a,b,n,ndim,reig,aeig,wk2,ierr)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: n,ndim,ierr
      complex(kind=p8),dimension(ndim,1):: a
      real(kind=p8),dimension(ndim,1):: b
      integer i,j
      real(kind=p8),dimension(1) :: reig,aeig,wk2
      real(kind=p8) :: zdum1,zdum2

      do 10 i = 1,n
      do 10 j = 1,n
      b(i,j) = real(a(i,j))
10    b(i,j+n) = aimag(a(i,j))
      call eigcg1(ndim,n,b(1,1),b(1,1+n),0,reig,aeig,zdum1,zdum2,a,wk2,ierr)
      return
      end
!
!
!
      subroutine solver(f,a,k,n,ndim,ir,ic)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: n,ndim,k
      real(kind=p8),dimension(ndim,1):: a,f
      real(kind=p8),dimension(ndim):: g
      integer,dimension(1):: ir, ic
!
!       the above mod is for solve which is called by eigrev which
!       is called by eigeset. soleve is also called directly by eigset.
!       if n= the length of the eigenvector is
!       greater than 100 change 100 to n.  in general n is the
!       number of radial collocation points.
!
!
      integer :: ni,i,j,kk,ici,i1,it,n1,iri
      real(p8) :: b

      n1 = n + 1
      do 1000 kk = 1,k
      do 100 i = 1,n
      iri = ir(i)
 100  g(i) = f(iri,kk)
      do 400 i = 2,n
      i1 = i - 1
      b = g(i)
      do 300 j = 1,i1
 300  b = b - a(i,j)*g(j)
 400  g(i) = b
      do 700 it = 1,n
      i = n1 - it
      i1 = i + 1
      b = g(i)
      if( i .eq. n)  go to 700
      do 600 j = i1,n
 600  b = b - a(i,j)*g(j)
 700  g(i) = b*a(i,i)
      do 900 i = 1,n
      ici = ic(i)
 900  f(ici,kk) = g(i)
 1000 continue
      return
      end
!
!
!
      subroutine solvec(f,a,k,n,ndim,ir,ic)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: k,n,ndim
      integer,dimension(1):: ir, ic
      complex(p8),dimension(ndim,1):: a,f
      complex(p8),dimension(1024):: g
      complex(p8):: b,c
      integer :: n1,kk,iri,i1,j,it,ici,i
!
!   program to solve k equations ax = f (dimensioned ndim)
!   assume previous call to routine lu has been made
!
!       this routine assumes ndim <= 256
!       this is used in solvec which is called by eigrevc which
!       is called by linear
!
!      if the length, n, of the eigenvector beign found in linear
!      is greater than 1000 than replace 1000 in the above card with n
!
      n1 = n + 1
      do 1000 kk = 1,k
      do 100 i = 1,n
      iri = ir(i)
 100  g(i) = f(iri,kk)
      do 400 i = 2,n
      i1 = i - 1
      b = g(i)
      do 300 j = 1,i1
 300  b = b - a(i,j)*g(j)
 400  g(i) = b
      do 700 it = 1,n
      i = n1 - it
      i1 = i + 1
      b = g(i)
      if( i .eq. n)  go to 700
      do 600 j = i1,n
 600  b = b - a(i,j)*g(j)
 700  g(i) = b*a(i,i)
      do 900 i = 1,n
      ici = ic(i)
 900  f(ici,kk) = g(i)
 1000 continue
      return
      end
!
!
!
      subroutine solvtrc(f,a,n,ndim,ir,ic)
      implicit none
      integer n,ndim
      integer, parameter :: p8=selected_real_kind(p=11)
      complex(p8),dimension(1):: f
      complex(p8),dimension(ndim,1):: a
      integer,dimension(1):: ir,ic
      complex(p8),dimension(1000):: g
      complex(p8):: b,c
      integer n1,i,j,ici,i1,iri,it
!
!       this routine assumes ndim <= 256
!
!       this is used in solvtrc which is called by eigrevc which
!       is called by linear
!
!      if the length, n, of the eigenvector beign found in linear
!      is greater than 1000 than replace 1000 in the above card with n
!
      n1=n+1
      do 100 i=1,n
      ici=ic(i)
100   g(i)=f(ici)
      g(1)=a(1,1)*g(1)
      do 200 i=2,n
      i1=i-1
      b=g(i)
      do 300 j=1,i1
300   b=b-a(j,i)*g(j)
200   g(i)=b*a(i,i)
      do 400 it=2,n
      i=n1-it
      i1=i+1
      b=g(i)
      if(i.eq.n) go to 400
      do 500 j=i1,n
500   b=b-a(j,i)*g(j)
400   g(i)=b
      do 600 i=1,n
      iri=ir(i)
600   f(iri)=g(i)
      return
      end
!
!
!
      subroutine solvtrn(f,a,n,ndim,ir,ic)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: n,ndim
      real(p8),dimension(ndim,1):: a
      real(p8),dimension(1):: f
      real(p8),dimension(ndim):: g
      integer,dimension(1):: ir,ic
      integer :: n1,i,j,it,i1,ici,iri
      real(p8):: b
!
!       the above mod is for solvtrn which is called by eigrev which
!       is called by eigeset. if n=the lensght of the eigenvector is
!       greater than 100 change 100 to n.  in general n is the
!       number of radial collocation points.
!
!
      n1=n+1
      do 100 i=1,n
      ici=ic(i)
100   g(i)=f(ici)
      g(1)=a(1,1)*g(1)
      do 200 i=2,n
      i1=i-1
      b=g(i)
      do 300 j=1,i1
300   b=b-a(j,i)*g(j)
200   g(i)=b*a(i,i)
      do 400 it=2,n
      i=n1-it
      i1=i+1
      b=g(i)
      if(i.eq.n) go to 400
      do 500 j=i1,n
500   b=b-a(j,i)*g(j)
400   g(i)=b
      do 600 i=1,n
      iri=ir(i)
600   f(iri)=g(i)
      return
      end
!
!
!
      function xnum(kr,jr,nr)
      implicit none
      integer, parameter :: p8=selected_real_kind(p=11)
      integer :: kr,jr,nr
      real(p8) :: pi,xnum,ff

      pi = 4._p8*atan(1._p8)
      if(kr .eq. 0) go to 20
      xnum = ff(kr+jr,nr) + ff(kr-jr,nr)
      return
20    if(jr .eq. 0 .or. jr .eq. nr) go to 30
      xnum = .5_p8*(nr+.5_p8)*( &
!      (nr+.5_p8) + (cot(pi*(jr)/(2*nr)))**2) &
       (nr+.5_p8) + (1/tan(pi*jr/(2*nr)))**2) &
       + 1/8._p8 - .25_p8/(sin(pi*jr/(2*nr))**2) &
       -.5_p8*(nr*nr)
      return
30    xnum = .5_p8
      if( mod(nr,2) .eq. 1 ) xnum = - xnum
      if(jr.eq.0) xnum = (1/3._p8)*(nr*nr) + 1/6._p8
      return
      end
