  Module Matrix_Mod
    use CFML_GlobalDeps,       only: sp,dp,cp,tpi
    use CFML_Math_General,     only: sort, trace, iminloc, SVDcmp
    use CFML_Rational_Arithmetic_test

    implicit none
    private
    public :: Identity_Matrix, Diagonal_Matrix, Vectorize, row_echelon, Kronecker_Product, Sylvester_Solver

    real(kind=dp), parameter :: epsd=1.0e-200_dp
    real(kind=dp), parameter :: epss=1.0e-30_dp

    interface Vectorize
    	module procedure Vectorize_int
    	module procedure Vectorize_sp
    	module procedure Vectorize_dp
    	module procedure Vectorize_complex
    	module procedure Vectorize_rational
    end interface Vectorize

    interface Kronecker_Product
    	module procedure Kronecker_Product_int
    	module procedure Kronecker_Product_sp
    	module procedure Kronecker_Product_dp
    	module procedure Kronecker_Product_complex
    	module procedure Kronecker_Product_rational
    end interface Kronecker_Product

    interface Identity_Matrix
    	module procedure Identity_int
    	module procedure Identity_sp
    	module procedure Identity_dp
    	module procedure Identity_complex
    	module procedure Identity_rational
    end interface Identity_Matrix

    interface Diagonal_Matrix
    	module procedure diagonal_int
    	module procedure diagonal_sp
    	module procedure diagonal_dp
    	module procedure diagonal_complex
    	module procedure diagonal_rational
    end interface Diagonal_Matrix

    interface row_echelon
    	module procedure row_echelon_sp
    	module procedure row_echelon_dp
    	module procedure row_echelon_rational
    end interface row_echelon

    interface Sylvester_Solver
    	module procedure Sylvester_Solver_sp
    	module procedure Sylvester_Solver_dp
    	module procedure Sylvester_Solver_rational
    end interface Sylvester_Solver

  Contains

    !!---  Set of procedures to be moved into Rational_Arithmetic

    Pure function diagonal_complex(v) result (vMat)
    	complex(kind=dp),dimension(:),    intent(in) :: v
    	complex(kind=dp),dimension(size(v),size(v))  :: vMat
    	integer :: n,i
    	n=size(v)
      vMat= (0.0_dp,0.0_dp)
      do i=1,n
      	vMat(i,i)=v(i)
      end do
    end function diagonal_complex

    Pure function diagonal_dp(v) result (vMat)
    	real(kind=dp),dimension(:),    intent(in) :: v
    	real(kind=dp),dimension(size(v),size(v))  :: vMat
    	integer :: n,i
    	n=size(v)
      vMat= 0.0_dp
      do i=1,n
      	vMat(i,i)=v(i)
      end do
    end function diagonal_dp

    Pure function diagonal_int(v) result (vMat)
    	integer,dimension(:),    intent(in) :: v
    	integer,dimension(size(v),size(v))  :: vMat
    	integer :: n,i
    	n=size(v)
      vMat= 0
      do i=1,n
      	vMat(i,i)=v(i)
      end do
    end function diagonal_int

    Pure function diagonal_sp(v) result (vMat)
    	real(kind=sp),dimension(:),    intent(in) :: v
    	real(kind=sp),dimension(size(v),size(v))  :: vMat
    	integer :: n,i
    	n=size(v)
      vMat= 0.0_sp
      do i=1,n
      	vMat(i,i)=v(i)
      end do
    end function diagonal_sp

    Pure function diagonal_rational(v) result (vMat)
    	type(rational),dimension(:),    intent(in) :: v
    	type(rational),dimension(size(v),size(v))  :: vMat
    	integer :: n,i
    	n=size(v)
      vMat= 0_ik//1_ik
      do i=1,n
      	vMat(i,i)=v(i)
      end do
    end function diagonal_rational

    Pure function vectorize_complex(A) result (vA)
    	complex(kind=dp),dimension(:,:),intent(in) :: A
    	complex(kind=dp),dimension(size(A,dim=1)*size(A,dim=2))  :: vA
    	integer :: m,n,i,j,k
    	n=size(A,dim=1) ; m=size(A,dim=2)
    	k=1
    	do j=1,m
          vA(k:k+n-1) = A(:,j)
          k=k+n
      end do
    end function vectorize_complex

    Pure function vectorize_dp(A) result (vA)
    	real(kind=dp),dimension(:,:),intent(in) :: A
    	real(kind=dp),dimension(size(A,dim=1)*size(A,dim=2))  :: vA
    	integer :: m,n,i,j,k
    	n=size(A,dim=1) ; m=size(A,dim=2)
    	k=1
    	do j=1,m
          vA(k:k+n-1) = A(:,j)
          k=k+n
      end do
    end function vectorize_dp

    Pure function vectorize_int(A) result (vA)
    	integer,dimension(:,:),intent(in) :: A
    	integer,dimension(size(A,dim=1)*size(A,dim=2))  :: vA
    	integer :: m,n,i,j,k
    	n=size(A,dim=1) ; m=size(A,dim=2)
    	k=1
    	do j=1,m
          vA(k:k+n-1) = A(:,j)
          k=k+n
      end do
    end function vectorize_int

    Pure function vectorize_rational(A) result (vA)
    	type(rational),dimension(:,:),intent(in) :: A
    	type(rational),dimension(size(A,dim=1)*size(A,dim=2))  :: vA
    	integer :: m,n,i,j,k
    	n=size(A,dim=1) ; m=size(A,dim=2)
    	k=1
    	do j=1,m
          vA(k:k+n-1) = A(:,j)
          k=k+n
      end do
    end function vectorize_rational

    Pure function vectorize_sp(A) result (vA)
    	real(kind=sp),dimension(:,:),intent(in) :: A
    	real(kind=sp),dimension(size(A,dim=1)*size(A,dim=2))  :: vA
    	integer :: m,n,i,j,k
    	n=size(A,dim=1) ; m=size(A,dim=2)
    	k=1
    	do j=1,m
          vA(k:k+n-1) = A(:,j)
          k=k+n
      end do
    end function vectorize_sp

    Pure function Kronecker_Product_complex(A,B) result (C)
    	complex(kind=dp),dimension(:,:),intent(in) :: A,B
    	complex(kind=dp),dimension(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))  :: C
    	integer :: m,n,p,q,i,j,k,l,ip,jp
    	n=size(A,dim=1) ; m=size(A,dim=2) ; p=size(B,dim=1) ; q=size(B,dim=2)
      forall ( i = 1:n , j = 1:m ) C(p*(i-1)+1 : p*i , q*(j-1)+1 : q*j) = A(i,j)*B
    end function Kronecker_Product_complex

    Pure function Kronecker_Product_dp(A,B) result (C)
    	real(kind=dp),dimension(:,:),intent(in) :: A,B
    	real(kind=dp),dimension(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))  :: C
    	integer :: m,n,p,q,i,j,k,l,ip,jp
    	n=size(A,dim=1) ; m=size(A,dim=2) ; p=size(B,dim=1) ; q=size(B,dim=2)
      forall ( i = 1:n , j = 1:m ) C(p*(i-1)+1 : p*i , q*(j-1)+1 : q*j) = A(i,j)*B
    end function Kronecker_Product_dp

    Pure function Kronecker_Product_int(A,B) result (C)
    	integer,dimension(:,:),intent(in) :: A,B
    	integer,dimension(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))  :: C
    	integer :: m,n,p,q,i,j,k,l,ip,jp
    	n=size(A,dim=1) ; m=size(A,dim=2) ; p=size(B,dim=1) ; q=size(B,dim=2)
      forall ( i = 1:n , j = 1:m ) C(p*(i-1)+1 : p*i , q*(j-1)+1 : q*j) = A(i,j)*B
    end function Kronecker_Product_int

    Pure function Kronecker_Product_sp(A,B) result (C)
    	real(kind=sp),dimension(:,:),intent(in) :: A,B
    	real(kind=sp),dimension(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))  :: C
    	integer :: m,n,p,q,i,j,k,l,ip,jp
    	n=size(A,dim=1) ; m=size(A,dim=2) ; p=size(B,dim=1) ; q=size(B,dim=2)
      forall ( i = 1:n , j = 1:m ) C(p*(i-1)+1 : p*i , q*(j-1)+1 : q*j) = A(i,j)*B
    end function Kronecker_Product_sp

    Pure function Kronecker_Product_rational(A,B) result (C)
    	type(rational),dimension(:,:),intent(in) :: A,B
    	type(rational),dimension(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))  :: C
    	integer :: m,n,p,q,i,j,k,l,ip,jp
    	n=size(A,dim=1) ; m=size(A,dim=2) ; p=size(B,dim=1) ; q=size(B,dim=2)
      forall ( i = 1:n , j = 1:m ) C(p*(i-1)+1 : p*i , q*(j-1)+1 : q*j) = A(i,j)*B
    end function Kronecker_Product_rational

    subroutine identity_complex(d, identity)                       ! Returns the complex dxd identity matrix
      integer,                             intent(in) :: d         ! Dimension of the identity matrix
      complex(kind=dp),dimension(1:d,1:d),intent(out) :: identity  ! Identity matrix
      integer :: j, k
      forall ( j = 1:d, k = 1:d, j /= k ) identity(j,k) = (0.0_sp,0.0_sp)
      forall ( j = 1:d, k = 1:d, j == k ) identity(j,k) = (1.0_sp,0.0_sp)
    end subroutine identity_complex

    subroutine identity_sp(d, identity)
      integer,                            intent(in) :: d         ! Dimension of the identity matrix
      real(kind=sp), dimension(1:d,1:d), intent(out) :: identity  ! Identity matrix
      integer :: j, k
      forall ( j = 1:d, k = 1:d, j /= k ) identity(j,k) = 0.0_sp
      forall ( j = 1:d, k = 1:d, j == k ) identity(j,k) = 1.0_sp
    end subroutine identity_sp

    subroutine identity_dp(d, identity)
      integer,                            intent(in) :: d         ! Dimension of the identity matrix
      real(kind=dp), dimension(1:d,1:d), intent(out) :: identity  ! Identity matrix
      integer :: j, k
      forall ( j = 1:d, k = 1:d, j /= k ) identity(j,k) = 0.0_dp
      forall ( j = 1:d, k = 1:d, j == k ) identity(j,k) = 1.0_dp
    end subroutine identity_dp

    subroutine identity_int(d, identity)
      integer,                      intent(in) :: d         ! Dimension of the identity matrix
      integer, dimension(1:d,1:d), intent(out) :: identity  ! Identity matrix
      integer :: j, k
      forall ( j = 1:d, k = 1:d, j /= k ) identity(j,k) = 0
      forall ( j = 1:d, k = 1:d, j == k ) identity(j,k) = 1
    end subroutine identity_int

    subroutine identity_rational(d, identity)
      integer,                             intent(in) :: d         ! Dimension of the identity matrix
      type(rational), dimension(1:d,1:d), intent(out) :: identity  ! Identity matrix
      integer :: j, k
      forall ( j = 1:d, k = 1:d, j /= k ) identity(j,k) = 0_ik//1_ik
      forall ( j = 1:d, k = 1:d, j == k ) identity(j,k) = 1_ik//1_ik
    end subroutine identity_rational

    !!Solve the Sylvester square matrix equation AX-XB=C to look for transformation matrices relating a set of Symmetry Operators
    !!The matrices A(n x n) and B(m x m) are square but they can be of different dimension. The matrix X(n x m) has the same
    !!shape as C(n x m)
    !!---- The resolution of the Sylvester equation for square matrices pass by the conversion to a linear equation
    !!----  (Im .kron. A - Bt .kron. In) vec(X) = vec(C)
    !!----
    Subroutine Sylvester_Solver_sp(A,B,C,X,nsol,ok)
    	real(kind=sp), dimension(:,:), intent(in) :: A,B,C
    	real(kind=sp), dimension(:,:), intent(out):: X !The declared dimensions of X and C should be at least (n*n,m)
    	integer,                       intent(out):: nsol
    	logical,                       intent(out):: ok
    	!--- Local variables ---!
    	real(kind=sp), dimension(size(A,dim=1)*size(B,dim=1),size(A,dim=1)*size(B,dim=1)) :: U,mw,V
    	real(kind=sp), dimension(size(A,dim=1)*size(B,dim=1)) :: vX,vC,w
    	real(kind=sp), dimension(size(A,dim=1),size(A,dim=1)) :: Idn
    	real(kind=sp), dimension(size(B,dim=1),size(B,dim=1)) :: Idm

    	integer :: i,j,k,n,m
    	logical :: homogeneous,singular

    	n=size(A,dim=1) ; m=size(B,dim=1); nsol=0
    	ok=.true.
    	call Identity_Matrix(n,Idn)
    	call Identity_Matrix(m,Idm)
    	homogeneous=.false.
    	if(sum(abs(C)) < epss) homogeneous=.true.
    	U=Kronecker_Product(Idm,A)-Kronecker_Product(transpose(B),Idn)
    	vC=vectorize(C)
      X=0.0_sp
      call SVDcmp(U,w,V)
      singular=.false.
      if(any(abs(w) < epss)) singular=.true.

      if(homogeneous .and. .not. singular) then

      	nsol=1
      	ok=.false.  ! The only solution is X=0.0_sp

      else if(homogeneous .and. singular) then !We store in the columns of X directly the
        nsol=count(abs(w) <= epss)
        nsol=0
        do i=1,n
        	if(abs(w(i)) <= epss) then
        		 nsol=nsol+1
        		 X(:,nsol)=V(:,i)
        	end if
        end do

      else if(.not. homogeneous .and. singular) then

      	nsol=0
      	ok=.false.  !no solution

      else if(.not. homogeneous .and. .not. singular) then

      	where(abs(w) > epss )
      		w=1/w
      	elsewhere
      		w=0
      	end where
      	nsol=1
      	vX=matmul(matmul(matmul(V,diagonal_matrix(w)),transpose(U)),vC)
      	X(1:n,1:m)=reshape(vX,[n,m])

      end if

    End Subroutine Sylvester_Solver_sp

    Subroutine Sylvester_Solver_dp(A,B,C,X,nsol,ok)
    	real(kind=dp), dimension(:,:), intent(in) :: A,B,C
    	real(kind=dp), dimension(:,:), intent(out):: X !The declared dimensions of X and C should be at least (n*n,m)
    	integer,                       intent(out):: nsol
    	logical,                       intent(out):: ok
    	!--- Local variables ---!
    	real(kind=dp), dimension(size(A,dim=1)*size(B,dim=1),size(A,dim=1)*size(B,dim=1)) :: U,mw,V
    	real(kind=dp), dimension(size(A,dim=1)*size(B,dim=1)) :: vX,vC,w
    	real(kind=dp), dimension(size(A,dim=1),size(A,dim=1)) :: Idn
    	real(kind=dp), dimension(size(B,dim=1),size(B,dim=1)) :: Idm

    	integer :: i,j,k,n,m
    	logical :: homogeneous,singular

    	n=size(A,dim=1) ; m=size(B,dim=1); nsol=0
    	ok=.true.
    	call Identity_Matrix(n,Idn)
    	call Identity_Matrix(m,Idm)
    	homogeneous=.false.
    	if(sum(abs(C)) < epsd) homogeneous=.true.
    	U=Kronecker_Product(Idm,A)-Kronecker_Product(transpose(B),Idn)
    	vC=vectorize(C)
      X=0.0_sp
      call SVDcmp(U,w,V)
      singular=.false.
      if(any(abs(w) < epsd)) singular=.true.

      if(homogeneous .and. .not. singular) then

      	nsol=1
      	ok=.false.  ! The only solution is X=0.0_sp

      else if(homogeneous .and. singular) then !We store in the columns of X directly the
        nsol=count(abs(w) <= epsd)
        nsol=0
        do i=1,n
        	if(abs(w(i)) <= epsd) then
        		 nsol=nsol+1
        		 X(:,nsol)=V(:,i)
        	end if
        end do

      else if(.not. homogeneous .and. singular) then

      	nsol=0
      	ok=.false.  !no solution

      else if(.not. homogeneous .and. .not. singular) then

      	where(abs(w) > epsd )
      		w=1/w
      	elsewhere
      		w=0
      	end where
      	nsol=1
      	vX=matmul(matmul(matmul(V,diagonal_matrix(w)),transpose(U)),vC)
      	X(1:n,1:m)=reshape(vX,[n,m])

      end if

    End Subroutine Sylvester_Solver_dp

    Subroutine Sylvester_Solver_rational(A,B,C,X,nsol,ok)
    	type(rational), dimension(:,:), intent(in) :: A,B,C
    	type(rational), dimension(:,:), intent(out):: X !The declared dimensions of X and C should be at least (n*n,m)
    	integer,                        intent(out):: nsol
    	logical,                        intent(out):: ok

    	!--- Local variables ---!
    	real(kind=dp), dimension(size(A,dim=1),size(A,dim=2)):: Ar !The real version of A
    	real(kind=dp), dimension(size(B,dim=1),size(B,dim=2)):: Br !The real version of B
    	real(kind=dp), dimension(size(C,dim=1),size(C,dim=2)):: Cr !The real version of C
    	real(kind=dp), dimension(size(X,dim=1),size(X,dim=2)):: Xr !The real version of X

    	Ar=A; Br=B; Cr=C  !converting to double precision the rational matrices
      call Sylvester_Solver_dp(Ar,Br,Cr,Xr,nsol,ok)
      X=Xr
    End Subroutine Sylvester_Solver_rational

    Subroutine row_echelon_sp(matrix)
      real(kind=sp), dimension(:,:), intent(in out) :: matrix
      !--- Local variables ---!
      integer :: pivot, n_row, n_column
      integer :: r, i
      real(kind=sp), dimension(:), allocatable :: trow

      pivot = 1
      n_row = size(matrix, 1)
      n_column = size(matrix, 2)

      allocate(trow(n_column))

      do r = 1, n_row
         if ( n_column <= pivot ) exit
         i = r
         do while ( matrix(i, pivot) == 0 )
            i = i + 1
            if ( n_row == i ) then
               i = r
               pivot = pivot + 1
               if ( n_column == pivot ) return
            end if
         end do
         trow = matrix(i, :)
         matrix(i, :) = matrix(r, :)
         matrix(r, :) = trow
         matrix(r, :) = matrix(r, :) / matrix(r, pivot)
         do i = 1, n_row
            if ( i /= r ) matrix(i, :) = matrix(i, :) - matrix(r, :) * matrix(i, pivot)
         end do
         pivot = pivot + 1
      end do
    End Subroutine row_echelon_sp

    Subroutine row_echelon_dp(matrix)
      real(kind=dp), dimension(:,:), intent(in out) :: matrix
      !--- Local variables ---!
      integer :: pivot, n_row, n_column
      integer :: r, i
      real(kind=dp), dimension(:), allocatable :: trow

      pivot = 1
      n_row = size(matrix, 1)
      n_column = size(matrix, 2)

      allocate(trow(n_column))

      do r = 1, n_row
         if ( n_column <= pivot ) exit
         i = r
         do while ( matrix(i, pivot) == 0 )
            i = i + 1
            if ( n_row == i ) then
               i = r
               pivot = pivot + 1
               if ( n_column == pivot ) return
            end if
         end do
         trow = matrix(i, :)
         matrix(i, :) = matrix(r, :)
         matrix(r, :) = trow
         matrix(r, :) = matrix(r, :) / matrix(r, pivot)
         do i = 1, n_row
            if ( i /= r ) matrix(i, :) = matrix(i, :) - matrix(r, :) * matrix(i, pivot)
         end do
         pivot = pivot + 1
      end do
    End Subroutine row_echelon_dp

    Subroutine row_echelon_rational(matrix)
      type(rational), dimension(:,:), intent(in out) :: matrix
      !--- Local variables ---!
      integer :: pivot, n_row, n_column
      integer :: r, i
      type(rational), dimension(:), allocatable :: trow

      pivot = 1
      n_row = size(matrix, 1)
      n_column = size(matrix, 2)

      allocate(trow(n_column))

      do r = 1, n_row
         if ( n_column <= pivot ) exit
         i = r
         do while ( matrix(i, pivot) == 0_ik//1_ik )
            i = i + 1
            if ( n_row == i ) then
               i = r
               pivot = pivot + 1
               if ( n_column == pivot ) return
            end if
         end do
         trow = matrix(i, :)
         matrix(i, :) = matrix(r, :)
         matrix(r, :) = trow
         matrix(r, :) = matrix(r, :) / matrix(r, pivot)
         do i = 1, n_row
            if ( i /= r ) matrix(i, :) = matrix(i, :) - matrix(r, :) * matrix(i, pivot)
         end do
         pivot = pivot + 1
      end do
    End Subroutine row_echelon_rational

    ! Matrix example (Gaussian elimination). This brings the given mxn matrix A
    ! into (column) echelon form. The algorithm uses partial pivoting. The
    ! permutation of the columns is returned in the index array.

    ! Note that this algorithm is prepared to work on transposed matrices, as Pure
    ! matrices are stored in row-major order. Hence it computes a column echelon
    ! form. In Pure land this becomes a row echelon form which is what we want.

    Subroutine gauss(A, indexp)
    	real(kind=dp),dimension(:,:), intent(in out) :: A
    	integer,      dimension(:),   intent(out)    :: indexp
      integer       :: i, j, k, p, q, n, m
      real(kind=dp) :: pivot, x, y

      n=size(A,dim=2)
      indexp(1:n) = [(i,i=1,n)]
      do i = 1, n
         ! partial pivoting
         k = 0; pivot = 0.0_dp
         do j = i, n
            x = A(i, indexp(j))
            if (abs(x) > abs(pivot)) then
               k = j; pivot = x
            end if
         end do
         x = pivot
         if (abs(x) <= epsd) exit ! zero pivot, bail out
         ! the pivot column
         p = indexp(k)
         if (i /= k) then
            indexp(k) = indexp(i); indexp(i) = p
         end if
         ! normalize the pivot column
         A(:, p) = A(:, p) / x
         ! subtract multiples of the pivot column from the remaining columns
         do k = i+1, n
            q = indexp(k); y = A(i, q)
            A(:, q) = A(:, q) - y*A(:, p)
         end do
      end do
    End Subroutine gauss

  End Module Matrix_Mod