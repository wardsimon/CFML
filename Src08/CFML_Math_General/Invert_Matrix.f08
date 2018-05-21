!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Invert_Matrix
 Contains
 
    !!---- SUBROUTINE INVERT_MATRIX
    !!--<<
    !!----    Subroutine to invert a real matrix using LU decomposition.
    !!----    In case of singular matrix (singular=.true.) instead of the inverse
    !!----    matrix, the subroutine provides the LU decomposed matrix as used
    !!----    in Numerical Recipes.
    !!----    The input matrix is preserved and its inverse (or its LU decomposition)
    !!----    is provided in "b". The optional argument "perm" holds the row permutation
    !!----    performed to obtain the LU decomposition.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine Invert_Matrix(a,b,singular,perm)    
       !---- Arguments ----!
       real(kind=cp), dimension(:,:),  intent(in ) :: a         ! Input Array
       real(kind=cp), dimension(:,:),  intent(out) :: b         ! Output array
       logical,                        intent(out) :: singular
       integer, dimension(:),optional, intent(out) :: perm

       !---- Local variables ----!
       integer                                       :: i,n
       integer,       dimension(size(a,1))           :: indx
       real(kind=cp)                                 :: d, det
       real(kind=cp), dimension(size(a,1),size(a,1)) :: lu

       n=size(a,1)
       lu=a(1:n,1:n)

       call LU_Decomp(lu,d,singular,indx)
       if (present(perm)) perm(1:n)=indx(1:n)

       if (singular) then
          b=lu
          return
       else
          det=0.0
          do i=1,n
             d=d*sign(1.0_cp,lu(i,i))
             det=det + log(abs(lu(i,i)))
          end do
          det=d*exp(det)
          if (abs(det) <= 1.0e-36) then
             singular=.true.
             b=lu
             return
          end if
       end if

       b=0.0
       do i=1,n
          b(i,i)=1.0
          call LU_backsub(lu,indx,b(:,i))
       end do

       return
    End Subroutine Invert_Matrix
    
    !!---- SUBROUTINE MATINV
    !!----
    !!----  Subroutine for inverting a real square matrix.
    !!----  The input matrix is replaced in output with its inverse.
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine Matinv(a,n)    
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       integer     ,                  intent(in)     :: n

       !---- Local variables ----!
       real(kind=cp)                 :: amax,savec
       integer, dimension(size(a,1)) :: ik,jk
       integer                       :: i,j,k,l

       !---- Subroutine to invert a real matrix ----!
       do k=1,n
          amax=0.0
          do
             do
                do i=k,n
                   do j=k,n
                      if (abs(amax)-abs(a(i,j)) > 0.0) cycle
                      amax=a(i,j)
                      ik(k)=i
                      jk(k)=j
                   end do
                end do
                i=ik(k)
                if (i-k < 0) cycle
                exit
             end do

             if (i-k /= 0) then
                do j=1,n
                   savec=a(k,j)
                   a(k,j)=a(i,j)
                   a(i,j)=-savec
                end do
             end if

             j=jk(k)
             if (j-k < 0) cycle
             exit
          end do

          if (j-k /= 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=a(i,j)
                a(i,j)=-savec
             end do
          end if

          do i=1,n
             if (i-k /= 0)  then
                a(i,k)=-a(i,k)/amax
             end if
          end do
          do i=1,n
             do j=1,n
                if (i-k == 0 .or. j-k == 0) cycle
                a(i,j)=a(i,j)+a(i,k)*a(k,j)
             end do
          end do
          do j=1,n
             if (j-k == 0)   cycle
             a(k,j)=a(k,j)/amax
          end do
          a(k,k)=1.0/amax
       end do     !k

       do l=1,n
          k=n-l+1
          j=ik(k)
          if (j-k > 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=-a(i,j)
                a(i,j)=savec
             end do
          end if
          i=jk(k)
          if (i-k > 0) then
             do j=1,n
                savec=a(k,j)
                a(k,j)=-a(i,j)
                a(i,j)=savec
             end do
          end if
       end do

       return
    End Subroutine Matinv
    
    !!---- SUBROUTINE LU_BACKSUB
    !!--<<
    !!----    Adapted from Numerical Recipes.
    !!----    Solves the set of N linear equations A  X = B. Here the N x N matrix A is input,
    !!----    not as the original matrix A, but rather as its LU decomposition, determined
    !!----    by the routine LU_DECOMP. INDX is input as the permutation vector of length N
    !!----    returned by LU_DECOMP. B is input as the right-hand-side vector B,
    !!----    also of length N, and returns with the solution vector X.
    !!----    A and INDX are not modified by this routine and can be left in place for successive calls
    !!----    with different right-hand sides B. This routine takes into account the possibility that B will
    !!----    begin with many zero elements, so it is efficient for use in matrix inversion.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine LU_Backsub(a,indx,b)    
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in)     :: a
       integer,         dimension(:), intent(in)     :: indx
       real(kind=cp),   dimension(:), intent(in out) :: b

       !---- Local Variables ----!
       integer       :: i,ii,ll,n
       real(kind=cp) :: summ

       n=size(a,1)
       ii=0              !When ii is set to a positive value, it will become the index
       do i=1,n          !of the first nonvanishing element of b. We now do
          ll=indx(i)     !the forward substitution. The only new wrinkle is to
          summ=b(ll)     !unscramble the permutation as we go.
          b(ll)=b(i)
          if (ii /= 0) then
             summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
          else if(summ /= 0.0) then   !A nonzero element was encountered, so from now on
             ii=i                       !we will have to do the dot product above.
          end if
          b(i)=summ
       end do

       do i=n,1,-1       !Now we do the backsubstitution
          b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
       end do

       return
    End Subroutine LU_Backsub
    
    !!---- SUBROUTINE LU_DECOMP
    !!--<<
    !!----    Subroutine to make the LU decomposition of an input matrix A.
    !!----    The input matrix is destroyed and replaced by a matrix containing
    !!----    in its upper triangular part (plus diagonal) the matrix U. The
    !!----    lower triangular part contains the nontrivial part (Lii=1) of matrix L.
    !!----    The output is rowwise permutation of the initial matrix. The vector INDX
    !!----    recording the row permutation. D is output as +/-1 depending on whether
    !!----    the number of row interchanges was even or odd, respectively.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine LU_Decomp(a,d,singular,indx)    
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       real(kind=cp),                 intent(out)    :: d
       logical,                       intent(out)    :: singular
       integer,  dimension(:), intent(out), optional :: indx

       !---- Local variables ----!
       real(kind=cp), dimension(size(a,1)):: vv  !vv stores the implicit scaling of each row.
       real(kind=cp), parameter           :: vtiny = 1.0e-7_sp !A small number.
       integer                            :: j,jj,imax,n

       singular=.false.
       n=size(a,1)
       if (present(indx)) then
          do j=1,n
             indx(j)=j
          end do
       end if
       d=1.0                      !No row interchanges yet.
       vv=maxval(abs(a),dim=2)    !Loop over rows to get the implicit scaling information.
       if (any(abs(vv) <= vtiny)) then   !There is a row of zeros.
          singular=.true.
          return
       end if
       vv=1.0_sp/vv     !Save the scaling.

       do j=1,n
          jj=maxloc(vv(j:n)*abs(a(j:n,j)),dim=1)
          imax=(j-1)+jj                               !Find the pivot row.
          if (j /= imax) then                         !Do we need to interchange rows?
             call swap(a(imax,:),a(j,:))              !Yes, do so...
             d=-d                                     !...and change the parity of d.
             vv(imax)=vv(j)                           !Also interchange the scale factor.
          end if
          if (present(indx)) indx(j)=imax
          if (abs(a(j,j)) <= vtiny) then !If the pivot element is zero the matrix is singular.
             a(j,j)=vtiny                !(at least to the precision of the algorithm)
             singular=.true.             !For some applications on singular matrices,
             !return                      !it is desirable to substitute vtiny for zero.
          end if                         !This is actually the present case
          a(j+1:n,j)=a(j+1:n,j)/a(j,j)                                    !Divide by the pivot element.
          a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))  !Reduce remaining submatrix.
       end do

       return
    End Subroutine LU_Decomp
    
    
 
End Submodule Invert_Matrix
