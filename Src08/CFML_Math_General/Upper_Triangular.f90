!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Upper_Triangular
 Contains
 
    !!---- FUNCTION UPPER_TRIANGULAR_I
    !!----
    !!----   Updated: October - 2014
    !!----
    Module Function Upper_Triangular_I(A,n) Result (T)    
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: A     ! Input array
       integer,                 intent(in) :: n     ! Dimension
       integer, dimension(n,n)             :: T

       !---- Local Variables ----!
       integer :: i,j,p,q,m

       m=n
       p=size(A(:,1)); q=size(A(1,:))
       if(n > p .or. n > q) m=min(p,q)

       T=0
       do j=1,m
          do i=1,j
             T(i,j)=A(i,j)
          end do
       end do

       return
    End Function Upper_Triangular_I
 
    !!---- FUNCTION UPPER_TRIANGULAR_R
    !!----
    !!----   Updated: October - 2014
    !!----
    Module Function Upper_Triangular_R(A,n) Result (T)    
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: A   ! Input array
       integer,                       intent(in) :: n   ! Dimension
       real(kind=cp), dimension(n,n)             :: T

       !---- Local Variables ----!
       integer :: i,j,p,q,m

       m=n
       p=size(A(:,1)); q=size(A(1,:))
       if(n > p .or. n > q) m=min(p,q)

       T=0.0_cp
       do j=1,m
          do i=1,j
             T(i,j)=A(i,j)
          end do
       end do

       return
    End Function  Upper_Triangular_R
 
   
End Submodule Upper_Triangular
