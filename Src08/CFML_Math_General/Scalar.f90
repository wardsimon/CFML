!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Scalar
 Contains
 
    !!
    Module Function Scalar_I(X,Y,G) Result(R)    
       !---- Arguments ----!
       integer, dimension(:),         intent(in) :: x,y     ! Input vectors
       real(kind=cp), dimension(:,:), intent(in) :: g       ! Metrics
       real(kind=cp)                             :: r       ! Scalar

       if (size(x)/= size(y) .or. size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=dot_product(real(x), matmul(g,real(y)))
       end if

       return
    End Function Scalar_I
 
    !!--++ FUNCTION SCALAR_R
    !!--++
    !!--++    Scalar Product including metrics
    !!--++
    !!--++ Update: April - 2009
    !!
    Module Function Scalar_R(X,Y,G) Result(R)    
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: x,y    ! Input vectors
       real(kind=cp), dimension(:,:), intent(in) :: g      ! Metrics
       real(kind=cp)                             :: r      ! Scalar

       if (size(x)/= size(y) .or. size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=dot_product(x, matmul(g,y))
       end if

       return
    End Function Scalar_R
 
   
End Submodule Scalar
