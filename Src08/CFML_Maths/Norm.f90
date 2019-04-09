!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Maths) Norm
 Contains
 
    !!---- 
    !!---- NORM_I
    !!----    Calculate the Norm of a vector
    !!----
    !!---- 28/03/2019 
    !!
    Module Pure Function Norm_I(X,G) Result(R)    
       !---- Arguments ----!
       integer,       dimension(:),   intent(in) :: x    ! Input vector
       real(kind=cp), dimension(:,:), intent(in) :: g    ! Metric array
       real(kind=cp)                             :: r    ! Norm of the input vector

       if (size(x)*size(x) /= size(g)) then
          r=tiny(0.0_cp)
       else
          r=sqrt(dot_product(real(x), matmul(g,real(x))))
       end if

       return
    End Function Norm_I
 
    !!---- 
    !!---- NORM_R
    !!----    Calculate the Norm of a vector
    !!----
    !!---- 28/03/2019 
    !!
    Module Pure Function Norm_R(X,G) Result(R)    
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: x   ! Input vector
       real(kind=cp), dimension(:,:), intent(in) :: g   ! Metrics
       real(kind=cp)                             :: r   ! Norm of the vector

       if (size(x)*size(x) /= size(g)) then
          r=tiny(0.0_cp)
       else
          r=sqrt(dot_product(x, matmul(g,x)))
       end if

       return
    End Function Norm_R
   
End Submodule Norm
