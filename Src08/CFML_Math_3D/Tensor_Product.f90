!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Tensor_Product
 Contains
 
    !!---- FUNCTION TENSOR_PRODUCT_CMPL_DP
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Module Function Tensor_Product_cmpl_dp(Vec1,Vec2) Result(w)    
       !---- Argument ----!
       complex(kind=dp), dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2
       complex(kind=dp), dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2

       !---- Local Arguments ----!
       complex(kind=dp), dimension(3,3)            :: mu,mv

       mu=0.0_dp;  mv=0.0_dp
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)

       return
    End Function Tensor_Product_cmpl_dp
 
    !!---- FUNCTION TENSOR_PRODUCT_CMPL_SP
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Module Function Tensor_Product_cmpl_sp(Vec1,Vec2) Result(w)    
       !---- Argument ----!
       complex(kind=sp), dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2
       complex(kind=sp), dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2

       !---- Local Arguments ----!
       complex(kind=sp), dimension(3,3)            :: mu,mv

       mu=0.0_sp;  mv=0.0_sp
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)
       return
    End Function Tensor_Product_cmpl_sp
 
    !!---- FUNCTION TENSOR_PRODUCT_DP
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Module Function Tensor_Product_dp(Vec1,Vec2) Result(w)    
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: Vec1,Vec2
       real(kind=dp), dimension(3,3)            :: w

       !---- Local Variables ----!
       real(kind=dp), dimension(3,3)            :: mu,mv

       mu=0.0_dp;  mv=0.0_dp
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)

       return
    End Function Tensor_Product_dp
 
    !!---- FUNCTION TENSOR_PRODUCT_IN
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Module Function Tensor_Product_in(Vec1,Vec2) Result(w)    
       !---- Argument ----!
       integer, dimension(3), intent( in) :: Vec1,Vec2
       integer, dimension(3,3)            :: w

       !---- Local Variables ----!
       integer, dimension(3,3)            :: mu,mv
       mu=0;  mv=0
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)

       return
    End Function Tensor_Product_in
 
    !!---- FUNCTION TENSOR_PRODUCT_SP
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Module Function Tensor_Product_sp(Vec1,Vec2) Result(w)    
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: Vec1,Vec2
       real(kind=sp), dimension(3,3)            :: w

       !---- Local Variables ----!
       real(kind=sp), dimension(3,3)            :: mu,mv

       mu=0.0_sp;  mv=0.0_sp
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)

       return
    End Function Tensor_Product_sp
 
   
End Submodule Tensor_Product
