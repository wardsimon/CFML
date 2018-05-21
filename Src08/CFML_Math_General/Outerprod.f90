!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Outerprod
 Contains
 
    !!--++ FUNCTION OUTERPROD_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = a(i)*b(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Outerprod_dp(a,b)  Result(c)    
       !---- Arguments ----!
       real(kind=dp),dimension(:),intent(in)    :: a,b
       real(kind=dp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod_dp
 
    !!--++ FUNCTION OUTERPROD_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = a(i)*b(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Outerprod_sp(a,b)  Result(c)    
       !---- Arguments ----!
       real(kind=sp),dimension(:),intent(in)    :: a,b
       real(kind=sp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod_sp
 
   
End Submodule Outerprod
