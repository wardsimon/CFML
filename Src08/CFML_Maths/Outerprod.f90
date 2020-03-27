!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Maths) CFML_Math_005
 Contains

    !!----
    !!---- OUTERPROD_R
    !!----    Computes the outer product (tensorial product) of two
    !!----    vectors to give a tensor (matrix) as the result:
    !!----                   c(i,j) = a(i)*b(j).
    !!----
    !!----    It uses the intrinsic Fortran 90 function SPREAD.
    !!----    Taken from Numerical Recipes.
    !!----
    !!---- 28/03/2019
    !!
    Module Pure Function Outerprod(a,b)  Result(c)
       !---- Arguments ----!
       real(kind=cp),dimension(:),intent(in)    :: a,b
       real(kind=cp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod

End Submodule CFML_Math_005
