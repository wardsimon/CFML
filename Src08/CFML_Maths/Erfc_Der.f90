!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Maths) CFML_Maths_012
 Contains
 
    !!----
    !!---- Erfc_Deriv(X)
    !!----    Derivative of the complementary error function
    !!----
    !!---- 09/04/2019 
    !!
    Module Elemental Function Erfc_Deriv(X) Result(Der)
       !---- Argument ----!
       real(kind=cp), intent(in)    :: x
       real(kind=cp)                :: der

       der = - 1.128379167096_cp * exp(-x*x)   !from Abramovitz

       return
    End Function Erfc_Deriv

End Submodule CFML_Maths_012
