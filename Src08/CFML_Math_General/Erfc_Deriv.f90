!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Erfc_Deriv
 Contains
 
    !!--++
    !!--++ Function Erfc_Deriv_DP(X)
    !!--++
    !!--++    (Overloaded)
    !!--++    Derivative of the complementary error function
    !!--++
    !!--++ Update: October - 2005
    !!
    Module Function Erfc_Deriv_DP(X) Result(Der)
       !---- Argument ----!
       real(kind=dp), intent(in)    :: x
       real(kind=dp)                :: der

       der = - 1.128379167096_dp * exp(-x*x)   !from Abramovitz

       return
    End Function Erfc_Deriv_DP

    !!--++
    !!--++ Function Erf_Deriv_SP(X)
    !!--++
    !!--++    (Overloaded)
    !!--++    Derivative of the complementary error function
    !!--++
    !!--++ Update: October - 2005
    !!
    Module Function Erfc_Deriv_SP(X)  Result(Der)
       !---- Argument ----!
       real(kind=sp), intent(in)    :: x
       real(kind=sp)                :: der

       der = - 1.128379167096 * exp(-x*x)

       return
    End Function Erfc_Deriv_SP
 
   
End Submodule Erfc_Deriv
