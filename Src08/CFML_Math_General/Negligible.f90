!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Negligible
 Contains
 
    !!--++ FUNCTION NEGLIGIBLEC
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if a complex number is negligible (abs < EPSS)
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Elemental Function Negligiblec(v) Result(Neglig)    
       !---- Argument ----!
       complex, intent( in) :: v         ! Complex number
       logical              :: Neglig

       !> Init
       Neglig=.false.

       if (abs(v) > epss) return
       Neglig=.true.

       return
    End Function Negligiblec
 
    !!--++ FUNCTION NEGLIGIBLER
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is negligible (abs < EPSS)
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Elemental Function Negligibler(v) Result(neglig)    
       !---- Argument ----!
       real(kind=cp), intent( in) :: v          ! Real number
       logical                    :: Neglig

       !> Init
       Neglig=.false.

       if (abs(v) > epss) return
       Neglig=.true.

       return
    End Function Negligibler
   
End Submodule Negligible
