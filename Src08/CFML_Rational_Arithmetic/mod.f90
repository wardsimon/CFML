Submodule (CFML_Rational_Arithmetic) Ov_Mod

 Contains
    !!----
    !!---- FUNCTION RATIONAL_MODULO
    !!----
    !!----
    !!
    Module Elemental Function Rational_Mod(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       integer                     :: res
      
       res = mod(r%numerator, r%denominator)
       
       return
    End Function Rational_Mod

    !!----
    !!---- FUNCTION RATIONAL_MODULO_INT
    !!----
    !!----
    !!
    Module Elemental Function Rational_Mod_Int(R,I) Result(Res)
       !---- Arguments ----!
      type(rational),   intent (in) :: r
      integer(kind=il), intent (in) :: i
      type(rational)                :: res
      
      !---- Local Variables ----!
      real(kind=dp) :: val
      
      val = mod(real(r%numerator,kind=dp) / real(r%denominator,kind=dp),real(i,kind=dp))
      res = val
      
      return
    End Function Rational_Mod_Int
 
End Submodule Ov_Mod