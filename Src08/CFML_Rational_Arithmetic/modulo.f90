Submodule (CFML_Rational_Arithmetic) Ov_Modulo

 Contains
    !!----
    !!---- FUNCTION RATIONAL_MODULO
    !!----
    !!----
    !!
    Module Elemental Function Rational_Modulo(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       integer                     :: res
      
       res = modulo (r%numerator, r%denominator)
       
       return
    End Function Rational_Modulo

    !!----
    !!---- FUNCTION RATIONAL_MODULO_INT
    !!----
    !!----
    !!
    Module Elemental Function Rational_Modulo_Int(R,I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent(in)  :: i
       type(rational)               :: res
      
       !---- Local Variables ----!
       real(kind=dp) :: val
      
       val = modulo (real(r%numerator,kind=dp) / real(r%denominator,kind=dp),real(i,kind=dp))
       res = val
       
       return
    End Function Rational_Modulo_Int
 
End Submodule Ov_Modulo 