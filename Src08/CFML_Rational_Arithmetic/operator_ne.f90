Submodule (CFML_Rational_Arithmetic) Operator_NE

 Contains
    !!----
    !!---- FUNCTION RATIONAL_NE
    !!----
    !!----
    !!
    Module Pure Function Rational_Ne(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
      
       res = r%numerator * s%denominator /= s%numerator * r%denominator
       
       return
    End Function Rational_Ne

    !!----
    !!---- FUNCTION RATIONAL_NE_INTEGER
    !!----
    !!----
    !!
    Module Elemental Function Rational_Ne_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       logical                      :: res
      
       res = r%numerator /= i * r%denominator
       
       return
    End Function Rational_Ne_Integer

    !!----
    !!---- FUNCTION INTEGER_NE_RATIONAL
    !!----
    !!----
    !!
    Module Elemental Function Integer_Ne_Rational(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
      
       res = r%numerator /= i * r%denominator
       
       return
    End Function Integer_Ne_Rational
    
End Submodule Operator_NE