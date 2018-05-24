Submodule (CFML_Rational_Arithmetic) Operator_Multiply

 Contains
    !!----
    !!---- FUNCTION RATIONAL_MULTIPLY 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Multiply(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       type(rational)              :: res
      
       res = r%numerator * s%numerator // (r%denominator * s%denominator)
       res=rational_simplify(res)
       
       return
    End Function Rational_Multiply

    !!----
    !!---- FUNCTION RATIONAL_MULTIPLY 
    !!----
    !!----
    !!
    Module Elemental Function Integer_Rational_Multiply(I, S) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res
      
       res = i * s%numerator // s%denominator
       res=rational_simplify(res)
       
       return
    End Function Integer_Rational_Multiply

    !!----
    !!---- FUNCTION RATIONAL_INTEGER_MULTIPLY 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Integer_Multiply(S,I) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res
      
       res = i * s%numerator // s%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Integer_Multiply
 
End Submodule Operator_Multiply 