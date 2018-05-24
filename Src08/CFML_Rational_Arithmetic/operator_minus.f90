Submodule (CFML_Rational_Arithmetic) Operator_Minus

 Contains
    !!----
    !!---- FUNCTION RATIONAL_MINUS 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Minus(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational)              :: res
      
       res = - r%numerator // r%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Minus
    
    !!----
    !!---- FUNCTION RATIONAL_SUBTRACT 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Subtract(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       type(rational)              :: res
      
       res = r%numerator * s%denominator - r%denominator * s%numerator // &
           & r%denominator * s%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Subtract
    
    !!----
    !!---- FUNCTION INTEGER_RATIONAL_SUBTRACT 
    !!----
    !!----
    !!
    Module Elemental Function Integer_Rational_Subtract(I, S) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res
      
       res = (i * s%denominator -  s%numerator) //  s%denominator
       res=rational_simplify(res)
       
       return
    End Function Integer_Rational_Subtract
    
    !!----
    !!---- FUNCTION RATIONAL_INTEGER_SUBTRACT 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Integer_Subtract(S,I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: s
       integer(kind=il),intent (in) :: i
       type(rational)               :: res
      
       res = (s%numerator - i * s%denominator) //  s%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Integer_Subtract
    
    
 
End Submodule Operator_Minus 