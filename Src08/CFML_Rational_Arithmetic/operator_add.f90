Submodule (CFML_Rational_Arithmetic) Operator_Add 

 Contains
    !!----
    !!---- FUNCTION RATIONAL_ADD 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Add(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       type(rational)              :: res
      
       res = r%numerator * s%denominator + r%denominator * s%numerator // &
           & r%denominator * s%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Add
    
    !!----
    !!---- FUNCTION RATIONAL_INTEGER_ADD 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Integer_Add(S, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: s
       integer(kind=il),intent (in) :: i
       type(rational)               :: res
      
       res = (i * s%denominator +  s%numerator) // s%denominator
       res=rational_simplify(res)
    
       return   
    End Function Rational_Integer_Add
    
    !!----
    !!---- FUNCTION INTEGER_RATIONAL_ADD 
    !!----
    !!----
    !!
    Module Elemental Function Integer_Rational_Add(I, S) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res
      
       res = (i * s%denominator +  s%numerator) // s%denominator
       res=rational_simplify(res)
       
       return
    End Function Integer_Rational_Add
 
End Submodule Operator_Add 