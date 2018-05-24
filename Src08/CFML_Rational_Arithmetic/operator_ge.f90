Submodule (CFML_Rational_Arithmetic) Operator_GE

 Contains
    !!----
    !!---- FUNCTION RATIONAL_GE 
    !!----
    !!----
    !!
    Module Pure Function Rational_GE(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple, s_simple
 
       r_simple = rational_simplify(r)
       s_simple = rational_simplify(s)
      
       res = r_simple%numerator * s_simple%denominator >= &
           & s_simple%numerator * r_simple%denominator
           
       return    
    End Function Rational_GE

    !!----
    !!---- FUNCTION RATIONAL_GE_INTEGER 
    !!----
    !!----
    !!
    Module Pure Function Rational_GE_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
      
       r_simple = rational_simplify(r)
       res = r_simple%numerator >= i * r_simple%denominator
       
       return
    End Function Rational_GE_Integer

    !!----
    !!---- FUNCTION INTEGER_GE_RATIONAL
    !!----
    !!----
    !!
    Module Pure Function Integer_GE_Rational(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
       
       r_simple = rational_simplify(r)
       res = i * r_simple%denominator >= r_simple%numerator
       
       return
    End Function Integer_GE_Rational
 
End Submodule Operator_GE