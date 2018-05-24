Submodule (CFML_Rational_Arithmetic) Operator_GT

 Contains
    !!----
    !!---- FUNCTION RATIONAL_GT 
    !!----
    !!----
    !!
    Module Pure Function Rational_GT(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple, s_simple

       r_simple = rational_simplify(r)
       s_simple = rational_simplify(s)
      
       res = r_simple%numerator * s_simple%denominator > &
           & s_simple%numerator * r_simple%denominator
           
       return    
    End Function Rational_GT

    !!----
    !!---- FUNCTION RATIONAL_GT_INTEGER 
    !!----
    !!----
    !!
    Module Pure Function Rational_Gt_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),   intent (in) :: r
       integer(kind=il), intent (in) :: i
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
      
       r_simple = rational_simplify(r)
      
       res = r_simple%numerator > i * r_simple%denominator
       
       return
    End Function Rational_Gt_Integer

    !!----
    !!---- FUNCTION INTEGER_GT_RATIONAL 
    !!----
    !!----
    !!
    Module Pure Function Integer_Gt_Rational(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
      
       r_simple = rational_simplify(r)
       res = i * r_simple%denominator > r_simple%numerator
       
       return
    End Function Integer_Gt_Rational
    
End Submodule Operator_GT