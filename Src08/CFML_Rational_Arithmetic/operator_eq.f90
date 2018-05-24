Submodule (CFML_Rational_Arithmetic) Operator_EQ

 Contains
    !!----
    !!---- FUNCTION RATIONAL_EQ
    !!----
    !!----
    !!
    Module Elemental Function Rational_EQ(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
      
       res = r%numerator * s%denominator == s%numerator * r%denominator
       
       return
    End Function Rational_EQ

    !!----
    !!---- FUNCTION RATIONAL_EQ_INTEGER
    !!----
    !!----
    !!
    Module Elemental Function Rational_EQ_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       logical                      :: res
      
       res = r%denominator == 1_il .and. r%numerator == i
       
       return
    End Function Rational_EQ_Integer
    
    !!----
    !!---- FUNCTION INTEGER_EQ_RATIONAL
    !!----
    !!----
    !!
    Module Elemental Function Integer_Eq_Rational(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
       
       res = r%denominator == 1_il .and. r%numerator == i
       
       return
    End Function Integer_EQ_Rational
    
End Submodule Operator_EQ