!!----
!!----
!!----
!!
Submodule (CFML_Rational) Operator_EQ

 Contains
    !!----
    !!---- RATIONAL_EQ
    !!----
    !!---- 08/04/2019
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
    !!---- RATIONAL_INTEGER_EQ
    !!----
    !!---- 08/04/2019
    !!
    Module Elemental Function Rational_Integer_EQ(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=LI),intent (in) :: i
       logical                      :: res

       res = r%denominator == 1_LI .and. r%numerator == i

       return
    End Function Rational_Integer_EQ

    !!----
    !!---- INTEGER_RATIONAL_EQ
    !!----
    !!---- 08/04/2019
    !!
    Module Elemental Function Integer_Rational_EQ(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=LI),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res

       res = r%denominator == 1_LI .and. r%numerator == i

       return
    End Function Integer_Rational_EQ

End Submodule Operator_EQ
