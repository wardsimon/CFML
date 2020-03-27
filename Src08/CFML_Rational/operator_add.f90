!!----
!!----
!!----
!!
Submodule (CFML_Rational) Operator_Add

 Contains
    !!----
    !!---- RATIONAL_ADD
    !!----
    !!---- 08/04/2019
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
    !!---- RATIONAL_INTEGER_ADD
    !!----
    !!---- 08/04/2019
    !!
    Module Elemental Function Rational_Integer_Add(S, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: s
       integer(kind=LI),intent (in) :: i
       type(rational)               :: res

       res = (i * s%denominator +  s%numerator) // s%denominator
       res=rational_simplify(res)

       return
    End Function Rational_Integer_Add

    !!----
    !!---- INTEGER_RATIONAL_ADD
    !!----
    !!---- 08/04/2019
    !!
    Module Elemental Function Integer_Rational_Add(I, S) Result(Res)
       !---- Arguments ----!
       integer(kind=LI),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res

       res = (i * s%denominator +  s%numerator) // s%denominator
       res=rational_simplify(res)

       return
    End Function Integer_Rational_Add

End Submodule Operator_Add