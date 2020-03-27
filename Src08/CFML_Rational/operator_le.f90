!!----
!!----
!!----
!!
Submodule (CFML_Rational) Operator_LE

 Contains
    !!----
    !!---- RATIONAL_LE
    !!----
    !!---- 08/04/2019
    !!
    Module Elemental Function Rational_LE(R, S) Result(Res)
       !---- Argument ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res

       !---- Local Variables ----!
       type(rational) :: r_simple, s_simple

       r_simple = rational_simplify(r)
       s_simple = rational_simplify(s)

       res = r_simple%numerator * s_simple%denominator <= &
           & s_simple%numerator * r_simple%denominator

       return
    End Function Rational_LE

    !!----
    !!---- RATIONAL_INTEGER_LE
    !!----
    !!---- 08/04/2019
    !!
    Module Elemental Function Rational_Integer_LE(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=LI),intent (in) :: i
       logical                      :: res

       !---- Local Variables ----!
       type(rational) :: r_simple

       r_simple = rational_simplify(r)
       res = r_simple%numerator <= i * r_simple%denominator

       return
    End Function Rational_Integer_LE

    !!----
    !!---- INTEGER_LE_RATIONAL
    !!----
    !!---- 08/04/2019
    !!
    Module Pure Function Integer_Rational_LE(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=LI),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res

       !---- Local Variables ----!
       type(rational) :: r_simple

       r_simple = rational_simplify(r)
       res = i * r_simple%denominator  <= r_simple%numerator

       return
    End Function Integer_Rational_LE

End Submodule Operator_LE

