Submodule (CFML_Rational_Arithmetic) Constructor

 Contains
    !!----
    !!---- FUNCTION MAKE_RATIONAL
    !!----
    !!----
    !!
    Module Elemental Function Make_Rational(Numerator, Denominator) Result(Res)
       !---- Arguments ----!
       integer(kind=il), intent (in) :: numerator
       integer(kind=il), intent (in) :: denominator
       type(rational)                :: res
       
       res = rational(numerator, denominator)
       
       return
    End Function Make_Rational

    !!----
    !!---- FUNCTION MAKE_RATIONAL_INT
    !!----
    !!----
    !!
    Module Elemental Function Make_Rational_Int(Numerator, Denominator) Result(Res)
       !---- Arguments ----! 
       integer, intent (in) :: numerator
       integer, intent (in) :: denominator
       type(rational)      :: res
      
       res = rational(int(numerator,kind=il), int(denominator,kind=il))
       
       return
    End Function Make_Rational_Int
 
End Submodule Constructor 