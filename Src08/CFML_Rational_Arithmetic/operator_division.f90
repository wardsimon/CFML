Submodule (CFML_Rational_Arithmetic) Operator_Division

 Contains
    !!----
    !!---- FUNCTION RATIONAL_DIVIDE 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Divide(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       type(rational)              :: res
       
       !---- Local Variables ----!
       integer(kind=il) :: denom
      
       denom = r%denominator * s%numerator
       if (denom /= 0_il) then
          res = r%numerator * s%denominator // denom
          res=rational_simplify(res)
       else
          res=0_il//1_il
       end if
       
       return
    End Function Rational_Divide
    
    !!----
    !!---- FUNCTION RATIONAL_DIVIDE_INT 
    !!----
    !!----
    !!
    Module Elemental Function Rational_Divide_Int(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       type(rational)               :: res
      
       if (i /= 0_il) then
          res = r%numerator // (r%numerator*i)
          res=rational_simplify(res)
       else
          res=0_il//1_il
       end if
       
       return
    End Function Rational_Divide_Int
 
End Submodule Operator_Division 