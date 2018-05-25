Submodule (CFML_Rational_Arithmetic) Overloads

 Contains
    !!----
    !!---- FUNCTION RATIONAL_ABS
    !!----
    !!----
    !!
    Module Elemental Function Rational_Abs(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational)              :: res
       
       res = sign (r%numerator, r%denominator) // r%denominator
       
       return
    End Function Rational_Abs
    
    !!----
    !!---- FUNCTION RATIONAL_INT
    !!----
    !!----
    !!
    Module Elemental Function Rational_Int(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       integer(kind=il)            :: res
       
       res = r%numerator / r%denominator
       
       return
    End Function Rational_Int

    !!----
    !!---- FUNCTION NINT_RATIONAL
    !!----
    !!----
    !!
    Module Elemental Function Nint_Rational(R) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in)  :: r
       integer(kind=il)              :: res
      
       res = nint(real(r%numerator,kind=dp)/real(r%denominator,kind=dp),kind=il)
       
       return
    End Function Nint_Rational
    
    !!----
    !!---- FUNCTION RATIONAL_DOT_PRODUCT
    !!----
    !!----
    !!
    Module Function Rational_Dot_Product(R1,R2) Result(Res)
       !---- Arguments ----!
       type(rational), dimension(:), intent (in) :: r1
       type(rational), dimension(:), intent (in) :: r2
       type(rational)                            :: res
      
       !---- Local Variables ----!
       integer :: n1,n2,i
      
       n1=size(r1); n2=size(r2)
       res=0_il//1_il
      
       if (n1 == n2) then
          do i=1,n1
            res = rational_simplify(res + r1(i)*r2(i))
          end do
       else
         Err_CFML%state=.true.
         err_cfml%flag=2
         write(unit=err_cfml%msg,fmt="(a,2i4)") &
              "Error in DOT_PRODUCT: the dimensions of the arguments are different",n1,n2
       end if
       
       return
    End Function Rational_Dot_Product
    
    !!----
    !!---- FUNCTION RATIONAL_SUM_VEC
    !!----
    !!----
    !!
    Module Pure Function Rational_Sum_Vec(Vec) Result(Suma)
       !---- Arguments ----!
       type(rational), dimension(:), intent(in) :: vec
       type(rational)                           :: suma
      
       !---- Local variables ----!
       real(kind=dp), dimension(size(vec)) :: rvec
      
       rvec=vec
       suma=sum(rvec)
       suma=rational_simplify(suma)
       
       return
    End Function Rational_Sum_Vec
    
 
End Submodule Overloads