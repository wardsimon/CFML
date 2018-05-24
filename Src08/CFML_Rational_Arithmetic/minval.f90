Submodule (CFML_Rational_Arithmetic) Minval

 Contains
    !!----
    !!---- FUNCTION RATIONAL_MINVAL_VECTOR
    !!----
    !!----
    !!
    Module Function Rational_Minval_Vector(R) Result(Res)
       !---- Arguments ----!
       type(rational), dimension(:), intent (in) :: r
       type(rational)                            :: res
       
       !---- Local Variables ----!
       integer :: i,n
       
       n=size(r)
       res=huge(1_il)//1_il
       do i=1,n
          if (r(i) < res) res=r(i)
       end do
       
       return
    End Function Rational_Minval_Vector

    !!----
    !!---- FUNCTION RATIONAL_MINVAL_MATRIX
    !!----
    !!----
    !!
    Module Function Rational_Minval_Matrix(R) Result(Res)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent (in) :: r
       type(rational)                              :: res
       
       !---- Local Variables ----!
       integer :: i,j,n1,n2
      
       n1=size(r,dim=1);  n2=size(r,dim=2)
       res=huge(1_il)//1_il
       do j=1,n2
          do i=1,n1
             if (r(i,j) < res) res=r(i,j)
          end do
       end do
       
       return
    End Function Rational_Minval_Matrix
 
End Submodule Minval