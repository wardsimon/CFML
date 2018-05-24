Submodule (CFML_Rational_Arithmetic) Maxloc

 Contains
    !!----
    !!---- FUNCTION RATIONAL_MAXLOC_MAT
    !!----
    !!----
    !!
    Module Pure Function Rational_Maxloc_Mat(Mat) Result(Pos_Max)
       !---- Arguments ----!
       type(rational),  dimension(:,:), intent(in) :: Mat
       integer(kind=il),dimension(2)               :: pos_max
      
       !---- Local variables ----!
       integer        :: nu1,nl1,nu2,nl2,i,j
       type(rational) :: res
      
       nu1=ubound(Mat,dim=1); nu2=ubound(Mat,dim=2)
       nl1=lbound(Mat,dim=1); nl2=lbound(Mat,dim=2)

       res=-huge(1_il)//1_il
       do j=nl2,nu2
          do i=nl1,nu1
             if (mat(i,j) > res) then
                res=mat(i,j)
                pos_max=[i,j]
             end if
          end do
       end do
       
       return
    End Function Rational_Maxloc_Mat

    !!----
    !!---- FUNCTION RATIONAL_MAXLOC_VECT
    !!----
    !!----
    !!
    Module Pure Function Rational_Maxloc_Vect(Vec) Result(Pos_Max)
       !---- Arguments ----!
       type(rational), dimension(:), intent(in) :: vec
       integer                                  :: pos_max
      
      
       !---- Local variables ----!
       integer:: nu,nl,i
       type(rational) :: res
      
       nu=ubound(vec,dim=1)
       nl=lbound(vec,dim=1)
       res=-huge(1_il)//1_il
       do i=nl,nu
          if (vec(i) > res) then
             res=vec(i)
             pos_max=i
          end if
       end do
       
       return
    End Function Rational_Maxloc_Vect
    
End Submodule Maxloc