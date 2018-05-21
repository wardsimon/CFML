!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Invert_Array3x3
 Contains
 
    !!--++  FUNCTION INVERT_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate de inverse of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Invert_Dp(array) Result(b)    
       !---- Arguments ----!
       real(kind=dp),dimension(3,3), intent(in) :: array
       real(kind=dp),dimension(3,3)             :: b

       !---- Local variables ----!
       real(kind=dp)  :: dmat

       b(1,1) =   array(2,2)*array(3,3)-array(2,3)*array(3,2)
       b(2,1) = -(array(2,1)*array(3,3)-array(2,3)*array(3,1))
       b(3,1) =   array(2,1)*array(3,2)-array(2,2)*array(3,1)
       b(1,2) = -(array(1,2)*array(3,3)-array(1,3)*array(3,2))
       b(2,2) =   array(1,1)*array(3,3)-array(1,3)*array(3,1)
       b(3,2) = -(array(1,1)*array(3,2)-array(1,2)*array(3,1))
       b(1,3) =   array(1,2)*array(2,3)-array(1,3)*array(2,2)
       b(2,3) = -(array(1,1)*array(2,3)-array(1,3)*array(2,1))
       b(3,3) =   array(1,1)*array(2,2)-array(1,2)*array(2,1)

       !> Determinant
       dmat = array(1,1)*b(1,1)+array(1,2)*b(2,1)+array(1,3)*b(3,1)

       if (abs(dmat) > tiny(dmat)) then
          b= b/dmat
       else
          b=0.0_dp
       end if

       return
    End Function Invert_Dp
 
End Submodule Invert_Array3x3
