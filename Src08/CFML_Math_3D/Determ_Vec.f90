!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Determ_Vec
 Contains
 
    !!--++ FUNCTION DETERM_V_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of the components of three vectors
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Determ_V_I(Vec1,Vec2,Vec3) Result(det)    
       !---- Arguments ----!
       integer, dimension(3), intent(in) :: Vec1,Vec2,Vec3
       integer                           :: det

       !---- Local variables ----!
       integer :: i,j,k

       det = 0
       do i = 1,3
          j = i+1
          if (j == 4) j = 1
          k = 6-i-j
          det = det+Vec1(i)*(Vec2(j)*Vec3(k)-Vec2(k)*Vec3(j))
       end do

       return
    End Function Determ_V_I
    
    !!--++ FUNCTION DETERM_V_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of the components of three vectors
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Determ_V_R(Vec1,Vec2,Vec3) Result(det)    
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: Vec1,Vec2,Vec3
       real(kind=cp)                           :: det

       !---- Local variables ----!
       integer :: i,j,k

       det = 0.0_cp
       do i = 1,3
          j = i+1
          if (j == 4) j = 1
          k = 6-i-j
          det = det+Vec1(i)*(Vec2(j)*Vec3(k)-Vec2(k)*Vec3(j))
       end do

       return
    End Function Determ_V_R
 
End Submodule Determ_Vec
