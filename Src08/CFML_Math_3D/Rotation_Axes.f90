!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Rotation_Axes
 Contains
 
    !!---- FUNCTION ROTATE_OX
    !!----
    !!----    X Rotation. Positive rotation is counter-clockwise
    !!----
    !!---- Update: February - 2005
    !!
    Module Function Rotate_OX(Vec,Angle) Result(Rvec)    
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: Vec      ! Vector
       real(kind=cp),               intent(in) :: angle    ! Angle
       real(kind=cp), dimension(3)             :: Rvec     ! Vector rotated

       !---- Variables locales ----!
       real(kind=cp), dimension(3,3)           :: rot

       rot(1,1)=  1.0_cp
       rot(2,1)=  0.0_cp
       rot(3,1)=  0.0_cp

       rot(1,2)=  0.0_cp
       rot(2,2)=  cosd(angle)
       rot(3,2)=  sind(angle)

       rot(1,3)=  0.0_cp
       rot(2,3)=  -sind(angle)
       rot(3,3)=  cosd(angle)

       Rvec=matmul(rot,vec)

       return
    End Function Rotate_OX
    
    !!----  FUNCTION ROTATE_OY
    !!----
    !!----    Y Rotation. Positive rotation is counter-clockwise
    !!----
    !!---- Update: February - 2005
    !!
    Module Function Rotate_OY(Vec,Angle) Result(Rvec)    
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: Vec     ! Vector
       real(kind=cp),               intent(in) :: angle   ! Angle
       real(kind=cp), dimension(3)             :: Rvec    ! Vector rotated

       !---- Variables locales ----!
       real(kind=cp), dimension(3,3)           :: rot

       rot(1,1)=  cosd(angle)
       rot(2,1)=  0.0_cp
       rot(3,1)=  -sind(angle)

       rot(1,2)=  0.0_cp
       rot(2,2)=  1.0_cp
       rot(3,2)=  0.0_cp

       rot(1,3)= sind(angle)
       rot(2,3)= 0.0_cp
       rot(3,3)= cosd(angle)

       Rvec=matmul(rot,vec)

       return
    End Function Rotate_OY
    
    !!----  FUNCTION ROTATE_OZ
    !!----
    !!----    Z Rotation
    !!----
    !!---- Update: February - 2005
    !!
    Module Function Rotate_OZ(Vec,Angle) Result(Rvec)    
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: Vec      ! Vector
       real(kind=cp),               intent(in) :: angle    ! Angle
       real(kind=cp), dimension(3)             :: Rvec     ! Vector rotated

       !---- Variables locales ----!
       real(kind=cp), dimension(3,3)           :: rot

       rot(1,1)=  cosd(angle)
       rot(2,1)=  sind(angle)
       rot(3,1)=  0.0_cp

       rot(1,2)=  -sind(angle)
       rot(2,2)=  cosd(angle)
       rot(3,2)=  0.0_cp

       rot(1,3)=  0.0_cp
       rot(2,3)=  0.0_cp
       rot(3,3)=  1.0_cp

       Rvec=matmul(rot,Vec)

       return
    End Function Rotate_OZ
 
End Submodule Rotation_Axes
