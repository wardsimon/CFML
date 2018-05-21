!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Get_Cart_from_Spher
 Contains
 
    !!--++  SUBROUTINE GET_CART_FROM_SPHER_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from spherical coordinates.
    !!--++    Theta is the azimutal angle
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Get_Cart_from_Spher_dp(r,Theta,Phi,Xo,Mode)    
       !---- Arguments ----!
       real(kind=dp),              intent( in)           :: r       ! Coordinates (R,Theta;Phi)
       real(kind=dp),              intent( in)           :: Theta
       real(kind=dp),              intent( in)           :: phi
       real(kind=dp), dimension(3),intent(out)           :: xo      ! Cartesian coordinates
       character(len=*),           intent( in), optional :: mode    ! If "D" the angles are in degrees, otherwise radians is considered

       !---- Local Variables ----!
       real(kind=dp) :: ph,th

       ph=Phi
       th=Theta
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             ph=Phi*to_rad
             th=Theta*to_rad
          end if
       end if
       xo(1)=r*cos(ph)*sin(th)
       xo(2)=r*sin(ph)*sin(th)
       xo(3)=r*cos(th)

       return
    End Subroutine Get_Cart_from_Spher_dp
 
    !!--++  SUBROUTINE GET_CART_FROM_SPHER_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from spherical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Get_Cart_from_Spher_sp(r,Theta,Phi,Xo,Mode)    
       !---- Arguments ----!
       real(kind=sp),              intent( in)           :: r       ! Coordinates R, Theta, Phi
       real(kind=sp),              intent( in)           :: Theta
       real(kind=sp),              intent( in)           :: phi
       real(kind=sp), dimension(3),intent(out)           :: xo      ! Cartesian Coordinates
       character(len=*),           intent( in), optional :: mode    ! If "D" then angles are given in degrees.

       !---- Local Variables ----!
       real(kind=sp) :: ph,th

       ph=Phi
       th=Theta
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             ph=Phi*to_rad
             th=Theta*to_rad
          end if
       end if
       xo(1)=r*cos(ph)*sin(th)
       xo(2)=r*sin(ph)*sin(th)
       xo(3)=r*cos(th)

       return
    End Subroutine Get_Cart_from_Spher_sp
 
   
End Submodule Get_Cart_from_Spher
