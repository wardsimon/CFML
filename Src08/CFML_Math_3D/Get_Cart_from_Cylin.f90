!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Get_Cart_from_Cylin
 Contains
 
    !!--++ SUBROUTINE GET_CART_FROM_CYLIN_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from cylindrical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Get_Cart_from_Cylin_dp(rho,Phi,z,Xo,Mode)    
       !---- Arguments ----!
       real(kind=dp),              intent( in)           ::  rho    ! Coordinates rho,phi,zeta
       real(kind=dp),              intent( in)           ::  phi
       real(kind=dp),              intent( in)           ::  z
       real(kind=dp), dimension(3),intent(out)           ::  xo     ! Cartesian coordinates
       character(len=*),           intent( in), optional ::  mode   ! "D" angles in degrees, otherwise in radians

       !---- Local Variables ----!
       real(kind=dp) :: ph

       ph=phi
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") ph=phi*to_rad
       end if
       xo(1)=rho*cos(ph)
       xo(2)=rho*sin(ph)
       xo(3)=z

       return
    End Subroutine Get_Cart_from_Cylin_dp
 
    !!--++   SUBROUTINE GET_CART_FROM_CYLIN_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from cylindrical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Get_Cart_from_Cylin_sp(rho,Phi,z,Xo,Mode)    
       real(kind=sp),              intent( in)           ::  rho  ! Coordinates rho,phi,zeta
       real(kind=sp),              intent( in)           ::  phi
       real(kind=sp),              intent( in)           ::  z
       real(kind=sp), dimension(3),intent(out)           ::  xo   ! Cartesian coordinates
       character(len=*),           intent( in), optional ::  mode ! "D" angles in degrees, otherwise in radians

       !---- Local Variables ----!
       real(kind=sp) :: ph

       ph=phi
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") ph=phi*to_rad
       end if
       xo(1)=rho*cos(ph)
       xo(2)=rho*sin(ph)
       xo(3)=z

       return
    End Subroutine Get_Cart_from_Cylin_sp
 
   
End Submodule Get_Cart_from_Cylin
