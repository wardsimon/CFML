!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Get_Spher_from_Cart
 Contains
 
    !!--++ SUBROUTINE GET_SPHER_FROM_CART_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the spheric coordinates from rectangular coordinates
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Get_Spher_from_Cart_dp(xo,r,theta,phi,mode)    
       !---- Arguments ----!
       real(kind=dp), intent( in), dimension(3)   :: xo      ! Cartesian
       real(kind=dp), intent(out)                 :: r      ! Spherical
       real(kind=dp), intent(out)                 :: theta
       real(kind=dp), intent(out)                 :: phi
       character(len=*), intent(in), optional     :: mode

       !---- Local Variables ----!
       integer :: j

       r=0.0_dp
       do j=1,3
          r=r+xo(j)*xo(j)
       end do
       r=sqrt(r)
       if (r > 0.0_dp) then
          theta=xo(3)/r
          if (abs(theta) > 1.0_dp) then
             theta=sign(1.0_dp,theta)
          end if
          theta=acos(theta)
          if (abs(theta) < eps .or. abs(theta-pi) < eps) then
             phi=0.0_dp
          else
             phi=atan2(xo(2),xo(1))
          end if
       else
          theta=0.0_dp
          phi=0.0_dp
       end if
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             theta=theta*to_deg
             phi=phi*to_deg
          end if
       end if

       return
    End Subroutine Get_Spher_from_Cart_dp
 
    !!--++  SUBROUTINE GET_SPHER_FROM_CART_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the spheric coordinates from rectangular coordinates
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Get_Spher_from_Cart_sp(xo,r,theta,phi,mode)    
       !---- Arguments ----!
       real(kind=sp), intent( in), dimension(3)   :: xo      ! Cartesian coordinates
       real(kind=sp), intent(out)                 :: r       ! Spherical
       real(kind=sp), intent(out)                 :: theta
       real(kind=sp), intent(out)                 :: phi
       character(len=*), intent(in), optional     :: mode

       !---- Local Variables ----!
       integer :: j

       r=0.0_sp
       do j=1,3
          r=r+xo(j)*xo(j)
       end do
       r=sqrt(r)
       if (r > 0.0_sp) then
          theta=xo(3)/r
          if (abs(theta) > 1.0_sp) then
             theta=sign(1.0_sp,theta)
          end if
          theta=acos(theta)
          if (abs(theta) < eps .or. abs(theta-pi) < eps) then
             phi=0.0_sp
          else
             phi=atan2(xo(2),xo(1))
          end if
       else
          theta=0.0_sp
          phi=0.0_sp
       end if
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             theta=theta*to_deg
             phi=phi*to_deg
          end if
       end if

       return
    End Subroutine Get_Spher_from_Cart_sp
 
   
End Submodule Get_Spher_from_Cart
