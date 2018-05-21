!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Get_Cylin_from_Cart
 Contains
 
    !!--++   SUBROUTINE GET_CYLIN_FROM_CART_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the cylindrical coordinates from Cartesian coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Get_Cylin_from_Cart_dp(Xo,rho,Phi,z,Mode)    
       !---- Arguments ----!
       real(kind=dp), dimension(3),intent( in)           ::  xo   ! Cartesian coordinatates
       real(kind=dp),              intent(out)           ::  rho  ! Cylindrical coordinates
       real(kind=dp),              intent(out)           ::  phi
       real(kind=dp),              intent(out)           ::  z
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       integer :: j

       z=xo(3)
       if( abs(xo(2)) > eps .or. abs(xo(1)) > eps) then
          phi=atan2(xo(2),xo(1))
       else
          phi= 0.0_dp
       end if

       rho=0.0_dp
       do j=1,2
          rho=rho+xo(j)*xo(j)
       end do
       rho=sqrt(rho)

       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") phi=phi*to_deg
       end if

       return
    End Subroutine Get_Cylin_from_Cart_dp
 
    !!--++   SUBROUTINE GET_CYLIN_FROM_CART_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the cylindrical coordinates from Cartesian coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Get_Cylin_from_Cart_sp(Xo,rho,Phi,z,Mode)    
       !---- Arguments ----!
       real(kind=sp), dimension(3),intent( in)           ::  xo    ! Cartesian
       real(kind=sp),              intent(out)           ::  rho   ! Cylindrical
       real(kind=sp),              intent(out)           ::  phi
       real(kind=sp),              intent(out)           ::  z
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       integer :: j

       z=xo(3)
       if( abs(xo(2)) > eps .or. abs(xo(1)) > eps) then
          phi=atan2(xo(2),xo(1))
       else
          phi= 0.0_sp
       end if
       rho=0.0_sp
       do j=1,2
          rho=rho+xo(j)*xo(j)
       end do
       rho=sqrt(rho)

       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") phi=phi*to_deg
       end if

       return
    End Subroutine Get_Cylin_from_Cart_sp
 
   
End Submodule Get_Cylin_from_Cart
