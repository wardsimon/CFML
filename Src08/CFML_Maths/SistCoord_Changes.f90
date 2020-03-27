!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Maths) CFML_Math_011
 Contains
    !!----
    !!---- GET_CART_FROM_SPHER
    !!----    Determine the Cartesian coordinates from spherical coordinates.
    !!----    Theta is the azimutal angle
    !!----
    !!---- 04/04/2019
    !!
    Module Pure Function Get_Cart_from_Spher(SphCoord,Mode) Result(CarCoord)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent( in) :: SphCoord ! Coordinates (R,Theta;Phi)
       character(len=*), optional,  intent( in) :: mode     ! If "D" the angles are in degrees, otherwise radians is considered
       real(kind=cp), dimension(3)              :: CarCoord ! Cartesian coordinates

       !---- Local Variables ----!
       real(kind=cp) :: ph,th

       th=SphCoord(2)
       ph=SphCoord(3)
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             th=SphCoord(2)*TO_RAD
             ph=SphCoord(3)*TO_RAD
          end if
       end if
       CarCoord(1)=SphCoord(1)*cos(ph)*sin(th)
       CarCoord(2)=SphCoord(1)*sin(ph)*sin(th)
       CarCoord(3)=SphCoord(1)*cos(th)

       return
    End Function Get_Cart_from_Spher

    !!----
    !!---- GET_CART_FROM_CYLIN
    !!----    Determine the Cartesian coordinates from cylindrical coordinates.
    !!----
    !!---- 04/04/2019
    !!
    Module Pure Function Get_Cart_from_Cylin(CilCoord,Mode) Result(CarCoord)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent( in) ::  CilCoord ! Coordinates rho,phi,zeta
       character(len=*), optional,  intent( in) ::  mode     ! "D" angles in degrees, otherwise in radians
       real(kind=cp), dimension(3)              ::  CarCoord ! Cartesian coordinates

       !---- Local Variables ----!
       real(kind=cp) :: ph

       ph=CilCoord(2)
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") ph=CilCoord(2)*TO_RAD
       end if
       CarCoord(1)=CilCoord(1)*cos(ph)
       CarCoord(2)=CilCoord(1)*sin(ph)
       CarCoord(3)=CilCoord(3)

       return
    End Function Get_Cart_from_Cylin

    !!----
    !!---- GET_CYLIN_FROM_CART
    !!----    Determine the cylindrical coordinates from Cartesian coordinates.
    !!----
    !!---- 04/04/2019
    !!
    Module Pure Function Get_Cylin_from_Cart(CarCoord, Mode) Result(CilCoord)
       !---- Arguments ----!
       real(kind=cp), dimension(3),intent(in) ::  CarCoord   ! Cartesian coordinatates
       character(len=*), optional, intent(in) ::  mode
       real(kind=cp), dimension(3)            ::  CilCoord   ! Cylindrical coordinates

       !---- Local Variables ----!
       integer :: j

       CilCoord(3)=CarCoord(3)
       if( abs(CarCoord(2)) > eps .or. abs(CarCoord(1)) > eps) then
          CilCoord(2)=atan2(CarCoord(2),CarCoord(1))
       else
          CilCoord(2)= 0.0_cp
       end if

       CilCoord(1)=0.0_cp
       do j=1,2
          CilCoord(1)=CilCoord(1)+CarCoord(j)*CarCoord(j)
       end do
       CilCoord(1)=sqrt(CilCoord(1))

       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") CilCoord(2)=CilCoord(2)*TO_DEG
       end if

       return
    End Function Get_Cylin_from_Cart

    !!----
    !!---- GET_SPHER_FROM_CART
    !!----    Determine the spheric coordinates from rectangular coordinates
    !!----
    !!---- 04/04/2019
    !!
    Module Pure Function Get_Spher_from_Cart(CarCoord,mode) Result(SphCoord)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: CarCoord ! Cartesian
       character(len=*), optional,  intent(in) :: mode
       real(kind=cp), dimension(3)             :: SphCoord ! Spherical

       !---- Local Variables ----!
       integer :: j

       SphCoord(1)=0.0_cp
       do j=1,3
          SphCoord(1)=SphCoord(1)+CarCoord(j)*CarCoord(j)
       end do
       SphCoord(1)=sqrt(SphCoord(1))

       if (SphCoord(1) > 0.0_cp) then
          SphCoord(2)=CarCoord(3)/SphCoord(1)
          if (abs(SphCoord(2)) > 1.0_cp) then
             SphCoord(2)=sign(1.0_cp,SphCoord(2))
          end if
          SphCoord(2)=acos(SphCoord(2))
          if (abs(SphCoord(2)) < eps .or. abs(SphCoord(2)-pi) < eps) then
             SphCoord(3)=0.0_cp
          else
             SphCoord(3)=atan2(CarCoord(2),CarCoord(1))
          end if
       else
          SphCoord(2)=0.0_cp
          SphCoord(3)=0.0_cp
       end if
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             SphCoord(2)=SphCoord(2)*TO_DEG
             SphCoord(3)=SphCoord(3)*TO_DEG
          end if
       end if

       return
    End Function Get_Spher_from_Cart

    !!----
    !!---- GET_SPHER_FROM_CYLIN
    !!----    Determine the spheric coordinates from cylinder coordinates
    !!----
    !!---- 04/04/2019
    !!
    Module Pure Function Get_Spher_from_Cylin(CilCoord,mode) Result(SphCoord)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: CilCoord ! Cylinder
       character(len=*), optional,  intent(in) :: mode
       real(kind=cp), dimension(3)             :: SphCoord ! Spherical

       !---- Local Variables ----!
       real(kind=cp), dimension(3) :: CarCoord

       if (present(mode)) then
          CarCoord=Get_Cart_from_Cylin(CilCoord,mode)
          SphCoord=Get_Spher_from_Cart(CarCoord,mode)
       else
          CarCoord=Get_Cart_from_Cylin(CilCoord)
          SphCoord=Get_Spher_from_Cart(CarCoord,mode)
       end if

       return
    End Function Get_Spher_from_Cylin

    !!----
    !!---- GET_CYLIN_FROM_SPHER
    !!----    Determine the spheric coordinates from cylinder coordinates
    !!----
    !!---- 04/04/2019
    !!
    Module Pure Function Get_Cylin_from_Spher(SphCoord,mode) Result(CilCoord)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: SphCoord ! Cylinder
       character(len=*), optional,  intent(in) :: mode
       real(kind=cp), dimension(3)             :: CilCoord ! Spherical

       !---- Local Variables ----!
       real(kind=cp), dimension(3) :: CarCoord

       if (present(mode)) then
          CarCoord=Get_Cart_from_Spher(SphCoord,mode)
          CilCoord=Get_Cylin_from_Cart(CarCoord,mode)
       else
          CarCoord=Get_Cart_from_Spher(SphCoord,mode)
          CilCoord=Get_Cylin_from_Cart(CarCoord,mode)
       end if

       return
    End Function Get_Cylin_from_Spher

End Submodule CFML_Math_011
