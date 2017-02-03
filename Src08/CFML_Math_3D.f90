!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!----
!!---- MODULE: CFML_Math_3D
!!----   INFO: Simple mathematics general utilities for 3D Systems
!!----
!!---- HISTORY
!!----    Update: 04/03/2011
!!----
!!
 Module CFML_Math_3D
    !---- Use Modules ----!
    Use CFML_DefPar,       only: CP, SP, DP, PI, TO_RAD, TO_DEG, Err_CFML, Err_CFML_Mess, &
                                 Init_Err_CFML

    Use CFML_Math_General, only: CosD, SinD

    implicit none

    private

    !---- List of public functions ----!
    public :: Cross_Product, Determ_3x3, Determ_Vec, Invert_Array3x3, Mat_Cross,  &
              Rotate_OX, Rotate_OY, Rotate_OZ, Tensor_Product, Vec_Length

    !---- List of public subroutines ----!
    public :: Get_Cart_From_Cylin, Get_Cart_From_Spher, Get_Cylin_from_Cart, Get_Spher_from_Cart,   &
              Matrix_DiagEigen, Resolv_Sist_1X2, Resolv_Sist_1X3, Resolv_Sist_2X2, Resolv_Sist_2X3, &
              Resolv_Sist_3X3, Set_Eps

    !-------------------!
    !---- Variables ----!
    !-------------------!
    real(kind=cp)  ::  eps=0.00001_cp    ! Epsilon value

    !---- Interfaces - Overlapp ----!
    Interface  Cross_Product
       Module Procedure Cross_product_sp
       Module Procedure Cross_product_dp
       Module Procedure Cross_product_in
       Module Procedure Cross_product_cmpl_sp
       Module Procedure Cross_product_cmpl_dp
    End Interface

    Interface  Determ_3x3
       Module Procedure Determ_A_I
       Module Procedure Determ_A_R
    End Interface

    Interface  Determ_Vec
       Module Procedure Determ_V_I
       Module Procedure Determ_V_R
    End Interface

    Interface  Invert_Array3x3
       Module Procedure Invert_sp
       Module Procedure Invert_dp
    End Interface

    Interface  Mat_Cross
       Module Procedure Mat_Cross_sp
       Module Procedure Mat_Cross_dp
       Module Procedure Mat_Cross_in
       Module Procedure Mat_Cross_cmpl_sp
       Module Procedure Mat_Cross_cmpl_dp
    End Interface

    Interface  Tensor_Product
       Module Procedure Tensor_product_sp
       Module Procedure Tensor_product_dp
       Module Procedure Tensor_product_in
       Module Procedure Tensor_product_cmpl_sp
       Module Procedure Tensor_product_cmpl_dp
    End Interface

    Interface  Get_Cart_from_Cylin
       Module Procedure Get_Cart_from_Cylin_dp
       Module Procedure Get_Cart_from_Cylin_sp
    End Interface

    Interface  Get_Cart_from_Spher
       Module Procedure Get_Cart_from_Spher_dp
       Module Procedure Get_Cart_from_Spher_sp
    End Interface

    Interface  Get_Cylin_from_Cart
       Module Procedure Get_Cylin_from_Cart_dp
       Module Procedure Get_Cylin_from_Cart_sp
    End Interface

    Interface  Get_Spher_from_Cart
       Module Procedure Get_Spher_from_Cart_dp
       Module Procedure Get_Spher_from_Cart_sp
    End Interface

 Contains

    !!--++
    !!--++ Function  Cross_Product_cmpl_dp(U,V) Result(W)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of the complex vectors u and v
    !!--++    Vectors, w = u x v, are given in cartesian components.
    !!--++
    !!--++ Update: June - 2012
    !!
    Function Cross_Product_cmpl_dp(u,v) Result(w)
       !---- Argument ----!
       complex(kind=dp), dimension(3), intent( in) :: u  ! Vector 1
       complex(kind=dp), dimension(3), intent( in) :: v  ! Vector 2
       complex(kind=dp), dimension(3)              :: w  ! u x v

       w(1)=u(2)*v(3)-u(3)*v(2)
       w(2)=u(3)*v(1)-u(1)*v(3)
       w(3)=u(1)*v(2)-u(2)*v(1)

       return
    End Function Cross_Product_cmpl_dp

    !!--++
    !!--++ Function  Cross_Product_cmpl_sp(U,V) Result(W)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of the complex vectors u and v
    !!--++    Vectors, w = u x v, are given in cartesian components.
    !!--++
    !!--++ Update: June - 2012
    !!
    Function Cross_Product_cmpl_sp(u,v) Result(w)
       !---- Argument ----!
       complex(kind=sp), dimension(3), intent( in) :: u   ! Vector 1
       complex(kind=sp), dimension(3), intent( in) :: v   ! Vector 2
       complex(kind=sp), dimension(3)              :: w   ! u x v

       w(1)=u(2)*v(3)-u(3)*v(2)
       w(2)=u(3)*v(1)-u(1)*v(3)
       w(3)=u(1)*v(2)-u(2)*v(1)

       return
    End Function Cross_Product_cmpl_sp

    !!--++
    !!--++ Function  Cross_Product_dp(U,V) Result(W)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of vectors u and v
    !!--++    Vectors, w= u x v, are given in cartesian components.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Cross_Product_dp(u,v) Result(w)
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: u   ! Vector 1
       real(kind=dp), dimension(3), intent( in) :: v   ! Vector 2
       real(kind=dp), dimension(3)              :: w   ! u x v

       w(1)=u(2)*v(3)-u(3)*v(2)
       w(2)=u(3)*v(1)-u(1)*v(3)
       w(3)=u(1)*v(2)-u(2)*v(1)

       return
    End Function Cross_Product_dp

    !!--++
    !!--++ Function  Cross_Product_in(U,V) Result(W)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of integer vectors u and v
    !!--++    In the indices are givent w.r.t the direct lattice, the cross product
    !!--++    are indices w.r.t. reciprocal lattice and viceversa.
    !!--++
    !!--++ Update: November - 2008
    !!
    Function Cross_Product_in(u,v) Result(w)
       !---- Argument ----!
       integer, dimension(3), intent( in) :: u    ! Vector 1
       integer, dimension(3), intent( in) :: v    ! Vector 2
       integer, dimension(3)              :: w    ! u x v

       w(1)=u(2)*v(3)-u(3)*v(2)  ! i  j   k !
       w(2)=u(3)*v(1)-u(1)*v(3)  !u1  u2  u3! = (u2.v3 - u3.v2)i + (v1.u3 - u1.v3)j + (u1.v2-u2.v1)k
       w(3)=u(1)*v(2)-u(2)*v(1)  !v1  v2  v3!

       return
    End Function Cross_Product_in

    !!--++
    !!--++ Function  Cross_Product_sp(U,V) Result(W)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of vectors u and v
    !!--++    Vectors, w= u x v, are given in cartesian components.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Cross_Product_sp(u,v) Result(w)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: u   ! Vector 1
       real(kind=sp), dimension(3), intent( in) :: v   ! Vector 2
       real(kind=sp), dimension(3)              :: w   ! u x v

       w(1)=u(2)*v(3)-u(3)*v(2)
       w(2)=u(3)*v(1)-u(1)*v(3)
       w(3)=u(1)*v(2)-u(2)*v(1)

       return
    End Function Cross_Product_sp

    !!--++
    !!--++ Function Determ_A_I(Array)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of an integer 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_A_I(Array) Result(determ)
       !---- Argument ----!
       integer, dimension(3,3), intent(in) :: array
       integer                             :: determ

       determ=array(1,1)*array(2,2)*array(3,3)+ &
              array(2,1)*array(3,2)*array(1,3)+ &
              array(1,2)*array(2,3)*array(3,1)- &
              array(1,3)*array(2,2)*array(3,1)- &
              array(1,1)*array(3,2)*array(2,3)- &
              array(1,2)*array(2,1)*array(3,3)

       return
    End Function Determ_A_I

    !!--++
    !!--++ Function Determ_A_R(Array)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_A_R(Array) Result (determ)
       !---- Argument ----!
       real(kind=cp), dimension(3,3), intent(in) :: array
       real(kind=cp)                             :: determ

       determ=array(1,1)*array(2,2)*array(3,3)+ &
              array(2,1)*array(3,2)*array(1,3)+ &
              array(1,2)*array(2,3)*array(3,1)- &
              array(1,3)*array(2,2)*array(3,1)- &
              array(1,1)*array(3,2)*array(2,3)- &
              array(1,2)*array(2,1)*array(3,3)

       return
    End Function Determ_A_R

    !!--++
    !!--++ Function Determ_V_I(A,B,C)
    !!--++    integer, dimension(3), intent(in) :: a,b,c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of the components of three vectors
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_V_I(Vec1,Vec2,Vec3) Result(det)
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

    !!--++
    !!--++ Function Determ_V_R(A,B,C)
    !!--++    real(kin=cp), dimension(3), intent(in) :: a,b,c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of the components of three vectors
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_V_R(Vec1,Vec2,Vec3) Result(det)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: Vec1,Vec2,Vec3
       real(kind=cp)                           :: det

       !---- Local variables ----!
       integer :: i,j,k

       det = 0.0
       do i = 1,3
          j = i+1
          if (j == 4) j = 1
          k = 6-i-j
          det = det+Vec1(i)*(Vec2(j)*Vec3(k)-Vec2(k)*Vec3(j))
       end do

       return
    End Function Determ_V_R

    !!--++
    !!--++ Funcion Invert_Dp(Array) Result(b)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate de inverse of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Invert_Dp(array) Result(b)
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

    !!--++
    !!--++ Funcion Invert_Sp(Array) Result(b)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate de inverse of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Invert_Sp(array) Result(b)
       !---- Arguments ----!
       real(kind=sp),dimension(3,3), intent(in) :: array
       real(kind=sp),dimension(3,3)             :: b

       !---- Local variables ----!
       real(kind=sp)  :: dmat

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
          b=0.0
       end if

       return
    End Function Invert_Sp

    !!--++
    !!--++ Function  Mat_Cross_cmpl_dp(Vec) Result(M)
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Function Mat_Cross_cmpl_dp(Vec) Result(M)
       !---- Argument ----!
       complex(kind=dp), dimension(3), intent( in) :: Vec
       complex(kind=dp), dimension(3,3)            :: M

       M = reshape( (/  (0.0_dp,0.0_dp),   -vec(3),         vec(2),  &
                            vec(3),   (0.0_dp,0.0_dp),     -vec(1),  &
                           -vec(2),           vec(1),   (0.0_dp,0.0_dp)/),(/3,3/))
       return
    End Function Mat_Cross_cmpl_dp

    !!--++
    !!--++ Function  Mat_Cross_cmpl_sp(Vec) Result(M)
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Function Mat_Cross_cmpl_sp(Vec) Result(M)
       !---- Argument ----!
       complex(kind=sp), dimension(3), intent( in) :: Vec
       complex(kind=sp), dimension(3,3)            :: M

       M = reshape( (/  (0.0_sp,0.0_sp),   -vec(3),         vec(2),  &
                            vec(3),   (0.0_sp,0.0_sp),     -vec(1),  &
                           -vec(2),           vec(1),   (0.0_sp,0.0_sp)/),(/3,3/))
       return
    End Function Mat_Cross_cmpl_sp

    !!--++
    !!--++ Function  Mat_Cross_dp(Vec) Result(M)
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Function Mat_Cross_dp(Vec) Result(M)
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: Vec
       real(kind=dp), dimension(3,3)            :: M

       M = reshape( (/ 0.0_dp,   -vec(3),     vec(2),  &
                        vec(3),    0.0_dp,   -vec(1),  &
                       -vec(2),     vec(1),    0.0_dp/),(/3,3/))
       return
    End Function Mat_Cross_dp

    !!--++
    !!--++ Function  Mat_Cross_in(Vec) Result(M)
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Function Mat_Cross_in(Vec) Result(M)
       !---- Argument ----!
       integer, dimension(3), intent( in) :: Vec
       integer, dimension(3,3)            :: M

       M = reshape( (/   0,    -vec(3),    vec(2),  &
                        vec(3),    0,     -vec(1),  &
                       -vec(2),   vec(1),     0 /),(/3,3/))
       return
    End Function Mat_Cross_in

    !!--++
    !!--++ Function  Mat_Cross_sp(Vec) Result(M)
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Function Mat_Cross_sp(Vec) Result(M)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: Vec
       real(kind=sp), dimension(3,3)            :: M

       M = reshape( (/ 0.0_sp, -vec(3),    vec(2),  &
                        vec(3),  0.0_sp,  -vec(1),  &
                       -vec(2),   vec(1),   0.0_sp/),(/3,3/))
       return
    End Function Mat_Cross_sp

    !!----
    !!---- Function Rotate_OX(Vec,Angle) Result (Vec)
    !!----
    !!----    X Rotation. Positive rotation is counter-clockwise
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OX(Vec,Angle) Result(Rvec)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: Vec      ! Vector
       real(kind=cp),               intent(in) :: angle    ! Angle
       real(kind=cp), dimension(3)             :: Rvec     ! Vector rotated

       !---- Variables locales ----!
       real(kind=cp), dimension(3,3)           :: rot

       rot(1,1)=  1.0
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

    !!----
    !!---- Function Rotate_OY(Vec,Angle) Result (Rec)
    !!----
    !!----    Y Rotation. Positive rotation is counter-clockwise
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OY(Vec,Angle) Result(Rvec)
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

    !!----
    !!---- Function Rotate_OZ(Vec,Angle) Result (RVec)
    !!----
    !!----    Z Rotation
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OZ(Vec,Angle) Result(Rvec)
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

    !!----
    !!---- Function  Tensor_Product_cmpl_dp(Vec1,Vec2) Result(W)
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Function Tensor_Product_cmpl_dp(Vec1,Vec2) Result(w)
       !---- Argument ----!
       complex(kind=dp), dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2
       complex(kind=dp), dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2

       !---- Local Arguments ----!
       complex(kind=dp), dimension(3,3)            :: mu,mv

       mu=0.0_dp;  mv=0.0_dp
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)

       return
    End Function Tensor_Product_cmpl_dp

    !!----
    !!---- Function  Tensor_Product_cmpl_sp(Vec1,Vec2) Result(W)
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Function Tensor_Product_cmpl_sp(Vec1,Vec2) Result(w)
       !---- Argument ----!
       complex(kind=sp), dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2
       complex(kind=sp), dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2

       !---- Local Arguments ----!
       complex(kind=sp), dimension(3,3)            :: mu,mv

       mu=0.0_sp;  mv=0.0_sp
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)
       return
    End Function Tensor_Product_cmpl_sp

    !!----
    !!---- Function  Tensor_Product_dp(Vec1,Vec2) Result(W)
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Function Tensor_Product_dp(Vec1,Vec2) Result(w)
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: Vec1,Vec2
       real(kind=dp), dimension(3,3)            :: w

       !---- Local Variables ----!
       real(kind=dp), dimension(3,3)            :: mu,mv

       mu=0.0_dp;  mv=0.0_dp
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)

       return
    End Function Tensor_Product_dp

    !!----
    !!---- Function  Tensor_Product_in(Vec1,Vec2) Result(W)
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Function Tensor_Product_in(Vec1,Vec2) Result(w)
       !---- Argument ----!
       integer, dimension(3), intent( in) :: Vec1,Vec2
       integer, dimension(3,3)            :: w

       !---- Local Variables ----!
       integer, dimension(3,3)            :: mu,mv
       mu=0;  mv=0
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)

       return
    End Function Tensor_Product_in

    !!----
    !!---- Function  Tensor_Product_sp(Vec1,Vec2) Result(W)
    !!----
    !!----    Calculates the tensor product of vectors Vec1 and Vec2
    !!----
    !!---- Updated: June - 2012
    !!
    Function Tensor_Product_sp(Vec1,Vec2) Result(w)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: Vec1,Vec2
       real(kind=sp), dimension(3,3)            :: w

       !---- Local Variables ----!
       real(kind=sp), dimension(3,3)            :: mu,mv

       mu=0.0_sp;  mv=0.0_sp
       mu(:,1)=Vec1
       mv(1,:)=Vec2
       w=matmul(mu,mv)

       return
    End Function Tensor_Product_sp

    !!----
    !!---- Function Vec_Length(G,Vec) Result(c)
    !!----
    !!----    Length of vector B when A is the Crystallographic
    !!----    to orthogonal matrix length=c
    !!----
    !!---- Update: February - 2005
    !!
    Function Vec_Length(G,Vec) Result(c)
       !---- Arguments ----!
       real(kind=cp), intent(in)  , dimension(3,3)       :: G      ! Metric array
       real(kind=cp), intent(in)  , dimension(3  )       :: Vec    ! Vector
       real(kind=cp)                                     :: c      ! Length of Vector

       !---- Local variables ----!
       integer                     :: i,j
       real(kind=cp), dimension(3) :: v

       v=0.0
       do i = 1,3
          do j = 1,3
             v(i) = v(i)+G(i,j)*Vec(j)
          end do
       end do

       c = sqrt(v(1)**2+v(2)**2+v(3)**2)

       return
    End Function Vec_Length

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!--++
    !!--++ Subroutine  Get_Cart_from_Cylin_dp(rho,Phi,z,Xo,Mode)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from cylindrical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Cylin_dp(rho,Phi,z,Xo,Mode)
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

    !!--++
    !!--++ Subroutine  Get_Cart_from_Cylin_sp(rho,Phi,z,Xo,Mode)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from cylindrical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Cylin_sp(rho,Phi,z,Xo,Mode)
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

    !!--++
    !!--++ Subroutine Get_Cart_from_Spher_dp(r,Theta,Phi,Xo,Mode)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from spherical coordinates.
    !!--++    Theta is the azimutal angle
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Spher_dp(r,Theta,Phi,Xo,Mode)
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

    !!--++
    !!--++ Subroutine Get_Cart_from_Spher_sp(r,Theta,Phi,Xo,Mode)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from spherical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Spher_sp(r,Theta,Phi,Xo,Mode)
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

    !!--++
    !!--++ Subroutine  Get_Cylin_from_Cart_dp(Xo,rho,Phi,zeta,Mode)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the cylindrical coordinates from Cartesian coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cylin_from_Cart_dp(Xo,rho,Phi,z,Mode)
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

    !!--++
    !!--++ Subroutine  Get_Cylin_from_Cart_sp(Xo,rho,Phi,z,Mode)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the cylindrical coordinates from Cartesian coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cylin_from_Cart_sp(Xo,rho,Phi,z,Mode)
       !---- Arguments ----!
       real(kind=sp), dimension(3),intent( in)           ::  xo
       real(kind=sp),              intent(out)           ::  rho
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

    !!--++
    !!--++ Subroutine Get_Spher_from_Cart_dp(Xo,Ss,Theta,Phi,Mode)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the spheric coordinates from rectangular coordinates
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Spher_from_Cart_dp(xo,r,theta,phi,mode)
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

    !!--++
    !!--++ Subroutine Get_Spher_from_Cart_sp(Xo,R,Theta,Phi,Mode)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the spheric coordinates from rectangular coordinates
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Spher_from_Cart_sp(xo,r,theta,phi,mode)
       !---- Arguments ----!
       real(kind=sp), intent( in), dimension(3)   :: xo
       real(kind=sp), intent(out)                 :: r
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

    !!----
    !!---- Subroutine Matrix_DiagEigen(Array, EigenValues, EigenVec)
    !!----
    !!----    Diagonalize the matrix Array, put eigenvalues in EigenValues and
    !!----    eigenvectors in EigenVec
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Matrix_DiagEigen(Array,EigenValues,EigenVec)
       !---- Arguments ----!
       real(kind=cp), intent(in)  , dimension(3,3)    :: array
       real(kind=cp), intent(out) , dimension(3)      :: EigenValues
       real(kind=cp), intent(out) , dimension(3,3)    :: EigenVec

       !---- Local Variables ----!
       integer, parameter            :: n=3
       integer                       :: i, j, k, itmax, nm1, ip1, iter
       real(kind=cp), dimension(3)   :: u
       real(kind=cp), dimension(3,3) :: e
       real(kind=cp), parameter      :: eps1=1.e-7 , eps2=1.e-7 , eps3=1.e-7
       real(kind=cp)                 :: sigma1, offdsq, p, q, spq, csa, sna
       real(kind=cp)                 :: holdik, holdki, sigma2

       !> Init
       call init_Err_CFML()

       nm1=n-1
       itmax=50
       do i=1,n
          do j=1,n
             e(i,j)=array(i,j)
             EigenVec(i,j)=0.0
             if (j < i) e(i,j)=0.0
          end do
       end do
       sigma1=0.0
       offdsq=0.0

       do i=1,n
          sigma1=sigma1+e(i,i)**2
          EigenVec(i,i)=1.0
          ip1=i+1
          if (i >= n) exit
          do j=ip1,n
             offdsq=offdsq+e(i,j)**2
          end do
       end do

       do iter=1,itmax
          do i=1,nm1
             ip1=i+1
             do j=ip1,n
                q=abs(e(i,i)-e(j,j))
                if (q <= eps1) then
                   csa=1.0/sqrt(2.0)
                   sna=csa
                else
                   if (abs(e(i,j)) <= eps2) then
                      e(i,j)=0.0
                      cycle
                   end if
                   p=2.0*e(i,j)*q/(e(i,i)-e(j,j))
                   spq=sqrt(p*p+q*q)
                   csa=sqrt((1.0+q/spq)/2.0)
                   sna=p/(2.0*csa*spq)
                end if
                do k=1,n
                   holdki=EigenVec(k,i)
                   EigenVec(k,i)=holdki*csa + EigenVec(k,j)*sna
                   EigenVec(k,j)=holdki*sna - EigenVec(k,j)*csa
                end do
                do k=i,n
                   if (k > j) then
                      holdik=e(i,k)
                      e(i,k)=csa*holdik + sna*e(j,k)
                      e(j,k)=sna*holdik - csa*e(j,k)
                   else
                      u(k)=e(i,k)
                      e(i,k)=csa*u(k)+sna*e(k,j)
                      if (k /= j) cycle
                      e(j,k)=sna*u(k)-csa*e(j,k)
                   end if
                end do
                u(j)=sna*u(i)-csa*u(j)
                do k=1,j
                   if (k <= i)  then
                      holdki=e(k,i)
                      e(k,i)=csa*holdki+sna*e(k,j)
                      e(k,j)=sna*holdki-csa*e(k,j)
                   else
                      e(k,j)=sna*u(k)-csa*e(k,j)
                   end if
                end do
                e(i,j)=0.0
             end do
          end do
          sigma2=0.0
          do i=1,n
             EigenValues(i)=e(i,i)
             sigma2=sigma2+EigenValues(i)*EigenValues(i)
          end do
          if (1.0-sigma1/sigma2 <= eps3) return
          sigma1=sigma2
       end do

       Err_CFML =.true.
       Err_CFML_Mess=" Convergence not reached in diagonalization "

       return
    End Subroutine Matrix_DiagEigen

    !!----
    !!---- Subroutine Resolv_Sist_1X2(W,T,Ts,X,Ix)
    !!--<<
    !!----              w11 x1 + w12 x2  = t1
    !!----              x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_1x2(w,t,ts,x,ix)
       !---- Arguments ----!
       integer, dimension(2),         intent( in) :: w      ! Input vector
       real(kind=cp),                 intent( in) :: t      ! Input value
       real(kind=cp), dimension(2),   intent(out) :: ts     ! Fixed value solution
       real(kind=cp), dimension(2),   intent(out) :: x      ! Fixed value for x,y
       integer, dimension(2),         intent(out) :: ix     ! 1: x, 2: y, 3: z

       !> Init
       call init_Err_CFML()
       ts = 0.0
       x  = 1.0
       ix = 0

       !> Both are zeros
       if ( all(w == 0)) then
          if (abs(t) < eps) then
             ix(1)=1
             ix(2)=2
          else
             Err_CFML=.true.
             Err_CFML_Mess="Inconsistent solution (1x2)"
          end if
          return
       end if

       !> Any is zero
       if (any(w == 0)) then
          if ( w(1) == 0 ) then
             ix(1)=1
             ts(2)=t/real(w(2))
              x(2)=0.0
          else
             ts(1)=t/real(w(1))
              x(1)=0.0
             ix(2)=2
          end if
       else
          ix(1)=1
          ts(2)=t/real(w(2))
           x(2)=-real(w(1))/real(w(2))
          ix(2)=1
       end if

       return
    End Subroutine Resolv_Sist_1x2

    !!----
    !!---- Subroutine Resolv_Sist_1X3(W,T,Ts,X,Ix)
    !!--<<
    !!----               w11 x1 + w12 x2 + w13 x3 = t1
    !!----               x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_1x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer, dimension(3),         intent( in) :: w      ! Input vector
       real(kind=cp),                 intent( in) :: t      ! Input value
       real(kind=cp), dimension(3),   intent(out) :: ts     ! Fixed value solution
       real(kind=cp), dimension(3),   intent(out) :: x      ! Fixed value for x,y,z
       integer, dimension(3),         intent(out) :: ix     ! 1: x, 2: y, 3: z

       !---- Local Variables ----!
       integer               :: i, zeros
       integer, dimension(2) :: w1
       integer, dimension(2) :: ix1
       real(kind=cp), dimension(2)    :: ts1
       real(kind=cp), dimension(2)    :: x1

       !> Init
       call init_Err_CFML()
       ts = 0.0
       x  = 1.0
       ix = 0

       zeros=0
       do i=1,3
          if (w(i) == 0) zeros=zeros+1
       end do
       select case (zeros)
          case (3)
             if (abs(t) < eps) then
                do i=1,3
                   ix(i)=i
                end do
             else
                Err_CFML=.true.
                Err_CFML_Mess="Inconsistent solution (1 x 3)"
             end if

          case (2)
             do i=1,3
                if (w(i) /= 0) then
                   ts(i)=t/real(w(i))
                   x(i) =0.0
                else
                   ix(i)=i
                end if
             end do

          case (1)
             do i=1,3
                if (w(i) == 0) exit
             end do
             select case (i)
                case (1)
                   w1=w(2:3)

                case (2)
                   w1(1)=w(1)
                   w1(2)=w(3)

                case (3)
                   w1=w(1:2)
             end select
             call resolv_sist_1x2(w1,t,ts1,x1,ix1)
             select case (i)
                case (1)
                   ix(1)  = 1
                   ts(2:3)= ts1
                   x(2:3) = x1
                   if (ix1(1)==1) ix(2)=2
                   if (ix1(1)==2) ix(2)=3
                   if (ix1(2)==1) ix(3)=2
                   if (ix1(2)==2) ix(3)=3

                  case (2)
                     ix(2)= 2
                     ts(1)= ts1(1)
                     ts(3)= ts1(2)
                     x(1) = x1(1)
                     x(3) = x1(2)
                     if (ix1(1)==1) ix(1)=1
                     if (ix1(1)==2) ix(1)=3
                     if (ix1(2)==1) ix(3)=1
                     if (ix1(2)==2) ix(3)=3

                  case (3)
                     ix(3)  = 3
                     ts(1:2)= ts1
                     x(1:2) = x1
                     ix(1:2)= ix1
               end select

          case (0)
             Err_CFML=.true.
             Err_CFML_Mess="Inconsistent case ax+by+cz=t (1x3)"
       end select

       return
    End Subroutine Resolv_Sist_1x3

    !!----
    !!---- Subroutine Resolv_Sist_2X2(W,T,Ts,X,Ix)
    !!--<<
    !!----                 w11 x1 + w12 x2  = t1
    !!----                 w21 x1 + w22 x2  = t2
    !!----                 x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_2x2(w,t,ts,x,ix)
       !---- Arguments ----!
       integer, dimension(2,2),       intent( in) :: w       ! Input vector
       real(kind=cp),dimension(2),    intent( in) :: t       ! Input value
       real(kind=cp),dimension(2),    intent(out) :: ts      ! Fixed value solution
       real(kind=cp),dimension(2),    intent(out) :: x       ! Fixed value for x,y
       integer, dimension(2),         intent(out) :: ix      ! 1: x, 2: y, 3: z

       !---- Local Variables ----!
       integer                 :: i,deter
       integer, dimension(2)   :: zeros,colum
       real(kind=cp)           :: rden, rnum

       !> Init
       call init_Err_CFML()
       ts    = 0.0
       x     = 1.0
       ix    = 0

       deter = w(1,1)*w(2,2) - w(1,2)*w(2,1)
       rden=real(deter)
       if (deter /= 0) then
          !---- X(1) ----!
          rnum=t(1)*w(2,2) - w(1,2)*t(2)
          ts(1)=rnum/rden

          !---- X(2) ----!
          rnum=w(1,1)*t(2) - t(1)*w(2,1)
          ts(2)=rnum/rden

          x =0.0

       else                        ! Singular Matrix
          !---- Are there zero rows? ----!
          zeros=0
          do i=1,2
             if (w(i,1) == 0 .and. w(i,2) == 0 )  zeros(i)=1
          end do
          select case (sum(zeros))
             case (2)
                if (abs(t(1)) <= eps .and. abs(t(2)) <= eps) then
                   ix(1)=1
                   ix(2)=2
                else
                   Err_CFML=.true.
                   Err_CFML_Mess="Inconsistent solution (2x2)"
                end if

             case (1)
                do i=1,2
                   if (zeros(i) == 0) exit
                end do
                call resolv_sist_1x2(w(i,:),t(i),ts,x,ix)

             case (0)
                !---- Are there zero columns? ----!
                colum=0
                do i=1,2
                   if (w(1,i) == 0 .and. w(2,i) == 0 ) colum(i)=1
                end do
                select case (sum(colum))
                   case (1)
                      do i=1,2
                         if (colum(i) == 0) exit
                      end do
                      if (w(1,i) /= 0) then
                         ts(i)=t(1)/real(w(1,i))
                      else
                         ts(i)=t(2)/real(w(2,i))
                      end if
                      x(i)=0.0
                      if (i == 1) then
                         ix(2)=2
                      else
                         ix(1)=1
                      end if

                   case (0)
                      call resolv_sist_1x2(w(1,:),t(1),ts,x,ix)

                end select
          end select
       end if

       return
    End Subroutine Resolv_Sist_2x2

    !!----
    !!---- Subroutine Resolv_Sist_2X3(W,T,Ts,X,Ix)
    !!----               w11 x1 + w12 x2 + w13 x3 = t1
    !!----               w21 x1 + w22 x2 + w23 x3 = t2
    !!----               x_sol(i)= ts(i) + x(i) ix(i)
    !!----
    !!----   Update: February - 2005
    !!
    Subroutine Resolv_Sist_2x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer, dimension(2,3),          intent( in) :: w     ! Input vector
       real(kind=cp), dimension(2),      intent( in) :: t     ! Input value
       real(kind=cp), dimension(3),      intent(out) :: ts    ! Fixed value solution
       real(kind=cp), dimension(3),      intent(out) :: x     ! Fixed value for x,y
       integer, dimension(3),            intent(out) :: ix    ! 1: x, 2: y, 3: z

       !---- Local Variables ----!
       integer                 :: i, j
       integer, dimension(2)   :: fila
       integer, dimension(2)   :: ix1
       integer, dimension(3)   :: colum
       integer, dimension(2,2) :: w1
       integer, dimension(2,3) :: wm
       integer, dimension(2)   :: wc
       real(kind=cp)                    :: tc
       real(kind=cp), dimension(2)      :: tm
       real(kind=cp), dimension(2)      :: ts1, x1

       !---- Initialize ----!
       call init_Err_CFML()
       ts    = 0.0
       x     = 1.0
       ix    = 0

       !---- Are there zero columns? ----!
       colum=0
       do i=1,3
            if (all(w(:,i) == 0)) colum(i)=1
       end do
       select case (sum(colum))
          case (3)
             if (abs(t(1)) <= eps .and. abs(t(2)) <= eps) then
                do i=1,3
                   ix(i)=i
                end do
             else
                Err_CFML=.true.
                Err_CFML_Mess="Inconsistent solution in (2x3)"
             end if

          case (2)
             do i=1,3
                if (colum(i) == 0) exit
             end do
             if (w(1,i) /= 0) then
                ts(i)=t(1)/real(w(1,i))
             else
                ts(i)=t(2)/real(w(2,i))
             end if
             x(i)=0.0
             select case (i)
                case (1)
                   ix(2)=2
                   ix(3)=3

                case (2)
                   ix(1)=1
                   ix(3)=3

                case (3)
                   ix(1)=1
                   ix(2)=2
             end select

          case (1)
             do i=1,3
                if (colum(i) == 1) exit
             end do
             select case (i)
                case (1)
                   w1=w(:,2:3)

                case (2)
                   w1(1,1)=w(1,1)
                   w1(1,2)=w(1,3)
                   w1(2,1)=w(2,1)
                   w1(2,2)=w(2,3)

                case (3)
                   w1=w(:,1:2)
             end select
             call resolv_sist_2x2(w1,t,ts1,x1,ix1)
             select case (i)
                case (1)
                   ix(1)  = 1
                   ts(2:3)= ts1
                   x (2:3)= x1
                   if (ix1(1) == 1) ix(2)=2
                   if (ix1(1) == 2) ix(2)=3
                   if (ix1(2) == 1) ix(3)=2
                   if (ix1(2) == 2) ix(3)=3

                case (2)
                   ix(2)=2
                   ts(1)=ts1(1)
                   ts(3)=ts1(2)
                   x(1) = x1(1)
                   x(3) = x1(2)
                   if (ix1(1) == 1) ix(1)=1
                   if (ix1(1) == 2) ix(1)=3
                   if (ix1(2) == 1) ix(3)=1
                   if (ix1(2) == 2) ix(3)=3

                case (3)
                   ix(3)  = 3
                   ts(1:2)= ts1
                   x (1:2)= x1
                   ix(1:2)= ix1
             end select

          case (0)
             !---- Are there zeros in any element of rows? ----!
             fila = 0
             do i=1,2
                if (all(w(i,:)==0)) fila(i)=1
             end do
             select case (sum(fila))
                case (1)
                   if (w(1,1) /= 0) then
                      call resolv_sist_1x3(w(1,:),t(1),ts,x,ix)
                   else
                      call resolv_sist_1x3(w(2,:),t(2),ts,x,ix)
                   end if

                case (0)
                   fila = 0
                   wm   = w
                   tm   = t
                   !---- Are there zeros in any element of rows? ----!
                   do i=1,2
                      do j=1,3
                         if (w(i,j)==0) fila(i)=fila(i)+1
                      end do
                   end do
                   if ( fila(2) > fila(1) ) then
                      wm(1,:)=w(2,:)
                      wm(2,:)=w(1,:)
                      tm(1)  =t(2)
                      tm(2)  =t(1)
                          j  =fila(1)
                      fila(1)=fila(2)
                      fila(2)=j
                   end if
                   select case (fila(1))
                      case (2)
                         do i=1,3
                            if (wm(1,i) /= 0) exit
                         end do
                         ts(i)=tm(1)/real(wm(1,i))
                         x(i)=0.0
                         select case (i)
                            case (1)
                               wc(1)=wm(2,2)
                               wc(2)=wm(2,3)
                               tc=tm(2)-(wm(2,1)*ts(i))

                            case (2)
                               wc(1)=wm(2,1)
                               wc(2)=wm(2,3)
                               tc=tm(2)-(wm(2,2)*ts(i))

                            case (3)
                               wc(1)=wm(2,1)
                               wc(2)=wm(2,2)
                               tc=tm(2)-(wm(2,3)*ts(i))
                         end select
                         call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                         select case(i)
                            case (1)
                               ts(2:3)=ts1
                                x(2:3)=x1
                                if (ix1(1)==1) ix(2)=2
                                if (ix1(1)==2) ix(2)=3
                                if (ix1(2)==1) ix(3)=2
                                if (ix1(2)==2) ix(3)=3

                            case (2)
                               ts(1)=ts1(1)
                               ts(3)=ts1(2)
                                x(1)=x1(1)
                                x(3)=x1(2)
                                if (ix1(1)==1) ix(1)=1
                                if (ix1(1)==2) ix(1)=3
                                if (ix1(2)==1) ix(3)=1
                                if (ix1(2)==2) ix(3)=3

                            case (3)
                               ts(1:2)=ts1
                                x(1:2)=x1
                               ix(1:2)=ix1
                         end select

                      case (1)
                         do i=1,3
                            if (wm(1,i) == 0) exit
                         end do
                         select case (fila(2))
                            case (1)
                               do j=1,3
                                  if (wm(2,j) == 0) exit
                               end do
                               select case (i)
                                  case (1)             ! 0 en w(1,1)
                                     select case (j)
                                        case (2)
                                           wc(1)=-wm(2,1)/wm(2,3)
                                           wc(2)= wm(1,2)/wm(1,3)
                                           tc=tm(1)/real(wm(1,3)) - tm(2)/real(wm(2,3))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1:2)=ts1
                                           x(1:2) =x1
                                           ix(1:2)=ix1
                                           if (ix(1) == 0) then
                                              ts(3)=tm(2)/real(wm(2,3)) - ts(1)*wm(2,1)/real(wm(2,3))
                                              x(3)=0.0
                                           else
                                              if (ix(2) == 0) then
                                                 ts(3)=tm(1)/real(wm(1,3)) - ts(2)*wm(1,2)/real(wm(1,3))
                                                 x(3)=0.0
                                              else
                                                 ts(3)=tm(2)/real(wm(2,3))
                                                 x(3)=-real(wm(2,1))/real(wm(2,3))
                                                 ix(3)=1

                                                 ts(2)=tc/real(wc(2))
                                                 x(2) =-real(wc(1))/real(wc(2))
                                                 ix(2)=1
                                              end if
                                           end if

                                        case (3)
                                           wc(1)=-wm(2,1)/wm(2,2)
                                           wc(2)= wm(1,3)/wm(1,2)
                                           tc=tm(1)/real(wm(1,2)) - tm(2)/real(wm(2,2))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1)=ts1(1)
                                           ts(3)=ts1(2)
                                           x(1) =x1(1)
                                           x(3) =x1(2)
                                           if (ix1(1) == 1) ix(1)=1
                                           if (ix1(1) == 2) ix(1)=3
                                           if (ix1(2) == 1) ix(3)=1
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(1) == 0) then
                                              ts(2)=tm(2)/real(wm(2,2)) - ts(1)*wm(2,1)/real(wm(2,2))
                                              x(2)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(2)=tm(1)/real(wm(1,2)) - ts(3)*wm(1,3)/real(wm(1,2))
                                                 x(2)=0.0
                                              else
                                                 ts(2)=tm(2)/real(wm(2,2))
                                                 x(3)=-real(wm(2,1))/real(wm(2,2))
                                                 ix(2)=1

                                                 ts(3)=tc/real(wc(2))
                                                 x(3) =-real(wc(1))/real(wc(2))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select

                                  case (2)             ! 0 en w(1,2)
                                     select case (j)
                                        case (1)
                                           wc(1)= wm(1,1)/wm(1,3)
                                           wc(2)=-wm(2,2)/wm(2,3)
                                           tc=tm(1)/real(wm(1,3)) - tm(2)/real(wm(2,3))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1:2)=ts1
                                           x(1:2) =x1
                                           ix(1:2)=ix1
                                           if (ix(1) == 0) then
                                              ts(3)=tm(1)/real(wm(1,3)) - ts(1)*wm(1,1)/real(wm(1,3))
                                              x(3)=0.0
                                           else
                                              if (ix(2) == 0) then
                                                 ts(3)=tm(2)/real(wm(2,3)) - ts(2)*wm(2,2)/real(wm(2,3))
                                                 x(3)=0.0
                                              else
                                                 ts(3)=tm(1)/real(wm(1,3))
                                                 x(3)=-real(wm(1,1))/real(wm(1,3))
                                                 ix(3)=1

                                                 ts(2)=tc/real(wc(2))
                                                 x(2) = -real(wc(1))/real(wc(2))
                                                 ix(2)= 1
                                              end if
                                           end if

                                        case (3)
                                           wc(1)=-wm(2,2)/wm(2,1)
                                           wc(2)= wm(1,3)/wm(1,1)
                                           tc=tm(1)/real(wm(1,1)) - tm(2)/real(wm(2,1))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(2:3)=ts1
                                           x(2:3) =x1
                                           if (ix1(1) == 1) ix(2)=2
                                           if (ix1(1) == 2) ix(2)=3
                                           if (ix1(2) == 1) ix(3)=2
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(2) == 0) then
                                              ts(1)=tm(2)/real(wm(2,1)) - ts(2)*wm(2,2)/real(wm(2,1))
                                              x(1)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(1)=tm(1)/real(wm(1,1)) - ts(3)*wm(1,3)/real(wm(1,1))
                                                 x(1)=0.0
                                              else
                                                 ix(1)=1

                                                 ts(2)=tm(2)/real(wm(2,2))
                                                 x(2) =-real(wm(2,1))/real(wm(2,2))
                                                 ix(2)=1

                                                 ts(3)=tm(1)/real(wm(1,3))
                                                 x(3) =-real(wm(1,1))/real(wm(1,3))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select

                                  case (3)             ! 0 en w(1,3)
                                     select case (j)
                                        case (1)
                                           wc(1)= wm(1,1)/wm(1,2)
                                           wc(2)=-wm(2,3)/wm(2,2)
                                           tc=tm(1)/real(wm(1,2)) - tm(2)/real(wm(2,2))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1)=ts1(1)
                                           ts(3)=ts1(2)
                                           x(1) =x1(1)
                                           x(3) =x1(2)
                                           if (ix1(1) == 1) ix(1)=1
                                           if (ix1(1) == 2) ix(1)=3
                                           if (ix1(2) == 1) ix(3)=1
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(1) == 0) then
                                              ts(2)=tm(1)/real(wm(1,2)) - ts(1)*wm(1,1)/real(wm(1,2))
                                              x(2)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(2)=tm(2)/real(wm(2,2)) - ts(3)*wm(2,3)/real(wm(2,2))
                                                 x(2)=0.0
                                              else
                                                 ts(2)=tm(1)/real(wm(1,2))
                                                 x(2) =-real(wm(1,1))/real(wm(1,2))
                                                 ix(2)=1

                                                 ts(3)=tc/real(wc(2))
                                                 x(3) =-real(wc(1))/real(wc(2))
                                                 ix(3)=1
                                              end if
                                           end if

                                        case (2)
                                           wc(1)= wm(1,2)/wm(1,1)
                                           wc(2)=-wm(2,3)/wm(2,1)
                                           tc=tm(1)/real(wm(1,1)) - tm(2)/real(wm(2,1))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(2:3)=ts1
                                           x(2:3) =x1
                                           if (ix1(1) == 1) ix(2)=2
                                           if (ix1(1) == 2) ix(2)=3
                                           if (ix1(2) == 1) ix(3)=2
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(2) == 0) then
                                              ts(1)=tm(1)/real(wm(1,1)) - ts(2)*wm(1,2)/real(wm(1,1))
                                              x(1)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(1)=tm(2)/real(wm(2,1)) - ts(3)*wm(2,3)/real(wm(2,1))
                                                 x(1)=0.0
                                              else
                                                 ix(1)=1

                                                 ts(2)=tm(1)/real(wm(1,2))
                                                 x(2) =-real(wm(1,1))/real(wm(1,2))
                                                 ix(2)=1

                                                 ts(3)=tm(2)/real(wm(2,3))
                                                 x(3) =-real(wm(2,1))/real(wm(2,3))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select
                               end select

                            case (0)
                               select case (i)
                                  case (1)
                                     wc(1)=wm(2,1)
                                     wc(2)=wm(2,2)- wm(2,3)*wm(1,2)/wm(1,3)
                                     tc=tm(2)-tm(1)*wm(2,3)/real(wm(1,3))
                                     call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                     ts(1:2)=ts1
                                     x(1:2)=x1
                                     ix(1:2)=ix1
                                     if (ix(2) == 0) then
                                        ts(3)=tm(1)/real(wm(1,3)) - ts(2)*real(wm(1,2))/real(wm(1,3))
                                        x(3)=0.0
                                     else
                                        ix(1)=1

                                        ts(2)=(tm(2) - tm(1)*wm(2,3)/real(wm(1,3))) / &
                                              (real(wm(2,2)) - real(wm(2,3)*wm(1,2))/real(wm(1,3)) )
                                        x(2) =-real(wm(2,1)) / &
                                              (real(wm(2,2)) - real(wm(2,3)*wm(1,2))/real(wm(1,3)) )
                                        ix(2)=1

                                        ts(3)= tm(1)/real(wm(1,3)) - (real(wm(1,2))/real(wm(1,3)))*ts(2)
                                        x(3) =- (real(wm(1,2))/real(wm(1,3)))*x(2)
                                        ix(3)=1
                                     end if

                                  case (2)
                                     wc(1)=wm(2,1)-wm(2,3)*wm(1,1)/wm(1,3)
                                     wc(2)=wm(2,2)
                                     tc=tm(2)-tm(1)*wm(2,3)/real(wm(1,3))
                                     call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                    ts(1:2)=ts1
                                    x(1:2)=x1
                                    ix(1:2)=ix1
                                    if (ix(1) == 0) then
                                       ts(3)=tm(1)/real(wm(1,3)) - ts(1)*real(wm(1,1))/real(wm(1,3))
                                       x(3)=0.0
                                    else
                                       ix(1)=1

                                       ts(2)=(tm(2) - tm(1)*wm(2,3)/real(wm(1,3)))/real(wm(2,2))
                                       x(2) =(real(wm(1,1)*wm(2,3))/real(wm(1,3)) - real(wm(2,1)))/real(wm(2,2))
                                       ix(2)=1

                                       ts(3)=tm(1)/real(wm(1,3))
                                       x(3) =-real(wm(1,1))/real(wm(1,3))
                                       ix(3)=1
                                    end if

                                 case (3)
                                    wc(1)=wm(2,1)-wm(1,1)*wm(2,2)/wm(1,2)
                                    wc(2)=wm(2,3)
                                    tc=tm(2)-tm(1)*wm(2,2)/real(wm(1,2))
                                    call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                    ts(1)=ts1(1)
                                    ts(3)=ts1(2)
                                    x(1)=x1(1)
                                    x(3)=x1(2)
                                    if (ix1(1) == 1) ix(1)=1
                                    if (ix1(1) == 2) ix(1)=3
                                    if (ix1(2) == 1) ix(3)=1
                                    if (ix1(2) == 2) ix(3)=3
                                    if (ix(1) == 0) then
                                       ts(2)=tm(1)/real(wm(1,2)) - ts(1)*real(wm(1,1))/real(wm(1,2))
                                       x(2)=0.0
                                    else
                                       ix(1) =1

                                       ts(2)=tm(1)/real(wm(1,2))
                                       x(2) =-real(wm(1,1))/real(wm(1,2))
                                       ix(2)=1

                                       ts(3)=(tm(2) - tm(1)*wm(2,2)/real(wm(1,2)))/real(wm(2,3))
                                       x(3) =(real(wm(1,1)*wm(2,2))/real(wm(1,2)) - real(wm(2,1)))/real(wm(2,3))
                                       ix(3)=1
                                    end if
                               end select
                         end select

                      case (0)
                         call resolv_sist_1x3(wm(1,:),tm(1),ts,x,ix)
                   end select

             end select
       end select

       return
    End Subroutine Resolv_Sist_2x3

    !!----
    !!---- Subroutine Resolv_Sist_3X3(W,T,Ts,X,Ix)
    !!--<<
    !!----              w11 x1 + w12 x2 + w13 x3 = t1
    !!----              w21 x1 + w22 x2 + w23 x3 = t2
    !!----              w31 x1 + w32 x2 + w33 x3 = t3
    !!----              x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_3x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer, dimension(3,3),          intent(in) :: w       ! Input vector
       real(kind=cp), dimension(3),      intent(in) :: t       ! Input value
       real(kind=cp), dimension(3),      intent(out):: ts      ! Fixed value solution
       real(kind=cp), dimension(3),      intent(out):: x       ! Fixed value for x,y
       integer, dimension(3),            intent(out):: ix      ! 1: x, 2: y, 3: z

       !---- Local variables ----!
       integer                 :: i,j,deter
       integer, dimension(3)   :: fila
       integer, dimension(3,3) :: w1
       integer, dimension(2,3) :: wm
       real(kind=cp)                    :: rnum, rden
       real(kind=cp), dimension(3)      :: t1
       real(kind=cp), dimension(2)      :: tm
       real(kind=cp),dimension(3,3)     :: rw

       !---- Initialize ----!
       call init_Err_CFML()
       ts  = 0.0
       x   = 1.0
       ix  = 0

       deter=Determ_3x3(w)
       rden=real(deter)

       if (deter /= 0) then
          !---- X(1) ----!
          rw=real(w)
          rw(:,1)=t
          rnum=Determ_3x3(rw)
          ts(1)=rnum/rden

          !---- X(2) ----!
          rw=real(w)
          rw(:,2)=t
          rnum=Determ_3x3(rw)
          ts(2)=rnum/rden

          !---- X(3) ----!
          rw=real(w)
          rw(:,3)=t
          rnum=Determ_3x3(rw)
          ts(3)=rnum/rden

          x=0.0

       else                     !  Singular Matrix
          !---- Are there zero rows? ----!
          fila=0
          do i=1,3
             if (all(w(i,:) == 0)) fila(i)=1
          end do
          select case (sum(fila))
             !---- All values are zeros ----!
             case (3)
                if (all(abs(t) < eps)) then
                   do i=1,3
                      ix(i)=i
                   end do
                else
                   Err_CFML=.true.
                   Err_CFML_Mess="Inconsistent system (3 x 3)"
                end if

             !---- Two rows with zeroes ----!
             case (2)
                do i=1,3
                   if (fila(i) == 0) exit
                end do
                call resolv_sist_1x3(w(i,:),t(i),ts,x,ix)

             !---- One row with zeroes ----!
             case (1)
                do i=1,3
                   if (fila(i) == 1) exit
                end do
                select case(i)
                   case (1)
                      wm(1,:)=w(2,:)
                      wm(2,:)=w(3,:)
                      tm=t(2:3)

                   case (2)
                      wm(1,:)=w(1,:)
                      wm(2,:)=w(3,:)
                      tm(1)=t(1)
                      tm(2)=t(3)

                   case (3)
                      wm(1,:)=w(1,:)
                      wm(2,:)=w(2,:)
                      tm=t(1:2)

                end select
                call resolv_sist_2x3(wm,tm,ts,x,ix)

             !---- Non zero rows ----!
             case (0)
                w1=w
                t1=t

                !---- Are there 2 rows proportional? ----!
                do i=1,3
                   if ( abs(w1(1,i)) > abs(w1(2,i)) ) then
                      if (w1(2,i) /= 0) then
                         j=w1(1,i)/w1(2,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(2,1) == w1(1,1) .and. j*w1(2,2) == w1(1,2) .and. &
                             j*w1(2,3) == w1(1,3) ) then
                            w1(1,:)=w1(2,:)
                            t1(1)  =t1(2)
                            exit
                         end if
                      end if
                   else
                      if (w1(1,i) /= 0) then
                         j=w1(2,i)/w1(1,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(1,1) == w1(2,1) .and. j*w1(1,2) == w1(2,2) .and. &
                             j*w1(1,3) == w1(2,3) ) then
                            w1(2,:)=w1(1,:)
                            t1(2)  =t1(1)
                            exit
                         end if
                      end if
                   end if
                end do

                do i=1,3
                   if ( abs(w1(1,i)) > abs(w1(3,i)) ) then
                      if (w1(3,i) /= 0) then
                         j=w1(1,i)/w1(3,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(3,1) == w1(1,1) .and. j*w1(3,2) == w1(1,2) .and. &
                             j*w1(3,3) == w1(1,3) ) then
                            w1(1,:)=w1(3,:)
                            t1(1)  =t1(3)
                            exit
                         end if
                      end if
                   else
                      if (w1(1,i) /= 0) then
                         j=w1(3,i)/w1(1,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(1,1) == w1(3,1) .and. j*w1(1,2) == w1(3,2) .and. &
                             j*w1(1,3) == w1(3,3) ) then
                            w1(3,:)=w1(1,:)
                            t1(3)  =t1(1)
                            exit
                         end if
                      end if
                   end if
                end do

                do i=1,3
                   if ( abs(w1(2,i)) > abs(w1(3,i)) ) then
                      if (w1(3,i) /= 0) then
                         j=w1(2,i)/w1(3,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(3,1) == w1(2,1) .and. j*w1(3,2) == w1(2,2) .and. &
                             j*w1(3,3) == w1(2,3) ) then
                            w1(2,:)=w1(3,:)
                            t1(2)  =t1(3)
                            exit
                         end if
                      end if
                   else
                      if (w1(2,i) /= 0) then
                         j=w1(3,i)/w1(2,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(2,1) == w1(3,1) .and. j*w1(2,2) == w1(3,2) .and. &
                             j*w1(2,3) == w1(3,3) ) then
                            w1(3,:)=w1(2,:)
                            t1(3)  =t1(2)
                            exit
                         end if
                      end if
                   end if
                end do

                !---- Are there 3 rows equal? ----!
                if ( (w1(1,1) == w1(2,1)) .and. (w1(1,1) == w1(3,1)) .and. &
                     (w1(1,2) == w1(2,2)) .and. (w1(1,2) == w1(3,2)) .and. &
                     (w1(1,3) == w1(2,3)) .and. (w1(1,3) == w1(3,3)) ) then

                   call resolv_sist_1x3(w1(1,:),t1(1),ts,x,ix)

                !---- Are there 2 rows equal? ----!
                elseif( (w1(1,1) == w1(2,1)) .and. (w1(1,2) == w1(2,2)) .and. &
                        (w1(1,3) == w1(2,3)) ) then

                   call resolv_sist_2x3(w1(2:3,:),t1(2:3),ts,x,ix)

                elseif( (w1(1,1) == w1(3,1)) .and. (w1(1,2) == w1(3,2)) .and. &
                        (w1(1,3) == w1(3,3)) ) then

                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                elseif( (w1(2,1) == w1(3,1)) .and. (w1(2,2) == w1(3,2)) .and. &
                        (w1(2,3) == w1(3,3)) ) then

                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                !---- Are linear combinations? ----!
                else
                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                end if

          end select
       end if

       return
    End Subroutine Resolv_Sist_3x3

    !!----
    !!---- Subroutine Set_Eps(Neweps)
    !!----
    !!----    Sets global EPS to the value "neweps"
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Eps(Neweps)
       !---- Arguments ----!
       real(kind=cp), optional, intent( in) :: neweps  ! Sets global EPS to the value "neweps"

       if (present(neweps)) then
          eps=neweps
       else
          eps=0.00001_cp
       end if

       return
    End Subroutine Set_Eps

 End Module CFML_Math_3D
