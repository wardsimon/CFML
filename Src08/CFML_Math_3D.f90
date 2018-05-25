!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2018  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!----               Ross Angel         (University of Pavia)
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
!!----
!!
 Module CFML_Math_3D
    !---- Use Modules ----!
    Use CFML_GlobalDeps, only: SP, DP, CP, PI, TO_RAD, TO_DEG, &
                               Err_CFML 


    implicit none

    private

    !---- List of public functions ----!
    public :: Rotate_OX, Rotate_OY, Rotate_OZ, Vec_length

    !---- List of public overloaded procedures: functions ----!
    public :: Cross_Product, Determ_3x3, Determ_Vec, Invert_Array3x3, Mat_Cross, &
              Tensor_Product

    !---- List of public subroutines ----!
    public :: Set_Eps,  Matrix_DiagEigen, Resolv_Sist_1X2, Resolv_Sist_1X3, Resolv_Sist_2X2, &
              Resolv_Sist_2X3, Resolv_Sist_3X3

    !---- List of public overloaded procedures: subroutines ----!
    public :: Get_Cart_From_Cylin, Get_Cart_From_Spher, Get_Cylin_from_Cart, Get_Spher_from_Cart


    !---------------------!
    !---- Definitions ----!
    !---------------------!
    real(kind=cp),  private  ::  eps=0.00001_cp   ! Epsilon value


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

 
    !> Interface Zone 
    Interface
 
       Module Function Cross_Product_cmpl_dp(u,v) Result(w)    
          !---- Argument ----!       
          complex(kind=dp), dimension(3), intent( in) :: u  ! Vector 1       
          complex(kind=dp), dimension(3), intent( in) :: v  ! Vector 2       
          complex(kind=dp), dimension(3)              :: w  ! u x v       
       End Function Cross_Product_cmpl_dp    
 
       Module Function Cross_Product_cmpl_sp(u,v) Result(w)    
          !---- Argument ----!       
          complex(kind=sp), dimension(3), intent( in) :: u   ! Vector 1       
          complex(kind=sp), dimension(3), intent( in) :: v   ! Vector 2       
          complex(kind=sp), dimension(3)              :: w   ! u x v       
       End Function Cross_Product_cmpl_sp    
 
       Module Function Cross_Product_dp(u,v) Result(w)    
          !---- Argument ----!       
          real(kind=dp), dimension(3), intent( in) :: u   ! Vector 1       
          real(kind=dp), dimension(3), intent( in) :: v   ! Vector 2       
          real(kind=dp), dimension(3)              :: w   ! u x v       
       End Function Cross_Product_dp    
 
       Module Function Cross_Product_in(u,v) Result(w)    
          !---- Argument ----!       
          integer, dimension(3), intent( in) :: u    ! Vector 1       
          integer, dimension(3), intent( in) :: v    ! Vector 2       
          integer, dimension(3)              :: w    ! u x v       
       End Function Cross_Product_in    
 
       Module Function Cross_Product_sp(u,v) Result(w)    
          !---- Argument ----!       
          real(kind=sp), dimension(3), intent( in) :: u   ! Vector 1       
          real(kind=sp), dimension(3), intent( in) :: v   ! Vector 2       
          real(kind=sp), dimension(3)              :: w   ! u x v       
       End Function Cross_Product_sp    
 
       Module Function Determ_A_I(Array) Result(determ)    
          !---- Argument ----!       
          integer, dimension(3,3), intent(in) :: array       
          integer                             :: determ       
       End Function Determ_A_I    
 
       Module Function Determ_A_R(Array) Result (determ)    
          !---- Argument ----!       
          real(kind=cp), dimension(3,3), intent(in) :: array       
          real(kind=cp)                             :: determ       
       End Function Determ_A_R    
 
       Module Function Determ_V_I(Vec1,Vec2,Vec3) Result(det)    
          !---- Arguments ----!       
          integer, dimension(3), intent(in) :: Vec1,Vec2,Vec3       
          integer                           :: det       
       End Function Determ_V_I    
 
       Module Function Determ_V_R(Vec1,Vec2,Vec3) Result(det)    
          !---- Arguments ----!       
          real(kind=cp), dimension(3), intent(in) :: Vec1,Vec2,Vec3       
          real(kind=cp)                           :: det       
       End Function Determ_V_R    
 
       Module Function Invert_Dp(array) Result(b)    
          !---- Arguments ----!       
          real(kind=dp),dimension(3,3), intent(in) :: array       
          real(kind=dp),dimension(3,3)             :: b       
       End Function Invert_Dp    
 
       Module Function Invert_Sp(array) Result(b)    
          !---- Arguments ----!       
          real(kind=sp),dimension(3,3), intent(in) :: array       
          real(kind=sp),dimension(3,3)             :: b       
       End Function Invert_Sp    
 
       Module Function Mat_Cross_cmpl_dp(Vec) Result(M)    
          !---- Argument ----!       
          complex(kind=dp), dimension(3), intent( in) :: Vec       
          complex(kind=dp), dimension(3,3)            :: M       
       End Function Mat_Cross_cmpl_dp    
 
       Module Function Mat_Cross_cmpl_sp(Vec) Result(M)    
          !---- Argument ----!       
          complex(kind=sp), dimension(3), intent( in) :: Vec       
          complex(kind=sp), dimension(3,3)            :: M       
       End Function Mat_Cross_cmpl_sp    
 
       Module Function Mat_Cross_dp(Vec) Result(M)    
          !---- Argument ----!       
          real(kind=dp), dimension(3), intent( in) :: Vec       
          real(kind=dp), dimension(3,3)            :: M       
       End Function Mat_Cross_dp    
 
       Module Function Mat_Cross_in(Vec) Result(M)    
          !---- Argument ----!       
          integer, dimension(3), intent( in) :: Vec       
          integer, dimension(3,3)            :: M       
       End Function Mat_Cross_in    
 
       Module Function Mat_Cross_sp(Vec) Result(M)    
          !---- Argument ----!       
          real(kind=sp), dimension(3), intent( in) :: Vec       
          real(kind=sp), dimension(3,3)            :: M       
       End Function Mat_Cross_sp    
 
       Module Function Rotate_OX(Vec,Angle) Result(Rvec)    
          !---- Arguments ----!       
          real(kind=cp), dimension(3), intent(in) :: Vec      ! Vector       
          real(kind=cp),               intent(in) :: angle    ! Angle       
          real(kind=cp), dimension(3)             :: Rvec     ! Vector rotated       
       End Function Rotate_OX    
 
       Module Function Rotate_OY(Vec,Angle) Result(Rvec)    
          !---- Arguments ----!       
          real(kind=cp), dimension(3), intent(in) :: Vec     ! Vector       
          real(kind=cp),               intent(in) :: angle   ! Angle       
          real(kind=cp), dimension(3)             :: Rvec    ! Vector rotated       
       End Function Rotate_OY    
 
       Module Function Rotate_OZ(Vec,Angle) Result(Rvec)    
          !---- Arguments ----!       
          real(kind=cp), dimension(3), intent(in) :: Vec      ! Vector       
          real(kind=cp),               intent(in) :: angle    ! Angle       
          real(kind=cp), dimension(3)             :: Rvec     ! Vector rotated       
       End Function Rotate_OZ    
 
       Module Function Tensor_Product_cmpl_dp(Vec1,Vec2) Result(w)    
          !---- Argument ----!       
          complex(kind=dp), dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2       
          complex(kind=dp), dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2       
       End Function Tensor_Product_cmpl_dp    
 
       Module Function Tensor_Product_cmpl_sp(Vec1,Vec2) Result(w)    
          !---- Argument ----!       
          complex(kind=sp), dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2       
          complex(kind=sp), dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2       
       End Function Tensor_Product_cmpl_sp    
 
       Module Function Tensor_Product_dp(Vec1,Vec2) Result(w)    
          !---- Argument ----!       
          real(kind=dp), dimension(3), intent( in) :: Vec1,Vec2       
          real(kind=dp), dimension(3,3)            :: w       
       End Function Tensor_Product_dp    
 
       Module Function Tensor_Product_in(Vec1,Vec2) Result(w)    
          !---- Argument ----!       
          integer, dimension(3), intent( in) :: Vec1,Vec2       
          integer, dimension(3,3)            :: w       
       End Function Tensor_Product_in    
 
       Module Function Tensor_Product_sp(Vec1,Vec2) Result(w)    
          !---- Argument ----!       
          real(kind=sp), dimension(3), intent( in) :: Vec1,Vec2       
          real(kind=sp), dimension(3,3)            :: w       
       End Function Tensor_Product_sp    
 
       Module Function Vec_Length(G,Vec) Result(c)    
          !---- Arguments ----!       
          real(kind=cp), intent(in)  , dimension(3,3)       :: G      ! Metric array       
          real(kind=cp), intent(in)  , dimension(3  )       :: Vec    ! Vector       
          real(kind=cp)                                     :: c      ! Length of Vector       
       End Function Vec_Length    
 
       Module Subroutine Set_Eps(Neweps)    
          !---- Arguments ----!       
          real(kind=cp), optional, intent( in) :: neweps  ! Sets global EPS to the value "neweps"       
       End Subroutine Set_Eps    
 
       Module Subroutine Get_Cart_from_Cylin_dp(rho,Phi,z,Xo,Mode)    
          !---- Arguments ----!       
          real(kind=dp),              intent( in)           ::  rho    ! Coordinates rho,phi,zeta       
          real(kind=dp),              intent( in)           ::  phi       
          real(kind=dp),              intent( in)           ::  z       
          real(kind=dp), dimension(3),intent(out)           ::  xo     ! Cartesian coordinates       
          character(len=*),           intent( in), optional ::  mode   ! "D" angles in degrees, otherwise in radians       
       End Subroutine Get_Cart_from_Cylin_dp    
 
       Module Subroutine Get_Cart_from_Cylin_sp(rho,Phi,z,Xo,Mode)    
          real(kind=sp),              intent( in)           ::  rho  ! Coordinates rho,phi,zeta       
          real(kind=sp),              intent( in)           ::  phi       
          real(kind=sp),              intent( in)           ::  z       
          real(kind=sp), dimension(3),intent(out)           ::  xo   ! Cartesian coordinates       
          character(len=*),           intent( in), optional ::  mode ! "D" angles in degrees, otherwise in radians       
       End Subroutine Get_Cart_from_Cylin_sp    
 
       Module Subroutine Get_Cart_from_Spher_dp(r,Theta,Phi,Xo,Mode)    
          !---- Arguments ----!       
          real(kind=dp),              intent( in)           :: r       ! Coordinates (R,Theta;Phi)       
          real(kind=dp),              intent( in)           :: Theta       
          real(kind=dp),              intent( in)           :: phi       
          real(kind=dp), dimension(3),intent(out)           :: xo      ! Cartesian coordinates       
          character(len=*),           intent( in), optional :: mode    ! If "D" the angles are in degrees, otherwise radians is considered       
       End Subroutine Get_Cart_from_Spher_dp    
 
       Module Subroutine Get_Cart_from_Spher_sp(r,Theta,Phi,Xo,Mode)    
          !---- Arguments ----!       
          real(kind=sp),              intent( in)           :: r       ! Coordinates R, Theta, Phi       
          real(kind=sp),              intent( in)           :: Theta       
          real(kind=sp),              intent( in)           :: phi       
          real(kind=sp), dimension(3),intent(out)           :: xo      ! Cartesian Coordinates       
          character(len=*),           intent( in), optional :: mode    ! If "D" then angles are given in degrees.       
       End Subroutine Get_Cart_from_Spher_sp    
 
       Module Subroutine Get_Cylin_from_Cart_dp(Xo,rho,Phi,z,Mode)    
          !---- Arguments ----!       
          real(kind=dp), dimension(3),intent( in)           ::  xo   ! Cartesian coordinatates       
          real(kind=dp),              intent(out)           ::  rho  ! Cylindrical coordinates       
          real(kind=dp),              intent(out)           ::  phi       
          real(kind=dp),              intent(out)           ::  z       
          character(len=*),           intent( in), optional ::  mode       
       End Subroutine Get_Cylin_from_Cart_dp    
 
       Module Subroutine Get_Cylin_from_Cart_sp(Xo,rho,Phi,z,Mode)    
          !---- Arguments ----!       
          real(kind=sp), dimension(3),intent( in)           ::  xo    ! Cartesian       
          real(kind=sp),              intent(out)           ::  rho   ! Cylindrical       
          real(kind=sp),              intent(out)           ::  phi       
          real(kind=sp),              intent(out)           ::  z       
          character(len=*),           intent( in), optional ::  mode       
       End Subroutine Get_Cylin_from_Cart_sp    
 
       Module Subroutine Get_Spher_from_Cart_dp(xo,r,theta,phi,mode)    
          !---- Arguments ----!       
          real(kind=dp), intent( in), dimension(3)   :: xo      ! Cartesian       
          real(kind=dp), intent(out)                 :: r      ! Spherical       
          real(kind=dp), intent(out)                 :: theta       
          real(kind=dp), intent(out)                 :: phi       
          character(len=*), intent(in), optional     :: mode       
       End Subroutine Get_Spher_from_Cart_dp    
 
       Module Subroutine Get_Spher_from_Cart_sp(xo,r,theta,phi,mode)    
          !---- Arguments ----!       
          real(kind=sp), intent( in), dimension(3)   :: xo      ! Cartesian coordinates       
          real(kind=sp), intent(out)                 :: r       ! Spherical       
          real(kind=sp), intent(out)                 :: theta       
          real(kind=sp), intent(out)                 :: phi       
          character(len=*), intent(in), optional     :: mode       
       End Subroutine Get_Spher_from_Cart_sp    
 
       Module Subroutine Matrix_DiagEigen(Array,EigenValues,EigenVec)    
          !---- Arguments ----!       
          real(kind=cp), intent(in)  , dimension(3,3)    :: array       
          real(kind=cp), intent(out) , dimension(3)      :: EigenValues       
          real(kind=cp), intent(out) , dimension(3,3)    :: EigenVec       
       End Subroutine Matrix_DiagEigen    
 
       Module Subroutine Resolv_Sist_1x2(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(2),         intent( in) :: w      ! Input vector       
          real(kind=cp),                 intent( in) :: t      ! Input value       
          real(kind=cp), dimension(2),   intent(out) :: ts     ! Fixed value solution       
          real(kind=cp), dimension(2),   intent(out) :: x      ! Fixed value for x,y       
          integer, dimension(2),         intent(out) :: ix     ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_1x2    
 
       Module Subroutine Resolv_Sist_1x3(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(3),         intent( in) :: w      ! Input vector       
          real(kind=cp),                 intent( in) :: t      ! Input value       
          real(kind=cp), dimension(3),   intent(out) :: ts     ! Fixed value solution       
          real(kind=cp), dimension(3),   intent(out) :: x      ! Fixed value for x,y,z       
          integer, dimension(3),         intent(out) :: ix     ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_1x3    
 
       Module Subroutine Resolv_Sist_2x2(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(2,2),       intent( in) :: w       ! Input vector       
          real(kind=cp),dimension(2),    intent( in) :: t       ! Input value       
          real(kind=cp),dimension(2),    intent(out) :: ts      ! Fixed value solution       
          real(kind=cp),dimension(2),    intent(out) :: x       ! Fixed value for x,y       
          integer, dimension(2),         intent(out) :: ix      ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_2x2    
 
       Module Subroutine Resolv_Sist_2x3(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(2,3),          intent( in) :: w     ! Input vector       
          real(kind=cp), dimension(2),      intent( in) :: t     ! Input value       
          real(kind=cp), dimension(3),      intent(out) :: ts    ! Fixed value solution       
          real(kind=cp), dimension(3),      intent(out) :: x     ! Fixed value for x,y       
          integer, dimension(3),            intent(out) :: ix    ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_2x3    
 
       Module Subroutine Resolv_Sist_3x3(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(3,3),          intent(in) :: w       ! Input vector       
          real(kind=cp), dimension(3),      intent(in) :: t       ! Input value       
          real(kind=cp), dimension(3),      intent(out):: ts      ! Fixed value solution       
          real(kind=cp), dimension(3),      intent(out):: x       ! Fixed value for x,y       
          integer, dimension(3),            intent(out):: ix      ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_3x3    
 
    End Interface
    
 Contains
    
    !!---- SUBROUTINE SET_EPS
    !!----
    !!----    Sets global EPS to the value "neweps"
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine Set_Eps(Neweps)    
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
