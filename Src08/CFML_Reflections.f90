!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2019  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_Reflections_Utilities
!!----   INFO: Series of procedures handling operation with
!!----         Bragg reflections
!!----
!!
Module CFML_Reflections
   !---- Use Modules ----!
   Use CFML_GlobalDeps,                only: CP, PI, TPI, Err_CFML, Clear_Error
   Use CFML_gSpaceGroups,              only: Spg_Type
   Use CFML_Maths,                     only: Trace, Sort, Equal_vector
   Use CFML_Metrics,                   only: Cell_G_Type
   Use CFML_Strings,                   only: l_case
   Use CFML_Rational

   !---- Variables ----!
   implicit none

   private

   !---- List of public functions ----!
   public :: H_Absent, mH_Absent , H_Equal, H_Latt_Absent , H_Equiv, H_Mult, H_S, &
             Get_MaxNumRef, Get_Asymm_Unit_H

   !---- List of public subroutines ----!
   public :: H_Equiv_List, Initialize_RefList, Gener_Reflections, Search_Extinctions, &
             Write_Info_RefList

   !---- Parameters ----!
   real(kind=cp), parameter :: EPS_REF  = 0.0002_cp   ! Epsilon for comparisons within this module

   !---- Types ----!

   !!----
   !!---- TYPE :: REFL_TYPE
   !!--..
   !!
   Type, public :: Refl_Type
      integer,dimension(:), allocatable :: H            ! H
      integer                           :: Mult  = 0     ! Mutiplicity
      real(kind=cp)                     :: S     = 0.0   ! Sin(Theta)/lambda=1/2d
      integer                           :: Imag  = 0     ! 0: nuclear reflection, 1:magnetic, 2=both
      integer                           :: Pcoeff= 0     ! Pointer to the harmonic q_coeff
   End Type Refl_Type

   !!----
   !!---- TYPE :: SREFL_TYPE
   !!--..
   !!
   Type, public, extends(Refl_Type) :: SRefl_Type
      real(kind=cp)        :: Fo    =0.0_cp  ! Observed Structure Factor
      real(kind=cp)        :: Fc    =0.0_cp  ! Calculated Structure Factor
      real(kind=cp)        :: SFo   =0.0_cp  ! Sigma of  Fo
      real(kind=cp)        :: Phase =0.0_cp  ! Phase in degrees
      real(kind=cp)        :: A     =0.0_cp  ! real part of the Structure Factor
      real(kind=cp)        :: B     =0.0_cp  ! Imaginary part of the Structure Factor
      real(kind=cp)        :: W     =1.0_cp  ! Weight factor
   End Type SRefl_Type

   !!----
   !!---- TYPE :: MREFL_TYPE
   !!--..
   !!
   Type, Public, extends(SRefl_Type) :: MRefl_Type
      real(kind=cp)                    :: mIvo  =0.0_cp                 ! Observed modulus of the Magnetic Interaction vector
      real(kind=cp)                    :: smIvo =0.0_cp                 ! Sigma of observed modulus of the Magnetic Interaction vector
      complex(kind=cp),dimension(3)    :: msF   =cmplx(0.0_cp,0.0_cp)   ! Magnetic Structure Factor
      complex(kind=cp),dimension(3)    :: mIv   =cmplx(0.0_cp,0.0_cp)   ! Magnetic Interaction vector
   End Type MRefl_Type

   !!----
   !!---- TYPE :: REFLIST_TYPE
   !!--..
   !!
   Type, public :: RefList_Type
      integer                                     :: NRef=0 ! Number of Reflections
      class(refl_type), dimension(:), allocatable :: Ref    ! Reflection List
   End Type RefList_Type

   !!----
   !!---- TYPE :: KVECT_INFO_TYPE
   !!--..
   !!
   Type, public :: Kvect_Info_Type
      integer                                      :: nk=0        ! Number of independent k-vectors
      real(kind=cp),allocatable,dimension(:,:)     :: kv          ! k-vectors (3,nk)
      real(kind=cp),allocatable,dimension(:)       :: sintlim     ! sintheta/lambda limits (nk)
      integer,allocatable,dimension(:)             :: nharm       ! number of harmonics along each k-vector
      integer                                      :: nq=0        ! number of effective set of Q_coeff > nk
      integer,allocatable,dimension(:,:)           :: q_coeff     ! number of q_coeff(nk,nq)
   End Type kvect_info_type


   !---- Private Variables ----!
   logical :: hkl_ref_cond_ini=.false.                ! reflection conditions array has been initialized

   !---- Public Variables ----!
   character(len=80), dimension(58),  public :: Hkl_Ref_Conditions=" "  ! Reflection conditions



   !---- Overload  Zone ----!
   Interface Search_Extinctions
       Module Procedure Search_Extinctions_Iunit
       Module Procedure Search_Extinctions_File
   End Interface Search_Extinctions

   Interface H_S
       Module Procedure H_S_int
       Module Procedure H_S_real
   End Interface H_S


   !---- Interface Zone ----!
   Interface
      Module Function H_Equal(H,K) Result (Info)
         !---- Arguments ----!
         integer, dimension(:), intent(in) :: H
         integer, dimension(:), intent(in) :: K
         logical                           :: info
      End Function H_Equal

      Module Function H_Absent(H, SpG) Result(Info)
         !---- Arguments ----!
         integer, dimension(:), intent (in) :: H
         class(SpG_Type),       intent (in) :: SpG
         logical                            :: info
      End Function H_Absent

      Module Function mH_Absent(H,SpG) Result(Info)
         !---- Arguments ----!
         integer, dimension(:), intent (in) :: H
         class(SpG_Type),       intent (in) :: SpG
         logical                            :: info
      End Function mH_Absent

      Module Function H_Latt_Absent(H, Latt, n) Result(info)
         !---- Arguments ----!
         integer,        dimension(:),  intent (in) :: h
         type(rational), dimension(:,:),intent (in) :: Latt
         integer,                       intent (in) :: n
         logical                                    :: info
      End Function H_Latt_Absent

      Module Function H_Equiv(H, K, SpG, Friedel) Result (Info)
         !---- Arguments ----!
         integer, dimension(:),        intent(in)  :: H
         integer, dimension(:),        intent(in)  :: K
         class(Spg_Type),              intent(in)  :: SpG
         logical, optional,            intent(in)  :: Friedel
         logical                                   :: info
      End Function H_Equiv

      Module Subroutine H_Equiv_List(H, SpG, Friedel, Mult, H_List)
         !---- Arguments ----!
         integer, dimension(:),    intent (in) :: H
         class(SpG_Type),          intent (in) :: SpG
         Logical,                  intent (in) :: Friedel
         integer,                  intent(out) :: mult
         integer, dimension(:,:),  intent(out) :: h_list
      End Subroutine H_Equiv_List

      Module Function H_Mult(H, SpG, Friedel) Result(N)
         !---- Arguments ----!
         integer, dimension(:),  intent (in)  :: H
         class(SpG_Type),        intent (in)  :: SpG
         Logical,                intent (in)  :: Friedel
         integer                              :: N
      End Function H_Mult

      Module Function H_S_real(H, Cell)  Result(S)
         !---- Arguments ----!
         real(kind=cp), dimension(3),            intent(in) :: h
         class(Cell_G_Type),                     intent(in) :: Cell
         real(kind=cp)                                      :: S
      End Function H_S_real

      Module Function H_S_int(H, Cell, Nk, Kv) Result(S)
         !---- Arguments ----!
         integer, dimension(:),                  intent(in) :: h
         class(Cell_G_Type),                     intent(in) :: Cell
         integer, optional,                      intent(in) :: Nk
         real(kind=cp),dimension(:,:), optional, intent(in) :: Kv
         real(kind=cp)                                      :: S
      End Function H_S_int

      Module Function Get_MaxNumRef(SinTLMax, VolCell, SinTLMin, Mult) Result(numref)
         !---- Arguments ----!
         real(kind=cp),           intent(in) :: SinTLMax
         real(kind=cp),           intent(in) :: VolCell
         real(kind=cp), optional, intent(in) :: SinTLMin
         integer,       optional, intent(in) :: Mult
         integer                             :: numref
      End Function Get_MaxNumRef

      Module Function Unitary_Vector_H(H, Cell) Result (U)
         !---- Arguments ----!
         integer, dimension(3), intent(in) :: H
         class(Cell_G_Type),    intent(in) :: Cell
         real(kind=cp), dimension(3)       :: U
      End Function Unitary_Vector_H

      Module Function Asu_H(H, SpG) Result(K)
         !---- Arguments ----!
         integer, dimension (3),  intent(in) :: h
         class(SpG_Type),         intent(in) :: SpG
         integer, dimension(3)               :: k
      End Function Asu_H

      Module Function Asu_H_Cubic(H, Laue) Result(K)
         !---- Argument ----!
         integer, dimension(3), intent(in) :: H
         character(len=*),      intent(in) :: Laue
         integer, dimension(3)             :: K
      End Function Asu_H_Cubic

      Module Function Asu_H_Hexagonal(H, Laue) Result(K)
         !---- Argument ----!
         integer, dimension(3), intent(in) :: h
         character(len=*),      intent(in) :: Laue
         integer, dimension(3)             :: k
      End Function Asu_H_Hexagonal

      Module Function Asu_H_Monoclinic(H, Axis) Result(K)
         !---- Argument ----!
         integer, dimension(3),      intent(in) :: h
         character(len=*), optional, intent(in) :: Axis
         integer, dimension(3)                  :: k
      End Function Asu_H_Monoclinic

      Module Function Asu_H_Orthorhombic(H) Result(K)
         !---- Argument ----!
         integer, dimension(3), intent(in) :: h
         integer, dimension(3)             :: k
      End Function Asu_H_Orthorhombic

      Module Function Asu_H_Tetragonal(H,Laue) Result(K)
         !---- Argument ----!
         integer, dimension(3), intent(in) :: h
         character(len=*),      intent(in) :: Laue
         integer, dimension(3)             :: k
      End Function Asu_H_Tetragonal

      Module Function Asu_H_Triclinic(H) Result(K)
         !---- Argument ----!
         integer, dimension(3), intent(in) :: h
         integer, dimension(3)             :: k
      End Function Asu_H_Triclinic

      Module Function Asu_H_Trigonal(H, Laue) Result(K)
         !---- Argument ----!
         integer, dimension(3), intent(in) :: h
         character(len=*),      intent(in) :: Laue
         integer, dimension(3)             :: k
      End Function Asu_H_Trigonal

      Module Function Get_Asymm_Unit_H(H,SpG) Result(k)
         !---- Arguments ----!
         integer, dimension (3),  intent(in) :: h
         class(SpG_Type),         intent(in) :: SpG
         integer, dimension(3)               :: k
      End Function Get_Asymm_Unit_H

      Module Subroutine Gener_Reflections(Cell,Sintlmax,Mag,Num_Ref,Reflex,SpG,kinfo,order,powder,mag_only)
         !---- Arguments ----!
         class(Cell_G_Type),                          intent(in)     :: Cell
         real(kind=cp),                               intent(in)     :: Sintlmax
         logical,                                     intent(in)     :: Mag
         integer,                                     intent(out)    :: Num_ref
         class(Refl_Type), dimension(:), allocatable, intent(out)    :: Reflex
         class(Spg_Type) ,              optional,     intent(in)     :: SpG
         type(kvect_info_type),         optional,     intent(in)     :: Kinfo
         character(len=*),              optional,     intent(in)     :: Order
         logical,                       optional,     intent(in)     :: Powder
         logical,                       optional,     intent(in)     :: Mag_only
      End Subroutine Gener_Reflections

      Module Subroutine Init_Refl_Conditions()
         !---- Arguments ----!
      End Subroutine Init_Refl_Conditions

      Module Subroutine Write_Integral_Conditions(SpG,iunit)
         !---- Arguments ----!
         class(SpG_Type),    intent(in)  :: SpG
         integer, optional,  intent(in)  :: iunit
      End Subroutine Write_Integral_Conditions

      Module Subroutine Write_Glide_Planes_Conditions(SpG,Iunit)
         !---- Arguments ----!
         class(SpG_Type),    intent(in) :: SpG
         integer, optional,  intent(in) :: Iunit
      End Subroutine Write_Glide_Planes_Conditions

      Module Subroutine Write_Screw_Axis_Conditions(SpG ,Iunit)
         !---- Arguments ----!
         class(SpG_Type),    intent(in) :: SpG
         integer, optional,  intent(in) :: Iunit
      End Subroutine Write_Screw_Axis_Conditions

      Module Subroutine Search_Extinctions_Iunit(SpG, Iunit)
         !---- Arguments ----!
         class(SpG_Type),   intent(in) :: SpG
         integer, optional, intent(in) :: Iunit
      End Subroutine Search_Extinctions_Iunit

      Module Subroutine Search_Extinctions_File(SpG, nlines, filevar)
         !---- Arguments ----!
         class(SpG_Type),                intent(in)   :: SpG
         integer,                        intent(out)  :: nlines
         character(len=*), dimension(:), intent(out)  :: filevar
      End Subroutine Search_Extinctions_File

      Module Subroutine Initialize_RefList(N, Reflex)
         !---- Arguments ----!
         integer,             intent(in)    :: N
         type(RefList_Type),  intent(in out) :: Reflex
      End Subroutine Initialize_RefList

      Module Subroutine Write_Info_RefList(Reflex, Iunit, Mode)
         !---- Arguments ----!
         type(RefList_Type),         intent(in) :: Reflex
         integer,          optional, intent(in) :: Iunit
         character(len=*), optional, intent(in) :: Mode
      End Subroutine Write_Info_RefList

   End Interface



End Module CFML_Reflections

