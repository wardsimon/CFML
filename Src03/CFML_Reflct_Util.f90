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
!!---- MODULE: CFML_Reflections_Utilities
!!----   INFO: Series of procedures handling operation with
!!----         Bragg reflections
!!----
!!---- HISTORY
!!----    Update: 06/03/2011
!!----
!!----
!!----
!!
 Module CFML_Reflections_Utilities

    !---- Use Modules ----!
    Use CFML_GlobalDeps,                only: sp, cp, pi
    Use CFML_Math_General,              only: sort
    Use CFML_String_Utilities,          only: l_case
    Use CFML_Crystallographic_Symmetry, only: Sym_Oper_Type, Space_Group_Type, Lattice_Centring_Type, &
                                              Allocate_Lattice_Centring
    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type

    !---- Definitions ----!
    implicit none

    private

    !---- Public procedures ----!
    public :: Asu_Hkl,Get_MaxNumRef, Hkl_Absent, Hkl_Equal, Hkl_Equiv, Hkl_Mult,     &
              Get_Hequiv_Asu,Hkl_R, Hkl_S, Unit_Cart_Hkl, Hkl_Lat_Absent

    public :: Hkl_Equiv_List, Hkl_Gen, Hkl_Rp, Hkl_Uni, Init_Err_Refl, Init_RefList, &
              Search_Extinctions, Write_Asu, Write_RefList_Info, Hkl_Gen_Sxtal

    !--------------------!
    !---- PARAMETERS ----!
    !--------------------!
    real(kind=cp), parameter :: EPS_REF  = 0.0002_cp  ! Epsilon for comparisons within this module

    !---------------!
    !---- TYPES ----!
    !---------------!

    !!----
    !!---- TYPE :: REFLECT_TYPE
    !!--..
    !!----
    !!---- Update: 20/07/2015
    !!
    Type, public :: Reflect_Type
       integer,dimension(3) :: H     ! H
       integer              :: Mult  ! mutiplicity
       real(kind=cp)        :: S     ! Sin(Theta)/lambda=1/2d
    End Type Reflect_Type

    !!----
    !!---- TYPE :: REFLECTION_TYPE
    !!--..
    !!---- Update: 20/07/2015
    !!
    Type, public :: Reflection_Type
       integer,dimension(3) :: H     ! H
       integer              :: Mult  ! mutiplicity
       real(kind=cp)        :: Fo    ! Observed Structure Factor
       real(kind=cp)        :: Fc    ! Calculated Structure Factor
       real(kind=cp)        :: SFo   ! Sigma of  Fo
       real(kind=cp)        :: S     ! Sin(Theta)/lambda
       real(kind=cp)        :: W     ! Weight
       real(kind=cp)        :: Phase ! Phase in degrees
       real(kind=cp)        :: A     ! real part of the Structure Factor
       real(kind=cp)        :: B     ! Imaginary part of the Structure Factor
       real(kind=cp)        :: AA    ! Free
       real(kind=cp)        :: BB    ! Free
    End Type Reflection_Type

    !!----
    !!---- TYPE :: REFLECTION_LIST_TYPE
    !!--..
    !!---- Update: 20/07/2015
    !!
    Type, public :: Reflection_List_Type
       integer                                         :: NRef  ! Number of Reflections
       type(reflection_type),allocatable, dimension(:) :: Ref   ! Reflection List
    End Type Reflection_List_Type

    !-------------------!
    !---- VARIABLES ----!
    !-------------------!
    logical                                   :: HKL_Ref_Cond_Ini=.false. ! reflection conditions array has been initialized

    logical,                           public :: ERR_Refl=.false.         ! Logical Variable indicating an error in CFML_Reflections_Utilities module
    character(len=256),                public :: ERR_Refl_Mess=" "        ! String containing information about the last error
    character(len=80), dimension(58),  public :: Hkl_Ref_Conditions       ! Reflection conditions for Lattices, glide planes, screw axes
                                                                          ! International tables vol. A, Fourth, revised edition (1996) p27-28-29

    !-----------------------------!
    !---- Interfaces/Overload ----!
    !-----------------------------!

    Interface Hkl_Absent
       Module Procedure hkl_AbsentI
       Module Procedure hkl_AbsentR
    End Interface Hkl_Absent

    Interface Hkl_Equal
       Module Procedure Hkl_EqualI
       Module Procedure Hkl_EqualR
    End Interface Hkl_Equal

    Interface Hkl_Equiv
       Module Procedure Hkl_EquivI
       Module Procedure Hkl_EquivR
    End Interface Hkl_Equiv

    Interface Hkl_Mult
       Module Procedure Hkl_MultI
       Module Procedure Hkl_MultR
    End Interface Hkl_Mult

    Interface Hkl_R
       Module Procedure HR_I
       Module Procedure HR_R
    End Interface Hkl_R

    Interface Hkl_S
       Module Procedure HS_I
       Module Procedure HS_R
    End Interface Hkl_S

    Interface Hkl_Equiv_List
       Module Procedure Hkl_Equiv_ListI
       Module Procedure Hkl_Equiv_ListR
    End Interface Hkl_Equiv_List

    Interface Hkl_Rp
       Module Procedure Hkl_RpI
       Module Procedure Hkl_RpR
    End Interface Hkl_Rp

    Interface HKL_GEN_SXTAL
       Module Procedure HKL_GEN_SXTAL_reflection
       Module Procedure HKL_GEN_SXTAL_list
    End Interface HKL_GEN_SXTAL

    Interface Hkl_Uni
       Module Procedure Hkl_Uni_reflect
       Module Procedure Hkl_Uni_reflection
       Module Procedure Hkl_Uni_ReflList
    End Interface Hkl_Uni

    Interface Search_Extinctions
       Module Procedure Search_Extinctions_Iunit
       Module Procedure Search_Extinctions_File
    End Interface Search_Extinctions

    Interface Unit_Cart_Hkl
       Module Procedure Unit_Cart_HklI
       Module Procedure Unit_Cart_HklR
    End Interface Unit_Cart_Hkl

 Contains

    !-------------------!
    !---- FUNCTIONS ----!
    !-------------------!

    !!----
    !!---- FUNCTION ASU_HKL
    !!----
    !!----    Obtain an equivalent reflection in asymmetric unit using
    !!----    simple transformation rules for each crystal system.
    !!----    When these rules are not satisfied the output is the
    !!----    (0,0,0) reflection. For obtaining a reflection within
    !!----    the asymmetric unit given an input reflection the best
    !!----    is to use the function: Get_Hequiv_Asu
    !!----
    !!--<<
    !!----    We assumed that F(hkl)=F(-h -k -l).
    !!-->>
    !!----    If and error occurs, the function returns also (0,0,0).
    !!----
    !!---- Update: 20/07/2015
    !!
    Function Asu_Hkl(H,SpG) Result(K)
       !---- Arguments ----!
       integer, dimension (3),  intent(in) :: h    ! Input reflection
       type (Space_Group_Type), intent(in) :: SpG  ! Spacegroup
       integer, dimension(3)               :: k

       !---- Local Variables ----!
       character(len=2)  :: inf

       !> Init
       k=0

       if (SpG%NumSpg > 0 .and. SpG%NumSpg <= 231) then
          select case (SpG%NumSpg)
             case (1:2)
                k=asu_hkl_triclinic(h)

             case (3:15)
                inf(1:2)=adjustl(SpG%info(1:2))
                if(inf(1:1) == "-") inf(1:1)=inf(2:2)
                select case (inf(1:1))
                   case ("b")
                      k=asu_hkl_monoclinic(h,"b")
                   case ("c")
                      k=asu_hkl_monoclinic(h,"c")
                   case ("a")
                      k=asu_hkl_monoclinic(h,"a")
                   case default
                      k=asu_hkl_monoclinic(h,"b")
                end select

             case (16:74)
                k=asu_hkl_orthorhombic(h)

             case (75:88)
                k=asu_hkl_tetragonal(h,"4/m  ")

             case (89:142)
                k=asu_hkl_tetragonal(h,"4/mmm")

             case (143:148)
                k=asu_hkl_trigonal(h,"-3  ")

             case (149,151,153,157,159,162,163)
                k=asu_hkl_trigonal(h,"-31m")

             case (150,152,154,155,156,158,160,161,164,165,166,167)
                k=asu_hkl_trigonal(h,"-3m")

             case (168:176)
                k=asu_hkl_hexagonal(h,"6/m  ")

             case (177:194)
                k=asu_hkl_hexagonal(h,"6/mmm")

             case (195:206)
                k=asu_hkl_cubic(h,"m-3 ")

             case (207:230)
                k=asu_hkl_cubic(h,"m-3m")

          end select

       else

          !---- General ----!
          select case(trim(SpG%Laue))
             case("-1")
                k=asu_hkl_triclinic(h)
             case("2/m")
                k=asu_hkl_monoclinic(h,"b")
             case("mmm")
                k=asu_hkl_orthorhombic(h)
             case("4/m")
                k=asu_hkl_tetragonal(h,"4/m  ")
             case("4/mmm")
                k=asu_hkl_tetragonal(h,"4/mmm")
             case("-3")
                k=asu_hkl_trigonal(h,"-3  ")
             case("-3m")
                k=asu_hkl_trigonal(h,"-3m")
             case("6/m")
                k=asu_hkl_hexagonal(h,"6/m  ")
             case("6/mmm")
                k=asu_hkl_hexagonal(h,"6/mmm")
             case("m-3")
                k=asu_hkl_cubic(h,"m-3 ")
             case("m-3m")
                k=asu_hkl_cubic(h,"m-3m")
             case default
               return
          end select

       end if

       return
    End Function Asu_Hkl

    !!--++
    !!--++ FUNCTION ASU_HKL_CUBIC
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtain a reflection in asymmetric unit for Cubic
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Asu_Hkl_Cubic(H,Mode) Result(K)
       !---- Argument ----!
       integer, dimension(3), intent(in) :: h     ! Reflection
       character(len=*),      intent(in) :: Mode  ! Laue information
       integer, dimension(3)             :: k

       !---- Local Variable ----!
       character(len=4)      :: mod_laue
       integer, dimension(3) :: hh

       !> Init
       k=0

       mod_laue=l_case(adjustl(Mode))
       if (len_trim(mod_laue) == 0) return

       select case(trim(mod_laue))
          case("m-3")
             !---- Laue: m-3 ----!
             !---- hkl: h>l, k>l, l>=0 ; hkk: k>=0 h>=k ----!
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select
             if (hh(3) >=0 .and. hh(1) >= hh(3) .and. hh(2) == hh(3)) k=hh
             if (hh(3) >=0 .and. hh(1) >  hh(3) .and. hh(2) >  hh(3)) k=hh

          case("m-3m")
             !---- Laue: m-3m ----!
             !---- hkl: h >=0, k >=0, l >=0, h >=k, k >=l ----!
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select
             if (hh(3) >= 0 .and. hh(2) >= hh(3) .and. hh(1) >= hh(2)) k=hh

          case default
             return
       end select

       return
    End Function Asu_Hkl_Cubic

    !!--++
    !!--++ FUNCTION ASU_HKL_HEXAGONAL
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtain a reflection in asymmetric unit for Hexagonal
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Asu_Hkl_Hexagonal(H,Mode) Result(K)
       !---- Argument ----!
       integer, dimension(3), intent(in) :: h      ! Input reflection
       character(len=*),      intent(in) :: Mode   ! Laue
       integer, dimension(3)             :: k

       !---- Local Variable ----!
       character(len=5)      :: mod_laue
       integer, dimension(3) :: hh

       k=0
       mod_laue=l_case(adjustl(Mode))
       if (len_trim(mod_laue) == 0) return

       select case(trim(mod_laue))
          case("6/m")
             !---- Laue: 6/m ----!
             !---- hkl: h>0,k>0,l>=0;  0kl k>=0,l>=0 ----!
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select
             if (hh(1) > 0 .and. hh(2) > 0 .and. hh(3) >= 0) k=hh
             if (hh(1) == 0 .and. hh(2) >= 0 .and. hh(3) >= 0) k=hh

          case("6/mmm")
             !---- Laue: 6/mmm ----!
             !---- hkl: h >=0, k >=0, l >=0, h >=k ----!
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select
             if (hh(2) >=0 .and. hh(1) >= hh(2) .and. hh(3) >= 0) k=hh

          case default
             return
       end select

       return
    End Function Asu_Hkl_Hexagonal

    !!--++
    !!--++ FUNCTION ASU_HKL_MONOCLINIC
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtain a reflection in asymmetric unit for Monoclinic
    !!--++    Unique axis b: hkl: k >=0, l >=0    hk0: h >=0
    !!--++    Unique axis c: hkl: k >=0, l >=0    h0l: h >=0
    !!--++    Unique axis a: hkl: h >=0, l >=0    0kl: l >=0
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Asu_Hkl_Monoclinic(H,Mode) Result(K)
       !---- Argument ----!
       integer, dimension(3),      intent(in) :: h     ! Input reflection
       character(len=*), optional, intent(in) :: mode  ! Axis
       integer, dimension(3)                  :: k

       !---- Local Variable ----!
       character(len=1)     :: ax
       integer,dimension(3) :: hh

       !> Init
       k=0

       if (present(mode)) then
          ax=l_case(adjustl(mode))
          if (ax ==" ") ax="b"
       else
          ax="b"
       end if

       select case (ax)
          !---- Laue: 2/m     Unique Axis: b ----!
          !---- hkl: k >=0, l >=0    hk0: h >=0
          case ("b")
             select case (h(3))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(1) >=0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select

             if (hh(3) == 0) then
                if (hh(1) >=0 ) k=hh
             else
                if (hh(2) >=0 .and. hh(3) >=0) k=hh
             end if

          !---- Laue: 2/m     Unique Axis: c ----!
          !---- hkl: k >=0, l >=0    h0l: h >=0
          case ("c")
             select case (h(3))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(1) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select

             if (hh(2) == 0) then
                if (hh(1) >= 0) k=hh
             else
                if (hh(2) >=0 .and. hh(3) >=0) k=hh
             end if

          !---- Laue: 2/m     Unique Axis: c ----!
          !---- hkl: h >=0, l >=0    0kl: l >=0
          case ("a")
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(3))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(2) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select

             if (hh(1) == 0) then
                if (hh(2) >= 0) k=hh
             else
                if (hh(1) >=0 .and. hh(3) >=0) k=hh
             end if

       end select

       return
    End Function Asu_Hkl_Monoclinic

    !!--++
    !!--++ FUNCTION ASU_HKL_ORTHORHOMBIC
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtain a reflection in asymmetric unit for Orthorhombic
    !!--++    hkl: h >=0, k >=0, l >=0
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Asu_Hkl_Orthorhombic(H) Result(K)
       !---- Argument ----!
       integer, dimension(3), intent(in) :: h  ! Input reflection
       integer, dimension(3)             :: k

       !---- Local Variable ----!
       integer, dimension(3) :: hh

       !> Init
       k=0

       !---- Laue: mmm ----!
       !---- hkl: h >=0, k >=0, l >=0 ----!
       select case (h(1))
          case (:-1)
             hh=-h
          case (0)
             select case (h(2))
                case (:-1)
                   hh=-h
                case (0)
                   if (h(3) >= 0) then
                      hh=h
                   else
                      hh=-h
                   end if
                case (1:)
                   hh=h
             end select
          case (1:)
             hh=h
       end select

       if (hh(1) >= 0 .and. hh(2) >= 0 .and. hh(3) >= 0) k=hh

       return
    End Function Asu_Hkl_Orthorhombic

    !!--++
    !!--++ FUNCTION ASU_HKL_TETRAGONAL
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtain a reflection in asymmetric unit for Tetragonal
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Asu_Hkl_Tetragonal(H,Mode) Result(K)
       !---- Argument ----!
       integer, dimension(3), intent(in) :: h     ! Input reflection
       character(len=*),      intent(in) :: Mode  ! Laue
       integer, dimension(3)             :: k

       !---- Local Variable ----!
       character(len=5)     :: mod_laue
       integer,dimension(3) :: hh

       k=0
       mod_laue=l_case(adjustl(Mode))
       if (len_trim(mod_laue) == 0) return

       select case(trim(mod_laue))
          case("4/m")
             !---- Laue: 4/m ----!
             !---- hkl: h >=0, l >=0, k >=0 if h = 0 ----!
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select
             if (hh(1) == 0 .and. hh(2) >= 0 .and. hh(3) >=0) k=hh
             if (hh(1)  > 0 .and. hh(2) >  0 .and. hh(3) >=0) k=hh

          case("4/mmm")
             !---- Laue: 4/mmm ----!
             !---- hkl: h >=0, l >=0, h >=k   (k >=0) ----!
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select
             if (hh(1) >=0 .and. hh(2) >=0 .and. hh(3) >=0 .and. hh(1) >= hh(2)) k=hh

          case default
             return
       end select

       return
    End Function Asu_Hkl_Tetragonal

    !!--++
    !!--++ FUNCTION ASU_HKL_TRICLINIC
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtain a reflection in asymmetric unit for Triclinic
    !!--++    hkl: l >=0    hk0: h >=0    0k0: k >=0
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Asu_Hkl_Triclinic(H) Result(K)
       !---- Argument ----!
       integer, dimension(3), intent(in) :: h  ! Input reflection
       integer, dimension(3)             :: k

       !> Init
       k=0

       !---- Laue: -1 ----!
       !---- hkl: l >=0    hk0: h >=0    0k0: k >=0
       select case (h(3))
          case (:-1)
             k=-h
          case (0)
             select case (h(1))
                case (:-1)
                   k=-h
                case (0)
                   if (h(2) < 0) then
                      k=-h
                   else
                      k=h
                   end if
                case (1:)
                   k=h
             end select
          case (1:)
             k=h
       end select

       return
    End Function Asu_Hkl_Triclinic

    !!--++
    !!--++ FUNCTION ASU_HKL_TRIGONAL
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtain a reflection in asymmetric unit for Trigonal
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Asu_Hkl_Trigonal(H,Mode) Result(K)
       !---- Argument ----!
       integer, dimension(3), intent(in) :: h     ! Input reflection
       character(len=*),      intent(in) :: Mode  ! Laue
       integer, dimension(3)             :: k

       !---- Local Variable ----!
       character(len=4)      :: mod_laue
       integer, dimension(3) :: hh

       k=0
       mod_laue=l_case(adjustl(Mode))
       if (len_trim(mod_laue) == 0) then
          return
       end if

       select case(trim(mod_laue))
          case("-3")
             !---- Laue: -3 ----!
             !---- hkl: h+k>0, l>0 ; hk0:h>0, k>=0
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h

                   end select
                case (1:)
                   hh=h
             end select
             if (hh(1) == 0 .and. hh(2) == 0 .and. hh(3) > 0) k=hh
             if (hh(1)+hh(2) > 0 .and. hh(3) > 0 ) k=hh
             if (hh(1) > 0  .and. hh(2) >= 0  .and. hh(3) == 0) k=hh

          case("-3m")
             !---- Laue: -3m ----!
             !---- hkl: h>=0, h>=k ; hhl: h>=0,l>=0 ----!
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                         hh=h
                   end select
                case (1:)
                   hh=h
             end select
             if (hh(1) >= hh(2) .and.  hh(2) >= 0 ) k=hh
             if (hh(1) >= 0 .and. hh(2) > 0 .and. hh(3) > 0 ) k=hh
             if (hh(1) >= 0 .and. hh(2) == hh(1) .and. hh(3) >=0) k=hh

          case("-31m")
             !---- Laue: -31m ----!
             !---- hkl: h>=0,h>=k>0 ; h0l: h>=0,l>=0 ----!
             select case (h(1))
                case (:-1)
                   hh=-h
                case (0)
                   select case (h(2))
                      case (:-1)
                         hh=-h
                      case (0)
                         if (h(3) >= 0) then
                            hh=h
                         else
                            hh=-h
                         end if
                      case (1:)
                           hh=h
                   end select
                case (1:)
                   hh=h
             end select
             if (hh(1) >= hh(2) .and. hh(1) >=0 .and. hh(2) > 0) k=hh
             if (hh(1) >= 0 .and. hh(2) ==0 .and. hh(3) >= 0) k=hh

          case default
             return
       end select

       return
    End Function Asu_Hkl_Trigonal

    !!----
    !!---- FUNCTION GET_HEQUIV_ASU
    !!----
    !!----    Provides a reflection equivalent to the input one but
    !!----    within the asymmetric unit
    !!----
    !!---- Update: 20/07/2015
    !!
    Function Get_Hequiv_Asu(H,SpG) Result(k)
       !---- Arguments ----!
       integer, dimension (3),  intent(in) :: h    ! Input reflection
       type (Space_Group_Type), intent(in) :: SpG  ! Spacegroup
       integer, dimension(3)               :: k

       !---- Local Variables ----!
       integer                             :: i
       integer, dimension(3)               :: kk,nul

       !> Init
       k=h

       nul=(/0,0,0/)
       do i=1,SpG%NumOps
          k=matmul(h,SpG%SymOp(i)%RotI)
          kk=asu_hkl(k,SpG)
          if (hkl_equal(kk,nul)) cycle
          k=kk
          exit
       end do

       return
    End Function Get_Hequiv_Asu

    !!----
    !!---- FUNCTION GET_MAXNUMREF
    !!----
    !!----    Provides un upper limit of the expected maximum number of
    !!----    reflections up to SinTLMax for a volume VolCell of the
    !!----    primitive cell. If the optional argument SinTLMin is given,
    !!----    the result is the number of reflections in the interval (SinTLMin,SinTLMax).
    !!----
    !!----    If Mult is provided the result is divided by half this multiplicity
    !!----    so we obtain an estimation of the expected mumber of unique reflections.
    !!----
    !!---- Update: 20/07/2015
    !!
    Function Get_MaxNumRef(SinTLMax, VolCell, SinTLMin, Mult) Result(numref)
       !---- Arguments ----!
       real(kind=cp),           intent(in) :: SinTLMax  ! Maximum sinTheta/Lambda
       real(kind=cp),           intent(in) :: VolCell   ! Direct Cell Volume
       real(kind=cp), optional, intent(in) :: SinTLMin  ! Minimum sinTheta/Lambda
       integer,       optional, intent(in) :: Mult      ! General Multiplicity
       integer                             :: numref

       !---- Local Variables ----!
       real(kind=cp) :: r3

       r3=8.0*SinTLMax*SinTLMax*SinTLMax*1.05
       if (present(SinTLMin)) r3= r3-8.0*SinTLMin*SinTLMin*SinTLMin

       numref=4.0*pi*r3*VolCell/3.0

       !> The factor 2 is given because, for high symmetry, sometimes the obtained number is
       !> not enough for allocating the real number of reflections
       if (present(Mult)) numref=2*numref/max(1,Mult)

       return
    End Function Get_MaxNumRef

    !!----
    !!---- FUNCTION HKL_ABSEN
    !!----
    !!----    Returns the value ".true." if the reflection is absent.
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ LOGICAL FUNCTION  HKL_ABSENTI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if the reflection is an absence
    !!--++
    !!--++  Update: 20/07/2015
    !!
    Function Hkl_AbsentI(H,SpG) Result(Info)
       !---- Arguments ----!
       integer, dimension(3),   intent (in) :: h     ! Input reflection
       Type (Space_Group_Type), intent (in) :: SpG   ! Spacegroup
       logical                              :: info

       !---- Local Variables ----!
       integer, dimension(3)  :: k
       integer                :: i
       real(kind=cp)          :: r1,r2

       info=.false.

       do i=1,SpG%multip
          k = hkl_r(h,SpG%SymOp(i))
          if (hkl_equal(h,k)) then
             r1=dot_product(SpG%SymOp(i)%Tr,real(h))
             r2=real(nint(r1))
             if (abs(r1-r2) > eps_ref) then
                info=.true.
                exit
             end if
          end if
       end do

       return
    End Function Hkl_AbsentI

    !!--++
    !!--++ LOGICAL FUNCTION HKL_ABSENTR
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if the reflection is an absence
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Hkl_AbsentR(H,SpG) Result(Info)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent (in) :: h     ! Input reflection
       Type (Space_Group_Type),     intent (in) :: SpG   ! SpaceGroup
       logical                                  :: info

       !---- Local Variables ----!
       integer                      :: i
       real(kind=cp), dimension(3)  :: k
       real(kind=cp)                :: r1,r2

       info=.false.

       do i=1,SpG%multip
          k = hkl_r(h,SpG%SymOp(i))
          if (hkl_equal(h,k)) then
             r1=dot_product(SpG%SymOp(i)%Tr,h)
             r2=real(nint(r1))
             if (abs(r1-r2) > eps_ref) then
                info=.true.
                exit
             end if
          end if
       end do

       return
    End Function Hkl_AbsentR

    !!----
    !!---- LOGICAL FUNCTION HKL_EQUAL
    !!----
    !!----    Calculate if two reflections are equal
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ LOGICAL FUNCTION HKL_EQUALI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    True if 2 reflections are equal
    !!--++
    !!--++  Update: 20/07/2015
    !!
    Function Hkl_EqualI(H,K) Result (Info)
       !---- Arguments ----!
       integer, dimension(3), intent(in) :: h,k   ! Reflections
       logical                           :: info

       info=.false.
       if (h(1)==k(1) .and. h(2)==k(2) .and. h(3)==k(3)) info=.true.

       return
    End Function Hkl_EqualI

    !!--++
    !!--++ LOGICAL FUNCTION HKL_EQUALR
    !!--++
    !!--++    (OVERLOADED)
    !!--++    True if 2 reflections are equal
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Hkl_EqualR(H,K) Result (Info)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: h,k  ! Reflections
       logical                                 :: info

       info=.false.
       if (abs(h(1)-k(1)) <= eps_ref .and. abs(h(2)-k(2)) <= eps_ref .and. &
           abs(h(3)-k(3)) <= eps_ref) info=.true.

       return
    End Function Hkl_EqualR

    !!----
    !!---- LOGICAL FUNCTION HKL_EQUIV
    !!----
    !!----    Calculate if two reflections are equivalent
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ LOGICAL FUNCTION HKL_EQUIVI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if the reflections are equivalent
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Hkl_EquivI(H,K,SpG,Friedel) Result (Info)
       !---- Arguments ----!
       integer, dimension(3),    intent(in)  :: h,k          ! Reflections
       Type (Space_Group_Type),  intent (in) :: SpG          ! SpG
       logical, optional,        intent(in)  :: Friedel      ! Apply Friedel law
       logical                               :: info

       !---- Local Variables ----!
       integer                           :: i, nops
       integer, dimension(3)             :: hh

       info=.false.

       nops= SpG%numops*max(SpG%centred,1)
       do i=1,nops
          hh = hkl_r(h,SpG%SymOp(i))
          if (hkl_equal(k,hh)) then
             info=.true.
             exit
          end if
          if (present(Friedel)) then
             if (Friedel) then
                if (hkl_equal(k,-hh)) then
                   info=.true.
                   exit
                end if
             end if
          end if
       end do

       return
    End Function Hkl_EquivI

    !!--++
    !!--++ LOGICAL FUNCTION HKL_EQUIVR
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if the reflections are equivalent
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Hkl_EquivR(H,K,SpG,Friedel) Result (Info)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: h,k     ! Reflections
       Type (Space_Group_Type),     intent(in) :: SpG     ! Space group
       logical, optional,           intent(in) :: Friedel ! Aply Friedel law
       logical                                 :: info

       !---- Local Variables ----!
       integer                            :: i, nops
       real(kind=cp), dimension(3)        :: hh

       info=.false.

       nops= SpG%numops*max(SpG%centred,1)
       do i=1, nops
          hh = hkl_r(h,SpG%SymOp(i))
          if (hkl_equal(k,hh)) then
             info=.true.
             exit
          end if
          if (present(Friedel)) then
             if (Friedel) then
                if (hkl_equal(k,-hh)) then
                   info=.true.
                   exit
                end if
             end if
          end if
       end do

       return
    End Function Hkl_EquivR

    !!--++
    !!--++ LOGICAL FUNCTION HKL_LAT_ABSENT
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if the reflection is a lattice absence
    !!--++
    !!--++  Update: 20/07/2015
    !!
    Function Hkl_Lat_Absent(H,Latt) Result(Info)
       !---- Arguments ----!
       integer, dimension(3),        intent (in) :: h      ! Reflection
       Type (Lattice_Centring_Type), intent (in) :: Latt   ! Lattice centring
       logical                                   :: info

       !---- Local Variables ----!
       integer               :: k,i
       logical               :: tinv
       real(kind=cp)         :: r1,r2

       info=.false.
       if (.not. Latt%set) return

       tinv=.false.
       if (ubound(Latt%Ltr,1) == 4) tinv=.true.

       do i=1,Latt%n_lat
          r1=dot_product(Latt%Ltr(1:3,i),real(h))
          r2=real(nint(r1))
          k=nint(2.0*r1)
          if (tinv) then  !Time inversion is considered
             if (Latt%Ltr(4,i) > 0.0) then  !No time inversion, lattice centring
                if (mod(k,2) /= 0) info=.true.
                exit
             else !now time inversion is associated with the translation (Anti-translation)
                if (abs(r1-r2) < eps_ref) info=.true.
                exit
             end if
          else  !No time inversion is considered only normal lattice centring vectors
             if (mod(k,2) /= 0) info=.true.
             exit
          end if
       end do

       return
    End Function Hkl_Lat_Absent

    !!----
    !!---- FUNCTION HKL_MULT
    !!----
    !!----    Calculate the multiplicity of the reflection
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ FUNCTION HKL_MULTI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate the multiplicity of the reflection
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Hkl_MultI(H,SpG,Friedel) Result(N)
       !---- Arguments ----!
       integer, dimension(3),   intent (in) :: h       ! Reflection
       Type (Space_Group_Type), intent (in) :: SpG     ! Space group
       Logical,                 intent (in) :: Friedel ! Apply Friedel law
       integer                              :: N

       !---- Local Variables ----!
       logical                         :: esta
       integer, dimension(3)           :: k
       integer                         :: i,j,ng
       integer, dimension(3,SpG%numops):: klist

       ng=SpG%numops
       n=1

       if (ng > 1) then
          klist(:,1)=h(:)

          do i=2,ng
             k = hkl_r(h,SpG%SymOp(i))
             esta=.false.
             do j=1,n
                if (hkl_equal(k,klist(:,j)) .or. hkl_equal(-k,klist(:,j))) then
                   esta=.true.
                   exit
                end if
             end do
             if (esta) cycle
             n=n+1
             klist(:,n) = k
          end do
       end if

       if (Friedel .or. SpG%centred == 2) then
           n=n*2
       end if

       return
    End Function Hkl_MultI

    !!--++
    !!--++ FUNCTION HKL_MULTR
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate the multiplicity of the reflection
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Hkl_MultR(H,SpG,Friedel) Result(N)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent (in) :: h        ! Reflection
       Type (Space_Group_Type),     intent (in) :: SpG      ! Spacegroup
       Logical,                     intent (in) :: Friedel  ! Friedel law
       integer                                  :: n

       !---- Local Variables ----!
       logical                               :: esta
       integer                               :: i,j,ng
       real(kind=cp), dimension(3)           :: k
       real(kind=cp), dimension(3,SpG%numops):: klist

       ng=SpG%numops
       n=1
       klist(:,1)=h

       do i=2,ng
          k = hkl_r(h,SpG%SymOp(i))
          esta=.false.
          do j=1,n
             if (hkl_equal(k,klist(:,j)) .or. hkl_equal(-k,klist(:,j))) then
                esta=.true.
                exit
             end if
          end do
          if (esta) cycle
          n=n+1
          klist(:,n) = k
       end do
       if (Friedel .or. SpG%centred == 2) n=n*2

       return
    End Function Hkl_MultR

    !!----
    !!---- FUNCTION HKL_R
    !!----
    !!----    Calculate the equivalent reflection
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ FUNCTION HR_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate the equivalent reflection
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function HR_I(H,Op) Result(K)
       !---- Arguments ----!
       integer, dimension(3), intent(in) :: h    ! Reflection
       Type(Sym_Oper_Type),   intent(in) :: Op   ! S.O
       integer, dimension(3)             :: k

       k = matmul(h,Op%RotI)

    End Function HR_I

    !!--++
    !!--++ FUNCTION HR_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate the equivalent reflection
    !!--++
    !!--++ Update: February - 2005
    !!
    Function HR_R(H,Op) Result(K)
       !---- Arguments ----!
       real(kind=cp), dimension(3),  intent(in) :: h    ! Reflection
       Type(Sym_Oper_Type),          intent(in) :: Op   ! S.O
       real(kind=cp), dimension(3)              :: k

       k = matmul(h,Op%RotI)

       return
    End Function HR_R

    !!----
    !!---- Function  Hkl_S
    !!--<<
    !!----    Calculates: sin_theta/lamda = 1/(2d)
    !!-->>
    !!----
    !!----  Update: 20/07/2015
    !!

    !!--++
    !!--++ FUNCTION HS_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate sin_theta/lamda = 1/ (2d)
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function HS_I(H,Crystalcell) Result(S)
       !---- Arguments ----!
       integer, dimension(3),    intent(in)  :: h            ! Reflection
       type (Crystal_Cell_Type), intent (in) :: CrystalCell  ! Cell object
       real(kind=cp)                         :: s

       s= 0.5*sqrt( h(1)*h(1)*CrystalCell%GR(1,1) +     h(2)*h(2)*CrystalCell%GR(2,2) + &
                    h(3)*h(3)*CrystalCell%GR(3,3) + 2.0*h(1)*h(2)*CrystalCell%GR(1,2) + &
                2.0*h(1)*h(3)*CrystalCell%GR(1,3) + 2.0*h(2)*h(3)*CrystalCell%GR(2,3) )

       return
    End Function HS_I

    !!--++
    !!--++ FUNCTION HS_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate sin_theta/lamda = 1/ (2d)
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function HS_R(H,Crystalcell) Result(S)
       !---- Arguments ----!
       real(kind=cp), dimension(3),intent(in)  :: h            ! Reflection
       type (Crystal_Cell_Type),   intent (in) :: CrystalCell  ! Cell object
       real(kind=cp)                           :: s

       s= 0.5*sqrt( h(1)*h(1)*CrystalCell%GR(1,1) +     h(2)*h(2)*CrystalCell%GR(2,2) + &
                    h(3)*h(3)*CrystalCell%GR(3,3) + 2.0*h(1)*h(2)*CrystalCell%GR(1,2) + &
                2.0*h(1)*h(3)*CrystalCell%GR(1,3) + 2.0*h(2)*h(3)*CrystalCell%GR(2,3) )

       return

    End Function HS_R

    !!----
    !!----  FUNCTION UNIT_CART_HKL
    !!----
    !!----     Calculate a unitary vector in the cartesian crystal frame
    !!----     along a reciprocal vector hkl (reciprocal lattice)
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ FUNCTION UNIT_CART_HKLI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate a unitary vector in the cartesian crystal
    !!--++    frame along a reciprocal vector hkl (reciprocal lattice)
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Unit_Cart_HklI(H, Crystalcell) Result (U)
       !---- Arguments ----!
       integer, dimension(3),    intent(in)  :: h            ! Reflection
       type (Crystal_Cell_Type), intent (in) :: CrystalCell  ! Cell object
       real(kind=cp), dimension(3)           :: u

       !---- Local Variables ----!
       real(kind=cp), dimension(3)           :: v

       v=matmul(CrystalCell%GR,real(h))     ![L-2]
       u=matmul(CrystalCell%Cr_Orth_cel,v)  ![L-1]
       u=u/sqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))

       return
    End Function Unit_Cart_HklI


    !!--++
    !!--++ FUNCTION UNIT_CART_HKLR
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate a unitary vector in the cartesian crystal
    !!--++    frame along a reciprocal vector hkl (reciprocal lattice)
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Function Unit_Cart_HklR(H, Crystalcell) Result (U)
       !---- Arguments ----!
       real(kind=cp), dimension(3),intent(in)  :: h            ! Reflection
       type (Crystal_Cell_Type),   intent (in) :: CrystalCell  ! Cell object
       real(kind=cp), dimension(3)             :: u

       !---- Local Variables ----!
       real(kind=cp), dimension(3)             :: v

       v=matmul(CrystalCell%GR,h)
       u=matmul(CrystalCell%Cr_Orth_cel,v)
       u=u/sqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))

       return
    End Function Unit_Cart_HklR

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!--++
    !!--++ SUBROUTINE GLIDE_PLANES_CONDITIONS
    !!--++
    !!--++    Reflections Conditions according with I.T. Table 2.2.13.2
    !!--++    space.
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Glide_Planes_Conditions(SpG,Iunit)
       !---- Arguments ----!
       type (Space_Group_Type), intent(in)     :: SpG
       integer, optional,       intent(in)     :: iunit

       !---- Local variables ----!
       integer               :: h, k,l, m
       integer               :: n, n_ext
       integer, dimension(3) :: hh
       integer               :: num_exti
       logical               :: zonal_condition

       zonal_condition   = .false.

       if (present(iunit) ) then
          write(unit=iunit,fmt=*) " "
          write(unit=iunit,fmt=*) " >>> Zonal reflections conditions for glide planes:"
          write(unit=iunit,fmt=*) "---------------------------------------------------"
          write(unit=iunit,fmt=*) " "
       end if

       !GLIDE PLANES and screw axes: table 2.13.2
       !-------------
       !
       !        0 k l:    k=2n    b/2             monoclinic, orthorhombic, tetragonal and cubic
       !        0 k l:    l=2n    c/2             monoclinic, orthorhombic, tetragonal and cubic
       !        0 k l:  k+l=2n    b/2 +  c/2      monoclinic, orthorhombic, tetragonal and cubic
       !        0 k l:  k+l=4n    b/4 +- c/4      orthorhombic and cubic
       !
       !
       !        h 0 l:    h=2n    a/2             monoclinic, orthorhombic, tetragonal and cubic
       !        h 0 l:    l=2n    c/2             monoclinic, orthorhombic, tetragonal and cubic
       !        h 0 l:  l+h=2n    c/2 +  a/2      monoclinic, orthorhombic, tetragonal and cubic
       !        h 0 l:  l+h=4n    c/4 +- a/4      orthorhombic and cubic
       !
       !        h k 0:    h=2n    a/2             monoclinic, orthorhombic, tetragonal and cubic
       !        h k 0:    k=2n    b/2             monoclinic, orthorhombic, tetragonal and cubic
       !        h k 0:  h+k=2n    a/2 +  b/2      monoclinic, orthorhombic, tetragonal and cubic
       !        h k 0:  h+k=4n    a/4 +- b/4      monoclinic, orthorhombic, tetragonal and cubic

       if (SpG%CrystalSys(1:10) == "Monoclinic"   .or. SpG%CrystalSys(1:10) == "Tetragonal" .or.     &
           SpG%CrystalSys(1:12) == "Orthorhombic" .or. SpG%CrystalSys(1:5)  == "Cubic") then

          !---- glide plane b/2:
          ! Hkl_Ref_Conditions(7)  =   "(0 k l)      k=2n : 0yz glide plane with b/2 translation"
          num_exti = 7
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do k=-6, 6
             do l=-6, 6
                hh(1)=0
                hh(2)=k
                hh(3)=l
                m =  k
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane c/2:
          ! Hkl_Ref_Conditions(8)  =   "(0 k l)      l=2n : 0yz glide plane with c/2 translation"
          num_exti = 8
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do k=-6, 6
             do l=-6, 6
                hh(1)=0
                hh(2)=k
                hh(3)=l
                m =  l
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane b/2 + c/2:
          !Hkl_Ref_Conditions(9)  =   "(0 k l)    k+l=2n : 0yz glide plane with b/2 + c/2 translation"
          num_exti = 9
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do k=-6, 6
             do l=-6, 6
                hh(1)=0
                hh(2)=k
                hh(3)=l
                m =  k+l
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if    ! fin de la condition "if monoclinic, tetragonal, ortho, cubic


       if (SpG%CrystalSys(1:12) == "Orthorhombic" .or. SpG%CrystalSys(1:5)  == "Cubic") then
          !---- glide plane b/4 + c/4:
          ! Hkl_Ref_Conditions(10)  =   "(0 k l)    k+l=4n : 0yz glide plane with b/4 +- c/4 translation"
          num_exti = 10
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do k=-6, 6, 1
             do l=-6, 6, 1
                hh(1)=0
                hh(2)=k
                hh(3)=l
                m =  k+l
                if (m /= int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if ! fin de la condition "if ortho, cubic

       if (SpG%CrystalSys(1:10) == "Monoclinic"   .or. SpG%CrystalSys(1:10) == "Tetragonal" .or.     &
          SpG%CrystalSys(1:12) == "Orthorhombic" .or. SpG%CrystalSys(1:5)  == "Cubic") then

          !---- glide plane a/2:
          !  Hkl_Ref_Conditions(11)  =   "(h 0 l)      h=2n : x0z glide plane with a/2 translation"
          num_exti = 11
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=0
                hh(3)=l
                m =  h
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane c/2:
          ! Hkl_Ref_Conditions(12) =   "(h 0 l)      l=2n : x0z glide plane with c/2 translation"
          num_exti = 12
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=0
                hh(3)=l
                m =  l
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane c/2 + a/2:
          ! Hkl_Ref_Conditions(13) =   "(h 0 l)    l+h=2n : x0z glide plane with a/2 + c/2 translations"
          num_exti = 13
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=0
                hh(3)=l
                m =  h+l
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if  ! fin de la condition "if monoclinic, tetragonal, ortho, cubic

       if (SpG%CrystalSys(1:12) == "Orthorhombic" .or. SpG%CrystalSys(1:5)  == "Cubic") then

          !---- glide plane c/4 + a/4:
          ! Hkl_Ref_Conditions(14) =   "(h 0 l)    l+h=4n : x0z glide plane with a/4 +- c/4 translations"
          num_exti = 14
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=0
                hh(3)=l
                m =  h+l
                if (m /= int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if ! fin de la condition "if ortho, cubic

       if (SpG%CrystalSys(1:10) == "Monoclinic"   .or. SpG%CrystalSys(1:10) == "Tetragonal" .or.     &
          SpG%CrystalSys(1:12) == "Orthorhombic" .or. SpG%CrystalSys(1:5)  == "Cubic") then

          !---- glide plane a/2:
          ! Hkl_Ref_Conditions(15) =   "(h k 0)      h=2n : xy0 glide plane with a/2 translation"
          num_exti = 15
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6
             do k=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=0
                m =  h
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do    ! k loop
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane b/2:
          !Hkl_Ref_Conditions(16) =   "(h k 0)      k=2n : xy0 glide plane with b/2 translation"
          num_exti = 16
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6
             do k=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=0
                m =  k
                if (m /= int(m/2)*2) then
                   n=n+1
                  if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do    ! k loop
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane a/2 + b/2:
          ! Hkl_Ref_Conditions(17) =   "(h k 0)    h+k=2n : xy0 glide plane with a/2 + b/2 translations"
          num_exti = 17
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6
             do k=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=0
                m =  h+k
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do    ! k loop
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if  ! fin de la condition "if monoclinic, tetragonal, ortho, cubic

       if (SpG%CrystalSys(1:12) == "Orthorhombic" .or. SpG%CrystalSys(1:5)  == "Cubic") then
          !---- glide plane a/4 + b/4:
          ! Hkl_Ref_Conditions(18) =   "(h k 0)    h+k=4n : xy0 glide plane with a/4 +- b/4 translations"
          num_exti = 18
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6, 1
             do k=-6, 6, 1
                hh(1)=h
                hh(2)=k
                hh(3)=0
                m =  h+k
                if (m /= int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do    ! k loop
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if  ! fin de la condition "if ortho, cubic

       if (SpG%SPG_Latsy(1:1) == "h") then
          !---- glide plane with c/2 translation: hexagonal
          !  Hkl_Ref_Conditions(19) =   "(  h  -h   0 l) l=2n : (11-20) glide plane with c/2 translation (c)"
          num_exti = 19
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=h
                hh(2)=-h
                hh(3)=l
                m=l
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane with c/2 translation: hexagonal
          !  Hkl_Ref_Conditions(20) =   "(  0   k  -k l) l=2n : (-2110) glide plane with c/2 translation (c)"
          num_exti = 20
          n = 0
          n_ext = 0
          do k=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=0
                hh(2)=k
                hh(3)=l
                m=l
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane with c/2 translation: hexagonal
          !Hkl_Ref_Conditions(21) =   "( -h   0   h l) l=2n : (1-210) glide plane with c/2 translation (c)"
          num_exti = 21
          n = 0
          n_ext = 0
          do h=-6, 6
             do l=-6, 6
                hh(1)=-h
                hh(2)=0
                hh(3)=l
                m=l
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane with c/2 translation: hexagonal
          ! Hkl_Ref_Conditions(22) =   "(  h   h -2h l) l=2n : (1-100) glide plane with c/2 translation (c)"
          num_exti = 22
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=h
                hh(2)=h
                hh(3)=l
                m=l
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane with c/2 translation: hexagonal
          !  Hkl_Ref_Conditions(23) =   "(-2h   h   h l) l=2n : (01-10) glide plane with c/2 translation (c)"
          num_exti = 23
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=-2*h
                hh(2)=h
                hh(3)=l
                m=l
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane with c/2 translation: hexagonal
          !  Hkl_Ref_Conditions(24) =   "(  h -2h   h l) l=2n : (-1010) glide plane with c/2 translation (c)"
          num_exti = 24
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=h
                hh(2)=-2*h
                hh(3)=l
                m=l
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if ! fin de la condition if hexagonal

       !25: glide plane with c/2 translation: rhomboedral
       !  Hkl_Ref_Conditions(25) =  "(  h  h  l) l=2n : (1-10) glide plane with c/2 translation (c,n)"
       num_exti = 25
       n = 0
       n_ext = 0
       do h=-6, +6, 1
          do l=-6, +6, 1
             hh(1)=h
             hh(2)=h
             hh(3)=l
             m=l
             if (m /=int(m/2)*2) then
                n=n+1
                if (hkl_absent(hh, SpG)) n_ext=n_ext+1
             end if
          end do  ! l loop
       end do   ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          zonal_condition = .true.
       end if

       !---- glide plane with c/2 translation: rhomboedral
       !  Hkl_Ref_Conditions(26) =  "(  h  k  k) h=2n : (01-1) glide plane with a/2 translation (a,n)"
       num_exti = 26
       n = 0
       n_ext = 0
       do h=-6, +6, 1
          do k=-6, +6, 1
             hh(1)=h
             hh(2)=k
             hh(3)=k
             m=h
             if (m /=int(m/2)*2) then
                n=n+1
                if (hkl_absent(hh, SpG)) n_ext=n_ext+1
             end if
          end do  ! l loop
       end do   ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          zonal_condition = .true.
       end if

       !27: glide plane with c/2 translation: rhomboedral
       !  Hkl_Ref_Conditions(27) =  "(  h  k  h) k=2n : (-101) glide plane with b/2 translation (b,n)"
       num_exti = 27
       n = 0
       n_ext = 0
       do h=-6, +6, 1
          do k=-6, +6, 1
             hh(1)=h
             hh(2)=k
             hh(3)=h
             m=k
             if (m /=int(m/2)*2) then
                n=n+1
                if (hkl_absent(hh, SpG)) n_ext=n_ext+1
             end if
          end do  ! l loop
       end do   ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          zonal_condition = .true.
       end if

       if (SpG%CrystalSys(1:10) == "Tetragonal" .or. SpG%CrystalSys(1:5) == "Cubic") then
          !---- glide plane with c/2 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(28) =  "(  h  h  l)    l=2n : (1-10) glide plane with c/2 translation (c,n)"
          num_exti = 28
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=h
                hh(2)=h
                hh(3)=l
                m=l
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane with a/4 +- b/4 +- c/4 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(29) =  "(  h  h  l) 2h+l=4n : (1-10) glide plane with a/4 +- b/4 +- c/4 translation (d)"
          num_exti = 29
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=h
                hh(2)=h
                hh(3)=l
                m=2*h+l
                if (m /=int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !30: glide plane with c/2 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(30) =  "(  h -h  l)    l=2n : (110)  glide plane with c/2 translation (c,n)"
          num_exti = 30
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=h
                hh(2)=-h
                hh(3)=l
                m=l
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          ! 31: glide plane with a/4 +- b/4 +- c/4 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(31) = "(  h -h  l) 2h+l=4n : (110)  glide plane with a/4 +- b/4 +- c/4 translation (d)"
          num_exti = 31
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do l=-6, +6, 1
                hh(1)=h
                hh(2)=-h
                hh(3)=l
                m=2*h+l
                if (m /=int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if   ! fin de la condition "if tetragonal .or. cubic

       if (SpG%CrystalSys(1:5) == "Cubic") then
          !---- glide plane with a/2 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(32) = "(  h  k  k)    h=2n : (01-1) glide plane with a/2 translation (a,n)"
          num_exti = 32
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do k=-6, +6, 1
                hh(1)=h
                hh(2)=k
                hh(3)=k
                m=h
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !---- glide plane with +-a/4 +- b/4 +- c/4 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(33) = "(  h  k  k) 2k+h=4n : (01-1) glide plane with +-a/4 + b/4 +- c/4 translation (d)"
          num_exti = 33
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do k=-6, +6, 1
                hh(1)=h
                hh(2)=k
                hh(3)=k
                m=2*k+h
                if (m /=int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !34: glide plane with a/2 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(34) =  "(  h  k -k)    h=2n : (011)  glide plane with a/2 translation (a,n)"
          num_exti = 34
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do k=-6, +6, 1
                hh(1)=h
                hh(2)=k
                hh(3)=-k
                m=h
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          ! 35: glide plane with a/4 +- b/4 +- c/4 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(351) = "(  h  k -k) 2k+h=4n : (011)  glide plane with +-a/4 + b/4 +- c/4 translation (d)"
          num_exti = 35
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do k=-6, +6, 1
                hh(1)=h
                hh(2)=k
                hh(3)=-k
                m=2*k+h
                if (m /=int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !36: glide plane with b/2 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(36) = "(  h  k  h)    k=2n : (-101) glide plane with b/2 translation (b,n)"
          num_exti = 36
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do k=-6, +6, 1
                hh(1)=h
                hh(2)=k
                hh(3)=h
                m=k
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !37: glide plane with +-a/4 +- b/4 +- c/4 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(33) = "(  h  k  h) 2h+k=4n : (-101) glide plane with +-a/4 + b/4 +- c/4 translation (d)"
          num_exti = 37
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do k=-6, +6, 1
                hh(1)=h
                hh(2)=k
                hh(3)=h
                m=2*h+k
                if (m /=int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          !38: glide plane with b/2 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(38) = "( -h  k  h)    k=2n : (101)  glide plane with b/2 translation (b,n)"
          num_exti = 38
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do k=-6, +6, 1
                hh(1)=-h
                hh(2)=k
                hh(3)=h
                m=k
                if (m /=int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if

          ! 39: glide plane with a/4 +- b/4 +- c/4 translation: tetragonal + cubic
          !  Hkl_Ref_Conditions(39) = "( -h  k  h) 2h+k=4n : (101)  glide plane with +-a/4 + b/4 +- c/4 translation (d)"
          num_exti = 39
          n = 0
          n_ext = 0
          do h=-6, +6, 1
             do k=-6, +6, 1
                hh(1)=-h
                hh(2)=k
                hh(3)=h
                m=2*h+k
                if (m /=int(m/4)*4) then
                   n=n+1
                   if (hkl_absent(hh, SpG)) n_ext=n_ext+1
                end if
             end do  ! l loop
          end do   ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             zonal_condition = .true.
          end if
       end if  ! fin de la condition "if cubic

       if (.not. zonal_condition)   then
          if (present(iunit)) write(unit=iunit,fmt=*) "     =====>>> no zonal reflection condition"
       end if

       return
    End Subroutine Glide_Planes_Conditions

    !!----
    !!---- SUBROUTINE HKL_EQUIV_LIST
    !!----
    !!----    Calculate the multiplicity of the reflection and the list of all
    !!----    equivalent reflections. Friedel law assumed if Friedel=.true.
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE HKL_EQUIV_LISTI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate the multiplicity of the reflection and the list of all
    !!--++    equivalent reflections. Friedel law assumed if Friedel=.true.
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Hkl_Equiv_Listi(H,SpG,Friedel,Mul,Hlist)
       !---- Arguments ----!
       integer, dimension(3),                     intent (in) :: h        ! Reflection
       Type (Space_Group_Type),                   intent (in) :: SpG      ! Space group
       Logical,                                   intent (in) :: Friedel  ! Friedel law
       integer,                                   intent(out) :: mul      ! Multiplicity
       integer, dimension(3,SpG%numops*2),        intent(out) :: hlist    ! List of equivalent reflections

       !---- Local Variables ----!
       logical              :: esta
       integer, dimension(3):: k
       integer              :: i,j,ng

       !> Init
       hlist = 0

       ng=SpG%numops
       mul=1
       hlist(:,1)=h
       do i=2,ng
          k = hkl_r(h,SpG%SymOp(i))
          esta=.false.
          do j=1,mul
             if (hkl_equal(k,hlist(:,j)) .or. (hkl_equal(-k,hlist(:,j)) .and. Friedel)) then
                esta=.true.
                exit
             end if
          end do
          if (esta) cycle
          mul=mul+1
          hlist(:,mul) = k
       end do

       if (Friedel .or. SpG%centred == 2) then
          j=mul
          mul=mul*2
          do i=j+1,mul
             hlist(:,i)=-hlist(:,i-j)
          end do
       end if

       return
    End Subroutine Hkl_Equiv_Listi

    !!--++
    !!--++ SUBROUTINE HKL_EQUIV_LISTR
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate the multiplicity of the reflection and the list of all
    !!--++    equivalent reflections. Friedel law assumed if Friedel=.true.
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Hkl_Equiv_ListR(H,SpG,Friedel,Mul,Hlist)
       !---- Arguments ----!
       real(kind=cp), dimension(3),                     intent (in)  :: h         ! Reflection
       Type (Space_Group_Type),                         intent (in)  :: SpG       ! Space group
       Logical,                                         intent (in)  :: Friedel   ! Friedel law
       integer,                                         intent (out) :: mul       ! Multiplicity
       real(kind=cp), dimension(3,SpG%numops*2),        intent (out) :: hlist     ! List of equivalent reflections

       !---- Local Variables ----!
       logical                    :: esta
       real(kind=cp), dimension(3):: k
       integer                    :: i,j,ng

       !> Init
       hlist = 0.0

       ng=SpG%numops
       mul=1
       hlist(:,1)=h
       do i=2,ng
          k = hkl_r(h,SpG%SymOp(i))
          esta=.false.
          do j=1,mul
             if (hkl_equal(k,hlist(:,j)) .or. (hkl_equal(-k,hlist(:,j)) .and. Friedel)) then
                esta=.true.
                exit
             end if
          end do
          if (esta) cycle
          mul=mul+1
          hlist(:,mul) = k
       end do
       if (Friedel .or. SpG%centred == 2) then
          j=mul
          mul=mul*2
          do i=j+1,mul
             hlist(:,i)=-hlist(:,i-j)
          end do
       end if

       return
    End Subroutine Hkl_Equiv_Listr

    !!----
    !!---- SUBROUTINE HKL_GEN
    !!----
    !!----    Calculate unique reflections between two values of
    !!----    sin_theta/lambda.  The output is not ordered.
    !!----
    !!---- Update: 20/07/2015
    !!
    Subroutine Hkl_Gen(Crystalcell,SpG,Friedel,Value1,Value2,Num_Ref,Reflex)
       !---- Arguments ----!
       type (Crystal_Cell_Type),          intent(in)     :: crystalcell      ! Cell object
       type (Space_Group_Type) ,          intent(in)     :: SpG              ! Space group object
       Logical,                           intent(in)     :: Friedel          ! If true, Friedel law applied
       real(kind=cp),                     intent(in)     :: value1,value2    ! Range in SinTheta/Lambda
       integer,                           intent(out)    :: num_ref          ! Number of generated reflections
       type (Reflect_Type), dimension(:), intent(out)    :: reflex           ! List of generated hkl,mult, s

       !---- Local variables ----!
       real(kind=cp)         :: vmin,vmax,sval
       integer               :: h,k,l,hmin,kmin,lmin,hmax,kmax,lmax, maxref
       integer, dimension(3) :: hh,kk,nulo
       character(len=2)      :: inf

       nulo=0
       maxref=size(reflex)
       vmin=min(value1,value2)
       vmax=max(value1,value2)
       hmax=nint(CrystalCell%cell(1)*2.0*vmax+1.0)
       kmax=nint(CrystalCell%cell(2)*2.0*vmax+1.0)
       lmax=nint(CrystalCell%cell(3)*2.0*vmax+1.0)
       lmin= 0  ! l positive or zero except for -3 1 m (see below)

       !---- Select approximate region to generate reflections depending
       !---- on the space group. This allows a faster generation.
       select Case(SpG%NumSpg)
          case (1:2)                 ! -1    -> hkl: l >=0; hk0: h >=0; 0k0: k >=0
             hmin=-hmax
             kmin=-kmax

          case (3:15)                ! 2/m
             inf(1:2)=adjustl(SpG%info(1:2))
             if(inf(1:1) == "-") inf(1:1)=inf(2:2)
             select case (inf(1:1))
                case ("b")     !       -> hkl: k >=0, l >=0; hk0: h >=0
                   hmin=-hmax
                   kmin=0
                case ("c")     !       -> hkl: k >=0, l >=0; h0l: h >=0
                   hmin=-hmax
                   kmin=0
                case ("a")     !       -> hkl: h >=0, l >=0; 0kl: l >=0  Provisional (to be tested)
                   kmin=-kmax
                   hmin=0
                case default
                   hmin=-hmax
                   kmin=0
             end select

          case (16:74)         ! mmm   -> hkl: h >=0, k >=0, l >=0
             hmin=0
             kmin=0

          case (75:88)         ! 4/m   -> hkl: h >=0, l >=0, k >=0 if h = 0
                               !                             k > 0 if h > 0
             hmin=0
             kmin=0

          case (89:142)        ! 4/mmm -> hkl: h >=0, k>=0, l>=0, h >=k
             hmin=0
             kmin=0

          case (143:148)       ! -3    -> hkl: h+k>0, l>0 ;  hk0: h>0, k>=0
             hmin=0
             kmin=-kmax

          case (149,151,153,157,159,162,163) ! -3 1 m  -> hkl: h>=0,h>=k>0 ; h0l: h>=0,l>=0
             hmin=0
             kmin=0
             lmin=-lmax

          case (150,152,154,155,156,158,160,161,164,165,166,167)
                              ! -3 m   -> hkl: h>=0 h>=k ; hhl: h>=0,l>=0
             hmin=0
             kmin=0

          case (168:176)    ! 6/m   -> hkl: h>0,k>0,l>=0;  0kl k>=0,l>=0
             hmin=0
             kmin=0

          case (177:194)    ! 6/mmm -> hkl: h >=0, k >=0, l >=0, h >=k
             hmin=0
             kmin=0

          case (195:206)    ! m-3   -> hkl: h > l, k > l, l >=0 ; hkk: k>=0 h>=k
             hmin=0
             kmin=0

          case (207:230)    ! m-3m  -> hkl: h >=0, k >=0, l >=0, h >=k, k >=l
             hmin=0
             kmin=0

          case default      ! Assumed -1
             hmin=-hmax
             kmin=-kmax
       end Select

       num_ref=0
       ext_do: do h=hmin,hmax
          do k=kmin,kmax
             do l=lmin,lmax

                hh(1)=h
                hh(2)=k
                hh(3)=l

                if (hkl_equal(hh,nulo)) cycle
                sval=hkl_s(hh,crystalcell)
                if (sval > vmax .or. sval < vmin) cycle
                if (hkl_absent(hh,SpG)) cycle

                kk=asu_hkl(hh,SpG)
                if (hkl_equal(kk,nulo)) cycle
                if (hkl_equal(kk,-hh)) cycle

                num_ref=num_ref+1
                if(num_ref > maxref) then
                   num_ref=maxref
                   exit ext_do
                end if
                reflex(num_ref)%h    = kk
                reflex(num_ref)%mult = hkl_mult(kk,SpG,friedel)
                reflex(num_ref)%S    = sval
             end do
          end do
       end do ext_do

       return
    End Subroutine Hkl_Gen

    !!----
    !!---- SUBROUTINE HKL_GEN_SXTAL
    !!----
    !!----    Calculate all allowed reflections between a minimum and a maximum value of sin_theta/lambda.
    !!----    If the limits of indices are provided in hlim, only the reflections verifying the prescription
    !!----    are finally kept. hlim(:,1) and hlim(:,2) contain the minimum and maximum values respectively.
    !!----    The output is not ordered but the user can obtain the reflections generated
    !!----    in a particular way by providing the integer vector "ord", containing a permutation
    !!----    of the three numbers 1,2,3. By default the loop generating the hkl-indices uses
    !!----    the vector ord=(/3,2,1/), this means that the inner loop (more rapidly changing index)
    !!----    is the l-index, then the k-index and finally the h-index.
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE HKL_GEN_SXTAL_REFLECTION
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate all allowed reflections between a minimum and a maximum value of sin_theta/lambda.
    !!--++    If the limits of indices are provided in hlim, only the reflections verifying the prescription
    !!--++    are finally kept. hlim(:,1) and hlim(:,2) contain the minimum and maximum values respectively.
    !!--++    The reflections are stored in the array Reflex, with components of type: Reflect_Type
    !!--++    The output is not ordered but the user can obtain the reflections generated
    !!--++    in a particular way by providing the integer vector "ord", containing a permutation
    !!--++    of the three numbers 1,2,3. By default the loop generating the hkl-indices uses
    !!--++    the vector ord=(/3,2,1/), this means that the inner loop (more rapidly changing index)
    !!--++    is the l-index, then the k-index and finally the h-index.
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Hkl_Gen_Sxtal_Reflection(Crystalcell,SpG,stlmin,stlmax,Num_Ref,Reflex,ord,hlim)
       !---- Arguments ----!
       type (Crystal_Cell_Type),          intent(in)  :: crystalcell       ! Cell object
       type (Space_Group_Type) ,          intent(in)  :: SpG               ! Space group object
       real(kind=cp),                     intent(in)  :: stlmin,stlmax     ! Minimum and Maximum SinTheta/Lambda
       integer,                           intent(out) :: num_ref           ! Number of generated reflections
       type (Reflect_Type), dimension(:), intent(out) :: reflex            ! List of generated hkl,mult, s
       Integer, dimension(3),   optional, intent(in)  :: ord               ! Order for loop of hkl-indices
       Integer, dimension(3,2), optional, intent(in)  :: hlim              ! hkl-limits

       !---- Local variables ----!
       real(kind=cp)         :: sval
       integer               :: h,k,l,hmax,kmax,lmax, maxref
       integer, dimension(3) :: hh,nulo,od,imin,imax

       nulo=0
       maxref=size(reflex)
       hmax=nint(CrystalCell%cell(1)*2.0*stlmax+1.0)
       kmax=nint(CrystalCell%cell(2)*2.0*stlmax+1.0)
       lmax=nint(CrystalCell%cell(3)*2.0*stlmax+1.0)
       if (present(hlim)) then
          imin=hlim(:,1)
          imax=hlim(:,2)
       else
          imin=(/-hmax,-kmax,-lmax/)
          imax=(/ hmax, kmax, lmax/)
       end if
       od=(/3,2,1/)
       if (present(ord)) od=ord

       num_ref=0
       ext_do: do h=imin(od(3)),imax(od(3))
          do k=imin(od(2)),imax(od(2))
             do l=imin(od(1)),imax(od(1))
                hh(od(3))=h
                hh(od(2))=k
                hh(od(1))=l
                if (hkl_equal(hh,nulo)) cycle
                sval=hkl_s(hh,crystalcell)
                if (sval > stlmax .or. sval < stlmin) cycle
                if (hkl_absent(hh,SpG)) cycle
                num_ref=num_ref+1
                if(num_ref > maxref) then
                   num_ref=maxref
                   exit ext_do
                end if
                reflex(num_ref)%h    = hh
                reflex(num_ref)%mult = hkl_mult(hh,SpG,.false.)
                reflex(num_ref)%S    = sval
            end do
          end do
       end do ext_do

       return
    End Subroutine Hkl_Gen_Sxtal_Reflection

    !!--++
    !!--++ SUBROUTINE HKL_GEN_SXTAL_LIST
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate all allowed reflections between a minimum and a maximum value of sin_theta/lambda.
    !!--++    If the limits of indices are provided in hlim, only the reflections verifying the prescription
    !!--++    are finally kept. hlim(:,1) and hlim(:,2) contain the minimum and maximum values respectively.
    !!--++    The reflections are stored in the scalar object Reflex of type: Reflection_List_Type
    !!--++    The output is not ordered but the user can obtain the reflections generated
    !!--++    in a particular way by providing the integer vector "ord", containing a permutation
    !!--++    of the three numbers 1,2,3. By default the loop generating the hkl-indices uses
    !!--++    the vector ord=(/3,2,1/), this means that the inner loop (more rapidly changing index)
    !!--++    is the l-index, then the k-index and finally the h-index.
    !!--++
    !!--++ Update: May - 2006
    !!
    Subroutine Hkl_Gen_Sxtal_List(Crystalcell,SpG,stlmin,stlmax,Num_Ref,Reflex,ord,hlim)
       !---- Arguments ----!
       type (Crystal_Cell_Type),          intent(in)  :: crystalcell     ! Unit cell object
       type (Space_Group_Type) ,          intent(in)  :: SpG             ! Space group
       real(kind=cp),                     intent(in)  :: stlmin,stlmax   ! Minimum and Maximum SinTheta/Lambda
       integer,                           intent(out) :: num_ref         ! Number of generated reflections
       Type(Reflection_List_Type),        intent(out) :: reflex          ! Generated set of reflections
       Integer, dimension(3),   optional, intent(in)  :: ord             ! Order for loop of hkl-indices
       Integer, dimension(3,2), optional, intent(in)  :: hlim            ! hkl-limits

       !---- Local variables ----!
       real(kind=cp)         :: sval
       integer               :: h,k,l,hmax,kmax,lmax, maxref,i
       integer, dimension(3) :: hh,nulo,od,imin,imax
       Type(Reflection_Type), dimension(:), allocatable :: tmp_reflex

       nulo=0

       hmax=nint(CrystalCell%cell(1)*2.0*stlmax+1.0)
       kmax=nint(CrystalCell%cell(2)*2.0*stlmax+1.0)
       lmax=nint(CrystalCell%cell(3)*2.0*stlmax+1.0)
       if(present(hlim)) then
         imin=hlim(:,1)
         imax=hlim(:,2)
       else
         imin=(/-hmax,-kmax,-lmax/)
         imax=(/ hmax, kmax, lmax/)
       end if
       od=(/3,2,1/)
       if(present(ord)) od=ord

       maxref=(2*hmax+1)*(2*kmax+1)*(2*lmax+1)
       if(allocated(tmp_reflex)) deallocate(tmp_reflex)
       allocate(tmp_reflex(maxref))

       num_ref=0
       ext_do: do h=imin(od(3)),imax(od(3))
          do k=imin(od(2)),imax(od(2))
             do l=imin(od(1)),imax(od(1))
                hh(od(3))=h
                hh(od(2))=k
                hh(od(1))=l
                if (hkl_equal(hh,nulo)) cycle
                sval=hkl_s(hh,crystalcell)
                if (sval > stlmax .or. sval < stlmin) cycle
                if (hkl_absent(hh,SpG)) cycle
                num_ref=num_ref+1
                if(num_ref > maxref) then
                   num_ref=maxref
                   exit ext_do
                end if

                tmp_reflex(num_ref)%h    = hh
                tmp_reflex(num_ref)%mult = hkl_mult(hh,SpG,.false.)
                tmp_reflex(num_ref)%S    = sval

             end do
          end do
       end do ext_do

       if(allocated(reflex%ref)) deallocate(reflex%ref)
       allocate(reflex%ref(num_ref))
       reflex%nref= num_ref

       do i=1,num_ref
          reflex%Ref(i)%h    = tmp_reflex(i)%h
          reflex%Ref(i)%mult = tmp_reflex(i)%mult
          reflex%Ref(i)%S    = tmp_reflex(i)%S
          reflex%Ref(i)%fo    =0.0
          reflex%Ref(i)%sfo   =0.0
          reflex%Ref(i)%fc    =0.0
          reflex%Ref(i)%w     =0.0
          reflex%Ref(i)%phase =0.0
          reflex%Ref(i)%a     =0.0
          reflex%Ref(i)%b     =0.0
          reflex%Ref(i)%aa    =0.0
          reflex%Ref(i)%bb    =0.0
       end do

       return
    End Subroutine Hkl_Gen_Sxtal_list


    !!----
    !!---- SUBROUTINE  HKL_RP
    !!----
    !!----    Calculate the equivalent reflection and Phase
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE HKL_RPI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate the equivalent reflection and new phase
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Hkl_RpI(H, Phase, Op, K, Phasen)
       !---- Arguments ----!
       integer, dimension(3), intent (in) :: h        ! Reflection
       real(kind=cp),         intent (in) :: phase    ! Phase
       Type(Sym_Oper_Type),   intent (in) :: Op       ! S.O
       integer, dimension(3), intent (out):: k        ! Equivalent reflection
       real(kind=cp),         intent (out):: phasen   ! Phase

       k = matmul(h,Op%RotI)
       phasen= phase - 360.0_cp*dot_product(Op%Tr,real(h))
       phasen=mod(phasen+3600.0_cp,360.0_cp)

       return
    End Subroutine Hkl_RpI

    !!--++
    !!--++ SUBROUTINE HKL_RPR
    !!--++
    !!--++     (OVERLOADED)
    !!--++     Calculate the equivalent reflection and new phase
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Hkl_RpR(h, phase, Op, k, phasen)
       !---- Arguments ----!
       real(kind=cp), dimension(3),intent (in) :: h            ! Reflection
       real(kind=cp),              intent (in) :: phase        ! Phase
       Type(Sym_Oper_Type),        intent(in)  :: Op           ! S.O
       real(kind=cp), dimension(3),intent (out):: k            ! Equivalent reflection
       real(kind=cp),              intent (out):: phasen       ! Phase

       k = matmul(h,Op%RotI)
       phasen= phase - 360.0_cp*dot_product(Op%Tr,h)
       phasen=mod(phasen+3600.0_cp,360.0_cp)

       return
    End Subroutine Hkl_RpR

    !!----
    !!---- SUBROUTINE HKL_UNI
    !!----
    !!----    Calculate unique reflections between two values (value1,value2)
    !!----    of sin_theta/lambda. If no_order is present and .true. the sort subroutine
    !!----    is not invoked.
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE HKL_UNI_REFLECT
    !!--++
    !!--++    (Overloaded)
    !!--++    Calculate unique reflections between two values (value1,value2)
    !!--++    of sin_theta/lambda. If no_order is present and .true. the sort subroutine
    !!--++    is not invoked.
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Hkl_Uni_Reflect(Crystalcell,SpG,Friedel,Value1,Value2,Code,Num_Ref,Reflex,no_order)
       !---- Arguments ----!
       type (Crystal_Cell_Type),             intent(in)     :: crystalcell     ! Cell object
       type (Space_Group_Type) ,             intent(in)     :: SpG             ! Space group object
       Logical,                              intent(in)     :: Friedel         ! If true Friedel law is applied
       real(kind=cp),                        intent(in)     :: value1,value2   ! Range in sintheta/Lambda
       character(len=1),                     intent(in)     :: code            ! If code="r", d-spacing are input in value1/value2
       integer,                              intent(out)    :: num_ref         ! Number of generated reflections
       type (Reflect_Type),    dimension(:), intent(out)    :: reflex          ! Reflections
       logical,                   optional,  intent(in)     :: no_order

       !---- Local variables ----!
       real(kind=cp)                         :: vmin,vmax,sval
       integer                               :: h,k,l,hmin,kmin,lmin,hmax,kmax,lmax, i, maxref
       integer, dimension(3)                 :: hh,kk,nulo
       integer,  dimension(  size(reflex))   :: ind
       integer,  dimension(  size(reflex))   :: mul
       integer,  dimension(3,size(reflex))   :: hkl
       real(kind=cp),dimension(size(reflex)) :: sv
       character(len=2)                      :: inf

       nulo=0
       maxref=size(reflex)
       vmin=min(value1,value2)
       vmax=max(value1,value2)
       if (code =="r" .or. code=="R") then
          vmin=1.0/(2.0*max(value1,value2))
          vmax=1.0/(2.0*min(value1,value2))
       end if

       hmax=nint(CrystalCell%cell(1)*2.0*vmax+1.0)
       kmax=nint(CrystalCell%cell(2)*2.0*vmax+1.0)
       lmax=nint(CrystalCell%cell(3)*2.0*vmax+1.0)
       lmin= 0  !l positive or zero except for -3 1 m (see below)

       !---- Select approximate region to generate reflections depending
       !---- on the space group. This allows a faster generation.
       Select Case(SpG%NumSpg)
          case (1:2)                   ! -1    -> hkl: l >=0; hk0: h >=0; 0k0: k >=0
             hmin=-hmax
             kmin=-kmax
             if(SpG%NumSpg == 1 .and. .not. Friedel) lmin=-lmax

          case (3:15)                  ! 2/m
             inf=SpG%info(1:2)
             if (inf(1:1) == "-") inf(1:1)=inf(2:2)
             select case (inf(1:1))
                case ("b")     !       -> hkl: k >=0, l >=0; hk0: h >=0
                   hmin=-hmax
                   kmin=0
                case ("c")     !       -> hkl: k >=0, l >=0; h0l: h >=0
                   hmin=-hmax
                   kmin=0
                case ("a")     !       -> hkl: h >=0, l >=0; 0kl: l >=0  Provisional (to be tested)
                   kmin=-kmax
                   hmin=0
             end select

          case (16:74)         ! mmm   -> hkl: h >=0, k >=0, l >=0
             hmin=0
             kmin=0

          case (75:88)         ! 4/m   -> hkl: h >=0, l >=0, k >=0 if h = 0
                               !                             k > 0 if h > 0
             hmin=0
             kmin=0

          case (89:142)        ! 4/mmm -> hkl: h >=0, k>=0, l>=0, h >=k
             hmin=0
             kmin=0

          case (143:148)       ! -3    -> hkl: h+k>0, l>0 ;  hk0: h>0, k>=0
             hmin=0
             kmin=-kmax

          case (149,151,153,157,159,162,163) ! -3 1 m  -> hkl: h>=0,h>=k>0 ; h0l: h>=0,l>=0
             hmin=0
             kmin=0
             lmin=-lmax

          case (150,152,154,155,156,158,160,161,164,165,166,167)
                              ! -3 m   -> hkl: h>=0 h>=k ; hhl: h>=0,l>=0
             hmin=0
             kmin=0

          case (168:176)    ! 6/m   -> hkl: h>0,k>0,l>=0;  0kl k>=0,l>=0
             hmin=0
             kmin=0

          case (177:194)    ! 6/mmm -> hkl: h >=0, k >=0, l >=0, h >=k
             hmin=0
             kmin=0

          case (195:206)    ! m-3   -> hkl: h > l, k > l, l >=0 ; hkk: k>=0 h>=k
             hmin=0
             kmin=0

          case (207:230)    ! m-3m  -> hkl: h >=0, k >=0, l >=0, h >=k, k >=l
             hmin=0
             kmin=0

          case default      ! Assumed -1
             hmin=-hmax
             kmin=-kmax
       end select

       num_ref=0
       ext_do: do h=hmin,hmax
          do k=kmin,kmax
             do l=lmin,lmax

                hh(1)=h
                hh(2)=k
                hh(3)=l

                if (hkl_equal(hh,nulo)) cycle
                sval=hkl_s(hh,crystalcell)
                if (sval > vmax .or. sval < vmin) cycle
                if (hkl_absent(hh,SpG)) cycle

                kk=asu_hkl(hh,SpG)
                if (hkl_equal(kk,nulo)) cycle
                if (hkl_equal(kk,-hh) .and. Friedel) cycle

                num_ref=num_ref+1
                if(num_ref > maxref) then
                   num_ref=maxref
                   exit ext_do
                end if
                hkl(:,num_ref)= kk
                mul(num_ref)  = hkl_mult(kk,SpG,friedel)
                sv(num_ref)   = sval
             end do
          end do
       end do ext_do

       if (present(no_order)) then
          if (no_order) then
             ind=(/(i,i=1,num_ref)/)
          else
             call sort(sv,num_ref,ind)
          end if
       else
          call sort(sv,num_ref,ind)
       end if

       do i=1,num_ref
          reflex(i)%h   = hkl(:,ind(i))
          reflex(i)%mult= mul(ind(i))
          reflex(i)%S   = sv(ind(i))
       end do

       return
    End Subroutine Hkl_Uni_reflect

    !!--++
    !!--++ SUBROUTINE HKL_UNI_REFLECTION
    !!--++
    !!--++    (Overloaded)
    !!--++    Calculate unique reflections between two values (value1,value2)
    !!--++    of sin_theta/lambda. If no_order is present and .true. the sort subroutine
    !!--++    is not invoked.
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Hkl_Uni_Reflection(Crystalcell,SpG,Friedel,Value1,Value2,Code,Num_Ref,Reflex,no_order)
       !---- Arguments ----!
       type (Crystal_Cell_Type),             intent(in)     :: crystalcell     ! Unit cell object
       type (Space_Group_Type) ,             intent(in)     :: SpG             ! Spacegroup object
       Logical,                              intent(in)     :: Friedel         ! If true Friedel law is applied
       real(kind=cp),                        intent(in)     :: value1,value2   ! Range in sintheta/Lambda
       character(len=1),                     intent(in)     :: code            ! If code="r", d-spacing are input
       integer,                              intent(out)    :: num_ref         ! Number of reflections generated
       type (Reflection_Type), dimension(:), intent(out)    :: reflex          ! Reflections
       logical,                   optional,  intent(in)     :: no_order

       !---- Local variables ----!
       real(kind=cp)                         :: vmin,vmax,sval
       integer                               :: h,k,l,hmin,kmin,lmin,hmax,kmax,lmax, i, maxref
       integer, dimension(3)                 :: hh,kk,nulo
       integer,  dimension(  size(reflex))   :: ind
       integer,  dimension(  size(reflex))   :: mul
       integer,  dimension(3,size(reflex))   :: hkl
       real(kind=cp),dimension(size(reflex)) :: sv
       character(len=2)                      :: inf

       nulo=0
       maxref=size(reflex)
       vmin=min(value1,value2)
       vmax=max(value1,value2)
       if (code =="r" .or. code=="R") then
          vmin=1.0/(2.0*max(value1,value2))
          vmax=1.0/(2.0*min(value1,value2))
       end if

       hmax=nint(CrystalCell%cell(1)*2.0*vmax+1.0)
       kmax=nint(CrystalCell%cell(2)*2.0*vmax+1.0)
       lmax=nint(CrystalCell%cell(3)*2.0*vmax+1.0)
       lmin= 0  !l positive or zero except for -3 1 m (see below)

       !---- Select approximate region to generate reflections depending
       !---- on the space group. This allows a faster generation.
       Select Case(SpG%NumSpg)
          case (1:2)                 ! -1    -> hkl: l >=0; hk0: h >=0; 0k0: k >=0
             hmin=-hmax
             kmin=-kmax
             if(SpG%NumSpg == 1 .and. .not. Friedel) lmin=-lmax

          case (3:15)                ! 2/m
             inf=SpG%info(1:2)
             if (inf(1:1) == "-") inf(1:1)=inf(2:2)
             select case (inf(1:1))
                case ("b")     !       -> hkl: k >=0, l >=0; hk0: h >=0
                   hmin=-hmax
                   kmin=0
                case ("c")     !       -> hkl: k >=0, l >=0; h0l: h >=0
                   hmin=-hmax
                   kmin=0
                case ("a")     !       -> hkl: h >=0, l >=0; 0kl: l >=0  Provisional (to be tested)
                   kmin=-kmax
                   hmin=0
             end select

          case (16:74)         ! mmm   -> hkl: h >=0, k >=0, l >=0
             hmin=0
             kmin=0

          case (75:88)         ! 4/m   -> hkl: h >=0, l >=0, k >=0 if h = 0
                               !                             k > 0 if h > 0
             hmin=0
             kmin=0

          case (89:142)        ! 4/mmm -> hkl: h >=0, k>=0, l>=0, h >=k
             hmin=0
             kmin=0

          case (143:148)       ! -3    -> hkl: h+k>0, l>0 ;  hk0: h>0, k>=0
             hmin=0
             kmin=-kmax

          case (149,151,153,157,159,162,163) ! -3 1 m  -> hkl: h>=0,h>=k>0 ; h0l: h>=0,l>=0
             hmin=0
             kmin=0
             lmin=-lmax

          case (150,152,154,155,156,158,160,161,164,165,166,167)
                              ! -3 m   -> hkl: h>=0 h>=k ; hhl: h>=0,l>=0
             hmin=0
             kmin=0

          case (168:176)    ! 6/m   -> hkl: h>0,k>0,l>=0;  0kl k>=0,l>=0
             hmin=0
             kmin=0

          case (177:194)    ! 6/mmm -> hkl: h >=0, k >=0, l >=0, h >=k
             hmin=0
             kmin=0

          case (195:206)    ! m-3   -> hkl: h > l, k > l, l >=0 ; hkk: k>=0 h>=k
             hmin=0
             kmin=0

          case (207:230)    ! m-3m  -> hkl: h >=0, k >=0, l >=0, h >=k, k >=l
             hmin=0
             kmin=0

          case default      ! Assumed -1
             hmin=-hmax
             kmin=-kmax
       end select

       num_ref=0
       ext_do: do h=hmin,hmax
          do k=kmin,kmax
             do l=lmin,lmax

                hh(1)=h
                hh(2)=k
                hh(3)=l

                if (hkl_equal(hh,nulo)) cycle
                sval=hkl_s(hh,crystalcell)
                if (sval > vmax .or. sval < vmin) cycle
                if (hkl_absent(hh,SpG)) cycle

                kk=asu_hkl(hh,SpG)
                if (hkl_equal(kk,nulo)) cycle
                if (hkl_equal(kk,-hh) .and. Friedel) cycle

                num_ref=num_ref+1
                if(num_ref > maxref) then
                   num_ref=maxref
                   exit ext_do
                end if
                hkl(:,num_ref)=kk
                mul(num_ref)  =hkl_mult(kk,SpG,friedel)
                sv(num_ref)   = sval
             end do
          end do
       end do ext_do

       if (present(no_order)) then
          if (no_order) then
             ind=(/(i,i=1,num_ref)/)
          else
             call sort(sv,num_ref,ind)
          end if
       else
          call sort(sv,num_ref,ind)
       end if

       do i=1,num_ref
          reflex(i)%h= hkl(:,ind(i))
          reflex(i)%mult= mul(ind(i))
          reflex(i)%S   = sv(ind(i))
       end do

       return
    End Subroutine Hkl_Uni_Reflection

    !!--++
    !!--++ SUBROUTINE HKL_UNI_REFLLIST
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate unique reflections between two values (value1,value2)
    !!--++    of sin_theta/lambda. If no_order is present and .true. the sort subroutine
    !!--++    is not invoked.
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Hkl_Uni_ReflList(Crystalcell,SpG,Friedel,Value1,Value2,Code,MaxRef,Reflex,no_order)
       !---- Arguments ----!
       type (Crystal_Cell_Type),       intent(in)     :: crystalcell      ! Unit cell object
       type (Space_Group_Type) ,       intent(in)     :: SpG              ! Spacegroup object
       Logical,                        intent(in)     :: Friedel          ! If true Friedel law is applied
       real(kind=cp),                  intent(in)     :: value1,value2    ! Range in sintheta/Lambda
       character(len=1),               intent(in)     :: code             ! If code="r", d-spacing are input
       integer,                        intent(in)     :: MaxRef           ! Maximum Number of reflections to be generated
       type (Reflection_List_Type),    intent(out)    :: reflex           ! set of reflections
       logical,             optional,  intent(in)     :: no_order

       !---- Local variables ----!
       real(kind=cp)                   :: vmin,vmax,sval
       integer                         :: h,k,l,hmin,kmin,lmin,hmax,kmax,lmax, i, num_ref
       integer, dimension(3)           :: hh,kk,nulo
       integer,  dimension(  MaxRef)   :: ind
       integer,  dimension(  MaxRef)   :: mul
       integer,  dimension(3,MaxRef)   :: hkl
       real(kind=cp),dimension(MaxRef) :: sv
       character(len=2)                :: inf

       nulo=0
       vmin=min(value1,value2)
       vmax=max(value1,value2)
       if (code =="r" .or. code=="R") then
          vmin=1.0/(2.0*max(value1,value2))
          vmax=1.0/(2.0*min(value1,value2))
       end if

       hmax=nint(CrystalCell%cell(1)*2.0*vmax+1.0)
       kmax=nint(CrystalCell%cell(2)*2.0*vmax+1.0)
       lmax=nint(CrystalCell%cell(3)*2.0*vmax+1.0)
       lmin= 0  !l positive or zero except for -3 1 m (see below)

       !---- Select approximate region to generate reflections depending
       !---- on the space group. This allows a faster generation.
       Select Case(SpG%NumSpg)
          case (1:2)                 ! -1    -> hkl: l >=0; hk0: h >=0; 0k0: k >=0
             hmin=-hmax
             kmin=-kmax
             if(SpG%NumSpg == 1 .and. .not. Friedel) lmin=-lmax

          case (3:15)                ! 2/m
             inf=SpG%info(1:2)
             if (inf(1:1) == "-") inf(1:1)=inf(2:2)
             select case (inf(1:1))
                case ("b")     !       -> hkl: k >=0, l >=0; hk0: h >=0
                   hmin=-hmax
                   kmin=0
                case ("c")     !       -> hkl: k >=0, l >=0; h0l: h >=0
                   hmin=-hmax
                   kmin=0
                case ("a")     !       -> hkl: h >=0, l >=0; 0kl: l >=0  Provisional (to be tested)
                   kmin=-kmax
                   hmin=0
             end select

          case (16:74)         ! mmm   -> hkl: h >=0, k >=0, l >=0
             hmin=0
             kmin=0

          case (75:88)         ! 4/m   -> hkl: h >=0, l >=0, k >=0 if h = 0
                               !                             k > 0 if h > 0
             hmin=0
             kmin=0

          case (89:142)        ! 4/mmm -> hkl: h >=0, k>=0, l>=0, h >=k
             hmin=0
             kmin=0

          case (143:148)       ! -3    -> hkl: h+k>0, l>0 ;  hk0: h>0, k>=0
             hmin=0
             kmin=-kmax

          case (149,151,153,157,159,162,163) ! -3 1 m  -> hkl: h>=0,h>=k>0 ; h0l: h>=0,l>=0
             hmin=0
             kmin=0
             lmin=-lmax

          case (150,152,154,155,156,158,160,161,164,165,166,167)
                              ! -3 m   -> hkl: h>=0 h>=k ; hhl: h>=0,l>=0
             hmin=0
             kmin=0

          case (168:176)    ! 6/m   -> hkl: h>0,k>0,l>=0;  0kl k>=0,l>=0
             hmin=0
             kmin=0

          case (177:194)    ! 6/mmm -> hkl: h >=0, k >=0, l >=0, h >=k
             hmin=0
             kmin=0

          case (195:206)    ! m-3   -> hkl: h > l, k > l, l >=0 ; hkk: k>=0 h>=k
             hmin=0
             kmin=0

          case (207:230)    ! m-3m  -> hkl: h >=0, k >=0, l >=0, h >=k, k >=l
             hmin=0
             kmin=0

          case default      ! Assumed -1
             hmin=-hmax
             kmin=-kmax
       end select

       num_ref=0
       ext_do: do h=hmin,hmax
          do k=kmin,kmax
             do l=lmin,lmax

                hh(1)=h
                hh(2)=k
                hh(3)=l

                if (hkl_equal(hh,nulo)) cycle
                sval=hkl_s(hh,crystalcell)
                if (sval > vmax .or. sval < vmin) cycle
                if (hkl_absent(hh,SpG)) cycle

                kk=asu_hkl(hh,SpG)
                if (hkl_equal(kk,nulo)) cycle
                if (hkl_equal(kk,-hh) .and. Friedel) cycle

                num_ref=num_ref+1
                if(num_ref > maxref) then
                   num_ref=maxref
                   exit ext_do
                end if
                hkl(:,num_ref)=kk
                mul(num_ref)  =hkl_mult(kk,SpG,friedel)
                sv(num_ref)   = sval
             end do
          end do
       end do ext_do

       if (present(no_order)) then
          if (no_order) then
             ind=(/(i,i=1,num_ref)/)
          else
             call sort(sv,num_ref,ind)
          end if
       else
          call sort(sv,num_ref,ind)
       end if

       if (allocated(reflex%ref)) deallocate(reflex%ref)
       allocate(reflex%ref(num_ref))
       reflex%nref= num_ref

       do i=1,num_ref
          reflex%Ref(i)%h    = hkl(:,ind(i))
          reflex%Ref(i)%mult = mul(ind(i))
          reflex%Ref(i)%S    = sv(ind(i))
          reflex%Ref(i)%fo    =0.0
          reflex%Ref(i)%sfo   =0.0
          reflex%Ref(i)%fc    =0.0
          reflex%Ref(i)%w     =0.0
          reflex%Ref(i)%phase =0.0
          reflex%Ref(i)%a     =0.0
          reflex%Ref(i)%b     =0.0
          reflex%Ref(i)%aa    =0.0
          reflex%Ref(i)%bb    =0.0
       end do

       return
    End Subroutine Hkl_Uni_ReflList


    !!----
    !!---- SUBROUTINE INIT_ERR_REFL
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: 20/07/2015
    !!
    Subroutine Init_Err_Refl()

       err_refl=.false.
       ERR_Refl_Mess=" "

       return
    End Subroutine Init_Err_Refl

    !!----
    !!---- SUBROUTINE INIT_REFLIST
    !!----
    !!----    Initialize the Reflection List Variable
    !!----
    !!---- Update: 20/07/2015
    !!
    Subroutine Init_RefList(Reflex,N)
       !---- Arguments ----!
       type(reflection_list_type), intent(out):: Reflex
       integer, optional,          intent(in) :: N

       !---- Local Variables ----!
       integer :: i

       !> Init
       if (allocated(reflex%ref)) deallocate(reflex%ref)

       if (present(n)) then
          reflex%nref=n
          if (n > 0) then
             allocate(reflex%ref(n))
             do i=1,n
                reflex%ref(i)%h     =0
                reflex%ref(i)%mult  =0
                reflex%ref(i)%fo    =0.0
                reflex%ref(i)%sfo   =0.0
                reflex%ref(i)%fc    =0.0
                reflex%ref(i)%w     =0.0
                reflex%ref(i)%phase =0.0
                reflex%ref(i)%a     =0.0
                reflex%ref(i)%b     =0.0
                reflex%ref(i)%aa    =0.0
                reflex%ref(i)%bb    =0.0
             end do
          end if
       else
          reflex%nref=0
       end if

       return
    End Subroutine Init_RefList


    !!--++
    !!--++ SUBROUTINE INIT_REF_COND
    !!--++
    !!--++    Initialize the Reflection conditions information array
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Init_Ref_Cond()

       Hkl_Ref_Conditions(1:20)(1:80)   = (/  &
             "(h k l)    h+k=2n : xy0 centred face (C)                                        " , &
             "(h k l)    k+l=2n : 0yz centred face (A)                                        " , &
             "(h k l)    h+l=2n : x0z centred face (B)                                        " , &
             "(h k l)  h+k+l=2n : body centred (I)                                            " , &
             "(h k l)  h,k,l same parity: all-face centred (F)                                " , &
             "(h k l) -h+k+l=3n : rhombohedrally centred (R)                                  " , &
             "(  0  k  l)     k=2n : (100) glide plane with b/2 translation (b)               " , &
             "(  0  k  l)     l=2n : (100) glide plane with c/2 translation (c)               " , &
             "(  0  k  l)   k+l=2n : (100) glide plane with b/2 + c/2 translations (n)        " , &
             "(  0  k  l)   k+l=4n : (100) glide plane with b/4 +- c/4 translations (d)       " , &
             "(  h  0  l)     h=2n : (010) glide plane with a/2 translation (a)               " , &
             "(  h  0  l)     l=2n : (010) glide plane with c/2 translation (c)               " , &
             "(  h  0  l)   l+h=2n : (010) glide plane with c/2 + a/2 translations (n)        " , &
             "(  h  0  l)   l+h=4n : (010) glide plane with c/4 +- a/4 translations (d)       " , &
             "(  h  k  0)     h=2n : (001) glide plane with a/2 translation (a)               " , &
             "(  h  k  0)     k=2n : (001) glide plane with b/2 translation (b)               " , &
             "(  h  k  0)   h+k=2n : (001) glide plane with a/2 + b/2 translations (n)        " , &
             "(  h  k  0)   h+k=4n : (001) glide plane with a/4 +- b/4 translations (d)       " , &
             "(  h  -h   0 l) l=2n : (11-20) glide plane with c/2 translation (c)             " , &
             "(  0   k  -k l) l=2n : (-2110) glide plane with c/2 translation (c)             " /)
       Hkl_Ref_Conditions(21:39)(1:80)   = (/  &
             "( -h   0   h l) l=2n : (1-210) glide plane with c/2 translation (c)             " , &
             "(  h   h -2h l) l=2n : (1-100) glide plane with c/2 translation (c)             " , &
             "(-2h   h   h l) l=2n : (01-10) glide plane with c/2 translation (c)             " , &
             "(  h -2h   h l) l=2n : (-1010) glide plane with c/2 translation (c)             " , &
             "(  h  h  l)     l=2n : (1-10) glide plane with c/2 translation (c,n)            " , &
             "(  h  k  k)     h=2n : (01-1) glide plane with a/2 translation (a,n)            " , &
             "(  h  k  h)     k=2n : (-101) glide plane with b/2 translation (b,n)            " , &
             "(  h  h  l)     l=2n : (1-10) glide plane with c/2 translation (c,n)            " , &
             "(  h  h  l)  2h+l=4n : (1-10) glide plane with a/4 +- b/4 +- c/4 translation (d)" , &
             "(  h -h  l)     l=2n : (110)  glide plane with c/2 translation (c,n)            " , &
             "(  h -h  l)  2h+l=4n : (110)  glide plane with a/4 +- b/4 +- c/4 translation (d)" , &
             "(  h  k  k)     h=2n : (01-1) glide plane with a/2 translation (a,n)            " , &
             "(  h  k  k)  2k+h=4n : (01-1) glide plane with +-a/4 + b/4 +- c/4 translation(d)" , &
             "(  h  k -k)     h=2n : (011)  glide plane with a/2 translation (a,n)            " , &
             "(  h  k -k)  2k+h=4n : (011)  glide plane with +-a/4 + b/4 +- c/4 translation(d)" , &
             "(  h  k  h)     k=2n : (-101) glide plane with b/2 translation (b,n)            " , &
             "(  h  k  h)  2h+k=4n : (-101) glide plane with +-a/4 +- b/4 + c/4 translation(d)" , &
             "( -h  k  h)     k=2n : (101)  glide plane with b/2 translation (b,n)            " , &
             "( -h  k  h)  2h+k=4n : (101)  glide plane with +-a/4 +- b/4 + c/4 translation(d)" /)
         Hkl_Ref_Conditions(40:58)(1:80)   = (/  &
             "(h 0 0)      h=2n : screw axis // [100] with  a/2 translation (21)              " , & ! monoclinic, ortho., tetra and cubic
             "(h 0 0)      h=2n : screw axis // [100] with 2a/4 translation (42)              " , & ! cubic
             "(h 0 0)      h=4n : screw axis // [100] with  a/4 translation (41)              " , & ! cubic
             "(h 0 0)      h=4n : screw axis // [100] with 3a/4 translation (43)              " , & ! cubic
             "(0 k 0)      k=2n : screw axis // [010] with  b/2 translation (21)              " , & ! monoclinic, ortho., tetra and cubic
             "(0 k 0)      k=2n : screw axis // [010] with 2b/4 translation (42)              " , & ! cubic
             "(0 k 0)      k=4n : screw axis // [010] with  b/4 translation (41)              " , & ! cubic
             "(0 k 0)      k=4n : screw axis // [010] with 3b/4 translation (43)              " , & ! cubic
             "(0 0 l)      l=2n : screw axis // [00l] with  c/2 translation (21)              " , & ! monoclinic, ortho., tetra and cubic
             "(0 0 l)      l=2n : screw axis // [00l] with 2c/4 translation (42)              " , & ! tetragonal and cubic
             "(0 0 l)      l=4n : screw axis // [00l] with  c/4 translation (41)              " , & ! tetragonal and cubic
             "(0 0 l)      l=4n : screw axis // [00l] with 3c/4 translation (43)              " , & ! tetragonal and cubic
             "(0 0 0 l)    l=2n : screw axis // [00l] axis with 3c/6 translation (63)         " , &
             "(0 0 0 l)    l=3n : screw axis // [00l] axis with  c/3 translation (31)         " , &
             "(0 0 0 l)    l=3n : screw axis // [00l] axis with 2c/3 translation (32)         " , &
             "(0 0 0 l)    l=3n : screw axis // [00l] axis with 2c/6 translation (62)         " , &
             "(0 0 0 l)    l=3n : screw axis // [00l] axis with 4c/6 translation (64)         " , &
             "(0 0 0 l)    l=6n : screw axis // [00l] axis with  c/6 translation (61)         " , &
             "(0 0 0 l)    l=6n : screw axis // [00l] axis with 5c/6 translation (65)         " /)

       return
    End Subroutine Init_Ref_Cond

    !!--++
    !!--++ SUBROUTINE INTEGRAL_CONDITIONS
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integral Conditions according with I.T. Table 2.2.13.1
    !!--++    space.
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Integral_Conditions(SpG,iunit)
       !---- Arguments ----!
       type (Space_Group_Type),  intent(in)     :: SpG  ! Space group object
       integer, optional,        intent(in)     :: iunit

       !---- local variables ----!
       integer               :: h, k,l, m
       integer               :: n, n_ext
       integer, dimension(3) :: hh
       integer               :: num_exti
       logical               :: integral_condition

       integral_condition   = .false.

       ! 1.       h+k   = 2n                   C-face centred                      C
       ! 2.       k+l   = 2n                   A-face centred                      A
       ! 3.       h+l   = 2n                   B-face centred                      B
       ! 4.       h+k+l = 2n                   Body centred                        I
       !
       ! 5.       h+k   = 2n
       !      and k+l   = 2n
       !      and h+l   = 2n                   All-face centred                    F
       !     or h,k,l all odd
       !     or h,k,l all even
       !
       ! 6.      -h+k+l = 3n                   Rhombohedrally centred,             R
       !                                     obverse setting

       if (present(iunit)) then
          write(unit=iunit,fmt=*) " "
          write(unit=iunit,fmt=*) " >>> Integral reflections conditions for centred lattices:"
          write(unit=iunit,fmt=*) "----------------------------------------------------------"
          write(unit=iunit,fmt=*) " "
       end if

       !---- C-face centred ----!
       !  Hkl_Ref_Conditions(1) =   "(h k l)  h+k=2n           : xy0 centered base"
       num_exti = 1
       n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
       n_ext = 0   ! nombre de reflecions obeissant a la regle
       do h=-6, 6
          do k=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=l
                m =  h+k
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
       end do     ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          integral_condition = .true.
       end if

       !---- A-face centred ----!
       !   Hkl_Ref_Conditions(2) =   "(h k l)  k+l=2n           : 0yz centered base"
       num_exti = 2
       n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
       n_ext = 0   ! nombre de reflecions obeissant a la regle
       do h=-6, 6
          do k=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=l
                m =  k+l
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
       end do     ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          integral_condition = .true.
       end if

       !---- B-face centred ----!
       !  Hkl_Ref_Conditions(3) =   "(h k l)  h+l=2n           : x0z centered base"
       num_exti = 3
       n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
       n_ext = 0   ! nombre de reflecions obeissant a la regle
       do h=-6, 6
          do k=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=l
                m =  h+l
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
       end do     ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          integral_condition = .true.
       end if

       !---- Body centred (I) ----!
       !  Hkl_Ref_Conditions(4) =   "(h k l)  h+k+l=2n         : body centred"
       num_exti = 4
       n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
       n_ext = 0   ! nombre de reflecions obeissant a la regle
       do h=-6, 6
          do k=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=l
                m =  h+k+l
                if (m /= int(m/2)*2) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
       end do     ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          integral_condition = .true.
       end if

       !---- all-face centred (F) ----!
       ! Hkl_Ref_Conditions(5) =   "(h k l)  h,k,l same parity: all-face centred"
       num_exti = 5
       n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
       n_ext = 0   ! nombre de reflecions obeissant a la regle
       do h=-6, 6
          do k=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=l
                if (h /= int(h/2)*2 .and.  k /= int(k/2)*2 .and. l == int(l/2)*2 ) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1

                else if(h /= int(h/2)*2 .and.  k == int(k/2)*2 .and. l /= int(l/2)*2 ) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1

                else if(h == int(h/2)*2 .and.  k /= int(k/2)*2 .and. l /= int(l/2)*2 ) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1

                else if(h == int(h/2)*2 .and.  k == int(k/2)*2 .and. l /= int(l/2)*2 ) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1

                else if(h == int(h/2)*2 .and.  k /= int(k/2)*2 .and. l == int(l/2)*2 ) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1

                else if(h /= int(h/2)*2 .and.  k == int(k/2)*2 .and. l == int(l/2)*2 ) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
       end do     ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          integral_condition = .true.
       end if

       !---- R network ----!
       !  Hkl_Ref_Conditions(6) =   "(h k l) -h+k+l=3n         : Rhombohedrally centred (R)"
       num_exti = 6
       n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
       n_ext = 0   ! nombre de reflecions obeissant a la regle
       do h=-6, 6
          do k=-6, 6
             do l=-6, 6
                hh(1)=h
                hh(2)=k
                hh(3)=l
                m =  -h+k+l
                if (m /= int(m/3)*3) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do   ! l loop
          end do    ! k loop
       end do     ! h loop
       if (n==n_ext) then
          if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
          integral_condition = .true.
       end if

       if (.not. integral_condition)   then
          if (present(iunit)) write(unit=iunit,fmt=*) "     =====>>> no general reflection condition"
       end if

       return
    End Subroutine Integral_Conditions

    !!--++
    !!--++ SUBROUTINE SCREW_AXIS_CONDITIONS
    !!--++
    !!--++    (PRIVATE)
    !!--++    Reflections conditions for Screw axes Table 2.2.13.2
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Screw_Axis_Conditions(SpG,Iunit)
       !---- Arguments ----!
       type (Space_Group_Type),       intent(in)     :: SpG   ! Space group object
       integer, optional,             intent(in)     :: iunit

       !---- Local variables ----!
       integer               :: h, k,l
       integer               :: n, n_ext
       integer, dimension(3) :: hh
       integer               :: num_exti
       logical               :: serial_condition

       serial_condition   = .false.

       if (present(iunit)) then
          write(unit=iunit,fmt=*) " "
          write(unit=iunit,fmt=*) " >>> Serial reflections conditions for screw axes:"
          write(unit=iunit,fmt=*) "---------------------------------------------------"
          write(unit=iunit,fmt=*) " "
       end if

       !SCREW AXES:      33 extinctions

       if (SpG%CrystalSys(1:10) == "Monoclinic" .or. SpG%CrystalSys(1:12) == "Orthorhombic" .or.   &
          SpG%CrystalSys(1:10) == "Tetragonal" .or. SpG%CrystalSys(1:5)  == "Cubic" ) then

          ! Hkl_Ref_Conditions(40) =   "(h 0 0)      h=2n : screw axis // [100] with  a/2 translation (21)"   ! monoclinic, ortho., tetra or cubic
          num_exti = 40
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6
             hh(1)=h
             hh(2)=0
             hh(3)=0
             if (h /= int(h/2)*2) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             serial_condition = .true.
          end if
       end if ! fin de la condition "if monoclinic, ortho, tetragonal, cubic

       if (SpG%CrystalSys(1:5) == "Cubic") then
          ! 41
          ! Hkl_Ref_Conditions(41) =   "(h 0 0)      h=2n : screw axis // [100] with  2a/4 translation (42)"   !  cubic
          num_exti = 41
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6, 1
             hh(1)=h
             hh(2)=0
             hh(3)=0
             if (h /= int(h/2)*2) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             serial_condition = .true.
          end if

          ! Hkl_Ref_Conditions(42) =   "(h 0 0)      h=4n : screw axis // [100] with  a/4 translation (41)"   ! cubic
          ! Hkl_Ref_Conditions(43) =   "(h 0 0)      h=4n : screw axis // [100] with 3a/4 translation (43)"   ! cubic
          num_exti = 42
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do h=-6, 6, 1
             hh(1)=h
             hh(2)=0
             hh(3)=0
             if (h /= int(h/4)*4) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti+1,": ", Hkl_Ref_Conditions(num_exti+1)
             serial_condition = .true.
          end if
       end if ! fin de la condition "if cubic

       if (SpG%CrystalSys(1:10) == "Monoclinic" .or. SpG%CrystalSys(1:12) == "Orthorhombic" .or.   &
          SpG%CrystalSys(1:10) == "Tetragonal" .or. SpG%CrystalSys(1:5)  == "Cubic" ) then
          ! Hkl_Ref_Conditions(44) =   "(0 k 0)      k=2n : screw axis // [010] with  b/2 translation (21)"   ! monoclinic, ortho., tetra and cubic
          num_exti = 44
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do k=-6, 6, 1
             hh(1)=0
             hh(2)=k
             hh(3)=0
             if (k /= int(k/2)*2) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             serial_condition = .true.
          end if
       end if   ! fin de la condition "if mono, ortho, tetra, cubic

       if (SpG%CrystalSys(1:5) == "Cubic") then
          ! Hkl_Ref_Conditions(45) =   "(0 k 0)      k=2n : screw axis // [010] with  2b/4 translation (42)"   ! cubic
          num_exti = 45
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do k=-6, 6, 1
             hh(1)=0
             hh(2)=k
             hh(3)=0
             if (k /= int(k/2)*2) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             serial_condition = .true.
          end if

          ! Hkl_Ref_Conditions(46) =   "(0 k 0)      k=4n : screw axis // [010] with  b/4 translation (41)"   ! cubic
          ! Hkl_Ref_Conditions(47) =   "(0 k 0)      k=4n : screw axis // [010] with 3b/4 translation (43)"   ! cubic
          num_exti = 46
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do k=-6, 6, 1
             hh(1)=0
             hh(2)=k
             hh(3)=0
             if (k /= int(k/4)*4) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti+1,": ", Hkl_Ref_Conditions(num_exti+1)
             serial_condition = .true.
          end if
       end if ! fin de la condition "if cubic

       if (SpG%CrystalSys(1:10) == "Monoclinic" .or. SpG%CrystalSys(1:12) == "Orthorhombic" .or.   &
          SpG%CrystalSys(1:5)  == "Cubic" ) then
          ! Hkl_Ref_Conditions(48) =   "(0 0 l)      l=2n : screw axis // [00l] with  c/2 translation (21)"   ! monoclinic, ortho. and cubic
          num_exti = 48
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do l=-6, 6, 1
             hh(1)=0
             hh(2)=0
             hh(3)=l
             if (l /= int(l/2)*2) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             serial_condition = .true.
          end if
       end if  ! fin de la condition mono, ortho, tetra, cubic

       if (SpG%CrystalSys(1:5) == "Cubic" .or. SpG%CrystalSys(1:10) == "Tetragonal") then
          ! 49
          ! Hkl_Ref_Conditions(49) =   "(0 0 l)      l=2n : screw axis // [00l] with  c/2 translation (21)"   ! monoclinic, ortho. and cubic
          num_exti = 49
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do l=-6, 6, 1
             hh(1)=0
             hh(2)=0
             hh(3)=l
             if (l /= int(l/2)*2) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             serial_condition = .true.
          end if

          ! 50 51
          ! Hkl_Ref_Conditions(50) =  "(0 0 l)      l=4n : screw axis // [00l] with  c/4 translation (41)"   ! tetragonal and cubic
          ! Hkl_Ref_Conditions(51) =  "(0 0 l)      l=4n : screw axis // [00l] with 3c/4 translation (43)"   ! tetragonal and cubic
          num_exti = 50
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do l=-6, 6, 1
             hh(1)=0
             hh(2)=0
             hh(3)=l
             if (l /= int(l/4)*4) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti+1,": ", Hkl_Ref_Conditions(num_exti+1)
             serial_condition = .true.
          end if
       end if ! fin de la condition "if cubic

       if (SpG%SPG_Latsy(1:1) == "h") then

          !52:
          ! Hkl_Ref_Conditions(52) =   "(0 0 0 l)    l=2n : screw axis // [00l] axis with 3c/6 translation (63)"
          num_exti = 52
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do l=-6, 6, 1
             hh(1)=0
             hh(2)=0
             hh(3)=l
             if (l /= int(l/2)*2) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop
          if (n==n_ext) then
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             serial_condition = .true.
          end if

          !53 54 55 56
          ! Hkl_Ref_Conditions(53) =   "(0 0 0 l)    l=3n : screw axis // [00l] axis with  c/3 translation (31)"
          ! Hkl_Ref_Conditions(54) =   "(0 0 0 l)    l=3n : screw axis // [00l] axis with 2c/3 translation (32)"
          ! Hkl_Ref_Conditions(55) =   "(0 0 0 l)    l=3n : screw axis // [00l] axis with 2c/6 translation (62)"
          ! Hkl_Ref_Conditions(56) =   "(0 0 0 l)    l=3n : screw axis // [00l] axis with 4c/6 translation (64)"
          num_exti = 53
          n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
          n_ext = 0   ! nombre de reflecions obeissant a la regle
          do l=-6, 6, 1
             hh(1)=0
             hh(2)=0
             hh(3)=l
             if (l /= int(l/3)*3) then
                n=n+1
                if (hkl_absent(hh,SpG)) n_ext=n_ext+1
             end if
          end do     ! h loop

          if (n==n_ext) then
                 if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
             if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti+1,": ", Hkl_Ref_Conditions(num_exti+1)
             if (SpG%laue(1:3) == "6/m") then
                  if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti+2,": ", Hkl_Ref_Conditions(num_exti+2)
                if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti+3,": ", Hkl_Ref_Conditions(num_exti+3)
             end if ! fin de la condition "6/m
             serial_condition   = .true.
          end if

          if (SpG%laue(1:3) == "6/m") then
             !57 58:
             ! Hkl_Ref_Conditions(57) =   "(0 0 0 l)    l=6n : screw axis // [00l] axis with  c/6 translation (61)"
             ! Hkl_Ref_Conditions(58) =   "(0 0 0 l)    l=6n : screw axis // [00l] axis with 5c/6 translation (65)"
             num_exti = 57
             n = 0       ! nombre de reflections pouvant obeir a la regle d"extinction
             n_ext = 0   ! nombre de reflecions obeissant a la regle
             do l=-6, 6, 1
                hh(1)=0
                hh(2)=0
                hh(3)=l
                if (l /= int(l/6)*6) then
                   n=n+1
                   if (hkl_absent(hh,SpG)) n_ext=n_ext+1
                end if
             end do     ! h loop
             if (n==n_ext) then
                  if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti,": ", Hkl_Ref_Conditions(num_exti)
                if (present(iunit)) write(unit=iunit,fmt="(tr5,a,i2,2a)")  "#",num_exti+1,": ", Hkl_Ref_Conditions(num_exti+1)
                serial_condition   = .true.
             end if
          end if ! fin de la condition "6/m
       end if  ! fin de la condition "if hexagonal

       if (.not. serial_condition)   then
          if (present(iunit)) write(unit=iunit,fmt=*) "     =====>>> no serial reflection condition"
       end if

       return
    End Subroutine Screw_Axis_Conditions

    !!----
    !!---- SUBROUTINE SEARCH_EXTINTIONS
    !!----
    !!----    Write information about the Reflections Extintion for SpaceGroup
    !!----
    !!---- Update: 20/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE SEARCH_EXTINTIONS_IUNIT
    !!--++
    !!--++    (Overloaded)
    !!--++    Write information about the Reflections Extintion for SpaceGroup
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Search_Extinctions_Iunit(SpG, Iunit)
       !---- Arguments ----!
       type (Space_Group_Type), intent(in)     :: SpG     ! Space group object
       integer,                 intent(in)     :: Iunit

       if (.not. hkl_ref_cond_ini) then
          call init_ref_cond()
          hkl_ref_cond_ini=.true.
       end if
       call integral_conditions(SpG,iunit)
       call glide_planes_conditions(SpG,iunit)
       call screw_axis_conditions(SpG,iunit)

       return
    End Subroutine Search_Extinctions_Iunit

    !!--++
    !!--++ SUBROUTINE SEARCH_EXTINCTIONS_FILE
    !!--++
    !!--++    (Overloaded)
    !!--++    Write information about the Reflections Extintion for SpaceGroup
    !!--++    in filevar variable
    !!--++
    !!--++ Update: 20/07/2015
    !!
    Subroutine Search_Extinctions_File(SpG, nlines, filevar)
       !---- Arguments ----!
       type (Space_Group_Type), intent(in)          :: SpG      ! Space group object
       integer,                 intent(out)         :: nlines   ! Number of lines
       character(len=*), dimension(:), intent(out)  :: filevar

       !---- Local Variables ----!
       integer            :: iunit,ierr
       character(len=132) :: line

       !> Init
       nlines=0
       filevar=' '

       !> Load Information
       if (.not. hkl_ref_cond_ini) then
          call init_ref_cond()
          hkl_ref_cond_ini=.true.
       end if

       open(newunit=iunit,file='search_extin_xyx.tmp')

       call integral_conditions(SpG,iunit)
       call glide_planes_conditions(SpG,iunit)
       call screw_axis_conditions(SpG,iunit)

       rewind(unit=iunit)
       do
          read(unit=iunit,fmt='(a)', iostat=ierr) line
          if (ierr /=0) exit
          nlines=nlines+1
          filevar(nlines)=trim(line)
       end do
       close(unit=iunit, status='delete')

       return
    End Subroutine Search_Extinctions_File

    !!----
    !!---- SUBROUTINE WRITE_ASU
    !!----
    !!----    Write information about the asymmetric unit for reciprocal
    !!----    space.
    !!----
    !!---- Update: 20/07/2015
    !!
    Subroutine Write_Asu(SpG, iunit)
       !---- Arguments ----!
       type (space_group_type), intent(in) :: SpG     ! Space group object
       integer,optional,        intent(in) :: iunit

       !---- Local Variables ----!
       character(len=120)                  :: line
       character(len=2)                    :: inf
       integer                             :: lun

       lun=6
       if (present(iunit)) lun=iunit

       line(1:10)=" [     ]  "
       line(3:7)=SpG%laue

       if (SpG%numspg > 0 .and. SpG%numspg <= 231) then
          select case (SpG%numspg)
             case (1:2)       ! -1
                line(11:)="hkl: l >=0    hk0: h >=0    0k0: k >=0"

             case (3:15)      ! 2/m
                inf(1:2)=adjustl(SpG%info(1:2))
                if(inf(1:1) == "-") inf(1:1)=inf(2:2)
                select case (inf(1:1))
                   case ("b")    ! 1 2/m 1
                      line(11:)="hkl: k >=0, l >=0   hk0: h >=0"
                   case ("c")    ! 1 1 2/m
                      line(11:)="hkl: k >=0, l >=0   h0l: h >=0"
                   case ("a")    ! 2/m 1 1
                      line(11:)="hkl: h >=0, l >=0   0kl: l >=0" !  Provisional (to be tested)
                end select

             case (16:74)      ! mmm
                line(11:)="hkl: h >=0, k >=0, l >=0"

             case (75:88)      ! 4/m
                line(11:)="hkl: h >=0, l >=0 with k >=0 if h =0 and k >0 if h >0"

             case (89:142)     ! 4/mmm
                line(11:)="hkl: h >=0, k >=0, l >=0 and h >=k"

             case (143:148)    ! -3
                line(11:)="hkl: h+k>0, l>0    hk0: h>0, k>=0"

             case (149,151,153,157,159,162,163)  ! -3 1 m
                line(11:)="hkl: h >=0, h >=k >0   and  h0l: h >=0, l >=0"

             case (150,152,154,155,156,158,160,161,164,165,166,167)   ! -3 m
                line(11:)="hkl: h >=0  h >=k  and   hhl: h >=0, l >=0 "

             case (168:176)  ! 6/m
                line(11:)="hkl: h > 0, k > 0, l >=0   and  0kl: k >=0, l >=0 "

             case (177:194)  ! 6/mmmm
                line(11:)="hkl: h >=0, k >=0, l >=0 with h >=k"

             case (195:206)  ! m-3
                line(11:)="hkl: h > l, k > l, l >=0  and   hkk: k >=0, h >=k"

             case (207:230)  ! m-3m
                line(11:)="hkl: h >=0, k >=0, l >=0 with h >=k  and k >=l"

          end select
       else
          select case(SpG%Laue)
             case("-1   ")
                line(11:)="hkl: l >=0    hk0: h >=0    0k0: k >=0"
             case("2/m  ")
                line(11:)="hkl: k >=0, l >=0   hk0: h >=0"
             case("mmm  ")
                line(11:)="hkl: h >=0, k >=0, l >=0"
             case("4/m  ")
                line(11:)="hkl: h >=0, l >=0 with k >=0 if h =0 and k >0 if h >0"
             case("4/mmm")
                 line(11:)="hkl: h >=0, k >=0, l >=0 and h >=k"
             case("-3   ")
                line(11:)="hkl: h+k>0, l>0    hk0: h>0, k>=0"
             case("-3m  ")
                line(11:)="hkl: h >=0  h >=k  and   hhl: h >=0, l >=0 "
             case("-31m ")
                line(11:)="hkl: h >=0, h >=k >0   and  h0l: h >=0, l >=0"
             case("6/m  ")
                line(11:)="hkl: h > 0, k > 0, l >=0   and  0kl: k >=0, l >=0 "
             case("6/mmm")
                line(11:)="hkl: h >=0, k >=0, l >=0 with h >=k"
             case("m-3  ")
                line(11:)="hkl: h > l, k > l, l >=0  and   hkk: k >=0, h >=k"
             case("m-3m ")
                line(11:)="hkl: h >=0, k >=0, l >=0 with h >=k  and k >=l"
             case default
                err_refl=.true.
                ERR_Refl_Mess=" SpG Laue Wrong"
                return
          end select
       end if

       write(unit=lun,fmt="(a)") " => Reciprocal Asymmetric Unit "
       write(unit=lun,fmt="(a)") "   "//line

       return
    End Subroutine Write_Asu

    !!----
    !!---- SUBROUTINE WRITE_REFLIST_INFO
    !!----
    !!----    Write information about the Reflection List
    !!----
    !!---- Update: 20/07/2015
    !!
    Subroutine Write_RefList_Info(Rfl, Iunit, Mode)
       !---- Arguments ----!
       type (Reflection_List_Type), intent(in) :: Rfl    ! Reflections
       integer, optional,           intent(in) :: iunit  ! Unit
       character(len=*), optional,  intent(in) :: Mode   ! Nuclear or Xray

       !---- Local variables ----!
       integer :: i,lun
       integer :: hmax,kmax,lmax,hmin,kmin,lmin
       real    :: delta

       lun=6
       if (present(iunit)) lun=iunit

       if (present(mode)) then
          Select Case (mode(1:3))
             Case("NUC","nuc")
                write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(NEUTRONS)"
                write(unit=lun,fmt="(a)")     "    ==================================================="
             Case default
                write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(X-RAYS)"
                write(unit=lun,fmt="(a)")     "    ================================================="
          End Select

       else
          write(unit=lun,fmt="(a)")   "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(X-RAYS)"
          write(unit=lun,fmt="(a)")   "    ================================================="
       end if

       hmax=maxval(rfl%ref%h(1))
       kmax=maxval(rfl%ref%h(2))
       lmax=maxval(rfl%ref%h(3))

       hmin=minval(rfl%ref%h(1))
       kmin=minval(rfl%ref%h(2))
       lmin=minval(rfl%ref%h(3))

       write(unit=lun,fmt="(/,a,/)") "   H   K   L   Mult  SinTh/Lda    |Fobs|      SFobs        |Fc|       Delta"
       do i=1,rfl%Nref
          delta=rfl%ref(i)%Fo-rfl%ref(i)%Fc
          write(unit=lun,fmt="(3i4,i5,5f12.5)") rfl%ref(i)%h, rfl%ref(i)%mult,     &
              rfl%ref(i)%S, rfl%ref(i)%Fo,rfl%ref(i)%SFo, rfl%ref(i)%Fc, delta
       end do

       write(unit=lun,fmt="(a)") " "
       write(unit=lun,fmt="(a)") " "
       write(unit=lun,fmt="(a,i6)") " => Number of Reflections: ", rfl%nref
       write(unit=lun,fmt="(a,i4,tr3,a,i4,tr3,a,i4)") " => H_max: ",hmax," K_max: ",kmax," L_max: ",lmax
       write(unit=lun,fmt="(a,i4,tr3,a,i4,tr3,a,i4)") " => H_min: ",hmin," K_min: ",kmin," L_min: ",lmin
       write(unit=lun,fmt="(a)") " "

       return
    End Subroutine Write_RefList_Info

 End Module CFML_Reflections_Utilities


