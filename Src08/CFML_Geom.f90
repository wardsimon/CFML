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
!!----          Nebil Ayape Katcho      (ILL)
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
!!---- MODULE: CFML_Geometry_Calc
!!----   INFO: Routines for Geometry Calculations
!!----
!!---- HISTORY
!!----    Update: 06/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!--++    CFML_Math_3D:  Cross_Product
!!--++    CFML_GlobalDeps: Eps, Pi, Cp, Sp, To_Rad, To_Deg
!!--++    CFML_Math_General: Acosd, Cosd, Sind
!!--++    CFML_Crystal_Metrics: Cell_G_Type
!!----
!!---- VARIABLES
!!----    COORDINATION_TYPE
!!----    COORD_INFO
!!--++    EPSI
!!----    ERR_GEOM
!!----    ERR_GEOM_MESS
!!----    POINT_LIST_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       ANGLE_DIHEDRAL
!!--++       ANGLE_DIHEDRAL_IJKN            [Overloaded]
!!--++       ANGLE_DIHEDRAL_UVW             [Overloaded]
!!----       ANGLE_MOD
!!--++       ANGLE_MODN                     [Overloaded]
!!--++       ANGLE_MODV                     [Overloaded]
!!----       ANGLE_UV
!!--++       ANGLE_UVI                      [Overloaded]
!!--++       ANGLE_UVR                      [Overloaded]
!!----       COORD_MOD
!!--++       COORD_MODN                     [Overloaded]
!!--++       COORD_MODV                     [Overloaded]
!!----       DISTANCE
!!--++       DISTANCE_FR                    [Overloaded]
!!--++       DISTANCE_FR_DP                 [Overloaded]
!!--++       DISTANCE_SC                    [Overloaded]
!!----       MATRIX_PHITHECHI
!!----       MATRIX_RX
!!----       MATRIX_RY
!!----       MATRIX_RZ
!!----
!!----    Subroutines:
!!----       ALLOCATE_COORDINATION_TYPE
!!----       ALLOCATE_POINT_LIST
!!----       ANGLE_AND_SIGMA
!!----       CALC_DIST_ANGLE
!!----       CALC_DIST_ANGLE_SIGMA
!!----       DEALLOCATE_COORDINATION_TYPE
!!----       DEALLOCATE_POINT_LIST
!!----       DISTANCE_AND_SIGMA
!!----       GET_ANGLEN_AXIS_FROM_ROTMAT
!!----       GET_EULER_FROM_FRACT
!!----       GET_MATRIX_MOVING_V_TO_U
!!----       GET_OMEGACHIPHI
!!----       GET_PHITHECHI
!!----       GET_TRANSF_LIST
!!----       INIT_ERR_GEOM
!!----       P1_DIST
!!----       PRINT_DISTANCES
!!----       SET_NEW_ASYMUNIT
!!----       SET_ORBITS_INLIST
!!----       SET_ROTATION_MATRIX
!!----       SET_TDIST_COORDINATION
!!----       SET_TDIST_PARTIAL_COORDINATION
!!----       TORSION_AND_SIGMA
!!----
!!
 Module CFML_Geom

    !---- Use Modules ----!
    use CFML_GlobalDeps
    use CFML_Rational
    use CFML_Maths,         only: Modulo_Lat, Cross_Product, Inverse_Matrix, Determ3D
    use CFML_Trigonometry,  only: acosd, cosd, sind
    use CFML_Strings,       only: Frac_Trans_1Dig, L_Case,U_Case,pack_string,String_NumStd
    use CFML_Metrics,       only: Cell_G_Type, Get_Deriv_Orth_Cell,Rot_Gibbs_Matrix
    use CFML_Atoms,         only: AtList_Type,Atm_Cell_Type,Equiv_Atm, Wrt_Lab, Atom_Equiv_List_Type, &
                                  Allocate_Atom_List, Atm_Std_Type
    use CFML_gSpaceGroups,  only: SpG_Type, Apply_OP, Get_Multip_Pos, Is_Lattice_Vec, &
                                  searchop, Write_SymTrans_Code, Get_Orbit

    implicit none

    private

    !---- List of public functions ----!

    !---- List of public overloaded procedures: functions ----!
    public :: Angle_Dihedral, Angle_Mod, Angle_Uv, Distance, Matrix_PhiTheChi, Matrix_Rx, &
              Matrix_Ry, Matrix_Rz

    !---- List of public subroutines ----!
    public :: Allocate_Coordination_Type, Allocate_Point_List, Calc_Dist_Angle, Calc_Dist_Angle_Sigma, &
              Deallocate_Coordination_Type, Deallocate_Point_List, Distance_and_Sigma, Get_Euler_From_Fract, &
              Get_PhiTheChi, P1_Dist, Print_Distances, Set_Orbits_InList, Set_TDist_Coordination, &
              Get_Transf_List, Set_TDist_Partial_Coordination, Get_Anglen_Axis_From_RotMat, Get_Matrix_moving_v_to_u, &
              Get_OmegaChiPhi, Set_Rotation_Matrix, Set_New_AsymUnit,Angle_and_Sigma, Torsion_and_Sigma

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!
    private :: Angle_Dihedral_Uvw,  Angle_Dihedral_Ijkn, Angle_Uvi, Angle_Uvr, Angle_Modn, Angle_Modv, &
               Distance_fr, Distance_fr_dp, Distance_sc

    !---- Definitions ----!

    !!----
    !!---- TYPE :: COORDINATION_TYPE
    !!--..
    !!---- Type, public :: Coordination_Type
    !!----    integer                                      :: Natoms    ! number of atoms
    !!----    integer                                      :: Max_Coor  ! Maximum number of connected atoms to a given one
    !!----    integer,       dimension(:),     allocatable :: Coord_Num ! Counter of distances connected to the current atom
    !!----    integer,       dimension(:,:),   allocatable :: N_Cooatm  ! Pointer to the ordinal number in the list of the attached
    !!----                                                              ! atom to the atom given by the first index
    !!----    integer,       dimension(:,:),   allocatable :: N_Sym     !
    !!----    real(kind=cp), dimension(:,:),   allocatable :: Dist      ! List of distances related to an atom
    !!----    real(kind=cp), dimension(:,:),   allocatable :: S_Dist    ! List of Sigma(distances)
    !!----    real(kind=cp), dimension(:,:,:), allocatable :: Tr_coo    !
    !!---- End type Coordination_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Coordination_Type
       integer                                      :: Natoms    ! number of atoms
       integer                                      :: Max_Coor  ! Maximum number of connected atoms to a given one
       integer,       dimension(:),     allocatable :: Coord_Num ! Counter of distances connected to the current atom
       integer,       dimension(:,:),   allocatable :: N_Cooatm  ! Pointer to the ordinal number in the list of the attached
                                                                 ! atom to the atom given by the first index
       integer,       dimension(:,:),   allocatable :: N_Sym     ! Number of symmetry operator to apply to N_Cooatm
       real(kind=cp), dimension(:,:),   allocatable :: Dist      ! List of distances related to an atom
       real(kind=cp), dimension(:,:),   allocatable :: S_Dist    ! List of Sigma(distances)
       real(kind=cp), dimension(:,:,:), allocatable :: Tr_coo
    End type Coordination_Type

    !!----
    !!---- COORD_INFO
    !!----    type(Coordination_Type), public :: coord_info
    !!----
    !!----    Coordination Information
    !!----
    !!---- Update: March - 2005
    !!
    type(Coordination_Type), public :: coord_info

    !!--++
    !!--++ EPSI
    !!--++    real(kind=cp), parameter :: epsi=0.001
    !!--++
    !!--++    (PRIVATE)
    !!--++    Epsilon for roughly comparing distances
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), parameter, private :: epsi=0.001


    !!----
    !!---- TYPE :: POINT_LIST_TYPE
    !!--..
    !!---- Type, public :: Point_List_Type
    !!----    integer                                       :: np   !number of points in list
    !!----    character(len=20), dimension(:),  allocatable :: nam  !name/label associated to each point
    !!----    integer,           dimension(:),  allocatable :: p    !integer pointer for various purposes
    !!----    real(kind=cp)      dimension(:,:),allocatable :: x    !fractional coordinates of points
    !!---- End type Point_List_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: point_list_type
       integer                                       :: np   !number of points in list
       character(len=20), dimension(:),  allocatable :: nam  !name/label associated to each point
       integer,           dimension(:),  allocatable :: p    !integer pointer for various purposes
       real(kind=cp),     dimension(:,:),allocatable :: x    !fractional coordinates of points
    End type point_list_type


    !---- Interfaces - Overlapp ----!
    Interface  Angle_Dihedral
       Module Procedure Angle_Dihedral_Ijkn
       Module Procedure Angle_Dihedral_Uvw
    End Interface

    Interface  Angle_Uv
       Module Procedure Angle_UvI
       Module Procedure Angle_UvR
    End Interface

    Interface  Angle_Mod
       Module Procedure Angle_ModN
       Module Procedure Angle_ModV
    End Interface

    Interface  Distance
       Module Procedure Distance_FR_DP
       Module Procedure Distance_FR
       Module Procedure Distance_SC
    End Interface

    !------------------------!
    !---- Interface Zone ----!
    !------------------------!
    Interface

      Pure Module Function Angle_Dihedral_Ijkn(ri,rj,rk,rn) result(angle)
         !---- Arguments ----!
         real(kind=cp), dimension(3), intent( in) :: ri,rj,rk,rn
         real(kind=cp)                            :: angle
      End Function Angle_Dihedral_Ijkn

      Pure Module Function Angle_Dihedral_Uvw(u,v,w) result(angle)
         !---- Argument ----!
         real(kind=cp), dimension(3), intent( in) :: u,v,w
         real(kind=cp)                            :: angle
      End Function Angle_Dihedral_Uvw

      Pure Module Function Angle_ModN(Angle) Result(AngMod)
         !---- Arguments ----!
         real(kind=cp), intent(in) :: Angle
         real(kind=cp)             :: AngMod
      End Function Angle_ModN

      Pure Module Function Angle_ModV(V_Angle) Result(VAngMod)
         !---- Arguments ----!
         real(kind=cp), dimension(:),intent(in) :: V_Angle
         real(kind=cp), dimension(size(V_Angle)):: VAngMod
      End Function Angle_ModV

      Pure Module Function Angle_UvI(Ui,Vi,G) Result(Angle)
         !---- Argument ----!
         integer, dimension(:),   intent( in)                 :: ui
         integer, dimension(:),   intent( in)                 :: vi
         real(kind=cp), dimension(:,:), intent( in), optional :: g   !metric tensor
         real(kind=cp)                                        :: angle
      End Function Angle_uvi

      Pure Module Function Angle_UvR(u,v,g) result(angle)
         !---- Argument ----!
         real(kind=cp), dimension(:),   intent( in)           :: u
         real(kind=cp), dimension(:),   intent( in)           :: v
         real(kind=cp), dimension(:,:), intent( in), optional :: g   !metric tensor
         real(kind=cp)                                        :: angle
      End Function Angle_uvr

      Pure Module Function Distance_Fr(X0,X1,Celda) Result(Dis)
         !---- Arguments ----!
         real(kind=cp), dimension(3), intent(in) :: x0,x1
         type (Cell_G_Type),          intent(in) :: Celda
         real(kind=cp)                           :: dis
      End Function Distance_Fr

      Pure Module Function Distance_Fr_dp(X0,X1,Celda) Result(Dis)
         !---- Arguments ----!
         real(kind=dp), dimension(3), intent(in) :: x0,x1
         type (Cell_G_Type),    intent(in) :: Celda
         real(kind=dp)                           :: dis
      End Function Distance_Fr_dp

      Pure Module Function Distance_SC(X0,X1,Code) Result(Dis)
         !---- Arguments ----!
         real(kind=cp), dimension(3), intent(in) :: x0,x1
         character(len=*), optional,  intent(in) :: Code
         real(kind=cp)                           :: dis
      End Function Distance_SC

      Pure Module Function Matrix_Phithechi(Phi,Theta,Chi,Code) Result(Mt)
         !---- Arguments ----!
         real(kind=cp),                intent(in) :: Phi
         real(kind=cp),                intent(in) :: Theta
         real(kind=cp),                intent(in) :: Chi
         character(len=*), optional,   intent(in) :: Code
         real(kind=cp), dimension(3,3)            :: Mt
      End Function Matrix_Phithechi

      Pure Module Function Matrix_Rx(Ang,Code) Result(Mt)
         !---- Arguments ----!
         real(kind=cp),               intent(in) :: Ang
         character(len=*), optional,  intent(in) :: Code
         real(kind=cp), dimension(3,3)           :: Mt
      End Function Matrix_Rx

      Pure Module Function Matrix_Ry(Ang,Code) Result(Mt)
         !---- Arguments ----!
         real(kind=cp),               intent(in) :: Ang
         character(len=*), optional,  intent(in) :: Code
         real(kind=cp), dimension(3,3)           :: Mt
      End Function Matrix_Ry

      Pure Module Function Matrix_Rz(Ang,Code) Result(Mt)
         !---- Arguments ----!
         real(kind=cp),               intent(in) :: Ang
         character(len=*), optional,  intent(in) :: Code
         real(kind=cp), dimension(3,3)           :: Mt
      End Function Matrix_Rz

      Pure Module Function Set_Rotation_Matrix(ang) Result(Rot)
        real(kind=cp), dimension(3),   intent( in) :: ang
        real(kind=cp), dimension(3,3)              :: Rot
      End Function Set_Rotation_Matrix

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

      Module Subroutine Allocate_Coordination_Type(nasu,numops,dmax,Max_Coor)
         !---- Arguments ----!
         integer,       intent(in) :: nasu
         integer,       intent(in) :: numops
         real(kind=cp), intent(in) :: dmax
         integer,      intent(out) :: Max_Coor
      End Subroutine Allocate_Coordination_Type

      Module Subroutine Allocate_Point_List(n,Pl,Ier)
         !---- Arguments ----!
         integer,               intent(in)     :: n
         type(point_list_type), intent(in out) :: pl
         integer,               intent(out)    :: ier
      End subroutine Allocate_Point_List

      Module Subroutine Angle_and_Sigma(Cellp,DerM,x1,x0,x2,s1,s0,s2,ang,s)
         !---- Arguments ----!
         Type(Cell_G_Type),         intent(in)  :: Cellp         ! Cell object
         real(kind=cp), dimension(3,3,6), intent(in)  :: DerM          ! Matrix of derivatives of Cellp%Cr_Orth_cel
         real(kind=cp), dimension(3),     intent(in)  :: x0,x1,x2      ! Three points in fractional coordinates and sigmas, X0 is central
         real(kind=cp), dimension(3),     intent(in)  :: s0,s1,s2      ! Sigmas of the three points
         real(kind=cp),                   intent(out) :: ang,s         ! Angle and sigma
      End Subroutine Angle_and_Sigma

      Module Subroutine Calc_Dist_Angle(Dmax, Dangl, Cell, Spg, A, Lun)
         !---- Arguments ----!
         real(kind=cp),         intent(in)   :: Dmax, Dangl
         type (Cell_G_Type),    intent(in)   :: Cell
         Class(SpG_Type),       intent(in)   :: SpG
         type (AtList_Type),    intent(in)   :: A
         integer, optional,     intent(in)   :: lun
      End Subroutine Calc_Dist_Angle

      Module Subroutine Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spg, A, Lun, Lun_cons, Lun_cif,filrest,rdmax,ramin)
         !---- Arguments ----!
         real(kind=cp),             intent(in)   :: dmax, dangl
         type (Cell_G_Type),        intent(in)   :: Cell
         Class(SpG_Type),           intent(in)   :: SpG
         type (AtList_Type),        intent(in)   :: A
         integer, optional,         intent(in)   :: lun
         integer, optional,         intent(in)   :: lun_cons
         integer, optional,         intent(in)   :: lun_cif
         character(len=*), optional,intent(in)   :: filrest
         real(kind=cp),    optional,intent(in)   :: rdmax, ramin
      End Subroutine Calc_Dist_Angle_Sigma

      Module Subroutine Deallocate_Coordination_Type()
      End Subroutine Deallocate_Coordination_Type

      Module Subroutine Deallocate_Point_List(Pl)
         !---- Arguments ----!
         type(point_list_type), intent(in out) :: pl
      End Subroutine Deallocate_Point_List

      Module Subroutine Distance_and_Sigma(Cellp,DerM,x0,x1,s0,s1,dis,s)
         !---- Arguments ----!
         Type(Cell_G_Type),         intent(in)  :: Cellp         ! Cell object
         real(kind=cp), dimension(3,3,6), intent(in)  :: DerM          ! Matrix of derivatives of Cellp%Cr_Orth_cel
         real(kind=cp), dimension(3),     intent(in)  :: x0,x1,s0,s1   ! Two points in fractional coordinates and sigmas
         real(kind=cp),                   intent(out) :: dis,s         ! Distance and sigma
      End Subroutine Distance_and_Sigma

      Module Subroutine Get_Anglen_Axis_From_RotMat(R,axis,angle)
        Real(kind=cp), dimension(3,3), intent(in) :: R
        Real(kind=cp), dimension(3),   intent(out):: axis
        Real(kind=cp),                 intent(out):: angle
      End Subroutine Get_Anglen_Axis_From_RotMat

      Module Subroutine Get_Euler_From_Fract(X1,X2,X3,Mt,Phi,Theta,Chi,Eum,Code)
         !---- Arguments ----!
         real(kind=cp),           dimension(3),   intent (in) :: x1,x2,x3
         real(kind=cp),           dimension(3,3), intent (in) :: Mt
         real(kind=cp),                           intent(out) :: theta,phi,chi
         real(kind=cp), optional, dimension(3,3), intent(out) :: EuM
         character(len=*), optional,              intent (in) :: Code
      End Subroutine Get_Euler_From_Fract

      Module Subroutine Get_Matrix_moving_v_to_u(v,u,R,w,ang)
        real(kind=cp), dimension(3),           intent(in)  :: v,u
        real(kind=cp), dimension(3,3),         intent(out) :: R
        real(kind=cp), optional,               intent(out) :: ang
        real(kind=cp), optional,dimension(3),  intent(out) :: w
      End Subroutine Get_Matrix_moving_v_to_u

      Module Subroutine Get_OmegaChiPhi(Mt,Omega,Chi,Phi,Code)  !Conventional Euler angles of diffractometry
         !---- Arguments ----!
         real(kind=cp), dimension(3,3),intent(in)  :: Mt
         real(kind=cp),                intent(out) :: Omega
         real(kind=cp),                intent(out) :: Chi
         real(kind=cp),                intent(out) :: Phi
         character(len=*), optional,   intent(in)  :: Code
      End Subroutine Get_OmegaChiPhi

      Module Subroutine Get_PhiTheChi(Mt,Phi,Theta,Chi,Code)
         !---- Arguments ----!
         real(kind=cp), dimension(3,3),intent(in)  :: Mt
         real(kind=cp),                intent(out) :: Phi
         real(kind=cp),                intent(out) :: Theta
         real(kind=cp),                intent(out) :: Chi
         character(len=*), optional,   intent(in)  :: Code
      End Subroutine Get_PhiTheChi

      Module Subroutine Get_Transf_List(trans,ox,pl,npl)
         !---- Arguments ----!
         real(kind=cp),         dimension(3,3), intent(in)     :: trans
         real(kind=cp),         dimension(3  ), intent(in)     :: ox
         type(point_list_type),                 intent(in)     :: pl
         type(point_list_type),                 intent(in out) :: npl
      End Subroutine Get_Transf_List

      Module Subroutine P1_Dist(Dmax, Cell, Spg, Ac, Lun)
         !---- Arguments ----!
         real(kind=cp),         intent(in)       :: dmax
         type (Cell_G_Type),    intent(in)       :: Cell
         Class(SpG_Type),       intent(in)       :: SpG
         type (Atm_Cell_Type),  intent(in out)   :: Ac
         integer, optional,     intent(in)       :: lun
      End Subroutine P1_Dist

      Module Subroutine Print_Distances(Lun, Dmax, Cell, Spg, A)
         !-- Arguments --!
         integer,               intent(in)   :: lun
         real(kind=cp),         intent(in)   :: dmax
         type (Cell_G_Type),    intent(in)   :: Cell
         Class(SpG_Type),       intent(in)   :: SpG
         type (AtList_Type),    intent(in)   :: A
      End Subroutine Print_Distances

      Module Subroutine Set_New_AsymUnit(SpGn,Ate,Mat,orig,A_n,matkind,debug)
         Class(SpG_Type) ,              intent(in    ) :: SpGn
         type (Atom_Equiv_List_Type),   intent(in    ) :: Ate !In old group
         real(kind=cp), dimension (3,3),intent(in    ) :: Mat
         real(kind=cp), dimension (  3),intent(in    ) :: orig
         type (AtList_Type),            intent(out   ) :: A_n
         character (len=*), optional,   intent(in    ) :: matkind
         character (len=*), optional,   intent(in    ) :: debug
      End Subroutine Set_New_AsymUnit


      Module Subroutine Set_Orbits_Inlist(Spg,Pl)
         !---- Arguments ----!
         Class(SpG_Type),        intent(in)     :: SpG
         type(point_list_type),  intent(in out) :: pl
      End Subroutine Set_Orbits_Inlist

      Module Subroutine Set_TDist_Coordination(max_coor,Dmax, Cell, Spg, A)
         !---- Arguments ----!
         integer,                  intent(in)   :: max_coor
         real(kind=cp),            intent(in)   :: dmax
         type (cell_G_Type),       intent(in)   :: Cell
         Class(SpG_Type),          intent(in)   :: SpG
         type (AtList_Type),       intent(in)   :: A
      End Subroutine Set_TDist_Coordination

      Module Subroutine Set_TDist_Partial_Coordination(List,max_coor,Dmax, Cell, Spg, A)
         !---- Arguments ----!
         integer,             intent(in)   :: List
         integer,             intent(in)   :: max_coor
         real(kind=cp),       intent(in)   :: dmax
         type (Cell_G_Type),  intent(in)   :: Cell
         Class(SpG_Type),     intent(in)   :: SpG
         type (AtList_Type),  intent(in)   :: A
      End Subroutine Set_TDist_Partial_Coordination

      Module Subroutine Torsion_and_Sigma(Cellp, x1,x2,x3,x4,sx1,sx2,sx3,sx4,tor,s)
         !---- Arguments ----!
         Type(Cell_G_Type),         intent(in)  :: Cellp         ! Cell object
         real(kind=cp), dimension(3),     intent(in)  :: x1,x2,x3,x4       ! Three points in fractional coordinates and sigmas, X0 is central
         real(kind=cp), dimension(3),     intent(in)  :: sx1,sx2,sx3,sx4   ! Sigmas of the three points
         real(kind=cp),                   intent(out) :: tor,s             ! Angle and sigma
      End Subroutine Torsion_and_Sigma

    End interface

 End Module CFML_Geom
