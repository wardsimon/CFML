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
!!---- MODULE: CFML_Atoms
!!----   INFO: Subroutines related to Atoms definitions
!!----
!!---- HISTORY
!!----    Update: 06/03/2011
!!----
!!
 Module CFML_Atoms

    !---- Use Modules ----!
    Use CFML_GlobalDeps
    Use CFML_Maths,        only: modulo_lat, equal_vector
    Use CFML_Metrics,      only: Cell_G_Type
    Use CFML_Strings,      only: u_case,l_case
    Use CFML_gSpaceGroups, only: spg_type, apply_op, SuperSpaceGroup_Type

    !---- Variables ----!
    implicit none

    private

    !---- List of public procedures ----!
    public :: Allocate_Atom_List, Extend_Atom_List, Init_Atom_Type, Read_Bin_Atom_List, &
              Write_Bin_atom_List, Write_Atom_List
    public :: Equiv_Atm, Wrt_Lab


    !---- Parameters ----!
    real(kind=cp), parameter :: R_ATOM=1.1_cp      ! Average atomic radius

    integer,public, parameter :: max_mod=8

    !---- Types ----!

    !!----
    !!---- TYPE :: ATM_TYPE
    !!----
    !!----  Simple Atom type containing the minimum set of model parameters
    !!----  to describe the atom within a periodic crystal structure. The
    !!----  simplest atom is characterized by its labels (Lab,ChemSymb,SfacSymb, ThType),
    !!----  position X, isotropic displacement parameter U_iso, occupancy factor Occ,
    !!----  magnetic moment, anisotropic displacement parameters U, charge, atomic number,
    !!----  maximum module of magnetic moment and additional information.
    !!----
    Type, Public :: Atm_Type
       character(len=20)             :: Lab     =" "     ! Label (identifier) of the atom
       character(len=2)              :: ChemSymb=" "     ! Chemical symbol
       character(len=4)              :: SfacSymb=" "     ! Scattering Factor Symbol
       integer                       :: Z       = 0      ! Atomic number
       integer                       :: Mult    = 0      ! Multiplicity of the Wyckoff position
       integer                       :: Charge  = 0      ! Charge, ionic state
       real(kind=cp), dimension(3)   :: X       = 0.0_cp ! Fractional Coordinates
       real(kind=cp)                 :: U_iso   = 0.0_cp ! Biso, Uiso or Ueq (if iso =Biso)
       real(kind=cp)                 :: Occ     = 1.0_cp ! Occupancy factor
       character(len=4)              :: UType   ="beta"  ! Options: U, B, beta
       character(len=3)              :: ThType  ="iso"   ! Options: iso, ani
       real(kind=cp), dimension(6)   :: U       = 0.0_cp ! Anisotropic thermal factors
       logical                       :: Magnetic=.false. ! Flag indication if the atom is magnetic or not.
       real(kind=cp)                 :: Mom     = 0.0_cp ! Maximum Module of Magnetic moment
       real(kind=cp), dimension(3)   :: Moment  = 0.0_cp ! Magnetic moment
       integer, dimension(3)         :: Ind_ff  = 0      ! Index of form factor (Xray, b, Magff)
       character(len=40)             :: AtmInfo = " "    ! Information string for different purposes
       character(len=5)              :: wyck    = " "    ! Wyckoff position label if known
    End Type Atm_Type

    !!----
    !!---- TYPE :: ATM_STD_TYPE
    !!----
    !!----  Simple Atom type + Standard deviation values
    !!----
    Type, Public, Extends(Atm_Type) :: Atm_Std_Type
       real(kind=cp), dimension(3)  :: X_Std      = 0.0_cp     ! standard deviations
       real(kind=cp)                :: Occ_Std    = 0.0_cp
       real(kind=cp)                :: U_iso_Std  = 0.0_cp
       real(kind=cp), dimension(6)  :: U_Std      = 0.0_cp
       real(kind=cp), dimension(3)  :: Moment_std = 0.0_cp
    End Type Atm_Std_Type

    !!----
    !!---- TYPE :: MATM_STD_TYPE
    !!----
    !!---- This type Modulated Atom type extends Atm_Std_Type by adding modulation
    !!---- Cosine(c) and Sine amplitudes(s) to each model parameter characterizing
    !!---- normal atoms. Up to max_mod harmonic numbers (Q_coeffs) are allowed
    !!----
    Type, Public, Extends(Atm_Std_Type)    :: MAtm_Std_Type
       integer                             :: n_oc   = 0       ! Number of occupation amplitudes
       integer                             :: n_mc   = 0       ! Number of moment amplitudes
       integer                             :: n_dc   = 0       ! Number of static displacement amplitudes
       integer                             :: n_uc   = 0       ! Number of thermal displacement amplitudes
       integer,      dimension(max_mod)    :: poc_q  = 0       ! Pointer to Q_coeffs of occupatiom amplitudes
       integer,      dimension(max_mod)    :: pmc_q  = 0       ! Pointer to Q_coeffs of moment amplitudes
       integer,      dimension(max_mod)    :: pdc_q  = 0       ! Pointer to Q_coeffs of displacement amplitudes
       integer,      dimension(max_mod)    :: puc_q  = 0       ! Pointer to Q_coeffs of thermal displacement amplitudes
       real(kind=cp),dimension(2, max_mod) :: Ocs    = 0.0_cp  ! Ocos,Osin up to 8  (Oc, Os)
       real(kind=cp),dimension(2, max_mod) :: Ocs_std= 0.0_cp  !
       real(kind=cp),dimension(6, max_mod) :: Mcs    = 0.0_cp  ! Mcos,Msin up to 8  (Mcx Mcy  Mcz , Msx  Msy  Msz)
       real(kind=cp),dimension(6, max_mod) :: Mcs_std= 0.0_cp  !
       real(kind=cp),dimension(6, max_mod) :: Dcs    = 0.0_cp  ! Dcos,Dsin up to 8  (Dcx Dcy  Dcz , Dsx  Dsy  Dsz)
       real(kind=cp),dimension(6, max_mod) :: Dcs_std= 0.0_cp  !
       real(kind=cp),dimension(12,max_mod) :: Ucs    = 0.0_cp  ! Ucos,Usin up to 8  (Dcx Dcy  Dcz , Dsx  Dsy  Dsz)
       real(kind=cp),dimension(12,max_mod) :: Ucs_std= 0.0_cp  !
       real(kind=cp),dimension(:),   allocatable :: Xs      ! Position in superspace
       real(kind=cp),dimension(:),   allocatable :: Moms    ! Moment in superspace
       real(kind=cp),dimension(:,:), allocatable :: Us      ! Thermal factos in superspace
    End Type MAtm_Std_Type

    !!----
    !!---- TYPE :: ATM_REF_TYPE
    !!----
    !!----   Refinement Atom type: This type extends Atm_Std_Type with the numbers
    !!----   of each model paramenter within the list of LSQ free parameters, as well
    !!----   as the corresponding multipliers for constraints.
    !!----
    Type, Public, Extends(Atm_Std_Type) :: Atm_Ref_Type
       integer,      dimension(3)               :: L_X       =0      ! Code number of free parameters
       integer                                  :: L_Occ     =0
       integer                                  :: L_U_iso   =0
       integer,      dimension(3)               :: L_moment  =0
       integer,      dimension(6)               :: L_U       =0
       real(kind=cp),dimension(3)               :: M_X       =0.0_cp ! Multipliers
       real(kind=cp)                            :: M_Occ     =0.0_cp
       real(kind=cp)                            :: M_U_iso   =0.0_cp
       real(kind=cp),dimension(3)               :: M_moment  =0.0_cp
       real(kind=cp),dimension(6)               :: M_U       =0.0_cp
    End Type Atm_Ref_Type

    !!----
    !!---- TYPE :: MATM_REF_TYPE
    !!----
    !!----   Refinement Atom type: This type extends MAtm_Std_Type with the numbers
    !!----   of each model paramenter within the list of LSQ free parameters, as well
    !!----   as the corresponding multipliers for constraints.
    !!----
    Type, Public, Extends(MAtm_Std_Type) :: MAtm_Ref_Type
       integer,      dimension(3)               :: L_X       =0       ! Code Numbers of parameter
       integer                                  :: L_Occ     =0
       integer                                  :: L_U_iso   =0
       integer,      dimension(6)               :: L_U       =0
       real(kind=cp),dimension(3)               :: M_X       =0.0_cp  ! Multipliers of refinement
       real(kind=cp)                            :: M_Occ     =0.0_cp
       real(kind=cp)                            :: M_U_iso   =0.0_cp
       real(kind=cp),dimension(6)               :: M_U       =0.0_cp
       integer,      dimension(2, max_mod)      :: L_Ocs    = 0       ! Code Numbers of parameter
       integer,      dimension(6, max_mod)      :: L_Mcs    = 0       !
       integer,      dimension(6, max_mod)      :: L_Dcs    = 0       !
       integer,      dimension(12,max_mod)      :: L_Ucs    = 0       !
       real(kind=cp),dimension(2, max_mod)      :: M_Ocs    = 0.0_cp  ! Multipliers
       real(kind=cp),dimension(6, max_mod)      :: M_Mcs    = 0.0_cp  !
       real(kind=cp),dimension(6, max_mod)      :: M_Dcs    = 0.0_cp  !
       real(kind=cp),dimension(12,max_mod)      :: M_Ucs    = 0.0_cp  !
    End Type MAtm_Ref_Type


    !!----
    !!---- TYPE ::ATM_CELL_TYPE
    !!--..
    !!---- This type is mostly used for distance-angle and Bond-valence calculations.
    !!---- It holds the position and coordination of all the atoms in the conventional
    !!---- unit cell as well as their distances to neighbours atoms.
    !!----
    !!---- 13/06/2019
    !!
    Type, public :: Atm_Cell_Type
       integer                                            :: nat            ! Total number of atoms
       character(len=20),       dimension(:), allocatable :: Lab            ! Labels for atoms (dimension Nat)
       real(kind=cp),         dimension(:,:), allocatable :: xyz            ! Fractional coordinates (3,nat)
       real(kind=cp),           dimension(:), allocatable :: charge
       real(kind=cp),           dimension(:), allocatable :: moment
       real(kind=cp),         dimension(:,:), allocatable :: var_free       ! Free variables (10,nat)
       integer,                 dimension(:), allocatable :: neighb         ! Number of neighbours (nat)
       integer,               dimension(:,:), allocatable :: neighb_atom    ! Ptr.->neighbour (# in list)(nat,idp)
       real(kind=cp),         dimension(:,:), allocatable :: distance       ! Corresponding distances (nat,idp)
       real(kind=cp),       dimension(:,:,:), allocatable :: trans          ! Lattice translations   (3,nat,idp)
       integer                                            :: ndist          ! Number of distinct distances
       real(kind=cp),           dimension(:), allocatable :: ddist          ! List of distinct distances(nat*idp)
       character (len=20),      dimension(:), allocatable :: ddlab          ! Labels of atoms at ddist (nat*idp)
    End Type Atm_Cell_Type

    !!---- Type, Public :: Atom_Equiv_Type
    !!----    integer                                        :: mult
    !!----    character(len=2)                               :: ChemSymb
    !!----    character(len=10),allocatable, dimension(:)    :: Lab
    !!----    real(kind=sp),    allocatable, dimension(:,:)  :: x
    !!---- End Type Atom_Equiv_Type
    !!----
    !!----  Updated: January 2014
    !!
    Type, Public :: Atom_Equiv_Type
       integer                                        :: mult
       character(len=2)                               :: ChemSymb
       character(len=20),allocatable, dimension(:)    :: Lab
       real(kind=cp),    allocatable, dimension(:,:)  :: x
    End Type Atom_Equiv_Type

    !!---- Type, Public :: Atom_Equiv_List_Type
    !!----    integer                                           :: nauas
    !!----    type (Atom_Equiv_Type), allocatable, dimension(:) :: atm
    !!---- End Type Atom_Equiv_List_Type
    !!----
    !!----  Updated: January 2014
    !!
    Type, Public :: Atom_Equiv_List_Type
       integer                                           :: nauas
       type (Atom_Equiv_Type), allocatable, dimension(:) :: atm
    End Type Atom_Equiv_List_Type

    !!----
    !!---- TYPE :: ALIST_TYPE
    !!--..
    !!
    Type, Public :: AtList_Type
       integer                                    :: natoms=0        ! Number of atoms in the list
       character(len=9)                           :: mcomp="Crystal" ! For magnetic moments and modulation functions Mcs and Dcs It may be also "Cartesian" or "Spherical"
       logical                                    :: symm_checked=.false.
       logical,         dimension(:), allocatable :: Active          ! Flag for active or not
       class(Atm_Type), dimension(:), allocatable :: Atom            ! Atoms
    End type AtList_Type

    !Overload

    Interface Extend_Atom_List
      Module Procedure Extend_List              !Creating a new AtList_Type with all atoms in unit cell
      Module Procedure Set_Atom_Equiv_List      !Creating a an Atom_Equiv_List_Type from AtList_Type in asymmetric unit
    End Interface Extend_Atom_List


    !---- Interface Zone ----!
    Interface

       Pure Module Function Equiv_Atm(Nam1,Nam2,NameAt) Result(Equiv_Atom)
          !---- Arguments ----!
          character (len=*), intent (in) :: nam1,nam2
          character (len=*), intent (in) :: NameAt
          logical                        :: equiv_atom
       End Function Equiv_Atm

       Pure Module Function Wrt_Lab(Nam1,Nam2) Result(Bilabel)
          !---- Arguments ----!
          character (len=*), intent (in) :: nam1,nam2
          character (len=8)              :: bilabel
       End Function Wrt_Lab

       Module Subroutine Init_Atom_Type(Atm,d)
          !---- Arguments ----!
          class(Atm_Type), intent(in out)   :: Atm
          integer,         intent(in)       :: d     ! Number of k-vectors
       End Subroutine Init_Atom_Type

       Module Subroutine Allocate_Atom_List(N, A,Type_Atm, d)
          !---- Arguments ----!
          integer,             intent(in)       :: n    ! Atoms in the List
          type(Atlist_type),   intent(in out)   :: A    ! Objet to be allocated
          character(len=*),    intent(in)       :: Type_Atm !Atomic type: Atm, Atm_Std, MAtm_Std, Atm_Ref, MAtm_Ref
          integer,             intent(in)       :: d    ! Number of k-vectors
       End Subroutine Allocate_Atom_List

       Module Subroutine Read_Bin_Atom_List(filename, A, Type_Atm)
          !---- Arguments ----!
          character(len=*),   intent(in)     :: filename
          type(atlist_type),  intent(in out) :: A
          character(len=*),   intent(in)     :: Type_Atm
       End Subroutine Read_Bin_Atom_List

       Module Subroutine Write_Bin_Atom_List(filename, A)
          !---- Arguments ----!
          character(len=*),   intent(in) :: filename
          type(atlist_type),  intent(in) :: A
       End Subroutine Write_Bin_Atom_List

       Module Subroutine Write_Atom_List(A, Iunit, SpG)
          !---- Arguments ----!
          type(atlist_type),                   intent(in) :: A        ! Atom list object
          integer, optional,                   intent(in) :: IUnit    ! Logical unit
          type(SuperSpaceGroup_type),optional, intent(in) :: SpG
       End Subroutine Write_Atom_List

       Module Subroutine Extend_List(A, B, Spg, Type_Atm,Conven)
          !---- Arguments ----!
          type(atlist_type),   intent(in)     :: A         ! Atom list (asymmetric unit)
          type(atlist_type),   intent(in out) :: B         ! Atom list into the unit cell
          type(SpG_Type),      intent(in)     :: SpG       ! SpaceGroup
          character(len=*),    intent(in)     :: Type_Atm  ! !Atomic type: Atm, Atm_Std, MAtm_Std, Atm_Ref, MAtm_Ref
          logical, optional,   intent(in)     :: Conven    ! If present and .true. using the whole conventional unit cell
       End Subroutine Extend_List

       Module Subroutine Set_Atom_Equiv_List(SpG,cell,A,Ate,lun)
         type(SpG_Type),             intent(in) :: SpG
         type(Cell_G_Type),          intent(in) :: Cell
         type(Atlist_Type),          intent(in) :: A
         type(Atom_Equiv_List_Type), intent(out):: Ate
         integer, optional,          intent(in) :: lun
       End Subroutine Set_Atom_Equiv_List

    End Interface

 End Module CFML_Atoms
