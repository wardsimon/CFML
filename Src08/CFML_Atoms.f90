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
    Use CFML_Strings,      only: u_case  
    Use CFML_gSpaceGroups, only: spg_type, apply_op           

    !---- Variables ----!
    implicit none

    private

    !---- List of public procedures ----!

    !---- Parameters ----!
    real(kind=cp), parameter :: R_ATOM=1.1           ! Average atomic radius
    
    !---- Types ----!
    
    !!----
    !!---- TYPE :: ATM_TYPE
    !!----
    !!----      Simple Atom type
    !!----
    Type, Public :: Atm_Type
       character(len=20)             :: Lab     =" "     ! Label for atom
       character(len=2)              :: ChemSymb=" "     ! Chemical symbol
       integer                       :: Z       = 0      ! Atomic number
       integer                       :: Mult    = 0      ! Multiplicity
       real(kind=cp), dimension(3)   :: X       = 0.0_cp ! Coordinates
       real(kind=cp)                 :: Occ     = 1.0_cp ! Occupancy factor
       character(len=4)              :: UType   ="beta"  ! Options: U, B, beta
       character(len=3)              :: ThType  ="iso"   ! Options: iso, ani
       real(kind=cp)                 :: U_iso   = 0.0_cp ! Uiso or Ueq
       real(kind=cp), dimension(6)   :: U       = 0.0_cp ! Anisotropic thermal factors
       logical                       :: Magnetic=.false.
       real(kind=cp), dimension(3)   :: Moment  = 0.0_cp ! Magnetic moment
       character(len=4)              :: SfacSymb=" "     ! SFac symbol
       integer, dimension(3)         :: Ind_ff  = 0      ! Index of form factor (Xray, b, Magff)
    End Type Atm_Type 
    
    !!----
    !!---- TYPE :: ATM_STD_TYPE
    !!----
    !!----      Simple Atom type + Standard deviation values
    !!----
    Type, Public, Extends(Atm_Type) :: Atm_Std_Type
       real(kind=cp), dimension(3)  :: X_Std      = 0.0_cp     ! standard deviations
       real(kind=cp)                :: Occ_Std    = 0.0_cp     
       real(kind=cp)                :: U_iso_Std  = 0.0_cp     
       real(kind=cp)                :: U_Std      = 0.0_cp   
       real(kind=cp), dimension(3)  :: Moment_std = 0.0_cp  
    End Type Atm_Std_Type 
    
    !!----
    !!---- TYPE :: MATM_STD_TYPE
    !!----
    !!----      
    !!----
    Type, Public, Extends(Atm_Std_Type) :: MAtm_Std_Type
       character(len=3)                         :: wyck   = " "
       integer                                  :: n_mc   = 0       ! up to 8
       integer                                  :: n_dc   = 0
       real(kind=cp),dimension(6,8)             :: Mcs    = 0.0_cp  ! Mcos,Msin up to 8  (Mcx Mcy  Mcz , Msx  Msy  Msz)
       real(kind=cp),dimension(6,8)             :: Mcs_std= 0.0_cp  ! Mcos,Msin up to 8  (Mcx Mcy  Mcz , Msx  Msy  Msz)
       real(kind=cp),dimension(6,8)             :: Dcs    = 0.0_cp  ! Dcos,Dsin up to 8  (Dcx Dcy  Dcz , Dsx  Dsy  Dsz)
       real(kind=cp),dimension(6,8)             :: Dcs_std= 0.0_cp  ! Dcos,Dsin up to 8  (Dcx Dcy  Dcz , Dsx  Dsy  Dsz)
    End Type MAtm_Std_Type
    
    !!----
    !!---- TYPE :: ATM_REF_TYPE
    !!----
    !!----      Refinement Atom type
    !!----
    Type, Public, Extends(Atm_Std_Type) :: Atm_Ref_Type
       integer,      dimension(3)               :: LX       =0      ! Code for parameters
       integer                                  :: LOcc     =0
       integer                                  :: LU_iso   =0 
       integer,      dimension(6)               :: LU       =0
       real(kind=cp),dimension(3)               :: MX       =0.0_cp ! Factor of refinement
       real(kind=cp)                            :: MOcc     =0.0_cp
       real(kind=cp)                            :: MU_iso   =0.0_cp
       real(kind=cp),dimension(6)               :: MU       =0.0_cp
    End Type Atm_Ref_Type 
    
    !!----
    !!---- TYPE :: MATM_REF_TYPE
    !!----
    !!----      Refinement Atom type
    !!----
    Type, Public, Extends(MAtm_Std_Type) :: MAtm_Ref_Type
       integer,      dimension(3)               :: LX       =0      ! Code for parameters
       integer                                  :: LOcc     =0
       integer                                  :: LU_iso   =0 
       integer,      dimension(6)               :: LU       =0
       real(kind=cp),dimension(3)               :: MX       =0.0_cp ! Factor of refinement
       real(kind=cp)                            :: MOcc     =0.0_cp
       real(kind=cp)                            :: MU_iso   =0.0_cp
       real(kind=cp),dimension(6)               :: MU       =0.0_cp
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
       integer,              dimension( :,:), allocatable :: neighb_atom    ! Ptr.->neighbour (# in list)(nat,idp)
       real(kind=cp),        dimension( :,:), allocatable :: distance       ! Corresponding distances (nat,idp)
       real(kind=cp),      dimension(:, :,:), allocatable :: trans          ! Lattice translations   (3,nat,idp)
       integer                                            :: ndist          ! Number of distinct distances
       real(kind=cp),           dimension(:), allocatable :: ddist          ! List of distinct distances(nat*idp)
       character (len=20),      dimension(:), allocatable :: ddlab          ! Labels of atoms at ddist (nat*idp)
    End Type Atm_Cell_Type
    
    !!----
    !!---- TYPE :: ALIST_TYPE
    !!--..
    !!
    Type, Public :: AtList_Type
       integer                                    :: natoms=0   ! Number of atoms in the list
       logical,         dimension(:), allocatable :: Active     ! Flag for active or not
       class(Atm_Type), dimension(:), allocatable :: Atom       ! Atoms
    End type AtList_Type
    
    !---- Interface Zone ----!
    Interface
       Module Subroutine Init_Atom_Type(Atm)
          !---- Arguments ----!
          class(Atm_Type), intent(in out)   :: Atm
       End Subroutine Init_Atom_Type
       
       Module Subroutine Allocate_Atom_List(N, A)
          !---- Arguments ----!
          integer,             intent(in)       :: n    
          class(atlist_type),  intent(in out)   :: A    
       End Subroutine Allocate_Atom_List
       
       Module Subroutine Read_Bin_Atom_List(filename, A)
          !---- Arguments ----!
          character(len=*),   intent(in)    :: filename
          class(atlist_type), intent(in out) :: A
       End Subroutine Read_Bin_Atom_List  
       
       Module Subroutine Write_Bin_Atom_List(filename, A)
          !---- Arguments ----!
          character(len=*),   intent(in) :: filename
          class(atlist_type), intent(in) :: A 
       End Subroutine Write_Bin_Atom_List   
       
       Module Subroutine Write_Info_Atom_List(A, Iunit)
          !---- Arguments ----!
          class(atlist_type),              intent(in) :: A        
          integer,              optional, intent(in) :: IUnit    
       End Subroutine Write_Info_Atom_List
       
       Module Subroutine Extend_List(A, B, Spg, Conven)
          !---- Arguments ----!
          class(atlist_type),   intent(in)     :: A         
          class(atlist_type),   intent(in out) :: B         
          type(SpG_Type),      intent(in)     :: SpG       
          logical, optional,   intent(in)     :: Conven    
       End Subroutine Extend_List
          
    End Interface
    
 Contains
 
End Module CFML_Atoms
