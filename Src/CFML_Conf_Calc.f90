!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2015  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_BVS_Energy_Calc
!!----   INFO: Subroutines related to calculations of energy or
!!----         configuration properties depending on the crystal structure: BVS, Energy,....
!!----
!!---- HISTORY
!!----    Updated: 16/12/2014
!!----
!!----
!!---- DEPENDENCIES
!!----
!!----
!!---- VARIABLES
!!----    ATOM_CONF_LIST
!!----    BVEL_ANIONS_N
!!----    BVEL_ANIONS
!!----    BVEL_ANIONS_RION
!!----    BVEL_SPECIES_N
!!----    BVEL_TABLE
!!----    BVS_ANIONS_N
!!----    BVS_ANIONS
!!----    BVS_ANIONS_RION
!!----    BVS_PAR_TYPE
!!----    BVS_SPECIES_N
!!----    BVS_TABLE
!!----    ERR_CONF
!!----    ERR_CONF_MESS
!!----    SSBVS_SPECIES_N
!!----    SBVS_TABLE
!!--++    TABLE_Alpha
!!--++    TABLE_Avcoor
!!--++    TABLE_B
!!--++    TABLE_D0
!!--++    TABLE_Dzero
!!--++    TABLE_Rcutoff
!!--++    TABLE_Ref
!!--++    TABLE_Rmin
!!--++    TABLE_Rzero
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ALLOCATE_ATOMS_CONF_LIST
!!--++       BOND_VALENCE                  [Private]
!!----       CALC_BVS
!!----       CALC_MAP_BVS
!!--++       COMPLETE_TABLE                [Private]
!!----       COST_BVS
!!----       COST_BVS_COULOMBREP
!!----       DEALLOCATE_ATOMS_CONF_LIST
!!----       DEALLOCATE_BVS_TABLE
!!----       INIT_ERR_CONF
!!----       SET_BVEL_TABLE
!!----       SET_BVS_TABLE
!!----       SET_SBVS_TABLE
!!----       SET_TABLE_D0_B
!!----       SPECIES_ON_LIST
!!----
!!
 Module CFML_BVS_Energy_Calc
    !---- Use Files ----!
    Use CFML_GlobalDeps,                 only: Sp,Cp
    Use CFML_Math_General,               only: Sort_Strings,Sort
    use CFML_String_Utilities,           only: Getword, U_Case,pack_string, get_logunit
    Use CFML_Scattering_Chemical_Tables, only: Get_Ionic_Radius, Get_Chemsymb, Get_Covalent_Radius
    use CFML_Crystal_Metrics,            only: Crystal_Cell_Type
    use CFML_Crystallographic_Symmetry,  only: Space_Group_Type
    use CFML_Atom_TypeDef,               only: Atom_type, Init_Atom_type, Write_Atom_List, Atom_list_Type, Allocate_Atom_List, &
                                               Deallocate_Atom_List, AtList1_ExtenCell_AtList2
    use CFML_Geometry_Calc,              only: Coord_Info, Distance
    use CFML_Export_VTK,                 only: write_grid_VESTA

    !---- Variables ----!
    implicit none

    private

    !---- List of public subroutines ----!
    public :: Allocate_Atoms_Conf_List, Calc_BVS, Deallocate_Atoms_Conf_List,    &
              Init_Err_Conf, Set_BVS_Table, Set_Table_d0_b, Species_on_List,     &
              Calc_Map_BVS, Cost_BVS, Cost_BVS_CoulombRep, Deallocate_BVS_Table, &
              Set_BVEL_Table, Deallocate_BVEL_Table, Set_Table_BVEL_Params,      &
              Set_SBVS_Table

    !---- List of public private ----!
    private :: Bond_Valence, Complete_Table, Complete_Table_BVEL

    !---- Definitions ----!

    !!----
    !!---- TYPE :: ATOMS_CONF_LIST_TYPE
    !!--..
    !!---- Type, public :: Atoms_Conf_List_Type
    !!----    integer                                     :: natoms    ! Total number of atoms in the list
    !!----    integer                                     :: N_Spec    ! Number of different species in the list
    !!----    integer                                     :: N_Anions  ! Number of anions in the list
    !!----    integer                                     :: N_Cations ! Number of cations in the list
    !!----    real(kind=cp)                               :: Tol       ! Tolerance(%) for sum of radii conditions
    !!----    real(kind=cp)                               :: totatoms  ! Total number of atoms in the unit cell
    !!----    character(len=4), dimension(:), allocatable :: Species   ! Symbol + valence
    !!----    real(kind=cp),    dimension(:), allocatable :: Radius    !ionic/atomic radius of species
    !!----    type(Atom_Type),  dimension(:), allocatable :: atom
    !!---- End Type Atoms_Conf_List_Type
    !!----
    !!---- Update: March - 2005
    !!
    Type, public :: Atoms_Conf_List_Type
       integer                                     :: natoms    ! Total number of atoms in the list
       integer                                     :: N_Spec    ! Number of different species in the list
       integer                                     :: N_Anions  ! Number of anions in the list
       integer                                     :: N_Cations ! Number of cations in the list
       real(kind=cp)                               :: Tol       ! Tolerance(%) for sum of radii conditions
       real(kind=cp)                               :: totatoms  ! Total number of atoms in the unit cell
       character(len=4), dimension(:), allocatable :: Species
       real(kind=cp),    dimension(:), allocatable :: Radius    !ionic/atomic radius of species
       type(Atom_Type),  dimension(:), allocatable :: Atom
    End type Atoms_Conf_List_Type

    !!----
    !!---- BVEL_ANIONS_N
    !!----    integer, parameter, public :: bvel_anions_n
    !!----
    !!----    Number of anions known in BVEL Table in: Stefan Adams and R. Prasada Rao
    !!----
    !!---- Created: December - 2014
    !!
    integer, parameter, public :: bvel_anions_n=1

    !!----
    !!---- BVEL_Anions
    !!----    character(len=4), parameter, dimension(bvel_anions_n) :: bvel_anions
    !!----
    !!----    Anions known from Stefan Adams and R. Prasada Rao
    !!----
    !!---- Created: December - 2014
    !!
    character(len=*), parameter, dimension(bvel_anions_n), public :: bvel_anions = &
                     (/"O-2 "/)

    !!----
    !!---- BVEL_Anions_Rion
    !!----    real(kind=cp), parameter, dimension(bvs_anions_n) :: bvs_anions_rion
    !!----
    !!----    Radii Ionic for Anions in BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp), parameter, dimension(bvel_anions_n), public :: bvel_anions_rion = &
                      (/1.40/)

    !!----
    !!---- TYPE :: BVEL_PAR_TYPE
    !!--..
    !!---- Type, public :: Bvel_Par_Type
    !!---   character(len=5)                       :: symb     !Symbol of the cation
    !!---   real(kind=cp),dimension(bvel_anions_n) :: Avcoor   !Average cation coordination number
    !!---   real(kind=cp),dimension(bvel_anions_n) :: Rzero    !Modified Bond-Valence parameter R0
    !!---   real(kind=cp),dimension(bvel_anions_n) :: Rcutoff  !Cutoff distance in Angstroms
    !!---   real(kind=cp),dimension(bvel_anions_n) :: Dzero    !First Morse potential parameter (eV)
    !!---   real(kind=cp),dimension(bvel_anions_n) :: Rmin     !Second Morse potential parameter (Angstroms)
    !!---   real(kind=cp),dimension(bvel_anions_n) :: alpha    !Third Morse potential parameter (1/b) (Angstroms^-1)
    !!---   integer      ,dimension(bvel_anions_n) :: refnum   !Pointer to reference paper
    !!---- End Type Bvel_Par_Type
    !!----
    !!----    Type Definition for BVEL Parameters
    !!----
    !!---- Created: December - 2014
    !!
    Type, public :: Bvel_Par_Type
       character(len=5)                       :: symb
       real(kind=cp),dimension(bvel_anions_n) :: Avcoor
       real(kind=cp),dimension(bvel_anions_n) :: Rzero
       real(kind=cp),dimension(bvel_anions_n) :: Rcutoff
       real(kind=cp),dimension(bvel_anions_n) :: Dzero
       real(kind=cp),dimension(bvel_anions_n) :: Rmin
       real(kind=cp),dimension(bvel_anions_n) :: alpha
       integer      ,dimension(bvel_anions_n) :: refnum
    End Type Bvel_Par_Type

    !!----
    !!---- BVEL_SPECIES_N
    !!----    integer, parameter, public :: bvel_species_n
    !!----
    !!----    Maximum Number of species in BVEL_Table
    !!----
    !!---- Created: December - 2014
    !!
    integer, parameter, public :: bvel_species_n=132

    !!----
    !!---- BVEL_TABLE
    !!----    Type(Bvel_Par_Type), allocatable, dimension(:), public :: BVEL_Table
    !!----
    !!----    BVEL Parameters for calculations
    !!----
    !!---- Created: December - 2014
    !!
    Type(Bvel_Par_Type), allocatable, dimension(:), public :: BVEL_Table


    !!----
    !!---- BVS_ANIONS_N
    !!----    integer, parameter, public :: bvs_anions_n
    !!----
    !!----    Number of anions known in BVS Table in O"Keefe, Breese, Brown
    !!----
    !!---- Update: March - 2005
    !!
    integer, parameter, public :: bvs_anions_n=14

    !!----
    !!---- BVS_Anions
    !!----    character(len=4), parameter, dimension(bvs_anions_n) :: bvs_anions
    !!----
    !!----    Anions known from O'Keefe, Bresse, Brown
    !!----
    !!---- Update: March - 2005
    !!
    character(len=*), parameter, dimension(bvs_anions_n), public :: bvs_anions = &
                     (/"O-2 ","F-1 ","CL-1","BR-1","I-1 ","S-2 ","SE-2","TE-2",  &
                       "N-3 ","P-3 ","AS-3","H-1 ","O-1 ","SE-1"/)

    !!----
    !!---- BVS_Anions_Rion
    !!----    real(kind=cp), parameter, dimension(bvs_anions_n) :: bvs_anions_rion
    !!----
    !!----    Ionic Radii for Anions
    !!----
    !!---- Update: March - 2005
    !!
    real(kind=cp), parameter, dimension(bvs_anions_n), public :: bvs_anions_rion = &
                      (/1.40,1.19,1.67,1.95,2.16,1.84,1.98,2.21,1.71,2.12,2.22,2.08,1.35,1.80/)



    !!---- TYPE :: BVS_PAR_TYPE
    !!--..
    !!---- Type, public :: Bvs_Par_Type
    !!----    character (len=4)                      :: Symb      ! Chemical symbol
    !!----    real(kind=cp), dimension(bvs_anions_n) :: D0        ! D0 Parameter
    !!----    real(kind=cp), dimension(bvs_anions_n) :: B_Par     ! B Parameter
    !!----    integer,       dimension(bvs_anions_n) :: refnum    ! Integer pointing to the reference paper
    !!---- End Type Bvs_Par_Type
    !!----
    !!----    Definition for BVS Parameters
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Bvs_Par_Type
       character (len=4)                     :: Symb
       real(kind=cp),dimension(bvs_anions_n) :: d0
       real(kind=cp),dimension(bvs_anions_n) :: b_par
       integer      ,dimension(bvs_anions_n) :: refnum
    End Type Bvs_Par_Type

    !!----
    !!---- BVS_SPECIES_N
    !!----    integer, parameter, public :: bvs_species_n
    !!----
    !!----    Maximum Number of species in BVS_Table
    !!----
    !!---- Update: March - 2005
    !!
    integer, parameter, public :: bvs_species_n=247

    !!----
    !!---- BVS_TABLE
    !!----    Type(Bvs_Par_Type), allocatable, dimension(:), public :: BVS_Table
    !!----
    !!----    BVS Parameters for calculations
    !!----
    !!---- Update: March - 2005
    !!
    Type(Bvs_Par_Type), allocatable, dimension(:), public :: BVS_Table

    !!----
    !!---- ERR_CONF
    !!----    logical, public  :: err_conf
    !!----
    !!----    Logical Variable taking the value .true. if an error in the module
    !!----    CONFIGURATIONS_CALCULATIONS occurs.
    !!----
    !!---- Update: March - 2005
    !!
    logical, public  :: err_conf

    !!----
    !!---- ERR_CONF_MESS
    !!----    character(len=150), public:: ERR_Conf_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: March - 2005
    !!
    character(len=150), public :: ERR_Conf_Mess

    !!----
    !!---- SBVS_SPECIES_N
    !!----    integer, parameter, public :: sbvs_species_n=5 !only alkali chalcogenides
    !!----
    !!----    Maximum Number of species in SBVS_Table
    !!----
    !!---- Update: December - 2014
    !!
    integer, parameter, public :: sbvs_species_n=5

    !!----
    !!---- SBVS_TABLE
    !!----    Type(Bvs_Par_Type), allocatable, dimension(:), public :: sBVS_Table
    !!----
    !!----    SBVS Parameters for calculations (only alkali chalcogenides are available)
    !!----
    !!---- Created: December - 2014
    !!
    Type(Bvs_Par_Type), allocatable, dimension(:), public :: sBVS_Table

    !!----
    !!---- Table_Alpha
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Alpha
    !!----
    !!----    Matrix N_Species x N_Species of Alpha (equivalent to 1/b in BVS) parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, private :: Table_Alpha

    !!----
    !!---- Table_Avcoor
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Avcoor
    !!----
    !!----    Matrix N_Species x N_Species of Average coordination parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, private :: Table_Avcoor

    !!----
    !!---- TABLE_B
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_b
    !!----
    !!----    Matrix N_Species x N_Species of B parameters for BVS
    !!----
    !!---- Update: March - 2005
    !!
    real(kind=cp),dimension(:,:), allocatable, public :: Table_b

    !!----
    !!---- TABLE_D0
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_d0
    !!----
    !!----    Matrix N_Species x N_Species of D0 for BVS
    !!----
    !!---- Update: March - 2005
    !!
    real(kind=cp),dimension(:,:), allocatable, public :: Table_d0

    !!----
    !!---- Table_Dzero
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Dzero
    !!----
    !!----    Matrix N_Species x N_Species of Dzero parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, private :: Table_Dzero

    !!----
    !!---- Table_Rcutoff
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rcutoff
    !!----
    !!----    Matrix N_Species x N_Species of Rcutoff parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rcutoff

    !!----
    !!---- TABLE_Ref
    !!----    integer,dimension(:,:), allocatable, private :: Table_ref
    !!----
    !!----    Matrix N_Species x N_Species with references for BVS parameters
    !!----
    !!---- Update: March - 2005
    !!
    integer,dimension(:,:), allocatable, private :: Table_ref

    !!----
    !!---- Table_Rmin
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rmin
    !!----
    !!----    Matrix N_Species x N_Species of Rmin parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rmin

     !!----
    !!---- Table_Rzero
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rzero
    !!----
    !!----    Matrix N_Species x N_Species of Rzero (equivalent to D0 in BVS) parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rzero

    !!----
    !!----  Reference list for BVS parameters
    !!----
    !!----
    !!----    List of Reference for BVS Data
    !!----
    !!---- Update: March - 2005
    !!
    character(len=*),dimension(0:34),parameter, private :: references = (/  &
         "Unknown                                                                         ", &
         "Brown and Altermatt, (1985), Acta Cryst. B41, 244-247 (empirical)               ", &
         "Brese and O'Keeffe, (1991), Acta Cryst. B47, 192-197 (extrapolated)             ", &
         "Adams, 2001, Acta Cryst. B57, 278-287 (includes second neighbours)              ", &
         "Hu et al. (1995) Inorg. Chim. Acta, 232, 161-165.                               ", &
         "I.D.Brown Private communication                                                 ", &
         "Brown et al. (1984) Inorg. Chem. 23, 4506-4508                                  ", &
         "Palenik (1997) Inorg. Chem. 36 4888-4890                                        ", &
         "Kanowitz and Palenik (1998) Inorg. Chem. 37 2086-2088                           ", &
         "Wood and Palenik (1998) Inorg. Chem. 37 4149-4151                               ", &
         "Liu and Thorp (1993) Inorg. Chem. 32 4102-4105                                  ", &
         "Palenik (1997) Inorg. Chem. 36 3394-3397                                        ", &
         "Shields, Raithby, Allen and Motherwell (1999) Acta Cryst.B56, 455-465           ", &
         "Chen, Zhou and Hu (2002) Chinese Sci. Bul. 47, 978-980.                         ", &
         "Kihlbourg (1963) Ark. Kemi 21 471; Schroeder 1975 Acta Cryst. B31, 2294         ", &
         "Allmann (1975) Monatshefte Chem. 106, 779                                       ", &
         "Zachariesen (1978) J.Less Common Metals 62, 1                                   ", &
         "Krivovichev and Brown (2001) Z. Krist. 216, 245                                 ", &
         "Burns, Ewing and Hawthorne (1997) Can. Miner. 35,1551-1570                      ", &
         "Garcia-Rodriguez, et al. (2000) Acta Cryst. B56, 565-569                        ", &
         "Mahapatra et al. (1996) J. Amer.Chem. Soc. 118, 11555                           ", &
         "Wood and Palenik (1999) Inorg. Chem. 38, 1031-1034                              ", &
         "Wood and Palenik (1999) Inorg. Chem. 38, 3926-3930                              ", &
         "Wood, Abboud, Palenik and Palenik (2000) Inorg. Chem. 39, 2065-2068             ", &
         "Tytko, Mehnike and Kurad (1999) Structure and Bonding 93, 1-66                  ", &
         "Gundemann, et al.(1999) J. Phys. Chem. A 103, 4752-4754                         ", &
         "Zocchi (2000) Solid State Sci. 2 383-387                                        ", &
         "Jensen, Palenik and Tiekiak (2001) Polyhedron 20, 2137                          ", &
         "Roulhac and Palenik (2002) Inorg. Chem. 42, 118-121                             ", &
         "Holsa et al.(2002) J.Solid State Chem 165, 48-55                                ", &
         "Trzesowska, Kruszynski & Bartezak (2004) Acta Cryst. B60, 174-178               ", &
         "Locock & Burns (2004) Z.Krist. 219, 267-271                                     ", &
         "J.Rodriguez-Carvajal Private communication                                      ", &
         "S. Adams and R. Prasada Rao,Phys. Status Solidi A 208, No. 8, 1746–1753 (2011)  ", &
         "S. Adams, Acta Crystallographica B57, 278-287 (2001)                            "/)

 Contains

    !!----
    !!---- Subroutine Allocate_Atoms_Conf_List(N,A)
    !!----    integer, intent(in)                         :: n    !  In -> Atoms in asymmetric unit
    !!----    type (Atoms_Conf_List_Type), intent(in out) :: A    !  In -> Objet to be allocated
    !!----
    !!----    Allocation of objet A of type atom_list_Conf. This subroutine
    !!----    should be called before using an object of type atom_list_Conf.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_Atoms_Conf_List(n,A)
       !---- Arguments ----!
       integer,                     intent(in)       :: N  !N. atoms in asymmetric unit
       type (Atoms_Conf_List_Type), intent(in out)   :: A  !Objet to be allocated

       !---- Local variables ----!
       integer :: i

       A%natoms   = n
       A%n_spec   = 0
       A%n_anions = 0
       A%n_cations= 0
       A%tol      = 20.0 ! 20% tolerance by default
       A%totatoms = 0.0
       if(.not. allocated(A%atom)) allocate (A%atom(n))
       do i=1,n
          call init_atom_type(A%atom(i))
       end do

       return
    End Subroutine Allocate_Atoms_Conf_List

    !!--++
    !!--++ Subroutine Bond_Valence(d0,b0,d,sd,bv,sbv)
    !!--++    real(kind=cp), intent(in)  :: d0    ! BVS parameter
    !!--++    real(kind=cp), intent(in)  :: b0    ! BVS parameter
    !!--++    real(kind=cp), intent(in)  :: d     ! Distance
    !!--++    real(kind=cp), intent(in)  :: sd    ! Sigma(d)
    !!--++    real(kind=cp), intent(out) :: bv    ! BVS
    !!--++    real(kind=cp), intent(out) :: sbv   ! Sigma(bv)
    !!--++
    !!--++    (Private)
    !!--++    Zachariasen exponential expression of Bond Valence
    !!--++
    !!--++ Update: Dic - 2007
    !!
    Subroutine Bond_Valence(D0,B0,D,Sd,Bv,Sbv)
       !---- Arguments ----!
       real(kind=cp),  intent(in)  :: d0,b0  !Bond-valence parameters
       real(kind=cp),  intent(in)  :: d,sd   !Bond distance and sigma
       real(kind=cp),  intent(out) :: bv,sbv !Bond-valence and sigma

       bv=EXP((d0-d)/b0)
       sbv=bv*sd/b0

       return
    End Subroutine Bond_Valence

    !!----
    !!---- Subroutine Calc_BVS(A, Ipr, N_BVSm, BVS_M, Filecod)
    !!----    type (Atoms_Conf_List_type),              intent(in)   :: A            !  In -> Object of Atoms_Conf_List_type
    !!----    integer,                        optional, intent(in)   :: Ipr
    !!----    integer,                        optional, intent(in)   :: n_bvsm       !  In -> Number of modifications
    !!----    character(len=*), dimension(:), optional, intent(in)   :: BVS_M        ! In -> Text with BVS parameters
    !!----    character(len=*),               optional, intent(in)   :: Filecod
    !!----
    !!----    Subroutine to calculate Bond-Valence sums.
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Calc_Dist_Angle_Sigma" in order
    !!----    to update the internal private variables related to distance/angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    Control for error is present.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Calc_BVS(A, Ipr, N_Bvsm, Bvs_M, Filecod)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),            intent(in)  :: A      !  In -> Object of Atoms_Conf_List_type
       integer,                      optional, intent(in)  :: Ipr
       integer,                      optional, intent(in)  :: N_bvsm
       character(len=*),dimension(:),optional, intent(in)  :: bvs_m
       character(len=*),             optional, intent(in)  :: Filecod

       !---- Local variables ----!
       integer           :: i,j,k,n1,n2,icm,l,icn,isoc,kk,is0, &
                            isstr,isdav,isigt,ibvs,sig1,sig2
       real(kind=cp)     :: tol,fact,del2,s0,q2,  &
                            dd,sigtot,efcn,sums,dav,sdav,q1,d2,  &
                            str2,sstr,ric,r_2,del,perc,spred,disp,  &
                            str,rg1,dist,gii_a,gii_b,gii_c
       character(len=4)  :: rnatom

       call init_err_conf()

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20


       if (present(ipr)) then
          write(unit=ipr,fmt="(/,a)")   "  ------------------------------------------------"
          write(unit=ipr,fmt="(a)")     "  {--- BOND-VALENCE AND POLYHEDRA DISTORTIONS ---}"
          write(unit=ipr,fmt="(a,/)")   "  ------------------------------------------------"
       end if


       if (present(n_bvsm).and. present(ipr)) then
          write(unit=ipr,fmt="(a,/,a,/)")  &
          " Bond-Valence parameters (d0,B0) for Zachariasen formula:  s= exp{(d0-d)/B0}", &
                "   (read from external conditions)"
       else
          write(unit=ipr,fmt="(a,/,a,/)")  &
          " Bond-Valence parameters (d0,B0) for Zachariasen formula:  s= exp{(d0-d)/B0}", &
                "   (data read from internal table)"
       end if

       !---- Each line: Cation Anion d0 b
       if (present(n_bvsm) .and. present(bvs_m)) then
          call Complete_Table(A,N_bvsm,bvs_m)
       end if

       do n1=1,A%N_Cations
          do j=1,A%N_Anions
             n2=A%N_Cations+j
             if (present(ipr)) then
                write(unit=ipr,fmt="(2(a,i3,a,a4),2(a,f6.3),a)") &
                      "   Type",n1,": ",A%Species(n1)," with type",n2,": ",A%Species(n2),&
                      " d0=",Table_d0(n1,n2),"    B0=",Table_b(n1,n2), "   => Reference: "//trim(references(Table_ref(n1,n2)))
                write(unit=ipr,fmt="(2(a,a,a,f6.3,a),/)") &
                      "   Cation (Eff. radius): ",A%Species(n1),"(",A%Radius(n1),")   ", &
                      "   Anion  (Eff. radius): ",A%Species(n2),"(",A%Radius(n2),")"
             end if
          end do
       end do

       del2=0.0
       call Get_LogUnit(ibvs)

       if (present(filecod)) then
          open(unit=ibvs,file=trim(filecod)//"_sum.bvs",status="replace",action="write", position="append")
          write(unit=ibvs,fmt="(a)") "  Subroutine Calc_BVS (JRC-LLB, version: March-2005)"
          write(unit=ibvs,fmt="(a)") "  Title: Summary of Bond-Valence calculations for file: "//trim(filecod)//".cfl"
          write(unit=ibvs,fmt="(a)") "   Atom      Coord  D_aver Sigm   Distort(x10-4)    Valence    BVSum(Sigma)"
       else
          open(unit=ibvs,file="summary.bvs",status="replace",action="write", position="append")
          write(unit=ibvs,fmt="(a)") "  Subroutine Calc_BVS (JRC-LLB, version: March-2005)"
          write(unit=ibvs,fmt="(a)") "  Title: Summary of Bond-Valence calculations "
          write(unit=ibvs,fmt="(a)") "   Atom      Coord  D_aver Sigm   Distort(x10-4)    Valence    BVSum(Sigma)"
       end if

       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii_a=0.0
       gii_b=0.0
       gii_c=0.0
       do i=1,A%natoms
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0_cp,q1)
          icn=0
          efcn=0.0
          sums=0.0
          sigtot=0.0
          dav=0.0
          sdav=0.0
          str2=0.0
          d2=0.0
          isoc=INT(A%atom(i)%VarF(2)*1000.0+0.5)
          if (present(ipr)) then
             write(unit=ipr,fmt="(/,/,a,/,a,a4,a,f6.3,a,i3,a,/,a,/)") &
                  "    -------------------------------------------------------------------",  &
                  " => Bond-valence and coordination of atom: ",A%atom(i)%lab ," occupancy: ",A%atom(i)%VarF(1),"(",isoc,")",  &
                  "    -------------------------------------------------------------------"
          end if
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(j,i))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(j,i))%charge
             sig2= SIGN(1.0_cp,q2)
             if (sig1 == sig2) cycle
             dd=coord_info%dist(j,i)
             if (dd > (A%Radius(l)+A%Radius(k))*(1.0+tol)) cycle
             icn=icn+1
             efcn=efcn+A%Atom(coord_info%n_cooatm(j,i))%VarF(1) !/A%Atom(i)%VarF(1)
             kk=k
             d2=d2+dd*dd
             s0=coord_info%s_dist(j,i)
             dav=dav+dd
             sdav=sdav+s0*s0
             rnatom=A%atom(coord_info%n_cooatm(j,i))%lab
             call Bond_valence(Table_d0(l,k),Table_b(l,k),dd,s0,str,sstr)
             str=str*A%Atom(coord_info%n_cooatm(j,i))%VarF(1)             !/A%Atom(i)%VarF(1)
             sstr=sstr*A%Atom(coord_info%n_cooatm(j,i))%VarF(1)           !/A%Atom(i)%VarF(1)
             sums=sums+str
             str2=str2+str*str
             sigtot=sigtot+sstr*sstr
             is0=nint(s0*10000)
             isstr=nint(sstr*1000)
             if (present(ipr)) then
                write(unit=ipr,fmt="(a,a4,a,a4,a,f8.4,a,i4,a,f6.3,a,i3,a)")  &
                     " (",A%atom(i)%lab ,")-(",rnatom,") :",dd,"(",is0,")  ",str,"(",isstr,")"
             end if
          end do

          ric=real(icn)
          if (icn == 0) then
            if (present(ipr) ) then
                write(unit=ipr,fmt=*) " => Warning!! atom: ",A%atom(i)%lab ," is non-coordinated"
                write(unit=ipr,fmt=*) "    Increase the tolerance factor for ionic radii"
             end if
             cycle
          end if
          d2=d2/ric
          sigtot=SQRT(sigtot)
          dav=dav/ric
          sdav=SQRT(sdav)/ric

          isdav=INT(sdav*10000+0.5)
          isigt=INT(sigtot*1000+0.5)
          dist=10000.0*(d2/(dav*dav)-1.0)
          r_2=sums/ric
          r_2=SQRT(ABS(str2/ric-r_2*r_2))
          del=sums-ABS(q1)
          del2=del2+del*del
          perc=100.0*ABS(del/q1)

          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii_a=gii_a+abs(del)*fact
          gii_b=gii_b+perc*fact
          gii_c=gii_c+del*del*fact

          !efcn=efcn/A%Atom(i)%VarF(1)                    !Division by the site occupancy
          spred=ABS(q1)/efcn                              !Predicted valence of a single bond
          disp=table_d0(l,kk)-table_b(l,kk)*LOG(spred)    !Pred. distance
          if (present(ipr)) then
             write(unit=ipr,fmt="(/,a,i5,a,f5.2,a,a4)") " Coordination number: ",icn, &
                 "      Eff.Coor. number: ",efcn,"  for atom: ",A%atom(i)%lab
             write(unit=ipr,fmt="(a,f8.4,a,i4,a,f8.3,a)")  &
                  " Average distance  :",dav,"(",isdav,")  Distortion:  ",dist," xE-04"
             write(unit=ipr,fmt="(a,f8.4,a,f6.3)") " Predicted distance:",disp, &
                  "        Single bond-valence S=",spred
             write(unit=ipr,fmt="(a,f8.3,/,a,f8.3,a,i3,a)") &
                  "                                      Valence: ",q1,  &
                  "                                         Sums: ",sums, "(",isigt,")"
             write(unit=ipr,fmt="(a,2f8.3,/,a)") &
                  " Deviation from the Valence Sum Rule (r1,%dev):",del,perc,  &
                  " {r1=Sumj(sij)-Vi, %dev=100abs(r1)/Vi} "
             write(unit=ipr,fmt="(a,f8.3,/,a)")  &
                  " Deviation from the Equal Valence Rule    (r2):",r_2,  &
                  " {r2=<sij-<sij>>rms}"
             write(unit=ibvs,fmt="(tr4,a4,tr4,f6.2,f8.4,a,i4,a,f14.3,2f12.3,a,i3,a)") &
                  A%atom(i)%lab,efcn,dav,"(",isdav,")",dist,q1,sums, "(",isigt,")"

          end if
       end do

       rg1=SQRT(del2/real(A%natoms))*100.0
       gii_a=gii_a*100.0/A%totatoms
       gii_b=gii_b/A%totatoms  !*100.0 already multiplied
       gii_c=sqrt(gii_c/A%totatoms)*100.0

       if (present(ipr)) then
          write(unit=ipr,fmt="(/,6(a,/))")  &
             " => Lines concerning predicted average distances and single",  &
             "    bond-valence values, as well as the deviations from the",  &
             "    Equal Valence Rule,  apply only to those central atoms",  &
             "    having N coordination-atoms of the same chemical species. ",  &
             "    (The term 'single bond-valence' refers to the valence value",  &
             "    of a single bond for a regular polyhedron, so S=Valence/N)"
           write(unit=ipr,fmt="(/,4(a,/))")  &
             " => The Old Global Instability Index (GII) is calculated with the atoms of the asymetric unit (Num_Atoms).",&
             "    The normalized GII(a,b,c) below are calculated using the sum over asymmetric unit but multiplying ",&
             "    differences by the multiplicity of the site. N_Atoms_UCell is the total number of atoms in the ", &
             "    conventional unit cell. In all cases the result of the different expressions is multiplied by 100.0"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Old Global Instability Index (  GII=SQRT{SUM{|BVS-abs(q)|^2}/Num_Atoms} ) =", &
               rg1," /100"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(a)=       SUM {|BVS-abs(q)|  *mult}       /N_Atoms_UCell =", &
               gii_a," /100"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(b)=       SUM {|BVS-abs(q)|  *mult/abs(q)}/N_Atoms_UCell =", &
               gii_b," %"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(c)= SQRT{ SUM {|BVS-abs(q)|^2*mult}       /N_Atoms_UCell}=", &
               gii_c," /100"
       end if

       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Old Global Instability Index (  GII=SQRT{SUM{|BVS-abs(q)|^2}/Num_Atoms} ) =", &
            rg1," /100"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(a)=       SUM {|BVS-abs(q)|  *mult}       /N_Atoms_UCell =", &
            gii_a," /100"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(b)=       SUM {|BVS-abs(q)|  *mult/abs(q)}/N_Atoms_UCell =", &
            gii_b," %"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(c)= SQRT{ SUM {|BVS-abs(q)|^2*mult}       /N_Atoms_UCell}=", &
            gii_c," /100"

       call flush(ibvs)
       close (unit=ibvs)

       return
    End Subroutine Calc_BVS

    !!----
    !!---- Subroutine Calc_Map_BVS(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta)
    !!----    type (Atoms_Conf_List_type), intent(in) :: A
    !!----    type (Space_Group_Type),     intent(in) :: SpG
    !!----    Type (Crystal_Cell_Type),    intent(in) :: Cell
    !!----    character(len=*),            intent(in) :: Filecod
    !!----    integer,                     intent(in) :: ndimx
    !!----    integer,                     intent(in) :: ndimy
    !!----    integer,                     intent(in) :: ndimz
    !!----    character(len=*),            intent(in) :: atname
    !!----    real(kind=cp),               intent(in) :: drmax
    !!----    real(kind=cp), optional,     intent(in) :: delta !Tolerance in v.u. for output the map
    !!----
    !!----    Calculate a map of BVS values where each point of the grid is determined
    !!----    by a species representative defined in atname. The BVS value is evaluated
    !!----    for distances below drmax value. If delta is present only the points with valence
    !!----    close ( q-delta < BVS < q+delta) to that of the atname are output in the map.
    !!----
    !!---- Update: December - 2007, December 2014 ( JRC,change the order of indices and VESTA output)
    !!
    Subroutine Calc_Map_BVS(A,Spg,Cell,Filecod,ndimx,ndimy,ndimz,atname,drmax,delta)
       !---- Arguments ----!
       type (Atoms_Conf_List_type), intent(in) :: A
       type (Space_Group_Type),     intent(in) :: SpG
       Type (Crystal_Cell_Type),    intent(in) :: Cell
       character(len=*),            intent(in) :: Filecod
       integer,                     intent(in) :: ndimx
       integer,                     intent(in) :: ndimy
       integer,                     intent(in) :: ndimz
       character(len=*),            intent(in) :: atname  !Species of the probe atom e.g. LI+1
       real(kind=cp),               intent(in) :: drmax
       real(kind=cp), optional,     intent(in) :: delta

       !---- Local variables ----!
       character(len=4)                             :: car,atm,chem
       integer                                      :: i,j,k,n,n1,n2,np,L,nL
       integer                                      :: nx1,nx2,ny1,ny2,nz1,nz2,nn
       integer                                      :: i1,j1,k1
       integer                                      :: jbvs
       integer,      dimension(:),     allocatable  :: ind
       real(kind=cp)                                :: rx1,ry1,rz1,qval,dif,q1,q2,sig1,sig2,rep
       real(kind=cp)                                :: sbvs, dd, occ, radius,sig,dmin
       real(kind=cp), dimension(3)                  :: pto,pta,step,extend
       real(kind=cp), dimension(:,:,:), allocatable :: map_bvs
       real(kind=cp), dimension(:,:),   allocatable :: peaks
       real(kind=cp), dimension(:),     allocatable :: VD_peaks
       type (Atom_list_Type)                        :: At1, At2
       logical                                      :: new_peak,anion,cation

       !---- Initial conditions ----!
       if (A%natoms <= 0) return
       if (ndimx <= 0 .or. ndimy <= 0 .or. ndimz <= 0) return

       !---- Preparing the Variables and calculations ----!
       call Allocate_Atom_List(A%natoms,At1)
       At1%atom=A%atom  !Atom list A%Atom(:)%ind(1) contains the pointer to the species
                        !coming from A (Atoms_Conf_List_type, set in the main program)
       atm=u_case(atname)
       n1=0
       do i=1,At1%natoms  !This exclude all atoms with the same species as Atname
          n2=A%Atom(i)%ind(1)
          if (trim(u_case(A%species(n2))) == trim(atm)) then
             At1%atom(i)%active=.false.
             q1=At1%atom(i)%charge
             sig1=SIGN(1.0_cp,q1)
             n1=n2
             radius=A%radius(n1)
             car=A%Species(n1)
          end if
       end do
       if (n1 ==0) then
          err_conf=.true.
          ERR_Conf_Mess="The point atom "//atname//" is not in the Species Table"
          return
       end if

       call AtList1_ExtenCell_AtList2(Spg,At1,At2,.true.)
       call Deallocate_atom_list(At1)
       !check that all species are well set in the list
       do n=1,At2%natoms
           n2=At2%Atom(n)%ind(1)
           if (n2 ==0) then
              err_conf=.true.
              ERR_Conf_Mess="The atom "//trim(At2%Atom(n)%lab)//" is not in the Species Table"
              return
           end if
       end do

       allocate(map_bvs(ndimx,ndimy,ndimz))
       map_bvs=0.0

       if(.not. present(delta)) then
          allocate(peaks(4,ndimx)) ! A maximum of ndimx peaks will be stored
          peaks=0.0
          allocate(VD_peaks(ndimx)) ! Vector holding the differences of bond-valences
          VD_peaks=0.0
          allocate(ind(ndimx))      ! Vector pointer for ordering the peaks
          ind=0
       end if

       ! Determination of the valence expected for a good position
       i=index(car,"+")
       cation=.false.
       anion=.false.
       if( i /= 0) then
         read(unit=car(i+1:),fmt=*) qval
         cation=.true.
       else
         i=index(car,"-")
         read(unit=car(i+1:),fmt=*) qval
         anion=.true.
       end if
       step=(/ 1.0/real(ndimx),  1.0/real(ndimy), 1.0/real(ndimz) /)
       extend=(/ drmax/cell%cell(1),  drmax/cell%cell(2), drmax/cell%cell(3) /)

       np=0
       !---- Map Calculation ----!
       do k=1,ndimz
          pto(3)=(k-1)*step(3)
          do j=1,ndimy
             pto(2)=(j-1)*step(2)
             do i=1,ndimx
                pto(1)=(i-1)*step(1)

                rx1=pto(1)-extend(1)
                if (rx1 <= 0.0) then
                    nx1=int(rx1)-1
                else
                    nx1=int(rx1)
                end if
                nx2=int(pto(1)+extend(1))

                ry1=pto(2)-extend(2)
                if (ry1 <= 0.0) then
                    ny1=int(ry1)-1
                else
                    ny1=int(ry1)
                end if
                ny2=int(pto(2)+extend(2))

                rz1=pto(3)-extend(3)
                if (rz1 <= 0.0) then
                    nz1=int(rz1)-1
                else
                    nz1=int(rz1)
                end if
                nz2=int(pto(3)+extend(3))

                sbvs=0.0
                rep=0.0
                do n=1,At2%natoms
                   n2=At2%Atom(n)%ind(1)
                   !if((cation .and. n2 <= A%N_cations)  .or. (anion .and. n2 > A%N_cations) ) cycle
                   q2=At2%Atom(n)%charge
                   sig2= SIGN(1.0_cp,q2)
                   sig=(radius+A%radius(n2))*0.99
                   do k1=nz1,nz2
                      do j1=ny1,ny2
                         do i1=nx1,nx2
                            pta=At2%Atom(n)%x+real((/i1,j1,k1/))
                            occ=At2%Atom(n)%VarF(1)
                            dd=Distance(pto,pta,Cell)
                            if (dd > drmax) cycle
                            if (sig1 == sig2) then
                               rep=rep + (sig/dd)**18
                               cycle
                            end if
                            sbvs=sbvs+occ*exp((table_d0(n1,n2)-dd)/table_b(n1,n2))
                         end do
                      end do
                   end do
                end do
                dif=abs(sbvs-qval)
                if(present(delta)) then
                  if(dif > delta .or. rep > 0.01) sbvs=0.0_cp
                else
                  !Algorithm for selecting the optimum positions
                  if(dif < 0.15 .and. np < ndimx) then
                    new_peak=.true.
                    do L=1,np
                       dd=Distance(peaks(1:3,L),pto,Cell)
                       if( dd < 0.8 ) then
                         new_peak=.false.
                         nL=L
                         exit
                       end if
                    end do
                    if(new_peak) then
                      np=np+1
                      peaks(1:3,np)= pto
                      peaks(4,np)  = sbvs
                      VD_peaks(np) = dif
                    else  !now compare with the peak stored at nL and interchange them if sbvs is more favourable
                      if( dif < abs(qval-peaks(4,nL)) ) then
                        peaks(1:3,nL) = pto
                        peaks(  4,nL) = sbvs
                        VD_peaks(nL)  = dif
                      end if
                    end if
                  end if
                end if
                !end of peaks construction
                map_bvs(i,j,k)=sbvs
             end do
          end do
       end do
       !Sorting the favourable peak positions for inserting the species of atom Atname

       if(.not. present(delta)) call sort(VD_peaks,np,ind)
       !---- Export a File ----!
       !call Get_LogUnit(jbvs)
       open(newunit=jbvs,file=trim(filecod)//".map",status="replace",action="write")

       write (unit=jbvs, fmt='(a)') "BVS Map Calculations using Bond_STR Program"
       write (unit=jbvs, fmt='(a)') "BVS Map for species "//trim(car)
       write (unit=jbvs, fmt='(a,3f12.4,3x,3f8.3)') "CELL ",cell%cell,cell%ang
       write (unit=jbvs, fmt='(a)')     "SPGR  "//trim(SpG%spg_symb)
       write (unit=jbvs, fmt='(a)')     "! List ot atoms  "
       do i=1,A%natoms
         write(unit=jbvs, fmt='(a,t20,5f12.5)')"Atom "//trim(A%Atom(i)%lab)//"  "//A%Atom(i)%SfacSymb, &
                                                A%Atom(i)%x, A%Atom(i)%biso, A%Atom(i)%occ
       end do
       if(.not. present(delta)) then
         write (unit=jbvs, fmt='(a)')     "! List ot favourable positions for inserting the species "//trim(car)
         do i=1,np
           j=ind(i)
           write(unit=jbvs, fmt='(a,i4,a,3f10.5,a,f8.3)')"#",i,"  Position: (",peaks(1:3,j),")  Valence: ",peaks(4,j)
         end do
       end if
       write (unit=jbvs, fmt='(a)')     "! Grid values: ndimx,ndimy,ndimz [0,ndimx-1] [0,ndimy-1] [0,ndimz-1]  "
       write (unit=jbvs, fmt='(a,9i6)') "GRID ",ndimx,ndimy,ndimz,0,ndimx-1,0,ndimy-1,0,ndimz-1
       write (unit=jbvs, fmt='(a)')     "DENSITY_MAP"
       write (unit=jbvs,fmt='(8g12.5)') map_bvs
       close(unit=jbvs)
       call write_grid_VESTA(map_bvs,cell,"Bond Valence Map",trim(filecod)//"_bvs","P")

       !---- End Procedure ----!
       if (allocated(map_bvs)) deallocate(map_bvs)
       call deallocate_atom_list(At2)

       return
    End Subroutine Calc_Map_BVS

    !!----
    !!---- Subroutine Complete_Table(A,N_bvsm,bvs_m)
    !!----    type (Atoms_Conf_List_type),   intent(in) :: A      !  In -> Object of Atoms_Conf_List_type
    !!----    Integer,                       intent(in) :: N_bvsm
    !!----    Character(len=*),dimension(:), intent(in) :: bvs_m
    !!----
    !!----    Sets up the table of BV parameters (add provided external values)
    !!----    completing the table when the user gives his/her own BV parameters
    !!----
    !!---- Update: January - 2008
    !!
    Subroutine Complete_Table(A,N_bvsm,bvs_m)
       type (Atoms_Conf_List_type),   intent(in) :: A      !  In -> Object of Atoms_Conf_List_type
       Integer,                       intent(in) :: N_bvsm
       Character(len=*),dimension(:), intent(in) :: bvs_m
       !---- Local variables ----!
       Integer       :: i,j,k,ic,icat,ian
       real(kind=cp) :: d0_n, b_n
       character(len=10), dimension(5)        :: dire

          do k=1,N_bvsm
                 dire=" "
                 call getword(bvs_m(k),dire,ic)
                 if (ic <= 0) cycle
                 icat=0
                 ian=0
                 if (ic < 3 ) then
                          err_conf=.true.
                          ERR_Conf_Mess="Cation-Anion D0,b parameters must be provided"
                          return
                 end if

                 do i=1,A%N_Spec
                          if (u_case(dire(1)(1:4)) /= A%Species(i)) cycle
                          icat=i
                 end do
                 do i=1,A%N_Spec
                          if (u_case(dire(2)(1:4)) /= A%Species(i)) cycle
                          ian=i
                 end do
                 if (icat == 0 .or. ian == 0) then
                    err_conf=.true.
                    ERR_Conf_Mess="The given Cation or Anion cannot be found in the atom list"
                    return
                 end if
                 if (icat > ian) then
                          j=icat
                          icat=ian
                          ian=j
                 end if
                 if (icat > A%N_Cations) then
                    err_conf=.true.
                    ERR_Conf_Mess="One of the given cations is not found in the atom list"
                    return
                 end if

                 read(unit=dire(3),fmt=*) d0_n
                 if (ic > 3) read(unit=dire(4),fmt=*) b_n

                 Table_d0(icat,ian)=d0_n
                 Table_b(icat,ian)=b_n
                 Table_ref(icat,ian)=0

                 Table_d0(ian,icat)=d0_n
                 Table_b(ian,icat)=b_n
                 Table_ref(ian,icat)=0

          end do

    End Subroutine Complete_Table

    !!----
    !!---- Subroutine Complete_Table_BVEL(A,N_bvel,bvel)
    !!----    type (Atoms_Conf_List_type),   intent(in) :: A      !  In -> Object of Atoms_Conf_List_type
    !!----    Integer,                       intent(in) :: N_bvel
    !!----    Character(len=*),dimension(:), intent(in) :: bvel
    !!----
    !!----    Sets up the table of BVEL parameters (add provided external values)
    !!----    Completing the table when the user gives his/her own BVEL parameters
    !!----
    !!---- Created: December - 2014
    !!
    Subroutine Complete_Table_BVEL(A,N_bvel,bvel)
       type (Atoms_Conf_List_type),   intent(in) :: A      !  In -> Object of Atoms_Conf_List_type
       Integer,                       intent(in) :: N_bvel
       Character(len=*),dimension(:), intent(in) :: bvel
       !---- Local variables ----!
       Integer       :: i,j,k,ic,icat,ian
       real(kind=cp) :: Avcoor,Rzero,Rcutoff,Dzero,Rmin,alpha
       character(len=12), dimension(9)  :: dire

          do k=1,N_bvel
                 dire=" "
                 call getword(bvel(k),dire,ic)
                 if (ic <= 0) cycle
                 icat=0
                 ian=0
                 if (ic < 8 ) then
                          err_conf=.true.
                          ERR_Conf_Mess="Cation-Anion {Nc,R0,D0,Rmin,alpha} parameters must be provided"
                          return
                 end if

                 do i=1,A%N_Spec
                          if (u_case(dire(1)(1:4)) /= A%Species(i)) cycle
                          icat=i
                 end do
                 do i=1,A%N_Spec
                          if (u_case(dire(2)(1:4)) /= A%Species(i)) cycle
                          ian=i
                 end do
                 if (icat == 0 .or. ian == 0) then
                    err_conf=.true.
                    ERR_Conf_Mess="The given Cation or Anion cannot be found in the atom list"
                    return
                 end if
                 if (icat > ian) then
                          j=icat
                          icat=ian
                          ian=j
                 end if
                 if (icat > A%N_Cations) then
                    err_conf=.true.
                    ERR_Conf_Mess="One of the given cations is not found in the atom list"
                    return
                 end if

                 read(unit=dire(3),fmt=*) Avcoor
                 read(unit=dire(4),fmt=*) Rzero
                 read(unit=dire(5),fmt=*) Rcutoff
                 read(unit=dire(6),fmt=*) Dzero
                 read(unit=dire(7),fmt=*) Rmin
                 read(unit=dire(8),fmt=*) alpha

                 Table_Avcoor  (icat,ian)=Avcoor
                 Table_Rzero   (icat,ian)=Rzero
                 Table_Rcutoff (icat,ian)=Rcutoff
                 Table_Dzero   (icat,ian)=Dzero
                 Table_Rmin    (icat,ian)=Rmin
                 Table_alpha   (icat,ian)=alpha
                 Table_ref     (icat,ian)=0

                 Table_Avcoor  (ian,icat)=Avcoor
                 Table_Rzero   (ian,icat)=Rzero
                 Table_Rcutoff (ian,icat)=Rcutoff
                 Table_Dzero   (ian,icat)=Dzero
                 Table_Rmin    (ian,icat)=Rmin
                 Table_alpha   (ian,icat)=alpha
                 Table_ref     (ian,icat)=0

          end do

    End Subroutine Complete_Table_BVEL
    !!----
    !!---- Subroutine Cost_BVS(A, GII, ERep, gic)
    !!----    type (Atoms_Conf_List_type),  intent(in out) :: A    !  In  -> Object of Atoms_Conf_List_type
    !!----    real(kind=cp),                intent(out)    :: GII  !  OUT -> Global instability index
    !!----    real(kind=cp),      optional, intent(out)    :: ERep !  OUT -> Repulstion term from soft spheres
    !!----    character(len=*),   optional, intent(in)     :: gic  ! If present GII_c is put in GII
    !!----
    !!----    Subroutine to calculate the Global Instability Index.
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Set_TDist_Coordination" in order
    !!----    to update the internal private variables related to distance/angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    All items corresponding to the bond-valence parameters contained in A have to
    !!----    be properly set before calling this procedure.
    !!----
    !!---- Update: April - 2005, July -2104 (JRC, calculation only for the subset having VarF(5) = 0)
    !!
    Subroutine Cost_BVS(A, GII, ERep,gic)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),  intent(in out)  :: A    !  In -> Object of Atoms_Conf_List_type
       real(kind=cp),                intent(out)     :: GII  !  GII_a
       real(kind=cp),      optional, intent(out)     :: ERep !  Repulsion term due to shoft spheres
       character(len=*),   optional, intent(in)      :: gic  !  If present GII_c is put in GII

       !---- Local variables ----!
       integer       :: i,j,k,icm,l,sig1,sig2
       real(kind=cp) :: tol,fact,q2,dd,sums,q1, del, bv,gii_a,gii_c,efcn,sig,rep

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20
       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii_a=0.0
       gii_c=0.0
       rep=0.0
       do i=1,A%natoms
          if(A%Atom(i)%VarF(5) > 0.01) cycle
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0_cp,q1)
          sums=0.0
          efcn=0.0
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(j,i))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(j,i))%charge
             sig2= SIGN(1.0_cp,q2)
             sig=A%radius(l)+A%radius(k)
             dd=coord_info%dist(j,i)
             if (sig1 == sig2) then
                rep=rep + (0.8*sig/dd)**18
                cycle
             end if
             if (dd > sig*(1.0+tol)) cycle
             efcn=efcn+A%Atom(coord_info%n_cooatm(j,i))%VarF(1)
             bv=EXP((Table_d0(l,k)-dd)/Table_b(l,k))
             bv=bv*A%Atom(coord_info%n_cooatm(j,i))%VarF(1) !Occupacy
             sums=sums+bv
          end do
          A%Atom(i)%varF(3)=efcn
          del=sums-ABS(q1)

          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii_a=gii_a+abs(del)*fact
          gii_c=gii_c+del*del*fact

       end do
       gii_a=gii_a*100.0/A%totatoms
       gii_c=sqrt(gii_c/A%totatoms)*100.0
       GII=gii_a
       if(present(ERep)) ERep=rep
       if(present(gic)) GII=gii_c
       return
    End Subroutine Cost_BVS

    !!----
    !!---- Subroutine Cost_BVS_CoulombRep(A, GII, ERep)
    !!----    type (Atoms_Conf_List_type),  intent(in out):: A     !  In  -> Object of Atoms_Conf_List_type
    !!----    real(kind=cp),                intent(out)   :: GII   !  OUT -> Global instability index
    !!----    real(kind=cp),                intent(out)   :: ERep   !  Pseudo Repulsion Coulomb "energy"
    !!----
    !!----
    !!----    Subroutine to calculate the Global Instability Index Gii_a and
    !!----    a pseudo Coulomb repulsion energy useful to avoid cation-cation and
    !!----    anion-anion overlap when using this cost function for predicting or
    !!----    solving a ionic crystal structure. It was used in the old program PiXSA,
    !!----    by J. Pannetier, J. Bassas-Alsina, J.Rodriguez-Carvajal and V. Caignaert,
    !!----    in "Prediction of Crystal Structures from Crystal Chemistry Rules by
    !!----    Simulated Annealing", Nature 346, 343-345 (1990).
    !!----    An additional repulsive term (soft sphere) of the form Esph=(sig/d)**18,
    !!----    with sig equal to the sum of ionic radii of the pair, has been introduced.
    !!----
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Set_TDist_Coordination" in order
    !!----    to update the internal Coord_Info variables related to distance and
    !!----    angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    All items corresponding to the bond-valence parameters contained in
    !!----    "A" have to be properly set before calling this procedure.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Cost_BVS_CoulombRep(A, GII, ERep)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),  intent(in out):: A      !  In -> Object of Atoms_Conf_List_type
       real(kind=cp),                intent(out)   :: GII    !  GII_a
       real(kind=cp),                intent(out)   :: ERep   !  Pseudo Repulsion Coulomb "energy"

       !---- Local variables ----!
       integer        :: i,j,k,icm,l,sig1,sig2
       real(kind=cp)  :: tol,fact,q2,dd,sums,q1, del, bv,efcn,sig

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20

       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii=0.0
       Erep=0.0
       do i=1,A%natoms
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0_cp,q1)
          sums=0.0
          efcn=0.0
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(j,i))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(j,i))%charge
             sig2= SIGN(1.0_cp,q2)
             dd=coord_info%dist(j,i)
             sig=A%radius(l)+A%radius(k) !sum of ionic radii of the current pair
             if (sig1 == sig2) then
                Erep=Erep+ q1*q2/dd + (sig/dd)**18  !Coulomb potential + soft sphere repulsion (avoid short distances)
                cycle
             end if
             if (dd > sig*(1.0+tol)) cycle
             efcn=efcn+A%Atom(coord_info%n_cooatm(j,i))%VarF(1)
             bv=EXP((Table_d0(l,k)-dd)/Table_b(l,k))
             bv=bv*A%Atom(coord_info%n_cooatm(j,i))%VarF(1)
             sums=sums+bv
          end do

          del=sums-ABS(q1)
          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii=gii+abs(del)*fact
          A%Atom(i)%varF(3)=efcn
       end do
       gii=gii*100.0/A%totatoms

       return
    End Subroutine Cost_BVS_CoulombRep

    !!----
    !!---- Subroutine Deallocate_Atoms_Conf_List(A)
    !!----    type (Atoms_Conf_List_Type), intent(in out)   :: A  ! In/ Out -> Objet to be deallocated
    !!----
    !!----    De-allocation of objet A of type Atoms_Conf_List. This subroutine
    !!----    should be after using an object of type Atoms_Conf_List that is no
    !!----    more needed.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Deallocate_Atoms_Conf_List(A)
       !---- Arguments ----!
       type (Atoms_Conf_List_Type), intent(in out)   :: A  !Objet to be deallocated

       if (allocated(A%atom)) deallocate (A%atom)
       A%natoms   = 0
       A%n_spec   = 0
       A%n_anions = 0
       A%n_cations= 0

       return
    End Subroutine Deallocate_Atoms_Conf_List

    !!----
    !!---- Subroutine Deallocate_BVEL_Table()
    !!----
    !!----    Deallocating BVEL_Table
    !!----
    !!---- Created: December - 2014
    !!
    Subroutine Deallocate_BVEL_Table()

       if (allocated(BVEL_Table)) deallocate(BVEL_Table)

       return
    End Subroutine Deallocate_BVEL_Table

    !!----
    !!---- Subroutine Deallocate_BVS_Table()
    !!----
    !!----    Deallocating BVS_Table
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Deallocate_BVS_Table()

       if (allocated(BVS_Table)) deallocate(BVS_Table)

       return
    End Subroutine Deallocate_BVS_Table

    !!----
    !!---- Subroutine Init_Err_Conf()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Init_Err_Conf()

       err_conf=.false.
       ERR_Conf_Mess=" "

       return
    End Subroutine Init_Err_Conf

    !!----
    !!----
    !!---- Subroutine Set_BVEL_Table()
    !!----
    !!----  Fills the parameters for BVEL from
    !!----  Stefan Adams and R. Prasada Rao
    !!---- "High power lithium ion battery materials by computational design"
    !!----  Phys. Status Solidi A 208, No. 8, 1746–1753 (2011) / DOI 10.1002/pssa.201001116
    !!----
    !!---- Only on anion is available (O-2) for the moment
    !!----
    !!---- Update: December - 2014
    !!
    Subroutine Set_BVEL_Table()

       if (.not. allocated(BVEL_Table)) allocate(BVEL_Table(bvel_species_n))

        BVEL_Table(  1)=BVEL_Par_Type("H+1",  (/  1.923000/),(/  0.870450/),(/  4.000000/),(/  1.885800/),(/  1.127680/),(/  2.188184/),(/33/))
        BVEL_Table(  2)=BVEL_Par_Type("Li+1", (/  5.021000/),(/  1.170960/),(/  5.500000/),(/  0.988160/),(/  1.940010/),(/  1.937984/),(/33/))
        BVEL_Table(  3)=BVEL_Par_Type("Be+2", (/  4.000000/),(/  1.209030/),(/  5.500000/),(/  2.768820/),(/  1.522170/),(/  1.848429/),(/33/))
        BVEL_Table(  4)=BVEL_Par_Type("B+3",  (/  3.417000/),(/  1.357610/),(/  4.500000/),(/  2.389240/),(/  1.340030/),(/  2.597403/),(/33/))
        BVEL_Table(  5)=BVEL_Par_Type("C+4",  (/  3.000000/),(/  1.398260/),(/  5.000000/),(/  4.791870/),(/  1.200890/),(/  2.237136/),(/33/))
        BVEL_Table(  6)=BVEL_Par_Type("C+2",  (/  1.000000/),(/  1.413680/),(/  5.000000/),(/  2.405530/),(/  1.030980/),(/  2.409639/),(/33/))
        BVEL_Table(  7)=BVEL_Par_Type("N+5",  (/  3.000000/),(/  1.462670/),(/  5.000000/),(/  6.276770/),(/  1.161420/),(/  2.222222/),(/33/))
        BVEL_Table(  8)=BVEL_Par_Type("N+3",  (/  2.000000/),(/  1.407950/),(/  5.000000/),(/  3.810890/),(/  1.137580/),(/  2.232143/),(/33/))
        BVEL_Table(  9)=BVEL_Par_Type("NH4+1",(/  3.467000/),(/  2.033800/),(/  6.000000/),(/  0.405370/),(/  2.453640/),(/  2.262443/),(/33/))
        BVEL_Table( 10)=BVEL_Par_Type("Na+1", (/  6.520000/),(/  1.562250/),(/  6.000000/),(/  0.575230/),(/  2.374330/),(/  2.074689/),(/33/))
        BVEL_Table( 11)=BVEL_Par_Type("Mg+2", (/  5.897000/),(/  1.483980/),(/  5.500000/),(/  1.575540/),(/  1.956270/),(/  1.953125/),(/33/))
        BVEL_Table( 12)=BVEL_Par_Type("Al+3", (/  5.327000/),(/  1.599010/),(/  5.000000/),(/  1.803460/),(/  1.758060/),(/  2.358491/),(/33/))
        BVEL_Table( 13)=BVEL_Par_Type("Si+4", (/  4.100000/),(/  1.608170/),(/  5.000000/),(/  2.857200/),(/  1.535940/),(/  2.314815/),(/33/))
        BVEL_Table( 14)=BVEL_Par_Type("P+5",  (/  4.000000/),(/  1.620380/),(/  5.000000/),(/  3.896350/),(/  1.440660/),(/  2.288330/),(/33/))
        BVEL_Table( 15)=BVEL_Par_Type("P+3",  (/  3.000000/),(/  1.515550/),(/  4.500000/),(/  2.020620/),(/  1.410510/),(/  2.487562/),(/33/))
        BVEL_Table( 16)=BVEL_Par_Type("S+6",  (/  4.000000/),(/  1.642200/),(/  5.000000/),(/  4.967260/),(/  1.381020/),(/  2.267574/),(/33/))
        BVEL_Table( 17)=BVEL_Par_Type("S+4",  (/  3.000000/),(/  1.642820/),(/  5.500000/),(/  3.036720/),(/  1.411880/),(/  2.341920/),(/33/))
        BVEL_Table( 18)=BVEL_Par_Type("Cl+7", (/  4.000000/),(/  1.679460/),(/  5.000000/),(/  5.991000/),(/  1.348010/),(/  2.257336/),(/33/))
        BVEL_Table( 19)=BVEL_Par_Type("Cl+5", (/  3.000000/),(/  1.695520/),(/  5.500000/),(/  4.290890/),(/  1.356530/),(/  2.247191/),(/33/))
        BVEL_Table( 20)=BVEL_Par_Type("Cl+3", (/  2.000000/),(/  1.722650/),(/  5.500000/),(/  3.071190/),(/  1.384410/),(/  2.036660/),(/33/))
        BVEL_Table( 21)=BVEL_Par_Type("K+1",  (/  8.846000/),(/  1.941170/),(/  6.000000/),(/  0.349850/),(/  2.766360/),(/  2.293578/),(/33/))
        BVEL_Table( 22)=BVEL_Par_Type("Ca+2", (/  7.544000/),(/  1.795190/),(/  5.500000/),(/  0.994290/),(/  2.320320/),(/  2.100840/),(/33/))
        BVEL_Table( 23)=BVEL_Par_Type("Sc+3", (/  6.255000/),(/  1.732200/),(/  5.500000/),(/  2.156100/),(/  1.996150/),(/  2.024291/),(/33/))
        BVEL_Table( 24)=BVEL_Par_Type("Ti+4", (/  6.000000/),(/  1.723940/),(/  5.500000/),(/  2.813330/),(/  1.831440/),(/  1.988072/),(/33/))
        BVEL_Table( 25)=BVEL_Par_Type("Ti+3", (/  6.000000/),(/  1.697660/),(/  5.500000/),(/  1.978510/),(/  1.886190/),(/  2.173913/),(/33/))
        BVEL_Table( 26)=BVEL_Par_Type("V+5",  (/  4.166000/),(/  1.794450/),(/  5.500000/),(/  3.695330/),(/  1.602580/),(/  1.960784/),(/33/))
        BVEL_Table( 27)=BVEL_Par_Type("V+4",  (/  5.738000/),(/  1.749320/),(/  5.000000/),(/  2.080470/),(/  1.776380/),(/  2.347418/),(/33/))
        BVEL_Table( 28)=BVEL_Par_Type("V+3",  (/  6.000000/),(/  1.677990/),(/  5.500000/),(/  1.829360/),(/  1.857970/),(/  2.277904/),(/33/))
        BVEL_Table( 29)=BVEL_Par_Type("Cr+6", (/  4.000000/),(/  1.824710/),(/  5.500000/),(/  3.687510/),(/  1.532510/),(/  2.100840/),(/33/))
        BVEL_Table( 30)=BVEL_Par_Type("Cr+5", (/  4.000000/),(/  1.767810/),(/  5.500000/),(/  2.365510/),(/  1.555460/),(/  2.487562/),(/33/))
        BVEL_Table( 31)=BVEL_Par_Type("Cr+4", (/  5.429000/),(/  1.760950/),(/  5.500000/),(/  1.933290/),(/  1.762090/),(/  2.444988/),(/33/))
        BVEL_Table( 32)=BVEL_Par_Type("Cr+3", (/  6.000000/),(/  1.661980/),(/  5.500000/),(/  1.773350/),(/  1.838870/),(/  2.325581/),(/33/))
        BVEL_Table( 33)=BVEL_Par_Type("Mn+7", (/  4.000000/),(/  1.873620/),(/  6.500000/),(/  4.916300/),(/  1.481710/),(/  1.923077/),(/33/))
        BVEL_Table( 34)=BVEL_Par_Type("Mn+6", (/  4.000000/),(/  1.820180/),(/  5.500000/),(/  2.822360/),(/  1.529310/),(/  2.403846/),(/33/))
        BVEL_Table( 35)=BVEL_Par_Type("Mn+5", (/  4.000000/),(/  1.788790/),(/  5.500000/),(/  2.464560/),(/  1.575770/),(/  2.421308/),(/33/))
        BVEL_Table( 36)=BVEL_Par_Type("Mn+4", (/  5.923000/),(/  1.732720/),(/  5.000000/),(/  1.858860/),(/  1.770450/),(/  2.487562/),(/33/))
        BVEL_Table( 37)=BVEL_Par_Type("Mn+3", (/  5.862000/),(/  1.689930/),(/  5.500000/),(/  1.812830/),(/  1.857860/),(/  2.288330/),(/33/))
        BVEL_Table( 38)=BVEL_Par_Type("Mn+2", (/  5.910000/),(/  1.627580/),(/  5.500000/),(/  1.641430/),(/  2.029690/),(/  2.079002/),(/33/))
        BVEL_Table( 39)=BVEL_Par_Type("Fe+4", (/  6.000000/),(/  1.765590/),(/  5.500000/),(/  1.872850/),(/  1.827860/),(/  2.439024/),(/33/))
        BVEL_Table( 40)=BVEL_Par_Type("Fe+3", (/  5.733000/),(/  1.708400/),(/  5.000000/),(/  1.666810/),(/  1.866470/),(/  2.380952/),(/33/))
        BVEL_Table( 41)=BVEL_Par_Type("Fe+2", (/  5.743000/),(/  1.579110/),(/  5.500000/),(/  1.692690/),(/  1.960050/),(/  2.083333/),(/33/))
        BVEL_Table( 42)=BVEL_Par_Type("Ni+3", (/  6.000000/),(/  1.648880/),(/  5.500000/),(/  1.661910/),(/  1.818870/),(/  2.415459/),(/33/))
        BVEL_Table( 43)=BVEL_Par_Type("Ni+2", (/  5.933000/),(/  1.559200/),(/  5.500000/),(/  1.468410/),(/  1.924520/),(/  2.257336/),(/33/))
        BVEL_Table( 44)=BVEL_Par_Type("Co+3", (/  6.000000/),(/  1.592340/),(/  5.500000/),(/  1.870240/),(/  1.776200/),(/  2.304147/),(/33/))
        BVEL_Table( 45)=BVEL_Par_Type("Co+2", (/  5.506000/),(/  1.597730/),(/  5.500000/),(/  1.514760/),(/  1.933620/),(/  2.217295/),(/33/))
        BVEL_Table( 46)=BVEL_Par_Type("Cu+3", (/  4.000000/),(/  1.709640/),(/  5.000000/),(/  1.882420/),(/  1.708230/),(/  2.341920/),(/33/))
        BVEL_Table( 47)=BVEL_Par_Type("Cu+2", (/  2.560000/),(/  1.574220/),(/  5.000000/),(/  1.853410/),(/  1.566330/),(/  2.227171/),(/33/))
        BVEL_Table( 48)=BVEL_Par_Type("Cu+1", (/  2.560000/),(/  1.587300/),(/  5.000000/),(/  0.664170/),(/  1.782690/),(/  2.932551/),(/33/))
        BVEL_Table( 49)=BVEL_Par_Type("Zn+2", (/  4.718000/),(/  1.653440/),(/  5.000000/),(/  1.240310/),(/  1.885570/),(/  2.481390/),(/33/))
        BVEL_Table( 50)=BVEL_Par_Type("Ga+3", (/  4.905000/),(/  1.716060/),(/  5.000000/),(/  1.184560/),(/  1.793910/),(/  2.680965/),(/33/))
        BVEL_Table( 51)=BVEL_Par_Type("Ge+4", (/  4.305000/),(/  1.739390/),(/  5.000000/),(/  1.913750/),(/  1.668720/),(/  2.525253/),(/33/))
        BVEL_Table( 52)=BVEL_Par_Type("As+5", (/  4.029000/),(/  1.766890/),(/  5.000000/),(/  2.719340/),(/  1.581270/),(/  2.433090/),(/33/))
        BVEL_Table( 53)=BVEL_Par_Type("As+3", (/  3.000000/),(/  1.767060/),(/  5.000000/),(/  1.514930/),(/  1.645540/),(/  2.475248/),(/33/))
        BVEL_Table( 54)=BVEL_Par_Type("Se+6", (/  4.000000/),(/  1.798660/),(/  5.500000/),(/  3.448650/),(/  1.532870/),(/  2.403846/),(/33/))
        BVEL_Table( 55)=BVEL_Par_Type("Se+4", (/  3.000000/),(/  1.800950/),(/  5.500000/),(/  2.380820/),(/  1.559570/),(/  2.341920/),(/33/))
        BVEL_Table( 56)=BVEL_Par_Type("Br+7", (/  4.000000/),(/  1.836580/),(/  5.500000/),(/  4.243390/),(/  1.502740/),(/  2.364066/),(/33/))
        BVEL_Table( 57)=BVEL_Par_Type("Rb+1", (/ 10.020000/),(/  2.085970/),(/  6.500000/),(/  0.268130/),(/  2.896830/),(/  2.421308/),(/33/))
        BVEL_Table( 58)=BVEL_Par_Type("Sr+2", (/  9.400000/),(/  1.953110/),(/  5.500000/),(/  0.743510/),(/  2.535890/),(/  2.197802/),(/33/))
        BVEL_Table( 59)=BVEL_Par_Type("Y+3",  (/  7.285000/),(/  1.903840/),(/  5.500000/),(/  1.627010/),(/  2.215230/),(/  2.092050/),(/33/))
        BVEL_Table( 60)=BVEL_Par_Type("Zr+4", (/  6.765000/),(/  1.845050/),(/  5.500000/),(/  2.191030/),(/  1.996020/),(/  2.040816/),(/33/))
        BVEL_Table( 61)=BVEL_Par_Type("Nb+5", (/  6.044000/),(/  1.865880/),(/  5.500000/),(/  2.723260/),(/  1.854590/),(/  2.008032/),(/33/))
        BVEL_Table( 62)=BVEL_Par_Type("Nb+4", (/  6.000000/),(/  1.785430/),(/  6.000000/),(/  2.709600/),(/  1.859890/),(/  1.901141/),(/33/))
        BVEL_Table( 63)=BVEL_Par_Type("Nb+3", (/  6.000000/),(/  1.745810/),(/  6.000000/),(/  2.028480/),(/  1.951900/),(/  1.996008/),(/33/))
        BVEL_Table( 64)=BVEL_Par_Type("Mo+6", (/  4.764000/),(/  1.909340/),(/  5.000000/),(/  1.991500/),(/  1.712540/),(/  2.557545/),(/33/))
        BVEL_Table( 65)=BVEL_Par_Type("W+5",  (/  6.000000/),(/  1.819750/),(/  6.000000/),(/  2.615700/),(/  1.762610/),(/  2.008032/),(/33/))
        BVEL_Table( 66)=BVEL_Par_Type("W+4",  (/  6.000000/),(/  1.745580/),(/  6.000000/),(/  2.471140/),(/  1.819450/),(/  1.923077/),(/33/))
        BVEL_Table( 67)=BVEL_Par_Type("Re+7", (/  4.098000/),(/  1.977920/),(/  6.000000/),(/  3.555930/),(/  1.596340/),(/  1.968504/),(/33/))
        BVEL_Table( 68)=BVEL_Par_Type("Re+6", (/  5.500000/),(/  1.910070/),(/  6.000000/),(/  2.950990/),(/  1.711470/),(/  2.008032/),(/33/))
        BVEL_Table( 69)=BVEL_Par_Type("Re+5", (/  6.000000/),(/  1.826640/),(/  6.000000/),(/  2.410990/),(/  1.769140/),(/  2.087683/),(/33/))
        BVEL_Table( 70)=BVEL_Par_Type("Re+3", (/  6.000000/),(/  2.207100/),(/  6.000000/),(/  0.810670/),(/  2.332180/),(/  2.493766/),(/33/))
        BVEL_Table( 71)=BVEL_Par_Type("Os+8", (/  5.333000/),(/  1.977280/),(/  6.000000/),(/  3.710190/),(/  1.661460/),(/  1.953125/),(/33/))
        BVEL_Table( 72)=BVEL_Par_Type("Os+7", (/  6.000000/),(/  1.957750/),(/  5.500000/),(/  2.919480/),(/  1.728690/),(/  2.087683/),(/33/))
        BVEL_Table( 73)=BVEL_Par_Type("Os+6", (/  6.000000/),(/  1.931920/),(/  5.500000/),(/  2.448710/),(/  1.782800/),(/  2.159827/),(/33/))
        BVEL_Table( 74)=BVEL_Par_Type("Os+4", (/  6.000000/),(/  1.753020/),(/  6.000000/),(/  2.275240/),(/  1.812440/),(/  2.008032/),(/33/))
        BVEL_Table( 75)=BVEL_Par_Type("Ir+5", (/  6.000000/),(/  1.897910/),(/  6.000000/),(/  2.324760/),(/  1.834760/),(/  2.087683/),(/33/))
        BVEL_Table( 76)=BVEL_Par_Type("Ir+4", (/  6.000000/),(/  1.832330/),(/  5.500000/),(/  1.686670/),(/  1.874020/),(/  2.293578/),(/33/))
        BVEL_Table( 77)=BVEL_Par_Type("Pt+4", (/  6.000000/),(/  1.821980/),(/  5.500000/),(/  2.038250/),(/  1.871740/),(/  2.087683/),(/33/))
        BVEL_Table( 78)=BVEL_Par_Type("Pt+2", (/  4.000000/),(/  1.512050/),(/  5.500000/),(/  2.149990/),(/  1.801790/),(/  1.742160/),(/33/))
        BVEL_Table( 79)=BVEL_Par_Type("Au+3", (/  4.000000/),(/  1.817610/),(/  5.500000/),(/  1.969670/),(/  1.813120/),(/  2.008032/),(/33/))
        BVEL_Table( 80)=BVEL_Par_Type("Au+1", (/  2.000000/),(/  1.718190/),(/  5.500000/),(/  0.853040/),(/  1.895430/),(/  2.267574/),(/33/))
        BVEL_Table( 81)=BVEL_Par_Type("Hg+2", (/  6.966000/),(/  1.812760/),(/  5.500000/),(/  1.128520/),(/  2.252750/),(/  2.150538/),(/33/))
        BVEL_Table( 82)=BVEL_Par_Type("Hg+1", (/  4.786000/),(/  1.812800/),(/  5.500000/),(/  0.739310/),(/  2.431550/),(/  2.150538/),(/33/))
        BVEL_Table( 83)=BVEL_Par_Type("Tl+3", (/  5.220000/),(/  2.062970/),(/  5.000000/),(/  0.676370/),(/  2.106420/),(/  2.958580/),(/33/))
        BVEL_Table( 84)=BVEL_Par_Type("Tl+1", (/  8.030000/),(/  1.917520/),(/  6.000000/),(/  0.349990/),(/  2.770860/),(/  2.070393/),(/33/))
        BVEL_Table( 85)=BVEL_Par_Type("Pb+4", (/  5.740000/),(/  2.032930/),(/  5.000000/),(/  1.027190/),(/  2.028570/),(/  2.824859/),(/33/))
        BVEL_Table( 86)=BVEL_Par_Type("Pb+2", (/  7.541000/),(/  2.018250/),(/  5.500000/),(/  0.638330/),(/  2.441910/),(/  2.309469/),(/33/))
        BVEL_Table( 87)=BVEL_Par_Type("Bi+5", (/  6.000000/),(/  2.044980/),(/  5.000000/),(/  1.440500/),(/  1.985990/),(/  2.695418/),(/33/))
        BVEL_Table( 88)=BVEL_Par_Type("Bi+3", (/  6.058000/),(/  2.036770/),(/  5.500000/),(/  0.979040/),(/  2.183210/),(/  2.415459/),(/33/))
        BVEL_Table( 89)=BVEL_Par_Type("Mo+5", (/  5.980000/),(/  1.847600/),(/  5.500000/),(/  2.648020/),(/  1.786700/),(/  2.074689/),(/33/))
        BVEL_Table( 90)=BVEL_Par_Type("Mo+4", (/  6.000000/),(/  1.723900/),(/  6.500000/),(/  3.108070/),(/  1.850990/),(/  1.779359/),(/33/))
        BVEL_Table( 91)=BVEL_Par_Type("Mo+3", (/  5.700000/),(/  1.789330/),(/  5.500000/),(/  1.428260/),(/  1.929740/),(/  2.392344/),(/33/))
        BVEL_Table( 92)=BVEL_Par_Type("Ru+6", (/  4.500000/),(/  1.925790/),(/  5.500000/),(/  2.421090/),(/  1.664310/),(/  2.352941/),(/33/))
        BVEL_Table( 93)=BVEL_Par_Type("Ru+5", (/  6.000000/),(/  1.874420/),(/  5.500000/),(/  2.132080/),(/  1.815710/),(/  2.293578/),(/33/))
        BVEL_Table( 94)=BVEL_Par_Type("Ru+4", (/  6.000000/),(/  1.793630/),(/  5.500000/),(/  1.995130/),(/  1.840530/),(/  2.227171/),(/33/))
        BVEL_Table( 95)=BVEL_Par_Type("Rh+4", (/  6.000000/),(/  1.776750/),(/  5.500000/),(/  1.627250/),(/  1.817930/),(/  2.481390/),(/33/))
        BVEL_Table( 96)=BVEL_Par_Type("Rh+3", (/  6.000000/),(/  1.670130/),(/  5.500000/),(/  1.928260/),(/  1.869150/),(/  2.092050/),(/33/))
        BVEL_Table( 97)=BVEL_Par_Type("Pd+4", (/  5.333000/),(/  1.805000/),(/  5.500000/),(/  2.042180/),(/  1.798130/),(/  2.227171/),(/33/))
        BVEL_Table( 98)=BVEL_Par_Type("Pd+2", (/  4.000000/),(/  1.623590/),(/  5.500000/),(/  1.739100/),(/  1.836710/),(/  2.008032/),(/33/))
        BVEL_Table( 99)=BVEL_Par_Type("Ag+1", (/  4.438000/),(/  1.782390/),(/  5.000000/),(/  0.635190/),(/  2.225780/),(/  2.538071/),(/33/))
        BVEL_Table(100)=BVEL_Par_Type("Cd+2", (/  6.176000/),(/  1.839260/),(/  5.500000/),(/  0.983460/),(/  2.169400/),(/  2.457002/),(/33/))
        BVEL_Table(101)=BVEL_Par_Type("Ta+4", (/  5.500000/),(/  1.756320/),(/  6.000000/),(/  2.756550/),(/  1.798260/),(/  1.831502/),(/33/))
        BVEL_Table(102)=BVEL_Par_Type("In+3", (/  6.024000/),(/  1.903050/),(/  5.000000/),(/  0.840760/),(/  2.024710/),(/  2.832861/),(/33/))
        BVEL_Table(103)=BVEL_Par_Type("Sn+4", (/  6.069000/),(/  1.890190/),(/  5.000000/),(/  1.352680/),(/  1.934220/),(/  2.638522/),(/33/))
        BVEL_Table(104)=BVEL_Par_Type("Sn+2", (/  3.325000/),(/  1.874990/),(/  5.500000/),(/  0.972610/),(/  1.964200/),(/  2.183406/),(/33/))
        BVEL_Table(105)=BVEL_Par_Type("Sb+5", (/  6.000000/),(/  1.897680/),(/  5.500000/),(/  1.955230/),(/  1.863180/),(/  2.500000/),(/33/))
        BVEL_Table(106)=BVEL_Par_Type("Sb+3", (/  6.000000/),(/  1.920360/),(/  5.000000/),(/  1.177860/),(/  2.075260/),(/  2.364066/),(/33/))
        BVEL_Table(107)=BVEL_Par_Type("I+7",  (/  5.800000/),(/  1.922740/),(/  5.500000/),(/  3.214240/),(/  1.741050/),(/  2.386635/),(/33/))
        BVEL_Table(108)=BVEL_Par_Type("I+5",  (/  3.100000/),(/  1.977750/),(/  6.000000/),(/  2.489470/),(/  1.644210/),(/  2.358491/),(/33/))
        BVEL_Table(109)=BVEL_Par_Type("Te+6", (/  6.000000/),(/  1.913430/),(/  5.500000/),(/  2.564060/),(/  1.808760/),(/  2.427184/),(/33/))
        BVEL_Table(110)=BVEL_Par_Type("Te+4", (/  3.396000/),(/  1.952900/),(/  5.500000/),(/  1.671690/),(/  1.752080/),(/  2.493766/),(/33/))
        BVEL_Table(111)=BVEL_Par_Type("Cs+1", (/ 11.790000/),(/  2.258990/),(/  6.500000/),(/  0.233070/),(/  3.131210/),(/  2.386635/),(/33/))
        BVEL_Table(112)=BVEL_Par_Type("Ba+2", (/ 10.320000/),(/  2.159980/),(/  6.000000/),(/  0.579940/),(/  2.737690/),(/  2.288330/),(/33/))
        BVEL_Table(113)=BVEL_Par_Type("La+3", (/  9.830000/),(/  2.063920/),(/  5.500000/),(/  1.185870/),(/  2.469890/),(/  2.217295/),(/33/))
        BVEL_Table(114)=BVEL_Par_Type("Ce+4", (/  7.867000/),(/  2.028210/),(/  5.500000/),(/  1.484120/),(/  2.198720/),(/  2.257336/),(/33/))
        BVEL_Table(115)=BVEL_Par_Type("Ce+3", (/  9.147000/),(/  2.031180/),(/  5.500000/),(/  1.220480/),(/  2.378610/),(/  2.227171/),(/33/))
        BVEL_Table(116)=BVEL_Par_Type("Pr+3", (/  9.067000/),(/  2.036520/),(/  5.500000/),(/  1.170410/),(/  2.371130/),(/  2.277904/),(/33/))
        BVEL_Table(117)=BVEL_Par_Type("Nd+3", (/  8.647000/),(/  2.024250/),(/  5.500000/),(/  1.132050/),(/  2.330160/),(/  2.336449/),(/33/))
        BVEL_Table(118)=BVEL_Par_Type("Sm+3", (/  8.119000/),(/  2.011680/),(/  5.500000/),(/  1.176220/),(/  2.295360/),(/  2.309469/),(/33/))
        BVEL_Table(119)=BVEL_Par_Type("Eu+3", (/  7.743000/),(/  2.004690/),(/  5.500000/),(/  1.195450/),(/  2.268880/),(/  2.304147/),(/33/))
        BVEL_Table(120)=BVEL_Par_Type("Eu+2", (/ 10.111000/),(/  1.891580/),(/  6.000000/),(/  1.130320/),(/  2.538460/),(/  2.024291/),(/33/))
        BVEL_Table(121)=BVEL_Par_Type("Gd+3", (/  8.052000/),(/  1.996540/),(/  5.500000/),(/  1.091610/),(/  2.271900/),(/  2.409639/),(/33/))
        BVEL_Table(122)=BVEL_Par_Type("Tb+4", (/  6.000000/),(/  1.962440/),(/  6.000000/),(/  1.701320/),(/  2.385060/),(/  2.024291/),(/33/))
        BVEL_Table(123)=BVEL_Par_Type("Tb+3", (/  7.958000/),(/  1.956750/),(/  5.500000/),(/  1.207640/),(/  2.235630/),(/  2.309469/),(/33/))
        BVEL_Table(124)=BVEL_Par_Type("Dy+3", (/  7.828000/),(/  1.960290/),(/  5.500000/),(/  1.173500/),(/  2.226890/),(/  2.347418/),(/33/))
        BVEL_Table(125)=BVEL_Par_Type("Ho+3", (/  7.500000/),(/  1.970990/),(/  5.500000/),(/  1.121570/),(/  2.211220/),(/  2.409639/),(/33/))
        BVEL_Table(126)=BVEL_Par_Type("Er+3", (/  7.135000/),(/  1.956080/),(/  5.500000/),(/  1.123940/),(/  2.174770/),(/  2.427184/),(/33/))
        BVEL_Table(127)=BVEL_Par_Type("Tm+3", (/  6.912000/),(/  1.949010/),(/  5.500000/),(/  1.181380/),(/  2.160420/),(/  2.375297/),(/33/))
        BVEL_Table(128)=BVEL_Par_Type("Yb+3", (/  6.875000/),(/  1.928720/),(/  5.500000/),(/  1.219890/),(/  2.142200/),(/  2.347418/),(/33/))
        BVEL_Table(129)=BVEL_Par_Type("Lu+3", (/  6.830000/),(/  1.917280/),(/  5.500000/),(/  1.194880/),(/  2.136000/),(/  2.375297/),(/33/))
        BVEL_Table(130)=BVEL_Par_Type("Hf+4", (/  7.105000/),(/  1.833610/),(/  6.000000/),(/  1.899920/),(/  1.999640/),(/  2.092050/),(/33/))
        BVEL_Table(131)=BVEL_Par_Type("Ta+5", (/  6.090000/),(/  1.868160/),(/  5.500000/),(/  2.366690/),(/  1.855320/),(/  2.057613/),(/33/))
        BVEL_Table(132)=BVEL_Par_Type("W+6",  (/  5.688000/),(/  1.906410/),(/  5.000000/),(/  1.842670/),(/  1.777130/),(/  2.493766/),(/33/))

    End Subroutine Set_BVEL_Table

    !!----
    !!----
    !!---- Subroutine Set_BVS_Table()
    !!----
    !!----    Fills the parameters for BVS from O'Keefe, Bresse, Brown
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Set_BVS_Table()

       if (.not. allocated(BVS_Table)) allocate(BVS_Table(bvs_species_n))

     BVS_Table(  1)=BVS_Par_Type("AC+3", &
                   (/ 2.240, 2.130, 2.630, 2.750, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  2)=BVS_Par_Type("AG+1", &
                   (/ 1.842, 1.800, 2.090, 0.000, 0.000, 2.119, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  3)=BVS_Par_Type("AG+2", &
                   (/ 0.000, 1.790, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  4)=BVS_Par_Type("AG+3", &
                   (/ 0.000, 1.830, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  5)=BVS_Par_Type("AG+9", &
                   (/ 0.000, 0.000, 0.000, 2.220, 2.380, 0.000, 2.260, 2.510, 1.850, 2.220, 2.300, 1.500, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     0,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(  6)=BVS_Par_Type("AL+3", &
                   (/ 1.620, 1.545, 2.032, 2.200, 2.410, 2.210, 2.270, 2.480, 1.790, 2.240, 2.300, 1.450, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     1,     1,     2,     2,     5,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(  7)=BVS_Par_Type("AM+3", &
                   (/ 2.110, 2.000, 2.480, 2.590, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  8)=BVS_Par_Type("AM+4", &
                   (/ 2.080, 1.960, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  9)=BVS_Par_Type("AM+5", &
                   (/ 2.070, 1.950, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 10)=BVS_Par_Type("AM+6", &
                   (/ 2.050, 1.950, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 11)=BVS_Par_Type("AS+2", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 2.240, 2.380, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 12)=BVS_Par_Type("AS+3", &
                   (/ 1.789, 1.700, 2.160, 2.350, 2.580, 2.272, 2.400, 2.650, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     1,     5,     5,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 13)=BVS_Par_Type("AS+5", &
                   (/ 1.767, 1.620, 0.000, 0.000, 0.000, 2.280, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 14)=BVS_Par_Type("AU+1", &
                   (/ 0.000, 0.000, 2.020, 0.000, 2.350, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 15)=BVS_Par_Type("AU+3", &
                   (/ 1.890, 1.890, 2.170, 2.320, 2.540, 2.390, 0.000, 0.000, 1.940, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     2,     5,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 16)=BVS_Par_Type("AU+5", &
                   (/ 0.000, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 17)=BVS_Par_Type("AU+9", &
                   (/ 0.000, 0.000, 0.000, 2.120, 2.340, 2.030, 2.180, 2.410, 1.720, 2.140, 2.220, 1.370, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 18)=BVS_Par_Type("B+3 ", &
                   (/ 1.371, 1.281, 1.740, 1.880, 2.100, 1.770, 1.950, 2.200, 1.470, 1.880, 1.970, 1.140, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     5,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 19)=BVS_Par_Type("BA+2", &
                   (/ 2.285, 2.188, 2.690, 2.880, 3.130, 2.769, 2.880, 3.080, 2.470, 2.880, 2.960, 2.220, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 20)=BVS_Par_Type("BE+2", &
                   (/ 1.381, 1.281, 1.760, 1.900, 2.100, 1.830, 1.970, 2.210, 1.500, 1.950, 2.000, 1.110, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 21)=BVS_Par_Type("BI+2", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.700, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 22)=BVS_Par_Type("BI+3", &
                   (/ 2.094, 1.990, 2.480, 2.590, 2.820, 2.570, 0.000, 0.000, 2.020, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     1,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 23)=BVS_Par_Type("BI+5", &
                   (/ 2.060, 1.970, 2.440, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 24)=BVS_Par_Type("BI+9", &
                   (/ 0.000, 0.000, 0.000, 2.620, 2.840, 2.550, 2.720, 2.870, 2.240, 2.630, 2.720, 1.970, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 25)=BVS_Par_Type("BK+3", &
                   (/ 2.080, 1.960, 2.350, 2.560, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     5,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 26)=BVS_Par_Type("BK+4", &
                   (/ 2.070, 1.930, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 27)=BVS_Par_Type("BR+3", &
                   (/ 1.900, 1.750, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 28)=BVS_Par_Type("BR+5", &
                   (/ 1.840, 1.760, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 29)=BVS_Par_Type("BR+7", &
                   (/ 1.810, 1.720, 2.190, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 30)=BVS_Par_Type("C+2 ", &
                   (/ 1.366, 0.000, 1.410, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 31)=BVS_Par_Type("C+4 ", &
                   (/ 1.390, 1.320, 1.760, 1.910, 0.000, 1.800, 0.000, 0.000, 1.442, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     5,     0,     0,     1,     0,     0,     0,     0,     0 /) )
     BVS_Table( 32)=BVS_Par_Type("C+9 ", &
                   (/ 0.000, 0.000, 0.000, 1.900, 2.120, 1.820, 1.970, 2.210, 1.470, 1.890, 1.990, 1.100, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 33)=BVS_Par_Type("CA+2", &
                   (/ 1.967, 1.842, 2.370, 2.507, 2.720, 2.450, 2.560, 2.760, 2.140, 2.550, 2.620, 1.830, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 34)=BVS_Par_Type("CD+2", &
                   (/ 1.904, 1.811, 2.212, 2.350, 2.570, 2.304, 2.400, 2.590, 1.960, 2.340, 2.430, 1.660, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     1,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 35)=BVS_Par_Type("CE+3", &
                   (/ 2.151, 2.036, 2.520, 2.650, 2.870, 2.650, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 36)=BVS_Par_Type("CE+4", &
                   (/ 2.028, 1.995, 0.000, 0.000, 0.000, 2.650, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 37)=BVS_Par_Type("CE+9", &
                   (/ 0.000, 0.000, 2.410, 2.690, 2.920, 2.620, 2.740, 2.920, 2.340, 2.700, 2.780, 2.040, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 38)=BVS_Par_Type("CF+3", &
                   (/ 2.070, 1.950, 2.450, 2.550, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 39)=BVS_Par_Type("CF+4", &
                   (/ 2.060, 1.920, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 40)=BVS_Par_Type("CL+3", &
                   (/ 1.710, 1.690, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 41)=BVS_Par_Type("CL+5", &
                   (/ 1.670, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 42)=BVS_Par_Type("CL+7", &
                   (/ 1.632, 1.550, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 43)=BVS_Par_Type("CF+3", &
                   (/ 0.000, 0.000, 2.450, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 44)=BVS_Par_Type("CM+3", &
                   (/ 2.230, 2.120, 2.620, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 45)=BVS_Par_Type("CM+4", &
                   (/ 2.080, 1.940, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 46)=BVS_Par_Type("CO+1", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     5,     0,     0 /) )
     BVS_Table( 47)=BVS_Par_Type("CO+2", &
                   (/ 1.692, 1.640, 2.033, 0.000, 0.000, 1.940, 0.000, 0.000, 1.650, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     1,     0,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 48)=BVS_Par_Type("CO+3", &
                   (/ 1.637, 1.620, 2.050, 0.000, 0.000, 2.020, 0.000, 0.000, 1.750, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     9,     2,     2,     0,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 49)=BVS_Par_Type("CO+4", &
                   (/ 1.720, 1.550, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 50)=BVS_Par_Type("CO+9", &
                   (/ 1.655, 0.000, 0.000, 2.180, 2.370, 2.060, 2.240, 2.460, 1.840, 2.210, 2.280, 1.440, 0.000, 0.000 /), &
                   (/ 0.420, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 51)=BVS_Par_Type("CR+2", &
                   (/ 1.730, 1.670, 2.090, 2.260, 2.480, 0.000, 0.000, 0.000, 1.830, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     5,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 52)=BVS_Par_Type("CR+3", &
                   (/ 1.724, 1.657, 2.080, 2.280, 0.000, 2.162, 0.000, 0.000, 1.810, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 53)=BVS_Par_Type("CR+4", &
                   (/ 1.810, 1.560, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 54)=BVS_Par_Type("CR+5", &
                   (/ 1.760, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    23,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 55)=BVS_Par_Type("CR+6", &
                   (/ 1.794, 1.740, 2.120, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 56)=BVS_Par_Type("CR+9", &
                   (/ 1.790, 0.000, 0.000, 2.260, 2.450, 2.180, 2.290, 2.520, 1.850, 2.270, 2.340, 1.520, 0.000, 0.000 /), &
                   (/ 0.340, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 57)=BVS_Par_Type("CS+1", &
                   (/ 2.417, 2.330, 2.791, 2.950, 3.180, 2.890, 2.980, 3.160, 2.830, 2.930, 3.040, 2.440, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     1,     2,     2,     2,     2,     2,     5,     2,     2,     2,     0,     0 /) )
     BVS_Table( 58)=BVS_Par_Type("CU+1", &
                   (/ 1.610, 1.600, 1.858, 2.030, 2.108, 1.898, 1.900, 0.000, 1.520, 1.774, 1.856, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,    12,     5,     1,     1,    12,     0,    12,    12,    12,     0,     0,     0 /) )
     BVS_Table( 59)=BVS_Par_Type("CU+2", &
                   (/ 1.679, 1.594, 2.000, 1.990, 2.160, 2.054, 2.020, 2.270, 1.751, 1.970, 2.080, 1.210, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     1,     2,     2,    10,     2,     2,     2,     0,     0 /) )
     BVS_Table( 60)=BVS_Par_Type("CU+3", &
                   (/ 1.735, 1.580, 2.078, 0.000, 0.000, 0.000, 0.000, 0.000, 1.768, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    20,     5,    12,     0,     0,     0,     0,     0,    12,     0,     0,     0,     0,     0 /) )
     BVS_Table( 61)=BVS_Par_Type("DY+2", &
                   (/ 1.900, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 62)=BVS_Par_Type("DY+3", &
                   (/ 2.001, 1.922, 2.410, 2.530, 2.760, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 63)=BVS_Par_Type("DY+9", &
                   (/ 0.000, 0.000, 0.000, 2.560, 2.770, 2.470, 2.610, 2.800, 2.180, 2.570, 2.640, 1.890, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 64)=BVS_Par_Type("ER+2", &
                   (/ 1.880, 0.000, 0.000, 0.000, 0.000, 2.520, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 65)=BVS_Par_Type("ER+3", &
                   (/ 1.988, 1.904, 2.390, 2.510, 2.750, 2.520, 2.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,    16,    16,     5,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 66)=BVS_Par_Type("ER+9", &
                   (/ 0.000, 0.000, 0.000, 2.540, 2.750, 2.460, 2.590, 2.780, 2.160, 2.550, 2.630, 1.860, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 67)=BVS_Par_Type("ES+3", &
                   (/ 2.080, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 68)=BVS_Par_Type("EU+2", &
                   (/ 2.130, 2.040, 2.530, 2.670, 2.900, 2.584, 0.000, 0.000, 2.340, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,     2,     5,     5,     1,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 69)=BVS_Par_Type("EU+3", &
                   (/ 2.074, 1.961, 2.480, 2.570, 2.790, 2.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     5,    16,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 70)=BVS_Par_Type("EU+9", &
                   (/ 0.000, 0.000, 0.000, 2.610, 2.830, 2.530, 2.660, 2.850, 2.240, 2.620, 2.700, 1.950, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 71)=BVS_Par_Type("FE+2", &
                   (/ 1.734, 1.650, 2.060, 2.210, 2.470, 2.120, 0.000, 0.000, 1.769, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     5,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table( 72)=BVS_Par_Type("FE+3", &
                   (/ 1.759, 1.679, 2.090, 0.000, 0.000, 2.149, 0.000, 0.000, 1.815, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     0,     0,     1,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table( 73)=BVS_Par_Type("FE+4", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 2.230, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 74)=BVS_Par_Type("FE+6", &
                   (/ 1.760, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 75)=BVS_Par_Type("FE+9", &
                   (/ 1.740, 0.000, 0.000, 2.260, 2.470, 2.160, 2.280, 2.530, 1.860, 2.270, 2.350, 1.530, 0.000, 0.000 /), &
                   (/ 0.380, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 76)=BVS_Par_Type("GA+1", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.550 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     5 /) )
     BVS_Table( 77)=BVS_Par_Type("GA+3", &
                   (/ 1.730, 1.620, 2.070, 2.200, 2.460, 2.163, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 78)=BVS_Par_Type("GA+9", &
                   (/ 0.000, 0.000, 0.000, 2.240, 2.450, 2.170, 2.300, 2.540, 1.840, 2.260, 2.340, 1.510, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 79)=BVS_Par_Type("GD+2", &
                   (/ 2.010, 2.400, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 80)=BVS_Par_Type("GD+3", &
                   (/ 2.065, 1.950, 2.445, 2.560, 2.780, 2.530, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 81)=BVS_Par_Type("GD+9", &
                   (/ 0.000, 0.000, 0.000, 2.600, 2.820, 2.530, 2.650, 2.840, 2.220, 2.610, 2.680, 1.930, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 82)=BVS_Par_Type("GE+4", &
                   (/ 1.748, 1.660, 2.140, 0.000, 0.000, 2.217, 2.350, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     1,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 83)=BVS_Par_Type("GE+9", &
                   (/ 0.000, 0.000, 0.000, 2.300, 2.500, 2.230, 2.350, 2.560, 1.880, 2.320, 2.430, 1.550, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 84)=BVS_Par_Type("H+1 ", &
                   (/ 0.569, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.940, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 85)=BVS_Par_Type("HF+3", &
                   (/ 0.000, 2.620, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 86)=BVS_Par_Type("HF+4", &
                   (/ 1.923, 1.850, 2.240, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 87)=BVS_Par_Type("HF+9", &
                   (/ 0.000, 0.000, 0.000, 2.470, 2.680, 2.390, 2.520, 2.720, 2.090, 2.480, 2.560, 1.780, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 88)=BVS_Par_Type("HG+1", &
                   (/ 1.900, 1.810, 2.280, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 89)=BVS_Par_Type("HG+2", &
                   (/ 1.972, 2.170, 2.280, 2.380, 2.620, 2.308, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     5,     5,     5,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 90)=BVS_Par_Type("HG+9", &
                   (/ 0.000, 0.000, 0.000, 2.400, 2.590, 2.320, 2.470, 2.610, 2.020, 2.420, 2.500, 1.710, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 91)=BVS_Par_Type("HO+3", &
                   (/ 2.025, 1.908, 2.401, 2.520, 2.760, 2.490, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,    16,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 92)=BVS_Par_Type("HO+9", &
                   (/ 0.000, 0.000, 0.000, 2.550, 2.770, 2.480, 2.610, 2.800, 2.180, 2.560, 2.640, 1.880, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 93)=BVS_Par_Type("I+1 ", &
                   (/ 0.000, 2.320, 2.470, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 94)=BVS_Par_Type("I+3 ", &
                   (/ 2.020, 1.900, 2.390, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 95)=BVS_Par_Type("I+5 ", &
                   (/ 2.003, 1.840, 2.380, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 96)=BVS_Par_Type("I+7 ", &
                   (/ 1.930, 1.830, 2.310, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 97)=BVS_Par_Type("IN+1", &
                   (/ 0.000, 0.000, 2.560, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 98)=BVS_Par_Type("IN+3", &
                   (/ 1.902, 1.792, 2.280, 2.510, 2.630, 2.370, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     5,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 99)=BVS_Par_Type("IN+9", &
                   (/ 0.000, 0.000, 0.000, 2.410, 2.630, 2.360, 2.470, 2.690, 2.030, 2.430, 2.510, 1.720, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(100)=BVS_Par_Type("IR+4", &
                   (/ 1.870, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(101)=BVS_Par_Type("IR+5", &
                   (/ 1.916, 1.820, 2.300, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(102)=BVS_Par_Type("IR+9", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.660, 2.380, 2.510, 2.710, 2.060, 2.460, 2.540, 1.760, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(103)=BVS_Par_Type("K+1 ", &
                   (/ 2.132, 1.992, 2.519, 2.660, 2.880, 2.590, 2.720, 2.930, 2.260, 2.640, 2.830, 2.100, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(104)=BVS_Par_Type("KR+2", &
                   (/ 0.000, 1.890, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(105)=BVS_Par_Type("LA+3", &
                   (/ 2.172, 2.020, 2.545, 2.720, 2.930, 2.643, 2.740, 2.940, 2.340, 2.730, 2.800, 2.060, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,    16,     2,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(106)=BVS_Par_Type("LI+1", &
                   (/ 1.466, 1.360, 1.910, 2.020, 2.220, 1.940, 2.090, 2.300, 1.610, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     2,     2,     2,     2,     0,     0,     0,     0,     0 /) )
     BVS_Table(107)=BVS_Par_Type("LU+3", &
                   (/ 1.971, 1.876, 2.361, 2.500, 2.730, 2.430, 2.560, 2.750, 2.110, 2.510, 2.590, 1.820, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(108)=BVS_Par_Type("MG+2", &
                   (/ 1.693, 1.578, 2.080, 2.280, 2.460, 2.180, 2.320, 2.530, 1.850, 2.290, 2.380, 1.530, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(109)=BVS_Par_Type("MN+2", &
                   (/ 1.790, 1.698, 2.133, 2.340, 0.000, 2.220, 0.000, 0.000, 1.849, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     1,     5,     0,     5,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(110)=BVS_Par_Type("MN+3", &
                   (/ 1.760, 1.660, 2.140, 0.000, 0.000, 0.000, 0.000, 0.000, 1.837, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(111)=BVS_Par_Type("MN+4", &
                   (/ 1.753, 1.710, 2.130, 0.000, 0.000, 0.000, 0.000, 0.000, 1.822, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(112)=BVS_Par_Type("MN+6", &
                   (/ 1.790, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(113)=BVS_Par_Type("MN+7", &
                   (/ 1.827, 1.720, 2.170, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(114)=BVS_Par_Type("MN+9", &
                   (/ 1.754, 0.000, 0.000, 2.260, 2.490, 2.200, 0.000, 2.550, 1.870, 2.240, 2.360, 1.550, 0.000, 2.320 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     7,     0,     0,     2,     2,     2,     0,     2,     2,     2,     2,     2,     0,     2 /) )
     BVS_Table(115)=BVS_Par_Type("MO+3", &
                   (/ 1.834, 1.760, 2.220, 2.340, 0.000, 0.000, 0.000, 0.000, 1.960, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    13,     5,     5,     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(116)=BVS_Par_Type("MO+4", &
                   (/ 1.886, 1.800, 2.170, 0.000, 0.000, 2.235, 0.000, 0.000, 2.043, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    10,     5,     5,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(117)=BVS_Par_Type("MO+5", &
                   (/ 1.907, 0.000, 2.260, 0.000, 0.000, 2.288, 0.000, 0.000, 2.009, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    10,     0,     5,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(118)=BVS_Par_Type("MO+6", &
                   (/ 1.907, 1.810, 2.280, 0.000, 0.000, 2.331, 0.000, 0.000, 2.009, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(119)=BVS_Par_Type("MO+9", &
                   (/ 1.879, 0.000, 0.000, 2.430, 2.640, 2.350, 2.490, 2.690, 2.040, 2.440, 2.520, 1.730, 0.000, 0.000 /), &
                   (/ 0.300, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    26,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(120)=BVS_Par_Type("N+3 ", &
                   (/ 1.361, 1.370, 1.750, 0.000, 0.000, 1.730, 0.000, 0.000, 1.440, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(121)=BVS_Par_Type("N+5 ", &
                   (/ 1.432, 1.360, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(122)=BVS_Par_Type("NA+1", &
                   (/ 1.803, 1.677, 2.150, 2.330, 2.560, 2.300, 2.410, 2.640, 1.930, 2.360, 2.530, 1.680, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(123)=BVS_Par_Type("NB+3", &
                   (/ 1.910, 1.710, 2.200, 2.350, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(124)=BVS_Par_Type("NB+4", &
                   (/ 1.880, 1.900, 2.260, 2.620, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(125)=BVS_Par_Type("NB+5", &
                   (/ 1.911, 1.870, 2.270, 0.000, 2.770, 0.000, 0.000, 0.000, 2.010, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     5,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(126)=BVS_Par_Type("NB+9", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.680, 2.370, 2.510, 2.700, 2.060, 2.460, 2.540, 1.750, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(127)=BVS_Par_Type("ND+2", &
                   (/ 1.950, 0.000, 0.000, 0.000, 0.000, 2.600, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(128)=BVS_Par_Type("ND+3", &
                   (/ 2.105, 2.008, 2.492, 2.660, 2.870, 2.590, 2.710, 2.890, 2.300, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0,     0,     0,     0 /) )
     BVS_Table(129)=BVS_Par_Type("NH+1", &
                   (/ 2.226, 2.129, 2.619, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    19,    19,    19,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(130)=BVS_Par_Type("NI+2", &
                   (/ 1.654, 1.596, 2.020, 2.200, 2.400, 1.980, 0.000, 0.000, 1.700, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(131)=BVS_Par_Type("NI+3", &
                   (/ 1.686, 1.580, 0.000, 0.000, 0.000, 2.040, 0.000, 0.000, 1.731, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    32,     5,     0,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(132)=BVS_Par_Type("NI+4", &
                   (/ 1.780, 1.610, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(133)=BVS_Par_Type("NI+9", &
                   (/ 0.000, 0.000, 0.000, 2.160, 2.340, 2.040, 2.140, 2.430, 1.750, 2.170, 2.240, 1.400, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(134)=BVS_Par_Type("NP+3", &
                   (/ 0.000, 2.000, 2.480, 2.620, 2.850, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.400, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,    16,    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(135)=BVS_Par_Type("NP+4", &
                   (/ 2.180, 2.020, 2.460, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(136)=BVS_Par_Type("NP+5", &
                   (/ 2.090, 1.970, 2.420, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(137)=BVS_Par_Type("NP+6", &
                   (/ 2.070, 1.970, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(138)=BVS_Par_Type("NP+7", &
                   (/ 2.060, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(139)=BVS_Par_Type("O+2 ", &
                   (/ 1.500, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(140)=BVS_Par_Type("OS+4", &
                   (/ 1.811, 1.720, 2.190, 2.370, 0.000, 2.210, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(141)=BVS_Par_Type("OS+5", &
                   (/ 0.000, 1.810, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(142)=BVS_Par_Type("OS+6", &
                   (/ 2.030, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(143)=BVS_Par_Type("OS+8", &
                   (/ 1.920, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(144)=BVS_Par_Type("P+3 ", &
                   (/ 1.630, 1.530, 0.000, 0.000, 0.000, 2.120, 2.240, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(145)=BVS_Par_Type("P+4 ", &
                   (/ 1.640, 1.660, 0.000, 0.000, 0.000, 2.130, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(146)=BVS_Par_Type("P+5 ", &
                   (/ 1.617, 1.540, 2.020, 2.170, 0.000, 2.145, 0.000, 0.000, 1.704, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     5,     5,     0,     1,     0,     0,     1,     0,     0,     0,     0,     0 /) )
     BVS_Table(147)=BVS_Par_Type("P+9 ", &
                   (/ 0.000, 0.000, 0.000, 2.150, 2.400, 2.110, 2.260, 2.440, 1.730, 2.190, 2.250, 1.410, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(148)=BVS_Par_Type("PA+4", &
                   (/ 2.150, 2.020, 2.490, 2.660, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(149)=BVS_Par_Type("PA+5", &
                   (/ 2.090, 2.040, 2.450, 2.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(150)=BVS_Par_Type("PB+2", &
                   (/ 1.963, 2.030, 2.530, 2.680, 2.830, 2.541, 2.690, 0.000, 2.180, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.490, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    17,     2,     2,     5,     5,     1,     5,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(151)=BVS_Par_Type("PB+4", &
                   (/ 2.042, 1.940, 2.430, 3.040, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(152)=BVS_Par_Type("PB+9", &
                   (/ 0.000, 0.000, 0.000, 2.640, 2.780, 2.550, 2.670, 2.840, 2.220, 2.640, 2.720, 1.970, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(153)=BVS_Par_Type("PD+2", &
                   (/ 1.792, 1.740, 2.050, 2.200, 2.360, 2.090, 0.000, 0.000, 1.820, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(154)=BVS_Par_Type("PD+4", &
                   (/ 0.000, 1.660, 0.000, 0.000, 0.000, 2.300, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(155)=BVS_Par_Type("PD+9", &
                   (/ 0.000, 0.000, 0.000, 2.190, 2.380, 2.100, 2.220, 2.480, 1.810, 2.220, 2.300, 1.470, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(156)=BVS_Par_Type("PM+3", &
                   (/ 0.000, 1.960, 2.450, 2.590, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(157)=BVS_Par_Type("PO+4", &
                   (/ 2.190, 2.380, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(158)=BVS_Par_Type("PR+3", &
                   (/ 2.138, 2.022, 2.500, 2.670, 2.890, 2.600, 0.000, 2.900, 2.300, 2.680, 2.750, 2.020, 0.000, 2.720 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     0,     2,     2,     2,     2,     2,     0,     2 /) )
     BVS_Table(159)=BVS_Par_Type("PT+2", &
                   (/ 1.768, 1.680, 2.050, 2.200, 0.000, 2.160, 0.000, 0.000, 1.810, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(160)=BVS_Par_Type("PT+3", &
                   (/ 1.870, 0.000, 2.300, 2.470, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(161)=BVS_Par_Type("PT+4", &
                   (/ 1.879, 1.759, 2.170, 2.600, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(162)=BVS_Par_Type("PT+9", &
                   (/ 0.000, 0.000, 0.000, 2.180, 2.370, 2.080, 2.190, 2.450, 1.770, 2.190, 2.260, 1.400, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(163)=BVS_Par_Type("PU+3", &
                   (/ 2.110, 2.000, 2.480, 2.600, 2.840, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(164)=BVS_Par_Type("PU+4", &
                   (/ 2.090, 1.970, 2.440, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(165)=BVS_Par_Type("PU+5", &
                   (/ 2.110, 1.960, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(166)=BVS_Par_Type("PU+6", &
                   (/ 2.060, 1.960, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(167)=BVS_Par_Type("PU+7", &
                   (/ 2.050, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(168)=BVS_Par_Type("RB+1", &
                   (/ 2.263, 2.160, 2.652, 2.780, 3.010, 2.700, 2.810, 3.000, 2.620, 2.760, 2.870, 2.260, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     1,     2,     2,     2,     2,     2,     5,     2,     2,     2,     0,     0 /) )
     BVS_Table(169)=BVS_Par_Type("RE+1", &
                   (/ 0.000, 0.000, 2.620, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(170)=BVS_Par_Type("RE+3", &
                   (/ 1.900, 0.000, 2.230, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(171)=BVS_Par_Type("RE+4", &
                   (/ 0.000, 1.810, 2.230, 2.350, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(172)=BVS_Par_Type("RE+5", &
                   (/ 1.860, 0.000, 2.240, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(173)=BVS_Par_Type("RE+6", &
                   (/ 0.000, 1.790, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(174)=BVS_Par_Type("RE+7", &
                   (/ 1.970, 1.860, 2.230, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(175)=BVS_Par_Type("RE+9", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.610, 2.370, 2.500, 2.700, 2.060, 2.460, 2.540, 1.750, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(176)=BVS_Par_Type("RH+3", &
                   (/ 1.793, 1.710, 2.080, 2.270, 0.000, 0.000, 0.000, 0.000, 1.820, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     5,     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(177)=BVS_Par_Type("RH+4", &
                   (/ 0.000, 1.590, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(178)=BVS_Par_Type("RH+5", &
                   (/ 0.000, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(179)=BVS_Par_Type("RH+9", &
                   (/ 0.000, 0.000, 0.000, 2.250, 2.480, 2.150, 0.000, 2.550, 1.880, 2.290, 2.370, 1.550, 0.000, 2.330 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     0,     2,     2,     2,     2,     2,     0,     2 /) )
     BVS_Table(180)=BVS_Par_Type("RU+2", &
                   (/ 0.000, 1.840, 0.000, 0.000, 0.000, 0.000, 2.110, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(181)=BVS_Par_Type("RU+3", &
                   (/ 1.770, 2.120, 2.250, 0.000, 0.000, 2.200, 0.000, 0.000, 1.820, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     5,     5,     0,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(182)=BVS_Par_Type("RU+4", &
                   (/ 1.834, 1.740, 2.210, 0.000, 0.000, 2.210, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(183)=BVS_Par_Type("RU+5", &
                   (/ 1.900, 1.820, 2.230, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(184)=BVS_Par_Type("RU+6", &
                   (/ 1.870, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(185)=BVS_Par_Type("RU+7", &
                   (/ 1.990, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(186)=BVS_Par_Type("RU+9", &
                   (/ 0.000, 0.000, 0.000, 2.260, 2.480, 2.160, 2.330, 2.540, 1.880, 2.290, 2.360, 1.610, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(187)=BVS_Par_Type("S+2 ", &
                   (/ 1.740, 0.000, 0.000, 0.000, 0.000, 2.030, 0.000, 0.000, 1.682, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     5,     0,     0,     1,     0,     0,     0,     0,     0 /) )
     BVS_Table(188)=BVS_Par_Type("S+4 ", &
                   (/ 1.644, 1.600, 2.020, 0.000, 0.000, 0.000, 0.000, 0.000, 1.762, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     1,     0,     0,     0,     0,     0 /) )
     BVS_Table(189)=BVS_Par_Type("S+6 ", &
                   (/ 1.624, 1.560, 2.030, 0.000, 0.000, 0.000, 0.000, 0.000, 1.720, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(190)=BVS_Par_Type("S+9 ", &
                   (/ 0.000, 0.000, 0.000, 2.170, 2.360, 2.070, 2.210, 2.450, 1.740, 2.150, 2.250, 1.380, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(191)=BVS_Par_Type("SB+3", &
                   (/ 1.973, 1.883, 2.350, 2.510, 2.760, 2.474, 2.600, 0.000, 2.108, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     5,     1,     5,     0,     4,     0,     0,     0,     0,     0 /) )
     BVS_Table(192)=BVS_Par_Type("SB+5", &
                   (/ 1.942, 1.797, 2.300, 2.480, 0.000, 0.000, 0.000, 0.000, 1.990, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(193)=BVS_Par_Type("SB+9", &
                   (/ 0.000, 0.000, 0.000, 2.500, 2.720, 2.450, 2.570, 2.780, 2.120, 2.520, 2.600, 2.770, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(194)=BVS_Par_Type("SC+3", &
                   (/ 1.849, 1.760, 2.360, 2.380, 2.590, 2.321, 2.440, 2.640, 1.980, 2.400, 2.480, 1.680, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     5,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(195)=BVS_Par_Type("SE+2", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 2.210, 2.330, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(196)=BVS_Par_Type("SE+4", &
                   (/ 1.811, 1.730, 2.220, 2.430, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(197)=BVS_Par_Type("SE+6", &
                   (/ 1.788, 1.690, 2.160, 0.000, 0.000, 0.000, 0.000, 0.000, 1.900, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(198)=BVS_Par_Type("SE+9", &
                   (/ 0.000, 0.000, 0.000, 2.330, 2.540, 2.250, 2.360, 2.550, 0.000, 2.340, 2.420, 1.540, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     0,     2,     2,     2,     0,     0 /) )
     BVS_Table(199)=BVS_Par_Type("SI+4", &
                   (/ 1.624, 1.580, 2.030, 2.200, 2.410, 2.126, 2.260, 2.490, 1.724, 2.230, 2.310, 1.470, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     2,     2,     1,     2,     2,     1,     2,     2,     2,     0,     0 /) )
     BVS_Table(200)=BVS_Par_Type("SM+3", &
                   (/ 2.088, 1.940, 1.977, 2.660, 2.840, 2.550, 2.670, 2.860, 2.240, 2.630, 2.700, 1.960, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,    16,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(201)=BVS_Par_Type("SN+2", &
                   (/ 1.940, 1.925, 2.410, 2.530, 2.810, 2.440, 0.000, 0.000, 2.030, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     1,     5,     4,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(202)=BVS_Par_Type("SN+3", &
                   (/ 0.000, 0.000, 2.360, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(203)=BVS_Par_Type("SN+4", &
                   (/ 1.905, 1.843, 2.276, 2.400, 0.000, 2.399, 2.510, 0.000, 2.030, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     1,     5,     0,     1,     5,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(204)=BVS_Par_Type("SN+9", &
                   (/ 0.000, 0.000, 0.000, 2.550, 2.760, 2.390, 2.590, 2.760, 2.060, 2.450, 2.620, 1.850, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,    27,     2,     2,    27,     2,     2,     2,     0,     0 /) )
     BVS_Table(205)=BVS_Par_Type("SR+2", &
                   (/ 2.118, 2.019, 2.510, 2.680, 2.880, 2.590, 2.720, 2.870, 2.230, 2.670, 2.760, 2.010, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(206)=BVS_Par_Type("TA+4", &
                   (/ 2.290, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(207)=BVS_Par_Type("TA+5", &
                   (/ 1.920, 1.880, 2.300, 0.000, 0.000, 2.470, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(208)=BVS_Par_Type("TA+9", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.660, 2.390, 2.510, 2.700, 2.010, 2.470, 2.550, 1.760, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(209)=BVS_Par_Type("TB+3", &
                   (/ 2.032, 1.936, 2.427, 2.580, 2.800, 2.510, 2.630, 2.820, 2.200, 2.590, 2.660, 1.910, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(210)=BVS_Par_Type("TC+4", &
                   (/ 0.000, 1.880, 2.210, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(211)=BVS_Par_Type("TC+7", &
                   (/ 1.900, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(212)=BVS_Par_Type("TE+4", &
                   (/ 1.977, 1.870, 2.370, 2.550, 2.787, 2.440, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(213)=BVS_Par_Type("TE+6", &
                   (/ 1.917, 1.820, 2.300, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(214)=BVS_Par_Type("TE+9", &
                   (/ 0.000, 0.000, 0.000, 2.530, 2.760, 2.450, 2.530, 2.760, 2.120, 2.520, 2.600, 1.830, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(215)=BVS_Par_Type("TH+4", &
                   (/ 2.167, 2.068, 2.550, 2.710, 2.930, 2.640, 2.760, 2.940, 2.340, 2.730, 2.800, 2.070, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(216)=BVS_Par_Type("TI+2", &
                   (/ 0.000, 2.150, 2.310, 2.490, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(217)=BVS_Par_Type("TI+3", &
                   (/ 1.791, 1.723, 2.220, 0.000, 2.520, 2.110, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     5,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(218)=BVS_Par_Type("TI+4", &
                   (/ 1.815, 1.760, 2.190, 2.360, 0.000, 2.290, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(219)=BVS_Par_Type("TI+9", &
                   (/ 1.790, 0.000, 2.184, 2.320, 2.540, 2.240, 2.380, 2.600, 1.930, 2.360, 2.420, 1.610, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    11,     0,    11,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(220)=BVS_Par_Type("TL+1", &
                   (/ 2.124, 2.150, 2.560, 2.690, 2.822, 2.545, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     1,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(221)=BVS_Par_Type("TL+3", &
                   (/ 2.003, 1.880, 2.320, 2.650, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(222)=BVS_Par_Type("TL+9", &
                   (/ 0.000, 0.000, 0.000, 2.700, 2.910, 2.630, 2.700, 2.930, 2.290, 2.710, 2.790, 2.050, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(223)=BVS_Par_Type("TM+3", &
                   (/ 2.000, 1.842, 2.380, 2.530, 2.740, 2.450, 2.580, 2.770, 2.140, 2.530, 2.620, 1.850, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(224)=BVS_Par_Type("U+2 ", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.080, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     5,     0 /) )
     BVS_Table(225)=BVS_Par_Type("U+3 ", &
                   (/ 0.000, 2.020, 2.490, 2.640, 2.870, 2.540, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.400, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,    16,    16,    16,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(226)=BVS_Par_Type("U+4 ", &
                   (/ 2.112, 2.038, 2.470, 2.600, 2.880, 2.550, 0.000, 0.000, 2.180, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     1,    16,    16,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(227)=BVS_Par_Type("U+5 ", &
                   (/ 2.075, 1.966, 2.460, 2.700, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(228)=BVS_Par_Type("U+6 ", &
                   (/ 2.051, 1.980, 2.420, 0.000, 0.000, 0.000, 0.000, 0.000, 1.930, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.519, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    18,    16,    16,     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(229)=BVS_Par_Type("U+9 ", &
                   (/ 0.000, 0.000, 0.000, 2.630, 2.840, 2.560, 2.700, 2.860, 2.240, 2.640, 2.720, 1.970, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(230)=BVS_Par_Type("V+1 ", &
                   (/ 1.880, 0.000, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(231)=BVS_Par_Type("V+2 ", &
                   (/ 1.700, 2.160, 2.440, 0.000, 0.000, 2.110, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(232)=BVS_Par_Type("V+3 ", &
                   (/ 1.743, 1.702, 2.190, 2.330, 0.000, 2.170, 0.000, 0.000, 1.813, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     5,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(233)=BVS_Par_Type("V+4 ", &
                   (/ 1.784, 1.700, 2.160, 0.000, 0.000, 2.226, 0.000, 0.000, 1.875, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(234)=BVS_Par_Type("V+5 ", &
                   (/ 1.803, 1.700, 2.160, 0.000, 0.000, 2.250, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     2,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(235)=BVS_Par_Type("V+9 ", &
                   (/ 1.810, 0.000, 0.000, 2.300, 2.510, 2.230, 2.330, 2.570, 1.860, 2.310, 2.390, 1.580, 0.000, 0.000 /), &
                   (/ 0.340, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(236)=BVS_Par_Type("W+5 ", &
                   (/ 1.890, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(237)=BVS_Par_Type("W+6 ", &
                   (/ 1.917, 1.830, 2.270, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(238)=BVS_Par_Type("W+9 ", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.660, 2.390, 2.510, 2.710, 2.060, 2.460, 2.540, 1.760, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(239)=BVS_Par_Type("XE+2", &
                   (/ 2.050, 2.020, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(240)=BVS_Par_Type("XE+4", &
                   (/ 0.000, 1.930, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(241)=BVS_Par_Type("XE+6", &
                   (/ 2.000, 1.890, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(242)=BVS_Par_Type("XE+8", &
                   (/ 1.940, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(243)=BVS_Par_Type("Y+3 ", &
                   (/ 2.019, 1.904, 2.400, 2.550, 2.770, 2.480, 2.610, 2.800, 2.170, 2.570, 2.640, 1.860, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(244)=BVS_Par_Type("YB+3", &
                   (/ 1.965, 1.875, 2.371, 2.451, 2.720, 2.430, 2.560, 2.760, 2.120, 2.530, 2.590, 1.820, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(245)=BVS_Par_Type("ZN+2", &
                   (/ 1.704, 1.620, 2.010, 2.150, 2.360, 2.090, 2.220, 2.450, 1.720, 2.150, 2.240, 1.420, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     5,     2,     2,     2,     0,     0 /) )
     BVS_Table(246)=BVS_Par_Type("ZR+2", &
                   (/ 2.340, 2.240, 2.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(247)=BVS_Par_Type("ZR+4", &
                   (/ 1.928, 1.846, 2.330, 2.480, 2.690, 2.410, 2.530, 2.670, 2.110, 2.520, 2.570, 1.790, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     return
    End Subroutine Set_BVS_Table


    !!----
    !!----
    !!---- Subroutine Set_SBVS_Table()
    !!----
    !!----    Fills the parameters for soft-BVS from Stefan Adams
    !!----    Acta Cryst B57, 278-287 (2001)
    !!----
    !!---- Update: December - 2014
    !!
    Subroutine Set_SBVS_Table()

       if (.not. allocated(sBVS_Table)) allocate(sBVS_Table(sbvs_species_n))
                     !  O      F     Cl     Br      I      S     Se     Te
     sBVS_Table(  1)=BVS_Par_Type("LI+1", &
                   (/1.1725,1.1011,1.3418,1.5340,1.6733,1.5070,1.5296,1.7340, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/ 0.515, 0.501, 0.661, 0.665, 0.723, 0.632, 0.735, 0.717, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/    34,    34,    34,    34,    34,    34,    34,    34,     0,     0,     0,     0,     0,     0 /) )
     sBVS_Table(  2)=BVS_Par_Type("NA+1", &
                   (/1.5602,1.4262,1.6940,0.0000,1.9694,1.8311,1.8787,2.0518, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/ 0.483, 0.475, 0.603, 0.000, 0.688, 0.621, 0.660, 0.684, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/    34,    34,    34,    34,    34,    34,    34,    34,     0,     0,     0,     0,     0,     0 /) )
     sBVS_Table(  3)=BVS_Par_Type("K+1", &
                   (/1.9729,1.8472,2.0866,2.1001,2.3202,2.1711,2.2569,2.3926, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/ 0.422, 0.422, 0.552, 0.625, 0.641, 0.571, 0.624, 0.662, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/    34,    34,    34,    34,    34,    34,    34,    34,     0,     0,     0,     0,     0,     0 /) )
     sBVS_Table(  4)=BVS_Par_Type("RB+1", &
                   (/2.0573,1.9572,2.2443,2.3272,2.4667,2.3011,2.4015,2.4600, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/ 0.425, 0.418, 0.540, 0.579, 0.631, 0.552, 0.581, 0.615, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/    34,    34,    34,    34,    34,    34,    34,    34,     0,     0,     0,     0,     0,     0 /) )
     sBVS_Table(  5)=BVS_Par_Type("CS+1", &
                   (/2.2985,2.1955,2.5046,2.5152,2.6951,2.5147,2.6568,2.7360, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/ 0.403, 0.411, 0.481, 0.538, 0.608, 0.522, 0.546, 0.617, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /),&
                   (/    34,    34,    34,    34,    34,    34,    34,    34,     0,     0,     0,     0,     0,     0 /) )
    End Subroutine Set_SBVS_Table



    !!----
    !!---- Subroutine Set_Table_BVEL_Params(A,N_bvel,bvel)
    !!---- type (Atoms_Conf_List_type),            intent(in)  :: A
    !!---- integer,                      optional, intent(in)  :: N_bvel   !Number of bvel strings with externally provided values
    !!---- character(len=*),dimension(:),optional, intent(in)  :: bvel     !bvs strings with externally provided values
    !!----
    !!----
    !!----
    !!---- Created: December - 2014
    !!
    Subroutine Set_Table_BVEL_Params(A,N_bvel,bvel)
    type (Atoms_Conf_List_type),            intent(in)  :: A
    integer,                      optional, intent(in)  :: N_bvel
    character(len=*),dimension(:),optional, intent(in)  :: bvel

       !---- Local Variables ----!
       integer :: i,j,k,ia,ic

       if (A%N_Spec == 0) then
          err_conf=.true.
          ERR_Conf_Mess=" The number of different species is zero, tables cannot be set"
          return
       end if

       if (allocated(Table_Avcoor))   deallocate(Table_Avcoor)
       if (allocated(Table_Rzero))    deallocate(Table_Rzero)
       if (allocated(Table_Rcutoff))  deallocate(Table_Rcutoff)
       if (allocated(Table_Dzero))    deallocate(Table_Dzero)
       if (allocated(Table_Rmin))     deallocate(Table_Rmin)
       if (allocated(Table_Alpha))    deallocate(Table_Alpha)
       if (allocated(Table_Ref))      deallocate(Table_Ref)

       allocate(Table_Avcoor(A%N_Spec,A%N_Spec))  ; Table_Avcoor  =0.0
       allocate(Table_Rzero(A%N_Spec,A%N_Spec))   ; Table_Rzero   =0.0
       allocate(Table_Rcutoff(A%N_Spec,A%N_Spec)) ; Table_Rcutoff =0.0
       allocate(Table_Dzero(A%N_Spec,A%N_Spec))   ; Table_Dzero   =0.0
       allocate(Table_Rmin(A%N_Spec,A%N_Spec))    ; Table_Rmin    =0.0
       allocate(Table_Alpha(A%N_Spec,A%N_Spec))   ; Table_Alpha   =0.0
       allocate(Table_ref(A%N_Spec,A%N_Spec))     ; Table_ref     =0

       call Set_BVEL_Table()

       do i=1,A%N_Cations
          ic=0
          do j=1,bvel_species_n
             if (A%Species(i) == BVEL_Table(j)%Symb) then
                ic=j
                exit
             end if
          end do
          if (ic == 0) then
             if(.not. present(N_bvel)) then
               err_conf=.true.
               ERR_Conf_Mess=" Cation not found on the internal list: "//A%Species(i)
               return
             else
                call Complete_Table_BVEL(A,N_bvel,bvel)
                if(err_conf) then
                     return
                else
                     cycle
                end if
             end if
          end if

          do k=1,A%N_Anions
             ia=0
             do j=1,bvs_anions_n
                if (A%Species(A%N_Cations+k) == bvel_anions(j)) then
                   ia=j
                   exit
                end if
             end do
             if (ia == 0) then
                err_conf=.true.
                ERR_Conf_Mess=" Anion not found on the internal list: "//A%Species(A%N_Cations+k)
                return
             end if
             Table_Avcoor (i,A%N_Cations+k)=bvel_table(ic)%Avcoor(ia)
             Table_Rzero  (i,A%N_Cations+k)=bvel_table(ic)%Rzero(ia)
             Table_Rcutoff(i,A%N_Cations+k)=bvel_table(ic)%Rcutoff(ia)
             Table_Dzero  (i,A%N_Cations+k)=bvel_table(ic)%Dzero(ia)
             Table_Rmin   (i,A%N_Cations+k)=bvel_table(ic)%Rmin(ia)
             Table_Alpha  (i,A%N_Cations+k)=bvel_table(ic)%Alpha(ia)
             Table_ref    (i,A%N_Cations+k)=bvel_table(ic)%refnum(ia)

             Table_Avcoor (A%N_Cations+k,i)=bvel_table(ic)%Avcoor(ia)
             Table_Rzero  (A%N_Cations+k,i)=bvel_table(ic)%Rzero(ia)
             Table_Rcutoff(A%N_Cations+k,i)=bvel_table(ic)%Rcutoff(ia)
             Table_Dzero  (A%N_Cations+k,i)=bvel_table(ic)%Dzero(ia)
             Table_Rmin   (A%N_Cations+k,i)=bvel_table(ic)%Rmin(ia)
             Table_Alpha  (A%N_Cations+k,i)=bvel_table(ic)%Alpha(ia)
             Table_ref    (A%N_Cations+k,i)=bvel_table(ic)%refnum(ia)

          end do
       end do

       call Deallocate_BVEL_Table()

       return
    End Subroutine Set_Table_BVEL_Params

    !!----
    !!---- Subroutine Set_Table_d0_b(A,N_bvsm,bvs_m)
    !!---- type (Atoms_Conf_List_type),            intent(in)  :: A
    !!---- integer,                      optional, intent(in)  :: N_bvsm   !Number of bvs strings with externally provided values
    !!---- character(len=*),dimension(:),optional, intent(in)  :: bvs_m    !bvs strings with externally provided values
    !!----
    !!----
    !!----
    !!---- Updated: March - 2005, December-2014
    !!
    Subroutine Set_Table_d0_b(A,N_bvsm,bvs_m,soft)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),            intent(in)  :: A
       integer,                      optional, intent(in)  :: N_bvsm
       character(len=*),dimension(:),optional, intent(in)  :: bvs_m
       logical,                      optional, intent(in)  :: soft
       !---- Local Variables ----!
       integer :: i,j,k,ia,ic
       logical :: found

       if (A%N_Spec == 0) then
          err_conf=.true.
          ERR_Conf_Mess=" The number of different species is zero, tables cannot be set"
          return
       end if

       if (allocated(Table_d0))   deallocate(Table_d0)
       if (allocated(Table_b))    deallocate(Table_b)
       if (allocated(Table_ref))  deallocate(Table_ref)

       allocate(Table_d0(A%N_Spec,A%N_Spec))
       allocate(Table_b(A%N_Spec,A%N_Spec))
       allocate(Table_ref(A%N_Spec,A%N_Spec))

       Table_d0=0.0
       Table_b=0.37
       Table_ref = 0

       call Set_BVS_Table()
       if(present(soft)) call set_SBVS_Table()

       do i=1,A%N_Cations
          ic=0; found=.false.
          if(present(soft)) then
             do j=1,sbvs_species_n
                if (A%Species(i) == sBVS_Table(j)%Symb) then
                   ic=j
                   found=.true.
                   exit
                end if
             end do
          end if
          if(ic == 0) then
            do j=1,bvs_species_n
               if (A%Species(i) == BVS_Table(j)%Symb) then
                  ic=j
                  exit
               end if
            end do
          end if
          if (ic == 0) then
             if(.not. present(N_bvsm)) then
               err_conf=.true.
               ERR_Conf_Mess=" Cation not found on the internal list: "//A%Species(i)
               return
             else
                call Complete_Table(A,N_bvsm,bvs_m)
                if(err_conf) then
                     return
                else
                     cycle
                end if
             end if
          end if

          do k=1,A%N_Anions
             ia=0
             do j=1,bvs_anions_n
                if (A%Species(A%N_Cations+k) == bvs_anions(j)) then
                   ia=j
                   exit
                end if
             end do
             if (ia == 0) then
                err_conf=.true.
                ERR_Conf_Mess=" Anion not found on the internal list: "//A%Species(A%N_Cations+k)
                return
             end if
             if(found) then
                Table_d0 (i,A%N_Cations+k)=sbvs_table(ic)%d0(ia)
                Table_b  (i,A%N_Cations+k)=sbvs_table(ic)%b_par(ia)
                Table_ref(i,A%N_Cations+k)=sbvs_table(ic)%refnum(ia)

                Table_d0 (A%N_Cations+k,i)=sbvs_table(ic)%d0(ia)
                Table_b  (A%N_Cations+k,i)=sbvs_table(ic)%b_par(ia)
                Table_ref(A%N_Cations+k,i)=sbvs_table(ic)%refnum(ia)
             else
                Table_d0 (i,A%N_Cations+k)=bvs_table(ic)%d0(ia)
                Table_b  (i,A%N_Cations+k)=bvs_table(ic)%b_par(ia)
                Table_ref(i,A%N_Cations+k)=bvs_table(ic)%refnum(ia)

                Table_d0 (A%N_Cations+k,i)=bvs_table(ic)%d0(ia)
                Table_b  (A%N_Cations+k,i)=bvs_table(ic)%b_par(ia)
                Table_ref(A%N_Cations+k,i)=bvs_table(ic)%refnum(ia)
             end if

          end do
       end do

       call Deallocate_BVS_Table()
       if(present(soft)) call Deallocate_BVS_Table()

       return
    End Subroutine Set_Table_d0_b


    !!----
    !!---- Subroutine Species_on_List(A,MulG,tol, covalent)
    !!----    type (Atoms_Conf_List_Type), intent(in out) :: A
    !!----    Integer, optional,           intent(in)     :: MulG
    !!----    real(kind=cp), optional,     intent(in)     :: tol
    !!----    logical,       optional,     intent(in)     :: covalent
    !!----
    !!----    Determines the different species in the List and,
    !!----    optionally, sets the tolerance factor for ionic radii
    !!----    conditions and provides "corrected" occupation factors
    !!----    (mult/MulG) when the user is using a multiplier. The
    !!----    general multiplicity of the space group MulG must be
    !!----    provided in such a case. This first free variable of the
    !!----    Atom-type A%Atom%VFree(1) is set to the corrected
    !!----    occupation. The first atom in the list must completely
    !!----    occupy its site. If covalent is provided the covalent
    !!----    radius is used instead of the ionic radius.
    !!----
    !!---- Update: March - 2005, December 2014
    !!
    Subroutine Species_on_List(A,MulG, tol, covalent)
       !---- Arguments ----!
       type (Atoms_Conf_List_Type), intent(in out) :: A
       Integer, optional,           intent(in)     :: MulG
       real(kind=cp), optional,     intent(in)     :: tol
       logical,       optional,     intent(in)     :: covalent

       !---- Local variables ----!
       character(len=4), dimension(50) :: cation,anion,spec
       character(len=2)                :: car,cv
       character(len=4)                :: canio
       integer                         :: i,im,j,v,ns,nc,na
       real(kind=cp)                   :: fac1,fact


       if (A%natoms == 0) return

       ns=0
       spec  = " "
       nc=0
       cation= " "
       na=0
       anion = " "

       if(present(tol)) A%tol=tol

       if(present(MulG)) then

         fac1=A%atom(1)%Occ*real(MulG)/real(A%atom(1)%mult)
         fac1=1.0/fac1
         A%totatoms=0.0
         do i=1,a%natoms
            fact=real(MulG)/real(a%atom(i)%mult)
            A%Atom(i)%VarF(1)=A%atom(i)%occ*fact*fac1      !Site Occupancy (=1, full occupation)
            A%Atom(i)%VarF(2)=A%atom(i)%occ_std*fact*fac1  !standard deviation of Site Occupancy
            A%totatoms=A%totatoms + A%Atom(i)%VarF(1)*real(a%atom(i)%mult) !total number of atoms/conventional unit cell
         end do

       else

         A%totatoms=0.0
         do i=1,a%natoms
            A%Atom(i)%VarF(1)=A%atom(i)%occ  !The user has given site occupancy
            A%Atom(i)%VarF(2)=A%atom(i)%occ_std
            A%totatoms=A%totatoms + A%Atom(i)%VarF(1)*real(a%atom(i)%mult) !total number of atoms/conventional unit cell
         end do

       end if

       loop1:do i=1, A%natoms
          car=u_case(a%atom(i)%ChemSymb)
          v=nint(a%atom(i)%charge)
          if (v == 0) then
             err_conf=.true.
             ERR_Conf_Mess=" The Atom "//a%atom(i)%lab//"has not charge"
             return
          end if
          write(unit=cv,fmt="(i2)") v
          if (v > 0) cv(1:1)="+"
          canio=car//cv
          canio=pack_string(canio)

          if (v > 0) then
             do j=1,nc
                if (canio == cation(j)) cycle loop1
             end do
             nc=nc+1
             cation(nc)=canio
          else
             do j=1,na
                if (canio == anion(j)) cycle loop1
             end do
             na=na+1
             anion(na)=canio
          end if
          if (na+nc == 50) exit
       end do loop1

       ns=nc+na
       A%N_Spec    = ns
       A%N_Anions  = na
       A%N_Cations = nc

       !---- Order the Species vector ----!
       call sort_strings(cation(1:nc))
       call sort_strings(anion(1:na))
       spec(1:nc)=cation(1:nc)
       spec(nc+1:nc+na)=anion(1:na)

       if (allocated(A%Species)) deallocate(A%Species)
       allocate(A%Species(ns))
       A%Species=spec(1:ns)

       if (allocated(A%Radius)) deallocate(A%Radius)
       allocate(A%Radius(ns))

       do i=1,nc
          im=index(A%Species(i),"+")
          car=A%Species(i)(im+1:im+2)
          read(unit=car,fmt="(i1)") j
          car=A%Species(i)(1:im-1)
          if(present(covalent)) then
             call get_covalent_radius(car,A%Radius(i))
          else
             call get_ionic_radius(car,j,A%Radius(i))
          end if
          if (A%Radius(i) < 0.01) A%Radius(i)=0.8
       end do

       do i=1,A%N_Anions
          do j=1,bvs_anions_n
             if (A%Species(nc+i) == bvs_anions(j)) then
                if(present(covalent)) then
                    call get_covalent_radius(A%Species(nc+i),A%Radius(nc+i))
                else
                    A%Radius(nc+i) = bvs_anions_rion(j)
                end if
                exit
             end if
          end do
       end do

       !---- Fix the index on Atom_type pointing to the Species vector ----!
       do i=1, A%natoms
          do j=1,ns
             im=index(A%Species(j),"+")
             if (im == 0) im=index(A%Species(j),"-")
             car=A%Species(j)(1:im-1)
             cv=A%Species(j)(im:im+1)
             if (cv(1:1)=="+") cv(1:1)=" "
                read(unit=cv,fmt="(i2)") v
                if (u_case(A%Atom(i)%ChemSymb) == car .and. nint(A%Atom(i)%charge) == v) then
                A%atom(i)%ind(1)=j
                exit
             end if
          end do
       end do

       return
    End Subroutine Species_on_List

 End Module CFML_BVS_Energy_Calc
