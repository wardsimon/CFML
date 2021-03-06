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
!!----    AP_SPECIES_N
!!----    AP_TABLE
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
!!----    SBVS_PAR_TYPE
!!----    SBVS_SPECIES_N
!!----    SBVS_TABLE
!!----    TABLE_Alpha
!!----    TABLE_Avcoor
!!----    TABLE_B
!!----    TABLE_D0
!!----    TABLE_Dzero
!!----    TABLE_Rcutoff
!!----    TABLE_Ref
!!----    TABLE_Rmin
!!----    TABLE_Rzero
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       DEALLOCATE_AP_TABLE
!!----       DEALLOCATE_BVEL_TABLE
!!----       DEALLOCATE_BVS_TABLE
!!----       SET_ATOMIC_PROPERTIES
!!----       SET_BVEL_TABLE
!!----       SET_BVS_TABLE
!!----       SET_SBVS_TABLE
!!----       SET_TABLE_D0_B
!!----       SET_COMMON_OXIDATION_STATES_TABLE
!!----       SET_OXIDATION_STATES_TABLE
!!----       SET_PAULING_ELECTRONEGATIVITY
!!----
!!
   Module CFML_BVSpar
    !---- Use Files ----!
    Use CFML_GlobalDeps,  only: Cp,Dp

    !---- Variables ----!
    implicit none

    private

    !---- List of public subroutines ----!

    public :: Deallocate_Bvs_Table, Deallocate_BVEL_Table, Deallocate_Ap_Table, Deallocate_sBVS_Table,&
              Set_Atomic_Properties, Set_Bvel_Table, Set_Bvs_Table, Set_Sbvs_Table, Set_Pauling_Electronegativity, &
              Set_Common_Oxidation_States_Table, Set_Oxidation_States_Table

    !!----
    !!---- Ap_species_n
    !!----    integer, parameter, public :: Ap_species_n=183
    !!----
    !!----    Number of atomic properties species in:  Ap_Table
    !!----
    !!---- Created: January - 2015
    !!
    integer, parameter, public :: Ap_species_n=183

    !!----
    !!---- TYPE :: Atomic_Properties_Type
    !!--..
    !!---- Type, public :: Atomic_Properties_Type
    !!----     integer          :: Z        ! Atomic number
    !!----     character(len=2) :: Symb     ! Element
    !!----     integer          :: oxs      ! Nominal oxidation state
    !!----     integer          :: dox      ! Default oxidation state
    !!----     real             :: Mass     ! Atomic mass in atomic units
    !!----     integer          :: n        ! Principal Quantum number (period)
    !!----     integer          :: g        ! Group in the periodic table
    !!----     integer          :: b        ! Block (s:0, p:1, d:2, f:3)
    !!----     real             :: Rc       ! Covalent radius
    !!----     real             :: sigma    ! Softness
    !!---- End Type Atomic_Properties_Type
    !!----
    !!----    Type Definition for single atomic properties
    !!----
    !!---- Created: January - 2015
    !!
    Type, public :: Atomic_Properties_Type
        integer          :: Z        ! Atomic number
        character(len=4) :: Symb     ! Element with charge
        integer          :: oxs      ! Nominal oxidation state
        integer          :: dox      ! Default oxidation state
        real             :: Mass     ! Atomic mass in atomic units
        integer          :: n        ! Principal Quantum number (period)
        integer          :: g        ! Group in the periodic table
        integer          :: b        ! Block (s:0, p:1, d:2, f:3)
        real             :: Rc       ! Covalent radius
        real             :: sigma    ! Softness
    End Type Atomic_Properties_Type

    Type(Atomic_Properties_Type), allocatable,  dimension(:),public :: Ap_Table

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
    !!----    real(kind=cp), parameter, dimension(bvel_anions_n) :: bvel_anions_rion
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

    !!---- TYPE :: sBVS_PAR_TYPE
    !!--..
    !!---- Type, public :: sBvs_Par_Type
    !!----    character (len=4)                      :: Symb      ! Chemical symbol
    !!----    real(kind=cp), dimension(bvs_anions_n) :: D0        ! D0 Parameter
    !!----    real(kind=cp), dimension(bvs_anions_n) :: B_Par     ! B Parameter
    !!----    real(kind=cp),dimension(bvs_anions_n)  :: cn        ! Preferred Coordination
    !!----    real(kind=cp),dimension(bvs_anions_n)  :: ctoff     ! Cutoff distance
    !!----    integer,       dimension(bvs_anions_n) :: refnum    ! Integer pointing to the reference paper
    !!---- End Type sBvs_Par_Type
    !!----
    !!----    Definition for sBVS Parameters
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: sBvs_Par_Type
       character (len=4)                     :: Symb
       real(kind=cp),dimension(bvs_anions_n) :: d0
       real(kind=cp),dimension(bvs_anions_n) :: b_par
       real(kind=cp),dimension(bvs_anions_n) :: cn
       real(kind=cp),dimension(bvs_anions_n) :: ctoff
       integer      ,dimension(bvs_anions_n) :: refnum
    End Type sBvs_Par_Type

    !!----
    !!---- SBVS_SPECIES_N
    !!----    integer, parameter, public :: sbvs_species_n=5 !only alkali chalcogenides
    !!----
    !!----    Maximum Number of species in SBVS_Table
    !!----
    !!---- Created: December - 2014, Updated: January 2015
    !!
    integer, parameter, public :: sbvs_species_n=168

    !!----
    !!---- SBVS_TABLE
    !!----    Type(sBvs_Par_Type), allocatable, dimension(:), public :: sBVS_Table
    !!----
    !!----    SBVS Parameters for calculations (only alkali chalcogenides are available)
    !!----
    !!---- Created: December - 2014
    !!
    Type(sBvs_Par_Type), allocatable, dimension(:), public :: sBVS_Table

    !!----
    !!---- Table_Alpha
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Alpha
    !!----
    !!----    Matrix N_Species x N_Species of Alpha (equivalent to 1/b in BVS) parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, public :: Table_Alpha

    !!----
    !!---- Table_Avcoor
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Avcoor
    !!----
    !!----    Matrix N_Species x N_Species of Average coordination parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, public :: Table_Avcoor

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
    real(kind=cp),dimension(:,:), allocatable, public :: Table_Dzero

    !!----
    !!---- Table_Rcutoff
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rcutoff
    !!----
    !!----    Matrix N_Species x N_Species of Rcutoff parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, public :: Table_Rcutoff

    !!----
    !!---- TABLE_Ref
    !!----    integer,dimension(:,:), allocatable, private :: Table_ref
    !!----
    !!----    Matrix N_Species x N_Species with references for BVS parameters
    !!----
    !!---- Update: March - 2005
    !!
    integer,dimension(:,:), allocatable, public :: Table_ref

    !!----
    !!---- Table_Rmin
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rmin
    !!----
    !!----    Matrix N_Species x N_Species of Rmin parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, public :: Table_Rmin

     !!----
    !!---- Table_Rzero
    !!----    real(kind=cp),dimension(:,:), allocatable, private :: Table_Rzero
    !!----
    !!----    Matrix N_Species x N_Species of Rzero (equivalent to D0 in BVS) parameters for BVEL
    !!----
    !!---- Created: December - 2014
    !!
    real(kind=cp),dimension(:,:), allocatable, public :: Table_Rzero

    !!----
    !!----  Reference list for BVS parameters
    !!----
    !!----
    !!----    List of Reference for BVS Data
    !!----
    !!---- Update: March - 2005
    !!
    character(len=*),dimension(0:35),parameter, public :: references = (/  &
         "Unknown  (provided by the user)                                                 ", &  !0
         "Brown and Altermatt, (1985), Acta Cryst. B41, 244-247 (empirical)               ", &  !1
         "Brese and O'Keeffe, (1991), Acta Cryst. B47, 192-197 (extrapolated)             ", &  !2
         "Adams, 2001, Acta Cryst. B57, 278-287 (includes second neighbours)              ", &  !3
         "Hu et al. (1995) Inorg. Chim. Acta, 232, 161-165.                               ", &  !4
         "I.D.Brown Private communication                                                 ", &  !5
         "Brown et al. (1984) Inorg. Chem. 23, 4506-4508                                  ", &  !6
         "Palenik (1997) Inorg. Chem. 36 4888-4890                                        ", &  !7
         "Kanowitz and Palenik (1998) Inorg. Chem. 37 2086-2088                           ", &  !8
         "Wood and Palenik (1998) Inorg. Chem. 37 4149-4151                               ", &  !9
         "Liu and Thorp (1993) Inorg. Chem. 32 4102-4105                                  ", &  !10
         "Palenik (1997) Inorg. Chem. 36 3394-3397                                        ", &  !11
         "Shields, Raithby, Allen and Motherwell (1999) Acta Cryst.B56, 455-465           ", &  !12
         "Chen, Zhou and Hu (2002) Chinese Sci. Bul. 47, 978-980.                         ", &  !13
         "Kihlbourg (1963) Ark. Kemi 21 471; Schroeder 1975 Acta Cryst. B31, 2294         ", &  !14
         "Allmann (1975) Monatshefte Chem. 106, 779                                       ", &  !15
         "Zachariesen (1978) J.Less Common Metals 62, 1                                   ", &  !16
         "Krivovichev and Brown (2001) Z. Krist. 216, 245                                 ", &  !17
         "Burns, Ewing and Hawthorne (1997) Can. Miner. 35,1551-1570                      ", &  !18
         "Garcia-Rodriguez, et al. (2000) Acta Cryst. B56, 565-569                        ", &  !19
         "Mahapatra et al. (1996) J. Amer.Chem. Soc. 118, 11555                           ", &  !20
         "Wood and Palenik (1999) Inorg. Chem. 38, 1031-1034                              ", &  !21
         "Wood and Palenik (1999) Inorg. Chem. 38, 3926-3930                              ", &  !22
         "Wood, Abboud, Palenik and Palenik (2000) Inorg. Chem. 39, 2065-2068             ", &  !23
         "Tytko, Mehnike and Kurad (1999) Structure and Bonding 93, 1-66                  ", &  !24
         "Gundemann, et al.(1999) J. Phys. Chem. A 103, 4752-4754                         ", &  !25
         "Zocchi (2000) Solid State Sci. 2 383-387                                        ", &  !26
         "Jensen, Palenik and Tiekiak (2001) Polyhedron 20, 2137                          ", &  !27
         "Roulhac and Palenik (2002) Inorg. Chem. 42, 118-121                             ", &  !28
         "Holsa et al.(2002) J.Solid State Chem 165, 48-55                                ", &  !29
         "Trzesowska, Kruszynski & Bartezak (2004) Acta Cryst. B60, 174-178               ", &  !30
         "Locock & Burns (2004) Z.Krist. 219, 267-271                                     ", &  !31
         "J.Rodriguez-Carvajal, Private communication                                     ", &  !32
         "S. Adams and R. Prasada Rao, (2011) Phys. Status Solidi A 208, No. 8, 1746-1753 ", &  !33
         "S. Adams (2013),  Structure and Bonding (eds. Brown & Poeppelmeier) 158, 91-128 ", &  !34
         "Adams S, Moretsky O and Canadell E (2004) Solid State Ionics 168, 281-290       "/)   !35

    !!----
    !!---- Common_OxStates_Table
    !!----    Integer, Dimension(:,:), Allocatable, public :: Common_OxStates_Table
    !!----
    !!----    Tables of Common Oxidation States
    !!----
    !!---- Created: December - 2016
    !!
    Integer, Dimension(:,:), Allocatable, public :: Common_OxStates_Table

    !!----
    !!---- OxStates_Table
    !!----    Integer, Dimension(:,:), Allocatable, public :: OxStates_Table
    !!----
    !!----    Tables of Oxidation States
    !!----
    !!---- Created: December - 2016
    !!
    Integer, Dimension(:,:), Allocatable, public :: OxStates_Table

    !!----
    !!---- PaulingX
    !!----    Real(kind=dp), Dimension(:), Allocatable, public :: PaulingX
    !!----
    !!----    Table of electronegativities
    !!----
    !!---- Created: December - 2016
    !!
    Real(kind=dp), Dimension(:), Allocatable, public :: PaulingX

    contains

    !!----
    !!---- Subroutine Deallocate_Ap_Table()
    !!----
    !!----    Deallocating Ap_Table
    !!----
    !!---- Created: January - 2015
    !!
    Subroutine Deallocate_Ap_Table()

       if (allocated(Ap_Table)) deallocate(Ap_Table)

       return
    End Subroutine Deallocate_Ap_Table

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
    !!---- Subroutine Deallocate_sBVS_Table()
    !!----
    !!----    Deallocating sBVS_Table
    !!----
    !!---- Updated: January - 2015
    !!
    Subroutine Deallocate_sBVS_Table()

       if (allocated(sBVS_Table)) deallocate(sBVS_Table)

       return
    End Subroutine Deallocate_sBVS_Table

    !!----
    !!----
    !!---- Subroutine Set_Atomic_Properties()
    !!----
    !!----  Fills the parameters for single atomic properties, including softness,
    !!----  for determining D0 and Rmin parameters of the Morse potential for non-tabulated
    !!----  values. The values are from supplementary material of the reference:
    !!----  Stefan Adams
    !!---- "Practical considerations in determining bond-valence parameters"
    !!----  Structure and Bonding (2014) 158, 91-128 / DOI 10.1007/430_2013_96
    !!----
    !!---- Created: January - 2015
    !!
    Subroutine Set_Atomic_Properties()

       if (.not. allocated(Ap_Table)) allocate(Ap_Table(Ap_species_n))

        Ap_Table(  1)=Atomic_Properties_Type( 1,"H+1 ", 1, 1,  1.0079, 1, 1, 0, 0.370, 0.00000)
        Ap_Table(  2)=Atomic_Properties_Type( 1,"HX+1", 1, 1,  1.0080, 1, 1, 0, 0.370, 0.00000)
        Ap_Table(  3)=Atomic_Properties_Type( 1,"D+1 ", 1, 1,  2.0141, 1, 1, 0, 0.370, 0.00000)
        Ap_Table(  4)=Atomic_Properties_Type( 3,"LI+1", 1, 1,  6.9410, 1, 1, 1, 1.310, 0.02847)
        Ap_Table(  5)=Atomic_Properties_Type( 4,"BE+2", 2, 2,  9.0122, 2, 2, 1, 0.910, 0.01474)
        Ap_Table(  6)=Atomic_Properties_Type( 5,"B+3 ", 3, 3, 10.8110, 1, 3, 1, 0.710, 0.00903)
        Ap_Table(  7)=Atomic_Properties_Type( 6,"C+4 ", 4, 4, 12.0107, 2, 4, 1, 0.770, 0.00611)
        Ap_Table(  8)=Atomic_Properties_Type( 6,"C+2 ", 2, 4, 12.0107, 2, 4, 1, 0.770, 0.08509)
        Ap_Table(  9)=Atomic_Properties_Type( 7,"N+5 ", 5, 5, 14.0067, 2, 5, 1, 0.740, 0.00440)
        Ap_Table( 10)=Atomic_Properties_Type( 7,"N+3 ", 3, 5, 14.0067, 2, 5, 1, 0.740, 0.06661)
        Ap_Table( 11)=Atomic_Properties_Type( 7,"N-3 ",-3, 5, 14.0067, 2, 5, 1, 1.320, 0.19700)
        Ap_Table( 12)=Atomic_Properties_Type( 8,"O-2 ",-2,-2, 15.9994, 2, 6, 1, 1.330, 0.14752)
        Ap_Table( 13)=Atomic_Properties_Type( 8,"OX-2",-2,-2, 15.9994, 2, 6, 1, 1.330, 0.14752)
        Ap_Table( 14)=Atomic_Properties_Type( 8,"O+6 ", 6,-2, 15.9994, 2, 6, 1, 0.740, 0.00333)
        Ap_Table( 15)=Atomic_Properties_Type( 9,"F-1 ",-1,-1, 18.9994, 2, 7, 1, 1.260, 0.14262)
        Ap_Table( 16)=Atomic_Properties_Type(11,"NA+1", 1, 1, 22.9898, 3, 1, 1, 1.660, 0.04746)
        Ap_Table( 17)=Atomic_Properties_Type(12,"MG+2", 2, 2, 24.3050, 3, 2, 1, 1.360, 0.03072)
        Ap_Table( 18)=Atomic_Properties_Type(13,"AL+3", 3, 3, 26.9815, 3, 3, 1, 1.110, 0.02185)
        Ap_Table( 19)=Atomic_Properties_Type(14,"SI+4", 4, 4, 28.0855, 3, 4, 1, 1.170, 0.01644)
        Ap_Table( 20)=Atomic_Properties_Type(15,"P+5 ", 5, 5, 30.9738, 3, 5, 1, 1.100, 0.01287)
        Ap_Table( 21)=Atomic_Properties_Type(15,"P+3 ", 3, 5, 30.9738, 3, 5, 1, 1.100, 0.09416)
        Ap_Table( 22)=Atomic_Properties_Type(16,"S-2 ",-2,-2, 32.0660, 3, 6, 1, 1.770, 0.22433)
        Ap_Table( 23)=Atomic_Properties_Type(16,"S+6 ", 6,-2, 32.0660, 3, 6, 1, 1.040, 0.01037)
        Ap_Table( 24)=Atomic_Properties_Type(16,"S+4 ", 4,-2, 32.0660, 3, 6, 1, 1.040, 0.07882)
        Ap_Table( 25)=Atomic_Properties_Type(17,"CL+7", 7,-1, 35.4527, 3, 7, 1, 0.990, 0.00854)
        Ap_Table( 26)=Atomic_Properties_Type(17,"CL+5", 5,-1, 35.4527, 3, 7, 1, 0.990, 0.07909)
        Ap_Table( 27)=Atomic_Properties_Type(17,"CL+3", 3,-1, 35.4527, 3, 7, 1, 0.990, 0.14437)
        Ap_Table( 28)=Atomic_Properties_Type(17,"CL-1",-1,-1, 35.4527, 3, 7, 1, 1.740, 0.21389)
        Ap_Table( 29)=Atomic_Properties_Type(18,"AR-0", 0, 0, 39.9480, 3, 8, 1, 0.000, 0.99999)
        Ap_Table( 30)=Atomic_Properties_Type(19,"K+1 ", 1, 1, 39.0983, 4, 1, 1, 2.060, 0.07328)
        Ap_Table( 31)=Atomic_Properties_Type(20,"CA+2", 2, 2, 40.0780, 4, 2, 1, 1.660, 0.05123)
        Ap_Table( 32)=Atomic_Properties_Type(21,"SC+3", 3, 3, 44.9559, 4, 3, 2, 1.460, 0.04104)
        Ap_Table( 33)=Atomic_Properties_Type(22,"TI+4", 4, 4, 47.8670, 4, 4, 2, 1.260, 0.03569)
        Ap_Table( 34)=Atomic_Properties_Type(22,"TI+3", 3, 4, 47.8670, 4, 4, 2, 1.260, 0.12678)
        Ap_Table( 35)=Atomic_Properties_Type(22,"TI+2", 2, 4, 47.8670, 4, 4, 2, 1.260, 0.14372)
        Ap_Table( 36)=Atomic_Properties_Type(23,"V+5 ", 5, 5, 50.9415, 4, 5, 2, 1.210, 0.03182)
        Ap_Table( 37)=Atomic_Properties_Type(23,"V+4 ", 4, 5, 50.9415, 4, 5, 2, 1.210, 0.10770)
        Ap_Table( 38)=Atomic_Properties_Type(23,"V+3 ", 3, 5, 50.9415, 4, 5, 2, 1.210, 0.11507)
        Ap_Table( 39)=Atomic_Properties_Type(23,"V+2 ", 2, 5, 50.9415, 4, 5, 2, 1.210, 0.13628)
        Ap_Table( 40)=Atomic_Properties_Type(24,"CR+6", 6, 6, 51.9961, 4, 6, 2, 1.210, 0.02876)
        Ap_Table( 41)=Atomic_Properties_Type(24,"CR+4", 4, 6, 51.9961, 4, 6, 2, 1.210, 0.09850)
        Ap_Table( 42)=Atomic_Properties_Type(24,"CR+3", 3, 6, 51.9961, 4, 6, 2, 1.210, 0.10989)
        Ap_Table( 43)=Atomic_Properties_Type(24,"CR+2", 2, 6, 51.9961, 4, 6, 2, 1.210, 0.13819)
        Ap_Table( 44)=Atomic_Properties_Type(24,"CR+5", 5, 6, 51.9961, 4, 6, 2, 1.210, 0.09446)
        Ap_Table( 45)=Atomic_Properties_Type(25,"MN+7", 7, 7, 54.9380, 4, 7, 2, 1.260, 0.02654)
        Ap_Table( 46)=Atomic_Properties_Type(25,"MN+6", 6, 7, 54.9380, 4, 7, 2, 1.260, 0.08464)
        Ap_Table( 47)=Atomic_Properties_Type(25,"MN+5", 5, 7, 54.9380, 4, 7, 2, 1.260, 0.08653)
        Ap_Table( 48)=Atomic_Properties_Type(25,"MN+4", 4, 7, 54.9380, 4, 7, 2, 1.260, 0.09413)
        Ap_Table( 49)=Atomic_Properties_Type(25,"MN+3", 3, 7, 54.9380, 4, 7, 2, 1.260, 0.11405)
        Ap_Table( 50)=Atomic_Properties_Type(25,"MN+2", 2, 7, 54.9380, 4, 7, 2, 1.260, 0.11097)
        Ap_Table( 51)=Atomic_Properties_Type(26,"FE+3", 3, 3, 55.8450, 4, 8, 2, 1.260, 0.08271)
        Ap_Table( 52)=Atomic_Properties_Type(26,"FE+2", 2, 3, 55.8450, 4, 8, 2, 1.260, 0.13832)
        Ap_Table( 53)=Atomic_Properties_Type(26,"FE+4", 4, 3, 55.8450, 4, 8, 2, 1.260, 0.04856)
        Ap_Table( 54)=Atomic_Properties_Type(27,"CO+3", 3, 3, 58.9332, 4, 8, 2, 1.210, 0.11232)
        Ap_Table( 55)=Atomic_Properties_Type(27,"CO+4", 4, 3, 58.9332, 4, 8, 2, 1.210, 0.07094)
        Ap_Table( 56)=Atomic_Properties_Type(27,"CO+2", 2, 3, 58.9332, 4, 8, 2, 1.210, 0.12182)
        Ap_Table( 57)=Atomic_Properties_Type(27,"CO+1", 1, 3, 58.9332, 4, 8, 2, 1.210, 0.21741)
        Ap_Table( 58)=Atomic_Properties_Type(28,"NI+2", 2, 2, 58.6934, 4, 8, 2, 1.210, 0.11752)
        Ap_Table( 59)=Atomic_Properties_Type(28,"NI+4", 4, 2, 58.6934, 4, 8, 2, 1.210, 0.09464)
        Ap_Table( 60)=Atomic_Properties_Type(28,"NI+3", 3, 2, 58.6934, 4, 8, 2, 1.210, 0.10130)
        Ap_Table( 61)=Atomic_Properties_Type(29,"CU+2", 2, 2, 63.5460, 4, 1, 2, 1.210, 0.12091)
        Ap_Table( 62)=Atomic_Properties_Type(29,"CU+3", 3, 2, 63.5460, 4, 1, 2, 1.210, 0.10890)
        Ap_Table( 63)=Atomic_Properties_Type(29,"CU+1", 1, 2, 63.5460, 4, 1, 2, 1.210, 0.15914)
        Ap_Table( 64)=Atomic_Properties_Type(30,"ZN+2", 2, 2, 65.3900, 4, 2, 2, 1.210, 0.09190)
        Ap_Table( 65)=Atomic_Properties_Type(31,"GA+3", 3, 3, 69.7230, 4, 3, 1, 1.160, 0.05998)
        Ap_Table( 66)=Atomic_Properties_Type(31,"GA+1", 1, 3, 69.7230, 4, 3, 1, 1.160, 0.13779)
        Ap_Table( 67)=Atomic_Properties_Type(32,"GE+4", 4, 4, 72.6100, 4, 4, 1, 1.220, 0.04187)
        Ap_Table( 68)=Atomic_Properties_Type(32,"GE+2", 2, 4, 72.6100, 4, 4, 1, 1.220, 0.10936)
        Ap_Table( 69)=Atomic_Properties_Type(33,"AS+5", 5, 5, 74.9216, 4, 5, 1, 1.210, 0.03079)
        Ap_Table( 70)=Atomic_Properties_Type(33,"AS+3", 3, 5, 74.9220, 4, 5, 1, 1.210, 0.09180)
        Ap_Table( 71)=Atomic_Properties_Type(34,"SE+6", 6, 6, 78.9600, 4, 6, 1, 1.170, 0.02714)
        Ap_Table( 72)=Atomic_Properties_Type(34,"SE+4", 4, 6, 78.9600, 4, 6, 1, 1.170, 0.07889)
        Ap_Table( 73)=Atomic_Properties_Type(34,"SE-2",-2, 6, 78.9600, 4, 6, 1, 1.910, 0.24167)
        Ap_Table( 74)=Atomic_Properties_Type(35,"BR+7", 7,-1, 79.9040, 4, 7, 1, 1.140, 0.02228)
        Ap_Table( 75)=Atomic_Properties_Type(35,"BR+5", 5,-1, 79.9040, 4, 7, 1, 1.140, 0.06916)
        Ap_Table( 76)=Atomic_Properties_Type(35,"BR-1",-1,-1, 79.9040, 4, 7, 1, 1.890, 0.23669)
        Ap_Table( 77)=Atomic_Properties_Type(36,"KR-0", 0, 0, 83.8000, 4, 8, 1, 0.000, 0.99999)
        Ap_Table( 78)=Atomic_Properties_Type(37,"RB+1", 1, 1, 85.4678, 5, 1, 1, 2.210, 0.08653)
        Ap_Table( 79)=Atomic_Properties_Type(38,"SR+2", 2, 2, 87.6200, 5, 2, 1, 1.860, 0.06278)
        Ap_Table( 80)=Atomic_Properties_Type(39,"Y+3 ", 3, 3, 88.9059, 5, 3, 2, 1.660, 0.04990)
        Ap_Table( 81)=Atomic_Properties_Type(40,"ZR+4", 4, 4, 91.2240, 5, 4, 2, 1.410, 0.04347)
        Ap_Table( 82)=Atomic_Properties_Type(40,"ZR+3", 3, 4, 91.2240, 5, 4, 2, 1.410, 0.17623)
        Ap_Table( 83)=Atomic_Properties_Type(40,"ZR+2", 2, 4, 91.2240, 5, 4, 2, 1.410, 0.20355)
        Ap_Table( 84)=Atomic_Properties_Type(41,"NB+5", 5, 5, 92.9064, 5, 5, 2, 1.310, 0.03883)
        Ap_Table( 85)=Atomic_Properties_Type(41,"NB+4", 4, 5, 92.9064, 5, 5, 2, 1.310, 0.16395)
        Ap_Table( 86)=Atomic_Properties_Type(41,"NB+3", 3, 5, 92.9064, 5, 5, 2, 1.310, 0.15029)
        Ap_Table( 87)=Atomic_Properties_Type(41,"NB+2", 2, 5, 92.9064, 5, 5, 2, 1.310, 0.18626)
        Ap_Table( 88)=Atomic_Properties_Type(42,"MO+6", 6, 6, 95.9400, 5, 6, 2, 1.310, 0.03519)
        Ap_Table( 89)=Atomic_Properties_Type(42,"MO+5", 5, 6, 95.9400, 5, 6, 2, 1.310, 0.13945)
        Ap_Table( 90)=Atomic_Properties_Type(42,"MO+4", 4, 6, 95.9400, 5, 6, 2, 1.310, 0.24835)
        Ap_Table( 91)=Atomic_Properties_Type(42,"MO+3", 3, 6, 95.9400, 5, 6, 2, 1.310, 0.10364)
        Ap_Table( 92)=Atomic_Properties_Type(42,"MO+2", 2, 6, 95.9400, 5, 6, 2, 1.310, 0.18239)
        Ap_Table( 93)=Atomic_Properties_Type(43,"TC+7", 7, 7, 98.0000, 5, 7, 2, 1.210, 0.02969)
        Ap_Table( 94)=Atomic_Properties_Type(43,"TC+4", 4, 7, 98.0000, 5, 7, 2, 1.210, 0.12061)
        Ap_Table( 95)=Atomic_Properties_Type(44,"RU+6", 6, 6,101.0700, 5, 8, 2, 1.160, 0.10721)
        Ap_Table( 96)=Atomic_Properties_Type(44,"RU+5", 5, 6,101.0700, 5, 8, 2, 1.160, 0.11351)
        Ap_Table( 97)=Atomic_Properties_Type(44,"RU+4", 4, 6,101.0700, 5, 8, 2, 1.160, 0.12061)
        Ap_Table( 98)=Atomic_Properties_Type(44,"RU+3", 3, 6,101.0700, 5, 8, 2, 1.160, 0.11008)
        Ap_Table( 99)=Atomic_Properties_Type(45,"RH+4", 4, 4,102.9060, 5, 8, 2, 1.210, 0.09189)
        Ap_Table(100)=Atomic_Properties_Type(45,"RH+3", 3, 4,102.9060, 5, 8, 2, 1.210, 0.13754)
        Ap_Table(101)=Atomic_Properties_Type(46,"PD+2", 2, 2,106.4200, 5, 8, 2, 1.260, 0.14821)
        Ap_Table(102)=Atomic_Properties_Type(46,"PD+4", 4, 2,106.4200, 5, 8, 2, 1.260, 0.12061)
        Ap_Table(103)=Atomic_Properties_Type(47,"AG+1", 1, 1,107.8680, 5, 1, 2, 1.460, 0.14379)
        Ap_Table(104)=Atomic_Properties_Type(47,"AG+3", 3, 1,107.8680, 5, 1, 2, 1.460, 0.11774)
        Ap_Table(105)=Atomic_Properties_Type(47,"AG+2", 2, 1,107.8680, 5, 1, 2, 1.460, 0.14982)
        Ap_Table(106)=Atomic_Properties_Type(48,"CD+2", 2, 2,112.4110, 5, 2, 2, 1.410, 0.09723)
        Ap_Table(107)=Atomic_Properties_Type(49,"IN+3", 3, 3,114.8180, 5, 3, 1, 1.410, 0.07700)
        Ap_Table(108)=Atomic_Properties_Type(49,"IN+1", 1, 3,114.8180, 5, 3, 1, 1.410, 0.15286)
        Ap_Table(109)=Atomic_Properties_Type(50,"SN+4", 4, 4,118.7100, 5, 4, 1, 1.400, 0.05473)
        Ap_Table(110)=Atomic_Properties_Type(50,"SN+2", 2, 4,118.7000, 5, 4, 1, 1.400, 0.12603)
        Ap_Table(111)=Atomic_Properties_Type(51,"SB+5", 5, 5,121.7600, 5, 5, 1, 1.410, 0.03859)
        Ap_Table(112)=Atomic_Properties_Type(51,"SB+3", 3, 5,121.7600, 5, 5, 1, 1.410, 0.10603)
        Ap_Table(113)=Atomic_Properties_Type(52,"TE+6", 6, 6,127.6000, 5, 6, 1, 1.370, 0.03025)
        Ap_Table(114)=Atomic_Properties_Type(52,"TE+4", 4, 6,127.6000, 5, 6, 1, 1.370, 0.09377)
        Ap_Table(115)=Atomic_Properties_Type(52,"TE+2", 2, 6,127.6000, 5, 6, 1, 1.370, 0.21252)
        Ap_Table(116)=Atomic_Properties_Type(52,"TE-2",-2, 6,127.6000, 5, 6, 1, 2.140, 0.26716)
        Ap_Table(117)=Atomic_Properties_Type(53,"I+7 ", 7,-1,126.9040, 5, 7, 1, 1.330, 0.02506)
        Ap_Table(118)=Atomic_Properties_Type(53,"I+5 ", 5,-1,126.9040, 5, 7, 1, 1.330, 0.08040)
        Ap_Table(119)=Atomic_Properties_Type(53,"I-1 ",-1,-1,126.9040, 5, 7, 1, 2.130, 0.27057)
        Ap_Table(120)=Atomic_Properties_Type(54,"XE+6", 6, 0,131.2900, 5, 8, 1, 1.300, 0.07147)
        Ap_Table(121)=Atomic_Properties_Type(54,"XE+4", 4, 0,131.2900, 5, 8, 1, 1.300, 0.16081)
        Ap_Table(122)=Atomic_Properties_Type(54,"XE+2", 2, 0,131.2900, 5, 8, 1, 1.300, 0.18361)
        Ap_Table(123)=Atomic_Properties_Type(55,"CS+1", 1, 1,132.9050, 6, 1, 1, 2.460, 0.10383)
        Ap_Table(124)=Atomic_Properties_Type(56,"BA+2", 2, 2,137.3270, 6, 2, 1, 2.010, 0.07324)
        Ap_Table(125)=Atomic_Properties_Type(57,"LA+3", 3, 3,138.9060, 6, 0, 3, 1.810, 0.06500)
        Ap_Table(126)=Atomic_Properties_Type(58,"CE+4", 4, 4,140.1160, 6, 0, 3, 1.710, 0.06946)
        Ap_Table(127)=Atomic_Properties_Type(58,"CE+3", 3, 4,140.1160, 6, 0, 3, 1.710, 0.12076)
        Ap_Table(128)=Atomic_Properties_Type(59,"PR+3", 3, 3,140.9080, 6, 0, 3, 1.710, 0.11521)
        Ap_Table(129)=Atomic_Properties_Type(60,"ND+3", 3, 3,144.2400, 6, 0, 3, 1.710, 0.10902)
        Ap_Table(130)=Atomic_Properties_Type(62,"SM+3", 3, 3,150.3600, 6, 0, 3, 1.710, 0.11154)
        Ap_Table(131)=Atomic_Properties_Type(63,"EU+3", 3, 3,151.9640, 6, 0, 3, 1.710, 0.11245)
        Ap_Table(132)=Atomic_Properties_Type(63,"EU+2", 2, 3,151.9640, 6, 0, 3, 1.710, 0.14630)
        Ap_Table(133)=Atomic_Properties_Type(64,"GD+3", 3, 3,157.2500, 6, 0, 3, 1.660, 0.08538)
        Ap_Table(134)=Atomic_Properties_Type(65,"TB+3", 3, 3,158.9254, 6, 0, 3, 1.610, 0.11187)
        Ap_Table(135)=Atomic_Properties_Type(65,"TB+4", 4, 3,158.9254, 6, 0, 3, 1.610, 0.99999)
        Ap_Table(136)=Atomic_Properties_Type(66,"DY+3", 3, 3,162.5000, 6, 0, 3, 1.610, 0.10780)
        Ap_Table(137)=Atomic_Properties_Type(67,"HO+3", 3, 3,164.9303, 6, 0, 3, 1.610, 0.10178)
        Ap_Table(138)=Atomic_Properties_Type(68,"ER+3", 3, 3,167.2590, 6, 0, 3, 1.610, 0.10019)
        Ap_Table(139)=Atomic_Properties_Type(69,"TM+3", 3, 3,168.9342, 6, 0, 3, 1.610, 0.10516)
        Ap_Table(140)=Atomic_Properties_Type(70,"YB+3", 3, 3,173.0545, 6, 0, 3, 1.610, 0.10805)
        Ap_Table(141)=Atomic_Properties_Type(70,"YB+2", 2, 3,173.0545, 6, 0, 3, 1.610, 0.15535)
        Ap_Table(142)=Atomic_Properties_Type(71,"LU+3", 3, 3,174.9670, 6, 3, 2, 1.610, 0.08220)
        Ap_Table(143)=Atomic_Properties_Type(72,"HF+4", 4, 4,178.4900, 6, 4, 2, 1.410, 0.05000)
        Ap_Table(144)=Atomic_Properties_Type(73,"TA+5", 5, 5,180.9480, 6, 5, 2, 1.310, 0.04545)
        Ap_Table(145)=Atomic_Properties_Type(73,"TA+4", 4, 5,112.4110, 6, 5, 2, 1.310, 0.17543)
        Ap_Table(146)=Atomic_Properties_Type(73,"TA+3", 3, 5,180.9490, 6, 5, 2, 1.310, 0.17543)
        Ap_Table(147)=Atomic_Properties_Type(73,"TA+2", 2, 5,180.9490, 6, 5, 2, 1.310, 0.32162)
        Ap_Table(148)=Atomic_Properties_Type(74,"W+6 ", 6, 6,183.8400, 6, 6, 2, 1.210, 0.03846)
        Ap_Table(149)=Atomic_Properties_Type(74,"W+5 ", 5, 6,183.8400, 6, 6, 2, 1.210, 0.14843)
        Ap_Table(150)=Atomic_Properties_Type(74,"W+4 ", 4, 6,183.8400, 6, 6, 2, 1.210, 0.16080)
        Ap_Table(151)=Atomic_Properties_Type(74,"W+3 ", 3, 6,183.8400, 6, 6, 2, 1.210, 0.17544)
        Ap_Table(152)=Atomic_Properties_Type(74,"W+2 ", 2, 6,183.8400, 6, 6, 2, 1.210, 0.32154)
        Ap_Table(153)=Atomic_Properties_Type(75,"RE+7", 7, 7,186.2070, 6, 7, 2, 1.210, 0.03333)
        Ap_Table(154)=Atomic_Properties_Type(75,"RE+6", 6, 7,186.2070, 6, 7, 2, 1.210, 0.14844)
        Ap_Table(155)=Atomic_Properties_Type(75,"RE+5", 5, 7,186.2070, 6, 7, 2, 1.210, 0.13784)
        Ap_Table(156)=Atomic_Properties_Type(75,"RE+4", 4, 7,186.2070, 6, 7, 2, 1.210, 0.15315)
        Ap_Table(157)=Atomic_Properties_Type(75,"RE+3", 3, 7,186.2070, 6, 7, 2, 1.210, 0.17077)
        Ap_Table(158)=Atomic_Properties_Type(76,"OS+8", 8, 8,190.2300, 6, 8, 2, 1.160, 0.03030)
        Ap_Table(159)=Atomic_Properties_Type(76,"OS+7", 7, 8,190.2300, 6, 8, 2, 1.160, 0.13784)
        Ap_Table(160)=Atomic_Properties_Type(76,"OS+6", 6, 8,190.2300, 6, 8, 2, 1.160, 0.12865)
        Ap_Table(161)=Atomic_Properties_Type(76,"OS+5", 5, 8,190.2300, 6, 8, 2, 1.160, 0.13784)
        Ap_Table(162)=Atomic_Properties_Type(76,"OS+4", 4, 8,190.2300, 6, 8, 2, 1.160, 0.14844)
        Ap_Table(163)=Atomic_Properties_Type(77,"IR+5", 5, 5,192.2170, 6, 8, 2, 1.210, 0.13784)
        Ap_Table(164)=Atomic_Properties_Type(77,"IR+4", 4, 5,192.2170, 6, 8, 2, 1.210, 0.11351)
        Ap_Table(165)=Atomic_Properties_Type(77,"IR+3", 3, 5,192.2170, 6, 8, 2, 1.210, 0.16081)
        Ap_Table(166)=Atomic_Properties_Type(78,"PT+2", 2, 2,195.0780, 6, 8, 2, 1.210, 0.19125)
        Ap_Table(167)=Atomic_Properties_Type(78,"PT+5", 5, 2,195.0780, 6, 8, 2, 1.210, 0.10156)
        Ap_Table(168)=Atomic_Properties_Type(78,"PT+4", 4, 2,195.0780, 6, 8, 2, 1.210, 0.13784)
        Ap_Table(169)=Atomic_Properties_Type(78,"PT+3", 3, 2,195.0780, 6, 8, 2, 1.210, 0.17543)
        Ap_Table(170)=Atomic_Properties_Type(79,"AU+3", 3, 3,196.9670, 6, 1, 2, 1.210, 0.14844)
        Ap_Table(171)=Atomic_Properties_Type(79,"AU+1", 1, 3,196.9670, 6, 1, 2, 1.210, 0.17705)
        Ap_Table(172)=Atomic_Properties_Type(80,"HG+2", 2, 2,200.5900, 6, 2, 2, 1.360, 0.12951)
        Ap_Table(173)=Atomic_Properties_Type(80,"HG+1", 1, 2,200.5900, 6, 2, 2, 1.360, 0.24034)
        Ap_Table(174)=Atomic_Properties_Type(81,"TL+3", 3, 3,204.3830, 6, 3, 1, 1.760, 0.09582)
        Ap_Table(175)=Atomic_Properties_Type(81,"TL+1", 1, 3,204.3830, 6, 3, 1, 1.760, 0.13967)
        Ap_Table(176)=Atomic_Properties_Type(82,"PB+4", 4, 4,207.2000, 6, 4, 1, 1.660, 0.07547)
        Ap_Table(177)=Atomic_Properties_Type(82,"PB+2", 2, 4,207.2000, 6, 4, 1, 1.660, 0.11831)
        Ap_Table(178)=Atomic_Properties_Type(83,"BI+5", 5, 5,208.9800, 6, 5, 1, 1.460, 0.06185)
        Ap_Table(179)=Atomic_Properties_Type(83,"BI+3", 3, 5,208.9800, 6, 5, 1, 1.460, 0.10135)
        Ap_Table(180)=Atomic_Properties_Type(90,"TH+4", 4, 4,232.0400, 7, 0, 3, 1.660, 0.06900)
        Ap_Table(181)=Atomic_Properties_Type(92,"U+6 ", 6, 6,238.0300, 7, 0, 3, 1.610, 0.03000)
        Ap_Table(182)=Atomic_Properties_Type(-1,"NH+1", 1, 1, 18.0400, 4, 1, 0, 2.100, 0.06993)
        Ap_Table(183)=Atomic_Properties_Type(-1,"DU-1",-1,-1,  1.0079, 1, 1, 0, 0.000, 0.00000)
    End Subroutine Set_Atomic_Properties

    !!----
    !!----
    !!---- Subroutine Set_BVEL_Table()
    !!----
    !!----  Fills the parameters for BVEL from
    !!----  Stefan Adams and R. Prasada Rao
    !!---- "High power lithium ion battery materials by computational design"
    !!----  Phys. Status Solidi A 208, No. 8, 1746, 1753 (2011) / DOI 10.1002/pssa.201001116
    !!----
    !!---- Only one anion is available (O-2) for the moment
    !!----
    !!---- Created: December - 2014
    !!
    Subroutine Set_BVEL_Table()

       if (.not. allocated(BVEL_Table)) allocate(BVEL_Table(bvel_species_n))

        BVEL_Table(  1)=BVEL_Par_Type("H+1",  (/  1.923000/),(/  0.870450/),(/  4.000000/),(/  1.885800/),(/  1.127680/),(/  2.188184/),(/33/))
        BVEL_Table(  2)=BVEL_Par_Type("LI+1", (/  5.021000/),(/  1.170960/),(/  5.500000/),(/  0.988160/),(/  1.940010/),(/  1.937984/),(/33/))
        BVEL_Table(  3)=BVEL_Par_Type("BE+2", (/  4.000000/),(/  1.209030/),(/  5.500000/),(/  2.768820/),(/  1.522170/),(/  1.848429/),(/33/))
        BVEL_Table(  4)=BVEL_Par_Type("B+3",  (/  3.417000/),(/  1.357610/),(/  4.500000/),(/  2.389240/),(/  1.340030/),(/  2.597403/),(/33/))
        BVEL_Table(  5)=BVEL_Par_Type("C+4",  (/  3.000000/),(/  1.398260/),(/  5.000000/),(/  4.791870/),(/  1.200890/),(/  2.237136/),(/33/))
        BVEL_Table(  6)=BVEL_Par_Type("C+2",  (/  1.000000/),(/  1.413680/),(/  5.000000/),(/  2.405530/),(/  1.030980/),(/  2.409639/),(/33/))
        BVEL_Table(  7)=BVEL_Par_Type("N+5",  (/  3.000000/),(/  1.462670/),(/  5.000000/),(/  6.276770/),(/  1.161420/),(/  2.222222/),(/33/))
        BVEL_Table(  8)=BVEL_Par_Type("N+3",  (/  2.000000/),(/  1.407950/),(/  5.000000/),(/  3.810890/),(/  1.137580/),(/  2.232143/),(/33/))
        BVEL_Table(  9)=BVEL_Par_Type("NH4+1",(/  3.467000/),(/  2.033800/),(/  6.000000/),(/  0.405370/),(/  2.453640/),(/  2.262443/),(/33/))
        BVEL_Table( 10)=BVEL_Par_Type("NA+1", (/  6.520000/),(/  1.562250/),(/  6.000000/),(/  0.575230/),(/  2.374330/),(/  2.074689/),(/33/))
        BVEL_Table( 11)=BVEL_Par_Type("MG+2", (/  5.897000/),(/  1.483980/),(/  5.500000/),(/  1.575540/),(/  1.956270/),(/  1.953125/),(/33/))
        BVEL_Table( 12)=BVEL_Par_Type("AL+3", (/  5.327000/),(/  1.599010/),(/  5.000000/),(/  1.803460/),(/  1.758060/),(/  2.358491/),(/33/))
        BVEL_Table( 13)=BVEL_Par_Type("SI+4", (/  4.100000/),(/  1.608170/),(/  5.000000/),(/  2.857200/),(/  1.535940/),(/  2.314815/),(/33/))
        BVEL_Table( 14)=BVEL_Par_Type("P+5",  (/  4.000000/),(/  1.620380/),(/  5.000000/),(/  3.896350/),(/  1.440660/),(/  2.288330/),(/33/))
        BVEL_Table( 15)=BVEL_Par_Type("P+3",  (/  3.000000/),(/  1.515550/),(/  4.500000/),(/  2.020620/),(/  1.410510/),(/  2.487562/),(/33/))
        BVEL_Table( 16)=BVEL_Par_Type("S+6",  (/  4.000000/),(/  1.642200/),(/  5.000000/),(/  4.967260/),(/  1.381020/),(/  2.267574/),(/33/))
        BVEL_Table( 17)=BVEL_Par_Type("S+4",  (/  3.000000/),(/  1.642820/),(/  5.500000/),(/  3.036720/),(/  1.411880/),(/  2.341920/),(/33/))
        BVEL_Table( 18)=BVEL_Par_Type("CL+7", (/  4.000000/),(/  1.679460/),(/  5.000000/),(/  5.991000/),(/  1.348010/),(/  2.257336/),(/33/))
        BVEL_Table( 19)=BVEL_Par_Type("CL+5", (/  3.000000/),(/  1.695520/),(/  5.500000/),(/  4.290890/),(/  1.356530/),(/  2.247191/),(/33/))
        BVEL_Table( 20)=BVEL_Par_Type("CL+3", (/  2.000000/),(/  1.722650/),(/  5.500000/),(/  3.071190/),(/  1.384410/),(/  2.036660/),(/33/))
        BVEL_Table( 21)=BVEL_Par_Type("K+1",  (/  8.846000/),(/  1.941170/),(/  6.000000/),(/  0.349850/),(/  2.766360/),(/  2.293578/),(/33/))
        BVEL_Table( 22)=BVEL_Par_Type("CA+2", (/  7.544000/),(/  1.795190/),(/  5.500000/),(/  0.994290/),(/  2.320320/),(/  2.100840/),(/33/))
        BVEL_Table( 23)=BVEL_Par_Type("SC+3", (/  6.255000/),(/  1.732200/),(/  5.500000/),(/  2.156100/),(/  1.996150/),(/  2.024291/),(/33/))
        BVEL_Table( 24)=BVEL_Par_Type("TI+4", (/  6.000000/),(/  1.723940/),(/  5.500000/),(/  2.813330/),(/  1.831440/),(/  1.988072/),(/33/))
        BVEL_Table( 25)=BVEL_Par_Type("TI+3", (/  6.000000/),(/  1.697660/),(/  5.500000/),(/  1.978510/),(/  1.886190/),(/  2.173913/),(/33/))
        BVEL_Table( 26)=BVEL_Par_Type("V+5",  (/  4.166000/),(/  1.794450/),(/  5.500000/),(/  3.695330/),(/  1.602580/),(/  1.960784/),(/33/))
        BVEL_Table( 27)=BVEL_Par_Type("V+4",  (/  5.738000/),(/  1.749320/),(/  5.000000/),(/  2.080470/),(/  1.776380/),(/  2.347418/),(/33/))
        BVEL_Table( 28)=BVEL_Par_Type("V+3",  (/  6.000000/),(/  1.677990/),(/  5.500000/),(/  1.829360/),(/  1.857970/),(/  2.277904/),(/33/))
        BVEL_Table( 29)=BVEL_Par_Type("CR+6", (/  4.000000/),(/  1.824710/),(/  5.500000/),(/  3.687510/),(/  1.532510/),(/  2.100840/),(/33/))
        BVEL_Table( 30)=BVEL_Par_Type("CR+5", (/  4.000000/),(/  1.767810/),(/  5.500000/),(/  2.365510/),(/  1.555460/),(/  2.487562/),(/33/))
        BVEL_Table( 31)=BVEL_Par_Type("CR+4", (/  5.429000/),(/  1.760950/),(/  5.500000/),(/  1.933290/),(/  1.762090/),(/  2.444988/),(/33/))
        BVEL_Table( 32)=BVEL_Par_Type("CR+3", (/  6.000000/),(/  1.661980/),(/  5.500000/),(/  1.773350/),(/  1.838870/),(/  2.325581/),(/33/))
        BVEL_Table( 33)=BVEL_Par_Type("MN+7", (/  4.000000/),(/  1.873620/),(/  6.500000/),(/  4.916300/),(/  1.481710/),(/  1.923077/),(/33/))
        BVEL_Table( 34)=BVEL_Par_Type("MN+6", (/  4.000000/),(/  1.820180/),(/  5.500000/),(/  2.822360/),(/  1.529310/),(/  2.403846/),(/33/))
        BVEL_Table( 35)=BVEL_Par_Type("MN+5", (/  4.000000/),(/  1.788790/),(/  5.500000/),(/  2.464560/),(/  1.575770/),(/  2.421308/),(/33/))
        BVEL_Table( 36)=BVEL_Par_Type("MN+4", (/  5.923000/),(/  1.732720/),(/  5.000000/),(/  1.858860/),(/  1.770450/),(/  2.487562/),(/33/))
        BVEL_Table( 37)=BVEL_Par_Type("MN+3", (/  5.862000/),(/  1.689930/),(/  5.500000/),(/  1.812830/),(/  1.857860/),(/  2.288330/),(/33/))
        BVEL_Table( 38)=BVEL_Par_Type("MN+2", (/  5.910000/),(/  1.627580/),(/  5.500000/),(/  1.641430/),(/  2.029690/),(/  2.079002/),(/33/))
        BVEL_Table( 39)=BVEL_Par_Type("FE+4", (/  6.000000/),(/  1.765590/),(/  5.500000/),(/  1.872850/),(/  1.827860/),(/  2.439024/),(/33/))
        BVEL_Table( 40)=BVEL_Par_Type("FE+3", (/  5.733000/),(/  1.708400/),(/  5.000000/),(/  1.666810/),(/  1.866470/),(/  2.380952/),(/33/))
        BVEL_Table( 41)=BVEL_Par_Type("FE+2", (/  5.743000/),(/  1.579110/),(/  5.500000/),(/  1.692690/),(/  1.960050/),(/  2.083333/),(/33/))
        BVEL_Table( 42)=BVEL_Par_Type("NI+3", (/  6.000000/),(/  1.648880/),(/  5.500000/),(/  1.661910/),(/  1.818870/),(/  2.415459/),(/33/))
        BVEL_Table( 43)=BVEL_Par_Type("NI+2", (/  5.933000/),(/  1.559200/),(/  5.500000/),(/  1.468410/),(/  1.924520/),(/  2.257336/),(/33/))
        BVEL_Table( 44)=BVEL_Par_Type("CO+3", (/  6.000000/),(/  1.592340/),(/  5.500000/),(/  1.870240/),(/  1.776200/),(/  2.304147/),(/33/))
        BVEL_Table( 45)=BVEL_Par_Type("CO+2", (/  5.506000/),(/  1.597730/),(/  5.500000/),(/  1.514760/),(/  1.933620/),(/  2.217295/),(/33/))
        BVEL_Table( 46)=BVEL_Par_Type("CU+3", (/  4.000000/),(/  1.709640/),(/  5.000000/),(/  1.882420/),(/  1.708230/),(/  2.341920/),(/33/))
        BVEL_Table( 47)=BVEL_Par_Type("CU+2", (/  2.560000/),(/  1.574220/),(/  5.000000/),(/  1.853410/),(/  1.566330/),(/  2.227171/),(/33/))
        BVEL_Table( 48)=BVEL_Par_Type("CU+1", (/  2.560000/),(/  1.587300/),(/  5.000000/),(/  0.664170/),(/  1.782690/),(/  2.932551/),(/33/))
        BVEL_Table( 49)=BVEL_Par_Type("ZN+2", (/  4.718000/),(/  1.653440/),(/  5.000000/),(/  1.240310/),(/  1.885570/),(/  2.481390/),(/33/))
        BVEL_Table( 50)=BVEL_Par_Type("GA+3", (/  4.905000/),(/  1.716060/),(/  5.000000/),(/  1.184560/),(/  1.793910/),(/  2.680965/),(/33/))
        BVEL_Table( 51)=BVEL_Par_Type("GE+4", (/  4.305000/),(/  1.739390/),(/  5.000000/),(/  1.913750/),(/  1.668720/),(/  2.525253/),(/33/))
        BVEL_Table( 52)=BVEL_Par_Type("AS+5", (/  4.029000/),(/  1.766890/),(/  5.000000/),(/  2.719340/),(/  1.581270/),(/  2.433090/),(/33/))
        BVEL_Table( 53)=BVEL_Par_Type("AS+3", (/  3.000000/),(/  1.767060/),(/  5.000000/),(/  1.514930/),(/  1.645540/),(/  2.475248/),(/33/))
        BVEL_Table( 54)=BVEL_Par_Type("SE+6", (/  4.000000/),(/  1.798660/),(/  5.500000/),(/  3.448650/),(/  1.532870/),(/  2.403846/),(/33/))
        BVEL_Table( 55)=BVEL_Par_Type("SE+4", (/  3.000000/),(/  1.800950/),(/  5.500000/),(/  2.380820/),(/  1.559570/),(/  2.341920/),(/33/))
        BVEL_Table( 56)=BVEL_Par_Type("BR+7", (/  4.000000/),(/  1.836580/),(/  5.500000/),(/  4.243390/),(/  1.502740/),(/  2.364066/),(/33/))
        BVEL_Table( 57)=BVEL_Par_Type("RB+1", (/ 10.020000/),(/  2.085970/),(/  6.500000/),(/  0.268130/),(/  2.896830/),(/  2.421308/),(/33/))
        BVEL_Table( 58)=BVEL_Par_Type("SR+2", (/  9.400000/),(/  1.953110/),(/  5.500000/),(/  0.743510/),(/  2.535890/),(/  2.197802/),(/33/))
        BVEL_Table( 59)=BVEL_Par_Type("Y+3",  (/  7.285000/),(/  1.903840/),(/  5.500000/),(/  1.627010/),(/  2.215230/),(/  2.092050/),(/33/))
        BVEL_Table( 60)=BVEL_Par_Type("ZR+4", (/  6.765000/),(/  1.845050/),(/  5.500000/),(/  2.191030/),(/  1.996020/),(/  2.040816/),(/33/))
        BVEL_Table( 61)=BVEL_Par_Type("NB+5", (/  6.044000/),(/  1.865880/),(/  5.500000/),(/  2.723260/),(/  1.854590/),(/  2.008032/),(/33/))
        BVEL_Table( 62)=BVEL_Par_Type("NB+4", (/  6.000000/),(/  1.785430/),(/  6.000000/),(/  2.709600/),(/  1.859890/),(/  1.901141/),(/33/))
        BVEL_Table( 63)=BVEL_Par_Type("NB+3", (/  6.000000/),(/  1.745810/),(/  6.000000/),(/  2.028480/),(/  1.951900/),(/  1.996008/),(/33/))
        BVEL_Table( 64)=BVEL_Par_Type("MO+6", (/  4.764000/),(/  1.909340/),(/  5.000000/),(/  1.991500/),(/  1.712540/),(/  2.557545/),(/33/))
        BVEL_Table( 65)=BVEL_Par_Type("W+5",  (/  6.000000/),(/  1.819750/),(/  6.000000/),(/  2.615700/),(/  1.762610/),(/  2.008032/),(/33/))
        BVEL_Table( 66)=BVEL_Par_Type("W+4",  (/  6.000000/),(/  1.745580/),(/  6.000000/),(/  2.471140/),(/  1.819450/),(/  1.923077/),(/33/))
        BVEL_Table( 67)=BVEL_Par_Type("RE+7", (/  4.098000/),(/  1.977920/),(/  6.000000/),(/  3.555930/),(/  1.596340/),(/  1.968504/),(/33/))
        BVEL_Table( 68)=BVEL_Par_Type("RE+6", (/  5.500000/),(/  1.910070/),(/  6.000000/),(/  2.950990/),(/  1.711470/),(/  2.008032/),(/33/))
        BVEL_Table( 69)=BVEL_Par_Type("RE+5", (/  6.000000/),(/  1.826640/),(/  6.000000/),(/  2.410990/),(/  1.769140/),(/  2.087683/),(/33/))
        BVEL_Table( 70)=BVEL_Par_Type("RE+3", (/  6.000000/),(/  2.207100/),(/  6.000000/),(/  0.810670/),(/  2.332180/),(/  2.493766/),(/33/))
        BVEL_Table( 71)=BVEL_Par_Type("OS+8", (/  5.333000/),(/  1.977280/),(/  6.000000/),(/  3.710190/),(/  1.661460/),(/  1.953125/),(/33/))
        BVEL_Table( 72)=BVEL_Par_Type("OS+7", (/  6.000000/),(/  1.957750/),(/  5.500000/),(/  2.919480/),(/  1.728690/),(/  2.087683/),(/33/))
        BVEL_Table( 73)=BVEL_Par_Type("OS+6", (/  6.000000/),(/  1.931920/),(/  5.500000/),(/  2.448710/),(/  1.782800/),(/  2.159827/),(/33/))
        BVEL_Table( 74)=BVEL_Par_Type("OS+4", (/  6.000000/),(/  1.753020/),(/  6.000000/),(/  2.275240/),(/  1.812440/),(/  2.008032/),(/33/))
        BVEL_Table( 75)=BVEL_Par_Type("IR+5", (/  6.000000/),(/  1.897910/),(/  6.000000/),(/  2.324760/),(/  1.834760/),(/  2.087683/),(/33/))
        BVEL_Table( 76)=BVEL_Par_Type("IR+4", (/  6.000000/),(/  1.832330/),(/  5.500000/),(/  1.686670/),(/  1.874020/),(/  2.293578/),(/33/))
        BVEL_Table( 77)=BVEL_Par_Type("PT+4", (/  6.000000/),(/  1.821980/),(/  5.500000/),(/  2.038250/),(/  1.871740/),(/  2.087683/),(/33/))
        BVEL_Table( 78)=BVEL_Par_Type("PT+2", (/  4.000000/),(/  1.512050/),(/  5.500000/),(/  2.149990/),(/  1.801790/),(/  1.742160/),(/33/))
        BVEL_Table( 79)=BVEL_Par_Type("AU+3", (/  4.000000/),(/  1.817610/),(/  5.500000/),(/  1.969670/),(/  1.813120/),(/  2.008032/),(/33/))
        BVEL_Table( 80)=BVEL_Par_Type("AU+1", (/  2.000000/),(/  1.718190/),(/  5.500000/),(/  0.853040/),(/  1.895430/),(/  2.267574/),(/33/))
        BVEL_Table( 81)=BVEL_Par_Type("HG+2", (/  6.966000/),(/  1.812760/),(/  5.500000/),(/  1.128520/),(/  2.252750/),(/  2.150538/),(/33/))
        BVEL_Table( 82)=BVEL_Par_Type("HG+1", (/  4.786000/),(/  1.812800/),(/  5.500000/),(/  0.739310/),(/  2.431550/),(/  2.150538/),(/33/))
        BVEL_Table( 83)=BVEL_Par_Type("TL+3", (/  5.220000/),(/  2.062970/),(/  5.000000/),(/  0.676370/),(/  2.106420/),(/  2.958580/),(/33/))
        BVEL_Table( 84)=BVEL_Par_Type("TL+1", (/  8.030000/),(/  1.917520/),(/  6.000000/),(/  0.349990/),(/  2.770860/),(/  2.070393/),(/33/))
        BVEL_Table( 85)=BVEL_Par_Type("PB+4", (/  5.740000/),(/  2.032930/),(/  5.000000/),(/  1.027190/),(/  2.028570/),(/  2.824859/),(/33/))
        BVEL_Table( 86)=BVEL_Par_Type("PB+2", (/  7.541000/),(/  2.018250/),(/  5.500000/),(/  0.638330/),(/  2.441910/),(/  2.309469/),(/33/))
        BVEL_Table( 87)=BVEL_Par_Type("BI+5", (/  6.000000/),(/  2.044980/),(/  5.000000/),(/  1.440500/),(/  1.985990/),(/  2.695418/),(/33/))
        BVEL_Table( 88)=BVEL_Par_Type("BI+3", (/  6.058000/),(/  2.036770/),(/  5.500000/),(/  0.979040/),(/  2.183210/),(/  2.415459/),(/33/))
        BVEL_Table( 89)=BVEL_Par_Type("MO+5", (/  5.980000/),(/  1.847600/),(/  5.500000/),(/  2.648020/),(/  1.786700/),(/  2.074689/),(/33/))
        BVEL_Table( 90)=BVEL_Par_Type("MO+4", (/  6.000000/),(/  1.723900/),(/  6.500000/),(/  3.108070/),(/  1.850990/),(/  1.779359/),(/33/))
        BVEL_Table( 91)=BVEL_Par_Type("MO+3", (/  5.700000/),(/  1.789330/),(/  5.500000/),(/  1.428260/),(/  1.929740/),(/  2.392344/),(/33/))
        BVEL_Table( 92)=BVEL_Par_Type("RU+6", (/  4.500000/),(/  1.925790/),(/  5.500000/),(/  2.421090/),(/  1.664310/),(/  2.352941/),(/33/))
        BVEL_Table( 93)=BVEL_Par_Type("RU+5", (/  6.000000/),(/  1.874420/),(/  5.500000/),(/  2.132080/),(/  1.815710/),(/  2.293578/),(/33/))
        BVEL_Table( 94)=BVEL_Par_Type("RU+4", (/  6.000000/),(/  1.793630/),(/  5.500000/),(/  1.995130/),(/  1.840530/),(/  2.227171/),(/33/))
        BVEL_Table( 95)=BVEL_Par_Type("RH+4", (/  6.000000/),(/  1.776750/),(/  5.500000/),(/  1.627250/),(/  1.817930/),(/  2.481390/),(/33/))
        BVEL_Table( 96)=BVEL_Par_Type("RH+3", (/  6.000000/),(/  1.670130/),(/  5.500000/),(/  1.928260/),(/  1.869150/),(/  2.092050/),(/33/))
        BVEL_Table( 97)=BVEL_Par_Type("PD+4", (/  5.333000/),(/  1.805000/),(/  5.500000/),(/  2.042180/),(/  1.798130/),(/  2.227171/),(/33/))
        BVEL_Table( 98)=BVEL_Par_Type("PD+2", (/  4.000000/),(/  1.623590/),(/  5.500000/),(/  1.739100/),(/  1.836710/),(/  2.008032/),(/33/))
        BVEL_Table( 99)=BVEL_Par_Type("AG+1", (/  4.438000/),(/  1.782390/),(/  5.000000/),(/  0.635190/),(/  2.225780/),(/  2.538071/),(/33/))
        BVEL_Table(100)=BVEL_Par_Type("CD+2", (/  6.176000/),(/  1.839260/),(/  5.500000/),(/  0.983460/),(/  2.169400/),(/  2.457002/),(/33/))
        BVEL_Table(101)=BVEL_Par_Type("TA+4", (/  5.500000/),(/  1.756320/),(/  6.000000/),(/  2.756550/),(/  1.798260/),(/  1.831502/),(/33/))
        BVEL_Table(102)=BVEL_Par_Type("IN+3", (/  6.024000/),(/  1.903050/),(/  5.000000/),(/  0.840760/),(/  2.024710/),(/  2.832861/),(/33/))
        BVEL_Table(103)=BVEL_Par_Type("SN+4", (/  6.069000/),(/  1.890190/),(/  5.000000/),(/  1.352680/),(/  1.934220/),(/  2.638522/),(/33/))
        BVEL_Table(104)=BVEL_Par_Type("SN+2", (/  3.325000/),(/  1.874990/),(/  5.500000/),(/  0.972610/),(/  1.964200/),(/  2.183406/),(/33/))
        BVEL_Table(105)=BVEL_Par_Type("SB+5", (/  6.000000/),(/  1.897680/),(/  5.500000/),(/  1.955230/),(/  1.863180/),(/  2.500000/),(/33/))
        BVEL_Table(106)=BVEL_Par_Type("SB+3", (/  6.000000/),(/  1.920360/),(/  5.000000/),(/  1.177860/),(/  2.075260/),(/  2.364066/),(/33/))
        BVEL_Table(107)=BVEL_Par_Type("I+7",  (/  5.800000/),(/  1.922740/),(/  5.500000/),(/  3.214240/),(/  1.741050/),(/  2.386635/),(/33/))
        BVEL_Table(108)=BVEL_Par_Type("I+5",  (/  3.100000/),(/  1.977750/),(/  6.000000/),(/  2.489470/),(/  1.644210/),(/  2.358491/),(/33/))
        BVEL_Table(109)=BVEL_Par_Type("TE+6", (/  6.000000/),(/  1.913430/),(/  5.500000/),(/  2.564060/),(/  1.808760/),(/  2.427184/),(/33/))
        BVEL_Table(110)=BVEL_Par_Type("TE+4", (/  3.396000/),(/  1.952900/),(/  5.500000/),(/  1.671690/),(/  1.752080/),(/  2.493766/),(/33/))
        BVEL_Table(111)=BVEL_Par_Type("CS+1", (/ 11.790000/),(/  2.258990/),(/  6.500000/),(/  0.233070/),(/  3.131210/),(/  2.386635/),(/33/))
        BVEL_Table(112)=BVEL_Par_Type("BA+2", (/ 10.320000/),(/  2.159980/),(/  6.000000/),(/  0.579940/),(/  2.737690/),(/  2.288330/),(/33/))
        BVEL_Table(113)=BVEL_Par_Type("LA+3", (/  9.830000/),(/  2.063920/),(/  5.500000/),(/  1.185870/),(/  2.469890/),(/  2.217295/),(/33/))
        BVEL_Table(114)=BVEL_Par_Type("CE+4", (/  7.867000/),(/  2.028210/),(/  5.500000/),(/  1.484120/),(/  2.198720/),(/  2.257336/),(/33/))
        BVEL_Table(115)=BVEL_Par_Type("CE+3", (/  9.147000/),(/  2.031180/),(/  5.500000/),(/  1.220480/),(/  2.378610/),(/  2.227171/),(/33/))
        BVEL_Table(116)=BVEL_Par_Type("PR+3", (/  9.067000/),(/  2.036520/),(/  5.500000/),(/  1.170410/),(/  2.371130/),(/  2.277904/),(/33/))
        BVEL_Table(117)=BVEL_Par_Type("ND+3", (/  8.647000/),(/  2.024250/),(/  5.500000/),(/  1.132050/),(/  2.330160/),(/  2.336449/),(/33/))
        BVEL_Table(118)=BVEL_Par_Type("SM+3", (/  8.119000/),(/  2.011680/),(/  5.500000/),(/  1.176220/),(/  2.295360/),(/  2.309469/),(/33/))
        BVEL_Table(119)=BVEL_Par_Type("EU+3", (/  7.743000/),(/  2.004690/),(/  5.500000/),(/  1.195450/),(/  2.268880/),(/  2.304147/),(/33/))
        BVEL_Table(120)=BVEL_Par_Type("EU+2", (/ 10.111000/),(/  1.891580/),(/  6.000000/),(/  1.130320/),(/  2.538460/),(/  2.024291/),(/33/))
        BVEL_Table(121)=BVEL_Par_Type("GD+3", (/  8.052000/),(/  1.996540/),(/  5.500000/),(/  1.091610/),(/  2.271900/),(/  2.409639/),(/33/))
        BVEL_Table(122)=BVEL_Par_Type("TB+4", (/  6.000000/),(/  1.962440/),(/  6.000000/),(/  1.701320/),(/  2.385060/),(/  2.024291/),(/33/))
        BVEL_Table(123)=BVEL_Par_Type("TB+3", (/  7.958000/),(/  1.956750/),(/  5.500000/),(/  1.207640/),(/  2.235630/),(/  2.309469/),(/33/))
        BVEL_Table(124)=BVEL_Par_Type("DY+3", (/  7.828000/),(/  1.960290/),(/  5.500000/),(/  1.173500/),(/  2.226890/),(/  2.347418/),(/33/))
        BVEL_Table(125)=BVEL_Par_Type("HO+3", (/  7.500000/),(/  1.970990/),(/  5.500000/),(/  1.121570/),(/  2.211220/),(/  2.409639/),(/33/))
        BVEL_Table(126)=BVEL_Par_Type("ER+3", (/  7.135000/),(/  1.956080/),(/  5.500000/),(/  1.123940/),(/  2.174770/),(/  2.427184/),(/33/))
        BVEL_Table(127)=BVEL_Par_Type("TM+3", (/  6.912000/),(/  1.949010/),(/  5.500000/),(/  1.181380/),(/  2.160420/),(/  2.375297/),(/33/))
        BVEL_Table(128)=BVEL_Par_Type("YB+3", (/  6.875000/),(/  1.928720/),(/  5.500000/),(/  1.219890/),(/  2.142200/),(/  2.347418/),(/33/))
        BVEL_Table(129)=BVEL_Par_Type("LU+3", (/  6.830000/),(/  1.917280/),(/  5.500000/),(/  1.194880/),(/  2.136000/),(/  2.375297/),(/33/))
        BVEL_Table(130)=BVEL_Par_Type("HF+4", (/  7.105000/),(/  1.833610/),(/  6.000000/),(/  1.899920/),(/  1.999640/),(/  2.092050/),(/33/))
        BVEL_Table(131)=BVEL_Par_Type("TA+5", (/  6.090000/),(/  1.868160/),(/  5.500000/),(/  2.366690/),(/  1.855320/),(/  2.057613/),(/33/))
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
                   (/ 0.569, 0.70779, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.940, 0.558,   0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,    34,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
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
    !!----    Practical Considerations in Determining Bond-Valence Parameters
    !!----    Structure and Bonding (2014) 158: 91-128
    !!----    DOI: 10.1007/430_2013_96
    !!----
    !!---- Created: December - 2014, Updated: January 2015 (JRC)
    !!
    Subroutine Set_SBVS_Table()

       if (.not. allocated(sBVS_Table)) allocate(sBVS_Table(sbvs_species_n))

     sBVS_Table(  1)=sBVS_Par_Type( "AG+1", &
                    (/ 1.78239, 1.58298, 1.98819, 2.03699, 2.08036, 2.11006, 2.17409, 2.34487, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.394,   0.498,   0.429,   0.470,   0.530,   0.365,   0.385,   0.420,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.438,   9.000,   5.667,   5.000,   4.692,   3.308,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.000,   5.000,   5.500,   6.000,   5.000,   5.000,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(  2)=sBVS_Par_Type( "AG+2", &
                    (/ 1.72209, 1.63780, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.501,   0.509,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   6.812,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(  3)=sBVS_Par_Type( "AG+3", &
                    (/ 1.84687, 1.71485, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.444,   0.452,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.250,   8.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(  4)=sBVS_Par_Type( "AL+3", &
                    (/ 1.59901, 1.41970, 1.89795, 2.04045, 2.24648, 2.05965, 2.17392, 2.38658, 1.69564, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.424,   0.519,   0.646,   0.687,   0.747,   0.550,   0.582,   0.632,   0.502,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.327,   0.000,   0.000,   0.000,   0.000,   4.077,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   6.500,   7.000,   7.000,   6.000,   6.500,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table(  5)=sBVS_Par_Type( "AS+3", &
                    (/ 1.76706, 1.68300, 2.06175, 0.00000, 0.00000, 2.20841, 2.34035, 2.48812, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.404,   0.406,   0.522,   0.000,   0.000,   0.540,   0.571,   0.616,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.000,   6.600,   3.000,   0.000,   0.000,   2.960,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.000,   6.000,   0.000,   0.000,   6.000,   6.500,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(  6)=sBVS_Par_Type( "AS+5", &
                    (/ 1.76689, 1.60379, 0.00000, 0.00000, 0.00000, 2.27986, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.411,   0.503,   0.000,   0.000,   0.000,   0.534,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.029,   6.000,   0.000,   0.000,   0.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   6.000,   0.000,   0.000,   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(  7)=sBVS_Par_Type( "AU+1", &
                    (/ 1.71819, 0.00000, 0.00000, 0.00000, 0.00000, 2.06404, 2.17613, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.441,   0.000,   0.000,   0.000,   0.000,   0.342,   0.346,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   2.000,   0.000,   0.000,   0.000,   0.000,   2.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   4.500,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(  8)=sBVS_Par_Type( "AU+3", &
                    (/ 1.81761, 1.69700, 2.14532, 2.27911, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.498,   0.507,   0.421,   0.461,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,  10.375,   4.250,   5.429,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(  9)=sBVS_Par_Type( "B+3 ", &
                    (/ 1.35761, 1.19909, 1.67500, 1.82100, 2.00700, 1.78195, 1.91870, 2.09820, 1.46869, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.385,   0.542,   0.669,   0.709,   0.770,   0.573,   0.607,   0.659,   0.315,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.417,   0.000,   0.000,   0.000,   0.000,   3.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.500,   5.000,   6.000,   6.500,   7.000,   6.000,   6.000,   6.500,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      33,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 10)=sBVS_Par_Type( "BA+2", &
                    (/ 2.15998, 2.06475, 2.41041, 2.48600, 2.64950, 2.44239, 2.51975, 2.68807, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.437,   0.428,   0.519,   0.595,   0.655,   0.573,   0.649,   0.649,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/  10.320,   0.000,   0.000,   0.000,   0.000,   8.703,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.500,   7.000,   7.500,   7.000,   7.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 11)=sBVS_Par_Type( "BE+2", &
                    (/ 1.20903, 1.14986, 1.52729, 1.65810, 1.80080, 1.57840, 1.70430, 1.84970, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.541,   0.532,   0.659,   0.699,   0.760,   0.677,   0.708,   0.753,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   7.000,   7.000,   6.500,   7.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 12)=sBVS_Par_Type( "BI+3", &
                    (/ 2.03677, 1.93046, 2.34105, 2.46708, 2.63989, 2.41824, 2.54545, 2.82304, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.414,   0.423,   0.505,   0.545,   0.605,   0.523,   0.554,   0.599,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.058,   7.400,   6.500,   6.667,   6.000,   7.917,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   7.000,   7.500,   6.500,   7.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 13)=sBVS_Par_Type( "BI+5", &
                    (/ 2.04498, 1.86984, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.371,   0.448,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 14)=sBVS_Par_Type( "BR+7", &
                    (/ 1.83658, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.423,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 15)=sBVS_Par_Type( "C+2 ", &
                    (/ 1.41368, 1.15907, 0.00000, 0.00000, 2.03000, 0.00000, 0.00000, 0.00000, 1.48523, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.415,   0.407,   0.000,   0.000,   0.634,   0.000,   0.000,   0.000,   0.524,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   1.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 16)=sBVS_Par_Type( "C+4 ", &
                    (/ 1.39826, 1.37249, 0.00000, 0.00000, 2.09679, 1.84354, 0.00000, 0.00000, 1.55787, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.447,   0.547,   0.000,   0.000,   0.775,   0.579,   0.000,   0.000,   0.529,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.000,   8.333,   0.000,   0.000,   0.000,   3.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   6.000,   0.000,   0.000,   7.000,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 17)=sBVS_Par_Type( "CA+2", &
                    (/ 1.79519, 1.71038, 2.06262, 2.16860, 2.32420, 2.11148, 2.20285, 2.38860, 2.09085, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.476,   0.467,   0.594,   0.634,   0.695,   0.612,   0.643,   0.689,   0.454,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.544,   0.000,   0.000,   0.000,   0.000,   7.250,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   7.000,   7.000,   7.000,   7.000,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 18)=sBVS_Par_Type( "CD+2", &
                    (/ 1.83926, 1.73948, 2.06950, 2.16913, 2.32440, 2.12789, 2.16209, 2.36683, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.407,   0.416,   0.512,   0.553,   0.613,   0.531,   0.561,   0.607,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.176,   7.125,   6.000,   5.273,   5.500,   4.783,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   6.000,   7.000,   5.500,   6.500,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 19)=sBVS_Par_Type( "CE+3", &
                    (/ 2.03118, 1.92656, 0.00000, 0.00000, 0.00000, 2.50404, 2.52663, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.449,   0.458,   0.000,   0.000,   0.000,   0.489,   0.520,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   9.147,  12.750,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 20)=sBVS_Par_Type( "CE+4", &
                    (/ 2.02821, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.443,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.867,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 21)=sBVS_Par_Type( "CL+3", &
                    (/ 1.72265, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.491,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   2.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 22)=sBVS_Par_Type( "CL+5", &
                    (/ 1.69552, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.445,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 23)=sBVS_Par_Type( "CL+7", &
                    (/ 1.67946, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.443,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 24)=sBVS_Par_Type( "CO+1", &
                    (/ 1.29501, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.621,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   2.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 25)=sBVS_Par_Type( "CO+2", &
                    (/ 1.59773, 1.51442, 1.93137, 0.00000, 2.19837, 1.99705, 2.10951, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.451,   0.460,   0.468,   0.000,   0.569,   0.487,   0.518,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.506,   6.000,   4.625,   0.000,   4.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   0.000,   6.000,   5.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 26)=sBVS_Par_Type( "CO+3", &
                    (/ 1.59234, 1.56345, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.434,   0.443,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 27)=sBVS_Par_Type( "CO+4", &
                    (/ 1.79444, 1.70035, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.441,   0.441,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 28)=sBVS_Par_Type( "CR+2", &
                    (/ 1.59356, 1.54102, 2.01499, 2.12678, 2.28404, 1.90930, 0.00000, 2.18257, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.480,   0.489,   0.439,   0.480,   0.540,   0.458,   0.000,   0.534,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.600,   6.778,   6.000,   6.000,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   6.000,   6.500,   5.500,   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,       0,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 29)=sBVS_Par_Type( "CR+3", &
                    (/ 1.66198, 1.58975, 2.00510, 2.13165, 0.00000, 2.04397, 2.14649, 2.30242, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.430,   0.438,   0.489,   0.530,   0.000,   0.508,   0.539,   0.584,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   6.500,   0.000,   5.500,   6.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,       0,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 30)=sBVS_Par_Type( "CR+4", &
                    (/ 1.76095, 1.64899, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.409,   0.418,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.429,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 31)=sBVS_Par_Type( "CR+5", &
                    (/ 1.76781, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.402,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 32)=sBVS_Par_Type( "CR+6", &
                    (/ 1.82471, 0.00000, 2.00700, 2.22400, 2.53730, 2.10300, 2.30000, 2.60900, 1.96859, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.476,   0.000,   0.634,   0.674,   0.735,   0.652,   0.683,   0.729,   0.605,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   6.500,   6.500,   7.000,   6.500,   7.000,   8.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 33)=sBVS_Par_Type( "CS+1", &
                    (/ 2.25899, 2.15318, 2.46037, 2.50818, 2.71194, 2.52283, 2.64858, 2.78226, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.419,   0.428,   0.500,   0.541,   0.601,   0.519,   0.550,   0.595,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/  11.790,   0.000,   0.000,   0.000,   0.000,   7.585,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.500,   6.000,   6.500,   7.000,   8.000,   6.500,   7.000,   7.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 34)=sBVS_Par_Type( "CU+1", &
                    (/ 1.58730, 0.00000, 1.79822, 1.86810, 1.93067, 1.80744, 1.89561, 1.89529, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.341,   0.000,   0.402,   0.442,   0.503,   0.420,   0.451,   0.497,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   2.560,   0.000,   3.929,   3.889,   4.692,   3.308,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   0.000,   5.000,   5.500,   6.000,   5.000,   5.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 35)=sBVS_Par_Type( "CU+2", &
                    (/ 1.57422, 1.49530, 1.90175, 0.00000, 0.00000, 1.94997, 2.02888, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.449,   0.457,   0.470,   0.000,   0.000,   0.488,   0.519,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   2.560,   6.135,   5.692,   0.000,   0.000,   3.880,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   5.500,   0.000,   0.000,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 36)=sBVS_Par_Type( "CU+3", &
                    (/ 1.70964, 1.54272, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.427,   0.437,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 37)=sBVS_Par_Type( "D+1 ", &
                    (/ 0.87141, 0.70779, 1.01040, 1.20000, 1.29000, 1.19164, 1.29750, 1.40390, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.457,   0.558,   0.685,   0.725,   0.786,   0.591,   0.625,   0.678,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   1.923,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   4.500,   5.500,   6.000,   6.500,   5.500,   6.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 38)=sBVS_Par_Type( "DY+3", &
                    (/ 1.96029, 1.85000, 2.28810, 0.00000, 0.00000, 2.34124, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.426,   0.435,   0.493,   0.000,   0.000,   0.512,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.828,   0.000,   7.000,   0.000,   0.000,   7.429,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 39)=sBVS_Par_Type( "ER+3", &
                    (/ 1.95608, 1.84663, 2.25982, 0.00000, 0.00000, 2.32366, 2.43420, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.412,   0.421,   0.507,   0.000,   0.000,   0.525,   0.556,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.135,   8.125,   6.333,   0.000,   0.000,   6.375,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 40)=sBVS_Par_Type( "EU+2", &
                    (/ 1.89158, 1.81420, 2.41447, 2.50011, 2.66292, 2.47173, 0.00000, 0.00000, 2.18740, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.494,   0.503,   0.425,   0.465,   0.525,   0.443,   0.000,   0.000,   0.430,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/  10.111,   7.000,   7.000,   7.333,   7.500,   7.579,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   5.500,   5.500,   6.000,   6.500,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 41)=sBVS_Par_Type( "EU+3", &
                    (/ 2.00469, 1.87763, 2.31998, 0.00000, 0.00000, 2.39179, 0.00000, 0.00000, 2.16270, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.434,   0.443,   0.485,   0.000,   0.000,   0.503,   0.000,   0.000,   0.455,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.743,   8.000,   8.600,   0.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,       0,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 42)=sBVS_Par_Type( "FE+2", &
                    (/ 1.57911, 1.50766, 2.00003, 2.09517, 2.25785, 2.01111, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.480,   0.489,   0.439,   0.479,   0.540,   0.457,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.743,   6.077,   6.000,   6.000,   5.000,   5.400,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.000,   5.500,   6.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 43)=sBVS_Par_Type( "FE+3", &
                    (/ 1.70840, 1.63674, 1.98266, 2.13627, 2.37658, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.420,   0.411,   0.538,   0.578,   0.523,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.733,   5.992,   5.429,   4.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.000,   6.000,   6.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 44)=sBVS_Par_Type( "FE+4", &
                    (/ 1.76559, 1.76220, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.410,   0.417,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 45)=sBVS_Par_Type( "GA+1", &
                    (/ 0.00000, 0.00000, 2.32040, 0.00000, 2.38769, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.440,   0.000,   0.541,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   9.000,   0.000,   9.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   6.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 46)=sBVS_Par_Type( "GA+3", &
                    (/ 1.71606, 1.56039, 1.97794, 2.11352, 2.31400, 2.12458, 2.24723, 2.44491, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.373,   0.451,   0.578,   0.619,   0.679,   0.483,   0.513,   0.558,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.905,   0.000,   0.000,   0.000,   0.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   6.500,   6.500,   7.000,   6.000,   6.500,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 47)=sBVS_Par_Type( "GD+3", &
                    (/ 1.99654, 1.90718, 2.27126, 0.00000, 0.00000, 2.37028, 0.00000, 0.00000, 2.09519, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.415,   0.406,   0.533,   0.000,   0.000,   0.552,   0.000,   0.000,   0.504,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   8.052,   8.333,   7.750,   0.000,   0.000,   7.200,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,       0,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 48)=sBVS_Par_Type( "GE+2", &
                    (/ 0.00000, 0.00000, 2.07494, 0.00000, 2.29129, 0.00000, 0.00000, 2.26937, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.490,   0.000,   0.591,   0.000,   0.000,   0.585,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.500,   0.000,   7.000,   0.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,       0,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 49)=sBVS_Par_Type( "GE+4", &
                    (/ 1.73939, 1.57864, 2.07619, 0.00000, 0.00000, 2.19183, 2.33264, 2.57461, 1.84758, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.396,   0.484,   0.610,   0.000,   0.000,   0.514,   0.545,   0.593,   0.469,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.305,   6.000,   5.000,   0.000,   0.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   6.000,   0.000,   0.000,   6.000,   6.000,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 50)=sBVS_Par_Type( "H+1 ", &
                    (/ 0.87045, 0.70779, 1.01040, 1.20000, 1.29000, 1.19164, 1.29750, 1.40390, 0.91281, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.457,   0.558,   0.685,   0.725,   0.786,   0.591,   0.625,   0.678,   0.540,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   1.923,   2.000,   6.333,   0.000,   0.000,   1.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   4.500,   5.500,   6.000,   6.500,   5.500,   6.000,   6.500,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      35,      34,      34,      34,      35,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 51)=sBVS_Par_Type( "HF+4", &
                    (/ 1.83361, 1.77520, 2.24440, 0.00000, 2.57677, 2.23550, 2.37893, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.478,   0.469,   0.483,   0.000,   0.584,   0.615,   0.645,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.105,   7.250,   6.000,   0.000,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   5.500,   6.500,   0.000,   7.000,   7.000,   7.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 52)=sBVS_Par_Type( "HG+1", &
                    (/ 1.81280, 1.76507, 2.17341, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.465,   0.473,   0.455,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.786,   5.000,   5.333,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 53)=sBVS_Par_Type( "HG+2", &
                    (/ 1.81276, 1.73380, 2.17009, 2.25138, 2.38163, 2.20501, 2.25179, 2.37917, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.465,   0.473,   0.455,   0.495,   0.555,   0.473,   0.504,   0.549,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.966,   6.286,   6.857,   5.167,   4.000,   5.222,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   6.000,   6.000,   5.500,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 54)=sBVS_Par_Type( "HO+3", &
                    (/ 1.97099, 1.84485, 0.00000, 0.00000, 0.00000, 2.33798, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.415,   0.424,   0.000,   0.000,   0.000,   0.522,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.500,   8.333,   0.000,   0.000,   0.000,   6.875,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 55)=sBVS_Par_Type( "HX+1", &
                    (/ 0.78150, 0.66653, 0.96000, 1.15000, 1.24000, 1.09778, 1.21280, 1.31710, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.457,   0.558,   0.685,   0.725,   0.786,   0.591,   0.625,   0.678,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   1.923,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   4.500,   5.500,   6.000,   6.500,   5.500,   6.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 56)=sBVS_Par_Type( "I+5 ", &
                    (/ 1.97775, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.424,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.100,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 57)=sBVS_Par_Type( "I+7 ", &
                    (/ 1.92274, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.419,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.800,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 58)=sBVS_Par_Type( "IN+1", &
                    (/ 0.00000, 0.00000, 2.43863, 2.46665, 2.50769, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.413,   0.454,   0.514,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   9.000,   8.700,   7.833,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   6.000,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 59)=sBVS_Par_Type( "IN+3", &
                    (/ 1.90305, 1.75728, 2.13474, 2.28936, 2.47773, 2.30797, 2.43051, 2.62278, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.353,   0.421,   0.548,   0.589,   0.667,   0.456,   0.484,   0.534,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.024,   6.941,   6.000,   4.400,   4.000,   5.453,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   6.500,   6.500,   7.000,   6.000,   6.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 60)=sBVS_Par_Type( "IR+3", &
                    (/ 0.00000, 0.00000, 2.05823, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.402,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 61)=sBVS_Par_Type( "IR+4", &
                    (/ 1.83233, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.436,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 62)=sBVS_Par_Type( "IR+5", &
                    (/ 1.89791, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.479,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 63)=sBVS_Par_Type( "K+1 ", &
                    (/ 1.94117, 1.83331, 2.07673, 2.17206, 2.28880, 2.16803, 2.29712, 2.41802, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.436,   0.428,   0.555,   0.595,   0.655,   0.573,   0.604,   0.649,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   8.846,   0.000,   0.000,   0.000,   0.000,   8.400,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   5.500,   6.500,   7.000,   7.000,   6.500,   7.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 64)=sBVS_Par_Type( "LA+3", &
                    (/ 2.06392, 1.96577, 2.43099, 0.00000, 0.00000, 2.40330, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.451,   0.443,   0.459,   0.000,   0.000,   0.588,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   9.830,  10.286,   9.000,   0.000,   0.000,   7.944,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 65)=sBVS_Par_Type( "LI+1", &
                    (/ 1.17096, 1.08674, 1.39892, 1.51288, 1.64829, 1.46652, 1.62021, 1.71028, 1.15430, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.516,   0.508,   0.634,   0.675,   0.735,   0.653,   0.684,   0.729,   0.637,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.021,   5.538,   5.867,   5.846,   5.818,   4.385,   4.400,   4.670,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   7.000,   7.000,   6.000,   7.000,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 66)=sBVS_Par_Type( "LU+3", &
                    (/ 1.91728, 1.80738, 2.28400, 0.00000, 0.00000, 2.28119, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.421,   0.412,   0.433,   0.000,   0.000,   0.557,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.830,   7.500,   6.000,   0.000,   0.000,   6.077,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 67)=sBVS_Par_Type( "MG+2", &
                    (/ 1.48398, 1.40956, 1.78276, 1.89059, 2.04500, 1.82372, 1.70430, 1.94970, 1.63729, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.512,   0.504,   0.630,   0.671,   0.731,   0.649,   0.708,   0.680,   0.601,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.897,   0.000,   0.000,   0.000,   0.000,   5.556,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   7.000,   7.000,   6.500,   7.000,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 68)=sBVS_Par_Type( "MN+2", &
                    (/ 1.62758, 1.56491, 1.99681, 2.08701, 2.26045, 2.15468, 2.24751, 2.43113, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.481,   0.474,   0.488,   0.528,   0.588,   0.406,   0.431,   0.470,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.910,   6.242,   5.905,   5.750,   6.000,   5.308,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   6.000,   6.500,   5.500,   5.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 69)=sBVS_Par_Type( "MN+3", &
                    (/ 1.68993, 1.59633, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.437,   0.446,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.862,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 70)=sBVS_Par_Type( "MN+4", &
                    (/ 1.73272, 1.61606, 2.06541, 0.00000, 0.00000, 2.30057, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.402,   0.410,   0.517,   0.000,   0.000,   0.536,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.923,   6.000,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.000,   6.500,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 71)=sBVS_Par_Type( "MN+5", &
                    (/ 1.78879, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.413,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 72)=sBVS_Par_Type( "MN+6", &
                    (/ 1.82018, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.416,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 73)=sBVS_Par_Type( "MN+7", &
                    (/ 1.87362, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.520,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 74)=sBVS_Par_Type( "MO+2", &
                    (/ 0.00000, 0.00000, 2.05175, 2.22218, 2.35705, 2.07169, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.441,   0.401,   0.461,   0.422,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.000,   5.000,   5.000,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.500,   5.500,   6.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      35,      34,      34,      35,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 75)=sBVS_Par_Type( "MO+3", &
                    (/ 1.78933, 1.73819, 2.08941, 2.19147, 0.00000, 2.06172, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.418,   0.427,   0.501,   0.541,   0.000,   0.519,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.700,   6.000,   6.286,   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      35,      35,      35,       0,      35,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 76)=sBVS_Par_Type( "MO+4", &
                    (/ 1.72390, 0.00000, 2.12830, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.562,   0.000,   0.558,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.500,   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      35,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 77)=sBVS_Par_Type( "MO+5", &
                    (/ 1.84760, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.482,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.980,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 78)=sBVS_Par_Type( "MO+6", &
                    (/ 1.90934, 0.00000, 2.08880, 2.30370, 2.61910, 2.18700, 2.38240, 2.66000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.391,   0.000,   0.622,   0.663,   0.723,   0.641,   0.672,   0.717,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.764,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   0.000,   6.500,   6.500,   7.000,   6.500,   6.500,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 79)=sBVS_Par_Type( "N+3 ", &
                    (/ 1.40795, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.448,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   2.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 80)=sBVS_Par_Type( "N+5 ", &
                    (/ 1.46267, 1.40717, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.450,   0.550,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 81)=sBVS_Par_Type( "NA+1", &
                    (/ 1.56225, 1.42885, 1.69489, 1.78228, 1.94353, 1.83088, 1.89880, 2.03011, 1.60725, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.482,   0.474,   0.601,   0.641,   0.701,   0.619,   0.650,   0.695,   0.571,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.520,   0.000,   0.000,   0.000,   0.000,   5.772,   5.210,   5.520,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   5.500,   6.500,   7.000,   7.000,   6.500,   7.000,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 82)=sBVS_Par_Type( "NB+2", &
                    (/ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 2.39087, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.448,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,       0,       0,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 83)=sBVS_Par_Type( "NB+3", &
                    (/ 1.74581, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.501,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 84)=sBVS_Par_Type( "NB+4", &
                    (/ 1.78543, 1.75844, 0.00000, 0.00000, 0.00000, 2.28012, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.526,   0.534,   0.000,   0.000,   0.000,   0.412,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   0.000,   0.000,   0.000,   5.750,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   0.000,   0.000,   0.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 85)=sBVS_Par_Type( "NB+5", &
                    (/ 1.86588, 1.76034, 2.20986, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 2.02880, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.498,   0.489,   0.616,   0.000,   0.000,   0.000,   0.000,   0.000,   0.480,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.044,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,       0,       0,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 86)=sBVS_Par_Type( "ND+3", &
                    (/ 2.02425, 1.91099, 0.00000, 0.00000, 0.00000, 2.44624, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.428,   0.437,   0.000,   0.000,   0.000,   0.510,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   8.647,  10.000,   0.000,   0.000,   0.000,   7.700,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 87)=sBVS_Par_Type( "NH+1", &
                    (/ 2.03380, 1.97933, 2.15362, 2.23280, 2.40576, 2.24264, 2.34800, 2.48100, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.442,   0.434,   0.561,   0.601,   0.661,   0.579,   0.610,   0.655,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.467,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   5.500,   6.500,   6.500,   7.000,   7.000,   7.000,   7.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 88)=sBVS_Par_Type( "NI+2", &
                    (/ 1.55920, 1.49655, 1.87666, 1.96729, 0.00000, 1.85695, 1.89863, 2.11175, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.443,   0.452,   0.476,   0.516,   0.000,   0.494,   0.525,   0.571,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.933,   6.000,   6.000,   6.000,   0.000,   4.700,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.000,   5.500,   6.000,   0.000,   5.500,   6.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,       0,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 89)=sBVS_Par_Type( "NI+3", &
                    (/ 1.64888, 1.56869, 0.00000, 0.00000, 0.00000, 1.95763, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.414,   0.423,   0.000,   0.000,   0.000,   0.523,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.000,   0.000,   0.000,   0.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 90)=sBVS_Par_Type( "NI+4", &
                    (/ 1.72159, 1.61746, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.402,   0.411,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 91)=sBVS_Par_Type( "OS+4", &
                    (/ 1.75302, 0.00000, 2.15572, 2.30615, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.498,   0.000,   0.421,   0.461,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 92)=sBVS_Par_Type( "OS+5", &
                    (/ 0.00000, 1.75815, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.488,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 93)=sBVS_Par_Type( "OS+6", &
                    (/ 1.93192, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.463,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 94)=sBVS_Par_Type( "OS+7", &
                    (/ 1.95775, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.479,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 95)=sBVS_Par_Type( "OS+8", &
                    (/ 1.97728, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.512,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.333,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 96)=sBVS_Par_Type( "P+3 ", &
                    (/ 1.51555, 0.00000, 0.00000, 0.00000, 0.00000, 1.96089, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.402,   0.000,   0.000,   0.000,   0.000,   0.536,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.000,   0.000,   0.000,   0.000,   0.000,   3.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.500,   0.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 97)=sBVS_Par_Type( "P+5 ", &
                    (/ 1.62038, 1.47610, 2.01093, 2.21050, 2.49792, 2.14917, 2.43695, 0.00000, 1.73267, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.437,   0.535,   0.662,   0.703,   0.763,   0.566,   0.600,   0.000,   0.517,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   6.000,   4.000,   0.000,   4.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   7.000,   7.000,   7.000,   7.500,   6.000,   7.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table( 98)=sBVS_Par_Type( "PB+2", &
                    (/ 2.01825, 1.90916, 2.34641, 2.43160, 2.56880, 2.38758, 2.46755, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.433,   0.453,   0.474,   0.515,   0.575,   0.493,   0.524,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.541,  10.588,   7.833,   7.714,   6.429,   7.551,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   6.000,   5.500,   6.000,   7.000,   6.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table( 99)=sBVS_Par_Type( "PB+4", &
                    (/ 2.03293, 1.93428, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.354,   0.424,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.740,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(100)=sBVS_Par_Type( "PD+2", &
                    (/ 1.62359, 1.56885, 1.99818, 2.10183, 0.00000, 2.02559, 2.13063, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.498,   0.506,   0.421,   0.462,   0.000,   0.439,   0.471,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   6.897,   5.250,   4.667,   0.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.000,   5.500,   0.000,   5.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(101)=sBVS_Par_Type( "PD+4", &
                    (/ 1.80500, 1.68582, 2.12662, 0.00000, 0.00000, 0.00000, 2.31115, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.449,   0.457,   0.470,   0.000,   0.000,   0.000,   0.520,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.333,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(102)=sBVS_Par_Type( "PR+3", &
                    (/ 2.03652, 1.91752, 2.48619, 0.00000, 0.00000, 2.46496, 2.52676, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.439,   0.448,   0.387,   0.000,   0.000,   0.499,   0.529,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   9.067,   9.000,   8.000,   0.000,   0.000,   8.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   0.000,   0.000,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(103)=sBVS_Par_Type( "PT+2", &
                    (/ 1.51205, 0.00000, 1.97699, 0.00000, 0.00000, 2.02450, 2.16525, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.574,   0.000,   0.456,   0.000,   0.000,   0.438,   0.407,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   4.600,   0.000,   0.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   5.500,   0.000,   0.000,   5.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      34,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(104)=sBVS_Par_Type( "PT+3", &
                    (/ 1.66559, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.546,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(105)=sBVS_Par_Type( "PT+4", &
                    (/ 1.82198, 1.69132, 2.14110, 0.00000, 0.00000, 2.19999, 2.30943, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.479,   0.488,   0.440,   0.000,   0.000,   0.458,   0.489,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(106)=sBVS_Par_Type( "PT+5", &
                    (/ 0.00000, 1.78573, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.424,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(107)=sBVS_Par_Type( "RB+1", &
                    (/ 2.08597, 1.99137, 2.26585, 2.34362, 2.46445, 2.30817, 2.40364, 2.43463, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.413,   0.404,   0.531,   0.572,   0.632,   0.550,   0.580,   0.626,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/  10.020,   0.000,   0.000,   0.000,   0.000,   6.552,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.500,   5.500,   6.500,   7.000,   7.000,   6.500,   7.000,   7.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(108)=sBVS_Par_Type( "RE+3", &
                    (/ 0.00000, 0.00000, 2.15781, 2.30419, 0.00000, 2.20710, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.420,   0.422,   0.000,   0.401,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   9.000,   4.750,   0.000,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.500,   5.500,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(109)=sBVS_Par_Type( "RE+4", &
                    (/ 1.78237, 0.00000, 2.18513, 0.00000, 0.00000, 2.26800, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.507,   0.000,   0.412,   0.000,   0.000,   0.431,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   6.000,   0.000,   0.000,   5.400,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   5.500,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(110)=sBVS_Par_Type( "RE+5", &
                    (/ 1.82664, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.479,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(111)=sBVS_Par_Type( "RE+6", &
                    (/ 1.91007, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.498,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(112)=sBVS_Par_Type( "RE+7", &
                    (/ 1.97792, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.508,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.098,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(113)=sBVS_Par_Type( "RH+3", &
                    (/ 1.67013, 1.63692, 2.03722, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.478,   0.488,   0.440,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(114)=sBVS_Par_Type( "RH+4", &
                    (/ 1.77675, 1.68686, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.403,   0.406,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(115)=sBVS_Par_Type( "RU+3", &
                    (/ 1.72066, 0.00000, 1.99520, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.430,   0.000,   0.489,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(116)=sBVS_Par_Type( "RU+4", &
                    (/ 1.79363, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.449,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(117)=sBVS_Par_Type( "RU+5", &
                    (/ 1.87442, 1.77634, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.436,   0.445,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(118)=sBVS_Par_Type( "RU+6", &
                    (/ 1.92579, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.425,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(119)=sBVS_Par_Type( "S+4 ", &
                    (/ 1.64282, 0.00000, 2.05269, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.427,   0.000,   0.545,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(120)=sBVS_Par_Type( "S+6 ", &
                    (/ 1.64220, 1.55000, 0.00000, 0.00000, 0.00000, 2.20000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.441,   0.441,   0.000,   0.000,   0.000,   0.571,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.000,   0.000,   0.000,   0.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(121)=sBVS_Par_Type( "SB+3", &
                    (/ 1.92036, 1.83448, 2.30831, 2.38165, 2.59905, 2.36947, 2.47586, 2.73139, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.423,   0.431,   0.496,   0.537,   0.597,   0.515,   0.546,   0.591,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   7.000,   6.400,   7.500,   6.200,   4.337,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   6.000,   7.000,   7.500,   6.000,   6.500,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(122)=sBVS_Par_Type( "SB+5", &
                    (/ 1.89768, 1.76603, 2.22320, 0.00000, 0.00000, 2.44254, 2.58978, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.400,   0.489,   0.616,   0.000,   0.000,   0.520,   0.551,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   0.000,   0.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   6.000,   6.500,   0.000,   0.000,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,      34,       0,       0,      33,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(123)=sBVS_Par_Type( "SC+3", &
                    (/ 1.73220, 1.65325, 2.11330, 0.00000, 0.00000, 2.09450, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.494,   0.485,   0.612,   0.000,   0.000,   0.631,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.255,   6.250,   6.000,   0.000,   0.000,   6.048,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   6.000,   6.500,   0.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(124)=sBVS_Par_Type( "SE+4", &
                    (/ 1.80095, 0.00000, 2.15849, 2.32723, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.427,   0.000,   0.545,   0.585,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.000,   0.000,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   6.500,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(125)=sBVS_Par_Type( "SE+6", &
                    (/ 1.79866, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.416,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(126)=sBVS_Par_Type( "SI+4", &
                    (/ 1.60817, 1.45654, 0.00000, 0.00000, 0.00000, 2.09975, 2.23551, 0.00000, 1.69319, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.432,   0.529,   0.000,   0.000,   0.000,   0.560,   0.529,   0.000,   0.511,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.100,   0.000,   0.000,   0.000,   0.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   0.000,   0.000,   0.000,   6.500,   5.500,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,       0,       0,       0,      34,      34,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table(127)=sBVS_Par_Type( "SM+3", &
                    (/ 2.01168, 1.88844, 0.00000, 0.00000, 0.00000, 2.41460, 0.00000, 2.71979, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.433,   0.441,   0.000,   0.000,   0.000,   0.505,   0.000,   0.581,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   8.119,   8.429,   0.000,   0.000,   0.000,   7.667,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   6.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,       0,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(128)=sBVS_Par_Type( "SN+2", &
                    (/ 1.87499, 1.80141, 2.28757, 2.37369, 2.54376, 2.32497, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.458,   0.467,   0.461,   0.501,   0.562,   0.479,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.325,   0.000,   0.000,   0.000,   0.000,   4.600,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   6.000,   5.500,   6.000,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(129)=sBVS_Par_Type( "SN+4", &
                    (/ 1.89019, 1.75558, 2.17643, 0.00000, 0.00000, 2.36256, 2.50725, 2.73112, 1.99361, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.379,   0.461,   0.588,   0.000,   0.000,   0.492,   0.522,   0.568,   0.449,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.069,   6.000,   5.889,   0.000,   0.000,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   6.500,   0.000,   0.000,   6.000,   6.000,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table(130)=sBVS_Par_Type( "SR+2", &
                    (/ 1.95311, 1.86817, 2.20197, 2.31250, 2.32420, 2.26765, 2.36718, 2.54050, 2.09477, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.455,   0.446,   0.573,   0.614,   0.695,   0.591,   0.622,   0.668,   0.544,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   9.400,   0.000,   0.000,   0.000,   0.000,   8.074,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   7.000,   7.000,   7.000,   7.000,   7.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table(131)=sBVS_Par_Type( "TA+2", &
                    (/ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 2.20689, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.593,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,       0,       0,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(132)=sBVS_Par_Type( "TA+3", &
                    (/ 0.00000, 1.52314, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 1.76100, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.555,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.457,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,      34,       0,       0,       0,       0,       0,       0,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table(133)=sBVS_Par_Type( "TA+4", &
                    (/ 1.75632, 0.00000, 0.00000, 0.00000, 0.00000, 2.28935, 2.36128, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.546,   0.000,   0.000,   0.000,   0.000,   0.410,   0.422,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   5.667,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   5.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(134)=sBVS_Par_Type( "TA+5", &
                    (/ 1.86816, 1.79011, 2.21442, 0.00000, 0.00000, 2.27860, 2.51263, 2.71824, 2.00405, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.486,   0.477,   0.604,   0.000,   0.000,   0.622,   0.654,   0.699,   0.463,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.090,   6.400,   6.000,   0.000,   0.000,   7.167,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   0.000,   0.000,   6.500,   7.000,   7.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,      34,      34,       0,       0,       0,       0,       0/) )
     sBVS_Table(135)=sBVS_Par_Type( "TB+3", &
                    (/ 1.95675, 1.83798, 2.29526, 0.00000, 0.00000, 2.36384, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.433,   0.442,   0.486,   0.000,   0.000,   0.505,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.958,   9.333,   9.000,   0.000,   0.000,   7.600,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(136)=sBVS_Par_Type( "TB+4", &
                    (/ 1.96244, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.494,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(137)=sBVS_Par_Type( "TC+4", &
                    (/ 0.00000, 0.00000, 2.15861, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.470,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(138)=sBVS_Par_Type( "TC+7", &
                    (/ 1.97036, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.514,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(139)=sBVS_Par_Type( "TE+2", &
                    (/ 1.39168, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.612,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   2.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(140)=sBVS_Par_Type( "TE+4", &
                    (/ 1.95290, 0.00000, 2.30563, 2.46108, 2.67663, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.401,   0.000,   0.518,   0.559,   0.619,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   3.396,   0.000,   6.300,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   6.500,   7.000,   7.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(141)=sBVS_Par_Type( "TE+6", &
                    (/ 1.91343, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.412,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(142)=sBVS_Par_Type( "TH+4", &
                    (/ 2.04983, 2.40121, 0.00000, 0.00000, 0.00000, 2.45709, 2.54850, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.486,   0.561,   0.000,   0.000,   0.000,   0.579,   0.611,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   8.941,   0.000,   0.000,   0.000,   0.000,   8.667,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   7.000,   0.000,   0.000,   0.000,   7.000,   7.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(143)=sBVS_Par_Type( "TI+2", &
                    (/ 0.00000, 0.00000, 2.03427, 2.15487, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.429,   0.470,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.667,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(144)=sBVS_Par_Type( "TI+3", &
                    (/ 1.69766, 1.63180, 2.13089, 0.00000, 0.00000, 2.11813, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.460,   0.468,   0.459,   0.000,   0.000,   0.478,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(145)=sBVS_Par_Type( "TI+4", &
                    (/ 1.72394, 1.66170, 2.12577, 0.00000, 0.00000, 2.17452, 0.00000, 2.41793, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.503,   0.449,   0.507,   0.000,   0.000,   0.640,   0.000,   0.716,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   0.000,   0.000,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.000,   0.000,   0.000,   7.000,   0.000,   7.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,       0,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(146)=sBVS_Par_Type( "TL+1", &
                    (/ 1.91752, 1.82737, 2.36792, 2.42772, 2.48689, 2.38260, 2.39935, 2.50518, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.483,   0.491,   0.491,   0.477,   0.537,   0.455,   0.508,   0.531,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   8.030,  12.625,   8.875,   8.400,   8.200,   7.341,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   6.000,   7.000,   6.000,   6.500,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(147)=sBVS_Par_Type( "TL+3", &
                    (/ 2.06297, 1.85700, 2.22567, 2.35992, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.338,   0.325,   0.514,   0.555,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.220,   7.250,   6.222,   4.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.500,   6.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(148)=sBVS_Par_Type( "TM+3", &
                    (/ 1.94901, 1.81795, 0.00000, 0.00000, 0.00000, 2.31070, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.421,   0.430,   0.000,   0.000,   0.000,   0.516,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.912,   8.000,   0.000,   0.000,   0.000,   6.312,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(149)=sBVS_Par_Type( "U+6 ", &
                    (/ 2.02317, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.527,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.987,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(150)=sBVS_Par_Type( "V+2 ", &
                    (/ 1.59976, 1.53303, 2.03267, 2.10218, 0.00000, 1.89815, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.477,   0.485,   0.442,   0.483,   0.000,   0.461,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,      34,      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(151)=sBVS_Par_Type( "V+3 ", &
                    (/ 1.67799, 1.61581, 2.06334, 0.00000, 0.00000, 2.02622, 2.13686, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.439,   0.448,   0.480,   0.000,   0.000,   0.499,   0.530,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   6.000,   6.000,   0.000,   0.000,   5.920,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   5.500,   0.000,   0.000,   6.500,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(152)=sBVS_Par_Type( "V+4 ", &
                    (/ 1.74932, 1.65629, 0.00000, 0.00000, 0.00000, 2.16760, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.426,   0.430,   0.000,   0.000,   0.000,   0.512,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.738,   6.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.000,   0.000,   0.000,   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,      34,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(153)=sBVS_Par_Type( "V+5 ", &
                    (/ 1.79445, 0.00000, 0.00000, 0.00000, 0.00000, 2.29252, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.510,   0.000,   0.000,   0.000,   0.000,   0.647,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.166,   0.000,   0.000,   0.000,   0.000,   4.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   0.000,   0.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(154)=sBVS_Par_Type( "W+2 ", &
                    (/ 0.00000, 0.00000, 0.00000, 1.95717, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.000,   0.647,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(155)=sBVS_Par_Type( "W+3 ", &
                    (/ 0.00000, 0.00000, 2.12237, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.428,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.667,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(156)=sBVS_Par_Type( "W+4 ", &
                    (/ 1.74558, 0.00000, 0.00000, 2.35170, 0.00000, 2.22672, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.520,   0.000,   0.000,   0.439,   0.000,   0.417,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   0.000,   6.500,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,       0,      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(157)=sBVS_Par_Type( "W+5 ", &
                    (/ 1.81975, 0.00000, 2.24489, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.498,   0.000,   0.420,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(158)=sBVS_Par_Type( "W+6 ", &
                    (/ 1.90641, 0.00000, 2.21497, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.401,   0.000,   0.617,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.688,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(159)=sBVS_Par_Type( "XE+2", &
                    (/ 0.00000, 1.77407, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.569,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   2.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(160)=sBVS_Par_Type( "XE+4", &
                    (/ 0.00000, 1.85567, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.529,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   5.286,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   6.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(161)=sBVS_Par_Type( "XE+6", &
                    (/ 0.00000, 1.86157, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.431,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   8.231,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(162)=sBVS_Par_Type( "Y+3 ", &
                    (/ 1.90384, 1.89262, 2.28867, 2.40922, 0.00000, 2.25923, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.478,   0.379,   0.483,   0.522,   0.000,   0.615,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   7.285,   7.824,   6.400,   6.000,   0.000,   6.417,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   6.500,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(163)=sBVS_Par_Type( "YB+2", &
                    (/ 1.63254, 0.00000, 2.29464, 0.00000, 0.00000, 2.40079, 2.62618, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.510,   0.000,   0.409,   0.000,   0.000,   0.427,   0.458,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   8.000,   0.000,   6.500,   0.000,   0.000,   7.750,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   0.000,   5.500,   0.000,   0.000,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      34,       0,      34,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(164)=sBVS_Par_Type( "YB+3", &
                    (/ 1.92872, 1.80974, 0.00000, 0.00000, 0.00000, 2.31245, 2.42113, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.426,   0.435,   0.000,   0.000,   0.000,   0.511,   0.542,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.875,   7.500,   0.000,   0.000,   0.000,   6.125,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   0.000,   0.000,   0.000,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,       0,       0,       0,      34,      34,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(165)=sBVS_Par_Type( "ZN+2", &
                    (/ 1.65344, 1.57447, 1.88195, 1.98741, 2.15151, 1.94373, 2.00075, 2.15126, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.403,   0.406,   0.521,   0.562,   0.622,   0.540,   0.571,   0.616,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   4.718,   6.269,   4.000,   4.000,   4.500,   4.481,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.000,   5.000,   6.000,   6.000,   6.500,   6.000,   6.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,      34,      34,      34,      34,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(166)=sBVS_Par_Type( "ZR+2", &
                    (/ 0.00000, 0.00000, 2.12072, 0.00000, 2.49296, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.478,   0.000,   0.424,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   5.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   6.000,   0.000,   5.500,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(167)=sBVS_Par_Type( "ZR+3", &
                    (/ 0.00000, 0.00000, 0.00000, 0.00000, 2.56625, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   0.472,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   0.000,   0.000,   0.000,   0.000,   6.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/       0,       0,       0,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0,       0/) )
     sBVS_Table(168)=sBVS_Par_Type( "ZR+4", &
                    (/ 1.84505, 1.83174, 2.18437, 2.32815, 0.00000, 2.26265, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/), &
                    (/   0.490,   0.388,   0.608,   0.648,   0.000,   0.626,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   6.765,   7.146,   6.000,   6.000,   0.000,   6.333,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/   5.500,   5.500,   6.500,   7.000,   0.000,   7.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000/), &
                    (/      33,      34,      34,      34,       0,      34,       0,       0,       0,       0,       0,       0,       0,       0/) )
    End Subroutine Set_SBVS_Table

    Subroutine Set_Common_Oxidation_States_Table()

      If (.Not. Allocated(Common_OxStates_Table)) Allocate(Common_OxStates_Table(8,108))

      Common_OxStates_Table(:,  1) = (/ -1 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,  2) = (/  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,  3) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,  4) = (/  2 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,  5) = (/  3 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,  6) = (/ -4 , -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 /)
      Common_OxStates_Table(:,  7) = (/ -3 ,  3 ,  5 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,  8) = (/ -2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,  9) = (/ -1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 10) = (/  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 11) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 12) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 13) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 14) = (/ -4 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 15) = (/ -3 ,  3 ,  5 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 16) = (/ -2 ,  2 ,  4 ,  6 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 17) = (/ -1 ,  1 ,  3 ,  5 ,  7 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 18) = (/  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 19) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 20) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 21) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 22) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 23) = (/  5 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 24) = (/  3 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 25) = (/  2 ,  4 ,  7 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 26) = (/  2 ,  3 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 27) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 28) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 29) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 30) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 31) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 32) = (/ -4 ,  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 33) = (/ -3 ,  3 ,  5 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 34) = (/ -2 ,  2 ,  4 ,  6 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 35) = (/ -1 ,  1 ,  3 ,  5 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 36) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 37) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 38) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 39) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 40) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 41) = (/  5 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 42) = (/  4 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 43) = (/  4 ,  7 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 44) = (/  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 45) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 46) = (/  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 47) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 48) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 49) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 50) = (/ -4 ,  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 51) = (/ -3 ,  3 ,  5 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 52) = (/ -2 ,  2 ,  4 ,  6 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 53) = (/ -1 ,  1 ,  3 ,  5 ,  7 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 54) = (/  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 55) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 56) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 57) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 58) = (/  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 59) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 60) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 61) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 62) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 63) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 64) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 65) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 66) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 67) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 68) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 69) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 70) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 71) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 72) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 73) = (/  5 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 74) = (/  4 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 75) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 76) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 77) = (/  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 78) = (/  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 79) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 80) = (/  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 81) = (/  1 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 82) = (/  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 83) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 84) = (/ -2 ,  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 85) = (/ -1 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 86) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 87) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 88) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 89) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 90) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 91) = (/  5 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 92) = (/  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 93) = (/  5 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 94) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 95) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 96) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 97) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 98) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:, 99) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,100) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,101) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,102) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,103) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,104) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,105) = (/  5 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,106) = (/  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,107) = (/  7 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      Common_OxStates_Table(:,108) = (/  8 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)

    End Subroutine Set_Common_Oxidation_States_Table

    Subroutine Set_Oxidation_States_Table()

      If (.Not. Allocated(OxStates_Table)) Allocate(OxStates_Table(11,108))

      OxStates_Table(:,  1) = (/ -1 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,  2) = (/  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,  3) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,  4) = (/  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,  5) = (/ -5 , -1 ,  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,  6) = (/ -4 , -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,  7) = (/ -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,  8) = (/ -2 , -1 ,  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,  9) = (/ -1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 10) = (/  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 11) = (/ -1 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 12) = (/  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 13) = (/ -2 , -1 ,  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 14) = (/ -4 , -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 15) = (/ -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 16) = (/ -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 17) = (/ -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 18) = (/  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 19) = (/ -1 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 20) = (/ -1 ,  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 21) = (/  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 22) = (/ -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 23) = (/ -3 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 24) = (/ -4 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 /)
      OxStates_Table(:, 25) = (/ -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  0 /)
      OxStates_Table(:, 26) = (/ -4 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 /)
      OxStates_Table(:, 27) = (/ -3 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 28) = (/ -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 29) = (/ -2 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 30) = (/ -2 ,  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 31) = (/ -5 , -4 , -2 , -1 ,  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 32) = (/ -4 , -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 33) = (/ -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 34) = (/ -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 35) = (/ -1 ,  1 ,  3 ,  4 ,  5 ,  7 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 36) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 37) = (/ -1 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 38) = (/  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 39) = (/  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 40) = (/ -2 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 41) = (/ -3 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 42) = (/ -4 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 /)
      OxStates_Table(:, 43) = (/ -3 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  0 ,  0 /)
      OxStates_Table(:, 44) = (/ -4 , -2 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  0 /)
      OxStates_Table(:, 45) = (/ -3 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 46) = (/  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 47) = (/ -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 48) = (/ -2 ,  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 49) = (/ -5 , -2 , -1 ,  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 50) = (/ -4 , -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 51) = (/ -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 52) = (/ -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 53) = (/ -1 ,  1 ,  3 ,  4 ,  5 ,  6 ,  7 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 54) = (/  2 ,  4 ,  6 ,  8 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 55) = (/ -1 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 56) = (/  1 ,  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 57) = (/  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 58) = (/  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 59) = (/  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 60) = (/  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 61) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 62) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 63) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 64) = (/  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 65) = (/  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 66) = (/  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 67) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 68) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 69) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 70) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 71) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 72) = (/ -2 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 73) = (/ -3 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 74) = (/ -4 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 /)
      OxStates_Table(:, 75) = (/ -3 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  0 ,  0 /)
      OxStates_Table(:, 76) = (/ -4 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 /)
      OxStates_Table(:, 77) = (/ -3 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  9 /)
      OxStates_Table(:, 78) = (/ -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 /)
      OxStates_Table(:, 79) = (/ -3 , -2 , -1 ,  1 ,  2 ,  3 ,  5 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 80) = (/ -2 ,  1 ,  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 81) = (/ -5 , -2 , -1 ,  1 ,  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 82) = (/ -4 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 83) = (/ -3 , -2 , -1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 84) = (/ -2 ,  2 ,  4 ,  5 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 85) = (/ -1 ,  1 ,  3 ,  5 ,  7 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 86) = (/  2 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 87) = (/  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 88) = (/  2 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 89) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 90) = (/  1 ,  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 91) = (/  2 ,  3 ,  4 ,  5 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 92) = (/  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 93) = (/  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 94) = (/  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 95) = (/  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 96) = (/  2 ,  3 ,  4 ,  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 97) = (/  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 98) = (/  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:, 99) = (/  2 ,  3 ,  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,100) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,101) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,102) = (/  2 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,103) = (/  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,104) = (/  4 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,105) = (/  5 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,106) = (/  6 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,107) = (/  7 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)
      OxStates_Table(:,108) = (/  8 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 /)

    End Subroutine Set_Oxidation_States_Table

    Subroutine Set_Pauling_Electronegativity()

      If (.Not. Allocated(PaulingX)) Allocate(PaulingX(108))

      PaulingX(  1) = 2.20
      PaulingX(  2) = 0.00
      PaulingX(  3) = 0.98
      PaulingX(  4) = 1.57
      PaulingX(  5) = 2.04
      PaulingX(  6) = 2.55
      PaulingX(  7) = 3.04
      PaulingX(  8) = 3.44
      PaulingX(  9) = 3.98
      PaulingX( 10) = 0.00
      PaulingX( 11) = 0.93
      PaulingX( 12) = 1.31
      PaulingX( 13) = 1.61
      PaulingX( 14) = 1.90
      PaulingX( 15) = 2.19
      PaulingX( 16) = 2.58
      PaulingX( 17) = 3.16
      PaulingX( 18) = 0.00
      PaulingX( 19) = 0.82
      PaulingX( 20) = 1.00
      PaulingX( 21) = 1.36
      PaulingX( 22) = 1.54
      PaulingX( 23) = 1.63
      PaulingX( 24) = 1.66
      PaulingX( 25) = 1.55
      PaulingX( 26) = 1.83
      PaulingX( 27) = 1.88
      PaulingX( 28) = 1.91
      PaulingX( 29) = 1.90
      PaulingX( 30) = 1.65
      PaulingX( 31) = 1.81
      PaulingX( 32) = 2.01
      PaulingX( 33) = 2.18
      PaulingX( 34) = 2.55
      PaulingX( 35) = 2.96
      PaulingX( 36) = 3.00
      PaulingX( 37) = 0.82
      PaulingX( 38) = 0.95
      PaulingX( 39) = 1.22
      PaulingX( 40) = 1.33
      PaulingX( 41) = 1.60
      PaulingX( 42) = 2.16
      PaulingX( 43) = 1.90
      PaulingX( 44) = 2.20
      PaulingX( 45) = 2.28
      PaulingX( 46) = 2.20
      PaulingX( 47) = 1.93
      PaulingX( 48) = 1.69
      PaulingX( 49) = 1.78
      PaulingX( 50) = 1.96
      PaulingX( 51) = 2.05
      PaulingX( 52) = 2.10
      PaulingX( 53) = 2.66
      PaulingX( 54) = 2.60
      PaulingX( 55) = 0.79
      PaulingX( 56) = 0.89
      PaulingX( 57) = 1.10
      PaulingX( 58) = 1.12
      PaulingX( 59) = 1.13
      PaulingX( 60) = 1.14
      PaulingX( 61) = 1.13
      PaulingX( 62) = 1.17
      PaulingX( 63) = 1.20
      PaulingX( 64) = 1.20
      PaulingX( 65) = 1.10
      PaulingX( 66) = 1.22
      PaulingX( 67) = 1.23
      PaulingX( 68) = 1.24
      PaulingX( 69) = 1.25
      PaulingX( 70) = 1.10
      PaulingX( 71) = 1.27
      PaulingX( 72) = 1.30
      PaulingX( 73) = 1.50
      PaulingX( 74) = 2.36
      PaulingX( 75) = 1.90
      PaulingX( 76) = 2.20
      PaulingX( 77) = 2.20
      PaulingX( 78) = 2.28
      PaulingX( 79) = 2.54
      PaulingX( 80) = 2.00
      PaulingX( 81) = 1.62
      PaulingX( 82) = 1.87
      PaulingX( 83) = 2.02
      PaulingX( 84) = 2.00
      PaulingX( 85) = 2.20
      PaulingX( 86) = 2.20
      PaulingX( 87) = 0.70
      PaulingX( 88) = 0.90
      PaulingX( 89) = 1.10
      PaulingX( 90) = 1.30
      PaulingX( 91) = 1.50
      PaulingX( 92) = 1.38
      PaulingX( 93) = 1.36
      PaulingX( 94) = 1.28
      PaulingX( 95) = 1.13
      PaulingX( 96) = 1.28
      PaulingX( 97) = 1.30
      PaulingX( 98) = 1.30
      PaulingX( 99) = 1.30
      PaulingX(100) = 1.30
      PaulingX(101) = 1.30
      PaulingX(102) = 1.30
      PaulingX(103) = 1.30
      PaulingX(104) = 0.00
      PaulingX(105) = 0.00
      PaulingX(106) = 0.00
      PaulingX(107) = 0.00
      PaulingX(108) = 0.00

    End Subroutine Set_Pauling_Electronegativity



   End Module CFML_BVSpar
