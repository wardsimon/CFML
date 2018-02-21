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
!!---- MODULE: CFML_Keywords_Code_Parser
!!----   INFO: Refinable Codes for Parameters
!!----   The procedures gathered in this module have as purpose the control of refinable
!!----   parameters in optimization techniques (least squares, simulated annealing, etc).
!!----   Many derived types in CrysFML have as component the code of a particular parameter
!!----   indicating that it may be refined in an optimization procedure with an associated
!!----   multiplier. In this module we associate a name to the particular parameter and we
!!----   add the parameter value, the code, the multiplier, the name, the range, step and
!!----   boundary conditions to general arrays of this module called V_vec, V_vec_std, V_name,
!!----   V_list, V_bounds and V_shift. The code of each derived type points to a particular
!!----   element of the V-arrays. The last refined parameter corresponds to the global variable
!!----   NP_refi that is updated each time a VARY or GVARY directive is processed.
!!----   There are procedures for deleting refined parameters by processing FIX or GFIX
!!----   instruction in the input CFL file.
!!----   For non-atomic parameters there are no specific derived types in CrysFML. In this
!!----   module we provide a simple derived type called "NonAtomic_Parameter_Type" containing
!!----   two real values, an integer code, a multiplier and a name. We provide also a derived
!!----   type "NonAtomic_Parameter_List_Type" containing an allocatable array of objects with
!!----   "NonAtomic_Parameter_Type", "par" and an integer "npar" with the effective number
!!----   of allocated elements. It is responsibility of the user of this module to allocate
!!----   an fill up an object of derived type "NonAtomic_Parameter_List_Type", before use of
!!----   the procedures of the present module. A call to Allocate_VParam with the maximum number
!!----   of expected refinable parameters is needed before the interpretation of VARY/GVARY (or
!!----   FIX/GFIX) instructions if the input CFL file(s)
!!----
!!----
!!---- HISTORY
!!----    Update: 07/03/2011
!!----    Modifications: 12/2011
!!----    Modifications: 02/2012
!!----    Modifications: 10/2013 (Introduction of non-atomic parameters, correction of bugs,
!!----                            and merging of some procedures, JRC)
!!----
!!----
!!---- DEPENDENCIES
!!----
!!----
!!---- VARIABLES
!!--..    Types
!!----    ANGLE_RESTRAINT_TYPE
!!----    DISTANCE_RESTRAINT_TYPE
!!----    NONATOMIC_PARAMETER_TYPE
!!----    NONATOMIC_PARAMETER_LIST_TYPE
!!----    TORSION_RESTRAINT_TYPE
!!--..
!!----    ANG_REST
!!----    CODE_NAM
!!----    mCODE_NAM                    !NEW 12/11
!!----    DIS_REST
!!----    ERR_REFCODES
!!----    ERR_REFCODES_MESS
!!--++    NCODE                        [Private]
!!--++    NKEY                         [Private]
!!--++    mNCODE                       [Private] !NEW 12/11
!!--++    mNKEY                        [Private] !NEW 12/11
!!----    KEY_CODE
!!----    KEY_mCODE                    !NEW 12/11
!!----    NP_CONS
!!----    NP_MAX
!!----    NP_REFI
!!----    NP_REST_ANG
!!----    NP_REST_DIS
!!----    NP_REST_TOR
!!----    TORSION_REST
!!----    V_BCON
!!----    V_BOUNDS
!!----    V_LIST
!!----    V_NAME
!!----    V_VEC
!!----    V_VEC_STD                    !NEW 11/13 (JRC)
!!----    V_SHIFT
!!----
!!---- PUBLIC PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ALLOCATE_RESTPARAM
!!----       ALLOCATE_VPARAM
!!----       DELETE_ELEMENT_IN_VARRAYS   [Private]    !NEW 11/13
!!--++       DELETE_REFCODES             [Private]
!!--++       DELETE_REFCODES_FATOM       [Overloaded]
!!--++       DELETE_REFCODES_FmATOM      [Overloaded] !NEW 12/11
!!--++       DELETE_REFCODES_MOLCRYS     [Overloaded]
!!--++       DELETE_REFCODES_MOLEC       [Overloaded]
!!--++       DELETE_REFCODES_MAGDOM                   !NEW 02/12
!!--++
!!--++       FILL_REFCODES_FATOM         [Overloaded]
!!--++       FILL_REFCODES_FmATOM        [Overloaded] !NEW 12/11
!!--++       FILL_REFCODES_MOLCRYS       [Overloaded]
!!--++       FILL_REFCODES_MOLEC         [Overloaded]
!!--++       FILL_REFCODES_MAGDOM                     !NEW 02/12
!!--++       FILL_REFGCODES              [Private]    !NEW 11/13
!!--++       GET_ATOMBET_CTR             [Private]
!!--++       GET_ATOMPOS_CTR             [Private]
!!--++       GET_CONCODES_LINE           [Private]
!!--++       GET_CONCODES_LINE_FATOM     [Overloaded]
!!--++       GET_CONCODES_LINE_FmATOM    [Overloaded] !NEW 12/11
!!--++       GET_CONCODES_LINE_MOLCRYS   [Overloaded]
!!--++       GET_CONCODES_LINE_MOLEC     [Overloaded]
!!--++       GET_CONCODES_LINE_MAGDOM                 !NEW 02/12
!!--++       GET_REFCODES_LINE           [Private]
!!--++       GET_REFCODES_LINE_FATOM     [Overloaded]
!!--++       GET_REFCODES_LINE_FmATOM    [Overloaded] !NEW 12/11
!!--++       GET_REFCODES_LINE_MOLCRYS   [Overloaded]
!!--++       GET_REFCODES_LINE_MOLEC     [Overloaded]
!!--++       GET_REFCODES_LINE_MAGDOM                 !NEW 02/12
!!--++       GET_REFGCODES_LINE          [Private]    !NEW 11/13
!!----       GET_RESTANG_LINE
!!----       GET_RESTDIS_LINE
!!----       GET_RESTTOR_LINE
!!----       INIT_ERR_REFCODES
!!----       INIT_REFCODES                            !NEW 11/13 (merged with optional arguments, JRC)
!!----       READ_REFCODES_FILE
!!--++       READ_REFCODES_FILE_FATOM    [Overloaded]
!!--++       READ_REFCODES_FILE_FmATOM   [Overloaded] !NEW 12/11 !replaced 02/12
!!--++       READ_REFCODES_FILE_MagStr                !NEW 02/12
!!--++       READ_REFCODES_FILE_MOLCRYS  [Overloaded]
!!--++       READ_REFCODES_FILE_MOLEC    [Overloaded]
!!----       READ_REFGCODES_FILE                      !NEW 11/13 (JRC, non atomic parameters)
!!--++       SPLIT_OPERATIONS            [Private]
!!--++       SPLIT_mOPERATIONS           [Private]    !NEW 12/11
!!----       VSTATE_TO_ATOMSPAR
!!--++       VSTATE_TO_ATOMSPAR_FATOM    [Overloaded]
!!--++       VSTATE_TO_ATOMSPAR_FmATOM   [Overloaded] !NEW 12/11 !modified 02/12 to include MagDom
!!--++       VSTATE_TO_ATOMSPAR_MOLCRYS  [Overloaded]
!!--++       VSTATE_TO_ATOMSPAR_MOLEC    [Overloaded]
!!----       VSTATE_TO_MODELPAR                       !NEW 11/13
!!----       WRITE_INFO_REFCODES
!!--++       WRITE_INFO_REFCODES_FATOM   [Overloaded]
!!--++       WRITE_INFO_REFCODES_FmATOM  [Overloaded] !NEW 12/11 !replaced 02/12
!!--++       WRITE_INFO_REFCODES_Magstr               !NEW 02/12
!!--++       WRITE_INFO_REFCODES_MOLCRYS [Overloaded]
!!--++       WRITE_INFO_REFCODES_MOLEC   [Overloaded]
!!----       WRITE_INFO_REFGCODES                     !NEW 11/13
!!----       WRITE_INFO_REFPARAMS
!!----       WRITE_RESTRAINTS_OBSCALC
!!----
!!
 Module CFML_Keywords_Code_Parser
    !---- Modules ----!
    Use CFML_GlobalDeps,                only: cp
    Use CFML_Math_General,              only: Sort
    Use CFML_String_Utilities,          only: Cutst, U_Case, L_Case, Getword, GetNum
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Get_Stabilizer, Symmetry_Symbol,   &
                                              Sym_B_Relations, Read_SymTrans_Code, Get_SymSymb,Read_Xsym
    Use CFML_Atom_TypeDef,              only: Atom_list_Type, mAtom_list_Type
    Use CFML_Molecular_Crystals,        only: Molecule_Type, Molecular_Crystal_Type
    Use CFML_IO_Formats,                only: File_List_Type
    Use CFML_Magnetic_Symmetry
    Use CFML_Math_3D,                   only: Get_Spheric_Coord

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public :: Allocate_VParam, Init_RefCodes, Read_RefCodes_File, VState_to_AtomsPar,  &
              Write_Info_RefCodes, Get_RestAng_Line, Get_RestDis_Line, Get_RestTor_Line, &
              Allocate_RestParam, Write_Restraints_ObsCalc, Init_Err_RefCodes, &
              Write_Info_RefParams, Read_RefGCodes_File, Write_Info_RefGCodes, &
              VState_to_ModelPar

    !---- List of private functions ----!

    !---- List of private subroutines ----!
    private :: Delete_RefCodes,   &
               Delete_RefCodes_FAtom, Delete_RefCodes_FmAtom, Delete_RefCodes_Molcrys,          &
               Delete_RefCodes_Molec, Delete_RefCodes_Magdom, &
               Fill_RefCodes,Fill_RefGCodes,     &
               Fill_RefCodes_FAtom, Fill_RefCodes_FmAtom, Fill_RefCodes_Molcrys,                &
               Fill_RefCodes_Molec, Fill_RefCodes_Magdom,  &
               Get_AtomBet_Ctr, Get_Atompos_Ctr,                                                &
               Get_ConCodes_Line, &
               Get_ConCodes_Line_FAtom, Get_ConCodes_Line_FmAtom, Get_ConCodes_Line_Molcrys,    &
               Get_ConCodes_Line_Molec, Get_ConCodes_Line_Magdom, &
               Get_RefCodes_Line,Get_RefGCodes_Line, &
               Get_RefCodes_Line_FAtom, Get_RefCodes_Line_FmAtom, Get_RefCodes_Line_Molcrys,    &
               Get_RefCodes_Line_Molec, Get_RefCodes_Line_Magdom, &
               Read_RefCodes_File_FAtom, Read_RefCodes_File_MagStr, Read_RefCodes_File_Molcrys, &
               Read_RefCodes_File_Molec, &
               Split_Operations, Split_mOperations, &
               VState_to_AtomsPar_FAtom, VState_to_AtomsPar_FmAtom, VState_to_AtomsPar_Molcrys, &
               VState_to_AtomsPar_Molec, &
               Write_Info_RefCodes_FAtom,Write_Info_RefCodes_MagStr,Write_Info_RefCodes_Molcrys,&
               Write_Info_RefCodes_Molec

    !---- Definitions ----!

    !!----
    !!---- TYPE :: ANGLE_RESTRAINT_TYPE
    !!--..
    !!---- Type, public :: Angle_Restraint_Type
    !!----    real(kind=cp)                 :: AObs
    !!----    real(kind=cp)                 :: ACalc
    !!----    real(kind=cp)                 :: Sigma
    !!----    integer,dimension(3)          :: P
    !!----    character(len=8),dimension(2) :: STCode
    !!---- End Type Angle_Restraint_Type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: Angle_Restraint_Type
      real(kind=cp)                 :: AObs
      real(kind=cp)                 :: ACalc
      real(kind=cp)                 :: Sigma
      integer,dimension(3)          :: P
      character(len=8),dimension(2) :: STCode
    End Type Angle_Restraint_Type

    !!----
    !!---- TYPE :: DISTANCE_RESTRAINT_TYPE
    !!--..
    !!---- Type, public :: Distance_Restraint_Type
    !!----    real(kind=cp)        :: DObs
    !!----    real(kind=cp)        :: DCalc
    !!----    real(kind=cp)        :: Sigma
    !!----    integer,dimension(2) :: P
    !!----    character(len=8)     :: STCode    ! _N.ABC
    !!---- End Type Distance_Restraint_Type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: Distance_Restraint_Type
      real(kind=cp)        :: DObs
      real(kind=cp)        :: DCalc
      real(kind=cp)        :: Sigma
      integer,dimension(2) :: P
      character(len=8)     :: STCode
    End Type Distance_Restraint_Type

    !!----
    !!---- TYPE :: NONATOMIC_PARAMETER_TYPE
    !!--..
    !!---- Type, public :: Nonatomic_Parameter_Type
    !!----    real(kind=cp)        :: Value
    !!----    real(kind=cp)        :: Sigma
    !!----    integer              :: Lcode
    !!----    real(kind=cp)        :: multip
    !!----    character(len=20)    :: Nam
    !!---- End Type Nonatomic_Parameter_Type
    !!----
    !!---- Update: November 1 - 2013
    !!
    Type, public :: Nonatomic_Parameter_Type
       real(kind=cp)        :: Value
       real(kind=cp)        :: Sigma
       integer              :: Lcode
       real(kind=cp)        :: multip
       character(len=20)    :: Nam
    End Type Nonatomic_Parameter_Type

    !!----
    !!---- TYPE :: NONATOMIC_PARAMETER_LIST_TYPE
    !!--..
    !!---- Type, public :: Nonatomic_Parameter_List_Type
    !!----    Integer                                                       :: npar
    !!----    Type(Nonatomic_Parameter_Type),dimension(:),allocatable :: par
    !!---- End Type Nonatomic_Parameter_List_Type
    !!----
    !!---- The user of this derived type must allocate an fill a type of this kind for his(her) own
    !!---- his(her) own problem. An object of this type is needed to fill the V_vec and
    !!---- accompanying arrays.
    !!----
    !!---- Update: November 1 - 2013
    !!
    Type, public :: Nonatomic_Parameter_List_Type
       Integer                                                 :: npar
       Type(Nonatomic_Parameter_Type),dimension(:),allocatable :: par
    End Type Nonatomic_Parameter_List_Type

    !!----
    !!---- TYPE :: TORSION_RESTRAINT_TYPE
    !!--..
    !!---- Type, public :: Torsion_Restraint_Type
    !!----    real(kind=cp)                 :: TObs
    !!----    real(kind=cp)                 :: TCalc
    !!----    real(kind=cp)                 :: Sigma
    !!----    integer,dimension(4)          :: P
    !!----    character(len=8),dimension(3) :: STCode
    !!---- End Type Torsion_Restraint_Type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: Torsion_Restraint_Type
      real(kind=cp)                 :: TObs
      real(kind=cp)                 :: TCalc
      real(kind=cp)                 :: Sigma
      integer,dimension(4)          :: P
      character(len=8),dimension(3) :: STCode
    End Type Torsion_Restraint_Type

    !!----
    !!---- ANG_REST
    !!----    type(Angle_Restraint_Type), public, dimension(:), allocatable :: Ang_rest
    !!----
    !!---- Relations for Angle Restraints
    !!----
    !!---- Update: March - 2005
    !!
    type(Angle_Restraint_Type), public, dimension(:), allocatable :: Ang_Rest

    !!--++
    !!--++ NCODE
    !!--++    integer, private :: NCode
    !!--++
    !!--++    Number of Code variables
    !!--++
    !!--++ Update: March - 2005
    !!
    integer, private, parameter :: NCode=21

    !!----
    !!---- CODE_NAM
    !!----    character(len=*), dimension(NCode), public, parameter :: Code_Nam
    !!----
    !!----    Variable for treatement Codes
    !!--..    21 debe ser igual a NCode
    !!----
    !!---- Update: March - 2005
    !!
    character(len=*), dimension(NCode), public, parameter :: Code_Nam=(/ "X_    ","Y_    ","Z_    ",       &
                                                                         "B_    ","Occ_  ","B11_  ",       &
                                                                         "B22_  ","B33_  ","B12_  ",       &
                                                                         "B13_  ","B23_  ","Bns_  ",       &
                                                                         "Xc_   ","Yc_   ","Zc_   ",       &
                                                                         "Theta_","Phi_  ","Chi_  ",       &
                                                                         "Th_L_ ","Th_T_ ","Th_S_ "/)
    !!--++
    !!--++ NCODE
    !!--++    integer, private :: NGCode
    !!--++
    !!--++    Number of Codes for non atomic variables
    !!--++
    !!--++ Update: November - 2013
    !!
    integer, private, parameter :: NGCode=87

    !!----
    !!---- GCODE_NAM
    !!----    character(len=*), dimension(NGCode), public, parameter :: GCode_Nam
    !!----
    !!----    Variable for treatement of GCodes
    !!--..
    !!----
    !!---- Update: November - 2013, February 2014
    !!
    character(len=*), dimension(NGCode), public, parameter :: &
                      GCode_Nam=(/"Scalef ","Cell-a ","Cell-b ", &
                                  "Cell-c ","C-alpha","C-beta ", &
                                  "C-gamma","Cell   ","Up     ","Vp     ", &
                                  "Wp     ","Xp     ","Yp     ","Xfract ", &
                                  "Size   ","Gsize  ","Strain ", &
                                  "LStrain","Qbroad ","Qdamp  ", &
                                  "delta1 ","delta2 ","Lratio ","Sharpf ", &
                                  "Sratio ","Zero   ","Diff1  ", &
                                  "Diff2  ","eta    ","sig2   ", &
                                  "sig1   ","sig0   ","gamm2  ", &
                                  "gamm1  ","gamm0  ","alph0  ", &
                                  "alph1  ","beta0  ","beta1  ", &
                                  "kappa  ","Bover  ","Abs1   ", &
                                  "Abs2   ","Extinct","BCExt1 ", &
                                  "BCExt2 ","Ext11  ","Ext22  ", &
                                  "Ext33  ","Ext12  ","Ext13  ", &
                                  "Ext23  ","S400   ","S040   ", &
                                  "S004   ","S220   ","S202   ", &
                                  "S022   ","S211   ","S121   ", &
                                  "S112   ","S310   ","S301   ", &
                                  "S130   ","S103   ","S013   ", &
                                  "S031   ","Sizh2  ","Sizk2  ", &
                                  "Sizl2  ","Siz2hk ","Siz2hl ", &
                                  "Siz2kl ","kx     ","ky     ", &
                                  "kz     ","bkg    ","bkg_   ", &
                                  "kx_    ","ky_    ","kz_    ", &
                                  "Sc_    ","sycos  ","sysin  ", &
                                  "Dom_   ","sycos_ ","sysin_ "/)
    integer, private, parameter :: mNCode=25

    !!----
    !!---- mCODE_NAM
    !!----    character(len=*), dimension(mNCode), public, parameter :: mCode_Nam
    !!----
    !!----    Variable for treatement Codes
    !!----    mag clone of CODE_NAM
    !!---- Created: December - 2011
    !!
    character(len=*),dimension(mNCode), public, parameter :: mCode_Nam=(/"Rx_   ","Ry_   ","Rz_   ",&
                                                                         "Ix_   ","Iy_   ","Iz_   ",&
                                                                         "Rm_   ","Rphi_ ","Rth_  ",&
                                                                         "Im_   ","Iphi_ ","Ith_  ",&
                                                                         "MagPh_",                  &
                                                                         "C1_   ","C2_   ","C3_   ",&
                                                                         "C4_   ","C5_   ","C6_   ",&
                                                                         "C7_   ","C8_   ","C9_   ",&
                                                                         "C10_  ","C11_  ","C12_  "/)

    !!----
    !!---- DIS_REST
    !!----    type(Distance_Restraint_Type), public, dimension(:), allocatable :: Dis_Rest
    !!----
    !!---- Relations for Angle Restraints
    !!----
    !!---- Update: March - 2005
    !!
    type(Distance_Restraint_Type), public, dimension(:), allocatable :: Dis_Rest

    !!----
    !!---- ERR_REFCODES
    !!----    logical, public :: Err_RefCodes
    !!----
    !!----    Error variable
    !!----
    !!---- Update: March - 2005
    !!
    logical, public :: Err_RefCodes = .false.

    !!----
    !!---- ERR_REFCODES_MESS
    !!----    character(len=150), public :: ERR_RefCodes_Mess
    !!----
    !!----    Error variable messages
    !!----
    !!---- Update: March - 2005
    !!
    character(len=150), public :: ERR_RefCodes_Mess = " "


    !!--++
    !!--++ NKEY
    !!--++    integer, public :: NKey
    !!--++
    !!--++    Number of Keys variables
    !!--++
    !!--++ Update: March - 2005
    !!
    integer, private, parameter :: NKey=8

    !!----
    !!---- KEY_CODE
    !!----    character(len=*), dimension(NKey), public, parameter :: key_Code
    !!----
    !!----     Key codes defined in the module
    !!----
    !!---- Update: March - 2005
    !!
    character(len=*), dimension(Nkey), public, parameter :: Key_Code=(/ "XYZ ","OCC ","BIS ","BAN ", &
                                                                        "ALL ","CEN ","ORI ","THE "/)
    integer, private, parameter :: mNKey=4

    !!----
    !!---- KEY_mCODE
    !!----    character(len=*), dimension(mNKey), public, parameter :: key_mCode
    !!----
    !!----    Key codes defined in the module
    !!----    mag clone of KEY_CODE
    !!---- Created: December - 2011
    !!
    character(len=*), dimension(mNkey), public, parameter :: Key_mCode=(/"Rxyz","Ixyz","Mxyz","Magd"/)

    !!----
    !!---- NP_CONS
    !!----    integer, public :: NP_Cons
    !!----
    !!----    Number of Constraints relations
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Cons

    !!----
    !!---- NP_MAX
    !!----    integer, public :: NP_Max
    !!----
    !!----    Number of Maximum Parameters to Refine
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Max

    !!----
    !!---- NP_REFI
    !!----    integer, public :: NP_Refi
    !!----
    !!----    Number of Refinable Parameters
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Refi

    !!----
    !!---- NP_REST_ANG
    !!----    integer, save, public :: NP_Rest_Ang=0
    !!----
    !!----    Number of Angle Restraints relations
    !!----
    !!---- Update: March-2005, May-2015
    !!
    integer, save, public :: NP_Rest_Ang=0

    !!----
    !!---- NP_REST_DIS
    !!----    integer, save, public :: NP_Rest_Dis=0
    !!----
    !!----    Number of Distance Restraints relations
    !!----
    !!---- Update: March-2005, May-2015
    !!
    integer, save, public :: NP_Rest_Dis=0

    !!----
    !!---- NP_REST_TOR
    !!----    iinteger, save, public :: NP_Rest_Tor=0
    !!----
    !!----    Number of Torsion Restraints relations
    !!----
    !!---- Update: March - 2005, May-2015
    !!
    integer, save, public :: NP_Rest_Tor=0

    !!----
    !!---- TOR_REST
    !!----    type(Torsion_Restraint_Type), public, dimension(:), allocatable :: Tor_Rest
    !!----
    !!---- Relations for Torsion Angle Restraints
    !!----
    !!---- Update: March - 2005
    !!
    type(Torsion_Restraint_Type), public, dimension(:), allocatable :: Tor_Rest

    !!----
    !!---- V_BCON
    !!----    integer, public, dimension(:), allocatable :: V_BCon
    !!----
    !!----    Vector of Boundary Conditions
    !!----
    !!---- Update: March - 2005
    !!
    integer, public, dimension(:), allocatable :: V_BCon

    !!----
    !!---- V_BOUNDS
    !!----    real(kind=cp), public, dimension(:,:), allocatable :: V_Bounds
    !!----
    !!----    Vector of Lower, Upper limits and Step for Parameters
    !!----
    !!---- Update: March - 2005
    !!
    real(kind=cp), public, dimension(:,:),  allocatable :: V_Bounds

    !!----
    !!---- V_LIST
    !!----    integer, public, dimension(:), allocatable :: V_List
    !!----
    !!----    Vector of indices pointing to the parameter number
    !!----
    !!---- Update: March - 2005
    !!
    integer, public, dimension(:),  allocatable :: V_List

    !!----
    !!---- V_NAME
    !!----    character(len=20), public, dimension(:), allocatable :: V_Name
    !!----
    !!----    Vector of  Name of Refinable Parameters
    !!----
    !!---- Update: March - 2005
    !!
    character(len=20), public, dimension(:), allocatable :: V_Name

    !!----
    !!---- V_VEC
    !!----    real(kind=cp), public, dimension(:), allocatable :: V_Vec
    !!----
    !!----    Vector of Parameters
    !!----
    !!---- Update: March - 2005
    !!

    real(kind=cp), public, dimension(:),    allocatable :: V_Vec
    !!----
    !!---- V_SAVE
    !!----    real(kind=cp), public, dimension(:), allocatable :: V_Save
    !!----
    !!----    Vector of Parameters saving previous values of parameters.
    !!----    This vector is handled by the calling program. It is only automatically
    !!----    allocated in this module by a call to Allocate_VParam.
    !!----
    !!---- Update: April - 2014
    !!

    real(kind=cp), public, dimension(:),    allocatable :: V_Save
    !!----
    !!---- V_VEC_std
    !!----    real(kind=cp), public, dimension(:), allocatable :: V_Vec_std
    !!----
    !!----    Standard deviations of the parameters
    !!----
    !!---- Update: November - 2013 (JRC)
    !!
    real(kind=cp), public, dimension(:),    allocatable :: V_Vec_std

    !!----
    !!---- V_SHIFT
    !!----    real(kind=cp), public, dimension(:), allocatable :: V_Shift
    !!----
    !!----    Vector of holding the shift of parameters
    !!----
    !!---- Update: March - 2005
    !!
    real(kind=cp), public, dimension(:),    allocatable :: V_Shift

    !---- Interfaces - Overloaded ----!
    Interface Allocate_RestParam
       Module Procedure Allocate_RestParam_Single
       Module Procedure Allocate_RestParam_Mult
    End Interface

    Interface Delete_RefCodes
       Module Procedure Delete_RefCodes_FAtom
       Module Procedure Delete_RefCodes_FmAtom
       Module Procedure Delete_RefCodes_Molcrys
       Module Procedure Delete_RefCodes_Molec
       Module Procedure Delete_RefCodes_Magdom
    End Interface

    Interface Fill_RefCodes
       Module Procedure Fill_RefCodes_FAtom
       Module Procedure Fill_RefCodes_FmAtom
       Module Procedure Fill_RefCodes_Molcrys
       Module Procedure Fill_RefCodes_Molec
       Module Procedure Fill_RefCodes_Magdom
    End Interface

    Interface Get_ConCodes_Line
       Module Procedure Get_ConCodes_Line_FAtom
       Module Procedure Get_ConCodes_Line_FmAtom
       Module Procedure Get_ConCodes_Line_Molcrys
       Module Procedure Get_ConCodes_Line_Molec
       Module Procedure Get_ConCodes_Line_Magdom
    End Interface

    Interface Get_RefCodes_Line
       Module Procedure Get_RefCodes_Line_FAtom
       Module Procedure Get_RefCodes_Line_FmAtom
       Module Procedure Get_RefCodes_Line_Molcrys
       Module Procedure Get_RefCodes_Line_Molec
       Module Procedure Get_RefCodes_Line_Magdom
    End Interface

    Interface Read_RefCodes_File
       Module Procedure Read_RefCodes_File_FAtom
       Module Procedure Read_RefCodes_File_MagStr
       Module Procedure Read_RefCodes_File_Molcrys
       Module Procedure Read_RefCodes_File_Molec
    End Interface

    Interface VState_to_AtomsPar
       Module Procedure VState_to_AtomsPar_FAtom
       Module Procedure VState_to_AtomsPar_FmAtom
       Module Procedure VState_to_AtomsPar_Molcrys
       Module Procedure VState_to_AtomsPar_Molec
    End Interface

    Interface Write_Info_RefCodes
       Module Procedure Write_Info_RefCodes_FAtom
       Module Procedure Write_Info_RefCodes_MagStr
       Module Procedure Write_Info_RefCodes_Molcrys
       Module Procedure Write_Info_RefCodes_Molec
    End Interface

 Contains

    !!--++
    !!--++ Subroutine Allocate_RestParam_Single(file_dat)
    !!--++    Type(file_list_type),     intent( in)    :: file_dat
    !!--++
    !!--++    Allocate vectors Ang_Rest, Dist_Rest, Tor_Rest. It is supposed
    !!--++    that a single phase is at work.
    !!--++
    !!--++ Update: March-2005, May-2015
    !!
    Subroutine Allocate_RestParam_Single(file_dat)
       !---- Arguments ----!
       Type(file_list_type),     intent( in) :: file_dat

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=15),dimension(40) :: car
       integer                         :: i,nc,nr

       if (allocated(Ang_Rest)) deallocate(Ang_Rest)
       if (allocated(Dis_Rest)) deallocate(Dis_Rest)
       if (allocated(Tor_Rest)) deallocate(Tor_Rest)

       NP_Rest_Ang=0
       NP_Rest_Dis=0
       NP_Rest_Tor=0

       !---- Dimension for AFIX ----!
       nr=0
       do i=1,file_dat%nlines
          line=adjustl(file_dat%line(i))
          if (u_case(line(1:4)) /= "AFIX") cycle
          call cutst(line)
          call getword(line,car,nc)
          nr=nr+nc/3
       end do
       if (nr >0) then
          allocate(Ang_Rest(nr))
          ang_rest%aobs =0.0
          ang_rest%acalc=0.0
          ang_rest%sigma=0.0
          ang_rest%p(1) = 0
          ang_rest%p(2) = 0
          ang_rest%STCode(1)=" "
          ang_rest%STCode(2)=" "
       end if

       !---- Dimension for DFIX ----!
       nr=0
       do i=1,file_dat%nlines
          line=adjustl(file_dat%line(i))
          if (u_case(line(1:4)) /= "DFIX") cycle
          call cutst(line)
          call getword(line,car,nc)
          nr=nr+nc/2
          if (modulo(nc,2) == 0) nr=nr-1
       end do
       if (nr > 0) then
          allocate(Dis_Rest(nr))
          dis_rest%dobs =0.0
          dis_rest%dcalc=0.0
          dis_rest%sigma=0.0
          dis_rest%p(1) = 0
          dis_rest%p(2) = 0
          dis_rest%STCode=" "
       end if

       !---- Dimension for TFIX ----!
       nr=0
       do i=1,file_dat%nlines
          line=adjustl(file_dat%line(i))
          if (u_case(line(1:4)) /= "TFIX") cycle
          call cutst(line)
          call getword(line,car,nc)
          nr=nr+nc/4
       end do
       if (nr > 0) then
          allocate(Tor_Rest(nr))
          tor_rest%tobs=0.0
          tor_rest%tcalc=0.0
          tor_rest%sigma=0.0
          tor_rest%p(1) = 0
          tor_rest%p(2) = 0
          tor_rest%STCode(1)=" "
          tor_rest%STCode(2)=" "
          tor_rest%STCode(3)=" "
       end if

       return
    End Subroutine Allocate_RestParam_Single

    !!--++
    !!--++ Subroutine Allocate_RestParam_Mult(file_dat,ndis_rest,nang_rest,ntor_rest)
    !!--++   Type(file_list_type),dimension(:),  intent( in) :: file_dat
    !!==++   integer,              dimension(:), intent(out) :: ndis_rest
    !!==++   integer,              dimension(:), intent(out) :: nang_rest
    !!==++   integer,              dimension(:), intent(out) :: ntor_rest
    !!--++
    !!--++    Allocate vectors Ang_Rest, Dist_Rest, Tor_Rest.
    !!--++    This is for a context of multiple phases. All restrainst
    !!--++    are stored in single array Ang_Rest, Dist_Rest, Tor_Rest, however
    !!--++    the information about the phases is stored in xxx_rest integer
    !!--++    arrays, from which one can know to which phase belong a particular
    !!--++    restrain. This works only if there is a single line per restraint
    !!--++    in the restraint files. The arrays xxx_rest stores the number of
    !!--++    restraints of each type in each phase.
    !!--++
    !!--++ Update: March-2005, May-2015
    !!
    Subroutine Allocate_RestParam_Mult(file_dat,ndis_rest,nang_rest,ntor_rest)
       !---- Arguments ----!
       Type(file_list_type), dimension(:), intent( in) :: file_dat
       integer,              dimension(:), intent(out) :: ndis_rest
       integer,              dimension(:), intent(out) :: nang_rest
       integer,              dimension(:), intent(out) :: ntor_rest

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=15),dimension(40) :: car
       integer                         :: i,nc,nr,nfiles,nf

       if (allocated(Ang_Rest)) deallocate(Ang_Rest)
       if (allocated(Dis_Rest)) deallocate(Dis_Rest)
       if (allocated(Tor_Rest)) deallocate(Tor_Rest)

       NP_Rest_Ang=0
       NP_Rest_Dis=0
       NP_Rest_Tor=0
       nfiles=size(ndis_rest)  !number of phases
       ndis_rest=0; nang_rest=0; ntor_rest=0
       !---- Dimension for AFIX ----!
       nr=0
       do nf=1,nfiles
         nang_rest(nf)=0
         do i=1,file_dat(nf)%nlines
            line=adjustl(file_dat(nf)%line(i))
            if (u_case(line(1:4)) /= "AFIX") cycle
            call cutst(line)
            call getword(line,car,nc)
            nr=nr+nc/3
            nang_rest(nf)=nang_rest(nf)+1
         end do
       end do
       if (nr > 0) then
          allocate(Ang_Rest(nr))
          ang_rest%aobs =0.0
          ang_rest%acalc=0.0
          ang_rest%sigma=0.0
          ang_rest%p(1) = 0
          ang_rest%p(2) = 0
          ang_rest%STCode(1)=" "
          ang_rest%STCode(2)=" "
       end if

       !---- Dimension for DFIX ----!
       nr=0
       do nf=1,nfiles
         ndis_rest(nf)=0
         do i=1,file_dat(nf)%nlines
            line=adjustl(file_dat(nf)%line(i))
            if (u_case(line(1:4)) /= "DFIX") cycle
            call cutst(line)
            call getword(line,car,nc)
            nr=nr+nc/2
            if (modulo(nc,2) == 0) nr=nr-1
            ndis_rest(nf)=ndis_rest(nf)+1
         end do
       end do
       if (nr > 0) then
          allocate(Dis_Rest(nr))
          dis_rest%dobs =0.0
          dis_rest%dcalc=0.0
          dis_rest%sigma=0.0
          dis_rest%p(1) = 0
          dis_rest%p(2) = 0
          dis_rest%STCode=" "
       end if

       !---- Dimension for TFIX ----!
       nr=0
       do nf=1,nfiles
         ntor_rest(nf)=0
         do i=1,file_dat(nf)%nlines
            line=adjustl(file_dat(nf)%line(i))
            if (u_case(line(1:4)) /= "TFIX") cycle
            call cutst(line)
            call getword(line,car,nc)
            nr=nr+nc/4
            ntor_rest(nf)=ntor_rest(nf)+1
         end do
       end do
      if (nr > 0) then
          allocate(Tor_Rest(nr))
          tor_rest%tobs=0.0
          tor_rest%tcalc=0.0
          tor_rest%sigma=0.0
          tor_rest%p(1) = 0
          tor_rest%p(2) = 0
          tor_rest%STCode(1)=" "
          tor_rest%STCode(2)=" "
          tor_rest%STCode(3)=" "
       end if

       return
    End Subroutine Allocate_RestParam_Mult
    !!---- Subroutine Allocate_VParam(N)
    !!----    integer, intent(in) :: N
    !!----
    !!----    Allocate vectors V_Vec, V_Bounds, V_Name, V_Bcon, V_Shift, V_list
    !!----    If N is equal zero it deallocates the vectors
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_VParam(N)
       !---- Arguments ----!
       integer, intent(in) :: N

       if (allocated(V_Vec))    deallocate(V_Vec)
       if (allocated(V_Save))   deallocate(V_Save)
       if (allocated(V_Vec_Std))deallocate(V_Vec_Std)
       if (allocated(V_Name))   deallocate(V_Name)
       if (allocated(V_Bounds)) deallocate(V_Bounds)
       if (allocated(V_BCon))   deallocate(V_BCon)
       if (allocated(V_Shift))  deallocate(V_Shift)
       if (allocated(V_List))   deallocate(V_List)

       if (N > 0) then
          allocate(V_Vec(n))
          V_Vec=0.0
          allocate(V_Save(n))
          V_Save=0.0
          allocate(V_Vec_Std(n))
          V_Vec_Std=0.0
          allocate(V_Name(n))
          V_Name=" "
          allocate(V_Bounds(3,n))
          V_Bounds=0.0
          allocate(V_BCon(n))
          V_BCon=0
          allocate(V_Shift(n))
          V_Shift=0.0
          allocate(V_List(n))
          V_list=0
          np_max=n
       else
          np_max=0
       end if

       return
    End Subroutine Allocate_VParam

    !!--++
    !!--++ Subroutine Delete_RefCodes(N, FAtom/FmAtom/MolCrys/Molec/Mag_Dom)
    !!--++    integer,                      intent(in)     :: N
    !!--++    type(Atom_List_Type),         intent(in out) :: FAtom
    !!--++    or
    !!--++    type(mAtom_List_Type),        intent(in out) :: FmAtom
    !!--++    or
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++    or
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    or
    !!--++    type(Magnetic_Domain_type),intent(in out) :: Mag_Dom
    !!--++
    !!--++    (Private)
    !!--++    Delete the number of Refinable Parameter (N) on the list
    !!--++
    !!--++ Update: March - 2005
    !!

    !!--++
    !!--++ Subroutine Delete_RefCodes_FAtom(N, FAtom)
    !!--++    integer,              intent(in)     :: N
    !!--++    type(Atom_List_Type), intent(in out) :: FAtom
    !!--++
    !!--++ Overloaded
    !!--++ Delete the number of Refinable Parameter (N) on the list
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Delete_RefCodes_FAtom(N, FAtom)
       !---- Arguments ----!
       integer,              intent(in)     :: N
       type(Atom_List_Type), intent(in out) :: FAtom

       !---- Local Variables ----!
       logical :: deleted
       integer :: i,j

       deleted=.false.

       !---- Eliminate N Parameter ----!
       do i=1,FAtom%natoms
          do j=1,3
             if (FAtom%atom(i)%lx(j) == N) then
                 FAtom%atom(i)%lx(j)=0
                 FAtom%atom(i)%mx(j)=0.0
                 deleted=.true.
             end if
          end do

          if (FAtom%atom(i)%lbiso == N) then
              FAtom%atom(i)%lbiso=0
              FAtom%atom(i)%mbiso=0.0
              deleted=.true.
          end if

          if (FAtom%atom(i)%locc == N) then
              FAtom%atom(i)%locc=0
              FAtom%atom(i)%mocc=0.0
              deleted=.true.
          end if

          do j=1,6
             if (FAtom%atom(i)%lu(j) == N) then
                 FAtom%atom(i)%lu(j)=0
                 FAtom%atom(i)%mu(j)=0.0
                 deleted=.true.
             end if
          end do
       end do

       !---- Updating Variables ----!
       do i=1,FAtom%natoms
          do j=1,3
             if (FAtom%atom(i)%lx(j) > N) then
                FAtom%atom(i)%lx(j)=FAtom%atom(i)%lx(j)-1
             end if
          end do

          if (FAtom%atom(i)%lbiso > N) then
             FAtom%atom(i)%lbiso=FAtom%atom(i)%lbiso-1
          end if

          if (FAtom%atom(i)%locc > N) then
             FAtom%atom(i)%locc=FAtom%atom(i)%locc-1
          end if

          do j=1,6
             if (FAtom%atom(i)%lu(j) > N) then
                FAtom%atom(i)%lu(j)=FAtom%atom(i)%lu(j)-1
             end if
          end do
       end do

       !---- Updating V_Vectors ----!
       if (deleted) call Delete_element_in_Varrays(N)

       return
    End Subroutine Delete_RefCodes_FAtom

    !!--++
    !!--++ Subroutine Delete_RefCodes_FmAtom(N, FmAtom)
    !!--++    integer,              intent(in)     :: N
    !!--++    type(mAtom_List_Type),intent(in out) :: FmAtom
    !!--++
    !!--++ Delete the number of Refinable Parameters (N) on the list
    !!--++ magnetic clone of Delete_RefCodes_FAtom
    !!--++ Created: December - 2011
    !!
    Subroutine Delete_RefCodes_FmAtom(N, FmAtom)
       !---- Arguments ----!
       integer,              intent(in)     :: N
       type(mAtom_List_Type), intent(in out) :: FmAtom

       !---- Local Variables ----!
       logical :: deleted
       integer :: i,j,k,ik

       deleted=.false.

       !---- Eliminate N Parameter ----!
       do i=1,FmAtom%natoms
       ik=FmAtom%atom(i)%nvk
          do j=1,3
           do k=1,ik
             if (FmAtom%atom(i)%lSkR(j,k) == N) then
                 FmAtom%atom(i)%mSkR(j,k)=0.0
                 FmAtom%atom(i)%lskr(j,k)=0
                 deleted=.true.
             end if
           end do
          end do

          do j=1,3
           do k=1,ik
             if (FmAtom%atom(i)%lSkI(j,k) == N) then
                 FmAtom%atom(i)%mSkI(j,k)=0.0
                 FmAtom%atom(i)%lski(j,k)=0
                 deleted=.true.
             end if
           end do
          end do

          do k=1,ik
             if (FmAtom%atom(i)%lmphas(k) == N) then
                 FmAtom%atom(i)%mmphas(k)=0.0
                 FmAtom%atom(i)%lmphas(k)=0
                 deleted=.true.
             end if
          end do

          do j=1,12
           do k=1,ik
             if (FmAtom%atom(i)%lbas(j,k) == N) then
                 FmAtom%atom(i)%mbas(j,k)=0.0
                 FmAtom%atom(i)%lbas(j,k)=0
                 deleted=.true.
             end if
           end do
          end do
       end do

       !---- Updating Variables ----!
       do i=1,FmAtom%natoms
          do j=1,3
           do k=1,ik
             if (FmAtom%atom(i)%lSkR(j,k) > N) then
                 FmAtom%atom(i)%lSkR(j,k)=FmAtom%atom(i)%lSkR(j,k)-1
             end if
           end do
          end do

          do j=1,3
           do k=1,ik
             if (FmAtom%atom(i)%lSkI(j,k) > N) then
                 FmAtom%atom(i)%lSkI(j,k)=FmAtom%atom(i)%lSkI(j,k)-1
             end if
           end do
          end do

          do k=1,ik
             if (FmAtom%atom(i)%lmphas(k) > N) then
                 FmAtom%atom(i)%lmphas(k)=FmAtom%atom(i)%lmphas(k)-1
             end if
          end do

          do j=1,12
           do k=1,ik
             if (FmAtom%atom(i)%lbas(j,k) > N) then
                 FmAtom%atom(i)%lbas(j,k)=FmAtom%atom(i)%lbas(j,k)-1
             end if
           end do
          end do
       end do

       !---- Updating V_Vectors ----!
       if (deleted) call Delete_element_in_Varrays(N)

       return
    End Subroutine Delete_RefCodes_FmAtom

    !!--++
    !!--++ Subroutine Delete_RefCodes_MolCrys(N,MolCrys)
    !!--++    integer,                      intent(in)     :: N
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++
    !!--++ Overloaded
    !!--++ Delete the number of Refinable Parameter (N) on the list
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Delete_RefCodes_MolCrys(N,Molcrys)
       !---- Arguments ----!
       integer,                      intent(in)     :: N
       type(molecular_Crystal_type), intent(in out) :: MolCrys

       !---- Local Variables ----!
       logical :: deleted
       integer :: i,j,k

       deleted=.false.

       if (MolCrys%N_Free > 0 ) then
          do i=1,MolCrys%N_Free
             do j=1,3
                if (MolCrys%Atm(i)%lx(j) == N) then
                    MolCrys%Atm(i)%lx(j)=0
                    MolCrys%Atm(i)%mx(j)=0.0
                    deleted=.true.
                end if
             end do

             if (MolCrys%Atm(i)%lbiso == N) then
                 MolCrys%Atm(i)%lbiso=0
                 MolCrys%Atm(i)%mbiso=0.0
                 deleted=.true.
             end if

             if (MolCrys%Atm(i)%locc == N) then
                 MolCrys%Atm(i)%locc=0
                 MolCrys%Atm(i)%mocc=0.0
                 deleted=.true.
             end if

             do j=1,6
                if (MolCrys%Atm(i)%lu(j) == N) then
                    MolCrys%Atm(i)%lu(j)=0
                    MolCrys%Atm(i)%mu(j)=0.0
                    deleted=.true.
                end if
             end do
          end do

          do i=1,MolCrys%N_Free
             do j=1,3
                if (MolCrys%Atm(i)%lx(j) > N) then
                    MolCrys%Atm(i)%lx(j)=MolCrys%Atm(i)%lx(j)-1
                end if
             end do

             if (MolCrys%Atm(i)%lbiso > N) then
                MolCrys%Atm(i)%lbiso=MolCrys%Atm(i)%lbiso-1
             end if

             if (MolCrys%Atm(i)%locc > N) then
                MolCrys%Atm(i)%locc=MolCrys%Atm(i)%locc-1
             end if

             do j=1,6
                if (MolCrys%Atm(i)%lu(j) > N) then
                   MolCrys%Atm(i)%lu(j)=MolCrys%Atm(i)%lu(j)-1
                end if
             end do
          end do
       end if

       if (MolCrys%N_Mol > 0 ) then

          do k=1,MolCrys%N_Mol
             do j=1,3
                if (Molcrys%Mol(k)%lxcentre(j) == N) then
                   Molcrys%Mol(k)%lxcentre(j)=0
                   Molcrys%Mol(k)%mxcentre(j)=0.0
                   deleted=.true.
                end if
             end do

             do j=1,3
                if (Molcrys%Mol(k)%lOrient(j) == N) then
                   Molcrys%Mol(k)%lOrient(j)=0
                   Molcrys%Mol(k)%mOrient(j)=0.0
                   deleted=.true.
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lT_TLS(j) == N) then
                   Molcrys%Mol(k)%lT_TLS(j)=0
                   Molcrys%Mol(k)%mT_TLS(j)=0.0
                   deleted=.true.
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lL_TLS(j) == N) then
                   Molcrys%Mol(k)%lL_TLS(j)=0
                   Molcrys%Mol(k)%mL_TLS(j)=0.0
                   deleted=.true.
                end if
             end do

             do i=1,3
                do j=1,3
                   if (Molcrys%Mol(k)%lS_TLS(i,j) == N) then
                      Molcrys%Mol(k)%lS_TLS(i,j)=0
                      Molcrys%Mol(k)%mS_TLS(i,j)=0.0
                      deleted=.true.
                   end if
                end do
             end do

             !---- Updating ----!
             do j=1,3
                if (Molcrys%Mol(k)%lxcentre(j) > N) then
                   Molcrys%Mol(k)%lxcentre(j)=Molcrys%Mol(k)%lxcentre(j)-1
                end if
             end do

             do j=1,3
                if (Molcrys%Mol(k)%lOrient(j) > N) then
                   Molcrys%Mol(k)%lOrient(j)=Molcrys%Mol(k)%lOrient(j)-1
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lT_TLS(j) > N) then
                   Molcrys%Mol(k)%lT_TLS(j)=Molcrys%Mol(k)%lT_TLS(j)-1
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lL_TLS(j) > N) then
                   Molcrys%Mol(k)%lL_TLS(j)=Molcrys%Mol(k)%lL_TLS(j)-1
                end if
             end do

             do i=1,3
                do j=1,3
                   if (Molcrys%Mol(k)%lS_TLS(i,j) > N) then
                      Molcrys%Mol(k)%lS_TLS(i,j)=Molcrys%Mol(k)%lS_TLS(i,j)-1
                   end if
                end do
             end do

             if (Molcrys%Mol(k)%natoms <=0) cycle

             do i=1,Molcrys%Mol(k)%natoms
                do j=1,3
                   if (MolCrys%Mol(k)%lI_coor(j,i) == N) then
                      MolCrys%Mol(k)%lI_coor(j,i)=0
                      MolCrys%Mol(k)%mI_coor(j,i)=0.0
                      deleted=.true.
                   end if
                end do

                if (MolCrys%Mol(k)%lbiso(i) == N) then
                   MolCrys%Mol(k)%lbiso(i)=0
                   MolCrys%Mol(k)%mbiso(i)=0.0
                   deleted=.true.
                end if

                if (MolCrys%Mol(k)%locc(i) == N) then
                   MolCrys%Mol(k)%locc(i)=0
                   MolCrys%Mol(k)%mocc(i)=0.0
                   deleted=.true.
                end if
             end do

             do i=1,Molcrys%Mol(k)%natoms
                do j=1,3
                   if (MolCrys%Mol(k)%lI_coor(j,i) > N) then
                      MolCrys%Mol(k)%lI_coor(j,i)=MolCrys%Mol(k)%lI_coor(j,i)-1
                   end if
                end do

                if (MolCrys%Mol(k)%lbiso(i) > N) then
                   MolCrys%Mol(k)%lbiso(i)=MolCrys%Mol(k)%lbiso(i)-1
                end if

                if (MolCrys%Mol(k)%locc(i) > N) then
                   MolCrys%Mol(k)%locc(i)=MolCrys%Mol(k)%locc(i)-1
                end if
             end do

          end do
       end if

       !---- Updating V_Vectors ----!
       if (deleted) call Delete_element_in_Varrays(N)

       return
    End Subroutine Delete_RefCodes_MolCrys

    !!--++
    !!--++ Subroutine Delete_RefCodes_Molec(N,Molec)
    !!--++    integer,             intent(in)     :: N
    !!--++    type(molecule_type), intent(in out) :: Molec
    !!--++
    !!--++ Overloaded
    !!--++ Delete the number of Refinable Parameter (N) on the list
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Delete_RefCodes_Molec(N,Molec)
       !---- Arguments ----!
       integer,             intent(in)     :: N
       type(molecule_type), intent(in out) :: Molec

       !---- Local Variables ----!
       logical :: deleted
       integer :: i,j

       deleted=.false.

       do j=1,3
          if (Molec%lxcentre(j) == N) then
             Molec%lxcentre(j)=0
             Molec%mxcentre(j)=0.0
             deleted=.true.
          end if
       end do

       do j=1,3
          if (Molec%lOrient(j) == N) then
             Molec%lOrient(j)=0
             Molec%mOrient(j)=0.0
             deleted=.true.
          end if
       end do

       do j=1,6
          if (Molec%lT_TLS(j) == N) then
             Molec%lT_TLS(j)=0
             Molec%mT_TLS(j)=0.0
             deleted=.true.
          end if
       end do

       do j=1,6
          if (Molec%lL_TLS(j) == N) then
             Molec%lL_TLS(j)=0
             Molec%mL_TLS(j)=0.0
             deleted=.true.
          end if
       end do

       do i=1,3
          do j=1,3
             if (Molec%lS_TLS(i,j) == N) then
                Molec%lS_TLS(i,j)=0
                Molec%mS_TLS(i,j)=0.0
                deleted=.true.
             end if
          end do
       end do

       !---- Updating ----!
       do j=1,3
          if (Molec%lxcentre(j) > N) then
             Molec%lxcentre(j)=Molec%lxcentre(j)-1
          end if
       end do

       do j=1,3
          if (Molec%lOrient(j) > N) then
             Molec%lOrient(j)=Molec%lOrient(j)-1
          end if
       end do

       do j=1,6
          if (Molec%lT_TLS(j) > N) then
             Molec%lT_TLS(j)=Molec%lT_TLS(j)-1
          end if
       end do

       do j=1,6
          if (Molec%lL_TLS(j) > N) then
             Molec%lL_TLS(j)=Molec%lL_TLS(j)-1
          end if
       end do

       do i=1,3
          do j=1,3
             if (Molec%lS_TLS(i,j) > N) then
                Molec%lS_TLS(i,j)=Molec%lS_TLS(i,j)-1
             end if
          end do
       end do

       if (molec%natoms <=0) return

       do i=1,Molec%Natoms
          do j=1,3
             if (Molec%lI_coor(j,i) == N) then
                Molec%lI_coor(j,i)=0
                Molec%mI_coor(j,i)=0.0
                deleted=.true.
             end if
          end do

          if (Molec%lbiso(i) == N) then
             Molec%lbiso(i)=0
             Molec%mbiso(i)=0.0
             deleted=.true.
          end if

          if (Molec%locc(i) == N) then
             Molec%locc(i)=0
             Molec%mocc(i)=0.0
             deleted=.true.
          end if
       end do

       !---- Updating ----!
       do i=1,Molec%Natoms
          do j=1,3
             if (Molec%lI_coor(j,i) > N) then
                Molec%lI_coor(j,i)=Molec%lI_coor(j,i)-1
             end if
          end do

          if (Molec%lbiso(i) > N) then
             Molec%lbiso(i)=Molec%lbiso(i)-1
          end if

          if (Molec%locc(i) > N) then
             Molec%locc(i)=Molec%locc(i)-1
          end if
       end do

       !---- Updating V_Vectors ----!
       if (deleted) call Delete_element_in_Varrays(N)

       return
    End Subroutine Delete_RefCodes_Molec

    !!--++
    !!--++ Subroutine Delete_RefCodes_Magdom(N, Mag_Dom)
    !!--++    integer,              intent(in)          :: N
    !!--++    type(Magnetic_Domain_type),intent(in out) :: Mag_Dom
    !!--++
    !!--++ Delete the number of Refinable Parameters (N) on the list
    !!--++ related to magnetic domains
    !!--++ Created: February - 2012
    !!
    Subroutine Delete_RefCodes_Magdom(N, Mag_Dom)
       !---- Arguments ----!
       integer,              intent(in)           :: N
       type(Magnetic_Domain_type), intent(in out) :: Mag_Dom

       !---- Local Variables ----!
       logical :: deleted
       integer :: i,j,ich

       deleted=.false.

       !---- Check is chirality is present ----!
       if (Mag_Dom%chir) then
        ich=2
       else
        ich=1
       end if

       !---- Eliminate N Parameter ----!
        do i=1,Mag_Dom%nd
         do j=1,ich
          if (Mag_Dom%Lpop(j,i) == N) then
              Mag_Dom%Mpop(j,i)=0.0
              Mag_Dom%Lpop(j,i)=0
              deleted=.true.
          end if
         end do
        end do

       !---- Updating Variables ----!
        do i=1,Mag_Dom%nd
         do j=1,ich
          if (Mag_Dom%Lpop(j,i) > N) then
              Mag_Dom%Lpop(j,i)=Mag_Dom%Lpop(j,i)-1
          end if
         end do
        end do

       !---- Updating V_Vectors ----!
       if (deleted) call Delete_element_in_Varrays(N)

       return
    End Subroutine Delete_RefCodes_Magdom

    !!----
    !!---- Subroutine Delete_RefGCodes(N, model)
    !!----   integer,                             intent(in)     :: N
    !!----   type(Nonatomic_Parameter_List_Type), intent(in out) :: model
    !!----
    !!---- Delete the number of Refinable Parameter (N) on the list of
    !!---- non atomic parameters
    !!----
    !!---- Update: November 2 - 2013
    !!
    Subroutine Delete_RefGCodes(N, model)
       !---- Arguments ----!
       integer,                             intent(in)     :: N
       type(Nonatomic_Parameter_List_Type), intent(in out) :: model

       !---- Local Variables ----!
       logical :: deleted
       integer :: i

       deleted=.false.

       !---- Eliminate N Parameter ----!
       do i=1,model%npar
          if (model%par(i)%Lcode == N) then
              model%par(i)%Lcode=0
              model%par(i)%multip=0.0
              deleted=.true.
          end if
       end do

       !---- Updating Variables ----!
       do i=1,model%npar
          if (model%par(i)%Lcode > N) then
              model%par(i)%Lcode=model%par(i)%Lcode-1
          end if
       end do

       !---- Updating V_Vectors ----!
       if (deleted) call Delete_element_in_Varrays(N)
       return
    End Subroutine Delete_RefGCodes

    Subroutine Delete_element_in_Varrays(N)
       integer, intent(in) :: N
       integer             :: i
       do i=N+1,Np_refi
          V_Vec(i-1)=V_Vec(i)
          V_Name(i-1)=V_Name(i)
          V_Bounds(:,i-1)=V_Bounds(:,i)
          V_BCon(i-1)=V_BCon(i)
          V_List(i-1)=V_List(i)
       end do
       V_Vec(np_refi)=0.0
       V_Name(np_refi)=" "
       V_Bounds(:,np_refi)=0.0
       V_BCon(np_refi)=0
       V_List(np_refi)=0
       np_refi=np_refi-1
       return
    End Subroutine Delete_element_in_Varrays
    !!--++
    !!--++ Subroutine Fill_RefCodes(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,FAtom/MolCrys/Molec,Spg)
    !!--++    integer,                       intent(in)     :: Key
    !!--++    character(len=*),              intent(in)     :: Dire
    !!--++    integer,                       intent(in)     :: Na
    !!--++    integer,                       intent(in)     :: Nb
    !!--++    real(kind=cp),                 intent(in)     :: Xl
    !!--++    real(kind=cp),                 intent(in)     :: Xu
    !!--++    real(kind=cp),                 intent(in)     :: Xs
    !!--++    integer,                       intent(in)     :: Ic
    !!--++    type(Atom_List_Type),          intent(in out) :: FAtom
    !!--++    or
    !!--++    type(molecular_Crystal_type),  intent(in out) :: MolCrys
    !!--++    or
    !!--++    type(molecule_type),           intent(in out) :: Molec
    !!--++    type(space_group_type),        intent(in)     :: Spg
    !!--++
    !!--++    (Private)
    !!--++    Write on Vectors the Information for Free Atoms
    !!--++
    !!--++ Update: March - 2005
    !!

    !!--++
    !!--++ Subroutine Fill_RefCodes_FAtom(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,FAtom,Spg)
    !!--++    integer,                       intent(in)     :: Key  !0=> nb as below,
    !!--++    character(len=*),              intent(in)     :: Dire !Var of Fix
    !!--++    integer,                       intent(in)     :: Na   !Number of the atom in the asymmetric unit
    !!--++    integer,                       intent(in)     :: Nb   !Type of individual parameter x=1,y=2,z=3,biso=4, etc
    !!--++    real(kind=cp),                 intent(in)     :: Xl   !Lower bound of parameter
    !!--++    real(kind=cp),                 intent(in)     :: Xu   !Upper bound of parameter
    !!--++    real(kind=cp),                 intent(in)     :: Xs   !Step of parameter
    !!--++    integer,                       intent(in)     :: Ic   !Boundary condition (0:fixed or 1:periodic)
    !!--++    type(Atom_List_Type),          intent(in out) :: FAtom
    !!--++    type(space_group_type),        intent(in)     :: Spg
    !!--++
    !!--++ Overloaded
    !!--++ Write on Vectors the Information for Free Atoms
    !!--++  Key=0 -> Provide information on individual parameter for atom na (nb should be given)
    !!--++  Key=1  XYZ -> fix or vary all coordinates of atom na
    !!--++  Key=2  OCC -> fix or vary occupation factor of atom na
    !!--++  Key=3  BIS -> fix or vary isotropic temperature factor of atom na
    !!--++  Key=4  BAN -> fix or vary anisotropic temperature factors of atom na
    !!--++  Key=5  ALL -> fix or vary all parameters of atom na
    !!--++
    !!--++  nb=1:3   X_ Y_ Z_
    !!--++  nb=4     Biso_
    !!--++  nb=5     Occ_
    !!--++  nb=6:11  B11_  B22_  B33_  B12_  B13_  B23_
    !!--++  nb=12    Bns_ (all Bij ...)
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Fill_RefCodes_FAtom(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,FAtom,Spg)
       !---- Arguments ----!
       integer,                       intent(in)     :: Key
       character(len=*),              intent(in)     :: Dire
       integer,                       intent(in)     :: Na
       integer,                       intent(in)     :: Nb
       real(kind=cp),                 intent(in)     :: Xl
       real(kind=cp),                 intent(in)     :: Xu
       real(kind=cp),                 intent(in)     :: Xs
       integer,                       intent(in)     :: Ic
       type(Atom_List_Type),          intent(in out) :: FAtom
       type(space_group_type),        intent(in)     :: Spg

       !---- Local variables ----!
       integer           :: j,nc,np_ini

       call init_err_refcodes()
       if (Na <= 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Number of atom no defined"
          return
       end if

       select case (dire)
          !---- FIX Directive ----!
          case ("fix")

             select case (key)
                case (0) !Key=0 -> Provide information on individual parameter for atom na (nb should be given)
                   !---- nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined"
                         return

                      case ( 1:3) ! key=0, nb=1 to 3. Fix particular coordinate x(nb) for atom na
                         !---- X_, Y_, Z_ ----!
                         if (FAtom%atom(na)%lx(nb) /=0) then
                            nc=FAtom%atom(na)%lx(nb)
                            call Delete_RefCodes(nc,FAtom)
                         end if

                      case ( 4) ! key=0, nb=4. Fix Biso for atom na
                         !---- Biso_ ----!
                         if (FAtom%atom(na)%lbiso /=0) then
                            nc=FAtom%atom(na)%lbiso
                            call Delete_RefCodes(nc,fatom)
                         end if

                      case ( 5) ! key=0, nb=5. Fix Occ for atom na
                         !---- Occ_ ----!
                         if (FAtom%atom(na)%locc /=0) then
                            nc=FAtom%atom(na)%locc
                            call Delete_RefCodes(nc,fatom)
                         end if

                      case ( 6:11) ! key=0, nb=6 to 11. Fix particula Baniso(nb) for atom na
                         !---- B11_,...,B23_ ----!
                         if (FAtom%atom(na)%lu(nb-5) /=0) then
                            nc=FAtom%atom(na)%lu(nb-5)
                            call Delete_RefCodes(nc,fatom)
                         end if

                      case (12) ! key=0, nb=12. Fix all B for atom na
                         !---- Banis_ or Bns_ ----!
                         do j=1,6
                            if (FAtom%atom(na)%lu(j) /=0) then
                               nc=FAtom%atom(na)%lu(j)
                               call Delete_RefCodes(nc,fatom)
                            end if
                         end do

                      case (13:)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined for this type of variable "
                         return
                   end select ! nb

                case (1) ! Key=1  XYZ -> Fix all coordinates of atom na
                   !---- XYZ ----!
                   do j=1,3
                      if (FAtom%atom(na)%lx(j) /=0) then
                         nc=FAtom%atom(na)%lx(j)
                         call Delete_RefCodes(nc,fatom)
                      end if
                   end do

                case (2) !  Key=2  OCC -> Fix occupation factor of atom na
                   !---- OCC ----!
                   if (FAtom%atom(na)%locc /=0) then
                      nc=FAtom%atom(na)%locc
                      call Delete_RefCodes(nc,fatom)
                   end if

                case (3) !  Key=3  BIS -> Fix isotropic temperature factor of atom na
                   !---- BIS ----!
                   if (FAtom%atom(na)%lbiso /=0) then
                      nc=FAtom%atom(na)%lbiso
                      call Delete_RefCodes(nc,fatom)
                   end if

                case (4) !  Key=4  BAN -> Fix anisotropic temperature factors of atom na
                  !---- BAN ----!
                  do j=1,6
                     if (FAtom%atom(na)%lu(j) /=0) then
                        nc=FAtom%atom(na)%lu(j)
                        call Delete_RefCodes(nc,fatom)
                     end if
                  end do

                case (5) !  Key=5  ALL -> Fix all parameters of atom na
                   !---- ALL ----!
                   do j=1,3
                      if (FAtom%atom(na)%lx(j) /=0) then
                         nc=FAtom%atom(na)%lx(j)
                         call Delete_RefCodes(nc,fatom)
                      end if
                   end do
                   if (FAtom%atom(na)%locc /=0) then
                      nc=FAtom%atom(na)%locc
                      call Delete_RefCodes(nc,fatom)
                   end if
                   if (FAtom%atom(na)%lbiso /=0) then
                      nc=FAtom%atom(na)%lbiso
                      call Delete_RefCodes(nc,fatom)
                   end if
                   do j=1,6
                      if (FAtom%atom(na)%lu(j) /=0) then
                         nc=FAtom%atom(na)%lu(j)
                         call Delete_RefCodes(nc,fatom)
                      end if
                   end do

                case (6:)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Incompatible information for this type of variable "
                   return
             end select

          !---- VARY Directive ----!
          case ("var")

             select case (key)
                case (0)

                   !---- nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined"
                         return

                      case ( 1:3) ! key=0, nb=1 to 3. Vary particular coordinate x(nb) for atom na
                         !---- X_, Y_, Z_ ----!
                         if (FAtom%atom(na)%lx(nb) ==0) then
                            FAtom%atom(na)%mx(nb)=1.0
                            call get_atompos_ctr(FAtom%atom(na)%x, Spg, np_refi,   &
                                                 FAtom%atom(na)%lx, FAtom%atom(na)%mx)
                            if (FAtom%atom(na)%lx(nb) == np_refi) then
                               V_Vec(np_refi)=FAtom%atom(na)%x(nb)
                               V_Name(np_refi)=trim(code_nam(nb))//trim(FAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if

                      case ( 4) ! key=0, nb=4. Vary Biso for atom na
                         !---- Biso_ ----!
                         if (FAtom%atom(na)%lbiso ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=FAtom%atom(na)%biso
                            V_Name(np_refi)=trim(code_nam(nb))//trim(FAtom%atom(na)%lab)
                            FAtom%atom(na)%mbiso=1.0
                            FAtom%atom(na)%lbiso=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case ( 5)  ! key=0, nb=5. Vary Occ for atom na
                         !---- Occ_ ----!
                         if (FAtom%atom(na)%locc ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=FAtom%atom(na)%occ
                            V_Name(np_refi)=trim(code_nam(nb))//trim(FAtom%atom(na)%lab)
                            FAtom%atom(na)%mocc=1.0
                            FAtom%atom(na)%locc=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case ( 6:11) ! key=0, nb=6 to 11. Vary particular Baniso(nb) for atom na
                         !---- B11_,...,B23_ ----!
                         if (FAtom%atom(na)%lu(nb-5) ==0) then
                            FAtom%atom(na)%mu(nb-5)=1.0
                            call get_atombet_ctr(FAtom%atom(na)%x,FAtom%atom(na)%u,Spg, &
                                                 np_refi,FAtom%atom(na)%lu,FAtom%atom(na)%mu)
                            if (FAtom%atom(na)%lu(nb-5) == np_refi) then
                               V_Vec(np_refi)=FAtom%atom(na)%u(nb-5)
                               V_Name(np_refi)=trim(code_nam(nb))//trim(FAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if

                      case (12)  ! key=0, nb=12. Vary all Baniso for atom na
                         !---- Banis_ ----!
                         np_ini=np_refi
                         FAtom%atom(na)%mu(:)=1.0
                         call get_atombet_ctr(FAtom%atom(na)%x,FAtom%atom(na)%u,Spg, &
                                              np_refi,FAtom%atom(na)%lu,FAtom%atom(na)%mu)
                         !write(*,"(a,2i4,3x,6i4)")"BANIS: np_ini,np_refi,LU",np_ini,np_refi,FAtom%atom(na)%lu
                         np_refi=np_ini
                         do j=1,6
                            if (FAtom%atom(na)%lu(j) > np_refi) then
                               np_refi=np_refi+1
                               FAtom%atom(na)%lu(j)=np_refi
                                  V_Vec(np_refi)=FAtom%atom(na)%u(j)
                                  V_Name(np_refi)=trim(code_nam(5+j))//trim(FAtom%atom(na)%lab)
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                            end if
                         end do

                      case (13:)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined by this type of variables"
                         return

                   end select ! nb

                case (1)  ! Key=1  XYZ -> Vary all coordinates of atom na
                   !---- XYZ ----!
                   do j=1,3
                      if (FAtom%atom(na)%lx(j) ==0) then
                         FAtom%atom(na)%mx(j)=1.0
                         call get_atompos_ctr(FAtom%atom(na)%x, Spg, np_refi,   &
                                              FAtom%atom(na)%lx, FAtom%atom(na)%mx)
                         if (FAtom%atom(na)%lx(j) == np_refi) then
                            V_Vec(np_refi)=FAtom%atom(na)%x(j)
                            V_Name(np_refi)=trim(code_nam(j))//trim(FAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case (2) !  Key=2  OCC -> Vary occupation factor of atom na
                   !---- OCC ----!
                   if (FAtom%atom(na)%locc ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=FAtom%atom(na)%occ
                      V_Name(np_refi)=trim(code_nam(5))//trim(FAtom%atom(na)%lab)
                      FAtom%atom(na)%mocc=1.0
                      FAtom%atom(na)%locc=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if

                case (3) !  Key=3  BIS -> Vary isotropic temperature factor of atom na
                   !---- BIS ----!
                   if (FAtom%atom(na)%lbiso ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=FAtom%atom(na)%biso
                      V_Name(np_refi)=trim(code_nam(4))//trim(FAtom%atom(na)%lab)
                      FAtom%atom(na)%mbiso=1.0
                      FAtom%atom(na)%lbiso=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if

                case (4) !  Key=4  BAN -> Vary all anisotropic temperature factors of atom na
                   !---- BAN ----!
                   np_ini=np_refi
                   FAtom%atom(na)%mu(:)=1.0
                   call get_atombet_ctr(FAtom%atom(na)%x,FAtom%atom(na)%u,Spg, &
                                        np_refi,FAtom%atom(na)%lu,FAtom%atom(na)%mu)
                   np_refi=np_ini
                   do j=1,6
                      if (FAtom%atom(na)%lu(j) > np_refi) then
                         np_refi=np_refi+1
                         FAtom%atom(na)%lu(j)=np_refi
                            V_Vec(np_refi)=FAtom%atom(na)%u(j)
                            V_Name(np_refi)=trim(code_nam(5+j))//trim(FAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                      end if
                   end do

                case (5) !  Key=5  ALL -> Vary all parameters of atom na
                   !---- ALL ----!
                   do j=1,3
                      if (FAtom%atom(na)%lx(j) ==0) then
                         FAtom%atom(na)%mx(j)=1.0
                         call get_atompos_ctr(FAtom%atom(na)%x, Spg, np_refi,   &
                                              FAtom%atom(na)%lx, FAtom%atom(na)%mx)
                         if (FAtom%atom(na)%lx(j) == np_refi) then
                            V_Vec(np_refi)=FAtom%atom(na)%x(j)
                            V_Name(np_refi)=trim(code_nam(j))//trim(FAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do
                   if (FAtom%atom(na)%locc ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=FAtom%atom(na)%occ
                      V_Name(np_refi)=trim(code_nam(5))//trim(FAtom%atom(na)%lab)
                      FAtom%atom(na)%mocc=1.0
                      FAtom%atom(na)%locc=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if
                   if (FAtom%atom(na)%lbiso ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=FAtom%atom(na)%biso
                      V_Name(np_refi)=trim(code_nam(4))//trim(FAtom%atom(na)%lab)
                      FAtom%atom(na)%mbiso=1.0
                      FAtom%atom(na)%lbiso=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if
                   np_ini=np_refi
                   FAtom%atom(na)%mu(:)=1.0
                   call get_atombet_ctr(FAtom%atom(na)%x,FAtom%atom(na)%u,Spg, &
                                        np_refi,FAtom%atom(na)%lu,FAtom%atom(na)%mu)
                   np_refi=np_ini
                   do j=1,6
                      if (FAtom%atom(na)%lu(j) > np_refi) then
                         np_refi=np_refi+1
                         FAtom%atom(na)%lu(j)=np_refi
                            V_Vec(np_refi)=FAtom%atom(na)%u(j)
                            V_Name(np_refi)=trim(code_nam(5+j))//trim(FAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                      end if
                   end do


                case(6:)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option Not defined by this type of variables"
                   return
             end select
       end select

       return
    End Subroutine Fill_RefCodes_FAtom

    !!--++
    !!--++ Subroutine Fill_RefCodes_FmAtom(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,FmAtom)
    !!--++    integer,                       intent(in)     :: Key
    !!--++    character(len=*),              intent(in)     :: Dire
    !!--++    integer,                       intent(in)     :: Na
    !!--++    integer,                       intent(in)     :: Nb
    !!--++    real(kind=cp),                 intent(in)     :: Xl
    !!--++    real(kind=cp),                 intent(in)     :: Xu
    !!--++    real(kind=cp),                 intent(in)     :: Xs
    !!--++    integer,                       intent(in)     :: Ic
    !!--++    type(mAtom_List_Type),         intent(in out) :: FmAtom
    !!--++
    !!--++ Write on Vectors the Information for Free Magnetic Atoms
    !!--++ magnetic clone of Fill_RefCodes_FAtom
    !!--++ Created: December - 2011
    !!--++ Updated: February - 2012
    !!--++
    Subroutine Fill_RefCodes_FmAtom(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,FmAtom)
       !---- Arguments ----!
       integer,                       intent(in)     :: Key
       character(len=*),              intent(in)     :: Dire
       integer,                       intent(in)     :: Na
       integer,                       intent(in)     :: Nb
       real(kind=cp),                 intent(in)     :: Xl
       real(kind=cp),                 intent(in)     :: Xu
       real(kind=cp),                 intent(in)     :: Xs
       integer,                       intent(in)     :: Ic
       type(mAtom_List_Type),         intent(in out) :: FmAtom

       !---- Local variables ----!
       integer           :: j,nc,ik

       call init_err_refcodes()
       if (Na <= 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Number of atom no defined"
          return
       end if

       !---- Get im, number of the magnetic matrices/irrep set
       ik=FmAtom%atom(na)%nvk

       select case (dire)
          !---- FIX Directive ----!
          case ("fix")

             select case (key)
                case (0)
                   !---- nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined"
                         return

                      case ( 1:3) ! key=0, nb=1 to 3, particular coordinate x(nb) for atom na
                         !---- Rx_, Ry_, Rz_ ----!
                         if (FmAtom%atom(na)%lSkR(nb,ik) /=0) then
                            nc=FmAtom%atom(na)%lSkR(nb,ik)
                            call Delete_RefCodes(nc,FmAtom)
                         end if

                      case ( 4:6)
                         !---- Ix_, Iy_, Iz_ ----!
                         if (FmAtom%atom(na)%lSkI(nb-3,ik) /=0) then
                            nc=FmAtom%atom(na)%lSkI(nb-3,ik)
                            call Delete_RefCodes(nc,FmAtom)
                         end if

                      case ( 7:9)
                         !---- Rm_, Rphi_, Rth_ ----!
                         if (FmAtom%atom(na)%lSkR(nb-6,ik) /=0) then
                            nc=FmAtom%atom(na)%lSkR(nb-6,ik)
                            call Delete_RefCodes(nc,FmAtom)
                         end if

                      case (10:12)
                         !---- Im_, Iphi_, Ith_ ----!
                         if (FmAtom%atom(na)%lSkI(nb-9,ik) /=0) then
                            nc=FmAtom%atom(na)%lSkI(nb-9,ik)
                            call Delete_RefCodes(nc,FmAtom)
                         end if

                      case (13)
                         !---- MagPh_ ----!
                            if (FmAtom%atom(na)%lmphas(ik) /=0) then
                               nc=FmAtom%atom(na)%lmphas(ik)
                               call Delete_RefCodes(nc,FmAtom)
                            end if

                      case (14:22)
                         !---- C1_,..C12_ ----!
                            if (FmAtom%atom(na)%lbas(nb-13,ik) /=0) then
                               nc=FmAtom%atom(na)%lbas(nb-13,ik)
                               call Delete_RefCodes(nc,FmAtom)
                            end if

                      case (23)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined for this type of variable "
                         return
                   end select ! nb

                case (1)
                   !---- Rxyz ----!
                   do j=1,3
                      if (FmAtom%atom(na)%lSkR(j,ik) /=0) then
                         nc=FmAtom%atom(na)%lSkR(j,ik)
                         call Delete_RefCodes(nc,FmAtom)
                      end if
                   end do

                case (2)
                   !---- Ixyz ----!
                   do j=1,3
                      if (FmAtom%atom(na)%lSkI(j,ik) /=0) then
                         nc=FmAtom%atom(na)%lSkI(j,ik)
                         call Delete_RefCodes(nc,FmAtom)
                      end if
                   end do

                case (3)
                   !---- Mxyz ----!
                   do j=1,3
                      if (FmAtom%atom(na)%lSkR(j,ik) /=0) then
                         nc=FmAtom%atom(na)%lSkR(j,ik)
                         call Delete_RefCodes(nc,FmAtom)
                      end if
                   end do
                   do j=1,3
                      if (FmAtom%atom(na)%lSkI(j,ik) /=0) then
                         nc=FmAtom%atom(na)%lSkI(j,ik)
                         call Delete_RefCodes(nc,FmAtom)
                      end if
                   end do

                case (4:)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Incompatible information for this type of variable "
                   return
             end select

          !---- VARY Directive ----!
          case ("var")

              select case (key)
                case (0)

                   !---- nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined"
                         return

                      case ( 1:3)
                         !---- Rx_, Ry_, Rz_ ----!
                         if (FmAtom%atom(na)%lSkR(nb,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lSkR(nb,ik)=np_refi
                            FmAtom%atom(na)%mSkR(nb,ik)=1.0

                            if (FmAtom%atom(na)%lSkR(nb,ik) == np_refi) then
                               V_Vec(np_refi)=FmAtom%atom(na)%SkR(nb,ik)
                               V_Name(np_refi)=trim(mcode_nam(nb))//trim(FmAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if

                      case ( 4:6)
                         !---- Ix_, Iy_, Iz_ ----!
                         if (FmAtom%atom(na)%lSkI(nb-3,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lSkI(nb-3,ik)=np_refi
                            FmAtom%atom(na)%mSkI(nb-3,ik)=1.0

                            if (FmAtom%atom(na)%lSkI(nb-3,ik) == np_refi) then
                               V_Vec(np_refi)=FmAtom%atom(na)%SkI(nb-3,ik)
                               V_Name(np_refi)=trim(mcode_nam(nb))//trim(FmAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                             else
                               np_refi=np_refi-1
                            end if
                         end if

                      case ( 7:9)
                         !---- Rm_, Rphi_, Rth_ ----!
                         if (FmAtom%atom(na)%lSkR(nb-6,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lSkR(nb-6,ik)=np_refi
                            FmAtom%atom(na)%mSkR(nb-6,ik)=1.0

                            if (FmAtom%atom(na)%lSkR(nb-6,ik) == np_refi) then
                               V_Vec(np_refi)=FmAtom%atom(na)%Spher_SkR(nb-6,ik)
                               V_Name(np_refi)=trim(mcode_nam(nb))//trim(FmAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if

                      case (10:12)
                         !---- Im_, Iphi_, Ith_ ----!
                         if (FmAtom%atom(na)%lSkI(nb-9,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lSkI(nb-9,ik)=np_refi
                            FmAtom%atom(na)%mSkI(nb-9,ik)=1.0

                            if (FmAtom%atom(na)%lSkI(nb-9,ik) == np_refi) then
                               V_Vec(np_refi)=FmAtom%atom(na)%Spher_SkI(nb-9,ik)
                               V_Name(np_refi)=trim(mcode_nam(nb))//trim(FmAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                             else
                               np_refi=np_refi-1
                            end if
                         end if

                      case (13)
                         !---- MagPh_ ----!
                         if (FmAtom%atom(na)%lmphas(ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%mmphas(ik)=1.0
                            FmAtom%atom(na)%lmphas(ik)=np_refi

                            V_Vec(np_refi)=FmAtom%atom(na)%mphas(ik)
                            V_Name(np_refi)=trim(mcode_nam(nb))//trim(FmAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case (14:22)
                         !---- C1_,..C12_ ----!
                         if (FmAtom%atom(na)%lbas(nb-13,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lbas(nb-13,ik)=np_refi
                            FmAtom%atom(na)%mbas(nb-13,ik)=1.0

                            if (FmAtom%atom(na)%lbas(nb-13,ik) == np_refi) then
                               V_Vec(np_refi)=FmAtom%atom(na)%cbas(nb-13,ik)
                               V_Name(np_refi)=trim(mcode_nam(nb))//trim(FmAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                             else
                               np_refi=np_refi-1
                            end if
                         end if

                      case (23:)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined by this type of variables"
                         return

                   end select ! nb

                case (1)
                   !---- Rxyz ----!
                   do j=1,3
                      if (FmAtom%atom(na)%lSkR(j,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lSkR(j,ik)=np_refi
                            FmAtom%atom(na)%mSkR(j,ik)=1.0

                         if (FmAtom%atom(na)%lSkR(j,ik) == np_refi) then
                            V_Vec(np_refi)=FmAtom%atom(na)%SkR(j,ik)
                            V_Name(np_refi)=trim(mcode_nam(j))//trim(FmAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case (2)
                   !---- Ixyz ----!
                   do j=1,3
                      if (FmAtom%atom(na)%lSkI(j,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lSkI(j,ik)=np_refi
                            FmAtom%atom(na)%mSkI(j,ik)=1.0

                         if (FmAtom%atom(na)%lSkI(j,ik) == np_refi) then
                            V_Vec(np_refi)=FmAtom%atom(na)%SkI(j,ik)
                            V_Name(np_refi)=trim(mcode_nam(j+3))//trim(FmAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case (3)
                    !---- Mxyz ----!
                   do j=1,3
                      if (FmAtom%atom(na)%lSkR(j,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lSkR(j,ik)=np_refi
                            FmAtom%atom(na)%mSkR(j,ik)=1.0

                         if (FmAtom%atom(na)%lSkR(j,ik) == np_refi) then
                            V_Vec(np_refi)=FmAtom%atom(na)%SkR(j,ik)
                            V_Name(np_refi)=trim(mcode_nam(j))//trim(FmAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                   do j=1,3
                      if (FmAtom%atom(na)%lSkI(j,ik) ==0) then

                            np_refi=np_refi+1
                            FmAtom%atom(na)%lSkI(j,ik)=np_refi
                            FmAtom%atom(na)%mSkI(j,ik)=1.0

                         if (FmAtom%atom(na)%lSkI(j,ik) == np_refi) then
                            V_Vec(np_refi)=FmAtom%atom(na)%SkI(j,ik)
                            V_Name(np_refi)=trim(mcode_nam(j+3))//trim(FmAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case (4:)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option Not defined by this type of variables"
                   return
             end select
       end select

       return
    End Subroutine Fill_RefCodes_FmAtom

    !!--++
    !!--++ Subroutine Fill_RefCodes_Molcrys(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Molcrys,NMol)
    !!--++    integer,                      intent(in)     :: Key
    !!--++    character(len=*),             intent(in)     :: Dire
    !!--++    integer,                      intent(in)     :: Na
    !!--++    integer,                      intent(in)     :: Nb
    !!--++    real(kind=cp),                intent(in)     :: Xl
    !!--++    real(kind=cp),                intent(in)     :: Xu
    !!--++    real(kind=cp),                intent(in)     :: Xs
    !!--++    integer,                      intent(in)     :: Ic
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++    integer,                      intent(in)     :: NMol
    !!--++
    !!--++ Overloaded
    !!--++ Write on Vectors the Information for Free Atoms
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Fill_RefCodes_Molcrys(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Molcrys,Nmol)
       !---- Arguments ----!
       integer,                      intent(in)     :: Key
       character(len=*),             intent(in)     :: Dire
       integer,                      intent(in)     :: Na
       integer,                      intent(in)     :: Nb
       real(kind=cp),                intent(in)     :: Xl
       real(kind=cp),                intent(in)     :: Xu
       real(kind=cp),                intent(in)     :: Xs
       integer,                      intent(in)     :: Ic
       type(molecular_Crystal_type), intent(in out) :: MolCrys
       integer,                      intent(in)     :: NMol

       !---- Local variables ----!
       character(len=2) :: car
       integer          :: i, j, nc, naa

       call init_err_refcodes()
       if (Na <= 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Number of atom no defined"
          return
       end if

       select case (dire)
          !---- FIX Directive ----!
          case ("fix")

             select case (key)
                case (0)
                   !---- Nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case ( 1:3)
                         !---- X_, Y_, Z_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lx(nb) /=0) then
                               nc=molcrys%atm(na)%lx(nb)
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%lI_coor(nb,naa) /=0) then
                                  nc=molcrys%mol(i)%lI_coor(nb,naa)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                            end do
                         end if

                      case ( 4) ! key=0, nb=4 Biso for atom na
                         !---- Biso_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lbiso /=0) then
                               nc=molcrys%atm(na)%lbiso
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%lbiso(naa) /=0) then
                                  nc=molcrys%mol(i)%lbiso(naa)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                            end do
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%locc /=0) then
                               nc=molcrys%atm(na)%locc
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%locc(naa) /=0) then
                                  nc=molcrys%mol(i)%locc(naa)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                            end do
                         end if

                      case ( 6:11)
                         !---- B11_, ..., B23_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lu(nb-5) /=0) then
                               nc=molcrys%atm(na)%lu(nb-5)
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         else
                            err_refcodes=.true.
                            ERR_RefCodes_Mess="Option Not defined"
                            return
                         end if

                      case (12)
                         !---- Banis_ ----!
                         if (na <= molcrys%n_free) then
                            do j=1,6
                               if (molcrys%atm(na)%lu(j) /=0) then
                                  nc=molcrys%atm(na)%lu(j)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                            end do
                         else
                            err_refcodes=.true.
                            ERR_RefCodes_Mess="Option Not defined"
                            return
                         end if

                      case (13:15)
                         !---- Xc_, Yc_, Zc_ ----!
                         select case (nmol)
                            case (-1)
                               err_refcodes=.true.
                               ERR_RefCodes_Mess="Option Not defined"
                               return

                            case (0)
                               do i=1,molcrys%n_mol
                                  if (molcrys%mol(i)%lxcentre(nb-12) /=0) then
                                     nc=molcrys%mol(i)%lxcentre(nb-12)
                                     call Delete_RefCodes(nc,MolCrys)
                                  end if
                               end do

                            case (1:)
                               if (molcrys%mol(nmol)%lxcentre(nb-12) /=0) then
                                  nc=molcrys%mol(nmol)%lxcentre(nb-12)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                         end select

                      case (16:18)
                         !---- Theta_, Phi_, Chi_ ----!
                         select case (nmol)
                            case (-1)
                               err_refcodes=.true.
                               ERR_RefCodes_Mess="Option Not defined"
                               return

                            case (0)
                               do i=1,molcrys%n_mol
                                  if (molcrys%mol(i)%lOrient(nb-15) /=0) then
                                     nc=molcrys%mol(i)%lOrient(nb-15)
                                     call Delete_RefCodes(nc,MolCrys)
                                  end if
                               end do

                            case (1:)
                               if (molcrys%mol(nmol)%lOrient(nb-15) /=0) then
                                  nc=molcrys%mol(nmol)%lOrient(nb-15)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                         end select

                      case (19:21)
                         !!! Not yet implemented !!!

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   if (na <= molcrys%n_free) then
                      do j=1,3
                         if (molcrys%atm(na)%lx(j) /=0) then
                            nc=molcrys%atm(na)%lx(j)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         do j=1,3
                            if (molcrys%mol(i)%lI_coor(j,naa) /=0) then
                               nc=molcrys%mol(i)%lI_coor(j,naa)
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         end do
                      end do
                   end if

                case (2)
                   !---- OCC ----!
                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%locc /=0) then
                         nc=molcrys%atm(na)%locc
                         call Delete_RefCodes(nc,MolCrys)
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%locc(naa) /=0) then
                            nc=molcrys%mol(i)%locc(naa)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                case (3)
                   !---- BIS ----!
                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%lbiso /=0) then
                         nc=molcrys%atm(na)%lbiso
                         call Delete_RefCodes(nc,MolCrys)
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%lbiso(naa) /=0) then
                            nc=molcrys%mol(i)%lbiso(naa)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                case (4)
                   !---- BAN ----!
                   if (na <= molcrys%n_free) then
                      do j=1,6
                         if (molcrys%atm(na)%lu(j) /=0) then
                            nc=molcrys%atm(na)%lu(j)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   else
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Option Not defined"
                      return
                   end if

                case (5)
                   !---- ALL ----!
                   if (na <= molcrys%n_free) then
                      do j=1,3
                         if (molcrys%atm(na)%lx(j) /=0) then
                            nc=molcrys%atm(na)%lx(j)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         do j=1,3
                            if (molcrys%mol(i)%lI_coor(j,naa) /=0) then
                               nc=molcrys%mol(i)%lI_coor(j,naa)
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         end do
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%locc /=0) then
                         nc=molcrys%atm(na)%locc
                         call Delete_RefCodes(nc,MolCrys)
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%locc(naa) /=0) then
                            nc=molcrys%mol(i)%locc(naa)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%lbiso /=0) then
                         nc=molcrys%atm(na)%lbiso
                         call Delete_RefCodes(nc,MolCrys)
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%lbiso(naa) /=0) then
                            nc=molcrys%mol(i)%lbiso(naa)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      do j=1,6
                         if (molcrys%atm(na)%lu(j) /=0) then
                            nc=molcrys%atm(na)%lu(j)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            do j=1,3
                               if (molcrys%mol(i)%lxcentre(j) /=0) then
                                  nc=molcrys%mol(i)%lxcentre(j)
                                  call Delete_RefCodes(nc,molcrys)
                               end if
                            end do

                            do j=1,3
                               if (molcrys%mol(i)%lOrient(j) /=0) then
                                  nc=molcrys%mol(i)%lOrient(j)
                                  call Delete_RefCodes(nc,molcrys)
                               end if
                            end do
                         end do

                      case (1:)
                         do j=1,3
                            if (molcrys%mol(nmol)%lxcentre(j) /=0) then
                               nc=molcrys%mol(nmol)%lxcentre(j)
                               call Delete_RefCodes(nc,molcrys)
                            end if
                         end do

                         do j=1,3
                            if (molcrys%mol(nmol)%lOrient(j) /=0) then
                               nc=molcrys%mol(nmol)%lOrient(j)
                               call Delete_RefCodes(nc,molcrys)
                            end if
                         end do
                   end select

                   !!! Not yet Implemented !!!

                case (6)
                   !---- CEN ----!
                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            do j=1,3
                               if (molcrys%mol(i)%lxcentre(j) /=0) then
                                  nc=molcrys%mol(i)%lxcentre(j)
                                  call Delete_RefCodes(nc,molcrys)
                               end if
                            end do
                         end do

                      case (1:)
                         do j=1,3
                            if (molcrys%mol(nmol)%lxcentre(j) /=0) then
                               nc=molcrys%mol(nmol)%lxcentre(j)
                               call Delete_RefCodes(nc,molcrys)
                            end if
                         end do
                   end select

                case (7)
                   !---- ORI ----!
                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            do j=1,3
                               if (molcrys%mol(i)%lOrient(j) /=0) then
                                  nc=molcrys%mol(i)%lOrient(j)
                                  call Delete_RefCodes(nc,molcrys)
                               end if
                            end do
                         end do

                      case (1:)
                         do j=1,3
                            if (molcrys%mol(nmol)%lOrient(j) /=0) then
                              nc=molcrys%mol(nmol)%lOrient(j)
                               call Delete_RefCodes(nc,molcrys)
                            end if
                         end do
                   end select

                case (8)
                   !---- THE ----!
                   !!! Not yet Implemented !!!!
             end select

          !---- VARY Directive ----!
          case ("var")

             select case (key)
                case (0)
                   !---- Nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case ( 1:3)
                         !---- X_, Y_, Z_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lx(nb) ==0) then
                               molcrys%atm(na)%mx(nb)=1.0
                               call get_atompos_ctr(molcrys%atm(na)%x,   &
                                                    molcrys%Spg,np_refi, &
                                                    molcrys%atm(na)%lx,  &
                                                    molcrys%atm(na)%mx)
                               if (molcrys%atm(na)%lx(nb) == np_refi) then
                                  V_Vec(np_refi)=molcrys%atm(na)%x(nb)
                                  V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%atm(na)%lab)
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if

                               if (molcrys%mol(i)%lI_coor(nb,naa) ==0) then
                                  molcrys%mol(i)%mI_coor(nb,naa)=1.0
                                  call get_atompos_ctr(molcrys%mol(i)%I_Coor(:,naa),  &
                                                       molcrys%Spg, np_refi,          &
                                                       molcrys%mol(i)%lI_coor(:,naa), &
                                                       molcrys%mol(i)%mI_coor(:,naa))
                                  if (molcrys%mol(i)%lI_coor(nb,naa) == np_refi) then
                                     V_Vec(np_refi)=molcrys%mol(i)%I_Coor(nb,naa)
                                     V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%mol(i)%AtName(naa))
                                     V_Bounds(1,np_refi)=xl
                                     V_Bounds(2,np_refi)=xu
                                     V_Bounds(3,np_refi)=xs
                                     V_BCon(np_refi)=ic
                                     V_list(np_refi)=na
                                  else
                                     np_refi=np_refi-1
                                  end if
                               end if
                            end do
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lbiso ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%atm(na)%biso
                               V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%atm(na)%lab)
                               molcrys%atm(na)%mbiso=1.0
                               molcrys%atm(na)%lbiso=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%lbiso(naa) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%biso(naa)
                                  V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%mol(i)%AtName(naa))
                                  molcrys%mol(i)%mbiso(naa)=1.0
                                  molcrys%mol(i)%lbiso(naa)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               end if
                            end do
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%locc ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%atm(na)%occ
                               V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%atm(na)%lab)
                               molcrys%atm(na)%mocc=1.0
                               molcrys%atm(na)%locc=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%locc(naa) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%occ(naa)
                                  V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%mol(i)%AtName(naa))
                                  molcrys%mol(i)%mocc(naa)=1.0
                                  molcrys%mol(i)%locc(naa)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               end if
                            end do
                         end if

                      case ( 6:11)
                         !---- B11_, ..., B23_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lu(nb-5) ==0) then
                               molcrys%atm(na)%mu(nb-5)=1.0
                               call get_atombet_ctr(molcrys%atm(na)%x,molcrys%atm(na)%u,molcrys%Spg, &
                                                   np_refi,molcrys%atm(na)%lu,molcrys%atm(na)%mu)
                               if (molcrys%atm(na)%lu(nb-5) == np_refi) then
                                  V_Vec(np_refi)=molcrys%atm(na)%u(nb-5)
                                  V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%atm(na)%lab)
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         else
                            err_refcodes=.true.
                            ERR_RefCodes_Mess="Option Not defined"
                            return
                         end if

                      case (12)
                         !---- Banis_ ----!
                         if (na <= molcrys%n_free) then
                            do j=1,6
                               if (molcrys%atm(na)%lu(j) ==0) then
                                  molcrys%atm(na)%mu(j)=1.0
                                  call get_atombet_ctr(molcrys%atm(na)%x,molcrys%atm(na)%u,molcrys%Spg, &
                                                       np_refi,molcrys%atm(na)%lu,molcrys%atm(na)%mu)
                                  if (molcrys%atm(na)%lu(j) == np_refi) then
                                     V_Vec(np_refi)=molcrys%atm(na)%u(j)
                                     V_Name(np_refi)=trim(code_nam(5+j))//trim(molcrys%atm(na)%lab)
                                     V_Bounds(1,np_refi)=xl
                                     V_Bounds(2,np_refi)=xu
                                     V_Bounds(3,np_refi)=xs
                                     V_BCon(np_refi)=ic
                                     V_list(np_refi)=na
                                  else
                                     np_refi=np_refi-1
                                  end if
                               end if
                            end do
                         else
                            err_refcodes=.true.
                            ERR_RefCodes_Mess="Option Not defined"
                            return
                         end if

                      case (13:15)
                         !---- Xc_, Yc_, Zc_ ----!
                         select case (nmol)
                            case (-1)
                               err_refcodes=.true.
                               ERR_RefCodes_Mess="Option Not defined"
                               return

                            case (0)
                               do i=1,molcrys%n_mol
                                  write(unit=car,fmt="(i2)") i
                                  car=adjustl(car)
                                  if (molcrys%mol(i)%lxcentre(nb-12) ==0) then
                                     np_refi=np_refi+1
                                     V_Vec(np_refi)=molcrys%mol(i)%xcentre(nb-12)
                                     V_Name(np_refi)=trim(code_nam(nb))//"Mol"//trim(car)
                                     molcrys%mol(i)%mxcentre(nb-12)=1.0
                                     molcrys%mol(i)%lxcentre(nb-12)=np_refi
                                     V_Bounds(1,np_refi)=xl
                                     V_Bounds(2,np_refi)=xu
                                     V_Bounds(3,np_refi)=xs
                                     V_BCon(np_refi)=ic
                                     V_list(np_refi)=-i
                                  end if
                               end do

                            case (1:)
                               write(unit=car,fmt="(i2)") nmol
                               car=adjustl(car)
                               if (molcrys%mol(nmol)%lxcentre(nb-12) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(nmol)%xcentre(nb-12)
                                  V_Name(np_refi)=trim(code_nam(nb))//"entre_Mol"//trim(car)
                                  molcrys%mol(nmol)%mxcentre(nb-12)=1.0
                                  molcrys%mol(nmol)%lxcentre(nb-12)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-nmol
                               end if
                         end select

                      case (16:18)
                         !---- Theta_, Phi_, Chi_ ----!
                         select case (nmol)
                            case (-1)
                               err_refcodes=.true.
                               ERR_RefCodes_Mess="Option Not defined"
                               return

                            case (0)
                               do i=1,molcrys%n_mol
                                  write(unit=car,fmt="(i2)") i
                                  car=adjustl(car)
                                  if (molcrys%mol(i)%lOrient(nb-15) ==0) then
                                     np_refi=np_refi+1
                                     V_Vec(np_refi)=molcrys%mol(i)%Orient(nb-15)
                                     V_Name(np_refi)=trim(code_nam(nb))//"Orient_Mol"//trim(car)
                                     molcrys%mol(i)%mOrient(nb-15)=1.0
                                     molcrys%mol(i)%lOrient(nb-15)=np_refi
                                     V_Bounds(1,np_refi)=xl
                                     V_Bounds(2,np_refi)=xu
                                     V_Bounds(3,np_refi)=xs
                                     V_BCon(np_refi)=ic
                                     V_list(np_refi)=-i
                                  end if
                               end do

                            case (1:)
                               write(unit=car,fmt="(i2)") nmol
                               car=adjustl(car)
                               if (molcrys%mol(nmol)%lOrient(nb-15) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(nmol)%Orient(nb-15)
                                  V_Name(np_refi)=trim(code_nam(nb))//"Orient_Mol"//trim(car)
                                  molcrys%mol(nmol)%mOrient(nb-15)=1.0
                                  molcrys%mol(nmol)%lOrient(nb-15)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-nmol
                               end if
                         end select

                      case (19:21)
                         !!! Not Yet Implemented !!!

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   if (na <= molcrys%n_free) then
                      do j=1,3
                         if (molcrys%atm(na)%lx(j) ==0) then
                            molcrys%atm(na)%mx(j)=1.0
                            call get_atompos_ctr(molcrys%atm(na)%x,   &
                                                 molcrys%Spg,np_refi, &
                                                 molcrys%atm(na)%lx,  &
                                                 molcrys%atm(na)%mx)
                            if (molcrys%atm(na)%lx(j) == np_refi) then
                               V_Vec(np_refi)=molcrys%atm(na)%x(j)
                               V_Name(np_refi)=trim(code_nam(j))//trim(molcrys%atm(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if
                      end do
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                             naa=naa-molcrys%mol(i)%natoms
                             cycle
                         end if
                         do j=1,3
                            if (molcrys%mol(i)%lI_coor(j,naa) ==0) then
                               molcrys%mol(i)%mI_coor(nb,naa)=1.0
                               call get_atompos_ctr(molcrys%mol(i)%I_Coor(:,naa), &
                                                    molcrys%Spg,np_refi,   &
                                                    molcrys%mol(i)%lI_coor(:,naa), &
                                                    molcrys%mol(i)%mI_coor(:,naa))
                               if (molcrys%mol(i)%lI_coor(j,naa) == np_refi) then
                                  V_Vec(np_refi)=molcrys%mol(i)%I_Coor(nb,naa)
                                  V_Name(np_refi)=trim(code_nam(j))//trim(molcrys%mol(i)%AtName(naa))
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         end do
                      end do
                   end if

                case (2)
                   !---- OCC ----!
                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%locc ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molcrys%atm(na)%occ
                         V_Name(np_refi)=trim(code_nam(5))//trim(molcrys%atm(na)%lab)
                         molcrys%atm(na)%mocc=1.0
                         molcrys%atm(na)%locc=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=na
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if

                         if (molcrys%mol(i)%locc(naa) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molcrys%mol(i)%Occ(naa)
                            V_Name(np_refi)=trim(code_nam(5))//trim(molcrys%mol(i)%AtName(naa))
                            molcrys%mol(i)%mocc(naa)=1.0
                            molcrys%mol(i)%locc(naa)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if
                      end do
                   end if

                case (3)
                   !---- BIS ----!
                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%lbiso ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molcrys%atm(na)%biso
                         V_Name(np_refi)=trim(code_nam(4))//trim(molcrys%atm(na)%lab)
                         molcrys%atm(na)%mbiso=1.0
                         molcrys%atm(na)%lbiso=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=na
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%lbiso(naa) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molcrys%mol(i)%biso(naa)
                            V_Name(np_refi)=trim(code_nam(4))//trim(molcrys%mol(i)%AtName(naa))
                            molcrys%mol(i)%mbiso(naa)=1.0
                            molcrys%mol(i)%lbiso(naa)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if
                      end do
                   end if

                case (4)
                   !---- BAN ----!
                   if (na <= molcrys%n_free) then
                      do j=1,6
                         if (molcrys%atm(na)%lu(j) ==0) then
                            molcrys%atm(na)%mu(j)=1.0
                            call get_atombet_ctr(molcrys%atm(na)%x,molcrys%atm(na)%u,molcrys%Spg, &
                                                 np_refi,molcrys%atm(na)%lu,molcrys%atm(na)%mu)
                            if (molcrys%atm(na)%lu(j) == np_refi) then
                               V_Vec(np_refi)=molcrys%atm(na)%u(j)
                               V_Name(np_refi)=trim(code_nam(5+j))//trim(molcrys%atm(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if
                      end do
                   else
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Option Not defined"
                      return
                   end if

                case (5)
                   !---- ALL ----!
                   if (na <= molcrys%n_free) then
                      do j=1,3
                         if (molcrys%atm(na)%lx(j) ==0) then
                            molcrys%atm(na)%mx(j)=1.0
                            call get_atompos_ctr(molcrys%atm(na)%x,   &
                                                 molcrys%Spg,np_refi, &
                                                 molcrys%atm(na)%lx,  &
                                                 molcrys%atm(na)%mx)
                            if (molcrys%atm(na)%lx(j) == np_refi) then
                               V_Vec(np_refi)=molcrys%atm(na)%x(j)
                               V_Name(np_refi)=trim(code_nam(j))//trim(molcrys%atm(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if
                      end do
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                             naa=naa-molcrys%mol(i)%natoms
                             cycle
                         end if

                         do j=1,3
                            if (molcrys%mol(i)%lI_coor(j,naa) ==0) then
                               molcrys%mol(i)%mI_coor(nb,naa)=1.0
                               call get_atompos_ctr(molcrys%mol(i)%I_Coor(:,naa), &
                                                    molcrys%Spg,np_refi,   &
                                                    molcrys%mol(i)%lI_coor(:,naa), &
                                                    molcrys%mol(i)%mI_coor(:,naa))
                               if (molcrys%mol(i)%lI_coor(j,naa) == np_refi) then
                                  V_Vec(np_refi)=molcrys%mol(i)%I_Coor(nb,naa)
                                  V_Name(np_refi)=trim(code_nam(j))//trim(molcrys%mol(i)%AtName(naa))
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         end do
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%locc ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molcrys%atm(na)%occ
                         V_Name(np_refi)=trim(code_nam(5))//trim(molcrys%atm(na)%lab)
                         molcrys%atm(na)%mocc=1.0
                         molcrys%atm(na)%locc=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=na
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%locc(naa) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molcrys%mol(i)%Occ(naa)
                            V_Name(np_refi)=trim(code_nam(5))//trim(molcrys%mol(i)%AtName(naa))
                            molcrys%mol(i)%mocc(naa)=1.0
                            molcrys%mol(i)%locc(naa)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%lbiso ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molcrys%atm(na)%biso
                         V_Name(np_refi)=trim(code_nam(4))//trim(molcrys%atm(na)%lab)
                         molcrys%atm(na)%mbiso=1.0
                         molcrys%atm(na)%lbiso=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=na
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if

                         if (molcrys%mol(i)%lbiso(naa) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molcrys%mol(i)%biso(naa)
                            V_Name(np_refi)=trim(code_nam(4))//trim(molcrys%mol(i)%AtName(naa))
                            molcrys%mol(i)%mbiso(naa)=1.0
                            molcrys%mol(i)%lbiso(naa)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      do j=1,6
                         if (molcrys%atm(na)%lu(j) ==0) then
                            molcrys%atm(na)%mu(j)=1.0
                            call get_atombet_ctr(molcrys%atm(na)%x,molcrys%atm(na)%u,molcrys%Spg, &
                                                 np_refi,molcrys%atm(na)%lu,molcrys%atm(na)%mu)
                            if (molcrys%atm(na)%lu(j) == np_refi) then
                               V_Vec(np_refi)=molcrys%atm(na)%u(j)
                               V_Name(np_refi)=trim(code_nam(5+j))//trim(molcrys%atm(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if
                      end do
                   end if

                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            write(unit=car,fmt="(i2)") i
                            car=adjustl(car)
                            do j=1,3
                               if (molcrys%mol(i)%lxcentre(j) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%xcentre(j)
                                  V_Name(np_refi)=trim(code_nam(12+j))//"entre_Mol"//trim(car)
                                  molcrys%mol(i)%mxcentre(j)=1.0
                                  molcrys%mol(i)%lxcentre(j)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-i
                               end if
                            end do
                         end do

                      case (1:)
                         write(unit=car,fmt="(i2)") nmol
                         car=adjustl(car)
                         do j=1,3
                            if (molcrys%mol(nmol)%lxcentre(j) ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%mol(nmol)%xcentre(j)
                               V_Name(np_refi)=trim(code_nam(12+j))//"entre_Mol"//trim(car)
                               molcrys%mol(nmol)%mxcentre(j)=1.0
                               molcrys%mol(nmol)%lxcentre(j)=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=-nmol
                            end if
                         end do
                   end select

                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            write(unit=car,fmt="(i2)") i
                            car=adjustl(car)
                            do j=1,3
                               if (molcrys%mol(i)%lOrient(j) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%Orient(j)
                                  V_Name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"//trim(car)
                                  molcrys%mol(i)%mOrient(j)=1.0
                                  molcrys%mol(i)%lOrient(j)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-i
                               end if
                            end do
                         end do

                      case (1:)
                         write(unit=car,fmt="(i2)") nmol
                         car=adjustl(car)
                         do j=1,3
                            if (molcrys%mol(nmol)%lOrient(j) ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%mol(nmol)%Orient(j)
                               V_Name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"//trim(car)
                               molcrys%mol(nmol)%mOrient(j)=1.0
                               molcrys%mol(nmol)%lOrient(j)=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=-nmol
                            end if
                         end do
                   end select

                case (6)
                   !---- CEN ----!
                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            write(unit=car,fmt="(i2)") i
                            car=adjustl(car)
                            do j=1,3
                               if (molcrys%mol(i)%lxcentre(j) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%xcentre(j)
                                  V_Name(np_refi)=trim(code_nam(12+j))//"entre_Mol"//trim(car)
                                  molcrys%mol(i)%mxcentre(j)=1.0
                                  molcrys%mol(i)%lxcentre(j)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-i
                               end if
                            end do
                         end do

                      case (1:)
                         write(unit=car,fmt="(i2)") nmol
                         car=adjustl(car)
                         do j=1,3
                            if (molcrys%mol(nmol)%lxcentre(j) ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%mol(nmol)%xcentre(j)
                               V_Name(np_refi)=trim(code_nam(12+j))//"entre_Mol"//trim(car)
                               molcrys%mol(nmol)%mxcentre(j)=1.0
                               molcrys%mol(nmol)%lxcentre(j)=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=-nmol
                            end if
                         end do
                   end select

                case (7)
                   !---- ORI  ----!
                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            write(unit=car,fmt="(i2)") i
                            car=adjustl(car)
                            do j=1,3
                               if (molcrys%mol(i)%lOrient(j) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%Orient(j)
                                  V_Name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"//trim(car)
                                  molcrys%mol(i)%mOrient(j)=1.0
                                  molcrys%mol(i)%lOrient(j)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-i
                               end if
                            end do
                         end do

                      case (1:)
                         write(unit=car,fmt="(i2)") nmol
                         car=adjustl(car)
                         do j=1,3
                            if (molcrys%mol(nmol)%lOrient(j) ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%mol(nmol)%Orient(j)
                               V_Name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"//trim(car)
                               molcrys%mol(nmol)%mOrient(j)=1.0
                               molcrys%mol(nmol)%lOrient(j)=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=-nmol
                            end if
                         end do
                   end select

                case (8)
                   !!! Not yet implemented !!!

             end select
       end select

       return
    End Subroutine Fill_RefCodes_Molcrys

    !!--++
    !!--++ Subroutine Fill_RefCodes_Molec(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Molec,Spg)
    !!--++    integer,                      intent(in)     :: Key
    !!--++    character(len=*),             intent(in)     :: Dire
    !!--++    integer,                      intent(in)     :: Na
    !!--++    integer,                      intent(in)     :: Nb
    !!--++    real(kind=cp),                intent(in)     :: Xl
    !!--++    real(kind=cp),                intent(in)     :: Xu
    !!--++    real(kind=cp),                intent(in)     :: Xs
    !!--++    integer,                      intent(in)     :: Ic
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    type(space_group_type),       intent(in)     :: Spg
    !!--++
    !!--++ Overloaded
    !!--++ Write on Vectors the Information for Molecule_Type
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Fill_RefCodes_Molec(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Molec,Spg)
       !---- Arguments ----!
       integer,                      intent(in)     :: Key
       character(len=*),             intent(in)     :: Dire
       integer,                      intent(in)     :: Na
       integer,                      intent(in)     :: Nb
       real(kind=cp),                intent(in)     :: Xl
       real(kind=cp),                intent(in)     :: Xu
       real(kind=cp),                intent(in)     :: Xs
       integer,                      intent(in)     :: Ic
       type(molecule_type),          intent(in out) :: Molec
       type(space_group_type),       intent(in)     :: Spg

       !---- Local variables ----!
       integer :: j, nc

       call init_err_refcodes()
       if (Na <= 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Number of atom no defined"
          return
       end if

       select case (dire)
          !---- FIX Directive ----!
          case ("fix")

             select case (key)
                case (0)

                   !---- nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined"
                         return

                      case ( 1:3)
                         !---- X_, Y_, Z_----!
                         if (molec%lI_coor(nb,na) /=0) then
                            nc=molec%lI_coor(nb,na)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (molec%lbiso(na) /=0) then
                            nc=molec%lbiso(na)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (molec%locc(na) /=0) then
                            nc=molec%locc(na)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case ( 6:12)
                         !---- B11_, ..., B23_ ----!
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined"
                         return

                      case (13:15)
                         !---- Xc_, Yc_, Zc_ ----!
                         if (molec%lxcentre(nb-12) /=0) then
                            nc=molec%lxcentre(nb-12)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case (16:18)
                         !---- Theta_, Phi_, Chi_ ----!
                         if (molec%lOrient(nb-15) /=0) then
                            nc=molec%lOrient(nb-15)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case (19:21)
                         !!! Not yet implement !!!

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   do j=1,3
                      if (molec%lI_coor(j,na) /=0) then
                         nc=molec%lI_coor(j,na)
                         call Delete_RefCodes(nc,molec)
                      end if
                   end do

                case (2)
                   !---- OCC ----!
                   if (molec%locc(na) /=0) then
                      nc=molec%locc(na)
                      call Delete_RefCodes(nc,molec)
                   end if

                case (3)
                   !---- BIS ----!
                   if (molec%lbiso(na) /=0) then
                      nc=molec%lbiso(na)
                      call Delete_RefCodes(nc,molec)
                   end if

                case (4)
                   !---- BAN ----!
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option not defined"
                   return

                case (5)
                  !---- ALL ----!
                  do j=1,3
                     if (molec%lI_coor(j,na) /=0) then
                        nc=molec%lI_coor(j,na)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do
                  if (molec%locc(na) /=0) then
                     nc=molec%locc(na)
                     call Delete_RefCodes(nc,molec)
                  end if
                  if (molec%lbiso(na) /=0) then
                     nc=molec%lbiso(na)
                     call Delete_RefCodes(nc,molec)
                  end if
                  do j=1,3
                     if (molec%lxcentre(j) /=0) then
                        nc=molec%lxcentre(j)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do
                  do j=1,3
                     if (molec%lOrient(j) /=0) then
                        nc=molec%lOrient(j)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do
                  !!! Falta THermal TLS !!!

               case (6)
                  !---- CEN ----!
                  do j=1,3
                     if (molec%lxcentre(j) /=0) then
                        nc=molec%lxcentre(j)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do

               case (7)
                  !---- ORI ----!
                  do j=1,3
                     if (molec%lOrient(j) /=0) then
                        nc=molec%lOrient(j)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do

               case (8)
                  !---- THE ----!
                  !!! Not Yet Implemented !!!
             end select

          !---- VARY Directive ----!
          case ("var")

             select case (key)
                case (0)

                   !---- Nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined"
                         return

                      case ( 1:3)
                         !--- X_, Y_, Z_ ----!
                         if (molec%lI_coor(nb,na) ==0) then
                            molec%mI_coor(nb,na)=1.0
                            call get_atompos_ctr(molec%I_Coor(:,na),  &
                                                 Spg, np_refi,  &
                                                 molec%lI_coor(:,na), &
                                                 molec%mI_coor(:,na))
                            if (molec%lI_coor(nb,na) == np_refi) then
                               V_Vec(np_refi)=molec%I_Coor(nb,na)
                               V_Name(np_refi)=trim(code_nam(nb))//trim(molec%AtName(na))
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (molec%lbiso(na) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molec%biso(na)
                            V_Name(np_refi)=trim(code_nam(nb))//trim(molec%AtName(na))
                            molec%mbiso(na)=1.0
                            molec%lbiso(na)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (molec%locc(na) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molec%occ(na)
                            V_Name(np_refi)=trim(code_nam(nb))//trim(molec%AtName(na))
                            molec%mocc(na)=1.0
                            molec%locc(na)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case ( 6:12)
                         !---- B11_, ..., B23_ ----!
                         err_refcodes=.true.
                         ERR_RefCodes_Mess="Option not defined"
                         return

                      case (13:15)
                         !---- Xc_, Yc_, Zc_ ----!
                         if (molec%lxcentre(nb-12) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molec%xcentre(nb-12)
                            V_Name(np_refi)=trim(code_nam(nb))//"Mol"
                            molec%mxcentre(nb-12)=1.0
                            molec%lxcentre(nb-12)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=0
                         end if

                      case (16:18)
                         !---- Theta_, Phi_, Chi_ ----!
                         if (molec%lOrient(nb-15) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molec%orient(nb-15)
                            V_Name(np_refi)=trim(code_nam(nb))//"Mol"
                            molec%mOrient(nb-15)=1.0
                            molec%lOrient(nb-15)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=0
                            V_list(np_refi)=0
                         end if

                      case (19:21)
                         !!! Not yet implement !!!

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   do j=1,3
                      if (molec%lI_coor(j,na) ==0) then
                         molec%mI_coor(j,na)=1.0
                         call get_atompos_ctr(molec%I_Coor(:,na),  &
                                              Spg, np_refi,  &
                                              molec%lI_coor(:,na), &
                                              molec%mI_coor(:,na))
                         if (molec%lI_coor(j,na) == np_refi) then
                            V_Vec(np_refi)=molec%I_Coor(j,na)
                            V_Name(np_refi)=trim(code_nam(j))//trim(molec%AtName(na))
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case (2)
                   !---- OCC ----!
                   if (molec%locc(na) ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=molec%occ(na)
                      V_Name(np_refi)=trim(code_nam(5))//trim(molec%AtName(na))
                      molec%mocc(na)=1.0
                      molec%locc(na)=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if

                case (3)
                   !---- BIS ----!
                   if (molec%lbiso(na) ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=molec%biso(na)
                      V_Name(np_refi)=trim(code_nam(4))//trim(molec%AtName(na))
                      molec%mbiso(na)=1.0
                      molec%lbiso(na)=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if

                case (4)
                   !---- BAN ----!
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option not defined"
                   return

                case (5)
                   !---- ALL ----!
                   do j=1,3
                      if (molec%lI_coor(j,na) ==0) then
                         molec%mI_coor(j,na)=1.0
                         call get_atompos_ctr(molec%I_Coor(:,na),  &
                                              Spg, np_refi,  &
                                              molec%lI_coor(:,na), &
                                              molec%mI_coor(:,na))
                         if (molec%lI_coor(j,na) == np_refi) then
                            V_Vec(np_refi)=molec%I_Coor(j,na)
                            V_Name(np_refi)=trim(code_nam(j))//trim(molec%AtName(na))
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do
                   if (molec%locc(na) ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=molec%occ(na)
                      V_Name(np_refi)=trim(code_nam(5))//trim(molec%AtName(na))
                      molec%mocc(na)=1.0
                      molec%locc(na)=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if
                   if (molec%lbiso(na) ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=molec%biso(na)
                      V_Name(np_refi)=trim(code_nam(4))//trim(molec%AtName(na))
                      molec%mbiso(na)=1.0
                      molec%lbiso(na)=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if
                   do j=1,3
                      if (molec%lxcentre(j) ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molec%xcentre(j)
                         V_name(np_refi)=trim(code_nam(12+j))//"entre_Mol"
                         molec%mxcentre(j)=1.0
                         molec%lxcentre(j)=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=0
                      end if
                   end do
                   do j=1,3
                      if (molec%lOrient(j) ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molec%orient(j)
                         V_name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"
                         molec%mOrient(j)=1.0
                         molec%lOrient(j)=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=0
                         V_list(np_refi)=0
                      end if
                   end do

                   !!! Falta THE !!!

                case (6)
                   !---- CEN ----!
                   do j=1,3
                      if (molec%lxcentre(j) ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molec%xcentre(j)
                         V_name(np_refi)=trim(code_nam(12+j))//"Mol"
                         molec%mxcentre(j)=1.0
                         molec%lxcentre(j)=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=0
                      end if
                   end do

                case (7)
                   !---- ORI ----!
                   do j=1,3
                      if (molec%lOrient(j) ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molec%orient(j)
                         V_name(np_refi)=trim(code_nam(15+j))//"Mol"
                         molec%mOrient(j)=1.0
                         molec%lOrient(j)=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=0
                      end if
                   end do

                case (8)
                   !---- THE ----!

                   !!! Not yet implemented !!!
             end select
       end select

       return
    End Subroutine Fill_RefCodes_Molec

    !!--++
    !!--++ Subroutine Fill_RefCodes_Magdom(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Mag_dom)
    !!--++    integer,                       intent(in)     :: Key
    !!--++    character(len=*),              intent(in)     :: Dire
    !!--++    integer,                       intent(in)     :: Na
    !!--++    integer,                       intent(in)     :: Nb
    !!--++    real(kind=cp),                 intent(in)     :: Xl
    !!--++    real(kind=cp),                 intent(in)     :: Xu
    !!--++    real(kind=cp),                 intent(in)     :: Xs
    !!--++    integer,                       intent(in)     :: Ic
    !!--++    type(Magnetic_Domain_type),    intent(in out) :: Mag_dom
    !!--++
    !!--++ Write on Vectors the Information for Magnetic Domains
    !!--++ Created: February - 2012
    !!--++
    Subroutine Fill_RefCodes_Magdom(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Mag_dom)
       !---- Arguments ----!
       integer,                       intent(in)     :: Key
       character(len=*),              intent(in)     :: Dire
       integer,                       intent(in)     :: Na
       integer,                       intent(in)     :: Nb
       real(kind=cp),                 intent(in)     :: Xl
       real(kind=cp),                 intent(in)     :: Xu
       real(kind=cp),                 intent(in)     :: Xs
       integer,                       intent(in)     :: Ic
       type(Magnetic_Domain_type),    intent(in out) :: Mag_dom

       !---- Local variables ----!
       integer           :: nc

       call init_err_refcodes()
       if (Na <= 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Number of domains no defined"
          return
       end if

       select case (dire)
          !---- FIX Directive ----!
          case ("fix")

             select case (key)
                case (0)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option not defined"
                   return

                case (4)
                   !---- Magdomain ----!
                      if (Mag_Dom%Lpop(nb,na) /=0) then
                         nc=Mag_Dom%Lpop(nb,na)
                         call Delete_RefCodes(nc,Mag_dom)
                      end if

                case (2:3)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option not defined"
                   return

                case (5:)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option not defined"
                   return
             end select

          !---- VARY Directive ----!
          case ("var")

              select case (key)
                case (0)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option not defined"
                   return

                case (4)
                   !---- Magd ----!
                      if (Mag_Dom%Lpop(nb,na) ==0) then
                            np_refi=np_refi+1
                            Mag_Dom%Lpop(nb,na)=np_refi
                            Mag_Dom%Mpop(nb,na)=1.0

                         if (Mag_Dom%Lpop(nb,na) == np_refi) then
                            V_Vec(np_refi)=Mag_Dom%pop(nb,na)
                            V_Name(np_refi)=trim(Mag_dom%lab(nb,na))
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=0
                         else
                            np_refi=np_refi-1
                         end if
                      end if

                case (2:3)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option not defined"
                   return

                case (5:)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option not defined"
                   return
              end select
       end select

       return
    End Subroutine Fill_RefCodes_Magdom

    !!----
    !!---- Subroutine Fill_RefGCodes(Key,Dire,namp,Xl,Xu,Xs,Ic,model,sys,Iphas)
    !!----    integer,                             intent(in)     :: Key  !0=> nb as below,
    !!----    character(len=*),                    intent(in)     :: Dire !GVar of GFix
    !!----    character(len=*),                    intent(in)     :: namp !Name of the parameter to be refined or fixed
    !!----    real(kind=cp),                       intent(in)     :: Xl   !Lower bound of parameter
    !!----    real(kind=cp),                       intent(in)     :: Xu   !Upper bound of parameter
    !!----    real(kind=cp),                       intent(in)     :: Xs   !Step of parameter
    !!----    integer,                             intent(in)     :: Ic   !Boundary condition (0:fixed or 1:periodic)
    !!----    type(NonAtomic_Parameter_List_Type), intent(in out) :: model
    !!----    character(len=*), optional,          intent( in)    :: sys
    !!----    integer,          optional,          intent(in)     :: Iphas
    !!----
    !!---- Write on Vectors the Information for Non atomic parameters
    !!----  Key=0 -> Provide information on individual parameter for atom na (nb should be given)
    !!----  Key=1  BKG -> fix or vary all background parameters
    !!----  Key=2  CELL -> fix or vary cell parameters
    !!             For the constraints on cell parameters the optional argument sys should be
    !!----         provided. It contains the crystal system and the (in the case of Monoclinic)
    !!----         the setting "a", "b" or "c" for the twofold axis. For instance the content
    !!----         of Sys may be Sys="Monoclinic c". Between the crystal system and the information
    !!----         about the setting a space must exist.
    !!----         The proper application of the constraints supposes that the cell parameters
    !!----         are contiguous in the list of non-atomic parameters.
    !!----  Key=3  UVW -> fix or vary u,v,w parameters
    !!----  Key=4  ASIZE -> fix or vary size parameters
    !!----  Key=5  ASTRAIN -> fix or vary strain parameters
    !!----  Key=6  EXTINCT -> fix or vary extinction parameters
    !!----  Key=7  SCALEFS -> fix or vary scale factors
    !!----  Key=9  ALL -> fix or vary all parameters
    !!----
    !!----
    !!---- Updated: November 3 - 2013
    !!
    Subroutine Fill_RefGCodes(Key,Dire,Namp,Xl,Xu,Xs,Ic,model,sys,Iphas)
       integer,                             intent(in)     :: Key  !0 => nb as below,
       character(len=*),                    intent(in)     :: Dire !GVar of GFix
       character(len=*),                    intent(in)     :: Namp !Name of the parameter to be refined or fixed
       real(kind=cp),                       intent(in)     :: Xl   !Lower bound of parameter
       real(kind=cp),                       intent(in)     :: Xu   !Upper bound of parameter
       real(kind=cp),                       intent(in)     :: Xs   !Step of parameter
       integer,                             intent(in)     :: Ic   !Boundary condition (0:fixed or 1:periodic)
       !integer,                             intent(in)     :: Imstart   !Starting index for mode search
       type(NonAtomic_Parameter_List_Type), intent(in out) :: model
       character(len=*), optional,          intent(in)     :: sys
       integer,          optional,          intent(in)     :: Iphas

       !---- Local variables ----!
       integer                  :: i,k,nc !j,,np_ini
       character(len=len(namp)) :: name_par
       character(len=15)        :: c_system
       character(len=5)         :: info
       character(len=2)         :: phase

       call init_err_refcodes()
       if (len_trim(namp) == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="No name for a parameter to be refined"
          return
       end if
       if(present(Iphas)) then
         write(unit=phase,fmt="(i2.2)") Iphas
       else
         phase="  "
       end if

       select case (trim(dire))
          !---- GFIX Directive ----!
          case ("gfix")

             select case (key)
                case (0) !Key=0 -> Provide information on individual parameter for type of parameter na (nb should be given)

                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par /= l_case(namp)) cycle
                      if ( model%par(i)%lcode /= 0) then
                         nc=model%par(i)%lcode
                         call Delete_RefGCodes(nc,model)
                      end if
                      exit
                   end do

                case (1) ! Key=1  BKG -> Fix all background parameters
                   !---- BKG ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:3) /= "bkg") cycle
                      if ( model%par(i)%lcode /= 0) then
                         nc=model%par(i)%lcode
                         call Delete_RefGCodes(nc,model)
                      end if
                      exit
                   end do

                case (2) !  Key=2  CELL -> Fix all cell parameters
                   !---- CELL ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if(index(name_par,"cell") /= 0 ) then
                         if(present(Iphas)) then
                           if(index(name_par,phase) == 0) cycle
                         end if
                         if ( model%par(i)%lcode /= 0) then
                            nc=model%par(i)%lcode
                            call Delete_RefGCodes(nc,model)
                         end if
                      end if
                   end do

                case (3) !  Key=3  UVW -> Fix U,V,W parameters
                   !---- UVW ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( trim(name_par) == "up" .or. trim(name_par) == "vp" .or. trim(name_par) == "wp" ) then
                         if ( model%par(i)%lcode /= 0) then
                            nc=model%par(i)%lcode
                            call Delete_RefGCodes(nc,model)
                         end if
                      end if
                   end do

                case (4) !  Key=4  ASIZE -> Fix all size parameters
                  !---- ASIZE ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:3) == "siz" .or. name_par(1:4) == "gsiz") then
                         if ( model%par(i)%lcode /= 0) then
                            nc=model%par(i)%lcode
                            call Delete_RefGCodes(nc,model)
                         end if
                      end if
                   end do

                case (5) !  Key=5  ASTRAIN -> Fix all strain parameters
                   !---- ASTRAIN ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:3) == "str" .or. name_par(1:4) == "lstr") then
                         if ( model%par(i)%lcode /= 0) then
                            nc=model%par(i)%lcode
                            call Delete_RefGCodes(nc,model)
                         end if
                      end if
                   end do

                case (6)  !  Key=6  EXTINCT -> Fix all extinction parameters
                   !---- EXTINCT ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:3) == "ext" .or. name_par(1:5) == "bcext") then
                         if ( model%par(i)%lcode /= 0) then
                            nc=model%par(i)%lcode
                            call Delete_RefGCodes(nc,model)
                         end if
                      end if
                   end do

                case (7)  !  Key=7  SCALEFS -> Fix all scale factors

                    !---- SCALEFS ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:2) == "sc" .or. name_par(1:5) == "scale") then
                         if ( model%par(i)%lcode /= 0) then
                            nc=model%par(i)%lcode
                            call Delete_RefGCodes(nc,model)
                         end if
                      end if
                   end do

                case (8)  !  Key=8  ALL -> Fix all non atomic parameters
                    !---- ALL ----!
                   do i=1,model%npar
                      if ( model%par(i)%lcode /= 0) then
                         nc=model%par(i)%lcode
                         call Delete_RefGCodes(nc,model)
                      end if
                   end do

                case (9:)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Incompatible information for this Non-atomic parameter "//trim(Namp)
                   return

             end select

          !---- GVARY Directive ----!
          case ("gvar")

             select case (key)
                case (0) !Key=0 -> Provide information on individual model parameter
                         !Vary a single model parameter
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par /= l_case(namp)) cycle
                      if ( model%par(i)%lcode == 0) then
                         call update_vect(i)
                         exit
                      end if
                   end do

                case (1) ! Key=1  BKG -> Vary all background parameters

                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:3) /= "bkg") cycle
                      if ( model%par(i)%lcode == 0) call update_vect(i)
                   end do

                case (2) !  Key=2  CELL -> Vary all cell parameters
                   !---- CELL ----!
                   if(present(sys)) then

                      do i=1,model%npar
                         name_par=l_case(model%par(i)%nam)
                         if ( index(name_par,"cell") /= 0) then
                            if(present(Iphas)) then
                              if(index(name_par,phase) == 0) cycle
                            end if
                            k=i
                            exit
                         end if
                      end do
                      i=index(sys," ")
                      c_system=sys(1:i-1)
                      info=sys(i+1:)
                      Select Case(c_system)

                          case("Triclinic")
                             do i=k,k+5
                               call update_vect(i)
                             end do

                          case("Monoclinic")

                             do i=k,k+2
                               call update_vect(i)
                             end do

                            if(index(Info,"b") /= 0) then
                               model%par(k+3)%value=90.0;  model%par(k+5)%value=90.0
                               call update_vect(k+4)

                            else if(index(Info,"c") /= 0) then
                               model%par(k+3)%value=90.0;  model%par(k+4)%value=90.0
                               call update_vect(k+5)

                            else if(index(Info,"a") /= 0) then
                               model%par(k+4)%value=90.0;  model%par(k+5)%value=90.0
                               call update_vect(k+3)

                            end if

                          case("Orthorhombic")
                            model%par(k+3:k+5)%value=90.0
                            do i=k,k+2
                              call update_vect(i)
                            end do

                          case("Tetragonal")
                            model%par(k+1)%value=model%par(k)%value
                            model%par(k+3:k+5)%value=90.0
                            call update_vect(k)
                            call update_vect(k+1,.false.)
                            call update_vect(k+2)

                          case("Trigonal","Hexagonal")
                            model%par(k+1)%value=model%par(k)%value
                            model%par(k+3:k+4)%value=90.0
                            model%par(k+5)%value=120.0
                            call update_vect(k)
                            call update_vect(k+1,.false.)
                            call update_vect(k+2)

                          case("Cubic")
                            model%par(k+3:k+5)%value=90.0
                            model%par(k+1)%value=model%par(k)%value
                            model%par(k+2)%value=model%par(k)%value
                            call update_vect(k)
                            call update_vect(k+1,.false.)
                            call update_vect(k+2,.false.)
                     End Select

                   else

                     do i=1,model%npar
                        name_par=l_case(model%par(i)%nam)
                        if ( index(name_par,"cell") /= 0) then
                            if(present(Iphas)) then
                              if(index(name_par,phase) == 0) cycle
                            end if
                           if ( model%par(i)%lcode == 0) call update_vect(i)
                        end if
                     end do

                   end if

                case (3) !  Key=3  UVW -> Vary U,V,W parameters
                   !---- UVW ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( trim(name_par) == "up" .or. trim(name_par) == "vp" .or. trim(name_par) == "wp" ) then
                         if ( model%par(i)%lcode == 0) call update_vect(i)
                      end if
                   end do

                case (4) !  Key=4  ASIZE -> Vary all size parameters
                  !---- ASIZE ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:3) == "siz" .or. name_par(1:4) == "gsiz") then
                         if ( model%par(i)%lcode == 0) call update_vect(i)
                      end if
                   end do

                case (5) !  Key=5  ASTRAIN -> Vary all strain parameters
                   !---- ASTRAIN ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:3) == "str" .or. name_par(1:4) == "lstr") then
                         if ( model%par(i)%lcode == 0) call update_vect(i)
                      end if
                   end do

                case (6)  !  Key=6  EXTINCT -> Vary all extinction parameters
                   !---- EXTINCT ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:3) == "ext") then
                         if ( model%par(i)%lcode == 0) call update_vect(i)
                      end if
                   end do

                case (7)  !  Key=7  SCALEFS -> Vary all scale factors
                   !---- SCALEFS ----!
                   do i=1,model%npar
                      name_par=l_case(model%par(i)%nam)
                      if ( name_par(1:2) == "sc" .or. name_par(1:5) == "scale") then
                         if ( model%par(i)%lcode == 0) call update_vect(i)
                      end if
                   end do

                case (8)  !  Key=8  ALL -> Vary all non atomic parameters
                  !---- ALL ----!
                   do i=1,model%npar
                      if ( model%par(i)%lcode == 0) call update_vect(i)
                   end do


                case (9:)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Incompatible information for this Non-atomic parameter "//trim(Namp)
                   return
             end select

       end select

       return

       contains

         Subroutine update_vect(n,up_np_refi)  !Internal subroutine to avoid copying the same text
           integer, intent(in):: n             !every time we need to do the same thing!
           logical, optional,intent(in):: up_np_refi
           logical :: local_up

           local_up=.true.
           if(present(up_np_refi)) local_up=up_np_refi
           if(local_up) np_refi=np_refi+1        !Here V_Vec_std should not be updated
           model%par(n)%lcode=np_refi
           model%par(n)%multip=1.0
           if(local_up) then
              V_Vec(np_refi)=model%par(n)%value
              V_Name(np_refi)=model%par(n)%nam
              V_Bounds(1,np_refi)=xl
              V_Bounds(2,np_refi)=xu
              V_Bounds(3,np_refi)=xs
              V_BCon(np_refi)=ic
              V_list(np_refi)=n
           end if
           return
         End Subroutine update_vect

    End Subroutine Fill_RefGCodes


    !!--++
    !!--++  Subroutine Get_Atombet_Ctr(X,Betas,Spgr,Codini,Icodes,Multip,Ord,Ss,Ipr)
    !!--++     real(kind=cp), dimension(3),             intent(in    ) :: X         !Atom position (fractional coordinates)
    !!--++     real(kind=cp), dimension(6),             intent(in out) :: Betas     !Anisotropic temperature factors
    !!--++     type(Space_Group_type),                  intent(in    ) :: Spgr      !Space Group
    !!--++     Integer,                                 intent(in out) :: Codini    !Last attributed parameter
    !!--++     Integer, dimension(6),                   intent(in out) :: Icodes    !codewords for betas only number
    !!--++     real(kind=cp), dimension(6),             intent(in out) :: Multip    !Multipliers
    !!--++     integer,                       optional, intent(in    ) :: Ord       !Order of the stabilizer
    !!--++     integer, dimension(:),         optional, intent(in    ) :: Ss        !Pointer to SymmOp. of stabilizer
    !!--++     integer,                       optional, intent(in    ) :: Ipr       !Printing unit for debug
    !!--++
    !!--++  Subroutine to get the appropriate constraints in the refinement codes of
    !!--++  anisotropic atomic displacement(thermal) parameters.
    !!--++  New algorithm based in the Wigner theorem.
    !!--++  The matrix Bet = Sum { R Beta RT} displays the symmetry constraints to be
    !!--++  applied to the anisotropic temperature factors. The sum runs over all rotational
    !!--++  symmetry operators of the stabilizer of the particular atom position in the given
    !!--++  space group.
    !!--++  There are a total of 29 kind of relations that may appear in the Bet matrix:
    !!--++
    !!--++     1    A A A 0   0   0  -> m-3m, -43m, 432, m-3,23, 3[111].2[001]
    !!--++     2    A A C 0   0   0  -> 4/mmm, -42m, 4mm, 422, 4/m, -4,4, 4[001]
    !!--++     3    A B A 0   0   0  -> 4[010]
    !!--++     4    A B B 0   0   0  -> 4[100]
    !!--++     5    A A A D   D   D  -> -3m, 3m, 32, -3, 3   3[111]
    !!--++     6    A A A D  -D  -D  -> 3[11-1]
    !!--++     7    A A A D  -D   D  -> 3[1-11]
    !!--++     8    A A A D   D  -D  -> 3[-111]
    !!--++     9    A A C A/2 0   0  -> 6/mmm, -6m2, 6mm, 622, 6/m, 6,-6,-3m, 32,-3, 3:  h 3[001]
    !!--++    10    A B C 0   0   0  -> mmm, mm2, 222  2[001] 2[100]
    !!--++    11    A A C D   0   0  -> 2[001], 2[110]    w
    !!--++    12    A B A 0   E   0  -> 2[010], 2[101]
    !!--++    13    A B B 0   0   F  -> 2[100], 2[011]
    !!--++    14    A B C B/2 0   0  -> 2[001], 2[100]    h
    !!--++    15    A B C A/2 0   0  -> 2[001], 2[010]    h
    !!--++    16    A B C D   0   0  -> 2/m, m, 2: 2[001] w
    !!--++    17    A B C 0   E   0  -> 2[010]
    !!--++    18    A B C 0   0   F  -> 2[100]
    !!--++    19    A A C D   E  -E  -> 2[110]            w
    !!--++    20    A A C D   E   E  -> 2[1-10]           w
    !!--++    21    A B A D   E  -D  -> 2[101]
    !!--++    22    A B A D   E   D  -> 2[10-1]
    !!--++    23    A B B D  -D   F  -> 2[011]
    !!--++    24    A B B D   D   F  -> 2[01-1]
    !!--++    25    A B C B/2 F/2 F  -> 2[100]            h
    !!--++    26    A B C A/2 0   F  -> 2[210]            h
    !!--++    27    A B C B/2 E   0  -> 2[120]            h
    !!--++    28    A B C A/2 E   E/2-> 2[010]            h
    !!--++    29    A B C D   E   F  -> 1, -1
    !!--++
    !!--++   Updated: 14 July 2011
    !!--++
    !!

    Subroutine Get_Atombet_Ctr(X,Betas,Spgr,Codini,Icodes,Multip,Ord,Ss,Ipr)
       real(kind=cp), dimension(3),             intent(in    ) :: X
       real(kind=cp), dimension(6),             intent(in out) :: Betas
       type(Space_Group_type),                  intent(in    ) :: Spgr
       Integer,                                 intent(in out) :: Codini
       Integer, dimension(6),                   intent(in out) :: Icodes
       real(kind=cp), dimension(6),             intent(in out) :: Multip
       integer,                       optional, intent(in    ) :: Ord
       integer, dimension(:),         optional, intent(in    ) :: Ss
       integer,                       optional, intent(in    ) :: Ipr


       ! Local variables
       character (len=1), dimension(6)   :: cdd
       integer                           :: i,j,order
       integer,           dimension(48)  :: ss_ptr
       integer,           dimension(6)   :: codd
       integer,           dimension(3,3) :: Rsym
       real(kind=cp),     dimension(3,3) :: bet,bett,Rs
       real(kind=cp),     dimension(6)   :: cod
       real(kind=cp),     dimension(3,48):: atr
       real(kind=cp),     parameter      :: epss=0.01_cp

       cod=real(icodes)

       do j=1,6
          if (cod(j) < 1.0 .and. abs(multip(j)) > epss)  then
             codini=codini+1
             cod(j) = real(codini)
          end if
       end do

       if (present(ord) .and. present(ss)) then
          order=ord
          ss_ptr(1:order) = ss(1:ord)
       else
          call get_stabilizer(x,Spgr,order,ss_ptr,atr)
       end if

       bet=reshape((/17.0, 7.0,3.0,  &
                     7.0,13.0,5.0,  &
                     3.0, 5.0,11.0/),(/3,3/))
       bett=bet
       if (order > 1 ) then
          do j=2,order
             Rsym=Spgr%SymOp(ss_ptr(j))%Rot
             Rs=real(Rsym)
             bett=bett+ matmul(Rs,matmul(bet,transpose(Rs)))
          end do
       end if
       Rsym=nint(1000.0*bett)
       codd=(/Rsym(1,1),Rsym(2,2),Rsym(3,3),Rsym(1,2),Rsym(1,3),Rsym(2,3)/)
       cdd=(/'a','b','c','d','e','f'/)
       multip=1.0
       !Search systematically all the possible constraints

       if(codd(1) == codd(2) .and. codd(1) == codd(3)) then ! a a a
         if(codd(4) == codd(5) .and. codd(4) == codd(6) ) then ! a a a d d d
           if(codd(4) == 0) then
             cdd=(/'a','a','a','0','0','0'/)     ! 1 A A A 0   0   0
             multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
             betas(4:6)=0.0
             betas(2:3)=betas(1)
             cod(2:3)=cod(1); cod(4:6)=0.0
           else
             cdd=(/'a','a','a','d','d','d'/)     ! 5 A A A D   D   D
             multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
             betas(5:6)=betas(4)
             betas(2:3)=betas(1)
             cod(2:3)=cod(1); cod(5:6)=cod(4)
           end if
         else if(codd(4) == -codd(5) .and. codd(4) == -codd(6) ) then !a a a d -d -d
           cdd=(/'a','a','a','d','d','d'/)       ! 6 A A A D  -D  -D
           multip=(/1.0,1.0,1.0,1.0,-1.0,-1.0/)
           betas(5:6)=-betas(4)
           betas(2:3)=betas(1)
           cod(2:3)=cod(1); cod(5:6)=cod(4)
         else if(codd(4) == -codd(5) .and. codd(4) ==  codd(6) ) then !a a a d -d  d
           cdd=(/'a','a','a','d','d','d'/)       ! 7 A A A D  -D   D
           multip=(/1.0,1.0,1.0,1.0,-1.0, 1.0/)
           betas(5)=-betas(4); betas(6)=betas(4)
           betas(2:3)=betas(1)
           cod(2:3)=cod(1); cod(5:6)= cod(4)
         else if(codd(4) ==  codd(5) .and. codd(4) == -codd(6) ) then !a a a d  d -d
           cdd=(/'a','a','a','d','d','d'/)       ! 8 A A A D   D  -D
           multip=(/1.0,1.0,1.0,1.0, 1.0,-1.0/)
           betas(6)=-betas(4); betas(5)=betas(4)
           betas(2:3)=betas(1)
           cod(2:3)=cod(1); cod(5:6)= cod(4)
         end if

       else if(codd(1) == codd(2)) then ! a a c
         if(codd(4) == codd(5) .and. codd(4) == codd(6) .and. codd(4) == 0) then ! a a c 0 0 0
             cdd=(/'a','a','c','0','0','0'/)     ! 2 A A C 0   0   0
             multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
             betas(4:6)=0.0
             betas(2)=betas(1)
             cod(2)=cod(1); cod(4:6)= 0.0
         else if(codd(5) == codd(6) .and. codd(5) == 0) then ! a a c x 0 0
             if(codd(4) == codd(1)/2) then
               cdd=(/'a','a','c','a','0','0'/)     ! 9 A A C A/2 0   0
               multip=(/1.0,1.0,1.0,0.5,0.0,0.0/)
               betas(5:6)=0.0; betas(4)=betas(1)*0.5
               betas(2)=betas(1)
               cod(2)=cod(1); cod(4)= cod(1); cod(5:6)=0.0
             else
               cdd=(/'a','a','c','d','0','0'/)     !11 A A C D   0   0
               multip=(/1.0,1.0,1.0,1.0,0.0,0.0/)
               betas(5:6)=0.0
               betas(2)=betas(1)
               cod(2)=cod(1); cod(5:6)=0.0
             end if
         else
             if(codd(5) == codd(6)) then  ! a a c d e e
               cdd=(/'a','a','c','d','e','e'/)     !20 A A C D   E   E
               multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
               betas(6)=betas(5)
               betas(2)=betas(1)
               cod(2)=cod(1); cod(6)=cod(5)
             else if(codd(5) == -codd(6)) then  ! a a c d e -e
               cdd=(/'a','a','c','d','e','e'/)     !19 A A C D   E  -E
               multip=(/1.0,1.0,1.0,1.0,1.0,-1.0/)
               betas(6)=-betas(5)
               betas(2)=betas(1)
               cod(2)=cod(1); cod(6)=cod(5)
             end if
         end if

       else if(codd(1) == codd(3)) then ! a b a
         if(codd(4) == codd(6)) then    ! a b a d x d
           if(codd(4) == 0) then  ! a b a 0 x 0
             if(codd(5) == 0) then ! a b a 0 0 0
               cdd=(/'a','b','a','0','0','0'/)     ! 3 A B A 0   0   0
               multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
               betas(4:6)=0.0
               betas(3)=betas(1)
               cod(3)=cod(1); cod(4:6)=0.0
             else                  ! a b a 0 e 0
               cdd=(/'a','b','a','0','e','0'/)     !12 A B A 0   E   0
               multip=(/1.0,1.0,1.0,0.0,1.0,0.0/)
               betas(4)=0.0;  betas(6)=0.0
               betas(3)=betas(1)
               cod(3)=cod(1); cod(4)=0.0;  cod(6)=0.0
             end if
           else  !! a b a d e d
             cdd=(/'a','b','a','d','e','d'/)       !22 A B A D   E   D
             multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
             betas(6)=betas(4)
             betas(3)=betas(1)
             cod(3)=cod(1); cod(6)=cod(4)
          end if

         else if(codd(4) == -codd(6)) then ! a b a d e -d
           cdd=(/'a','b','a','d','e','d'/)         !21 A B A D   E  -D
           multip=(/1.0,1.0,1.0,1.0,1.0,-1.0/)
           betas(6)=-betas(4)
           betas(3)=betas(1)
           cod(3)=cod(1); cod(6)=cod(4)
         end if

       else if(codd(2) == codd(3)) then ! a b b
         if(codd(4) == codd(5)) then    ! a b b d d x
           if(codd(4) == 0) then  ! a b b 0 0 x
             if(codd(6) == 0) then ! a b b 0 0 0
               cdd=(/'a','b','b','0','0','0'/)     ! 4 A B B 0   0   0
               multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
               betas(4:6)=0.0
               betas(3)=betas(2)
               cod(3)=cod(2); cod(4:6)=0.0
             else                  ! a b b 0 0 f
               cdd=(/'a','b','b','0','0','f'/)     !13 A B B 0   0   F
               multip=(/1.0,1.0,1.0,0.0,0.0,1.0/)
               betas(4:5)=0.0
               betas(3)=betas(2)
               cod(3)=cod(2); cod(4:5)=0.0
             end if
           else  !! a b b d d f
             cdd=(/'a','b','b','d','d','f'/)       !24 A B B D   D   F
             multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
             betas(5)=betas(4)
             betas(3)=betas(2)
             cod(3)=cod(2); cod(5)=cod(4)
           end if
         else if(codd(4) == -codd(5)) then ! a b b d -d e
           cdd=(/'a','b','b','d','d','f'/)         !23 A B B D  -D   F
           multip=(/1.0,1.0,1.0,1.0,-1.0,1.0/)
           betas(5)=-betas(4)
           betas(3)=betas(2)
           cod(3)=cod(2); cod(5)=cod(4)
         end if

       else !Now a /= b /= c

         if(codd(4) == codd(5) .and. codd(4) == 0) then ! a b c 0 0 x
           if(codd(6) == 0) then ! a b c 0 0 0
             cdd=(/'a','b','c','0','0','0'/)          !10 A B C 0   0   0
             multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
             betas(4:6)=0.0
             cod(4:6)=0.0
           else
             cdd=(/'a','b','c','0','0','f'/)          !18 A B C 0   0   F
             multip=(/1.0,1.0,1.0,0.0,0.0,1.0/)
             betas(4:5)=0.0
             cod(4:5)=0.0
           end  if
         else if(codd(5) == codd(6) .and. codd(5) == 0) then  ! a b c x 0 0
           if(codd(4) == codd(1)/2) then ! a b c a/2 0 0
             cdd=(/'a','b','c','a','0','0'/)          !15 A B C A/2 0   0
             multip=(/1.0,1.0,1.0,0.5,0.0,0.0/)
             betas(5:6)=0.0; betas(4)=betas(1)*0.5
             cod(4)=cod(1); cod(5:6)=0.0
           else if(codd(4) == codd(2)/2) then    !a b c b/2 0 0
             cdd=(/'a','b','c','b','0','0'/)          !14 A B C B/2 0   0
             multip=(/1.0,1.0,1.0,0.5,0.0,0.0/)
             betas(5:6)=0.0; betas(4)=betas(2)*0.5
             cod(4)=cod(2); cod(5:6)=0.0
           else
             cdd=(/'a','b','c','d','0','0'/)          !16 A B C D   0   0
             multip=(/1.0,1.0,1.0,1.0,0.0,0.0/)
             betas(5:6)=0.0
             cod(5:6)=0.0
           end  if
         else if(codd(4) == codd(6) .and. codd(4) == 0) then !a b c 0 e 0
           cdd=(/'a','b','c','0','e','0'/)            !17 A B C 0   E   0
           multip=(/1.0,1.0,1.0,0.0,1.0,0.0/)
           betas(4)=0.0; betas(6)=0.0
           cod(4)=0.0; cod(6)=0.0
         else if(codd(4) == codd(1)/2 .and. codd(5) == 0) then !a b c a/2 0 f
           cdd=(/'a','b','c','a','0','f'/)            !26 A B C A/2 0   F
           multip=(/1.0,1.0,1.0,0.5,0.0,1.0/)
           betas(4)=betas(1)*0.5; betas(5)=0.0
           cod(4)=cod(1); cod(5)=0.0
         else if(codd(4) == codd(2)/2 .and. codd(6) == 0) then !a b c b/2 e 0
           cdd=(/'a','b','c','b','e','0'/)            !27 A B C B/2 E   0
           multip=(/1.0,1.0,1.0,0.5,1.0,0.0/)
           betas(4)=betas(2)*0.5; betas(6)=0.0
           cod(4)=cod(2); cod(6)=0.0
         else if(codd(4) == codd(2)/2 .and. codd(5) == codd(6)/2) then !a b c b/2 f/2 f
           cdd=(/'a','b','c','b','f','f'/)            !25 A B C B/2 F/2 F
           multip=(/1.0,1.0,1.0,0.5,0.5,1.0/)
           betas(4)=betas(2)*0.5; betas(5)=betas(6)*0.5
           cod(4)=cod(2); cod(5)=cod(6)
         else if(codd(4) == codd(1)/2 .and. codd(6) == codd(5)/2) then !a b c a/2 e e/2
           cdd=(/'a','b','c','a','e','e'/)            !28 A B C A/2 E   E/2
           multip=(/1.0,1.0,1.0,0.5,1.0,0.5/)
           betas(4)=betas(1)*0.5; betas(6)=betas(5)*0.5
           cod(4)=cod(1); cod(6)=cod(5)
         else
           cdd=(/'a','b','c','d','e','f'/)            !29 A B C D   E   F
           multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
         end if
       end if

       do j=1,6
          if (multip(j) < epss .or. cdd(j) == "0" ) then
             icodes(j) = 0
          else
             icodes(j) = nint(cod(j))
          end if
       end do

       if(present(Ipr)) then
         Write(Ipr,'(a,6i5)')           '     Codes on Betas       :  ',Icodes
         Write(Ipr,'(a,6(a,1x),6f7.3)') '     Codes and multipliers:  ',cdd,multip
         Write(Ipr,'(a)')               '     Beta_TOT matrix:  '
         Do I=1,3
          Write(Ipr,'(a,3f12.4)')       '                      ',bett(i,:)
         End Do
       end if
       return
    End Subroutine Get_Atombet_Ctr

    !!--++
    !!--++  Subroutine Get_Atompos_Ctr(X,Spgr,Codini,Icodes,Multip,Ord,Ss,Att,Ipr)
    !!--++     real(kind=cp), dimension(3),       intent(in    ) :: x      !Atom position (fractional coordinates)
    !!--++     type(Space_Group_type),            intent(in    ) :: Spgr   !Space Group
    !!--++     Integer,                           intent(in out) :: codini !Last attributed parameter
    !!--++     integer,       dimension(3),       intent(in out) :: Icodes
    !!--++     real(kind=cp), dimension(3),       intent(in out) :: Multip
    !!--++     integer,                 optional, intent(in    ) :: Ord
    !!--++     integer, dimension(:),   optional, intent(in    ) :: Ss
    !!--++     integer, dimension(:,:), optional, intent(in    ) :: Atr
    !!--++     integer,                 optional, intent(in    ) :: Ipr
    !!--++
    !!--++     (Private)
    !!--++     Subroutine to get the appropriate constraints in the refinement codes of
    !!--++     atoms positions. The algorithm is based in an analysis of the symbol generated
    !!--++     for the symmetry elements of the operators belonging to the stabilizer of the
    !!--++     atom position. This routine operates with splitted codes in the sense of
    !!--++     FullProf rules
    !!--++
    !!--++     Updated: July - 2011 (JRC, the old subroutine has been completely changed)
    !!
    Subroutine Get_Atompos_Ctr(X,Spgr,Codini,ICodes,Multip,Ord,Ss,Att,Ipr)
       real(kind=cp), dimension(3),            intent(in)     :: X
       type(Space_Group_type),                 intent(in)     :: Spgr
       Integer,                                intent(in out) :: Codini
       Integer,       dimension(3),            intent(in out) :: ICodes
       real(kind=cp), dimension(3),            intent(in out) :: Multip
       integer,                       optional,intent(in)     :: Ord
       integer, dimension(:),         optional,intent(in)     :: Ss
       real(kind=cp), dimension(:,:), optional,intent(in)     :: Att
       integer,                       optional,intent(in)     :: Ipr

       ! Local variables
       integer                          :: i,j,k,order,L,L1,L2,ipar,j1
       integer,          dimension(3,3) :: RSym
       integer,          dimension(48)  :: ss_ptr
       real(kind=cp),    dimension(3,48):: atr
       real(kind=cp),    dimension(3)   :: tr

       character(len=40)                :: symbol,tsymbol,sym_symb
       Character(len=10), dimension(3)  :: nsymb
       Character(len=3),  dimension(3)  :: ssymb
       real(kind=cp),     parameter     :: epss=0.001

       if(present(ord) .and. present(ss) .and. present(att)) then
         order=ord
         ss_ptr(1:order) = ss(1:ord)
         atr(:,1:order)  = att(:,1:ord)
       else
         call get_stabilizer(x,Spgr,order,ss_ptr,atr)
       end if

       !If codes were not assigned with explicit number
       !attribute numbers bigger than initial Codini

       do j=1,3
          if(Icodes(j) < 1  .and. abs(multip(j)) > epss)  then
               codini = codini+1
               Icodes(j) = codini
          end if
       end do

       ssymb=(/"  x","  y","  z"/)

       if(present(Ipr)) write(unit=Ipr,fmt='(/a,3f10.5)')  ' => Atom Position:',x

       if(order > 1 ) then  !A constraint in position must exist

          if(present(Ipr)) write(unit=Ipr,fmt='(a)')   ' => List of symmetry element of the stabilizer without identity:'

          do k=2,order
           symbol=" "
           Rsym=Spgr%SymOp(ss_ptr(k))%Rot
           tr=Spgr%SymOp(ss_ptr(k))%tr + atr(:,k)
           call Get_SymSymb(Rsym,tr,Sym_Symb)
           call symmetry_symbol(Sym_Symb,tsymbol)
             i=index(tsymbol,";")
             if(i /= 0) then
               symbol=tsymbol(1:i-1)
               call Read_Xsym(tsymbol(i+1:),1,Rsym,Tr,.false.)
               if(sum(abs(x-tr)) < epss) then
                  ssymb=(/"  0","  0","  0"/)
                  if(present(Ipr)) then
                    write(unit=Ipr,fmt="(a,i2,a,t20,a,t55,a,t90,4a)") "     Operator ",k,": ", &
                    trim(Sym_Symb),trim(tsymbol),"  ssymb:" ,(ssymb(j)//"  ",j=1,3)
                  end if
                  cycle
               end if
             else
               symbol=tsymbol
             end if
             ipar=index(symbol,")")              !Translation element appears before position
             L =index(symbol(ipar+1:)," ")+ipar  !Position of the first blank after translation
             L1=index(symbol(ipar+1:),",")+ipar  !Position of the first comma after translation
             L2=index(symbol(L1+1:),",")+L1      !Position of the second comma
             if(L1 == 0) L1=1
             if(L2 == 0) L2=1
             if(L  == 0) L=1

             !Construct a new symbol that estabish automatically the constraints
             nsymb = (/symbol(L+1:L1-1),symbol(L1+1:L2-1),symbol(L2+1:)/)

             do i=1,3
                 do j=1,10  !Delete unwanted symbols (keep only x,y,z,2 and -
                    if(nsymb(i)(j:j) == " ") cycle
                    if(nsymb(i)(j:j) /= "x" .and. nsymb(i)(j:j) /= "y" .and. &
                       nsymb(i)(j:j) /= "z" .and. nsymb(i)(j:j) /= "-" .and. &
                       nsymb(i)(j:j) /= "2" ) nsymb(i)(j:j)=" "
                 end do
                 if(len_trim(nsymb(i))  == 0 .or. (index(nsymb(i),"x") == 0 .and. &
                    index(nsymb(i),"y") == 0 .and. index(nsymb(i),"z") == 0  ) ) then
                    ssymb(i)="  0"
                    cycle
                 end if
                 !Now remove 2s on the right of x,y, or z
                 j1=index(nsymb(i),"2")
                 if( j1 /= 0) then
                    if(len_trim(nsymb(i)) == j1) nsymb(i)=nsymb(i)(1:j1-1)
                 end if
                 !Now remove -s on the right of x,y, or z
                 j1=index(nsymb(i),"-")
                 if( j1 /= 0) then
                    if(len_trim(nsymb(i)) == j1) nsymb(i)=nsymb(i)(1:j1-1)
                 end if
                 nsymb(i)= adjustl(nsymb(i))
             end do

             if(ssymb(1) /= "  0" .and. ssymb(1) /= "  a") then
                ssymb(1)= nsymb(1)
                ssymb(1)= adjustr(ssymb(1))
             end if

             if(ssymb(2) /= "  0" .and. ssymb(2) /= "  a" .and. ssymb(2) /= "  b" .and. &
                ssymb(2) /= " -a" .and. ssymb(2) /= " 2a"   ) then
                ssymb(2) = nsymb(2)
                ssymb(2) = adjustr(ssymb(2))
             end if

             if(ssymb(3) /= "  0" .and. ssymb(3) /= "  a" .and. ssymb(3) /= "  b" .and. &
                ssymb(3) /= "  c" .and. ssymb(3) /= " 2a" .and. ssymb(3) /= " 2b" .and. &
                ssymb(3) /= " -a" .and. ssymb(3) /= " -b") then
                ssymb(3) = nsymb(3)
                ssymb(3) = adjustr(ssymb(3))
             end if

             do i=1,3
                if(ssymb(i)(3:3) == "x")  ssymb(i)(3:3) = "a"
             end do
             do i=1,3
                if(ssymb(i)(3:3) == "y")  ssymb(i)(3:3) = "b"
             end do
             do i=1,3
                if(ssymb(i)(3:3) == "z")  ssymb(i)(3:3) = "c"
             end do
             if(present(Ipr)) then
                write(unit=Ipr,fmt="(a,i2,a,t20,a,t55,a,t90,4a)") "     Operator ",k,": ", &
                trim(Sym_Symb),trim(tsymbol),"  Ssymb:" ,(ssymb(j)//"  ",j=1,3)
             end if

          end do !do k=1,order  over operators of the stabilizer

       else
         ssymb=(/"  a","  b","  c"/)

       end if  !order > 1

       do i=1,3                  !Fixing codes
         if(ssymb(i)=="  0") then
           Icodes(i)=0
           multip(i)=0.0
         end if
       end do

       if(index(ssymb(1),"a") /= 0) then

         do i=2,3  !Fixing codes
           if(index(ssymb(i),"-a") /= 0) then
             Icodes(i)=Icodes(1)
             multip(i)=-multip(1)
           else if(index(ssymb(i),"a") /= 0) then
             Icodes(i)=Icodes(1)
             multip(i)=multip(1)

             if(index(ssymb(i),"2") /= 0) then
               multip(i)=2.0* multip(1)
             else if(index(ssymb(1),"2") /= 0) then
               multip(i)=0.5* multip(1)
             end if

           end if
         end do
       else  !the x-coordinate is fixed, analyse y and z
         if(index(ssymb(2),"b") /= 0 .and. index(ssymb(3),"b") /= 0) then
           Icodes(3)=Icodes(2)
           if(ssymb(2) == ssymb(3)) then
             multip(3)= multip(2)
           else if(ssymb(3) == " -b" .and. ssymb(2) == "  b") then
             multip(3)= -multip(2)
           else if(ssymb(3) == "  b" .and. ssymb(2) == " -b") then
             multip(3)= -multip(2)
           end if
         end if

       end if !if(index(ssymb(1),"a") /= 0)

       do j=1,3
         if(abs(multip(j)) < epss) then
           Icodes(j) = 0
         end if
       end do
       if(present(Ipr)) then
         write(unit=Ipr,fmt="(a,3i5)")    "     Codes positions: ",Icodes
         write(unit=Ipr,fmt="(a,3f5.1)")  "     Multipliers    : ",multip
         write(unit=Ipr,fmt="(5a)")       "     Codes   string : ( ",(ssymb(j),j=1,3) ," )"
       end if
       return
    End Subroutine Get_Atompos_Ctr

    !!--++
    !!--++ Subroutine Get_ConCodes_Line(Line,FAtom/FmAtom/MolCrys/Molec/MagDom)
    !!--++    character(len=*),             intent(in)     :: Line
    !!--++    integer,                      intent(in)     :: Nat
    !!--++    type(Atom_List_Type),         intent(in out) :: FAtom
    !!--++    or
    !!--++    type(mAtom_List_Type),        intent(in out) :: FmAtom
    !!--++    or
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++    or
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    or
    !!--++    type(Magnetic_Domain_type),   intent(in out) :: Mag_Dom
    !!--++
    !!--++    (Private)
    !!--++    Get the Constraints relations
    !!--++
    !!--++ Update: March - 2005
    !!

    !!--++
    !!--++ Subroutine Get_ConCodes_Line_FAtom(Line,FAtom)
    !!--++    character(len=*),         intent(in)     :: Line
    !!--++    integer,                  intent(in)     :: Nat
    !!--++    type(Atom_List_Type),     intent(in out) :: FAtom
    !!--++
    !!--++ Overloaded
    !!--++ Get the Constraints relations
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_ConCodes_Line_FAtom(Line,FAtom)
       !---- Arguments ----!
       character(len=*),     intent(in)     :: Line
       type(Atom_List_Type), intent(in out) :: FAtom

       !---- Local variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: j,ic,n,na,nb,nc,nd,npos
       integer                          :: nl,nl2,iv
       integer, dimension(1)            :: ivet
       real(kind=cp)                    :: fac_0,fac_1
       real(kind=cp),dimension(1)       :: vet

       call init_err_refcodes()

       nl=0
       call getword(line,label,ic)
       if (ic < 2) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="EQUAL keyword needs two labels: "//trim(line)
          return
       end if

       !---- Set Futher Information----!
       !---- Na is the number of atom on List
       !---- Nb is the key (X,Y,Z,Occ,...)
       !---- Fac0 is the multiplier
       !---- Nl is the number of refinement parameter

       npos=index(label(1),"_")
       if (npos == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="The name "//trim(label(1))//" does not fit any known code-name of CrysFML "
          return
       end if

       na=0
       do j=1,FAtom%Natoms
          if (u_case(FAtom%atom(j)%lab) == u_case(label(1)(npos+1:npos+6))) then
             na=j
             exit
          end if
       end do
       if (na == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess=" Atom label not found for "//trim(line)
          return
       end if

       nb=0
       do j=1,ncode
          if (u_case(label(1)(1:npos))==u_case(trim(code_nam(j)))) then
             nb=j
             exit
          end if
       end do
       if (nb == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Code-name not found for parameter name: "//trim(label(1))
          return
       end if

       select case (nb)
          case ( 1:3)
             fac_0=FAtom%atom(na)%mx(nb)
                nl=FAtom%atom(na)%lx(nb)
          case ( 4)
             fac_0=FAtom%atom(na)%mbiso
                nl=FAtom%atom(na)%lbiso
          case ( 5)
             fac_0=FAtom%atom(na)%mocc
                nl=FAtom%atom(na)%locc
          case ( 6:11)
             fac_0=FAtom%atom(na)%mu(nb-5)
                nl=FAtom%atom(na)%lu(nb-5)
          case (12)
             fac_0=FAtom%atom(na)%mu(1)
          case (13:)
             err_refcodes=.true.
             ERR_RefCodes_Mess="Incompatible Code-name for parameter name: "//trim(label(1))
             return
       end select ! nb

       if (nb < ncode) then
          if (nl == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No refinable parameter was selected for "//trim(label(1))
             return
          end if
       end if

       !---- Set the rest elements in Constraints ----!
       n=1
       do
          n=n+1
          if (n > ic) exit

          npos=index(label(n),"_")
          if (npos ==0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No CrysFML code-name was found for "//trim(label(n))
             return
          end if

          nc=0
          do j=1,FAtom%Natoms
             if (u_case(FAtom%atom(j)%lab) == u_case(label(n)(npos+1:npos+6))) then
                nc=j
                exit
             end if
          end do
          if (nc == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess=" Atom label not found for "//trim(label(n))
             return
          end if

          nd=0
          do j=1,ncode
             if (u_case(label(n)(1:npos))==u_case(trim(code_nam(j)))) then
                nd=j
                exit
             end if
          end do
          if (nd == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Code-name not found for "//trim(label(n))
             return
          end if

          !---- Is there a new multiplier?
          n=n+1
          call getnum(label(n),vet,ivet,iv)
          if (iv == 1) then
             fac_1=vet(1)
          else
             fac_1=fac_0
             n=n-1
          end if

          select case (nd)
             case ( 1:3)
                nl2=FAtom%atom(nc)%lx(nd)
                call Delete_refCodes(nl2,FAtom)
                FAtom%atom(nc)%mx(nd)=fac_1
                FAtom%atom(nc)%lx(nd)=nl
             case ( 4)
                nl2=FAtom%atom(nc)%lbiso
                call Delete_refCodes(nl2,FAtom)
                FAtom%atom(nc)%mbiso=fac_1
                FAtom%atom(nc)%lbiso=nl
             case ( 5)
                nl2=FAtom%atom(nc)%locc
                call Delete_refCodes(nl2,FAtom)
                FAtom%atom(nc)%mocc=fac_1
                FAtom%atom(nc)%locc=nl
             case ( 6:11)
                nl2=FAtom%atom(nc)%lu(nd-5)
                call Delete_refCodes(nl2,FAtom)
                FAtom%atom(nc)%mu(nd-5)=fac_1
                FAtom%atom(nc)%lu(nd-5)=nl
             case (12)
                do j=1,6
                   nl2=FAtom%atom(nc)%lu(j)
                   call Delete_refCodes(nl2,FAtom)
                   FAtom%atom(nc)%mu(j)=fac_1
                   FAtom%atom(nc)%lu(j)=FAtom%atom(na)%lu(j)
                end do
                np_cons=np_cons+5
             case (13:)
                err_refcodes=.true.
                ERR_RefCodes_Mess="Incompatible Code-name for parameter name: "//trim(label(1))
                return
          end select ! nb

          np_cons=np_cons+1

       end do

       return
    End Subroutine Get_ConCodes_Line_FAtom

    !!--++
    !!--++ Subroutine Get_ConCodes_Line_FmAtom(Line,FmAtom)
    !!--++    character(len=*),         intent(in)     :: Line
    !!--++    integer,                  intent(in)     :: Nat
    !!--++    type(mAtom_List_Type),    intent(in out) :: FmAtom
    !!--++
    !!--++ Get the Magnetic Constraints Relations: Presently only 'equal'
    !!--++ mag clone of Get_ConCodes_Line_FAtom
    !!--++ Created: December - 2011
    !!
    Subroutine Get_ConCodes_Line_FmAtom(Line,FmAtom)
       !---- Arguments ----!
       character(len=*),     intent(in)     :: Line
       type(mAtom_List_Type),intent(in out) :: FmAtom

       !---- Local variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: j,ic,n,na,nb,nc,nd,npos
       integer                          :: nl,nl2,iv,ik
       integer, dimension(1)            :: ivet
       real(kind=cp)                    :: fac_0,fac_1
       real(kind=cp),dimension(1)       :: vet

       call init_err_refcodes()

       nl=0
       call getword(line,label,ic)
       if (ic < 2) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="EQUAL keyword needs two labels: "//trim(line)
          return
       end if

       !---- Set Futher Information----!
       !---- Na is the number of atom on List
       !---- Nb is the number of keys (Rx,Ry,Rz,Ix,Iy,Iz,MagPh...)
       !---- Fac0 is the multiplier
       !---- Nl is the number of refinement parameter

       npos=index(label(1),"_")
       if (npos ==0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="The name "//trim(label(1))//" does not fit any known code-name of CrysFML "
          return
       end if

       na=0
       do j=1,FmAtom%Natoms
          if (u_case(FmAtom%atom(j)%lab) == u_case(label(1)(npos+1:npos+6))) then
             na=j
             exit
          end if
       end do
       if (na == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess=" Atom label not found for "//trim(line)
          return
       end if

       nb=0
       do j=1,mNcode
          if (u_case(label(1)(1:npos))==u_case(trim(mcode_nam(j)))) then
             nb=j
             exit
          end if
       end do
       if (nb == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Code-name not found for parameter name: "//trim(label(1))
          return
       end if

       !---- Get im, number of the magnetic matrices/irrep set
       ik=FmAtom%atom(na)%nvk

       select case (nb)
          case ( 1:3)
             fac_0=FmAtom%atom(na)%mSkR(nb,ik)
                nl=FmAtom%atom(na)%lSkR(nb,ik)
          case ( 4:6)
             fac_0=FmAtom%atom(na)%mSkI(nb-3,ik)
                nl=FmAtom%atom(na)%lSkI(nb-3,ik)
          case ( 7:9)
             fac_0=FmAtom%atom(na)%mSkR(nb-6,ik)
                nl=FmAtom%atom(na)%lSkR(nb-6,ik)
          case (10:12)
             fac_0=FmAtom%atom(na)%mSkI(nb-9,ik)
                nl=FmAtom%atom(na)%lSkI(nb-9,ik)
          case (13)
             fac_0=FmAtom%atom(na)%mmphas(ik)
                nl=FmAtom%atom(na)%lmphas(ik)
           case (14:22)
             fac_0=FmAtom%atom(na)%mbas(nb-13,ik)
                nl=FmAtom%atom(na)%lbas(nb-13,ik)
         case (23:)
             err_refcodes=.true.
             ERR_RefCodes_Mess="Incompatible Code-name for parameter name: "//trim(label(1))
             return
       end select ! nb

       if (nb < mNcode) then
          if (nl == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No refinable parameter was selected for "//trim(label(1))
             return
          end if
       end if

       !---- Set the rest elements in Constraints ----!
       n=1
       do
          n=n+1
          if (n > ic) exit

          npos=index(label(n),"_")
          if (npos ==0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No CrysFML code-name was found for "//trim(label(n))
             return
          end if

          nc=0
          do j=1,FmAtom%Natoms
             if (u_case(FmAtom%atom(j)%lab) == u_case(label(n)(npos+1:npos+6))) then
                nc=j
                exit
             end if
          end do
          if (nc == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess=" Atom label not found for "//trim(label(n))
             return
          end if

          nd=0
          do j=1,ncode
             if (u_case(label(n)(1:npos))==u_case(trim(mcode_nam(j)))) then
                nd=j
                exit
             end if
          end do
          if (nd == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Code-name not found for "//trim(label(n))
             return
          end if

          !---- Is there a new multiplier?
          n=n+1
          call getnum(label(n),vet,ivet,iv)
          if (iv == 1) then
             fac_1=vet(1)
          else
             fac_1=fac_0
             n=n-1
          end if

          select case (nd)
             case ( 1:3)
                nl2=FmAtom%atom(nc)%lSkR(nd,ik)
                call Delete_refCodes(nl2,FmAtom)
                FmAtom%atom(nc)%mSkR(nd,ik)=fac_1
                FmAtom%atom(nc)%lSkR(nd,ik)=nl
             case ( 4:6)
                nl2=FmAtom%atom(nc)%lSkI(nd-3,ik)
                call Delete_refCodes(nl2,FmAtom)
                FmAtom%atom(nc)%mSkI(nd-3,ik)=fac_1
                FmAtom%atom(nc)%lSkI(nd-3,ik)=nl
             case ( 7:9)
                nl2=FmAtom%atom(nc)%lSkR(nd-6,ik)
                call Delete_refCodes(nl2,FmAtom)
                FmAtom%atom(nc)%mSkR(nd-6,ik)=fac_1
                FmAtom%atom(nc)%lSkR(nd-6,ik)=nl
             case (10:12)
                nl2=FmAtom%atom(nc)%lSkI(nd-9,ik)
                call Delete_refCodes(nl2,FmAtom)
                FmAtom%atom(nc)%mSkI(nd-9,ik)=fac_1
                FmAtom%atom(nc)%lSkI(nd-9,ik)=nl
             case (13)
                nl2=FmAtom%atom(nc)%lmphas(ik)
                call Delete_refCodes(nl2,FmAtom)
                FmAtom%atom(nc)%mmphas(ik)=fac_1
                FmAtom%atom(nc)%lmphas(ik)=nl
             case (14:22)
                nl2=FmAtom%atom(nc)%lbas(nd-13,ik)
                call Delete_refCodes(nl2,FmAtom)
                FmAtom%atom(nc)%mbas(nd-13,ik)=fac_1
                FmAtom%atom(nc)%lbas(nd-13,ik)=nl
             case (23:)
                err_refcodes=.true.
                ERR_RefCodes_Mess="Incompatible Code-name for parameter name: "//trim(label(1))
                return
          end select ! nb

          np_cons=np_cons+1

       end do

       return
    End Subroutine Get_ConCodes_Line_FmAtom

    !!--++
    !!--++ Subroutine Get_ConCodes_Line_Molcrys(Line,Molcrys)
    !!--++    character(len=*),             intent(in)     :: Line
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++
    !!--++ Overloaded
    !!--++ Get the Constraints relations
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_ConCodes_Line_Molcrys(Line,Molcrys)
       !---- Arguments ----!
       character(len=*),             intent(in)     :: Line
       type(molecular_Crystal_type), intent(in out) :: MolCrys

       !---- Loval variables ----!
       character(len=5)                 :: car
       character(len=20), dimension(30) :: label
       integer                          :: i,j,ic,n,na,naa,nb,nc,ncc,nd
       integer                          :: npos, nposm, nmol1,nmol2
       integer                          :: nl,nl2,iv
       integer, dimension(1)            :: ivet
       real(kind=cp)                    :: fac_0,fac_1
       real(kind=cp),dimension(1)       :: vet

       call init_err_refcodes()

       call getword(line,label,ic)
       if (ic < 2) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="EQUAL keyword needs two labels: "//trim(line)
          return
       end if

       !---- Set Father Information ----!
       !---- Na is the number of atom on List
       !---- Nb is the key (X,Y,Z,Occ,...)
       !---- Fac0 is the multiplier
       !---- Nl is the number of refinement parameter
       npos=index(label(1),"_")
       if (npos ==0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="The name "//trim(label(1))//" does not fit any known code-name of CrysFML "
          return
       end if
       nposm=index(u_case(label(1)),"MOL")
       if (nposm /= 0) then
          car=adjustl(label(1)(nposm+3:))
          if (car(1:2) == "  ") then
             nmol1=0
          else
             read(unit=car,fmt="(i2)") nmol1
          end if
       else
          nmol1=-1
       end if

       na=0
       do j=1,molcrys%n_free
          if (u_case(molcrys%atm(j)%lab) == u_case(label(1)(npos+1:npos+6))) then
             na=j
             exit
          end if
       end do
       if (na == 0) then
          do i=1,molcrys%n_mol
             do j=1,molcrys%mol(i)%natoms
                if (u_case(molcrys%mol(i)%Atname(j)) == u_case(label(1)(npos+1:npos+6))) then
                   if (j > 1) then
                      na=molcrys%n_free+sum(molcrys%mol(1:j-1)%natoms)+j
                   else
                      na=molcrys%n_free+j
                   end if
                   exit
                end if
             end do
          end do
       end if

       nb=0
       do j=1,ncode
          if (u_case(label(1)(1:npos))==u_case(trim(code_nam(j)))) then
             nb=j
             exit
          end if
       end do
       if (nb == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Code-name not found for parameter name: "//trim(label(1))
          return
       end if

       !---- Checking ----!
       if (nb <= 12 .and. na==0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Incompatible option for parameter name: "//trim(label(1))
          return
       end if

       nl=0
       select case (nb)
          case ( 1:3)
             !---- X_, Y_, Z_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mx(nb)
                   nl=molcrys%atm(na)%lx(nb)
             else
                naa=na-molcrys%n_free
                do i=1,molcrys%n_mol
                   if (naa > molcrys%mol(i)%natoms) then
                      naa=naa-molcrys%mol(i)%natoms
                      cycle
                   end if
                   fac_0=molcrys%mol(i)%mI_coor(nb,naa)
                      nl=molcrys%mol(i)%lI_coor(nb,naa)
                end do
             end if

          case ( 4)
             !---- Biso_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mbiso
                   nl=molcrys%atm(na)%lbiso
             else
                naa=na-molcrys%n_free
                do i=1,molcrys%n_mol
                   if (naa > molcrys%mol(i)%natoms) then
                      naa=naa-molcrys%mol(i)%natoms
                      cycle
                   end if
                   fac_0=molcrys%mol(i)%mbiso(naa)
                      nl=molcrys%mol(i)%lbiso(naa)
                end do
             end if

          case ( 5)
             !---- Occ_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mocc
                   nl=molcrys%atm(na)%locc
             else
                naa=na-molcrys%n_free
                do i=1,molcrys%n_mol
                   if (naa > molcrys%mol(i)%natoms) then
                      naa=naa-molcrys%mol(i)%natoms
                      cycle
                   end if
                   fac_0=molcrys%mol(i)%mocc(naa)
                      nl=molcrys%mol(i)%locc(naa)
                end do
             end if

          case ( 6:11)
             !---- B11_, ..., B23_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mu(nb-5)
                   nl=molcrys%atm(na)%lu(nb-5)
             else
                err_refcodes=.true.
                ERR_RefCodes_Mess="Option no valid"
                return
             end if

          case (12)
             !---- Banis_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mu(1)
             else
                err_refcodes=.true.
                ERR_RefCodes_Mess="Option no valid"
                return
             end if

          case (13:15)
             !---- Xc_, Yc_, Zc_ ----!
             select case (nmol1)
                case (-1)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option no valid"
                   return

                case (0)
                   fac_0=molcrys%mol(1)%mxcentre(nb-12)
                   nl=molcrys%mol(1)%lxcentre(nb-12)

                case (1:)
                   fac_0=molcrys%mol(nmol1)%mxcentre(nb-12)
                   nl=molcrys%mol(nmol1)%lxcentre(nb-12)
             end select

          case (16:18)
             !---- Theta_ , Phi_, Chi_ ----!
             select case (nmol1)
                case (-1)
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option no valid"
                   return

                case (0)
                   fac_0=molcrys%mol(1)%mOrient(nb-15)
                   nl=molcrys%mol(1)%lOrient(nb-15)

                case (1)
                   fac_0=molcrys%mol(nmol1)%mOrient(nb-15)
                   nl=molcrys%mol(nmol1)%lOrient(nb-15)
             end select

          case (19:21)
             !!! Not yet Implemented !!!
       end select ! nb

       if (nb < ncode) then
          if (nl == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No refinable parameter was selected for "//trim(label(1))
             return
          end if
       end if

       !---- Set the rest elements in Constraints ----!
       n=1
       do
          n=n+1
          if (n > ic) exit

          npos=index(label(n),"_")
          if (npos ==0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No CrysFML code-name was found for "//trim(label(n))
             return
          end if
          nposm=index(u_case(label(n)),"MOL")
          if (nposm /= 0) then
             car=adjustl(label(n)(nposm+3:))
             if (car(1:2) == "  ") then
                nmol2=0
             else
                read(unit=car,fmt="(i2)") nmol2
             end if
          else
             nmol2=-1
          end if

          nc=0
          do j=1,molcrys%n_free
             if (u_case(molcrys%atm(j)%lab) == u_case(label(n)(npos+1:npos+6))) then
                nc=j
                exit
             end if
          end do
          if (nc == 0) then
             do i=1,molcrys%n_mol
                do j=1,molcrys%mol(i)%natoms
                   if (u_case(molcrys%mol(i)%Atname(j)) == u_case(label(n)(npos+1:npos+6))) then
                      if (j > 1) then
                         nc=molcrys%n_free+sum(molcrys%mol(1:j-1)%natoms)+j
                      else
                         nc=molcrys%n_free+j
                      end if
                      exit
                   end if
                end do
             end do
          end if

          nd=0
          do j=1,ncode
             if (u_case(label(n)(1:npos))==u_case(trim(code_nam(j)))) then
                nd=j
                exit
             end if
          end do
          if (nd == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Code-name not found for "//trim(label(n))
             return
          end if

          !---- Checking ----!
          if (nd <= 12 .and. nc==0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Incompatible option for parameter name: "//trim(label(n))
             return
          end if

          !---- Is there a new multiplier? ----!
          n=n+1
          call getnum(label(n),vet,ivet,iv)
          if (iv == 1) then
             fac_1=vet(1)
          else
             fac_1=fac_0
             n=n-1
          end if

          select case (nd)
             case ( 1:3)
                !---- X_, Y_, Z_ ----!
                if (nc <= molcrys%n_free) then
                   nl2=molcrys%atm(nc)%lx(nd)
                   call Delete_RefCodes(nl2,molcrys)
                   molcrys%atm(nc)%mx(nd)=fac_1
                   molcrys%atm(nc)%lx(nd)=nl
                else
                   ncc=nc-molcrys%n_free
                   do i=1,molcrys%n_mol
                      if (ncc > molcrys%mol(i)%natoms) then
                         ncc=ncc-molcrys%mol(i)%natoms
                         cycle
                      end if
                      nl2=molcrys%mol(i)%lI_coor(nd,ncc)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(i)%mI_coor(nd,ncc)=fac_1
                      molcrys%mol(i)%lI_coor(nd,ncc)=nl
                   end do
                end if

             case ( 4)
                !---- Biso_ ----!
                if (nc <= molcrys%n_free) then
                   nl2=molcrys%atm(nc)%lbiso
                   call Delete_RefCodes(nl2,molcrys)
                   molcrys%atm(nc)%mbiso=fac_1
                   molcrys%atm(nc)%lbiso=nl
                else
                   ncc=nc-molcrys%n_free
                   do i=1,molcrys%n_mol
                      if (ncc > molcrys%mol(i)%natoms) then
                         ncc=ncc-molcrys%mol(i)%natoms
                         cycle
                      end if
                      nl2=molcrys%mol(i)%lbiso(ncc)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(i)%mbiso(ncc)=fac_1
                      molcrys%mol(i)%lbiso(ncc)=nl
                   end do
                end if

             case ( 5)
                !---- Occ_ ----!
                if (nc <= molcrys%n_free) then
                   nl2=molcrys%atm(nc)%locc
                   call Delete_RefCodes(nl2,molcrys)
                   molcrys%atm(nc)%mocc=fac_1
                   molcrys%atm(nc)%locc=nl
                else
                   ncc=nc-molcrys%n_free
                   do i=1,molcrys%n_mol
                      if (ncc > molcrys%mol(i)%natoms) then
                         ncc=ncc-molcrys%mol(i)%natoms
                         cycle
                      end if
                      nl2=molcrys%mol(i)%locc(ncc)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(i)%mocc(ncc)=fac_1
                      molcrys%mol(i)%locc(ncc)=nl
                   end do
                end if

             case ( 6:11)
                !---- B11_, ...., B23_ ----!
                if (nc <= molcrys%n_free) then
                   nl2=molcrys%atm(nc)%lu(nd-5)
                   call Delete_RefCodes(nl2,molcrys)
                   molcrys%atm(nc)%mu(nd-5)=fac_1
                   molcrys%atm(nc)%lu(nd-5)=nl
                else
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option no valid"
                   return
                end if

             case (12)
                !---- Banis_ ----!
                if (nc <= molcrys%n_free) then
                   do j=1,6
                      nl2=molcrys%atm(nc)%lu(j)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%atm(nc)%mu(j)=fac_1
                      molcrys%atm(nc)%lu(j)=molcrys%atm(na)%lu(j)
                   end do
                   np_cons=np_cons+5
                else
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Option no valid"
                   return
                end if

             case (13:15)
                select case (nmol2)
                   case (-1)
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Option no valid"
                      return

                   case (0)
                      do i=1,molcrys%n_mol
                         nl2=molcrys%mol(i)%lxcentre(nc-12)
                         call Delete_RefCodes(nl2,molcrys)
                         molcrys%mol(i)%mxcentre(nc-12)=fac_1
                         molcrys%mol(i)%lxcentre(nc-12)=nl
                      end do
                      np_cons=np_cons+(i-1)

                   case (1:)
                      nl2=molcrys%mol(nmol2)%lxcentre(nc-12)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(nmol2)%mxcentre(nc-12)=fac_1
                      molcrys%mol(nmol2)%lxcentre(nc-12)=nl
                end select

             case (16:18)
                select case (nmol2)
                   case (-1)
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Option no valid"
                      return

                   case (0)
                      do i=1,molcrys%n_mol
                         nl2=molcrys%mol(i)%lorient(nc-15)
                         call Delete_RefCodes(nl2,molcrys)
                         molcrys%mol(i)%morient(nc-15)=fac_1
                         molcrys%mol(i)%lorient(nc-15)=nl
                      end do
                      np_cons=np_cons+(i-1)

                   case (1:)
                      nl2=molcrys%mol(nmol2)%lorient(nc-15)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(nmol2)%morient(nc-15)=fac_1
                      molcrys%mol(nmol2)%lorient(nc-15)=nl
                end select

             case (19:21)
                !!! Not yet implemented !!!

          end select ! nb
          np_cons=np_cons+1

       end do

       return
    End Subroutine Get_ConCodes_Line_Molcrys

    !!--++
    !!--++ Subroutine Get_ConCodes_Line_Molec(Line,Molec)
    !!--++    character(len=*),    intent(in)     :: Line
    !!--++    type(molecule_type), intent(in out) :: Molec
    !!--++
    !!--++ Overloaded
    !!--++ Get the Constraints relations
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_ConCodes_Line_Molec(Line,Molec)
       !---- Arguments ----!
       character(len=*),    intent(in)     :: Line
       type(molecule_type), intent(in out) :: Molec

       !---- Loval variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: j,ic,na,nb,nc,nd,npos!,i,naa,ncc
       integer                          :: n,nl,nl2,iv
       integer, dimension(1)            :: ivet
       real(kind=cp)                    :: fac_0,fac_1
       real(kind=cp),dimension(1)       :: vet

       call init_err_refcodes()

       call getword(line,label,ic)
       if (ic < 2) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="EQUAL keyword needs two labels: "//trim(line)
          return
       end if

       !---- Set Father Information ----!
       !---- Na is the number of atom on List
       !---- Nb is the key (X,Y,Z,Occ,...)
       !---- Fac0 is the multiplier
       !---- Nl is the number of refinement parameter
       npos=index(label(1),"_")
       if (npos ==0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="The name "//trim(label(1))//" does not fit any known code-name of CrysFML "
          return
       end if

       na=0
       do j=1,molec%natoms
          if (u_case(molec%Atname(j)) == u_case(label(1)(npos+1:npos+6))) then
             na=j
             exit
          end if
       end do

       nb=0
       do j=1,ncode
          if (u_case(label(1)(1:npos))==u_case(trim(code_nam(j)))) then
             nb=j
             exit
          end if
       end do
       if (nb == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Code-name not found for parameter name: "//trim(label(1))
          return
       end if

       !---- Checking ----!
       if (nb < 6 .and. na == 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Incompatible relation: "//trim(label(1))
          return
       end if

       nl=0
       select case (nb)
          case ( 1:3)
             !---- X_, Y_, Z_ ----!
             fac_0=molec%mI_coor(nb,na)
                nl=molec%lI_coor(nb,na)
          case ( 4)
             !---- Biso_ ----!
             fac_0=molec%mbiso(na)
                nl=molec%lbiso(na)
          case ( 5)
             !---- Occ_ ----!
             fac_0=molec%mocc(na)
                nl=molec%locc(na)
          case ( 6:12)
             !---- Anisotropic Parameters ----!
             err_refcodes=.true.
             ERR_RefCodes_Mess="Incompatible Code-name for parameter name: "//trim(label(1))
             return
          case (13:15)
             !---- Xc_, Yc_, Zc_ ----!
             fac_0=molec%mxcentre(nb-12)
                nl=molec%lxcentre(nb-12)
          case (16:18)
             !---- Theta_, Phi_, Chi_ ----!
             fac_0=molec%mxcentre(nb-15)
                nl=molec%lxcentre(nb-15)
          case (19:21)
             !!! Not Yet Implemented !!!
       end select ! nb

       if (nb < ncode) then
          if (nl == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No refinable parameter was selected for "//trim(label(1))
             return
          end if
       end if

       !---- Set Others ----!
       n=1
       do
          n=n+1
          if (n > ic) exit

          npos=index(label(n),"_")
          if (npos ==0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No CrysFML code-name was found for "//trim(label(n))
             return
          end if

          nc=0
          do j=1,molec%natoms
             if (u_case(molec%Atname(j)) == u_case(label(n)(npos+1:npos+6))) then
                nc=j
                exit
             end if
          end do

          nd=0
          do j=1,ncode
             if (u_case(label(n)(1:npos))==u_case(trim(code_nam(j)))) then
                nd=j
                exit
             end if
          end do
          if (nd == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Code-name not found for "//trim(label(n))
             return
          end if

          n=n+1
          call getnum(label(n),vet,ivet,iv)
          if (iv == 1) then
             fac_1=vet(1)
          else
             fac_1=fac_0
             n=n-1
          end if

          !---- Checking ----!
          if (nd < 6 .and. nc == 0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Incompatible relation: "//trim(label(n))
             return
          end if

          select case (nd)
             case ( 1:3)
                !---- X_, Y_, Z_ ----!
                nl2=molec%lI_coor(nd,nc)
                call Delete_RefCodes(nl2,molec)
                molec%mI_coor(nd,nc)=fac_1
                molec%lI_coor(nd,nc)=nl
             case ( 4)
                !---- Biso_ ----!
                nl2=molec%lbiso(nc)
                call Delete_RefCodes(nl2,molec)
                molec%mbiso(nc)=fac_1
                molec%lbiso(nc)=nl
             case ( 5)
                !---- Occ_ ----!
                nl2=molec%locc(nc)
                call Delete_RefCodes(nl2,molec)
                molec%mocc(nc)=fac_1
                molec%locc(nc)=nl
             case ( 6:12)
                err_refcodes=.true.
                ERR_RefCodes_Mess="Incompatible Code-name for "//trim(label(n))
                return
             case (13:15)
                !---- Xc_, Yc_, Zc_ ----!
                nl2=molec%lxcentre(nc-12)
                call Delete_RefCodes(nl2,molec)
                molec%mxcentre(nc-12)=fac_1
                molec%lxcentre(nc-12)=nl
             case (16:18)
                !---- Theta_, Phi_, Chi_ ----!
                nl2=molec%lorient(nc-15)
                call Delete_RefCodes(nl2,molec)
                molec%morient(nc-15)=fac_1
                molec%lorient(nc-15)=nl
             case (19:21)
                !!! not yet implemented !!!
          end select ! nb
          np_cons=np_cons+1

       end do

       return
    End Subroutine Get_ConCodes_Line_Molec

    !!--++
    !!--++ Subroutine Get_ConCodes_Line_Magdom(Line,Mag_dom)
    !!--++    character(len=*),         intent(in)          :: Line
    !!--++    type(Magnetic_Domain_type),    intent(in out) :: Mag_dom
    !!--++
    !!--++ Get the Magnetic Constraints Relations: Presently only 'equal'
    !!--++ related to magnetic domains
    !!--++ Created: February - 2012
    !!
    Subroutine Get_ConCodes_Line_Magdom(Line,Mag_dom)
       !---- Arguments ----!
       character(len=*),     intent(in)          :: Line
       type(Magnetic_Domain_type),intent(in out) :: Mag_dom

       !---- Local variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: ic,n,na,nb,ina,inb,npos,ich
       integer                          :: nl,nl2,iv
       integer, dimension(1)            :: ivet
       real(kind=cp)                    :: fac_0,fac_1
       real(kind=cp),dimension(1)       :: vet

       call init_err_refcodes()

       !---- Check is chirality is present ----!
       if (Mag_Dom%chir) then
        ich=2
       else
        ich=1
       end if
       nl=0
       call getword(line,label,ic)
       if (ic < 2) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="EQUAL keyword needs two labels: "//trim(line)
          return
       end if

       !---- Set Futher Information----!
       !---- Na is the number of S-domains
       !---- Nb is the number of chiral domains
       !---- Fac0 is the multiplier
       !---- Nl is the number of refinement parameter

       npos=index(label(1),"d") ! identifying magDom word
       if (npos ==0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="The name "//trim(label(1))//" does not fit any known code-name of CrysFML "
          return
       end if
          do na=1,Mag_Dom%nd
           do nb=1,ich ! npos+8 as magdom0N should have length magD + om0N
             if (u_case(Mag_Dom%Lab(nb,na)) == u_case(label(1)(1:npos+8))) then
             fac_0=Mag_Dom%Mpop(nb,na)
                nl=Mag_Dom%Lpop(nb,na)
                exit
             end if

           end do
          end do

       !---- Set the rest elements in Contsraints ----!
       n=1
       do
          n=n+1
          if (n > ic) exit

          npos=index(label(n),"d")! identifying magDom word
          if (npos ==0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="No CrysFML code-name was found for "//trim(label(n))
             return
          end if

          do na=1,Mag_Dom%nd
           do nb=1,ich ! npos+8 as magdom0N should have length magD + om0N
             if (u_case(trim(Mag_Dom%Lab(nb,na))) == u_case(label(n)(1:npos+8))) then
             ina=na
             inb=nb
                exit
             end if
           end do
          end do

          !---- Is there a new multiplier?
          n=n+1
          call getnum(label(n),vet,ivet,iv)
          if (iv == 1) then
             fac_1=vet(1)
          else
             fac_1=fac_0
             n=n-1
          end if

            nl2=Mag_Dom%Lpop(inb,ina)
            call Delete_refCodes(nl2,Mag_dom)
            Mag_Dom%Mpop(inb,ina)=fac_1
            Mag_Dom%Lpop(inb,ina)=nl

          np_cons=np_cons+1

       end do

       return
    End Subroutine Get_ConCodes_Line_Magdom

    !!--++
    !!--++ Subroutine Get_RefCodes_Line(Key,Dire,Line,FAtom/FmAtom/MolCrys/Molec/MagDom,Spg)
    !!--++    integer,                      intent(in)     :: Key
    !!--++    character(len=*),             intent(in)     :: Dire
    !!--++    character(len=*),             intent(in)     :: Line
    !!--++    type(Atom_List_Type),         intent(in out) :: FAtom
    !!--++    or
    !!--++    type(mAtom_List_Type),        intent(in out) :: FmAtom
    !!--++    or
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++    or
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    or
    !!--++    type(Magnetic_Domain_type),   intent(in out) :: Magdom
    !!--++
    !!--++    type(space_group_type),       intent(in)     :: Spg
    !!--++
    !!--++    (Private)
    !!--++    Get Refinement Codes for Free atoms type
    !!--++
    !!--++ Update: March - 2005
    !!

    !!--++
    !!--++ Subroutine Get_RefCodes_Line_FAtom(Key,Dire,Line,FAtom,Spg)
    !!--++    integer,                 intent(in)     :: Key
    !!--++    character(len=*),        intent(in)     :: Dire
    !!--++    character(len=*),        intent(in)     :: Line
    !!--++    type(Atom_List_Type),    intent(in out) :: FAtom
    !!--++    type(space_group_type),  intent(in)     :: Spg
    !!--++
    !!--++ Overloaded
    !!--++ Get Refinement Codes for Free atoms type
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RefCodes_Line_FAtom(Key,Dire,Line,FAtom,Spg)
       !---- Arguments ----!
       integer,                 intent(in)     :: Key
       character(len=*),        intent(in)     :: Dire
       character(len=*),        intent(in)     :: Line
       type(Atom_List_Type),    intent(in out) :: FAtom
       type(space_group_type),  intent(in)     :: Spg

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: i,j,n,na,nb,ndir,npos,nlong,ic !,k,nc
       integer                          :: icond,iv,n_ini,n_end
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real(kind=cp)                    :: x_low,x_up,x_step
       real(kind=cp),dimension(5)       :: vet

       call init_err_refcodes()
       nlong=len_trim(line)

       if (nlong ==0) then
          !---- Default Values ----!
          do i=1,FAtom%natoms
             call Fill_RefCodes(Key,Dire,i,0,0.0_cp,0.0_cp,0.0_cp,0,Fatom,Spg)
          end do

       else
          !---- VARY/FIX Line: [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1)

          if (ndir <=0) then
             !--- [INF,[SUP,[STEP,[COND]]]] ----!
             call getnum(line,vet,ivet,iv)
             select case (iv)
                case (1)
                   x_low=vet(1)
                    x_up=vet(1)
                  x_step=0.0
                   icond=0

                case (2)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=0.0
                   icond=0

                case (3)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=0

                case (4)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=ivet(4)

                case default
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Only numbers in "//trim(line)
                   return
             end select

             nb=0
             do na=1,FAtom%natoms
                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,fatom,spg)
             end do
             if (err_refcodes) return

          else
             !---- [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
             ! If Ilabel(i) == 0 then is a label otherwise is a number
             n_ini=minloc(ilabel,dim=1)
             ilabel(n_ini)=2
             n_end=minloc(ilabel,dim=1)-1

             do n=1,ndir
                na=0
                nb=0

                !---- Default values ----!
                x_low =0.0
                x_up  =0.0
                x_step=0.0
                icond =0

                !---- Label ----!
                npos=index(label(n_ini),"_")
                if (npos >0) then
                   do j=1,ncode
                      if (u_case(label(n_ini)(1:npos))==u_case(trim(code_nam(j)))) then
                         nb=j
                         exit
                      end if
                   end do
                   if (nb == 0) then
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Code-name not found for "//trim(label(n_ini))
                      return
                   end if
                end if

                do j=1,FAtom%natoms
                   if (u_case(fatom%atom(j)%lab) == u_case(label(n_ini)(npos+1:npos+6))) then
                      na=j
                      exit
                   end if
                end do
                if (na == 0) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess=" Atom label not found for "//trim(line)
                   return
                end if

                !---- Valu List: Inf,Sup,Step,Cond ----!
                i=0
                do j=n_ini+1,n_end
                   i=i+1
                   call getnum(label(j),vet,ivet,iv)
                   select case (i)
                      case (1)
                         x_low=vet(1)
                         x_up =vet(1)

                      case (2)
                         x_up =vet(1)

                      case (3)
                         x_step =vet(1)

                      case (4)
                         icond = ivet(1)
                   end select
                end do

                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,fatom,spg)
                if (err_refcodes) return

                n_ini=minloc(ilabel,dim=1)
                ilabel(n_ini)=2
                n_end=minloc(ilabel,dim=1)-1

             end do
          end if

       end if

       return
    End Subroutine Get_RefCodes_Line_FAtom

    !!--++
    !!--++ Subroutine Get_RefCodes_Line_FmAtom(Key,Dire,Line,FmAtom)
    !!--++    integer,                 intent(in)     :: Key
    !!--++    character(len=*),        intent(in)     :: Dire
    !!--++    character(len=*),        intent(in)     :: Line
    !!--++    type(mAtom_List_Type),   intent(in out) :: FmAtom
    !!--++
    !!--++ Get Refinement Codes for Free Magnetic Atom type
    !!--++ magnetic clone of Get_RefCodes_Line_FAtom
    !!--++ Created: December - 2011
    !!--++ Updated: February - 2012
    !!
    Subroutine Get_RefCodes_Line_FmAtom(Key,Dire,Line,FmAtom)
       !---- Arguments ----!
       integer,                 intent(in)     :: Key
       character(len=*),        intent(in)     :: Dire
       character(len=*),        intent(in)     :: Line
       type(mAtom_List_Type),   intent(in out) :: FmAtom

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: i,j,n,na,nb,ndir,npos,nlong,ic
       integer                          :: icond,iv,n_ini,n_end
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real(kind=cp)                    :: x_low,x_up,x_step
       real(kind=cp),dimension(5)       :: vet


       call init_err_refcodes()

       nlong=len_trim(line)
       if (nlong ==0) then
          !---- Default Values ----!
          do i=1,FmAtom%natoms
             call Fill_Refcodes(Key,Dire,i,0,0.0_cp,0.0_cp,0.0_cp,0,FmAtom)
          end do

       else
          !---- VARY/FIX Line: [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1)

          if (ndir <=0) then
             !--- [INF,[SUP,[STEP,[COND]]]] ----!
             call getnum(line,vet,ivet,iv)
             select case (iv)
                case (1)
                   x_low=vet(1)
                    x_up=vet(1)
                  x_step=0.0
                   icond=0

                case (2)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=0.0
                   icond=0

                case (3)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=0

                case (4)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=ivet(4)

                case default
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Only numbers in "//trim(line)
                   return
             end select

             nb=0
             do na=1,FmAtom%natoms
                call Fill_Refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,FmAtom)
             end do
             if (err_refcodes) return

          else

             !---- [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
             ! If Ilabel(i) == 0 then is a label otherwise is a number
             n_ini=minloc(ilabel,dim=1)
             ilabel(n_ini)=2
             n_end=minloc(ilabel,dim=1)-1

             do n=1,ndir
                na=0
                nb=0

                !---- Default values ----!
                x_low =0.0
                x_up  =0.0
                x_step=0.0
                icond =0

                !---- Label ----!
                npos=index(label(n_ini),"_")
                if (npos >0) then
                   do j=1,mncode
                      if (u_case(label(n_ini)(1:npos))==u_case(trim(mcode_nam(j)))) then
                         nb=j
                         exit
                      end if
                   end do
                   if (nb == 0) then
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Code-name not found for "//trim(label(n_ini))
                      return
                   end if
                end if

                do j=1,FmAtom%natoms
                   if (u_case(FmAtom%atom(j)%lab) == u_case(label(n_ini)(npos+1:npos+6))) then
                      na=j
                      exit
                   end if
                end do
                if (na == 0) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess=" Atom label not found for "//trim(line)
                   return
                end if

                !---- Valu List: Inf,Sup,Step,Cond ----!
                i=0
                do j=n_ini+1,n_end
                   i=i+1
                   call getnum(label(j),vet,ivet,iv)
                   select case (i)
                      case (1)
                         x_low=vet(1)
                         x_up =vet(1)

                      case (2)
                         x_up =vet(1)

                      case (3)
                         x_step =vet(1)

                      case (4)
                         icond = ivet(1)
                   end select
                end do

                call Fill_Refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,FmAtom)
                if (err_refcodes) return

                n_ini=minloc(ilabel,dim=1)
                ilabel(n_ini)=2
                n_end=minloc(ilabel,dim=1)-1

             end do
          end if

       end if

       return
    End Subroutine Get_RefCodes_Line_FmAtom

    !!--++
    !!--++ Subroutine Get_RefCodes_Line_Molcrys(Key,Dire,Line,Molcrys,NMol)
    !!--++   integer,                      intent(in)     :: Key
    !!--++   character(len=*),             intent(in)     :: Dire
    !!--++   character(len=*),             intent(in)     :: Line
    !!--++   type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++   integer,                      intent(in)     :: NMol
    !!--++
    !!--++ Overloaded
    !!--++ Get Refinement Codes for Free atoms type
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RefCodes_Line_Molcrys(Key,Dire,Line,Molcrys,NMol)
       !---- Arguments ----!
       integer,                      intent(in)     :: Key
       character(len=*),             intent(in)     :: Dire
       character(len=*),             intent(in)     :: Line
       type(molecular_Crystal_type), intent(in out) :: MolCrys
       integer,                      intent(in)     :: NMol

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: i,j,k,n,na,nb,ndir,npos,nlong,ic !,nc
       integer                          :: icond,iv,n_ini,n_end
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real(kind=cp)                    :: x_low,x_up,x_step
       real(kind=cp),dimension(5)       :: vet


       call init_err_refcodes()

       nlong=len_trim(line)
       if (nlong ==0) then
          !---- Default values ----!
          select case (NMol)
             case (-1)
                !---- No Molecule Information ----!
                do i=1,molcrys%n_free
                   call Fill_RefCodes(Key,Dire,i,0,0.0_cp,0.0_cp,0.0_cp,0,Molcrys,NMol)
                end do

             case (0)
                do k=1,molcrys%n_mol
                   do i=1,molcrys%mol(k)%natoms
                      if (k > 1) then
                         na=molcrys%n_free+sum(molcrys%mol(1:k-1)%natoms)+i
                      else
                         na=molcrys%n_free+i
                      end if
                      call Fill_RefCodes(Key,Dire,na,0,0.0_cp,0.0_cp,0.0_cp,0,Molcrys,NMol)
                   end do
                end do

             case (1:)
                do i=1,molcrys%mol(nmol)%natoms
                   if (nmol > 1) then
                      na=molcrys%n_free+sum(molcrys%mol(1:nmol-1)%natoms)+i
                   else
                      na=molcrys%n_free+i
                   end if
                   call Fill_RefCodes(Key,Dire,na,0,0.0_cp,0.0_cp,0.0_cp,0,Molcrys,NMol)
                end do
          end select

       else
          !---- VARY/FIX Line: [LABEL,[INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1)

          if (ndir <=0) then
             !---- [INF,[SUP,[STEP,[COND]]]] ----!
             call getnum(line,vet,ivet,iv)
             select case (iv)
                case (1)
                   x_low=vet(1)
                    x_up=vet(1)
                  x_step=0.0
                   icond=0

                case (2)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=0.0
                   icond=0

                case (3)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=0

                case (4)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=ivet(4)

                case default
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Only numbers in "//trim(line)
                   return
             end select

             nb=0
             select case (nmol)
                case (-1)
                   !---- No Molecule Information ----!
                   do na=1,molcrys%n_free
                      call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,Molcrys,NMol)
                   end do

                case ( 0)
                   !---- For all Molecules Defined ----!
                   do n=1,molcrys%n_mol
                      do i=1,molcrys%mol(n)%natoms
                         if (n > 1) then
                            na=molcrys%n_free+sum(molcrys%mol(1:n-1)%natoms)+i
                         else
                            na=molcrys%n_free+i
                         end if
                      end do
                      call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,Molcrys,NMol)
                   end do

                case (1:)
                   !---- Particular molecule ----!
                   do i=1,molcrys%mol(nmol)%natoms
                      if (nmol > 1) then
                         na=molcrys%n_free+sum(molcrys%mol(1:nmol-1)%natoms)+i
                      else
                         na=molcrys%n_free+i
                      end if
                      call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,Molcrys,NMol)
                   end do
             end select
             if (err_refcodes) return

          else
             !---- [LABEL,[INF,[SUP,[STEP,[COND]]]]] ----!
             ! If Ilabel(i) == 0 then is a label otherwise is a number
             n_ini=minloc(ilabel,dim=1)
             ilabel(n_ini)=2
             n_end=minloc(ilabel,dim=1)-1

             do n=1,ndir
                na=0
                nb=0

                !---- Deafult values ----!
                x_low =0.0
                x_up  =0.0
                x_step=0.0
                icond =0

                !---- Label ----!
                npos=index(label(n_ini),"_")
                if (npos >0) then
                   do j=1,ncode
                      if (u_case(label(n_ini)(1:npos))==u_case(trim(code_nam(j)))) then
                         nb=j
                         exit
                      end if
                   end do
                   if (nb == 0) then
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Code-name not found for "//trim(label(n_ini))
                      return
                   end if
                end if

                select case (nmol)
                   case (-1)
                      !---- No molecule Information ----!
                      do j=1,molcrys%n_free
                         if (u_case(molcrys%atm(j)%lab) == u_case(label(n_ini)(npos+1:npos+6))) then
                            na=j
                            exit
                         end if
                      end do

                   case ( 0)
                      !---- For all molecules defined ----!
                      do i=1,molcrys%n_mol
                         do j=1,molcrys%mol(i)%natoms
                            if (u_case(molcrys%mol(i)%Atname(j)) == u_case(label(n_ini)(npos+1:npos+6))) then
                               if (j > 1) then
                                  na=molcrys%n_free+sum(molcrys%mol(1:j-1)%natoms)+j
                               else
                                  na=molcrys%n_free+j
                               end if
                               exit
                            end if
                         end do
                      end do

                   case (1:)
                      !---- Particular Molecule ----!
                      do j=1,molcrys%mol(nmol)%natoms
                         if (u_case(molcrys%mol(nmol)%Atname(j)) == u_case(label(n_ini)(npos+1:npos+6))) then
                            if (j > 1) then
                               na=molcrys%n_free+sum(molcrys%mol(1:j-1)%natoms)+j
                            else
                               na=molcrys%n_free+j
                            end if
                            exit
                         end if
                      end do
                end select
                if (na == 0) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess=" Atom label not found for "//trim(line)
                   return
                end if

                !---- Valu List: Inf,Sup,Step,Cond ----!
                i=0
                do j=n_ini+1,n_end
                   i=i+1
                   call getnum(label(j),vet,ivet,iv)
                   select case (i)
                      case (1)
                         x_low=vet(1)
                         x_up =vet(1)

                      case (2)
                         x_up =vet(1)

                      case (3)
                         x_step =vet(1)

                      case (4)
                         icond = ivet(1)
                   end select
                end do

                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,molcrys,NMol)
                if (err_refcodes) return

                n_ini=minloc(ilabel,dim=1)
                ilabel(n_ini)=2
                n_end=minloc(ilabel,dim=1)-1

             end do
          end if
       end if

       return
    End Subroutine Get_RefCodes_Line_Molcrys

    !!--++
    !!--++ Subroutine Get_RefCodes_Line_Molec(Key,Dire,Line,Molec,Spg)
    !!--++    integer,                      intent(in)     :: Key
    !!--++    character(len=*),             intent(in)     :: Dire
    !!--++    character(len=*),             intent(in)     :: Line
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    type(space_group_type),       intent(in)     :: Spg
    !!--++
    !!--++ Overloaded
    !!--++ Get Refinement Codes for Free atoms type
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RefCodes_Line_Molec(Key,Dire,Line,Molec,Spg)
       !---- Arguments ----!
       integer,                      intent(in)     :: Key
       character(len=*),             intent(in)     :: Dire
       character(len=*),             intent(in)     :: Line
       type(molecule_type),          intent(in out) :: Molec
       type(space_group_type),       intent(in)     :: Spg

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: i,j,n,na,nb,ndir,npos,nlong,ic !,k,nc
       integer                          :: icond,iv,n_ini,n_end
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real(kind=cp)                    :: x_low,x_up,x_step
       real(kind=cp),dimension(5)       :: vet

       call init_err_refcodes()

       nlong=len_trim(line)

       if (nlong ==0) then
          !---- Default values ----!
          do i=1,molec%natoms
             call Fill_RefCodes(key,dire,i,0,0.0_cp,0.0_cp,0.0_cp,0,molec,spg)
          end do

       else
          !---- VARY/FIX Line: [LABEL,[INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1)

          if (ndir <=0) then
             !---- [INF,[SUP,[STEP,[COND]]]] ----!
             call getnum(line,vet,ivet,iv)
             select case (iv)
                case (1)
                   x_low=vet(1)
                    x_up=vet(1)
                  x_step=0.0
                   icond=0

                case (2)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=0.0
                   icond=0

                case (3)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=0

                case (4)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=ivet(4)

                case default
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Only numbers in "//trim(line)
                   return
             end select

             nb=0
             do na=1,molec%natoms
                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,molec,spg)
             end do
             if (err_refcodes) return

          else
             !---- [LABEL,[INF,[SUP,[STEP,[COND]]]]] ----!
             ! If Ilabel(i) == 0 then is a label otherwise is a number
             n_ini=minloc(ilabel,dim=1)
             ilabel(n_ini)=2
             n_end=minloc(ilabel,dim=1)-1

             do n=1,ndir
                na=0
                nb=0

                !---- Default values ----!
                x_low =0.0
                x_up  =0.0
                x_step=0.0
                icond =0

                !---- Label ----!
                npos=index(label(n_ini),"_")
                if (npos >0) then
                   do j=1,ncode
                      if (u_case(label(n_ini)(1:npos))==u_case(trim(code_nam(j)))) then
                         nb=j
                         exit
                      end if
                   end do
                   if (nb == 0) then
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Code-name not found for "//trim(label(n_ini))
                      return
                   end if
                end if

                do j=1,molec%natoms
                   if (u_case(molec%AtName(j)) == u_case(label(n_ini)(npos+1:npos+6))) then
                      na=j
                      exit
                   end if
                end do
                if (na == 0) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess=" Atom label not found for "//trim(line)
                   return
                end if

                !---- Value List: Inf,Sup,Step,Cond ----!
                i=0
                do j=n_ini+1,n_end
                   i=i+1
                   call getnum(label(j),vet,ivet,iv)
                   select case (i)
                      case (1)
                         x_low=vet(1)
                         x_up =vet(1)

                      case (2)
                         x_up =vet(1)

                      case (3)
                         x_step =vet(1)

                      case (4)
                         icond = ivet(1)
                   end select
                end do

                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,molec,spg)
                if (err_refcodes) return

                n_ini=minloc(ilabel,dim=1)
                ilabel(n_ini)=2
                n_end=minloc(ilabel,dim=1)-1

             end do
          end if

       end if

       return
    End Subroutine Get_RefCodes_Line_Molec

    !!--++
    !!--++ Subroutine Get_RefCodes_Line_Magdom(Key,Dire,Line,Mag_dom)
    !!--++    integer,                   intent(in)     :: Key
    !!--++    character(len=*),          intent(in)     :: Dire
    !!--++    character(len=*),          intent(in)     :: Line
    !!--++    type(Magnetic_Domain_type),intent(in out) :: Mag_dom
    !!--++
    !!--++ Get Refinement Codes for Magnetic domain
    !!--++ Created: February - 2012
    !!
    !!
    Subroutine Get_RefCodes_Line_Magdom(Key,Dire,Line,Mag_dom)
       !---- Arguments ----!
       integer,                   intent(in)     :: Key
       character(len=*),          intent(in)     :: Dire
       character(len=*),          intent(in)     :: Line
       type(Magnetic_Domain_type),intent(in out) :: Mag_dom

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: i,na,nb,ndir,nlong,ic,ich
       integer                          :: icond,iv
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real(kind=cp)                    :: x_low,x_up,x_step
       real(kind=cp),dimension(5)       :: vet

       !---- Check is chirality is present ----!
        if (Mag_Dom%chir) then
         ich=2
        else
         ich=1
        end if

       call init_err_refcodes()

       nlong=len_trim(line)

       if (nlong ==0) then

          !---- Default Values ----!
        do i=1,Mag_Dom%nd*ich
             call Fill_Refcodes(Key,Dire,i,0,0.0_cp,0.0_cp,0.0_cp,0,Mag_dom)
        end do

       else
          !---- VARY/FIX Line: [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1) !number of refined domains

                !---- Default values ----!
                x_low =0.0
                x_up  =1.0
                x_step=0.01
                icond =1

          do na=1,Mag_Dom%nd
           do nb=1,ich
            do i=1,ic
             if ( u_case(Mag_dom%lab(nb,na)) == u_case(label(i)) ) &
              call Fill_Refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,Mag_dom)
            end do
           end do
          end do

             return
       end if

       return
    End Subroutine Get_RefCodes_Line_Magdom

    !!----
    !!---- Subroutine Get_RefGCodes_Line(Key,Dire,Line,namp,model,Sys,Iphas)
    !!----    integer,                             intent(in)     :: Key
    !!----    character(len=*),                    intent(in)     :: Dire ! "gvary" or "gfix"
    !!----    character(len=*),                    intent(in)     :: Line ! directive after gvary or gfix
    !!----    character(len=*),                    intent(in)     :: namp ! name of the parameter to be fixed or refined
    !!----    type(Nonatomic_Parameter_List_Type), intent(in out) :: model
    !!----    character(len=*), optional,          intent( in)    :: sys
    !!----    integer,          optional,          intent( in)    :: Iphas
    !!----
    !!---- Get Refinement Codes for non-atomic parameters in the current line
    !!----
    !!---- Update: November 1 - 2013
    !!
    Subroutine Get_RefGCodes_Line(Key,Dire,Line,namp,model,sys,Iphas)
       integer,                             intent(in)     :: Key
       character(len=*),                    intent(in)     :: Dire
       character(len=*),                    intent(in)     :: Line
       character(len=*),                    intent(in)     :: namp
       type(Nonatomic_Parameter_List_Type), intent(in out) :: model
       character(len=*), optional,          intent(in)     :: sys
       integer,          optional,          intent( in)    :: Iphas

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       character(len=20)                :: new_label,aux_string
       character(len=2)                 :: phase
       integer                          :: i,j,n,na,nb,ndir,nlong,ic !,k,nc
       integer                          :: icond,iv,n_ini,n_end
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real(kind=cp)                    :: x_low,x_up,x_step
       real(kind=cp),dimension(5)       :: vet

       call init_err_refcodes()
       nlong=len_trim(line)  !length of the directive, eg. for "cell_a  5 6 0.1 0", nlong=17

       if (nlong ==0) then  !In this case "key" cannot be zero and "namp" does not contain numbers
          !---- Default Values ----!
          if(present(sys)) then
            if(present(Iphas)) then
               call Fill_RefGCodes(Key,Dire,namp,0.0_cp,0.0_cp,0.0_cp,0,model,sys,Iphas)
            else
               call Fill_RefGCodes(Key,Dire,namp,0.0_cp,0.0_cp,0.0_cp,0,model,sys)
            end if
          else
            if(present(Iphas)) then
               call Fill_RefGCodes(Key,Dire,namp,0.0_cp,0.0_cp,0.0_cp,0,model,Iphas=Iphas)
            else
               call Fill_RefGCodes(Key,Dire,namp,0.0_cp,0.0_cp,0.0_cp,0,model)
            end if
          end if

       else
          !---- GVARY/GFIX Line: [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1)  !Counts the number of zeroes  (label(1) contains always the name of the parameter)
                                        !This means the number of keywords

          if (ndir <= 0) then  !label(i) contain only numbers (no more keywords!)
             !--- [INF,[SUP,[STEP,[COND]]]] ----! It corresponds to a generic key /= 0 keyword plus numbers
             call getnum(line,vet,ivet,iv)

             select case (iv)
                case (1)
                   x_low=vet(1)
                    x_up=vet(1)
                  x_step=0.0
                   icond=0

                case (2)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=0.0
                   icond=0

                case (3)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=0

                case (4)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=ivet(4)

                case default
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Maximum of 4 numbers in "//trim(namp)
                   return
             end select

             if(present(sys)) then
               if(present(Iphas)) then
                  call Fill_refGcodes(key,dire,namp,x_low,x_up,x_step,icond,model,sys,Iphas)
               else
                   call Fill_refGcodes(key,dire,namp,x_low,x_up,x_step,icond,model,sys) !dire="gvary" or "gfix", namp may contain numbers
               end if
             else
               if(present(Iphas)) then
                  call Fill_refGcodes(key,dire,namp,x_low,x_up,x_step,icond,model,Iphas=Iphas)
               else
                  call Fill_refGcodes(key,dire,namp,x_low,x_up,x_step,icond,model)
               end if
             end if
             if (err_refcodes) return

          else  !Now there are more keywords and eventually numbers
             !---- [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
             ! If Ilabel(i) == 0 then is a label otherwise is a number
             n_ini=minloc(ilabel,dim=1)
             ilabel(n_ini)=2
             n_end=minloc(ilabel,dim=1)-1

             do n=1,ndir  !This runs over the number of directives (keywords, names of parameters)
                nb=0
                !---- Default values ----!
                x_low =0.0
                x_up  =0.0
                x_step=0.0
                icond =0
                aux_string=u_case(label(n_ini))

                if(present(Iphas)) then
                  write(unit=phase,fmt="(i2.2)") iphas
                  if(index(aux_string,phase) == 0) then
                    aux_string=trim(aux_string)//"_"//phase
                  end if
                end if

                do j=1,model%npar
                   if (index(u_case(model%Par(j)%nam),trim(aux_string)) /= 0 ) then
                      na=j
                      exit
                   end if
                end do

                if (na == 0) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess=" NonAtomic label not found for "//trim(line)
                   return
                else
                   new_label=model%Par(na)%nam  !this is the name of the parameter to be fixed or refined
                end if

                !---- Valu List: Inf,Sup,Step,Cond ----!
                i=0
                do j=n_ini+1,n_end
                   i=i+1
                   call getnum(label(j),vet,ivet,iv)
                   select case (i)
                      case (1)
                         x_low=vet(1)
                         x_up =vet(1)

                      case (2)
                         x_up =vet(1)

                      case (3)
                         x_step =vet(1)

                      case (4)
                         icond = ivet(1)
                   end select
                end do

                if(present(sys)) then
                  if(present(Iphas)) then
                     call Fill_refGcodes(key,dire,new_label,x_low,x_up,x_step,icond,model,sys,Iphas=Iphas)
                  else
                     call Fill_refGcodes(key,dire,new_label,x_low,x_up,x_step,icond,model,sys)
                  end if
                else
                  if(present(Iphas)) then
                     call Fill_refGcodes(key,dire,new_label,x_low,x_up,x_step,icond,model,Iphas=Iphas)
                  else
                     call Fill_refGcodes(key,dire,new_label,x_low,x_up,x_step,icond,model)
                  end if
                end if
                if (err_refcodes) return

                n_ini=minloc(ilabel,dim=1)
                ilabel(n_ini)=2
                n_end=minloc(ilabel,dim=1)-1

             end do
          end if

       end if

       return
    End Subroutine Get_RefGCodes_Line


    !!----
    !!---- Subroutine Get_RestAng_Line(Line, FAtom)
    !!----    character(len=*),        intent(in)     :: Line
    !!----    type(Atom_List_Type),    intent(in out) :: FAtom
    !!----
    !!----     Get Distance Restraints relations for Free atoms type
    !!----     Line: Angle [sig] At1a At1b At1c At2a At2b At2c....
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Get_RestAng_Line(Line, FAtom)
       !---- Arguments ----!
       character(len=*),        intent(in) :: Line
       type(Atom_List_Type),    intent(in) :: FAtom

       !---- Local variables ----!
       integer, parameter                :: np=30
       character(len=30),dimension(np)   :: dire
       character(len=8), dimension(2,np) :: symtrans
       integer,dimension(3,np)           :: p
       real                              :: ang,sig

       character(len=8), dimension(2)  :: car
       integer                         :: i,j,iv,nc,nr,n_ini,n_end,npos
       integer, dimension(3)           :: ivet
       real(kind=cp), dimension(3)     :: vet


       if (len_trim(line) == 0 .or. .not. allocated(ang_rest)) return

       !---- Description for each word ----!
       call getword(line,dire,nc)

       !---- Get Angle ----!
       call getnum(dire(1),vet,ivet,iv)
       if (iv /= 1) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Error in AFIX line: "//trim(line)
          return
       end if
       ang=vet(1)

       !---- Get Sigma ----!
       call getnum(dire(2),vet,ivet,iv)
       if (iv /= 1) then
          sig=0.2
          n_ini=2
       else
          sig=max(vet(1),0.001_cp)
          n_ini=3
       end if

       nr=0
       symtrans=" "
       do i=n_ini,nc,3
          ivet=0
          car=" "
          npos=index(dire(i),"_")
          if (npos /=0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess=" The first atom in AFIX command must belong to the asymmetric unit: "//trim(Line)
             return
          end if
          npos=index(dire(i+1),"_")
          if (npos /=0) then
             car(1)=dire(i+1)(npos:)
             dire(i+1)=dire(i+1)(1:npos-1)
          end if
          npos=index(dire(i+2),"_")
          if (npos /=0) then
             car(2)=dire(i+2)(npos:)
             dire(i+2)=dire(i+2)(1:npos-1)
          end if

          do j=1,FAtom%natoms
             if (trim(u_case(dire(i))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(1)=j
             end if
             if (trim(u_case(dire(i+1))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(2)=j
             end if
             if (trim(u_case(dire(i+2))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(3)=j
             end if
             if (all(ivet > 0) ) exit
          end do
          if (any(ivet == 0)) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="  Some atom names in "//trim(line)//" not found in the asymmetric unit"
             return
          end if

          !---- New Relation ----!
          nr=nr+1
          p(:,nr)=ivet
          symtrans(:,nr)=car
       end do
       if (nr <= 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Illegal AFIX command  "//trim(line)
          return
       end if

       !---- Adding relations ----!
       n_ini=np_rest_ang+1
       n_end=np_rest_ang+nr
       ang_rest(n_ini:n_end)%aobs=ang
       ang_rest(n_ini:n_end)%acalc=0.0
       ang_rest(n_ini:n_end)%sigma=sig
       ang_rest(n_ini:n_end)%p(1) = p(1,1:nr)
       ang_rest(n_ini:n_end)%p(2) = p(2,1:nr)
       ang_rest(n_ini:n_end)%p(3) = p(3,1:nr)
       ang_rest(n_ini:n_end)%STCode(1)=symtrans(1,1:nr)
       ang_rest(n_ini:n_end)%STCode(2)=symtrans(2,1:nr)
       np_rest_ang=n_end

       return
    End Subroutine Get_RestAng_Line


    !!----
    !!---- Subroutine Get_RestDis_Line(Line, FAtom)
    !!----    character(len=*),        intent(in)     :: Line
    !!----    type(Atom_List_Type),    intent(in out) :: FAtom
    !!----
    !!----    Get Distance Restraints relations for Free atoms type
    !!----    Line: Dist [sig] At1a At1b At2a At2b......
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Get_RestDis_Line(Line, FAtom)
       !---- Arguments ----!
       character(len=*),        intent(in) :: Line
       type(Atom_List_Type),    intent(in) :: FAtom

       !---- Local variables ----!
       integer, parameter              :: np=20
       character(len=30),dimension(np) :: dire
       character(len=8), dimension(np) :: symtrans
       integer,dimension(2,np)         :: p
       real                            :: dis,sig

       character(len=8)                :: car
       integer                         :: i,j,iv,nc,nr,n_ini,n_end,npos
       integer, dimension(2)           :: ivet
       real(kind=cp), dimension(2)     :: vet


       if (len_trim(line) == 0 .or. .not. allocated(dis_rest) ) return

       !---- Description for each word ----!
       call getword(line,dire,nc)

       !---- Get Dist ----!
       call getnum(dire(1),vet,ivet,iv)
       if (iv /= 1) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Error in DFIX line: "//trim(line)
          return
       end if
       dis=vet(1)

       !---- Get Sigma ----!
       call getnum(dire(2),vet,ivet,iv)
       if (iv /= 1) then
          sig=0.02
          n_ini=2
       else
          sig=max(vet(1),0.0001_cp)
          n_ini=3
       end if

       nr=0
       symtrans=" "
       do i=n_ini,nc,2
          ivet=0
          car=" "
          npos=index(dire(i),"_")
          if (npos /=0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess=" The first atom in DFIX command must belong to the asymmetric unit: "//trim(Line)
             return
          end if
          npos=index(dire(i+1),"_")
          if (npos /=0) then
             car=dire(i+1)(npos:)
             dire(i+1)=dire(i+1)(1:npos-1)
          end if

          do j=1,FAtom%natoms
             if (trim(u_case(dire(i))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(1)=j
             end if
             if (trim(u_case(dire(i+1))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(2)=j
             end if
             if (all(ivet > 0) ) exit
          end do
          if (any(ivet == 0)) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="  Some atom names in"//trim(line)//" not found in the asymmetric unit"
             return
          end if

          !---- New Relation ----!
          nr=nr+1
          p(:,nr)=ivet
          symtrans(nr)=car
       end do
       if (nr <= 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Illegal DFIX command  "//trim(line)
          return
       end if

       !---- Adding relations ----!
       n_ini=np_rest_dis+1
       n_end=np_rest_dis+nr
       dis_rest(n_ini:n_end)%dobs=dis
       dis_rest(n_ini:n_end)%dcalc=0.0
       dis_rest(n_ini:n_end)%sigma=sig
       dis_rest(n_ini:n_end)%p(1) = p(1,1:nr)
       dis_rest(n_ini:n_end)%p(2) = p(2,1:nr)
       dis_rest(n_ini:n_end)%STCode=symtrans(1:nr)
       np_rest_dis=n_end

       return
    End Subroutine Get_RestDis_Line

    !!----
    !!---- Subroutine Get_RestTor_Line(Line, FAtom)
    !!----    character(len=*),        intent(in)     :: Line
    !!----    type(Atom_List_Type),    intent(in out) :: FAtom
    !!----
    !!----    Get Torsion Restraints relations for Free atoms type
    !!----    Line: Torsion_Angle [sig] At1a At1b At1c At1d ...
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Get_RestTor_Line(Line, FAtom)
       !---- Arguments ----!
       character(len=*),        intent(in) :: Line
       type(Atom_List_Type),    intent(in) :: FAtom

       !---- Local variables ----!
       integer, parameter                :: np=30
       character(len=30),dimension(np)   :: dire
       character(len=8), dimension(3,np) :: symtrans
       integer,dimension(4,np)           :: p
       real                              :: tor,sig

       character(len=8), dimension(3)  :: car
       integer                         :: i,j,iv,nc,nr,n_ini,n_end,npos
       integer, dimension(4)           :: ivet
       real(kind=cp), dimension(4)     :: vet


       if (len_trim(line) == 0 .or. .not. allocated(tor_rest)) return

       !---- Description for each word ----!
       call getword(line,dire,nc)

       !---- Get Angle ----!
       call getnum(dire(1),vet,ivet,iv)
       if (iv /= 1) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Error in TFIX line: "//trim(line)
          return
       end if
       tor=vet(1)

       !---- Get Sigma ----!
       call getnum(dire(2),vet,ivet,iv)
       if (iv /= 1) then
          sig=0.5
          n_ini=2
       else
          sig=max(vet(1),0.02_cp)
          n_ini=3
       end if

       nr=0
       symtrans=" "
       do i=n_ini,nc,4
          ivet=0
          car=" "
          npos=index(dire(i),"_")
          if (npos /=0) then
             err_refcodes=.true.
             ERR_RefCodes_Mess=" The first atom in TFIX must belong to the asymmetric unit: "//trim(Line)
             return
          end if
          npos=index(dire(i+1),"_")
          if (npos /=0) then
             car(1)=dire(i+1)(npos:)
             dire(i+1)=dire(i+1)(1:npos-1)
          end if
          npos=index(dire(i+2),"_")
          if (npos /=0) then
             car(2)=dire(i+2)(npos:)
             dire(i+2)=dire(i+2)(1:npos-1)
          end if
          npos=index(dire(i+3),"_")
          if (npos /=0) then
             car(3)=dire(i+3)(npos:)
             dire(i+3)=dire(i+3)(1:npos-1)
          end if
          do j=1,FAtom%natoms
             if (trim(u_case(dire(i))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(1)=j
             end if
             if (trim(u_case(dire(i+1))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(2)=j
             end if
             if (trim(u_case(dire(i+2))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(3)=j
             end if
             if (trim(u_case(dire(i+3))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(4)=j
             end if
             if (all(ivet > 0) ) exit
          end do
          if (any(ivet == 0)) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="  Some atom names in"//trim(line)//" not found in the asymmetric unit"
             return
          end if

          !---- New Relation ----!
          nr=nr+1
          p(:,nr)=ivet
          symtrans(:,nr)=car
       end do
       if (nr <= 0) then
          err_refcodes=.true.
          ERR_RefCodes_Mess="Illegal TFIX command  "//trim(line)
          return
       end if

       !---- Adding relations ----!
       n_ini=np_rest_tor+1
       n_end=np_rest_tor+nr
       tor_rest(n_ini:n_end)%tobs=tor
       tor_rest(n_ini:n_end)%tcalc=0.0
       tor_rest(n_ini:n_end)%sigma=sig
       tor_rest(n_ini:n_end)%p(1) = p(1,1:nr)
       tor_rest(n_ini:n_end)%p(2) = p(2,1:nr)
       tor_rest(n_ini:n_end)%p(3) = p(3,1:nr)
       tor_rest(n_ini:n_end)%p(4) = p(4,1:nr)
       tor_rest(n_ini:n_end)%STCode(1)=symtrans(1,1:nr)
       tor_rest(n_ini:n_end)%STCode(2)=symtrans(2,1:nr)
       tor_rest(n_ini:n_end)%STCode(3)=symtrans(3,1:nr)
       np_rest_tor=n_end

       return
    End Subroutine Get_RestTor_Line

    !!----
    !!---- Subroutine Init_Err_RefCodes()
    !!----
    !!----    Initialize the errors flags in RefCodes
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Init_Err_RefCodes()

       err_refcodes=.false.
       ERR_RefCodes_Mess=" "

       return
    End Subroutine Init_Err_RefCodes

    !!----
    !!---- Subroutine Init_RefCodes(FAtom,FmAtom,Mag_dom,MolCrys,Molec,Model)
    !!----    type(Atom_List_Type),               optional,intent(in out) :: FAtom   ! Free Atom Object
    !!----    type(mAtom_List_Type),              optional,intent(in out) :: FmAtom  ! Magnetic Atom Object
    !!----    type(Magnetic_Domain_type),         optional,intent(in out) :: Mag_dom ! Magnetic domain object
    !!----    type(Molecular_Crystal_Type),       optional,intent(in out) :: MolCrys ! Molecular Crystal Object
    !!----    type(Molecule_Type),                optional,intent(in out) :: Molec   ! Molecule Object
    !!----    type(Nonatomic_Parameter_List_Type),optional,intent(in out) :: Model   !Non atomic parameter object
    !!----
    !!----
    !!---- Initialize all refinement codes. This subroutine has been merged with individual ones
    !!---- into a single subroutine with optional arguments. It seems more simple because the
    !!---- global set of parameters is nullified and initialized.
    !!----
    !!---- Update: November 2 - 2012
    !!
    Subroutine Init_RefCodes(FAtom,FmAtom,Mag_dom,MolCrys,Molec,Model)
       !---- Arguments ----!
       type(Atom_List_Type),               optional,intent(in out) :: FAtom   ! Free Atom Object
       type(mAtom_List_Type),              optional,intent(in out) :: FmAtom  ! Magnetic Atom Object
       type(Magnetic_Domain_type),         optional,intent(in out) :: Mag_dom ! Magnetic domain object
       type(Molecular_Crystal_Type),       optional,intent(in out) :: MolCrys ! Molecular Crystal Object
       type(Molecule_Type),                optional,intent(in out) :: Molec   ! Molecule Object
       type(Nonatomic_Parameter_List_Type),optional,intent(in out) :: Model   !Non atomic parameter object

       !---- Local variables ----!
       integer :: i

       !NP_Refi=0  !This has been removed because the initialization of V-arrays
       !NP_Cons=0  !is done in Allocate_VParam
       !
       !V_Vec   =0.0
       !V_Name=" "
       !V_Bounds=0.0
       !V_BCon  =0

       if(present(FAtom)) then
          do i=1,FAtom%Natoms
             FAtom%atom(i)%mx=0.0
             FAtom%atom(i)%lx=0

             FAtom%atom(i)%mbiso=0.0
             FAtom%atom(i)%lbiso=0

             FAtom%atom(i)%mocc=0.0
             FAtom%atom(i)%locc=0

             FAtom%atom(i)%mu=0.0
             FAtom%atom(i)%lu=0
          end do
       end if

       if(present(FmAtom)) then
          do i=1,FmAtom%Natoms

             FmAtom%atom(i)%mSkR=0.0
             FmAtom%atom(i)%lskr=0

             FmAtom%atom(i)%mSkI=0.0
             FmAtom%atom(i)%lski=0

             FmAtom%atom(i)%mmphas=0.0
             FmAtom%atom(i)%lmphas=0

             FmAtom%atom(i)%mbas=0.0
             FmAtom%atom(i)%lbas=0

          end do
       end if !present FmAtom

       if(present(Mag_dom)) then
          do i=1,Mag_Dom%nd

             Mag_Dom%Mpop(1:2,i)=0.0
             Mag_Dom%Lpop(1:2,i)=0

          end do
       end if  !present Mag_Dom

       if (present(MolCrys)) then

          if (MolCrys%N_Free > 0 .and. allocated(MolCrys%Atm)) then
             do i=1,MolCrys%N_Free
                MolCrys%Atm(i)%mx=0.0
                MolCrys%Atm(i)%lu=0

                MolCrys%Atm(i)%mbiso=0.0
                MolCrys%Atm(i)%lbiso=0

                MolCrys%Atm(i)%mocc=0.0
                MolCrys%Atm(i)%locc=0

                MolCrys%Atm(i)%mu=0.0
                MolCrys%Atm(i)%lu=0
             end do
          end if

          if (MolCrys%N_Mol > 0 .and. allocated(MolCrys%Mol)) then
             do i=1,MolCrys%N_Mol
                MolCrys%mol(i)%mxcentre=0.0
                MolCrys%mol(i)%lxcentre=0

                MolCrys%mol(i)%morient=0.0
                MolCrys%mol(i)%lorient=0

                MolCrys%mol(i)%mT_TLS=0.0
                MolCrys%mol(i)%lT_TLS=0

                MolCrys%mol(i)%mL_TLS=0.0
                MolCrys%mol(i)%lL_TLS=0

                MolCrys%mol(i)%mS_TLS=0.0
                MolCrys%mol(i)%lS_TLS=0

                if (MolCrys%Mol(i)%natoms > 0) then
                   MolCrys%mol(i)%mI_coor=0.0
                   MolCrys%mol(i)%lI_coor=0

                   MolCrys%mol(i)%mbiso=0.0
                   MolCrys%mol(i)%lbiso=0

                   MolCrys%mol(i)%mocc=0.0
                   MolCrys%mol(i)%locc=0
                end if
             end do
          end if
       end if !present MolCrys

       if(present(Molec)) then
          Molec%mxcentre=0.0
          Molec%lxcentre=0

          Molec%morient=0.0
          Molec%lorient=0

          Molec%mT_TLS=0.0
          Molec%lT_TLS=0

          Molec%mL_TLS=0.0
          Molec%lL_TLS=0

          Molec%mS_TLS=0.0
          Molec%lS_TLS=0

          do i=1,molec%natoms
             Molec%mI_coor(:,i)=0.0
             Molec%lI_coor(:,i)=0

             Molec%mbiso(i)=0.0
             Molec%lbiso(i)=0

             Molec%mocc(i)=0.0
             Molec%locc(i)=0
          end do
       end if !if present Molec

       if(present(Model)) then
         do i=1,Model%npar
            Model%par(i)%multip=0.0
            Model%par(i)%Lcode=0
         end do
       end if

       return
    End Subroutine Init_RefCodes

    !!----
    !!---- Subroutine Read_RefCodes_File(file_dat,n_ini,n_end,FAtom/MolCrys/Molec/MagStr,Spg)
    !!----    Type(file_list_type),         intent( in)    :: file_dat
    !!----    integer,                      intent( in)    :: n_ini
    !!----    integer,                      intent( in)    :: n_end
    !!----    type(Atom_List_Type),         intent(in out) :: FAtom
    !!----    or
    !!----    type(mAtom_List_Type),        intent(in out) :: FmAtom
    !!----    or
    !!----    type(molecular_crystal_type), intent(in out) :: molcrys
    !!----    or
    !!----    type(molecule_type),          intent(in out) :: molec
    !!----    type(Atom_List_Type),         intent(in out) :: FAtom
    !!----    or
    !!----    type(mAtom_List_Type),        intent(in out) :: FmAtom
    !!----    and type(Magnetic_Domain_type),intent(in out) :: Mag_dom
    !!----
    !!----    type(space_group_type),       intent(in)     :: Spg
    !!----
    !!----    Subroutine for treatment of Codes controls taken from FAtom/Molcrys/Molec
    !!----
    !!---- Update: March - 2005
    !!---- Update: February - 2012
    !!

    !!--++
    !!--++ Subroutine Read_RefCodes_File_FAtom(file_dat,n_ini,n_end,FAtom,Spg)
    !!--++    Type(file_list_type),    intent( in)    :: file_dat
    !!--++    integer,                 intent( in)    :: n_ini
    !!--++    integer,                 intent( in)    :: n_end
    !!--++    type(Atom_List_Type),    intent(in out) :: FAtom
    !!--++    type(space_group_type),  intent(in)     :: Spg
    !!--++
    !!--++    (Overloaded)
    !!--++    Subroutine for treatment of Codes controls
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Read_RefCodes_File_FAtom(file_dat,n_ini,n_end,FAtom,Spg)
       !---- Arguments ----!
       Type(file_list_type),     intent( in)    :: file_dat
       integer,                  intent( in)    :: n_ini
       integer,                  intent( in)    :: n_end
       type(Atom_List_Type),     intent(in out) :: fatom
       type(space_group_type),   intent(in)     :: Spg

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=132),dimension(5) :: dire
       integer                         :: i,k,nlong
       integer                         :: nop_in_line,key

       call init_err_refcodes()

       do i=n_ini,n_end
          line=adjustl(file_dat%line(i))
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle
          k=index(line,"!")
          if( k /= 0) line=trim(line(1:k-1))

          select case (l_case(line(1:4)))

             !---- Main Directive: FIX ----!
             case ("fix ", "fixe")
                call cutst(line,nlong)
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                         return       !!!!
                      case ("orie")
                         key=7
                         return       !!!!
                      case ("ther")
                         key=8
                         return       !!!!
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"fix",dire(k),fatom,spg)
                   if (err_refcodes) return
                end do

             !---- Main Directive: VARY ----!
             case ("vary")
                call cutst(line,nlong)
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                         return     !!!!
                      case ("orie")
                         key=7
                         return     !!!!
                      case ("ther")
                         key=8
                         return     !!!!
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"var",dire(k),fatom,spg)
                   if (err_refcodes) return
                end do

             !---- Main Directive: EQUAL ----!
             case ("equa")
                call cutst(line,nlong)
                call get_concodes_line(line,fatom)

             !---- Restraints Cases ----!
             case ("aequ") ! AEQU sigma        (Equal Angles restraints)

             case ("afix") ! AFIX ang sigma    (Angles restraints)
                call cutst(line,nlong)
                call get_restang_line(line,fatom)

             case ("dequ") ! DEQU sigma        (Equal Distance restraints)

             case ("dfix") ! DFIX d sigma      (Distance restraints)
                call cutst(line,nlong)
                call get_restdis_line(line,fatom)

             case ("flat") ! FLAT

             case ("tequ") ! TEQU sigma        (Equal Torsion angle restraints)

             case ("tfix") ! TFIX ang sigma    (Torsion angle restraints)
                call cutst(line,nlong)
                call get_resttor_line(line,fatom)

          end select
       end do

       return
    End Subroutine Read_RefCodes_File_FAtom

    !!--++
    !!--++ Subroutine Read_RefCodes_File_MagStr(file_dat,n_ini,n_end,FmAtom,Mag_dom)
    !!--++    Type(file_list_type),    intent( in)    :: file_dat
    !!--++    integer,                 intent( in)    :: n_ini
    !!--++    integer,                 intent( in)    :: n_end
    !!--++    type(mAtom_List_Type),   intent(in out) :: FmAtom
    !!--++    type(Magnetic_Domain_type),intent(in out) :: Mag_dom
    !!--++
    !!--++    Subroutine for treatment of magnetic Codes controls
    !!--++    magnetic clone of Read_RefCodes_File_FmAtom + Magdom
    !!--++ Created: February - 2012
    !!
    Subroutine Read_RefCodes_File_MagStr(file_dat,n_ini,n_end,FmAtom,Mag_dom)
       !---- Arguments ----!
       Type(file_list_type),              intent( in)    :: file_dat
       integer,                           intent( in)    :: n_ini
       integer,                           intent( in)    :: n_end
       type(mAtom_List_Type),             intent(in out) :: FmAtom
       type(Magnetic_Domain_type),optional,intent(in out):: Mag_dom

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=132),dimension(5) :: dire
       integer                         :: i,k,nlong
       integer                         :: nop_in_line,key

       call init_err_refcodes()

       do i=n_ini,n_end
          line=adjustl(file_dat%line(i))
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle
          k=index(line,"!")
          if( k /= 0) line=trim(line(1:k-1))

          select case (l_case(line(1:4)))

             !---- Main Directive: FIX ----!
             case ("fix ", "fixe")
                call cutst(line,nlong)
                call split_moperations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("rxyz")
                         key=1
                      case ("ixyz")
                         key=2
                      case ("mxyz")
                         key=3
                      case ("magd")
                         key=4
                      case default
                         key=0
                   end select
                   if (key /=0.and.key /=4) call cutst(dire(k),nlong)
                   if (key ==4) then
                    call get_refcodes_line(key,"fix",dire(k),Mag_dom)
                   else
                    call get_refcodes_line(key,"fix",dire(k),FmAtom)
                   end if
                   if (err_refcodes) return
                end do

             !---- Main Directive: VARY ----!
             case ("vary")
                call cutst(line,nlong)
                call split_moperations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("rxyz")
                         key=1
                      case ("ixyz")
                         key=2
                      case ("mxyz")
                         key=3
                      case ("magd")
                         key=4
                      case default
                         key=0
                   end select

                   if (key /=0.and.key /=4) call cutst(dire(k),nlong)
                   if (key ==4) then
                    call get_refcodes_line(key,"var",dire(k),Mag_dom)
                   else
                    call get_refcodes_line(key,"var",dire(k),FmAtom)
                   end if
                   if (err_refcodes) return
                end do

             !---- Main Directive: EQUAL ----! under construction
             case ("equa")
                call cutst(line,nlong)
                   if (key ==4) then
                    call get_concodes_line(line,Mag_dom)
                   else
                    call get_concodes_line(line,FmAtom)
                   end if

          end select
       end do

       return
    End Subroutine Read_RefCodes_File_MagStr

    !!--++
    !!--++ Subroutine Read_RefCodes_File_Molcrys(file_dat,n_ini,n_end,molcrys)
    !!--++    Type(file_list_type),         intent( in)    :: file_dat
    !!--++    integer,                      intent( in)    :: n_ini
    !!--++    integer,                      intent( in)    :: n_end
    !!--++    type(molecular_crystal_type), intent(in out) :: molcrys
    !!--++
    !!--++    (Overloaded)
    !!--++    Subroutine for treatment of Codes controls
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Read_RefCodes_File_Molcrys(file_dat,n_ini,n_end,molcrys)
       !---- Arguments ----!
       Type(file_list_type),         intent( in)    :: file_dat
       integer,                      intent( in)    :: n_ini
       integer,                      intent( in)    :: n_end
       type(molecular_crystal_type), intent(in out) :: molcrys

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=132),dimension(5) :: dire
       integer                         :: i,k,npos,nlong,iv
       integer                         :: nop_in_line,key,nmol
       integer, dimension(1)           :: ivet
       real(kind=cp), dimension(1)     :: vet

       call init_err_refcodes()

       do i=n_ini,n_end
          line=adjustl(file_dat%line(i))
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle
          k=index(line,"!")
          if( k /= 0) line=trim(line(1:k-1))

          nmol=-1
          select case (l_case(line(1:4)))

             !---- Main Directive: FIX ----!
             case ("fix ", "fixe") ! FIX
                call cutst(line,nlong)

                !---- Molecule Information ----!
                if (u_case(line(1:3)) =="MOL") then
                   npos=index(line," ")
                   call getnum(line(4:npos),vet,ivet,iv)
                   if (iv /= 1) then
                      nmol=0
                   else
                      nmol=ivet(1)
                   end if
                end if
                call cutst(line,nlong)

                !---- Splitting Line ----!
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                      case ("orie")
                         key=7
                      case ("ther")
                         key=8
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"fix",dire(k),molcrys,nmol)
                   if (err_refcodes) return
                end do

             !---- Main Directive: VARY ----!
             case ("vary") ! VARY
                call cutst(line,nlong)
                !---- Molecule Information ----!
                if (u_case(line(1:3)) =="MOL") then
                   npos=index(line," ")
                   call getnum(line(4:npos),vet,ivet,iv)
                   if (iv /= 1) then
                      nmol=0
                   else
                      nmol=ivet(1)
                   end if
                end if
                call cutst(line,nlong)

                !---- Splitting Line ----!
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                      case ("orie")
                         key=7
                      case ("ther")
                         key=8
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"var",dire(k),molcrys,nmol)
                   if (err_refcodes) return
                end do

             !---- Main Directive: EQUAL ----!
             case ("equa")
                call cutst(line,nlong)
                call get_concodes_line(line,molcrys)

             !---- Restraints Cases ----!
             case ("aequ") ! AEQU sigma        (Angles restraints)
             case ("afix") ! AFIX ang sigma    (Angles restraints)
             case ("dequ") ! DEQU sigma        (Distance restraints)
             case ("dfix") ! DFIX d sigma      (Distance restraints)
             case ("flat") ! FLAT
             case ("tequ") ! TEQU sigma        (Torsion angle restraints)
             case ("tfix") ! TFIX ang sigma    (Torsion angle restraints)

          end select
       end do

       return
    End Subroutine Read_RefCodes_File_Molcrys

    !!--++
    !!--++ Subroutine Read_RefCodes_File_Molec(file_dat,n_ini,n_end,molec,spg)
    !!--++    Type(file_list_type),   intent( in)    :: file_dat
    !!--++    integer,                intent( in)    :: n_ini
    !!--++    integer,                intent( in)    :: n_end
    !!--++    type(molecule_type),    intent(in out) :: molec
    !!--++    type(space_group_type), intent(in)     :: Spg
    !!--++
    !!--++    (Overloaded)
    !!--++    Subroutine for treatment of Codes controls
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Read_RefCodes_File_Molec(file_dat,n_ini,n_end,molec,spg)
       !---- Arguments ----!
       Type(file_list_type),   intent( in)    :: file_dat
       integer,                intent( in)    :: n_ini
       integer,                intent( in)    :: n_end
       type(molecule_type),    intent(in out) :: molec
       type(space_group_type), intent(in)     :: Spg

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=132),dimension(5) :: dire
       integer                         :: i,k,nlong
       integer                         :: nop_in_line,key

       call init_err_refcodes()

       do i=n_ini,n_end
          line=adjustl(file_dat%line(i))
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle
          k=index(line,"!")
          if( k /= 0) line=trim(line(1:k-1))

          select case (l_case(line(1:4)))

             !---- Main Directive: FIX ----!
             case ("fix ", "fixe")
                call cutst(line,nlong)

                !---- Splitting Line ----!
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                         return     !!!!
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                      case ("orie")
                         key=7
                      case ("ther")
                         key=8
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"fix",dire(k),molec,spg)
                   if (err_refcodes) return
                end do

             !---- Main Directive: VARY ----!
             case ("vary")
                call cutst(line,nlong)

                !---- Splitting Line ----!
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                         return     !!!!
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                      case ("orie")
                         key=7
                      case ("ther")
                         key=8
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"var",dire(k),molec,spg)
                   if (err_refcodes) return
                end do

             !---- Main Directive: EQUAL ----!
             case ("equa")
                call cutst(line,nlong)
                call get_concodes_line(line,molec)

             !---- Restraints Cases ----!
             case ("aequ") ! AEQU sigma        (Angles restraints)
             case ("afix") ! AFIX ang sigma    (Angles restraints)
             case ("dequ") ! DEQU sigma        (Distance restraints)
             case ("dfix") ! DFIX d sigma      (Distance restraints)
             case ("flat") ! FLAT
             case ("tequ") ! TEQU sigma        (Torsion angle restraints)
             case ("tfix") ! TFIX ang sigma    (Torsion angle restraints)

          end select
       end do

       return
    End Subroutine Read_RefCodes_File_Molec

    !!----
    !!---- Subroutine Read_RefGCodes_File(file_dat,n_ini,n_end,model,sys,Iphas)
    !!----    Type(file_list_type),                intent( in)    :: file_dat
    !!----    integer,                             intent( in)    :: n_ini
    !!----    integer,                             intent( in)    :: n_end
    !!----    type(Nonatomic_Parameter_List_Type), intent(in out) :: model
    !!----    character(len=*), optional,          intent( in)    :: sys
    !!----    integer, optional,                   intent( in)    :: Iphas
    !!----
    !!----    Subroutine for treatment of Codes for non-atomic parameters.
    !!----    The optional argument Sys containts the crystal system and the
    !!----    setting in case of the monoclinic, e.g. Sys="Monoclinic b".
    !!----
    !!---- Update: November - 2013
    !!
    Subroutine Read_RefGCodes_File(file_dat,n_ini,n_end,model,sys,Iphas)
       Type(file_list_type),                intent( in)    :: file_dat
       integer,                             intent( in)    :: n_ini
       integer,                             intent( in)    :: n_end
       type(Nonatomic_Parameter_List_Type), intent(in out) :: model
       character(len=*), optional,          intent( in)    :: sys
       integer, optional,                   intent( in)    :: Iphas
       !---- Local variables ----!
       character(len=132)              :: line
       character(len=20)               :: namp
       character(len=132),dimension(5) :: dire
       integer                         :: i,k,nlong
       integer                         :: nop_in_line,key

       call init_err_refcodes()

       do i=n_ini,n_end
          line=adjustl(file_dat%line(i))
          if (line(1:1) =="!" .or. len_trim(line) == 0) cycle
          namp=" "
          k=index(line,"!")
          if( k /= 0) line=trim(line(1:k-1))

          select case (l_case(line(1:5)))

             !---- Main Directive: GFIX ----!
             case ("gfix ")
                call cutst(line,nlong)
                call split_goperations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   namp=trim(l_case(dire(k)))  !keep the name of the parameter in the directive
                   select case (trim(namp))    !Simple general directives key /= 0
                      case ("bkg")
                         key=1
                      case ("cell")
                         key=2
                      case ("uvw")
                         key=3
                      case ("asize")
                         key=4
                      case ("astrain")
                         key=5
                      case ("extinct")
                         key=6
                      case ("scalefs")
                         key=7
                      case ("all")
                         key=8
                      case default  !other parameter names or the directive contains numbers or other options
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)  !in case key=0  dire(k)=namp, otherwise dire(k) does not contain the generic name of the parameter
                   if(present(iphas)) then
                     call get_refGcodes_line(key,"gfix",dire(k),namp,model,Iphas=iphas)
                   else
                     call get_refGcodes_line(key,"gfix",dire(k),namp,model)
                   end if
                   if (err_refcodes) return
                end do

             !---- Main Directive: GVARY ----!
             case ("gvary")
                call cutst(line,nlong)
                call split_goperations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   namp=trim(l_case(dire(k)))  !keep the name of the parameter in the directive
                   select case (trim(namp))    !Simple general directives key /= 0
                      case ("bkg")
                         key=1
                      case ("cell")
                         key=2
                      case ("uvw")
                         key=3
                      case ("asize")
                         key=4
                      case ("astrain")
                         key=5
                      case ("extinct")
                         key=6
                      case ("scalefs")
                         key=7
                      case ("all")
                         key=8
                      case default  !other parameter names or the directive contains numbers or other options
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong) !in case key=0  dire(k)=namp,otherwise dire(k) does not contain the generic name of the parameter
                   if(present(Sys)) then
                      if(present(iphas)) then
                        call get_refGcodes_line(key,"gvar",dire(k),namp,model,sys,Iphas=iphas)
                      else
                        call get_refGcodes_line(key,"gvar",dire(k),namp,model,sys)
                      end if
                   else
                      if(present(iphas)) then
                        call get_refGcodes_line(key,"gvar",dire(k),namp,model,Iphas=iphas)
                      else
                        call get_refGcodes_line(key,"gvar",dire(k),namp,model)
                      end if
                   end if
                   if (err_refcodes) return
                end do

             !---- Main Directive: EQUAL ----!
              !case ("equa")
              !   call cutst(line,nlong)
              !   call get_congcodes_line(line,model)

          end select
       end do

       return
    End Subroutine Read_RefGCodes_File

    !!--++
    !!--++ Subroutine Split_Operations(Line, Ni, S_Lines)
    !!--++    character(len=*),              intent( in) :: line
    !!--++    integer,                       intent(out) :: ni
    !!--++    character(len=*),dimension(:), intent(out) :: s_lines
    !!--++
    !!--++    (Private)
    !!--++    Split diferent directives according to KEY_CODE Variable
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Split_Operations(Line, Ni, S_Lines)
       !---- Arguments ----!
       character(len=*),              intent( in) :: line
       integer,                       intent(out) :: ni
       character(len=*),dimension(:), intent(out) :: s_lines

       !---- Local variables ----!
       character(len=150)      :: linec
       character(len=10)       :: car
       integer                 :: i,j,npos
       integer,dimension(nkey) :: ip,ipc

       ni=0
       s_lines=" "

       linec=u_case(line)

       !---- Search Subkeys: XYZ,OCC,BIS... ----!
       ip=0
       do i=1,nkey
          npos=index(linec,key_code(i))
          if (npos > 0) then
             car=linec(npos:)
             j=index(car," ")
             if (j > 0) car=car(:j-1)
             j=index(car,"_")
             if (j > 0) cycle
          end if
          ip(i)=npos
       end do

       npos=count(ip > 0)
       if (npos == 0) then
          ni=1
          s_lines(1)=adjustl(line)
       else
          call sort(ip,nkey,ipc)
          do i=1,nkey
             if (ip(ipc(i)) == 0) then
                if (ip(ipc(i+1)) <= 1) cycle
                ni=ni+1
                s_lines(ni)=adjustl(line(1:ip(ipc(i+1))-1))
             else
                ni=ni+1
                if (i < nkey) then
                   s_lines(ni)=adjustl(line(ip(ipc(i)):ip(ipc(i+1))-1))
                else
                   s_lines(ni)=adjustl(line(ip(ipc(i)):))
                end if
             end if
          end do
       end if

       return
    End Subroutine Split_Operations
    !!----
    !!---- Subroutine Split_GOperations(Line, Ni, S_Lines)
    !!----    character(len=*),              intent( in) :: line
    !!----    integer,                       intent(out) :: ni
    !!----    character(len=*),dimension(:), intent(out) :: s_lines
    !!----
    !!----    (Private)
    !!----    Split diferent directives according to GCODE_NAM Variable
    !!----
    !!---- Update: November - 2013
    !!
    Subroutine Split_GOperations(Line, Ni, S_Lines)
       !---- Arguments ----!
       character(len=*),              intent( in) :: line
       integer,                       intent(out) :: ni
       character(len=*),dimension(:), intent(out) :: s_lines

       !---- Local variables ----!
       character(len=150)        :: linec
       character(len=10)         :: up_gcode_nam
       integer                   :: i,npos
       integer,dimension(NGCode) :: ip,ipc

       ni=0
       s_lines=" "

       linec=u_case(line)

       !---- Search Subkeys: scale,cell,delta... ----!
       ip=0
       do i=1,NGCode
          up_gcode_nam=u_case(Gcode_Nam(i))
          ip(i)=index(linec,up_gcode_nam)
       end do

       npos=count(ip > 0)
       if (npos == 0) then
          ni=1
          s_lines(1)=adjustl(line)
       else
          call sort(ip,NGCode,ipc)
          do i=1,NGCode
             if (ip(ipc(i)) == 0) then
                if (ip(ipc(i+1)) <= 1) cycle
                ni=ni+1
                s_lines(ni)=adjustl(line(1:ip(ipc(i+1))-1))
             else
                ni=ni+1
                if (i < NGCode) then
                   s_lines(ni)=adjustl(line(ip(ipc(i)):ip(ipc(i+1))-1))
                else
                   s_lines(ni)=adjustl(line(ip(ipc(i)):))
                end if
             end if
          end do
       end if

       return
    End Subroutine Split_GOperations
    !!--++
    !!--++ Subroutine Split_mOperations(Line, Ni, S_Lines)
    !!--++    character(len=*),              intent( in) :: line
    !!--++    integer,                       intent(out) :: ni
    !!--++    character(len=*),dimension(:), intent(out) :: s_lines
    !!--++
    !!--++    (Private)
    !!--++    Split diferent directives according to KEY_mCODE Variable
    !!--++    magnetic clone of Subroutine Split_Operations
    !!--++ Created: December - 2011
    !!
    Subroutine Split_mOperations(Line, Ni, S_Lines)
       !---- Arguments ----!
       character(len=*),              intent( in) :: line
       integer,                       intent(out) :: ni
       character(len=*),dimension(:), intent(out) :: s_lines

       !---- Local variables ----!
       character(len=150)      :: linec
       character(len=10)       :: car
       integer                 :: i,j,npos
       integer,dimension(mnkey) :: ip,ipc

       ni=0
       s_lines=" "

       linec=u_case(line)

       !---- Search Subkeys: Rxyz,Ixyz,Mxyz ----!
       ip=0
       do i=1,mnkey
          npos=index(linec,key_mcode(i))
          if (npos > 0) then
             car=linec(npos:)
             j=index(car," ")
             if (j > 0) car=car(:j-1)
             j=index(car,"_")
             if (j > 0) cycle
          end if
          ip(i)=npos
       end do

       npos=count(ip > 0)
       if (npos == 0) then
          ni=1
          s_lines(1)=adjustl(line)
       else
          call sort(ip,mnkey,ipc)
          do i=1,mnkey
             if (ip(ipc(i)) == 0) then
                if (ip(ipc(i+1)) <= 1) cycle
                ni=ni+1
                s_lines(ni)=adjustl(line(1:ip(ipc(i+1))-1))
             else
                ni=ni+1
                if (i < mnkey) then
                   s_lines(ni)=adjustl(line(ip(ipc(i)):ip(ipc(i+1))-1))
                else
                   s_lines(ni)=adjustl(line(ip(ipc(i)):))
                end if
             end if
          end do
       end if

       return
    End Subroutine Split_mOperations

    !!----
    !!---- Subroutine VState_to_AtomsPar(FAtom/MolCrys/Molec/Magdom,Mode)
    !!----    type(Atom_List_Type),         intent(in out) :: FAtom
    !!----    or
    !!----    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!----    or
    !!----    type(molecule_type),          intent(in out) :: Molec
    !!----    character(len=*), optional,   intent(in)     :: Mode
    !!----
    !!----    Update the values to the variable FAtom/MolCrys/Molec from Vector
    !!----
    !!---- Update: November 22 - 2013
    !!

    !!--++
    !!--++ Subroutine VState_to_AtomsPar_FAtom(FAtom,Mode)
    !!--++    type(Atom_List_Type),       intent(in out) :: FAtom
    !!--++    character(len=*), optional, intent(in)     :: Mode
    !!--++
    !!--++    (Overloaded)
    !!--++
    !!--++ Update: November 22 - 2013.
    !!--++ Modified to include standard deviations, November 3, 2013 (JRC)
    !!
    Subroutine VState_to_AtomsPar_FAtom(FAtom,Mode)
       !---- Arguments ----!
       type(Atom_List_Type),       intent(in out) :: FAtom
       character(len=*), optional, intent(in)     :: Mode

       !---- Local Variables ----!
       integer          :: i,j,l
       character(len=1) :: car

       call init_err_refcodes()

       car="s"
       if (present(mode)) car=adjustl(mode)

       do i=1,FAtom%natoms
          !---- XYZ ----!
          do j=1,3
             l=FAtom%atom(i)%lx(j)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   FAtom%atom(i)%x(j)=v_vec(l)*FAtom%atom(i)%mx(j)

                case ("s","S") ! Passing Shift
                   FAtom%atom(i)%x(j)=FAtom%atom(i)%x(j)+v_shift(l)*FAtom%atom(i)%mx(j)
             end select
             FAtom%atom(i)%x_std(j)=v_vec_std(l)*FAtom%atom(i)%mx(j)
          end do

          !---- BISO ----!
          l=FAtom%atom(i)%lbiso
          if (l > 0) then
              if (l > np_refi) then
                  err_refcodes=.true.
                  ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                  return
              end if
              select case (car)
                 case ("v","V") ! Passing Value
                    FAtom%atom(i)%biso=v_vec(l)*FAtom%atom(i)%mbiso

                 case ("s","S") ! Passing Shift
                    FAtom%atom(i)%biso=FAtom%atom(i)%biso+v_shift(l)*FAtom%atom(i)%mbiso
              end select
              FAtom%atom(i)%biso_std=v_vec_std(l)*FAtom%atom(i)%mbiso
          end if


          !---- OCC ----!
          l=FAtom%atom(i)%locc
          if (l > 0) then
              if (l > np_refi) then
                 err_refcodes=.true.
                 ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                 return
              end if
              select case (car)
                 case ("v","V") ! Passing Value
                    FAtom%atom(i)%occ=v_vec(l)*FAtom%atom(i)%mocc
                 case ("s","S") ! Passing Shift
                    FAtom%atom(i)%occ=FAtom%atom(i)%occ+v_shift(l)*FAtom%atom(i)%mocc

              end select
              FAtom%atom(i)%occ_std=v_vec_std(l)*FAtom%atom(i)%mocc
          end if


          !---- BANIS ----!
          do j=1,6
             l=FAtom%atom(i)%lu(j)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   FAtom%atom(i)%u(j)=v_vec(l)*FAtom%atom(i)%mu(j)

                case ("s","S") ! Passing Shift
                   FAtom%atom(i)%u(j)=FAtom%atom(i)%u(j)+v_shift(l)*FAtom%atom(i)%mu(j)
             end select
             FAtom%atom(i)%u_std(j)=v_vec_std(l)*FAtom%atom(i)%mu(j)
          end do

       end do

       return
    End Subroutine VState_to_AtomsPar_FAtom

    !!--++
    !!--++ Subroutine VState_to_AtomsPar_FmAtom(FmAtom,Mode,MGp,Mag_dom)
    !!--++    type(mAtom_List_Type),      intent(in out)           :: FmAtom
    !!--++    character(len=*), optional, intent(in)               :: Mode
    !!--++    type(MagSymm_k_Type), optional, intent(in)           :: MGp
    !!--++    type(Magnetic_Domain_type), optional, intent(in out) :: Mag_dom
    !!--++
    !!--++ magnetic clone of VState_to_AtomsPar_FAtom
    !!--++ Created: December - 2011
    !!--++ Updated: February - 2012, November 3, 2013 (standard deviations,JRC)
    !!
    Subroutine VState_to_AtomsPar_FmAtom(FmAtom,Mode,MGp,Mag_dom)
       !---- Arguments ----!
       type(mAtom_List_Type),                intent(in out) :: FmAtom
       character(len=*),           optional, intent(in)     :: Mode
       type(MagSymm_k_Type),       optional, intent(in)     :: MGp
       type(Magnetic_Domain_type), optional, intent(in out) :: Mag_dom
       !---- Local Variables ----!
       integer          :: i,j,l,ik,ich
       character(len=1) :: car

       call init_err_refcodes()
       car="s"
       if (present(mode)) car=adjustl(mode)

       do i=1,FmAtom%natoms
          ik=FmAtom%atom(i)%nvk

          !---- Rxyz ----!
          do j=1,3
             l=FmAtom%atom(i)%lSkR(j,ik)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if

             select case (car)

                case ("v","V") ! Passing Value
                 if(MGp%Sk_type == "Spherical_Frame") then
                   FmAtom%atom(i)%Spher_SkR(j,ik)=v_vec(l)*FmAtom%atom(i)%mSkR(j,ik)
                   FmAtom%atom(i)%Spher_SkR_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkR(j,ik)
                 else
                   FmAtom%atom(i)%SkR(j,ik)=v_vec(l)*FmAtom%atom(i)%mSkR(j,ik)
                   FmAtom%atom(i)%SkR_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkR(j,ik)
                 end if

                case ("s","S") ! Passing Shift
                 if(MGp%Sk_type == "Spherical_Frame") then
                   FmAtom%atom(i)%Spher_SkR(j,ik)=v_vec(l)*FmAtom%atom(i)%mSkR(j,ik)
                   FmAtom%atom(i)%Spher_SkR_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkR(j,ik)
                 else
                   FmAtom%atom(i)%SkR(j,ik)=FmAtom%atom(i)%SkR(j,ik)+v_shift(l)*FmAtom%atom(i)%mSkR(j,ik)
                   FmAtom%atom(i)%SkR_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkR(j,ik)
                 end if

             end select
          end do

          !---- Ixyz ----!
          do j=1,3
             l=FmAtom%atom(i)%lSkI(j,ik)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                 if(MGp%Sk_type == "Spherical_Frame") then
                   FmAtom%atom(i)%Spher_SkI(j,ik)=v_vec(l)*FmAtom%atom(i)%mSkI(j,ik)
                   FmAtom%atom(i)%Spher_SkI_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkI(j,ik)
                 else
                   FmAtom%atom(i)%SkI(j,ik)=v_vec(l)*FmAtom%atom(i)%mSkI(j,ik)
                   FmAtom%atom(i)%SkI_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkI(j,ik)
                 end if

                case ("s","S") ! Passing Shift
                 if(MGp%Sk_type == "Spherical_Frame") then
                   FmAtom%atom(i)%Spher_SkI(j,ik)=v_vec(l)*FmAtom%atom(i)%mSkI(j,ik)
                   FmAtom%atom(i)%Spher_SkI_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkI(j,ik)
                 else
                   FmAtom%atom(i)%SkI(j,ik)=FmAtom%atom(i)%SkI(j,ik)+v_shift(l)*FmAtom%atom(i)%mSkI(j,ik)
                   FmAtom%atom(i)%SkI_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkI(j,ik)
                 end if

             end select
          end do

          !---- Mxyz ----!
          do j=1,6
            if(1 <= j .and. j <= 3) l=FmAtom%atom(i)%lSkR(j,ik)
            if(4 <= j .and. j <= 6) l=FmAtom%atom(i)%lSkI(j-3,ik)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                 if(MGp%Sk_type == "Spherical_Frame") then
                  if(1 <= j .and. j <= 3) then
                    FmAtom%atom(i)%Spher_SkR(j,ik)=v_vec(l)*FmAtom%atom(i)%mSkR(j,ik)
                    FmAtom%atom(i)%Spher_SkR_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkR(j,ik)
                  end if
                  if(4 <= j .and. j <= 6) then
                    FmAtom%atom(i)%Spher_SkI(j-3,ik)=v_vec(l)*FmAtom%atom(i)%mSkI(j-3,ik)
                    FmAtom%atom(i)%Spher_SkI_std(j-3,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkI(j-3,ik)
                  end if
                 else
                  if(1 <= j .and. j <= 3) then
                    FmAtom%atom(i)%SkR(j,ik)=v_vec(l)*FmAtom%atom(i)%mSkR(j,ik)
                    FmAtom%atom(i)%SkR_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkR(j,ik)
                  end if
                  if(4 <= j .and. j <= 6) then
                    FmAtom%atom(i)%SkI(j-3,ik)=v_vec(l)*FmAtom%atom(i)%mSkI(j-3,ik)
                    FmAtom%atom(i)%SkI_std(j-3,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkI(j-3,ik)
                  end if
                 end if

                case ("s","S") ! Passing Shift
                 if(MGp%Sk_type == "Spherical_Frame") then
                  if(1 <= j .and. j <= 3) then
                    FmAtom%atom(i)%Spher_SkR(j,ik)=FmAtom%atom(i)%Spher_SkR(j,ik)+ &
                                               v_shift(l)*FmAtom%atom(i)%mSkR(j,ik)
                    FmAtom%atom(i)%Spher_SkR_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkR(j,ik)
                  end if
                  if(4 <= j .and. j <= 6) then
                    FmAtom%atom(i)%Spher_SkI(j-3,ik)=FmAtom%atom(i)%Spher_SkI(j-3,ik)+ &
                                                v_shift(l)*FmAtom%atom(i)%mSkI(j-3,ik)
                    FmAtom%atom(i)%Spher_SkI_std(j-3,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkI(j-3,ik)
                  end if
                 else
                  if(1 <= j .and. j <= 3) then
                    FmAtom%atom(i)%SkR(j,ik)=FmAtom%atom(i)%SkR(j,ik)+ &
                                             v_shift(l)*FmAtom%atom(i)%mSkR(j,ik)
                    FmAtom%atom(i)%SkR_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkR(j,ik)
                  end if
                  if(4 <= j .and. j <= 6) then
                    FmAtom%atom(i)%SkI(j-3,ik)=FmAtom%atom(i)%SkI(j-3,ik)+ &
                                               v_shift(l)*FmAtom%atom(i)%mSkI(j-3,ik)
                    FmAtom%atom(i)%SkI_std(j-3,ik)=v_vec_std(l)*FmAtom%atom(i)%mSkI(j-3,ik)
                  end if
                 end if
             end select
          end do

          !---- MagPh ----!
          l=FmAtom%atom(i)%lmphas(ik)
          if (l /= 0) then
            if (l > np_refi) then
               err_refcodes=.true.
               ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
               return
            end if
            select case (car)
               case ("v","V") ! Passing Value
                  FmAtom%atom(i)%mphas(ik)=v_vec(l)*FmAtom%atom(i)%mmphas(ik)

               case ("s","S") ! Passing Shift
                  FmAtom%atom(i)%mphas(ik)=FmAtom%atom(i)%mphas(ik)+v_shift(l)*FmAtom%atom(i)%mmphas(ik)
            end select
            FmAtom%atom(i)%mphas_std(ik)=v_vec_std(l)*FmAtom%atom(i)%mmphas(ik)
          end if
          !---- C1-12 ----!
          do j=1,12
             l=FmAtom%atom(i)%lbas(j,ik)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   FmAtom%atom(i)%cbas(j,ik)=v_vec(l)*FmAtom%atom(i)%mbas(j,ik)

                case ("s","S") ! Passing Shift
                   FmAtom%atom(i)%cbas(j,ik)=FmAtom%atom(i)%cbas(j,ik)+v_shift(l)*FmAtom%atom(i)%mbas(j,ik)
             end select
             FmAtom%atom(i)%cbas_std(j,ik)=v_vec_std(l)*FmAtom%atom(i)%mbas(j,ik)
          end do
       end do !on atoms

       !---- Check is chirality is present ----!
       if (Mag_Dom%chir) then
        ich=2
       else
        ich=1
       end if

       do i=1,Mag_Dom%nd
         do j=1,ich
             l=Mag_Dom%Lpop(j,i)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   Mag_Dom%pop(j,i)=v_vec(l)*Mag_Dom%Mpop(j,i)

                case ("s","S") ! Passing Shift
                   Mag_Dom%pop(j,i)=Mag_Dom%pop(j,i)+v_shift(l)*Mag_Dom%Mpop(j,i)
             end select
             Mag_Dom%pop_std(j,i)=v_vec_std(l)*Mag_Dom%Mpop(j,i)
         end do
        end do

       return
    End Subroutine VState_to_AtomsPar_FmAtom

    !!--++
    !!--++ Subroutine VState_to_AtomsPar_Molcrys(Molcrys,Mode)
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++    character(len=*), optional,   intent(in)     :: Mode
    !!--++
    !!--++    (Overloaded)
    !!--++
    !!--++ Update: November 22 - 2013
    !!
    Subroutine VState_to_AtomsPar_Molcrys(Molcrys,Mode)
       !---- Arguments ----!
       type(molecular_Crystal_type), intent(in out) :: MolCrys
       character(len=*), optional,    intent(in)     :: Mode

       !---- Local variables ----!
       integer          :: i,j,jj,k,l
       character(len=1) :: car

       call init_err_refcodes()

       car="s"
       if (present(mode)) car=adjustl(mode)

       do i=1,molcrys%n_free
          !---- XYZ ----!
          do j=1,3
             l=molcrys%atm(i)%lx(j)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molcrys%atm(i)%x(j)=v_vec(l)*molcrys%atm(i)%mx(j)

                case ("s","S") ! Passing Shift
                   molcrys%atm(i)%x(j)=molcrys%atm(i)%x(j)+v_shift(l)*molcrys%atm(i)%mx(j)
             end select
          end do

          !---- BISO ----!
          l=molcrys%atm(i)%lbiso
          if (l > 0) then
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molcrys%atm(i)%biso=v_vec(l)*molcrys%atm(i)%mbiso

                case ("s","S") ! Passing Shift
                   molcrys%atm(i)%biso=molcrys%atm(i)%biso+v_shift(l)*molcrys%atm(i)%mbiso
             end select
          end if

          !---- OCC ----!
          l=molcrys%atm(i)%locc
          if (l > 0) then
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molcrys%atm(i)%occ=v_vec(l)*molcrys%atm(i)%mocc

                case ("s","S") ! Passing Shift
                   molcrys%atm(i)%occ=molcrys%atm(i)%occ+v_shift(l)*molcrys%atm(i)%mocc
             end select
          end if

          !---- BANIS ----!
          do j=1,6
             l=molcrys%atm(i)%lu(i)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molcrys%atm(i)%u(j)=v_vec(l)*molcrys%atm(i)%mu(j)

                case ("s","S") ! Passing Shift
                   molcrys%atm(i)%u(j)=molcrys%atm(i)%u(j)+v_shift(l)*molcrys%atm(i)%mu(j)
             end select
          end do
       end do

       do k=1,molcrys%n_mol
          do i=1,molcrys%mol(k)%natoms
             !---- Coordinates ----!
             do j=1,3
                l=molcrys%mol(k)%lI_Coor(j,i)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(k)%I_Coor(j,i)=v_vec(l)*molcrys%mol(k)%mI_Coor(j,i)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(k)%I_Coor(j,i)=molcrys%mol(k)%I_Coor(j,i)+v_shift(l)*molcrys%mol(k)%mI_Coor(j,i)
                end select
             end do

             !---- Biso ----!
             l=molcrys%mol(k)%lbiso(i)
             if (l > 0) then
                if (l > np_refi) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(k)%biso(i)=v_vec(l)*molcrys%mol(k)%mbiso(i)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(k)%biso(i)=molcrys%mol(k)%biso(i)+v_shift(l)*molcrys%mol(k)%mbiso(i)
                end select
             end if

             !---- Occ ----!
             l=molcrys%mol(k)%locc(i)
             if (l > 0) then
                if (l > np_refi) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(k)%occ(i)=v_vec(l)*molcrys%mol(k)%mocc(i)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(k)%occ(i)=molcrys%mol(k)%occ(i)+v_shift(l)*molcrys%mol(k)%mocc(i)
                end select
             end if

             !---- Centre ----!
             do j=1,3
                l=molcrys%mol(i)%lxcentre(j)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(i)%xcentre(j)=v_vec(l)*molcrys%mol(i)%mxcentre(j)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(i)%xcentre(j)=molcrys%mol(i)%xcentre(j)+v_shift(l)*molcrys%mol(i)%mxcentre(j)
                end select
             end do

             !---- Orient ----!
             do j=1,3
                l=molcrys%mol(i)%lorient(j)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(i)%orient(j)=v_vec(l)*molcrys%mol(i)%morient(j)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(i)%orient(j)=molcrys%mol(i)%orient(j)+v_shift(l)*molcrys%mol(i)%morient(j)
                end select
             end do

             !---- T_TLS ----!
             do j=1,6
                l=molcrys%mol(i)%lT_TLS(j)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(i)%T_TLS(j)=v_vec(l)*molcrys%mol(i)%mT_TLS(j)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(i)%T_TLS(j)=molcrys%mol(i)%T_TLS(j)+v_shift(l)*molcrys%mol(i)%mT_TLS(j)
                end select
             end do

             !---- L_TLS ----!
             do j=1,6
                l=molcrys%mol(i)%lL_TLS(j)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(i)%L_TLS(j)=v_vec(l)*molcrys%mol(i)%mL_TLS(j)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(i)%L_TLS(j)=molcrys%mol(i)%L_TLS(j)+v_shift(l)*molcrys%mol(i)%mL_TLS(j)
                end select
             end do

             !---- S_TLS ----!
             do j=1,3
                do jj=1,3
                   l=molcrys%mol(i)%lS_TLS(j,jj)
                   if (l == 0) cycle
                   if (l > np_refi) then
                      err_refcodes=.true.
                      ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                      return
                   end if
                   select case (car)
                      case ("v","V") ! Passing Value
                         molcrys%mol(i)%S_TLS(j,jj)=v_vec(l)*molcrys%mol(i)%mS_TLS(j,jj)

                      case ("s","S") ! Passing Shift
                         molcrys%mol(i)%S_TLS(j,jj)=molcrys%mol(i)%S_TLS(j,jj)+v_shift(l)*molcrys%mol(i)%mS_TLS(j,jj)
                   end select
                end do
             end do
          end do
       end do

       return
    End Subroutine VState_to_AtomsPar_Molcrys

    !!--++
    !!--++ Subroutine VState_to_AtomsPar_Molec(Molec,Mode)
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    character(len=*), optional,   intent(in)     :: Mode
    !!--++
    !!--++    (Overloaded)
    !!--++
    !!--++ Update: November 22 - 2013
    !!
    Subroutine VState_to_AtomsPar_Molec(Molec,Mode)
       !---- Arguments ----!
       type(molecule_type),          intent(in out) :: Molec
       character(len=*), optional,   intent(in)     :: Mode

       !---- Local variables ----!
       integer          :: i,j,l
       character(len=1) :: car

       call init_err_refcodes()

       car="s"
       if (present(mode)) car=adjustl(mode)

       do i=1,molec%natoms
          !---- Coordinates ----!
          do j=1,3
             l=molec%lI_Coor(j,i)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molec%I_Coor(j,i)=v_vec(l)*molec%mI_Coor(j,i)

                case ("s","S") ! Passing Shift
                   molec%I_Coor(j,i)=molec%I_Coor(j,i)+v_shift(l)*molec%mI_Coor(j,i)
             end select
          end do

          !---- Biso ----!
          l=molec%lbiso(i)
          if (l > 0) then
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                 case ("v","V") ! Passing Value
                    molec%biso(i)=v_vec(l)*molec%mbiso(i)

                 case ("s","S") ! Passing Shift
                    molec%biso(i)=molec%biso(i)+v_shift(l)*molec%mbiso(i)
             end select
          end if

          !---- Occ ----!
          l=molec%locc(i)
          if (l > 0) then
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                 case ("v","V") ! Passing Value
                    molec%occ(i)=v_vec(l)*molec%mocc(i)

                 case ("s","S") ! Passing Shift
                    molec%occ(i)=molec%occ(i)+v_shift(l)*molec%mocc(i)
             end select
          end if
       end do

       !---- Centre ----!
       do j=1,3
          l=molec%lxcentre(j)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%xcentre(j)=v_vec(l)*molec%mxcentre(j)

              case ("s","S") ! Passing Shift
                 molec%xcentre(j)=molec%xcentre(j)+v_shift(l)*molec%mxcentre(j)
          end select
       end do

       !---- Orient ----!
       do j=1,3
          l=molec%lorient(j)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%orient(j)=v_vec(l)*molec%morient(j)

              case ("s","S") ! Passing Shift
                 molec%orient(j)=molec%orient(j)+v_shift(l)*molec%morient(j)
          end select
       end do

       !---- T_TLS ----!
       do j=1,6
          l=molec%lT_TLS(j)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%T_TLS(j)=v_vec(l)*molec%mT_TLS(j)

              case ("s","S") ! Passing Shift
                 molec%T_TLS(j)=molec%T_TLS(j)+v_shift(l)*molec%mT_TLS(j)
          end select
       end do

       !---- L_TLS ----!
       do j=1,6
          l=molec%lL_TLS(j)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%L_TLS(j)=v_vec(l)*molec%mL_TLS(j)

              case ("s","S") ! Passing Shift
                 molec%L_TLS(j)=molec%L_TLS(j)+v_shift(l)*molec%mL_TLS(j)
          end select
       end do

       !---- S_TLS ----!
       do i=1,3
          do j=1,3
             l=molec%lS_TLS(i,j)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                 case ("v","V") ! Passing Value
                    molec%S_TLS(i,j)=v_vec(l)*molec%mS_TLS(i,j)

                 case ("s","S") ! Passing Shift
                    molec%S_TLS(i,j)=molec%S_TLS(i,j)+v_shift(l)*molec%mS_TLS(i,j)
             end select
          end do
       end do

       return
    End Subroutine VState_to_AtomsPar_Molec

    !!----
    !!---- Subroutine VState_to_ModelPar(Model,Mode)
    !!----    type(Nonatomic_Parameter_List_Type), intent(in out) :: model
    !!----    character(len=*), optional,          intent(in)     :: Mode
    !!----
    !!----    (Overloaded)
    !!----
    !!---- Update: November 2 - 2013
    !!
    Subroutine VState_to_ModelPar(Model,Mode)
       !---- Arguments ----!
       type(Nonatomic_Parameter_List_Type), intent(in out) :: model
       character(len=*), optional,          intent(in)     :: Mode

       !---- Local variables ----!
       integer          :: i,l
       character(len=1) :: car

       call init_err_refcodes()

       car="s"
       if (present(mode)) car=adjustl(mode)

       do i=1,model%npar
          l=model%par(i)%Lcode
          !write(*,"(a,i5,a,i5)") " Parameter: "//trim(model%par(i)%nam),i,"  Code-number:",l
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             ERR_RefCodes_Mess="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
             case ("v","V") ! Passing Value
                model%par(i)%value=v_vec(l)*model%par(i)%multip

             case ("s","S") ! Passing Shift
                model%par(i)%value=model%par(i)%value+v_shift(l)*model%par(i)%multip
          end select
          model%par(i)%sigma=v_vec_std(l)*model%par(i)%multip
          !write(*,"(a,2f14.5)") " New value and sigma: "//trim(model%par(i)%nam),model%par(i)%value,model%par(i)%sigma
       end do
       return
    End Subroutine VState_to_ModelPar

    !!----
    !!---- Subroutine Write_Info_RefCodes(FAtom/FmAtom/MolCrys/Molec/MagStr, Iunit)
    !!----    type(Atom_List_Type),         intent(in) :: FAtom
    !!----    or
    !!----    type(mAtom_List_Type),        intent(in) :: FmAtom
    !!----    or
    !!----    type(molecular_crystal_type), intent(in) :: Molcrys
    !!----    or
    !!----    type(molecule_type),          intent(in) :: Molec
    !!----    or
    !!--++    type(mAtom_List_Type),               intent(in) :: FmAtom
    !!--++    and type(Magnetic_Domain_type),optional, intent(in) :: Mag_dom
    !!----    integer, optional,            intent(in) :: Iunit
    !!----
    !!----    Write the Information about Refinement Codes
    !!----
    !!---- Update: March - 2005
    !!---- Update: February - 2012
    !!

    !!--++
    !!--++ Subroutine Write_Info_RefCodes_FAtom(FAtom, Spg, Iunit)
    !!--++    type(Atom_List_Type),  intent(in) :: FAtom
    !!--++    type(Space_Group_Type), intent(in) :: Spg
    !!--++    integer, optional,     intent(in) :: Iunit
    !!--++
    !!--++    Overloaded
    !!--++    Write the Information about Refinement Codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Write_Info_RefCodes_FAtom(FAtom, Spg, Iunit)
       !---- Arguments ----!
       type(Atom_List_Type),   intent(in) :: FAtom
       type(Space_Group_Type), intent(in) :: Spg
       integer, optional,      intent(in) :: Iunit

       !---- Local variables ----!
       character(len=20)              :: car
       character(len=60)              :: fmt1,fmt2,fmt3,fmt4,fmt5
       Character(len=25),dimension(3) :: symcar
       integer                        :: i,j,k,n,na,np,lun,p1,p2,p3,p4
       real(kind=cp)                  :: mu
       real(kind=cp),dimension(3)     :: tr

       !---- Format Zone ----!
       fmt1="(t5,a,t16,i3,t27,a,t33,4(tr6,f8.4),tr8,i2,tr6,f8.3,i9)"
       fmt2="(t10,i3,t21,a,t31,a,t39,f8.4)"

       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Refinable Parameters: ",np_refi
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,a)")"    Name      N.Param   Code-Name       Value        L.Bound       ",&
                                      "U.Bound         Step    Condition   Multiplier    N.Atom"
          write(unit=lun, fmt="(a,a)")" ==================================================================",&
                                      "========================================================="
          do i=1,FAtom%natoms
             do j=1,3
                if (FAtom%atom(i)%lx(j) > 0) then
                   na=FAtom%atom(i)%lx(j)
                   mu=FAtom%atom(i)%mx(j)
                   car=trim(code_nam(j))//trim(FAtom%atom(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             if (FAtom%atom(i)%locc > 0) then
                na=FAtom%atom(i)%locc
                mu=FAtom%atom(i)%mocc
                car=trim(code_nam(5))//trim(FAtom%atom(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             if (FAtom%atom(i)%lbiso > 0) then
                na=FAtom%atom(i)%lbiso
                mu=FAtom%atom(i)%mbiso
                car=trim(code_nam(4))//trim(FAtom%atom(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             do j=1,6
                if (FAtom%atom(i)%lu(j) > 0) then
                   na=FAtom%atom(i)%lu(j)
                   mu=FAtom%atom(i)%mu(j)
                   car=trim(code_nam(5+j))//trim(FAtom%atom(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do
          end do
       end if

       if (NP_Cons > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Constraints relations: ",np_cons
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") "       N.Constr     Name      Father    Factor"
          write(unit=lun, fmt="(a)") "    ============================================="

          np=0
          do i=1,NP_Refi
             n=0
             do j=1,FAtom%natoms
                n=n+count(FAtom%atom(j)%lx ==i)
                n=n+count(FAtom%atom(j)%lu==i)
                if (FAtom%atom(j)%locc==i) n=n+1
                if (FAtom%atom(j)%lbiso==i) n=n+1
             end do
             if ( n > 1) then
                do j=1,FAtom%natoms
                   do k=1,3
                      if (FAtom%atom(j)%lx(k) == i) then
                        car=trim(code_nam(i))//trim(FAtom%atom(j)%lab)
                        if (trim(car)==trim(V_Name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  np, trim(car), &
                              trim(V_Name(i)),FAtom%atom(j)%mx(k)
                      end if
                   end do

                   if (FAtom%atom(j)%lbiso == i) then
                      car=trim(code_nam(4))//trim(FAtom%atom(j)%lab)
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np, trim(car), &
                           trim(V_Name(i)),FAtom%atom(j)%mbiso
                   end if

                   if (FAtom%atom(j)%locc == i) then
                      car=trim(code_nam(5))//trim(FAtom%atom(j)%lab)
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np, trim(car), &
                           trim(V_Name(i)),FAtom%atom(j)%mocc
                   end if

                   do k=1,6
                      if (FAtom%atom(j)%lu(k) == i) then
                         car=trim(code_nam(5+k))//trim(FAtom%atom(j)%lab)
                         if (trim(car)==trim(V_Name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  np, trim(car), &
                              trim(V_Name(i)),FAtom%atom(j)%mu(k)
                      end if
                   end do
                end do
             end if

          end do
       end if

       if (NP_Rest_Dis > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Distance Restraints relations: ",np_rest_dis
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") "   N.Rest  Distance    Sigma    Atom1     Atom2: Symmetry_Op"
          write(unit=lun, fmt="(a)") " ==============================================================="
          fmt3="(i7,tr3,f8.4,tr4,f6.3,t34,a,t43,a)"
          do i=1,np_rest_dis
             p1=dis_rest(i)%p(1)
             p2=dis_rest(i)%p(2)
             n=0
             tr=0.0
             symcar(1)=" "
             if (len_trim(dis_rest(i)%stcode) > 0) then
                call Read_SymTrans_Code(dis_rest(i)%stcode,n,tr)
                tr=tr+Spg%Symop(n)%Tr
                call get_symSymb(Spg%Symop(n)%Rot,Tr,symcar(1))
             end if
             symcar(1)=": "//symcar(1)
             write(unit=lun, fmt=fmt3) i,dis_rest(i)%dobs,dis_rest(i)%sigma,trim(FAtom%Atom(p1)%Lab), &
                  trim(FAtom%Atom(p2)%Lab)//trim(symcar(1))
          end do
       end if

       if (NP_Rest_Ang > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Angle Restraints relations: ",np_rest_ang
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") &
          "   N.Rest   Atom1        Atom2: Symmetry_Op      Atom3: Symmetry_Op              Angle   Sigma"
          write(unit=lun, fmt="(a)") &
          " =============================================================================================="
          fmt4="(i7,tr5,a,t22,a,t50,a,t75,f12.4,tr3,f7.4)"
          do i=1,np_rest_ang
             p1=ang_rest(i)%p(1)
             p2=ang_rest(i)%p(2)
             p3=ang_rest(i)%p(3)
             do j=1,2
                n=0
                tr=0.0
                symcar(j)=" "
                if (len_trim(Ang_rest(i)%stcode(j)) > 0) then
                   call Read_SymTrans_Code(Ang_rest(i)%stcode(j),n,tr)
                   tr=tr+Spg%Symop(n)%Tr
                   call get_symSymb(Spg%Symop(n)%Rot,Tr,symcar(j))
                end if
                symcar(j)=": "//symcar(j)
             end do
             write(unit=lun, fmt=fmt4) i,FAtom%Atom(p1)%Lab,trim(FAtom%Atom(p2)%Lab)//trim(symcar(1)), &
                  trim(FAtom%Atom(p3)%Lab)//trim(symcar(2)),ang_rest(i)%aobs,ang_rest(i)%sigma
          end do
       end if

       if (NP_Rest_Tor > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Torsion Angle Restraints relations: ",np_rest_tor
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") &
          "   N.Rest   Atom1        Atom2: Symmetry_Op      Atom3: Symmetry_Op      Atom4: Symmetry_Op              Angle   Sigma"
          write(unit=lun, fmt="(a)") &
          " ======================================================================================================================"
          fmt5="(i7,tr5,a,t26,a,t34,a,t52,a,t75,f12.4,tr3,f7.4)"
          do i=1,np_rest_tor
             p1=tor_rest(i)%p(1)
             p2=tor_rest(i)%p(2)
             p3=tor_rest(i)%p(3)
             p4=tor_rest(i)%p(4)
             do j=1,3
                n=0
                tr=0.0
                symcar(j)=" "
                if (len_trim(Tor_rest(i)%stcode(j)) > 0) then
                   call Read_SymTrans_Code(Tor_rest(i)%stcode(j),n,tr)
                   tr=tr+Spg%Symop(n)%Tr
                   call get_symSymb(Spg%Symop(n)%Rot,Tr,symcar(j))
                end if
                symcar(j)=": "//symcar(j)
             end do
             write(unit=lun, fmt=fmt5) i,trim(FAtom%Atom(p1)%Lab),trim(FAtom%Atom(p2)%Lab)//trim(symcar(1)), &
                  trim(FAtom%Atom(p3)%Lab)//trim(symcar(2)),trim(FAtom%Atom(p4)%Lab)//trim(symcar(3)),       &
                  tor_rest(i)%tobs,tor_rest(i)%sigma
          end do
       end if

       return
    End Subroutine Write_Info_RefCodes_FAtom

    !!--++
    !!--++ Subroutine Write_Info_RefCodes_Molcrys(MolCrys, iunit)
    !!--++    type(molecular_crystal_type), intent(in) :: molcrys
    !!--++    integer, optional,            intent(in) :: iunit
    !!--++
    !!--++ Overloaded
    !!--++ Write the Information about Refinement Codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Write_Info_RefCodes_Molcrys(MolCrys,iunit)
       !---- Arguments ----!
       type(molecular_crystal_type), intent(in) :: molcrys
       integer, optional,            intent(in) :: iunit

       !---- Local variables ----!
       character(len=20) :: car
       character(len=60) :: fmt1,fmt2
       integer           :: i,j,k,kk,n,na,np, lun
       real              :: mu

       fmt1="(t5,a,t16,i3,t27,a,t33,4(tr6,f8.4),tr8,i2,tr6,f8.3,i9)"
       fmt2="(t10,i3,t21,a,t31,a,t39,f8.4)"

       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " =>  Refinable Parameters: ",np_refi
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,a)")"    Name      N.Param   Code-Name       Value        L.Bound       ",&
                                      "U.Bound         Step    Condition   Multiplier    N.Atom"
          write(unit=lun, fmt="(a,a)")" ==================================================================",&
                                      "========================================================="

          do i=1,Molcrys%N_Free
             do j=1,3
                if (molcrys%atm(i)%lx(j) /=0) then
                   na=molcrys%atm(i)%lx(j)
                   mu=molcrys%atm(i)%mx(j)
                   car=trim(code_nam(j))//trim(molcrys%atm(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             if (molcrys%atm(i)%locc /=0) then
                na=molcrys%atm(i)%locc
                mu=molcrys%atm(i)%mocc
                car=trim(code_nam(5))//trim(molcrys%atm(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             if (molcrys%atm(i)%lbiso /=0) then
                na=molcrys%atm(i)%lbiso
                mu=molcrys%atm(i)%mbiso
                car=trim(code_nam(4))//trim(molcrys%atm(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             do j=1,6
                if (molcrys%atm(i)%lu(j) /=0) then
                   na=molcrys%atm(i)%lu(j)
                   mu=molcrys%atm(i)%mu(j)
                   car=trim(code_nam(5+j))//trim(molcrys%atm(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do
          end do

          do k=1,Molcrys%N_Mol

             do i=1,Molcrys%Mol(k)%natoms
                do j=1,3
                   if (Molcrys%Mol(k)%lI_coor(j,i) /=0) then
                      na=Molcrys%Mol(k)%lI_coor(j,i)
                      mu=Molcrys%Mol(k)%mI_coor(j,i)
                      car=trim(code_nam(j))//trim(molcrys%mol(k)%AtName(i))
                      write(unit=lun,fmt=fmt1)  &
                           trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                   end if
                end do

                if (Molcrys%Mol(k)%locc(i) /=0) then
                   na=Molcrys%Mol(k)%locc(i)
                   mu=Molcrys%Mol(k)%mocc(i)
                   car=trim(code_nam(5))//trim(molcrys%mol(k)%AtName(i))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if

                if (Molcrys%Mol(k)%lbiso(i) /=0) then
                   na=Molcrys%Mol(k)%lbiso(i)
                   mu=Molcrys%Mol(k)%mbiso(i)
                   car=trim(code_nam(4))//trim(molcrys%mol(k)%AtName(i))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             write(unit=lun, fmt="(a)") " "

             do j=1,3
                if (Molcrys%Mol(k)%lxcentre(j) /=0) then
                   na=Molcrys%Mol(k)%lxcentre(j)
                   mu=Molcrys%Mol(k)%mxcentre(j)
                   car=trim(code_nam(12+j))//"entre"
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             do j=1,3
                if (Molcrys%Mol(k)%lOrient(j) /=0) then
                   na=Molcrys%Mol(k)%lOrient(j)
                   mu=Molcrys%Mol(k)%mOrient(j)
                   car=trim(code_nam(15+j))//"Orient"
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lT_TLS(j) /=0) then
                   na=Molcrys%Mol(k)%lT_TLS(j)
                   mu=Molcrys%Mol(k)%mT_TLS(j)
                   car=trim(code_nam(19))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lL_TLS(j) /=0) then
                   na=Molcrys%Mol(k)%lL_TLS(j)
                   mu=Molcrys%Mol(k)%mL_TLS(j)
                   car=trim(code_nam(10))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             do i=1,3
                do j=1,3
                   if (Molcrys%Mol(k)%lS_TLS(i,j) /=0) then
                      na=Molcrys%Mol(k)%lS_TLS(i,j)
                      mu=Molcrys%Mol(k)%mS_TLS(i,j)
                      car=trim(code_nam(21))
                      write(unit=lun,fmt=fmt1)  &
                           trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                   end if
                end do
             end do
          end do
       end if

       if (NP_Cons > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Constraints relations: ",np_cons
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") "       N.Constr     Name      Father    Factor"
          write(unit=lun, fmt="(a)") "    ============================================="

          np=0
          do i=1,NP_Refi
             n=0
             do j=1,Molcrys%N_Free
                n=n+count(molcrys%atm(j)%lx ==i)
                n=n+count(molcrys%atm(j)%lu==i)
                if (molcrys%atm(j)%locc==i) n=n+1
                if (molcrys%atm(j)%lbiso==i) n=n+1
             end do
             if ( n > 1) then
                do j=1,Molcrys%N_Free
                   do k=1,3
                      if (molcrys%atm(j)%lx(k) == i) then
                         car=trim(code_nam(k))//trim(molcrys%atm(j)%lab)
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car),trim(V_Name(i)),molcrys%atm(j)%mx(k)
                      end if
                   end do

                   if (molcrys%atm(j)%lbiso == i) then
                      car=trim(code_nam(4))//trim(molcrys%atm(j)%lab)
                      if (trim(car)==trim(V_name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  &
                           np, trim(car), trim(V_Name(i)),molcrys%atm(j)%mbiso
                   end if

                   if (molcrys%atm(j)%locc == i) then
                      car=trim(code_nam(5))//trim(molcrys%atm(j)%lab)
                      if (trim(car)==trim(V_name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  &
                           np, trim(car), trim(V_Name(i)),molcrys%atm(j)%mocc
                   end if

                   do k=1,6
                      if (molcrys%atm(j)%lu(k) == i) then
                         car=trim(code_nam(5+k))//trim(molcrys%atm(j)%lab)
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car), trim(V_Name(i)),molcrys%atm(j)%mu(k)
                      end if
                   end do
                end do
             end if
          end do

          do i=1,NP_Refi
             do k=1,Molcrys%N_Mol
                n=count(molcrys%mol(k)%lxcentre  == i)
                n=n+count(molcrys%mol(k)%lorient == i)
                n=n+count(molcrys%mol(k)%lT_TLS  == i)
                n=n+count(molcrys%mol(k)%lL_TLS  == i)
                n=n+count(molcrys%mol(k)%lS_TLS  == i)
                if (n > 1) then
                   do j=1,3
                      if (molcrys%mol(k)%lxcentre(j) == i) then
                         car=trim(code_nam(12+j))//"entre"
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car), trim(V_Name(i)),molcrys%mol(k)%mxcentre(j)
                      end if
                   end do

                   do j=1,3
                      if (molcrys%mol(k)%lOrient(j) == i) then
                         car=trim(code_nam(15+j))//"Orient"
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car), trim(V_Name(i)),molcrys%mol(k)%morient(j)
                      end if
                   end do

                   do j=1,6
                      if (molcrys%mol(k)%lT_TLS(j) == i) then
                         car=trim(code_nam(19))
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np,trim(car), trim(V_Name(i)),molcrys%mol(k)%mT_TLS(j)
                      end if
                   end do

                   do j=1,6
                      if (molcrys%mol(k)%lL_TLS(j) == i) then
                         car=trim(code_nam(20))
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car), trim(V_Name(i)),molcrys%mol(k)%mL_TLS(j)
                      end if
                   end do

                   do j=1,3
                      do kk=1,3
                         if (molcrys%mol(k)%lS_TLS(j,kk) == i) then
                            car=trim(code_nam(21))
                            if (trim(car)==trim(V_name(i))) cycle
                            np=np+1
                            write(unit=lun,fmt=fmt2)  &
                                 np,trim(car),trim(V_Name(i)),molcrys%mol(k)%mS_TLS(j,kk)
                         end if
                      end do
                   end do

                end if
             end do
          end do
       end if

       return
    End Subroutine Write_Info_RefCodes_Molcrys

    !!--++
    !!--++ Subroutine Write_Info_RefCodes_Molec(Molec, iunit)
    !!--++    type(molecule_type), intent(in) :: molec
    !!--++    integer, optional,   intent(in) :: iunit
    !!--++
    !!--++ Overloaded
    !!--++ Write the Information about Refinement Codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Write_Info_RefCodes_Molec(Molec,iunit)
       !---- Arguments ----!
       type(molecule_type), intent(in) :: molec
       integer, optional,   intent(in) :: iunit

       !---- Local variables ----!
       character(len=60) :: fmt1,fmt2
       character(len=20) :: car
       integer           :: i,j,k,n,na,np,lun
       real              :: mu

       fmt1="(t5,a,t16,i3,t27,a,t33,4(tr6,f8.4),tr8,i2,tr6,f8.3,i9)"
       fmt2="(t10,i3,t21,a,t31,a,t39,f8.4)"

       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Refinable Parameters: ",np_refi
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,a)")"    Name      N.Param   Code-Name       Value        L.Bound       ",&
                                      "U.Bound         Step    Condition   Multiplier    N.Atom"
          write(unit=lun, fmt="(a,a)")" ==================================================================",&
                                      "========================================================="

          do i=1,Molec%natoms
             do j=1,3
                if (molec%lI_coor(j,i) /=0) then
                   na=molec%lI_coor(j,i)
                   mu=molec%mI_coor(j,i)
                   car=trim(code_nam(j))//trim(molec%AtName(i))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             if (molec%locc(i) /=0) then
                na=molec%locc(i)
                mu=molec%mocc(i)
                car=trim(code_nam(4))//trim(molec%AtName(i))
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             if (molec%lbiso(i) /=0) then
                na=molec%lbiso(i)
                mu=molec%mbiso(i)
                car=trim(code_nam(5))//trim(molec%AtName(i))
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          write(unit=lun, fmt="(a)") " "

          do j=1,3
             if (molec%lxcentre(j) /=0) then
                na=molec%lxcentre(j)
                mu=molec%mxcentre(j)
                car=trim(code_nam(12+j))//"entre"
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          do j=1,3
             if (molec%lOrient(j) /=0) then
                na=molec%lOrient(j)
                mu=molec%mOrient(j)
                car=trim(code_nam(15+j))//"Orient"
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          do j=1,6
             if (molec%lT_TLS(j) /=0) then
                na=molec%lT_TLS(j)
                mu=molec%mT_TLS(i)
                car=trim(code_nam(19))
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          do j=1,6
             if (molec%lL_TLS(j) /=0) then
                na=molec%lL_TLS(j)
                mu=molec%mL_TLS(i)
                car=trim(code_nam(20))
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          do i=1,3
             do j=1,3
                if (molec%lS_TLS(i,j) /=0) then
                   na=molec%lS_TLS(i,j)
                   mu=molec%mS_TLS(i,j)
                   car=trim(code_nam(21))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do
          end do
       end if

       if (NP_Cons > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Constraints relations: ",np_cons
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") "       N.Constr     Name      Father    Factor"
          write(unit=lun, fmt="(a)") "    ============================================="

          np=0
          do i=1,NP_Refi
             n=0
             do j=1,Molec%natoms
                n=n+count(molec%lI_coor(:,j) ==i)
                if (molec%locc(j)==i) n=n+1
                if (molec%lbiso(j)==i) n=n+1
             end do
             if ( n > 1) then
                do j=1,Molec%natoms
                   do k=1,3
                      if (molec%lI_coor(k,j) == i) then
                         car=trim(code_nam(k))//trim(molec%AtName(j))
                         if (trim(car)==trim(V_Name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  np, trim(car),trim(V_Name(i)),molec%mI_coor(k,j)
                      end if
                   end do

                   if (molec%lbiso(j) == i) then
                      car=trim(code_nam(5))//trim(molec%AtName(j))
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2) np, trim(car), trim(V_Name(i)),molec%mbiso(j)
                   end if

                   if (molec%locc(j) == i) then
                      car=trim(code_nam(4))//trim(molec%AtName(j))
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np, trim(car),trim(V_Name(i)),molec%mocc(j)
                   end if
                end do
             end if
          end do

          do i=1,NP_Refi
             n=count(molec%lxcentre  == i)
             n=n+count(molec%lorient == i)
             n=n+count(molec%lT_TLS  == i)
             n=n+count(molec%lL_TLS  == i)
             n=n+count(molec%lS_TLS  == i)
             if (n > 1) then
                do j=1,3
                   if (molec%lxcentre(j) == i) then
                      car=trim(code_nam(12+j))//"entre"
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np, trim(car),trim(V_Name(i)),molec%mxcentre(j)
                   end if
                end do

                do j=1,3
                   if (molec%lOrient(j) == i) then
                      car=trim(code_nam(15+j))//"Orient"
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np,trim(car),trim(V_Name(i)),molec%morient(j)
                   end if
                end do

                do j=1,6
                   if (molec%lT_TLS(j) == i) then
                      car=trim(code_nam(19))
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2) np,trim(car), trim(V_Name(i)),molec%mT_TLS(j)
                  end if
                end do

                do j=1,6
                   if (molec%lL_TLS(j) == i) then
                      car=trim(code_nam(20))
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2) np,trim(car), trim(V_Name(i)),molec%mL_TLS(j)
                   end if
                end do

                do j=1,3
                   do k=1,3
                      if (molec%lS_TLS(j,k) == i) then
                         car=trim(code_nam(21))
                         if (trim(car)==trim(V_Name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2) np,trim(car), trim(V_Name(i)),molec%mS_TLS(j,k)
                      end if
                   end do
                end do
             end if
          end do

       end if

       return
    End Subroutine Write_Info_RefCodes_Molec

    !!--++
    !!--++ Subroutine Write_Info_RefCodes_MagStr(FmAtom, Mag_dom, MGp, Iunit)
    !!--++    type(mAtom_List_Type),               intent(in) :: FmAtom
    !!--++    type(Magnetic_Domain_type),optional, intent(in) :: Mag_dom
    !!--++    type(MagSymm_k_Type),   intent(in)    :: MGp
    !!--++    integer, optional,      intent(in)    :: Iunit
    !!--++
    !!--++ Write the Information about the Magnetic Refinement Codes
    !!--++    magnetic domains
    !!--++ Created: February - 2012
    !!
    Subroutine Write_Info_RefCodes_MagStr(FmAtom, Mag_dom, MGp, Iunit)
       !---- Arguments ----!
       type(mAtom_List_Type),               intent(in) :: FmAtom
       type(Magnetic_Domain_type),optional, intent(in) :: Mag_dom
       type(MagSymm_k_Type),      intent(in) :: MGp
       integer, optional,         intent(in) :: Iunit

       !---- Local variables ----!
       character(len=20)              :: car
       character(len=60)              :: fmt1,fmt2
       integer                        :: i,j,na,lun,ik,ich
       real(kind=cp)                  :: mu

       !---- Format Zone ----!
       fmt1="(t5,a,t16,i3,t27,a,t33,4(tr6,f8.4),tr8,i2,tr6,f8.3,i9)"
       fmt2="(t10,i3,t21,a,t31,a,t39,f8.4)"
       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Refinable Parameters: ",np_refi
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,a)")"    Name      N.Param   Code-Name       Value        L.Bound       ",&
                                      "U.Bound         Step    Condition   Multiplier    N.Atom"
          write(unit=lun, fmt="(a,a)")" ==================================================================",&
                                      "========================================================="

          do i=1,FmAtom%natoms
             !---- Get im, number of the magnetic matrices/irrep set
             ik=FmAtom%atom(i)%nvk

             !----Real components
             do j=1,3
                if (FmAtom%atom(i)%lSkR(j,ik) /= 0) then
                   na=FmAtom%atom(i)%lSkR(j,ik)
                   mu=FmAtom%atom(i)%mSkR(j,ik)
                 if(MGp%Sk_type == "Spherical_Frame") then
                   car=trim(mcode_nam(j+6))//trim(FmAtom%atom(i)%lab)
                 else
                   car=trim(mcode_nam(j  ))//trim(FmAtom%atom(i)%lab)
                 end if
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)

                end if
             end do

             !----Imaginary components
             do j=1,3
                if (FmAtom%atom(i)%lSkI(j,ik) /= 0) then
                   na=FmAtom%atom(i)%lSkI(j,ik)
                   mu=FmAtom%atom(i)%mSkI(j,ik)
                 if(MGp%Sk_type == "Spherical_Frame") then
                   car=trim(mcode_nam(j+9))//trim(FmAtom%atom(i)%lab)
                 else
                   car=trim(mcode_nam(j+3))//trim(FmAtom%atom(i)%lab)
                 end if
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             if (FmAtom%atom(i)%lmphas(ik) /=0) then
                na=FmAtom%atom(i)%lmphas(ik)
                mu=FmAtom%atom(i)%mmphas(ik)
                car=trim(mcode_nam(13))//trim(FmAtom%atom(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          !----C1-12 coefficients
          do i=1,FmAtom%natoms
             do j=1,12
                if (FmAtom%atom(i)%lbas(j,ik) /= 0) then
                   na=FmAtom%atom(i)%lbas(j,ik)
                   mu=FmAtom%atom(i)%mbas(j,ik)
                   car=trim(mcode_nam(j+13))//trim(FmAtom%atom(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do
          end do

       !---- Check is chirality is present ----!
        if (Mag_Dom%chir) then
         ich=2
        else
         ich=1
        end if

          do i=1,Mag_Dom%nd
             do j=1,ich
                if (Mag_Dom%Lpop(j,i) /= 0) then
                   na=Mag_Dom%Lpop(j,i)
                   mu=Mag_Dom%Mpop(j,i)
                   car=trim(Mag_Dom%Lab(j,i))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do
          end do

       end if
       return

    End Subroutine Write_Info_RefCodes_MagStr

    !!----
    !!---- Subroutine Write_Info_RefGCodes(model, Iunit)
    !!----    type(Nonatomic_Parameter_List_Type), intent(in) :: model
    !!----    integer, optional,                   intent(in) :: Iunit
    !!----
    !!----    Write the Information about Refinement Codes of non-atomic parameters
    !!----
    !!---- Update: November 2 - 2013
    !!
    Subroutine Write_Info_RefGCodes(model, Iunit)
       !---- Arguments ----!
       type(Nonatomic_Parameter_List_Type), intent(in) :: model
       integer, optional,                   intent(in) :: Iunit

       !---- Local variables ----!
       character(len=20)              :: car
       character(len=60)              :: fmt1,fmt2 !,fmt3,fmt4,fmt5
       !Character(len=25),dimension(3) :: symcar
       integer                        :: i,na,lun !,n,j,k,np,p1,p2,p3,p4
       real(kind=cp)                  :: mu
       !real(kind=cp),dimension(3)     :: tr

       !---- Format Zone ----!
       fmt1="(t5,a,t16,i3,t27,a,t33,4(tr6,f8.4),tr8,i2,tr6,f8.3,i9)"
       fmt2="(t10,i3,t21,a,t31,a,t39,f8.4)"

       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Non-atomic refinable Parameters: ",model%npar
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,a)")"    Name      N.Param   Code-Name       Value        L.Bound       ",&
                                      "U.Bound         Step    Condition   Multiplier   Order in Model"
          write(unit=lun, fmt="(a,a)")" ==================================================================",&
                                      "==========================================================="
          do i=1,model%npar
             if (model%par(i)%lcode > 0) then
                na=model%par(i)%lcode
                mu=model%par(i)%multip
                car=trim(model%par(i)%nam)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                !write(unit=*,fmt="(2i6,a)") i,na, "    "//trim(car)
             !else
                !write(unit=*,fmt="(2i6,a)") i,model%par(i)%lcode, "    "//trim(model%par(i)%nam)
             end if
          end do
       end if


       return
    End Subroutine Write_Info_RefGCodes



    !!----
    !!---- Subroutine Write_Info_RefParams(iunit)
    !!----    integer, optional,  intent(in) :: iunit
    !!----
    !!----    Write the Information about Refinement parameters in file associated with
    !!----    logical unit "iunit". If no argument is passed the standard output (iunit=6)
    !!----    is used.
    !!----
    !!---- Update: August - 2007
    !!
    Subroutine Write_Info_RefParams(iunit)
       !---- Arguments ----!
       integer, optional,   intent(in) :: iunit

       !---- Local variables ----!
       integer           :: i,lun

       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Number of Refinable Parameters: ",np_refi
          write(unit=lun, fmt="(a)") " "
          do i=1,nP_refi
            write(unit=lun,fmt="(i6,tr5,a20,4f14.5, i4, i6)") i,V_Name(i),V_Vec(i),V_Bounds(:,i),V_BCon(i),V_List(i)
          end do
       end if

       return
    End Subroutine Write_Info_RefParams

    !!----
    !!---- Subroutine Write_Restraints_ObsCalc(A,iunit)
    !!----    type(Atom_List_Type),intent(in) :: A
    !!----    integer, optional,   intent(in) :: iunit
    !!----
    !!----    Write the current values of the "observed" and calculated
    !!----    restraints, as well as the corresponding cost value.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Write_Restraints_ObsCalc(A,iunit)
       !---- Arguments ----!
       type(Atom_List_Type),intent(in) :: A
       integer, optional,   intent(in) :: iunit

       !---- Local variables ----!
       character(len=14) :: car1,car2,car3,car4
       integer           :: i,i1,i2,i3,i4,lun
       real              :: disto,distc,ango,angc,sigm, cost,w, delta

       lun=6
       if (present(iunit)) lun=iunit

       if( NP_Rest_Dis > 0) then
          Write(unit=lun,fmt="(/,a)") " ============================================================"
          Write(unit=lun,fmt="(a)")   "   Distance Restraints: Atoms, Dobs, Dcalc, Sigma, delt/Sigma"
          Write(unit=lun,fmt="(a,/)") " ============================================================"
          Write(unit=lun,fmt="(a)") " Rest#    Atom1         Atom2              Dobs        Dcalc       Sigma   (Do-Dc)/Sigma"
          cost=0.0
          do i=1,NP_Rest_Dis
             i1=Dis_rest(i)%p(1)
             i2=Dis_rest(i)%p(2)
             car1=trim(A%Atom(i1)%lab)
             car2=trim(A%Atom(i2)%lab)//dis_rest(i)%stcode
             disto=Dis_rest(i)%dobs
             distc=Dis_rest(i)%dcalc
             delta=disto-distc
             sigm=Dis_rest(i)%sigma
             w= 1.0/(sigm*sigm)
             cost= cost+delta*delta*w
             Write(unit=lun,fmt="(i6,tr4,2a,4f12.5)") i,car1,car2,disto,distc,sigm,delta/sigm
          end do

          Write(unit=lun,fmt="(/,a,f12.5)") "   Distance Restraints Cost = Sum{[(dobs-dcalc)/Sigma]^2} = ",cost
       end if


       if( NP_Rest_Ang > 0) then
          Write(unit=lun,fmt="(/,a)") " ============================================================="
          Write(unit=lun,fmt="(a)")   "   Angle Restraints: Atoms, Angobs, Angcalc, Sigma, delt/Sigma"
          Write(unit=lun,fmt="(a,/)") " ============================================================="
          Write(unit=lun,fmt="(a)") &
          " Rest#    Atom1         Atom2          Atom3            Ang_obs    Ang_calc      Sigma   (Ao-Ac)/Sigma"

          cost=0.0
          do i=1,NP_Rest_Ang
             i1=Ang_rest(i)%p(1)
             i2=Ang_rest(i)%p(2)
             i3=Ang_rest(i)%p(3)
             car1=trim(A%Atom(i1)%lab)
             car2=trim(A%Atom(i2)%lab)//ang_rest(i)%stcode(1)
             car3=trim(A%Atom(i3)%lab)//ang_rest(i)%stcode(2)
             ango=Ang_rest(i)%Aobs
             angc=Ang_rest(i)%Acalc
             delta=ango-angc
             sigm=Ang_rest(i)%sigma
             w= 1.0/(sigm*sigm)
             cost= cost+delta*delta*w
             Write(unit=lun,fmt="(i6,tr4,3a,4f12.5)") i,car1,car2,car3,ango,angc,sigm,delta/sigm
          end do

          Write(unit=lun,fmt="(/,a,f12.5)") "   Angle Restraints Cost = Sum{[(Ang_obs-Ang_calc)/Sigma]^2} = ",cost
       End If



       if( NP_Rest_tor > 0) then
          Write(unit=lun,fmt="(/,a)") " ====================================================================="
          Write(unit=lun,fmt="(a)")   "   Torsion Angle Restraints: Atoms, Angobs, Angcalc, Sigma, delt/Sigma"
          Write(unit=lun,fmt="(a,/)") " ====================================================================="
          Write(unit=lun,fmt="(a)") " Rest#    Atom1         Atom2          Atom3          Atom4            "//&
               "Ang_obs    Ang_calc      Sigma   (Ao-Ac)/Sigma"

          cost=0.0
          do i=1,NP_Rest_tor
             i1=Tor_rest(i)%p(1)
             i2=Tor_rest(i)%p(2)
             i3=Tor_rest(i)%p(3)
             i4=Tor_rest(i)%p(4)
             car1=trim(A%Atom(i1)%lab)
             car2=trim(A%Atom(i2)%lab)//tor_rest(i)%stcode(1)
             car3=trim(A%Atom(i3)%lab)//tor_rest(i)%stcode(2)
             car4=trim(A%Atom(i4)%lab)//tor_rest(i)%stcode(3)
             ango=Ang_rest(i)%Aobs
             angc=Ang_rest(i)%Acalc
             delta=ango-angc
             sigm=Ang_rest(i)%sigma
             w= 1.0/(sigm*sigm)
             cost= cost+delta*delta*w
             Write(unit=lun,fmt="(i6,tr4,4a,4f12.5)") i,car1,car2,car3,car4,ango,angc,sigm,delta/sigm
          end do

          Write(unit=lun,fmt="(/,a,f12.5)") "   Torsion Angle Restraints Cost = Sum{[(Ang_obs-Ang_calc)/Sigma]^2} = ",cost
       End If

       return
    End Subroutine Write_Restraints_ObsCalc

 End Module CFML_Keywords_Code_Parser
