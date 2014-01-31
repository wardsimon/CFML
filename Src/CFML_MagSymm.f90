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
!!---- MODULE: CFML_Magnetic_Symmetry
!!----   INFO: Series of procedures handling operations with Magnetic Symmetry
!!----         and Magnetic Structures
!!----
!!---- HISTORY
!!----
!!----    Update: 07/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,                only: cp, tpi
!!--++    Use CFML_Math_General,              only: Modulo_Lat
!!--++    Use CFML_Math_3D,                   only: Get_Cart_From_Spher
!!--++    Use CFML_Symmetry_Tables,           only: ltr_a,ltr_b,ltr_c,ltr_i,ltr_r,ltr_f
!!--++    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Read_Xsym, Get_SymSymb, &
!!--++                                              Sym_Oper_Type, Set_SpaceGroup,read_msymm, symmetry_symbol, &
!!--++                                              err_symm,err_symm_mess, set_SpG_Mult_Table,ApplySO,   &
!!--++                                              Lattice_Trans
!!--++    Use CFML_String_Utilities,          only: u_case, l_case, Frac_Trans_1Dig, Get_Separator_Pos,Pack_String, &
!!--++                                              Frac_Trans_2Dig, Get_Mat_From_Symb, getnum_std, Err_String,     &
!!--++                                              Err_String_Mess,setnum_std, getword
!!--++    Use CFML_IO_Formats,                only: file_list_type, File_To_FileList
!!--++    Use CFML_Atom_TypeDef,              only: Allocate_mAtom_list, mAtom_List_Type
!!--++    Use CFML_Scattering_Chemical_Tables,only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
!!--++                                              Magnetic_Form
!!--++    Use CFML_Propagation_Vectors,       only: K_Equiv_Minus_K
!!--++    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Set_Crystal_Cell
!!----
!!---- VARIABLES
!!--..    Types
!!----    MSYM_OPER_TYPE
!!----    MAGNETIC_DOMAIN_TYPE
!!----    MAGNETIC_GROUP_TYPE
!!----    MAGSYMM_K_TYPE
!!--..
!!----    ERR_MAGSYM
!!----    ERR_MAGSYM_MESS
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       APPLYMSO
!!----
!!----    Subroutines:
!!----       INIT_ERR_MAGSYM
!!----       INIT_MAGSYMM_K_TYPE             !OZ made it public to use in Read_Refcodes_Magnetic_Structure
!!----       READN_SET_MAGNETIC_STRUCTURE
!!--++       READN_SET_MAGNETIC_STRUCTURE_CFL    [Overloaded]
!!--++       READN_SET_MAGNETIC_STRUCTURE_MCIF   [Overloaded]
!!----       SET_SHUBNIKOV_GROUP
!!----       WRITE_MAGNETIC_STRUCTURE
!!----       WRITE_SHUBNIKOV_GROUP
!!----
!!
 Module CFML_Magnetic_Symmetry

    !---- Use Modules ----!
    Use CFML_GlobalDeps,                only: cp, tpi,Write_Date_Time
    Use CFML_Math_General,              only: Modulo_Lat
    Use CFML_Math_3D,                   only: Get_Cart_From_Spher
    Use CFML_Symmetry_Tables,           only: ltr_a,ltr_b,ltr_c,ltr_i,ltr_r,ltr_f
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Read_Xsym, Get_SymSymb, &
                                              Sym_Oper_Type, Set_SpaceGroup,read_msymm, symmetry_symbol, &
                                              err_symm,err_symm_mess, set_SpG_Mult_Table,ApplySO,   &
                                              Lattice_Trans
    Use CFML_String_Utilities,          only: u_case, l_case, Frac_Trans_1Dig, Get_Separator_Pos,Pack_String, &
                                              Frac_Trans_2Dig, Get_Mat_From_Symb, getnum_std, Err_String,     &
                                              Err_String_Mess,setnum_std, getword
    Use CFML_IO_Formats,                only: file_list_type, File_To_FileList
    Use CFML_Atom_TypeDef,              only: Allocate_mAtom_list, mAtom_List_Type
    Use CFML_Scattering_Chemical_Tables,only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
                                              Magnetic_Form
    Use CFML_Propagation_Vectors,       only: K_Equiv_Minus_K
    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Set_Crystal_Cell

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public :: ApplyMSO

    !---- List of public subroutines ----!
    public :: Readn_Set_Magnetic_Structure, Write_Magnetic_Structure, Set_Shubnikov_Group, &
              Write_Shubnikov_Group, Init_MagSymm_k_Type, Write_MCIF


    !---- Definitions ----!

    !!----
    !!---- TYPE :: MSYM_OPER_TYPE
    !!--..
    !!---- Type, public :: MSym_Oper_Type
    !!----    integer, dimension(3,3) :: Rot     !  Rotational Part of Symmetry Operator
    !!----    real(kind=cp)           :: Phas    !  Phase in fraction of 2pi
    !!---- End Type  MSym_Oper_Type
    !!----
    !!----  Definition of Magnetic symmetry operator type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: MSym_Oper_Type
       integer, dimension(3,3) :: Rot
       real(kind=cp)           :: Phas
    End Type MSym_Oper_Type

    !!----
    !!---- TYPE :: MAGNETIC_DOMAIN_TYPE
    !!--..
    !!---- Type, public :: Magnetic_Domain_type
    !!----    integer                           :: nd=0          !Number of rotational domains (not counting chiral domains)
    !!----    logical                           :: Chir=.false.  !True if chirality domains exist
    !!----    logical                           :: trans=.false. !True if translations are associated to matrix domains
    !!----    logical                           :: Twin=.false.  !True if domains are to be interpreted as twins
    !!----    integer,dimension(3,3,24)         :: DMat=0        !Domain matrices to be applied to Fourier Coefficients
    !!----    real(kind=cp), dimension (2,24)   :: Dt=0.0        !Translations associated to rotation matrices
    !!----    real(kind=cp), dimension (2,24)   :: pop=0.0       !Populations of domains (sum=1,
    !!----                                                       !the second value is /=0 for chir=.true.)
    !!----    real(kind=cp), dimension (2,24)   :: pop_std=0.0   !Standard deviations of Populations of domains
    !!----    integer,dimension (2,24)          :: Lpop=0        !Number of the refined parameter
    !!----    real(kind=cp), dimension (2,24)   :: Mpop=0.0      !Refinement codes for populations
    !!----    character(len=10),dimension (2,24):: Lab           !Label of domain
    !!---- End type Magnetic_Domain_type
    !!----
    !!----
    !!--<<
    !!----  Magnetic S-domains corresponds to a different magnetic structure obtained from
    !!----  the domain 1 (actual model) by applying a rotational operator to the Fourier
    !!----  coefficients of magnetic moments. This rotational operator corresponds to a
    !!----  symmetry operator of the paramagnetic group that is lost in the ordered state.
    !!----  Chirality domains are simply obtained by changing the sign of the imaginary
    !!----  components of the Fourier coefficients. For each rotational domain two chiralities
    !!----  domains exist.
    !!-->>
    !!----
    !!---- Updated: October - 2006, July-2012 (JRC, more type of domains), November 2013 (standard deviations)
    !!
    Type, public :: Magnetic_Domain_type
       integer                           :: nd=0          !Number of rotational domains (not counting chiral domains)
       logical                           :: Chir=.false.  !True if chirality domains exist
       logical                           :: trans=.false. !True if translations are associated to matrix domains
       logical                           :: Twin=.false.  !True if domains are to be interpreted as twins
       integer,dimension(3,3,24)         :: DMat=0        !Domain matrices to be applied to Fourier Coefficients
       real(kind=cp), dimension (3,24)   :: Dt=0.0        !Translations associated to rotation matrices
       real(kind=cp), dimension (2,24)   :: pop=0.0       !Populations of domains (sum=1,
                                                          !the second value is /=0 for chir=.true.)
       integer      , dimension (2,24)   :: Lpop=0        !Number of the refined parameter
       real(kind=cp), dimension (2,24)   :: Mpop=0.0      !Multipliers for population
       real(kind=cp), dimension (2,24)   :: pop_std=0.0   !Standard deviations of Populations of domains
       character(len=10),dimension (2,24):: Lab           !Label of domain
    End type Magnetic_Domain_type

    !!----
    !!---- TYPE :: MAGNETIC_SPACE_GROUP_TYPE
    !!--..
    !!---- Type, Public :: Magnetic_Space_Group_Type
    !!----    Integer                              :: number
    !!----    character(len=15)                    :: BNS_number
    !!----    character(len=15)                    :: OG_number
    !!----    Character(len=34)                    :: BNS_symbol
    !!----    Character(len=34)                    :: OG_symbol
    !!----    Integer                              :: MagType
    !!----    Integer                              :: Parent_num
    !!----    Character(len=20)                    :: Parent_spg
    !!----    logical                              :: standard_setting  !true or false ! 'yes' or 'no'
    !!----    logical                              :: mcif !true if mx,my,mz notation is used , false is u,v,w notation is used
    !!----    logical                              :: m_cell !true if magnetic cell is used for symmetry operators
    !!----    logical                              :: m_constr !true if constraints have been provided, strings are in atom types
    !!----    Character(len=40)                    :: trn_from_parent
    !!----    Character(len=40)                    :: trn_to_standard
    !!----    Integer                              :: n_sym
    !!----    Integer                              :: n_wyck  !Number of Wyckoff positions of the magnetic group
    !!----    Integer                              :: n_kv
    !!----    Integer                              :: n_irreps
    !!----    Integer,             dimension(:),allocatable  :: irrep_dim       !Dimension of the irreps
    !!----    Integer,             dimension(:),allocatable  :: small_irrep_dim !Dimension of the small irrep
    !!----    Character(len=15),   dimension(:),allocatable  :: irrep_id        !Labels for the irreps
    !!----    Character(len=20),   dimension(:),allocatable  :: irrep_direction !Irrep direction in representation space
    !!----    Character(len=20),   dimension(:),allocatable  :: irrep_action    !Irrep character primary or secondary
    !!----    Character(len=15),   dimension(:),allocatable  :: kv_label
    !!----    real(kind=cp),     dimension(:,:),allocatable  :: kv
    !!----    character(len=40),   dimension(:),allocatable  :: Wyck_Symb  ! Alphanumeric Symbols for first representant of Wyckoff positions
    !!----    character(len=40),   dimension(:),allocatable  :: SymopSymb  ! Alphanumeric Symbols for SYMM
    !!----    type(Sym_Oper_Type), dimension(:),allocatable  :: SymOp      ! Crystallographic symmetry operators
    !!----    character(len=40),   dimension(:),allocatable  :: MSymopSymb ! Alphanumeric Symbols for MSYMM
    !!----    type(MSym_Oper_Type),dimension(:),allocatable  :: MSymOp     ! Magnetic symmetry operators
    !!---- End Type Magnetic_Space_Group_Type
    !!----
    !!--<<
    !!----    The magnetic group type defined here satisfy all the needs for working with
    !!----    standard data bases for BNS and OG notations and also for working with
    !!----    simplified methods with the crystallographic cell and propagation vectors
    !!----    The component Phas in MSym_Oper_Type is used for time inversion Phas=+1 no time inversion
    !!----    and Phas=-1 if time inversion is associated with the operator (Not needed for real calculations).
    !!-->>
    !!----
    !!----  Created: January - 2014
    !!
    Type, Public :: Magnetic_Space_Group_Type
       Integer                              :: Sh_number
       character(len=15)                    :: BNS_number
       character(len=15)                    :: OG_number
       Character(len=34)                    :: BNS_symbol
       Character(len=34)                    :: OG_symbol
       Integer                              :: MagType
       Integer                              :: Parent_num
       Character(len=20)                    :: Parent_spg
       logical                              :: standard_setting  !true or false
       logical                              :: mcif !true if mx,my,mz notation is used , false is u,v,w notation is used
       logical                              :: m_cell !true if magnetic cell is used for symmetry operators
       logical                              :: m_constr !true if constraints have been provided
       Character(len=40)                    :: trn_from_parent
       Character(len=40)                    :: trn_to_standard
       Integer                              :: n_sym
       Integer                              :: n_wyck   !Number of Wyckoff positions of the magnetic group
       Integer                              :: n_kv
       Integer                              :: n_irreps
       Integer,             dimension(:),allocatable  :: irrep_dim       !Dimension of the irreps
       Integer,             dimension(:),allocatable  :: small_irrep_dim !Dimension of the small irrep
       Integer,             dimension(:),allocatable  :: irrep_modes_number !Number of the mode of the irrep
       Character(len=15),   dimension(:),allocatable  :: irrep_id        !Labels for the irreps
       Character(len=20),   dimension(:),allocatable  :: irrep_direction !Irrep direction in representation space
       Character(len=20),   dimension(:),allocatable  :: irrep_action    !Irrep character primary or secondary
       Character(len=15),   dimension(:),allocatable  :: kv_label
       real(kind=cp),     dimension(:,:),allocatable  :: kv
       character(len=40),   dimension(:),allocatable  :: Wyck_Symb  ! Alphanumeric Symbols for first representant of Wyckoff positions
       character(len=40),   dimension(:),allocatable  :: SymopSymb  ! Alphanumeric Symbols for SYMM
       type(Sym_Oper_Type), dimension(:),allocatable  :: SymOp      ! Crystallographic symmetry operators
       character(len=40),   dimension(:),allocatable  :: MSymopSymb ! Alphanumeric Symbols for MSYMM
       type(MSym_Oper_Type),dimension(:),allocatable  :: MSymOp     ! Magnetic symmetry operators
    End Type Magnetic_Space_Group_Type

    !!----
    !!---- TYPE :: MAGNETIC_GROUP_TYPE
    !!--..
    !!---- Type, Public :: Magnetic_Group_Type
    !!----    Character(len=30)           :: Shubnikov !Shubnikov symbol (Hermman-Mauguin + primes)
    !!----    type(Space_Group_Type)      :: SpG       !Crystallographic space group
    !!----    integer, dimension(192)     :: tinv      !When a component is +1 no time inversion is associated
    !!---- End Type Magnetic_Group_Type                !If tinv(i)=-1, the time inversion is associated to operator "i"
    !!----
    !!--<<
    !!----    A magnetic group type is adequate when k=(0,0,0). It contains as the second
    !!----    component the crystallographic space group. The first component is
    !!----    the Shubnikov Group symbol and the third component is an integer vector with
    !!----    values -1 or 1 when time inversion is associated (-1) with the corresponding
    !!----    crystallographic symmetry operator o not (1).
    !!-->>
    !!----
    !!---- Update: April - 2005
    !!
    Type, Public :: Magnetic_Group_Type
       Character(len=30)           :: Shubnikov
       type(Space_Group_Type)      :: SpG
       integer, dimension(192)     :: tinv
    End Type Magnetic_Group_Type
    !!----
    !!---- TYPE :: MAGSYMM_K_TYPE
    !!--..
    !!---- Type, Public :: MagSymm_k_Type
    !!----    character(len=31)                        :: MagModel   ! Name to characterize the magnetic symmetry
    !!----    character(len=10)                        :: Sk_type    ! If Sk_type="Spherical_Frame" the input Fourier coefficients are in spherical components
    !!----    character(len=1)                         :: Latt       ! Symbol of the crystallographic lattice
    !!----    integer                                  :: nirreps    ! Number of irreducible representations (max=4, if nirreps /= 0 => nmsym=0)
    !!----    integer                                  :: nmsym      ! Number of magnetic operators per crystallographic operator (max=8)
    !!----    integer                                  :: centred    ! =0 centric centre not at origin, =1 acentric, =2 centric (-1 at origin)
    !!----    integer                                  :: mcentred   ! =1 Anti/a-centric Magnetic symmetry, = 2 centric magnetic symmetry
    !!----    integer                                  :: nkv        ! Number of independent propagation vectors
    !!----    real(kind=cp),dimension(3,12)            :: kvec       ! Propagation vectors
    !!----    integer                                  :: NumLat     ! Number of centring lattice vectors
    !!----    real(kind=cp), dimension(3,4)            :: Ltr        ! Centring translations
    !!----    integer                                  :: Numops     ! Reduced number of crystallographic Symm. Op.
    !!----    integer                                  :: Multip     ! General multiplicity of the space group
    !!----    integer,             dimension(4)        :: nbas       ! Number of basis functions per irrep (if nbas < 0, the corresponding basis is complex).
    !!----    integer,             dimension(12,4)     :: icomp      ! Indicator (0 pure real/ 1 pure imaginary) for coefficients of basis fucntions
    !!----    character(len=40),   dimension(48)       :: SymopSymb  ! Alphanumeric Symbols for SYMM
    !!----    type(Sym_Oper_Type), dimension(48)       :: SymOp      ! Crystallographic symmetry operators
    !!----    character(len=40),   dimension(48,8)     :: MSymopSymb ! Alphanumeric Symbols for MSYMM
    !!----    type(MSym_Oper_Type),dimension(48,8)     :: MSymOp     ! Magnetic symmetry operators
    !!----    Complex(kind=cp),    dimension(3,12,48,4):: basf       ! Basis functions of the irreps of Gk
    !!---- End Type MagSymm_k_Type
    !!----
    !!----  Definition of the MagSymm_k_type derived type, encapsulating the information
    !!----  concerning the crystallographic symmetry, propagation vectors and magnetic matrices.
    !!----  Needed for calculating magnetic structure factors.
    !!----
    !!---- Created: April - 2005
    !!---- Updated: April - 2005
    !!
    Type, Public :: MagSymm_k_Type
       character(len=31)                        :: MagModel
       character(len=15)                        :: Sk_type
       character(len=1)                         :: Latt
       integer                                  :: nirreps
       integer                                  :: nmsym
       integer                                  :: centred
       integer                                  :: mcentred
       integer                                  :: nkv
       real(kind=cp),dimension(3,12)            :: kvec
       integer                                  :: NumLat
       real(kind=cp), dimension(3,4)            :: Ltr
       integer                                  :: Numops
       integer                                  :: Multip
       integer,             dimension(4)        :: nbas
       integer,             dimension(12,4)     :: icomp
       character(len=40),   dimension(48)       :: SymopSymb
       type( Sym_Oper_Type),dimension(48)       :: SymOp
       character(len=40),   dimension(48,8)     :: MSymopSymb
       type(MSym_Oper_Type),dimension(48,8)     :: MSymOp
       Complex(kind=cp),    dimension(3,12,48,4):: basf
    End Type MagSymm_k_Type

    !!----
    !!---- ERR_MAGSYM
    !!----    logical, public :: err_MagSym
    !!----
    !!----    Logical Variable indicating an error in CFML_Magnetic_Symmetry
    !!----
    !!---- Update: April - 2005
    !!
    logical, public :: err_MagSym

    !!----
    !!---- ERR_MAGSYM_MESS
    !!----    character(len=150), public :: ERR_MagSym_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: April - 2005
    !!
    character(len=150), public :: ERR_MagSym_Mess

    Interface  Readn_Set_Magnetic_Structure
       Module Procedure Readn_Set_Magnetic_Structure_CFL
       Module Procedure Readn_Set_Magnetic_Structure_MCIF
    End Interface  Readn_Set_Magnetic_Structure

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!----
    !!---- Function ApplyMso(Op,Sk) Result(Skp)
    !!----    Type(MSym_Oper_Type),   intent(in) :: Op        !  Magnetic Symmetry Operator Type
    !!----    complex, dimension(3) , intent(in) :: Sk        !  Complex vector
    !!----    complex, dimension(3)              :: Skp       !  Transformed complex vector
    !!----
    !!----    Apply a magnetic symmetry operator to a complex vector:  Skp = ApplyMSO(Op,Sk)
    !!----
    !!---- Update: April - 2005
    !!
    Function ApplyMSO(Op,Sk) Result(Skp)
       !---- Arguments ----!
       Type(MSym_Oper_Type), intent(in) :: Op
       Complex, dimension(3),intent(in) :: Sk
       Complex, dimension(3)            :: Skp

       Skp = matmul(Op%Rot,Sk) * cmplx(cos(tpi*Op%Phas),sin(tpi*Op%Phas))

       return
    End Function ApplyMSO


    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Init_Err_MagSym()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Init_Err_MagSym()

       err_magsym=.false.
       ERR_MagSym_Mess=" "

       return
    End Subroutine Init_Err_MagSym

    !!----
    !!---- subroutine Init_MagSymm_k_Type(MGp)
    !!----   type(MagSymm_k_Type),  intent (in out) :: MGp
    !!----
    !!----   Subroutine to initialize the MagSymm_k_Type variable MGp.
    !!----   It is called inside Readn_set_Magnetic_Structure
    !!----
    !!----  Update: April 2005
    !!
    Subroutine Init_MagSymm_k_Type(MGp)
       !---- Arguments ----!
       type(MagSymm_k_Type),  intent (in out) :: MGp

       !---- Local variables ----!
       integer :: i,j


       MGp%MagModel="Unnamed Model"
       MGp%Sk_Type="Crystal_Frame"       ! "Spherical_Frame"
       MGp%Latt="P"
       MGp%nmsym=0
       MGp%nirreps=0
       MGp%centred=1    !By default the crystal structure is acentric
       MGp%mcentred=1   !By default the magnetic structure is anti-centric (if there is -1 it is combined with time inversion)
       MGp%nkv=0
       MGp%kvec=0.0
       MGp%NumLat=1
       MGp%Ltr=0.0
       MGp%Numops=0
       MGp%Multip=0
       MGp%nbas=0
       MGp%icomp=0
       MGp%basf=cmplx(0.0,0.0)
       do i=1,48
          MGp%SymopSymb(i)=" "
          MGp%SymOp(i)%Rot(:,:)=0
          MGp%SymOp(i)%tr(:)=0.0
          do j=1,8
             MGp%MSymopSymb(i,j)=" "
             MGp%MSymOp(i,j)%Rot(:,:)=0
             MGp%MSymOp(i,j)%Phas=0.0
          end do
       end do

       return
    End Subroutine Init_MagSymm_k_Type

    !!---- Subroutine Init_Magnetic_Space_Group_Type(MGp)
    !!----   type(Magnetic_Space_Group_Type),  intent (in out) :: MGp
    !!----
    !!----   Initialize the non-allocatle parts of Magnetic_Space_Group_Type MGp
    !!----   It is called inside Readn_set_Magnetic_Structure
    !!----
    !!----   Updated: January-2014
    !!
    Subroutine Init_Magnetic_Space_Group_Type(MGp)
       !---- Arguments ----!
       type(Magnetic_Space_Group_Type),  intent (in out) :: MGp

       !---- Local variables ----!

       MGp%Sh_number=0
       MGp%BNS_number=" "
       MGp%OG_number=" "
       MGp%BNS_symbol=" "
       MGp%OG_symbol=" "
       MGp%MagType=0
       MGp%Parent_num=0
       MGp%Parent_spg=" "
       MGp%standard_setting=.false.
       MGp%mcif=.true.
       MGp%m_cell=.false.
       MGp%m_constr=.false.
       MGp%trn_from_parent=" "
       MGp%trn_to_standard=" "
       MGp%n_sym=0
       MGp%n_wyck=0
       MGp%n_kv=0
       return
    End Subroutine Init_Magnetic_Space_Group_Type

    Subroutine MagSymm_k_Type_to_Magnetic_Space_Group_Type(MG_Symk,MSpG)
       Type(MagSymm_k_Type),              intent(in)  :: MG_Symk
       Type(Magnetic_Space_Group_Type),   intent(out) :: MSpG
       !---- Local variables ----!
       Type(Space_Group_Type) :: SpG
       integer :: i,j,k,L,m,n, ngen
       integer,      dimension(5)    :: pos
       real(kind=cp)                 :: ph
       character(len=40),dimension(:), allocatable   :: gen
       character(len=132)   :: lowline,line
       character(len=30)    :: magmod, shubk
       character(len=2)     :: lattice, chardom
       character(len=4)     :: symbcar

       !
       call Init_Magnetic_Space_Group_Type(MSpG)

       !Verify the crystal structure information contained in MG_Symk by constructing the full Space group
       n=MG_Symk%Numops
       m=MG_Symk%Numops*MG_Symk%centred*MG_Symk%NumLat
       ngen=n-1
       allocate(gen(ngen))
       ngen=0
       do i=2,MG_Symk%Numops
         ngen=ngen+1
         gen(ngen)=MG_Symk%SymopSymb(i)
       end do
       if(MG_Symk%centred == 2) then
         ngen=ngen+1
         gen(ngen)="-x,-y,-z"
       end if
       Select Case(MG_Symk%Latt)
       End Select
       if(MG_Symk%Latt == "A") then
           ngen=ngen+1
           call Get_SymSymb(MG_Symk%SymOp(1)%rot,MG_Symk%Ltr(:,i),gen(ngen))
           gen(ngen)="-x,-y,-z"
       end if
       !Still to be finished
       call Set_SpaceGroup(" ",SpG,gen,ngen,"gen")


       return
    End Subroutine MagSymm_k_Type_to_Magnetic_Space_Group_Type

    !!----
    !!---- Subroutine Readn_Set_Magnetic_Structure_CFL(file_cfl,n_ini,n_end,MGp,Am,SGo,Mag_dom,Cell)
    !!----    type(file_list_type),                intent (in)     :: file_cfl
    !!----    integer,                             intent (in out) :: n_ini, n_end
    !!----    type(MagSymm_k_Type),                intent (out)    :: MGp
    !!----    type(mAtom_List_Type),               intent (out)    :: Am
    !!----    type(Magnetic_Group_Type), optional, intent (out)    :: SGo
    !!----    type(Magnetic_Domain_type),optional, intent (out)    :: Mag_dom
    !!----    type(Crystal_Cell_type),   optional, intent (in)     :: Cell
    !!----
    !!----    Subroutine for reading and construct the MagSymm_k_Type variable MGp.
    !!----    It is supposed that the CFL file is included in the file_list_type
    !!----    variable file_cfl. On output n_ini, n_end hold the lines with the
    !!----    starting and ending lines with information about a magnetic phase.
    !!----    Optionally the Magnetig space group (Shubnikov group) may be obtained
    !!----    separately for further use.
    !!----    Magnetic S-domains are also read in case of providing the optional variable Mag_dom.
    !!----
    !!---- Updates: November-2006, December-2011, July-2012 (JRC)
    !!
    Subroutine Readn_Set_Magnetic_Structure_CFL(file_cfl,n_ini,n_end,MGp,Am,SGo,Mag_dom,Cell)
       !---- Arguments ----!
       type(file_list_type),                intent (in)     :: file_cfl
       integer,                             intent (in out) :: n_ini, n_end
       type(MagSymm_k_Type),                intent (out)    :: MGp
       type(mAtom_List_Type),               intent (out)    :: Am
       type(Magnetic_Group_Type), optional, intent (out)    :: SGo
       type(Magnetic_Domain_type),optional, intent (out)    :: Mag_dom
       type(Crystal_Cell_type),   optional, intent (in)     :: Cell

       !---- Local Variables ----!
       integer :: i,no_iline,no_eline, num_k, num_xsym, num_irrep, num_dom, num_defdom, &
                  num_msym, ier, j, m, n, num_matom, num_skp, ik,im, ip, ncar
       integer,      dimension(5)    :: pos
       real(kind=cp)                 :: ph
       real(kind=cp),dimension(3)    :: rsk,isk,car,side
       real(kind=cp),dimension(3,12) :: br,bi
       real(kind=cp),dimension(3,3)  :: cart_to_cryst
       real(kind=cp),dimension(12)   :: coef
       character(len=132)   :: lowline,line
       character(len=30)    :: magmod, shubk
       character(len=2)     :: lattice, chardom
       character(len=4)     :: symbcar
       character(len=30)    :: msyr
       logical              :: msym_begin, kvect_begin, skp_begin, shub_given, irreps_given, &
                               irreps_begin, bfcoef_begin, magdom_begin
       type(Magnetic_Group_Type)  :: SG

       call init_err_MagSym()

       if(n_ini == 0) n_ini=1
       if(n_end == 0) n_end= file_cfl%nlines

       no_iline=0
       no_eline=0

       if(present(Cell)) then
         side(:)=Cell%cell
         cart_to_cryst=Cell%Orth_Cr_cel
       end if

       do i=n_ini,n_end
          ! Read comment
          if (index(file_cfl%line(i)(1:1),"!")/=0 .or. index(file_cfl%line(i)(1:1),"#")/=0) cycle
          lowline=adjustl(l_case(file_cfl%line(i)))

          if (lowline(1:13) == "mag_structure") then
             no_iline=i
          end if
          if (lowline(1:7) =="end_mag" ) then
             no_eline=i
             exit
          end if
       end do
       n_ini=no_iline
       n_end=no_eline

       if (n_ini == 0 .or. n_end == 0) then
          err_magsym=.true.
          ERR_MagSym_Mess=" No magnetic phase found in file!"
          return
       end if

       num_matom=0
       do i=n_ini,n_end
          lowline=l_case(adjustl(file_cfl%line(i)))
          if (index(lowline(1:5),"matom") ==0 ) cycle
          num_matom=num_matom+1
       end do

       Call Allocate_mAtom_list(num_matom,Am)  !Am contains Am%natoms = num_matom
       num_matom=0

       num_k=0
       num_dom=0
       num_defdom=0
       num_xsym=0
       kvect_begin=.true.
       magdom_begin=.true.
       call Init_MagSymm_k_Type(MGp)
       i=n_ini
       shub_given  =.false.
       irreps_given=.false.
       irreps_begin=.false.
       msym_begin  =.false.
       skp_begin   =.false.
       bfcoef_begin=.false.
       if (present(mag_dom)) then  !Initialise Mag_dom
          Mag_dom%nd=1
          Mag_dom%Chir=.false.
          Mag_dom%Twin=.false.
          Mag_dom%trans=.false.
          Mag_dom%DMat=0
          do j=1,3
           Mag_dom%DMat(j,j,1)=1
          end do
          Mag_dom%Dt=0.0
          Mag_dom%pop=0.0
          Mag_dom%pop(1,1)=1.0 !one domain is always present
          Mag_dom%Lab=" "
       end if

       do
          i=i+1
          if(i >= n_end) exit

          ! Read comment
          if( len_trim(file_cfl%line(i)) == 0) cycle
          lowline=adjustl(l_case(file_cfl%line(i)))
          if (lowline(1:1) == "!" .or. lowline(1:1)=="#") cycle

          ! Detect keywords

          ! Read magnetic model
          ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
          if (lowline(1:6) == "magmod") then
             read(unit=lowline(8:),fmt=*,iostat=ier) magmod
             if (ier /= 0) then
                err_magsym=.true.
                ERR_MagSym_Mess=" Error reading magnetic model name in magnetic phase"
                return
             end if
             MGp%MagModel= adjustl(magmod)
             cycle
          end if

          ! Read lattice
          if (lowline(1:7) == "lattice") then
             read(unit=lowline(9:),fmt=*,iostat=ier) lattice
             if (ier /= 0) then
                err_magsym=.true.
                ERR_MagSym_Mess=" Error reading lattice type in magnetic phase"
                return
             end if
             lattice=adjustl(lattice)
             if (lattice(1:1)=="-") then
                MGp%centred = 2
                MGp%latt=u_case(lattice(2:2))
             else
                MGp%centred = 1
                MGp%Latt= u_case(lattice(1:1))
             end if
             cycle
          end if

          ! Read type of Fourier coefficients
          if (lowline(1:9) == "spherical") then
             if(.not. present(Cell)) then
               err_magsym=.true.
               ERR_MagSym_Mess=" Cell argument is needed when Spherical components are used for Fourier Coefficients!"
             end if
             MGp%Sk_type = "Spherical_Frame"
             cycle
          end if

          ! Read magnetic centrig
          if (lowline(1:7) == "magcent") then
             MGp%mcentred = 2
             cycle
          end if

          ! Read propagation vectors
          if (lowline(1:5) == "kvect" .and. kvect_begin) then
             num_k=num_k+1
             read(unit=lowline(6:),fmt=*,iostat=ier) MGp%kvec(:,num_k)
             if (ier /= 0) then
                err_magsym=.true.
                ERR_MagSym_Mess=" Error reading propagation vectors"
                return
             end if
             do !repeat reading until continuous KVECT lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:5) == "kvect") then
                   num_k=num_k+1
                   read(unit=lowline(6:),fmt=*,iostat=ier) MGp%kvec(:,num_k)
                   if (ier /= 0) then
                      err_magsym=.true.
                      ERR_MagSym_Mess=" Error reading propagation vectors"
                      return
                   end if
                else
                   i=i-1
                   kvect_begin=.false.
                   exit
                end if
             end do
             cycle
          end if

          ! Read magnetic S-domains
          if (present(mag_dom)) then
             if (lowline(1:6) == "magdom" .and. magdom_begin) then
                num_dom=num_dom+1
                num_defdom=num_defdom+1
                if(index(lowline,"twin") /= 0) Mag_Dom%twin=.true.
                ip=index(lowline,":")
                if(index(lowline,"magdomt") == 0) then
                  msyr=lowline(8:ip-1)
                  call read_msymm(msyr,Mag_Dom%Dmat(:,:,num_dom),ph)
                  Mag_Dom%Dt(:,num_dom)=0.0
                  Mag_Dom%trans=.false.
                else
                  msyr=lowline(9:ip-1)
                  Call Get_Separator_Pos(msyr,",",pos,ncar)
                  if(ncar == 3) then
                    read(unit=msyr(pos(3)+1:),fmt=*,iostat=ier) ph
                    if(ier /= 0) ph=0.0
                    msyr=msyr(1:pos(3)-1)
                  else
                    ph=0.0
                  end if
                  call read_xsym(msyr,1,Mag_Dom%Dmat(:,:,num_dom),Mag_Dom%Dt(:,num_dom))
                  Mag_Dom%Dt(:,num_dom)=0.0
                  Mag_Dom%trans=.true.
                end if
                if (ph > 0.001) then
                  Mag_Dom%chir=.true.
                else
                  Mag_Dom%chir=.false.
                end if
                 if (Mag_Dom%chir) then
                   read(unit=lowline(ip+1:),fmt=*, iostat=ier) Mag_Dom%Pop(1:2,num_dom)
                   write(chardom,"(i2.2)") num_defdom
                   Mag_Dom%Lab(1,num_dom)="magdom"//chardom
                   num_defdom=num_defdom+1
                   write(chardom,"(i2.2)") num_defdom
                   Mag_Dom%Lab(2,num_dom)="magdom"//chardom
                else
                   read(unit=lowline(ip+1:),fmt=*, iostat=ier) Mag_Dom%Pop(1,num_dom)  !, Mag_Dom%MPop(1,num_dom)
                   write(chardom,"(i2.2)") num_defdom
                   Mag_Dom%Lab(1,num_dom)="magdom"//chardom
                end if
                if (ier /= 0) then
                   err_magsym=.true.
                   ERR_MagSym_Mess=" Error reading magnetic S-domains"
                   return
                end if
                Mag_Dom%nd = num_dom

                do  !repeat reading until continuous MAGDOM lines are exhausted
                   i=i+1
                   lowline=adjustl(l_case(file_cfl%line(i)))
                   ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                   if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                   if (lowline(1:6) == "magdom") then
                      if(index(lowline,"twin") /= 0) Mag_Dom%twin=.true.
                      num_dom=num_dom+1
                      num_defdom=num_defdom+1
                      ip=index(lowline,":")
                      if(index(lowline,"magdomt") == 0) then
                        msyr=lowline(8:ip-1)
                        call read_msymm(msyr,Mag_Dom%Dmat(:,:,num_dom),ph)
                        Mag_Dom%Dt(:,num_dom)=0.0
                        Mag_Dom%trans=.false.
                      else
                        msyr=lowline(9:ip-1)
                        Call Get_Separator_Pos(msyr,",",pos,ncar)
                        if(ncar == 3) then
                          read(unit=msyr(pos(3)+1:),fmt=*,iostat=ier) ph
                          if(ier /= 0) ph=0.0
                          msyr=msyr(1:pos(3)-1)
                        else
                          ph=0.0
                        end if
                        call read_xsym(msyr,1,Mag_Dom%Dmat(:,:,num_dom),Mag_Dom%Dt(:,num_dom))
                        Mag_Dom%Dt(:,num_dom)=0.0
                        Mag_Dom%trans=.true.
                      end if
                      if (ph > 0.001) then
                        Mag_Dom%chir=.true.
                      else
                         Mag_Dom%chir=.false.
                      end if
                      if (Mag_Dom%chir) then
                         read(unit=lowline(ip+1:),fmt=*, iostat=ier) Mag_Dom%Pop(1:2,num_dom) !, Mag_Dom%MPop(1:2,num_dom)
                         write(chardom,"(i2.2)") num_defdom
                         Mag_Dom%Lab(1,num_dom)="magdom"//chardom
                         num_defdom=num_defdom+1
                         write(chardom,"(i2.2)") num_defdom
                         Mag_Dom%Lab(2,num_dom)="magdom"//chardom
                      else
                         read(unit=lowline(ip+1:),fmt=*, iostat=ier) Mag_Dom%Pop(1,num_dom) !, Mag_Dom%MPop(1,num_dom)
                         write(chardom,"(i2.2)") num_defdom
                         Mag_Dom%Lab(1,num_dom)="magdom"//chardom
                      end if
                      if (ier /= 0) then
                         err_magsym=.true.
                         ERR_MagSym_Mess=" Error reading magnetic S-domains"
                         return
                      end if
                      Mag_Dom%nd = num_dom
                   else
                      i=i-1
                      magdom_begin=.false.
                      exit
                   end if
                end do
                cycle
             end if
          end if

          ! Read number of irreducible representations and number of basis functions for each
          if (lowline(1:6) == "irreps") then
             read(unit=lowline(7:),fmt=*,iostat=ier) MGp%nirreps
             if (ier /= 0) then
                err_magsym=.true.
                ERR_MagSym_Mess=" Error reading number of irreducible representations"
                return
             end if
             read(unit=lowline(7:),fmt=*,iostat=ier) n, (MGp%nbas(j),j=1,MGp%nirreps)
             if (ier /= 0) then
                err_magsym=.true.
                ERR_MagSym_Mess=" Error reading number of basis functions of irreducible representations"
                return
             end if
             irreps_given=.true.
             cycle
          end if

          ! Read the indicator real(0)/imaginary(1) of coefficients for basis functions of
          ! irreducible representations
          if (lowline(1:5) == "icomp" .and. irreps_given) then
             num_irrep=1
             n=MGp%nbas(num_irrep)
             read(unit=lowline(6:),fmt=*,iostat=ier) MGp%icomp(1:abs(n),num_irrep)
             if (ier /= 0) then
                err_magsym=.true.
                ERR_MagSym_Mess=" Error reading real/imaginary indicators of BF coeff. of irreducible representations"
                return
             end if
             do  !repeat reading until continuous icoebf lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:5) == "icomp") then
                   num_irrep=num_irrep+1
                   n=MGp%nbas(num_irrep)
                   read(unit=lowline(6:),fmt=*,iostat=ier) MGp%icomp(1:abs(n),num_irrep)
                   if (ier /= 0) then
                      err_magsym=.true.
                      ERR_MagSym_Mess=" Error reading real/imaginary indicators of BF coeff. of irreducible representations"
                      return
                   end if
                else
                   i=i-1
                   irreps_given=.false.
                   exit
                end if
             end do
             cycle
          end if

          ! Read Shubnikov group
          if (lowline(1:9) == "shubnikov") then
             shubk=adjustl(file_cfl%line(i)(10:))
             Call Set_Shubnikov_Group(shubk,SG,MGp)
             if (err_magsym) return
             shub_given=.true.
          end if

          ! Read SYMM operators
          if (lowline(1:4) == "symm" .and. .not. shub_given) then
             num_xsym=num_xsym+1
             num_msym=0
             num_irrep=0
             read(unit=lowline(5:),fmt="(a)") MGp%SymopSymb(num_xsym)
             msym_begin=.true.
             irreps_begin=.true.
          end if

          ! Read MSYM operators
          if (lowline(1:4) == "msym" .and. msym_begin .and. .not. shub_given) then
             num_msym=num_msym+1
             read(unit=lowline(5:),fmt="(a)") MGp%MSymopSymb(num_xsym,num_msym)
             do  !repeat reading until continuous MSYM lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:4) == "msym") then
                   num_msym=num_msym+1
                   read(unit=lowline(5:),fmt="(a)") MGp%MSymopSymb(num_xsym,num_msym)
                else
                   i=i-1
                   msym_begin=.false.
                   exit
                end if
             end do
             cycle
          end if

          ! Read basis functions of irreducible representations
          if (lowline(1:4) == "basr" .and. irreps_begin .and. .not. shub_given) then
             num_irrep=num_irrep+1
             n=MGp%nbas(num_irrep)
             br=0.0; bi=0.0
             read(unit=lowline(5:),fmt=*,iostat=ier) (br(:,j),j=1,abs(n))
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=ERR_MagSym_Mess,fmt="(2(a,i3))")" Error reading basis fuctions (BASR) of irrep ",num_irrep,&
                                                           " for symmetry operator # ",num_xsym
                return
             end if
             if (n < 0) then  !Read the imaginary part of the basis functions
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:4) == "basi") then
                   read(unit=lowline(5:),fmt=*,iostat=ier) (bi(:,j),j=1,abs(n))
                   if (ier /= 0) then
                      err_magsym=.true.
                      write(unit=ERR_MagSym_Mess,fmt="(2(a,i3))")" Error reading basis fuctions (BASI) of irrep ",num_irrep,&
                                                                 " for symmetry operator # ",num_xsym
                      return
                   end if
                else
                   err_magsym=.true.
                   write(unit=ERR_MagSym_Mess,fmt="(2(a,i3))")" Lacking BASI keyword of irrep ",num_irrep,&
                                                               " for symmetry operator # ",num_xsym
                   return
                end if
             end if
             do j=1,abs(n)
                MGp%basf(:,j,num_xsym,num_irrep)=cmplx( br(:,j),bi(:,j) )
             end do

             do  !repeat reading until continuous BASR or BASI lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:4) == "basr") then
                   num_irrep=num_irrep+1
                   n=MGp%nbas(num_irrep)
                   br=0.0; bi=0.0
                   read(unit=lowline(5:),fmt=*,iostat=ier) (br(:,j),j=1,abs(n))
                   if (ier /= 0) then
                      err_magsym=.true.
                      write(unit=ERR_MagSym_Mess,fmt="(2(a,i3))")" Error reading basis fuctions (BASR) of irrep ",num_irrep,&
                                                                 " for symmetry operator # ",num_xsym
                      return
                   end if
                   if (n < 0) then  !Read the imaginary part of the basis functions
                      i=i+1
                      lowline=adjustl(l_case(file_cfl%line(i)))
                      !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                      if (lowline(1:4) == "basi") then
                         read(unit=lowline(5:),fmt=*,iostat=ier) (bi(:,j),j=1,abs(n))
                         if (ier /= 0) then
                            err_magsym=.true.
                            write(unit=ERR_MagSym_Mess,fmt="(2(a,i3))")" Error reading basis fuctions (BASI) of irrep ",num_irrep,&
                                                                       " for symmetry operator # ",num_xsym
                            return
                         end if
                      else
                         err_magsym=.true.
                         write(unit=ERR_MagSym_Mess,fmt="(2(a,i3))")" Lacking BASI keyword of irrep ",num_irrep,&
                                                                    " for symmetry operator # ",num_xsym
                         return
                      end if
                   end if
                   do j=1,abs(n)
                      MGp%basf(:,j,num_xsym,num_irrep)=cmplx( br(:,j),bi(:,j) )
                   end do
                else
                   i=i-1
                   irreps_begin=.false.
                   exit
                end if
             end do
             cycle
          end if

          ! Read magnetic atoms:  label, magnetic form factor label,x,y,z,Biso,occ
          if (lowline(1:5) == "matom") then
             num_matom=num_matom+1
             num_skp=0
             line=adjustl(file_cfl%line(i))
             read(unit=line(6:),fmt=*,iostat=ier) Am%atom(num_matom)%lab,      & !Label
                                                  Am%atom(num_matom)%SfacSymb, & !Formfactor label
                                                  Am%atom(num_matom)%x,        & !Fract. coord.
                                                  Am%atom(num_matom)%Biso,     & !Is. Temp. Fact.
                                                  Am%atom(num_matom)%occ         !occupation
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=ERR_MagSym_Mess,fmt="(a,i4)")" Error reading magnetic atom #",num_matom
                return
             end if
             skp_begin=.true.
             bfcoef_begin=.true.
             cycle
          end if

          ! Read Fourier coefficients in cryst. axes and phase
          if (lowline(1:3) == "skp" .and. skp_begin) then
             num_skp=num_skp+1
             read(unit=lowline(4:),fmt=*,iostat=ier) ik,im,rsk,isk,ph
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=ERR_MagSym_Mess,fmt="(a,i3)") " Error reading Fourier Coefficient #", num_skp
                return
             end if
               Am%atom(num_matom)%nvk= num_skp
               Am%atom(num_matom)%imat(ik)= im
               Am%atom(num_matom)%mphas(ik)= ph

             if(MGp%Sk_type == "Spherical_Frame") then
               Am%atom(num_matom)%Spher_Skr(:,ik)= rsk(:)
               Am%atom(num_matom)%Spher_Ski(:,ik)= isk(:)
               !Transform from Cartesian coordinates to unitary Crystallographic frame
               call Get_Cart_from_Spher(rsk(1),rsk(3),rsk(2),car,"D")
               Am%atom(num_matom)%Skr(:,ik)=matmul(cart_to_cryst,car)*side(:)
               call Get_Cart_from_Spher(isk(1),isk(3),isk(2),car,"D")
               Am%atom(num_matom)%Ski(:,ik)=matmul(cart_to_cryst,car)*side(:)
             else  !In this case, the Cell argument may be not given
                   !so no transformation is done. This can be done in other parts of the calling program
               Am%atom(num_matom)%Skr(:,ik)= rsk(:)
               Am%atom(num_matom)%Ski(:,ik)= isk(:)
             end if

             do  !repeat reading until continuous SPK lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:3) == "skp") then
                   num_skp=num_skp+1
                   if (num_skp > 12) then
                      err_magsym=.true.
                      ERR_MagSym_Mess= " Too many Fourier Coefficients, the maximum allowed is 12! "
                      return
                   end if
                   read(unit=lowline(4:),fmt=*,iostat=ier) ik,im,rsk,isk,ph
                   if (ier /= 0) then
                      err_magsym=.true.
                      write(unit=ERR_MagSym_Mess,fmt="(a,i3)") " Error reading Fourier Coefficient #", num_skp
                      return
                   end if
                   Am%atom(num_matom)%nvk= num_skp
                   Am%atom(num_matom)%imat(ik)= im
                   Am%atom(num_matom)%Skr(:,ik)= rsk(:)
                   Am%atom(num_matom)%Ski(:,ik)= isk(:)
                   Am%atom(num_matom)%mphas(ik)= ph
                else
                   i=i-1
                   skp_begin=.false.
                   Am%atom(num_matom)%nvk= num_skp
                   exit
                end if
             end do
          end if

          if (lowline(1:6) == "bfcoef" .and. bfcoef_begin) then
             num_skp=num_skp+1
             read(unit=lowline(7:),fmt=*,iostat=ier) ik,im
             n=abs(MGp%nbas(im))
             read(unit=lowline(7:),fmt=*,iostat=ier) ik,im,coef(1:n),ph
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=ERR_MagSym_Mess,fmt="(a,i3)") " Error reading Coefficient of Basis Functions #", num_skp
                return
             end if
             Am%atom(num_matom)%nvk= num_skp
             Am%atom(num_matom)%imat(ik)= im
             Am%atom(num_matom)%cbas(1:n,ik)= coef(1:n)
             Am%atom(num_matom)%mphas(ik)= ph

             do  !repeat reading until continuous bfcoef lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:6) == "bfcoef" ) then
                   num_skp=num_skp+1
                   if (num_skp > 12) then
                      err_magsym=.true.
                      ERR_MagSym_Mess= " Too many sets of Coefficients, the maximum allowed is 12! "
                      return
                   end if
                   read(unit=lowline(7:),fmt=*,iostat=ier) ik,im
                   n=abs(MGp%nbas(im))
                   read(unit=lowline(7:),fmt=*,iostat=ier) ik,im,coef(1:n),ph
                   if (ier /= 0) then
                      err_magsym=.true.
                      write(unit=ERR_MagSym_Mess,fmt="(a,i3)") " Error reading Coefficient of Basis Functions #", num_skp
                      return
                   end if
                   Am%atom(num_matom)%nvk= num_skp
                   Am%atom(num_matom)%imat(ik)= im
                   Am%atom(num_matom)%cbas(1:n,ik)= coef(1:n)
                   Am%atom(num_matom)%mphas(ik)= ph
                else
                   i=i-1
                   bfcoef_begin=.false.
                   Am%atom(num_matom)%nvk= num_skp
                   exit
                end if
             end do
          end if
       end do

       !Arriving here we have exhausted reading magnetic phase

       !Get pointers to the magnetic form factors
       !Stored for each atom in the component ind(1)
       call Set_Magnetic_Form()

       !---- Find Species in Magnetic_Form ----!
       do i=1,Am%natoms
          symbcar=u_case(Am%atom(i)%SfacSymb)
          do j=1,num_mag_form
             if (symbcar /= Magnetic_Form(j)%Symb) cycle
             Am%atom(i)%ind(1)=j
             exit
          end do
       end do

       !Now construct the rest of magnetic symmetry type variable MGp
       MGp%nmsym =num_msym
       MGp%Numops=num_xsym
       MGp%nkv   =num_k

       !Construct the numerical symmetry operators
       do i=1,MGp%Numops
          Call Read_Xsym(MGp%SymopSymb(i),1,MGp%Symop(i)%Rot,MGp%Symop(i)%tr)
          do j=1,MGp%nmsym
             Call Read_Msymm(MGp%MSymopSymb(i,j),MGp%MSymop(i,j)%Rot,MGp%MSymop(i,j)%Phas)
          end do
       end do
       if (err_symm) then
          err_magsym=.true.
          write(unit=ERR_MagSym_Mess,fmt="(a)") " Error reading symmetry: "//trim(err_symm_mess)
          return
       end if

       !Complete the set of symmetry operators with the centre of symmetry
       m=MGp%Numops
       if (MGp%centred == 2) then
          do i=1,MGp%Numops
             m=m+1
             MGp%Symop(m)%Rot(:,:) = -MGp%Symop(i)%Rot(:,:)
             MGp%Symop(m)%tr(:)    =  modulo_lat(-MGp%Symop(m)%tr(:))
             call Get_SymSymb(MGp%Symop(m)%Rot(:,:), &
                              MGp%Symop(m)%tr(:), MGp%SymopSymb(m))
             if (Mgp%mcentred == 1) then  !Anticentre in the magnetic structure
                do j=1,MGp%nmsym
                   MGp%MSymop(m,j)%Rot(:,:) = -MGp%MSymop(i,j)%Rot(:,:)
                   MGp%MSymop(m,j)%Phas     = -MGp%MSymop(i,j)%Phas
                end do
             else if(Mgp%mcentred == 2) then
                do j=1,MGp%nmsym
                   MGp%MSymop(m,j)%Rot(:,:) =  MGp%MSymop(i,j)%Rot(:,:)
                   MGp%MSymop(m,j)%Phas     =  MGp%MSymop(i,j)%Phas
                end do
             end if
          end do
       end if

       !Get the centring lattice translations of the crystallographic structure
       !and calculate the general multiplicity of the group.
       Mgp%NumLat=1
       MGp%Ltr(:,:) = 0.0
       Select Case(MGp%Latt)
          case ("A")
             Mgp%NumLat=2
             MGp%Ltr(:,1:2)=Ltr_a(:,1:2)
          case ("B")
             Mgp%NumLat=2
             MGp%Ltr(:,1:2)=Ltr_b(:,1:2)
          case ("C")
             Mgp%NumLat=2
             MGp%Ltr(:,1:2)=Ltr_c(:,1:2)
          case ("I")
             Mgp%NumLat=2
             MGp%Ltr(:,1:2)=Ltr_i(:,1:2)
          case ("R")
             Mgp%NumLat=3
             MGp%Ltr(:,1:3)=Ltr_r(:,1:3)
          case ("F")
             Mgp%NumLat=4
             MGp%Ltr(:,1:4)=Ltr_f(:,1:4)
       End Select

       select case (MGp%centred)
          case (1)
             MGp%Multip =   MGp%Numops * Mgp%NumLat
          case (2)
             MGp%Multip = 2 * MGp%Numops * Mgp%NumLat
       end select

       if (present(SGo)) then
          if (shub_given) then
             SGo=SG
          else
             err_magsym=.true.
             ERR_MagSym_Mess=" Shubnikov Group has not been provided "
          end if
       end if

       return
    End Subroutine Readn_Set_Magnetic_Structure_CFL


    !!----
    !!---- Subroutine Readn_Set_Magnetic_Structure_MCIF(file_mcif,mCell,MGp,Am)
    !!----    character(len=*),        intent (in)     :: file_mcif
    !!----    type(Crystal_Cell_type), intent (out)    :: mCell
    !!----    type(MagSymm_k_Type),    intent (out)    :: MGp
    !!----    type(mAtom_List_Type),   intent (out)    :: Am
    !!----
    !!----    Subroutine for reading and construct the MagSymm_k_Type variable MGp.
    !!----    The magnetic atom list and the magnetic cell reading an mCIF file.
    !!----
    !!----  Created: January-2014 (JRC)
    !!
    Subroutine Readn_Set_Magnetic_Structure_MCIF(file_mcif,mCell,MGp,Am)
       character(len=*),               intent (in)  :: file_mcif
       type(Crystal_Cell_type),        intent (out) :: mCell
       type(Magnetic_Space_Group_Type),intent (out) :: MGp
       type(mAtom_List_Type),          intent (out) :: Am

       !---- Local Variables ----!
       integer :: i,num_sym, num_constr, num_kvs,num_msym,num_matom, num_mom,   &
                  ier, j, m, n, k, ncar,mult,nitems,iv, num_irreps, nitems_irreps
       integer,   dimension(9)             :: lugar
       integer,   dimension(6)             :: irrep_pos
       integer,   dimension(5)             :: pos
       integer,   dimension(3)             :: code
       real(kind=cp)                       :: ph
       real(kind=cp),dimension(3)          :: cel,ang,cel_std,ang_std
       real(kind=cp),dimension(6)          :: values,std
       real(kind=cp),dimension(3,3)        :: matr
       real(kind=cp),dimension(3,384)      :: orb
       character(len=132)                  :: lowline,keyword,line
       character(len=132),dimension(384)   :: sym_strings
       character(len=132),dimension(384)   :: atm_strings
       character(len=132),dimension(384)   :: mom_strings
       character(len=132),dimension(30)    :: constr_strings
       character(len=132),dimension(30)    :: irreps_strings
       character(len=132),dimension(30)    :: kv_strings
       character(len=20), dimension(15)    :: lab_items
       character(len=40)    :: shubk
       character(len=2)     :: chars
       character(len=10)    :: label
       character(len=4)     :: symbcar
       logical              :: ktag

       type(Magnetic_Group_Type)  :: SG
       type(file_list_type)       :: mcif

       call init_err_MagSym()

       call File_To_FileList(file_mcif,mcif)
       !Remove all possible tabs and non-ASCII characters in the CIF
       do i=1,mcif%nlines
         do j=1,len_trim(mcif%line(i))
           if(mcif%line(i)(j:j) == char(9)) mcif%line(i)(j:j)=" "
         end do
       end do
       num_constr=0; num_kvs=0; num_matom=0; num_mom=0; num_sym=0
       cel=0.0; ang=0.0
       i=0
       call Init_Magnetic_Space_Group_Type(MGp)
       ktag=.false.
       do
          i=i+1
          if(i > mcif%nlines) exit
          if (index(mcif%line(i)(1:1),"!")/=0 .or. index(mcif%line(i)(1:1),"#")/=0 .or. len_trim(mcif%line(i)) == 0) cycle
          line=adjustl(mcif%line(i))
          lowline=l_case(line)
          j=index(lowline," ")
          keyword=lowline(1:j-1)

          Select Case (trim(keyword))

             Case("_magnetic_space_group_standard_setting")
                chars=adjustl(line(j+1:))
                if(chars(2:2) == "y" .or. chars(2:2) == "Y") MGp%standard_setting=.true.

             Case("_parent_space_group.name_h-m")
                shubk=adjustl(line(j+1:))
                m=len_trim(shubk)
                MGp%Parent_spg=shubk(2:m-1)

             Case("_parent_space_group.it_number")
                read(unit=lowline(j:),fmt=*,iostat=ier) m
                if(ier /= 0) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the number of the parent space group"
                  return
                end if
                MGp%Parent_num=m

             Case("_magnetic_space_group_bns_number")
                shubk=adjustl(line(j+1:))
                MGp%BNS_number=shubk

             Case("_magnetic_space_group_bns_name")
                shubk=adjustl(line(j+1:))
                m=len_trim(shubk)
                MGp%BNS_symbol=shubk(2:m-1)

             Case("_magnetic_space_group_og_number")
                shubk=adjustl(line(j+1:))
                MGp%OG_number=shubk

             Case("_magnetic_space_group_og_name")
                shubk=adjustl(line(j+1:))
                m=len_trim(shubk)
                MGp%OG_symbol=shubk(2:m-1)

             Case("_magnetic_space_group.transform_from_parent_pp_abc")
                shubk=adjustl(line(j+1:))
                m=len_trim(shubk)
                MGp%trn_from_parent=shubk(2:m-1)

             Case("_magnetic_space_group.transform_to_standard_pp_abc")
                shubk=adjustl(line(j+1:))
                m=len_trim(shubk)
                MGp%trn_to_standard=shubk(2:m-1)

             Case("_magnetic_cell_length_a")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'a' -> "//trim(err_string_mess)
                  return
                end if
                cel(1)=values(1)
                cel_std(1)=std(1)
                MGp%m_cell=.true.

             Case("_magnetic_cell_length_b")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'b' -> "//trim(err_string_mess)
                  return
                end if
                cel(2)=values(1)
                cel_std(2)=std(1)

             Case("_magnetic_cell_length_c")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'c' -> "//trim(err_string_mess)
                  return
                end if
                cel(3)=values(1)
                cel_std(3)=std(1)

             Case("_magnetic_cell_angle_alpha")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'alpha' -> "//trim(err_string_mess)
                  return
                end if
                ang(1)=values(1)
                ang_std(1)=std(1)

             Case("_magnetic_cell_angle_beta")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'beta' -> "//trim(err_string_mess)
                  return
                end if
                ang(2)=values(1)
                ang_std(2)=std(1)

             Case("_magnetic_cell_angle_gamma")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'gamma' -> "//trim(err_string_mess)
                  return
                end if
                ang(3)=values(1)
                ang_std(3)=std(1)

             Case("loop_")
                 i=i+1
                 line=adjustl(mcif%line(i))
                 lowline=l_case(line)
                 j=index(lowline," ")
                 keyword=lowline(1:j-1)

                 Select Case(trim(keyword))

                   Case("_irrep_id")
                      irrep_pos=0
                      irrep_pos(1)=1
                      j=1
                      do k=1,6
                         i=i+1
                         if(index(mcif%line(i),"_irrep_dimension") /= 0) then
                            j=j+1
                            irrep_pos(2)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_small_irrep_dimension") /= 0) then
                            j=j+1
                            irrep_pos(3)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_direction_type") /= 0) then
                            j=j+1
                            irrep_pos(4)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_action") /= 0) then
                            j=j+1
                            irrep_pos(5)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_modes_number") /= 0) then
                            j=j+1
                            irrep_pos(6)=j
                            cycle
                         end if
                         exit
                      end do

                      i=i-1
                      nitems_irreps=count(irrep_pos > 0)

                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        irreps_strings(k)=mcif%line(i)
                      end do
                      num_irreps=k
                      !Treat later the list of irreps

                   Case("_magnetic_propagation_vector_seq_id")
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_magnetic_propagation_vector") == 0) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the propagation vector loop"
                          return
                        end if
                        if(index(mcif%line(i),"_magnetic_propagation_vector_kxkykz") /= 0) then
                          ktag=.true.  !new format for k-vector klabel '0,1/2,0'
                          exit
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        kv_strings(k)=mcif%line(i)
                      end do
                      num_kvs=k
                      MGp%n_kv=k
                      if(allocated(Mgp%kv)) deallocate(Mgp%kv)
                      allocate(Mgp%kv(3,k))
                      if(allocated(Mgp%kv_label)) deallocate(Mgp%kv_label)
                      allocate(Mgp%kv_label(k))
                      !Treat later the propagation vectors

                   Case("_magnetic_atom_site_moment_symmetry_constraints_label")
                      i=i+1
                      if(index(mcif%line(i),"_atom_site_magnetic_moment_symmetry_constraints_mxmymz") == 0) then
                        err_magsym=.true.
                        ERR_MagSym_Mess=" Error reading the magnetic_atom_site_moment_symmetry_constraints loop"
                        return
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        constr_strings(k)=mcif%line(i)
                      end do
                      num_constr=k
                      MGp%m_constr=.true.
                      !Treat later the constraints

                   Case("_magnetic_space_group_symop_id")
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_magnetic_space_group_symop_operation") == 0) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the magnetic_space_group_symop_operation loop"
                          return
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do
                      !now allocate the list of symmetry operators
                      num_sym=k
                      MGp%n_sym=k
                      if(allocated(Mgp%SymopSymb)) deallocate(Mgp%SymopSymb)
                      allocate(Mgp%SymopSymb(k))
                      if(allocated(Mgp%Symop)) deallocate(Mgp%Symop)
                      allocate(Mgp%Symop(k))
                      if(allocated(Mgp%MSymopSymb)) deallocate(Mgp%MSymopSymb)
                      allocate(Mgp%MSymopSymb(k))
                      if(allocated(Mgp%MSymop)) deallocate(Mgp%MSymop)
                      allocate(Mgp%MSymop(k))

                   Case("_magnetic_atom_site_label")
                      lugar=0
                      lugar(1)=1
                      j=1
                      do k=1,9
                         i=i+1
                         if(index(mcif%line(i),"_atom_site_type_symbol") /= 0) then
                            j=j+1
                            lugar(2)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_x") /= 0) then
                            j=j+1
                            lugar(3)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_y") /= 0) then
                            j=j+1
                            lugar(4)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_z") /= 0) then
                            j=j+1
                            lugar(5)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_U_iso_or_equiv") /= 0) then
                            j=j+1
                            lugar(6)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_occupancy") /= 0) then
                            j=j+1
                            lugar(7)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_symmetry_multiplicity") /= 0) then
                            j=j+1
                            lugar(8)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_Wyckoff_label") /= 0) then
                            j=j+1
                            lugar(9)=j
                            cycle
                         end if
                         exit
                      end do

                      if (any(lugar(3:5) == 0)) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the asymmetric unit of magnetic atoms"
                          return
                      end if

                      i=i-1
                      nitems=count(lugar > 0)

                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        atm_strings(k)=mcif%line(i)
                      end do
                      num_matom=k
                      !Treat late the list atoms

                   Case("_magnetic_atom_site_moment_label")
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_magnetic_atom_site_moment_crystalaxis") == 0) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the magnetic_atom_site_moment loop"
                          return
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        mom_strings(k)=mcif%line(i)
                      end do
                      num_mom=k
                      !Treat later the magnetic moment of the atoms
                 End Select
          End Select
       end do

       if(MGp%m_cell) then
         call Set_Crystal_Cell(cel,ang,mCell)
         mCell%cell_std=cel_std
         mCell%ang_std=ang_std
       end if

       !Treat symmetry operators
       if(num_sym == 0) then
          err_magsym=.true.
          ERR_MagSym_Mess=" No symmetry operators have been provided in the MCIF file "//trim(file_mcif)
          return
       else  !Decode the symmetry operators
         do i=1,num_sym
           line=adjustl(sym_strings(i))
           j=index(line," ")
           line=adjustl(line(j+1:))
           j=index(line," ")
           MGp%SymopSymb(i)=line(1:j-1)
           line=adjustl(line(j+1:))
           j=index(line," ")
           MGp%MSymopSymb(i)=line(1:j-1)
           read(unit=line(j:),fmt=*,iostat=ier) n
           if(ier /= 0) then
              err_magsym=.true.
              ERR_MagSym_Mess=" Error reading the time inversion in line: "//trim(sym_strings(i))
              return
           else
              MGp%MSymOp(i)%phas=real(n)
           end if
           call Read_Xsym(MGp%SymopSymb(i),1,MGp%Symop(i)%Rot,MGp%Symop(i)%tr)
           line=MGp%MSymopSymb(i)
           do k=1,len_trim(line)
             if(line(k:k) == "m") line(k:k)=" "
           end do
           line=Pack_String(line)
           call Read_Xsym(line,1,MGp%MSymop(i)%Rot)
         end do
       end if
       ! Symmetry operators treatment done!


       !Treating irreps

       if(num_irreps == 0) then

          MGp%n_irreps=0

        else

          MGp%n_irreps=num_irreps
          if(allocated(MGp%irrep_dim))          deallocate(MGp%irrep_dim)
          if(allocated(MGp%small_irrep_dim))    deallocate(MGp%small_irrep_dim)
          if(allocated(MGp%irrep_id))           deallocate(MGp%irrep_id)
          if(allocated(MGp%irrep_direction))    deallocate(MGp%irrep_direction)
          if(allocated(MGp%irrep_action))       deallocate(MGp%irrep_action)
          if(allocated(MGp%irrep_modes_number)) deallocate(MGp%irrep_modes_number)
          allocate(MGp%irrep_dim(num_irreps),MGp%small_irrep_dim(num_irreps),MGp%irrep_id(num_irreps), &
                   MGp%irrep_direction(num_irreps),MGp%irrep_action(num_irreps),MGp%irrep_modes_number(num_irreps))

          MGp%irrep_dim=0; MGp%small_irrep_dim=0; MGp%irrep_id=" "; MGp%irrep_direction=" "; MGp%irrep_action=" "
          MGp%irrep_modes_number=0

          do i=1,MGp%n_irreps

            call getword(irreps_strings(i),lab_items,iv)

            !if(iv /= nitems_irreps) write(*,"(2(a,i2))") " => Warning irreps_nitems=",nitems_irreps," /= items read=",iv

            MGp%irrep_id(i)=lab_items(irrep_pos(1))
            if(MGp%irrep_id(i) == "?") then
               MGp%n_irreps=0
               exit
            end if

            if (irrep_pos(2) /= 0) then
               read(unit=lab_items(irrep_pos(2)),fmt=*,iostat=ier) MGp%irrep_dim(i)
               if(ier /= 0) MGp%irrep_dim(i)=0
            end if

            if (irrep_pos(3) /= 0) then
               read(unit=lab_items(irrep_pos(3)),fmt=*,iostat=ier) MGp%small_irrep_dim(i)
               if(ier /= 0) MGp%small_irrep_dim(i)=0
            end if

            if (irrep_pos(4) /= 0) then
               MGp%irrep_direction(i)=lab_items(irrep_pos(4))
            end if

            if (irrep_pos(5) /= 0) then
               MGp%irrep_action(i)=lab_items(irrep_pos(5))
            end if

            if (irrep_pos(6) /= 0) then
               read(unit=lab_items(irrep_pos(6)),fmt=*,iostat=ier) MGp%irrep_modes_number(i)
               if(ier /= 0) MGp%irrep_modes_number(i)=0
            end if

          end do
       end if
       ! End treatment of irreps

       ! Treating propagation vectors
       if(num_kvs == 0) then
         MGp%n_kv=0
       else
         do i=1,MGp%n_kv
            line=adjustl(kv_strings(i))
            j=index(line," ")
            MGp%kv_label(i)=line(1:j-1)
            line=adjustl(line(j+1:))
            n=len_trim(line)
            if(ktag) then
              line=line(2:n-1)
              n=n-2
              Call Get_Separator_Pos(line,",",pos,ncar)
            else
              Call Get_Separator_Pos(line," ",pos,ncar)
            end if
            keyword=line(1:pos(1)-1)//"a,"//line(pos(1)+1:pos(2)-1)//"b,"//trim(line(pos(2)+1:))//"c"
            keyword=Pack_String(keyword)
            call Get_Mat_From_Symb(keyword,Matr, (/"a","b","c"/) )
            do k=1,3
               MGp%kv(k,i)=Matr(k,k)
            end do
         end do
       end if
       ! Propagation vectors treatment done!

       !Treating magnetic atoms
       if(num_matom == 0) then
          Am%natoms = 0
          return
       else
          Call Allocate_mAtom_list(num_matom,Am)

          do i=1,Am%natoms

            call getword(atm_strings(i),lab_items,iv)
            !if(iv /= nitems) write(*,"(2(a,i2))") " => Warning nitems=",nitems," /= items read=",iv
            Am%atom(i)%lab=lab_items(lugar(1))
            if (lugar(2) /= 0) then
               Am%atom(i)%SfacSymb=lab_items(lugar(2))(1:4)
               if(index("1234567890+-",lab_items(lugar(2))(2:2)) /= 0 ) then
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(2))(1:1))
               else
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(2))(1:1))//L_case(lab_items(lugar(2))(2:2))
               end if
            else
               if(index("1234567890+-",lab_items(lugar(1))(2:2)) /= 0 ) then
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(1))(1:1))
               else
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(1))(1:1))//L_case(lab_items(lugar(1))(2:2))
               end if
               Am%atom(i)%SfacSymb=Am%atom(i)%chemSymb
            end if
            call getnum_std(lab_items(lugar(3)),values,std,iv)    ! _atom_site_fract_x
            Am%atom(i)%x(1)=values(1)
            Am%atom(i)%x_std(1)=std(1)
            call getnum_std(lab_items(lugar(4)),values,std,iv)    ! _atom_site_fract_y
            Am%atom(i)%x(2)=values(1)
            Am%atom(i)%x_std(2)=std(1)
            call getnum_std(lab_items(lugar(5)),values,std,iv)    ! _atom_site_fract_z
            Am%atom(i)%x(3)=values(1)
            Am%atom(i)%x_std(3)=std(1)

            if (lugar(6) /= 0) then  ! _atom_site_Uiso_or_equiv
               call getnum_std(lab_items(lugar(6)),values,std,iv)
            else
               values=0.0
               std=0.0
            end if
            Am%atom(i)%ueq=values(1)
            Am%atom(i)%Biso=values(1)*78.95683521     !If anisotropic they
            Am%atom(i)%Biso_std=std(1)*78.95683521    !will be put to zero
            Am%atom(i)%utype="u_ij"

            if (lugar(7) /= 0) then ! _atom_site_occupancy
               call getnum_std(lab_items(lugar(7)),values,std,iv)
            else
               values=1.0
               std=0.0
            end if
            Am%atom(i)%occ=values(1)
            Am%atom(i)%occ_std=std(1)

            if(lugar(8) /= 0) then
              read(unit=lab_items(lugar(8)),fmt=*) Mult
              Am%atom(i)%mult=Mult
            else
              Call Get_mOrbit(Am%atom(i)%x,MGp,Mult,orb)
              Am%atom(i)%mult=Mult
            end if
            !Conversion from occupancy to occupation factor
            Am%atom(i)%occ=Am%atom(i)%occ*real(Mult)/real(MGp%n_sym)

            if(lugar(9) /= 0) then
               Am%atom(i)%wyck=adjustl(trim(lab_items(lugar(9))))
            end if

          end do
       end if

       !Treating moments of magnetic atoms
       if(num_mom /= 0) then
          do i=1,num_mom
            call getword(mom_strings(i),lab_items,iv)
            !write(*,"(4(a,tr3))")lab_items(1:iv)
            if(iv /= 4) then
               err_magsym=.true.
               write(unit=ERR_MagSym_Mess,fmt="(a,i4)")" Error reading magnetic moment #",i
               ERR_MagSym_Mess=trim(ERR_MagSym_Mess)//" -> 4 items expected in this line: 'Label mx my mz', read: "// &
                                                      trim(mom_strings(i))
               return
            end if
            label=Lab_items(1)
            do j=1,Am%natoms
               if(label == Am%Atom(j)%lab) then
                 do k=1,3
                     call getnum_std(lab_items(1+k),values,std,iv)
                     Am%Atom(j)%SkR(k,1)=values(1)
                     Am%Atom(j)%SkR_std(k,1)=std(1)
                 end do
               end if
            end do
          end do
       end if

       if(num_constr /= 0) then

         do i=1,num_constr
           line=adjustl(constr_strings(i))
           j=index(line," ")
           label=line(1:j-1)
           keyword=adjustl(line(j+1:))
           Call Get_Separator_Pos(keyword,",",pos,ncar)
           if(ncar == 0) then !There are no ","
             j=index(keyword," ")
             shubk=keyword(1:j-1)//","
             keyword=adjustl(keyword(j+1:))
             j=index(keyword," ")
             shubk=trim(shubk)//keyword(1:j-1)//","
             keyword=trim(shubk)//trim(adjustl(keyword(j+1:)))
           end if
           do j=1,len_trim(keyword)
             if(keyword(j:j) == "m") keyword(j:j) = " "
           end do
           keyword=Pack_String(keyword)
           !write(*,"(a)") "  constr_string: "//trim(line)
           !write(*,"(a)") "        keyword: "//trim(keyword)
           call Get_Mat_From_Symb(keyword,Matr, (/"x","y","z"/) )
           !write(*,"(9f10.3)") Matr
           do j=1,Am%natoms
             if(label == Am%Atom(j)%lab) then
                Am%Atom(j)%SkR=matmul(Matr,Am%Atom(j)%SkR)
                Am%Atom(j)%AtmInfo=constr_strings(i)
                exit
             end if
           end do
           !The treatment of the codes will be done in the future
         end do
       end if

       !Get pointers to the magnetic form factors
       !Stored for each atom in the component ind(1)
       call Set_Magnetic_Form()

       !---- Find Species in Magnetic_Form ----!
       do i=1,Am%natoms
          symbcar=u_case(Am%atom(i)%SfacSymb)
          do j=1,num_mag_form
             if (symbcar /= Magnetic_Form(j)%Symb) cycle
             Am%atom(i)%ind(1)=j
             exit
          end do
       end do

       return
    End Subroutine Readn_Set_Magnetic_Structure_MCIF

    Subroutine Get_mOrbit(x,Spg,Mult,orb,ptr)
       !---- Arguments ----!
       real(kind=cp), dimension(3),    intent (in) :: x
       type(Magnetic_Space_Group_type),intent (in) :: spg
       integer,                        intent(out) :: mult
       real(kind=cp),dimension(:,:),   intent(out) :: orb
       integer,dimension(:),optional,  intent(out) :: ptr

       !---- Local variables ----!
       integer                                :: j, nt
       real(kind=cp), dimension(3)            :: xx,v
       character(len=1)                       :: laty

       laty="P"
       mult=1
       orb(:,1)=x(:)
       if(present(ptr)) ptr(mult) = 1
       ext: do j=2,Spg%n_sym
          xx=ApplySO(Spg%SymOp(j),x)
          xx=modulo_lat(xx)
          do nt=1,mult
             v=orb(:,nt)-xx(:)
             if (Lattice_trans(v,laty)) cycle ext
          end do
          mult=mult+1
          orb(:,mult)=xx(:)
          if(present(ptr)) ptr(mult) = j   !Effective symop
       end do ext
       return
    End Subroutine Get_mOrbit
    !!----
    !!---- Subroutine Set_Shubnikov_Group(shubk,SG,MGp)
    !!----    character (len=*),         intent (in)    :: Shubk
    !!----    type(Magnetic_Group_Type), intent (out)   :: SG
    !!----    type(MagSymm_k_Type),      intent (in out):: MGp
    !!----
    !!----  This subroutined is not completed ... it is still in development
    !!---- Update: April 2008
    !!
    Subroutine Set_Shubnikov_Group(shubk,SG,MGp)
       !---- Arguments ----!
       character (len=*),         intent (in)    :: Shubk
       type(Magnetic_Group_Type), intent (out)   :: SG
       type(MagSymm_k_Type),      intent (in out):: MGp

       !---- Local Variables ----!
       !character (len=132) :: line
       character (len=20)  :: symb
       character (len=4)   :: gn
       character (len=4),dimension(10) :: gen
       logical,          dimension(10) :: found
       integer :: i,j, ng, k,m,n  !,nbl
       integer,              dimension(3)   :: bl
       integer,              dimension(10)  :: syp, numop
       !integer, allocatable, dimension(:,:) :: tab
       !character(len=*),parameter, dimension(26) :: oper = &
       !(/"1 ","-1","m ","2 ","21","3 ","31","32","-3","4 ","41","42","43",&
       !  "-4","6 ","61","62","63","64","65","-6","a ","b ","c ","d ","n "/)
       character(len=40),allocatable, dimension(:) :: ope


       SG%Shubnikov=" "
       SG%Shubnikov=adjustl(Shubk)
       gen = " "

       ! Generate the space group
       j=0
       bl=len_trim(SG%Shubnikov)
       !numop=0
       do i=1,len_trim(SG%Shubnikov)
          if (SG%Shubnikov(i:i) == " ") then
             j=j+1
             bl(j)=i
          end if
       end do

       SG%Shubnikov(bl(1):) = l_case( Shubk(bl(1):))   !Ensures lower case for symmetry elements

       !nbl=j
       j=0
       ng=0
       symb=" "
       syp=0
       do i=1,len_trim(SG%Shubnikov)
          j=j+1
          symb(j:j)=SG%Shubnikov(i:i)
          if (symb(j:j) == "'") then
             ng=ng+1
             k=5
             gn=" "
             do m=j-1,1,-1
                if (symb(m:m) == " ") exit
                if (symb(m:m) == "/") exit
                k=k-1
                gn(k:k)= symb(m:m)
             end do
             gen(ng)=adjustl(gn)
             if (i > bl(1)) syp(ng) = 1
             if (i > bl(2)) syp(ng) = 2
             if (i > bl(3)) syp(ng) = 3
             symb(j:j)=" "
             j= j-1
          end if
       end do
       i=index(symb," ")
       if ( i > 2) then
          symb=symb(1:1)//symb(i:)
       end if
       !write(*,*) " Space group symbol: ", trim(symb)
       !write(*,*) "  Primed Generators: ", (gen(i),i=1,ng), " in positions: ",(syp(i),i=1,ng)

       call Set_SpaceGroup(symb, SG%SpG)

                  !Determine the vector tinv from the information given for the generators
       SG%tinv=1  !by default the magnetic group is identical to the crystallographic group

       m=SG%SpG%Multip
       if (allocated(ope)) deallocate(ope)
       allocate(ope(m))

       found=.false.
       do j=1,ng
          if (gen(j) == "-1") then
             SG%tinv(SG%SpG%Numops+1) = -1
               found(SG%SpG%Numops+1) =.true.
             numop(j)= SG%SpG%Numops+1
          end if
       end do

       !         "Triclinic   ","Monoclinic  ","Orthorhombic","Tetragonal  ",    &
       !  "Trigonal","Hexagonal   ","Cubic       " /)
       n=1
       gn=SG%SpG%CrystalSys(1:4)
       do i=2,m  !over all symmetry operators of Space Group
          if(n == 0) exit  !all operators have been found
          call Symmetry_Symbol(SG%SpG%SymopSymb(i),ope(i))
          n=0
          do j=1,ng
             if (found(j)) cycle
             n=n+1
             if (gen(j)(1:1) == "-") then           !Search for roto-inversion axes
                k=index(ope(i),gen(j)(1:2)//"+")
                if (k /= 0) then    !Operator found
                   if (gn == "Cubi") then
                      k=index(ope(i),"x,x,x")
                      if (k /= 0) then
                         found(j)=.true.
                         SG%tinv(i)=-1
                         numop(j)= i
                         exit
                      end if
                   else
                      found(j)=.true.
                      SG%tinv(i)=-1
                      numop(j)= i
                      exit
                   end if
                end if
             else

                Select Case (gn)
                   case("Mono")
                      k=index(ope(i),gen(j)(1:1))    !Valid for all operators
                      if (k /= 0) then    !Operator found
                         found(j)=.true.
                         SG%tinv(i)=-1
                         numop(j)= i
                         exit
                      end if

                   case("Orth")
                      Select Case (gen(j))
                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)
                                     k=index(ope(i),"x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)
                                     k=index(ope(i),"y")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)
                                     k=index(ope(i),"z")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","a ","b ","c ","d ","n ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)
                                     k=index(ope(i),"x")
                                     if (k == 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)
                                     k=index(ope(i),"y")
                                     if (k == 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)
                                     k=index(ope(i),"z")
                                     if (k == 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                      End Select

                   case("Tetr")
                      Select Case (gen(j))
                         Case("4 ","41","42","43")           ! Look for 4-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")
                            if (k /= 0) then
                               found(j)=.true.
                               SG%tinv(i)=-1
                               numop(j)= i
                               exit
                            end if

                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(2)                ! along [100]
                                     k=index(ope(i),"x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! along [1-10]
                                     k=index(ope(i),"-x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","a ","b ","c ","d ","n ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)               ! Plane perp. to z (x,y,..) plane
                                     k=index(ope(i),"z")
                                     if (k == 0) then    !z should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)                ! perp. to [100]
                                     k=index(ope(i),"x")
                                     if (k == 0) then    !x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! perp. to [1-100]
                                     k=index(ope(i),"-x")
                                     if (k == 0) then    !-x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                      End Select

                   case("Rhom")
                      Select Case (gen(j))
                         Case("3 ","31","32")           ! Look for 3-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")
                            if (k /= 0) then
                               found(j)=.true.
                               SG%tinv(i)=-1
                               numop(j)= i
                               exit
                            end if

                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(2)                ! along [100]
                                     k=index(ope(i),"x,0")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! along [1-10]
                                     k=index(ope(i),"-x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","c ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)               ! Plane perp. to z (x,y,..) plane
                                     k=index(ope(i),"z")
                                     if (k == 0) then    !z should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)                ! perp. to [100]
                                     k=index(ope(i),"x")
                                     if (k == 0) then    !x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! perp. to [1-100]
                                     k=index(ope(i),"-x")
                                     if (k == 0) then    !-x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                      End Select

                   case("Hexa")
                      Select Case (gen(j))
                         Case("6 ","61","62","63","64","65")    ! Look for 6-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")    !only along z
                            if (k /= 0) then
                               found(j)=.true.
                               SG%tinv(i)=-1
                               numop(j)= i
                               exit
                            end if

                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(2)                ! along [100]
                                     k=index(ope(i),"x,0")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! along [1-10]
                                     k=index(ope(i),"-x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","c ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)               ! Plane perp. to z (x,y,..) plane
                                     k=index(ope(i),"z")
                                     if (k == 0) then    !z should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)                ! perp. to [100]
                                     k=index(ope(i),"x")
                                     if (k == 0) then    !x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! perp. to [1-100]
                                     k=index(ope(i),"-x")
                                     if (k == 0) then    !-x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                      End Select

                   case("Cubi")
                      Select Case (gen(j))
                         Case("4 ","41","42","43")    ! Look for 4-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")    !only along z
                            if (k /= 0) then
                               k=index(ope(i),"z")
                               if (k /= 0) then    !Operator found
                                  found(j)=.true.
                                  SG%tinv(i)=-1
                                  numop(j)= i
                                  exit
                               end if
                            end if

                         Case("3 ")    ! Look for 3-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")    !only along [111]
                            if (k /= 0) then
                               k=index(ope(i),"x,x,x")
                               if (k /= 0 ) then    !Operator found
                                  found(j)=.true.
                                  SG%tinv(i)=-1
                                  numop(j)= i
                                  exit
                               end if
                            end if

                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)                ! along [001]
                                     k=index(ope(i),"z")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! along [110]
                                     k=index(ope(i),"-x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","a ","b ","c ","d ","n ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)               ! Plane perp. to z (x,y,..) plane
                                     k=index(ope(i),"x,y,0")
                                     if (k /= 0) then
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! perp. to [1-10]
                                     k=index(ope(i),"x,x,z")
                                     if (k /= 0) then
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                     End Select

                End Select
             end if

          end do  !j=1,ng over all primed symmetry
       end do    !i=2,m over all symmetry operators of Space Group

       !write(*,*) "  Primed Generators: ", (gen(i),i=1,ng), " Correspond to operators: ",(numop(i),i=1,ng)

       !if(allocated(tab)) deallocate(tab)
       !allocate(tab(m,m))
       !call  Set_SpG_Mult_Table(SG%SpG,tab,.true.)

       !Construct MGp from the Shubnikov group
       !Just a dummy construction of MGp for avoiding warning or missbehaviour of compilers (provisional)
       call Init_MagSymm_k_Type(MGp)
       return
    End Subroutine Set_Shubnikov_Group

    !!----
    !!---- Subroutine Write_Magnetic_Structure(Ipr,MGp,Am,Mag_Dom)
    !!----    Integer,                    intent(in)           :: Ipr
    !!----    type(MagSymm_k_Type),       intent(in)           :: MGp
    !!----    type(mAtom_List_Type),      intent(in)           :: Am
    !!----    type(Magnetic_Domain_Type), intent(in), optional :: Mag_Dom
    !!----
    !!----    Subroutine to write out the information about the magnetic symmetry
    !!----    and mangnetic structure in unit Ipr.
    !!----
    !!---- Update: November 2006
    !!
    Subroutine Write_Magnetic_Structure(Ipr,MGp,Am,Mag_Dom)
       !---- Arguments ----!
       Integer,                    intent(in)           :: Ipr
       type(MagSymm_k_Type),       intent(in)           :: MGp
       type(mAtom_List_Type),      intent(in)           :: Am
       type(Magnetic_Domain_Type), intent(in), optional :: Mag_Dom

       !---- Local Variables ----!
       character (len=100), dimension( 4):: texto
       character (len=40)                :: aux
       integer :: i,j,k,l, nlines,n,m
       real(kind=cp)                  :: x
       complex                        :: ci
       real(kind=cp), dimension(3)    :: xp,xo
       complex, dimension(3)          :: Sk


       Write(unit=ipr,fmt="(/,a)")  "==================================="
       Write(unit=ipr,fmt="(  a)")  "== Magnetic Symmetry Information =="
       Write(unit=ipr,fmt="(a,/)")  "==================================="

       write(unit=ipr,fmt="(a)")    " => Magnetic  model name: "//trim(MGp%MagModel)
       write(unit=ipr,fmt="(a)")    " => Crystal lattice type: "//MGp%Latt
       if (MGp%nirreps == 0) then
          write(unit=ipr,fmt="(a,i2)") " => Number of Magnetic operators/Crystallographic operator: ",MGp%nmsym
       else
          write(unit=ipr,fmt="(a,i2)") " => Number of Irreducible Representations: ",MGp%nirreps
          do i=1,MGp%nirreps
             write(unit=ipr,fmt="(2(a,i3),a,12i2)") " => Number of basis functions of Irreducible Representation #",i," :", &
                                 MGp%nbas(i),"  Indicators for real(0)/imaginary(1): ", MGp%icomp(1:abs(MGp%nbas(i)),i)
          end do
       end if

       If(Am%Natoms > 0) then
         If (MGp%Centred == 2) then
            write(unit=ipr,fmt="(a)")    " => The crystallographic structure is centric (-1 at origin) "
         else
            write(unit=ipr,fmt="(a)")    " => The crystallographic structure is acentric  "
         End if
         if (MGp%MCentred == 2) then
            write(unit=ipr,fmt="(a)")    " => The magnetic structure is centric "
         else
            if (MGp%Centred == 2) then
               write(unit=ipr,fmt="(a)")    " => The magnetic structure is anti-centric  "
            else
               write(unit=ipr,fmt="(a)")    " => The magnetic structure is acentric  "
            end if
         End if
       End if
       write(unit=ipr,fmt="(a,i2)") " => Number of propagation vectors: ",MGp%nkv
       do i=1,MGp%nkv
          write(unit=ipr,fmt="(a,i2,a,3f8.4,a)") " => Propagation vectors #",i," = (",MGp%kvec(:,i)," )"
       end do
       if (MGp%Numlat > 1) then
          write(unit=ipr,fmt="(a,i3)")  " => Centring vectors:",MGp%Numlat-1
          nlines=1
          texto(:) (1:100) = " "
          do i=2,MGp%Numlat
             call Frac_Trans_2Dig(MGp%Ltr(:,i),aux)
             if (mod(i-1,2) == 0) then
                write(unit=texto(nlines)(51:100),fmt="(a,i2,a,a)") " => Latt(",i-1,"): ",trim(aux)
                nlines=nlines+1
             else
                write(unit=texto(nlines)( 1:50),fmt="(a,i2,a,a)") " => Latt(",i-1,"): ",trim(aux)
             end if
          end do
          do i=1,nlines
             write(unit=ipr,fmt="(a)") texto(i)
          end do
       end if

       If(MGp%Numops > 0) then
         write(unit=ipr,fmt="(/,a,/)")        " => List of all Symmetry Operators and Symmetry Symbols"

         do i=1,MGp%Numops
            texto(1)=" "
            call Symmetry_Symbol(MGp%SymopSymb(i),texto(1))
            write(unit=ipr,fmt="(a,i3,2a,t50,2a)") " => SYMM(",i,"): ",trim(MGp%SymopSymb(i)), &
                                                            "Symbol: ",trim(texto(1))
            if (MGp%nirreps == 0) then
              do j=1,MGp%NMSym
                 write(unit=ipr,fmt="(a,2(i2,a))")      "    MSYMM(",i,",",j,"): "//trim(MGp%MSymopSymb(i,j))
              end do
            else
              do j=1,MGp%nirreps
                write(unit=ipr,fmt="(a,2(i2,a),12(3f9.4,tr2))")"    BASR(",i,",",j,"): ",real(MGp%Basf(:,1:abs(MGp%nbas(j)),i,j))
                if (MGp%nbas(j) < 0) &
                write(unit=ipr,fmt="(a,2(i2,a),12(3f9.4,tr2))")"    BASI(",i,",",j,"): ",AIMAG(MGp%Basf(:,1:abs(MGp%nbas(j)),i,j))
              end do
            end if
         end do
       End if  !MGp%Numops > 0

       If(Am%Natoms > 0) then
         Write(unit=ipr,fmt="(/,a)")  "===================================="
         Write(unit=ipr,fmt="(  a)")  "== Magnetic Structure Information =="
         Write(unit=ipr,fmt="(a,/)")  "===================================="

         Write(unit=ipr,fmt="(a)")    " "
         Write(unit=ipr,fmt="(  a)")  "== Magnetic Asymmetric Unit Data =="
         Write(unit=ipr,fmt="(a,/)")  " "

         if (MGp%nirreps == 0) then
            Write(unit=ipr,fmt="(a)")  &
            "  The Fourier coefficients are of the form: Sk(j) = 1/2 { Rk(j) + i Ik(j) } exp {-2pi i Mphask(j)}"
            Write(unit=ipr,fmt="(a)")  &
            "  They are written for each atom j as Sk( j)= 1/2 {(Rx Ry Rz) + i ( Ix Iy Iz)} exp {-2pi i Mphask} -> MagMatrix # imat"
            Write(unit=ipr,fmt="(a)")  "  In case of k=2H (H reciprocal lattice vector) Sk(j)= (Rx Ry Rz)"

            do i=1,Am%Natoms
               Write(unit=ipr,fmt="(a,a,5f10.5)")  &
                 "   Atom "//Am%Atom(i)%Lab, Am%Atom(i)%SfacSymb, Am%Atom(i)%x,Am%Atom(i)%Biso,Am%Atom(i)%occ
               do j=1,Am%Atom(i)%nvk
                  if (K_Equiv_Minus_K(MGp%kvec(:,j),MGp%latt)) then
                     Write(unit=ipr,fmt="(a,i2,a,3f10.5,a,i4)")  &
                     "     Sk(",j,") =  (", Am%Atom(i)%Skr(:,j),")  -> MagMatrix #", Am%Atom(i)%imat(j)
                  else
                     Write(unit=ipr,fmt="(a,i2,a,2(3f10.5,a),f9.5,a,i4)")  &
                     "     Sk(",j,") = 1/2 {(", Am%Atom(i)%Skr(:,j),") + i (",Am%Atom(i)%Ski(:,j),")}  exp { -2pi i ",&
                     Am%Atom(i)%MPhas(j),"}  -> MagMatrix #", Am%Atom(i)%imat(j)
                  end if
               end do
            end do

         else

            Write(unit=ipr,fmt="(a)")  &
            "  The Fourier coefficients are of the form: Sk(j) = 1/2 Sum(i){Ci* Basf(i,imat)} exp {-2pi i Mphask(j)}"
            Write(unit=ipr,fmt="(a)")  &
            "  Where Ci are coefficients given below, Basf are the basis functions given above -> Irrep# imat"

            do i=1,Am%Natoms
               Write(unit=ipr,fmt="(a,a,5f10.5)")  &
                 "   Atom "//Am%Atom(i)%Lab, Am%Atom(i)%SfacSymb, Am%Atom(i)%x,Am%Atom(i)%Biso,Am%Atom(i)%occ
               do j=1,Am%Atom(i)%nvk
                  m=Am%Atom(i)%imat(j)
                  n=abs(MGp%nbas(m))
                  !1234567890123456789012345678
                  aux="(a,i2,a,  f10.5,a,f9.5,a,i4)"
                  write(unit=aux(9:10),fmt="(i2)") n
                  Write(unit=ipr,fmt=aux)  &
                     "  Coef_BasF(",j,") = 1/2 {(", Am%Atom(i)%cbas(1:n,j),")}  exp { -2pi i ",&
                  Am%Atom(i)%MPhas(j),"}  -> Irrep #", m
               end do
            end do
         end if

         ! Complete list of all atoms per primitive cell
         Write(unit=ipr,fmt="(/,a)")  " "
         Write(unit=ipr,fmt="(  a)")  "== List of all atoms and Fourier coefficients in the primitive cell =="
         Write(unit=ipr,fmt="(a,/)")  " "

         ! Construct the Fourier coefficients in case of Irreps
         if (MGp%nirreps /= 0 ) then
            do i=1,Am%natoms
               xo=Am%Atom(i)%x
               do k=1,MGp%NumOps
                  xp=ApplySO(MGp%SymOp(k),xo)
                  Write(unit=ipr,fmt="(a,i2,a,3f8.5)") " =>  Atom "//Am%Atom(i)%lab//"(",k,") :",xp
                  do j=1,Am%Atom(i)%nvk
                     m=Am%Atom(i)%imat(j)
                     n=abs(MGp%nbas(m))
                     Sk(:) = cmplx(0.0,0.0)
                     do l=1,n !cannot be greater than 12 at present
                        x=real(MGp%icomp(l,m))
                        ci=cmplx(1.0-x,x)
                        Sk(:)=Sk(:)+ Am%atom(i)%cbas(l,m)*ci* MGp%basf(:,l,k,m)
                     end do
                     x=-tpi*Am%atom(i)%mphas(j)
                     Sk=Sk*cmplx(cos(x),sin(x))
                     Write(unit=ipr,fmt="(a,i2,a,2(3f10.5,a),f9.5,a)")  &
                      "     Sk(",j,") = 1/2 {(", real(Sk),") + i (",aimag(Sk),")}"
                  end do
               end do  !Ops
               Write(unit=ipr,fmt="(a)") "  "
            end do  !atoms

         else
            do i=1,Am%natoms
               xo=Am%Atom(i)%x
               do k=1,MGp%NumOps
                  xp=ApplySO(MGp%SymOp(k),xo)
                  Write(unit=ipr,fmt="(a,i2,a,3f8.5)") " =>  Atom "//Am%Atom(i)%lab//"(",k,") :",xp
                  do j=1,Am%Atom(i)%nvk
                     m=Am%Atom(i)%imat(j)
                     n=abs(MGp%nbas(m))
                     x=-tpi*Am%atom(i)%mphas(j)
                     Sk=cmplx(Am%Atom(i)%Skr(:,j),Am%Atom(i)%Ski(:,j))
                     Sk= ApplyMSO(MGp%MSymOp(k,j),Sk)*cmplx(cos(x),sin(x))
                     Write(unit=ipr,fmt="(a,i2,a,2(3f10.5,a),f9.5,a)")  &
                      "     Sk(",j,") = 1/2 {(", real(Sk),") + i (",aimag(Sk),")}"
                  end do
               end do  !Ops
               Write(unit=ipr,fmt="(a)") "  "
            end do  !atoms
         end if

       End If !Am%Natoms > 0

       ! Writing information about domains (like in FullProf)
       if (present(Mag_Dom)) then
          write(unit=ipr,fmt="(a)") " => Magnetic S-Domains are present"
          if(Mag_Dom%chir) write(unit=ipr,fmt="(a)")"    Chirality domains are also present                     Chir-1      Chir-2"
          do i=1,Mag_Dom%nd
             if (Mag_Dom%chir) then
                write(unit=ipr,fmt="(a,i2,1(a,2f12.4))")"      Matrix of Magnetic Domain #:",i, &
                   " -> Populations: ",Mag_Dom%Pop(1:2,i) !,'  Codes:',MagDom(iom)%MPop(1:2,i)
             else
                write(unit=ipr,fmt="(a,i2,1(a,f12.4))")"      Matrix of Magnetic Domain #:",i,  &
                   " -> Population: ",Mag_Dom%Pop(1,i) !,'  Code:',MagDom(iom)%MPop(1,i)
             end if
             do j=1,3
                write(unit=ipr,fmt="(a,3i4)")  "                    ",Mag_Dom%Dmat(j,:,i)
            end do
          end do
       end if

       return
    End Subroutine Write_Magnetic_Structure

    Subroutine Write_MCIF(Ipr,mCell,MSGp,Am,Cell)
       Integer,                         intent(in)           :: Ipr
       type(Magnetic_Space_Group_Type), intent(in)           :: MSGp
       type(Crystal_Cell_Type),         intent(in)           :: mCell
       type(mAtom_List_Type),           intent(in)           :: Am
       type(Crystal_Cell_Type),optional,intent(in)           :: Cell
       !
       Character(len=132)             :: line
       character(len=40),dimension(6) :: text
       character(len=2)               :: invc
       real                           :: occ,occ_std,uiso,uiso_std
       integer :: i,j,k

       write(unit=Ipr,fmt="(a)") "#  --------------------------------------"
       write(unit=Ipr,fmt="(a)") "#  Magnetic CIF file generated by CrysFML"
       write(unit=Ipr,fmt="(a)") "#  --------------------------------------"
       write(unit=Ipr,fmt="(a)") "# https://forge.epn-campus.eu/projects/crysfml/repository"
       call Write_Date_Time(dtim=line)
       write(unit=Ipr,fmt="(a)") trim(line)
       write(unit=Ipr,fmt="(a)") " "

       write(unit=Ipr,fmt="(a)") "data_"
       write(unit=Ipr,fmt="(a)") "_citation_journal_abbrev ?"
       write(unit=Ipr,fmt="(a)") "_citation_journal_volume ?"
       write(unit=Ipr,fmt="(a)") "_citation_page_first     ?"
       write(unit=Ipr,fmt="(a)") "_citation_page_last      ?"
       write(unit=Ipr,fmt="(a)") "_citation_article_id     ?"
       write(unit=Ipr,fmt="(a)") "_citation_year           ?"
       write(unit=Ipr,fmt="(a)") "_loop "
       write(unit=Ipr,fmt="(a)") "_citation_author_name"
       write(unit=Ipr,fmt="(a)") "?"
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "_atomic_positions_source_database_code_ICSD  ?"
       write(unit=Ipr,fmt="(a)") "_atomic_positions_source_other    .  "
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "_Neel_temperature  ?"
       write(unit=Ipr,fmt="(a)") "_magn_diffrn_temperature  ?"
       write(unit=Ipr,fmt="(a)") "_exptl_crystal_magnetic_properties_details"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") "_active_magnetic_irreps_details"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") " "
       if(MSGp%standard_setting) then
          write(unit=Ipr,fmt="(a)") "_magnetic_space_group_standard_setting  'yes'"
       else
          write(unit=Ipr,fmt="(a)") "_magnetic_space_group_standard_setting  'no'"
       end if
       write(unit=Ipr,fmt="(a)")    '_parent_space_group.name_H-M  "'//trim(MSGp%Parent_spg)//'"'
       write(unit=Ipr,fmt="(a,i3)") "_parent_space_group.IT_number  ",MSGp%Parent_num
       write(unit=Ipr,fmt="(a)")    "_magnetic_space_group.transform_from_parent_Pp_abc  '"//trim(MSGp%trn_from_parent)//"'"
       write(unit=Ipr,fmt="(a)")    "_magnetic_space_group.transform_to_standard_Pp_abc  '"//trim(MSGp%trn_to_standard)//"'"
       write(unit=Ipr,fmt="(a)")
       if(len_trim(MSGp%BNS_number) /= 0) &
       write(unit=Ipr,fmt="(a)") "_magnetic_space_group_BNS_number  "//trim(MSGp%BNS_number)
       if(len_trim(MSGp%BNS_symbol) /= 0) &
       write(unit=Ipr,fmt="(a)") '_magnetic_space_group_BNS_name  "'//trim(MSGp%BNS_symbol)//'"'
       if(len_trim(MSGp%OG_number) /= 0) &
       write(unit=Ipr,fmt="(a)") '_magnetic_space_group_OG_number '//trim(MSGp%OG_number)
       if(len_trim(MSGp%OG_symbol) /= 0) &
       write(unit=Ipr,fmt="(a)") '_magnetic_space_group_OG_name  "'//trim(MSGp%OG_symbol)//'"'
       write(unit=Ipr,fmt="(a)")

       if(MSGp%n_irreps /= 0) then
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_irrep_id"
          write(unit=Ipr,fmt="(a)") "_irrep_dimension"
          if( any(MSGp%small_irrep_dim > 0) ) write(unit=Ipr,fmt="(a)") "_small_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_irrep_direction_type"
          write(unit=Ipr,fmt="(a)") "_irrep_action"
          if( any(MSGp%irrep_modes_number > 0) ) write(unit=Ipr,fmt="(a)") "_irrep_modes_number"
          do i=1,MSGp%n_irreps
            if(MSGp%small_irrep_dim(i) > 0) then
               write(unit=line,fmt=("(2i4)"))  MSGp%irrep_dim(i), MSGp%small_irrep_dim(i)
            else
               write(unit=line,fmt=("(i4)"))  MSGp%irrep_dim(i)
            end if
            line= trim(MSGp%irrep_id(i))//"  "//trim(line)//"   "// &
                                      trim(MSGp%irrep_direction(i))//"  "//trim(MSGp%irrep_action(i))
            if( MSGp%irrep_modes_number(i) > 0) then
               j=len_trim(line)
              write(unit=line(j+1:),fmt="(i4)") MSGp%irrep_modes_number(i)
            end if
            write(unit=Ipr,fmt="(a)") trim(line)
          end do
          write(unit=Ipr,fmt="(a)")
       else
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_irrep_id"
          write(unit=Ipr,fmt="(a)") "_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_small_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_irrep_direction_type"
          write(unit=Ipr,fmt="(a)") "_irrep_action"
          write(unit=Ipr,fmt="(a)") "_irrep_modes_number"
          write(unit=Ipr,fmt="(a)") " ?  ?  ?  ?  ?  ?"
          write(unit=Ipr,fmt="(a)")
       end if

       if(MSGp%m_cell) then
          do i=1,3
            call setnum_std(mCell%Cell(i),mCell%cell_std(i),text(i))
            call setnum_std(mCell%ang(i),mCell%ang_std(i),text(i+3))
          end do
          write(unit=Ipr,fmt="(a)") "_magnetic_cell_length_a    "//trim(text(1))
          write(unit=Ipr,fmt="(a)") "_magnetic_cell_length_b    "//trim(text(2))
          write(unit=Ipr,fmt="(a)") "_magnetic_cell_length_c    "//trim(text(3))
          write(unit=Ipr,fmt="(a)") "_magnetic_cell_angle_alpha "//trim(text(4))
          write(unit=Ipr,fmt="(a)") "_magnetic_cell_angle_beta  "//trim(text(5))
          write(unit=Ipr,fmt="(a)") "_magnetic_cell_angle_gamma "//trim(text(6))
          write(unit=Ipr,fmt="(a)")
       else
          if(present(Cell)) then
             do i=1,3
               call setnum_std(Cell%Cell(i),Cell%cell_std(i),text(i))
               call setnum_std(Cell%ang(i),Cell%ang_std(i),text(i+3))
             end do
             write(unit=Ipr,fmt="(a)") "_cell_length_a    "//trim(text(1))
             write(unit=Ipr,fmt="(a)") "_cell_length_b    "//trim(text(2))
             write(unit=Ipr,fmt="(a)") "_cell_length_c    "//trim(text(3))
             write(unit=Ipr,fmt="(a)") "_cell_angle_alpha "//trim(text(4))
             write(unit=Ipr,fmt="(a)") "_cell_angle_beta  "//trim(text(5))
             write(unit=Ipr,fmt="(a)") "_cell_angle_gamma "//trim(text(6))
             write(unit=Ipr,fmt="(a)")
          end if
       end if
       if(MSGp%n_kv > 0) then
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_magnetic_propagation_vector_seq_id"
          write(unit=Ipr,fmt="(a)") "_magnetic_propagation_vector_kxkykz"
          do i=1,MSGp%n_kv
            call Frac_Trans_2Dig(MSGp%kv(:,i),line)
            line=line(2:len_trim(line)-1)
            write(unit=Ipr,fmt="(a)") trim(MSGp%kv_label(i))//"  '"//trim(line)//"'"
          end do
       end if
       if(MSGp%m_constr) then
          write(unit=Ipr,fmt="(a)")
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_moment_symmetry_constraints_label"
          write(unit=Ipr,fmt="(a)") "_atom_site_magnetic_moment_symmetry_constraints_mxmymz"
          do i=1,Am%natoms
            line=Am%Atom(i)%AtmInfo
            if(len_trim(line) < 8) cycle
            write(unit=Ipr,fmt="(a)")trim(line)
          end do
       end if
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)")  "loop_"
       write(unit=Ipr,fmt="(a)")  "_magnetic_space_group_symop_id"
       write(unit=Ipr,fmt="(a)")  "_magnetic_space_group_symop_operation_xyz"
       write(unit=Ipr,fmt="(a)")  "_magnetic_space_group_symop_operation_mxmymz"
       write(unit=Ipr,fmt="(a)")  "_magnetic_space_group_symop_operation_timereversal"
       do i=1,MSGp%n_sym
          write(unit=invc,fmt="(i2)") nint(MSgp%MSymop(i)%Phas)
          if(invc(1:1) == " ") invc(1:1)="+"
          write(unit=Ipr,fmt="(i3,a)") i," "//trim(MSgp%SymopSymb(i))//" "//trim(MSgp%MSymopSymb(i))//" "//invc
       end do
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "loop_"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_label"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_type_symbol"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_fract_x"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_fract_y"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_fract_z"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_U_iso_or_equiv"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_occupancy"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_symmetry_multiplicity"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_Wyckoff_label"
       line=" "
       do i=1,Am%natoms
          do j=1,3
            call setnum_std(Am%atom(i)%x(j),Am%atom(i)%x_std(j),text(j))
          end do
          occ=real(MSgp%n_sym)/real(Am%atom(i)%Mult)*Am%atom(i)%occ
          occ_std=real(MSgp%n_sym)/real(Am%atom(i)%Mult)*Am%atom(i)%occ_std
          call setnum_std(occ,occ_std,text(5))
          uiso=Am%atom(i)%biso/78.95683521
          uiso_std=Am%atom(i)%biso_std/78.95683521
          call setnum_std(uiso,uiso_std,text(4))
          write(unit=Ipr,fmt="(a6,a6,3a13,2a11,i4,a)") Am%Atom(i)%lab, Am%atom(i)%SfacSymb,(text(j),j=1,5),&
                                                       Am%atom(i)%Mult," "//Am%atom(i)%wyck
       end do
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "loop_"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_moment_label"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_moment_crystalaxis_mx"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_moment_crystalaxis_my"
       write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_moment_crystalaxis_mz"
       do i=1,Am%natoms
          if(sum(abs(Am%Atom(i)%Skr(:,1))) < 0.0001) cycle
          do j=1,3
            call setnum_std(Am%atom(i)%Skr(j,1),Am%atom(i)%Skr_std(j,1),text(j))
          end do
          write(unit=Ipr,fmt="(a8,3a12)") Am%Atom(i)%lab,(text(j),j=1,3)
       end do
       write(unit=Ipr,fmt="(a)")
       return
    End Subroutine Write_MCIF

    !!----
    !!---- Subroutine Write_Shubnikov_Group(SG,Iunit)
    !!----    type (Magnetic_Group_Type),intent(in) :: SG
    !!----    integer,   optional,       intent(in) :: iunit
    !!----
    !!----    Subroutine to write out the information about the Shubnikov_Group
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Write_Shubnikov_Group(SG,Iunit)
       !---- Arguments ----!
       type (Magnetic_Group_Type),intent(in) :: SG
       integer,   optional,       intent(in) :: iunit

       !---- Local variables ----!
       character (len=100), dimension(24):: texto
       character (len=40)                :: aux
       integer                           :: lun
       integer                           :: i, nlines
       logical                           :: print_latt

       !---- Initializing variables ----!
       lun=6
       if (present(iunit)) lun=iunit
       print_latt=.true.

       !---- Printing ----!
       write(unit=lun,fmt="(/,/,a)")          "        Information on Space Group: "
       write(unit=lun,fmt="(a,/ )")           "        --------------------------- "
       write(unit=lun,fmt="(a,a )")          " =>       Shubnikov Symbol: ", SG%Shubnikov
       write(unit=lun,fmt="(a,i3)")          " =>  Number of Space group: ", SG%SpG%NumSpg
       write(unit=lun,fmt="(a,a)")           " => Hermann-Mauguin Symbol: ", SG%SpG%SPG_Symb
       write(unit=lun,fmt="(a,a)")           " =>            Hall Symbol: ", SG%SpG%Hall
       write(unit=lun,fmt="(a,a)")           " =>   Table Setting Choice: ", SG%SpG%info
       write(unit=lun,fmt="(a,a)")           " =>           Setting Type: ", SG%SpG%SG_setting
       write(unit=lun,fmt="(a,a)")           " =>         Crystal System: ", SG%SpG%CrystalSys
       write(unit=lun,fmt="(a,a)")           " =>             Laue Class: ", SG%SpG%Laue
       write(unit=lun,fmt="(a,a)")           " =>            Point Group: ", SG%SpG%Pg
       write(unit=lun,fmt="(a,a)")           " =>        Bravais Lattice: ", SG%SpG%SPG_Lat
       write(unit=lun,fmt="(a,a)")           " =>         Lattice Symbol: ", SG%SpG%SPG_Latsy
       write(unit=lun,fmt="(a,i3)")          " => Reduced Number of S.O.: ", SG%SpG%NumOps
       write(unit=lun,fmt="(a,i3)")          " =>   General multiplicity: ", SG%SpG%Multip
       write(unit=lun,fmt="(a,a)")           " =>         Centrosymmetry: ", SG%SpG%Centre
       write(unit=lun,fmt="(a,i3)")          " => Generators (exc. -1&L): ", SG%SpG%num_gen
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") " =>        Asymmetric unit: ", SG%SpG%R_Asym_Unit(1,1), &
                                                                " <= x <= ", SG%SpG%R_Asym_Unit(1,2)
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") "                            ", SG%SpG%R_Asym_Unit(2,1), &
                                                                " <= y <= ", SG%SpG%R_Asym_Unit(2,2)
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") "                            ", SG%SpG%R_Asym_Unit(3,1), &
                                                                " <= z <= ", SG%SpG%R_Asym_Unit(3,2)

       if (SG%SpG%centred == 0) then
          call Frac_Trans_1Dig(SG%SpG%Centre_coord,texto(1))
          write(unit=lun,fmt="(a,a)")          " =>              Centre at: ", trim(texto(1))
       end if
       if (SG%SpG%SPG_Lat == "Z" .or. print_latt) then
          texto(:) (1:100) = " "
          if (SG%SpG%SPG_Lat == "Z") then
             write(unit=lun,fmt="(a,i3)")          " => Non-conventional Centring vectors:",SG%SpG%Numlat
          else
             write(unit=lun,fmt="(a,i3)")          " => Centring vectors:",SG%SpG%Numlat-1
          end if
          nlines=1
          do i=2,SG%SpG%Numlat
             call Frac_Trans_1Dig(SG%SpG%Latt_trans(:,i),aux)
             if (mod(i-1,2) == 0) then
                write(unit=texto(nlines)(51:100),fmt="(a,i2,a,a)") &
                                           " => Latt(",i-1,"): ",trim(aux)
                nlines=nlines+1
             else
                write(unit=texto(nlines)( 1:50),fmt="(a,i2,a,a)")  &
                                           " => Latt(",i-1,"): ",trim(aux)
             end if
          end do
          do i=1,nlines
             write(unit=lun,fmt="(a)") texto(i)
          end do
       end if

       !---- Symmetry Operators ----!
       write(unit=lun,fmt="(/,a,/)")        " => List of all Symmetry Operators and Symmetry Symbols"

       do i=1,SG%SpG%Multip
          texto(1)=" "
          call Symmetry_Symbol(SG%SpG%SymopSymb(i),texto(1))
          write(unit=lun,fmt="(a,i3,2a,t50,2a)") " => SYMM(",i,"): ",trim(SG%SpG%SymopSymb(i)), &
                                                    "Symbol: ",trim(texto(1))
       end do

       return
    End Subroutine Write_Shubnikov_Group

 End Module CFML_Magnetic_Symmetry

