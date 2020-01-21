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
!!----    Update: January 2020 (JRC: many subroutines have been removed from this module and moved to CFML_Crystallographic_Symmety
!!---                                or to CFML_IO_Formats, so no magnetic groups remains in this module except for obsolete Shubnikov group type)
!!----
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,                only: cp, dp, tpi
!!--++    Use CFML_Math_General,              only: Modulo_Lat
!!--++    Use CFML_Math_3D,                   only: Get_Cart_From_Spher, matrix_inverse, Veclength
!!--++    Use CFML_Symmetry_Tables,           only: ltr_a,ltr_b,ltr_c,ltr_i,ltr_r,ltr_f
!!--++    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Read_Xsym, Get_SymSymb, axes_rotation, &
!!--++                                              Sym_Oper_Type, Set_SpaceGroup,read_msymm, symmetry_symbol, &
!!--++                                              err_symm,err_symm_mess, set_SpG_Mult_Table,ApplySO,   &
!!--++                                              Lattice_Trans, Get_SO_from_Gener,Get_Centring_Vectors
!!--++    Use CFML_String_Utilities,          only: u_case, l_case, Frac_Trans_1Dig, Get_Separator_Pos,Pack_String, &
!!--++                                              Frac_Trans_2Dig, Get_Mat_From_Symb, getnum_std, Err_String,     &
!!--++                                              Err_String_Mess,setnum_std, getword
!!--++    Use CFML_IO_Formats,                only: file_list_type, File_To_FileList
!!--++    Use CFML_Atom_TypeDef,              only: Allocate_mAtom_list, mAtom_List_Type, Get_Atom_2nd_Tensor_Ctr
!!--++    Use CFML_Scattering_Chemical_Tables,only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
!!--++                                              Magnetic_Form
!!--++    Use CFML_Propagation_Vectors,       only: K_Equiv_Minus_K
!!--++    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Set_Crystal_Cell
!!--++    Use CFML_Magnetic_Groups
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
!!----       CALC_INDUCED_SK
!!----       INIT_ERR_MAGSYM
!!----       MAGNETIC_SPACE_GROUP_TYPE_TO_MAGSYMM_K_TYPE
!!----       MAGSYMM_K_TYPE_TO_MAGNETIC_SPACE_GROUP_TYPE
!!----       READN_SET_MAGNETIC_STRUCTURE
!!--++       READN_SET_MAGNETIC_STRUCTURE_CFL    [Overloaded]
!!--++       READN_SET_MAGNETIC_STRUCTURE_MCIF   [Overloaded]
!!----       SET_SHUBNIKOV_GROUP
!!----       WRITE_MAGNETIC_STRUCTURE
!!----       WRITE_MCIF
!!----       WRITE_SHUBNIKOV_GROUP
!!----
!!
 Module CFML_Magnetic_Symmetry

    !---- Use Modules ----!
    Use CFML_GlobalDeps,                only: cp, dp, tpi, Write_Date_Time
    Use CFML_Math_General,              only: Trace, Zbelong, Modulo_Lat, equal_matrix,             &
                                              Equal_Vector
    Use CFML_Math_3D,                   only: Get_Cart_From_Spher,Determ_A, matrix_inverse, Veclength
    Use CFML_Symmetry_Tables,           only: ltr_a,ltr_b,ltr_c,ltr_i,ltr_r,ltr_f,Sys_cry,LATT
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Read_Xsym, Get_SymSymb,axes_rotation, &
                                              Sym_Oper_Type, Set_SpaceGroup,read_msymm, symmetry_symbol, &
                                              err_symm,err_symm_mess, set_SpG_Mult_Table,ApplySO,   &
                                              Lattice_Trans, Get_SO_from_Gener, Get_Centring_Vectors, &
                                              Get_Shubnikov_Operator_Symbol, MSym_Oper_Type, LatSym,&
                                              Magnetic_Space_Group_Type, Get_Stabilizer,            &
                                              Init_Magnetic_Space_Group_Type, Get_mOrbit
    Use CFML_String_Utilities,          only: u_case, l_case, Frac_Trans_1Dig, Get_Separator_Pos,Pack_String, &
                                              Frac_Trans_2Dig, Get_Mat_From_Symb, getnum_std, Err_String,     &
                                              Err_String_Mess,setnum_std, getword, Get_Transf,ucase
    Use CFML_IO_Formats,                only: file_list_type, File_To_FileList, cleanup_symmetry_operators
    Use CFML_Atom_TypeDef,              only: Allocate_mAtom_list, mAtom_List_Type, Get_Atom_2nd_Tensor_Ctr
    Use CFML_Scattering_Chemical_Tables,only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
                                              Magnetic_Form, get_magnetic_form_factor
    Use CFML_Propagation_Vectors,       only: K_Equiv_Minus_K
    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Set_Crystal_Cell
    Use CFML_Magnetic_Groups

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public :: ApplyMSO

    !---- List of public subroutines ----!
    public :: Readn_Set_Magnetic_Structure, Write_Magnetic_Structure, Set_Shubnikov_Group, &
              Write_Shubnikov_Group, Init_MagSymm_k_Type, &
              Calc_Induced_Sk

    !---- Definitions ----!


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
    !!---- TYPE :: MAGNETIC_GROUP_TYPE
    !!--..
    !!---- Type, Public :: Magnetic_Group_Type
    !!----    Character(len=50)           :: Shubnikov !Shubnikov symbol (Hermman-Mauguin + primes)
    !!----    type(Space_Group_Type)      :: SpG       !Crystallographic space group
    !!----    integer, dimension(192)     :: tinv      !When a component is +1 no time inversion is associated
    !!---- End Type Magnetic_Group_Type                !If tinv(i)=-1, the time inversion is associated to operator "i"
    !!----
    !!--<<
    !!----    A magnetic group type is adequate when k=(0,0,0). It contains as the second
    !!----    component the crystallographic space group. The first component is
    !!----    the Shubnikov Group symbol and the third component is an integer vector with
    !!----    values -1 or 1 when time inversion is associated (-1) with the corresponding
    !!----    crystallographic symmetry operator o not (1). This is now OBSOLETE thanks to the
    !!----    availability of Shubnikov tables and the Magnetic_Space_Group_Type defined in
    !!----    Crystallographic Symmetry module.
    !!-->>
    !!----
    !!---- Update: April - 2005
    !!
    Type, Public :: Magnetic_Group_Type
       Character(len=50)           :: Shubnikov
       type(Space_Group_Type)      :: SpG
       integer, dimension(192)     :: tinv
    End Type Magnetic_Group_Type
    !!----
    !!---- TYPE :: MAGSYMM_K_TYPE
    !!--..
    !!---- Type, Public :: MagSymm_k_Type
    !!----    character(len=31)                        :: MagModel   ! Name to characterize the magnetic symmetry
    !!----    character(len=10)                        :: Sk_type    ! If Sk_type="Spherical_Frame" the input Fourier coefficients are in spherical components
    !!----    character(len=15)                        :: BNS_number ! Added for keeping the same information
    !!----    character(len=15)                        :: OG_number  ! as in Magnetic_Space_Group_Type
    !!----    Character(len=34)                        :: BNS_symbol !             "
    !!----    Character(len=34)                        :: OG_symbol  !             "
    !!----    Integer                                  :: MagType    !             "
    !!----    Integer                                  :: Parent_num !             "
    !!----    Character(len=20)                        :: Parent_spg !             "
    !!----    character(len=1)                         :: Latt       ! Symbol of the crystallographic lattice
    !!----    integer                                  :: nirreps    ! Number of irreducible representations (max=4, if nirreps /= 0 => nmsym=0)
    !!----    Integer,             dimension(4)        :: irrep_dim       !Dimension of the irreps
    !!----    Integer,             dimension(4)        :: small_irrep_dim !Dimension of the small irrep
    !!----    Integer,             dimension(4)        :: irrep_modes_number !Number of the mode of the irrep
    !!----    Character(len=15),   dimension(4)        :: irrep_id        !Labels for the irreps
    !!----    Character(len=20),   dimension(4)        :: irrep_direction !Irrep direction in representation space
    !!----    Character(len=20),   dimension(4)        :: irrep_action    !Irrep character primary or secondary
    !!----    integer                                  :: nmsym      ! Number of magnetic operators per crystallographic operator (max=8)
    !!----    integer                                  :: centred    ! =0 centric centre not at origin, =1 acentric, =2 centric (-1 at origin)
    !!----    integer                                  :: mcentred   ! =1 Anti/a-centric Magnetic symmetry, = 2 centric magnetic symmetry
    !!----    integer                                  :: nkv        ! Number of independent propagation vectors
    !!----    real(kind=cp),       dimension(3,12)     :: kvec       ! Propagation vectors
    !!----    Character(len=15),   dimension(12)       :: kv_label
    !!----    integer                                  :: Num_Lat    ! Number of centring lattice vectors
    !!----    real(kind=cp), dimension(3,4)            :: Ltr        ! Centring translations
    !!----    integer                                  :: Numops     ! Reduced number of crystallographic Symm. Op.
    !!----    integer                                  :: Multip     ! General multiplicity of the space group
    !!----    integer,             dimension(4)        :: nbas       ! Number of basis functions per irrep (if nbas < 0, the corresponding basis is complex).
    !!----    integer,             dimension(12,4)     :: icomp      ! Indicator (0 pure real/ 1 pure imaginary) for coefficients of basis fucntions
    !!----    Complex(kind=cp),    dimension(3,12,48,4):: basf       ! Basis functions of the irreps of Gk
    !!----    character(len=40),   dimension(:),   allocatable :: SymopSymb  ! Alphanumeric Symbols for SYMM
    !!----    type(Sym_Oper_Type), dimension(:),   allocatable :: SymOp      ! Crystallographic symmetry operators (48)
    !!----    character(len=40),   dimension(:,:), allocatable :: MSymopSymb ! Alphanumeric Symbols for MSYMM (48,8)
    !!----    type(MSym_Oper_Type),dimension(:,:), allocatable :: MSymOp     ! Magnetic symmetry operators (48,8)
    !!---- End Type MagSymm_k_Type
    !!----
    !!----  Definition of the MagSymm_k_type derived type, encapsulating the information
    !!----  concerning the crystallographic symmetry, propagation vectors and magnetic matrices.
    !!----  Needed for calculating magnetic structure factors.
    !!----
    !!---- Created: April   - 2005
    !!---- Updated: January - 2014
    !!
    Type, Public :: MagSymm_k_Type
       character(len=31)                        :: MagModel
       character(len=15)                        :: Sk_type
       character(len=15)                        :: BNS_number ! Added for keeping the same information
       character(len=15)                        :: OG_number  ! as in Magnetic_Space_Group_Type
       Character(len=34)                        :: BNS_symbol !             "
       Character(len=34)                        :: OG_symbol  !             "
       Integer                                  :: MagType    !             "
       Integer                                  :: Parent_num !             "
       Character(len=20)                        :: Parent_spg !             "
       character(len=1)                         :: Latt
       integer                                  :: nirreps
       Integer,             dimension(4)        :: irrep_dim          !Dimension of the irreps
       Integer,             dimension(4)        :: small_irrep_dim    !Dimension of the small irrep
       Integer,             dimension(4)        :: irrep_modes_number !Number of the mode of the irrep
       Character(len=15),   dimension(4)        :: irrep_id           !Labels for the irreps
       Character(len=20),   dimension(4)        :: irrep_direction    !Irrep direction in representation space
       Character(len=20),   dimension(4)        :: irrep_action       !Irrep character primary or secondary
       integer                                  :: nmsym
       integer                                  :: centred
       integer                                  :: mcentred
       integer                                  :: nkv
       real(kind=cp),dimension(3,12)            :: kvec
       integer                                  :: Num_Lat
       real(kind=cp), dimension(3,4)            :: Ltr
       integer                                  :: Numops
       integer                                  :: Multip
       integer,             dimension(4)        :: nbas
       integer,             dimension(12,4)     :: icomp
       Complex(kind=cp),    dimension(3,12,48,4):: basf
       character(len=40),   dimension(:),   allocatable :: SymopSymb  ! Alphanumeric Symbols for SYMM
       type(Sym_Oper_Type), dimension(:),   allocatable :: SymOp      ! Crystallographic symmetry operators (48)
       character(len=40),   dimension(:,:), allocatable :: MSymopSymb ! Alphanumeric Symbols for MSYMM (48,8)
       type(MSym_Oper_Type),dimension(:,:), allocatable :: MSymOp     ! Magnetic symmetry operators (48,8)
    End Type MagSymm_k_Type


    real(kind=cp), parameter, private:: eps_symm  = 0.0002_cp

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

    character (len=1), dimension(26),parameter, private   :: &
    cdd=(/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r', &
        's','t','u','v','w','x','y','z'/)
    real(kind=dp), parameter, private :: epps=0.000001_dp

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

    !!---- Subroutine Calc_Induced_Sk(cell,SpG,MField,dir_MField,Atm)
    !!----    !---- Arguments ----!
    !!----   type(Crystal_Cell_type),    intent(in)     :: Cell
    !!----   type(Space_Group_Type),     intent(in)     :: SpG
    !!----   real(kind=cp),              intent(in)     :: MField
    !!----   real(kind=cp),dimension(3), intent(in)     :: dir_MField
    !!----   type(Matom_list_type),    intent(in out)   :: Atm
    !!----
    !!----  This subroutine completes the object Am by calculating the
    !!----  induced magnetic moments of the representant atoms in the asymmetric unit.
    !!----  It modifies also the Chi tensor according to the symmetry constraints of
    !!----  the crystallographic site.
    !!----
    !!----  Created: June 2014 (JRC)
    !!----
    Subroutine Calc_Induced_Sk(cell,SpG,MField,dir_MField,Atm,ipr)
       !---- Arguments ----!
       type(Crystal_Cell_type),    intent(in)     :: Cell
       type(Space_Group_Type),     intent(in)     :: SpG
       real(kind=cp),              intent(in)     :: MField
       real(kind=cp),dimension(3), intent(in)     :: dir_MField
       type(Matom_list_type),    intent(in out)   :: Atm
       integer, optional,          intent(in)     :: ipr
       !--- Local variables ---!
       integer                          :: i,codini
       integer, dimension(6)            :: icodes
       real(kind=cp),    dimension(6)   :: multip
       real(kind=cp),    dimension(3)   :: u_vect,x
       real(kind=cp),    dimension(3,3) :: chi

       !
       u_vect=MField * dir_MField / Veclength(Cell%Cr_Orth_cel,dir_MField)
       if(present(ipr)) write(unit=ipr,fmt="(a,3f8.4)") " => Applied Magnetic Field: ",u_vect
       icodes=(/1,2,3,4,5,6/); multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
       codini=1
       do i=1,Atm%natoms
          x=atm%atom(i)%x
          call Get_Atom_2nd_Tensor_Ctr(x,atm%atom(i)%chi,SpG,Codini,Icodes,Multip)
          chi=reshape((/atm%atom(i)%chi(1),atm%atom(i)%chi(4), atm%atom(i)%chi(5), &
                        atm%atom(i)%chi(4),atm%atom(i)%chi(2), atm%atom(i)%chi(6), &
                        atm%atom(i)%chi(6),atm%atom(i)%chi(6), atm%atom(i)%chi(3) /),(/3,3/))
          Atm%atom(i)%SkR(:,1)=matmul(Chi,u_vect)
          Atm%atom(i)%SkI(:,1)=0.0
          if(present(ipr)) then
             write(unit=ipr,fmt="(a,i3,a,6f8.4)")     " Atom # ",i," Chi      values: ",atm%atom(i)%chi
             write(unit=ipr,fmt="(a,6i4,6f6.2)")      "            Chi constraints: ",Icodes,multip
             write(unit=ipr,fmt="(a,3f8.4)")          "            Induced  Moment: ",Atm%atom(i)%SkR(:,1)
          end if
       end do ! Atoms

       return
    End Subroutine Calc_Induced_Sk

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
    !!----  Update: April 2005, January 2014
    !!
    Subroutine Init_MagSymm_k_Type(MGp)
       !---- Arguments ----!
       type(MagSymm_k_Type),  intent (in out) :: MGp

       MGp%MagModel="Unnamed Model"
       MGp%Sk_Type="Crystal_Frame"       ! "Spherical_Frame"
       MGp%Latt="P"
       MGp%BNS_number=" "
       MGp%OG_number=" "
       MGp%BNS_symbol=" "
       MGp%OG_symbol=" "
       MGp%MagType=0
       MGp%Parent_num=0
       MGp%Parent_spg=" "
       MGp%nmsym=0
       MGp%nirreps=0
       MGp%irrep_dim=0          !Dimension of the irreps
       MGp%small_irrep_dim=0    !Dimension of the small irrep
       MGp%irrep_modes_number=0 !Number of the mode of the irrep
       MGp%irrep_id=" "         !Labels for the irreps
       MGp%irrep_direction=" "  !Irrep direction in representation space
       MGp%irrep_action=" "     !Irrep character primary or secondary
       MGp%centred=1    !By default the crystal structure is acentric
       MGp%mcentred=1   !By default the magnetic structure is anti-centric (if there is -1 it is combined with time inversion)
       MGp%nkv=0
       MGp%kvec=0.0
       MGp%Num_Lat=1
       MGp%Ltr=0.0
       MGp%Numops=0
       MGp%Multip=0
       MGp%nbas=0
       MGp%icomp=0
       MGp%basf=cmplx(0.0,0.0)
       return
    End Subroutine Init_MagSymm_k_Type



    Subroutine Magnetic_Space_Group_Type_to_MagSymm_k_Type(MSpG,mode,MG_Symk)
       Type(Magnetic_Space_Group_Type),   intent(in)  :: MSpG
       character(len=*),                  intent(in)  :: mode
       Type(MagSymm_k_Type),              intent(out) :: MG_Symk
       !---- Local variables ----!
       integer :: i,k
       logical :: full_convertion
       !real(kind=cp)        :: ph
       character(len=132)    :: line !lowline,
       !character(len=30)    :: magmod, shubk
       !character(len=2)     :: lattice, chardom
       !character(len=4)     :: symbcar
       character(Len=*),dimension(4),parameter :: cod=(/"a","b","c",";"/)
       real(kind=cp), dimension(3,3) :: Mat  !Matrix from parent space group to actual setting in magnetic cell
       real(kind=cp), dimension(3)   :: v    !Change of origin with respect to the standard
       !
       call Init_MagSymm_k_Type(MG_Symk)

       full_convertion=.false.
       if(MSpG%parent_num > 0 .or. len_trim(MSpG%Parent_spg) /= 0) Then
         if(len_trim(MSpG%trn_from_parent) /= 0) then !The transformation from the parent space group to the
            call Get_Transf(MSpG%trn_from_parent,mat,v,cod)  !actual given cell is read from this item
            if(Err_String) then
               full_convertion=.false.
            else
               full_convertion=.true.
            end if
         end if
       end if

       Select Case (l_case(Mode(1:2)))
         Case("mc")  !Use magnetic cell is equivalent to use a k=0 in MG_Symk
             if(MSpG%m_cell) then
                MG_Symk%MagModel  ="Using the magnetic cell "
                MG_Symk%Sk_Type   ="Crystal_Frame"
                MG_Symk%BNS_number=MSpG%BNS_number
                MG_Symk%OG_number =MSpG%OG_number
                MG_Symk%BNS_symbol=MSpG%BNS_symbol
                MG_Symk%OG_symbol =MSpG%OG_symbol
                MG_Symk%MagType   =MSpG%MagType
                MG_Symk%Parent_num=MSpG%Parent_num
                MG_Symk%Parent_spg=MSpG%Parent_spg
                MG_Symk%Latt="P"
                MG_Symk%nmsym=1
                MG_Symk%nirreps=0
                MG_Symk%centred=1    !By default the crystal structure is acentric
                MG_Symk%mcentred=1   !By default the magnetic structure is anti-centric (if there is -1 it is combined with time inversion)
                MG_Symk%nkv=1        !always a propagation vector even if not provided in MSpG
                MG_Symk%kvec=0.0     !The propagation vector is assumed to be (0,0,0) w.r.t. "magnetic cell"
                MG_Symk%Num_Lat=1    !No lattice centring are considered (all of them should be included in the list of operators)
                MG_Symk%Ltr=0.0
                MG_Symk%Numops=MSpG%Multip  !Reduced number of symmetry operators is equal to the multiplicity in this
                MG_Symk%Multip=MSpG%Multip  !case ...
                allocate(MG_Symk%Symop(MSpG%Multip))
                allocate(MG_Symk%SymopSymb(MSpG%Multip))
                allocate(MG_Symk%MSymop(MSpG%Multip,1))
                allocate(MG_Symk%MSymopSymb(MSpG%Multip,1))
                MG_Symk%Symop=MSpG%Symop
                do i=1,MG_Symk%Multip
                  MG_Symk%MSymop(i,1)=MSpG%MSymop(i)
                end do
                if(MSpG%mcif) then
                    do i=1,MG_Symk%Multip
                       line=MSpG%MSymopSymb(i)
                       do k=1,len_trim(line)
                         if(line(k:k) == "m") line(k:k)=" "
                         if(line(k:k) == "x") line(k:k)="u"
                         if(line(k:k) == "y") line(k:k)="v"
                         if(line(k:k) == "z") line(k:k)="w"
                       end do
                       line=Pack_String(line)//", 0.0"
                       MG_Symk%MSymopSymb(i,1)=trim(line)
                    end do
                else
                    do i=1,MG_Symk%Multip
                      MG_Symk%MSymopSymb(i,1)=MSpG%MSymopSymb(i)
                    end do
                end if
             else
                Err_Magsym=.true.
                Err_Magsym_Mess=" The magnetic cell in the Magnetic_Space_Group_Type is not provided! Use mode CC!"
                return
             end if

         Case("cc")
             !This only possible if full_convertion is true
             if(full_convertion) then
                MG_Symk%MagModel  = "Using the crystal cell and k-vectors "
                MG_Symk%Sk_Type   = "Crystal_Frame"
                MG_Symk%BNS_number= MSpG%BNS_number
                MG_Symk%OG_number = MSpG%OG_number
                MG_Symk%BNS_symbol= MSpG%BNS_symbol
                MG_Symk%OG_symbol = MSpG%OG_symbol
                MG_Symk%MagType   = MSpG%MagType
                MG_Symk%Parent_num= MSpG%Parent_num
                MG_Symk%Parent_spg= MSpG%Parent_spg
                MG_Symk%Latt      = MSpG%Parent_spg(1:1)
                MG_Symk%nmsym     = 1
                MG_Symk%nirreps=0

                MG_Symk%centred=1    !By default the crystal structure is acentric
                MG_Symk%mcentred=1   !By default the magnetic structure is anti-centric (if there is -1 it is combined with time inversion)
                MG_Symk%nkv=1        !always a propagation vector even if not provided in MSpG
                MG_Symk%kvec=0.0     !The propagation vector is assumed to be (0,0,0) w.r.t. "magnetic cell"
                MG_Symk%Num_Lat=1    !No lattice centring are considered (all of them should be included in the list of operators)
                MG_Symk%Ltr=0.0
                MG_Symk%Numops=MSpG%Multip  !Reduced number of symmetry operators is equal to the multiplicity in this
                MG_Symk%Multip=MSpG%Multip  !case ...
                allocate(MG_Symk%Symop(MSpG%Multip))
                allocate(MG_Symk%SymopSymb(MSpG%Multip))
                allocate(MG_Symk%MSymop(MSpG%Multip,1))
                allocate(MG_Symk%MSymopSymb(MSpG%Multip,1))
                MG_Symk%Symop=MSpG%Symop
                do i=1,MG_Symk%Multip
                  MG_Symk%MSymop(i,1)=MSpG%MSymop(i)
                end do
                if(MSpG%mcif) then
                    do i=1,MG_Symk%Multip
                       line=MSpG%MSymopSymb(i)
                       do k=1,len_trim(line)
                         if(line(k:k) == "m") line(k:k)=" "
                         if(line(k:k) == "x") line(k:k)="u"
                         if(line(k:k) == "y") line(k:k)="v"
                         if(line(k:k) == "z") line(k:k)="w"
                       end do
                       line=Pack_String(line)//", 0.0"
                       MG_Symk%MSymopSymb(i,1)=trim(line)
                    end do
                else
                    do i=1,MG_Symk%Multip
                      MG_Symk%MSymopSymb(i,1)=MSpG%MSymopSymb(i)
                    end do
                end if
             else
                Err_Magsym=.true.
                Err_Magsym_Mess=&
                " This option is available only if the parent group and the transformation has ben provided! Use mode MC!"
                return
             end if

       End Select
       return
    End Subroutine Magnetic_Space_Group_Type_to_MagSymm_k_Type

    !!  WARNING this subroutine is not operational!!!
    Subroutine MagSymm_k_Type_to_Magnetic_Space_Group_Type(MG_Symk,MSpG)
       Type(MagSymm_k_Type),              intent(in)  :: MG_Symk
       Type(Magnetic_Space_Group_Type),   intent(out) :: MSpG
       !---- Local variables ----!
       Type(Space_Group_Type) :: SpG
       integer :: i,m,n, ngen
       character(len=40),dimension(:), allocatable   :: gen
       !character(len=132)   :: lowline
       !character(len=30)    :: magmod, shubk
       !character(len=2)     :: lattice, chardom
       !character(len=4)     :: symbcar

       !
       call Init_Magnetic_Space_Group_Type(MSpG)

       !Verify the crystal structure information contained in MG_Symk by constructing the full Space group
       n=MG_Symk%Numops
       m=MG_Symk%Numops*MG_Symk%centred*MG_Symk%Num_Lat
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
       real(kind=cp), parameter      :: epsi=0.00001
       real(kind=cp)                 :: ph
       real(kind=cp),dimension(3)    :: rsk,isk,car,side
       real(kind=cp),dimension(3,12) :: br,bi
       real(kind=cp),dimension(3,3)  :: cart_to_cryst
       real(kind=cp),dimension(12)   :: coef
       character(len=132)            :: lowline,line
       character(len=50)             :: magmod, shubk
       character(len=2)              :: lattice, chardom
       character(len=4)              :: symbcar
       character(len=50)             :: msyr
       logical                       :: msym_begin, kvect_begin, skp_begin, shub_given, irreps_given, &
                                        irreps_begin, bfcoef_begin, magdom_begin, done, spg_given, lattice_given, &
                                        kvec_given
       type(Magnetic_Group_Type)     :: SG
       type(Space_Group_Type)        :: SpG

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
       call Init_MagSymm_k_Type(MGp)

       !Determine the number of symmetry operators existing the magnetic part of the file
       !This is for allocating the dimension of allocatable arrays in MagSymm_k_Type object
       !We will allocate the double for taking into account the possible centring of the magnetic structure
       n=0
       done=.false.
       do i=n_ini,n_end
          lowline=l_case(adjustl(file_cfl%line(i)))
          if (index(lowline(1:4),"symm") == 0 ) cycle
          n=n+1

          !determine now the number of msym cards per symm card
          if(.not. done) then
            m=0
            do j=i+1,i+8
               lowline=l_case(adjustl(file_cfl%line(j)))
               if (index(lowline(1:4),"msym") /= 0 ) then
                 m=m+1
                 cycle
               end if
               done=.true.
               exit
            end do
          end if

       end do
       n=2*n !if it is centred we will need this space
       if(n > 0) then
          !Allocate the allocatable components of MagSymm_k_Type
          if(allocated(MGp%Symop))      deallocate(MGp%Symop)
          if(allocated(MGp%SymopSymb))  deallocate(MGp%SymopSymb)
          if(allocated(MGp%MSymop))     deallocate(MGp%MSymop)
          if(allocated(MGp%MSymopSymb)) deallocate(MGp%MSymopSymb)
          allocate(MGp%Symop(n))
          allocate(MGp%SymopSymb(n))
          if(m > 0) then
             allocate(MGp%MSymop(n,m))
             allocate(MGp%MSymopSymb(n,m))
          end if
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
       i=n_ini
       shub_given  =.false.
       irreps_given=.false.
       irreps_begin=.false.
       msym_begin  =.false.
       skp_begin   =.false.
       bfcoef_begin=.false.
       spg_given   =.false.
       lattice_given=.false.
       kvec_given=.false.
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

          ! Read magnetic field for paramagnetic-induced magnetic moments
          ! The first item is the strength of the magnetic field in Tesla and the three other
          ! items correspond to the vector (in crystallographic space) of the direction of applied field.
          if (lowline(1:9) == "mag_field") then
             read(unit=lowline(10:),fmt=*,iostat=ier) Am%MagField, Am%dir_MField
             if (ier /= 0) then
                err_magsym=.true.
                ERR_MagSym_Mess=" Error reading magnetic field in magnetic phase"
                return
             end if
             if( Am%MagField > 0.0001) Am%suscept=.true.
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
             lattice_given=.true.
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
             kvec_given=.true.
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

          ! Read Local Susceptibility coefficients in cryst. axes
          if (lowline(1:3) == "chi") then
             read(unit=lowline(4:),fmt=*,iostat=ier) coef(1:6)
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=ERR_MagSym_Mess,fmt="(a,i3)") " Error reading Local Susceptibility Coefficient for atom #", num_matom
                return
             end if
               Am%atom(num_matom)%chi= coef(1:6)
               if(abs(coef(1)-coef(2)) < epsi .and. abs(coef(2)-coef(3)) < epsi .and. sum(abs(coef(4:6))) < epsi ) then
                 Am%atom(num_matom)%chitype="isotr"
               else
                 Am%atom(num_matom)%chitype="aniso"
               end if
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
                !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
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
       if(num_matom == 0 .and. .not. kvec_given .and. .not. shub_given) then ! No information on magnetic structure is really provided
         err_magsym=.true.
         ERR_MagSym_Mess= " No information on magnetic structure is really provided! "
         return
       end if

       !Check if it is an induced paramagnetic magnetic structure due to an applied magnetic field
       !In such a case use the crystal space group to construct the magnetic matrices. If the symbol
       !of the space group is not provided it is supposed that the symmetry operators have been provided
       !together with th SYMM and MSYM matrices
       if(Am%suscept) then
         do i=1,file_cfl%nlines
            lowline=adjustl(l_case(file_cfl%line(i)))
            if (lowline(1:4) == "spgr" .or. lowline(1:3) == "spg" .or. lowline(1:6) == "spaceg") then
               lowline=adjustl(file_cfl%line(i))
               j=index(lowline," ")
               lowline=lowline(j+1:)
               call Set_SpaceGroup(trim(lowline),SpG)
               spg_given   =.true.
               exit
            end if
         end do
         if(spg_given) then
            n=SpG%Numops * SpG%Centred
            MGp%Centred=SpG%Centred
            MGp%MCentred=1  !Same rotation matrices as that of the space group
            MGp%Latt=SpG%SPG_Symb(1:1)
            num_xsym=SpG%Numops
            num_msym=1
            num_k=1
            if(allocated(MGp%Symop))      deallocate(MGp%Symop)
            if(allocated(MGp%SymopSymb))  deallocate(MGp%SymopSymb)
            if(allocated(MGp%MSymop))     deallocate(MGp%MSymop)
            if(allocated(MGp%MSymopSymb)) deallocate(MGp%MSymopSymb)
            allocate(MGp%Symop(n))
            allocate(MGp%SymopSymb(n))
            allocate(MGp%MSymop(n,1))
            allocate(MGp%MSymopSymb(n,1))

            do i=1,num_xsym
              lowline=" "
              MGp%SymopSymb(i)=SpG%SymopSymb(i)
              call Get_SymSymb(Spg%Symop(i)%Rot,(/0.0,0.0,0.0/),lowline)
              do j=1,len_trim(lowline)
                 if(lowline(j:j) == "x") lowline(j:j) = "u"
                 if(lowline(j:j) == "y") lowline(j:j) = "v"
                 if(lowline(j:j) == "z") lowline(j:j) = "w"
              end do
              MGp%MSymopSymb(i,1) = trim(lowline)//", 0.0"
            end do
         else
           !No action to be taken ... the symmetry operators are read in SYMM cards
         end if
       end if

       !Get pointers to the magnetic form factors
       !Stored for each atom in the component ind(1)
       call Set_Magnetic_Form()

       !---- Find Species in Magnetic_Form ----!
       do i=1,Am%natoms
          symbcar=u_case(Am%atom(i)%SfacSymb)
          do j=1,num_mag_form
             if (symbcar /= Magnetic_Form(j)%Symb) cycle
             Am%atom(i)%ind(2)=j
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
       Mgp%Num_Lat=1
       MGp%Ltr(:,:) = 0.0
       Select Case(MGp%Latt)
          case ("A")
             Mgp%Num_Lat=2
             MGp%Ltr(:,1:2)=Ltr_a(:,1:2)
          case ("B")
             Mgp%Num_Lat=2
             MGp%Ltr(:,1:2)=Ltr_b(:,1:2)
          case ("C")
             Mgp%Num_Lat=2
             MGp%Ltr(:,1:2)=Ltr_c(:,1:2)
          case ("I")
             Mgp%Num_Lat=2
             MGp%Ltr(:,1:2)=Ltr_i(:,1:2)
          case ("R")
             Mgp%Num_Lat=3
             MGp%Ltr(:,1:3)=Ltr_r(:,1:3)
          case ("F")
             Mgp%Num_Lat=4
             MGp%Ltr(:,1:4)=Ltr_f(:,1:4)
       End Select

       select case (MGp%centred)
          case (1)
             MGp%Multip =   MGp%Numops * Mgp%Num_Lat
          case (2)
             MGp%Multip = 2 * MGp%Numops * Mgp%Num_Lat
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
    !!----  Updated: August-2014 (JRC)
    !!
    Subroutine Readn_Set_Magnetic_Structure_MCIF(file_mcif,mCell,MGp,Am)
       character(len=*),               intent (in)  :: file_mcif
       type(Crystal_Cell_type),        intent (out) :: mCell
       type(Magnetic_Space_Group_Type),intent (out) :: MGp
       type(mAtom_List_Type),          intent (out) :: Am

       !---- Local Variables ----!
       integer :: i,num_sym, num_constr, num_kvs,num_matom, num_mom, num_magscat, ier, j, m, n, k, L,   &
                  ncar,mult,nitems,iv, num_irreps, nitems_irreps, num_rsym, num_centering,det,kfin
       integer,          dimension(10)     :: lugar
       integer,          dimension(7)      :: irrep_pos
       integer,          dimension(5)      :: pos
       integer,          dimension(3,3)    :: Rot
       real(kind=cp),    dimension(3)      :: cel,ang,cel_std,ang_std,tr,v
       real(kind=cp),    dimension(6)      :: values,std
       real(kind=cp),    dimension(3,3)    :: matr
       real(kind=cp),    dimension(3,384)  :: orb
       character(len=132)                  :: lowline,keyword,line, mxmymz_op,linat
       character(len=132),dimension(384)   :: sym_strings, cent_strings
       character(len=132),dimension(384)   :: atm_strings
       character(len=132),dimension(384)   :: mom_strings
       character(len=132),dimension(30)    :: constr_strings, mag_scatt_string
       character(len=132),dimension(30)    :: irreps_strings
       character(len=132),dimension(30)    :: kv_strings
       character(len=20), dimension(15)    :: lab_items
       character(len=50)                   :: shubk
       character(len=2)                    :: chars
       character(len=10)                   :: label
       character(len=4)                    :: symbcar
       logical                             :: ktag,no_symop_mxmymz,no_cent_mxmymz,mom_symmform

       !type(Magnetic_Group_Type)  :: SG
       type(file_list_type)       :: mcif

       call init_err_MagSym()
       call File_To_FileList(file_mcif,mcif)
       !Remove all possible tabs and non-ASCII characters in the CIF
       do i=1,mcif%nlines
         do j=1,len_trim(mcif%line(i))
           if(mcif%line(i)(j:j) == char(9)) mcif%line(i)(j:j)=" "
         end do
       end do
       num_constr=0; num_kvs=0; num_matom=0; num_mom=0; num_sym=0; num_magscat=0; num_rsym=0; num_centering=0
       cel=0.0; ang=0.0; num_irreps=0; nitems_irreps=0
       i=0
       call Init_Magnetic_Space_Group_Type(MGp)
       ktag=.false.
       no_symop_mxmymz=.false.
       no_cent_mxmymz=.false.
       mom_symmform=.false.

       do
          i=i+1
          if(i > mcif%nlines) exit
          if (index(mcif%line(i)(1:1),"!")/=0 .or. index(mcif%line(i)(1:1),"#")/=0 .or. len_trim(mcif%line(i)) == 0) cycle
          line=adjustl(mcif%line(i))
          lowline=l_case(line)
          j=index(lowline," ")
          keyword=lowline(1:j-1)
          !write(*,"(a)") " Keyword: "//trim(keyword)

          Select Case (trim(keyword))

             Case("_magnetic_space_group_standard_setting","_magnetic_space_group.standard_setting")
                chars=adjustl(line(j+1:))
                if(chars(2:2) == "y" .or. chars(2:2) == "Y") MGp%standard_setting=.true.
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_standard_setting -> "//trim(chars)

             Case("_parent_space_group.name_h-m", "_parent_space_group_name_h-m","_parent_space_group.name_h-m_alt")
                shubk=adjustl(line(j+1:))
                m=len_trim(shubk)
                MGp%Parent_spg=shubk(2:m-1)
                !write(unit=*,fmt="(a)") "  Treating item: _parent_space_group_name_h-m -> "// MGp%Parent_spg

             Case("_parent_space_group.it_number","_parent_space_group_it_number")
                read(unit=lowline(j:),fmt=*,iostat=ier) m
                if(ier /= 0) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the number of the parent space group"
                  return
                end if
                MGp%Parent_num=m
                !write(unit=*,fmt="(a,i4)") "  Treating item: _parent_space_group_it_number -> ", MGp%Parent_num

             Case("_magnetic_space_group_bns_number","_space_group.magn_number_bns","_space_group_magn.number_bns")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%BNS_number=shubk
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_number_bns -> "//trim(MGp%BNS_number)

             Case("_magnetic_space_group_bns_name","_space_group_magn.name_bns","_space_group.magn_name_bns")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%BNS_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_name_bns -> "//trim(MGp%BNS_symbol)

             Case("_magnetic_space_group_og_number","_space_group_magn.number_og","_space_group.magn_number_og")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%OG_number=shubk
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_number_og -> "//trim(MGp%OG_number)

             Case("_magnetic_space_group_point_group","_space_group_magn.point_group","_space_group.magn_point_group")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%PG_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group_magn.point_group -> "//trim(MGp%PG_symbol)

             Case("_magnetic_space_group_og_name","_space_group_magn.name_og","_space_group.magn_name_og")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%OG_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_name_og -> "//trim(MGp%OG_symbol)

             Case("_magnetic_space_group.transform_from_parent_pp_abc","_magnetic_space_group_transform_from_parent_pp_abc", &
                   "_parent_space_group.child_transform_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_from_parent=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_transform_from_parent_pp_abc -> "//trim(MGp%trn_from_parent)

             Case("_parent_space_group.transform_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_parent=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_transform_from_parent_pp_abc -> "//trim(MGp%trn_from_parent)

             Case("_magnetic_space_group.transform_to_standard_pp_abc","_magnetic_space_group_transform_to_standard_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_standard=pack_string(shubk)

             Case("_space_group_magn.transform_bns_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_standard=pack_string(shubk)

             Case("_magnetic_cell_length_a","_cell_length_a")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'a' -> "//trim(err_string_mess)
                  return
                end if
                cel(1)=values(1)
                cel_std(1)=std(1)
                MGp%m_cell=.true.
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_a"

             Case("_magnetic_cell_length_b","_cell_length_b")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'b' -> "//trim(err_string_mess)
                  return
                end if
                cel(2)=values(1)
                cel_std(2)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_b"

             Case("_magnetic_cell_length_c","_cell_length_c")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'c' -> "//trim(err_string_mess)
                  return
                end if
                cel(3)=values(1)
                cel_std(3)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_c"

             Case("_magnetic_cell_angle_alpha","_cell_angle_alpha")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'alpha' -> "//trim(err_string_mess)
                  return
                end if
                ang(1)=values(1)
                ang_std(1)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_alpha"

             Case("_magnetic_cell_angle_beta","_cell_angle_beta")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'beta' -> "//trim(err_string_mess)
                  return
                end if
                ang(2)=values(1)
                ang_std(2)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_beta"

             Case("_magnetic_cell_angle_gamma","_cell_angle_gamma")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  err_magsym=.true.
                  ERR_MagSym_Mess=" Error reading the magnetic unit cell parameter 'gamma' -> "//trim(err_string_mess)
                  return
                end if
                ang(3)=values(1)
                ang_std(3)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_gamma"

             Case("loop_")
                 i=i+1
                 line=adjustl(mcif%line(i))
                 lowline=l_case(line)
                 j=index(lowline," ")
                 keyword=lowline(1:j-1)
                 !write(*,"(a)") "         Loop_Keyword: "//trim(keyword)
                 Select Case(trim(keyword))

                   Case("_space_group_magn_transforms.id")
                      !write(*,"(a)") "         Loop_Keyword: "//trim(keyword)

                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_magn_transforms") == 0) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading _space_group_magn_transforms in loop"
                          return
                        end if
                      end do
                      i=i+1
                      call getword(mcif%line(i),lab_items,iv)
                      !write(unit=*,fmt="(3a)")  (lab_items(k),k=1,3)
                      if(lab_items(3)(1:3) == "BNS") then
                        MGp%trn_to_standard=lab_items(2)
                      end if
                      i=i+1
                      call getword(mcif%line(i),lab_items,iv)
                      !write(unit=*,fmt="(3a)")  (lab_items(k),k=1,3)
                      if(lab_items(3)(1:2) == "OG") then
                        !nothing to do
                      end if

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
                         if(index(mcif%line(i),"_small_irrep_dimension") /= 0 .or.  &
                            index(mcif%line(i),"_irrep_small_dimension") /= 0) then
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
                         if(index(mcif%line(i),"_irrep_presence") /= 0) then
                            j=j+1
                            irrep_pos(7)=j
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

                   Case("_atom_type_symbol")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_type_symbol"
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_atom_type_symbol") == 0) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the _atom_type_symbol in loop"
                          return
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        mag_scatt_string(k)=mcif%line(i)
                      end do
                      num_magscat=k
                      !Treat later the scattering factor

                   Case("_magnetic_atom_site_moment_symmetry_constraints_label")
                      !write(unit=*,fmt="(a)") "  Treating item: _magnetic_atom_site_moment_symmetry_constraints_label"
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
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_operation.id"
                      do k=1,3
                        i=i+1
                        j=index(mcif%line(i),"_magnetic_space_group_symop_operation")
                        if(j == 0 ) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the _magnetic_space_group_symop_operation loop"
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
                      MGp%Multip=k

                   Case("_space_group_symop_magn_operation.id","_space_group_symop_magn.id") !The second item is added to be compatible with BCS error
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_operation.id"

                      i=i+1
                      j=index(mcif%line(i),"_space_group_symop_magn_operation.xyz")
                      if(j == 0 ) then
                        err_magsym=.true.
                        ERR_MagSym_Mess=" Error reading the _space_group_symop_magn_operation loop"
                        return
                      end if

                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do
                      !now allocate the list of symmetry operators
                      num_sym=k
                      MGp%Multip=k

                   Case("_space_group_symop.magn_id")
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_id"
                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_symop.magn_operation") == 0 .and. &
                           index(mcif%line(i),"_space_group_symop_magn_operation") == 0) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the _space_group_symop_magn_operation loop"
                          return
                        end if
                      end do
                      if(index(mcif%line(i),"_space_group_symop.magn_operation_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation_mxmymz") == 0) then
                         i=i-1
                         no_symop_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do

                      num_rsym=k

                   Case("_space_group_symop_magn_id")   !here the symmetry operators are separated from the translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_id"
                      i=i+1
                      if(index(mcif%line(i),"_space_group_symop.magn_operation") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation") == 0) then
                        err_magsym=.true.
                        ERR_MagSym_Mess=" Error reading the _space_group_symop.magn_operation loop"
                        return
                      end if
                      if(index(mcif%line(i),"_space_group_symop.magn_operation_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation_mxmymz") == 0) then
                         no_symop_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do

                      num_rsym=k

                   Case("_space_group_symop.magn_centering_id")   !here we read the translations and anti-translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_centering_id"
                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_symop.magn_centering") == 0 .and. &
                           index(mcif%line(i),"_space_group_symop_magn_centering") == 0 ) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the _space_group_symop_magn_centering loop"
                          return
                        end if
                      end do
                      if(index(mcif%line(i),"_space_group_symop.magn_centering_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_centering_mxmymz") == 0 ) then
                         i=i-1
                         no_cent_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        cent_strings(k)=mcif%line(i)
                      end do
                      num_centering=k

                   Case("_space_group_symop_magn_centering.id")   !here we read the translations and anti-translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_centering_id"
                      i=i+1
                      if(index(mcif%line(i),"_space_group_symop_magn_centering.xyz") == 0) then
                        err_magsym=.true.
                        ERR_MagSym_Mess=" Error reading the _space_group_symop_magn_centering.xyz loop"
                        return
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        cent_strings(k)=mcif%line(i)
                      end do
                      num_centering=k

                   Case("_magnetic_atom_site_label","_atom_site_label")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_site_label"
                      !Count the number of keywords following the _loop
                      do k=1,10
                        linat=adjustl(mcif%line(i+k))
                        if(linat(1:1) /=  "_") then
                          kfin=k+1
                          iv=i+k
                          exit
                        end if
                      end do
                      lugar=0
                      lugar(1)=1
                      j=1
                      do k=1,kfin
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
                         if (index(mcif%line(i),"_atom_site_B_iso_or_equiv") /= 0) then
                            j=j+1
                            lugar(10)=j
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

                      i=iv-1
                      nitems=count(lugar > 0)

                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        atm_strings(k)=adjustl(mcif%line(i))
                      end do
                      num_matom=k
                      !Treat late the list atoms

                   Case("_magnetic_atom_site_moment_label","_atom_site_moment_label","_atom_site_moment.label")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_site_moment_label"
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_atom_site_moment_crystalaxis") == 0 .and. &
                           index(mcif%line(i),"_atom_site_moment.crystalaxis") == 0) then
                          err_magsym=.true.
                          ERR_MagSym_Mess=" Error reading the magnetic_atom_site_moment loop"
                          return
                        end if
                      end do
                      i=i+1
                      if(index(mcif%line(i),"_atom_site_moment.symmform") /= 0) then
                        !write(*,*) " _atom_site_moment.symmform FOUND"
                        mom_symmform=.true.
                      else
                        i=i-1
                      end if
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
       !write(unit=*,fmt="(a,2i4)") " num_sym, num_rsym :",num_sym,num_rsym
       if(num_sym == 0 .and. num_rsym == 0) then
          err_magsym=.true.
          ERR_MagSym_Mess=" No symmetry operators have been provided in the MCIF file "//trim(file_mcif)
          return
       else
          if(no_cent_mxmymz) then  !Full number of symmetry operators is not separated from the centering

            if(allocated(Mgp%SymopSymb)) deallocate(Mgp%SymopSymb)
            allocate(Mgp%SymopSymb(num_sym))
            if(allocated(Mgp%Symop)) deallocate(Mgp%Symop)
            allocate(Mgp%Symop(num_sym))
            if(allocated(Mgp%MSymopSymb)) deallocate(Mgp%MSymopSymb)
            allocate(Mgp%MSymopSymb(num_sym))
            if(allocated(Mgp%MSymop)) deallocate(Mgp%MSymop)
            allocate(Mgp%MSymop(num_sym))
            !write(unit=*,fmt="(a)") "  Decoding symmetry operators 1"

            ! Decode the symmetry operators
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

          else

            if( num_rsym == 0) num_rsym=num_sym
            ! First allocate the full number of symmetry operators after decoding if centering lattice
            ! have been provided and if the group is centred or not
            if(num_centering == 0) then
               MGp%Multip=num_rsym
            else
               MGp%Multip=num_rsym*num_centering
            end if

            num_sym=MGp%Multip
            if(allocated(Mgp%SymopSymb)) deallocate(Mgp%SymopSymb)
            allocate(Mgp%SymopSymb(num_sym))
            if(allocated(Mgp%Symop)) deallocate(Mgp%Symop)
            allocate(Mgp%Symop(num_sym))
            if(allocated(Mgp%MSymopSymb)) deallocate(Mgp%MSymopSymb)
            allocate(Mgp%MSymopSymb(num_sym))
            if(allocated(Mgp%MSymop)) deallocate(Mgp%MSymop)
            allocate(Mgp%MSymop(num_sym))
            ! Decode the symmetry operators
            !write(unit=*,fmt="(a)") "  Decoding symmetry operators 2"
            do i=1,num_rsym
              line=adjustl(sym_strings(i))
              j=index(line," ")
              line=adjustl(line(j+1:))
              j=index(line," ")
              MGp%SymopSymb(i)=line(1:j-1)
              k=index(MGp%SymopSymb(i),",",back=.true.)
              read(unit=MGp%SymopSymb(i)(k+1:),fmt=*,iostat=ier) n
              if(ier /= 0) then
                 err_magsym=.true.
                 ERR_MagSym_Mess=" Error reading the time inversion in line: "//trim(sym_strings(i))
                 return
              else
                 MGp%MSymOp(i)%phas=real(n)
              end if
              MGp%SymopSymb(i)=MGp%SymopSymb(i)(1:k-1)
              call Read_Xsym(MGp%SymopSymb(i),1,MGp%Symop(i)%Rot,MGp%Symop(i)%tr)

              !Now construc the magnetic rotation symbols
              line=adjustl(line(j+1:))
              if(len_trim(line) /= 0) then
                j=index(line," ")
                MGp%MSymopSymb(i)=line(1:j-1)
                line=MGp%MSymopSymb(i)
                do k=1,len_trim(line)
                  if(line(k:k) == "m") line(k:k)=" "
                end do
                line=Pack_String(line)
                call Read_Xsym(line,1,MGp%MSymop(i)%Rot)
              else
                det=determ_a(MGp%Symop(i)%Rot)
                MGp%MSymop(i)%Rot=MGp%Symop(i)%Rot*det*nint(MGp%MSymOp(i)%phas)
                call Get_Symsymb(MGp%MSymOp(i)%Rot,(/0.0,0.0,0.0/),line)
                !Expand the operator "line" to convert it to mx,my,mz like
                mxmymz_op=" "
                do j=1,len_trim(line)
                  Select Case(line(j:j))
                    case("x")
                       mxmymz_op=trim(mxmymz_op)//"mx"
                    case("y")
                       mxmymz_op=trim(mxmymz_op)//"my"
                    case("z")
                       mxmymz_op=trim(mxmymz_op)//"mz"
                    case default
                       mxmymz_op=trim(mxmymz_op)//line(j:j)
                  End Select
                end do
                MGp%MSymopSymb(i)=trim(mxmymz_op)
              end if


            end do
            !Decode lattice translations and anti-translations

            !write(unit=*,fmt="(a)") "  Decoding lattice translations and anti-translations"
            m=num_rsym
            do L=2,num_centering
              line=adjustl(cent_strings(L))
              j=index(line," ")
              line=adjustl(line(j+1:))
              j=index(line," ")
              line=line(1:j-1)
              k=index(line,",",back=.true.)
              read(unit=line(k+1:),fmt=*,iostat=ier) n
              if(ier /= 0) then
                 err_magsym=.true.
                 ERR_MagSym_Mess=" Error reading the time inversion in line: "//trim(cent_strings(i))
                 return
              end if
              line=line(1:k-1)
              call Read_Xsym(line,1,Rot,tr)

              do j=1,num_rsym
                m=m+1
                v=MGp%SymOp(j)%tr(:) + tr
                MGp%SymOp(m)%Rot  = MGp%SymOp(j)%Rot
                MGp%SymOp(m)%tr   = modulo_lat(v)
                MGp%MSymOp(m)%Rot = n*MGp%MSymOp(j)%Rot
                MGp%MSymOp(m)%phas= n*MGp%MSymOp(j)%phas
                call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
                call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
                !Expand the operator "line" to convert it to mx,my,mz like
                mxmymz_op=" "
                do i=1,len_trim(line)
                  Select Case(line(i:i))
                    case("x")
                       mxmymz_op=trim(mxmymz_op)//"mx"
                    case("y")
                       mxmymz_op=trim(mxmymz_op)//"my"
                    case("z")
                       mxmymz_op=trim(mxmymz_op)//"mz"
                    case default
                       mxmymz_op=trim(mxmymz_op)//line(i:i)
                  End Select
                end do
                MGp%MSymopSymb(m)=trim(mxmymz_op)
              end do
            end do
          end if
       end if
       ! Symmetry operators treatment done
       Call cleanup_symmetry_operators(MGp)
       if(err_MagSym) then
          return
          !write(unit=*,fmt="(a)") " => "//trim(err_MagSym)
       end if

       !Treating irreps

       if(num_irreps == 0) then

          MGp%n_irreps=0

       else
          !write(*,"(a,i3)") " Treating irreps: ",num_irreps
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
         !write(*,"(a,i3)") " Treating propagation vectors: ",num_kvs
         do i=1,MGp%n_kv
            line=adjustl(kv_strings(i))
            j=index(line," ")
            MGp%kv_label(i)=line(1:j-1)
            line=adjustl(line(j+1:))
            n=len_trim(line)
            if(ktag) then
              line=adjustl(line(2:n-1))
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
          !write(*,"(a,i4)") " Treating magnetic atoms:  ",num_matom
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

            if (lugar(6) /= 0) then  ! _atom_site_U_iso_or_equiv
               call getnum_std(lab_items(lugar(6)),values,std,iv)
               Am%atom(i)%ueq=values(1)
               Am%atom(i)%Biso=values(1)*78.95683521     !If anisotropic they
               Am%atom(i)%Biso_std=std(1)*78.95683521    !will be put to zero
            else if (lugar(10) /= 0) then    ! _atom_site_B_iso_or_equiv
               call getnum_std(lab_items(lugar(10)),values,std,iv)
               Am%atom(i)%ueq=values(1)/78.95683521
               Am%atom(i)%Biso=values(1)     !If anisotropic they
               Am%atom(i)%Biso_std=std(1)    !will be put to zero
            else
               Am%atom(i)%ueq=0.0
               Am%atom(i)%Biso=0.0
               Am%atom(i)%Biso_std=0.0
            end if
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
            Am%atom(i)%occ=Am%atom(i)%occ*real(Mult)/real(MGp%Multip)

            if(lugar(9) /= 0) then
               Am%atom(i)%wyck=adjustl(trim(lab_items(lugar(9))))
            end if

          end do
       end if

       !Treating moments of magnetic atoms
       if(num_mom /= 0) then
          !write(*,"(a,i4)") " Treating magnetic moments:  ",num_mom
          m=4
          if(mom_symmform) m=5
          do i=1,num_mom
            call getword(mom_strings(i),lab_items,iv)
            !write(*,"(2i6,tr4,5(a,tr3))") k,iv,lab_items(1:iv)
            if(iv /= m) then
               err_magsym=.true.
               write(unit=ERR_MagSym_Mess,fmt="(a,i4)")" Error reading magnetic moment #",i
               ERR_MagSym_Mess=trim(ERR_MagSym_Mess)//" -> 4-5 items expected in this line: 'Label mx my mz', read: "// &
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
                     Am%Atom(j)%SkI(k,1)=0.0
                     Am%Atom(j)%SkI_std(k,1)=0.0
                 end do
                 Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
               end if
            end do
          end do
       end if

       if(num_constr /= 0) then

         !write(*,"(a,i4)") " Treating constraints:  ",num_constr
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
                Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
                exit
             end if
           end do
           !The treatment of the codes will be done in the future
         end do
       end if

       if(num_magscat > 0) then !Reading the valence for determining the magnetic form factor
         do i=1,num_magscat
           call getword(mag_scatt_string(i),lab_items,iv)
           do j=1,Am%natoms
             if(Am%atom(j)%chemSymb == lab_items(1)) then
               Am%atom(j)%SfacSymb=lab_items(2)
               if(lab_items(2) /= ".") then !magnetic atoms
                  Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
               end if
             end if
           end do
         end do
       end if

       !Get pointers to the magnetic form factors
       !Stored for each atom in the component ind(1)
       call Set_Magnetic_Form()

       !---- Find Species in Magnetic_Form ----!
       do i=1,Am%natoms
          symbcar=get_magnetic_form_factor(Am%atom(i)%SfacSymb)
          do j=1,num_mag_form
             if (symbcar /= Magnetic_Form(j)%Symb) cycle
             Am%atom(i)%ind(2)=j
             exit
          end do
       end do

       return
    End Subroutine Readn_Set_Magnetic_Structure_MCIF

    !!
    !!----  Subroutine get_moment_ctr(xnr,moment,Spgr,codini,codes,ord,ss,att,Ipr)
    !!----     real(kind=cp), dimension(3),            intent(in    ) :: xnr    !Atom position (fractional coordinates)
    !!----     real(kind=cp), dimension(3),            intent(in out) :: moment !Moment at position xnr
    !!----     type(Magnetic_Space_Group_type),        intent(in    ) :: Spgr   !Magnetic Space Group
    !!----     Integer,                                intent(in out) :: codini !Last attributed parameter
    !!----     real(kind=cp), dimension(3),            intent(in out) :: codes  !codewords for positions
    !!----     integer,                       optional,intent(in)     :: ord
    !!----     integer, dimension(:),         optional,intent(in)     :: ss
    !!----     real(kind=cp), dimension(:,:), optional,intent(in)     :: att
    !!----     integer,                       optional,intent(in)     :: Ipr
    !!----
    !!----  Subroutine to get the appropriate constraints in the refinement codes of
    !!----  magnetic moment parameters.
    !!----  Algorithm based in the Wigner theorem.
    !!----  The vector Mom = Sum { R Moment} displays the symmetry constraints to be
    !!----  applied to the magnetic moments. The sum runs over all magnetic
    !!----  matrices of the stabilizer of the particular atom position in the given
    !!----  space group.
    !!----
    !!----   Updated: 16 April 2016
    !!----
    !!
    Subroutine get_moment_ctr_old(xnr,moment,Spgr,codini,codes,ord,ss,att,Ipr)
       real(kind=cp), dimension(3),            intent(in)     :: xnr
       real(kind=cp), dimension(3),            intent(in out) :: moment
       type(Magnetic_Space_Group_type),        intent(in)     :: Spgr
       Integer,                                intent(in out) :: codini
       real(kind=cp), dimension(3),            intent(in out) :: codes
       integer,                       optional,intent(in)     :: ord
       integer, dimension(:),         optional,intent(in)     :: ss
       real(kind=cp), dimension(:,:), optional,intent(in)     :: att
       integer,                       optional,intent(in)     :: Ipr

       ! Local variables
       character (len=4), dimension(3)   :: cdd
       character (len=4)                 :: cditem
       real(kind=cp),     dimension(3)   :: multip
       integer                           :: j,order
       real(kind=cp)                     :: suma,dif
       integer,           dimension(48)  :: ss_ptr
       integer,           dimension(3)   :: codd,msym
       real(kind=cp),     dimension(3,3) :: Rs
       real(kind=cp),     dimension(3)   :: x,cod,multi,mom,mome,Rsym
       real(kind=cp),     dimension(3,48):: atr
       real(kind=cp),     parameter      :: epss=0.01_cp

       suma=0.0
       do j=1,3
          suma=suma+abs(codes(j))
          cod(j)=int(abs(codes(j))/10.0_cp)             !Input Parameter number with sign
          multi(j)=mod(codes(j),10.0_cp)                !Input Multipliers
          if(cod(j) < 1.0 .and. abs(multi(j)) > epss)  then
               codini=codini+1
               cod(j) = real(codini)
          end if
       end do

       if(suma < epss) return  !No refinement is required
       x=modulo_lat(xnr)

       if(present(ord) .and. present(ss) .and. present(att)) then
         order=ord
         ss_ptr(1:order) = ss(1:ord)
         atr(:,1:order)  = att(:,1:ord)
       else
         call get_stabilizer(x,Spgr,order,ss_ptr,atr)
       end if

       mom=(/17.0, 7.0,5.0/)
       mome=mom
       if(present(ipr)) Write(unit=ipr,fmt="(a,i3)") " => Magnetic stabilizer without identity, order:",order
       if (order > 1 ) then
          do j=2,order
             Rs=real(Spgr%MSymOp(ss_ptr(j))%Rot)
             Rsym=matmul(Rs,mom)
             mome=mome+ Rsym
             if(present(ipr)) then
               write(unit=ipr,fmt='(a,i2,a,t20,a,t55,a,t75,3f8.1)') '     Operator ',j,": ",trim(Spgr%SymopSymb(ss_ptr(j))), &
                trim(Spgr%MSymopSymb(ss_ptr(j))), Rsym
             end if
          end do
          mome=mome/real(order)
       end if
       msym=nint(1000.0*mome)
       codd=msym
       cdd=(/'a','b','c'/)
       multip=1.0

       !Search systematically all the possible constraints

       if(codd(1) == codd(2) .and. codd(1) == codd(3)) then ! a a a
         cdd=(/'a','a','a'/)     ! 1 A A A
         multip=(/1.0,1.0,1.0/)
         moment(2:3)=moment(1)
         cod(2:3)=cod(1)
         if(codd(1) == 0) then !No magnetic moment allowed for this site
           cod=0
           moment=0.0
           multip=0.0
           cdd=(/'0','0','0'/)
         end if

       else if(codd(1) == codd(2)) then ! a a c
         cdd=(/'a','a','c'/)     ! 2  A A C
         multip=(/1.0,1.0,1.0/)
         moment(2)=moment(1)
         cod(2)=cod(1)
         if(codd(1) == 0) then ! 0 0 c
           cod(1:2)=0
           moment(1:2)=0.0
           multip(1:2)=0.0
           cdd=(/'0','0','c'/)
         else if(codd(3) == 0) then  ! a a 0
           cod(3)=0
           moment(3)=0.0
           multip(3)=0.0
           cdd=(/'a','a','0'/)
         else if(codd(3) == -codd(1)) then  ! a a -a
           cod(3)=cod(1)
           moment(3)=-moment(1)
           multip(3)=-1.0
           cdd=(/'a ','a ','-a'/)
         end if

       else if(codd(1) == codd(3)) then ! a b a
         cdd=(/'a','b','a'/)     ! 3  A B A
         multip=(/1.0,1.0,1.0/)
         moment(3)=moment(1)
         cod(3)=cod(1)
         if(codd(1) == 0) then !0 b 0
           cod(1)=0; cod(3)=0
           moment(1)=0.0; moment(3)=0.0
           multip(1)=0.0; multip(3)=0.0
           cdd=(/'0','b','0'/)
         else if(codd(2) == 0) then  ! a 0 a
           cod(2)=0
           moment(2)=0.0
           multip(2)=0.0
           cdd=(/'a','0','a'/)
         else if(codd(2) == -codd(1)) then  ! a -a a
           cod(2)=cod(1)
           moment(2)=-moment(1)
           multip(2)=-1.0
           cdd=(/'a ','-a','a '/)
         end if

       else if(codd(2) == codd(3)) then ! a b b
         cdd=(/'a','b','b'/)     ! 4  A B B
         multip=(/1.0,1.0,1.0/)
         moment(3)=moment(2)
         cod(3)=cod(2)
         if(codd(2) == 0) then !a 0 0
           cod(2:3)=0
           moment(2:3)=0.0
           multip(2:3)=0.0
           cdd=(/'a','0','0'/)
         else if(codd(1) == 0) then  ! 0 b b
           cod(1)=0
           moment(1)=0.0
           multip(1)=0.0
           cdd=(/'0','b','b'/)
         else if(codd(1) == -codd(2)) then  ! -b b b
           cod(1)=cod(2)
           moment(1)=-moment(2)
           multip(1)=-1.0
           cdd=(/'-b','b ','b '/)
         end if

       else !Now a /= b /= c

         if(codd(1) == 0) then  !0 b c
           cod(1)=0
           moment(1)=0.0
           multip(1)=0.0
           cdd=(/'0','b','c'/)
         end if
         if(codd(2) == 0) then  !a 0 c
           cod(2)=0
           moment(2)=0.0
           multip(2)=0.0
           cdd=(/'a','0','c'/)
         end if
         if(codd(3) == 0) then  !a b 0
           cod(3)=0
           moment(3)=0.0
           multip(3)=0.0
           cdd=(/'a','b','0'/)
         end if
         !Comparison a,b
         if(codd(1) /= 0 .and. codd(2)/=0) then
           suma=real(codd(1))/real(codd(2))
           if(abs(suma) < 1.0) then
             suma=1.0/suma
             order=codd(2)/codd(1)
             dif=abs(suma-real(order))
             if(dif < epss) then
               cod(2)=cod(1)
               multip(2)=suma
               moment(2)=suma*moment(1)
               write(unit=cditem,fmt="(i2,a)") order,"a"
               !cdd=(/'a',cditem,'c'/)  !incompatible with Lahey compiler
               cdd(1)='a'
               cdd(2)=cditem
               cdd(3)='c'
             end if
           else
             order=codd(1)/codd(2)
             dif=abs(suma-real(order))
             if(dif < epss) then
               cod(1)=cod(2)
               multip(1)=suma
               moment(1)=suma*moment(2)
               write(unit=cditem,fmt="(i2,a)") order,"b"
               !cdd=(/cditem,'b','c'/)
               cdd(1)=cditem
               cdd(2)='b'
               cdd(3)='c'
             end if
            end if
         end if
         !Comparison a,c
         if(codd(1) /= 0 .and. codd(3)/=0) then
           suma=real(codd(1))/real(codd(3))
           if(abs(suma) < 1.0) then
             suma=1.0/suma
             order=codd(3)/codd(1)
             dif=abs(suma-real(order))
             if(dif < epss) then
               cod(3)=cod(1)
               multip(3)=suma
               moment(3)=suma*moment(1)
               write(unit=cditem,fmt="(i2,a)") order,"a"
               !cdd=(/'a','b',cditem/)
               cdd(1)='a'
               cdd(2)='b'
               cdd(3)=cditem
             end if
           else
             order=codd(1)/codd(3)
             dif=abs(suma-real(order))
             if(dif < epss) then
               cod(1)=cod(3)
               multip(1)=suma
               moment(1)=suma*moment(3)
               write(unit=cditem,fmt="(i2,a)") order,"c"
               !cdd=(/cditem,'b','c'/)
               cdd(1)=cditem
               cdd(2)='b'
               cdd(3)='c'
             end if
            end if
         end if
         !Comparison b,c
         if(codd(2) /= 0 .and. codd(3)/=0) then
           suma=real(codd(2))/real(codd(3))
           if(abs(suma) < 1.0) then
             suma=1.0/suma
             order=codd(3)/codd(2)
             dif=abs(suma-real(order))
             if(dif < epss) then
               cod(3)=cod(2)
               multip(3)=suma
               moment(3)=suma*moment(2)
               write(unit=cditem,fmt="(i2,a)") order,"b"
               !cdd=(/'a','b',cditem/)
               cdd(1)='a'
               cdd(2)='b'
               cdd(3)=cditem
             end if
           else
             order=codd(2)/codd(3)
             dif=abs(suma-real(order))
             if(dif < epss) then
               cod(2)=cod(3)
               multip(2)=suma
               moment(2)=suma*moment(3)
               write(unit=cditem,fmt="(i2,a)") order,"c"
               !cdd=(/'a',cditem,'c'/)
               cdd(1)='a'
               cdd(2)=cditem
               cdd(3)='c'
             end if
            end if
         end if

       end if

       do j=1,3
         if(abs(multi(j)) < epss .or. cdd(j) == '0' ) then
           codes(j) = 0.0_cp
         else if(multi(j) < 0) then
           codes(j) = sign(1.0_cp, multi(j))*(abs(cod(j))*10.0_cp + abs(multip(j)) )
         else
           codes(j) = sign(1.0_cp, multip(j))*(abs(cod(j))*10.0_cp + abs(multip(j)) )
         end if
       end do

       if(present(Ipr)) then
         write(Ipr,'(a,3f10.4)')        '     Codes on Moments     : ',codes
         Write(Ipr,'(a,3(a,1x),6f7.3)') '     Codes and multipliers: ',cdd,multip
         Write(Ipr,'(a,3f12.4)')        '     Moment_TOT vector    : ',mome
       end if
       return
    End Subroutine get_moment_ctr_old


    !!----
    !!---- Subroutine Set_Shubnikov_Group(shubk,SG,MGp)
    !!----    character (len=*),         intent (in)    :: Shubk
    !!----    type(Magnetic_Group_Type), intent (out)   :: SG
    !!----    type(MagSymm_k_Type),      intent (in out):: MGp
    !!----
    !!----  This subroutined is not completed and it is now OBSOLETE because the
    !!----  availability of magnetic space groups databases
    !!----
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
    !!---- Updated: November 2006, June 2014
    !!
    Subroutine Write_Magnetic_Structure(Ipr,MGp,Am,Mag_Dom,cell)
       !---- Arguments ----!
       Integer,                    intent(in)           :: Ipr
       type(MagSymm_k_Type),       intent(in)           :: MGp
       type(mAtom_List_Type),      intent(in)           :: Am
       type(Magnetic_Domain_Type), intent(in), optional :: Mag_Dom
       type(Crystal_Cell_type),    intent(in), optional :: cell

       !---- Local Variables ----!
       character (len=100), dimension( 4):: texto
       character (len=80)                :: aux
       integer :: i,j,k,l, nlines,n,m,mult,nt
       real(kind=cp)                  :: x
       complex                        :: ci
       real(kind=cp), dimension(3)    :: xp,xo,u_vect,Mom,v
       real(kind=cp), dimension(3,3)  :: chi,chit
       real(kind=cp), dimension(3,48) :: orb
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
       if (MGp%Num_lat > 1) then
          write(unit=ipr,fmt="(a,i3)")  " => Centring vectors:",MGp%Num_lat-1
          nlines=1
          texto(:) (1:100) = " "
          do i=2,MGp%Num_lat
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

            if(Am%suscept) then
               Write(unit=ipr,fmt="(a,f8.3,a)")  &
               "  The magnetic structure is induced by an applied magnetic field of ",Am%MagField," Tesla"
               Write(unit=ipr,fmt="(a,3f8.3,a)")  &
               "  The direction of the applied magnetic field is: [",Am%dir_MField," ] in crystal space"
               do i=1,Am%Natoms
                  Write(unit=ipr,fmt="(a,a,5f10.5)")  &
                    "   Atom "//Am%Atom(i)%Lab, Am%Atom(i)%SfacSymb, Am%Atom(i)%x,Am%Atom(i)%Biso,Am%Atom(i)%occ
                  Write(unit=ipr,fmt="(a,6f10.5,a)")  &
                        "     Chi-Tensor( Chi11,Chi22,Chi33,Chi12,Chi13,Chi23) =  (", Am%Atom(i)%chi(:),")"
               end do

            else

               Write(unit=ipr,fmt="(a)")  &
               "  The Fourier coefficients are of the form: Sk(j) = 1/2 { Rk(j) + i Ik(j) } exp {-2pi i Mphask(j)}"
               Write(unit=ipr,fmt="(a)")  &
               "  They are written for each atom j as Sk( j)= 1/2 {(Rx Ry Rz)+i( Ix Iy Iz)} exp{-2pi i Mphask} -> MagMatrix # imat"
               Write(unit=ipr,fmt="(a)")  "  In case of k=2H (H reciprocal lattice vector) Sk(j)= (Rx Ry Rz)"

               do i=1,Am%Natoms
                  Write(unit=ipr,fmt="(a,a,5f10.5)")  &
                    "   Atom "//Am%Atom(i)%Lab, Am%Atom(i)%SfacSymb, Am%Atom(i)%x,Am%Atom(i)%Biso,Am%Atom(i)%occ
                  do j=1,Am%Atom(i)%nvk
                     if (K_Equiv_Minus_K(MGp%kvec(:,j),MGp%latt)) then
                        Write(unit=ipr,fmt="(a,i2,a,3f11.5,a,i4)")  &
                        "     Sk(",j,") =  (", Am%Atom(i)%Skr(:,j),")  -> MagMatrix #", Am%Atom(i)%imat(j)
                     else
                        Write(unit=ipr,fmt="(a,i2,a,2(3f11.5,a),f9.5,a,i4)")  &
                        "     Sk(",j,") = 1/2 {(", Am%Atom(i)%Skr(:,j),") + i (",Am%Atom(i)%Ski(:,j),")}  exp { -2pi i ",&
                        Am%Atom(i)%MPhas(j),"}  -> MagMatrix #", Am%Atom(i)%imat(j)
                     end if
                  end do
               end do
            end if

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
                  aux="(a,i2,a,  f11.5,a,f9.5,a,i4)"
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
               mult=0
               orb=0.0
               SOps: do k=1,MGp%NumOps
                  xp=ApplySO(MGp%SymOp(k),xo)
                  do nt=1,mult
                    v=orb(:,nt)-xp(:)
                    if (Lattice_trans(v,MGp%latt)) cycle SOps
                  end do
                  mult=mult+1
                  orb(:,mult)=xp(:)
                  Write(unit=ipr,fmt="(a,i2,a,3f9.5)") " =>  Atom "//Am%Atom(i)%lab//"(",k,") :",xp
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
                     Write(unit=ipr,fmt="(a,i2,a,2(3f11.5,a),f9.5,a)")  &
                      "     Sk(",j,") = 1/2 {(", real(Sk),") + i (",aimag(Sk),")}"
                  end do
               end do  SOps !Ops
               Write(unit=ipr,fmt="(a)") "  "
            end do  !atoms

         else !MGp%nirreps == 0

            if(Am%suscept .and. present(cell)) then
                u_vect=Am%MagField * Am%dir_MField / Veclength(Cell%Cr_Orth_cel,Am%dir_MField)
                do i=1,Am%natoms
                  xo=Am%Atom(i)%x
                  xo=modulo_lat(xo)
                  chi=reshape((/am%atom(i)%chi(1),am%atom(i)%chi(4), am%atom(i)%chi(5), &
                                am%atom(i)%chi(4),am%atom(i)%chi(2), am%atom(i)%chi(6), &
                                am%atom(i)%chi(6),am%atom(i)%chi(6), am%atom(i)%chi(3) /),(/3,3/))
                  mult=0
                  orb=0.0
                  sym: do k=1,MGp%Numops
                     xp=ApplySO(MGp%SymOp(k),xo)
                     xp=modulo_lat(xp)
                     do nt=1,mult
                       v=orb(:,nt)-xp(:)
                       if (Lattice_trans(v,MGp%latt)) cycle sym
                     end do
                     mult=mult+1
                     orb(:,mult)=xp(:)
                     chit=matmul(MGp%SymOp(k)%Rot,chi)
                     chit=matmul(chit,transpose(MGp%SymOp(k)%Rot))
                     Mom=matmul(Chit,u_vect)

                     Write(unit=ipr,fmt="(a,i2,2(a,3f11.5),a)") " =>  Atom "//Am%Atom(i)%lab//"(",k,") :",xp, &
                                                                "   Induced moment: [",Mom," ]"
                     Write(unit=ipr,fmt="(a)")            "             Local Susceptibility Tensor: "
                     do j=1,3
                        Write(unit=ipr,fmt="(a,3f14.5)")  "                            ",chit(j,:)
                     end do
                  end do sym ! symmetry
                end do ! Atoms

            else !suscept

              do i=1,Am%natoms
                 xo=Am%Atom(i)%x
                 mult=0
                 orb=0.0
                 Ops: do k=1,MGp%NumOps
                    xp=ApplySO(MGp%SymOp(k),xo)
                    do nt=1,mult
                      v=orb(:,nt)-xp(:)
                      if (Lattice_trans(v,MGp%latt)) cycle Ops
                    end do
                    mult=mult+1
                    orb(:,mult)=xp(:)
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
                 end do Ops
                 Write(unit=ipr,fmt="(a)") "  "
              end do  !atoms
            end if !suscept
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

