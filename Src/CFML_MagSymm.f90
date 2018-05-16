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
!!----       CLEANUP_SYMMETRY_OPERATORS
!!----       GET_MOMENT_CTR
!!--++       GET_MORBIT
!!--++       GET_MORBIT_MOMENT
!!----       GET_STABILIZERM
!!----       INIT_ERR_MAGSYM
!!--++       INIT_MAGNETIC_SPACE_GROUP_TYPE
!!----       INIT_MAGSYMM_K_TYPE             !OZ made it public to use in Read_Refcodes_Magnetic_Structure
!!----       MAGNETIC_SPACE_GROUP_TYPE_TO_MAGSYMM_K_TYPE
!!----       MAGSYMM_K_TYPE_TO_MAGNETIC_SPACE_GROUP_TYPE
!!----       READN_SET_MAGNETIC_SPACE_GROUP
!!----       READN_SET_MAGNETIC_STRUCTURE
!!--++       READN_SET_MAGNETIC_STRUCTURE_CFL    [Overloaded]
!!--++       READN_SET_MAGNETIC_STRUCTURE_MCIF   [Overloaded]
!!----       SET_MAGNETIC_SPACE_GROUP
!!----       SET_SHUBNIKOV_GROUP
!!----       SETTING_CHANGE_MAGGROUP
!!----       WRITE_MAGNETIC_STRUCTURE
!!----       WRITE_MCIF
!!----       WRITE_SHUBNIKOV_GROUP
!!----
!!
 Module CFML_Magnetic_Symmetry

    !---- Use Modules ----!
    Use CFML_GlobalDeps,                only: cp, tpi,Write_Date_Time
    Use CFML_Math_General,              only: Trace, Zbelong, Modulo_Lat, equal_matrix,             &
                                              Equal_Vector,Sort
    Use CFML_Math_3D,                   only: Get_Cart_From_Spher,Determ_A, matrix_inverse, Veclength
    Use CFML_Symmetry_Tables,           only: ltr_a,ltr_b,ltr_c,ltr_i,ltr_r,ltr_f,Sys_cry,LATT
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Read_Xsym, Get_SymSymb,axes_rotation, &
                                              Sym_Oper_Type, Set_SpaceGroup,read_msymm, symmetry_symbol, &
                                              err_symm,err_symm_mess, set_SpG_Mult_Table,ApplySO,   &
                                              Lattice_Trans, Get_SO_from_Gener, Get_Centring_Vectors, &
                                              Get_Shubnikov_Operator_Symbol, MSym_Oper_Type, LatSym,&
                                              Magnetic_Space_Group_Type
    Use CFML_String_Utilities,          only: u_case, l_case, Frac_Trans_1Dig, Get_Separator_Pos,Pack_String, &
                                              Frac_Trans_2Dig, Get_Mat_From_Symb, getnum_std, Err_String,     &
                                              Err_String_Mess,setnum_std, getword, Get_Transf,ucase, &
                                              Get_Symb_From_Mat
    Use CFML_IO_Formats,                only: file_list_type, File_To_FileList
    Use CFML_Atom_TypeDef,              only: Allocate_mAtom_list, mAtom_List_Type, Get_Atom_2nd_Tensor_Ctr
    Use CFML_Scattering_Chemical_Tables,only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
                                              Magnetic_Form
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
              Write_Shubnikov_Group, Init_MagSymm_k_Type, Write_MCIF, get_magnetic_form_factor, &
              Calc_Induced_Sk, Readn_Set_Magnetic_Space_Group,Cleanup_Symmetry_Operators, &
              Set_Magnetic_Space_Group, Get_mOrbit_mom, get_moment_ctr, get_stabilizerm, &
              Setting_Change_MagGroup

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
    !!----    crystallographic symmetry operator o not (1).
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

    !!---- Function is_Lattice_vec(V,Ltr,nlat,nl) Result(Lattice_Transl)
    !!----    !---- Argument ----!
    !!----    real(kind=cp), dimension(3),   intent( in) :: v
    !!----    real(kind=cp), dimension(:,:), intent( in) :: Ltr
    !!----    integer,                       intent( in) :: nlat
    !!----    integer,                       intent(out) :: nl
    !!----    logical                                    :: Lattice_Transl
    !!----
    !!----  Logical function that provides the value .true. if the vector V is a
    !!----  lattice vector.
    !!----
    !!----  Created: February 2014 (JRC)
    !!----
    Function is_Lattice_vec(V,Ltr,nlat,nl) Result(Lattice_Transl)
       !---- Argument ----!
       real(kind=cp), dimension(3),   intent( in) :: v
       real(kind=cp), dimension(:,:), intent( in) :: Ltr
       integer,                       intent( in) :: nlat
       integer,                       intent(out) :: nl
       logical                                    :: Lattice_Transl

       !---- Local variables ----!
       real(kind=cp)   , dimension(3) :: vec
       integer                        :: i

       Lattice_Transl=.false.
       nl=0

       if (Zbelong(v)) then       ! if v is an integral vector =>  v is a lattice vector
          Lattice_Transl=.true.
       else                       ! if not look for lattice type
          do i=1,nlat
            vec=Ltr(:,i)-v
            if (Zbelong(vec)) then
              Lattice_Transl=.true.
              nl=i
              exit
            end if
          end do
       end if
       return
    End Function is_Lattice_vec

    !!---- Function get_magnetic_form_factor(element) result(formf)
    !!----   character(len=*),intent(in) :: element
    !!----   character(len=6)            :: formf
    !!----
    !!----   Function to get the symbol for the magnetic scattering vector corresponding to the
    !!----   input symbol (element + valence state). Useful for transforming magCIF files to PCR.
    !!----
    !!----  Created: February 2014 (JRC)
    !!----
    Function get_magnetic_form_factor(element) result(formf)
      character(len=*),intent(in) :: element
      character(len=6)            :: formf
      ! Local variables
      logical :: is_re
      integer :: i,valence,ier
      character(len=6)   :: melem,aux
      integer, parameter :: n_re =12
      character(len=*), parameter, dimension(n_re) :: re=(/"ce","pr","nd","sm","eu","gd","tb","dy","ho","er","tm","yb"/)

      melem=l_case(element)
      is_re=.false.
      do i=1,n_re
        if(index(melem,re(i)) /= 0) then
          is_re=.true.
           exit
        end if
      end do
      if(is_re) then
        aux=melem(3:)
        i=index(aux,"+")
        if(i /= 0) then
          aux(i:i)=" "
          read(unit=aux,fmt=*,iostat=ier) valence
          if(ier /= 0) valence=3
        else
           valence=3
        end if
        write(unit=formf,fmt="(a,i1)") "J"//melem(1:2),valence
      else
        i=index(melem,"+")
        if(i /= 0) then
          melem(i:i)=" "
          aux=melem(i-1:i-1)
          read(unit=aux,fmt=*,iostat=ier) valence
          if(ier /= 0) valence=3
          melem(i-1:i-1)=" "
        else
           valence=2
        end if
        write(unit=formf,fmt="(a,i1)") "M"//trim(melem),valence
      end if
      formf=u_case(formf)
      return
    End Function get_magnetic_form_factor

    Function Get_MagMatSymb(line,mcif) result(mxmymz_op)
      character(len=*), intent(in) :: line
      logical,          intent(in) :: mcif
      character(len=len(line))     :: mxmymz_op
      !--- Local variables ---!
      integer :: i
      mxmymz_op=" "
      if(mcif) then
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
      else
         do i=1,len_trim(line)
           Select Case(line(i:i))
             case("x")
                mxmymz_op=trim(mxmymz_op)//"u"
             case("y")
                mxmymz_op=trim(mxmymz_op)//"v"
             case("z")
                mxmymz_op=trim(mxmymz_op)//"w"
             case default
                mxmymz_op=trim(mxmymz_op)//line(i:i)
           End Select
         end do
      end if
    End Function Get_MagMatSymb

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

    !!---- Subroutine Cleanup_Symmetry_Operators(MSpG)
    !!----   Type(Magnetic_Space_Group_Type), intent(in out) :: MSpG
    !!----
    !!----  Subroutine to re-organize symmetry operators extracting lattice translations
    !!----  and anti-translations and reordering the whole set of operators.
    !!----  (Still under development). It is supposed that the identity symmetry operator
    !!----  is provided in the input MSpG object, otherwise ok=.false. and
    !!----  no re-order is done.
    !!----
    !!----  Created: February 2014 (JRC), July 2016 (JRC)
    !!----
    Subroutine Cleanup_Symmetry_Operators(MSpG)
      Type(Magnetic_Space_Group_Type), intent(in out) :: MSpG
      !--- Local variables ---!
      integer,      dimension(    MSpG%Multip) :: ip,inp
      real,         dimension(    MSpG%Multip) :: tr
      logical,      dimension(    MSpG%Multip) :: nul
      real(kind=cp),dimension(3,192)           :: Lat_tr
      real(kind=cp),dimension(3,192)           :: aLat_tr
      integer :: i,j,k,kp,L,m, Ng,num_lat, num_alat,invt,nl,i_centre,centri
      integer,    dimension(3,3) :: identity, nulo, inver,mat,imat
      real(kind=cp),dimension(3) :: v
      character(len=80)          :: ShOp_symb !
      character(len=4)           :: ired !
      logical                    :: centrosymm
      Type(MSym_Oper_Type),dimension(MSpG%Multip) :: MSymOp
      Type(Sym_Oper_Type), dimension(MSpG%Multip) :: SymOp
      character(len=80),   dimension(MSpG%Multip) :: SymbSymOp,SymbMSymOp
      character (len=*),dimension(0:2), parameter  :: Centro = &
                                         (/"Centric (-1 not at origin)", &
                                           "Acentric                  ", &
                                           "Centric (-1 at origin)    "/)

      identity=0; nulo=0
      do i=1,3
        identity(i,i)=1
      end do
      inver=-identity
      num_lat=0; num_alat=0
      centrosymm=.false.
      nul=.false.
      MSpG%MagType=1
      centri=1 !Default value for non-centrosymmetric groups or for those having
               !the centre of symmetry out of the origin.

      !The code below is to re-order the symmetry operators provided in the input MSpG object
      !----Start re-ordering
      do i=1,MSpG%Multip
        tr(i)=sum(abs(MSpG%SymOp(i)%tr))
      end do
      ip=0
      call sort(tr,MSpG%Multip,ip)
      do i=1,MSpG%Multip
        j=ip(i)
        SymOp(i) = MSpG%SymOp(j)
        MSymOp(i)= MSpG%MSymOp(j)
        SymbSymOp(i)=MSpG%SymOpSymb(j)
        SymbMSymOp(i)=MSpG%MSymOpSymb(j)
      end do
      MSpG%SymOp(:)=SymOp(:)
      MSpG%MSymOp(:)=MSymOp(:)
      MSpG%SymOpSymb(:) = SymbSymOp(:)
      MSpG%MSymOpSymb(:)= SymbMSymOp(:)

      !Reorder again the operators in case the identity is not given as the first operator
      j=0
      imat=MSpG%SymOp(1)%Rot(:,:)
      if(.not. ( equal_matrix(imat,identity,3) .and.  sum(abs(MSpG%SymOp(1)%tr))  < 0.0001)) then
        do i=2,MSpG%Multip
          imat=MSpG%SymOp(i)%Rot(:,:)
          if(equal_matrix(imat,identity,3) .and. sum(abs(MSpG%SymOp(i)%tr))  < 0.0001) then
            j=i
            exit
          end if
        end do
        if(j == 0) then
          err_MagSym=.true.
          err_MagSym_Mess="The identity operator is not provided in the mCIF file"
          return
        end if
        MSpG%SymOp(j)=MSpG%SymOp(1)
        MSpG%MSymOp(j)=MSpG%MSymOp(1)
        MSpG%SymOpSymb(j)=MSpG%SymOpSymb(1)
        MSpG%MSymOpSymb(j)=MSpG%MSymOpSymb(1)
        MSpG%SymOp(1)%Rot=identity
        MSpG%SymOp(1)%tr=0.0
        MSpG%MSymOp(1)%Rot=identity
        MSpG%MSymOp(1)%phas=1.0
        MSpG%SymOpSymb(1)="x,y,z"
        MSpG%MSymOpSymb(1)="mx,my,mz"
      end if

      !Now look for centre of symmetry associated with time inversion and promote
      !it to the second position
      j=0
      do i=2,MSpG%Multip
        imat=MSpG%SymOp(i)%Rot(:,:)
        if(equal_matrix(imat,inver,3) .and. sum(abs(MSpG%SymOp(i)%tr))  < 0.0001 .and. MSpG%MSymOp(i)%phas < 0.0) then
          j=i
          exit
        end if
      end do
      if(j /= 0) then
        MSpG%SymOp(j)=MSpG%SymOp(2)
        MSpG%MSymOp(j)=MSpG%MSymOp(2)
        MSpG%SymOpSymb(j)=MSpG%SymOpSymb(2)
        MSpG%MSymOpSymb(j)=MSpG%MSymOpSymb(2)
        MSpG%SymOp(2)%Rot=inver
        MSpG%SymOp(2)%tr=0.0
        MSpG%MSymOp(2)%Rot=inver
        MSpG%MSymOp(2)%phas=-1.0
        MSpG%SymOpSymb(2)="-x,-y,-z"
        MSpG%MSymOpSymb(2)="-mx,-my,-mz"
      end if

      !----End re-ordering

      err_MagSym=.false.
      ip=0

      !Determine the lattice translations and anti-translations
      !Eliminate lattice translations and anti-translations
      do j=2,MSpG%Multip
        if(nul(j)) cycle
        invt= nint(MSpG%MSymOp(j)%phas)
        if(invt < 0) MSpG%MagType=3
        if(equal_matrix(identity,MSpG%SymOp(j)%Rot(:,:),3)) then
           if(invt == 1) then
              num_lat=num_lat+1
              Lat_tr(:,num_lat)=MSpG%SymOp(j)%tr(:)
              nul(j)=.true.   !Nullify centring translations
           else
              num_alat=num_alat+1
              aLat_tr(:,num_alat)=MSpG%SymOp(j)%tr(:)
              nul(j)=.true.  !Nullify anti-centring translations
           end if
        end if
      end do  !j=2,MSpG%Multip

      if(allocated(MSpG%Latt_trans)) deallocate(MSpG%Latt_trans)
      allocate(MSpG%Latt_trans(3,num_lat+1))
      MSpG%Latt_trans=0.0
      m=1
      do j=1,num_lat
        m=m+1
        MSpG%Latt_trans(:,m) = Lat_tr(:,j)
      end do
      MSpG%Num_Lat=num_lat+1

      if(num_alat > 0) then
        MSpG%MagType=4
        if(allocated(MSpG%aLatt_trans)) deallocate(MSpG%aLatt_trans)
        allocate(MSpG%aLatt_trans(3,num_alat))
        MSpG%aLatt_trans = aLat_tr(:,1:num_alat)
        MSpG%Num_aLat=num_alat
      end if

      !Eliminate centre of symmetry {-1|t} and select that having
      !t=0 if it exist
      k=0; kp=0
      do j=2,MSpG%Multip
          invt= nint(MSpG%MSymOp(j)%phas)
          imat=MSpG%SymOp(j)%Rot(:,:)
          if(equal_matrix(imat,inver,3)) then
            if(invt == 1) then
              kp=kp+1
              ip(kp)=j
            else
              k=k+1
              inp(k)=j
            end if
          end if
      end do

      i_centre=0
      if(kp > 0) then  !Centre of symmetry exist!, select that without translations
         i_centre=ip(1)
         do j=1,kp
           i=ip(j)
           if(sum(abs(MSpG%SymOp(i)%tr))  < 0.0001) then
             i_centre=i    !localization of the -x,-y,-z,+1 operator within the list
             centri=2      !Now this value indicates that the operor -x,-y,-z,+1 exists
             centrosymm=.true.
             nul(i)=.true.
             exit
           end if
         end do
      end if

      !Nullify the operators of inversion centres associated with time inversion
      !and have a translation corresponding to a centring or anticentring vector

      do i=1,k
         j=inp(i)
         v=MSpG%SymOp(j)%tr(:)
         if(sum(abs(v)) < 0.0001) cycle !Maintain the operaror -x,-y,-z,-1
         if(is_Lattice_vec(V,Lat_tr,num_lat,nl)) then
            nul(j)=.true.
            cycle
         end if

         if(is_Lattice_vec(V,aLat_tr,num_alat,nl)) then
            nul(j)=.true.
            cycle
         end if
      end do

      !Nullify the operators that can be deduced from others by applying translations,
      !anti-translations and centre of symmetry

      ip=0
      do j=2,MSpG%Multip-1
         do i=j+1,MSpG%Multip
           if(nul(i)) cycle
           mat=MSpG%SymOp(i)%Rot(:,:)-MSpG%SymOp(j)%Rot(:,:)
           if(equal_matrix(mat,nulo,3) ) then  !Pure lattice translation or antitranslation
              v=MSpG%SymOp(i)%tr(:)-MSpG%SymOp(j)%tr(:)

              if(is_Lattice_vec(V,Lat_tr,num_lat,nl)) then
                 nul(i)=.true.
                 cycle
              end if

              if(is_Lattice_vec(V,aLat_tr,num_alat,nl)) then
                 nul(i)=.true.
                 cycle
              end if

           end if

           if(centrosymm) then
              imat=MSpG%SymOp(i)%Rot(:,:)+MSpG%SymOp(j)%Rot(:,:)
              k=nint(MSpG%MSymOp(i)%phas)
              invt=nint(MSpG%MSymOp(j)%phas)

              if(equal_matrix(imat,nulo,3) .and. k == invt) then
                 v=MSpG%SymOp(i_centre)%tr(:)-MSpG%SymOp(i)%tr(:)-MSpG%SymOp(j)%tr(:)
                 if(is_Lattice_vec(V,Lat_tr,num_lat,nl)) then
                    nul(i)=.true.
                    cycle
                 end if
              end if

              if(equal_matrix(imat,nulo,3) .and. k /= invt) then
                 if(is_Lattice_vec(V,aLat_tr,num_alat,nl)) then
                    nul(i)=.true.
                    cycle
                 end if
              end if

           end if
         end do
      end do
      j=0

      ! => This is the reduced set of symmetry operators"
      do i=1,MSpG%Multip
        !write(*,"(a,i4,2x,L)") "  "//trim(MSpG%SymOpSymb(i))//"   "//trim(MSpG%MSymOpSymb(i)), nint(MSpG%MSymOp(i)%phas), nul(i)
        if(nul(i)) cycle
        j=j+1
        SymOp(j) = MSpG%SymOp(i)
        MSymOp(j)= MSpG%MSymOp(i)
      end do

      m=j*centri*(num_alat+num_lat+1)
      if( m /= MSpG%Multip) then
        write(unit=err_MagSym_Mess,fmt="(2(a,i4))") " Warning! Multip=",MSpG%Multip, " Calculated Multip: ",m
        err_MagSym=.true.
        return
      end if
      !Promote the reduced set of symmetry operators to the top of the list
      MSpG%SymOp(1:j)=SymOp(1:j)
      MSpG%MSymOp(1:j)=MSymOp(1:j)
      MSpG%Numops=j

      !Re-Construct, in an ordered way, all the symmetry operators in MSpG
      !starting with the reduced set
      m=MSpG%Numops

      if(centrosymm) then   !First apply the centre of symmetry
        MSpG%Centred=2
        MSpG%centre= centro(MSpG%Centred)
        do i=1,MSpG%Numops
          m=m+1
          MSpG%SymOp(m)%Rot  = -MSpG%SymOp(i)%Rot
          MSpG%SymOp(m)%tr   =  modulo_lat(-MSpG%SymOp(i)%tr)
          MSpG%MSymOp(m)= MSpG%MSymOp(i)
        end do
      else
        if(i_centre /= 0) then
          MSpG%Centred      = 0
          MSpG%centre       = centro(MSpG%Centred)
          MSpG%Centre_coord = MSpG%SymOp(i_centre)%tr(:)/2.0
        else
          MSpG%Centred=1
          MSpG%centre= centro(MSpG%Centred)
        end if
      end if

      ng=m

      if(MSpG%Num_aLat > 0) then   !Second apply the lattice centring anti-translations
        do L=1,MSpG%Num_aLat
           do i=1,ng
             m=m+1
             v=MSpG%SymOp(i)%tr(:) + MSpG%aLatt_trans(:,L)
             MSpG%SymOp(m)%Rot  = MSpG%SymOp(i)%Rot
             MSpG%SymOp(m)%tr   = modulo_lat(v)
             MSpG%MSymOp(m)%Rot = -MSpG%MSymOp(i)%Rot
             MSpG%MSymOp(m)%phas= -MSpG%MSymOp(i)%phas
           end do
        end do
      end if

      if(MSpG%Num_Lat > 1) then  !Third apply the lattice centring translations
        do L=2,MSpG%Num_Lat
           do i=1,ng
             m=m+1
             v=MSpG%SymOp(i)%tr(:) + MSpG%Latt_trans(:,L)
             MSpG%SymOp(m)%Rot  = MSpG%SymOp(i)%Rot
             MSpG%SymOp(m)%tr   = modulo_lat(v)
             MSpG%MSymOp(m)     = MSpG%MSymOp(i)
           end do
        end do
      end if


      !Normally here the number of operators should be equal to multiplicity
      !Test that everything is OK
      ng=m
      if(ng /= MSpG%Multip) then
        write(unit=err_MagSym_Mess,fmt="(2(a,i3))") " => Problem! the multiplicity ",MSpG%Multip," has not been recovered, value of ng=",ng
        err_MagSym=.true.
        return
      end if
      !Now re-generate all symbols from the symmetry operators and magnetic matrices
      ired=" => "
      do i=1,MSpG%Multip
         if(i > MSpG%Numops) ired="    "
         call Get_Shubnikov_Operator_Symbol(MSpG%SymOp(i)%Rot,MSpG%MSymOp(i)%Rot,MSpG%SymOp(i)%tr,ShOp_symb,.true.)
         j=index(ShOp_symb," ")
         MSpG%SymOpSymb(i)=ShOp_symb(1:j-1)
         ShOp_symb=adjustl(ShOp_symb(j:))
         j=index(ShOp_symb," ")
         MSpG%MSymOpSymb(i)=ShOp_symb(1:j-1)
         !write(*,"(a,i6,a,i4)") ired,i, "  "//trim(MSpG%SymOpSymb(i))//"   "//trim(MSpG%MSymOpSymb(i)), nint(MSpG%MSymOp(i)%phas)
      end do
      return
    End Subroutine Cleanup_Symmetry_Operators

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
       MGp%trn_to_parent=" "
       MGp%trn_from_parent=" "
       MGp%trn_to_standard=" "
       MGp%trn_from_standard=" "
       MGp%Multip=0
       MGp%n_wyck=0
       MGp%n_kv=0
       return
    End Subroutine Init_Magnetic_Space_Group_Type

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

    !!---- Subroutine Readn_Set_Magnetic_Space_Group(file_line,n_ini,n_end,MGp,mode,uvw)
    !!----    character(len=*),dimension(:),  intent (in)  :: file_line
    !!----    integer,                        intent (in)  :: n_ini,n_end
    !!----    type(Magnetic_Space_Group_Type),intent (out) :: MGp
    !!----    character(len=*),               intent (in)  :: mode
    !!----    character(len=*), optional,     intent (in)  :: uvw
    !!----
    !!----  This subroutine reads the lattice centring and anti-centring vectors
    !!----  as well as the symmetry operators of a magnetic space group in an
    !!----  arbitrary BNS setting. It construct the relevant magnetic space group
    !!----  components necessary for magnetic structure factor calculations.
    !!----  It may be used for reading from a CFL or a PCR file.
    !!----
    !!----  Created: March 2016.
    !!----
    !!----
    Subroutine Readn_Set_Magnetic_Space_Group(file_line,n_ini,n_end,MGp,mode,uvw)
       character(len=*),dimension(:),  intent (in)  :: file_line
       integer,                        intent (in)  :: n_ini,n_end
       type(Magnetic_Space_Group_Type),intent (out) :: MGp
       character(len=*),               intent (in)  :: mode
       character(len=*), optional,     intent (in)  :: uvw
       !
       ! --- Local variables ---!
       character(len=8)                 :: typ
       character(len=40)                :: Symbol
       integer                          :: i,j,ind,Nsym, Cen, N_Clat, N_Ant,ini, &
                                           num_sym,m,nop,k,n,L,ier,icount
       integer, dimension(3,3)          :: isim,msim
       real(kind=cp)                    :: p_mag
       real(kind=cp), dimension(3)      :: tr,v
       character(len=132)               :: line,ShOp_symb,setting,Parent
       character(len=40),dimension(10)  :: words
       logical                          :: u_type,m_type,inv_type, ttst,nonmag

       typ=l_case(adjustl(mode))

       call Init_Magnetic_Space_Group_Type(MGp)

       !Check if the database has to be read.
       nonmag=.false.; ttst=.false.
       if(typ /= "database") then
          do i=n_ini,n_end
           line=l_case(adjustl(file_line(i)))
           ind=index(line,"transform to standard:")
           if(ind /= 0) ttst=.true.
           ind=index(line,"<--nonmagnetic")
           if(ind /= 0) nonmag=.true.
           if(nonmag .and. ttst) then
             typ="database"
             exit
           end if
          end do
       end if

       Select Case(trim(typ))

          Case("pcr")
             line=adjustl(file_line(n_ini))
             ind=index(line,"Magnetic Space")
             if(ind == 0) then
               Err_Magsym=.true.
               Err_Magsym_Mess=" The Magnetic Space Group symbol is not provided in the PCR file! "
               return
             else
               j=index(line,"number:")
               MGp%BNS_symbol=trim(line(1:j-1))
               MGp%BNS_number=trim(line(j+7:ind-4))
             end if
             ini=n_ini+1
             Nsym=0; Cen=0; N_Clat=0;  N_Ant=0
             do i=ini,N_end
               line=adjustl(file_line(i))
               ind=index(line,"Transform to standard:")
               if(ind /= 0) then
                 MGp%trn_to_standard=adjustl(line(ind+22:))
               end if
               ind=index(line,"Parent Space Group:")
               if(ind /= 0) then
                 j=index(line,"IT_number:")
                 if( j /= 0) then
                   MGp%Parent_spg=adjustl(line(ind+19:j-1))
                   read(unit=line(j+10:),fmt=*,iostat=ier) MGp%Parent_num
                   if(ier /= 0) MGp%Parent_num=0
                 else
                   MGp%Parent_spg=adjustl(line(ind+19:))
                 end if
               end if
               ind=index(line,"Transform from Parent:")
               if(ind /= 0) then
                 MGp%trn_from_parent=adjustl(line(ind+22:))
               end if
               ind=index(line,"N_Clat")
               if(ind == 0) cycle
               read(unit=file_line(i+1),fmt=*) Nsym, Cen, N_Clat, N_Ant
               ini=i+2
               exit
             end do
             if(Nsym == 0) then
               Err_Magsym=.true.
               Err_Magsym_Mess=" The number of symmetry operators is not provided in the PCR file! "
               return
             end if
             !Allocate components of the magnetic space group
             MGp%Num_aLat=0
             allocate(MGp%Latt_trans(3,N_Clat+1))
             MGp%Latt_trans=0.0
             MGp%Num_Lat=N_Clat+1
             if(N_Ant > 0) then
               allocate(MGp%aLatt_trans(3,N_Ant))
               MGp%aLatt_trans=0.0
               MGp%Num_aLat=N_Ant
               MGp%MagType=4
             end if
             MGp%Numops = Nsym
             MGp%Centred= max(1,Cen)
             MGp%Multip = MGp%Numops * MGp%Centred * (MGp%Num_Lat + MGp%Num_aLat)
             num_sym=MGp%Multip
             allocate(Mgp%SymopSymb(num_sym))
             allocate(Mgp%Symop(num_sym))
             allocate(Mgp%MSymopSymb(num_sym))
             allocate(Mgp%MSymop(num_sym))
             if(N_Clat > 0) then
               do i=ini,N_end
                 line=adjustl(file_line(i))
                 ind=index(line,"Centring vectors")
                 if(ind == 0) cycle
                 ini=i+1
                 exit
               end do
               if(ind == 0) then
                 Err_Magsym=.true.
                 Err_Magsym_Mess=" 'Centring vectors' line is not provided in the PCR file! "
                 return
               end if
               m=1
               do i=ini,ini+N_Clat-1
                 m=m+1
                 read(unit=file_line(i),fmt=*) MGp%Latt_trans(:,m)
               end do
               ini=ini+N_Clat
             end if
             if(N_Ant > 0) then
               do i=ini,N_end
                 line=adjustl(file_line(i))
                 ind=index(line,"Anti-Centring vectors")
                 if(ind == 0) cycle
                 ini=i+1
                 exit
               end do
               if(ind == 0) then
                 Err_Magsym=.true.
                 Err_Magsym_Mess=" 'Anti-Centring vectors' line is not provided in the PCR file! "
                 return
               end if
               m=0
               do i=ini,ini+N_Ant-1
                 m=m+1
                 read(unit=file_line(i),fmt=*) MGp%aLatt_trans(:,m)
               end do
               ini=ini+N_Ant
             end if
             !Check the type of symmetry operators given
             do i=ini,N_end
                line=adjustl(file_line(i))
                if(line(1:1) == "!") cycle
                j=index(line,"!")
                if( j > 1) line=line(1:j-1)  !remove comments
                call Getword(line, words, icount)
                ! Icount=2 => SHSYM  x,-y,z+1/2,-1    <= This type
                ! Icount=3 => SHSYM  x,-y,z+1/2  -1   <= This type or these types => SHSYM x,-y,z+1/2  -u,v,-w  or SHSYM x,-y,z+1/2  -mx,my,-mz
                ! Icount=4 => SHSYM x,-y,z+1/2  -u,v,-w -1    <= This type or this type =>  SHSYM  x,-y,z+1/2  -mx,my,-mz  -1
                if( icount < 2 .or. icount > 4) then
                 Err_Magsym=.true.
                 Err_Magsym_Mess=" Error in Shubnikov operator: "//trim(line)
                 return
                end if
                u_type=(index(line,"u") /= 0) .or. (index(line,"U") /= 0)
                m_type=(index(line,"mx") /= 0) .or. (index(line,"MX") /= 0)
                if(.not. (u_type .or. m_type)) inv_type=.true.
                exit
             end do

             !Reading reduced set of symmetry operators
             m=0
             do i=ini,N_end
               line=adjustl(file_line(i))
               if(line(1:1) == "!") cycle
               j=index(line,"!")
               if( j > 1) line=line(1:j-1)  !remove comments
               j=index(line," ")
               line=adjustl(line(j:))
               m=m+1
               if(m > Nsym) exit
               call Getword(line, words, j)
               Select Case (icount)
                 Case(2)
                    j=index(line,",",back=.true.)
                    MGp%SymopSymb(m)=line(1:j-1)
                    read(unit=line(j+1:),fmt=*,iostat=ier) n
                    if(ier /= 0) then
                       err_magsym=.true.
                       ERR_MagSym_Mess=" Error reading the time inversion in line: "//trim(file_line(i))
                       return
                    else
                       MGp%MSymOp(m)%phas=real(n)
                    end if
                    !write(*,"(a,i3)") trim(MGp%SymopSymb(m)),n
                 Case(3)
                    MGp%SymopSymb(m)=words(1)
                    MGp%MSymopSymb(m)=words(2)  !u,v,w or mx,my,mz or +/-1

                 Case(4)
                    MGp%SymopSymb(m)=words(1)
                    MGp%MSymopSymb(m)=words(2)  !u,v,w or mx,my,mz
                    read(unit=words(3),fmt=*,iostat=ier) n
                    if(ier /= 0) then
                       err_magsym=.true.
                       ERR_MagSym_Mess=" Error reading the time inversion in line: "//trim(file_line(i))
                       return
                    else
                       MGp%MSymOp(m)%phas=real(n)
                    end if

               End Select
               call Read_Xsym(MGp%SymopSymb(m),1,isim,tr)
               MGp%Symop(m)%Rot=isim
               MGp%Symop(m)%tr=tr
               if(inv_type) then
                 j=determ_a(isim)
                 msim=nint(MGp%MSymOp(m)%phas)*j*isim
               else if (u_type) then
                 line=trim(MGp%MSymopSymb(m))//",0.0"
                 CALL read_msymm(line,msim,p_mag)
               else !should be mx,my,mz
                 line=trim(MGp%MSymopSymb(m))
                 do j=1,len_trim(line)
                    if(line(j:j) == "m" .or. line(j:j) == "M") line(j:j)=" "
                    if(line(j:j) == "x" .or. line(j:j) == "X") line(j:j)="u"
                    if(line(j:j) == "y" .or. line(j:j) == "Y") line(j:j)="v"
                    if(line(j:j) == "z" .or. line(j:j) == "Z") line(j:j)="w"
                 end do
                 line=pack_string(line)//",0.0"
                 CALL read_msymm(line,msim,p_mag)
               end if
               MGp%MSymop(m)%Rot=msim
               if(m_type .and. .not. present(uvw)) then
                 call Get_Shubnikov_Operator_Symbol(isim,msim,tr,ShOp_symb,.true.,invt=j)
                 MGp%mcif=.true.
               else
                 call Get_Shubnikov_Operator_Symbol(isim,msim,tr,ShOp_symb,invt=j)
                 MGp%mcif=.false.
               end if
               !write(*,"(a,i3)") trim(ShOp_symb),j
               MGp%MSymOp(m)%phas=j
               if(m_type .and. .not. present(uvw)) then
                 call Getword(ShOp_symb, words, j)
                 MGp%MSymopSymb(m)=words(2)
               else
                 j=index(ShOp_symb,";")
                 k=index(ShOp_symb,")")
                 MGp%MSymopSymb(m)=ShOp_symb(j+1:k-1)
               end if
             end do

          Case("cfl") !to be implemented
             write(unit=*,fmt="(a)") " => CFL file not yet implemented!"

          Case("database")
             line=adjustl(file_line(n_ini))
             ind=index(line,"Magnetic Space")
             if(ind == 0) then
               Err_Magsym=.true.
               Err_Magsym_Mess=" The Magnetic Space Group symbol is not provided in the PCR file! "
               return
             else
               j=index(line," ")
               symbol=trim(line(1:j-1))
             end if
             ini=n_ini+1
             line=adjustl(file_line(ini))
                       !     123456789012345678901234567890
             ind=index(line,"Transform to standard:")
             if(ind == 0) then
               Err_Magsym=.true.
               Err_Magsym_Mess=" The transformation to standard is needed even if it is: a,b,c;0,0,0 "
               return
             else
               ind=index(line,"<--")
               if( ind == 0) then
                 line=adjustl(line(23:))
                 ind=index(line," ")
                 setting=line(23:ind)
               else
                 setting=line(23:ind-1)
               end if
             end if
             ini=ini+1
             line=adjustl(file_line(ini))
             Parent=" "      !12345678901234567890
             ind= index(line,"Parent space group:")
             j  = index(line,"IT_number:")
             if(ind /= 0 .and. j /= 0) then
               Parent= adjustl(line(20:j-1))
               ind=index(line,"<--")
               Parent=trim(Parent)//" "//line(j+10:ind-1)
             end if
             ini=ini+1
             line=adjustl(file_line(ini))
             ind= index(line,"Transform from Parent:")
             if(ind /= 0) then
               j=index(line,"<--")
               Parent=trim(Parent)//"  "//line(23:j-1)
             end if
!C_ac  number: "9.41"                           <--Magnetic Space Group (BNS symbol and number)
!Transform to standard:  c,-b,a;0,0,0           <--Basis transformation from current setting to standard BNS
!Parent Space Group: Pna2_1  IT_number:   33    <--Non-magnetic Parent Group
!123456789012345678901234567890
!Transform from Parent:   a,2b,2c;0,0,0         <--Basis transformation from parent to current setting
             !write(*,"(a)") trim(symbol)//" "//trim(setting)//" "//trim(parent)
             ! trn_to=.true. always because magCIF considers thre transformation from the current
             ! setting to the standard setting
             if(len_trim(Parent) /= 0) then
               call Set_Magnetic_Space_Group(symbol,setting,MGp,parent,trn_to=.true.)
             else
               call Set_Magnetic_Space_Group(symbol,setting,MGp,trn_to=.true.)
             end if
             return       !The clean-up of operators is not needed
       End Select

       !Expand symmetry operators if Cen=2 (centre of symmetry at the origin)
       m=MGp%Numops
       if(Cen == 2) then
          do i=1,MGp%Numops
            m=m+1
            MGp%SymOp(m)%Rot  = -MGp%SymOp(i)%Rot
            MGp%SymOp(m)%tr   =  modulo_lat(-MGp%SymOp(i)%tr)
            MGp%MSymOp(m)%phas= MGp%MSymOp(i)%phas
            MGp%MSymOp(m)%Rot = MGp%MSymOp(i)%Rot
            call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
            call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
            !Expand the operator "line" to convert it to mx,my,mz like
            MGp%MSymopSymb(m)=Get_MagMatSymb(line,MGp%mcif)
          end do
       end if
       nop=m
       !Expand symmetry operators for lattice centrings
       do L=1,N_clat
         tr=MGp%Latt_trans(:,L+1)
         do j=1,nop
           m=m+1
           v=MGp%SymOp(j)%tr(:) + tr
           MGp%SymOp(m)%Rot  = MGp%SymOp(j)%Rot
           MGp%SymOp(m)%tr   = modulo_lat(v)
           MGp%MSymOp(m)%Rot = MGp%MSymOp(j)%Rot
           MGp%MSymOp(m)%phas= MGp%MSymOp(j)%phas
           call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
           call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
           !Expand the operator "line" to convert it to mx,my,mz like
           MGp%MSymopSymb(m)=Get_MagMatSymb(line,MGp%mcif)
         end do
       end do
       !Expand symmetry operators for lattice anti-centrings
       do L=1,N_Ant
         tr=MGp%aLatt_trans(:,L)
         do j=1,nop
           m=m+1
           v=MGp%SymOp(j)%tr(:) + tr
           MGp%SymOp(m)%Rot  = MGp%SymOp(j)%Rot
           MGp%SymOp(m)%tr   = modulo_lat(v)
           MGp%MSymOp(m)%Rot = -MGp%MSymOp(j)%Rot
           MGp%MSymOp(m)%phas= -MGp%MSymOp(j)%phas
           call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
           call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
           !Expand the operator "line" to convert it to mx,my,mz like
           MGp%MSymopSymb(m)=Get_MagMatSymb(line,MGp%mcif)
         end do
       end do
       ! Symmetry operators treatment done!

    End Subroutine Readn_Set_Magnetic_Space_Group

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
                                        irreps_begin, bfcoef_begin, magdom_begin, done, spg_given
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
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
                MGp%BNS_number=shubk
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_number_bns -> "//trim(MGp%BNS_number)

             Case("_magnetic_space_group_bns_name","_space_group_magn.name_bns","_space_group.magn_name_bns")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
                MGp%BNS_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_name_bns -> "//trim(MGp%BNS_symbol)

             Case("_magnetic_space_group_og_number","_space_group_magn.number_og","_space_group.magn_number_og")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
                MGp%OG_number=shubk
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_number_og -> "//trim(MGp%OG_number)

             Case("_magnetic_space_group_point_group","_space_group_magn.point_group","_space_group.magn_point_group")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
                MGp%PG_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group_magn.point_group -> "//trim(MGp%PG_symbol)

             Case("_magnetic_space_group_og_name","_space_group_magn.name_og","_space_group.magn_name_og")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
                MGp%OG_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_name_og -> "//trim(MGp%OG_symbol)

             Case("_magnetic_space_group.transform_from_parent_pp_abc","_magnetic_space_group_transform_from_parent_pp_abc", &
                   "_parent_space_group.child_transform_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
                MGp%trn_from_parent=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_transform_from_parent_pp_abc -> "//trim(MGp%trn_from_parent)

             Case("_parent_space_group.transform_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
                MGp%trn_to_parent=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_transform_from_parent_pp_abc -> "//trim(MGp%trn_from_parent)

             Case("_magnetic_space_group.transform_to_standard_pp_abc","_magnetic_space_group_transform_to_standard_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
                MGp%trn_to_standard=pack_string(shubk)

             Case("_space_group_magn.transform_bns_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=shubk(2:k-1)
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
       ext: do j=2,Spg%Multip
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

    Subroutine Get_mOrbit_mom(x,mom,Spg,Mult,orb,morb,ptr)
       !---- Arguments ----!
       real(kind=cp), dimension(3),    intent (in) :: x,mom
       type(Magnetic_Space_Group_type),intent (in) :: spg
       integer,                        intent(out) :: mult
       real(kind=cp),dimension(:,:),   intent(out) :: orb,morb
       integer,dimension(:),optional,  intent(out) :: ptr

       !---- Local variables ----!
       integer                                :: j, nt
       real(kind=cp), dimension(3)            :: xx,v,mmom,w
       character(len=1)                       :: laty

       laty="P"
       mult=1
       orb(:,1)=x(:)
       morb(:,1)=mom(:)
       if(present(ptr)) ptr(mult) = 1
       ext: do j=2,Spg%Multip
          xx=ApplySO(Spg%SymOp(j),x)
          xx=modulo_lat(xx)
          mmom=matmul(Spg%MSymOp(j)%Rot,mom)  !*Spg%MSymOp(j)%Phas <= This is already contained in Rot
          do nt=1,mult
             v=orb(:,nt)-xx(:)
             w=morb(:,nt)-mmom(:)
             if (Lattice_trans(v,laty)) then
              if(abs(sum(w)) > eps_symm) mmom=0.0
              cycle ext
             end if
          end do
          mult=mult+1
          orb(:,mult)=xx(:)
          morb(:,mult)=mmom(:)
          if(present(ptr)) ptr(mult) = j   !Effective symop
       end do ext
       return
    End Subroutine Get_mOrbit_mom

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
    Subroutine get_moment_ctr(xnr,moment,Spgr,codini,codes,ord,ss,att,Ipr)
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
         call get_stabilizerm(x,Spgr,order,ss_ptr,atr)
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
               write(unit=ipr,fmt='(a,i2,a,t20,a,t55,a,3f8.1)') '     Operator ',j,": ",trim(Spgr%SymopSymb(ss_ptr(j))), &
                trim(Spgr%MSymopSymb(ss_ptr(j))), Rsym
             end if
          end do
          mome=mome/real(order)
       end if
       do j=1,3
         if(abs(mome(j)-mom(j)) > epss) mome(j)=0.0
       end do
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
    End Subroutine get_moment_ctr

    !!---- Subroutine get_stabilizerm(x,Spg,order,ptr,atr)
    !!----    !---- Arguments ----!
    !!----    real(kind=cp), dimension(3),    intent (in)  :: x     ! real space position (fractional coordinates)
    !!----    type(Magnetic_Space_Group_type),intent (in)  :: Spg   ! Magnetic Space group
    !!----    integer,                        intent(out)  :: order ! Number of sym.op. keeping invariant the position x
    !!----    integer, dimension(:),          intent(out)  :: ptr   ! Array pointing to the symmetry operators numbers
    !!----                                                          ! of the stabilizer of x
    !!----    real(kind=cp), dimension(:,:),  intent(out)  :: atr   ! Associated additional translation to the symmetry operator
    !!----
    !!----
    !!----
    Subroutine get_stabilizerm(x,Spg,order,ptr,atr)
       !---- Arguments ----!
       real(kind=cp), dimension(3),    intent (in)  :: x     ! real space position (fractional coordinates)
       type(Magnetic_Space_Group_type),intent (in)  :: Spg   ! Space group
       integer,                        intent(out)  :: order ! Number of sym.op. keeping invariant the position x
       integer, dimension(:),          intent(out)  :: ptr   ! Array pointing to the symmetry operators numbers
                                                             ! of the stabilizer of x
       real(kind=cp), dimension(:,:),  intent(out)  :: atr   ! Associated additional translation to the symmetry operator
       !---- Local variables ----!
       real(kind=cp), dimension(3)    :: xx, tr

       integer                        :: j,n1,n2,n3

       order    = 1    !Identity belongs always to the stabilizer
       ptr(:)   = 0
       atr(:,:) = 0.0
       ptr(1)   = 1

       do n1=-1,1
        do n2=-1,1
          do n3=-1,1
            tr=real((/n1,n2,n3/))
             do j=2,Spg%multip
                xx=ApplySO(Spg%SymOp(j),x)+tr-x
                if (sum(abs(xx)) > 0.001) cycle
                order=order+1
                ptr(order)=j
                atr(:,order)=tr
             end do
          end do
        end do
       end do

       return
    End Subroutine get_stabilizerm

    !!----
    !!---- Subroutine Set_Magnetic_Space_Group(symb,setting,MSpg,parent,mcif,keepd,trn_to)
    !!----    character (len=*),                intent(in) :: symb        !  In -> String with the BNS symbol of the Shubnikov Group
    !!----    character (len=*),                intent(in ):: setting     !  In -> setting in the form -a,c,2b;1/2,0,0 (if empty no transformation is performed)
    !!----    Type (Magnetic_Space_Group_Type), intent(out):: MGp         ! Out -> Magnetic Space Group object
    !!----    character (len=*), optional,      intent(in ):: Parent      !  In -> Parent crystallographic group
    !!----    logical,  optional,               intent(in ):: mcif        !  In -> True if one wants to store the symbols as mx,my,mz
    !!----    logical,  optional,               intent(in ):: keepd       !  In -> True if one wants to keep the database allocated
    !!----    logical,  optional,               intent(in ):: trn_to      !  In -> True if the setting is from current TO standard setting
    !!----
    !!----    Subroutine constructing the object MGp from the BNS symbol by
    !!----    reading the database compiled by Harold T. Stokes and Branton J. Campbell
    !!----
    !!---- Created: November - 2016 (JRC)
    !!
    Subroutine Set_Magnetic_Space_Group(symb,setting,MSpg,parent,mcif,keepd,trn_to)
      character(len=*),               intent (in)  :: symb,setting
      type(Magnetic_Space_Group_Type),intent (out) :: MSpg
      character(len=*),optional,      intent (in)  :: parent
      logical,         optional,      intent (in)  :: mcif
      logical,         optional,      intent (in)  :: keepd
      logical,         optional,      intent (in)  :: trn_to
      !--- Local variables ---!
      integer                          :: i,j,m,k,n,L,ier,num,idem !,inv_time
      real(kind=cp)                    :: det
      !real(kind=cp), dimension(3)      :: orig
      real(kind=cp), dimension(3,3)    :: e !,S,Sinv
      integer, dimension(3,3)          :: identity
      character(len=256)               :: line,ShOp_symb
      logical                          :: change_setting,centring
      type(Magnetic_Space_Group_Type)  :: MGp

      call Init_Err_MagSym()
      call Allocate_DataBase()
      call read_magnetic_data()
      identity=0
      do i=1,3
        identity(i,i)=1
      end do
      e=identity
      !write(*,"(a)") trim(symb)//"  "//trim(setting)
      !if(present(parent)) write(*,"(a)") trim(Parent)
      !Check if the number of the magnetic group has been given
      !instead of the symbol
      read(unit=symb,fmt=*,iostat=ier) num
      if(ier /= 0) then
        num=0 !It is supposed that a symbol has been provided
        do i=1,magcount
          !write(*,"(i5,tr5,a)") i, spacegroup_label_bns(i)
          if(trim(symb) == trim(spacegroup_label_bns(i)) .or. &
             trim(symb) == trim(spacegroup_label_og(i))) then
            num=i
            exit
          end if
        end do
        if(num == 0) then
           write(unit=Err_MagSym_Mess,fmt="(a)") " => The BNS symbol: "//trim(symb)//" is illegal! "
           Err_MagSym=.true.
           if(.not. present(keepd)) call deAllocate_DataBase()
           return
        end if
      else
        if(num < 1 .or. num > magcount) then !magcount=1651
           write(unit=Err_MagSym_Mess,fmt="(a,i4,a)") " => The number of the Shubnikov group: ",num," is illegal!"
           Err_MagSym=.true.
           if(.not. present(keepd)) call deAllocate_DataBase()
           return
        end if
      end if
      if(len_trim(setting) == 0 .or. setting =='a,b,c;0,0,0') then
        change_setting=.false.
      else
        change_setting=.true.
      end if

      MGp%Sh_number=num
      MGp%BNS_number=nlabel_bns(num)
      MGp%OG_number= nlabel_og(num)
      MGp%BNS_symbol=spacegroup_label_bns(num)
      MGp%OG_symbol=spacegroup_label_og(num)
      MGp%MagType=magtype(num)
      !Setting the magnetic point group symbol from the BNS label
      m=0
      Select Case (MGp%MagType)
         Case(1,2,3)
           MGp%PG_Symbol=MGp%BNS_symbol(2:) !Remove the type of lattice
           do i=2,len_trim(MGp%BNS_symbol)
             m=m+1
             if(MGp%BNS_symbol(i:i) == "a" .or. MGp%BNS_symbol(i:i) == "b"  &
           .or. MGp%BNS_symbol(i:i) == "c" .or. MGp%BNS_symbol(i:i) == "d"  &
           .or. MGp%BNS_symbol(i:i) == "e" .or. MGp%BNS_symbol(i:i) == "g"  &
           .or. MGp%BNS_symbol(i:i) == "n") MGp%PG_Symbol(m:m)="m"
             if(MGp%BNS_symbol(i:i) == "_") MGp%PG_Symbol(m:m+1)=" "
           end do
           MGp%PG_Symbol=pack_string(MGp%PG_Symbol)

         Case(4)
           MGp%PG_Symbol=MGp%BNS_symbol(4:) !Remove the type of lattice
           do i=4,len_trim(MGp%BNS_symbol)
             m=m+1
             if(MGp%BNS_symbol(i:i) == "a" .or. MGp%BNS_symbol(i:i) == "b"  &
           .or. MGp%BNS_symbol(i:i) == "c" .or. MGp%BNS_symbol(i:i) == "d"  &
           .or. MGp%BNS_symbol(i:i) == "e" .or. MGp%BNS_symbol(i:i) == "g"  &
           .or. MGp%BNS_symbol(i:i) == "n") MGp%PG_Symbol(m:m)="m"
             if(MGp%BNS_symbol(i:i) == "_") MGp%PG_Symbol(m:m+1)=" "
           end do
           MGp%PG_Symbol=pack_string(MGp%PG_Symbol//"1'")
      End Select

      if(len_trim(setting) == 0 .or. setting =='a,b,c;0,0,0') then
        MGp%standard_setting=.true.
      else
        MGp%standard_setting=.false.
      end if
      MGp%mcif=.false.     !true if mx,my,mz notation is used , false is u,v,w notation is used
      if(present(mcif)) MGp%mcif=mcif
      MGp%m_cell=.true.    !true if magnetic cell is used for symmetry operators
      MGp%m_constr=.false. !true if constraints have been provided
      MGp%trn_from_parent=" "
      MGp%trn_to_standard="a,b,c;0,0,0"
      MGp%trn_from_standard="a,b,c;0,0,0"
      !Info about Parent Crystallographic Space Group
      if(present(parent)) then
        !Parent should be of the form  Xnnn  num  trn_from_parent
        line=adjustl(parent)
        i=index(line," ")
        MGp%Parent_spg=parent(1:i-1)
        line=adjustl(line(i:))
        i=index(line," ")
        read(unit=line(1:i),fmt=*,iostat=ier) MGp%Parent_num
        if(ier /= 0) then
           MGp%Parent_num=0
           MGp%trn_from_parent=line(1:i)
        else
           line=adjustl(line(i:))
           i=index(line," ")
           MGp%trn_from_parent=line(1:i-1)
        end if
      else
        !Try to deduce the parent space group from the BNS/OG numbers
        line=MGp%BNS_number
        i=index(line,".")
        line=line(1:i-1)
        read(unit=line,fmt=*) MGp%Parent_num
        if(MGp%MagType < 4) then
          MGp%Parent_spg=MGp%BNS_symbol
          if(MGp%MagType == 2) MGp%Parent_spg=MGp%Parent_spg(1:len_trim(MGp%Parent_spg)-2)
          do i=1,len_trim(MGp%Parent_spg)
            if(MGp%Parent_spg(i:i) == "'") MGp%Parent_spg(i:i) = " "
          end do
          MGp%Parent_spg=Pack_String(MGp%Parent_spg)
        else
          line=MGp%OG_number
          i=index(line,".")
          line=line(1:i-1)
          read(unit=line,fmt=*) MGp%Parent_num
          MGp%Parent_spg=MGp%OG_symbol
          if(MGp%Parent_spg(3:3) == "2") then
             MGp%Parent_spg(2:4)=" "
          else
             MGp%Parent_spg(2:3)=" "
          end if
          do i=1,len_trim(MGp%Parent_spg)
            if(MGp%Parent_spg(i:i) == "'") MGp%Parent_spg(i:i) = " "
          end do
          MGp%Parent_spg=Pack_String(MGp%Parent_spg)
        end if
      end if
      MGp%standard_setting = .true.
      ! Crystal system
      Select Case (num)
        case(1:7)
          MGp%CrystalSys="Triclinic"
        case(8:98)
          MGp%CrystalSys="Monoclinic"
        case(99:660)
          MGp%CrystalSys="Orthorhombic"
        case(661:1230)
          MGp%CrystalSys="Tetragonal"
        case(1231:1338)
          MGp%CrystalSys="Trigonal"
        case(1339:1502)
          MGp%CrystalSys="Hexagonal"
        case(1503:1651)
          MGp%CrystalSys="Cubic"
        case default
          MGp%CrystalSys="Unknown"
      End Select
      if(MGp%MagType == 4) then
        MGp%SPG_lat=spacegroup_label_bns(num)(1:3)
      else
        MGp%SPG_lat=spacegroup_label_bns(num)(1:1)
      end if
      MGp%SPG_latsy=MGp%SPG_lat !provisional before knowing s crystal system

      MGp%Num_Lat=lattice_bns_vectors_count(num)-2         ! Number of lattice points in a cell
      if(allocated(MGp%Latt_trans)) deallocate(MGp%Latt_trans)
      allocate(MGp%Latt_trans(3,MGp%Num_Lat))
      MGp%Latt_trans=0.0
      centring=.false.
      if(MGp%Num_Lat > 1) centring=.true.
      m=1
      do j=4,lattice_bns_vectors_count(num)
         m=m+1
         MGp%Latt_trans(:,m)= real(lattice_bns_vectors(:,j,num))/real(lattice_bns_vectors_denom(j,num))
      end do

      j=1
      MGp%Multip=wyckoff_mult(j,num)
      if(allocated(MGp%SymopSymb)) deallocate(MGp%SymopSymb)  ! Alphanumeric Symbols for SYMM
      if(allocated(MGp%SymOp))     deallocate(MGp%SymOp)      ! Crystallographic symmetry operators
      if(allocated(MGp%MSymopSymb))deallocate(MGp%MSymopSymb) ! Alphanumeric Symbols for MSYMM
      if(allocated(MGp%MSymOp))    deallocate(MGp%MSymOp)     ! Magnetic symmetry operators
      allocate(MGp%SymOp(MGp%Multip))
      allocate(MGp%SymopSymb(MGp%Multip))
      allocate(MGp%MSymOp(MGp%Multip))
      allocate(MGp%MSymopSymb(MGp%Multip))

      m=0
      !write(*,"(2(a,i5))") "Shubnikov number: ",num,"Wyckoff position count: ",wyckoff_pos_count(j,num)
      Do k=1,wyckoff_pos_count(j,num)
        idem=wyckoff_bns_fract_denom(k,j,num)
        MGp%SymOp(k)%tr=real(wyckoff_bns_fract(:,k,j,num))/real(idem)
        MGp%SymOp(k)%Rot = wyckoff_bns_xyz(:,:,k,j,num)
        MGp%MSymOp(k)%Rot = wyckoff_bns_mag(:,:,k,j,num)
        !inv_time=ops_bns_timeinv(k,num)  !Errors in the Database ... to be explored
        !MGp%MSymOp(k)%Phas=inv_time
        det=determ_a(MGp%SymOp(k)%Rot)
        if(det > 0.0) then
           if(equal_matrix(MGp%MSymOp(k)%Rot,MGp%SymOp(k)%Rot,3)) then
              MGp%MSymOp(k)%Phas=1.0
           else
              MGp%MSymOp(k)%Phas=-1.0
           end if
        else
           if(equal_matrix(MGp%MSymOp(k)%Rot,-MGp%SymOp(k)%Rot,3)) then
              MGp%MSymOp(k)%Phas=1.0
           else
              MGp%MSymOp(k)%Phas=-1.0
           end if
        end if
        if(MGp%mcif) then
           Call Get_Shubnikov_Operator_Symbol(MGp%SymOp(k)%Rot,MGp%MSymOp(k)%Rot,MGp%SymOp(k)%tr,ShOp_symb,MGp%mcif)
        else
           Call Get_Shubnikov_Operator_Symbol(MGp%SymOp(k)%Rot,MGp%MSymOp(k)%Rot,MGp%SymOp(k)%tr,ShOp_symb)
        end if
        !write(*,"(a)") trim(ShOp_symb)
        i=index(ShOp_symb,";")
        MGp%SymopSymb(k)=ShOp_symb(2:i-1)
        MGp%MSymopSymb(k)=ShOp_symb(i+1:len_trim(ShOp_symb)-1)
        if(MGp%MagType == 2) cycle
        if(equal_matrix(MGp%SymOp(k)%Rot,identity,3) .and. MGp%MSymOp(k)%Phas < 0.0) m=m+1 !counting anti-translations
      End Do

      if(centring) then
        n=wyckoff_pos_count(j,num)
        m=m*MGp%Num_Lat
        do L=2,MGp%Num_Lat
         do k=1,wyckoff_pos_count(j,num)
           MGp%SymOp(k+n)%Rot=MGp%SymOp(k)%Rot
           MGp%SymOp(k+n)%tr=Modulo_Lat(MGp%SymOp(k)%tr+MGp%Latt_trans(:,L))
           MGp%MSymOp(k+n)%Rot=MGp%MSymOp(k)%Rot
           MGp%MSymOp(k+n)%Phas=MGp%MSymOp(k)%Phas
           MGp%MSymopSymb(k+n)=MGp%MSymopSymb(k)
           Call Get_Shubnikov_Operator_Symbol(MGp%SymOp(k+n)%Rot,MGp%MSymOp(k+n)%Rot,MGp%SymOp(k+n)%tr,ShOp_symb)
           i=index(ShOp_symb,";")
           MGp%SymopSymb(k+n)=ShOp_symb(2:i-1)
         end do
         n=n+wyckoff_pos_count(j,num)
        end do
      end if

      MGp%Num_aLat=m       ! Number of anti-lattice points in a cell
      if(allocated(MGp%aLatt_trans)) deallocate(MGp%aLatt_trans)
      allocate(MGp%aLatt_trans(3,m))     ! Lattice anti-translations

      m=0
      if(MGp%MagType /= 2) then
        do k=1,MGp%multip
          if(equal_matrix(MGp%SymOp(k)%Rot,identity,3) .and. MGp%MSymOp(k)%Phas < 0) then
            m=m+1
            MGp%aLatt_trans(:,m) = MGp%SymOp(k)%tr
          end if
        end do
      end if
      MGp%Centred=0        ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
      MGp%Centre_coord=0.0 ! Fractional coordinates of the inversion centre
      do k=1,wyckoff_pos_count(j,num)
        if(equal_matrix(MGp%SymOp(k)%Rot,-identity,3) .and. MGp%MSymOp(k)%Phas > 0) then
          m=k
          MGp%Centred=max(MGp%Centred,1)
          if(sum(abs(MGp%SymOp(k)%tr)) < 0.001) then
            MGp%Centred=2
            exit
          end if
        end if
      end do
      MGp%NumOps=wyckoff_pos_count(j,num)
      MGp%Centre="Non-Centrosymmetric"       ! Alphanumeric information about the center of symmetry
      if(MGp%Centred == 1) then
        MGp%Centre="Centrosymmetric, -1 not @the origin "       ! Alphanumeric information about the center of symmetry
        MGp%Centre_coord=0.5*MGp%SymOp(m)%tr
      else if(MGp%Centred == 2) then
        MGp%Centre="Centrosymmetric, -1@the origin "       ! Alphanumeric information about the center of symmetry
        MGp%NumOps=MGp%NumOps/2
      end if
      if(change_setting) then
        if(present(trn_to)) then
          call Setting_Change_MagGroup(setting,MGp,MSpg,trn_to)
        else
          call Setting_Change_MagGroup(setting,MGp,MSpg)
        end if
        if(Err_MagSym) then
          if(.not. present(keepd)) call deAllocate_DataBase()
          return
        end if
      else
        MSpg=MGp !everything is allocated in the assignement (Fortran 2003)
      end if
      if(.not. present(keepd)) call deAllocate_DataBase()
    End Subroutine Set_Magnetic_Space_Group

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
    !!---- Subroutine Setting_Change_MagGroup(setting,MSpg,MSpgn,trn_to)
    !!----   character(len=*),                 intent(in) :: setting ! New origing in the old basis
    !!----   type (Magnetic_Space_Group_Type), intent(in) :: MSpG    ! Input space group
    !!----   type (Magnetic_Space_Group_Type), intent(out):: MSpGn   ! New space group in the new setting.
    !!----   logical, optional,                intent(in) :: trn_to  ! True if the input setting is "to standard"
    !!----                                                           ! If not present it is intrepreded as "from standard"
    !!----
    !!----    Transform the symmetry operators of the magnetic space group to
    !!----    a new basis given by the matrix "mat" and vector "orig", extracted
    !!----    from the string "setting" that is of the form:  a,b+c,c;1/2,0,0
    !!----
    !!---- Created: November - 2016 (JRC)
    !!
    Subroutine Setting_Change_MagGroup(setting,MSpg,MSpgn,trn_to)
       !---- Arguments ----!
       character(len=*),                 intent(in) :: setting
       type (Magnetic_Space_Group_Type), intent(in) :: MSpG
       type (Magnetic_Space_Group_Type),intent(out) :: MSpGn
       logical, optional,                intent(in) :: trn_to
       !--- Local variables ---!
       integer                 :: i, j, k, L, im, m, ngm,n, &
                                  nalat,nlat
       real(kind=cp)           :: det
       !character(len=40)       :: transla
       character(len=1)        :: LatSymb
       real(kind=cp), dimension (3,3), parameter :: e = reshape ((/1.0,0.0,0.0,  &
                                                                   0.0,1.0,0.0,  &
                                                                   0.0,0.0,1.0/),(/3,3/))
       real(kind=cp), dimension (3,192):: newlat = 0.0,alat=0.0 !big enough number of centring tranlations
       real(kind=cp), dimension (3,3)  :: S, Sinv, rot, rotn, mat, Matinv  !S is the ITC matrix P.
       integer,       dimension (3,3)  :: identity
       real(kind=cp), dimension (  3)  :: tr, trn, v, orig, iorig
       logical                         :: lattl
       character(len=80)               :: symbtr
       !character(len=60),dimension(15) :: gen
       character(len=120)              :: isetting
       character(len=132)              :: ShOp_symb
       real(kind=cp), allocatable, dimension(:,:,:) :: sm
       real(kind=cp), allocatable, dimension(:,:)   :: tm
       integer,       allocatable, dimension(:)     :: inv_time

       call Init_Err_MagSym()
       identity=nint(e)
       if(len_trim(setting) == 0) then !Check the error in the calling program
         Err_MagSym=.true.
         Err_MagSym_Mess=" => The string containing the setting change is empty!"
         return
       else
         call Get_Transf(setting,mat,orig)
         det=determ_a(mat)
         if( det < 0.0) then
           Err_MagSym_Mess=" => The transformation matrix should have a positive determinant!"
           Err_MagSym=.true.
           return
         end if
         S=transpose(mat)
         call matrix_inverse(S,Sinv,i)
         if (i /= 0) then
            Err_MagSym=.true.
            Err_MagSym_Mess= " => Wrong setting! Inversion Matrix Failed in Set_Magnetic_Space_Group!"
            return
         end if
         !----  A'= M A,  origin O =>  X'=inv(Mt)(X-O)   O'=-inv(Mt)O
         !----  A=inv(M)A'             X= Mt X'+ O       O= - Mt O'
         call matrix_inverse(Mat,Matinv,i)
         iorig=-matmul(Sinv,Orig)
         !Invers transformation -> to standard
         call Frac_Trans_2Dig(iorig,symbtr)
         i=len_trim(symbtr)
         symbtr=symbtr(2:i-1)
         call Get_Symb_From_Mat(Matinv,isetting,(/"a","b","c"/))
         isetting=trim(isetting)//";"//trim(symbtr)
         if(present(trn_to)) then
            if(trn_to) then
              rot=S
              rotn=Sinv
              S=rotn
              Sinv=rot
              v=orig
              tr=iorig
              orig=tr
              iorig=v
              write(unit=MSpGn%trn_from_standard,fmt="(a,f8.4)") adjustl(trim(isetting)//" -> det: "),1.0/det
              write(unit=MSpGn%trn_to_standard,fmt="(a,f8.4)") adjustl(trim(setting)//" -> det: "),det
            else
              write(unit=MSpGn%trn_to_standard,fmt="(a,f8.4)") adjustl(trim(isetting)//" -> det: "),1.0/det
              write(unit=MSpGn%trn_from_standard,fmt="(a,f8.4)") adjustl(trim(setting)//" -> det: "),det
            end if
         else
            write(unit=MSpGn%trn_to_standard,fmt="(a,f8.4)") adjustl(trim(isetting)//" -> det: "),1.0/det
         end if
       end if

       L=0
       if (MSpG%Num_Lat > 1) then  !Original lattice is centered
          do i=2,MSpG%Num_Lat      !Transform the centring vectors to the new lattice
             v=Modulo_Lat(matmul(Sinv,MSpG%Latt_trans(:,i)))
             if (sum(v) < eps_symm) cycle
             L=L+1
             newlat(:,L)=v
          end do
       end if
       do i=1,3  !Test the basis vectors of the original setting
         rot(:,i)=Modulo_Lat(Sinv(:,i))
         if (sum(rot(:,i)) < eps_symm) cycle
         L=L+1
         newlat(:,L)=rot(:,i)
       end do

       if (det > 1 ) then  !The new lattice is centred
          im=nint(det)-1   !Determine the new lattice translations
          ngm=L+im
          doi: do i=0,im
             v(1) = i
             do j=0,im
                v(2) = j
                do k=0,im
                   v(3) = k
                   if (nint(sum(v)) == 0) cycle
                   tr=Modulo_Lat(matmul(Sinv,v))
                   if (sum(tr) < eps_symm) cycle
                   lattl =.true.
                   do m=1,L
                      if (sum(abs(tr-newlat(:,m))) < eps_symm) then
                         lattl =.false.
                         exit
                      end if
                   end do
                   if (lattl) then ! new lattice translation
                      L=L+1
                      newlat(:,L) = tr(:)
                      if (L == ngm) exit doi
                   end if
                end do !k
             end do !j
          end do doi !i
       end if
       !The new multiplicity is obtained by multiplying the old one by the determinant
       MSpGn%PG_Symbol=MSpG%PG_Symbol
       MSpGn%multip=nint(MSpG%multip*det)
       allocate(sm(3,3,MSpGn%multip),tm(3,MSpGn%multip),inv_time(MSpGn%multip))

       call get_centring_vectors(L,newlat,LatSymb)  !Complete the centring vectors
       !Now we have L centring translations
       call LatSym(LatSymb,L,newlat)  !provides the value of the global variable inlat: index of the type of lattice
       MSpGn%SPG_lat     = LatSymb
       MSpGn%SPG_latsy   = LatSymb
       nlat=L+1
       MSpGn%Num_Lat=nlat
       allocate(MSpGn%Latt_trans(3,nlat))
       MSpGn%Latt_trans=0.0
       MSpGn%Latt_trans(:,2:nlat)= newlat(:,1:L)

       !---- Change of symmetry operator under a change of basis and origin
       !----  A'= M A,  origin O =>  X'=inv(Mt)(X-O)   O'=-inv(Mt)O
       !----  A=inv(M)A'             X= Mt X'+ O       O= - Mt O'
       !----  Symmetry operator C = (R,T)  -> C' = (R',T')
       !----   R' = inv(Mt) R Mt                 ITC:    R'= inv(P) R P
       !----   T' = inv(Mt) (T -(E-R)O)                  T'= inv(P) (T-(E-R)O)
       sm=0.0
       tm=0.0
       inv_time=0
       sm(:,:,1)=MSpG%SymOp(1)%Rot
       tm(:,1)=MSpG%SymOp(1)%tr
       inv_time(1)=nint(MSpG%MSymOp(1)%Phas)
       n=1
       !Transforming all the previous operators to the new cell,
       !including the previous lattice centring
       do_i:do i=2,MSpG%Multip
          Rot=MSpG%SymOp(i)%rot
          Rotn=matmul(matmul(Sinv,Rot),S)
          tr=MSpG%SymOp(i)%tr
          trn=Modulo_Lat(matmul(Sinv,tr-matmul(e-Rot,orig)))
          L=nint(MSpG%MSymOp(i)%Phas)
          do k=n,1,-1
            if(equal_matrix(Rotn,sm(:,:,k),3) .and. equal_vector(tm(:,k),trn,3) .and. &
            inv_time(k) == L)  cycle do_i
          end do
          n=n+1
          sm(:,:,n)=Rotn
          tm(:,n)=trn
          inv_time(n)=L
       end do do_i

       !Now complete the total set of operators by adding the new found
       !lattice translations to the previous set.
       n=MSpG%Multip
       do L=MSpG%Num_Lat,nlat-1
         do_im: do i=1,MSpG%Multip
           trn=modulo_lat(tm(:,i)+newlat(:,L))
           Rotn=sm(:,:,i)
           im=inv_time(i)
           do k=n,1,-1
             if(equal_matrix(Rotn,sm(:,:,k),3) .and. equal_vector(tm(:,k),trn,3) .and. &
             inv_time(k) == im)  cycle do_im
           end do
           n=n+1
           sm(:,:,n)=Rotn
           tm(:,n)=trn
           inv_time(n)=im
         end do do_im
       end do
       if( n /= MSpGn%multip) then
         Err_MagSym=.true.
         Err_MagSym_Mess=" => Error! The total multiplicity has not been recovered"
         return
       end if
       !Now we have the full set of operators sm,tm,inv_time
       !Construct the new magnetic space group
       allocate(MSpGn%SymOp(MSpGn%multip), MSpGn%SymOpSymb(MSpGn%multip))
       allocate(MSpGn%MSymOp(MSpGn%multip), MSpGn%MSymOpSymb(MSpGn%multip))
       MSpGn%NumOps=MSpGn%Multip/MSpGn%Num_Lat
       nalat=0
       do i=1,MSpGn%multip
         MSpGn%SymOp(i)%Rot=nint(sm(:,:,i))
         MSpGn%SymOp(i)%tr=tm(:,i)
         im=nint(determ_a(sm(:,:,i)))*inv_time(i)
         MSpGn%MSymOp(i)%Rot=nint(sm(:,:,i))*im
         MSpGn%MSymOp(i)%Phas=inv_time(i)
         if(equal_matrix(MSpGn%SymOp(i)%Rot,identity,3) .and. inv_time(i)==-1) then
           nalat=nalat+1
           alat(:,nalat)=MSpGn%SymOp(i)%tr
         end if
       end do
       MSpGn%Num_aLat=nalat
       allocate(MSpGn%aLatt_trans(3,nalat))
       MSpGn%aLatt_trans=alat

       MSpGn%Sh_number        = MSpG%Sh_number
       MSpGn%BNS_number       = MSpG%BNS_number
       MSpGn%OG_number        = MSpG%OG_number
       MSpGn%BNS_symbol       = MSpG%BNS_symbol
       MSpGn%OG_symbol        = MSpG%OG_symbol
       MSpGn%MagType          = MSpG%MagType
       MSpGn%mcif             = MSpG%mcif
       MSpGn%Parent_num       = MSpG%Parent_num
       MSpGn%Parent_spg       = MSpG%Parent_spg
       MSpGn%standard_setting = .false.
       MSpGn%CrystalSys       = MSpG%CrystalSys
       MSpGn%Centred=0
       do k=1,MSpGn%multip
         if(equal_matrix(MSpGn%SymOp(k)%Rot,-identity,3) .and. MSpGn%MSymOp(k)%Phas > 0) then
           m=k
           MSpGn%Centred=max(MSpGn%Centred,1)
           if(sum(abs(MSpGn%SymOp(k)%tr)) < 0.001) then
             MSpGn%Centred=2
             exit
           end if
         end if
       end do
       MSpGn%NumOps=MSpG%NumOps
       MSpGn%Centre="Non-Centrosymmetric"       ! Alphanumeric information about the center of symmetry
       if(MSpGn%Centred == 1) then
         MSpGn%Centre="Centrosymmetric, -1 not @the origin "
         MSpGn%Centre_coord=0.5*MSpGn%SymOp(m)%tr
       else if(MSpGn%Centred == 2) then
         MSpGn%Centre="Centrosymmetric, -1@the origin "
         MSpGn%NumOps=MSpGn%NumOps/2
       end if

       if(MSpG%mcif) then
         do i=1,MSpGn%multip
            call Get_Shubnikov_Operator_Symbol(MSpGn%Symop(i)%Rot,  &
                             MSpGn%MSymop(i)%Rot,MSpGn%Symop(i)%tr, &
                             ShOp_symb,.true.)
            j=index(ShOp_symb,";")
            MSpGn%SymopSymb(i)=ShOp_symb(2:j-1)
            MSpGn%MSymopSymb(i)=ShOp_symb(j+1:len_trim(ShOp_symb)-1)
         end do
       else
         do i=1,MSpGn%multip
            call Get_Shubnikov_Operator_Symbol(MSpGn%Symop(i)%Rot,  &
                             MSpGn%MSymop(i)%Rot,MSpGn%Symop(i)%tr, &
                             ShOp_symb)
            j=index(ShOp_symb,";")
            MSpGn%SymopSymb(i)=ShOp_symb(2:j-1)
            MSpGn%MSymopSymb(i)=ShOp_symb(j+1:len_trim(ShOp_symb)-1)
         end do
       end if
       return
    End Subroutine Setting_Change_MagGroup

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

    Subroutine Write_MCIF(Ipr,mCell,MSGp,Am,Cell)
       Integer,                         intent(in)           :: Ipr
       type(Magnetic_Space_Group_Type), intent(in)           :: MSGp
       type(Crystal_Cell_Type),         intent(in)           :: mCell
       type(mAtom_List_Type),           intent(in)           :: Am
       type(Crystal_Cell_Type),optional,intent(in)           :: Cell
       !
       Character(len=132)             :: line
       character(len=80),dimension(6) :: text
       character(len=2)               :: invc
       real(kind=cp)                  :: occ,occ_std,uiso,uiso_std
       integer :: i,j

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
       write(unit=Ipr,fmt="(a)") "_space_group.magn_number_BNS  "//trim(MSGp%BNS_number)
       if(len_trim(MSGp%BNS_symbol) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_name_BNS  "'//trim(MSGp%BNS_symbol)//'"'
       if(len_trim(MSGp%OG_number) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_number_OG '//trim(MSGp%OG_number)
       if(len_trim(MSGp%OG_symbol) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_name_OG  "'//trim(MSGp%OG_symbol)//'"'
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
          write(unit=Ipr,fmt="(a)") "_cell_length_a    "//trim(text(1))
          write(unit=Ipr,fmt="(a)") "_cell_length_b    "//trim(text(2))
          write(unit=Ipr,fmt="(a)") "_cell_length_c    "//trim(text(3))
          write(unit=Ipr,fmt="(a)") "_cell_angle_alpha "//trim(text(4))
          write(unit=Ipr,fmt="(a)") "_cell_angle_beta  "//trim(text(5))
          write(unit=Ipr,fmt="(a)") "_cell_angle_gamma "//trim(text(6))
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
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.id"
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.xyz"
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.mxmymz"
       do i=1,MSGp%Multip            !New mCIF format
          write(unit=invc,fmt="(i2)") nint(MSgp%MSymop(i)%Phas)
          if(invc(1:1) == " ") invc(1:1)="+"
          write(unit=Ipr,fmt="(i3,a)") i," "//trim(MSgp%SymopSymb(i))//","//invc//" "//trim(MSgp%MSymopSymb(i))
       end do
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "loop_"
       write(unit=Ipr,fmt="(a)") "_atom_site_label"
       write(unit=Ipr,fmt="(a)") "_atom_site_type_symbol"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_x"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_y"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_z"
       write(unit=Ipr,fmt="(a)") "_atom_site_U_iso_or_equiv"
       write(unit=Ipr,fmt="(a)") "_atom_site_occupancy"
       write(unit=Ipr,fmt="(a)") "_atom_site_symmetry_multiplicity"
       write(unit=Ipr,fmt="(a)") "_atom_site_Wyckoff_label"
       line=" "
       do i=1,Am%natoms
          do j=1,3
            call setnum_std(Am%atom(i)%x(j),Am%atom(i)%x_std(j),text(j))
          end do
          occ=real(MSgp%Multip)/real(Am%atom(i)%Mult)*Am%atom(i)%occ
          occ_std=real(MSgp%Multip)/real(Am%atom(i)%Mult)*Am%atom(i)%occ_std
          call setnum_std(occ,occ_std,text(5))
          uiso=Am%atom(i)%biso/78.95683521
          uiso_std=Am%atom(i)%biso_std/78.95683521
          call setnum_std(uiso,uiso_std,text(4))
          write(unit=Ipr,fmt="(a6,a6,3a13,2a11,i4,a)") Am%Atom(i)%lab, Am%atom(i)%SfacSymb,(text(j),j=1,5),&
                                                       Am%atom(i)%Mult," "//Am%atom(i)%wyck
       end do
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "loop_"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_label"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_x"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_y"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_z"
       do i=1,Am%natoms
          !if(sum(abs(Am%Atom(i)%Skr(:,1))) < 0.0001) cycle
          if(Am%Atom(i)%moment < 0.01) cycle
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

