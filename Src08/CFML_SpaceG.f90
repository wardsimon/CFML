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
!!----
!!---- MODULE: CFML_Groups
!!----         Space Groups and their algebra
!!----
!!--.. Rational matrix of special type dimension (3+d+1,3+d+1). The matrix of the 
!!--.. symmetry operator is extended with a column containing the translation in 
!!--.. the 3+d space plus a zero 3+d+1 row and +1 at position (3+d+1,3+d+1).
!!--..
!!--.. In order to limit the operators to the factor group w.r.t. traslations, a
!!--.. modulo-1 is applied in the multiplication of two operators.
!!----
!!
Module CFML_SpaceG
    !---- Use Modules ----!
    Use CFML_GlobalDeps
    Use CFML_Rational
    Use CFML_Strings,    only: u_case, pack_string

    !---- Variables ----!
    implicit none

    private
    
    !---- List of public operators ----!
    public :: operator (*) 
    public :: operator (==) 
    
    !---- List of public subroutines ----!
    public :: Init_SpaceG, &
              Read_Magnetic_Data, Read_Magnetic_Binary, &
              Set_Conditions_NumOP_EPS, &
              Write_SpaceG_Info
              
    !---- Parameters ----!
    integer,          dimension(0:2), parameter :: CENT=[2,1,2]       ! Multiplier for calculating the total multiplicity
    character(len=1), dimension(10),  parameter :: XYZ=["x","y","z","t","u","v","w","p","q","r"]
    character(len=1), dimension(10),  parameter :: ABC=["a","b","c","d","e","f","g","h","i","j"]
    character(len=3), dimension(10),  parameter :: X1X2X3=["x1 ","x2 ","x3 ","x4 ","x5 ","x6 ","x7 ","x8 ","x9 ","x10"]

    integer,                          parameter :: MAGCOUNT=1651      ! Magnetic Groups
    !---- Types ----!
    
    Type, public :: Symm_Oper_Type
       integer                                     :: Time_Inv =1         ! Time inversion for Magnetic groups
       integer                                     :: Dt       =1         ! determinant of the submatrix (1:3,1:3), it should be 1 or -1
       type(rational), dimension(:,:), allocatable :: Mat                 ! Matrix operator
    End Type Symm_Oper_Type
    
    Type, public :: Group_Type
       integer                                         :: Multip =0       ! Multiplicity of the Group 
       integer                                         :: D=0             ! Dimension of operator matrices (2:2D, 3:D3,...)
       type(Symm_Oper_Type), dimension(:), allocatable :: Op              ! Symmetry operator
       character(len=80),    dimension(:), allocatable :: Symb_Op         ! Symbol operator
    End Type Group_Type
    
    Type, public, extends(Group_Type) :: SPG_Type
       integer                                    :: numspg = 0
       integer                                    :: numshu = 0
       integer                                    :: numops = 0
       integer                                    :: centred= 0           !0 Centric(-1 no at origin), =1 Acentric, =2 Centric(-1 at origin)
       integer                                    :: anticentred=0        !=0 Centric(-1 no at origin), =1 Acentric, =2 Centric(-1' at origin)
       integer                                    :: mag_type = 0
       integer                                    :: num_lat  = 0
       integer                                    :: num_alat = 0
       character(len=1)                           :: spg_lat  =" "
       character(len=1), dimension(2)             :: shu_lat  =" "
       character(len=:),              allocatable :: spg_symb
       character(len=:),              allocatable :: shu_symb
       character(len=:),              allocatable :: pg
       character(len=:),              allocatable :: laue
       character(len=:),              allocatable :: mat2std
       character(len=:),              allocatable :: mat2std_shu
       character(len=:),              allocatable :: generators_list
       type(rational),dimension(:),   allocatable :: centre_coord
       type(rational),dimension(:),   allocatable :: anticentre_coord
       type(rational),dimension(:,:), allocatable :: Lat_tr
       type(rational),dimension(:,:), allocatable :: aLat_tr
    End Type SPG_Type
    
    !---- Private Variables ----!
    integer                                     :: MaxNum_OP=2048     ! Maximum number of Operators
    character(len=120), dimension(230)          :: it_spg_gen=" "     ! Generator of space groups in the standard setting
    type(rational), dimension(:,:), allocatable :: Identity_Matrix    ! Identity matrix
    
    logical                                     :: Magnetic_DBase_allocated=.false.
    logical                                     :: mcif=.false.
    
    !> For the ith nonhexagonal point operator:
    Character(Len=8),  dimension(:), public, allocatable :: point_op_label   ! point_op_label(i): point operator symbol (from Litvin)
    Character(Len=10), dimension(:), public, allocatable :: point_op_xyz
    Integer,       dimension(:,:,:), public, allocatable :: point_op_matrix  ! point_op_matrix(i): point operator matrix
   
    !> For the ith hexagonal point operator:
    Character(Len=8),  dimension(:), public, allocatable :: point_op_hex_label  ! point_op_hex_label(i): point operator symbol (from Litvin)
    Character(Len=10), dimension(:), public, allocatable :: point_op_hex_xyz    ! point_op_hex_xyz(i): point operator in x,y,z notation
    Integer,       dimension(:,:,:), public, allocatable :: point_op_hex_matrix ! point_op_hex_matrix(i): point operator matrix

    !> For the ith magnetic space group
    Character(Len=12), dimension(  :), public, allocatable :: nlabel_bns           ! nlabel_bns(i): numerical label in BNS setting
    Integer,           dimension(:,:), public, allocatable :: nlabelparts_bns      ! nlabel_parts_bns(j,i): jth part of nlabel_bns
    Character(Len=14), dimension(  :), public, allocatable :: spacegroup_label_bns ! label_bns(i): group symbol
    Character(Len=12), dimension(  :), public, allocatable :: nlabel_og            ! nlabel_og(i): numerical label in OG setting
    Integer,           dimension(:,:), public, allocatable :: nlabelparts_og       ! nlabel_parts_og(j,i): jth part of nlabel_og
    Character(Len=14), dimension(  :), public, allocatable :: spacegroup_label_og  ! label_og(i): group symbol
    Integer,           dimension(  :), public, allocatable :: magtype              ! magtype(i): type of magnetic space group (1-4)
   
    !> BNS-OG transformation (if type-4)
    Integer,         dimension(:,:,:), public, allocatable :: bnsog_point_op     ! bnsog_point_op(j,k,i): 3x3 point operator part of transformation
    Integer,         dimension(:,  :), public, allocatable :: bnsog_origin       ! bnsog_origin(j,i): translation part of transformation
    Integer,         dimension(    :), public, allocatable :: bnsog_origin_denom ! bnsog_point_origin(i): common denominator
    Integer,         dimension(    :), public, allocatable :: ops_count          ! iops_count(i): number of point operators
    Integer,         dimension(    :), public, allocatable :: wyckoff_site_count ! wyckoff_count(i): number of wyckoff sites
    Integer,         dimension(: , :), public, allocatable :: wyckoff_pos_count  ! wyckoff_pos_count(j,i): number of positions in jth wyckoff site
    Integer,         dimension(: , :), public, allocatable :: wyckoff_mult       ! wyckoff_mult(j,i): multiplicity for jth wyckoff site
    Character(Len=1),dimension(: , :), public, allocatable :: wyckoff_label      ! wyckoff_label(j,i): symbol (a,b,c,...,z,alpha) for jth wyckoff site

    !> For BNS setting
    Integer, dimension(    :), public, allocatable :: lattice_bns_vectors_count  ! number of lattice vectors defining the lattice
    Integer, dimension(:,:,:), public, allocatable :: lattice_bns_vectors        ! (k,j,i): kth component of the jth lattice vector
    Integer, dimension(:,  :), public, allocatable :: lattice_bns_vectors_denom  !(j,i): common denominator
   
    !> For jth operator
    Integer, dimension(  :,:), public, allocatable :: ops_bns_point_op    ! ops_bns_point_op(j,i): point operator part
    Integer, dimension(:,:,:), public, allocatable :: ops_bns_trans       ! ops_bns_trans(k,j,i): kth component of translation part
    Integer, dimension(  :,:), public, allocatable :: ops_bns_trans_denom ! ops_bns_trans_denom(j,i): common denominator
    Integer, dimension(  :,:), public, allocatable :: ops_bns_timeinv     ! ops_bns_timeinv(j,i): 1=no time inversion, -1=time inversion
   
    !> For jth wyckoff site
    Integer, dimension(:,  :,:,:), public, allocatable :: wyckoff_bns_fract       ! wyckoff_bns_fract(k,j,i): kth component of fractional part of wyckoff position
    Integer, dimension(    :,:,:), public, allocatable :: wyckoff_bns_fract_denom ! wyckoff_bns_fract_denom(j,i): common denominator
    Integer, dimension(:,:,:,:,:), public, allocatable :: wyckoff_bns_xyz         ! wyckoff_bns_xyz(m,k,j,i): mth component to coeffcient of kth parameter (x,y,z)
    Integer, dimension(:,:,:,:,:), public, allocatable :: wyckoff_bns_mag  ! wyckoff_bns_mag(m,k,j,i): mth component to coeffcient of kth magnetic parameter (mx,my,mz)
   
    !> For OG setting (for type-4 groups)
    Integer, dimension(    :), public, allocatable :: lattice_og_vectors_count  ! lattice_og_vectors_count(i): number of lattice vectors defining the lattice
    Integer, dimension(:,:,:), public, allocatable :: lattice_og_vectors   ! lattice_og_vectors(k,j,i): kth component of the jth lattice vector
    Integer, dimension(:,  :), public, allocatable :: lattice_og_vectors_denom  ! lattice_og_vectors_denom(j,i): common denominator
   
    !> For jth operator
    Integer, dimension(  :,:), public, allocatable :: ops_og_point_op    ! ops_og_point_op(j,i): point operator part
    Integer, dimension(:,:,:), public, allocatable :: ops_og_trans       ! ops_og_trans(k,j,i): kth component of translation part
    Integer, dimension(  :,:), public, allocatable :: ops_og_trans_denom ! ops_og_trans_denom(j,i): common denominator
    Integer, dimension(  :,:), public, allocatable :: ops_og_timeinv     ! ops_og_timeinv(j,i): 1=no time inversion, -1=time inversion
   
    !> For jth wyckoff site
    Integer, dimension(:,  :,:,:), public, allocatable :: wyckoff_og_fract        ! wyckoff_og_fract(k,j,i): kth component of fractional part of wyckoff position
    Integer, dimension(    :,:,:), public, allocatable :: wyckoff_og_fract_denom  ! wyckoff_og_fract_denom(j,i): common denominator
    Integer, dimension(:,:,:,:,:), public, allocatable :: wyckoff_og_xyz          ! wyckoff_og_xyz(m,k,j,i): mth component to coefficient of kth parameter (x,y,z)
    Integer, dimension(:,:,:,:,:), public, allocatable :: wyckoff_og_mag          ! wyckoff_og_mag(m,k,j,i): mth component to coefficient of kth magnetic parameter (mx,my,mz)

    
    !-------------------!
    !---- Operators ----!
    !-------------------!
    
    Interface operator (*)
       module procedure Multiply_Symm_Oper
    End interface
    
    Interface operator (==)
       module procedure Equal_Symm_Oper
       module procedure Equal_Group
    End interface
    
    !------------------!
    !---- Overload ----!
    !------------------!
    
    Interface Group_Constructor
       module procedure SpaceG_Constructor_GenV
       module procedure SpaceG_Constructor_Str
    End Interface Group_Constructor
    
    Interface Get_Lattice_Type
       module procedure Get_Lattice_Type_L
       module procedure Get_Lattice_Type_from_Mat
    End Interface Get_Lattice_Type


    !------------------------!
    !---- Interface Zone ----!
    !------------------------!
    Interface
       Module Subroutine Allocate_Magnetic_DBase()
          !---- Arguments ----!
       End Subroutine Allocate_Magnetic_DBase
          
       Module Subroutine Allocate_Operators(D, NMax, Op)
          !---- Arguments ----!
          integer,                                         intent(in)     :: d       ! Dimension
          integer,                                         intent(in)     :: NMax    ! is the expected maximum number of operators
          type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
       End Subroutine Allocate_Operators
       
       Module Subroutine Allocate_SpaceG(D, Multip, Grp)
          !---- Arguments ----!
          integer,         intent(in)     :: d
          integer,         intent(in)     :: multip
          class(Spg_Type), intent(in out) :: Grp
       End Subroutine Allocate_SpaceG
       
       Module Subroutine Allocate_Symm_Op(d, Op)
          !---- Arguments ----!
          integer,              intent(in)    :: d
          type(Symm_Oper_Type), intent(inout) :: Op
       End Subroutine Allocate_Symm_Op
       
       Module Subroutine Check_Gener(Gen_in, Gen_out)
          !---- Arguments ----!
          character(len=*), dimension(:),              intent(in)  :: gen_in
          character(len=*), dimension(:), allocatable, intent(out) :: gen_out
       End Subroutine Check_Gener 
       
       Module Subroutine Deallocate_Magnetic_DBase()
          !---- Arguments ----!
       End Subroutine Deallocate_Magnetic_DBase
       
       Module Function Equal_Group(Gr1, Gr2) Result(info)
          !---- Arguments ----!
          type(Spg_Type), intent(in) :: Gr1,Gr2
          logical                    :: info
       End Function Equal_Group
       
       Module Function Equal_Symm_Oper(Op1, Op2) Result(info)
          !---- Arguments ----!
          type(Symm_Oper_Type), intent(in) :: Op1,Op2
          logical                          :: info  
       End Function Equal_Symm_Oper
       
       Module Subroutine Get_Cosets(G,H, cosets)
          !---- Arguments ----!
          type(Spg_Type),                     intent(in)  :: G  
          type(Spg_Type),                     intent(in)  :: H  
          integer, dimension(:), allocatable, intent(out) :: cosets
       End Subroutine Get_Cosets
       
       Module Function Get_Dimension_Gener(Symb) Result(d)
          !---- Arguments ----! 
          character(len=*), intent(in) :: Symb
          integer                      :: d
       End Function Get_Dimension_Gener 
       
       Module Subroutine Get_Gener_From_Str(StrGen, d, ngen, gen)
          !---- Arguments ----!
          character(len=*),                            intent(in)  :: StrGen
          integer,                                     intent(out) :: d
          integer,                                     intent(out) :: ngen
          character(len=*), dimension(:), allocatable, intent(out) :: gen
       End Subroutine Get_Gener_From_Str  
       
       Module Function Get_Lattice_Type_from_MAT(M) Result(lattyp)
          !---- Arguments ----!
          type(rational), dimension(3,3), intent(in)  :: M
          character(len=1)                            :: lattyp
       End Function Get_Lattice_Type_from_MAT
       
       Module Function Get_Lattice_Type_L(L, Latc) Result(lattyp)
          !---- Arguments ----!
          integer,                        intent( in) :: L
          type(rational), dimension(:,:), intent( in) :: Latc
          character(len=1)                            :: lattyp
       End Function Get_Lattice_Type_L
       
       Module Subroutine Get_Magnetic_Lattice_Type(G)
          !---- Arguments ----!
          type(spg_type), intent(in out) :: G
       End Subroutine Get_Magnetic_Lattice_Type
       
       Module Subroutine Get_Mat_From_Symb(Symb, Mat, Invt) 
          !---- Arguments ----!
          character(len=*),                intent(in)  :: Symb
          type(rational), dimension(:,:),  intent(out) :: Mat
          integer, optional,               intent(out) :: invt
       End Subroutine Get_Mat_From_Symb
       
       Module Subroutine Get_Multip_OP_Table(Op,Table)
          !---- Arguments ----!
          type(Symm_Oper_Type), dimension(:), intent(in) :: Op
          integer, dimension(:,:),allocatable,intent(out):: Table
       End Subroutine Get_Multip_OP_Table
       
       Module Function Get_Op_from_Symb(Symb) Result(Op)
          !---- Arguments ----!
          character(len=*),     intent(in) :: symb
          type(Symm_Oper_Type)             :: Op
       End Function Get_Op_from_Symb
       
       Module Subroutine Get_OPS_From_Gener(Ngen, Ops, Multip, Table)
          !---- Arguments ----!
          integer,                                        intent(in)     :: ngen
          type(Symm_Oper_Type), dimension(:),             intent(in out) :: Ops
          integer,                                        intent(out)    :: multip
          integer, dimension(:,:), allocatable, optional, intent(out)    :: table
       End Subroutine Get_OPS_From_Gener 
       
       Module Subroutine Get_SubGroups(SpG, SubG, nsg, indexg, point)
          !---- Arguments ----!
          type(Spg_Type),                   intent( in) :: SpG
          type(Spg_Type),dimension(:),      intent(out) :: SubG
          integer,                          intent(out) :: nsg
          integer,                 optional,intent(in)  :: indexg
          logical, dimension(:,:), optional,intent(out) :: point
       End Subroutine Get_SubGroups
       
       Module Subroutine Get_SubGroups_Subgen(SpG, SubG, nsg, indexg)
          !---- Arguments ----!
          type(Spg_Type),                   intent( in) :: SpG
          type(Spg_Type),dimension(:),      intent(out) :: SubG
          integer,                          intent(out) :: nsg
          integer,                 optional,intent(in)  :: indexg
       End Subroutine Get_SubGroups_Subgen
       
       Module Function Get_Symb_from_Mat(Mat, Strcode, Invt) Result(Symb)
          !---- Arguments ----!
          type(rational), dimension(:,:), intent(in) :: Mat
          character(len=*), optional,     intent(in) :: strcode
          integer,          optional,     intent(in) :: invt
          character(len=80)                          :: symb
       End Function Get_Symb_from_Mat
       
       Module Function Get_Symb_from_Op(Op, Strcode) Result(symb)
          !---- Arguments ----!
          type(Symm_Oper_Type),       intent(in) :: Op
          character(len=*), optional, intent(in) :: Strcode
          character(len=80)                      :: symb
       End Function Get_Symb_from_Op
       
       Module Subroutine SpaceG_Constructor_GenV(GenV, Spg, StrCode)
          !---- Arguments ----!
          character(len=*),dimension(:),intent(in)     :: GenV
          class(Spg_Type),              intent(in out) :: Spg
          character(len=*),optional,    intent(in)     :: StrCode
       End Subroutine SpaceG_Constructor_GenV 
       
       Module Subroutine SpaceG_Constructor_Str(ListGen, Spg, Strcode)
          !---- Arguments ----!
          character(len=*),           intent(in)     :: ListGen
          class(Spg_Type),            intent(in out) :: Spg
          character(len=*), optional, intent(in)     :: Strcode  
       End Subroutine SpaceG_Constructor_Str 
       
       Module Subroutine Init_SpaceG(Grp)
          !---- Arguments ----!
          class(Group_type),  intent(in out) :: Grp 
       End Subroutine Init_SpaceG  
       
       Module Function Is_Inversion_Centre(Op) Result(Info)
          !---- Arguments ----!
          type(Symm_Oper_Type),intent(in) :: Op
          logical                         :: info
       End Function Is_Inversion_Centre
       
       Module Function Is_Lattice_Centring(Op) Result(Info)
          !---- Arguments ----!
          type(Symm_Oper_Type), intent(in) :: Op
          logical                          :: info
       End Function Is_Lattice_Centring
       
       Module Function Is_Lattice_Vec(V,Ltr,Nlat) Result(Lattice)
          !---- Argument ----!
          type(rational), dimension(:),   intent( in) :: v
          type(rational), dimension(:,:), intent( in) :: Ltr
          integer,                        intent( in) :: nlat
          logical                                     :: Lattice
       End Function Is_Lattice_Vec
       !
       Module Pure Function Multiply_Symm_Oper(Op1, Op2) Result (Op3)
          !---- Arguments ----!
          type(Symm_Oper_Type), intent(in) :: Op1,Op2
          type(Symm_Oper_Type)             :: Op3  
       End Function Multiply_Symm_Oper  
       
       Module Subroutine Read_Magnetic_Binary()
          !---- Arguments ----!
       End Subroutine Read_Magnetic_Binary 
       
       Module Subroutine Read_Magnetic_Data()
          !---- Arguments ----!
       End Subroutine Read_Magnetic_Data   
       
       Module Subroutine Reduced_Translation(Mat)
          !---- Arguments ----!
          type(rational), dimension(:,:), intent(in out) :: Mat
       End Subroutine Reduced_Translation  
       
       Module Subroutine Reorder_Operators(Multip, Op, Centred, Centre_Coord, Anticentred, &
                                           Anticentre_Coord, Numops, Num_Lat, Num_Alat,    &
                                           Lat_Tr, Alat_Tr, Mag_Type)
          !---- Arguments ----!
          integer,                            intent(in)     :: multip
          type(Symm_Oper_Type), dimension(:), intent(in out) :: Op
          integer,                            intent(out)    :: num_lat
          integer,                            intent(out)    :: num_alat
          integer,                            intent(out)    :: Numops
          integer,                            intent(out)    :: centred
          integer,                            intent(out)    :: anticentred
          integer,                            intent(out)    :: mag_type
          type(rational),dimension(:,:),      intent(out)    :: Lat_tr
          type(rational),dimension(:,:),      intent(out)    :: aLat_tr
          type(rational),dimension(:),        intent(out)    :: centre_coord
          type(rational),dimension(:),        intent(out)    :: anticentre_coord
       End Subroutine Reorder_Operators
       
       Module Subroutine Set_Conditions_NumOP_EPS(maxop,epsg)
          !---- Arguments ----!
          integer,       optional, intent(in) :: maxop
          real(kind=cp), optional, intent(in) :: epsg
       End Subroutine Set_Conditions_NumOP_EPS 
       
       Module Subroutine Set_Identity_Matrix(D)
          !---- Arguments ----! 
          integer, intent(in) :: D    
       End Subroutine Set_Identity_Matrix  
       
       Module Subroutine Smallest_Integral_Vector(v)
          !---- Arguments ----!
          type(rational), dimension(:), intent(inout) :: v
       End Subroutine Smallest_Integral_Vector
       
       Module Subroutine Sort_Oper(N, Op, Cod)
          !---- Arguments ----! 
          integer,                            intent(in)    :: n
          type(Symm_Oper_Type) ,dimension(n), intent(inout) :: Op
          character(len=*),                   intent(in)    :: cod   
       End Subroutine Sort_Oper  
       
       Module Subroutine Write_SpaceG_Info(Grp,Lun)
          !---- Arguments ----!
          class(Spg_Type),    intent(in)   :: Grp
          integer, optional,  intent(in)   :: lun 
       End Subroutine Write_SpaceG_Info   
      
    End Interface

End Module CFML_SpaceG
