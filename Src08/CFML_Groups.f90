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
Module CFML_Groups
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
    public :: Initialize_Conditions_Group, Initialize_Group, &
              Get_Subgroups_Subgen, Group_Constructor

    !---- Parameters ----!
    integer,          dimension(0:2), parameter :: CENT=[2,1,2]       ! Multiplier for calculating the total multiplicity
    character(len=1), dimension(10),  parameter :: XYZ=["x","y","z","t","u","v","w","p","q","r"]
    character(len=1), dimension(10),  parameter :: ABC=["a","b","c","d","e","f","g","h","i","j"]
    character(len=3), dimension(10),  parameter :: X1X2X3=["x1 ","x2 ","x3 ","x4 ","x5 ","x6 ","x7 ","x8 ","x9 ","x10"]

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
       module procedure Group_Constructor_GenV
       module procedure Group_Constructor_Str
    End Interface Group_Constructor


    !------------------------!
    !---- Interface Zone ----!
    !------------------------!
    Interface
       Module Subroutine Allocate_Group(D, Multip, Grp)
          !---- Arguments ----!
          integer,         intent(in)     :: d
          integer,         intent(in)     :: multip
          class(Spg_Type), intent(in out) :: Grp
       End Subroutine Allocate_Group
       
       Module Subroutine Allocate_Operator(d, Op)
          !---- Arguments ----!
          integer,              intent(in)    :: d
          type(Symm_Oper_Type), intent(inout) :: Op
       End Subroutine Allocate_Operator
       
       Module Subroutine Allocate_Operators(D, Multip, Op)
          !---- Arguments ----!
          integer,                                         intent(in)     :: d       ! Dimension
          integer,                                         intent(in)     :: multip  ! is the expected maximum number of operators
          type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
       End Subroutine Allocate_Operators
       
       Module Subroutine Check_Generators(Gen_in, Gen_out)
          !---- Arguments ----!
          character(len=*), dimension(:),              intent(in)  :: gen_in
          character(len=*), dimension(:), allocatable, intent(out) :: gen_out
       End Subroutine Check_Generators 
       
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
          type(Spg_Type),                     intent(in)  :: G  ! Group G > H
          type(Spg_Type),                     intent(in)  :: H  ! Subgroup of G
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
       
       Module Subroutine Get_Group_From_Gener(Ngen, Op, Multip, Table)
          !---- Arguments ----!
          integer,                                        intent(in)     :: ngen
          type(Symm_Oper_Type), dimension(:),             intent(in out) :: Op
          integer,                                        intent(out)    :: multip
          integer, dimension(:,:), allocatable, optional, intent(out)    :: table
       End Subroutine Get_Group_From_Gener 
       
       Module Subroutine Get_Group_From_Table(N, G, Table, Ord)  
          !---- Arguments ----! 
          integer,                 intent(in)     :: n   !Number of initial elements in G
          integer, dimension(:),   intent(in out) :: G   !Vector containing the different elements of the groups
          integer, dimension(:,:), intent(in)     :: Table
          integer,                 intent(out)    :: ord !order or the final group
       End Subroutine Get_Group_From_Table
       
       Module Subroutine Get_Mat_From_Symb(Symb,Mat,Invt)
          !---- Arguments ----!
          character(len=*),                intent(in)  :: Symb
          type(rational), dimension(:,:),  intent(out) :: Mat
          integer, optional,               intent(out) :: invt
       End Subroutine Get_Mat_From_Symb
       
       Module Subroutine Get_Multiplication_Table(Op,Table)
          !---- Arguments ----!
          type(Symm_Oper_Type), dimension(:), intent(in) :: Op
          integer, dimension(:,:),allocatable,intent(out):: Table
       End Subroutine Get_Multiplication_Table
       
       Module Function Get_Oper_from_Symb(Symb) Result(Op)
          !---- Arguments ----!
          character(len=*),     intent(in) :: symb
          type(Symm_Oper_Type)             :: Op
       End Function Get_Oper_from_Symb
       
       Module Function Get_Symb_from_Mat(Mat, Strcode, Invt) Result(Symb)
          !---- Arguments ----!
          type(rational), dimension(:,:), intent(in) :: Mat
          character(len=*), optional,     intent(in) :: strcode
          integer,          optional,     intent(in) :: invt
          character(len=80)                          :: symb
       End Function Get_Symb_from_Mat
       
       Module Function Get_Symb_from_Oper(Op, Strcode) Result(symb)
          !---- Arguments ----!
          type(Symm_Oper_Type),       intent(in) :: Op
          character(len=*), optional, intent(in) :: Strcode
          character(len=80)                      :: symb
       End Function Get_Symb_from_Oper
       
       Module Subroutine Get_SubGroups(SpG,SubG,nsg,indexg,point)
          !---- Arguments ----!
          type(Spg_Type),                   intent( in) :: SpG
          type(Spg_Type),dimension(:),      intent(out) :: SubG
          integer,                          intent(out) :: nsg
          integer,                 optional,intent(in)  :: indexg
          logical, dimension(:,:), optional,intent(out) :: point
       End Subroutine Get_SubGroups
       
       Module Subroutine Get_SubGroups_Subgen(SpG,SubG,nsg,indexg)
          !---- Arguments ----!
          type(Spg_Type),                   intent( in) :: SpG
          type(Spg_Type),dimension(:),      intent(out) :: SubG
          integer,                          intent(out) :: nsg
          integer,                 optional,intent(in)  :: indexg
       End Subroutine Get_SubGroups_Subgen
       
       Module Subroutine Group_Constructor_GenV(GenV, Grp, StrCode)
          !---- Arguments ----!
          character(len=*),dimension(:),intent(in)     :: GenV
          class(Spg_Type),              intent(in out) :: Grp
          character(len=*),optional,    intent(in)     :: StrCode
       End Subroutine Group_Constructor_GenV 
       
       Module Subroutine Group_Constructor_Str(ListGen, Grp, Strcode)
          !---- Arguments ----!
          character(len=*),           intent(in)     :: ListGen
          class(Spg_Type),            intent(in out) :: Grp
          character(len=*), optional, intent(in)     :: Strcode  
       End Subroutine Group_Constructor_Str 
       
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
       
       Module Subroutine Initialize_Conditions_Group(maxop,epsg)
          !---- Arguments ----!
          integer,       optional, intent(in) :: maxop
          real(kind=cp), optional, intent(in) :: epsg
       End Subroutine Initialize_Conditions_Group  
       
       Module Subroutine Initialize_Group(Grp)
          !---- Arguments ----!
          class(Group_type),  intent(inout) :: Grp 
       End Subroutine Initialize_Group  
       
       Module Pure Function Multiply_Symm_Oper(Op1, Op2) Result (Op3)
          !---- Arguments ----!
          type(Symm_Oper_Type), intent(in) :: Op1,Op2
          type(Symm_Oper_Type)             :: Op3  
       End Function Multiply_Symm_Oper  
       
       Module Subroutine Reduced_Translation(Mat)
          !---- Arguments ----!
          type(rational), dimension(:,:), intent(inout) :: Mat
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
       
       Module Subroutine Set_Identity_Matrix(d)
          !---- Arguments ----! 
          integer, intent(in) :: D    
       End Subroutine Set_Identity_Matrix  
       
       Module Subroutine Sort_Oper(N, Op,Cod)
          !---- Arguments ----! 
          integer,                            intent(in)    :: n
          type(Symm_Oper_Type) ,dimension(n), intent(inout) :: Op
          character(len=*),                   intent(in)    :: cod   
       End Subroutine Sort_Oper  
       
       Module Subroutine Write_Group_Info(Grp,Lun)
          !---- Arguments ----!
          class(Spg_Type),    intent(in)   :: Grp
          integer, optional,  intent(in)   :: lun 
       End Subroutine Write_Group_Info   
      
    End Interface

 Contains
   
   !!----
   !!---- RDET
   !!----
   !!---- 20/04/19
   !!
   Pure Recursive Function Rdet(A) Result(Acc)
      !---- Arguments ----!
      type(rational), dimension(:,:), intent(in) :: a
      type(rational)                             :: acc
      
      !---- Local variables ----!
      type(rational), dimension(size(a,dim=1)-1, size(a,dim=1)-1) :: b
      type(rational) :: sgn
      integer        :: i, n
      
      n=size(a,dim=1)
      if (n == 1) then
         acc = a(1,1)
      else
         acc = 0_LI//1_LI
         sgn = 1_LI/1_LI
         do i=1,n
            b(:, :(i-1)) = a(2:, :i-1)
            b(:, i:) = a(2:, i+1:)
            acc = acc + sgn * a(1, i) * rdet(b)
            sgn = sgn * (-1_LI/1_LI)
         end do
      end if
      
      return
   End Function rdet
   
   !!----
   !!---- LU_DESCOMPOSITION
   !!----
   !!---- 19/04/2019
   !!
   Pure Subroutine LU_Descomposition(a,p)
      !---- Arguments ----!
      real(kind=dp), intent(in out) :: a(:,:)
      integer,       intent(   out) :: p(:)
      
      !---- Local Variables ----!
      integer :: n,i,j,k,kmax
      
      n=size(a,1)
      p=[ ( i, i=1,n ) ]
      do k = 1,n-1
         kmax = maxloc(abs(a(p(k:),k)),1) + k-1
         if (kmax /= k ) p([k, kmax]) = p([kmax, k])
         a(p(k+1:),k) = a(p(k+1:),k) / a(p(k),k)
         forall (j=k+1:n) a(p(k+1:),j) = a(p(k+1:),j) - a(p(k+1:),k) * a(p(k),j)
      end do
      
      return
   End Subroutine LU_Descomposition
    
End Module CFML_Groups
