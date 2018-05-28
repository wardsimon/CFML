!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2018  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_Crystal_Metrics
!!----   INFO: Module to define crystallographic types and to provide
!!----         automatic crystallographic operations.
!!----
!!----
!!--.. INFORMATION
!!--..
!!--..    List Of Matrix Relationships For Crystallographic Applications
!!--..
!!--..    Small "t" is for transpose, inv(F) is the inverse of matrix F
!!--..
!!--..    Basis vectors as symbolic matrices
!!--..       At = (a,b,c)  At'=(a',b',c') ;  At* = (a*,b*,c*)  At*'=(a*',b*',c*')
!!--..
!!--..    Direct and reciprocal metric tensors: G, G*=inv(G)
!!--..    X  column vector in     direct space, referred to basis A
!!--..    X* column vector in reciprocal space, referred to basis A*
!!--..
!!--..       A'  = M  A           X'  = inv(Mt) X
!!--..       A*  = G* A           X*  =   G     X
!!--..       A*' = inv(Mt) A*     X*' =   M     X*
!!--..
!!--..       G' = M G Mt          G*' = inv(Mt) G* inv(M)
!!--..
!!--..   Symmetry operator defined in bases: A, A', A*, A*'
!!--..       C = (R,T), C'= (R',T'), C*= (R*,T*), C*'= (R*',T*')
!!--..
!!--..       R'  = inv(Mt) R Mt  ; T' = inv(Mt) T
!!--..       R*' =  M  R* inv(M) ; T*' = M T*
!!--..       R*  = G R G*  = inv(Rt)
!!--..
!!--..   If a change of origin is performed the positions are changed
!!--..   Ot=(o1,o2,o3) origin of the new basis A' w.r.t. old basis A
!!--..
!!--..       X' = inv(Mt) (X-O)
!!--..
!!--..   Changing just the origin   Xn  = C  X  = R  X  + T
!!--..                              Xn' = C' X' = R' X' + T'
!!--..          R=R'                X'  = X -O
!!--..                              Xn' = Xn-O
!!--..                  Xn-O = R' (X-O) + T' = R X + T - O
!!--..                   R X - R O + T' = R X + T - O
!!--..                               T' = T - (O - R O) = T - (E-R)O
!!--..
!!--..   Changing the basis (A,o) -> (A',o')
!!--..                  Xn  = C  X  = R  X  + T
!!--..                  Xn' = C' X' = R' X' + T'
!!--..                  X'= inv(Mt) (X-O), Xn' = inv(Mt) (Xn-O)
!!--..
!!--..            inv(Mt) (Xn-O) = R' inv(Mt) (X-O) + T'
!!--..            inv(Mt) (R  X  + T -O) = R' inv(Mt) (X-O) + T'
!!--..            inv(Mt) R X + inv(Mt)(T-O) = R' inv(Mt) X - R' inv(Mt) O + T'
!!--..            inv(Mt) R = R' inv(Mt)  => R' = inv(Mt) R Mt
!!--..            inv(Mt) (T-O)  = - R' inv(Mt) O + T'
!!--..            T' = R' inv(Mt) O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R Mt inv(Mt) O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R  O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R  O + inv(Mt) T - inv(Mt) O
!!--..            T' = inv(Mt)( R  O + T -  O) = inv(Mt) (T -(E-R)O)
!!--..
!!--..
!!--..                       R' = inv(Mt) R Mt
!!--..
!!--..                       T' = inv(Mt) (T -(E-R)O)
!!--..
!!--..
!!--..   A symmetry operator does not change the modulus of vectors and
!!--..   the angles between vectors (dot product is invariant):
!!--..
!!--..      X' = R X ,  Y' = R Y  =>  Xt' = Xt Rt,  Yt' = Yt Rt
!!--..
!!--..      Xt' G Y' = Xt Rt G R Y = Xt G Y  =>  G = Rt G R
!!--..
!!--..
!!--..   Second rank tensor Q and Q* defined in bases A and A*.
!!--..
!!--..      Q' = M Q Mt      Q* = G* Q G*     Q*'= inv(Mt) Q* inv(M)
!!--..
!!--..   A symmetry operator R is equivalent to a transformation
!!--..   M = inv(Rt) acting on basis vectors => G' = inv(Rt) G inv(R) = G
!!--..   The anisotropic temperature factors Beta is defined in reciprocal
!!--..   space: is a tensor like Q*, the transformation of beta under
!!--..   a symmetry operator is then :
!!--..
!!--..           Beta' = Inv(Mt) Beta inv(M) = R Beta Rt
!!--..
!!----
!!----
!!
 Module CFML_Crystal_Metrics

    !---- Use files ----!
    Use CFML_GlobalDeps,       only: CP, eps, PI, TO_RAD, err_CFML, clear_error
    Use CFML_Math_3D,          only: Invert_Array3x3, Determ_Vec, determ_3x3, Cross_Product
    Use CFML_Math_General,     only: co_linear, sort, co_prime, swap
    Use CFML_String_Utilities, only: U_Case

    implicit none

    private
    
    !---- Public Functions ----!
    public :: Volume_Cell, SigmaVolume, Set_Crystal_Cell
    public :: Write_Crystal_Cell, Write_Bin_Crystal_Cell, Read_Bin_Crystal_Cell

    
    !---- Definitions ----!
    
    !!----
    !!----  TYPE :: CRYSCELL_TYPE
    !!--..
    Type, public :: CrysCell_Type
       real(kind=cp),dimension(3)   :: cell  =0.0_cp      ! Cell parameters
       real(kind=cp),dimension(3)   :: ang   =0.0_cp      !
       real(kind=cp),dimension(3)   :: scell =0.0_cp      ! Standard deviation for cell paramaters
       real(kind=cp),dimension(3)   :: sang  =0.0_cp      !
       real(kind=cp)                :: vol   =0.0_cp      ! Volume and sig(V)
       real(kind=cp)                :: svol  =0.0_cp      !
    End Type CrysCell_Type
    
    !!----
    !!----  TYPE :: CRYSCELL_M_TYPE
    !!--..
    Type, public, extends(CrysCell_Type) :: CrysCell_M_Type
       real(kind=cp),dimension(3)   :: rcell  =0.0_cp      ! Reciprocal Cell parameters
       real(kind=cp),dimension(3)   :: rang   =0.0_cp      !
       real(kind=cp)                :: rvol   =0.0_cp
       real(kind=cp),dimension(3,3) :: GD     =0.0_cp      ! Direct Metric Tensor
       real(kind=cp),dimension(3,3) :: GR     =0.0_cp      ! Reciprocal Metric Tensor
       real(kind=cp),dimension(3,3) :: Cr_Orth_cel=0.0_cp  ! Fractional to Cartesian
       real(kind=cp),dimension(3,3) :: Orth_Cr_cel=0.0_cp  ! Cartesian to Fractional
       real(kind=cp),dimension(3,3) :: BL_M   =0.0_cp      ! Busing-Levy B-matrix 
       real(kind=cp),dimension(3,3) :: BL_Minv=0.0_cp      ! Inverse of Busing-Levy B-matrix
       character(len=1)             :: CartType=" "        ! if "A" x// a
    End Type CrysCell_M_Type
    
    !!----
    !!----  TYPE :: CRYSCELL_LS_TYPE
    !!--..
    Type, public, extends(CrysCell_Type) :: CrysCell_LS_Type
       integer, dimension(3) :: lcell=0   ! code for Refinements
       integer, dimension(3) :: lang =0   ! code for Refinements 
    End Type CrysCell_LS_Type
    
    !!----
    !!----  TYPE :: CRYSCELL_MLS_TYPE
    !!--..
    Type, public, extends(CrysCell_M_Type) :: CrysCell_MLS_Type
       integer, dimension(3) :: lcell=0   ! code for Refinements
       integer, dimension(3) :: lang =0   ! code for Refinements 
    End Type CrysCell_MLS_Type
    
    !!----
    !!----  TYPE :: TWOFOLD_AXES_TYPE
    !!--..
    !!
    Type, public :: Twofold_Axes_Type
       integer                        :: ntwo    =0                  ! Number of two-fold axes                           
       real(kind=cp)                  :: tol     =3.0_cp             ! Angular tolerance (ca 3 degrees)                  
       real(kind=cp) ,dimension(3,12) :: caxes   =0.0_cp             ! Cartesian components of two-fold axes             
       integer,dimension(3,12)        :: dtwofold=0                  ! Direct indices of two-fold axes                   
       integer,dimension(3,12)        :: rtwofold=0                  ! Reciprocal indices of two-fold axes               
       integer,dimension(12)          :: dot     =0                  ! Scalar product of reciprocal and direct indices   
       real(kind=cp), dimension(12)   :: cross   =0.0_cp             ! Angle between direct and reciprocal axes ( < tol) 
       real(kind=cp), dimension(12)   :: maxes   =0.0_cp             ! Modulus of the zone axes (two-fold axes) vectors  
       real(kind=cp), dimension(3)    :: a       =0.0_cp             ! Cartesian components of direct cell parameters    
       real(kind=cp), dimension(3)    :: b       =0.0_cp
       real(kind=cp), dimension(3)    :: c       =0.0_cp
    End Type Twofold_Axes_Type

    !!----
    !!----  TYPE :: ZONE_AXIS_TYPE
    !!--..
    !!
    Type, public :: Zone_Axis_Type
      Integer               :: nlayer =0   ! number of the reciprocal layer considered normally nlayer=0
      Integer, dimension(3) :: uvw    =0   ! Indices of the zone axis                                   
      Integer, dimension(3) :: rx     =0   ! Indices (reciprocal vector) of the basis vector 1          
      Integer, dimension(3) :: ry     =0   ! Indices (reciprocal vector) of the basis vector 2          
    End Type Zone_Axis_Type


    !> Parameters 
    real(kind=cp), parameter                 :: TPI2=2.0*PI*PI
    real(kind=cp), dimension(3,3), parameter :: IDENTITY= &
                   reshape ([1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0],[3,3])

    
    !---- Interfaces - Overloaded ----!
    Interface SigmaVolume
       Module Procedure SigmaV_CellType
       Module Procedure SigmaV_Cell
    End Interface SigmaVolume
    
    Interface  Niggli_Cell                   ! The first(s) argument(s) is(are)
      Module Procedure Niggli_Cell_abc       ! List of cell parameters passed as a 6D vector
      Module Procedure Niggli_Cell_nigglimat ! Niggli matrix passed as a 2x3 matrix (ultimately applying the algorithm)
      Module Procedure Niggli_Cell_Params    ! List of cell parameters a,b,c,alpha,beta,gamma
      Module Procedure Niggli_Cell_type      ! The object Cell is passed as argument
      Module Procedure Niggli_Cell_Vect      ! Input three vectors in Cartesian components
    End Interface  Niggli_Cell

    Interface
       Module Subroutine ReciprocalCell(cell,ang,rcell,rang,rVol)
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in ) :: cell,ang
          real(kind=cp), dimension(3), intent(out) :: rcell,rang
          real(kind=cp),               intent(out) :: rvol
       End Subroutine ReciprocalCell
       
       Module Subroutine Get_Cryst_Orthog_Matrix(Cell,Ang, Mat,CarType)
          !---- Arguments ----!
          real(kind=cp), dimension(3  ), intent (in ) :: cell,ang   ! Cell Parameters
          real(kind=cp), dimension(3,3), intent (out) :: Mat        ! Convsersion matrix
          character(len=*), optional,    intent (in ) :: CarType    ! Type of Cartesian axes  
       End Subroutine Get_Cryst_Orthog_Matrix   
       
       Module Function Metrics(cell,ang) Result(G)
          !---- Arguments ----!
          real(kind=cp), dimension(3)  , intent(in ) :: cell  ! Cell Parameters
          real(kind=cp), dimension(3)  , intent(in ) :: ang
          real(kind=cp), dimension(3,3)              :: G     ! Metric Tensor
       End Function Metrics   
       
       Module Pure Function Volume_Cell(cell,ang) Result(Vol)
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: cell
          real(kind=cp), dimension(3), intent(in) :: ang
          real(kind=cp)                           :: vol
       End Function Volume_Cell
       
       Module Pure Function SigmaV_CellType(Cell) Result(sigma)
          !---- Arguments ----!
          class(CrysCell_Type), intent(in) :: Cell      ! Cell Parameters
          real(kind=cp)                    :: sigma     ! Sigma
       End Function SigmaV_CellType
       
       Module Pure Function SigmaV_Cell(cell,ang,scell,sang) Result(sigma)
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: cell      ! Cell parameters
          real(kind=cp), dimension(3), intent(in) :: ang 
          real(kind=cp), dimension(3), intent(in) :: scell     ! standard deviation for cell parameters
          real(kind=cp), dimension(3), intent(in) :: sang
          real(kind=cp)                           :: sigma     ! Sigma
       End Function SigmaV_Cell 
       
       Module Subroutine Write_Crystal_Cell(Cell, Iunit)
          !---- Arguments ----!
          class(CrysCell_Type),  intent(in) :: Cell         ! Cell object
          Integer, optional,     intent(in) :: Iunit 
       End Subroutine Write_Crystal_Cell 
       
       Module Subroutine Write_Bin_Crystal_Cell(Cell,Iunit)
          !---- Arguments ----!
          class(CrysCell_Type),  intent(in) :: Cell       ! Cell object
          Integer,               intent(in) :: Iunit  
       End Subroutine Write_Bin_Crystal_Cell  
       
       Module Subroutine Read_Bin_Crystal_Cell(Cell,Iunit)
          !---- Arguments ----!
          class(CrysCell_Type),  intent(out) :: Cell       ! Cell object
          Integer,               intent(in)  :: Iunit  
       End Subroutine Read_Bin_Crystal_Cell 
       
       Module Subroutine Set_Crystal_Cell(VCell,VAng,Cell,Cartype,Vscell,Vsang)
          !---- Arguments ----!
          real(kind=cp), dimension(3),         intent(in ) :: Vcell, Vang    ! Cell parameters
          class(CrysCell_Type),                intent(out) :: Cell           ! Cell Object
          character (len=*),          optional,intent(in ) :: CarType        ! Orientation in Cartesian
          real(kind=cp), dimension(3),optional,intent(in ) :: Vscell, Vsang ! Standard deviations  
       End Subroutine Set_Crystal_Cell 
       
       Module Subroutine Get_Cryst_Family(Cell, Str_Family, Str_Symbol, Str_System)
          !---- Arguments ----!
          class(CrysCell_Type),   intent(in ) :: Cell
          character(len=*),       intent(out) :: Str_Family
          character(len=*),       intent(out) :: Str_Symbol
          character(len=*),       intent(out) :: Str_System 
       End Subroutine Get_Cryst_Family  
       
       Module Subroutine Get_Deriv_Orth_Cell(Cell,De_Orthcell,Cartype)
          !---- Arguments ----!
          class(CrysCell_type),            intent(in ) :: cell
          real(kind=cp), dimension(3,3,6), intent(out) :: de_Orthcell
          character (len=1), optional,     intent(in ) :: CarType  
       End Subroutine Get_Deriv_Orth_Cell 
       
       Module Subroutine Change_Setting_Cell(Cell,Mat,Celln,Matkind)
          !---- Arguments ----!
          class(CrysCell_Type),          intent( in)     :: Cell
          real(kind=cp), dimension (3,3),intent( in)     :: Mat
          class(CrysCell_Type),          intent(out)     :: Celln
          character(len=*),  optional,   intent (in)     :: Matkind  
       End Subroutine Change_Setting_Cell
       
       Module Subroutine Get_Primitive_Cell(Lat_Type,Centred_Cell,Primitive_Cell,Transfm)
          !---- Arguments ----!
          character(len=*),              intent(in)  :: lat_type          
          class(CrysCell_Type),          intent(in)  :: centred_cell      
          class(CrysCell_Type),          intent(out) :: primitive_cell    
          real(kind=cp), dimension(3,3), intent(out) :: transfm   
       End Subroutine Get_Primitive_Cell
       
       Module Subroutine Get_Transfm_Matrix(Cella,Cellb,Trm,tol)
          !---- Arguments ----!
          class(CrysCell_Type),          intent(in) :: cella  ! Cell object
          class(CrysCell_Type),          intent(in) :: cellb  ! Cell object
          real(kind=cp), dimension(3,3), intent(out):: trm    ! Transformation matrix
          real(kind=cp), optional,       intent(in) :: tol    ! Tolerance   
       End Subroutine Get_Transfm_Matrix 
       
       Module Subroutine Get_basis_from_uvw(dmin,u,cell,ZoneB,mode)
          !--- Arguments ---!
          real(kind=cp),              intent(in) :: dmin      ! Minimum d-spacing (smax=1/dmin)
          integer, dimension(3),      intent(in) :: u         ! Zone axis indices
          class(CrysCell_Type),       intent(in) :: cell      ! Cell object
          type (Zone_Axis_Type),      intent(out):: ZoneB     ! !Object containing u and basis vector in the plane
          character(len=*), optional, intent(in) :: mode         
       End Subroutine Get_basis_from_uvw
       
       Module Subroutine Get_TwoFold_Axes(Cell,Tol,Twofold)
          !---- Arguments ----!
          class(CrysCell_M_Type),   intent (in) :: Cell     ! Cell object
          real(kind=cp),           intent (in) :: Tol       ! angular tolerance in degrees
          Type(twofold_axes_type), intent(out) :: Twofold
       End Subroutine Get_TwoFold_Axes
       
       Module Subroutine Get_Conventional_Cell(Twofold,Cell,Tr,Message,told)
          !---- Arguments ----!
          Type(Twofold_Axes_Type), intent(in)  :: Twofold
          class(CrysCell_Type),    intent(out) :: Cell
          integer, dimension(3,3), intent(out) :: tr
          character(len=*),        intent(out) :: message
          real(kind=cp), optional, intent(in)  :: told  
       End Subroutine Get_Conventional_Cell 
       
       Module Function Cart_Vector(Mode,V,Cell) Result(Vc)
          !---- Arguments ----!
          character(len=*),            intent(in) :: mode      !  D: Direct, R: Reciprocal, BL or BLD 
          real(kind=cp), dimension(3), intent(in) :: v         !  Vector                   
          class(CrysCell_M_Type),      intent(in) :: Cell      !  Cell object           
          real(kind=cp), dimension(3)             :: vc        !                
       End Function Cart_Vector
       
       Module Function Cart_U_Vector(Mode,V,Cell) Result(Vc)
          !---- Arguments ----!
          character (len=*),           intent(in) :: Mode   ! Options, D, R, BL, BLD
          real(kind=cp), dimension(3), intent(in) :: v      ! Vector
          class(CrysCell_M_Type),      intent(in) :: Cell   ! Cell object
          real(kind=cp), dimension(3)             :: vc
       End Function Cart_U_Vector
       
       Module Function Rot_MetricalMatrix(V,Phi,Cell) Result(Mat)
          !---- Argument ----!
          real(kind=cp), dimension(3),      intent(in) :: V     ! Direction vector
          real(kind=cp),                    intent(in) :: phi   ! Degree of rotarion around V
          class(CrysCell_M_Type), optional, intent(in) :: cell  ! Cell object
          real(kind=cp), dimension(3,3)                :: Mat   ! Metrical Matrix rotated
       End Function Rot_MetricalMatrix  
       
       Module Pure Function U_Equiv(Cell, Th_U) Result(Uequi)
          !---- Arguments ----!
          class(CrysCell_M_Type),      intent(in)  :: Cell    ! Cell object
          real(kind=cp), dimension(6), intent(in)  :: Th_U    ! U thermal parameters
          real(kind=cp)                            :: Uequi   ! Uequiv
       End Function U_Equiv
       
       Module Pure Function Get_Betas_from_B(B,Cell) Result(Beta)
          !---- Arguments ----!
          real(kind=cp),dimension(6), intent(in)  :: B
          class(CrysCell_M_Type),     intent(in)  :: Cell
          real(kind=cp),dimension(6)              :: Beta
       End Function Get_Betas_from_B
       
       Module Pure Function Get_U_from_B(B) Result(U)
          !---- Arguments ----!
          real(kind=cp),dimension(6),  intent(in)  :: B
          real(kind=cp),dimension(6)               :: U
       End Function Get_U_from_B
       
       Module Pure Function Get_B_from_Betas(Beta,Cell) Result(B)
          !---- Arguments ----!
          real(kind=cp),dimension(6), intent(in)  :: Beta
          class(CrysCell_M_Type),     intent(in)  :: Cell
          real(kind=cp),dimension(6)              :: B   
       End Function Get_B_from_Betas
       
       Module Pure Function Get_Betas_from_U(U,Cell) Result(Beta)
          !---- Arguments ----!
          real(kind=cp),dimension(6),intent(in)  :: U
          class(CrysCell_M_Type),    intent(in)  :: Cell
          real(kind=cp),dimension(6)             :: Beta    
       End Function Get_Betas_from_U
       
       Module Pure Function Get_Betas_from_Biso(Biso,Cell) Result(Betas)
          !--- Argument ----!
          real(kind=cp),           intent(in)  :: Biso
          class(CrysCell_M_Type),  intent(in)  :: Cell
          real(kind=cp), dimension(6)          :: Betas   
       End Function Get_Betas_from_Biso
       
       Module Pure Function Get_U_from_Betas(Beta,Cell) Result(U)
          !---- Arguments ----!
          real(kind=cp),dimension(6),intent(in)  :: Beta
          class(CrysCell_M_Type),    intent(in)  :: Cell
          real(kind=cp),dimension(6)             :: U  
       End Function Get_U_from_Betas
       
       Module Pure Function Get_B_from_U(U) Result(B)
          !---- Arguments ----!
          real(kind=cp),dimension(6), intent(in)  :: U
          real(kind=cp),dimension(6)              :: B
       End Function Get_B_from_U  
       
       Module Subroutine Niggli_Cell_ABC(CellV,Niggli_Point,Cell,Trans)    
          !---- Arguments ----!
          real(kind=cp),dimension(6),             intent(in out) :: CellV         ! Cell parameters in a vector
          real(kind=cp),dimension(5),   optional, intent(out)    :: Niggli_Point  ! Niggli points
          class(CrysCell_Type),         optional, intent(out)    :: cell          ! Cell Object
          real(kind=cp), dimension(3,3),optional, intent(out)    :: trans         ! Transformation matrix 
       End Subroutine Niggli_Cell_ABC
       
       Module Subroutine Niggli_Cell_Nigglimat(N_Mat,Niggli_Point,Cell,Trans)    
          !---- Arguments ----!
          real(kind=cp),dimension(2,3),              intent(in out) :: n_mat
          real(kind=cp),dimension(5),      optional, intent(out)    :: Niggli_Point
          class(CrysCell_Type),            optional, intent(out)    :: cell
          real(kind=cp), dimension(3,3),   optional, intent(out)    :: trans  
       End Subroutine Niggli_Cell_Nigglimat
       
       Module Subroutine Niggli_Cell_Params(A,B,C,Al,Be,Ga,Niggli_Point,Cell,Trans)
          !---- Arguments ----!
          real(kind=cp),                           intent (in out)  :: a,b,c,al,be,ga
          real(kind=cp),dimension(5),    optional, intent(out)      :: Niggli_Point
          class(CrysCell_Type),          optional, intent(out)      :: cell
          real(kind=cp), dimension(3,3), optional, intent(out)      :: trans
       End Subroutine Niggli_Cell_Params
       
       Module Subroutine Niggli_Cell_Type(Cell,Niggli_Point,Celln,Trans)
          !---- Arguments ----!
          class(CrysCell_M_Type),                  intent(in out ) :: cell
          real(kind=cp),dimension(5),    optional, intent(out)     :: Niggli_Point
          class(CrysCell_Type),          optional, intent(out)     :: celln
          real(kind=cp), dimension(3,3), optional, intent(out)     :: trans   
       End Subroutine Niggli_Cell_Type
       
       Module Subroutine Niggli_Cell_Vect(A,B,C,Niggli_Point,Cell,Trans)
          !---- Arguments ----!
          real(kind=cp),dimension(3),                intent(in)     :: a,b,c
          real(kind=cp),dimension(5),      optional, intent(out)    :: Niggli_Point
          class(CrysCell_Type),            optional, intent(out)    :: cell
          real(kind=cp), dimension(3,3),   optional, intent(out)    :: trans
       End Subroutine Niggli_Cell_Vect
            
    End Interface
     

 End Module CFML_Crystal_Metrics
