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
    Use CFML_GlobalDeps,   only : CP, eps, PI, TO_RAD, err_CFML, clear_error
    Use CFML_Math_General, only : Co_Prime, swap, Sort, Co_Linear
    Use CFML_Math_3D,      only : Invert_Array3x3, determ_3x3, determ_Vec, Cross_Product

    implicit none

    private

    
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

    Interface
       Module Pure Function Metrics(cell,ang) Result(G)
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
             
    End Interface
     
 Contains

    

 End Module CFML_Crystal_Metrics
