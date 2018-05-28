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
!!---- MODULE: CFML_Rational_Arithmetic
!!----   INFO: Extension of the module: "module_rational" from
!!----         https://rosettacode.org/wiki/Arithmetic/Rational#Fortran
!!----         It is intented for making algebra with not too big numerators
!!----         and denominators. The implementation of error control against
!!----         oveflows is not yet done.
!!----
!!
 Module CFML_Rational_Arithmetic
    !---- Use Modules ----!
    Use CFML_GlobalDeps,       only : cp, dp, il, err_cfml
    Use CFML_String_Utilities, only : Pack_String
    Use CFML_Math_general,     only : determinant,invert_matrix

    !---- Variables ----!
    implicit none
    
    private

    !> Public operators and assignment
    public :: assignment (=)
    public :: operator (//)
    public :: operator (+)
    public :: operator (-)
    public :: operator (*)
    public :: operator (/)
    public :: operator (<)
    public :: operator (<=)
    public :: operator (>)
    public :: operator (>=)
    public :: operator (==)
    public :: operator (/=)

    !> Public functions
    public :: rational_simplify, rational_determinant, &       !Calculation Procedures
              string_rational, &                               !transform a rational type to a string like xxx/yyy
              equal_rational_matrix,equal_rational_vector

    !> Public subroutines
    public :: rational_inv_matrix, rational_modulo_lat


    !> Public overloaded intrinsic functions (transpose is not needed)
    public :: abs, int, nint, modulo, mod, dot_product, maxval, minval, &
              maxloc,minloc, matmul, sum

    !> Parameters
    integer(kind=IL), parameter :: MAXIMUM_DENOMINATOR=999_IL
    
    !> Types
    Type, Public :: Rational
       integer(kind=IL) :: numerator
       integer(kind=IL) :: denominator
    End Type Rational

    Interface assignment (=)
       module procedure assign_rational_int
       module procedure assign_rational_int_il
       module procedure assign_int_rational
       module procedure assign_intil_rational
       module procedure assign_rational_real_cp
       module procedure assign_rational_real_dp
       module procedure assign_real_rational_cp
       module procedure assign_real_rational_dp
    end Interface

    !> Constructor of a rational from two integers
    Interface operator (//)
       module procedure make_rational
       module procedure make_rational_int
    end Interface

    Interface operator (+)
       module procedure rational_add
       module procedure rational_integer_add
       module procedure integer_rational_add
    end Interface

    Interface operator (-)
       module procedure rational_minus
       module procedure rational_subtract
       module procedure integer_rational_subtract
       module procedure rational_integer_subtract
    end Interface

    Interface operator (*)
       module procedure rational_multiply
       module procedure integer_rational_multiply
       module procedure rational_integer_multiply
    end Interface

    Interface operator (/)
       module procedure rational_divide
       module procedure rational_divide_int
    end Interface

    Interface operator (<)
       module procedure rational_lt
       module procedure rational_lt_integer
       module procedure integer_lt_rational
    end Interface

    Interface operator (<=)
       module procedure rational_le
       module procedure rational_le_integer
       module procedure integer_le_rational
    end Interface

    Interface operator (>)
       module procedure rational_gt
       module procedure rational_gt_integer
       module procedure integer_gt_rational
    end Interface

    Interface operator (>=)
       module procedure rational_ge
       module procedure rational_ge_integer
       module procedure integer_ge_rational
    end Interface

    Interface operator (==)
       module procedure rational_eq
       module procedure rational_eq_integer
       module procedure integer_eq_rational
    end Interface

    Interface operator (/=)
       module procedure rational_ne
       module procedure rational_ne_integer
       module procedure integer_ne_rational
    end Interface

    Interface abs
       module procedure rational_abs
    end Interface

    Interface int
       module procedure rational_int
    end Interface

    Interface nint
       module procedure nint_rational
    end Interface

    Interface modulo
       module procedure rational_modulo
       module procedure rational_modulo_int
    end Interface

    Interface mod
       module procedure rational_mod
       module procedure rational_mod_int
    end Interface

    Interface dot_product
       module procedure rational_dot_product
    end Interface

    Interface maxval
       module procedure rational_maxval_vector
       module procedure rational_maxval_matrix
    end Interface

    Interface minval
       module procedure rational_minval_vector
       module procedure rational_minval_matrix
    end Interface

    Interface matmul
       module procedure rational_matmul_matvec
       module procedure rational_matmul_matmat
    end Interface

    Interface maxloc
       module procedure rational_maxloc_vect
       module procedure rational_maxloc_mat
    end Interface

    Interface minloc
       module procedure rational_minloc_vect
       module procedure rational_minloc_mat
    end Interface

    Interface sum
       module procedure rational_sum_vec
    end Interface
    
    Interface
       Module Pure Subroutine Assign_Rational_Int(Res, I)
          !---- Arguments ----!
          type(rational),  intent (out) :: res  ! volatile
          integer,         intent (in)  :: i
       End Subroutine Assign_Rational_Int
       
       Module Pure Subroutine Assign_Rational_Int_IL(Res, I)
          !---- Arguments ----!
          type(rational),   intent (out) :: res  ! volatile
          integer(kind=il), intent (in)  :: i
       End Subroutine Assign_Rational_Int_IL
       
       Module Elemental Subroutine Assign_Int_Rational(I, Res)
          !---- Arguments ----!
          type(rational), intent (in)   :: res  !, volatile
          integer,        intent (out)  :: i
       End Subroutine Assign_Int_Rational
       
       Module Elemental Subroutine Assign_Intil_Rational(I, Res)
          !---- Arguments ----!
          type(rational),  intent (in)   :: res  !, volatile
          integer(kind=il),intent (out)  :: i
       End Subroutine Assign_Intil_Rational
       
       Module Elemental Subroutine Assign_Rational_Real_Cp(Res, Xr)
          !---- Arguments ----!
          type(rational), intent(out) :: res  ! volatile
          real(kind=cp),  intent (in) :: xr
       End Subroutine Assign_Rational_Real_Cp
       
       Module Elemental Subroutine Assign_Rational_Real_Dp(Res,Xr)
          !---- Arguments ----!
          type(rational), intent(out) :: res
          real(kind=dp),  intent (in) :: xr
       End Subroutine Assign_Rational_Real_Dp
       
       Module Elemental Subroutine Assign_Real_Rational_Cp(X, Res)
          !---- Arguments ----!
          type(rational), intent(in)   :: res
          real(kind=cp),  intent (out) :: x
       End Subroutine Assign_Real_Rational_Cp
       
       Module Elemental Subroutine Assign_Real_Rational_Dp(X, Res)
          !---- Arguments ----!
          type(rational), intent(in)   :: res
          real(kind=dp),  intent (out) :: x
       End Subroutine Assign_Real_Rational_Dp
       
       Module Elemental Function Make_Rational(Numerator, Denominator) Result(Res)
          !---- Arguments ----!
          integer(kind=il), intent (in) :: numerator
          integer(kind=il), intent (in) :: denominator
          type(rational)                :: res
       End Function Make_Rational
       
       Module Elemental Function Make_Rational_Int(Numerator, Denominator) Result(Res)
          !---- Arguments ----! 
          integer, intent (in) :: numerator
          integer, intent (in) :: denominator
          type(rational)      :: res
       End Function Make_Rational_Int
       
       Module Elemental Function Rational_Add(R, S) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          type(rational)              :: res
       End Function Rational_Add
       
       Module Elemental Function Rational_Integer_Add(S, I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: s
          integer(kind=il),intent (in) :: i
          type(rational)               :: res
       End Function Rational_Integer_Add
       
       Module Elemental Function Integer_Rational_Add(I, S) Result(Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: s
          type(rational)               :: res
       End Function Integer_Rational_Add
       
       Module Elemental Function Rational_Minus(R) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational)              :: res
       End Function Rational_Minus
       
       Module Elemental Function Rational_Subtract(R, S) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          type(rational)              :: res
       End Function Rational_Subtract
       
       Module Elemental Function Integer_Rational_Subtract(I, S) Result(Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: s
          type(rational)               :: res
       End Function Integer_Rational_Subtract

       Module Elemental Function Rational_Integer_Subtract(S,I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: s
          integer(kind=il),intent (in) :: i
          type(rational)               :: res
       End Function Rational_Integer_Subtract
       
       Module Elemental Function Rational_Multiply(R, S) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          type(rational)              :: res
       End Function Rational_Multiply
       
       Module Elemental Function Integer_Rational_Multiply(I, S) Result(Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: s
          type(rational)               :: res
       End Function Integer_Rational_Multiply
       
       Module Elemental Function Rational_Integer_Multiply(S,I) Result(Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: s
          type(rational)               :: res
       End Function Rational_Integer_Multiply
       
       Module Elemental Function Rational_Divide(R, S) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          type(rational)              :: res
       End Function Rational_Divide
       
       Module Elemental Function Rational_Divide_Int(R, I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: r
          integer(kind=il),intent (in) :: i
          type(rational)               :: res
       End Function Rational_Divide_Int
       
       Module Pure Function Rational_Lt(R, S) Result(Res)
          !---- Arguments ----! 
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          logical                     :: res
       End Function Rational_Lt
       
       Module Pure Function Rational_Lt_Integer(R, I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: r
          integer(kind=il),intent (in) :: i
          logical                      :: res
       End Function Rational_Lt_Integer
       
       Module Pure Function Integer_Lt_Rational(I,R) Result (Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: r
          logical                      :: res
       End Function Integer_Lt_Rational
       
       Module Pure Function Rational_LE(R, S) Result(Res)
          !---- Argument ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          logical                     :: res
       End Function Rational_LE
       
       Module Pure Function Rational_Le_Integer(R, I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: r
          integer(kind=il),intent (in) :: i
          logical                      :: res
       End Function Rational_Le_Integer
       
       Module Pure Function Integer_Le_Rational(I, R) Result(Res)
          !---- Arguments ----! 
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: r
          logical                      :: res
       End Function Integer_Le_Rational
       
       Module Pure Function Rational_GT(R, S) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          logical                     :: res
       End Function Rational_GT
       
       Module Pure Function Rational_Gt_Integer(R, I) Result(Res)
          !---- Arguments ----!
          type(rational),   intent (in) :: r
          integer(kind=il), intent (in) :: i
          logical                      :: res
       End Function Rational_Gt_Integer
       
       Module Pure Function Integer_Gt_Rational(I, R) Result(Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: r
          logical                      :: res
       End Function Integer_Gt_Rational
       
       Module Pure Function Rational_GE(R, S) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          logical                     :: res
       End Function Rational_GE
       
       Module Pure Function Rational_GE_Integer(R, I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: r
          integer(kind=il),intent (in) :: i
          logical                      :: res
       End Function Rational_GE_Integer
       
       Module Pure Function Integer_GE_Rational(I, R) Result(Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: r
          logical                      :: res
       End Function Integer_GE_Rational
       
       Module Elemental Function Rational_EQ(R, S) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          logical                     :: res
       End Function Rational_EQ
       
       Module Elemental Function Rational_EQ_Integer(R, I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: r
          integer(kind=il),intent (in) :: i
          logical                      :: res
       End Function Rational_EQ_Integer
       
       Module Elemental Function Integer_Eq_Rational(I, R) Result(Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: r
          logical                      :: res
       End Function Integer_EQ_Rational
       
       Module Pure Function Rational_Ne(R, S) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational), intent (in) :: s
          logical                     :: res
       End Function Rational_Ne
       
       Module Elemental Function Rational_Ne_Integer(R, I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: r
          integer(kind=il),intent (in) :: i
          logical                      :: res
       End Function Rational_Ne_Integer
       
       Module Elemental Function Integer_Ne_Rational(I, R) Result(Res)
          !---- Arguments ----!
          integer(kind=il),intent (in) :: i
          type(rational),  intent (in) :: r
          logical                      :: res
       End Function Integer_Ne_Rational
       
       Module Elemental Function Rational_Abs(R) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          type(rational)              :: res
       End Function Rational_Abs
       
       Module Elemental Function Rational_Int(R) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          integer(kind=il)            :: res
       End Function Rational_Int
       
       Module Elemental Function Nint_Rational(R) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in)  :: r
          integer(kind=il)              :: res
       End Function Nint_Rational
       
       Module Function Rational_Dot_Product(R1,R2) Result(Res)
          !---- Arguments ----!
          type(rational), dimension(:), intent (in) :: r1
          type(rational), dimension(:), intent (in) :: r2
          type(rational)                            :: res
       End Function Rational_Dot_Product
       
       Module Pure Function Rational_Sum_Vec(Vec) Result(Suma)
          !---- Arguments ----!
          type(rational), dimension(:), intent(in) :: vec
          type(rational)                           :: suma
       End Function Rational_Sum_Vec
       
       Module Elemental Function Rational_Modulo(R) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          integer                     :: res
       End Function Rational_Modulo
       
       Module Elemental Function Rational_Modulo_Int(R,I) Result(Res)
          !---- Arguments ----!
          type(rational),  intent (in) :: r
          integer(kind=il),intent(in)  :: i
          type(rational)               :: res
       End Function Rational_Modulo_Int
       
       Module Elemental Function Rational_Mod(R) Result(Res)
          !---- Arguments ----!
          type(rational), intent (in) :: r
          integer                     :: res
       End Function Rational_Mod
       
       Module Elemental Function Rational_Mod_Int(R,I) Result(Res)
          !---- Arguments ----!
         type(rational),   intent (in) :: r
         integer(kind=il), intent (in) :: i
         type(rational)                :: res
       End Function Rational_Mod_Int
       
       Module Function Rational_Minval_Vector(R) Result(Res)
          !---- Arguments ----!
          type(rational), dimension(:), intent (in) :: r
          type(rational)                            :: res
       End Function Rational_Minval_Vector
       
       Module Function Rational_Minval_Matrix(R) Result(Res)
          !---- Arguments ----!
          type(rational), dimension(:,:), intent (in) :: r
          type(rational)                              :: res
       End Function Rational_Minval_Matrix
       
       Module Function Rational_Maxval_Vector(R) Result(Res)
          !---- Arguments ----!
          type(rational), dimension(:), intent (in) :: r
          type(rational)                            :: res
       End Function Rational_Maxval_Vector
       
       Module Function Rational_Maxval_Matrix(R) Result(Res)
          !---- Arguments ----!
          type(rational), dimension(:,:), intent (in) :: r
          type(rational)                              :: res
       End Function Rational_Maxval_Matrix
       
       Module Function Rational_Matmul_Matvec(Mat,Vec) Result(Vec_Out)
          !---- Arguments ----!
          type(rational), dimension(:,:), intent (in) :: mat
          type(rational), dimension(:),   intent (in) :: vec
          type(rational), dimension(size(vec))        :: vec_out
       End Function Rational_Matmul_Matvec
       
       Module Function Rational_Matmul_Matmat(Mat1,Mat2) Result(Mat_Out)
          !---- Arguments ----! 
          type(rational), dimension(:,:), intent (in)                 :: mat1
          type(rational), dimension(:,:), intent (in)                 :: mat2
          type(rational),dimension(size(mat1,dim=1),size(mat2,dim=2)) :: mat_out
       End Function Rational_Matmul_Matmat
       
       Module Pure Function Rational_Maxloc_Mat(Mat) Result(Pos_Max)
          !---- Arguments ----!
          type(rational),  dimension(:,:), intent(in) :: Mat
          integer(kind=il),dimension(2)               :: pos_max
       End Function Rational_Maxloc_Mat
       
       Module Pure Function Rational_Maxloc_Vect(Vec) Result(Pos_Max)
          !---- Arguments ----!
          type(rational), dimension(:), intent(in) :: vec
          integer                                  :: pos_max
       End Function Rational_Maxloc_Vect
       
       Module Pure Function Rational_Minloc_Mat(Mat) Result(Pos_Min)
          !---- Arguments ----!
          type(rational),  dimension(:,:), intent(in) :: Mat
          integer(kind=il),dimension(2)               :: pos_min
       End Function Rational_Minloc_Mat
       
       Module Pure Function Rational_Minloc_Vect(Vec) Result(Pos_Min)
          !---- Arguments ----!
          type(rational), dimension(:), intent(in) :: vec
          integer                                  :: pos_min
       End Function Rational_Minloc_Vect
    
    End Interface

 Contains

    !!----
    !!---- FUNCTION RATIONAL_NINT
    !!----
    !!----
    !!
    Elemental Function Rational_Nint(R) Result(Res)
       !---- Arguments ----!
       real(kind=cp), intent (in)  :: r
       type(rational)              :: res
       
       res = nint(r,kind=il) // 1_il
       
       return
    End Function Rational_Nint

    !!----
    !!---- FUNCTION GCD
    !!----
    !!----
    !!
    Pure Recursive Function Gcd(I, J) Result(Res)
       !---- Arguments ----!
       integer(kind=il), intent (in) :: i
       integer(kind=il), intent (in) :: j
       integer(kind=il)              :: res
      
       if (j == 0_il) then
          res = i
       else
          res = gcd (j, modulo (i, j))
       end if
       
       return
    End Function Gcd

    !!----
    !!---- FUNCTION RATIONAL_SIMPLIFY
    !!----
    !!----
    !! 
    Elemental Function Rational_Simplify(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational)              :: res
      
       !---- Local Variables ----!
       integer(kind=il) :: g
      
       g = gcd(r%numerator, r%denominator)
       if (g /= 0_il) then
          res = (r%numerator / g) // (r%denominator / g)
       else
          res= r
       end if
       
       return
    End Function Rational_Simplify

    !!----
    !!---- SUBROUTINE RATIONAL_MODULO_LAT
    !!----
    !!----
    !!
    Elemental Subroutine Rational_Modulo_Lat(R)
       !---- Arguments ----!
       type(rational), intent (in out) :: r

       do
          if (r < 0_il//1_il) then
             r=r+1_il
          else
             exit
          end if
       end do

       do
          if (r >= 1_il//1_il) then
             r=r-1_il
          else
             exit
          end if
       end do

       return
    End Subroutine Rational_Modulo_Lat

    !!----
    !!---- FUNCTION STRING_RATIONAL
    !!----
    !!----
    !!
    Elemental Function String_Rational(R) Result (Str)
       !---- Arguments ----!
       type(rational),  intent(in)  :: r
       character(len=50)            :: Str
      
       !---- Local Variable ----!
       type(rational) :: sr
      
       if (r%denominator /= 0_il) then
          sr = rational_simplify(r)
       else
          sr = r
       end if
      
       if (sr%denominator == 1_il) then
          write(unit=str,fmt="(i20)") sr%numerator
       else
          write(unit=str,fmt="(i20,a,i20)") sr%numerator,"/",sr%denominator
       end if
      
       str=adjustl(Pack_String(str))
       
       return
    End Function String_Rational

    !!----
    !!---- FUNCTION RATIONAL_DETERMINANT
    !!----
    !!----
    !!
    Function Rational_Determinant(Mat) Result(Det)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent(in) :: Mat
       type(rational)                             :: det
      
       !---- Local variables ----!
       real(kind=cp), dimension(size(Mat,dim=1),size(Mat,dim=2)) :: A
       real(kind=cp) :: determ
       integer:: n1,n2

       n1=size(Mat,dim=1); n2=size(Mat,dim=2)
       if (n1 == n2) then
          A=Mat
          call Determinant(A,n1,determ)
          det=determ
       else
          Err_CFML%state=.true.
          err_cfml%flag=2
          write(unit=err_cfml%msg,fmt="(a)") &
               "Error in Determinant: the provided matrix is not square!"
       end if
       
       return
    End Function Rational_Determinant
    
    !!----
    !!---- SUBROUTINE RATIONAL_INV_MATRIX
    !!----
    !!----
    !!
    Subroutine Rational_Inv_Matrix(Mat,Invmat)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent(in)  :: Mat
       type(rational), dimension(:,:), intent(out) :: invMat
      
      
       !---- Local variables ----!
       real(kind=cp), dimension(size(Mat,dim=1),size(Mat,dim=2)) :: A,invA
       integer :: n1,n2
       logical :: singular

       singular=.false.
       n1=size(Mat,dim=1); n2=size(Mat,dim=2)
       if (n1 == n2) then
          A=Mat
          call Invert_Matrix(A,invA,singular)
          if (singular) then
             Err_CFML%state=.true.
             err_cfml%flag=2
             write(unit=err_cfml%msg,fmt="(a)") "Singular Matrix!"
          else
             invMat=invA
          end if
       else
          Err_CFML%state=.true.
          err_cfml%flag=2
          write(unit=err_cfml%msg,fmt="(a)") &
               "Error in Determinant: the provided matrix is not square!"
       end if
       
       return
    End Subroutine Rational_Inv_Matrix

    !!----
    !!---- FUNCTION EQUAL_RATIONAL_VECTOR
    !!----
    !!----
    !!
    Pure Function Equal_Rational_Vector(Vec1,Vec2) Result(Eq)
       !---- Arguments ----!
       type(rational), dimension(:), intent(in) :: vec1
       type(rational), dimension(:), intent(in) :: vec2
       logical                                  :: eq      
       
       !---- Local Variables ----!
       integer:: i,n1,n2

       n1=size(vec1); n2=size(vec2)
       eq=.false.
       if (n1 /= n2) return
      
       do i=1,n1
          if (vec1(i) /= vec2(i)) return
       end do
       eq=.true.
       
       return
    End Function Equal_Rational_Vector

    !!----
    !!---- FUNCTION EQUAL_RATIONAL_MATRIX
    !!----
    !!----
    !!
    Pure Function Equal_Rational_Matrix(Mat1,Mat2) Result(Eq)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent(in) :: Mat1
       type(rational), dimension(:,:), intent(in) :: Mat2
       logical                                    :: eq
      
       !---- Local Variables ----!
       integer:: i,j,n11,n12,n21,n22

       n11=size(Mat1,dim=1); n12=size(Mat1,dim=2)
       n21=size(Mat2,dim=1); n22=size(Mat2,dim=2)
       eq=.false.
       if (n11 /= n21 .or. n12 /= n22) return
      
       do j=1,n12
          do i=1,n11
             if (Mat1(i,j) /= Mat2(i,j)) return
          end do
       end do
       eq=.true.
       
       return
    End Function Equal_Rational_Matrix

  End Module CFML_Rational_Arithmetic