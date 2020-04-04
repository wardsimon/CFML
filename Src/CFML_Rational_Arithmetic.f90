!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2017  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_Rational_Arithmetic
!!----   INFO: Extension of the module: "module_rational" from
!!----         https://rosettacode.org/wiki/Arithmetic/Rational#Fortran
!!----         It is intented for making algebra with not too big numerators
!!----         and denominators. The implementation of error control against
!!----         oveflows is not yet done. To avoid partially this the module
!!----         uses integer of 8 bytes (integer kind ik=8). Subroutines
!!----         for inverting matrices are provided: a special subroutine working
!!----         strictly with rational numbers 'Matinv_rational' and the subroutine
!!----         'rational_inv_matrix' that uses LU decomposition through double
!!----         precision and reconvesion to an approximate rational matrix.
!!----         The control of the maximum denominator is done through a global
!!----         variable called 'maximum_denominator'.
!!----         For small rational number everything seems to work, specially for
!!----         matrices representing symmetry operators in superspace and
!!----         transformation of them by not too complex change of basis.
!!----
!!---- HISTORY
!!----    Created: 01/02/2017
!!----    Last updated: 03/11/2018 (JRC)
!!----
!!---- DEPENDENCIES
!!----
!!----    CFML_GlobalDeps, CFML_String_Utilities, CFML_Math_general
!!----
!!----
!!
  Module CFML_Rational_Arithmetic

    Use CFML_GlobalDeps,       only : sp,cp,dp
    Use CFML_String_Utilities, only : Pack_String
    Use CFML_Math_general,     only : invert_matrix

    implicit none
    private

    ! Public operators and assignment
    public :: assignment (=)     ! Elemental Assignement rational<=>rational,real<=>rational, integer<=>rational of different types (sp,dp,int,int_ik)
    public :: operator (//)      ! Elemental Definition of rationals using numerator and denominator integers, e.g. 23_ik//34_ik
    public :: operator (+)       ! Elemental Addition of rationals and rational plus integers
    public :: operator (-)       ! Elemental subtraction of rationals and rational plus integers
    public :: operator (*)       ! Elemental Multiplication of rationals or rationals by integers and viceversa
    public :: operator (/)       ! Elemental division of rationals or rationals by integers and viceversa
    public :: operator (<)       ! Pure 'less than' operator comparing rational<=>rational rational<=>integer
    public :: operator (<=)      ! Pure 'less or equal than' operator comparing rational<=>rational rational<=>integer
    public :: operator (>)       ! Pure 'greater than' operator comparing rational<=>rational rational<=>integer
    public :: operator (>=)      ! Pure 'greater or equal than' operator comparing rational<=>rational rational<=>integer
    public :: operator (==)      ! Pure 'logical equal to' operator comparing rational<=>rational rational<=>integer
    public :: operator (/=)      ! Pure 'logical no-equal to' operator comparing rational<=>rational rational<=>integer

    !Public types
    public :: rational               !type definition. Constructor:  rational(numerator, denominator)

    !Public functions
    public :: rational_simplify, &   !Simplification of fractions, normally called in the majority of operations
              rational_determinant,& !Calculation Procedures using a recursive subroutine (valid only for small rationals)
              Short_Rational,  &     !Convert rational with big numerator and denominator to a limited denominator given by user
              print_rational,  &     !transform a rational type to a string like xxx/yyy
              IsInteger,       &     !Logical function applied to rationals: scalar, vector and matrix
              equal_rational_matrix,&!Logical function telling if two matrices are equal
              equal_rational_vector,&!Logical function telling if two vectors are equal
              recip,&                !Calculates the reciprocal of a rational  a/b -> b/a
              rational_colinear,&    !Logical function telling if two vectors are colinear
              rational_trace,&       !Returns the trace of a square matrix
              IsDiagonalMatrix,&     !Logical function telling if the matrix is diagonal
              IsNullVector           !Logical function telling if the vector is the null vector

    !Public subroutines
    public :: rational_inv_matrix, & !Calculates the inverse of a rational matrix using double precision arithmetic
              rational_modulo_lat, & !Reduces a translation vector to that with values in the interval [0_ik, 1_ik)
              Matinv_rational,&      !Uses rational arithmetic to invert a matrix of small rationals (no checking of overflow)
              rational_rank,&        !Computes the rank of a rational matrix
              rational_identity_matrix,& !Returns an identity matrix
              rational_rowechelonform, & !Put a matrix in a rowechelonform
              rational_smithnormalform

    !Public overloaded intrinsic functions (transpose is not needed)
    public :: abs, int, nint, modulo, mod, dot_product, maxval, minval, &
              maxloc,minloc, matmul, sum, real

    integer, public, parameter :: ik=8

    type :: rational
      integer(kind=ik) :: numerator
      integer(kind=ik) :: denominator
    end type rational

    logical,            public :: Err_Rational=.false.
    character(len=132), public :: Err_Rational_Mess
    integer(kind=ik),   public :: maximum_denominator=999_ik

    interface assignment (=)
      module procedure assign_rational_int
      module procedure assign_rational_intik
      module procedure assign_int_rational
      module procedure assign_intik_rational
      module procedure assign_rational_real_sp
      module procedure assign_rational_real_dp
      module procedure assign_real_rational_sp
      module procedure assign_real_rational_dp
    end interface

    !Constructor of a rational from two integers
    interface operator (//)
      module procedure make_rational       !input integers of 8bytes   4_ik//5_ik
      module procedure make_rational_int   !input integers of 4bytes   3//4
    end interface

    interface operator (+)
      module procedure rational_add
      module procedure rational_integer_add
      module procedure integer_rational_add
    end interface

    interface operator (-)
      module procedure rational_minus
      module procedure rational_subtract
      module procedure integer_rational_subtract
      module procedure rational_integer_subtract
    end interface

    interface operator (*)
      module procedure rational_multiply
      module procedure integer_rational_multiply
      module procedure rational_integer_multiply
    end interface

    interface operator (/)
      module procedure rational_divide
      module procedure rational_divide_int
    end interface

    interface operator (<)
      module procedure rational_lt
      module procedure rational_lt_integer
      module procedure integer_lt_rational
    end interface

    interface operator (<=)
      module procedure rational_le
      module procedure rational_le_integer
      module procedure integer_le_rational
    end interface

    interface operator (>)
      module procedure rational_gt
      module procedure rational_gt_integer
      module procedure integer_gt_rational
    end interface

    interface operator (>=)
      module procedure rational_ge
      module procedure rational_ge_integer
      module procedure integer_ge_rational
    end interface

    interface operator (==)
      module procedure rational_eq
      module procedure rational_eq_integer
      module procedure integer_eq_rational
    end interface

    interface operator (/=)
      module procedure rational_ne
      module procedure rational_ne_integer
      module procedure integer_ne_rational
    end interface

    interface abs
      module procedure rational_abs
    end interface

    interface int
      module procedure rational_int
    end interface

    interface nint
      module procedure nint_rational      !Produces an integer kind=ik
    end interface

    interface modulo
      module procedure rational_modulo
      module procedure rational_modulo_int
    end interface

    interface mod
      module procedure rational_mod
      module procedure rational_mod_int
    end interface

    interface dot_product
      module procedure rational_dot_product
    end interface

    interface maxval
      module procedure rational_maxval_vector
      module procedure rational_maxval_matrix
    end interface

    interface minval
      module procedure rational_minval_vector
      module procedure rational_minval_matrix
    end interface

    interface matmul
      module procedure rational_matmul_matvec
      !module procedure rational_matmul_vecmat
      module procedure rational_matmul_matmat
    end interface

    interface maxloc
      module procedure rational_maxloc_vect
      module procedure rational_maxloc_mat
    end interface

    interface minloc
      module procedure rational_minloc_vect
      module procedure rational_minloc_mat
    end interface

    interface real
      module procedure rational_real
    end interface

    interface IsInteger
      module procedure IsInteger_rational_scalar
      module procedure IsInteger_rational_vector
      module procedure IsInteger_rational_matrix
    end interface

    interface sum
      module procedure rational_sum_vec
    end interface

    interface Rational_RowEchelonForm
        module procedure RowEchelonForm_Rational
        module procedure RowEchelonFormT_Rational
    end interface

  contains

    elemental function make_rational (numerator, denominator) result (res)
      integer(kind=ik), intent (in) :: numerator
      integer(kind=ik), intent (in) :: denominator
      type (rational) :: res
      res = rational (numerator, denominator)
    end function make_rational

    elemental function make_rational_int (numerator, denominator) result (res)
      integer, intent (in) :: numerator
      integer, intent (in) :: denominator
      type (rational) :: res
      res = rational (int(numerator,kind=ik), int(denominator,kind=ik))
    end function make_rational_int

    pure recursive function gcd (i, j) result (res)
      integer(kind=ik), intent (in) :: i
      integer(kind=ik), intent (in) :: j
      integer(kind=ik) :: res
      if (j == 0) then
        res = i
      else
        res = gcd (j, modulo (i, j))
      end if
    end function gcd

    elemental function rational_simplify (r) result (res)
      type (rational), intent (in) :: r
      type (rational)              :: res
      integer (kind=ik) :: g
      g = gcd (r % numerator, r % denominator)
      if(g /= 0) then
        res = (r % numerator / g) // (r % denominator / g)
      else
        res= r
      end if
    end function rational_simplify

    elemental subroutine assign_rational_intik (res, i)
      type (rational),  intent (out) :: res  !, volatile
      integer(kind=ik), intent (in)  :: i
      res = i // 1_ik
    end subroutine assign_rational_intik

    elemental subroutine assign_rational_int (res, i)
      type (rational),  intent (out) :: res  !, volatile
      integer, intent (in)           :: i
      res = i // 1
    end subroutine assign_rational_int

    elemental subroutine assign_rational_real_sp (res,xr)
      type (rational),  intent(out) :: res  !, volatile
      real(kind=sp),    intent(in)  :: xr
      integer(kind=ik)              :: maxden,ai,t,si
      real(kind=sp)                 :: x,rai,eps !,startx,er1,er2
      integer(kind=ik), dimension(0:1,0:1) :: m
      maxden=maximum_denominator; eps=1.0e-6_sp
      m = 0; m(0,0)=1_ik; m(1,1)=1_ik
      si=sign(1.0,xr)
      x=abs(xr)
      !startx=x
      !First (more precise) option
      do
        ai=int(x)
        if( m(1,0)*ai+m(1,1) > maxden) exit
        t = m(0,0) * ai + m(0,1)
        m(0,1) = m(0,0); m(0,0) = t
        t = m(1,0) * ai + m(1,1)
        m(1,1) = m(1,0)
        m(1,0) =  t
        rai=real(ai,kind=sp)
        if( abs(x - rai) < eps) exit !division by zero
        x = 1.0_sp/(x - rai)
      end do
      res= si*m(0,0)// m(1,0)

      !er1=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)
      !Second option
      ! ai = (maxden - m(2,2)) / m(1,2)
      ! m(1,1) = m(1,1) * ai + mm(2,1)
      ! m(1,2) = mm(1,2) * ai + mm(2,2)
      ! er2=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)

    end subroutine assign_rational_real_sp

    elemental subroutine assign_rational_real_dp(res,xr)
      type (rational),  intent(out) :: res
      real(kind=dp),    intent (in) :: xr
      integer(kind=ik)              :: maxden,ai,t,si
      real(kind=dp)                 :: x,rai,eps
      integer(kind=ik), dimension(0:1,0:1)  :: m

      maxden=maximum_denominator; eps=1.0e-8_dp
      m = 0; m(0,0)=1_ik; m(1,1)=1_ik
      si=sign(1.0_dp,xr)
      x=abs(xr)

      do
        ai=int(x)
        if( m(1,0)*ai+m(1,1) > maxden) exit
        t = m(0,0) * ai + m(0,1)
        m(0,1) = m(0,0); m(0,0) = t
        t = m(1,0) * ai + m(1,1)
        m(1,1) = m(1,0)
        m(1,0) =  t
        rai=real(ai,kind=dp)
        if( abs(x - rai ) < eps) exit !division by zero
        x = 1.0_dp/(x - rai)
      end do
      res= si*m(0,0)// m(1,0)
    end subroutine assign_rational_real_dp

    elemental function recip(r) result (reciprocal)
      type(rational), intent (in) :: r
      type(rational)              :: reciprocal
      reciprocal= 0_ik
      if(r%numerator /= 0_ik) &
        reciprocal= r%denominator // r%numerator
    end function recip

    elemental subroutine assign_int_rational (i, res)
      type (rational), intent (in)   :: res  !, volatile
      integer,         intent (out)  :: i
      i= nint(real(res%numerator,kind=dp)/real(res%denominator,kind=dp))
    end subroutine assign_int_rational

    elemental subroutine assign_intik_rational (i, res)
      type (rational), intent (in)   :: res  !, volatile
      integer(kind=ik),intent (out)  :: i
      i= nint(real(res%numerator,kind=dp)/real(res%denominator,kind=dp))
    end subroutine assign_intik_rational

    elemental subroutine assign_real_rational_sp (x,res)
      type (rational), intent(in)   :: res
      real(kind=sp),   intent (out) :: x
      x=real(res%numerator,kind=sp)/real(res%denominator,kind=sp)
    end subroutine assign_real_rational_sp

    elemental subroutine assign_real_rational_dp (x,res)
      type (rational), intent(in)   :: res
      real(kind=dp),   intent (out) :: x
      x=real(res%numerator,kind=dp)/real(res%denominator,kind=dp)
    end subroutine assign_real_rational_dp

    elemental function rational_add (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % denominator + r % denominator * s % numerator // &
          & r % denominator * s % denominator
      res=rational_simplify(res)
    end function rational_add

    elemental function integer_rational_add (i, s) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: s
      type (rational) :: res
      res = (i * s % denominator +  s % numerator) // s % denominator
      res=rational_simplify(res)
    end function integer_rational_add

    elemental function rational_integer_add (s, i) result (res)
      type (rational), intent (in) :: s
      integer(kind=ik),intent (in) :: i
      type (rational) :: res
      res = (i * s % denominator +  s % numerator) // s % denominator
      res=rational_simplify(res)
    end function rational_integer_add

    elemental function rational_minus (r) result (res)
      type (rational), intent (in) :: r
      type (rational) :: res
      res = - r % numerator // r % denominator
      res=rational_simplify(res)
    end function rational_minus

    elemental function rational_subtract (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % denominator - r % denominator * s % numerator // &
          & r % denominator * s % denominator
      res=rational_simplify(res)
    end function rational_subtract

    elemental function integer_rational_subtract (i, s) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: s
      type (rational) :: res
      res = (i * s % denominator -  s % numerator) //  s % denominator
      res=rational_simplify(res)
    end function integer_rational_subtract

    elemental function rational_integer_subtract (s,i) result (res)
      type (rational), intent (in) :: s
      integer(kind=ik),intent (in) :: i
      type (rational) :: res
      res = (s % numerator - i * s % denominator) //  s % denominator
      res=rational_simplify(res)
    end function rational_integer_subtract

    elemental function rational_multiply (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % numerator // (r % denominator * s % denominator)
      res=rational_simplify(res)
    end function rational_multiply

    elemental function integer_rational_multiply (i, s) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: s
      type (rational) :: res
      res = i * s % numerator // s % denominator
      res=rational_simplify(res)
    end function integer_rational_multiply

    elemental function rational_integer_multiply (s,i) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: s
      type (rational) :: res
      res = i * s % numerator // s % denominator
      res=rational_simplify(res)
    end function rational_integer_multiply

    elemental function rational_divide (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational)  :: res
      integer(kind=ik) :: denom
      denom = r % denominator * s % numerator
      if(denom /= 0) then
        res = r % numerator * s % denominator // denom
        res=rational_simplify(res)
      else
        res=0_ik//1_ik
      end if
    end function rational_divide

    elemental function rational_divide_int (r, is) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),intent (in) :: is
      type (rational) :: res
      if(is /= 0) then
        res = r % numerator // (r % denominator*is)
        res=rational_simplify(res)
      else
        res=0_ik//1_ik
      end if
    end function rational_divide_int

    elemental function int_divide_rational (is,r) result (res)
      integer(kind=ik),intent (in) :: is
      type (rational), intent (in) :: r
      type (rational) :: res
      if(r /= 0_ik/1_ik) then
        res = r % denominator // (r % numerator*is)
        res=rational_simplify(res)
      else
        res=0_ik//1_ik
      end if
    end function int_divide_rational

    pure function rational_lt (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: r_simple
      type (rational) :: s_simple
      logical :: res
      r_simple = rational_simplify (r)
      s_simple = rational_simplify (s)
      res = r_simple % numerator * s_simple % denominator < &
          & s_simple % numerator * r_simple % denominator
    end function rational_lt

    pure function rational_lt_integer (r, i) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),intent (in) :: i
      logical                      :: res
      type (rational) :: r_simple

      r_simple = rational_simplify (r)
      res = r_simple % numerator < i * r_simple % denominator
    end function rational_lt_integer

    pure function integer_lt_rational (i,r) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: r
      logical                      :: res
      type (rational) :: r_simple
      r_simple = rational_simplify (r)
      res = i * r_simple % denominator < r_simple % numerator
    end function integer_lt_rational

    pure function rational_le (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: r_simple
      type (rational) :: s_simple
      logical :: res
      r_simple = rational_simplify (r)
      s_simple = rational_simplify (s)
      res = r_simple % numerator * s_simple % denominator <= &
          & s_simple % numerator * r_simple % denominator
    end function rational_le

    pure function rational_le_integer (r, i) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),intent (in) :: i
      logical                      :: res
      type (rational) :: r_simple
      r_simple = rational_simplify (r)
      res = r_simple % numerator <= i * r_simple % denominator
    end function rational_le_integer

    pure function integer_le_rational (i, r) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: r
      logical                      :: res
      type (rational) :: r_simple
      r_simple = rational_simplify (r)
      res = i * r_simple % denominator  <= r_simple % numerator
    end function integer_le_rational


    pure function rational_gt (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      logical                      :: res
      type (rational) :: r_simple
      type (rational) :: s_simple
      r_simple = rational_simplify (r)
      s_simple = rational_simplify (s)
      res = r_simple % numerator * s_simple % denominator > &
          & s_simple % numerator * r_simple % denominator
    end function rational_gt

    pure function rational_gt_integer (r, i) result (res)
      type (rational),  intent (in) :: r
      integer(kind=ik), intent (in) :: i
      logical                      :: res
      type (rational) :: r_simple
      r_simple = rational_simplify (r)
      res = r_simple % numerator > i * r_simple % denominator
    end function rational_gt_integer

    pure function integer_gt_rational (i, r) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: r
      logical                      :: res
      type (rational) :: r_simple
      r_simple = rational_simplify (r)
      res = i * r_simple % denominator > r_simple % numerator
    end function integer_gt_rational

    pure function rational_ge (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: r_simple
      type (rational) :: s_simple
      logical :: res
      r_simple = rational_simplify (r)
      s_simple = rational_simplify (s)
      res = r_simple % numerator * s_simple % denominator >= &
          & s_simple % numerator * r_simple % denominator
    end function rational_ge

    pure function rational_ge_integer (r, i) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),         intent (in) :: i
      logical                      :: res
      type (rational) :: r_simple
      r_simple = rational_simplify (r)
      res = r_simple % numerator >= i * r_simple % denominator
    end function rational_ge_integer

    pure function integer_ge_rational (i, r) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: r
      logical                      :: res
      type (rational) :: r_simple
      r_simple = rational_simplify (r)
      res = i * r_simple % denominator >= r_simple % numerator
    end function integer_ge_rational

    elemental function rational_eq (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      logical :: res
      res = r % numerator * s % denominator == s % numerator * r % denominator
    end function rational_eq

    Pure function rational_eq_integer (r, i) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),         intent (in) :: i
      logical :: res
      res = r % denominator == 1 .and. r % numerator == i
    end function rational_eq_integer

    pure function integer_eq_rational(i, r) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: r
      logical :: res
      res = r % denominator == 1 .and. r % numerator == i
    end function integer_eq_rational

    pure function rational_ne (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      logical :: res
      res = r % numerator * s % denominator /= s % numerator * r % denominator
    end function rational_ne

    pure function rational_ne_integer (r, i) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),         intent (in) :: i
      logical :: res
      res = r % numerator /= i * r % denominator
    end function rational_ne_integer

    pure function integer_ne_rational(i, r) result (res)
      integer(kind=ik),intent (in) :: i
      type (rational), intent (in) :: r
      logical :: res
      res = r % numerator /= i * r % denominator
    end function integer_ne_rational

    elemental function rational_abs (r) result (res)
      type (rational), intent (in) :: r
      type (rational) :: res
      res = sign (r % numerator, r % denominator) // r % denominator
    end function rational_abs

    elemental function rational_int (r) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik) :: res
      res = r % numerator / r % denominator
    end function rational_int

    elemental function rational_real(r) result (res)
      type (rational),  intent (in) :: r
      real(kind=dp) :: res
      res = real(r % numerator,kind=dp) / real(r % denominator,kind=dp)
    end function rational_real

    elemental function nint_rational (r) result(res)
      type (rational),  intent (in)  :: r
      integer(kind=ik)               :: res
      res = nint(real(r%numerator,kind=dp)/real(r%denominator,kind=dp),kind=ik)
    end function nint_rational

    elemental subroutine rational_modulo_lat (r)
      type (rational), intent (in out) :: r

      do
      	 if(r < 0_ik//1_ik) then
      	 	 r=r+1_ik
      	 else
      	 	 exit
      	 end if
      end do

      do
      	 if(r >= 1_ik//1_ik) then
      	 	 r=r-1_ik
      	 else
      	 	 exit
      	 end if
      end do

    end subroutine rational_modulo_lat

    elemental function rational_modulo (r) result (res)
      type (rational), intent (in) :: r
      integer :: res
      res = modulo (r % numerator, r % denominator)
    end function rational_modulo

    elemental function rational_modulo_int (r,i) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),intent(in)  :: i
      type (rational) :: res
      real(kind=dp) :: val
      val = modulo (real(r % numerator,kind=dp)/ real(r % denominator,kind=dp),real(i,kind=dp))
      res = val
    end function rational_modulo_int

    elemental function rational_mod (r) result (res)
      type (rational), intent (in) :: r
      integer :: res
      res = mod (r % numerator, r % denominator)
    end function rational_mod

    !Fortran MOD and MODULO
    !Syntax: MOD (a, p)  MODULO(a, p)
    !Arguments
    !a must be of type INTEGER or REAL.
    !p must be of the same type and kind as a. Its value must not be zero.
    !The result is the same type and kind as a
    !Results
    ! MOD -> Its value is a - INT(a / p) * p.
    ! MODULO ->
    ! If a is a REAL, the result value is: a - FLOOR(a / p) * p.
    ! If a is an integer, MODULO(a, p) has the value r such that
    ! a = q * p + r, where q is an INTEGER and r is nearer to zero than p.

    elemental function rational_mod_int (r,i) result (res)
      type (rational),  intent (in) :: r
      integer(kind=ik), intent (in) :: i
      type (rational)               :: res
      real(kind=dp) :: val
      val = mod (real(r % numerator,kind=dp)/ real(r % denominator,kind=dp),real(i,kind=dp))
      res = val
    end function rational_mod_int

    Pure function rational_dot_product (r1,r2) result (res)
      type (rational), dimension(:), intent (in) :: r1,r2
      type (rational) :: res
      integer :: n1,n2,i,nm
      n1=size(r1); n2=size(r2)
      res=0_ik//1_ik
      nm=min(n1,n2)
      do i=1,nm
        res = rational_simplify (res + r1(i)*r2(i))
      end do
    end function rational_dot_Product

    Pure function rational_matmul_matvec (mat,vec) result (vec_out)
      type (rational), dimension(:,:), intent (in) :: mat
      type (rational), dimension(:),   intent (in) :: vec
      type (rational), dimension(size(vec))        :: vec_out

      integer :: n1,n2,n3,i,nm
      n1=size(mat,dim=1); n2=size(mat,dim=2); n3=size(vec)
      nm=min(n2,n3)
      vec_out=0_ik//1_ik
      do i=1,nm
        vec_out(i) = rational_simplify(dot_product(mat(i,1:nm),vec(1:nm)))
      end do
    end function rational_matmul_matvec

    !function rational_matmul_vecmat (vec,mat) result (vec_out)
    !  type (rational), dimension(:),   intent (in) :: vec
    !  type (rational), dimension(:,:), intent (in) :: mat
    !  type (rational), dimension(size(vec))        :: vec_out
    !
    !  integer :: n1,n2,n3,i
    !  n1=size(mat,dim=1); n2=size(mat,dim=2); n3=size(vec)
    !  if(n1 == n2 .and. n2 == n3) then
    !    do i=1,n3
    !      vec_out(i) = rational_simplify (dot_product(vec,mat(:,i)))
    !    end do
    !  else
    !    Err_Rational=.true.
    !    write(unit=Err_Rational_Mess,fmt="(a,3i4)") "Error in MATMUL: the dimensions of the arguments are non-conformable",n1,n2,n3
    !  end if
    !end function rational_matmul_vecmat

    Function Short_Rational(r,maxdenom) result(sh)
      type(rational),   intent (in) :: r
      integer(kind=ik), intent (in) :: maxdenom
      type(rational)  :: sh
      real(kind=dp)    :: val
      integer(kind=ik) :: ival
      ival=maximum_denominator
      maximum_denominator=maxdenom
      val=real(r%numerator,kind=dp)/real(r%denominator,kind=dp)
      sh=val
      maximum_denominator=ival
    End Function Short_Rational

    Pure function rational_matmul_matmat(mat1,mat2) result (mat_out)
      type(rational), dimension(:,:), intent (in) :: mat1
      type(rational), dimension(:,:), intent (in) :: mat2
      type(rational), dimension(size(mat1,dim=1),size(mat2,dim=2)) :: mat_out

      integer :: n1,n2,n3,n4,i,j,nm
      n1=size(mat1,dim=1); n2=size(mat1,dim=2); n3=size(mat2,dim=1); n4=size(mat2,dim=2)
      forall ( i = 1:n1 , j = 1:n4 ) mat_out(i,j) = 0_ik/1_ik
      if(n2 == n3) then
        do j=1,n4
          do i=1,n1
            mat_out(i,j) = rational_simplify (dot_product(mat1(i,:),mat2(:,j)))
          end do
        end do
      else
        nm=min(n2,n3)
        do j=1,n4
          do i=1,n1
            mat_out(i,j) = rational_simplify (dot_product(mat1(i,1:nm),mat2(1:nm,j)))
          end do
        end do
      end if
    end function rational_matmul_matmat

    Pure function rational_maxval_vector (r) result (res)
      type (rational), dimension(:), intent (in) :: r
      type (rational) :: res
      integer :: i,n
      n=size(r)
      res=-huge(1_ik)//1_ik
      do i=1,n
        if(r(i) > res) res=r(i)
      end do
    end function rational_maxval_vector

    Pure function rational_maxval_matrix (r) result (res)
      type (rational), dimension(:,:), intent (in) :: r
      type (rational) :: res
      integer :: i,j,n1,n2
      n1=size(r,dim=1);  n2=size(r,dim=2)
      res=-huge(1_ik)//1_ik
      do j=1,n2
        do i=1,n1
          if(r(i,j) > res) res=r(i,j)
        end do
      end do
    end function rational_maxval_matrix

    Pure function rational_minval_vector (r) result (res)
      type (rational), dimension(:), intent (in) :: r
      type (rational) :: res
      integer :: i,n
      n=size(r)
      res=huge(1_ik)//1_ik
      do i=1,n
        if(r(i) < res) res=r(i)
      end do
    end function rational_minval_vector

    Pure function rational_minval_matrix (r) result (res)
      type (rational), dimension(:,:), intent (in) :: r
      type (rational) :: res
      integer :: i,j,n1,n2
      n1=size(r,dim=1);  n2=size(r,dim=2)
      res=huge(1_ik)//1_ik
      do j=1,n2
        do i=1,n1
          if(r(i,j) < res) res=r(i,j)
        end do
      end do
    end function rational_minval_matrix

    elemental function print_rational(r) result (rtxt)
      type(rational),  intent(in)  :: r
      character(len=81)            :: rtxt
      type(rational) :: sr
      if(r%denominator /= 0_ik) then
        sr=rational_simplify(r)
      else
        sr=r
      end if
      if(sr%denominator == 1_ik) then
        write(unit=rtxt,fmt="(i40)") sr%numerator
      else
        write(unit=rtxt,fmt="(i40,a,i40)") sr%numerator,"/",sr%denominator
      end if
      rtxt=adjustl(Pack_String(rtxt))
    end function print_rational

    !!---- Pure recursive Function rational_determinant(a) result(det)
    !!----   type(rational), dimension(:,:), intent(in) :: a
    !!----   type(rational) :: det
    !!----
    !!---- This function for calculating the determinant is not very efficient but
    !!---- largely enough for our needs. Will be probably replaced by another one
    !!---- using LU decomposition
    !!----
    Pure recursive Function rational_determinant(a) result(det)
      type(rational), dimension(:,:), intent(in) :: a
      type(rational) :: det
      !Local variables
      type(rational), dimension(size(a,dim=1)-1, size(a,dim=1)-1) :: b
      type(rational) :: sgn
      integer :: i, n
      n=size(a,dim=1)
      if (n == 1) then
        det = a(1,1)
      else
        det = 0_ik//1_ik
        sgn = 1_ik/1_ik
        do i=1,n
          b(:, :(i-1)) = a(2:, :i-1)
          b(:, i:) = a(2:, i+1:)
          det = det + sgn * a(1, i) * rational_determinant(b)
          sgn = sgn * (-1_ik/1_ik)
        end do
      end if
    End Function rational_determinant

    !Pure function rational_determinant(Mat) result(det)
    !  type(rational), dimension(:,:), intent(in) :: Mat
    !  type(rational) :: det
    !  !Local variables
    !  real(kind=cp), dimension(size(Mat,dim=1),size(Mat,dim=2)) :: A
    !  real(kind=cp) :: determ
    !  integer:: n1,n2,nm
    !
    !  n1=size(Mat,dim=1); n2=size(Mat,dim=2)
    !  nm=min(n1,n2)
    !  A=Mat
    !  call Determinant(A(1:nm,1:nm),nm,determ)
    !  det=determ
    !end function rational_determinant

    !!---- Subroutine rational_inv_matrix(Mat,invMat)
    !!----   type(rational), dimension(:,:), intent(in)  :: Mat
    !!----   type(rational), dimension(:,:), intent(out) :: invMat
    !!----
    !!----  This calculates the inverse of a matrix converting it previously to
    !!----  a double precision matrix. The final invMat is an approximation according
    !!----  to the value of the 'maximum_denominator' value
    !!----
    Subroutine rational_inv_matrix(Mat,invMat)
      type(rational), dimension(:,:), intent(in)  :: Mat
      type(rational), dimension(:,:), intent(out) :: invMat
      !Local variables
      real(kind=dp), dimension(size(Mat,dim=1),size(Mat,dim=2)) :: A,invA
      integer :: n1,n2
      logical :: singular

      singular=.false.
      err_rational=.false.
      n1=size(Mat,dim=1); n2=size(Mat,dim=2)
      if(n1 == n2) then
         A=Mat
         call Invert_Matrix(A,invA,singular)
         if(singular) then
            Err_Rational=.true.
            write(unit=Err_Rational_Mess,fmt="(a)") "Singular Matrix!"
         else
            invMat=invA  !This gives an approximation according to the 'maximum_denominator' value
         end if
      else
        Err_Rational=.true.
        write(unit=Err_Rational_Mess,fmt="(a)") "Error: the provided matrix is not square!"
      end if
    end Subroutine rational_inv_matrix

    !!---- Subroutine Matinv_rational(a,ainv)
    !!----    !---- Arguments ----!
    !!----    type(rational), dimension(:,:), intent(in) :: a
    !!----    type(rational), dimension(:,:), intent(out):: ainv
    !!----
    !!----    This subroutine uses rational types for calculating the inverse of
    !!----    a rational matrix. It works only with small rational and dimension
    !!----    of the matrices (the maximum dimension may be n=9-10)
    !!----
    Subroutine Matinv_rational(a,ainv)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent(in) :: a
       type(rational), dimension(:,:), intent(out):: ainv
       !---- Local variables ----!
       type(rational)                :: amax,savec
       integer, dimension(size(a,1)) :: ink,jnk
       integer                       :: i,j,k,l,n

       !---- Subroutine to invert a rational matrix ----!
       ainv=a
       n=size(a,1)
       do k=1,n
          amax=0_ik//1_ik
          do
             do
                do i=k,n
                   do j=k,n
                      if (abs(amax)-abs(ainv(i,j)) > 0_ik//1_ik) cycle
                      amax=ainv(i,j)
                      ink(k)=i
                      jnk(k)=j
                   end do
                end do
                i=ink(k)
                if (i-k < 0) cycle
                exit
             end do

             if (i-k /= 0) then
                do j=1,n
                   savec=ainv(k,j)
                   ainv(k,j)=ainv(i,j)
                   ainv(i,j)=-savec
                end do
             end if

             j=jnk(k)
             if (j-k < 0) cycle
             exit
          end do

          if (j-k /= 0) then
             do i=1,n
                savec=ainv(i,k)
                ainv(i,k)=ainv(i,j)
                ainv(i,j)=-savec
             end do
          end if

          do i=1,n
             if (i-k /= 0)  then
                ainv(i,k)=-ainv(i,k)*recip(amax)
             end if
          end do
          do i=1,n
             do j=1,n
                if (i-k == 0 .or. j-k == 0) cycle
                ainv(i,j)=ainv(i,j)+ainv(i,k)*ainv(k,j)
             end do
          end do
          do j=1,n
             if (j-k == 0)   cycle
             ainv(k,j)=ainv(k,j)*recip(amax)
          end do
          ainv(k,k)=recip(amax)
       end do     !k

       do l=1,n
          k=n-l+1
          j=ink(k)
          if (j-k > 0) then
             do i=1,n
                savec=ainv(i,k)
                ainv(i,k)=-ainv(i,j)
                ainv(i,j)=savec
             end do
          end if
          i=jnk(k)
          if (i-k > 0) then
             do j=1,n
                savec=ainv(k,j)
                ainv(k,j)=-ainv(i,j)
                ainv(i,j)=savec
             end do
          end if
       end do

       return
    End Subroutine Matinv_rational

    pure function rational_maxloc_mat(Mat) result(pos_max)
      type(rational),  dimension(:,:), intent(in) :: Mat
      integer(kind=ik),dimension(2)               :: pos_max
      !Local variables
      integer:: nu1,nl1,nu2,nl2,i,j
      type(rational) :: res
      nu1=ubound(Mat,dim=1); nu2=ubound(Mat,dim=2)
      nl1=lbound(Mat,dim=1); nl2=lbound(Mat,dim=2)

      res=-huge(1_ik)//1_ik
      do j=nl2,nu2
        do i=nl1,nu1
          if(mat(i,j) > res) then
            res=mat(i,j)
            pos_max=[i,j]
          end if
        end do
      end do
    end function rational_maxloc_mat

    pure function rational_maxloc_vect(vec) result(pos_max)
      type(rational), dimension(:), intent(in) :: vec
      integer                                  :: pos_max
      !Local variables
      integer:: nu,nl,i
      type(rational) :: res
      nu=ubound(vec,dim=1)
      nl=lbound(vec,dim=1)
      res=-huge(1_ik)//1_ik
      do i=nl,nu
        if(vec(i) > res) then
          res=vec(i)
          pos_max=i
        end if
      end do
    end function rational_maxloc_vect

    pure function rational_minloc_mat(Mat) result(pos_min)
      type(rational),  dimension(:,:), intent(in) :: Mat
      integer(kind=ik),dimension(2)               :: pos_min
      !Local variables
      integer:: nu1,nl1,nu2,nl2,i,j
      type(rational) :: res
      nu1=ubound(Mat,dim=1); nu2=ubound(Mat,dim=2)
      nl1=lbound(Mat,dim=1); nl2=lbound(Mat,dim=2)

      res=huge(1_ik)//1_ik
      do j=nl2,nu2
        do i=nl1,nu1
          if(mat(i,j) < res) then
            res=mat(i,j)
            pos_min=[i,j]
          end if
        end do
      end do

    end function rational_minloc_mat

    pure function rational_minloc_vect(vec) result(pos_min)
      type(rational), dimension(:), intent(in) :: vec
      integer                                  :: pos_min
      !Local variables
      integer:: nu,nl
      integer:: i
      type(rational) :: res
      nu=ubound(vec,dim=1)
      nl=lbound(vec,dim=1)
      res=huge(1_ik)//1_ik
      do i=nl,nu
        if(vec(i) < res) then
          res=vec(i)
          pos_min=i
        end if
      end do
    end function rational_minloc_vect

    pure function rational_sum_vec(vec) result(suma)
      type(rational), dimension(:), intent(in) :: vec
      type(rational)                           :: suma
      !Local variables
      integer :: i,n
      n=size(vec)
      suma=0_ik/1_ik
      do i=1,n
        suma=suma+vec(i)
      end do
    end function rational_sum_vec

    pure function equal_rational_vector(vec1,vec2) result(eq)
      type(rational), dimension(:), intent(in) :: vec1,vec2
      logical                                  :: eq      !
      integer:: i,n1,n2

      n1=size(vec1); n2=size(vec2)
      eq=.false.
      if(n1 /= n2) return
      do i=1,n1
        if(vec1(i) /= vec2(i)) return
      end do
      eq=.true.
    end function equal_rational_vector

    pure function equal_rational_matrix(Mat1,Mat2) result(eq)
      type(rational), dimension(:,:), intent(in) :: Mat1,Mat2
      logical                                    :: eq
      !
      integer:: i,j,n11,n12,n21,n22

      n11=size(Mat1,dim=1); n12=size(Mat1,dim=2)
      n21=size(Mat2,dim=1); n22=size(Mat2,dim=2)
      eq=.false.
      if(n11 /= n21 .or. n12 /= n22) return
      do j=1,n12
        do i=1,n11
          if(Mat1(i,j) /= Mat2(i,j)) return
        end do
      end do
      eq=.true.
    end function equal_rational_matrix

    pure function IsInteger_rational_matrix(Mat) result(is_int)
      type(rational), dimension(:,:), intent(in) :: Mat
      logical                                    :: is_int
      !
      integer:: i,j,n1,n2

      n1=size(Mat,dim=1); n2=size(Mat,dim=2)
      is_int=.false.
      do j=1,n1
        do i=1,n2
          if(Mat(i,j)%denominator /= 1_ik ) return
          !if (mod(Mat(i,j)%Numerator,Mat(i,j)%Denominator) /= 0) return
        end do
      end do
      is_int=.true.
    end function IsInteger_rational_matrix

    pure function IsInteger_rational_vector(vec) result(is_int)
      type(rational), dimension(:), intent(in) :: vec
      logical                                  :: is_int
      !
      integer:: i,n

      n=size(vec)
      is_int=.false.
      do i=1,n
        if(vec(i)%denominator /= 1_ik ) return
      end do
      is_int=.true.
    end function IsInteger_rational_vector

    pure function IsInteger_rational_scalar(rat) result(is_int)
      type(rational), intent(in) :: rat
      logical         :: is_int
      !
      is_int=.false.
      if(rat%denominator == 1_ik ) is_int=.true.
    end function IsInteger_rational_scalar

    function Rational_Colinear(a,b,n) Result(co_linear)

        !---- Argument ----!
        type(rational), dimension(n), intent(in) :: a
        type(rational), dimension(n), intent(in) :: b
        integer,                      intent(in) :: n
        logical                                  :: co_linear

        !---- Local variables ----!
        integer        :: i,ia,ib
        type(rational) :: c

        co_linear=.true.
        do i=1,n
            if (abs(a(i)%numerator) > 0) then
                ia=i
                exit
            end if
        end do
        do i=1,n
            if (abs(b(i)%numerator) > 0) then
                ib=i
                exit
            end if
        end do
        if (ia /= ib) then
            co_linear=.false.
            return
        else
            c = a(ia) / b(ib)
            do i = 1 , n
                if ((a(i) - c * b(i) /= (0//1))) then
                    co_linear=.false.
                    return
                end if
            end do
        end if

    end function Rational_Colinear

    function Rational_Trace(Mat) result(tr)

        !---- Arguments ----!
        type(rational), dimension(:,:), intent(in) :: Mat
        type(rational) :: tr

        !---- Local variables ----!
        integer:: n1,n2,i

        n1=size(Mat,dim=1); n2=size(Mat,dim=2)
        if(n1 == n2) then
            tr = 0
            do i = 1 , n1
                tr = tr + Mat(i,i)
            end do
        else
            Err_Rational=.true.
            write(unit=Err_Rational_Mess,fmt="(a)") "Error in Trace: matrix is not a square matrix!"
        end if

    end function Rational_Trace

    logical function IsDiagonalMatrix(A) Result(diagonal)

        !---- Arguments ----!
        type(rational), dimension(:,:), intent(in) :: A

        !---- Local variables ----!
        integer :: i,j

        diagonal = .true.
        !if (size(A,1) /= size(A,2)) then
        !    diagonal = .false.
        !    return
        !end if
        do i = 1 , size(A,1)
            do j = 1 , size(A,2)
                if (i /= j .and. A(i,j) /= (0//1)) then
                    diagonal = .false.
                    return
                end if
            end do
        end do

    end function IsDiagonalMatrix

    logical function IsNullVector(v) Result(null)

        !---- Arguments ----!
        type(rational), dimension(:), intent(in) :: v

        !---- Local variables ----!
        integer :: i

        null = .true.
        do i = 1 , size(v)
            if (v(i) /= (0//1)) then
                null = .false.
                return
            end if
        end do

    end function IsNullVector

    subroutine Rational_Rank(M,r)

        !---- Arguments ----!
        type(rational), dimension(:,:), intent(in)  :: M
        integer,                        intent(out) :: r

        !---- Local variables ----!
        integer :: i,nNull
        type(rational), dimension(:),   allocatable :: nullVector
        type(rational), dimension(:,:), allocatable :: U

        allocate(U(size(M,1),size(M,2)),nullVector(size(M,2)))

        U = M
        nNull = 0
        nullVector(:) = 0 // 1

        call Rational_RowEchelonForm(U)
        do i = 1 , size(U,1)
            if (equal_rational_vector(U(i,:),nullVector)) nNull = nNull + 1
        end do

        r = size(M,1) - nNull

    end subroutine Rational_Rank

    subroutine Rational_Identity_Matrix(n,I)

        !---- Arguments ----!
        integer, intent(in)            :: n
        type(rational), dimension(n,n) :: I

        !---- Local variables ----!
        integer :: j

        I(:,:) = 0 // 1
        do j = 1 , n
            I(j,j) = 1 // 1
        end do

    end subroutine Rational_Identity_Matrix

    !!---- Subroutine RowEchelonForm_Rational(M)
    !!----    integer, dimension(:,:), intent(inout) :: M
    !!----
    !!---- Fortran version of RowEchelonForm from the CrystGAP package.
    !!---- Adapted for rational matrices. The original source code
    !!---- can be found at:
    !!----      https://fossies.org/linux/gap/pkg/cryst/gap/common.gi
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine RowEchelonForm_Rational(M)

        !---- Arguments ----!
        type(rational), dimension(:,:), intent(inout) :: M

        !---- Local variables ----!
        integer :: r,c,i,j
        integer :: nrow,ncolumn
        logical :: cleared
        type(rational) :: a
        type(rational), dimension(:), allocatable :: row

        nrow    = size(M,1)
        ncolumn = size(M,2)
        allocate(row(ncolumn))
        r = 1  ! index for rows
        c = 1  ! index for columns

        do
            if (r > nrow .or. c > ncolumn) exit
            i = r
            do
                !if ( i > r .or. M(i,c) /= 0 ) exit
                if (i > nrow) exit
                if (M(i,c) /= (0//1)) exit
                i = i + 1
            end do

            if ( i <= nrow ) then
                row(:) = M(r,:)
                M(r,:) = M(i,:)
                M(i,:) = row(:)
                do j = i + 1 , nrow
                    a = abs(M(j,c))
                    if ( a /= (0//1) .and. a < abs(M(r,c)) ) then
                        row(:) = M(r,:)
                        M(r,:) = M(j,:)
                        M(j,:) = row(:)
                    end if
                end do
                if ( M(r,c) < (0//1) ) M(r,:) = -M(r,:)
                cleared = .true.
                do i = r + 1 , nrow
                    a = M(i,c)/M(r,c)
                    if ( a /= (0//1) ) M(i,:) = M(i,:) - a * M(r,:)
                    if ( M(i,c) /= (0//1) ) cleared = .false.
                end do
                if ( cleared ) then
                    r = r + 1
                    c = c + 1
                end if
            else
                c = c + 1
            end if
        end do

    end subroutine RowEchelonForm_Rational

    !!---- Subroutine RowEchelonFormT_Rational
    !!----      integer, dimension(:,:), intent(inout) :: M
    !!----      integer, dimension(:,:), intent(inout) :: T
    !!----
    !!---- Fortran version of RowEchelonFormT from the CrystGAP package
    !!---- The original source code can be found at:
    !!----         https://fossies.org/linux/gap/pkg/cryst/gap/common.gi
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine RowEchelonFormT_Rational(M,T)

        !---- Arguments ----!
        type(rational), dimension(:,:), intent(inout) :: M
        type(rational), dimension(:,:), intent(inout) :: T

        !---- Local variables ----!
        integer        :: r,c,i,j
        integer        :: nrow,ncolumn
        logical        :: cleared
        type(rational) :: a
        type(rational), dimension(:), allocatable :: row, Trow

        nrow    = size(M,1)
        ncolumn = size(M,2)
        allocate(row(ncolumn))
        allocate(Trow(nrow))

        r = 1  ! index for rows
        c = 1  ! index for columns

        do
            if (r > nrow .or. c > ncolumn) exit
            i = r

            do
                !if ( i > r .or. M(i,c) /= 0 ) exit
                if (i > nrow) exit
                if (M(i,c) /= (0//1)) exit
                i = i + 1
            end do

            if ( i <= nrow ) then
                row(:)  = M(r,:)
                M(r,:)  = M(i,:)
                M(i,:)  = row(:)
                Trow(:) = T(r,:)
                T(r,:)  = T(i,:)
                T(i,:)  = Trow(:)
                do j = i + 1 , nrow
                    a = abs(M(j,c))
                    if ( a /= (0//1) .and. a < abs(M(r,c)) ) then
                        row(:)  = M(r,:)
                        M(r,:)  = M(j,:)
                        M(j,:)  = row(:)
                        Trow(:) = T(r,:)
                        T(r,:)  = T(j,:)
                        T(j,:)  = Trow(:)
                    end if
                end do
                if ( M(r,c) < (0//1) ) then
                    M(r,:) = - M(r,:)
                    T(r,:) = - T(r,:)
                end if
                cleared = .true.
                do i = r + 1 , nrow
                    a = M(i,c)/M(r,c)
                    if ( a /= (0//1) ) then
                        M(i,:) = M(i,:) - a * M(r,:)
                        T(i,:) = T(i,:) - a * T(r,:)
                    end if
                    if ( M(i,c) /= (0//1) ) cleared = .false.
                end do
                if ( cleared ) then
                    r = r + 1
                    c = c + 1
                end if
            else
                c = c + 1
            end if
        end do

    end subroutine RowEchelonFormT_Rational

    subroutine Rational_SmithNormalForm(M,nr,nc,D,P,Q)

        !---- Arguments ----!
        type(rational), dimension(nr,nc), intent(in)  :: M
        integer,                          intent(in)  :: nr
        integer,                          intent(in)  :: nc
        type(rational), dimension(nr,nc), intent(out) :: D
        type(rational), dimension(nr,nr), intent(out) :: P
        type(rational), dimension(nc,nc), intent(out) :: Q

        !---- Local variables ----!
        integer                          :: ndiag
        type(rational), dimension(nc,nr) :: Dt

        ! P and Q must be initialized to the identity matrix
        call Rational_Identity_Matrix(nr,P)
        call Rational_Identity_Matrix(nc,Q)

        D = M
        ndiag = 0

        do
            if (mod(ndiag,2) == 0) then
                call Rational_RowEchelonForm(D,P)
                ndiag = ndiag + 1
                Dt = transpose(D)
            else
                call Rational_RowEchelonForm(Dt,Q)
                ndiag = ndiag + 1
                D = transpose(Dt)
            end if
            if (IsDiagonalMatrix(D)) exit
            if (ndiag > 100) then
                err_rational = .true.
                err_rational_mess = "Error in Rational_SmithNormalForm. Unable to diagonalize matrix."
                return
            end if
        end do

        Q = transpose(Q)

    end subroutine Rational_SmithNormalForm

  End Module CFML_Rational_Arithmetic