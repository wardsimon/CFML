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
!!----         oveflows is not yet done.
!!----
!!---- HISTORY
!!----    Created: 01/02/2017
!!----
!!---- DEPENDENCIES
!!----
!!----    CFML_GlobalDeps, CFML_String_Utilities, CFML_Math_general
!!----
!!----
!!
  Module CFML_Rational_Arithmetic_test

    Use CFML_GlobalDeps,       only : cp,dp
    Use CFML_String_Utilities, only : Pack_String
    Use CFML_Math_general,     only : determinant,invert_matrix

    implicit none
    private

    ! Public operators and assignment
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

    !Public types
    public :: rational

    !Public functions
    public :: rational_simplify, rational_determinant,& !Calculation Procedures
              print_rational, &  !transform a rational type to a string like xxx/yyy
              equal_rational_matrix,equal_rational_vector

    !Public subroutines
    public :: rational_inv_matrix, rational_modulo_lat


    !Public overloaded intrinsic functions (transpose is not needed)
    public :: abs, int, nint, modulo, mod, dot_product, maxval, minval, &
              maxloc,minloc, matmul, sum

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
      module procedure assign_rational_real_cp
      module procedure assign_rational_real_dp
      module procedure assign_real_rational_cp
      module procedure assign_real_rational_dp
    end interface

    !Constructor of a rational from two integers
    interface operator (//)
      module procedure make_rational
      module procedure make_rational_int
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
      module procedure nint_rational
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

    interface sum
      module procedure rational_sum_vec
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

    pure subroutine assign_rational_intik (res, i)
      type (rational),  intent (out) :: res  !, volatile
      integer(kind=ik), intent (in)  :: i
      res = i // 1_ik
    end subroutine assign_rational_intik
    
    pure subroutine assign_rational_int (res, i)
      type (rational),  intent (out) :: res  !, volatile
      integer, intent (in)           :: i
      res = i // 1 
    end subroutine assign_rational_int

    elemental subroutine assign_rational_real_cp (res,xr)
      type (rational), intent(out) :: res  !, volatile
      real(kind=cp),   intent (in) :: xr
      integer(kind=ik)             :: maxden,ai,t,si
      real(kind=cp)                :: x,rai,eps !,startx,er1,er2
      integer(kind=ik), dimension(0:1,0:1) :: m

      maxden=maximum_denominator; eps=1.0e-6_cp
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
        rai=real(ai,kind=cp)
        if( abs(x - rai) < eps) exit !division by zero
        x = 1.0_cp/(x - rai)
      end do
      res= si*m(0,0)// m(1,0)

      !er1=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)
      !Second option
      ! ai = (maxden - m(2,2)) / m(1,2)
      ! m(1,1) = m(1,1) * ai + mm(2,1)
      ! m(1,2) = mm(1,2) * ai + mm(2,2)
      ! er2=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)

    end subroutine assign_rational_real_cp

    elemental subroutine assign_rational_real_dp (res,xr)
      type (rational), intent(out) :: res
      real(kind=dp),   intent (in) :: xr
      integer(kind=ik)             :: maxden,ai,t,si
      real(kind=dp)                :: x,rai,eps
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

    elemental subroutine assign_real_rational_cp (x,res)
      type (rational), intent(in)   :: res
      real(kind=cp),   intent (out) :: x
      x=real(res%numerator,kind=cp)/real(res%denominator,kind=cp)
    end subroutine assign_real_rational_cp

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
        res = r % numerator // (r % numerator*is)
        res=rational_simplify(res)
      else
        res=0_ik//1_ik
      end if
    end function rational_divide_int

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

    elemental function rational_eq_integer (r, i) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),         intent (in) :: i
      logical :: res
      res = r % denominator == 1 .and. r % numerator == i
    end function rational_eq_integer

    elemental function integer_eq_rational(i, r) result (res)
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

    elemental function rational_ne_integer (r, i) result (res)
      type (rational), intent (in) :: r
      integer(kind=ik),         intent (in) :: i
      logical :: res
      res = r % numerator /= i * r % denominator
    end function rational_ne_integer

    elemental function integer_ne_rational(i, r) result (res)
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

    elemental function rational_nint (r) result(res)
      real(kind=cp),  intent (in)  :: r
      type (rational)              :: res
      res = nint(r,kind=ik) // 1_ik
    end function rational_nint

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

    function rational_dot_product (r1,r2) result (res)
      type (rational), dimension(:), intent (in) :: r1,r2
      type (rational) :: res
      integer :: n1,n2,i
      n1=size(r1); n2=size(r2)
      res=0_ik//1_ik
      if(n1 == n2) then
        do i=1,n1
          res = rational_simplify (res + r1(i)*r2(i))
        end do
      else
        Err_Rational=.true.
        write(unit=Err_Rational_Mess,fmt="(a,2i4)") "Error in DOT_PRODUCT: the dimensions of the arguments are different",n1,n2
      end if
    end function rational_dot_Product

    function rational_matmul_matvec (mat,vec) result (vec_out)
      type (rational), dimension(:,:), intent (in) :: mat
      type (rational), dimension(:),   intent (in) :: vec
      type (rational), dimension(size(vec))        :: vec_out

      integer :: n1,n2,n3,i
      n1=size(mat,dim=1); n2=size(mat,dim=2); n3=size(vec)
      if(n1 == n2 .and. n2 == n3) then
        do i=1,n3
          vec_out(i) = rational_simplify (dot_product(mat(i,:),vec))
        end do
      else
        Err_Rational=.true.
        write(unit=Err_Rational_Mess,fmt="(a,3i4)") "Error in MATMUL: the dimensions of the arguments are non-conformable",n1,n2,n3
      end if
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

    function rational_matmul_matmat (mat1,mat2) result (mat_out)
      type (rational), dimension(:,:), intent (in) :: mat1
      type (rational), dimension(:,:), intent (in) :: mat2
      type (rational),dimension(size(mat1,dim=1),size(mat2,dim=2)) :: mat_out

      integer :: n1,n2,n3,n4,i,j
      n1=size(mat1,dim=1); n2=size(mat1,dim=2); n3=size(mat2,dim=1); n4=size(mat2,dim=2)
      if(n2 == n3) then
        do j=1,n4
          do i=1,n1
            mat_out(i,j) = rational_simplify (dot_product(mat1(i,:),mat2(:,j)))
          end do
        end do
      else
        Err_Rational=.true.
        write(unit=Err_Rational_Mess,fmt="(a,4i4)") "Error in MATMUL: the dimensions of the arguments are non-conformable",n1,n2,n3,n4
      end if
    end function rational_matmul_matmat

    function rational_maxval_vector (r) result (res)
      type (rational), dimension(:), intent (in) :: r
      type (rational) :: res
      integer :: i,n
      n=size(r)
      res=-huge(1_ik)//1_ik
      do i=1,n
        if(r(i) > res) res=r(i)
      end do
    end function rational_maxval_vector

    function rational_maxval_matrix (r) result (res)
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

    function rational_minval_vector (r) result (res)
      type (rational), dimension(:), intent (in) :: r
      type (rational) :: res
      integer :: i,n
      n=size(r)
      res=huge(1_ik)//1_ik
      do i=1,n
        if(r(i) < res) res=r(i)
      end do
    end function rational_minval_vector

    function rational_minval_matrix (r) result (res)
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
      character(len=50)            :: rtxt
      type(rational) :: sr
      if(r%denominator /= 0_ik) then
        sr=rational_simplify(r)
      else
        sr=r
      end if
      if(sr%denominator == 1_ik) then
        write(unit=rtxt,fmt="(i20)") sr%numerator
      else
        write(unit=rtxt,fmt="(i20,a,i20)") sr%numerator,"/",sr%denominator
      end if
      rtxt=adjustl(Pack_String(rtxt))
    end function print_rational

    function rational_determinant(Mat) result(det)
      type(rational), dimension(:,:), intent(in) :: Mat
      type(rational) :: det
      !Local variables
      real(kind=cp), dimension(size(Mat,dim=1),size(Mat,dim=2)) :: A
      real(kind=cp) :: determ
      integer:: n1,n2

      n1=size(Mat,dim=1); n2=size(Mat,dim=2)
      if(n1 == n2) then
         A=Mat
         call Determinant(A,n1,determ)
         det=determ
      else
        Err_Rational=.true.
        write(unit=Err_Rational_Mess,fmt="(a)") "Error in Determinant: the provided matrix is not square!"
      end if
    end function rational_determinant

    Subroutine rational_inv_matrix(Mat,invMat)
      type(rational), dimension(:,:), intent(in)  :: Mat
      type(rational), dimension(:,:), intent(out) :: invMat
      !Local variables
      real(kind=cp), dimension(size(Mat,dim=1),size(Mat,dim=2)) :: A,invA
      integer :: n1,n2
      logical :: singular

      singular=.false.
      n1=size(Mat,dim=1); n2=size(Mat,dim=2)
      if(n1 == n2) then
         A=Mat
         call Invert_Matrix(A,invA,singular)
         if(singular) then
            Err_Rational=.true.
            write(unit=Err_Rational_Mess,fmt="(a)") "Singular Matrix!"
         else
            invMat=invA
         end if
      else
        Err_Rational=.true.
        write(unit=Err_Rational_Mess,fmt="(a)") "Error in Determinant: the provided matrix is not square!"
      end if
    end Subroutine rational_inv_matrix

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
      real(kind=dp), dimension(size(vec)) :: rvec
      rvec=vec
      suma=sum(rvec)
      suma=rational_simplify(suma)
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

  End Module CFML_Rational_Arithmetic_test