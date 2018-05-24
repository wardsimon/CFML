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
    Use CFML_GlobalDeps,       only : cp, dp, il, err_cfml, err_cfml_flag, err_cfml_msg
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

 Contains

    !!----
    !!---- FUNCTION MAKE_RATIONAL
    !!----
    !!----
    !!
    Elemental Function Make_Rational(Numerator, Denominator) Result(Res)
       !---- Arguments ----!
       integer(kind=il), intent (in) :: numerator
       integer(kind=il), intent (in) :: denominator
       type(rational)                :: res
       
       res = rational(numerator, denominator)
       
       return
    End Function Make_Rational

    !!----
    !!---- FUNCTION MAKE_RATIONAL_INT
    !!----
    !!----
    !!
    Elemental Function Make_Rational_Int(Numerator, Denominator) Result(Res)
       !---- Arguments ----! 
       integer, intent (in) :: numerator
       integer, intent (in) :: denominator
       type(rational)      :: res
      
       res = rational(int(numerator,kind=il), int(denominator,kind=il))
       
       return
    End Function Make_Rational_Int

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
    !!---- FUNCTION ASSIGN_RATIONAL_INT_IL 
    !!----
    !!----
    !! 
    Pure Subroutine Assign_Rational_Int_IL(Res, I)
       !---- Arguments ----!
       type(rational),   intent (out) :: res  ! volatile
       integer(kind=il), intent (in)  :: i
      
       res = i // 1_il
       
       return
    End Subroutine Assign_Rational_Int_IL
    
    !!----
    !!---- FUNCTION ASSIGN_RATIONAL_INT_IT 
    !!----
    !!----
    !! 
    Pure Subroutine Assign_Rational_Int(Res, I)
       !---- Arguments ----!
       type(rational),  intent (out) :: res  ! volatile
       integer,         intent (in)  :: i
      
       res = i // 1
       
       return 
    End Subroutine Assign_Rational_Int

    !!----
    !!---- SUBROUTINE ASSIGN_RATIONAL_REAL_CP 
    !!----
    !!----
    !! 
    Elemental Subroutine Assign_Rational_Real_Cp(Res, Xr)
       !---- Arguments ----!
       type(rational), intent(out) :: res  ! volatile
       real(kind=cp),  intent (in) :: xr
      
       !---- Local variables ----!      
       integer(kind=il)                     :: maxden,ai,t,si
       integer(kind=il), dimension(0:1,0:1) :: m
       real(kind=cp)                        :: x,rai,eps !,startx,er1,er2

       maxden=maximum_denominator; eps=1.0e-6_cp
       m = 0; m(0,0)=1_il; m(1,1)=1_il
       si=sign(1.0,xr)
       x=abs(xr)
      
       !startx=x
       !First (more precise) option
       do
          ai=int(x)
          if ( m(1,0)*ai+m(1,1) > maxden) exit
          t = m(0,0) * ai + m(0,1)
          m(0,1) = m(0,0); m(0,0) = t
          t = m(1,0) * ai + m(1,1)
          m(1,1) = m(1,0)
          m(1,0) =  t
          rai=real(ai,kind=cp)
          if ( abs(x - rai) < eps) exit !division by zero
          x = 1.0_cp/(x - rai)
       end do
      
       res= si*m(0,0)// m(1,0)

       !er1=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)
       !Second option
       ! ai = (maxden - m(2,2)) / m(1,2)
       ! m(1,1) = m(1,1) * ai + mm(2,1)
       ! m(1,2) = mm(1,2) * ai + mm(2,2)
       ! er2=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)

       return
    End Subroutine Assign_Rational_Real_Cp
    
    !!----
    !!---- SUBROUTINE ASSIGN_RATIONAL_REAL_DP 
    !!----
    !!----
    !! 
    Elemental Subroutine Assign_Rational_Real_Dp(Res,Xr)
       !---- Arguments ----!
       type(rational), intent(out) :: res
       real(kind=dp),  intent (in) :: xr
      
       !---- Local variables ----!
       integer(kind=il)                      :: maxden,ai,t,si
       integer(kind=il), dimension(0:1,0:1)  :: m
       real(kind=dp)                         :: x,rai,eps

       maxden=maximum_denominator; eps=1.0e-8_dp
       m = 0; m(0,0)=1_il; m(1,1)=1_il
       si=sign(1.0_dp,xr)
       x=abs(xr)

       do
          ai=int(x)
          if ( m(1,0)*ai+m(1,1) > maxden) exit
          t = m(0,0) * ai + m(0,1)
          m(0,1) = m(0,0); m(0,0) = t
          t = m(1,0) * ai + m(1,1)
          m(1,1) = m(1,0)
          m(1,0) =  t
          rai=real(ai,kind=dp)
          if ( abs(x - rai ) < eps) exit !division by zero
          x = 1.0_dp/(x - rai)
       end do
       res= si*m(0,0)// m(1,0)
       
       return
    End Subroutine Assign_Rational_Real_Dp
    
    !!----
    !!---- SUBROUTINE ASSIGN_INT_RATIONAL 
    !!----
    !!----
    !! 
    Elemental Subroutine Assign_Int_Rational(I, Res)
       !---- Arguments ----!
       type(rational), intent (in)   :: res  !, volatile
       integer,        intent (out)  :: i
      
       i= nint(real(res%numerator,kind=dp)/real(res%denominator,kind=dp))
       
       return
    End Subroutine Assign_Int_Rational
    
    !!----
    !!---- SUBROUTINE ASSIGN_INTIL_RATIONAL 
    !!----
    !!----
    !! 
    Elemental Subroutine Assign_Intil_Rational(I, Res)
       !---- Arguments ----!
       type(rational),  intent (in)   :: res  !, volatile
       integer(kind=il),intent (out)  :: i
      
       i= nint(real(res%numerator,kind=dp)/real(res%denominator,kind=dp))
       
       return
    End Subroutine Assign_Intil_Rational

    !!----
    !!---- SUBROUTINE ASSIGN_REAL_RATIONAL_CP 
    !!----
    !!----
    !! 
    Elemental Subroutine Assign_Real_Rational_Cp(X, Res)
       !---- Arguments ----!
       type(rational), intent(in)   :: res
       real(kind=cp),  intent (out) :: x
      
       x=real(res%numerator,kind=cp)/real(res%denominator,kind=cp)
       
       return
    End Subroutine Assign_Real_Rational_Cp
    
    !!----
    !!---- SUBROUTINE ASSIGN_REAL_RATIONAL_DP 
    !!----
    !!----
    !!
    Elemental Subroutine Assign_Real_Rational_Dp(X, Res)
       !---- Arguments ----!
       type(rational), intent(in)   :: res
       real(kind=dp),  intent (out) :: x
      
       x=real(res%numerator,kind=dp)/real(res%denominator,kind=dp)
    
       return   
    End Subroutine Assign_Real_Rational_Dp

    !!----
    !!---- FUNCTION RATIONAL_ADD 
    !!----
    !!----
    !!
    Elemental Function Rational_Add(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       type(rational)              :: res
      
       res = r%numerator * s%denominator + r%denominator * s%numerator // &
           & r%denominator * s%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Add

    !!----
    !!---- FUNCTION INTEGER_RATIONAL_ADD 
    !!----
    !!----
    !!
    Elemental Function Integer_Rational_Add(I, S) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res
      
       res = (i * s%denominator +  s%numerator) // s%denominator
       res=rational_simplify(res)
       
       return
    End Function Integer_Rational_Add

    !!----
    !!---- FUNCTION RATIONAL_INTEGER_ADD 
    !!----
    !!----
    !!
    Elemental Function Rational_Integer_Add(S, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: s
       integer(kind=il),intent (in) :: i
       type(rational)               :: res
      
       res = (i * s%denominator +  s%numerator) // s%denominator
       res=rational_simplify(res)
    
       return   
    End Function Rational_Integer_Add

    !!----
    !!---- FUNCTION RATIONAL_MINUS 
    !!----
    !!----
    !!
    Elemental Function Rational_Minus(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational)              :: res
      
       res = - r%numerator // r%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Minus
    
    !!----
    !!---- FUNCTION RATIONAL_SUBTRACT 
    !!----
    !!----
    !!
    Elemental Function Rational_Subtract(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       type(rational)              :: res
      
       res = r%numerator * s%denominator - r%denominator * s%numerator // &
           & r%denominator * s%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Subtract

    !!----
    !!---- FUNCTION INTEGER_RATIONAL_SUBTRACT 
    !!----
    !!----
    !!
    Elemental Function Integer_Rational_Subtract(I, S) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res
      
       res = (i * s%denominator -  s%numerator) //  s%denominator
       res=rational_simplify(res)
       
       return
    End Function Integer_Rational_Subtract

    !!----
    !!---- FUNCTION RATIONAL_INTEGER_SUBTRACT 
    !!----
    !!----
    !!
    Elemental Function Rational_Integer_Subtract(S,I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: s
       integer(kind=il),intent (in) :: i
       type(rational)               :: res
      
       res = (s%numerator - i * s%denominator) //  s%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Integer_Subtract

    !!----
    !!---- FUNCTION RATIONAL_MULTIPLY 
    !!----
    !!----
    !!
    Elemental Function Rational_Multiply(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       type(rational)              :: res
      
       res = r%numerator * s%numerator // (r%denominator * s%denominator)
       res=rational_simplify(res)
       
       return
    End Function Rational_Multiply

    !!----
    !!---- FUNCTION RATIONAL_MULTIPLY 
    !!----
    !!----
    !!
    Elemental Function Integer_Rational_Multiply(I, S) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res
      
       res = i * s%numerator // s%denominator
       res=rational_simplify(res)
       
       return
    End Function Integer_Rational_Multiply

    !!----
    !!---- FUNCTION RATIONAL_INTEGER_MULTIPLY 
    !!----
    !!----
    !!
    Elemental Function Rational_Integer_Multiply(S,I) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: s
       type(rational)               :: res
      
       res = i * s%numerator // s%denominator
       res=rational_simplify(res)
       
       return
    End Function Rational_Integer_Multiply

    !!----
    !!---- FUNCTION RATIONAL_DIVIDE 
    !!----
    !!----
    !!
    Elemental Function Rational_Divide(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       type(rational)              :: res
       
       !---- Local Variables ----!
       integer(kind=il) :: denom
      
       denom = r%denominator * s%numerator
       if (denom /= 0_il) then
          res = r%numerator * s%denominator // denom
          res=rational_simplify(res)
       else
          res=0_il//1_il
       end if
       
       return
    End Function Rational_Divide
    
    !!----
    !!---- FUNCTION RATIONAL_DIVIDE_INT 
    !!----
    !!----
    !!
    Elemental Function Rational_Divide_Int(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       type(rational)               :: res
      
       if (i /= 0_il) then
          res = r%numerator // (r%numerator*i)
          res=rational_simplify(res)
       else
          res=0_il//1_il
       end if
       
       return
    End Function Rational_Divide_Int

    !!----
    !!---- FUNCTION RATIONAL_LT 
    !!----
    !!----
    !!
    Pure Function Rational_Lt(R, S) Result(Res)
       !---- Arguments ----! 
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
       
       !---- Local Variables ----!
       type(rational) :: r_simple, s_simple

       r_simple = rational_simplify(r)
       s_simple = rational_simplify(s)
      
       res = r_simple%numerator * s_simple%denominator < &
           & s_simple%numerator * r_simple%denominator
           
       return    
    End Function Rational_Lt

    !!----
    !!---- FUNCTION RATIONAL_LT_INTEGER 
    !!----
    !!----
    !!
    Pure Function Rational_Lt_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       logical                      :: res
       
       !---- Local Variables ----!
       type(rational) :: r_simple

       r_simple = rational_simplify(r)
       res = r_simple%numerator < i * r_simple%denominator
      
       return
    End Function Rational_Lt_Integer

    !!----
    !!---- FUNCTION INTEGER_LT_RATIONAL 
    !!----
    !!----
    !!
    Pure Function Integer_Lt_Rational (I,R) Result (Res)
      integer(kind=il),intent (in) :: i
      type(rational), intent (in) :: r
      logical                      :: res
      type(rational) :: r_simple
      r_simple = rational_simplify (r)
      res = i * r_simple % denominator < r_simple % numerator
    End Function Integer_Lt_Rational

    !!----
    !!---- FUNCTION RATIONAL_LE 
    !!----
    !!----
    !!
    Pure Function Rational_LE(R, S) Result(Res)
       !---- Argument ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
       
       !---- Local Variables ----!
       type(rational) :: r_simple, s_simple

       r_simple = rational_simplify(r)
       s_simple = rational_simplify(s)
      
       res = r_simple%numerator * s_simple%denominator <= &
           & s_simple%numerator * r_simple%denominator
           
       return    
    End Function Rational_LE

    !!----
    !!---- FUNCTION RATIONAL_LE_INTEGER 
    !!----
    !!----
    !!
    Pure Function Rational_Le_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       logical                      :: res
       
       !---- Local Variables ----!
       type(rational) :: r_simple
      
       r_simple = rational_simplify(r)
       res = r_simple%numerator <= i * r_simple%denominator
       
       return
    End Function Rational_Le_Integer

    !!----
    !!---- FUNCTION INTEGER_LE_RATIONAL 
    !!----
    !!----
    !!
    Pure Function Integer_Le_Rational(I, R) Result(Res)
       !---- Arguments ----! 
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
      
       r_simple = rational_simplify(r)
       res = i * r_simple%denominator  <= r_simple%numerator
       
       return
    End Function Integer_Le_Rational

    !!----
    !!---- FUNCTION RATIONAL_GT 
    !!----
    !!----
    !!
    Pure Function Rational_GT(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple, s_simple

       r_simple = rational_simplify(r)
       s_simple = rational_simplify(s)
      
       res = r_simple%numerator * s_simple%denominator > &
           & s_simple%numerator * r_simple%denominator
           
       return    
    End Function Rational_GT

    !!----
    !!---- FUNCTION RATIONAL_GT_INTEGER 
    !!----
    !!----
    !!
    Pure Function Rational_Gt_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),   intent (in) :: r
       integer(kind=il), intent (in) :: i
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
      
       r_simple = rational_simplify(r)
      
       res = r_simple%numerator > i * r_simple%denominator
       
       return
    End Function Rational_Gt_Integer

    !!----
    !!---- FUNCTION INTEGER_GT_RATIONAL 
    !!----
    !!----
    !!
    Pure Function Integer_Gt_Rational(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
      
       r_simple = rational_simplify(r)
       res = i * r_simple%denominator > r_simple%numerator
       
       return
    End Function Integer_Gt_Rational
    
    !!----
    !!---- FUNCTION RATIONAL_GE 
    !!----
    !!----
    !!
    Pure Function Rational_GE(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple, s_simple
 
       r_simple = rational_simplify(r)
       s_simple = rational_simplify(s)
      
       res = r_simple%numerator * s_simple%denominator >= &
           & s_simple%numerator * r_simple%denominator
           
       return    
    End Function Rational_GE

    !!----
    !!---- FUNCTION RATIONAL_GE_INTEGER 
    !!----
    !!----
    !!
    Pure Function Rational_GE_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
      
       r_simple = rational_simplify(r)
       res = r_simple%numerator >= i * r_simple%denominator
       
       return
    End Function Rational_GE_Integer

    !!----
    !!---- FUNCTION INTEGER_GE_RATIONAL
    !!----
    !!----
    !!
    Pure Function Integer_GE_Rational(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
      
       !---- Local Variables ----!
       type(rational) :: r_simple
       
       r_simple = rational_simplify(r)
       res = i * r_simple%denominator >= r_simple%numerator
       
       return
    End Function Integer_GE_Rational

    !!----
    !!---- FUNCTION RATIONAL_EQ
    !!----
    !!----
    !!
    Elemental Function Rational_EQ(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
      
       res = r%numerator * s%denominator == s%numerator * r%denominator
       
       return
    End Function Rational_EQ

    !!----
    !!---- FUNCTION RATIONAL_EQ_INTEGER
    !!----
    !!----
    !!
    Elemental Function Rational_EQ_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       logical                      :: res
      
       res = r%denominator == 1_il .and. r%numerator == i
       
       return
    End Function Rational_EQ_Integer
    
    !!----
    !!---- FUNCTION INTEGER_EQ_RATIONAL
    !!----
    !!----
    !!
    Elemental Function Integer_Eq_Rational(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
       
       res = r%denominator == 1_il .and. r%numerator == i
       
       return
    End Function Integer_EQ_Rational

    !!----
    !!---- FUNCTION RATIONAL_NE
    !!----
    !!----
    !!
    Pure Function Rational_Ne(R, S) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational), intent (in) :: s
       logical                     :: res
      
       res = r%numerator * s%denominator /= s%numerator * r%denominator
       
       return
    End Function Rational_Ne

    !!----
    !!---- FUNCTION RATIONAL_NE_INTEGER
    !!----
    !!----
    !!
    Elemental Function Rational_Ne_Integer(R, I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent (in) :: i
       logical                      :: res
      
       res = r%numerator /= i * r%denominator
       
       return
    End Function Rational_Ne_Integer

    !!----
    !!---- FUNCTION INTEGER_NE_RATIONAL
    !!----
    !!----
    !!
    Elemental Function Integer_Ne_Rational(I, R) Result(Res)
       !---- Arguments ----!
       integer(kind=il),intent (in) :: i
       type(rational),  intent (in) :: r
       logical                      :: res
      
       res = r%numerator /= i * r%denominator
       
       return
    End Function Integer_Ne_Rational

    !!----
    !!---- FUNCTION RATIONAL_ABS
    !!----
    !!----
    !!
    Elemental Function Rational_Abs(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       type(rational)              :: res
       
       res = sign (r%numerator, r%denominator) // r%denominator
       
       return
    End Function Rational_Abs

    !!----
    !!---- FUNCTION RATIONAL_INT
    !!----
    !!----
    !!
    Elemental Function Rational_Int(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       integer(kind=il)            :: res
       
       res = r%numerator / r%denominator
       
       return
    End Function Rational_Int

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
    !!---- FUNCTION NINT_RATIONAL
    !!----
    !!----
    !!
    Elemental Function Nint_Rational(R) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in)  :: r
       integer(kind=il)              :: res
      
       res = nint(real(r%numerator,kind=dp)/real(r%denominator,kind=dp),kind=il)
       
       return
    End Function Nint_Rational

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
    !!---- FUNCTION RATIONAL_MODULO
    !!----
    !!----
    !!
    Elemental Function Rational_Modulo(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       integer                     :: res
      
       res = modulo (r%numerator, r%denominator)
       
       return
    End Function Rational_Modulo

    !!----
    !!---- FUNCTION RATIONAL_MODULO_INT
    !!----
    !!----
    !!
    Elemental Function Rational_Modulo_Int(R,I) Result(Res)
       !---- Arguments ----!
       type(rational),  intent (in) :: r
       integer(kind=il),intent(in)  :: i
       type(rational)               :: res
      
       !---- Local Variables ----!
       real(kind=dp) :: val
      
       val = modulo (real(r%numerator,kind=dp) / real(r%denominator,kind=dp),real(i,kind=dp))
       res = val
       
       return
    End Function Rational_Modulo_Int

    !!----
    !!---- FUNCTION RATIONAL_MODULO
    !!----
    !!----
    !!
    Elemental Function Rational_Mod(R) Result(Res)
       !---- Arguments ----!
       type(rational), intent (in) :: r
       integer                     :: res
      
       res = mod (r%numerator, r%denominator)
       
       return
    End Function Rational_Mod

    !!----
    !!---- FUNCTION RATIONAL_MODULO_INT
    !!----
    !!----
    !!
    Elemental Function Rational_Mod_Int(R,I) Result(Res)
       !---- Arguments ----!
      type(rational),   intent (in) :: r
      integer(kind=il), intent (in) :: i
      type(rational)                :: res
      
      !---- Local Variables ----!
      real(kind=dp) :: val
      
      val = mod(real(r%numerator,kind=dp) / real(r%denominator,kind=dp),real(i,kind=dp))
      res = val
      
      return
    End Function Rational_Mod_Int

    !!----
    !!---- FUNCTION RATIONAL_DOT_PRODUCT
    !!----
    !!----
    !!
    Function Rational_Dot_Product(R1,R2) Result(Res)
       !---- Arguments ----!
       type(rational), dimension(:), intent (in) :: r1
       type(rational), dimension(:), intent (in) :: r2
       type(rational)                            :: res
      
       !---- Local Variables ----!
       integer :: n1,n2,i
      
       n1=size(r1); n2=size(r2)
       res=0_il//1_il
      
       if (n1 == n2) then
          do i=1,n1
            res = rational_simplify(res + r1(i)*r2(i))
          end do
       else
         Err_CFML=.true.
         Err_CFML_Flag=2
         write(unit=Err_CFML_Msg,fmt="(a,2i4)") &
              "Error in DOT_PRODUCT: the dimensions of the arguments are different",n1,n2
       end if
       
       return
    End Function Rational_Dot_Product

    !!----
    !!---- FUNCTION RATIONAL_MATMUL_MATVEC
    !!----
    !!----
    !!
    Function Rational_Matmul_Matvec(Mat,Vec) Result(Vec_Out)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent (in) :: mat
       type(rational), dimension(:),   intent (in) :: vec
       type(rational), dimension(size(vec))        :: vec_out
       
       !---- Local Variables ----!
       integer :: n1,n2,n3,i
 
       n1=size(mat,dim=1); n2=size(mat,dim=2); n3=size(vec)
       if (n1 == n2 .and. n2 == n3) then
          do i=1,n3
             vec_out(i) = rational_simplify(dot_product(mat(i,:),vec))
          end do
       else
          Err_CFML=.true.
          Err_CFML_Flag=0
          write(unit=Err_CFML_Msg,fmt="(a,3i4)") &
               "Error in MATMUL: the dimensions of the arguments are non-conformable",n1,n2,n3
       end if
       
       return
    End Function Rational_Matmul_Matvec

    !!----
    !!---- FUNCTION RATIONAL_MATMUL_MATMAT
    !!----
    !!----
    !!
    Function Rational_Matmul_Matmat(Mat1,Mat2) Result(Mat_Out)
       !---- Arguments ----! 
       type(rational), dimension(:,:), intent (in)                 :: mat1
       type(rational), dimension(:,:), intent (in)                 :: mat2
       type(rational),dimension(size(mat1,dim=1),size(mat2,dim=2)) :: mat_out

       !---- Local Variables ----!
       integer :: n1,n2,n3,n4,i,j
      
       n1=size(mat1,dim=1); n2=size(mat1,dim=2); n3=size(mat2,dim=1); n4=size(mat2,dim=2)
       if (n2 == n3) then
          do j=1,n4
             do i=1,n1
                mat_out(i,j) = rational_simplify(dot_product(mat1(i,:),mat2(:,j)))
             end do
          end do
       else
          Err_CFML=.true.
          Err_CFML_Flag=2
          write(unit=Err_CFML_Msg,fmt="(a,4i4)") &
               "Error in MATMUL: the dimensions of the arguments are non-conformable",n1,n2,n3,n4
       end if
       
       return
    End Function Rational_Matmul_Matmat

    !!----
    !!---- FUNCTION RATIONAL_MAXVAL_VECTOR
    !!----
    !!----
    !!
    Function Rational_Maxval_Vector(R) Result(Res)
       !---- Arguments ----!
       type(rational), dimension(:), intent (in) :: r
       type(rational)                            :: res
       
       !---- Local Variables ----!
       integer :: i,n
      
       n=size(r)
       res=-huge(1_il)//1_il
       do i=1,n
          if (r(i) > res) res=r(i)
       end do
       
       return
    End Function Rational_Maxval_Vector

    !!----
    !!---- FUNCTION RATIONAL_MAXVAL_MATRIX
    !!----
    !!----
    !!
    Function Rational_Maxval_Matrix(R) Result(Res)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent (in) :: r
       type(rational)                              :: res
      
       !---- Local Variables ----!
       integer :: i,j,n1,n2
      
       n1=size(r,dim=1);  n2=size(r,dim=2)
       res=-huge(1_il)//1_il
       do j=1,n2
          do i=1,n1
             if (r(i,j) > res) res=r(i,j)
          end do
       end do
       
       return
    End Function Rational_Maxval_Matrix

    !!----
    !!---- FUNCTION RATIONAL_MINVAL_VECTOR
    !!----
    !!----
    !!
    Function Rational_Minval_Vector(R) Result(Res)
       !---- Arguments ----!
       type(rational), dimension(:), intent (in) :: r
       type(rational)                            :: res
       
       !---- Local Variables ----!
       integer :: i,n
       
       n=size(r)
       res=huge(1_il)//1_il
       do i=1,n
          if (r(i) < res) res=r(i)
       end do
       
       return
    End Function Rational_Minval_Vector

    !!----
    !!---- FUNCTION RATIONAL_MINVAL_MATRIX
    !!----
    !!----
    !!
    Function Rational_Minval_Matrix(R) Result(Res)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent (in) :: r
       type(rational)                              :: res
       
       !---- Local Variables ----!
       integer :: i,j,n1,n2
      
       n1=size(r,dim=1);  n2=size(r,dim=2)
       res=huge(1_il)//1_il
       do j=1,n2
          do i=1,n1
             if (r(i,j) < res) res=r(i,j)
          end do
       end do
       
       return
    End Function Rational_Minval_Matrix

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
          Err_CFML=.true.
          Err_CFML_Flag=2
          write(unit=Err_CFML_Msg,fmt="(a)") &
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
             Err_CFML=.true.
             Err_CFML_Flag=2
             write(unit=Err_CFML_Msg,fmt="(a)") "Singular Matrix!"
          else
             invMat=invA
          end if
       else
          Err_CFML=.true.
          Err_CFML_Flag=2
          write(unit=Err_CFML_Msg,fmt="(a)") &
               "Error in Determinant: the provided matrix is not square!"
       end if
       
       return
    End Subroutine Rational_Inv_Matrix

    !!----
    !!---- FUNCTION RATIONAL_MAXLOC_MAT
    !!----
    !!----
    !!
    Pure Function Rational_Maxloc_Mat(Mat) Result(Pos_Max)
       !---- Arguments ----!
       type(rational),  dimension(:,:), intent(in) :: Mat
       integer(kind=il),dimension(2)               :: pos_max
      
       !---- Local variables ----!
       integer        :: nu1,nl1,nu2,nl2,i,j
       type(rational) :: res
      
       nu1=ubound(Mat,dim=1); nu2=ubound(Mat,dim=2)
       nl1=lbound(Mat,dim=1); nl2=lbound(Mat,dim=2)

       res=-huge(1_il)//1_il
       do j=nl2,nu2
          do i=nl1,nu1
             if (mat(i,j) > res) then
                res=mat(i,j)
                pos_max=[i,j]
             end if
          end do
       end do
       
       return
    End Function Rational_Maxloc_Mat

    !!----
    !!---- FUNCTION RATIONAL_MAXLOC_VECT
    !!----
    !!----
    !!
    Pure Function Rational_Maxloc_Vect(Vec) Result(Pos_Max)
       !---- Arguments ----!
       type(rational), dimension(:), intent(in) :: vec
       integer                                  :: pos_max
      
      
       !---- Local variables ----!
       integer:: nu,nl,i
       type(rational) :: res
      
       nu=ubound(vec,dim=1)
       nl=lbound(vec,dim=1)
       res=-huge(1_il)//1_il
       do i=nl,nu
          if (vec(i) > res) then
             res=vec(i)
             pos_max=i
          end if
       end do
       
       return
    End Function Rational_Maxloc_Vect

    !!----
    !!---- FUNCTION RATIONAL_MINLOC_MAT
    !!----
    !!----
    !!
    Pure Function Rational_Minloc_Mat(Mat) Result(Pos_Min)
       !---- Arguments ----!
       type(rational),  dimension(:,:), intent(in) :: Mat
       integer(kind=il),dimension(2)               :: pos_min
      
       !---- Local variables ----!
       integer:: nu1,nl1,nu2,nl2,i,j
       type(rational) :: res
      
       nu1=ubound(Mat,dim=1); nu2=ubound(Mat,dim=2)
       nl1=lbound(Mat,dim=1); nl2=lbound(Mat,dim=2)

       res=huge(1_il)//1_il
       do j=nl2,nu2
          do i=nl1,nu1
             if (mat(i,j) < res) then
                res=mat(i,j)
                pos_min=[i,j]
             end if
          end do
       end do
       
       return
    End Function Rational_Minloc_Mat

    !!----
    !!---- FUNCTION RATIONAL_MAXLOC_VECT
    !!----
    !!----
    !!
    Pure Function Rational_Minloc_Vect(Vec) Result(Pos_Min)
       !---- Arguments ----!
       type(rational), dimension(:), intent(in) :: vec
       integer                                  :: pos_min
      
       !---- Local variables ----!
       integer        :: nu,nl,i
       type(rational) :: res
      
       nu=ubound(vec,dim=1)
       nl=lbound(vec,dim=1)
       res=huge(1_il)//1_il
       do i=nl,nu
          if (vec(i) < res) then
             res=vec(i)
             pos_min=i
          end if
       end do
       
       return
    End Function Rational_Minloc_Vect

    !!----
    !!---- FUNCTION RATIONAL_SUM_VEC
    !!----
    !!----
    !!
    Pure Function Rational_Sum_Vec(Vec) Result(Suma)
       !---- Arguments ----!
       type(rational), dimension(:), intent(in) :: vec
       type(rational)                           :: suma
      
       !---- Local variables ----!
       real(kind=dp), dimension(size(vec)) :: rvec
      
       rvec=vec
       suma=sum(rvec)
       suma=rational_simplify(suma)
       
       return
    End Function Rational_Sum_Vec
    
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