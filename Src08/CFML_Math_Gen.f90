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
!!---- MODULE: CFML_Math_General
!!----   INFO: Mathematic general utilities for use in Crystallography and
!!----         Solid State Physics and Chemistry.
!!----
!!
 Module CFML_Math_General
    !---- Use Modules ----!
    Use CFML_GlobalDeps
    Use CFML_DefPar, only : Err_MathGen, Err_MathGen_Mess, Primes

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public :: Factorial, Factorial_SP, Factorial_DP, Negligible, Co_Linear, Co_prime, &
              Debye, Equal_Matrix, Equal_Vector, Euclidean_Norm, In_limits, Locate,   &
              Lower_Triangular, Modulo_Lat, Norm, Outerprod, Pgcd, Ppcm, Pythag,      &
              Scalar, Trace, Upper_Triangular, Zbelong, Determinant, MatInv,          &
              Linear_Dependent, Splint, Poly_Legendre

    !---- List of public subroutines ----!
    public :: Co_Prime_vector, Diagonalize_SH, First_Derivative, In_Sort,             &
              Init_Err_Mathgen, Invert_Matrix, LU_Decomp, LU_Backsub,                 &
              Sort_Strings, Median_QS, Points_in_Line2D, Rank, Second_Derivative,     &
              Set_Epsg, SmoothingVec, Sort, Spline, Svdcmp, Swap, Sph_Jn


    !---- Definitions ----!
    real(kind=cp), parameter :: EP_SS=1.0E-12_cp  ! Internal epsilon value used for comparison in matrix operations

    real(kind=cp)            :: epss=1.0E-5_cp    ! Internal epsilon value used for comparing reals to integers


    !---- Interfaces - Overloaded ----!
    Interface  Debye
       Module Procedure Debye_DP
       Module Procedure Debye_SP
    End Interface

    Interface  Negligible
       Module Procedure Negligibler
       Module Procedure Negligiblec
    End Interface

    Interface  Co_Linear
       Module Procedure Co_linear_C
       Module Procedure Co_linear_I
       Module Procedure Co_linear_R
    End Interface

    Interface  Determinant
       Module Procedure Determinant_c
       Module Procedure Determinant_r
    End Interface

    Interface  Equal_Matrix
       Module Procedure Equal_Matrix_I
       Module Procedure Equal_Matrix_R
    End Interface

    Interface  Equal_Vector
       Module Procedure Equal_Vector_I
       Module Procedure Equal_Vector_R
    End Interface

    Interface In_Limits
       Module Procedure In_Limits_dp
       Module Procedure In_Limits_int
       Module Procedure In_Limits_sp
    End Interface

    Interface  Linear_Dependent
       Module Procedure Linear_Dependentc
       Module Procedure Linear_Dependenti
       Module Procedure Linear_Dependentr
    End Interface

    Interface  Locate
       Module Procedure Locate_I
       Module Procedure Locate_R
       Module Procedure Locate_Ib
       Module Procedure Locate_Rb
    End Interface

    Interface  Lower_Triangular
       Module Procedure Lower_Triangular_I
       Module Procedure Lower_Triangular_R
    End Interface

    Interface Norm
       Module Procedure Norm_I
       Module Procedure Norm_R
    End Interface Norm

    Interface  Outerprod
       Module Procedure Outerprod_dp
       Module Procedure Outerprod_sp
    End Interface

    Interface  Pythag
       Module Procedure Pythag_dp
       Module Procedure Pythag_sp
    End Interface

    Interface Scalar
       Module Procedure Scalar_I
       Module Procedure Scalar_R
    End Interface Scalar

    Interface  Trace
       Module Procedure Trace_C
       Module Procedure Trace_I
       Module Procedure Trace_R
    End Interface

    Interface  Upper_Triangular
       Module Procedure Upper_Triangular_I
       Module Procedure Upper_Triangular_R
    End Interface

    Interface  Zbelong
       Module Procedure ZbelongM
       Module Procedure ZbelongN
       Module Procedure ZbelongV
    End Interface

    Interface  Diagonalize_SH
       Module Procedure Diagonalize_HERM
       Module Procedure Diagonalize_SYMM
    End Interface

    Interface  Rank
       Module Procedure Rank_dp
       Module Procedure Rank_sp
    End Interface

    Interface  Sort
       Module Procedure Sort_I
       Module Procedure Sort_R
    End Interface

    Interface  Svdcmp
       Module Procedure Svdcmp_dp
       Module Procedure Svdcmp_sp
    End Interface

    Interface Swap
        Module Procedure swap_c
        Module Procedure swap_cm
        Module Procedure swap_cv
        Module Procedure swap_i
        Module Procedure swap_im
        Module Procedure swap_iv
        Module Procedure swap_r
        Module Procedure swap_rm
        Module Procedure swap_rv
        Module Procedure masked_swap_r
        Module Procedure masked_swap_rm
        Module Procedure masked_swap_rv
    End interface

 Contains

    !---- Functions ----!

    !!----
    !!---- Elemental Function Factorial(n) Result(fact)
    !!----
    !!----    Factorial of N.
    !!----    This function works fine only for N <= 12 (Integer *4)
    !!----
    !!---- Update: January - 2017
    !!
    Elemental Function Factorial(n) Result(fact)
       !---- Argument ----!
       integer, intent(in) :: n      ! Argument
       integer             :: fact

       !---- Local variables ----!
       integer, parameter :: N_LIM = 12
       integer            :: i

       !> Check the current limits
       if ( n > N_LIM) then
          fact=0
          return
       end if

       if (n ==0) then
          fact=1
       else
          fact=1
          do i=1,n
             fact=fact*i
          end do
       end if

       return
    End Function Factorial

    !!----
    !!---- Elemental Function Factorial_DP(N) Result(Res)
    !!----
    !!----    Factorial of N but the value returned is a double number
    !!----
    !!---- Update: January - 2017
    !!
    Elemental Function Factorial_DP(N) Result(Res)
       !---- Arguments ----!
       integer, intent(in) :: N

       !---- Local Variables ----!
       integer       :: i
       real(kind=dp) :: Res

       if (n == 0) then
          res = 1.0_dp
       else
          res=1.0_dp
          do i=1,n
             res = res * dble(i)
          end do
       end if

       return
    End Function Factorial_DP

    !!----
    !!---- Elemental Function Factorial_SP(N) Result(Res)
    !!----
    !!----    Factorial of N but the value returned is a real number
    !!----
    !!---- Update: January - 2017
    !!
    Elemental Function Factorial_SP(N) Result(Res)
       !---- Arguments ----!
       integer, intent(in) :: N

       !---- Local Variables ----!
       integer       :: i
       real(kind=sp) :: Res

       if (n == 0) then
          res = 1.0_sp
       else
          res=1.0_sp
          do i=1,n
             res = res * float(i)
          end do
       end if

       return
    End Function Factorial_SP

    !!--++
    !!--++ Elemental Function Negligiblec(X)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if a complex number is negligible
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Negligiblec(X) Result(Neglig)
       !---- Argument ----!
       complex, intent( in) :: X        ! Argument
       logical              :: Neglig

       Neglig=.false.
       if (abs(X) > epss) return

       Neglig=.true.

       return
    End Function Negligiblec

    !!--++
    !!--++ Elemental Function Negligibler(X)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is negligible (abs < EPSS)
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Negligibler(X) Result(neglig)
       !---- Argument ----!
       real(kind=cp), intent( in) :: X     ! Vector
       logical                    :: Neglig

       Neglig=.false.
       if (abs(X) > epss) return

       Neglig=.true.

       return
    End Function Negligibler

    !!----
    !!---- FUNCTION PLGNDR
    !!----
    !!----    Compute the Associated Legendre Polynomial Pml(x).
    !!----
    !!---- Here m (order) and l (degree) are integers satisfying
    !!----    0 <= m <= l
    !!----   -1 <= x <= 1.
    !!----
    !!---- Update: 11/07/2015
    !!
    Elemental Function Poly_Legendre(l,m,x) result(plmx)
       !---- Arguments ----!
       integer,      intent (in) :: l
       integer,      intent (in) :: m
       real(kind=cp),intent (in) :: x
       real(kind=cp)             :: plmx

       !---- Local variables ----!
       integer       :: i, ll
       real(kind=cp) :: fact, pll, pmm, pmmp1, somx2


       !> Check limits
       if (m < 0 .or. m > l .or. abs(x) > 1.0) then
          plmx=0.0
          return
       end if

       !> Calculation
       pmm=1.0
       if (m > 0) then                  !Compute P(m,m)
          somx2=sqrt((1.0_cp-x)*(1.0_cp+x))
          fact=1.0_cp
          do i=1,m
             pmm=-pmm*fact*somx2
             fact=fact+2.0_cp
          end do
       end if

       if (l == m) then
          plmx=pmm
       else
          pmmp1=x*real(2*m+1,kind=cp)*pmm           !Compute P(m,m+1)
          if (l == m+1) then
             plmx=pmmp1                 !Compute P(m,L), L > m+1
          else
             do ll=m+2,l
                pll=(x*real(2*ll-1,kind=cp)*pmmp1-real(ll+m-1,kind=cp)*pmm)/real(ll-m,kind=cp)
                pmm=pmmp1
                pmmp1=pll
             end do
             plmx=pll
          end if
       end if

       return
    End Function Poly_Legendre

    !!--++
    !!--++ Function CHEVAL
    !!--++
    !!--++ PRIVATE (USED FOR DEBYE FUNCTIONS)
    !!--++ This function evaluates a Chebyshev series, using the Clenshaw method
    !!--++ with Reinsch modification, as analysed in the paper by Oliver in
    !!--++ J.I.M.A., vol. 20, 1977, pp379-391
    !!--++
    !!--++  Update:  January - 2017
    !!
    Function Cheval(n, a, t) Result(fval)
       !---- Arguments ----!
       integer,                       intent(in) :: N  ! The no. of terms in the sequence
       real(kind=dp), dimension(0:N), intent(in) :: A  ! The coefficients of the Chebyshev series
       real(kind=dp),                 intent(in) :: T  ! The value at which the series is to be evaluated

       !---- Local Variables ----!
       integer       :: i
       real(kind=dp) :: fval
       real(kind=dp) :: d1, d2, tt, u0, u1, u2
       real(kind=dp), parameter  :: HALF = 0.5_dp
       real(kind=dp), parameter  :: TWO  = 2.0_dp
       real(kind=dp), parameter  :: TEST = 0.6_dp


       !> Init
       u1 = 0.0_dp

       !>  If ABS ( T )  < 0.6 use the standard Clenshaw method
       if (abs(t) < test) then
          u0 = 0.0_dp
          tt = t + t
          do i = n, 0, -1
             u2 = u1
             u1 = u0
             u0 = tt * u1 + a(i) - u2
          end do
          fval = (u0-u2) / two

       else
          !> If ABS ( T )  > =  0.6 use the Reinsch modification
          d1 = 0.0_dp

          !>  T > =  0.6 code
          if (t > 0.0_dp) then
             tt = (t-half) - half
             tt = tt + tt
             do i = n, 0, -1
                d2 = d1
                u2 = u1
                d1 = tt * u2 + a(i) + d2
                u1 = d1 + u2
             end do
             fval = (d1+d2) / two

          else
             !> T < =  -0.6 code
             tt = (t+half) + half
             tt = tt + tt
             do i = n, 0, -1
                d2 = d1
                u2 = u1
                d1 = tt * u2 + a(i) - d2
                u1 = d1 - u2
             end do
             fval = (d1-d2) / two
          end if
       end if

       return
    End Function Cheval

    !!--++
    !!--++ Logical Function Co_Linear_C(Vec1, Vec2, N)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two complex vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_C(Vec1,Vec2,n) Result(co_linear)
       !---- Argument ----!
       complex, dimension(:), intent(in) :: Vec1  ! Complex vector
       complex, dimension(:), intent(in) :: Vec2  ! Complex vector
       integer,               intent(in) :: n     ! Dimension of the vector
       logical                           :: co_linear

       !---- Local variables ----!
       integer :: i,ia,ib
       complex :: c

       co_linear=.true.
       do i=1,n
          if (abs(Vec1(i)) > epss) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(Vec2(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=Vec1(ia)/Vec2(ib)
          do i=1,n
             if (abs(Vec1(i)-c*Vec2(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_C

    !!--++
    !!--++ Logical Function Co_Linear_I(Vec1, Vec2, N)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are co-linear
    !!--++
    !!--++ Update: October - 2008
    !!
    Function Co_linear_I(Vec1,Vec2,n) Result(co_linear)
       !---- Argument ----!
       integer, dimension(:), intent(in) :: Vec1  ! Integer vector
       integer, dimension(:), intent(in) :: Vec2  ! Integer vector
       integer,               intent(in) :: n     ! Dimension of the vector
       logical                           :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib
       real(kind=cp) :: c

       co_linear=.true.
       do i=1,n
          if (abs(vec1(i)) > 0) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(vec2(i)) > 0) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=real(vec1(ia))/real(vec2(ib))
          do i=1,n
             if (abs( real(vec1(i))-c*real(vec2(i)) ) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_I

    !!--++
    !!--++ Logical Function Co_Linear_R(Vec1, Vec2, N)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_R(Vec1,Vec2,n) Result(co_linear)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in) :: Vec1    ! Real vector
       real(kind=cp), dimension(:), intent(in) :: Vec2    ! Real vector
       integer,                     intent(in) :: n    ! Dimension of the vector
       logical                                 :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib
       real(kind=cp) :: c

       co_linear=.true.
       do i=1,n
          if (abs(vec1(i)) > epss) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(vec2(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=vec1(ia)/vec2(ib)
          do i=1,n
             if (abs(vec1(i)-c*vec2(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_R

    !!----
    !!---- Logical Function Co_Prime(Vec,Prime_max) result(cop)
    !!----
    !!---- Provides the value .TRUE. if the array Vec contains co-prime
    !!---- integers: there is no common divisor for all the integers.
    !!----
    !!---- Only the first 1000 prime numbers are stored in the module array "primes"
    !!---- Prime_max is the maximum prime number to be tested. It is calculated if not given.
    !!----
    !!
    Function Co_Prime(Vec,Prime_max) result(cop)
      integer, dimension(:), intent(in) :: Vec           ! Vector containing co-prime integers
      integer, optional,     intent(in) :: Prime_max     ! Maximum prime number to be tested
      Logical                           :: cop
      !---- Local variables ----!
      integer :: i,j,im,k,dimv,imaxv,maxv

      !> Init
      cop=.true.

      maxv=maxval(abs(vec))
      if (present(Prime_max)) then
         imaxv=Prime_max
      else
         imaxv=maxv
      end if

      !> If the maximum value of the indices is 1 they are coprimes
      if (maxv == 1) return
      if (maxv == 0) then
         cop=.false.
         return
      end if

      !> Search the maximum prime number to be tested
      if (imaxv > 7919) then
         im=1000
      else
         do i=1,1000
            if (imaxv > primes(i)) cycle
            im=i
            exit
         end do
      end if

      !---- Indices greater than 1
      dimv=size(vec)
      do_p: do i=1,im
         k=primes(i)
         do j=1,dimv
            if ( mod(vec(j),k) /= 0) cycle do_p
         end do
         cop=.false.
         exit
      end do do_p

      return
    End Function Co_Prime

    !!--++
    !!--++ FUNCTION DEBYE_DP(N,X)
    !!--++
    !!--++ OVERLOADED
    !!--++ Calculates the Debye function of order N
    !!--++
    !!--++ If X < 0.0 then limited to |x| < 2*Pi
    !!--++
    !!--++ Update: January 2017
    !!
    Function Debye_DP(N,X) Result(Fval)
       !---- Arguments ----!
       integer,       intent(in) :: N ! Order of the Debye function
       real(kind=dp), intent(in) :: X ! Argument Value

       !---- Local Variables ----!
       real(kind=dp) :: fval

       !> Init
       fval=0.0_dp

       !> Check
       if (n <= 0) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess="The order for Debye function was ZERO!"
          return
       end if

       if (x < 0.0_dp) then
          if (abs(x) > tpi) then
             ERR_MathGen=.true.
             ERR_MathGen_Mess="The argument is negative and less than 2Pi"
             return
          end if
          fval =DebyeN(n,x)
       else
          select case (n)
             case (1)
                fval=Debye1(x)
             case (2)
                fval=Debye2(x)
             case (3)
                fval=Debye3(x)
             case (4)
                fval=Debye4(x)
             case (5:)
                if (x > tpi) then
                   ERR_MathGen=.true.
                   ERR_MathGen_Mess="The argument was greater then 2Pi and the order >= 5!"
                  return
                end if
                fval=DebyeN(n,x)

          end select
       end if

       return
    End Function Debye_DP

    !!--++
    !!--++ FUNCTION DEBYE_SP(N,X)
    !!--++
    !!--++ OVERLOADED
    !!--++ Calculates the Debye function of order N
    !!--++
    !!--++ If X < 0.0 then limited to |x| < 2*Pi
    !!--++
    !!--++ Update: January 2017
    !!
    Function Debye_SP(N,X) Result(Fval)
       !---- Arguments ----!
       integer,       intent(in) :: N ! Order of the Debye function
       real(kind=sp), intent(in) :: X ! Argument value

       !---- Local Variables ----!
       real(kind=sp) :: fval
       real(kind=dp) :: xx,ff

       !> Init
       fval=0.0_sp
       xx=dble(x)
       ff=0.0_dp

       !> Check
       if (n <= 0) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess="The order for Debye function was ZERO!"
          return
       end if

       if (x < 0.0) then
          if (abs(x) > tpi) then
             ERR_MathGen=.true.
             ERR_MathGen_Mess="The argument is negative and less than 2Pi"
             return
          end if
          ff =DebyeN(n,xx)
       else
          select case (n)
             case (1)
                ff=Debye1(xx)
             case (2)
                ff=Debye2(xx)
             case (3)
                ff=Debye3(xx)
             case (4)
                ff=Debye4(xx)
             case (5:)
                if (x > tpi) then
                   ERR_MathGen=.true.
                   ERR_MathGen_Mess="The argument was greater then 2Pi and the order >= 5!"
                   return
                end if
                ff=DebyeN(n,xx)
          end select
       end if

       !> Result
       if (dble(huge(fval)) < ff) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess="The value is greater than huge value for real case! in Debye Function"
          return
       else
          fval=real(ff)
       end if

       return
    End Function Debye_SP


    !!--++
    !!--++ Function Debye1(x)
    !!--++
    !!--++ Calculates the Debye function of order 1, defined as
    !!--++    DEBYE1(x) = [Integral {0 to x} t/(exp(t)-1) dt] / x
    !!--++
    !!--++ The code uses Chebyshev series whose coefficients are given to 20
    !!--++ decimal places.
    !!--++
    !!--.. EXTRA INFORMATION
    !!--..
    !!--.. If X < 0.0 an error message is defined and the function returns the value 0.0
    !!--..
    !!--.. NTERMS: The no. of elements of the array ADEB1. The recommended value is such that
    !!--..         ABS(ADEB1(NTERMS)) < EPS/100 , with 1 <= NTERMS <= 18
    !!--..
    !!--..   XLOW: The value below which DEBYE1 = 1 - x/4 + x*x/36 to machine precision.
    !!--..         The recommended value is SQRT(8*EPSNEG)
    !!--..
    !!--.. XUPPER: The value above which DEBYE1 = (pi*pi/(6*x)) - exp(-x)(x+1)/x.
    !!--..         The recommended value is -LOG(2*EPS)
    !!--..
    !!--..   XLIM: The value above which DEBYE1 = pi*pi/(6*x)
    !!--..         The recommended value is -LOG(XMIN)
    !!--++
    !!--++ Update: January 2017
    !!
    Function Debye1(X) Result(Fval)
       !---- Arguments ----!
       real(kind=dp), intent(in) :: X  ! Argument

       !---- Local Variables ----!
       integer       :: i, nexp, nterms
       real(kind=dp) :: fval
       real(kind=dp) :: expmx, rk, suma, t, xx, xk, xlim, xlow, xupper

       real(kind=dp), parameter :: QUART = 0.25_dp
       real(kind=dp), parameter :: HALF  = 0.50_dp
       real(kind=dp), parameter :: ONE   = 1.00_dp
       real(kind=dp), parameter :: FOUR  = 4.00_dp
       real(kind=dp), parameter :: EIGHT = 8.00_dp
       real(kind=dp), parameter :: NINE  = 9.00_dp
       real(kind=dp), parameter :: THIRT6= 36.0_dp
       real(kind=dp), parameter :: ONEHUN=100.0_dp
       real(kind=dp), parameter :: DEBINF= 0.60792710185402662866_dp

       real(kind=dp), parameter :: ADEB1(0:18) = (/  &
                                                 2.40065971903814101941_dp, 0.19372130421893600885_dp,  &
                                                -0.623291245548957703D-2,   0.35111747702064800D-3,     &
                                                -0.2282224667012310D-4,     0.158054678750300D-5,       &
                                                -0.11353781970719D-6,       0.835833611875D-8,          &
                                                -0.62644247872D-9,          0.4760334890D-10,           &
                                                -0.365741540D-11,           0.28354310D-12,             &
                                                -0.2214729D-13,             0.174092D-14,               &
                                                -0.13759D-15,               0.1093D-16,                 &
                                                -0.87D-18,                  0.7D-19,                    &
                                                -0.1D-19 /)

       !> Start computation
       xx = x

       !> Check xx >= 0.0
       if (xx < 0.0_dp) then
          !> Error activated
          ERR_MathGen=.true.
          ERR_MathGen_Mess="DEBYE1 doesn't work with negative Argument"
          fval = 0.0_dp
          return
       end if

       !> Constants.
       t = EPSILON(0.0_dp)
       xlow = SQRT(t*eight)
       xupper = -LOG(t+t)
       xlim = -LOG(TINY(0.0_dp))
       t = t / onehun
       do nterms = 18, 0, -1
          if (abs(adeb1(nterms)) > t) exit
       end do

       !> Code for xx <= 4.0
       if (xx <= four) then
          if (xx < xlow) then
             fval = ((xx-nine)*xx+thirt6) / thirt6
          else
             t = ((xx*xx/eight)-half) - half
             fval = cheval(nterms,adeb1,t) - quart * xx
          end if

       else
          !> Code for xx > 4.0
          fval = one / (xx*debinf)
          if (xx < xlim) then
             expmx = EXP(-xx)
             if (xx > xupper) then
                fval = fval - expmx * (one+one/xx)
             else
                suma = 0.0_dp
                rk = AINT(xlim/xx)
                nexp = INT(rk)
                xk = rk * xx
                do i = nexp, 1, -1
                   t = (one+one/xk) / rk
                   suma = suma * expmx + t
                   rk = rk - one
                   xk = xk - xx
                end do
                fval = fval - suma * expmx
             end if
          end if
       end if

       return
    End Function Debye1

    !!--++
    !!--++ Function Debye2(x)
    !!--++
    !!--++ Calculates the Debye function of order 2, defined as
    !!--++    DEBYE2(x) = 2*[Integral {0 to x} t*t/(exp(t)-1) dt] / (x*x)
    !!--++
    !!--++ The code uses Chebyshev series whose coefficients are given to 20
    !!--++ decimal places.
    !!--++
    !!--.. EXTRA INFORMATION
    !!--..
    !!--.. If X < 0.0 an error message is defined and the function returns the value 0.0
    !!--..
    !!--.. NTERMS: The no. of elements of the array ADEB2. The recommended value is such that
    !!--..         ABS(ADEB2(NTERMS)) < EPS/100 , with 1 <= NTERMS <= 18
    !!--..
    !!--..   XLOW: The value below which DEBYE2 = 1 - x/3 + x*x/24 to machine precision.
    !!--..         The recommended value is SQRT(8*EPSNEG)
    !!--..
    !!--.. XUPPER: The value above which DEBYE2 = (4*zeta(3)/x^2) - 2*exp(-x)(x^2+2x+1)/x^2.
    !!--..         The recommended value is -LOG(2*EPS)
    !!--..
    !!--..  XLIM1: The value above which DEBYE2 = 4*zeta(3)/x^2
    !!--..         The recommended value is -LOG(XMIN)
    !!--..
    !!--..  XLIM2: The value above which DEBYE2 = 0.0 to machine precision.
    !!--..         The recommended value is SQRT(4.8/XMIN)
    !!--++
    !!--++ Update: January 2017
    !!
    Function Debye2(X) Result(Fval)
       !---- Argument ----!
       real(kind=dp), intent(in) :: X ! Argument

       !---- Local Variables ----!
       integer       :: i, nexp, nterms
       real(kind=dp) :: fval
       real(kind=dp) :: expmx, rk, suma, t, xx, xk, xlim1, xlim2, xlow, xupper

       real(kind=dp), parameter :: QUART = 0.25_dp
       real(kind=dp), parameter :: HALF  = 0.50_dp
       real(kind=dp), parameter :: ONE   = 1.00_dp
       real(kind=dp), parameter :: TWO   = 2.00_dp
       real(kind=dp), parameter :: THREE = 3.00_dp
       real(kind=dp), parameter :: FOUR  = 4.00_dp
       real(kind=dp), parameter :: EIGHT = 8.00_dp
       real(kind=dp), parameter :: TWENT4= 24.0_dp
       real(kind=dp), parameter :: ONEHUN=100.0_dp
       real(kind=dp), parameter :: DEBINF= 4.80822761263837714160_dp

       real(kind=dp), parameter :: ADEB2(0:18) = (/  &
                                                 2.59438102325707702826_dp, 0.28633572045307198337_dp,  &
                                                -0.1020626561580467129D-1,  0.60491097753468435D-3,     &
                                                -0.4052576589502104D-4,     0.286338263288107D-5,       &
                                                -0.20863943030651D-6,       0.1552378758264D-7,         &
                                                -0.117312800866D-8,         0.8973585888D-10,           &
                                                -0.693176137D-11,           0.53980568D-12,             &
                                                -0.4232405D-13,             0.333778D-14,               &
                                                -0.26455D-15,               0.2106D-16,                 &
                                                -0.168D-17,                 0.13D-18,                   &
                                                -0.1D-19 /)

       !> Start computation
       xx = x

       !> Check xx >= 0.0
       if (xx < 0.0_dp) then
          ! Error activated
          ERR_MathGen=.true.
          ERR_MathGen_Mess="DEBYE2 doesn't work with negative Argument"
          fval = 0.0_dp
          return
       end if

       !> Constants.
       t = TINY(0.0_dp)
       xlim1 = -LOG(t)
       xlim2 = SQRT(debinf) / SQRT(t)
       t = EPSILON(0.0_dp)
       xlow = SQRT(t*eight)
       xupper = -LOG(t+t)
       t = t / onehun
       do nterms = 18, 0, -1
          if (ABS(adeb2(nterms)) > t) exit
       end do

       !> Code for xx <= 4.0
       if (xx <= four) then
          if (xx < xlow) then
             fval = ((xx-eight)*xx+twent4) / twent4
          else
             t = ((xx*xx/eight)-half) - half
             fval = cheval(nterms,adeb2,t) - xx / three
          end if

       else
          !> Code for xx > 4.0
          if (xx > xlim2) then
             fval = 0.0_dp
          else
             fval = debinf / (xx*xx)
             if (xx < xlim1) then
                expmx = EXP(-xx)
                if (xx > xupper) then
                   suma = ((xx+two)*xx+two) / (xx*xx)
                else
                   suma = 0.0_dp
                   rk = AINT(xlim1/xx)
                   nexp = INT(rk)
                   xk = rk * xx
                   do i = nexp, 1, -1
                      t = (one+two/xk+two/(xk*xk)) / rk
                      suma = suma * expmx + t
                      rk = rk - one
                      xk = xk - xx
                   end do
                end if
                fval = fval - two * suma * expmx
             end if
          end if
       end if

       return
    End Function Debye2

    !!--++
    !!--++ Function Debye3(x)
    !!--++
    !!--++ Calculates the Debye function of order 3, defined as
    !!--++    DEBYE3(x) = 3*[Integral {0 to x} t^3/(exp(t)-1) dt] / (x^3)
    !!--++
    !!--++ The code uses Chebyshev series whose coefficients are given to 20
    !!--++ decimal places.
    !!--++
    !!--.. EXTRA INFORMATION
    !!--..
    !!--.. If X < 0.0 an error message is defined and the function returns the value 0.0
    !!--..
    !!--.. NTERMS: The no. of elements of the array ADEB3. The recommended value is such that
    !!--..         ABS(ADEB3(NTERMS)) < EPS/100 , with 1 <= NTERMS <= 18
    !!--..
    !!--..   XLOW: The value below which DEBYE3 = 1 - 3x/8 + x*x/20 to machine precision.
    !!--..         The recommended value is SQRT(8*EPSNEG)
    !!--..
    !!--.. XUPPER: The value above which DEBYE3 = (18*zeta(4)/x^3) - 3*exp(-x)(x^3+3x^2+6x+6)/x^3.
    !!--..         The recommended value is -LOG(2*EPS)
    !!--..
    !!--..  XLIM1: The value above which DEBYE3 = 18*zeta(4)/x^3
    !!--..         The recommended value is -LOG(XMIN)
    !!--..
    !!--..  XLIM2: The value above which DEBYE3 = 0.0 to machine precision.
    !!--..         The recommended value is CUBE ROOT(19/XMIN)
    !!--++
    !!--++ Update: January 2017
    !!
    Function Debye3(X) Result(Fval)
       !---- Argument ----!
       real(kind=dp), intent(in) :: X ! Argument

       !---- Local Variables ----!
       integer       :: i, nexp, nterms
       real(kind=dp) :: fval
       real(kind=dp) :: expmx, rk, suma, t, xx, xk, xki, xlim1, xlim2, xlow, xupper

       real(kind=dp), parameter :: PT375 = 0.375_dp
       real(kind=dp), parameter :: HALF  = 0.50_dp
       real(kind=dp), parameter :: ONE   = 1.00_dp
       real(kind=dp), parameter :: THREE = 3.00_dp
       real(kind=dp), parameter :: FOUR  = 4.00_dp
       real(kind=dp), parameter :: SIX   = 6.00_dp
       real(kind=dp), parameter :: SEVP5 = 7.50_dp
       real(kind=dp), parameter :: EIGHT = 8.00_dp
       real(kind=dp), parameter :: TWENTY= 20.0_dp
       real(kind=dp), parameter :: ONEHUN=100.0_dp
       real(kind=dp), parameter :: DEBINF= 0.51329911273421675946D-1

       real(kind=dp), parameter :: ADEB3(0:18) = (/  &
                                                 2.70773706832744094526_dp, 0.34006813521109175100_dp,  &
                                                -0.1294515018444086863D-1,  0.79637553801738164D-3,     &
                                                -0.5463600095908238D-4,     0.392430195988049D-5,       &
                                                -0.28940328235386D-6,       0.2173176139625D-7,         &
                                                -0.165420999498D-8,         0.12727961892D-9,           &
                                                -0.987963459D-11,           0.77250740D-12,             &
                                                -0.6077972D-13,             0.480759D-14,               &
                                                -0.38204D-15,               0.3048D-16,                 &
                                                -0.244D-17,                 0.20D-18,                   &
                                                -0.2D-19 /)

       !> Start computation
       xx = x

       !> Error test
       if (xx < 0.0_dp) then
          ! Error activated
          ERR_MathGen=.true.
          ERR_MathGen_Mess="DEBYE3 doesn't work with negative Argument"
          fval = 0.0_dp
          return
       end if

       !> Constants.
       t = TINY(0.0_dp)
       xlim1 = -LOG(t)
       xk = one / three
       xki = (one/debinf) ** xk
       rk = t ** xk
       xlim2 = xki / rk
       t = EPSILON(0.0_dp)
       xlow = SQRT(t*eight)
       xupper = -LOG(t+t)
       t = t / onehun
       do nterms = 18, 0, -1
          if (ABS(adeb3(nterms)) > t) exit
       end do

       !> Code for xx <= 4.0
       if (xx <= four) then
          if (xx < xlow) then
             fval = ((xx-sevp5)*xx+twenty) / twenty
          else
             t = ((xx*xx/eight)-half) - half
             fval = cheval(nterms,adeb3,t) - pt375 * xx
          end if

       else
          !> Code for xx > 4.0
          if (xx > xlim2) then
             fval = 0.0_dp
          else
             fval = one / (debinf*xx*xx*xx)
             if (xx < xlim1) then
                expmx = EXP(-xx)
                if (xx > xupper) then
                   suma = (((xx+three)*xx+six)*xx+six) / (xx*xx*xx)
                else
                   suma = 0.0_dp
                   rk = AINT(xlim1/xx)
                   nexp = INT(rk)
                   xk = rk * xx
                   do i = nexp, 1, -1
                      xki = one / xk
                      t = (((six*xki+six)*xki+three)*xki+one) / rk
                      suma = suma * expmx + t
                      rk = rk - one
                      xk = xk - xx
                   end do
                end if
                fval = fval - three * suma * expmx
             end if
          end if
       end if

       return
    End Function Debye3

    !!--++
    !!--++ Function Debye4(x)
    !!--++
    !!--++ Calculates the Debye function of order 4, defined as
    !!--++    DEBYE4(x) = 4*[Integral {0 to x} t^4/(exp(t)-1) dt] / (x^4)
    !!--++
    !!--++ The code uses Chebyshev series whose coefficients are given to 20
    !!--++ decimal places.
    !!--++
    !!--.. EXTRA INFORMATION
    !!--..
    !!--.. If X < 0.0 an error message is defined and the function returns the value 0.0
    !!--..
    !!--.. NTERMS: The no. of elements of the array ADEB4. The recommended value is such that
    !!--..         ABS(ADEB4(NTERMS)) < EPS/100 , with 1 <= NTERMS <= 18
    !!--..
    !!--..   XLOW: The value below which  DEBYE4 = 1 - 4x/10 + x*x/18 to machine precision.
    !!--..         The recommended value is SQRT(8*EPSNEG)
    !!--..
    !!--.. XUPPER: The value above which DEBYE4=(96*zeta(5)/x^4)-4*exp(-x)(x^4+4x^2+12x^2+24x+24)/x^4
    !!--..         The recommended value is -LOG(2*EPS)
    !!--..
    !!--..  XLIM1: The value above which DEBYE4 = 96*zeta(5)/x^4
    !!--..         The recommended value is -LOG(XMIN)
    !!--..
    !!--..  XLIM2: The value above which DEBYE4 = 0.0 to machine precision.
    !!--..         The recommended value is FOURTH ROOT(99/XMIN)
    !!--++
    !!--++ Update: January 2017
    !!
    Function Debye4(X) Result(FVal)
       !---- Argument ----!
       real(kind=dp), intent(in) :: X ! Argument

       !---- Local Variables ----!
       integer       :: i, nexp, nterms
       real(kind=dp) :: fval
       real(kind=dp) :: expmx, rk, suma, t, xx, xk, xki, xlim1, xlim2, xlow, xupper

       real(kind=dp), parameter :: HALF  = 0.50_dp
       real(kind=dp), parameter :: ONE   = 1.00_dp
       real(kind=dp), parameter :: TWOPT5= 2.50_dp
       real(kind=dp), parameter :: FOUR  = 4.00_dp
       real(kind=dp), parameter :: FIVE  = 5.00_dp
       real(kind=dp), parameter :: EIGHT = 8.00_dp
       real(kind=dp), parameter :: TWELVE= 12.0_dp
       real(kind=dp), parameter :: EIGHTN= 18.0_dp
       real(kind=dp), parameter :: TWENT4= 24.0_dp
       real(kind=dp), parameter :: FORTY5= 45.0_dp
       real(kind=dp), parameter :: ONEHUN=100.0_dp
       real(kind=dp), parameter :: DEBINF= 99.54506449376351292781_dp


       real(kind=dp), parameter :: ADEB4(0:18) = (/  &
                                                 2.78186941502052346008_dp, 0.37497678352689286364_dp,  &
                                                -0.1494090739903158326D-1,  0.94567981143704274D-3,     &
                                                -0.6613291613893255D-4,     0.481563298214449D-5,       &
                                                -0.35880839587593D-6,       0.2716011874160D-7,         &
                                                -0.208070991223D-8,         0.16093838692D-9,           &
                                                -0.1254709791D-10,          0.98472647D-12,             &
                                                -0.7772369D-13,             0.616483D-14,               &
                                                -0.49107D-15,               0.3927D-16,                 &
                                                -0.315D-17,                 0.25D-18,                   &
                                                -0.2D-19 /)

       !> Start computation
       xx = x

       !> Error test
       if (xx < 0.0_dp) then
          ! Error activated
          ERR_MathGen=.true.
          ERR_MathGen_Mess="DEBYE4 doesn't work with negative Argument"
          fval = 0.0_dp
          return
       end if

       !> Constants.
       t = TINY(0.0_dp)
       xlim1 = -LOG(t)
       rk = one / four
       xk = debinf ** rk
       xki = t ** rk
       xlim2 = xk / xki
       t = EPSILON(0.0_dp)
       xlow = SQRT(t*eight)
       xupper = -LOG(t+t)
       t = t / onehun
       do nterms = 18, 0, -1
          if (ABS(adeb4(nterms)) > t) exit
       end do

       !> Code for xx <= 4.0
       if (xx <= four) then
          if (xx < xlow) then
             fval = ((twopt5*xx-eightn)*xx+forty5) / forty5
          else
             t = ((xx*xx/eight)-half) - half
             fval = cheval(nterms,adeb4,t) - (xx+xx) / five
          end if

       else
          !> Code for xx > 4.0
          if (xx > xlim2) then
             fval = 0.0_dp
          else
             t = xx * xx
             fval = (debinf/t) / t
             if (xx < xlim1) then
                expmx = EXP(-xx)
                if (xx > xupper) then
                   suma = ((((xx+four)*xx+twelve)*xx+twent4)*xx+twent4) / (xx*xx*xx*xx )
                else
                   suma = 0.0_dp
                   rk = INT(xlim1/xx)
                   nexp = INT(rk)
                   xk = rk * xx
                   do i = nexp, 1, -1
                      xki = one / xk
                      t = ((((twent4*xki+twent4)*xki+twelve)*xki+four)*xki+one ) / rk
                      suma = suma * expmx + t
                      rk = rk - one
                      xk = xk - xx
                   end do
                end if
                fval = fval - four * suma * expmx
             end if
          end if
       end if

       return
    End Function Debye4

    !!--++
    !!--++ Function DebyeN(n,x)
    !!--++
    !!--++ Calculates the Debye function of order N of X.
    !!--++
    !!--++ Limitation: |x| < 2*Pi
    !!--++
    !!--++ Update: January 2017
    !!
    Function DebyeN(n,x) Result(Fval)
       !---- Arguments ----!
       integer,       intent(in) :: N ! Order of Debye function
       real(kind=dp), intent(in) :: X ! Argument

       !---- Local Variables ----!
       integer, parameter                        :: KMAX=12
       integer                                   :: k,i
       real(kind=dp)                             :: Fval
       real(kind=dp)                             :: den,t1,t2
       real(kind=dp), dimension(KMAX), parameter :: B2K=(/                                       &
                                                      1.0_dp/  6.0_dp,        -1.0_dp/  30.0_dp, &
                                                      1.0_dp/ 42.0_dp,        -1.0_dp/  30.0_dp, &
                                                      5.0_dp/ 66.0_dp,      -691.0_dp/2730.0_dp, &
                                                      7.0_dp/  6.0_dp,     -3617.0_dp/ 510.0_dp, &
                                                  43867.0_dp/798.0_dp,   -174611.0_dp/ 330.0_dp, &
                                                 854513.0_dp/138.0_dp,-236364091.0_dp/2730.0_dp/)

       !> Init
       fval=0.0_dp

       !> Check
       if (abs(x) > 2.0*pi) then
          ! Error Flag
          ERR_MathGen=.true.
          ERR_MathGen_Mess="The absolute value of argument for DEBYEN was greater than 2PI "
          return
       end if

       !> Calculation
       t1=dble(n)/dble(2*(n+1))

       t2=0.0_dp
       do k=1,kmax
          i=2*k
          den=dble(2*k+n)*factorial_dp(i)
          t2 = t2 + (B2K(k)/den)*(x**(2*k))
       end do
       fval = 1.0_dp - (t1*x) + (n*t2)

       return
    End Function DebyeN

    !!--++
    !!--++ Function Determinant_C(Array,n,determ)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real square matrix.
    !!--++    Calculates the pseudo-determinant of a complex square matrix.
    !!--++    determ=det(AR)^2 + det(AI)^2 (if complex)
    !!--++
    !!--++    The calculated value is only useful for linear dependency purposes.
    !!--++    It tell us if the complex matrix is singular or not.
    !!--++
    !!--++    P R O V I S I O N A L (The determinant of A is not calculated at present)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determinant_C(Array,n) Result (Determ)
       !---- Arguments ----!
       complex, dimension(:,:), intent( in) :: Array   ! Complex imput array
       integer,                 intent( in) :: n       ! Dimension of Array
       real(kind=cp)                        :: determ

       !---- local variables ----!
       real(kind=cp),    dimension(2*n,2*n) :: AC   !real square matrix
       real(kind=cp)                        :: d
       integer                              :: i,nn
       logical                              :: singular

       !> Init
       Determ=0.0_cp

       nn=2*n
       AC(  1:n ,  1:n ) =  real(Array(1:n ,1:n))
       AC(n+1:nn,  1:n ) = aimag(Array(1:n ,1:n))
       AC(n+1:nn,n+1:nn) =    AC(  1:n ,1:n)
       AC(  1:n ,n+1:nn) =   -AC(n+1:nn,1:n)

       call lu_decomp(ac(1:nn,1:nn),d,singular)

       if (.not. singular) then
          do i=1,nn
             d=d*sign(1.0_cp,ac(i,i))
             determ=determ+ log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Function Determinant_C

    !!--++
    !!--++ Function Determinant_R(A,n,determ)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real square matrix.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determinant_R(Array,n) Result(Determ)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent( in) :: Array  ! Real matrix (N,N)
       integer,                       intent( in) :: n      ! Dimension of Array
       real(kind=cp)                              :: determ ! Value

       !---- local variables ----!
       real(kind=cp),    dimension(n,n)  :: AC
       real(kind=cp)                     :: d
       integer                           :: i
       logical                           :: singular

       !> Init
       determ=0.0_cp

       ac=Array(1:n,1:n)
       call lu_decomp(ac,d,singular)

       if (.not. singular) then
          do i=1,n
             d=d*sign(1.0_cp,ac(i,i))
             determ=determ + log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Function Determinant_R

    !!--++
    !!--++ FUNCTION ENVJ
    !!--++
    !!--++    (PRIVATE)
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Pure Function Envj(N,X) Result(Y)
       !---- Arguments ----!
       integer,       intent(in) :: n
       real(Kind=dp), intent(in) :: x
       real(Kind=dp)             :: y

       y=0.5_dp*log10(6.28_dp*real(n,kind=dp))-n*log10(1.36_dp*x/real(n,kind=dp))

       return
    End Function Envj

    !!--++
    !!--++ Logical Function Equal_Matrix_I(Arr1, Arr2, N)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Matrix_I(Array1,Array2,n) result(info)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: Array1 ! First Array(N,N)
       integer, dimension(:,:), intent(in) :: Array2 ! Second array(N,N)
       integer                , intent(in) :: n      ! Dimension of arrays
       logical                             :: info

       !---- Local variables ----!
       integer :: i,j

       info=.false.
       do i=1,n
          do j=1,n
             if (array1(i,j) /= array2(i,j)) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_I

    !!--++
    !!--++ Logical Function Equal_Matrix_R(Array1, Array2, N)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Matrix_R(Array1,Array2,n) result(info)
       !---- Argument ----!
       real(kind=cp), dimension(:,:)   , intent(in) :: Array1 ! First Array(N,N)
       real(kind=cp), dimension(:,:)   , intent(in) :: Array2 ! Second Array(N,N)
       integer,                          intent(in) :: n      ! Dimension of arrays
       logical                                      :: info

       !---- Local variables ----!
       integer :: i,j

       info=.false.
       do i=1,n
          do j=1,n
             if (abs(array1(i,j) - array2(i,j)) > epss ) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_R

    !!--++
    !!--++ Logical Function Equal_Vector_I(Vec1, Vec2, N)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Vector_I(Vec1,Vec2,n) result(info)
       !---- Argument ----!
       integer, dimension(:),   intent(in) :: Vec1   ! First vector
       integer, dimension(:),   intent(in) :: Vec2   ! Second vector
       integer                , intent(in) :: n      ! Dimension of the vectors
       logical                             :: info

       !---- Local variables ----!
       integer :: i

       info=.false.
       do i=1,n
          if (vec1(i) /= vec2(i)) return
       end do
       info=.true.

       return
    End Function Equal_Vector_I

    !!--++
    !!--++ Logical Function Equal_Vector_R(Vec1, Vec2, N)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real(kind=sp) vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Vector_R(Vec1,Vec2,n) result(info)
       !---- Argument ----!
       real(kind=cp), dimension(:)   ,   intent(in) :: Vec1   ! First vector
       real(kind=cp), dimension(:)   ,   intent(in) :: Vec2   ! Second vector
       integer,                          intent(in) :: n
       logical                                      :: info

       !---- Local variables ----!
       integer :: i

       info=.false.
       do i=1,n
          if (abs(vec1(i) - vec2(i)) > epss ) return
       end do
       info=.true.

       return
    End Function Equal_Vector_R

    !!----
    !!----  Function Euclidean_Norm(Vec,n) Result(Fn_Val)
    !!----
    !!----  This function calculates safely the Euclidean norm of a vector.
    !!----  Intermediate overflows are avoided using this function. The original
    !!----  name "enorm" from MINPACK has been changed and the subroutine has
    !!----  been translated to Fortran 90.
    !!----
    !!----
    !!--..  Original documentation (from MINPACK):
    !!--..
    !!--..  Function enorm
    !!--..
    !!--..  Given an n-vector x, this function calculates the euclidean norm of x.
    !!--..
    !!--..  The euclidean norm is computed by accumulating the sum of squares in
    !!--..  three different sums.  The sums of squares for the small and large
    !!--..  components are scaled so that no overflows occur.  Non-destructive
    !!--..  underflows are permitted.  Underflows and overflows do not occur in the
    !!--..  computation of the unscaled sum of squares for the intermediate
    !!--..  components.  The definitions of small, intermediate and large components
    !!--..  depend on two constants, rdwarf and rgiant.  The main restrictions on
    !!--..  these constants are that rdwarf**2 not underflow and rgiant**2 not
    !!--..  overflow.  The constants given here are suitable for every known computer.
    !!--..
    !!--..  Argonne National Laboratory. MINPACK project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!----
    !!----  Update: August - 2009
    !!----
    Function Euclidean_Norm(Vec,n) Result(Fn_Val)
       !---- Arguments ----!
       real (kind=cp), dimension(:), intent(In)  :: Vec  ! Input vector
       integer,                      intent(In)  :: n    ! Dimension of the Vector
       real (kind=cp)                            :: Fn_Val

       !--- Local Variables ---!
       integer                   :: i
       real (Kind=cp)            :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max

       real (Kind=cp), Parameter :: ONE = 1.0_cp
       real (Kind=cp), Parameter :: ZERO= 0.0_cp
       real (Kind=cp), Parameter :: RDWARF = 3.834e-20_cp
       real (Kind=cp), Parameter :: RGIANT = 1.304e+19_cp

       !> Init
       s1 = zero
       s2 = zero
       s3 = zero
       x1max = zero
       x3max = zero
       floatn = real(n)
       agiant = rgiant/floatn

       do i = 1, n
          xabs = Abs(vec(i))
          if (.not. (xabs > rdwarf .AND. xabs < agiant)) then
             !> sum for large components.
             if (xabs > rdwarf) then
                if (xabs > x1max) then
                   s1 = one + s1*(x1max/xabs)**2
                   x1max = xabs
                   cycle
                end if
                s1 = s1 + (xabs/x1max)**2
                cycle
             end if

             !> sum for small components.
             If (xabs > x3max) Then
                s3 = one + s3*(x3max/xabs)**2
                x3max = xabs
                cycle
             end if

             if (xabs /= zero) s3 = s3 + (xabs/x3max)**2
             cycle
          end if

          !>  sum for intermediate components.
          s2 = s2 + xabs**2
       end do

       !> calculation of Norm.
       if (s1 /= zero) Then
          Fn_Val = x1max*Sqrt(s1 + (s2/x1max)/x1max)
          return
       end if

       if (s2 /= zero) Then
          if (s2 >= x3max) Fn_Val = Sqrt(s2*(one + (x3max/s2)*(x3max*s3)))
          if (s2 < x3max) Fn_Val = Sqrt(x3max*((s2/x3max) + (x3max*s3)))
          return
       end if

       Fn_Val = x3max*Sqrt(s3)

       Return
    End Function Euclidean_Norm

    !!--++
    !!--++ Function in_limits_dp(vec,n,limits) result(ok)
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: March - 2013
    !!
    Function IN_Limits_dp(vec,n,limits) result(ok)
       !---- Arguments ----!
       real(kind=dp), dimension(:),   intent(in) :: vec     ! Vector
       integer,                       intent(in) :: n       ! Dimension of the input vector
       real(kind=dp), dimension(:,:), intent(in) :: limits  ! Normally (2,n)
       logical                                   :: ok

       !---- Local Variables ----!
       integer :: i

       !> Init
       ok=.false.
       if (size(vec) /= n) return

       ok=.true.
       do i=1,n
          if(vec(i) >= limits(1,i) .and. vec(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_dp

    !!--++
    !!--++ Function in_limits_int(vec,n,limits) result(ok)
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: March - 2013
    !!
    Function IN_Limits_int(vec,n,limits) result(ok)
       !---- Arguments ----!
       integer, dimension(:),   intent(in) :: vec      ! Vector
       integer,                 intent(in) :: n        ! Dimension of Vector
       integer, dimension(:,:), intent(in) :: limits   ! Normally (2,n)
       logical                             :: ok

       !---- Local Variables ----!
       integer :: i

       !> Init
       ok=.false.
       if (size(vec) /= n) return

       ok=.true.
       do i=1,n
          if (vec(i) >= limits(1,i) .and. vec(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_int

    !!--++
    !!--++ Function in_limits_sp(vec,n,limits) result(ok)
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: March - 2013
    !!
    Function IN_Limits_sp(Vec,n,limits) result(ok)
       !---- Arguments ----!
       real(kind=sp), dimension(:),   intent(in) :: vec      ! Vector
       integer,                       intent(in) :: n        ! Dimension odf the vector
       real(kind=sp), dimension(:,:), intent(in) :: limits   ! Normally (2,n)
       logical                                   :: ok

       !---- Local Variables ----!
       integer :: i

       !> Init
       ok=.false.
       if (size(vec) /= n) return

       ok=.true.
       do i=1,n
          if (vec(i) >= limits(1,i) .and. vec(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_sp

    !!--++
    !!--++ Subroutine Linear_DependentC(a,na,b,nb,mb,info)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++
    !!--++    For the case of complex vectors in Cn the problem can be reduced to real vectors
    !!--++    of dimension R2n. Each complex vector contributes as two real vectors of dimension
    !!--++    2n: (R,I) and (-I,R). A complex vector V is linearly dependent on n complex vectors
    !!--++    if V can be written as: V = Sigma{j=1,n}(Cj.Vj), with Cj complex numbers and Vj
    !!--++    having n complex components. One may write:
    !!--++
    !!--++     V = Sigma{j=1,n}(Cj.Vj)
    !!--++     (R,I) = Sigma{j=1,n} (Cjr Vj + i Cji Vj) = Sigma{j=1,n} (Cjr (Rj,Ij) +  Cji (-Ij,Rj) )
    !!--++     (R,I) = Sigma{j=1,n} (aj (Rj,Ij) + bj (-Ij,Rj) )  = Sigma{j=1,2n} (Aj.Uj)
    !!--++     Were Uj=(Rj,Ij) and U(j+1)= (-Ij,Rj)
    !!--++
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Linear_DependentC(Vec,na,Array,nb,mb) Result(info)
       !---- Arguments ----!
       complex, dimension(:),   intent(in)  :: Vec        ! Input Vector
       integer,                 intent(in)  :: na         ! Dimension of Vec
       complex, dimension(:,:), intent(in)  :: Array      ! Array(NB,MB)
       integer,                 intent(in)  :: nb,mb      ! Dimension of Array
       logical                              :: info

       !---- Local variables ----!
       integer                                                     :: r,n1
       real(kind=dp), parameter                                    :: tol= 100.0_dp*deps
       real(kind=dp), dimension(2*max(nb+1,mb+1),2*max(nb+1,mb+1)) :: c

       !> Init
       call init_err_mathgen()
       info=.true.
       if (nb > size(array,1) .or. mb > size(array,2) .or. na > size(vec) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentC: Error in dimension of input matrix or vector"
          return
       end if

       c=0.0
       if ( na == mb) then
          n1=2*nb+1
          if(n1+1 > 2*mb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(array(1:nb,1:mb))
          c(1:nb,     mb+1:mb+na) = aimag(array(1:nb,1:mb))
          c(nb+1:2*nb,      1:mb) =-aimag(array(1:nb,1:mb))
          c(nb+1:2*nb,mb+1:mb+na) =  real(array(1:nb,1:mb))
          c(n1,             1:mb) =  real(vec(1:na))
          c(n1,      mb+1:mb+na ) = aimag(vec(1:na))
          c(n1+1,           1:mb) =-aimag(vec(1:na))
          c(n1+1,    mb+1:mb+na ) =  real(vec(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*mb)) info=.false.

       else if( na == nb) then
          n1=2*mb+1
          if(n1+1 > 2*nb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(array(1:nb,1:mb))
          c(nb+1:nb+na,     1:mb) = aimag(array(1:nb,1:mb))
          c(1:nb,      mb+1:2*mb) =-aimag(array(1:nb,1:mb))
          c(nb+1:nb+na,mb+1:2*mb) =  real(array(1:nb,1:mb))
          c(1:na,             n1) =  real(vec(1:na))
          c(nb+1:nb+na,       n1) = aimag(vec(1:na))
          c(1:na,           1+n1) =-aimag(vec(1:na))
          c(nb+1:nb+na,     1+n1) =  real(vec(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentC: input dimension of vector incompatible with matrix"
       end if

       return
    End Function Linear_DependentC

    !!--++
    !!--++ Function Linear_DependentI(Vec,na,Array,nb,mb) Result(info)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Linear_DependentI(Vec,na,Array,nb,mb) Result(info)
       !---- Arguments ----!
       integer, dimension(:),   intent(in)  :: Vec   ! Integer vector
       integer,                 intent(in)  :: na    ! Dimension of Vec
       integer, dimension(:,:), intent(in)  :: Array ! Array(NB,MB)
       integer,                 intent(in)  :: nb,mb ! Dimension of Array
       logical                              :: info

       !---- Local variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       !> Init
       call init_err_mathgen()
       info=.true.
       if (nb > size(array,1) .or. mb > size(array,2) .or. na > size(vec) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentI: Error in dimension of input matrix or vector"
          return
       end if

       c=0.0
       if ( na == mb) then
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(array(1:nb,1:mb))
          c(n1,  1:mb)=real(vec(1:na))      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) info=.false.

       else if( na == nb) then
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(array(1:nb,1:mb))
          c(1:nb,  n1)=real(vec(1:na))     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentI: input dimension of vector incompatible with matrix"
       end if

       return
    End Function Linear_DependentI

    !!--++
    !!--++ Function Linear_DependentR(Vec,na,Array,nb,mb) Result(info)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Linear_DependentR(Vec,na,Array,nb,mb) Result(info)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in)  :: Vec      ! Real vector
       integer,                       intent(in)  :: na       ! Dimension of Vec
       real(kind=cp), dimension(:,:), intent(in)  :: Array    ! Array(NB,MB)
       integer,                       intent(in)  :: nb,mb    ! Dimension of Array
       logical                                    :: info

       !---- Local Variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       !> Init
       call init_err_mathgen()
       info=.true.
       if (nb > size(array,1) .or. mb > size(array,2) .or. na > size(vec) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentR: Error in dimension of input matrix or vector"
          return
       end if

       c=0.0
       if ( na == mb) then    !Vector added as an additional row
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=array(1:nb,1:mb)
          c(n1,  1:mb)=vec(1:na)      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) info=.false.

       else if( na == nb) then   !Vector added as an additional column
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=array(1:nb,1:mb)
          c(1:nb,  n1)=vec(1:na)     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentR: input dimension of vector incompatible with matrix"
       end if

       return
    End Function Linear_DependentR

    !!--++
    !!--++ Function Locate_I(Vec, n, x) Result(j)
    !!--++
    !!--++    Procedure for locating the index J of a ordered vector Vec(N)
    !!--++    satisfying:
    !!--++
    !!--++               Vec(J) <= X < Vec(J+1)
    !!--++
    !!--++ Update: June - 2011
    !!
    Function Locate_I(Vec,n,x) Result(j)
       !---- Argument ----!
       integer, dimension(:), intent(in):: Vec  ! Vector
       integer ,              intent(in):: n    ! Dimension of Vector
       integer,               intent(in):: x    ! Value
       integer                          :: j    ! Index

       !---- Local Variables ----!
       integer :: jl, ju, jm

       !> Init
       j=0

       !> Check of limits
       if (x <= Vec(1)) then
          j=1
          return
       end if

       if (x >= Vec(n)) then
          j=n
          return
       end if

       !> Searching...
       jl=0
       ju=n+1
       do
          if (ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((Vec(n) > Vec(1)) .eqv. (x > Vec(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_I

    !!--++
    !!--++ Function Locate_Ib(Vec, x) Result(j)
    !!--++
    !!--++    Subroutine for locating the index J of a ordered Vector Vec(:)
    !!--++    satisfying:
    !!--++
    !!--++               Vec(J) <= X < Vec(J+1)
    !!--++
    !!--++ Update: June - 2011
    !!
    Function Locate_Ib(Vec,x) Result(j)
       !---- Argument ----!
       integer, dimension(:), intent(in):: Vec   ! Vector
       integer,               intent(in):: x     ! Value
       integer                          :: j     ! Index

       !---- Local Variables ----!
       integer :: jl, ju, jm, i1,i2


       !> Init
       j=0

       i1=lbound(Vec,dim=1)
       i2=ubound(Vec,dim=1)

       !> Checking
       if (x <= Vec(i1)) then
          j=i1
          return
       end if

       if (x >= Vec(i2)) then
          j=i2
          return
       end if

       !> Searching...
       jl=i1-1
       ju=i2+1
       do
          if (ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((Vec(i2) > Vec(i1)) .eqv. (x > Vec(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_Ib

    !!--++
    !!--++ Function Locate_R(Vec, n, x) Result(j)
    !!--++
    !!--++    Function for locating the index J of an ordered vector Vec(N)
    !!--++    satisfying:
    !!--++
    !!--++               Vec(J) <= X < Vec(J+1)
    !!--++
    !!--++ Update: June - 2011
    !!
    Function Locate_R(Vec,n,x) Result(j)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in):: Vec  ! Vector
       integer ,                    intent(in):: n    ! Dimension of Vector
       real(kind=cp),               intent(in):: x    ! Value
       integer                                :: j    ! Index

       !---- Local Variables ----!
       integer :: jl, ju, jm

       !> Init
       j=0

       !> Check
       if (x <= Vec(1)) then
          j=1
          return
       end if

       if (x >= Vec(n)) then
          j=n
          return
       end if

       !> Searching...
       jl=0
       ju=n+1
       do
          if (ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((Vec(n) > Vec(1)) .eqv. (x > Vec(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_R

    !!--++
    !!--++ Function Locate_Rb(Vec, x) Result(j)
    !!--++
    !!--++    Function for locating the index J of an ordered Vector Vec(:)
    !!--++    satisfying:
    !!--++
    !!--++               Vec(J) <= X < Vec(J+1)
    !!--++
    !!--++ Update: June - 2011
    !!
    Function Locate_Rb(Vec,x) Result(j)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in):: Vec  ! Vector
       real(kind=cp),               intent(in):: x    ! Value
       integer                                :: j    ! Index

       !---- Local Variables ----!
       integer :: jl, ju, jm, i1,i2

       i1=lbound(Vec,dim=1)
       i2=ubound(Vec,dim=1)

       !> Init
       j=0

       if (x <= Vec(i1)) then
          j=i1
          return
       end if

       if (x >= Vec(i2)) then
          j=i2
          return
       end if

       !> Searching...
       jl=i1-1
       ju=i2+1
       do
          if (ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((Vec(i2) > Vec(i1)) .eqv. (x > Vec(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_Rb

    !!----
    !!---- Function Lower_Triangular_I(Array,n) Result (T)
    !!----
    !!---- Calculate the Lower Triangular
    !!----
    !!---- Updated: October - 2014
    !!
    Function Lower_Triangular_I(Array,n) Result (T)
       !---- Arguments ----!
       integer, dimension(:,:), intent(in) :: Array  ! Array(N,N)
       integer,                 intent(in) :: n      ! Dimension
       integer, dimension(n,n)             :: T

       !---- Local Variables ----!
       integer :: i,j,p,q,m

       m=n
       p=size(Array(:,1))
       q=size(Array(1,:))

       if(n > p .or. n > q) m=min(p,q)

       !> Init
       T=0

       do j=1,m
          do i=j,m
             T(i,j)=Array(i,j)
          end do
       end do

       return
    End Function  Lower_Triangular_I

    !!----
    !!---- Function Lower_Triangular_R(Array,n) Result (T)
    !!----
    !!----   Updated: October - 2014
    !!----
    Function Lower_Triangular_R(Array,n) Result (T)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: Array ! Array (N,N)
       integer,                       intent(in) :: n     ! Dimension of Array
       real(kind=cp), dimension(n,n)             :: T

       !---- Local Variables ----!
       integer :: i,j,p,q,m


       m=n
       p=size(Array(:,1))
       q=size(Array(1,:))

       if(n > p .or. n > q) m=min(p,q)

       T=0
       do j=1,m
          do i=j,m
             T(i,j)=Array(i,j)
          end do
       end do

       return
    End Function Lower_Triangular_R

    !!----
    !!---- Function Matinv(Array,n) Result(a)
    !!----
    !!----  Inverting a real square matrix.
    !!----
    !!---- Update: February - 2005
    !!
    Function Matinv(Array,N) Result(A)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in) :: array
       integer,                       intent(in) :: n
       real(kind=cp), dimension(n,n)             :: a

       !---- Local variables ----!
       real(kind=cp)                 :: amax,savec
       integer, dimension(size(a,1)) :: ik,jk
       integer                       :: i,j,k,l

       !> Init
       a=0.0_cp
       a=array(1:n,1:n)

       !---- Subroutine to invert a real matrix ----!
       do k=1,n
          amax=0.0
          do
             do
                do i=k,n
                   do j=k,n
                      if (abs(amax)-abs(a(i,j)) > 0.0) cycle
                      amax=a(i,j)
                      ik(k)=i
                      jk(k)=j
                   end do
                end do
                i=ik(k)
                if (i-k < 0) cycle
                exit
             end do

             if (i-k /= 0) then
                do j=1,n
                   savec=a(k,j)
                   a(k,j)=a(i,j)
                   a(i,j)=-savec
                end do
             end if

             j=jk(k)
             if (j-k < 0) cycle
             exit
          end do

          if (j-k /= 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=a(i,j)
                a(i,j)=-savec
             end do
          end if

          do i=1,n
             if (i-k /= 0)  then
                a(i,k)=-a(i,k)/amax
             end if
          end do
          do i=1,n
             do j=1,n
                if (i-k == 0 .or. j-k == 0) cycle
                a(i,j)=a(i,j)+a(i,k)*a(k,j)
             end do
          end do
          do j=1,n
             if (j-k == 0)   cycle
             a(k,j)=a(k,j)/amax
          end do
          a(k,k)=1.0/amax
       end do     !k

       do l=1,n
          k=n-l+1
          j=ik(k)
          if (j-k > 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=-a(i,j)
                a(i,j)=savec
             end do
          end if
          i=jk(k)
          if (i-k > 0) then
             do j=1,n
                savec=a(k,j)
                a(k,j)=-a(i,j)
                a(i,j)=savec
             end do
          end if
       end do

       return
    End Function Matinv

    !!----
    !!---- Function Modulo_Lat(U)
    !!----
    !!----    Reduces a real vector to another with components in
    !!----    the interval [0,1)
    !!----
    !!---- Updated: February - 2005
    !!
    Function Modulo_Lat(Vec) result(v)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in) :: Vec
       real(kind=cp), dimension(1:size(Vec))   :: v

       v=mod(vec+10.0_cp,1.0_cp)

       return
    End Function  Modulo_Lat

    !!--++
    !!--++ Function Norm_I(Vec,G) Result(R)
    !!--++
    !!--++    Calculate the Norm of a vector
    !!--++
    !!--++ Update: April - 2009
    !!
    Function Norm_I(Vec,G) Result(R)
       !---- Arguments ----!
       integer,       dimension(:),   intent(in) :: Vec  ! Vector
       real(kind=cp), dimension(:,:), intent(in) :: G    ! Metric array
       real(kind=cp)                             :: r    ! Norm

       if (size(vec)*size(vec) /= size(g)) then
          r=tiny(0.0)
       else
          r=sqrt(dot_product(real(vec), matmul(g,real(vec))))
       end if

       return
    End Function Norm_I

    !!--++
    !!--++ Function Norm_R(X,G) Result(R)
    !!--++
    !!--++    Calculate the Norm of a vector
    !!--++
    !!--++ Update: April - 2009
    !!
    Function Norm_R(Vec,G) Result(R)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: Vec ! Vector
       real(kind=cp), dimension(:,:), intent(in) :: G   ! Metric array
       real(kind=cp)                             :: r   ! Norm

       if (size(vec)*size(vec) /= size(g)) then
          r=tiny(0.0)
       else
          r=sqrt(dot_product(vec, matmul(g,vec)))
       end if

       return
    End Function Norm_R

    !!--++
    !!--++ Function Outerprod_dp(Vec1,Vec2) Result(c)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = Vec1(i)*Vec2(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Outerprod_dp(Vec1,Vec2)  Result(c)
       !---- Arguments ----!
       real(kind=dp),dimension(:),intent(in)          :: Vec1   ! Vector 1
       real(kind=dp),dimension(:),intent(in)          :: Vec2   ! Vector 2
       real(kind=dp),dimension(size(Vec1),size(Vec2)) :: c      ! Tensor Matrix

       c =spread(vec1,dim=2,ncopies=size(vec2))*spread(vec2,dim=1,ncopies=size(vec1))

       return
    End Function Outerprod_dp

    !!--++
    !!--++ Function Outerprod_sp(Vec1,Vec2) Result(c)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = Vec1(i)*Vec2(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Outerprod_sp(Vec1,Vec2)  Result(c)
       !---- Arguments ----!
       real(kind=sp),dimension(:),intent(in)          :: Vec1  ! Vector 1
       real(kind=sp),dimension(:),intent(in)          :: Vec2  ! Vector 2
       real(kind=sp),dimension(size(Vec1),size(Vec2)) :: c

       c =spread(Vec1,dim=2,ncopies=size(Vec2))*spread(Vec2,dim=1,ncopies=size(Vec1))

       return
    End Function Outerprod_sp

    !!----
    !!---- Function Pgcd(i,j) Result(mcd)
    !!----
    !!----    Function calculating the maximum common divisor of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Function Pgcd(i,j) Result(mcd)
       !---- Arguments ----!
       integer, intent(in) :: i    ! Integer 1
       integer, intent(in) :: j    ! Integer2
       integer             :: mcd

       !---- Local variables ----!
       integer  :: u,v,m

       u=max(i,j)
       v=min(i,j)
       m=0
       do
          if (m == 1) exit
          m=mod(u,v)
          u=v
          v=m
       end do
       mcd=u

       return
    End Function Pgcd

    !!----
    !!---- Function Ppcm(i,j) result(mcm)
    !!----
    !!----    Function calculating the minimum common multiple of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Function Ppcm(i,j) result(mcm)
       !---- Arguments ----!
       integer, intent(in) :: i    ! Integer 1
       integer, intent(in) :: j    ! Integer 2
       integer             :: mcm

       !---- Local variables ----!
       integer :: u,v,w,ii

       u=max(i,j)
       v=min(i,j)
       mcm=1

       if (v <= 1) then
          mcm=u
          return
       end if

       w=int(sqrt(real(u)))+1
       do ii=2,w
          do
             if (.not. ((mod(u,ii)==0) .or. (mod(v,ii)==0)) ) exit
             mcm=mcm*ii
             if (modulo(u,ii) == 0) u=u/ii
             if (modulo(v,ii) == 0) v=v/ii
          end do
       end do

       return
    End Function Ppcm

    !!--++
    !!--++ Function Pythag_dp(a,b) Result (c)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes c=sqrt(a^2 +b^2 ) without destructive underflow or overflow.
    !!--++    Adapted from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Pythag_dp(a,b) Result (c)
       !---- Arguments ----!
       real(kind=dp),intent(in):: a    ! Real 1
       real(kind=dp),intent(in):: b    ! Real 2
       real(kind=dp)           :: c

       !---- Local variables ----!
       real(kind=dp)           :: absa,absb

       absa=abs(a)
       absb=abs(b)
       if (absa >absb)then
          c=absa*sqrt(1.0_dp+(absb/absa)**2)
       else
          if (absb < tiny(1.0_dp))then
             c=0.0
          else
             c=absb*sqrt(1.0_dp+(absa/absb)**2)
          end if
       end if

       return
    End Function Pythag_dp

    !!--++
    !!--++ Function Pythag_sp(a,b) result (c)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes c=sqrt(a^2 +b^2 ) without destructive underflow or overflow.
    !!--++    Adapted from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Pythag_sp(a,b) Result (c)
       !---- Arguments ----!
       real(kind=sp),intent(in):: a  ! Real 1
       real(kind=sp),intent(in):: b  ! Real 2
       real(kind=sp)           :: c

       !---- Local variables ----!
       real(kind=sp)           :: absa,absb

       absa=abs(a)
       absb=abs(b)
       if (absa > absb) then
          c=absa*sqrt(1.0_sp+(absb/absa)**2)
       else
          if (absb < tiny(1.0_sp)) then
             c=0.0
          else
             c=absb*sqrt(1.0_sp+(absa/absb)**2)
          end if
       end if

       return
    End Function Pythag_sp

    !!--++
    !!--++ Function Scalar_R(Vec1,Vec2,G) Result(R)
    !!--++
    !!--++    Scalar Product including metrics
    !!--++
    !!--++ Update: April - 2009
    !!
    Function Scalar_I(Vec1,Vec2,G) Result(R)
       !---- Arguments ----!
       integer, dimension(:),         intent(in) :: Vec1  ! Vector 1
       integer, dimension(:),         intent(in) :: Vec2  ! Vector 2
       real(kind=cp), dimension(:,:), intent(in) :: G     ! Metric array
       real(kind=cp)                             :: r

       if (size(Vec1)/= size(Vec2) .or. size(Vec1)*size(Vec1) /= size(g)) then
          r=tiny(0.0)
       else
          r=dot_product(real(Vec1), matmul(g,real(Vec2)))
       end if

       return
    End Function Scalar_I

    !!--++
    !!--++ Function Scalar_R(Vec1,Vec2,G) Result(R)
    !!--++
    !!--++    Scalar Product including metrics
    !!--++
    !!--++ Update: April - 2009
    !!
    Function Scalar_R(Vec1,Vec2,G) Result(R)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: Vec1  ! Vector 1
       real(kind=cp), dimension(:),   intent(in) :: Vec2  ! Vector 2
       real(kind=cp), dimension(:,:), intent(in) :: g     ! Metric array
       real(kind=cp)                             :: r

       if (size(Vec1)/= size(Vec2) .or. size(Vec1)*size(Vec1) /= size(g)) then
          r=tiny(0.0)
       else
          r=dot_product(Vec1, matmul(g,Vec2))
       end if

       return
    End Function Scalar_R

    !!----
    !!---- Function Splint(Xa, Ya, Y2a, n, x, y)
    !!----
    !!----    Cubic - spline interpolation at x value
    !!----
    !!---- Note: Y2a is obtained from spline procedure
    !!----
    !!---- Update: February - 2005
    !!
    Function Splint(Xa, Ya, Y2a, N, X) Result(y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: Xa   ! Vector of X points
       real(kind=cp), dimension(:), intent(in)  :: Ya   ! Vecot of Y=F(Xi) points
       real(kind=cp), dimension(:), intent(in)  :: Y2a  ! Second derivatives at Xi points (Spline)
       integer ,                    intent(in)  :: N    ! Number of points
       real(kind=cp),               intent(in)  :: X    ! X Value to be evaluated
       real(kind=cp)                            :: Y    ! The value

       !---- Local Variables ----!
       integer          :: klo, khi, k
       real(kind=cp)    :: h, a, b

       klo=1
       khi=n
       do
          if (khi-klo > 1) then
             k=(khi+klo)/2
             if (xa(k) > x) then
                khi=k
             else
                klo=k
             end if
             cycle
          else
             exit
          end if
       end do

       h=xa(khi)-xa(klo)
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h

       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)* y2a(khi))*(h**2)/6.0

       return
    End Function Splint

    !!--++
    !!--++ FUNCTION START1
    !!--++
    !!--++    (PRIVATE)
    !!--++    Determine the starting point for backward
    !!--++    recurrence such that the magnitude of Jn(x) at that point is
    !!--++    about 10^(-MP).
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Pure Function Start1(X,Mp) Result (Start)
       !---- Arguments ----!
       real(kind=dp), intent(in) :: x         ! Argument of Jn(x)
       integer, intent(in)       :: mp        ! Value of magnitude
       integer                   :: start

       !---- Local variables ----!
       integer      :: n1,n0,nn, it
       real(kind=dp):: f,f0,f1,a0

       a0=abs(x)
       n0=int(1.1_dp*a0)+1
       f0=envj(n0,a0)-mp
       n1=n0+5
       f1=envj(n1,a0)-mp
       do it=1,20
          nn=n1-(n1-n0)/(1.0_dp-f0/f1)
          f=envj(nn,a0)-mp
          if (abs(nn-n1) < 1) exit
          n0=n1
          f0=f1
          n1=nn
          f1=f
       end do
       start=nn

       return
    End Function Start1

    !!--++
    !!--++ FUNCTION START2
    !!--++
    !!--++    (PRIVATE)
    !!--++    Determine the starting point for backward
    !!--++    recurrence such that all Jn(x) has MP significants digits
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Pure Function Start2(X,N,Mp) Result(Start)
       !---- Arguments ----!
       real(kind=dp), intent(in) :: x     ! Argument of Jn(x)
       integer,       intent(in) :: n     ! Order of Jn(x)
       integer,       intent(in) :: mp    ! Significant digit
       integer                   :: start

       !---- Local variables ----!
       real(kind=dp) :: a0, hmp, ejn, obj,f,f0,f1
       integer       :: n0,n1,nn, it

       a0=abs(x)
       hmp=0.5_dp*mp
       ejn=envj(n,a0)
       if (ejn <= hmp) then
          obj=mp
          n0=int(1.1_dp*a0)+1  ! +1 was absent in the original version ... this solves the problem of
       else                    ! Intel, gfortran and g95 compilers ... Lahey was calculating well event if n0=0!
          obj=hmp+ejn
          n0=n
       end if
       f0=envj(n0,a0)-obj
       n1=n0+5
       f1=envj(n1,a0)-obj
       do it=1,20
          nn=n1-(n1-n0)/(1.0_dp-f0/f1)
          f=envj(nn,a0)-obj
          if (abs(nn-n1) < 1) exit
          n0=n1
          f0=f1
          n1=nn
          f1=f
       end do
       start=nn+10

       return
    End Function Start2

    !!--++
    !!--++ Function Trace_C(Array)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a complex nxn array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Trace_C(Array) Result(b)
       !---- Argument ----!
       complex, dimension(:,:), intent(in) :: array   ! Complex array
       complex                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=(0.0,0.0)
       imax=min(size(array,1),size(array,2))
       do i=1,imax
          b=b+array(i,i)
       end do

       return
    End Function Trace_C

    !!--++
    !!--++ Function Trace_I(Array)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of an integer 3x3 array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Trace_I(Array) Result(b)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: Array  ! Integer array
       integer                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0
       imax=min(size(array,1),size(array,2))
       do i=1,imax
          b=b+array(i,i)
       end do

       return
    End Function Trace_I

    !!--++
    !!--++ Function Trace_R(Array)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a real 3x3 array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Trace_R(Array) Result(b)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: Array  ! Real array
       real(kind=cp)                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0.0
       imax=min(size(array,1),size(array,2))
       do i=1,imax
          b=b+array(i,i)
       end do

       return
    End Function Trace_R

    !!----
    !!---- Function Upper_Triangular_I(Array,n) Result (T)
    !!----
    !!---- Calculate the Upper Triangular
    !!----
    !!----   Updated: October - 2014
    !!
    Function Upper_Triangular_I(Array,n) Result (T)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: Array   ! Integer array
       integer,                 intent(in) :: n       ! Dimension of Array
       integer, dimension(n,n)             :: T

       !---- Local Variables ----!
       integer :: i,j,p,q,m

       m=n
       p=size(Array(:,1))
       q=size(Array(1,:))
       if (n > p .or. n > q) m=min(p,q)

       T=0
       do j=1,m
          do i=1,j
             T(i,j)=Array(i,j)
          end do
       end do

       return
    End Function Upper_Triangular_I

    !!----
    !!---- Function Upper_Triangular_R(Array,n) Result (T)
    !!----
    !!---- Calculate the Upper Triangular
    !!----
    !!---- Updated: October - 2014
    !!
    Function Upper_Triangular_R(Array,n) Result (T)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: Array ! Real array
       integer,                       intent(in) :: n     ! Dimension of array
       real(kind=cp), dimension(n,n)             :: T

       !---- Local Variables ----!
       integer :: i,j,p,q,m


       m=n
       p=size(Array(:,1))
       q=size(Array(1,:))

       if (n > p .or. n > q) m=min(p,q)

       T=0
       do j=1,m
          do i=1,j
             T(i,j)=Array(i,j)
          end do
       end do

       return
    End Function Upper_Triangular_R

    !!--++
    !!--++ Logical Function ZbelongM(Array)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real array is an Integer matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function ZbelongM(Array) Result(belong)
       !---- Argument ----!
       real(kind=cp),   dimension(:,:), intent( in) :: Array  ! Real array
       logical                                      :: belong

       !---- Local variables ----!
       real(kind=cp),   dimension(size(array,1),size(array,2)) :: vec

       vec= abs(real(nint (array))-array)
       belong=.not. ANY(vec > epss)

       return
    End Function ZbelongM

    !!--++
    !!--++ Logical Function ZbelongN(A)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is an Integer
    !!--++
    !!--++ Update: February - 2005
   !!
    Function ZbelongN(a) Result(belong)
       !---- Argument ----!
       real(kind=cp), intent( in) :: a       ! Real number
       logical                    :: belong

       belong=.false.
       if (abs(real(nint (a))-a) > epss) return
       belong=.true.

       return
    End Function ZbelongN

    !!--++
    !!--++ Logical Function ZbelongV(Vec)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real vector is an Integer vector
    !!--++
    !!--++ Update: February - 2005
    !!
    Function ZbelongV(Vec) Result(belong)
       !---- Argument ----!
       real(kind=cp),   dimension(:), intent( in) :: Vec  ! Real vector
       logical                                    :: belong

       !---- Local variables ----!
       integer                               :: i
       real(kind=cp),   dimension(size(vec)) :: vect

       belong=.false.
       vect= abs(real(nint (vec))-vec)

       do i=1,size(vec)
          if (vect(i) > epss) return
       end do
       belong=.true.

       return
    End Function ZbelongV

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!---
    !!---- Subroutine Sort_Strings(Str)
    !!----
    !!----    Sort an array of string
    !!----
    !!---- Update: March - 2005
    !!
    Recursive Subroutine Sort_Strings(Str)
       !---- Argument ----!
       character(len=*), dimension(:), intent(in out) :: Str   ! String vector

       !---- Local variables ----!
       integer :: iq

       if (size(Str) > 1) then
          call Partition(Str, iq)
          call Sort_Strings(Str(:iq-1))
          call Sort_Strings(Str(iq:))
       end if

       return
    End Subroutine Sort_Strings

    !!----
    !!----  Subroutine Co_Prime_Vector(Vec,Cop,F)
    !!----
    !!----     Calculates the co-prime vector (cop) parallel to the input vector (v)
    !!----     It uses the list of the first thousand prime numbers.
    !!----
    !!----     copied from Nodal_Indices (Laue_Mod) in July 2013 (JRC)
    !!----
    !!----   Updated: January 2013
    !!
    Subroutine Co_Prime_Vector(Vec,Cop,f)
       !---- Arguments ----!
       integer, dimension(:), intent(in)  :: Vec    ! Input integer vector
       integer, dimension(:), intent(out) :: Cop    ! Ouput Co-prime vector
       integer,  optional,    intent(out) :: f      ! Common multiplicative factor

       !---- Local variables ----!
       integer                     :: i,j,max_ind,k,im,dimv,n

       cop=vec
       n=1
       if (present(f)) f=1
       max_ind=maxval(abs(cop))

       !> If the maximum value of the indices is 1 they are already coprimes
       if (max_ind <= 1) return

       !> Indices greater than 1
       dimv=size(vec)
       im=0
       do i=1,size(primes)
          if (primes(i) > max_ind) then  !primes is an array within this module
             im=i
             exit
          end if
       end do
       if (im == 0) return

       do_p: do i=1,im
          k=primes(i)
          do
             do j=1,dimv
                if( mod(cop(j),k) /= 0) cycle do_p
             end do
             n=n*k
             cop=cop/k
          end do
       end do do_p

       if (present(f)) f=n

       return
    End Subroutine Co_Prime_vector

    !!--++
    !!--++ Subroutine Diagonalize_Herm(Array,n,e_val,e_vect)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Diagonalize Hermitian matrices.
    !!--++    The eigen_values E_val are sorted in descending order. The columns
    !!--++    of E_vect are the corresponding eigenvectors.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Diagonalize_Herm(Array,N,E_val,E_vect)
       !---- Arguments ----!
       complex,           dimension(:,:), intent( in)  :: Array   ! Input array NxN
       integer,                           intent( in)  :: N       ! Dimension of array
       real(kind=cp),     dimension(:),   intent(out)  :: E_val   ! Eigenvalues
       complex, optional, dimension(:,:), intent(out)  :: E_vect  ! Eigenvectors

       !---- Local variables ----!
       real(kind=cp),  dimension(2*n,2*n)   :: aux
       real(kind=cp),  dimension(2*n)       :: e,d
       integer                              :: nn


       !> Init
       call init_err_mathgen()
       e_val=0.0_cp
       if (present(e_vect)) e_vect=0.0_cp

       !> Check
       if (n > size(Array,1) .or. n > size(Array,2)) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Diagonalize_HERM: Error in dimension of input matrix: Array(m,m) with m < n "
          return
       end if

       nn=2*n
       aux(  1:n ,  1:n ) =  real(array(1:n ,1:n))   !      (  U   V )
       aux(n+1:nn,n+1:nn) =  real(array(1:n ,1:n))   !   M=(          ),   A = U + i V
       aux(n+1:nn,  1:n ) = aimag(array(1:n ,1:n))   !      ( -V   U )
       aux(  1:n ,n+1:nn) =-aimag(array(1:n ,1:n))   !

       if (present(E_vect)) then
          call tred2(aux,nn,d,e)
          call tqli2(d,e,nn,aux)
          call eigsrt(d,aux,nn,1)
          e_vect(1:n,1:n)=cmplx(aux(1:n,1:nn:2),aux(n+1:nn,1:nn:2))
       else
          call tred1(aux,nn,d,e)
          call tqli1(d,e,nn)
          call eigsrt(d,aux,nn,0)
       end if
       e_val(1:n)=d(1:nn:2)

       return
    End Subroutine Diagonalize_Herm

    !!--++
    !!--++ Subroutine Diagonalize_Symm(Array,n,e_val,e_vect)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Diagonalize symmetric matrices
    !!--++    The eigen_values E_val are sorted in descending order. The columns
    !!--++    of E_vect are the corresponding eigenvectors.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Diagonalize_Symm(Array,N,E_Val,E_vect)
       !---- Arguments ----!
       real(kind=cp),           dimension(:,:), intent( in)  :: Array    ! Input symmetric array
       integer,                                 intent( in)  :: n        ! Dimension of Array
       real(kind=cp),           dimension(:),   intent(out)  :: E_val    ! Eigenvalues
       real(kind=cp), optional, dimension(:,:), intent(out)  :: E_vect   ! Eigenvectors

       !---- Local variables ----!
       real(kind=cp),        dimension(n,n)   :: aux
       real(kind=cp),        dimension(n)     :: e

       !> Init
       call init_err_mathgen()
       e_val=0.0_cp
       if (present(e_vect)) e_vect=0.0_cp

       !> Check
       if (n > size(Array,1) .or. n > size(Array,2)) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Diagonalize_SYMM: Error in dimension of input matrix: Array(m,m) with m < n "
          return
       end if

       aux=array(1:n,1:n)
       if (present(E_vect)) then
          call tred2(aux,n,E_val,e)
          call tqli2(E_val,e,n,aux)
          call eigsrt(E_val,aux,n,1)
          e_vect(1:n,1:n)=aux
       else
          call tred1(aux,n,E_val,e)
          call tqli1(E_val,e,n)
          call eigsrt(E_val,aux,n,0)
       end if

       return
    End Subroutine Diagonalize_Symm

    !!--++
    !!--++ Subroutine Eigsrt(d,v,n,io)
    !!--++    real(kind=cp), dimension(:),   intent(in out) :: d
    !!--++    real(kind=cp), dimension(:,:), intent(in out) :: v
    !!--++    integer,                       intent (in)    :: n
    !!--++    integer,                       intent (in)    :: io
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for sorting eigenvalues in d(n) and eigenvectors
    !!--++    in columns of v(n,n). Sorts d(n) in descending order and
    !!--++    rearranges v(n,n) correspondingly. The method is the straight
    !!--++    insertion. If io=0 order  only the eigenvalues are treated.
    !!--++    Adapted from Numerical Recipes. Valid for hermitian matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Eigsrt(d,v,n,io)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in out) :: d
       real(kind=cp), dimension(:,:), intent(in out) :: v
       integer,                       intent(in)     :: n
       integer,                       intent(in)     :: io

       !---- Local Variables ----!
       integer          :: i,j,k
       real(kind=cp)    :: p

       do i=1,n-1
          k=i
          p=d(i)
          do j=i+1,n
             if (d(j) >= p) then
                k=j
                p=d(j)
             end if
          end do
          if (k /= i) then
             d(k)=d(i)
             d(i)=p
             if (io == 1) then
                do j=1,n
                   p=v(j,i)
                   v(j,i)=v(j,k)
                   v(j,k)=p
                end do
             end if
          end if
       end do

       return
    End Subroutine Eigsrt

    !!----
    !!---- Subroutine First_Derivative(X, Y, N, D2Y, D1Y)
    !!----
    !!----    Calculate the First derivate values of the N points
    !!----
    !!---- Update: January - 2006
    !!
    Subroutine First_Derivative(X,Y,N,D2Y,D1Y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: X   ! Vector X
       real(kind=cp), dimension(:), intent(in)  :: Y   ! Vector Yi=F(Xi)
       integer ,                    intent(in)  :: N   ! Dimension of Vector N (Number of Points)
       real(kind=cp), dimension(:), intent(in)  :: d2y ! Vector containing the second derivative
       real(kind=cp), dimension(:), intent(out) :: d1y ! Vector contaning the First derivatives at the given points

       !---- Local Variables ----!
       integer       :: i
       real(kind=cp) :: step, x0, y0, y1, y2

       do i=1,n
          if (i /= n) then
             step = x(i+1)-x(i)
          end if
          x0 = x(i) - step/2.0
          y1=splint(x,y, d2y, n, x0)
          x0 = x(i) + step/2
          y2=splint(x,y, d2y, n, x0)
          y2 = y0
          d1y(i) = (y2 - y1) / step
       end do

       return
    End Subroutine First_Derivative

    !!----
    !!---- Subroutine In_Sort(Vec,N,Ind_i,Ind_f)
    !!--<<
    !!----    Subroutine to order in ascending mode the integer vector VEC.
    !!----    The input value N is the number of items to be ordered in VEC.
    !!----    The vector Ind_i is the initial pointer to Vec (coming from a previous call)
    !!----    The final pointer holding the order of items.
    !!-->>
    !!----
    !!---- Update: November - 2008
    !!
    Subroutine In_Sort(Vec,N,Ind_i,Ind_f)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: Vec    ! Integer array to be sorted
       integer,               intent(in) :: n      ! Number items in the array
       integer, dimension(:), intent(in) :: Ind_i  ! Initial pointer from a previous related call
       integer, dimension(:), intent(out):: Ind_f  ! Final pointer doing the sort of id

       !--- Local Variables ----!
       integer :: i,j,k,l,m
       integer, dimension(:),allocatable :: it

       !> Init
       Ind_f=0

       l=minval(Vec)
       m=maxval(Vec)
       l=l-1
       m=m-l

       allocate(it(m))
       it(1:m)=0

       do i=1,n
          j=vec(ind_i(i))-l
          it(j)=it(j)+1
       end do
       j=0
       do i=1,m
          k=j
          j=j+it(i)
          it(i)=k
       end do
       do i=1,n
          j=vec(ind_i(i))-l
          it(j)=it(j)+1
          j=it(j)
          ind_f(j)=ind_i(i)
       end do

       return
    End Subroutine In_Sort

    !!----
    !!---- Subroutine Init_Err_Mathgen()
    !!----
    !!----    Initialize the errors flags in CFML_Math_General
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_MathGen()

       ERR_MathGen=.false.
       ERR_MathGen_Mess=" "

       return
    End Subroutine Init_Err_MathGen

    !!----
    !!---- Subroutine Invert_Matrix(Array,Array_Inv,Singular,Perm)
    !!--<<
    !!----    Subroutine to invert a real matrix using LU decomposition.
    !!----    In case of singular matrix (singular=.true.) instead of the inverse
    !!----    matrix, the subroutine provides the LU decomposed matrix as used
    !!----    in Numerical Recipes.
    !!----    The input matrix is preserved and its inverse (or its LU decomposition)
    !!----    is provided in "b". The optional argument "perm" holds the row permutation
    !!----    performed to obtain the LU decomposition.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Invert_Matrix(Array,Array_Inv,Singular,Perm)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:),  intent(in ) :: Array      ! Input array
       real(kind=cp), dimension(:,:),  intent(out) :: Array_Inv  ! Inverse array
       logical,                        intent(out) :: singular   ! Flag for singular array
       integer, dimension(:),optional, intent(out) :: perm       ! Row permutation

       !---- Local variables ----!
       integer                                               :: i,n
       integer,       dimension(size(array,1))               :: indx
       real(kind=cp)                                         :: d, det
       real(kind=cp), dimension(size(array,1),size(array,1)) :: lu

       n=size(array,1)
       lu=array(1:n,1:n)

       call LU_Decomp(lu,d,singular,indx)
       if (present(perm)) perm(1:n)=indx(1:n)

       if (singular) then
          array_inv=lu
          return
       else
          det=0.0
          do i=1,n
             d=d*sign(1.0_cp,lu(i,i))
             det=det + log(abs(lu(i,i)))
          end do
          det=d*exp(det)
          if (abs(det) <= 1.0e-36) then
             singular=.true.
             array_inv=lu
             return
          end if
       end if

       array_inv=0.0
       do i=1,n
          array_inv(i,i)=1.0
          call LU_backsub(lu,indx,array_inv(:,i))
       end do

       return
    End Subroutine Invert_Matrix

    !!----
    !!---- Subroutine LU_Backsub(Array,indx,Vec)
    !!--<<
    !!----    Adapted from Numerical Recipes.
    !!----    Solves the set of N linear equations A  X = B. Here the N x N matrix A is input,
    !!----    not as the original matrix A, but rather as its LU decomposition, determined
    !!----    by the routine LU_DECOMP. INDX is input as the permutation vector of length N
    !!----    returned by LU_DECOMP. B is input as the right-hand-side vector B,
    !!----    also of length N, and returns with the solution vector X.
    !!----    A and INDX are not modified by this routine and can be left in place for successive calls
    !!----    with different right-hand sides B. This routine takes into account the possibility that B will
    !!----    begin with many zero elements, so it is efficient for use in matrix inversion.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine LU_Backsub(Array,indx,Vec)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in)     :: Array  ! Input array(N,N)
       integer,         dimension(:), intent(in)     :: indx   ! Permutation vector
       real(kind=cp),   dimension(:), intent(in out) :: Vec    ! Vector

       !---- Local Variables ----!
       integer       :: i,ii,ll,n
       real(kind=cp) :: summ

       n=size(array,1)
       ii=0              !When ii is set to a positive value, it will become the index
       do i=1,n          !of the first nonvanishing element of b. We now do
          ll=indx(i)     !the forward substitution. The only new wrinkle is to
          summ=vec(ll)     !unscramble the permutation as we go.
          vec(ll)=vec(i)
          if (ii /= 0) then
             summ=summ-dot_product(array(i,ii:i-1),vec(ii:i-1))
          else if(summ /= 0.0) then   !A nonzero element was encountered, so from now on
             ii=i                       !we will have to do the dot product above.
          end if
          vec(i)=summ
       end do

       do i=n,1,-1       !Now we do the backsubstitution
          vec(i) = (vec(i)-dot_product(array(i,i+1:n),vec(i+1:n)))/array(i,i)
       end do

       return
    End Subroutine LU_Backsub

    !!----
    !!---- Subroutine LU_Decomp(Array,d,singular,indx)
    !!--<<
    !!----    Subroutine to make the LU decomposition of an input matrix A.
    !!----    The input matrix is destroyed and replaced by a matrix containing
    !!----    in its upper triangular part (plus diagonal) the matrix U. The
    !!----    lower triangular part contains the nontrivial part (Lii=1) of matrix L.
    !!----    The output is rowwise permutation of the initial matrix. The vector INDX
    !!----    recording the row permutation. D is output as +/-1 depending on whether
    !!----    the number of row interchanges was even or odd, respectively.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine LU_Decomp(Array,d,singular,indx)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: array     ! Input array
       real(kind=cp),                 intent(out)    :: d         ! Information about permutation rows
       logical,                       intent(out)    :: singular
       integer,  dimension(:), intent(out), optional :: indx

       !---- Local variables ----!
       real(kind=cp), parameter               :: VTINY = 1.0e-7_sp    !A small number.
       real(kind=cp), dimension(size(array,1)):: vv                   !vv stores the implicit scaling of each row.
       integer                                :: j,imax,n

       !> Init
       d=0.0_cp
       singular=.false.

       n=size(array,1)
       if (present(indx)) then
          do j=1,n
             indx(j)=j
          end do
       end if

       d=1.0                      !No row interchanges yet.
       vv=maxval(abs(array),dim=2)    !Loop over rows to get the implicit scaling information.
       if (any(abs(vv) <= vtiny)) then   !There is a row of zeros.
          singular=.true.
          return
       end if
       vv=1.0_sp/vv     !Save the scaling.
       do j=1,n
          !imax=(j-1)+imaxloc(vv(j:n)*abs(array(j:n,j)))   !Find the pivot row.
          imax=(j-1)+maxloc(vv(j:n)*abs(array(j:n,j)),dim=1)
          if (j /= imax) then                            !Do we need to interchange rows?
             call swap(array(imax,:),array(j,:))         !Yes, do so...
             d=-d                                        !...and change the parity of d.
             vv(imax)=vv(j)                              !Also interchange the scale factor.
          end if
          if (present(indx)) indx(j)=imax
          if (abs(array(j,j)) <= vtiny) then !If the pivot element is zero the matrix is singular.
             array(j,j)=vtiny                !(at least to the precision of the algorithm)
             singular=.true.             !For some applications on singular matrices,
             !return                      !it is desirable to substitute vtiny for zero.
          end if                         !This is actually the present case
          array(j+1:n,j)=array(j+1:n,j)/array(j,j)                                    !Divide by the pivot element.
          array(j+1:n,j+1:n)=array(j+1:n,j+1:n)-outerprod(array(j+1:n,j),array(j,j+1:n))  !Reduce remaining submatrix.
       end do

       return
    End Subroutine LU_Decomp

    !!--++
    !!--++ Subroutine Masked_Swap_R(A,B,Mask)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b if mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_R(A,B,Mask)
       !---- Arguments ----!
       real(kind=cp), intent(in out) :: a       ! Value 1
       real(kind=cp), intent(in out) :: b       ! Value 2
       logical,           intent(in) :: mask    ! Mask to apply

       !---- Local Variables ----!
       real(kind=cp) :: swp

       if (mask) then
          swp=a
          a=b
          b=swp
       end if

       return
    End Subroutine Masked_Swap_R

    !!--++
    !!--++ Subroutine Masked_Swap_Rm(Array1,Array2,Mask)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b where mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_Rm(Array1,Array2,Mask)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: array1  ! Array 1
       real(kind=cp), dimension(:,:), intent(in out) :: array2  ! Array2
       logical,       dimension(:,:), intent(in)     :: mask    ! Mask to apply

       !---- Local variables ----!
       real(kind=cp), dimension(size(array1,1),size(array1,2)) :: swp

       where (mask)
          swp=array1
          array1=array2
          array2=swp
       end where

       return
    End Subroutine Masked_Swap_Rm

    !!--++
    !!--++ Subroutine Masked_Swap_Rv(Vec1,Vec2,Mask)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of VEC1 and VEC2 where mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_Rv(Vec1,Vec2,Mask)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out) :: Vec1  ! Vector 1
       real(kind=cp), dimension(:), intent(in out) :: Vec2  ! Vector 2
       logical,       dimension(:), intent(in)     :: mask  ! Mask to apply

       !---- Local variables ----!
       real(kind=cp), dimension(size(vec1)) :: swp

       where (mask)
          swp=vec1
          vec1=vec2
          vec2=swp
       end where

       return
    End Subroutine Masked_Swap_Rv

    !!----
    !!---- Subroutine Median_QS(Vec, n, Xmed)
    !!---- integer, intent(in)                :: n
    !!---- real, intent(in out), dimension(:) :: x
    !!---- real, intent(out)                  :: xmed
    !!----
    !!---- Subroutine calculating the median of a real array using a partially
    !!---- Quicsort method. On exit, the array x is partially ordered.
    !!---- Based in Alan Miller's median.f90 code.
    !!----
    !!
    Subroutine Median_QS(Vec, N, Xmed)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in out) :: Vec      ! Vector
       integer,                     intent(in)     :: n        ! Dimension of Vector
       real(kind=cp),               intent(out)    :: xmed     ! Mean value of Vector

       !---- Local Variables ----!
       integer       :: hi, lo, nby2, nby2p1, mid, i, j, k
       real(kind=cp) :: temp, xhi, xlo, xmax, xmin
       logical       :: odd

       nby2 = n / 2
       nby2p1 = nby2 + 1
       odd = .true.

       !>  HI & LO are position limits encompassing the median.
       if (n == 2 * nby2) odd = .false.
       lo = 1
       hi = n
       if (n < 3) then
          if (n < 1) then
             xmed = 0.0
             return
          end if
          xmed = vec(1)
          if (n == 1) return
          xmed = 0.5*(xmed + vec(2))
          return
       end if

       !>  Find median of 1st, middle & last values.
       do
          mid = (lo + hi)/2
          xmed = vec(mid)
          xlo  = vec(lo)
          xhi  = vec(hi)
          if (xhi < xlo) then          ! swap xhi & xlo
             temp = xhi
             xhi = xlo
             xlo = temp
          end if
          if (xmed > xhi) then
             xmed = xhi
          else if (xmed < xlo) then
             xmed = xlo
          end if

          !> The basic quicksort algorithm to move all values <= the sort key (XMED)
          !> to the left-hand end, and all higher values to the other end.
          i = lo
          j = hi
          do
              do
                 if (vec(i) >= xmed) exit
                 i = i + 1
              end do
              do
                 if (vec(j) <= xmed) exit
                 j = j - 1
              end do
              if (i < j) then
                 temp = vec(i)
                 vec(i) = vec(j)
                 vec(j) = temp
                 i = i + 1
                 j = j - 1
                 !> Decide which half the median is in.
                 if (i <= j) cycle
              end if
              exit
          end do

          if (.not. odd) then
             if (j == nby2 .and. i == nby2p1) then
                !> Special case, N even, J = N/2 & I = J + 1, so the median is
                !> between the two halves of the series.   Find max. of the first
                !> half & min. of the second half, then average.
                xmax = vec(1)
                do k = lo, j
                   xmax = max(xmax, vec(k))
                end do
                xmin = vec(n)
                do k = i, hi
                   xmin = min(xmin, vec(k))
                end do
                xmed = 0.5*(xmin + xmax)
                return
             end if

             if (j < nby2) lo = i
             if (i > nby2p1) hi = j
             if (i /= j) then
                if (lo < hi - 1) cycle
                exit
             end if
             if (i == nby2) lo = nby2
             if (j == nby2p1) hi = nby2p1

          else
             if (j < nby2p1) lo = i
             if (i > nby2p1) hi = j
             if (i /= j) then
                if (lo < hi - 1) cycle
                exit
             end if

             !> test whether median has been isolated.
             if (i == nby2p1) return
          end if

          if (lo < hi - 1) cycle
          exit
       end do

       if (.not. odd) then
          xmed = 0.5*(vec(nby2) + vec(nby2p1))
          return
       end if
       temp = vec(lo)
       if (temp > vec(hi)) then
          vec(lo) = vec(hi)
          vec(hi) = temp
       end if
       xmed = vec(nby2p1)

       return
    End Subroutine Median_QS

    !!--++
    !!--++ Subroutine Partition(A, marker)
    !!--++
    !!--++    (Private)
    !!--++    Utilised by Sort_Strings.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Partition(A, Marker)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in out) :: A         ! String
       integer,                        intent(   out) :: marker    ! Integer

       !---- Local variables ----!
       integer                  :: i, j
       character(len=len(A(1))) :: temp
       character(len=len(A(1))) :: x      ! pivot point

       x = A(1)
       i= 0
       j= size(A) + 1

       do
          j = j-1
          do
             if (A(j) <= x) exit
             j = j-1
          end do
          i = i+1
          do
             if (A(i) >= x) exit
             i = i+1
          end do
          if (i < j) then
             !> exchange A(i) and A(j)
             temp = A(i)
             A(i) = A(j)
             A(j) = temp
          else if (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          end if
       end do

       return
    End Subroutine Partition

    !!----
    !!---- Subroutine Points_In_Line2D(X1, XN, N, XP)
    !!----
    !!----    The routine calculate N points belonging to the line defined
    !!----    by X1 and Xn with equal distance between them.
    !!----    XP (2,N) contains X1,X2,.....,XN points.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Points_In_Line2D(X1, XN, N, XP)
       !---- Arguments ----!
       real(kind=cp), dimension(2),   intent(in)  :: X1   ! Point1 in 2D
       real(kind=cp), dimension(2),   intent(in)  :: XN   ! PointN in 2D
       integer,                       intent(in)  :: N    ! Number of Total points
       real(kind=cp), dimension(:,:), intent(out) :: XP   ! List of points

       !---- Local Variables ----!
       integer :: i
       real(kind=cp)    :: ml,bl,dl,t
       real(kind=cp)    :: a,b,c,d
       real(kind=cp)    :: xa,xb

       xp=0.0

       if (n <= 1) return

       !---- Calculating the distance between two points to
       !---- eliminate rare considerations as the same point
       dl=sqrt( (xn(1)-x1(1))**2 + (xn(2)-x1(2))**2 )
       if (dl <= 0.0001) return

       !---- When N=2 is trivial case ----!
       if (n == 2) then
          xp(:,1)=x1
          xp(:,2)=xn
          return
       end if

       !---- Case 1: Y=cte ----!
       !Xn(2) and X1(2) are equal, then we have a line  with Y=cte
       if (abs(xn(2)-x1(2)) <= 0.0001) then
          dl=abs(xn(1)-x1(1))
          d=dl/real(n-1)
          xp(:,1)=x1
          if (xn(1) > x1(1)) then
             do i=2,n-1
                xp(1,i)=xp(1,i-1)+d
                xp(2,i)=xp(2,1)
             end do
          else
             do i=2,n-1
                xp(1,i)=xp(1,i-1)-d
                xp(2,i)=xp(2,1)
             end do
          end if
          xp(:,n)=xn

          return
       end if

       !---- Case 2: X=cte ----!
       !Xn(1) - X1(1) are equal, then we have a line with X=cte
       if (abs(xn(1)-x1(1)) <= 0.0001) then
          dl=abs(xn(2)-x1(2))
          d=dl/real(n-1)
          xp(:,1)=x1
          if (xn(2) > x1(2)) then
             do i=2,n-1
                xp(1,i)=xp(1,1)
                xp(2,i)=xp(2,i-1)+d
             end do
          else
             do i=2,n-1
                xp(1,i)=xp(1,1)
                xp(2,i)=xp(2,i-1)-d
             end do
          end if
          xp(:,n)=xn

          return
       end if

       !---- Case 3: General case ----!
       ml=(x1(2)-xn(2))/(x1(1)-xn(1))
       bl=x1(2) - (ml * x1(1))

       !---- Distance between X1 and XN ----!
       dl=sqrt( (xn(1)-x1(1))**2 + (xn(2)-x1(2))**2 )

       !---- Creating the list ----!
       a=ml**2 + 1.0
       b=2.0 *( ml*(bl-x1(2)) -x1(1) )

       xp(:,1)=x1
       do i=2,n-1
          t=(dl**2)*((real(i-1)/real(n-1))**2)
          c=(x1(2)-bl)**2 + x1(1)**2 - t

          xa=(-b + sqrt(b**2 - 4.0*a*c))/(2.0*a)
          xb=(-b - sqrt(b**2 - 4.0*a*c))/(2.0*a)
          if (x1(1) <= xa .and. xa <= xn(1)) then
             xp(1,i)=xa
             xp(2,i)=ml*xa+bl
          else
             xp(1,i)=xb
             xp(2,i)=ml*xb+bl
          end if
       end do
       xp(:,n)=xn

       return
    End Subroutine Points_In_Line2D

    !!--++
    !!--++ Subroutine Rank_dp(a,tol,r)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rank_dp(Array,tol,r)
       !---- Arguments ----!
       real(kind=dp), dimension(:,:),intent( in)      :: Array     ! Rectangular Matrix
       real(kind=dp),                intent( in)      :: tol       ! Tolerance
       integer,                      intent(out)      :: r         ! Rank of Array

       !---- Local Variables ----!
       real(kind=dp), dimension(size(array,1),size(array,2))  :: u
       real(kind=dp), dimension(size(array,2))                :: w
       real(kind=dp), dimension(size(array,2),size(array,2))  :: v
       integer                                                :: i

       u=array
       call svdcmp(u,w,v)

       if (ERR_MathGen) then
          r=0
       else
          r=0
          do i=1,size(array,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_dp

    !!--++
    !!--++ Subroutine Rank_sp(a,tol,r)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rank_sp(Array,tol,r)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:),intent( in)      :: array   ! Rectagular matrix
       real(kind=sp),                intent( in)      :: tol     ! Tolerance
       integer,                      intent(out)      :: r       ! Rank of Array

       !---- Local variables ----!
       real(kind=sp), dimension(size(array,1),size(array,2))  :: u
       real(kind=sp), dimension(size(array,2))            :: w
       real(kind=sp), dimension(size(array,2),size(array,2))  :: v
       integer :: i

       u=array
       call svdcmp(u,w,v)
       if (ERR_MathGen) then
          r=0
       else
          r=0
          do i=1,size(array,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_sp

    !!----
    !!---- Subroutine Second_Derivative(x, y, n, d2y)
    !!----
    !!----    Calculate the second derivate of N Points
    !!----
    !!---- Update: January - 2006
    !!
    Subroutine Second_Derivative(x,y,n,d2y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: x     ! X array
       real(kind=cp), dimension(:), intent(in)  :: y     ! Y array Yi=F(Xi)
       integer ,                    intent(in)  :: n     ! Dimension of X and Y
       real(kind=cp), dimension(:), intent(out) :: d2y   ! Vector conytining the second derivatives en Xi

       !---- Local Variables ----!
       integer                     :: i, k
       real(kind=cp), dimension(n) :: u
       real(kind=cp)               :: yp1, ypn, sig, p, qn, un

       yp1=(y(2) - y(1))   / (x(2) - x(1))     ! derivative at point 1
       ypn=(y(n) - y(n-1)) / (x(n) - x(n-1))   ! derivative at point n

       d2y(1)=-0.5
       u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)

       do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*d2y(i-1)+2.0
          d2y(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))  &
               /(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do

       qn=0.5
       un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       d2y(n)=(un-qn*u(n-1))/(qn*d2y(n-1)+1.0)
       do k=n-1,1,-1
          d2y(k)=d2y(k)*d2y(k+1)+u(k)
       end do

       return
    End Subroutine Second_Derivative

    !!----
    !!---- Subroutine Set_Epsg(Neweps)
    !!----
    !!---- Sets global EPSS to the value "neweps".
    !!---- If Not argument is given, then EPSS is set to default value
    !!----
    !!---- EPSS (default) = 1.0E-5_sp
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Set_Epsg(Neweps)
       !---- Arguments ----!
       real(kind=cp), optional, intent( in) :: neweps

       if (present(neweps)) then
          epss=neweps
       else
          epss=1.0E-5_sp
       end if

       return
    End Subroutine Set_Epsg

    !!----
    !!---- Subroutine SmoothingVec(Y, N, NIter, Ys)
    !!----    real(kind=cp),    dimension(:),           intent(in out) :: Y      !  In Out-> Array to be smoothed
    !!----    integer,                                  intent(in)     :: N      !  In -> Number of points
    !!----    integer,                                  intent(in)     :: NIter  !  In -> Number of iterations
    !!----    real(kind=cp),    dimension(:), optional, intent(out)    :: datY   !  Out-> Array smoothed
    !!----
    !!----    Procedure to smooth a vector values
    !!----
    !!---- Update: January - 2006
    !!
    Subroutine SmoothingVec(Vec, N, Iter, VecS)
       !---- Arguments ----!
       real(kind=cp),dimension(:),            intent(in out) :: Vec    ! Vector to be smoothed
       integer,                               intent(in)     :: n      ! Number of Points
       integer,                               intent(in)     :: iter   ! Number of iterations
       real(kind=cp),dimension(:), optional,  intent(out)    :: VecS   ! Vector smoothed if present

       !---- Local Variables ----!
       integer                     :: n1, n2
       integer                     :: i, j
       real(kind=cp), dimension (n):: Y,datYs


       n1 = 4
       n2 = n-3

       y=vec(1:n)

       do j = 1 ,iter
          datYs(n1-1)=((Y(n1-2)+Y(n1))*10.0+(Y(n1-3)+Y(n1+1))*5.0+Y(n1+2))/31.0
          datYs(n1-2)=((Y(n1-3)+Y(n1-1))*10.0+Y(n1)*5.0+Y(n1+1))/26.0
          datYs(n1-3)=(Y(n1-2)*10.0+Y(n1-1)*5.0+Y(n1))/16.0

          do i=n1,n2
             datYs(i)=(Y(i-3)+Y(i+3)+5.0*(Y(i-2)+Y(i+2))+10.0*(Y(i-1)+Y(i+1)))/ 32.0
          end do

          datYs(n2+1)=((Y(n2+2)+Y(n2))*10.0+(Y(n2+3)+Y(n2-1))*5.0+Y(n2-2))/31.0
          datYs(n2+2)=((Y(n2+3)+Y(n2+1))*10.0+Y(n2)*5.0+Y(n2-1))/26.0
          datYs(n2+3)=(Y(n2+2)*10.0+Y(n2+1)*5.0+Y(n2))/16.0

          y=datYs
       end do

       !> Output procedure
       if (present(VecS)) then
          VecS(1:n) = datYs(1:n)
       else
          Vec(1:n) = datYs(1:n)
       end if

       return
    End Subroutine SmoothingVec

    !!--++
    !!--++ Subroutine Sort_I(Vec,N,Indx)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sort a vector such the vec(indx(j)) is in ascending
    !!--++    order for j=1,2,...,N.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sort_I(vec,n,indx)
       !---- Arguments ----!
       integer, dimension(:), intent(in ) :: vec     ! Vector to be sorted
       integer              , intent(in ) :: n       ! Dimension of the vector
       integer, dimension(:), intent(out) :: indx    ! Index for sort values

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer                      :: i,indxt,ir,itemp,j,jstack,k,l
       integer                      :: a

       call init_Err_MathGen()
       do j=1,n
          indx(j)=j
       end do

       istack=0
       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=Vec(indxt)
                do i=j-1,1,-1
                   if (Vec(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (Vec(indx(l+1)) > Vec(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (Vec(indx(l)) > Vec(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (Vec(indx(l+1)) > Vec(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=Vec(indxt)
             do
                i=i+1
                if (Vec(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (Vec(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                ERR_MathGen=.true.
                ERR_MathGen_Mess=" NSTACK too small in SORT"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Subroutine Sort_I

    !!--++
    !!--++ Subroutine Sort_R(vec,n,indx)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sort a vector such the vec(indx(j)) is in ascending
    !!--++    order for j=1,2,...,N.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sort_R(Vec,n,indx)
       !---- Arguments ----!
       real(kind=cp),dimension(:), intent(in) :: Vec    ! Vector to be sorted
       integer,                    intent(in) :: n      ! Dimension of Vector
       integer,      dimension(:), intent(out):: indx   ! Index

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer                      :: i,indxt,ir,itemp,j,jstack,k,l
       real(kind=cp)                :: a

       call init_Err_MathGen()
       do j=1,n
          indx(j)=j
       end do

       istack=0
       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=Vec(indxt)
                do i=j-1,1,-1
                   if (Vec(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (Vec(indx(l+1)) > Vec(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (Vec(indx(l)) > Vec(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (Vec(indx(l+1)) > Vec(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=Vec(indxt)
             do
                i=i+1
                if (Vec(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (Vec(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                ERR_MathGen=.true.
                ERR_MathGen_Mess=" NSTACK too small in SORT"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Subroutine Sort_R

    !!----
    !!---- SUBROUTINE SPH_JN
    !!----
    !!----    Compute spherical Bessel functions jn(x) and their derivatives
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Sph_Jn(n,x,nm,jn,djn)
       !---- Arguments ----!
       integer,                       intent(in)  :: n   !Order of jn(x) (n=0,1,2,3,...)
       real(kind=dp),                 intent(in)  :: x   !Argument of jn(x)
       integer,                       intent(out) :: nm  !Highest order computed
       real(kind=dp), dimension(0:n), intent(out) :: jn  !array with spherical Bessel functions jn(x)
       real(kind=dp), dimension(0:n), intent(out) :: djn !array with derivatives jn'(x)

       !---- Local variables ----!
       integer       :: k,m
       real(kind=dp) :: sa,sb, f,f1,f0, cs

       nm=n
       if (abs(x) <= 1.0e-30_dp) then
          do k=0,n
             jn(k) = 0.0_dp
             djn(k)= 0.0_dp
          end do
          jn(0)=1.0_dp
          djn(1)=1.0_dp/3.0_dp
          return
       end if

       jn(0)=sin(x)/x
       jn(1)=(jn(0)-cos(x))/x
       if (n >= 2) then
          sa=jn(0)
          sb=jn(1)
          m=start1(x,200)
          if (m < n) then
             nm=m
          else
             m=start2(x,n,15)
          end if
          f0=0.0_dp
          f1=1.0e-100_dp
          do k=m,0,-1
             f=(2.0_dp*k+3.0_dp)*f1/x-f0
             if (k <= nm) jn(k)=f
             f0=f1
             f1=f
          end do
          if (abs(sa) > abs(sb)) cs=sa/f
          if (abs(sa) <= abs(sb)) cs=sb/f0
          do k=0,nm
             jn(k)=cs*jn(k)
          end do
       end if

       djn(0)=(cos(x)-sin(x)/x)/x
       do k=1,nm
          djn(k)=jn(k-1)-(k+1.0_dp)*jn(k)/x
       end do

       return
    End Subroutine Sph_Jn

    !!----
    !!---- Subroutine Spline(x, y, n, yp1, ypn, y2)
    !!----
    !!----    Spline  N points at the given points
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Spline(x,y,n,yp1,ypn,y2)
       !---- Arguments ----!!
       real(kind=cp), dimension(:), intent(in)  :: x    ! Vector X
       real(kind=cp), dimension(:), intent(in)  :: y    ! Vector Y=F(Xi)
       integer ,                    intent(in)  :: n    ! Dimension of X
       real(kind=cp),               intent(in)  :: yp1  ! Derivate on Point 1
       real(kind=cp),               intent(in)  :: ypn  ! Derivate on Point N
       real(kind=cp), dimension(:), intent(out) :: y2   ! Vector containing Second derivatives at given points

       !---- Local Variables ----!
       integer                     :: i, k
       real(kind=cp), dimension(n) :: u
       real(kind=cp)               :: sig, p, qn, un

       if (yp1 > 1.0e+30) then
          y2(1)=0.0
          u(1)=0.0
       else
          y2(1)=-0.5
          u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       end if

       do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.0
          y2(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))  &
               /(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do
       if (ypn > 1.0e+30) then
          qn=0.0
          un=0.0
       else
          qn=0.5
          un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       end if
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
       do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
       end do

       return
    End Subroutine Spline

    !!--++
    !!--++ Subroutine Svdcmp_dp(a,w,v)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Given an M x N matrix A ,this routine computes its singular value decomposition,
    !!--++    A = U  W  VT . The matrix U replaces A on output. The diagonal matrix of
    !!--++    singular values W is output as the N-dimensional vector w. The NxN matrix V
    !!--++    (not the transpose VT )is output as v .
    !!--++
    !!--++    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Svdcmp_dp(a,w,v)
       !---- Arguments ----!
       real(kind=dp),dimension(:,:),intent(in out) ::a   ! A(m,n)
       real(kind=dp),dimension(:),  intent(   out) ::w   ! W(n)
       real(kind=dp),dimension(:,:),intent(   out) ::v   ! V(n,n)

       !---- Local variables ----!
       integer, parameter                  :: num_its=500
       integer                             :: i,its,j,k,l,m,n,nm
       real(kind=dp)                       :: anorm,c,f,g,h,s,scal,x,y,z
       real(kind=dp),dimension(size(a,1))  :: tempm
       real(kind=dp),dimension(size(a,2))  :: rv1,tempn

       m=size(a,1)
       n=size(a,2)
       call init_err_mathgen()
       if ( .not. (size(v,1) == n .and. size(v,2) == n .and. size(w) == n)) then
          ERR_MathGen = .true.
          ERR_MathGen_Mess = " => Physical dimensions of arguments in SVDcmp_dp are not compatible "
          return
       end if
       g=0.0_dp
       scal=0.0_dp
       do i=1,n
          l=i+1
          rv1(i)=scal*g
          g=0.0_dp
          scal=0.0_dp
          if (i <=m)then
             scal=sum(abs(a(i:m,i)))
             if ( abs(scal) > tiny(1.0_dp) ) then
                a(i:m,i)=a(i:m,i)/scal
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scal*a(i:m,i)
             end if
          end if
          w(i)=scal*g
          g=0.0_dp
          scal=0.0_dp
          if ((i <=m).and.(i /=n))then
             scal=sum(abs(a(i,l:n)))
             if ( abs(scal) > tiny(1.0_dp) ) then
                a(i,l:n)=a(i,l:n)/scal
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scal*a(i,l:n)
             end if
          end if
       end do
       anorm=maxval(abs(w)+abs(rv1))
       do i=n,1,-1
          if (i <n) then
             if ( abs(g) > tiny(1.0_dp) ) then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
             end if
             v(i,l:n)=0.0_dp
             v(l:n,i)=0.0_dp
          end if
          v(i,i)=1.0_dp
          g=rv1(i)
          l=i
       end do
       do i=min(m,n),1,-1
          l=i+1
          g=w(i)
          a(i,l:n)=0.0_dp
          if ( abs(g) > tiny(1.0_dp) ) then
             g=1.0_dp/g
             tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=a(i:m,i)*g
          else
             a(i:m,i)=0.0_dp
          end if
          a(i,i)=a(i,i)+1.0_dp
       end do
       do k=n,1,-1
          do its=1,num_its
             do l=k,1,-1
                nm=l-1
                if ((abs(rv1(l))+anorm)==anorm) exit
                if ((abs(w(nm))+anorm)==anorm) then
                   c=0.0_dp
                   s=1.0_dp
                   do i=l,k
                      f=s*rv1(i)
                      rv1(i)=c*rv1(i)
                      if ((abs(f)+anorm)==anorm)exit
                      g=w(i)
                      h=pythag(f,g)
                      w(i)=h
                      h=1.0_dp/h
                      c=(g*h)
                      s=-(f*h)
                      tempm(1:m)=a(1:m,nm)
                      a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                   end do
                   exit
                end if
             end do
             z=w(k)
             if (l ==k)then
                if (z <0.0_dp)then
                   w(k)=-z
                   v(1:n,k)=-v(1:n,k)
                end if
                exit
             end if
             if (its == num_its) then
                ERR_MathGen = .true.
                ERR_MathGen_Mess = " => SVDcmp_dp: convergence not reached ! "
                return
             end if
             x=w(l)
             nm=k-1
             y=w(nm)
             g=rv1(nm)
             h=rv1(k)
             f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dp*h*y)
             g=pythag(f,1.0_dp)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.0_dp
             s=1.0_dp
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=pythag(f,h)
                w(j)=z
                if ( abs(z) > tiny(1.0_dp) ) then
                   z=1.0_dp/z
                   c=f*z
                   s=h*z
                end if
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
             end do
             rv1(l)=0.0_dp
             rv1(k)=f
             w(k)=x
          end do
       end do

       return
    End Subroutine Svdcmp_dp

    !!--++
    !!--++ Subroutine Svdcmp_sp(a,w,v)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Given an M x N matrix A ,this routine computes its singular value decomposition,
    !!--++    A = U W VT . The matrix U replaces A on output. The diagonal matrix of
    !!--++    singular values W is output as the N-dimensional vector w. The N X N matrix V
    !!--++    (not the transpose VT )is output as v.
    !!--++    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Svdcmp_sp(a,w,v)
       !---- Arguments ----!
       real(kind=sp),dimension(:,:),intent(in out) :: a  ! A(m,n)
       real(kind=sp),dimension(:),  intent(   out) :: w  ! W(n)
       real(kind=sp),dimension(:,:),intent(   out) :: v  ! V(n,n)

       !---- Local variables ----!
       integer, parameter                 :: num_its=500
       integer                            :: i,its,j,k,l,m,n,nm
       real(kind=sp)                      :: anorm,c,f,g,h,s,scala,x,y,z
       real(kind=sp),dimension(size(a,1)) :: tempm
       real(kind=sp),dimension(size(a,2)) :: rv1,tempn


       m=size(a,1)
       n=size(a,2)
       call init_err_mathgen()
       if ( .not. (size(v,1) == n .and. size(v,2) == n .and. size(w) == n)) then
          ERR_MathGen = .true.
          ERR_MathGen_Mess = " => Physical dimensions of arguments in SVDcmp_sp are not compatible "
          return
       end if
       g=0.0
       scala=0.0
       do i=1,n                        !Householder reduction to bidiagonal form.
          l=i+1
          rv1(i)=scala*g
          g=0.0
          scala=0.0
          if (i <=m)then
             scala=sum(abs(a(i:m,i)))
             if (abs(scala) > tiny(1.0_sp))then
                a(i:m,i)=a(i:m,i)/scala
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scala*a(i:m,i)
             end if
          end if
          w(i)=scala*g
          g=0.0
          scala=0.0
          if ((i <=m).and.(i /=n))then
             scala=sum(abs(a(i,l:n)))
             if (abs(scala) > tiny(1.0_sp))then
                a(i,l:n)=a(i,l:n)/scala
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scala*a(i,l:n)
             end if
          end if
       end do
       anorm=maxval(abs(w)+abs(rv1))
       do i=n,1,-1                    ! Accumulation of right-hand transformations.
          if (i <n)then
             if (abs(g) > tiny(1.0_sp))then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g   !Double division to avoid possible underflow.
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
             end if
             v(i,l:n)=0.0
             v(l:n,i)=0.0
          end if
          v(i,i)=1.0
          g=rv1(i)
          l=i
       end do
       do i=min(m,n),1,-1  !Accumulation of left-hand transformations.
          l=i+1
          g=w(i)
          a(i,l:n)=0.0
          if (abs(g) > tiny(1.0_sp))then
             g=1.0_sp/g
             tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=a(i:m,i)*g
          else
             a(i:m,i)=0.0
          end if
          a(i,i)=a(i,i)+1.0_sp
       end do
       do k=n,1,-1            !Diagonalization of the idiagonal form:Loop over
          do its=1,num_its    !singular values,and over allowed iterations.
             do l=k,1,-1      !Test for splitting.
                nm=l-1        !Note that rv1(1)is always zero,so can never fall through bottom of loop.
                if ((abs(rv1(l))+anorm)==anorm) exit
                if ((abs(w(nm))+anorm)==anorm) then
                   c=0.0       ! Cancellation of rv1(l),if l >1 .
                   s=1.0
                   do i=l,k
                      f=s*rv1(i)
                      rv1(i)=c*rv1(i)
                      if ((abs(f)+anorm)==anorm)exit
                      g=w(i)
                      h=pythag(f,g)
                      w(i)=h
                      h=1.0_sp/h
                      c=(g*h)
                      s=-(f*h)
                      tempm(1:m)=a(1:m,nm)
                      a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                   end do
                   exit
                end if
             end do
             z=w(k)
             if (l ==k) then    !Convergence.
                if (z <0.0)then !Singular value is made nonnegative.
                   w(k)=-z
                   v(1:n,k)=-v(1:n,k)
                end if
                exit
             end if
             if (its == num_its) then
                ERR_MathGen = .true.
                ERR_MathGen_Mess = " => SVDcmp_sp: convergence not reached ! "
                return
             end if
             x=w(l)             !Shift from ottom 2-y-2 minor.
             nm=k-1
             y=w(nm)
             g=rv1(nm)
             h=rv1(k)
             f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
             g=pythag(f,1.0_sp)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.0  ! Next QR transformation:
             s=1.0
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=pythag(f,h)
                w(j)=z                 !Rotation can e arbitrary if z =0 .
                if (abs(z) > tiny(1.0_sp) )then
                   z=1.0_sp/z
                   c=f*z
                   s=h*z
                end if
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
             end do
             rv1(l)=0.0
             rv1(k)=f
             w(k)=x
          end do
       end do

       return
    End Subroutine Svdcmp_sp

    !!--++
    !!--++ Subroutine Swap_C(a,b)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_C(a,b)
       !---- Arguments ----!
       complex, intent(in out) :: a  ! Complex value
       complex, intent(in out) :: b  ! Complex value

       !---- Local variables ----!
       complex :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_C

    !!--++
    !!--++ Subroutine Swap_Cm(array1,array2)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of array1 and array2
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Cm(array1,array2)
       !---- Arguments ----!
       complex, dimension(:,:), intent(in out) :: array1
       complex, dimension(:,:), intent(in out) :: array2

       !---- Local variables ----!
       complex, dimension(size(array1,1),size(array1,2)) :: dum

       dum=array1
       array1=array2
       array2=dum

       return
    End Subroutine Swap_Cm

    !!--++
    !!--++ Subroutine Swap_Cv(Vec1,Vec2)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of Vec1 and Vec2
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Cv(Vec1,Vec2)
       !---- Arguments ----!
       complex, dimension(:), intent(in out) :: Vec1
       complex, dimension(:), intent(in out) :: Vec2

       !---- Local variables ----!
       complex, dimension(size(vec1)) :: dum

       dum=vec1
       vec1=vec2
       vec2=dum

       return
    End Subroutine Swap_Cv

    !!--++
    !!--++ Subroutine Swap_I(i,j)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of i and j
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_I(i,j)
       !---- Arguments ----!
       integer , intent(in out) :: i ! Integer
       integer , intent(in out) :: j ! Integer

       !---- Local variables ----!
       integer  :: dum

       dum=i
       i=j
       j=dum

       return
    End Subroutine Swap_I

    !!--++
    !!--++ Subroutine Swap_Im(array1,array2)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Im(Array1,Array2)
       !---- Arguments ----!
       integer, dimension(:,:), intent(in out) :: array1
       integer, dimension(:,:), intent(in out) :: array2

       !---- Local Variables ----!
       integer, dimension(size(array1,1),size(array1,2)) :: dum

       dum=array1
       array1=array2
       array2=dum

       return
    End Subroutine Swap_Im

    !!--++
    !!--++ Subroutine Swap_Iv(Vec1,Vec2)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of Vec1 and Vec2
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Iv(Vec1,Vec2)
       !---- Arguments ----!
       integer, dimension(:), intent(in out) :: Vec1
       integer, dimension(:), intent(in out) :: Vec2

       !---- Local Variables ----!
       integer, dimension(size(Vec1)) :: dum

       dum=Vec1
       Vec1=Vec2
       Vec2=dum

       return
    End Subroutine Swap_Iv

    !!--++
    !!--++ Subroutine Swap_R(A,B)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_R(A,B)
       !---- Arguments ----!
       real(kind=cp), intent(in out) :: a
       real(kind=cp), intent(in out) :: b

       !---- Local variables ----!
       real(kind=cp) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_R

    !!--++
    !!--++ Subroutine Swap_Rm(Array1,Array2)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of array1 and array2
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Rm(Array1,Array2)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: array1
       real(kind=cp), dimension(:,:), intent(in out) :: array2

       !---- Local variables ----!
       real(kind=cp), dimension(size(array1,1),size(array1,2)) :: dum

       dum=array1
       array1=array2
       array2=dum

       return
    End Subroutine Swap_Rm

    !!--++
    !!--++ Subroutine Swap_Rv(Vec1,Vec2)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of Vec1 and Vec2
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Rv(Vec1,Vec2)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out) :: Vec1
       real(kind=cp), dimension(:), intent(in out) :: Vec2

       !---- Local variables ----!
       real(kind=cp), dimension(size(Vec1)) :: dum

       dum=vec1
       vec1=vec2
       vec2=dum

       return
    End Subroutine Swap_Rv

    !!--++
    !!--++ Subroutine Tqli1(d,e,n)
    !!--++
    !!--++    (PRIVATE)
    !!--++    QL-algorithm with implicit shifts, to determine the eigenvalues
    !!--++    and eigenvectors of a real tridiagonal symmetric matrix, or of
    !!--++    a real symmetric matrix previously reduced by tred. D is a vector
    !!--++    with the diagonal elements of the tridiagonal matrix. on output
    !!--++    it returns the eigenvalues. the vector e inputs the subdiagonal
    !!--++    elements of the tridiagonal matrix, with E(1) arbitrary. on
    !!--++    output e is destroyed.
    !!--++    In TLQ1 only the eigenvalues are calculated
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tqli1(d,e,n)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out):: d, e ! d(np),e(np)
       integer,                     intent(in )   :: n

       !---- Local variables ----!
       integer, parameter :: NMAX_ITER=500
       integer      :: i, iter, l, m, mv
       real(kind=cp):: b, c, dd, f, g, p, r, s, comp

       call init_Err_MathGen()
       do i=2,n
          e(i-1)=e(i)
       end do
       e(n)=0.0
       do l=1,n
          iter=0
          do_g : do
             mv=n
             do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                comp= abs(e(m))+dd
                if (abs(comp-dd) <= ep_ss) then
                   mv=m
                   exit
                end if
             end do
             m=mv

             if (m /= l) then
                if (iter == NMAX_ITER) then
                   ERR_MathGen=.true.
                   ERR_MathGen_Mess=" Too many iterations in TQLI1"
                   exit
                end if

                iter=iter+1
                g=(d(l+1)-d(l))/(2.0*e(l))
                r=sqrt(g*g+1.0)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                   f=s*e(i)
                   b=c*e(i)
                   r=sqrt(f*f+g*g)
                   e(i+1)=r
                   if (abs(r)  <= ep_ss) then
                      d(i+1)=d(i+1)-p
                      e(m)=0.0
                      cycle do_g
                   end if
                   s=f/r
                   c=g/r
                   g=d(i+1)-p
                   r=(d(i)-g)*s+2.0*c*b
                   p=s*r
                   d(i+1)=g+p
                   g=c*r-b
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
                cycle do_g
             end if
             exit
          end do do_g
       end do

       return
    End Subroutine Tqli1

    !!--++
    !!--++ Subroutine Tqli2(d,e,n,z)
    !!--++
    !!--++    (PRIVATE)
    !!--++    QL-algorithm with implicit shifts, to determine the eigenvalues
    !!--++    and eigenvectors of a real tridiagonal symmetric matrix, or of
    !!--++    a real symmetric matrix previously reduced by tred. D is a vector
    !!--++    with the diagonal elements of the tridiagonal matrix. on output
    !!--++    it returns the eigenvalues. the vector e inputs the subdiagonal
    !!--++    elements of the tridiagonal matrix, with E(1) arbitrary. on
    !!--++    output e is destroyed.
    !!--++    The eigenvectors of the tridiagonal matrix are calculated in TLQ2
    !!--++    by providing the matrix Z  as the identity matrix on input. if the
    !!--++    eigenvectors of the matrix reduced by tred are required, then Z
    !!--++    is input as the matrix output of tred. in either cased, the k-th
    !!--++    column of Z returns the mormalized eigenvector corresponding to
    !!--++    D(k).
    !!--++
    !!--++  Update: February - 2005
    !!
    Subroutine Tqli2(d,e,n,z)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)
       integer,                       intent(in )    :: n
       real(kind=cp), dimension(:,:), intent(in out) :: z    ! z(np,np)

       !---- Local Variables ----!
       integer, parameter :: NMAX_ITER=500
       integer       :: i, iter, k, l, m, mv
       real(kind=cp) :: b, c, dd, f, g, p, r, s, comp

       call init_Err_MathGen()
       do i=2,n
          e(i-1)=e(i)
       end do

       e(n)=0.0
       do l=1,n
          iter=0
          do_g: do
             mv=n
             do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                comp= abs(e(m))+dd
                if (abs(comp-dd) <= ep_ss) then
                   mv=m
                   exit
                end if
             end do
             m=mv
             if (m /= l) then
                if (iter == NMAX_ITER) then
                   ERR_MathGen=.true.
                   ERR_MathGen_Mess=" Too many iterations in TQLI2"
                   exit
                end if

                iter=iter+1
                g=(d(l+1)-d(l))/(2.0*e(l))
                r=sqrt(g*g+1.0)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                   f=s*e(i)
                   b=c*e(i)
                   r=sqrt(f*f+g*g)
                   e(i+1)=r
                   if (abs(r) <= ep_ss) then
                      d(i+1)=d(i+1)-p
                      e(m)=0.0
                      cycle do_g
                   end if
                   s=f/r
                   c=g/r
                   g=d(i+1)-p
                   r=(d(i)-g)*s+2.0*c*b
                   p=s*r
                   d(i+1)=g+p
                   g=c*r-b

                   !---- omit lines from here ...
                   do k=1,n
                      f=z(k,i+1)
                      z(k,i+1)=s*z(k,i)+c*f
                      z(k,i)=c*z(k,i)-s*f
                   end do

                   !---- ... to here when finding only eigenvalues.
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
                cycle do_g
             end if
             exit
          end do do_g
       end do

       return
    End Subroutine Tqli2

    !!--++
    !!--++ Subroutine Tred1(a,n,d,e)
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for preparing the matrix to find only eigenvalues
    !!--++    Householder reduction of a real symetric nxn matrix A.
    !!--++    On output A is replaced by the orthogonal matrix Q effecting
    !!--++    the transformation. D returns the diagonal elements of the tri-
    !!--++    diagonal matrix and E the off-diagonal elements with E(1)=0.
    !!--++    In tred1 several lines have been deleted and A contains no
    !!--++    useful information on output.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tred1(a,n,d,e)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)
       integer,                       intent(in)     :: n
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)

       !---- Local Variables ----!
       integer :: i, j, k, l
       real(kind=cp)    :: f, g, h, hh, scala

       do i=n,2,-1
          l=i-1
          h=0.0
          scala=0.0
          if (l > 1)then
             do k=1,l
                scala=scala+abs(a(i,k))
             end do
             if (abs(scala) <= ep_ss) then
                e(i)=a(i,l)
             else
                do k=1,l
                   a(i,k)=a(i,k)/scala
                   h=h+a(i,k)**2
                end do
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scala*g
                h=h-f*g
                a(i,l)=f-g
                f=0.0
                do j=1,l
                   g=0.0
                   do k=1,j
                      g=g+a(j,k)*a(i,k)
                   end do
                   do k=j+1,l
                      g=g+a(k,j)*a(i,k)
                   end do
                   e(j)=g/h
                   f=f+e(j)*a(i,j)
                end do
                hh=f/(h+h)
                do j=1,l
                   f=a(i,j)
                   g=e(j)-hh*f
                   e(j)=g
                   do k=1,j
                      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                   end do
                end do
             end if
          else
             e(i)=a(i,l)
          end if
          d(i)=h
       end do

       e(1)=0.0
       do i=1,n
          d(i)=a(i,i)
       end do

       return
    End Subroutine Tred1

    !!--++
    !!--++ Subroutine Tred2(a,n,d,e)
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for preparing the matrix to find the complete set
    !!--++    of eigenvectors.
    !!--++    Householder reduction of a real symetric nxn matrix A.
    !!--++    On output A is replaced by the orthogonal matrix Q effecting
    !!--++    the transformation. D returns the diagonal elements of the tri-
    !!--++    diagonal matrix and E the off-diagonal elements with E(1)=0.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tred2(a,n,d,e)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)
       integer,                       intent(in)     :: n
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)

       !---- Local variables ----!
       integer :: i, j, k, l
       real(kind=cp)    :: f, g, h, hh, scala

       do i=n,2,-1
          l=i-1
          h=0.0
          scala=0.0
          if (l > 1)then
             do k=1,l
                scala=scala+abs(a(i,k))
             end do
             if (abs(scala) <= ep_ss) then
                e(i)=a(i,l)
             else
                do k=1,l
                   a(i,k)=a(i,k)/scala
                   h=h+a(i,k)**2
                end do
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scala*g
                h=h-f*g
                a(i,l)=f-g
                f=0.0
                do j=1,l
                   !---- omit following line if finding only eigenvalues
                   a(j,i)=a(i,j)/h
                   g=0.0
                   do k=1,j
                      g=g+a(j,k)*a(i,k)
                   end do
                   do k=j+1,l
                      g=g+a(k,j)*a(i,k)
                   end do
                   e(j)=g/h
                   f=f+e(j)*a(i,j)
                end do
               hh=f/(h+h)
                do j=1,l
                   f=a(i,j)
                   g=e(j)-hh*f
                   e(j)=g
                   do k=1,j
                      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                   end do
                end do
             end if
          else
             e(i)=a(i,l)
          end if
          d(i)=h
       end do

       !---- omit following line if finding only eigenvalues.
       d(1)=0.0
       e(1)=0.0
       do i=1,n
          !---- delete lines from here ...
          l=i-1
          if (abs(d(i)) > ep_ss)then
             do j=1,l
                g=0.0
                do k=1,l
                   g=g+a(i,k)*a(k,j)
                end do
                do k=1,l
                   a(k,j)=a(k,j)-g*a(k,i)
                end do
             end do
          end if
          !---- ... to here when finding only eigenvalues.
          d(i)=a(i,i)
          !---- also delete lines from here ...
          a(i,i)=1.0
          do j=1,l
             a(i,j)=0.0
             a(j,i)=0.0
          end do
          !---- ... to here when finding only eigenvalues.
       end do

       return
    End Subroutine Tred2

 End Module CFML_Math_General

