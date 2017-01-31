!!----  Module CFML_Rational_Arithmetic
!!----  Extension of the module: module_rational from
!!----  Rosetta Code: https://rosettacode.org/wiki/Arithmetic/Rational#Fortran
!!----
!!----
!!----
  Module CFML_Rational_Arithmetic

    implicit none
    private
    public :: rational
    public :: rational_simplify
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
    public :: abs
    public :: int
    public :: modulo
    public :: mod
    public :: dot_product
    public :: maxval
    public :: minval
    public :: maxloc
    public :: minloc

    type :: rational
      integer :: numerator
      integer :: denominator
    end type rational

    logical           :: Error_Rational
    character(len=80) :: Error_Rational_Mess

    interface assignment (=)
      module procedure assign_rational_int, assign_rational_real
    end interface

    interface operator (//)
      module procedure make_rational
    end interface

    interface operator (+)
      module procedure rational_add
    end interface

    interface operator (-)
      module procedure rational_minus, rational_subtract
    end interface

    interface operator (*)
      module procedure rational_multiply
    end interface

    interface operator (/)
      module procedure rational_divide
    end interface

    interface operator (<)
      module procedure rational_lt
    end interface

    interface operator (<=)
      module procedure rational_le
    end interface

    interface operator (>)
      module procedure rational_gt
    end interface

    interface operator (>=)
      module procedure rational_ge
    end interface

    interface operator (==)
      module procedure rational_eq
    end interface

    interface operator (/=)
      module procedure rational_ne
    end interface

    interface abs
      module procedure rational_abs
    end interface

    interface int
      module procedure rational_int
    end interface

    interface modulo
      module procedure rational_modulo
    end interface

    interface mod
      module procedure rational_mod
    end interface

    interface dot_product
      module procedure rational_dot_product
    end interface

    interface maxval
      module procedure rational_maxval
    end interface

    interface minval
      module procedure rational_minval
    end interface

   ! interface maxloc
   !   module procedure rational_maxloc
   ! end interface
   !
   ! interface minloc
   !   module procedure rational_minloc
   ! end interface

  contains

    recursive function gcd (i, j) result (res)
      integer, intent (in) :: i
      integer, intent (in) :: j
      integer :: res
      if (j == 0) then
        res = i
      else
        res = gcd (j, modulo (i, j))
      end if
    end function gcd

    function rational_simplify (r) result (res)
      type (rational), intent (in) :: r
      type (rational) :: res
      integer :: g
      g = gcd (r % numerator, r % denominator)
      res = r % numerator / g // r % denominator / g
    end function rational_simplify

    elemental function make_rational (numerator, denominator) result (res)
      integer, intent (in) :: numerator
      integer, intent (in) :: denominator
      type (rational) :: res
      res = rational (numerator, denominator)
    end function make_rational

    subroutine assign_rational_int (res, i)
      type (rational), intent (out), volatile :: res
      integer, intent (in) :: i
      res = i // 1
    end subroutine assign_rational_int

    elemental subroutine assign_rational_real (res,x)
      real,            intent (in)           :: x
      type (rational), intent(out), volatile :: res
      integer :: x_floor
      real :: x_frac
      x_floor = floor (x)
      x_frac = x - x_floor
      if (x_frac == 0) then
        res = x_floor // 1
      else
        res = (x_floor // 1) + (1 // floor (1 / x_frac))
      end if
    end subroutine assign_rational_real

    elemental function rational_add (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % denominator + r % denominator * s % numerator // &
          & r % denominator * s % denominator
    end function rational_add

    elemental function rational_minus (r) result (res)
      type (rational), intent (in) :: r
      type (rational) :: res
      res = - r % numerator // r % denominator
    end function rational_minus

    elemental function rational_subtract (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % denominator - r % denominator * s % numerator // &
          & r % denominator * s % denominator
    end function rational_subtract

    elemental function rational_multiply (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % numerator // r % denominator * s % denominator
    end function rational_multiply

    elemental function rational_divide (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % denominator // r % denominator * s % numerator
    end function rational_divide

    function rational_lt (r, s) result (res)
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

    function rational_le (r, s) result (res)
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

    function rational_gt (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: r_simple
      type (rational) :: s_simple
      logical :: res
      r_simple = rational_simplify (r)
      s_simple = rational_simplify (s)
      res = r_simple % numerator * s_simple % denominator > &
          & s_simple % numerator * r_simple % denominator
    end function rational_gt

    function rational_ge (r, s) result (res)
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

    function rational_eq (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      logical :: res
      res = r % numerator * s % denominator == s % numerator * r % denominator
    end function rational_eq

    function rational_ne (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      logical :: res
      res = r % numerator * s % denominator /= s % numerator * r % denominator
    end function rational_ne

    elemental function rational_abs (r) result (res)
      type (rational), intent (in) :: r
      type (rational) :: res
      res = sign (r % numerator, r % denominator) // r % denominator
    end function rational_abs

    elemental function rational_int (r) result (res)
      type (rational), intent (in) :: r
      integer :: res
      res = r % numerator / r % denominator
    end function rational_int

    elemental function rational_modulo (r) result (res)
      type (rational), intent (in) :: r
      integer :: res
      res = modulo (r % numerator, r % denominator)
    end function rational_modulo

    elemental function rational_mod (r) result (res)
      type (rational), intent (in) :: r
      integer :: res
      res = mod (r % numerator, r % denominator)
    end function rational_mod

    function rational_dot_product (r1,r2) result (res)
      type (rational), dimension(:), intent (in) :: r1,r2
      type (rational) :: res
      integer :: n1,n2,i
      n1=size(r1); n2=size(r2)
      res=0//1
      if(n1 == n2) then
        do i=1,n1
          res = rational_simplify (res + r1(i)*r2(i))
        end do
      else
        Error_Rational=.true.
        write(unit=Error_Rational_Mess,fmt="(a,2i4)") "Error in DOT_PRODUCT: the dimensions of the argument are different",n1,n2
      end if
    end function rational_dot_Product

    function rational_maxval (r) result (res)
      type (rational), dimension(*), intent (in) :: r
      type (rational) :: res
      integer :: i,n
      n=size(r)
      res=-(2**32-1)//1
      do i=1,n
        if(r(i) > res) res=r(i)
      end do
    end function rational_maxval

    function rational_minval (r) result (res)
      type (rational), dimension(*), intent (in) :: r
      type (rational) :: res
      integer :: i,n
      n=size(r)
      res=(2**32-1)//1
      do i=1,n
        if(r(i) < res) res=r(i)
      end do
    end function rational_minval

  End Module CFML_Rational_Arithmetic



program perfect_numbers

  use CFML_Rational_Arithmetic
  implicit none
  character(len=8)   :: key
  integer, parameter :: n_min = 2
  integer, parameter :: n_max = 2 ** 19 - 1
  integer :: n
  integer :: factor
  type (rational) :: sum
  type (rational), dimension(n_max) :: rat_vec

  write(*,"(2(a,i8))") " Writing Perfect Numbers between ",n_min, " and ",n_max
  do n = n_min, n_max
    sum = 1 // n
    factor = 2
    do
      if (factor * factor >= n) then
        exit
      end if
      if (modulo (n, factor) == 0) then
        sum = rational_simplify (sum + (1 // factor) + (factor // n))
      end if
      factor = factor + 1
    end do
    if (sum % numerator == 1 .and. sum % denominator == 1) then
      write (*, '(i0)') n
    end if
    rat_vec(n)=sum
  end do
  !do n = n_min, n_max
  !  write(*,"(i8,a, i12,a,i12)") n," -> ",rat_vec(n)%numerator," / ",rat_vec(n)%denominator
  !end do
  do
    write(*,"(a)", advance="no") " => Enter Option: "
    read(*,"(a)") key
    if(len_trim(key) == 0) exit
  end do

end program perfect_numbers