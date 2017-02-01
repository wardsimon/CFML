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
  Module CFML_Rational_Arithmetic

    Use CFML_GlobalDeps,       only : cp,dp
    Use CFML_String_Utilities, only : Pack_String
    Use CFML_Math_general,     only : determinant,invert_matrix
    implicit none
    private
    public :: rational  !Type
    public :: rational_simplify, rational_determinant,& !Calculation Procedures
              rational_inv_matrix
    public :: print_rational !transform a rational type to a string like xxx/yyy
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
    public :: nint
    public :: modulo
    public :: mod
    public :: dot_product
    public :: maxval
    public :: minval
    public :: maxloc
    public :: minloc
    public :: matmul

    type :: rational
      integer :: numerator
      integer :: denominator
    end type rational

    logical,           public :: Err_Rational=.false.
    character(len=80), public :: Err_Rational_Mess

    interface assignment (=)
      module procedure assign_rational_int, assign_rational_real
      module procedure assign_int_rational, assign_real_rational
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

    interface nint
      module procedure nint_rational
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
      module procedure rational_maxval_vector
      module procedure rational_maxval_matrix
    end interface

    interface minval
      module procedure rational_minval_vector
      module procedure rational_minval_matrix
    end interface

    interface matmul
      module procedure rational_matmul_matvec
      module procedure rational_matmul_vecmat
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

  contains

    pure recursive function gcd (i, j) result (res)
      integer, intent (in) :: i
      integer, intent (in) :: j
      integer :: res
      if (j == 0) then
        res = i
      else
        res = gcd (j, modulo (i, j))
      end if
    end function gcd

    elemental function rational_simplify (r) result (res)
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
      type (rational), intent (out) :: res  !, volatile
      integer,         intent (in)  :: i
      res = i // 1
    end subroutine assign_rational_int

    elemental subroutine assign_rational_real (res,xr)
      type (rational), intent(out) :: res  !, volatile
      real(kind=cp),   intent (in) :: xr
      integer                     :: maxden,ai,t,si
      real(kind=cp)               :: x,startx !,er1,er2
      integer, dimension(0:1,0:1) :: m

      maxden=99
      m = 0; m(0,0)=1; m(1,1)=1
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
        if( x == real(ai,kind=cp)) exit !division by zero
        x = 1/(x - real(ai,kind=cp))
      end do
      res= si*m(0,0)// m(1,0)

      !er1=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)
      !Second option
      ! ai = (maxden - m(2,2)) / m(1,2)
      ! m(1,1) = m(1,1) * ai + mm(2,1)
      ! m(1,2) = mm(1,2) * ai + mm(2,2)
      ! er2=startx - real(m(1,1),kind=8) / real(m(1,2),kind=8)

    end subroutine assign_rational_real

    elemental subroutine assign_int_rational (i, res)
      type (rational), intent (in)   :: res  !, volatile
      integer,         intent (out)  :: i
      i= nint(real(res%numerator,kind=dp)/real(res%denominator,kind=dp))
    end subroutine assign_int_rational

    elemental subroutine assign_real_rational (x,res)
      type (rational), intent(in)   :: res
      real(kind=cp),   intent (out) :: x
      x=real(res%numerator,kind=dp)/real(res%denominator,kind=dp)
    end subroutine assign_real_rational

    elemental function rational_add (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % denominator + r % denominator * s % numerator // &
          & r % denominator * s % denominator
      res=rational_simplify(res)
    end function rational_add

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

    elemental function rational_multiply (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      res = r % numerator * s % numerator // r % denominator * s % denominator
      res=rational_simplify(res)
    end function rational_multiply

    elemental function rational_divide (r, s) result (res)
      type (rational), intent (in) :: r
      type (rational), intent (in) :: s
      type (rational) :: res
      integer :: denom
      denom = r % denominator * s % numerator
      if(denom /= 0) then
        res = r % numerator * s % denominator // denom
        res=rational_simplify(res)
      else
        res=0//0
      end if
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

    elemental function rational_nint (r) result(res)
      real(kind=cp),  intent (in)  :: r
      type (rational)              :: res
      res = nint(r) // 1
    end function rational_nint

    elemental function nint_rational (r) result(res)
      type (rational),  intent (in)  :: r
      integer                        :: res
      res = nint(real(r%numerator,kind=dp)/real(r%denominator,kind=dp))
    end function nint_rational


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

    function rational_matmul_vecmat (vec,mat) result (vec_out)
      type (rational), dimension(:),   intent (in) :: vec
      type (rational), dimension(:,:), intent (in) :: mat
      type (rational), dimension(size(vec))        :: vec_out

      integer :: n1,n2,n3,i
      n1=size(mat,dim=1); n2=size(mat,dim=2); n3=size(vec)
      if(n1 == n2 .and. n2 == n3) then
        do i=1,n3
          vec_out(i) = rational_simplify (dot_product(vec,mat(:,i)))
        end do
      else
        Err_Rational=.true.
        write(unit=Err_Rational_Mess,fmt="(a,3i4)") "Error in MATMUL: the dimensions of the arguments are non-conformable",n1,n2,n3
      end if
    end function rational_matmul_vecmat

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
      res=-huge(1)//1
      do i=1,n
        if(r(i) > res) res=r(i)
      end do
    end function rational_maxval_vector

    function rational_maxval_matrix (r) result (res)
      type (rational), dimension(:,:), intent (in) :: r
      type (rational) :: res
      integer :: i,j,n1,n2
      n1=size(r,dim=1);  n2=size(r,dim=2)
      res=-huge(1)//1
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
      res=huge(1)//1
      do i=1,n
        if(r(i) < res) res=r(i)
      end do
    end function rational_minval_vector

    function rational_minval_matrix (r) result (res)
      type (rational), dimension(:,:), intent (in) :: r
      type (rational) :: res
      integer :: i,j,n1,n2
      n1=size(r,dim=1);  n2=size(r,dim=2)
      res=huge(1)//1
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
      if(r%denominator /= 0) then
        sr=rational_simplify(r)
      else
        sr=r
      end if
      if(sr%denominator == 1) then
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
      integer:: n1,n2
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

    function rational_maxloc_mat(Mat) result(pos_max)
      type(rational), dimension(:,:), intent(in) :: Mat
      integer,        dimension(2)               :: pos_max
      !Local variables
      integer:: nu1,nl1,nu2,nl2,i,j
      type(rational) :: res
      nu1=ubound(Mat,dim=1); nu2=ubound(Mat,dim=2)
      nl1=lbound(Mat,dim=1); nl2=lbound(Mat,dim=2)

      res=-huge(1)//1
      do j=nl2,nu2
        do i=nl1,nu1
          if(mat(i,j) > res) then
            res=mat(i,j)
            pos_max=[i,j]
          end if
        end do
      end do
    end function rational_maxloc_mat

    function rational_maxloc_vect(vec) result(pos_max)
      type(rational), dimension(:), intent(in) :: vec
      integer                                  :: pos_max
      !Local variables
      integer:: nu,nl,i
      type(rational) :: res
      nu=ubound(vec,dim=1)
      nl=lbound(vec,dim=1)
      res=-huge(1)//1
      do i=nl,nu
        if(vec(i) > res) then
          res=vec(i)
          pos_max=i
        end if
      end do
    end function rational_maxloc_vect

    function rational_minloc_mat(Mat) result(pos_min)
      type(rational), dimension(:,:), intent(in) :: Mat
      integer,        dimension(2)               :: pos_min
      !Local variables
      integer:: nu1,nl1,nu2,nl2,i,j
      type(rational) :: res
      nu1=ubound(Mat,dim=1); nu2=ubound(Mat,dim=2)
      nl1=lbound(Mat,dim=1); nl2=lbound(Mat,dim=2)

      res=huge(1)//1
      do j=nl2,nu2
        do i=nl1,nu1
          if(mat(i,j) < res) then
            res=mat(i,j)
            pos_min=[i,j]
          end if
        end do
      end do

    end function rational_minloc_mat

    function rational_minloc_vect(vec) result(pos_min)
      type(rational), dimension(:), intent(in) :: vec
      integer                                  :: pos_min
      !Local variables
      integer:: nu,nl
      integer:: i
      type(rational) :: res
      nu=ubound(vec,dim=1)
      nl=lbound(vec,dim=1)
      res=huge(1)//1
      do i=nl,nu
        if(vec(i) < res) then
          res=vec(i)
          pos_min=i
        end if
      end do
    end function rational_minloc_vect

  End Module CFML_Rational_Arithmetic