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
!!---- MODULE: CFML_FFT
!!----   INFO: FFT Calculations Routines
!!----
!!---- HISTORY
!!--..    Update: 02/03/2011
!!--..
!!--..    Multivariate Fast Fourier Transform
!!--..    Fortran 90 (ELF90) Implementation of Singleton's mixed-radix
!!--..    algorithm, RC Singleton, Stanford Research Institute, Sept. 1968.
!!--..    Adapted from fftn.c, translated from Fortran 66 to C by Mark Olesen
!!--..    and John Beale.
!!--..    Fourier transforms can be computed either in place, using assumed
!!--..    size arguments, or by generic function, using assumed shape arguments.
!!--..
!!--..    FFT(array, dim, inv)                 generic transform function
!!--..    COMPLEX(fftkind), DIMENSION(:,...,:), INTENT(IN)           :: array
!!--..    INTEGER,          DIMENSION(:),       INTENT(IN),  OPTIONAL:: dim
!!--..    LOGICAL,                              INTENT(IN),  OPTIONAL:: inv
!!--..
!!--..    Formal Parameters:
!!--..
!!--..    array    The complex array to be transformed. array can be of arbitrary
!!--..             rank (i.e. up to seven).
!!--..
!!--..    shape    With subroutine fftn, the shape of the array to be transformed
!!--..             has to be passed separately, since fftradix - the internal trans-
!!--..             formation routine - will treat array always as one dimensional.
!!--..             The product of elements in shape must be the number of
!!--..             elements in array.
!!--..             Although passing array with assumed shape would have been nicer,
!!--..             I prefered assumed size in order to prevent the compiler from
!!--..             using a copy-in-copy-out mechanism. That would generally be
!!--..             necessary with fftn passing array to fftradix and with fftn
!!--..             being prepared for accepting non consecutive array sections.
!!--..             Using assumed size, it's up to the user to pass an array argu-
!!--..             ment, that can be addressed as continous one dimensional array
!!--..             without copying. Otherwise, transformation will not really be
!!--..             performed in place.
!!--..             On the other hand, since the rank of array and the size of
!!--..             shape needn't match, fftn is appropriate for handling more than
!!--..             seven dimensions.
!!--..             As far as function fft is concerned all this doesn't matter,
!!--..             because the argument will be copied anyway. Thus no extra
!!--..             shape argument is needed for fft.
!!--..
!!--..    Optional Parameters:
!!--..
!!--..    dim      One dimensional integer array, containing the dimensions to be
!!--..             transformed. Default is (/1,...,N/) with N being the rank of
!!--..             array, i.e. complete transform. dim can restrict transformation
!!--..             to a subset of available dimensions. Its size must not exceed the
!!--..             rank of array or the size of shape respectivly.
!!--..
!!--..    inv      If .true., inverse transformation will be performed. Default is
!!--..             .false., i.e. forward transformation.
!!--..
!!--..    stat     If present, a system dependent nonzero status value will be
!!--..             returned in stat, if allocation of temporary storage failed.
!!--..             For functions, the integer variable status is used.
!!--..
!!--..    Scaling:
!!--..             Transformation results will always be scaled by the square
!!--..             root of the product of sizes of each dimension in dim.
!!--..             (See examples below)
!!--..
!!--..    Examples:
!!--..             Let A be a L*M*N three dimensional complex array. Then
!!--..             result = fft(A) will produce a three dimensional transform,
!!--..             scaled by sqrt(L*M*N), while call fftn(A, SHAPE(A)) will do
!!--..             the same in place.
!!--..
!!--..             result = fft(A, dim=(/1,3/)) will transform with respect to
!!--..             the first and the third dimension, scaled by sqrt(L*N).
!!--..
!!--..             result = fft(fft(A), inv=.true.) should (approximately)
!!--..             reproduce A. With B having the same shape as A
!!--..
!!--..             result = fft(fft(A) * CONJG(fft(B)), inv=.true.) will
!!--..             correlate A and B.
!!--..
!!--..     Remarks:
!!--..
!!--..             Following changes have been introduced with respect to fftn.c:
!!--..             - complex arguments and results are of type complex, rather
!!--..               than real an imaginary part separately.
!!--..             - increment parameter (magnitude of isign) has been dropped,
!!--..               inc is always one, direction of transform is given by inv.
!!--..             - maxf and maxp have been dropped. The amount of temporary
!!--..               storage needed is determined by the fftradix routine.
!!--..               Both fftn and fft can handle any size of array. (Maybe they
!!--..               take a lot of time and memory, but they will do it)
!!--..
!!--..             Redesigning fftradix in a way, that it handles assumed shape
!!--..             arrays would have been desirable. However, I found it rather
!!--..             hard to do this in an efficient way. Problems were:
!!--..             - to prevent stride multiplications when indexing arrays. At
!!--..               least our compiler was not clever enough to discover that
!!--..               in fact additions would do the job as well. On the other
!!--..               hand, I haven't been clever enough to find an implementation
!!--..               using array operations.
!!--..             - fftradix is rather large and different versions would be
!!--..               necessaray for each possible rank of array.
!!--..             Consequently, in place transformation still needs the
!!--..             argument stored in a consecutive bunch of memory and can't be
!!--..             performed on array sections like A(100:199:-3, 50:1020).
!!--..             Calling fftn with such sections will most probably imply
!!--..             copy-in-copy-out. However, the function fft works with
!!--..             everything it gets and should be convenient to use.
!!--..
!!--..             To enable this module to be used with ELF90 it appears to be
!!--..             necessary to allocate a 1-D work array into which the
!!--..             multi-dimensional array is copied, and then to copy the
!!--..             results back from the 1-D array to the multi-dimensional
!!--..             array ft.
!!--..
!!--..             Unfortunately, ELF90 will not allow a function to return more
!!--..             than one output variable.   The variable `stat' has been
!!--..             dropped from the function arguments. Users should examine the
!!--..             value of the variable `StatusF' instead. This is a PUBLIC
!!--..             variable declared in this module.
!!--..
!!--..             Michael Steffens, 09.12.96,
!!--..             <Michael.Steffens@mbox.muk.uni-hannover.de>
!!--..
!!--..             ELF90-compatible version by Alan Miller, 29 April 1997
!!--..             Alan.Miller @ mel.dms.csiro.au
!!--..
!!---- DEPENDENCIES
!!----    CFML_GlobalDeps, only: cp
!!----
!!---- VARIABLES
!!--++    FFTKIND                  [Private]
!!--++    COS72                    [Private]
!!--++    PI                       [Private]
!!--++    SIN60                    [Private]
!!--++    SIN72                    [Private]
!!--++    STATUSF                  [Private]
!!--..
!!----    POINTS_INTERVAL_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       CONVOL
!!----       CONVOL_PEAKS
!!----       F_FFT
!!----       FFT
!!--++       FFT1D                 [Overloaded]
!!--++       FFT2D                 [Overloaded]
!!--++       FFT3D                 [Overloaded]
!!--++       FFT4D                 [Overloaded]
!!--++       FFT5D                 [Overloaded]
!!--++       FFT6D                 [Overloaded]
!!--++       FFT7D                 [Overloaded]
!!----
!!----    Subroutines:
!!--++       FFTN                  [Private]
!!--++       FFTRADIX
!!--..       FACTORIZE
!!--..       PERMUTE
!!--..       TRANSFORM
!!----       HFFT
!!----       SFFT
!!----
!!
Module CFML_FFT
   !---- Use Modules ----!
   use CFML_GlobalDeps, only : cp

   !---- Local Variables ----!
   implicit none

   private

   !---- List of public variables ----!

   !---- List of public functions ----!
   public :: Convol, Convol_Peaks, F_FFT

   !---- List of public overloaded procedures: functions ----!
   public :: FFT

   !---- List of public subroutines ----!
   public :: SFFT, HFFT

   !---- List of public overloaded procedures: subroutines ----!

   !---- List of private functions ----!
   private :: Fft1D, Fft2D, Fft3D, Fft4D, Fft5D, Fft6D, FFt7D

   !---- List of private subroutines ----!
   private :: FFTN

   !---- Definitions ----!

   !!--++
   !!--++ FFTKIND
   !!--++    integer, private, parameter:: Fftkind = Kind(0.0)
   !!--++
   !!--++    (PRIVATE)
   !!--++    Default Precicion for FFT Variables
   !!--++
   !!--++ Update: February - 2005
   !!
   integer, private, parameter:: Fftkind = Kind(0.0) !--- Adjust Here For Other Precisions

   !!--++
   !!--++ COS72
   !!--++    Real(Fftkind), Parameter:: Cos72
   !!--++
   !!--++    (PRIVATE)
   !!--++    Cos72 = 0.30901699437494742_Fftkind
   !!--++
   !!--++ Update: February - 2005
   !!
   Real(Fftkind), Parameter:: Cos72 = 0.30901699437494742_Fftkind

   !!--++
   !!--++ PI
   !!--++    Real(Fftkind), Parameter:: Pi
   !!--++
   !!--++    (PRIVATE)
   !!--++    Pi    = 3.14159265358979323_Fftkind
   !!--++
   !!--++ Update: February - 2005
   !!
   Real(Fftkind), Parameter:: Pi    = 3.14159265358979323_Fftkind

   !!--++
   !!--++ SIN60
   !!--++    Real(Fftkind), Parameter:: Sin60
   !!--++
   !!--++    (PRIVATE)
   !!--++    Sin60 = 0.86602540378443865_Fftkind
   !!--++
   !!--++ Update: February - 2005
   !!
   Real(Fftkind), Parameter:: Sin60 = 0.86602540378443865_Fftkind

   !!--++
   !!--++ SIN72
   !!--++    Real(Fftkind), Parameter:: Sin72
   !!--++
   !!--++    (PRIVATE)
   !!--++    Sin72 = 0.95105651629515357_Fftkind
   !!--++
   !!--++ Update: February - 2005
   !!
   Real(Fftkind), Parameter:: Sin72 = 0.95105651629515357_Fftkind

   !!--++
   !!--++ STATUSF
   !!--++    integer, private, Save     :: StatusF
   !!--++
   !!--++    Information on FFT Routines
   !!--++
   !!--++ Update: February - 2005
   !!
   integer, private, Save     :: StatusF    !--- Shifted To Here As Elf90 Does Not Allow
                                            !    Arguments To Be Intent(Out)
   !!----
   !!---- TYPE, public :: Points_Interval_Type
   !!--..
   !!---- Type, public :: Points_Interval_Type
   !!----   integer       :: Np                 ! Number of Points
   !!----   real(kind=cp) :: Low                ! Lower range value
   !!----   real(kind=cp) :: High               ! Upper range value
   !!---- End Type Points_Interval_Type
   !!----
   !!---- Update: April 2005
   !!----
   Type, public :: Points_Interval_Type
     integer       :: Np
     real(kind=cp) :: Low
     real(kind=cp) :: High
   End Type Points_Interval_Type

   !---- Interfaces - Overlapp ----!
   Interface Fft
      Module Procedure Fft1D
      Module Procedure Fft2D
      Module Procedure Fft3D
      Module Procedure Fft4D
      Module Procedure Fft5D
      Module Procedure Fft6D
      Module Procedure Fft7D
   End Interface

 Contains

    !---- Functions ----!

    !!----
    !!---- Pure Function Convol(F,Pf,G,Pg,Interval)  Result(Conv)
    !!----    real(kind=cp),dimension(:),          intent(in) :: Pf
    !!----    real(kind=cp),dimension(:),          intent(in) :: Pg
    !!----    type(Points_Interval_Type), intent(in) :: Interval
    !!----    real(kind=cp), dimension(interval%np)           :: conv
    !!----    Interface F_Function
    !!----       Pure Function f(x,parf)  result (vf)
    !!----          real(kind=cp),              intent(in) :: x
    !!----          real(kind=cp), dimension(:),intent(in) :: parf
    !!----          real(kind=cp)                          :: vf
    !!----       End Function f
    !!----    End Interface F_Function
    !!----    Interface G_Function
    !!----       Pure function g(x,parg)  result (vg)
    !!----          real(kind=cp), intent(in)              :: x
    !!----          real(kind=cp), dimension(:),intent(in) :: parg
    !!----       End Function G
    !!----    End Interface G_Function
    !!----
    !!--..
    !!--<<   Convolution of the user-provided centred (x=0) peak functions
    !!----   f and g. The characteristic parameters of the functions f and
    !!----   g are provided in vectors pf and pg.
    !!----   The intent-in Points_Interval_Type variable "Interval" gives
    !!----   the number of points and the limits of the interval
    !!----         Number of points:  Interval%np
    !!----     Range of calculation: [ Interval%low, Interval%high ]
    !!----                    step : (Interval%high-Interval%low)/(Interval%np-1)
    !!----   The convolution function is normalized to unit area .
    !!----
    !!----   Example of use:
    !!----      h = convol(Pseudo_Voigt,P_PV, hat, P_hat, my_interval)
    !!----
    !!----   generates my_interval%np values  h(i), i=1,my_interval%np corresponding
    !!-->>   to the convolution of a pseudo-Voigt function with a hat function
    !!----
    !!---- Update: April - 2005
    !!
    Pure Function Convol(F,Pf,G,Pg,Interval) Result(Conv)
       !---- Arguments ----!
       real(kind=cp),dimension(:),          intent(in) :: pf
       real(kind=cp),dimension(:),          intent(in) :: pg
       type(Points_Interval_Type),          intent(in) :: interval
       real(kind=cp), dimension(interval%np)           :: conv

       Interface F_Function
          Pure function f(x,parf)  result (vf)
             use CFML_GlobalDeps, only: cp
             real(kind=cp),              intent(in) :: x
             real(kind=cp), dimension(:),intent(in) :: parf
             real(kind=cp)                          :: vf
          End Function F
       End Interface F_Function

       Interface G_Function
          Pure Function G(X,Parg)  Result (Vg)
             use CFML_GlobalDeps, only: cp
             real(kind=cp),              intent(in) :: x
             real(kind=cp), dimension(:),intent(in) :: parg
             real(kind=cp)                          :: vg
          End Function G
       End Interface G_Function

       !---- Local variables ----!
       integer                         :: i, n, nd2
       real(kind=cp)                   :: step,xvar, value, area
       complex, dimension(interval%np) :: fx,gx,convo

       n=interval%np-1
       nd2=interval%np/2
       step = (interval%high - interval%low)/real(n)
       conv = 0.0

       do i=1,nd2
          xvar=interval%low+real(i-1)*step
          value =  f(xvar,pf)
          fx(nd2+i) =  cmplx(value)
          value =  g(xvar,pg)
          gx(i) =  cmplx(value)
       end do

       do i=nd2+1,interval%np
          xvar=interval%low+real(i-1)*step
          value =  f(xvar,pf)
          fx(i-nd2) =  cmplx(value)
          value =  g(xvar,pg)
          gx(i) =  cmplx(value)
       end do

       call Sfft(fx)
       call Sfft(gx)

       convo = fx * gx

       call Sfft(convo,"INV")

       !> correct for a displacement of 1 step
       !> To recalculate the array for the same input points
       !> one has to interpolate between succesive values

       conv(1) = real(convo(1))
       do i=2,interval%np
          conv(i) = 0.5 * ( real(convo(i-1)) + real(convo(i)) )
       end do

       area=sum(conv)*step
       conv=conv/area

       return
    End Function Convol

    !!----
    !!---- Pure Function Convol_Peaks(F,Pf,G,Pg,Wd,Np)  Result(Conv)
    !!----    real(kind=cp),dimension(:),          intent(in) :: pf !Parameters of the function f (starting with FWHM)
    !!----    real(kind=cp),dimension(:),          intent(in) :: pg !Parameters of the function g (starting with FWHM)
    !!----    real(kind=cp),                       intent(in) :: wd !Number of times a FWHM of the f-function to calculate range
    !!----    integer,                             intent(in) :: np !Number of points (it is increased internally up to the closest power of 2)
    !!----    Interface F_Function
    !!----       Pure function f(x,parf)  result (vf)
    !!----          real(kind=cp),              intent(in) :: x
    !!----          real(kind=cp), dimension(:),intent(in) :: parf
    !!----          real(kind=cp)                          :: vf
    !!----       End Function F
    !!----    End Interface F_Function
    !!----    Interface G_Function
    !!----       Pure function g(x,parg)  result (vg)
    !!----          real(kind=cp), intent(in)              :: x
    !!----          real(kind=cp), dimension(:),intent(in) :: parg
    !!----       End Function G
    !!----    End Interface G_Function
    !!----
    !!----   Convolution of the user-provided centred (x=0) peak functions
    !!----   f and g. The characteristic parameters of the functions f and
    !!----   g are provided in vectors pf and pg. The first component should
    !!----   be the value of the parameter related to the FWHM.
    !!----   The parameter wd controls the range of the calculation in terms
    !!----   of the FWHM of peaks. The definition interval [a,b] of the peaks
    !!----   is calculated as: a=-b, with b=wd*FWHM=wd*pf(1).
    !!----   The number of points to calculate the convolution is "np".
    !!----   Internally, the actual number of points "n". Corresponding to
    !!----   and increased value of np up to the nearest higher power of 2.
    !!----   The convolution function is normalized to unit area .
    !!----   The internal step is:  step=(b-a)/(n-1)
    !!----   Example of use:
    !!----      h = convol_peaks (Pseudo_Voigt,P_PV, hat, P_hat, 15.0, 150)
    !!----   generates 150 values  h(i), i=1,150 corresponding to the convolution
    !!----   of a pseudo-Voigt function with a hat function
    !!----
    !!---- Update: April - 2005
    !!
    Pure Function Convol_Peaks(F,Pf,G,Pg,Wd,Np) Result(Conv)
       !---- Arguments ----!
       real(kind=cp),dimension(:),          intent(in) :: pf !Parameters of the function f (starting with FWHM)
       real(kind=cp),dimension(:),          intent(in) :: pg !Parameters of the function g (starting with FWHM)
       real(kind=cp),                       intent(in) :: wd !Number of times a FWHM of the f-function to calculate range
       integer,                             intent(in) :: np !Number of points (it is increased internally up to the closest power of 2)
       real(kind=cp), dimension(np)                    :: conv

       Interface F_Function
          Pure Function F(X,Parf)  Result (Vf)
             use CFML_GlobalDeps, only: cp
             real(kind=cp),              intent(in) :: x
             real(kind=cp), dimension(:),intent(in) :: parf
             real(kind=cp)                          :: vf
          End Function F
       End Interface F_Function

       Interface G_Function
          Pure Function G(X,Parg) Result (Vg)
             use CFML_GlobalDeps, only: cp
             real(kind=cp),              intent(in) :: x
             real(kind=cp), dimension(:),intent(in) :: parg
             real(kind=cp)                          :: vg
          End Function G
       End Interface G_Function

       !---- Local variables ----!
       integer                                 :: i,j, n, nd2, m
       real(kind=cp)                           :: step, xvar, value, area, a,b, nstep
       logical                                 :: is_powerof2
       complex, dimension(:),allocatable       :: fx,gx,convo
       real(kind=cp), dimension(:),allocatable :: xv

       ! m will be the effective number of points used in FFT
       do i=1,20
          m=2**i
          if (m >= np) exit
       end do
       is_powerof2=.false.
       if( m == np) is_powerof2=.true.

       n=m-1
       nd2=m/2
       b=wd*pf(1)
       a=-b
       step = (b-a)/real(n)
       conv = 0.0

       if(allocated(fx)) deallocate(fx)
       allocate(fx(m))

       if(allocated(gx)) deallocate(gx)
       allocate(gx(m))

       if(allocated(convo)) deallocate(convo)
       allocate(convo(m))

       if(allocated(xv)) deallocate(xv)
       allocate(xv(m))

       do i=1,nd2
          xvar  = a+real(i-1)*step
          xv(i) =  xvar
          value =  f(xvar,pf)
          fx(nd2+i) =  cmplx(value)
          value =  g(xvar,pg)
          gx(i) =  cmplx(value)
       end do

       do i=nd2+1,m
          xvar  = a+real(i-1)*step
          xv(i) = xvar
          value =  f(xvar,pf)
          fx(i-nd2) =  cmplx(value)
          value =  g(xvar,pg)
          gx(i) =  cmplx(value)
       end do

       call Sfft(fx)     !Forward FFT
       call Sfft(gx)

       convo = fx * gx

       call Sfft(convo,"INV")  !Backward FFT

       !correct for a displacement of 1 step
       !To recalculate the array for the same input points
       !one has to interpolate between succesive values

       fx(1) = convo(1)
       do i=2,m
          fx(i) = 0.5 * ( convo(i-1) + convo(i) )
       end do
       area=sum(real(fx))*step
       fx=fx/area

       ! Calculate an interpolated array for the number of points supplied by the user
       if (.not. is_powerof2) then
          nstep = (b-a)/real(np-1) !New step
          n=1
          conv(np) = real(fx(m))

          do i=1,np-1
             xvar=a+real(i-1)*nstep
             do j=n,m-1
                if (xv(j) >= xvar) then
                   n=j
                   exit
                end if
             end do
             conv(i) = real(fx(n))+(xvar-xv(n))*(real(fx(n+1))- real(fx(n)))/step
          end do

       else
          conv=real(fx)
       end if

       return
    End Function Convol_Peaks

    !!----
    !!---- Pure Function F_FFT(Array,Mode) result(fft)
    !!----    complex, dimension(:),     intent (in) :: Array  !  In -> real array containing real parts of transform
    !!----    character(len=*),optional, intent (in) :: Mode   !  In -> type="INV"    : backward transform
    !!----                                                              Absent or whatever : forward transform
    !!----
    !!----    This Function is a slight modification of a complex split
    !!----    radix FFT routine presented by C.S. Burrus. There is no control
    !!----    of the error consisting in giving a dimension that is not a power
    !!----    of two. It is the responsibility of the user to provide a complex
    !!----    array of dimension equal to a power of 2.
    !!----    The function is similar to subroutine SFFT and it is useful only when
    !!----    one is interested in conserving the original array.
    !!----    Example of use:
    !!----
    !!----    FX = F_fft(X)
    !!----     Y = F_fft(FY,"INV")
    !!----
    !!---- Update: February - 2005
    !!
    Pure Function F_FFT(Array,Mode ) result (fft)
       !---- Arguments ----!
       complex, dimension(:),      intent(in) :: Array
       character(len=*), optional, intent(in) :: Mode
       complex, dimension(size(Array))           :: fft

       !---- Local variables ----!
       integer                               :: i, j, k, n, m, n1, n2, n4, is, id, i0, i1, i2, i3
       real(kind=cp)                         :: r1, r2, s1, s2, s3, xt
       real(kind=cp)                         :: e, a, a3, cc1, ss1, cc3, ss3
       real(kind=cp), parameter              :: twopi=6.2831853071795864769
       real(kind=cp), dimension(size(Array)) :: x
       real(kind=cp), dimension(size(Array)) :: y

       n=size(Array)
       m=0
       do i=1,20
          if (2**i /= n) cycle
          m=i
          exit
       end do

       do i=1,n
          x(i)=real(Array(i))
          y(i)=aimag(Array(i))
       end do
       if(present(Mode)) then
        if (Mode == "INV") y=-y
       end if

       n2 = 2 * n
       do k = 1, m-1
          n2 = n2 / 2
          n4 = n2 / 4
          e = twopi / n2
          a = 0.0
          do j = 1, n4
             a3 = 3 * a
             cc1 = cos( a )
             ss1 = sin( a )
             cc3 = cos( a3 )
             ss3 = sin( a3 )
             a = j * e
             is = j
             id = 2 * n2
             do
                do i0 = is, n-1, id
                   i1 = i0 + n4
                   i2 = i1 + n4
                   i3 = i2 + n4
                   r1 = x(i0) - x(i2)
                   x(i0) = x(i0) + x(i2)
                   r2 = x(i1) - x(i3)
                   x(i1) = x(i1) + x(i3)
                   s1 = y(i0) - y(i2)
                   y(i0) = y(i0) + y(i2)
                   s2 = y(i1) - y(i3)
                   y(i1) = y(i1) + y(i3)
                   s3 = r1 - s2
                   r1 = r1 + s2
                   s2 = r2 - s1
                   r2 = r2 + s1
                   x(i2) = r1 * cc1 - s2 * ss1
                   y(i2) = - s2 * cc1 - r1 * ss1
                   x(i3) = s3 * cc3 + r2 * ss3
                   y(i3) = r2 * cc3 - s3 * ss3
                end do
                is = 2 * id - n2 + j
                id = 4 * id
                if (is >= n) exit
             end do
          end do
       end do

       !- LAST STAGE, LENGTH-2 BUTTERFLY
       is = 1
       id = 4
       do
          do i0 = is, n, id
             i1 = i0 + 1
             r1 = x(i0)
             x(i0) = r1 + x(i1)
             x(i1) = r1 - x(i1)
             r1 = y(i0)
             y(i0) = r1 + y(i1)
             y(i1) = r1 - y(i1)
          end do
          is = 2 * id - 1
          id = 4 * id
          if (is >= n) exit
       end do

       ! BIT REVERSE COUNTER
       j = 1
       n1 = n - 1
       do i = 1, n1
          if (i < j) then
             xt = x(j)
             x(j) = x(i)
             x(i) = xt
             xt = y(j)
             y(j) = y(i)
             y(i) = xt
          end if
          k = n / 2

          do
             if (k >=j) exit
             j = j - k
             k = k / 2
          end do
          j = j + k
       end do

       if(present(Mode)) then
        if (Mode == "INV") then
          x=x/n
          y=-y/n
        end if
       end if

       do i=1,n
          fft(i)=cmplx(x(i),y(i))
       end do

       return
    End Function F_FFT

    !!----
    !!---- Function Fft(Array, Dim, Inv) Result(Ft)
    !!----    complex(fftkind), dimension(:), intent(in)            :: array  !  In -> Complex array
    !!----    integer,          dimension(:), intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!----                                                                             to be transformed
    !!----    logical,                        intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!----                                                                             Default is .false., i.e. forward transformation.
    !!----    Calculation of FFT from 1 to up 7 dimensions
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Fft1D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:), intent(in)            :: array  !  In -> Complex array
    !!--++    integer,          dimension(:), intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                             to be transformed
    !!--++    logical,                        intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                             Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT one dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft1D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:), intent(in)           :: array
       integer,          dimension(:), intent(in),  optional:: dim
       logical,                        intent(in),  optional:: inv

       !--- function result
       complex(fftkind), dimension(size(array, 1)):: ft
       ft = array
       call fftn(ft, shape(array), dim, inv = inv, stat = StatusF)

       return
    End Function Fft1D

    !!--++
    !!--++ Function Fft2D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:), intent(in)            :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),   intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                               to be transformed
    !!--++    logical,                          intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                               Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT two dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft2D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:), intent(in)           :: array
       integer,          dimension(:),   intent(in),  optional:: dim
       logical,                          intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension(size(array, 1), size(array, 2)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft2D

    !!--++
    !!--++ Function Fft3D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:), intent(in)            :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),     intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                 to be transformed
    !!--++    logical,                            intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                 Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT three dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft3D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:), intent(in)           :: array
       integer,          dimension(:),     intent(in),  optional:: dim
       logical,                            intent(in),  optional:: inv
       !--- function result
       complex(fftkind), &
            dimension(size(array, 1), size(array, 2), size(array, 3)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft3D

    !!--++
    !!--++ Function Fft4D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:,:), intent(in)            :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),       intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                   to be transformed
    !!--++    logical,                              intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                   Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT four dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft4D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:,:), intent(in)           :: array
       integer,          dimension(:),       intent(in),  optional:: dim
       logical,                              intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3), &
                             size(array, 4) /))
       if (allocated(work)) deallocate(work)

       return
    End Function Fft4D

    !!--++
    !!--++ Function Fft5D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:,:,:), intent(in)          :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),       intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                   to be transformed
    !!--++    logical,                              intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                   Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT five dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft5D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:,:,:), intent(in)           :: array
       integer,          dimension(:),         intent(in),  optional:: dim
       logical,                                intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3), &
                             size(array, 4), size(array, 5) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft5D

    !!--++
    !!--++ Function Fft6D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:,:,:,:), intent(in)        :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),       intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                   to be transformed
    !!--++    logical,                              intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                   Default is .false., i.e. forward transformation.
    !!--++    (Overloaded)
    !!--++    FFT six dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft6D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:,:,:,:), intent(in)           :: array
       integer,          dimension(:),           intent(in),  optional:: dim
       logical,                                  intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5), size(array, 6)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3), &
                             size(array, 4), size(array, 5), size(array, 6) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft6D

    !!--++
    !!--++ Function Fft7D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:,:,:,:,:), intent(in)      :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),       intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                   to be transformed
    !!--++    logical,                              intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                   Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT seven dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft7D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:,:,:,:,:), intent(in)           :: array
       integer,          dimension(:),             intent(in),  optional:: dim
       logical,                                    intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5), size(array, 6), size(array, 7)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3), &
                             size(array, 4), size(array, 5), size(array, 6), &
                             size(array, 7) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft7D

    !---- Subroutines ----!

    !!--++
    !!--++ Subroutine Fftn(Array, Shape, Dim, Inv, Stat)
    !!--++    complex(fftkind), dimension(:), intent(in out)       :: array
    !!--++    integer,          dimension(:), intent(in)           :: shape
    !!--++    integer,          dimension(:), intent(in), optional :: dim
    !!--++    logical,                        intent(in), optional :: inv
    !!--++    integer,                        intent(in), optional :: stat
    !!--++
    !!--++    (PRIVATE)
    !!--++    General routine for FFT calculations
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Fftn(Array, Shape, Dim, Inv, Stat)
       !--- formal parameters
       complex(fftkind), dimension(:), intent(in out)       :: array
       integer,          dimension(:), intent(in)           :: shape
       integer,          dimension(:), intent(in),  optional:: dim
       logical,                        intent(in),  optional:: inv
       integer,                        intent(out), optional:: stat

       !--- local arrays
       integer, dimension(size(shape)):: d
       !--- local scalars
       logical      :: inverse
       integer      :: i, ndim, ntotal
       real(fftkind):: scal

       !--- optional parameter settings
       if (present(inv)) then
          inverse = inv
       else
          inverse = .false.
       end if
       if (present(dim)) then
          ndim = min(size(dim), size(d))
          d(1:ndim) = dim(1:ndim)
       else
          ndim = size(d)
          d = (/(i, i = 1, size(d))/)
       end if
       ntotal = product(shape)
       scal = sqrt(1.0_fftkind / product(shape(d(1:ndim))))
       array(1:ntotal) = array(1:ntotal) * scal
       do i = 1, ndim
          call fftradix(array, ntotal, shape(d(i)), product(shape(1:d(i))), &
               inverse, stat)
          if (present(stat)) then
             if (stat /=0) return
          end if
       end do

       return
    End Subroutine Fftn

    !!--++
    !!--++ Subroutine Fftradix(Array, Ntotal, Npass, Nspan, Inv, Stat)
    !!--++    complex(fftkind), dimension(:), intent(in out)       :: array
    !!--++    integer,                        intent(in)           :: ntotal
    !!--++    integer,                        intent(in)           :: npass
    !!--++    integer,                        intent(in)           :: nspan
    !!--++    logical,                        intent(in)           :: inv
    !!--++    integer,                        intent(in), optional :: stat
    !!--++
    !!--++    (PRIVATE)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Fftradix(Array, Ntotal, Npass, Nspan, Inv, Stat)
       !--- formal parameters
       INTEGER,                        INTENT(IN)           :: ntotal, npass, nspan
       COMPLEX(fftkind), DIMENSION(:), INTENT(IN OUT)       :: array
       LOGICAL,                        INTENT(IN)           :: inv
       INTEGER,                        INTENT(OUT), OPTIONAL:: stat

       !--- local arrays
       INTEGER,          DIMENSION(BIT_SIZE(0))     :: factor
       COMPLEX(fftkind), DIMENSION(:), ALLOCATABLE  :: ctmp
       REAL(fftkind),    DIMENSION(:), ALLOCATABLE  :: sine, cosine
       INTEGER,          DIMENSION(:), ALLOCATABLE  :: perm

       !--- local scalars
       INTEGER         :: ii, kspan, ispan
       INTEGER         :: j, jc, jf, jj, k, k1, k2, k3, k4, kk, kt, nn, ns, nt
       INTEGER         :: maxfactor, nfactor, nperm
       REAL(fftkind)   :: s60, c72, s72, pi2
       REAL(fftkind)   :: radf
       REAL(fftkind)   :: c1, c2, c3, cd, ak
       REAL(fftkind)   :: s1, s2, s3, sd
       COMPLEX(fftkind):: cc, cj, ck, cjp, cjm, ckp, ckm

       IF (npass <= 1) RETURN
       c72 = cos72
       IF (inv) THEN
          s72 = sin72
          s60 = sin60
          pi2 = pi
       ELSE
          s72 = -sin72
          s60 = -sin60
          pi2 = -pi
       END IF
       nt = ntotal
       ns = nspan
       kspan = ns
       nn = nt - 1
       jc = ns / npass
       radf = pi2 * jc
       pi2 = pi2 * 2.0_fftkind !-- use 2 PI from here on
       CALL factorize()
       maxfactor = MAXVAL(factor (:nfactor))
       IF (nfactor - ISHFT(kt, 1) > 0) THEN
          nperm = MAX(nfactor + 1, PRODUCT(factor(kt+1: nfactor-kt)) - 1)
       ELSE
          nperm = nfactor + 1
       END IF
       IF (PRESENT(stat)) THEN
          ALLOCATE(ctmp(maxfactor), sine(maxfactor), cosine(maxfactor), STAT=stat)
          IF (stat /= 0) RETURN
          CALL transform()
          DEALLOCATE(sine, cosine, STAT=stat)
          IF (stat /= 0) RETURN
          ALLOCATE(perm(nperm), STAT=stat)
          IF (stat /= 0) RETURN
          CALL permute()
          DEALLOCATE(perm, ctmp, STAT=stat)
          IF (stat /= 0) RETURN
       ELSE
          ALLOCATE(ctmp(maxfactor), sine(maxfactor), cosine(maxfactor))
          CALL transform()
          DEALLOCATE(sine, cosine)
          ALLOCATE(perm(nperm))
          CALL permute()
          DEALLOCATE(perm, ctmp)
       END IF
       RETURN

    Contains

       !!--++
       !!--++ Subroutine factorize()
       !!--++
       !!--++    (PRIVATE)
       !!--++
       !!--++ Update: February - 2005
       !!
       SUBROUTINE factorize()

          nfactor = 0
          k = npass
          DO WHILE (MOD(k, 16) == 0)
             nfactor = nfactor + 1
             factor (nfactor) = 4
             k = k / 16
          END DO

          j = 3
          jj = 9
          DO
             DO WHILE (MOD(k, jj) == 0)
                nfactor=nfactor + 1
                factor (nfactor) = j
                k = k / jj
             END DO
             j = j + 2
             jj = j * j
             IF (jj > k) EXIT
          END DO

          IF (k <= 4) THEN
             kt = nfactor
             factor (nfactor + 1) = k
             IF (k /= 1) nfactor = nfactor + 1
          ELSE
             IF (k - ISHFT(k / 4, 2) == 0) THEN
                nfactor = nfactor + 1
                factor (nfactor) = 2
                k = k / 4
             END IF
             kt = nfactor
             j = 2
             DO
                IF (MOD(k, j) == 0) THEN
                   nfactor = nfactor + 1
                   factor (nfactor) = j
                   k = k / j
                END IF
                j = ISHFT((j + 1)/2, 1) + 1
                IF (j > k) EXIT
             END DO
          END IF

          IF (kt > 0) THEN
             j = kt
             DO
                nfactor = nfactor + 1
                factor (nfactor) = factor (j)
                j = j - 1
                IF (j==0) EXIT
             END DO
          END IF

          RETURN
       END SUBROUTINE factorize

       !!--++
       !!--++ Subroutine transform()
       !!--++
       !!--++    compute fourier transform
       !!--++    (PRIVATE)
       !!--++
       !!--++ Update: February - 2005
       !!
       SUBROUTINE transform()

          ii = 0
          jf = 0
          DO
             sd = radf / kspan
             cd = SIN(sd)
             cd = 2.0_fftkind * cd * cd
             sd = SIN(sd + sd)
             kk = 1
             ii = ii + 1
             SELECT CASE (factor (ii))
                CASE (2)
                   !-- transform for factor of 2 (including rotation factor)
                   kspan = kspan / 2
                   k1 = kspan + 2
                   DO
                      DO
                         k2 = kk + kspan
                         ck = array(k2)
                         array(k2) = array(kk) - ck
                         array(kk) = array(kk) + ck
                         kk = k2 + kspan
                         IF (kk > nn) EXIT
                      END DO
                      kk = kk - nn
                      IF (kk > jc) EXIT
                   END DO

                   IF (kk > kspan) RETURN
                   DO
                      c1 = 1.0_fftkind - cd
                      s1 = sd
                      DO
                         DO
                            DO
                               k2 = kk + kspan
                               ck = array(kk) - array(k2)
                               array(kk) = array(kk) + array(k2)
                               array(k2) = ck * CMPLX(c1, s1, kind=fftkind)
                               kk = k2 + kspan
                               IF (kk >= nt) EXIT
                            END DO
                            k2 = kk - nt
                            c1 = -c1
                            kk = k1 - k2
                            IF (kk <= k2) EXIT
                         END DO
                         ak = c1 - (cd * c1 + sd * s1)
                         s1 = sd * c1 - cd * s1 + s1
                         c1 = 2.0_fftkind - (ak * ak + s1 * s1)
                         s1 = s1 * c1
                         c1 = c1 * ak
                         kk = kk + jc
                         IF (kk >= k2) EXIT
                      END DO

                      k1 = k1 + 1 + 1
                      kk = (k1 - kspan) / 2 + jc
                      IF (kk > jc + jc) EXIT
                   END DO

                CASE (4) !-- transform for factor of 4
                   ispan = kspan
                   kspan = kspan / 4
                   DO
                      c1 = 1.0_fftkind
                      s1 = 0.0_fftkind
                      DO
                         DO
                            k1 = kk + kspan
                            k2 = k1 + kspan
                            k3 = k2 + kspan
                            ckp = array(kk) + array(k2)
                            ckm = array(kk) - array(k2)
                            cjp = array(k1) + array(k3)
                            cjm = array(k1) - array(k3)
                            array(kk) = ckp + cjp
                            cjp = ckp - cjp
                            IF (inv) THEN
                               ckp = ckm + CMPLX(-AIMAG(cjm), REAL(cjm), kind=fftkind)
                               ckm = ckm + CMPLX(AIMAG(cjm), -REAL(cjm), kind=fftkind)
                            ELSE
                               ckp = ckm + CMPLX(AIMAG(cjm), -REAL(cjm), kind=fftkind)
                               ckm = ckm + CMPLX(-AIMAG(cjm), REAL(cjm), kind=fftkind)
                            END IF
                            !-- avoid useless multiplies
                            IF (s1 == 0.0_fftkind) THEN
                               array(k1) = ckp
                               array(k2) = cjp
                               array(k3) = ckm
                            ELSE
                               array(k1) = ckp * CMPLX(c1, s1, kind=fftkind)
                               array(k2) = cjp * CMPLX(c2, s2, kind=fftkind)
                               array(k3) = ckm * CMPLX(c3, s3, kind=fftkind)
                            END IF
                            kk = k3 + kspan
                            IF (kk > nt) EXIT
                         END DO
                         c2 = c1 - (cd * c1 + sd * s1)
                         s1 = sd * c1 - cd * s1 + s1
                         c1 = 2.0_fftkind - (c2 * c2 + s1 * s1)
                         s1 = s1 * c1
                         c1 = c1 * c2
                         !-- values of c2, c3, s2, s3 that will get used next time
                         c2 = c1 * c1 - s1 * s1
                         s2 = 2.0_fftkind * c1 * s1
                         c3 = c2 * c1 - s2 * s1
                         s3 = c2 * s1 + s2 * c1
                         kk = kk - nt + jc
                         IF (kk > kspan) EXIT
                      END DO
                      kk = kk - kspan + 1
                      IF (kk > jc) EXIT
                   END DO
                   IF (kspan == jc) RETURN

                CASE default
                   !-- transform for odd factors
                   k = factor (ii)
                   ispan = kspan
                   kspan = kspan / k
                   SELECT CASE (k)
                      CASE (3) !-- transform for factor of 3 (optional code)
                         DO
                            DO
                               k1 = kk + kspan
                               k2 = k1 + kspan
                               ck = array(kk)
                               cj = array(k1) + array(k2)
                               array(kk) = ck + cj
                               ck = ck - 0.5_fftkind * cj
                               cj = (array(k1) - array(k2)) * s60
                               array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), kind=fftkind)
                               array(k2) = ck + CMPLX(AIMAG(cj), -REAL(cj), kind=fftkind)
                               kk = k2 + kspan
                               IF (kk >= nn) EXIT
                            END DO
                            kk = kk - nn
                            IF (kk > kspan) EXIT
                         END DO

                      CASE (5) !-- transform for factor of 5 (optional code)
                         c2 = c72 * c72 - s72 * s72
                         s2 = 2.0_fftkind * c72 * s72
                         DO
                            DO
                               k1 = kk + kspan
                               k2 = k1 + kspan
                               k3 = k2 + kspan
                               k4 = k3 + kspan
                               ckp = array(k1) + array(k4)
                               ckm = array(k1) - array(k4)
                               cjp = array(k2) + array(k3)
                               cjm = array(k2) - array(k3)
                               cc = array(kk)
                               array(kk) = cc + ckp + cjp
                               ck = ckp * c72 + cjp * c2 + cc
                               cj = ckm * s72 + cjm * s2
                               array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), kind=fftkind)
                               array(k4) = ck + CMPLX(AIMAG(cj), -REAL(cj), kind=fftkind)
                               ck = ckp * c2 + cjp * c72 + cc
                               cj = ckm * s2 - cjm * s72
                               array(k2) = ck + CMPLX(-AIMAG(cj), REAL(cj), kind=fftkind)
                               array(k3) = ck + CMPLX(AIMAG(cj), -REAL(cj), kind=fftkind)
                               kk = k4 + kspan
                               IF (kk >= nn) EXIT
                            END DO
                            kk = kk - nn
                            IF (kk > kspan) EXIT
                         END DO

                      CASE default
                         IF (k /= jf) THEN
                            jf = k
                            s1 = pi2 / k
                            c1 = COS(s1)
                            s1 = SIN(s1)
                            cosine (jf) = 1.0_fftkind
                            sine (jf) = 0.0_fftkind
                            j = 1
                            DO
                               cosine (j) = cosine (k) * c1 + sine (k) * s1
                               sine (j) = cosine (k) * s1 - sine (k) * c1
                               k = k-1
                               cosine (k) = cosine (j)
                               sine (k) = -sine (j)
                               j = j + 1
                               IF (j >= k) EXIT
                            END DO
                         END IF

                         DO
                            DO
                               k1 = kk
                               k2 = kk + ispan
                               cc = array(kk)
                               ck = cc
                               j = 1
                               k1 = k1 + kspan
                               DO
                                  k2 = k2 - kspan
                                  j = j + 1
                                  ctmp(j) = array(k1) + array(k2)
                                  ck = ck + ctmp(j)
                                  j = j + 1
                                  ctmp(j) = array(k1) - array(k2)
                                  k1 = k1 + kspan
                                  IF (k1 >= k2) EXIT
                               END DO
                               array(kk) = ck
                               k1 = kk
                               k2 = kk + ispan
                               j = 1
                               DO
                                  k1 = k1 + kspan
                                  k2 = k2 - kspan
                                  jj = j
                                  ck = cc
                                  cj = (0.0_fftkind, 0.0_fftkind)
                                  k = 1
                                  DO
                                     k = k + 1
                                     ck = ck + ctmp(k) * cosine (jj)
                                     k = k + 1
                                     cj = cj + ctmp(k) * sine (jj)
                                     jj = jj + j
                                     IF (jj > jf) jj = jj - jf
                                     IF (k >= jf) EXIT
                                  END DO
                                  k = jf - j
                                  array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), kind=fftkind)
                                  array(k2) = ck + CMPLX(AIMAG(cj), -REAL(cj), kind=fftkind)
                                  j = j + 1
                                  IF (j >= k) EXIT
                               END DO
                               kk = kk + ispan
                               IF (kk > nn) EXIT
                            END DO
                            kk = kk - nn
                            IF (kk > kspan) EXIT
                         END DO
                   END SELECT !k

                   !--  multiply by rotation factor (except for factors of 2 and 4)
                   IF (ii == nfactor) RETURN
                   kk = jc + 1
                   DO
                      c2 = 1.0_fftkind - cd
                      s1 = sd
                      DO
                         c1 = c2
                         s2 = s1
                         kk = kk + kspan
                         DO
                            DO
                               array(kk) = CMPLX(c2, s2, kind=fftkind) * array(kk)
                               kk = kk + ispan
                               IF (kk > nt) EXIT
                            END DO
                            ak = s1 * s2
                            s2 = s1 * c2 + c1 * s2
                            c2 = c1 * c2 - ak
                            kk = kk - nt + kspan
                            IF (kk > ispan) EXIT
                         END DO
                         c2 = c1 - (cd * c1 + sd * s1)
                         s1 = s1 + sd * c1 - cd * s1
                         c1 = 2.0_fftkind - (c2 * c2 + s1 * s1)
                         s1 = s1 * c1
                         c2 = c2 * c1
                         kk = kk - ispan + jc
                         IF (kk > kspan) EXIT
                      END DO
                      kk = kk - kspan + jc + 1
                      IF (kk > jc + jc) EXIT
                   END DO
             END SELECT ! Factor
          END DO

          RETURN
       END SUBROUTINE transform

       !!--++
       !!--++ Subroutine permute()
       !!--++
       !!--++    permute the results to normal order---done in two stages
       !!--++    permutation for square factors of n
       !!--++
       !!--++    (PRIVATE)
       !!--++
       !!--++ Update: February - 2005
       !!
       SUBROUTINE permute()

          perm (1) = ns
          IF (kt > 0) THEN
             k = kt + kt + 1
             IF (nfactor < k) k = k - 1
             j = 1
             perm (k + 1) = jc
             DO
                perm (j + 1) = perm (j) / factor (j)
                perm (k) = perm (k + 1) * factor (j)
                j = j + 1
                k = k - 1
                IF (j >= k) EXIT
             END DO
             k3 = perm (k + 1)
             kspan = perm (2)
             kk = jc + 1
             k2 = kspan + 1
             j = 1
             IF (npass /= ntotal) THEN
                permute_multi: DO
                   DO
                      DO
                         k = kk + jc
                         DO
                            !-- swap array(kk) <> array(k2)
                            ck = array(kk)
                            array(kk) = array(k2)
                            array(k2) = ck
                            kk = kk + 1
                            k2 = k2 + 1
                            IF (kk >= k) EXIT
                         END DO
                         kk = kk + ns - jc
                         k2 = k2 + ns - jc
                         IF (kk >= nt) EXIT
                      END DO
                      kk = kk - nt + jc
                      k2 = k2 - nt + kspan
                      IF (k2 >= ns) EXIT
                   END DO

                   DO
                      DO
                         k2 = k2 - perm (j)
                         j = j + 1
                         k2 = perm (j + 1) + k2
                         IF (k2 <= perm (j)) EXIT
                      END DO
                      j = 1

                      DO
                         IF (kk < k2) CYCLE permute_multi
                         kk = kk + jc
                         k2 = k2 + kspan
                         IF (k2 >= ns) EXIT
                      END DO
                      IF (kk >= ns) EXIT
                   END DO

                   EXIT
                END DO permute_multi
             ELSE
                permute_single: DO
                   DO
                      !-- swap array(kk) <> array(k2)
                      ck = array(kk)
                      array(kk) = array(k2)
                      array(k2) = ck
                      kk = kk + 1
                      k2 = k2 + kspan
                      IF (k2 >= ns) EXIT
                   END DO

                   DO
                      DO
                         k2 = k2 - perm (j)
                         j = j + 1
                         k2 = perm (j + 1) + k2
                         IF (k2 <= perm (j)) EXIT
                      END DO
                      j = 1
                      DO
                         IF (kk < k2) CYCLE permute_single
                         kk = kk + 1
                         k2 = k2 + kspan
                         IF (k2 >= ns) EXIT
                      END DO
                      IF (kk >= ns) EXIT
                   END DO
                   EXIT
                END DO permute_single
             END IF
             jc = k3
          END IF

          IF (ISHFT(kt, 1) + 1 >= nfactor) RETURN
          ispan = perm (kt + 1)

          !-- permutation for square-free factors of n
          j = nfactor - kt
          factor (j + 1) = 1
          DO
             factor(j) = factor(j) * factor(j+1)
             j = j - 1
             IF (j == kt) EXIT
          END DO
          kt = kt + 1
          nn = factor(kt) - 1
          j = 0
          jj = 0
          DO
             k = kt + 1
             k2 = factor(kt)
             kk = factor(k)
             j = j + 1
             IF (j > nn) EXIT !-- exit infinite loop
             jj = jj + kk
             DO WHILE (jj >= k2)
                jj = jj - k2
                k2 = kk
                k = k + 1
                kk = factor(k)
                jj = jj + kk
             END DO
             perm (j) = jj
          END DO

          !--  determine the permutation cycles of length greater than 1
          j = 0
          DO
             DO
                j = j + 1
                kk = perm(j)
                IF (kk >= 0) EXIT
             END DO
             IF (kk /= j) THEN
                DO
                   k = kk
                   kk = perm (k)
                   perm (k) = -kk
                   IF (kk == j) EXIT
                END DO
                k3 = kk
             ELSE
                perm (j) = -j
                IF (j == nn) EXIT !-- exit infinite loop
             END IF
          END DO

          !--  reorder a and b, following the permutation cycles
          DO
             j = k3 + 1
             nt = nt - ispan
             ii = nt - 1 + 1
             IF (nt < 0) EXIT !-- exit infinite loop
             DO
                DO
                   j = j-1
                   IF (perm(j) >= 0) EXIT
                END DO
                jj = jc

                DO
                   kspan = jj
                   IF (jj > maxfactor) kspan = maxfactor
                   jj = jj - kspan
                   k = perm(j)
                   kk = jc * k + ii + jj
                   k1 = kk + kspan
                   k2 = 0
                   DO
                      k2 = k2 + 1
                      ctmp(k2) = array(k1)
                      k1 = k1 - 1
                      IF (k1 == kk) EXIT
                   END DO
                   DO
                      k1 = kk + kspan
                      k2 = k1 - jc * (k + perm(k))
                      k = -perm(k)
                      DO
                         array(k1) = array(k2)
                         k1 = k1 - 1
                         k2 = k2 - 1
                         IF (k1 == kk) EXIT
                      END DO
                      kk = k2
                      IF (k == j) EXIT
                   END DO
                   k1 = kk + kspan
                   k2 = 0
                   DO
                      k2 = k2 + 1
                      array(k1) = ctmp(k2)
                      k1 = k1 - 1
                      IF (k1 == kk) EXIT
                   END DO
                   IF (jj == 0) EXIT
                END DO
                IF (j == 1) EXIT
             END DO
          END DO

          RETURN
       END SUBROUTINE permute

    End Subroutine Fftradix

    !!----
    !!---- Pure Subroutine Hfft(Aa,Ifset,Iferr)
    !!----    complex, dimension(:),    intent (in out) :: AA      In -> Contains the complex 3D array to be transformed
    !!----    integer,                  intent (in    ) :: IFSET   In ->  1, 2 Inverse Fourier Transform
    !!----                                                               -1,-2 Fourier Transform
    !!----    integer,                  intent (   out) :: IFERR  Out -> Flags to error. 0 No error
    !!----
    !!----    Performs discrete complex fourier transforms on a complex
    !!----    three dimensional array.
    !!----    This subroutine is to be used for complex, 3-dimensional
    !!----    arrays in which each dimension is a power of 2.  the
    !!----    maximum m(i) must not be less than 3 or greater than 20,
    !!----    i = 1,2,3
    !!----
    !!--..    Translation:
    !!--..               Unknown author
    !!--..               (http://risc2.numis.nwu.edu/fft/transforms/harm.f)
    !!--..               Translation to F by Javier Gonzalez-Platas
    !!--<<
    !!----  Method
    !!----
    !!----     For IFSET = -1, or -2, the fourier transform of complex
    !!----     array a is obtained.
    !!----
    !!----                  N1-1   N2-1   N3-1                L1   L2   L3
    !!----     X(J1,J2,J3)=SUM    SUM    SUM    A(K1,K2,K3)*W1  *W2  *W3
    !!----                  K1=0   K2=0   K3=0
    !!----
    !!----     where wi is the n(i) root of unity and L1=K1*J1,L2=K2*J2, L3=K3*J3
    !!----
    !!----
    !!----     For IFSET = +1, or +2, the inverse fourier transform a of
    !!----     complex array x is obtained.
    !!----
    !!----     A(K1,K2,K3)=
    !!----
    !!----               1      N1-1   N2-1   N3-1                -L1  -L2  -L3
    !!----           -------- *SUM    SUM    SUM    X(J1,J2,J3)*W1  *W2  *W3
    !!----           N1*N2*N3   J1=0   J2=0   J3=0
    !!-->>
    !!--..
    !!--..  Reference
    !!--..
    !!--..     See J.W. COOLEY and J.W. TUKEY, "an algorithm for the
    !!--..     machine calculation of complex fourier series",
    !!--..     Mathematics of Computations, Vol. 19 (apr. 1965), p. 297.
    !!--..
    !!--..
    !!---- Update: February - 2005
    !!
    Pure Subroutine Hfft(Array,Ifset,Iferr)
       !---- Arguments ----!
       complex, dimension(0:,0:,0:), intent(in out) :: Array
       integer,                      intent(in    ) :: IfSet
       integer,                      intent(   out) :: IFerr

       !---- Local variables ----!
       integer          :: n1,n2,n3,m1,m2,m3,k1,k2,k3,nx,fn,nnn
       integer          :: i,j,k,l,ii,jj,i2,i3, jj1,jj2,jj3,ip1,ip2,ip3, jp1,jp2,jp3
       integer          :: mm,mt,mtt,nt,ntv2,jstep,jstep2,idif,jdif,jlast,mtlexp,lm1exp
       Integer          :: kbit,mev,klast,jjdif,ilast,lfirst
       integer          :: jc,jc1,jd,id,il,il1,mi,kl
       integer          :: ntsq,n3vnt,n2vnt,n1vnt
       integer          :: ipp1,ipp2,ipp3,jpp1,jpp2,jpp3,jjd1,jjd2,jjd3,igo1,igo2,igo3
       integer          :: m1mt,m2mt,m3mt,minn1,minn2,minn3,ntvn1,ntvn2,ntvn3
       integer, dimension(3)              :: m,ngr,np
       integer, dimension(:), allocatable :: inv

       real(kind=cp)                               :: t,r,theta,root2,awr,awi
       real(kind=cp)                               :: ak0_0,ak0_1,ak1_0,ak1_1,ak2_0,ak2_1,ak3_0,ak3_1
       real(kind=cp), dimension(2*size(Array))     :: a
       real(kind=cp), dimension(2)                 :: w, w2, w3
       real(kind=cp), dimension(:), allocatable    :: s

       !---- Init ----!
       iferr = 0
       ngr(1)=size(Array,1)
       ngr(2)=size(Array,2)
       ngr(3)=size(Array,3)

       m=0
       do i=3,20
          j=2**i
          if (ngr(1) == j) m(1)=i
          if (ngr(2) == j) m(2)=i
          if (ngr(3) == j) m(3)=i
          if (all(m /= 0) ) exit
       end do
       if (any(m ==0)) then
          iferr=1
          return
       end if
       mm=maxval(ngr)/4

       allocate(inv(mm))
       allocate(s(mm))

       do i=0,ngr(1)-1
          do j=0,ngr(2)-1
             do k=0,ngr(3)-1
                ii=2*(k*ngr(1)*ngr(2) + j*ngr(1) + i) + 1
                a(ii)=real(Array(i,j,k))
                a(ii+1)=aimag(Array(i,j,k))
             end do
          end do
       end do

       !----
       !---- The following are parameters which may be overwritten.
       !---- The original coding in fact assumed that these variables
       !---- would be retained
       mt    =maxval(m)-2
       mt    =max(2,mt)
       if (mt >= 20) then
          iferr=1
          return
       end if

       nt=2**mt

       if (abs(ifset) <= 1) then
          !---- Computing the sin and inv tables. ----!
          ntv2=nt/2

          !---- SET UP SIN TABLE ----!
          Theta=asin(1.0)/2.0

          !---- JSTEP=2**(MT-L+1) FOR L=1 ----!
          jstep=nt

          !---- JDIF=2**(MT-L) FOR L=1 ----!
          jdif=ntv2
          s(jdif)=sin(theta)
          do l=2,mt
             theta=theta/2.0
             jstep2=jstep
             jstep=jdif
             jdif=jstep/2
             s(jdif)=sin(theta)
             jc1=nt-jdif
             s(jc1)=cos(theta)
             jlast=nt-jstep2
             if (jlast < jstep) cycle
             do j=jstep,jlast,jstep
                jc=nt-j
                jd=j+jdif
                s(jd)=s(j)*s(jc1)+s(jdif)*s(jc)
             end do
          end do

          !---- SET UP INV(J) TABLE ----!
          mtlexp=ntv2

          !---- MTLEXP=2**(MT-L). FOR L=1
          lm1exp=1

          !---- LM1EXP=2**(L-1). FOR L=1
          inv(1)=0
          do l=1,mt
             inv(lm1exp+1) = mtlexp
             do j=2,lm1exp
                jj=j+lm1exp
                inv(jj)=inv(j)+mtlexp
             end do
             mtlexp=mtlexp/2
             lm1exp=lm1exp*2
          end do

          if (ifset == 0) return
       end if

       mtt=maxval(m)-2
       root2=sqrt(2.0)
       if (mtt > mt) then
          iferr=1
          return
       end if

       m1=m(1)
       m2=m(2)
       m3=m(3)
       n1=2**m1
       n2=2**m2
       n3=2**m3

       if (ifset < 0) then
          nx= n1*n2*n3
          fn = nx
          !---- may be faster than * ----!
          nnn=nx*2
          do i = 1,nnn,2
             a(i) = a(i)/fn
             a(i+1) = -a(i+1)/fn
          end do
       end if

       np(1)= n1*2
       np(2)= np(1)*n2
       np(3)= np(2)*n3

       do id=1,3
          il = np(3)-np(id)
          il1 = il+1
          mi = m(id)
          if (mi <=0) cycle
          idif=np(id)
          kbit=np(id)
          mev = 2*(mi/2)

          if (mi > mev) then
             !---- m is odd. do l=1 case
             kbit=kbit/2
             kl=kbit-2
             do i=1,il1,idif
                klast=kl+i
                do k=i,klast,2
                   k1=k+kbit
                   ak0_0=a(k)
                   ak0_1=a(k+1)
                   ak1_0=a(k1)
                   ak1_1=a(k1+1)
                   a(k)=ak0_0+ak1_0
                   a(k+1)=ak0_1+ak1_1
                   a(k1)=ak0_0-ak1_0
                   a(k1+1)=ak0_1-ak1_1
                end do
             end do

             if (mi > 1) then
                lfirst=3
                jlast=1
             else
                cycle
             end if
          else
             !---- m is even ----!
             lfirst = 2
             jlast=0
          end if

          do l=lfirst,mi,2
             jjdif=kbit
             kbit=kbit/4
             kl=kbit-2

             do i=1,il1,idif
                klast=i+kl
                do k=i,klast,2
                   k1=k+kbit
                   k2=k1+kbit
                   k3=k2+kbit
                   ak0_0=a(k)
                   ak0_1=a(k+1)
                   ak1_0=a(k1)
                   ak1_1=a(k1+1)
                   ak2_0=a(k2)
                   ak2_1=a(k2+1)
                   ak3_0=a(k3)
                   ak3_1=a(k3+1)

                   t=ak2_0
                   ak2_0=ak0_0-t
                   ak0_0=ak0_0+t
                   t=ak2_1
                   ak2_1=ak0_1-t
                   ak0_1=ak0_1+t

                   t=ak3_0
                   ak3_0=ak1_0-t
                   ak1_0=ak1_0+t
                   t=ak3_1
                   ak3_1=ak1_1-t
                   ak1_1=ak1_1+t
                   a(k)=ak0_0+ak1_0
                   a(k+1)=ak0_1+ak1_1
                   a(k1)=ak0_0-ak1_0
                   a(k1+1)=ak0_1-ak1_1
                   a(k2)=ak2_0-ak3_1
                   a(k2+1)=ak2_1+ak3_0
                   a(k3)=ak2_0+ak3_1
                   a(k3+1)=ak2_1-ak3_0
                end do
             end do

             if (jlast > 0) then
                jj=jjdif   +1
                ilast= il +jj
                do i = jj,ilast,idif
                   klast = kl+i
                   do k=i,klast,2
                      k1 = k+kbit
                      k2 = k1+kbit
                      k3 = k2+kbit
                      ak0_0=a(k)
                      ak0_1=a(k+1)
                      ak1_0=a(k1)
                      ak1_1=a(k1+1)
                      ak2_0=a(k2)
                      ak2_1=a(k2+1)
                      ak3_0=a(k3)
                      ak3_1=a(k3+1)
                      r =-ak2_1
                      t = ak2_0
                      ak2_0 = ak0_0-r
                      ak0_0 = ak0_0+r
                      ak2_1=ak0_1-t
                      ak0_1=ak0_1+t

                      awr=ak1_0-ak1_1
                      awi = ak1_1+ak1_0
                      r=-ak3_0-ak3_1
                      t=ak3_0-ak3_1
                      ak3_0=(awr-r)/root2
                      ak3_1=(awi-t)/root2
                      ak1_0=(awr+r)/root2
                      ak1_1=(awi+t)/root2
                      a(k)=ak0_0+ak1_0
                      a(k+1)=ak0_1+ak1_1
                      a(k1)=ak0_0-ak1_0
                      a(k1+1)=ak0_1-ak1_1
                      a(k2)=ak2_0-ak3_1
                      a(k2+1)=ak2_1+ak3_0
                      a(k3)=ak2_0+ak3_1
                      a(k3+1)=ak2_1-ak3_0
                   end do
                end do

                if (jlast-1 > 0) then
                   jj= jj + jjdif
                   do j=2,jlast
                      i=inv(j+1)
                      w(1)=s(nt-i)
                      w(2)=s(i)
                      k=i+i
                      k1=nt-k

                      select case (k1)
                         case (:-1)
                            !---- 2*i is in second quadrant ----!
                            k2 = k1+nt
                            k1=-k1
                            w2(1)=-s(k1)
                            w2(2)=s(k2)
                         case (0)
                            w2(1)=0.0
                            w2(2)=1.0
                         case (1:)
                            !---- 2*i is in first quadrant ----!
                            w2(1)=s(k1)
                            w2(2)=s(k)
                      end select

                      k1=i+k
                      k2=nt-k1

                      select case (k2)
                         case (:-1)
                            k=k2+nt
                            select case (k)
                               case (:-1)
                                  !---- 3*i in third quadrant ----!
                                  k1=nt+k
                                  k = -k
                                  w3(1)=-s(k1)
                                  w3(2)=-s(k)
                               case (0)
                                  w3(1)=-1.0
                                  w3(2)=0.0
                               case (1:)
                                  !---- i3 in second quadrant ----!
                                  k2=-k2
                                  w3(1)=-s(k2)
                                  w3(2)=s(k)
                            end select
                         case (0)
                            w3(1)=0.0
                            w3(2)=1.0
                         case (1:)
                            !---- i3 in first quadrant ----!
                            w3(1)=s(k2)
                            w3(2)=s(k1)
                      end select

                      ilast=il+jj
                      do i=jj,ilast,idif
                         klast=kl+i
                         do k=i,klast,2
                            k1=k+kbit
                            k2=k1+kbit
                            ak0_0=a(k)
                            ak0_1=a(k+1)
                            ak1_0=a(k1)
                            ak1_1=a(k1+1)
                            ak2_0=a(k2)
                            ak2_1=a(k2+1)
                            ak3_0=a(k2+kbit)
                            ak3_1=a(k2+kbit+1)
                            r=ak2_0*w2(1)-ak2_1*w2(2)
                            t=ak2_0*w2(2)+ak2_1*w2(1)
                            ak2_0=ak0_0-r
                            ak0_0=ak0_0+r
                            ak2_1=ak0_1-t
                            ak0_1=ak0_1+t

                            r=ak3_0*w3(1)-ak3_1*w3(2)
                            t=ak3_0*w3(2)+ak3_1*w3(1)
                            awr=ak1_0*w(1)-ak1_1*w(2)
                            awi=ak1_0*w(2)+ak1_1*w(1)
                            ak3_0=awr-r
                            ak3_1=awi-t
                            ak1_0=awr+r
                            ak1_1=awi+t
                            a(k)=ak0_0+ak1_0
                            a(k+1)=ak0_1+ak1_1
                            a(k1)=ak0_0-ak1_0
                            a(k1+1)=ak0_1-ak1_1
                            a(k2)=ak2_0-ak3_1
                            a(k2+1)=ak2_1+ak3_0
                            a(k2+kbit)=ak2_0+ak3_1
                            a(k2+kbit+1)=ak2_1-ak3_0
                         end do
                      end do
                      jj=jjdif+jj
                   end do ! j
                end if
             end if
             jlast=4*jlast+3

          end do ! end l
       end do ! id

       !---- We now have the complex fourier sums but their addresses are
       !---- bit-reversed. The following routine puts them in order.
       ntsq=nt*nt
       m3mt=m3-mt

       if (m3mt >=0) then
          !---- m3 gr. or eq. mt ----!
          igo3=1
          n3vnt=n3/nt
          minn3=nt
       else
           !---- m3 less than mt ----!
          igo3=2
          n3vnt=1
          ntvn3=nt/n3
          minn3=n3
       end if
       jjd3 = ntsq/n3
       m2mt=m2-mt

       if (m2mt >=0) then
          !---- m2 gr. or eq. mt ----!
          igo2=1
          n2vnt=n2/nt
          minn2=nt
       else
          !---- m2 less than mt ----!
          igo2 = 2
          n2vnt=1
          ntvn2=nt/n2
          minn2=n2
       end if
       jjd2=ntsq/n2
       m1mt=m1-mt

       if (m1mt >= 0) then
          !---- m1 gr. or eq. mt ----!
          igo1=1
          n1vnt=n1/nt
          minn1=nt
       else
          !---- m1 less than mt ----!
          igo1=2
          n1vnt=1
          ntvn1=nt/n1
          minn1=n1
       end if

       jjd1=ntsq/n1
       jj3=1
       j=1

       do jpp3=1,n3vnt
          ipp3=inv(jj3)
          do jp3=1,minn3
             select case (igo3)
                case (1)
                   ip3=inv(jp3)*n3vnt
                case (2)
                   ip3=inv(jp3)/ntvn3
             end select
             i3=(ipp3+ip3)*n2
             jj2=1
             do jpp2=1,n2vnt
                ipp2=inv(jj2)+i3
                do jp2=1,minn2
                   select case (igo2)
                      case (1)
                         ip2=inv(jp2)*n2vnt
                      case (2)
                         ip2=inv(jp2)/ntvn2
                   end select
                   i2=(ipp2+ip2)*n1
                   jj1=1
                   do jpp1=1,n1vnt
                      ipp1=inv(jj1)+i2
                      do jp1=1,minn1
                         select case (igo1)
                            case (1)
                               ip1=inv(jp1)*n1vnt
                            case (2)
                               ip1=inv(jp1)/ntvn1
                         end select
                         i=2*(ipp1+ip1)+1
                         if (j < i) then
                            t=a(i)
                            a(i)=a(j)
                            a(j)=t
                            t=a(i+1)
                            a(i+1)=a(j+1)
                            a(j+1)=t
                         end if
                         j=j+2
                      end do
                      jj1=jj1+jjd1
                   end do !jpp1
                end do ! jp2
                jj2=jj2+jjd2
             end do ! jpp2
          end do ! jp3
          jj3 = jj3+jjd3
       end do ! jpp3

       if (ifset < 0) then
          nnn=nx*2
          do i=2,nnn,2
             a(i)=-a(i)
          end do
       end if

       do i=0,ngr(1)-1
          do j=0,ngr(2)-1
             do k=0,ngr(3)-1
                ii=2*(k*ngr(1)*ngr(2) + j*ngr(1) + i) + 1
                Array(i,j,k)=cmplx(a(ii),a(ii+1))
             end do
          end do
       end do

       return
    End Subroutine HFFT

    !!----
    !!---- Pure Subroutine Sfft(Array,Mode, Iferr)
    !!----    complex, dimension(:),     intent (in out) :: Array  !  In -> real array containing real parts of transform
    !!----    character(len=*),optional, intent (in)     :: Mode   !  In -> type="INV"    : backward transform
    !!----                                                                  Absent or whatever : forward transform
    !!----    integer, optional,         intent (   out) :: IFERR  ! Out -> Flags to error. 0 No error
    !!----
    !!----    This routine is a slight modification of a complex split
    !!----    radix FFT routine presented by C.S. Burrus.
    !!--..    The original program header is shown below.
    !!--<<
    !!----    The forward transform computes
    !!----        X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)
    !!----
    !!----    The backward transform computes
    !!----        x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)
    !!-->>
    !!--.. Authors:
    !!--..
    !!--..      Steve Kifowit, 9 July 1997
    !!--..      Traslation to Fortran90 by Javier Gonzalez-Platas
    !!--..
    !!--.. Refrences:
    !!--..
    !!--..  A Duhamel-Hollman Split-Radix DIF FFT
    !!--..  Electronics Letters, January 5, 1984
    !!----
    !!---- Update: February - 2005
    !!
    Pure Subroutine Sfft( Array, Mode,Iferr )
       !---- Arguments ----!
       complex, dimension(:),      intent(in out) :: Array
       character(len=*), optional, intent(in)     :: Mode
       integer,          optional, intent(   out) :: Iferr

       !---- Local variables ----!
       integer                            :: i, j, k, n, m, n1, n2, n4, is, id, i0, i1, i2, i3
       real(kind=cp)                      :: r1, r2, s1, s2, s3, xt
       real(kind=cp)                      :: e, a, a3, cc1, ss1, cc3, ss3
       real(kind=cp), parameter           :: twopi=6.2831853071795864769
       real(kind=cp), dimension(size(Array)) :: x
       real(kind=cp), dimension(size(Array)) :: y

       n=size(Array)
       m=0
       do i=1,20
          if (2**i /= n) cycle
          m=i
          exit
       end do

       if(present(iferr)) then
         iferr=0
         if (n == 1 .or. m ==0) then
            iferr=1
            return
         end if
       end if

       do i=1,n
          x(i)=real(Array(i))
          y(i)=aimag(Array(i))
       end do

       if(present(Mode)) then
         if (Mode == "INV") y=-y
       end if

       n2 = 2 * n
       do k = 1, m-1
          n2 = n2 / 2
          n4 = n2 / 4
          e = twopi / n2
          a = 0.0
          do j = 1, n4
             a3 = 3 * a
             cc1 = cos( a )
             ss1 = sin( a )
             cc3 = cos( a3 )
             ss3 = sin( a3 )
             a = j * e
             is = j
             id = 2 * n2
             do
                do i0 = is, n-1, id
                   i1 = i0 + n4
                   i2 = i1 + n4
                   i3 = i2 + n4
                   r1 = x(i0) - x(i2)
                   x(i0) = x(i0) + x(i2)
                   r2 = x(i1) - x(i3)
                   x(i1) = x(i1) + x(i3)
                   s1 = y(i0) - y(i2)
                   y(i0) = y(i0) + y(i2)
                   s2 = y(i1) - y(i3)
                   y(i1) = y(i1) + y(i3)
                   s3 = r1 - s2
                   r1 = r1 + s2
                   s2 = r2 - s1
                   r2 = r2 + s1
                   x(i2) = r1 * cc1 - s2 * ss1
                   y(i2) = - s2 * cc1 - r1 * ss1
                   x(i3) = s3 * cc3 + r2 * ss3
                   y(i3) = r2 * cc3 - s3 * ss3
                end do
                is = 2 * id - n2 + j
                id = 4 * id
                if (is >= n) exit
             end do
          end do
       end do

       ! LAST STAGE, LENGTH-2 BUTTERFLY
       is = 1
       id = 4
       do
          do i0 = is, n, id
             i1 = i0 + 1
             r1 = x(i0)
             x(i0) = r1 + x(i1)
             x(i1) = r1 - x(i1)
             r1 = y(i0)
             y(i0) = r1 + y(i1)
             y(i1) = r1 - y(i1)
          end do
          is = 2 * id - 1
          id = 4 * id
          if (is >= n) exit
       end do

       ! BIT REVERSE COUNTER
       j = 1
       n1 = n - 1
       do i = 1, n1
          if (i < j) then
             xt = x(j)
             x(j) = x(i)
             x(i) = xt
             xt = y(j)
             y(j) = y(i)
             y(i) = xt
          end if
          k = n / 2

          do
             if (k >=j) exit
             j = j - k
             k = k / 2
          end do
          j = j + k
       end do

       if(present(Mode)) then
         if (Mode == "INV") then
          x=x/n
          y=-y/n
         end if
       end if

       do i=1,n
          Array(i)=cmplx(x(i),y(i))
       end do

       return
    End Subroutine SFFT

 End Module CFML_FFT
