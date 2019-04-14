!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2019  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_FFT
!!----   INFO: FFT Calculations Routines
!!----
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
!!----
!!
Module CFML_FFT
   !---- Use Modules ----!
   use CFML_GlobalDeps, only : CP, SP, DP, ERR_CFML

   !---- Local Variables ----!
   implicit none

   private

   !---- List of public functions ----!
   public :: Convol, Convol_Peaks, &
             FFT, F_FFT

   !---- List of public subroutines ----!
   public :: SFFT, HFFT

   !---- Definitions ----!
   
   !--------------------!
   !---- PARAMETERS ----!
   !--------------------!
   integer, parameter:: FFTKIND = Kind(0.0_cp)                             ! Default Precicion for FFT Variables.
                                                                           ! Adjust Here For Other Precisions

   real(kind=FFTKIND), parameter:: Cos72 = 0.30901699437494742_Fftkind  ! COS72
   real(kind=FFTKIND), parameter:: Pi    = 3.14159265358979323_Fftkind  ! Pi
   real(kind=FFTKIND), parameter:: Sin60 = 0.86602540378443865_Fftkind  ! SIN60
   real(kind=FFTKIND), parameter:: Sin72 = 0.95105651629515357_Fftkind  ! SIN72


   !-------------------!
   !---- VARIABLES ----!
   !-------------------!
   !integer, save :: StatusF                  ! Information on FFT Routines

   
   !!----
   !!---- TYPE, public :: Points_Interval_Type
   !!----
   Type, public :: Points_Interval_Type
      integer       :: Np                      ! Number of Point
      real(kind=cp) :: Low                     ! Lower range value
      real(kind=cp) :: High                    ! Upper range value
   End Type Points_Interval_Type

   !---- Overlapp ----!
   Interface Fft
      Module Procedure Fft1D
      Module Procedure Fft2D
      Module Procedure Fft3D
      Module Procedure Fft4D
      Module Procedure Fft5D
      Module Procedure Fft6D
      Module Procedure Fft7D
   End Interface
   
   !---- General Procedures ----!
   Interface
      Module Pure Function Convol(F, Pf, G, Pg, Interval) Result(Conv)
         !---- Arguments ----!
         real(kind=cp),dimension(:),          intent(in) :: pf
         real(kind=cp),dimension(:),          intent(in) :: pg
         type(Points_Interval_Type),          intent(in) :: interval
         real(kind=cp), dimension(interval%np)           :: conv
         
         Interface 
            Pure function f(x,parf)  result (vf)
               use CFML_GlobalDeps, only: cp
               real(kind=cp),              intent(in) :: x
               real(kind=cp), dimension(:),intent(in) :: parf
               real(kind=cp)                          :: vf
            End Function F
         
            Pure Function G(X,Parg)  Result (Vg)
               use CFML_GlobalDeps, only: cp
               real(kind=cp),              intent(in) :: x
               real(kind=cp), dimension(:),intent(in) :: parg
               real(kind=cp)                          :: vg
            End Function G
          End Interface 
      End Function Convol
      
      Module Pure Function Convol_Peaks(F,Pf,G,Pg,Wd,Np) Result(Conv)
         !---- Arguments ----!
         real(kind=cp),dimension(:),          intent(in) :: pf !Parameters of the function f (starting with FWHM)
         real(kind=cp),dimension(:),          intent(in) :: pg !Parameters of the function g (starting with FWHM)
         real(kind=cp),                       intent(in) :: wd !Number of times a FWHM of the f-function to calculate range
         integer,                             intent(in) :: np !Number of points (it is increased internally up to the closest power of 2)
         real(kind=cp), dimension(np)                    :: conv

         Interface 
            Pure Function F(X,Parf)  Result (Vf)
               use CFML_GlobalDeps, only: cp
               real(kind=cp),              intent(in) :: x
               real(kind=cp), dimension(:),intent(in) :: parf
               real(kind=cp)                          :: vf
            End Function F

            Pure Function G(X,Parg) Result (Vg)
               use CFML_GlobalDeps, only: cp
               real(kind=cp),              intent(in) :: x
               real(kind=cp), dimension(:),intent(in) :: parg
               real(kind=cp)                          :: vg
            End Function G
         End Interface 
      End Function Convol_Peaks 
      
      Module Pure Function F_FFT(Array, Mode ) Result(fft)
         !---- Arguments ----!
         complex(kind=cp), dimension(:),      intent(in) :: Array  !  In -> real array containing real parts of transform
         character(len=*),          optional, intent(in) :: Mode   !  In -> type="INV" : backward transform. Absent or whatever -> forward transform
         complex(kind=cp), dimension(size(Array))        :: fft 
      End Function F_FFT
      
      Module Function Fft1D(Array, Dim, Inv) Result(Ft)
         !--- formal parameters
         complex(fftkind), dimension(:), intent(in)             :: array     ! In -> Complex array
         integer,          dimension(:), intent(in),  optional  :: dim       ! In -> array containing the dimensions to be transformed
         logical,                        intent(in),  optional  :: inv       ! In -> If .true., inverse transformation will be performed.
         complex(fftkind), dimension(size(array, 1))            :: ft
      End Function Fft1D
      
      Module Function Fft2D(Array, Dim, Inv) Result(Ft)
         !--- formal parameters
         complex(fftkind), dimension(:,:), intent(in)                :: array         ! In -> Complex array
         integer,          dimension(:),   intent(in),  optional     :: dim           ! In -> array containing the dimensions to be transformed
         logical,                          intent(in),  optional     :: inv 
         complex(fftkind), dimension(size(array, 1), size(array, 2)) :: ft
      End Function Fft2D
      
      Module Function Fft3D(Array, Dim, Inv) Result(Ft)
         !--- formal parameters
         complex(fftkind), dimension(:,:,:), intent(in)               :: array        ! In -> Complex array
         integer,          dimension(:),     intent(in),  optional    :: dim          ! In -> array containing the dimensions to be transformed
         logical,                            intent(in),  optional    :: inv 
         complex(fftkind), &
            dimension(size(array, 1), size(array, 2), size(array, 3)) :: ft
      End Function Fft3D
      
      Module Function Fft4D(Array, Dim, Inv) Result(Ft)
         !--- formal parameters
         complex(fftkind), dimension(:,:,:,:), intent(in)                   :: array         ! In -> Complex array
         integer,          dimension(:),       intent(in),  optional        :: dim           ! In -> array containing the dimensions to be transformed
         logical,                              intent(in),  optional        :: inv  
         complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4)) :: ft
      End Function Fft4D 
      
      Module Function Fft5D(Array, Dim, Inv) Result(Ft)
         !--- formal parameters
         complex(fftkind), dimension(:,:,:,:,:), intent(in)                   :: array           ! In -> Complex array
         integer,          dimension(:),         intent(in),  optional        :: dim             ! In -> array containing the dimensions to be transformed
         logical,                                intent(in),  optional        :: inv
         complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5))                                                   :: ft
      End Function Fft5D
      
      Module Function Fft6D(Array, Dim, Inv) Result(Ft)
         !--- formal parameters
         complex(fftkind), dimension(:,:,:,:,:,:), intent(in)                 :: array         ! In -> Complex array
         integer,          dimension(:),           intent(in),  optional      :: dim           ! In -> array containing the dimensions to be transformed
         logical,                                  intent(in),  optional      :: inv
         complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5), size(array, 6))                                   :: ft
      End Function Fft6D 
      
      Module Function Fft7D(Array, Dim, Inv) Result(Ft)
         !--- formal parameters
         complex(fftkind), dimension(:,:,:,:,:,:,:), intent(in)               :: array      ! In -> Complex array
         integer,          dimension(:),             intent(in),  optional    :: dim        ! In -> array containing the dimensions to be transformed
         logical,                                    intent(in),  optional    :: inv 
         complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5), size(array, 6), size(array, 7))                   :: ft
      End Function Fft7D
      
      Module Subroutine Fftn(Array, Shape, Dim, Inv, Stat)
         !--- formal parameters
         complex(fftkind), dimension(:), intent(in out)       :: array
         integer,          dimension(:), intent(in)           :: shape
         integer,          dimension(:), intent(in),  optional:: dim
         logical,                        intent(in),  optional:: inv
         integer,                        intent(out), optional:: stat  
      End Subroutine Fftn
      
      Module Pure Subroutine Hfft(Array,Ifset,Iferr)
         !---- Arguments ----!
         complex(kind=cp), dimension(0:,0:,0:), intent(in out) :: Array       ! In -> Contains the complex 3D array to be transformed
         integer,                               intent(in    ) :: IfSet       ! In -> 1,2 Inverse Fourier Transform; -1,-2 Fourier Transform
         integer,                               intent(   out) :: IFerr       ! Out -> Flags to error. 0 No error
      End Subroutine Hfft
      
      Module Pure Subroutine Sfft( Array, Mode,Iferr )
         !---- Arguments ----!
         complex(kind=cp), dimension(:),      intent(in out) :: Array   ! In -> real array containing real parts of transform
         character(len=*), optional,          intent(in)     :: Mode    ! In -> type="INV" : backward transform; Absent or whatever : forward transform
         integer,          optional,          intent(   out) :: Iferr   ! Out -> Flags to error. 0 No error
      End Subroutine Sfft
      
   End Interface
   
 Contains
    !!--++
    !!--++ SUBROUTINE FFTRADIX
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
       !!--++ SUBROUTINE FACTORIZE
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
       !!--++ SUBROUTINE TRANSFORM
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
       !!--++ SUBROUTINE PERMUTE
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
   
 End Module CFML_FFT
