!!----
!!---- Copyleft(C) 1999-2010,              Version: 4.1
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: CFML_GlobalDeps  (Mac OS version)
!!----   INFO: Precision for for CrysFML library and Operating System information
!!----         All the global variables defined in this module are implicitly public.
!!----
!!---- HISTORY
!!--..    Update: January - 2009
!!--..
!!---- VARIABLES
!!--..
!!--..    Operating system
!!--..
!!----    OPS
!!----    OPS_NAME
!!----    OPS_SEP
!!--..
!!--..    Precision Data
!!--..
!!----    SP
!!----    DP
!!----    CP
!!--..
!!--..    Trigonometric
!!--..
!!----    PI
!!----    TO_DEG
!!----    TO_RAD
!!----    TPI
!!--..
!!--..    Numeric
!!--..
!!----    DEPS
!!----    EPS
!!--..
!!--..    Compiler
!!--..
!!----    COMPNAME
!!--..
!!---- FUNCTIONS
!!--..
!!----    DIRECTORY_EXISTS
!!----
!!
Module CFML_GlobalDeps

   !---- Variables ----!
   implicit None

   public

   !------------------------------------!
   !---- Operating System variables ----!
   !------------------------------------!

   !!----
   !!---- OPS
   !!----   Integer variable 1: Windows, 2: Mac, 3: Linux, ....
   !!----   This is a variable set by the user of the library for the case
   !!----   that there is no external library with a procedure for getting
   !!----   the operating system.
   !!----
   !!---- Update: March 2009
   !!
   integer, parameter :: OPS= 3    ! MacOS

   !!----
   !!---- OPS_NAME
   !!----   Character variable containing the name of the operating system
   !!----   This is a variable set by the user of the library for the case
   !!----   that there is no external library with a procedure for getting
   !!----   the operating system.
   !!----
   !!---- Update: March 2009
   !!
   character(len=*), parameter :: OPS_NAME="MacOS"

   !!----
   !!---- OPS_SEP
   !!----   ASCII code of directory separator character
   !!----   Here it is written explicitly as a character variable
   !!----
   !!---- Update: March 2009
   !!
   character(len=*), parameter :: OPS_SEP="/"

   !------------------------------!
   !---- Precision Parameters ----!
   !------------------------------!

   !!----
   !!---- SP
   !!----    SP: Single precision ( sp = selected_real_kind(6,30) )
   !!----
   !!---- Update: January - 2009
   !!
   integer, parameter :: sp = selected_real_kind(6,30)

   !!----
   !!---- DP
   !!----    DP: Double precision ( dp = selected_real_kind(14,150) )
   !!----
   !!---- Update: January - 2009
   !!
   integer, parameter :: dp = selected_real_kind(14,150)

   !!----
   !!---- CP
   !!----    CP: Current precision
   !!----
   !!---- Update: January - 2009
   !!
   integer, parameter :: cp = sp

   !----------------------------------!
   !---- Trigonometric Parameters ----!
   !----------------------------------!

   !!----
   !!---- PI
   !!----    real(kind=dp), parameter ::  pi = 3.141592653589793238463_dp
   !!----
   !!----    Pi value
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter ::  pi = 3.141592653589793238463_dp

   !!----
   !!---- TO_DEG
   !!----    real(kind=dp), parameter ::  to_DEG = 180.0_dp/pi
   !!----
   !!----    Conversion from Radians to Degrees
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter ::  to_DEG  = 180.0_dp/pi

   !!----
   !!---- TO_RAD
   !!----    real(kind=dp), parameter ::  to_RAD  = pi/180.0_dp
   !!----
   !!----    Conversion from Degrees to Radians
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter ::  to_RAD  = pi/180.0_dp

   !!----
   !!---- TPI
   !!----  real(kind=dp), parameter ::  tpi = 6.283185307179586476925_dp
   !!----
   !!----  2.0*Pi value
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter ::  tpi = 6.283185307179586476925_dp

   !----------------------------!
   !---- Numeric Parameters ----!
   !----------------------------!

   !!----
   !!---- DEPS
   !!----    real(kind=dp), parameter :: deps=0.00000001_dp
   !!----
   !!----    Epsilon value use for comparison of real numbers
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter, public :: deps=0.00000001_dp

   !!----
   !!----  EPS
   !!----     real(kind=cp), public ::  eps=0.00001_cp
   !!----
   !!----     Epsilon value use for comparison of real numbers
   !!----
   !!----  Update: January - 2009
   !!
   real(kind=cp),  parameter, public  ::  eps=0.00001_cp

 Contains

   !-------------------!
   !---- Functions ----!
   !-------------------!

   !!----
   !!---- Function Directory_Exists(Dirname) Result(info)
   !!----    character(len=*), intent(in) :: Dirname
   !!----    logical                      :: info
   !!----
   !!---- Generic function dependent of the compiler that return
   !!---- a logical value if a directory exists or not.
   !!----
   !!---- Update: April - 2009
   !!
   Function Directory_Exists(Dirname) Result(info)
      !---- Argument ----!
      character(len=*), intent(in) :: Dirname
      logical                      :: info

      !---- Local Variables ----!
      character(len=512) :: linea
      integer            :: nlong

      ! Init value
      info=.false.

      linea=trim(dirname)
      nlong=len_trim(linea)
      if (nlong ==0) return

      if (linea(nlong:nlong) /= ops_sep) linea=trim(linea)//ops_sep

      ! All compilers except Intel
      !inquire(file=trim(linea)//'.' , exist=info)

      ! Intel
      inquire(directory=trim(linea), exist=info)

      return
   End Function Directory_Exists


End Module CFML_GlobalDeps
