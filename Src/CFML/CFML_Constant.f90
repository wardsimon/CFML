!!----
!!---- Copyleft(C) 1999-2009,              Version: 4.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: CFML_Constants
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
!!
Module CFML_Constants

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
   integer, parameter :: OPS= 1    ! Windows
   !integer, parameter :: OPS= 2    ! Linux
   !integer, parameter :: OPS= 3    ! MacOS
   !integer, parameter :: OPS= 4    ! Solaris

   !!----
   !!---- OPS_NAME
   !!----   Character variable containing the name of the operating system
   !!----   This is a variable set by the user of the library for the case
   !!----   that there is no external library with a procedure for getting
  !!----   the operating system.
   !!----
   !!---- Update: March 2009
   !!
   character(len=*), parameter :: OPS_NAME="Windows"
   !character(len=*), parameter :: OPS_NAME="Linux"
   !character(len=*), parameter :: OPS_NAME="MacOS"
   !character(len=*), parameter :: OPS_NAME="Solaris"

   !!----
   !!---- OPS_SEP
   !!----   ASCII code of directory separator character
   !!----   Here it is written explicitly as a character variable
   !!----
   !!---- Update: March 2009
   !!
   character(len=*), parameter :: OPS_SEP="\"
   !character(len=*), parameter :: OPS_SEP="/"
   !character(len=*), parameter :: OPS_SEP="/"
   !character(len=*), parameter :: OPS_SEP=":"  !Old MacOS

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

End Module CFML_Constants
