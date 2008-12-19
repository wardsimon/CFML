!!----
!!---- Copyleft(C) 1999-2009,              Version: 4.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: CFML_Constants
!!----   INFO: Precision for for CrysFML library
!!----
!!---- HISTORY
!!--..    Update: January - 2009
!!--..
!!---- VARIABLES
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
   !!----    Conversion from Rad to Degree
   !!----
   !!---- Update: January - 2009 
   !!
   real(kind=dp), parameter ::  to_DEG  = 180.0_dp/pi

   !!----
   !!---- TO_RAD
   !!----    real(kind=dp), parameter ::  to_RAD  = pi/180.0_dp
   !!----
   !!----    Conversion from Degree to Rad
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
   !!----    Epsilon value
   !!----
   !!---- Update: January - 2009 
   !!
   real(kind=dp), parameter, public :: deps=0.00000001_dp

   !!----
   !!----  EPS
   !!----     real(kind=cp), public ::  eps=0.00001_cp
   !!----
   !!----     Epsilon value
   !!----
   !!----  Update: January - 2009 
   !!
   real(kind=cp),  parameter, public  ::  eps=0.00001_cp
   
End Module CFML_Constants
