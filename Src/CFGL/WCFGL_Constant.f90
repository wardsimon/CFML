module WCFGL_constant
!------------------------------------------------
! Written by Laurent C.Chapon 
! March 2004. 
! Updated : 
!------------------------------------------------ 
implicit none 
  integer, parameter, public :: sp = selected_real_kind(6,30)
  integer, parameter, public :: dp = selected_real_kind(14,80)
  real, parameter :: PI=3.1415926535897932384626433832795_dp
  real, parameter :: deg2rad=0.017453292519943295769236907684886_dp
  real, parameter :: rad2deg=57.295779513082320876798154814105_dp

end module WCFGL_constant 