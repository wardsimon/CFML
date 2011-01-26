!     Last change:  TR    9 May 2006    4:15 pm
 module  Neutrons_data

 use atome_module, ONLY : atom
  implicit none


  contains

  subroutine definition_data_neutrons()
    implicit none

 ! bcoh (b coherent):
 !    =>  neutron coherent scattering length (in 10-12 cm)
 ! SEDinc (Section Efficace de Diffusion Incoherente)
 !    => Incoherent scattering neutron cross-section (in barns, i.e. 10-24 cm2)
 ! SEA (Section Efficace d'Absorption
 !    => neutron absorption cross-section (in barns), for v = 2200 m/s (v = neutrons velocity)
 !            rk: l (A) = 3.95 / v (km/s)
 !
 ! extrait de: V.F. Sears
 !             Neutron News, vol.3 n°3, 1992, 26-37

! coherent scattering length (in barns)
      atom(1:201)%bcoh = 0.
      atom(  1: 10)%bcoh = (/ -.3739,  .326,   -.19,    .779,    .53,   .6646,  .936,   .5803,  .5654,  .4566 /)
      atom( 11: 20)%bcoh = (/  .363,   .5375,   .3449,  .41491,  .513,  .2847,  .9577,  .1909,  .367,   .47   /)
      atom( 21: 30)%bcoh = (/ 1.229,  -.3438,  -.03824, .3635,  -.373,  .945,   .249,  1.03,    .7718,  .568  /)
      atom( 31: 40)%bcoh = (/  .7288,  .8185,   .658,   .797,    .6795, .781,   .709,   .702,   .775,   .716  /)
      atom( 41: 50)%bcoh = (/  .7054,  .6715,   .68,    .703,    .588 , .591,   .5922 , .487 ,  .4065,  .6225 /)
      atom( 51: 60)%bcoh = (/  .557,   .58,     .528,   .492,    .542,  .507,   .824,   .484,   .458,   .769  /)
      atom( 61: 70)%bcoh = (/ 1.26,    .08,     .722,   .65,     .738, 1.69,    .801,   .779,   .707,  1.243  /)
      atom( 71: 80)%bcoh = (/  .721,   .777,    .691,   .486,    .92,  1.07,   1.06,    .96,    .763,  1.2692 /)
      atom( 81: 90)%bcoh = (/  .8776,  .9405,   .8532,  .0,      .0,    .0,     .0,    1.,      .0,    1.031  /)
      atom( 91:100)%bcoh = (/  .91,    .8417,  1.055,  1.41,    0.83,   .95,    .0,     .0,     .0,    .0     /)
      atom(101:200)%bcoh = 0.
      atom(201)%bcoh     = .6671

!incoherent scattering cross_section:
      atom(1:201)%SEDinc = 0.
      atom(  1: 10)%SEDinc = (/ 80.26,  0.,    .92,    .0018,  1.7,     .001,   .5,  0.,     .0008, .008  /)
      atom( 11: 20)%SEDinc = (/  1.62,   .08,  .0082,  .004,    .005,   .007,  5.3,   .225,  .27,   .05   /)
      atom( 21: 30)%SEDinc = (/  4.5,   2.87, 5.08,   1.837,    .4,     .4,    4.8,  5.2,    .55,   .077  /)
      atom( 31: 40)%SEDinc = (/   .16,   .17,  .06,    .32,     .1,     .01,    .5,   .06,   .15,   .02   /)
      atom( 41: 50)%SEDinc = (/   .0024, .04,  .5,     .4,      .3,     .093,   .58, 3.46,   .54,   .022  /)
      atom( 51: 60)%SEDinc = (/  0.,     .09,  .31,   0.,       .21,    .15,   1.13, 0.,     .015, 9.2    /)
      atom( 61: 70)%SEDinc = (/  1.3,  39.,   2.5,    0.,       .004, 54.4,     .36, 1.1,    .1,   4.     /)
      atom( 71: 80)%SEDinc = (/   .7,   2.6,   .01,   1.63,     .9,     .3,    0.,    .13,   .43,  6.6    /)
      atom( 81: 90)%SEDinc = (/  .21,    .003, .0084, 0.,      0.,     0.,     0.,   0.,    0.,    0.     /)
      atom( 91: 95)%SEDinc = (/  .1,     .005, .5,    0.,       .3/)
      atom(96:200)%SEDinc = 0.
      atom(201)%SEDinc = 2.05

!absorption cross sections:

      atom(1:210)%SEA = 0.
      atom(  1: 10)%SEA = (/    .3326,    .00747,  70.5,       .0076,   767.,      .0035,   1.9,      .00019,    .0096,    .0039 /)
      atom( 11: 20)%SEA = (/    .53,      .063,      .231,     .171,       .172,   .53,    33.5,      .675,     2.1,       .43   /)
      atom( 21: 30)%SEA = (/  27.5,      6.09,      5.08,     3.05,      13.3,    2.56 ,   37.18,    4.49,      3.78,     1.11   /)
      atom( 31: 40)%SEA = (/   2.75,     2.2,       4.5,     11.7,        6.9,   25.,        .38,    1.28,      1.28,      .185  /)
      atom( 41: 50)%SEA = (/   1.15,     2.48,     20.,       2.56,     144.8,    6.9,     63.3,  2520.,      193.8,       .626  /)
      atom( 51: 60)%SEA = (/   4.91,     4.7,       6.15,    23.9,       29.,     1.1,      8.97,     .63,     11.5,     50.5    /)
      atom( 61: 70)%SEA = (/ 168.4,   5922.,     4530.,   49700.,        23.4,  994.,      64.7,   159.,      100.,      34.8    /)
      atom( 71: 80)%SEA = (/  74.,      74.,       20.6,     18.3,       89.7,   16.,     425.,     10.3,      98.65,   372.3    /)
      atom( 81: 90)%SEA = (/   3.43,      .171,      .0338,   0.,         0.,     0.,       0.,     12.8,       0.,        7.37  /)
      atom( 91: 96)%SEA = (/ 200.6,      7.57,    175.9,    558.,        75.3,   16.2/)
      atom(200)%SEA = 0.
      atom(201)%SEA = .00052


     return

  end subroutine definition_data_neutrons
end module Neutrons_data
