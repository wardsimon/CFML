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
!!----
!!---- MODULE: CFML_GlobalDeps (Windows version)
!!----   INFO: Precision for CrysFML library and Operating System information
!!----         All the global variables defined in this module are implicitly public.
!!----
!!----
!!
Module CFML_DefPar
   !---- Use Modules ----!
   use CFML_GlobalDeps

   !---- Local Variables ----!
   implicit none

   !--------------------!
   !---- PARAMETERS ----!
   !--------------------!
   integer, parameter :: AP_SPECIES_N  = 183   ! Number of atomic properties species in:  Ap_Table
   integer, parameter :: BVEL_ANIONS_N =   1   ! Number of anions known in BVEL Table in: Stefan Adams and R. Prasada Rao
   integer, parameter :: BVEL_SPECIES_N= 132   ! Maximum Number of species in BVEL_Table
   integer, parameter :: BVS_ANIONS_N  =  14   ! Number of anions known in BVS Table in O"Keefe, Breese, Brown
   integer, parameter :: BVS_SPECIES_N = 247   ! Maximum Number of species in BVS_Table
   integer, parameter :: MAX_FREE_PAR  =3000   ! Maximum number of free parameters
   integer, parameter :: NUM_CHEM_INFO = 108   ! Number of total Chem_info Data
   integer, parameter :: NUM_DELTA_FP  =  98   ! Number of total Delta (Fp,Fpp) Data
   integer, parameter :: NUM_MAG_FORM  = 119   ! Number of total Magnetic_Form Data
   integer, parameter :: NUM_MAG_J2    =  97   ! Number of <j2> Magnetic_Form Data
   integer, parameter :: NUM_MAG_J4    =  97   ! Number of <j4> Magnetic_Form Data
   integer, parameter :: NUM_MAG_J6    =  39   ! Number of <j6> Magnetic_Form Data
   integer, parameter :: NUM_XRAY_FORM = 214   ! Number of total Xray_Form Data
   integer, parameter :: SBVS_SPECIES_N= 168   ! Maximum Number of species in SBVS_Table

   integer, parameter, dimension(36,3,3)  :: MOD6 = reshape (  (/              &  ! Matrix Types For Rotational Operators In Conventional Basis
       1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,                   &  ! 1->24 Oh, 25->36 D6h
      -1, 1, 0, 0, 0, 0, 1, 0,-1,-1, 0, 1, 0, 1,-1, 0,-1, 1,                   &
       0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 1,-1, 0,-1, 1, 1, 0,-1,-1, 0, 1,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0,                   &
       0, 0,-1, 1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 1,-1, 1,-1, 0,-1, 1, 0,                   &
       1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 1,-1, 1,-1, 1,-1, 0,-1, 1, 0, 0,-1, 1, 0, 1,-1,                   &
       0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1,                   &
      -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1, 1,                   &
      -1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1, 1, 1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1 /), (/36,3,3/) )

   integer, parameter, dimension(1000) :: PRIMES =    & ! List of the first 1000 prime numbers
           (/ 2,      3,      5,      7,     11,     13,     17,     19,     23,     29,  &
             31,     37,     41,     43,     47,     53,     59,     61,     67,     71,  &
             73,     79,     83,     89,     97,    101,    103,    107,    109,    113,  &
            127,    131,    137,    139,    149,    151,    157,    163,    167,    173,  &
            179,    181,    191,    193,    197,    199,    211,    223,    227,    229,  &
            233,    239,    241,    251,    257,    263,    269,    271,    277,    281,  &
            283,    293,    307,    311,    313,    317,    331,    337,    347,    349,  &
            353,    359,    367,    373,    379,    383,    389,    397,    401,    409,  &
            419,    421,    431,    433,    439,    443,    449,    457,    461,    463,  &
            467,    479,    487,    491,    499,    503,    509,    521,    523,    541,  &
            547,    557,    563,    569,    571,    577,    587,    593,    599,    601,  &
            607,    613,    617,    619,    631,    641,    643,    647,    653,    659,  &
            661,    673,    677,    683,    691,    701,    709,    719,    727,    733,  &
            739,    743,    751,    757,    761,    769,    773,    787,    797,    809,  &
            811,    821,    823,    827,    829,    839,    853,    857,    859,    863,  &
            877,    881,    883,    887,    907,    911,    919,    929,    937,    941,  &
            947,    953,    967,    971,    977,    983,    991,    997,   1009,   1013,  &
           1019,   1021,   1031,   1033,   1039,   1049,   1051,   1061,   1063,   1069,  &
           1087,   1091,   1093,   1097,   1103,   1109,   1117,   1123,   1129,   1151,  &
           1153,   1163,   1171,   1181,   1187,   1193,   1201,   1213,   1217,   1223,  &
           1229,   1231,   1237,   1249,   1259,   1277,   1279,   1283,   1289,   1291,  &
           1297,   1301,   1303,   1307,   1319,   1321,   1327,   1361,   1367,   1373,  &
           1381,   1399,   1409,   1423,   1427,   1429,   1433,   1439,   1447,   1451,  &
           1453,   1459,   1471,   1481,   1483,   1487,   1489,   1493,   1499,   1511,  &
           1523,   1531,   1543,   1549,   1553,   1559,   1567,   1571,   1579,   1583,  &
           1597,   1601,   1607,   1609,   1613,   1619,   1621,   1627,   1637,   1657,  &
           1663,   1667,   1669,   1693,   1697,   1699,   1709,   1721,   1723,   1733,  &
           1741,   1747,   1753,   1759,   1777,   1783,   1787,   1789,   1801,   1811,  &
           1823,   1831,   1847,   1861,   1867,   1871,   1873,   1877,   1879,   1889,  &
           1901,   1907,   1913,   1931,   1933,   1949,   1951,   1973,   1979,   1987,  &
           1993,   1997,   1999,   2003,   2011,   2017,   2027,   2029,   2039,   2053,  &
           2063,   2069,   2081,   2083,   2087,   2089,   2099,   2111,   2113,   2129,  &
           2131,   2137,   2141,   2143,   2153,   2161,   2179,   2203,   2207,   2213,  &
           2221,   2237,   2239,   2243,   2251,   2267,   2269,   2273,   2281,   2287,  &
           2293,   2297,   2309,   2311,   2333,   2339,   2341,   2347,   2351,   2357,  &
           2371,   2377,   2381,   2383,   2389,   2393,   2399,   2411,   2417,   2423,  &
           2437,   2441,   2447,   2459,   2467,   2473,   2477,   2503,   2521,   2531,  &
           2539,   2543,   2549,   2551,   2557,   2579,   2591,   2593,   2609,   2617,  &
           2621,   2633,   2647,   2657,   2659,   2663,   2671,   2677,   2683,   2687,  &
           2689,   2693,   2699,   2707,   2711,   2713,   2719,   2729,   2731,   2741,  &
           2749,   2753,   2767,   2777,   2789,   2791,   2797,   2801,   2803,   2819,  &
           2833,   2837,   2843,   2851,   2857,   2861,   2879,   2887,   2897,   2903,  &
           2909,   2917,   2927,   2939,   2953,   2957,   2963,   2969,   2971,   2999,  &
           3001,   3011,   3019,   3023,   3037,   3041,   3049,   3061,   3067,   3079,  &
           3083,   3089,   3109,   3119,   3121,   3137,   3163,   3167,   3169,   3181,  &
           3187,   3191,   3203,   3209,   3217,   3221,   3229,   3251,   3253,   3257,  &
           3259,   3271,   3299,   3301,   3307,   3313,   3319,   3323,   3329,   3331,  &
           3343,   3347,   3359,   3361,   3371,   3373,   3389,   3391,   3407,   3413,  &
           3433,   3449,   3457,   3461,   3463,   3467,   3469,   3491,   3499,   3511,  &
           3517,   3527,   3529,   3533,   3539,   3541,   3547,   3557,   3559,   3571,  &
           3581,   3583,   3593,   3607,   3613,   3617,   3623,   3631,   3637,   3643,  &
           3659,   3671,   3673,   3677,   3691,   3697,   3701,   3709,   3719,   3727,  &
           3733,   3739,   3761,   3767,   3769,   3779,   3793,   3797,   3803,   3821,  &
           3823,   3833,   3847,   3851,   3853,   3863,   3877,   3881,   3889,   3907,  &
           3911,   3917,   3919,   3923,   3929,   3931,   3943,   3947,   3967,   3989,  &
           4001,   4003,   4007,   4013,   4019,   4021,   4027,   4049,   4051,   4057,  &
           4073,   4079,   4091,   4093,   4099,   4111,   4127,   4129,   4133,   4139,  &
           4153,   4157,   4159,   4177,   4201,   4211,   4217,   4219,   4229,   4231,  &
           4241,   4243,   4253,   4259,   4261,   4271,   4273,   4283,   4289,   4297,  &
           4327,   4337,   4339,   4349,   4357,   4363,   4373,   4391,   4397,   4409,  &
           4421,   4423,   4441,   4447,   4451,   4457,   4463,   4481,   4483,   4493,  &
           4507,   4513,   4517,   4519,   4523,   4547,   4549,   4561,   4567,   4583,  &
           4591,   4597,   4603,   4621,   4637,   4639,   4643,   4649,   4651,   4657,  &
           4663,   4673,   4679,   4691,   4703,   4721,   4723,   4729,   4733,   4751,  &
           4759,   4783,   4787,   4789,   4793,   4799,   4801,   4813,   4817,   4831,  &
           4861,   4871,   4877,   4889,   4903,   4909,   4919,   4931,   4933,   4937,  &
           4943,   4951,   4957,   4967,   4969,   4973,   4987,   4993,   4999,   5003,  &
           5009,   5011,   5021,   5023,   5039,   5051,   5059,   5077,   5081,   5087,  &
           5099,   5101,   5107,   5113,   5119,   5147,   5153,   5167,   5171,   5179,  &
           5189,   5197,   5209,   5227,   5231,   5233,   5237,   5261,   5273,   5279,  &
           5281,   5297,   5303,   5309,   5323,   5333,   5347,   5351,   5381,   5387,  &
           5393,   5399,   5407,   5413,   5417,   5419,   5431,   5437,   5441,   5443,  &
           5449,   5471,   5477,   5479,   5483,   5501,   5503,   5507,   5519,   5521,  &
           5527,   5531,   5557,   5563,   5569,   5573,   5581,   5591,   5623,   5639,  &
           5641,   5647,   5651,   5653,   5657,   5659,   5669,   5683,   5689,   5693,  &
           5701,   5711,   5717,   5737,   5741,   5743,   5749,   5779,   5783,   5791,  &
           5801,   5807,   5813,   5821,   5827,   5839,   5843,   5849,   5851,   5857,  &
           5861,   5867,   5869,   5879,   5881,   5897,   5903,   5923,   5927,   5939,  &
           5953,   5981,   5987,   6007,   6011,   6029,   6037,   6043,   6047,   6053,  &
           6067,   6073,   6079,   6089,   6091,   6101,   6113,   6121,   6131,   6133,  &
           6143,   6151,   6163,   6173,   6197,   6199,   6203,   6211,   6217,   6221,  &
           6229,   6247,   6257,   6263,   6269,   6271,   6277,   6287,   6299,   6301,  &
           6311,   6317,   6323,   6329,   6337,   6343,   6353,   6359,   6361,   6367,  &
           6373,   6379,   6389,   6397,   6421,   6427,   6449,   6451,   6469,   6473,  &
           6481,   6491,   6521,   6529,   6547,   6551,   6553,   6563,   6569,   6571,  &
           6577,   6581,   6599,   6607,   6619,   6637,   6653,   6659,   6661,   6673,  &
           6679,   6689,   6691,   6701,   6703,   6709,   6719,   6733,   6737,   6761,  &
           6763,   6779,   6781,   6791,   6793,   6803,   6823,   6827,   6829,   6833,  &
           6841,   6857,   6863,   6869,   6871,   6883,   6899,   6907,   6911,   6917,  &
           6947,   6949,   6959,   6961,   6967,   6971,   6977,   6983,   6991,   6997,  &
           7001,   7013,   7019,   7027,   7039,   7043,   7057,   7069,   7079,   7103,  &
           7109,   7121,   7127,   7129,   7151,   7159,   7177,   7187,   7193,   7207,  &
           7211,   7213,   7219,   7229,   7237,   7243,   7247,   7253,   7283,   7297,  &
           7307,   7309,   7321,   7331,   7333,   7349,   7351,   7369,   7393,   7411,  &
           7417,   7433,   7451,   7457,   7459,   7477,   7481,   7487,   7489,   7499,  &
           7507,   7517,   7523,   7529,   7537,   7541,   7547,   7549,   7559,   7561,  &
           7573,   7577,   7583,   7589,   7591,   7603,   7607,   7621,   7639,   7643,  &
           7649,   7669,   7673,   7681,   7687,   7691,   7699,   7703,   7717,   7723,  &
           7727,   7741,   7753,   7757,   7759,   7789,   7793,   7817,   7823,   7829,  &
           7841,   7853,   7867,   7873,   7877,   7879,   7883,   7901,   7907,   7919 /)

   character(len=*), parameter, dimension(24)            :: BC_D6h =(/                &      ! Bradley & Cracknell Notation
             "  E  "," C+_3"," C-_3"," C_2 "," C-_6"," C+_6","C'_23","C'_21","C'_22", &
             "C`_23","C`_21","C`_22","  I  "," S-_6"," S+_6"," s_h "," S+_3"," S-_3", &
             " s_v3"," s_v1"," s_v2"," s_d3"," s_d1"," s_d2" /)
   character(len=*), parameter, dimension(48)            :: BC_Oh =(/                 &      ! Bradley & Cracknell Notation
             "  E  "," C_2z"," C_2y"," C_2x","C+_31","C+_34","C+_33","C+_32","C-_31", &
             "C-_33","C-_32","C-_34"," C_2a"," C_2b","C-_4z","C+_4z","C-_4x"," C_2d", &
             " C_2f","C+_4x","C+_4y"," C_2c","C-_4y"," C_2e","  I  "," s_z "," s_y ", &
             " s_x ","S-_61","S-_64","S-_63","S-_62","S+_61","S+_63","S+_62","S+_64", &
             " s_da"," s_db","S+_4z","S-_4z","S+_4x"," s_dd"," s_df","S-_4x","S-_4y", &
             " s_dc","S+_4y"," s_de"  /)
   character(len=*), parameter, dimension(BVEL_ANIONS_N) :: BVEL_ANIONS = (/"O-2 "/)          ! Anions known from Stefan Adams and R. Prasada Rao
   character(len=*), parameter, dimension(BVS_ANIONS_N)  :: BVS_ANIONS  = (/"O-2 ","F-1 ", &  ! Anions known from O'Keefe, Bresse, Brown
             "CL-1","BR-1","I-1 ","S-2 ","SE-2","TE-2","N-3 ","P-3 ","AS-3","H-1 ","O-1 ", &
             "SE-1"/)
   character(len=*), parameter, dimension(72)            :: DEPMAT = (/           &           ! Magnetic array
             "( Dx, Dy, Dz)      ","(-Dx,-Dy, Dz)      ","(-Dx, Dy,-Dz)      ",   &
             "( Dx,-Dy,-Dz)      ","( Dz, Dx, Dy)      ","( Dz,-Dx,-Dy)      ",   &
             "(-Dz,-Dx, Dy)      ","(-Dz, Dx,-Dy)      ","( Dy, Dz, Dx)      ",   &
             "(-Dy, Dz,-Dx)      ","( Dy,-Dz,-Dx)      ","(-Dy,-Dz, Dx)      ",   &
             "( Dy, Dx,-Dz)      ","(-Dy,-Dx,-Dz)      ","( Dy,-Dx, Dz)      ",   &
             "(-Dy, Dx, Dz)      ","( Dx, Dz,-Dy)      ","(-Dx, Dz, Dy)      ",   &
             "(-Dx,-Dz,-Dy)      ","( Dx,-Dz, Dy)      ","( Dz, Dy,-Dx)      ",   &
             "( Dz,-Dy, Dx)      ","(-Dz, Dy, Dx)      ","(-Dz,-Dy,-Dx)      ",   &
             "(-Dx,-Dy,-Dz)      ","( Dx, Dy,-Dz)      ","( Dx,-Dy, Dz)      ",   &
             "(-Dx, Dy, Dz)      ","(-Dz,-Dx,-Dy)      ","(-Dz, Dx, Dy)      ",   &
             "( Dz, Dx,-Dy)      ","( Dz,-Dx, Dy)      ","(-Dy,-Dz,-Dx)      ",   &
             "( Dy,-Dz, Dx)      ","(-Dy, Dz, Dx)      ","( Dy, Dz,-Dx)      ",   &
             "(-Dy,-Dx, Dz)      ","( Dy, Dx, Dz)      ","(-Dy, Dx,-Dz)      ",   &
             "( Dy,-Dx,-Dz)      ","(-Dx,-Dz, Dy)      ","( Dx,-Dz,-Dy)      ",   &
             "( Dx, Dz, Dy)      ","(-Dx, Dz,-Dy)      ","(-Dz,-Dy, Dx)      ",   &
             "(-Dz, Dy,-Dx)      ","( Dz,-Dy,-Dx)      ","( Dz, Dy, Dx)      ",   &
             "( Dx   ,    Dy, Dz)","(   -Dy, Dx-Dy, Dz)","(-Dx+Dy,-Dx   , Dz)",   &
             "(-Dx   ,   -Dy, Dz)","(    Dy,-Dx+Dy, Dz)","( Dx-Dy, Dx   , Dz)",   &
             "(    Dy, Dx   ,-Dz)","( Dx-Dy,   -Dy,-Dz)","(-Dx   ,-Dx+Dy,-Dz)",   &
             "(   -Dy,-Dx   ,-Dz)","(-Dx+Dy,    Dy,-Dz)","( Dx   , Dx-Dy,-Dz)",   &
             "(-Dx   ,   -Dy,-Dz)","(    Dy,-Dx+Dy,-Dz)","( Dx-Dy, Dx   ,-Dz)",   &
             "( Dx   ,    Dy,-Dz)","(   -Dy, Dx-Dy,-Dz)","(-Dx+Dy,-Dx   ,-Dz)",   &
             "(   -Dy,-Dx   , Dz)","(-Dx+Dy,    Dy, Dz)","( Dx   , Dx-Dy, Dz)",   &
             "(    Dy, Dx   , Dz)","( Dx-Dy,   -Dy, Dz)","(-Dx   ,-Dx+Dy, Dz)"   /)
   character(len=* ), parameter, dimension(24)           :: INTSYMD6H =(/          &     ! International Symbols For Point Group Elements Of 6/mmm (D6h)
             "  1           "," 3+ ( 0, 0, z)"," 3- ( 0, 0, z)","  2 ( 0, 0, z)",  &
             " 6- ( 0, 0, z)"," 6+ ( 0, 0, z)","  2 ( x, x, 0)","  2 ( x, 0, 0)",  &
             "  2 ( 0, y, 0)","  2 ( x,-x, 0)","  2 ( x,2x, 0)","  2 (2x, x, 0)",  &
             " -1           ","-3+ ( 0, 0, z)","-3- ( 0, 0, z)","  m ( x, y, 0)",  &
             "-6- ( 0, 0, z)","-6+ ( 0, 0, z)","  m ( x,-x, z)","  m ( x,2x, z)",  &
             "  m (2x, x, z)","  m ( x, x, z)","  m ( x, 0, z)","  m ( 0, y, z)"   /)
   character(len=* ), parameter, dimension(48)           :: INTSYMOH = (/          &
             "  1           ","  2 ( 0, 0, z)","  2 ( 0, y, 0)","  2 ( x, 0, 0)",  &     ! International Symbols For Point Group Elements Of M3M (Oh)
             " 3+ ( x, x, x)"," 3+ (-x, x,-x)"," 3+ ( x,-x,-x)"," 3+ (-x,-x, x)",  &
             " 3- ( x, x, x)"," 3- ( x,-x,-x)"," 3- (-x,-x, x)"," 3- (-x, x,-x)",  &
             "  2 ( x, x, 0)","  2 ( x,-x, 0)"," 4- ( 0, 0, z)"," 4+ ( 0, 0, z)",  &
             " 4- ( x, 0, 0)","  2 ( 0, y, y)","  2 ( 0, y,-y)"," 4+ ( x, 0, 0)",  &
             " 4+ ( 0, y, 0)","  2 ( x, 0, x)"," 4- ( 0, y, 0)","  2 (-x, 0, x)",  &
             " -1           ","  m ( x, y, 0)","  m ( x, 0, z)","  m ( 0, y, z)",  &
             "-3+ ( x, x, x)","-3+ (-x, x,-x)","-3+ ( x,-x,-x)","-3+ (-x,-x, x)",  &
             "-3- ( x, x, x)","-3- ( x,-x,-x)","-3- (-x,-x, x)","-3- (-x, x,-x)",  &
             "  m ( x,-x, z)","  m ( x, x, z)","-4- ( 0, 0, z)","-4+ ( 0, 0, z)",  &
             "-4- ( x, 0, 0)","  m ( x, y,-y)","  m ( x, y, y)","-4+ ( x, 0, 0)",  &
             "-4+ ( 0, y, 0)","  m (-x, y, x)","-4- ( 0, y, 0)","  m ( x, y, x)"   /)
   character(len=*), parameter, dimension(24)            :: KOV_D6H=(/           &       ! Kovalev Notation
             " h1"," h3"," h5"," h4"," h6"," h2","h11"," h9"," h7"," h8","h12",  &
             "h10","h13","h15","h17","h16","h18","h14","h23",                    &
             "h21","h19","h20","h24","h22"/)
   character(len=*), parameter, dimension(48)            :: KOV_OH=(/                   &  ! Kovalev Notation
             " h1"," h4"," h3"," h2"," h9","h10","h12","h11"," h5"," h7"," h6"," h8",   &
             "h16","h13","h15","h14","h20","h18","h17","h19","h24","h23",               &
             "h22","h21","h25","h28","h27","h26","h33","h34","h36","h35",               &
             "h29","h31","h30","h32","h40","h37","h39","h38","h44","h42",               &
             "h41","h43","h48","h47","h46","h45"/)
   character(len=* ), parameter, dimension( 8)           :: LATT =(/         &   ! Lattice Traslations
             "  P: { 000 }                                       ",          &
             "  A: { 000;  0  1/2 1/2 }+                         ",          &
             "  B: { 000; 1/2  0  1/2 }+                         ",          &
             "  C: { 000; 1/2 1/2  0  }+                         ",          &
             "  I: { 000; 1/2 1/2 1/2 }+                         ",          &
             "  R: { 000; 2/3 1/3 1/3; 1/3 2/3 2/3   }+          ",          &
             "  F: { 000;  0  1/2 1/2; 1/2  0  1/2; 1/2 1/2  0 }+",          &
             "  Z: { 000;  Unconventional Z-centering vectors  }+"   /)
   character(len=*), parameter, dimension(16)            :: LAUE_CLASS=(/     &   ! Laue Class
             "-1   ","2/m  ","mmm  ","4/m  ","4/mmm","-3 R ","-3m R","-3   ", &
             "-3m1 ","-31m ","6/m  ","6/mmm","m-3  ","m-3m ","m3   ","m3m  "/)
   character(len=*), parameter, dimension(48)            :: LITVIN_POINT_OP_LABEL=(/                  &
             "1       ","2x      ","2y      ","2z      ","3xyz-1  ","3xy-z   ","3-xyz   ","3x-yz   ", &
             "3xyz    ","3x-yz-1 ","3xy-z-1 ","3-xyz-1 ","2-xy    ","4z      ","4z-1    ","2xy     ", &
             "2-yz    ","2yz     ","4x      ","4x-1    ","2-xz    ","4y-1    ","2xz     ","4y      ", &
             "-1      ","mx      ","my      ","mz      ","-3xyz-1 ","-3xy-z  ","-3-xyz  ","-3x-yz  ", &
             "-3xyz   ","-3x-yz-1","-3xy-z-1","-3-xyz-1","m-xy    ","-4z     ","-4z-1   ","mxy     ", &
             "m-yz    ","myz     ","-4x     ","-4x-1   ","m-xz    ","-4y-1   ","mxz     ","-4y     "/)
   character(len=*), parameter, dimension(48)            :: LITVIN_POINT_OP=(/     &
             "x,y,z   ", "x,-y,-z ", "-x,y,-z ", "-x,-y,z ", "y,z,x   ",           &
             "y,-z,-x ", "-y,z,-x ", "-y,-z,x ", "z,x,y   ", "z,-x,-y ",           &
             "-z,x,-y ", "-z,-x,y ", "-y,-x,-z", "-y,x,z  ", "y,-x,z  ",           &
             "y,x,-z  ", "-x,-z,-y", "-x,z,y  ", "x,-z,y  ", "x,z,-y  ",           &
             "-z,-y,-x", "-z,y,x  ", "z,-y,x  ", "z,y,-x  ", "-x,-y,-z",           &
             "-x,y,z  ", "x,-y,z  ", "x,y,-z  ", "-y,-z,-x", "-y,z,x  ",           &
             "y,-z,x  ", "y,z,-x  ", "-z,-x,-y", "-z,x,y  ", "z,-x,y  ",           &
             "z,x,-y  ", "y,x,z   ", "y,-x,-z ", "-y,x,-z ", "-y,-x,z ",           &
             "x,z,y   ", "x,-z,-y ", "-x,z,-y ", "-x,-z,y ", "z,y,x   ",           &
             "z,-y,-x ", "-z,y,-x ", "-z,-y,x "/)
   character(len=*), parameter, dimension(24)            :: LITVIN_POINT_OP_HEX_LABEL=(/     &
             "1    ","6z   ","3z   ","2z   ","3z-1 ","6z-1 ","2x   ","21   ",                &
             "2xy  ","22   ","2y   ","23   ","-1   ","-6z  ","-3z  ","mz   ",                &
             "-3z-1","-6z-1","mx   ","m1   ","mxy  ","m2   ","my   ","m3   "/)
   character(len=*), parameter, dimension(24)            :: LITVIN_POINT_OP_HEX=(/          &
             "x,y,z     ","x-y,x,z   ","-y,x-y,z  ","-x,-y,z   ","-x+y,-x,z ","y,-x+y,z  ", &
             "x-y,-y,-z ","x,x-y,-z  ","y,x,-z    ","-x+y,y,-z ","-x,-x+y,-z","-y,-x,-z  ", &
             "-x,-y,-z  ","-x+y,-x,-z","y,-x+y,-z ","x,y,-z    ","x-y,x,-z  ","-y,x-y,-z ", &
             "-x+y,y,z  ","-x,-x+y,z ","-y,-x,z   ","x-y,-y,z  ","x,x-y,z   ","y,x,z     "/)
   character(len=*), parameter, dimension(72)            :: MAGMAT = (/           &
             "( Mx, My, Mz)      ","(-Mx,-My, Mz)      ","(-Mx, My,-Mz)      ",   &
             "( Mx,-My,-Mz)      ","( Mz, Mx, My)      ","( Mz,-Mx,-My)      ",   &
             "(-Mz,-Mx, My)      ","(-Mz, Mx,-My)      ","( My, Mz, Mx)      ",   &
             "(-My, Mz,-Mx)      ","( My,-Mz,-Mx)      ","(-My,-Mz, Mx)      ",   &
             "( My, Mx,-Mz)      ","(-My,-Mx,-Mz)      ","( My,-Mx, Mz)      ",   &
             "(-My, Mx, Mz)      ","( Mx, Mz,-My)      ","(-Mx, Mz, My)      ",   &
             "(-Mx,-Mz,-My)      ","( Mx,-Mz, My)      ","( Mz, My,-Mx)      ",   &
             "( Mz,-My, Mx)      ","(-Mz, My, Mx)      ","(-Mz,-My,-Mx)      ",   &
             "(-Mx,-My,-Mz)      ","( Mx, My,-Mz)      ","( Mx,-My, Mz)      ",   &
             "(-Mx, My, Mz)      ","(-Mz,-Mx,-My)      ","(-Mz, Mx, My)      ",   &
             "( Mz, Mx,-My)      ","( Mz,-Mx, My)      ","(-My,-Mz,-Mx)      ",   &
             "( My,-Mz, Mx)      ","(-My, Mz, Mx)      ","( My, Mz,-Mx)      ",   &
             "(-My,-Mx, Mz)      ","( My, Mx, Mz)      ","(-My, Mx,-Mz)      ",   &
             "( My,-Mx,-Mz)      ","(-Mx,-Mz, My)      ","( Mx,-Mz,-My)      ",   &
             "( Mx, Mz, My)      ","(-Mx, Mz,-My)      ","(-Mz,-My, Mx)      ",   &
             "(-Mz, My,-Mx)      ","( Mz,-My,-Mx)      ","( Mz, My, Mx)      ",   &
             "( Mx   ,    My, Mz)","(   -My, Mx-My, Mz)","(-Mx+My,-Mx   , Mz)",   &
             "(-Mx   ,   -My, Mz)","(    My,-Mx+My, Mz)","( Mx-My, Mx   , Mz)",   &
             "(    My, Mx   ,-Mz)","( Mx-My,   -My,-Mz)","(-Mx   ,-Mx+My,-Mz)",   &
             "(   -My,-Mx   ,-Mz)","(-Mx+My,    My,-Mz)","( Mx   , Mx-My,-Mz)",   &
             "(-Mx   ,   -My,-Mz)","(    My,-Mx+My,-Mz)","( Mx-My, Mx   ,-Mz)",   &
             "( Mx   ,    My,-Mz)","(   -My, Mx-My,-Mz)","(-Mx+My,-Mx   ,-Mz)",   &
             "(   -My,-Mx   , Mz)","(-Mx+My,    My, Mz)","( Mx   , Mx-My, Mz)",   &
             "(    My, Mx   , Mz)","( Mx-My,   -My, Mz)","(-Mx   ,-Mx+My, Mz)"   /)
   character(len=*), parameter, dimension(24)            :: ML_D6H=(/                   &   ! Miller & Love Notation
             " 1"," 3"," 5"," 4"," 6"," 2"," 9"," 7","11","12","10"," 8","13","15","17",&
             "16","18","14","21","19","23","24","22","20"/)
   character(len=*), parameter, dimension(48)            :: ML_Oh=(/                    &   ! Miller & Love Notation
             " 1"," 4"," 3"," 2"," 9","10","12","11"," 5"," 7"," 6"," 8","16","13","15",&
             "14","20","18","17","19","24","23","22","21","25","28","27","26","33","34",&
             "36","35","29","31","30","32","40","37","39","38","44","42","41","43","48",&
             "47","46","45"/)
   character(len=*), parameter, dimension(41)            :: POINT_GROUP=(/      &           ! Point Group Symbols
             "1    ","-1   ","2    ","m    ","2/m  ","222  ","mm2  ","m2m  ",   &
             "2mm  ","mmm  ","4    ","-4   ","4/m  ","422  ","4mm  ","-42m ",   &
             "-4m2 ","4/mmm","3    ","-3   ","32   ","3m   ","-3m  ","312  ",   &
             "31m  ","-31m ","6    ","-6   ","6/m  ","622  ","6mm  ","-62m ",   &
             "-6m2 ","6/mmm","23   ","m-3  ","432  ","-43m ","m-3m ","m3   ",   &
             "m3m  "/)
   character(len=*), parameter, dimension(0:35)          :: REFERENCES  = (/                     &  ! List of Reference for BVS Data
             "Unknown                                                                         ", &  !0
             "Brown and Altermatt, (1985), Acta Cryst. B41, 244-247 (empirical)               ", &  !1
             "Brese and O'Keeffe, (1991), Acta Cryst. B47, 192-197 (extrapolated)             ", &  !2
             "Adams, 2001, Acta Cryst. B57, 278-287 (includes second neighbours)              ", &  !3
             "Hu et al. (1995) Inorg. Chim. Acta, 232, 161-165.                               ", &  !4
             "I.D.Brown Private communication                                                 ", &  !5
             "Brown et al. (1984) Inorg. Chem. 23, 4506-4508                                  ", &  !6
             "Palenik (1997) Inorg. Chem. 36 4888-4890                                        ", &  !7
             "Kanowitz and Palenik (1998) Inorg. Chem. 37 2086-2088                           ", &  !8
             "Wood and Palenik (1998) Inorg. Chem. 37 4149-4151                               ", &  !9
             "Liu and Thorp (1993) Inorg. Chem. 32 4102-4105                                  ", &  !10
             "Palenik (1997) Inorg. Chem. 36 3394-3397                                        ", &  !11
             "Shields, Raithby, Allen and Motherwell (1999) Acta Cryst.B56, 455-465           ", &  !12
             "Chen, Zhou and Hu (2002) Chinese Sci. Bul. 47, 978-980.                         ", &  !13
             "Kihlbourg (1963) Ark. Kemi 21 471; Schroeder 1975 Acta Cryst. B31, 2294         ", &  !14
             "Allmann (1975) Monatshefte Chem. 106, 779                                       ", &  !15
             "Zachariesen (1978) J.Less Common Metals 62, 1                                   ", &  !16
             "Krivovichev and Brown (2001) Z. Krist. 216, 245                                 ", &  !17
             "Burns, Ewing and Hawthorne (1997) Can. Miner. 35,1551-1570                      ", &  !18
             "Garcia-Rodriguez, et al. (2000) Acta Cryst. B56, 565-569                        ", &  !19
             "Mahapatra et al. (1996) J. Amer.Chem. Soc. 118, 11555                           ", &  !20
             "Wood and Palenik (1999) Inorg. Chem. 38, 1031-1034                              ", &  !21
             "Wood and Palenik (1999) Inorg. Chem. 38, 3926-3930                              ", &  !22
             "Wood, Abboud, Palenik and Palenik (2000) Inorg. Chem. 39, 2065-2068             ", &  !23
             "Tytko, Mehnike and Kurad (1999) Structure and Bonding 93, 1-66                  ", &  !24
             "Gundemann, et al.(1999) J. Phys. Chem. A 103, 4752-4754                         ", &  !25
             "Zocchi (2000) Solid State Sci. 2 383-387                                        ", &  !26
             "Jensen, Palenik and Tiekiak (2001) Polyhedron 20, 2137                          ", &  !27
             "Roulhac and Palenik (2002) Inorg. Chem. 42, 118-121                             ", &  !28
             "Holsa et al.(2002) J.Solid State Chem 165, 48-55                                ", &  !29
             "Trzesowska, Kruszynski & Bartezak (2004) Acta Cryst. B60, 174-178               ", &  !30
             "Locock & Burns (2004) Z.Krist. 219, 267-271                                     ", &  !31
             "J.Rodriguez-Carvajal, Private communication                                     ", &  !32
             "S. Adams and R. Prasada Rao, (2011) Phys. Status Solidi A 208, No. 8, 1746-1753 ", &  !33
             "S. Adams (2013),  Structure and Bonding (eds. Brown & Poeppelmeier) 158, 91-128 ", &  !34
             "Adams S, Moretsky O and Canadell E (2004) Solid State Ionics 168, 281-290       "/)   !35
   character(len=*), parameter, dimension(7)             :: SYS_CRY =(/      &
             "Triclinic   ","Monoclinic  ","Orthorhombic","Tetragonal  ",    &
             "Trigonal    ","Hexagonal   ","Cubic       " /)
   character(len=*), parameter, dimension(24)            :: X_D6H = (/             &
             "( x  ,   y, z)","(  -y, x-y, z)","(-x+y,-x  , z)","(-x  ,  -y, z)",  &
             "(   y,-x+y, z)","( x-y, x  , z)","(   y, x  ,-z)","( x-y,  -y,-z)",  &
             "(-x  ,-x+y,-z)","(  -y,-x  ,-z)","(-x+y,   y,-z)","( x  , x-y,-z)",  &
             "(-x  ,  -y,-z)","(   y,-x+y,-z)","( x-y, x  ,-z)","( x  ,   y,-z)",  &
             "(  -y, x-y,-z)","(-x+y,-x  ,-z)","(  -y,-x  , z)","(-x+y,   y, z)",  &
             "( x  , x-y, z)","(   y, x  , z)","( x-y,  -y, z)","(-x  ,-x+y, z)"   /)
   character(len=*), parameter, dimension(48)            :: X_OH = (/                       &
             "( x, y, z)","(-x,-y, z)","(-x, y,-z)","( x,-y,-z)","( z, x, y)","( z,-x,-y)", &
             "(-z,-x, y)","(-z, x,-y)","( y, z, x)","(-y, z,-x)","( y,-z,-x)","(-y,-z, x)", &
             "( y, x,-z)","(-y,-x,-z)","( y,-x, z)","(-y, x, z)","( x, z,-y)","(-x, z, y)", &
             "(-x,-z,-y)","( x,-z, y)","( z, y,-x)","( z,-y, x)","(-z, y, x)","(-z,-y,-x)", &
             "(-x,-y,-z)","( x, y,-z)","( x,-y, z)","(-x, y, z)","(-z,-x,-y)","(-z, x, y)", &
             "( z, x,-y)","( z,-x, y)","(-y,-z,-x)","( y,-z, x)","(-y, z, x)","( y, z,-x)", &
             "(-y,-x, z)","( y, x, z)","(-y, x,-z)","( y,-x,-z)","(-x,-z, y)","( x,-z,-y)", &
             "( x, z, y)","(-x, z,-y)","(-z,-y, x)","(-z, y,-x)","( z,-y,-x)","( z, y, x)"  /)
   character(len=*), parameter, dimension(24)            :: ZAK_D6H =(/               &         ! Zak Notation
             "   E   "," C(z)_3","C(2z)_3","  C_2  ","C(5z)_6"," C(z)_6","  U(xy)",   &
             "  U(x) ","  U(y) ","  U(3) ","  U(2) ","  U(1) ","   I   ","S(5z)_6",   &
             " S(z)_6","  s(z) "," S(z)_3","S(2z)_3"," s(xy) ","  s(x) ","  s(y) ",   &
             "  s(3) ","  s(2) ","  s(1) " /)
   character(len=*), parameter, dimension(48)            :: ZAK_OH =(/                &         ! Zak Notation
             "     E     ","    U(z)   ","    U(y)   ","    U(x)   ","  C(xyz)_3 ",   &
             " C(-xy-z)_3"," C(x-y-z)_3"," C(-x-yz)_3"," C(2xyz)_3 ","C(2x-y-z)_3",   &
             " C(2x-yz)_3","C(-2xy-z)_3","    U(xy)  ","   U(-xy)  ","   C(3z)_4 ",   &
             "   C(z)_4  ","   C(3x)_4 ","    U(yz)  ","   U(y-z)  ","   C(x)_4  ",   &
             "   C(y)_4  ","    U(xz)  ","   C(3y)_4 ","   U(x-z)  ","      I    ",   &
             "    s(z)   ","    s(y)   ","    s(x)   "," S(5xyz)_6 ","S(-5xy-z)_6",   &
             "S(5x-y-z)_6","S(-5x-yz)_6","  S(xyz)_6 "," S(x-y-z)_6"," S(-x-yz)_6",   &
             " S(-xy-z)_6","    s(xy)  ","   s(-xy)  ","   S(z)_4  ","  S(3z)_4  ",   &
             "   S(x)_4  ","    s(yz)  ","   s(y-z)  ","  S(3x)_4  ","  S(3y)_4  ",   &
             "    s(xz)  ","   S(y)_4  ","   s(x-z)  " /)


   real(kind=cp), parameter, dimension(BVEL_ANIONS_N) :: BVEL_ANIONS_RION = (/1.40/)            ! Radii Ionic for Anions in BVEL
   real(kind=cp), parameter, dimension(BVS_ANIONS_N)  :: BVS_ANIONS_RION =  (/1.40, 1.19, &     ! Ionic Radii for Anions
                  1.67, 1.95, 2.16, 1.84, 1.98, 2.21, 1.71, 2.12, 2.22, 2.08, 1.35, 1.80/)
   real(kind=cp), parameter, dimension(3,2)           :: LTR_A =reshape ( (/0.0,0.0,0.0, 0.0,0.5,0.5/), (/3,2/) )
   real(kind=cp), parameter, dimension(3,2)           :: LTR_B =reshape ( (/0.0,0.0,0.0, 0.5,0.0,0.5/), (/3,2/) )
   real(kind=cp), parameter, dimension(3,2)           :: LTR_C =reshape ( (/0.0,0.0,0.0, 0.5,0.5,0.0/), (/3,2/) )
   real(kind=cp), parameter, dimension(3,4)           :: LTR_F =reshape ( (/0.0,0.0,0.0, 0.0,0.5,0.5, &
                                                                            0.5,0.0,0.5, 0.5,0.5,0.0 /),(/3,4/) )
   real(kind=cp), parameter, dimension(3,2)           :: LTR_I =reshape ( (/0.0,0.0,0.0, 0.5,0.5,0.5/), (/3,2/) )
   real(kind=cp), parameter, dimension(3,3)           :: LTR_R =reshape ( (/0.0,0.0,0.0, 2.0/3.0,1.0/3.0,1.0/3.0, &
                                                                              1.0/3.0,2.0/3.0,2.0/3.0/),(/3,3/) )


   !---------------!
   !---- TYPES ----!
   !---------------!

   !!----
   !!---- TYPE, PUBLIC :: ANOMALOUS_SC_TYPE
   !!--..
   !!---- Update: February - 2005
   !!
   Type :: Anomalous_Sc_Type
      character(len= 2)            :: Symb = " "        ! Symbol of the Chemical species
      real(kind=cp), dimension(5)  :: Fp   = 0.0        ! Delta Fp
      real(kind=cp), dimension(5)  :: Fpp  = 0.0        ! Delta Fpp
   End Type Anomalous_Sc_Type

   !!----
   !!---- TYPE :: Atomic_Properties_Type
   !!--..
   !!----    Type Definition for single atomic properties
   !!----
   !!---- Created: January - 2015
   !!
   Type :: Atomic_Properties_Type
       integer          :: Z      = 0   ! Atomic number
       character(len=4) :: Symb   = " " ! Element with charge
       integer          :: oxs    = 0   ! Nominal oxidation state
       integer          :: dox    = 0   ! Default oxidation state
       real             :: Mass   = 0.0 ! Atomic mass in atomic units
       integer          :: n      = 0   ! Principal Quantum number (period)
       integer          :: g      = 0   ! Group in the periodic table
       integer          :: b      = 0   ! Block (s:0, p:1, d:2, f:3)
       real             :: Rc     = 0.0 ! Covalent radius
       real             :: sigma  = 0.0 ! Softness
   End Type Atomic_Properties_Type

   !!----
   !!---- TYPE :: BVEL_PAR_TYPE
   !!--..
   !!----    Type Definition for BVEL Parameters
   !!----
   !!---- Created: December - 2014
   !!
   Type :: Bvel_Par_Type
      character(len=5)                       :: symb    = " "    !Symbol of the cation
      real(kind=cp),dimension(bvel_anions_n) :: Avcoor  = 0.0    !Average cation coordination number
      real(kind=cp),dimension(bvel_anions_n) :: Rzero   = 0.0    !Modified Bond-Valence parameter R0
      real(kind=cp),dimension(bvel_anions_n) :: Rcutoff = 0.0    !Cutoff distance in Angstroms
      real(kind=cp),dimension(bvel_anions_n) :: Dzero   = 0.0    !First Morse potential parameter (eV)
      real(kind=cp),dimension(bvel_anions_n) :: Rmin    = 0.0    !Second Morse potential parameter (Angstroms)
      real(kind=cp),dimension(bvel_anions_n) :: alpha   = 0.0    !Third Morse potential parameter (1/b) (Angstroms^-1)
      integer      ,dimension(bvel_anions_n) :: refnum  = 0      !Pointer to reference paper
   End Type Bvel_Par_Type

   !!---- TYPE :: BVS_PAR_TYPE
   !!--..
   !!----    Definition for BVS Parameters
   !!----
   !!---- Update: February - 2005
   !!
   Type :: Bvs_Par_Type
      character (len=4)                     :: Symb    = " "      ! Chemical symbol
      real(kind=cp),dimension(bvs_anions_n) :: d0      = 0.0      ! D0 Parameter
      real(kind=cp),dimension(bvs_anions_n) :: b_par   = 0.0      ! B Parameter
      integer      ,dimension(bvs_anions_n) :: refnum  = 0        ! Integer pointing to the reference paper
   End Type Bvs_Par_Type

   !!----
   !!---- TYPE, PUBLIC :: CHEM_INFO_TYPE
   !!--..
   !!---- Update: February - 2005
   !!
   Type :: Chem_Info_Type
      character (len= 2)         :: Symb   = " "      ! Symbol of the Element
      character (len=12)         :: Name   = " "      ! Name of the Element
      integer                    :: Z      = 0        ! Atomic Number
      real(kind=cp)              :: AtWe   = 0.0      ! Atomic weight
      real(kind=cp)              :: RCov   = 0.0      ! Covalent Radius
      real(kind=cp)              :: RWaals = 0.0      ! van der Waals Radius
      real(kind=cp)              :: VAtm   = 0.0      ! Atomic volumen
      integer, dimension(5)      :: Oxid   = 0        ! Oxidation State
      real(kind=cp), dimension(5):: Rion   = 0.0      ! Ionic Radius (depending of the oxidation)
      real(kind=cp)              :: SctF   = 0.0      ! Fermi length [10**(-12) cm]
      real(kind=cp)              :: SedInc = 0.0      ! Incoherent Scattering Neutron cross-section (barns -> [10**(-24) cm**2] )
      real(kind=cp)              :: Sea    = 0.0      ! Neutron Absorption cross-section ( barns, for v= 2200m/s, l(A)=3.95/v (km/s) )
   End Type Chem_Info_Type

   !!----
   !!---- TYPE :: DERIV_TOF_TYPE
   !!--..
   !!---- Type Definition for TOF Profiles
   !!----
   !!---- Update: 11/07/2015
   !!
   Type :: Deriv_TOF_Type
      real(kind=cp) :: alfa  = 0.0    ! omega_a  DOmega/Dalpha
      real(kind=cp) :: beta  = 0.0    ! omega_b  DOmega/Dbeta
      real(kind=cp) :: dt    = 0.0    ! omega_t  DOmega/Ddt      (dt=TOFi-TOF(Bragg))
      real(kind=cp) :: sigma = 0.0    ! omega_s  DOmega/Dsigma   (for tof_Jorgensen function)
      real(kind=cp) :: gamma = 0.0    ! omega_g  DOmega/Dgamma   (for tof_Jorgensen_VonDreele function)
      real(kind=cp) :: eta   = 0.0    ! omega_e  DOmega/Deta                     "
      real(kind=cp) :: kappa = 0.0    ! omega_e  DOmega/kappa    (for tof_Carpenter function)
   End Type Deriv_TOF_Type

   !!----
   !!---- TYPE :: ERR_TEXT_TYPE
   !!--..
   !!---- Update: February - 2005
   !!
   Type :: Err_Text_Type
      integer :: nlines
      character (len=132), dimension(5) :: txt
   End Type Err_Text_Type

   !!----
   !!---- TYPE :: LSQ_CONDITIONS_TYPE
   !!--..
   !!----  Derived type encapsulating all necessary conditions for running the LSQ algorithm
   !!----
   !!
   Type :: LSQ_Conditions_Type
      logical          :: constr  =.false.   ! if true box constraint of percent% are applied to parameters
      logical          :: reached =.false.   ! if true convergence was reached in the algorithm
      integer          :: corrmax =50        ! value of correlation in % to output
      integer          :: nfev    =0         ! number of function evaluations (output component, useful for assessing LM algorithm)
      integer          :: njev    =0         ! number of Jacobian evaluations                 "
      integer          :: icyc    =0         ! number of cycles of refinement or maximum number of function evaluations in LM
                                             ! In LM procedures the default value is icyc = maxfev = 100(npvar+1)
      integer          :: npvar   =0         ! number of effective free parameters of the model
      integer          :: iw      =0         ! indicator for weighting scheme (if iw=1 => w=1/yc)
      integer          :: nprint  =0         ! indicator for printing during iterations, if nprint > 0 printing each nprint iterations
      real(kind=cp)    :: tol     =0.0       ! tolerance value for applying stopping criterion in LM algorithm
      real(kind=cp)    :: percent =0.0       ! %value of maximum variation of a parameter w.r.t.
                                             ! the intial value before fixing it
   End Type LSQ_Conditions_Type

   !!----
   !!---- TYPE :: LSQ_DATA_TYPE
   !!--..
   !!----
   !!----  Derived type encapsulating the observed and calculated data as well as the
   !!----  weighting factors, a variable related with each observed value and the
   !!----  total number of observations. It is responsibility of the calling program to
   !!----  allocate the components before calling the Marquardt_fit procedure.
   !!----
   !!---- Update: 11/07/2015
   !!
   Type :: LSQ_Data_Type
      integer                                 :: nobs=0  !total number of observations
      integer                                 :: iw  =0  !Indicator for type of values contained in component sw
      real(kind=cp), dimension(:),allocatable :: x       !Vector containing a relevant quantity for each observation (x-coordinate ...)
      real(kind=cp), dimension(:),allocatable :: y       !Vector containing the observed values
      real(kind=cp), dimension(:),allocatable :: sw      !Vector containing the standard deviation of observations if iw=0
                                                         !or the weight factors for least squares refinement if iw=1
      real(kind=cp), dimension(:),allocatable :: yc      !Vector containing the calculated values
   End Type LSQ_Data_type

   !!----
   !!---- TYPE :: LSQ_STATE_VECTOR_TYPE
   !!--..
   !!----  Derived type encapsulating the vector state defining a set of parameter
   !!----  for calculating the model function and running the LSQ algorithm.
   !!----  Now, with the introduction of code_comp and mul, the codes may be also interpreted
   !!----  as the ordinal number in the LSQ list of independent parameters. Depending on the
   !!----  the way the user program attributes codes and constraints a call to the subroutine
   !!----  Modify_Codes_State_Vector (see below)
   !!----
   !!---- Update: 11/07/2015
   !!
   Type :: LSQ_State_Vector_Type
      integer                                    :: np        = 0       !total number of model parameters <= Max_Free_Par
      logical                                    :: code_comp = .false. !If .true. the codes are interpreted as number in the LSQ list
      integer(kind=2)                            :: code_max  = 0       !Maximum code number (used in case of code_comp=.true.)
      real(kind=cp),     dimension(Max_Free_Par) :: mul       = 1.0     !Vector of multipliers (used in case of code_comp=.true.)
      real(kind=cp),     dimension(Max_Free_Par) :: pv        = 0.0     !Vector of parameters
      real(kind=cp),     dimension(Max_Free_Par) :: spv       = 0.0     !Vector of standard deviations
      real(kind=cp),     dimension(Max_Free_Par) :: dpv       = 0.0     !Vector of derivatives at a particular point
      integer(kind=2),   dimension(Max_Free_Par) :: code      = 0       !pointer for selecting variable parameters
      character(len=40), dimension(Max_Free_Par) :: nampar    =" "      !Names of parameters
   End Type LSQ_State_Vector_type

   !!----
    !!---- TYPE :: MAGNETIC_FORM_TYPE
    !!--..
    !!---- Update: February - 2005
    !!
    Type :: Magnetic_Form_Type
       character (len= 4)         :: Symb = " "         ! Symbol of the Chemical species
       real(kind=cp), dimension(7):: SctM = 0.0        ! Scattering Factors
    End Type Magnetic_Form_Type

   !!----
   !!---- TYPE :: OPT_CONDITIONS_TYPE
   !!----
   !!----    This TYPE has been introduced to simplify the call to optimization
   !!----    procedures. It contains the optimization parameters useful for different
   !!----    algorithms.
   !!----
   !!----    All integer components are initialized to zero and the real components
   !!----    are initilized as indicated below.
   !!----    A variable of this type should be defined by the user and all their
   !!----    input parameters (in) must be provided before calling the procedures.
   !!----    On output from the procedure the (out) items are provided for checking.
   !!----
   !!
   Type :: Opt_Conditions_Type
      character(len=20) :: method = " "    ! Name of the method
      integer           :: nmeth  = 0      ! 0= conjugate gradient, 1= BFGS method
      integer           :: npar   = 0      ! Number of free parameters
      integer           :: mxfun  = 0      ! Maximum number function calls
      integer           :: loops  = 0      ! Useful for SIMPLEX method: = 0
      integer           :: iquad  = 0      ! For SIMPLEX, if iquad/= 0 fitting to a quadratic
      integer           :: iout   = 0      ! =0 no printing for Quasi_Newton & Conjugate Gradient. Partial printing for Simplex (<0 no printing).
                                           ! > 0 printing each iout iterations/evaluations
      integer           :: nflag  = 0      ! Flag value states which condition caused the exit of the optimization subroutine
                                           !       If NFLAG=0, the algorithm has converged.
                                           !       If NFLAG=1, the maximum number of function
                                           !          evaluations have been used.
                                           !       If NFLAG=2, the linear search has failed to
                                           !          improve the function value. This is the
                                           !          usual exit if either the function or the
                                           !          gradient is incorrectly coded.
                                           !       If NFLAG=3, The search vector was not
                                           !          a descent direction. This can only be caused
                                           !          by roundoff,and may suggest that the
                                           !          convergence criterion is too strict.
      integer           :: ifun   = 0      ! Total number of function and gradient evaluations
      integer           :: iter   = 0      ! Total number of search directions used in the algorithm
      real(kind=cp)     :: eps    = 0.0    ! Convergence occurs when the norm of the gradient is less than or equal to EPS times the maximum
                                           ! of one and the norm of the vector X. Initialized to 1.0E-6
      real(kind=cp)     :: acc    = 0.0    ! ACC is a user supplied estimate of machine accuracy ACC=10.0E-20 has proved satisfactory
                                           ! For simplex method this should be changed to 1.0e-6
   End Type Opt_Conditions_Type

   !!----
   !!---- TYPE :: Points_Interval_Type
   !!--..
   !!---- Type used for FFT routines
   !!----
   !!
   Type :: Points_Interval_Type
      integer       :: Np   = 0              ! Number of Points
      real(kind=cp) :: Low  = 0.0            ! Lower range value
      real(kind=cp) :: High = 0.0            ! Upper range value
   End Type Points_Interval_Type

   !!---- TYPE :: sBVS_PAR_TYPE
   !!--..
   !!----    Definition for sBVS Parameters
   !!----
   !!---- Update: February - 2005
   !!
   Type :: sBvs_Par_Type
      character (len=4)                     :: Symb    = " "      ! Chemical symbol
      real(kind=cp),dimension(bvs_anions_n) :: d0      = 0.0      ! D0 Parameter
      real(kind=cp),dimension(bvs_anions_n) :: b_par   = 0.0      ! B Parameter
      real(kind=cp),dimension(bvs_anions_n) :: cn      = 0.0      ! Preferred Coordination
      real(kind=cp),dimension(bvs_anions_n) :: ctoff   = 0.0      ! Cutoff distance
      integer      ,dimension(bvs_anions_n) :: refnum  = 0        ! Integer pointing to the reference paper
   End Type sBvs_Par_Type

   !!----
   !!---- TYPE :: SPGR_INFO_TYPE
   !!--..
   !!----    Definition for General Info about Space Groups
   !!----
   !!---- Update: February - 2005
   !!
   Type :: Spgr_Info_Type
      integer                 :: N         = 0        ! Number of the Spacegroup
      character (len=12)      :: HM        = " "      ! Hermann-Mauguin
      character (len=16)      :: Hall      = " "      ! Hall
      integer                 :: Laue      = 0        ! Laue Group
      integer                 :: Pg        = 0        ! Point group
      integer, dimension(6)   :: Asu       = 0        ! Asymmetric unit * 24
      character (len= 5)      :: Inf_Extra = " "      ! Extra information
   End Type Spgr_Info_Type

   !!----
   !!---- TYPE :: TABLE_EQUIV_TYPE
   !!--..
   !!----    Definition for Equivalences on a Table
   !!----
   !!---- Update: February - 2005
   !!
   Type :: Table_Equiv_Type
      character(len= 6)      :: SC  = " "           ! Schoenflies
      character(len=17)      :: ML  = " "           ! Miller & Love
      character(len=18)      :: KO  = " "           ! Kovalev
      character(len=32)      :: BC  = " "           ! Bradley & Cracknell
      character(len=18)      :: ZA  = " "           ! Zak
   End Type Table_Equiv_Type

   !!----
    !!---- TYPE :: WYCK_INFO_TYPE
    !!--..
    !!----    Definition for Wyckoff Positions acording to IT
    !!----
    !!---- Update: February - 2005
    !!
    Type :: Wyck_Info_Type
       character (len=12)               :: HM     = " "    ! Hermann-Mauguin
       integer                          :: Norbit = 0      ! Number of orbites
       character (len=15),dimension(26) :: Corbit = " "    ! Generator of the orbit
    End Type Wyck_Info_Type

   !!----
   !!---- TYPE :: XRAY_FORM_TYPE
   !!--..
   !!---- Update: February - 2005
   !!
   Type :: Xray_Form_Type
      character (len= 4)         :: Symb = " "      ! Symbol of the Chemical species
      integer                    :: Z    = 0        ! Atomic Number
      real(kind=cp), dimension(4):: a    = 0.0      ! Coefficients for calculating the X-ray scattering factors
      real(kind=cp), dimension(4):: b    = 0.0      ! f(s) = Sum_{i=1,4} { a(i) exp(-b(i)*s^2) } + c
      real(kind=cp)              :: c    = 0.0      ! s=sinTheta/Lambda
   End Type Xray_Form_Type

   !!----
   !!---- TYPE :: XRAY_WAVELENGTH_TYPE
   !!--..
   !!---- Update: February - 2005
   !!
   Type :: Xray_Wavelength_Type
      character (len= 2)         :: Symb  = " "       ! Symbol of the Chemical species
      real(kind=cp), dimension(2):: Kalfa = 0.0       ! K-Serie for X-ray
      real(kind=cp)              :: Kbeta = 0.0       ! K-Serie for X-ray
   End Type Xray_Wavelength_Type


   !-------------------!
   !---- VARIABLES ----!
   !-------------------!

   logical            :: ERR_Bond    =.false.          ! Error flag in CFML_Bond_Tables module
   logical            :: ERR_MathGen =.false.          ! Error flag in CFML_Math_General module
   logical            :: ERR_Random  =.false.          ! Error flag in CFML_Random_Generators module
   logical            :: ERR_Spher   =.false.          ! Error flag in CFML_Spherical_Harmonics module
   logical            :: ERR_String  =.false.          ! Error flag in CFML_String_Utilities module
   logical            :: ERR_Math3D  =.false.          ! Error flag in CFML_Math_3D module
   logical            :: ERR_Optim   =.false.          ! Error flag in CFML_Optimization_General module
   logical            :: ERR_LSQ     =.false.          ! Error flag in CFML_Optimization_LSQ module
   logical            :: ERR_Symtab  =.false.          ! Error flag in CFML_Symmetry_Tables module

   logical            :: Init_ProfVal=.false.
   logical            :: Lorcomp     =.false.          ! .true. if there are Lorentzian components

   character(len=256) :: ERR_Bond_Mess    = " "        ! String containing information about the last error
   character(len=256) :: ERR_MathGen_Mess = " "        ! String containing information about the last error
   character(len=256) :: ERR_Random_Mess  = " "        ! String containing information about the last error
   character(len=256) :: ERR_Spher_Mess   = " "        ! String containing information about the last error
   character(len=256) :: ERR_String_Mess  = " "        ! String containing information about the last error
   character(len=256) :: ERR_Math3D_Mess  = " "        ! String containing information about the last error
   character(len=256) :: ERR_Optim_Mess   = " "        ! String containing information about the last error
   character(len=256) :: ERR_LSQ_Mess     = " "        ! String containing information about the last error
   character(len=256) :: ERR_SymTab_Mess  = " "        ! String containing information about the last error

   character(len=150) :: Info_Lsq_Mess    = " "        ! Information in Levenberg_Marquardt_Fit procedure

   integer            :: iErr_fmt    = 0               ! Integer signaling if an error has occurred (/=0) in using the procedure findFMT
   integer            :: win_console = -1              ! Code number for Scroll Window (Variable only in use for Winteracter code)

   integer,                      allocatable, dimension(:,:) :: Table_ref    ! Matrix N_Species x N_Species with references for BVS parameters

   real(kind=cp),                allocatable, dimension(:,:,:) :: Bond_Length_Table ! Variable to hold the Bonds length between type of atoms. Order by Z
   real(kind=cp),                allocatable, dimension(:,:) :: Table_Alpha  ! Matrix N_Species x N_Species of Alpha (equivalent to 1/b in BVS) parameters for BVEL
   real(kind=cp),                allocatable, dimension(:,:) :: Table_Avcoor ! Matrix N_Species x N_Species of Average coordination parameters for BVEL
   real(kind=cp),                allocatable, dimension(:,:) :: Table_b      ! Matrix N_Species x N_Species of B parameters for BVS
   real(kind=cp),                allocatable, dimension(:,:) :: Table_d0     ! Matrix N_Species x N_Species of D0 for BVS
   real(kind=cp),                allocatable, dimension(:,:) :: Table_Dzero  ! Matrix N_Species x N_Species of Dzero parameters for BVEL
   real(kind=cp),                allocatable, dimension(:,:) :: Table_Rcutoff! Matrix N_Species x N_Species of Rcutoff parameters for BVEL
   real(kind=cp),                allocatable, dimension(:,:) :: Table_Rmin   ! Matrix N_Species x N_Species of Rmin parameters for BVEL
   real(kind=cp),                allocatable, dimension(:,:) :: Table_Rzero  ! Matrix N_Species x N_Species of Rzero (equivalent to D0 in BVS) parameters for BVEL



   Type(Anomalous_Sc_Type),      allocatable, dimension(:)   :: Anomalous_ScFac ! Table of Delta-Fp and Delta-Fpp for 5 common radiations.
   Type(Atomic_Properties_Type), allocatable, dimension(:)   :: AP_Table
   Type(Bvel_Par_Type),          allocatable, dimension(:)   :: BVEL_Table
   Type(Bvs_Par_Type),           allocatable, dimension(:)   :: BVS_Table
   Type(Chem_Info_Type),         allocatable, dimension(:)   :: Chem_Info       ! Tabulated chemical data
   Type(Magnetic_Form_Type),     allocatable, dimension(:)   :: Magnetic_Form   ! Tabulated magnetic form factor data
   Type(Magnetic_Form_Type),     allocatable, dimension(:)   :: Magnetic_j2     ! Tabulated magnetic form factor J2
   Type(Magnetic_Form_Type),     allocatable, dimension(:)   :: Magnetic_j4     ! Tabulated magnetic form factor J4
   Type(Magnetic_Form_Type),     allocatable, dimension(:)   :: Magnetic_j6     ! Tabulated magnetic form factor J6
   Type(sBvs_Par_Type),          allocatable, dimension(:)   :: sBVS_Table
   Type(Xray_Form_Type),         allocatable, dimension(:)   :: Xray_Form       ! Tabulated Xray scattering factor coefficients
   Type(Spgr_Info_Type),         allocatable, dimension(:)   :: Spgr_Info       ! General Info about Space Groups
   Type(Table_Equiv_Type),       allocatable, dimension(:)   :: System_Equiv
   Type(Wyck_Info_Type),         allocatable, dimension(:)   :: Wyckoff_Info    ! Wyckoff information

   Type(Err_Text_Type)                     :: Mess_FindFMT = &  ! Text composed of a maximum of 5 lines to inform about position or error (findFMT)
        Err_Text_Type(0,(/" "," "," "," "," "/))

   Type(Xray_Wavelength_Type), parameter, dimension(7) :: XRAY_WAVELENGTHS =(/ &  ! Tabulated K-Series for Xray
        Xray_Wavelength_type("CR",(/2.28988,2.29428/),2.08480), &
        Xray_Wavelength_type("FE",(/1.93631,1.94043/),1.75650), &
        Xray_Wavelength_type("CU",(/1.54059,1.54431/),1.39220), &
        Xray_Wavelength_type("MO",(/0.70932,0.71360/),0.63225), &
        Xray_Wavelength_type("AG",(/0.55942,0.56380/),0.49708), &
        Xray_Wavelength_type("CO",(/1.78919,1.79321/),1.62083), &
        Xray_Wavelength_type("NI",(/1.65805,1.66199/),1.50017)  /)


End Module CFML_DefPar
