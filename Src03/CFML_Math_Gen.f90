!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2015  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_Math_General
!!----   INFO: Mathematic general utilities for use in Crystallography and
!!----         Solid State Physics and Chemistry.
!!----
!!---- HISTORY
!!----    Updated: 19/03/2015
!!----
!!
 Module CFML_Math_General
    !---- Use Modules ----!
    Use CFML_GlobalDeps

    !---- Definitions ----!
    implicit none

    private

    !---- Public Procedures ----!
    public :: Factorial, Pgcd, Ppcm, Modulo_Lat, Co_Prime

    public :: Acosd, Asind, Atan2d, Atand, Cosd, Sind, Tand, Negligible,           &
              Co_Linear, Equal_Matrix, Equal_Vector, Locate, Outerprod, Trace,     &
              Zbelong, Norm, Scalar, In_limits, Lower_Triangular, Upper_Triangular

    public :: Init_Err_Mathgen, Invert_Matrix, LU_Decomp, LU_Backsub, Matinv,        &
              Sort_Strings, Spline, Splint, Set_Epsg, Set_Epsg_Default,In_Sort,      &
              First_Derivative, Second_Derivative, SmoothingVec, Points_in_Line2D,   &
              Co_Prime_vector

    public :: RTan, Determinant, Diagonalize_Sh, Linear_Dependent, Rank, Sort,   &
              Svdcmp, Swap


    !--------------------!
    !---- PARAMETERS ----!
    !--------------------!
    integer, parameter, dimension(1000) :: primes =                                       &  ! List of the first 1000 prime numbers.
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

    real(kind=cp), parameter :: EP_SS=1.0E-12_cp   ! Internal epsilon value used for comparison in matrix operations

    !-------------------!
    !---- VARIABLES ----!
    !-------------------!
    character(len=256), public :: ERR_MathGen_Mess     ! String containing information about the last error
    logical,            public :: err_mathgen          ! Logical Variable indicating an error in CFML_Math_General module
    real(kind=cp)              :: epss=1.0E-5_cp       ! Internal epsilon value used for comparing reals to integers

    !---------------------------------!
    !---- Interfaces - Overloaded ----!
    !---------------------------------!
    Interface  Acosd
       Module Procedure Acosd_dp
       Module Procedure Acosd_sp
    End Interface

    Interface  Asind
       Module Procedure Asind_dp
       Module Procedure Asind_sp
    End Interface

    Interface  Atan2d
       Module Procedure Atan2d_dp
       Module Procedure Atan2d_sp
    End Interface

    Interface  Atand
       Module Procedure Atand_dp
       Module Procedure Atand_sp
    End Interface

    Interface  Cosd
       Module Procedure Cosd_dp
       Module Procedure Cosd_sp
    End Interface

    Interface  Sind
       Module Procedure Sind_dp
       Module Procedure Sind_sp
    End Interface

    Interface  Tand
       Module Procedure Tand_dp
       Module Procedure Tand_sp
    End Interface

    Interface  Negligible
       Module Procedure Negligibler
       Module Procedure Negligiblec
    End Interface

    Interface  Co_Linear
       Module Procedure Co_linear_C
       Module Procedure Co_linear_I
       Module Procedure Co_linear_R
    End Interface

    Interface  Equal_Matrix
       Module Procedure Equal_Matrix_I
       Module Procedure Equal_Matrix_R
    End Interface

    Interface  Equal_Vector
       Module Procedure Equal_Vector_I
       Module Procedure Equal_Vector_R
    End Interface

    Interface  Locate
       Module Procedure Locate_I
       Module Procedure Locate_R
       Module Procedure Locate_Ib
       Module Procedure Locate_Rb
    End Interface

    Interface  Lower_Triangular
       Module Procedure Lower_Triangular_I
       Module Procedure Lower_Triangular_R
    End Interface

    Interface Norm
       Module Procedure Norm_I
       Module Procedure Norm_R
    End Interface Norm

    Interface  Outerprod
       Module Procedure Outerprod_dp
       Module Procedure Outerprod_sp
    End Interface

    Interface Scalar
       Module Procedure Scalar_I
       Module Procedure Scalar_R
    End Interface Scalar

    Interface  Trace
       Module Procedure Trace_C
       Module Procedure Trace_I
       Module Procedure Trace_R
    End Interface

    Interface  Upper_Triangular
       Module Procedure Upper_Triangular_I
       Module Procedure Upper_Triangular_R
    End Interface

    Interface  Zbelong
       Module Procedure ZbelongM
       Module Procedure ZbelongN
       Module Procedure ZbelongV
    End Interface

    Interface  Rtan
       Module Procedure Rtan_dp
       Module Procedure Rtan_sp
    End Interface

    Interface  Determinant
       Module Procedure Determinant_c
       Module Procedure Determinant_r
    End Interface

    Interface  Diagonalize_SH
       Module Procedure Diagonalize_HERM
       Module Procedure Diagonalize_SYMM
    End Interface

    Interface In_Limits
       Module Procedure In_Limits_int
       Module Procedure In_Limits_dp
       Module Procedure In_Limits_sp
    End Interface

    Interface  Linear_Dependent
       Module Procedure Linear_Dependentc
       Module Procedure Linear_Dependenti
       Module Procedure Linear_Dependentr
    End Interface

    Interface  Rank
       Module Procedure Rank_dp
       Module Procedure Rank_sp
    End Interface

    Interface  Sort
       Module Procedure Sort_I
       Module Procedure Sort_R
    End Interface

    Interface  Svdcmp
       Module Procedure Svdcmp_dp
       Module Procedure Svdcmp_sp
    End Interface

    Interface Swap
        Module Procedure swap_c
        Module Procedure swap_cm
        Module Procedure swap_cv
        Module Procedure swap_i
        Module Procedure swap_im
        Module Procedure swap_iv
        Module Procedure swap_r
        Module Procedure swap_rm
        Module Procedure swap_rv
        Module Procedure masked_swap_r
        Module Procedure masked_swap_rm
        Module Procedure masked_swap_rv
    End interface

 Contains

    !---- Functions ----!

    !!----
    !!---- ELEMENTAL FUNCTION ACOSD
    !!----
    !!----    Inverse cosine function -> output in Degrees
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ Elemental Function Acosd_dp
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse cosine function -> output in Degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Acosd_dp(x) Result(arc_cos)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: arc_cos

       if (abs(x) > 1.0_dp ) then
          if (x > 0.0_dp)  then
             arc_cos=0.0_dp
          else
             arc_cos=180.0_dp
          end if
       else
          arc_cos=acos(x)*to_DEG
       end if

       return
    End Function Acosd_dp

    !!--++
    !!--++ ELEMENTAL FUNCTION ACOSD_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse cosine function -> output in Degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Acosd_sp(x) Result(arc_cos)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: arc_cos

       if (abs(x) > 1.0_sp ) then
          if (x > 0.0_sp)  then
             arc_cos=0.0_sp
          else
             arc_cos=180.0_sp
          end if
       else
          arc_cos=acos(x)*to_DEG
       end if

       return
    End Function Acosd_sp

    !!----
    !!---- FUNCTION ASIND
    !!----
    !!----    Inverse sine function -> output in Degrees
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ ELEMENTAL FUNCTION ASIND_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse sine function -> output in Degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Asind_dp(x) Result(arc_sin)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: arc_sin

       if (abs(x) > 1.0_dp ) then
          if (x > 0.0_dp) then
             arc_sin=90.0_dp
          else
             arc_sin=-90.0_dp
          end if
       else
          arc_sin=asin(x)*to_DEG
       end if

       return
    End Function Asind_dp

    !!--++
    !!--++ ELEMENTAL FUNCTION ASIND_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse sine function -> output in Degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Asind_sp(x) Result(arc_sin)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: arc_sin

       if (abs(x) > 1.0_sp ) then
          if (x > 0.0_sp) then
             arc_sin=90.0_sp
          else
             arc_sin=-90.0_sp
          end if
       else
          arc_sin=asin(x)*to_DEG
       end if

       return
    End Function Asind_sp

    !!----
    !!---- ELEMENTAL FUNCTION ATAN2D
    !!----
    !!----    Inverse tangent function of y/x
    !!----    y,x have the same units -> output in Degrees
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ ELEMENTAL FUNCTION ATAN2D_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function of y/x
    !!--++    y,x have the same units -> output in Degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Atan2d_dp(y,x) Result(atand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: y,x
       real(kind=dp)             :: atand

       atand=atan2(y,x)*to_DEG

       return
    End Function Atan2d_dp

    !!--++
    !!--++ ELEMENTAL FUNCTION ATAN2D_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function of y/x
    !!--++    y,x have the same units -> output in Degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Atan2d_sp(y,x) Result(atande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: y,x
       real(kind=sp)             :: atande

       atande=atan2(y,x)*to_DEG

       return
    End Function Atan2d_sp

    !!----
    !!---- ELEMENTAL FUNCTION ATAND
    !!----
    !!----    Inverse tangent function, X no units -> output in Degrees
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ ELEMENTAL FUNCTION ATAND_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function, X no units -> output in Degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Atand_dp(x) Result(atand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: atand

       atand=atan(x)*to_DEG

       return
    End Function Atand_dp

    !!--++
    !!--++ FUNCTION ATAND_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function, X no units -> output in Degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Atand_sp(x) Result(atande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: atande

       atande=atan(x)*to_DEG

       return
    End Function Atand_sp

    !!----
    !!---- ELEMENTAL FUNCTION COSD
    !!----
    !!----    Cosine function, X in degrees
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ ELEMENTAL FUNCTION COSD_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Cosine function, X in degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Cosd_dp(x) Result(cosine)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: cosine

       cosine=cos(to_RAD*x)

       return
    End Function Cosd_dp

    !!--++
    !!--++ ELEMENTAL FUNCTION COSD_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Cosine function, X in degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Cosd_sp(x) Result(cosine)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: cosine

       cosine=cos(to_RAD*x)

       return
    End Function Cosd_sp

    !!----
    !!---- ELEMENTAL FUNCTION SIND
    !!----
    !!----    Sine function, X in degrees
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ ELEMENTAL FUNCTION SIND_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sine function, X in degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Sind_dp(x) Result(sine)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: sine

       sine=sin(to_RAD*x)

       return
    End Function Sind_dp

    !!--++
    !!--++ ELEMENTAL FUNCTION SIND_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sine function, X in degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Sind_sp(x) Result(sine)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: sine

       sine=sin(to_RAD*x)

       return
    End Function Sind_sp

    !!----
    !!---- ELEMENTAL FUNCTION TAND
    !!----
    !!----    Tangent function, X in degrees
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ ELEMENTAL FUNCTION TAND_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Tangent function, X in degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Tand_dp(x) Result(tand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: tand

       tand=tan(to_RAD*x)

       return
    End Function Tand_dp

    !!--++
    !!--++ ELEMENTAL FUNCTION TAND_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Tangent function, X in degrees
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Tand_sp(x) Result(tande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: tande

       tande=tan(to_RAD*x)

       return
    End Function Tand_sp

    !!----
    !!---- ELEMENTAL FUNCTION FACTORIAL
    !!----
    !!----    Factorial of N
    !!----
    !!---- Update: 11/07/2015
    !!
    Elemental Function Factorial(n) Result(fact)
       !---- Argument ----!
       integer, intent(in) :: n
       integer             :: fact

       !---- Local variables ----!
       integer   :: nt, np

       if (n ==0) then
          fact=1
       else
          nt=1
          np=abs(n)
          do
             nt=nt*np
             np=np-1
             if (np == 1) exit
          end do
          fact=nt
       end if

       return
    End Function Factorial

    !!----
    !!---- Elemental Function Negligible
    !!----
    !!----    Provides the value .TRUE. if the real/complex
    !!----    number V is less than EPS
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ Elemental Function Negligiblec
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if a complex number is negligible
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Negligiblec(v) Result(Neglig)
       !---- Argument ----!
       complex, intent( in) :: v
       logical              :: Neglig

       Neglig=.false.
       if (abs(v) > epss) return
       Neglig=.true.

       return
    End Function Negligiblec

    !!--++
    !!--++ ELEMENTAL FUNCTION NEGLIGIBLER
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is negligible (abs < EPSS)
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Elemental Function Negligibler(v) Result(neglig)
       !---- Argument ----!
       real(kind=cp), intent( in) :: v
       logical                    :: Neglig

       Neglig=.false.
       if (abs(v) > epss) return
       Neglig=.true.

       return
    End Function Negligibler

    !!----
    !!---- FUNCTION PGCD
    !!----
    !!----    Function calculating the maximum common divisor of two integers
    !!----
    !!---- Update: 11/07/2015
    !!
    Function Pgcd(a,b) Result(mcd)
       !---- Arguments ----!
       integer, intent(in) :: a,b
       integer             :: mcd

       !---- Local variables ----!
       integer  :: u,v,m

       u=max(a,b)
       v=min(a,b)
       m=0
       do
          if (m == 1) exit
          m=mod(u,v)
          u=v
          v=m
       end do
       mcd=u

       return
    End Function Pgcd

    !!----
    !!---- FUNCTION PPCM
    !!----
    !!----    Function calculating the minimum common multiple of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Function Ppcm(a,b) result(mcm)
       !---- Arguments ----!
       integer, intent(in) :: a,b
       integer             :: mcm

       !---- Local variables ----!
       integer :: u,v,w,i

       u=max(a,b)
       v=min(a,b)
       mcm=1
       if (v <= 1) then
          mcm=u
          return
       end if
       w=int(sqrt(real(u)))+1
       do i=2,w
          do
             if(.not. ((mod(u,i)==0) .or. (mod(v,i)==0)) ) exit
             mcm=mcm*i
             if (modulo(u,i) == 0) u=u/i
             if (modulo(v,i) == 0) v=v/i
          end do
       end do

       return
    End Function Ppcm

    !!----
    !!---- LOGICAL FUNCTION CO_LINEAR
    !!----
    !!----    Provides the value .TRUE. if the vectors A and B are co-linear
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ LOGICAL FUNCTION CO_LINEAR_C
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two complex vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_C(a,b,n) Result(co_linear)
       !---- Argument ----!
       complex, dimension(:), intent(in) :: a,b           ! Complex vectors
       integer,               intent(in) :: n             ! Dimension
       logical                           :: co_linear

       !---- Local variables ----!
       integer :: i,ia,ib
       complex :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > epss) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=a(ia)/b(ib)
          do i=1,n
             if (abs(a(i)-c*b(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_C

    !!--++
    !!--++ LOGICAL FUNCTION CO_LINEAR_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are co-linear
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Co_linear_I(a,b,n) Result(co_linear)
       !---- Argument ----!
       integer, dimension(:), intent(in) :: a,b            ! Integer vectors
       integer,               intent(in) :: n              ! Dimension
       logical                           :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib
       real(kind=cp) :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > 0) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > 0) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=real(a(ia))/real(b(ib))
          do i=1,n
             if (abs( real(a(i))-c*real(b(i)) ) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_I

    !!--++
    !!--++ LOGICAL FUNCTION CO_LINEAR_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_R(a,b,n) Result(co_linear)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in) :: a,b           ! Real vectors
       integer,                     intent(in) :: n             ! Dimension
       logical                                 :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib
       real(kind=cp) :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > epss) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=a(ia)/b(ib)
          do i=1,n
             if (abs(a(i)-c*b(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_R

    !!----
    !!---- FUNCTION CO_PRIME
    !!----
    !!----    Provides the value .TRUE. if the array V contains co-prime
    !!----    integers: there is no common divisor for all the integers.
    !!----    Only the first 1000 prime numbers are stored in the module array "primes"
    !!----    imax is the maximum prime number to be tested. It is calculated if not given.
    !!----    (imax argument made optional, really not needed)
    !!----
    !!---- Updated: 11/07/2015
    !!
    Function Co_Prime(v,imax) result(cop)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: v           ! Integer vector
       integer, optional,     intent(in) :: imax        ! Maximun prime number to be tested
       Logical                           :: cop

       !---- Local variables ----!
       integer :: i,j,im,k,dimv,imaxv,maxv

       cop=.true.
       maxv=maxval(abs(v))
       if (present(imax)) then
          imaxv=imax
       else
          imaxv=maxv
       end if

      !> If the maximum value of the indices is 1 they are coprimes
      if (maxv == 1) return
      if (maxv == 0) then
         cop=.false.
         return
      end if

      !> Search the maximum prime number to be tested
      if (imaxv > 7919) then
         im=1000
      else
         do i=1,1000
            if(imaxv > primes(i)) cycle
            im=i
            exit
         end do
      end if

      !> Indices greater than 1
      dimv=size(v)
      do_p: do i=1,im
         k=primes(i)
         do j=1,dimv
            if( mod(v(j),k) /= 0) cycle do_p
         end do
         cop=.false.
         exit
      end do do_p

      return
    End Function Co_Prime

    !!----
    !!---- LOGICAL FUNCTION EQUAL_MATRIX
    !!----
    !!----    Provides the value .TRUE. if the array A is equal to array B in NxN
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ LOGICAL FUNCTION EQUAL_MATRIX_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Equal_Matrix_I(a,b,n) result(info)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: a,b        ! Arrays
       integer                , intent(in) :: n          ! Dimension
       logical                             :: info

       !---- Local variables ----!
       integer :: i,j

       info=.false.
       do i=1,n
          do j=1,n
             if (a(i,j) /= b(i,j)) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_I

    !!--++
    !!--++ LOGICAL FUNCTION EQUAL_MATRIX_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Matrix_R(a,b,n) result(info)
       !---- Argument ----!
       real(kind=cp), dimension(:,:)   , intent(in) :: a,b    ! Arrays
       integer,                          intent(in) :: n      ! Dimension
       logical                                      :: info

       !---- Local variables ----!
       integer :: i,j

       info=.false.
       do i=1,n
          do j=1,n
             if (abs(a(i,j) - b(i,j)) > epss ) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_R

    !!----
    !!---- LOGICAL FUNCTION EQUAL_VECTOR
    !!----
    !!----    Provides the value .TRUE. if the vector A is equal to vector B
    !!----
    !!---- Update: 01/07/2015
    !!

    !!--++
    !!--++ Logical Function Equal_Vector_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Vector_I(a,b,n) result(info)
       !---- Argument ----!
       integer, dimension(:),   intent(in) :: a,b         ! Vectors
       integer                , intent(in) :: n           ! Dimension
       logical                             :: info

       !---- Local variables ----!
       integer :: i

       info=.false.
       do i=1,n
          if (a(i) /= b(i)) return
       end do
       info=.true.

       return
    End Function Equal_Vector_I

    !!--++
    !!--++ Logical Function Equal_Vector_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real(kind=sp) vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Vector_R(a,b,n) result(info)
       !---- Argument ----!
       real(kind=cp), dimension(:)   ,   intent(in) :: a,b      ! Vectors
       integer,                          intent(in) :: n        ! Dimension
       logical                                      :: info

       !---- Local variables ----!
       integer :: i

       info=.false.
       do i=1,n
          if (abs(a(i) - b(i)) > epss ) return
       end do
       info=.true.

       return
    End Function Equal_Vector_R

    !!----
    !!---- FUNCTION IN_LIMITS
    !!----
    !!----   Logical function that is true if all the components of the vector vect
    !!----   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!----
    !!----   Updated: 11/07/2015
    !!

    !!--++
    !!--++ FUNCTION IN_LIMITS_INT
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: 11/07/2015
    !!
    Function in_limits_int(n,limits,vect) result(ok)
       !---- Arguments ----!
       integer,                 intent(in) :: n
       integer, dimension(:,:), intent(in) :: limits   ! Normally (2,n)
       integer, dimension(n),   intent(in) :: vect
       logical                             :: ok

       !---- Local Variables ----!
       integer :: i

       !> Init
       ok=.true.

       do i=1,n
          if (vect(i) >= limits(1,i) .and. vect(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_int

    !!--++
    !!--++ FUNCTION IN_LIMITS_DP
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: 11/07/2015
    !!
    Function in_limits_dp(n,limits,vect) result(ok)
       !---- Arguments ----!
       integer,                       intent(in) :: n
       real(kind=dp), dimension(:,:), intent(in) :: limits   ! Normally (2,n)
       real(kind=dp), dimension(n),   intent(in) :: vect
       logical :: ok

       !---- Local Variables ----!
       integer :: i

       !> Init
       ok=.true.

       do i=1,n
          if (vect(i) >= limits(1,i) .and. vect(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_dp

    !!--++
    !!--++ FUNCTION IN_LIMITS_SP
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: 11/07/2015
    !!
    Function in_limits_sp(n,limits,vect) result(ok)
       !---- Arguments ----!
       integer,                       intent(in) :: n
       real(kind=sp), dimension(:,:), intent(in) :: limits   ! Normally (2,n)
       real(kind=sp), dimension(n),   intent(in) :: vect
       logical :: ok

       !---- Local Variables ----!
       integer :: i

       !> Init
       ok=.true.

       do i=1,n
          if (vect(i) >= limits(1,i) .and. vect(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_sp

    !!----
    !!---- FUNCTION LOCATE
    !!----
    !!----    Function for locating the index J of an array XX(N)
    !!----    satisfying:
    !!--<<
    !!----               XX(J) <= X < XX(J+1)
    !!-->>
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ FUNCTION LOCATE_I
    !!--++
    !!--++    Subroutine for locating the index J of an array XX(N)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Locate_I(xx,n,x) Result(j)
       !---- Argument ----!
       integer, dimension(:), intent(in):: xx    ! Vector
       integer ,              intent(in):: n     ! Dimension
       integer,               intent(in):: x     ! Value
       integer                          :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm

       if(x <= xx(1)) then
         j=1
         return
       end if
       if(x >= xx(n)) then
         j=n
         return
       end if
       jl=0
       ju=n+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_I

    !!--++
    !!--++ FUNCTION LOCATE_IB
    !!--++
    !!--++    Subroutine for locating the index J of an array XX(:)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Locate_Ib(xx,x) Result(j)
       !---- Argument ----!
       integer, dimension(:), intent(in):: xx     ! Vector
       integer,               intent(in):: x      ! Value
       integer                          :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm, i1,i2
       integer, dimension(1) :: mi

       mi=lbound(xx)
       i1=mi(1)
       mi=ubound(xx)
       i2=mi(1)

       if(x <= xx(i1)) then
         j=i1
         return
       end if
       if(x >= xx(i2)) then
         j=i2
         return
       end if
       jl=i1-1
       ju=i2+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(i2) > xx(i1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl
       return
    End Function Locate_Ib
    !!--++
    !!--++ FUNCTION LOCATE_R
    !!--++
    !!--++    Function for locating the index J of an array XX(N)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update:11/07/2015
    !!
    Function Locate_R(xx,n,x) Result(j)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in):: xx    ! Vector
       integer ,                    intent(in):: n     ! Dimension
       real(kind=cp),               intent(in):: x     ! Value
       integer                                :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm

       if(x <= xx(1)) then
         j=1
         return
       end if
       if(x >= xx(n)) then
         j=n
         return
       end if
       jl=0
       ju=n+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_R

    !!--++
    !!--++ FUNCTION LOCATE_RB
    !!--++
    !!--++    Function for locating the index J of an array XX(:)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Locate_Rb(xx,x) Result(j)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in):: xx
       real(kind=cp),               intent(in):: x
       integer                                :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm, i1,i2
       integer, dimension(1) :: mi

       mi=lbound(xx)
       i1=mi(1)
       mi=ubound(xx)
       i2=mi(1)

       if(x <= xx(i1)) then
         j=i1
         return
       end if
       if(x >= xx(i2)) then
         j=i2
         return
       end if
       jl=i1-1
       ju=i2+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(i2) > xx(i1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_Rb

    !!----
    !!---- FUNCTION LOWER_TRIANGULAR
    !!----
    !!----   Return a Lower triangular array
    !!----
    !!----   Updated: 01/07/2015
    !!

    !!--++
    !!--++ FUNCTION LOWER_TRIANGULAR_I
    !!--++
    !!--++   Updated: 11/07/2015
    !!
    Function Lower_Triangular_I(A,n) Result (T)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: A     ! Array
       integer,                 intent(in) :: n     ! Dimension
       integer, dimension(n,n)             :: T     ! Lower triangular array (n x n)

       !---- Local Variable ----!
       integer :: i,j,p,q,m

       m=n
       p=size(A(:,1))
       q=size(A(1,:))
       if (n > p .or. n > q) m=min(p,q)

       T=0
       do j=1,m
          do i=j,m
             T(i,j)=A(i,j)
          end do
       end do

       return
    End Function  Lower_Triangular_I

    !!--++
    !!--++ FUNCTION LOWER_TRIANGULAR_R
    !!--++
    !!--++   Updated: 11/07/2015
    !!
    Function Lower_Triangular_R(A,n) Result (T)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: A
       integer,                       intent(in) :: n
       real(kind=cp), dimension(n,n)             :: T

       !---- Local Variables ----!
       integer :: i,j,p,q,m

       m=n
       p=size(A(:,1))
       q=size(A(1,:))
       if (n > p .or. n > q) m=min(p,q)

       T=0
       do j=1,m
          do i=j,m
             T(i,j)=A(i,j)
          end do
       end do

       return
    End Function  Lower_Triangular_R

    !!----
    !!---- FUNCTION MODULO_LAT
    !!----
    !!----    Reduces a real vector to another with components in
    !!----    the interval [0,1)
    !!----
    !!---- Updated: 11/07/2015
    !!
    Function Modulo_Lat(u) result(v)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent( in) :: u   ! Vector
       real(kind=cp), dimension(1:size(u))      :: v

       v=mod(u+10.0_cp,1.0_cp)

       return
    End Function  Modulo_Lat

    !!----
    !!---- FUNCTION NORM
    !!----
    !!----    Calculate the Norm of a vector
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ FUNCTION NORM_I
    !!--++
    !!--++    Calculate the Norm of a vector
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Norm_I(X,G) Result(R)
       !---- Arguments ----!
       integer,       dimension(:),   intent(in) :: x        ! Vector
       real(kind=cp), dimension(:,:), intent(in) :: g        ! Metrics
       real(kind=cp)                             :: r

       if (size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=sqrt(dot_product(real(x), matmul(g,real(x))))
       end if

       return
    End Function Norm_I

    !!--++
    !!--++ FUNCTION NORM_R
    !!--++
    !!--++    Calculate the Norm of a vector
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Norm_R(X,G) Result(R)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: x
       real(kind=cp), dimension(:,:), intent(in) :: g
       real(kind=cp)                             :: r

       if (size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=sqrt(dot_product(x, matmul(g,x)))
       end if

       return
    End Function Norm_R


    !!----
    !!---- FUNCTION OUTERPROD
    !!----
    !!----    Computes the outer product (tensorial product) of two
    !!----    vectors to give a tensor (matrix) as the result:
    !!--<<
    !!----                   c(i,j) = a(i)*b(j).
    !!-->>
    !!--..    It uses the intrinsic Fortran 90 function SPREAD.
    !!--..    Function adapted from Numerical Recipes.
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ FUNCTION OUTERPROD_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = a(i)*b(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Outerprod_dp(a,b)  Result(c)
       !---- Arguments ----!
       real(kind=dp),dimension(:),intent(in)    :: a,b
       real(kind=dp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod_dp

    !!--++
    !!--++ FUNCTION OUTERPROD_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = a(i)*b(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Outerprod_sp(a,b)  Result(c)
       !---- Arguments ----!
       real(kind=sp),dimension(:),intent(in)    :: a,b
       real(kind=sp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod_sp

    !!----
    !!---- FUNCTION SCALAR
    !!----
    !!----    Scalar Product including metrics
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ FUNCTION SCALAR_R
    !!--++
    !!--++    Scalar Product including metrics
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Scalar_I(X,Y,G) Result(R)
       !---- Arguments ----!
       integer, dimension(:),         intent(in) :: x             ! Vector1
       integer, dimension(:),         intent(in) :: y             ! Vector2
       real(kind=cp), dimension(:,:), intent(in) :: g             ! Metrics
       real(kind=cp)                             :: r

       if (size(x)/= size(y) .or. size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=dot_product(real(x), matmul(g,real(y)))
       end if

       return
    End Function Scalar_I

    !!--++
    !!--++ FUNCTION SCALAR_R
    !!--++
    !!--++    Scalar Product including metrics
    !!--++
    !!--++ Update: 11/07/201
    !!
    Function Scalar_R(X,Y,G) Result(R)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: x      ! Vector1
       real(kind=cp), dimension(:),   intent(in) :: y      ! Vector2
       real(kind=cp), dimension(:,:), intent(in) :: g      ! Metrics
       real(kind=cp)                             :: r

       if (size(x)/= size(y) .or. size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=dot_product(x, matmul(g,y))
       end if

       return
    End Function Scalar_R

    !!----
    !!---- FUNCTION TRACE
    !!----
    !!----    Provides the trace of a complex/real or integer matrix
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ FUNCTION TRACE_C
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a complex nxn array
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Trace_C(a) Result(b)
       !---- Argument ----!
       complex, dimension(:,:), intent(in) :: a
       complex                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=(0.0,0.0)
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_C

    !!--++
    !!--++ FUNCTION TRACE_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of an integer 3x3 array
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Trace_I(a) Result(b)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: a
       integer                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_I

    !!--++
    !!--++ FUNCTION TRACE_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a real 3x3 array
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function Trace_R(a) Result(b)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: a
       real(kind=cp)                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0.0
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_R

    !!----
    !!---- FUNCTION UPPER_TRIANGULAR
    !!----
    !!----    Return the upper triangular matrix
    !!----
    !!----   Updated: 11/07/2015
    !!

    !!--++
    !!--++ FUNCTION UPPER_TRIANGULAR_I
    !!--++
    !!--++   Updated: 11/07/2015
    !!--++
    Function Upper_Triangular_I(A,n) Result (T)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: A
       integer,                 intent(in) :: n
       integer, dimension(n,n)             :: T
       integer :: i,j,p,q,m
       m=n
       p=size(A(:,1)); q=size(A(1,:))
       if(n > p .or. n > q) m=min(p,q)
       T=0
       do j=1,m
         do i=1,j
           T(i,j)=A(i,j)
         end do
       end do
    End Function  Upper_Triangular_I

    !!--++
    !!--++ FUNCTION UPPER_TRIANGULAR_R
    !!--++
    !!--++   Updated: 11/07/2015
    !!--++
    Function Upper_Triangular_R(A,n) Result (T)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: A
       integer,                       intent(in) :: n
       real(kind=cp), dimension(n,n)             :: T
       integer :: i,j,p,q,m
       m=n
       p=size(A(:,1)); q=size(A(1,:))
       if(n > p .or. n > q) m=min(p,q)
       T=0
       do j=1,m
         do i=1,j
           T(i,j)=A(i,j)
         end do
       end do
    End Function  Upper_Triangular_R

    !!----
    !!---- LOGICAL FUNCTION ZBELONG
    !!----
    !!----    Provides the value .TRUE. if the real number (or array) V is close enough
    !!----    (whithin EPS) to an integer.
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ LOGICAL FUNCTION ZBELONGM
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real array is an Integer matrix
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function ZbelongM(v) Result(belong)
       !---- Argument ----!
       real(kind=cp),   dimension(:,:), intent( in) :: v
       logical                                      :: belong

       !---- Local variables ----!
       real(kind=cp),   dimension(size(v,1),size(v,2)) :: vec

       vec= abs(real(nint (v))-v)
       belong=.not. ANY(vec > epss)

       return
    End Function ZbelongM

    !!--++
    !!--++ LOGICAL FUNCTION ZBELONGN
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is an Integer
    !!--++
    !!--++ Update: 11/07/2015
   !!
    Function ZbelongN(a) Result(belong)
       !---- Argument ----!
       real(kind=cp), intent( in) :: a
       logical                    :: belong

       belong=.false.
       if (abs(real(nint (a))-a) > epss) return
       belong=.true.

       return
    End Function ZbelongN

    !!--++
    !!--++ LOGICAL FUNCTION ZBELONGV
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real vector is an Integer vector
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Function ZbelongV(v) Result(belong)
       !---- Argument ----!
       real(kind=cp),   dimension(:), intent( in) :: v
       logical                                    :: belong

       !---- Local variables ----!
       integer                             :: i
       real(kind=cp),   dimension(size(v)) :: vec

       belong=.false.
       vec= abs(real(nint (v))-v)
       do i=1,size(v)
          if (vec(i) > epss) return
       end do
       belong=.true.

       return
    End Function ZbelongV

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- SUBROUTINE INIT_ERR_MATHGEN
    !!----
    !!----    Initialize the errors flags in CFML_Math_General
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Init_Err_MathGen()

       ERR_MathGen=.false.
       ERR_MathGen_Mess=" "

       return
    End Subroutine Init_Err_MathGen

    !!----
    !!---- SUBROUTINE SET_EPSG
    !!----
    !!----    Sets Global EPS to the value "neweps"
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Set_Epsg(Neweps)
       !---- Arguments ----!
       real(kind=cp), intent( in) :: neweps

       epss=neweps

       return
    End Subroutine Set_Epsg

    !!----
    !!---- SUBROUTINE SET_EPSG_DEFAULT
    !!----
    !!----    Sets Global EPS to the default value: epss=1.0E-5_sp
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Set_Epsg_Default()

       epss=1.0E-5_sp

       return
    End Subroutine Set_Epsg_Default

    !!----
    !!---- SUBROUTINE RTAN
    !!----
    !!----    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!----    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE RTAN_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!--++    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Rtan_dp(y,x,ang,deg)
       !---- Arguments ----!
       real(kind=dp),              Intent( In)   :: x,y
       real(kind=dp),              Intent(Out)   :: ang
       character(len=*), optional, Intent( In)   :: deg

       !---- Local variables ----!
       real(kind=dp):: abx,aby

       abx=abs(x)
       aby=abs(y)
       if ((abx < eps) .and. (aby < eps)) then
          ang = 0.0_dp
          return
       else if(abx < eps) then
          ang = pi/2.0_dp
       else if(aby < abx) then
          ang = atan(aby/abx)
          if(x < 0.0_dp) ang = pi-ang
       else
          ang = pi/2.0_dp - atan(abx/aby)
          if(x < 0.0_dp) ang = pi-ang
       end if
       if (y < 0.0_dp) ang = -ang
       if (present(deg)) ang = ang*to_deg

       return
    End Subroutine Rtan_dp

    !!--++
    !!--++ SUBROUTINE RTAN_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!--++    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Rtan_sp(y,x,ang,deg)
       !---- Arguments ----!
       real(kind=sp),              Intent( In)   :: x,y
       real(kind=sp),              Intent(Out)   :: ang
       character(len=*), optional, Intent( In)   :: deg

       !---- local variables ----!
       real(kind=sp):: abx,aby

       abx=abs(x)
       aby=abs(y)
       if ((abx < eps) .and. (aby < eps)) then
          ang = 0.0_sp
          return
       else if(abx < eps) then
          ang = pi/2.0_sp
       else if(aby < abx) then
          ang = atan(aby/abx)
          if(x < 0.0_sp) ang = pi-ang
       else
          ang = pi/2.0_sp - atan(abx/aby)
          if(x < 0.0_sp) ang = pi-ang
       end if
       if(y < 0.0_sp) ang = -ang
       if (present(deg)) ang = ang*to_deg

       return
    End Subroutine Rtan_sp

    !!----
    !!----  SUBROUTINE CO_PRIME_VECTOR
    !!----
    !!----     Calculates the co-prime vector (cop) parallel to the input vector (v)
    !!----     It uses the list of the first thousand prime numbers.
    !!----
    !!--..     copied from Nodal_Indices (Laue_Mod) in July 2013 (JRC)
    !!----
    !!----   Updated: 11/07/2015
    !!----
    Subroutine Co_Prime_Vector(V,Cop,f)
       !---- Arguments ----!
       integer, dimension(:), intent(in)  :: v              ! input integer vector
       integer, dimension(:), intent(out) :: cop            ! Output co-prime vector
       integer,  optional,    intent(out) :: f              ! Common multiplicative factor

       !---- Local variables ----!
       integer :: i,j,max_ind,k,im,dimv,n

       cop=v
       n=1
       if (present(f)) f=1
       max_ind=maxval(abs(cop))

       !> If the maximum value of the indices is 1 they are already coprimes
       if (max_ind <= 1) return

       !> Indices greater than 1
       dimv=size(v)
       im=0
       do i=1,size(primes)
          if(primes(i) > max_ind) then  !primes is an array within this module
             im=i
             exit
          end if
       end do
       if(im == 0) return
       do_p: do i=1,im
         k=primes(i)
         do
           do j=1,dimv
              if( mod(cop(j),k) /= 0) cycle do_p
           end do
           n=n*k
           cop=cop/k
         end do
       end do do_p

       if (present(f)) f=n

       return
    End Subroutine Co_Prime_vector

    !!----
    !!---- SUBROUTINE DETERMINANT
    !!----
    !!----    Calculates the determinant of a real square matrix.
    !!----    Calculates the pseudo-determinant of a complex square matrix.
    !!----    The calculated value is only useful for linear dependency purposes.
    !!----    It tell us if the complex matrix is singular or not.
    !!--..
    !!--..    Calculates the determinant of a complex square matrix selected from a rectangular
    !!--..    matrix A, n x m, where m >= n. determ=determinant_of_A(1:n,icol:icol+n-1)
    !!--..    If icol is absent, the calculation is performed as if icol=1.
    !!--..    If icol+n-1 > m, or m < n, determ is set to 0.0 and an error message is generated.
    !!----
    !!--..    P R O V I S I O N A L (The determinant of A is not calculated at present)
    !!----
    !!---- Update: 01/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE DETERMINANT_C
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real square matrix.
    !!--++    Calculates the pseudo-determinant of a complex square matrix.
    !!--++    The calculated value is only useful for linear dependency purposes.
    !!--++    It tell us if the complex matrix is singular or not.
    !!--++
    !!--++    P R O V I S I O N A L (The determinant of A is not calculated at present)
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Determinant_C(A,n,determ)
       !---- Arguments ----!
       complex, dimension(:,:), intent( in) :: A              ! Complex array
       integer,                 intent( in) :: n              ! Dimension
       real(kind=cp),           intent(out) :: determ         ! det(AR)^2 + det(AI)^2

       !---- local variables ----!
       real(kind=cp),    dimension(2*n,2*n) :: AC   !real square matrix
       real(kind=cp)                        :: d
       integer                              :: i,nn
       logical                              :: singular

       nn=2*n
       AC(  1:n ,  1:n ) =  real(A(1:n ,1:n))
       AC(n+1:nn,  1:n ) = aimag(A(1:n ,1:n))
       AC(n+1:nn,n+1:nn) =    AC(  1:n ,1:n)
       AC(  1:n ,n+1:nn) =   -AC(n+1:nn,1:n)

       call lu_decomp(ac(1:nn,1:nn),d,singular)

       if (singular) then
          determ=0.0
       else
          determ=0.0
          do i=1,nn
             d=d*sign(1.0_cp,ac(i,i))
             determ=determ+ log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Subroutine Determinant_C

    !!--++
    !!--++ SUBROUTINE DETERMINANT_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real square matrix.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Determinant_R(A,n,determ)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent( in) :: A          ! Real array
       integer,                       intent( in) :: n          ! Dimension
       real(kind=cp),                 intent(out) :: determ

       !---- local variables ----!
       real(kind=cp),    dimension(n,n)  :: AC
       real(kind=cp)                     :: d
       integer                           :: i
       logical                           :: singular

       ac=A(1:n,1:n)
       call lu_decomp(ac,d,singular)

       if (singular) then
          determ=0.0
       else
          determ=0.0
          do i=1,n
             d=d*sign(1.0_cp,ac(i,i))
             determ=determ + log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Subroutine Determinant_R

    !!----
    !!---- SUBROUTINE DIAGONALIZE_SH
    !!----
    !!----    Diagonalize Symmetric/Hermitian matrices.
    !!----    The eigen_values E_val are sorted in descending order. The columns
    !!----    of E_vect are the corresponding eigenvectors.
    !!----
    !!---- Update:11/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE DIAGONALIZE_HERM
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Diagonalize Hermitian matrices.
    !!--++    The eigen_values E_val are sorted in descending order. The columns
    !!--++    of E_vect are the corresponding eigenvectors.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Diagonalize_Herm(a,n,e_val,e_vect)
       !---- Arguments ----!
       complex,           dimension(:,:), intent( in)  :: A           ! Complex array
       integer,                           intent( in)  :: n           ! Dimension
       real(kind=cp),     dimension(:),   intent(out)  :: E_val       ! Autovalores
       complex, optional, dimension(:,:), intent(out)  :: E_vect      ! Autovectores

       !---- Local variables ----!
       real(kind=cp),        dimension(2*n,2*n)   :: aux
       real(kind=cp),        dimension(2*n)       :: e,d
       integer :: nn

       e_val=0.0
       call init_err_mathgen()
       if (n > size(A,1) .or. n > size(A,2)) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Diagonalize_HERM: Error in dimension of input matrix: A(m,m) with m < n "
          return
       end if

       nn=2*n
       aux(  1:n ,  1:n ) =  real(a(1:n ,1:n))   !      (  U   V )
       aux(n+1:nn,n+1:nn) =  real(a(1:n ,1:n))   !   M=(          ),   A = U + i V
       aux(n+1:nn,  1:n ) = aimag(a(1:n ,1:n))   !      ( -V   U )
       aux(  1:n ,n+1:nn) =-aimag(a(1:n ,1:n))   !

       if (present(E_vect)) then
          call tred2(aux,nn,d,e)
          call tqli2(d,e,nn,aux)
          call eigsrt(d,aux,nn,1)
          e_vect(1:n,1:n)=cmplx(aux(1:n,1:nn:2),aux(n+1:nn,1:nn:2))
       else
          call tred1(aux,nn,d,e)
          call tqli1(d,e,nn)
          call eigsrt(d,aux,nn,0)
       end if
       e_val(1:n)=d(1:nn:2)

       return
    End Subroutine Diagonalize_Herm

    !!--++
    !!--++ SUBROUTINE DIAGONALIZE_SYMM
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Diagonalize symmetric matrices
    !!--++    The eigen_values E_val are sorted in descending order. The columns
    !!--++    of E_vect are the corresponding eigenvectors.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Diagonalize_Symm(A,n,E_Val,E_vect)
       !---- Arguments ----!
       real(kind=cp),           dimension(:,:), intent( in)  :: A          ! Array
       integer,                                 intent( in)  :: n          ! Dimension
       real(kind=cp),           dimension(:),   intent(out)  :: E_val      ! Autovalores
       real(kind=cp), optional, dimension(:,:), intent(out)  :: E_vect     ! Autovectores

       !---- Local variables ----!
       real(kind=cp),        dimension(n,n)   :: aux
       real(kind=cp),        dimension(n)     :: e

       e_val=0.0
       call init_err_mathgen()
       if (n > size(A,1) .or. n > size(A,2)) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Diagonalize_SYMM: Error in dimension of input matrix: A(m,m) with m < n "
          return
       end if

       aux=a(1:n,1:n)
       if (present(E_vect)) then
          call tred2(aux,n,E_val,e)
          call tqli2(E_val,e,n,aux)
          call eigsrt(E_val,aux,n,1)
          e_vect(1:n,1:n)=aux
       else
          call tred1(aux,n,E_val,e)
          call tqli1(E_val,e,n)
          call eigsrt(E_val,aux,n,0)
       end if

       return
    End Subroutine Diagonalize_Symm

    !!--++
    !!--++ SUBROUTINE EIGSRT
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for sorting eigenvalues in d(n) and eigenvectors
    !!--++    in columns of v(n,n). Sorts d(n) in descending order and
    !!--++    rearranges v(n,n) correspondingly. The method is the straight
    !!--++    insertion. If io=0 order  only the eigenvalues are treated.
    !!--++    Adapted from Numerical Recipes. Valid for hermitian matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Eigsrt(d,v,n,io)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in out) :: d
       real(kind=cp), dimension(:,:), intent(in out) :: v
       integer,                       intent(in)     :: n
       integer,                       intent(in)     :: io

       !---- Local Variables ----!
       integer          :: i,j,k
       real(kind=cp)    :: p

       do i=1,n-1
          k=i
          p=d(i)
          do j=i+1,n
             if (d(j) >= p) then
                k=j
                p=d(j)
             end if
          end do
          if (k /= i) then
             d(k)=d(i)
             d(i)=p
             if (io == 1) then
                do j=1,n
                   p=v(j,i)
                   v(j,i)=v(j,k)
                   v(j,k)=p
                end do
             end if
          end if
       end do

       return
    End Subroutine Eigsrt

    !!----
    !!---- SUBROUTINE FIRST_DERIVATIVE
    !!----
    !!----    Calculate the First derivate values of the N points
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine First_Derivative(x,y,n,d2y,d1y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: x                      ! X Vector
       real(kind=cp), dimension(:), intent(in)  :: y                      ! Yi=F(Xi)
       integer ,                    intent(in)  :: n                      ! Dimension of X, Y
       real(kind=cp), dimension(:), intent(in)  :: d2y                    ! Second derivative
       real(kind=cp), dimension(:), intent(out) :: d1y                    ! First derivative

       !---- Local Variables ----!
       integer       :: i
       real(kind=cp) :: step, x0, y0, y1, y2

       do i=1,n
         if (i /= n) then
           step = x(i+1)-x(i)
         end if
         x0 = x(i) - step/2.0
         call splint(x,y, d2y, n, x0, y0)
         y1 = y0
         x0 = x(i) + step/2
         call splint(x,y, d2y, n, x0, y0)
         y2 = y0
         d1y(i) = (y2 - y1) / step
       end do

       return
    End Subroutine First_Derivative

    !!----
    !!---- SUBROUTINE IN_SORT
    !!--<<
    !!----    Subroutine to order in ascending mode the integer array "id".
    !!----    The input value "n" is the number of items to be ordered in "id".
    !!----    The array "p" is the initial pointer to "id" (coming from a previous call)
    !!----    The final pointer holding the order of items.
    !!-->>
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine In_Sort(id,n,p,q)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: id  !Integer array to be sorted
       integer,               intent(in) :: n   !Number items in the array
       integer, dimension(:), intent(in) :: p   !Initial pointer from a previous related call
       integer, dimension(:), intent(out):: q   !Final pointer doing the sort of id

       !--- Local Variables ----!
       integer :: i,j,k,l,m
       integer, dimension(:),allocatable :: it

       l=minval(id)
       m=maxval(id)
       l=l-1
       m=m-l
       allocate(it(m))
       it(1:m)=0
       do i=1,n
          j=id(p(i))-l
          it(j)=it(j)+1
       end do
       j=0
       do i=1,m
          k=j
          j=j+it(i)
          it(i)=k
       end do
       do i=1,n
          j=id(p(i))-l
          it(j)=it(j)+1
          j=it(j)
          q(j)=p(i)
       end do

       return
    End Subroutine In_Sort

    !!----
    !!---- Subroutine INVERT_MATRIX
    !!--<<
    !!----    Subroutine to invert a real matrix using LU decomposition.
    !!----    In case of singular matrix (singular=.true.) instead of the inverse
    !!----    matrix, the subroutine provides the LU decomposed matrix as used
    !!----    in Numerical Recipes.
    !!----    The input matrix is preserved and its inverse (or its LU decomposition)
    !!----    is provided in "b". The optional argument "perm" holds the row permutation
    !!----    performed to obtain the LU decomposition.
    !!-->>
    !!----
    !!---- Update:11/07/2015
    !!
    Subroutine Invert_Matrix(a,b,singular,perm)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:),  intent(in ) :: a
       real(kind=cp), dimension(:,:),  intent(out) :: b
       logical,                        intent(out) :: singular
       integer, dimension(:),optional, intent(out) :: perm

       !---- Local variables ----!
       integer                                       :: i,n
       integer,       dimension(size(a,1))           :: indx
       real(kind=cp)                                 :: d, det
       real(kind=cp), dimension(size(a,1),size(a,1)) :: lu

       n=size(a,1)
       lu=a(1:n,1:n)

       call LU_Decomp(lu,d,singular,indx)
       if (present(perm)) perm(1:n)=indx(1:n)

       if (singular) then
          b=lu
          return
       else
          det=0.0
          do i=1,n
             d=d*sign(1.0_cp,lu(i,i))
             det=det + log(abs(lu(i,i)))
          end do
          det=d*exp(det)
          if (abs(det) <= 1.0e-36) then
             singular=.true.
             b=lu
             return
          end if
       end if

       b=0.0
       do i=1,n
          b(i,i)=1.0
          call LU_backsub(lu,indx,b(:,i))
       end do

       return
    End Subroutine Invert_Matrix

    !!----
    !!---- SUBROUTINE LINEAR_DEPENDENT
    !!--<<
    !!----    Provides the value .TRUE. if the vector A is linear dependent of the
    !!----    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!----    are the number of rows and columns of B to be considered. The actual
    !!----    dimension of vector a should be na=max(nb,mb).
    !!----    The problem is equivalent to determine the rank (in algebraic sense)
    !!----    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!----    case it is supposed that na = mb and in the second na = nb.
    !!----    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!----    is generated. The function uses floating arithmetic for all types.
    !!-->>
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE LINEAR_DEPENDENTC
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++
    !!--++    For the case of complex vectors in Cn the problem can be reduced to real vectors
    !!--++    of dimension R2n. Each complex vector contributes as two real vectors of dimension
    !!--++    2n: (R,I) and (-I,R). A complex vector V is linearly dependent on n complex vectors
    !!--++    if V can be written as: V = Sigma{j=1,n}(Cj.Vj), with Cj complex numbers and Vj
    !!--++    having n complex components. One may write:
    !!--++
    !!--++     V = Sigma{j=1,n}(Cj.Vj)
    !!--++     (R,I) = Sigma{j=1,n} (Cjr Vj + i Cji Vj) = Sigma{j=1,n} (Cjr (Rj,Ij) +  Cji (-Ij,Rj) )
    !!--++     (R,I) = Sigma{j=1,n} (aj (Rj,Ij) + bj (-Ij,Rj) )  = Sigma{j=1,2n} (Aj.Uj)
    !!--++     Were Uj=(Rj,Ij) and U(j+1)= (-Ij,Rj)
    !!--++
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Linear_DependentC(A,na,B,nb,mb,info)
       !---- Arguments ----!
       complex, dimension(:),   intent(in)  :: a
       complex, dimension(:,:), intent(in)  :: b
       integer,                 intent(in)  :: na,nb,mb
       logical,                 intent(out) :: info

       !---- Local variables ----!
       integer                                                     :: r,n1
       real(kind=dp), parameter                                    :: tol= 100.0_dp*deps
       real(kind=dp), dimension(2*max(nb+1,mb+1),2*max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentC: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then
          n1=2*nb+1
          if(n1+1 > 2*mb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(b(1:nb,1:mb))
          c(1:nb,     mb+1:mb+na) = aimag(b(1:nb,1:mb))
          c(nb+1:2*nb,      1:mb) =-aimag(b(1:nb,1:mb))
          c(nb+1:2*nb,mb+1:mb+na) =  real(b(1:nb,1:mb))
          c(n1,             1:mb) =  real(a(1:na))
          c(n1,      mb+1:mb+na ) = aimag(a(1:na))
          c(n1+1,           1:mb) =-aimag(a(1:na))
          c(n1+1,    mb+1:mb+na ) =  real(a(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*mb)) info=.false.
       else if( na == nb) then
          n1=2*mb+1
          if(n1+1 > 2*nb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(b(1:nb,1:mb))
          c(nb+1:nb+na,     1:mb) = aimag(b(1:nb,1:mb))
          c(1:nb,      mb+1:2*mb) =-aimag(b(1:nb,1:mb))
          c(nb+1:nb+na,mb+1:2*mb) =  real(b(1:nb,1:mb))
          c(1:na,             n1) =  real(a(1:na))
          c(nb+1:nb+na,       n1) = aimag(a(1:na))
          c(1:na,           1+n1) =-aimag(a(1:na))
          c(nb+1:nb+na,     1+n1) =  real(a(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentC: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentC

    !!--++
    !!--++ SUBROUTINE LINEAR_DEPENDENTI
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Linear_DependentI(A,na,B,nb,mb,info)
       !---- Arguments ----!
       integer, dimension(:),   intent(in)  :: a
       integer, dimension(:,:), intent(in)  :: b
       integer,                 intent(in)  :: na,nb,mb
       logical,                 intent(out) :: info

       !---- Local variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentI: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(b(1:nb,1:mb))
          c(n1,  1:mb)=real(a(1:na))      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) info=.false.
       else if( na == nb) then
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(b(1:nb,1:mb))
          c(1:nb,  n1)=real(a(1:na))     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentI: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentI

    !!--++
    !!--++ SUBROUTINE LINEAR_DEPENDENTR
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Linear_DependentR(A,na,B,nb,mb,info)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in)  :: a
       real(kind=cp), dimension(:,:), intent(in)  :: b
       integer,                       intent(in)  :: na,nb,mb
       logical,                       intent(out) :: info

       !---- Local Variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentR: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then    !Vector added as an additional row
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=b(1:nb,1:mb)
          c(n1,  1:mb)=a(1:na)      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) info=.false.
       else if( na == nb) then   !Vector added as an additional column
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=b(1:nb,1:mb)
          c(1:nb,  n1)=a(1:na)     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentR: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentR

    !!----
    !!---- SUBROUTINE LU_BACKSUB
    !!--<<
    !!----    Adapted from Numerical Recipes.
    !!----    Solves the set of N linear equations A  X = B. Here the N x N matrix A is input,
    !!----    not as the original matrix A, but rather as its LU decomposition, determined
    !!----    by the routine LU_DECOMP. INDX is input as the permutation vector of length N
    !!----    returned by LU_DECOMP. B is input as the right-hand-side vector B,
    !!----    also of length N, and returns with the solution vector X.
    !!----    A and INDX are not modified by this routine and can be left in place for successive calls
    !!----    with different right-hand sides B. This routine takes into account the possibility that B will
    !!----    begin with many zero elements, so it is efficient for use in matrix inversion.
    !!-->>
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine LU_Backsub(a,indx,b)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in)     :: a
       integer,         dimension(:), intent(in)     :: indx
       real(kind=cp),   dimension(:), intent(in out) :: b

       !---- Local Variables ----!
       integer       :: i,ii,ll,n
       real(kind=cp) :: summ

       n=size(a,1)
       ii=0              !When ii is set to a positive value, it will become the index
       do i=1,n          !of the first nonvanishing element of b. We now do
          ll=indx(i)     !the forward substitution. The only new wrinkle is to
          summ=b(ll)     !unscramble the permutation as we go.
          b(ll)=b(i)
          if (ii /= 0) then
             summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
          else if(summ /= 0.0) then   !A nonzero element was encountered, so from now on
             ii=i                       !we will have to do the dot product above.
          end if
          b(i)=summ
       end do

       do i=n,1,-1       !Now we do the backsubstitution
          b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
       end do

       return
    End Subroutine LU_Backsub

    !!----
    !!---- SUBROUTINE LU_DECOMP
    !!--<<
    !!----    Subroutine to make the LU decomposition of an input matrix A.
    !!----    The input matrix is destroyed and replaced by a matrix containing
    !!----    in its upper triangular part (plus diagonal) the matrix U. The
    !!----    lower triangular part contains the nontrivial part (Lii=1) of matrix L.
    !!----    The output is rowwise permutation of the initial matrix. The vector INDX
    !!----    recording the row permutation. D is output as +/-1 depending on whether
    !!----    the number of row interchanges was even or odd, respectively.
    !!-->>
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine LU_Decomp(a,d,singular,indx)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       real(kind=cp),                 intent(out)    :: d
       logical,                       intent(out)    :: singular
       integer,  dimension(:), intent(out), optional :: indx

       !---- Local variables ----!
       real(kind=cp), dimension(size(a,1)):: vv  !vv stores the implicit scaling of each row.
       real(kind=cp), parameter           :: vtiny = 1.0e-20_sp !A small number.
       integer                            :: j,imax,n

       singular=.false.
       n=size(a,1)
       d=1.0                      !No row interchanges yet.
       vv=maxval(abs(a),dim=2)    !Loop over rows to get the implicit scaling information.
       if (any(abs(vv) <= vtiny)) then   !There is a row of zeros.
          singular=.true.
          return
       end if
       vv=1.0_sp/vv     !Save the scaling.
       do j=1,n
          imax=(j-1)+maxloc(vv(j:n)*abs(a(j:n,j)),dim=1)   !Find the pivot row.
          if (j /= imax) then                         !Do we need to interchange rows?
             call swap(a(imax,:),a(j,:))              !Yes, do so...
             d=-d                                     !...and change the parity of d.
             vv(imax)=vv(j)                           !Also interchange the scale factor.
          end if
          if (present(indx)) indx(j)=imax
          if (abs(a(j,j)) <= vtiny) then !If the pivot element is zero the matrix is singular.
             a(j,j)=vtiny                !(at least to the precision of the algorithm)
             singular=.true.             !For some applications on singular matrices,
             return                      !it is desirable to substitute vtiny for zero.
          end if                         !This is actually the present case
          a(j+1:n,j)=a(j+1:n,j)/a(j,j)                                    !Divide by the pivot element.
          a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))  !Reduce remaining submatrix.
       end do

       return
    End Subroutine LU_Decomp

    !!----
    !!---- SUBROUTINE MATINV
    !!----
    !!----  Subroutine for inverting a real square matrix.
    !!----  The input matrix is replaced in output with its inverse.
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Matinv(a,n)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       integer     ,                  intent(in)     :: n

       !---- Local variables ----!
       real(kind=cp)                 :: amax,savec
       integer, dimension(size(a,1)) :: ik,jk
       integer                       :: i,j,k,l

       !---- Subroutine to invert a real matrix ----!
       do k=1,n
          amax=0.0
          do
             do
                do i=k,n
                   do j=k,n
                      if (abs(amax)-abs(a(i,j)) > 0.0) cycle
                      amax=a(i,j)
                      ik(k)=i
                      jk(k)=j
                   end do
                end do
                i=ik(k)
                if (i-k < 0) cycle
                exit
             end do

             if (i-k /= 0) then
                do j=1,n
                   savec=a(k,j)
                   a(k,j)=a(i,j)
                   a(i,j)=-savec
                end do
             end if

             j=jk(k)
             if (j-k < 0) cycle
             exit
          end do

          if (j-k /= 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=a(i,j)
                a(i,j)=-savec
             end do
          end if

          do i=1,n
             if (i-k /= 0)  then
                a(i,k)=-a(i,k)/amax
             end if
          end do
          do i=1,n
             do j=1,n
                if (i-k == 0 .or. j-k == 0) cycle
                a(i,j)=a(i,j)+a(i,k)*a(k,j)
             end do
          end do
          do j=1,n
             if (j-k == 0)   cycle
             a(k,j)=a(k,j)/amax
          end do
          a(k,k)=1.0/amax
       end do     !k

       do l=1,n
          k=n-l+1
          j=ik(k)
          if (j-k > 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=-a(i,j)
                a(i,j)=savec
             end do
          end if
          i=jk(k)
          if (i-k > 0) then
             do j=1,n
                savec=a(k,j)
                a(k,j)=-a(i,j)
                a(i,j)=savec
             end do
          end if
       end do

       return
    End Subroutine Matinv

    !!--++
    !!--++ SUBROUTINE PARTITION
    !!--++    (Private)
    !!--++    Utilised by Sort_Strings.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Partition(A, Marker)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in out) :: A
       integer,                        intent(   out) :: marker

       !---- Local variables ----!
       integer                  :: i, j
       character(len=len(A(1))) :: temp
       character(len=len(A(1))) :: x      ! pivot point

       x = A(1)
       i= 0
       j= size(A) + 1

       do
          j = j-1
          do
             if (A(j) <= x) exit
             j = j-1
          end do
          i = i+1
          do
             if (A(i) >= x) exit
             i = i+1
          end do
          if (i < j) then
             !---- exchange A(i) and A(j)
             temp = A(i)
             A(i) = A(j)
             A(j) = temp
          else if (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          end if
       end do

       return
    End Subroutine Partition

    !!----
    !!---- SUBROUTINE POINTS_IN_LINE2D
    !!----
    !!----    The routine calculate N points belonging to the line defined
    !!----    by X1 and Xn with equal distance between them. XP contains
    !!----    X1,X2,.....,XN points.
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Points_In_Line2D(X1, XN, N, XP)
       !---- Arguments ----!
       real(kind=cp), dimension(2),   intent(in)  :: X1   ! Point1 in 2D
       real(kind=cp), dimension(2),   intent(in)  :: XN   ! PointN in 2D
       integer,                       intent(in)  :: N    ! Number of Total points
       real(kind=cp), dimension(:,:), intent(out) :: XP   ! List of points

       !---- Local Variables ----!
       integer :: i
       real(kind=cp)    :: ml,bl,dl,t
       real(kind=cp)    :: a,b,c,d
       real(kind=cp)    :: xa,xb

       xp=0.0

       if (n <= 1) return

       !---- Calculating the distance between two points to
       !---- eliminate rare considerations as the same point
       dl=sqrt( (xn(1)-x1(1))**2 + (xn(2)-x1(2))**2 )
       if (dl <= 0.0001) return

       !---- When N=2 is trivial case ----!
       if (n == 2) then
          xp(:,1)=x1
          xp(:,2)=xn
          return
       end if

       !---- Case 1: Y=cte ----!
       !Xn(2) and X1(2) are equal, then we have a line  with Y=cte
       if (abs(xn(2)-x1(2)) <= 0.0001) then
          dl=abs(xn(1)-x1(1))
          d=dl/real(n-1)
          xp(:,1)=x1
          if (xn(1) > x1(1)) then
             do i=2,n-1
                xp(1,i)=xp(1,i-1)+d
                xp(2,i)=xp(2,1)
             end do
          else
             do i=2,n-1
                xp(1,i)=xp(1,i-1)-d
                xp(2,i)=xp(2,1)
             end do
          end if
          xp(:,n)=xn

          return
       end if

       !---- Case 2: X=cte ----!
       !Xn(1) - X1(1) are equal, then we have a line with X=cte
       if (abs(xn(1)-x1(1)) <= 0.0001) then
          dl=abs(xn(2)-x1(2))
          d=dl/real(n-1)
          xp(:,1)=x1
          if (xn(2) > x1(2)) then
             do i=2,n-1
                xp(1,i)=xp(1,1)
                xp(2,i)=xp(2,i-1)+d
             end do
          else
             do i=2,n-1
                xp(1,i)=xp(1,1)
                xp(2,i)=xp(2,i-1)-d
             end do
          end if
          xp(:,n)=xn

          return
       end if

       !---- Case 3: General case ----!
       ml=(x1(2)-xn(2))/(x1(1)-xn(1))
       bl=x1(2) - (ml * x1(1))

       !---- Distance between X1 and XN ----!
       dl=sqrt( (xn(1)-x1(1))**2 + (xn(2)-x1(2))**2 )

       !---- Creating the list ----!
       a=ml**2 + 1.0
       b=2.0 *( ml*(bl-x1(2)) -x1(1) )

       xp(:,1)=x1
       do i=2,n-1
          t=(dl**2)*((real(i-1)/real(n-1))**2)
          c=(x1(2)-bl)**2 + x1(1)**2 - t

          xa=(-b + sqrt(b**2 - 4.0*a*c))/(2.0*a)
          xb=(-b - sqrt(b**2 - 4.0*a*c))/(2.0*a)
          if (x1(1) <= xa .and. xa <= xn(1)) then
             xp(1,i)=xa
             xp(2,i)=ml*xa+bl
          else
             xp(1,i)=xb
             xp(2,i)=ml*xb+bl
          end if
       end do
       xp(:,n)=xn

       return
    End Subroutine Points_In_Line2D

    !!----
    !!---- SUBROUTINE RANK
    !!----
    !!----    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE RANK_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Rank_dp(a,tol,r)
       !---- Arguments ----!
       real(kind=dp), dimension(:,:),intent( in)      :: a
       real(kind=dp),                intent( in)      :: tol
       integer,                      intent(out)      :: r

       !---- Arguments ----!
       real(kind=dp), dimension(size(a,1),size(a,2))  :: u
       real(kind=dp), dimension(size(a,2))            :: w
       real(kind=dp), dimension(size(a,2),size(a,2))  :: v
       integer                                        :: i

       u=a
       call svdcmp(u,w,v)
       if (ERR_MathGen) then
          r=0
       else
          r=0
          do i=1,size(a,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_dp

    !!--++
    !!--++ SUBROUTINE RANK_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Rank_sp(a,tol,r)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:),intent( in)      :: a
       real(kind=sp),                intent( in)      :: tol
       integer,                      intent(out)      :: r

       !---- Local variables ----!
       real(kind=sp), dimension(size(a,1),size(a,2))  :: u
       real(kind=sp), dimension(size(a,2))            :: w
       real(kind=sp), dimension(size(a,2),size(a,2))  :: v
       integer :: i

       u=a
       call svdcmp(u,w,v)
       if (ERR_MathGen) then
          r=0
       else
          r=0
          do i=1,size(a,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_sp

    !!----
    !!---- SUBROUTINE SECOND_DERIVATIVE
    !!----
    !!----    Calculate the second derivate of N Points
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Second_Derivative(x,y,n,d2y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: x      ! X array
       real(kind=cp), dimension(:), intent(in)  :: y      ! Yi=F(Xi)
       integer ,                    intent(in)  :: n      ! Dimension of X,Y
       real(kind=cp), dimension(:), intent(out) :: d2y    ! Second derivative

       !---- Local Variables ----!
       integer                     :: i, k
       real(kind=cp), dimension(n) :: u
       real(kind=cp)               :: yp1, ypn, sig, p, qn, un

       yp1=(y(2) - y(1))   / (x(2) - x(1))     ! derivative at point 1
       ypn=(y(n) - y(n-1)) / (x(n) - x(n-1))   ! derivative at point n

       d2y(1)=-0.5
       u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)

       do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*d2y(i-1)+2.0
          d2y(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))  &
               /(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do

       qn=0.5
       un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       d2y(n)=(un-qn*u(n-1))/(qn*d2y(n-1)+1.0)
       do k=n-1,1,-1
          d2y(k)=d2y(k)*d2y(k+1)+u(k)
       end do

       return
    End Subroutine Second_Derivative

    !!----
    !!---- SUBROUTINE SMOOTHINGVEC
    !!----
    !!----    Procedure to smooth the array values
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine SmoothingVec(Y, N, Niter, Ys)
       !---- Arguments ----!
       real(kind=cp),dimension(:),            intent(in out) :: Y          ! Array to be smoothed
       integer,                               intent(in)     :: n          ! Number of points
       integer,                               intent(in)     :: niter      ! Number of iterations
       real(kind=cp),dimension(:), optional,  intent(out)    :: Ys         ! Array smoothed

       !---- Local Variables ----!
       integer                     :: n1, n2
       integer                     :: i, iter
       real(kind=cp), dimension (n):: datYs


       n1 = 4
       n2 = n-3

       do iter = 1 ,niter
          datYs(n1-1)=((Y(n1-2)+Y(n1))*10.0+(Y(n1-3)+Y(n1+1))*5.0+Y(n1+2))/31.0
          datYs(n1-2)=((Y(n1-3)+Y(n1-1))*10.0+Y(n1)*5.0+Y(n1+1))/26.0
          datYs(n1-3)=(Y(n1-2)*10.0+Y(n1-1)*5.0+Y(n1))/16.0

          do i=n1,n2
             datYs(i)=(Y(i-3)+Y(i+3)+5.0*(Y(i-2)+Y(i+2))+10.0*(Y(i-1)+Y(i+1)))/ 32.0
          end do

          datYs(n2+1)=((Y(n2+2)+Y(n2))*10.0+(Y(n2+3)+Y(n2-1))*5.0+Y(n2-2))/31.0
          datYs(n2+2)=((Y(n2+3)+Y(n2+1))*10.0+Y(n2)*5.0+Y(n2-1))/26.0
          datYs(n2+3)=(Y(n2+2)*10.0+Y(n2+1)*5.0+Y(n2))/16.0

          if(present(Ys)) then
             Ys(1:n) = datYs(1:n)
          else
             Y(1:n) = datYs(1:n)
          end if
       end do

       return
    End Subroutine SmoothingVec

    !!---
    !!---- SUBROUTINE SORT
    !!----
    !!----    Sort an array such the a(indx(j)) is in ascending
    !!----    order for j=1,2,...,N.
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE SORT_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sort an array such the arr(indx(j)) is in ascending
    !!--++    order for j=1,2,...,N.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Sort_I(arr,n,indx)
       !---- Arguments ----!
       integer, dimension(:), intent(in ) :: arr
       integer              , intent(in ) :: n
       integer, dimension(:), intent(out) :: indx

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer                      :: i,indxt,ir,itemp,j,jstack,k,l
       integer                      :: a

       call init_Err_MathGen()
       do j=1,n
          indx(j)=j
       end do

       istack=0
       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if (arr(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l+1)) > arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
             do
                i=i+1
                if (arr(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (arr(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                ERR_MathGen=.true.
                ERR_MathGen_Mess=" NSTACK too small in SORT"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Subroutine Sort_I

    !!--++
    !!--++ SUBROUTINE SORT_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sort an array such the arr(indx(j)) is in ascending
    !!--++    order for j=1,2,...,N.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Sort_R(arr,n,indx)
       !---- Arguments ----!
       real(kind=cp),dimension(:), intent(in) :: arr
       integer,                    intent(in) :: n
       integer,      dimension(:), intent(out):: indx

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer :: i,indxt,ir,itemp,j,jstack,k,l
       real(kind=cp)    :: a

       call init_Err_MathGen()
       do j=1,n
          indx(j)=j
       end do

       istack=0
       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if (arr(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l+1)) > arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
             do
                i=i+1
                if (arr(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (arr(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                ERR_MathGen=.true.
                ERR_MathGen_Mess=" NSTACK too small in SORT"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Subroutine Sort_R

    !!----
    !!---- SUBROUTINE SORT_STRINGS
    !!----
    !!----    Sort an array of string
    !!----
    !!---- Update: 11/07/2015
    !!
    Recursive Subroutine Sort_Strings(VStr)
       !---- Argument ----!
       character(len=*), dimension(:), intent(in out) :: VStr

       !---- Local variables ----!
       integer :: iq

       if (size(VStr) > 1) then
          call Partition(VStr, iq)
          call Sort_Strings(VStr(:iq-1))
          call Sort_Strings(VStr(iq:))
       end if

       return
    End Subroutine Sort_Strings

    !!----
    !!---- SUBROUTINE SPLINE
    !!----
    !!----    Spline  N points
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Spline(x,y,n,yp1,ypn,y2)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: x       ! X array
       real(kind=cp), dimension(:), intent(in)  :: y       ! Yi=F(Xi)
       integer ,                    intent(in)  :: n       ! Dimension X,Y
       real(kind=cp),               intent(in)  :: yp1     ! Derivate of Point 1
       real(kind=cp),               intent(in)  :: ypn     ! Derivate of Point N
       real(kind=cp), dimension(:), intent(out) :: y2      ! Second derivatives at the given points

       !---- Local Variables ----!
       integer                     :: i, k
       real(kind=cp), dimension(n) :: u
       real(kind=cp)               :: sig, p, qn, un

       if (yp1 > 1.0e+30) then
          y2(1)=0.0
          u(1)=0.0
       else
          y2(1)=-0.5
          u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       end if

       do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.0
          y2(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))  &
               /(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do
       if (ypn > 1.0e+30) then
          qn=0.0
          un=0.0
       else
          qn=0.5
          un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       end if
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
       do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
       end do

       return
    End Subroutine Spline

    !!----
    !!---- SUBROUTINE SPLINT
    !!----
    !!----    Spline Interpolation
    !!----
    !!---- Update: 11/07/2015
    !!
    Subroutine Splint(xa,ya,y2a,n,x,y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: xa     ! Array X
       real(kind=cp), dimension(:), intent(in)  :: ya     ! Y=F(x)
       real(kind=cp), dimension(:), intent(in)  :: y2a    ! Second derivative
       integer ,                    intent(in)  :: n      ! Dimension X,Y,Y2
       real(kind=cp),               intent(in)  :: x      ! Point to evaluate
       real(kind=cp),               intent(out) :: y      ! Interpoled value

       !---- Local Variables ----!
       integer          :: klo, khi, k
       real(kind=cp)    :: h, a, b

       klo=1
       khi=n
       do
          if (khi-klo > 1) then
             k=(khi+klo)/2
             if (xa(k) > x) then
                khi=k
             else
                klo=k
             end if
             cycle
          else
             exit
          end if
       end do

       h=xa(khi)-xa(klo)
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)* y2a(khi))*(h**2)/6.0

       return
    End Subroutine Splint

    !!----
    !!---- SUBROUTINE SVDCMP
    !!--<<
    !!----    Given an M�N matrix A ,this routine computes its singular value decomposition,
    !!----    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!----    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!----    (not the transpose VT )is output as v .
    !!----    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!-->>
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE SVDCMP_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Given an M �N matrix A ,this routine computes its singular value decomposition,
    !!--++    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!--++    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!--++    (not the transpose VT )is output as v .
    !!--++    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Svdcmp_dp(a,w,v)
       !---- Arguments ----!
       real(kind=dp),dimension(:,:),intent(in out) ::a   ! A(m,n)
       real(kind=dp),dimension(:),  intent(   out) ::w   ! W(n)
       real(kind=dp),dimension(:,:),intent(   out) ::v   ! V(n,n)

       !---- Local variables ----!
       integer, parameter                          :: num_its=500
       integer                                     ::i,its,j,k,l,m,n,nm
       real(kind=dp)                               ::anorm,c,f,g,h,s,scal,x,y,z
       real(kind=dp),dimension(size(a,1))          ::tempm
       real(kind=dp),dimension(size(a,2))          ::rv1,tempn

       m=size(a,1)
       n=size(a,2)
       call init_err_mathgen()
       if ( .not. (size(v,1) == n .and. size(v,2) == n .and. size(w) == n)) then
          ERR_MathGen = .true.
          ERR_MathGen_Mess = " => Physical dimensions of arguments in SVDcmp_dp are not compatible "
          return
       end if
       g=0.0_dp
       scal=0.0_dp
       do i=1,n
          l=i+1
          rv1(i)=scal*g
          g=0.0_dp
          scal=0.0_dp
          if (i <=m)then
             scal=sum(abs(a(i:m,i)))
             if ( abs(scal) > tiny(1.0_dp) ) then
                a(i:m,i)=a(i:m,i)/scal
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scal*a(i:m,i)
             end if
          end if
          w(i)=scal*g
          g=0.0_dp
          scal=0.0_dp
          if ((i <=m).and.(i /=n))then
             scal=sum(abs(a(i,l:n)))
             if ( abs(scal) > tiny(1.0_dp) ) then
                a(i,l:n)=a(i,l:n)/scal
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scal*a(i,l:n)
             end if
          end if
       end do
       anorm=maxval(abs(w)+abs(rv1))
       do i=n,1,-1
          if (i <n) then
             if ( abs(g) > tiny(1.0_dp) ) then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
             end if
             v(i,l:n)=0.0_dp
             v(l:n,i)=0.0_dp
          end if
          v(i,i)=1.0_dp
          g=rv1(i)
          l=i
       end do
       do i=min(m,n),1,-1
          l=i+1
          g=w(i)
          a(i,l:n)=0.0_dp
          if ( abs(g) > tiny(1.0_dp) ) then
             g=1.0_dp/g
             tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=a(i:m,i)*g
          else
             a(i:m,i)=0.0_dp
          end if
          a(i,i)=a(i,i)+1.0_dp
       end do
       do k=n,1,-1
          do its=1,num_its
             do l=k,1,-1
                nm=l-1
                if ((abs(rv1(l))+anorm)==anorm) exit
                if ((abs(w(nm))+anorm)==anorm) then
                   c=0.0_dp
                   s=1.0_dp
                   do i=l,k
                      f=s*rv1(i)
                      rv1(i)=c*rv1(i)
                      if ((abs(f)+anorm)==anorm)exit
                      g=w(i)
                      h=hypot(f,g)
                      w(i)=h
                      h=1.0_dp/h
                      c=(g*h)
                      s=-(f*h)
                      tempm(1:m)=a(1:m,nm)
                      a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                   end do
                   exit
                end if
             end do
             z=w(k)
             if (l ==k)then
                if (z <0.0_dp)then
                   w(k)=-z
                   v(1:n,k)=-v(1:n,k)
                end if
                exit
             end if
             if (its == num_its) then
                ERR_MathGen = .true.
                ERR_MathGen_Mess = " => SVDcmp_dp: convergence not reached ! "
                return
             end if
             x=w(l)
             nm=k-1
             y=w(nm)
             g=rv1(nm)
             h=rv1(k)
             f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dp*h*y)
             g=hypot(f,1.0_dp)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.0_dp
             s=1.0_dp
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=hypot(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=hypot(f,h)
                w(j)=z
                if ( abs(z) > tiny(1.0_dp) ) then
                   z=1.0_dp/z
                   c=f*z
                   s=h*z
                end if
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
             end do
             rv1(l)=0.0_dp
             rv1(k)=f
             w(k)=x
          end do
       end do

       return
    End Subroutine Svdcmp_dp

    !!--++
    !!--++ SUBROUTINE SVDCMP_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Given an M �N matrix A ,this routine computes its singular value decomposition,
    !!--++    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!--++    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!--++    (not the transpose VT )is output as v .
    !!--++    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Svdcmp_sp(a,w,v)
       !---- Arguments ----!
       real(kind=sp),dimension(:,:),intent(in out) :: a
       real(kind=sp),dimension(:),  intent(   out) :: w
       real(kind=sp),dimension(:,:),intent(   out) :: v

       !---- Local variables ----!
       integer, parameter                          :: num_its=500
       integer                                     ::i,its,j,k,l,m,n,nm
       real(kind=sp)                               ::anorm,c,f,g,h,s,scala,x,y,z
       real(kind=sp),dimension(size(a,1))          ::tempm
       real(kind=sp),dimension(size(a,2))          ::rv1,tempn


       m=size(a,1)
       n=size(a,2)
       call init_err_mathgen()
       if ( .not. (size(v,1) == n .and. size(v,2) == n .and. size(w) == n)) then
          ERR_MathGen = .true.
          ERR_MathGen_Mess = " => Physical dimensions of arguments in SVDcmp_sp are not compatible "
          return
       end if
       g=0.0
       scala=0.0
       do i=1,n                        !Householder reduction to bidiagonal form.
          l=i+1
          rv1(i)=scala*g
          g=0.0
          scala=0.0
          if (i <=m)then
             scala=sum(abs(a(i:m,i)))
             if (abs(scala) > tiny(1.0_sp))then
                a(i:m,i)=a(i:m,i)/scala
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scala*a(i:m,i)
             end if
          end if
          w(i)=scala*g
          g=0.0
          scala=0.0
          if ((i <=m).and.(i /=n))then
             scala=sum(abs(a(i,l:n)))
             if (abs(scala) > tiny(1.0_sp))then
                a(i,l:n)=a(i,l:n)/scala
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scala*a(i,l:n)
             end if
          end if
       end do
       anorm=maxval(abs(w)+abs(rv1))
       do i=n,1,-1                    ! Accumulation of right-hand transformations.
          if (i <n)then
             if (abs(g) > tiny(1.0_sp))then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g   !Double division to avoid possible underflow.
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
             end if
             v(i,l:n)=0.0
             v(l:n,i)=0.0
          end if
          v(i,i)=1.0
          g=rv1(i)
          l=i
       end do
       do i=min(m,n),1,-1  !Accumulation of left-hand transformations.
          l=i+1
          g=w(i)
          a(i,l:n)=0.0
          if (abs(g) > tiny(1.0_sp))then
             g=1.0_sp/g
             tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=a(i:m,i)*g
          else
             a(i:m,i)=0.0
          end if
          a(i,i)=a(i,i)+1.0_sp
       end do
       do k=n,1,-1           !Diagonalization of the idiagonal form:Loop over
          do its=1,num_its    !singular values,and over allowed iterations.
             do l=k,1,-1      !Test for splitting.
                nm=l-1        !Note that rv1(1)is always zero,so can never fall through bottom of loop.
                if ((abs(rv1(l))+anorm)==anorm) exit
                if ((abs(w(nm))+anorm)==anorm) then
                   c=0.0       ! Cancellation of rv1(l),if l >1 .
                   s=1.0
                   do i=l,k
                      f=s*rv1(i)
                      rv1(i)=c*rv1(i)
                      if ((abs(f)+anorm)==anorm)exit
                      g=w(i)
                      h=hypot(f,g)
                      w(i)=h
                      h=1.0_sp/h
                      c=(g*h)
                      s=-(f*h)
                      tempm(1:m)=a(1:m,nm)
                      a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                   end do
                   exit
                end if
             end do
             z=w(k)
             if (l ==k) then    !Convergence.
                if (z <0.0)then !Singular value is made nonnegative.
                   w(k)=-z
                   v(1:n,k)=-v(1:n,k)
                end if
                exit
             end if
             if (its == num_its) then
                ERR_MathGen = .true.
                ERR_MathGen_Mess = " => SVDcmp_sp: convergence not reached ! "
                return
             end if
             x=w(l)             !Shift from ottom 2-y-2 minor.
             nm=k-1
             y=w(nm)
             g=rv1(nm)
             h=rv1(k)
             f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
             g=hypot(f,1.0_sp)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.0  ! Next QR transformation:
             s=1.0
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=hypot(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=hypot(f,h)
                w(j)=z                 !Rotation can e arbitrary if z =0 .
                if (abs(z) > tiny(1.0_sp) )then
                   z=1.0_sp/z
                   c=f*z
                   s=h*z
                end if
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
             end do
             rv1(l)=0.0
             rv1(k)=f
             w(k)=x
          end do
       end do

       return
    End Subroutine Svdcmp_sp

    !!----
    !!---- SUBROUTINE SWAP
    !!----
    !!----    Swap the contents of a and b, when mask (if given) is true.
    !!----
    !!----
    !!---- Update: 11/07/2015
    !!

    !!--++
    !!--++ SUBROUTINE SWAP_C
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_C(a,b)
       !---- Arguments ----!
       complex, intent(in out) :: a
       complex, intent(in out) :: b

       !---- Local variables ----!
       complex :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_C

    !!--++
    !!--++ SUBROUTINE SWAP_CM
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_Cm(a,b)
       !---- Arguments ----!
       complex, dimension(:,:), intent(in out) :: a
       complex, dimension(:,:), intent(in out) :: b

       !---- Local variables ----!
       complex, dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Cm

    !!--++
    !!--++ SUBROUTINE SWAP_CV
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_Cv(a,b)
       !---- Arguments ----!
       complex, dimension(:), intent(in out) :: a
       complex, dimension(:), intent(in out) :: b

       !---- Local variables ----!
       complex, dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Cv

    !!--++
    !!--++ SUBROUTINE SWAP_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_I(A,B)
       !---- Arguments ----!
       integer , intent(in out) :: a
       integer , intent(in out) :: b

       !---- Local variables ----!
       integer  :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_I

    !!--++
    !!--++ SUBROUTINE SWAP_IM
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_Im(A,B)
       !---- Arguments ----!
       integer, dimension(:,:), intent(in out) :: a
       integer, dimension(:,:), intent(in out) :: b

       !---- Local Variables ----!
       integer, dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Im

    !!--++
    !!--++ SUBROUTINE SWAP_IV
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_Iv(A,B)
       !---- Arguments ----!
       integer, dimension(:), intent(in out) :: a
       integer, dimension(:), intent(in out) :: b

       !---- Local Variables ----!
       integer, dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Iv

    !!--++
    !!--++ SUBROUTINE SWAP_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_R(A,B)
       !---- Arguments ----!
       real(kind=cp), intent(in out) :: a
       real(kind=cp), intent(in out) :: b

       !---- Local variables ----!
       real(kind=cp) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_R

    !!--++
    !!--++ SUBROUTINE SWAP_RM
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_Rm(A,B)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       real(kind=cp), dimension(:,:), intent(in out) :: b

       !---- Local variables ----!
       real(kind=cp), dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Rm

    !!--++
    !!--++ SUBROUTINE SWAP_RV
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Swap_Rv(A,B)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out) :: a
       real(kind=cp), dimension(:), intent(in out) :: b

       !---- Local variables ----!
       real(kind=cp), dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Rv

    !!--++
    !!--++ SUBROUTINE MASKED_SWAP_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b if mask=.true.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Masked_Swap_R(A,B,Mask)
       !---- Arguments ----!
       real(kind=cp), intent(in out) :: a
       real(kind=cp), intent(in out) :: b
       logical,           intent(in) :: mask

       !---- Local Variables ----!
       real(kind=cp) :: swp

       if (mask) then
          swp=a
          a=b
          b=swp
       end if

       return
    End Subroutine Masked_Swap_R

    !!--++
    !!--++ SUBROUTINE MASKED_SWAP_RM
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b where mask=.true.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Masked_Swap_Rm(A,B,Mask)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       real(kind=cp), dimension(:,:), intent(in out) :: b
       logical,       dimension(:,:), intent(in)     :: mask

       !---- Local variables ----!
       real(kind=cp), dimension(size(a,1),size(a,2)) :: swp

       where (mask)
          swp=a
          a=b
          b=swp
       end where

       return
    End Subroutine Masked_Swap_Rm

    !!--++
    !!--++ SUBROUTINE MASKED_SWAP_RV
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b where mask=.true.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Masked_Swap_Rv(A,B,Mask)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out) :: a
       real(kind=cp), dimension(:), intent(in out) :: b
       logical,       dimension(:), intent(in)     :: mask

       !---- Local variables ----!
       real(kind=cp), dimension(size(a))           :: swp

       where (mask)
          swp=a
          a=b
          b=swp
       end where

       return
    End Subroutine Masked_Swap_Rv

    !!--++
    !!--++ SUBROUTINE TQLI1
    !!--++
    !!--++    (PRIVATE)
    !!--++    QL-algorithm with implicit shifts, to determine the eigenvalues
    !!--++    and eigenvectors of a real tridiagonal symmetric matrix, or of
    !!--++    a real symmetric matrix previously reduced by tred. D is a vector
    !!--++    with the diagonal elements of the tridiagonal matrix. on output
    !!--++    it returns the eigenvalues. the vector e inputs the subdiagonal
    !!--++    elements of the tridiagonal matrix, with E(1) arbitrary. on
    !!--++    output e is destroyed.
    !!--++    In TLQ1 only the eigenvalues are calculated
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Tqli1(d,e,n)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out):: d, e ! d(np),e(np)
       integer,                     intent(in )   :: n

       !---- Local variables ----!
       integer      :: i, iter, l, m, mv
       real(kind=cp):: b, c, dd, f, g, p, r, s, comp

       call init_Err_MathGen()
       do i=2,n
          e(i-1)=e(i)
       end do
       e(n)=0.0
       do l=1,n
          iter=0
          do_g : do
             mv=n
             do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                comp= abs(e(m))+dd
                if (abs(comp-dd) <= ep_ss) then
                   mv=m
                   exit
                end if
             end do
             m=mv

             if (m /= l) then
                if (iter == 40) then
                   ERR_MathGen=.true.
                   ERR_MathGen_Mess=" Too many iterations in TQLI1"
                   exit
                end if

                iter=iter+1
                g=(d(l+1)-d(l))/(2.0*e(l))
                r=sqrt(g*g+1.0)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                   f=s*e(i)
                   b=c*e(i)
                   r=sqrt(f*f+g*g)
                   e(i+1)=r
                   if (abs(r)  <= ep_ss) then
                      d(i+1)=d(i+1)-p
                      e(m)=0.0
                      cycle do_g
                   end if
                   s=f/r
                   c=g/r
                   g=d(i+1)-p
                   r=(d(i)-g)*s+2.0*c*b
                   p=s*r
                   d(i+1)=g+p
                   g=c*r-b
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
                cycle do_g
             end if
             exit
          end do do_g
       end do

       return
    End Subroutine Tqli1

    !!--++
    !!--++ SUBROUTINE TQLI2
    !!--++
    !!--++    (PRIVATE)
    !!--++    QL-algorithm with implicit shifts, to determine the eigenvalues
    !!--++    and eigenvectors of a real tridiagonal symmetric matrix, or of
    !!--++    a real symmetric matrix previously reduced by tred. D is a vector
    !!--++    with the diagonal elements of the tridiagonal matrix. on output
    !!--++    it returns the eigenvalues. the vector e inputs the subdiagonal
    !!--++    elements of the tridiagonal matrix, with E(1) arbitrary. on
    !!--++    output e is destroyed.
    !!--++    The eigenvectors of the tridiagonal matrix are calculated in TLQ2
    !!--++    by providing the matrix Z  as the identity matrix on input. if the
    !!--++    eigenvectors of the matrix reduced by tred are required, then Z
    !!--++    is input as the matrix output of tred. in either cased, the k-th
    !!--++    column of Z returns the mormalized eigenvector corresponding to
    !!--++    D(k).
    !!--++
    !!--++  Update: 11/07/2015
    !!
    Subroutine Tqli2(d,e,n,z)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)
       integer,                       intent(in )    :: n
       real(kind=cp), dimension(:,:), intent(in out) :: z    ! z(np,np)

       !---- Local Variables ----!
       integer       :: i, iter, k, l, m, mv
       real(kind=cp) :: b, c, dd, f, g, p, r, s, comp

       call init_Err_MathGen()
       do i=2,n
          e(i-1)=e(i)
       end do

       e(n)=0.0
       do l=1,n
          iter=0
          do_g: do
             mv=n
             do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                comp= abs(e(m))+dd
                if (abs(comp-dd) <= ep_ss) then
                   mv=m
                   exit
                end if
             end do
             m=mv
             if (m /= l) then
                if (iter == 40) then
                   ERR_MathGen=.true.
                   ERR_MathGen_Mess=" Too many iterations in TQLI2"
                   exit
                end if

                iter=iter+1
                g=(d(l+1)-d(l))/(2.0*e(l))
                r=sqrt(g*g+1.0)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                   f=s*e(i)
                   b=c*e(i)
                   r=sqrt(f*f+g*g)
                   e(i+1)=r
                   if (abs(r) <= ep_ss) then
                      d(i+1)=d(i+1)-p
                      e(m)=0.0
                      cycle do_g
                   end if
                   s=f/r
                   c=g/r
                   g=d(i+1)-p
                   r=(d(i)-g)*s+2.0*c*b
                   p=s*r
                   d(i+1)=g+p
                   g=c*r-b

                   !---- omit lines from here ...
                   do k=1,n
                      f=z(k,i+1)
                      z(k,i+1)=s*z(k,i)+c*f
                      z(k,i)=c*z(k,i)-s*f
                   end do

                   !---- ... to here when finding only eigenvalues.
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
                cycle do_g
             end if
             exit
          end do do_g
       end do

       return
    End Subroutine Tqli2

    !!--++
    !!--++ SUBROUTINE TRED1
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for preparing the matrix to find only eigenvalues
    !!--++    Householder reduction of a real symetric nxn matrix A.
    !!--++    On output A is replaced by the orthogonal matrix Q effecting
    !!--++    the transformation. D returns the diagonal elements of the tri-
    !!--++    diagonal matrix and E the off-diagonal elements with E(1)=0.
    !!--++    In tred1 several lines have been deleted and A contains no
    !!--++    useful information on output.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Tred1(a,n,d,e)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)
       integer,                       intent(in)     :: n
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)

       !---- Local Variables ----!
       integer :: i, j, k, l
       real(kind=cp)    :: f, g, h, hh, scala

       do i=n,2,-1
          l=i-1
          h=0.0
          scala=0.0
          if (l > 1)then
             do k=1,l
                scala=scala+abs(a(i,k))
             end do
             if (abs(scala) <= ep_ss) then
                e(i)=a(i,l)
             else
                do k=1,l
                   a(i,k)=a(i,k)/scala
                   h=h+a(i,k)**2
                end do
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scala*g
                h=h-f*g
                a(i,l)=f-g
                f=0.0
                do j=1,l
                   g=0.0
                   do k=1,j
                      g=g+a(j,k)*a(i,k)
                   end do
                   do k=j+1,l
                      g=g+a(k,j)*a(i,k)
                   end do
                   e(j)=g/h
                   f=f+e(j)*a(i,j)
                end do
                hh=f/(h+h)
                do j=1,l
                   f=a(i,j)
                   g=e(j)-hh*f
                   e(j)=g
                   do k=1,j
                      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                   end do
                end do
             end if
          else
             e(i)=a(i,l)
          end if
          d(i)=h
       end do

       e(1)=0.0
       do i=1,n
          d(i)=a(i,i)
       end do

       return
    End Subroutine Tred1

    !!--++
    !!--++ SUBROUTINE TRED2
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for preparing the matrix to find the complete set
    !!--++    of eigenvectors.
    !!--++    Householder reduction of a real symetric nxn matrix A.
    !!--++    On output A is replaced by the orthogonal matrix Q effecting
    !!--++    the transformation. D returns the diagonal elements of the tri-
    !!--++    diagonal matrix and E the off-diagonal elements with E(1)=0.
    !!--++
    !!--++ Update: 11/07/2015
    !!
    Subroutine Tred2(a,n,d,e)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)
       integer,                       intent(in)     :: n
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)

       !---- Local variables ----!
       integer :: i, j, k, l
       real(kind=cp)    :: f, g, h, hh, scala

       do i=n,2,-1
          l=i-1
          h=0.0
          scala=0.0
          if (l > 1)then
             do k=1,l
                scala=scala+abs(a(i,k))
             end do
             if (abs(scala) <= ep_ss) then
                e(i)=a(i,l)
             else
                do k=1,l
                   a(i,k)=a(i,k)/scala
                   h=h+a(i,k)**2
                end do
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scala*g
                h=h-f*g
                a(i,l)=f-g
                f=0.0
                do j=1,l
                   !---- omit following line if finding only eigenvalues
                   a(j,i)=a(i,j)/h
                   g=0.0
                   do k=1,j
                      g=g+a(j,k)*a(i,k)
                   end do
                   do k=j+1,l
                      g=g+a(k,j)*a(i,k)
                   end do
                   e(j)=g/h
                   f=f+e(j)*a(i,j)
                end do
               hh=f/(h+h)
                do j=1,l
                   f=a(i,j)
                   g=e(j)-hh*f
                   e(j)=g
                   do k=1,j
                      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                   end do
                end do
             end if
          else
             e(i)=a(i,l)
          end if
          d(i)=h
       end do

       !---- omit following line if finding only eigenvalues.
       d(1)=0.0
       e(1)=0.0
       do i=1,n
          !---- delete lines from here ...
          l=i-1
          if (abs(d(i)) > ep_ss)then
             do j=1,l
                g=0.0
                do k=1,l
                   g=g+a(i,k)*a(k,j)
                end do
                do k=1,l
                   a(k,j)=a(k,j)-g*a(k,i)
                end do
             end do
          end if
          !---- ... to here when finding only eigenvalues.
          d(i)=a(i,i)
          !---- also delete lines from here ...
          a(i,i)=1.0
          do j=1,l
             a(i,j)=0.0
             a(j,i)=0.0
          end do
          !---- ... to here when finding only eigenvalues.
       end do

       return
    End Subroutine Tred2

 End Module CFML_Math_General

