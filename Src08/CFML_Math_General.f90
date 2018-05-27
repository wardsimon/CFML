!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2018  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_Math_General
!!----         Mathematic general utilities for use in Crystallography and
!!----         Solid State Physics and Chemistry.
!!----
!!
 Module CFML_Math_General
    !---- Use Modules ----!
    Use CFML_GlobalDeps
    
    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public :: Factorial, Factorial_SP, Factorial_DP, Co_Prime, Euclidean_Norm, &
              Pgcd, Ppcm, Modulo_Lat, Poly_Legendre

    !---- List of public overloaded procedures: functions ----!
    public :: Co_Linear, Debye, Negligible, Equal_Matrix, Equal_Vector, In_limits, &
              Locate, Lower_Triangular,  Outerprod, Trace, Zbelong,  Norm,  &
              Scalar, Upper_Triangular, Erfc_Deriv

    !---- List of public subroutines ----!
    public ::  Invert_Matrix, LU_Decomp, LU_Backsub, Matinv,        &
               Sort_Strings, Spline, Splint, Set_Epsg,In_Sort,      &
               First_Derivative, Second_Derivative, SmoothingVec, Points_in_Line2D,   &
               Co_Prime_vector, Median_QS, Sph_Jn

    !---- List of public overloaded procedures: subroutines ----!
    public ::  Determinant, Diagonalize_Sh, Linear_Dependent, Rank, Sort,   &
               Svdcmp, Swap


    !---- Parameters ----!
    integer, dimension(1000), parameter, public :: PRIMES =                               & ! List of the first 1000 prime numbers.
            [ 2,      3,      5,      7,     11,     13,     17,     19,     23,     29,  &
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
           7841,   7853,   7867,   7873,   7877,   7879,   7883,   7901,   7907,   7919 ]

    real(kind=cp), parameter :: EP_SS=1.0E-12_cp  ! Internal epsilon value used for comparison in matrix operations

    !---- Variables ----!
    real(kind=cp)            :: epss= 1.0E-5_cp   ! Internal epsilon value used for comparing reals to integers

    !---- Interfaces - Overloaded ----!
    Interface  Negligible
       Module Procedure Negligibler
       Module Procedure Negligiblec
    End Interface

    Interface  Co_Linear
       Module Procedure Co_linear_C
       Module Procedure Co_linear_I
       Module Procedure Co_linear_R
    End Interface

    Interface  Debye
       Module Procedure Debye_DP
       Module Procedure Debye_SP
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
        Module Procedure swap_i
        Module Procedure swap_r
        Module Procedure swap_masked_r
    End interface
    
    Interface Erfc_Deriv
        Module Procedure Erfc_Deriv_dp
        Module Procedure Erfc_Deriv_sp
    End interface

 
    !> Interface Zone 
    Interface
 
       Module Elemental Function Factorial(N) Result(Fact)    
          !---- Argument ----!       
          integer, intent(in) :: N        ! Factorial of N       
          integer             :: Fact       
       End Function Factorial    
 
       Module Elemental Function Factorial_SP(N) Result(Fact)    
          !---- Arguments ----!       
          integer, intent(in) :: N    ! Factorial of N       
          real(kind=sp)       :: Fact       
       End Function Factorial_SP    
 
       Module Elemental Function Factorial_DP(N) Result(Fact)    
          !---- Arguments ----!       
          integer, intent(in) :: N      ! Factorial of N       
          real(kind=dp)       :: Fact       
       End Function Factorial_DP    
 
       Module Elemental Function Negligiblec(v) Result(Neglig)    
          !---- Argument ----!       
          complex, intent( in) :: v         ! Complex number       
          logical              :: Neglig       
       End Function Negligiblec    
 
       Module Elemental Function Negligibler(v) Result(neglig)    
          !---- Argument ----!       
          real(kind=cp), intent( in) :: v          ! Real number       
          logical                    :: Neglig       
       End Function Negligibler    
 
       Module Elemental Function Poly_Legendre(l,m,x) result(plmx)    
          !---- Arguments ----!       
          integer,      intent (in) :: l       
          integer,      intent (in) :: m       
          real(kind=cp),intent (in) :: x       
          real(kind=cp)             :: plmx       
       End Function Poly_Legendre    
 
       Module Function Debye_PR_ChebyshevSeries(n, a, t) Result(fval)    
          !---- Arguments ----!       
          integer,                       intent(in) :: N    ! The no. of terms in the sequence       
          real(kind=dp), dimension(0:N), intent(in) :: A    ! The coefficients of the Chebyshev series       
          real(kind=dp),                 intent(in) :: T    ! The value at which the series is to be evaluated       
          real(kind=dp)                             :: fval ! Return value       
       End Function Debye_PR_ChebyshevSeries    
 
       Module Function Co_linear_C(a,b,n) Result(co_linear)    
          !---- Argument ----!       
          complex, dimension(:), intent(in)           :: a,b    ! Complex vectors       
          integer,               optional, intent(in) :: n      ! Dimension of the vectors       
          logical                                     :: co_linear       
       End Function Co_linear_C    
 
       Module Function Co_linear_I(a,b,n) Result(co_linear)    
          !---- Argument ----!       
          integer, dimension(:),           intent(in) :: a,b        ! Input vectors       
          integer,               optional, intent(in) :: n          ! Dimension of the vector       
          logical                                     :: co_linear       
       End Function Co_linear_I    
 
       Module Function Co_linear_R(a,b,n) Result(co_linear)    
          !---- Argument ----!       
          real(kind=cp), dimension(:),           intent(in) :: a,b        ! Input real vectors       
          integer,                     optional, intent(in) :: n          ! Dimension of the vectors       
          logical                                           :: co_linear       
       End Function Co_linear_R    
 
       Module Function Co_Prime(v,imax) result(cop)    
          !---- Arguments ----!       
          integer, dimension(:),           intent(in) :: v          ! Input vector of numbers       
          integer,               optional, intent(in) :: imax       ! Maximun prime number to be tested       
          Logical                                     :: cop       
       End Function Co_Prime    
 
       Module Function Debye_DP(N,X) Result(Fval)    
          !---- Arguments ----!       
          integer,       intent(in) :: N ! Order of the Debye function       
          real(kind=dp), intent(in) :: X ! Value       
          real(kind=dp)             :: fval       
       End Function Debye_DP    
 
       Module Function Debye_SP(N,X) Result(Fval)    
          !---- Arguments ----!       
          integer,       intent(in) :: N ! Order of the Debye function       
          real(kind=sp), intent(in) :: X ! Value       
          real(kind=sp)             :: fval       
       End Function Debye_SP    
 
       Module Function Debye1(X) Result(Fval)    
          !---- Arguments ----!       
          real(kind=dp), intent(in) :: X       
          real(kind=dp)             :: fval       
       End Function Debye1    
 
       Module Function Debye2(X) Result(Fval)    
          !---- Argument ----!       
          real(kind=dp), intent(in) :: X       
          real(kind=dp)             :: fval       
       End Function Debye2    
 
       Module Function Debye3(X) Result(Fval)    
          !---- Argument ----!       
          real(kind=dp), intent(in) :: X       
          real(kind=dp)             :: fval       
       End Function Debye3    
 
       Module Function Debye4(X) Result(FVal)    
          !---- Argument ----!       
          real(kind=dp), intent(in) :: X       
          real(kind=dp)             :: fval       
       End Function Debye4    
 
       Module Function DebyeN(n,x) Result(Fval)    
          !---- Arguments ----!       
          integer,       intent(in) :: N ! Order of Debye function       
          real(kind=dp), intent(in) :: X       
          real(kind=dp)             :: Fval       
       End Function DebyeN    
 
       Module Pure Function Sphjn_PR_Envj(N,X) Result(Y)    
          !---- Arguments ----!       
          integer,       intent(in) :: n       
          real(Kind=dp), intent(in) :: x       
          real(Kind=dp)             :: y       
       End Function Sphjn_PR_Envj    
 
       Module Function Equal_Matrix_I(a,b,n) result(info)    
          !---- Argument ----!       
          integer, dimension(:,:),          intent(in) :: a,b     ! Input arrays (NxN)       
          integer,                optional, intent(in) :: n       ! Dimension of Arrays       
          logical                                      :: info       
       End Function Equal_Matrix_I    
 
       Module Function Equal_Matrix_R(a,b,n) result(info)    
          !---- Argument ----!       
          real(kind=cp), dimension(:,:),          intent(in) :: a,b      ! Input arrays NxN       
          integer,                      optional, intent(in) :: n        ! Dimensions N       
          logical                                            :: info       
       End Function Equal_Matrix_R    
 
       Module Function Equal_Vector_I(a,b,n) result(info)    
          !---- Argument ----!       
          integer, dimension(:),           intent(in) :: a,b    ! Input vectors       
          integer,               optional, intent(in) :: n      ! Dimension of the vectors       
          logical                                     :: info       
       End Function Equal_Vector_I    
 
       Module Function Equal_Vector_R(a,b,n) result(info)    
          !---- Argument ----!       
          real(kind=cp), dimension(:),           intent(in) :: a,b      ! Input vectors       
          integer,                     optional, intent(in) :: n        ! Dimension of the vector       
          logical                                           :: info       
       End Function Equal_Vector_R    
 
       Module Function Euclidean_Norm(x,n) Result(Fn_Val)    
          !---- Arguments ----!       
          Real (Kind=cp), Dimension(:),           intent(In)  :: x      ! Input vector       
          Integer,                      optional, intent(In)  :: n      ! Dimension of the vector       
          Real (Kind=cp)                                      :: Fn_Val ! Return value       
       End Function Euclidean_Norm    
 
       Module Function in_limits_dp(v,limits,n) result(ok)    
          !---- Arguments ----!       
          real(kind=dp), dimension(:),             intent(in) :: v        ! Input vector       
          real(kind=dp), dimension(:,:),           intent(in) :: limits   ! Normally (2,n)       
          integer,                       optional, intent(in) :: n        ! Dimension of Vect       
          logical                                             :: ok       
       End Function in_limits_dp    
 
       Module Function in_limits_int(v,limits,n) result(ok)    
          !---- Arguments ----!       
          integer, dimension(:),             intent(in) :: v        ! Input Vector       
          integer, dimension(:,:),           intent(in) :: limits   ! Normally (2,n)       
          integer,                 optional, intent(in) :: n        ! Dimension of vect       
          logical                                       :: ok       
       End Function in_limits_int    
 
       Module Function in_limits_sp(v,limits,n) result(ok)    
          !---- Arguments ----!       
          real(kind=sp), dimension(:),             intent(in) :: v        ! Input vector       
          real(kind=sp), dimension(:,:),           intent(in) :: limits   ! Normally (2,n)       
          integer,                       optional, intent(in) :: n        ! Dimension of the vector V       
          logical                                             :: ok       
       End Function in_limits_sp    
 
       Module Elemental Function Pgcd(a,b) Result(mcd)    
          !---- Arguments ----!       
          integer, intent(in) :: a,b   ! Integers       
          integer             :: mcd   ! Maximum common divisor       
       End Function Pgcd    
 
       Module Elemental Function Ppcm(a,b) result(mcm)    
          !---- Arguments ----!       
          integer, intent(in) :: a,b    ! Integers       
          integer             :: mcm    ! Minimum common multiple       
       End Function Ppcm    
 
       Module Function Locate_I(V,x) Result(j)    
          !---- Argument ----!       
          integer, dimension(:), intent(in):: v  ! Input vector       
          integer,               intent(in):: x  ! Value       
          integer                          :: j       
       End Function Locate_I    
 
       Module Function Locate_R(V,x) Result(j)    
          !---- Argument ----!       
          real(kind=cp), dimension(:), intent(in):: v       
          real(kind=cp),               intent(in):: x       
          integer                                :: j       
       End Function Locate_R    
 
       Module Function Lower_Triangular_I(A,n) Result (T)    
          !---- Argument ----!       
          integer, dimension(:,:), intent(in) :: A    ! Input array       
          integer,                 intent(in) :: n    ! Dimension of array       
          integer, dimension(n,n)             :: T       
       End Function Lower_Triangular_I    
 
       Module Function Lower_Triangular_R(A,n) Result (T)    
          !---- Argument ----!       
          real(kind=cp), dimension(:,:), intent(in) :: A    ! Input Array       
          integer,                       intent(in) :: n    ! Dimension of A       
          real(kind=cp), dimension(n,n)             :: T       
       End Function Lower_Triangular_R    
 
       Module Function Modulo_Lat(v) result(u)    
          !---- Argument ----!       
          real(kind=cp), dimension(:), intent( in) :: v       
          real(kind=cp), dimension(1:size(v))      :: u       
       End Function Modulo_Lat    
 
       Module Function Norm_I(X,G) Result(R)    
          !---- Arguments ----!       
          integer,       dimension(:),   intent(in) :: x    ! Input vector       
          real(kind=cp), dimension(:,:), intent(in) :: g    ! Metric array       
          real(kind=cp)                             :: r    ! Norm of the input vector       
       End Function Norm_I    
 
       Module Function Norm_R(X,G) Result(R)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:),   intent(in) :: x   ! Input vector       
          real(kind=cp), dimension(:,:), intent(in) :: g   ! Metrics       
          real(kind=cp)                             :: r   ! Norm of the vector       
       End Function Norm_R    
 
       Module Function Outerprod_dp(a,b)  Result(c)    
          !---- Arguments ----!       
          real(kind=dp),dimension(:),intent(in)    :: a,b       
          real(kind=dp),dimension(size(a),size(b)) :: c       
       End Function Outerprod_dp    
 
       Module Function Outerprod_sp(a,b)  Result(c)    
          !---- Arguments ----!       
          real(kind=sp),dimension(:),intent(in)    :: a,b       
          real(kind=sp),dimension(size(a),size(b)) :: c       
       End Function Outerprod_sp    
 
       Module Function Scalar_I(X,Y,G) Result(R)    
          !---- Arguments ----!       
          integer, dimension(:),         intent(in) :: x,y     ! Input vectors       
          real(kind=cp), dimension(:,:), intent(in) :: g       ! Metrics       
          real(kind=cp)                             :: r       ! Scalar       
       End Function Scalar_I    
 
       Module Function Scalar_R(X,Y,G) Result(R)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:),   intent(in) :: x,y    ! Input vectors       
          real(kind=cp), dimension(:,:), intent(in) :: g      ! Metrics       
          real(kind=cp)                             :: r      ! Scalar       
       End Function Scalar_R    
 
       Module Pure Function Sphjn_PR_Start1(X,Mp) Result (Start)    
          !---- Arguments ----!       
          real(kind=dp), intent(in) :: x         ! Argument of Jn(x)       
          integer, intent(in)       :: mp        ! Value of magnitude       
          integer                   :: start       
       End Function Sphjn_PR_Start1    
 
       Module Pure Function Sphjn_PR_Start2(X,N,Mp) Result(Start)    
          !---- Arguments ----!       
          real(kind=dp), intent(in) :: x     ! Argument of Jn(x)       
          integer,       intent(in) :: n     ! Order of Jn(x)       
          integer,       intent(in) :: mp    ! Significant digit       
          integer                   :: start       
       End Function Sphjn_PR_Start2    
 
       Module Function Trace_C(a) Result(b)    
          !---- Argument ----!       
          complex, dimension(:,:), intent(in) :: a       
          complex                             :: b       
       End Function Trace_C    
 
       Module Function Trace_I(a) Result(b)    
          !---- Argument ----!       
          integer, dimension(:,:), intent(in) :: a       
          integer                             :: b       
       End Function Trace_I    
 
       Module Function Trace_R(a) Result(b)    
          !---- Argument ----!       
          real(kind=cp), dimension(:,:), intent(in) :: a       
          real(kind=cp)                             :: b       
       End Function Trace_R    
 
       Module Function Upper_Triangular_I(A,n) Result (T)    
          !---- Argument ----!       
          integer, dimension(:,:), intent(in) :: A     ! Input array       
          integer,                 intent(in) :: n     ! Dimension       
          integer, dimension(n,n)             :: T       
       End Function Upper_Triangular_I    
 
       Module Function Upper_Triangular_R(A,n) Result (T)    
          !---- Argument ----!       
          real(kind=cp), dimension(:,:), intent(in) :: A   ! Input array       
          integer,                       intent(in) :: n   ! Dimension       
          real(kind=cp), dimension(n,n)             :: T       
       End Function  Upper_Triangular_R    
 
       Module Function ZbelongM(v) Result(belong)    
          !---- Argument ----!       
          real(kind=cp),   dimension(:,:), intent( in) :: v        ! Input array       
          logical                                      :: belong       
       End Function ZbelongM    
 
       Module Function ZbelongN(a) Result(belong)    
          !---- Argument ----!       
          real(kind=cp), intent( in) :: a              ! Input number       
          logical                    :: belong       
       End Function ZbelongN    
 
       Module Function ZbelongV(v) Result(belong)    
          !---- Argument ----!       
          real(kind=cp),   dimension(:), intent( in) :: v      ! Input vector       
          logical                                    :: belong       
       End Function ZbelongV  
       
       Module Function Erfc_Deriv_DP(X) Result(Der)
          !---- Argument ----!
          real(kind=dp), intent(in)    :: x
          real(kind=dp)                :: der
       End Function Erfc_Deriv_DP 
       
       Module Function Erfc_Deriv_SP(X) Result(Der)
          !---- Argument ----!
          real(kind=sp), intent(in)    :: x
          real(kind=sp)                :: der
       End Function Erfc_Deriv_SP 
 
       Module Pure Subroutine Median_QS(x, n, xmed)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:), intent(in out) :: x      ! In: Vector  Out: Sorted vector       
          integer,                     intent(in)     :: n      ! Number of data in X       
          real(kind=cp),               intent(out)    :: xmed   ! Media of consided data       
       End Subroutine Median_QS    
 
       Module Subroutine Co_Prime_Vector(V,Cop,f)    
          !---- Arguments ----!       
          integer, dimension(:),           intent(in)  :: v              !input integer vector       
          integer, dimension(:),           intent(out) :: cop            !Output co-prime vector       
          integer,               optional, intent(out) :: f              !Common multiplicative factor       
       End Subroutine Co_Prime_vector    
 
       Module Subroutine Determinant_C(A,n,determ)    
          !---- Arguments ----!       
          complex, dimension(:,:), intent( in) :: A         !input square matrix (n,n)       
          integer,                 intent( in) :: n         !actual dimension of A       
          real(kind=cp),           intent(out) :: determ    !det(A) if real and det(AR)^2 + det(AI)^2 if complex       
       End Subroutine Determinant_C    
 
       Module Subroutine Determinant_R(A,n,determ)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:,:), intent( in) :: A      ! Input array NxN       
          integer,                       intent( in) :: n      ! Dimension of A       
          real(kind=cp),                 intent(out) :: determ ! Value       
       End Subroutine Determinant_R    
 
       Module Subroutine Diagonalize_Herm(a,n,e_val,e_vect,norder)    
          !---- Arguments ----!       
          complex,           dimension(:,:), intent( in)  :: A       
          integer,                           intent( in)  :: n       
          real(kind=cp),     dimension(:),   intent(out)  :: E_val   ! Eigenvalues       
          complex, optional, dimension(:,:), intent(out)  :: E_vect  ! Eigenvectors       
          logical, optional,                 intent(in)   :: norder  ! If present no ordering       
       End Subroutine Diagonalize_Herm    
 
       Module Subroutine Diagonalize_Symm(A,n,E_Val,E_vect,norder)    
          !---- Arguments ----!       
          real(kind=cp),           dimension(:,:), intent( in)  :: A       
          integer,                                 intent( in)  :: n       
          real(kind=cp),           dimension(:),   intent(out)  :: E_val    ! Eigenvalues       
          real(kind=cp), optional, dimension(:,:), intent(out)  :: E_vect   ! Eigenvectors       
          logical,       optional,                 intent(in)   :: norder   ! If present no ordering       
       End Subroutine Diagonalize_Symm    
 
       Module Subroutine Diagonalize_EigenvSort(d,v,n,io)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:),   intent(in out) :: d       
          real(kind=cp), dimension(:,:), intent(in out) :: v       
          integer,                       intent(in)     :: n       
          integer,                       intent(in)     :: io       
       End Subroutine Diagonalize_EigenvSort    
 
       Module Pure Subroutine First_Derivative(x,y,n,d2y,d1y)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:), intent(in)  :: x     ! Input vector       
          real(kind=cp), dimension(:), intent(in)  :: y     ! Yi=F(xi)       
          integer ,                    intent(in)  :: n     ! Dimension of X and Y       
          real(kind=cp), dimension(:), intent(in)  :: d2y   ! Vector containing the second derivative in each point       
          real(kind=cp), dimension(:), intent(out) :: d1y   ! Vector containing the first derivative in each point       
       End Subroutine First_Derivative    
 
       Module Subroutine In_Sort(id,n,p,q)    
          !---- Arguments ----!       
          integer, dimension(:), intent(in) :: id  !Integer array to be sorted       
          integer,               intent(in) :: n   !Number items in the array       
          integer, dimension(:), intent(in) :: p   !Initial pointer from a previous related call       
          integer, dimension(:), intent(out):: q   !Final pointer doing the sort of id       
       End Subroutine In_Sort    
 
       Module Subroutine Invert_Matrix(a,b,singular,perm)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:,:),  intent(in ) :: a         ! Input Array       
          real(kind=cp), dimension(:,:),  intent(out) :: b         ! Output array       
          logical,                        intent(out) :: singular       
          integer, dimension(:),optional, intent(out) :: perm       
       End Subroutine Invert_Matrix    
 
       Module Subroutine Linear_DependentC(A,na,B,nb,mb,info)    
          !---- Arguments ----!       
          complex, dimension(:),   intent(in)  :: a       
          complex, dimension(:,:), intent(in)  :: b       
          integer,                 intent(in)  :: na,nb,mb       
          logical,                 intent(out) :: info       
       End Subroutine Linear_DependentC    
 
       Module Subroutine Linear_DependentI(A,na,B,nb,mb,info)    
          !---- Arguments ----!       
          integer, dimension(:),   intent(in)  :: a       
          integer, dimension(:,:), intent(in)  :: b       
          integer,                 intent(in)  :: na,nb,mb       
          logical,                 intent(out) :: info       
       End Subroutine Linear_DependentI    
 
       Module Subroutine Linear_DependentR(A,na,B,nb,mb,info)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:),   intent(in)  :: a       
          real(kind=cp), dimension(:,:), intent(in)  :: b       
          integer,                       intent(in)  :: na,nb,mb       
          logical,                       intent(out) :: info       
       End Subroutine Linear_DependentR    
 
       Module Subroutine LU_Backsub(a,indx,b)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:,:), intent(in)     :: a       
          integer,         dimension(:), intent(in)     :: indx       
          real(kind=cp),   dimension(:), intent(in out) :: b       
       End Subroutine LU_Backsub    
 
       Module Subroutine LU_Decomp(a,d,singular,indx)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:,:), intent(in out) :: a       
          real(kind=cp),                 intent(out)    :: d       
          logical,                       intent(out)    :: singular       
          integer,  dimension(:), intent(out), optional :: indx       
       End Subroutine LU_Decomp    
 
       Module Subroutine Matinv(a,n)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:,:), intent(in out) :: a       
          integer     ,                  intent(in)     :: n       
       End Subroutine Matinv    
 
       Module Subroutine Sort_PR_Partition(A, Marker)    
          !---- Arguments ----!       
          character(len=*), dimension(:), intent(in out) :: A       
          integer,                        intent(   out) :: marker       
       End Subroutine Sort_PR_Partition    
 
       Module Subroutine Points_In_Line2D(X1, XN, N, XP)    
          !---- Arguments ----!       
          real(kind=cp), dimension(2),   intent(in)  :: X1   ! Point1 in 2D       
          real(kind=cp), dimension(2),   intent(in)  :: XN   ! PointN in 2D       
          integer,                       intent(in)  :: N    ! Number of Total points       
          real(kind=cp), dimension(:,:), intent(out) :: XP   ! List of points       
       End Subroutine Points_In_Line2D    
 
       Module Subroutine Rank_dp(a,tol,r)    
          !---- Arguments ----!       
          real(kind=dp), dimension(:,:),intent( in)      :: a     ! Input array       
          real(kind=dp),                intent( in)      :: tol   ! Tolerance       
          integer,                      intent(out)      :: r       
       End Subroutine Rank_dp    
 
       Module Subroutine Rank_sp(a,tol,r)    
          !---- Arguments ----!       
          real(kind=sp), dimension(:,:),intent( in)      :: a       
          real(kind=sp),                intent( in)      :: tol       
          integer,                      intent(out)      :: r       
       End Subroutine Rank_sp    
 
       Module Pure Subroutine Second_Derivative(x,y,n,d2y)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:), intent(in)  :: x    ! Input X vector       
          real(kind=cp), dimension(:), intent(in)  :: y    ! Yi=F(xi)       
          integer ,                    intent(in)  :: n    ! Number of Points       
          real(kind=cp), dimension(:), intent(out) :: d2y  ! Second derivative       
       End Subroutine Second_Derivative    
 
       Module Subroutine SmoothingVec(Y, N, Niter, Ys)    
          !---- Arguments ----!       
          real(kind=cp),dimension(:),            intent(in out) :: Y         !  In Out-> Array to be smoothed       
          integer,                               intent(in)     :: n         !  In -> Number of points       
          integer,                               intent(in)     :: niter     !  In -> Number of iterations       
          real(kind=cp),dimension(:), optional,  intent(out)    :: Ys        !  Out-> Array smoothed       
       End Subroutine SmoothingVec    
 
       Module Subroutine Sort_I(arr,n,indx)    
          !---- Arguments ----!       
          integer, dimension(:), intent(in ) :: arr       ! Vector       
          integer              , intent(in ) :: n         ! Dimension       
          integer, dimension(:), intent(out) :: indx      ! Index       
       End Subroutine Sort_I    
 
       Module Subroutine Sort_R(arr,n,indx)    
          !---- Arguments ----!       
          real(kind=cp),dimension(:), intent(in) :: arr       ! Input Vector       
          integer,                    intent(in) :: n         ! Dimension       
          integer,      dimension(:), intent(out):: indx       
       End Subroutine Sort_R    
 
       Module Recursive Subroutine Sort_Strings(Str)    
          !---- Argument ----!       
          character(len=*), dimension(:), intent(in out) :: Str       
       End Subroutine Sort_Strings    
 
       Module Pure Subroutine Spline(x,y,n,yp1,ypn,y2)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:), intent(in)  :: x               !  In -> Array X       
          real(kind=cp), dimension(:), intent(in)  :: y               !  In -> Array Yi=F(Xi)       
          integer ,                    intent(in)  :: n               !  In -> Dimension of X, Y       
          real(kind=cp),               intent(in)  :: yp1             !  In -> Derivate of Point 1       
          real(kind=cp),               intent(in)  :: ypn             !  In -> Derivate of Point N       
          real(kind=cp), dimension(:), intent(out) :: y2              ! Out -> array containing second derivatives       
       End Subroutine Spline    
 
       Module Pure Function Splint(xa,ya,y2a,n,x) Result(y)   
          !---- Arguments ----!       
          real(kind=cp), dimension(:), intent(in)  :: xa          !  In -> Array X       
          real(kind=cp), dimension(:), intent(in)  :: ya          !  In -> Array Y=F(X)       
          real(kind=cp), dimension(:), intent(in)  :: y2a         !  In -> Array Second Derivatives in X       
          integer ,                    intent(in)  :: n           !  In -> Dimension of XA,YA,Y2A       
          real(kind=cp),               intent(in)  :: x           !  In -> Point to evaluate       
          real(kind=cp)                            :: y           ! Out -> Interpoled value       
       End Function Splint    
 
       Module Subroutine Sph_Jn(n,x,nm,jn,djn)    
          !---- Arguments ----!       
          integer,                       intent(in)  :: n   !Order of jn(x) (n=0,1,2,3,...)       
          real(kind=dp),                 intent(in)  :: x   !Argument of jn(x)       
          integer,                       intent(out) :: nm  !Highest order computed       
          real(kind=dp), dimension(0:n), intent(out) :: jn  !array with spherical Bessel functions jn(x)       
          real(kind=dp), dimension(0:n), intent(out) :: djn !array with derivatives jn'(x)       
       End Subroutine Sph_Jn    
 
       Module Subroutine Svdcmp_dp(a,w,v)    
          !---- Arguments ----!       
          real(kind=dp),dimension(:,:),intent(in out) ::a   ! A(m,n)       
          real(kind=dp),dimension(:),  intent(   out) ::w   ! W(n)       
          real(kind=dp),dimension(:,:),intent(   out) ::v   ! V(n,n)       
       End Subroutine Svdcmp_dp    
 
       Module Subroutine Svdcmp_sp(a,w,v)    
          !---- Arguments ----!       
          real(kind=sp),dimension(:,:),intent(in out) :: a  ! A(M,n)       
          real(kind=sp),dimension(:),  intent(   out) :: w  ! W(n)       
          real(kind=sp),dimension(:,:),intent(   out) :: v  ! V(n,n)       
       End Subroutine Svdcmp_sp    
 
       Module Elemental Subroutine Swap_C(a,b)    
          !---- Arguments ----!       
          complex, intent(in out) :: a       
          complex, intent(in out) :: b       
       End Subroutine Swap_C    
 
       Module Elemental Subroutine Swap_I(A,B)    
          !---- Arguments ----!       
          integer , intent(in out) :: a       
          integer , intent(in out) :: b       
       End Subroutine Swap_I    
 
       Module Elemental Subroutine Swap_R(A,B)    
          !---- Arguments ----!       
          real(kind=cp), intent(in out) :: a       
          real(kind=cp), intent(in out) :: b       
       End Subroutine Swap_R    
 
       Module Elemental Subroutine Swap_Masked_R(A,B,Mask)    
          !---- Arguments ----!       
          real(kind=cp), intent(in out) :: a       
          real(kind=cp), intent(in out) :: b       
          logical,           intent(in) :: mask       
       End Subroutine Swap_Masked_R    
 
       Module Subroutine Diagonalize_PR_Tqli1(d,e,n)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:), intent(in out):: d, e ! d(np),e(np)       
          integer,                     intent(in )   :: n       
       End Subroutine Diagonalize_PR_Tqli1    
 
       Module Subroutine Diagonalize_PR_Tqli2(d,e,n,z)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)       
          integer,                       intent(in )    :: n       
          real(kind=cp), dimension(:,:), intent(in out) :: z    ! z(np,np)       
       End Subroutine Diagonalize_PR_Tqli2    
 
       Module Subroutine Diagonalize_PR_Tred1(a,n,d,e)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)       
          integer,                       intent(in)     :: n       
          real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)       
       End Subroutine Diagonalize_PR_Tred1    
 
       Module Subroutine Diagonalize_PR_Tred2(a,n,d,e)    
          !---- Arguments ----!       
          real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)       
          integer,                       intent(in)     :: n       
          real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)       
       End Subroutine Diagonalize_PR_Tred2    
 
    End Interface
    
 Contains
 
    !!---- SUBROUTINE SET_EPSG
    !!----
    !!----    Sets/Modify global EPSS.
    !!----    Calling without arguments set to default value
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Set_Epsg(Neweps)    
       !---- Arguments ----!
       real(kind=cp), optional, intent( in) :: neweps

       if (present(neweps)) then
          epss=neweps
       else
          epss=1.0E-5_cp
       end if

       return
    End Subroutine Set_Epsg   
 
End Module CFML_Math_General 
