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
!!---- MODULE: CFML_Maths
!!----         Mathematic general utilities for use in Crystallography and
!!----         Solid State Physics and Chemistry.
!!----
!!
 Module CFML_Maths
    !---- Use Modules ----!
    Use CFML_GlobalDeps

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public :: Co_Linear, Co_Prime, Cross_Product,                                & 
              Debye, Determ, Determ_V,       &
              Equal_Matrix, Equal_Vector,    &
              Factorial_I, Factorial_R, First_Derivative,                        &
              Gcd, Get_EPS_Math, Get_Cart_from_Cylin, Get_Cart_from_Spher,       &
              Get_Cylin_from_Cart, Get_Cylin_from_Spher, Get_Spher_from_Cart,    &
              Get_Spher_from_Cylin,          &
              Inverse_Array, In_Limits, Is_Diagonal_Matrix, Is_Null_Vector,      &
              Lcm, Linear_Dependent, Linear_Interpol, Locate, Lower_Triangular,  &
              Mat_Cross, Modulo_Lat,         &
              Negligible, Norm,              &
              Outerprod,                     &
              Poly_Legendre,                 &
              Rank, Rotation_OX, Rotation_OY, Rotation_OZ,                       &
              Scalar, Second_Derivative, Smoothing_Vec, Sort, Spline_Interpol, Spline_D2y,  &
              Tensor_Product, Trace,         &
              Upper_Triangular,              &
              Vec_Length,                    &
              Zbelong


    !---- List of public subroutines ----!
    public :: Co_Prime_Vector,               &
              Diagonalize_SH,                &
              Invert_Matrix_R,               &
              Points_In_Line2D,              &
              RowEchelonForm,                &
              Set_EPS_Math, SmithNormalForm, Sort_Strings, Svdcmp, Swap, &
              Resolv_Sist_1x2, Resolv_Sist_1x3, Resolv_Sist_2x2, Resolv_Sist_2x3,&
              Resolv_Sist_3x3
    


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

    real(kind=cp), parameter :: EPS_ARR=1.0E-12_cp  ! Internal epsilon value used for 
                                                    ! comparison in matrix operations
    real(kind=cp), parameter :: EPS_RTI= 1.0E-5_cp  ! Internal default epsilon for comparison
                                                    ! reals to integers

    !---- Variables ----!
    real(kind=cp) :: epss= EPS_RTI   

    !--------------------! 
    !---- Overloaded ----!
    !--------------------!
    Interface  Co_Linear
       Module Procedure Co_linear_C
       Module Procedure Co_linear_I
       Module Procedure Co_linear_R
    End Interface
    
    Interface  Cross_Product
       Module Procedure Cross_Product_C
       Module Procedure Cross_Product_I
       Module Procedure Cross_Product_R
    End Interface
    
    Interface  Determ
       Module Procedure Determinant_C
       Module Procedure Determinant_I
       Module Procedure Determinant_R
    End Interface
    
    Interface  Determ_V
       Module Procedure Determ_V_I
       Module Procedure Determ_V_R
    End Interface
    
    Interface  Diagonalize_SH
       Module Procedure Diagonalize_HERM
       Module Procedure Diagonalize_SYMM
    End Interface
    
    Interface  Equal_Matrix
       Module Procedure Equal_Matrix_C
       Module Procedure Equal_Matrix_I
       Module Procedure Equal_Matrix_R
    End Interface
    
    Interface  Equal_Vector
       Module Procedure Equal_Vector_C
       Module Procedure Equal_Vector_I
       Module Procedure Equal_Vector_R
    End Interface
    
    Interface Inverse_Array
       Module Procedure Inverse_Array_C
       Module Procedure Inverse_Array_I
       Module Procedure Inverse_Array_R
    End Interface
    
    Interface In_Limits
       Module Procedure In_Limits_I
       Module Procedure In_Limits_R
    End Interface
    
    interface Is_Diagonal_Matrix
        module procedure Is_Diagonal_Matrix_I
        module procedure Is_Diagonal_Matrix_R
    end interface
    
    interface Is_Null_Vector
        module procedure Is_Null_Vector_I
        module procedure Is_Null_Vector_R
    end interface
    
    Interface  Linear_Dependent
       Module Procedure Linear_Dependent_C
       Module Procedure Linear_Dependent_I
       Module Procedure Linear_Dependent_R
    End Interface
    
    Interface  Locate
       Module Procedure Locate_I
       Module Procedure Locate_R
    End Interface
    
    Interface  Lower_Triangular
       Module Procedure Lower_Triangular_I
       Module Procedure Lower_Triangular_R
    End Interface
    
    Interface  Mat_Cross
       Module Procedure Mat_Cross_C
       Module Procedure Mat_Cross_I
       Module Procedure Mat_Cross_R
    End Interface
    
    Interface  Negligible
       Module Procedure Negligible_C
       Module Procedure Negligible_R
    End Interface
    
    Interface Norm
       Module Procedure Norm_I
       Module Procedure Norm_R
    End Interface Norm
    
    Interface RowEchelonForm               
       Module Procedure RowEchelonFormM
       Module Procedure RowEchelonFormT
    End Interface RowEchelonForm
    
    Interface Scalar
       Module Procedure Scalar_I
       Module Procedure Scalar_R
    End Interface Scalar
    
    Interface Sort
       Module Procedure Sort_I
       Module Procedure Sort_R
    End Interface Sort
    
    Interface Swap
        Module Procedure Swap_C
        Module Procedure Swap_I
        Module Procedure Swap_R
        Module Procedure Swap_masked_C
        Module Procedure Swap_masked_I
        Module Procedure Swap_masked_R
    End interface
    
    Interface  Tensor_Product
       Module Procedure Tensor_product_C
       Module Procedure Tensor_product_I
       Module Procedure Tensor_product_R
    End Interface
    
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
       Module Procedure Zbelong_M
       Module Procedure Zbelong_R
       Module Procedure Zbelong_V
    End Interface
    

    !------------------------! 
    !---- Interface Zone ----!
    !------------------------!
    Interface
       
       Module Function Co_linear_C(a,b,n) Result(co_linear)
          !---- Argument ----!
          complex(kind=cp), dimension(:), intent(in) :: a,b    ! Complex vectors
          integer,              optional, intent(in) :: n      ! Dimension of the vectors
          logical                                    :: co_linear
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
       
       Module Pure Function Cross_Product_C(u,v) Result(w)    
          !---- Argument ----!
          complex(kind=cp), dimension(3), intent( in) :: u    ! Vector 1
          complex(kind=cp), dimension(3), intent( in) :: v    ! Vector 2
          complex(kind=cp), dimension(3)              :: w    ! u x v
       End Function Cross_Product_C
       
       Module Pure Function Cross_Product_I(u,v) Result(w)    
          !---- Argument ----!
          integer, dimension(3), intent( in) :: u    ! Vector 1
          integer, dimension(3), intent( in) :: v    ! Vector 2
          integer, dimension(3)              :: w    ! u x v
       End Function Cross_Product_I  
       
       Module Pure Function Cross_Product_R(u,v) Result(w)    
          !---- Argument ----!
          real(kind=cp), dimension(3), intent( in) :: u    ! Vector 1
          real(kind=cp), dimension(3), intent( in) :: v    ! Vector 2
          real(kind=cp), dimension(3)              :: w    ! u x v
       End Function Cross_Product_R
       
       Module Function Debye(N,X) Result(Fval)
          !---- Arguments ----!
          integer,       intent(in) :: N ! Order of the Debye function
          real(kind=CP), intent(in) :: X ! Value
          real(kind=CP)             :: fval
       End Function Debye
       
       Module Function Debye_DP(N,X) Result(Fval)
          !---- Arguments ----!
          integer,       intent(in) :: N ! Order of the Debye function
          real(kind=DP), intent(in) :: X ! Value
          real(kind=DP)             :: fval
       End Function Debye_DP

       Module Function Debye1(X) Result(Fval)
          !---- Arguments ----!
          real(kind=DP), intent(in) :: X
          real(kind=DP)             :: fval
       End Function Debye1

       Module Function Debye2(X) Result(Fval)
          !---- Argument ----!
          real(kind=DP), intent(in) :: X
          real(kind=DP)             :: fval
       End Function Debye2

       Module Function Debye3(X) Result(Fval)
          !---- Argument ----!
          real(kind=DP), intent(in) :: X
          real(kind=DP)             :: fval
       End Function Debye3

       Module Function Debye4(X) Result(FVal)
          !---- Argument ----!
          real(kind=DP), intent(in) :: X
          real(kind=DP)             :: fval
       End Function Debye4

       Module Function DebyeN(n,x) Result(Fval)
          !---- Arguments ----!
          integer,       intent(in) :: N ! Order of Debye function
          real(kind=DP), intent(in) :: X
          real(kind=DP)             :: Fval
       End Function DebyeN
       
       Module Function Debye_PR_ChebyshevSeries(n, a, t) Result(fval)
          !---- Arguments ----!
          integer,                       intent(in) :: N    ! The no. of terms in the sequence
          real(kind=DP), dimension(0:N), intent(in) :: A    ! The coefficients of the Chebyshev series
          real(kind=DP),                 intent(in) :: T    ! The value at which the series is to be evaluated
          real(kind=DP)                             :: fval ! Return value
       End Function Debye_PR_ChebyshevSeries
       
       Module Pure Function Determ_V_I(Vec1,Vec2,Vec3) Result(det)    
          !---- Arguments ----!
          integer, dimension(3), intent(in) :: Vec1,Vec2,Vec3
          integer                           :: det
       End Function Determ_V_I   
       
       Module Pure Function Determ_V_R(Vec1,Vec2,Vec3) Result(det)    
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: Vec1,Vec2,Vec3
          real(kind=cp)                           :: det
       End Function Determ_V_R
       
       Module Pure Subroutine Diagonalize_EigenvSort(d,v,n,io)    
          !---- Arguments ----!
          real(kind=cp), dimension(:),   intent(in out) :: d
          real(kind=cp), dimension(:,:), intent(in out) :: v
          integer,                       intent(in)     :: n
          integer,                       intent(in)     :: io
       End Subroutine Diagonalize_EigenvSort
       
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

       Module Pure Subroutine Diagonalize_PR_Tred1(a,n,d,e)
          !---- Arguments ----!
          real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)
          integer,                       intent(in)     :: n
          real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)
       End Subroutine Diagonalize_PR_Tred1

       Module Pure Subroutine Diagonalize_PR_Tred2(a,n,d,e)
          !---- Arguments ----!
          real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)
          integer,                       intent(in)     :: n
          real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)
       End Subroutine Diagonalize_PR_Tred2
       
       Module Subroutine Diagonalize_Herm(a,n,e_val,e_vect,norder)
          !---- Arguments ----!
          complex(kind=cp),           dimension(:,:), intent( in)  :: A
          integer,                                    intent( in)  :: n
          real(kind=cp),              dimension(:),   intent(out)  :: E_val   ! Eigenvalues
          complex(kind=cp), optional, dimension(:,:), intent(out)  :: E_vect  ! Eigenvectors
          logical, optional,                          intent(in)   :: norder  ! If present no ordering
       End Subroutine Diagonalize_Herm

       Module Subroutine Diagonalize_Symm(A,n,E_Val,E_vect,norder)
          !---- Arguments ----!
          real(kind=cp),           dimension(:,:), intent( in)  :: A
          integer,                                 intent( in)  :: n
          real(kind=cp),           dimension(:),   intent(out)  :: E_val    ! Eigenvalues
          real(kind=cp), optional, dimension(:,:), intent(out)  :: E_vect   ! Eigenvectors
          logical,       optional,                 intent(in)   :: norder   ! If present no ordering
       End Subroutine Diagonalize_Symm
       
       Module Pure Function Co_Prime(V,Imax) result(Cop)
          !---- Arguments ----!
          integer, dimension(:),           intent(in) :: V          ! Input vector of numbers
          integer,               optional, intent(in) :: Imax       ! Maximun prime number to be tested
          Logical                                     :: Cop
       End Function Co_Prime
       
       Module Subroutine Co_Prime_Vector(V,Cop,Ifact)
          !---- Arguments ----!
          integer, dimension(:),           intent(in)  :: V              ! input integer vector
          integer, dimension(:),           intent(out) :: Cop            ! Output co-prime vector
          integer,               optional, intent(out) :: Ifact          ! Common multiplicative factor
       End Subroutine Co_Prime_vector
       
       Module Function Determinant_C(A,n) result(det)
          !---- Arguments ----!
          complex(kind=cp), dimension(:,:), intent( in) :: A      ! Input array NxN
          integer,                          intent( in) :: n      ! Dimension of A
          complex(kind=cp)                              :: Det    ! Value
       End Function Determinant_C
       
       Module Function Determinant_I(A,n) result(det)
          !---- Arguments ----!
          integer, dimension(:,:), intent( in) :: A      ! Input array NxN
          integer,                 intent( in) :: n      ! Dimension of A
          integer                              :: Det    ! Value
       End Function Determinant_I
       
       Module Function Determinant_R(A,n) result(det)
          !---- Arguments ----!
          real(kind=cp), dimension(:,:), intent( in) :: A      ! Input array NxN
          integer,                       intent( in) :: n      ! Dimension of A
          real(kind=cp)                              :: Det    ! Value
       End Function Determinant_R
       
       Module Pure Function Deter2_C(A) Result(Det)
          !---- arguments ----!
          complex(kind=cp), dimension(2,2), intent(in) :: A   !! Matrix
          complex(kind=cp)                             :: Det      !! Determinant
       End Function Deter2_C
       
       Module Pure Function Deter2_I(A) Result(Det)
          !---- arguments ----!
          integer, dimension(2,2), intent(in) :: A   !! Matrix
          integer                             :: Det      !! Determinant
       End Function Deter2_I
       
       Module Pure Function Deter2_R(A) Result(Det)
          !---- arguments ----!
          real(kind=cp), intent(in) :: A(2,2)   !! Matrix
          real(kind=cp)             :: Det      !! Determinant
       End Function Deter2_R
       
       Module Pure Function Deter3_C(A) Result(Det)
          !---- arguments ----!
          complex(kind=cp), dimension(3,3), intent(in) :: A   !! Matrix
          complex(kind=cp)                             :: Det      !! Determinant
       End Function Deter3_C
       
       Module Pure Function Deter3_I(A) Result(Det)
          !---- arguments ----!
          integer, dimension(3,3), intent(in) :: A   !! Matrix
          integer                             :: Det      !! Determinant
       End Function Deter3_I
       
       Module Pure Function Deter3_R(A) Result(Det)
          !---- arguments ----!
          real(kind=cp), dimension(3,3), intent(in) :: A   !! Matrix
          real(kind=cp)                             :: Det      !! Determinant
       End Function Deter3_R
       
       Module Pure Function Deter4_C(A) Result(Det)
          !---- arguments ----!
          complex(kind=cp), dimension(4,4), intent(in) :: A   !! Matrix
          complex(kind=cp)                             :: Det      !! Determinant
       End Function Deter4_C
       
       Module Pure Function Deter4_I(A) Result(Det)
          !---- arguments ----!
          integer, dimension(4,4), intent(in) :: A   !! Matrix
          integer                             :: Det      !! Determinant
       End Function Deter4_I
       
       Module Pure Function Deter4_R(A) Result(Det)
          !---- arguments ----!
          real(kind=cp), dimension(4,4), intent(in) :: A   !! Matrix
          real(kind=cp)                             :: Det      !! Determinant
       End Function Deter4_R
       
       Module Function DeterN_C(A,n) result(det)
          !---- Arguments ----!
          complex(kind=cp), dimension(:,:), intent( in) :: A      ! Input array NxN
          integer,                          intent( in) :: n      ! Dimension of A
          complex(kind=cp)                              :: Det    ! Value
       End Function DeterN_C
       
       Module Function DeterN_I(A,n) result(det)
          !---- Arguments ----!
          integer, dimension(:,:), intent( in) :: A      ! Input array NxN
          integer,                 intent( in) :: n      ! Dimension of A
          integer                              :: Det    ! Value
       End Function DeterN_I
       
       Module Function DeterN_R(A,n) result(det)
          !---- Arguments ----!
          real(kind=cp), dimension(:,:), intent( in) :: A      ! Input array NxN
          integer,                       intent( in) :: n      ! Dimension of A
          real(kind=cp)                              :: Det    ! Value
       End Function DeterN_R
       
       Module Function Equal_Matrix_C(a,b,n) result(info)
          !---- Argument ----!
          complex(kind=cp), dimension(:,:),       intent(in) :: a,b      ! Input arrays NxN
          integer,                      optional, intent(in) :: n        ! Dimensions N
          logical                                            :: info
       End Function Equal_Matrix_C
       
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
       
       Module Function Equal_Vector_C(a,b,n) result(info)
          !---- Argument ----!
          complex(kind=cp), dimension(:),           intent(in) :: a,b      ! Input vectors
          integer,                        optional, intent(in) :: n        ! Dimension of the vector
          logical                                           :: info
       End Function Equal_Vector_C
       
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
       
       Module Elemental Function Factorial_I(N) Result(Fact)
          !---- Argument ----!
          integer, intent(in) :: N        ! Factorial of N
          integer             :: Fact
       End Function Factorial_I

       Module Elemental Function Factorial_R(N) Result(Fact)
          !---- Arguments ----!
          integer, intent(in) :: N    ! Factorial of N
          real(kind=cp)       :: Fact
       End Function Factorial_R
       
       Module Pure Function First_Derivative(x,y,n) Result(d1y)
          !---- Arguments ----!
          real(kind=cp), dimension(:), intent(in)  :: x    ! Vector containing Xi
          real(kind=cp), dimension(:), intent(in)  :: y    ! Vector containing Yi
          integer ,                    intent(in)  :: n    ! Dimension
          real(kind=cp), dimension(n)              :: d1y  ! Vector containing the first derivative
       End Function First_Derivative
       
       Module Elemental Function Gcd(a,b) Result(mcd)
          !---- Arguments ----!
          integer, intent(in) :: a,b   ! Integers
          integer             :: mcd   ! Maximum common divisor
       End Function Gcd
       
       Module Pure Function Get_Cart_from_Cylin(CilCoord,Mode) Result(CarCoord)   
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent( in) ::  CilCoord ! Coordinates rho,phi,zeta
          character(len=*), optional,  intent( in) ::  mode     ! "D" angles in degrees, otherwise in radians
          real(kind=cp), dimension(3)              ::  CarCoord ! Cartesian coordinates
       End Function Get_Cart_from_Cylin
       
       Module Pure Function Get_Cart_from_Spher(SphCoord,Mode) Result(CarCoord)   
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent( in) :: SphCoord ! Coordinates (R,Theta;Phi)
          character(len=*), optional,  intent( in) :: mode     ! If "D" the angles are in degrees, otherwise radians is considered
          real(kind=cp), dimension(3)              :: CarCoord ! Cartesian coordinates
       End Function Get_Cart_from_Spher
       
       Module Pure Function Get_Cylin_from_Cart(CarCoord, Mode) Result(CilCoord)    
          !---- Arguments ----!
          real(kind=cp), dimension(3),intent(in) ::  CarCoord   ! Cartesian coordinatates
          character(len=*), optional, intent(in) ::  mode
          real(kind=cp), dimension(3)            ::  CilCoord   ! Cylindrical coordinates
       End Function Get_Cylin_from_Cart
       
       Module Pure Function Get_Cylin_from_Spher(SphCoord,mode) Result(CilCoord)    
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: SphCoord ! Cylinder
          character(len=*), optional,  intent(in) :: mode
          real(kind=cp), dimension(3)             :: CilCoord ! Spherical
       End Function Get_Cylin_from_Spher
       
       Module Pure Function Get_Spher_from_Cart(CarCoord,mode) Result(SphCoord)    
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: CarCoord ! Cartesian
          character(len=*), optional,  intent(in) :: mode
          real(kind=cp), dimension(3)             :: SphCoord ! Spherical
       End Function Get_Spher_from_Cart
       
       Module Pure Function Get_Spher_from_Cylin(CilCoord,mode) Result(SphCoord)    
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: CilCoord ! Cylinder
          character(len=*), optional,  intent(in) :: mode
          real(kind=cp), dimension(3)             :: SphCoord ! Spherical
       End Function Get_Spher_from_Cylin

       Module Function Inverse_Array_C(A) Result(Ainv)
          !---- Arguments ----!
          complex(kind=cp), dimension(:,:), intent(in)     :: A
          complex(kind=cp), dimension(size(a,1),size(a,1)) :: Ainv 
       End Function Inverse_Array_C   
       
       Module Function Inverse_Array_I(A) Result(Ainv)
          !---- Arguments ----!
          integer, dimension(:,:),       intent(in)     :: A
          real(kind=cp), dimension(size(a,1),size(a,1)) :: Ainv 
       End Function Inverse_Array_I
       
       Module Function Inverse_Array_R(A) Result(Ainv)
          !---- Arguments ----!
          real(kind=cp), dimension(:,:), intent(in)     :: A
          real(kind=cp), dimension(size(a,1),size(a,1)) :: Ainv 
       End Function Inverse_Array_R
       
       Module Pure Function In_Limits_I(v,limits,n) result(ok)    
          !---- Arguments ----!
          integer, dimension(:),             intent(in) :: v        ! Input Vector
          integer, dimension(:,:),           intent(in) :: limits   ! Normally (2,n)
          integer,                 optional, intent(in) :: n        ! Dimension of vect
          logical                                       :: ok
       End Function In_Limits_I 
       
       Module Pure Function In_Limits_R(v,limits,n) result(ok)    
          !---- Arguments ----!
          real(kind=cp), dimension(:),             intent(in) :: v        ! Input Vector
          real(kind=cp), dimension(:,:),           intent(in) :: limits   ! Normally (2,n)
          integer,                       optional, intent(in) :: n        ! Dimension of vect
          logical                                             :: ok
       End Function In_Limits_R
       
       Module Subroutine Invert_Matrix_R(a,b,perm)    
          !---- Arguments ----!
          real(kind=cp), dimension(:,:),  intent(in ) :: a         ! Input Array
          real(kind=cp), dimension(:,:),  intent(out) :: b         ! Output array
          integer, dimension(:),optional, intent(out) :: perm
       End Subroutine Invert_Matrix_R
       
       Module Pure Function Is_Diagonal_Matrix_I(A) Result(info)
          !---- Arguments ----! 
          integer, dimension(:,:), intent(in)  :: A
          logical                              :: info
       End Function Is_Diagonal_Matrix_I
       
       Module Pure Function Is_Diagonal_Matrix_R(A) Result(info)
          !---- Arguments ----! 
          real(kind=cp), dimension(:,:), intent(in)  :: A
          logical                                    :: info
       End Function Is_Diagonal_Matrix_R
       
       Module Pure Function Is_Null_Vector_I(V) Result(info)
          !---- Arguments ----!
          integer,  dimension(:), intent(in)  :: V
          logical                             :: Info
       End Function Is_Null_Vector_I   
       
       Module Pure Function Is_Null_Vector_R(V) Result(info)
          !---- Arguments ----!
          real(kind=cp),  dimension(:), intent(in)  :: V
          logical                             :: Info
       End Function Is_Null_Vector_R
       
       Module Elemental Function Lcm(a,b) result(mcm)
          !---- Arguments ----!
          integer, intent(in) :: a,b    ! Integers
          integer             :: mcm    ! Minimum common multiple
       End Function Lcm
       
       Module Function Linear_Dependent_C(A,na,B,nb,mb) Result(info)
          !---- Arguments ----!
          complex(kind=cp), dimension(:),   intent(in)  :: a
          complex(kind=cp), dimension(:,:), intent(in)  :: b
          integer,                          intent(in)  :: na,nb,mb
          logical                                       :: info
       End Function Linear_Dependent_C

       Module Function Linear_Dependent_I(A,na,B,nb,mb) Result(info)
          !---- Arguments ----!
          integer, dimension(:),   intent(in)  :: a
          integer, dimension(:,:), intent(in)  :: b
          integer,                 intent(in)  :: na,nb,mb
          logical                              :: info
       End Function Linear_Dependent_I

       Module Function Linear_Dependent_R(A,na,B,nb,mb) Result(info)
          !---- Arguments ----!
          real(kind=cp), dimension(:),   intent(in)  :: a
          real(kind=cp), dimension(:,:), intent(in)  :: b
          integer,                       intent(in)  :: na,nb,mb
          logical                                    :: info
       End Function Linear_Dependent_R
       
       Module Pure Function Linear_Interpol(x,y,xi) Result(yi)
          !---- Arguments ----!
          real(kind=cp), dimension(:),intent(in)   :: x  ! Vector containing Xi points
          real(kind=cp), dimension(:),intent(in)   :: y  ! Vector Yi=F(xi)
          real(kind=cp),              intent(in)   :: xi ! X point to evaluate
          real(kind=cp)                            :: yi ! Output
       End Function Linear_Interpol   
       
       Module Pure Function Locate_I(V,x,n) Result(j)
          !---- Argument ----!
          integer, dimension(:), intent(in):: v  ! Input vector
          integer,               intent(in):: x  ! Value
          integer, optional,     intent(in):: n  ! Value
          integer                          :: j
       End Function Locate_I

       Module Pure Function Locate_R(V,x,n) Result(j)
          !---- Argument ----!
          real(kind=cp), dimension(:), intent(in):: v
          real(kind=cp),               intent(in):: x
          integer, optional,           intent(in):: n  ! Value
          integer                                :: j
       End Function Locate_R
       
       Module Pure Function Lower_Triangular_I(A,n) Result (T)
          !---- Argument ----!
          integer, dimension(:,:), intent(in) :: A    ! Input array
          integer,                 intent(in) :: n    ! Dimension of array
          integer, dimension(n,n)             :: T
       End Function Lower_Triangular_I

       Module Pure Function Lower_Triangular_R(A,n) Result (T)
          !---- Argument ----!
          real(kind=cp), dimension(:,:), intent(in) :: A    ! Input Array
          integer,                       intent(in) :: n    ! Dimension of A
          real(kind=cp), dimension(n,n)             :: T
       End Function Lower_Triangular_R
       
       Module Subroutine LU_Backsub(a,indx,b)
          !---- Arguments ----!
          real(kind=cp), dimension(:,:), intent(in)     :: a
          integer,         dimension(:), intent(in)     :: indx
          real(kind=cp),   dimension(:), intent(in out) :: b
       End Subroutine LU_Backsub

       Module Subroutine LU_Decomp(a,d,singular,indx)
          !---- Arguments ----!
          real(kind=cp), dimension(:,:),   intent(in out) :: a
          real(kind=cp),                   intent(out)    :: d
          logical,                         intent(out)    :: singular
          integer, dimension(:), optional, intent(out)    :: indx
       End Subroutine LU_Decomp
       
       Module Pure Function Mat_Cross_C(Vec) Result(M)    
          !---- Argument ----!
          complex(kind=cp), dimension(3), intent( in) :: Vec
          complex(kind=cp), dimension(3,3)            :: M
       End Function Mat_Cross_C  
       
       Module Pure Function Mat_Cross_I(Vec) Result(M)    
          !---- Argument ----!
          integer, dimension(3), intent( in) :: Vec
          integer, dimension(3,3)            :: M
       End Function Mat_Cross_I 
       
       Module Pure Function Mat_Cross_R(Vec) Result(M)    
          !---- Argument ----!
          real(kind=cp), dimension(3), intent( in) :: Vec
          real(kind=cp), dimension(3,3)            :: M
       End Function Mat_Cross_R 
       
       Module Pure Function MatInv2_C(A) Result(B)
          !---- arguments ----!
          complex(kind=cp), dimension(2,2), intent(in) :: A   
          complex(kind=cp), dimension(2,2)             :: B    
       End Function MatInv2_C
       
       Module Pure Function MatInv2_R(A) Result(B)
          !---- arguments ----!
          real(kind=cp), dimension(2,2), intent(in) :: A   
          real(kind=cp), dimension(2,2)             :: B   
       End Function MatInv2_R
       
       Module Pure Function MatInv3_C(A) Result(B)
          !---- arguments ----!
          complex(kind=cp), dimension(3,3), intent(in) :: A   
          complex(kind=cp), dimension(3,3)             :: B    
       End Function MatInv3_C
       
       Module Pure Function MatInv3_R(A) Result(B)
          !---- arguments ----!
          real(kind=cp), dimension(3,3), intent(in) :: A   
          real(kind=cp), dimension(3,3)             :: B   
       End Function MatInv3_R
       
       Module Pure Function MatInv4_C(A) Result(B)
          !---- arguments ----!
          complex(kind=cp), dimension(4,4), intent(in) :: A   
          complex(kind=cp), dimension(4,4)             :: B    
       End Function MatInv4_C
       
       Module Pure Function MatInv4_R(A) Result(B)
          !---- arguments ----!
          real(kind=cp), dimension(4,4), intent(in) :: A   
          real(kind=cp), dimension(4,4)             :: B   
       End Function MatInv4_R
       
       Module Function MatInvN_C(A,n) Result(Ainv)    
          !---- Arguments ----!
          complex(kind=cp), dimension(:,:), intent(in) :: a
          integer,                       intent(in) :: n
          complex(kind=cp), dimension(n,n)             :: Ainv
       End Function MatInvN_C 
       
       Module Function MatInvN_R(A,n) Result(Ainv)    
          !---- Arguments ----!
          real(kind=cp), dimension(:,:), intent(in) :: a
          integer,                       intent(in) :: n
          real(kind=cp), dimension(n,n)             :: Ainv
       End Function MatInvN_R   
       
       Module Pure Function Modulo_Lat(v) result(u)
          !---- Argument ----!
          real(kind=cp), dimension(:), intent( in) :: v
          real(kind=cp), dimension(1:size(v))      :: u
       End Function Modulo_Lat
       
       Module Elemental Function Negligible_C(C) Result(Neglig)
          !---- Argument ----!
          complex(kind=cp), intent( in) :: C         ! Complex number
          logical                       :: Neglig
       End Function Negligible_C

       Module Elemental Function Negligible_R(R) Result(neglig)
          !---- Argument ----!
          real(kind=cp), intent( in) :: R          ! Real number
          logical                    :: Neglig
       End Function Negligible_R
       
       Module Pure Function Norm_I(X,G) Result(R)
          !---- Arguments ----!
          integer,       dimension(:),   intent(in) :: x    ! Input vector
          real(kind=cp), dimension(:,:), intent(in) :: g    ! Metric array
          real(kind=cp)                             :: r    ! Norm of the input vector
       End Function Norm_I

       Module Pure Function Norm_R(X,G) Result(R)
          !---- Arguments ----!
          real(kind=cp), dimension(:),   intent(in) :: x   ! Input vector
          real(kind=cp), dimension(:,:), intent(in) :: g   ! Metrics
          real(kind=cp)                             :: r   ! Norm of the vector
       End Function Norm_R
       
       Module Pure Function Outerprod(a,b)  Result(c)
          !---- Arguments ----!
          real(kind=cp),dimension(:),intent(in)    :: a,b
          real(kind=cp),dimension(size(a),size(b)) :: c
       End Function Outerprod
       
       Module Pure Subroutine Points_In_Line2D(X1, XN, N, XP)    
          !---- Arguments ----!
          real(kind=cp), dimension(2),   intent(in)  :: X1   ! Point1 in 2D
          real(kind=cp), dimension(2),   intent(in)  :: XN   ! PointN in 2D
          integer,                       intent(in)  :: N    ! Number of Total points
          real(kind=cp), dimension(:,:), intent(out) :: XP   ! List of points
       End Subroutine Points_In_Line2D  
       
       Module Elemental Function Poly_Legendre(L,M,X) Result(Plmx)
          !---- Arguments ----!
          integer,      intent (in) :: L
          integer,      intent (in) :: M
          real(kind=cp),intent (in) :: X
          real(kind=cp)             :: Plmx
       End Function Poly_Legendre
       
       Module Function PseudoDeterm_C(A,n) result(det)
          !---- Arguments ----!
          complex(kind=cp), dimension(:,:), intent( in) :: A      ! Input array NxN
          integer,                          intent( in) :: n      ! Dimension of A
          real(kind=cp)                              :: Det    ! Value
       End Function PseudoDeterm_C
       
       Module Function Rank(a,tol) Result(r)    
          !---- Arguments ----!
          real(kind=cp), dimension(:,:),intent( in)      :: a     ! Input array
          real(kind=cp),                intent( in)      :: tol   ! Tolerance
          integer                                        :: r
       End Function Rank   
       
       Module Subroutine Resolv_Sist_1x2(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(2),         intent( in) :: w      ! Input vector       
          real(kind=cp),                 intent( in) :: t      ! Input value       
          real(kind=cp), dimension(2),   intent(out) :: ts     ! Fixed value solution       
          real(kind=cp), dimension(2),   intent(out) :: x      ! Fixed value for x,y       
          integer, dimension(2),         intent(out) :: ix     ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_1x2    
 
       Module Subroutine Resolv_Sist_1x3(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(3),         intent( in) :: w      ! Input vector       
          real(kind=cp),                 intent( in) :: t      ! Input value       
          real(kind=cp), dimension(3),   intent(out) :: ts     ! Fixed value solution       
          real(kind=cp), dimension(3),   intent(out) :: x      ! Fixed value for x,y,z       
          integer, dimension(3),         intent(out) :: ix     ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_1x3    
 
       Module Subroutine Resolv_Sist_2x2(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(2,2),       intent( in) :: w       ! Input vector       
          real(kind=cp),dimension(2),    intent( in) :: t       ! Input value       
          real(kind=cp),dimension(2),    intent(out) :: ts      ! Fixed value solution       
          real(kind=cp),dimension(2),    intent(out) :: x       ! Fixed value for x,y       
          integer, dimension(2),         intent(out) :: ix      ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_2x2    
 
       Module Subroutine Resolv_Sist_2x3(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(2,3),          intent( in) :: w     ! Input vector       
          real(kind=cp), dimension(2),      intent( in) :: t     ! Input value       
          real(kind=cp), dimension(3),      intent(out) :: ts    ! Fixed value solution       
          real(kind=cp), dimension(3),      intent(out) :: x     ! Fixed value for x,y       
          integer, dimension(3),            intent(out) :: ix    ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_2x3    
 
       Module Subroutine Resolv_Sist_3x3(w,t,ts,x,ix)    
          !---- Arguments ----!       
          integer, dimension(3,3),          intent(in) :: w       ! Input vector       
          real(kind=cp), dimension(3),      intent(in) :: t       ! Input value       
          real(kind=cp), dimension(3),      intent(out):: ts      ! Fixed value solution       
          real(kind=cp), dimension(3),      intent(out):: x       ! Fixed value for x,y       
          integer, dimension(3),            intent(out):: ix      ! 1: x, 2: y, 3: z       
       End Subroutine Resolv_Sist_3x3 
       
       Module Pure Function Rotation_OX(Vec,Angle) Result(Rvec)    
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: Vec      ! Vector
          real(kind=cp),               intent(in) :: angle    ! Angle
          real(kind=cp), dimension(3)             :: Rvec     ! Vector rotated
       End Function Rotation_OX  
       
       Module Pure Function Rotation_OY(Vec,Angle) Result(Rvec)    
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: Vec      ! Vector
          real(kind=cp),               intent(in) :: angle    ! Angle
          real(kind=cp), dimension(3)             :: Rvec     ! Vector rotated
       End Function Rotation_OY  
       
       Module Pure Function Rotation_OZ(Vec,Angle) Result(Rvec)    
          !---- Arguments ----!
          real(kind=cp), dimension(3), intent(in) :: Vec      ! Vector
          real(kind=cp),               intent(in) :: angle    ! Angle
          real(kind=cp), dimension(3)             :: Rvec     ! Vector rotated
       End Function Rotation_OZ 
       
       Module Pure Subroutine RowEchelonFormM(M)
          !---- Arguments ----!
          integer, dimension(:,:), intent(in out) :: M
       End Subroutine RowEchelonFormM
       
       Module Pure Subroutine RowEchelonFormT(M,T)
          !---- Arguments ----!
          integer, dimension(:,:), intent(in out) :: M
          integer, dimension(:,:), intent(in out) :: T 
       End Subroutine RowEchelonFormT     
       
       Module Pure Function Scalar_I(X,Y,G) Result(R)
          !---- Arguments ----!
          integer, dimension(:),         intent(in) :: x,y     ! Input vectors
          real(kind=cp), dimension(:,:), intent(in) :: g       ! Metrics
          real(kind=cp)                             :: r       ! Scalar
       End Function Scalar_I

       Module Pure Function Scalar_R(X,Y,G) Result(R)
          !---- Arguments ----!
          real(kind=cp), dimension(:),   intent(in) :: x,y    ! Input vectors
          real(kind=cp), dimension(:,:), intent(in) :: g      ! Metrics
          real(kind=cp)                             :: r      ! Scalar
       End Function Scalar_R
       
       Module Pure Function Second_Derivative(x,y,n) Result(d2y)
          !---- Arguments ----!
          real(kind=cp), dimension(:), intent(in)  :: x     ! Vector xi
          real(kind=cp), dimension(:), intent(in)  :: y     ! Vector Yi=F(xi) 
          integer ,                    intent(in)  :: n     ! Dimension
          real(kind=cp), dimension(n)              :: d2y   ! Second derivate
       End Function Second_Derivative   
       
       Module Pure Subroutine SmithNormalForm(M,D,P,Q)
          !---- Arguments ----!
          integer, dimension(:,:), intent(in)  :: M !(nr,nc)
          integer, dimension(:,:), intent(out) :: D !(nr,nc)
          integer, dimension(:,:), intent(out) :: P !(nr,nr)
          integer, dimension(:,:), intent(out) :: Q !(nc,nc)
       End Subroutine SmithNormalForm  
       
       Module Pure Function Smoothing_Vec(Y, N, Niter) Result(Ys)    
          !---- Arguments ----!
          real(kind=cp),dimension(:),            intent(in) :: Y         !  In Out-> Array to be smoothed
          integer,                               intent(in) :: n         !  In -> Number of points
          integer,                               intent(in) :: niter     !  In -> Number of iterations
          real(kind=cp),dimension(n)                        :: Ys        !  Out-> Array smoothed
       End Function Smoothing_Vec
       
       Module Function Sort_I(arr,n) Result(indx)    
          !---- Arguments ----!
          integer, dimension(:), intent(in ) :: arr       ! Vector
          integer              , intent(in ) :: n         ! Dimension
          integer, dimension(n)              :: indx      ! Index
       End Function Sort_I   
       
       Module Pure Subroutine Sort_PR_Partition(A, Marker)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in out) :: A
          integer,                        intent(   out) :: marker
       End Subroutine Sort_PR_Partition
       
       Module Function Sort_R(arr,n) Result(indx)    
          !---- Arguments ----!
          real(kind=cp), dimension(:), intent(in ) :: arr       ! Vector
          integer                    , intent(in ) :: n         ! Dimension
          integer, dimension(n)                    :: indx      ! Index
       End Function Sort_R 
       
       Module Recursive Subroutine Sort_Strings(Str)
          !---- Argument ----!
          character(len=*), dimension(:), intent(in out) :: Str
       End Subroutine Sort_Strings
       
       Module Pure Function Spline_D2Y(x,y,n,yp1,ypn) Result(ys)
          !---- Arguments ----!
          real(kind=cp), dimension(:), intent(in)  :: x               !  In -> Array X
          real(kind=cp), dimension(:), intent(in)  :: y               !  In -> Array Yi=F(Xi)
          integer ,                    intent(in)  :: n               !  In -> Dimension of X, Y
          real(kind=cp),               intent(in)  :: yp1             !  In -> Derivate of Point 1
          real(kind=cp),               intent(in)  :: ypn             !  In -> Derivate of Point N
          real(kind=cp), dimension(n)            :: ys              ! Out -> array containing second derivatives
       End Function Spline_D2Y
       
       Module Pure Function Spline_Interpol(x,y,d2y,n,xi) Result(yi)
          !---- Arguments ----!
          real(kind=cp), dimension(:), intent(in)  :: x    ! Vector Xi points
          real(kind=cp), dimension(:), intent(in)  :: y    ! Vector Yi points
          real(kind=cp), dimension(:), intent(in)  :: d2y  ! Vector Second derivate of Yi points
          integer ,                    intent(in)  :: n    ! Dimension of vectors
          real(kind=cp),               intent(in)  :: xi   ! X value for evaluation
          real(kind=cp)                            :: yi
       End Function Spline_Interpol
       
       Module Subroutine Svdcmp(a,w,v)    
          !---- Arguments ----!
          real(kind=cp),dimension(:,:),intent(in out) ::a   ! A(m,n)
          real(kind=cp),dimension(:),  intent(   out) ::w   ! W(n)
          real(kind=cp),dimension(:,:),intent(   out) ::v   ! V(n,n)
       End Subroutine Svdcmp
       
       Module Elemental Subroutine Swap_C(a,b)
          !---- Arguments ----!
          complex(kind=cp), intent(in out) :: a
          complex(kind=cp), intent(in out) :: b
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
       
       Module Elemental Subroutine Swap_Masked_C(A,B,Mask)
          !---- Arguments ----!
          complex(kind=cp), intent(in out) :: a
          complex(kind=cp), intent(in out) :: b
          logical,           intent(in) :: mask
       End Subroutine Swap_Masked_C
       
       Module Elemental Subroutine Swap_Masked_I(A,B,Mask)
          !---- Arguments ----!
          integer, intent(in out) :: a
          integer, intent(in out) :: b
          logical,           intent(in) :: mask
       End Subroutine Swap_Masked_I
       
       Module Elemental Subroutine Swap_Masked_R(A,B,Mask)
          !---- Arguments ----!
          real(kind=cp), intent(in out) :: a
          real(kind=cp), intent(in out) :: b
          logical,           intent(in) :: mask
       End Subroutine Swap_Masked_R
       
       Module Pure Function Tensor_Product_C(Vec1,Vec2) Result(w)    
          !---- Argument ----!
          complex(kind=cp), dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2
          complex(kind=cp), dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2
       End Function Tensor_Product_C
       
       Module Pure Function Tensor_Product_I(Vec1,Vec2) Result(w)    
          !---- Argument ----!
          integer, dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2
          integer, dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2
       End Function Tensor_Product_I
       
       Module Pure Function Tensor_Product_R(Vec1,Vec2) Result(w)    
          !---- Argument ----!
          real(kind=cp), dimension(3), intent( in) :: Vec1,Vec2  ! Vector 1, Vector 2
          real(kind=cp), dimension(3,3)            :: w          ! Tensor product Vector1 (o) Vector2
       End Function Tensor_Product_R
       
       Module Pure Function Trace_C(a) Result(b)
          !---- Argument ----!
          complex(kind=cp), dimension(:,:), intent(in) :: a
          complex(kind=cp)                             :: b
       End Function Trace_C

       Module Pure Function Trace_I(a) Result(b)
          !---- Argument ----!
          integer, dimension(:,:), intent(in) :: a
          integer                             :: b
       End Function Trace_I

       Module Pure Function Trace_R(a) Result(b)
          !---- Argument ----!
          real(kind=cp), dimension(:,:), intent(in) :: a
          real(kind=cp)                             :: b
       End Function Trace_R
       
       Module Pure Function Upper_Triangular_I(A,n) Result (T)
          !---- Argument ----!
          integer, dimension(:,:), intent(in) :: A     ! Input array
          integer,                 intent(in) :: n     ! Dimension
          integer, dimension(n,n)             :: T
       End Function Upper_Triangular_I

       Module Pure Function Upper_Triangular_R(A,n) Result (T)
          !---- Argument ----!
          real(kind=cp), dimension(:,:), intent(in) :: A   ! Input array
          integer,                       intent(in) :: n   ! Dimension
          real(kind=cp), dimension(n,n)             :: T
       End Function  Upper_Triangular_R
       
       Module Pure Function Vec_Length(G,Vec) Result(c)    
          !---- Arguments ----!
          real(kind=cp), intent(in)  , dimension(3,3)       :: G      ! Metric array
          real(kind=cp), intent(in)  , dimension(3  )       :: Vec    ! Vector
          real(kind=cp)                                     :: c      ! Length of Vector
       End Function Vec_Length
       
       Module Pure Function Zbelong_M(A) Result(belong)
          !---- Argument ----!
          real(kind=cp),   dimension(:,:), intent( in) :: A        ! Input array
          logical                                      :: belong
       End Function Zbelong_M

       Module Pure Function Zbelong_R(X) Result(belong)
          !---- Argument ----!
          real(kind=cp), intent( in) :: X              ! Input number
          logical                    :: belong
       End Function Zbelong_R

       Module Pure Function Zbelong_V(V) Result(belong)
          !---- Argument ----!
          real(kind=cp),   dimension(:), intent( in) :: v      ! Input vector
          logical                                    :: belong
       End Function Zbelong_V
       
    End Interface

 Contains

    !!----
    !!---- GET_EPS_MATH()
    !!----    Gets global EPSS
    !!----
    !!---- 27/03/2019 
    !!
    Function Get_Eps_Math() result(V)
       !---- Arguments ----!
       real(kind=cp) :: v

       v=epss
       
       return
    End Function Get_Eps_Math

    !!---- 
    !!---- SET_EPS_MATH
    !!----    Sets/Modify global EPSS.
    !!----    Calling without arguments set to default value
    !!----
    !!---- 27/03/2019
    !!
    Subroutine Set_Eps_Math(Neweps)
       !---- Arguments ----!
       real(kind=cp), optional, intent( in) :: neweps

       if (present(neweps)) then
          epss=neweps
       else
          epss=EPS_RTI
       end if

       return
    End Subroutine Set_Eps_Math

End Module CFML_Maths
