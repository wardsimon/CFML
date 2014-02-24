!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2013  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----                          Universita di Padova, Padova, ITALY
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL) (CrysFML)
!!----          Javier Gonzalez-Platas  (ULL) (CrysFML)
!!----          Ross John Angel               (EoS)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----               Ross John Angel    (Dpto. Geoscienze, Universita di Padova, Italy)
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
!!---- MODULE: CFML_EoS
!!----   INFO:
!!----
!!---- HISTORY
!!----    Update: 16/02/2013
!!----    Incorporates bug fix to get_volume: JGP 23/02/2013
!!----    Validation for non-thermal EOS and modifications: RJA 22/02/2013 .....
!!----
!!---- DEPENDENCIES
!!----    CFML_GLOBALDEPS,       only: cp
!!----    CFML_CRYSTAL_METRICS,  only: Crystal_Cell_Type,Get_Cryst_Family,Volume_Sigma_from_Cell
!!----    CFML_STRING_UTILITIES
!!----
!!---- PARAMETERS
!!----    N_EOSPAR
!!----    NCOL_DATA_MAX
!!----
!!---- VARIABLES
!!----    ERR_EOS
!!----    ERR_EOS_MESS
!!----    EOS_TYPE
!!----    EOS_LIST_TYPE
!!----    EOS_DATA_TYPE
!!----    EOS_DATA_LIST_TYPE
!!----
!!---- PROCEDURES
!!----
!!----    Functions:
!!----    ----------------------
!!----    ALPHA_CAL
!!----    DKDT_CAL
!!----    GET_K0_T                       [PRIVATE]
!!----    GET_PRESSURE
!!----    GET_PRESSURE_ESD
!!----    GET_TEMPERATURE
!!----    GET_V0_T                       [PRIVATE]
!!----    GET_VOLUME
!!----    GET_VOLUME_S
!!----    K_CAL
!!----    KP_CAL
!!----    KPP_CAL
!!----    NORMPRESSURE_EOS               [PRIVATE]
!!----    NORMPRESSURE_P                 [PRIVATE]
!!----    PRESSURE_F
!!----    PTHERMAL                       [PRIVATE]
!!----    STRAIN
!!----    STRAIN_EOS
!!----
!!----
!!----    Subroutines:
!!----    --------------------------
!!----    ALLOCATE_EOS_DATA_LIST
!!----    ALLOCATE_EOS_LIST
!!----    DEALLOCATE_EOS_DATA_LIST
!!----    DEALLOCATE_EOS_LIST
!!----    DEF_CRYSTAL_SYSTEM             [PRIVATE]
!!----    DERIV_PARTIAL_P
!!----    DERIV_PARTIAL_P_ANALYTIC       [PRIVATE]
!!----    DERIV_PARTIAL_P_NUMERIC        [PRIVATE]
!!----    EOS_CAL
!!----    EOS_CAL_ESD
!!----    EOS_TO_VEC                     [PRIVATE]
!!----    FFCAL_DAT
!!----    FFCAL_DAT_ESD
!!----    FFCAL_EOS
!!----    GET_TAIT                       [PRIVATE]
!!----    INIT_EOS_DATA_TYPE
!!----    INIT_EOS_THERMAL
!!----    INIT_EOS_TYPE
!!----    INIT_ERR_EOS
!!----    PHYSICAL_CHECK                 [PRIVATE]
!!----    READ_EOS_DATAFILE
!!----    SET_EOS_NAMES
!!----    SET_EOS_USE
!!----    SET_KP_KPP_COND
!!----    SET_VOLUME_FROM_CELL           [PRIVATE]
!!----    TRANSFORM_ESD                  [PRIVATE]
!!----    VEC_TO_EOS                     [PRIVATE]
!!----    WRITE_INFO_EOS
!!----    WRITE_INFO_EOS_THERMAL         [PRIVATE]
!!----
!!
Module CFML_EoS
   !---- Use Modules ----!
   Use CFML_GlobalDeps,       only: cp
   Use CFML_Crystal_Metrics,  only: Crystal_Cell_Type,Get_Cryst_Family,Volume_Sigma_from_Cell
   Use CFML_String_Utilities

   implicit none

   private

   !---- List of public functions ----!
   public :: Alpha_Cal, Dkdt_Cal,Get_Pressure, Get_Pressure_Esd, Get_Temperature, Get_Volume, Get_Volume_S, &
             K_Cal, Kp_Cal, Kpp_Cal, Pressure_F, Strain, Strain_EOS


   !---- List of public subroutines ----!
   public :: Allocate_EoS_Data_List, Allocate_EoS_List, Deallocate_EoS_Data_List, Deallocate_EoS_List, &
             Deriv_Partial_P, EoS_Cal,  EoS_Cal_Esd, FfCal_Dat, FfCal_Dat_Esd, FfCal_EoS, Init_EoS_Data_Type,  &
             Init_EoS_Type, Init_Err_EoS, Init_Eos_Thermal, Read_EoS_DataFile, Set_Eos_Names, &
             Set_Eos_Use, Set_Kp_Kpp_Cond, Write_Info_EoS

   !---- Definitions ----!

   !!----
   !!---- N_EOSPAR
   !!----    integer, public, parameter :: N_EosPAR=13
   !!----
   !!----    Specify the maximum number of Eos parameters allowed in Eos_type data type
   !!----
   !!---- Update: 08/07/2013 for thermal
   !!
   integer, public, parameter :: N_EOSPAR=13

   !!----
   !!---- NCOL_DATA_MAX
   !!----    integer, public, parameter :: ncol_data_max=22
   !!----
   !!----    Defines the maximum number of columns allowed in input data files
   !!----
   integer, public, parameter :: NCOL_DATA_MAX=22

   !!----
   !!---- ERR_EOS
   !!----    logical, public :: Err_EoS
   !!----
   !!----    Logical Variable indicating an error in CFML_EoS module
   !!----
   !!---- Update: January - 2013
   !!
   logical, public          :: ERR_EOS

   !!----
   !!---- ERR_EOS_MESS
   !!----    character(len=150), public :: ERR_EoS_Mess
   !!----
   !!----    String containing information about the last error
   !!----
   !!---- Update: January - 2013
   !!
   character(len=150), public :: ERR_EOS_MESS

   !!----
   !!----  TYPE :: EOS_TYPE
   !!--..
   !!----  Type, public :: EoS_Type
   !!----     character(len=80)                         :: Title     ! Descriptive title of EoS, set by user
   !!----     character(len=15)                         :: Model     ! Murnaghan, Birch-Murnaghan, Vinet, Natural-Strain
   !!----     character(len=15)                         :: TModel    ! Name for thermal model
   !!----     character(len=5), dimension(n_eospar)     :: ParName   ! Names of the Eos variables...init
   !!----     character(len=50),dimension(n_eospar)     :: Comment   ! Description of the Eos variables inclduing units...init
   !!----     integer                                   :: IModel    ! Index for Model
   !!----     integer                                   :: IOrder    ! Order for the Model
   !!----     logical                                   :: Linear    ! Flag for Linear EoS not volume
   !!----     integer                                   :: ITherm    ! Index for thermal expansion model, =0 for none
   !!----     integer,      dimension(n_eospar)         :: Iuse      ! Flags for parameters alllowed for a given EoS
   !!----     real(kind=cp)                             :: PRef      ! Pressure of Reference
   !!----     real(kind=cp)                             :: TRef      ! Temperature of Reference
   !!----     logical                                   :: TRef_fixed! If true, then Tref cannot be changed by user
   !!----     real(kind=cp), dimension(n_eospar)        :: Params    ! Eos Parameters
   !!----     real(kind=cp), dimension(n_eospar)        :: Esd       ! Sigma Eos Parameters
   !!----     real(kind=cp),                            :: X         ! Spare intensive variable, after P,T
   !!----     real(kind=cp)                             :: WChi2     ! weighted chi-squared for the refined parameters...set by ls
   !!----     real(kind=cp)                             :: DelPMax   ! Maximum misfit in P between obs and calc P...set by ls
   !!----     integer, dimension(4)                     :: IWt       ! Choice for weighting in LS: order is P,T,V,X, 1 = use sigmas for weights
   !!----     integer,      dimension(n_eospar)         :: IRef      ! Refinement switches
   !!----     real(kind=cp),dimension(n_eospar)         :: Factor    ! Scale factor for variables
   !!----     real(kind=cp),dimension(n_eospar)         :: Lastshift ! Shift applied in last LS cycle to parameters
   !!----     real(kind=cp),dimension(n_eospar,n_eospar):: VCV       ! Var-Covar matrix from refinement
   !!----  End Type EoS_Type
   !!----
   !!---- Update: January - 2013
   !!
   Type, public :: EoS_Type
      character(len=80)                         :: Title     ! Descriptive title of EoS, set by user
      character(len=15)                         :: Model     ! Murnaghan, Birch-Murnaghan, Vinet, Natural-Strain
      character(len=15)                         :: TModel    ! Name for thermal model
      character(len=5), dimension(n_eospar)     :: ParName   ! Names of the Eos variables...init
      character(len=50),dimension(n_eospar)     :: Comment   ! Description of the Eos variables inclduing units...init
      integer                                   :: IModel    ! Index for Model
      integer                                   :: IOrder    ! Order for the Model
      logical                                   :: Linear    ! Flag for Linear EoS not volume
      integer                                   :: ITherm    ! Index for thermal expansion model, =0 for none
      integer,      dimension(n_eospar)         :: Iuse      ! Flags for parameters allowed for a given EoS =0 (not), =1 (refineable), =2 (implied non-zero)
      real(kind=cp)                             :: PRef      ! Pressure of Reference
      real(kind=cp)                             :: TRef      ! Temperature of Reference
      logical                                   :: TRef_fixed! If true, then Tref cannot be changed by user
      real(kind=cp), dimension(n_eospar)        :: Params    ! EoS Parameters
      real(kind=cp), dimension(n_eospar)        :: Esd       ! Sigma EoS Parameters
      real(kind=cp)                             :: X         ! Spare intensive variable, after P,T
      real(kind=cp)                             :: WChi2     ! weighted chi-squared for the refined parameters...set by ls
      real(kind=cp)                             :: DelPMax   ! Maximum misfit in P between obs and calc P...set by ls
      integer, dimension(4)                     :: IWt       ! Choice for weighting in LS: order is P,T,V,X, 1 = use sigmas for weights
      integer,      dimension(n_eospar)         :: IRef      ! Refinement switches
      real(kind=cp),dimension(n_eospar)         :: Factor    ! Scale factor for variables
      real(kind=cp),dimension(n_eospar)         :: Lastshift ! Shift applied in last LS cycle to parameters
      real(kind=cp),dimension(n_eospar,n_eospar):: VCV       ! Var-Covar matrix from refinement
   End Type EoS_Type

   !!----
   !!---- TYPE :: EOS_LIST_TYPE
   !!--..
   !!---- Type, public :: EoS_List_Type
   !!----    integer                                   :: N   ! Number of EoS List
   !!----    type(EoS_Type),dimension(:),allocatable   :: EoS ! EoS Parameters
   !!---- End Type EoS_List_Type
   !!----
   !!---- Update: January - 2013
   !!
   Type, public :: EoS_List_Type
      integer                                 :: N      ! Number of EoS List
      type(EoS_Type),dimension(:),allocatable :: EoS    ! EoS Parameters
   End Type EoS_List_Type

   !!----
   !!----  TYPE :: EOS_DATA_TYPE
   !!--..
   !!----  Type, public :: EoS_Data_Type
   !!----     integer                     :: IUse    ! 0: No active, 1: active
   !!----     integer                     :: IGrp    ! Group
   !!----     real(kind=cp)               :: T       ! Temperature
   !!----     real(kind=cp)               :: P       ! Pressure
   !!----     real(kind=cp)               :: V       ! Volume
   !!----     real(kind=cp), dimension(3) :: cell    ! a,b,c parameters
   !!----     real(kind=cp), dimension(3) :: ang     ! alpha, beta, gamma parameters
   !!----     real(kind=cp)               :: SigT    ! Sigma Temperature
   !!----     real(kind=cp)               :: SigP    ! Sigma Pressure
   !!----     real(kind=cp)               :: SigV    ! Sigma Volume
   !!----     real(kind=cp), dimension (3):: sigC    ! Sigma for a, b and c parameters
   !!----     real(kind=cp), dimension (3):: sigA    ! Sigma for alpha, beta and gamma angles
   !!----  End Type EoS_Type
   !!----
   !!---- Update: January - 2013
   !!
   Type, public :: EoS_Data_Type
      integer                     :: IUse    ! 0=No active, 1= active
      integer                     :: IGrp    ! Group
      real(kind=cp)               :: T       ! Temperature
      real(kind=cp)               :: P       ! Pressure
      real(kind=cp)               :: V       ! Volume
      real(kind=cp), dimension(3) :: cell    ! a,b,c parameters
      real(kind=cp), dimension(3) :: ang     ! alpha, beta, gamma parameters
      real(kind=cp)               :: SigT    ! Sigma Temperature
      real(kind=cp)               :: SigP    ! Sigma Pressure
      real(kind=cp)               :: SigV    ! Sigma Volume
      real(kind=cp), dimension (3):: sigC    ! Sigma for a, b and c parameters
      real(kind=cp), dimension (3):: sigA    ! Sigma for alpha, beta and gamma angles
   End Type EoS_Data_Type

   !!----
   !!---- TYPE :: EOS_DATA_LIST_TYPE
   !!--..
   !!---- Type, public :: EoS_Data_List_Type
   !!----    character(len=80)                              :: Title     ! Title of dataset (normally from input datafile)
   !!----    character(len=40)                              :: System    ! Crystal System  (normally set by Def_Crystal_System)
   !!----    integer                                        :: N         ! Number of EoS Data List
   !!----    integer, dimension(ncol_data_max)              :: IC_Dat    ! Which values are input
   !!----    type(EoS_Data_Type),dimension(:),allocatable   :: EoSD      ! EoS Data Parameters
   !!---- End Type EoS_Data_List_Type
   !!----
   !!---- Update: January - 2013
   !!
   Type, public :: EoS_Data_List_Type
      character(len=80)                            :: Title     ! Title of dataset (normally from input datafile)
      character(len=40)                            :: System    ! Crystal System  (normally set by Def_Crystal_System)
      integer                                      :: N         ! Number of EoS Data List
      integer, dimension(ncol_data_max)            :: IC_Dat    ! Which values are input
      type(EoS_Data_Type),dimension(:),allocatable :: EoSD      ! EoS Data Parameters
   End Type EoS_Data_List_Type

Contains

   !!----
   !!---- FUNCTION Alpha_Cal(P,T,Eospar, DeltaT)
   !!----    real(kind=cp),            intent(in) :: P       ! Pressure
   !!----    real(kind=cp),            intent(in) :: T       ! Temperature
   !!----    type(Eos_Type),           intent(in) :: EoSPar  ! Eos Parameter
   !!----    real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T value
   !!----
   !!---- Calculate the alpha parameter in thermal case
   !!---- For linear case alpha is correct as 1/a da/dT
   !!----
   !!---- Date: 09/09/2013
   !!
   Function Alpha_Cal(P,T,Eospar, DeltaT) Result(alpha)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: P       ! Pressure
      real(kind=cp),            intent(in) :: T       ! Temperature
      type(Eos_Type),           intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T value

      !---- Local Variables ----!
      real(kind=cp)                  :: alpha,del,tt,delmin,tlimit
      real(kind=cp), dimension(-2:2) :: v       ! array for calc v values
      integer                        :: j

      !> Init
      alpha=0.0_cp
      tlimit=0.0_cp

      !> Need to trap numerical problems with Kroll, Salje, Pthermal at low T
      select case(eospar%itherm)
         case(0)                ! no thermal parameters
            return

         case(4,5,6)                ! Kroll and Pthermal: For T < 0.05T(Einstein) alpha=0. Same for Salje but test is 0.05Tsat
            if (t < 0.05_cp*eospar%params(11) ) return
      end select

      !> numerical solutions: set step size
      del =0.001_cp/eospar%params(10)      ! set step in T to get about 0.1% shift in V
      if (del > 30.0) del=30.0_cp          ! otherwise for small alpha, del is too big and alpha is inaccurate
      if (present(deltat)) del=deltat

      select case(eospar%itherm)           ! adjustment of step size
         case(2,3)                         ! T ok, but do not step down into invalid region
            if (abs(eospar%params(12)) > tiny(0.0_cp)) then
               delmin=abs(t-tlimit)/2.1_cp
               if (del > delmin) del=delmin
            end if

         case(4,5,6)
            ! if (t < 0.2_cp*eospar%params(11)) del=0.02_cp*eospar%params(11)
            delmin=(t-0.025_cp*eospar%params(11))/2.0_cp     ! do not allow step in to area where alpha=0
            if (del > delmin) del=delmin                     ! ensures T at all steps is positive
      end select

      !> do the numerical solution
      do j=-2,2,1
         tt=t+real(j)*del                                       ! apply shift to temp
         v(j)=get_volume(p,tt,eospar)                           ! calc resulting V
      end do

      alpha=(v(-2)+8.0_cp*(v(1)-v(-1))-v(2))/(12.0_cp*del)/v(0) ! derivative to second order approximation

      !> Trap non-physical results that arise due to rounding errors
      select case(eospar%itherm)
         case(4,5,6)
            if (alpha < tiny(0.0)) alpha=0.0_cp
      end select

      !> get_volume returns 'a' for linear, so alpha is correct as 1/a da/dT for linear

      return
   End Function Alpha_Cal

   !!----
   !!---- FUNCTION dKdT_Cal(P, T, Eospar, DeltaT)
   !!----    real(kind=cp),            intent(in) :: P       ! Pressure
   !!----    real(kind=cp),            intent(in) :: T       ! Temperature
   !!----    type(Eos_Type),           intent(in) :: EoSPar  ! Eos Parameter
   !!----    real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T value
   !!----
   !!---- Calculate the derivative dK/dt (or dM/dt in the linear case) at P and T
   !!----
   !!---- Date: 09/09/2013
   !!
   Function dKdT_Cal(P,T,Eospar,DeltaT) Result(dKdT)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: P       ! Pressure
      real(kind=cp),            intent(in) :: T       ! Temperature
      type(Eos_Type),           intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T

      !---- Local Variables ----!
      real(kind=cp) :: dKdT, vp, vm, del

      !> Init
      del=30.0_cp                      ! good number for accuracy
      if (present(deltat)) del=deltat

      !> Volume
      vp=get_volume(p,t+del,eospar)
      vm=get_volume(p,t-del,eospar)

      !> In the linear case: this function return dM/dT
      dKdT=(k_cal(vp,t+del,eospar)-k_cal(vm,t-del,eospar))/(2.0_cp*del)

      return
   End Function dKdT_Cal


   !!--++
   !!--++ FUNCTION Get_K0_T(T,EoSPar)
   !!--++    real(kind=cp),  intent(in) :: T       ! Temperature in K
   !!--++    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!--++
   !!--++ PRIVATE
   !!--++ Returns k0 needed for Eos calculations at T this means for pthermal,
   !!--++ k at Tref is returned.
   !!--++ In the linear case then  M(T, P=0) from params is returned
   !!--++
   !!--++ If k is calculated as being negative, an error message is set
   !!--++ and k0 is returned as the value at Tref for safety.
   !!--++
   !!--++ Date:09/09/2013
   !!
   Function Get_K0_T(T,Eospar) Result(k0)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: k0

      k0=eospar%params(2)

      select case(eospar%itherm)
          case(1:5)
              k0=k0+eospar%params(5)*(t-eospar%tref)
      end select

      return
   End Function Get_K0_T

   !!----
   !!---- FUNCTION Get_Pressure(V,T,EosPar)
   !!----    real(kind=cp),  intent(in) :: V       ! Volume
   !!----    real(kind=cp),  intent(in) :: T       ! Temperature
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Returns value of pressure (P) at (V,T) for EoS Model defined in EoSPar variable
   !!---- In any error occurs then the final value is set to 0.0
   !!----
   !!--.. Note:
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 22/02/2013
   !!--.. Changed 28/03/2013 RJA to use strain function to get f from Vo/V.
   !!--.. This will ensure better consistency
   !!----
   !!---- Date: 16/02/2013
   !!
   Function Get_Pressure(V,T,EosPar) Result(P)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                     :: p
      real(kind=cp)                     :: K0,Kp,kpp,vv0,x,u
      real(kind=cp)                     :: a,b,c,f
      real(kind=cp),dimension(n_eospar) :: ev
      real(kind=cp),dimension(3)        :: abc      ! for Tait parameters

      !> Init
      p=0.0_cp

      !> copy Eos parameters to local
      call EoS_to_Vec(eospar,ev)   ! Volume or linear case is covered

      !> Thermal case
      select case (eospar%itherm)
         case (0,6) ! No thermal case, or (6) thermal pressure which uses params at Tref
            vv0=v/eospar%params(1)      !vv0 is now V/V0 or a/a0
         case (1:5)
            vv0=v/Get_V0_T(T,EosPar)
      end select

      k0=Get_K0_T(T,eospar)             !handles thermal pressure case, returns K0 or M0
      if (eospar%linear) k0=k0/3.0_cp

      kp=ev(3)
      kpp=ev(4)
      f=strain(vv0,eospar)              ! use strain to get finite strain from v/v0 or a/a0

      !> Using volume dimensions
      vv0=1.0_cp/vv0                    ! vv0 now V0/V for ease of use in rest of subprogram
      if (eospar%linear) vv0=vv0**(3.0_cp)

      select case (eospar%imodel)
         case (1) !Murnaghan
            P=K0/Kp*(vv0**Kp - 1.0_cp)

         case (2) !Birch-Murnaghan
            a=f*(1.0_cp+2.0_cp*f)**(2.5_cp)     ! changed expressions to use only f 04/04/2013 RJA
            b=0.0_cp
            c=0.0_cp
            if (eospar%iorder > 2) b=1.5_cp*(Kp-4.0_cp)
            if (eospar%iorder == 4)c = 1.5_cp*(K0*Kpp + (Kp-4.0_cp)*(Kp-3.0_cp)+35.0_cp/9.0_cp)
            p=3.0_cp*K0*a*(1.0_cp + b*f + c*f*f)

         case (3) !Vinet
            x=vv0**(-1.0_cp/3.0_cp)
            u=1.0_cp -x
            p=3.0_cp*K0*u/x/x*exp(1.5_cp*(Kp-1.0_cp)*u)

         case (4) !Natural
            b=0.0_cp
            c=0.0_cp
            if (eospar%iorder > 2) b=1.5_cp*(Kp-2.0_cp)
            if (eospar%iorder == 4)c =1.5_cp*(K0*Kpp + 1.0_cp +(Kp-2.0_cp)+(Kp-2.0_cp)**2.0_cp)
            p=3.0_cp*vv0*K0*f*(1.0_cp + b*f + c*f*f)

         case (5) ! Tait
            call get_tait(eospar,t,abc)
            vv0=1.0_cp/vv0     ! back to vv0=v/v0
            p=(((vv0 +abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3)) - 1.0_cp)/abc(2)
      end select

      !> handle pthermal EoS
      if (eospar%itherm == 6 .and. eospar%imodel /= 0) p=p+pthermal(T,eospar)


      return
   End Function Get_Pressure

   !!----
   !!---- FUNCTION Get_Pressure_Esd(V,T,EosPar)
   !!----    real(kind=cp),  intent(in) :: V       ! Volume
   !!----    real(kind=cp),  intent(in) :: T       ! Temperature
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Gets the partial derivatives of P with respect to the EoS
   !!---- parameters at a given v,t point in order to calculate from vcv
   !!---- matrix the esd(p)
   !!----
   !!--.. Rewritten 09/04/2013 RJA to use deriv_partial
   !!----
   !!
   Function Get_Pressure_Esd(V,T,EosPar) Result(esd)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      integer                           :: i,j
      real(kind=cp)                     :: esd
      real(kind=cp),dimension(n_eospar) :: td
      real(kind=cp)                     :: vol,temp
      type(Eos_Type)                    :: E  ! Eos Parameter local copy

      !> Init
      esd=0.0_cp

      !> local copies
      vol=v
      temp=t
      e=eospar
      call Deriv_Partial_P(vol,temp,e,td)  ! gets all derivatives dP/dparam in array td

      !> esd
      do i=1,n_eospar
         do j=1,n_eospar
            esd=esd+eospar%vcv(i,j)*td(i)*td(j)
         end do
      end do

      !> Final
      esd=sqrt(esd)

      return
   End Function Get_Pressure_Esd

   !!----
   !!---- FUNCTION Get_Temperature(Vt,EosPar)
   !!----    real(kind=cp),  intent(in) :: Vt       ! Volume at temperature T
   !!----    type(Eos_Type), intent(in) :: EoSPar   ! Eos Parameter
   !!----
   !!---- Returns Temperature at P=0 for given V(T)
   !!----
   !!---- Date: 14/02/2013
   !!
   Function Get_Temperature(Vt,EosPar) Result(t)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: Vt      ! Volume at temperature T
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                     :: t
      real(kind=cp)                     :: V0,Tref,Alpha0,Alpha1,cvt
      real(kind=cp)                     :: c,delt,root
      real(kind=cp),dimension(n_eospar) :: ev

      !> Init
      t=0.0_cp
      Tref=eospar%tref

      !> Local copy Eospar
      call EoS_to_Vec(eospar,ev) ! Volume or linear case is covered

      !> This subroutine, if needed, needs complete rewriting!********
      v0=ev(1)
      alpha0=ev(10)
      alpha1=ev(11)

      !>
      cvt=vt
      if (eospar%linear) cvt=cvt**3.0_cp

      !> Check alpha0

      if (abs(alpha1) < 0.00001) then
         !> Only linear Thermal expansion
         delt=log(cvt/V0)/alpha0
         t=delt + TRef
      else
         c= -1.0_cp*( 0.5_cp*Alpha1*TRef**2.0_cp + Alpha0*TRef + Log(cvt/V0))
         Root=Alpha0*Alpha0 - 2.0_cp*Alpha1*c
         if(root < 0.0) root=0.0_cp
         t= (sqrt(root)-Alpha0)/Alpha1
      end if

      return
   End Function Get_Temperature

   !!--++
   !!--++ FUNCTION Get_V0_T(T,EosPar)
   !!--++    real(kind=cp),  intent(in) :: T        ! Temperature
   !!--++    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!--++
   !!--++ PRIVATE
   !!--++ Returns the volume at P=0 and T=T from V0 and Thermal expansion
   !!--++ Except for Pthermal, for which it returns V at P=0, T=Tref
   !!--++ Therfore this must remain PRIVATE
   !!--++
   !!--++ Date: 08/07/2013
   !!
   Function Get_V0_T(T,EosPar) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T        ! Temperature
      type(Eos_Type), intent(in) :: EoSPar   ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: V
      real(kind=cp)                      :: Tref,A,B,C,Tn,tt
      real(kind=cp)                      :: delt,delt2
      real(kind=cp), dimension(n_eospar) :: ev

      !> Init
      Tref=eospar%tref

      !> Local copy Eospar
      call EoS_to_Vec(eospar,ev) ! Volume or linear case is covered
                                 ! all equations written for volume
      delt=T-Tref
      select case(eospar%itherm)
         case(0)
            v=ev(1)                 ! no thermal eos: V is V0

         case(1)                     ! Berman, works at all T
            V=ev(1)*(1.0_cp +ev(10)*delt+0.5_cp*ev(11)*delt*delt)

         case(2)                     ! Fei,
            TT=T
            delt2=t*t-tref*tref
            if(tt < tiny(0._cp))tt=0.001        ! to prevent divide by zero

            if(abs(ev(12)) > tiny(0._cp))then
                V=ev(1)*exp(ev(10)*delt + 0.5_cp*ev(11)*delt2 - ev(12)*(1.0_cp/TT - 1.0_cp/Tref))
            else
                V=ev(1)*exp(ev(10)*delt + 0.5_cp*ev(11)*delt2) ! to protect from divide by zero when alpha2=0 and T=0
            endif

         case(3)                                    ! HP 1998 modified
            tt=t
            if(t < tiny(0._cp))tt=tiny(0._cp)      ! prevents sqrt(-ve T)
            V=ev(1)*(1.0_cp + ev(10)*(tt-tref) - 2.0_cp*(10.0_cp*ev(10)+ev(11))*(sqrt(tt) -sqrt(tref)))

         case(4)                                    ! Holland-Powell 2011 in the Kroll form
            Tn= ev(11)/Tref                        ! theta/Tref
            C=Tn*Tn*exp(Tn)/(exp(tn)-1.0_cp)**2.0_cp
            B=-1.0/ev(3)/(ev(3)+2.0_cp)
            if (t > 0.05_cp*ev(11)) then              ! to avoid numerical problems at T=0
               A=ev(10)*ev(11)/C *(1.0_cp/(exp(ev(11)/T)-1.0_cp) - 1.0_cp/(exp(Tn)-1.0_cp) )
            else
               A=ev(10)*ev(11)/C *(-1.0_cp/(exp(Tn)-1.0_cp) )          ! because when T=0 then 1.0/exp(Tein/T) = 0
            end if
            V=ev(1)*(-1.0_cp*ev(3) + (1.0_cp+ev(3))*(1.0_cp - ev(3)*(ev(3)+2.0_cp)*A/(ev(3)+1.0_cp))**B)

         case(5)                                     ! Salje, Tref fixed at zero
            A=T/ev(11)
            if (A < 0.001) then
               V=ev(1)                 ! ultra-low T: coth(theta_sat/T)=1.0
            else
               A=1.0_cp/tanh(1.0_cp/A) ! the coth(theta_sat/T)
               V=(ev(1)**(1.0_cp/3.0_cp) + ev(10)*ev(11)*(A-1.0_cp))**3.0_cp
            end if

         case(6)
            v=ev(1)         ! Pthermal needs V0 at Tref
      end select

      !> Linear
      if (eospar%linear) v=v**(1.0_cp/3.0_cp)

      return
   End Function Get_V0_T

   !!----
   !!---- FUNCTION Get_Volume(P,T,Eospar)
   !!----    real(kind=cp),  intent(in) :: P       ! Pressure
   !!----    real(kind=cp),  intent(in) :: T       ! Temperature
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Find volume from EoS at given P and T
   !!----
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 25/02/2013
   !!----
   !!---- Date: 16/02/13
   !!
   Function Get_Volume(P,T,EosPar) Result(v)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      integer                           :: nstep

      real(kind=cp), parameter          :: prec=0.000001_cp  !precision to find V.
      real(kind=cp)                     :: V
      real(kind=cp)                     :: V0,K0,Kp,k
      real(kind=cp)                     :: Vol, step,dp1,dp2
      real(kind=cp),dimension(n_eospar) :: ev
      real(kind=cp),dimension(3)        :: abc          ! Tait parameters
      real(kind=cp)                     :: pa           ! pa=p-pth

      !> Init
      v=0.0_cp
      pa=p            ! local copy p

      !> Local copy Eospar
      call EoS_to_Vec(eospar,ev) ! Volume or linear case is covered

      !> Set appropriate V0, and adjust for thermal pressure:
      select case (eospar%itherm)
         case (0)                                     ! 0=no thermal,
            v0=ev(1)                                  ! v0 is volume eos, (a0)^3 for linear

         case (1:5)
            v0=get_v0_t(t,eospar)                     ! returns a0 for linear
            if (eospar%linear) v0=v0**3.0_cp

         case (6)                                     ! thermal pressure
            v0=ev(1)
            pa=p-pthermal(t,eospar)                   ! adjust pressure to isothermal pressure for murn and tait
      end select

      !> set K0, for various purposes
      k0=Get_K0_T(T,eospar)                           !  for pthermal case tis returns k0 at Tref for volume,
      if (eospar%linear)k0=k0/3.0_cp                  ! correct M0 to K0 for linear
      kp=ev(3)

      !> Analytic solution for Murnaghan:  use this for first guess for other EoS except Tait
      v=v0*(1.0_cp + kp*pa/k0)**(-1.0_cp/kp)
      if (eospar%imodel ==1) then
         if (eospar%linear) v=v**(1.0_cp/3.0_cp)
         return
      end if

      !> analytic solution for Tait
      if (eospar%imodel ==5) then
         call get_tait(eospar,t,abc)                  ! get_tait returns volume-like parameters even for linear
         v=v0*(1.0_cp-abc(1)*(1.0_cp-(1.0_cp + abc(2)*pa)**(-1.0_cp*abc(3))))
         if (eospar%linear) v=v**(1.0_cp/3.0_cp)
         return
      end if

      vol=v
      if (eospar%linear) vol=vol**(1.0_cp/3.0_cp)

      !> Find iterative solution for the rest of functions: get_pressure includes the tehrmal pressure term
      dp1=p-get_pressure(vol,t,eospar)

      !> estimate the step to make in vol to get closer
      k=k0+p*kp                          ! Murn-like guess estimate to avoid recursive call back here when pthermal used
      if (eospar%linear)k=k*3.0_cp       ! The iteration is working with linear quantities
      step= -1.0_cp*vol*dp1/k            ! by definition of bulk modulus

      nstep=0
      do
         !> Trap infinite loops
         if (nstep > 10000) exit

         !> Increment Volume
         vol=vol+step
         dp2=p-get_pressure(vol,t,eospar)
         nstep=nstep+1

         !> Compare delta P's
         !if (dp1*dp2 < 0.0_cp) then
             !> overshot ptr:reverse and restimate step size
         !   step=-1.0_cp*dp2*vol/k
         !else
         !   if (abs(dp2) > abs(dp1)) step=-1.0*step      ! wrong direction: reverse
         !end if

         !> test for sufficient convergence
         !> switched to scaling by K because otherwise problems with small p
         !if (abs(dp1-dp2) < prec*abs(k) ) exit       ! dp1=dp2 so converged
         !if (abs(dp2) < prec*abs(k)) exit           ! dp smaller than required
         if (abs(step) < prec*Vol)exit          ! 1 part in 100,000 in volume

         !> not converged, so adjust step size
         if (dp1*dp2 < 0.0_cp) then
            !> overshot ptr:reverse step direction and make size smaller
            step=-0.5_cp*step
         else
            if (abs(dp2) > abs(dp1))then
               step=-1.0*step      ! wrong direction: reverse
            else
               step=0.9_cp*dp2/dp1*step   ! correct direction, should get a smaller step size
            end if                        ! the factor of 0.9 is a safety factor to ensure step gets smaller
         end if

         dp1=dp2        ! update delta-p values and go back for next cycle
      end do

      !> now set return value depending on success or not
      v=vol           ! success

      return
   End Function Get_Volume

   !!----
   !!---- FUNCTION Get_Volume_S(S,T,Eospar)
   !!----    real(kind=cp),  intent(in) :: S   ! Strain
   !!----    real(kind=cp),  intent(in) :: T   ! Temperature
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Returns the value of Volume obtained from Strain (S)
   !!----
   !!---- Date: 15/02/2013
   !!
   Function Get_Volume_S(S,T,Eospar) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S       ! Strain
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: V,v0
      real(kind=cp), dimension(n_eospar) :: ev

      !> Init
      v=0.0_cp

      !> Local Eos Copy
      call EoS_to_Vec(eospar,ev)

      !> Allow for thermal: s=function of V(P,T)/V(P=0,T)
      select case (eospar%itherm)
         case (0)
            v0=ev(1)

         case (1:6)
            v0=get_volume(0.0_cp,t,eospar)
            if (eospar%linear) v0=v0**3.0_cp
      end select

      select case (eospar%imodel)
         case (1,5) ! Murnaghan and Tait, no strain defined
            v=0.0_cp

         case (2) ! Birch-Murnaghan
            V=v0*(1.0_cp+2.0_cp*s)**(-1.5_cp)

         case (3) ! Vinet
            V=v0*(1.0_cp-s)**3.0_cp   !JGP

         case (4) ! Natural Strain
            V=v0*exp(-3.0_cp*s)
      end select

      !> Linear case
      if (eospar%linear) v=v**(1.0_cp/3.0_cp)

      return
   End Function Get_Volume_S

   !!----
   !!---- FUNCTION K_Cal(V,T,Eospar)
   !!----    real(kind=cp),  intent(in) :: V       ! Volume
   !!----    real(kind=cp),  intent(in) :: T       ! Temperature
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Returns value of K at this volume for EoS
   !!----
   !!--.. Changed code to use STRAIN function with v/v0: 27/02/2013
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- Date: 16/02/13
   !!
   Function K_Cal(V,T,Eospar) Result(kc)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: kc
      real(kind=cp) :: vv0,k0,kp,kpp,cvv0,p,vol1,f
      real(kind=cp) :: a,b,c,nu
      real(kind=cp), dimension(n_eospar) :: ev
      real(kind=cp), dimension(3)        :: abc     ! Tait parameters

      !> Init
      kc=0.0_cp

      !> Local EoS copy
      call EoS_to_Vec(eospar,ev)

      !> Get pressuure corrected for pthermal
      p=get_pressure(v,t,eospar)-pthermal(t,eospar)

      !> Set appropriate V/V0
      select case (eospar%itherm)
         case (0)                             ! No thermal
            vv0=v/eospar%params(1)            ! vv0 or aa0
         case (1:5)
            vv0=v/get_V0_T(t,eospar)          ! vv0 is  v(p,t)/v(p=0,t)
         case (6)
              vol1=get_volume(p,eospar%tref,eospar)    ! get volume at p-pth, at t=tref
              vv0=vol1/eospar%params(1)                ! vv0 is now for eoS at Tref
      end select

      !> Now get local copies of EoS parameters appropriate for EoS and thermal type
      k0=Get_K0_T(T,eospar)             ! returns M(T) for linear, handles pthermal
      if(err_eos)return                 ! exit with value eosparms(2) if k0 at P=0 calculated as negative

      if (eospar%linear) k0=k0/3.0_cp
      kp=ev(3)
      kpp=ev(4)

      !> Strain for BM, NS, Vinet EoS equations: this is f at p-pth if pthermal used
      f=strain(vv0,eospar)

      cvv0=vv0
      if (eospar%linear) cvv0=cvv0**3.0_cp

      select case (eospar%imodel)
         case (1)   !Murnaghan
            kc=k0+kp*p

         case (2)   !Birch-Murnaghan
            a=0.0_cp
            b=0.0_cp
            if (eospar%iorder > 2) a=1.5_cp*(kp-4.0_cp)
            if (eospar%iorder ==4) b =(9.0_cp*k0*kpp + (9.0_cp*kp*kp) - 63.0_cp*kp + 143.0_cp)/6.0_cp

            kc=k0*(1.0_cp+2.0_cp*f)**2.5_cp * (1.0_cp+ (7.0_cp+2.0_cp*a)*f + &
               (9.0_cp*a + 3.0_cp*b)*f*f + 11.0_cp*b*f*f*f)

         case (3) ! Vinet: from Schlosser & Ferrante PRB 37:4351; this form avoids pressure (RJA 26/2/2013)
                  ! new definition of f used 28/03/2013 RJA
            nu=1.5_cp*(kp-1)
            kc=k0*exp(nu*f) * (1.0_cp + (1.0_cp + nu)*f - nu*f*f) / (1.0_cp-f)**2.0_cp

         case (4) ! Natural strain: expanded from Poirier and Tarantola PEPI 109:1-8 (RJA 26/2/2013 - 28/03/2013)
            a=0.0_cp
            b=0.0_cp
            c=0.0_cp
            select case (eospar%iorder)
               case (2)
                  a=1.0_cp
               case (3)
                  a=kp-1.0_cp
                  b=1.5_cp*(kp-2.0_cp)
               case (4)
                  a=kp-1.0_cp
                  c=1.5_cp*(1.0_cp + k0*kpp + (kp-2.0_cp) + (kp-2.0_cp)**2.0_cp)
                  b=c+1.5_cp*(kp-2.0_cp)
            end select
            kc= 3.0_cp*k0/cvv0*( 1.0_cp/3.0_cp + a*f + b*f*f + c*f*f*f)

         case(5) ! Tait
             call get_tait(eospar,t,abc)              ! get_tait returns volume-like parameters even for linear
             kc=k0*(1.0_cp-abc(1)*(1.0_cp-(1.0_cp + abc(2)*p)**(-1.0_cp*abc(3))))*(1.0_cp + abc(2)*p)**(1.0_cp+abc(3))
      end select
      if (eospar%linear) kc=kc*3.0_cp

      return
   End Function K_Cal

   !!----
   !!---- FUNCTION Kp_Cal(V,T,EoSpar)
   !!----    real(kind=cp),  intent(in) :: V       ! Volume
   !!----    real(kind=cp),  intent(in) :: T       ! Temperature
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Returns value of kprime at this volume for EoS
   !!----
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- Date: 16/02/13
   !!
   Function Kp_Cal(V,T,EoSpar) Result (kpc)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: kpc
      real(kind=cp)                      :: vv0,k0,kp,kpp,vol1,p
      real(kind=cp)                      :: a,b,f,rkp_top, rkp_bot,nu
      real(kind=cp), dimension(n_eospar) :: ev
      real(kind=cp), dimension(3)        :: abc     ! Tait parameters

      !> Init
      kpc=0.0_cp

      !> Local copy
      call EoS_to_Vec(eospar,ev)

      !> Get pressuure corrected for pthermal
      p=get_pressure(v,t,eospar)-pthermal(t,eospar)

      !> Set appropriate V/V0  so that strain can be calculated
      select case (eospar%itherm)
         case (0)                             ! No thermal
            vv0=v/eospar%params(1)            ! vv0 or aa0
         case (1:5)
            vv0=v/get_V0_T(t,eospar)          ! vv0 is  v(p,t)/v(p=0,t)
         case (6)
            vol1=get_volume(p,eospar%tref,eospar)    ! get volume at p-pth, at t=tref
            vv0=vol1/eospar%params(1)                ! vv0 is now for eoS at Tref
      end select

      !> Now get local copies of EoS parameters appropriate for EoS and thermal type
      k0=Get_K0_T(T,eospar)             ! returns M(T) for linear, handles pthermal
      if (eospar%linear) k0=k0/3.0_cp
      kp=ev(3)
      kpp=ev(4)

      !> Strain for BM, NS, Vinet EoS equations: this is f at p-pth if pthermal used
      f=strain(vv0,eospar)

      select case(eospar%imodel)
         case (1) ! Murnaghan
            kpc=kp

         case (2) ! Birch-Murnaghan
            a=0.0_cp
            b=0.0_cp
            if (eospar%iorder > 2) a=1.5_cp*(kp-4.0_cp)
            if (eospar%iorder ==4) b = (9.0_cp*k0*kpp + (9.0_cp*kp*kp) - 63.0_cp*kp + 143.0_cp)/6.0_cp

            rkp_top= 12.0_cp + 2.0_cp*a + (49.0_cp+ 32.0_cp*a + 6.0_cp*b)*f + &
                     (81.0_cp*a + 60.0_cp*b)*f*f + 121.0_cp*b*f*f*f
            rkp_bot= 3.0_cp*(1.0_cp + (7.0_cp + 2.0_cp*a)*f + (7.0_cp*a + 3.0_cp*b)*f*f + 7.0_cp*b*f*f*f)
            kpc=rkp_top/rkp_bot

         case (3) ! Vinet: expression derived RJA 26-Feb-2013. Use new definition of f: RJA 28/03/2013
            nu=1.5_cp*(kp-1.0_cp)
            a= (1.0_cp + (1.0_cp + nu)*f - nu*f*f)
            kpc= (2.0_cp + nu*(1.0_cp-f) + (1.0_cp- f)*(1.0_cp + nu - 2.0_cp*nu*f)/a)/3.0_cp


         case (4) ! Natural strain: expanded from Poirier and Tarantola PEPI 109:1-8, derived RJA 26-Feb-2013
            a=0.0_cp
            b=0.0_cp
            if (eospar%iorder > 2) a=1.5_cp*(kp-2.0_cp)
            if (eospar%iorder ==4) b=1.5_cp*(1.0_cp + k0*kpp + (kp-2.0_cp) + (kp-2.0_cp)**2.0_cp)
            rkp_top=  1.0_cp + 2.0_cp/3.0_cp*a + 2.0_cp*(a+b)*f + 3.0_cp*b*f*f
            rkp_bot=  1.0_cp + (3.0_cp +2.0_cp*a)*f + 3.0_cp*(a+b)*f*f +3.0_cp*b*f*f*f
            kpc=1.0_cp+rkp_top/rkp_bot

          case(5) ! Tait
             call get_tait(eospar,t,abc)              ! get_tait returns volume-like parameters even for linear
             kpc=(kp+1.0_cp)*((1.0_cp + abc(2)*p)**abc(3)*(1.0_cp-abc(1)) + abc(1)) -1.0_cp
      end select
      if (eospar%linear) kpc=kpc*3.0_cp

      return
   End Function Kp_Cal

   !!----
   !!---- FUNCTION Kpp_Cal(V,T,EoSpar)
   !!----    real(kind=cp),  intent(in) :: V       ! Volume
   !!----    real(kind=cp),  intent(in) :: T       ! Temperature
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Returns value of kdouble-prime at this volume for EoS
   !!----
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- Date: 16/02/13
   !!
   Function Kpp_Cal(V,T,EoSpar) Result (kppc)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: kppc
      real(kind=cp) :: delv,vol1,vol2,rk1,rk2,p1,p2

      !> Init
      kppc=0.0_cp

      !> jump on equation type
      if (eospar%imodel == 1) return ! Murnaghan

      !> numerical solution
      delv=0.0001_cp*v
      vol1=v+delv
      vol2=v-delv
      rk1=kp_cal(vol1,t,eospar)
      rk2=kp_cal(vol2,t,eospar)
      p1=get_pressure(vol1,t,eospar)
      p2=get_pressure(vol2,t,eospar)

      !> trap case when p2=p1 because there is no eos loaded
      if(abs(p2-p1) .gt. tiny(0._cp))kppc=(rk2-rk1)/(p2-p1)


      !> No linear conversion is required because kp_cal returns values for "linear Kp" = Mp,
      !> so kppc is already dMp/dP = Mpp

      return
   End Function Kpp_Cal

   !!--++
   !!--++ FUNCTION NormPressure_Eos(S,T,EosPar)
   !!--++    real(kind=cp),  intent(in) :: S       ! Strain
   !!--++    real(kind=cp),  intent(in) :: T       ! Temperature
   !!--++    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!--++
   !!--++ PRIVATE
   !!--++ Returns the value of Normalised Pressure (F) at this Strain (S) using
   !!--++ the EoS parameters as K0, V0, Kp, Kpp
   !!--++
   !!--++ Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!--++ Modified for thermal EoS: requires T to be meaningful
   !!--++
   !---++ Date: 10/09/2013
   !!
   Function NormPressure_Eos(S, T, EosPar) Result(F)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S       ! Strain
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                     :: F
      real(kind=cp)                     :: k0,kp,kpp, b,c,v0
      real(kind=cp),dimension(n_eospar) :: ev

      !> Init
      f=0.0_cp

      !> Local copy
      call EoS_to_Vec(Eospar,Ev)

      !> Get correct parameters for this T
      select case(eospar%itherm)
         case(0)
            k0=ev(2)
            kp=ev(3)
            kpp=ev(4)

         case(1:5)     ! normal thermal expansion models with dK/dT
            k0=Get_K0_T(T,eospar)                     ! returns M(T) for linear,
            if (eospar%linear) k0=k0/3.0_cp
            kp=ev(3)
            kpp=ev(4)

         case(6)       ! Pthermal: cannot use get_V0_T or get_K0_T because they return V and K at Tref for pthermal
            v0=get_volume(0.0_cp,T,eospar)        ! determine V0 at this T
             k0=k_cal(v0,T,eospar)
             kp=kp_cal(v0,T,eospar)
             kpp=kpp_cal(v0,T,eospar)
             if (eospar%linear)then
                k0=k0/3.0_cp
                kp=kp/3.0_cp
                kpp=kpp/3.0_cp
             end if
      end select

      select case(eospar%imodel)
         case (1,5) ! Murnaghan, Tait
            f=0.0_cp

         case (2) ! Birch-Murnaghan
            b=0.0_cp
            c=0.0_cp
            if (eospar%iorder > 2)b=1.5_cp*(kp-4.0_cp)
            if (eospar%iorder ==4)c=1.5_cp*(k0*kpp + (kp-4.0_cp)*(kp-3.0_cp)+35.0_cp/9.0_cp)
            f=K0*(1.0_cp + b*s + c*s*s)

         case (3) ! Vinet: new definition RJA 28/03/2013
            f= K0*exp(1.5_cp*(Kp-1.0_cp)*s)

         case (4) ! Natural Strain
            b=0.0_cp
            c=0.0_cp
            if (eospar%iorder > 2)b=1.5_cp*(Kp - 2.0_cp)
            if (eospar%iorder ==4)c=1.5_cp*(K0*Kpp + 1.0_cp + (Kp-2.0_cp)+(Kp-2.0_cp)**2.0_cp)
            f=K0*(1.0_cp + b*s + c*s*s)
      end select

      return
   End Function NormPressure_Eos

   !!--++
   !!--++ FUNCTION NormPressure_P(S,P,imodel)
   !!--++    real(kind=cp),  intent(in) :: S      ! Strain
   !!--++    real(kind=cp),  intent(in) :: P      ! Pressure
   !!--++    integer      ,  intent(in) :: imodel ! Type of EoS
   !!--++
   !!--++ PRIVATE
   !!--++ Returns the value of Normalised Pressure (F)  at this Strain (S) and Pressure (P)
   !!--++ for the type of Eos indicated by imodel
   !!--++
   !!--++ NOTE:
   !!--++    The eos paramaters  are NOT used in this calculation
   !!--++    To get Normalised Pressure (F) from Strain (S) and Eos parameters,
   !!--++    use Function NormPressure_Eos(S,EosPar)
   !!--++
   !!--++ Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!--++
   !!--++ 13/3/2013  Changed default value of F to zero and trapped s=0 for which F is not defined: RJA
   !!--++ 05/09/2013 Changed argument from eospar to just imodel
   !!--++
   !!--++ Date: 15/02/2013
   !!
   Function NormPressure_P(S,P,imodel) Result(F)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S      ! Strain
      real(kind=cp),  intent(in) :: P      ! Presure
      integer      ,  intent(in) :: imodel ! type of eos

      !---- Local Variables ----!
      real(kind=cp) :: f

      !> Init: a zero value is returned for undefined values
      f=0.0_cp

      !> If the finite strain is zero, the normalised pressure is not defined
      if (s > tiny(0.0) ) then
         select case (imodel)
            case (1,5) ! Murnaghan, Tait
               f=0.0_cp

            case (2) ! Birch-Murnaghan
               f=p/3.0_cp/s/(1.0_cp+2.0_cp*s)**2.5_cp

            case (3) ! Vinet: new definition RJA 28/03/2013
               f= (p*(1.0_cp - s)**2.0_cp) / 3.0_cp / s

            case (4) ! Natural Strain
               f=p/3.0_cp/s * exp(-3.0_cp*s)
         end select
      end if

      return
   End Function NormPressure_P

   !!----
   !!---- FUNCTION Pressure_F(F,S,EosPar)
   !!----    real(kind=cp),  intent(in) :: F       ! Normalize Pressure
   !!----    real(kind=cp),  intent(in) :: S       ! Strain
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Returns the value of Pressure (P) from the input Normalized Pressure(F)
   !!---- and Strain (S) according to the EoS model in EosPar.
   !!----
   !!--.. NOTE:
   !!--..    The Eos parameters are NOT used in this calculation, on the
   !!--..    eospar%imodel (ie type of Eos)
   !!----
   !!---- Date: 21/02/2013
   !!
   Function Pressure_F(F,S,EosPar) Result(P)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: F       ! Normalized Pressure
      real(kind=cp),  intent(in) :: S       ! Strain
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: p, cF

      !> Init
      p=0.0_cp
      cf=f

      select case (eospar%imodel)
         case (1,5) ! Murnaghan, Tait

         case (2) ! Birch-Murnaghan
            p=3.0_cp*f*s*(1.0_cp +2.0_cp*s)**2.5_cp

         case (3) ! Vinet
            p=(3.0_cp*(1.0_cp-s)/s)*(10.0_cp)**f

         case (4) ! Natural
            p=3.0_cp*f*s*exp(3.0_cp*s)
      end select

      return
   End Function Pressure_F

   !!--++
   !!--++ FUNCTION Pthermal(T,EosPar)
   !!--++    real(kind=cp),  intent(in) :: T       ! Temperature
   !!--++    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!--++
   !!--++ PRIVATE
   !!--++  Calculate Pthermal from eosparameters at temperature T
   !!--++
   !!--++ Date: 10/09/2013
   !!
   Function Pthermal(T,EosPar) Result(Pth)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: pth,thtref,exp0,eta0
      real(kind=cp),dimension(n_eospar) :: ev

      !> Local copy
      call eos_to_vec(eospar,ev)    !handle linear case

      select case (eospar%itherm)
         case (6) ! Thermal pressure from Holland and Powelll
            thtref=ev(11)/eospar%tref         ! T_einstein/Tref
            exp0=exp(thtref)                  ! exp(T_Ein/Tref)
            eta0= thtref*thtref*exp0/(exp0-1.0_cp)**2.0_cp  ! eta at Tref

            pth = ev(10)*ev(2)*ev(11)/eta0*( 1.0_cp/(exp(ev(11)/t)-1.0_cp) -1.0_cp/(exp0 -1.0_cp))

         case default
            pth=0.0_cp
      end select

      return
   End Function Pthermal

   !!----
   !!---- FUNCTION Strain(VV0,EosPar)
   !!----    real(kind=cp),  intent(in) :: VVo   ! Volume divided by V0 at this temperature
   !!----    type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
   !!----
   !!---- Returns the value of Strain (S) at this V according with the EoS model
   !!----
   !!--.. NOTE:
   !!--..    The values of EosPar are NOT used. Only the type of Eos function is needed!
   !!--..    For linear eos, the strain is that of a^3
   !!--..    If V/V0 = 0. then an error is set
   !!----
   !!---- Date: 15/02/2013
   !!
   Function Strain(VV0,EosPar) Result(S)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: VV0     ! Volume divided by V0 at this temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: s
      real(kind=cp) :: cvv0

      !> Init
      s=0.0_cp
      if (vv0 <= 0.00001) then
         Err_EoS=.true.
         Err_EoS_Mess="Strain calculation called with invalid V/V0 =< 0: strain set zero"
         return
      end if

      !> Local copy
      cvv0=vv0
      if (eospar%linear) cvv0=vv0**3.0_cp

      select case (eospar%imodel)
         case (1,5) ! Murnaghan, Tait
            s=0.0_cp

         case (2) ! Birch-Murnaghan
            s=(cvv0**(-2.0_cp/3.0_cp) - 1.0_cp)/2.0_cp

         case (3) ! Vinet: new definition RJA 28/03/2013
            s= 1.0_cp - cvv0**(1.0_cp/3.0_cp)

         case (4) ! Natural Strain: the original definition of strain was wrong in v5 by a change in sign
            s= -1.0_cp*log(cvv0)/3.0_cp
      end select

      return
   End Function Strain

   !!----
   !!---- FUNCTION Strain_EoS(V,T,EosPar)
   !!----    real(kind=cp),  intent(in) :: V          ! Volume
   !!----    real(kind=cp),  intent(in) :: T          ! Temperature
   !!----    type(Eos_Type), intent(inout) :: EoSPar  ! Eos Parameters
   !!----
   !!---- Returns the value of Strain (S) at this V and T according with the EoS parameters
   !!----
   !!---- Date: 05/09/2013
   !!
   Function Strain_EOS(V,T,EosPar) Result(S)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameters

      !---- Local Variables ----!
      real(kind=cp) :: vvo,s
      type(Eos_Type):: E  ! Eos Parameters copy

      !> Init
      s=0.0_cp
      if (v <= 0.00001) then
         Err_EoS=.true.
         Err_EoS_Mess="Strain calculation called with invalid V =< 0: strain set zero"
         return
      end if
      e=eospar

      !> Calculate
      vvo=v/get_volume(0.0_cp,t,e)     ! vv0 is V(P,T)/V(0,T) or a(P,T)/a(0,T)
      s=strain(vvo,e)                  ! cubes vv0 on input if linear

      return
   End Function Strain_EOS

   !>----
   !>---- SUBROUTINE ZONE
   !>----

   !!----
   !!---- SUBROUTINE Allocate_EoS_Data_List(N, E)
   !!----    integer,                    intent(in)      :: N    !  In  -> Number of elements of E
   !!----    type (eos_data_list_type),  intent(in out)  :: E    !  In  -> Objet to be allocated
   !!----
   !!----    Allocation of objet E of eos_list_data_type.
   !!----    This subroutine should be called before using an object of type eos_data_list.
   !!----
   !!---- Update: 16/02/2013
   !!
   Subroutine Allocate_EoS_Data_List(N, E)
      !---- Arguments ----!
      integer,                    intent(in)       :: N  ! Number of elements of E
      type (eos_data_list_type),  intent(in out)   :: E  ! Objet to be allocated

      !---- Local Variables ----!
      integer :: i,ier

      !> Set dimension
      e%n = n
      if (allocated(e%eosd)) deallocate(e%eosd)
      allocate (e%eosd(n),stat=ier)
      if (ier /= 0) then
         e%n = 0
         err_eos=.true.
         err_eos_mess="Problems allocating memory for Eos_Data_List_Type variable"
         return
      end if

      !> Initializing
      do i=1,n
         call init_eos_data_type(e%eosd(i))
      end do

      return
   End Subroutine Allocate_EoS_Data_List

   !!----
   !!---- SUBROUTINE Allocate_EoS_List(N, E)
   !!----    integer,               intent(in)      :: N    !  In  -> Number of elements of E
   !!----    type (eos_list_type),  intent(in out)  :: E    !  In  -> Objet to be allocated
   !!----
   !!----    Allocation of objet E of eos_list_type.
   !!----    This subroutine should be called before using an object of type eos_list.
   !!----
   !!---- Update: 16/02/2013
   !!
   Subroutine Allocate_EoS_List(N, E)
      !---- Arguments ----!
      integer,               intent(in)       :: N  ! Number of elements of E
      type (eos_list_type),  intent(in out)   :: E  ! Objet to be allocated

      !---- Local Variables ----!
      integer :: i,ier

      !> Set dimension
      e%n = n
      if (allocated(e%eos)) deallocate(e%eos)
      allocate (e%eos(n),stat=ier)
      if (ier /= 0) then
         e%n = 0
         err_eos=.true.
         err_eos_mess="Problems allocating memory for Eos_List_Type variable"
         return
      end if

      !> Initializing
      do i=1,n
         call init_eos_type(e%eos(i))
      end do

      return
   End Subroutine Allocate_EoS_List

   !!----
   !!---- SUBROUTINE Deallocate_EoS_Data_List(E)
   !!----    type (eos_data_list_type), intent(in out)   :: E  ! In/ Out -> Objet to be deallocated
   !!----
   !!----    De-allocation of objet E of type eos_data_list.
   !!----    This subroutine should be after using an object of type eos_data_list that is no
   !!----    more needed.
   !!----
   !!---- Update: 16/02/2013
   !!
   Subroutine Deallocate_EoS_Data_List(E)
      !---- Arguments ----!
      type (eos_data_list_type), intent(in out)   :: E  ! Objet to be deallocated

      if (allocated(E%eosd)) deallocate (E%eosd)
      E%n=0

      return
   End Subroutine Deallocate_EoS_Data_List

   !!----
   !!---- SUBROUTINE Deallocate_EoS_List(E)
   !!----    type (eos_list_type), intent(in out)   :: E  ! In/ Out -> Objet to be deallocated
   !!----
   !!----    De-allocation of objet E of type eos_list.
   !!----    This subroutine should be after using an object of type eos_list that is no
   !!----    more needed.
   !!----
   !!---- Update: 16/02/2013 15:49:04
   !!
   Subroutine Deallocate_EoS_List(E)
      !---- Arguments ----!
      type (eos_list_type), intent(in out)   :: E  ! Objet to be deallocated

      if (allocated(E%eos)) deallocate (E%eos)
      E%n=0

      return
   End Subroutine Deallocate_EoS_List

   !!--++
   !!--++ SUBROUTINE Def_Crystal_System(dat)
   !!--++    type (eos_data_list_type),  intent(in out)   :: dat   !data structure
   !!--++
   !!--++ PRIVATE
   !!--++ Either sets cell parameters to conform to specifed system Or,
   !!--++ if no system specified, tries to determine system if all cell parameters
   !!--++ provided
   !!--++
   !!--++ Update: 18/04/2013
   !!
   Subroutine Def_Crystal_System(dat)
      !---- Arguments ----!
      type (eos_data_list_type),  intent(in out)   :: dat  ! data structure

      !---- Local Variables ----!
      character(len=40)       :: Family, SystemC, system
      character(len=1)        :: Symbol
      integer                 :: i,ndat
      type(Crystal_Cell_Type) :: ncell

      !> Check
      ndat=dat%n
      if (ndat <= 0) return

      !> local copies from data construct
      system=dat%system

      if (len_trim(system) > 0) then
         system=adjustl(system)
         select case (u_case(system(1:4)))
         case ('MONO')
               if(index(u_case(system),' C ') /= 0)then
               !> alpha = beta = 90º
                   dat%eosd(1:ndat)%ang(1)=90.0
                   dat%eosd(1:ndat)%ang(2)=90.0
                   dat%eosd(1:ndat)%siga(1)=0.0
                   dat%eosd(1:ndat)%siga(2)=0.0
               else
               !> alpha = gamma = 90º
                   dat%eosd(1:ndat)%ang(1)=90.0
                   dat%eosd(1:ndat)%ang(3)=90.0
                   dat%eosd(1:ndat)%siga(1)=0.0
                   dat%eosd(1:ndat)%siga(3)=0.0
               endif
            case ('ORTH')
               !> Angles =90º
               dat%eosd(1:ndat)%ang(1)=90.0
               dat%eosd(1:ndat)%ang(2)=90.0
               dat%eosd(1:ndat)%ang(3)=90.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

            case ('TETR')
               !> Angles =90º
               dat%eosd(1:ndat)%ang(1)=90.0
               dat%eosd(1:ndat)%ang(2)=90.0
               dat%eosd(1:ndat)%ang(3)=90.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

               !> a=b
               dat%eosd(1:ndat)%cell(2)=dat%eosd(1:ndat)%cell(1)
               dat%eosd(1:ndat)%sigc(2)=dat%eosd(1:ndat)%sigc(1)

            case ('TRIG','HEXA')
               !> Angles alpha=beta=90º, gamma=120º
               dat%eosd(1:ndat)%ang(1)= 90.0
               dat%eosd(1:ndat)%ang(2)= 90.0
               dat%eosd(1:ndat)%ang(3)  =120.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

               !> a=b
               dat%eosd(1:ndat)%cell(2)=dat%eosd(1:ndat)%cell(1)
               dat%eosd(1:ndat)%sigc(2)=dat%eosd(1:ndat)%sigc(1)

            case ('RHOM')
               !> alpha = beta = gamma
               dat%eosd(1:ndat)%ang(2)=dat%eosd(1:ndat)%ang(1)
               dat%eosd(1:ndat)%ang(3)=dat%eosd(1:ndat)%ang(1)
               dat%eosd(1:ndat)%siga(2)=dat%eosd(1:ndat)%siga(1)
               dat%eosd(1:ndat)%siga(3)=dat%eosd(1:ndat)%siga(1)

               !> a = b = c
               dat%eosd(1:ndat)%cell(2)=dat%eosd(1:ndat)%cell(1)
               dat%eosd(1:ndat)%cell(3)=dat%eosd(1:ndat)%cell(1)
               dat%eosd(1:ndat)%sigc(2)=dat%eosd(1:ndat)%sigc(1)
               dat%eosd(1:ndat)%sigc(3)=dat%eosd(1:ndat)%sigc(1)

            case ('CUBI')
               !> Angles =90º
               dat%eosd(1:ndat)%ang(1)=90.0
               dat%eosd(1:ndat)%ang(2)=90.0
               dat%eosd(1:ndat)%ang(3)=90.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

               !> a=b=c
               dat%eosd(1:ndat)%cell(2)=dat%eosd(1:ndat)%cell(1)
               dat%eosd(1:ndat)%cell(3)=dat%eosd(1:ndat)%cell(1)
               dat%eosd(1:ndat)%sigc(2)=dat%eosd(1:ndat)%sigc(1)
               dat%eosd(1:ndat)%sigc(3)=dat%eosd(1:ndat)%sigc(1)
         end select

      else
         !> checking cell parameters
         system=' '
         if (any(dat%eosd(1)%cell <= 0.1)) return
         if (any(dat%eosd(1)%ang <= 0.5)) return

         !> Obtaining system from cell parameters
         ncell%cell=dat%eosd(1)%cell
         ncell%ang =dat%eosd(1)%ang
         call Get_Cryst_Family(ncell,Family,Symbol,System)
         system=adjustl(system)
         do i=2,ndat
            ncell%cell=dat%eosd(i)%cell
            ncell%ang =dat%eosd(i)%ang
            call Get_Cryst_Family(ncell,Family,Symbol,SystemC)
            systemC=adjustl(systemc)
            if (systemC /= System) then
                err_eos=.true.
                Err_EoS_Mess="Cells values in the Data file are inconsistent"
                system=' '
            end if
         end do
      end if

      return
   End Subroutine Def_Crystal_System

   !!----
   !!---- SUBROUTINE Deriv_Partial_P(V,T,Eospar,td)
   !!----    real(kind=cp),                      intent(in)  :: V       ! Volume
   !!----    real(kind=cp),                      intent(in)  :: T       ! Temperature
   !!----    type(Eos_Type),                     intent(in)  :: Eospar  ! Eos Parameter
   !!----    real(kind=cp), dimension(n_eospar), intent(out) :: td      ! derivatives dP/d(param)
   !!----
   !!---- Calculate Partial derivates of Pressure respect to Params
   !!----
   !!---- Date: 10/09/2013
   !!
   Subroutine Deriv_Partial_P(V,T,Eospar,td)
      !---- Arguments ----!
      real(kind=cp),                      intent(in)  :: V       ! Volume
      real(kind=cp),                      intent(in)  :: T       ! Temperature
      type(Eos_Type),                     intent(in)  :: Eospar  ! Eos Parameter
      real(kind=cp), dimension(n_eospar), intent(out) :: td      ! derivatives dP/d(param)

      !---- Local Variables ----!
      real(kind=cp), dimension(n_eospar) :: tda,tdn                ! analytic and numeric derivatives

      !> Calculate derivatives by both methods: correct values are returned for linear
      call Deriv_Partial_P_Analytic(V,T,Eospar,tda)
      call Deriv_Partial_P_Numeric(V,T,Eospar,tdn)

      !> Init
      td=0.0_cp

      !> now decide which ones to use
      td(1:5)=tda(1:5)                   ! analytic for Vo and moduli terms because these are exact even at small P
      td(6:n_eospar)=tdn(6:n_eospar)

      return
   End Subroutine Deriv_Partial_P

   !!--++
   !!--++ SUBROUTINE Deriv_Partial_P_Analytic(V,T,Eospar,td)
   !!--++    real(kind=cp),                      intent(in)  :: V       ! Volume
   !!--++    real(kind=cp),                      intent(in)  :: T       ! Temperature
   !!--++    type(Eos_Type),                     intent(in)  :: Eospar  ! Eos Parameter
   !!--++    real(kind=cp), dimension(n_eospar), intent(out) :: td      ! derivatives dP/d(param)
   !!--++
   !!--++ PRIVATE
   !!--++ Calculates the partial derivatives of P with respect to the EoS
   !!--++ at a given v,t point, and returns  them in array td
   !!--++
   !!--.. 27-Feb-2013: RJA edits to use IREF as in original
   !!--++
   !!--++ Date: 10/09/2013
   !!
   Subroutine Deriv_Partial_P_Analytic(V,T,Eospar,td)
      !---- Arguments ----!
      real(kind=cp),                      intent(in)  :: V       ! Volume
      real(kind=cp),                      intent(in)  :: T       ! Temperature
      type(Eos_Type),                     intent(in)  :: Eospar  ! Eos Parameter
      real(kind=cp), dimension(n_eospar), intent(out) :: td      ! derivatives dP/d(param)

      !---- Local Variables ----!
      real(kind=cp), dimension(n_eospar) :: ev
      real(kind=cp)                      :: vv0,k0,delt, delt2,vt0,dpdvt,deltinv
      real(kind=cp)                      :: f,b,c,cv
      real(kind=cp)                      :: term1,term2,term3,vterm
      real(kind=cp)                      :: a,f52,da,db,dc,nu   ! new guys
      real(kind=cp),dimension(3)         :: abc                 !Tait parammeters

      !> Local copy
      call EoS_to_Vec(eospar,ev)    ! Volume or linear case is covered: ev contains volume-like parameters and
                                    ! all derivatives will be calculated as volume-like derivatives
      cv=v                          ! cv has the current volume or 'a' input. It is needed for dP/dVo

      select case (eospar%Itherm)
        case (0)
            vt0=eospar%params(1)                ! vt0 is the V (or a) at this T and P=0

         case (1:)
            vt0=get_V0_T(t,eospar)
      end select
      vv0=v/vt0                                 ! vv0 is V(P,T)/V(P,T=0), or the linear equivalent

      k0=Get_K0_T(T,eospar)                     ! returns M(T) for linear
      if (eospar%linear) k0=k0/3.0_cp

      !> now get the finite strain for vv0
      f=strain(vv0,eospar)                      ! use strain to get finite strain from v/v0 or a/a0

      !> Using volume dimensions for remainder of calculations
      vv0=1.0_cp/vv0                               ! vv0 now V0/V for ease of use in rest of subprogram
      if (eospar%linear)then
         vt0=vt0**3.0_cp           ! Vo at this T
         vv0=vv0**3.0_cp           ! V/Vo at this T
         cv=cv**3.0_cp             ! current volume at this P,T
      end if

      td=0.0_cp
      select case (eospar%imodel)
         case (1) ! Murnaghan: validated algebra, 10/04/2013
            td(1)=vv0**(ev(3)-1.0_cp)*k0/cv
            td(2)=(vv0**ev(3)-1.0_cp)/ev(3)
            td(3)=k0/ev(3)*(vv0**ev(3)*(log(vv0)-1.0_cp/ev(3)) +1.0_cp/ev(3))
            td(4)=0.0_cp

         case (2) ! Birch-Murnaghan: new code 10-11/04/2013 Validated
            !> Assign coefficients of f in expansion of P: use K0 as that is K at P=0
            !  f was calculated above at init
            a=K0                           ! P=3(1+2f)^(5/2) . (af + bf^2 +cf^3)
            b=0.0_cp
            c=0.0_cp
            f52=(1.0_cp+2.0*f)**2.5_cp
            if (eospar%iorder > 2)b=1.5_cp*K0*(ev(3)-4.0_cp)
            if (eospar%iorder == 4)c = 1.5_cp*K0*(K0*ev(4) + (ev(3)-4.0_cp)*(ev(3)-3.0_cp)+35.0_cp/9.0_cp)

            !> dP/dVo
            td(1)=f52/vt0 * ( a + (7.0_cp*a+2.0_cp*b)*f + (9.0_cp*b+3.0_cp*c)*f*f + 11.0_cp*c*f*f*f)

            !> dP/dKo:  da, db, dc are da/dparam etc for each param K0, Kp0, Kpp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder > 2)db = 1.5_cp*(ev(3)-4.0_cp)
            if (eospar%iorder == 4)dc = 3.0_cp*K0*ev(4) + 1.5_cp*((ev(3)-4.0_cp)*(ev(3)-3.0_cp)+35.0_cp/9.0_cp)

            td(2)= 3.0*f52*f* (1.0 + db*f + dc*f*f)

            !> dP/dKp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder > 2)db = 1.5_cp*K0
            if (eospar%iorder == 4)dc = 1.5_cp*K0*(2.0_cp*ev(3)-7.0_cp)

            td(3) = 3.0*f52* (db + dc*f) *f*f

            !> dP/dKpp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder == 4)then
                dc = 1.5_cp*K0*K0
                td(4) = 3.0*f52*dc*f*f*f
            else
                td(4)=0.0_cp
            endif

         case (3) ! Vinet: new code and validated 11/04/2013
            if (eospar%iorder == 2) then
               td(1) = k0/vt0 *(1.0_cp+f)/(1.0_cp-f)**2.0_cp
               td(2) =      3.0_cp * f /(1.0_cp-f)**2.0_cp
               td(3) = 0.0_cp
               td(4) = 0.0_cp

            else
               nu=1.5_cp*(ev(3)-1.0)
               td(1) = k0/vt0 * (1.0_cp + (1.0_cp+nu)*f - nu*f*f)/(1.0_cp-f)**2.0_cp * exp(nu*f)
               td(2) =      3.0_cp * f /(1.0_cp-f)**2.0_cp * exp(nu*f)
               td(3) = k0 * 4.5_cp * f*f/(1.0_cp-f)**2.0_cp * exp(nu*f)
               td(4) = 0.0_cp
            end if

         case (4) ! Natural Strain: new code and validated 11/04/2013
            !> Assign coefficients of f in expansion of P: use K0 as that is K at P=0,T=Tdata
            !  f was calculated above at init
            !
            !  the coefficients a and b are for P= 3 K0 (V0/V) f(1+ af + bf^2)
            a=0.0_cp
            b=0.0_cp
            if(eospar%iorder > 2) a=1.5_cp*(ev(3)-2.0_cp)
            if(eospar%iorder ==4) b=1.5_cp*(1.0_cp + K0*ev(4) + (ev(3)-2.0_cp) + (ev(3)-2.0_cp)**2.0_cp)

            !> dP/dVo
            td(1) = K0/cv * (1.0_cp  + (2.0_cp*a+3.0_cp)*f + 3.0_cp*(a+b)*f*f +3.0_cp*b*f*f*f)

            !> d(p)/d(k0)
            td(2) = 3.0_cp * vv0 * (f + a*f*f + (b + 1.5_cp*ev(4)*k0)*f*f*f)

            !> dp/dKp0
            da = 0.0_cp
            db = 0.0_cp
            if(eospar%iorder > 2) da=1.5_cp
            if (eospar%iorder == 4)db = 3.0_cp*ev(3) - 4.5_cp

            td(3) = 3.0_cp * k0 * vv0 * (da + db*f)*f*f

            !> dp/dKpp0
            td(4)=0.0_cp
            if (eospar%iorder == 4)td(4) = 4.5_cp* k0*k0 * vv0 *f*f*f

         case(5) ! Tait: new code and validated against Murnaghan and against finite diffs 15/06/2013
            call get_tait(eospar,t,abc)
            vv0=1.0/vv0        ! now vv0 is v/v0

            !> dP/dVo
            td(1)= k0 * vv0**2.0_cp/cv * ( (vv0 + abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3) -1.0_cp)

            !> d(p)/d(k0)
            td(2)=abc(1)*abc(3)*(( (vv0 + abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3)) -1.0_cp)

            !> dp/dKp0
            if (eospar%iorder > 2)then
               da= k0*ev(4)/(1.0_cp +ev(3)+k0*ev(4))**2.0_cp      ! da/dKp
               db= 1.0_cp/k0 + ev(4)/(1.0_cp+ev(3))**2.0_cp       ! db/dKp
               dc= 1.0_cp/(ev(3)*ev(3) + ev(3) -k0*ev(4)) - ((1.0_cp)+ev(3)+ &
                   k0*ev(4))*(2.0_cp*ev(3) +1.0_cp)/(ev(3)*ev(3) + ev(3) -k0*ev(4))**2.0_cp

               vterm = (vv0/abc(1) + 1.0_cp - 1.0_cp/abc(1))**(-1.0_cp/abc(3))
               term1=-1.0_cp/abc(2)/abc(2) * db * vterm
               term2= 1.0_cp/abc(2)*vterm/abc(3)/abc(3)*dc*log((vv0-1.0_cp)/abc(1) + 1.0_cp)
               term2=term2+((vv0-1.0_cp)/abc(1) + 1.0_cp)**(-1.0_cp/abc(3)-1.0_cp)/abc(2)/abc(3)/abc(1)/abc(1)*(vv0-1.0_cp)*da
               term3=db/abc(2)/abc(2)
               td(3)=term1+term2+term3
            end if

            !> dp/dKpp0
            if (abs(ev(4)) > tiny(0.) .and. eospar%iorder == 4)then          ! if Kpp0 = 0 then the derivative is undefined
               da=-1.0_cp*k0*abc(1)/(1.0_cp+ev(3)+k0*ev(4))
               db=-1.0_cp/(1.0+ev(3))
               dc=(k0*(ev(3)+1.0_cp)**2.0_cp)/(ev(3)*ev(3)+ev(3)-k0*ev(4))**2.0_cp

               term1= -1.0_cp/abc(2)/abc(2)*db*(vterm-1.0_cp)
               term2= 1.0_cp/abc(2)*vterm/abc(3)/abc(3)*dc*log((vv0-1.0_cp)/abc(1) + 1.0_cp)
               term3= ((vv0-1.0_cp)/abc(1) + 1.0_cp)**(-1.0_cp/abc(3)-1.0_cp)/abc(2)/abc(3)/abc(1)/abc(1)*(vv0-1.0_cp)*da
               td(4)=term1+term2+term3
            end if
      end select

      if (eospar%ITherm > 0) then
         !> First do the parameters common to all thermal
         delt=t-eospar%tref
         delt2=t*t-eospar%tref*eospar%tref
         deltinv=0.0_cp
         if(eospar%tref .gt. 1.0)deltinv=1.0_cp/t-1.0_cp/eospar%tref

         !> Adjust dp/dv(0,t) to dp/dv(0,0):
         dpdvt=td(1)
         td(1)=td(1)*vt0/ev(1)

         !> dp/dk(Tref,0) = dp/dk(T,0)

         !> d(k)/d(t) ...(the temperature dependence of k)
         td(5)=td(2)*delt

         !> Now do the specific derivatives for each equation, where possible. Program uses numeric derivs anyway
         select case(eospar%itherm)     !  dp/dalpha=dp/dV(0T) . dV(0T)/dalpha
            case(1)                ! Berman
               !> dp/dalph0:
               td(10)=dpdvt*delt*ev(1)

               !> dp/dalpha1
               td(11)=dpdvt*0.5_cp*delt*delt*ev(1)

            case(2)                ! Fei
               !> dp/dalph0:
               td(10)=dpdvt*delt*vt0

               !> dp/dalpha1
               td(11)=dpdvt*0.5_cp*delt2*vt0

               !> dp/dalpha2
               td(12)= -1.0_cp*dpdvt*deltinv*vt0
         end select
      end if

      !> fix derivatives for linear eos
      if (eospar%linear)then
         td(1) =td(1) *(3.0_cp*eospar%params(1)**2.0_cp)
         td(2:5)=td(2:5)/3.0_cp

         select case(eospar%itherm)                 ! thermal expansion terms
            case(1:2)                              ! Berman, Fei, alpha terms
               td(10:12)=td(10:12)*3.0_cp

            case(4)                                ! Holland-Powell alpha
               td(10)=td(10)*3.0_cp
         end select                     !other thermal equations have nothing to convert
      end if

      return
   End Subroutine Deriv_Partial_P_Analytic

   !!--++
   !!--++ SUBROUTINE Deriv_Partial_P_Numeric(V,T,Eospar,td)
   !!--++    real(kind=cp),                      intent(in)    :: V       ! Volume
   !!--++    real(kind=cp),                      intent(in)    :: T       ! Temperature
   !!--++    type(Eos_Type),                     intent(in)    :: Eospar  ! Eos Parameter
   !!--++    real(kind=cp), dimension(n_eospar), intent(out)   :: td      ! derivatives dP/d(param)
   !!--++
   !!--++ PRIVATE
   !!--++ Calculates the partial derivatives of P with respect to the EoS
   !!--++ at a given v,t point, and returns  them in array td
   !!--++
   !!--.. 27-Feb-2013: RJA edits to use IREF as in original
   !!--++
   !!--++ Date: 10/09/2013
   !!
   Subroutine Deriv_Partial_P_Numeric(V,T,Eospar,td)
      !---- Arguments ----!
      real(kind=cp),                      intent(in) :: V       ! Volume
      real(kind=cp),                      intent(in) :: T       ! Temperature
      type(Eos_Type),                     intent(in) :: Eospar  ! Eos Parameter
      real(kind=cp), dimension(n_eospar), intent(out):: td      ! derivatives dP/d(param)

      !---- Local Variables ----!
      type(Eos_Type)                 :: Eost                 ! Eos Parameter local copy
      real(kind=cp), dimension(-2:2) :: p                    ! array for calc p values
      real(kind=cp)                  :: delfactor,del,d_prev,delmin
      integer                        :: i,j,icycle

      !> Initialise
      td=0.0_cp
      eost=eospar
      call Set_Eos_Use(eost)

      !> Set the inital shift factor (fraction of parameter value)
      delfactor=0.01_cp

      do i=1,n_eospar
          if (eospar%iuse(i) == 1) then             ! refineable parameters only
             del=delfactor*eospar%params(i)         ! the initial shift estimate
             delmin=delfactor/eospar%factor(i)      ! scale required min shift by print factors
             if (abs(del) < delmin) del=delmin      ! trap param values that are zero
             icycle=0                               ! iteration count
             d_prev=0.0_cp

             iter:do                                ! top of loop over iterations
                do j=-2,2,1
                   eost=eospar                                             !reset eos params
                   eost%params(i)=eost%params(i)+float(j)*del              ! apply shift to a parameter
                   p(j)=get_pressure(v,t,eost)                             ! calc resulting P
                end do
                td(i)=(p(-2)+8.0_cp*(p(1)-p(-1))-p(2))/(12.0_cp*del)       ! derivative to second order approximation

                !write(6,'(''  Param # '',i2,'' Cycle '',i2,'': deriv = '',ES14.6,'' delp = '',f5.3,'' del= '',f9.4)')i,icycle,td(i),p(2)-p(-2),del

                if (abs(td(i)) < 1.0E-8) exit iter                         ! zero deriv
                if (icycle > 0 .and.&
                            abs(d_prev-td(i))/td(i) < 1.0E-4) exit iter    ! deriv converged to 1 part in 10^4

                d_prev=td(i)                ! store last deriv value
                del=2.0_cp*del              ! increase the shift
                icycle=icycle+1
                if (icycle > 5)exit iter    ! Do not allow 2*shift to exceed 64% of param value
             end do iter
          end if
      end do

      !> no need to fix derivatives for linear eos by this method

      return
   End Subroutine Deriv_Partial_P_Numeric

   !!----
   !!---- SUBROUTINE EoS_Cal(P,T,EoSpar,Parvals)
   !!----    real(kind=cp),               intent(in)  :: P       ! Pressure
   !!----    real(kind=cp),               intent(in)  :: T       ! Temperature
   !!----    type(Eos_Type),              intent(in)  :: EoSPar  ! Eos Parameter
   !!----    real(kind=cp),  dimension(:),intent(out) :: Parvals ! Output parameter values
   !!----
   !!---- Returns elastic properties (not the parameter values) at this P,T for EoS
   !!----
   !!----  Date: 25-Feb-2013 RJA
   !!
   Subroutine EoS_Cal(P,T,EoSpar,Parvals)
      !---- Arguments ----!
      real(kind=cp),                intent(in)  :: P       ! Pressure
      real(kind=cp),                intent(in)  :: T       ! Temperature
      type(Eos_Type),               intent(in)  :: EoSPar  ! Eos Parameter
      real(kind=cp),  dimension(:), intent(out) :: Parvals ! Output parameter values


      !> Init
      parvals=0.0_cp

      parvals(1)=get_volume(p,t,eospar)
      parvals(2)=k_cal(parvals(1),t,eospar)
      parvals(3)=kp_cal(parvals(1),t,eospar)
      parvals(4)=kpp_cal(parvals(1),t,eospar)
      parvals(5)=dKdT_cal(p,t,eospar)           ! dK/dT at this P,T
      parvals(6)=alpha_cal(p,t,eospar)          ! 1/V.dV/dT at this T

      call physical_check(p,t,eospar)           ! produce warnings based on P,T



      return
   End Subroutine EoS_Cal
   !!----
   !!---- SUBROUTINE EoS_Cal_Esd(P,T,EoSpar,Esd)
   !!----    real(kind=cp),               intent(in)  :: P       ! Pressure
   !!----    real(kind=cp),               intent(in)  :: T       ! Temperature
   !!----    type(Eos_Type),              intent(in)  :: EoSPar  ! Eos Parameter
   !!----    real(kind=cp),  dimension(:),intent(out) :: Esd     ! Output esd values
   !!----
   !!---- Returns esd's of the elastic properties (not the parameter values) at this P and T for EoS
   !!----
   !!----  Date: 18/09/2013 RJA
   !!
   Subroutine EoS_Cal_Esd(P,T,EoSpar,Esd)
      !---- Arguments ----!
      real(kind=cp),                intent(in)  :: P       ! Pressure
      real(kind=cp),                intent(in)  :: T       ! Temperature
      type(Eos_Type),               intent(in)  :: EoSPar  ! Eos Parameter
      real(kind=cp),  dimension(:), intent(out) :: Esd     ! Output esd values

      !---- Local Variables ----!
      real(kind=cp),dimension(n_eospar) :: esdfull

      !> Init
      esd=0.0_cp
      esdfull=0.0_cp

      !> calculate parameter esd's transformed to this P,T
      call transform_Esd(P,T,EoSpar,esdfull)
      esd(1:5)=esdfull(1:5)
      esd(6)=esdfull(10)

      !> make adjustment for only using alpha0
      select case(EoSpar%itherm)
          case(5)               ! Salje
            esd(6)=esd(6)/3.0_cp            ! alpha is 1/3 of p1
      end select



      return
   End Subroutine EoS_Cal_Esd

   !!--++
   !!--++ SUBROUTINE EoS_to_Vec(EosPar, Vec)
   !!--++    type(EoS_Type),              intent(in)  :: Eospar
   !!--++    real(kind=cp),  dimension(:),intent(out) :: Vec
   !!--++
   !!--++ PRIVATE
   !!--++ Copy parameters from EoS type to a vector
   !!--++
   !!--++ Date: 28/02/2013
   !!
   Subroutine EoS_to_Vec(EosPar, Vec)
      !---- Arguments ----!
      type(EoS_Type),              intent(in)  :: Eospar
      real(kind=cp), dimension(:), intent(out) :: Vec

      !> Init
      vec=0.0_cp

      if (.not. eospar%linear) then
         vec=eospar%params
      else
         vec(1)=eospar%params(1)**3.0_cp
         vec(2:5)=eospar%params(2:5)/3.0_cp
         select case(eospar%itherm)                 ! thermal expansion terms
             case(1:3)
                vec(10:12)=eospar%params(10:12)*3.0_cp

             case(4,5,6)
                vec(10)=eospar%params(10)*3.0_cp
                vec(11)=eospar%params(11)
         end select
      end if

      return
   End Subroutine EoS_to_Vec

   !!----
   !!---- SUBROUTINE FfCal_Dat(V,V0,P,Eospar,F,S)
   !!----    real(kind=cp),            intent(in)  :: V       ! Volume
   !!----    real(kind=cp),            intent(in)  :: V0      ! Volume at zero pressure
   !!----    real(kind=cp),            intent(in)  :: P       ! Pressure
   !!----    type(EoS_Type),           intent(in)  :: Eospar  ! Eos Parameter: only imodel and linear used
   !!----    real(kind=cp),  optional, intent(out) :: F       ! Normalised pressure
   !!----    real(kind=cp),  optional, intent(out) :: S       ! Strain
   !!----
   !!---- Return Normalized pressure and/or Strain at V and P
   !!----
   !!---- Date: 10/09/2013
   !!
   Subroutine FfCal_Dat(V,V0,P,Eospar,F,S)
      !---- Arguments ----!
      real(kind=cp),           intent(in)  :: V       ! Volume
      real(kind=cp),           intent(in)  :: V0      ! Volume at zero pressure
      real(kind=cp),           intent(in)  :: P       ! Pressure
      type(Eos_Type),          intent(in)  :: EoSPar  ! Eos Parameter: only imodel and linear used
      real(kind=cp), optional, intent(out) :: F       ! Normalised pressure
      real(kind=cp), optional, intent(out) :: S       ! Strain

      !---- Local Variables ----!
      real(kind=cp) :: vv0, sc, fc

      !> Check
      if (present(f)) f=0.0_cp
      if (present(s)) s=0.0_cp

      if (abs(v0) <= 0.00001) then
         err_eos=.true.
         err_eos_mess='V0 is zero - No strain can be calculated!'
         return
      end if

      !> Calculate VV0 or a/a0
      vv0=v/v0
      sc=Strain(VV0,EosPar)                      ! returns volume like strain for linear
      fc=NormPressure_P(sc,P,EoSPar%imodel)

      if (present(f)) f=fc
      if (present(s)) s=sc

      return
   End Subroutine FfCal_Dat

   !!----
   !!---- SUBROUTINE FfCal_Dat_Esd(V,SigV,V0,SigV0,P,SigP,EosPar,F,SigF,S,SigS)
   !!----    real(kind=cp),  intent(in)  :: V       ! Volume
   !!----    real(kind=cp),  intent(in)  :: SigV    ! Sigma (V)
   !!----    real(kind=cp),  intent(in)  :: V0      ! Vo at this temperature
   !!----    real(kind=cp),  intent(in)  :: SigV0   ! Sigma (Vo)
   !!----    real(kind=cp),  intent(in)  :: P       ! Pressure
   !!----    real(kind=cp),  intent(in)  :: SigP    ! Sigma (P)
   !!----    type(Eos_Type), intent(in)  :: EoSPar  ! Eos Parameter
   !!----    real(kind=cp),  intent(out) :: F       ! Normalised pressure
   !!----    real(kind=cp),  intent(out) :: SigF    ! Sigma (F)
   !!----    real(kind=cp),  intent(out) :: S       ! Strain
   !!----    real(kind=cp),  intent(out) :: SigS    ! Sigma (S)
   !!----
   !!---- Calculate the Normalised Pressure (F) and Strain (S) value at
   !!---- Volume (V) and Pressure (P) and also their ESD
   !!----
   !!--.. IMPORTANT: Eosparameters are not used in this calculation!
   !!--.. The only element of esopar that is used is the eos model type
   !!--.. and the linear flag when the strain function is called
   !!----
   !!---- Date: 12/03/2013
   !!
   Subroutine FfCal_Dat_Esd(V,SigV,V0,SigV0,P,SigP,EosPar,F,SigF,S,SigS)
      !---- Arguments ----!
      real(kind=cp),  intent(in)  :: V       ! Volume
      real(kind=cp),  intent(in)  :: SigV    ! Sigma (V)
      real(kind=cp),  intent(in)  :: V0      ! Vo
      real(kind=cp),  intent(in)  :: SigV0   ! Sigma (Vo)
      real(kind=cp),  intent(in)  :: P       ! Pressure
      real(kind=cp),  intent(in)  :: SigP    ! Sigma (P)
      type(Eos_Type), intent(in)  :: EoSPar  ! Eos Parameter
      real(kind=cp),  intent(out) :: F       ! Normalised pressure
      real(kind=cp),  intent(out) :: SigF    ! Sigma (F)
      real(kind=cp),  intent(out) :: S       ! Strain
      real(kind=cp),  intent(out) :: SigS    ! Sigma (S)

      !---- Local Variables ----!
      real(kind=cp) :: vv0,vv0_esd, sigmap,e

      !> Init
      f=0.0_cp
      sigf=0.0_cp
      s=0.0_cp
      sigs=0.0_cp

      if (abs(v0) <= 0.00001) then
         err_eos=.true.
         err_eos_mess='V0 is zero - No strain can be calculated!'
         return
      end if

      !> Note that V0 in call is the V0 at this temperature
      vv0=v/v0
      vv0_esd=vv0*sqrt( (sigv/abs(v))**2.0_cp + (sigv0/abs(v0))**2.0_cp )

      !> Strain
      s=Strain(vv0,EosPar)      ! input a/a0 or v/v0 to strain. It always returns volume strain

      !> Normalized Pressure
      f=NormPressure_P(s,p,EoSPar%imodel)

      !> Check Pressure values
      if (abs(p) <= 0.0001) then
         e=0.0_cp
      else
         e=sigp/p
      end if

      !> If the finite strain is zero, the esd of the normalised pressure is not defined
      if (s < tiny(0.0)) return      ! changed 23/04/2013 was ">"

      !> ESD Calculations: all done in volume terms
      if (eospar%linear)then
         vv0_esd=3.0*vv0**2.0_cp*vv0_esd
         vv0=vv0**3.0_cp
      end if

      select case (eospar%imodel)
         case (1,5) ! Murnaghan, Tait
            ! Nothing to do

         case (2) ! Birch-Murnaghan
            sigmap=((7.0_cp * vv0**(-2.0_cp/3.0_cp) -5.0_cp)*vv0_esd) / &           ! this is sigma_prime in texts
                   (3.0_cp*(1.0_cp - vv0**(-2.0_cp/3.0_cp))*vv0)

            sigs=vv0**(-5.0_cp/3.0_cp) * vv0_esd / 3.0_cp
            sigf=f*sqrt(e**2.0_cp + sigmap**2.0_cp)

         case (3) ! Vinet:
            sigs=vv0**(-2.0_cp/3.0_cp) * vv0_esd/3.0_cp
            sigmap= (s*s-1.0_cp)/s/(1.0_cp-s)**2.0_cp
            sigmap=sigmap*sigs

            sigf=f*sqrt(e**2.0_cp + sigmap**2.0_cp)

         case (4) ! Natural
            sigmap=(1.0_cp - 1.0_cp/log(vv0))*vv0_esd/vv0

            sigs=vv0_esd/3.0_cp/vv0
            sigf=f*sqrt(e**2.0_cp + sigmap**2.0_cp)
      end select

      return
   End Subroutine FfCal_Dat_Esd

   !!----
   !!----Subroutine FfCal_EoS(P,T,Eospar,F,S)
   !!----   real(kind=cp),           intent(in)  :: P       ! Pressure
   !!----   real(kind=cp),           intent(in)  :: T       ! Temperature
   !!----   type(Eos_Type),          intent(in)  :: E       ! Eos Parameters
   !!----   real(kind=cp),           intent(out) :: F       ! Normalised pressure
   !!----   real(kind=cp),           intent(out) :: S       ! Strain
   !!----
   !!----  Returns normalised pressure and Strain at this P,T,
   !!----  calculated from the eos parameters
   !!----
   !!---- Date: 05/09/2013
   !!
   Subroutine FfCal_EoS(P,T,Eospar,F,S)
      !---- Arguments ----!
      real(kind=cp),           intent(in)     :: P       ! Pressure
      real(kind=cp),           intent(in)     :: T       ! Temperature
      type(Eos_Type),          intent(in)     :: Eospar  ! Eos Parameter
      real(kind=cp),           intent(out)    :: F       ! Normalised pressure
      real(kind=cp),           intent(out)    :: S       ! Strain

      !---- Local Variables ----!
      real(kind=cp) :: v

      !> Init
      f=0.0_cp
      s=0.0_cp

      v=get_volume(p,t,eospar)        ! V at this P,T from eos parameters
      s=strain_eos(v,t,eospar)         ! finite strain at this P,T
      f=normpressure_eos(s,t,eospar)   ! Norm pressure

      return
   End Subroutine FfCal_EoS

   !!--++
   !!--++ Subroutine Get_Tait(Eospar,T,Vec)
   !!--++    type(Eos_Type),              intent(in)  :: Eospar     ! Eos Parameters
   !!--++    real(kind=cp),               intent(in)  :: T          ! Temperature
   !!--++    real(kind=cp), dimension(3), intent(out) :: Vec        ! Vector (a,b,c) of Tait parameters
   !!--++
   !!--++ Returns a,b,c Tait parameters in a vector.
   !!--++
   !!--++ Date: 11/07/2013 RJA
   !!
   Subroutine Get_Tait(Eospar,T,Vec)
     !---- Arguments ----!
     type(Eos_Type),            intent(in)  :: Eospar
     real(kind=cp),             intent(in)  :: T
     real(kind=cp),dimension(3),intent(out) :: Vec

     !---- Local Variables ----!
     real(kind=cp) :: k0,kp,kpp

     !> Init
     Vec=0.0_cp
     if (eospar%imodel /= 5)return

     select case(eospar%itherm)
         case default               ! includes no thermal model, also pthermal which requires params at Tref
             k0 =eospar%params(2)
             kp =eospar%params(3)
         case(1:5)                  ! normal thermal expansion models with dK/dT
             k0 =Get_K0_T(T,eospar)
             kp =eospar%params(3)

     end select

     if(eospar%iorder < 4)then
         kpp= -1.0_cp*kp/k0        ! implied value for Kpp except for 4th order
     else
         kpp=eospar%params(4)
     endif


     if (eospar%linear)then
        k0  =k0/3.0_cp
        kp  =kp/3.0_cp
        kpp =kpp/3.0_cp
     end if


     Vec(1)=(1.0_cp + kp)/(1.0_cp+kp+k0*kpp)

     Vec(2)= kp/k0
     if (abs(1_cp+kp) > tiny(0.0) ) Vec(2)=Vec(2)-kpp/(1.0_cp+kp)

     Vec(3)= (1.0_cp+kp+k0*kpp)/(kp*kp+kp-k0*kpp)

     return
   End Subroutine Get_Tait

   !!----
   !!---- SUBROUTINE Init_EoS_Data_Type(E)
   !!----    type (EoS_Data_Type),  intent(in out) :: E
   !!----
   !!---- Initialize EoS_Data_Type
   !!----
   !!---- Date: 16/02/2013
   !!
   Subroutine Init_EoS_Data_Type(E)
      !---- Arguments ----!
      type (EoS_Data_Type), intent(in out)   :: E

      E%IUse = 0
      E%IGrp = 0
      E%T    = 0.0
      E%P    = 0.0
      E%V    = 0.0
      E%cell = 0.0
      E%ang  = 0.0

      E%SigT = 0.0
      E%SigP = 0.0
      E%SigV = 0.0
      E%sigc = 0.0
      E%siga = 0.0

      return
   End Subroutine Init_EoS_Data_Type

   !!----
   !!---- SUBROUTINE Init_EoS_Thermal(Eospar)
   !!----    type (EoS_Type), intent(in out) :: Eospar
   !!----
   !!---- Initialize the EoS Type for Thermal case
   !!----
   !!---- Date: 10/09/2013
   !!
   Subroutine Init_EoS_Thermal(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      select case(eospar%itherm)
         case(1)
            eospar%tmodel     = 'Berman (1988)'
            eospar%parname(10:11) = (/'alph0','alph1'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Linear term thermal expansion x10^8 K^-2'
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp
            eospar%TRef    =298.0_cp       ! Simple thermal expansion,
            eospar%TRef_fixed = .false.    ! alpha terms can be safely 0.

         case(2)
            eospar%tmodel     = 'Fei (1995)'
            eospar%parname(10:12) = (/'alph0','alph1','alph2'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Linear term thermal expansion x10^8 K^-2'
            eospar%comment(12) = '1/T^2 term thermal expansion, K'
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp
            eospar%factor(12)  = 1.0_cp
            eospar%TRef    =298.0_cp       ! Simple thermal expansion,
            eospar%TRef_fixed = .false.

         case(3)
            eospar%tmodel     = 'Modified HP1998'
            eospar%parname(10:11) = (/'alph0','alph1'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Sqrt term of thermal expansion x10^4 K^-1/2'
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E4_cp
            eospar%TRef        =298.0_cp       ! Simple thermal expansion,
            eospar%TRef_fixed  = .false.

         case(4)
            eospar%tmodel     = 'Kroll'
            eospar%parname(10:11) = (/'alph0','Th_E '/)
            eospar%comment(10) = 'Constant of thermal expansion at Tref x10^5 K^-1'
            eospar%comment(11) = 'Einstein temperature in K'
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef    =298.0_cp     ! Holland and Powell thermal expansion
            eospar%TRef_fixed = .false.
            eospar%params(11) = 298.0_cp   ! Einstein temperature

         case(5)
            eospar%tmodel     = 'Salje low-T'
            eospar%parname(10:11) = (/'p1   ','T_sat'/)
            eospar%comment(10) = 'Approx 3x highT thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Saturation temperature in K'
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef    =0.0_cp     ! Salje thermal expansion
            eospar%TRef_fixed = .true.
            eospar%params(11) = 298.0_cp   ! Saturation temperature

         case(6)
            eospar%tmodel     = 'HP Thermal Pressure'
            eospar%parname(10:11) = (/'alph0','Th_E '/)
            eospar%comment(10) = 'Constant of thermal expansion at Tref x10^5 K^-1'
            eospar%comment(11) = 'Einstein temperature in K'
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef    =298.0_cp
            eospar%TRef_fixed = .false.
            eospar%params(11) = 298.0_cp   ! Einstein temperature

      end select

      call Set_Eos_Use(Eospar)                 ! update the use flags

      return
   End Subroutine Init_EoS_Thermal

   !!----
   !!---- SUBROUTINE Init_EoS_Type(Eospar,CLin,IThermal)
   !!----    type (EoS_Type),  intent(in out) :: Eospar
   !!----
   !!---- Initialize EoS_Type setting all parameters to sensible values
   !!----
   !!---- Date: 10/09/2013
   !!
   Subroutine Init_EoS_Type(Eospar,CLin,IThermal)
      !---- Arguments ----!
      type (EoS_Type),            intent(in out) :: Eospar    ! EoS Type
      character(len=*), optional, intent(in)     :: CLin      ! Character variable to indicate linear EoS or not
      integer,          optional, intent(in)     :: IThermal  ! integer to indicate ithermal type

      !> test for optional argument for linear
      Eospar%Linear  =.false.
      if (present(clin))then
         if (index(U_case(clin(1:3)),'LIN') > 0) Eospar%linear=.true.
      end if

      Eospar%IModel  =0
      Eospar%IOrder  =3
      Eospar%Title   =' '

      eospar%ParName=' '
      eospar%comment=' '
      call Set_Eos_Names(Eospar)         ! also sets the print factors for the pressure part

      Eospar%PRef    =0.0_cp

      Eospar%Iuse    =0
      Eospar%Iuse(1:4)=1                 ! Vo, Ko, Kpp

      Eospar%params  =0.0_cp
      Eospar%esd     =0.0_cp

      Eospar%params(1)=1.0_cp
      Eospar%params(2)=10.0_cp
      Eospar%params(3)=4.0_cp

      Eospar%X       =0.0_cp

      Eospar%WChi2   =0.0_cp
      Eospar%DelPMax =0.0_cp
      Eospar%IWt     =0

      Eospar%IRef    =0
      Eospar%factor  =1.0_cp
      Eospar%LastShift=0.0_cp
      Eospar%VCV     =0.0_cp

      !> test for optional argument for thermal
      Eospar%ITherm  =0
      Eospar%Tref = 298.0_cp         !sensible default
      if (present(ithermal) .and. ithermal /= 0)then
         Eospar%Itherm = ithermal
         call Init_Eos_Thermal(Eospar)              ! set up default values, names for specific thermal eqn
      end if

      return
   End Subroutine Init_EoS_Type

   !!----
   !!---- SUBROUTINE Init_Err_Eos()
   !!----
   !!---- Initialize Flags of Errors in this module
   !!----
   !!---- Date: 16/02/2013
   !!
   Subroutine Init_Err_EoS()

      Err_EoS=.false.
      Err_EoS_Mess=" "

      return
   End Subroutine Init_Err_EoS

   !!--++
   !!--++ SUBROUTINE Physical_Check(P,T,E)
   !!--++    real(kind=cp),    intent(in) :: p      ! Pressure
   !!--++    real(kind=cp),    intent(in) :: t      ! Temperature
   !!--++    type(Eos_Type),   intent(in) :: E      ! EoS type
   !!--++
   !!--++ PRIVATE
   !!--++ Check if the parameters have physical sense
   !!--++
   !!--++ Date: 16/02/2013
   !!
   Subroutine Physical_Check(P,T,E)
      !---- Arguments ----!
      real(kind=cp),    intent(in) :: p
      real(kind=cp),    intent(in) :: t
      type(Eos_Type),   intent(in) :: E

      !---- Local variables ----!
      character(len=20)   :: car
      real(kind=cp)       :: tlimit,v

      !> Init
      call Init_Err_EoS()

      !> Basic checks
      v=get_volume(p,t,e)
      if (v < tiny(0.0) ) then
         err_eos=.true.
         write(unit=car, fmt='(2f10.1)') p, t
         car=adjustl(car)
         err_eos_mess='Volume calculated as zero or negative at P,T = '//trim(car)
      end if

      if (e%imodel > 0)then
         if (.not. e%linear .and. k_cal(v,t,e) < tiny(0._cp)) then
            err_eos=.true.
            write(unit=car, fmt='(2f10.1)') p, t
            car=adjustl(car)
            err_eos_mess='Bulk modulus calculated as zero or negative at P,T = '//trim(car)
         end if
      end if

      !> Check validity of thermal model
      select case(e%itherm)
         case(2)                ! Fei:
            if (e%params(12) > tiny(0.0_cp)) then  ! non-physical V and divergent alpha at low T when alpha2 .ne. 0
               tlimit=(2.0_cp*e%params(12)/e%params(11))**(1.0_cp/3.0_cp)
               if (t < tlimit) then
                  err_eos=.true.
                  write(unit=car,fmt='(f5.1)')tlimit
                  car=adjustl(car)
                  err_eos_mess='Fei equation yields non-physical behaviour below T = '//trim(car)//'K'
               end if
            else if(e%params(12) < tiny(0.0_cp)) then  ! alpha2 < 0
               tlimit=sqrt(-1.0_cp*e%params(12)/e%params(10))
               if (t < tlimit)then
                  err_eos=.true.
                  write(unit=car,fmt='(f5.1)')tlimit
                  car=adjustl(car)
                  err_eos_mess='Fei equation yields non-physical behaviour below T = '//trim(car)//'K'
               end if
            end if

         case(3)               ! HP 1998: trap non-physical behaviour at low T
            tlimit=((10.0_cp*e%params(10)+e%params(11))/e%params(10))**2.0_cp
            if (t < tlimit) then
               err_eos=.true.
               write(unit=car,fmt='(f5.1)')tlimit
               car=adjustl(car)
               err_eos_mess='HP1998 equation yields non-physical behaviour below T = '//trim(car)//'K'
            end if

         case(1,4,5,6)
            if (t < 0.0_cp) then
               err_eos=.true.
               err_eos_mess='T  less than zero K'
            end if
      end select

      return
   End Subroutine Physical_Check

   !!----
   !!---- SUBROUTINE Read_EoS_DataFile(fname,dat)
   !!----    character(len=*),          intent(in)  :: fname
   !!----    type (eos_data_list_type), intent(out) :: dat  ! data structure
   !!----
   !!---- General routine to read data for Eos
   !!----
   !!---- Update: 18/04/2013
   !!
   Subroutine Read_EoS_DataFile(fname,dat)
      !---- Arguments ----!
      character(len=*),          intent(in)  :: fname
      type (eos_data_list_type), intent(out) :: dat  ! data structure

      !---- Local Variables ----!
      character(len=132), dimension(:), allocatable :: flines
      character(len=132)                            :: line
      character(len=5)                              :: car
      character(len=1)                              :: Ts
      character(len=30), dimension(ncol_data_max)   :: dire

      integer                                       :: nldata,ndat
      integer                                       :: i,j,kk,nl,nk,nlines,m
      integer                                       :: iv, inum
      integer, dimension(ncol_data_max)             :: ivet,iorden

      real(kind=cp), dimension(ncol_data_max)       :: vet,vetsd
      real(kind=cp), dimension(ncol_data_max)       :: rvet

      !---- Data Files ----!
      character(len=512)                            :: filedat=' '
      integer                                       :: NCDat         ! Total columns for data
      integer, dimension(ncol_data_max)             :: IC_Dat        ! Which values are input - local copy

      !> Eos init
      call Init_err_Eos()

      !> Init
      dire=' '                       ! Clear, otherwise a second read may have items in dire from 1st read
      filedat=trim(fname)

      !> Number of lines
      call number_lines (trim(filedat),nlines)
      if (nlines <= 0) then
         err_eos=.true.
         Err_EoS_Mess="Impossible to read the file "
         return
      end if

      !> Read lines
      if (allocated(flines)) deallocate(flines)
      allocate(flines(nlines))
      flines=' '
      call reading_lines(trim(filedat),nlines,flines)

      !> Title
      i=1
      j=nlines
      dat%title=''
      call Read_Key_Str(flines, i, j, 'Title', line,'#')
      if (len_trim(line) > 0) dat%title=trim(line)

      !> System
      i=1
      j=nlines
      dat%system=''
      call Read_Key_Str(flines, i, j, 'System', line,'#')
      if (len_trim(line) > 0) dat%system=trim(line)

      !> TScale
      Ts='K'
      i=1
      j=nlines
      call Read_Key_Str(flines, i, j, 'TScale', line,'#')
      if (len_trim(line) > 0) Ts=U_case(trim(adjustl(line)))

      !> Format
      i=1
      j=nlines
      call Read_Key_Str(flines, i, j, 'Format', line,'#')

      !> for compatability with eosfit v5, remove any commas
      !> from format directive before processing
      call Get_Separator_pos(line,',',ivet,iv)
      do m=1,iv
         line(ivet(m):ivet(m))=' '
      end do

      call getword(line,dire,iv)
      if (iv <= 1) then
         err_eos=.true.
         Err_EoS_Mess="The options for directive Format are wrong!"
         return
      end if

      nldata=i+1      ! Line where begin the data

      !> Is a number the first value on format line?
      !if 'yes' then we need to read dire(2:iv) and ncdat=iv-1
      !if 'no' then we need to read dire(1:iv) and ncdat=iv

      inum=0
      call getnum(dire(1),vet,ivet,m)
      if (m == 1) then
         inum=1
         ncdat=iv-1
      else
         ncdat=iv
      end if
      ic_dat=0
      ic_dat(1)=1     ! CheckBox (active data)
      ic_dat(22)=1    ! Group 1
      iorden=0

      do i=1,ncdat          ! this must loop from 1 to ncdat always so that iorden is filled correctly
         car=U_case(adjustl(dire(i+inum)))

         !> allow 'press, pressure, vol, volume, temp...' for full compatability with old files
         if (index(car,'PRESS') /= 0)  car='P'
         if (index(car,'VOL')   /= 0)  car='V'
         if (index(car,'TEM')   /= 0)  car='T'

         select case (trim(car))
            case ('T')
               ic_dat(4)=1
               iorden(i)=4
            case ('SIGT')
               ic_dat(5)=1
               iorden(i)=5
            case ('P')
               ic_dat(6)=1
               iorden(i)=6
            case ('SIGP')
               ic_dat(7)=1
               iorden(i)=7
            case ('V')
               ic_dat(8)=1
               iorden(i)=8
            case ('SIGV')
               ic_dat(9)=1
               iorden(i)=9
            case ('A')
               ic_dat(10)=1
               iorden(i)=10
            case ('SIGA')
               ic_dat(11)=1
               iorden(i)=11
            case ('B')
               ic_dat(12)=1
               iorden(i)=12
            case ('SIGB')
               ic_dat(13)=1
               iorden(i)=13
            case ('C')
               ic_dat(14)=1
               iorden(i)=14
            case ('SIGC')
               ic_dat(15)=1
               iorden(i)=15
            case ('ALPHA')
               ic_dat(16)=1
               iorden(i)=16
            case ('SIGAL')
               ic_dat(17)=1
               iorden(i)=17
            case ('BETA')
               ic_dat(18)=1
               iorden(i)=18
            case ('SIGBE')
               ic_dat(19)=1
               iorden(i)=19
            case ('GAMMA')
               ic_dat(20)=1
               iorden(i)=20
            case ('SIGGA')
               ic_dat(21)=1
               iorden(i)=21
         end select
      end do

      !> Now calculate ncdat per line
      !if (n_ini ==0) ncdat=ncdat+1

      call getnum(dire(1),vet,ivet,iv)
      if (iv <=0) then
         nl=1
      else
         nl=ivet(1)
      end if

      !> Estimating Number of points
      ndat=0
      do j=nldata,nlines
         line=adjustl(flines(j))
         if (len_trim(line) <= 0) cycle
         if (line(1:1) == '#') cycle
         ndat=ndat+1
      end do
      ndat=ndat/nl
      if (ndat <=0) then
         err_eos=.true.
         Err_EoS_Mess='Number of data points estimated in the file was zero!'
         return
      end if

      !> Allocating data
      call allocate_eos_data_list(ndat,dat)

      !> reading file
      ndat=0
      nk=0
      rvet=0.0
      do i=nldata,nlines
         line=adjustl(flines(i))
         if (len_trim(line) <= 0 .or. line(1:1) =='#') cycle

         !---- Data Points ----!
         !> Replace ',' by blank character
         call Get_Separator_pos(line,',',ivet,iv)
         do j=1,iv
            line(ivet(j):ivet(j))=' '
         end do

         !> Replace tab by blank character
         call Get_Separator_pos(line,char(9),ivet,iv)
         do j=1,iv
            line(ivet(j):ivet(j))=' '
         end do

         call getnum_std(line,vet,vetsd,iv)
         if (iv == 0) then
            err_eos=.true.
            Err_EoS_Mess='Error in the number of values in a data line!'
            dat%n=ndat
            return
         end if

         do j=1,iv
            kk=nk+j
            if (iorden(kk) < 1 .or. iorden(kk) > ncol_data_max)then
               err_eos=.true.
               write(unit=car,fmt='(i5)') i
               car=adjustl(car)
               Err_EoS_Mess='Error reading data at line '//trim(car)//'  Did you put sig in format line, but have sig in ()?'
               dat%n=ndat
               return
            end if
            rvet(iorden(kk))=vet(j)             ! if esd's listed as separate items, they are in vet

            select case (iorden(kk))            ! if esd's listed with () they are in vetsd
               case (4) ! T
                  if (vetsd(j) > 0.0)rvet(5)=vetsd(j)

               case (6) ! P
                  if (vetsd(j) > 0.0) rvet(7)=vetsd(j)

               case (8) ! V
                  if (vetsd(j) > 0.0) rvet(9)=vetsd(j)

               case (10) ! a
                  if (vetsd(j) > 0.0) rvet(11)=vetsd(j)

               case (12) ! b
                  if (vetsd(j) > 0.0) rvet(13)=vetsd(j)

               case (14) ! c
                  if (vetsd(j) > 0.0) rvet(15)=vetsd(j)

               case (16) ! alpha
                  if (vetsd(j) > 0.0) rvet(17)=vetsd(j)

               case (18) ! beta
                  if (vetsd(j) > 0.0) rvet(19)=vetsd(j)

               case (20)
                  if (vetsd(j) > 0.0) rvet(21)=vetsd(j)
            end select
         end do
         nk=nk+iv
         if (nk < ncdat) cycle

         !> Writting values on Type
         ndat=ndat+1
         dat%eosd(ndat)%iuse=1          ! Active data
         dat%eosd(ndat)%igrp=1          ! Group 1

         !> Convert to Kelvin
         select case (Ts)
            case ('C')
               dat%eosd(ndat)%t=rvet(4) + 273.15

            case ('F')
               dat%eosd(ndat)%t=(rvet(4) + 459.67)/1.8

            case default
               dat%eosd(ndat)%t=rvet(4)
         end select

         dat%eosd(ndat)%p=rvet(6)
         dat%eosd(ndat)%v=rvet(8)
         dat%eosd(ndat)%cell=(/rvet(10),rvet(12),rvet(14)/)
         dat%eosd(ndat)%ang =(/rvet(16),rvet(18),rvet(20)/)
         dat%eosd(ndat)%sigt=rvet(5)
         dat%eosd(ndat)%sigp=rvet(7)
         dat%eosd(ndat)%sigv=rvet(9)
         dat%eosd(ndat)%sigc=(/rvet(11),rvet(13),rvet(15)/)
         dat%eosd(ndat)%siga=(/rvet(17),rvet(19),rvet(21)/)

         !> New 02/08/2013 to set flags for esd's when they were provided in ()
         do j=1,21,2
            if (rvet(j) > 0.0)ic_dat(j)=1
         end do

         nk=0
         rvet=0.0
      end do

      dat%n=ndat
      dat%ic_dat=ic_dat         ! flags for original input

      !> Default values in function of system
      call Def_Crystal_System(dat)
      if (err_eos) return         ! trap error because of inconsistent data

      !> Volume calculation if have cell parameters
      call Set_Volume_from_Cell(dat)

      return
   End Subroutine Read_EoS_DataFile

   !!----
   !!---- SUBROUTINE Set_Eos_Names(Eospar)
   !!----    type (EoS_Type), intent(in out)      :: Eospar
   !!----
   !!---- Set the character variables in eos_type data structures
   !!---- to match the flags already set
   !!----
   !!---- Date: 10/09/2013
   !!
   Subroutine Set_Eos_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !> Set the Eos name
      select case(eospar%imodel)
         case(0)
            eospar%model='None'
         case(1)
            eospar%model='Murnaghan'
         case(2)
            eospar%model='Birch-Murnaghan'
         case(3)
            eospar%model='Vinet'
         case(4)
            eospar%model='Natural Strain'
         case(5)
            eospar%model='Tait'
      end select

      !> set the comments for parameters for volume or linear eos
      !eospar%ParName=' '
      !eospar%comment=' '
      if (.not. eospar%linear) then
         eospar%ParName(1:5) =(/'V0   ','K0   ','Kp   ','Kpp  ','dK/dT'/)

         eospar%comment(1) = 'Zero pressure volume: units as volume data'
         eospar%comment(2) = 'Bulk modulus: same units as pressure data'
         eospar%comment(3) = 'dK/dP: dimensionless'
         eospar%comment(4) = 'd2K/dP2 inverse pressure units'
         eospar%comment(5) = 'dK/dT: units are P units/K'
      else
         eospar%ParName(1:5) =(/'a0   ','M0   ','Mp   ','Mpp  ','dM/dT'/)

         eospar%comment(1) = 'Zero pressure length: units as cell parameters'
         eospar%comment(2) = 'Linear modulus: same units as pressure data'
         eospar%comment(3) = 'dM/dP: dimensionless'
         eospar%comment(4) = 'd2M/dP2: inverse pressure units'
         eospar%comment(5) = 'dM/dT: : units are P units/K'
      end if

      !> Thermal models are only set in init_EoS_thermal
      return
   End Subroutine Set_Eos_Names

   !!----
   !!---- SUBROUTINE Set_EoS_Use(Eospar)
   !!----    type (EoS_Type), intent(in out) :: Eospar
   !!----
   !!---- sets the 'use' flags for Eos type based on all current settings
   !!----
   !!----    Iuse     Comments
   !!----   ----------------------------------------------------------------------------
   !!----      0      parameter not used
   !!----      1      parameter is used, settable, refineable
   !!----      2      parameter is used and/or should be reported, settable, but cannot be refined
   !!----      3      parameter is used and/or should be reported, not settable, cannot be refined (includes implied values)
   !!----
   !!---- Date: 10/09/2013
   !!
   Subroutine Set_EoS_Use(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !---- Local Variables ----!
      integer :: i

      !> Init
      eospar%iuse=0

      select case(eospar%imodel)
      case(0)
          eospar%iuse(1)=1                                      ! None eg thermal only
      case(1)
         eospar%iuse(1:3)=1                                     ! Murnaghan

      case default
        eospar%iuse(1:eospar%iorder)=1                          ! other isothermal EoS
        if(eospar%iorder < 4) eospar%iuse(eospar%iorder+1:4)=3  !implied values
      end select

      select case(eospar%itherm)
         case(1)             ! Berman
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10:11)=1 ! alpha terms

         case(2)             ! Fei
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10:12)=1 ! alpha terms

         case(3)             ! HP 1998
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10:11)=1 ! alpha terms

         case(4)             ! Holland-Powell thermal expansion, in Kroll form
             if(eospar%iuse(3) == 0)eospar%iuse(3)=2     ! require Kprime_zero but not stable in refinement if no P data (added 27/01/2014 RJA)
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10)=1    ! alpha at Tref
             eospar%iuse(11)=2    ! Einstein T should be reported but cannot be refined

         case(5)             ! Salje
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10:11)=1

         case(6)             ! Thermal pressure in H&P form (no dK/dT): requires a eos model as well
             eospar%iuse(2:3)=2   ! K, Kp refinement is not stable
             eospar%iuse(5)=0     ! No dK/dT parameter:
             eospar%iuse(10)=1    ! alpha at Tref
             eospar%iuse(11)=2    ! Einstein T should be reported but cannot be refined
      end select

      !> JGP
      do i=1,n_eospar
         if (eospar%iuse(i) /=1) eospar%iref(i)=0
      end do

      return
   End Subroutine Set_Eos_Use

   !!----
   !!---- SUBROUTINE Set_Kp_Kpp_Cond(Eospar)
   !!----    type (EoS_Type), intent(in out) :: Eospar
   !!----
   !!---- Fix Kp and Kpp values from Model and Order of EoSpar
   !!----
   !!---- Date: 16/02/2013
   !!
   Subroutine Set_Kp_Kpp_Cond(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !---- Local Variables ----!
      real(kind=cp),dimension(N_EOSPar):: ev        ! local copies of room pressure parameter values

      !> Local copy
      call EoS_to_Vec(eospar,ev) !  ev contains volume-like parameters

      select case (eospar%imodel)
         case (1) ! Murnaghan
            ev(4)=0.0_cp

         case (2) ! Birch-Murnaghan
            if (eospar%iorder == 2) ev(3)=4.0_cp
            if (eospar%iorder == 2 .or. eospar%iorder == 3) then
               if (abs(ev(2)) > 0.0) ev(4)=-1.0_cp*((ev(3)-4.0_cp)*(ev(3)-3.0_cp)+35.0_cp/9.0_cp)/ev(2)  !for order 2 and 3
            end if

         case (3) ! Vinet
            if (eospar%iorder == 2) ev(3)=1.0_cp
            if (eospar%iorder == 2 .or. eospar%iorder == 3) then
               if (abs(ev(2)) > 0.0) ev(4)=-1.0_cp*((0.5_cp*ev(3))**2+0.5*ev(3)-19.0_cp/36.0_cp)/ev(2) !for order 2 and 3
            end if

         case (4) ! Natural
            if (eospar%iorder == 2) ev(3)=2.0_cp
            if (eospar%iorder == 2 .or. eospar%iorder == 3) then
               if (abs(ev(2)) > 0.0) ev(4)=-1.0_cp*(1.0_cp + (ev(3)-2.0_cp)+(ev(3)-2.0_cp)**2.0_cp)/ev(2) !for order 2 and 3
            end if

         case(5) ! Tait with definitions of order derived from Holland and Powell (2011)
            if (eospar%iorder == 2) ev(3)=4.0_cp
            if (eospar%iorder == 2 .or. eospar%iorder == 3)then
               if (abs(ev(2)) > 0.0)ev(4)=-1.0_cp*ev(3)/ev(2)
            endif
      end select

      if (.not. eospar%linear) then
         if (eospar%iorder == 2) eospar%params(3)=ev(3)
         if (eospar%iorder == 2 .or. eospar%iorder == 3) eospar%params(4)=ev(4)
      else
         if (eospar%iorder == 2) eospar%params(3)=ev(3)*3.0_cp
         if (eospar%iorder == 2 .or. eospar%iorder == 3) eospar%params(4)=ev(4)*3.0_cp
      end if

      !> take opportunity to update the use flags
      call Set_Eos_Use(eospar)

      return
   End Subroutine Set_Kp_Kpp_Cond

   !!--++
   !!--++ SUBROUTINE Set_Volume_from_Cell(dat)
   !!--++    type (eos_data_list_type),  intent(in out)   :: dat  ! data structure
   !!--++
   !!--++ PRIVATE
   !!--++ Sets V and esd(V) from cell parameter data for all data items in dat
   !!--++ If V is present in first data item, no esd is calculated
   !!--++
   !!--++ Update: 18/04/2013
   !!
   Subroutine Set_Volume_from_Cell(dat)
      !---- Arguments ----!
      type (eos_data_list_type),  intent(in out)  :: dat  ! data structure

      !---- Local Variables ----!
      integer :: i

      !> Check
      if (dat%eosd(1)%v > 0.0) return
      if (any(dat%eosd(1)%cell < 0.5)) return

      !> Calculations
      do i=1,dat%n
         call Volume_Sigma_from_Cell(dat%eosd(i)%cell,dat%eosd(i)%ang,dat%eosd(i)%sigc, &
                                     dat%eosd(i)%siga,dat%eosd(i)%v,dat%eosd(i)%sigv)
      end do

      return
   End Subroutine Set_Volume_from_Cell

   !!--++
   !!--++ SUBROUTINE Transform_Esd(P,T,Eospar,Esd)
   !!--++    real(kind=cp),                intent(in)      :: P       ! Pressure at which to calculate parameter esds
   !!--++    real(kind=cp),                intent(in)      :: T       ! Temperature at which to calculate parameter esds
   !!--++    type (EoS_Type),              intent(in)      :: Eospar  ! The EoS parameters and vcv
   !!--++    real(kind=cp),dimension(:),   intent(out)     :: Esd
   !!--++
   !!--++ PRIVATE
   !!--++ New generalised version for V,K,Kp,Kpp,dK/dT and alpha
   !!--++
   !!--..  NOTE:
   !!--..       the value of 0.001 was chosen to balance non-linearity of
   !!--..       equations (which implies small value) against the arithmatic
   !!--..       precision (wants big value). this value gives same esd's as
   !!--..       double-precision eosfit5
   !!--++
   !!--++ Date: 02/08/2013
   !!
   Subroutine Transform_Esd(P,T,Eospar,Esd)
      !---- Arguments ----!
      real(kind=cp),   intent(in)              :: p        ! Pressure at which to calculate parameter esds
      real(kind=cp),   intent(in)              :: t        ! Temperature at which to calculate parameter esds
      type (EoS_Type), intent(in)              :: eospar   ! The EoS parameters and vcv
      real(kind=cp),dimension(:), intent(out)  :: esd      ! The esd's of Eos parameters at this P and T

      !---- Local Variables ----!
      real(kind=cp),dimension(N_EOSPar,N_EOSPar)  :: d         ! cross derivatives of parameters
      real(kind=cp),dimension(-2:2,6)             :: par
      real(kind=cp), parameter                    :: factor=0.01        ! shift factor for calculating derivatives
      real(kind=cp)                               :: shift              ! shift applied to parameter
      type (eos_type)                             :: eost               ! local copy of input eos
      integer                                     :: i,j,k

      !> initialisation
      d=0.0_cp
      esd=0.0_cp

      !> loop to calculate d(param_i at P,T)/d(param_j at Pref,Tref)
      do j=1,10
          if (eospar%iuse(j) == 1)then
              do k=-2,2,1
                  eost=eospar                                                     !reset eos params
                  if(abs(eospar%params(j)) > tiny(0.0_cp))then
                      shift=factor*eospar%params(j)
                      if(j == 10)shift=10.0*shift                   ! alpha0
                  else
                      shift=1.0_cp/eospar%factor(j)                 ! shift to a parameter with zero value
                      if(j == 5)shift=0.002                         ! dK/dT has print factor 1, but typical value -.02
                  endif
                  eost%params(j)=eospar%params(j)+float(k)*shift    ! apply shift to a parameter
                  call EoS_Cal(P,T,Eost,Par(k,1:6))                 ! calc resulting parvals
              end do
              d(1:6,j)=(par(-2,1:6)+8.0_cp*(par(1,1:6)-par(-1,1:6))-par(2,1:6))/(12.0_cp*shift) ! derivative to second order approximation
          end if
      end do

      !> d(1:6,j) contains the derivatives in order V,K0,Kp,Kpp,dK/dT,alpha0 wrt param(j) at Pref,Tref
      !  now switch these to 10 for alpha0 and ignore any other alpha coeffs
      d(10,1:10)=d(6,1:10)
      d(6,1:10)=0.0_cp

      !> Now calculate esd-squared array = vcv(i,i). The input vcv is already linear for linear eos!
      do k=1,N_EOSPar
         do i=1,N_EOSPar
            if(eospar%iref(i) == 0)cycle ! because vcv(i,j) will be 0
            do j=1,N_EOSPar
            if(eospar%iref(j) == 0)cycle ! because vcv(i,j) will be 0
               esd(k)=esd(k)+eospar%vcv(i,j)*d(k,i)*d(k,j)
            end do
         end do
      end do

      !> Final
      esd=sqrt(esd)

      return
   End Subroutine Transform_Esd

   !!--++
   !!--++ SUBROUTINE Vec_to_EoS(Vec, EosPar)
   !!--++    real(kind=cp),  dimension(:),intent(in)     :: Vec
   !!--++    type(EoS_Type),              intent(in out) :: Eospar
   !!--++
   !!--++ PRIVATE
   !!--++ Pass values fron a vector to respective EoS parameter
   !!--++
   !!--++ Date: 28/02/2013
   !!
   Subroutine Vec_to_EoS(Vec,EosPar)
      !---- Arguments ----!
      real(kind=cp), dimension(:), intent(in)     :: Vec
      type(EoS_Type),              intent(in out) :: Eospar

      if (.not. eospar%linear) then
         eospar%params=vec
      else
         eospar%params(1)=vec(1)**(1.0_cp/3.0_cp)
         eospar%params(2:5)=vec(2:5)*3.0_cp
         select case(eospar%itherm)                 ! thermal expansion terms
            case(1,2,3)
               eospar%params(10:12)=vec(10:12)/3.0_cp

            case(4,5,6)
               eospar%params(10)=vec(10)/3.0_cp
               eospar%params(11)=vec(11)
         end select
      end if

      return
   End Subroutine Vec_to_EoS

   !!----
   !!---- SUBROUTINE Write_Info_Eos(Eospar, iout)
   !!----    type(Eos_Type),    intent(in) :: Eospar
   !!----    integer, optional, intent(in) :: iout
   !!----
   !!---- Subroutine that print information on iout unit
   !!----
   !!---- Date: 13/06/13
   !!
   Subroutine Write_Info_Eos(Eospar,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eospar
      integer, optional, intent(in) :: iout

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  EOS Information'
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' Title: '//trim(eospar%title)
      write(unit=lun,fmt='(a)') ' Model: '//trim(eospar%model)
      write(unit=lun,fmt='(a,i2)') ' Order: ',eospar%iorder
      if (eospar%linear) then
         write(unit=lun,fmt='(a)') ' Class: Linear'
      else
         write(unit=lun,fmt='(a)') ' Class: Volume'
      end if

      !> Pressure Parameters
      write(unit=lun,fmt='(a,f8.3)') ' Pressure of reference: ',eospar%pref
      write(unit=lun,fmt='(a)') ' '

      do i=1,4
         if (eospar%iuse(i) /= 0)then
            call setnum_std(eospar%params(i)*eospar%factor(i),eospar%esd(i)*eospar%factor(i),line)     ! include scaling
            string=' '
            if (eospar%iuse(i) == 2)string=' [IMPLIED VALUE]'
            write(unit=lun,fmt='(1x,a5,'': '',a,T30,'':'',a)') &
                 trim(eospar%parname(i)),trim(line),trim(eospar%comment(i))//trim(string)
         end if
      end do
      write(unit=lun,fmt='(a)') ' '

      !> Thermal EOS
      if (eospar%itherm > 0) call write_info_eos_thermal(eospar,lun)

      !> End
      write(unit=lun,fmt='(a)') ' '

      return
   End Subroutine Write_Info_Eos

   !!--++
   !!--++ SUBROUTINE Write_Info_Eos_Thermal(Eospar, iout)
   !!--++    type(Eos_Type),    intent(in) :: Eospar
   !!--++    integer, optional, intent(in) :: iout
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 13/06/13
   !!
   Subroutine Write_Info_Eos_Thermal(Eospar,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eospar
      integer, optional, intent(in) :: iout

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Init
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Thermal Expansion'
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' Model: '//trim(eospar%tmodel)
      write(unit=lun,fmt='(a)') ' '

      write(unit=lun,fmt='(a,f8.2,a)') ' Temperature of reference: ',eospar%tref,' K'
      write(unit=lun,fmt='(a)') ' '
      do i=5,n_eospar
         if (eospar%iuse(i) /= 0) then
             call setnum_std(eospar%params(i)*eospar%factor(i),eospar%esd(i)*eospar%factor(i),line)     ! include scaling
             string=' '
             if (eospar%iuse(i) > 1) string=' [FIXED VALUE]'
             write(unit=lun,fmt='(1x,a5,'': '',a,T30,'':'',a)') &
                  trim(eospar%parname(i)),trim(line),trim(eospar%comment(i))//trim(string)
         end if
      end do

      return
   End Subroutine Write_Info_Eos_Thermal

End Module CFML_EoS
