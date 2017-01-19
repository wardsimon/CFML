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
!!----    Validation for non-thermal EOS and modifications: RJA 22/02/2013 ..... 2014
!!----    Addition of code to handle Landau-type phase transitions: RJA 08/2014 - 02/2015
!!----    Modifications 06/2015 - 11/2015 for 2003 compatibility and Eos development RJA, JGP
!!----    Modifications 10/2015 - 07/2016 for curved phase boundaries and fitting moduli. RJA
!!
Module CFML_EoS
   !---- Use Modules ----!
   Use CFML_GlobalDeps,       only: cp, pi
   Use CFML_Crystal_Metrics,  only: Crystal_Cell_Type,Get_Cryst_Family,Volume_Sigma_from_Cell
   Use CFML_String_Utilities

   !---- Definitions ----!
   implicit none

   private

   !---- Public procedures ----!
   public :: Alpha_Cal, Dkdt_Cal, Get_GPT, Get_Grun_PT, Get_Grun_V, Get_K, Get_Pressure, Get_Pressure_Esd, Get_Temperature, &
             Get_Volume, Get_Volume_S, K_Cal, Kp_Cal, Kpp_Cal, Pressure_F, Strain, Strain_EOS, Get_Temperature_P0,          &
             Transition_phase, Get_Transition_Strain,Get_Transition_Temperature,Get_Transition_Pressure


   public :: Allocate_EoS_Data_List, Allocate_EoS_List, Calc_Conlev, Deallocate_EoS_Data_List, Deallocate_EoS_List,         &
             Deriv_Partial_P, EoS_Cal,  EoS_Cal_Esd, Eosparams_check, FfCal_Dat, FfCal_Dat_Esd, FfCal_EoS,                  &
             Init_EoS_Data_Type, Init_EoS_Type, Init_Err_EoS, Init_Eos_Thermal, Init_EoS_Transition, Init_Eos_Shear,        &
             Read_EoS_File, Read_EoS_DataFile, Set_Eos_Names, Set_Eos_Use, Set_Kp_Kpp_Cond, Write_Data_Conlev,              &
             Write_Info_Conlev, Write_EoS_File, Write_Eoscal,Write_EoS_DataFile, Write_Info_EoS, Read_Multiple_EoS_File


   !--------------------!
   !---- PARAMETERS ----!
   !--------------------!
   integer, public, parameter :: NCOL_DATA_MAX=22   ! Defines the maximum number of columns allowed in input data files

   integer, public, parameter :: N_EOSPAR=40        ! Specify the maximum number of Eos parameters allowed in Eos_type data type

   integer, public, parameter :: N_PRESS_MODELS=5   ! Number of possible pressure models
   integer, public, parameter :: N_THERM_MODELS=6   ! Number of possible Thermal models
   integer, public, parameter :: N_TRANS_MODELS=3   ! Number of possible Transition models
   integer, public, parameter :: N_SHEAR_MODELS=1   ! Number of possible Shear models
   integer, public, parameter :: N_DATA_TYPES=2     ! Number of possible data types (V,Kt,Ks etc)


   character(len=*), public, parameter, dimension(0:N_PRESS_MODELS) :: PMODEL_NAMES=(/    &    ! Name of the Pressure Models
                                                                       'None           ', &
                                                                       'Murnaghan      ', &
                                                                       'Birch-Murnaghan', &
                                                                       'Vinet          ', &
                                                                       'Natural Strain ', &
                                                                       'Tait           '/)

   character(len=*), public, parameter, dimension(0:N_THERM_MODELS) :: TMODEL_NAMES=(/        &  ! Name of the Thermal Models
                                                                       'None               ', &
                                                                       'Berman 1988        ', &
                                                                       'Fei 1995           ', &
                                                                       'Modified HP1998    ', &
                                                                       'Kroll              ', &
                                                                       'Salje low-T        ', &
                                                                       'HP Thermal Pressure'/)

   character(len=*), public, parameter, dimension(0:N_TRANS_MODELS) :: TRANMODEL_NAMES=(/ &    ! Name of Transition models
                                                                       'None           ', &
                                                                       'Landau P only  ', &
                                                                       'Landau T only  ', &
                                                                       'Landau PVT     '/)

   character(len=*), public, parameter, dimension(0:N_SHEAR_MODELS) :: SHEARMODEL_NAMES=(/ &    ! Name of Transition models
                                                                       'None           ', &
                                                                       'Polynomial     '/)


    character(len=*), public, parameter, dimension(0:N_DATA_TYPES) :: DATATYPE_NAMES=(/ &    ! Name of Data types
                                                                       'Cell parameters    ', &
                                                                       'Isothermal moduli  ', &
                                                                       'Adiabatic moduli   '/)

   real(kind=cp), public, parameter, dimension(6) :: DELCHI_LEVELS=(/68.30,90.00,95.40,99.00,99.73,99.99/) ! Confidence Levels
   real(kind=cp), public, parameter, dimension(6) :: DELCHI       =(/ 2.30, 4.61, 6.17, 9.21,11.80,18.40/) ! Delta Chi2 values

   !---------------!
   !---- TYPES ----!
   !---------------!

   !!----
   !!----  TYPE :: EOS_TYPE
   !!--..
   !!---- Update: 17/07/2015
   !!
   Type, public :: EoS_Type
      character(len=80)                         :: Title=" "         ! Descriptive title of EoS, set by user
      character(len=15)                         :: Model=" "         ! Murnaghan, Birch-Murnaghan, Vinet, Natural-Strain
      character(len=15)                         :: TModel=" "        ! Name for thermal model
      character(len=15)                         :: TranModel=" "     ! Name for phase transition model
      character(len=15)                         :: SModel=" "        ! Name for shear model
      character(len=5), dimension(N_EOSPAR)     :: ParName=" "       ! Names of the Eos variables...init
      character(len=50),dimension(N_EOSPAR)     :: Comment=" "       ! Description of the Eos variables inclduing units...init
      character(len=15)                         :: Pscale_name=" "   ! Description of the Pressure scale. Only used for output
      character(len=15)                         :: Vscale_name=" "   ! Description of the Volume scale. Only used for output
      character(len=120),dimension(10)          :: doc=" "           ! Documentation for eos: set by user
      integer                                   :: IModel=0          ! Index for Model
      integer                                   :: IOrder=0          ! Order for the Model
      logical                                   :: Linear=.false.    ! Flag for Linear EoS not volume
      integer                                   :: ITherm=0          ! Index for thermal expansion model, =0 for none
      integer                                   :: ITran=0           ! Index for phase transition model, =0 for none
      integer                                   :: IShear=0          ! Index for shear model, =0 for none
      integer,      dimension(N_EOSPAR)         :: Iuse=0            ! Flags for parameters allowed for a given EoS =0 (not), =1 (refineable), =2 (implied non-zero)
      real(kind=cp)                             :: PRef=0.0          ! Pressure of Reference
      real(kind=cp)                             :: TRef=298.0        ! Temperature of Reference
      real(kind=cp)                             :: Stoich=0.0        ! Stocihometry factor for multiple-phase calculations
      real(kind=cp)                             :: Density0=0.0      ! Density at reference conditions
      logical                                   :: TRef_fixed=.false.! If true, then Tref cannot be changed by user
      real(kind=cp), dimension(N_EOSPAR)        :: Params=0.0        ! EoS Parameters
      real(kind=cp), dimension(N_EOSPAR)        :: Esd=0.0           ! Sigma EoS Parameters
      real(kind=cp)                             :: X=0.0             ! Spare intensive variable, after P,T
      real(kind=cp)                             :: WChi2=0.0         ! weighted chi-squared for the refined parameters...set by ls
      real(kind=cp)                             :: DelPMax=0.0       ! Maximum misfit in P between obs and calc P...set by ls
      integer, dimension(4)                     :: IWt=0             ! Choice for weighting in LS: order is P,T,V,X, 1 = use sigmas for weights
      integer,      dimension(N_EOSPAR)         :: IRef=0            ! Refinement switches
      real(kind=cp),dimension(N_EOSPAR)         :: Factor=1.0        ! Scale factor to multiply variables for output (and divide on input)
      real(kind=cp),dimension(N_EOSPAR)         :: Lastshift=0.0     ! Shift applied in last LS cycle to parameters
      real(kind=cp),dimension(N_EOSPAR,N_EOSPAR):: VCV=0.0           ! Var-Covar matrix from refinement
   End Type EoS_Type

   !!----
   !!---- TYPE :: EOS_LIST_TYPE
   !!--..
   !!---- Update: 17/07/2015
   !!
   Type, public :: EoS_List_Type
      integer                                 :: N=0    ! Number of EoS List
      type(EoS_Type),dimension(:),allocatable :: EoS    ! EoS Parameters
   End Type EoS_List_Type

   !!----
   !!----  TYPE :: EOS_DATA_TYPE
   !!--..
   !!---- Update: 17/07/2015
   !!
   Type, public :: EoS_Data_Type
      integer                     :: IUse=0    ! 0=No active, 1= active
      integer                     :: IGrp=0    ! Group
      integer                     :: xtype     ! Indicates type of data in V,cell, etc xtype=0 default, xtype=1 isothermal moduli etc
      real(kind=cp)               :: T=298.0   ! Temperature
      real(kind=cp)               :: P=0.0     ! Pressure
      real(kind=cp)               :: V=0.0     ! Volume (xtype=0) or volume derivative as indicated by xtype
      real(kind=cp), dimension(3) :: cell=0.0  ! a,b,c parameters or derivatives as indicated by xtype
      real(kind=cp), dimension(3) :: ang=0.0   ! alpha, beta, gamma parameters
      real(kind=cp)               :: SigT=0.0  ! Sigma Temperature
      real(kind=cp)               :: SigP=0.0  ! Sigma Pressure
      real(kind=cp)               :: SigV=0.0  ! Sigma Volume
      real(kind=cp), dimension (3):: sigC=0.0  ! Sigma for a, b and c parameters
      real(kind=cp), dimension (3):: sigA=0.0  ! Sigma for alpha, beta and gamma angles
   End Type EoS_Data_Type

   !!----
   !!---- TYPE :: EOS_DATA_LIST_TYPE
   !!--..
   !!---- Update: January - 2013
   !!
   Type, public :: EoS_Data_List_Type
      character(len=80)                            :: Title=" "     ! Title of dataset (normally from input datafile)
      character(len=40)                            :: System=" "    ! Crystal System  (normally set by Def_Crystal_System)
      integer                                      :: N=0           ! Number of EoS Data List
      integer, dimension(ncol_data_max)            :: IC_Dat=0      ! Which values are input
      type(EoS_Data_Type),dimension(:),allocatable :: EoSD          ! EoS Data Parameters
   End Type EoS_Data_List_Type

   !-------------------!
   !---- VARIABLES ----!
   !-------------------!
   logical,            public :: Err_EOS=.false.  ! Logical Variable indicating an error in CFML_EoS module
   character(len=256), public :: Err_EOS_Mess=" " ! String containing information about the last error


Contains

   !!----
   !!---- FUNCTION ALPHA_CAL
   !!----
   !!---- Calculate the alpha parameter in thermal case
   !!---- For linear case alpha is correct as 1/a da/dT
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Alpha_Cal(P,T,Eospar, DeltaT) Result(alpha)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: P       ! Pressure
      real(kind=cp),            intent(in) :: T       ! Temperature
      type(Eos_Type),           intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T value

      !---- Local Variables ----!
      real(kind=cp)                  :: alpha,del,tt,delmin,tlimit,tr
      real(kind=cp), dimension(-2:2) :: v       ! array for calc v values
      integer                        :: j

      !> Init
      alpha=0.0_cp
      Tlimit=0.0_cp

      !> Need to trap numerical problems with Kroll, Salje, Pthermal at low T
      select case(eospar%itherm)
         case(0) ! no thermal parameters
            return

         case(4,5,6) ! Kroll and Pthermal: For T < 0.05T(Einstein) alpha=0. Same for Salje but test is 0.05Tsat
            if (t < 0.05_cp*eospar%params(11) ) return
      end select

      !> Numerical solutions: set step size
      del =abs(0.001_cp/eospar%params(10))      ! Set step in T to get about 0.1% shift in V
      if (del > 80.0) del=80.0_cp               ! otherwise for small alpha, del is too big and alpha is inaccurate
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

      !> Stop calculation going across a phase boundary
      if (eospar%itran > 0)then
         Tr=get_transition_temperature(p,eospar)
         if (transition_phase(P,T,eospar) .neqv. transition_phase(P,T+2.0*del,eospar)) del=abs(T-Tr)/2.1_cp
         if (transition_phase(P,T,eospar) .neqv. transition_phase(P,T-2.0*del,eospar)) del=abs(T-Tr)/2.1_cp
         if (del < 1.0) del=1.0
      end if

      !> Do the numerical solution
      do j=-2,2,1
         tt=t+real(j)*del                 ! apply shift to temp
         v(j)=get_volume(p,tt,eospar)     ! calc resulting V
      end do

      alpha=(v(-2)+8.0_cp*(v(1)-v(-1))-v(2))/(12.0_cp*del)/v(0)     ! Derivative to second order approximation

      !> Trap non-physical results that arise due to rounding errors
      select case(eospar%itherm)
         case(4,5,6)
            if (alpha < tiny(0.0)) alpha=0.0_cp
      end select

      !> get_volume returns 'a' for linear, so alpha is correct as 1/a da/dT for linear

      return
   End Function Alpha_Cal

   !!----
   !!---- FUNCTION DKDT_CAL
   !!----
   !!---- Calculate the derivative dK/dt (or dM/dt in the linear case) at P and T
   !!----
   !!---- Date: 17/07/2015
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
      del=30.0_cp                        ! good number for accuracy
      if (present(deltat)) del=deltat

      !> Volume
      vp=get_volume(p,t+del,eospar)
      vm=get_volume(p,t-del,eospar)

      !> In the linear case: this function return dM/dT
      dKdT=(k_cal(vp,t+del,eospar)-k_cal(vm,t-del,eospar))/(2.0_cp*del)

      return
   End Function dKdT_Cal

   !!----
   !!---- FUNCTION GET_GPT
   !!----
   !!---- Obtain the value of G (or Glinear) at P and T
   !!----
   !!---- Date: 11/07/2016
   !!
   Function Get_GPT(P,T,EoSPar) Result(gpt)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: p    ! Pressure
      real(kind=cp),  intent(in) :: t    ! Temperature
      type(Eos_Type), intent(in) :: EosPar    ! EoS Variable

      !---- Local Variables ----!
      integer                     :: i
      real(kind=cp)               :: gpt, delp

      !> default
      gpt=eospar%params(30)          ! default is g(Pref,Tref)

      !> T variation
      select case(eospar%ishear)       ! choice of model
         case(1)                  ! model 1 is polynomial in P and T
            gpt=gpt+(t-eospar%tref)*eospar%params(34)
            delp=p-eospar%pref
            do i=1,3
               if (eospar%params(i+30) < tiny(0.0)) exit
               gpt=gpt+eospar%params(i+30)*delp**i      ! eg linear term is (P-Pref)*dG/dP with params(31)
            end do

      end select

      return
   End Function Get_GPT

   !!----
   !!---- FUNCTION  GET_GRUN_PT
   !!----
   !!---- Returns Gruneisen parameter at this P,T as gamma0*(V/V0)
   !!----
   !!---- Date: 18/07/2016
   !!
   Function Get_Grun_PT(P,T,Eospar)  Result(G)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: T       ! Temperarture
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: V
      real(kind=cp) :: G

      V=get_volume(P,T,eospar)
      G=Get_Grun_V(V,Eospar)

      return
   End Function Get_Grun_PT

   !!----
   !!---- FUNCTION  GET_GRUN_V
   !!----
   !!---- Returns Gruneisen parameter at this volume as gamma0*(V/V0)
   !!---- If linear it calculates gamma0*(a/a0)^3
   !!----
   !!---- Date: 18/07/2016
   !!
   Function Get_Grun_V(V,Eospar)  Result(G)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume or length if linear
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: G
      real(kind=cp) :: V0,VV0

      !> Must be careful with transitions because eospar&params(1) is the high phase V0
      V0=get_volume(eospar%pref,eospar%tref,eospar)
      VV0=V/V0
      if (eospar%linear) VV0=VV0**3.0_cp
      if (abs(eospar%params(14)-1.0_cp) > 0.001_cp) VV0=VV0**eospar%params(14)

      G=eospar%params(13)*VV0

      return
   End Function Get_Grun_V

   !!----
   !!---- FUNCTION GET_K
   !!----
   !!---- Returns the value of K (or M if linear) at P and T
   !!----
   !!---- Works by using get_volume(P,T) and then using the V in k_cal
   !!----
   !!---- Date: 18/07/2016
   !!----
   Function Get_K(P,T,EosPar) Result(K)
      !---- Arguments ----!
      real(kind=cp),intent(in)        :: p       ! Pressure
      real(kind=cp),intent(in)        :: t       ! Tenperature
      type(Eos_Type),intent(in)       :: Eospar  ! Eos Parameters

      !---- Local Variables ----!
      real(kind=cp)                   :: k, v

      v=get_volume(p,t,eospar)
      k=k_cal(v,t,eospar)

      return
   End Function Get_K

   !!--++
   !!--++ FUNCTION Get_K0_T
   !!--++
   !!--++ PRIVATE
   !!--++ Returns k0 needed for Eos calculations at T this means for pthermal,
   !!--++ k at Tref is returned.
   !!--++ In the linear case then  M(T, P=0) from params is returned
   !!--++
   !!--++ If k is calculated as being negative, an error message is set
   !!--++ and k0 is returned as the value at Tref for safety.
   !!--++
   !!--++ Date:17/07/2015
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
   !!---- FUNCTION GET_PRESSURE
   !!----
   !!---- Returns value of pressure (P) at (V,T) for EoS Model defined in EoSPar variable
   !!---- In any error occurs then the final value is set to 0.0
   !!----
   !!--.. Note:
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 22/02/2013
   !!--.. Changed 28/03/2013 RJA to use strain function to get f from Vo/V.
   !!--.. This will ensure better consistency
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Get_Pressure(V,T,EosPar) Result(P)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      logical                           :: first
      integer                           :: i
      real(kind=cp)                     :: p
      real(kind=cp)                     :: K0,Kp,kpp,vv0,x,u,vol,plast,vs,ptr,vtr,difp,volp
      real(kind=cp)                     :: a,b,c,f
      real(kind=cp),dimension(n_eospar) :: ev
      real(kind=cp),dimension(3)        :: abc      ! for Tait parameters

      !> Init
      p=0.0_cp
      vol=v             ! needed to allow for transition strain: vol is the V of the bare phase
      vs=0.0
      plast=0.0
      first=.true.

      !> copy Eos parameters to local
      call EoS_to_Vec(eospar,ev)         ! Volume or linear case is covered

      !> These parameters depend only on T, not P or V
      k0=Get_K0_T(T,eospar)              ! Handles thermal pressure case, returns K0 or M0
      if (eospar%linear) k0=k0/3.0_cp
      kp=ev(3)
      kpp=ev(4)

      !> Start increment loop to get transition factor
      do
         !> Thermal case
         select case (eospar%itherm)
            case (0,6) ! No thermal case, or (6) thermal pressure which uses params at Tref
               vv0=vol/eospar%params(1)      !vv0 is now V/V0 or a/a0

            case (1:5)
               vv0=vol/Get_V0_T(T,EosPar)    ! In the case of Phase transition, Get_V0_T always returns V0 of high phase at this T: This is correct!
         end select

         f=strain(vv0,eospar)              ! use strain to get finite strain from v/v0 or a/a0
         if (err_eos)then
             err_eos_mess=trim(err_eos_mess)//' called from Get_Pressure'
             exit
         end if


         !> Using volume dimensions
         vv0=1.0_cp/vv0                    ! vv0 now V0/V for easy of use in rest of subprogram
         if (eospar%linear) vv0=vv0**(3.0_cp)

         select case (eospar%imodel)
            case (1) ! Murnaghan
               P=K0/Kp*(vv0**Kp - 1.0_cp)

            case (2) ! Birch-Murnaghan
               a=f*(1.0_cp+2.0_cp*f)**(2.5_cp)     ! changed expressions to use only f 04/04/2013 RJA
               b=0.0_cp
               c=0.0_cp
               if (eospar%iorder > 2)  b=1.5_cp*(Kp-4.0_cp)
               if (eospar%iorder == 4) c = 1.5_cp*(K0*Kpp + (Kp-4.0_cp)*(Kp-3.0_cp)+35.0_cp/9.0_cp)
               p=3.0_cp*K0*a*(1.0_cp + b*f + c*f*f)

            case (3) ! Vinet
               x=vv0**(-1.0_cp/3.0_cp)
               u=1.0_cp -x
               p=3.0_cp*K0*u/x/x*exp(1.5_cp*(Kp-1.0_cp)*u)

            case (4) ! Natural
               b=0.0_cp
               c=0.0_cp
               if (eospar%iorder > 2)  b=1.5_cp*(Kp-2.0_cp)
               if (eospar%iorder == 4) c =1.5_cp*(K0*Kpp + 1.0_cp +(Kp-2.0_cp)+(Kp-2.0_cp)**2.0_cp)
               p=3.0_cp*vv0*K0*f*(1.0_cp + b*f + c*f*f)

            case (5) ! Tait
               call get_tait(eospar,t,abc)
               vv0=1.0_cp/vv0     ! back to vv0=v/v0
               p=(((vv0 +abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3)) - 1.0_cp)/abc(2)
         end select

         !> Handle pthermal EoS
         if (eospar%itherm == 6 .and. eospar%imodel /= 0) p=p+pthermal(T,eospar)

         !> Iteration required if there is a phase transition
         if (eospar%itran == 0) exit

         !> Determine and set the transition pressure at this T
         if (first) then
            i=0
            if (eospar%itran == 1)then
               ptr=ev(21)       !PV transition, so Ptr is the eos%param(21)
            else
                !> ptr=(T-ev(21))/ev(22)
                ptr=get_transition_pressure(T,eospar)
            end if
            vtr=huge(0._cp)
         end if

         i=i+1

         if (abs(plast-p) < 0.0001) exit         ! iteration has converged
         if (i > 1000)exit

         !> Here if not converged:      ! works provided phase boundary is almost linear
         if (first .or. abs(plast-p) < difp) then
            vs=get_transition_strain(p,T,eospar)
            volp=vol
            vol=(v/(1.0+vs) +vol)/2.0
            difp=abs(plast-p)
            plast=p
         else            ! shift was too big
            vol=(vol+volp)/2.0
         end if

         first=.false.
      end do

      return
   End Function Get_Pressure

   !!----
   !!---- FUNCTION GET_PRESSURE_ESD
   !!----
   !!---- Gets the partial derivatives of P with respect to the EoS
   !!---- parameters at a given v,t point in order to calculate from vcv
   !!---- matrix the esd(p)
   !!----
   !!--.. Rewritten 09/04/2013 RJA to use deriv_partial
   !!----
   !!---- Update: 17/07/2015
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
   !!---- FUNCTION GET_TEMPERATURE
   !!----
   !!---- Returns Temperature at given P,V
   !!----
   !!---- Date: 01/04/2014
   !!
   Function Get_Temperature(P,V,EosPar) Result(T)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: V       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      integer                           :: nstep

      real(kind=cp)                     :: t
      real(kind=cp)                     :: pa,va,ta
      real(kind=cp)                     :: step,dp1,dp2

      !> Init
      t=eospar%tref
      pa=p               ! local copy p
      va=v
      call Init_Err_EoS()

      !> Check
      if (eospar%itherm == 0) return          ! no calcs possible without thermal model

      !> First estimate at P=0
      t=Get_Temperature_P0(va,EosPar)
      if (eospar%imodel ==0) return
      if (err_eos) return

      !> Use iterative solution, by minimising p-pcalc
      !> Init this part
      ta=t
      dp1=p-get_pressure(va,ta,eospar)
      nstep=0
      step=100.0_cp  ! For positive thermal expansion, should increase T from P=0

      !> Do iteration
      do
         !> Trap infinite loops
         if (nstep > 10000)then
             err_eos=.true.
             write(err_eos_mess,'("No convergence in get_temperature after ",i5," steps, dp= ",f10.4)')nstep,dp1
             exit
         end if

         !> Increment Temperature
         ta=ta+step
         dp2=p-get_pressure(va,ta,eospar)
         nstep=nstep+1

         !> test for sufficient convergence
         if (abs(step) < 1) exit          ! Go to 1 K precision

         !> not converged, so adjust step size
         if (dp1*dp2 < 0.0_cp) then
            !> overshot ptr:reverse step direction and make size smaller
            step=-0.5_cp*step
         else
            if (abs(dp2) > abs(dp1))then
               step=-1.0*step      ! wrong direction: reverse
            end if                 !(correct direction uses same step again)
         end if

         dp1=dp2        ! update delta-p values and go back for next cycle
      end do

      t=ta              ! success

      return
   End Function Get_Temperature

   !!----
   !!---- FUNCTION GET_TEMPERATURE_P0
   !!----
   !!---- Returns Temperature at P=0 for given V0T
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Get_Temperature_P0(Vin,EosPar,Tmin,Tmax) Result(Tk)
      !---- Arguments ----!
      real(kind=cp),           intent(in) :: Vin       ! Volume at temperature T or Pth (Pthermal case)
      type(Eos_Type),          intent(in) :: EoSPar    ! Eos Parameter
      real(kind=cp), optional, intent(in) :: Tmin      ! Range for solution in T
      real(kind=cp), optional, intent(in) :: Tmax

      !---- Local Variables ----!
      real(kind=cp)                      :: tk,tref
      real(kind=cp)                      :: v00,v0T
      real(kind=cp)                      :: a,b,c,d,t1,t2,t3,x1,x2,x3
      real(kind=cp)                      :: a1,a2,a3,q,r,s1,s2,dd,th
      real(kind=cp)                      :: y,p
      real(kind=cp)                      :: kp,th_e,eps0
      real(kind=cp)                      :: Tkmin,Tkmax
      real(kind=cp), dimension(n_eospar) :: ev

      !> Init
      Tk=eospar%tref

      !> Local copy Eospar to handle linear or volume
      call EoS_to_Vec(eospar,ev)                         ! all equations written for volume

      v00=ev(1)
      if (eospar%linear) then
         v0T=Vin**3.0_cp
      else
         v0T=Vin
      end if

      !> Optional arguments
      TKmin=0.0_cp
      if (present(Tmin)) TKmin=tmin

      TKmax=1.0e4
      if (present(Tmax)) TKmax=tmax

      !> Check
      Tref=eospar%tref
      if (v00 <= tiny(0.0)) return

      !> Init local variables
      t1=0.0
      t2=0.0
      t3=0.0

      select case (eospar%itherm)
         case (1) ! Berman
            ! modified RJA 31.03.2014 to fix case params(11)=0.
            if (abs(ev(11)) < tiny(0.0_cp)) then
               if (ev(10) > tiny(0.0_cp)) tk=Tref+(v0T/v00 -1.0_cp)/ev(10)
            else
               a= 0.5*ev(11)
               b= ev(10) - ev(11)*tref
               c= 1.0 - ev(10)*tref + 0.5*ev(11)*tref*tref - (v0t/v00)

               d= b*b - 4.0*a*c
               if (d >= 0.0) then
                  t1=(-b+sqrt(d))/(2.0*a)
                  t2=(-b-sqrt(d))/(2.0*a)
                  if (t1 > 0.0) then
                     if (t1 >=TKmin .and. t1 <=TKmax)then
                         tk=t1
                     else
                         err_eos=.true.
                         err_eos_mess='No valid solution for temperature'
                     end if

                  else if (t2 > 0.0) then
                     if (t2 >=TKmin .and. t2 <=TKmax) then
                         tk=t2
                     else
                         err_eos=.true.
                         err_eos_mess='No valid solution for temperature'
                     end if
                  else
                     err_eos=.true.
                     err_eos_mess='No valid solution for temperature'
                  end if
               end if
            end if

         case (2) ! Fei
            ! modified RJA 31.03.2014 to fix case params(11)=0.
            if (abs(ev(11)) < tiny(0.0_cp)) then
               if (ev(10) > tiny(0.0_cp))tk=Tref+log(v0T/v00)/ev(10)
            else
               a=0.5*ev(11)
               b=ev(10)
               c=-b*tref -a*tref*tref + ev(12)/tref -log(v0t/v00)
               d=-ev(12)

               a1=b/a
               a2=c/a
               a3=d/a
               q=(3.0*a2 - a1*a1)/9.0
               r=(9.0*a1*a2-27.0*a3-2.0*a1*a1*a1)/54.0

               dd=q**3+r**2
               if (dd >= 0.0) then
                  s1=(r + sqrt(q**3+r**2))**(1.0/3.0)
                  s2=(r - sqrt(q**3+r**2))**(1.0/3.0)

                  t1=s1+s2-a1/3.0
                  if (abs(s1-s2) <= tiny(0.0)) t2=-0.5*(s1+s2)-a1/3.0

                  if (t1 > 0.0) then
                     if (t1 >=TKmin .and. t1 <=TKmax)then
                         tk=t1
                     else
                         err_eos=.true.
                         err_eos_mess='No valid solution for temperature'
                     end if

                  elseif (t2 > 0.0) then
                     if (t2 >=TKmin .and. t2 <=TKmax)then
                         tk=t2
                     else
                         err_eos=.true.
                         err_eos_mess='No valid solution for temperature'
                     end if
                  else
                     err_eos=.true.
                     err_eos_mess='No valid solution for temperature'
                  end if

               else
                  th=acos(r/sqrt(-q**3))
                  t1=2.0*sqrt(-q)*cos(th/3.0)-a1/3.0
                  t2=2.0*sqrt(-q)*cos((th+2.0*pi)/3.0) -a1/3.0
                  t3=2.0*sqrt(-q)*cos((th+4.0*pi)/3.0) -a1/3.0
                  if (t1 > 0.0) then
                     if (t1 >=TKmin .and. t1 <=TKmax) tk=t1
                  elseif (t2 > 0.0) then
                     if (t2 >=TKmin .and. t2 <=TKmax) tk=t2
                  elseif (t3 > 0.0) then
                     if (t3 >=TKmin .and. t3 <=TKmax) tk=t3
                  else
                     err_eos=.true.
                     err_eos_mess='No valid solution for temperature'
                  end if
               end if
            end if

         case (3) ! HP
            ! modified RJA 31.03.2014 to fix coding errors
            a=ev(10)
            b=-2.0*(10.0*ev(10) + ev(11))
            c=1.0 - ev(10)*tref + 2.0*(10.0*ev(10) + ev(11))*sqrt(tref) - (V0t/V00)

            d=b*b - 4.0*a*c
            if (d >= 0.0) then
               x1=(-b+sqrt(d))/(2.0*a)
               x2=(-b-sqrt(d))/(2.0*a)
               if (x1 > 0.0) then
                  t1=x1*x1
                  if (t1 >=TKmin .and. t1 <=TKmax) tk=t1
               elseif (x2 > 0.0) then
                  t2=x2*x2
                  if (t2 >=TKmin .and. t2 <=TKmax) tk=t2
               end if
            else
               err_eos=.true.
               err_eos_mess='No valid solution for temperature'
            end if

         case (4) ! Kroll
            kp=ev(3)
            th_e=ev(11)
            a=th_e/Tref
            b=-1.0/(kp*(kp+2.0))
            eps0=(a**2)*exp(a)/(exp(a)-1.0)**2

            x1=(-b*(1.0+kp))*eps0/(ev(10)*th_e)
            x2=1.0 - (((v0t/v00)+kp)/(1.0+kp))**(1.0/b)
            x3=1.0/(exp(a)-1)

            c=1.0/(x1*x2+x3)
            if (c <= -1.0_cp)then
               err_eos=.true.
               err_eos_mess='Temperature calculated as negative'
            else
               d=log(1.0 +c)
               t1=th_e/d
               if (t1 >=TKmin .and. t1 <=TKmax)then
                  tk=t1
               else
                  err_eos=.true.
                  err_eos_mess='No valid solution for temperature'
               end if
            end if


         case (5) ! Salje
            x1=0.0
            x2=0.0
            x3=0.0

            y=ev(10)*ev(11)
            p=v00**(1.0/3.0) -y

            a=y**3
            b=3.0 * y*y*p
            c=3.0 * y * p*p
            d=p**3 - v0t

            a1=b/a
            a2=c/a
            a3=d/a
            q=(3.0*a2 - a1*a1)/9.0
            r=(9.0*a1*a2-27.0*a3-2.0*a1*a1*a1)/54.0

            dd=q**3+r**2
            if (dd >= 0.0) then
               s1=(r + sqrt(q**3+r**2))**(1.0/3.0)
               s2=(r - sqrt(q**3+r**2))**(1.0/3.0)

               x1=s1+s2-a1/3.0
               if (abs(s1-s2) <= tiny(0.0)) then
                  x2=-0.5*(s1+s2)-a1/3.0
               end if
            else
               th=acos(r/sqrt(-q**3))
               x1=2.0*sqrt(-q)*cos(th/3.0)-a1/3.0
               x2=2.0*sqrt(-q)*cos((th+2.0*pi)/3.0) -a1/3.0
               x3=2.0*sqrt(-q)*cos((th+4.0*pi)/3.0) -a1/3.0
            end if
            if (abs(x1) >= 1.0) t1=2.0*ev(11)/log((1.0+x1)/(x1-1.0))
            if (abs(x2) >= 1.0) t2=2.0*ev(11)/log((1.0+x2)/(x2-1.0))
            if (abs(x3) >= 1.0) t3=2.0*ev(11)/log((1.0+x3)/(x3-1.0))

            if (t1 > 0.0) then
               if (t1 >=TKmin .and. t1 <=TKmax) tk=t1
            end if
            if (t2 > 0.0 .and. tk < 0.0) then
               if (t2 >=TKmin .and. t2 <=TKmax) tk=t2
            end if
            if (t3 > 0.0 .and. tk < 0.0) then
               if (t3 >=TKmin .and. t3 <=TKmax) tk=t3
            end if

         case(6)
            th_e=ev(11)
            a=th_e/Tref
            eps0=(a**2)*exp(a)/(exp(a)-1.0)**2

            x1=eps0/(ev(10)*ev(2)*ev(11))
            x2=1.0/(exp(a)-1.0)

            c=1.0/(v0t*x1+x2)
            d=log(1.0 + c)

            t1=th_e/d
            if (t1 >=TKmin .and. t1 <=TKmax) tk=t1

      end select

      return
   End Function Get_Temperature_P0

   !!----
   !!---- FUNCTION GET_TRANSITION_PRESSURE
   !!----
   !!---- Returns the transition pressure at this T
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Get_Transition_Pressure(T,EosPar) Result(Ptr)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)              :: Ptr       ! The transition T at this P
      real(kind=cp)              :: sqroot

      !> init
      ptr=0._cp
      Err_EoS=.false.
      call Init_Err_EoS()

      !> Check for valid model number. If not valid, return with zero Tr
      if (eospar%itran < 1 .or. eospar%itran > N_TRANS_MODELS) return

      !> Now determine if in high-symmetry or low field
      select case(eospar%itran)
         case(1) ! Landau PV
            ptr=eospar%params(21)

         case(3) ! Landau PVT
            ! ptr = (T-eospar%params(21))/eospar%params(22) original linear
            if (abs(eospar%params(23)) < tiny(0.0)) then
               ptr = (T-eospar%params(21))/eospar%params(22)        ! linear
            else
               sqroot=eospar%params(22)*eospar%params(22)-4.0_cp*eospar%params(23)*(eospar%params(21)-T)
               if (sqroot > tiny(0.0)) then
                  ptr=2.0_cp*(T-eospar%params(21))/(eospar%params(22)+sqrt(sqroot))        ! Viet/Muller formula for root
               else if(sqroot > -1.0*tiny(0.0))then
                  ptr=-2.0_cp*(T-eospar%params(21))/eospar%params(22)
               else
                  Err_EoS=.true.
                  Write(Err_EoS_Mess,'(a,f8.3,a)')'No real solution for Ptr at T = ',T,'K'
               end if
            end if

      end select

      return
   End Function Get_Transition_Pressure

   !!----
   !!---- FUNCTION GET_TRANSITION_STRAIN
   !!----
   !!---- Returns the strain at P,T due to the transition, including any softening
   !!---- in the high-symm phase for Volume eos returns the volume strain term, for
   !!---- linear eos it returns the linear term!!
   !!---- Vs is defined relative to the 'bare' eos of the high phase (ie the high
   !!---- phase without transition effects)Returns the transition pressure at this T
   !!----
   !!---- Date: 16/02/2015
   !!
   Function Get_Transition_Strain(P,T,EosPar) Result(vs)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !----Local Variables ----!
      real(kind=cp)              :: Vs      ! The volume strain
      real(kind=cp)              :: Ttr , a ! transition temperature at this pressure

      !> init
      vs=0._cp

      !> Check for valid model number. If not valid, return with zero Vs
      if (eospar%itran < 1 .or. eospar%itran > N_TRANS_MODELS) return

      !> Now determine if in high-symmetry or low field
      if (Transition_phase(P,T,Eospar)) then
         !> This section for being in the low field
         select case(eospar%itran)
            case(1) ! Landau PV
               vs=eospar%params(24)*abs(eospar%params(21)-P)**eospar%params(25)

            case(2) ! Landau TV
               vs=eospar%params(24)*abs(eospar%params(21)-T)**eospar%params(25)

            case(3) ! Landau PVT
               Ttr = Get_Transition_Temperature(P,EosPar)
               a=eospar%params(24)
               vs=a*abs(Ttr-T)**eospar%params(25)            !abs function to handle highT being low sym
         end select

      else
         !> This section for being in the high field
         select case(eospar%itran)
            case(1) ! Landau PV
               vs=eospar%params(26)*abs(eospar%params(21)-P)**eospar%params(27)

            case(2) ! Landau TV
               vs=eospar%params(26)*abs(eospar%params(21)-T)**eospar%params(27)

            case(3) ! Landau PVT:  Note no da/dP for highP phase
               Ttr = Get_Transition_Temperature(P,EosPar)
               vs=eospar%params(26)*abs(Ttr-T)**eospar%params(27)            !abs function to handle highT being low sym
         end select
      end if

      return
   End Function Get_Transition_Strain

   !!----
   !!---- FUNCTION GET_TRANSITION_TEMPERATURE
   !!----
   !!---- Returns the transition temperature at this pressure
   !!----
   !!---- Date: 16/02/2015
   !!
   Function Get_Transition_Temperature(P,EosPar) Result(Tr)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)              :: Tr       ! The transition T at this P

      !>init
      tr=0._cp

      !> Check for valid model number. If not valid, return with zero Tr
      if (eospar%itran < 1 .or. eospar%itran > N_TRANS_MODELS) return

      !> Now determine if in high-symmetry or low field
      select case(eospar%itran)
         case(2) ! Landau TV
            Tr=eospar%params(21)

         case(3) ! Landau PVT: with a curved phase boundary
            Tr = eospar%params(21)+p*eospar%params(22)+p*p*eospar%params(23)

      end select

      return
   End Function Get_Transition_Temperature

   !!--++
   !!--++ FUNCTION Get_V0_T
   !!--++
   !!--++ PRIVATE
   !!--++ Returns the volume at P=0 and T=T from V0 and Thermal expansion
   !!--++ Except for Pthermal, for which it returns V at P=0, T=Tref
   !!--++ It calculates V0 only for thermal expansion, not including transition effects!
   !!--++ Therfore this must remain PRIVATE
   !!--++
   !!--++ Date: 17/07/2015
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

         case(1)                    ! Berman, works at all T
            V=ev(1)*(1.0_cp +ev(10)*delt+0.5_cp*ev(11)*delt*delt)

         case(2)                    ! Fei,
            TT=T
            delt2=t*t-tref*tref
            if (tt < tiny(0._cp)) tt=0.001        ! to prevent divide by zero

            if (abs(ev(12)) > tiny(0._cp)) then
               V=ev(1)*exp(ev(10)*delt + 0.5_cp*ev(11)*delt2 - ev(12)*(1.0_cp/TT - 1.0_cp/Tref))
            else
               V=ev(1)*exp(ev(10)*delt + 0.5_cp*ev(11)*delt2)    ! to protect from divide by zero when alpha2=0 and T=0
            end if

         case(3)                    ! HP 1998 modified
            tt=t
            if (t < tiny(0._cp)) tt=tiny(0._cp)      ! prevents sqrt(-ve T)
            V=ev(1)*(1.0_cp + ev(10)*(tt-tref) - 2.0_cp*(10.0_cp*ev(10)+ev(11))*(sqrt(tt) -sqrt(tref)))

         case(4)                    ! Holland-Powell 2011 in the Kroll form
            Tn= ev(11)/Tref                        ! theta/Tref
            C=Tn*Tn*exp(Tn)/(exp(tn)-1.0_cp)**2.0_cp
            B=-1.0/ev(3)/(ev(3)+2.0_cp)
            if (t > 0.05_cp*ev(11)) then                               ! to avoid numerical problems at T=0
               A=ev(10)*ev(11)/C *(1.0_cp/(exp(ev(11)/T)-1.0_cp) - 1.0_cp/(exp(Tn)-1.0_cp) )
            else
               A=ev(10)*ev(11)/C *(-1.0_cp/(exp(Tn)-1.0_cp) )          ! because when T=0 then 1.0/exp(Tein/T) = 0
            end if
            V=ev(1)*(-1.0_cp*ev(3) + (1.0_cp+ev(3))*(1.0_cp - ev(3)*(ev(3)+2.0_cp)*A/(ev(3)+1.0_cp))**B)

         case(5)                    ! Salje, Tref fixed at zero
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
   !!---- FUNCTION GET_VOLUME
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
      real(kind=cp), parameter          :: PREC=0.000001_cp  !precision to find V.

      integer                           :: nstep
      type(Eos_Type)                    :: EoS  ! Eos copy
      real(kind=cp)                     :: V
      real(kind=cp)                     :: V0,K0,Kp,k,strain
      real(kind=cp)                     :: Vol, step,dp1,dp2
      real(kind=cp),dimension(n_eospar) :: ev
      real(kind=cp),dimension(3)        :: abc          ! Tait parameters
      real(kind=cp)                     :: pa           ! pa=p-pth

      !> Init
      v=0.0_cp
      strain=0.       ! strain from transition: linear or volume to match eos type
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
      if (eospar%linear) k0=k0/3.0_cp                 ! correct M0 to K0 for linear
      kp=ev(3)

      !> Get the volume strain due to transition: only a function of P,T NOT V!!
      if (eospar%itran > 0) then
         strain=get_transition_strain(P,T,eospar)     ! returns the linear or volume strain as appropriate
      end if

      !> If there is no eos model, we are finished because V0 is the V at this T
      if (eospar%imodel == 0) then
         if (eospar%linear) v0=v0**(1.0_cp/3.0_cp)
         v=v0*(1.0_cp + strain)
         return
      end if

      !> Analytic solution for Murnaghan:  use this for first guess for other EoS except Tait
      v=v0*(1.0_cp + kp*pa/k0)**(-1.0_cp/kp)
      if (eospar%imodel ==1) then
         if (eospar%linear) v=v**(1.0_cp/3.0_cp)
         if (eospar%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
         return
      end if

      !> Analytic solution for Tait
      if (eospar%imodel ==5) then
         call get_tait(eospar,t,abc)                     ! get_tait returns volume-like parameters even for linear
         v=v0*(1.0_cp-abc(1)*(1.0_cp-(1.0_cp + abc(2)*pa)**(-1.0_cp*abc(3))))
         if (eospar%linear) v=v**(1.0_cp/3.0_cp)
         if (eospar%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
         return
      end if

      vol=v
      if (eospar%linear) vol=vol**(1.0_cp/3.0_cp)

      !> Find iterative solution for the rest of functions: get_pressure includes the thermal pressure term
      !> But if there is a transition, we only want the P/V for the bare high-symm phase without softening
      eos=eospar        ! copy
      eos%itran=0       ! turn off transition

      dp1=p-get_pressure(vol,t,eos)

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
         dp2=p-get_pressure(vol,t,eos)
         nstep=nstep+1

         !> test for sufficient convergence
         if (abs(step) < PREC*Vol) exit          ! 1 part in 100,000 in volume

         !> not converged, so adjust step size
         if (dp1*dp2 < 0.0_cp) then
            !> overshot ptr:reverse step direction and make size smaller
            step=-0.5_cp*step
         else
            if (abs(dp2) > abs(dp1)) then
               step=-1.0*step      ! wrong direction: reverse
            else
               step=0.9_cp*dp2/dp1*step   ! correct direction, should get a smaller step size
            end if                        ! the factor of 0.9 is a safety factor to ensure step gets smaller
         end if

         dp1=dp2        ! update delta-p values and go back for next cycle
      end do

      !> now set return value depending on success or not
      v=vol           ! success
      if (eospar%itran > 0) v=vol*(1.0_cp + strain)  ! apply transition strain ('vol' is actually linear if linear eos)

      return
   End Function Get_Volume

   !!----
   !!---- FUNCTION GET_VOLUME_S
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
            V=v0*(1.0_cp-s)**3.0_cp

         case (4) ! Natural Strain
            V=v0*exp(-3.0_cp*s)
      end select

      !> Linear case
      if (eospar%linear) v=v**(1.0_cp/3.0_cp)

      return
   End Function Get_Volume_S

   !!----
   !!---- FUNCTION K_CAL
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
      real(kind=cp)                      :: kc
      real(kind=cp)                      :: vv0,k0,kp,kpp,cvv0,p,vol1,f,vol
      real(kind=cp)                      :: a,b,c,nu
      real(kind=cp)                      :: vs,Ttr,dVs                   ! for transition calculations
      real(kind=cp), dimension(n_eospar) :: ev
      real(kind=cp), dimension(3)        :: abc     ! Tait parameters

      !> Init
      kc=0.0_cp

      !> Correct the volume to the high phase only if there is a transition
      if (eospar%itran > 0) then
         p=get_pressure(v,t,eospar)
         Vs=get_transition_strain(P,T,eospar)     ! returns the linear or volume strain as appropriate
         Vol=v/(1.0+vs)
      else
         Vol=v
      end if

      !> Get pressure corrected for pthermal
      p=get_pressure(v,t,eospar)-pthermal(t,eospar)

      !> Set appropriate V/V0
      select case (eospar%itherm)
         case (0)                               ! No thermal
            vv0=vol/eospar%params(1)            ! vv0 or aa0
         case (1:5)
            vv0=vol/get_V0_T(t,eospar)          ! vv0 is  v(p,t)/v(p=0,t)
         case (6)
            vol1=get_volume(p,eospar%tref,eospar)    ! get volume at p-pth, at t=tref
            vv0=vol1/eospar%params(1)                ! vv0 is now for eoS at Tref
      end select

      !> Now get local copies of EoS parameters appropriate for EoS and thermal type
      k0=Get_K0_T(T,eospar)               ! returns M(T) for linear, handles pthermal
      if (err_eos) return                 ! exit with value eosparms(2) if k0 at P=0 calculated as negative

      if (eospar%linear) k0=k0/3.0_cp

      !> Local EoS copy
      call EoS_to_Vec(eospar,ev)

      kp=ev(3)
      kpp=ev(4)

      !> Strain for BM, NS, Vinet EoS equations: this is f at p-pth if pthermal used
      f=strain(vv0,eospar)
         if (err_eos) then
            err_eos_mess=trim(err_eos_mess)//' called from Kcal'
            return
         endif

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

      !> To this point Kc is the bulk modulus of the 'bare' eos of the high phase
      !> Now handle phase transition

      if (eospar%itran > 0) then
         p=get_pressure(V,T,eospar)
         Vs=get_transition_strain(P,T,eospar)
         if (eospar%linear) Vs=3.0_cp*Vs    ! Here we want to work in volume-terms

         select case(eospar%itran)
            case(1) ! Landau P-V
               if (transition_phase(P,T,eospar)) then             ! in the low phase
                  dVs=ev(24)*ev(25)*abs(P-ev(21))**(ev(25)-1)     ! correct if lowP phase is highsym
                  if (nint(ev(20)) ==1) dVs=-1.0_cp*dVs            ! change sign for highp phase is highsym
                  kc=1.0_cp/(1.0_cp/kc - dVs/(1+Vs))
               else
                                                                  ! in the high phase
                  dVs=ev(26)*ev(27)*abs(P-ev(21))**(ev(27)-1)     ! correct if highP phase is highsym
                  if (nint(ev(20)) /= 1) dVs=-1.0_cp*dVs
                  kc=1.0_cp/(1.0_cp/kc - dVs/(1+Vs))
               end if

            case(2,3) ! Landau VT or PVT
               if (transition_phase(P,T,eospar)) then
                  ! low phase
                  Ttr = get_transition_temperature(P,eospar)
                  dVs=ev(24)*ev(25)*(abs(Ttr-T))**(ev(25)-1)*(ev(22)+2.0_cp*ev(23)*p) ! for curved phase boundary
                  if (nint(ev(20)) /= 1) dVs=-1.0_cp*dVs                              ! sign change when low phase is highT
                  kc=1.0_cp/(1.0_cp/kc - dVs/(1+Vs))
               else
                  ! in the highphase
                  Ttr = get_transition_temperature(P,eospar)
                  dVs=ev(26)*ev(27)*(abs(Ttr-T))**(ev(27)-1)*ev(22)                     ! simpler because no da/dP: expression if High phase is lowT   ! for straight boundary
                  dVs=ev(26)*ev(27)*(abs(Ttr-T))**(ev(27)-1)*(ev(22)+2.0_cp*ev(23)*p)   ! simpler because no da/dP: expression if High phase is lowT: curved boundary
                  if (nint(ev(20)) ==1) dVs=-1.0_cp*dVs                                 ! sign change when high phase is highT
                  kc=1.0_cp/(1.0_cp/kc - dVs/(1+Vs))
               end if
         end select

      end if

      if (eospar%linear) kc=kc*3.0_cp

      return
   End Function K_Cal

   !!----
   !!---- FUNCTION KP_CAL
   !!----
   !!---- Returns value of kprime at this volume for EoS
   !!----
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Kp_Cal(V,T,EoSpar) Result (kpc)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: kpc
      real(kind=cp)                      :: vv0,k0,kp,kpp,vol1,p,ptr
      real(kind=cp)                      :: a,b,f,rkp_top, rkp_bot,nu
      real(kind=cp), dimension(n_eospar) :: ev
      real(kind=cp), dimension(3)        :: abc     ! Tait parameters
      real(kind=cp) :: delv,vol2,rk1,rk2,p1,p2

      !> Init
      kpc=0.0_cp

      !> Local copy
      call EoS_to_Vec(eospar,ev)

      !> Get pressuure corrected for pthermal
      p=get_pressure(v,t,eospar)-pthermal(t,eospar)

      !> Set appropriate V/V0  so that strain can be calculated
      select case (eospar%itherm)
         case (0)                                    ! No thermal
            vv0=v/eospar%params(1)                   ! vv0 or aa0
         case (1:5)
            vv0=v/get_V0_T(t,eospar)                 ! vv0 is  v(p,t)/v(p=0,t)
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
      if (err_eos) then
         err_eos_mess=trim(err_eos_mess)//' called from Kp_cal'
         return
      end if

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

      !> Code to handle transition calculation, numerically: rewrote this 17/11/2014
      if (eospar%itran > 0) then
         delv=0.0001_cp*v

         !> Fix delv so that we do not step over into other phase
         Ptr=get_transition_pressure(t,eospar)
         p1=get_pressure(v,t,eospar)           ! The pressure at VT point of interest
         p2=get_pressure(v-delv,t,eospar)
         if (transition_phase(P1,T,eospar) .neqv. transition_phase(P2,T,eospar)) delv=0.9*abs(get_volume(ptr,t,eospar)-v)

         p2=get_pressure(v+delv,t,eospar)
         if(transition_phase(P1,T,eospar) .neqv. transition_phase(P2,T,eospar)) delv=0.9*abs(get_volume(ptr,t,eospar)-v)

         vol1=v+delv
         vol2=v-delv
         rk1=k_cal(vol1,t,eospar)
         rk2=k_cal(vol2,t,eospar)
         p1=get_pressure(vol1,t,eospar)
         p2=get_pressure(vol2,t,eospar)

         !> Trap case when p2=p1 because there is no eos loaded
         if(abs(p2-p1) > tiny(0.0_cp)) kpc=(rk2-rk1)/(p2-p1)
      end if

      return
   End Function Kp_Cal

   !!----
   !!---- FUNCTION KPP_CAL
   !!----
   !!---- Returns value of kdouble-prime at this volume for EoS
   !!----
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Kpp_Cal(V,T,EoSpar) Result (kppc)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: kppc,ptr
      real(kind=cp) :: delv,vol1,vol2,rk1,rk2,p1,p2

      !> Init
      kppc=0.0_cp

      !> jump on equation type
      if (eospar%imodel == 1 .and. eospar%itran == 0) return ! Murnaghan without transition

      !> numerical solution
      delv=0.0001_cp*v

      !> Code to prevent stepping across transition
      if (eospar%itran > 0)then
         Ptr=get_transition_pressure(t,eospar)
         p1=get_pressure(v,t,eospar)           ! The pressure at VT point of interest

         p2=get_pressure(v-delv,t,eospar)
         if (transition_phase(P1,T,eospar) .neqv. transition_phase(P2,T,eospar)) delv=0.9*abs(get_volume(ptr,t,eospar)-v)

         p2=get_pressure(v+delv,t,eospar)
         if (transition_phase(P1,T,eospar) .neqv. transition_phase(P2,T,eospar)) delv=0.9*abs(get_volume(ptr,t,eospar)-v)
      end if

      vol1=v+delv
      vol2=v-delv
      rk1=kp_cal(vol1,t,eospar)
      rk2=kp_cal(vol2,t,eospar)
      p1=get_pressure(vol1,t,eospar)
      p2=get_pressure(vol2,t,eospar)

      !> Trap case when p2=p1 because there is no eos loaded
      if (abs(p2-p1) > tiny(0.0_cp)) kppc=(rk2-rk1)/(p2-p1)

      !> No linear conversion is required because kp_cal returns values for "linear Kp" = Mp,
      !> so kppc is already dMp/dP = Mpp

      return
   End Function Kpp_Cal

   !!--++
   !!--++ FUNCTION NORMPRESSURE_EOS
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
      !> 17/11/2014: Before transitions, this used eosparams for normal thermal expansion
      !> but with transitions safer to do following:

      if (eospar%itherm == 0 .and. eospar%itran == 0) then
         ! simple PV eos
         k0=ev(2)
         kp=ev(3)
         kpp=ev(4)

      else if(eospar%itran == 0 .and. eospar%itherm /= 6) then
         ! normal thermal expansion models with dK/dT and no transition
         k0=Get_K0_T(T,eospar)                     ! returns M(T) for linear,
         if (eospar%linear) k0=k0/3.0_cp
         kp=ev(3)
         kpp=ev(4)

      else
         ! Transition model, or Pthermal: cannot use get_V0_T or get_K0_T because they return V and K at Tref for pthermal
         v0=get_volume(0.0_cp,T,eospar)        ! determine V0 at this T
         k0=k_cal(v0,T,eospar)
         kp=kp_cal(v0,T,eospar)
         kpp=kpp_cal(v0,T,eospar)
         if (eospar%linear) then
            k0=k0/3.0_cp
            kp=kp/3.0_cp
            kpp=kpp/3.0_cp
         end if
      end if

      select case(eospar%imodel)
         case (1,5) ! Murnaghan, Tait
            f=0.0_cp

         case (2) ! Birch-Murnaghan
            b=0.0_cp
            c=0.0_cp
            if (eospar%iorder > 2) b=1.5_cp*(kp-4.0_cp)
            if (eospar%iorder ==4) c=1.5_cp*(k0*kpp + (kp-4.0_cp)*(kp-3.0_cp)+35.0_cp/9.0_cp)
            f=K0*(1.0_cp + b*s + c*s*s)

         case (3) ! Vinet: new definition RJA 28/03/2013
            f= K0*exp(1.5_cp*(Kp-1.0_cp)*s)

         case (4) ! Natural Strain
            b=0.0_cp
            c=0.0_cp
            if (eospar%iorder > 2) b=1.5_cp*(Kp - 2.0_cp)
            if (eospar%iorder ==4) c=1.5_cp*(K0*Kpp + 1.0_cp + (Kp-2.0_cp)+(Kp-2.0_cp)**2.0_cp)
            f=K0*(1.0_cp + b*s + c*s*s)
      end select

      return
   End Function NormPressure_Eos

   !!--++
   !!--++ FUNCTION NORMPRESSURE_P
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
   !!---- FUNCTION PRESSURE_F
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
   !!--++ FUNCTION PTHERMAL
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
   !!---- FUNCTION STRAIN
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
      call init_err_eos
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
   !!---- FUNCTION STRAIN_EOS
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

   !!----
   !!---- LOGICAL FUNCTION TRANSITION_PHASE
   !!----
   !!---- Returns .true. if P and T are in the low phase stability field
   !!---- and .false. (default) if in the high-symm field, or exactly on
   !!---- the boundary.
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Transition_Phase(P,T,Eospar) Result(Ip)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      logical                    :: Ip

      !---- Local Variables ----!
      real(kind=cp)              :: Ttr      ! transition temperature at this pressure

      !> default to 'high' phase for safety
      ip=.false.

      !> Check for valid model number. If not valid, return (with high phase indicated).
      if (eospar%itran < 1 .or. eospar%itran > N_TRANS_MODELS) return

      !> Test P and T against the Tr, set ip=.true. if in the low phase field for highT =high symm
      select case(eospar%itran)
         case(1) ! Landau PV
            if (P < eospar%params(21)) ip=.true.

         case(2) ! Landau TV
            if (T < eospar%params(21)) ip=.true.

         case(3) ! Landau PVT
          !  Ttr = eospar%params(21)+p*eospar%params(22)  changed from this on 18/12/2015
            Ttr = Get_Transition_Temperature(P,EosPar)  ! general case
            if ( (T <  Ttr) ) ip=.true.

      end select

      !> Now invert result if the lowT phase is high symm phase:
      if (eospar%params(20) < 0) ip = .not. ip

      return
   End Function Transition_Phase

   !---------------------!
   !---- SUBROUTINES ----!
   !---------------------!

   !!----
   !!---- SUBROUTINE ALLOCATE_EOS_DATA_LIST
   !!----
   !!----    Allocation of objet E of eos_list_data_type.
   !!----    This subroutine should be called before using an object of type eos_data_list.
   !!----
   !!---- Update: 17/07/2015
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
   !!---- SUBROUTINE ALLOCATE_EOS_LIST
   !!----
   !!----    Allocation of objet E of eos_list_type.
   !!----    This subroutine should be called before using an object of type eos_list.
   !!----
   !!---- Update: 17/07/2015
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
   !!---- SUBROUTINE CALC_CONLEV
   !!----
   !!----    Performs the calculations of one confidence ellipse
   !!----    for a pair of EoS parameters refined by least-squares
   !!----
   !!--.. Added to cfml_eos_mod: 05/10/2015 RJA
   !!----
   !!---- Date: 18/07/2016
   !!
   Subroutine Calc_Conlev(Eos,ix,iy,isig,xyy,n)
      !---- Arguments ----!
      type(Eos_Type),              intent(in)  :: Eos          ! EoS with refined parameters
      integer,                     intent(in)  :: ix,iy        ! input pointers to two variables in the variance-covariance matrix
      integer,                     intent(in)  :: isig         ! requested confidence level, values 1-6
      real(kind=cp),dimension(:,:),intent(out) :: xyy          ! output points for plotting; array must by of dimension >100
      integer,                     intent(out) :: n            ! number of data output

      !---- Local Variables ----!
      integer                    :: ilast,nlimit
      real(kind=cp)              :: c11,c22,c12,det
      real(kind=cp)              :: cinv11,cinv22,cinv12,a,b,c,root
      real(kind=cp)              :: xe,xs,xinc,x,y1,y2

      !> Init error variables
      call Init_Err_Eos()

      !> init
      nlimit=size(xyy,dim=2)

      n=0
      xyy=0.0_cp

      !> Check
      if (isig < 1 .or. isig > 6) then
         err_eos=.true.
         err_eos_mess="Confidence level is out of range"
         return
      end if

      !> Copy over vcv values
      c11=eos%vcv(ix,ix)
      c22=eos%vcv(iy,iy)
      c12=eos%vcv(ix,iy)

      !> Invert matrix
      det=c11*c22-c12*c12
      if (abs(det) <=tiny(0.0_cp)) then
         err_eos=.true.
         err_eos_mess="Determinant value is zero in the Confidence ellipses calculation"
         return
      end if

      cinv22=c11/det
      cinv11=c22/det
      cinv12=-1.0_cp*c12/det

      xe=sqrt(cinv22*delchi(isig)/(cinv22*cinv11-cinv12**2.0_cp))    !x at the end
      xs=-1.0_cp*xe
      xinc=sqrt(c11)/10.0_cp

      !> at x=xs equal roots, calculate explicity
      x=xs
      y1= -1.0_cp*cinv12*xs/cinv22  ! gives equal roots
      y2=y1

      ilast=0
      n=n+1
      xyy(1,n)=x+eos%params(ix)
      xyy(2,n)=y1+eos%params(iy)
      xyy(3,n)=y2+eos%params(iy)

      !> loop over x starts here
      do
         if (abs(x-xe) < xinc .or. abs(x-xs) < xinc) then
            x=x+0.2*xinc     ! increment x
         else                ! small step close to end
            x=x+xinc
         end if

         if (x > xe) then               ! last of this delchi
            ilast=1                     ! calc for x=xe which
            y1= -1.0*cinv12*xe/cinv22   ! gives equal roots
            y2=y1
            x=xe
         else
            ! solve for y
            a=cinv22
            b=2.0*cinv12*x
            c=cinv11*x*x-delchi(isig)
            root=b*b-4.0*a*c
            if (root < 0.0) cycle                  ! error: try next x
            y1= (sqrt(root)-b)/2.0/a
            y2= -1.0*(sqrt(root)+b)/2.0/a
         end if

         n=n+1
         xyy(1,n)=x+eos%params(ix)
         xyy(2,n)=y1+eos%params(iy)
         xyy(3,n)=y2+eos%params(iy)
         if (ilast /= 0) exit

         if (n == nlimit) then
            err_eos=.true.
            err_eos_mess="Number of points arrived to the limit for Confidence ellipses"
            exit
         end if
      end do

      return
   End Subroutine Calc_Conlev

   !!----
   !!---- SUBROUTINE DEALLOCATE_EOS_DATA_LIST
   !!----
   !!----    De-allocation of objet E of type eos_data_list.
   !!----    This subroutine should be after using an object of type eos_data_list that is no
   !!----    more needed.
   !!----
   !!---- Update: 17/07/2015
   !!
   Subroutine Deallocate_EoS_Data_List(E)
      !---- Arguments ----!
      type (eos_data_list_type), intent(in out)   :: E  ! Objet to be deallocated

      if (allocated(E%eosd)) deallocate (E%eosd)
      E%n=0

      return
   End Subroutine Deallocate_EoS_Data_List

   !!----
   !!---- SUBROUTINE DEALLOCATE_EOS_LIST
   !!----
   !!----    De-allocation of objet E of type eos_list.
   !!----    This subroutine should be after using an object of type eos_list that is no
   !!----    more needed.
   !!----
   !!---- Update: 17/07/2015
   !!
   Subroutine Deallocate_EoS_List(E)
      !---- Arguments ----!
      type (eos_list_type), intent(in out)   :: E  ! Objet to be deallocated

      if (allocated(E%eos)) deallocate (E%eos)
      E%n=0

      return
   End Subroutine Deallocate_EoS_List

   !!--++
   !!--++ SUBROUTINE DEF_CRYSTAL_SYSTEM
   !!--++
   !!--++ PRIVATE
   !!--++ Either sets cell parameters to conform to specifed system Or,
   !!--++ if no system specified, tries to determine system if all cell parameters
   !!--++ provided
   !!--++
   !!--++ Update: 17/07/2015
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
               if (index(u_case(system),' C ') /= 0) then
                  !> alpha = beta = 90�
                  dat%eosd(1:ndat)%ang(1)=90.0
                  dat%eosd(1:ndat)%ang(2)=90.0
                  dat%eosd(1:ndat)%siga(1)=0.0
                  dat%eosd(1:ndat)%siga(2)=0.0
               else
                  !> alpha = gamma = 90�
                  dat%eosd(1:ndat)%ang(1)=90.0
                  dat%eosd(1:ndat)%ang(3)=90.0
                  dat%eosd(1:ndat)%siga(1)=0.0
                  dat%eosd(1:ndat)%siga(3)=0.0
               end if

            case ('ORTH')
               !> Angles =90�
               dat%eosd(1:ndat)%ang(1)=90.0
               dat%eosd(1:ndat)%ang(2)=90.0
               dat%eosd(1:ndat)%ang(3)=90.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

            case ('TETR')
               !> Angles =90�
               dat%eosd(1:ndat)%ang(1)=90.0
               dat%eosd(1:ndat)%ang(2)=90.0
               dat%eosd(1:ndat)%ang(3)=90.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

               !> b=a
               dat%eosd(1:ndat)%cell(2)=dat%eosd(1:ndat)%cell(1)
               dat%eosd(1:ndat)%sigc(2)=dat%eosd(1:ndat)%sigc(1)

            case ('TRIG','HEXA')
               !> Angles alpha=beta=90�, gamma=120�
               dat%eosd(1:ndat)%ang(1)= 90.0
               dat%eosd(1:ndat)%ang(2)= 90.0
               dat%eosd(1:ndat)%ang(3)  =120.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

               !> b=a
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
               !> Angles =90�
               dat%eosd(1:ndat)%ang(1)=90.0
               dat%eosd(1:ndat)%ang(2)=90.0
               dat%eosd(1:ndat)%ang(3)=90.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

               !> a=b=c: Modified RJA 14 Jan to handle case if V is supplied, but 'a' is not
               do i=1,ndat
                   if(dat%eosd(i)%cell(1) < tiny(0.))then
                       dat%eosd(i)%cell(1)=dat%eosd(i)%v**(1.0_cp/3.0_cp)
                       dat%eosd(i)%sigc(1)=dat%eosd(i)%sigc(1)/3.0_cp
                   endif
                   dat%eosd(i)%cell(2)=dat%eosd(i)%cell(1)
                   dat%eosd(i)%cell(3)=dat%eosd(i)%cell(1)
                   dat%eosd(i)%sigc(2)=dat%eosd(i)%sigc(1)
                   dat%eosd(i)%sigc(3)=dat%eosd(i)%sigc(1)
               enddo
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
         if (len_trim(system) <= 0) then
            err_eos=.true.
            Err_EoS_Mess="Cell values for first data in the data file are inconsistent"
            return
         end if
         !> Now check if all data the same. If they are not, set system=' ' because there may be a transition
         do i=2,ndat
            ncell%cell=dat%eosd(i)%cell
            ncell%ang =dat%eosd(i)%ang
            call Get_Cryst_Family(ncell,Family,Symbol,SystemC)
            systemC=adjustl(systemc)

            if (systemC /= System) then
!                err_eos=.true.
!                write(Err_EoS_Mess,'(a,i5,a)')"Cell values for data #",i," in the data file indicate a change in crystal system"
                dat%system=' '
                return
            end if
         end do
      end if

      return
   End Subroutine Def_Crystal_System

   !!----
   !!---- SUBROUTINE DERIV_PARTIAL_P
   !!----
   !!---- Calculate Partial derivates of Pressure respect to Params
   !!----
   !!---- Date: 17/07/2015
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
      if (eospar%itran ==0) then
         td(1:5)=tda(1:5)                   ! analytic for Vo and moduli terms because these are exact even at small P
         td(6:n_eospar)=tdn(6:n_eospar)
      else
         td(1:n_eospar)=tdn(1:n_eospar)
      end if

      return
   End Subroutine Deriv_Partial_P

   !!--++
   !!--++ SUBROUTINE DERIV_PARTIAL_P_ANALYTIC
   !!--++
   !!--++ PRIVATE
   !!--++ Calculates the partial derivatives of P with respect to the EoS
   !!--++ at a given v,t point, and returns  them in array td
   !!--++
   !!--.. 27-Feb-2013: RJA edits to use IREF as in original
   !!--++
   !!--++ Date: 17/07/2015
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
      if (eospar%linear) then
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
            if (eospar%iorder > 2) b=1.5_cp*K0*(ev(3)-4.0_cp)
            if (eospar%iorder == 4)c = 1.5_cp*K0*(K0*ev(4) + (ev(3)-4.0_cp)*(ev(3)-3.0_cp)+35.0_cp/9.0_cp)

            !> dP/dVo
            td(1)=f52/vt0 * ( a + (7.0_cp*a+2.0_cp*b)*f + (9.0_cp*b+3.0_cp*c)*f*f + 11.0_cp*c*f*f*f)

            !> dP/dKo:  da, db, dc are da/dparam etc for each param K0, Kp0, Kpp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder > 2) db = 1.5_cp*(ev(3)-4.0_cp)
            if (eospar%iorder == 4)dc = 3.0_cp*K0*ev(4) + 1.5_cp*((ev(3)-4.0_cp)*(ev(3)-3.0_cp)+35.0_cp/9.0_cp)

            td(2)= 3.0*f52*f* (1.0 + db*f + dc*f*f)

            !> dP/dKp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder > 2) db = 1.5_cp*K0
            if (eospar%iorder == 4)dc = 1.5_cp*K0*(2.0_cp*ev(3)-7.0_cp)

            td(3) = 3.0*f52* (db + dc*f) *f*f

            !> dP/dKpp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder == 4) then
               dc = 1.5_cp*K0*K0
               td(4) = 3.0*f52*dc*f*f*f
            else
               td(4)=0.0_cp
            end if

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
            if (eospar%iorder > 2) a=1.5_cp*(ev(3)-2.0_cp)
            if (eospar%iorder ==4) b=1.5_cp*(1.0_cp + K0*ev(4) + (ev(3)-2.0_cp) + (ev(3)-2.0_cp)**2.0_cp)

            !> dP/dVo
            td(1) = K0/cv * (1.0_cp  + (2.0_cp*a+3.0_cp)*f + 3.0_cp*(a+b)*f*f +3.0_cp*b*f*f*f)

            !> d(p)/d(k0)
            td(2) = 3.0_cp * vv0 * (f + a*f*f + (b + 1.5_cp*ev(4)*k0)*f*f*f)

            !> dp/dKp0
            da = 0.0_cp
            db = 0.0_cp
            if(eospar%iorder > 2)  da=1.5_cp
            if (eospar%iorder == 4)db = 3.0_cp*ev(3) - 4.5_cp

            td(3) = 3.0_cp * k0 * vv0 * (da + db*f)*f*f

            !> dp/dKpp0
            td(4)=0.0_cp
            if (eospar%iorder == 4) td(4) = 4.5_cp* k0*k0 * vv0 *f*f*f

         case(5) ! Tait: new code and validated against Murnaghan and against finite diffs 15/06/2013
            call get_tait(eospar,t,abc)
            vv0=1.0/vv0        ! now vv0 is v/v0

            !> dP/dVo
            td(1)= k0 * vv0**2.0_cp/cv * ( (vv0 + abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3) -1.0_cp)

            !> d(p)/d(k0)
            td(2)=abc(1)*abc(3)*(( (vv0 + abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3)) -1.0_cp)

            !> dp/dKp0
            if (eospar%iorder > 2) then
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
         if (eospar%tref > 1.0_cp) deltinv=1.0_cp/t-1.0_cp/eospar%tref

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
      if (eospar%linear) then
         td(1) =td(1) *(3.0_cp*eospar%params(1)**2.0_cp)
         td(2:5)=td(2:5)/3.0_cp

         select case(eospar%itherm)                 ! thermal expansion terms
            case(1:2)                               ! Berman, Fei, alpha terms
               td(10:12)=td(10:12)*3.0_cp

            case(4)                                 ! Holland-Powell alpha
               td(10)=td(10)*3.0_cp
         end select                                 !other thermal equations have nothing to convert
      end if

      return
   End Subroutine Deriv_Partial_P_Analytic

   !!--++
   !!--++ SUBROUTINE DERIV_PARTIAL_P_NUMERIC
   !!--++
   !!--++ PRIVATE
   !!--++ Calculates the partial derivatives of P with respect to the EoS
   !!--++ at a given v,t point, and returns  them in array td
   !!--++
   !!--.. 27-Feb-2013: RJA edits to use IREF as in original
   !!--++
   !!--++ Date: 17/07/2015
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
               if (icycle > 0 .and. &
                   abs(d_prev-td(i))/td(i) < 1.0E-4) exit iter            ! deriv converged to 1 part in 10^4

               d_prev=td(i)                ! store last deriv value
               del=2.0_cp*del              ! increase the shift
               icycle=icycle+1
               if (icycle > 5) exit iter   ! Do not allow 2*shift to exceed 64% of param value
            end do iter
         end if
      end do

      !> no need to fix derivatives for linear eos by this method

      return
   End Subroutine Deriv_Partial_P_Numeric

   !!----
   !!---- SUBROUTINE EOS_CAL
   !!----
   !!---- Returns elastic properties (not the parameter values) at this P,T for EoS
   !!----
   !!---- Date: 17/07/2015
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
   !!---- SUBROUTINE EOS_CAL_ESD
   !!----
   !!---- Returns esd's of the elastic properties (not the parameter values) at this P and T for EoS
   !!----
   !!---- Date: 17/07/2015
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

   !!----
   !!---- SUBROUTINE EOSPARAMS_CHECK
   !!----
   !!---- Check for Params that are invalid for all cases.
   !!----
   !!--.. NOTE:
   !!--..      Written 7-April-2014 to fix invalid eos parameters and report as error message
   !!--..      This is different from physical_check, because that checks properties at specific
   !!--..      p and T
   !!----
   !!---- Date: 11/07/2016
   !!
   Subroutine EoSParams_Check(EoSPar)
      !---- Argument ----!
      type (EoS_Type), intent(in out) :: EoSPar

      !---- Local Variables ----!
      real(kind=cp)     :: pinf
      character(len=40) :: text

      !> Init
      err_eos=.false.
      err_eos_mess=''

      !> Check that v0 is positive
      if (eospar%params(1) < tiny(0.0))then
         eospar%params(1)=1.0_cp
         err_eos=.true.
         if (eospar%linear) then
            err_eos_mess='*****WARNING: a0 was < 0. Not allowed! Reset to 1.00'
         else
            err_eos_mess='*****WARNING: V0 was < 0. Not allowed! Reset to 1.00'
         end if
      end if

      !> Check K0 is positive for V-P (-ve K ok for linear)
      if (.not. eospar%linear .and. eospar%imodel > 0)then
         if (eospar%params(2) < tiny(0.0_cp))then
            eospar%params(2)=10.0_cp
            err_eos=.true.
            if (len_trim(err_eos_mess) == 0) then
               err_eos_mess=' *****WARNING: K0 was < 0. Not allowed! Reset to 10.0'
            else
               err_eos_mess=trim(err_eos_mess)//' And K0 was < 0. Not allowed! Reset to 10.0'
            end if
         end if
      end if

      !> Check Tref is positive
      if (eospar%tref < -1.0_cp*tiny(0.0_cp)) then
         eospar%tref=0.0_cp
         err_eos=.true.
         if (len_trim(err_eos_mess) == 0) then
            err_eos_mess=' *****WARNING: Tref was < 0. Not allowed! Reset to 0 K'
         else
            err_eos_mess=trim(err_eos_mess)//' And Tref was < 0. Not allowed! Reset to 0 K'
         end if
      end if


      select case(eospar%itherm)  ! for specific thermal parameters
         case(4,6)    !>Kroll or Pthermal must have P_einstein > 0
            if (eospar%params(11) < 0.1) then
               eospar%params(11)=eospar%Tref
               if (eospar%Tref < 0.1) eospar%params=0.1
               err_eos=.true.
               if (len_trim(err_eos_mess) == 0) then
                  err_eos_mess=' *****WARNING: Einstein T was =< 0. Not allowed! Reset to Tref'
               else
                  err_eos_mess=trim(err_eos_mess)//' And Einstein T was =< 0. Not allowed! Reset to Tref'
               end if
            end if
      end select

      !> Produce warning for curved phase boundaries: Pinflection = a/-2b when Ttr=Tr0+aP+bP^2
      if (eospar%itran>0 .and. abs(eospar%params(24)) > tiny(0.)) then
         pinf=-1.0*eospar%params(22)/2.0/eospar%params(24)
         if (abs(pinf) < 10.0) then
            err_eos=.true.
            write(text,'(a,f4.1,1x,a)')'Phase boundary inflects at P ~ ',pinf,trim(eospar%pscale_name)
            if (len_trim(err_eos_mess) == 0) then
               err_eos_mess=' *****WARNING: '//trim(text)
            else
               err_eos_mess=trim(err_eos_mess)//' And '//trim(text)
            end if
         end if
      end if

      return
   End Subroutine EoSParams_Check

   !!--++
   !!--++ SUBROUTINE EOS_TO_VEC
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

      !> Copy everything 1 to 1 as default
      vec=eospar%params

      !> adjust any linearr values as necessary
      if (eospar%linear) then
         vec(1)=eospar%params(1)**3.0_cp
         vec(2:5)=eospar%params(2:5)/3.0_cp
         select case(eospar%itherm)                 ! thermal expansion terms
            case(1:3)
               vec(10:12)=eospar%params(10:12)*3.0_cp

            case(4,5,6)
               vec(10)=eospar%params(10)*3.0_cp
               vec(11)=eospar%params(11)
         end select

         select case(eospar%itran)          ! phase transition
            case(1,2)      ! V or T only models
               vec(24)=3.0_cp*eospar%params(24)
               vec(26)=3.0_cp*eospar%params(26)

            case(3)        ! PVT model
               vec(23:24)=3.0_cp*eospar%params(23:24)
               vec(26)=3.0_cp*eospar%params(26)
         end select

      end if

      return
   End Subroutine EoS_to_Vec

   !!----
   !!---- SUBROUTINE FFCAL_DAT
   !!----
   !!---- Return Normalized pressure and/or Strain at V and P
   !!----
   !!---- Date: 17/07/2015
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
   !!---- SUBROUTINE FFCAL_DAT_ESD
   !!----
   !!---- Calculate the Normalised Pressure (F) and Strain (S) value at
   !!---- Volume (V) and Pressure (P) and also their ESD
   !!----
   !!--.. IMPORTANT: Eosparameters are not used in this calculation!
   !!--.. The only element of esopar that is used is the eos model type
   !!--.. and the linear flag when the strain function is called
   !!----
   !!---- Date: 17/07/2015
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
   !!----Subroutine FFCAL_EOS
   !!----
   !!----  Returns normalised pressure and Strain at this P,T,
   !!----  calculated from the eos parameters
   !!----
   !!---- Date: 17/07/2015
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

   !!----
   !!---- SUBROUTINE GET_TAIT
   !!----
   !!---- Returns a,b,c Tait parameters in a vector.
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Get_Tait(Eospar,T,Vec)
      !---- Arguments ----!
      type(Eos_Type),            intent(in)  :: Eospar  ! Eos Parameters
      real(kind=cp),             intent(in)  :: T       ! Temperature
      real(kind=cp),dimension(3),intent(out) :: Vec     ! Vector (a,b,c) of Tait parameters

      !---- Local Variables ----!
      real(kind=cp) :: k0,kp,kpp

      !> Init
      Vec=0.0_cp
      if (eospar%imodel /= 5) return

      select case(eospar%itherm)
         case(1:5)                  ! normal thermal expansion models with dK/dT
            k0 =Get_K0_T(T,eospar)
            kp =eospar%params(3)

         case default               ! includes no thermal model, also pthermal which requires params at Tref
            k0 =eospar%params(2)
            kp =eospar%params(3)
      end select

      if (eospar%iorder < 4) then
         kpp= -1.0_cp*kp/k0        ! implied value for Kpp except for 4th order
      else
         kpp=eospar%params(4)
      end if

      if (eospar%linear) then
         k0  =k0/3.0_cp
         kp  =kp/3.0_cp
         kpp =kpp/3.0_cp
      end if

      Vec(1)=(1.0_cp + kp)/(1.0_cp+kp+k0*kpp)
      Vec(2)= kp/k0
      if (abs(1_cp+kp) > tiny(0.0)) Vec(2)=Vec(2)-kpp/(1.0_cp+kp)

      Vec(3)= (1.0_cp+kp+k0*kpp)/(kp*kp+kp-k0*kpp)

      return
   End Subroutine Get_Tait

   !!----
   !!---- SUBROUTINE Init_EoS_Data_Type
   !!----
   !!---- Initialize EoS_Data_Type
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Init_EoS_Data_Type(E)
      !---- Arguments ----!
      type (EoS_Data_Type), intent(in out)   :: E

      E%IUse = 0
      E%IGrp = 0
      E%T    = 298.0
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
   !!---- SUBROUTINE INIT_EOS_Shear
   !!----
   !!---- Initialize the EoS Type for Shear case
   !!----
   !!---- Date: 17/02/2015
   !!
   Subroutine Init_EoS_Shear(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !---- Variables ----!
      integer  :: n

      !> Check for valid model number. If not valid, set zero
      if (eospar%ishear < 0 .or. eospar%ishear > N_SHEAR_MODELS) eospar%ishear=0

      !> Set upper limit to parameter numbers
      n=34
      if (n > n_eospar) n=n_eospar

      select case(eospar%ishear)
         case(0)
            eospar%params(30)        = huge(0.0_cp)   ! default is infinitely stiff
            eospar%vcv(30,1:n)         = 0.0_cp
            eospar%vcv(30:n_eospar,1:n)= 0.0_cp
            eospar%vcv(1:n_eospar,30:n)= 0.0_cp
            eospar%factor(30:n)        = 1.0_cp

         case(1)        ! polynomial
            eospar%params(30) = 100.0_cp       ! G0 at Pref Tref
            eospar%params(31:34) = 0.0_cp      ! Polynomial coefficients
            eospar%factor(30:n) = 1.0_cp
      end select

      call Set_Shear_Names(Eospar)    ! Set the variable names
      call Set_Eos_Use(Eospar)        ! update the use flags

      return
   End Subroutine Init_EoS_Shear

   !!----
   !!---- SUBROUTINE INIT_EOS_THERMAL
   !!----
   !!---- Initialize the EoS Type for Thermal case
   !!----
   !!---- Date: 10/09/2013
   !!
   Subroutine Init_EoS_Thermal(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !---- Variables ----!
      integer    :: n

      !> Check for valid model number. If not valid, set zero
      if(eospar%itherm < 0 .or. eospar%itherm > N_THERM_MODELS) eospar%itherm=0

      !> Set upper limit to thermal parameter numbers
      n=19
      if (n > n_eospar)n=n_eospar

      select case(eospar%itherm)
         case(0)
            eospar%params(5)    = 0.0_cp
            eospar%params(10:n) = 0.0_cp
            eospar%vcv(5,1:n)   = 0.0_cp
            eospar%vcv(10:n,1:n)= 0.0_cp
            eospar%vcv(1:n,10:n)= 0.0_cp
            eospar%factor(10:n) = 1.0_cp
            eospar%TRef         = 298.0_cp            ! Simple thermal expansion,
            eospar%TRef_fixed   = .false.

         case(1)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp
            eospar%TRef        = 298.0_cp             ! Simple thermal expansion,
            eospar%TRef_fixed  = .false.              ! alpha terms can be safely 0.

         case(2)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp
            eospar%factor(12)  = 1.0_cp
            eospar%TRef        = 298.0_cp             ! Simple thermal expansion,
            eospar%TRef_fixed = .false.

         case(3)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E4_cp
            eospar%TRef        = 298.0_cp             ! Simple thermal expansion,
            eospar%TRef_fixed  = .false.

         case(4)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef        = 298.0_cp             ! Holland and Powell thermal expansion
            eospar%TRef_fixed  = .false.
            eospar%params(11)  = 298.0_cp             ! Einstein temperature

         case(5)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef        = 0.0_cp               ! Salje thermal expansion
            eospar%TRef_fixed  = .true.
            eospar%params(11)  = 298.0_cp             ! Saturation temperature

         case(6)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef        = 298.0_cp
            eospar%TRef_fixed  = .false.
            eospar%params(11)  = 298.0_cp             ! Einstein temperature

      end select

      !> Set the common terms for the Gruneisen model
      eospar%params(13)=1.0      ! gamma0
      eospar%params(14)=1.0      ! q

      call Set_Thermal_Names(Eospar)                  ! Set the variable names
      call Set_Eos_Use(Eospar)                        ! update the use flags

      return
   End Subroutine Init_EoS_Thermal

   !!----
   !!---- SUBROUTINE INIT_EOS_TRANSITION
   !!----
   !!---- Initialize the EoS Type for Transition case
   !!----
   !!---- Date: 17/02/2015
   !!
   Subroutine Init_EoS_Transition(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !---- Variables ----!
      integer  :: n

      !> Check for valid model number. If not valid, set zero
      if (eospar%itran < 0 .or. eospar%itran > N_TRANS_MODELS) eospar%itran=0

      !> Set upper limit to parameter numbers
      n=29
      if (n > n_eospar) n=n_eospar

      select case(eospar%itran)
         case(0)
            eospar%params(20:n)        = 0.0_cp
            eospar%vcv(20,1:n)         = 0.0_cp
            eospar%vcv(20:n_eospar,1:n)= 0.0_cp
            eospar%vcv(1:n_eospar,20:n)= 0.0_cp
            eospar%factor(20:n)        = 1.0_cp

         case(1)        ! Landau PV
            eospar%params(20) = 1.0_cp        ! high P is high sym
            eospar%params(21) = 5.0_cp        ! Safe default Ptr
            eospar%params(22) = 0.0_cp        ! dTr/dP: not used
            eospar%params(23) = 0.0_cp        ! d2Tr/dP2 not used
            eospar%params(24) = 0.0_cp        ! no excess V
            eospar%params(25) = 0.5_cp        ! power law
            eospar%params(26) = 0.0_cp        ! excess V high phase
            eospar%params(27) = 0.5_cp        ! power law high phase

            eospar%factor(20:n) = 1.0_cp
            eospar%factor(24)   = 1.0E3_cp    ! 1000 for aL
            eospar%factor(26)   = 1.0E3_cp    ! 1000 for aH

         case(2)        ! Landau TV
            eospar%params(20) =   1.0_cp        ! high T is high sym
            eospar%params(21) = 800.0_cp        ! Safe default Ttr
            eospar%params(22) = 100.0_cp        ! dTr/dP: dummy will be used in get_pressure
            eospar%params(23) =   0.0_cp        ! d2Tr/dP2 not used
            eospar%params(24) =   0.0_cp        ! no excess V
            eospar%params(25) =   0.5_cp        ! power law
            eospar%params(26) =   0.0_cp        ! excess V high phase
            eospar%params(27) =   0.5_cp        ! power law high phase

            eospar%factor(20:n) = 1.0_cp
            eospar%factor(24)   = 1.0E3_cp      ! 1000 for aL
            eospar%factor(26)   = 1.0E3_cp      ! 1000 for aH

         case(3)        ! Landau PVT
            eospar%params(20) =   1.0_cp        ! high T is high sym
            eospar%params(21) = 800.0_cp        ! Safe defaulat Tr
            eospar%params(22) =   1.0_cp
            eospar%params(23) =   0.0_cp
            eospar%params(24) =   0.0_cp        ! no excess V
            eospar%params(25) =   0.5_cp
            eospar%params(26) =   0.0_cp        ! excess V high phase
            eospar%params(27) =   0.5_cp        ! power law high phase

            eospar%factor(20:n) = 1.0_cp
            eospar%factor(24)   = 1.0E3_cp      ! 1000 for aL
            eospar%factor(23)   = 1.0E5_cp      ! for da/dP
            eospar%factor(26)   = 1.0E3_cp      ! 1000 for aH

      end select

      call Set_Transition_Names(Eospar)    ! Set the variable names
      call Set_Eos_Use(Eospar)             ! update the use flags

      return
   End Subroutine Init_EoS_Transition

   !!----
   !!---- SUBROUTINE INIT_EOS_TYPE
   !!----
   !!---- Initialize EoS_Type setting all parameters to sensible values
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Init_EoS_Type(Eospar,CLin,IThermal,ITransition,Ishear)
      !---- Arguments ----!
      type (EoS_Type),            intent(in out) :: Eospar       ! EoS Type
      character(len=*), optional, intent(in)     :: CLin         ! Character variable to indicate linear EoS or not
      integer,          optional, intent(in)     :: IThermal     ! integer to indicate ithermal type
      integer,          optional, intent(in)     :: ITransition  ! integer to indicate transition type
      integer,          optional, intent(in)     :: IShear       ! integer to indicate shear type

      !> test for optional argument for linear
      Eospar%Linear  =.false.
      if (present(clin)) then
         if (index(U_case(clin(1:3)),'LIN') > 0) Eospar%linear=.true.
      end if

      Eospar%Title   =' '
      Eospar%IModel  =0
      Eospar%IOrder  =3

      eospar%ParName=' '
      eospar%comment=' '
      eospar%doc=' '
      eospar%pscale_name=' '
      eospar%vscale_name=' '
      call Set_Eos_Names(Eospar)         ! also sets the print factors for the pressure part

      Eospar%PRef     = 0.0_cp
      Eospar%Density0 = 0.0_cp

      Eospar%Iuse     =0
      Eospar%Iuse(1:4)=1                 ! Vo, Ko, Kpp

      Eospar%params   = 0.0_cp
      Eospar%esd      = 0.0_cp

      Eospar%params(1)= 1.0_cp
      Eospar%params(2)=10.0_cp
      Eospar%params(3)= 4.0_cp

      Eospar%X        = 0.0_cp
      Eospar%stoich   = 1.0_cp

      Eospar%WChi2    = 0.0_cp
      Eospar%DelPMax  = 0.0_cp
      Eospar%IWt      = 0

      Eospar%IRef     = 0
      Eospar%factor   = 1.0_cp
      Eospar%LastShift= 0.0_cp
      Eospar%VCV      = 0.0_cp

      !> Test for optional argument for thermal
      Eospar%ITherm  = 0
      Eospar%Tref    = 298.0_cp         ! Sensible default
      if (present(ithermal) .and. ithermal > -1 )then
         Eospar%Itherm = ithermal
         call Init_Eos_Thermal(Eospar)              ! set up default values, names for specific thermal eqn
      end if

      !> Test for optional argument for transition
      Eospar%ITran  =0
      if (present(itransition) .and. itransition  > -1 )then
         Eospar%ITran = itransition
         call Init_Eos_Transition(Eospar)              ! set up default values, names for specific transition
      end if

      !> Test for optional argument for shear model
      Eospar%Ishear  =0
      if (present(ishear) .and. ishear  > -1 )then
         Eospar%Ishear = ishear
      end if
      call Init_Eos_Shear(Eospar)              ! set up default values, names for specific shear model

      return
   End Subroutine Init_EoS_Type

   !!----
   !!---- SUBROUTINE INIT_ERR_EOS
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
   !!--++ SUBROUTINE PHYSICAL_CHECK
   !!--++
   !!--++ PRIVATE
   !!--++ Check if the parameters have physical sense
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Physical_Check(P,T,E)
      !---- Arguments ----!
      real(kind=cp),    intent(in) :: p  ! Pressure
      real(kind=cp),    intent(in) :: t  ! Temperature
      type(Eos_Type),   intent(in) :: E  ! EoS object

      !---- Local variables ----!
      character(len=20)   :: car
      real(kind=cp)       :: tlimit,v,pinf

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

      !> Produce warning for curved phase boundaries: Pinflection = a/-2b when Ttr=Tr0+aP+bP^2
      if (e%itran>0 .and. abs(e%params(23)) > tiny(0.0))then
         pinf=abs(e%params(22)/2.0/e%params(23))
         if (abs(p/pinf -1.0) < 0.1) then
            err_eos=.true.
            err_eos_mess='P in region of boundary inflection P: PVT calculations may be inaccurate or wrong'
         end if
      end if

      return
   End Subroutine Physical_Check

   !!----
   !!---- SUBROUTINE READ_EOS_DATAFILE
   !!----
   !!---- General routine to read data for Eos
   !!----
   !!---- Update: 17/07/2015
   !!
   Subroutine Read_EoS_DataFile(fname,dat)
      !---- Arguments ----!
      character(len=*),          intent(in)  :: fname   ! File name
      type (eos_data_list_type), intent(out) :: dat     ! data structure

      !---- Local Variables ----!
      character(len=255), dimension(:), allocatable :: flines
      character(len=255)                            :: line
      character(len=5)                              :: car
      character(len=1)                              :: Ts
      character(len=30), dimension(ncol_data_max)   :: dire

      integer                                       :: nldata,ndat, npos
      integer                                       :: i,j,kk,nl,nk,nlines,m,idatatype,nlines_datum
      integer                                       :: iv, inum
      integer, dimension(ncol_data_max)             :: ivet,iorden

      real(kind=cp), dimension(ncol_data_max)       :: vet,vetsd
      real(kind=cp), dimension(ncol_data_max)       :: rvet
      logical                                       :: esd_as_num

      !---- Data Files ----!
      character(len=512)                            :: filedat=' '
      integer                                       :: NCDat         ! Total columns for data
      integer, dimension(ncol_data_max)             :: IC_Dat        ! Which values are input - local copy

      !> Eos init
      call Init_err_Eos()

      !> Init
      dire=' '
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

      !> DataType: sets idatatype for temporary use
      i=1
      j=nlines
      idatatype=0   ! default
      call Read_Key_Str(flines, i, j, 'DataType', line,'#')
      if (len_trim(line) > 0) then
         line=U_case(trim(adjustl(line)))
         if (index(line,'MODUL') > 0) then         ! allows modulus and moduli !!
            if (index(line,'ADIA') > 0) then
               idatatype=2
            else
               idatatype=1
            end if
         end if
      end if

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

      !> Replace tabs in the format line by blank character
      call Get_Separator_pos(line,char(9),ivet,iv)
      do j=1,iv
         line(ivet(j):ivet(j))=' '
      end do

      call getword(line,dire,iv)
      if (iv <= 1) then
         err_eos=.true.
         Err_EoS_Mess="No keywords found on format line in data file!"
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
      ic_dat(1)=1     ! Set all data read in to active
      ic_dat(22)=1    ! Set all data read in to Group 1
      iorden=0

      do i=1,ncdat          ! this must loop from 1 to ncdat always so that iorden is filled correctly
         car=U_case(adjustl(dire(i+inum)))

         !> allow 'press, pressure, vol, volume, temp...' for full compatability with v5.2 old files
         if (index(car,'PRESS') /= 0)  car='P'
         if (index(car,'VOL')   /= 0)  car='V'
         if (index(car,'TEM')   /= 0)  car='T'
         if (index(car,'LIN')   /= 0)  car='A'      ! Linear - set it as 'a' in v7
         if (index(car,'SIGL')  /= 0)  car='SIGA'

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
         Err_EoS_Mess='Number of data points estimated in the data file was zero!'
         return
      end if

      !> Allocating data
      call allocate_eos_data_list(ndat,dat)

      !> Set a flag for expecting esd's as separate numbers
      esd_as_num=.false.        ! if esd's then format line implies in brackets
      if (sum(ic_dat(5:21:2)) > 0) esd_as_num=.true.

      !> reading the data starts here
      ndat=0
      nk=0
      rvet=0.0
      nlines_datum=0

      do i=nldata,nlines                ! loop over reading lines: nlines per datum
         write(unit=car,fmt='(i5)') i   ! store the line number in car for error messages
         car=adjustl(car)
         line=adjustl(flines(i))
         if (len_trim(line) <= 0 .or. line(1:1) =='#') cycle

         !> After "#" character is a comment
         npos=index(line,'#',back=.true.)
         if (npos > 0) line=line(:npos-1)

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

         !> Check for brackets indicating esd's consistent with format
         if (esd_as_num) then
            if (index(line,'(') > 0) then
               err_eos=.true.
               Err_EoS_Mess='Format says esds as numbers but esd in bracket found on  line '//trim(car)//' of data file!'
               dat%n=ndat
               return
            end if
         end if

         !> Read the numbers on this line
         nlines_datum=nlines_datum+1
         call getnum_std(line,vet,vetsd,iv)
         if (iv == 0) then
            err_eos=.true.
            Err_EoS_Mess='No values found on line '//trim(car)//' of data file!'
            dat%n=ndat
            return

         else if (iv > ncdat) then        ! This catches too many data when 1 line/datum
            err_eos=.true.
            Err_EoS_Mess='Error reading data at line '//trim(car)//': Too many data items found'
            dat%n=ndat
            return
         end if

         do j=1,iv
            kk=nk+j
            if (iorden(kk) < 1)then ! This catches too many data when >1 line/datum
                err_eos=.true.
                Err_EoS_Mess='Error reading data at line '//trim(car)//' or previous line: Too many data items found'
                dat%n=ndat
                return

            else if(iorden(kk) > ncol_data_max) then
               err_eos=.true.
               Err_EoS_Mess='Error reading data at line '//trim(car)//':  Did you put sig in format line, but have sig in ()?'
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

         !> Test for the correct number of data found
         if (nlines_datum < nl) cycle         ! not read enough lines for this datum yet
         if (nk < ncdat) then
            err_eos=.true.
            Err_EoS_Mess='Error reading data at line '//trim(car)//': Not enough data items found'
            dat%n=ndat
            return

         else if (nk > ncdat) then
            err_eos=.true.
            Err_EoS_Mess='Error reading data at line '//trim(car)//': Too many data items found'
            dat%n=ndat
            return
         end if
         nlines_datum=0     ! re-init counter for neaxt datum

         !> Writting values on Type
         ndat=ndat+1
         dat%eosd(ndat)%iuse=1          ! Active data
         dat%eosd(ndat)%igrp=1          ! Group 1
         dat%eosd(ndat)%xtype=idatatype ! Data type

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
      ! Do not test for error here because we can still set the volume

      !> Volume calculation if have cell parameters
      call Set_Volume_from_Cell(dat)

      return
   End Subroutine Read_EoS_DataFile

   !!----
   !!---- SUBROUTINE READ_EOS_FILE
   !!----
   !!---- General routine a single Eos from a file
   !!----
   !!---- Created 12/10/2015 to allow for files with more than eos
   !!---- and return the first eos in the file.
   !!----
   !!
   Subroutine Read_Eos_File(Fname,Eos)
      !---- Arguments ----!
      character(len=*),intent(in)  :: fname  ! File name
      type (EoS_Type), intent(out) :: Eos    ! EoS object

      !---- Local Variables ----!
      type (EoS_List_Type)  :: Eoslist    ! EoS list object


      call Read_Multiple_Eos_File(Fname,Eoslist)

      if (eoslist%n < 1)then
         err_eos=.true.
         Err_EoS_Mess="No EoS found in EoS file "
      else
          EoS=Eoslist%eos(1)
      end if

      return
   End Subroutine Read_Eos_File

   !!--++
   !!--++ SUBROUTINE READ_EOS_IN
   !!--++
   !!--++ PRIVATE
   !!--++ General routine to read Eos from a file
   !!--++
   !!--++ Date: 12/10/2015
   !!--++
   !!
   Subroutine Read_Eos_In(Flines,Eos)
      !---- Arguments ----!
      character(len=*),dimension(:),intent(in)   :: flines
      type(Eos_type),               intent(out)  :: eos

      !---- Local Variables ----!
      integer       :: nl, imax,ierr,idoc,nlines,i,c,j
      character(len=255)                            :: text
      character(len=10)                             :: forma
      real                                          :: val

      !> initialisation
      call init_eos_type(eos)

      nl=0
      imax=0
      ierr=0
      idoc=0                      ! local counter
      nlines=size(flines)

      do
         nl=nl+1
         if (nl > nlines) exit
         text=adjustl(flines(nl))
         if (len_trim(text) <=0) cycle   ! blank line
         if (text(1:1) == '!') cycle     ! comment line
         if (ierr /=0) exit

         c=index(text,'=')+1             !  1 place after the =
         text(1:c-1)=U_case(text(1:c-1)) ! set keyword in caps

         if (index(text,'TITLE') /= 0)then
            eos%title=trim(text(c:))

         else if(index(text,'COMMENT') /= 0)then
            idoc=idoc+1
            eos%doc(idoc)=trim(text(c:))

         else if(index(text,'MODEL') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%imodel
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Model number"

         else if(index(text,'ORDER') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%iorder
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Order number"

         else if(index(text,'THERMAL') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%itherm
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Thermal model"

         else if(index(text,'TRANS') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%itran
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Transition model"

         else if(index(text,'SHEAR') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%ishear
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Shear model"


         else if(index(text,'PSCALE') /= 0)then
            read(text(c:),'(a)',iostat=ierr)eos%pscale_name
            if (ierr /=0) Err_EoS_Mess="Error reading the Pressure Scale info"

         else if(index(text,'VSCALE') /= 0)then
            read(text(c:),'(a)',iostat=ierr)eos%vscale_name
            if (ierr /=0) Err_EoS_Mess="Error reading the Volume Scale info"

         else if(index(text,'TYPE') /= 0)then
            if(index(U_case(text),'LINEAR') /= 0) eos%linear=.true.

         else if(index(text,'PREF') /= 0)then
            read(text(c:),'(f10.0)',iostat=ierr)eos%pref
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Pressure reference"

         else if(index(text,'TREF') /= 0)then
            read(text(c:),'(f10.0)',iostat=ierr)eos%tref
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Temperature reference"

         else if(index(text,'STOICH') /= 0)then
            read(text(c:),'(f10.0)',iostat=ierr)eos%stoich
            if (ierr /=0) Err_EoS_Mess="Error reading the Stochiometry"

         else if(index(text,'DENSITY0') /= 0)then
            read(text(c:),'(f10.0)',iostat=ierr)eos%density0
            if (ierr /=0) Err_EoS_Mess="Error reading the reference density"
            if(eos%density0 < 0.00005)eos%density0=0.0_cp                   ! test against min value allowed by format

         else if(index(text,'PARAM') /= 0)then
            read(text(c:),'(i2,f12.6)',iostat=ierr)i,val
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Parameters"

            if (i > 0 .and. i <= n_eospar)then
               eos%params(i)=val
               imax=max(imax,i)   ! use this for vcv reading: allows reading of old files when n_eospar is increased
            end if

         else if(index(text,'VARIANCE') /= 0)then
            ierr=0
            forma="(   e12.5)"
            write(unit=forma(2:4),fmt="(i3)") imax
            do i=1,imax
               ! The variable format <> is a VAX extension that is taken into account by Intel Fortran
               ! however Gfortran does not support it. It is better to use a general variable for writing
               ! dynamically whatever kind of format.
               !read(unit=flines(nl+i),fmt='(<imax>e12.5)',iostat=ierr)eos%vcv(i,1:imax)
               read(unit=flines(nl+i),fmt=forma,iostat=ierr) eos%vcv(i,1:imax)
               if (ierr /= 0)exit
            end do
            if (ierr /= 0)then
               Err_EoS_Mess="Error reading the EoS Variance information"
               exit
            end if
         end if
      end do
      !>If space in comments, add eos file creation date
      j=0
      do i=1,idoc
          if(index(eos%doc(i),'Eos saved on ') > 0)j=1
      end do
      if (j == 0 .and. idoc < size(eos%doc)) then
         idoc=idoc+1
         i=1
         call Read_Key_Str(flines, i, nl, 'Current', text)
         i=index(text,':')
         if (i > 0) then
            eos%doc(idoc)='Eos saved on '//trim(text(i+1:i+30))
         end if
      end if

      !> Error during reading
      if (ierr /=0) then
         err_eos=.true.
         return
      end if

      !> Now finish setting the other eos components
      call set_eos_names(eos)
      call set_thermal_names(eos)
      call Set_Transition_Names(eos)
      call Set_Shear_Names(eos)
      call Set_EoS_Use(eos)
      call set_eos_factors(eos)           ! sets the eos factors without resetting param values

      eos%params=eos%params/eos%factor    ! rescale the values
      do i=1,n_eospar                     ! rescale the vcv matrix
         do j=1,n_eospar
            eos%vcv(i,j)=eos%vcv(i,j)/eos%factor(i)/eos%factor(j)
         end do
         eos%esd(i)=sqrt(eos%vcv(i,i))   ! set the esd's from the vcv
      end do

      call Set_Kp_Kpp_Cond(eos)           ! set default values

      !> we have to also set the refinement flags to match vcv
      do i=1,n_eospar
         if (abs(eos%vcv(i,i)) > tiny(0.0)) eos%iref(i)=1
      end do

      return
   End Subroutine Read_Eos_In

   !!----
   !!---- SUBROUTINE READ_MULTIPLE_EOS_FILE
   !!----
   !!---- General routine to read Eos from a file
   !!----
   !!---- Created 12/10/2015 to read files with more than eos
   !!----
   !!
   Subroutine Read_Multiple_Eos_File(Fname,Eoslist)
      !---- Arguments ----!
      character(len=*),     intent(in)  :: fname      ! File name
      type (EoS_List_Type), intent(out) :: Eoslist    ! EoS list object

      !---- Variables ----!
      type (EoS_Type)                               :: Eos    ! EoS  object
      character(len=255), dimension(:), allocatable :: flines
      character(len=512)                            :: filedat

      integer                                       :: nlines,neos,i
      integer,dimension(10)                         :: istart


      !> Eos init
      call Init_err_Eos()

      !> Init
      filedat=' '
      filedat=trim(fname)
      istart=0

      !> Number of lines
      call number_lines (trim(filedat),nlines)
      if (nlines <= 0) then
         err_eos=.true.
         Err_EoS_Mess="Impossible to read the EoS file "
         return
      end if

      !> Read lines
      if (allocated(flines)) deallocate(flines)
      allocate(flines(nlines))
      flines=' '
      call reading_lines(trim(filedat),nlines,flines)

      !> Find how many eos in file, and split up flines
      neos=0
      do i=1,nlines
         if (index(U_case(flines(i)),'CURRENT DATE AND TIME') > 0)then      ! This allows the date and time to be transferred into the eos%doc
            neos=neos+1
            istart(neos)=i          ! line number of TITLE
         end if
      end do
      istart(neos+1)=nlines+1             ! last line number+1
      !> If the file is from GUI there is no 'Current data and time' line
      !> Therefore set default pointers for one eos if this line not detected in file:
      if (neos == 0) then
         neos=1
         istart(1)=1
         istart(2)=nlines+1
      end if

      call Allocate_EoS_List(neos, eoslist)

      do i=1,neos
         call read_eos_in(flines(istart(i):istart(i+1)-1),eos)
         eoslist%eos(i)=eos
      end do

      return
   End Subroutine Read_Multiple_Eos_File

   !!--++
   !!--++ SUBROUTINE SET_EOS_FACTORS
   !!--++
   !!--++ PRIVATE
   !!--++ Initialize the EoS Factors without change in the parameters values
   !!--++
   !!--++ Date: 17/02/2015
   !!
   Subroutine Set_EoS_Factors(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !> Init
      eospar%factor=1.0

      select case(eospar%itherm)
         case(0)
            eospar%factor(10:19) = 1.0_cp

         case(1)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp

         case(2)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp
            eospar%factor(12)  = 1.0_cp

         case(3)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E4_cp

         case(4)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp

         case(5)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp

         case(6)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp

      end select

      select case(eospar%itran)
         case(1:3)
            eospar%factor(20:n_eospar) = 1.0_cp
            eospar%factor(24) = 1.0E3_cp         ! 1000 for aL
            eospar%factor(26) = 1.0E3_cp         ! 1000 for aH
      end select

      return
   End Subroutine Set_Eos_Factors

   !!----
   !!---- SUBROUTINE SET_EOS_NAMES
   !!----
   !!---- Set the character variables in eos_type data structures
   !!---- to match the flags already set
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Set_Eos_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar   ! EoS object

      !---- Local Variables ----!
      character(len=50),dimension(5)  :: ptext    ! local variable to hold name of pressure scale

      !> Check for valid model number. If not valid, set zero
      if (eospar%imodel < 0 .or. eospar%imodel > N_PRESS_MODELS) eospar%imodel=0

      !> Set the Eos name
      eospar%model=Pmodel_names(eospar%imodel)

      !> set the comments for parameters for volume or linear eos
      !> set the pressure scale text first
      ptext=' '
      if (len_trim(eospar%pscale_name) > 0) then
         ptext(2)='units are '//trim(eospar%pscale_name)
         ptext(4)='units are inverse '//trim(eospar%pscale_name)
         ptext(5)='units are '//trim(eospar%pscale_name)//'/K'
      else
         ptext(2)='same units as pressure data'
         ptext(4)='inverse pressure units'
         ptext(5)='units are P units/K'
      end if

      !> Set the volume/linear scale name
      if (len_trim(eospar%vscale_name) > 0) then
         ptext(1)='units are '//trim(eospar%vscale_name)
      else
         if (.not. eospar%linear) then
            ptext(1)='units as volume data'
         else
            ptext(1)='units as cell parameters'
         end if
      end if

      if (.not. eospar%linear) then
         eospar%ParName(1:5) =(/'V0   ','K0   ','Kp   ','Kpp  ','dK/dT'/)

         eospar%comment(1) = 'Reference pressure volume: '//trim(ptext(1))
         eospar%comment(2) = 'Bulk modulus: '//trim(ptext(2))
         eospar%comment(3) = 'dK/dP: dimensionless'
         eospar%comment(4) = 'd2K/dP2: '//trim(ptext(4))
         eospar%comment(5) = 'dK/dT: '//trim(ptext(5))

      else
         eospar%ParName(1:5) =(/'L0   ','M0   ','Mp   ','Mpp  ','dM/dT'/)

         eospar%comment(1) = 'Reference pressure length: '//trim(ptext(1))
         eospar%comment(2) = 'Linear modulus: '//trim(ptext(2))
         eospar%comment(3) = 'dM/dP: dimensionless'
         eospar%comment(4) = 'd2M/dP2: '//trim(ptext(4))
         eospar%comment(5) = 'dM/dT: '//trim(ptext(5))
      end if

      !> Thermal models are only set in init_EoS_thermal

      return
   End Subroutine Set_Eos_Names

   !!----
   !!---- SUBROUTINE SET_EOS_USE
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
   !!---- Date: 17/07/2015
   !!
   Subroutine Set_EoS_Use(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

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
             eospar%iuse(13)=2    ! Grunesien parameter at Pref,Tref
             eospar%iuse(14)=2    ! Grunesien q power law parameter

         case(2)             ! Fei
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10:12)=1 ! alpha terms
             eospar%iuse(13)=2    ! Grunesien parameter at Pref,Tref
             eospar%iuse(14)=2    ! Grunesien q power law parameter

         case(3)             ! HP 1998
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10:11)=1 ! alpha terms
             eospar%iuse(13)=2
             eospar%iuse(14)=2    ! Grunesien q power law parameter

         case(4)             ! Holland-Powell thermal expansion, in Kroll form
             if (eospar%iuse(3) == 0)eospar%iuse(3)=2     ! require Kprime_zero but not stable in refinement if no P data (added 27/01/2014 RJA)
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10)=1    ! alpha at Tref
             eospar%iuse(11)=2    ! Einstein T should be reported but cannot be refined
             eospar%iuse(13)=2    ! Grunesien parameter at Pref,Tref
             eospar%iuse(14)=2    ! Grunesien q power law parameter

         case(5)             ! Salje
             if (eospar%imodel /=0) eospar%iuse(5)=1     ! allow dK/dT
             eospar%iuse(10:11)=1
             eospar%iuse(13)=2    ! Grunesien parameter at Pref,Tref
             eospar%iuse(14)=2    ! Grunesien q power law parameter

         case(6)             ! Thermal pressure in H&P form (no dK/dT): requires a eos model as well
             if (eospar%imodel==0)eospar%iuse(2:3)=2   ! K, Kp refinement is not stable without a pressure model
             eospar%iuse(5)=0     ! No dK/dT parameter:
             eospar%iuse(10)=1    ! alpha at Tref
             eospar%iuse(11)=2    ! Einstein T should be reported but cannot be refined
             eospar%iuse(13)=2    ! Grunesien parameter at Pref,Tref
             eospar%iuse(14)=2    ! Grunesien q power law parameter
      end select

      !> Phase transition model
      select case(eospar%itran)
         case(1,2)     ! Landau PV or TV
            eospar%iuse(20)=2           !settable, no refine: sense of transition,
            eospar%iuse(21)=1           !settable, allow refine: Ptr or Ttr
            eospar%iuse(24)=1           !settable, allow refine: aL
            eospar%iuse(25)=1           !settable, allow refine: betaL
            eospar%iuse(26)=1           !settable, allow refine: aH
            eospar%iuse(27)=1           !settable, allow refine: betaH

         case(3)     ! Landau PVT
            eospar%iuse(20:22)=2        !settable, no refine: sense of transition, T(tr), dT(Tr)/dP
            eospar%iuse(24)=1           !settable, allow refine: aL,
            eospar%iuse(23)=2           !settable, fixed d2Tr/dP2
            eospar%iuse(25:27)=1        !settable, allow refine:  betaL,aH,betaH
      end select

      !> Phase transition model
      select case(eospar%ishear)
         case(0)
            eospar%iuse(30)=3            ! No model, G0 set very large
         case(1)
            eospar%iuse(30:34)=2         ! Polynomial model: settable, not refineable
      end select

      !> Set the refine flags to be consistent with the use flags
      do i=1,n_eospar
         if (eospar%iuse(i) /=1) eospar%iref(i)=0
      end do

      return
   End Subroutine Set_Eos_Use

   !!----
   !!---- SUBROUTINE SET_KP_KPP_COND
   !!----
   !!---- Fix Kp and Kpp values from Model and Order of EoSpar
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Set_Kp_Kpp_Cond(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

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
      !call Set_Eos_Use(eospar)

      return
   End Subroutine Set_Kp_Kpp_Cond

   !!--++
   !!--++ SUBROUTINE SET_SHEAR_NAMES
   !!--++
   !!--++ PRIVATE
   !!--++ Set the character variables in eos_type data structures for SHEAR EoS
   !!--++
   !!--++
   !!--++ Date: 11/07/2016
   !!
   Subroutine Set_Shear_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

      !---- Local Variables ----!
      integer       :: n

      !> Check for valid model number. If not valid, set zero
      if(eospar%ishear < 0 .or. eospar%ishear > N_SHEAR_MODELS) eospar%ishear=0

      !> Set the Eos name
      eospar%smodel=shearmodel_names(eospar%ishear)

      !> Set upper limit to thermal parameter numbers
      n=34
      if(n >n_eospar)n=n_eospar

      select case(eospar%ishear)
         case(0)
            eospar%parname(30:n) = ' '
            eospar%comment(30:n) = ' '

         case(1)
            eospar%parname(30)='G0'
            eospar%parname(31)='dG/dP'
            eospar%parname(32)='d2G/'
            eospar%parname(33)='d3G/'
            eospar%parname(34)='dG/dT'
            eospar%comment(30)='Shear modulus at Pref, Tref, in pressure units'
            eospar%comment(31)='Pressure derivative of shear modulus: no units'
            eospar%comment(32)='2nd Pressure derivative of shear modulus: P^-1'
            eospar%comment(33)='3rd Pressure derivative of shear modulus: P^-2'
            eospar%comment(34)='Temperature derivative of shear modulus'
      end select

      return
   End Subroutine Set_Shear_Names

   !!--++
   !!--++ SUBROUTINE SET_THERMAL_NAMES
   !!--++
   !!--++ PRIVATE
   !!--++ Set the character variables in eos_type data structures for thermal EoS
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Set_Thermal_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

      !---- Local Variables ----!
      integer       :: n

      !> Check for valid model number. If not valid, set zero
      if(eospar%itherm < 0 .or. eospar%itherm > N_THERM_MODELS) eospar%itherm=0

      !> Set the Eos name
      eospar%tmodel=Tmodel_names(eospar%itherm)

      !> Set upper limit to thermal parameter numbers
      n=19
      if(n >n_eospar)n=n_eospar

      select case(eospar%itherm)
         case(0)
            eospar%parname(10:n) = ' '
            eospar%comment(10:n) = ' '

         case(1)
            eospar%parname(10:11) = (/'alph0','alph1'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Linear term thermal expansion x10^8 K^-2'

         case(2)
            eospar%parname(10:12) = (/'alph0','alph1','alph2'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Linear term thermal expansion x10^8 K^-2'
            eospar%comment(12) = '1/T^2 term thermal expansion, K'

         case(3)
            eospar%parname(10:11) = (/'alph0','alph1'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Sqrt term of thermal expansion x10^4 K^-1/2'

         case(4)
            eospar%parname(10:11) = (/'alph0','Th_E '/)
            eospar%comment(10) = 'Constant of thermal expansion at Tref x10^5 K^-1'
            eospar%comment(11) = 'Einstein temperature in K'


         case(5)
            eospar%parname(10:11) = (/'p1   ','T_sat'/)
            eospar%comment(10) = 'Approx 3x highT thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Saturation temperature in K'


         case(6)
            eospar%parname(10:11) = (/'alph0','Th_E '/)
            eospar%comment(10) = 'Constant of thermal expansion at Tref x10^5 K^-1'
            eospar%comment(11) = 'Einstein temperature in K'

      end select

      !> Common terms for all thermal
      eospar%parname(13) = 'Gamm0'
      eospar%comment(13) = 'Gruneisen parameter at Tref,Pref'
      eospar%parname(14) = 'q    '
      eospar%comment(14) = 'Gruneisen power law in V/V0'

      return
   End Subroutine Set_Thermal_Names

   !!--++
   !!--++ SUBROUTINE SET_TRANSITION_NAMES
   !!--++
   !!--++ PRIVATE
   !!--++ Set the character variables in eos_type data structures for Transition EoS
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Set_Transition_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

      !---- Local Variables ----!
      integer :: n

      !> Check for valid model number. If not valid, set zero
      if (eospar%itran < 0 .or. eospar%itran > N_TRANS_MODELS) eospar%itran=0

      !> Set the model name
      eospar%tranmodel=Tranmodel_names(eospar%itran)

      !> Set upper limit to parameter numbers
      n=29
      if (n >n_eospar)n=n_eospar

      select case(eospar%itran)
         case(0)
            eospar%parname(20:n) = ' '
            eospar%comment(20:n) = ' '

         case(1)       ! Landau power law P-V
            eospar%parname(20:27) = (/'High ','Ptr  ','     ','     ','aL   ','betaL','aH   ','betaH'/)
            eospar%comment(20) = 'Indicator = +1 if high P phase is high sym phase'
            eospar%comment(21) = 'Transition pressure'
            eospar%comment(22) = 'dTr/dP: slope of transition boundary'
            eospar%comment(23) = ''
            eospar%comment(24) = 'Scaling parameter, low phase x10^3'
            eospar%comment(25) = 'Power law term, low phase'
            eospar%comment(26) = 'Scaling parameter, high phase x10^3'
            eospar%comment(27) = 'Power law term, high phase'

         case(2)       ! Landau power law V-T
            eospar%parname(20:27) = (/'High ','Ttr  ','     ','     ','aL   ','betaL','aH   ','betaH'/)
            eospar%comment(20) = 'Indicator = +1 if high T phase is high sym phase'
            eospar%comment(21) = 'Transition temperature'
            eospar%comment(22) = ''
            eospar%comment(23) = ''
            eospar%comment(24) = 'Scaling parameter, low phase x10^3'
            eospar%comment(25) = 'Power law term, low phase'
            eospar%comment(26) = 'Scaling parameter, high phase x10^3'
            eospar%comment(27) = 'Power law term, high phase'

         case(3)       ! Landau power law PVT
            eospar%parname(20:27) = (/'High ','Ttr  ','Tr/dP','dTr2 ','aL   ','betaL','aH   ','betaH'/)
            eospar%comment(20) = 'Indicator = +1 if high T phase is high sym phase'
            eospar%comment(21) = 'Transition temperature'
            eospar%comment(22) = 'dTr/dP: slope of transition boundary'
            eospar%comment(23) = 'd2Tr/dP2: curvature of phase boundary'
            eospar%comment(24) = 'Scaling parameter, low phase x10^3'
            eospar%comment(25) = 'Power law term, low phase'
            eospar%comment(26) = 'Scaling parameter, high phase x10^3'
            eospar%comment(27) = 'Power law term, high phase'

      end select

      return
   End Subroutine Set_Transition_Names

   !!--++
   !!--++ SUBROUTINE SET_VOLUME_FROM_CELL
   !!--++
   !!--++ PRIVATE
   !!--++ Sets V and esd(V) from cell parameter data for all data items in dat
   !!--++ If V is present in first data item, no esd is calculated
   !!--++
   !!--++ Update: 17/07/2015 16:47:03
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
   !!--++ SUBROUTINE TRANSFORM_ESD
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
      real(kind=cp),dimension(N_EOSPar,N_EOSPar)  :: d                  ! cross derivatives of parameters
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
         if (eospar%iuse(j) == 1) then
            do k=-2,2,1
               eost=eospar                                                     !reset eos params
               if (abs(eospar%params(j)) > tiny(0.0_cp)) then
                  shift=factor*eospar%params(j)
                  if (j == 10) shift=10.0*shift                    ! alpha0
               else
                  shift=1.0_cp/eospar%factor(j)                    ! shift to a parameter with zero value
                  if (j == 5) shift=0.002                          ! dK/dT has print factor 1, but typical value -.02
               end if

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
            if (eospar%iref(i) == 0) cycle ! because vcv(i,j) will be 0
            do j=1,N_EOSPar
               if (eospar%iref(j) == 0) cycle ! because vcv(i,j) will be 0
               esd(k)=esd(k)+eospar%vcv(i,j)*d(k,i)*d(k,j)
            end do
         end do
      end do

      !> Final: extra trap June 2016 in case round off leaves esd(i)^2 < 0
      do i=1,n_eospar
          if(esd(i) > tiny(0._cp))then
            esd(i)=sqrt(esd(i))
          else
            esd(i)=0._cp
          endif
      enddo


      return
   End Subroutine Transform_Esd

   !!--++
   !!--++ SUBROUTINE VEC_TO_EOS
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

      !> Copy vec to eosparams as 1 to 1 as default
      eospar%params=vec

      if (eospar%linear) then
         eospar%params(1)=vec(1)**(1.0_cp/3.0_cp)
         eospar%params(2:5)=vec(2:5)*3.0_cp
         select case(eospar%itherm)                 ! thermal expansion terms
            case(1,2,3)
               eospar%params(10:12)=vec(10:12)/3.0_cp

            case(4,5,6)
               eospar%params(10)=vec(10)/3.0_cp
               eospar%params(11)=vec(11)
         end select

         select case(eospar%itran)          ! phase transition
            case(1,2)      ! V or T only models
               eospar%params(24)=vec(24)/3.0_cp
               eospar%params(26)=vec(26)/3.0_cp

            case(3)        ! PVT model
               eospar%params(23:24)=vec(23:24)/3.0_cp
               eospar%params(26)=vec(26)/3.0_cp
         end select
      end if

      return
   End Subroutine Vec_to_EoS

   !!----
   !!---- SUBROUTINE WRITE_DATA_CONLEV
   !!----
   !!---- Writes out confidence ellipse data to a file
   !!----
   !!---- Date: 05/12/2015
   !!
   Subroutine Write_Data_Conlev(xyy,n,iout)
      !---- Arguments ----!
      real(kind=cp),dimension(:,:),intent(in)   :: xyy  ! output points for plotting
      integer,                     intent(in)   :: n    ! Number of points
      integer, optional,           intent(in)   :: iout ! Iunit output

      !---- Local Variables ----!
      integer :: lun,i

      !> Check
      if (n < 1)return

      !> Unit to print the data
      lun=6
      if (present(iout)) lun=iout

      do i=1,n
         write(lun,'(3x,a,3x,a,3x,a)')trim(rformat(xyy(1,i),10)),trim(rformat(xyy(2,i),10)),trim(rformat(xyy(3,i),10))
      end do

      return
   End Subroutine Write_Data_Conlev

   !!----
   !!---- SUBROUTINE WRITE_EOS_DATAFILE
   !!----
   !!---- General routine to Write Data in a Lun iunit
   !!----
   !!---- Update: 17/07/2015
   !!
   Subroutine Write_EoS_DataFile(dat,lun)
      !---- Arguments ----!
      type (eos_data_list_type), intent(in) :: dat  ! data structure
      integer,                   intent(in) :: lun  ! Unit to write the information

      !---- Variables ----!
      integer              :: ierr,i,j,k
      character(len=128)   :: text

      !> set up the labels in order
      integer, parameter           :: ini=4,iend=21
      character(len=5),dimension(ini:iend) :: lab=(/'T    ','sigT ','P    ','sigP ','V    ','sigV ',  &
                                                    'A    ','sigA ','B    ','sigB ','C    ','sigC ',  &
                                                    'ALPHA','sigAL','BETA ','sigBE','GAMMA','sigGA'/)

      !>
      !> assume that unit is connected and open.
      !>

      !> Write header info
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'TITLE ',trim(dat%title)
      write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      write(unit=lun,fmt='(a)',iostat=ierr)    '#  Data file written by CFML eos module'
      write(unit=lun,fmt='(a)',iostat=ierr)    '#'

      !> Crystal system
      if (len_trim(dat%system) > 0)then
         write(unit=lun,fmt='(a,a)',iostat=ierr)  'SYSTEM ',trim(dat%system)
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if

      !> Original Tscale of data is not known, write out Tscale K if T data present
      if (dat%ic_dat(4) == 1)then
         write(unit=lun,fmt='(a)',iostat=ierr)  'TSCALE  K'
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if

      !> build format line
      text='FORMAT 1'
      do i=ini,iend
         if (dat%ic_dat(i) ==1) text=trim(text)//' '//lab(i)
      end do
      write(unit=lun,fmt='(a)',iostat=ierr)    trim(text)

      !> write the data: all data written out
      do i=1,dat%n        ! loop over data points
         if (dat%ic_dat(4) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%t
         if (dat%ic_dat(5) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%sigt
         if (dat%ic_dat(6) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%p
         if (dat%ic_dat(7) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%sigp
         if (dat%ic_dat(8) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%v
         if (dat%ic_dat(9) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%sigv

         !> small loop over cell edges
         do j=1,3
            k=8+2*j
            if (dat%ic_dat(k) == 1)   write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%cell(j)
            if (dat%ic_dat(k+1) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%sigc(j)
         end do

         !> small loop over cell angles
         do j=1,3
            k=14+2*j
            if (dat%ic_dat(k) == 1)   write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%ang(j)
            if (dat%ic_dat(k+1) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%siga(j)
         end do
         write(unit=lun,fmt='(1x)',iostat=ierr)          ! forces new line
      end do

      return
   End Subroutine Write_Eos_Datafile

   !!----
   !!---- SUBROUTINE WRITE_EOS_FILE
   !!----
   !!---- General routine to Write EoS in a Lun iunit
   !!----
   !!---- Update: 17/07/2015
   !!
   Subroutine Write_Eos_File(Eos,Lun)
      !---- Arguments ----!
      type (EoS_Type),intent(in)   :: Eos ! EoS object
      integer,intent(in)           :: lun ! Unit

      !---- Variables ----!
      character(len=12)            :: stext
      character(len=512)           :: text
      integer                      :: ierr,i,j

      !>
      !> assume that unit is connected and open.
      !>

      !> Header info written to file from calling program
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
      write(unit=lun,fmt='(a)',iostat=ierr) ' EosFit parameter file'
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
      write(unit=lun,fmt='(a)',iostat=ierr) ' This is a fixed-format file. If you edit this file, make sure you do not move '// &
                                            'anything!'
      write(unit=lun,fmt='(a)',iostat=ierr) ' It is safer to change the parameters by loading the file to EosFit, and saving '// &
                                            'the file after making changes'
      write(unit=lun,fmt='(a)',iostat=ierr) ' _______________________________________________________________________________'// &
                                            '_____________________________'
      write(unit=lun,fmt='(a)',iostat=ierr) ' '


      !> title and comments
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Title =',trim(eos%title)
      do i=1,10
         if (len_trim(eos%doc(i)) > 0)then
            write(unit=lun,fmt='(a)') 'Comment ='//trim(eos%doc(i))
         end if
      end do
      write(unit=lun,fmt='(a)',iostat=ierr) ' '

      !> eos type
      text=',  ('//trim(eos%model)//')'
      if (eos%imodel == 0) text=',  (none)'
      write(unit=lun,fmt='(a,i3,a)',iostat=ierr) 'Model =',eos%imodel,text
      write(unit=lun,fmt='(a,i3)',iostat=ierr) 'Order =',eos%iorder

      text=',  ('//trim(eos%tmodel)//')'
      if (eos%itherm == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Thermal =',eos%itherm,text

      text=',  ('//trim(eos%tranmodel)//')'
      if (eos%itran == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Trans =',eos%itran,text

      text=',  ('//trim(eos%smodel)//')'
      if (eos%ishear == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Shear =',eos%ishear,text

      if (eos%linear)then
         write(unit=lun,fmt='(a)',iostat=ierr) 'Type = Linear'
      else
         write(unit=lun,fmt='(a)',iostat=ierr) 'Type = Volume'
      end if
      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Pref =',eos%pref
      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Tref =',eos%tref
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Pscale =',trim(eos%pscale_name)
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Vscale =',trim(eos%vscale_name)
      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Stoich =',eos%stoich
      if (eos%density0 > tiny(0.)) then
         write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Density0 =',eos%density0
      end if

      write(unit=lun,fmt='(a)',iostat=ierr) ' '

      !> Eos parameters
      do i=1,n_eospar
         if (eos%iuse(i) == 0)then
            write(unit=lun,fmt='(a,i2,f12.6,5a)')'Param =',i,eos%params(i)*eos%factor(i)
         else
            write(unit=lun,fmt='(a,i2,f12.6,5a)')'Param =',i,eos%params(i)*eos%factor(i),'     (',eos%parname(i),',  ',&
                 trim(eos%comment(i)),')'
         end if
      end do

      !> VCV: stored as scaled values for precision
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
      write(unit=lun,fmt='(a)',iostat=ierr) '  Variance-Covariance matrix='
      do i=1,n_eospar
         text=''
         do j=1,n_eospar
            if (abs(eos%vcv(i,j)) > tiny(0.0))then
               write(stext,fmt='(e12.5)')eos%vcv(i,j )*eos%factor(i)*eos%factor(j)
               text=trim(text)//stext
            else              !   123456789012
               text=trim(text)//' 0.00000E+00'
            end if
         end do
         write(unit=lun,fmt='(a)',iostat=ierr)trim(text)
      end do
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
      write(unit=lun,fmt='(a)',iostat=ierr) ' '

      return
   End Subroutine Write_Eos_File

   !!----
   !!---- SUBROUTINE WRITE_EOSCAL
   !!----
   !!----   Subroutine to write the calculated parameters of an eos to file at a series of PT points
   !!----   NO  header info is printed here
   !!----   Therefore the program header and write_info_eos have to be called first before calling this routine
   !!----   Then write_eoscal_header is called from here
   !!----
   !!----   Change: 06/10/2015 to make write_eoscal_header private, and change name from write_eoscal_file
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Write_Eoscal(Pmin,Pmax,Pstep,Tmin,Tmax,Tstep,Tscale_In,Eos,Lun,Nprint)
      !---- Arguments ----!
      real(kind=cp),    intent(in)  ::  pmin, pmax, pstep   !P to calculate properties
      real(kind=cp),    intent(in)  ::  tmin,tmax,tstep     !T to calculate properties
      character(len=*), intent(in)  ::  tscale_in           ! Name of the Tscale for output, either C or K
                                                            ! If Pstep or Tstep  < tiny(0.0) then only Pmin (or Tmin) calculated
      type(EoS_Type),   intent(in)  ::  eos                 ! Eos
      integer,          intent(in)  :: lun                  ! logical unit for printing
      integer,          intent(out) :: nprint               ! Number of data printed

      !---- Local variable ----!
      real(kind=cp)           :: p,t      ! The p and T of each calculation
      real(kind=cp)           :: pst,tst  ! local copy of tstep and pstep
      character(len=255)      :: text     ! local text variable
      character(len=1)        :: tscale   ! local name of tscale
      logical                 :: loop_p   ! loop indicator .true. for inner loop of calcs over P
      integer,dimension(14)   :: ip=(/6,6,9,8,6,5,  5, 9, 7, 7,    5,  7, 7,7/) ! format for output
      integer                 :: i

      real(kind=cp),dimension(6) :: parvals(7)
      real(kind=cp),dimension(6) :: esd
      real(kind=cp),dimension(14):: parout,esdout
      real(kind=cp)              :: v0,fp,fs

      !> init
      nprint=0    ! output counter

      !> Tscale for output: C or K
      if (len_trim(tscale_in) == 0)then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      !> Write file header
      call write_eoscal_header(eos,lun,tscale)

      !> copy Pstep/Tstep
      tst=tstep
      pst=pstep

      !> set up loop control variables
      if (abs(pst) > tiny(0.))then       ! inner loop over P
         loop_p=.true.
         if (abs(tst) < tiny(0.))then   ! no outerloop
            tst=10.*max((tmax-tmin),1.0)       ! set tstep big enough to stop loop
         end if
      else
         loop_p=.false.                  ! inner loop over T
         if (abs(pst) < tiny(0.))then   ! no outerloop
            pst=10.*max((pmax-pmin),1.0)       ! set pstep big enough to stop loop
         end if
      end if

      !> Initialise loop variables
      p=pmin
      t=tmin

      !> Start of outer loop
      outer: do
         !> reset inner loop variable to start value
         if (loop_p)then
            p=pmin
         else
            t=tmin
         end if

         inner: do
            call init_err_eos()
            esd=0.
            esdout=0.
            parout=0.

            !> Now do the calculations at P,T
            call EoS_Cal(P,T,eos,Parvals)    ! GET V,K ETC
            if (sum(eos%vcv) > tiny(0.0)) CALL eos_cal_esd(P,T,eos,esd)

            !> build ouput value array
            V0=Get_Volume(0.0,T,Eos)

            parout(1)=p
            parout(2)=t
            if (tscale =='C')parout(2)=parout(2)-273.16
            parout(3)=parvals(1)*eos%factor(1)      ! v
            esdout(3)=esd(1)*eos%factor(1)
            parout(4)=parvals(1)/v0                 ! v/V0 at this T
            esdout(4)=esdout(3)/v0

            !> convert  V,K,Kp,Kpp to output values
            do i=2,4
               parout(i+3)=parvals(i)*eos%factor(i)
               esdout(i+3)=esd(i)*eos%factor(i)
            end do

            !>deal with f-F:
            if (abs(p) < tiny(0.))then
               CALL FFCAL_EOS(P,T,EOS,FP,FS)      ! because F not defined numerically at P=0
               parout(9)=FP
               esdout(9)=0.
            else
               Call FfCal_Dat_Esd(parvals(1),esd(1),V0,0.0_cp,P,0.0_cp,Eos, &          ! only esd input is esd(V) at this P
                    parout(9),esdout(9),parout(8),esdout(8))
            end if

            !> dK/dT
            parout(10)=parvals(5)*eos%factor(5)
            esdout(10)=esd(5)*eos%factor(5)

            !> handle alpha
            parout(11)=parvals(6)*eos%factor(10)
            esdout(11)=esd(6)*eos%factor(10)

            !> spon strain
            if (eos%itran > 0) parout(12)=Get_Transition_Strain(P,T,Eos)

            !> density
            if (eos%density0 > tiny(0.0)) then
               parout(13)=eos%density0*eos%params(1)/parvals(1)
               parout(14)=parout(13)*esd(1)/parvals(1)
            end if

            !> output this datum: dynamic formatting to text string
            nprint=nprint+1

            !>init
            text=''

            !> pressure (no esd)
            text=trim(rformat(parout(1),ip(1)))

            !> T value (no esd)
            text=trim(text)//'  '//trim(rformat((parout(2)),ip(2)))

            !> other params
            do i=3,11
                text=trim(text)//'  '//trim(rformat(parout(i),ip(i)))//' '//trim(rformat(esdout(i),ip(i)))
            end do

            if (eos%itran > 0) text=trim(text)//'  '//trim(rformat(parout(12),ip(12)))    ! spontaneous strain

            if (eos%density0 > tiny(0.0)) &
               text=trim(text)//'  '//trim(rformat(parout(13),ip(13)))//' '//trim(rformat(parout(14),ip(14)))

            write(lun,'(a)')trim(text)      ! This way we get to see the calculated values even if error
            if (err_eos) then
               text=text(1:14)//':   '//trim(err_eos_mess)
               write(lun,'(a)')trim(text)
            end if

            !> Now increment inner loop variable and test for completion
            if (loop_p) then
               p=p+pst                           ! inner loop over p
               if (p > pmax+0.99_cp*pst)exit inner
            else
               t=t+tst                           ! inner loop over t
               if (t > tmax+0.99_cp*tst)exit inner
            end if
         end do inner

         write(lun,'(/)')        ! blank line to help some plotting programs

         !> Now increment outer loop variable and test for completion
         if (loop_p)then
            t=t+tst                           ! outer loop over t
            if (t > tmax+0.99_cp*tst)exit outer
         else
            p=p+pst                           ! outer loop over p
            if (p > pmax+0.99_cp*pst)exit outer
         end if
      end do outer

      return
   End Subroutine Write_Eoscal

   !!--++
   !!--++ SUBROUTINE WRITE_EOSCAL_HEADER
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Write_Eoscal_Header(Eos,Lun,Tscale_In)
      !---- Arguments ----!
      type(EoS_Type),intent(in)   :: eos         ! Eos information
      integer,       intent(in)   :: lun         ! logical unit for printing
      character(len=*),intent(in) :: tscale_in   ! Scale for Temp

      !---- Local Variables ----!
      character(len=1)     :: tscale
      character(len=255)   :: head     ! local text variable for column headers

      !> Warning for Tait or Murnaghan
      write(lun,'(//)')
      if (eos%itherm /= 0)then
         write(lun,'("  Note that values of alpha are multiplied by a factor of ",f5.1,"x10^5"//)')eos%factor(10)/1.0E5
      end if

      if (eos%imodel == 1 .or. eos%imodel == 5) then
         write(lun,'("  Do not forget: Normalised Pressure and strain not defined for ",a," Eos")')trim(eos%model)
      else if(eos%itherm /= 0)then
         write(lun,'("  Normalised Pressure (NP) and finite strain (f) are defined relative to V at P=0 and same T")')
      else
         write(lun,'("  Normalised Pressure is NP and finite strain is f")')
      end if

      !> tscale for output: C or K
      if (len_trim(tscale_in) == 0)then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      !> create column header
      if (eos%linear) then
         write(head,'(" Press   Temp",a1,"     Length  esdL       L/L0T  esd(L/L0T)  M      esdM  Mprime ", &
             &"esdMp   Mpp  esdMpp  f         esdf       NP      esdNP   dM/dT    esddK   alpha  esda")' ) Tscale
      else
         write(head,'(" Press   Temp",a1,"     Volume  esdV       V/V0T  esd(V/V0T)  K      esdK  Kprime ", &
             &"esdKp   Kpp  esdKpp  f         esdf       NP      esdNP   dK/dT    esddK   alpha  esda")' ) Tscale
      end if
      if (eos%itran > 0)head=trim(head)//' spstrain'

      if (eos%density0 > tiny(0.0)) head=trim(head)//'  density  esdden'

      !> Write header
      write(lun,'(/a)')trim(head)

      return
   End Subroutine Write_Eoscal_Header

   !!----
   !!---- SUBROUTINE Write_Info_Conlev
   !!----
   !!---- Writes out header info specific to a confidence ellipse
   !!----
   !!---- Date: 05/12/2015
   !!
   Subroutine Write_Info_Conlev(Eos,ix,iy,isig,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos
      integer,           intent(in) :: ix,iy,isig  ! parameter numbers for the plot, and sigma level
      integer, optional, intent(in) :: iout

      !---- Local Variables ----!
      integer :: lun


      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      !> Info
      write(lun,'(//"Coordinates for confidence ellipse at confidence level of ",f6.2,"%")')delchi_levels(isig)

      write(lun,'(3x,"X axis as ",a," with variance = ",a)')trim(eos%parname(ix)),trim(rformat(eos%vcv(ix,ix),10))
      write(lun,'(3x,"Y axis as ",a," with variance = ",a)')trim(eos%parname(iy)),trim(rformat(eos%vcv(iy,iy),10))
      write(lun,'(3x,"The covariance of the two variables = ",a)')trim(rformat(eos%vcv(ix,iy),10))

      write(lun,'(//"Data points for plotting: first column is X", &
                " with two columns with the two Y values at that X"//)')


      return
   End Subroutine Write_Info_Conlev

   !!----
   !!---- SUBROUTINE WRITE_INFO_EOS
   !!----
   !!---- Subroutine that print information on iout unit
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Write_Info_Eos(Eospar,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eospar  ! EoS object
      integer, optional, intent(in) :: iout    ! Logical unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun


      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      !> Header / Title
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  EOS Information'
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '    Title: '//trim(eospar%title)

      !> Doc Information
      do i=1,10
          if (len_trim(eospar%doc(i)) > 0) then
             write(unit=lun,fmt='(a)') '  Comment: '//trim(eospar%doc(i))
          end if
      end do
      write(unit=lun,fmt='(a)') ' '
      if (eospar%imodel > 0) then
         if (len_trim(eospar%Pscale_name) > 0)write(unit=lun,fmt='(a)') '   Pscale: '//trim(eospar%Pscale_name)
         if (len_trim(eospar%Vscale_name) > 0)write(unit=lun,fmt='(a)') '   Vscale: '//trim(eospar%Vscale_name)
         write(unit=lun,fmt='(a,t27,f8.3)') '   Stoichiometry: ',eospar%stoich


         !> Reference Density
         if (eospar%density0 > tiny(0.0)) then
            write(unit=lun,fmt='(a,t27,f8.3)') '   Reference density: ',eospar%density0
         end if

         write(unit=lun,fmt='(a)') ' '
         write(unit=lun,fmt='(a)') '  Compressibility'
         write(unit=lun,fmt='(a)') '-------------------'
         write(unit=lun,fmt='(a)') ' '
         write(unit=lun,fmt='(a)') '   Model: '//trim(eospar%model)
         if (eospar%imodel > 0) then                ! no output if no p eos: all done in write_info_eos_thermal
            write(unit=lun,fmt='(a,i2)') '   Order: ',eospar%iorder
            if (eospar%linear) then
               write(unit=lun,fmt='(a)') '   Class: Linear'
            else
               write(unit=lun,fmt='(a)') '   Class: Volume'
            end if

            !> Pressure Parameters
            write(unit=lun,fmt='(a,t27,f8.3)') '   Pressure of reference: ',eospar%pref
            write(unit=lun,fmt='(a)') ' '

            do i=1,4
               if (eospar%iuse(i) /= 0) then
                  call setnum_std(eospar%params(i)*eospar%factor(i),eospar%esd(i)*eospar%factor(i),line)     ! include scaling
                  string=' '
                  if (eospar%iuse(i) == 2)string=' [IMPLIED VALUE]'
                  write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                       trim(eospar%parname(i)),trim(line),trim(eospar%comment(i))//trim(string)
               end if
            end do
            write(unit=lun,fmt='(a)') ' '
         end if
      end if

      !> Thermal EOS
      if (eospar%itherm > 0) call write_info_eos_thermal(eospar,lun)

      !> Transition
      if (eospar%itran > 0) call write_info_eos_transition(eospar,lun)

      !> Shear
      if (eospar%ishear > 0) call write_info_eos_shear(eospar,lun)

      !> End
      write(unit=lun,fmt='(a)') ' '

      return
   End Subroutine Write_Info_Eos

   !!--++
   !!--++ SUBROUTINE WRITE_INFO_EOS_SHEAR
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 11/07/2016
   !!
   Subroutine Write_Info_Eos_Shear(Eos,Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos
      integer, optional, intent(in) :: iout

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Check
      if(Eos%ishear ==0) return

      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Shear modulus Information'
      write(unit=lun,fmt='(a)') '---------------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(eos%smodel)

      do i=30,39
         if (eos%iuse(i) /= 0) then
            call setnum_std(eos%params(i)*eos%factor(i),eos%esd(i)*eos%factor(i),line)     ! include scaling
            string=' '
            if (eos%iuse(i) == 2)string=' [IMPLIED VALUE]'
            write(unit=lun,fmt='(3x,a5,'': '',a,T30,'':'',a)') &
                 trim(eos%parname(i)),trim(line),trim(eos%comment(i))//trim(string)
         end if
      end do
      write(unit=lun,fmt='(a)') ' '

      return
   End Subroutine Write_Info_Eos_Shear

   !!--++
   !!--++ SUBROUTINE WRITE_INFO_EOS_THERMAL
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Write_Info_Eos_Thermal(Eospar,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eospar  ! EoS object
      integer, optional, intent(in) :: iout    ! Logical unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun,is

      !> Init
      lun=6
      if (present(iout)) lun=iout

      is=5                             ! If pmodel present, it was already reported
      if (eospar%imodel == 0) is=1     ! if no pmodel report all params here


      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Thermal Expansion'
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(eospar%tmodel)
      write(unit=lun,fmt='(a)') ' '

      write(unit=lun,fmt='(a,f8.2,a)') '   Temperature of reference: ',eospar%tref,' K'
      write(unit=lun,fmt='(a)') ' '
      do i=is,19
         if (eospar%iuse(i) /= 0) then
            call setnum_std(eospar%params(i)*eospar%factor(i),eospar%esd(i)*eospar%factor(i),line)     ! include scaling
            string=' '
            if (eospar%iuse(i) > 1) string=' [FIXED VALUE]'
            write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                  trim(eospar%parname(i)),trim(line),trim(eospar%comment(i))//trim(string)
         end if
      end do

      return
   End Subroutine Write_Info_Eos_Thermal

   !!--++
   !!--++ SUBROUTINE WRITE_INFO_EOS_TRANSITION
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Write_Info_Eos_Transition(Eospar,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eospar  ! EoS object
      integer, optional, intent(in) :: iout    ! Logical Unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Check
      if(eospar%itran == 0)return

      !> Init
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Phase Transition '
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(eospar%tranmodel)
      write(unit=lun,fmt='(a)') ' '

      do i=20,29
         if (eospar%iuse(i) /= 0) then
             call setnum_std(eospar%params(i)*eospar%factor(i),eospar%esd(i)*eospar%factor(i),line)     ! include scaling
             string=' '
             if (eospar%iuse(i) > 1) string=' [FIXED VALUE]'
             write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                  trim(eospar%parname(i)),trim(line),trim(eospar%comment(i))//trim(string)
         end if
      end do

      return
   End Subroutine Write_Info_Eos_Transition

End Module CFML_EoS