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
!!----                          Universita di Pavia, Pavia, ITALY
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL) (CrysFML)
!!----          Javier Gonzalez-Platas  (ULL) (CrysFML)
!!----          Ross John Angel         (Pavia)   (EoS)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----               Ross John Angel    (Universita di Pavia, Italy)
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
!!----    Update: 16/03/2017
!!----
!!----    ...   2013: First version of the Module
!!----    23/02/2013: Incorporates bug fix to get_volume (JGP)
!!----    ...   2014: Validation for non-thermal EOS and modifications (RJA)
!!----    ...   2015: Addition of code to handle Landau-type phase transitions (RJA)
!!----    ...   2015: Modifications for 2003 compatibility and Eos development (RJA,  JGP)
!!----    ...   2016: Modifications for curved phase boundaries and fitting moduli (RJA)
!!----    ...   2016: Modifications to add PVT table to EoS structure as alternative to EoS parameters (RJA)
!!----    ...   2017: Split write_eoscal to improve error handling, new thermalP eos (RJA) 
!!----    ...   2018: New routines: physical_check, get_kp, added pscale,vscale,lscale to data_list_type (RJA)
!!
Module CFML_EoS
   !---- Use Modules ----!
   Use CFML_GlobalDeps
   Use CFML_Maths,     only: Debye,Second_Derivative,spline_interpol
   Use CFML_Metrics,   only: Cell_Type, Get_Cryst_Family, Set_Crystal_Cell
   Use CFML_Strings,   only: u_case, string_real, string_numstd

   !---- Definitions ----!
   implicit none

   private

   !---- Public procedures ----!
   public :: Alpha_Cal, Dkdt_Cal, Get_DebyeT, Get_GPT, Get_Grun_PT, Get_K, Get_Kp, Get_Pressure, Get_Pressure_Esd, &
             Get_Pressure_X, Get_Property_X, Get_Temperature, Get_Temperature_P0, Get_Transition_Pressure, &
             Get_Transition_Strain, Get_Transition_Temperature, Get_Volume, Get_Volume_S, K_Cal, Kp_Cal,   &
             Kpp_Cal, Pressure_F, Strain, Strain_EOS, Transition_phase, Get_volume_New

   public :: Allocate_EoS_Data_List, Allocate_EoS_List, Calc_Conlev, Check_scales, Deallocate_EoS_Data_List, Deallocate_EoS_List,    &
             Deriv_Partial_P, Deriv_Partial_P_Numeric, EoS_Cal, EoS_Cal_Esd, EosCal_text, EosParams_Check, FfCal_Dat, FfCal_Dat_Esd,&
             FfCal_EoS, Init_EoS_Cross, Init_EoS_Data_Type, Init_Eos_Shear, Init_Eos_Thermal, Init_EoS_Transition,     &
             Init_EoS_Type, Physical_check, Read_EoS_DataFile, Read_EoS_File, Read_Multiple_EoS_File,    &
             Set_Eos_Names, Set_Eos_Use, Set_Kp_Kpp_Cond, Write_Data_Conlev, Write_EoS_DataFile, Write_EoS_File,       &
             Write_Eoscal, Write_Eoscal_Header, Write_Info_Conlev, Write_Info_EoS, Define_Crystal_System


   !--------------------!
   !---- PARAMETERS ----!
   !--------------------!
   integer, public, parameter :: NCOL_DATA_MAX=22   ! Defines the maximum number of columns allowed in input data files

   integer, public, parameter :: N_EOSPAR=40        ! Specify the maximum number of Eos parameters allowed in Eos_type data type

   integer, public, parameter :: N_PRESS_MODELS=6   ! Number of possible pressure models
   integer, public, parameter :: N_THERM_MODELS=7   ! Number of possible Thermal models
   integer, public, parameter :: N_TRANS_MODELS=3   ! Number of possible Transition models
   integer, public, parameter :: N_SHEAR_MODELS=1   ! Number of possible Shear models
   integer, public, parameter :: N_CROSS_MODELS=2   ! Number of possible Cross-term models
   integer, public, parameter :: N_DATA_TYPES=2     ! Number of possible data types in addition to V (Kt,Ks etc)


   character(len=*), public, parameter, dimension(-1:N_PRESS_MODELS) :: PMODEL_NAMES=(/    &      ! Name of the Pressure Models
                                                                        'PTV Table      ', &
                                                                        'None           ', &
                                                                        'Murnaghan      ', &
                                                                        'Birch-Murnaghan', &
                                                                        'Vinet          ', &
                                                                        'Natural Strain ', &
                                                                        'Tait           ', &
                                                                        'APL2           '/)

   character(len=*), public, parameter, dimension(-1:N_THERM_MODELS) :: TMODEL_NAMES=(/        &  ! Name of the Thermal Models
                                                                        'PTV Table          ', &
                                                                        'None               ', &
                                                                        'Berman 1988        ', &
                                                                        'Fei 1995           ', &
                                                                        'Modified HP1998    ', &
                                                                        'Kroll              ', &
                                                                        'Salje low-T        ', &
                                                                        'HP Thermal Pressure', &
                                                                        'Mie-Gruneisen-Debye'/)

   character(len=*), public, parameter, dimension(-1:N_TRANS_MODELS) :: TRANMODEL_NAMES=(/ &      ! Name of Transition models
                                                                        'PTV Table      ', &
                                                                        'None           ', &
                                                                        'Landau P only  ', &
                                                                        'Landau T only  ', &
                                                                        'Landau PVT     '/)

   character(len=*), public, parameter, dimension(0:N_SHEAR_MODELS) :: SHEARMODEL_NAMES=(/ &      ! Name of Shear models
                                                                       'None           ',  &
                                                                       'Polynomial     '/)

   character(len=*), public, parameter, dimension(0:N_CROSS_MODELS) :: CROSSMODEL_NAMES=(/ &      ! Name of Cross-term models
                                                                       'None           ',  &
                                                                       'Linear dK/dT   ',  &
                                                                       'Hellfrich-Conn '/)

   character(len=*), public, parameter, dimension(0:N_DATA_TYPES) :: DATATYPE_NAMES=(/      &     ! Name of Data types
                                                                     'Cell parameters    ', &
                                                                     'Isothermal moduli  ', &
                                                                     'Adiabatic moduli   '/)

   real(kind=cp), public, parameter               :: AFERMIGAS    = 2337.0                                 ! Fermi Gas constant in GPa/A^5
   real(kind=cp), public, parameter, dimension(6) :: DELCHI       =(/ 2.30, 4.61, 6.17, 9.21,11.80,18.40/) ! Delta Chi2 values
   real(kind=cp), public, parameter, dimension(6) :: DELCHI_LEVELS=(/68.30,90.00,95.40,99.00,99.73,99.99/) ! Confidence Levels

   !---------------!
   !---- TYPES ----!
   !---------------!

   !!----
   !!----  TYPE :: PVT_Table
   !!--..
   !!---- Update: 23/09/2016
   !!
   Type, public :: PVT_Table
      integer                                      :: np       ! number of pressure lines
      integer                                      :: nt       ! number of temperature columns
      real(kind=cp)                                :: pmin     ! smallest pressure
      real(kind=cp)                                :: pmax     ! biggest pressure
      real(kind=cp)                                :: tmin     ! smallest temperature
      real(kind=cp)                                :: tmax     ! biggest temperature
      real(kind=cp), allocatable, dimension(:,:,:) :: ptv      ! The table, last index is 1=p, 2=t, 3=v
   End Type PVT_Table

   !!----
   !!----  TYPE :: EOS_TYPE
   !!--..
   !!---- Update: 23/09/2016
   !!
   Type, public :: EoS_Type
      character(len=80)                         :: Title=" "             ! Descriptive title of EoS, set by user
      character(len=15)                         :: Model=" "             ! Murnaghan, Birch-Murnaghan, Vinet, Natural-Strain
      character(len=15)                         :: TModel=" "            ! Name for thermal model
      character(len=15)                         :: TranModel=" "         ! Name for phase transition model
      character(len=15)                         :: SModel=" "            ! Name for shear model
      character(len=15)                         :: CModel=" "            ! Name for cross-terms model
      character(len=5), dimension(N_EOSPAR)     :: ParName=" "           ! Names of the Eos variables...init
      character(len=50), dimension(N_EOSPAR)    :: Comment=" "           ! Description of the Eos variables inclduing units...init
      character(len=15)                         :: Pscale_name=" "       ! Description of the Pressure scale. Only used for output
      character(len=15)                         :: Vscale_name=" "       ! Description of the Volume scale. Only used for output
      character(len=120), dimension(20)         :: doc=" "               ! Documentation for eos: set by user
      character(len=120)                        :: savedate=" "          ! Documentation on last save
      integer                                   :: IModel=0              ! Index for Model
      integer                                   :: IOrder=0              ! Order for the Model
      logical                                   :: Linear=.false.        ! Flag for Linear EoS not volume
      integer                                   :: ITherm=0              ! Index for thermal expansion model, =0 for none
      integer                                   :: ITran=0               ! Index for phase transition model, =0 for none
      integer                                   :: IShear=0              ! Index for shear model, =0 for none
      integer                                   :: ICross=0              ! Index for P-T cross-terms model, =0 for none or Pth
      integer, dimension(N_EOSPAR)              :: Iuse=0                ! Flags for parameters allowed for a given EoS =0 (not), =1 (refineable), =2 (implied non-zero)
      real(kind=cp)                             :: PRef=0.0              ! Pressure of Reference
      real(kind=cp)                             :: TRef=298.0            ! Temperature of Reference
      real(kind=cp)                             :: Stoich=0.0            ! Stocihometry factor for multiple-phase calculations
      real(kind=cp)                             :: Density0=0.0          ! Density at reference conditions
      logical                                   :: TRef_fixed=.false.    ! If true, then Tref cannot be changed by user
      logical                                   :: Pthermaleos=.false.   ! Indicates a Pthermal model if .true.
      real(kind=cp), dimension(N_EOSPAR)        :: Params=0.0            ! EoS Parameters
      real(kind=cp), dimension(N_EOSPAR)        :: Esd=0.0               ! Sigma EoS Parameters
      real(kind=cp)                             :: X=0.0                 ! Spare intensive variable, after P,T
      real(kind=cp)                             :: WChi2=0.0             ! weighted chi-squared for the refined parameters...set by ls
      real(kind=cp)                             :: DelPMax=0.0           ! Maximum misfit in P between obs and calc P...set by ls
      integer, dimension(4)                     :: IWt=0                 ! Choice for weighting in LS: order is P,T,V,X, 1 = use sigmas for weights
      integer, dimension(N_EOSPAR)              :: IRef=0                ! Refinement switches
      real(kind=cp),dimension(N_EOSPAR)         :: Factor=1.0            ! Scale factor to multiply variables for output (and divide on input)
      real(kind=cp)                             :: AlphaFactor=1.0E5_cp  ! Scale factor to multiply values of alpha (not parameters) for output
      real(kind=cp),dimension(N_EOSPAR)         :: Lastshift=0.0         ! Shift applied in last LS cycle to parameters
      real(kind=cp),dimension(N_EOSPAR,N_EOSPAR):: VCV=0.0               ! Var-Covar matrix from refinement
      Type(PVT_Table)                           :: Table                 ! A pvt table, used instead of eos parameters when imodel=-1
   End Type EoS_Type

   !!----
   !!---- TYPE :: EOS_LIST_TYPE
   !!--..
   !!---- Update: 17/07/2015
   !!
   Type, public :: EoS_List_Type
      integer                                   :: N=0    ! Number of EoS List
      character(len=30)                         :: system ! Crystal system name, including setting info (e.g. b-unique for mono)
      type(EoS_Type), allocatable, dimension(:) :: EoS    ! EoS Parameters
   End Type EoS_List_Type

   !!----
   !!----  TYPE :: EOS_DATA_TYPE
   !!--..
   !!---- Update: 17/07/2015
   !!
   Type, public :: EoS_Data_Type
      integer                     :: IUse=0    ! 0=No active, 1= active
      integer , dimension(5)      :: IGrp=0    ! Group
      integer                     :: xtype=0   ! Indicates type of data in V,cell, etc xtype=0 default, xtype=1 isothermal moduli etc
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
      character(len=80)                              :: Title=" "     ! Title of dataset (normally from input datafile)
      character(len=40)                              :: System=" "    ! Crystal System  (normally set by Def_Crystal_System)
      integer                                        :: N=0           ! Number of EoS Data List
      integer, dimension(NCOL_DATA_MAX)              :: IC_Dat=0      ! Which values are input
      character(len=15)                              :: Pscale_name=" "       ! Description of the Pressure scale of data (e.g. GPa)
      character(len=15)                              :: Vscale_name=" "       ! Description of the units of volume data (e.g. A3/cell)
      character(len=15)                              :: Lscale_name=" "       ! Description of the units of linear data  (e.g. A)
      type(EoS_Data_Type), allocatable, dimension(:) :: EoSD          ! EoS Data Parameters
   End Type EoS_Data_List_Type

   Interface
      Module Function Get_Pressure(V,T,EosPar) Result(P)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: V       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: P
      End Function Get_Pressure
      
      Module Function Get_Pressure_Esd(V,T,EosPar) Result(esd)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: V       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar 
         real(kind=cp)              :: EsD 
      End Function Get_Pressure_Esd
      
      Module Function Get_Pressure_K(K,T,Eos,Itype,Pest) Result(P)
         !---- Arguments ----!
         real(kind=cp),           intent(in) :: K       
         real(kind=cp),           intent(in) :: T       
         type(Eos_Type),          intent(in) :: EoS     
         integer,                 intent(in) :: itype
         real(kind=cp), optional, intent(in) :: Pest   
         real(kind=cp)                       :: P 
      End Function Get_Pressure_K
      
      Module Function Get_Pressure_X(X,T,Eos,Xtype,Pest) Result(P)
         !---- Arguments ----!
         real(kind=cp),           intent(in) :: x       
         real(kind=cp),           intent(in) :: T       
         type(Eos_Type),          intent(in) :: EoS     
         integer,                 intent(in) :: Xtype   
         real(kind=cp), optional, intent(in) :: Pest 
         real(kind=cp)                       :: P   
      End Function Get_Pressure_X
      
      Module Function Get_V0_T(T,EosPar) Result(V)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: T        
         type(Eos_Type), intent(in) :: EoSPar   
         real(kind=cp)              :: V   
      End Function Get_V0_T
      
      Module Function Get_Volume(P,T,EosPar) Result(v)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: P       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar 
         real(kind=cp)              :: V 
      End Function Get_Volume
      
      Module Function Get_Volume_New(P,T,EosPar) Result(v)
        !---- Arguments ----!
        real(kind=cp),  intent(in) :: P       
        real(kind=cp),  intent(in) :: T       
        type(Eos_Type), intent(in) :: EoSPar  
        real(kind=cp)              :: V
      End Function Get_Volume_New
      
      Module Function Get_Volume_K(K,T,E) Result(V)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: K      
         real(kind=cp),  intent(in) :: T      
         type(Eos_Type), intent(in) :: E    
         real(kind=cp)              :: V  
      End Function Get_Volume_K
      
      Module Function Get_Volume_K_old(K,T,E) Result(V)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: K      
         real(kind=cp),  intent(in) :: T      
         type(Eos_Type), intent(in) :: E 
         real(kind=cp)              :: V     
      End Function Get_Volume_K_old
      
      Module Function Get_Volume_S(S,T,Eospar) Result(V)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: S       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: V
      End Function Get_Volume_S
      
      Module Function Get_Transition_Pressure(T,EosPar) Result(Ptr)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: Ptr
      End Function Get_Transition_Pressure
      
      Module Function Get_Transition_Strain(P,T,EosPar) Result(vs)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: P       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar 
         real(kind=cp)              :: Vs 
      End Function Get_Transition_Strain
      
      Module Function Get_Transition_Temperature(P,EosPar) Result(Tr)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: P       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: Tr       
      End Function Get_Transition_Temperature
      
      Module Function Transition_Phase(P,T,Eospar) Result(Ip)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: P       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar 
         logical                    :: Ip 
      End Function Transition_Phase
      
      Module Function Get_Temperature(P,V,EosPar) Result(T)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: P       
         real(kind=cp),  intent(in) :: V       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: T
      End Function Get_Temperature
      
      Module Function Get_Temperature_P0(Vin,EosPar,Tmin,Tmax) Result(Tk)
         !---- Arguments ----!
         real(kind=cp),           intent(in) :: Vin       
         type(Eos_Type),          intent(in) :: EoSPar    
         real(kind=cp), optional, intent(in) :: Tmin      
         real(kind=cp), optional, intent(in) :: Tmax
         real(kind=cp)                       :: Tk
      End Function Get_Temperature_P0
      
      Module Function NormPressure_Eos(S,T,EosPar) Result(F)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: S       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: F
      End Function NormPressure_Eos
      
      Module Function NormPressure_P(S,P,imodel) Result(F)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: S      
         real(kind=cp),  intent(in) :: P      
         integer      ,  intent(in) :: imodel 
         real(kind=cp)              :: F
      End Function NormPressure_P
      
      Module Function Pressure_F(F,S,EosPar) Result(P)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: F       
         real(kind=cp),  intent(in) :: S       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: P
      End Function Pressure_F
      
      Module Function Strain(VV0,EosPar) Result(S)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: VV0     
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: S
      End Function Strain
      
      Module Function Strain_EOS(V,T,EosPar) Result(S)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: V       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar 
         real(kind=cp)              :: S 
      End Function Strain_EOS
      
      Module Function Get_K(P,T,EosPar) Result(K)
         !---- Arguments ----!
         real(kind=cp),intent(in)        :: p       
         real(kind=cp),intent(in)        :: t       
         type(Eos_Type),intent(in)       :: Eospar  
         real(kind=cp)                   :: K
      End Function Get_K
      
      Module Function Get_K0_T(T,Eospar) Result(k0)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: k0
      End Function Get_K0_T
      
      Module Function Get_Kp(P,T,EosPar) Result(Kp)
         !---- Arguments ----!
         real(kind=cp),intent(in)        :: p       
         real(kind=cp),intent(in)        :: t       
         type(Eos_Type),intent(in)       :: Eospar 
         real(kind=cp)                   :: kp 
      End  Function Get_Kp
      
      Module Function Get_Kp0_T(T,Eospar) Result(kp0)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar 
         real(kind=cp)              :: kp0  
      End Function Get_Kp0_T
      
      Module Function Get_Kpp0_T(T,Eospar) Result(kpp0)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar 
         real(kind=cp)              :: kpp0  
      End Function Get_Kpp0_T
      
      Module Function Get_GPT(P,T,EoSPar) Result(gpt)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: P      
         real(kind=cp),  intent(in) :: T      
         type(Eos_Type), intent(in) :: EosPar 
         real(kind=cp)              :: Gpt
      End Function Get_GPT
      
      Module Function Get_Grun_PT(P,T,Eospar) Result(G)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: P       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar 
         real(kind=cp)              :: G 
      End Function Get_Grun_PT
      
      Module Function Get_Grun_V(V,Eospar)  Result(Grun)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: V       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: Grun
      End Function Get_Grun_V
      
      Module Function K_Cal(V,T,Eospar,P) Result(kc)
         !---- Arguments ----!
         real(kind=cp),            intent(in) :: V       
         real(kind=cp),            intent(in) :: T       
         type(Eos_Type),           intent(in) :: EoSPar  
         real(kind=cp),  optional, intent(in) :: P  
         real(kind=cp)                        :: kc     
      End Function K_Cal
      
      Module Function Kp_Cal(V,T,EoSpar,P) Result (kpc)
         !---- Arguments ----!
         real(kind=cp),           intent(in) :: V       
         real(kind=cp),           intent(in) :: T       
         type(Eos_Type),          intent(in) :: EoSPar  
         real(kind=cp), optional, intent(in) :: P 
         real(kind=cp)                       :: kpc       
      End Function Kp_Cal
      
      Module Function Kpp_Cal(V,T,EoSpar) Result (kppc)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: V       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar
         real(kind=cp)              :: kppc    
      End Function Kpp_Cal
      
      Module Function Get_Props_PTVTable(P,T,V,EosPar,Res) Result(Val)
         !---- Arguments ----!
         real(kind=cp),    intent(in) :: P       
         real(kind=cp),    intent(in) :: T       
         real(kind=cp),    intent(in) :: V       
         type(Eos_Type),   intent(in) :: EoSPar  
         character(len=*), intent(in) :: Res 
         real(kind=cp)                :: Val    
      End Function Get_Props_PTVTable
      
      Module Function Murn_Interpolate_PTVTable(I,J,P,Eos) Result(V)
         !---- Arguments ----!
         integer,        intent(in) :: I      
         integer,        intent(in) :: J      
         real(kind=cp),  intent(in) :: P      
         type(Eos_Type), intent(in) :: EoS  
         real(kind=cp)              :: V
      End Function Murn_Interpolate_PTVTable
      
      Module Function EthDebye(T,Theta,Eospar) Result(Eth)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: T       
         real(kind=cp),  intent(in) :: Theta   
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: Eth
      End Function EthDebye
      
      Module Function Get_DebyeT(V, Eos) result(DebyeT)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: V       
         type(Eos_Type), intent(in) :: EoS 
         real(kind=cp)              :: DebyeT    
      End Function Get_DebyeT
      
      Module Function Pthermal(V,T,EosPar) Result(Pth)
         !---- Arguments ----!
         real(kind=cp),  intent(in) :: V       
         real(kind=cp),  intent(in) :: T       
         type(Eos_Type), intent(in) :: EoSPar  
         real(kind=cp)              :: PTh
      End Function Pthermal
      
      Module Function VscaleMGD(E) Result(MGD)
         !---- Arguments ----!
         type(Eos_Type),intent(in)  :: E          
         logical                    :: MGD        
      End Function VscaleMGD
      
      Module Function Get_Property_X(P,T,Eos,Xtype) Result(val)
         !---- Arguments ----!
         real(kind=cp),     intent(in) :: P       
         real(kind=cp),     intent(in) :: T       
         type(Eos_Type),    intent(in) :: EoS     
         integer, optional, intent(in) :: xtype   
         real(kind=cp)                 :: Val
      End Function Get_Property_X
      
      Module Function Alpha_Cal(P,T,Eospar, DeltaT) Result(Alpha)
         !---- Arguments ----!
         real(kind=cp),            intent(in) :: P       
         real(kind=cp),            intent(in) :: T       
         type(Eos_Type),           intent(in) :: EoSPar  
         real(kind=cp),  optional, intent(in) :: DeltaT  
      End Function Alpha_Cal
      
      Module Function dKdT_Cal(P,T,Eospar,DeltaT) Result(dKdT)
         !---- Arguments ----!
         real(kind=cp),            intent(in) :: P       
         real(kind=cp),            intent(in) :: T       
         type(Eos_Type),           intent(in) :: EoSPar  
         real(kind=cp),  optional, intent(in) :: DeltaT  
         real(kind=cp)                        :: dKdT
      End Function dKdT_Cal
      
      Module Subroutine Allocate_EoS_Data_List(N, E)
         !---- Arguments ----!
         integer,                    intent(in)       :: N  
         type (eos_data_list_type),  intent(in out)   :: E  
      End Subroutine Allocate_EoS_Data_List
      
      Module Subroutine Allocate_EoS_List(N, E)
         !---- Arguments ----!
         integer,               intent(in)       :: N  
         type (eos_list_type),  intent(in out)   :: E  
      End Subroutine Allocate_EoS_List
      
      Module Subroutine Deallocate_EoS_Data_List(E)
         !---- Arguments ----!
         type (eos_data_list_type), intent(in out)   :: E  
      End Subroutine Deallocate_EoS_Data_List
      
      Module Subroutine Deallocate_EoS_List(E)
         !---- Arguments ----!
         type (eos_list_type), intent(in out)   :: E  
      End Subroutine Deallocate_EoS_List
      
      Module Subroutine Init_EoS_Cross(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar
      End Subroutine Init_EoS_Cross
      
      Module Subroutine Init_EoS_Data_Type(E)
         !---- Arguments ----!
         type (EoS_Data_Type), intent(in out)   :: E
      End Subroutine Init_EoS_Data_Type
      
      Module Subroutine Init_EoS_Shear(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar
      End Subroutine Init_EoS_Shear
      
      Module Subroutine Init_EoS_Thermal(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar
      End Subroutine Init_EoS_Thermal
      
      Module Subroutine Init_EoS_Transition(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar
      End Subroutine Init_EoS_Transition
      
      Module Subroutine Init_EoS_Type(Eospar,CLin,IThermal,ITransition,Ishear,Icross)
         !---- Arguments ----!
         type (EoS_Type),            intent(out)    :: Eospar       
         character(len=*), optional, intent(in)     :: CLin         
         integer,          optional, intent(in)     :: IThermal     
         integer,          optional, intent(in)     :: ITransition  
         integer,          optional, intent(in)     :: IShear       
         integer,          optional, intent(in)     :: ICross       
      End Subroutine Init_EoS_Type
      
      Module Subroutine Set_Cross_Names(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar  
      End Subroutine Set_Cross_Names
      
      Module Subroutine Set_EoS_Factors(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar
      End Subroutine Set_EoS_Factors
      
      Module Subroutine Set_Eos_Names(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar   
      End Subroutine Set_Eos_Names
      
      Module Subroutine Set_EoS_Use(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  
      End Subroutine Set_EoS_Use
      
      Module Subroutine Set_Kp_Kpp_Cond(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar  
      End Subroutine Set_Kp_Kpp_Cond
      
      Module Subroutine Set_Shear_Names(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar  
      End Subroutine Set_Shear_Names
      
      Module Subroutine Set_Thermal_Names(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar  
      End Subroutine Set_Thermal_Names
      
      Module Subroutine Set_Transition_Names(Eospar)
         !---- Arguments ----!
         type (EoS_Type), intent(in out) :: Eospar  
      End Subroutine Set_Transition_Names

      Module Subroutine Read_EoS_DataFile(fname,dat)
         !---- Arguments ----!
         character(len=*),          intent(in)  :: fname   
         type (eos_data_list_type), intent(out) :: dat     
      End Subroutine Read_EoS_DataFile
      
      Module Subroutine Read_Eos_File(Fname,Eos)
         !---- Arguments ----!
         character(len=*),intent(in)  :: fname  
         type (EoS_Type), intent(out) :: Eos    
      End Subroutine Read_Eos_File
      
      Module Subroutine Read_Eos_In(Flines,Eos)
         !---- Arguments ----!
         character(len=*),dimension(:),intent(in)   :: flines
         type(Eos_type),               intent(out)  :: eos
      End Subroutine Read_Eos_In
      
      Module Subroutine Read_Multiple_Eos_File(Fname,Eoslist)
         !---- Arguments ----!
         character(len=*),     intent(in)  :: fname     
         type (EoS_List_Type), intent(out) :: Eoslist   
      End Subroutine Read_Multiple_Eos_File
      
      Module Subroutine Write_EoS_DataFile(dat,lun)
         !---- Arguments ----!
         type (eos_data_list_type), intent(in) :: dat  
         integer,                   intent(in) :: lun  
      End Subroutine Write_EoS_DataFile
      
      Module Subroutine Write_Eos_File(Eos,Lun)
         !---- Arguments ----!
         type (EoS_Type),intent(in)   :: Eos 
         integer,intent(in)           :: lun 
      End Subroutine Write_Eos_File
      
      Module Subroutine Write_Eoscal(Pmin,Pmax,Pstep,Tmin,Tmax,Tstep,Tscale_In,Eos,Lun,Nprint,eoscal_err)
         !---- Arguments ----!
         real(kind=cp),    intent(in)  ::  pmin, pmax, pstep   
         real(kind=cp),    intent(in)  ::  tmin,tmax,tstep     
         character(len=*), intent(in)  ::  tscale_in           
         type(EoS_Type),   intent(in)  ::  eos                 
         integer,          intent(in)  ::  lun                  
         integer,          intent(out) ::  nprint               
         logical,          intent(out) ::  eoscal_err           
      End Subroutine Write_Eoscal
      
      Module Subroutine Write_Eoscal_Header(Eos,Lun,Tscale_In)
         !---- Arguments ----!
         type(EoS_Type),intent(in)   :: eos         
         integer,       intent(in)   :: lun         
         character(len=*),intent(in) :: tscale_in   
      End Subroutine Write_Eoscal_Header
      
      Module Subroutine Write_Info_Eos(Eospar,iout)
         !---- Arguments ----!
         type(Eos_Type),    intent(in) :: Eospar  
         integer, optional, intent(in) :: iout    
      End Subroutine Write_Info_Eos
      
      Module Subroutine Write_Info_Eos_Cross(Eos,Iout)
         !---- Arguments ----!
         type(Eos_Type),    intent(in) :: Eos
         integer, optional, intent(in) :: iout
      End Subroutine Write_Info_Eos_Cross
      
      Module Subroutine Write_Info_Eos_Shear(Eos,Iout)
         !---- Arguments ----!
         type(Eos_Type),    intent(in) :: Eos
         integer, optional, intent(in) :: iout
      End Subroutine Write_Info_Eos_Shear
      
      Module Subroutine Write_Info_Eos_Thermal(Eospar,iout)
         !---- Arguments ----!
         type(Eos_Type),    intent(in) :: Eospar  
         integer, optional, intent(in) :: iout    
      End Subroutine Write_Info_Eos_Thermal
      
      Module Subroutine Write_Info_Eos_Transition(Eospar,iout)
         !---- Arguments ----!
         type(Eos_Type),    intent(in) :: Eospar  
         integer, optional, intent(in) :: iout    
      End Subroutine Write_Info_Eos_Transition
      
      Module Subroutine Calc_Conlev(Eos,ix,iy,isig,xyy,n)
         !---- Arguments ----!
         type(Eos_Type),              intent(in)  :: Eos         
         integer,                     intent(in)  :: ix       
         integer,                     intent(in)  :: iy       
         integer,                     intent(in)  :: isig        
         real(kind=cp),dimension(:,:),intent(out) :: xyy         
         integer,                     intent(out) :: n          
      End Subroutine Calc_Conlev
      
      Module Subroutine Write_Data_Conlev(xyy,n,iout)
         !---- Arguments ----!
         real(kind=cp),dimension(:,:),intent(in)   :: xyy  
         integer,                     intent(in)   :: n    
         integer, optional,           intent(in)   :: iout 
      End Subroutine Write_Data_Conlev
      
      Module Subroutine Write_Info_Conlev(Eos,ix,iy,isig,iout)
         !---- Arguments ----!
         type(Eos_Type),    intent(in) :: Eos
         integer,           intent(in) :: ix  
         integer,           intent(in) :: iy  
         integer,           intent(in) :: isig  
         integer, optional, intent(in) :: iout
      End Subroutine Write_Info_Conlev
      
      Module Function Get_Tait(Eospar,T) Result(Vec)
         !---- Arguments ----!
         type(Eos_Type),            intent(in)  :: Eospar  
         real(kind=cp),             intent(in)  :: T       
         real(kind=cp),dimension(3)             :: Vec     
      End Function Get_Tait
      
      Module Subroutine Deriv_Partial_P(V,T,Eospar,td,xtype,calc)
         !---- Arguments ----!
         real(kind=cp),                      intent(in)  :: V       
         real(kind=cp),                      intent(in)  :: T       
         type(Eos_Type),                     intent(in)  :: Eospar  
         real(kind=cp), dimension(n_eospar), intent(out) :: td      
         integer,optional,                   intent(in)  :: xtype   
         character(len=*),optional,          intent(in)  :: calc    
      End Subroutine Deriv_Partial_P
      
      Module Subroutine Deriv_Partial_P_Analytic(V,T,Eospar,td)
         !---- Arguments ----!
         real(kind=cp),                      intent(in)  :: V       
         real(kind=cp),                      intent(in)  :: T       
         type(Eos_Type),                     intent(in)  :: Eospar  
         real(kind=cp), dimension(n_eospar), intent(out) :: td      
      End Subroutine Deriv_Partial_P_Analytic
      
      Module Subroutine Deriv_Partial_P_Numeric(X1,X2,Eospar,td,xtype,calc)
         !---- Arguments ----!
         real(kind=cp),                      intent(in) :: X1,X2   
         type(Eos_Type),                     intent(in) :: Eospar  
         real(kind=cp), dimension(n_eospar), intent(out):: td      
         integer,         optional,          intent(in) :: xtype   
         character(len=*),optional,          intent(in) :: calc    
      End Subroutine Deriv_Partial_P_Numeric 
      
      Module Subroutine FfCal_Dat(V,V0,P,Eospar,F,S)
         !---- Arguments ----!
         real(kind=cp),           intent(in)  :: V       
         real(kind=cp),           intent(in)  :: V0      
         real(kind=cp),           intent(in)  :: P       
         type(Eos_Type),          intent(in)  :: EoSPar  
         real(kind=cp), optional, intent(out) :: F       
         real(kind=cp), optional, intent(out) :: S       
      End Subroutine FfCal_Dat
      
      Module Subroutine FfCal_Dat_Esd(V,SigV,V0,SigV0,P,SigP,EosPar,F,SigF,S,SigS)
         !---- Arguments ----!
         real(kind=cp),  intent(in)  :: V       
         real(kind=cp),  intent(in)  :: SigV    
         real(kind=cp),  intent(in)  :: V0      
         real(kind=cp),  intent(in)  :: SigV0   
         real(kind=cp),  intent(in)  :: P       
         real(kind=cp),  intent(in)  :: SigP    
         type(Eos_Type), intent(in)  :: EoSPar  
         real(kind=cp),  intent(out) :: F       
         real(kind=cp),  intent(out) :: SigF    
         real(kind=cp),  intent(out) :: S       
         real(kind=cp),  intent(out) :: SigS    
      End Subroutine FfCal_Dat_Esd
      
      Module Subroutine FfCal_EoS(P,T,Eospar,F,S)
         !---- Arguments ----!
         real(kind=cp),           intent(in)     :: P       
         real(kind=cp),           intent(in)     :: T       
         type(Eos_Type),          intent(in)     :: Eospar  
         real(kind=cp),           intent(out)    :: F       
         real(kind=cp),           intent(out)    :: S       
      End Subroutine FfCal_EoS
      
      Module Subroutine Murn_PTVTable(i,j,Eos,Eosm)
         !---- Arguments ----!
         integer,        intent(in)  :: i      
         integer,        intent(in)  :: j      
         type(Eos_Type), intent(in)  :: EoS    
         type(Eos_Type), intent(out) :: EoSm   
      End Subroutine Murn_PTVTable
      
      Module Subroutine Check_Scales(E,dat)
         !---- Arguments ----!
         type(Eos_Type),                       intent(in)  :: E     
         type (eos_data_list_type), optional,  intent(in)  :: dat   
      End  Subroutine Check_Scales
      
      Module Subroutine EoSParams_Check(EoSPar)
         !---- Argument ----!
         type (EoS_Type), intent(in out) :: EoSPar
      End Subroutine EoSParams_Check
      
      Module Subroutine Physical_Check(Ein,Pin,Tin,Vin)
         !---- Arguments ----!
         type(Eos_Type),        intent(in) :: Ein 
         real(kind=cp),optional,intent(in) :: pin 
         real(kind=cp),optional,intent(in) :: tin 
         real(kind=cp),optional,intent(in) :: vin 
      End Subroutine Physical_Check
      
      Module Subroutine Physical_Check_old(P,T,E)
         !---- Arguments ----!
         real(kind=cp),    intent(in) :: p  
         real(kind=cp),    intent(in) :: t  
         type(Eos_Type),   intent(in) :: E  
      End Subroutine Physical_Check_old
      
      Module Subroutine PVEos_Check(Pin,Vin,Tin,ein,vpresent)
         !---- Arguments ----!
         real(kind=cp),optional,intent(in) :: pin  
         real(kind=cp),optional,intent(in) :: vin  
         real(kind=cp),optional,intent(in) :: tin  
         type(Eos_Type),        intent(in) :: Ein  
         logical,               intent(in) :: vpresent 
      End Subroutine PVEos_Check
      
      Module Function Eoscal_Text(P,T,Tscale_In,Eos) Result(text)
         !---- Arguments ----!
         real(kind=cp),    intent(in)  :: p             
         real(kind=cp),    intent(in)  :: t             
         character(len=*), intent(in)  :: tscale_in     
         type(EoS_Type),   intent(in)  :: eos           
         character(len=:), allocatable :: text             
      End Function Eoscal_Text
      
      Module Subroutine Transform_Esd(P,T,Eospar,Esd)
         !---- Arguments ----!
         real(kind=cp),   intent(in)              :: p        
         real(kind=cp),   intent(in)              :: t        
         type (EoS_Type), intent(in)              :: eospar   
         real(kind=cp),dimension(:), intent(out)  :: esd      
      End Subroutine Transform_Esd
      
      Module Subroutine EoS_Cal(P,T,EoSpar,Parvals)
         !---- Arguments ----!
         real(kind=cp),                intent(in)  :: P       
         real(kind=cp),                intent(in)  :: T       
         type(Eos_Type),               intent(in)  :: EoSPar  
         real(kind=cp),  dimension(:), intent(out) :: Parvals 
      End Subroutine EoS_Cal
      
      Module Subroutine EoS_Cal_Esd(P,T,EoSpar,Esd)
         !---- Arguments ----!
         real(kind=cp),                intent(in)  :: P       
         real(kind=cp),                intent(in)  :: T       
         type(Eos_Type),               intent(in)  :: EoSPar  
         real(kind=cp),  dimension(:), intent(out) :: Esd     
      End Subroutine EoS_Cal_Esd
      
   End Interface
   
Contains
   !!--++
   !!--++ EOS_TO_VEC
   !!--++
   !!--++ PRIVATE
   !!--++ Copy parameters from EoS type to a vector
   !!--++
   !!--++ 28/02/2013
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
         select case(eospar%imodel)
            case (6)
               vec(2:3)=eospar%params(2:3)/3.0_cp

            case default
               vec(2:4)=eospar%params(2:4)/3.0_cp
         end select

         select case(eospar%icross)
            case (1)
               vec(5)=eospar%params(5)/3.0_cp

            case (2)
               vec(5:6)=eospar%params(5:6)
         end select

         select case(eospar%itherm)                 ! thermal expansion terms
            case(1:3)
               vec(10:12)=eospar%params(10:12)*3.0_cp

            case(4,5,6)
               vec(10)=eospar%params(10)*3.0_cp
               vec(11)=eospar%params(11)
         end select

         !> phase transition: no changes required for linear
      end if
   End Subroutine EoS_to_Vec
   
   !!--++
   !!--++ VEC_TO_EOS
   !!--++
   !!--++ PRIVATE
   !!--++ Pass values fron a vector to respective EoS parameter
   !!--++
   !!--++ 28/02/2013
   !!
   Subroutine Vec_to_EoS(Vec,EosPar)
      !---- Arguments ----!
      real(kind=cp), dimension(:), intent(in)     :: Vec
      type(EoS_Type),              intent(in out) :: Eospar

      !> Copy vec to eosparams as 1 to 1 as default
      eospar%params=vec

      if (eospar%linear) then
         eospar%params(1)=vec(1)**(1.0_cp/3.0_cp)
         select case(eospar%imodel)
            case (6)
               eospar%params(2:3)=vec(2:3)*3.0_cp

            case default
               eospar%params(2:5)=vec(2:5)*3.0_cp          ! recheck cross terms******
         end select

         select case(eospar%itherm)                 ! thermal expansion terms
            case (1,2,3)
               eospar%params(10:12)=vec(10:12)/3.0_cp

            case (4,5,6)
               eospar%params(10)=vec(10)/3.0_cp
               eospar%params(11)=vec(11)
         end select
      end if
   End Subroutine Vec_to_EoS

   !!--++
   !!--++ DEFINE_CRYSTAL_SYSTEM
   !!--++
   !!--++ PRIVATE
   !!--++ Either sets cell parameters to conform to specifed system Or,
   !!--++ if no system specified, tries to determine system if all cell parameters
   !!--++ provided
   !!--++
   !!--++ 17/07/2015
   !!
   Subroutine Define_Crystal_System(dat)
      !---- Arguments ----!
      type (eos_data_list_type),  intent(in out)   :: dat  ! data structure

      !---- Local Variables ----!
      character(len=40)       :: Family, SystemC, system
      character(len=1)        :: Symbol
      integer                 :: i,ndat
      type(Cell_Type)         :: ncell

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
               end if

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

               !> b=a
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
               !> Angles =90º
               dat%eosd(1:ndat)%ang(1)=90.0
               dat%eosd(1:ndat)%ang(2)=90.0
               dat%eosd(1:ndat)%ang(3)=90.0
               dat%eosd(1:ndat)%siga(1)=0.0
               dat%eosd(1:ndat)%siga(2)=0.0
               dat%eosd(1:ndat)%siga(3)=0.0

               !> a=b=c: Modified RJA 14 Jan to handle case if V is supplied, but 'a' is not
               do i=1,ndat
                  if (dat%eosd(i)%cell(1) < tiny(0.0) ) then
                     dat%eosd(i)%cell(1)=dat%eosd(i)%v**(1.0_cp/3.0_cp)
                     dat%eosd(i)%sigc(1)=dat%eosd(i)%sigv/3.0_cp/(dat%eosd(i)%cell(1)**2.0_cp)
                  end if
                  dat%eosd(i)%cell(2)=dat%eosd(i)%cell(1)
                  dat%eosd(i)%cell(3)=dat%eosd(i)%cell(1)
                  dat%eosd(i)%sigc(2)=dat%eosd(i)%sigc(1)
                  dat%eosd(i)%sigc(3)=dat%eosd(i)%sigc(1)
               end do
         end select

      else
         !> checking cell parameters
         system=' '
         if (any(dat%eosd(1)%cell <= 0.1)) return
         if (any(dat%eosd(1)%ang <= 0.5)) return

         !> Obtaining system from cell parameters
         ncell%cell=dat%eosd(1)%cell
         ncell%ang =dat%eosd(1)%ang
         call Get_Cryst_Family(ncell, Family, Symbol, System)
         system=adjustl(system)
         if (len_trim(system) <= 0) then
            err_CFML%IErr=1
            Err_CFML%Msg="Cell values for first data in the data file are inconsistent"
            return
         end if
         
         !> Now check if all data the same. 
         !> If they are not, set system=' ' because there may be a transition
         do i=2,ndat
            ncell%cell=dat%eosd(i)%cell
            ncell%ang =dat%eosd(i)%ang
            call Get_Cryst_Family(ncell, Family, Symbol, SystemC)
            systemC=adjustl(systemc)

            if (systemC /= System) then
               dat%system=' '
               return
            end if
         end do
      end if
   End Subroutine Define_Crystal_System

   !!--++
   !!--++ SET_VOLUME_FROM_CELL
   !!--++
   !!--++ PRIVATE
   !!--++ Sets V and esd(V) from cell parameter data for all data items in dat
   !!--++ If V is present in first data item, no esd is calculated
   !!--++
   !!--++ 17/07/2015
   !!
   Subroutine Set_Volume_from_Cell(dat)
      !---- Arguments ----!
      type (eos_data_list_type),  intent(in out)  :: dat  ! data structure

      !---- Local Variables ----!
      integer         :: i
      type(Cell_Type) :: cell

      !> Check
      if (dat%eosd(1)%v > 0.0) return
      if (any(dat%eosd(1)%cell < 0.5)) return

      !> Calculations
      do i=1,dat%n
         call Set_Crystal_Cell(dat%eosd(i)%cell,dat%eosd(i)%ang,cell, &
                               Vscell=dat%eosd(i)%sigc,Vsang=dat%eosd(i)%siga)
         dat%eosd(i)%v=cell%vol
         dat%eosd(i)%sigv=cell%svol
      end do
   End Subroutine Set_Volume_from_Cell

End Module CFML_EoS
