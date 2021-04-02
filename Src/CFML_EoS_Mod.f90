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
!!----                          Universita di Pavia, Pavia, ITALY
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL) (CrysFML)
!!----          Javier Gonzalez-Platas  (ULL) (CrysFML)
!!----          Ross John Angel         (Padova)  (EoS)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----               Ross John Angel    (IGG-CNR, Italy)
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
!!----    ...   2019: Bug fixes, addition of new thermal-pressure EoS types (RJA)
!!----    ...   2020: Jan: E%LinearDir to label linear EoS (RJA)
!!----    ...   2020: March: Addition of extra oscillators for thermal pressure, addition of scales block to eos%params (RJA)
!!----    ...   2020: Sept: Moved enquiry routines for groups into this module (RJA)
!!----    ...   2020: Sept: Moved calculation and management routines for eos of cells into this module (RJA)
Module CFML_EoS
   !---- Use Modules ----!
   Use CFML_GlobalDeps,       only: CP, PI, TO_RAD
   Use CFML_Math_General,     only: Debye,ERR_MathGen,ERR_MathGen_Mess,First_Derivative,Second_Derivative, &
                                    splint,Diagonalize_SH,tand,cosd,sind,asind
   Use CFML_Crystal_Metrics,  only: Crystal_Cell_Type,Get_Cryst_Family,Volume_Sigma_from_Cell,Strain_Tensor_type,&
                                    Fix_tensor,Set_Crystal_Cell,Orient_Eigenvectors,Calc_Paxes_Angles,Init_Strain_Tensor
   Use CFML_String_Utilities

   !---- Definitions ----!
   implicit none

   private

   !---- Public procedures ----!
   public :: Alpha_Cal, &
             Deriv_Partial_P, Deriv_Partial_P_Numeric, Deriv_Partial_P_Scales, dKdt_Cal, &
             EoS_Cal, EoS_Cal_Esd, &
             Get_Alpha_Cell, Get_Angle_Deriv, Get_Cp, Get_Cv, Get_DebyeT, Get_GPT, Get_Grun_PT, &
             Get_Grun_Th,Get_Grun_V, Get_K, Get_Kp, Get_Mod_Axis, Get_Mod_Cell,                 &
             Get_Modp_Axis, Get_Modp_Cell, Get_Press_Axis, Get_Press_Cell,                      &
             Get_Pressure, Get_Pressure_Esd, Get_Pressure_X, Get_Property_X, Get_Temperature,   &
             Get_Temperature_P0,                                                                &
             Get_Transition_Pressure, Get_Transition_Strain, Get_Transition_Temperature,        &
             Get_Volume, Get_Volume_Axis, Get_Volume_Cell, Get_Volume_S, Get_Params_Cell,       &
             Get_Props_General, Get_Props_Third, Isotropic_Cell,  &
             K_Cal, Kp_Cal, Kpp_Cal, &
             Linear_EoS_Allowed, &
             Pressure_F, Principal_Eos, Pthermal, &
             Set_XdataTypes, Strain, Strain_EOS, &
             Thermal_Pressure_Eos, Transition_Phase, &
             VscaleMGD


   public :: Allocate_EoS_Data_List, Allocate_EoS_List, &
             Calc_Conlev, Check_scales, Copy_Eos_Data_List, &
             Deallocate_EoS_Data_List, Deallocate_EoS_List, Def_Crystal_System, &
             EosCal_text, Eos_Cell_Loaded_Check, EosParams_Check,  &
             FfCal_Dat, FfCal_Dat_Esd, FfCal_EoS, &
             Get_Tensor_Eos, &
             Init_EoS_Angles, Init_Eos_Cell_Type, Init_EoS_Cross, Init_EoS_Data_Type, &
             Init_EoS_GroupScales, Init_EoS_Osc, Init_Eos_Shear, Init_Eos_Thermal,    &
             Init_EoS_Transition, Init_EoS_Type, Init_Err_EoS,                        &
             Physical_check, &
             Read_EoS_DataFile, Read_EoS_File, Read_Multiple_EoS_File,                &
             Set_Cell_Types, Set_Eos_Implied_Values, Set_Eos_Names, Set_Eos_Use,      &
             Write_Data_Conlev, Write_EoS_DataFile, Write_EoS_File, Write_Eoscal,     &
             Write_Eoscal_Header, Write_Info_Conlev, Write_Info_EoS, Write_Info_Eos_Cell_Type


   !--------------------!
   !---- PARAMETERS ----!
   !--------------------!
   integer, public, parameter :: NCOL_DATA_MAX=22   ! Defines the maximum number of columns allowed in input data files

   integer, public, parameter :: N_EOSPAR=59        ! Specify the maximum number of Eos parameters allowed in Eos_type data type
   integer, public, parameter :: N_ANGPOLY=3        ! Dimension of polynomial for angles

   integer, public, parameter :: N_PRESS_MODELS=7   ! Number of possible pressure models
   integer, public, parameter :: N_THERM_MODELS=8  ! Number of possible Thermal models
   integer, public, parameter :: N_TRANS_MODELS=3   ! Number of possible Transition models
   integer, public, parameter :: N_SHEAR_MODELS=1   ! Number of possible Shear models
   integer, public, parameter :: N_CROSS_MODELS=2   ! Number of possible Cross-term models
   integer, public, parameter :: N_OSC_MODELS=2     ! Number of possible extra oscillator model types.
   integer, public, parameter :: N_ANGLE_MODELS=1   ! Number of possible angle polynomial models
   integer, public, parameter :: N_DATA_TYPES=2     ! Number of possible data types in addition to V (Kt,Ks etc)


   character(len=*), public, parameter, dimension(-1:N_PRESS_MODELS) :: PMODEL_NAMES=(/    &      ! Name of the Pressure Models
                                                                        'PTV Table      ', &
                                                                        'None           ', &
                                                                        'Murnaghan      ', &
                                                                        'Birch-Murnaghan', &
                                                                        'Vinet          ', &
                                                                        'Natural Strain ', &
                                                                        'Tait           ', &
                                                                        'APL            ', &
                                                                        'Kumar          '/)

   character(len=*), public, parameter, dimension(-1:N_THERM_MODELS) :: TMODEL_NAMES=(/        &  ! Name of the Thermal Models
                                                                        'PTV Table          ', &
                                                                        'None               ', &
                                                                        'Berman 1988        ', &
                                                                        'Fei 1995           ', &
                                                                        'Modified HP1998    ', &
                                                                        'Kroll              ', &
                                                                        'Salje low-T        ', &
                                                                        'HP Thermal Pressure', &
                                                                        'Mie-Gruneisen-Debye', &
                                                                        'Einstein Oscillator'/)

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

   character(len=*), public, parameter, dimension(0:N_OSC_MODELS)   :: OSCMODEL_NAMES=(/ &      ! Name of extra oscillator types
                                                                       'None           ',  &
                                                                       'Debye          ',  &
                                                                       'Einstein       '/)

   character(len=*), public, parameter, dimension(0:N_DATA_TYPES) :: DATATYPE_NAMES=(/      &     ! Name of Data types
                                                                     'Cell parameters    ', &
                                                                     'Isothermal moduli  ', &
                                                                     'Adiabatic moduli   '/)

   real(kind=cp), public, parameter               :: AFERMIGAS    = 2337.0                                 ! Fermi Gas constant in GPa/A^5
   real(kind=cp), public, parameter, dimension(6) :: DELCHI       =(/ 2.30, 4.61, 6.17, 9.21,11.80,18.40/) ! Delta Chi2 values
   real(kind=cp), public, parameter, dimension(6) :: DELCHI_LEVELS=(/68.30,90.00,95.40,99.00,99.73,99.99/) ! Confidence Levels


   character(len=*), public, parameter, dimension(4:21) :: DATA_NAMES=(/     &    !Names of data variables in ic_dat in EoS_Data_List_Type
                                                                       'T    ', &
                                                                       'SIGT ', &
                                                                       'P    ', &
                                                                       'SIGP ', &
                                                                       'V    ', &
                                                                       'SIGV ', &
                                                                       'A    ', &
                                                                       'SIGA ', &
                                                                       'B    ', &
                                                                       'SIGB ', &
                                                                       'C    ', &
                                                                       'SIGC ', &
                                                                       'ALPHA', &
                                                                       'SIGAL', &
                                                                       'BETA ', &
                                                                       'SIGBE', &
                                                                       'GAMMA', &
                                                                       'SIGGA'/)

   character(len=*), public, parameter, dimension(1:7) :: CELLLABEL =  (/'a    ', &
                                                                         'b    ', &
                                                                         'c    ', &
                                                                         'alpha', &
                                                                         'beta ', &
                                                                         'gamma', &
                                                                         'Vol  '/)   !labels for unit-cell parameters

   character(len=*), public, parameter, dimension(0:6) :: AXISLABEL =  (/'V   ', &
                                                                         'a   ', &
                                                                         'b   ', &
                                                                         'c   ', &
                                                                         'd100', &
                                                                         'd010', &
                                                                         'd001'/)  !labels for axis eos
   !---------------!
   !---- TYPES ----!
   !---------------!

   !!----
   !!----  TYPE :: PVT_Table
   !!--..
   !!---- Update: 23/09/2016
   !!
   Type, public :: PVT_Table
      integer                                      :: np  =0        ! number of pressure lines
      integer                                      :: nt  =0        ! number of temperature columns
      real(kind=cp)                                :: pmin=0.0      ! smallest pressure
      real(kind=cp)                                :: pmax=1.0E6    ! biggest pressure
      real(kind=cp)                                :: tmin=0.0      ! smallest temperature
      real(kind=cp)                                :: tmax=1.0E6    ! biggest temperature
      real(kind=cp), allocatable, dimension(:,:,:) :: ptv           ! The table, last index is 1=p, 2=t, 3=v
   End Type PVT_Table

   !!----
   !!----  TYPE :: EOS_TYPE
   !!--..
   !!---- Update: 23/09/2016
   !!
   Type, public :: EoS_Type
      character(len=80)                         :: Title=" "             ! Descriptive title of EoS, set by user
      character(len=20)                         :: System=" "            ! Crystal system name
      character(len=15)                         :: Model=" "             ! Murnaghan, Birch-Murnaghan, Vinet, Natural-Strain
      character(len=20)                         :: TModel=" "            ! Name for thermal model
      character(len=15)                         :: TranModel=" "         ! Name for phase transition model
      character(len=15)                         :: SModel=" "            ! Name for shear model
      character(len=15)                         :: CModel=" "            ! Name for cross-terms model
      character(len=25),dimension(2)            :: OscModel=" "          ! Names for oscillator models
      character(len=5), dimension(N_EOSPAR)     :: ParName=" "           ! Names of the Eos variables...init
      character(len=50), dimension(N_EOSPAR)    :: Comment=" "           ! Description of the Eos variables inclduing units...init
      character(len=15)                         :: Pscale_name=" "       ! Description of the Pressure scale. Only used for output
      character(len=15)                         :: Vscale_name=" "       ! Description of the Volume scale. Only used for output
      character(len=120), dimension(20)         :: doc=" "               ! Documentation for eos: set by user
      character(len=120)                        :: savedate=" "          ! Documentation on last save
      character(len=32)                         :: LinearDir=" "         ! Label for linear direction; default blank
      integer                                   :: IModel=0              ! Index for Model
      integer                                   :: IOrder=0              ! Order for the Model
      logical                                   :: Linear=.false.        ! Flag for Linear EoS not volume
      integer                                   :: ITherm=0              ! Index for thermal expansion model, =0 for none
      integer                                   :: ITran=0               ! Index for phase transition model, =0 for none
      integer                                   :: IShear=0              ! Index for shear model, =0 for none
      integer                                   :: ICross=0              ! Index for P-T cross-terms model, =0 for none or Pth
      integer                                   :: IAngle=0              ! Index for angle polynomial and not eos
      integer,dimension(2)                      :: IOsc=0                ! Index for extra oscillator models
      integer, dimension(N_EOSPAR)              :: Iuse=0                ! Flags for parameters allowed for a given EoS =0 (not), =1 (refineable), =2 (implied non-zero)
      real(kind=cp)                             :: PRef=0.0              ! Pressure of Reference
      real(kind=cp)                             :: TRef=298.0            ! Temperature of Reference
      real(kind=cp)                             :: Stoich=0.0            ! Stocihometry factor for multiple-phase calculations
      real(kind=cp)                             :: Density0=0.0          ! Density at reference conditions
      logical                                   :: TRef_fixed=.false.    ! If true, then Tref cannot be changed by user
      logical                                   :: Pthermaleos=.false.   ! Indicates a Pthermal model if .true.
      logical                                   :: Osc_allowed=.false.   ! Indicates if extra oscillators can be used with current thermal model
      logical, dimension(2:4)                   :: allowed_orders        ! Indicates which orders (.true.) are allowed for the PV Eos
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
      real(kind=cp),dimension(3,0:3,N_ANGPOLY)  :: angpoly=0.0           ! Polynomial coefficients for unit-cell angles
      Type(PVT_Table)                           :: Table                 ! A pvt table, used instead of eos parameters when imodel=-1
   End Type EoS_Type

   !!----
   !!---- TYPE :: EOS_LIST_TYPE
   !!--..
   !!---- Update: 17/07/2015
   !!
   Type, public :: EoS_List_Type
      integer                                   :: N=0          ! Number of EoS List
      character(len=30)                         :: system=" "   ! Crystal system name, including setting info (e.g. b-unique for mono)
      type(EoS_Type), allocatable, dimension(:) :: EoS          ! EoS Parameters
   End Type EoS_List_Type

   !!----
   !!---- TYPE :: EOS_CELL_TYPE
   !!--..
   !!---- New 14/02/2020. Specific list of eos and pointers for full description of a unit cell. RJA
   !!
   Type, public :: EoS_Cell_Type
      integer                                     :: N=0              ! Max index of used EoS  in List - depends on crystal system:
      character(len=30)                           :: system=" "       ! Crystal system name, including setting info (e.g. b-unique for mono)
      type(EoS_Type),dimension(0:6)               :: EoS              ! EoS Parameters for V,a,b,c,d100,d010,d001
      character(len=1)                            :: unique_label=" " ! A,B, or C to indicate unqiue axis. Only used in monoclinic
      integer                                     :: unique=0         ! integer to indicate unique axis
      logical,dimension(3)                        :: obtuse=.false.   ! .true. if cell angle is obtuse. Only used in monoclinic and triclinic
      type(Eos_type)                              :: eosc             ! The common factors to all EoS in an EoS
      type(Eos_type)                              :: eosang           ! The unit cell angle information, if stored as polynomials
      integer,dimension(0:6)                      :: loaded = 0       ! 0 when absent, 1 when eos present, 2 set by symmetry, 3 when possible to calc, 4 monoclinic d_unique (set by set_cell_types)
      character(len=1),dimension(0:6,3)           :: cout = " "       ! output array for reporting PV, VT and PVT types of EoS
      character(len=30)                           :: inputlist= " "   ! List of allowed eos that can be selected, given the cell symmetry. Useful for i/o prompts
   End Type EoS_Cell_Type

   !!----
   !!---- TYPE :: AXIS_TYPE
   !!--..
   !!---- Update: 03/02/2021
   !!
   Type, public :: Axis_type
      real(kind=cp), dimension(3)  :: v=0.0        ! UVW or hkl
      character(len=1)             :: atype=' '    ! axis type, H=hkl, U=UVW
      integer                      :: Ieos=0       ! >0 if a primary axis. 0 = volume,  -1 error, if -2 axis vector in array axis
   End Type Axis_type

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
      character(len=80)                              :: Title=" "           ! Title of dataset (normally from input datafile)
      character(len=40)                              :: System=" "          ! Crystal System  (normally set by Def_Crystal_System)
      integer                                        :: N=0                 ! Number of EoS Data List
      integer, dimension(NCOL_DATA_MAX)              :: IC_Dat=0            ! Which values are input
      character(len=15)                              :: Pscale_name=" "     ! Description of the Pressure scale of data (e.g. GPa)
      character(len=15)                              :: Vscale_name=" "     ! Description of the units of volume data (e.g. A3/cell)
      character(len=15)                              :: Lscale_name=" "     ! Description of the units of linear data  (e.g. A)
      type(EoS_Data_Type), allocatable, dimension(:) :: EoSD                ! Data values
   End Type EoS_Data_List_Type

   !-------------------!
   !---- VARIABLES ----!
   !-------------------!
   logical,            public :: Err_EOS=.false.   ! Logical Variable indicating an error in CFML_EoS module
   character(len=256), public :: Err_EOS_Mess=" "  ! String containing information about the last error

   logical,            public :: Warn_EOS=.false.  ! Logical Variable indicating an error in CFML_EoS module
   character(len=256), public :: Warn_EOS_Mess=" " ! String containing information about the last error

   !---------------------------------!
   !---- Interfaces - Overloaded ----!
   !---------------------------------!
   Interface  Linear_EoS_Allowed
       Module Procedure Linear_EoS_Allowed_I
       Module Procedure Linear_EoS_Allowed_EoS
   End Interface

Contains

   !!----
   !!---- FUNCTION ALPHA_CAL
   !!----
   !!---- Calculate the alpha parameter in thermal case
   !!---- For linear case alpha is correct as 1/a da/dT
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Alpha_Cal(P, T, EoS, DeltaT) Result(Alpha)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: P       ! Pressure
      real(kind=cp),            intent(in) :: T       ! Temperature
      type(Eos_Type),           intent(in) :: EoS     ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T value
      real(kind=cp)                        :: alpha

      !---- Local Variables ----!
      integer                        :: j
      real(kind=cp)                  :: del, tt, delmin, tlimit, tr, alphaest, vol
      real(kind=cp), dimension(-2:2) :: v  ! array for calc v values

      !> Init
      alpha=0.0_cp
      Tlimit=0.0_cp

      !> PTV table
      if (eos%imodel == -1) then
         alpha=get_props_ptvtable(p, t, 0.0, eos, 'AL')
         return
      end if

      !>No thermal model
      if (eos%itherm < 0 .or. eos%itherm > N_THERM_MODELS) return

      !> Oscillator based eos (itherm=7,8), with no phase transition
      !> Alpha calculated direct from Cv based on equation 2.83 of Anderson (1995)
      !> Cv is returned in J/mol/K by direct algebra expression, and using R=8.314
      if (eos%Osc_allowed .and. eos%itran == 0) then
         vol=get_volume(p,t,eos)
         alpha=get_grun_v(vol,eos)*get_cv(p,t,eos)/k_cal(vol,t,eos,p=p)/vol

         !> scaling
         alpha=alpha*EPThermal_factor(eos)
         return
      end if

      !> Need to trap numerical problems with Kroll, Salje, Pthermal at low T
      select case(eos%itherm)
         case(0) ! no thermal parameters
            return

         case(4:6) ! Kroll, Salje, and HP Pthermal: For T < 0.05T(Einstein) alpha=0. Same for Salje but test is 0.05Tsat
            if (t < 0.05_cp*eos%params(11) ) return
      end select

      !> Numerical solutions: set step size (not used in MGD)
      del =abs(0.001_cp/eos%params(10))      ! Set step in T to get about 0.1% shift in V
      if (del > 80.0_cp) del=80.0_cp         ! otherwise for small alpha, del is too big and alpha is inaccurate
      if (present(deltaT)) del=deltaT

      select case(eos%itherm)           ! adjustment of step size
         !> Fei, HP98
         case(2,3)
            !> T ok, but do not step down into invalid region
            if (abs(eos%params(12)) > tiny(0.0_cp) ) then
               delmin=abs(t-tlimit)/2.1_cp
               if (del > delmin) del=delmin
            end if

         !> Kroll, Salje, HP and linear Thermal pressure
         case(4:6)
            delmin=(t-0.025_cp*eos%params(11))/2.0_cp     ! do not allow step in to area where alpha=0
            if (del > delmin) del=delmin                  ! ensures T at all steps is positive

         !> MGD Pthermal and q-compromise
         case(7:8)
            ! so no alpha available for estimation: changed from T+100 to T-100 to avoid going
            ! into area where large V invalid
            alphaest=(get_volume(p,t,eos)-get_volume(p,t-100._cp,eos))/get_volume(p,t-50._cp,eos)/100._cp
            del=abs(0.001_cp/alphaest)
            delmin=(t-0.025_cp*eos%params(11))/2.0_cp     ! do not allow step in to area where alpha=0
            if (del > delmin) del=delmin                  ! ensures T at all steps is positive

            ! now stop the del taking us into illegal area
            tlimit=t+2.0*del
            do
               ! search for positive K at this P
               if (get_K(p,tlimit,eos) > 0.0_cp .and. (.not. err_eos)) exit
               call init_err_eos()
               tlimit=tlimit-0.1_cp*(tlimit-t)
               if (tlimit < t) exit                    ! should never happen because P,T is valid
            end do
            del=0.4*abs(tlimit-t)
      end select

      !> Stop calculation going across a phase boundary
      if (eos%itran > 0) then
         Tr=get_transition_temperature(p,eos)
         if (transition_phase(P,T,eos) .neqv. transition_phase(P,T+2.0*del,eos)) del=abs(T-Tr)/2.1_cp
         if (transition_phase(P,T,eos) .neqv. transition_phase(P,T-2.0*del,eos)) del=abs(T-Tr)/2.1_cp
         if (del < 1.0) del=1.0
      end if

      !> Do the numerical solution
      do j=-2,2,1
         tt=t+real(j)*del              ! apply shift to temp
         v(j)=get_volume(p,tt,eos)     ! calc resulting V
      end do

      alpha=(v(-2)+8.0_cp*(v(1)-v(-1))-v(2))/(12.0_cp*del)/v(0)     ! Derivative to second order approximation

      return
   End Function Alpha_Cal

   !!--++
   !!--++ FUNCTION DERIV_PARTIAL_P_ANALYTIC
   !!--++
   !!--++ Calculates the partial derivatives of P with respect to the EoS
   !!--++ at a given v,t point, and returns  them in array td
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Function Deriv_Partial_P_Analytic(V,T,Eos) Result(TD)
      !---- Arguments ----!
      real(kind=cp),                      intent(in)  :: V       ! Volume
      real(kind=cp),                      intent(in)  :: T       ! Temperature
      type(Eos_Type),                     intent(in)  :: Eos     ! Eos Parameter
      real(kind=cp), dimension(N_EOSPAR)              :: td      ! derivatives dP/d(param)

      !---- Local Variables ----!
      real(kind=cp), dimension(N_EOSPAR) :: ev
      real(kind=cp)                      :: vv0,k0,delt, delt2,vt0,dpdvt,deltinv
      real(kind=cp)                      :: f,b,c,cv,kp,kpp
      real(kind=cp)                      :: term1,term2,term3,vterm
      real(kind=cp)                      :: a,f52,da,db,dc,nu   ! new guys
      real(kind=cp),dimension(3)         :: abc                 !Tait parammeters

      !> Init
      td=0.0_cp

      !> Local copy
      ev= EoS_to_Vec(eos)       ! Volume or linear case is covered: ev contains volume-like parameters and
                                    ! all derivatives will be calculated as volume-like derivatives
      cv=v                          ! cv has the current volume or 'a' input. It is needed for dP/dVo

      select case (eos%Itherm)
         case (0)
            vt0=eos%params(1)                ! vt0 is the V (or a) at this T and P=0

         case (1:)
            vt0=get_V0_T(t,eos)
      end select
      vv0=v/vt0                                 ! vv0 is V(P,T)/V(P,T=0), or the linear equivalent

      k0=Get_K0_T(T,eos)                     ! returns M(T) for linear
      if (eos%linear) k0=k0/3.0_cp
      kp=Get_Kp0_T(T,eos)
      if (eos%linear) kp=kp/3.0_cp
      kpp=Get_Kpp0_T(T,eos)
      if (eos%linear) kpp=kpp/4.0_cp

      !> now get the finite strain for vv0
      f=strain(vv0,eos)                      ! use strain to get finite strain from v/v0 or a/a0

      !> Using volume dimensions for remainder of calculations
      vv0=1.0_cp/vv0                               ! vv0 now V0/V for ease of use in rest of subprogram
      if (eos%linear) then
         vt0=vt0**3.0_cp           ! Vo at this T
         vv0=vv0**3.0_cp           ! V/Vo at this T
         cv=cv**3.0_cp             ! current volume at this P,T
      end if

      select case (eos%imodel)
         case (1) ! Murnaghan: validated algebra, 10/04/2013
            td(1)=vv0**(kp-1.0_cp)*k0/cv
            td(2)=(vv0**kp-1.0_cp)/kp
            td(3)=k0/kp*(vv0**kp*(log(vv0)-1.0_cp/kp) +1.0_cp/kp)
            td(4)=0.0_cp

         case (2) ! Birch-Murnaghan: new code 10-11/04/2013 Validated
            !> Assign coefficients of f in expansion of P: use K0 as that is K at P=0
            !  f was calculated above at init
            a=K0                           ! P=3(1+2f)^(5/2) . (af + bf^2 +cf^3)
            b=0.0_cp
            c=0.0_cp
            f52=(1.0_cp+2.0*f)**2.5_cp
            if (eos%iorder > 2) b=1.5_cp*K0*(kp-4.0_cp)
            if (eos%iorder == 4)c = 1.5_cp*K0*(K0*kpp + (kp-4.0_cp)*(kp-3.0_cp)+35.0_cp/9.0_cp)

            !> dP/dVo
            td(1)=f52/vt0 * ( a + (7.0_cp*a+2.0_cp*b)*f + (9.0_cp*b+3.0_cp*c)*f*f + 11.0_cp*c*f*f*f)

            !> dP/dKo:  da, db, dc are da/dparam etc for each param K0, Kp0, Kpp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eos%iorder > 2) db = 1.5_cp*(kp-4.0_cp)
            if (eos%iorder == 4)dc = 3.0_cp*K0*kpp + 1.5_cp*((kp-4.0_cp)*(kp-3.0_cp)+35.0_cp/9.0_cp)

            td(2)= 3.0*f52*f* (1.0 + db*f + dc*f*f)

            !> dP/dKp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eos%iorder > 2) db = 1.5_cp*K0
            if (eos%iorder == 4)dc = 1.5_cp*K0*(2.0_cp*kp-7.0_cp)

            td(3) = 3.0*f52* (db + dc*f) *f*f

            !> dP/dKpp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eos%iorder == 4) then
               dc = 1.5_cp*K0*K0
               td(4) = 3.0*f52*dc*f*f*f
            else
               td(4)=0.0_cp
            end if

         case (3) ! Vinet: new code and validated 11/04/2013
            if (eos%iorder == 2) then
               td(1) = k0/vt0 *(1.0_cp+f)/(1.0_cp-f)**2.0_cp
               td(2) =      3.0_cp * f /(1.0_cp-f)**2.0_cp
               td(3) = 0.0_cp
               td(4) = 0.0_cp

            else
               nu=1.5_cp*(kp-1.0)
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
            if (eos%iorder > 2) a=1.5_cp*(kp-2.0_cp)
            if (eos%iorder ==4) b=1.5_cp*(1.0_cp + K0*kpp + (kp-2.0_cp) + (kp-2.0_cp)**2.0_cp)

            !> dP/dVo
            td(1) = K0/cv * (1.0_cp  + (2.0_cp*a+3.0_cp)*f + 3.0_cp*(a+b)*f*f +3.0_cp*b*f*f*f)

            !> d(p)/d(k0)
            td(2) = 3.0_cp * vv0 * (f + a*f*f + (b + 1.5_cp*kpp*k0)*f*f*f)

            !> dp/dKp0
            da = 0.0_cp
            db = 0.0_cp
            if(eos%iorder > 2)  da=1.5_cp
            if (eos%iorder == 4)db = 3.0_cp*kp - 4.5_cp

            td(3) = 3.0_cp * k0 * vv0 * (da + db*f)*f*f

            !> dp/dKpp0
            td(4)=0.0_cp
            if (eos%iorder == 4) td(4) = 4.5_cp* k0*k0 * vv0 *f*f*f

         case (5) ! Tait: new code and validated against Murnaghan and against finite diffs 15/06/2013
            abc= get_tait(t,eos)
            vv0=1.0/vv0        ! now vv0 is v/v0

            !> dP/dVo
            td(1)= k0 * vv0**2.0_cp/cv * ( (vv0 + abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3) -1.0_cp)

            !> d(p)/d(k0)
            td(2)=abc(1)*abc(3)*(( (vv0 + abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3)) -1.0_cp)

            !> dp/dKp0
            if (eos%iorder > 2) then
               da= k0*kpp/(1.0_cp +kp+k0*kpp)**2.0_cp      ! da/dKp
               db= 1.0_cp/k0 + kpp/(1.0_cp+kp)**2.0_cp       ! db/dKp
               dc= 1.0_cp/(kp*kp + kp -k0*kpp) - ((1.0_cp)+kp+ &
                   k0*kpp)*(2.0_cp*kp +1.0_cp)/(kp*kp + kp -k0*kpp)**2.0_cp

               vterm = (vv0/abc(1) + 1.0_cp - 1.0_cp/abc(1))**(-1.0_cp/abc(3))
               term1=-1.0_cp/abc(2)/abc(2) * db * vterm
               term2= 1.0_cp/abc(2)*vterm/abc(3)/abc(3)*dc*log((vv0-1.0_cp)/abc(1) + 1.0_cp)
               term2=term2+((vv0-1.0_cp)/abc(1) + 1.0_cp)**(-1.0_cp/abc(3)-1.0_cp)/abc(2)/abc(3)/abc(1)/abc(1)*(vv0-1.0_cp)*da
               term3=db/abc(2)/abc(2)
               td(3)=term1+term2+term3
            end if

            !> dp/dKpp0
            if (abs(kpp) > tiny(0.) .and. eos%iorder == 4)then          ! if Kpp0 = 0 then the derivative is undefined
               da=-1.0_cp*k0*abc(1)/(1.0_cp+kp+k0*kpp)
               db=-1.0_cp/(1.0+kp)
               dc=(k0*(kp+1.0_cp)**2.0_cp)/(kp*kp+kp-k0*kpp)**2.0_cp

               term1= -1.0_cp/abc(2)/abc(2)*db*(vterm-1.0_cp)
               term2= 1.0_cp/abc(2)*vterm/abc(3)/abc(3)*dc*log((vv0-1.0_cp)/abc(1) + 1.0_cp)
               term3= ((vv0-1.0_cp)/abc(1) + 1.0_cp)**(-1.0_cp/abc(3)-1.0_cp)/abc(2)/abc(3)/abc(1)/abc(1)*(vv0-1.0_cp)*da
               td(4)=term1+term2+term3
            end if
      end select

      if (eos%ITherm > 0) then
         !> First do the parameters common to all thermal
         delt=t-eos%tref
         delt2=t*t-eos%tref*eos%tref
         deltinv=0.0_cp
         if (eos%tref > 1.0_cp) deltinv=1.0_cp/t-1.0_cp/eos%tref

         !> Adjust dp/dv(0,t) to dp/dv(0,0):
         dpdvt=td(1)
         td(1)=td(1)*vt0/ev(1)

         !> dp/dk(Tref,0) = dp/dk(T,0)

         !> d(k)/d(t) ...(the temperature dependence of k)
         td(5)=td(2)*delt

         !> Now do the specific derivatives for each equation, where possible. Program uses numeric derivs anyway
         select case(eos%itherm)     !  dp/dalpha=dp/dV(0T) . dV(0T)/dalpha
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
      if (eos%linear) then
         td(1) =td(1) *(3.0_cp*eos%params(1)**2.0_cp)
         td(2:5)=td(2:5)/3.0_cp

         select case(eos%itherm)                 ! thermal expansion terms
            case(1:2)                               ! Berman, Fei, alpha terms
               td(10:12)=td(10:12)*3.0_cp

            case(4)                                 ! Holland-Powell alpha
               td(10)=td(10)*3.0_cp
         end select                                 !other thermal equations have nothing to convert
      end if

      return
   End Function Deriv_Partial_P_Analytic

   !!--++
   !!--++ FUNCTION DERIV_PARTIAL_P_NUMERIC
   !!--++
   !!--++   Calculates the partial derivatives of P with respect to the EoS
   !!--++   at a given property and t point, and returns  them in array td
   !!--++
   !!--++
   !!--++ Date: 17/03/2017
   !!
   Function Deriv_Partial_P_Numeric(X1, X2, Eos, xtype, calc) Result(TD)
      !---- Arguments ----!
      real(kind=cp),                      intent(in) :: X1,X2   ! The two parameter values (eg V and T)
      type(Eos_Type),                     intent(in) :: Eos     ! Eos Parameters
      integer,         optional,          intent(in) :: xtype   ! =0 for V,T input, =1 for Kt,T =2 for Ks,T
      character(len=*),optional,          intent(in) :: calc    ! 'all' if all derivs required.
                                                                ! If absent, then just params with iref=1 calculated
      real(kind=cp), dimension(n_eospar)             :: td      ! derivatives dP/d(param)


      !---- Local Variables ----!
      type(Eos_Type)                 :: Eost                 ! Eos Parameter local copy
      real(kind=cp), dimension(-2:2) :: p                    ! array for calc p values
      real(kind=cp)                  :: delfactor,del,d_prev,delmin
      integer                        :: i,j,icycle,itype
      logical                        :: warn            ! local warn flag
      logical                        :: iall            ! .true. if all derivs required


      !> Init
      TD=0.0_cp
      call Init_err_eos()

      itype=0
      if(present(xtype))itype=xtype

      iall=.false.      ! default is calc params with iref=1
      if (present(calc) )then
         if (index(u_case(adjustl(calc)),'ALL') > 0) iall=.true.
      end if

      eost=eos
      call Set_Eos_Use(eost)

      warn=.false.

      !> Check
      if (itype < 0 .or. itype > n_data_types)then
         err_eos=.true.
         Err_EoS_Mess='No type set for deriv_partial_p'
         return
      end if

      !> Set the inital shift factor (fraction of parameter value)
      delfactor=0.01_cp

      do i=1,N_EOSPAR
         if (eos%iref(i) == 1 .or. iall) then    ! only refined parameters
            del=delfactor*eos%params(i)         ! the initial shift estimate
            delmin=delfactor/eos%factor(i)      ! scale required min shift by print factors
            if (abs(del) < delmin) del=delmin      ! trap param values that are zero
            icycle=0                               ! iteration count
            d_prev=0.0_cp

            iter:do                                ! top of loop over iterations
               do j=-2,2,1
                  eost=eos                                             !reset eos params
                  eost%params(i)=eost%params(i)+float(j)*del              ! apply shift to a parameter
                  p(j)=get_pressure_x(x1,x2,eost,itype)                  ! calc resulting P
               end do

               td(i)=(p(-2)+8.0_cp*(p(1)-p(-1))-p(2))/(12.0_cp*del)       ! derivative to second order approximation

               !write(6,'(''  Param # '',i2,'' Cycle '',i2,'': deriv = '',ES14.6,'' delp = '',f5.3,'' del= '',f9.4)')i,icycle,td(i),p(2)-p(-2),del

               ! to trap problems
               if (err_eos)then
                  td(i)=d_prev         ! previous cycle value
                  call init_err_eos()  ! clear errors
                  warn=.true.          ! warning flag
                  exit iter
               end if

               if (abs(td(i)) < 1.0E-8) exit iter                         ! zero deriv
               if (icycle > 0 .and. &
                            abs(d_prev-td(i))/td(i) < 1.0E-4) exit iter    ! deriv converged to 1 part in 10^4

               d_prev=td(i)                ! store last deriv value
               del=2.0_cp*del              ! increase the shift
               icycle=icycle+1
               if (icycle > 5) exit iter    ! Do not allow 2*shift to exceed 64% of param value
            end do iter
         end if
      end do

      if (warn) then
         warn_eos=.true.
         warn_eos_mess='Error calculating some derivatives in Least squares'
      end if

      !> no need to fix derivatives for linear eos by this method

      return
   End Function Deriv_Partial_P_Numeric

   !!--++
   !!--++ FUNCTION DERIV_PARTIAL_P_SCALES
   !!--++
   !!--++ Calculates the partial derivative of P with respect to the EoS scale parameters only
   !!--++ at a given property and t point, and it in td
   !!--++
   !!--++
   !!--++ Date: 24/03/2020
   !!
   Function Deriv_Partial_P_Scales(V, T, EoS, Xtype,Igp) Result(TD)
      !---- Arguments ----!
      real(kind=cp),                      intent(in) :: V,T     ! The two parameter values (eg V and T, used to make coding clearer)
      type(Eos_Type),                     intent(in) :: Eos     ! Eos Parameters
      integer,                            intent(in) :: xtype   ! =0 for V,T input, =1 for Kt,T =2 for Ks,T
      integer,                            intent(in) :: igp     ! igp is the group number of the V,so implies scale factor in params(50+igp)
      real(kind=cp)                                  :: td      ! derivative dP/d(param)

      !---- Local Variables ----!
      integer                        :: ip                   ! param number=50+igp
      integer                        :: icycle,j
      real(kind=cp), dimension(-2:2) :: p                    ! array for calc p values
      real(kind=cp)                  :: del,d_prev,vol
      logical                        :: warn            ! local warn flag


      !>Init
      td=0.0_cp
      if (igp < 1 .or. igp > 9)return

      ip=igp+50                         !the param number
      if (eos%iref(ip) /= 1 )return     ! this scale not refined

      !> init iteration
      del=0.01*eos%params(ip)
      icycle=0                               ! iteration count
      d_prev=0.0_cp

      !Note: the input V is on the same scale as the eospar
      iter:do                                ! top of loop over iterations
         do j=-2,2,1
            vol=V*eos%params(ip)/(eos%params(ip)+float(j)*del)       ! apply shift to the Vol by a shift to the scale
            p(j)=get_pressure_x(Vol,T,eos,xtype)                     ! calc resulting P
         end do

         td=(p(-2)+8.0_cp*(p(1)-p(-1))-p(2))/(12.0_cp*del)       ! derivative to second order approximation

         !write(6,'(''  Param # '',i2,'' Cycle '',i2,'': deriv = '',ES14.6,'' delp = '',f5.3,'' del= '',f9.4)')i,icycle,td(i),p(2)-p(-2),del

         ! to trap problems
         if (err_eos)then
            td=d_prev     ! previous cycle value
            call init_err_eos()  ! clear errors
            warn=.true.       ! warning flag
            exit iter
         end if

         if (abs(td) < 1.0E-8) exit iter                         ! zero deriv
         if (icycle > 0 .and. &
             abs(d_prev-td)/td < 1.0E-4) exit iter    ! deriv converged to 1 part in 10^4

         d_prev=td                ! store last deriv value
         del=2.0_cp*del              ! increase the shift
         icycle=icycle+1
         if (icycle > 5) exit iter    ! Do not allow 2*shift to exceed 64% of param value
      end do iter

      if (warn)then
         warn_eos=.true.
         warn_eos_mess='Error calculating some derivatives of scales in Least squares'
      end if

      return
   End Function Deriv_Partial_P_Scales

   !!----
   !!---- FUNCTION DERIV_PARTIAL_P
   !!----
   !!---- Calculate Partial derivates of Pressure respect to Params
   !!----
   !!---- Date: 21/03/2017
   !!
   Function Deriv_Partial_P(V, T, EoS, Xtype, Calc) Result(TD)
      !---- Arguments ----!
      real(kind=cp),                      intent(in)  :: V       ! Volume
      real(kind=cp),                      intent(in)  :: T       ! Temperature
      type(Eos_Type),                     intent(in)  :: EoS     ! Eos Parameter
      integer,optional,                   intent(in)  :: xtype   ! =0 for V,T input, =1 for Kt,T =2 for Ks,T
      character(len=*),optional,          intent(in)  :: calc    ! 'all' if all derivs required. If absent, then just params with iref=1 calculated
      real(kind=cp), dimension(N_EOSPAR)              :: td      ! derivatives dP/d(param)

      !---- Local Variables ----!
      integer                            :: itype
      real(kind=cp), dimension(N_EOSPAR) :: tda,tdn                ! analytic and numeric derivatives
      character(len=10)                  :: cstring

      !> Init
      td=0.0_cp

      itype=0
      if (present(xtype))itype=xtype
      cstring='ref'
      if (present(calc))cstring=u_case(adjustl(calc))

      !> Calculate derivatives by both methods if possible: correct values are returned for linear
      TDN=Deriv_Partial_P_Numeric(V,T,Eos,itype,cstring)

      !> Default to numeric, because they are always available:
      td(1:N_EOSPAR)=tdn(1:N_EOSPAR)

      if (itype == 0 .and. .not. Eos%pthermaleos)then
         TDA=Deriv_Partial_P_Analytic(V,T,Eos)
         if (eos%itran ==0 .and. eos%imodel /=6) then ! imodel=6 is APL, not yet coded
            td(1:4)=tda(1:4)                          ! analytic for Vo and moduli terms because these are exact even at small P
         end if
      end if

      return
   End Function Deriv_Partial_P

   !!----
   !!---- FUNCTION DKDT_CAL
   !!----
   !!---- Calculate the derivative dK/dt (or dM/dt in the linear case) at P and T
   !!----
   !!---- Date: 17/07/2015
   !!
   Function dKdT_Cal(P, T, EoS, DeltaT) Result(dKdT)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: P       ! Pressure
      real(kind=cp),            intent(in) :: T       ! Temperature
      type(Eos_Type),           intent(in) :: EoS     ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T
      real(kind=cp)                        :: dKdT

      !---- Local Variables ----!
      integer                      :: j
      real(kind=cp)                :: del,tlimit,tcal,Ttr !,vlimit
      real(kind=cp),dimension(-2:2):: kpt


      !> Init: set step size for numerical differentiation
      del=30.0_cp                        ! good number for accuracy
      if (present(deltat)) del=deltat
      if (t-2.0*del < 0.0_cp) del=t/2.1  ! prevents going to negative T

      !> Code to prevent crossing a phase boundary
      if (eos%itran > 0) then
         Ttr = Get_Transition_Temperature(P,eos)
         if (transition_phase(P,T+2.0*del,eos) .neqv. transition_phase(P,T,eos)) del=0.4*abs(Ttr-T)
         if (transition_phase(P,T-2.0*del,eos) .neqv. transition_phase(P,T,eos)) del=0.4*abs(Ttr-T)
      end if

      !> Code to stop some Pthermal EoS going into illegal large volume above T
      if (eos%itherm == 7 .or. eos%itherm == 8) then
         tlimit=t+2.0*del
         do                                        ! search for positive K at this P
            if (get_K(p,tlimit,eos) > 0.0_cp .and. (.not. err_eos)) exit
            call init_err_eos()
            tlimit=tlimit-0.1_cp*(tlimit-t)
            if (tlimit < t)exit                    ! should never happen because P,T is valid
         end do
         del=0.4*abs(tlimit-t)
      end if

      !> Trap close to zero K
      if (t < 1.0_cp) then
         dKdT=Get_K(P,1.0,eos)-Get_K(P,0.0,eos)

      else
         do j=-2,2,1
            tcal=t+real(j)*del              ! apply shift to t
            kpt(j)=Get_K(P,tcal,eos)        ! calc resulting K
         end do
         dKdT=(kpt(-2)+8.0_cp*(kpt(1)-kpt(-1))-kpt(2))/(12.0_cp*del)     ! Derivative to second order approximation
      end if

      !> No linear conversion is required because get_K returns values for "linear Kp" = Mp,
      !> so kppc is already dMp/dP = Mpp

      return
   End Function dKdT_Cal

   !!----
   !!---- FUNCTION EOS_CAL
   !!----
   !!---- Returns elastic properties (not the parameter values) at this P,T for EoS
   !!----
   !!---- Date: 17/07/2015
   !!
   Function EoS_Cal(P,T,EoS) Result(Parvals)
      !---- Arguments ----!
      real(kind=cp),                intent(in)  :: P       ! Pressure
      real(kind=cp),                intent(in)  :: T       ! Temperature
      type(Eos_Type),               intent(in)  :: EoS     ! Eos Parameter
      real(kind=cp), dimension(6)               :: Parvals ! Output parameter values

      !> Init
      parvals=0.0_cp

      call physical_check(eos,Pin=p,Tin=t)           ! produce warnings based on P,T
      if (Err_eos)return

      parvals(1)=get_volume(p,t,eos)
      parvals(2)=k_cal(parvals(1),t,eos,P=p)
      parvals(3)=kp_cal(parvals(1),t,eos,P=p)
      parvals(4)=kpp_cal(parvals(1),t,eos)
      parvals(5)=dKdT_cal(p,t,eos)           ! dK/dT at this P,T
      parvals(6)=alpha_cal(p,t,eos)          ! 1/V.dV/dT at this T

      return
   End Function EoS_Cal

   !!----
   !!---- FUNCTION EOS_CAL_ESD
   !!----
   !!---- Returns esd's of the elastic properties (not the parameter values) at this P and T for EoS
   !!----
   !!---- Date: 17/07/2015
   !!
   Function EoS_Cal_Esd(P,T,EoS) Result(Esd)
      !---- Arguments ----!
      real(kind=cp),                intent(in)  :: P       ! Pressure
      real(kind=cp),                intent(in)  :: T       ! Temperature
      type(Eos_Type),               intent(in)  :: EoS     ! Eos Parameter
      real(kind=cp),  dimension(6)              :: Esd     ! Output esd values

      !---- Local Variables ----!
      real(kind=cp),dimension(N_EOSPAR) :: esdfull

      !> Init
      esd=0.0_cp
      esdfull=0.0_cp

      !> calculate parameter esd's transformed to this P,T
      esdfull= transform_Esd(P,T,EoS)
      esd(1:5)=esdfull(1:5)
      esd(6)=esdfull(10)

      !> make adjustment for only using alpha0
      select case(EoS%itherm)
          case(5)               ! Salje
            esd(6)=esd(6)/3.0_cp            ! alpha is 1/3 of p1
      end select

      return
   End Function EoS_Cal_Esd

   !!--++
   !!--++ FUNCTION EOS_TO_VEC
   !!--++
   !!--++ Copy parameters from EoS type to a vector
   !!--++
   !!--++ Date: 28/02/2013
   !!
   Function EoS_to_Vec(Eos) Result(Vec)
      !---- Arguments ----!
      type(EoS_Type),              intent(in)  :: Eos
      real(kind=cp), dimension(N_EOSPAR)       :: Vec

      !> Init
      vec=0.0_cp

      !> Copy everything 1 to 1 as default
      vec=eos%params

      !> adjust any linearr values as necessary
      if (eos%linear) then
         vec(1)=eos%params(1)**3.0_cp
         vec(2:4)=eos%params(2:4)/3.0_cp

         select case(eos%icross)
            case (1)
               vec(8)=eos%params(8)/3.0_cp

            case (2)
               vec(8:9)=eos%params(8:9)
         end select

         select case(eos%itherm)                 ! thermal expansion terms
            case(1:3)
               vec(10:12)=eos%params(10:12)*3.0_cp

            case(4,5,6,8)
               vec(10)=eos%params(10)*3.0_cp
               vec(11)=eos%params(11)
         end select

         !> phase transition: no changes required for linear

      end if

      return
   End Function EoS_to_Vec

   !!--++
   !!--++ FUNCTION EPTHERMAL_FACTOR
   !!--++
   !!--++  Calculate Scale factor for E(thermal) to P(thermal) on basis of units for P and V
   !!--++  Scale should be used to multiply the Pthermal calculated from eos parameters
   !!--++
   !!--++ Date: 03/09/2020
   !!
   Function EPthermal_factor(Eos) Result(scale)
      !---- Arguments ----!
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameters
      real(kind=cp)              :: scale

      !---- Local Variables ----!
      character(len=len(eos%vscale_name)) :: vname

      !>init
      scale=1.0_cp

      !>if the thermal energy was from EthDebye or EthEinstein, it is in J/mol pth
      !>Then if V in m3/mol  Pth=Eth/V is in J/m3=Pa

      !> Pressure scales
      if (index(U_case(eos%pscale_name),'GPA') > 0)  scale=1.0E-9
      if (index(U_case(eos%pscale_name),'KBAR') > 0) scale=1.0E-8

      !> Volume
      vname=adjustl(U_case(eos%vscale_name))
      if (len_trim(vname) == 0)return

      !> test for cm3/mol or equivalent
      if (index(vname,'CM') > 0 .and. index(vname,'3') > 0 .and. index(vname,'MOL') > 0)scale=scale*1.0E+6

      return
   End Function EPthermal_factor

   !!--++
   !!--++ FUNCTION ETHDEBYE
   !!--++
   !!--++  Calculates the Debye thermal Energy in Jmol(-1)
   !!--++  because R is given in Jmol(-1)K(-1)
   !!--++
   !!--++ Date: 18/11/2015
   !!
   Function EthDebye(T, Theta, Natom) Result(Eth)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      real(kind=cp),  intent(in) :: Theta   ! Debye T
      real(kind=cp),  intent(in) :: Natom   ! Number of atoms in formula unit
      real(kind=cp)              :: Eth

      !---- Local Variables ----!
      real(kind=8)  :: x

      if (T < 0.1) then
         Eth=0.0_cp

      else
         x=theta/t
         Eth=debye(3,x)
         Eth=3.0_cp*Natom*8.314_cp*T*Eth
      end if

      if (err_Mathgen)then
         err_eos=.true.
         err_eos_mess=trim(err_Mathgen_mess)
      end if

      return
   End Function EthDebye

   !!--++
   !!--++ FUNCTION ETHEINSTEIN
   !!--++
   !!--++  Calculates the thermal Energy in Jmol(-1) of N Einstein oscillators
   !!--++  because R is given in Jmol(-1)K(-1)
   !!--++
   !!--++ Date: 3/2020
   !!
   Function EthEinstein(T, Theta, Natom) Result(Eth)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      real(kind=cp),  intent(in) :: Theta   ! Einstein T
      real(kind=cp),  intent(in) :: Natom   ! Number of atoms in formula unit
      real(kind=cp)              :: Eth

      !---- Local Variables ----!

      if (T < 0.1) then
         Eth=0.0_cp
      else
         Eth=3.0_cp*Natom*8.314_cp*Theta/(exp(Theta/T)-1.0_cp)
      end if

      return
   End Function EthEinstein

   !!--++
   !!--++ FUNCTION GET_ALPHA_AXIS
   !!--++
   !!--++ Returns the value of alpha of principal axis (ieos) in unit cell
   !!--++ in cell_eos at P,T
   !!--++ Call this Function directly when the calling routine  knows that the direction is a principal axis
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Alpha_Axis(P, T, Cell_EoS, Ieos) Result(Alpha)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: P, T
      type(eos_cell_type),intent(in)  :: Cell_EoS
      integer,            intent(in)  :: Ieos      !axis indicator, as in axis_type%ieos
      real(kind=cp)                   :: Alpha     !returned thermal expansion

      !---- Local Variables ----!

      !> init
      Alpha=0.0_cp

      select case(cell_eos%loaded(ieos))
         case(1)
            Alpha=alpha_cal(p,t,cell_eos%eos(ieos))

         case(2) ! sym equiv. Always uses eos(1) for a-axis
            Alpha=alpha_cal(p,t,cell_eos%eos(1))

         case(3)
            Alpha=get_alpha_third(p,T,cell_eos,ieos)

         case(4)
            alpha=alpha_cal(p,t,cell_eos%eos(cell_eos%unique))
      end select

      return
   End Function Get_Alpha_Axis

   !!----
   !!---- FUNCTION GET_ALPHA_CELL
   !!----
   !!---- Returns the value of alpha of any axis in unit cell in cell_eos at P,T
   !!----
   !!---- Call this Function when the calling routine does not know if the direction is a principal axis or not
   !!---- If a principal direction is requested, only axis%ieos is required
   !!---- axis%v and axis%atype only used if axis%ieos=-2
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Alpha_Cell(P, T, Cell_EoS, Axis) Result(Alpha)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: alpha !returned thermal expansion

      !---- Local Variables ----!

      !> init
      alpha=0.0_cp

      select case(axis%ieos)
         case(0:6)
            !> principal direction for which eos exists, or can be calculated
            alpha=get_alpha_axis(p,t,cell_eos,axis%ieos)

         case(-2)   !general direction
            alpha=get_alpha_general(p,T,cell_eos,axis)
      end select

      return
   End Function Get_Alpha_Cell

   !!--++
   !!--++ FUNCTION GET_ALPHA_GENERAL
   !!--++
   !!--++ Returns the value of alpha of any axis in unit cell in cell_eos at P,T
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Alpha_General(P, T, Cell_Eos, Axis) Result(Alpha)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: alpha

      !---- Local Variables ----!
      integer         :: i
      real(kind=cp)   :: tstep,tcal

      !> for spline
      integer, parameter             :: NSTEP=21   !must be odd
      integer                        :: imid
      real(kind=cp), dimension(NSTEP):: x,y,d2y,dy

      tstep=20.
      tcal=t-int(NSTEP/2)*tstep
      do i=1,NSTEP
         x(i)=tcal
         y(i)=get_Volume_general(P,Tcal,cell_eos,axis)
         tcal=tcal+tstep
      end do

      call Second_Derivative(x, y, NSTEP, d2y)
      call First_Derivative(x, y, NSTEP, d2y, dy)

      imid=int(NSTEP/2) + 1
      alpha=dy(imid)/y(imid)

      return
   End Function Get_Alpha_General

   !!--++
   !!--++ FUNCTION GET_ALPHA_THIRD
   !!--++
   !!--++ Returns the value of alpha of a principal axis ieos in unit cell in cell_eos at P,T
   !!--++ when it can be calculated from others
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Alpha_Third(P, T, Cell_Eos, Ieos) Result(Alpha)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: p,T
      type(eos_cell_type),intent(in) :: cell_eos
      integer,            intent(in) :: ieos     ! the alpha of the axis (1,2,3) or V (0) to be calculated
      real(kind=cp)                  :: alpha

      !---- Local Variables ----!
      integer         :: i
      real(kind=cp)   :: alpha_ang

      !> init
      alpha=0._cp

      !> safety check: if mono or triclinic, should only be called if angle poly used
      if (U_case(cell_eos%system(1:4)) == 'TRIC' .or. U_case(cell_eos%system(1:3)) == 'MONO')then
         if (cell_eos%eosang%iangle == 0)then
            err_eos=.true.
            err_eos_mess='Get_Alpha_Third called for mono or triclinic, without angle poly set'
         end if
      end if

      select case(U_case(cell_eos%system(1:4)))
         case('TRIC','MONO','ORTH')
            alpha_ang=Get_Angle_Volfactor_Deriv(P,T,cell_eos,'T')/Get_Angle_Volfactor(P,T,cell_eos)  !1/A dA/dT

            select case(ieos)
               case(0)
                  alpha=Alpha_Cal(P,T,cell_eos%eos(1))+ Alpha_Cal(P,T,cell_eos%eos(2)) + &
                        Alpha_Cal(P,T,cell_eos%eos(3)) + alpha_ang

               case default
                  alpha=Alpha_Cal(P,T,cell_eos%eos(0))
                  do i=1,3
                     if (i == ieos)cycle
                     alpha=alpha- Alpha_Cal(P,T,cell_eos%eos(i))
                  end do
                  alpha=alpha-alpha_ang
            end select

         case('TRIG','HEXA','TETR')
            select case(ieos)
               case(0)     ! calc volumes from a and c
                  alpha= 2.0_cp*Alpha_Cal(P,T,cell_eos%eos(1))+Alpha_Cal(P,T,cell_eos%eos(3))

               case(1)     ! a from V and c
                  alpha= (Alpha_Cal(P,T,cell_eos%eos(0))-Alpha_Cal(P,T,cell_eos%eos(3)))/2.0_cp

               case(3)     ! c from a and V
                  alpha=  Alpha_Cal(P,T,cell_eos%eos(0))-2.0_cp*Alpha_Cal(P,T,cell_eos%eos(1))
            end select

         case('CUBI','ISOT')
            select case(ieos)
               case(0)     ! calc volume from a
                  alpha=Alpha_Cal(P,T,cell_eos%eos(1))*3.0_cp

               case(1,2,3)     ! a,b, or c from V
                  alpha=Alpha_Cal(P,T,cell_eos%eos(0))/3.0_cp
            end select
      end select

      return
   End Function Get_Alpha_Third

   !!----
   !!---- FUNCTION GET_ANGLE_DERIV
   !!----
   !!----
   !!---- Date: 04/02/2021
   !!
   Function Get_Angle_Deriv(P, T, Cell_Eos, Ia, RealAng, Dx) Result(Da)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos     ! The cell eos
      integer,            intent(in)  :: ia           ! The angle number
      logical,            intent(in)  :: realang      !.true. for real angle, .false. for recip
      character(len=1),   intent(in)  :: dx           ! either T or P
      real(kind=cp)                   :: da           ! the resulting angle derivative w.r.t. variable dx

      !---- Local Variables ----!

      !> init
      da=0.0_cp

      !> check that ia is valid and calculations required because da /= 0
      select case(U_case(cell_eos%system(1:4)))
         case('MONO')
            if (ia /=cell_eos%unique)return

         case('TRIC')
            if (ia < 1 .and. ia > 3)return

         case default
            return
      end select

      !> calculate: outer loop over angles
      if (realang .and. cell_eos%eosang%iangle > 0)then
         ! polynomial model for angles
         da=get_angle_poly_deriv(p,t,cell_eos%eosang,ia,dx)
      else
         !all other cases
         da=get_angle_eos_deriv(p,t,cell_eos,ia,realang,dx)
      end if

      return
   End Function Get_Angle_Deriv

   !!--++
   !!--++ FUNCTION GET_ANGLE_EOS_DERIV
   !!--++
   !!--++ Date: 04/02/2021
   !!
   Function Get_Angle_Eos_Deriv(Pin, Tin, Cell_EoS, Ia, RealAng, Xl) Result(D)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: pin,tin
      type(eos_cell_type),intent(in)  :: cell_eos
      integer,            intent(in)  :: ia       ! The angle number
      logical,            intent(in)  :: realang  !.true. for real angle, .false. for recip
      character(len=1),   intent(in)  :: xl       ! either T or P

      !---- Local Variables ----!
      integer           :: i
      real(kind=cp)     :: d        !the derivative d(angle_i)/dx in radians per x
      real(kind=cp)     :: p,t

      !> for spline
      integer,parameter              :: nstep=11   !must be odd
      integer                        :: imid
      real(kind=cp),dimension(nstep) :: x,y,d2y,dy
      type(crystal_cell_type)        :: c

      !>init
      d=0.0_cp
      imid=int(nstep/2) + 1

      do i=1,nstep
          if (U_case(xl) == 'P')then
             T=Tin
             P=Pin+(i-imid)*0.1_cp       !step P in 0.1
             x(i)=P

          else
             T=tin+(i-imid)*10.0_cp
             P=Pin
             x(i)=T
          end if

          c= get_params_cell(P,T,cell_eos)
          if (realang)then
             y(i)=c%ang(ia)
          else
             y(i)=c%rang(ia)
          end if
      end do

      call Second_Derivative(x, y, nstep, d2y)
      call First_Derivative(x, y, nstep, d2y, dy)

      d=dy(imid)*to_rad

      return
   End Function Get_Angle_Eos_Deriv

   !!--++
   !!--++ FUNCTION GET_ANGLE_POLY
   !!--++
   !!--++ Date: 04/02/2021
   !!
   Function Get_Angle_Poly(P, T, EoS, Ia) Result(Ang)
      !---- Arguments ----!
      real(kind=cp), intent(in) :: p,t
      type(eos_type),intent(in) :: eos
      integer,       intent(in) :: ia            ! The angle number

      !---- Local Variables ----!
      integer           :: i
      real(kind=cp)     :: ang,dt

      if (eos%angpoly(ia,0,1) < tiny(0.0_cp))then
         ang=90.0_cp
         return
      else
         ang=eos%angpoly(ia,0,1)
      end if

      dt=t-eos%tref
      do i=1,N_ANGPOLY
         ang=ang+eos%angpoly(ia,1,i)*p**i
         ang=ang+eos%angpoly(ia,2,i)*dt**i
      end do

      ang=ang+eos%angpoly(ia,3,1)*p*dt
      ang=ang+eos%angpoly(ia,3,2)*p*p*dt
      ang=ang+eos%angpoly(ia,3,3)*p*dt*dt

      if (ang < 0.0_cp .or. ang > 180.0_cp)then
         err_eos=.true.
         err_eos_mess='Angle polynomial predicted cell angle <0 or >180: reset to 90deg'
         ang=90._cp
      end if

      return
   End Function Get_Angle_Poly

   !!--++
   !!--++ FUNCTION GET_ANGLE_POLY_DERIV
   !!--++
   !!--++ Date: 04/02/2021
   !!
   Function Get_Angle_Poly_Deriv(P, T, EoS, Ia, X) Result(D)
      !---- Arguments ----!
      real(kind=cp),   intent(in) :: p,t
      type(eos_type),  intent(in) :: eos
      integer,         intent(in) :: ia ! The angle number
      character(len=1),intent(in) :: x  ! either T or P

      !---- Local Variables ----!
      integer           :: i
      real(kind=cp)     :: d  !the derivative d(angle_i)/dx in radians per x
      real(kind=cp)     :: dt !T-Tref

      !>init
      d=0._cp
      dt=t-eos%tref

      select case(U_case(x))
         case('P')
            do i=1,N_ANGPOLY
               d=d+eos%angpoly(ia,1,i)*i*p**(i-1)
            end do
            d=d+eos%angpoly(ia,3,1)*dt
            d=d+eos%angpoly(ia,3,2)*2.0_cp*p*dt
            d=d+eos%angpoly(ia,3,3)*dt*dt

         case('T')
            do i=1,N_ANGPOLY
               d=d+eos%angpoly(ia,2,i)*i*dt**(i-1)
            end do
            d=d+eos%angpoly(ia,3,1)*p
            d=d+eos%angpoly(ia,3,2)*p*p
            d=d+eos%angpoly(ia,3,3)*2.0_cp*p*dt

      end select
      d=d*to_rad

      return
   End Function Get_Angle_Poly_Deriv

   !!--++
   !!--++ FUNCTION GET_ANGLE_VOLFACTOR
   !!--++
   !!--++ Date: 04/02/2021
   !!
   Function Get_Angle_VolFactor(P, T, E) Result(vf)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: p,T
      type(eos_cell_type),intent(in) :: e
      real(kind=cp)                  :: vf !returned volume  factor

      !---- Local Variables ----!
      real(kind=cp)   :: aprod, cosang
      integer         :: i

      select case(U_case(e%system(1:4)))
         case default
            vf=1.0_cp

         case('MONO')
            i=2     !default b-unique
            if (e%unique > 0 .and. e%unique < 4)i=e%unique
            vf=sind(get_angle_poly(P,T,e%eosang,i))

         case('TRIC')
            aprod=2._cp
            vf=1.0_cp
            do i=1,3
               cosang=cosd(get_angle_poly(P,T,e%eosang,i) )
               vf=vf-cosang**2._cp
               aprod=aprod*cosang
            end do
            vf=vf+aprod
            if (vf > tiny(0.0_cp))then
               vf=sqrt(vf)
            else
               vf=1.0_cp
            end if
      end select

      return
   End Function Get_Angle_VolFactor

   !!--++
   !!--++ FUNCTION GET_ANGLE_VOLFACTOR_DERIV
   !!--++
   !!--++ Date: 04/02/2021
   !!
   Function Get_Angle_VolFactor_Deriv(P, T, E, X) Result(D)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: p,T
      type(eos_cell_type),intent(in) :: e
      character(len=1),   intent(in) :: x        !P or T
      real(kind=cp)                  :: d !returned volume  factor derivative

      !---- Local Variables ----!
      real(kind=cp)                 ::  pi,ti,del
      real(kind=cp),dimension(-2:2) :: a
      integer                       :: ia,j


      select case(U_case(e%system(1:4)))
         case default
            d=0.0_cp

         case('MONO')    !calculates cos(beta). d(beta)/dx
            ia=2     !default b-unique
            if (e%unique > 0 .and. e%unique < 4)ia=e%unique

            if (U_case(x) =='P')then
               Ti=T
               del=0.1
               do j=-2,2,1
                  pi=p+real(j)*del
                  a(j)=get_angle_poly(Pi,Ti,e%eosang,ia)
               end do

            else
               Pi=P
               del=10._cp
               do j=-2,2,1
                  Ti=T+real(j)*del
                  a(j)=get_angle_poly(Pi,Ti,e%eosang,ia)
               end do
            end if

            d=(a(-2)+8.0_cp*(a(1)-a(-1))-a(2))/(12.0_cp*del)     ! Derivative to second order approximation
            d=cosd(get_angle_poly(Pi,Ti,e%eosang,ia))*d*to_rad

         case('TRIC')
            !> calculate Volfactor as function of x, then direct deriv
            if (U_case(x) =='P')then
               Ti=T
               del=0.1
               do j=-2,2,1
                  pi=p+real(j)*del
                  a(j)=get_angle_volfactor(Pi,Ti,e)
               end do

            else
               Pi=P
               del=10._cp
               do j=-2,2,1
                  Ti=T+real(j)*del
                  a(j)=get_angle_volfactor(Pi,Ti,e)
               end do
            end if

            d=(a(-2)+8.0_cp*(a(1)-a(-1))-a(2))/(12.0_cp*del)     ! Derivative to second order approximation

      end select

      return
   End Function Get_Angle_VolFactor_Deriv

   !!--++
   !!--++ FUNCTION GET_ANGLE_VOLFACTOR_DERIV2
   !!--++
   !!--++ Date: 04/02/2021
   !!
   Function Get_Angle_Volfactor_Deriv2(P, T, E, X) Result(D)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: p,T
      type(eos_cell_type),intent(in) :: e
      character(len=*),   intent(in) :: x        !P or T or PT
      real(kind=cp)                  :: d !returned volume  factor 2nd derivative

      !---- Local Variables ----!
      real(kind=cp)                 ::  pi,ti,del
      real(kind=cp),dimension(-2:2) :: a
      integer                       :: j

      !> Init
      d=0.0_cp
      if (index(U_case(e%system),'TRIC') == 0 .and. index(U_case(e%system),'MONO') == 0) return

      !> calculate Volfactor as function of x, then direct 2nd deriv
      if (index(U_case(x),'P') > 0 .and. index(U_case(x),'T') > 0)then
         !> cross derivative
         del=40._cp
         do j=-2,2,1
            ti=t+real(j)*del
            a(j)=get_angle_volfactor_deriv(P,Ti,e,'P')
         end do

      else if(index(U_case(x),'P') > 0 .and. index(U_case(x),'T') == 0)then
         Ti=T
         del=0.4
         do j=-2,2,1
            pi=p+real(j)*del
            a(j)=get_angle_volfactor_deriv(Pi,Ti,e,'P')
         end do

      else
         Pi=P
         del=40._cp
         do j=-2,2,1
            Ti=T+real(j)*del
            a(j)=get_angle_volfactor_deriv(Pi,Ti,e,'T')
         end do
      end if

      d=(a(-2)+8.0_cp*(a(1)-a(-1))-a(2))/(12.0_cp*del)     ! Derivative to second order approximation

      return
   End Function Get_Angle_VolFactor_Deriv2

   !!--++
   !!--++ FUNCTION GET_APL
   !!--++
   !!--++ Returns a,b,c APL parameters and their derivatives in vectors.
   !!--++
   !!--++ output values:
   !!--++      a(1,j) are a,b,c
   !!--++      a(2,j) are first derivs
   !!--++      a(3,j) are second derivs of the a,b,c
   !!--++
   !!--++ Explicit values of VV0 etc used as input, because this depends on thermal model etc
   !!--++ not just on the PV model parameters
   !!--++
   !!--++ Date: 20/11/2019 RJA
   !!
   Function Get_APL(VV0, V0, K0, Kp, Kpp, Z, Iorder) Result(A)
      !---- Arguments ----!
      real(kind=cp),               intent(in)  :: VV0, V0, K0, Kp, Kpp, Z   !input parameters: VV0 is V/V0
      integer,                     intent(in)  :: iorder
      real(kind=cp),dimension(3,3)             :: a

      !---- Local Variables ----!
      real(kind=cp) :: x,c0,c2,c3,pFG0

      !> init
      a=0.0_cp

      x=vv0**0.333333_cp
      pFG0=AFERMIGAS*(Z/v0)**1.66666667_cp
      c0=-1.0_cp*log(3.0_cp*K0/pFG0)         ! assumes V in A^3

      select case(iorder)
         case(2) !AP1
            c2=0._cp
            c3=0._cp

         case(3) !AP2
            c2=1.5_cp*(Kp-3.0_cp)-c0
            c3=0._cp

         case(4) !AP3
            c2=1.5_cp*(Kp-3.0_cp)-c0
            !c3 in steps, using the Holzapfel expression for Kpp0
            c3=-9._cp*kpp*k0
            c3=(20._cp + 12._cp*c0 + c0*c0 + 2.0_cp*c2*(9.0_cp+c0) + 4.0_cp*c2*c2 - c3)/6.0_cp
      end select

      !> terms in pressure expression up to AP3
      a(1,1)=1.0_cp/x**5.0_cp*(1.0_cp-x)
      a(1,2)=exp(c0*(1.0_cp-x))
      a(1,3)=1.0_cp+c2*x*(1.0_cp-x) + c3*x*(1.0_cp-x)**2.0_cp     !Only one with extra term for AP3

      !> first derivatives wrt x up to AP2
      a(2,1)= -5.0_cp/x**6.0_cp +4.0_cp/x**5.0_cp
      a(2,2)= -1.0_cp*c0*a(1,2)
      a(2,3)= c2*(1.0_cp-2.0_cp*x) + c3*(1.0_cp-4.0_cp*x+3.0_cp*x*x) !Only one with extra term for AP3

      !> second derivatives wrt x up to AP2
      a(3,1)= 30.0_cp/x**7.0_cp - 20.0_cp/x**6.0_cp
      a(3,2)=c0*c0*a(1,2)
      a(3,3)=-2.0_cp*c2 + c3*(-4.0_cp*x+6.0_cp*x) !Only one with extra term for AP3

      return
   End Function Get_APL

   !!----
   !!---- FUNCTION GET_Cp
   !!----
   !!---- Calculates the heat capacity at constant pressure by the Gruenesien relation
   !!----
   !!---- If problem, like no gamma0, or  linear eos, default return is Cp=0.
   !!----
   !!---- Date: 01/04/2020
   !!
   Function Get_Cp(P, T, Eos) result(C)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P,T
      type(Eos_Type), intent(in) :: EoS     ! Eos Parameters

      !---- Local Variables ----!
      real(kind=cp) :: C           !Heat capacity in J/mol/K if V is molar
      real(kind=cp) :: v,k,al,gamma2,cv,gamma

      !> Default
      C=0.0_cp

      !> checks
      if (abs(eos%params(18)) < tiny(0.0_cp)) return          ! gamma=0. therefore Cp undefined

      !> get V and K and into SI units
      v=get_volume(p,t,eos)

      k=k_cal(v,t,eos,p=p)
      if (index(U_case(eos%pscale_name),'GPA') > 0)  k=k*1.0E9
      if (index(U_case(eos%pscale_name),'KBAR') > 0) k=k*1.0E8

      select case(eos%itherm)
         case(6,7,8)                           !Calculation from Cv, avoiding alpha, which would be recursive
            v=get_volume(p,t,eos)
            gamma2=get_grun_th(P,T,Eos)**2.0_cp
            !          gamma2=get_grun_V(V,Eos)**2.0_cp
            if (VscaleMGD(eos)) v=v*1.0E-6     !V now in m3/mol
            cv=get_cv(p,t,eos)
            c=(1.0_cp + gamma2*Cv*T/k/v)*Cv

         case default   !from Cp=(1+alpha.gamma.T)Cv
            al=Alpha_Cal(P,T,eos)
            gamma=Get_Grun_V(v,eos)
            c=(1.0_cp + al*gamma*T)*get_cv(p,t,eos)
      end select

      return
   End Function Get_Cp

   !!----
   !!---- FUNCTION GET_Cv
   !!----
   !!---- Calculates the heat capacity at constant volume by the Gruenesien relation for most EoS
   !!----  for those that are oscillator-based, does explicit calculation from oscillator equations
   !!---- If problem, like no gamma0, or  linear eos, default return is Cv=0.
   !!----
   !!---- Date: 01/04/2020
   !!
   Function Get_Cv(P, T, Eos,j) result(Cv)
      !---- Arguments ----!
      real(kind=cp),     intent(in) :: P,T
      type(Eos_Type),    intent(in) :: EoS   ! Eos Parameters
      integer, optional, intent(in) :: j     ! which oscillator: then Pthermal only calculates for this one

      !---- Local Variables ----!
      integer                      :: i,jo
      real(kind=cp)                :: Cv           !Heat capacity in J/mol/K if V is molar
      real(kind=cp),dimension(0:2) :: Cvpart
      real(kind=cp)                :: v,k,thetaD,gammaV,x

      !> init
      jo=-1     !calculate all:
      if (present(j))then
         if (j > -1 .and. j < 3)jo=j
      end if
      cvpart=0.0_cp
      Cv=0.0_cp

      v=get_volume(p,t,eos)

      select case(eos%itherm)
         case(7)                     !MGD
            thetaD=get_DebyeT(V,Eos)
            x=thetaD/t
            if (x < huge(0.))  &
               cvpart(0)=3.0_cp*eos%params(13)*8.314_cp * (4.0_cp*debye(3,x) -3.0_cp*x/(exp(x)-1))
               ! no scaling, units are in R=8.314 J/mol/K



         case(6,8)  !Einstein: includes Holland-Powell
            x=get_DebyeT(V,Eos)/t
            if (x < 20)     &        !corresponds to 0.05ThetaE where Cv < 0.00002 J/mol/K
               cvpart(0)=3.0_cp*eos%params(13)*8.314_cp * x**2._cp * exp(x)/(exp(x)-1)**2._cp
               ! no scaling, units are in R=8.314 J/mol/K

         case default    ! Cv = alpha.Kt/gamma/V
            if (abs(eos%params(18)) < tiny(0.)) return          ! gamma=0. therefore Cv undefined
            k=k_cal(v,t,eos,p=p)
            cv=Alpha_Cal(P,T,eos)*k/Get_Grun_V(v,eos) * v
            ! scaling when getting Cv from other params
            ! if (index(U_case(eos%pscale_name),'GPA') > 0)  factor=1.0E9
            ! if (index(U_case(eos%pscale_name),'KBAR') > 0) factor=1.0E8
            ! if (VscaleMGD(eos)) factor=factor*1.0E-6     !test for cm3/mol or equivalent in eos%vscale_name
            cv=cv/EPthermal_factor(eos)
            return        ! this approach not compatible with mode calculations
      end select

      !> Extra oscillators: only allowed in combination with models 6,7,8
      if (eos%osc_allowed .and. sum(eos%iosc) > 0)then
         cvpart(0)=(1._cp-eos%params(40)-eos%params(45))*cvpart(0)     ! partial contribution main oscillator

         do i=1,2
            select case(eos%iosc(i))
               case(0)
                  cycle

               case(1)  !DEBYE
                  thetaD=get_DebyeT(V,Eos,i)
                  gammaV=get_grun_V(V,Eos,i)
                  x=thetaD/t
                  if (x < huge(0.))  &
                     cvpart(i)=eos%params(35+5*i)*3.0_cp*eos%params(13)*8.314_cp * (4.0_cp*debye(3,x) -3.0_cp*x/(exp(x)-1))

               case(2)     ! Einstein
                  x=get_DebyeT(V,Eos,i)/t
                  if (x < 20)     &        !corresponds to 0.05ThetaE where Cv < 0.00002 J/mol/K
                     cvpart(i)=eos%params(35+5*i)*3.0_cp*eos%params(13)*8.314_cp * x**2._cp * exp(x)/(exp(x)-1)**2._cp

            end select
         end do
      end if

      !> Now return requested part of pth:
      if (jo == -1)then
         cv=sum(cvpart)
      else
         cv=cvpart(jo)
      end if

      return
   End Function Get_Cv

   !!----
   !!---- FUNCTION GET_DEBYET
   !!----
   !!---- Calculate the Debye Temperature at V
   !!----
   !!---- Date: 16/03/2017
   !!
   Function Get_DebyeT(V, Eos,i) result(DebyeT)
      !---- Arguments ----!
      real(kind=cp),    intent(in) :: V       ! Volume or length
      type(Eos_Type),   intent(in) :: EoS     ! Eos Parameters
      integer,optional, intent(in) :: i       ! which oscillator
      real(kind=cp)                :: DebyeT

      !---- Local Variables ----!
      integer       :: io
      real(kind=cp) :: gammaV,V0V

      !> Default
      DebyeT=eos%tref

      if (.not. eos%Osc_allowed)return

      V0V=eos%params(1)/V
      if (eos%linear) V0V=V0V**3.0_cp

      !> local copy of pointer to oscillator
      io=0
      if (present(i))then
         if (i > 0 .and. i <= N_OSC_MODELS)io=i
      end if

      if (io == 0)then                    !main thermal model
                                          !> For linear uses the same parameter values, no factor of 3
         if (eos%params(14) > 0.5_cp)then
            !> q-compromise
            DebyeT=eos%params(11)
         else
            !> normal
            if (abs(eos%params(19)) < tiny(0._cp)) then
               DebyeT=eos%params(11)*V0V**eos%params(18)                  ! when q=0, gamma=gamma0
            else
               gammaV=get_Grun_v(v,eos)               ! Get_grun knows about linear/volume
               DebyeT=eos%params(11)*exp((eos%params(18)-gammaV)/eos%params(19))  ! if q=0 then this gives DebyeT=nan
            end if
         end if

      else                              !extra oscillator: Debye and einstein
         if (eos%iosc(io) == 0)return

         if (eos%params(39+5*io) > 0.5_cp)then
            !> q-compromise
            DebyeT=eos%params(36+5*io)
         else
            !> normal
            if (abs(eos%params(38+5*io)) < tiny(0._cp)) then
               DebyeT=eos%params(36+5*io)*V0V**eos%params(37+5*io)                 ! when q=0, gamma=gamma0
            else
               gammaV=get_Grun_v(v,eos,io)               ! Get_grun knows about linear/volume
               DebyeT=eos%params(36+5*io)*exp((eos%params(37+5*io)-gammaV)/eos%params(38+5*io))  ! if q=0 then this gives DebyeT=nan
            end if
         end if
      end if

      return
   End Function Get_DebyeT

   !!--++
   !!--++ FUNCTION GET_DMDT_AXIS
   !!--++
   !!--++ Returns the value of temperature derivative of the modulus of principal axis (ieos) in unit
   !!--++ cell in cell_eos at P,T
   !!--++
   !!--++ Call this Function directly when the calling routine  knows that the direction is a principal axis
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_DmDt_Axis(P, T, Cell_eos, Ieos) Result(dMdT)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: P,T
      type(eos_cell_type),intent(in)  :: cell_eos
      integer,            intent(in)  :: ieos      !axis indicator, as in axis_type%ieos
      real(kind=cp)                   :: dMdT

      !---- Local Variables ----!

      !> init
      dmdt=0._cp

      select case(cell_eos%loaded(ieos))
         case(1)
            dmdt=dKdT_Cal(p,t,cell_eos%eos(ieos))

         case(2) ! sym equiv. Always uses eos(1) for a-axis
            dmdt=dKdT_Cal(p,t,cell_eos%eos(1))

         case(3)
            dmdt=get_dmdt_third(p,T,cell_eos,ieos)

         case(4)
            dmdt=dKdt_cal(p,t,cell_eos%eos(cell_eos%unique))
      end select

      return
   End Function Get_DmDt_Axis

   !!--++
   !!--++ FUNCTION GET_DMDT_CELL
   !!--++
   !!--++ Returns the value of temperature derivative of modulus of any axis in unit cell in cell_eos at P,T
   !!--++ Call this Function when the calling routine does not know if the direction is a principal axis or not
   !!--++ If a principal direction is requested, only axis%ieos is required
   !!--++ axis%v and axis%atype only used if axis%ieos=-2
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_DmDt_Cell(P, T, Cell_eos, Axis) Result(dMdT)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: dMdT

      !---- Local Variables ----!

      !> init
      dmdt=0.0_cp

      select case(axis%ieos)      !invalid numbers just return
         case(0:6)   !principal direction for which eos exists, or can be calculated
            dmdt=get_dmdt_axis(p,t,cell_eos,axis%ieos)

         case(-2)   !general direction
            dmdt=get_dmdt_general(p,T,cell_eos,axis)

      end select

      return
   End Function Get_DmDt_Cell

   !!--++
   !!--++ FUNCTION GET_DMDT_GENERAL
   !!--++
   !!--++ Returns the value of temperature derivative of modulus of any axis in unit
   !!--++ cell in cell_eos at P,T
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_DmDt_General(P, T, Cell_eos, Axis) Result(dMdT)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: dMdT

      !---- Local Variables ----!
      real(kind=cp)   :: tstep,tcal
      integer         :: i

      !> for spline
      integer,parameter :: nstep=21   !must be odd
      integer           :: imid
      real(kind=cp),dimension(nstep):: x,y,d2y,dy

      tstep=20.
      tcal=t-int(nstep/2)*tstep
      do i=1,nstep
         x(i)=tcal
         y(i)=get_mod_general(P,Tcal,cell_eos,axis)
         tcal=tcal+tstep
      end do

      call Second_Derivative(x, y, nstep, d2y)
      call First_Derivative(x, y, nstep, d2y, dy)

      imid=int(nstep/2) + 1
      dMdT=dy(imid)

      return
   End Function Get_DmDt_General

   !!--++
   !!--++ FUNCTION GET_DMDT_THIRD
   !!--++
   !!--++ Returns the value of temperature derivative of modulus of a principal axis ieos in unit
   !!--++ cell in cell_eos at P,T when it can be calculated from others
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_DmDt_Third(P, T, Cell_eos, Ieos) Result(modp)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: p,T
      type(eos_cell_type),intent(in) :: cell_eos
      integer,            intent(in) :: ieos     ! the modulus of the axis (1,2,3) or V (0) to be calculated
      real(kind=cp)                  :: modp

      !---- Local Variables ----!
      integer        :: i
      real(kind=cp)  :: Kp,M1p,M2p,M3p,Mangp,vf

      !> init
      modp=0.0_cp

      !> safety check: if mono or triclinic, should only be called if angle poly used
      if (U_case(cell_eos%system(1:4)) == 'TRIC' .or. U_case(cell_eos%system(1:3)) == 'MONO')then
         if (cell_eos%eosang%iangle == 0)then
            err_eos=.true.
            err_eos_mess='Get_DmDt_Third called for mono or triclinic, without angle poly set'
         end if
      end if

      select case(U_case(cell_eos%system(1:4)))
         case('TRIC','MONO','ORTH')
            if (U_case(cell_eos%system(1:4)) == 'ORTHO')then
               Mangp=0._cp
            else
               vf=Get_Angle_Volfactor(P,T,cell_eos)
               Mangp=(Get_Angle_Volfactor_Deriv(P,T,cell_eos,'P') * &
                      Get_Angle_Volfactor_Deriv(P,T,cell_eos,'T')/vf - &
                      Get_Angle_Volfactor_Deriv2(P,T,cell_eos,'PT'))/vf
            end if

            select case(ieos)
               case(0)
                  M1p= dKdT_cal(P,T,cell_eos%eos(1))/Get_K(P,T,cell_eos%eos(1))**2.0_cp
                  M2p= dKdT_cal(P,T,cell_eos%eos(2))/Get_K(P,T,cell_eos%eos(2))**2.0_cp
                  M3p= dKdT_cal(P,T,cell_eos%eos(3))/Get_K(P,T,cell_eos%eos(3))**2.0_cp
                  modp=get_mod_third(P,T,cell_eos,0)**2.0_cp*(M1p+M2p+M3p-Mangp)

               case default
                  modp= dKdT_cal(P,T,cell_eos%eos(0))/Get_K(P,T,cell_eos%eos(0))**2.0_cp
                  do i=1,3
                     if (i == ieos)cycle
                     modp = modp - dKdT_cal(P,T,cell_eos%eos(i))/Get_K(P,T,cell_eos%eos(i))**2.0_cp
                  end do
                  modp=modp+Mangp
                  modp = get_mod_third(P,T,cell_eos,ieos)**2.0_cp * modp
            end select

         case('TRIG','HEXA','TETR')
            select case(ieos)
               case(0)     ! calc V from a and c
                  M1p= dKdT_cal(P,T,cell_eos%eos(1))/Get_K(P,T,cell_eos%eos(1))**2.0_cp
                  M3p= dKdT_cal(P,T,cell_eos%eos(3))/Get_K(P,T,cell_eos%eos(3))**2.0_cp
                  modp=get_mod_third(P,T,cell_eos,0)**2.0_cp*(2.0*M1p+M3p)

               case(1)     ! a from V and c
                  Kp=  dKdT_cal(P,T,cell_eos%eos(0))/Get_K(P,T,cell_eos%eos(0))**2.0_cp
                  M3p= dKdT_cal(P,T,cell_eos%eos(3))/Get_K(P,T,cell_eos%eos(3))**2.0_cp
                  modp=get_mod_third(P,T,cell_eos,1)**2.0_cp*(Kp-M3p)/2.0

               case(3)     ! c from a and V
                  Kp= dKdT_cal(P,T,cell_eos%eos(0))/Get_K(P,T,cell_eos%eos(0))**2.0_cp
                  M1p= dKdT_cal(P,T,cell_eos%eos(1))/Get_K(P,T,cell_eos%eos(1))**2.0_cp
                  modp=get_mod_third(P,T,cell_eos,3)**2.0_cp*(Kp-2.0*M1p)
            end select

         case('CUBI','ISOT')
            select case(ieos)
               case(0)     ! calc volume from a
                  modp=dKdT_cal(P,T,cell_eos%eos(1))/3.0_cp

               case(1,2,3)     ! a,b, or c from V
                  modp=dKdT_cal(P,T,cell_eos%eos(0))*3.0_cp
            end select
      end select

      return
   End Function Get_DmDt_Third

   !!----
   !!---- FUNCTION GET_GPT
   !!----
   !!---- Obtain the value of G (or Glinear) at P and T
   !!----
   !!---- Date: 11/07/2016
   !!
   Function Get_GPT(P, T, EoS) Result(gpt)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P   ! Pressure
      real(kind=cp),  intent(in) :: T   ! Temperature
      type(Eos_Type), intent(in) :: Eos ! EoS Variable

      !---- Local Variables ----!
      integer         :: i
      real(kind=cp)   :: gpt, delp

      !> default
      gpt=eos%params(30)          ! default is g(Pref,Tref)

      !> T variation
      select case(eos%ishear)       ! choice of model
         case(1)                    ! model 1 is polynomial in P and T
            gpt=gpt+(t-eos%tref)*eos%params(34)
            delp=p-eos%pref
            do i=1,3
               if (eos%params(i+30) < tiny(0.0)) exit
               gpt=gpt+eos%params(i+30)*delp**i      ! eg linear term is (P-Pref)*dG/dP with params(31)
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
   Function Get_Grun_PT(P, T, Eos, i) Result(G)
      !---- Arguments ----!
      real(kind=cp),     intent(in) :: P    ! Pressure
      real(kind=cp),     intent(in) :: T    ! Temperarture
      type(Eos_Type),    intent(in) :: EoS  ! Eos Parameter
      integer, optional, intent(in) :: i   ! which oscillator

      !---- Local Variables ----!
      real(kind=cp) :: v,g

      v=get_volume(P,T,eos)
      G=Get_Grun_V(v,Eos,i)

      return
   End Function Get_Grun_PT

   !!----
   !!---- FUNCTION  GET_GRUN_Th
   !!----
   !!---- Returns thermal Gruneisen parameter at this volume
   !!----
   !!---- Date: 26/05/2020 Restructured from previous version of 13/05/2020 for special cases
   !!
   Function Get_Grun_Th(P, T, Eos)  Result(G)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P,T
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter

      !---- Local Variables ----!
      integer              :: i
      integer,dimension(1) :: ii
      real(kind=cp)        :: g,v,sumc,sumb,sumt,cvf,cv
      real(kind=cp),dimension(0:2) :: Tc,Cvi

      !> init
      G=0._cp

      v=get_volume(P,T,eos)

      !> No extra oscillators
      if (sum(eos%IOsc) == 0 .or. eos%itran > 0)then
         G=Get_Grun_V(v,Eos)
         return
      end if

      !> get Cv
      Cvi=0._cp
      do i=0,2
         Cvi(i)=Get_Cv(P, T, Eos,i)
      end do
      Cv=Cvi(0)+Cvi(1)+Cvi(2)

      if (Cv > tiny(0._cp))then
         ! normal case at finite T
         sumc=0._cp
         do i=0,2
            sumc=sumc+Get_Grun_V(v,Eos,i)*Get_Cv(P, T, Eos,i)
         end do
         G=sumc/cv

      else
         ! here for Cv = 0 : Either all einstein at low T or T=0
         if (eos%itherm == 7 .or. eos%iosc(1) == 1 .or. eos%iosc(2) == 1)then
            !At least one Debye: only Debye contribute to gamma_th
            sumb=0._cp
            sumt=0._cp

            !main oscillator
            if (eos%itherm == 7)then
               cvf=(1._cp-eos%params(40)-eos%params(45))/(Get_DebyeT(V, Eos,0)**3._cp)
               sumb=cvf
               sumt=cvf*Get_Grun_V(v,Eos,0)
            end if

            !extra osc1
            if (eos%iosc(1) == 1)then
               cvf=eos%params(40)/(Get_DebyeT(V, Eos,1)**3._cp)
               sumb=sumb+cvf
               sumt=sumt+cvf*Get_Grun_V(v,Eos,1)
            end if

            !extra osc2
            if (eos%iosc(2) == 1)then
               cvf=eos%params(45)/(Get_DebyeT(V, Eos,2)**3._cp)
               sumb=sumb+cvf
               sumt=sumt+cvf*Get_Grun_V(v,Eos,2)
            end if
            G=sumt/sumb
         else
            !All Einstein at very low or zero T
            tc=huge(0._cp)
            tc(0)=Get_DebyeT(V, Eos,0)
            if (eos%iosc(1) == 2)tc(1)=Get_DebyeT(V, Eos,1)
            if (eos%iosc(2) == 2)tc(2)=Get_DebyeT(V, Eos,2)
            ii=minloc(tc)
            i=ii(1)-1   !minloc returns absolute location, not its label
            g=Get_Grun_V(v,Eos,i)
        end if
      end if

      return
   End Function Get_Grun_Th

   !!--++
   !!--++ FUNCTION  GET_GRUN_Th_old
   !!--++
   !!--++ Returns thermal Gruneisen parameter at this volume
   !!--++
   !!--++ Date: 13/05/2020
   !!
   Function Get_Grun_Th_old(P, T, Eospar)  Result(G)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P,T
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      integer       :: i
      integer,dimension(1) :: ii
      real(kind=cp) :: g,v,sumc,sumb,sumt,cvf
      real(kind=cp),dimension(0:2) :: Tc

      !init
      G=0._cp

      v=get_volume(P,T,eospar)

      if (sum(eospar%IOsc) == 0 .or. eospar%itran > 0)then
         G=Get_Grun_V(v,Eospar)
      else
         !here for multi-oscillator models
         if (T > tiny(0.))then
            sumc=0._cp
            do i=0,2
               sumc=sumc+Get_Grun_V(v,Eospar,i)*Get_Cv(P, T, Eospar,i)
            end do
            cvf=Get_Cv(P, T, Eospar)
            if (cvf > tiny(0.))then
               G=sumc/cvf
            else
               ! must have all Einstein oscillators, and if Cv = 0 then dominated by lowest theta mode
               tc=huge(0._cp)
               tc(0)=Get_DebyeT(V, Eospar,0)
               if(eospar%iosc(1) == 2)tc(1)=Get_DebyeT(V, Eospar,1)
               if(eospar%iosc(2) == 2)tc(2)=Get_DebyeT(V, Eospar,2)
               ii=minloc(tc)
               i=ii(1)-1   !minloc returns absolute location, not its label
               g=Get_Grun_V(v,Eospar,i)
            end if

         else  !at lowT limit where all Cv are zero: then Einstein do not contribute to gamma_th
            sumb=0._cp
            sumt=0._cp
            !main oscillator
            if (eospar%itherm == 7)then
               cvf=(1._cp-eospar%params(40)-eospar%params(45))/(Get_DebyeT(V, Eospar,0)**3._cp)
               sumb=cvf
               sumt=cvf*Get_Grun_V(v,Eospar,0)
            end if

            !extra osc1
            if (eospar%iosc(1) == 1)then
               cvf=eospar%params(40)/(Get_DebyeT(V, Eospar,1)**3._cp)
               sumb=sumb+cvf
               sumt=sumt+cvf*Get_Grun_V(v,Eospar,1)
            end if

            !extra osc2
            if (eospar%iosc(2) == 1)then
               cvf=eospar%params(45)/(Get_DebyeT(V, Eospar,2)**3._cp)
               sumb=sumb+cvf
               sumt=sumt+cvf*Get_Grun_V(v,Eospar,2)
            end if
            G=sumt/sumb
         end if
      end if

      return
   End Function Get_Grun_Th_old

   !!----
   !!---- FUNCTION  GET_GRUN_V
   !!----
   !!---- Returns Gruneisen parameter at this volume as gamma0*(V/V0)
   !!---- If linear it calculates gamma0*(a/a0)^3
   !!----
   !!---- Date: 18/07/2016
   !!
   Function Get_Grun_V(V, Eos, i) Result(Grun)
      !---- Arguments ----!
      real(kind=cp),     intent(in) :: V    ! Volume or length if linear
      type(Eos_Type),    intent(in) :: EoS  ! Eos Parameter
      integer, optional, intent(in) :: i    ! which oscillator

      !---- Local Variables ----!
      integer       :: io
      real(kind=cp) :: Grun                 !The resulting gruneisen gamma
      real(kind=cp) :: VV0,q

      !>local copy of pointer to oscillator
      io=0
      if (present(i))then
         if (i > 0 .and. i <= N_OSC_MODELS)io=i
      end if

      !Init
      grun=0._cp

      !> Must be careful with transitions because eospar%params(1) is the high phase V0
      !> V0=get_volume(eospar%pref,eospar%tref,eospar) (Nov 2016)
      !> no I don't think so. Grun is a property of the high phase

      VV0=V/eos%params(1)
      if (eos%linear) VV0=VV0**3.0_cp

      if (io == 0)then                    !main thermal model
         if (eos%osc_allowed .and. eos%params(14) > 0.5_cp )then
            !q-compromise: gamma/V is constant: MGD and Einstein are allowed q-comp, Holland-Powell is always q-comp
            Grun=eos%params(18)*VV0

         else
            !Normal
            if (abs(eos%params(19)) > 0.00001_cp) then
               VV0=VV0**eos%params(19)
            else
               VV0=1.0_cp
            end if
            Grun=eos%params(18)*VV0
         end if
      else if(eos%iosc(io) > 0) then          !additonal Debye or Einstein oscillator
         if (eos%params(39+5*io) > 0.5_cp)then
            Grun=eos%params(37+5*io)*VV0

         else
            !Normal
            q =eos%params(38+5*io)      ! 43 or 48
            if (abs(q) > tiny(0._cp))then
               VV0=VV0**q
            else
               VV0=1.0_cp
            end if
            Grun=eos%params(37+5*io)*VV0
         end if
      end if

      return
   End Function Get_Grun_V

   !!----
   !!---- FUNCTION GET_K
   !!----
   !!---- Returns the value of K (or M if linear) at P and T
   !!---- Works by using get_volume(P,T) and then using the V in k_cal
   !!----
   !!---- Date: 18/07/2016
   !!
   Function Get_K(P, T, Eos) Result(K)
      !---- Arguments ----!
      real(kind=cp), intent(in)    :: p    ! Pressure
      real(kind=cp), intent(in)    :: t    ! Tenperature
      type(Eos_Type),intent(in)    :: Eos  ! Eos Parameters
      real(kind=cp)                :: k

      !---- Local Variables ----!
      real(kind=cp) :: v

      v=get_volume(p,t,eos)
      k=k_cal(v,t,eos,p=p)

      return
   End Function Get_K

   !!--++
   !!--++ FUNCTION GET_K0_T
   !!--++
   !!--++ Returns k0 needed for Eos calculations at T this means for Pthermal,
   !!--++ k at Tref is returned.
   !!--++ In the linear case then  M(T, P=0) from params is returned
   !!--++
   !!--++ If k is calculated as being negative, an error message is set
   !!--++ and k0 is returned as the value at Tref for safety.
   !!--++
   !!--++ Date:17/07/2015
   !!
   Function Get_K0_T(T, Eos) Result(k0)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: k0

      !---- Local Variables ----!
      real(kind=cp) :: vr

      !> Init
      k0=eos%params(2) !default (and correct if no thermal model)

      if (.not. eos%Pthermaleos)then
         select case(eos%icross)
            case(1)
               k0=eos%params(2)+eos%params(8)*(t-eos%tref)  !Old linear variation of K with T

            case(2)
               vr=eos%params(1)/Get_V0_T(T,Eos)          ! Get_V0_T returns a0 for linear
               if (eos%linear) vr=vr**3.0_cp
               if (vr > 0.001 .and. vr < huge(0.0_cp)) k0=eos%params(2)*vr**eos%params(8)   !Anderson Gruneisen approach using params(8) as delta
         end select
      end if

      return
   End Function Get_K0_T

   !!----
   !!---- FUNCTION GET_Kp
   !!----
   !!---- Returns the value of K (or M if linear) at P and T
   !!---- Works by using get_volume(P,T) and then using the V in k_cal
   !!----
   !!---- Date: 1/03/2018
   !!
   Function Get_Kp(P, T, Eos) Result(Kp)
      !---- Arguments ----!
      real(kind=cp),intent(in)        :: p    ! Pressure
      real(kind=cp),intent(in)        :: t    ! Tenperature
      type(Eos_Type),intent(in)       :: Eos  ! Eos Parameters
      real(kind=cp)                   :: kp

      !---- Local Variables ----!
      real(kind=cp) :: v

      v=get_volume(p,t,eos)
      kp=kp_cal(v,t,eos,p=p)

      return
   End Function Get_Kp

   !!--++
   !!--++ FUNCTION GET_KP0_T
   !!--++
   !!--++ Returns kp0 needed for Eos calculations at T this means for Pthermal,
   !!--++ kp at Tref is returned.
   !!--++
   !!--++ In the linear case then  Mp(T, P=0) from params is returned
   !!--++
   !!--++ Date:17/11/2016
   !!
   Function Get_Kp0_T(T, Eos) Result(kp0)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: kp0

      !---- Local Variables ----!
      real(kind=cp) :: vr

      !> Init
      kp0=eos%params(3) !default (and correct if no thermal model)

      if (.not. eos%Pthermaleos)then
         select case(eos%icross)
            case(1)  !Old linear variation of K with T, no Kp variation

            case(2)
               vr=Get_V0_T(T,Eos)/eos%params(1)            !normally Vr > 0, but if negative thermal expansion, not
               if (eos%linear) vr=vr**3.0_cp
               if (vr > 0.001 .and. vr < huge(0._cp) ) kp0=eos%params(3)*vr**eos%params(9)   ! params(9) is delta-prime
         end select
      end if

      return
   End Function Get_Kp0_T

   !!--++
   !!--++ FUNCTION GET_KPP0_T
   !!--++
   !!--++ Returns kpp0 needed for Eos calculations at T this means for pthermal,
   !!--++ kpp at Tref is returned.
   !!--++
   !!--++ In the linear case then  Mp(T, P=0) from params is returned
   !!--++
   !!--++ Date:17/11/2016
   !!
   Function Get_Kpp0_T(T, Eos) Result(kpp0)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: kpp0

      !---- Local Variables ----!
      type(Eos_Type) :: eost     ! workspace

      !> Init
      kpp0=eos%params(4) !default (and correct if no thermal model, or if Pthermal model)

      select case(eos%itherm)
         case(1:5)  ! Normal isothermal eos at T
            eost=eos           ! Modify eost to hold K0 and Kp at the T of interest

            eost%params(2)=Get_K0_T(T,Eos)
            eost%params(3)=Get_Kp0_T(T,Eos)
            call set_eos_implied_values(Eost)

            kpp0=eost%params(4)
      end select

      return
   End Function Get_Kpp0_T



   !!----
   !!---- FUNCTION GET_MOD_AXIS
   !!----
   !!---- Returns the value of modulus of principal axis (ieos) in unit
   !!---- cell in cell_eos at P,T
   !!----
   !!---- Call this Function directly when the calling routine  knows that
   !!---- the direction is a principal axis
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Mod_Axis(P, T, Cell_eos, Ieos) result(Modu)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      integer,            intent(in)  :: ieos      !axis indicator, as in axis_type%ieos
      real(kind=cp)                   :: modu      !returned modulus

      !---- Local Variables ----!

      !>init
      modu=0.0_cp

      select case(cell_eos%loaded(ieos))
         case(1)
            modu=get_k(p,t,cell_eos%eos(ieos))

         case(2) ! sym equiv. Always uses eos(1) for a-axis
            modu=get_k(p,t,cell_eos%eos(1))

         case(3)
            modu=get_mod_third(p,T,cell_eos,ieos)

         case(4)
            modu=get_k(p,t,cell_eos%eos(cell_eos%unique))
      end select

      return
   End Function Get_Mod_Axis

   !!----
   !!---- FUNCTION GET_MOD_CELL
   !!----
   !!---- Returns the value of modulus of any axis in unit cell in cell_eos at P,T
   !!---- Call this Function when the calling routine does not know if the
   !!---- direction is a principal axis or not
   !!----
   !!---- If a principal direction is requested, only axis%ieos is required
   !!---- axis%v and axis%atype only used if axis%ieos=-2
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Mod_Cell(P, T, Cell_eos, Axis) Result(Modu)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: modu !returned modulus

      !---- Local Variables ----!

      !> init
      modu=0.0_cp

      select case(axis%ieos)      !invalid numbers just return
         case(0:6)   !principal direction for which eos exists, or can be calculated
            modu=get_mod_axis(p,t,cell_eos,axis%ieos)

         case(-2)   !general direction
            modu=get_mod_general(p,T,cell_eos,axis)
      end select

      return
   End Function Get_Mod_Cell

   !!--++
   !!--++ FUNCTION GET_MOD_GENERAL
   !!--++
   !!--++ Returns the value of modulus of any axis in unit cell in cell_eos at P,T
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Mod_General(P, T, Cell_eos, Axis) Result(Modu)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: modu

      !---- Local Variables ----!
      integer         :: i
      real(kind=cp)   :: k,kmax,pstep,pcal,vm,vp

      !> for spline
      integer,parameter              :: NSTEP=11   !must be odd
      integer                        :: imid
      real(kind=cp),dimension(NSTEP) :: x,y,d2y,dy

      !> find largest linear modulus of axes
      kmax=tiny(0._cp)
      do i=1,3
         k=get_mod_axis(P,T,cell_eos,i)
         if (k > kmax) kmax=k
      end do

      !> now do pstep as kmax/1000.  Seems good compromise from tests - also depends on nstep
      pstep=kmax/500.

      !> approximate modulus
      vm=get_Volume_cell(P-0.5_cp,T,cell_eos,axis)
      vp=get_Volume_cell(P+0.5_cp,T,cell_eos,axis)
      K=(vm+vp)/2.0_cp/(vm-vp)
      pstep=K/500._cp  !should give delv of ca 1%

      pcal=p-int(NSTEP/2)*pstep
      do i=1,NSTEP
         x(i)=pcal
         y(i)=get_Volume_cell(Pcal,T,cell_eos,axis)
         pcal=pcal+pstep
      end do

      call Second_Derivative(x, y, NSTEP, d2y)
      call First_Derivative(x, y, NSTEP, d2y, dy)

      imid=int(NSTEP/2) + 1
      Modu=-1.0*y(imid)/dy(imid)

      return
   End Function  Get_Mod_General

   !!--++
   !!--++ FUNCTION GET_MOD_THIRD
   !!--++
   !!--++ Returns the value of modulus of a principal axis ieos in unit
   !!--++ cell in cell_eos at P,T when it can be calculated from others
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Mod_Third(P, T, Cell_eos, Ieos) Result(Modu)
      !---- Arguments ----!
      real(kind=cp),       intent(in) :: p,T
      type(eos_cell_type), intent(in) :: cell_eos
      integer,             intent(in) :: ieos     ! the modulus of the axis (1,2,3) or V (0) to be calculated
      real(kind=cp)                   :: modu

      !---- Local Variables ----!
      integer       :: i
      real(kind=cp) :: beta_ang

      !>init
      modu=0.0_cp

      !>safety check: if mono or triclinic, should only be called if angle poly used
      if (U_case(cell_eos%system(1:4)) == 'TRIC' .or. U_case(cell_eos%system(1:3)) == 'MONO')then
         if (cell_eos%eosang%iangle == 0)then
            err_eos=.true.
            err_eos_mess='Get_Mod_Third called for mono or triclinic, without angle poly set'
         end if
      end if

      select case(U_case(cell_eos%system(1:4)))
         case('TRIC','MONO','ORTH')
            beta_ang=Get_Angle_Volfactor_Deriv(P,T,cell_eos,'P')/Get_Angle_Volfactor(P,T,cell_eos)  !1/A dA/dP
            select case(ieos)
               case(0)
                  modu=  1.0_cp/(1.0_cp/Get_K(P,T,cell_eos%eos(1))+1.0_cp/Get_K(P,T,cell_eos%eos(2))+ &
                         1.0_cp/Get_K(P,T,cell_eos%eos(3))- beta_ang)

               case default
                  modu=  1.0_cp/Get_K(P,T,cell_eos%eos(0))
                  do i=1,3
                     if (i == ieos)cycle
                     modu= modu - 1.0_cp/Get_K(P,T,cell_eos%eos(i))
                  end do
                  modu=modu+beta_ang
                  modu=1.0_cp/modu
            end select

         case('TRIG','HEXA','TETR')
            select case(ieos)
               case(0)     ! calc V from a and c
                  modu= 1.0_cp/(2.0_cp/Get_K(P,T,cell_eos%eos(1))+1.0_cp/Get_K(P,T,cell_eos%eos(3)))

               case(1)     ! a from V and c
                  modu= 2.0_cp/(1.0_cp/Get_K(P,T,cell_eos%eos(0))-1.0_cp/Get_K(P,T,cell_eos%eos(3)))

               case(3)     ! c from a and V
                  modu= 1.0_cp/(1.0_cp/Get_K(P,T,cell_eos%eos(0))-2.0_cp/Get_K(P,T,cell_eos%eos(1)))
            end select

         case('CUBI','ISOT')
            select case(ieos)
               case(0)     ! calc volume from a
                  modu= Get_K(P,T,cell_eos%eos(1))/3.0_cp

               case(1,2,3)     ! a,b, or c from V
                  modu= Get_K(P,T,cell_eos%eos(0))*3.0_cp
            end select
      end select

      return
   End Function Get_Mod_Third

   !!----
   !!---- FUNCTION GET_MODP_AXIS
   !!----
   !!---- Returns the value of pressure derivative of the modulus of principal axis (ieos)
   !!---- in unit cell in cell_eos at P,T
   !!---- Call this Function directly when the calling routine  knows that the direction is a principal axis
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Modp_Axis(P, T, Cell_eos, Ieos) Result(Modp)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      integer,            intent(in)  :: ieos      !axis indicator, as in axis_type%ieos
      real(kind=cp)                   :: modp !returned modulus derivative

      !---- Local Variables ----!

      !> init
      modp=0.0_cp

      select case(cell_eos%loaded(ieos))
         case(1)
            modp=get_kp(p,t,cell_eos%eos(ieos))

         case(2) ! sym equiv. Always uses eos(1) for a-axis
            modp=get_kp(p,t,cell_eos%eos(1))

         case(3)
            modp=get_modp_third(p,T,cell_eos,ieos)

         case(4)
            modp=get_kp(p,t,cell_eos%eos(cell_eos%unique))
      end select

      return
   End Function Get_Modp_Axis

   !!----
   !!---- FUNCTION GET_MODP_CELL
   !!----
   !!---- Returns the value of pressure derivative of modulus of any axis in unit
   !!---- cell in cell_eos at P,T
   !!---- Call this Function when the calling routine does not know if the direction
   !!---- is a principal axis or not. If a principal direction is requested, only
   !!---- axis%ieos is required axis%v and axis%atype only used if axis%ieos=-2
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Modp_Cell(P, T, Cell_eos, Axis) Result(Modp)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: modp !returned modulus  derivative

      !---- Local Variables ----!

      !>init
      modp=0.0_cp

      select case(axis%ieos)      !invalid numbers just return
         case(0:6)   !principal direction for which eos exists, or can be calculated
            modp=get_modp_axis(p,t,cell_eos,axis%ieos)

         case(-2)   !general direction
            modp=get_modp_general(p,T,cell_eos,axis)
      end select

      return
   End Function Get_Modp_Cell

   !!--++
   !!--++ FUNCTION GET_MODP_GENERAL
   !!--++
   !!--++ Returns the value of pressure derivative of modulus of any axis in unit cell in cell_eos at P,T
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Modp_General(P, T, Cell_eos, Axis) Result(Mp)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: mp

      !---- Local Variables ----!
      integer         :: i
      real(kind=cp)   :: k,kmax,pstep,pcal

      !> for spline
      integer,parameter                 :: NSTEP=11   !must be odd
      integer                           :: imid
      real(kind=cp),dimension(NSTEP)    :: x,y,d2y,dy

      !> find largest linear modulus of axes
      kmax=tiny(0.0_cp)
      do i=1,3
         k=get_mod_axis(P,T,cell_eos,i)
         if (k > kmax) kmax=k
      end do

      !> Seems good compromise from tests - also depends on nstep
      pstep=kmax/500.
      pcal=p-int(NSTEP/2)*pstep

      do i=1,NSTEP
         x(i)=pcal
         y(i)=get_mod_general(Pcal,T,cell_eos,axis)
         pcal=pcal+pstep
      end do

      call Second_Derivative(x, y, NSTEP, d2y)
      call First_Derivative(x, y, NSTEP, d2y, dy)
      imid=int(NSTEP/2) + 1
      Mp=dy(imid)

      return
   End Function Get_Modp_General

   !!--++
   !!--++ FUNCTION GET_MODP_THIRD
   !!--++
   !!--++ Returns the value of pressure derivative of modulus of a principal axis
   !!--++ ieos in unit cell in cell_eos at P,T when it can be calculated from others
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Modp_Third(P, T, Cell_Eos, Ieos) Result(Modp)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: p,T
      type(eos_cell_type),intent(in) :: cell_eos
      integer,            intent(in) :: ieos     ! the modulus of the axis (1,2,3) or V (0) to be calculated
      real(kind=cp)                  :: modp

      !---- Local Variables ----!
      real(kind=cp) :: Kp,M1p,M2p,M3p,Mangp,Vf
      integer       :: i

      !> init
      modp=0.0_cp

      !> safety check: if mono or triclinic, should only be called if angle poly used
      if (U_case(cell_eos%system(1:4)) == 'TRIC' .or. U_case(cell_eos%system(1:3)) == 'MONO')then
         if (cell_eos%eosang%iangle == 0)then
            err_eos=.true.
            err_eos_mess='Get_Modp_Third called for mono or triclinic, without angle poly set'
         end if
      end if

      select case(U_case(cell_eos%system(1:4)))
         case('TRIC','MONO','ORTH')
            if (U_case(cell_eos%system(1:4)) == 'ORTHO')then
               Mangp=0._cp
            else
               vf=Get_Angle_Volfactor(P,T,cell_eos)
               Mangp=(Get_Angle_Volfactor_Deriv(P,T,cell_eos,'P')**2.0_cp/vf - &
                      Get_Angle_Volfactor_Deriv2(P,T,cell_eos,'P'))/vf
            end if

            select case(ieos)
               case(0)
                  M1p= get_Kp(P,T,cell_eos%eos(1))/Get_K(P,T,cell_eos%eos(1))**2.0_cp
                  M2p= get_Kp(P,T,cell_eos%eos(2))/Get_K(P,T,cell_eos%eos(2))**2.0_cp
                  M3p= get_Kp(P,T,cell_eos%eos(3))/Get_K(P,T,cell_eos%eos(3))**2.0_cp
                  modp=get_mod_third(P,T,cell_eos,0)**2.0_cp*(M1p+M2p+M3p-Mangp)

               case default
                  modp= get_Kp(P,T,cell_eos%eos(0))/Get_K(P,T,cell_eos%eos(0))**2.0_cp
                  do i=1,3
                     if (i == ieos)cycle
                     modp = modp - get_Kp(P,T,cell_eos%eos(i))/Get_K(P,T,cell_eos%eos(i))**2.0_cp
                  end do
                  modp=modp+Mangp
                  modp = get_mod_third(P,T,cell_eos,ieos)**2.0_cp * modp
            end select

         case('TRIG','HEXA','TETR')
            select case(ieos)
               case(0)     ! calc V from a and c
                  M1p= get_Kp(P,T,cell_eos%eos(1))/Get_K(P,T,cell_eos%eos(1))**2.0_cp
                  M3p= get_Kp(P,T,cell_eos%eos(3))/Get_K(P,T,cell_eos%eos(3))**2.0_cp
                  modp=get_mod_third(P,T,cell_eos,0)**2.0_cp*(2.0*M1p+M3p)

               case(1)     ! a from V and c
                  Kp= get_Kp(P,T,cell_eos%eos(0))/Get_K(P,T,cell_eos%eos(0))**2.0_cp
                  M3p= get_Kp(P,T,cell_eos%eos(3))/Get_K(P,T,cell_eos%eos(3))**2.0_cp
                  modp=get_mod_third(P,T,cell_eos,1)**2.0_cp*(Kp-M3p)/2.0

               case(3)     ! c from a and V
                  Kp= get_Kp(P,T,cell_eos%eos(0))/Get_K(P,T,cell_eos%eos(0))**2.0_cp
                  M1p= get_Kp(P,T,cell_eos%eos(1))/Get_K(P,T,cell_eos%eos(1))**2.0_cp
                  modp=get_mod_third(P,T,cell_eos,3)**2.0_cp*(Kp-2.0*M1p)
            end select

         case('CUBI','ISOT')
            select case(ieos)
               case(0)     ! calc volume from a
                  modp=get_Kp(P,T,cell_eos%eos(1))/3.0_cp

               case(1:3)     ! a,b, or c from V
                  modp=get_Kp(P,T,cell_eos%eos(0))*3.0_cp
            end select
      end select

      return
   End Function Get_Modp_Third



   !!----
   !!---- FUNCTION GET_PARAMS_CELL
   !!----    returns full cell parameters in crystal_cell_type for input P and T
   !!----
   !!---- Date: 04/02/2021
   !!
   Function Get_Params_Cell(P, T, Cell_Eos, Cartype) Result(Cell)
      !---- Arguments ----!
      real(kind=cp),              intent(in) :: p,t
      type(eos_cell_type),        intent(in) :: cell_eos
      character(len=2), optional, intent(in) :: cartype    ! orientation
      type(crystal_cell_type)                :: cell  !cell params, metric tensor at this P,T

      !---- Local Variables ----!
      integer                     :: i,j,k
      real(kind=cp)               :: v,arg
      real(kind=cp),dimension(3)  :: abc,ang
      character(len=2)            :: ctype


      !> init
      abc=10._cp
      ctype='  '
      if (present(cartype))ctype=U_case(cartype)

      !> check if all needed eos are loaded
      call Eos_Cell_Loaded_Check(cell_eos)
      if (warn_eos)return

      !> Get cell edges
      do i=1,3
         call init_err_eos()

         abc(i)=get_Volume_axis(P,T,cell_eos,i)
         if (Err_EoS)return
      end do

      !> get the angles
      ang=90._cp
      select case(U_case(cell_eos%system(1:4)))
         case('TRIG','HEXA')
            ang(3)=120._cp

         case('MONO')
            if (cell_eos%eosang%iangle == 0)then
               v=get_Volume_axis(P,T,cell_eos,0)
               i=cell_eos%unique
               arg=V/product(abc)
               if (arg > 0.999999_cp)then
                  ang(i)=90._cp
               else
                  ang(i)=asind(arg)
                  if (cell_eos%obtuse(i))ang(i)=180.0_cp-ang(i)
               end if

            else
               ang(cell_eos%unique)=get_angle_poly(p,t,cell_eos%eosang,cell_eos%unique)
            end if

         case('TRIC')
            if (cell_eos%eosang%iangle == 0)then
               !all eos for V,abc, and d's present
               v=get_volume(p,t,cell_eos%eos(0))
               do i=1,3
                  j=mod(i,3)+1
                  k=mod(j,3)+1
                  arg=V/abc(j)/abc(k)/get_volume(p,t,cell_eos%eos(i+3))
                  if (arg > 0.999999_cp)then
                     ang(i)=90._cp
                  else
                     ang(i)=asind(arg)
                     if (cell_eos%obtuse(i))ang(i)=180.0_cp-ang(i)
                  end if
               end do
            else
               do i=1,3
                  ang(i)=get_angle_poly(p,t,cell_eos%eosang,i)
               end do
            end if

      end select

      !> Set the metric tensor
      call Set_Crystal_Cell (abc, Ang, Cell,cartype=ctype)

      return
   End Function Get_Params_Cell

   !!----
   !!---- FUNCTION GET_PRESS_AXIS
   !!----
   !!---- Returns the value of pressure at input length or volume of principal axis
   !!---- ioes in unit cell in cell_eos
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Press_Axis(A, T, Cell_eos, Ieos)  Result(P)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: a,T       ! a is volume or linear value, as appropriate
      type(eos_cell_type),intent(in) :: cell_eos  !the eos for the cell axes
      integer,            intent(in) :: ieos      ! the axis (1,2,3) or V (0) to be calculated
      real(kind=cp)                  :: p

      !---- Local Variables ----!

      !>init
      p=0.0_cp

      select case(cell_eos%loaded(ieos))
         case(1)
            p=get_pressure(a,t,cell_eos%eos(ieos))

         case(2) ! sym equiv. Always uses eos(1) for a-axis
            p=get_pressure(a,t,cell_eos%eos(1))

         case(3)
            p=get_press_third(a,T,cell_eos,ieos)

         case(4)
            p=get_pressure(p,t,cell_eos%eos(cell_eos%unique))
      end select

      return
   End Function Get_Press_Axis

   !!----
   !!---- FUNCTION GET_PRESS_CELL
   !!----
   !!---- Returns the value of pressure at input length or volume of axis in unit cell in cell_eos
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Press_Cell(A, T, Cell_eos, Axis) Result(P)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: a,t      !a is in length not length cubed
      type(eos_cell_type),intent(in) :: cell_eos !the eos for the cell axes
      type(axis_type),    intent(in) :: axis  ! The direction. if a cell axis then only axis%ieos is required
      real(kind=cp)                  :: p

      !---- Local Variables ----!

      !>init
      p=0.0_cp

      select case(axis%ieos)      !invalid numbers just return
         case(0:6)   !principal direction for which eos exists, or can be calculated
            P=get_press_axis(a,t,cell_eos,axis%ieos)

         case(-2)   !general direction
            P=get_press_general(a,T,cell_eos,axis)

      end select

      return
   End Function Get_Press_Cell

   !!--++
   !!--++ FUNCTION GET_PRESS_GENERAL
   !!--++
   !!--++ Returns the value of pressure at input length of a general direction
   !!--++ (axis) in unit cell in cell_eos
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Press_General(A, T, Cell_Eos, Axis) Result(P)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: a,t      !a is in length not length cubed
      type(eos_cell_type),intent(in) :: cell_eos
      type(axis_type),    intent(in) :: axis
      real(kind=cp)                  :: p

      !---- Local Variables ----!
      real(kind=cp) :: del,delprev,step,tol,acalc
      integer       :: ic

      !> init
      p=0.0_cp

      del=0.1
      delprev=del
      step=-1.0_cp*del
      tol=0.0005_cp*a/get_mod_general(P,T,cell_eos,axis) !tolerance in V scaled by M to give 0.0005 error in P

      ic=0
      do while (ic <= 1000)
         ic=ic+1

         !> calc the ratio at the current p
         acalc=get_Volume_general(P,T,cell_eos,axis)
         del=acalc-a

         if (ic > 1)then                      ! need to get two calcs before adjusting step size and dirn
            if(abs(del) .lt. tol) exit

            if (del*delprev < 0._cp)then     ! over-stepped: reverse with half the step size
               step=-0.5_cp*step

            else                            ! same signs
               if (abs(del) > abs(delprev))then ! going the wrong direction
                  step=-1.0_cp*step
               end if
            end if
            if (abs(step) < 0.000001) exit  ! step in p got too small
         end if
         delprev=del
         p=p+step
      end do

      return
   End Function Get_Press_General

   !!--++
   !!--++ FUNCTION GET_PRESS_THIRD
   !!--++
   !!--++ Returns the value of pressure at input length of a principal axis (ieos) in unit cell in cell_eos
   !!--++   when it must be calculated from other eos
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Press_Third(A, T, Cell_Eos, Ieos, Pguess) Result(P)
      !---- Arguments ----!
      real(kind=cp),         intent(in) :: a,T      ! v is volume or linear value, as appropriate
      type(eos_cell_type),   intent(in) :: cell_eos
      integer,               intent(in) :: ieos     ! the axis (1,2,3) or V (0) to be calculated
      real(kind=cp),optional,intent(in) :: pguess   !initial guess to P
      real(kind=cp)                     :: p

      !---- Local Variables ----!
      real(kind=cp) :: del,delprev,step,tol,vcalc
      integer       :: ic

      !> init
      p=0.0_cp
      if (present(pguess))p=pguess
      call init_err_eos()

      del=0.1
      delprev=del
      step=-1.0_cp*del

      tol=0.0005_cp*a/get_mod_third(P,T,cell_eos,ieos)   !tolerance in V scaled by M to give 0.0005 error in P

      ic=0
      do while (ic <= 1000)
         ic=ic+1

         ! calc the ratio at the current p
         vcalc=get_volume_third(P,T,cell_eos,ieos)
         del=vcalc-a

         if (ic > 1)then                      ! need to get two calcs before adjusting step size and dirn
            if (abs(del) < tol)return

            if (del*delprev < 0._cp)then     ! over-stepped: reverse with half the step size
               step=-0.5_cp*step

            else                            ! same signs
               if (abs(del) > abs(delprev))then ! going the wrong direction
                  step=-1.0_cp*step
               end if
            end if
            if (abs(step) < 0.000001)return  ! step in p got too small
         end if
         delprev=del
         p=p+step
      end do

      return
   End Function Get_Press_Third

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
   Function Get_Pressure(V, T, Eos) Result(P)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V    ! Volume
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: p

      !---- Local Variables ----!
      logical                           :: first
      integer                           :: i
      real(kind=cp)                     :: K0,Kp,kpp,vv0,x,u,vol,plast,vs,ptr,vtr,difp,volp
      real(kind=cp)                     :: a,b,c,f
      real(kind=cp),dimension(N_EOSPAR) :: ev
      real(kind=cp),dimension(3)        :: abc      ! for Tait parameters
      real(kind=cp),dimension(3,3)      :: apl      ! for APL parameters

      !> Init
      p=0.0_cp
      vol=v             ! needed to allow for transition strain: vol is the V of the bare phase
      vs=0.0_cp
      plast=0.0_cp
      first=.true.

      !> If PTV table, go directly
      if (EoS%imodel == -1) then
         p=get_props_ptvtable(0.0,t,v,EoS,'P')     ! get_props_ptvtable works in length if linear
         return
      end if



      !> copy Eos parameters to local
      ev= EoS_to_Vec(EoS)         ! Volume or linear case is covered

      !> These parameters depend only on T, not P or V
      k0=Get_K0_T(T,EoS)              ! Handles thermal pressure case, returns K0 or M0
      if (EoS%linear) k0=k0/3.0_cp

      kp=Get_Kp0_T(T,EoS)
      if (EoS%linear) kp=kp/3.0_cp
      kpp=Get_Kpp0_T(T,EoS)
      if (EoS%linear) kpp=kpp/3.0_cp

      !> Start increment loop to get transition factor
      do
         !> Thermal case
         select case (EoS%itherm)
            case (0,6,7,8) ! No thermal case, or  thermal pressure which uses params at Tref
               vv0=vol/EoS%params(1)      !vv0 is now V/V0 or a/a0

            case (1:5)
               vv0=vol/Get_V0_T(T,EoS)    ! In the case of Phase transition, Get_V0_T always returns V0 of high phase at this T: This is correct!
         end select

         f=strain(vv0,EoS)                ! use strain to get finite strain from v/v0 or a/a0
         if (err_eos)then
            err_eos_mess=trim(err_eos_mess)//' called from Get_Pressure'
            exit
         end if

         !> Using volume dimensions
         vv0=1.0_cp/vv0                    ! vv0 now V0/V for easy of use in rest of subprogram
         if (EoS%linear) vv0=vv0**(3.0_cp)

         select case (EoS%imodel)
            case (1) ! Murnaghan
               P=K0/Kp*(vv0**Kp - 1.0_cp)

            case (2) ! Birch-Murnaghan
               a=f*(1.0_cp+2.0_cp*f)**(2.5_cp)     ! changed expressions to use only f 04/04/2013 RJA
               b=0.0_cp
               c=0.0_cp
               if (EoS%iorder > 2)  b=1.5_cp*(Kp-4.0_cp)
               if (EoS%iorder == 4) c = 1.5_cp*(K0*Kpp + (Kp-4.0_cp)*(Kp-3.0_cp)+35.0_cp/9.0_cp)
               p=3.0_cp*K0*a*(1.0_cp + b*f + c*f*f)

            case (3) ! Vinet
               x=vv0**(-1.0_cp/3.0_cp)
               u=1.0_cp -x
               p=3.0_cp*K0*u/x/x*exp(1.5_cp*(Kp-1.0_cp)*u)

            case (4) ! Natural
               b=0.0_cp
               c=0.0_cp
               if (EoS%iorder > 2)  b=1.5_cp*(Kp-2.0_cp)
               if (EoS%iorder == 4) c =1.5_cp*(K0*Kpp + 1.0_cp +(Kp-2.0_cp)+(Kp-2.0_cp)**2.0_cp)
               p=3.0_cp*vv0*K0*f*(1.0_cp + b*f + c*f*f)

            case (5) ! Tait
               abc= get_tait(t,EoS)
               vv0=1.0_cp/vv0     ! back to vv0=v/v0
               p=(((vv0 +abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3)) - 1.0_cp)/abc(2)

            case(6) ! APL
               apl= Get_APL(1.0_cp/VV0,vv0*vol,K0,Kp,Kpp,ev(5),EoS%iorder)
               p=3.0_cp*k0*product(apl(1,:))

            case(7) !Kumar
               a=kp+1.0_cp
               vv0=1.0_cp/vv0     ! back to vv0=v/v0
               p=k0/a * (exp(a*(1.0_cp-vv0))-1)

         end select

         !> Handle pthermal EoS
         if (EoS%pthermaleos .and. EoS%imodel /= 0) p=p+pthermal(Vol,T,EoS)        !Vol is a or V

         !> Iteration required if there is a phase transition
         if (EoS%itran == 0) exit

         !> Determine and set the transition pressure at this T
         if (first) then
            i=0
            if (EoS%itran == 1)then
               ptr=ev(21)       !PV transition, so Ptr is the eos%param(21)

            else
               !> ptr=(T-ev(21))/ev(22)
               ptr=get_transition_pressure(T,EoS)
            end if
            vtr=huge(0.0_cp)
         end if

         i=i+1

         if (abs(plast-p) < 0.0001) exit         ! iteration has converged
         if (i > 1000)then
            err_eos=.true.
            err_eos_mess='No convergence in Transition loop in Get_Pressure'
            return
         endif

         !> Here if not converged:               ! works provided phase boundary is almost linear
         if (first .or. abs(plast-p) < difp) then
            vs=get_transition_strain(p,T,EoS)
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
   Function Get_Pressure_Esd(V, T, Eos) Result(Esd)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V    ! Volume
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: esd

      !---- Local Variables ----!
      integer                           :: i,j
      real(kind=cp),dimension(n_eospar) :: td
      real(kind=cp)                     :: vol,temp
      type(Eos_Type)                    :: E  ! Eos Parameter local copy

      !> Init
      esd=0.0_cp

      !> local copies
      vol=v
      temp=t
      e=eos
      td= Deriv_Partial_P(vol,temp,e,xtype=0,calc='ALL')  ! gets all derivatives dP/dparam in array td

      !> esd
      do i=1,N_EOSPAR
         do j=1,N_EOSPAR
            esd=esd+eos%vcv(i,j)*td(i)*td(j)
         end do
      end do

      !> Final
      esd=sqrt(esd)

      return
   End Function Get_Pressure_Esd

   !!--++
   !!--++ FUNCTION GET_PRESSURE_K
   !!--++   returns the pressure for a known bulk modulus K and temperature T
   !!--++
   !!--++ Date: 16/03/2017
   !!
   Function Get_Pressure_K(K, T, Eos, Itype, Pest) Result(P)
      !---- Arguments ----!
      real(kind=cp),           intent(in) :: K       ! Bulk modulus
      real(kind=cp),           intent(in) :: T       ! Temperature
      type(Eos_Type),          intent(in) :: EoS     ! Eos Parameter
      integer,                 intent(in) :: itype
      real(kind=cp), optional, intent(in) :: Pest    ! approx pressure: needed if transitions
      real(kind=cp)                       :: p

      !---- Local Variables ----!
      integer            :: ic,irev
      real(kind=cp)      :: temp,ktarget,kcalc,ktol,step,del,delprev,pprev,v0,step_prev
      real(kind=cp)      :: pmindiff
      real(kind=cp)      :: k0,steptest
      !character(len=255) :: ltext

      !> init
      ic=0
      irev=0      ! reversal counter
      temp=t
      ktarget=k

      !> initial guess...use Pest if present
      if (present(Pest)) then
         p=Pest

      else
         !> Murngahan guess based on T=Tref
         V0=get_volume(0.0_cp,eos%Tref,eos)
         K0=k_cal(V0,eos%Tref,eos)
         p=(ktarget-K0)/eos%params(3)
      end if

      pmindiff=p                  !> record p of best fit
      step=eos%params(2)/100.0    !> make initial step 0.01 K0 of first phase: not critical
      ktol=0.001*ktarget

      do
          !> counter
          ic=ic+1

          kcalc=get_property_x(p,t,eos,itype)
          del=ktarget-kcalc
          if (ic > 1) then                      ! need to get two calcs before adjusting step size and dirn
             if (abs(del) < ktol .or. abs(p-pprev) < 0.0005*P) return
             if (abs(step) < 0.0001) return                             ! step size getting down to precision

             if (abs(del) < abs(delprev) ) pmindiff=p
             if (abs(del-delprev) > 10.0*tiny(1.0_cp) ) then
                steptest= del/(del-delprev) *(pprev-p)  ! linear approx: newstep= dx/dy * (ytarget-yc)
                if (steptest < step) then
                   step=steptest
                else
                   step=sign(step,steptest)
                end if
             end if
          end if

          pprev=p          ! this x stored
          delprev=del      ! ytarget-yc stored
          if (eos%itran > 0) then
             !> trap jumping back and forth across phase boundary
             if (transition_phase(P,Temp,eos) .neqv. transition_phase(P+step,Temp,eos) ) then
                step=0.9*Get_Transition_Pressure(Temp,Eos)-p
                !write(ltext,'(''  Step limited by transition pressure of '',f7.4)')Get_Transition_Pressure(Temp,Eos)
                !if (idmode == 2)call write_out(ltext)
             end if
          end if

          if (ic > 2 .and. abs(step) > abs(step_prev) ) then
             step=sign(0.9*step_prev,step)
          end if

          !write(ltext,'('' P ='',f7.4,'' Step = '',f7.4,''  Kcalc ='',f10.5,'' Ktarget = '',f10.5,'' Ktol ='',f10.5)')p,step,kcalc,ktarget,ktol
          !if (idmode == 2)call write_out(ltext)

          p=p+step
          step_prev=step

          if (ic > 100) then
             !if(idmode == 2)call write_out('********Convergence failure in get_pressure_K')
             !write(ltext,'('' P ='',f7.4,'' Step = '',f7.4,''  Kcalc ='',f10.5,'' Ktarget = '',f10.5,'' Ktol ='',f10.5)')p,step,kcalc,ktarget,ktol

             !if(idmode == 2)call write_out(ltext)
             !if(idmode == 2)call write_out_blank(1)
             warn_eos=.true.
             warn_eos_mess='Convergence failure in get_pressure_K'
             if (eos%itran > 0) warn_eos_mess=trim(warn_eos_mess)//' probably due to transition giving max in K'
             p=pmindiff
             return
          end if
      end do

      return
   End Function Get_Pressure_K

   !!----
   !!---- FUNCTION GET_PRESSURE_X
   !!----
   !!---- Gets ...
   !!----
   !!---- Date: 16/03/2017
   !!
   Function Get_Pressure_X(X, T, Eos, Xtype, Pest) Result(P)
      !---- Arguments ----!
      real(kind=cp),           intent(in) :: x       ! V or Bulk modulus
      real(kind=cp),           intent(in) :: T       ! Temperature
      type(Eos_Type),          intent(in) :: EoS     ! Eos Parameter
      integer,                 intent(in) :: Xtype   ! =0 when X=V, =1 for X=K (isothermal)  =2 for adiabatic
      real(kind=cp), optional, intent(in) :: Pest    ! Approx pressure: needed if transitions
      real(kind=cp)                       :: p

      !---- Local Variables ----!
      integer :: itype

      !> Init
      p=0.0_cp

      itype=xtype
      if (itype < 0 .or. itype > N_DATA_TYPES) itype=0               ! default

      select case(itype) ! Get calc pressure: depends on data type
         case(0)                 ! X is V
            p=get_pressure(x,t,eos)

         case(1)                 ! X is KT
            p=get_pressure_K(x,t,eos,itype,Pest)

         case(2)                 ! X is Ks:
            p=get_pressure_K(x,t,eos,itype,Pest)

         case default
            p=0.0_cp
      end select

      return
   End Function Get_Pressure_X

   !!----
   !!---- FUNCTION GET_PROPERTY_X
   !!----
   !!---- Returns the ...
   !!----
   !!---- Date: 16/03/2017
   !!
   Function Get_Property_X(P, T, Eos, Xtype) Result(Val)
      !---- Arguments ----!
      real(kind=cp),     intent(in) :: P       ! Pressure
      real(kind=cp),     intent(in) :: T       ! Temperature
      type(Eos_Type),    intent(in) :: EoS     ! Eos Parameter
      integer, optional, intent(in) :: xtype   ! =0 when X=V, =1 for X=K (isothermal)
      real(kind=cp)                 :: Val

      !---- Local Variables ----!
      integer       :: itype
      real(kind=cp) :: vol, agt

      !> Init
      itype=0               ! default
      if (present(xtype)) itype=xtype

      select case(itype)
         case(0)                 ! Volume
            val=Get_Volume(P,T,Eos)

         case(1)                 ! Isothermal modulus
            vol=Get_Volume(P,T,Eos)
            val=K_Cal(Vol,T,Eos,p=p)

         case(2)                 ! Adiabatic modulus
            vol=Get_Volume(P,T,Eos)
            agt=Alpha_Cal(P,T,Eos)*get_grun_th(P,T,eos)*T        ! Get_Grun knows about linear
            if (eos%linear) agt=3.0_cp*agt
            val=(1.0_cp+agt)*K_Cal(Vol,T,Eos,p=p)

         case default
            val=0.0
      end select

      return
   End Function Get_Property_X

   !!----
   !!---- FUNCTION GET_PROPS_GENERAL
   !!----
   !!---- Returns elastic properties input P,T for a general direction axis in cell_eos
   !!---- Use Eos_Cal when the eos for an axis is known
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Props_General(P, T, Cell_Eos, Axis) Result(Parvals)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp), dimension(6)     :: parvals

      !---- Local Variables ----!

      !>init
      parvals=0.0_cp

      parvals(1)=get_Volume_general(P,T,cell_eos,axis)      ! length
      parvals(2)=get_mod_general(P,T,cell_eos,axis)         ! M
      parvals(3)=get_modp_general(P,T,cell_eos,axis)        ! Kp
      parvals(5)=get_dMdT_general(P,T,cell_eos,axis)        ! dK/dT
      parvals(6)=get_alpha_general(P,T,cell_eos,axis)       ! alpha

      return
   End Function Get_Props_General

   !!--++
   !!--++ FUNCTION GET_PROPS_PTVTABLE
   !!--++
   !!--++ Returns the requested property at two of P,T,V from an eos loaded as a table of values
   !!--++ The volume values in the table are expected to be scaled to V=1.0 at Pref, Tref
   !!--++ and the true V0 is in eospar%params(1)
   !!--++
   !!--++ Date: 14/10/16
   !!
   Function Get_Props_PTVTable(P, T, V, Eos, Res) Result(Val)
      !---- Arguments ----!
      real(kind=cp),    intent(in) :: P       ! Pressure
      real(kind=cp),    intent(in) :: T       ! Temperature
      real(kind=cp),    intent(in) :: V       ! Volume
      type(Eos_Type),   intent(in) :: EoS     ! Eos Parameter
      character(len=*), intent(in) :: Res     ! Parameter requested for calculation (P,T, or V)
      real(kind=cp)                :: Val

      !---- Local Variables ----!
      integer                   :: i,j
      real(kind=cp)             :: VV0,tt,vm,vp,km,kp,va,vb,dvdp
      type(Eos_Type)            :: eosm,eosp,eosv     ! Eos for local Murn at Tminus and Tplus
      character(len=10)         :: Var
      !character(len=60)         :: text

      !>Init
      Val=0.0_cp

      !> Check valid request
      if (EoS%imodel/= -1) then
         err_eos=.true.
         write(err_eos_mess,'(''Request to get_props_pvttable with invalid imodel #'',i5)')EoS%imodel
         return
      end if

      VV0=V/EoS%params(1)        ! table values are all V/V0

      !> Determine what is the request: P,T,V, K, KP, AL
      Var=Res
      Var=U_Case(adjustl(Var))

      !> Find lower-corner coords for P,T  : works because P and T always ascend
      if (Var(1:1) /= 'P')then                         ! Only test if P supplied in argument, not if P requested from T and V
         do i=1,EoS%table%np
            if (p < EoS%table%ptv(i,1,1)) exit
         end do

         if (i < 2)then
            warn_eos=.true.
            write(warn_eos_mess, &
                 '(''PTV table request for P = '',f6.2,'' smaller than Pmin = '',f6.2,'' at T ='',f6.1,'': Linear guess made'')') &
                 p,EoS%table%pmin,t
            i=2
         end if

         if (i > EoS%table%np-1) then
            warn_eos=.true.
            write(warn_eos_mess,&
                 '(''PTV table request for P = '',f6.2,'' bigger than Pmax = '',f6.2,'' at T ='',f6.1,'': Linear guess made'')') &
                 p,EoS%table%pmax,t
            i=EoS%table%np-1
         end if
      end if

      !> Check the T limits
      do j=1,EoS%table%nt
         if (t < EoS%table%ptv(1,j,2) ) exit
      end do

      if (j < 2 ) then
         warn_eos=.true.
         write(warn_eos_mess, &
              '(''PTV table request for T = '',f6.1,'' smaller than Tmin = '',f6.1,'' at P ='',f6.2,'': Linear guess made'')') &
              t,EoS%table%tmin,p
         j=2
      end if

      if (j > EoS%table%nt-1) then
         warn_eos=.true.
         write(warn_eos_mess, &
              '(''PTV table request for T = '',f6.1,'' bigger than Tmax = '',f6.1,'' at P ='',f6.1,'': Linear guess made'')') &
              t,EoS%table%tmax,p
         j = EoS%table%nt-1
      end if

      !> Get the Murngahan EoS from Pminus at Tminus and Tplus, but set as if a T=Tref
      !> Must then call Eos functions with Tref as argument
      if (Var(1:1) /= 'P')then       ! cannot do this if P being requested
         EoSM= Murn_Ptv_table(i-1,j-1,EoS)
         EoSp= Murn_Ptv_table(i-1,j,EoS)
      end if

      Select case(Var(1:3))
         case('V  ')
            !> Murnaghan interpolation on P axis, followed by linear in T
            Vm=get_volume(p-EoS%table%ptv(i-1,j-1,1),eosm%tref,eosm)
            Vp=get_volume(p-EoS%table%ptv(i-1,j,1),eosp%tref,eosp)
            val=Vm+(Vp-Vm)*(t-EoS%table%ptv(1,j-1,2))/(EoS%table%ptv(1,j,2)-EoS%table%ptv(1,j-1,2))  ! linear T
            val=val*EoS%params(1)

         case('P  ')
            tt=(t-EoS%table%ptv(1,j-1,2))/(EoS%table%ptv(1,j,2)-EoS%table%ptv(1,j-1,2))
            do i=1,EoS%table%np
               vb=va
               va=tt*(EoS%table%ptv(i,j,3)-EoS%table%ptv(i,j-1,3)) +EoS%table%ptv(i,j-1,3)      ! volume at each P row for target T
               if (va < vv0) exit
            end do

            if (i == 1) then
               warn_eos=.true.
               write(warn_eos_mess, &
                    '(''PTV table request for V = '',f6.1,''at T = '',f6.1,'' is below Pmin: Linear guess made'')') v,t

               !> linear guess beyond table bottom
               dvdp=(tt*(EoS%table%ptv(2,j,3)-EoS%table%ptv(2,j-1,3)) +EoS%table%ptv(2,j-1,3)-va)/(EoS%table%ptv(2,j,1)- &
                    EoS%table%ptv(1,j,1))
               val=EoS%table%ptv(1,j,1)- (va-vv0)/dvdp
               return
            end if

            if (i >= EoS%table%np) then
               warn_eos=.true.
               write(warn_eos_mess, &
                    '(''PTV table request for V = '',f6.1,''at T = '',f6.1,'' is above Pmax: Linear guess made'')') v,t
               !> linear guess beyond table bottom
               dvdp=(va-vb)/(EoS%table%ptv(EoS%table%np,j,1)-EoS%table%ptv(EoS%table%np-1,j,1))
               val=EoS%table%ptv(EoS%table%np,j,1)+ (vv0-va)/dvdp
               return
            end if

            !> Set up Murn in Eosm
            call Init_EoS_Type(Eosv)
            eosv%imodel=1
            call set_eos_names(eosv)      ! sets names for params
            eosv%title='Murnaghan for interpolation'
            eosv%params(1)=vb
            EoSM= Murn_Ptv_table(i-1,j-1,EoS)
            EoSP= Murn_Ptv_table(i-1,j,EoS)
            eosv%params(2)=eosm%params(2)+(eosp%params(2)-eosm%params(2))*tt
            eosv%params(3)=eosm%params(3)+(eosp%params(3)-eosm%params(3))*tt
            val=get_pressure(vv0,eosv%tref,eosv)+EoS%table%ptv(i-1,j,1)

         case('K  ') !input must be p and t
            km=get_k(p-EoS%table%ptv(i-1,j-1,1),eosm%tref,eosm)      ! K at Tminus
            kp=get_k(p-EoS%table%ptv(i-1,j,1),eosp%tref,eosp)          ! K at Tplus
            tt=(t-EoS%table%ptv(1,j-1,2))/(EoS%table%ptv(1,j,2)-EoS%table%ptv(1,j-1,2))
            val=(kp-km)*tt + km

         case('KP ')
            tt=(t-EoS%table%ptv(1,j-1,2))/(EoS%table%ptv(1,j,2)-EoS%table%ptv(1,j-1,2))
            val=eosm%params(3)+(eosp%params(3)-eosm%params(3))*tt

         case('AL ')
            !> Murnaghan interpolation on P axis, followed by linear in T
            Vm=get_volume(p-EoS%table%ptv(i-1,j-1,1),eosm%tref,eosm)
            Vp=get_volume(p-EoS%table%ptv(i-1,j,1),eosp%tref,eosp)
            val=2.0_cp*(Vp-Vm)/(EoS%table%ptv(i-1,j,2)-EoS%table%ptv(i-1,j-1,2))/(Vm+Vp)

         case default
            err_eos=.true.
            write(err_eos_mess,'(''Request for '',a1,'' to get_props_pvttable not valid'')')Var(1:1)

      end select

      return
   End Function Get_Props_PTVTable

   !!----
   !!---- FUNCTION GET_PROPS_THIRD
   !!----
   !!---- Returns elastic properties input P,T for principal axis ieos in cell_eos
   !!---- when it can calculated from other eos
   !!---- Use Eos_Cal when the eos for an axis is known
   !!----
   !!---- Date: 09/09/2020
   !!----
   Function Get_Props_Third(P, T, Cell_Eos, Ieos) Result(Parvals)
      !---- Arguments ----!
      real(kind=cp),             intent(in)  :: p,T
      type(eos_cell_type),       intent(in)  :: cell_eos
      integer,                   intent(in)  :: ieos     ! the axis (1,2,3) or V (0) to be calculated
      real(kind=cp), dimension(6)            :: parvals

      !---- Local Variables ----!

      !>init
      parvals=0.0_cp

      parvals(1)=get_volume_third(P,T,cell_eos,ieos)       ! V
      parvals(2)=get_mod_third(P,T,cell_eos,ieos)          ! K
      parvals(3)=get_modp_third(P,T,cell_eos,ieos)         ! Kp
      parvals(5)=get_dMdT_third(P,T,cell_eos,ieos)         ! dK/dT
      parvals(6)=get_alpha_third(P,T,cell_eos,ieos)        ! alpha

      return
   End Function Get_Props_Third

   !!--++
   !!--++ SUBROUTINE GET_TAIT
   !!--++
   !!--++ Returns a,b,c Tait parameters in a vector.
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Function Get_Tait(T, Eos) Result(Vec)
      !---- Arguments ----!
      real(kind=cp),             intent(in)  :: T    ! Temperature
      type(Eos_Type),            intent(in)  :: Eos  ! Eos Parameters
      real(kind=cp), dimension(3)            :: Vec  ! Vector (a,b,c) of Tait parameters

      !---- Local Variables ----!
      real(kind=cp) :: k0,kp,kpp

      !> Init
      Vec=0.0_cp
      if (eos%imodel /= 5) return

      select case(eos%itherm)
         case (1:5)                  ! normal thermal expansion models with dK/dT
            k0 =Get_K0_T(T,eos)
            kp =Get_Kp0_T(T,eos)

         case default               ! includes no thermal model, also pthermal which requires params at Tref
            k0 =eos%params(2)
            kp =eos%params(3)
      end select

      if (eos%iorder < 4) then
         kpp= -1.0_cp*kp/k0        ! implied value for Kpp except for 4th order
      else
         kpp=eos%params(4)
      end if

      if (eos%linear) then
         k0  =k0/3.0_cp
         kp  =kp/3.0_cp
         kpp =kpp/3.0_cp
      end if

      Vec(1)=(1.0_cp + kp)/(1.0_cp+kp+k0*kpp)
      Vec(2)= kp/k0
      if (abs(1_cp+kp) > tiny(0.0)) Vec(2)=Vec(2)-kpp/(1.0_cp+kp)

      Vec(3)= (1.0_cp+kp+k0*kpp)/(kp*kp+kp-k0*kpp)

      return
   End Function Get_Tait

   !!----
   !!---- FUNCTION GET_TEMPERATURE
   !!----
   !!---- Returns Temperature at given P,V
   !!----
   !!---- Date: 01/04/2014
   !!
   Function Get_Temperature(P, V, Eos) Result(T)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P    ! Pressure
      real(kind=cp),  intent(in) :: V    ! Volume
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter

      !---- Local Variables ----!
      integer                           :: nstep
      real(kind=cp)                     :: t
      real(kind=cp)                     :: pa,va,ta
      real(kind=cp)                     :: step,dp1,dp2

      !> Init
      t=EoS%tref
      pa=p               ! local copy p
      va=v

      !> Check
      if (EoS%itherm == 0) return          ! no calcs possible without thermal model

      !> First estimate at P=0
      t=Get_Temperature_P0(va,EoS)
      if (EoS%imodel ==0) return
      if (err_eos) return

      !> Use iterative solution, by minimising p-pcalc
      !> Init this part
      ta=t
      dp1=p-get_pressure(va,ta,EoS)
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
         dp2=p-get_pressure(va,ta,EoS)
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

   !!--++
   !!--++ FUNCTION GET_TEMPERATURE_P0
   !!--++
   !!--++ Returns Temperature at P=0 for given V0T
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Function Get_Temperature_P0(V, EoS, Tmin, Tmax) Result(Tk)
      !---- Arguments ----!
      real(kind=cp),           intent(in) :: V       ! Volume at temperature T or Pth (Pthermal case)
      type(Eos_Type),          intent(in) :: EoS     ! Eos Parameter
      real(kind=cp), optional, intent(in) :: Tmin    ! Range for solution in T
      real(kind=cp), optional, intent(in) :: Tmax
      real(kind=cp)                       :: tk

      !---- Local Variables ----!
      real(kind=cp)                      :: tref
      real(kind=cp)                      :: v00,v0T
      real(kind=cp)                      :: a,b,c,d,t1,t2,t3,x1,x2,x3
      real(kind=cp)                      :: a1,a2,a3,q,r,s1,s2,dd,th
      real(kind=cp)                      :: y,p
      real(kind=cp)                      :: kp,th_e,eps0
      real(kind=cp)                      :: Tkmin,Tkmax
      real(kind=cp), dimension(N_EOSPAR) :: ev

      !> Init
      Tk=EoS%tref

      !> Local copy EoS to handle linear or volume
      ev= EoS_to_Vec(EoS)

      !> all equations written for volume
      v00=ev(1)
      if (EoS%linear) then
         v0T=V**3.0_cp
      else
         v0T=V
      end if

      !> Optional arguments
      TKmin=0.0_cp
      if (present(Tmin)) TKmin=tmin

      TKmax=1.0e4
      if (present(Tmax)) TKmax=tmax

      !> Check
      Tref=EoS%tref
      if (v00 <= tiny(0.0)) return

      !> Init local variables
      t1=0.0
      t2=0.0
      t3=0.0

      select case (EoS%itherm)
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

                  elseif (t2 > 0.0) then
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
               else if (x2 > 0.0) then
                  t2=x2*x2
                  if (t2 >=TKmin .and. t2 <=TKmax) tk=t2
               end if
            else
               err_eos=.true.
               err_eos_mess='No valid solution for temperature'
            end if

         case (4) ! Kroll
            !>>>>> kp=ev(3) : version before 11/11/2016
            if (EoS%icross == 2)then
               kp=ev(8)                    !uses Anderson delta
            else
               kp=ev(3)
            end if

            th_e=ev(11)
            a=th_e/Tref
            b=-1.0/(kp*(kp+2.0))
            eps0=(a**2)*exp(a)/(exp(a)-1.0)**2

            x1=(-b*(1.0+kp))*eps0/(ev(10)*th_e)
            x2=1.0 - (((v0t/v00)+kp)/(1.0+kp))**(1.0/b)
            x3=1.0/(exp(a)-1)

            c=1.0/(x1*x2+x3)
            if (c <= -1.0_cp) then
               err_eos=.true.
               err_eos_mess='Temperature calculated as negative'
            else
               d=log(1.0 +c)
               t1=th_e/d
               if (t1 >=TKmin .and. t1 <=TKmax) then
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

         case(7,8)  !For MGD Pthermal, no real meaningful value of T to return, so set as Tref
            tk=EoS%tref

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
   Function Get_Transition_Pressure(T, Eos) Result(Ptr)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: Ptr  ! The transition T at this P

      !---- Local Variables ----!
      real(kind=cp) :: sqroot

      !> Init
      ptr=0.0_cp

      !> Check for valid model number. If not valid, return with zero Tr
      if (EoS%itran < 1 .or. EoS%itran > N_TRANS_MODELS) return

      !> Now determine if in high-symmetry or low field
      select case(EoS%itran)
         case(1) ! Landau PV
            ptr=EoS%params(21)

         case(3) ! Landau PVT
            !ptr = (T-EoS%params(21))/EoS%params(22) original linear
            if (abs(EoS%params(23)) < tiny(0.0) ) then
               ptr = (T-EoS%params(21))/EoS%params(22)        ! linear
            else
               sqroot=EoS%params(22)*EoS%params(22)-4.0_cp*EoS%params(23)*(EoS%params(21)-T)
               if (sqroot > tiny(0.0) ) then
                  ptr=2.0_cp*(T-EoS%params(21))/(EoS%params(22)+sqrt(sqroot))        ! Viet/Muller formula for root
               else if(sqroot > -1.0*tiny(0.0) ) then
                  ptr=-2.0_cp*(T-EoS%params(21))/EoS%params(22)
               else
                  Err_EoS=.true.
                  write(Err_EoS_Mess,'(a,f8.3,a)')'No real solution for Ptr at T = ',T,'K'
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
   Function Get_Transition_Strain(P, T, Eos) Result(VS)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P    ! Pressure
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: Vs   ! The volume strain

      !----Local Variables ----!
      real(kind=cp) :: Ttr , a ! transition temperature at this pressure

      !> init
      vs=0._cp

      !> Check for valid model number. If not valid, return with zero Vs
      if (EoS%itran < 1 .or. EoS%itran > N_TRANS_MODELS) return

      !> Now determine if in high-symmetry or low field
      if (Transition_phase(P,T,EoS)) then
         !> This section for being in the low field
         select case(EoS%itran)
            case(1) ! Landau PV
               vs=EoS%params(24)*abs(EoS%params(21)-P)**EoS%params(25)

            case(2) ! Landau TV
               vs=EoS%params(24)*abs(EoS%params(21)-T)**EoS%params(25)

            case(3) ! Landau PVT
               Ttr = Get_Transition_Temperature(P,EoS)
               a=EoS%params(24)
               vs=a*abs(Ttr-T)**EoS%params(25)            !abs function to handle highT being low sym
         end select

      else
         !> This section for being in the high field
         select case(EoS%itran)
            case(1) ! Landau PV
               vs=EoS%params(26)*abs(EoS%params(21)-P)**EoS%params(27)

            case(2) ! Landau TV
               vs=EoS%params(26)*abs(EoS%params(21)-T)**EoS%params(27)

            case(3) ! Landau PVT:  Note no da/dP for highP phase
               Ttr = Get_Transition_Temperature(P,EoS)
               vs=EoS%params(26)*abs(Ttr-T)**EoS%params(27)            !abs function to handle highT being low sym
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
   Function Get_Transition_Temperature(P, Eos) Result(Tr)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P    ! Pressure
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: Tr   ! The transition T at this P

      !---- Local Variables ----!

      !>init
      tr=0._cp

      !> Check for valid model number. If not valid, return with zero Tr
      if (eos%itran < 1 .or. eos%itran > N_TRANS_MODELS) return

      !> Now determine if in high-symmetry or low field
      select case(eos%itran)
         case(2) ! Landau TV
            Tr=eos%params(21)

         case(3) ! Landau PVT: with a curved phase boundary
            Tr = eos%params(21)+p*eos%params(22)+p*p*eos%params(23)
      end select

      return
   End Function Get_Transition_Temperature

   !!--++
   !!--++ FUNCTION GET_V0_AXIS
   !!--++
   !!--++ Returns the value of volume or length of principal axis (ieos) in unit cell
   !!--++ in cell_eos at Pref,Tref
   !!--++
   !!--++ Call this Function directly when the calling routine  knows that the direction
   !!--++ is a principal axis
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_V0_Axis(Cell_eos, Ieos) result(L)
      !---- Arguments ----!
      type(eos_cell_type), intent(in)  :: cell_eos
      integer,             intent(in)  :: ieos      !axis indicator, as in axis_type%ieos
      real(kind=cp)                    :: L         !returned length or volume

      !---- Local Variables ----!

      !> init
      l=10.0_cp

      select case(cell_eos%loaded(ieos))
         case(1)
            L=get_volume(cell_eos%eosc%pref,cell_eos%eosc%tref,cell_eos%eos(ieos))

         case(2) ! sym equiv. Always uses eos(1) for a-axis
            L=get_volume(cell_eos%eosc%pref,cell_eos%eosc%tref,cell_eos%eos(1))

         case(3)
            L=get_volume_third(cell_eos%eosc%pref,cell_eos%eosc%tref,cell_eos,ieos)
      end select

      return
   End Function Get_V0_Axis

   !!--++
   !!--++ FUNCTION GET_V0_CELL
   !!--++
   !!--++ Returns the value of volume or length of any axis in unit cell in cell_eos
   !!--++ at Pref,Tref
   !!--++ Call this Function when the calling routine does not know if the direction
   !!--++ is a principal axis or not
   !!--++ If a principal direction is requested, only axis%ieos is required
   !!--++ axis%v and axis%atype only used if axis%ieos=-2
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_V0_Cell(Cell_eos, Axis) result(L)
      !---- Arguments ----!
      type(eos_cell_type), intent(in)  :: cell_eos
      type(axis_type),     intent(in)  :: axis
      real(kind=cp)                    :: L !returned length or volume

      !---- Local Variables ----!

      !> init
      l=10.0_cp

      select case(axis%ieos)      !invalid numbers just return
         case(0:6)   !principal direction for which eos exists, or can be calculated
            L=get_v0_axis(cell_eos,axis%ieos)

         case(-2)   !general direction
            L=get_volume_general(cell_eos%eosc%pref,cell_eos%eosc%tref,cell_eos,axis)
      end select

      return
   End Function Get_V0_Cell

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
   Function Get_V0_T(T, EoS) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T        ! Temperature
      type(Eos_Type), intent(in) :: EoS   ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: V,kp,AK
      real(kind=cp)                      :: Tref,A,B,C,Tn,tt
      real(kind=cp)                      :: delt,delt2
      real(kind=cp), dimension(n_eospar) :: ev

      !> Init
      Tref=EoS%tref

      !> Local copy EoS
      ev= EoS_to_Vec(EoS) ! Volume or linear case is covered
                                 ! all equations written for volume
      delt=T-Tref
      select case(EoS%itherm)
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
            !>>>>>kp=ev(3)   : version before 11/11/2016
            if (EoS%icross == 2) then
               kp=ev(8)
            else
               kp=ev(3)
            end if
            if(abs(kp-1) < 0.0001 .or. abs(kp/(kp+2.0_cp)) < 0.0001)then
                V=ev(1)                             ! In these cases algebra shows V=V0
            else
                Tn= ev(11)/Tref                        ! theta/Tref
                C=Tn*Tn*exp(Tn)/(exp(tn)-1.0_cp)**2.0_cp
                B=-1.0/kp/(kp+2.0_cp)
                if (t > 0.05_cp*ev(11)) then                               ! to avoid numerical problems at T=0
                   A=ev(10)*ev(11)/C *(1.0_cp/(exp(ev(11)/T)-1.0_cp) - 1.0_cp/(exp(Tn)-1.0_cp) )
                else
                   A=ev(10)*ev(11)/C *(-1.0_cp/(exp(Tn)-1.0_cp) )          ! because when T=0 then 1.0/exp(Tein/T) = 0
                end if

                !  V=ev(1)*(-1.0_cp*kp + (1.0_cp+kp)*(1.0_cp - kp*(kp+2.0_cp)*A/(kp+1.0_cp))**B)
                AK=1.0_cp - kp*(kp+2.0_cp)*A/(kp+1.0_cp)
                if (AK < tiny(0._cp))then
                   V=ev(1)       ! for safe return
                   err_eos=.true.
                   err_eos_mess='T exceeds valid limit for Kroll expansion in get_V0_T'
                else
                   V=ev(1)*(-1.0_cp*kp + (1.0_cp+kp)*AK**B)
                end if
            end if

         case(5)                    ! Salje, Tref fixed at zero
            A=T/ev(11)
            if (A < 0.001) then
               V=ev(1)                 ! ultra-low T: coth(theta_sat/T)=1.0
            else
               A=1.0_cp/tanh(1.0_cp/A) ! the coth(theta_sat/T)
               V=(ev(1)**(1.0_cp/3.0_cp) + ev(10)*ev(11)*(A-1.0_cp))**3.0_cp
            end if

         case(6,7,8)
            v=ev(1)         ! Pthermal needs V0 at Tref
      end select

      !> Linear
      if (EoS%linear) v=v**(1.0_cp/3.0_cp)

      return
   End Function Get_V0_T

   !!----
   !!---- FUNCTION GET_VOLUME
   !!----
   !!---- Find volume from EoS at given P and T
   !!----
   !!---- Uses a spline to solve for volume if no analytical approach available
   !!---- Written 3/2019 RJA
   !!---- Under test
   !!----
   !!---- Date: 03/02/2021
   !!
   Function Get_Volume(P, T, Eos) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P    ! Pressure
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: V

      !---- Local Variables ----!
      real(kind=cp), parameter          :: PREC=0.000001_cp  !precision to find V.

      integer,parameter                 :: nstep=30
      integer                           :: i,ic
      real(kind=cp),dimension(nstep)    :: x,y,d2y

      type(Eos_Type)                    :: EoSc  ! Eos copy
      real(kind=cp)                     :: V0,K0,Kp,strain,vfactor !,k
      real(kind=cp)                     :: Vol, vstep, delp_prev,delp,v_prev,a,logterm
      real(kind=cp),dimension(N_EOSPAR) :: ev
      real(kind=cp),dimension(3)        :: abc          ! Tait parameters
      real(kind=cp)                     :: pa           ! pa=p-pth

      logical                           :: reverse


      !> Init
      v=0.0_cp
      strain=0.0_cp   ! strain from transition: linear or volume to match eos type
      pa=p            ! local copy p

      !> If PTV table, go directly
      if (EoS%imodel == -1) then
         v=get_props_ptvtable(pa,t,0.0,EoS,'V')     ! get_props_ptvtable returns length if linear
         return
      end if

      !> Local copy EoS
      ev= EoS_to_Vec(EoS) ! Volume or linear case is covered

      !> Set appropriate V0, and adjust for thermal pressure:
      select case (EoS%itherm)
         case (0)                                   ! 0=no thermal,
            v0=ev(1)                                  ! v0 is volume eos, (a0)^3 for linear

         case (1:5)
            v0=get_v0_t(t,EoS)                     ! returns a0 for linear
            if (EoS%linear) v0=v0**3.0_cp

         case (6)                                     ! HP or linear thermal pressure
            v0=ev(1)
            pa=p-pthermal(0.0,t,EoS)               ! adjust pressure to isothermal pressure for murn and tait estimates

         case(7,8)                                    ! MGD - do a guess on basis parallel isochors of pa
            v0=ev(1)
            pa=p - EoS%params(2)*(t-EoS%tref)/100000.         ! have to guess an alpha because that requires V !!!
      end select

      !> set K0, for various purposes
      k0=Get_K0_T(T,EoS)              ! Handles thermal pressure case, returns K0 or M0
      if (EoS%linear) k0=k0/3.0_cp

      kp=Get_Kp0_T(T,EoS)
      if (EoS%linear) kp=kp/3.0_cp

      !> Get the volume strain due to transition: only a function of P,T NOT V!!
      if (EoS%itran > 0) then
         strain=get_transition_strain(P,T,EoS)     ! returns the linear or volume strain as appropriate
      end if

      !> If there is no eos model, we are finished because V0 is the V at this T
      if (EoS%imodel == 0) then
         if (EoS%linear) v0=v0**(1.0_cp/3.0_cp)
         v=v0*(1.0_cp + strain)
         return
      end if

      !> Analytic solution for Murnaghan:  use this for first guess for other EoS except Tait
      vfactor=(1.0_cp + kp*pa/k0)
      if (vfactor < 0.0)then
         v=v0        ! safe value for when 1+kp*pa/k0 is negative
      else
         v=v0*(vfactor)**(-1.0_cp/kp)
      end if

      !> Cannot do the following if MGD pthermal
      if (EoS%itherm /=7  .and. EoS%itherm /=8) then
         if (EoS%imodel ==1) then
            !> Exact solution for Murnaghan
            if (EoS%linear) v=v**(1.0_cp/3.0_cp)
            if (EoS%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
            !if (vfactor < 0.0 .and. kp > 0.0) then
            if (vfactor < 0.0) then
               err_eos=.true.
               err_eos_mess='Pressure < -K0/Kp: yields meaningless volumes for Murnaghhan EoS'
            end if
            return
         end if

         !> Analytic solution for Tait
         if (EoS%imodel ==5) then
            abc= get_tait(t,EoS)                     ! get_tait returns volume-like parameters even for linear
            if (abc(2)*pa < -0.999999_cp) then
               err_eos=.true.
               err_eos_mess='Pressure yields infinite volume for Tait Eos'
               v=9999.0        ! safe value return
            else
               v=v0*(1.0_cp-abc(1)*(1.0_cp-(1.0_cp + abc(2)*pa)**(-1.0_cp*abc(3))))
               if (EoS%linear) v=v**(1.0_cp/3.0_cp)
               if (EoS%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
            end if
            return
         end if

         !> Analytic solution for Kumar
         if (EoS%imodel ==7) then
             a=kp+1.0_cp
             logterm=a*pa/k0 +1.0_cp    !This is for safety: we should have checked with physical check
             if (logterm < tiny(0._cp))then
                err_eos=.true.
                return
            else
               v=v0*(1.0_cp - log(logterm)/a)
               if (EoS%linear) v=v**(1.0_cp/3.0_cp)
               if (EoS%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
            end if
            return
         end if
      end if

      !> Find iterative solution for the rest of functions: get_pressure includes the thermal pressure term
      !> But if there is a transition, we only want the P/V for the bare high-symm phase without softening

      !From here work with Vol, starting from Murnaghan estimate
      vol=v
      if (EoS%linear) vol=vol**(1.0_cp/3.0_cp)
      eosc=EoS        ! copy
      eosc%itran=0       ! turn off transition

      !> initial simple hunt
      delp_prev=huge(0._cp)
      ic = 0
      reverse=.false.
      Vstep=eosc%params(1)/100._cp
      do
         ic=ic+1
         if (ic > 1000)then
            err_eos=.true.
            err_eos_mess=' *****No solution found in get_volume after 1000 cycles'
            return
         end if

         call init_err_eos()                   ! have to clear the previous errors, otherwise get_pressure will return 0
         delp=p-get_pressure(Vol,T,eosc)
         if (abs(delp) < 0.000001_cp)then  ! hit correct vol by accident. Happens if P=0 at Tref for MGD
            v=vol
            if (EoS%itran > 0) v=vol*(1.0_cp + strain)
            return
         end if

         if (delp*delp_prev < 0._cp .and. ic > 1)then     ! over-stepped solution: solution between v_prev and v
            vol=vol-delp*Vstep/(delp-delp_prev)                               ! best guess
            exit
         end if

         if (abs(delp) > abs(delp_prev))then ! delta-pressure getting bigger
            if (reverse)then               ! found a minimum between v_prev-vstep and v
               err_eos=.true.
               err_eos_mess=' *****No volume found in get_volume'
               return
            else
               reverse=.true.            ! just going the wrong way
               vstep=-1.0_cp*vstep
            end if
         end if
         v_prev=vol         ! this volume stored
         delp_prev=delp  ! store delp
         Vol=Vol+Vstep
      end do

      ! now calculate PV around the solution: we want increasing P, so this means vstep < 0
      Vstep=-1.0_cp*abs(Vstep)
      Vol=Vol-2.0*Vstep
      Vstep=4.0*Vstep/nstep

      do i=1,nstep
         x(i)=get_pressure(vol,t,eosc)
         y(i)=vol
         vol=vol+vstep
      end do
      call Second_Derivative(x, y, nstep, d2y)
      call splint(x,y,d2y,nstep,p,vol)

      v=vol
      if (EoS%itran > 0) v=vol*(1.0_cp + strain)  ! apply transition strain ('vol' is actually linear if linear eos)

      return
   End Function Get_Volume

   !!----
   !!---- FUNCTION GET_VOLUME_AXIS
   !!----
   !!---- Returns the value of volume or length of principal axis (ieos) in unit cell
   !!---- in cell_eos at P,T
   !!---- Call this Function directly when the calling routine  knows that the direction is a principal axis
   !!----
   !!---- Date: 09/09/2020
   !!
   Function Get_Volume_Axis(P, T, Cell_Eos, Ieos) Result(L)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      integer,            intent(in)  :: ieos      !axis indicator, as in axis_type%ieos
      real(kind=cp)                   :: L !returned length or volume

      !---- Local Variables ----!

      !> init
      l=10.0_cp

      select case(cell_eos%loaded(ieos))
         case(1)
            L=get_volume(p,t,cell_eos%eos(ieos))

         case(2) ! sym equiv. Always uses eos(1) for a-axis
            L=get_volume(p,t,cell_eos%eos(1))

         case(3)
            L=get_volume_third(p,T,cell_eos,ieos)

         case(4)     ! only in mono  ...the d-sapcing of the unique axis
            L=get_volume(p,t,cell_eos%eos(cell_eos%unique))
      end select

      return
   End Function Get_Volume_Axis

   !!----
   !!---- FUNCTION GET_VOLUME_CELL
   !!----
   !!---- Returns the value of volume or length of any axis in unit cell in
   !!---- cell_eos at P,T
   !!---- Call this Function when the calling routine does not know if the direction
   !!---- is a principal axis or not
   !!---- If a principal direction is requested, only axis%ieos is required
   !!---- axis%v and axis%atype only used if axis%ieos=-2
   !!----
   !!---- Date:  09/09/2020
   !!
   Function Get_Volume_Cell(P, T, Cell_Eos, Axis) Result(L)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: L !returned length or volume

      !---- Local Variables ----!

      !> init
      l=10.0_cp

      select case(axis%ieos)      !invalid numbers just return
         case(0:6)   !principal direction for which eos exists, or can be calculated
            L=get_volume_axis(p,t,cell_eos,axis%ieos)

         case(-2)   !general direction
            L=get_volume_general(p,T,cell_eos,axis)
      end select

      return
   End Function Get_Volume_Cell

   !!--++
   !!--++ FUNCTION GET_VOLUME_GENERAL
   !!--++
   !!--++ Returns the value of volume or length of any axis in unit cell in cell_eos at P,T
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Volume_General(P, T, Cell_Eos, Axis) Result(L)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  :: p,T
      type(eos_cell_type),intent(in)  :: cell_eos
      type(axis_type),    intent(in)  :: axis
      real(kind=cp)                   :: L !returned length

      !---- Local Variables ----!
      type(crystal_cell_type) :: cell

      !> get the unit cell, metric tensors
      cell= get_params_cell(P,T,cell_eos)

      !> calculate the distance
      select case(U_case(axis%atype))
         case('U')
            L= dot_product(axis%v,matmul(cell%gd,axis%v))
            L=sqrt(abs(L))

         case('H')
            L= 1.0_cp/dot_product(axis%v,matmul(cell%gr,axis%v))
            L=sqrt(abs(L))

         case default
            L=10.0_cp
      end select

      return
   End Function Get_Volume_General

   !!--++
   !!--++ FUNCTION GET_VOLUME_K
   !!--++
   !!--++ Returns the value of Volume for a given K  at T without using pressure
   !!--++ This has limited precision when Kp is small, so do not use except to
   !!--++ obtain approximate V (eg for limits to eos)
   !!--++
   !!--++ Date: 06/03/2019
   !!
   Function Get_Volume_K(K, T, EoS) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: K       ! Bulk modulus
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoS     ! Eos Parameter
      real(kind=cp)              :: V
      !---- Local Variables ----!
      integer                            :: ic
      real(kind=cp)                      :: vprev,Kprev,kc,vnew,delv,delvprev,Vstep

      !> Init
      v=0.0_cp

      Vprev=EoS%params(1)             !This is Vo
      Kprev=K_cal(Vprev,T,EoS)        !This is Ko

      !> Initial search
      Vstep=0.001_cp*EoS%params(1)            !If K is smaller than Ko, must go up in volume
      if (K > Kprev) Vstep=-1.0_cp*Vstep

      do
         Vprev=Vprev+Vstep
         kc=K_cal(Vprev,T,EoS)
         if (K < Kprev)then                   !going to volumes > 1.0
            if (Kc < K)exit
            if (Vprev > 2.0_cp*EoS%params(1) )then
               V=2.0_cp*EoS%params(1)            !stop infinite looping
               return
            end if

         else
            if (Kc > K)exit
            if (Vprev < 0.001_cp*EoS%params(1) )then
               V=0.001_cp*EoS%params(1)      !stop infinite looping
               return
            end if
        end if
      end do

      !> set-up for Newton-raphson
      ic=0
      Kprev=Kc
      V=Vprev-Vstep

    !     write(unit=6,fmt='(a,f10.3,a,f10.3)') 'Got to NR, Kprev = ',Kprev,'      V = ',V
    !     write(unit=6,fmt='(a,f10.3,a,f10.3)') 'Got to NR, Vprev = ',Vprev,'  Vstep = ',Vstep

      do     ! does a newton-raphson search
         ic=ic+1
         if (ic > 10)exit
         err_eos=.false.
         !  write(unit=6,fmt='(a,f10.3,a,f10.3)') 'In NR, T =  ',T,'      V = ',V
         kc=K_cal(V,T,EoS)

         !   write(unit=6,fmt='(a,f10.3)') 'In  NR, Kc = ',Kc
         if (abs(kc-k) < 0.001*k) exit                                 !new limit May 2019: was 0.001 now 0.001*K
         if (ic > 1)then                                              !introduced if(ic > 1) May 2019: otherwise delVprev is not initialised and delV becomes nan
            delV= (k-kc)*(V-Vprev)/(kc-Kprev)
            !   write(unit=6,fmt='(a,f10.3)') 'In  NR, delv= ',delv
            if (abs(delV) > abs(delVprev))delV=sign(delVprev,delV)       !prevents step getting bigger

         else
            delV=Vstep
         end if

         Vnew= V + delV
         if (Vnew < 0._cp) Vnew=0.99*V          ! stops V going negative
         Kprev=Kc
         Vprev=V
         delVprev=delV
         V=vnew
      end do

      return
   End Function Get_Volume_K

   !!--++
   !!--++ FUNCTION GET_VOLUME_K_OLD
   !!--++
   !!--++
   !!
   Function Get_Volume_K_old(K,T,E) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: K       ! Bulk modulus
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: E      ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: V,vprev,Kprev,kc,vnew,delv,delvprev

      !> Init
      v=0.0_cp

      Vprev=e%params(1)             !This is Vo
      Kprev=K_cal(Vprev,T,E)        !This is Ko

      V=1.01_cp*e%params(1)

      do  ! does a newton-raphson search
         err_eos=.false.
         kc=K_cal(V,T,E)
         if (abs(kc-k) < 0.001)exit
         delV= (k-kc)*(V-Vprev)/(kc-Kprev)
         if (abs(delV) > abs(delVprev))delV=sign(delVprev,delV)       !prevents step getting bigger
         Vnew= V + delV
         if (Vnew < 0._cp)Vnew=0.99*V          ! stops V going negative
         Kprev=Kc
         Vprev=V
         delVprev=delV
         V=vnew
      end do

      !> Linear case
      if (e%linear) v=v**(1.0_cp/3.0_cp)

      return
   End Function Get_Volume_K_old

   !!----
   !!---- FUNCTION GET_VOLUME_S
   !!----
   !!---- Returns the value of Volume obtained from Strain (S)
   !!----
   !!---- Date: 15/02/2013
   !!
   Function Get_Volume_S(S, T, Eos) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S    ! Strain
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: V

      !---- Local Variables ----!
      real(kind=cp)                      :: v0
      real(kind=cp), dimension(N_EOSPAR) :: ev

      !> Init
      v=0.0_cp

      !> Local Eos Copy
      ev= EoS_to_Vec(EoS)

      !> Allow for thermal: s=function of V(P,T)/V(P=0,T)
      select case (EoS%itherm)
         case (0)
            v0=ev(1)

         case (1:n_therm_models)
            v0=get_volume(0.0_cp,t,EoS)
            if (EoS%linear) v0=v0**3.0_cp
      end select

      select case (EoS%imodel)
         case (1,5,7) ! Murnaghan and Tait, no strain defined
            v=0.0_cp

         case (2) ! Birch-Murnaghan
            V=v0*(1.0_cp+2.0_cp*s)**(-1.5_cp)

         case (3) ! Vinet
            V=v0*(1.0_cp-s)**3.0_cp

         case (4) ! Natural Strain
            V=v0*exp(-3.0_cp*s)
      end select

      !> Linear case
      if (EoS%linear) v=v**(1.0_cp/3.0_cp)

      return
   End Function Get_Volume_S

   !!--++
   !!--++ FUNCTION GET_VOLUME_THIRD
   !!--++
   !!--++ Returns the value of volume or length of a principal axis ieos in unit
   !!--++ cell in cell_eos at P,T when it can be calculated from others
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Get_Volume_Third(P, T, Cell_Eos, Ieos) Result(L)
      !---- Arguments ----!
      real(kind=cp),      intent(in) :: p,T
      type(eos_cell_type),intent(in) :: cell_eos
      integer,            intent(in) :: ieos     ! the axis (1,2,3) or V (0) to be calculated
      real(kind=cp)                  :: l        !returned length or volume

      !---- Local Variables ----!
      integer       :: i
      real(kind=cp) :: vfactor

      !> init
      l=10.0_cp

      !> safety check: if mono or triclinic, should only be called if angle poly used
      if (U_case(cell_eos%system(1:4)) == 'TRIC' .or. U_case(cell_eos%system(1:3)) == 'MONO')then
         if (cell_eos%eosang%iangle == 0)then
            err_eos=.true.
            err_eos_mess='Get_Volume_Third called for mono or triclinic, without angle poly set'
         end if
      end if

      !> Factor for unit cell volume V = a.b.c.vfactor
      vfactor=1.0
      if (U_case(cell_eos%system(1:4)) == 'TRIG' .or. U_case(cell_eos%system(1:3)) == 'HEX ') &
         vfactor=sqrt(3.0_cp)/2.0_cp

      call init_err_eos()

      select case(U_case(cell_eos%system(1:4)))
         case('ORTH','MONO','TRIC')
            vfactor=Get_Angle_Volfactor(P,T,cell_eos)

            select case(ieos)
               case(0)
                  l=Get_Volume(P,T,cell_eos%eos(1))*Get_Volume(P,T,cell_eos%eos(2))* &
                    Get_Volume(P,T,cell_eos%eos(3))*vfactor

               case default
                  l=Get_Volume(P,T,cell_eos%eos(0))/vfactor
                  do i=1,3
                     if (i == ieos)cycle
                     l=l/Get_Volume(P,T,cell_eos%eos(i))
                  end do
            end select

         case('TRIG','HEXA','TETR')
            select case(ieos)
               case(0)     ! calc volumes from a and c
                  l=Get_Volume(P,T,cell_eos%eos(1))**2.0_cp*Get_Volume(P,T,cell_eos%eos(3))*vfactor

               case(1)     ! a from V and c
                  l=sqrt(Get_Volume(P,T,cell_eos%eos(0))/vfactor/Get_Volume(P,T,cell_eos%eos(3)))

               case(3)     ! c from a and V
                  l=Get_Volume(P,T,cell_eos%eos(0))/vfactor/Get_Volume(P,T,cell_eos%eos(1))**2.0_cp
            end select

         case('CUBI','ISOT')
            select case(ieos)
               case(0)     ! calc volume from a
                  l=Get_Volume(P,T,cell_eos%eos(1))**3.0_cp

               case(1:3)     ! a,b, or c from V
                  l=Get_Volume(P,T,cell_eos%eos(0))**(1.0_cp/3.0_cp)
            end select
      end select

      return
   End Function Get_Volume_Third

   !!----
   !!---- FUNCTION ISOTROPIC_CELL
   !!----
   !!---- returns .true. if crystal system is isotropic or cubic
   !!----
   !!---- Date: 23/02/2021
   !!
   Function Isotropic_Cell(cell_eos) result(isotropic)

      !---- Arguments ----!
    type(eos_cell_type),intent(in)      :: cell_eos
    logical                             :: isotropic

    !---- Local Variables ----!
    character(len=len(cell_eos%system)) :: sys


    isotropic=.true.
    sys=U_case(cell_eos%system)
    if(index(sys,'ISOT') > 0)return
    if(index(sys,'CUB')  > 0)return
    isotropic=.false.

    return
    end function Isotropic_Cell

   !!----
   !!---- FUNCTION K_CAL
   !!----
   !!---- Returns value of K at this volume for EoS
   !!----
   !!--.. Changed code to use STRAIN function with v/v0: 27/02/2013
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- Date: 10/02/2017
   !!
   Function K_Cal(V, T, Eos, P) Result(Kc)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: V    ! Volume
      real(kind=cp),            intent(in) :: T    ! Temperature
      type(Eos_Type),           intent(in) :: EoS  ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: P    ! Pressure if known
      real(kind=cp)                        :: kc

      !---- Local Variables ----!
      integer                            :: j
      real(kind=cp)                      :: vv0,k0,kp,kpp,f,vol
      real(kind=cp)                      :: a,b,c,nu,vt,dPdV,delv
      real(kind=cp)                      :: vs,Ttr,dVs                         ! for transition calculations
      real(kind=cp), dimension(N_EOSPAR) :: ev
      real(kind=cp), dimension(3)        :: abc      ! Tait parameters
      real(kind=cp), dimension(-2:2)     :: pcal     ! array for calc p values numeric solutions
      real(Kind=cp)                      :: Ptrue    ! Input pressure if present. The true pressure without Pthermal subtracted
      real(Kind=cp)                      :: Pcorr    ! Pressure minus Pthermal.
      real(kind=cp),dimension(3,3)       :: apl      ! for APL parameters
      type(Eos_Type)                     :: EoST     ! Copy of eos


      !> Init
      kc=EoS%params(2)               ! safe default, K(P,T)=K0

      !> Pressure is needed for Tait, Murngahan, and for spontaneous strain at transitions
      !> This is the 'true' pressure, uncorrected for Pthermal
      if (present(P) )then
         Ptrue=P
      else
         Ptrue=get_pressure(v,t,EoS)
      end if

      !> PTV table
      if (EoS%imodel == -1) then
         if (.not. present(P) ) Ptrue=get_props_ptvtable(0.0,T,V,EoS,'P')
         kc=get_props_ptvtable(Ptrue,T,V,EoS,'K')
         return
      end if

      !> Isothermal EoS: get K0 and Kp0 at this T, then do algebra for normal EoS
      !>Or pthermal at Tref and do additional d(Pth)/dV contribution afterwards
      !>If transition, do calculation for bare phase, and add correction for transition afterwards

      !> Correct the volume to the high phase only if there is a transition
      if (EoS%itran > 0) then
         Vs=get_transition_strain(Ptrue,T,EoS)     ! returns the linear or volume strain as appropriate
         Vol=v/(1.0+vs)
      else
         Vol=v
      end if

      select case (EoS%itherm)
         case (0,6,7,8)                      ! No thermal model, or we have pthermal, so need params at Tref
            vv0=vol/EoS%params(1)           ! vv0 or aa0
            k0=EoS%params(2)
            kp=EoS%params(3)
            kpp=EoS%params(4)

         case (1:5)
            vv0=vol/get_V0_T(t,EoS)          ! vv0 is  v(p,t)/v(p=0,t): in transition, highphase
            k0=Get_K0_T(T,EoS)               ! returns M(T) for linear,
            if (err_eos) return                 ! exit with value eosparms(2) if k0 at P=0 calculated as negative

            kp=Get_Kp0_T(T,EoS)
            kpp=Get_Kpp0_T(T,EoS)
      end select

      !> Strain for BM, NS, Vinet EoS equations
      f=strain(vv0,EoS)
      if (err_eos) then
         err_eos_mess=trim(err_eos_mess)//' called from Kcal'
         return
      end if

      !> adjust for linear
      if (EoS%linear) then
         vol=vol**3.0_cp
         vv0=vv0**3.0_cp
         k0=k0/3.0_cp
         kp=kp/3.0_cp
         kpp=kpp/3.0_cp
      end if

      !> Now do the calculation for an isothermal eos at this T
      select case (EoS%imodel)
         case (1)   !Murnaghan
            Pcorr=Ptrue
            if (EoS%pthermaleos) Pcorr=Ptrue - pthermal(Vol,T,EoS)
            kc=k0+kp*Pcorr

         case (2)   !Birch-Murnaghan
            a=0.0_cp
            b=0.0_cp
            if (EoS%iorder > 2) a=1.5_cp*(kp-4.0_cp)
            if (EoS%iorder ==4) b =(9.0_cp*k0*kpp + (9.0_cp*kp*kp) - 63.0_cp*kp + 143.0_cp)/6.0_cp

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
            select case (EoS%iorder)
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
            kc= 3.0_cp*k0/vv0*( 1.0_cp/3.0_cp + a*f + b*f*f + c*f*f*f)

         case(5) ! Tait
            abc= get_tait(t,EoS)              ! get_tait returns volume-like parameters even for linear
            Pcorr=Ptrue
            if (EoS%pthermaleos) Pcorr=Ptrue - pthermal(Vol,T,EoS)
            if (abc(2)*Pcorr < -0.999999_cp) then
               err_eos=.true.
               err_eos_mess='Pressure yields zero bulk modulus for Tait Eos'  !Kc will be returned as eosparams(2)
               return
            else
               kc=k0*(1.0_cp-abc(1)*(1.0_cp-(1.0_cp + abc(2)*Pcorr)**(-1.0_cp*abc(3))))*(1.0_cp + abc(2)*Pcorr)** &
                  (1.0_cp+abc(3))
            end if

         case(6) ! APL
            apl= Get_APL(VV0,vol/vv0,K0,Kp,Kpp,EoS%params(5),EoS%iorder)
            kc=-1.0_cp*k0*vv0**0.333333_cp*(apl(2,1)*apl(1,2)*apl(1,3) + apl(1,1)*apl(2,2)*apl(1,3) + &
                apl(1,1)*apl(1,2)*apl(2,3))

         case(7) ! Kumar
            kc=k0*vv0*exp((kp+1.0_cp)*(1.0_cp-vv0))

         case default
            err_eos=.true.
            err_eos_mess='Invalid number for eos%imodel in K_cal'
            return
      end select

      !> To this point Kc is the bulk modulus of the 'bare' eos of the high phase if it was an isothermal eos
      !> Or for thermal-pressure EoS it is the bulk modulus at Tref and V

      !> Now correct thermal-pressure EoS for d(Pth)/dV contribution to bulk modulus
      if (EoS%Pthermaleos) then
         select case(EoS%itherm)
            case(7,8)           !MGD or Einstein EoS: Do this numerically,
               eost=EoS
               eost%itran=0    ! clear any transition terms
               delv=0.01_cp*vol
               do j=-2,2,1
                  vt=vol+real(j)*delv             ! apply shift to v
                  Pcal(j)=pthermal(vt,t,eost)     ! calc resulting thermal pressure
               end do
               dPdV=(Pcal(-2)+8.0_cp*(Pcal(1)-Pcal(-1))-Pcal(2))/(12.0_cp*delv)   !dPth/dV
               kc= kc -1.0_cp*vol*dPdV

            case default      !includes HP thermal-pressure which has dK = 0 along isochor

         end select
      end if

      !>Now correct from volume to linear if required, because transition code works in linear
      if (EoS%linear) kc=kc*3.0_cp

      !> Now handle phase transition
      if (EoS%itran > 0 ) then
         !> Local EoS copy
         ev= EoS_to_Vec(EoS)  ! but no scaling done

         select case(EoS%itran)
            case (1) ! Landau P-V
               if (abs(Ptrue-ev(21)) < 0.0001) then
                  !> at the transition pressure
                  if (abs(ev(25)-1.0_cp) < 0.0001) then
                     dVs=ev(24)*ev(25)           ! exact second order
                  else
                     dvs=huge(0._cp)             ! effectively infinite
                  end if
               else if (transition_phase(Ptrue,T,EoS) ) then             ! in the low phase
                  dVs=ev(24)*ev(25)*abs(Ptrue-ev(21))**(ev(25)-1.0_cp)      ! correct if lowP phase is highsym
                  if (nint(ev(20)) ==1) dVs=-1.0_cp*dVs                     ! change sign for highp phase is highsym
               else                                                         ! in the high phase
                  dVs=ev(26)*ev(27)*abs(Ptrue-ev(21))**(ev(27)-1)           ! correct if highP phase is highsym
                  if (nint(ev(20)) /= 1) dVs=-1.0_cp*dVs
               end if

            case (2,3) ! Landau VT or PVT
               if (transition_phase(Ptrue,T,EoS) ) then
                  !> low phase
                  Ttr = get_transition_temperature(Ptrue,EoS)
                  if (abs(Ttr-T) < 0.0001) then
                     dvs=0.0_cp
                  else
                     dVs=ev(24)*ev(25)*(abs(Ttr-T))**(ev(25)-1)*(ev(22)+2.0_cp*ev(23)*Ptrue) ! for curved phase boundary
                  end if
                  if (nint(ev(20)) /= 1) dVs=-1.0_cp*dVs                              ! sign change when low phase is highT
               else
                  !> in the highphase
                  Ttr = get_transition_temperature(Ptrue,EoS)
                  if (abs(Ttr-T) < 0.0001) then
                     dvs=0.0_cp
                  else
                     dVs=ev(26)*ev(27)*(abs(Ttr-T))**(ev(27)-1)*(ev(22)+2.0_cp*ev(23)*Ptrue)   !  curved boundary
                  end if
                  if (nint(ev(20)) ==1) dVs=-1.0_cp*dVs                                 ! sign change when high phase is highT
               end if
         end select

         !>Apply the correction to Kc due to dVs/dP:
         Vs=get_transition_strain(Ptrue,T,EoS) ! Vs  if volume eos, linear strain if linear eos
         kc=1.0_cp/(1.0_cp/kc - dVs/(1+Vs))
      end if

      return
   End Function K_Cal

   !!----
   !!---- FUNCTION KP_CAL
   !!----
   !!---- Returns value of kprime at this volume for EoS
   !!----
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- Date: 07/02/2017
   !!
   Function Kp_Cal(V, T, EoS, P) Result (Kpc)
      !---- Arguments ----!
      real(kind=cp),           intent(in) :: V       ! Volume
      real(kind=cp),           intent(in) :: T       ! Temperature
      type(Eos_Type),          intent(in) :: EoS  ! Eos Parameter
      real(kind=cp), optional, intent(in) :: P       ! Pressure if known
      real(kind=cp)                       :: kpc

      !---- Local Variables ----!
      integer                            :: j
      type(Eos_Type)                     :: eosbare
      real(kind=cp)                      :: vv0,k0,kp,kpp,ptr,vol,vs,ttr,betap_tran,betap_bare,kc,k
      real(kind=cp)                      :: a,b,f,rkp_top, rkp_bot,nu,dkdv,vt
      real(kind=cp), dimension(N_EOSPAR) :: ev
      real(kind=cp), dimension(3)        :: abc     ! Tait parameters
      real(kind=cp),dimension(3,3)       :: apl     ! for APL parameters
      real(kind=cp)                      :: dVsdP,d2VsdP2
      real(kind=cp),dimension(-2:2)      :: kpt
      real(kind=cp)                      :: delv,group,dgroup
      real(Kind=cp)                      :: Ptrue    ! Input pressure if present. The true pressure without Pthermal subtracted
      real(Kind=cp)                      :: Pcorr    ! Pressure minus Pthermal.

      !> Init
      kpc=0.0_cp

      !> Pressure is needed for Tait, Murngahan, and for spontaneous strain at transitions
      !> This is the 'true' pressure, uncorrected for Pthermal
      !> Using Pin avoids rounding errors building up when get_pressure is called, with transition eos.
      if (present(P) ) then
         Ptrue=P
      else
         Ptrue=get_pressure(v,t,EoS)
      end if

      !> PTV table
      if (EoS%imodel == -1) then
         if (.not. present(P) ) Ptrue=get_props_ptvtable(0.,T,V,EoS,'P')
         kpc=get_props_ptvtable(Ptrue,T,V,EoS,'KP')
         return
      end if

      !> Isothermal EoS: get K0 and Kp0 at this T, then do algebra for normal EoS
      !>Or pthermal at Tref and do additional d(Pth)/dV contribution afterwards
      !>If transition, do calculation for bare phase, and add correction for transition afterwards

      !> Correct the volume to the high phase only if there is a transition
      if (EoS%itran > 0) then
         Vs=get_transition_strain(Ptrue,T,EoS)     ! returns the linear or volume strain as appropriate
         Vol=v/(1.0+vs)
      else
         Vol=v
      end if

      select case (EoS%itherm)
         case (0,6,7,8)                           ! No thermal model, or we have pthermal, so need params at Tref
            vv0=vol/EoS%params(1)            ! vv0 or aa0
            k0=EoS%params(2)
            kp=EoS%params(3)
            kpp=EoS%params(4)

         case (1:5)
            vv0=vol/get_V0_T(t,EoS)          ! vv0 is  v(p,t)/v(p=0,t): in transition, highphase
            k0=Get_K0_T(T,EoS)               ! returns M(T) for linear,
            if (err_eos) return                 ! exit with value eosparms(2) if k0 at P=0 calculated as negative

            kp=Get_Kp0_T(T,EoS)
            kpp=Get_Kpp0_T(T,EoS)
      end select

      !> Strain for BM, NS, Vinet EoS equations
      f=strain(vv0,EoS)
      if (err_eos) then
         err_eos_mess=trim(err_eos_mess)//' called from Kpcal'
         return
      end if

      !>adjust for linear
      if (EoS%linear) then
         vol=vol**3.0_cp
         vv0=vv0**3.0_cp
         k0=k0/3.0_cp
         kp=kp/3.0_cp
         kpp=kpp/3.0_cp
      end if

      !>Now do the calculation for an isothermal eos at this T
      select case(EoS%imodel)
         case (1) ! Murnaghan
            kpc=kp

         case (2) ! Birch-Murnaghan
            a=0.0_cp
            b=0.0_cp
            if (EoS%iorder > 2) a=1.5_cp*(kp-4.0_cp)
            if (EoS%iorder ==4) b = (9.0_cp*k0*kpp + (9.0_cp*kp*kp) - 63.0_cp*kp + 143.0_cp)/6.0_cp

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
            if (EoS%iorder > 2) a=1.5_cp*(kp-2.0_cp)
            if (EoS%iorder ==4) b=1.5_cp*(1.0_cp + k0*kpp + (kp-2.0_cp) + (kp-2.0_cp)**2.0_cp)
            rkp_top=  1.0_cp + 2.0_cp/3.0_cp*a + 2.0_cp*(a+b)*f + 3.0_cp*b*f*f
            rkp_bot=  1.0_cp + (3.0_cp +2.0_cp*a)*f + 3.0_cp*(a+b)*f*f +3.0_cp*b*f*f*f
            kpc=1.0_cp+rkp_top/rkp_bot

          case(5) ! Tait
            abc= get_tait(t,EoS)              ! get_tait returns volume-like parameters even for linear
            Pcorr=Ptrue
            if (EoS%pthermaleos) pcorr=ptrue - pthermal(Vol,T,EoS)   ! correct P for pthermal, needed only for Tait
            if (abc(2)*Pcorr < -0.999999_cp) then
               err_eos=.true.
               err_eos_mess='Pressure yields infinite Kp for Tait Eos'
               kpc=9999.0        ! safe value return
            else
               kpc=(kp+1.0_cp)*((1.0_cp + abc(2)*Pcorr)**abc(3)*(1.0_cp-abc(1)) + abc(1)) -1.0_cp
            end if

         case(6) ! APL
            apl= Get_APL(VV0,vol/vv0,K0,Kp,Kpp,EoS%params(5),EoS%iorder)
            group=apl(2,1)*apl(1,2)*apl(1,3) + apl(1,1)*apl(2,2)*apl(1,3) + apl(1,1)*apl(1,2)*apl(2,3)

            dgroup=apl(3,1)*apl(1,2)*apl(1,3) + apl(1,1)*apl(3,2)*apl(1,3) + apl(1,1)*apl(1,2)*apl(3,3) &
                  +2.0_cp*(apl(1,1)*apl(2,2)*apl(2,3) + apl(2,1)*apl(1,2)*apl(2,3) + apl(2,1)*apl(2,2)*apl(1,3))
            kpc= -1.0_cp/3.0_cp - vv0**0.333333_cp*dgroup/3.0_cp/group

         case(7) ! Kumar
            kpc=(kp+1.0_cp)*vv0 -1

         case default
            err_eos=.true.
            err_eos_mess='Invalid number for eos%imodel in K_cal'
            return

      end select

      !> To this point Kpc is the bulk modulus of the 'bare' eos of the high phase if it was an isothermal eos
      !> Or for thermal-pressure EoS it is the bulk modulus at Tref and V

      !> Now correct thermal-pressure EoS for d(Pth)/dV contribution to bulk modulus, when possible
      if (EoS%Pthermaleos) then
         select case(EoS%itherm)
            case(7,8)           !MGD EoS: Do this numerically on complete EoS without transition
               eosbare=EoS
               eosbare%itran=0    ! clear any transition terms
               delv=0.01_cp*vol
               do j=-2,2,1             ! loop calculates dK/dV
                  vt=vol+real(j)*delv  ! apply shift to v
                  kpt(j)=k_cal(vt,t,eosbare)     ! calc resulting k
               end do
               dkdV=(Kpt(-2)+8.0_cp*(Kpt(1)-Kpt(-1))-Kpt(2))/(12.0_cp*delv)
               kpc= -1.0_cp*dkdV*vol/kpt(0)

            case default      !includes HP thermal-pressure which has dK = 0 along isochor


         end select
      end if

      !> Now correct from volume to linear if required, because transition code works in linear
      if (EoS%linear) kpc=kpc*3.0_cp

      !> And the transition effects
      if (EoS%itran > 0 ) then
         !> Local EoS copy, but no scaling of transition parameters for linear
         ev= EoS_to_Vec(EoS)

         select case(EoS%itran)
            case(1) ! Landau P-V
               if (transition_phase(Ptrue,T,EoS)) then             ! in the low phase
                  dVsdP=ev(24)*ev(25)*abs(Ptrue-ev(21))**(ev(25)-1)   ! d(Vs)/dP correct if lowP phase is highsym
                  if (nint(ev(20)) /=1) dVsdP=-1.0_cp*dVsdP           ! change sign for lowp phase is highsym
                  d2VsdP2=dVsdP*(ev(25)-1)/abs(Ptrue-ev(21))
                  if (nint(ev(20)) /=1) d2VsdP2=-1.0_cp*d2VsdP2
               else                                                                      ! in the high phase
                  dVsdP=ev(26)*ev(27)*abs(Ptrue-ev(21))**(ev(27)-1)     ! correct if highP phase is highsym
                  if (nint(ev(20)) == 1) dVsdP=-1.0_cp*dVsdP
                  d2VsdP2=dVsdP*(ev(27)-1)/abs(Ptrue-ev(21))
                  if (nint(ev(20)) ==1) d2VsdP2=-1.0_cp*d2VsdP2
               end if

            case(2,3) ! Landau VT or PVT
               Ttr = get_transition_temperature(Ptrue,EoS)
               Ptr = get_transition_pressure(T,EoS)
               if (transition_phase(Ptrue,T,EoS)) then
                  if (abs(Ttr-T) < 0.0001) then
                     dVsdP  = 0._cp            ! on the boundary
                     d2VsdP2= 0._cp
                  else                          ! in the low phase
                     dVsdP=ev(24)*ev(25)*(abs(Ttr-T))**(ev(25)-1)*(ev(22)+2.0_cp*ev(23)*Ptrue) ! for curved phase boundary
                     if (nint(ev(20)) /= 1) dVsdP=-1.0_cp*dVsdP                              ! sign change when low phase is highT

                     d2VsdP2= dVsdP*(2.0_cp*ev(23))/(ev(22)+2.0_cp*Ptrue*ev(23))      ! This works 15-Feb
                     if (nint(ev(20)) == 1) then
                        d2VsdP2= d2VsdP2 + dVsdP*(ev(25)-1)/abs(Ttr-T)*(ev(22)+2.0_cp*Ptrue*ev(23))
                     else
                        d2VsdP2= d2VsdP2 - dVsdP*(ev(25)-1)/abs(Ttr-T)*(ev(22)+2.0_cp*Ptrue*ev(23))
                     end if
                  end if

               else
                  !> in the highphase
                  if (abs(Ttr-T) < 0.0001) then
                     dVsdP  = 0._cp
                     d2VsdP2= 0._cp
                  else
                     dVsdP=ev(26)*ev(27)*(abs(Ttr-T))**(ev(27)-1)*(ev(22)+2.0_cp*ev(23)*Ptrue)   !  curved boundary
                     if (nint(ev(20)) ==1) dVsdP=-1.0_cp*dVsdP                                 ! sign change when high phase is highT
                     d2VsdP2=dVsdP*(2.0_cp*ev(23))/(ev(22)+2.0_cp*Ptrue*ev(23))
                     if (nint(ev(20)) == 1) then
                        d2VsdP2=d2VsdP2-1.0_cp*dVsdP*(ev(27)-1)/abs(Ttr-T)*(ev(22)+2.0_cp*Ptrue*ev(23))  !high phase is highT
                     else
                        d2VsdP2=d2VsdP2+1.0_cp*dVsdP*(ev(27)-1)/abs(Ttr-T)*(ev(22)+2.0_cp*Ptrue*ev(23)) !high phase is lowT
                     end if
                  end if
               end if
         end select

         !> Now apply the transition effect
         Vs=get_transition_strain(Ptrue,T,EoS) ! Vs  if volume eos, linear strain if linear eos

         betap_tran=   d2VsdP2/(1.0_cp+Vs) -(dVsdP/(1.0_cp+Vs))**2.0_cp         ! transition contribution

         eosbare=EoS
         eosbare%itran=0     ! create a bare phase eos
         kc=k_cal(vol,t,eosbare,p=Ptrue)   ! k_cal returns linear moduli if linear
         k=k_cal(v,t,EoS,p=Ptrue)
         betap_bare= -1.0_cp*kpc/kc**2.0_cp
         kpc= (betap_tran - betap_bare)*(k**2.0_cp)
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
   !!---- Date: 25/012/2016
   !!
   Function Kpp_Cal(V, T, EoS) Result (Kppc)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: kppc

      !---- Local Variables ----!
      integer                      :: j
      real(kind=cp)                :: ptr,vlimit,plimit,vlimitk
      real(kind=cp)                :: delp,p0,p,vt
      real(kind=cp),dimension(-2:2):: kpt

      !> Init
      kppc=0.0_cp

      !> Trap for pvt table....do calculation later
      if(EoS%imodel == -1)return

      !> jump on equation type
      if (EoS%imodel == 1 .and. EoS%itran == 0) return ! Murnaghan without transition

      !> numerical solution
      ! To get reasonable Kprime numerically, need delP as 1 GPa for K0=160 GPa because this generates 1% change in Kp
      ! Therefore need as big del for Kpp, if not bigger

      delp=K_cal(V,T,EoS)/100.        !delp=K/100
      p0=get_pressure(V,T,EoS)         ! current p

      !> Code to prevent stepping across transition
      if (EoS%itran > 0) then
         Ptr=get_transition_pressure(t,EoS)
         if (transition_phase(P0+2.0*delp,T,EoS) .neqv. transition_phase(P0,T,EoS)) delp=0.2*abs(ptr-p0)
         if (transition_phase(P0-2.0*delp,T,EoS) .neqv. transition_phase(P0,T,EoS)) delp=0.2*abs(ptr-p0)
      end if

      !> Code to stop MGD EoS going into illegal large volume at negative delp
      if (EoS%itherm == 7 .or. EoS%itherm == 8)then
         vlimitk=get_volume_K(0._cp,t,EoS)
         vlimit=get_volume_K(EoS%params(2)/2.0_cp,EoS%tref,EoS)
         if (vlimitk > tiny(0.0_cp) .and. vlimitk < vlimit)vlimit=vlimitk
         plimit=get_pressure(vlimit,t,EoS)
         if (p0-2.0*delp < plimit)delp=0.2*abs(p0-plimit)
      end if

      do j=-2,2,1
         p=p0+real(j)*delp                 ! apply shift to p
         vt=get_volume(p,t,EoS)
         kpt(j)=kp_cal(vt,t,EoS,p=p)     ! calc resulting Kp
      end do

      kppc=(kpt(-2)+8.0_cp*(kpt(1)-kpt(-1))-kpt(2))/(12.0_cp*delp)     ! Derivative to second order approximation

      !> No linear conversion is required because kp_cal returns values for "linear Kp" = Mp,
      !> so kppc is already dMp/dP = Mpp

      return
   End Function Kpp_Cal

   !!--++
   !!--++ FUNCTION Linear_EoS_Allowed
   !!--++
   !!--++ Date: 03/02/2021
   !!
   Function Linear_EoS_Allowed_Eos(Eos) Result(Allowed)
      !---- Arguments ----!
      type(Eos_Type), intent(in) :: EoS    ! Eos Parameter
      logical                    :: allowed

      !---- Local Variables ----!

      !> init
      allowed=.true.

      !> Model
      if (eos%imodel == 6) then
         allowed=.false.    !APL
         return
      end if

      !> Thermal model
      select case (eos%itherm)
         case (6:8)
            allowed=.false.
      end select

      return
   End Function Linear_EoS_Allowed_Eos

   !!--++
   !!--++ FUNCTION Linear_EoS_Allowed
   !!--++
   !!--++ Date: 03/02/2021
   !!
   Function Linear_EoS_Allowed_I(Imodel, Itherm) Result(Allowed)
      !---- Arguments ----!
      integer, optional, intent(in) :: imodel      !number of eos PV model
      integer, optional, intent(in) :: itherm      !number of thermal model
      logical                       :: allowed

      !---- Local Variables ----!

      !> init
      allowed=.true.

      !> Model
      if (present(imodel)) then
         if (imodel == 6) allowed=.false.    !APL
      end if

      !> Thermal model
      if (present(itherm)) then
         if (itherm == 6 .or. itherm == 7 .or. itherm == 8) allowed=.false.  !MGD & Einstein
      end if

      return
   End Function Linear_EoS_Allowed_I

   !!--++
   !!--++ FUNCTION MURN_INTERPOLATE_PTVTABLE
   !!--++
   !!--++ Returns ...
   !!--++
   !!--++ Date: 25/012/2016 RJA
   !!
   Function Murn_Interpolate_PTV_Table(I, J, P, Eos) Result(V)
      !---- Arguments ----!
      integer,        intent(in) :: I      ! pointers to table point just beyond P,T required
      integer,        intent(in) :: J      ! pointers to table point just beyond P,T required
      real(kind=cp),  intent(in) :: P      ! pressure required
      type(Eos_Type), intent(in) :: EoS    ! Eos Parameter

      !---- Local Variables ----!
      integer         :: is
      real(kind=cp)   :: k12,K23,KP,V0,K0,V

      !> Set up pointer
      is=i-2
      if (is < 1) is=1

      !> Estimate K and Kp
      k12=-0.5_cp*(eos%table%ptv(is+1,j,3) + eos%table%ptv(is,j,3))*(eos%table%ptv(is+1,j,1) - &
           eos%table%ptv(is,j,1))/(eos%table%ptv(is+1,j,3)-eos%table%ptv(is,j,3))
      k23=-0.5_cp*(eos%table%ptv(is+2,j,3) + eos%table%ptv(is+1,j,3))*(eos%table%ptv(is+2,j,1) - &
           eos%table%ptv(is+1,j,1))/(eos%table%ptv(is+2,j,3)-eos%table%ptv(is+1,j,3))

      kp=2.0_cp*(k23-k12)/(eos%table%ptv(is+2,j,1)-eos%table%ptv(is,j,1))       ! Kp at point is+1
      if (abs(kp) < 0.0001) kp=4.0_cp
      k0=k12+kp*(eos%table%ptv(is+1,j,1)-eos%table%ptv(is+1,j,1))/2.0_cp          ! K at point is+1
      v0=eos%table%ptv(is+1,j,3)                                                  ! V0 at point is+1

      !> Value
      v=v0*(1.0_cp+kp/k0*(p-eos%table%ptv(is+1,j,1)))**(-1.0_cp/kp)

      return
   End Function Murn_Interpolate_PTV_Table

   !!--++
   !!--++ SUBROUTINE MURN_PTV_TABLE
   !!--++
   !!--++ Calculate ....
   !!--++
   !!--++ Date: 17/03/2017
   !!
   Function Murn_PTV_Table(I, J, Eos) Result(Eosm)
      !---- Arguments ----!
      integer,        intent(in)  :: i      ! pointers to table point
      integer,        intent(in)  :: j      ! pointers to table point
      type(Eos_Type), intent(in)  :: EoS    ! Eos Parameter in (with ptvtable)
      type(Eos_Type)              :: EoSm   ! Eos Parameter of Murn output for point i,j

      !---- Local Variables ----!
      integer         :: im                   ! normally i but if edge point it will be shifted
      real(kind=cp)   :: k2,k4

      !> Set up Murn in Eosm
      call Init_EoS_Type(Eosm)

      eosm%imodel=1
      call set_eos_names(eosm)      ! sets names for params
      eosm%title='Murnaghan for interpolation'

      !> Set up pointer
      im=i
      if (i < 3) im=3
      if (i > eos%table%np-2) im=eos%table%np-2

      !>Vo
      eosm%params(1)=eos%table%ptv(i,j,3)

      !>Estimate K0 at point i,j
      eosm%params(2)=-1.0_cp*eos%table%ptv(im,j,3)*(eos%table%ptv(im+1,j,1)-eos%table%ptv(im-1,j,1))/ &
                     (eos%table%ptv(im+1,j,3)-eos%table%ptv(im-1,j,3))

      !> Kp
      k2=-1.0_cp*eos%table%ptv(im-1,j,3)*(eos%table%ptv(im,j,1)-eos%table%ptv(im-2,j,1))/ &
         (eos%table%ptv(im,j,3)-eos%table%ptv(im-2,j,3))

      k4=-1.0_cp*eos%table%ptv(im+1,j,3)*(eos%table%ptv(im+2,j,1)-eos%table%ptv(im,j,1))/ &
         (eos%table%ptv(im+2,j,3)-eos%table%ptv(im,j,3))

      eosm%params(3)=(k4-k2)/(eos%table%ptv(im+1,j,1)-eos%table%ptv(im-1,j,1))       ! Kp at point im

      if (abs(eosm%params(3)) < 0.0001) eosm%params(3)=4.0_cp

      !> adjust K now if im /= i; Kp should be constant for a Murn
      if (i /= im) eosm%params(2)=eosm%params(2)+eosm%params(3)*(eos%table%ptv(i,j,1)-eos%table%ptv(im,j,1))

      return
   End Function Murn_PTV_Table

   !!--++
   !!--++ FUNCTION NORMPRESSURE_EOS
   !!--++
   !!--++ Returns the value of Normalised Pressure (F) at this Strain (S) using
   !!--++ the EoS parameters as K0, V0, Kp, Kpp
   !!--++
   !!--++ Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!--++ Modified for thermal EoS: requires T to be meaningful
   !!--++
   !---++ Date: 10/09/2013
   !!
   Function NormPressure_Eos(S, T, Eos) Result(F)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S    ! Strain
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: F

      !---- Local Variables ----!
      real(kind=cp)                     :: k0,kp,kpp, b,c,v0
      real(kind=cp),dimension(n_eospar) :: ev

      !> Init
      f=0.0_cp

      !> Local copy
      ev= EoS_to_Vec(EoS)

      !> Get correct parameters for this T
      !> 17/11/2014: Before transitions, this used eosparams for normal thermal expansion
      !> but with transitions safer to do following:

      if (EoS%itherm == 0 .and. EoS%itran == 0) then
         ! simple PV eos
         k0=ev(2)
         kp=ev(3)
         kpp=ev(4)

      else if(EoS%itran == 0 .and. .not. EoS%pthermaleos) then
         ! normal thermal expansion models with dK/dT and no transition
         k0=Get_K0_T(T,EoS)                     ! returns M(T) for linear,
         if (EoS%linear) k0=k0/3.0_cp
         kp=ev(3)
         kpp=ev(4)

      else
         ! Transition model, or Pthermal: cannot use get_V0_T or get_K0_T because they return V and K at Tref for pthermal
         v0=get_volume(0.0_cp,T,EoS)        ! determine V0 at this T
         k0=k_cal(v0,T,EoS)
         kp=kp_cal(v0,T,EoS)
         kpp=kpp_cal(v0,T,EoS)
         if (EoS%linear) then
            k0=k0/3.0_cp
            kp=kp/3.0_cp
            kpp=kpp/3.0_cp
         end if
      end if

      select case(EoS%imodel)
         case (1,5,7) ! Murnaghan, Tait, Kumar
            f=0.0_cp

         case (2) ! Birch-Murnaghan
            b=0.0_cp
            c=0.0_cp
            if (EoS%iorder > 2) b=1.5_cp*(kp-4.0_cp)
            if (EoS%iorder ==4) c=1.5_cp*(k0*kpp + (kp-4.0_cp)*(kp-3.0_cp)+35.0_cp/9.0_cp)
            f=K0*(1.0_cp + b*s + c*s*s)

         case (3) ! Vinet: new definition RJA 28/03/2013
            f= K0*exp(1.5_cp*(Kp-1.0_cp)*s)

         case (4) ! Natural Strain
            b=0.0_cp
            c=0.0_cp
            if (EoS%iorder > 2) b=1.5_cp*(Kp - 2.0_cp)
            if (EoS%iorder ==4) c=1.5_cp*(K0*Kpp + 1.0_cp + (Kp-2.0_cp)+(Kp-2.0_cp)**2.0_cp)
            f=K0*(1.0_cp + b*s + c*s*s)
      end select

      return
   End Function NormPressure_Eos

   !!--++
   !!--++ FUNCTION NORMPRESSURE_P
   !!--++
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
   Function NormPressure_P(S, P, Imodel) Result(F)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S      ! Strain
      real(kind=cp),  intent(in) :: P      ! Presure
      integer      ,  intent(in) :: imodel ! type of eos
      real(kind=cp)              :: f

      !---- Local Variables ----!

      !> Init: a zero value is returned for undefined values
      f=0.0_cp

      !> If the finite strain is zero, the normalised pressure is not defined
      if (s > tiny(0.0) ) then
         select case (imodel)
            case (1,5,6,7) ! Murnaghan, Tait, APL, Kumar
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
   Function Pressure_F(F, S, Eos) Result(P)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: F    ! Normalized Pressure
      real(kind=cp),  intent(in) :: S    ! Strain
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: p

      !---- Local Variables ----!
      real(kind=cp) :: cF

      !> Init
      p=0.0_cp
      cf=f

      select case (eos%imodel)
         case (1,5,6,7) ! Murnaghan, Tait, APL, Kumar

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
   !!--++ FUNCTION PRINCIPAL_EOS
   !!--++
   !!--++  For use with eos_cell_type: Converts request for recip axis to real axis, if allowed
   !!--++
   !!--++ Date: 09/09/2020
   !!
   Function Principal_EoS(Cell_Eos, I) Result(Ieos)
      !---- Arguments ----!
      type(eos_cell_type),intent(in) :: cell_eos !eos for cells
      integer,            intent(in) :: i   !proposed axis number
      integer                        :: ieos

      !---- Local Variables ----!

      select case(i)
         case(:-3,-1,7:)
            ieos=-1                   !error

         case(-2,0:3)                    ! general dir, V a b c
            ieos=i

         case(4:6)
            select case(U_case(cell_eos%system(1:4)))
               case('ISOT','CUBI','TETR','ORTH')
                  ieos=i-3

               case('HEXA','TRIG')       !onl d(001) is equivalent to an axis in length
                  if (i == 6)then
                     ieos=3
                  else
                     ieos=-2
                  end if

               case('MONO')
                  if (i-3 == cell_eos%unique)then
                     ieos=cell_eos%unique
                  else
                     ieos=-2
                  end if

               case('TRIC')
                  ieos=i
            end select
      end select

      return
   End Function Principal_EoS

   !!----
   !!---- FUNCTION PSCALEMGD
   !!----
   !!---- Date: 31/03/2021
   !!
   Function PscaleMGD(EoS) Result(MGD)
      !---- Arguments ----!
      type(Eos_Type), intent(in)  :: EoS          ! EoS
      logical                     :: MGD        ! .true. if e%pscale_name iskbar or Gpa

      !---- Local Variables ----!
      character(len=len(eos%vscale_name)) :: vname

      MGD=.false.

      vname=adjustl(U_case(eos%pscale_name))
      if (len_trim(vname) == 0)return

      if (index(vname,'GPA') > 0 ) MGD=.true.
      if (index(vname,'KBAR') > 0 ) MGD=.true.


      return
   End Function PscaleMGD

   !!----
   !!---- FUNCTION PTHERMAL
   !!----
   !!----  Calculate Pthermal from eosparameters at temperature T
   !!----
   !!---- Date: 10/09/2013
   !!
   Function Pthermal(V, T, Eos, J) Result(Pth)
      !---- Arguments ----!
      real(kind=cp),   intent(in) :: V       ! Volume: not needed for HP2011 of linear-thermalpressure, needed for MGD, q-comp and Einstein
                                            !assumed this is  'a' if linear
      real(kind=cp),   intent(in) :: T       ! Temperature
      type(Eos_Type),  intent(in) :: EoS     ! Eos Parameter
      integer,optional,intent(in) :: j      ! which oscillator: then Pthermal only calculates for this one
      real(kind=cp)               :: pth
      !---- Local Variables ----!
      integer                           :: i,jo
      real(kind=cp)                     :: thtref,exp0,eta0,eth,eth0
      real(kind=cp)                     :: gammaV, thetaD,thetaE
      real(kind=cp),dimension(0:2)      :: pthp  !contributions to pth
      real(kind=cp),dimension(N_EOSPAR) :: ev


      !> Local copies
      ev= eos_to_vec(EoS)    !handle linear case
      jo=-1     !calculate all
      if(present(j))then
          if(j > -1 .and. j < 3)jo=j
      endif
      pthp=0._cp
      pth=0._cp

      select case (EoS%itherm)
         case (6) ! Thermal pressure from Holland and Powell 2011
            thtref=ev(11)/EoS%tref         ! T_einstein/Tref
            exp0=exp(thtref)                  ! exp(T_Ein/Tref)
            eta0= thtref*thtref*exp0/(exp0-1.0_cp)**2.0_cp  ! eta at Tref

            pthp(0) = ev(10)*ev(2)*ev(11)/eta0*( 1.0_cp/(exp(ev(11)/t)-1.0_cp) -1.0_cp/(exp0 -1.0_cp))
            !no scaling required because pth is scaled by Ko

         case (7) ! MGD in the form of Kroll et al (2012)
            thetaD=get_DebyeT(V,EoS)                 ! Both get_Debye and get_grun expect 'a' value if linear
            gammaV=get_grun_V(V,EoS)
            pthp(0)=gammaV/v*(EthDebye(T,thetaD,EoS%params(13))-EthDebye(EoS%tref,thetaD,EoS%params(13)))

         case(8)   !Einstein oscillator
            gammaV=get_grun_V(V,EoS)
            thetaE=get_DebyeT(V,EoS)
            eth=EthEinstein(T,thetaE,EoS%params(13))
            eth0=EthEinstein(EoS%tref,thetaE,EoS%params(13))
            pthp(0)=gammaV/v*(eth-eth0)

         case default
            pthp(0)=0.0_cp
      end select

      !>Extra oscillators: only allowed in combination with models 6,7,8
      if (EoS%osc_allowed .and. sum(EoS%iosc) > 0)then
         pthp(0)=(1._cp-EoS%params(40)-EoS%params(45))*pthp(0)     ! partial contribution main oscillator

         do i=1,2
            select case(EoS%iosc(i))
               case(0)
                  cycle

               case(1)  !DEBYE
                  thetaD=get_DebyeT(V,EoS,i)
                  gammaV=get_grun_V(V,EoS,i)
                  pthp(i)=gammaV/v*EoS%params(35+5*i)*  &
                      (EthDebye(T,thetaD,EoS%params(13))-EthDebye(EoS%tref,thetaD,EoS%params(13)))

               case(2)     ! Einstein
                  gammaV=get_grun_V(V,EoS,i)
                  thetaE=get_DebyeT(V,EoS,i)
                  pthp(i)=gammaV/v*EoS%params(35+5*i)*  &
                      (EthEinstein(T,thetaE,EoS%params(13))-EthEinstein(EoS%Tref,thetaE,EoS%params(13)))
            end select
         end do
      end if

      ! if the thermal energy was from EthDebye or EthEinstein, it is in J/mol pth
      ! Then if V in m3/mol  Eth/V is in J/m3=Pa
      select case(EoS%itherm)
         case(7,8)
            pthp=pthp*EPthermal_factor(EoS)
      end select

      ! Now return requested part of pth:
      if (jo == -1)then
         pth=sum(pthp)
      else
         pth=pthp(jo)
      end if

      return
   End Function Pthermal

   !!----
   !!---- FUNCTION SET_XDATATYPES
   !!----
   !!----  returns array result xdatatypes(i)=1 if datatype i is present in dataset gdat
   !!----  if used=.true. then only = 1 if at least one datum of the type is used
   !!----
   !!---- Date: 18/02/2021
   !!
   Function Set_XdataTypes(Gdat, Used) Result(XdataTypes)
      !---- Arguments ----!
      type (EoS_Data_List_Type), intent(in) :: Gdat  ! the data list
      logical,                   intent(in) :: Used  ! .true. to request
      integer, dimension(0:N_DATA_TYPES)    :: XdataTypes

      !---- Local Variables ----!
      integer :: i

      !> Init
      xdatatypes=0

      if (used)then
         do i=1,gdat%n
            if (gdat%eosd(i)%iuse == 1)xdatatypes(gdat%eosd(i)%xtype) = 1
         end do

      else        ! all data, used or not
         do i=1,gdat%n
            xdatatypes(gdat%eosd(i)%xtype) = 1
         end do
      end if

      return
   End Function Set_XdataTypes

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
   Function Strain(VV0, Eos) Result(S)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: VV0  ! Volume divided by V0 at this temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      real(kind=cp)              :: s

      !---- Local Variables ----!
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
      if (eos%linear) cvv0=vv0**3.0_cp

      select case (eos%imodel)
         case (1,5,6,7) ! Murnaghan, Tait, APL, Kumar
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
   Function Strain_EoS(V, T,Eos) Result(S)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V    ! Volume
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameters
      real(kind=cp)              :: s

      !---- Local Variables ----!
      real(kind=cp) :: vvo

      !> Init
      s=0.0_cp
      if (v <= 0.00001) then
         Err_EoS=.true.
         Err_EoS_Mess="Strain calculation called with invalid V =< 0: strain set zero"
         return
      end if

      !> Calculate
      vvo=v/get_volume(0.0_cp,t,eos)     ! vv0 is V(P,T)/V(0,T) or a(P,T)/a(0,T)
      s=strain(vvo,eos)                  ! cubes vv0 on input if linear

      return
   End Function Strain_EOS

   !!----
   !!---- FUNCTION THERMAL_PRESSURE_EOS
   !!----
   !!----
   !!---- Date: 03/02/2021
   !!
   Function Thermal_Pressure_Eos(I) Result(Pth)
      !---- Arguments ----!
      integer, intent(in) :: i        !a thermal model number
      logical             :: pth

      !---- Local Variables ----!

      select case(i)
         case(1:5)
            pth=.false.

         case(6:10)
            pth=.true.

         case default
            pth=.false.
      end select

      return
   End Function Thermal_Pressure_Eos

   !!--++
   !!--++ FUNCTION TRANSFORM_ESD
   !!--++
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
   Function Transform_Esd(P, T, Eos) Result(Esd)
      !---- Arguments ----!
      real(kind=cp),   intent(in)              :: p        ! Pressure at which to calculate parameter esds
      real(kind=cp),   intent(in)              :: t        ! Temperature at which to calculate parameter esds
      type (EoS_Type), intent(in)              :: eos      ! The EoS parameters and vcv
      real(kind=cp),dimension(N_EOSPAR)        :: esd      ! The esd's of Eos parameters at this P and T

      !---- Local Variables ----!
      real(kind=cp),dimension(N_EOSPAR,N_EOSPAR)  :: d                  ! cross derivatives of parameters
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
         if (EoS%iuse(j) == 1) then
            do k=-2,2,1
               eost=EoS                                                     !reset eos params
               if (abs(EoS%params(j)) > tiny(0.0_cp)) then
                  shift=factor*EoS%params(j)
                  if (j == 10) shift=10.0*shift                    ! alpha0
               else
                  shift=1.0_cp/EoS%factor(j)                    ! shift to a parameter with zero value
                  if (j == 5) shift=0.002                          ! dK/dT has print factor 1, but typical value -.02
               end if

               eost%params(j)=EoS%params(j)+float(k)*shift    ! apply shift to a parameter
               Par(k,1:6)= EoS_Cal(P,T,Eost)                 ! calc resulting parvals
            end do
            d(1:6,j)=(par(-2,1:6)+8.0_cp*(par(1,1:6)-par(-1,1:6))-par(2,1:6))/(12.0_cp*shift) ! derivative to second order approximation
         end if
      end do

      !> d(1:6,j) contains the derivatives in order V,K0,Kp,Kpp,dK/dT,alpha0 wrt param(j) at Pref,Tref
      !  now switch these to 10 for alpha0 and ignore any other alpha coeffs
      d(10,1:10)=d(6,1:10)
      d(6,1:10)=0.0_cp

      !> Now calculate esd-squared array = vcv(i,i). The input vcv is already linear for linear eos!
      do k=1,N_EOSPAR
         do i=1,N_EOSPAR
            if (EoS%iref(i) == 0) cycle ! because vcv(i,j) will be 0
            do j=1,N_EOSPar
               if (EoS%iref(j) == 0) cycle ! because vcv(i,j) will be 0
               esd(k)=esd(k)+EoS%vcv(i,j)*d(k,i)*d(k,j)
            end do
         end do
      end do

      !> Final: extra trap June 2016 in case round off leaves esd(i)^2 < 0
      do i=1,N_EOSPAR
          if (esd(i) > tiny(0.0_cp) ) then
             esd(i)=sqrt(esd(i))
          else
             esd(i)=0.0_cp
          end if
      end do

      return
   End Function Transform_Esd

   !!----
   !!---- LOGICAL FUNCTION TRANSITION_PHASE
   !!----
   !!---- Returns .true. if P and T are in the low phase stability field
   !!---- and .false. (default) if in the high-symm field, or exactly on
   !!---- the boundary.
   !!----
   !!---- Date: 17/07/2015
   !!
   Function Transition_Phase(P, T, Eos) Result(Ip)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P    ! Pressure
      real(kind=cp),  intent(in) :: T    ! Temperature
      type(Eos_Type), intent(in) :: EoS  ! Eos Parameter
      logical                    :: Ip

      !---- Local Variables ----!
      real(kind=cp) :: Ttr      ! transition temperature at this pressure

      !> default to 'high' phase for safety
      ip=.false.

      !> Check for valid model number. If not valid, return (with high phase indicated).
      if (eos%itran < 1 .or. eos%itran > N_TRANS_MODELS) return

      !> Test P and T against the Tr, set ip=.true. if in the low phase field for highT =high symm
      select case(eos%itran)
         case (1) ! Landau PV
            if (P < eos%params(21)) ip=.true.

         case (2) ! Landau TV
            if (T < eos%params(21)) ip=.true.

         case (3) ! Landau PVT
            !Ttr = eos%params(21)+p*eos%params(22)   changed from this on 18/12/2015
            Ttr = Get_Transition_Temperature(P,Eos)  ! general case
            if ( T < Ttr ) ip=.true.
      end select

      !> Now invert result if the lowT phase is high symm phase:
      if (eos%params(20) < 0) ip = .not. ip

      return
   End Function Transition_Phase

   !!----
   !!---- FUNCTION VSCALEMGD
   !!----
   !!---- Date: 03/02/2021
   !!
   Function VscaleMGD(EoS) Result(MGD)
      !---- Arguments ----!
      type(Eos_Type), intent(in)  :: EoS          ! EoS
      logical                     :: MGD        ! .true. if e%vscale_name is cm3/mol

      !---- Local Variables ----!
      character(len=len(eos%vscale_name)) :: vname

      MGD=.false.

      vname=adjustl(U_case(eos%vscale_name))
      if (len_trim(vname) == 0)return

      if (index(vname,'CM') > 0 .and. index(vname,'3') > 0 .and. index(vname,'MOL') > 0) MGD=.true.

      return
   End Function VscaleMGD



   !---------------------!
   !---- SUBROUTINES ----!
   !---------------------!

   !!----
   !!---- SUBROUTINE ALLOCATE_EOS_DATA_LIST
   !!----
   !!----    Allocation of objet E of eos_list_data_type.
   !!----    This subroutine should be called before using an object of type eos_data_list.
   !!----
   !!---- Date: 17/07/2015
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
   !!---- Date: 17/07/2015
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

   !!--++
   !!--++ SUBROUTINE CHECK_AXIS
   !!--++
   !!--++ Returns .true. if length of axis can be calculated from eos in cell
   !!--++
   !!--++ Date: 15/12/2020  (does not appear to be used???)
   !!
   Subroutine Check_Axis(Cell, Axis, OK)
      !---- Arguments ----!
      type(eos_cell_type), intent(in out) :: cell
      type(axis_type),     intent(in)     :: axis
      logical                             :: ok

      !---- Local Variables ----!

      !> init
      ok=.true.

      !> set flags
      call Set_Cell_Types(cell)

      !> Test
      select case(axis%ieos)
         case (0:6)
            if (cell%loaded(axis%ieos) > 0) return

         case(-2)
            call Eos_Cell_Loaded_Check(cell)
            if (.not. warn_eos) return
      end select

      ok=.false.

      return
   End Subroutine Check_Axis

   !!----
   !!---- SUBROUTINE CHECK_SCALES
   !!----
   !!----  Check Scales definitions
   !!----
   !!---- Date: 03/02/2021
   !!
   Subroutine Check_Scales(EoS, Dat)
      !---- Arguments ----!
      type(Eos_Type),                     intent(in)  :: EoS     ! EoS
      type (eos_data_list_type),optional, intent(in)  :: Dat   ! data structure

      !---- Local Variables ----!
      character(len=40) :: name

      !> Init
      Call Init_Err_EoS()

      !> APL
      if (eos%imodel == 6)then
         if (len_trim(eos%pscale_name) == 0)then
            Warn_EoS=.true.
            if (eos%linear)then
               Warn_Eos_Mess='APL EoS must have a Pscale (and M0) in GPa'
            else
               Warn_Eos_Mess='APL EoS must have a Pscale (and K0) in GPa'
            end if
         end if

         if (len_trim(eos%vscale_name) == 0 .or. index(U_case(eos%Vscale_name),'A') == 0)then
            Warn_EoS=.true.
            if (len_trim(Warn_EoS_Mess) == 0)then
               if (eos%linear)then
                  Warn_Eos_Mess='APL EoS must have a Vscale and L0 in A'
               else
                  Warn_Eos_Mess='APL EoS must have a Vscale and V0 in A^3'
               end if
            else
               if (eos%linear)then
                  Warn_Eos_Mess=trim(Warn_Eos_Mess)//' and a Vscale and L0 in A'
               else
                  Warn_Eos_Mess=trim(Warn_Eos_Mess)//' and a Vscale and V0 in A^3'
               end if
            end if
         end if
      end if

      !> If MGD or q-compromise type thermal EoS, must have eos%pscale_name and eos%_Vscale_name
      if (eos%itherm == 7 .or. eos%itherm == 8)then
         if ( .not. pscaleMGD(Eos))then
            Warn_EoS=.true.
            Warn_Eos_Mess='EoS must have a Pscale in kbar or GPa'
         end if

         if ( .not. VscaleMGD(EoS))then
            Warn_EoS=.true.
            if (len_trim(Warn_EoS_Mess) == 0)then
               Warn_Eos_Mess='EoS must have a Vscale in cm3/mol'
            else
               Warn_Eos_Mess=trim(Warn_Eos_Mess)//' and a Vscale in cm3/mol'
            end if
         end if
         if (len_trim(Warn_EoS_Mess) /= 0)Warn_Eos_Mess=trim(Warn_Eos_Mess)//' set to get correct results. '
      end if

      !> End checks here if only eos present
      if (.not. present(dat))return

      !>For all EoS compare data and eos scales
      if (len_trim(eos%pscale_name) /= 0 .and. len_trim(dat%Pscale_name) /=0)then
         if (trim(u_case(adjustl(eos%pscale_name))) /= trim(u_case(adjustl(dat%Pscale_name))))then
            Warn_EoS=.true.
            if (len_trim(Warn_EoS_Mess) > 0) Warn_Eos_Mess=trim(Warn_Eos_Mess)//' And'
            Warn_Eos_Mess=trim(Warn_Eos_Mess)//' Pscales of data and EoS are different.'
         end if
      end if

      if (len_trim(eos%vscale_name) /= 0 )then
         if (eos%linear)then
            name=trim(u_case(adjustl(dat%Lscale_name)))
         else
            name=trim(u_case(adjustl(dat%Vscale_name)))
         end if

         if (len_trim(name) /= 0)then
            if (trim(u_case(adjustl(eos%vscale_name))) /= trim(name))then
               Warn_EoS=.true.
               if (len_trim(Warn_EoS_Mess) > 0) Warn_Eos_Mess=trim(Warn_Eos_Mess)//' And'
               Warn_Eos_Mess=trim(Warn_Eos_Mess)//' Vscales of data and EoS are different'
            end if
         end if
      end if

      return
   End Subroutine Check_Scales

   !!----
   !!---- SUBROUTINE COPY_EOS_DATA_LIST
   !!----
   !!----
   !!---- Date: 03/02/2021
   !!
   Subroutine Copy_Eos_Data_List(Dat1, Dat2)
      !---- Arguments ----!
      type (eos_data_list_type), intent(in)      :: Dat1  ! Object to be copied
      type (eos_data_list_type), intent(in out)  :: Dat2  ! Output copy

      if (allocated(dat2%eosd))then
         call Deallocate_EoS_Data_List(dat2)
      end if

      call Allocate_EoS_Data_List(dat1%N, dat2)
      dat2=dat1

      return
   End Subroutine Copy_Eos_Data_List

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

   !!----
   !!---- SUBROUTINE DEF_CRYSTAL_SYSTEM
   !!----
   !!----   Either sets cell parameters to conform to specifed system Or, if no system specified,
   !!----   tries to determine system if all cell parameters provided
   !!----
   !!---- Update: 30/10/2020
   !!
   Subroutine Def_Crystal_System(Dat)
      !---- Arguments ----!
      type (eos_data_list_type),  intent(in out)   :: dat  ! data structure

      !---- Local Variables ----!
      character(len=40)       :: Family, SystemC, system
      character(len=1)        :: Symbol,U
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
               i=index(system,'-un')
               U='B'        !default is b-unique
               if (i > 0)U=U_case(system(i-1:i-1))
               select case(U)
                  case('A')
                     !> beta = gamma = 90º
                     dat%eosd(1:ndat)%ang(2)=90.0
                     dat%eosd(1:ndat)%ang(3)=90.0
                     dat%eosd(1:ndat)%siga(2)=0.0
                     dat%eosd(1:ndat)%siga(3)=0.0
                  case('B')
                     !> alpha = gamma = 90º
                     dat%eosd(1:ndat)%ang(1)=90.0
                     dat%eosd(1:ndat)%ang(3)=90.0
                     dat%eosd(1:ndat)%siga(1)=0.0
                     dat%eosd(1:ndat)%siga(3)=0.0
                  case('C')
                     !> alpha = beta = 90º
                     dat%eosd(1:ndat)%ang(1)=90.0
                     dat%eosd(1:ndat)%ang(2)=90.0
                     dat%eosd(1:ndat)%siga(1)=0.0
                     dat%eosd(1:ndat)%siga(2)=0.0
               end select

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
               ! err_eos=.true.
               ! write(Err_EoS_Mess,'(a,i5,a)')"Cell values for data #",i," in the data file indicate a change in crystal system"
               dat%system=' '
               return
            end if
         end do
      end if

      return
   End Subroutine Def_Crystal_System
   !!----
   !!---- SUBROUTINE EOS_CELL_LOADED_CHECK
   !!----
   !!---- Subroutine to check whether required eos are available for cell calculations
   !!---- If not, issues a warning
   !!----
   !!---- Date: 09/09/2020
   !!
   Subroutine Eos_Cell_Loaded_Check(Cell_Eos)
      !---- Arguments ----!
      type(eos_cell_type), intent(in) :: cell_eos

      !---- Variables ----!
      integer :: i, isum

      !>Init
      Warn_EOS=.true.

      select case(U_case(cell_eos%system(1:4)))
         case('CUBI','ISOT')
            if (cell_eos%loaded(0) == 0)then
               Warn_EOS_Mess='No eos loaded for cubic system, so no calculations possible'
               return
            end if

         case('ORTH','TRIG','HEXA','TETR')
            !>check for all eos present in some form
            do i = 0,3
               if (cell_eos%loaded(i) > 0)cycle
               if (index(U_case(cell_eos%system(1:4)),'ORTH') == 1)then
                  Warn_EOS_Mess='Orthorhombic crystal system, but not enough EoS loaded to do calculations'

               else
                  Warn_EOS_Mess='Uniaxial crystal system, but not enough EoS loaded to do calculations'
               end if

               return
            end do

         case('MONO')
            !> check for all eos present in some form
            isum=0
            do i = 0,3
               if (cell_eos%loaded(i) > 0)isum=isum+1
            end do

            if (isum < 4)then        !This takes care of angle poly because if present and 3 eos, the 4th is set as calculated
               Warn_EOS_Mess='Monoclinic crystal system, but not enough EoS loaded to do calculations'
               return
            end if

            !> check for unique flag
            if (cell_eos%unique == 0)then
               Warn_EOS_Mess='Monoclinic crystal system, but unique axis not defined'
               return
            end if

            if (cell_eos%eosang%iangle > 0)then
               isum=0
               do i = 0,3
                  if (cell_eos%loaded(i) == 1)isum=isum+1
               end do
               if (isum == 4)then
                  Warn_EOS_Mess='More than required EoS are loaded with angles: either delete angle polynomial or 1 EoS'
                  return
               end if
            end if

         case('TRIC')
            if (cell_eos%eosang%iangle == 0)then
               !> require all
               if (sum(cell_eos%loaded(0:6)) < 7)then
                  Warn_EOS_Mess='Triclinic crystal system, but not enough EoS loaded to do calculations'
                  return
               end if

            else
               isum=0
               do i = 0,3
                  if (cell_eos%loaded(i) == 1)isum=isum+1
               end do
               if (isum == 4)then
                  Warn_EOS_Mess='More than required EoS are loaded with angles: either delete angle polynomial or 1 EoS'
                  return

               else if(isum < 3)then
                  Warn_EOS_Mess='Triclinic crystal system, but not enough EoS loaded to do calculations'
                  return
               end if
            end if

         case default
            Warn_EOS_Mess='Unrecognised Crystal System'
            return

      end select

      Warn_EOS=.false.

      return
   End Subroutine Eos_Cell_Loaded_Check
   !!----
   !!---- SUBROUTINE EOSCAL_TEXT
   !!----
   !!----   Subroutine to write the calculated parameters of an eos to file at one P T point
   !!----   NO  header info is printed here.
   !!----
   !!----   Normally called from Write_eoscal
   !!----
   !!---- Date: 16/02/2021
   !!
   Subroutine Eoscal_Text(P, T, Tscale_In, Eos, text)
      !---- Arguments ----!
      real(kind=cp),     intent(in)  ::  p                   ! P to calculate properties
      real(kind=cp),     intent(in)  ::  t                   ! T to calculate properties
      character(len=*),  intent(in)  ::  tscale_in           ! Name of the Tscale for output, either C or K
      type(EoS_Type),    intent(in)  ::  eos                 ! Eos
      character(len=255),intent(out) :: text                 ! character string with results

      !---- Local variable ----!
      integer,parameter       :: nout=21
      character(len=1)        :: tscale   ! local name of tscale
      integer,dimension(nout) :: ip=(/6,6,9,8,6,5,  5, 9, 7, 7,    5,  9, 7,7,6,6,6,6,6,6,6/) ! format for output
      integer                 :: i

      real(kind=cp),dimension(6)   :: parvals
      real(kind=cp),dimension(6)   :: esd
      real(kind=cp),dimension(nout):: parout,esdout
      real(kind=cp)                :: v0,fp,fs,agt


      !> Tscale for output: C or K
      if (len_trim(tscale_in) == 0)then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      call init_err_eos()
      esd=0.0_cp
      esdout=0.0_cp
      parout=0.0_cp

      !> Now do the calculations at P,T
      Parvals= EoS_Cal(P,T,eos)
      if (sum(eos%vcv) > tiny(0.0_cp)) esd= eos_cal_esd(P,T,eos)

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
      if (abs(p) < tiny(0.0) ) then
         call ffcal_eos(p,t,eos,fp,fs)      ! because F not defined numerically at P=0
         parout(9)=FP
         esdout(9)=0.0_cp
      else
         call ffcal_dat_esd(parvals(1),esd(1),V0,0.0_cp,P,0.0_cp,Eos, &          ! only esd input is esd(V) at this P
              parout(9),esdout(9),parout(8),esdout(8))
      end if

      !> dK/dT
      parout(10)=parvals(5)*eos%factor(8)
      esdout(10)=esd(5)*eos%factor(8)

      !> handle alpha
      parout(11)=parvals(6)*eos%alphafactor
      esdout(11)=esd(6)*eos%alphafactor

      !> spon strain
      if (eos%itran > 0) parout(12)=Get_Transition_Strain(P,T,Eos)

      !> density
      if (eos%density0 > tiny(0.0)) then
         parout(13)=eos%density0*eos%params(1)/parvals(1)
         parout(14)=parout(13)*esd(1)/parvals(1)
      end if

      !> Thermal pressure
      if (eos%pthermaleos .and. eos%itran ==0) parout(15)=p-get_pressure(parvals(1),eos%tref,eos)

      !> Report adiabatic properties
      if (eos%itherm > 0 .and. abs(eos%params(18)) > tiny(0.0)) then
         parout(19)=Get_Grun_th(p,t,eos)             !Gruneisen for Kt--> Ks
         agt=parvals(6)*parout(19)*T        ! Get_Grun knows about linear
         if (eos%linear) agt=3.0_cp*agt
         parout(18)=(1.0_cp+agt)*parvals(2) ! Ks/Ms
      end if

      !> MGD EoS parameters:
      if (eos%itherm == 7 .or. eos%itherm == 8) then
          parout(17)=get_DebyeT(parvals(1),eos)      !Debye or Einstein T
      end if

      !> Cp and CV - write these provided there is a thermal and eos model and non-zero gamma0
          if(eos%itherm /= 0 .and. eos%imodel /= 0 .and. (abs(eos%params(18)) >tiny(0._cp) .or. eos%osc_allowed))then
              if(VscaleMGD(Eos) .and. .not. eos%linear)then
                 parout(20)=get_cp(P,T,Eos)
                 parout(21)=get_cv(P,T,Eos)
              endif
      end if

      !> output this datum: dynamic formatting to text string
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

      if (eos%pthermaleos .and. eos%itran ==0)text=trim(text)//'  '//trim(rformat(parout(15),ip(15)))

      if (eos%itherm > 0 .and. abs(eos%params(18)) > tiny(0.)) &
         text=trim(text)//'  '//trim(rformat(parout(18),ip(18)))//'  '//trim(rformat(parout(19),ip(19)))
         !>Cp and CV - write these provided there is a thermal and eos model and non-zero gamma0

      if (eos%itherm == 7 .or.  eos%itherm == 8) &
         text=trim(text)//'  '//trim(rformat(parout(17),ip(17)))

      if(eos%itherm /= 0 .and. eos%imodel /= 0 .and. (abs(eos%params(18)) >tiny(0._cp) .or. eos%osc_allowed))then
              if(VscaleMGD(Eos) .and. .not. eos%linear) &
                text=trim(text)//'  '//trim(rformat(parout(20),ip(20)))//'  '//trim(rformat(parout(21),ip(21)))
      endif


      return
   End Subroutine Eoscal_Text

   !!--++
   !!--++ SUBROUTINE EOSCAL_TEXT_DIRECTION
   !!--++
   !!--++   Subroutine to write the calculated parameters of calculated direction to file at one PT point
   !!--++   NO  header info is printed here.
   !!--++
   !!--++   Normally called from Write_eoscal
   !!--++
   !!--++   Date: 09/2020
   !!
   Subroutine Eoscal_Text_Direction(P, T, Tscale_In, Cell_eos, Axis, Text)
      !---- Arguments ----!
      real(kind=cp),      intent(in)  ::  p                   !P to calculate properties
      real(kind=cp),      intent(in)  ::  t                   !T to calculate properties
      character(len=*),   intent(in)  ::  tscale_in           ! Name of the Tscale for output, either C or K
      type(axis_type),    intent(in)  :: axis                 ! The direction
      type(EoS_Cell_Type),intent(in)  :: cell_eos             ! The eos for all the cell
      character(len=255), intent(out) :: text                 ! character string with results

      !---- Local variable ----!
      integer,parameter         :: nout=8
      character(len=1)          :: tscale   ! local name of tscale
      integer,dimension(nout)   :: ip=(/6,6,8,8,6,6,7,7/) ! format for output
      integer                   :: i

      real(kind=cp),dimension(6)   :: parvals
      real(kind=cp),dimension(nout):: parout
      real(kind=cp)                :: v0

      !> Tscale for output: C or K
      if (len_trim(tscale_in) == 0)then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      call init_err_eos()
      parout=0.0_cp

      !> Now do the calculations at P,T
      if (axis%ieos == -2)then
         Parvals= get_props_general(P,T,cell_eos,axis)
         V0=get_Volume_general(0.0_cp,T,cell_eos,axis)
      else
         Parvals= get_props_third(P,T,cell_eos,axis%ieos)
         V0=get_volume_third(0.0_cp,T,cell_eos,axis%ieos)
      end if

      !> build ouput value array
      parout(1)=p
      parout(2)=t
      if (tscale =='C')parout(2)=parout(2)-273.16
      parout(3)=parvals(1)*cell_eos%eosc%factor(1)      ! v
      parout(4)=parvals(1)/v0                 ! v/V0 at this T

      parout(5)=parvals(2)*cell_eos%eosc%factor(2)   !K or M
      parout(6)=parvals(3)*cell_eos%eosc%factor(3)   !Kp or Mp
      parout(7)=parvals(5)*cell_eos%eosc%factor(5)   !dK/dT or dM/dT
      parout(8)=parvals(6)*cell_eos%eosc%factor(10)   !alpha

      !> output this datum: dynamic formatting to text string
      !>init
      text=''

      !> No esd's
      do i=1,nout
         text=trim(text)//'  '//trim(rformat(parout(i),ip(i)))
      end do

      return
   End Subroutine Eoscal_Text_Direction

   !!----
   !!---- SUBROUTINE EOSPARAMS_CHECK
   !!----
   !!---- Check for Params that are invalid for all cases.
   !!----
   !!--.. NOTE:
   !!--..      Written 7-April-2014 to fix invalid eos parameters and report as error message
   !!--..      This is different from physical_check, because that checks properties at specific
   !!--..      p and T
   !!----      Later added warning state for non-fatal problems
   !!---- Date: 11/07/2016
   !!
   Subroutine EoSParams_Check(EoS)
      !---- Argument ----!
      type (EoS_Type), intent(in out) :: EoS

      !---- Local Variables ----!
      integer           :: i
      real(kind=cp)     :: pinf
      character(len=80) :: text

      !> Init
      call init_err_eos()

      !> Check for valid model numbers
      if (EoS%imodel < -1 .and. EoS%imodel > N_PRESS_MODELS)then
         err_eos=.true.
         err_eos_mess=' Invalid number for type of compressional eos'
      end if

      if (EoS%itherm < -1 .and. EoS%itherm > N_THERM_MODELS)then
         err_eos=.true.
         err_eos_mess=' Invalid number for type of thermal model'
      end if

      if (EoS%itran < -1 .and. EoS%itran > N_TRANS_MODELS)then
         err_eos=.true.
         err_eos_mess=' Invalid number for type of phase transition model'
      end if

      if (EoS%ishear < 0 .and. EoS%ishear > N_SHEAR_MODELS)then
         err_eos=.true.
         err_eos_mess=' Invalid number for type of shear modulus model'
      end if

      if (EoS%icross < 0 .and. EoS%icross > N_CROSS_MODELS)then
         err_eos=.true.
         err_eos_mess=' Invalid number for type of PT cross-terms model'
      end if

      if (EoS%iangle < 0 .and. EoS%iangle > N_ANGLE_MODELS)then
         err_eos=.true.
         err_eos_mess=' Invalid number for type of angle polynomial'
      end if

      !> Check that v0 is positive
      if (EoS%params(1) < tiny(0.0) .and. EoS%iangle == 0) then
         EoS%params(1)=1.0_cp
         err_eos=.true.
         if (EoS%linear) then
            err_eos_mess=' a0 was < 0. Not allowed! Reset to 1.00'
         else
            err_eos_mess=' V0 was < 0. Not allowed! Reset to 1.00'
         end if
      end if

      !> Check K0 is positive for V-P (-ve K ok for linear)
      if (.not. EoS%linear .and. EoS%imodel > 0) then
         if (EoS%params(2) < tiny(0.0_cp) ) then
            EoS%params(2)=10.0_cp
            err_eos=.true.
            if (len_trim(err_eos_mess) == 0) then
               err_eos_mess=' K0 was < 0. Not allowed! Reset to 10.0'
            else
               err_eos_mess=trim(err_eos_mess)//' And K0 was < 0. Not allowed! Reset to 10.0'
            end if
         end if
      end if


      !> Check that Z> 0 for APL EoS
      if (EoS%imodel ==6) then
         if (EoS%params(5) < tiny(0.0_cp) ) then
            EoS%params(5)=1.0_cp
            err_eos=.true.
            if (len_trim(err_eos_mess) == 0) then
               err_eos_mess=' Z was < 0. Not allowed! Reset to 1.0'
            else
               err_eos_mess=trim(err_eos_mess)//' Z was < 0. Not allowed! Reset to 1.0'
            end if
         end if
      end if



      !> Check Tref is positive
      if (EoS%tref < -1.0_cp*tiny(0.0_cp)) then
         EoS%tref=0.0_cp
         err_eos=.true.
         if (len_trim(err_eos_mess) == 0) then
            err_eos_mess=' Tref was < 0. Not allowed! Reset to 0 K'
         else
            err_eos_mess=trim(err_eos_mess)//' And Tref was < 0. Not allowed! Reset to 0 K'
         end if
      end if

      !> Thermal cases
      select case(EoS%itherm)  ! for specific thermal parameters
         case (4,6,7,8)    !>Kroll orPthermal must have characteristic T > 0.
            if (EoS%params(11) < 0.1) then
               EoS%params(11)=EoS%Tref
               if (EoS%Tref < 0.1) EoS%params=0.1
               err_eos=.true.
               if (len_trim(err_eos_mess) == 0) then
                  err_eos_mess=trim(EoS%comment(11))//' was =< 0. Not allowed! Reset to Tref'
               else
                  err_eos_mess=trim(err_eos_mess)//' And '//trim(EoS%comment(11))//' was =< 0. Not allowed! Reset to Tref'
               end if
            end if

            if(eos%itherm == 7 .or. eos%itherm ==8)then !thermal P require Natom, only Cv and alpha need this in HP2011
               if(eos%params(13) < 1.0)then
                   err_eos=.true.
                   if (len_trim(err_eos_mess) == 0) then
                        err_eos_mess='Natom < 1.0 not valid: reset it!'
                    else
                        err_eos_mess=trim(err_eos_mess)//' And Natom < 1.0 not valid: reset it!'
                    end if
               endif
               if ( .not. pscaleMGD(Eos))then
                   err_eos=.true.
                   if (len_trim(err_eos_mess) == 0) then
                        err_eos_mess='Pscale must be GPa or kbar'
                    else
                        err_eos_mess=trim(err_eos_mess)//' And Pscale must be GPa or kbar'
                    end if
               endif
               if ( .not. vscaleMGD(Eos))then
                   err_eos=.true.
                   if (len_trim(err_eos_mess) == 0) then
                        err_eos_mess='Vscale must be cm^3/mol'
                    else
                        err_eos_mess=trim(err_eos_mess)//' And Vscale must be cm^3/mol'
                    end if
               endif
             elseif(eos%itherm == 6)then !only Cv and alpha need this in HP2011
               if(eos%params(13) < 1.0)then
                   warn_eos=.true.
                   if (len_trim(warn_eos_mess) == 0) then
                        warn_eos_mess='Natom =0 so PVT  correct, but not heat capacities or calculated alpha'
                    else
                        warn_eos_mess=trim(warn_eos_mess)//' And Natom =0 so PVT  correct, but not heat capacities or calculated alpha'
                    endif
               end if

            endif


         end select

      !> Check q-comp switch
      if (EoS%itherm == 7 .or. EoS%itherm == 8)then
         if (EoS%params(14) < 0._cp)then
            EoS%params(14)=0.0_cp

         else if(EoS%params(14) > 1._cp)then
                 EoS%params(14)=1.0_cp
         end if
      end if

      !> Extra oscillator models
      if (EoS%iosc(1) >0 .and. EoS%params(41) < 0.1_cp)then
         EoS%params(41)=EoS%Tref
         err_eos=.true.
         if (len_trim(err_eos_mess) == 0) then
             err_eos_mess=trim(EoS%parname(41))//'for 2nd oscillator was =< 0. Not allowed! Reset to Tref'
         else
            err_eos_mess=trim(err_eos_mess)//' And '//trim(EoS%parname(41))//'for 2nd oscillator was =< 0. Not allowed! Reset to Tref'
         end if
      end if

      if (EoS%iosc(2) >0 .and. EoS%params(46) < 0.1_cp)then
         EoS%params(46)=EoS%Tref
         err_eos=.true.
         if (len_trim(err_eos_mess) == 0) then
            err_eos_mess=trim(EoS%parname(46))//'for 3rd oscillator was =< 0. Not allowed! Reset to Tref'
         else
            err_eos_mess=trim(err_eos_mess)//' And '//trim(EoS%parname(46))//'for 3rd oscillator was =< 0. Not allowed! Reset to Tref'
         end if
      end if

      if (EoS%iosc(1) > 0 .and. EoS%params(40)+EoS%params(45) > 1.0)then
         err_eos=.true.
         if (len_trim(err_eos_mess) == 0) then
             err_eos_mess='Mode fractions of extra oscillators reset to sum to 1'
         else
            err_eos_mess=trim(err_eos_mess)//' And mode fractions of extra oscillators reset to sum to 1'
         end if

         if (EoS%params(40) > 1.0_cp)then
            EoS%params(40)=1.0_cp
            EoS%params(45)=0.0_cp
         else
            EoS%params(45)=1.0_cp-EoS%params(40)
         end if
      end if

      !> Scale factors for data must be positive
      text=' '
      do i=51,59
         if (EoS%iuse(i) > 0 .and. EoS%params(i) < tiny(0._cp))then
            err_eos=.true.
            text='Data scale factor =< 0., reset to 1.0'
            EoS%params(i)=1.0_cp
         end if
      end do

      if (len_trim(text) /= 0)then
         if (len_trim(err_eos_mess) == 0)then
            err_eos_mess=trim(text)
         else
            err_eos_mess=trim(err_eos_mess)//' and '//trim(text)
         end if
      end if

      !> Produce warning for curved phase boundaries: Pinflection = a/-2b when Ttr=Tr0+aP+bP^2
      if (EoS%itran>1 .and. abs(EoS%params(24)) > tiny(0.0) ) then
         pinf=-1.0*EoS%params(22)/2.0/EoS%params(24)
         if (abs(pinf) < 10.0) then
            err_eos=.true.
            write(text,'(a,f4.1,1x,a)')'Phase boundary inflects at P ~ ',pinf,trim(EoS%pscale_name)
            if (len_trim(err_eos_mess) == 0) then
               err_eos_mess=trim(text)
            else
               err_eos_mess=trim(err_eos_mess)//' And '//trim(text)
            end if
         end if
      end if

      !> If MGD and linear warn that this is not generally valid
      if ((EoS%itherm == 6 .or. EoS%itherm == 7 .or. EoS%itherm == 8) .and. EoS%linear) then
         err_eos=.true.
         text='Linear EoS only has valid parameters if the material is cubic'
         if (len_trim(err_eos_mess) == 0) then
            err_eos_mess=' *****WARNING: '//trim(text)
        else
           err_eos_mess=trim(err_eos_mess)//' And '//trim(text)
        end if
      end if

      return
   End Subroutine EoSParams_Check

   !!----
   !!---- SUBROUTINE FFCAL_DAT
   !!----
   !!---- Return Normalized pressure and/or Strain at V and P
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine FfCal_Dat(V, V0, P, EoS, F, S)
      !---- Arguments ----!
      real(kind=cp),           intent(in)  :: V    ! Volume
      real(kind=cp),           intent(in)  :: V0   ! Volume at zero pressure
      real(kind=cp),           intent(in)  :: P    ! Pressure
      type(Eos_Type),          intent(in)  :: EoS  ! Eos Parameter: only imodel and linear used
      real(kind=cp), optional, intent(out) :: F    ! Normalised pressure
      real(kind=cp), optional, intent(out) :: S    ! Strain

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
      sc=Strain(VV0,Eos)                      ! returns volume like strain for linear
      fc=NormPressure_P(sc,P,EoS%imodel)

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
   Subroutine FfCal_Dat_Esd(V, SigV, V0, SigV0, P, SigP, Eos, F, SigF, S, SigS)
      !---- Arguments ----!
      real(kind=cp),  intent(in)  :: V       ! Volume
      real(kind=cp),  intent(in)  :: SigV    ! Sigma (V)
      real(kind=cp),  intent(in)  :: V0      ! Vo
      real(kind=cp),  intent(in)  :: SigV0   ! Sigma (Vo)
      real(kind=cp),  intent(in)  :: P       ! Pressure
      real(kind=cp),  intent(in)  :: SigP    ! Sigma (P)
      type(Eos_Type), intent(in)  :: EoS     ! Eos Parameter
      real(kind=cp),  intent(out) :: F       ! Normalised pressure
      real(kind=cp),  intent(out) :: SigF    ! Sigma (F)
      real(kind=cp),  intent(out) :: S       ! Strain
      real(kind=cp),  intent(out) :: SigS    ! Sigma (S)

      !---- Local Variables ----!
      real(kind=cp) :: vv0,vv0_esd, sigmap,e

      !> Init
      f=0.0_cp; sigf=0.0_cp
      s=0.0_cp; sigs=0.0_cp

      if (abs(v0) <= 0.00001) then
         err_eos=.true.
         err_eos_mess='V0 is zero - No strain can be calculated!'
         return
      end if

      !> Note that V0 in call is the V0 at this temperature
      vv0=v/v0
      vv0_esd=vv0*sqrt( (sigv/abs(v))**2.0_cp + (sigv0/abs(v0))**2.0_cp )

      !> Strain
      s=Strain(vv0,Eos)      ! input a/a0 or v/v0 to strain. It always returns volume strain

      !> Normalized Pressure
      f=NormPressure_P(s,p,EoS%imodel)

      !> Check Pressure values
      if (abs(p) <= 0.0001) then
         e=0.0_cp
      else
         e=sigp/p
      end if

      !> If the finite strain is zero, the esd of the normalised pressure is not defined
      if (s < tiny(0.0)) return

      !> ESD Calculations: all done in volume terms
      if (eos%linear) then
         vv0_esd=3.0*vv0**2.0_cp*vv0_esd
         vv0=vv0**3.0_cp
      end if

      select case (eos%imodel)
         case (1,5,6,7) ! Murnaghan, Tait, APL, Kumar
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
   !!---- SUBROUTINE FFCAL_EOS
   !!----
   !!----  Returns normalised pressure and Strain at this P,T,
   !!----  calculated from the eos parameters
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine FfCal_EoS(P, T, EoS, F, S)
      !---- Arguments ----!
      real(kind=cp),           intent(in)     :: P    ! Pressure
      real(kind=cp),           intent(in)     :: T    ! Temperature
      type(Eos_Type),          intent(in)     :: Eos  ! Eos Parameter
      real(kind=cp),           intent(out)    :: F    ! Normalised pressure
      real(kind=cp),           intent(out)    :: S    ! Strain

      !---- Local Variables ----!
      real(kind=cp) :: v

      !> Init
      f=0.0_cp
      s=0.0_cp

      v=get_volume(p,t,eos)         ! V at this P,T from eos parameters
      s=strain_eos(v,t,eos)         ! finite strain at this P,T
      f=normpressure_eos(s,t,eos)   ! Norm pressure

      return
   End Subroutine FfCal_EoS



   !!----
   !!---- SUBROUTINE GET_TENSOR_EOS
   !!----
   !!---- Returns the value of alpha or beta tensor calculated from eos of cell and
   !!---- its eigen vectors etc
   !!----
   !!---- Date: 22/10/2020
   !!
   Subroutine Get_Tensor_Eos(P, T, Cell_Eos, dx, x)
      !---- Arguments ----!
      real(kind=cp),           intent(in)     :: p,T
      type(eos_cell_type),     intent(in)     :: cell_eos
      character(len=1),        intent(in)     :: dx  ! ='T' for alpha or 'P' for beta
      type(Strain_Tensor_Type),intent(in out) :: x   ! only 'in' is cartype, this routine loads system and property and paxis angles

      !---- Local Variables ----!
      integer                     :: i
      real(kind=cp)               :: dr,cotbs,cotgs
      real(kind=cp),dimension(3)  :: d,da ! for derivatives 1/a. da/dP and dangle/dP
      type(axis_type)             :: axis
      type(crystal_cell_type)     :: c  !cell params, metric tensor at this P,T
      character(len=2)            :: cartype  !local copy
      character(len=2)            :: dtype  !local copy of P or T

      !> Init
      cartype=x%cartype
      call init_strain_tensor(x)
      x%cartype=U_case(cartype)
      x%system=cell_eos%system
      dtype='P'
      if (U_case(dx) == 'T')dtype='T'

      !> Calculate cell edge compressibilities
      do i = 1,3
         axis%ieos=i
         if (dtype == 'P')then
            d(i)=-1.0_cp/Get_Mod_Cell(P,T,cell_eos,axis)        !d(1) is 1/a . da/dP: the negative of the compressibility
         else
            d(i)=Get_alpha_cell(P,T,cell_eos,axis)
         end if
      end do

      !> now do triclinic or monoclinic
      if (index(U_case(x%system),'TRIC') > 0 .or. index(U_case(x%system),'MONO') > 0)then
         !> First calculate the cell params and recip cell at this point
         c=  get_params_cell(P,T,cell_eos)
         do i=1,3
            da(i)=get_angle_deriv(P,T,cell_eos,i,.true.,dtype)
         end do

         select case(x%cartype)
            case('BC')   !cartype=2: Redfern & Carpenter Y // b and Z //c*
               dr=get_angle_deriv(P,T,cell_eos,2,.false.,dtype)   !d(beta*)/dP;  needed if beta changing but is 90.0
               !tensor coeffs if beta*=90
               x%ep(1,1)=d(1) +da(3)/tand(c%ang(3))
               x%ep(2,2)=d(2)
               x%ep(3,3)=d(3) +da(1)/tand(c%ang(1))
               x%ep(1,3)=0.5_cp*dr
               x%ep(2,3)=0.5_cp*((d(3)-d(2))/tand(c%ang(1))/sind(c%rang(2)) - da(1)/sind(c%rang(2)))
               x%ep(1,2)=0.5_cp*((d(1)-d(2))/tand(c%ang(3)) - da(3))

               if (abs(c%rang(2)-90.0) > 0.01)then
                  !add in terms when beta* /=90
                  cotbs=1.0_cp/tand(c%rang(2))        !cot(beta*)
                  x%ep(3,3)=x%ep(3,3) + dr*cotbs
                  x%ep(2,3)=x%ep(2,3) + 0.5_cp*cotbs*((d(1)-d(2))/tand(c%ang(3))     -da(3))
                  x%ep(1,3)=x%ep(1,3) + 0.5_cp*cotbs*(d(1)-d(3) -da(1)*cosd(c%ang(1)) + da(3)/tand(c%ang(3)))
               end if

            case('BA')  !cartype=3: Brown and Angel, Equations from Tribaudino et al (2011)
               dr=get_angle_deriv(P,T,cell_eos,2,.false.,dtype)   !d(be*)/dP;  needed if beta changing but is 90.0
               !tensor coeffs if beta*=90
               x%ep(1,1)=d(1) +da(3)/tand(c%ang(3))
               x%ep(2,2)=d(2)
               x%ep(3,3)=d(3) + da(1)/tand(c%ang(1))
               x%ep(1,2)=0.5_cp*((d(1)-d(2))/tand(c%ang(3))/sind(c%rang(2)) - da(3)/sind(c%rang(2)))
               x%ep(1,3)=0.5_cp*dr
               x%ep(2,3)=0.5_cp*((d(3)-d(2))/tand(c%ang(1)) - da(1))
               if (abs(c%rang(2)-90.0) > 0.01)then
                  !add in terms when beta* /=90
                  cotbs=1.0_cp/tand(c%rang(2))        !cot(beta*)
                  x%ep(1,1)=x%ep(1,1) + dr*cotbs
                  x%ep(1,2)=x%ep(1,2) + 0.5_cp*cotbs*((d(3)-d(2))/tand(c%ang(1))     -da(1))
                  x%ep(1,3)=x%ep(1,3) + 0.5_cp*cotbs*(d(3)-d(1) -da(3)*cosd(c%ang(3)) + da(1)/tand(c%ang(1)))
               end if

            case('CB')      ! cartype=4: Neumann (1861) Equations from Pauffler and Weber (1999)
               dr=get_angle_deriv(P,T,cell_eos,3,.false.,dtype)   !d(ga*)/dP;  needed if gamma changing but is 90.0
               !tensor coeffs if gamma*=90
               x%ep(1,1)=d(1) +da(2)/tand(c%ang(2))
               x%ep(2,2)=d(2) +da(1)/tand(c%ang(1))
               x%ep(3,3)=d(3)
               x%ep(1,2)=0.5_cp*dr
               x%ep(1,3)=0.5_cp*(d(1)-d(3))/tand(c%ang(2)) - 0.5_cp*da(2)
               x%ep(2,3)=0.5_cp*((d(2)-d(3))/tand(c%ang(1))/sind(c%rang(3)) - da(1)/sind(c%rang(3)))
               if (abs(c%rang(3)-90.0) > 0.01)then
                  !add in terms when gamma* /=90
                  cotgs=1.0_cp/tand(c%rang(3))        !cot(gamma*)
                  x%ep(2,2)=x%ep(2,2) + dr*cotgs
                  x%ep(1,2)=x%ep(1,2) + 0.5_cp*cotgs*(d(1)-d(2) -da(1)*cosd(c%ang(1)) + da(2)/tand(c%ang(2)))
                  x%ep(2,3)=x%ep(2,3) + 0.5_cp*cotgs*((d(1)-d(3))/tand(c%ang(2))     -da(2))
               end if

            case default ! CA cartype=1:  This is Z //C X//A*: IRE convention
               ! Invalid orientation code defaults to this one
               ! These equations derived by RJA, October 2020
               dr=get_angle_deriv(P,T,cell_eos,3,.false.,dtype)  !d(ga*)/dP;  needed if gamma changing but is 90.0
               !tensor coeffs if gamma*=90
               x%ep(2,2)=d(2) +da(1)/tand(c%ang(1))
               x%ep(1,1)=d(1) +da(2)/tand(c%ang(2))
               x%ep(3,3)=d(3)
               x%ep(1,2)=0.5_cp*dr
               x%ep(2,3)=0.5_cp*((d(2)-d(3))/tand(c%ang(1)) - da(1))
               x%ep(1,3)=0.5_cp*((d(1)-d(3))/tand(c%ang(2))/sind(c%rang(3)) - da(2)/sind(c%rang(3)))
               if (abs(c%rang(3)-90.0) > 0.01)then
                  !add in terms when gamma* /=90
                  cotgs=1.0_cp/tand(c%rang(3))        !cot(gamma*)
                  x%ep(1,1)=x%ep(1,1) + dr*cotgs
                  x%ep(1,2)=x%ep(1,2) + 0.5_cp*cotgs*(d(2)-d(1) -da(2)*cosd(c%ang(2)) + da(1)/tand(c%ang(1)))
                  x%ep(1,3)=x%ep(1,3) + 0.5_cp*cotgs*((d(2)-d(3))/tand(c%ang(1))     -da(1))
               end if

         end select

         !finish
         x%ep(2,1)=x%ep(1,2)
         x%ep(3,1)=x%ep(1,3)
         x%ep(3,2)=x%ep(2,3)
         if (dtype == 'P')then
            x%ep=-1000.0_cp*x%ep       ! because compressibilities are negative of 1/a da/dP etc
            if (len_trim(cell_eos%eosc%Pscale_name) > 0)then
               x%property='Compressibility in units of inverse '//trim(cell_eos%eosc%Pscale_name)//' x 10^3'
            else
               x%property='Compressibility in units of inverse pressure units x 10^3'
            end if

         else
            x%property='Thermal expansion x 10^5'
            x%ep=100000.0_cp*x%ep
         end if

         call fix_tensor(x%ep,x%system)   ! make strain conform to crystal system, and thus eliminate round-off error

         !> for monoclinic or triclinic calculate Eigenvalues and vectors from tensor of properties x%ep
         call Diagonalize_SH (X%Ep, 3, X%evalp, X%Evec)
         call orient_eigenvectors(X%evalp,X%evec)       !sort the eigen vectors so that #1 is close to +X etc
         call calc_Paxes_angles(x,c,3)

      else        !higher symmetries
         do i=1,3
            x%ep(i,i)=d(i)
         end do
         if (dtype == 'P')then
            x%ep=-1000.0_cp*x%ep       ! because compressibilities are negative of 1/a da/dP etc
            if (len_trim(cell_eos%eosc%Pscale_name) > 0)then
               x%property='Compressibility in units of inverse '//trim(cell_eos%eosc%Pscale_name)//' x 10^3'
            else
               x%property='Compressibility in units of inverse pressure units x 10^3'
            end if

         else
            x%property='Thermal expansion x 10^5'
            x%ep=100000.0_cp*x%ep
         end if
         call fix_tensor(x%ep,x%system)   ! make strain conform to crystal system, and thus eliminate round-off error
      end if

      return
   End Subroutine Get_Tensor_Eos

   !!----
   !!---- SUBROUTINE INIT_EOS_Angles
   !!----
   !!---- Initialize the EoS Type for Angles polynomial
   !!----
   !!---- Date: 07/10/2020
   !!
   Subroutine Init_EoS_Angles(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS

      !> Check for valid model number. If not valid, set zero
      if (eos%iangle < 0 .or. eos%iangle > N_ANGLE_MODELS) eos%iangle=0

      Eos%AngPoly  = 0.0_cp
      Eos%angpoly(1:3,0,1)=90._cp

      return
   End Subroutine Init_EoS_Angles

   !!----
   !!---- SUBROUTINE INIT_EOS_CELL_TYPE
   !!----
   !!---- Subroutine to initialise eos_cell_type and set to default orthorhombic
   !!----
   !!---- Date: 09/09/2020
   !!
   Subroutine Init_Eos_Cell_Type(Cell_Eos)
      !---- Arguments ----!
      type(eos_cell_type),intent(in out) :: cell_eos

      !---- Local Variables ----!
      integer i

      !> clear eos
      do i=0,6
         call Init_EoS_Type(cell_eos%eos(i))
      end do

      !>reset to default orthorhombic
      cell_eos%n=3
      cell_eos%system='ORTHORHOMBIC'
      cell_eos%obtuse = .true.
      cell_eos%unique_label=' '
      cell_eos%unique=0
      call Init_EoS_Type(cell_eos%eosang)
      call set_cell_types(cell_eos)      ! sets the cout array PV VT etc

      return
    end subroutine Init_Eos_Cell_Type

   !!----
   !!---- SUBROUTINE INIT_EOS_CROSS
   !!----
   !!---- Initialize the EoS Type for P-T cross-terms
   !!----
   !!---- Date: 11/11/2016
   !!
   Subroutine Init_EoS_Cross(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS

      !> Check for valid model number. If not valid, set zero
      if (EoS%icross < 0 .or. EoS%icross > N_CROSS_MODELS) EoS%icross=0
      if (EoS%pthermaleos) EoS%icross=0

      EoS%params(8)           = 0.0_cp
      EoS%params(9)           = 0.0_cp
      EoS%vcv(8:9,1:N_EOSPAR) = 0.0_cp
      EoS%vcv(1:N_EOSPAR,5:6) = 0.0_cp
      EoS%factor(8)           = 1.0_cp
      EoS%factor(9)           = 1.0_cp

      call Set_Cross_Names(EoS)    ! Set the variable names
      call Set_Eos_Use(EoS)        ! update the use flags

      return
   End Subroutine Init_EoS_Cross

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
   !!---- SUBROUTINE INIT_EOS_GROUPSCALES
   !!----
   !!---- Initialize the EoS Type for extra oscillators
   !!----
   !!---- Date: 23/03/2020
   !!
   Subroutine Init_EoS_GroupScales(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS

      !---- Variables ----!
      integer  :: i

      !>Default values
      EoS%params(50:59) = 1.0_cp
      EoS%factor(50:59) = 1.0_cp
      EoS%vcv(50:59,1:N_EOSPAR)= 0.0_cp
      EoS%vcv(1:N_EOSPAR,50:59)= 0.0_cp

      !>Names
      EoS%ParName(1) ='     '
      EoS%comment(1) ='Not used'
      do i=1,9
        write(EoS%ParName(50+i),'(''Sca '',i1)')i
        write(EoS%comment(50+i),'(''Scale factor for data group '',i1)')i
      enddo

      !>Use flags default to 0. The values depend on the groups of data present
      ! and therefore cannot be set in cfml_eos_mod but must be reset by main programs
      EoS%iuse(50:59)=0

      return
   End Subroutine Init_EoS_GroupScales

   !!----
   !!---- SUBROUTINE INIT_EOS_OSC
   !!----
   !!---- Initialize the EoS Type for extra oscillators
   !!----
   !!---- Date: 09/03/2020
   !!
   Subroutine Init_EoS_Osc(EoS, i)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS
      integer,         intent(in)     :: i     !=1 for 1st oscillator, =2 for second, = 3 for both


      if (i == 1 .or. i == 3) then
         if (eos%iosc(1) < 0 .or. eos%iosc(1) > N_OSC_MODELS) eos%iosc(1)=0

         !> initial values same for all models
         eos%params(40)        = 0.0_cp     !fraction
         eos%params(41)        = eos%tref  ! Characteristic T
         eos%params(42)        = 1.0_cp     ! gamma
         eos%params(43)        = 0.0_cp     ! q
         eos%params(44)        = 0.0_cp     ! not used

         eos%vcv(40:44,1:N_EOSPAR)= 0.0_cp
         eos%vcv(1:N_EOSPAR,40:44)= 0.0_cp

         eos%factor(40:44)        = 1.0_cp
      end if

      if (i == 2 .or. i == 3)then
         if (eos%iosc(2) < 0 .or. eos%iosc(2) > N_OSC_MODELS) eos%iosc(2)=0

         !> initial values same for all models
         eos%params(45)        = 0.0_cp     !fraction
         eos%params(46)        = eos%tref  ! Characteristic T
         eos%params(47)        = 1.0_cp     ! gamma
         eos%params(48)        = 0.0_cp     ! q
         eos%params(49)        = 0.0_cp     ! not used

         eos%vcv(45:49,1:N_EOSPAR)= 0.0_cp
         eos%vcv(1:N_EOSPAR,45:49)= 0.0_cp

         eos%factor(45:49)        = 1.0_cp
      end if

      call Set_Osc_Names(Eos)    ! Set the variable names
      call Set_Eos_Use(Eos)        ! update the use flags

      return
   End Subroutine Init_EoS_Osc

   !!----
   !!---- SUBROUTINE INIT_EOS_Shear
   !!----
   !!---- Initialize the EoS Type for Shear case
   !!----
   !!---- Date: 17/02/2015
   !!
   Subroutine Init_EoS_Shear(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS

      !---- Variables ----!
      integer  :: n

      !> Check for valid model number. If not valid, set zero
      if (EoS%ishear < 0 .or. EoS%ishear > N_SHEAR_MODELS) EoS%ishear=0

      !> Set upper limit to parameter numbers
      n=34
      if (n > N_EOSPAR) n=N_EOSPAR

      select case(EoS%ishear)
         case (0)
            EoS%params(30)          = huge(0.0_cp)   ! default is infinitely stiff
            EoS%vcv(30,1:n)         = 0.0_cp
            EoS%vcv(30:N_EOSPAR,1:n)= 0.0_cp
            EoS%vcv(1:N_EOSPAR,30:n)= 0.0_cp
            EoS%factor(30:n)        = 1.0_cp

         case (1)        ! polynomial
            EoS%params(30)    = 100.0_cp      ! G0 at Pref Tref
            EoS%params(31:34) =   0.0_cp      ! Polynomial coefficients
            EoS%factor(30:n)  =   1.0_cp
      end select

      call Set_Shear_Names(EoS)    ! Set the variable names
      call Set_Eos_Use(EoS)        ! update the use flags

      return
   End Subroutine Init_EoS_Shear

   !!----
   !!---- SUBROUTINE INIT_EOS_THERMAL
   !!----
   !!---- Initialize the EoS Type for Thermal case
   !!----
   !!---- Date: 10/09/2013
   !!
   Subroutine Init_EoS_Thermal(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS

      !---- Variables ----!
      integer    :: n

      !> Check for valid model number. If not valid, set zero
      if(EoS%itherm < -1 .or. EoS%itherm > N_THERM_MODELS) EoS%itherm=0

      !> Set upper limit to thermal parameter numbers
      n=19
      if (n > N_EOSPAR)n=N_EOSPAR

      EoS%alphafactor=1.0E5_cp                      ! Normal scale factor for printing values of alpha

      select case(EoS%itherm)
         case (-1)           ! PTV table
            EoS%factor(10)   = 1.0E5_cp            ! factor to multiply alpha values on printing
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (0)
            EoS%params(8)    = 0.0_cp
            EoS%params(10:n) = 0.0_cp
            EoS%vcv(8,1:n)   = 0.0_cp
            EoS%vcv(10:n,1:n)= 0.0_cp
            EoS%vcv(1:n,10:n)= 0.0_cp
            EoS%factor(10:n) = 1.0_cp
            EoS%TRef         = 298.0_cp
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (1)
            EoS%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            EoS%factor(11)  = 1.0E8_cp
            EoS%TRef        = 298.0_cp             ! Simple thermal expansion,
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (2)
            EoS%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            EoS%factor(11)  = 1.0E8_cp
            EoS%factor(12)  = 1.0_cp
            EoS%TRef        = 298.0_cp             ! Simple thermal expansion,
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (3)
            EoS%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            EoS%factor(11)  = 1.0E4_cp
            EoS%TRef        = 298.0_cp             ! Simple thermal expansion,
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (4)
            EoS%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            EoS%factor(11)  = 1.0_cp
            EoS%TRef        = 298.0_cp             ! Holland and Powell thermal expansion without P
            EoS%TRef_fixed  = .false.
            EoS%params(11)  = 298.0_cp             ! Einstein temperature
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (5)
            EoS%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            EoS%factor(11)  = 1.0_cp
            EoS%TRef        = 0.0_cp               ! Salje thermal expansion
            EoS%TRef_fixed  = .true.
            EoS%params(11)  = 298.0_cp             ! Saturation temperature
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (6)
            EoS%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            EoS%factor(11)  = 1.0_cp
            EoS%TRef        = 298.0_cp
            EoS%TRef_fixed  = .false.
            EoS%params(11)  = 298.0_cp             ! Einstein temperature default
            EoS%params(14)  = 1.0                   ! q-compromise
            EoS%pthermaleos  =.true.
            EoS%Osc_allowed  =.true.

         case(7,8)                                 ! MGD,  Einstein Osc:
            EoS%factor(10:14)  = 1.0_cp            ! plus gamma0 and q as (18),(19)
            EoS%TRef           = 298.0_cp
            EoS%TRef_fixed     = .false.

            EoS%params(11)     = 298.0_cp          ! Debye/Einstein temperature default

            EoS%params(13)     = 1.0               ! Natoms/molecule for MGD
            EoS%params(14)     = 0.0               ! flag to use full q
            EoS%pthermaleos    =.true.
            EoS%Osc_allowed    =.true.

      end select

      !> Set the common terms for Ks to Kt conversion: also used in  thermal pressure with oscillator
      EoS%params(18)=1.0_cp      ! gamma0
      EoS%params(19)=0.0_cp      ! q

      call Init_EoS_Cross(EoS)                     ! init the cross-terms
      if (.not. EoS%osc_allowed)call Init_EoS_Osc(EoS,3)  ! clear extra oscillators
      call Set_Thermal_Names(EoS)                  ! Set the variable names
      call Set_Eos_Use(EoS)                        ! update the use flags and other pointers

      return
   End Subroutine Init_EoS_Thermal

   !!----
   !!---- SUBROUTINE INIT_EOS_TRANSITION
   !!----
   !!---- Initialize the EoS Type for Transition case
   !!----
   !!---- Date: 17/02/2015
   !!
   Subroutine Init_EoS_Transition(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS

      !---- Variables ----!
      integer  :: n

      !> Check for valid model number. If not valid, set zero
      if (EoS%itran < 0 .or. EoS%itran > N_TRANS_MODELS) EoS%itran=0

      !> Set upper limit to parameter numbers
      n=29
      if (n > N_EOSPAR) n=N_EOSPAR

      select case(EoS%itran)
         case (0)
            EoS%params(20:n)        = 0.0_cp
            EoS%vcv(20,1:n)         = 0.0_cp
            EoS%vcv(20:n_eospar,1:n)= 0.0_cp
            EoS%vcv(1:n_eospar,20:n)= 0.0_cp
            EoS%factor(20:n)        = 1.0_cp

         case (1)        ! Landau PV
            EoS%params(20) = 1.0_cp        ! high P is high sym
            EoS%params(21) = 5.0_cp        ! Safe default Ptr
            EoS%params(22) = 0.0_cp        ! dTr/dP: not used
            EoS%params(23) = 0.0_cp        ! d2Tr/dP2 not used
            EoS%params(24) = 0.0_cp        ! no excess V
            EoS%params(25) = 0.5_cp        ! power law
            EoS%params(26) = 0.0_cp        ! excess V high phase
            EoS%params(27) = 0.5_cp        ! power law high phase

            EoS%factor(20:n) = 1.0_cp
            EoS%factor(24)   = 1.0E3_cp    ! 1000 for aL
            EoS%factor(26)   = 1.0E3_cp    ! 1000 for aH

         case (2)        ! Landau TV
            EoS%params(20) =   1.0_cp        ! high T is high sym
            EoS%params(21) = 800.0_cp        ! Safe default Ttr
            EoS%params(22) = 100.0_cp        ! dTr/dP: dummy will be used in get_pressure
            EoS%params(23) =   0.0_cp        ! d2Tr/dP2 not used
            EoS%params(24) =   0.0_cp        ! no excess V
            EoS%params(25) =   0.5_cp        ! power law
            EoS%params(26) =   0.0_cp        ! excess V high phase
            EoS%params(27) =   0.5_cp        ! power law high phase

            EoS%factor(20:n) = 1.0_cp
            EoS%factor(24)   = 1.0E3_cp      ! 1000 for aL
            EoS%factor(26)   = 1.0E3_cp      ! 1000 for aH

         case (3)        ! Landau PVT
            EoS%params(20) =   1.0_cp        ! high T is high sym
            EoS%params(21) = 800.0_cp        ! Safe defaulat Tr
            EoS%params(22) =   1.0_cp
            EoS%params(23) =   0.0_cp
            EoS%params(24) =   0.0_cp        ! no excess V
            EoS%params(25) =   0.5_cp
            EoS%params(26) =   0.0_cp        ! excess V high phase
            EoS%params(27) =   0.5_cp        ! power law high phase

            EoS%factor(20:n) = 1.0_cp
            EoS%factor(24)   = 1.0E3_cp      ! 1000 for aL
            EoS%factor(23)   = 1.0E5_cp      ! for da/dP
            EoS%factor(26)   = 1.0E3_cp      ! 1000 for aH

      end select

      call Set_Transition_Names(EoS)    ! Set the variable names
      call Set_Eos_Use(EoS)             ! update the use flags

      return
   End Subroutine Init_EoS_Transition

   !!----
   !!---- SUBROUTINE INIT_EOS_TYPE
   !!----
   !!---- Initialize EoS_Type setting all parameters to sensible values
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Init_EoS_Type(Eos, CLin, IThermal, ITransition, Ishear, Icross)
      !---- Arguments ----!
      type (EoS_Type),            intent(out)    :: Eos          ! EoS Type
      character(len=*), optional, intent(in)     :: CLin         ! Character variable to indicate linear EoS or not
      integer,          optional, intent(in)     :: IThermal     ! integer to indicate ithermal type
      integer,          optional, intent(in)     :: ITransition  ! integer to indicate transition type
      integer,          optional, intent(in)     :: IShear       ! integer to indicate shear type
      integer,          optional, intent(in)     :: ICross       ! integer to indicate cross-terms type

      !> test for optional argument for linear
      EoS%Linear  =.false.
      if (present(clin)) then
         if (index(U_case(clin(1:3)),'LIN') > 0) EoS%linear=.true.
      end if

      EoS%Title   =' '
      EoS%IModel  =0
      EoS%IOrder  =3
      EoS%IAngle  =0

      EoS%ParName=' '
      EoS%comment=' '
      EoS%doc=' '
      EoS%savedate=' '
      EoS%lineardir=' '
      EoS%pscale_name=' '
      EoS%vscale_name=' '
      call Set_Eos_Names(EoS)         ! also sets the print factors for the pressure part

      EoS%PRef     = 0.0_cp
      EoS%Density0 = 0.0_cp

      EoS%Iuse     =0
      EoS%Iuse(1:4)=1                 ! Vo, Ko, Kpp

      EoS%params   = 0.0_cp
      EoS%esd      = 0.0_cp

      EoS%params(1)= 1.0_cp
      EoS%params(2)=10.0_cp
      EoS%params(3)= 4.0_cp

      EoS%params(5)= 1.0_cp          ! Z for APL, set non-zero for safety. params(5) not used by any other EoS
      
      
      EoS%X        = 0.0_cp
      EoS%stoich   = 1.0_cp

      EoS%WChi2    = 0.0_cp
      EoS%DelPMax  = 0.0_cp
      EoS%IWt      = 0

      EoS%IRef     = 0
      EoS%factor   = 1.0_cp
      EoS%LastShift= 0.0_cp
      EoS%VCV      = 0.0_cp

      !> Test for optional argument for thermal
      EoS%ITherm  = 0
      EoS%Tref    = 298.0_cp         ! Sensible default
      if (present(ithermal) .and. ithermal > -2 )then
         EoS%Itherm = ithermal
         call Init_Eos_Thermal(EoS)              ! set up default values, names for specific thermal eqn
      end if

      !> Test for optional argument for transition
      EoS%ITran  =0
      if (present(itransition) .and. itransition  > -1 )then
         EoS%ITran = itransition
         call Init_Eos_Transition(EoS)              ! set up default values, names for specific transition
      end if

      !> Test for optional argument for cross terms
      EoS%ICross  =0
      if (present(icross) .and. icross  > -1 )then
         EoS%Icross = icross
         call Init_Eos_Cross(EoS)              ! set up default values, names for specific crosstrem model
      end if

      !> Test for optional argument for shear model
      EoS%Ishear  =0
      if (present(ishear) .and. ishear  > -1 )then
         EoS%Ishear = ishear
      end if

      call Init_Eos_Shear(EoS)              ! set up default values, names for specific shear model

      call Init_EoS_Osc(EoS,3)

      call Init_EoS_Groupscales(EoS)

      call Init_EoS_Angles(Eos)

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
      Warn_EoS=.false.
      Warn_Eos_Mess=" "

      return
   End Subroutine Init_Err_EoS



   !!----
   !!---- SUBROUTINE PHYSICAL_CHECK
   !!----
   !!---- Check if the parameters have physical sense
   !!---- New routine with new logic
   !!---- Returns on first error
   !!----
   !!---- Date: 19/07/2018
   !!
   Subroutine Physical_Check(EoS, Pin, Tin, Vin)
      !---- Arguments ----!
      type(Eos_Type),        intent(in) :: EoS  ! EoS object
      real(kind=cp),optional,intent(in) :: pin  ! Pressure
      real(kind=cp),optional,intent(in) :: vin  ! volume
      real(kind=cp),optional,intent(in) :: tin  ! Temperature

      !---- Local variables ----!
      integer             :: n
      character(len=100)  :: car
      real(kind=cp)       :: tlimit,pinf,p,v,t,vmin !,pmin
      type(eos_type)      :: e,eiso
      logical             :: vpresent

      !>local copies
      E=EoS
      T=e%tref

      !> check PVT present
      n=0
      if (present(Tin))then
         T=Tin
         n=n+1
      end if
      P=0._cp
      if (present(Pin))then
         P=Pin
         n=n+1
      end if

      !> Volume : This is needed for most tests of most EoS
      V=0._cp
      Vpresent=.false.
      if (present(Vin))then
         if (Vin < 0._cp)then
            err_eos=.true.
            err_eos_mess='Volume is negative'
            return
         end if
         V=Vin
         n=n+1
         Vpresent=.true.
      end if
      if (n == 0)return      !no arguments
      if (e%imodel > 0 .and. e%itherm > 0 .and. n < 2)return   ! not enough arguments for PT eos

      !> Positive T
      if (t < 0.0_cp) then
         err_eos=.true.
         err_eos_mess='T is less than zero K'
         return
      end if

      !> Now check for valid parameters at reference
      call EoSParams_Check(E)
      if (err_eos)return

      !> Now check pthermal and isothermal seperately: Pthermal is first
      if (e%pthermaleos)then
          if (e%params(3) > 0._cp)then  ! K limit does not occur if Kp or Mp negative
             ! FIRST find the V at which K=K0/2 at Tref, WITHOUT using pressure
             eiso=e
             eiso%pthermaleos=.false.
             eiso%itherm=0
             vmin=get_volume_K(eiso%params(2)/2.0_cp,eiso%tref,eiso)
             if (vpresent)then
                if (v > vmin)then
                   err_eos=.true.
                   err_eos_mess='Thermal pressure EoS not valid at this V and T: the V is too big so the compressional part of the EoS at Tref is not valid'
                   return
                end if
                if (k_cal(v,t,e) < tiny(0._cp))then
                   err_eos=.true.
                   err_eos_mess='Thermal pressure EoS not valid at this V and T: the K is negative (maybe because of q large?)'
                   return
                end if

             else
                ! No volume input. So calculate the isochor Pressure of Vmin at the input T, and compare to input P
                if (get_k(p,t,e) < tiny(0._cp))then
                   err_eos=.true.
                   err_eos_mess='Thermal pressure EoS not valid at this P and T: the K is negative (maybe because of q large?)'
                   return
                end if
                if (get_volume(p,t,e) > Vmin)then
                   err_eos=.true.
                   err_eos_mess='Thermal pressure EoS not valid at this P and T: the V is too big so the compressional part of the EoS at Tref is not valid'
                   return
                end if
             end if
          end if

      else  !isothermal or no thermal: check thermal part first for T being valid
         !> Check validity of normal-type thermal model: only needs T
         select case(e%itherm)
            case (2)                ! Fei:
               if (e%params(12) > tiny(0.0_cp)) then  ! non-physical V and divergent alpha at low T when alpha2 .ne. 0
                  tlimit=(2.0_cp*e%params(12)/e%params(11))**(1.0_cp/3.0_cp)
                  if (t < tlimit) then
                     err_eos=.true.
                     write(unit=car,fmt='(f5.1)')tlimit
                     car=adjustl(car)
                     err_eos_mess='Fei equation yields non-physical behaviour below T = '//trim(car)//'K'
                     return
                  end if

               else if(e%params(12) < tiny(0.0_cp)) then  ! alpha2 < 0
                  tlimit=sqrt(-1.0_cp*e%params(12)/e%params(10))
                  if (t < tlimit) then
                     err_eos=.true.
                     write(unit=car,fmt='(f5.1)')tlimit
                     car=adjustl(car)
                     err_eos_mess='Fei equation yields non-physical behaviour below T = '//trim(car)//'K'
                     return
                  end if
               end if

            case(3)               ! HP 1998: trap non-physical behaviour at low T
               tlimit=((10.0_cp*e%params(10)+e%params(11))/e%params(10))**2.0_cp
               if (t < tlimit) then
                  err_eos=.true.
                  write(unit=car,fmt='(f5.1)')tlimit
                  car=adjustl(car)
                  err_eos_mess='HP1998 equation yields non-physical behaviour below T = '//trim(car)//'K'
                  return
               end if
         end select

         !> Now check the validity of Eos params at T
         call pveos_check(e,P,V,T,vpresent)
         if (err_eos)then
            err_eos_mess='Compressional EoS not valid at this PV: '//trim(err_eos_mess)
            return
         end if
      end if

      !If got to here, now check that properties at P,T,V valid of Full EoS
      ! because  checks  above are for the PV part and the TV part, without transitions.
      ! all must be valid for the Eos to be valid

      if (e%itherm /=7 .and. e%itherm /=8 .and. .not. vpresent )then        !only done if V not provided at start
         v=get_volume(p,t,e)
         if (err_eos)then         ! added 22/05/2017
            write(unit=car, fmt='(2f10.1)') p, t
            car=adjustl(car)
            err_eos_mess='Volume cannot be calculated at P,T = '//trim(car)
            return
         end if

         if (v < tiny(0.0) ) then
            err_eos=.true.
            write(unit=car, fmt='(2f10.1)') p, t
            car=adjustl(car)
            err_eos_mess='Volume calculated as zero or negative at P,T = '//trim(car)
            return
         end if
      end if

      if (.not. e%linear .and.  V > tiny(0._cp))then
         if (K_cal(V,T,E,P) < tiny(0._cp))then
            write(unit=car, fmt='(2f10.1)') p, t
            car=adjustl(car)
            err_eos_mess='Bulk modulus calculated as zero or negative at P,T = '//trim(car)
            return
         end if
      end if

      !> Produce warning for curved phase boundaries: Pinflection = a/-2b when Ttr=Tr0+aP+bP^2
      if (e%itran>0 .and. abs(e%params(23)) > tiny(0.0) )then
         pinf=abs(e%params(22)/2.0/e%params(23))
         if (abs(p/pinf -1.0) < 0.1) then
            err_eos=.true.
            err_eos_mess='P in region of boundary inflection P: PVT calculations may be inaccurate or wrong'
         end if
      end if

      return
   End Subroutine Physical_Check

   !!----
   !!---- SUBROUTINE PVEOS_CHECK
   !!----
   !!---- Checks compressional part of Eos at P,T for validity (normally that K > 0, or K > K0/2)
   !!---- for volume
   !!---- Does not do transition part
   !!----
   !!---- Update: 17/12/2018
   !!
   Subroutine PVEoS_Check(EoS, Pin, Vin, Tin, Vpresent)
      !---- Arguments ----!
      type(Eos_Type),          intent(in) :: Eos  ! EoS object
      real(kind=cp), optional, intent(in) :: pin  ! Pressure
      real(kind=cp), optional, intent(in) :: vin  ! volume
      real(kind=cp), optional, intent(in) :: tin  ! Temperature
      logical,                 intent(in) :: vpresent ! .true. when the Vin is meaningful

      !---- Local variables ----!
      real(kind=cp)       :: p,v,t
      !real(kind=cp),dimension(3)       :: abc
      real(kind=cp)       :: plim,klim,logterm,vv0,k0,kp !kc,bp,step,kprev,Vnew,Vprev
      type(eos_type)      :: e


      if (EoS%linear)return

      !>local copies
      e=EoS
      p=pin
      t=tin
      v=vin

      !set no transitions
      e%itran=0

      !This routine is private and only called from physical_check
      ! therefore if pthermaleos then T will always be Tref. But set it to be safe, and suppress all thermal part
      if (e%pthermaleos)then
         t=e%tref
         e%pthermaleos=.false.
         e%itherm=0
      end if

      k0=Get_K0_T(T,e)              ! Handles thermal pressure case, returns K0 or M0
      if (k0 < 0._cp)then
         err_eos_mess='K is negative at P=0 and this T'
         err_eos=.true.
         return
      end if
      kp=Get_Kp0_T(T,e)

      ! now do further tests dependening on Vpresent
      ! When V is present, calculate K from V,T
      ! And error state when K <  K(P=0,T)/2, except for Murnaghan which is stable to K=0
      if (vpresent)then
         if (v > Get_V0_T(T,E))then
            select case(e%imodel)
               case(1) ! Murngahan: limit is when K=0
                  if (p < -1.0_cp*k0/kp)err_eos=.true.

               case(2,3,4,5,6)   ! BM, Vinet, NS, Tait, APL
                  if (K_cal(V,T,E) < get_K0_T(T,E)/2.0)err_eos=.true.

               case(7)       !Kumar
                  vv0=v/Get_V0_T(T,E)
                  if (vv0*exp((kp+1)*(1-vv0)) < 0.5_cp)err_eos=.true.

            end select
         end if

      else      ! V was not given, but p was
         if (p < 0._cp)then
            select case(e%imodel)
               case(1) ! Murngahan
                  if (p + 1.0_cp*k0/kp < tiny(0.))err_eos=.true.

               case(2,3,4,5,6) ! find V that gives K = K(P=0,T)/2, by iteration
                  Klim=get_K0_T(T,E)/2.0_cp
                  V=get_volume_K(Klim,e%tref,e)
                  plim=get_pressure(V,T,e)
                  if (p < plim)then
                     err_eos=.true.
                  end if

               case(7)       !Kumar
                  logterm=(kp+1)*p/k0 +1
                  if (logterm < tiny(0._cp))then
                     err_eos=.true.
                     return
                  else
                     if ((1-log(logterm)/(kp+1))*logterm < 0.5_cp)err_eos=.true.
                  end if
            end select

         end if

      end if

      if (err_eos)err_eos_mess='K < K0/2'

      return
   End Subroutine PVEoS_Check

   !!----
   !!---- SUBROUTINE READ_EOS_DATAFILE
   !!----
   !!---- General routine to read data for Eos
   !!----
   !!---- Update: 06/12/2018  RJA: reading of Vscale and Pscale from file
   !!
   Subroutine Read_EoS_DataFile(FName, Dat)
      !---- Arguments ----!
      character(len=*),          intent(in)  :: fname   ! File name
      type (eos_data_list_type), intent(out) :: dat     ! data structure

      !---- Local Variables ----!
      character(len=255), dimension(:), allocatable :: flines
      character(len=255)                            :: line
      character(len=5)                              :: car
      character(len=1)                              :: Ts
      character(len=30), dimension(NCOL_DATA_MAX)   :: dire

      integer                                       :: nldata,ndat, npos
      integer                                       :: i,j,kk,nl,nk,nlines,m,idatatype,nlines_datum
      integer                                       :: iv, inum
      integer, dimension(NCOL_DATA_MAX)             :: ivet,iorden

      real(kind=cp), dimension(NCOL_DATA_MAX)       :: vet,vetsd
      real(kind=cp), dimension(NCOL_DATA_MAX)       :: rvet
      logical                                       :: esd_as_num

      !---- Data Files ----!
      character(len=512)                            :: filedat=' '
      integer                                       :: NCDat         ! Total columns for data
      integer, dimension(NCOL_DATA_MAX)             :: IC_Dat        ! Which values are input - local copy

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

      !> PScale
      i=1
      j=nlines
      call Read_Key_Str(flines, i, j, 'PScale', line,'#')
      if (len_trim(line) > 0)then
          line=adjustl(U_case(line))
          j=len_trim(line)
          if(j-i > 15)j=15
          dat%Pscale_name=line(1:j)
      endif

      !> VScale
      i=1
      j=nlines
      call Read_Key_Str(flines, i, j, 'VScale', line,'#')
      if (len_trim(line) > 0)then
          line=adjustl(U_case(line))
          j=len_trim(line)
          if(j-i > 15)j=15
          dat%Vscale_name=line(1:j)
      endif

      !> LScale
      i=1
      j=nlines
      call Read_Key_Str(flines, i, j, 'LScale', line,'#')
      if (len_trim(line) > 0)then
          line=adjustl(U_case(line))
          j=len_trim(line)
          if(j-i > 15)j=15
          dat%Lscale_name=line(1:j)
      endif

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
         dat%eosd(ndat)%igrp(1)=1       ! Group 1
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
   !!---- Date: 12/10/2015
   !!
   Subroutine Read_Eos_File(FName, Eos)
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
   !!--++ General routine to read one Eos from a file
   !!--++
   !!--++ Date: 12/10/2015
   !!--++
   !!
   Subroutine Read_Eos_In(Flines,Eos)
      !---- Arguments ----!
      character(len=*),dimension(:),intent(in)   :: flines
      type(Eos_type),               intent(out)  :: eos

      !---- Local Variables ----!
      integer       :: nl, imax,ierr,idoc,nlines,i,c,j,jc,kc
      character(len=255)                            :: text
      character(len=10)                             :: forma
      real(kind=cp)                                 :: val
      character(len=255)                            :: wtext    !local warning message, transferred to Warn_Eos_Mess on exit
      logical                                       :: warn     !local flag: needed to prevent potential over-writing by other routines

      !> initialisation
      call init_eos_type(eos)

      nl=0
      imax=0
      ierr=0
      idoc=0                      ! local counter
      nlines=size(flines)
      warn=.false.
      wtext=' '

      do
         nl=nl+1
         if (nl > nlines) exit
         text=adjustl(flines(nl))
         if (len_trim(text) <=0) cycle   ! blank line
         if (text(1:1) == '!') cycle     ! comment line
         if (ierr /=0) exit

         c=index(text,'=')+1             !  1 place after the =
         text(1:c-1)=U_case(text(1:c-1)) ! set keyword in caps

         if (index(text,'TITLE') /= 0) then
            eos%title=trim(text(c:))

         else if(index(text,'SAVEDATE') /= 0)then
            eos%savedate=trim(text(c:))

         else if(index(text,'COMMENT') /= 0)then
            idoc=idoc+1
            if(idoc <= size(eos%doc))eos%doc(idoc)=trim(text(c:))

         else if(index(text,'SYSTEM') /= 0)then
             eos%system=trim(adjustl(text(c:)))

         else if(index(text,'MODEL') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%imodel
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Model number"

         else if(index(text,'ORDER') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%iorder
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Order number"

         else if(index(text,'THERMAL') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%itherm
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Thermal model"

         else if(index(text,'CROSS') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%icross
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Cross-terms model"

         else if(index(text,'TRANS') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%itran
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Transition model"

         else if(index(text,'SHEAR') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%ishear
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Shear model"


         else if(index(text,'OSC2') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%iosc(1)
            if (ierr /=0) Err_EoS_Mess="Error reading the Oscillator2 model"

         else if(index(text,'OSC3') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%iosc(2)
            if (ierr /=0) Err_EoS_Mess="Error reading the Oscillator3 model"

         else if(index(text,'PSCALE') /= 0)then
            read(text(c:),'(a)',iostat=ierr)eos%pscale_name
            if (ierr /=0) Err_EoS_Mess="Error reading the Pressure Scale info"

         else if(index(text,'VSCALE') /= 0)then
            read(text(c:),'(a)',iostat=ierr)eos%vscale_name
            if (ierr /=0) Err_EoS_Mess="Error reading the Volume Scale info"

         else if(index(text,'TYPE') /= 0)then
            if(index(U_case(text),'LINEAR') /= 0) eos%linear=.true.

         else if(index(text,'DIRECTION') /= 0)then
            read(text(c:),'(a)',iostat=ierr)eos%LinearDir
            if (ierr /=0) Err_EoS_Mess="Error reading the direction info"

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

         else if(index(text,'ANGLES') /= 0)then
            read(text(c:),'(i5)',iostat=ierr)eos%iangle
            if (ierr /=0) Err_EoS_Mess="Error reading the Angle Model number"

         !> Parameter values: one at a time
         else if(index(text,'PARAM') /= 0)then
            if (index(U_case(text(c:)),'INF')> 0)then
               read(text(c:),'(i2)',iostat=ierr)i
               val=huge(0._cp)
            else
               read(text(c:),'(i2,f12.6)',iostat=ierr)i,val
               if (i > 50 .and. i < 60)then  !read scale factor name
                  jc=index(U_case(text),',')+1
                  kc=index(U_case(text),')')-1
                  if (jc > 1 .and. kc > jc)eos%comment(i)=trim(adjustl(text(jc:kc)))
               end if
            end if
            if (ierr /=0) Err_EoS_Mess="Error reading the EoS Parameters"

            if (i > 0 .and. i <= N_EOSPAR) then
               eos%params(i)=val
               imax=max(imax,i)   ! use this for vcv reading: allows reading of old files when n_eospar is increased
            end if

         else if(index(text,'VARIANCE') /= 0)then
            ierr=0
            forma="(   e12.5)"
            write(unit=forma(2:4),fmt="(i3)") imax
            do i=1,imax
               read(unit=flines(nl+i),fmt=forma,iostat=ierr) eos%vcv(i,1:imax)
               if (ierr /= 0)exit
            end do
            if (ierr /= 0)then
               Err_EoS_Mess="Error reading the EoS Variance information"
               exit
            end if
         end if
      end do

      !> Error during reading
      if (ierr /=0) then
         err_eos=.true.
         return
      end if

      if (eos%iangle > 0 .and. eos%imodel+eos%itherm+eos%itran+eos%icross+eos%ishear /=0)then
         Err_EoS_Mess="EoS file contains both angle and eos models: not allowed"
         err_eos=.true.
         return
      end if

      !> Do stuff to allow for old files made prior to Nov 2016 not having icross:
      if (eos%icross == 0 .and. abs(eos%params(8)) > 0.000001_cp) eos%icross=1
          !>0.000001 is the smallest non-zero number in the eos file format
          !>Old files cannot be icross=2, only =1

      !> Trap move of cross terms from params(5:6) to params(8:9) in Nov 2019
      !> old files will return icross > 0 but params(5) /= 0, while params(8:9) = 0
      !> new files wil have params(8) and/or (9) non zero

      if (eos%icross > 0 .and. abs(eos%params(8)) < 0.000001_cp .and. eos%imodel /=6)then
         eos%params(8:9) = eos%params(5:6)
         eos%params(5:6) = 0._cp
         eos%vcv(8:9,1:imax)=eos%vcv(5:6,1:imax)
         eos%vcv(1:imax,8:9)=eos%vcv(1:imax,5:6)
         eos%vcv(5:6,1:imax)=0._cp
         eos%vcv(1:imax,5:6)=0._cp
         warn=.true.
         Wtext="EoS file in old format for cross terms: check parameter values carefully"
      end if

      !> Trap move of Z for APL from params(4) to (5): do it after icross so no conflict with cross terms
      if (eos%imodel == 6 .and. abs(eos%params(5)) < 0.000001_cp)then
         eos%params(5) = eos%params(4)
         eos%params(4) = 0._cp     ! Kpp0 will be reset by implied values
         eos%vcv(5,1:imax)=eos%vcv(4,1:imax)
         eos%vcv(1:imax,5)=eos%vcv(1:imax,4)
         eos%vcv(4,1:imax)=0._cp
         eos%vcv(1:imax,4)=0._cp
         warn=.true.
         Wtext=trim(Wtext)//"  EoS file in old format for APL EoS: check parameter values carefully"
      end if

      !> Move angle polynomial values from params(1:30) and clear params(1:30)
      if (eos%iangle > 0)then
         do i=1,3
            j=10*i-9        ! j=1,11,21
            eos%angpoly(i,0,1)=eos%params(j)
            eos%angpoly(i,1,1:3)=eos%params(j+1:j+3)      !P terms
            eos%angpoly(i,2,1:3)=eos%params(j+4:j+6)      !T terms
            eos%angpoly(i,3,1:3)=eos%params(j+7:j+9)      !PT terms
         end do
         eos%params(1:30)=0._cp
      end if

      !>Trap Natom=0 in HP thermal pressure model
      if(eos%itherm == 6 .and. eos%params(13) == 0)then
         warn=.true.
         Wtext=trim(Wtext)//"  Natom was read as zero. PVT will be correct, but not heat capacities or calculated alpha"
      endif

      !>Warn if extra oscillators
      if(sum(eos%iosc) > 0)then
         warn=.true.
         Wtext=trim(Wtext)//"  Extra oscillators in eos file. These are not supported in this version"
      endif





      !> Now finish setting the other eos components
      call set_eos_names(eos)
      call set_thermal_names(eos)
      call Set_Transition_Names(eos)
      call Set_Shear_Names(eos)
      call Set_Cross_Names(eos)
      call Set_Osc_Names(eos)
      call Set_EoS_Use(eos)
      call set_eos_factors(eos)           ! sets the eos factors without resetting param values

      eos%params=eos%params/eos%factor    ! rescale the values
      do i=1,n_eospar                     ! rescale the vcv matrix
         do j=1,n_eospar
            eos%vcv(i,j)=eos%vcv(i,j)/eos%factor(i)/eos%factor(j)
         end do
         eos%esd(i)=sqrt(eos%vcv(i,i))   ! set the esd's from the vcv
      end do

      call set_eos_implied_values(eos)           ! set implied values. Calls from this routine may clear Warn_Eos_Mess
      if(warn)then
          Warn_Eos_Mess=trim(Warn_Eos_Mess)//' '//trim(wtext)
          Warn_eos=.true.
      endif


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
   Subroutine Read_Multiple_Eos_File(FName, Eoslist)
      !---- Arguments ----!
      character(len=*),     intent(in)  :: fname      ! File name
      type (EoS_List_Type), intent(out) :: Eoslist    ! EoS list object

      !---- Variables ----!
      type (EoS_Type)                                :: Eos    ! EoS  object
      character(len=1024), dimension(:), allocatable :: flines
      character(len=512)                             :: filedat

      integer                                        :: nlines,neos,i
      integer,dimension(10)                          :: istart


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
         if (index(U_case(flines(i)),'EOSFIT PARAMETER FILE') > 0)then      ! Use this because it is generated by cfml_eos and not the header
            neos=neos+1
            istart(neos)=i          ! line number of TITLE
         end if
      end do
      istart(neos+1)=nlines+1             ! last line number+1

      !> Set default pointers for one eos if this line not detected in file:
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

   !!----
   !!---- SUBROUTINE SET_CELL_TYPES
   !!----
   !!---- Subroutine to set info, eosc and flags inside eos_cell_type
   !!----
   !!---- Date: 25/09/2020
   !!
   Subroutine Set_Cell_Types(Cell_Eos)
      !---- Arguments ----!
      type(eos_cell_type), intent(in out) :: cell_eos

      !---- Local Variables ----!
      integer :: i

      !> Init
      cell_eos%cout='N'
      cell_eos%loaded=0        ! clear and reset

      select case(U_case(cell_eos%system(1:4)))
         case('TRIC')
            if (cell_eos%eosang%iangle >  0)then
               cell_eos%n=3
               cell_eos%inputlist='(a,b,c,V,Ang)'
            else
               cell_eos%n=6
               cell_eos%inputlist='(a,b,c,d100,d010,d001,V)'
            end if

         case('MONO')
            cell_eos%n=3
            cell_eos%inputlist='(a,b,c,V)'
            if (cell_eos%eosang%iangle >  0)cell_eos%inputlist='(a,b,c,V,Ang)'

         case('ISOT')
            cell_eos%n= 0
            cell_eos%inputlist='V'

         case default
            cell_eos%n=3
            cell_eos%inputlist='(a,b,c,V)'
      end select

      !> set all eos%system to cell_eos%system
      if (len_trim(cell_eos%system) > 0)then
         cell_eos%eos(0:6)%system=cell_eos%system
         cell_eos%eosang%system=cell_eos%system
      end if

      !> initial testing of eos to see if present
      do i=0,cell_eos%n
         if (cell_eos%eos(i)%imodel /= 0 .or.  cell_eos%eos(i)%itherm /= 0)cell_eos%loaded(i)=1
      end do

      !>now check to see if we can calculate missing
      call init_err_eos()
      select case(U_case(cell_eos%system(1:4)))
         case('TRIC')
            !requires all 7 eos or only 3 of 4 if angle poly
            if (cell_eos%eosang%iangle > 0)then        !angle poly
               !angle poly should be loaded: only need 3 of 4 eos Vabc
               if (sum(cell_eos%loaded(0:3)) == 3)then
                  do i=0,3
                     if (cell_eos%loaded(i) == 0)then
                        cell_eos%loaded(i)=3
                        exit
                     end if
                  end do
               end if
            end if

         case('MONO')
            !Set the unique axis flags
            if (index('ABC',U_case(cell_eos%unique_label)) > 0)then
               cell_eos%unique=index('ABC',U_case(cell_eos%unique_label))
            else
               i=index(cell_eos%system,'-')-1
               if (i > 0)cell_eos%unique=index('ABC',U_case(cell_eos%system(i:i)))
            end if

            ! no symmetry equivs, so loaded is either 1 or 0, until here. V,a,b,c are all required unless angle poly set
            if (cell_eos%eosang%iangle > 0)then        !angle poly
               !angle poly should be loaded: only need 3 of 4 eos Vabc
               if (sum(cell_eos%loaded(0:3)) == 3)then
                  do i=0,3
                     if (cell_eos%loaded(i) == 0)then
                        cell_eos%loaded(i)=3
                        exit
                     end if
                  end do
               end if
            end if

         case('ORTH')
            ! no symmetry equivs, so loaded is either 1 or 0, until here
            if (sum(cell_eos%loaded(0:3)) == 3)then
               do i=0,3
                  if (cell_eos%loaded(i) == 0)then
                     cell_eos%loaded(i)=3
                     exit
                  end if
               end do
            end if

         case('TETR','TRIG','HEXA')
            if (cell_eos%loaded(1) == 0 .and. cell_eos%loaded(2) == 1)then     ! move b to a
               cell_eos%eos(1)=cell_eos%eos(2)
               cell_eos%loaded(1)=1
            end if

            if (cell_eos%loaded(1) == 1)then                          ! make b = a
               cell_eos%loaded(2)=2
               cell_eos%eos(2)=cell_eos%eos(1)
            end if

            if (cell_eos%loaded(0) == 0 .and. cell_eos%loaded(1) == 1 .and. cell_eos%loaded(3) == 1) cell_eos%loaded(0)=3   !V from a and c
            if (cell_eos%loaded(1) == 0 .and. cell_eos%loaded(0) == 1 .and. cell_eos%loaded(3) == 1) cell_eos%loaded(1)=3   ! a from V and c
            if (cell_eos%loaded(3) == 0 .and. cell_eos%loaded(0) == 1 .and. cell_eos%loaded(1) == 1) cell_eos%loaded(3)=3   ! C from a and V

         case('CUBI','ISOT') !
            !> first move b or c loaded to a
            if (cell_eos%loaded(2) == 1)then
               cell_eos%eos(1)=cell_eos%eos(2)
               cell_eos%loaded(1)=1

            else if(cell_eos%loaded(3) == 1)then
               cell_eos%eos(1)=cell_eos%eos(3)
               cell_eos%loaded(1)=1
            end if

            !> now set dependencies
            if (cell_eos%loaded(0) == 1)then              ! V loaded
               cell_eos%loaded(1:3)=3                   ! calc a,b,c from V

            else if(cell_eos%loaded(1) == 1)then          ! a loaded
               cell_eos%loaded(0)=3
               cell_eos%loaded(2:3)=2
               cell_eos%eos(2)=cell_eos%eos(1)
               cell_eos%eos(3)=cell_eos%eos(1)
            end if

      end select

      !> Update flags and eosc
      call init_eos_type(cell_eos%eosc)      ! clears eosc
      do i=0,cell_eos%n
         if (cell_eos%loaded(i) == 1)then
            cell_eos%eosc=cell_eos%eos(i)
            exit
         end if
      end do
      cell_eos%eosc%imodel=1                 ! dummies for i/o control
      cell_eos%eosc%itherm=1
      do i=0,cell_eos%n
         if (cell_eos%loaded(i) == 1)then
            if (cell_eos%eos(i)%imodel /= 0)then
               cell_eos%cout(i,1)='Y'

            else
               cell_eos%eosc%imodel=0
            end if

            if (cell_eos%eos(i)%itherm /= 0)then
               cell_eos%cout(i,2)='Y'

            else
               cell_eos%eosc%itherm=0
            end if

            if (cell_eos%cout(i,1) == 'Y' .and. cell_eos%cout(i,2) == 'Y')cell_eos%cout(i,3)='Y'
         end if
      end do

      return
   End Subroutine Set_Cell_Types

   !!--++
   !!--++ SUBROUTINE SET_CROSS_NAMES
   !!--++
   !!--++ Set the character variables in eos_type data structures for Cross-terms
   !!--++
   !!--++ Date: 11/11/2016
   !!
   Subroutine Set_Cross_Names(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS  ! EoS object

      !---- Local Variables ----!
      character(len=50) :: ptext

      !> Check for valid model number. If not valid, set zero
      if (EoS%icross < 0 .or. EoS%icross > N_CROSS_MODELS) EoS%icross=0

      EoS%cmodel=crossmodel_names(EoS%icross)

      select case(EoS%icross)
         case (0)
            EoS%parname(8:9) = ' '

         case (1)
            if (len_trim(EoS%vscale_name) > 0)then
               ptext='units are '//trim(EoS%pscale_name)//'/K'
            else
               ptext='units are P units/K'
            end if
            if (EoS%linear)then
               EoS%parname(8) = 'dM/dT'
               EoS%comment(8) = 'dM/dT '//trim(ptext)
            else
               EoS%parname(8) = 'dK/dT'//trim(ptext)
               EoS%comment(8) = 'dK/dT '//trim(ptext)
            end if

         case (2)
            EoS%parname(8) = 'delta'
            EoS%comment(8) = 'Anderson delta_T, without units'
            EoS%parname(9) = 'delPr'
            EoS%comment(9) = 'delta_prime for Kprime power law, without units'
      end select

      return
   End Subroutine Set_Cross_Names

   !!--++
   !!--++ SUBROUTINE SET_EOS_FACTORS
   !!--++
   !!--++ Initialize the EoS Factors without change in the parameters values
   !!--++
   !!--++ Date: 17/02/2015
   !!
   Subroutine Set_EoS_Factors(Eos)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eos

      !> Init
      EoS%factor=1.0
      EoS%alphafactor=1.0E5_cp

      select case(EoS%itherm)
         case (-1)       ! pvt table
            EoS%factor(10)  = 1.0E5_cp                ! factor to multiply alpha values on printing

         case (0)
            EoS%factor(10:19) = 1.0_cp

         case (1)
            EoS%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            EoS%factor(11)  = 1.0E8_cp

         case (2)
            EoS%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            EoS%factor(11)  = 1.0E8_cp
            EoS%factor(12)  = 1.0_cp

         case (3)
            EoS%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            EoS%factor(11)  = 1.0E4_cp

         case (4)
            EoS%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            EoS%factor(11)  = 1.0_cp

         case (5)
            EoS%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            EoS%factor(11)  = 1.0_cp

         case (6)
            EoS%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            EoS%factor(11)  = 1.0_cp

      end select

      select case(EoS%itran)
         case (1:3)
            EoS%factor(20:n_eospar) = 1.0_cp
            EoS%factor(24)          = 1.0E3_cp         ! 1000 for aL
            EoS%factor(26)          = 1.0E3_cp         ! 1000 for aH
      end select

      return
   End Subroutine Set_Eos_Factors

   !!----
   !!---- SUBROUTINE SET_EOS_IMPLIED_VALUES
   !!----
   !!---- Fix Kp and Kpp values from Model and Order of EoSpar
   !!---- And other implied values
   !!---- Date: 17/07/2015
   !!
   Subroutine Set_Eos_Implied_Values(Eos)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eos  ! EoS object

      !---- Local Variables ----!
      real(kind=cp),dimension(N_EOSPAR):: ev           ! local copies of room pressure parameter values
      real(kind=cp)                    :: pfg0,c0,c2   ! variables for APL

      !> Local copy
      ev= EoS_to_Vec(eos) !  ev contains volume-like parameters

      select case (EoS%imodel)
         case (1) ! Murnaghan
            return      ! no defaults

         case (2) ! Birch-Murnaghan
            if (EoS%iorder == 2) ev(3)=4.0_cp
            if (EoS%iorder == 2 .or. EoS%iorder == 3) then
               if (abs(ev(2)) > 0.0) ev(4)=-1.0_cp*((ev(3)-4.0_cp)*(ev(3)-3.0_cp)+35.0_cp/9.0_cp)/ev(2)  !for order 2 and 3
            end if

         case (3) ! Vinet
            if (EoS%iorder == 2) ev(3)=1.0_cp
            if (EoS%iorder == 2 .or. EoS%iorder == 3) then
               if (abs(ev(2)) > 0.0) ev(4)=-1.0_cp*((0.5_cp*ev(3))**2+0.5*ev(3)-19.0_cp/36.0_cp)/ev(2) !for order 2 and 3
            end if

         case (4) ! Natural
            if (EoS%iorder == 2) ev(3)=2.0_cp
            if (EoS%iorder == 2 .or. EoS%iorder == 3) then
               if (abs(ev(2)) > 0.0) ev(4)=-1.0_cp*(1.0_cp + (ev(3)-2.0_cp)+(ev(3)-2.0_cp)**2.0_cp)/ev(2) !for order 2 and 3
            end if

         case (5) ! Tait with definitions of order derived from Holland and Powell (2011)
            if (EoS%iorder == 2) ev(3)=4.0_cp
            if (EoS%iorder == 2 .or. EoS%iorder == 3)then
               if (abs(ev(2)) > 0.0)ev(4)=-1.0_cp*ev(3)/ev(2)
            end if

         case (6) !APL
            pFG0=AFERMIGAS*(ev(5)/ev(1))**1.66666667_cp
            c0=-1.0_cp*log(3.0_cp*ev(2)/pFG0)           ! assumes V in A^3
            if (EoS%iorder == 2)ev(3)=3.0_cp+2.0_cp*c0/3.0_cp

            if (EoS%iorder < 4)then
               c2=1.5_cp*(ev(3)-3.0_cp)-c0
               ev(4)=(20._cp + 12._cp*c0 + c0*c0 + 2.0_cp*c2*(9.0_cp+c0) + 4.0_cp*c2*c2)
               ev(4)=-1.0_cp*ev(4)/9._cp/ev(2)
            end if

         case(7) !Kumar
            if (EoS%iorder == 2)ev(3)=4.0_cp
         end select

         !> Thermal models
         select case(eos%itherm)
         case(6)        ! Holland-Powell thermal pressure
             eos%params(18)=eos%params(1)*eos%params(10)*eos%params(2)/get_cv(eos%pref,eos%tref,eos)
             eos%params(18)=eos%params(18)/EPthermal_factor(Eos)
         end select



      !> Handle linear or volume

         if (.not. EoS%linear) then
         if (EoS%iorder == 2) EoS%params(3)=ev(3)
         if (EoS%iorder == 2 .or. EoS%iorder == 3) EoS%params(4)=ev(4)

      else
         if (EoS%iorder == 2) EoS%params(3)=ev(3)*3.0_cp
         if (EoS%iorder == 2 .or. EoS%iorder == 3) EoS%params(4)=ev(4)*3.0_cp
      end if

      return
   End Subroutine Set_Eos_Implied_Values

   !!----
   !!---- SUBROUTINE SET_EOS_NAMES
   !!----
   !!---- Set the character variables in eos_type data structures
   !!---- to match the flags already set
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Set_Eos_Names(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS   ! EoS object

      !---- Local Variables ----!
      character(len=50),dimension(5)  :: ptext    ! local variable to hold name of pressure scale

      !> Check for valid model number. If not valid, set zero
      if (EoS%imodel < -1 .or. EoS%imodel > N_PRESS_MODELS) EoS%imodel=0

      !> Set the Eos name
      EoS%model=Pmodel_names(EoS%imodel)
      if (EoS%imodel <= 0) return

      !> set the comments for parameters for volume or linear eos
      !> set the pressure scale text first
      ptext=' '
      if (len_trim(EoS%pscale_name) > 0) then
         ptext(2)='units are '//trim(EoS%pscale_name)
         ptext(4)='units are inverse '//trim(EoS%pscale_name)

      else
         ptext(2)='same units as pressure data'
         ptext(4)='inverse pressure units'
      end if

      !> Set the volume/linear scale name
      if (len_trim(EoS%vscale_name) > 0) then
         ptext(1)='units are '//trim(EoS%vscale_name)
      else
         if (.not. EoS%linear) then
            ptext(1)='units as volume data'
         else
            ptext(1)='units as length data'
         end if
      end if

      if (.not. EoS%linear) then
         EoS%ParName(1:4) =(/'V0   ','K0   ','Kp   ','Kpp  '/)

         EoS%comment(1) = 'Reference pressure volume: '//trim(ptext(1))
         EoS%comment(2) = 'Bulk modulus: '//trim(ptext(2))
         EoS%comment(3) = 'dK/dP: dimensionless'
         EoS%comment(4) = 'd2K/dP2: '//trim(ptext(4))
         EoS%LinearDir  = ' '

         select case(EoS%imodel)
            case (6)                !APL
               EoS%ParName(5) = 'Z   '
               EoS%comment(5) = 'N(electrons) in V0'
         end select

      else
         EoS%ParName(1:4) =(/'L0   ','M0   ','Mp   ','Mpp  '/)

         EoS%comment(1) = 'Reference pressure length: '//trim(ptext(1))
         EoS%comment(2) = 'Linear modulus: '//trim(ptext(2))
         EoS%comment(3) = 'dM/dP: dimensionless'
         EoS%comment(4) = 'd2M/dP2: '//trim(ptext(4))

         select case(EoS%imodel)
            case (6)
               EoS%ParName(5) = 'Z   '
               EoS%comment(5) = 'N(electrons) in V0'

         end select

      end if

      !> Thermal models are only set in init_EoS_thermal

      !>Adjustl all scale names so that they can be compared in index
      EoS%pscale_name=trim(adjustl(EoS%pscale_name))
      EoS%vscale_name=trim(adjustl(EoS%vscale_name))

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
   !!----      3      parameter is used and/or should be reported, not settable, cannot be refined
   !!----             (includes implied values)
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Set_EoS_Use(Eos)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eos  ! EoS object

      !---- Local Variables ----!
      integer                     :: i
      integer,dimension(N_EOSPAR) :: useflags       ! local copy: useful for avoiding some resets

      !> Init
      useflags=EoS%iuse
      EoS%iuse=0
      EoS%allowed_orders=.true.

      !> EoS Model
      select case(EoS%imodel)
         case (0)
            EoS%iuse(1)=1                                    ! None eg thermal only
            EoS%allowed_orders=.false.

         case (1,7)
            EoS%iuse(1:3)=1                                  ! Murnaghan, Kumar
            EoS%allowed_orders(2)=.false.
            EoS%allowed_orders(4)=.false.

         case (2,4,5)                                           ! other isothermal EoS
            EoS%iuse(1:EoS%iorder)=1
            if (EoS%iorder < 4) EoS%iuse(EoS%iorder+1:4)=3  !implied values

         case (3)                                           ! Vinet
            EoS%iuse(1:EoS%iorder)=1
            if (EoS%iorder < 4) EoS%iuse(EoS%iorder+1:4)=3  !implied values
            EoS%allowed_orders(4)=.false.

         case (6)                                               !APL only
            EoS%iuse(1:EoS%iorder)=1
            if (EoS%iorder < 4) EoS%iuse(EoS%iorder+1:4)=3  !implied values
            EoS%iuse(5)=2

      end select

      !> Thermal Model
      select case(EoS%itherm)
         case (1)             ! Berman
            EoS%iuse(10:11)=1 ! alpha terms
            EoS%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            EoS%iuse(19)=2    ! Grunesien q power law parameter
            EoS%TRef_fixed   = .false.
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (2)             ! Fei
            EoS%iuse(10:12)=1 ! alpha terms
            EoS%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            EoS%iuse(19)=2    ! Grunesien q power law parameter
            EoS%TRef_fixed   = .false.
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (3)             ! HP 1998
            EoS%iuse(10:11)=1 ! alpha terms
            EoS%iuse(18)=2
            EoS%iuse(19)=2    ! Grunesien q power law parameter
            EoS%TRef_fixed   = .false.
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (4)             ! Holland-Powell thermal expansion, in Kroll form
            if (EoS%imodel ==0) EoS%iuse(3)=2     ! require Kprime_zero but not stable in refinement if no P data (added 27/01/2014 RJA)
            EoS%iuse(10)=1    ! alpha at Tref
            EoS%iuse(11)=1    ! Einstein T set refineable 2 Sept 2020
            EoS%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            EoS%iuse(19)=2    ! Grunesien q power law parameter
            EoS%TRef_fixed   = .false.
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (5)             ! Salje
            EoS%iuse(10:11)=1
            EoS%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            EoS%iuse(19)=2    ! Grunesien q power law parameter
            EoS%TRef_fixed   = .true.
            EoS%pthermaleos  =.false.
            EoS%Osc_allowed  =.false.

         case (6)             ! Thermal pressure in H&P form (no dK/dT): requires a eos model as well
            if (EoS%imodel==0) EoS%iuse(2:4)=2   ! K, Kp refinement is not stable without a pressure model
            EoS%iuse(8:9)=0     ! No dK/dT parameter:
            EoS%iuse(10)=1    ! alpha at Tref
            EoS%iuse(11)=1    ! Einstein T should be reported
            EoS%iuse(13)=2    ! Natoms per formula unit
            EoS%iuse(14)=0    ! Flag for q-compromise: does not appear to user
            EoS%iuse(18)=3    ! Grunesien parameter at Pref,Tref. This is implied by alpha0
            EoS%iuse(19)=0    ! Grunesien q power law parameter: this is a q-comp model
            EoS%TRef_fixed   = .false.
            EoS%pthermaleos  =.true.
            EoS%Osc_allowed  =.true.

         case (7,8)             ! Thermal pressure in MGD or Einstein form,
            if (EoS%imodel==0) EoS%iuse(2:4)=2   ! K, Kp refinement is not stable without a pressure model
            EoS%iuse(8:9)=0     ! No dK/dT parameter:

            EoS%iuse(11)=1    ! Debye T

            EoS%iuse(13)=2    ! Natoms per formula unit
            EoS%iuse(14)=0    ! Flag for q-compromise: does not appear to user
            EoS%iuse(18)=1    ! Grunesien parameter at Pref,Tref
            EoS%iuse(19)=1    ! Grunesien q power law parameter
            if (EoS%params(14) > 0.5_cp)then
               EoS%params(19)=0._cp
               EoS%iuse(19)=0
            end if

            EoS%TRef_fixed   = .false.
            EoS%pthermaleos  = .true.
            EoS%osc_allowed  = .true.

      end select

      !> Phase transition model
      select case(EoS%itran)
         case (1,2)     ! Landau PV or TV
            EoS%iuse(20)=2           !settable, no refine: sense of transition,
            EoS%iuse(21)=1           !settable, allow refine: Ptr or Ttr
            EoS%iuse(24)=1           !settable, allow refine: aL
            EoS%iuse(25)=1           !settable, allow refine: betaL
            EoS%iuse(26)=1           !settable, allow refine: aH
            EoS%iuse(27)=1           !settable, allow refine: betaH

         case (3)     ! Landau PVT
            EoS%iuse(20:22)=2        !settable, no refine: sense of transition, T(tr), dT(Tr)/dP
            EoS%iuse(24)=1           !settable, allow refine: aL,
            EoS%iuse(23)=2           !settable, fixed d2Tr/dP2
            EoS%iuse(25:27)=1        !settable, allow refine:  betaL,aH,betaH
      end select

      !> Shear model
      select case(EoS%ishear)
         case (0)
            EoS%iuse(30)=0            ! No model, G0 set very large

         case (1)
            EoS%iuse(30:34)=2         ! Polynomial model: settable, not refineable
      end select

      !> Cross terms model
      if (EoS%pthermaleos) then
         EoS%icross=0                ! Kill the cross-terms
         EoS%iuse(8:9)=0
         EoS%params(8:9)=0.

      else
         select case(EoS%icross)
            case (0)
               EoS%iuse(8)=0

            case(1)
               EoS%iuse(8)=1

            case(2)
               EoS%iuse(8)=1
               EoS%iuse(9)=2  !settable, no refine: del-prime refinement always unstable
         end select
      end if

      !> Additional oscillators
      EoS%iuse(40:49)=0
      if (EoS%iosc(1) > 0)then
         EoS%iuse(40:43)=1
         if (EoS%params(44) > 0.5_cp)EoS%iuse(43)=0
      end if

      if (EoS%iosc(2) > 0)then
         EoS%iuse(45:48)=1
         if (EoS%params(49) > 0.5_cp)EoS%iuse(48)=0
      end if

      !> Use flags for data scales: leave them unchanged, because we do not have the dataset
      EoS%iuse(50:59)=useflags(50:59)

      !> Set the refine flags to be consistent with the use flags
      do i=1,N_EOSPAR
         if (EoS%iuse(i) /=1) EoS%iref(i)=0
      end do

      return
   End Subroutine Set_Eos_Use


   !!--++
   !!--++ SUBROUTINE SET_OSC_NAMES
   !!--++
   !!--++ Set the character variables in eos_type data structures for 2nd oscillators
   !!--++
   !!--++ Date: 09/03/2020
   !!
   Subroutine Set_Osc_Names(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS  ! EoS object

      !---- Local Variables ----!
      integer :: i,n0,n1

      !> Check for valid model number. If not valid, set zero
      if (EoS%iosc(1) < 0 .or. EoS%iosc(1) > N_OSC_MODELS) EoS%iosc(1)=0
      if (EoS%iosc(2) < 0 .or. EoS%iosc(2) > N_OSC_MODELS) EoS%iosc(2)=0

      !> Set the name
      EoS%oscmodel(1)=oscmodel_names(EoS%iosc(1))
      EoS%oscmodel(2)=oscmodel_names(EoS%iosc(2))

      !> Set all the parameter names again, no harm in doing so

      do i=1,2
          n0=35+5*i
          n1=39+5*i

          select case(EoS%iosc(i))
             case (0)
                EoS%parname(n0:n1) = ' '
                EoS%comment(n0:n1) = ' '

             case (1)
                EoS%parname(n0)  ='mfrac'
                EoS%parname(n0+1)='ThD'
                EoS%parname(n0+2)='gamma'
                EoS%parname(n0+3)='q'
                EoS%parname(n0+4) ='qcomp'
                EoS%comment(n0)='Fraction of modes with this oscillator'
                EoS%comment(n0+1)='Debye Temperature in K'
                EoS%comment(n0+2)='Gruenesien mode gamma for this oscillator'
                EoS%comment(n0+3)='Gruneisen power law in V/V0  for this oscillator'
                EoS%comment(n0+4)='Switch for q-compromise model, +1 for compromise'

             case (2)
                EoS%parname(n0)  ='mfrac'
                EoS%parname(n0+1)='Th_E '
                EoS%parname(n0+2)='gamma'
                EoS%parname(n0+3)='q    '
                EoS%parname(n0+4) ='qcomp'
                EoS%comment(n0)='Fraction of modes with this oscillator'
                EoS%comment(n0+1)='Einstein Temperature in K'
                EoS%comment(n0+2)='Gruenesien mode gamma for this oscillator'
                EoS%comment(n0+3)='Gruneisen power law in V/V0  for this oscillator'
                EoS%comment(n0+4)='Switch for q-compromise model, +1 for compromise'

          end select
      end do

      return
   End Subroutine Set_Osc_Names

   !!--++
   !!--++ SUBROUTINE SET_SHEAR_NAMES
   !!--++
   !!--++ PRIVATE
   !!--++ Set the character variables in eos_type data structures for SHEAR EoS
   !!--++
   !!--++
   !!--++ Date: 11/07/2016
   !!
   Subroutine Set_Shear_Names(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS  ! EoS object

      !---- Local Variables ----!
      integer :: n

      !> Check for valid model number. If not valid, set zero
      if (EoS%ishear < 0 .or. EoS%ishear > N_SHEAR_MODELS) EoS%ishear=0

      !> Set the Eos name
      EoS%smodel=shearmodel_names(EoS%ishear)

      !> Set upper limit to thermal parameter numbers
      n=34
      if (n > N_EOSPAR) n=N_EOSPAR

      select case(EoS%ishear)
         case (0)
            EoS%parname(30:n) = ' '
            EoS%comment(30:n) = ' '

         case (1)
            EoS%parname(30)='G0'
            EoS%parname(31)='dG/dP'
            EoS%parname(32)='d2G/'
            EoS%parname(33)='d3G/'
            EoS%parname(34)='dG/dT'
            EoS%comment(30)='Shear modulus at Pref, Tref, in pressure units'
            EoS%comment(31)='Pressure derivative of shear modulus: no units'
            EoS%comment(32)='2nd Pressure derivative of shear modulus: P^-1'
            EoS%comment(33)='3rd Pressure derivative of shear modulus: P^-2'
            EoS%comment(34)='Temperature derivative of shear modulus'
      end select

      return
   End Subroutine Set_Shear_Names

   !!--++
   !!--++ SUBROUTINE SET_THERMAL_NAMES
   !!--++
   !!--++ Set the character variables in eos_type data structures for thermal EoS
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Set_Thermal_Names(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS  ! EoS object

      !---- Local Variables ----!
      integer :: n

      !> Check for valid model number. If not valid, set zero
      if (EoS%itherm < -1 .or. EoS%itherm > N_THERM_MODELS) EoS%itherm=0

      !> Set the Eos name
      EoS%tmodel=Tmodel_names(EoS%itherm)

      !> Set upper limit to thermal parameter numbers
      n=19
      if (n > n_eospar) n=n_eospar

      !> Set the V0 name and comment here, in case eos is thermal only
      if (.not. EoS%linear) then
         EoS%ParName(1) ='V0   '
         EoS%comment(1) = 'Reference pressure volume:'
         if (len_trim(EoS%vscale_name) > 0) then
            EoS%comment(1) = trim(EoS%comment(1))//' units are '//trim(EoS%vscale_name)
         else
            EoS%comment(1) = trim(EoS%comment(1))//' units as volume data'
         end if
      else          !Linear
         EoS%ParName(1) ='L0   '
         EoS%comment(1) ='Reference pressure length:'
         if (len_trim(EoS%vscale_name) > 0) then
            EoS%comment(1) = trim(EoS%comment(1))//' units are '//trim(EoS%vscale_name)
         else
            EoS%comment(1) = trim(EoS%comment(1))//' units as length data'
         end if
      end if

      select case(EoS%itherm)
         case (0)
            EoS%parname(10:n) = ' '
            EoS%comment(10:n) = ' '

         case (1)
            EoS%parname(10:11) = (/'alph0','alph1'/)
            EoS%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            EoS%comment(11) = 'Linear term thermal expansion x10^8 K^-2'

         case (2)
            EoS%parname(10:12) = (/'alph0','alph1','alph2'/)
            EoS%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            EoS%comment(11) = 'Linear term thermal expansion x10^8 K^-2'
            EoS%comment(12) = '1/T^2 term thermal expansion, K'

         case (3)
            EoS%parname(10:11) = (/'alph0','alph1'/)
            EoS%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            EoS%comment(11) = 'Sqrt term of thermal expansion x10^4 K^-1/2'

         case (4)    ! Kroll needs Kp as well (in case no pressure eos)
            EoS%parname(3) = 'Kp   '
            EoS%comment(3) = 'dK/dP: dimensionless'
            if (EoS%linear)then
               EoS%parname(3) = 'Mp   '
               EoS%comment(3) = 'dM/dP: dimensionless'
            end if
            EoS%parname(10:11) = (/'alph0','Th_E '/)
            EoS%comment(10) = 'Constant of thermal expansion at Tref x10^5 K^-1'
            EoS%comment(11) = 'Einstein temperature in K'

         case (5)
            EoS%parname(10:11) = (/'p1   ','T_sat'/)
            EoS%comment(10) = 'Approx 3x highT thermal expansion x10^5 K^-1'
            EoS%comment(11) = 'Saturation temperature in K'

         case (6)
            EoS%parname(10:11) = (/'alph0','Th_E '/)
            EoS%comment(10) = 'Constant of thermal expansion at Tref x10^5 K^-1'
            EoS%comment(11) = 'Einstein temperature in K'
            EoS%parname(13) = 'Natom'
            EoS%comment(13) = 'Number of atoms per formula unit'
            EoS%parname(14) = 'qcomp'
            EoS%comment(14) = 'Switch for q-compromise model, +1 for compromise'

         case (7)
            EoS%parname(11) = 'ThMGD'
            EoS%comment(11) = 'Debye temperature in K'
            EoS%parname(13) = 'Natom'
            EoS%comment(13) = 'Number of atoms per formula unit'
            EoS%parname(14) = 'qcomp'
            EoS%comment(14) = 'Switch for q-compromise model, +1 for compromise'

         case (8)
            EoS%parname(11) = 'Th_E'
            EoS%comment(11) = 'Einstein temperature in K'
            EoS%parname(13) = 'Natom'
            EoS%comment(13) = 'Number of atoms per formula unit'
            EoS%parname(14) = 'qcomp'
            EoS%comment(14) = 'Switch for q-compromise model, +1 for compromise'
      end select

      !> Common terms for all thermal
      EoS%parname(18) = 'Gamm0'
      EoS%comment(18) = 'Gruneisen parameter at Tref,Pref'
      EoS%parname(19) = 'q    '
      EoS%comment(19) = 'Gruneisen power law in V/V0'

      return
   End Subroutine Set_Thermal_Names

   !!--++
   !!--++ SUBROUTINE SET_TRANSITION_NAMES
   !!--++
   !!--++ Set the character variables in eos_type data structures for Transition EoS
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Set_Transition_Names(EoS)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: EoS  ! EoS object

      !---- Local Variables ----!
      integer :: n

      !> Check for valid model number. If not valid, set zero
      if (EoS%itran < 0 .or. EoS%itran > N_TRANS_MODELS) EoS%itran=0

      !> Set the model name
      EoS%tranmodel=Tranmodel_names(EoS%itran)

      !> Set upper limit to parameter numbers
      n=29
      if (n > n_eospar) n=n_eospar

      select case(EoS%itran)
         case (0)
            EoS%parname(20:n) = ' '
            EoS%comment(20:n) = ' '

         case (1)       ! Landau power law P-V
            EoS%parname(20:27) = (/'High ','Ptr  ','     ','     ','aL   ','betaL','aH   ','betaH'/)
            EoS%comment(20) = 'Indicator = +1 if high P phase is high sym phase'
            EoS%comment(21) = 'Transition pressure'
            EoS%comment(22) = ''
            EoS%comment(23) = ''
            EoS%comment(24) = 'Scaling parameter, low phase x10^3'
            EoS%comment(25) = 'Power law term, low phase'
            EoS%comment(26) = 'Scaling parameter, high phase x10^3'
            EoS%comment(27) = 'Power law term, high phase'

         case (2)       ! Landau power law V-T
            EoS%parname(20:27) = (/'High ','Ttr  ','     ','     ','aL   ','betaL','aH   ','betaH'/)
            EoS%comment(20) = 'Indicator = +1 if high T phase is high sym phase'
            EoS%comment(21) = 'Transition temperature'
            EoS%comment(22) = ''
            EoS%comment(23) = ''
            EoS%comment(24) = 'Scaling parameter, low phase x10^3'
            EoS%comment(25) = 'Power law term, low phase'
            EoS%comment(26) = 'Scaling parameter, high phase x10^3'
            EoS%comment(27) = 'Power law term, high phase'

         case (3)       ! Landau power law PVT
            EoS%parname(20:27) = (/'High ','Ttr  ','Ttr-P','TtrP2','aL   ','betaL','aH   ','betaH'/)
            EoS%comment(20) = 'Indicator = +1 if high T phase is high sym phase'
            EoS%comment(21) = 'Transition temperature at P=0'
            EoS%comment(22) = 'Ttr=Ttr0 + uP + vP^2: P coeff'
            EoS%comment(23) = 'Ttr=Ttr0 + uP + vP^2: P^2 coeff'
            EoS%comment(24) = 'Scaling parameter, low phase x10^3'
            EoS%comment(25) = 'Power law term, low phase'
            EoS%comment(26) = 'Scaling parameter, high phase x10^3'
            EoS%comment(27) = 'Power law term, high phase'
      end select

      return
   End Subroutine Set_Transition_Names

   !!--++
   !!--++ SUBROUTINE SET_VOLUME_FROM_CELL
   !!--++
   !!--++ Sets V and esd(V) from cell parameter data for all data items in dat
   !!--++ If V is present in first data item, no esd is calculated
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Set_Volume_from_Cell(Dat)
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
   !!--++ SUBROUTINE VEC_TO_EOS
   !!--++
   !!--++ Pass values fron a vector to respective EoS parameter
   !!--++
   !!--++ Date: 28/02/2013
   !!
   Subroutine Vec_to_EoS(Vec, Eos)
      !---- Arguments ----!
      real(kind=cp), dimension(:), intent(in)     :: Vec
      type(EoS_Type),              intent(in out) :: Eos

      !> Copy vec to eosparams as 1 to 1 as default
      EoS%params=vec

      if (EoS%linear) then
         EoS%params(1)=vec(1)**(1.0_cp/3.0_cp)

         EoS%params(2:4)=vec(2:4)*3.0_cp          ! recheck cross terms******


         select case(EoS%itherm)                 ! thermal expansion terms
            case (1,2,3)
               EoS%params(10:12)=vec(10:12)/3.0_cp

            case (4,5,6)
               EoS%params(10)=vec(10)/3.0_cp
         end select

         select case(EoS%icross)
            case (1)
               EoS%params(8)=vec(8)*3.0_cp

            case (2)
               EoS%params(8:9)=vec(8:9)
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
   Subroutine Write_Data_Conlev(Xyy, N, Iout)
      !---- Arguments ----!
      real(kind=cp), dimension(:,:), intent(in)   :: xyy  ! output points for plotting
      integer,                       intent(in)   :: n    ! Number of points
      integer, optional,             intent(in)   :: iout ! Iunit output

      !---- Local Variables ----!
      integer :: lun,i

      !> Check
      if (n < 1) return

      !> Unit to print the data
      lun=6
      if (present(iout)) lun=iout

      do i=1,n
         write(lun,'(3x,a,3x,a,3x,a)') trim(rformat(xyy(1,i),10)), trim(rformat(xyy(2,i),10)), &
                                       trim(rformat(xyy(3,i),10))
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
   Subroutine Write_EoS_DataFile(Dat, Lun)
      !---- Arguments ----!
      type (eos_data_list_type), intent(in) :: dat  ! data structure
      integer,                   intent(in) :: lun  ! Unit to write the information

      !---- Variables ----!
      integer              :: ierr,i,j,k
      character(len=128)   :: text

      !> set up the labels in order
      integer, parameter           :: INI=4, IEND=21
      character(len=5),dimension(INI:IEND) :: lab=(/'T    ','sigT ','P    ','sigP ','V    ','sigV ',  &
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

      !> Scales
      if (len_trim(dat%Pscale_name) > 0)then
         write(unit=lun,fmt='(a,a)',iostat=ierr)  'PSCALE ',trim(dat%Pscale_name)
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if

      if (len_trim(dat%Vscale_name) > 0)then
         write(unit=lun,fmt='(a,a)',iostat=ierr)  'VSCALE ',trim(dat%Vscale_name)
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if

      if (len_trim(dat%Lscale_name) > 0)then
         write(unit=lun,fmt='(a,a)',iostat=ierr)  'LSCALE ',trim(dat%Lscale_name)
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if

      !> Datatype: we assume that all data are the same type: responsibility of calling program
      select case(dat%eosd(1)%xtype)
         case(1)
            write(unit=lun,fmt='(a)',iostat=ierr)  'DATATYPE MODULI ISOTHERMAL'
            write(unit=lun,fmt='(a)',iostat=ierr)    '#'

         case(2)
            write(unit=lun,fmt='(a)',iostat=ierr)  'DATATYPE MODULI ADIABATIC'
            write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end select

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
   Subroutine Write_Eos_File(Eos, Lun)
      !---- Arguments ----!
      type (EoS_Type),intent(in)   :: Eos ! EoS object
      integer,        intent(in)   :: lun ! Unit

      !---- Variables ----!
      character(len=12)            :: stext
      character(len=1024)          :: text
      integer                      :: ierr,i,j,k
      real(kind=cp)                :: valp
      real(kind=cp),dimension(10)  :: p

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
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Savedate =',trim(eos%savedate)

      do i=1,size(eos%doc)
         if (len_trim(eos%doc(i)) > 0)then
            write(unit=lun,fmt='(a)') 'Comment ='//trim(eos%doc(i))
         end if
      end do
      write(unit=lun,fmt='(a)',iostat=ierr) ' '

      !> Crystal system
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'System =',trim(eos%system)

      !> eos type
      text=',  ('//trim(eos%model)//')'
      if (eos%imodel == 0) text=',  (none)'
      write(unit=lun,fmt='(a,i3,a)',iostat=ierr) 'Model =',eos%imodel,trim(text)
      write(unit=lun,fmt='(a,i3)',iostat=ierr) 'Order =',eos%iorder

      text=',  ('//trim(eos%tmodel)//')'
      if (eos%itherm == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Thermal =',eos%itherm,trim(text)

      text=',  ('//trim(eos%cmodel)//')'
      if (eos%icross == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Cross =',eos%icross,trim(text)

      text=',  ('//trim(eos%tranmodel)//')'
      if (eos%itran == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Trans =',eos%itran,trim(text)

      text=',  ('//trim(eos%smodel)//')'
      if (eos%ishear == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Shear =',eos%ishear,trim(text)

      write(unit=lun,fmt='(a,i3)',iostat=ierr) 'Angles =',eos%iangle

      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Pref =',eos%pref
      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Tref =',eos%tref
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Pscale =',trim(eos%pscale_name)
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Vscale =',trim(eos%vscale_name)

      do i=1,2
         text=',  ('//trim(eos%oscmodel(i))//')'
         if (eos%iosc(i) == 0)text=',  (none)'
         write(unit=lun,fmt='(a,i1,a,i3,a)',iostat=ierr) 'Osc',i+1,' =',eos%iosc(i),trim(text)
      end do

      if (eos%linear)then
         write(unit=lun,fmt='(a)',iostat=ierr) 'Type = Linear'
         if (len_trim(eos%LinearDir) == 0)then
            write(unit=lun,fmt='(a)') 'Direction = Unknown'
         else
            write(unit=lun,fmt='(a,a)') 'Direction =',trim(adjustl(eos%LinearDir))
         end if

      else if(eos%iangle == 0)then
         write(unit=lun,fmt='(a)',iostat=ierr) 'Type = Volume'

      else
         write(unit=lun,fmt='(a)',iostat=ierr) 'Type = Angles'
      end if

      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Stoich =',eos%stoich
      if (eos%density0 > tiny(0.)) then
         write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Density0 =',eos%density0
      end if

      write(unit=lun,fmt='(a)',iostat=ierr) ' '

      !> Eos parameters
      if (eos%iangle == 0)then
         !> Normal Eos parameters
         do i=1,n_eospar
            valp=eos%params(i)*eos%factor(i)
            if (abs(valp) < 1.0E7_cp)then
               text=rformat(valp,precision(valp)+2)
            else
               write(text,'(''    Inf'')')
            end if

            if (eos%iuse(i) == 0) then
               write(unit=lun,fmt='(a,i2,a12,5a)')'Param =',i,text(1:12)
            else
               write(unit=lun,fmt='(a,i2,a12,5a)')'Param =',i,text(1:12),'     (',eos%parname(i),',  ',&
                     trim(eos%comment(i)),')'
            end if
         end do

      else
         !> angle polynomial to be written into space for params(1:30)
         do i=1,3      ! loop over angles
            p(1)=eos%angpoly(i,0,1)
            p(2:4)=eos%angpoly(i,1,1:3)      !P terms
            p(5:7)=eos%angpoly(i,2,1:3)      !T terms
            p(8:10)=eos%angpoly(i,3,1:3)      !PT terms
            do k=1,10
               j=10*(i-1)+k
               text=rformat(p(k),precision(p(k))+2)
               write(unit=lun,fmt='(a,i2,a12)')'Param =',j,text(1:12)
            end do
         end do

      end if

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
   !!----   NO  header info is printed here.
   !!----
   !!----   Therefore the program header and write_info_eos have to be called first before calling this routine
   !!----   Then write_eoscal_header is called from here
   !!----
   !!----   Change: 06/10/2015 to make write_eoscal_header private, and change name from write_eoscal_file
   !!----   Change: 12/12/2017 created eoscal_text so that errors and values are printed when error state
   !!----   Change: 19/12/2018 added error flag to return to calling program, if warning or error on at least one calc
   !!----   Change: 9/2020 added handling of calculated directions in unit cell
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Write_Eoscal(Pmin, Pmax, Pstep, Tmin, Tmax, Tstep, Tscale_In, Eos, Lun, &
                           Nprint, Eoscal_ERR, Cell_Eos, Axis)
      !---- Arguments ----!
      real(kind=cp),                 intent(in)     ::  pmin, pmax, pstep   ! P to calculate properties
      real(kind=cp),                 intent(in)     ::  tmin,tmax,tstep     ! T to calculate properties
      character(len=*),              intent(in)     ::  tscale_in           ! Name of the Tscale for output, either C or K
                                                                            ! If Pstep or Tstep  < tiny(0.0) then only Pmin (or Tmin) calculated
      type(EoS_Type),                intent(in)     ::  eos                 ! Eos
      integer,                       intent(in)     :: lun                  ! logical unit for printing
      integer,                       intent(out)    :: nprint               ! Number of data printed
      logical,                       intent(out)    :: eoscal_err           ! error flag
      type(EoS_Cell_Type), optional, intent(in out) :: cell_eos
      type(axis_type),     optional, intent(in out) :: axis

      !---- Local variable ----!
      real(kind=cp)           :: p,t      ! The p and T of each calculation
      real(kind=cp)           :: pst,tst  ! local copy of tstep and pstep
      character(len=255)      :: text     ! local text variable
      character(len=1)        :: tscale   ! local name of tscale
      logical                 :: loop_p   ! loop indicator .true. for inner loop of calcs over P
      integer                 :: cellcase ! indicator for type of eos or cell_eos

      !> init
      nprint=0    ! output counter
      eoscal_err=.false.
      cellcase=0

      if (present(axis) .and. present(cell_eos))then
         cell_eos%eosc%linear = .true.       !safety for using write_eoscal_header
         if (axis%ieos == -2)then
            cellcase = 1        !general drection
         else
            cellcase = 3        !third-axis from others
         end if
      end if

      !> Tscale for output: C or K
      if (len_trim(tscale_in) == 0)then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      !> Write file header
      if (cellcase ==0)then
         call write_eoscal_header(eos,lun,tscale)
      else
         call write_eoscal_cell_header(axis,lun,tscale)
      end if

      !> copy Pstep/Tstep
      tst=tstep
      pst=pstep

      !> set up loop control variables
      if (abs(pst) > tiny(0.0))then       ! inner loop over P
         loop_p=.true.
         if (abs(tst) < tiny(0.0))then   ! no outerloop
            tst=10.*max((tmax-tmin),1.0)       ! set tstep big enough to stop loop
         end if
      else
         loop_p=.false.                  ! inner loop over T
         if (abs(pst) < tiny(0.0))then   ! no outerloop
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
            if (cellcase == 0)then
               call init_err_eos()
               call physical_check(eos,Pin=p,Tin=T)
               if (err_eos) then
                  text=trim(rformat(p,6))//'  '//trim(rformat(t,6))//' :   '//trim(err_eos_mess)
                  write(lun,'(a)')trim(text)
                  eoscal_err=.true.

               else
                  call eoscal_text(p,t,Tscale_In,Eos,text)
                  write(lun,'(a)')trim(text)      ! This way we get to see the calculated values even if error in calcs with valid eos
                  if (err_eos)then
                     write(lun,'(a)')'   *****WARNING:   '//trim(err_eos_mess)
                     eoscal_err=.true.
                  end if
               end if

            else
               call init_err_eos()
               call eoscal_text_direction(P,T,Tscale_in,cell_eos,axis,text)
               write(lun,'(a)')trim(text)      ! This way we get to see the calculated values even if error in calcs with valid eos
               if (err_eos)then
                  write(lun,'(a)')'   *****WARNING:   '//trim(err_eos_mess)
                  eoscal_err=.true.
               end if

            end if
            nprint=nprint+1

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
   !!--++ SUBROUTINE WRITE_EOSCAL_CELL_HEADER
   !!--++
   !!--++ Subroutine that print information on iout unit for calculated directions without own eos
   !!--++
   !!--++ Date: 15/09/2020
   !!
   Subroutine Write_Eoscal_Cell_Header(Axis, Lun, Tscale_In)
      !---- Arguments ----!
      type(axis_type),  intent(in) :: axis
      integer,          intent(in) :: lun         ! logical unit for printing
      character(len=*), intent(in) :: tscale_in   ! Scale for Temp

      !---- Local Variables ----!
      character(len=1)     :: tscale
      character(len=255)   :: head     ! local text variable for column headers

      !> alpha scale
      write(lun,'(//)')
      write(lun,'("  Note that values of alpha are multiplied by a factor of 10^5")')

      !> tscale for output: C or K
      if (len_trim(tscale_in) == 0) then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      !> create column header
      if (axis%ieos == 0) then
         write(head,'("   Press   Temp",a1,"   Volume    V/V0T       K    Kprime    dK/dT    alpha")' ) Tscale
      else
         write(head,'("   Press   Temp",a1,"   Length    L/L0T       M    Mprime    dM/dT    alpha")' ) Tscale
      end if

      !> Write header
      write(lun,'(/a)')trim(head)

      return
   End Subroutine Write_Eoscal_Cell_Header

   !!----
   !!---- SUBROUTINE WRITE_EOSCAL_HEADER
   !!----
   !!---- Subroutine that print information on iout unit
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Write_Eoscal_Header(Eos, Lun, Tscale_In)
      !---- Arguments ----!
      type(EoS_Type),   intent(in) :: eos         ! Eos information
      integer,          intent(in) :: lun         ! logical unit for printing
      character(len=*), intent(in) :: tscale_in   ! Scale for Temp

      !---- Local Variables ----!
      character(len=1)     :: tscale
      character(len=255)   :: head     ! local text variable for column headers

      !> alpha scale
      write(lun,'(//)')
      if (eos%itherm /= 0) then
         write(lun,'("  Note that values of alpha are multiplied by a factor of ",f5.1,"x10^5")')eos%alphafactor/1.0E5
      end if

      !> Cp, Cv units
      if (eos%imodel > 0 .and. eos%itherm > 0 .and. .not. eos%linear)then
         if (VscaleMGD(eos) )then
            write(lun,'("  Heat capacities in J/mol/K provided K0 in kbar or GPa")')

         else
            write(lun,'("  Heat capacities are in units that depend on the volume and pressure units")')
         end if
      end if
      write(lun,'(//)')

      if (eos%imodel == 1 .or. eos%imodel == 5 .or. eos%imodel == 6) then
         write(lun,'("  Do not forget: Normalised Pressure and strain not defined for ",a," Eos")')trim(eos%model)

      else if(eos%itherm /= 0)then
         write(lun,'("  Normalised Pressure (NP) and finite strain (f) are defined relative to V at P=0 and same T")')

      else
         write(lun,'("  Normalised Pressure is NP and finite strain is f")')
      end if

      !> tscale for output: C or K
      if (len_trim(tscale_in) == 0) then
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
      if (eos%itran > 0)head=trim(head)//'   spstrain'

      if (eos%density0 > tiny(0.0)) head=trim(head)//'  density  esdden'

      if (eos%pthermaleos .and. eos%itran ==0)head=trim(head)//'  Ptherm'

      if (eos%itherm > 0 .and. abs(eos%params(18)) > tiny(0.0)) then
         if (eos%linear)then
            head=trim(head)//'    Ms    Gamma'
         else
            head=trim(head)//'    Ks    Gamma'
         end if
      end if

      if (eos%itherm == 7)head=trim(head)//' DebyeT'
      if (eos%itherm == 8)head=trim(head)//'    EinT'

      if(eos%itherm /= 0 .and. eos%imodel /= 0 .and. (abs(eos%params(18)) >tiny(0._cp) .or. eos%osc_allowed))then
            if(VscaleMGD(Eos) .and. .not. eos%linear)head=trim(head)//'    Cp      Cv'
      endif


      !> Write header
      write(lun,'(/a)')trim(head)

      return
   End Subroutine Write_Eoscal_Header

   !!--++
   !!--++ Subroutine Write_Info_Angle_Poly
   !!--++
   !!--++ Date: 08/02/2021
   !!
   Subroutine Write_Info_Angle_Poly(EoS, Iang, Iout)
      !---- Arguments ----!
      type(eos_type),    intent(in) :: eos     ! eos type with angles polynomial
      integer, optional, intent(in) :: iang    ! The angle number to print, If missing prints all
      integer, optional, intent(in) :: iout    ! Logical unit

      !---- Local Variables ----!
      integer           :: lun,i,ia,angflag
      character(len=180):: ltext
      character(len=20) :: sign,var

      !> check for angle to print
      angflag=0
      if (present(iang))then
         if (iang < 1 .or. iang > 3)return        !illegal angl number
         angflag=iang
      end if

      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      if (eos%iangle == 0) return              !no polynomial

      do ia=1,3
         if (angflag > 0 .and. ia /= angflag)cycle
         ltext='     '//celllabel(ia+3)//'='//trim(rformat(eos%angpoly(ia,0,1),7))

         if (eos%iangle == 1)then         !polynomial
            do i=1,N_ANGPOLY            ! P terms
               if (abs(eos%angpoly(ia,1,i)) > tiny(0._cp))then
                  sign='+'
                  if (eos%angpoly(ia,1,i) < 0._cp) sign=' '     !rformat sets sign if <0
                  var='P'
                  if (i > 1)write(unit=var,fmt='("P^",i1)')i
                  ltext=trim(ltext)//' '//trim(sign)//trim(adjustl(rformat(eos%angpoly(ia,1,i),8)))//trim(var)
               end if
            end do

            do i=1,N_ANGPOLY
               if (abs(eos%angpoly(ia,2,i)) > tiny(0._cp))then
                  sign='+'
                  if (eos%angpoly(ia,2,i) < 0._cp)sign=' '     !rformat sets sign if <0
                  var='T'
                  if (i > 1)write(unit=var,fmt='("T^",i1)')i
                  ltext=trim(ltext)//' '//trim(sign)//trim(adjustl(rformat(eos%angpoly(ia,2,i),8)))//trim(var)
               end if
            end do

            do i=1,N_ANGPOLY
               if (abs(eos%angpoly(ia,3,i)) > tiny(0._cp))then
                  sign='+'
                  if (eos%angpoly(ia,3,i) < 0._cp)sign=' '     !rformat sets sign if <0
                  var='PT'
                  if (i == 2)var='P^2T'
                  if (i == 3)var='PT^2'
                  ltext=trim(ltext)//' '//trim(sign)//trim(adjustl(rformat(eos%angpoly(ia,3,i),8)))//trim(var)
               end if
            end do
         end if

         !>Write out the out the info
         write(unit=lun,fmt='(a)')trim(ltext)
      end do

      return
   End Subroutine Write_Info_Angle_Poly

   !!----
   !!---- SUBROUTINE Write_Info_Conlev
   !!----
   !!---- Writes out header info specific to a confidence ellipse
   !!----
   !!---- Date: 05/12/2015
   !!
   Subroutine Write_Info_Conlev(Eos, Ix, Iy, Isig, Iout)
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
      write(lun,'(//a,f6.2,a)') "Coordinates for confidence ellipse at confidence level of ",delchi_levels(isig),"%"

      write(lun,'(3x,"X axis as ",a," with variance = ",a)')trim(eos%parname(ix)),trim(rformat(eos%vcv(ix,ix),10))
      write(lun,'(3x,"Y axis as ",a," with variance = ",a)')trim(eos%parname(iy)),trim(rformat(eos%vcv(iy,iy),10))
      write(lun,'(3x,"The covariance of the two variables = ",a)')trim(rformat(eos%vcv(ix,iy),10))

      write(lun,'(//,a,//)') &
                "Data points for plotting: first column is X with two columns with the two Y values at that X"


      return
   End Subroutine Write_Info_Conlev

   !!----
   !!---- SUBROUTINE WRITE_INFO_EOS
   !!----
   !!---- Subroutine that print information on iout unit
   !!----
   !!---- Date: 17/07/2015
   !!
   Subroutine Write_Info_Eos(EoS, Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos   ! EoS object
      integer, optional, intent(in) :: iout  ! Logical unit

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
      write(unit=lun,fmt='(a)') '    Title: '//trim(EoS%title)
      write(unit=lun,fmt='(a)') 'Eos Saved: '//trim(EoS%savedate)

      !> Doc Information
      do i=1,size(EoS%doc)
         if (len_trim(EoS%doc(i)) > 0) then
            write(unit=lun,fmt='(a)') '  Comment: '//trim(EoS%doc(i))
         end if
      end do
      write(unit=lun,fmt='(a)') ' '

      !> Crystal system
      write(unit=lun,fmt='(a)') '   System: '//trim(EoS%system)
      if (EoS%linear) then
         write(unit=lun,fmt='(a)') '    Class: Linear'
         if (len_trim(EoS%LinearDir) == 0)then
            write(unit=lun,fmt='(a)') 'Direction: Unknown'

         else
            write(unit=lun,fmt='(a,a)') 'Direction: ',trim(adjustl(EoS%LinearDir))
         end if

      else if (EoS%iangle == 0)then
         write(unit=lun,fmt='(a)') '    Class: Volume'

      else
         write(unit=lun,fmt='(a)') '    Class: Angles'
      end if

      if (EoS%imodel /= 0 .or. EoS%iangle /= 0) then
         if (len_trim(EoS%Pscale_name) > 0)write(unit=lun,fmt='(a)') '   Pscale: '//trim(EoS%Pscale_name)
      end if

      if (EoS%imodel /= 0) then
         if (len_trim(EoS%Vscale_name) > 0)write(unit=lun,fmt='(a)') '   Vscale: '//trim(EoS%Vscale_name)
         write(unit=lun,fmt='(a,t27,f8.3)') '   Stoichiometry: ',EoS%stoich

         !> Reference Density
         if (EoS%density0 > tiny(0.0)) then
            write(unit=lun,fmt='(a,t27,f8.3)') '   Reference density: ',EoS%density0
         end if

         write(unit=lun,fmt='(a)') ' '
         write(unit=lun,fmt='(a)') '  Compressibility'
         write(unit=lun,fmt='(a)') '-------------------'
         write(unit=lun,fmt='(a)') ' '
         write(unit=lun,fmt='(a)') '    Model: '//trim(EoS%model)
         if (EoS%imodel > 0) then                ! no output if no p eos: all done in write_info_eos_thermal
            write(unit=lun,fmt='(a,i2)') '    Order: ',EoS%iorder

            !> Pressure Parameters
            write(unit=lun,fmt='(a,t27,f8.3)') '   Pressure of reference: ',EoS%pref
            write(unit=lun,fmt='(a)') ' '

            do i=1,5
               if (EoS%iuse(i) /= 0) then
                  call setnum_std(EoS%params(i)*EoS%factor(i),EoS%esd(i)*EoS%factor(i),line)     ! include scaling
                  string=' '
                  select case(EoS%iuse(i))
                     case (2)
                        string=' [FIXED VALUE]'
                     case (3)
                        string=' [IMPLIED VALUE]'
                  end select
                  write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                       trim(EoS%parname(i)),trim(line),trim(EoS%comment(i))//trim(string)
               end if
            end do
            write(unit=lun,fmt='(a)') ' '
         end if
      end if

      !> Thermal EOS
      if (EoS%itherm > 0) call write_info_eos_thermal(EoS,lun)

      !>Extra oscillators
      if (EoS%iosc(1) > 0) call write_info_eos_oscillator(EoS,lun)

      !> Cross terms
      if (EoS%icross > 0) call write_info_eos_cross(EoS,lun)

      !> Transition
      if (EoS%itran > 0) call write_info_eos_transition(EoS,lun)

      !> Shear
      if (EoS%ishear > 0) call write_info_eos_shear(EoS,lun)

      !> Group scales
      if(sum(EoS%iuse(51:59)) > 0) call write_info_eos_groupscales(EoS,lun)

      !> Angle polynomial
      if (EoS%iangle > 0) call Write_Info_Angle_Poly(EoS,iout=lun)

      !> End
      write(unit=lun,fmt='(a)') ' '

      return
   End Subroutine Write_Info_Eos

   !!----
   !!---- SUBROUTINE WRITE_INFO_EOS_CELL
   !!----
   !!---- Subroutine that prints information about the list of eos in cell_eos on iout unit
   !!----
   !!---- Date: 09/09/2020
   !!
   Subroutine Write_Info_Eos_Cell_Type(Cell_Eos, Iout)
      !---- Arguments ----!
      type(eos_cell_type), intent(in out) :: cell_eos   !must be inout to allow the flags to be changed
      integer, optional,   intent(in)     :: iout       ! Logical unit

      !---- Local Variables ----!
      character(len=110) :: ltext
      integer            :: i,lun

      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      !>Update all the flags and pointers
      call set_cell_types(cell_eos)

      write(ltext,'(110(''_''))')
      write(unit=lun,fmt='(a)')ltext
      write(unit=lun,fmt='(a)')'     The crystal system is '//trim(cell_eos%system)

      !>Angle information
      select case(U_case(cell_eos%system(1:4)))
         case('MONO')
            if (cell_eos%unique < 1 .or. cell_eos%unique > 3)then
               write(unit=lun,fmt='(''     Monoclinic unique axis has not been set'')')
            else
               if (cell_eos%eosang%iangle > 0)then
                  write(unit=lun,fmt='(''     Monoclinic angle defined by polynomial:'')')
                  call write_info_angle_poly(cell_eos%eosang,cell_eos%unique,lun)

               else
                  write(unit=lun,fmt='(''     Monoclinic angle defined eos of cell edges and volume:'')')
                  if (cell_eos%obtuse(cell_eos%unique))then
                     write(unit=lun,fmt='(''     Monoclinic angle '',a,'' is set obtuse (>90deg)'')')celllabel(cell_eos%unique+3)

                  else
                     write(unit=lun,fmt='(''     Monoclinic angle '',a,'' is set acute (<90deg)'')')celllabel(cell_eos%unique+3)
                  end if
               end if
            end if

         case('TRIC')
            if (cell_eos%eosang%iangle > 0)then
               write(unit=lun,fmt='(''     Triclinic angles defined by polynomials:'')')
               do i=1,3
                  call write_info_angle_poly(cell_eos%eosang,i,lun)
               end do

            else
               write(unit=lun,fmt='(''     Triclinic angles defined by eos of cell edges, volume, and d-spacings:'')')
               do i=1,3
                  if (cell_eos%obtuse(i))then
                     write(unit=lun,fmt='(''     Triclinic angle '',a,'' is set obtuse (>90deg)'')')celllabel(i+3)

                  else
                     write(unit=lun,fmt='(''     Triclinic angle '',a,'' is set acute  (<90deg)'')')celllabel(i+3)
                  end if
               end do
            end if

         case default
            if (cell_eos%eosang%iangle> 0)then
               write(unit=lun,fmt='(''     Polynomials for angles have been loaded but the crystal system has fixed angles'')')
               write(unit=lun,fmt='(''     The angle polynomials will be ignored unless you change the system to triclinic or monoclinic'')')
            end if

      end select

      write(unit=lun,fmt='(a)')ltext

      write(unit=lun,fmt='(a)')'  The loaded EoS are:                                         PV    VT    PVT  PSCALE      VSCALE'
      do i=0,cell_eos%n
         select case(cell_eos%loaded(i))
            case default
               write(unit=lun,fmt='(2x,a4,a)')axislabel(i),':  No eos loaded'

            case(1)
               write(unit=lun,fmt='(2x,a4,a,a,3x,3(5x,a1),t80,a,t91,a)')axislabel(i),': ',cell_eos%eos(i)%title(1:47), &
                     cell_eos%cout(i,1:3),trim(cell_eos%eos(i)%Pscale_name),trim(' '//cell_eos%eos(i)%Vscale_name)

            case(2)
               write(unit=lun,fmt='(2x,a4,a)')axislabel(i),':  No eos loaded but set by symmetry'

            case(3)
               write(unit=lun,fmt='(2x, a4,a)')axislabel(i),':  No eos loaded but will be calculated from others'

            case(4)
               write(unit=lun,fmt='(2x,a4,a)')axislabel(i),':  Monoclinic unique axis'
         end select
      end do

      write(unit=lun,fmt='(a)')ltext

      write(unit=lun,fmt='(a)') ' '

      return
   End Subroutine Write_Info_Eos_Cell_Type

   !!--++
   !!--++ SUBROUTINE WRITE_INFO_EOS_CROSS
   !!--++
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 11/07/2016
   !!
   Subroutine Write_Info_Eos_Cross(EoS, Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos
      integer, optional, intent(in) :: iout

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Check
      if (Eos%icross ==0) return

      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '      P-T cross-terms'
      write(unit=lun,fmt='(a)') '---------------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(eos%cmodel)

      do i=8,9
         if (eos%iuse(i) /= 0) then
            call setnum_std(eos%params(i)*eos%factor(i),eos%esd(i)*eos%factor(i),line)     ! include scaling
            string=' '
            select case(eos%iuse(i))
               case(2)
                  string=' [FIXED VALUE]'

               case(3)
                  string=' [IMPLIED VALUE]'
            end select
            write(unit=lun,fmt='(3x,a5,'': '',a,T30,'':'',a)') &
                 trim(eos%parname(i)),trim(line),trim(eos%comment(i))//trim(string)
         end if
      end do
      write(unit=lun,fmt='(a)') ' '

      return
   End Subroutine Write_Info_Eos_Cross

   !!--++
   !!--++ SUBROUTINE WRITE_INFO_EOS_GROUPSCALES
   !!--++
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 23/03/2020 RJA
   !!
   Subroutine Write_Info_Eos_GroupScales(EoS, Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos   ! EoS object
      integer, optional, intent(in) :: iout  ! Logical unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Init
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Scale factors for data groups'
      write(unit=lun,fmt='(a)') '---------------------------------'

      do i=51,59
         if (EoS%iuse(i) /= 0) then
            call setnum_std(EoS%params(i)*EoS%factor(i),EoS%esd(i)*EoS%factor(i),line)     ! include scaling
            string=' '
            select case(EoS%iuse(i))
               case(2)
                  string=' [FIXED VALUE]'
            end select
            write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                  trim(EoS%parname(i)),trim(line),trim(EoS%comment(i))//trim(string)
         end if
      end do

      return
   End Subroutine Write_Info_Eos_GroupScales

   !!--++
   !!--++ SUBROUTINE WRITE_INFO_EOS_OSCILLATOR
   !!--++
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 09/03/2020 RJA
   !!
   Subroutine Write_Info_Eos_Oscillator(EoS, Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos     ! EoS object
      integer, optional, intent(in) :: iout    ! Logical unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun,j

      !> Check
      if (EoS%iosc(1) == 0) return  !if no first oscillator, there is no second
      if (.not. EoS%osc_allowed)return

      !> Init
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Additional Oscillators'
      write(unit=lun,fmt='(a)') '--------------------------'

      do j=1,2
         write(unit=lun,fmt='(a)') ' '
         write(unit=lun,fmt='(a,i1,a)') 'Oscillator ',j+1,': '//trim(EoS%oscmodel(j))
         write(unit=lun,fmt='(a)') ' '

         do i=35+5*j,39+5*j
            if (EoS%iuse(i) /= 0) then
               call setnum_std(EoS%params(i)*EoS%factor(i),EoS%esd(i)*EoS%factor(i),line)     ! include scaling
               string=' '
               select case(EoS%iuse(i))
                  case(2)
                     string=' [FIXED VALUE]'

                  case(3)
                     string=' [IMPLIED VALUE]'
               end select
               write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                     trim(EoS%parname(i)),trim(line),trim(EoS%comment(i))//trim(string)
            end if
         end do
      end do

      return
   End Subroutine Write_Info_Eos_Oscillator

   !!--++
   !!--++ SUBROUTINE WRITE_INFO_EOS_SHEAR
   !!--++
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 11/07/2016
   !!
   Subroutine Write_Info_Eos_Shear(EoS, Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos
      integer, optional, intent(in) :: iout

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Check
      if (Eos%ishear ==0) return

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
            select case(eos%iuse(i))
               case(2)
                  string=' [FIXED VALUE]'

               case(3)
                  string=' [IMPLIED VALUE]'
            end select
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
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Write_Info_Eos_Thermal(EoS, Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos   ! EoS object
      integer, optional, intent(in) :: iout  ! Logical unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun,is

      !> Init
      lun=6
      if (present(iout)) lun=iout

      is=10                            ! If pmodel present, it was already reported
      if (EoS%imodel == 0) is=1     ! if no pmodel report all params here


      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Thermal Expansion'
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(EoS%tmodel)

      if (EoS%osc_allowed .and. EoS%params(14) > 0.5_cp)then
          if(eos%itherm == 6)then
            write(unit=lun,fmt='(a)') '           (this is a q-compromise EoS)'
          else
            write(unit=lun,fmt='(a)') '           with q-compromise'
          endif
      endif

      write(unit=lun,fmt='(a)') ' '

      write(unit=lun,fmt='(a,f8.2,a)') '   Temperature of reference: ',EoS%tref,' K'
      write(unit=lun,fmt='(a)') ' '

      do i=is,19
         if (EoS%iuse(i) /= 0) then
            call setnum_std(EoS%params(i)*EoS%factor(i),EoS%esd(i)*EoS%factor(i),line)     ! include scaling
            string=' '
            select case(EoS%iuse(i))
               case(2)
                  string=' [FIXED VALUE]'

               case(3)
                  string=' [IMPLIED VALUE]'
            end select
            write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                  trim(EoS%parname(i)),trim(line),trim(EoS%comment(i))//trim(string)
         end if
      end do

      !> extra stuff if additional oscilators
      if (sum(EoS%iosc) > 0. .and. EoS%osc_allowed)then
         write(unit=lun,fmt='(a)') ' '
         write(unit=lun,fmt='(a,f6.3,a)') '   This oscillator models ',1.0_cp-EoS%params(40)-EoS%params(45),' of the total modes'
      end if

      return
   End Subroutine Write_Info_Eos_Thermal

   !!--++
   !!--++ SUBROUTINE WRITE_INFO_EOS_TRANSITION
   !!--++
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 17/07/2015
   !!
   Subroutine Write_Info_Eos_Transition(EoS, Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: EoS  ! EoS object
      integer, optional, intent(in) :: iout    ! Logical Unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Check
      if(EoS%itran == 0)return

      !> Init
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Phase Transition '
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(EoS%tranmodel)
      write(unit=lun,fmt='(a)') ' '

      do i=20,29
         if (EoS%iuse(i) /= 0) then
            call setnum_std(EoS%params(i)*EoS%factor(i),EoS%esd(i)*EoS%factor(i),line)     ! include scaling
            string=' '
            select case(EoS%iuse(i))
               case (2)
                  string=' [FIXED VALUE]'
               case (3)
                  string=' [IMPLIED VALUE]'
            end select
            write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                  trim(EoS%parname(i)),trim(line),trim(EoS%comment(i))//trim(string)
         end if
      end do

      return
   End Subroutine Write_Info_Eos_Transition

End Module CFML_EoS
