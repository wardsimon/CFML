!!----
!!----
!!----
SubModule (CFML_Eos) EoS_022
   Contains
   
   !!--++
   !!--++ SET_CROSS_NAMES
   !!--++
   !!--++ PRIVATE
   !!--++ Set the character variables in eos_type data structures for Cross-terms
   !!--++
   !!--++ 11/11/2016
   !!
   Module Subroutine Set_Cross_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

      !---- Local Variables ----!
      character(len=50) :: ptext

      !> Check for valid model number. If not valid, set zero
      if (eospar%icross < 0 .or. eospar%icross > N_CROSS_MODELS) eospar%icross=0

      eospar%cmodel=crossmodel_names(eospar%icross)

      select case(eospar%icross)
         case (0)
            eospar%parname(5:6) = ' '

         case (1)
            if (len_trim(eospar%vscale_name) > 0)then
               ptext='units are '//trim(eospar%pscale_name)//'/K'
            else
               ptext='units are P units/K'
            end if
            if (eospar%linear)then
               eospar%parname(5) = 'dM/dT'
               eospar%comment(5) = 'dM/dT '//trim(ptext)
            else
               eospar%parname(5) = 'dK/dT'//trim(ptext)
               eospar%comment(5) = 'dK/dT '//trim(ptext)
            end if

         case (2)
            eospar%parname(5) = 'delta'
            eospar%comment(5) = 'Anderson delta_T, without units'
            eospar%parname(6) = 'delPr'
            eospar%comment(6) = 'delta_prime for Kprime power law, without units'
      end select
   End Subroutine Set_Cross_Names

   !!--++
   !!--++ SET_EOS_FACTORS
   !!--++
   !!--++ PRIVATE
   !!--++ Initialize the EoS Factors without change in the parameters values
   !!--++
   !!--++ 17/02/2015
   !!
   Module Subroutine Set_EoS_Factors(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !> Init
      eospar%factor=1.0
      eospar%alphafactor=1.0E5_cp

      select case(eospar%itherm)
         case (-1)       ! pvt table
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply alpha values on printing

         case (0)
            eospar%factor(10:19) = 1.0_cp

         case (1)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp

         case (2)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp
            eospar%factor(12)  = 1.0_cp

         case (3)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E4_cp

         case (4)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp

         case (5)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp

         case (6)
            eospar%factor(10)  = 1.0E5_cp                ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp

      end select

      select case(eospar%itran)
         case (1:3)
            eospar%factor(20:n_eospar) = 1.0_cp
            eospar%factor(24)          = 1.0E3_cp         ! 1000 for aL
            eospar%factor(26)          = 1.0E3_cp         ! 1000 for aH
      end select
   End Subroutine Set_Eos_Factors

   !!----
   !!---- SET_EOS_NAMES
   !!----    Set the character variables in eos_type data structures
   !!----    to match the flags already set
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine Set_Eos_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar   ! EoS object

      !---- Local Variables ----!
      character(len=50),dimension(5)  :: ptext    ! local variable to hold name of pressure scale

      !> Check for valid model number. If not valid, set zero
      if (eospar%imodel < -1 .or. eospar%imodel > N_PRESS_MODELS) eospar%imodel=0

      !> Set the Eos name
      eospar%model=Pmodel_names(eospar%imodel)
      if (eospar%imodel <= 0) return

      !> set the comments for parameters for volume or linear eos
      !> set the pressure scale text first
      ptext=' '
      if (len_trim(eospar%pscale_name) > 0) then
         ptext(2)='units are '//trim(eospar%pscale_name)
         ptext(4)='units are inverse '//trim(eospar%pscale_name)

      else
         ptext(2)='same units as pressure data'
         ptext(4)='inverse pressure units'
      end if

      !> Set the volume/linear scale name
      if (len_trim(eospar%vscale_name) > 0) then
         ptext(1)='units are '//trim(eospar%vscale_name)
      else
         if (.not. eospar%linear) then
            ptext(1)='units as volume data'
         else
            ptext(1)='units as length data'
         end if
      end if

      if (.not. eospar%linear) then
         eospar%ParName(1:4) =(/'V0   ','K0   ','Kp   ','Kpp  '/)

         eospar%comment(1) = 'Reference pressure volume: '//trim(ptext(1))
         eospar%comment(2) = 'Bulk modulus: '//trim(ptext(2))
         eospar%comment(3) = 'dK/dP: dimensionless'

         select case(eospar%imodel)
            case (6)
               eospar%ParName(4) = 'Z   '
               eospar%comment(4) = 'N(electrons) in V0'

            case default
               eospar%comment(4) = 'd2K/dP2: '//trim(ptext(4))
         end select

      else
         eospar%ParName(1:4) =(/'L0   ','M0   ','Mp   ','Mpp  '/)

         eospar%comment(1) = 'Reference pressure length: '//trim(ptext(1))
         eospar%comment(2) = 'Linear modulus: '//trim(ptext(2))
         eospar%comment(3) = 'dM/dP: dimensionless'

         select case(eospar%imodel)
            case (6)
               eospar%ParName(4) = 'Z   '
               eospar%comment(4) = 'N(electrons) in V0'

            case default
               eospar%comment(4) = 'd2M/dP2: '//trim(ptext(4))
         end select
      end if

      !> Thermal models are only set in init_EoS_thermal

      !>Adjustl all scale names so that they can be compared in index
      eospar%pscale_name=trim(adjustl(eospar%pscale_name))
      eospar%vscale_name=trim(adjustl(eospar%vscale_name))
   End Subroutine Set_Eos_Names

   !!----
   !!---- SET_EOS_USE
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
   !!---- 17/07/2015
   !!
   Module Subroutine Set_EoS_Use(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

      !---- Local Variables ----!
      integer :: i

      !> Init
      eospar%iuse=0

      !> EoS Model
      select case(eospar%imodel)
         case (0)
            eospar%iuse(1)=1                                      ! None eg thermal only

         case (1)
            eospar%iuse(1:3)=1                                     ! Murnaghan

         case (6)
            eospar%iuse(1:3)=1
            eospar%iuse(4)=2

         case default
            eospar%iuse(1:eospar%iorder)=1                          ! other isothermal EoS
            if (eospar%iorder < 4) eospar%iuse(eospar%iorder+1:4)=3  !implied values
      end select

      !> Thermal Model
      select case(eospar%itherm)
         case (1)             ! Berman
            eospar%iuse(10:11)=1 ! alpha terms
            eospar%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            eospar%iuse(19)=2    ! Grunesien q power law parameter
            eospar%TRef_fixed   = .false.
            eospar%pthermaleos  =.false.

         case (2)             ! Fei
            eospar%iuse(10:12)=1 ! alpha terms
            eospar%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            eospar%iuse(19)=2    ! Grunesien q power law parameter
            eospar%TRef_fixed   = .false.
            eospar%pthermaleos  =.false.

         case (3)             ! HP 1998
            eospar%iuse(10:11)=1 ! alpha terms
            eospar%iuse(18)=2
            eospar%iuse(19)=2    ! Grunesien q power law parameter
            eospar%TRef_fixed   = .false.
            eospar%pthermaleos  =.false.

         case (4)             ! Holland-Powell thermal expansion, in Kroll form
            if (eospar%imodel ==0) eospar%iuse(3)=2     ! require Kprime_zero but not stable in refinement if no P data (added 27/01/2014 RJA)
            eospar%iuse(10)=1    ! alpha at Tref
            eospar%iuse(11)=2    ! Einstein T should be reported but cannot be refined
            eospar%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            eospar%iuse(19)=2    ! Grunesien q power law parameter
            eospar%TRef_fixed   = .false.
            eospar%pthermaleos  =.false.

         case (5)             ! Salje
            eospar%iuse(10:11)=1
            eospar%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            eospar%iuse(19)=2    ! Grunesien q power law parameter
            eospar%TRef_fixed   = .true.
            eospar%pthermaleos  =.false.

         case (6)             ! Thermal pressure in H&P form (no dK/dT): requires a eos model as well
            if (eospar%imodel==0) eospar%iuse(2:3)=2   ! K, Kp refinement is not stable without a pressure model
            eospar%iuse(5:6)=0     ! No dK/dT parameter:
            eospar%iuse(10)=1    ! alpha at Tref
            eospar%iuse(11)=2    ! Einstein T should be reported but cannot be refined
            eospar%iuse(18)=2    ! Grunesien parameter at Pref,Tref
            eospar%iuse(19)=2    ! Grunesien q power law parameter
            eospar%TRef_fixed   = .false.
            eospar%pthermaleos  =.true.

         case (7)             ! Thermal pressure in MGD form
            eospar%iuse(5:6)=0     ! No dK/dT parameter:
            !eospar%iuse(10)=1    ! MGD Gamma0
            eospar%iuse(11)=1    ! Debye T
            !eospar%iuse(12)=1    ! MGD q at Pref,Tref
            eospar%iuse(13)=2    ! Natoms per formula unit
            eospar%iuse(18)=1    ! Grunesien parameter at Pref,Tref for Ks to Kt
            eospar%iuse(19)=1    ! Grunesien q power law parameter for Ks to Kt
            eospar%TRef_fixed   = .false.
            eospar%pthermaleos  = .true.
      end select

      !> Phase transition model
      select case(eospar%itran)
         case (1,2)     ! Landau PV or TV
            eospar%iuse(20)=2           !settable, no refine: sense of transition,
            eospar%iuse(21)=1           !settable, allow refine: Ptr or Ttr
            eospar%iuse(24)=1           !settable, allow refine: aL
            eospar%iuse(25)=1           !settable, allow refine: betaL
            eospar%iuse(26)=1           !settable, allow refine: aH
            eospar%iuse(27)=1           !settable, allow refine: betaH

         case (3)     ! Landau PVT
            eospar%iuse(20:22)=2        !settable, no refine: sense of transition, T(tr), dT(Tr)/dP
            eospar%iuse(24)=1           !settable, allow refine: aL,
            eospar%iuse(23)=2           !settable, fixed d2Tr/dP2
            eospar%iuse(25:27)=1        !settable, allow refine:  betaL,aH,betaH
      end select

      !> Shear model
      select case(eospar%ishear)
         case (0)
            eospar%iuse(30)=0            ! No model, G0 set very large

         case (1)
            eospar%iuse(30:34)=2         ! Polynomial model: settable, not refineable
      end select

      !> Cross terms model
      if (eospar%pthermaleos) then
         eospar%icross=0                ! Kill the cross-terms
         eospar%iuse(5:6)=0
         eospar%params(5:6)=0.

      else
         select case(eospar%icross)
            case (0)
               eospar%iuse(5)=0

            case(1)
               eospar%iuse(5)=1

            case(2)
               eospar%iuse(5)=1
               eospar%iuse(6)=2              !settable, no refine: del-prime refinement always unstable
         end select
      end if

      !> Set the refine flags to be consistent with the use flags
      do i=1,n_eospar
         if (eospar%iuse(i) /=1) eospar%iref(i)=0
      end do
   End Subroutine Set_Eos_Use

   !!----
   !!---- SET_KP_KPP_COND
   !!----    Fix Kp and Kpp values from Model and Order of EoSpar
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine Set_Kp_Kpp_Cond(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

      !---- Local Variables ----!
      real(kind=cp),dimension(N_EOSPar):: ev        ! local copies of room pressure parameter values

      !> Local copy
      call EoS_to_Vec(eospar,ev) !  ev contains volume-like parameters

      select case (eospar%imodel)
         case (1) ! Murnaghan
            return      ! no defaults

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

         case (5) ! Tait with definitions of order derived from Holland and Powell (2011)
            if (eospar%iorder == 2) ev(3)=4.0_cp
            if (eospar%iorder == 2 .or. eospar%iorder == 3)then
               if (abs(ev(2)) > 0.0)ev(4)=-1.0_cp*ev(3)/ev(2)
            endif

         case (6) !No defaults
            return      ! no defaults

      end select

      !> Handle linear or volume
      if (.not. eospar%linear) then
         if (eospar%iorder == 2) eospar%params(3)=ev(3)
         if (eospar%iorder == 2 .or. eospar%iorder == 3) eospar%params(4)=ev(4)

      else
         if (eospar%iorder == 2) eospar%params(3)=ev(3)*3.0_cp
         if (eospar%iorder == 2 .or. eospar%iorder == 3) eospar%params(4)=ev(4)*3.0_cp
      end if
   End Subroutine Set_Kp_Kpp_Cond

   !!--++
   !!--++ SET_SHEAR_NAMES
   !!--++
   !!--++ PRIVATE
   !!--++ Set the character variables in eos_type data structures for SHEAR EoS
   !!--++
   !!--++ 11/07/2016
   !!
   Module Subroutine Set_Shear_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

      !---- Local Variables ----!
      integer :: n

      !> Check for valid model number. If not valid, set zero
      if(eospar%ishear < 0 .or. eospar%ishear > N_SHEAR_MODELS) eospar%ishear=0

      !> Set the Eos name
      eospar%smodel=shearmodel_names(eospar%ishear)

      !> Set upper limit to thermal parameter numbers
      n=34
      if (n > n_eospar) n=n_eospar

      select case(eospar%ishear)
         case (0)
            eospar%parname(30:n) = ' '
            eospar%comment(30:n) = ' '

         case (1)
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
   End Subroutine Set_Shear_Names

   !!--++
   !!--++ SET_THERMAL_NAMES
   !!--++
   !!--++ PRIVATE
   !!--++ Set the character variables in eos_type data structures for thermal EoS
   !!--++
   !!--++ 17/07/2015
   !!
   Module Subroutine Set_Thermal_Names(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar  ! EoS object

      !---- Local Variables ----!
      integer :: n

      !> Check for valid model number. If not valid, set zero
      if (eospar%itherm < -1 .or. eospar%itherm > N_THERM_MODELS) eospar%itherm=0

      !> Set the Eos name
      eospar%tmodel=Tmodel_names(eospar%itherm)

      !> Set upper limit to thermal parameter numbers
      n=19
      if (n > n_eospar) n=n_eospar

      !> Set the V0 name and comment here, in case eos is thermal only
      if (.not. eospar%linear) then
         eospar%ParName(1) ='V0   '
         eospar%comment(1) = 'Reference pressure volume:'
         if (len_trim(eospar%vscale_name) > 0) then
            eospar%comment(1) = trim(eospar%comment(1))//' units are '//trim(eospar%vscale_name)
         else
            eospar%comment(1) = trim(eospar%comment(1))//' units as volume data'
         end if
      else          !Linear
         eospar%ParName(1) ='L0   '
         eospar%comment(1) ='Reference pressure length:'
         if (len_trim(eospar%vscale_name) > 0) then
            eospar%comment(1) = trim(eospar%comment(1))//' units are '//trim(eospar%vscale_name)
         else
            eospar%comment(1) = trim(eospar%comment(1))//' units as length data'
         end if
      end if

      select case(eospar%itherm)
         case (0)
            eospar%parname(10:n) = ' '
            eospar%comment(10:n) = ' '

         case (1)
            eospar%parname(10:11) = (/'alph0','alph1'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Linear term thermal expansion x10^8 K^-2'

         case (2)
            eospar%parname(10:12) = (/'alph0','alph1','alph2'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Linear term thermal expansion x10^8 K^-2'
            eospar%comment(12) = '1/T^2 term thermal expansion, K'

         case (3)
            eospar%parname(10:11) = (/'alph0','alph1'/)
            eospar%comment(10) = 'Constant of thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Sqrt term of thermal expansion x10^4 K^-1/2'

         case (4)    ! Kroll needs Kp as well (in case no pressure eos)
            eospar%parname(3) = 'Kp   '
            eospar%comment(3) = 'dK/dP: dimensionless'
            if (eospar%linear) then
                eospar%parname(3) = 'Mp   '
                eospar%comment(3) = 'dM/dP: dimensionless'
            end if
            eospar%parname(10:11) = (/'alph0','Th_E '/)
            eospar%comment(10) = 'Constant of thermal expansion at Tref x10^5 K^-1'
            eospar%comment(11) = 'Einstein temperature in K'

         case (5)
            eospar%parname(10:11) = (/'p1   ','T_sat'/)
            eospar%comment(10) = 'Approx 3x highT thermal expansion x10^5 K^-1'
            eospar%comment(11) = 'Saturation temperature in K'

         case (6)
            eospar%parname(10:11) = (/'alph0','Th_E '/)
            eospar%comment(10) = 'Constant of thermal expansion at Tref x10^5 K^-1'
            eospar%comment(11) = 'Einstein temperature in K'

         case (7)
            eospar%parname(10:12) = (/'Gamm0','ThMGD','qMGD '/)
       !     eospar%comment(10) = 'MGD Gamma0, dimensionless'  ! no longer used - uses 18 and 19
            eospar%comment(11) = 'Debye temperature in K'
       !     eospar%comment(12) = 'Volume scaling of MGD gamma, dimensionless'
            eospar%parname(13) = 'Natom'
            eospar%comment(13) = 'Number of atoms per formula unit'
      end select

      !> Common terms for all thermal
      eospar%parname(18) = 'Gamm0'
      eospar%comment(18) = 'Gruneisen parameter at Tref,Pref'
      eospar%parname(19) = 'q    '
      eospar%comment(19) = 'Gruneisen power law in V/V0'
   End Subroutine Set_Thermal_Names

   !!--++
   !!--++ SET_TRANSITION_NAMES
   !!--++
   !!--++ PRIVATE
   !!--++ Set the character variables in eos_type data structures for Transition EoS
   !!--++
   !!--++ 17/07/2015
   !!
   Module Subroutine Set_Transition_Names(Eospar)
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
      if (n > n_eospar) n=n_eospar

      select case(eospar%itran)
         case (0)
            eospar%parname(20:n) = ' '
            eospar%comment(20:n) = ' '

         case (1)       ! Landau power law P-V
            eospar%parname(20:27) = (/'High ','Ptr  ','     ','     ','aL   ','betaL','aH   ','betaH'/)
            eospar%comment(20) = 'Indicator = +1 if high P phase is high sym phase'
            eospar%comment(21) = 'Transition pressure'
            eospar%comment(22) = 'dTr/dP: slope of transition boundary'
            eospar%comment(23) = ''
            eospar%comment(24) = 'Scaling parameter, low phase x10^3'
            eospar%comment(25) = 'Power law term, low phase'
            eospar%comment(26) = 'Scaling parameter, high phase x10^3'
            eospar%comment(27) = 'Power law term, high phase'

         case (2)       ! Landau power law V-T
            eospar%parname(20:27) = (/'High ','Ttr  ','     ','     ','aL   ','betaL','aH   ','betaH'/)
            eospar%comment(20) = 'Indicator = +1 if high T phase is high sym phase'
            eospar%comment(21) = 'Transition temperature'
            eospar%comment(22) = ''
            eospar%comment(23) = ''
            eospar%comment(24) = 'Scaling parameter, low phase x10^3'
            eospar%comment(25) = 'Power law term, low phase'
            eospar%comment(26) = 'Scaling parameter, high phase x10^3'
            eospar%comment(27) = 'Power law term, high phase'

         case (3)       ! Landau power law PVT
            eospar%parname(20:27) = (/'High ','Ttr  ','Ttr-P','TtrP2','aL   ','betaL','aH   ','betaH'/)
            eospar%comment(20) = 'Indicator = +1 if high T phase is high sym phase'
            eospar%comment(21) = 'Transition temperature at P=0'
            eospar%comment(22) = 'Ttr=Ttr0 + uP + vP^2: P coeff'
            eospar%comment(23) = 'Ttr=Ttr0 + uP + vP^2: P^2 coeff'
            eospar%comment(24) = 'Scaling parameter, low phase x10^3'
            eospar%comment(25) = 'Power law term, low phase'
            eospar%comment(26) = 'Scaling parameter, high phase x10^3'
            eospar%comment(27) = 'Power law term, high phase'
      end select
   End Subroutine Set_Transition_Names
   
End SubModule EoS_022   