!!----
!!----
!!----
SubModule (CFML_EoS) EoS_021
   Contains

   !!----
   !!---- INIT_EOS_CROSS
   !!----    Initialize the EoS Type for P-T cross-terms
   !!----
   !!---- 11/11/2016
   !!
   Module Subroutine Init_EoS_Cross(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !> Check for valid model number. If not valid, set zero
      if (eospar%icross < 0 .or. eospar%icross > N_CROSS_MODELS) eospar%icross=0
      if (eospar%pthermaleos) eospar%icross=0

      eospar%params(5)           = 0.0_cp
      eospar%params(6)           = 0.0_cp
      eospar%vcv(5:6,1:N_EOSPAR) = 0.0_cp
      eospar%vcv(1:N_EOSPAR,5:6) = 0.0_cp
      eospar%factor(5)           = 1.0_cp
      eospar%factor(6)           = 1.0_cp

      call Set_Cross_Names(Eospar)    ! Set the variable names
      call Set_Eos_Use(Eospar)        ! update the use flags
   End Subroutine Init_EoS_Cross

   !!----
   !!---- INIT_EOS_DATA_TYPE
   !!----    Initialize EoS_Data_Type
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine Init_EoS_Data_Type(E)
      !---- Arguments ----!
      type (EoS_Data_Type), intent(in out)   :: E

      !> Init
      E%IUse = 0
      E%IGrp = 0
      E%T    = 298.0
      E%P    = 0.0_cp
      E%V    = 0.0_cp
      E%cell = 0.0_cp
      E%ang  = 0.0_cp

      E%SigT = 0.0_cp
      E%SigP = 0.0_cp
      E%SigV = 0.0_cp
      E%sigc = 0.0_cp
      E%siga = 0.0_cp
   End Subroutine Init_EoS_Data_Type

   !!----
   !!---- INIT_EOS_SHEAR
   !!----    Initialize the EoS Type for Shear case
   !!----
   !!---- 17/02/2015
   !!
   Module Subroutine Init_EoS_Shear(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !---- Variables ----!
      integer  :: n

      !> Check for valid model number. If not valid, set zero
      if (eospar%ishear < 0 .or. eospar%ishear > N_SHEAR_MODELS) eospar%ishear=0

      !> Set upper limit to parameter numbers
      n=34
      if (n > N_EOSPAR) n=N_EOSPAR

      select case(eospar%ishear)
         case (0)
            eospar%params(30)          = huge(0.0_cp)   ! default is infinitely stiff
            eospar%vcv(30,1:n)         = 0.0_cp
            eospar%vcv(30:N_EOSPAR,1:n)= 0.0_cp
            eospar%vcv(1:N_EOSPAR,30:n)= 0.0_cp
            eospar%factor(30:n)        = 1.0_cp

         case (1)        ! polynomial
            eospar%params(30)    = 100.0_cp      ! G0 at Pref Tref
            eospar%params(31:34) =   0.0_cp      ! Polynomial coefficients
            eospar%factor(30:n)  =   1.0_cp
      end select

      call Set_Shear_Names(Eospar)    ! Set the variable names
      call Set_Eos_Use(Eospar)        ! update the use flags
   End Subroutine Init_EoS_Shear

   !!----
   !!---- INIT_EOS_THERMAL
   !!----    Initialize the EoS Type for Thermal case
   !!----
   !!---- 10/09/2013
   !!
   Module Subroutine Init_EoS_Thermal(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !---- Variables ----!
      integer    :: n

      !> Check for valid model number. If not valid, set zero
      if(eospar%itherm < -1 .or. eospar%itherm > N_THERM_MODELS) eospar%itherm=0

      !> Set upper limit to thermal parameter numbers
      n=19
      if (n > N_EOSPAR)n=N_EOSPAR

      eospar%alphafactor=1.0E5_cp                      ! Normal scale factor for printing values of alpha

      select case(eospar%itherm)
         case (-1)           ! PTV table
            eospar%factor(10)   = 1.0E5_cp            ! factor to multiply alpha values on printing
            eospar%pthermaleos  =.false.

         case (0)
            eospar%params(5)    = 0.0_cp
            eospar%params(10:n) = 0.0_cp
            eospar%vcv(5,1:n)   = 0.0_cp
            eospar%vcv(10:n,1:n)= 0.0_cp
            eospar%vcv(1:n,10:n)= 0.0_cp
            eospar%factor(10:n) = 1.0_cp
            eospar%TRef         = 298.0_cp
            eospar%pthermaleos  =.false.

         case (1)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp
            eospar%TRef        = 298.0_cp             ! Simple thermal expansion,
            eospar%pthermaleos  =.false.

         case (2)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E8_cp
            eospar%factor(12)  = 1.0_cp
            eospar%TRef        = 298.0_cp             ! Simple thermal expansion,
            eospar%pthermaleos  =.false.

         case (3)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0E4_cp
            eospar%TRef        = 298.0_cp             ! Simple thermal expansion,
            eospar%pthermaleos  =.false.

         case (4)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef        = 298.0_cp             ! Holland and Powell thermal expansion without P
            eospar%TRef_fixed  = .false.
            eospar%params(11)  = 298.0_cp             ! Einstein temperature
            eospar%pthermaleos  =.false.

         case (5)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef        = 0.0_cp               ! Salje thermal expansion
            eospar%TRef_fixed  = .true.
            eospar%params(11)  = 298.0_cp             ! Saturation temperature
            eospar%pthermaleos  =.false.

         case (6)
            eospar%factor(10)  = 1.0E5_cp             ! factor to multiply values on printing
            eospar%factor(11)  = 1.0_cp
            eospar%TRef        = 298.0_cp
            eospar%TRef_fixed  = .false.
            eospar%params(11)  = 298.0_cp             ! Einstein temperature default
            eospar%pthermaleos  =.true.

         case(7)                                      ! MGD Pthermal (only uses params (11) and (13)
            eospar%factor(10:14)  = 1.0_cp
            eospar%TRef           = 298.0_cp
            eospar%TRef_fixed     = .false.
            eospar%params(11)     = 298.0_cp          ! Debye temperature default
            eospar%params(13)     = 1.0               ! Natoms/molecule for MGD
            eospar%pthermaleos  =.true.
      end select

      !> Set the common terms for the Gruneisen model/Debye model
      eospar%params(18)=1.0_cp      ! gamma0
      eospar%params(19)=1.0_cp      ! q

      call Init_EoS_Cross(Eospar)                     ! init the cross-terms
      call Set_Thermal_Names(Eospar)                  ! Set the variable names
      call Set_Eos_Use(Eospar)                        ! update the use flags and other pointers
   End Subroutine Init_EoS_Thermal

   !!----
   !!---- INIT_EOS_TRANSITION
   !!----    Initialize the EoS Type for Transition case
   !!----
   !!---- 17/02/2015
   !!
   Module Subroutine Init_EoS_Transition(Eospar)
      !---- Arguments ----!
      type (EoS_Type), intent(in out) :: Eospar

      !---- Variables ----!
      integer  :: n

      !> Check for valid model number. If not valid, set zero
      if (eospar%itran < 0 .or. eospar%itran > N_TRANS_MODELS) eospar%itran=0

      !> Set upper limit to parameter numbers
      n=29
      if (n > N_EOSPAR) n=N_EOSPAR

      select case(eospar%itran)
         case (0)
            eospar%params(20:n)        = 0.0_cp
            eospar%vcv(20,1:n)         = 0.0_cp
            eospar%vcv(20:n_eospar,1:n)= 0.0_cp
            eospar%vcv(1:n_eospar,20:n)= 0.0_cp
            eospar%factor(20:n)        = 1.0_cp

         case (1)        ! Landau PV
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

         case (2)        ! Landau TV
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

         case (3)        ! Landau PVT
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
   End Subroutine Init_EoS_Transition

   !!----
   !!---- INIT_EOS_TYPE
   !!----   Initialize EoS_Type setting all parameters to sensible values
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine Init_EoS_Type(Eospar,CLin,IThermal,ITransition,Ishear,Icross)
      !---- Arguments ----!
      type (EoS_Type),            intent(out)    :: Eospar       ! EoS Type
      character(len=*), optional, intent(in)     :: CLin         ! Character variable to indicate linear EoS or not
      integer,          optional, intent(in)     :: IThermal     ! integer to indicate ithermal type
      integer,          optional, intent(in)     :: ITransition  ! integer to indicate transition type
      integer,          optional, intent(in)     :: IShear       ! integer to indicate shear type
      integer,          optional, intent(in)     :: ICross       ! integer to indicate cross-terms type

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
      eospar%savedate=' '
      eospar%pscale_name=' '
      eospar%vscale_name=' '
      call Set_Eos_Names(Eospar)         ! also sets the print factors for the pressure part

      Eospar%PRef     = 0.0_cp
      Eospar%Density0 = 0.0_cp

      Eospar%Iuse     =0
      Eospar%Iuse(1:4)=1                 ! Vo, Ko, Kp, Kpp

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
      if (present(ithermal) .and. ithermal > -2 )then
         Eospar%Itherm = ithermal
         call Init_Eos_Thermal(Eospar)              ! set up default values, names for specific thermal eqn
      end if

      !> Test for optional argument for transition
      Eospar%ITran  =0
      if (present(itransition) .and. itransition  > -1 )then
         Eospar%ITran = itransition
         call Init_Eos_Transition(Eospar)              ! set up default values, names for specific transition
      end if

      !> Test for optional argument for cross terms
      Eospar%ICross  =0
      if (present(icross) .and. icross  > -1 )then
         Eospar%Icross = icross
         call Init_Eos_Cross(Eospar)              ! set up default values, names for specific crosstrem model
      end if

      !> Test for optional argument for shear model
      Eospar%Ishear  =0
      if (present(ishear) .and. ishear  > -1 )then
         Eospar%Ishear = ishear
      end if

      call Init_Eos_Shear(Eospar)              ! set up default values, names for specific shear model
   End Subroutine Init_EoS_Type

End SubModule EoS_021