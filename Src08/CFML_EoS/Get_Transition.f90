!!----
!!----
!!----
SubModule (CFML_EoS) EoS_003
   Contains
   !!----
   !!---- GET_TRANSITION_PRESSURE
   !!----    Returns the transition pressure at this T
   !!----
   !!---- 17/07/2015
   !!
   Module Function Get_Transition_Pressure(T,EosPar) Result(Ptr)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: Ptr     ! The transition T at this P

      !---- Local Variables ----!
      real(kind=cp)              :: sqroot

      !> Init
      ptr=0.0_cp

      !> Check for valid model number. If not valid, return with zero Tr
      if (eospar%itran < 1 .or. eospar%itran > N_TRANS_MODELS) return

      !> Now determine if in high-symmetry or low field
      select case(eospar%itran)
         case(1) ! Landau PV
            ptr=eospar%params(21)

         case(3) ! Landau PVT
            !ptr = (T-eospar%params(21))/eospar%params(22) original linear
            if (abs(eospar%params(23)) < tiny(0.0) ) then
               ptr = (T-eospar%params(21))/eospar%params(22)        ! linear
            else
               sqroot=eospar%params(22)*eospar%params(22)-4.0_cp*eospar%params(23)*(eospar%params(21)-T)
               if (sqroot > tiny(0.0) ) then
                  ptr=2.0_cp*(T-eospar%params(21))/(eospar%params(22)+sqrt(sqroot))        ! Viet/Muller formula for root
               else if(sqroot > -1.0*tiny(0.0) ) then
                  ptr=-2.0_cp*(T-eospar%params(21))/eospar%params(22)
               else
                  err_CFML%Ierr=1 
                  write(err_CFML%Msg,'(a,f8.3,a)')'No real solution for Ptr at T = ',T,'K'
               end if
            end if
      end select
   End Function Get_Transition_Pressure

   !!----
   !!---- GET_TRANSITION_STRAIN
   !!----    Returns the strain at P,T due to the transition, including any softening
   !!----    in the high-symm phase for Volume eos returns the volume strain term, for
   !!----    linear eos it returns the linear term!!
   !!----    Vs is defined relative to the 'bare' eos of the high phase (ie the high
   !!----    phase without transition effects)Returns the transition pressure at this T
   !!----
   !!---- 16/02/2015
   !!
   Module Function Get_Transition_Strain(P,T,EosPar) Result(vs)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: Vs      ! The volume strain

      !----Local Variables ----!
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
   End Function Get_Transition_Strain

   !!----
   !!---- GET_TRANSITION_TEMPERATURE
   !!----    Returns the transition temperature at this pressure
   !!----
   !!---- 16/02/2015
   !!
   Module Function Get_Transition_Temperature(P,EosPar) Result(Tr)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: Tr       ! The transition T at this P
      
      !---- Local Variables ----!

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
   End Function Get_Transition_Temperature
   
   !!----
   !!---- TRANSITION_PHASE
   !!----    Returns .true. if P and T are in the low phase stability field
   !!----    and .false. (default) if in the high-symm field, or exactly on
   !!----    the boundary.
   !!----
   !!---- 17/07/2015
   !!
   Module Function Transition_Phase(P,T,Eospar) Result(Ip)
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
         case (1) ! Landau PV
            if (P < eospar%params(21)) ip=.true.

         case (2) ! Landau TV
            if (T < eospar%params(21)) ip=.true.

         case (3) ! Landau PVT
            !Ttr = eospar%params(21)+p*eospar%params(22)  changed from this on 18/12/2015
            Ttr = Get_Transition_Temperature(P,EosPar)  ! general case
            if ( T < Ttr ) ip=.true.
      end select

      !> Now invert result if the lowT phase is high symm phase:
      if (eospar%params(20) < 0) ip = .not. ip
   End Function Transition_Phase
   
End SubModule EoS_003   