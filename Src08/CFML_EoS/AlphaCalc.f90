!!----
!!----
!!----
SubModule (CFML_EoS) EoS_013
   Contains
   
   !!----
   !!---- ALPHA_CAL
   !!----   Calculate the alpha parameter in thermal case
   !!----   For linear case alpha is correct as 1/a da/dT
   !!----
   !!---- 17/07/2015
   !!
   Module Function Alpha_Cal(P,T,Eospar, DeltaT) Result(Alpha)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: P       ! Pressure
      real(kind=cp),            intent(in) :: T       ! Temperature
      type(Eos_Type),           intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T value
      real(kind=cp)                        :: Alpha

      !---- Local Variables ----!
      integer                        :: j
      real(kind=cp), dimension(-2:2) :: v  ! array for calc v values
      real(kind=cp)                  :: del, tt, delmin, tlimit, tr, alphaest, vlimit

      !> Init
      alpha=0.0_cp
      Tlimit=0.0_cp

      !> PTV table
      if (eospar%imodel == -1) then
         alpha=get_props_ptvtable(p, t, 0.0, eospar, 'AL')
         return
      end if

      !> Need to trap numerical problems with Kroll, Salje, Pthermal at low T
      select case(eospar%itherm)
         case(0) ! no thermal parameters
            return

         case(4:7) ! Kroll and Pthermal: For T < 0.05T(Einstein) alpha=0. Same for Salje but test is 0.05Tsat
            if (t < 0.05_cp*eospar%params(11) ) return
      end select

      !> Numerical solutions: set step size
      del =abs(0.001_cp/eospar%params(10))      ! Set step in T to get about 0.1% shift in V
      if (del > 80.0) del=80.0_cp               ! otherwise for small alpha, del is too big and alpha is inaccurate
      if (present(deltaT)) del=deltaT

      select case(eospar%itherm)           ! adjustment of step size
         !> Fei, HP98
         case(2,3)                         ! T ok, but do not step down into invalid region
            if (abs(eospar%params(12)) > tiny(0.0_cp) ) then
               delmin=abs(t-tlimit)/2.1_cp
               if (del > delmin) del=delmin
            end if

         !> Kroll, Salje, Thermal pressure
         case(4,5,6)
            delmin=(t-0.025_cp*eospar%params(11))/2.0_cp     ! do not allow step in to area where alpha=0
            if (del > delmin) del=delmin                     ! ensures T at all steps is positive

         !> MGD Pthermal
         case(7)  ! so no alpha available for estimation: changed from T+100 to T-100 to avoid going into area where large V invalid
            alphaest=(get_volume(p,t,eospar)-get_volume(p,t-100._cp,eospar))/get_volume(p,t-50._cp,eospar)/100._cp
            del=abs(0.001_cp/alphaest)
            delmin=(t-0.025_cp*eospar%params(11))/2.0_cp     ! do not allow step in to area where alpha=0
            if (del > delmin) del=delmin                     ! ensures T at all steps is positive
            ! now stop the del taking us into illegal area
            tlimit=t+2.0*del
            do                                        ! search for positive K at this P
               if (get_K(p,tlimit,eospar) > 0.0_cp .and. err_CFML%IErr==0)exit
               call clear_error()
               tlimit=tlimit-0.1_cp*(tlimit-t)
               if (tlimit < t) exit                    ! should never happen because P,T is valid    
            end do
            del=0.4*abs(tlimit-t)
      end select

      !> Stop calculation going across a phase boundary
      if (eospar%itran > 0)then
         Tr=get_transition_temperature(p,eospar)
         if (transition_phase(P,T,eospar) .neqv. transition_phase(P,T+2.0*del,eospar)) del=abs(T-Tr)/2.1_cp
         if (transition_phase(P,T,eospar) .neqv. transition_phase(P,T-2.0*del,eospar)) del=abs(T-Tr)/2.1_cp
         if (del < 1.0_cp) del=1.0_cp
      end if

      !> Do the numerical solution
      do j=-2,2,1
         tt=t+real(j)*del                 ! apply shift to temp
         v(j)=get_volume(p,tt,eospar)     ! calc resulting V
      end do

      alpha=(v(-2)+8.0_cp*(v(1)-v(-1))-v(2))/(12.0_cp*del)/v(0)     ! Derivative to second order approximation

      !> Trap non-physical results that arise due to rounding errors
      !> But negative thermal expansion is allowed, so original test alpha < 0 not valid
      !> Removed 1-Feb-2019 because small T is tested for these itherm eos at start
      !select case(eospar%itherm)
      !    case(4,5,6,7)
      !      if (alpha < tiny(0.0)) alpha=0.0_cp
      !end select

      !> get_volume returns 'a' for linear, so alpha is correct as 1/a da/dT for linear

   End Function Alpha_Cal
   
End SubModule EoS_013   