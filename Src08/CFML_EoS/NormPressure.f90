!!----
!!----
!!----
SubModule (CFML_EoS) EoS_005
   Contains
   !!--++
   !!--++ NORMPRESSURE_EOS
   !!--++
   !!--++ PRIVATE
   !!--++ Returns the value of Normalised Pressure (F) at this Strain (S) using
   !!--++ the EoS parameters as K0, V0, Kp, Kpp
   !!--++
   !!--++ Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!--++ Modified for thermal EoS: requires T to be meaningful
   !!--++
   !---++ 10/09/2013
   !!
   Module Function NormPressure_Eos(S,T,EosPar) Result(F)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S       ! Strain
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)                     :: F

      !---- Local Variables ----!
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
         !> simple PV eos
         k0=ev(2)
         kp=ev(3)
         kpp=ev(4)

      else if(eospar%itran == 0 .and. eospar%itherm /= 6) then
         !> normal thermal expansion models with dK/dT and no transition
         k0=Get_K0_T(T,eospar)                     ! returns M(T) for linear,
         if (eospar%linear) k0=k0/3.0_cp
         kp=ev(3)
         kpp=ev(4)

      else
         !> Transition model, or Pthermal: cannot use get_V0_T or get_K0_T because they return V and K at Tref for pthermal
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
   End Function NormPressure_Eos

   !!--++
   !!--++ NORMPRESSURE_P
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
   !!--++ 15/02/2013
   !!
   Module Function NormPressure_P(S,P,imodel) Result(F)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S      ! Strain
      real(kind=cp),  intent(in) :: P      ! Presure
      integer      ,  intent(in) :: imodel ! type of eos
      real(kind=cp)              :: F

      !---- Local Variables ----!

      !> Init: a zero value is returned for undefined values
      f=0.0_cp

      !> If the finite strain is zero, the normalised pressure is not defined
      if (s > tiny(0.0) ) then
         select case (imodel)
            case (1,5,6) ! Murnaghan, Tait, APL
               f=0.0_cp

            case (2) ! Birch-Murnaghan
               f=p/3.0_cp/s/(1.0_cp+2.0_cp*s)**2.5_cp

            case (3) ! Vinet: new definition RJA 28/03/2013
               f= (p*(1.0_cp - s)**2.0_cp) / 3.0_cp / s

            case (4) ! Natural Strain
               f=p/3.0_cp/s * exp(-3.0_cp*s)
         end select
      end if
   End Function NormPressure_P

   !!----
   !!---- PRESSURE_F
   !!----    Returns the value of Pressure (P) from the input Normalized Pressure(F)
   !!----    and Strain (S) according to the EoS model in EosPar.
   !!----
   !!--.. NOTE:
   !!--..    The Eos parameters are NOT used in this calculation, on the
   !!--..    eospar%imodel (ie type of Eos)
   !!----
   !!---- 21/02/2013
   !!
   Module Function Pressure_F(F,S,EosPar) Result(P)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: F       ! Normalized Pressure
      real(kind=cp),  intent(in) :: S       ! Strain
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: P

      !---- Local Variables ----!
      real(kind=cp) :: cF

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
   End Function Pressure_F

End SubModule EoS_005
