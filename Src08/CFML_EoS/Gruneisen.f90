!!----
!!----
!!----
SubModule (CFML_EoS) EoS_008
   Contains
   !!----
   !!---- GET_GPT
   !!----    Obtain the value of G (or Glinear) at P and T
   !!----
   !!---- 11/07/2016
   !!
   Module Function Get_GPT(P,T,EoSPar) Result(gpt)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P      ! Pressure
      real(kind=cp),  intent(in) :: T      ! Temperature
      type(Eos_Type), intent(in) :: EosPar ! EoS Variable
      real(kind=cp)              :: GPT

      !---- Local Variables ----!
      integer         :: i
      real(kind=cp)   :: delp

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
   End Function Get_GPT

   !!----
   !!---- GET_GRUN_PT
   !!----    Returns Gruneisen parameter at this P,T as gamma0*(V/V0)
   !!----
   !!---- 18/07/2016
   !!
   Module Function Get_Grun_PT(P,T,Eospar) Result(G)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: T       ! Temperarture
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: G

      !---- Local Variables ----!
      real(kind=cp) :: v

      v=get_volume(P,T,eospar)
      G=Get_Grun_V(v,Eospar)
   End Function Get_Grun_PT

   !!----
   !!---- GET_GRUN_V
   !!----    Returns Gruneisen parameter at this volume as gamma0*(V/V0)
   !!----    If linear it calculates gamma0*(a/a0)^3
   !!----
   !!---- 18/07/2016
   !!
   Module Function Get_Grun_V(V,Eospar)  Result(Grun)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume or length if linear
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: GRun    ! The resulting gruneisen gamma

      !---- Local Variables ----!
      real(kind=cp) :: VV0

      !> Must be careful with transitions because eospar%params(1) is the high phase V0
      !> V0=get_volume(eospar%pref,eospar%tref,eospar) (Nov 2016)
      !> no I don't think so. Grun is a property of the high phase
      VV0=V/eospar%params(1)
      if (eospar%linear) VV0=VV0**3.0_cp

      if (abs(eospar%params(19)) > 0.00001_cp) then
         VV0=VV0**eospar%params(19)
      else
         VV0=1.0_cp
      end if

      Grun=eospar%params(18)*VV0
   End Function Get_Grun_V
   
End SubModule EoS_008   