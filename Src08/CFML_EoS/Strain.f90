!!----
!!----
!!----
SubModule (CFML_EoS) EoS_006
   Contains
   !!----
   !!---- STRAIN
   !!----   Returns the value of Strain (S) at this V according with the EoS model
   !!----
   !!--.. NOTE:
   !!--..    The values of EosPar are NOT used. Only the type of Eos function is needed!
   !!--..    For linear eos, the strain is that of a^3
   !!--..    If V/V0 = 0. then an error is set
   !!----
   !!---- 15/02/2013
   !!
   Module Function Strain(VV0,EosPar) Result(S)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: VV0     ! Volume divided by V0 at this temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: S

      !---- Local Variables ----!
      real(kind=cp) :: cvv0

      !> Init
      s=0.0_cp
      if (vv0 <= 0.00001) then
         err_CFML%IErr=1
         err_CFML%Msg="Strain calculation called with invalid V/V0 =< 0: strain set zero"
         return
      end if

      !> Local copy
      cvv0=vv0
      if (eospar%linear) cvv0=vv0**3.0_cp

      select case (eospar%imodel)
         case (1,5,6) ! Murnaghan, Tait, APL
            s=0.0_cp

         case (2) ! Birch-Murnaghan
            s=(cvv0**(-2.0_cp/3.0_cp) - 1.0_cp)/2.0_cp

         case (3) ! Vinet: new definition RJA 28/03/2013
            s= 1.0_cp - cvv0**(1.0_cp/3.0_cp)

         case (4) ! Natural Strain: the original definition of strain was wrong in v5 by a change in sign
            s= -1.0_cp*log(cvv0)/3.0_cp
      end select
   End Function Strain

   !!----
   !!---- STRAIN_EOS
   !!----   Returns the value of Strain (S) at this V and T according with the EoS parameters
   !!----
   !!---- 05/09/2013
   !!
   Module Function Strain_EOS(V,T,EosPar) Result(S)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameters
      real(kind=cp)              :: S

      !---- Local Variables ----!
      real(kind=cp) :: vvo
      type(Eos_Type):: E  ! Eos Parameters copy

      !> Init
      s=0.0_cp
      if (v <= 0.00001) then
         err_CFML%IErr=1
         err_CFML%Msg="Strain calculation called with invalid V =< 0: strain set zero"
         return
      end if
      e=eospar

      !> Calculate
      vvo=v/get_volume(0.0_cp,t,e)     ! vv0 is V(P,T)/V(0,T) or a(P,T)/a(0,T)
      s=strain(vvo,e)                  ! cubes vv0 on input if linear
   End Function Strain_EOS
   
End SubModule EoS_006   