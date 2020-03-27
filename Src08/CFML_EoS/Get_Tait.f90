!!----
!!----
!!----
SubModule (CFML_EoS) EoS_127
   Contains

   !!--++
   !!--++ GET_TAIT
   !!--++
   !!--++ PRIVATE
   !!--++ Returns a,b,c Tait parameters in a vector.
   !!--++
   !!--++ 17/07/2015
   !!
   Module Function Get_Tait(Eospar,T) Result(Vec)
      !---- Arguments ----!
      type(Eos_Type),            intent(in)  :: Eospar  ! Eos Parameters
      real(kind=cp),             intent(in)  :: T       ! Temperature
      real(kind=cp),dimension(3)             :: Vec     ! Vector (a,b,c) of Tait parameters

      !---- Local Variables ----!
      real(kind=cp) :: k0,kp,kpp

      !> Init
      Vec=0.0_cp
      if (eospar%imodel /= 5) return

      select case(eospar%itherm)
         case (1:5)                  ! normal thermal expansion models with dK/dT
            k0 =Get_K0_T(T,eospar)
            kp =Get_Kp0_T(T,eospar)

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
   End Function Get_Tait

End SubModule EoS_127
