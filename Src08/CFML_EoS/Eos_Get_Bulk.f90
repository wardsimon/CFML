!!----
!!----
!!----
SubModule (CFML_EoS) EoS_Bulk
   implicit none
   Contains

   !!----
   !!---- GET_K
   !!----    Returns the value of K (or M if linear) at P and T
   !!----    Works by using get_volume(P,T) and then using the V in k_cal
   !!----
   !!---- 18/07/2016
   !!
   Module Function Get_K(P,T,EosPar) Result(K)
      !---- Arguments ----!
      real(kind=cp), intent(in)        :: p       ! Pressure
      real(kind=cp), intent(in)        :: t       ! Tenperature
      type(Eos_Type),intent(in)        :: Eospar  ! Eos Parameters
      real(kind=cp)                    :: k

      !---- Local Variables ----!
      real(kind=cp)                   :: v

      v=get_volume(p,t,eospar)
      k=k_cal(v,t,eospar,p=p)
   End Function Get_K

   !!----
   !!---- GET_Kp
   !!----    Returns the value of K (or M if linear) at P and T
   !!----    Works by using get_volume(P,T) and then using the V in k_cal
   !!----
   !!---- 1/03/2018
   !!
   Module Function Get_Kp(P,T,EosPar) Result(Kp)
      !---- Arguments ----!
      real(kind=cp),intent(in)        :: p       ! Pressure
      real(kind=cp),intent(in)        :: t       ! Tenperature
      type(Eos_Type),intent(in)       :: Eospar  ! Eos Parameters
      real(kind=cp)                   :: kp

      !---- Local Variables ----!
      real(kind=cp)                   :: v

      v=get_volume(p,t,eospar)
      kp=kp_cal(v,t,eospar,p=p)
   End Function Get_Kp

   !!--++
   !!--++ GET_K0_T
   !!--++
   !!--++ PRIVATE
   !!--++ Returns k0 needed for Eos calculations at T this means for Pthermal,
   !!--++ k at Tref is returned.
   !!--++ In the linear case then  M(T, P=0) from params is returned
   !!--++
   !!--++ If k is calculated as being negative, an error message is set
   !!--++ and k0 is returned as the value at Tref for safety.
   !!--++
   !!--++ 17/07/2015
   !!
   Module Function Get_K0_T(T,Eospar) Result(k0)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type),intent(in)  :: Eospar  ! Eos Parameters
      real(kind=cp)              :: k0

      !---- Local Variables ----!
      real(kind=cp) :: vr

      !> Init
      k0=eospar%params(2) !default (and correct if no thermal model)

      select case(eospar%itherm)
         case(1:5)  !> Normal isothermal EoS at T
            select case(eospar%icross)
               case(1)
                  k0=eospar%params(2)+eospar%params(5)*(t-eospar%tref)  !Old linear variation of K with T

               case(2)
                  vr=eospar%params(1)/Get_V0_T(T,EosPar)          ! Get_V0_T returns a0 for linear
                  if (eospar%linear) vr=vr**3.0_cp
                  if (vr > 0.001 .and. vr < huge(0.0_cp)) k0=eospar%params(2)*vr**eospar%params(5)   !Anderson Gruneisen approach using params(5) as delta
            end select

         case(6,7) !> Thermal pressure model
            !k0=eospar%params(2)   ! Is set in the beginning of procedure
      end select
   End Function Get_K0_T

   !!--++
   !!--++ GET_KP0_T
   !!--++ Returns kp0 needed for Eos calculations at T this means for Pthermal,
   !!--++ kp at Tref is returned.
   !!--++ In the linear case then  Mp(T, P=0) from params is returned
   !!--++
   !!--++ 17/11/2016
   !!
   Module Function Get_Kp0_T(T,Eospar) Result(kp0)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: kp0

      !---- Local Variables ----!
      real(kind=cp) :: vr

      !> Init
      kp0=eospar%params(3) !default (and correct if no thermal model)

      select case(eospar%itherm)
         case(1:5)  ! Normal isothermal eos at T
            select case(eospar%icross)
               case(1)  !Old linear variation of K with T, no Kp variation

               case(2)
                  vr=Get_V0_T(T,EosPar)/eospar%params(1)            !normally Vr > 0, but if negative thermal expansion, not
                  if (eospar%linear) vr=vr**3.0_cp
                  if (vr > 0.001 .and. vr < huge(0._cp) ) kp0=eospar%params(3)*vr**eospar%params(6)   ! params(6) is delta-prime
            end select

         case(6,7)        ! Thermal pressure model: value at Tref

      end select
   End Function Get_Kp0_T

   !!--++
   !!--++ GET_KPP0_T
   !!--++  Returns kpp0 needed for Eos calculations at T this means for pthermal,
   !!--++  kpp at Tref is returned.
   !!--++  In the linear case then  Mp(T, P=0) from params is returned
   !!--++
   !!--++ 17/11/2016
   !!
   Module Function Get_Kpp0_T(T,Eospar) Result(kpp0)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: kpp0

      !---- Local Variables ----!
      type(Eos_Type) :: eost     ! workspace

      !> Init
      select case(eospar%imodel)
         case(6)
            kpp0=0.0_cp

         case default
            kpp0=eospar%params(4) !default (and correct if no thermal model, or if Pthermal model)
      end select

      select case(eospar%itherm)
         case(1:5)  ! Normal isothermal eos at T
            eost=eospar           ! Modify eost to hold K0 and Kp at the T of interest

            eost%params(2)=Get_K0_T(T,Eospar)
            eost%params(3)=Get_Kp0_T(T,Eospar)
            !call Set_Kp_Kpp_Cond(Eost)
            call Set_Eos_Implied_Values(Eost)
            kpp0=eost%params(4)
      end select
   End Function Get_Kpp0_T

End SubModule EoS_Bulk