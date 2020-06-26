!!----
!!----
!!----
SubModule (CFML_EoS) EoS_Property
   implicit none
   Contains
   !!----
   !!---- GET_PROPERTY_X
   !!----    Returns the ...
   !!----
   !!---- 16/03/2017
   !!
   Module Function Get_Property_X(P,T,Eos,Xtype) Result(val)
      !---- Arguments ----!
      real(kind=cp),     intent(in) :: P       ! Pressure
      real(kind=cp),     intent(in) :: T       ! Temperature
      type(Eos_Type),    intent(in) :: EoS     ! Eos Parameter
      integer, optional, intent(in) :: xtype   ! =0 when X=V, =1 for X=K (isothermal)
      real(kind=cp)                 :: val

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
            agt=Alpha_Cal(P,T,Eos)*get_grun_v(vol,eos)*T        ! Get_Grun knows about linear
            if (eos%linear) agt=3.0_cp*agt
            val=(1.0_cp+agt)*K_Cal(Vol,T,Eos,p=p)

         case default
            val=0.0
      end select
   End Function Get_Property_X

End SubModule EoS_Property