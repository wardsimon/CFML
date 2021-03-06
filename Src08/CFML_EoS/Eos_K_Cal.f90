!!----
!!----
!!----
SubModule (CFML_EoS) EoS_K_Cal
   implicit none
   Contains

   !!----
   !!---- K_CAL
   !!----    Returns value of K at this volume for EoS
   !!----
   !!--..    Changed code to use STRAIN function with v/v0: 27/02/2013
   !!--..    Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- 10/02/2017
   !!
   Module Function K_Cal(V,T,Eospar,P) Result(kc)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: V       ! Volume
      real(kind=cp),            intent(in) :: T       ! Temperature
      type(Eos_Type),           intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: P       ! Pressure if known
      real(kind=cp)                        :: kc

      !---- Local Variables ----!
      integer                            :: j
      real(kind=cp)                      :: vv0,k0,kp,kpp,f,vol
      real(kind=cp)                      :: a,b,c,nu,vt,dPdV,delv
      real(kind=cp)                      :: x,pFG0,c0,c2,da,db,dc              ! For APL
      real(kind=cp)                      :: vs,Ttr,dVs                         ! for transition calculations
      real(kind=cp), dimension(N_EOSPAR) :: ev
      real(kind=cp), dimension(3)        :: abc      ! Tait parameters
      real(kind=cp), dimension(-2:2)     :: pcal     ! array for calc p values numeric solutions
      real(Kind=cp)                      :: Ptrue    ! Input pressure if present. The true pressure without Pthermal subtracted
      real(Kind=cp)                      :: Pcorr    ! Pressure minus Pthermal.
      type(Eos_Type)                     :: EoST     ! Copy of eos


      !> Init
      kc=eospar%params(2)               ! safe default, K(P,T)=K0

      !> Pressure is needed for Tait, Murngahan, and for spontaneous strain at transitions
      !> This is the 'true' pressure, uncorrected for Pthermal
      if (present(P) )then
         Ptrue=P
      else
         Ptrue=get_pressure(v,t,eospar)
      end if

      !> PTV table
      if (eospar%imodel == -1) then
         if (.not. present(P) ) Ptrue=get_props_ptvtable(0.0,T,V,eospar,'P')
         kc=get_props_ptvtable(Ptrue,T,V,eospar,'K')
         return
      end if

      !> Isothermal EoS: get K0 and Kp0 at this T, then do algebra for normal EoS
      !>Or pthermal at Tref and do additional d(Pth)/dV contribution afterwards
      !>If transition, do calculation for bare phase, and add correction for transition afterwards

      !> Correct the volume to the high phase only if there is a transition
      if (eospar%itran > 0) then
         Vs=get_transition_strain(Ptrue,T,eospar)     ! returns the linear or volume strain as appropriate
         Vol=v/(1.0+vs)
      else
         Vol=v
      end if

      select case (eospar%itherm)
         case (0,6,7)                          ! No thermal model, or we have pthermal, so need params at Tref
            vv0=vol/eospar%params(1)           ! vv0 or aa0
            k0=eospar%params(2)
            kp=eospar%params(3)
            kpp=eospar%params(4)

         case (1:5)
            vv0=vol/get_V0_T(t,eospar)          ! vv0 is  v(p,t)/v(p=0,t): in transition, highphase
            k0=Get_K0_T(T,eospar)               ! returns M(T) for linear,
            if (err_CFML%IErr==1) return        ! exit with value eosparms(2) if k0 at P=0 calculated as negative

            kp=Get_Kp0_T(T,eospar)
            kpp=Get_Kpp0_T(T,eospar)
      end select

      !> Strain for BM, NS, Vinet EoS equations
      f=strain(vv0,eospar)
      if (err_CFML%IErr==1) then
         err_CFML%Msg=trim(err_CFML%Msg)//' called from Kcal'
         return
      end if

      !> adjust for linear
      if (eospar%linear) then
         vol=vol**3.0_cp
         vv0=vv0**3.0_cp
         k0=k0/3.0_cp
         kp=kp/3.0_cp
         kpp=kpp/3.0_cp
      end if

      !> Now do the calculation for an isothermal eos at this T
      select case (eospar%imodel)
         case (1)   !Murnaghan
            Pcorr=Ptrue
            if (eospar%pthermaleos) Pcorr=Ptrue - pthermal(Vol,T,eospar)
            kc=k0+kp*Pcorr

         case (2)   !Birch-Murnaghan
            a=0.0_cp
            b=0.0_cp
            if (eospar%iorder > 2) a=1.5_cp*(kp-4.0_cp)
            if (eospar%iorder ==4) b =(9.0_cp*k0*kpp + (9.0_cp*kp*kp) - 63.0_cp*kp + 143.0_cp)/6.0_cp

            kc=k0*(1.0_cp+2.0_cp*f)**2.5_cp * (1.0_cp+ (7.0_cp+2.0_cp*a)*f + &
               (9.0_cp*a + 3.0_cp*b)*f*f + 11.0_cp*b*f*f*f)

         case (3) ! Vinet: from Schlosser & Ferrante PRB 37:4351; this form avoids pressure (RJA 26/2/2013)
                  ! new definition of f used 28/03/2013 RJA
            nu=1.5_cp*(kp-1)
            kc=k0*exp(nu*f) * (1.0_cp + (1.0_cp + nu)*f - nu*f*f) / (1.0_cp-f)**2.0_cp

         case (4) ! Natural strain: expanded from Poirier and Tarantola PEPI 109:1-8 (RJA 26/2/2013 - 28/03/2013)
            a=0.0_cp
            b=0.0_cp
            c=0.0_cp
            select case (eospar%iorder)
               case (2)
                  a=1.0_cp
               case (3)
                  a=kp-1.0_cp
                  b=1.5_cp*(kp-2.0_cp)
               case (4)
                  a=kp-1.0_cp
                  c=1.5_cp*(1.0_cp + k0*kpp + (kp-2.0_cp) + (kp-2.0_cp)**2.0_cp)
                  b=c+1.5_cp*(kp-2.0_cp)
            end select
            kc= 3.0_cp*k0/vv0*( 1.0_cp/3.0_cp + a*f + b*f*f + c*f*f*f)

         case(5) ! Tait
            abc=get_tait(eospar,t)              ! get_tait returns volume-like parameters even for linear
            Pcorr=Ptrue
            if (eospar%pthermaleos) Pcorr=Ptrue - pthermal(Vol,T,eospar)
            if (abc(2)*Pcorr < -0.999999_cp) then
               err_CFML%IErr=1
               err_CFML%Msg='Pressure yields zero bulk modulus for Tait Eos'  !Kc will be returned as eosparams(2)
               return
            else
               kc=k0*(1.0_cp-abc(1)*(1.0_cp-(1.0_cp + abc(2)*Pcorr)**(-1.0_cp*abc(3))))*(1.0_cp + abc(2)*Pcorr)**(1.0_cp+abc(3))
            end if

         case(6) ! APL
            x=vv0**0.333333_cp
            pFG0=AFERMIGAS*(eospar%params(4)*vv0/vol)**1.66666667_cp
            c0=-1.0_cp*log(3.0_cp*K0/pFG0)           ! assumes V in A^3
            c2=1.5_cp*(Kp-3.0_cp)-c0

            !> terms in pressure expression
            a=1.0_cp/x**5.0_cp*(1.0_cp-x)
            b=exp(c0*(1.0_cp-x))
            c=1.0_cp+c2*x*(1.0_cp-x)

            !> derivatives
            da= -5.0_cp/x**6.0_cp +4.0_cp/x**5.0_cp
            db= -1.0_cp*c0*b
            dc= c2*(1.0_cp-2.0_cp*x)
            kc= -1.0_cp*k0*x*(da*b*c+a*db*c+a*b*dc)
      end select

      !> To this point Kc is the bulk modulus of the 'bare' eos of the high phase if it was an isothermal eos
      !> Or for thermal-pressure EoS it is the bulk modulus at Tref and V

      !> Now correct thermal-pressure EoS for d(Pth)/dV contribution to bulk modulus
      if (eospar%Pthermaleos) then
         select case(eospar%itherm)
            case(7)           !MGD EoS: Do this numerically
               eost=eospar
               eost%itran=0    ! clear any transition terms
               delv=0.01_cp*vol
               do j=-2,2,1
                  vt=vol+real(j)*delv             ! apply shift to v
                  Pcal(j)=pthermal(vt,t,eost)     ! calc resulting thermal pressure
               end do
               dPdV=(Pcal(-2)+8.0_cp*(Pcal(1)-Pcal(-1))-Pcal(2))/(12.0_cp*delv)   !dPth/dV
               kc= kc -1.0_cp*vol*dPdV

            case default      !includes HP thermal-pressure which has dK = 0 along isochor

         end select
      end if

      !>Now correct from volume to linear if required, because transition code works in linear
      if (eospar%linear) kc=kc*3.0_cp

      !> Now handle phase transition
      if (eospar%itran > 0 ) then
         !> Local EoS copy
         call EoS_to_Vec(eospar,ev)  ! but no scaling done

         select case(eospar%itran)
            case (1) ! Landau P-V
               if (abs(Ptrue-ev(21)) < 0.0001) then
                  !> at the transition pressure
                  if (abs(ev(25)-1.0_cp) < 0.0001) then
                     dVs=ev(24)*ev(25)           ! exact second order
                  else
                     dvs=huge(0._cp)             ! effectively infinite
                  end if
               else if (transition_phase(Ptrue,T,eospar) ) then             ! in the low phase
                  dVs=ev(24)*ev(25)*abs(Ptrue-ev(21))**(ev(25)-1.0_cp)      ! correct if lowP phase is highsym
                  if (nint(ev(20)) ==1) dVs=-1.0_cp*dVs                     ! change sign for highp phase is highsym
               else                                                         ! in the high phase
                  dVs=ev(26)*ev(27)*abs(Ptrue-ev(21))**(ev(27)-1)           ! correct if highP phase is highsym
                  if (nint(ev(20)) /= 1) dVs=-1.0_cp*dVs
               end if

            case (2,3) ! Landau VT or PVT
               if (transition_phase(Ptrue,T,eospar) ) then
                  !> low phase
                  Ttr = get_transition_temperature(Ptrue,eospar)
                  if (abs(Ttr-T) < 0.0001) then
                     dvs=0.0_cp
                  else
                     dVs=ev(24)*ev(25)*(abs(Ttr-T))**(ev(25)-1)*(ev(22)+2.0_cp*ev(23)*Ptrue) ! for curved phase boundary
                  end if
                  if (nint(ev(20)) /= 1) dVs=-1.0_cp*dVs                              ! sign change when low phase is highT
               else
                  !> in the highphase
                  Ttr = get_transition_temperature(Ptrue,eospar)
                  if (abs(Ttr-T) < 0.0001) then
                     dvs=0.0_cp
                  else
                     dVs=ev(26)*ev(27)*(abs(Ttr-T))**(ev(27)-1)*(ev(22)+2.0_cp*ev(23)*Ptrue)   !  curved boundary
                  end if
                  if (nint(ev(20)) ==1) dVs=-1.0_cp*dVs                                 ! sign change when high phase is highT
               end if
         end select

         !>Apply the correction to Kc due to dVs/dP:
         Vs=get_transition_strain(Ptrue,T,eospar) ! Vs  if volume eos, linear strain if linear eos
         kc=1.0_cp/(1.0_cp/kc - dVs/(1+Vs))
      end if
   End Function K_Cal

   !!----
   !!---- KP_CAL
   !!----    Returns value of kprime at this volume for EoS
   !!----
   !!--..    Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- 07/02/2017
   !!
   Module Function Kp_Cal(V,T,EoSpar,P) Result (kpc)
      !---- Arguments ----!
      real(kind=cp),           intent(in) :: V       ! Volume
      real(kind=cp),           intent(in) :: T       ! Temperature
      type(Eos_Type),          intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp), optional, intent(in) :: P       ! Pressure if known
      real(kind=cp)                       :: kpc

      !---- Local Variables ----!
      integer                            :: j
      type(Eos_Type)                     :: eosbare
      real(kind=cp)                      :: vv0,k0,kp,kpp,ptr,vol,vs,ttr,betap_tran,betap_bare,k
      real(kind=cp)                      :: a,b,f,rkp_top, rkp_bot,nu,dkdv,vt
      real(kind=cp)                      :: c,c0,c2,da,db,dc,d2a,d2b,d2c,group,dgroup,x,PFG0,kc   ! APL
      real(kind=cp), dimension(N_EOSPAR) :: ev
      real(kind=cp), dimension(3)        :: abc     ! Tait parameters
      real(kind=cp)                      :: dVsdP,d2VsdP2
      real(kind=cp),dimension(-2:2)      :: kpt
      real(kind=cp)                      :: delv
      real(Kind=cp)                      :: Ptrue    ! Input pressure if present. The true pressure without Pthermal subtracted
      real(Kind=cp)                      :: Pcorr    ! Pressure minus Pthermal.

      !> Init
      kpc=0.0_cp

      !> Pressure is needed for Tait, Murngahan, and for spontaneous strain at transitions
      !> This is the 'true' pressure, uncorrected for Pthermal
      !> Using Pin avoids rounding errors building up when get_pressure is called, with transition eos.
      if (present(P) ) then
         Ptrue=P
      else
         Ptrue=get_pressure(v,t,eospar)
      end if

      !> PTV table
      if (eospar%imodel == -1) then
         if (.not. present(P) ) Ptrue=get_props_ptvtable(0.,T,V,eospar,'P')
         kpc=get_props_ptvtable(Ptrue,T,V,eospar,'KP')
         return
      end if

      !> Isothermal EoS: get K0 and Kp0 at this T, then do algebra for normal EoS
      !>Or pthermal at Tref and do additional d(Pth)/dV contribution afterwards
      !>If transition, do calculation for bare phase, and add correction for transition afterwards

      !> Correct the volume to the high phase only if there is a transition
      if (eospar%itran > 0) then
         Vs=get_transition_strain(Ptrue,T,eospar)     ! returns the linear or volume strain as appropriate
         Vol=v/(1.0+vs)
      else
         Vol=v
      end if

      select case (eospar%itherm)
         case (0,6,7)                           ! No thermal model, or we have pthermal, so need params at Tref
            vv0=vol/eospar%params(1)            ! vv0 or aa0
            k0=eospar%params(2)
            kp=eospar%params(3)
            kpp=eospar%params(4)

         case (1:5)
            vv0=vol/get_V0_T(t,eospar)          ! vv0 is  v(p,t)/v(p=0,t): in transition, highphase
            k0=Get_K0_T(T,eospar)               ! returns M(T) for linear,
            if (err_CFML%IErr == 1) return                 ! exit with value eosparms(2) if k0 at P=0 calculated as negative

            kp=Get_Kp0_T(T,eospar)
            kpp=Get_Kpp0_T(T,eospar)
      end select

      !> Strain for BM, NS, Vinet EoS equations
      f=strain(vv0,eospar)
      if (err_CFML%IErr==1) then
         err_CFML%Msg=trim(err_CFML%Msg)//' called from Kpcal'
         return
      end if

      !>adjust for linear
      if (eospar%linear) then
         vol=vol**3.0_cp
         vv0=vv0**3.0_cp
         k0=k0/3.0_cp
         kp=kp/3.0_cp
         kpp=kpp/3.0_cp
      end if

      !>Now do the calculation for an isothermal eos at this T
      select case(eospar%imodel)
         case (1) ! Murnaghan
            kpc=kp

         case (2) ! Birch-Murnaghan
            a=0.0_cp
            b=0.0_cp
            if (eospar%iorder > 2) a=1.5_cp*(kp-4.0_cp)
            if (eospar%iorder ==4) b = (9.0_cp*k0*kpp + (9.0_cp*kp*kp) - 63.0_cp*kp + 143.0_cp)/6.0_cp

            rkp_top= 12.0_cp + 2.0_cp*a + (49.0_cp+ 32.0_cp*a + 6.0_cp*b)*f + &
                     (81.0_cp*a + 60.0_cp*b)*f*f + 121.0_cp*b*f*f*f
            rkp_bot= 3.0_cp*(1.0_cp + (7.0_cp + 2.0_cp*a)*f + (7.0_cp*a + 3.0_cp*b)*f*f + 7.0_cp*b*f*f*f)
            kpc=rkp_top/rkp_bot

         case (3) ! Vinet: expression derived RJA 26-Feb-2013. Use new definition of f: RJA 28/03/2013
            nu=1.5_cp*(kp-1.0_cp)
            a= (1.0_cp + (1.0_cp + nu)*f - nu*f*f)
            kpc= (2.0_cp + nu*(1.0_cp-f) + (1.0_cp- f)*(1.0_cp + nu - 2.0_cp*nu*f)/a)/3.0_cp

         case (4) ! Natural strain: expanded from Poirier and Tarantola PEPI 109:1-8, derived RJA 26-Feb-2013
            a=0.0_cp
            b=0.0_cp
            if (eospar%iorder > 2) a=1.5_cp*(kp-2.0_cp)
            if (eospar%iorder ==4) b=1.5_cp*(1.0_cp + k0*kpp + (kp-2.0_cp) + (kp-2.0_cp)**2.0_cp)
            rkp_top=  1.0_cp + 2.0_cp/3.0_cp*a + 2.0_cp*(a+b)*f + 3.0_cp*b*f*f
            rkp_bot=  1.0_cp + (3.0_cp +2.0_cp*a)*f + 3.0_cp*(a+b)*f*f +3.0_cp*b*f*f*f
            kpc=1.0_cp+rkp_top/rkp_bot

          case(5) ! Tait
            abc=get_tait(eospar,t)              ! get_tait returns volume-like parameters even for linear
            Pcorr=Ptrue
            if (eospar%pthermaleos) pcorr=ptrue - pthermal(Vol,T,eospar)   ! correct P for pthermal, needed only for Tait
            if (abc(2)*Pcorr < -0.999999_cp) then
               err_CFML%IErr=1
               err_CFML%Msg='Pressure yields infinite Kp for Tait Eos'
               kpc=9999.0        ! safe value return
            else
               kpc=(kp+1.0_cp)*((1.0_cp + abc(2)*Pcorr)**abc(3)*(1.0_cp-abc(1)) + abc(1)) -1.0_cp
            end if

         case(6) ! APL
            x=vv0**0.333333_cp
            pFG0=AFERMIGAS*(eospar%params(4)*vv0/vol)**1.66666667_cp
            c0=-1.0_cp*log(3.0_cp*K0/pFG0)         ! assumes V in A^3
            c2=1.5_cp*(Kp-3.0_cp)-c0

            !> terms in pressure expression
            a=1.0_cp/x**5.0_cp*(1.0_cp-x)
            b=exp(c0*(1.0_cp-x))
            c=1.0_cp+c2*x*(1.0_cp-x)

            !> derivatives
            da= -5.0_cp/x**6.0_cp +4.0_cp/x**5.0_cp
            db= -1.0_cp*c0*b
            dc= c2*(1.0_cp-2.0_cp*x)
            d2a= 30.0_cp/x**7.0_cp - 20.0_cp/x**6.0_cp
            d2b=c0*c0*b
            d2c=-2.0_cp*c2

            !> collect terms: checked ok 25-Feb-2017
            group=(da*b*c+a*db*c+a*b*dc)
            dgroup=a*b*d2c + a*d2b*c + d2a*b*c + 2.0_cp*(a*db*dc + b*da*dc + c*da*db)
            kpc= -1.0_cp/3.0_cp - x*dgroup/3.0_cp/group
      end select

      !> To this point Kpc is the bulk modulus of the 'bare' eos of the high phase if it was an isothermal eos
      !> Or for thermal-pressure EoS it is the bulk modulus at Tref and V

      !> Now correct thermal-pressure EoS for d(Pth)/dV contribution to bulk modulus, when possible
      if (eospar%Pthermaleos) then
         select case(eospar%itherm)
            case(7)           !MGD EoS: Do this numerically on complete EoS without transition
               eosbare=eospar
               eosbare%itran=0    ! clear any transition terms
               delv=0.01_cp*vol
               do j=-2,2,1             ! loop calculates dK/dV
                  vt=vol+real(j)*delv  ! apply shift to v
                  kpt(j)=k_cal(vt,t,eosbare)     ! calc resulting k
               end do
               dkdV=(Kpt(-2)+8.0_cp*(Kpt(1)-Kpt(-1))-Kpt(2))/(12.0_cp*delv)
               kpc= -1.0_cp*dkdV*vol/kpt(0)

            case default      !includes HP thermal-pressure which has dK = 0 along isochor


         end select
      end if

      !> Now correct from volume to linear if required, because transition code works in linear
      if (eospar%linear) kpc=kpc*3.0_cp

      !> And the transition effects
      if (eospar%itran > 0 ) then
         !> Local EoS copy, but no scaling of transition parameters for linear
         call EoS_to_Vec(eospar,ev)

         select case(eospar%itran)
            case(1) ! Landau P-V
               if (transition_phase(Ptrue,T,eospar)) then             ! in the low phase
                  dVsdP=ev(24)*ev(25)*abs(Ptrue-ev(21))**(ev(25)-1)   ! d(Vs)/dP correct if lowP phase is highsym
                  if (nint(ev(20)) /=1) dVsdP=-1.0_cp*dVsdP           ! change sign for lowp phase is highsym
                  d2VsdP2=dVsdP*(ev(25)-1)/abs(Ptrue-ev(21))
                  if (nint(ev(20)) /=1) d2VsdP2=-1.0_cp*d2VsdP2
               else                                                                      ! in the high phase
                  dVsdP=ev(26)*ev(27)*abs(Ptrue-ev(21))**(ev(27)-1)     ! correct if highP phase is highsym
                  if (nint(ev(20)) == 1) dVsdP=-1.0_cp*dVsdP
                  d2VsdP2=dVsdP*(ev(27)-1)/abs(Ptrue-ev(21))
                  if (nint(ev(20)) ==1) d2VsdP2=-1.0_cp*d2VsdP2
               end if

            case(2,3) ! Landau VT or PVT
               Ttr = get_transition_temperature(Ptrue,eospar)
               Ptr = get_transition_pressure(T,eospar)
               if (transition_phase(Ptrue,T,eospar)) then
                  if (abs(Ttr-T) < 0.0001) then
                     dVsdP  = 0._cp            ! on the boundary
                     d2VsdP2= 0._cp
                  else                          ! in the low phase
                     dVsdP=ev(24)*ev(25)*(abs(Ttr-T))**(ev(25)-1)*(ev(22)+2.0_cp*ev(23)*Ptrue) ! for curved phase boundary
                     if (nint(ev(20)) /= 1) dVsdP=-1.0_cp*dVsdP                              ! sign change when low phase is highT

                     d2VsdP2= dVsdP*(2.0_cp*ev(23))/(ev(22)+2.0_cp*Ptrue*ev(23))      ! This works 15-Feb
                     if (nint(ev(20)) == 1) then
                        d2VsdP2= d2VsdP2 + dVsdP*(ev(25)-1)/abs(Ttr-T)*(ev(22)+2.0_cp*Ptrue*ev(23))
                     else
                        d2VsdP2= d2VsdP2 - dVsdP*(ev(25)-1)/abs(Ttr-T)*(ev(22)+2.0_cp*Ptrue*ev(23))
                     end if
                  end if

               else
                  !> in the highphase
                  if (abs(Ttr-T) < 0.0001) then
                     dVsdP  = 0._cp
                     d2VsdP2= 0._cp
                  else
                     dVsdP=ev(26)*ev(27)*(abs(Ttr-T))**(ev(27)-1)*(ev(22)+2.0_cp*ev(23)*Ptrue)   !  curved boundary
                     if (nint(ev(20)) ==1) dVsdP=-1.0_cp*dVsdP                                 ! sign change when high phase is highT
                     d2VsdP2=dVsdP*(2.0_cp*ev(23))/(ev(22)+2.0_cp*Ptrue*ev(23))
                     if (nint(ev(20)) == 1) then
                        d2VsdP2=d2VsdP2-1.0_cp*dVsdP*(ev(27)-1)/abs(Ttr-T)*(ev(22)+2.0_cp*Ptrue*ev(23))  !high phase is highT
                     else
                        d2VsdP2=d2VsdP2+1.0_cp*dVsdP*(ev(27)-1)/abs(Ttr-T)*(ev(22)+2.0_cp*Ptrue*ev(23)) !high phase is lowT
                     end if
                  end if
               end if
         end select

         !> Now apply the transition effect
         Vs=get_transition_strain(Ptrue,T,eospar) ! Vs  if volume eos, linear strain if linear eos

         betap_tran=   d2VsdP2/(1.0_cp+Vs) -(dVsdP/(1.0_cp+Vs))**2.0_cp         ! transition contribution

         eosbare=eospar
         eosbare%itran=0     ! create a bare phase eos
         kc=k_cal(vol,t,eosbare,p=Ptrue)   ! k_cal returns linear moduli if linear
         k=k_cal(v,t,eospar,p=Ptrue)
         betap_bare= -1.0_cp*kpc/kc**2.0_cp
         kpc= (betap_tran - betap_bare)*(k**2.0_cp)
      end if
   End Function Kp_Cal

   !!----
   !!---- KPP_CAL
   !!----    Returns value of kdouble-prime at this volume for EoS
   !!----
   !!--..    Validated code against Eosfit v5.2 for non-thermal EoS: RJA 27/02/2013
   !!----
   !!---- 25/012/2016
   !!
   Module Function Kpp_Cal(V,T,EoSpar) Result (kppc)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: kppc

      !---- Local Variables ----!
      integer                      :: j
      real(kind=cp)                :: ptr,vlimit,plimit,vlimitk
      real(kind=cp)                :: delp,p0,p,vt
      real(kind=cp),dimension(-2:2):: kpt

      !> Init
      kppc=0.0_cp

      !> Trap for pvt table....do calculation later
      if(eospar%imodel == -1)return

      !> jump on equation type
      if (eospar%imodel == 1 .and. eospar%itran == 0) return ! Murnaghan without transition

      !> numerical solution
      ! To get reasonable Kprime numerically, need delP as 1 GPa for K0=160 GPa because this generates 1% change in Kp
      ! Therefore need as big del for Kpp, if not bigger

      delp=K_cal(V,T,eospar)/100.        !delp=K/100
      p0=get_pressure(V,T,eospar)         ! current p

      !> Code to prevent stepping across transition
      if (eospar%itran > 0) then
         Ptr=get_transition_pressure(t,eospar)
         if (transition_phase(P0+2.0*delp,T,eospar) .neqv. transition_phase(P0,T,eospar)) delp=0.2*abs(ptr-p0)
         if (transition_phase(P0-2.0*delp,T,eospar) .neqv. transition_phase(P0,T,eospar)) delp=0.2*abs(ptr-p0)
      end if
      !> Code to stop MGD EoS going into illegal large volume at negative delp
      if(eospar%itherm == 7)then
          vlimitk=get_volume_K(0._cp,t,eospar)
          vlimit=get_volume_K(eospar%params(2)/2.0_cp,eospar%tref,eospar)
          if(vlimitk > tiny(0.0_cp) .and. vlimitk < vlimit)vlimit=vlimitk
          plimit=get_pressure(vlimit,t,eospar)
          if(p0-2.0*delp < plimit)delp=0.2*abs(p0-plimit)
      endif

      do j=-2,2,1
         p=p0+real(j)*delp                 ! apply shift to p
         vt=get_volume(p,t,eospar)
         kpt(j)=kp_cal(vt,t,eospar,p=p)     ! calc resulting Kp
      end do

      kppc=(kpt(-2)+8.0_cp*(kpt(1)-kpt(-1))-kpt(2))/(12.0_cp*delp)     ! Derivative to second order approximation

      !> No linear conversion is required because kp_cal returns values for "linear Kp" = Mp,
      !> so kppc is already dMp/dP = Mpp

   End Function Kpp_Cal

End SubModule EoS_K_Cal