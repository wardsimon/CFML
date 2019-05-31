!!----
!!----
!!----
SubModule (CFML_EoS) EoS_002
   Contains
   
   !!--++
   !!--++ GET_V0_T
   !!--++
   !!--++ PRIVATE
   !!--++    Returns the volume at P=0 and T=T from V0 and Thermal expansion
   !!--++    Except for Pthermal, for which it returns V at P=0, T=Tref
   !!--++    It calculates V0 only for thermal expansion, not including transition effects!
   !!--++    Therfore this must remain PRIVATE
   !!--++
   !!--++ 17/07/2015
   !!
   Module Function Get_V0_T(T,EosPar) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T        ! Temperature
      type(Eos_Type), intent(in) :: EoSPar   ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: V,kp,AK
      real(kind=cp)                      :: Tref,A,B,C,Tn,tt
      real(kind=cp)                      :: delt,delt2
      real(kind=cp), dimension(n_eospar) :: ev

      !> Init
      Tref=eospar%tref

      !> Local copy Eospar
      call EoS_to_Vec(eospar,ev) ! Volume or linear case is covered
                                 ! all equations written for volume
      delt=T-Tref
      select case(eospar%itherm)
         case(0)
            v=ev(1)                 ! no thermal eos: V is V0

         case(1)                    ! Berman, works at all T
            V=ev(1)*(1.0_cp +ev(10)*delt+0.5_cp*ev(11)*delt*delt)

         case(2)                    ! Fei,
            TT=T
            delt2=t*t-tref*tref
            if (tt < tiny(0._cp)) tt=0.001        ! to prevent divide by zero

            if (abs(ev(12)) > tiny(0._cp)) then
               V=ev(1)*exp(ev(10)*delt + 0.5_cp*ev(11)*delt2 - ev(12)*(1.0_cp/TT - 1.0_cp/Tref))
            else
               V=ev(1)*exp(ev(10)*delt + 0.5_cp*ev(11)*delt2)    ! to protect from divide by zero when alpha2=0 and T=0
            end if

         case(3)                    ! HP 1998 modified
            tt=t
            if (t < tiny(0._cp)) tt=tiny(0._cp)      ! prevents sqrt(-ve T)
            V=ev(1)*(1.0_cp + ev(10)*(tt-tref) - 2.0_cp*(10.0_cp*ev(10)+ev(11))*(sqrt(tt) -sqrt(tref)))

         case(4)                    ! Holland-Powell 2011 in the Kroll form
            !>>>>>kp=ev(3)   : version before 11/11/2016
            if (eospar%icross == 2) then
               kp=ev(5)
            else
               kp=ev(3)
            end if
            if (abs(kp-1) < 0.0001 .or. abs(kp/(kp+2.0_cp)) < 0.0001)then
               V=ev(1)                             ! In these cases algebra shows V=V0
            else
               Tn= ev(11)/Tref                        ! theta/Tref
               C=Tn*Tn*exp(Tn)/(exp(tn)-1.0_cp)**2.0_cp
               B=-1.0/kp/(kp+2.0_cp)
               if (t > 0.05_cp*ev(11)) then                               ! to avoid numerical problems at T=0
                  A=ev(10)*ev(11)/C *(1.0_cp/(exp(ev(11)/T)-1.0_cp) - 1.0_cp/(exp(Tn)-1.0_cp) )
               else
                  A=ev(10)*ev(11)/C *(-1.0_cp/(exp(Tn)-1.0_cp) )          ! because when T=0 then 1.0/exp(Tein/T) = 0
               end if
               ! V=ev(1)*(-1.0_cp*kp + (1.0_cp+kp)*(1.0_cp - kp*(kp+2.0_cp)*A/(kp+1.0_cp))**B)
               AK=1.0_cp - kp*(kp+2.0_cp)*A/(kp+1.0_cp)
               if (AK < tiny(0.0_cp) )then
                  V=ev(1)       ! for safe return
                  err_CFML%Ierr=1
                  err_CFML%Msg='T exceeds valid limit for Kroll expansion in get_V0_T'
               else
                  V=ev(1)*(-1.0_cp*kp + (1.0_cp+kp)*AK**B)
               end if
            end if

         case(5)                    ! Salje, Tref fixed at zero
            A=T/ev(11)
            if (A < 0.001) then
               V=ev(1)                 ! ultra-low T: coth(theta_sat/T)=1.0
            else
               A=1.0_cp/tanh(1.0_cp/A) ! the coth(theta_sat/T)
               V=(ev(1)**(1.0_cp/3.0_cp) + ev(10)*ev(11)*(A-1.0_cp))**3.0_cp
            end if

         case(6,7)
            v=ev(1)         ! Pthermal needs V0 at Tref
      end select

      !> Linear
      if (eospar%linear) v=v**(1.0_cp/3.0_cp)
   End Function Get_V0_T

   !!----
   !!---- GET_VOLUME
   !!----    Find volume from EoS at given P and T
   !!----
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 25/02/2013
   !!----
   !!---- 16/02/13
   !!
   Module Function Get_Volume(P,T,EosPar) Result(v)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp), parameter          :: PREC=0.000001_cp  !precision to find V.

      integer                           :: nstep
      type(Eos_Type)                    :: EoS  ! Eos copy
      real(kind=cp)                     :: V
      real(kind=cp)                     :: V0,K0,Kp,k,strain,vfactor
      real(kind=cp)                     :: Vol, step,dp1,dp2
      real(kind=cp),dimension(N_EOSPAR) :: ev
      real(kind=cp),dimension(3)        :: abc          ! Tait parameters
      real(kind=cp)                     :: pa           ! pa=p-pth

      !> Init
      v=0.0_cp
      strain=0.0_cp   ! strain from transition: linear or volume to match eos type
      pa=p            ! local copy p

      !> If PTV table, go directly
      if (eospar%imodel == -1) then
         v=get_props_ptvtable(pa,t,0.0,eospar,'V')     ! get_props_ptvtable returns length if linear
         return
      end if

      !> Local copy Eospar
      call EoS_to_Vec(eospar,ev) ! Volume or linear case is covered

      !> Set appropriate V0, and adjust for thermal pressure:
      select case (eospar%itherm)
         case (0)                                   ! 0=no thermal,
            v0=ev(1)                                  ! v0 is volume eos, (a0)^3 for linear

         case (1:5)
            v0=get_v0_t(t,eospar)                     ! returns a0 for linear
            if (eospar%linear) v0=v0**3.0_cp

         case (6)                                     ! HP thermal pressure
            v0=ev(1)
            pa=p-pthermal(0.0,t,eospar)               ! adjust pressure to isothermal pressure for murn and tait estimates
            
         case(7)                                    ! MGD - do a guess on basis parallel isochors of pa
            v0=ev(1)
            pa=p - eospar%params(2)*(t-eospar%tref)/100000.         ! have to guess an alpha because that requires V !!!  
      end select

      !> set K0, for various purposes
      k0=Get_K0_T(T,eospar)              ! Handles thermal pressure case, returns K0 or M0
      if (eospar%linear) k0=k0/3.0_cp

      kp=Get_Kp0_T(T,eospar)
      if (eospar%linear) kp=kp/3.0_cp

      !> Get the volume strain due to transition: only a function of P,T NOT V!!
      if (eospar%itran > 0) then
         strain=get_transition_strain(P,T,eospar)     ! returns the linear or volume strain as appropriate
      end if

      !> If there is no eos model, we are finished because V0 is the V at this T
      if (eospar%imodel == 0) then
         if (eospar%linear) v0=v0**(1.0_cp/3.0_cp)
         v=v0*(1.0_cp + strain)
         return
      end if

      !> Analytic solution for Murnaghan:  use this for first guess for other EoS except Tait
      vfactor=(1.0_cp + kp*pa/k0)
  
      !if (vfactor < 0.0 .and. kp > 0.0)then
      if (vfactor < 0.0)then
         v=v0        ! safe value for when 1+kp*pa/k0 is negative
      else
         v=v0*(vfactor)**(-1.0_cp/kp)
      end if

      !> Cannot do the following if MGD pthermal
      if (eospar%itherm /=7) then
         if (eospar%imodel ==1) then
            !> Exact solution for Murnaghan
            if (eospar%linear) v=v**(1.0_cp/3.0_cp)
            if (eospar%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
            !if (vfactor < 0.0 .and. kp > 0.0) then
            if (vfactor < 0.0) then
               err_CFML%Ierr=1
               err_CFML%Msg='Pressure < -K0/Kp: yields meaningless volumes for Murnaghhan EoS'
            end if
            return
         end if

         !> Analytic solution for Tait
         if (eospar%imodel ==5) then
            call get_tait(eospar,t,abc)                     ! get_tait returns volume-like parameters even for linear
            if (abc(2)*pa < -0.999999_cp) then
               err_CFML%Ierr=1
               err_CFML%Msg='Pressure yields infinite volume for Tait Eos'
               v=9999.0        ! safe value return
            else
               v=v0*(1.0_cp-abc(1)*(1.0_cp-(1.0_cp + abc(2)*pa)**(-1.0_cp*abc(3))))
               if (eospar%linear) v=v**(1.0_cp/3.0_cp)
               if (eospar%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
            end if
            return
         end if
      end if

      !> Find iterative solution for the rest of functions: get_pressure includes the thermal pressure term
      !> But if there is a transition, we only want the P/V for the bare high-symm phase without softening
      vol=v
      if (eospar%linear) vol=vol**(1.0_cp/3.0_cp)
      eos=eospar        ! copy
      eos%itran=0       ! turn off transition
      
      dp1=p-get_pressure(vol,t,eos)
      if (err_CFML%Ierr==1)then
         if(eospar%linear)v=v**(1.0_cp/3.0_cp)           ! to ensure reasonable return value if linear
         return
      end if

      !> estimate the step to make in vol to get closer
      k=k0+p*kp                          ! Murn-like guess estimate to avoid recursive call back here when pthermal used
      if (eospar%linear)k=k*3.0_cp       ! The iteration is working with linear quantities
      if(abs(k) > abs(dp1/1000.))then
          step= -1.0_cp*vol*dp1/k            ! by definition of bulk modulus
      else
          step = -1.0_cp*vol*dp1/10.        ! trap for k being close to zero, to prevent step becoming NaN
      endif
      if(abs(step) < eos%params(1)/1000.)step=eos%params(1)/1000.

      nstep=0
      iter: do
         !> Trap infinite loops
         if (nstep > 10000)then
            err_CFML%Ierr=2
            err_CFML%Msg='No convergence in Get_Volume'
            exit iter
         endif

         !> Increment Volume
         vol=vol+step
         dp2=p-get_pressure(vol,t,eos)
         if (err_CFML%Ierr==1)exit iter
         nstep=nstep+1

         !> test for sufficient convergence
         !if (abs(step) < PREC*Vol)exit  ! 1 part in 100,000 in volume
         if (abs(dp2) < abs(prec*k0))exit    ! If dP2 very small, finished; from dP= -K.dV/V and prec=dV/V
         if (abs(step) < PREC*Vol) then  ! 1 part in 100,000 in volume
            if (abs(abs(dP2)-abs(10.0*step*k0/vol)) > abs(max(0.005,p/1000.)))then   ! was 0.005,1000
               err_CFML%Ierr=2
               write(err_CFML%Msg,'(''Convergence problem in Get_Volume: dP = '',f6.3,''>> K0*dV/V'')')dP2
            end if
            exit iter
         end if

         !> not converged, so adjust step size
         if (dp1*dp2 < 0.0_cp) then
            !> overshot ptr:reverse step direction and make size smaller
            step=-0.5_cp*step
         else
            if (abs(dp2) > abs(dp1)) then
               step=-1.0*step      ! wrong direction: reverse
            else
        !       step=0.9_cp*dp2/dp1*step   ! correct direction, should get a smaller step size  
        !    end if                        ! the factor of 0.9 is a safety factor to ensure step gets smaller
                step= -1.0*dp2*step/(dp2-dp1)    ! new version Dec 2018: adjust step size in Newton-Raphson manner
            end if
         end if

         dp1=dp2        ! update delta-p values and go back for next cycle
      end do iter
      
      !> now set return value depending on success or not
      v=vol

      if (eospar%itran > 0) v=vol*(1.0_cp + strain)  ! apply transition strain ('vol' is actually linear if linear eos)
   End Function Get_Volume
   
   !!----
   !!---- GET_VOLUME
   !!----    Find volume from EoS at given P and T
   !!----
   !!--.. Validated code against Eosfit v5.2 for non-thermal EoS: RJA 25/02/2013
   !!----
   !!---- 16/02/13
   !!
   Module Function Get_Volume_New(P,T,EosPar) Result(v)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp), parameter          :: PREC=0.000001_cp  !precision to find V.

      integer,parameter                 :: nstep=30
      integer                           :: i,ic
      real(kind=cp),dimension(nstep):: x,y,d2y

      
      type(Eos_Type)                    :: EoS  ! Eos copy
      real(kind=cp)                     :: V
      real(kind=cp)                     :: V0,K0,Kp,k,strain,vfactor
      real(kind=cp)                     :: Vol, vstep, delp_prev,delp,v_prev
      real(kind=cp),dimension(N_EOSPAR) :: ev
      real(kind=cp),dimension(3)        :: abc          ! Tait parameters
      real(kind=cp)                     :: pa           ! pa=p-pth

      logical                           :: reverse

      !> Init
      v=0.0_cp
      strain=0.0_cp   ! strain from transition: linear or volume to match eos type
      pa=p            ! local copy p

      !> If PTV table, go directly
      if (eospar%imodel == -1) then
         v=get_props_ptvtable(pa,t,0.0,eospar,'V')     ! get_props_ptvtable returns length if linear
         return
      end if

      !> Local copy Eospar
      call EoS_to_Vec(eospar,ev) ! Volume or linear case is covered

      !> Set appropriate V0, and adjust for thermal pressure:
      select case (eospar%itherm)
         case (0)                                   ! 0=no thermal,
            v0=ev(1)                                  ! v0 is volume eos, (a0)^3 for linear

         case (1:5)
            v0=get_v0_t(t,eospar)                     ! returns a0 for linear
            if (eospar%linear) v0=v0**3.0_cp

         case (6)                                     ! HP thermal pressure
            v0=ev(1)
            pa=p-pthermal(0.0,t,eospar)               ! adjust pressure to isothermal pressure for murn and tait estimates
            
         case(7)                                    ! MGD - do a guess on basis parallel isochors of pa
            v0=ev(1)
            pa=p - eospar%params(2)*(t-eospar%tref)/100000.         ! have to guess an alpha because that requires V !!!  
      end select

      !> set K0, for various purposes
      k0=Get_K0_T(T,eospar)              ! Handles thermal pressure case, returns K0 or M0
      if (eospar%linear) k0=k0/3.0_cp

      kp=Get_Kp0_T(T,eospar)
      if (eospar%linear) kp=kp/3.0_cp

      !> Get the volume strain due to transition: only a function of P,T NOT V!!
      if (eospar%itran > 0) then
         strain=get_transition_strain(P,T,eospar)     ! returns the linear or volume strain as appropriate
      end if

      !> If there is no eos model, we are finished because V0 is the V at this T
      if (eospar%imodel == 0) then
         if (eospar%linear) v0=v0**(1.0_cp/3.0_cp)
         v=v0*(1.0_cp + strain)
         return
      end if

      !> Analytic solution for Murnaghan:  use this for first guess for other EoS except Tait
      vfactor=(1.0_cp + kp*pa/k0)
      if (vfactor < 0.0)then
         v=v0        ! safe value for when 1+kp*pa/k0 is negative
      else
         v=v0*(vfactor)**(-1.0_cp/kp)
      end if

      !> Cannot do the following if MGD pthermal
      if (eospar%itherm /=7) then
         if (eospar%imodel ==1) then
         !> Exact solution for Murnaghan
            if (eospar%linear) v=v**(1.0_cp/3.0_cp)
            if (eospar%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
            !if (vfactor < 0.0 .and. kp > 0.0) then
            if (vfactor < 0.0) then
               err_CFML%IErr=1
               err_CFML%Msg='Pressure < -K0/Kp: yields meaningless volumes for Murnaghhan EoS'
            end if
            return
         end if

         !> Analytic solution for Tait
         if (eospar%imodel ==5) then
            call get_tait(eospar,t,abc)                     ! get_tait returns volume-like parameters even for linear
            if (abc(2)*pa < -0.999999_cp) then
               err_CFML%IErr=1
               err_CFML%Msg='Pressure yields infinite volume for Tait Eos'
               v=9999.0        ! safe value return
            else
               v=v0*(1.0_cp-abc(1)*(1.0_cp-(1.0_cp + abc(2)*pa)**(-1.0_cp*abc(3))))
               if (eospar%linear) v=v**(1.0_cp/3.0_cp)
               if (eospar%itran > 0) v=v*(1.0_cp + strain)     ! apply transition strain (already converted if linear)
            end if
            return
         end if
      end if

      !> Find iterative solution for the rest of functions: get_pressure includes the thermal pressure term
      !> But if there is a transition, we only want the P/V for the bare high-symm phase without softening
      
      !From here work with Vol, starting from Murnaghan estimate
      vol=v
      if (eospar%linear) vol=vol**(1.0_cp/3.0_cp)
      eos=eospar        ! copy
      eos%itran=0       ! turn off transition          
      
      !> initial simple hunt
      delp_prev=huge(0._cp)
      ic = 0 
      reverse=.false.
      Vstep=eos%params(1)/100._cp
      do 
         ic=ic+1
         if (ic > 1000)then
            err_CFML%Ierr=1
            err_CFML%Msg=' *****No solution found in get_volume after 1000 cycles'
            return
         end if
        
         !> ! have to clear the previous errors, otherwise get_pressure will return 0   
         call clear_error()
         delp=p-get_pressure(Vol,T,eos) 
          
         if (delp*delp_prev < 0._cp .and. ic > 1)then     ! over-stepped solution: solution between v_prev and v 
            vol=vol-delp*Vstep/(delp-delp_prev)                               ! best guess
            exit
         end if
   
         if (abs(delp) > abs(delp_prev))then ! delta-pressure getting bigger
            if (reverse)then               ! found a minimum between v_prev-vstep and v 
               err_CFML%Ierr=1
               err_CFML%Msg=' *****No volume found in get_volume'
               return
            else
               reverse=.true.            ! just going the wrong way
               vstep=-1.0_cp*vstep
            end if
         end if        
         v_prev=vol         ! this volume stored
         delp_prev=delp  ! store delp
         Vol=Vol+Vstep
      end do 
      
      !> now calculate PV around the solution: we want increasing P, so this means vstep < 0
      Vstep=-1.0_cp*abs(Vstep)
      Vol=Vol-2.0*Vstep
      Vstep=4.0*Vstep/nstep 

      do i=1,nstep
         x(i)=get_pressure(vol,t,eos)
         y(i)=vol
         vol=vol+vstep
      end do
      d2y=Second_Derivative(x, y, nstep)
      vol=spline_interpol(x,y,d2y,nstep,p)
      
      v=vol
      if (eospar%itran > 0) v=vol*(1.0_cp + strain)  ! apply transition strain ('vol' is actually linear if linear eos)
   End Function Get_Volume_New

   !!----
   !!---- GET_VOLUME_K
   !!----    Returns the value of Volume for a given K  at T without using pressure
   !!----    This has limited precision when Kp is small, so do not use except to obatin 
   !!----    approximate V (eg for limits to eos)
   !!----
   !!---- 06/03/2019
   !!
   Module Function Get_Volume_K(K,T,E) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: K       ! Bulk modulus
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: E      ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: V,vprev,Kprev,kc,vnew,delv,delvprev,Vstep

      !> Init
      v=0.0_cp
      
      Vprev=e%params(1)             !This is Vo
      Kprev=K_cal(Vprev,T,E)        !This is Ko

      !> Initial search
      Vstep=0.001_cp*e%params(1)            !If K is smaller than Ko, must go up in volume
      if (K > Kprev)Vstep=-1.0_cp*Vstep

      do 
         Vprev=Vprev+Vstep
         kc=K_cal(Vprev,T,E)
         if (K < Kprev) then                   !going to volumes > 1.0
            if (Kc < K) exit
            if (Vprev > 2.0_cp*e%params(1))then
               V=2.0_cp*e%params(1)            !stop infinite looping
               return
            end if
         else
            if (Kc > K)exit
            if (Vprev < 0.001_cp*e%params(1))then
               V=0.001_cp*e%params(1)      !stop infinite looping
               return
            end if
         end if
      end do
      
      !> set-up for Newton-raphson
      Kprev=Kc
      V=Vprev-Vstep

      do     ! does a newton-raphson search
         call clear_error()
         kc=K_cal(V,T,E)
         if (abs(kc-k) < 0.001)exit
         delV= (k-kc)*(V-Vprev)/(kc-Kprev)
         if (abs(delV) > abs(delVprev))delV=sign(delVprev,delV)       !prevents step getting bigger
         Vnew= V + delV
         if (Vnew < 0._cp)Vnew=0.99*V          ! stops V going negative
         Kprev=Kc
         Vprev=V
         delVprev=delV
         V=vnew
      end do
      
      !> Linear case
      if (e%linear) v=v**(1.0_cp/3.0_cp)
   End Function Get_Volume_K
   
   !!----
   !!---- GET_VOLUME_K_OLD
   !!----
   !!----
   !!
   Module Function Get_Volume_K_old(K,T,E) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: K       ! Bulk modulus
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: E       ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: V,vprev,Kprev,kc,vnew,delv,delvprev

      !> Init
      v=0.0_cp
      Vprev=e%params(1)             !This is Vo
      Kprev=K_cal(Vprev,T,E)        !This is Ko

      V=1.01_cp*e%params(1)
      
      do     ! does a newton-raphson search
         call clear_error()
         kc=K_cal(V,T,E)
         if (abs(kc-k) < 0.001)exit
         
         delV= (k-kc)*(V-Vprev)/(kc-Kprev)
         if (abs(delV) > abs(delVprev)) delV=sign(delVprev,delV)       !prevents step getting bigger
         Vnew= V + delV
         if (Vnew < 0.0_cp) Vnew=0.99*V          ! stops V going negative
         Kprev=Kc
         Vprev=V
         delVprev=delV
         V=vnew
      end do
     
      !> Linear case
      if (e%linear) v=v**(1.0_cp/3.0_cp)
   End Function Get_Volume_K_old

   !!----
   !!---- GET_VOLUME_S
   !!----    Returns the value of Volume obtained from Strain (S)
   !!----
   !!---- 15/02/2013
   !!
   Module Function Get_Volume_S(S,T,Eospar) Result(V)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: S       ! Strain
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp)                      :: V,v0
      real(kind=cp), dimension(n_eospar) :: ev

      !> Init
      v=0.0_cp

      !> Local Eos Copy
      call EoS_to_Vec(eospar,ev)

      !> Allow for thermal: s=function of V(P,T)/V(P=0,T)
      select case (eospar%itherm)
         case (0)
            v0=ev(1)

         case (1:7)
            v0=get_volume(0.0_cp,t,eospar)
            if (eospar%linear) v0=v0**3.0_cp
      end select

      select case (eospar%imodel)
         case (1,5) ! Murnaghan and Tait, no strain defined
            v=0.0_cp

         case (2) ! Birch-Murnaghan
            V=v0*(1.0_cp+2.0_cp*s)**(-1.5_cp)

         case (3) ! Vinet
            V=v0*(1.0_cp-s)**3.0_cp

         case (4) ! Natural Strain
            V=v0*exp(-3.0_cp*s)
      end select

      !> Linear case
      if (eospar%linear) v=v**(1.0_cp/3.0_cp)
   End Function Get_Volume_S
   
End SubModule (CFML_EoS) EoS_002