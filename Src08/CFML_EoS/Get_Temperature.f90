!!----
!!----
!!----
SubModule (CFML_EoS) EoS_004
   Contains
   !!----
   !!---- GET_TEMPERATURE
   !!----    Returns Temperature at given P,V
   !!----
   !!---- 01/04/2014
   !!
   Module Function Get_Temperature(P,V,EosPar) Result(T)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: P       ! Pressure
      real(kind=cp),  intent(in) :: V       ! Volume
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp)              :: t

      !---- Local Variables ----!
      integer                           :: nstep
      real(kind=cp)                     :: pa,va,ta
      real(kind=cp)                     :: step,dp1,dp2

      !> Init
      t=eospar%tref
      pa=p               ! local copy p
      va=v

      !> Check
      if (eospar%itherm == 0) return          ! no calcs possible without thermal model

      !> First estimate at P=0
      t=Get_Temperature_P0(va,EosPar)
      if (eospar%imodel ==0) return
      if (err_CFML%IErr ==1 ) return

      !> Use iterative solution, by minimising p-pcalc
      !> Init this part
      ta=t
      dp1=p-get_pressure(va,ta,eospar)
      nstep=0
      step=100.0_cp  ! For positive thermal expansion, should increase T from P=0

      !> Do iteration
      do
         !> Trap infinite loops
         if (nstep > 10000) then
            err_CFML%IErr=1
            write(err_CFML%Msg,'("No convergence in get_temperature after ",i5," steps, dp= ",f10.4)')nstep,dp1
            exit
         end if

         !> Increment Temperature
         ta=ta+step
         dp2=p-get_pressure(va,ta,eospar)
         nstep=nstep+1

         !> test for sufficient convergence
         if (abs(step) < 1) exit          ! Go to 1 K precision

         !> not converged, so adjust step size
         if (dp1*dp2 < 0.0_cp) then
            !> overshot ptr:reverse step direction and make size smaller
            step=-0.5_cp*step
         else
            if (abs(dp2) > abs(dp1))then
               step=-1.0_cp*step      ! wrong direction: reverse
            end if                 !(correct direction uses same step again)
         end if

         dp1=dp2        ! update delta-p values and go back for next cycle
      end do

      t=ta              ! success
   End Function Get_Temperature

   !!----
   !!---- GET_TEMPERATURE_P0
   !!----    Returns Temperature at P=0 for given V0T
   !!----
   !!---- 17/07/2015
   !!
   Module Function Get_Temperature_P0(Vin,EosPar,Tmin,Tmax) Result(Tk)
      !---- Arguments ----!
      real(kind=cp),           intent(in) :: Vin       ! Volume at temperature T or Pth (Pthermal case)
      type(Eos_Type),          intent(in) :: EoSPar    ! Eos Parameter
      real(kind=cp), optional, intent(in) :: Tmin      ! Range for solution in T
      real(kind=cp), optional, intent(in) :: Tmax
      real(kind=cp)                       :: tk

      !---- Local Variables ----!
      real(kind=cp)                      :: tref
      real(kind=cp)                      :: v00,v0T
      real(kind=cp)                      :: a,b,c,d,t1,t2,t3,x1,x2,x3
      real(kind=cp)                      :: a1,a2,a3,q,r,s1,s2,dd,th
      real(kind=cp)                      :: y,p
      real(kind=cp)                      :: kp,th_e,eps0
      real(kind=cp)                      :: Tkmin,Tkmax
      real(kind=cp), dimension(N_EOSPAR) :: ev

      !> Init
      Tk=eospar%tref

      !> Local copy Eospar to handle linear or volume
      call EoS_to_Vec(eospar,ev)

      !> all equations written for volume
      v00=ev(1)
      if (eospar%linear) then
         v0T=Vin**3.0_cp
      else
         v0T=Vin
      end if

      !> Optional arguments
      TKmin=0.0_cp
      if (present(Tmin)) TKmin=tmin

      TKmax=1.0e4
      if (present(Tmax)) TKmax=tmax

      !> Check
      Tref=eospar%tref
      if (v00 <= tiny(0.0)) return

      !> Init local variables
      t1=0.0
      t2=0.0
      t3=0.0

      select case (eospar%itherm)
         case (1) ! Berman
            ! modified RJA 31.03.2014 to fix case params(11)=0.
            if (abs(ev(11)) < tiny(0.0_cp)) then
               if (ev(10) > tiny(0.0_cp)) tk=Tref+(v0T/v00 -1.0_cp)/ev(10)
            else
               a= 0.5*ev(11)
               b= ev(10) - ev(11)*tref
               c= 1.0 - ev(10)*tref + 0.5*ev(11)*tref*tref - (v0t/v00)

               d= b*b - 4.0*a*c
               if (d >= 0.0) then
                  t1=(-b+sqrt(d))/(2.0*a)
                  t2=(-b-sqrt(d))/(2.0*a)
                  if (t1 > 0.0) then
                     if (t1 >=TKmin .and. t1 <=TKmax)then
                        tk=t1
                     else
                        err_CFML%IErr=1
                        err_CFML%Msg='No valid solution for temperature'
                     end if

                  else if (t2 > 0.0) then
                     if (t2 >=TKmin .and. t2 <=TKmax) then
                        tk=t2
                     else
                        err_CFML%IErr=1
                        err_CFML%Msg='No valid solution for temperature'
                     end if
                  else
                     err_CFML%IErr=1
                     err_CFML%Msg='No valid solution for temperature'
                  end if
               end if
            end if

         case (2) ! Fei
            !> modified RJA 31.03.2014 to fix case params(11)=0.
            if (abs(ev(11)) < tiny(0.0_cp)) then
               if (ev(10) > tiny(0.0_cp))tk=Tref+log(v0T/v00)/ev(10)
            else
               a=0.5*ev(11)
               b=ev(10)
               c=-b*tref -a*tref*tref + ev(12)/tref -log(v0t/v00)
               d=-ev(12)

               a1=b/a
               a2=c/a
               a3=d/a
               q=(3.0*a2 - a1*a1)/9.0
               r=(9.0*a1*a2-27.0*a3-2.0*a1*a1*a1)/54.0

               dd=q**3+r**2
               if (dd >= 0.0) then
                  s1=(r + sqrt(q**3+r**2))**(1.0/3.0)
                  s2=(r - sqrt(q**3+r**2))**(1.0/3.0)

                  t1=s1+s2-a1/3.0
                  if (abs(s1-s2) <= tiny(0.0)) t2=-0.5*(s1+s2)-a1/3.0

                  if (t1 > 0.0) then
                     if (t1 >=TKmin .and. t1 <=TKmax)then
                        tk=t1
                     else
                        err_CFML%IErr=1
                        err_CFML%Msg='No valid solution for temperature'
                     end if

                  elseif (t2 > 0.0) then
                     if (t2 >=TKmin .and. t2 <=TKmax)then
                        tk=t2
                     else
                        err_CFML%IErr=1
                        err_CFML%Msg='No valid solution for temperature'
                     end if
                  else
                     err_CFML%IErr=1
                     err_CFML%Msg='No valid solution for temperature'
                  end if

               else
                  th=acos(r/sqrt(-q**3))
                  t1=2.0*sqrt(-q)*cos(th/3.0)-a1/3.0
                  t2=2.0*sqrt(-q)*cos((th+2.0*pi)/3.0) -a1/3.0
                  t3=2.0*sqrt(-q)*cos((th+4.0*pi)/3.0) -a1/3.0
                  if (t1 > 0.0) then
                     if (t1 >=TKmin .and. t1 <=TKmax) tk=t1
                  elseif (t2 > 0.0) then
                     if (t2 >=TKmin .and. t2 <=TKmax) tk=t2
                  elseif (t3 > 0.0) then
                     if (t3 >=TKmin .and. t3 <=TKmax) tk=t3
                  else
                     err_CFML%IErr=1
                     err_CFML%Msg='No valid solution for temperature'
                  end if
               end if
            end if

         case (3) ! HP
            ! modified RJA 31.03.2014 to fix coding errors
            a=ev(10)
            b=-2.0*(10.0*ev(10) + ev(11))
            c=1.0 - ev(10)*tref + 2.0*(10.0*ev(10) + ev(11))*sqrt(tref) - (V0t/V00)

            d=b*b - 4.0*a*c
            if (d >= 0.0) then
               x1=(-b+sqrt(d))/(2.0*a)
               x2=(-b-sqrt(d))/(2.0*a)
               if (x1 > 0.0) then
                  t1=x1*x1
                  if (t1 >=TKmin .and. t1 <=TKmax) tk=t1
               else if (x2 > 0.0) then
                  t2=x2*x2
                  if (t2 >=TKmin .and. t2 <=TKmax) tk=t2
               end if
            else
               err_CFML%IErr=1
               err_CFML%Msg='No valid solution for temperature'
            end if

         case (4) ! Kroll
            !>>>>> kp=ev(3) : version before 11/11/2016
            if (eospar%icross == 2)then
               kp=ev(5)                    !uses Anderson delta
            else
               kp=ev(3)
            end if

            th_e=ev(11)
            a=th_e/Tref
            b=-1.0/(kp*(kp+2.0))
            eps0=(a**2)*exp(a)/(exp(a)-1.0)**2

            x1=(-b*(1.0+kp))*eps0/(ev(10)*th_e)
            x2=1.0 - (((v0t/v00)+kp)/(1.0+kp))**(1.0/b)
            x3=1.0/(exp(a)-1)

            c=1.0/(x1*x2+x3)
            if (c <= -1.0_cp) then
               err_CFML%IErr=1
               err_CFML%Msg='Temperature calculated as negative'
            else
               d=log(1.0 +c)
               t1=th_e/d
               if (t1 >=TKmin .and. t1 <=TKmax) then
                  tk=t1
               else
                  err_CFML%IErr=1
                  err_CFML%Msg='No valid solution for temperature'
               end if
            end if

         case (5) ! Salje
            x1=0.0
            x2=0.0
            x3=0.0

            y=ev(10)*ev(11)
            p=v00**(1.0/3.0) -y

            a=y**3
            b=3.0 * y*y*p
            c=3.0 * y * p*p
            d=p**3 - v0t

            a1=b/a
            a2=c/a
            a3=d/a
            q=(3.0*a2 - a1*a1)/9.0
            r=(9.0*a1*a2-27.0*a3-2.0*a1*a1*a1)/54.0

            dd=q**3+r**2
            if (dd >= 0.0) then
               s1=(r + sqrt(q**3+r**2))**(1.0/3.0)
               s2=(r - sqrt(q**3+r**2))**(1.0/3.0)

               x1=s1+s2-a1/3.0
               if (abs(s1-s2) <= tiny(0.0)) then
                  x2=-0.5*(s1+s2)-a1/3.0
               end if
            else
               th=acos(r/sqrt(-q**3))
               x1=2.0*sqrt(-q)*cos(th/3.0)-a1/3.0
               x2=2.0*sqrt(-q)*cos((th+2.0*pi)/3.0) -a1/3.0
               x3=2.0*sqrt(-q)*cos((th+4.0*pi)/3.0) -a1/3.0
            end if
            if (abs(x1) >= 1.0) t1=2.0*ev(11)/log((1.0+x1)/(x1-1.0))
            if (abs(x2) >= 1.0) t2=2.0*ev(11)/log((1.0+x2)/(x2-1.0))
            if (abs(x3) >= 1.0) t3=2.0*ev(11)/log((1.0+x3)/(x3-1.0))

            if (t1 > 0.0) then
               if (t1 >=TKmin .and. t1 <=TKmax) tk=t1
            end if
            if (t2 > 0.0 .and. tk < 0.0) then
               if (t2 >=TKmin .and. t2 <=TKmax) tk=t2
            end if
            if (t3 > 0.0 .and. tk < 0.0) then
               if (t3 >=TKmin .and. t3 <=TKmax) tk=t3
            end if

         case(6)
            th_e=ev(11)
            a=th_e/Tref
            eps0=(a**2)*exp(a)/(exp(a)-1.0)**2

            x1=eps0/(ev(10)*ev(2)*ev(11))
            x2=1.0/(exp(a)-1.0)

            c=1.0/(v0t*x1+x2)
            d=log(1.0 + c)

            t1=th_e/d
            if (t1 >=TKmin .and. t1 <=TKmax) tk=t1

         case(7)  !For MGD Pthermal, no real meaningful value of T to return, so set as Tref
            tk=eospar%tref
      end select
   End Function Get_Temperature_P0

End SubModule EoS_004