!!----
!!----
!!----
SubModule (CFML_EoS) EoS_Partial_P
   implicit none
   Contains

   !!----
   !!---- DERIV_PARTIAL_P
   !!----    Calculate Partial derivates of Pressure respect to Params
   !!----
   !!---- 21/03/2017
   !!
   Module Subroutine Deriv_Partial_P(V,T,Eospar,td,xtype,calc)
      !---- Arguments ----!
      real(kind=cp),                      intent(in)  :: V       ! Volume
      real(kind=cp),                      intent(in)  :: T       ! Temperature
      type(Eos_Type),                     intent(in)  :: Eospar  ! Eos Parameter
      real(kind=cp), dimension(n_eospar), intent(out) :: td      ! derivatives dP/d(param)
      integer,optional,                   intent(in)  :: xtype   ! =0 for V,T input, =1 for Kt,T =2 for Ks,T
      character(len=*),optional,          intent(in)  :: calc    ! 'all' if all derivs required. If absent, then just params with iref=1 calculated

      !---- Local Variables ----!
      integer                                         :: itype
      real(kind=cp), dimension(n_eospar)              :: tda,tdn                ! analytic and numeric derivatives
      character(len=10)                               :: cstring

      !> Init
      td=0.0_cp
      itype=0
      if(present(xtype))itype=xtype
      cstring='ref'
      if(present(calc))cstring=u_case(adjustl(calc))

      !> Calculate derivatives by both methods if possible: correct values are returned for linear
      call Deriv_Partial_P_Numeric(V,T,Eospar,tdn,itype,cstring)

      !> Default to numeric, because they are always available:
      td(1:n_eospar)=tdn(1:n_eospar)

      if (itype == 0 .and. .not. Eospar%pthermaleos)then
         call Deriv_Partial_P_Analytic(V,T,Eospar,tda)
         if (eospar%itran ==0 .and. eospar%imodel /=6) then ! imodel=6 is APL, not yet coded
            td(1:4)=tda(1:4)                                ! analytic for Vo and moduli terms because these are exact even at small P
         end if
      end if
   End Subroutine Deriv_Partial_P

   !!--++
   !!--++ DERIV_PARTIAL_P_ANALYTIC
   !!--++
   !!--++ PRIVATE
   !!--++ Calculates the partial derivatives of P with respect to the EoS
   !!--++ at a given v,t point, and returns  them in array td
   !!--++
   !!--.. 27-Feb-2013: RJA edits to use IREF as in original
   !!--++
   !!--++ 17/07/2015
   !!
   Module Subroutine Deriv_Partial_P_Analytic(V,T,Eospar,td)
      !---- Arguments ----!
      real(kind=cp),                      intent(in)  :: V       ! Volume
      real(kind=cp),                      intent(in)  :: T       ! Temperature
      type(Eos_Type),                     intent(in)  :: Eospar  ! Eos Parameter
      real(kind=cp), dimension(n_eospar), intent(out) :: td      ! derivatives dP/d(param)

      !---- Local Variables ----!
      real(kind=cp), dimension(n_eospar) :: ev
      real(kind=cp)                      :: vv0,k0,delt, delt2,vt0,dpdvt,deltinv
      real(kind=cp)                      :: f,b,c,cv,kp,kpp
      real(kind=cp)                      :: term1,term2,term3,vterm
      real(kind=cp)                      :: a,f52,da,db,dc,nu   ! new guys
      real(kind=cp),dimension(3)         :: abc                 !Tait parammeters

      !> Local copy
      call EoS_to_Vec(eospar,ev)    ! Volume or linear case is covered: ev contains volume-like parameters and
                                    ! all derivatives will be calculated as volume-like derivatives
      cv=v                          ! cv has the current volume or 'a' input. It is needed for dP/dVo

      select case (eospar%Itherm)
         case (0)
            vt0=eospar%params(1)                ! vt0 is the V (or a) at this T and P=0

         case (1:)
            vt0=get_V0_T(t,eospar)
      end select
      vv0=v/vt0                                 ! vv0 is V(P,T)/V(P,T=0), or the linear equivalent

      k0=Get_K0_T(T,eospar)                     ! returns M(T) for linear
      if (eospar%linear) k0=k0/3.0_cp
      kp=Get_Kp0_T(T,eospar)
      if (eospar%linear) kp=kp/3.0_cp
      kpp=Get_Kpp0_T(T,eospar)
      if (eospar%linear) kpp=kpp/4.0_cp

      !> now get the finite strain for vv0
      f=strain(vv0,eospar)                      ! use strain to get finite strain from v/v0 or a/a0

      !> Using volume dimensions for remainder of calculations
      vv0=1.0_cp/vv0                               ! vv0 now V0/V for ease of use in rest of subprogram
      if (eospar%linear) then
         vt0=vt0**3.0_cp           ! Vo at this T
         vv0=vv0**3.0_cp           ! V/Vo at this T
         cv=cv**3.0_cp             ! current volume at this P,T
      end if

      td=0.0_cp
      select case (eospar%imodel)
         case (1) ! Murnaghan: validated algebra, 10/04/2013
            td(1)=vv0**(kp-1.0_cp)*k0/cv
            td(2)=(vv0**kp-1.0_cp)/kp
            td(3)=k0/kp*(vv0**kp*(log(vv0)-1.0_cp/kp) +1.0_cp/kp)
            td(4)=0.0_cp

         case (2) ! Birch-Murnaghan: new code 10-11/04/2013 Validated
            !> Assign coefficients of f in expansion of P: use K0 as that is K at P=0
            !  f was calculated above at init
            a=K0                           ! P=3(1+2f)^(5/2) . (af + bf^2 +cf^3)
            b=0.0_cp
            c=0.0_cp
            f52=(1.0_cp+2.0*f)**2.5_cp
            if (eospar%iorder > 2) b=1.5_cp*K0*(kp-4.0_cp)
            if (eospar%iorder == 4)c = 1.5_cp*K0*(K0*kpp + (kp-4.0_cp)*(kp-3.0_cp)+35.0_cp/9.0_cp)

            !> dP/dVo
            td(1)=f52/vt0 * ( a + (7.0_cp*a+2.0_cp*b)*f + (9.0_cp*b+3.0_cp*c)*f*f + 11.0_cp*c*f*f*f)

            !> dP/dKo:  da, db, dc are da/dparam etc for each param K0, Kp0, Kpp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder > 2) db = 1.5_cp*(kp-4.0_cp)
            if (eospar%iorder == 4)dc = 3.0_cp*K0*kpp + 1.5_cp*((kp-4.0_cp)*(kp-3.0_cp)+35.0_cp/9.0_cp)

            td(2)= 3.0*f52*f* (1.0 + db*f + dc*f*f)

            !> dP/dKp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder > 2) db = 1.5_cp*K0
            if (eospar%iorder == 4)dc = 1.5_cp*K0*(2.0_cp*kp-7.0_cp)

            td(3) = 3.0*f52* (db + dc*f) *f*f

            !> dP/dKpp0
            db = 0.0_cp
            dc = 0.0_cp
            if (eospar%iorder == 4) then
               dc = 1.5_cp*K0*K0
               td(4) = 3.0*f52*dc*f*f*f
            else
               td(4)=0.0_cp
            end if

         case (3) ! Vinet: new code and validated 11/04/2013
            if (eospar%iorder == 2) then
               td(1) = k0/vt0 *(1.0_cp+f)/(1.0_cp-f)**2.0_cp
               td(2) =      3.0_cp * f /(1.0_cp-f)**2.0_cp
               td(3) = 0.0_cp
               td(4) = 0.0_cp

            else
               nu=1.5_cp*(kp-1.0)
               td(1) = k0/vt0 * (1.0_cp + (1.0_cp+nu)*f - nu*f*f)/(1.0_cp-f)**2.0_cp * exp(nu*f)
               td(2) =      3.0_cp * f /(1.0_cp-f)**2.0_cp * exp(nu*f)
               td(3) = k0 * 4.5_cp * f*f/(1.0_cp-f)**2.0_cp * exp(nu*f)
               td(4) = 0.0_cp
            end if

         case (4) ! Natural Strain: new code and validated 11/04/2013
            !> Assign coefficients of f in expansion of P: use K0 as that is K at P=0,T=Tdata
            !  f was calculated above at init
            !
            !  the coefficients a and b are for P= 3 K0 (V0/V) f(1+ af + bf^2)
            a=0.0_cp
            b=0.0_cp
            if (eospar%iorder > 2) a=1.5_cp*(kp-2.0_cp)
            if (eospar%iorder ==4) b=1.5_cp*(1.0_cp + K0*kpp + (kp-2.0_cp) + (kp-2.0_cp)**2.0_cp)

            !> dP/dVo
            td(1) = K0/cv * (1.0_cp  + (2.0_cp*a+3.0_cp)*f + 3.0_cp*(a+b)*f*f +3.0_cp*b*f*f*f)

            !> d(p)/d(k0)
            td(2) = 3.0_cp * vv0 * (f + a*f*f + (b + 1.5_cp*kpp*k0)*f*f*f)

            !> dp/dKp0
            da = 0.0_cp
            db = 0.0_cp
            if(eospar%iorder > 2)  da=1.5_cp
            if (eospar%iorder == 4)db = 3.0_cp*kp - 4.5_cp

            td(3) = 3.0_cp * k0 * vv0 * (da + db*f)*f*f

            !> dp/dKpp0
            td(4)=0.0_cp
            if (eospar%iorder == 4) td(4) = 4.5_cp* k0*k0 * vv0 *f*f*f

         case (5) ! Tait: new code and validated against Murnaghan and against finite diffs 15/06/2013
            abc=get_tait(eospar,t)
            vv0=1.0/vv0        ! now vv0 is v/v0

            !> dP/dVo
            td(1)= k0 * vv0**2.0_cp/cv * ( (vv0 + abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3) -1.0_cp)

            !> d(p)/d(k0)
            td(2)=abc(1)*abc(3)*(( (vv0 + abc(1) -1.0_cp)/abc(1))**(-1.0_cp/abc(3)) -1.0_cp)

            !> dp/dKp0
            if (eospar%iorder > 2) then
               da= k0*kpp/(1.0_cp +kp+k0*kpp)**2.0_cp      ! da/dKp
               db= 1.0_cp/k0 + kpp/(1.0_cp+kp)**2.0_cp       ! db/dKp
               dc= 1.0_cp/(kp*kp + kp -k0*kpp) - ((1.0_cp)+kp+ &
                   k0*kpp)*(2.0_cp*kp +1.0_cp)/(kp*kp + kp -k0*kpp)**2.0_cp

               vterm = (vv0/abc(1) + 1.0_cp - 1.0_cp/abc(1))**(-1.0_cp/abc(3))
               term1=-1.0_cp/abc(2)/abc(2) * db * vterm
               term2= 1.0_cp/abc(2)*vterm/abc(3)/abc(3)*dc*log((vv0-1.0_cp)/abc(1) + 1.0_cp)
               term2=term2+((vv0-1.0_cp)/abc(1) + 1.0_cp)**(-1.0_cp/abc(3)-1.0_cp)/abc(2)/abc(3)/abc(1)/abc(1)*(vv0-1.0_cp)*da
               term3=db/abc(2)/abc(2)
               td(3)=term1+term2+term3
            end if

            !> dp/dKpp0
            if (abs(kpp) > tiny(0.) .and. eospar%iorder == 4)then          ! if Kpp0 = 0 then the derivative is undefined
               da=-1.0_cp*k0*abc(1)/(1.0_cp+kp+k0*kpp)
               db=-1.0_cp/(1.0+kp)
               dc=(k0*(kp+1.0_cp)**2.0_cp)/(kp*kp+kp-k0*kpp)**2.0_cp

               term1= -1.0_cp/abc(2)/abc(2)*db*(vterm-1.0_cp)
               term2= 1.0_cp/abc(2)*vterm/abc(3)/abc(3)*dc*log((vv0-1.0_cp)/abc(1) + 1.0_cp)
               term3= ((vv0-1.0_cp)/abc(1) + 1.0_cp)**(-1.0_cp/abc(3)-1.0_cp)/abc(2)/abc(3)/abc(1)/abc(1)*(vv0-1.0_cp)*da
               td(4)=term1+term2+term3
            end if
      end select

      if (eospar%ITherm > 0) then
         !> First do the parameters common to all thermal
         delt=t-eospar%tref
         delt2=t*t-eospar%tref*eospar%tref
         deltinv=0.0_cp
         if (eospar%tref > 1.0_cp) deltinv=1.0_cp/t-1.0_cp/eospar%tref

         !> Adjust dp/dv(0,t) to dp/dv(0,0):
         dpdvt=td(1)
         td(1)=td(1)*vt0/ev(1)

         !> dp/dk(Tref,0) = dp/dk(T,0)

         !> d(k)/d(t) ...(the temperature dependence of k)
         td(5)=td(2)*delt

         !> Now do the specific derivatives for each equation, where possible. Program uses numeric derivs anyway
         select case(eospar%itherm)     !  dp/dalpha=dp/dV(0T) . dV(0T)/dalpha
            case(1)                ! Berman
               !> dp/dalph0:
               td(10)=dpdvt*delt*ev(1)

               !> dp/dalpha1
               td(11)=dpdvt*0.5_cp*delt*delt*ev(1)

            case(2)                ! Fei
               !> dp/dalph0:
               td(10)=dpdvt*delt*vt0

               !> dp/dalpha1
               td(11)=dpdvt*0.5_cp*delt2*vt0

               !> dp/dalpha2
               td(12)= -1.0_cp*dpdvt*deltinv*vt0
         end select
      end if

      !> fix derivatives for linear eos
      if (eospar%linear) then
         td(1) =td(1) *(3.0_cp*eospar%params(1)**2.0_cp)
         td(2:5)=td(2:5)/3.0_cp

         select case(eospar%itherm)                 ! thermal expansion terms
            case(1:2)                               ! Berman, Fei, alpha terms
               td(10:12)=td(10:12)*3.0_cp

            case(4)                                 ! Holland-Powell alpha
               td(10)=td(10)*3.0_cp
         end select                                 !other thermal equations have nothing to convert
      end if
   End Subroutine Deriv_Partial_P_Analytic

   !!--++
   !!--++ DERIV_PARTIAL_P_NUMERIC
   !!--++    Calculates the partial derivatives of P with respect to the EoS
   !!--++    at a given property and t point, and returns  them in array td
   !!--++
   !!--++
   !!--++ 17/03/2017
   !!
   Module Subroutine Deriv_Partial_P_Numeric(X1,X2,Eospar,td,xtype,calc)
      !---- Arguments ----!
      real(kind=cp),                      intent(in) :: X1,X2   ! The two parameter values (eg P and T)
      type(Eos_Type),                     intent(in) :: Eospar  ! Eos Parameters
      real(kind=cp), dimension(n_eospar), intent(out):: td      ! derivatives dP/d(param)
      integer,         optional,          intent(in) :: xtype   ! =0 for V,T input, =1 for Kt,T =2 for Ks,T
      character(len=*),optional,          intent(in) :: calc    ! 'all' if all derivs required.
                                                                ! If absent, then just params with iref=1 calculated


      !---- Local Variables ----!
      type(Eos_Type)                 :: Eost                 ! Eos Parameter local copy
      real(kind=cp), dimension(-2:2) :: p                    ! array for calc p values
      real(kind=cp)                  :: delfactor,del,d_prev,delmin
      integer                        :: i,j,icycle,itype
      logical                        :: warn            ! local warn flag
      logical                        :: iall            ! .true. if all derivs required


      !> Initialise
      call clear_error()
      itype=0
      if(present(xtype))itype=xtype

      iall=.false.      ! default is calc params with iref=1
      if (present(calc) )then
         if (index(u_case(adjustl(calc)),'ALL') > 0) iall=.true.
      end if

      td=0.0_cp
      eost=eospar
      call Set_Eos_Use(eost)

      warn=.false.

      !> Check
      if (itype < 0 .or. itype > n_data_types) then
         err_CFML%IErr=1
         Err_CFML%Msg='No type set for deriv_partial_p'
         return
      end if

      !> Set the inital shift factor (fraction of parameter value)
      delfactor=0.01_cp

      do i=1,n_eospar
         if (eospar%iref(i) == 1 .or. iall) then    ! only refined parameters
            del=delfactor*eospar%params(i)         ! the initial shift estimate
            delmin=delfactor/eospar%factor(i)      ! scale required min shift by print factors
            if (abs(del) < delmin) del=delmin      ! trap param values that are zero
            icycle=0                               ! iteration count
            d_prev=0.0_cp

       iter:do                                ! top of loop over iterations
               do j=-2,2,1
                  eost=eospar                                             !reset eos params
                  eost%params(i)=eost%params(i)+float(j)*del              ! apply shift to a parameter
                  p(j)=get_pressure_x(x1,x2,eost,itype)                  ! calc resulting P
               end do

               td(i)=(p(-2)+8.0_cp*(p(1)-p(-1))-p(2))/(12.0_cp*del)       ! derivative to second order approximation

               !write(6,'(''  Param # '',i2,'' Cycle '',i2,'': deriv = '',ES14.6,'' delp = '',f5.3,'' del= '',f9.4)')i,icycle,td(i),p(2)-p(-2),del

               ! to trap problems
               if (err_CFML%IErr==1) then
                  td(i)=d_prev        ! previous cycle value
                  call clear_error()  ! clear errors
                  warn=.true.       ! warning flag
                  exit iter
               end if

               if (abs(td(i)) < 1.0E-8) exit iter                         ! zero deriv
               if (icycle > 0 .and. &
                  abs(d_prev-td(i))/td(i) < 1.0E-4) exit iter    ! deriv converged to 1 part in 10^4

               d_prev=td(i)                ! store last deriv value
               del=2.0_cp*del              ! increase the shift
               icycle=icycle+1
               if (icycle > 5) exit iter    ! Do not allow 2*shift to exceed 64% of param value
            end do iter
         end if
      end do

      if (warn) then
         err_CFML%IErr=2
         err_CFML%Msg='Error calculating some derivatives in Least squares'
      end if

      !> no need to fix derivatives for linear eos by this method
   End Subroutine Deriv_Partial_P_Numeric

End SubModule EoS_Partial_P