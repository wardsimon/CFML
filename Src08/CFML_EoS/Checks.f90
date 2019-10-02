!!----
!!----
!!----
SubModule (CFML_EoS) EoS_028
   Contains
   
   !!----
   !!---- CHECK_SCALES
   !!----
   !!----
   !!---- 31/05/2019 
   !!
   Module Subroutine Check_Scales(E,dat)
      !---- Arguments ----!
      type(Eos_Type),                       intent(in)  :: E          ! EoS 
      type (eos_data_list_type), optional,  intent(in)  :: dat        ! data structure
      
      !---- Local Variables ----!
      character(len=40)       :: name
      
      !> Init
      call clear_error()
      
      !>Checks of EoS only
      
      !> APL
      if (e%imodel == 6)then
         if (len_trim(E%pscale_name) == 0) then
            err_CFML%IErr=2
            if (e%linear) then
               err_CFML%Msg='APL EoS must have a Pscale (and M0) in GPa'
            else
               err_CFML%Msg='APL EoS must have a Pscale (and K0) in GPa'
            end if
         end if
        
         if (len_trim(E%vscale_name) == 0 .or. index(U_case(E%Vscale_name),'A') == 0)then 
            err_CFML%IErr=2
            if (len_trim(err_CFML%Msg) == 0)then
               if (e%linear)then
                  err_CFML%Msg='APL EoS must have a Vscale and L0 in A'
               else
                  err_CFML%Msg='APL EoS must have a Vscale and V0 in A^3'
               end if
            else
               if (e%linear)then
                  err_CFML%Msg=trim(err_CFML%Msg)//' and a Vscale and L0 in A'
               else
                  err_CFML%Msg=trim(err_CFML%Msg)//' and a Vscale and V0 in A^3'
               end if
            end if
         end if
      end if
      
      !> If MGD type thermal EoS, must have eos%pscale_name and eos%_Vscale_name
      if (e%itherm == 7) then
         if (len_trim(E%pscale_name) == 0)then
            err_CFML%IErr=2
            err_CFML%Msg='MGD EoS must have a Pscale in kbar or GPa'
         end if
         if (len_trim(E%Vscale_name) == 0 .or. .not. VscaleMGD(E))then
            err_CFML%IErr=2
            if (len_trim(err_CFML%Msg) == 0)then
               err_CFML%Msg='MGD EoS must have a Vscale in cm3/mol'
            else
               err_CFML%Msg=trim(err_CFML%Msg)//' and a Vscale in cm3/mol'
            end if
         end if
         if (len_trim(err_CFML%Msg) /= 0) err_CFML%Msg=trim(err_CFML%Msg)//' set to get correct results. '       
      end if
      if (.not. present(dat))return
      
      !>For all EoS compare data and eos scales
      if (len_trim(E%pscale_name) /= 0 .and. len_trim(dat%Pscale_name) /=0)then
         if (trim(u_case(adjustl(E%pscale_name))) /= trim(u_case(adjustl(dat%Pscale_name))))then
            err_CFML%IErr=2
            if (len_trim(err_CFML%Msg) > 0) err_CFML%Msg=trim(err_CFML%Msg)//' And'
            err_CFML%Msg=trim(err_CFML%Msg)//' Pscales of data and EoS are different.'
         end if
      end if
   
      if (len_trim(E%vscale_name) /= 0 )then
         if (e%linear)then
            name=trim(u_case(adjustl(dat%Lscale_name)))
         else
            name=trim(u_case(adjustl(dat%Vscale_name)))
         end if
         if (len_trim(name) /= 0)then
            if (trim(u_case(adjustl(E%vscale_name))) /= trim(name))then
               err_CFML%IErr=2
               if (len_trim(err_CFML%Msg) > 0) err_CFML%Msg=trim(err_CFML%Msg)//' And'
               err_CFML%Msg=trim(err_CFML%Msg)//' Vscales of data and EoS are different'
            end if
         end if
           
      end if    
   End Subroutine Check_Scales
   
   !!--++
   !!--++ PHYSICAL_CHECK
   !!--++    Check if the parameters have physical sense
   !!--++
   !!--++ 22/05/2017
   !!
   Module Subroutine Physical_Check_old(P,T,E)
      !---- Arguments ----!
      real(kind=cp),    intent(in) :: p  ! Pressure
      real(kind=cp),    intent(in) :: t  ! Temperature
      type(Eos_Type),   intent(in) :: E  ! EoS object

      !---- Local variables ----!
      character(len=20)   :: car
      real(kind=cp)       :: tlimit,v,pinf

      !> Basic checks
      v=get_volume(p,t,e)
      if (err_CFML%IErr==1)then         ! added 22/05/2017
         write(unit=car, fmt='(2f10.1)') p, t
         car=adjustl(car)
         err_CFML%Msg='Volume cannot be calculated at P,T = '//trim(car)
         return             ! have to return if error, because other tests cannot be done
      end if

      if (v < tiny(0.0) ) then
         err_CFML%IErr=1
         write(unit=car, fmt='(2f10.1)') p, t
         car=adjustl(car)
         err_CFML%Msg='Volume calculated as zero or negative at P,T = '//trim(car)
      end if

      if (e%imodel > 0) then
         if (.not. e%linear .and. k_cal(v,t,e) < tiny(0.0_cp)) then
            err_CFML%IErr=1
            write(unit=car, fmt='(2f10.1)') p, t
            car=adjustl(car)
            err_CFML%Msg='Bulk modulus calculated as zero or negative at P,T = '//trim(car)
         end if
      end if

      !> Check validity of thermal model
      select case(e%itherm)
         case (2)                ! Fei:
            if (e%params(12) > tiny(0.0_cp)) then  ! non-physical V and divergent alpha at low T when alpha2 .ne. 0
               tlimit=(2.0_cp*e%params(12)/e%params(11))**(1.0_cp/3.0_cp)
               if (t < tlimit) then
                  err_CFML%IErr=1
                  write(unit=car,fmt='(f5.1)')tlimit
                  car=adjustl(car)
                  err_CFML%Msg='Fei equation yields non-physical behaviour below T = '//trim(car)//'K'
               end if
            else if(e%params(12) < tiny(0.0_cp)) then  ! alpha2 < 0
               tlimit=sqrt(-1.0_cp*e%params(12)/e%params(10))
               if (t < tlimit) then
                  err_CFML%IErr=1
                  write(unit=car,fmt='(f5.1)')tlimit
                  car=adjustl(car)
                  err_CFML%Msg='Fei equation yields non-physical behaviour below T = '//trim(car)//'K'
               end if
            end if

         case(3)               ! HP 1998: trap non-physical behaviour at low T
            tlimit=((10.0_cp*e%params(10)+e%params(11))/e%params(10))**2.0_cp
            if (t < tlimit) then
               err_CFML%IErr=1
               write(unit=car,fmt='(f5.1)')tlimit
               car=adjustl(car)
               err_CFML%Msg='HP1998 equation yields non-physical behaviour below T = '//trim(car)//'K'
            end if

         case(1,4,5,6,7)
            if (t < 0.0_cp) then
               err_CFML%IErr=1
               err_CFML%Msg='T  less than zero K'
            end if
      end select

      !> Produce warning for curved phase boundaries: Pinflection = a/-2b when Ttr=Tr0+aP+bP^2
      if (e%itran>0 .and. abs(e%params(23)) > tiny(0.0) )then
         pinf=abs(e%params(22)/2.0/e%params(23))
         if (abs(p/pinf -1.0) < 0.1) then
            err_CFML%IErr=1
            err_CFML%Msg='P in region of boundary inflection P: PVT calculations may be inaccurate or wrong'
         end if
      end if
   End Subroutine Physical_Check_old

   !!--++
   !!--++ PHYSICAL_CHECK
   !!--++    Check if the parameters have physical sense
   !!--++
   !!--++ 31/05/2019 
   Module Subroutine Physical_Check(Ein,Pin,Tin,Vin)
      !---- Arguments ----!
      type(Eos_Type),        intent(in) :: Ein  ! EoS object
      real(kind=cp),optional,intent(in) :: pin  ! Pressure
      real(kind=cp),optional,intent(in) :: tin  ! Temperature
      real(kind=cp),optional,intent(in) :: vin  ! volume

      !---- Local variables ----!
      integer             :: n
      character(len=100)   :: car
      real(kind=cp)       :: tlimit,pinf,p,v,t,pmin,vmin
      type(eos_type)      :: e,eiso
      logical             :: vpresent


      !>local copies
      E=Ein
      T=e%tref
      
      !> check PVT present
      n=0
      if (present(Tin))then
         T=Tin
         n=n+1
      end if
      P=0.0_cp
      if (present(Pin))then
         P=Pin
         n=n+1
      end if

      !> Volume : This is needed for most tests of most EoS
      V=0.0_cp
      Vpresent=.false.
      if (present(Vin))then
         if (Vin < 0.0_cp) then
            err_CFML%IErr=1
            err_CFML%Msg='Volume is negative'
            return
         end if
         V=Vin
         n=n+1
         Vpresent=.true.
      end if
      if (n == 0)return      !no arguments
      if (e%imodel > 0 .and. e%itherm > 0 .and. n < 2)return   ! not enough arguments for PT eos

      !> Positive T
      if (t < 0.0_cp) then
         err_CFML%IErr=1
         err_CFML%Msg='T is less than zero K'
         return
      end if

      !> Now check for valid parameters at reference
      call EoSParams_Check(E)
      if (err_CFML%IErr==1)return

      !> Now check pthermal and isothermal seperately: Pthermal is first
      if (e%pthermaleos)then
         if (e%params(3) > 0._cp)then  ! K limit does not occur if Kp or Mp negative
            !> FIRST find the V at which K=K0/2 at Tref, WITHOUT using pressure
            eiso=e
            eiso%pthermaleos=.false.
            eiso%itherm=0
            vmin=get_volume_K(eiso%params(2)/2.0_cp,eiso%tref,eiso)
            if (vpresent) then
               if (v > vmin) then
                  err_CFML%IErr=1
                  err_CFML%Msg='Thermal pressure EoS not valid at this V and T: the V is too big so the compressional part of the EoS at Tref is not valid'
                  return
               end if
               if (k_cal(v,t,e) < tiny(0._cp)) then
                  err_CFML%IErr=1
                  err_CFML%Msg='Thermal pressure EoS not valid at this V and T: the K is negative (maybe because of q large?)'
                  return
               end if    
            else
               !> No volume input. So calculate the isochor Pressure of Vmin at the input T, and compare to input P
               if (get_k(p,t,e) < tiny(0._cp))then
                  err_CFML%IErr=1
                  err_CFML%Msg='Thermal pressure EoS not valid at this P and T: the K is negative (maybe because of q large?)'
                  return
               end if  
               if (get_volume(p,t,e) > Vmin) then
                  err_CFML%IErr=1
                  err_CFML%Msg='Thermal pressure EoS not valid at this P and T: the V is too big so the compressional part of the EoS at Tref is not valid'
                  return
               end if 
            end if
         end if  
         
      else  !isothermal or no thermal: check thermal part first for T being valid
         !> Check validity of normal-type thermal model: only needs T
         select case(e%itherm)
            case (2)                ! Fei:
               if (e%params(12) > tiny(0.0_cp)) then  ! non-physical V and divergent alpha at low T when alpha2 .ne. 0
                  tlimit=(2.0_cp*e%params(12)/e%params(11))**(1.0_cp/3.0_cp)
                  if (t < tlimit) then
                     err_CFML%IErr=1
                     write(unit=car,fmt='(f5.1)')tlimit
                     car=adjustl(car)
                     err_CFML%Msg='Fei equation yields non-physical behaviour below T = '//trim(car)//'K'
                     return
                  end if
               else if (e%params(12) < tiny(0.0_cp)) then  ! alpha2 < 0
                  tlimit=sqrt(-1.0_cp*e%params(12)/e%params(10))
                  if (t < tlimit) then
                     err_CFML%IErr=1
                     write(unit=car,fmt='(f5.1)')tlimit
                     car=adjustl(car)
                     err_CFML%Msg='Fei equation yields non-physical behaviour below T = '//trim(car)//'K'
                     return
                  end if
               end if

            case(3)               ! HP 1998: trap non-physical behaviour at low T
               tlimit=((10.0_cp*e%params(10)+e%params(11))/e%params(10))**2.0_cp
               if (t < tlimit) then
                  err_CFML%IErr=1
                  write(unit=car,fmt='(f5.1)')tlimit
                  car=adjustl(car)
                  err_CFML%Msg='HP1998 equation yields non-physical behaviour below T = '//trim(car)//'K'
                  return
               end if
         end select

         !> Now check the validity of Eos params at T
         call pveos_check(P,V,T,e,vpresent)
         if (err_CFML%IErr==1)then
            err_CFML%Msg='Compressional EoS not valid at this PV: '//trim(err_CFML%Msg)
            return
         end if
      end if

      !> If got to here, now check that properties at P,T,V valid of Full EoS
      !> because  checks  above are for the PV part and the TV part, without transitions.
      !> all must be valid for the Eos to be valid
      if (.not. vpresent .and. e%itherm /=7)then        !only done if V not provided at start
         v=get_volume(p,t,e)
         if (err_CFML%IErr==1)then         ! added 22/05/2017
            write(unit=car, fmt='(2f10.1)') p, t
            car=adjustl(car)
            err_CFML%Msg='Volume cannot be calculated at P,T = '//trim(car)
            return
         end if

         if (v < tiny(0.0) ) then
            err_CFML%IErr=1
            write(unit=car, fmt='(2f10.1)') p, t
            car=adjustl(car)
            err_CFML%Msg='Volume calculated as zero or negative at P,T = '//trim(car)
            return
         end if
      end if

      if (.not. e%linear .and.  V > tiny(0._cp))then
         if (K_cal(V,T,E,P) < tiny(0._cp))then
            err_CFML%IErr=1
            write(unit=car, fmt='(2f10.1)') p, t
            car=adjustl(car)
            err_CFML%Msg='Bulk modulus calculated as zero or negative at P,T = '//trim(car)
            return
         end if
      end if

      !> Produce warning for curved phase boundaries: Pinflection = a/-2b when Ttr=Tr0+aP+bP^2
      if (e%itran>0 .and. abs(e%params(23)) > tiny(0.0) )then
         pinf=abs(e%params(22)/2.0/e%params(23))
         if (abs(p/pinf -1.0) < 0.1) then
            err_CFML%IErr=1
            err_CFML%Msg='P in region of boundary inflection P: PVT calculations may be inaccurate or wrong'
         end if
      end if
   End Subroutine Physical_Check

   !!----
   !!---- PVEOS_CHECK
   !!----    Checks compressional part of Eos at P,T for validity (normally that K > 0, or K > K0/2) for volume
   !!----    Does not do transition part
   !!----
   !!---- 17/12/2018
   !!
   Module Subroutine PVEos_Check(Pin,Vin,Tin,ein,vpresent)
      !---- Arguments ----!
      real(kind=cp),optional,intent(in) :: pin  ! Pressure
      real(kind=cp),optional,intent(in) :: vin  ! volume
      real(kind=cp),optional,intent(in) :: tin  ! Temperature
      type(Eos_Type),   intent(in)      :: Ein  ! EoS object
      logical,          intent(in)      :: vpresent ! .true. when the Vin is meaningful

      !---- Local variables ----!
      real(kind=cp)       :: p,v,t
      real(kind=cp),dimension(3)       :: abc
      real(kind=cp)       :: bp,kc,step,plim,kprev,Vnew,Vprev,klim
      type(eos_type)      :: e

      !> Init
      if(ein%linear)return

      !> local copies
      e=ein
      p=pin
      t=tin
      v=vin

      !> set no transitions
      e%itran=0
      
      !> This routine is private and only called from physical_check
      !> therefore if pthermaleos then T will always be Tref. But set it to be safe, and suppress all thermal part
      if (e%pthermaleos)then
         t=e%tref
         e%pthermaleos=.false.
         e%itherm=0
      end if

      !>prelim stuff for each type of eos, used if vpresent or not

      ! now do further tests dependening on Vpresent
      ! When V is present, calculate K from V,T
      ! And error state when K < K = K(P=0,T)/2, except for Murnaghan which is stable to K=0
      if (vpresent)then
         if (v > Get_V0_T(T,E))then
            select case(e%imodel)
               case(1) ! Murngahan: limit is when K=0
                  if (p < -1.0_cp*e%params(2)/e%params(3)) err_CFML%IErr=1

              case(2,3,4,5,6)   ! BM, Vinet, NS, Tait, APL
                 if (K_cal(V,T,E) < get_K0_T(T,E)/2.0) err_CFML%IErr=1
            end select
         end if

      else      ! V was not given, but p was
         if (p < 0.0_cp)then
            select case(e%imodel)
               case(1) ! Murngahan
                  if (p + 1.0_cp*e%params(2)/e%params(3) < tiny(0.0)) err_CFML%IErr=1

               case(2:6) ! find V that gives K = K(P=0,T)/2, by iteration
                  Klim=get_K0_T(T,E)/2.0_cp
                  V=get_volume_K(Klim,e%tref,e)
                  plim=get_pressure(V,T,e) 
                  if (p < plim) then
                     err_CFML%IErr=1
                     return
                  end if
            end select
         end if
      end if
   End Subroutine pveos_check
   
   !!----
   !!---- EOSPARAMS_CHECK
   !!----    Check for Params that are invalid for all cases.
   !!----
   !!--.. NOTE:
   !!--..      Written 7-April-2014 to fix invalid eos parameters and report as error message
   !!--..      This is different from physical_check, because that checks properties at specific
   !!--..      p and T
   !!----
   !!---- 11/07/2016
   !!
   Module Subroutine EoSParams_Check(EoSPar)
      !---- Argument ----!
      type (EoS_Type), intent(in out) :: EoSPar

      !---- Local Variables ----!
      real(kind=cp)     :: pinf
      character(len=80) :: text

      !> Init
      call clear_error()

      !> Check that v0 is positive
      if (eospar%params(1) < tiny(0.0)) then
         eospar%params(1)=1.0_cp
         err_CFML%IErr=1
         if (eospar%linear) then
            err_CFML%Msg='*****WARNING: a0 was < 0. Not allowed! Reset to 1.00'
         else
            err_CFML%Msg='*****WARNING: V0 was < 0. Not allowed! Reset to 1.00'
         end if
      end if

      !> Check K0 is positive for V-P (-ve K ok for linear)
      if (.not. eospar%linear .and. eospar%imodel > 0) then
         if (eospar%params(2) < tiny(0.0_cp) ) then
            eospar%params(2)=10.0_cp
            err_CFML%IErr=1
            if (len_trim(err_CFML%Msg) == 0) then
               err_CFML%Msg=' *****WARNING: K0 was < 0. Not allowed! Reset to 10.0'
            else
               err_CFML%Msg=trim(err_CFML%Msg)//' And K0 was < 0. Not allowed! Reset to 10.0'
            end if
         end if
      end if

      !> Check Tref is positive
      if (eospar%tref < -1.0_cp*tiny(0.0_cp)) then
         eospar%tref=0.0_cp
         err_CFML%IErr=1
         if (len_trim(err_CFML%Msg) == 0) then
            err_CFML%Msg=' *****WARNING: Tref was < 0. Not allowed! Reset to 0 K'
         else
            err_CFML%Msg=trim(err_CFML%Msg)//' And Tref was < 0. Not allowed! Reset to 0 K'
         end if
      end if

      !> Thermal cases
      select case(eospar%itherm)  ! for specific thermal parameters
         case (4,6,7)    !>Kroll or Pthermal must have characteristic T > o
            if (eospar%params(11) < 0.1) then
               eospar%params(11)=eospar%Tref
               if (eospar%Tref < 0.1) eospar%params=0.1
               err_CFML%IErr=1
               if (len_trim(err_CFML%Msg) == 0) then
                  err_CFML%Msg=' *****WARNING: '//trim(eospar%comment(11))//' was =< 0. Not allowed! Reset to Tref'
               else
                  err_CFML%Msg=trim(err_CFML%Msg)//' And '//trim(eospar%comment(11))//' was =< 0. Not allowed! Reset to Tref'
               end if
            end if
      end select

      !> Produce warning for curved phase boundaries: Pinflection = a/-2b when Ttr=Tr0+aP+bP^2
      if (eospar%itran>1 .and. abs(eospar%params(24)) > tiny(0.0) ) then
         pinf=-1.0*eospar%params(22)/2.0/eospar%params(24)
         if (abs(pinf) < 10.0) then
            err_CFML%IErr=1
            write(text,'(a,f4.1,1x,a)')'Phase boundary inflects at P ~ ',pinf,trim(eospar%pscale_name)
            if (len_trim(err_CFML%Msg) == 0) then
               err_CFML%Msg=' *****WARNING: '//trim(text)
            else
               err_CFML%Msg=trim(err_CFML%Msg)//' And '//trim(text)
            end if
         end if
      end if

      !>If MGD and linear warn that this is not generally valid
      if (eospar%itherm == 7 .and. eospar%linear) then
         err_CFML%IErr=1
         text='Linear MGD EoS only has valid parameters if the material is cubic'
         if (len_trim(err_CFML%Msg) == 0) then
            err_CFML%Msg=' *****WARNING: '//trim(text)
         else
            err_CFML%Msg=trim(err_CFML%Msg)//' And '//trim(text)
         end if
      end if
   End Subroutine EoSParams_Check
   
End SubModule EoS_028   