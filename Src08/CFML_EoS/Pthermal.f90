!!----
!!----
!!----
SubModule (CFML_EoS) EoS_011
   Contains
   
   !!--++
   !!--++ PTHERMAL
   !!--++
   !!--++ PRIVATE
   !!--++  Calculate Pthermal from eosparameters at temperature T
   !!--++
   !!--++ 10/09/2013
   !!
   Module Function Pthermal(V,T,EosPar) Result(Pth)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume: not needed for HP2011, needed for MGD, is 'a' if linear
      real(kind=cp),  intent(in) :: T       ! Temperature
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: pth,thtref,exp0,eta0,vlocal
      real(kind=cp) :: gammaV, thetaD,factor
      real(kind=cp),dimension(n_eospar) :: ev

      !> Local copy
      call eos_to_vec(eospar,ev)    !handle linear case

      select case (eospar%itherm)
         case (6) ! Thermal pressure from Holland and Powell 2011
            thtref=ev(11)/eospar%tref         ! T_einstein/Tref
            exp0=exp(thtref)                  ! exp(T_Ein/Tref)
            eta0= thtref*thtref*exp0/(exp0-1.0_cp)**2.0_cp  ! eta at Tref

            pth = ev(10)*ev(2)*ev(11)/eta0*( 1.0_cp/(exp(ev(11)/t)-1.0_cp) -1.0_cp/(exp0 -1.0_cp))

         case (7) ! MGD in the form of Kroll et al (2012)
            vlocal=v
            thetaD=get_DebyeT(V,EosPar)                 ! Both get_Debye and get_grun expect 'a' value if linear
            gammaV=get_grun_V(V,Eospar)
            if (eospar%linear) vlocal=v**3.0_cp          ! cube length to get volume for Pthermal
            pth=gammaV/Vlocal*(EthDebye(T,thetaD,eospar)-EthDebye(eospar%tref,thetaD,eospar))

            !>handle scaling: EthDebye returns energy in Jmol-1. Then if V in m3/mol  Eth/V is in J/m3=Pa
            !>Try to convert volume to m3/mol
            factor=1.0
            if (index(U_case(eospar%pscale_name),'GPA') > 0)  factor=1.0E-9
            if (index(U_case(eospar%pscale_name),'KBAR') > 0) factor=1.0E-8
            
            if (VscaleMGD(eospar)) factor=factor*1.0E+6     !test for cm3/mol or equivalent in eos%vscale_name

            pth=pth*factor

         case default
            pth=0.0_cp
      end select
   End Function Pthermal
   
   !!----
   !!---- VSCALEMGD
   !!----
   !!---- 31/05/2019 
   !!
   Module Function VscaleMGD(E) Result(MGD)
      !---- Arguments ----!
      type(Eos_Type),intent(in)  :: E          ! EoS    
      logical                    :: MGD        ! .true. if e%vscale_name is cm3/mol           
   
      !---- Local Variables ----!
      character(len=len(e%vscale_name)) :: vname

      !> Init
      MGD=.false.  
      
      vname=adjustl(U_case(e%vscale_name))
      if(len_trim(vname) == 0)return

      if (index(vname,'CM') > 0 .and. index(vname,'3') > 0 .and. index(vname,'MOL') > 0) MGD=.true.
   End Function VscaleMGD
   
   !!--++
   !!--++ ETHDEBYE
   !!--++
   !!--++ PRIVATE
   !!--++  Calculates the Debye thermal Energy in Jmol(-1)
   !!--++  because R is given in Jmol(-1)K(-1)
   !!--++
   !!--++ 18/11/2015
   !!
   Module Function EthDebye(T,Theta,Eospar) Result(Eth)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: T       ! Temperature
      real(kind=cp),  intent(in) :: Theta   ! Debye T
      type(Eos_Type), intent(in) :: EoSPar  ! Eos Parameter

      !---- Local Variables ----!
      real(kind=cp) :: Eth
      real(kind=8)  :: x

      if (T < 0.1) then
         Eth=0.0_cp

      else
         x=theta/t
         Eth=3.0_cp*eospar%params(13)*8.314*T*debye(3,x)
      end if
   End Function EthDebye

   !!----
   !!---- GET_DEBYET
   !!----    Calculate the Debye Temperature at V
   !!----
   !!---- 16/03/2017
   !!
   Module Function Get_DebyeT(V, Eos) result(DebyeT)
      !---- Arguments ----!
      real(kind=cp),  intent(in) :: V       ! Volume or length
      type(Eos_Type), intent(in) :: EoS     ! Eos Parameters

      !---- Local Variables ----!
      real(kind=cp) :: DebyeT,gammaV

      !> Default
      DebyeT=eos%tref

      select case(eos%itherm)
         case(7) ! MGD Pthermal
            !> For linear uses the same parameter values, no factor of 3
            if (abs(eos%params(19)) < 0.0001) then
               DebyeT=eos%params(11)                  ! when q=0 DebyeT is fixed at reference value
            else
               gammaV=get_Grun_v(v,eos)               ! Get_grun knows about linear/volume
               DebyeT=eos%params(11)*exp((eos%params(18)-gammaV)/eos%params(19))  ! if q=0 then this gives DebyeT=nan
            end if
      end select
   End Function Get_DebyeT
   
End SubModule EoS_011   