!!----
!!----
!!----
SubModule (CFML_EoS) EoS_014
   Contains
   
   !!----
   !!---- DKDT_CAL
   !!----    Calculate the derivative dK/dt (or dM/dt in the linear case) at P and T
   !!----
   !!---- 17/07/2015
   !!
   Module Function dKdT_Cal(P,T,Eospar,DeltaT) Result(dKdT)
      !---- Arguments ----!
      real(kind=cp),            intent(in) :: P       ! Pressure
      real(kind=cp),            intent(in) :: T       ! Temperature
      type(Eos_Type),           intent(in) :: EoSPar  ! Eos Parameter
      real(kind=cp),  optional, intent(in) :: DeltaT  ! Delta T
      real(kind=cp)                        :: dKdT
           
      !---- Local Variables ----!
      integer       :: j
      real(kind=cp) :: del,vlimit,tlimit,tcal,Ttr
      real(kind=cp),dimension(-2:2):: kpt
      
      !> Init: set step size for numerical differentiation
      del=30.0_cp                        ! good number for accuracy
      if (present(deltat)) del=deltat
      if (t-2.0*del < 0.0_cp) del=t/2.1          ! prevents going to negative T
      
      !> Code to prevent crossing a phase boundary
      if (eospar%itran > 0) then      
         Ttr = Get_Transition_Temperature(P,EosPar)
         if (transition_phase(P,T+2.0*del,eospar) .neqv. transition_phase(P,T,eospar)) del=0.4*abs(Ttr-T)
         if (transition_phase(P,T-2.0*del,eospar) .neqv. transition_phase(P,T,eospar)) del=0.4*abs(Ttr-T)
      end if
      
      !> Code to stop MGD EoS going into illegal large volume above T
      if (eospar%itherm == 7)then    
         tlimit=t+2.0*del
         do                                        ! search for positive K at this P
            if (get_K(p,tlimit,eospar) > 0._cp .and. err_CFML%Ierr==0)exit
            call clear_error() 
            tlimit=tlimit-0.1_cp*(tlimit-t)
            if (tlimit < t)exit                    ! should never happen because P,T is valid    
         end do
         del=0.4*abs(tlimit-t)
      end if      

      !> Trap close to zero K
      if (t < 1.0_cp)then
         dKdT=Get_K(P,1.0,EosPar)-Get_K(P,0.0,EosPar)
      else
         do j=-2,2,1
            tcal=t+real(j)*del                 ! apply shift to t
            kpt(j)=Get_K(P,tcal,EosPar)        ! calc resulting K
         end do
         dKdT=(kpt(-2)+8.0_cp*(kpt(1)-kpt(-1))-kpt(2))/(12.0_cp*del)     ! Derivative to second order approximation
      end if
      
      !> No linear conversion is required because get_K returns values for "linear Kp" = Mp,
      !> so kppc is already dMp/dP = Mpp

   End Function dKdT_Cal
   
End SubModule EoS_014   