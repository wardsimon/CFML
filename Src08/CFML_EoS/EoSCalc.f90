!!----
!!----
!!----
SubModule (CFML_EoS) EoS_030
   Contains
   !!----
   !!---- EOS_CAL
   !!----    Returns elastic properties (not the parameter values) at this P,T for EoS
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine EoS_Cal(P,T,EoSpar,Parvals)
      !---- Arguments ----!
      real(kind=cp),                intent(in)  :: P       ! Pressure
      real(kind=cp),                intent(in)  :: T       ! Temperature
      type(Eos_Type),               intent(in)  :: EoSPar  ! Eos Parameter
      real(kind=cp),  dimension(:), intent(out) :: Parvals ! Output parameter values

      !> Init
      parvals=0.0_cp

      call physical_check(eospar,Pin=p,Tin=t)           ! produce warnings based on P,T
      if (err_CFML%IErr==1)return

      parvals(1)=get_volume(p,t,eospar)
      parvals(2)=k_cal(parvals(1),t,eospar,P=p)
      parvals(3)=kp_cal(parvals(1),t,eospar,P=p)
      parvals(4)=kpp_cal(parvals(1),t,eospar)
      parvals(5)=dKdT_cal(p,t,eospar)           ! dK/dT at this P,T
      parvals(6)=alpha_cal(p,t,eospar)          ! 1/V.dV/dT at this T

   End Subroutine EoS_Cal

   !!----
   !!---- EOS_CAL_ESD
   !!----    Returns esd's of the elastic properties (not the parameter values) at this P and T for EoS
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine EoS_Cal_Esd(P,T,EoSpar,Esd)
      !---- Arguments ----!
      real(kind=cp),                intent(in)  :: P       ! Pressure
      real(kind=cp),                intent(in)  :: T       ! Temperature
      type(Eos_Type),               intent(in)  :: EoSPar  ! Eos Parameter
      real(kind=cp),  dimension(:), intent(out) :: Esd     ! Output esd values

      !---- Local Variables ----!
      real(kind=cp),dimension(n_eospar) :: esdfull

      !> Init
      esd=0.0_cp
      esdfull=0.0_cp

      !> calculate parameter esd's transformed to this P,T
      call transform_Esd(P,T,EoSpar,esdfull)
      esd(1:5)=esdfull(1:5)
      esd(6)=esdfull(10)

      !> make adjustment for only using alpha0
      select case(EoSpar%itherm)
          case(5)               ! Salje
            esd(6)=esd(6)/3.0_cp            ! alpha is 1/3 of p1
      end select
   End Subroutine EoS_Cal_Esd

   !!--++
   !!--++ TRANSFORM_ESD
   !!--++   New generalised version for V,K,Kp,Kpp,dK/dT and alpha
   !!--++
   !!--..  NOTE:
   !!--..       the value of 0.001 was chosen to balance non-linearity of
   !!--..       equations (which implies small value) against the arithmatic
   !!--..       precision (wants big value). this value gives same esd's as
   !!--..       double-precision eosfit5
   !!--++
   !!--++ 02/08/2013
   !!
   Module Subroutine Transform_Esd(P,T,Eospar,Esd)
      !---- Arguments ----!
      real(kind=cp),   intent(in)              :: p        ! Pressure at which to calculate parameter esds
      real(kind=cp),   intent(in)              :: t        ! Temperature at which to calculate parameter esds
      type (EoS_Type), intent(in)              :: eospar   ! The EoS parameters and vcv
      real(kind=cp),dimension(:), intent(out)  :: esd      ! The esd's of Eos parameters at this P and T

      !---- Local Variables ----!
      real(kind=cp),dimension(N_EOSPar,N_EOSPar)  :: d                  ! cross derivatives of parameters
      real(kind=cp),dimension(-2:2,6)             :: par
      real(kind=cp), parameter                    :: factor=0.01_cp     ! shift factor for calculating derivatives
      real(kind=cp)                               :: shift              ! shift applied to parameter
      type (eos_type)                             :: eost               ! local copy of input eos
      integer                                     :: i,j,k

      !> initialisation
      d=0.0_cp
      esd=0.0_cp

      !> loop to calculate d(param_i at P,T)/d(param_j at Pref,Tref)
      do j=1,10
         if (eospar%iuse(j) == 1) then
            do k=-2,2,1
               eost=eospar                                                     !reset eos params
               if (abs(eospar%params(j)) > tiny(0.0_cp)) then
                  shift=FACTOR*eospar%params(j)
                  if (j == 10) shift=10.0*shift                    ! alpha0
               else
                  shift=1.0_cp/eospar%factor(j)                    ! shift to a parameter with zero value
                  if (j == 5) shift=0.002                          ! dK/dT has print factor 1, but typical value -.02
               end if

               eost%params(j)=eospar%params(j)+float(k)*shift    ! apply shift to a parameter
               call EoS_Cal(P,T,Eost,Par(k,1:6))                 ! calc resulting parvals
            end do
            d(1:6,j)=(par(-2,1:6)+8.0_cp*(par(1,1:6)-par(-1,1:6))-par(2,1:6))/(12.0_cp*shift) ! derivative to second order approximation
         end if
      end do

      !> d(1:6,j) contains the derivatives in order V,K0,Kp,Kpp,dK/dT,alpha0 wrt param(j) at Pref,Tref
      !  now switch these to 10 for alpha0 and ignore any other alpha coeffs
      d(10,1:10)=d(6,1:10)
      d(6,1:10)=0.0_cp

      !> Now calculate esd-squared array = vcv(i,i). The input vcv is already linear for linear eos!
      do k=1,N_EOSPar
         do i=1,N_EOSPar
            if (eospar%iref(i) == 0) cycle ! because vcv(i,j) will be 0
            do j=1,N_EOSPar
               if (eospar%iref(j) == 0) cycle ! because vcv(i,j) will be 0
               esd(k)=esd(k)+eospar%vcv(i,j)*d(k,i)*d(k,j)
            end do
         end do
      end do

      !> Final: extra trap June 2016 in case round off leaves esd(i)^2 < 0
      do i=1,n_eospar
          if (esd(i) > tiny(0.0_cp) ) then
             esd(i)=sqrt(esd(i))
          else
             esd(i)=0.0_cp
          end if
      end do
   End Subroutine Transform_Esd
End SubModule EoS_030