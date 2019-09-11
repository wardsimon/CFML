!!----
!!----
!!----
SubModule (CFML_EoS) EoS_010
   Contains
   
   !!--++
   !!--++ GET_PROPS_PTVTABLE
   !!--++   Returns the requested property at two of P,T,V from an eos loaded as a table of values
   !!--++   The volume values in the table are expected to be scaled to V=1.0 at Pref, Tref
   !!--++   and the true V0 is in eospar%params(1)
   !!--++
   !!--++ 14/10/16  
   !!
   Module Function Get_Props_PTVTable(P,T,V,EosPar,Res) Result(Val)
      !---- Arguments ----!
      real(kind=cp),    intent(in) :: P       ! Pressure
      real(kind=cp),    intent(in) :: T       ! Temperature
      real(kind=cp),    intent(in) :: V       ! Volume
      type(Eos_Type),   intent(in) :: EoSPar  ! Eos Parameter
      character(len=*), intent(in) :: Res     ! Parameter requested for calculation (P,T, or V)
      real(kind=cp)                :: Val

      !---- Local Variables ----!
      integer                   :: i,j
      real(kind=cp)             :: VV0,tt,vm,vp,km,kp,va,vb,dvdp
      type(Eos_Type)            :: eosm,eosp,eosv     ! Eos for local Murn at Tminus and Tplus
      character(len=10)         :: Var

      !>Init
      Val=0.0_cp

      !> Check valid request
      if (eospar%imodel/= -1) then
         err_CFML%IErr=1
         write(err_CFML%Msg,'(''Request to get_props_pvttable with invalid imodel #'',i5)')eospar%imodel
         return
      end if

      VV0=V/eospar%params(1)        ! table values are all V/V0

      !> Determine what is the request: P,T,V, K, KP, AL
      Var=Res
      Var=U_Case(adjustl(Var))

      !> Find lower-corner coords for P,T  : works because P and T always ascend
      if (Var(1:1) /= 'P')then                         ! Only test if P supplied in argument, not if P requested from T and V
         do i=1,eospar%table%np
            if (p < eospar%table%ptv(i,1,1)) exit
         end do

         if (i < 2) then
            err_CFML%IErr=2
            write(err_CFML%Msg, &
                 '(''PTV table request for P = '',f6.2,'' smaller than Pmin = '',f6.2,'' at T ='',f6.1,'': Linear guess made'')') &
                 p,eospar%table%pmin,t
            i=2
         end if

         if (i > eospar%table%np-1) then
            err_CFML%IErr=2
            write(err_CFML%Msg,&
                 '(''PTV table request for P = '',f6.2,'' bigger than Pmax = '',f6.2,'' at T ='',f6.1,'': Linear guess made'')') &
                 p,eospar%table%pmax,t
            i=eospar%table%np-1
         end if
      end if

      !> Check the T limits
      do j=1,eospar%table%nt
         if (t < eospar%table%ptv(1,j,2) ) exit
      end do

      if (j < 2 ) then
         err_CFML%IErr=2
         write(err_CFML%Msg, &
              '(''PTV table request for T = '',f6.1,'' smaller than Tmin = '',f6.1,'' at P ='',f6.2,'': Linear guess made'')') &
              t,eospar%table%tmin,p
         j=2
      end if

      if (j > eospar%table%nt-1) then
         err_CFML%IErr=2
         write(err_CFML%Msg, &
              '(''PTV table request for T = '',f6.1,'' bigger than Tmax = '',f6.1,'' at P ='',f6.1,'': Linear guess made'')') &
              t,eospar%table%tmax,p
         j = eospar%table%nt-1
      end if

      !> Get the Murngahan EoS from Pminus at Tminus and Tplus, but set as if a T=Tref
      !> Must then call Eos functions with Tref as argument
      if (Var(1:1) /= 'P')then       ! cannot do this if P being requested
         Call Murn_Ptvtable(i-1,j-1,Eospar,Eosm)
         Call Murn_Ptvtable(i-1,j,Eospar,Eosp)
      end if

      Select case(Var(1:3))
         case('V  ')
            !> Murnaghan interpolation on P axis, followed by linear in T
            Vm=get_volume(p-eospar%table%ptv(i-1,j-1,1),eosm%tref,eosm)
            Vp=get_volume(p-eospar%table%ptv(i-1,j,1),eosp%tref,eosp)
            val=Vm+(Vp-Vm)*(t-eospar%table%ptv(1,j-1,2))/(eospar%table%ptv(1,j,2)-eospar%table%ptv(1,j-1,2))  ! linear T
            val=val*eospar%params(1)

         case('P  ')
            tt=(t-eospar%table%ptv(1,j-1,2))/(eospar%table%ptv(1,j,2)-eospar%table%ptv(1,j-1,2))
            do i=1,eospar%table%np
               vb=va
               va=tt*(eospar%table%ptv(i,j,3)-eospar%table%ptv(i,j-1,3)) +eospar%table%ptv(i,j-1,3)      ! volume at each P row for target T
               if (va < vv0) exit
            end do

            if (i == 1) then
               err_CFML%IErr=2
               write(err_CFML%Msg, &
                    '(''PTV table request for V = '',f6.1,''at T = '',f6.1,'' is below Pmin: Linear guess made'')') v,t

               !> linear guess beyond table bottom
               dvdp=(tt*(eospar%table%ptv(2,j,3)-eospar%table%ptv(2,j-1,3)) +eospar%table%ptv(2,j-1,3)-va)/(eospar%table%ptv(2,j,1)- &
                    eospar%table%ptv(1,j,1))
               val=eospar%table%ptv(1,j,1)- (va-vv0)/dvdp
               return
            end if

            if (i >= eospar%table%np) then
               err_CFML%IErr=2
               write(err_CFML%Msg, &
                    '(''PTV table request for V = '',f6.1,''at T = '',f6.1,'' is above Pmax: Linear guess made'')') v,t
               !> linear guess beyond table bottom
               dvdp=(va-vb)/(eospar%table%ptv(eospar%table%np,j,1)-eospar%table%ptv(eospar%table%np-1,j,1))
               val=eospar%table%ptv(eospar%table%np,j,1)+ (vv0-va)/dvdp
               return
            end if

            !> Set up Murn in Eosm
            call Init_EoS_Type(Eosv)
            eosv%imodel=1
            call set_eos_names(eosv)      ! sets names for params
            eosv%title='Murnaghan for interpolation'
            eosv%params(1)=vb
            call Murn_Ptvtable(i-1,j-1,Eospar,Eosm)
            call Murn_Ptvtable(i-1,j,Eospar,Eosp)
            eosv%params(2)=eosm%params(2)+(eosp%params(2)-eosm%params(2))*tt
            eosv%params(3)=eosm%params(3)+(eosp%params(3)-eosm%params(3))*tt
            val=get_pressure(vv0,eosv%tref,eosv)+eospar%table%ptv(i-1,j,1)

         case('K  ') !input must be p and t
            km=get_k(p-eospar%table%ptv(i-1,j-1,1),eosm%tref,eosm)      ! K at Tminus
            kp=get_k(p-eospar%table%ptv(i-1,j,1),eosp%tref,eosp)          ! K at Tplus
            tt=(t-eospar%table%ptv(1,j-1,2))/(eospar%table%ptv(1,j,2)-eospar%table%ptv(1,j-1,2))
            val=(kp-km)*tt + km

         case('KP ')
            tt=(t-eospar%table%ptv(1,j-1,2))/(eospar%table%ptv(1,j,2)-eospar%table%ptv(1,j-1,2))
            val=eosm%params(3)+(eosp%params(3)-eosm%params(3))*tt

         case('AL ')
            !> Murnaghan interpolation on P axis, followed by linear in T
            Vm=get_volume(p-eospar%table%ptv(i-1,j-1,1),eosm%tref,eosm)
            Vp=get_volume(p-eospar%table%ptv(i-1,j,1),eosp%tref,eosp)
            val=2.0_cp*(Vp-Vm)/(eospar%table%ptv(i-1,j,2)-eospar%table%ptv(i-1,j-1,2))/(Vm+Vp)

         case default
            err_CFML%IErr=1
            write(err_CFML%Msg,'(''Request for '',a1,'' to get_props_pvttable not valid'')')Var(1:1)
      end select
   End Function Get_Props_PTVTable

   !!--++
   !!--++ MURN_INTERPOLATE_PTVTABLE
   !!--++
   !!--++ PRIVATE
   !!--++ Returns ...
   !!--++
   !!--++ 25/012/2016
   !!
   Module Function Murn_Interpolate_PTVTable(I,J,P,Eos) Result(V)
      !---- Arguments ----!
      integer,        intent(in) :: I      ! pointers to table point just beyond P,T required
      integer,        intent(in) :: J      ! pointers to table point just beyond P,T required
      real(kind=cp),  intent(in) :: P      ! pressure required
      type(Eos_Type), intent(in) :: EoS    ! Eos Parameter
      real(kind=cp)              :: V

      !---- Local Variables ----!
      integer         :: is
      real(kind=cp)   :: k12,K23,KP,V0,K0

      !> Set up pointer
      is=i-2
      if (is < 1) is=1

      !>Estimate K and Kp
      k12=-0.5_cp*(eos%table%ptv(is+1,j,3) + eos%table%ptv(is,j,3))*(eos%table%ptv(is+1,j,1) - &
           eos%table%ptv(is,j,1))/(eos%table%ptv(is+1,j,3)-eos%table%ptv(is,j,3))
      k23=-0.5_cp*(eos%table%ptv(is+2,j,3) + eos%table%ptv(is+1,j,3))*(eos%table%ptv(is+2,j,1) - &
           eos%table%ptv(is+1,j,1))/(eos%table%ptv(is+2,j,3)-eos%table%ptv(is+1,j,3))

      kp=2.0_cp*(k23-k12)/(eos%table%ptv(is+2,j,1)-eos%table%ptv(is,j,1))         ! Kp at point is+1
      if (abs(kp) < 0.0001) kp=4.0_cp
      k0=k12+kp*(eos%table%ptv(is+1,j,1)-eos%table%ptv(is+1,j,1))/2.0_cp          ! K at point is+1
      v0=eos%table%ptv(is+1,j,3)                                                  ! V0 at point is+1

      v=v0*(1.0_cp+kp/k0*(p-eos%table%ptv(is+1,j,1)))**(-1.0_cp/kp)
   End Function Murn_Interpolate_PTVTable
   
   !!--++
   !!--++ MURN_PTVTABLE
   !!--++   Calculate ....
   !!--++
   !!--++ 17/03/2017
   !!
   Module Subroutine Murn_PTVTable(i,j,Eos,Eosm)
      !---- Arguments ----!
      integer,        intent(in)  :: i      ! pointers to table point
      integer,        intent(in)  :: j      ! pointers to table point
      type(Eos_Type), intent(in)  :: EoS    ! Eos Parameter in (with ptvtable)
      type(Eos_Type), intent(out) :: EoSm   ! Eos Parameter of Murn output for point i,j

      !---- Local Variables ----!
      integer         :: im                   ! normally i but if edge point it will be shifted
      real(kind=cp)   :: k2,k4

      !> Set up Murn in Eosm
      call Init_EoS_Type(Eosm)

      eosm%imodel=1
      call set_eos_names(eosm)      ! sets names for params
      eosm%title='Murnaghan for interpolation'

      !> Set up pointer
      im=i
      if (i < 3) im=3
      if (i > eos%table%np-2) im=eos%table%np-2

      !>Vo
      eosm%params(1)=eos%table%ptv(i,j,3)

      !>Estimate K0 at point i,j
      eosm%params(2)=-1.0_cp*eos%table%ptv(im,j,3)*(eos%table%ptv(im+1,j,1)-eos%table%ptv(im-1,j,1))/ &
                     (eos%table%ptv(im+1,j,3)-eos%table%ptv(im-1,j,3))

      !> Kp
      k2=-1.0_cp*eos%table%ptv(im-1,j,3)*(eos%table%ptv(im,j,1)-eos%table%ptv(im-2,j,1))/ &
         (eos%table%ptv(im,j,3)-eos%table%ptv(im-2,j,3))

      k4=-1.0_cp*eos%table%ptv(im+1,j,3)*(eos%table%ptv(im+2,j,1)-eos%table%ptv(im,j,1))/ &
         (eos%table%ptv(im+2,j,3)-eos%table%ptv(im,j,3))

      eosm%params(3)=(k4-k2)/(eos%table%ptv(im+1,j,1)-eos%table%ptv(im-1,j,1))       ! Kp at point im

      if (abs(eosm%params(3)) < 0.0001) eosm%params(3)=4.0_cp

      !> adjust K now if im /= i; Kp should be constant for a Murn
      if (i /= im) eosm%params(2)=eosm%params(2)+eosm%params(3)*(eos%table%ptv(i,j,1)-eos%table%ptv(im,j,1))
   End Subroutine Murn_PTVTable
   
End SubModule EoS_010   