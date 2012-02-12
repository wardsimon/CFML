  Module cost_magfunctions 
      use CFML_GlobalDeps,                only: cp,sp,eps
      use CFML_crystallographic_symmetry, only: space_group_type

      use CFML_IO_Formats,                only: file_list_type
      use CFML_string_utilities,          only: l_case, setnum_std
      use CFML_Atom_TypeDef,              only: mAtom_List_Type
      use CFML_crystal_Metrics,           only: Crystal_Cell_Type
      use CFML_Keywords_Code_Parser,      only: NP_Max,NP_Refi,v_Vec,v_Shift,v_Bounds,v_BCon,v_Name,v_List, &
                                                VState_to_AtomsPar
      use CFML_Magnetic_Structure_Factors
      use CFML_Magnetic_Symmetry

      use prep_input

      implicit none
      private

      public::  General_cost_function, Readn_Set_CostFunctPars, Write_CostFunctPars, Write_FinalCost, &
                Write_SOL_mCFL

      logical,                      public :: err_cost=.false.
      character(len=132),           public :: err_mess_cost=" "
      type (space_group_type),      public :: SpG

      Integer, parameter,             public :: N_costf=1
      Integer,dimension(0:N_costf),   public :: Icost
      real,   dimension(0:N_costf),   public :: Wcost
      real,   dimension(0:N_costf),   public :: P_cost !Partial cost

      real,            public :: wavel
      character(len=3),public :: diff_mode="NUC"   ! XRA for x-rays, ELE for electrons

    Contains

!******************************************!
    Subroutine Readn_Set_CostFunctPars(file_cfl)
!******************************************!
       !---- Arguments ----!
       Type(file_list_type),   intent( in)  :: file_cfl
       !---- Local variables ----!
       character(len=132)   :: line
       real                 :: w,tol
       integer              :: i,ier,j,nrel

       Icost=0
       wcost=0.0

       i=0

       do j=1,file_cfl%nlines
          line=adjustl(file_cfl%line(j))
          line=l_case(line)
          if (len_trim(line) == 0) cycle
          if (line(1:1) =="!" .or. line(1:1) =="#") cycle
          i=index(line,"!")
          if( i /= 0) line=trim(line(1:i-1))

          select case (line(1:4))

              case ("opti")    !Optimization

              i=index(line,"f2mag")
              if(i == 0) then
                   Icost(0)=0; wcost(0)=0.0
                 else
                   Icost(0)=1
                   read(unit=line(i+9:),fmt=*,iostat=ier) w
                   if(ier /= 0) then
                      wcost(0)=1.0
                   else
                      wcost(0)= w
                   end if
                 end if

          end select

       end do

       return
    End Subroutine Readn_Set_CostFunctPars

!******************************************!
    Subroutine Write_CostFunctPars(lun)
!******************************************!
       !---- Arguments ----!
       integer,   intent( in)    :: lun
       !---- Local variables ----!
       character(len=132)   :: line
       real                 :: w,tol
       integer              :: i,ier


       Write(unit=lun,fmt="(/,a)")    "=================================="
       Write(unit=lun,fmt="(a)")      "= Cost functions to be minimized ="
       Write(unit=lun,fmt="(a,/)")    "=================================="

       do i=0,n_costf

          if(icost(i) == 0) cycle
          select case (i)

              case (0)    !Optimization "f2mag"
                 Write(unit=lun,fmt="(a,f8.4)") &
                 "  => Cost(F2mag): Optimization of Sum(|Obs-Scale*Calc|^2) / Sum (sigma^2)) / (Nobs-Nref), with weight: ",wcost(i)
          end select

       end do

       return
    End Subroutine Write_CostFunctPars

!******************************************!
    Subroutine Write_FinalCost(lun)
!******************************************!
       !---- Arguments ----!
       integer,   intent( in)    :: lun
       !---- Local variables ----!
       character(len=132)   :: line
       real                 :: w,tol
       integer              :: i


       Write(unit=lun,fmt="(/,a)")    "====================================="
       Write(unit=lun,fmt="(a)")      "= Final Cost of minimized functions ="
       Write(unit=lun,fmt="(a,/)")    "====================================="

       do i=0,n_costf

          if(icost(i) == 0) cycle

          select case (i)

              case (0)    !Optimization "F2mag"

                 Write(unit=lun,fmt="(a,f8.4,a,f12.4)") &
                 "  => Cost(F2mag): Optimization of Sum(|Obs-Scale*Calc|^2) / Sum(sigma^2)) / (Nobs-Nref), with weight: ",wcost(i),&
                 "  Final Cost: ",P_cost(0)

          end select

       end do

       return
    End Subroutine Write_FinalCost

!******************************************!
    Subroutine General_Cost_function(v,cost)
!******************************************!
      real,dimension(:),    intent( in):: v
      real,                 intent(out):: cost
      
       !---- Arguments ----!

      !---- Local variables ----!
      integer :: i,ic, nop, nlist=1, numv
      integer, dimension(1) :: List

      v_shift(1:NP_Refi)=v(1:NP_Refi)-v_vec(1:NP_Refi)  !Calculate the shifts w.r.t. old configuration
      v_vec(1:NP_Refi)=v(1:NP_Refi)

      numv=count(abs(v_shift) > 0.00001, dim=1)
      if(numv == 1) then
         i=maxloc(abs(v_shift),dim=1)
         List(1)=v_list(i)
         call VState_to_AtomsPar(mA,mode="S",MGp=MGp,Mag_dom=Mag_dom) !Update Atomic parameters with the proper constraints, 
         cost=0.0

         do ic=0,N_costf

            if(Icost(ic) == 0) cycle

            Select Case(ic)

               case(0)      !F2mag
                     call Calc_sqMiV_Data
                     call Cost_sqMiV(P_cost(0),Scale)
                     cost=cost+ P_cost(0)* WCost(0)

            End Select
         end do

      else            !New configuration

         call VState_to_AtomsPar(mA,mode="V",MGp=MGp,Mag_dom=Mag_dom)   !Update Atomic parameters with the proper constraints
         cost=0.0

         do ic=0,N_costf

            if(Icost(ic) == 0) cycle

            Select Case(ic)

               case(0)      !F2mag=sqMiV
                     call Calc_sqMiV_Data
                     call Cost_sqMiV(P_cost(0),Scale)
                     cost=cost+ P_cost(0)* WCost(0)

            End Select
         end do
      end if

      return
    End Subroutine General_Cost_function

!******************************************!
    Subroutine Cost_sqMiV(cost,Scale)
!******************************************!
       real,                 intent(out):: cost
       real,                 intent(in out):: Scale
       !---- Local variables ----!
       integer              :: j,n

       n=Oblist%Nobs
       Scale=sum( [(MhList%Mh(j)%sqMiV*Oblist%Ob(j)%Gobs*Oblist%Ob(j)%wGobs,j=1,n)])/ &
             sum( [(MhList%Mh(j)%sqMiV**2 * Oblist%Ob(j)%wGobs,j=1,n)] )
       cost=sum(([(Oblist%Ob(j)%wGobs* (Oblist%Ob(j)%Gobs-Scale*MhList%Mh(j)%sqMiV)**2,j=1,n)]))/(n-NP_Refi)
       return
    End Subroutine Cost_sqMiV

!******************************************!
Subroutine Write_SOL_mCFL(lun,file_cfl,mA,Mag_dom,comment)
!******************************************!
    !!----
    !!---- Subroutine Write_SOL_mCFL(lun,file_cfl,Cel,SpG,mA,Mag_dom,comment)
    !!----    integer,                  intent(in)      :: lun
    !!----    type(file_list_type),     intent (in out) :: file_cfl
    !!----    type (atom_list_type),    intent(in)      :: mA
    !!----    type (Magnetic_Domain_type),optional,intent(in)    :: Mag_dom
    !!----    character(len=*),optional,intent(in)      :: comment
    !!----
    !!----    Write a file mCFL
    !!----
    !!---- Created: January - 2012
    !!

       !---- Arguments ----!
       integer,                  intent(in)           :: lun
       type(file_list_type),     intent (in out)      :: file_cfl
       type (mAtom_list_type),   intent(in)           :: mA
       type (Magnetic_Domain_type),optional,intent(in):: Mag_dom
       character(len=*),optional,intent(in)           :: comment

       !----- Local variables -----!
       integer                         :: j,i,n,ier
       integer                         :: num_matom,num_skp,num_dom,ik,im,ip
       real,dimension(3)               :: Rsk,Isk
       real(kind=cp),dimension(12)     :: coef
       real(kind=cp)                   :: Ph
       character(len=132)              :: lowline,line
       character(len=30)               :: forma
       character(len=14)               :: pop
       logical                         :: skp_begin, bfcoef_begin, magdom_begin=.true.

       if(present(comment)) write(unit=lun,fmt="(a)") "TITLE "//trim(comment)
       write(unit=lun,fmt="(a)") "!  Automatically generated CFL file (Write_SOL_mCFL)"
       write(unit=lun,fmt="(a)") "!  "

      num_matom=0
      num_dom=0
      i=1

      do 
      i=i+1
      if(i >= file_cfl%nlines) exit

       lowline=l_case(adjustl(file_cfl%line(i)))
 
       if(lowline(1:6) == "magdom".and.magdom_begin) then
        num_dom=num_dom+1 
        ip=index(lowline,":")
        forma="(a  ,a14)"
        write(unit=forma(3:4),fmt="(i2)") ip
        write(pop,"(2f7.4)") Mag_Dom%Pop(1:2,num_dom)/sum(Mag_dom%Pop) !normalization to Sum_dom=1
        write(unit=file_cfl%line(i),fmt=forma) lowline(1:ip),pop

        do
         i=i+1
         lowline=adjustl(l_case(file_cfl%line(i)))         
         if(lowline(1:6) == "magdom") then
          num_dom=num_dom+1 
          ip=index(lowline,":")
          forma="(a  ,a14)"
          write(unit=forma(3:4),fmt="(i2)") ip
          write(pop,"(2f7.4)") Mag_Dom%Pop(1:2,num_dom)/sum(Mag_dom%Pop) !normalization to Sum_dom=1
          write(unit=file_cfl%line(i),fmt=forma) lowline(1:ip),pop
         else
          i=i-1
          magdom_begin=.false.
          exit
         endif
        enddo
       endif! end magdom

       if(lowline(1:5) == "matom") then
          num_matom=num_matom+1 !max NMatom=mA%natoms
          num_skp=0
          skp_begin=.true.
          bfcoef_begin=.true.
          write(unit=lun,fmt="(a)") trim(file_cfl%line(i))
          cycle
       endif

       if(lowline(1:3) == "skp".and.skp_begin) then
        num_skp=num_skp+1 !max mA%atom(num_matom)%nvk 
        read(unit=lowline(4:),fmt=*,iostat=ier) ik,im,Rsk,Isk,Ph
        if(MGp%Sk_type == "Spherical_Frame") then
         Rsk(:)=mA%atom(num_matom)%Spher_Skr(:,ik)
         Isk(:)=mA%atom(num_matom)%Spher_Ski(:,ik)
         Ph=mA%atom(num_matom)%mphas(ik)
        else
         Rsk(:)=mA%atom(num_matom)%Skr(:,ik)
         Isk(:)=mA%atom(num_matom)%Ski(:,ik)
         Ph=mA%atom(num_matom)%mphas(ik)        
        endif
        
        write(unit=file_cfl%line(i),fmt='(a,i8,i3,7f8.3)') 'skp',ik,im,Rsk,Isk,Ph
        do
         i=i+1
         lowline=adjustl(l_case(file_cfl%line(i)))         
         if(lowline(1:3) == "skp") then
          num_skp=num_skp+1
          read(unit=lowline(4:),fmt=*,iostat=ier) ik,im,Rsk,Isk,Ph
          if(MGp%Sk_type == "Spherical_Frame") then
           Rsk(:)=mA%atom(num_matom)%Spher_Skr(:,ik)
           Isk(:)=mA%atom(num_matom)%Spher_Ski(:,ik)
           Ph=mA%atom(num_matom)%mphas(ik)
          else
           Rsk(:)=mA%atom(num_matom)%Skr(:,ik)
           Isk(:)=mA%atom(num_matom)%Ski(:,ik)
           Ph=mA%atom(num_matom)%mphas(ik)        
          endif
          write(unit=file_cfl%line(i),fmt='(a,i8,i3,7f8.3)') 'skp',ik,im,Rsk,Isk,Ph
         else
          i=i-1
          skp_begin=.false.
          exit
         endif
        enddo

       endif! end Rsk,Isk,Ph

       forma="(a6,i8,i3,  f8.3)"

       if(lowline(1:6) == "bfcoef".and.bfcoef_begin) then
        num_skp=num_skp+1 !max mA%atom(num_matom)%nvk 
        read(unit=lowline(7:),fmt=*,iostat=ier) ik,im
        n=abs(MGp%nbas(im))
        write(unit=forma(11:12),fmt="(i2)") n+1
        read(unit=lowline(7:),fmt=*,iostat=ier) ik,im,coef(1:n),ph
        coef(1:n)= mA%atom(num_matom)%cbas(1:n,ik) 
        Ph=mA%atom(num_matom)%mphas(ik)        
        write(unit=file_cfl%line(i),fmt=forma) 'bfcoef',ik,im,coef(1:n),Ph        
        do
         i=i+1
         lowline=adjustl(l_case(file_cfl%line(i)))         
         if(lowline(1:6) == "bfcoef") then
          num_skp=num_skp+1
          read(unit=lowline(7:),fmt=*,iostat=ier) ik,im
          n=abs(MGp%nbas(im))
          write(unit=forma(11:12),fmt="(i2)") n+1
          read(unit=lowline(7:),fmt=*,iostat=ier) ik,im,coef(1:n),ph
          coef(1:n)= mA%atom(num_matom)%cbas(1:n,ik) 
          Ph=mA%atom(num_matom)%mphas(ik)        
          write(unit=file_cfl%line(i),fmt=forma) 'bfcoef',ik,im,coef(1:n),Ph
         else
          i=i-1
          bfcoef_begin=.false.
          exit
        endif
        enddo
       endif! end bfcoef

       write(unit=lun,fmt="(a)") trim(file_cfl%line(i))

      end do

    End Subroutine Write_SOL_mCFL

  End Module cost_magfunctions