!!----
!!----
!!----
SubModule (CFML_gSpaceGroups) Spg_060
   Contains
   !!----
   !!---- SET_SPACEGROUP
   !!----    For non magnetic Space group 
   !!----
   !!---- 07/05/19
   !!
   Module Subroutine Set_SpaceGroup(Str, SpaceG, NGen, Gen)
      !---- Arguments ----!
      character(len=*),                          intent(in ) :: Str  
      class(spg_type),                           intent(out) :: SpaceG
      integer,                         optional, intent(in ) :: NGen
      character(len=*),  dimension(:), optional, intent(in ) :: Gen
      
      !---- Local Variables ----!
      integer            :: i,j,n_gen, n_it, d, ier
      integer            :: n_laue, n_pg, nfin
      character(len=40), dimension(:), allocatable :: l_gen
      character(len=20)                            :: str_HM, str_HM_std, str_Hall, str_CHM
      character(len=5)                             :: car
      character(len=256)                           :: gList
      real(kind=cp), dimension(3)                  :: co
      
      integer                          :: n
      type(rational), dimension(3,3)   :: P,Mp,Mc,M
      type(rational), dimension(3,3,6) :: A
      
      logical :: by_Gen=.false.
             
      !> Init
      n_gen=0
      gList=" "
      n_it=0
      str_HM=" "
      str_Hall=" "
      str_CHM=" "
      str_HM_std=" "
      n_laue=0
      n_pg=0
      
      call Init_SpaceGroup(SpaceG)
      if (present(ngen) .and. present(gen)) n_gen=ngen
      
      !> Check
      if (len_trim(Str) <= 0 .and. n_gen <= 0) then 
         err_CFML%Ierr=1
         Err_CFML%Msg="Set_SpaceGroup@GSPACEGROUPS: Argument not valid for define Space Group!"
         return
      end if  
      
      !> Spacegroup from Str argument
      if (n_gen ==0) then
         
         !> Init Spgr_Info
         call Set_Spgr_Info()
      
         !> Is a IT  number?
         read(unit=str, fmt=*, iostat=ier) n_it
         if (ier == 0) then
            gList=get_IT_Generators(str)  ! IT take our default choice for SpaceGroup
            call Get_SpaceGroup_Symbols(Str, str_HM, str_Hall)
            do i=1,NUM_SPGR_INFO
               if (trim(str_HM) /= trim(spgr_info(i)%hm)) cycle
               n_laue=spgr_info(i)%laue
               n_pg=  spgr_info(i)%pg
               exit
            end do
            str_HM_std=trim(str_HM)
         end if   
            
         !> Is HM symbol defined in SPGR_INFO?
         if (n_it ==0) then   
            str_HM=u_case(trim(Str))
            do i=1,NUM_SPGR_INFO
               if (trim(str_HM) /= trim(spgr_info(i)%hm)) cycle
               n_it=spgr_info(i)%n
               str_hall=spgr_info(i)%hall
               n_laue=spgr_info(i)%laue
               n_pg=  spgr_info(i)%pg
               
               !> Get generators list from standard
               write(unit=car, fmt='(i3)') n_it
               call Get_SpaceGroup_Symbols(car, str_HM_std)
               if (trim(str_HM_std)==trim(str_HM)) then
                  gList=get_IT_Generators(car)
               end if    
               
               exit
            end do 
         end if   
         
         !> Is Hall symbol defined in SPGR_INFO?
         if (n_it ==0) then
            Str_Hall=l_case(trim(Str))
            if (str_hall(1:1) == '-') then
               str_hall(2:2)=u_case(str_hall(2:2))
            else
               str_hall(1:1)=u_case(str_hall(1:1))
            end if 
            do i=1,NUM_SPGR_INFO
               if (trim(str_Hall) /= trim(spgr_info(i)%hall)) cycle
               n_it=spgr_info(i)%n
               str_hm=spgr_info(i)%hm
               n_laue=spgr_info(i)%laue
               n_pg=  spgr_info(i)%pg
               
               !> Get generators list from standard
               write(unit=car,fmt=*) n_it
               call Get_SpaceGroup_Symbols(car, str_HM_std)
               if (trim(str_HM_std)==trim(str_HM)) then
                  gList=get_IT_Generators(car)
               end if
               exit
            end do  
         end if  
         
         !> Is compact HM symbol?
         if (n_it ==0) then
            Str_CHM=u_case(trim(Str))
            call Get_SpaceGroup_Symbols(str_CHM, str_HM)
            do i=1,NUM_SPGR_INFO
               if (trim(str_HM) /= trim(spgr_info(i)%hm)) cycle
               n_it=spgr_info(i)%n
               str_hall=spgr_info(i)%hall
               n_laue=spgr_info(i)%laue
               n_pg=  spgr_info(i)%pg
               
               !> Get generators list from standard
               write(unit=car,fmt=*) n_it
               call Get_SpaceGroup_Symbols(car, str_HM_std)
               if (trim(str_HM_std)==trim(str_HM)) then
                  gList=get_IT_Generators(car)
               end if
               exit
            end do 
         end if 

         
         if (len_trim(gList) == 0) then
            call Get_Generators_from_Hall(str_hall,n_gen, l_gen)
            if (Err_CFML%Ierr /= 0) return
         else
            call Get_Generators_from_Str(gList, d, n_gen, l_gen)   
         end if   
         
      else
         by_Gen=.true.
         allocate(l_gen(n_gen))
         l_gen=gen(1:n_gen)
      end if
      
      !> Generate the spacegroup
      call Group_Constructor(l_gen,SpaceG)
      if (Err_CFML%Ierr /= 0) return
      
      !> More info
      if (.not. by_Gen) then
         SpaceG%Numspg=n_it
         SpaceG%spg_symb=str_HM(1:1)//l_case(str_HM(2:))
      end if

      !> Identify Group
      call Identify_Group(SpaceG)
      
   End Subroutine Set_SpaceGroup
   
End SubModule Spg_060   
