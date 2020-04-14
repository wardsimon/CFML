!!----
!!---- Test 230 Standard Groups
!!----
!!----
!!
Program test230_groups
   !---- Use Modules ----!
   use CFML_GlobalDeps
   use CFML_Symmetry_Tables
   use CFML_Strings
   use CFML_gSpaceGroups

   !---- Local Variables ----!
   character(len=256)                  :: generatorList
   character(len=5)                    :: aux
   character(len=15)                   :: forma,hm
   integer                             :: i,j,L,nsg,lun,nc,inig,fing,ier,minuts
   integer, dimension(:), allocatable  :: cosets
   real(kind=cp)                       :: start, fin,par,secs

   type(Spg_Type)                      :: Grp
   type(Spg_Type), dimension(4096)     :: sGrp

   !> Init
   call Set_Conditions_NumOp_EPS(2048) ! Maximum admissible multiplicity
   write(*,"(a)",advance="no") " => Enter the numbers of the space groups to be analysed (two integers <=230): "
   read(*,*,iostat=ier) inig,fing
   if(ier /= 0) then
     inig=1
     fing=230
   end if

   !> Start Time
   call CPU_TIME(start)

   !> Output results
   open(newunit=lun,file="test_230groups.out",status="replace",action="write")
   do i=inig,fing
      call CPU_TIME(par)
      write(unit=aux,fmt="(i3)") i
      generatorList=" "
      generatorList=adjustl(Get_IT_Generators(aux))
      if (len_trim(generatorList) == 0) exit
      call Group_Constructor(generatorList,Grp)
      Grp%numspg=i
      call Get_SpaceGroup_Symbols(aux,HM=hm)
      Grp%spg_symb=hm(1:1)//l_case(hm(2:))
      if (Err_CFML%Ierr /= 0) then
          write(lun,'(i4,a)') i," "//trim(Err_CFML%Msg)
          cycle
      else
          call Identify_Group(Grp)
          if (Err_CFML%Ierr /= 0) then
              write(lun,'(/,a,i4,a)') "  Group number: ",i, &
                                      " => Error in the identification of the group: "//trim(Err_CFML%Msg)
          end if
          write(lun,"(/,a)") "  ------------------------------------------------------------------------------"
          write(lun,"(a,i4, a20,a)") "  Group number: ",i,"  Symb: "//trim(Grp%Spg_symb),&
                                     "  Transf. to standard: "//trim(Grp%mat2std_shu)
          write(lun,"(a)")   "  ------------------------------------------------------------------------------"
      end if

      call get_subgroups_subgen(Grp, sGrp, nsg)

      if (nsg > 0) Then
         write(lun,"(a,i4)") "  Total number of subgroups: ",nsg
         do L=1,nsg
            call Identify_Group(sGrp(L))
            if (Err_CFML%Ierr /= 0) then
               write(lun,'(/,a,i4,a)') " => Error in the identification of the sub-group: ",L,"  => "//trim(Err_CFML%Msg)
               write(*,'(a,i4,a)') " => Error in the identification of the sub-group: ",L,"  => "//trim(Err_CFML%Msg)
            end if !else
               write(lun,"(/,2(a,i4),a20,a)") "  Sub-Group Number #",L, " of index: ",Grp%multip/sGrp(L)%multip, &
               "  Symb: "//trim(sGrp(L)%BNS_symb),"  Transf. to standard: "//trim(sGrp(L)%mat2std_shu)//"    Generators: "//trim(sGrp(L)%generators_list)

               !> Write the coset decomposition of Grp with respect to the current subgroup
               call Get_Cosets(Grp,sGrp(L),cosets)
               nc=size(cosets)
               if (nc > 0) then
                  !     12345678901234
                  forma="(a,i3,a,    a)"
                  write(forma(9:12),"(i4)") nc
                  write(lun,forma) "  Coset decomposition of  G: "//trim(Grp%Spg_symb)//&
                                   "(",nc+1,") =  H("//trim(sGrp(L)%BNS_symb)//")  + ",("{"//&
                                   trim(Grp%Symb_Op(cosets(j)))//"} H + ",j=1,nc-1), &
                                   "{"//trim(Grp%Symb_Op(cosets(nc)))//"} H"
               end if
            !end if
         end do
      end if
      call CPU_TIME(fin)
      write(*,"(a,i4,t25,a,t45,a,t95,i4,a,t115,a,f8.3,a)") " Group number: ",i,"Symb: "//trim(Grp%Spg_symb),        &
                                                       "Transf. to standard: "//trim(Grp%mat2std_shu)//" with ",&
                                                       nsg," subgroups","CPU-time: ",fin-par," seconds"
   end do
   call CPU_TIME(fin)
   minuts=nint((fin-start)/60.0)
   secs=((fin-start)/60.0-real(minuts))*60.0
   write(*,"(/,a,i3,a,f9.5,a)")   "  CPU_Time for all calculations: ",minuts," minutes, ",secs," seconds"
   write(lun,"(/,a,i3,a,f9.5,a)") "  CPU_Time for all calculations: ",minuts," minutes, ",secs," seconds"

End Program test230_groups