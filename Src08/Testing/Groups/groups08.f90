 Module Get_gSpG
   use CFML_gSpaceGroups
   implicit none
   private
   Public :: Set_gSpG
   contains
     Subroutine Set_gSpG(str,SpG,Mode,Setting)
       character(len=*),          intent(in) :: str
       Class(SpG_Type),           intent(out):: SpG
       Character(len=*),          intent(in) :: Mode
       Character(len=*),optional, intent(in) :: Setting
       ! --- Local Variables ---!
       integer :: i
       character(len=len(str)) :: st

       Select Case (trim(Mode))
         Case("SHUBN","SUPER")
           if(present(setting)) then
             call Set_SpaceGroup(Str,Mode,SpG,Setting=Setting)
           else
             call Set_SpaceGroup(Str,Mode,SpG)
           end if
         Case default
           call Set_SpaceGroup(Str,SpG)
       End Select
     End Subroutine Set_gSpG

 End Module Get_gSpG

!!----
!!----
!!----
!!----
 Program Test_Groups
    !---- Use Modules ----!
    use CFML_Globaldeps, only: cp, err_cfml
    use CFML_Symmetry_Tables
    use CFML_gSpaceGroups
    use Get_gSpG

    character(len=256)                  :: generatorList
    character(len=80)                   :: setting
    character(len=25)                   :: forma="(i5,tr2,a,   i4,a,i8)"
    character(len=5)                    :: aux
    type(Spg_Type)                      :: Grp
    type(Spg_Type), dimension(512)      :: sGrp
    integer, dimension(:,:),allocatable :: table
    integer, dimension(:,:),allocatable :: G
    integer, dimension(:),  allocatable :: ord
    integer :: i, j, L, nsg, ind, indexg, num_group, ier
    real(kind=cp) :: start, fin
    logical :: set_given, sup_given, datb_given

    !> Init
    call Set_Conditions_NumOp_EPS(4096) !Maximum admissible multiplicity

    do
       set_given=.false.; sup_given=.false.; datb_given=.false.
       write(*,'(/,a,/)') " => Examples of input to the next question:"
       write(*,'(a)') "1236   :: a,b,c;0,0,0                  <-- Shubnikov group type 1236 in standard setting"
       write(*,'(a)') "123    :: a,c,b;1/2,0,0                <-- Space group type 123 in non-standard setting"
       write(*,'(a)') "Pn'ma                                  <-- Shubnikov group type Pn'ma in standard setting"
       write(*,'(a)') "Pn'ma  :: -c,b,a;0,0,0                 <-- Shubnikov group type Pn'ma in the non-standard setting Pcmn'"
       write(*,'(a)') "B 2 C B                                <-- space #41 of standard symbol Aba2"
       write(*,'(a)') "13232  :: sup                          <-- SuperSpace group 13232 in standard setting"
       write(*,'(a)') "221    :: a2,-a1,a3,a4;0,0,1/2,0  sup  <-- SuperSpace group #221 in non-standard setting"
       write(*,'(a)') "Pnma(0,0,g)ss0                         <-- SuperSpace group of standard symbol"
       write(*,'(a)') "x,-y,z,t,1;x,y,z,t+1/2,-1              <-- The generators of a Magnetic SuperSpace group"
       write(*,'(a)') "x1,-x2,x3,x4,1;x1,x2,x3,x4+1/2,-1      <-- The generators of a Magnetic SuperSpace group"
       write(*,'(/,a)',advance='no') " => Introduce generators, HM or Hall symbol, or the number of a standard group as indicated above: "
       read(*,'(a)') generatorList
       if (len_trim(generatorList) == 0) exit

       !> Determine if it is a number
       aux=" "
       read(unit=generatorList,fmt=*,iostat=ier) num_group
       if (ier == 0) then
          !Now determine if setting is appearing in which case it uses data bases
          i=index(generatorList,"::")
          if(i /= 0) then
             j=index(generatorList,"sup")  !Superspace
             if(j /= 0) then
               sup_given=.true.
               if(index(generatorList,"a1") /= 0) then
                 setting= adjustl(generatorlist(i+2:j-1))
                 set_given=.true.
               end if
             else
               setting= adjustl(generatorlist(i+2:))
               set_given=.true.
             end if
          end if
       else
          i=index(generatorList,"'")
          j=index(generatorList,"_")
          if(i /= 0 .or. j /= 0) datb_given=.true.
       end if

       call CPU_TIME(start)
       if(sup_given) then
          if(set_given) then
             call Set_gSpG(generatorList,Grp,"SUPER",Setting)
          else
             call Set_gSpG(generatorList,Grp,"SUPER")
          end if
       else
          if(set_given ) then
             call Set_gSpG(generatorList,Grp,"SHUBN",Setting)
          else if (datb_given) then
             call Set_gSpG(generatorList,Grp,"SHUBN")
          else
             call Set_gSpG(generatorList,Grp,"GEN")
          end if
       end if
       if (Err_CFML%Ierr /= 0) then
          write(*,'(/,4x,a)') trim(Err_CFML%Msg)
          cycle
       else
          call Write_SpaceGroup_Info(Grp)
       end if

       !call Group_Constructor(generatorList,Grp)
       !if (Err_CFML%Ierr /= 0) then
       !   write(*,'(/,4x,a)') trim(Err_CFML%Msg)
       !   cycle
       !else
       !   call Identify_Group(Grp) !.false.
       !   if (Err_CFML%Ierr /= 0) then
       !      write(*,'(/,4x,"=> Error in the identification of the group: ",a)') trim(Err_CFML%Msg)
       !   end if
       !   call Write_SpaceGroup_Info(Grp)
       !end if

       do
          write(*,'(/,a)',advance='no') "Introduce the index of subgroups (if = 0, no restriction, if < 0 no calculation): "
          read(*,*,iostat=ier) indexg
          if(ier == 0) exit
       end do
       !> Testing Get_subgroups_cosets
       if (indexg == 0) then
          call get_subgroups_subgen(Grp,sGrp,nsg)

       else if(indexg > 0) then
          call get_subgroups_subgen(Grp,sGrp,nsg,indexg)
       else
          cycle
       end if

       if (nsg > 0) Then
          do L=1,nsg
             write(*,"(/2(a,i3))") "  SUB-GROUP NUMBER #",L, " of index: ",Grp%multip/sGrp(L)%multip
             call Identify_Group(sGrp(L)) !.false.
             if (Err_CFML%Ierr /= 0) then
                write(*,'(/,4x,"=> Error in the identification of the group: ",a)') trim(Err_CFML%Msg)
             end if
             call Write_SpaceGroup_Info(sGrp(L))
          end do
       end if

       call CPU_TIME(fin)
       write(*,"(a,f12.3,a)") "CPU_TIME for this calculation: ",fin-start," seconds"
    end do

End Program Test_Groups