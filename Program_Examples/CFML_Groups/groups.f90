   Program test_groups
     use CFML_Rational_Groups
     character(len=256)   :: generatorList
     type(Spg_Type)                   :: Grp
     type(Spg_Type), dimension(512)   :: sGrp
     integer :: i,j,L,nsg
     real :: start, fin
     call Init_Group(2048) !Maximum admissible multiplicity
     do
         write(*,'(/,a)',advance='no') "Introduce generators: "
         read(*,'(a)') generatorList
         if (len_trim(generatorList) == 0) exit
         call CPU_TIME(start)
         call Group_Constructor(generatorList,Grp)
         if (Err_group) then
            write(*,'(/,4x,a)') Err_group_Mess
         else
            call print_Group(Grp)
         end if
         call get_subgroups(Grp,sGrp,nsg)
         if(nsg > 0) Then
           do L=1,nsg
              write(*,"(/a,i3)") "  SUB-GROUP NUMBER #",L
              call print_Group(sGrp(L))
           end do
         end if
         call CPU_TIME(fin)
         write(*,"(a,f12.3,a)") "CPU_TIME for this calculation: ",fin-start," seconds"
     end do
   End Program test_groups