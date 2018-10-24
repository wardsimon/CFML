   Program test_groups
     use CFML_Groups
     character(len=256)   :: generatorList
     !type(raw_spg_type)        :: Grp
     type(rational_spg_type)                   :: Grp
     !type(rational_spg_type), dimension(512)   :: sGrp
     integer :: i,j,L,nsg
     call Init_Group(2048) !Maximum admissible multiplicity
     do
         write(*,'(/,a)',advance='no') "Introduce generators: "
         read(*,'(a)') generatorList
         if (len_trim(generatorList) == 0) exit
         call Group_Constructor(generatorList,Grp)
         if (Err_group) then
            write(*,'(/,4x,a)') Err_group_Mess
         else
            call print_Group(Grp)
         end if
         !call get_subgroups(Grp,sGrp,nsg)
         !if(nsg > 0) Then
         !  do L=1,nsg
         !     call print_Group(sGrp(L))
         !  end do
         !end if
     end do
   End Program test_groups