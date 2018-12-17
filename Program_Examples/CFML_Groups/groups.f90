   Program test_groups
     use CFML_Rational_Groups
     character(len=256)                  :: generatorList
                                                  !123456789012345678
     character(len=25)                   :: forma="(i5,tr2,a,   i4,a,i8)"
     type(Spg_Type)                      :: Grp
     type(Spg_Type), dimension(512)      :: sGrp
     integer, dimension(:,:),allocatable :: table
     integer, dimension(:,:),allocatable :: G
     integer, dimension(:),  allocatable :: ord
     integer :: i,j,L,nsg,ind,indexg
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
         write(*,'(/,a)',advance='no') "Introduce the index of subgroups (if = 0, no restriction): "
         read(*,*) indexg
         !call Get_Multiplication_Table(Grp%Op,table)
         !!do i=1,Grp%multip
         !!   write(*,"(128i4)") (Table(i,j),j=1,Grp%multip)
         !!end do
         !Call Get_subgroups_from_Table(Table(1:Grp%Numops,1:Grp%Numops),G,ord,nsg)
         !do i=1,nsg
         !   write(forma(11:13),"(i3)") ord(i)
         !   write(*,forma) i," {",G(1:ord(i),i),"   }",ord(i)
         !end do
         if(indexg == 0) then
           call get_subgroups(Grp,sGrp,nsg)
         else
           call get_subgroups(Grp,sGrp,nsg,indexg)
         end if
         if(nsg > 0) Then
           do L=1,nsg
              write(*,"(/2(a,i3))") "  SUB-GROUP NUMBER #",L, " of index: ",Grp%multip/sGrp(L)%multip
              call print_Group(sGrp(L))
           end do
         end if
         call CPU_TIME(fin)
         write(*,"(a,f12.3,a)") "CPU_TIME for this calculation: ",fin-start," seconds"
     end do
   End Program test_groups