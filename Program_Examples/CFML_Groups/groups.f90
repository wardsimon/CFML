Program test_groups
    use CFML_Symmetry_tables, only: Get_Generators
    use CFML_Rational_Groups
    use CFML_Standard_Settings

    character(len=256)                  :: generatorList
                                              !123456789012345678
    character(len=25)                   :: forma="(i5,tr2,a,   i4,a,i8)"
    character(len=5)                    :: aux
    type(Spg_Type)                      :: Grp
    type(Spg_Type), dimension(512)      :: sGrp
    integer, dimension(:,:),allocatable :: table
    integer, dimension(:,:),allocatable :: G
    integer, dimension(:),  allocatable :: ord
    integer :: i,j,L,nsg,ind,indexg, num_group,ier
    real :: start, fin

    call Init_Group(4096) !Maximum admissible multiplicity
    do
        write(*,'(/,a)',advance='no') "Introduce generators or number of a standard group: "
        read(*,'(a)') generatorList
        if (len_trim(generatorList) == 0) exit
        !Determine if it is a number
        read(unit=generatorList,fmt=*,iostat=ier) num_group
        if(ier == 0) then
          write(aux,"(i3)") num_group
          generatorList=" "
          call Get_Generators(aux,generatorList)
        end if
        call CPU_TIME(start)

        call Group_Constructor(generatorList,Grp)
        if (Err_group) then
            write(*,'(/,4x,a)') Err_group_Mess
            cycle
        else
            call Identify_Group(Grp)
            if (err_std) then
                write(*,'(/,4x,"=> Error in the identification of the group: ",a)') trim(err_std_mess)
                !stop
            end if
            call print_Group(Grp)
        end if
        write(*,'(/,a)',advance='no') "Introduce the index of subgroups (if = 0, no restriction, if < 0 no calculation): "
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

        !Testing Get_subgroups_cosets

        if(indexg == 0) then
          call get_subgroups_cosets(Grp,sGrp,nsg)
        else if(indexg > 0) then
          call get_subgroups_cosets(Grp,sGrp,nsg,indexg)
        else
          cycle
        end if
        !cycle
        !if(indexg == 0) then
        !    call get_subgroups(Grp,sGrp,nsg)
        !else if(indexg > 0) then
        !    call get_subgroups(Grp,sGrp,nsg,indexg)
        !else
        !    cycle
        !end if

        if(nsg > 0) Then
            do L=1,nsg
                write(*,"(/2(a,i3))") "  SUB-GROUP NUMBER #",L, " of index: ",Grp%multip/sGrp(L)%multip
                call Identify_Group(sGrp(L))
                if (err_std) then
                    write(*,'(/,4x,"=> Error in the identification of the group: ",a)') trim(err_std_mess)
                else
                    call print_Group(sGrp(L))
                end if
            end do
        end if
        call CPU_TIME(fin)
        write(*,"(a,f12.3,a)") "CPU_TIME for this calculation: ",fin-start," seconds"
    end do

End Program test_groups