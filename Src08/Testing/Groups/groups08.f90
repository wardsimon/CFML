!!----
!!----
!!----
!!----
 Program Test_Groups
    !---- Use Modules ----!
    use CFML_Globaldeps, only: cp, err_cfml
    use CFML_Symmetry_Tables   
    use CFML_gSpaceGroups

    character(len=256)                  :: generatorList
    character(len=25)                   :: forma="(i5,tr2,a,   i4,a,i8)"
    character(len=5)                    :: aux
    type(Spg_Type)                      :: Grp
    type(Spg_Type), dimension(512)      :: sGrp
    integer, dimension(:,:),allocatable :: table
    integer, dimension(:,:),allocatable :: G
    integer, dimension(:),  allocatable :: ord
    integer :: i, j, L, nsg, ind, indexg, num_group, ier
    real(kind=cp) :: start, fin

    !> Init
    call Set_Conditions_NumOp_EPS(4096) !Maximum admissible multiplicity
    
    do
       write(*,'(/,a)',advance='no') "Introduce generators or number of a standard group: "
       read(*,'(a)') generatorList
       if (len_trim(generatorList) == 0) exit
        
       !> Determine if it is a number
       aux=" "
       read(unit=generatorList,fmt=*,iostat=ier) num_group
       if (ier == 0) then
          write(aux,"(i3)") num_group
          generatorList=" "
          generatorlist=Get_IT_Generators(aux)
       end if
       
       call CPU_TIME(start)

       call Group_Constructor(generatorList,Grp)
       if (Err_CFML%Ierr /= 0) then
          write(*,'(/,4x,a)') trim(Err_CFML%Msg)
          cycle
       
       else
          call Identify_Group(Grp) !.false.
          if (Err_CFML%Ierr /= 0) then
             write(*,'(/,4x,"=> Error in the identification of the group: ",a)') trim(Err_CFML%Msg)
          end if
          call Write_SpaceGroup_Info(Grp)
       end if
        
       write(*,'(/,a)',advance='no') "Introduce the index of subgroups (if = 0, no restriction, if < 0 no calculation): "
       read(*,*) indexg
       
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