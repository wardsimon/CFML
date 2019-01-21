Program test230_groups

    use CFML_Symmetry_tables, only: Get_Generators
    use CFML_Rational_Groups
    use CFML_Standard_Settings

    character(len=256)                  :: generatorList
    character(len=5)                    :: aux
    type(Spg_Type)                      :: Grp
    type(Spg_Type), dimension(4096)     :: sGrp
    integer :: i,j,L,nsg,ind,indexg,lun
    real :: start, fin,par

    call Init_Group(2048) !Maximum admissible multiplicity
    call CPU_TIME(start)
    open(newunit=lun,file="test_230groups.out",status="replace",action="write")
    do i=2,230
        call CPU_TIME(par)
        write(aux,"(i3)") i
        call Get_Generators(aux,generatorList)
        if (len_trim(generatorList) == 0) exit
        call Group_Constructor(generatorList,Grp)
        if (Err_group) then
            write(lun,'(i4,a)') i," "//Err_group_Mess
            cycle
        else
            call Identify_Group(Grp)
            if (err_std) then
                write(lun,'(/,a,i4,a)') "  Group number: ",i, " => Error in the identification of the group: "//trim(err_std_mess)
            end if
            write(lun,"(/,a)") "  ------------------------------------------------------------------------------"
            write(lun,"(a,i4, a20,a)") "  Group number: ",i,"  Symb: "//trim(Grp%spg_symb),"  Transf. to standard: "//trim(Grp%mat2std)
            write(lun,"(a)")   "  ------------------------------------------------------------------------------"
        end if
        call get_subgroups_cosets(Grp,sGrp,nsg)
        if(nsg > 0) Then
            write(lun,"(a,i4)") "  Total number of subgroups: ",nsg
            do L=1,nsg
              call Identify_Group(sGrp(L))
              if (err_std) then
                write(lun,'(a,i4,a)') " => Error in the identification of the sub-group: ",L,trim(err_std_mess)
              else
                write(lun,"(2(a,i4),a20,a)") "  Sub-Group Number #",L, " of index: ",Grp%multip/sGrp(L)%multip, &
                "  Symb: "//trim(sGrp(L)%spg_symb),"  Transf. to standard: "//trim(sGrp(L)%mat2std)
              end if
            end do
        end if
        call CPU_TIME(fin)
        write(*,"(a,i4, a20,a,i4,a,f12.3,a)") "  Group number: ",i,"  Symb: "//trim(Grp%spg_symb),"  Transf. to standard: "//trim(Grp%mat2std)//" with ",&
                                 nsg," subgroups CPU-time: ",fin-par," seconds"
        !write(*,"(a,f12.5,a)") "  CPU_Time for this calculation: ",fin-par," seconds"
    end do
    call CPU_TIME(fin)
    write(*,"(a,f12.5,a)") "  CPU_Time for all calculations: ",fin-start," seconds"
    write(lun,"(a,f12.5,a)") "  CPU_Time for all calculations: ",fin-start," seconds"

End Program test230_groups