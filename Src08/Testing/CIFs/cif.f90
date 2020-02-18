!!----
!!----
!!----
!!----
 Program Test_CIF_CFL
    !---- Use Modules ----!
    use CFML_Globaldeps
    !use CFML_Symmetry_Tables
    use CFML_Metrics
    use CFML_gSpaceGroups
    !use CFML_Rational
    use CFML_IOForm
    use CFML_Atoms
    !use CFML_Strings,only: pack_string
    implicit none

    character(len=256)                  :: fname
    character(len=256)                  :: setting
    character(len=25)                   :: forma
    character(len=5)                    :: aux
    type(Cell_G_Type)                   :: Cell,Celln
    !type(Spg_Type)                      :: Grp
    type(AtList_Type)                   :: Atm
    type(File_type)                     :: flist
    type(SuperSpaceGroup_Type)          :: Grp
    !type(rational),   dimension(:,:),allocatable :: Mat
    !character(len=40),dimension(:,:),allocatable :: matrix
    integer :: i, j, L,k, d,Dd,nsg, ind, indexg, num_group, ier
    real(kind=cp) :: start, fin

    do

       write(*,'(/,a)',advance='no') " => Introduce the name of the CFL file: "
       read(*,"(a)") fname
       if(len_trim(fname) == 0) exit
       call CPU_TIME(start)
       call Readn_Set_Xtal_Structure(fname,Cell,Grp,Atm,"MAtm_std","CFL")!,file_list=flist) !,Iphase,Job_Info,file_list,CFrame)
       if(Err_CFML%Ierr == 0) then
          !write(*,"(/,a,/)")  " => Content of the CFL-file: "//flist%Fname
          !do i=1,flist%nlines
          !   write(*,"(i6,a)") i,"    "//flist%line(i)%Str
          !end do
          call Write_Crystal_Cell(Cell)
          if(len_trim(Grp%setting) /= 0) then
            write(*,"(/,a)") " => Transformed Cell"
            if(Grp%D > 4) then
              i=index(Grp%setting,"d")
              setting=Grp%setting(1:d-2)//";0,0,0"
            else
              setting=Grp%setting
            end if
            call Change_Setting_Cell(Cell,setting,Celln)
            call Write_Crystal_Cell(Celln)
          end if
          call Write_SpaceGroup_Info(Grp)
          !if(Atm%natoms > 0) call Write_Info_Atom_List(Atm)

          if(Atm%natoms > 0) then
             write(*,"(a,i5)") "  Number of atoms:",Atm%natoms
             call Write_Atom_List(Atm)
          end if
       else
          write(*,'(/,a)') " => ERROR: "//trim(Err_CFML%Msg)
       end if
       call CPU_TIME(fin)
       write(*,"(/,a,f12.3,a)") "CPU_TIME for this calculation: ",fin-start," seconds"
   end do

End Program Test_CIF_CFL