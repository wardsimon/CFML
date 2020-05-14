!!----
!!----
!!----
!!----
 Program Test_SHX_CIF_CFL
    !---- Use Modules ----!
    use CFML_Globaldeps
    use CFML_Strings,      only: File_type, u_case, Get_extension
    use CFML_Metrics,      only: Cell_G_Type, Write_Crystal_Cell
    use CFML_gSpaceGroups, only: Spg_Type, SuperSpaceGroup_Type, Write_SpaceGroup_Info
    use CFML_Atoms,        only: AtList_Type, Write_Atom_List
    use CFML_IOForm

    !---- Local Variables ----!
    implicit none

    character(len=:), allocatable       :: ext
    character(len=512)                  :: fname,cmdline
    integer                             :: nlong,narg
    real(kind=cp)                       :: start, fin

    type(Cell_G_Type)                   :: Cell
    !type(Spg_Type)                      :: SpG
    type(SuperSpaceGroup_Type)          :: SpG
    type(AtList_Type)                   :: Atm

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !character(len=256)                  :: setting,ctr_code
    !character(len=256),dimension(26)    :: tctr_code
    !type(Cell_G_Type)                   :: Cell,Celln
    !!type(Spg_Type)                      :: Grp
    !type(File_type)                     :: flist
    !type(SuperSpaceGroup_Type)          :: Grp
    !!type(rational),   dimension(:,:),allocatable :: Mat
    !!character(len=40),dimension(:,:),allocatable :: matrix
    !integer :: i, j, L,k, d,Dd,nsg, ind, indexg, num_group, ier,mult,codini,narg,len_cmdline
    !real(kind=cp), dimension(:,:),allocatable :: orb,morb
    !real(kind=cp), dimension(3)               :: codes=1.0
    !real(kind=cp), dimension(:,:),allocatable :: codeT
    !integer,       dimension(:),  allocatable :: ptr

    !> Init
    narg=COMMAND_ARGUMENT_COUNT()
    cmdline=" "; nlong=0
    if (narg ==0) then
       write(unit=*,fmt='(/,a)',advance='no') " => Introduce the name of the file: "
       read(unit=*,fmt='(a)') fname
       if (len_trim(fname) <=0 ) call CloseProgram()
       cmdline=trim(fname)

    else
       call GET_COMMAND_ARGUMENT(1, cmdline)
    end if
    nlong=len_trim(cmdline)

    !> Start
    call CPU_TIME(start)

    !> Type of Files
    ext=get_extension(cmdline)
    call Read_Xtal_Structure(trim(cmdline), Cell, Spg, Atm)

    !> Print Information
    if (Err_CFML%Ierr == 0) then
       !> Cell Info
       call Write_Crystal_Cell(Cell)
       write(unit=*,fmt='(a)') " "

       !> SpaceGroup
       call Write_SpaceGroup_Info(SpG)
       write(unit=*,fmt='(a)') " "

       !> Atoms
       call Write_Atom_List(Atm)
       write(unit=*,fmt='(a)') " "
    else
       write(unit=*,fmt='(/,a)') " => ERROR: "//trim(Err_CFML%Msg)
    end if
    call CPU_TIME(fin)
    write(unit=*,fmt="(/,a,f12.3,a)") "CPU_TIME for this calculation: ",fin-start," seconds"


    !!!! TEST JRC
    !call Set_Eps_Math(0.0002_cp)
    !
    !call Readn_Set_Xtal_Structure(fname,Cell,Grp,Atm,"MAtm_std","CFL")!,file_list=flist) !,Iphase,Job_Info,file_list,CFrame)
    !if(Err_CFML%Ierr == 0) then
    !   !write(*,"(/,a,/)")  " => Content of the CFL-file: "//flist%Fname
    !   !do i=1,flist%nlines
    !   !   write(*,"(i6,a)") i,"    "//flist%line(i)%Str
    !   !end do
    !   call Write_Crystal_Cell(Cell)
    !   if(len_trim(Grp%setting) /= 0) then
    !     write(*,"(/,a)") " => Transformed Cell"
    !     if(Grp%D > 4) then
    !       i=index(Grp%setting,"d")
    !       setting=Grp%setting(1:d-2)//";0,0,0"
    !     else
    !       setting=Grp%setting
    !     end if
    !     call Change_Setting_Cell(Cell,setting,Celln)
    !     call Write_Crystal_Cell(Celln)
    !   end if
    !   call Write_SpaceGroup_Info(Grp)
    !
    !   i=index(fname,".")
    !   filename=fname(1:i)//"cif"
    !   call Write_Cif_Template(filename, Cell, Grp, Atm, 2, "Testing WriteCIF")
    !
    !   if(Atm%natoms > 0) then
    !      !First Check symmetry constraints in magnetic moments and Fourier coefficients
    !      !call Check_Symmetry_Constraints(Grp,Atm)
    !      write(*,"(//a,i5)") "  Number of atoms:",Atm%natoms
    !      call Write_Atom_List(Atm,SpG=Grp)
    !      !Calculate all atoms in the unit cell
    !      forma="(i5, f10.5,tr8, f10.5,i8)"
    !      formb="(a, i3,a,6f10.5,a)"
    !      write(unit=formb(4:4),fmt="(i1)") Grp%nk
    !      write(forma(5:5),"(i1)") Grp%d-1
    !      write(forma(16:16),"(i1)") Grp%d-1
    !      write(*,"(//a)") "  Orbits of atoms after applying constraints on moments:"
    !      write(*,"(  a)") "  ======================================================"
    !
    !
    !      do i=1,Atm%natoms
    !        !codini=1; codes=1.0
    !        call Get_moment_ctr(Atm%Atom(i)%x,Atm%Atom(i)%moment,Grp,codini,codes,ctr_code=ctr_code)!,Ipr=6)
    !        write(*,"(a,3f10.5,a)") " => Moment of atom "//trim(Atm%Atom(i)%Lab)//": ",Atm%Atom(i)%moment,"    CtrCode: "//trim(ctr_code)
    !        call Get_Orbit(Atm%Atom(i)%x,Grp,Mult,orb,Atm%Atom(i)%moment,morb,ptr)
    !        write(*,"(a)") " => Orbit of atom: "//trim(Atm%Atom(i)%Lab)
    !
    !        Select Case(Grp%d-1)
    !          Case(3)
    !            write(*,"(a)") "    N      X         Y         Z                 Mx        My       Mz      PointoOP"
    !          Case(4)
    !            write(*,"(a)") "    N     X1        X2        X3        X4                 M1        M2         M3        M4      PointoOP"
    !          Case(5)
    !            write(*,"(a)") "    N     X1        X2        X3        X4        X5                 M1        M2        M3        M4        M5      PointoOP"
    !          Case(6)
    !            write(*,"(a)") "    N     X1        X2        X3        X4        X5        X6                 M1        M2        M3        M4        M5        M6      PointoOP"
    !        End Select
    !
    !        do j=1,Mult
    !            write(*,forma) j,orb(:,j),morb(:,j),ptr(j)
    !        end do
    !       Select Type(at => Atm%Atom(i))
    !         class is (MAtm_Std_Type)
    !           write(*,"(a)") " => Modulation amplitudes of atom: "//trim(Atm%Atom(i)%Lab)
    !           if(allocated(CodeT)) deallocate(CodeT)
    !           allocate(CodeT(6,at%n_mc))
    !           CodeT=1.0
    !           call Get_TFourier_Ctr(At%x,At%Mcs(:,1:at%n_mc),codeT,Grp,codini,"M",ctr_code=tctr_code)
    !           do j=1,At%n_mc
    !             write(*,formb) "     Mcs: [",Grp%Q_coeff(:,j),"]",At%Mcs(:,j),"    CtrCode: "//trim(tctr_code(j))
    !           end do
    !           if(allocated(CodeT)) deallocate(CodeT)
    !           allocate(CodeT(6,at%n_dc))
    !           CodeT=1.0
    !           call Get_TFourier_Ctr(At%x,At%Dcs(:,1:at%n_dc),codeT,Grp,codini,"D",ctr_code=tctr_code)
    !           do j=1,At%n_dc
    !             write(*,formb) "     Dcs: [",Grp%Q_coeff(:,j),"]",At%Dcs(:,j),"    CtrCode: "//trim(tctr_code(j))
    !           end do
    !       end select
    !      end do
    !   end if
    !end if

 contains
    !!----
    !!---- CLOSEPROGRAM
    !!----
    !!---- 09/05/2020
    Subroutine CloseProgram()
       !---- Local Variables ----!
       character(len=1) :: ans

       write(unit=*,fmt="(a)")   " "
       write(unit=*,fmt="(a)")   " => Press <cr> to finish ..."
       read(unit=*,fmt="(a)") ans

       stop
    End Subroutine CloseProgram

End Program Test_SHX_CIF_CFL