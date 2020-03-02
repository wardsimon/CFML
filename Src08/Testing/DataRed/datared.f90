!! --- New version of the program DataRed based in CrysFML08
!! ---

  Program DataRed

    Use CFML_GlobalDeps
    Use CFML_Maths,   only: Set_Eps_Math
    Use CFML_gSpaceGroups
    Use CFML_Metrics
    Use CFML_Reflections
    Use CFML_strings, only: File_Type, Reading_File
    Use Twin_Mod
    Use DataRed_Mod
    Use DataRed_Reflections
    Use CFML_Propagation_Vectors
    implicit none

    integer :: narg,len_cmdline,ier,lenf,idot,lun,i
    character(len=:),allocatable :: cmdline,filered
    character(len=256)           :: line, title
    logical :: esta

    type(RefList_Type)          :: Ref
    type(File_Type)             :: cfl_file
    type(Twin_Type)             :: Tw
    Type(SPG_Type)              :: SpG, twSpG
    Type(SuperSpaceGroup_Type)  :: SSpG
    Type(Cell_G_Type)           :: cell
    Type(Conditions_Type)       :: cond
    Type(kvect_info_Type)       :: kinf
    Type(Group_k_Type),dimension(:),allocatable :: Gk


    narg=COMMAND_ARGUMENT_COUNT()
    cmdline="                          "
    len_cmdline=0
    if(narg > 0) then
        call GET_COMMAND_ARGUMENT(1,cmdline)
        len_cmdline=len_trim(cmdline)
    else
        write(unit=*,fmt="(a)") " => The program DataRed should be invoked as: DataRed xxx.red or xxx.cfl ! "
        call CloseProgram()
    end if
    !Attempt to read the input *.red file
    call Header_Output()
    if(len_cmdline /= 0) then
        idot=index(cmdline,".")
        if(idot == 0) then
           filered=cmdline//".red"
        else
           filered=cmdline
        end if
        inquire(file=filered,exist=esta)

        if(.not. esta) then
          write(unit=*,fmt="(a)") " => File: "//filered//" not found! "
          call CloseProgram()

        else

          cfl_file=Reading_File(filered)
          if(Err_CFML%Ierr /= 0) then
            write(unit=*,fmt="(a)")  Err_CFML%Msg
            call CloseProgram()
          end if

          Call Read_DataRed_File(cfl_file,cond,kinf,tw,SpG,Cell)

          if(Err_CFML%Ierr /= 0) then
            write(unit=*,fmt="(a)")  Err_CFML%Msg
            call CloseProgram()
          else
            if(cond%eps_given) then
              call Set_Eps_Math(cond%epsg)
            else
              call Set_Eps_Math(0.001_cp)
            end if
          end if

          if(cond%prop .and. SpG%mag_type == 1) then !Setting up the propagation vectors groups
            allocate(Gk(kinf%nk))
            do i=1,kinf%nk
              call K_Star(kinf%kv(:,i),SpG,Gk(i),.true.)
            end do
          end if

          if(len_trim(cond%filhkl) /= 0) then  !Read Reflections and contruct RefList_Type object Ref
            if(cond%prop .and. SpG%mag_type == 1) then
              call Read_Reflections_File(cond%filhkl,cond,kinf,Cell,Ref,Gk)
            else
              call Read_Reflections_File(cond%filhkl,cond,kinf,Cell,Ref)
            end if
            if(Err_CFML%Ierr /= 0) then
              write(unit=*,fmt="(a)")  Err_CFML%Msg
              call CloseProgram()
            end if
          end if

          if(len_trim(cond%fileout) /= 0) then
             open(newunit=lun,file=trim(cond%fileout)//".out",status="replace",action="write")
             !Call Write_Conditions(cfl_file,cond,Cell,SpG,kinf,tw,cond%forma)
             Call Header_Output(lun)
             if(cond%prop .and. SpG%mag_type == 1) then
               Call Write_Conditions(cfl_file,cond,Cell,SpG,kinf,tw,cond%forma,lun,Gk)
             else
               Call Write_Conditions(cfl_file,cond,Cell,SpG,kinf,tw,cond%forma,lun)
             end if
             Call Write_Reflections(Ref,kinf,lun)
             Call Treat_Reflections(Ref,cond,cell,SpG,kinf,tw,lun)
          end if
        end if
    else
        write(unit=*,fmt="(a)") " => An input (xxx.red or xxx.cfl) datared file is needed !"
        call CloseProgram()
    End if
    !Now Read the Reflections File
    !call

    Contains
      Subroutine CloseProgram()
         character(len=1) :: ans
         write(unit=*,fmt="(a)")   " "
         write(unit=*,fmt="(a)")   " => Press <cr> to finish ...."
         read(unit=*,fmt="(a)") ans
         stop
      End Subroutine CloseProgram

  End Program DataRed
