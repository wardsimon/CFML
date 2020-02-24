!! --- New version of the program DataRed based in CrysFML08
!! ---

  Program DataRed
    Use Twin_Mod
    Use CFML_Reflections
    Use CFML_strings, only: File_Type, Reading_File
    Use DataRed_Mod
    Use Reading_Reflections

    integer :: narg,len_cmdline,ier,lenf,idot
    character(len=:),allocatable :: cmdline
    character(len=256)           :: line, title
    logical :: esta
    type(File_Type) :: cfl_file
    type(Twin_Type) :: Tw

    narg=COMMAND_ARGUMENT_COUNT()
    len_cmdline=0
    if(narg > 0) then
        call GET_COMMAND_ARGUMENT(1,cmdline)
        len_cmdline=len_trim(cmdline)
    else
        write(unit=*,fmt="(a)") " => The program DataRed should be invoked as: DataRed xxx.red ! "
        call CloseProgram()
    end if

    write(unit=*,fmt="(a)") "       ==============================="
    write(unit=*,fmt="(a)") "       DATA REDUCTION PROGRAM: DataRed"
    write(unit=*,fmt="(a)") "       ==============================="
    write(unit=*,fmt="(a)") "          JRC-ILL version:30-2-2020"
    write(unit=*,fmt="(a)") "       "

    !Attempt to read the input *.red file
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
          Call Process_DataRed_File(cfl_file,cond,tw)
        end if
    else
        write(unit=*,fmt="(a)") " => An input (xxx.red) datared file is needed !"
        call CloseProgram()
    End if
    !Now Read the Reflections File
    !call

    Contains
      Subroutine CloseProgram()
         write(unit=*,fmt="(a)")   " "
         write(unit=*,fmt="(a)")   "  => Press <cr> to finish ...."
         read(unit=*,fmt="(a)") ans
         stop
      End Subroutine CloseProgram

  End Program DataRed
