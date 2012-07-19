  !!---- Example of program reading powder diffraction patterns and making
  !!---- "corrections" of bad behaving channels (cells) in the detector.
  !!---- The program assumes that the outliers are isolated in the sense that
  !!---- there are no more than two consecutive cells to be corrected.
  !!---- The program interpolates linearly between two proper adjacent points
  !!---- and generates a random number with a Poisson distribution to replace
  !!---- the bad intensity of the cell.

  Program correction_cells
    Use CFML_Diffraction_Patterns
    Use CFML_Math_General,     only: locate
    Use CFML_String_Utilities, only: getword
    Use CFML_Random_Generators,only: Random_Poisson
    type (diffraction_pattern_type) :: diffpat
    character(len=256)              :: dfile,filebuf,line
    character(len=20)               :: fmode
    character(len=*),parameter,dimension(22) :: inst_mode=(/"FREE","D1B","D20","NLS", &
      "G41","INSTRM5","D1A","D2B","3T2","G42","D1AOLD","D2BOLD","OLDD1A", "OLDD2B","DMC","HRPT", &
      "SOCABIM","XYSIGMA","GSAS","GSASTOF","PANALYTICAL","TIMEVARIABLE"/)
    integer, parameter, dimension(22) :: inst=(/0,3,3,4,5,5,6,6,6,6,1,1,1,1,8,8,9,10,12,14,13,11/)

    integer, parameter                :: maxp=200,i_buff=1,i_dat=2,i_out=3
    integer                           :: i, narg, nzone, instr, im,ip, genpoi, ier, &
                                         linum, nfiles
    integer,          dimension(maxp) :: ncl
    real                              :: ang, rint1,rint2,det,rinn
    character(len=20),dimension(maxp) :: argv
    logical :: buffer_given=.false., angles=.false.

    linum=0; nfiles=0
    narg=command_argument_count()
    Select Case(narg)

      Case(1)  !A buffer file with all the information is given
        call get_command_argument(1,filebuf)
        buffer_given=.true.
        open(unit=i_buff,file=trim(filebuf),status="old",action="read",position="rewind",iostat=ier)
        if(ier /= 0) then
          write(unit=*,fmt="(a)") " => Error opening the file: "//trim(filebuf)
          stop
        end if

      Case(3:)  !The name of the file to correct, instrument type and list of cell values
        call get_command_argument(1,dfile)
        call get_command_argument(2,fmode)
        nzone=narg-2
        call Check_cells()
        do i=1,nzone
          call get_command_argument(i+2,argv(i))
        end do

      Case default
        call usage()
        stop

    End Select

    write(unit=*,fmt='(a)') '   ----------------------------------------------------'
    write(unit=*,fmt='(a)') '   Program to correct bad cells in diffraction patterns'
    write(unit=*,fmt='(a)') '   ----------------------------------------------------'

    do
       if(buffer_given) then
         read(unit=i_buff,fmt="(a)",iostat=ier) line
         if(ier /= 0) exit
         linum = linum + 1
         line=adjustl(line)
         if(line(1:1)=="!" .or. line(1:1)=="#") cycle
         call Getword(line,argv,narg)
         if(narg < 3) then
           write(unit=*,fmt="(a,i4)") " => Error in buffer file at line: ", linum
           write(unit=*,fmt="(a)") " => The format of the lines should be: filedat instr_type n1 n2 .."
           stop
         end if
         dfile=argv(1); fmode=argv(2)
         nzone=narg-2
         call Check_cells()
       end if
       !Check if fmode contains a number or a keyword
       read(unit=fmode,fmt=*,iostat=ier) instr
       if(ier == 0) then  !Assign the proper keyword to fmode
         fmode=" "
         do i=1,22
           if(instr == inst(i)) then
             fmode=inst_mode(i)
             exit
           end if
         end do
         if(len_trim(fmode) == 0) then
           write(unit=*,fmt="(a,i4)") " => Illegal instrument format: ",instr
           stop
         end if
       end if
       call read_pattern(dfile,diffpat,fmode)
       if(Err_diffpatt) then
         write(unit=*,fmt="(a)") trim(err_diffpatt_mess)
         stop
       end if
       nfiles=nfiles + 1
       !reading the cells
       do i=1,nzone
         if(index(argv(i),".") /= 0) then
           angles=.true.
           exit
         end if
       end do
       if(angles) then
         do i=1,nzone
           read(unit=argv(i),fmt=*) ang
           ncl(i)=locate(diffpat%x,diffpat%npts,ang)+1
           !write(*,*) ncl(i),ang
         end do
       else
         do i=1,nzone
           read(unit=argv(i),fmt=*) ncl(i)
         end do
       end if
       do i=1,nzone
         j=ncl(i)
         ip=j+1
         im=j-1
         if(i /= 1) then   !Check that the points do not coincide with excluded cells
           if(ncl(i-1) == im) im = im-1
           if(ncl(i+1) == ip) ip = ip+1
         end if
         rint1=diffpat%y(im)
         rint2=diffpat%y(ip)
         det=rint1/diffpat%sigma(im)
         det=0.5*(det+rint2/diffpat%sigma(ip))
         rin=0.5*(rint2+rint1)
         rin=rin*det
         call Random_Poisson(rin, Genpoi)
         diffpat%y(j)=real(genpoi)/det
       end do
       line="corr_"//trim(dfile)
       write(unit=*,fmt="(a,50i5)") " => Corrected Cells: ",ncl(1:nzone)
       write(unit=*,fmt="(a)") " => Writing corrected file: "//trim(line)
       !call Write_Pattern_INSTRM5(line,diffpat,"var")
       call Write_Pattern_XYSig(line,diffpat)
       if(buffer_given) cycle
       exit
    end do
    if(buffer_given) then
      write(unit=*,fmt="(a,i4,a)") " => The buffer file treated contained ",linum," lines"
    end if
    write(unit=*,fmt="(a,i4)") " => Number of treated files: ",nfiles

  Contains

    Subroutine Check_cells()
      if( nzone > maxp) then
        write(unit=*,fmt="(a,i4)") " => Too many bad cells, maximum: ",maxp
        stop
      end if
      return
    End Subroutine Check_cells

    Subroutine usage()
      write(unit=*,fmt="(a)") " => Program usage:"
      write(unit=*,fmt="(a)") "  > corr_cells  filebuf  <cr>  "
      write(unit=*,fmt="(a)") "            or  "
      write(unit=*,fmt="(a)") "  > corr_cells  filedat instr_type n1 n2 ..  <cr>  "
      write(unit=*,fmt="(a)") "   where n1, n2, ... are the number of the cells or their angular values"
      write(unit=*,fmt="(a)") "   "
      write(unit=*,fmt="(a)") "  If a buffer file is given, each line should be of the form: filedat instr_time n1 n2 ... "
      return
    End Subroutine usage

  End Program correction_cells

