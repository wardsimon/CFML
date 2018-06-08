SubModule (CFML_Diffraction_Patterns) ReadPat_ILL
 
 Contains
    !!--++
    !!--++ Subroutine Read_Pattern_D1A_D2B
    !!--++
    !!--++    Read a pattern for D1A, D2B
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Read_Pattern_D1A_D2B(Filename,Pat)
       !---- Arguments ----!
       character(len=*),        intent(in)  :: Filename      ! Path+Filename
       class(DiffPat_ILL_Type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=180)                           :: txt1
       integer                                      :: i, nlines, j, no, ier,i_dat
       integer, dimension(:), allocatable           :: iww
       real(kind=cp)                                :: rmoni, rmoniold, cnorm,step
       logical                                      :: info

       !> Init
       call clear_error()

       !> File exists?
       inquire(file=trim(filename),exist=info)
       if ( .not. info) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" The file "//trim(filename)//" doesn't exist"
          return
       end if
          
       !> Open File
       open(newunit=i_dat,file=trim(filename),status="old",action="read",position="rewind",iostat=ier)
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error opening the file: "//trim(filename)//" for reading!"
          return
       end if
       
       !> Title
       read(unit=i_dat,fmt="(a)",iostat=ier) txt1
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if
       pat%title=trim(txt1)

       !> Step
       read(unit=i_dat,fmt="(tr16,F8.3)",iostat=ier) step
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt="(F8.3)",iostat=ier)pat%xmin
       if (ier /= 0)then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML_Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt="(2f8.0)",iostat=ier) rmoni,rmoniold
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML_Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if
       pat%monitor=rmoni

       if (rmoniold < 1.0) then
          cnorm=1.00
          rmoniold=rmoni
       else
          cnorm=rmoni/rmoniold
       end if

       nlines = nint(18.0/step)
       pat%npts  = 10*nlines
       if (pat%npts <= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML_Msg=" Error in Intensity file, Number of points negative or zero!"
          
          close(unit=i_dat)
          return
       end if

       !> Allocating Pattern Object
       call Allocate_Diffraction_Pattern(pat)

       if(allocated(iww) ) deallocate(iww)
       allocate(iww(pat%npts))

       j=0
       do i=1,nlines
          read(unit=i_dat,fmt="(10(i2,f6.0))",iostat=ier)(iww(j+no),pat%y(j+no),no=1,10)
          if (ier /= 0) then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML_Msg=" Error in Intensity file, check your instr parameter!"
             
             close(unit=i_dat)
             return
          end if
          if(abs(pat%y(j+1)+1000.0) < 1.0e-03) exit
          j = j+10
       end do
       j=j-10
       pat%npts=j
       pat%xmax = pat%xmin+(pat%npts-1)*step
       do i=1,pat%npts
          !if (pat%y(i) <= 0.00001) pat%y(i) = 1.0
          if (iww(i) == 0) iww(i) = 1
          pat%sigma(i) = cnorm*abs(pat%y(i))/real(iww(i))
          pat%x(i)= pat%xmin+(i-1)*step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       
       close(unit=i_dat)
       
       return
    End Subroutine Read_Pattern_D1A_D2B
    
    !!--++
    !!--++ Subroutine Read_Pattern_D1A_D2B_OLD
    !!--++
    !!--++    Read a pattern for D1A, D2B (Old Format)
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Read_Pattern_D1A_D2B_OLD(Filename,Pat)
       !---- Arguments ----!
       character(len=*),        intent(in)  :: Filename      ! Path+Filename
       class(DiffPat_ILL_Type), intent(out) :: Pat
       
       !---- Local Variables ----!
       integer                                      :: ier,i,i_dat
       integer, dimension(:), allocatable           :: iww
       logical                                      :: info
       real                                         :: step

       !> Init
       call clear_error()

       !> File exists?
       inquire(file=trim(filename),exist=info)
       if ( .not. info) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" The file "//trim(filename)//" doesn't exist"
          return
       end if
          
       !> Open File
       open(newunit=i_dat,file=trim(filename),status="old",action="read",position="rewind",iostat=ier)
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error opening the file: "//trim(filename)//" for reading!"
          return
       end if
       

       read(unit=i_dat,fmt=*,iostat=ier)pat%xmin,step,pat%xmax
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if
       pat%title=" No title: data format -> old D1A"

       pat%npts = (pat%xmax-pat%xmin)/step+1.5
       if (pat%npts <= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, Number of points negative or zero!"
          
          close(unit=i_dat)
          return
       end if

       !> Allocating
       call Allocate_Diffraction_Pattern(pat)

       if(allocated(iww) ) deallocate(iww)
       allocate(iww(pat%npts))

       read(unit=i_dat,fmt="(10(i2,f6.0))",iostat=ier)(iww(i),pat%y(i),i=1,pat%npts)
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if

       do i=1,pat%npts
          !if (pat%y(i) <= 0.00001) pat%y(i) = 1.0
          if (iww(i) == 0) iww(i) = 1
          pat%sigma(i) = abs(pat%y(i))/real(iww(i))
          pat%x(i)= pat%xmin+(i-1)*step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       
       close(unit=i_dat)

       return
    End Subroutine Read_Pattern_D1A_D2B_Old
    
    !!--++
    !!--++ Subroutine Read_Pattern_D1B_D20
    !!--++
    !!--++    Read a pattern for D1B or D20
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Read_Pattern_D1B_D20(Filename,Pat)
       !---- Arguments ----!
       character(len=*),        intent(in)  :: Filename      ! Path+Filename
       class(DiffPat_ILL_Type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=180)                           :: line
       integer                                      :: i,j,ier,i_dat
       integer, dimension(:), allocatable           :: iww
       real(kind=cp)                                :: aux,step
       logical                                      :: info

       !> Init
       call clear_error()

       !> File exists?
       inquire(file=trim(filename),exist=info)
       if ( .not. info) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" The file "//trim(filename)//" doesn't exist"
          return
       end if
          
       !> Open File
       open(newunit=i_dat,file=trim(filename),status="old",action="read",position="rewind",iostat=ier)
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error opening the file: "//trim(filename)//" for reading!"
          return
       end if
       
       do i=1,3
          read(unit=i_dat,fmt="(a)", iostat=ier)line
          if (ier /= 0 )then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"
             
             close(unit=i_dat)
             return
          end if
          if( i == 2) pat%title=line
       end do
       pat%title=trim(pat%title)//" "//trim(line)

       read(unit=i_dat,fmt="(f13.0,tr10,f8.3,tr45,4f9.3)  ",iostat=ier) pat%monitor,pat%xmin,step,pat%tset,aux,pat%tsamp
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt="(i4)",iostat=ier) pat%npts
       if (ier /= 0 )then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if

       if (pat%npts <= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, Number of points negative or zero!"
          
          close(unit=i_dat)
          return
       end if
  
       !> Allocating memory
       call Allocate_Diffraction_Pattern(pat)

       if(allocated(iww)) deallocate(iww)
       allocate(iww(pat%npts))

       read(unit=i_dat,fmt="(10(i2,f8.0))",iostat=ier) (iww(j),pat%y(j),j=1,pat%npts)

       if (ier /= 0 )then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
       end if

       pat%xmax = pat%xmin+(pat%npts-1)*step

       read(unit=i_dat,fmt=*,iostat=ier)line
       if (ier /= 0 )then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter!"
          
          close(unit=i_dat)
          return
      end if

       do i=1,pat%npts
          !if (pat%y(i) <= 0.00001) pat%y(i) = 1.0
          if (iww(i) <= 0) iww(i) = 1
          pat%sigma(i) = abs(pat%y(i))/REAL(iww(i))
          pat%x(i)= pat%xmin+(i-1)*step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       
       close(unit=i_dat)
       
       return
    End Subroutine Read_Pattern_D1B_D20
    
End SubModule ReadPat_ILL 