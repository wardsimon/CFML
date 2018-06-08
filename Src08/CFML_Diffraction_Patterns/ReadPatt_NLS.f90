SubModule (CFML_Diffraction_Patterns) ReadPat_NLS

  Contains
   !!--++
   !!--++ Subroutine Read_Pattern_Nls
   !!--++
   !!--++    Read a pattern for NLS
   !!--++
   !!--++ Update: February - 2005
   !!
   Module Subroutine Read_Pattern_NLS(Filename,Pat)
      !---- Arguments ----!
      character(len=*),    intent(in)  :: Filename      ! Path+Filename
      class(DiffPat_Type), intent(out) :: Pat
       
      !---- Local Variables ----!
      character(len=132)                           :: aline
      integer                                      :: nlines,j,i,ier, no, i_dat
      logical                                      :: title_given, info
      real(kind=cp)                                :: step

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
       
      title_given=.false.

      do
         read(unit=i_dat,fmt="(a)") aline
         aline=adjustl(aline)
         if (.not. title_given) then
            Pat%title=trim(aline)
            title_given=.true.
         end if
         if (aline(1:1) == "!") cycle
         
         read(unit=aline,fmt=*,iostat=ier) pat%xmin, step, pat%xmax
         if (ier /= 0 ) then
            Err_CFML%state=.true.
            Err_CFML%Flag=2
            Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"
            
            close(unit=i_dat)
            return
         end if
         exit
      end do

      pat%npts = (pat%xmax-pat%xmin)/step+1.5
      if (pat%npts <= 0) then
         Err_CFML%state=.true.
         Err_CFML%Flag=2
         Err_CFML%Msg=" Error in Intensity file, Number of points negative or zero!"
         
         close(unit=i_dat)
         return
      end if

      nlines = pat%npts/10-1

      !> Allocating
      call Allocate_Diffraction_Pattern(pat)

      j = 0
      do i=1,nlines
         read(unit=i_dat,fmt="(10F8.0)",iostat=ier)(pat%y(j+no),no=1,10)
         if (ier /= 0 ) then
            Err_CFML%state=.true.
            Err_CFML%Flag=2
            Err_CFML%Msg=" Error in (NLS) Intensity file, check your instr parameter!1"
            
            close(unit=i_dat)
            return
         end if
         
         read(unit=i_dat,fmt="(10F8.0)",iostat=ier)(pat%sigma(j+no),no=1,10)
         if (ier /= 0 ) then
            Err_CFML%state=.true.
            Err_CFML%Flag=2
            Err_CFML%Msg=" Error in (NLS) Intensity file, check your instr parameter!2"
            
            close(unit=i_dat)
            return
         end if
         j = j+10
      end do

      pat%sigma(1) = pat%sigma(1)**2
      pat%x(1)=pat%xmin

      do i=2,pat%npts
         if (pat%sigma(i) < 0.00001) pat%sigma(i) = 1.0
         pat%sigma(i) = pat%sigma(i)**2.0
         pat%x(i)= pat%xmin+(i-1)*step
      end do
      pat%ymin=minval(pat%y(1:pat%npts))
      pat%ymax=maxval(pat%y(1:pat%npts))
      
      close(unit=i_dat)
      return
   End Subroutine Read_Pattern_NLS
    
End SubModule ReadPat_NLS