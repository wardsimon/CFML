SubModule (CFML_Diffraction_Patterns) ReadPat_LLB
 
 Contains
    !!--++
    !!--++ Subroutine Read_Pattern_G41(Filename,Pat)
    !!--++
    !!--++    Read a pattern for G41
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Read_Pattern_G41(Filename,Pat)
       !---- Arguments ----!
       character(len=*),        intent(in)  :: Filename      ! Path+Filename
       class(DiffPat_ILL_Type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=180)                           :: txt1, txt2, txt3
       integer                                      :: i, ier, ivari,i_dat
       real(kind=cp)                                :: cnorm,step
       real(kind=cp)                                :: rmon1, rmon2
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
       
       read(unit=i_dat,fmt="(A)",iostat=ier) txt1                  !1
       pat%title=txt1
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter! "
          
          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier) txt2                  !2
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter! "
          
          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier) txt3                  !3
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter! "
          
          close(unit=i_dat)
          return
       end if
       
       do
         read(unit=i_dat,fmt="(A)",iostat=ier) txt3
         txt3=adjustl(txt3)
         if(txt3(1:1) /= "!") exit
       end do
       read(unit=txt3, fmt=*, iostat=ier)  pat%npts,pat%tsamp,pat%tset,ivari,rmon1,rmon2
       if (ier /= 0 )then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter! "
          
          close(unit=i_dat)
          return
       end if
       pat%monitor=rmon1
       if (pat%npts <= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, Number of points negative or zero!"
          
          close(unit=i_dat)
          return
       end if

       !> Allocating
       call Allocate_Diffraction_Pattern(pat)


       read(unit=i_dat,fmt=*,iostat=ier)pat%xmin,step,pat%xmax              !5
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter! "
          
          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt=*,iostat=ier)(pat%y(i),i=1, pat%npts)
       if (ier /= 0 )then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, check your instr parameter! "
          
          close(unit=i_dat)
          return
       end if

       if (ivari /= 0) then          !IVARI
          read(unit=i_dat,fmt=*,iostat=ier)(pat%sigma(i),i=1, pat%npts)
          if (ier /= 0 ) then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML%Msg=" Error in Intensity file, check your instr parameter! "
          
             close(unit=i_dat)
             return
          end if
          
          cnorm=0.0
          do i=1,pat%npts
             pat%sigma(i)=pat%sigma(i)*pat%sigma(i)
             if (pat%sigma(i) < 0.000001) pat%sigma(i)=1.0
             pat%x(i)= pat%xmin+(i-1)*step
             cnorm=cnorm+pat%sigma(i)/pat%y(i)
          end do
          cnorm=cnorm/real(pat%npts)
       else                         !ivari
          if (rmon1 > 1.0 .and. rmon2 > 1.0) then
             cnorm=rmon1/rmon2
          else
             cnorm=1.0
          end if
          do i=1,pat%npts
             pat%sigma(i)=pat%y(i)*cnorm
             pat%x(i)= pat%xmin+(i-1)*step
          end do
       end if                        !IVARI
       
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       
       !> Close unit
       close(unit=i_dat)
       
       return
    End Subroutine Read_Pattern_G41
 
End SubModule ReadPat_LLB 