Submodule (CFML_Diffraction_Patterns) WritePatterns

 Contains
    !!----
    !!---- SUBROUTINE WRITE_PATTERN_XYSIG
    !!----
    !!----    Write a pattern in X,Y,Sigma format
    !!----
    !!---- Update: March - 2007
    !!
    Module Subroutine Write_Pattern_XYSig(Filename,Pat,excl,xmin,xmax)
       !---- Arguments ----!
       character(len=*),               intent(in) :: filename     ! Path+Filename
       class(DiffPat_Type),            intent(in) :: Pat          ! Pattern object
       logical, dimension(:),optional, intent(in) :: excl         ! Exclusion zones
       real(kind=cp),        optional, intent(in) :: xmin         ! Limits
       real(kind=cp),        optional, intent(in) :: xmax    

       !---- Local Variables ----!
       integer                :: i,n, ini,ifin,npoi, ier, i_dat
       character(len=256)     :: excluded
       character(len=6)       :: cellexc
       
       !> Init
       call clear_error()
       
       open(newunit=i_dat,file=trim(filename),status="replace",action="write",iostat=ier)
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error opening the file: "//trim(filename)//" for writing!"
          return
       end if
       
       ini=1
       if (present(xmin)) then
          do i=1,Pat%npts
             if (pat%x(i) >= xmin) then
                ini=i
                exit
             end if
          end do
       end if
       
       ifin=Pat%npts
       if (present(xmax)) then
          do i=Pat%npts,1,-1
             if (pat%x(i) <= xmax) then
                ifin=i 
                exit
             end if
          end do
       end if
       
       npoi=ifin-ini+1
       excluded=" "
       if (present(excl)) then
          n=0
          do i=1,size(excl)
             if (excl(i)) then
                n=n+1
                write(unit=cellexc,fmt="(i6)") i
                excluded=trim(excluded)//cellexc
             end if
          end do
          npoi=npoi-n
       end if
       
       write(unit=i_dat,fmt="(a)")"XYDATA"
       write(unit=i_dat,fmt="(a)")"TITLE "//trim(pat%title)
       write(unit=i_dat,fmt="(a)")"COND: "//trim(pat%kindRad)//"-"//trim(pat%ScatVar)
       
       select type (Pat)
          type is (DiffPat_Type)
             write(unit=i_dat,fmt="(a)")"FILE: "//trim(filename)
             write(unit=i_dat,fmt="(a,i8)") "!N POINTS ",  npoi
             
             if (index(U_case(pat%ScatVar),"THET") /= 0) then
                write(unit=i_dat,fmt="(a,3f9.5)")"FILE: "//trim(filename)//"   Wavelengths: ",pat%wave(1:3)
       
             else if (index(U_case(pat%ScatVar),"TOF") /= 0) then
                write(unit=i_dat,fmt="(a,2f9.5)")"FILE: "//trim(filename)//"   TOF Dtt1, Dtt2: ",pat%wave(1:2)
             
             else
                write(unit=i_dat,fmt="(a)")"FILE: "//trim(filename)
             end if
             
          class is (DiffPat_ILL_Type)  
             write(unit=i_dat,fmt="(a,2f10.3)") "TEMP", pat%tsample, pat%tset
             if (pat%ct_step) then
                write(unit=i_dat,fmt="(a,2f8.4,i3,f8.5,a)") &
                      "INTER ", 1.0,1.0,2,pat%x(ini+1)-pat%x(ini)," <- internal multipliers for X, Y-Sigma, Interpol, StepIn"
             else
                write(unit=i_dat,fmt="(a,2f8.4,i3,f8.5,a)") &
                      "INTER ", 1.0,1.0,0,0.0," <- internal multipliers for X, Y-Sigma, Interpol, StepIn"
             end if
             
             write(unit=i_dat,fmt="(a,f12.2,i8)") "! MONITOR & N POINTS ", pat%monitor, npoi
       end select
       
       write(unit=i_dat,fmt="(a)") "! Scatt. Var., Profile Intensity, Standard Deviation "
       if (present(excl)) write(unit=i_dat,fmt="(a)") "! Excluded points (absent in the file):"//trim(excluded)
       write(unit=i_dat,fmt="(a,a10,a)") "!     ",pat%ScatVar,"        Y          Sigma "

       if (present(excl)) then
          if (pat%sigVar) then
             do i=ini,ifin
                if (excl(i)) cycle
                write(unit=i_dat,fmt="(3f14.5)") pat%x(i), pat%y(i),sqrt(pat%sigma(i))
             end do
         
          else
             do i=ini,ifin
                if (excl(i)) cycle
                write(unit=i_dat,fmt="(3f14.5)") pat%x(i), pat%y(i), pat%sigma(i)
             end do
          end if
       
       else
          if (pat%sig_var) then
             do i=ini,ifin
                write(unit=i_dat,fmt="(3f14.5)") pat%x(i),pat%y(i),sqrt(pat%sigma(i))
             end do
          else
             do i=ini,ifin
                write(unit=i_dat,fmt="(3f14.5)") pat%x(i),pat%y(i),pat%sigma(i)
             end do
          end if
       end if
       close(unit=i_dat)

       return
    End Subroutine Write_Pattern_XYSig
    
    !!----
    !!---- Subroutine Write_Pattern_FreeFormat(Filename,Pat,excl,xmin,xmax)
    !!----    character (len=*),               intent(in)     :: Filename
    !!----    type (diffraction_pattern_type), intent(in out) :: Pat
    !!----    logical, dimension(:),optional,  intent(in)     :: excl
    !!----    real,                 optional,  intent(in)     :: xmin,xmax
    !!----
    !!----    Write a pattern in Free Format (Instrm=0)
    !!----
    !!---- Update: 21/03/2011
    !!
    Module Subroutine Write_Pattern_FreeFormat(Filename,Pat,excl,xmin,xmax)
       !---- Arguments ----!
       character (len=*),               intent(in)     :: Filename      ! Path+Filename
       class(DiffPat_Type),             intent(in out) :: Pat           ! Pat object
       logical, dimension(:),optional,  intent(in)     :: excl
       real,                 optional,  intent(in)     :: xmin,xmax

       !---- Local Variables ----!
       integer   :: i,j,k,nl,ier,i_dat,ini,ifin, npoi,jmin,jmax

       !> Init
       call clear_error()
       
       open(newunit=i_dat,file=trim(filename),status="replace",action="write",iostat=ier)
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error opening the file: "//trim(filename)//" for writing!"
          return
       end if
       
       ini=1
       ifin=Pat%npts
       if (present(xmin)) then
          do i=1,Pat%npts
             if (pat%x(i) >= xmin) then
                ini=i
                exit
             end if
          end do
       end if
       if (present(xmax)) then
          do i=Pat%npts,1,-1
             if (pat%x(i) <= xmax) then
                ifin=i
                exit
             end if
          end do
       end if
       
       npoi=ifin-ini+1
       if (present(excl)) then !Replace the excluded points by the average of adjacent non-excluded points
          do i=ini+1,ifin-1
             if (excl(i)) then
                do j=max(i-1,ini),1,-1
                   if (.not. excl(j)) then
                      jmin=j
                      exit
                   end if
                end do
                do j=i+1,ifin
                   if (.not. excl(j)) then
                      jmax=j
                      exit
                   end if
                end do
                Pat%y(i)=0.5*(Pat%y(jmin)+Pat%y(jmax))
            end if
         end do
       end if

       write(unit=i_dat,fmt='(3(1x,f14.6),2x,a)') Pat%x(ini), Pat%x(ini+1)-Pat%x(ini), Pat%x(ifin), trim(pat%Title)
       nl=npoi/10
       if (mod(npoi,10) /= 0) nl=nl+1
       j=1
       do i=1,nl
          if (i /= nl) then
             write(unit=i_dat,fmt='(10i8)') nint(pat%y(j:j+9))
          else
             k=pat%npts - 10*(i-1)
             write(unit=i_dat,fmt='(10i8)') nint(pat%y(j:j+k-1))
          end if
          j=j+10
       end do

       close(unit=i_dat)
       return
    End Subroutine Write_Pattern_FreeFormat
    
    !!----
    !!---- SUBROUTINE WRITE_PATTERN_INSTRM5
    !!----
    !!----    Write a pattern in 2-axis format with fixed step (Instrm=5)
    !!----    The pattern Pat is modified on output if excl is present.
    !!----    Only the points starting with the next value higher than
    !!----    xmin (if present) are written
    !!----    If var is present the standard deviations are also provided,
    !!----    otherwise they are calculated from the number of counts and the
    !!----    values of the normalisation monitor and the used monitor.
    !!----
    !!---- Updated: 29/04/2011, 18/07/2012, 25/10/2015 (JRC)
    !!
    Module Subroutine Write_Pattern_INSTRM5(Filename,Pat,excl,xmin,xmax,var)
       !---- Arguments ----!
       character (len=*),               intent(in)     :: Filename
       class(DiffPat_ILL_Type),         intent(in out) :: Pat
       logical, dimension(:),optional,  intent(in)     :: excl
       real,                 optional,  intent(in)     :: xmin,xmax
       character (len=*), optional,     intent(in)     :: var

       !---- Local Variables ----!
       integer   :: i,j,ier,i_dat,ini,ifin, npoi,jmin,jmax
  
       !> Init
       call clear_error()
       
       opennew(unit=i_dat,file=trim(filename),status="replace",action="write",iostat=ier)
       if (ier /= 0 ) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          ERR_CFML%Msg=" Error opening the file: "//trim(filename)//" for writing!"
          return
       end if
       
       ifin=Pat%npts
       ini=1
       if (present(xmin)) then
          do i=1,ifin
             if (pat%x(i) >= xmin) then
                ini=i
                exit
             end if
          end do
       end if
       if (present(xmax)) then
          do i=Pat%npts,1,-1
             if (pat%x(i) <= xmax) then
                ifin=i
                exit
             end if
          end do
       end if
       
       npoi=ifin-ini+1
       if (present(excl)) then !Replace the excluded points by the average of adjacent non-excluded points
          do i=ini+1,ifin-1
             if (excl(i)) then
                do j=max(i-1,ini),1,-1
                   if (.not. excl(j)) then
                      jmin=j
                      exit
                   end if
                end do
                do j=i+1,ifin
                   if (.not. excl(j)) then
                      jmax=j
                      exit
                   end if
                end do
                Pat%y(i)=0.5*(Pat%y(jmin)+Pat%y(jmax))
             end if
          end do
       end if

       Write(unit=i_dat,fmt='(a)') trim(Pat%Title)
       Write(unit=i_dat,fmt="(a,f10.5)") trim(pat%kindRad)//" "//trim(pat%ScatVar)//", Wavelength (angstroms): ",pat%wave(1)
       Write(unit=i_dat,fmt="(a)")  "! Npoints TSample  TSetting Variance    Norm-Monitor           Monitor"
       if (Pat%ymax > 999999999.0) then
          write(unit=i_dat,fmt="(a,f16.1)") "!  Warning! Maximum number of counts above the allowed format ",Pat%ymax
          Err_CFML%state=.true.
          Err_CFML%Flag=1
          Err_CFML%Msg=" Too high counts ... format error in the file: "//trim(filename)//" at writing!"
       end if
       if (present(var)) then
          Write(unit=i_dat,fmt="(i6,tr1,2F10.3,i5,2f18.1)")  npoi, Pat%tsamp,Pat%tset, 1,&
                                                             Pat%Norm_Mon, Pat%Monitor
       else
          Write(unit=i_dat,fmt="(i6,tr1,2F10.3,i5,2f18.1)")  npoi, Pat%tsamp,Pat%tset, 0,&
                                                             Pat%Norm_Mon, Pat%Monitor
       end if
       Write(unit=i_dat,fmt="(3F12.5)")  Pat%x(ini),Pat%x(ini+1)-Pat%x(ini),Pat%x(ifin)
       Write(unit=i_dat,fmt="(8F14.2)")  Pat%y(ini:ifin)
       if (present(var)) then
          if (pat%sigvar) then
             Write(unit=i_dat,fmt="(8F14.2)")  sqrt(Pat%sigma(ini:ifin))
          else
             Write(unit=i_dat,fmt="(8F14.2)")  Pat%sigma(ini:ifin)
          end if
       end if
       
       close(unit=i_dat)
       
       return
    End Subroutine Write_Pattern_INSTRM5
 
End Submodule WritePatterns 