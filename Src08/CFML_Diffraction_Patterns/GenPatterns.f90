Submodule (CFML_Diffraction_Patterns) GenPatterns

 Contains
    !!----
    !!---- FUNCTION CALC_FWHM_PEAK
    !!----
    !!---- Function that calculate the FHWM of a peak situated on (xi,yi). Then
    !!---- the routine search the Ym value in the range (xi-rlim, xi+rlim) to
    !!---- obtain the FWHM.
    !!----
    !!---- The function return a negative values if an error is ocurred during calculation
    !!----
    !!---- Update: April - 2009
    !!
    Module Function Calc_FWHM_Peak(Pat, Xi, Yi, Ybi, RLim) Result(v)
       !---- Arguments ----!
       class(DiffPat_Type),       intent(in) :: Pat      ! Pattern object
       real(kind=cp),             intent(in) :: Xi       ! (Xi,Yi) for point i
       real(kind=cp),             intent(in) :: Yi       !
       real(kind=cp),             intent(in) :: Ybi      ! Background at point i
       real(kind=cp),optional,    intent(in) :: RLim     ! Limit range in X units to search the point
       real(kind=cp)                         :: V

       !---- Local variables ----!
       integer        :: j, i1, j1,n,nlim
       real(kind=cp)  :: xml, xmr, ym, x1, x2, y1, y2
       real(kind=cp)  :: difx


       !> Init value
       v=-1.0

       !> Y value for FHWM
       ym=0.5*(yi-ybi) + ybi

       !> Limit to search
       difx=pat%x(2)-pat%x(1)
       if (present(rlim)) then
          nlim=nint(rlim/difx)
       else
          nlim=nint(0.5/difx)     ! 0.5
       end if

       !> Locating the index that X(i1) <= Xi < X(i1+1)
       i1=0
       i1=locate(Pat%x,xi)
       if (i1 <=0 .or. i1 > Pat%npts) return  ! Error in the search

       !> Searching on Left side: Y(j1) <= ym < Y(j1+1)
       n=max(1,i1-nlim)
       j1=0
       do j=i1,n,-1
          if (pat%y(j) < ym) then
             j1=j
             exit
          end if
       end do
       if (j1 <= 0) j1=i1-1

       x1=Pat%x(j1)
       y1=Pat%y(j1)
       x2=Pat%x(j1+1)
       y2=Pat%y(j1+1)
       xml= x1 + ((ym-y1)/(y2-y1) )*(x2-x1)

       !> Searching on Right side: Y(j1) <= yn < Y(j1+1)
       n=min(i1+nlim,pat%npts)
       j1=0
       do j=i1,n
          if (pat%y(j) < ym) then
             j1=j
             exit
          end if
       end do
       if (j1 ==0) j1=i1

       x1=Pat%x(j1-1)
       y1=Pat%y(j1-1)
       x2=Pat%x(j1)
       y2=Pat%y(j1)
       xmr= x1 + ((ym-y1)/(y2-y1) )*(x2-x1)

       v=xmr-xml

       return
    End Function Calc_FWHM_Peak

    !!----
    !!---- Subroutine Add_Patterns
    !!----
    !!---- Add Patterns
    !!----
    !!---- Date: 25/03/2011
    !!
    Module Subroutine Add_Patterns(Patterns, N, Active, Pat, VNorm)
        !---- Arguments ----!
        class(DiffPat_Type), dimension(:), intent(in)  :: Patterns
        integer,                           intent(in)  :: N
        logical,             dimension(:), intent(in)  :: Active
        class(DiffPat_Type),               intent(out) :: Pat
        real(kind=cp), optional,           intent(in)  :: VNorm

        !---- Local Variables ----!
        integer                           :: i,j,k,npts,nc,np
        real(kind=cp)                     :: xmin,xmax,step,x1,x2,y,cnorm,fac
        real(kind=cp), dimension(:,:), allocatable :: d2y

        !> Init
        Pat%npts=0

        !> Checking
        if (N <= 0) return

        !> if (all(active) == .false.) return
        if (all(active) .eqv. .false.) return

        !> Initial values
        xmin=minval(Patterns(1:N)%xmin, mask= (active .eqv. .true.) )
        xmax=maxval(Patterns(1:N)%xmax, mask= (active .eqv. .true.) )

        npts=maxval(Patterns(1:N)%npts, mask= (active .eqv. .true.) )
        if (npts <= 1) then
           Err_CFML%state=.true.
           Err_CFML%Flag=2
           Err_CFML%Msg="Number of Points in the new Pattern was zero! "

           return
        end if

        step=(xmax-xmin)/real(npts-1)
        if (abs(step) <= 0.00001) then
           Err_CFML%state=.true.
           Err_CFML%Flag=2
           Err_CFML%Msg="Step size in the new Pattern was close to zero! "

           return
        end if

        !> Second Derivative
        if (allocated(d2y)) deallocate(d2y)
        allocate(d2y(npts,n))
        d2y=0.0

        do i=1,n
           if (.not. active(i)) cycle
           call second_derivative(Patterns(i)%x,Patterns(i)%y,Patterns(i)%npts,d2y(:,i))
        end do

        np=nint((xmax-xmin)/step)+1

        !> Allocating New Pat
        call Allocate_Diffraction_Pattern (Pat, np)

        if (present(vnorm)) then
           cnorm=vnorm
        else
           cnorm=maxval(Patterns(1:N)%ymax, mask= (active .eqv. .true.) )
        end if

        do j=1,np
           Pat%x(j)=xmin + (j-1)*step
           nc=0
           do i=1,N
              if (.not. active(i) ) cycle

              x1=minval(Patterns(i)%x)
              x2=maxval(Patterns(i)%x)
              k=locate(Patterns(i)%x,Pat%x(j))
              if (k == 0) cycle
              nc=nc+1
              y=splint(Patterns(i)%x,Patterns(i)%y,d2y(:,i),Patterns(i)%npts,Pat%x(j))
              fac=cnorm/Patterns(i)%ymax
              Pat%y(j)=Pat%y(j)+ y*fac
              Pat%sigma(j)=Pat%sigma(j)+ Patterns(i)%sigma(k)
           end do

           !> control
           if (nc > 0) then
              Pat%y(i)=Pat%y(i)/real(nc)
              Pat%sigma(i)=abs(Pat%sigma(i))/real(nc*nc)  ! No lo tengo muy claro
              !Pat%nd(i)=nc
           else
              Pat%y(i)=0.0
              Pat%sigma(i)=1.0
              !Pat%nd(i)=0
           end if
        end do

        !Pat%Monitor=cnorm
        Pat%xmin=xmin
        Pat%xmax=xmax
        !Pat%step=step
        Pat%ymin=minval(Pat%y)
        Pat%ymax=maxval(Pat%y)

        return
    End Subroutine Add_Patterns

    !!----
    !!---- SUBROUTINE CALC_BACKGROUND
    !!----
    !!----    Calculate a Background using an iterative process according
    !!----    to Bruckner, S. (2000). J. Appl. Cryst., 33, 977-979.
    !!----
    !!----    Np is the extension of Np points at Left and Right. Normally it could be
    !!----    around 10-40 points.
    !!----
    !!---- Update: December - 2008
    !!
    Module Subroutine Calc_BackGround(Pat, Ncyc, Np, Xmin, Xmax)
       !---- Arguments ----!
       class(DiffPat_E_Type),   intent(in out) :: Pat        ! Pattern object
       integer,                   intent(in)     :: NCyc       ! Number of Cycles to apply
       integer,                   intent(in)     :: Np         ! Number of extension points at L and R.
       real(kind=cp), optional,   intent(in)     :: Xmin       ! Min, max values in X
       real(kind=cp), optional,   intent(in)     :: Xmax

       !---- Variables ----!
       integer                                 :: n,n_ini,n_fin
       integer                                 :: i,j,k,ind1,ind2,nt
       real(kind=cp)                           :: x_ini,x_fin, yc_min, yc_max, yc_ave
       real(kind=cp),dimension(:), allocatable :: yc, yb

       !> Init
       call clear_error()
       pat%bgr=0.0_cp

       !> Check Pattern
       if (pat%npts <= 1) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          err_CFML%Msg="There aren't points in the Pattern used for Background calculation!"
          return
       end if

       !> Number of points into the range
       x_ini=pat%xmin
       x_fin=pat%xmax
       if (present(xmin)) x_ini=xmin
       if (present(xmax)) x_fin=xmax

       nt=0
       do i=1,pat%npts
          if (pat%x(i) <= x_ini) cycle
          if (pat%x(i) > x_fin) cycle
          nt=nt+1
       end do
       if (nt < 1) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          err_CFML%Msg="There aren't background points into the range"
          return
       end if

       !> Locating index that define the range to study
       ind1=0
       if (abs(x_ini-pat%xmin) <= eps) then
          ind1=1
       else
          ind1=locate(pat%x,x_ini)
          ind1=max(ind1,1)
          ind1=min(ind1,pat%npts)
       end if

       ind2=0
       if (abs(x_fin-pat%xmax) <= eps) then
          ind2=pat%npts
       else
          ind2=locate(pat%x,x_fin)
          ind2=min(ind2,pat%npts)
          ind2=max(ind2,1)
       end if

       if (ind1 == ind2) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          err_CFML%Msg="Same value for Xmin and Xmax for Background calculation"
          return
       end if
       if (ind1 > ind2) then
          i=ind1
          ind1=ind2
          ind2=i
       end if


       !> Allocating arrays
       if (ind2-ind1+1 > nt) nt=ind2-ind1+1
       if (allocated(yc)) deallocate(yc)
       if (allocated(yb)) deallocate(yb)
       allocate(yc(nt+2*np))
       allocate(yb(nt+2*np))

       yc=0.0_cp

       !> Load initial values
       n_ini=np+1
       n_fin=np+nt
       yc(1:np)=pat%y(ind1)
       yc(n_ini:n_fin)=pat%y(ind1:ind2)
       yc(n_fin+1:n_fin+np)=pat%y(ind2)

       yc_min=minval(pat%y(ind1:ind2))
       yc_ave=sum(pat%y(ind1:ind2))/real(nt)
       yc_max=yc_ave+2.0_cp*(yc_ave-yc_min)
       where(yc > yc_max) yc=yc_max

       !> Main cycles
       do n=1,Ncyc
          yb=0.0_cp
          do k=n_ini,n_fin ! Points Observed
             do i=-np,np
                if (i == 0) cycle
                j=k+i
                yb(k)=yb(k)+yc(j)
             end do
             yb(k)=yb(k)/real(2*np)
          end do
          do k=n_ini,n_fin
             j=k-np+ind1-1
             if (yb(k) > pat%y(j)) yb(k)=pat%y(j)
          end do
          yb(1:np)=yb(n_ini)
          yb(n_fin+1:n_fin+np)=yb(n_fin)
          yc=yb
       end do

       !> Save the result
       pat%bgr(ind1:ind2)=yc(n_ini:n_fin)

       return
    End Subroutine Calc_BackGround

    !!----
    !!---- SUBROUTINE DEL_NOISYPOINTS
    !!----
    !!---- Delete noisy points in a Pattern.
    !!----
    !!---- If FileInfo is .true. then a file is created containing
    !!---- information about the elimination of noisy points
    !!----
    !!---- Date: 26/03/2011
    !!
    Module Subroutine Del_NoisyPoints(Pat, NoisyP, FileInfo)
        !---- Arguments ----!
        class(DiffPat_Type),  intent(in out) :: Pat        ! Pattern object
        integer,              intent(out)    :: NoisyP     ! Noisy points
        logical, optional,    intent(in)     :: FileInfo   ! .true. For create an information file

        !---- Local Variables ----!
        logical                                  :: info
        integer                                  :: i,j,nomo1,nomo2,lun
        real(kind=cp), dimension(5)              :: cc
        real(kind=cp)                            :: suma,sc,dif1,dif2
        real(kind=cp)                            :: ci2,ci1,c,cd1,cd2,cn
        real(kind=cp), dimension(:), allocatable :: yc

        !> Initializing errors
        call clear_error()

        info=.false.
        if (present(FileInfo)) info=FileInfo

        NoisyP=0

        !> Check Pattern
        if (pat%npts < 1) then
           err_CFML%state=.true.
           err_CFML%Flag=2
           err_CFML%Msg=" There aren't points in the Pattern passed to Del_NoisyPoints"
           return
        end if

        if (info) then
           open(newunit=lun, file='Noisy_Points_Information.txt')
           write(unit=lun,fmt='(a/)')  " => Analysis of Noisy points of Pattern "//trim(Pat%title)
           write(unit=lun,fmt='(/a/)') " => A Noisy point means the following:"
           write(unit=lun,fmt='(a/)')  "        NoMono .and. Iosci = .true. "
           write(unit=lun,fmt='(a/)')  " where:"
           write(unit=lun,fmt='(a)')   "     ci2 : counts at          left-left position"
           write(unit=lun,fmt='(a)')   "     ci1 : counts at               left position"
           write(unit=lun,fmt='(a)')   "     cc  : counts at            current position"
           write(unit=lun,fmt='(a)')   "     cd1 : counts at              right position"
           write(unit=lun,fmt='(a)')   "     cd2 : counts at              right position"
           write(unit=lun,fmt='(a)')   "     sc  : 8.0*sqrt((ci1+ci2+cd1+cd2)/4.0)"
           write(unit=lun,fmt='(a)')   "     dif1: cc -2.0*ci1+ci2"
           write(unit=lun,fmt='(a)')   "     dif2: cc -2.0*cd1+cd2"
           write(unit=lun,fmt='(a)')   "    Iosci: .not.(dif1 < sc .or. dif2 < sc)"
           write(unit=lun,fmt='(a)')   "   NoMono: Non monotonic ci2,ci1,cc,cd1,cd2"
           write(unit=lun,fmt='(a/)')  "  cc(new): 0.5*(ci1+cd1)"
        end if

        !> Copy Y values
        if (allocated(yc)) deallocate(yc)
        allocate(yc(Pat%npts))
        yc=Pat%y

        cyc_1: do j=3,Pat%npts-2
           suma=0.0
           do i=1,5
              cc(i)=yc(i+j-3)
              if (cc(i) <= 0.0) cycle cyc_1
              if (i /= 3) suma=suma+cc(i)
           end do
           nomo1=0
           nomo2=0
           do i=2,5
              if (cc(i) > cc(i-1)) nomo1=nomo1+1
              if (cc(i) < cc(i-1)) nomo2=nomo2+1
           end do
           if (nomo1 == 4 .or. nomo2 == 4) cycle cyc_1

           sc=4.0_cp*sqrt(suma)
           dif1=cc(3)-2.0_cp*cc(1)+cc(2)
           dif2=cc(3)-2.0_cp*cc(4)+cc(5)
           if (.not. (dif1 <= sc .or. dif2 <= sc) ) then
              if (info) then
                 ci2=cc(1)
                 ci1=cc(2)
                 c=cc(3)
                 cd1=cc(4)
                 cd2=cc(5)
                 cn=0.5*(ci1+cd1)
                 write(unit=lun,fmt='(a,2i6,2(a,i6),a,a,2i6)') "  Counts-left: ",nint(ci2),nint(ci1), &
                                                               " Counts: ",nint(c)," (",nint(cn),")",  &
                                                               " Counts-right: ",nint(cd1),nint(cd2)
              end if
              noisyp=noisyp+1
              yc(j)=0.5*(yc(j-1)+yc(j+1))
           end if
        end do cyc_1

        if (info) then
           select case (noisyP)
              case (0)
                 write(unit=lun,fmt='(/a)')  " => No noisy points were found for this Pattern!"
              case (1)
                 write(unit=lun,fmt='(/a)')  " => Only one noisy point was found for this Pattern!"
              case (2:)
                 write(unit=lun,fmt='(/a,i3,a)')  " => A ",noisyP," noisy points were found for this Pattern!"
           end select
           close(unit=lun)
        end if

        Pat%y=yc

        return
    End Subroutine Del_NoisyPoints

    !!----
    !!---- SUBROUTINE READ_BACKGOUND_FILE
    !!----
    !!----    Read background pattern from an external file.
    !!----
    !!----    Mode:
    !!----         Poly | Inter
    !!----
    !!----
    !!---- Update: February - 2005
    !!
    Module Subroutine Read_Background_File(Bck_File, Bck_Mode, Pat)
       !---- Arguments ----!
       character(len=*),         intent(in   )    :: bck_file      ! Path+Filename of Background file
       character(len=*),         intent(in   )    :: bck_mode
       class(DiffPat_E_Type),  intent(in out)   :: Pat

       !---- local variables ----!
       logical                                       :: esta
       character(len=80)                             :: line
       character(len=3)                              :: car
       integer                                       :: bck_points
       integer                                       :: i,j,i_bck
       integer                                       :: ier, alloc_error
       real(kind=cp), dimension (:), allocatable     :: bck_v
       real(kind=cp), dimension (:), allocatable     :: bck_p

       !> Init
       call clear_error()

       inquire(file=trim(bck_file), exist =esta)
       if (.not. esta) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg="The file "//trim(bck_file)//" doesn't exist"
          return
       end if

       !> Open Backgriund file
       open(newunit=i_bck,file=trim(bck_file),status="old",action="read",position="rewind",iostat=ier)
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg="Error opening the Background file: "//trim(bck_file)
          return
       end if

       !> Number of background points
       i=0
       do
          read(unit=i_bck,fmt="(a)",iostat=ier) line
          if (ier /= 0) exit

          if (len_trim(line) == 0) cycle
          if (index(line,"!") /= 0) cycle
          if (index(line,"#") /= 0) cycle
          i=i+1
       end do
       bck_points=i

       if (bck_points <=0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg="Was impossible to read any background point in the file: "//trim(bck_file)

          close (unit=i_bck)
          return
       end if

       !> Allocating variables
       if (allocated(bck_v)) deallocate(bck_v)
       allocate(bck_v(bck_points+1))
       bck_v=0.0_cp

       if (allocated(bck_p)) deallocate(bck_p)
       allocate(bck_p(bck_points+1))
       bck_p=0.0_cp

       rewind(unit=i_bck)
       do j=1, bck_points
          read(unit=i_bck,fmt="(a)",iostat=ier) line
          if (ier /= 0) exit

          if (len_trim(line) == 0 .or. line(1:1) == "!" .or. line(1:1)=="#") cycle

          read(unit=line, fmt=*, iostat=ier)  bck_p(j), bck_v(j)
          if (ier /= 0) then
             bck_points=j-1
             close(unit=i_bck)

             Err_CFML%state=.true.
             Err_CFML%Flag=1
             ERR_CFML%Msg=" WARNING: The reading of Background points was incomplete. Please, check it!"
             exit
          end if
       end do
       if (bck_points <=0) then
          Err_CFML%Flag=2
          Err_CFML%Msg="Was impossible to read any background point in the file: "//trim(bck_file)
          return
       end if

       close (unit=i_bck)
       car=adjustl(u_case(bck_mode))

       select case (car)
          case ("POL") ! Polynomial
             call set_background_poly (Pat, 50.0_cp, bck_p, bck_points )

          case ("INT") ! Interpolation
             call  set_background_inter (Pat, bck_v, bck_p, bck_points )

          case default
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             ERR_CFML%Msg=" Set Polynomial or Interpolation in the MODE for Read_Background procedure"
             return
       end select

       return
    End Subroutine Read_Background_File

    !!--++
    !!--++ SUBROUTINE SET_BACKGROUND_POLY
    !!--++
    !!--++    (PRIVATE)
    !!--++    Define a n-polynomial with constanat vlue at bkpos Background
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Set_Background_Poly(Pat, Bkpos, Bckx, N)
       !---- Arguments ----!
       class(DiffPat_E_Type),         intent(in out) :: Pat
       real (kind=cp),                intent(in    ) :: bkpos
       real (kind=cp), dimension(:),  intent(in    ) :: bckx
       integer,                       intent(in    ) :: n

       !---- Local Variables ----!
       integer :: i,j

       if (allocated(pat%bgr)) deallocate(pat%bgr)
       allocate(pat%bgr(pat%npts))
       pat%bgr=0.0_cp
       al_bgr=.true.

       do i=1, pat%npts
          do j=1,n
             pat%bgr(i)= pat%bgr(i) + bckx(j)*((pat%x(i)/bkpos-1.0)**(j-1))
          end do
       end do

       return
    End Subroutine Set_Background_Poly

    !!--++
    !!--++ SUBROUTINE SET_BACKGROUND_INTER
    !!--++
    !!--++    (PRIVATE)
    !!--++    Define a Background
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Set_Background_Inter(Pat, Bcky, Bckx, N)
       !---- Arguments ----!
       class(DiffPat_E_Type),       intent(in out) :: Pat
       real (kind=cp), dimension(:),  intent(in out) :: bcky
       real (kind=cp), dimension(:),  intent(in out) :: bckx
       integer,                       intent(in    ) :: n

       !---- Local variables ----!
       integer        :: nbx, nbac1 , i , j  , nxx
       real(kind=cp)  :: difl, difr , thx , delt, slope, bstep,p,step

       nbx=1
       nbac1=n

       difl=bckx(1)-pat%xmin
       difr=bckx(n)-pat%xmax

       if (difl >= 0) then
          if (pat%ct_step) then
             step=pat%x(2)-pat%x(1)
             nbx=difl/step + 1.5
          else
             nbx=locate(pat%x,bckx(1))
             if (nbx <= 0) nbx=1
          end if
          do i=1,nbx
             pat%bgr(i)=bcky(1)
          end do
       end if

       if (difr <= 0) then
          nbac1=n+1
          bckx(nbac1)=pat%xmax
          bcky(nbac1)=bcky(n)
       end if

       nxx=2
       do_i: do i=nbx,pat%npts
          thx=pat%x(i)
          do j=nxx,nbac1
             delt=bckx(j)-thx
             if (delt > 0.0) then
                p=bckx(j)-bckx(j-1)
                if (abs(p) > eps) then
                   slope=(bcky(j)-bcky(j-1))/p
                else
                   slope=0.0
                end if
                bstep=(thx-bckx(j-1))*slope
                pat%bgr(i)=bcky(j-1)+bstep
                nxx=j-1
                cycle do_i
             end if
          end do
       end do  do_i

       return
    End Subroutine Set_Background_Inter

    !!--++
    !!--++ Subroutine Read_Pattern_CIF
    !!--++
    !!--++    (PRIVATE)
    !!--++    Read a pattern from a CIF file
    !!--++
    !!--++ Update: November - 2017
    !!
    Module Subroutine Read_Pattern_CIF(Filename,Pat)
       !---- Arguments ----!
       character(len=*),    intent(in)  :: Filename      ! Path+Filename
       class(DiffPat_Type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=120),dimension(:),allocatable  :: file_lines
       character(len=1)                             :: aux
       integer                                      :: i, n, nlines, ic, ier,i_dat
       real(kind=cp), dimension(6)                  :: values, std
       real(kind=cp)                                :: chi2
       integer,       dimension(6)                  :: pos
       integer                                      :: line_block_id, line_probe, line_loop, line_point_id, line_npoints, line_start_data
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

       !> Init
       nlines=0
       do
         read(unit=i_dat,fmt="(a)",iostat=ier) aux
         if (ier /= 0) exit
         nlines=nlines+1
       end do
       if (nlines < 10) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Number of lines too small to hold a diffraction pattern, check your CIF file!"

          close(unit=i_dat)
          return
       end if

       allocate(file_lines(nlines))
       line_block_id=0; line_probe=0; line_loop=0; line_point_id=0; line_npoints=0; line_start_data=0

       rewind(unit=i_dat)
       do i=1,nlines
          read(unit=i_dat,fmt="(a)",iostat=ier) file_lines(i)

          if (index(file_lines(i),"_pd_block_id") /= 0) line_block_id=i
          if (index(file_lines(i),"_diffrn_radiation_probe") /= 0) line_probe=i
          if (index(file_lines(i),"loop_") /= 0 .and. line_loop == 0) line_loop=i
          if (index(file_lines(i),"_pd_proc_number_of_points") /= 0 ) line_npoints=i
       end do
       close(unit=i_dat)

       if (line_npoints == 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" No line with the number of points in the pattern, check your CIF file!"

          return
       else
          read(unit=file_lines(line_npoints)(27:),fmt=*,iostat=ier) pat%npts
          if (ier /= 0) then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML%Msg=" Error reading the number of points in the pattern, check your CIF file!"

             return
          end if
       end if

       !> Allocating
       call Allocate_Diffraction_Pattern(pat)

       if (line_block_id > 0) pat%title=file_lines(line_block_id)(13:)
       if (line_probe > 0)    pat%kindRad = adjustl(file_lines(line_probe)(24:))

       n=0
       pos=0
       do i=line_loop+1, line_loop+10
          n=n+1
          if (index(file_lines(i),"_pd_proc_point_id") /= 0) then
             pos(1)=n
             cycle
          end if
          if (index(file_lines(i),"_pd_proc_2theta_corrected") /= 0) then
             pat%scatvar =  "2theta"
             !pat%xax_text =  "2theta(degrees)"
             pos(2)=n
             cycle
          end if
          if (index(file_lines(i),"_pd_proc_d_spacing") /= 0) then
             pat%scatvar =  "d-spacing"
             !pat%xax_text =  "d-spacing(Angstroms)"
             pat%kindRad=  "Neutrons (TOF)"
             pos(2)=n
             cycle
          end if
          if (index(file_lines(i),"_pd_proc_intensity_total") /= 0) then
             pos(3)=n
             cycle
          end if
          if (index(file_lines(i),"_pd_calc_intensity_total") /= 0) then
             pos(4)=n
             cycle
          end if
          if (index(file_lines(i),"_pd_proc_intensity_bkg_calc") /= 0) then
             pos(5)=n
             cycle
          end if
          if (len_trim(file_lines(i)) == 0) then
             line_start_data=i+1
             exit
          end if
       end do

       if (pos(2) == 0 .or. pos(3) == 0) then  !Error experimental data are lacking
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error: experimental data are lacking, check your CIF file!"

          return
       end if

       n=0
       do i=line_start_data, nlines
          n=n+1
          if (n > pat%npts) exit
          call Get_NumStd(file_lines(i), values, std, ic)
          if (ic < 2) exit

          pat%x(n) =  values(pos(2))
          pat%y(n) =  values(pos(3))
          pat%sigma(n) = std(pos(3))
       end do

       pat%xmin=pat%x(1)
       pat%xmax=pat%x(pat%npts)
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       select type (Pat)
          class is (DiffPat_E_Type)
             n=0
             chi2=0.0
             do i=line_start_data, nlines
                n=n+1
                if (n > pat%npts) exit
                call Get_NumStd(file_lines(i), values, std, ic)
                if (ic < 2) exit

                if (pos(4) /= 0) pat%ycalc(n) = values(pos(4))
                if (pos(5) /= 0) pat%bgr(n) = values(pos(5))

                if (pat%sigma(n) > 0.001) then
                   chi2=chi2+((pat%y(n)-pat%ycalc(n))/pat%sigma(n))**2
                end if
             end do
             chi2=chi2/real(pat%npts)

             i=len_trim(pat%title)
             write(unit=pat%title(i+2:),fmt="(a,g12.4)") "  Chi2(free) = ",chi2

       end select

       return
    End Subroutine Read_Pattern_CIF

    !!--++
    !!--++ Subroutine Read_Pattern_Free
    !!--++
    !!--++    Read a pattern for Free
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Read_Pattern_Free(Filename,Pat,ext)
       !---- Arguments ----!
       character(len=*),           intent(in)  :: Filename      ! Path+Filename
       class(DiffPat_Type),        intent(out) :: Pat
       character(len=*), optional, intent(in)  :: ext

       !---- Local Variables ----!
       integer                                      :: i,no,ier,inum,nc,iv,i_dat
       integer, dimension(3)                        :: ivet
       character(len=180)                           :: aline
       character(len=20), dimension(10)             :: dire
       real(kind=cp), dimension(3)                  :: vet
       real(kind=cp)                                :: step
       logical                                      :: title_given,ext_given, rigaku,info


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
       ext_given=.false.
       rigaku=.false.

       no=0
       pat%ScatVar="2theta"

       if (present(ext)) ext_given=.true.
       if (ext_given) then
          if (trim(ext) == ".MDI") rigaku=.true.
       end if

       do
          read(unit=i_dat,fmt="(a)",iostat=ier) aline
          if (ier /= 0) then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML%Msg=" End of file *.dat"

             close(unit=i_dat)
             return
          end if

          aline=adjustl(aline)

          !> Comment lines using ! or #
          if (aline(1:1) == "!" .or. aline(1:1) == "#") then
             i=index(aline,"Legend_X")
             if (i /= 0) then
                !pat%xax_text=adjustl(aline(i+9:))
             end if

             i=index(aline,"Legend_Y")
             if (i /= 0) then
                !pat%yax_text=adjustl(aline(i+9:))
             end if

             i=index(aline,"Scattering variable:")
             if (i /= 0) then
                pat%ScatVar=adjustl(aline(i+20:))
             end if

             i=index(aline,"TSAMP")
             if (i /= 0) then
                !read(unit=aline(i+5:),fmt=*,iostat=ier) pat%tsamp
                !if (ier /= 0) pat%tsamp = 0.0
             end if

             i=index(aline,"TITLE")
             if (i /= 0) then
                Pat%title=trim(aline(i+6:))
                title_given=.true.
             end if

             i=index(aline,"Title:")
             if (i /= 0) then
                Pat%title=trim(aline(i+7:))
                title_given=.true.
             end if

             cycle
          end if

          !> BANK Information
          if (aline(1:4) == "BANK") then
             read(unit=aline(5:41),fmt=*) inum, pat%npts
             read(unit=aline(47:90),fmt=*) pat%xmin,step
             pat%xmax=pat%xmin+(pat%npts-1)*step
             exit
          end if

          if (rigaku) then
             if (.not. title_given) then
                title_given=.true.
                Pat%title=trim(adjustl(aline))
                cycle

             else
                call get_word(aline,dire,nc)
                call get_num(trim(dire(1))//' '//trim(dire(2))//' '//trim(dire(6)),vet,ivet,iv)
                pat%xmin=vet(1)
                step=vet(2)
                pat%xmax=vet(3)
                pat%npts = nint((pat%xmax-pat%xmin)/step+1.0)
                exit
             end if
          end if

          !> Reading Xmin, Step, Xmax, Title (optional)
          call get_word(aline,dire,nc)
          if (nc > 2) then
             call get_num(trim(dire(1))//' '//trim(dire(2))//' '//trim(dire(3)),vet,ivet,iv)
             if (iv == 3) then
                pat%xmin=vet(1)
                step=vet(2)
                pat%xmax=vet(3)

                if (step <= 1.0e-6 ) then
                   Err_CFML%state=.true.
                   Err_CFML%Flag=2
                   Err_CFML%Msg=" Error in Intensity file, Step value was zero!"

                   close(unit=-i_dat)
                   return
                end if

                pat%npts = nint((pat%xmax-pat%xmin)/step+1.0)

                !> Title?
                i=index(aline,trim(dire(3)))
                nc=len_trim(dire(3))

                if (len_trim(aline(i+nc+1:)) > 0 .and. .not. title_given) then
                   Pat%title=trim(aline(i+nc+1:))
                   title_given=.true.
                end if

                exit  ! Salida del Bucle
             end if

             !> TSAMP
             i=index(aline,"TSAMP")
             if (i /= 0) then
                !read(unit=aline(i+5:),fmt=*,iostat=ier) pat%tsamp
                !if (ier /= 0) pat%tsamp = 0.0
             end if
          end if

          !> Probably Coment line or Title
          if (.not. title_given) then
             Pat%title=trim(aline)
             title_given=.true.
          end if

          no=no+1
          if (no > 7)then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML%Msg=" Error on Intensity file, Number of Comment lines was exceeded ( > 7) !"

             close(unit=i_dat)
             return
          else
             cycle
          end if
       end do

       !> Aditional checks
       if (pat%npts <= 10 .or. pat%xmax <  pat%xmin  .or. step > pat%xmax) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, Problems reading 2Theta_ini, Step, 2Theta_end !"

          close(unit=i_dat)
          return
       end if

       ! Allocating memory
       call Allocate_Diffraction_Pattern(pat)

       ! Reading intensities values
       read(unit=i_dat,fmt=*,iostat=ier)(pat%y(i),i=1,pat%npts)
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity file, Number of intensities values is wrong!!"

          close(unit=i_dat)
          return
       end if

       do i=1,pat%npts
          pat%sigma(i) = pat%y(i)
          pat%x(i)= pat%xmin+(i-1)*step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       if (pat%scatvar == "2theta" .and. pat%xmax > 180.0 ) pat%Scatvar="TOF"

       close(unit=i_dat)
       return
    End Subroutine Read_Pattern_Free

    !!--++
    !!--++ Subroutine Read_Pattern_Time_Variable
    !!--++
    !!--++    Read a pattern for Time Variable
    !!--++
    !!--++ Update: January - 2005
    !!
    Module Subroutine Read_Pattern_Time_Variable(Filename, Pat)
       !---- Arguments ----!
       character(len=*),    intent(in)  :: Filename      ! Path+Filename
       class(DiffPat_Type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=180)                           :: txt1
       character(len=132)                           :: txt2
       character(len=132)                           :: txt3
       real(kind=cp), dimension(:), allocatable     :: bk
       real(kind=cp)                                :: cnorma,step
       integer                                      :: i,ier,i_dat
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
       read(unit=i_dat,fmt="(A)",iostat=ier) txt1   !1
       pat%title=trim(txt1)
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"

          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier) txt2   !2
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"

          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier) txt3   !3
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"

          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier) txt3   !4
       if (ier /= 0)then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"

          close(unit=i_dat)
          return
       end if

       read(unit=i_dat,fmt=*,iostat=ier)pat%xmin, step, pat%xmax
       if (ier /= 0)then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"

          close(unit=i_dat)
          return
       end if

       pat%npts = (pat%xmax-pat%xmin)/step+1.5
       if (pat%npts <=0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"

          close(unit=i_dat)
          return
       end if

       !> Allocating
       call Allocate_Diffraction_Pattern(pat)

       if(allocated(bk) ) deallocate(bk)
       allocate(bk(pat%npts))

       read(unit=i_dat,fmt="(5(F6.0,F10.0))",iostat=ier)(bk(i),pat%y(i),i=1,pat%npts)
       if (ier /= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in  Intensity file, check your instr parameter!"

          close(unit=i_dat)
          return
       end if

       !> Normalize data to constant time
       cnorma=0.0
       DO i=1,pat%npts
          IF (bk(i) < 1.0E-06) then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML%Msg=" Zero time in *.DAT "

             close(unit=i_dat)
             return
          end if
          cnorma=cnorma+bk(i)
          pat%x(i)= pat%xmin+(i-1)*step
       end do
       cnorma=cnorma/real(pat%npts)
       do i=1,pat%npts
          pat%y(i)=pat%y(i)*cnorma/bk(i)
          pat%sigma(i)=pat%y(i)
          bk(i)=0.0
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       close(unit=i_dat)

       return
    End subroutine Read_Pattern_Time_Variable

    !!--++
    !!--++ Subroutine Read_Pattern_XYSigma
    !!--++
    !!--++    Read a pattern for X,Y,Sigma.
    !!--++    Adding (2014) the possibility to read a calculated pattern
    !!--++    in a fouth column. If gr is present a PDFGUI pattern is read.
    !!--++    If header is present the full header of the file is stored in
    !!--++    the hopefully long string: header
    !!--++
    !!--++ Updated: January - 2014, Nov-2015 (JRC)
    !!
    Module Subroutine Read_Pattern_XYSigma(Filename, Pat, PDF, Header)
       !---- Arguments ----!
       character(len=*),           intent(in)  :: Filename      ! Path+Filename
       class(DiffPat_Type),        intent(out) :: Pat
       logical,          optional, intent(in)  :: PDF
       character(len=*), optional, intent (out):: Header

       !---- Local Variables ----!
       real(kind=cp), parameter                     :: EPS1=1.0E-6

       character(len=180)                           :: txt1, aline, fmtfields, fmtformat
       character (len=5)                            :: date1
       integer                                      :: line_da, ntt, interpol, i, j,ier,npp,lenhead,i_dat
       real(kind=cp)                                :: fac_x, fac_y,  yp1, sumavar, cnorm, step
       real(kind=cp)                                :: ycor, xt, stepin, ypn, dum
       real(kind=cp), dimension(:), allocatable     :: yc, bk
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

       !---- Or X,Y sigma data ----!
       fac_x=1.0
       fac_y=1.0
       interpol=0
       line_da=1
       npp=0
       ntt=0
       if (present(header)) then
          lenhead=len(header)
          header=" "
       end if

       do
          read(unit=i_dat,fmt="(a)", iostat=ier) txt1
          if (ier /= 0)   exit
          npp=npp+1
          if (index(txt1,"#L r(A)") /= 0) then
             line_da=npp
             npp=0
          end if
       end do

       pat%npts=npp
       if (pat%npts <= 0) then
          Err_CFML%state=.true.
          Err_CFML%Flag=2
          Err_CFML%Msg=" Error in Intensity/Profile file, Number of points negative or zero!"

          close(unit=i_dat)
          return
       end if
       rewind(unit=i_dat)

       if (present(PDF)) then
          do i=1, line_da
             read(unit=i_dat,fmt="(a)") txt1
             if (present(header)) header=trim(header)//trim(txt1)//char(0)
             j= index(txt1,"title=")
             if ( j /= 0 ) then !Title given
                pat%title=adjustl(txt1(j+6:))
             end if
          end do
          !pat%xax_text="Distance R(Angstroms)"
          !pat%yax_text="G(R)x1000"
          pat%scatvar="Distance"

       else
          read(unit=i_dat,fmt="(a)") txt1
          IF (txt1(1:6) /= "XYDATA") THEN
             pat%title=trim(txt1)
             do
                read(unit=i_dat,fmt="(a)", iostat=ier) txt1
                if (ier /= 0) then
                   Err_CFML%state=.true.
                   Err_CFML%Flag=2
                   Err_CFML%Msg=" Error reading a profile DATA file of XYSigma format"

                   close(unit=i_dat)
                   return
                end if
                txt1=adjustl(txt1)
                if(present(header)) header=trim(header)//trim(txt1)//char(0)
                j= index(txt1,"TITLE")
                if ( j /= 0 ) then !Title given
                   pat%title=adjustl(txt1(j+5:))
                   cycle
                end if

                i=index(txt1,"Scattering variable:")
                if (i /= 0) then
                   pat%scatvar=adjustl(txt1(i+20:))
                   cycle
                end if

                i=index(txt1,"Legend_X")
                if (i /= 0) then
                   !pat%xax_text = adjustl(txt1(i+9:))
                   cycle
                end if

                i=index(txt1,"Legend_Y")
                if (i /= 0) then
                   !pat%yax_text=adjustl(txt1(i+9:))
                   cycle
                end if

                if (txt1(1:1) == "!" .or. txt1(1:1) == "#") cycle
                read(unit=txt1, fmt=*, iostat=ier) yp1, ypn !This is to detect the beginning of numerical values
                if ( ier /= 0) cycle
                backspace (unit=i_dat)

                call init_findfmt(line_da)
                exit
             end do
             if (present(header)) then
                j=len_trim(header)-1
                i=index(header(1:j),char(0),back=.true.)
                header=header(1:i)
             end if
          else
             pat%title=trim(txt1)
             i=0
             do
                i=i+1
                if (i == 6) exit
                read(unit=i_dat,fmt="(a)", iostat=ier) txt1
                if (ier /= 0) then
                   Err_CFML%state=.true.
                   Err_CFML%Flag=2
                   Err_CFML%Msg=" Error reading a profile DATA file of XYSigma format"

                   close(unit=i_dat)
                   return
                end if
                if (txt1(1:1) == "#") i=i-1
                if (present(header)) header=trim(header)//trim(txt1)//char(0)

                line_da=line_da+1
                j= index(txt1,"TITLE")
                if ( j /= 0 ) then !Title given
                   pat%title=adjustl(txt1(j+5:))
                end if

                j=index(txt1,"Scattering variable:")
                if (j /= 0) then
                   pat%scatvar=adjustl(txt1(j+20:))
                end if

                j=index(txt1,"Legend_X")
                if (j /= 0) then
                   !pat%xax_text=adjustl(txt1(j+9:))
                end if

                j=index(txt1,"Legend_Y")
                if (j /= 0) then
                   !pat%yax_text=adjustl(txt1(j+9:))
                end if

                if (txt1(1:5) == "INTER") then !Interpolation possible!
                   backspace (unit=i_dat)
                   line_da=line_da-2
                   call init_findfmt(line_da)
                   fmtfields = "5ffif"
                   call findfmt(i_dat,aline,fmtfields,fmtformat)
                   if (ierr_fmt /= 0) then
                      Err_CFML%state=.true.
                      Err_CFML%Flag=2
                      Err_CFML%Msg=" Error reading"

                      close(unit=i_dat)
                      return
                   end if

                   read(unit=aline,fmt=fmtformat) date1,fac_x,fac_y,interpol,stepin
                   if (fac_x <= 0.0) fac_x=1.0
                   if (fac_y <= 0.0) fac_y=1.0
                end if

                if (txt1(1:4) == "TEMP") then
                   !read(unit=txt1(5:80),fmt=*, iostat=ier) pat%tsamp
                   !if(ier == 0) then
                   !  pat%tset=pat%tsamp
                   !else
                   !  pat%tsamp=0.0
                   !  pat%tset=0.0
                   !end if
                end if
             end do
          end if

          if (interpol == 0) then
             !pat%ct_step = .false.
          else if(interpol == 1) then
             !pat%ct_step = .true.
          else if(interpol == 2) then
             !pat%ct_step = .true.
          end if
       end if

       !> Allocating
       call Allocate_Diffraction_Pattern(pat)

       fmtfields = "ffff"  !Now four columns are read in order to incorporate the calcualted pattern
       sumavar=0.0
       cnorm=0.0
       i=0
       do j=1,pat%npts
          call findfmt(i_dat,aline,fmtfields,fmtformat)
          if (ierr_fmt == -1) exit
          if (ierr_fmt /= 0) then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML%Msg=" Error reading X,Y, Sigma Ycalc in profile DATA file"

             close(unit=i_dat)
             return
          end if
          if(aline(1:1) == "!" .or. aline(1:1) == "#") cycle

          i=i+1
          if (present(PDF)) then
             read(unit=aline,fmt = fmtformat, iostat=ier ) pat%x(i),pat%y(i),dum,pat%sigma(i)
          else
             !read(unit=aline,fmt = fmtformat, iostat=ier ) pat%x(i),pat%y(i),pat%sigma(i),pat%ycalc(i)
          end if
          if (ier /=0) then
             Err_CFML%state=.true.
             Err_CFML%Flag=2
             Err_CFML%Msg=" Error in Intensity file, check your instr parameter!"

             close(unit=i_dat)
             return
          end if
          if (i > 10 .and. ABS(pat%x(i)) < EPS1 .AND. pat%y(i) < EPS1 .AND.  pat%sigma(i) < EPS1) exit

          pat%x(i)=pat%x(i)*fac_x
          pat%y(i)=pat%y(i)*fac_y

          !pat%ycalc(i)=pat%ycalc(i)*fac_y

          pat%sigma(i)=pat%sigma(i)*fac_y
          pat%sigma(i)=pat%sigma(i)*pat%sigma(i)
          sumavar=sumavar+pat%sigma(i)
          if (pat%sigma(i) < EPS1) pat%sigma(i) =1.0_cp
          cnorm=cnorm+pat%sigma(i)/MAX(abs(pat%y(i)),0.001_cp)
          !if (i > 1) then
          !   step=step+pat%x(i)-pat%x(i-1)
          !end if
       end do

       ntt=i-1
       pat%xmin=pat%x(1)
       pat%xmax=pat%x(ntt)
       cnorm=cnorm/real(ntt)
       if (sumavar < EPS1) then
          do i=1,ntt
             pat%sigma(i)=abs(pat%y(i))
          end do
          cnorm=1.0
       end if

       if (interpol == 0 .or. interpol == 2) then      !if interpol
          step=pat%x(2)-pat%x(1)
          pat%npts=ntt

       else                        !else interpol
          step=stepin
          j=(pat%x(ntt)-pat%x(1))/step+1.05
          if( j > pat%npts) then
             step=(pat%x(ntt)-pat%x(1))/(ntt-1)
             pat%npts=ntt
          else
             pat%npts=j
          end if
          if(allocated(bk) ) deallocate(bk)
          allocate(bk(pat%npts))
          if(allocated(yc) ) deallocate(yc)
          allocate(yc(pat%npts))

          !>First interpolate the raw intensities ----!
          yp1=huge(1.0)
          ypn=huge(1.0)
          call spline(pat%x(:),pat%y(:),ntt,yp1,ypn,bk(:))
          do i=1,pat%npts
             xt=pat%x(1)+(i-1)*step
             ycor=splint(pat%x(:),pat%y(:),bk(:),ntt,xt)
             yc(i)=ycor !max(1.0_cp,ycor)
          end do
          do i=1,pat%npts
             pat%y(i)=yc(i)
             yc(i)=0.0
             bk(i)=0.0
          end do

          !> Second interpolate the sigmas
          call spline(pat%x(:),pat%sigma(:),ntt,yp1,ypn,bk(:))
          do i=1,pat%npts
             xt=pat%x(1)+(i-1)*step
             ycor=splint(pat%x(:),pat%sigma(:),bk(:),ntt,xt)
             yc(i)=ycor !max(1.0_cp,ycor)
          end do
          do i=1,pat%npts
             pat%sigma(i)=abs(yc(i))
             yc(i)=0.0
             bk(i)=0.0
          end do
          pat%xmax=pat%xmin+step*(pat%npts-1)
       end if                       !End If interpol
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       close(unit=i_dat)

       return
    End Subroutine Read_Pattern_XYSigma

End Submodule GenPatterns