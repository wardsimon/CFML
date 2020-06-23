Submodule (CFML_DiffPatt) Bckg_Patt

 implicit none

 real(kind=cp), parameter :: eps=0.0001_cp

 Contains

    !!----
    !!---- CALC_BACKGROUND
    !!----
    !!----    Calculate a Background using an iterative process according
    !!----    to Bruckner, S. (2000). J. Appl. Cryst., 33, 977-979.
    !!----
    !!----    Np is the extension of Np points at Left and Right. Normally it could be
    !!----    around 10-40 points.
    !!----
    !!---- 30/04/2019
    !!
    Module Subroutine Calc_BackGround(Pat, Ncyc, Np, Xmin, Xmax)
       !---- Arguments ----!
       class(DiffPat_E_Type),     intent(in out) :: Pat        ! Pattern object
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
       pat%bgr=0.0_cp

       !> Check
       call clear_error()
       if (pat%npts <= 1) then
          err_CFML%Ierr=1
          err_CFML%Msg="Calc_BackGround@DIFFPAT: There aren't points in the Pattern used for Background calculation!"
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
          err_CFML%Ierr=1
          err_CFML%Msg="Calc_BackGround@DIFFPAT: There aren't background points into the range"
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
          err_CFML%IErr=1
          err_CFML%Msg="Calc_BackGround@DIFFPAT: Same value for Xmin and Xmax! "
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

    End Subroutine Calc_BackGround

    !!----
    !!---- READ_BACKGOUND_FILE
    !!----
    !!----    Read background pattern from an external file.
    !!----    Mode:
    !!----         Poly | Inter
    !!----
    !!----
    !!---- 30/04/2019
    !!
    Module Subroutine Read_Background_File(Bck_File, Bck_Mode, Pat)
       !---- Arguments ----!
       character(len=*),         intent(in   )    :: bck_file      ! Path+Filename of Background file
       character(len=*),         intent(in   )    :: bck_mode
       class(DiffPat_E_Type),    intent(inout)    :: Pat

       !---- local variables ----!
       logical                                       :: esta
       character(len=80)                             :: line
       character(len=3)                              :: car
       integer                                       :: bck_points
       integer                                       :: i,j,i_bck
       integer                                       :: ier
       real(kind=cp), dimension (:), allocatable     :: bck_v
       real(kind=cp), dimension (:), allocatable     :: bck_p

       !> Init
       call clear_error()

       inquire(file=trim(bck_file), exist =esta)
       if (.not. esta) then
          Err_CFML%IErr=1
          Err_CFML%Msg="Read_Background_File@DIFFPATT: The file "//trim(bck_file)//" doesn't exist"
          return
       end if

       !> Open Backgriund file
       open(newunit=i_bck,file=trim(bck_file),status="old",action="read",position="rewind",iostat=ier)
       if (ier /= 0) then
          Err_CFML%Ierr=1
          Err_CFML%Msg="Read_Background_File@DIFFPATT: Problems opening the file: "//trim(bck_file)
          return
       end if

       !> Number of background points
       i=0
       do
          read(unit=i_bck,fmt="(a)",iostat=ier) line
          if (ier /= 0) exit

          if (len_trim(line) == 0)  cycle
          if (index(line,"!") /= 0) cycle
          if (index(line,"#") /= 0) cycle
          i=i+1
       end do
       bck_points=i

       if (bck_points <=0) then
          Err_CFML%IErr=1
          Err_CFML%Msg="Read_Background_File@DIFFPATT: Impossible to read any background point in the file: "//trim(bck_file)
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
       j=0
       do
          read(unit=i_bck,fmt="(a)",iostat=ier) line
          if (ier /= 0) exit

          if (len_trim(line) == 0 .or. line(1:1) == "!" .or. line(1:1)=="#") cycle

          j=j+1
          read(unit=line, fmt=*, iostat=ier)  bck_p(j), bck_v(j)
          if (ier /= 0) then
             bck_points=j-1
             close(unit=i_bck)

             Err_CFML%Ierr=1
             ERR_CFML%Msg="Read_Background_File@DIFFPATT: Problems during loading the background points!"
             exit
          end if
       end do
       close (unit=i_bck)

       bck_points=j
       if (bck_points <=0) then
          Err_CFML%IErr=1
          Err_CFML%Msg="Read_Background_File@DIFFPATT: Zero background points in the file: "//trim(bck_file)
          return
       end if

       car=u_case(adjustl(bck_mode))
       select case (car)
          case ("POL") ! Polynomial
             call set_background_poly (Pat, 50.0_cp, bck_p, bck_points )

          case ("INT") ! Interpolation
             call  set_background_inter (Pat, bck_v, bck_p, bck_points )

          case default
             Err_CFML%IErr=1
             ERR_CFML%Msg="Read_Background_File@DIFFPATT: Define the mode: Polynomial or Interpolation"
             return
       end select
    End Subroutine Read_Background_File

    !!--++
    !!--++ SET_BACKGROUND_POLY
    !!--++
    !!--++    (PRIVATE)
    !!--++    Define a n-polynomial with constanat vlue at bkpos Background
    !!--++
    !!--++ 30/04/2019
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
       pat%al_bgr=.true.

       do i=1, pat%npts
          do j=1,n
             pat%bgr(i)= pat%bgr(i) + bckx(j)*((pat%x(i)/bkpos-1.0)**(j-1))
          end do
       end do
    End Subroutine Set_Background_Poly

    !!--++
    !!--++ SET_BACKGROUND_INTER
    !!--++
    !!--++    (PRIVATE)
    !!--++    Define a Background
    !!--++
    !!--++ 30/04/2019
    !!
    Module Subroutine Set_Background_Inter(Pat, Bcky, Bckx, N)
       !---- Arguments ----!
       class(DiffPat_E_Type),         intent(in out) :: Pat
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
    End Subroutine Set_Background_Inter

End Submodule Bckg_Patt