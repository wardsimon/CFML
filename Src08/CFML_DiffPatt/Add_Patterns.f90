Submodule (CFML_DiffPatt) AddPatt

 Contains

    !!----
    !!---- ADD_PATTERNS
    !!----
    !!---- Add Patterns
    !!----
    !!---- 30/04/2019 
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
        if (all(active) .eqv. .false.) return

        !> Initial values
        xmin=minval(Patterns(1:N)%xmin, mask= (active .eqv. .true.) )
        xmax=maxval(Patterns(1:N)%xmax, mask= (active .eqv. .true.) )

        npts=maxval(Patterns(1:N)%npts, mask= (active .eqv. .true.) )
        if (npts <= 1) then
           Err_CFML%IErr=1
           Err_CFML%Msg="Add_Patterns@DIFFPAT: Number of Points in the new Pattern was zero! "
           return
        end if

        step=(xmax-xmin)/real(npts-1)
        if (abs(step) <= epsilon(1.0_cp)) then
           Err_CFML%IErr=1
           Err_CFML%Msg="Add_Patterns@DIFFPAT: Step size in the new Pattern was close to zero! "
           return
        end if

        !> Second Derivative
        if (allocated(d2y)) deallocate(d2y)
        allocate(d2y(npts,n))
        d2y=0.0_cp

        do i=1,n
           if (.not. active(i)) cycle
           d2y(:,i)=second_derivative(Patterns(i)%x,Patterns(i)%y,Patterns(i)%npts)
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
              y=spline_interpol(Patterns(i)%x,Patterns(i)%y,d2y(:,i),Patterns(i)%npts,Pat%x(j))
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
              Pat%y(i)=0.0_cp
              Pat%sigma(i)=1.0_cp
              !Pat%nd(i)=0
           end if
        end do

        !Pat%Monitor=cnorm
        Pat%xmin=xmin
        Pat%xmax=xmax
        !Pat%step=step
        Pat%ymin=minval(Pat%y)
        Pat%ymax=maxval(Pat%y)
    End Subroutine Add_Patterns

End Submodule AddPatt