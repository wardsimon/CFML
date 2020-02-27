  Module DataRed_Reflections

    Use CFML_GlobalDeps
    Use CFML_Rational
    Use CFML_Maths,   only: zbelong,epss,sort,set_eps_math
    Use CFML_Metrics, only: Cell_Type, Cell_G_Type
    Use CFML_Strings, only: number_lines
    Use CFML_gSpaceGroups
    Use CFML_Reflections
    Use CFML_Propagation_Vectors
    Use Datared_Mod
    Use Twin_Mod

    Implicit None
    Private

    Public :: Read_Reflections_File, Write_Reflections, Treat_Reflections

    integer, parameter :: nref=400000

    character(len=20)                    :: line1,wmess
    integer, dimension(nref)             :: idomain
    !integer, dimension(nref)             :: itreat, iord, nequv, ini, fin, numor, ivk, icod, warn
    integer, dimension(nmax)             :: contr

    real(kind=cp),    dimension(nref)    :: intav, sigmav, sigstat,lambda_laue, &
                                            tbar,twotet,absorpt
    real(kind=cp),    dimension(4)       :: angles
    real(kind=cp),    dimension(256)     :: weight  ! A maximum of 256 reflections equivalent to one of them can be treated
    real(kind=cp),    dimension(3,nmax)  :: hkln
    character(len=*),parameter,dimension(0:1) :: warn_mess=["                      ",  &
                                                             " <- Dubious reflection"]
    real    :: total, sig, suma, suman, sumaw, sumanw, Rint, Rwint, aver_sig, &
               wavel,sigg, int_rej, epsg, delt, warning, scal_fact, q2, aver_int
    integer :: i,j,k, ier, nr=0, iv, ns, rej, ival,  nkv, nn, &
               lenf, nin, hkl_type, ivp, nk, nequiv, L, drej, Lmin,i_refout, a,b

    Type, extends(Refl_Type)    :: ObsRef
      real(kind=cp),dimension(3):: hr      = 0.0  ! Real indices in reciprocal space
      real(kind=cp)             :: intens  = 0.0  ! Observed intensity
      real(kind=cp)             :: sigma   = 0.0  ! Estimated standard deviation
      real(kind=cp)             :: twtheta = 0.0  ! Scattering angle
      real(kind=cp)             :: omega   = 0.0  ! Angle of orienting device
      real(kind=cp)             :: chi     = 0.0  ! Angle of orienting device
      real(kind=cp)             :: phi     = 0.0  ! Angle of orienting device
      real(kind=cp)             :: tbar    = 0.0  ! Weigthed path for absortion
      real(kind=cp)             :: absorpt = 0.0  ! Transmission factor
      real(kind=cp)             :: lambda_laue = 0.0 !lambda of reflection in a Laue or TOF experiment
      integer                   :: idomain = 1    ! Indicator of the domain to which reflection refers
      integer                   :: icod    = 0    ! Code for treating the reflection
      integer                   :: pfn     = 0    ! Indicator of problem with the reflection
      integer                   :: numor   = 0    ! Number order in the data collection
    End Type ObsRef


    contains

    Subroutine Read_Reflections_File(filhkl,cond,kinfo,Cell,Ref,Gk)
      character(len=*),                        intent(in)     :: filhkl
      type(Conditions_Type),                   intent(in out) :: cond
      type(kvect_Info_Type),                   intent(in)     :: kinfo
      type(Cell_G_type),                       intent(in)     :: Cell
      class(RefList_Type),                     intent(out)    :: Ref
      type(Group_k_Type),optional,dimension(:),intent(in out) :: Gk
      !local variables
      integer                                :: i,nr,inp,ier,nkv,ivk,numor
      real(kind=cp)                          :: wavel
      character(len=132)                     :: line
      character(len=6)                       :: keyw
      integer,             dimension(3)      :: hkl
      real(kind=cp),       dimension(3)      :: h1,h2,h3
      type(ObsRef), dimension(:),allocatable :: Ob
      integer,      dimension(:),allocatable :: ind,pnf !pnf points to not found integer indices
      logical :: found

      nr=Number_lines(trim(filhkl))+1
      allocate(Ob(nr))
      do i=1,nr
         allocate(Ob(i)%h(3+kinfo%nk))
      end do

      Open(newunit=inp,file=trim(filhkl),status="old",action="read",position="rewind")

      call clear_error()

      nr=0

      Select Case(cond%hkl_type)

        Case(0)       !Shelx-like input file (3i4,2f8.2)

            if(.not. cond%cell_given) then
              Err_CFML%Msg=" => UNIT CELL not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if
            if(.not. cond%wave_given) then
              Err_CFML%Msg=" => WAVELENGTH not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if
            cond%forma="(3i4,2f8.2)"
            do
               nr=nr+1
               read(unit=inp,fmt="(3i4,2f8.2)",iostat=ier)  Ob(nr)%h, Ob(nr)%intens, Ob(nr)%sigma
               if(ier /= 0 .or. (Ob(nr)%h(1)==0 .and. Ob(nr)%h(2)==0 .and. Ob(nr)%h(3) == 0) ) then
                nr=nr-1
                exit
               end if
               Ob(nr)%s=h_s(Ob(nr)%h,cell)
               h1(:)=real(Ob(nr)%h(:))
               if(cond%transf_ind) then
                 h2=matmul(cond%transhkl,h1)
               else
                 h2=h1
               end if
               Ob(nr)%h=nint(h2)
               Ob(nr)%hr=Ob(nr)%h
            end do

        Case(10)       !Shelx-like HKLF5 input file (3i4,2f8.2,i4)

             if(.not. cond%cell_given) then
              Err_CFML%Msg=" => UNIT CELL not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if
            if(.not. cond%wave_given) then
              Err_CFML%Msg=" => WAVELENGTH not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if

            cond%forma="(3i4,2f8.2,i4)"
            do
               nr=nr+1
               read(unit=inp,fmt="(3i4,2f8.2,i4)",iostat=ier)  Ob(nr)%h, Ob(nr)%intens, Ob(nr)%sigma, Ob(nr)%idomain
               if(ier /= 0 .or. (hkl(1)==0 .and. hkl(2)==0 .and. hkl(3) == 0) ) then
                nr=nr-1
                exit
               end if
               h1(:)=real(Ob(nr)%h(:))
               if(cond%transf_ind) then
                 h2=matmul(cond%transhkl,h1)
               else
                 h2=h1
               end if
               Ob(nr)%h=nint(h2)
               Ob(nr)%hr=Ob(nr)%h
            end do

        Case(1)           !Free format  h,k,l,int,sigma (h,k,l real)

            if(.not. cond%cell_given) then
              Err_CFML%Msg=" => UNIT CELL not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if
            if(.not. cond%wave_given) then
              Err_CFML%Msg=" => WAVELENGTH not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if

            cond%forma="*"
            do
               nr=nr+1
               read(unit=inp,fmt=*,iostat=ier) h1, Ob(nr)%intens, Ob(nr)%sigma
               if(ier /= 0) then
                 nr=nr-1
                 exit
               end if
               if(.not. zbelong(h1(:))) then
                  if(cond%prop) then
                    cond%hkl_real=.true.
                  else
                    if(cond%transf_ind) then
                       h1=matmul(cond%transhkl,h1)
                       if(.not. zbelong(h1(:))) then
                         write(unit=*,fmt="(a,3f8.3)")" => Removed reflection (non-integer): ",h1(:)
                         nr= nr-1
                         cycle
                       end if
                    else
                       write(unit=*,fmt="(a,3f8.3)")" => Removed reflection (non-integer): ",h1(:)
                       nr= nr-1
                       cycle
                    end if
                  end if
               else
                  if(cond%transf_ind) then
                     h1=matmul(cond%transhkl,h1)
                  else
                     h1=real(nint(h1))
                  end if
               end if
               Ob(nr)%hr=h1
            end do

        Case(2)     ! Jane's format for .fsq:   Two lines per reflection

            if(.not. cond%cell_given) then
              Err_CFML%Msg=" => UNIT CELL not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if
            if(.not. cond%wave_given) then
              Err_CFML%Msg=" => WAVELENGTH not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if

            line="(i6,3f7.3,2f10.2,3f8.2)"
            cond%forma="(i6,3f7.3,2f10.2,3f8.2)"
            do
               nr=nr+1
               read(unit=inp,fmt=line,iostat=ier) &
                   Ob(nr)%numor,h1, Ob(nr)%intens, Ob(nr)%sigma, Ob(nr)%omega, Ob(nr)%chi, Ob(nr)%phi
               if(ier /= 0) then
                 nr=nr-1
                 exit
               end if
               if(.not. zbelong(h1)) then
                 if(cond%prop) then
                   cond%hkl_real=.true.
                 else
                   if(cond%transf_ind) then
                      h1=matmul(cond%transhkl,h1)
                   else
                      write(unit=*,fmt="(a,3f8.3)")" => Removed reflection (non-integer): ",h1
                      nr= nr-1
                      cycle
                   end if
                 end if
               else
                   if(cond%transf_ind) then
                      h1=matmul(cond%transhkl,h1)
                   else
                      h1=real(nint(h1))
                   end if
               end if
               Ob(nr)%hr=h1
               read(unit=inp,fmt="(a)",iostat=ier) keyw
               if(ier /= 0) exit
            end do

        Case(3)                     !FullProf format

            read(unit=inp,fmt="(a)") line
            read(unit=inp,fmt="(a)") line
            read(unit=inp,fmt=*) wavel

            cond%forma=line

            do
               nr=nr+1
               read(unit=inp,fmt=line,iostat=ier) &
                    hkl(:), Ob(nr)%intens, Ob(nr)%sigma,i,Ob(nr)%twtheta, Ob(nr)%omega, Ob(nr)%chi, Ob(nr)%phi
               if(ier /= 0) then
                 nr=nr-1
                 exit
               end if
               h1=hkl
               if(cond%transf_ind) then
                  h1=matmul(cond%transhkl,h1)
               end if
               Ob(nr)%hr=h1
               Ob(nr)%h =hkl
            end do

        Case(4)                  !COLL 5 hkl integer

            cond%forma="(i6,3i4,2f10.2,4f8.2)"
            do
               nr=nr+1
               read(unit=inp,fmt="(i6,3i4,2f10.2,4f8.2)",iostat=ier) &
                    Ob(nr)%numor, Ob(nr)%h, Ob(nr)%intens, Ob(nr)%sigma,Ob(nr)%twtheta, Ob(nr)%omega, Ob(nr)%chi, Ob(nr)%phi
               if(ier /= 0 .or. Ob(nr)%numor == 0) then
                 nr=nr-1
                 exit
               end if
               Ob(nr)%hr=Ob(nr)%h
            end do

        Case(5) !COLL 5 hkl real numbers (old and new format with autodetection)

            read(unit=inp,fmt="(a)" ,iostat=ier) line  !Detect one of the two formats

            if(line(10:10) == ".") then  !If this position is . the old format is active
               if(line(16:16) == ".") then
                 cond%forma="(i6,3f6.2,f8.0,f4.0,4f8.2)"
               else
                 cond%forma="(i6,3f7.3,2f10.2,4f8.2)"
               end if
            else ! autodetection of the new format
               a = index(line, ".")
               b = index(line(a+1:), ".")
               if (b < 10) then
                 write(cond%forma,"(a,i1,a,i1,a)") '(i6,3f', b, '.', b+6-a, ',2f10.2,4f8.2)'
               else
                 write(cond%forma,"(a,i2,a,i1,a)") '(i6,3f', b, '.', b+6-a, ',2f10.2,4f8.2)'
               end if
               write(unit=*,fmt=*)"=> Will use format ",cond%forma
            end if
            rewind(unit=inp)

            do
               nr=nr+1
               read(unit=inp,fmt=cond%forma ,iostat=ier) &
                    Ob(nr)%numor,  Ob(nr)%hr, Ob(nr)%intens, Ob(nr)%sigma,Ob(nr)%twtheta, Ob(nr)%omega, Ob(nr)%chi, Ob(nr)%phi
               if(ier /= 0) then
                 nr=nr-1
                 exit
               end if
               h1=Ob(nr)%hr
               if(.not. zbelong(h1)) then
                 if(cond%prop) then
                   cond%hkl_real=.true.
                 else
                   if(cond%transf_ind) then
                      h1=matmul(cond%transhkl,h1)
                   else
                      write(unit=*,fmt="(a,3f8.3)")" => Removed reflection (non-integer): ",h1
                      nr= nr-1
                      cycle
                   end if
                 end if
               else
                 if(cond%transf_ind) then
                    h1=matmul(cond%transhkl,h1)
                 else
                    h1=real(nint(h1))
                 end if
               end if
               Ob(nr)%hr=h1
            end do

        Case(6)    !COLL 5 hkl real numbers - 6T2 ?  Only difference with 4: i4 instead of i6 for the numor
                   !and 3f6.2 instead of 3i4 for reading the hkl-indices
            cond%forma="(i4,3f6.2,2f10.2,4f8.2)"
            do
               nr=nr+1
               read(unit=inp,fmt="(i4,3f6.2,2f10.2,4f8.2)" ,iostat=ier) &
                    Ob(nr)%numor, h1, Ob(nr)%intens, Ob(nr)%sigma,Ob(nr)%twtheta, Ob(nr)%omega, Ob(nr)%chi, Ob(nr)%phi
               if(ier /= 0) then
                 nr=nr-1
                 exit
               end if
               if(.not. zbelong(h1)) then
                 if(cond%prop) then
                   cond%hkl_real=.true.
                 else
                   if(cond%transf_ind) then
                      h1=matmul(cond%transhkl,h1)
                   else
                      write(unit=*,fmt="(a,3f8.3)")" => Removed reflection (non-integer): ",h1
                      nr= nr-1
                      cycle
                   end if
                 end if
               else
                 if(cond%transf_ind) then
                    h1=matmul(cond%transhkl,h1)
                 else
                    h1=real(nint(h1))
                 end if
               end if
               Ob(nr)%hr=h1
            end do

        Case(7)    ! SXD data for FullProf

            if(.not. cond%cell_given) then
              Err_CFML%Msg=" => UNIT CELL not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if

            if(cond%prop) cond%hkl_real=.true.

            do
               read(unit=inp,fmt="(a)",iostat=ier) line
               if( ier /= 0) exit
               if(line(1:1) == "#" .or. line(1:1) == "!") cycle
               if(line(1:1) == "(") then
                cond%forma=trim(line)
                read(unit=inp,fmt="(a)",iostat=ier) line
                if(cond%prop) then
                    read(unit=inp,fmt="(a)",iostat=ier) line
                    read(unit=line,fmt=*,iostat=ier) nkv
                    do i=1,nkv
                       read(unit=inp,fmt="(a)",iostat=ier) line
                    end do
                end if
                exit
               end if
            end do

            do
               read(unit=inp,fmt="(a)",iostat=ier) line
               if( ier /= 0) exit
               if(line(1:1) == "#" .or. line(1:1) == "!")  cycle
               nr=nr+1
               if(cond%prop) then
                 read(unit=line,fmt=cond%forma,iostat=ier) hkl(:), ivk,Ob(nr)%intens, Ob(nr)%sigma, Ob(nr)%icod, Ob(nr)%Lambda_Laue, &
                                                      Ob(nr)%twtheta,Ob(nr)%absorpt,Ob(nr)%tbar
                  if(ier /= 0 .or. (hkl(1)==0 .and. hkl(2)==0 .and. hkl(3) == 0 .and. ivk == 0) ) then
                   nr=nr-1
                   exit
                  end if
               else
                 read(unit=line,fmt=cond%forma,iostat=ier) hkl(:), Ob(nr)%intens, Ob(nr)%sigma, Ob(nr)%icod, Ob(nr)%Lambda_Laue, &
                                                      Ob(nr)%twtheta,Ob(nr)%absorpt,Ob(nr)%tbar
                  if(ier /= 0 .or. (hkl(1)==0 .and. hkl(2)==0 .and. hkl(3) == 0) ) then
                   nr=nr-1
                   exit
                  end if
               end if
               h1(:)=real(hkl(:))
               if(cond%transf_ind) then
                 h2=matmul(cond%transhkl,h1)
                 if(cond%prop) then
                   h3=matmul(cond%transhkl,kinfo%kv(:,ivk))
                   h2=h2+h3
                 end if
               else
                 h2=h1
                 if(cond%prop) h2=h2+kinfo%kv(:,ivk)
               end if
               Ob(nr)%hr=h2
            end do

        Case(8)       ! Flipping ratios from D3

            if(.not. cond%cell_given) then
              Err_CFML%Msg=" => UNIT CELL not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if
            if(.not. cond%wave_given) then
              Err_CFML%Msg=" => WAVELENGTH not GIVEN! Modify your input file."
              Err_CFML%Ierr=1
              return
            end if

            cond%forma="(i8,6f8.3,2f10.2)"
            do
               nr=nr+1      ! gamma  omega  nu  => in fact it is: omega, gamma, nu!!!!
               read(unit=inp,fmt=cond%forma,iostat=ier) &
                   Ob(nr)%numor,h1, Ob(nr)%omega, Ob(nr)%phi,  Ob(nr)%chi,  Ob(nr)%intens, Ob(nr)%sigma
                   !                   ^omega     ^gamma    ^nu
               if(ier /= 0) then
                 nr=nr-1
                 exit
               end if
               if(.not. zbelong(h1)) then
                   if(cond%transf_ind) then
                      h1=matmul(cond%transhkl,h1)
                   else
                      write(unit=*,fmt="(a,3f8.3)")" => Removed reflection (non-integer): ",h1
                      nr= nr-1
                      cycle
                   end if
               else
                   if(cond%transf_ind) then
                      h1=matmul(cond%transhkl,h1)
                   else
                      h1=real(nint(h1))
                   end if
               end if
               Ob(nr)%hr=h1
            end do

        Case (9)  !i6,3f10.4,2f10.2,4f8.2  Suggested by Juerg Schefer
            cond%forma="(i6,3f10.4,2f10.2,4f8.2)"
            do
               nr=nr+1
               read(unit=inp,fmt="(i6,3f10.4,2f10.2,4f8.2)" ,iostat=ier) &
                    Ob(nr)%numor, h1, Ob(nr)%intens, Ob(nr)%sigma,Ob(nr)%twtheta, Ob(nr)%omega, Ob(nr)%chi, Ob(nr)%phi
               if(ier /= 0) then
                 nr=nr-1
                 exit
               end if

               if(.not. zbelong(h1)) then
                 if(cond%prop) then
                   cond%hkl_real=.true.
                 else
                   if(cond%transf_ind) then
                      h1=matmul(cond%transhkl,h1)
                   else
                      write(unit=*,fmt="(a,3f8.3)")" => Removed reflection (non-integer): ",h1
                      nr= nr-1
                      cycle
                   end if
                 end if
               else
                 if(cond%transf_ind) then
                    h1=matmul(cond%transhkl,h1)
                 else
                    h1=real(nint(h1))
                 end if
               end if
               Ob(nr)%hr=h1
            end do


      End Select

      do i=1,nr
         h1=Ob(i)%hr
         if(cond%prop) then
           Ob(i)%s=h_s(Ob(i)%hr,cell)
         else
           Ob(i)%h=nint(h1)
           Ob(i)%s=h_s(Ob(i)%h,cell)
         end if
         if(cond%hkl_type == 7) then
           Ob(i)%twtheta=2.0* asind(Ob(i)%s*Ob(i)%Lambda_Laue)
         else
           Ob(i)%twtheta=2.0* asind(Ob(i)%s*cond%wavel)
         end if
         if(Ob(i)%twtheta < 0.0001)  ier=1
         if(Ob(i)%intens  < 0.0001)  Ob(i)%intens=0.0001
         if(Ob(i)%sigma   < 0.00001) Ob(i)%sigma=sqrt(abs(Ob(i)%intens))
      end do

      Ref%nref=nr
      allocate(ind(nr))
      ind=sort(Ob(:)%s,nr)
      allocate(Ref%Ref(nr),source=Ob)
      do i=1,nr
        j=ind(i)
        !Now calculate the integer indices
        if(present(Gk)) then
          call get_integer_indices(Ob(j)%hr,kinfo,Ob(j)%h,found,Gk)
        else
          call get_integer_indices(Ob(j)%hr,kinfo,Ob(j)%h,found)
        end if
        if(.not. found) then
           Ob(j)%pfn=1
           write(unit=*,fmt="(a,3f8.4)") " => Warning: no combination of the given Q_coefficients have been found for reflection: ",Ob(j)%hr
        end if
        Ref%Ref(i)=Ob(j)
      end do

    End Subroutine Read_Reflections_File

    Subroutine get_integer_indices(hr,kinfo,h,found,Gk)
      real(kind=cp), dimension(3),intent(in) :: hr
      type(kvect_Info_Type),      intent(in) :: kinfo
      integer, dimension(:),      intent(out):: h
      logical,                    intent(out):: found
      type(Group_k_Type),dimension(:),optional,intent(in) :: Gk
      !
      integer                        :: i,j,k,L,m,nk,ia
      real(kind=cp)                  :: rv
      real(kind=cp), dimension(3)    :: ha,hc,kv
      integer,       dimension(3)    :: ih
      real(kind=cp), dimension(3,48,kinfo%nk) :: kvec

      ha=nint(hr)
      h=0; found=.false.;kvec=0.0
      if(sum(abs(ha-hr)) < epss) then  !Fundamental reflection
        h(1:3)=nint(hr)
        found=.true.
        return
      end if

      do i=1,3
        j=int(hr(i))
        rv=hr(i)-real(j,kind=cp)
        if(rv > 0.50001_cp) then
          ih(i)=nint(hr(i))
        else if(rv <= -0.50001_cp) then
          ih(i)=int(hr(i))
        else
          ih(i)=j
        end if
      end do

      ha=ih
      h(1:3)=ih
      nk=kinfo%nk

      if(present(Gk)) then
        m=0
        do_j:do j=1,kinfo%nk
          do k=1,Gk(j)%ngk
            do ia=-1,1,2
              m=m+1
              kvec(:,m,j)=ia*Gk(j)%Eqv_k(:,k)
              write(*,"(i8,3f10.4)") m, kvec(:,m,j)
            end do
          end do
        end do do_j

        do_i:do i=1,kinfo%nq
          do ia=-1,1,2
            h(4:3+nk)=ia*kinfo%q_coeff(1:nk,i)
            do L=1,m
              hc=0.0
              do j=1,kinfo%nk
                  hc=hc+h(3+j)*kvec(:,L,j)
              end do
              hc=ha+hc
              if(sum(abs(hc-hr)) < epss) then !The good q_coeff are found!
                found=.true.
                exit do_i
              end if
            end do
            write(*,"(2i3,tr5,3f8.4,tr5,3f8.4,tr5,10i4)") i , ia, hr, hc, h
          end do !ia
        end do do_i

      else

        do_ii:do i=1,kinfo%nq
          do ia=-1,1,2
              h(4:3+nk)=ia*kinfo%q_coeff(1:nk,i)
              hc=0.0
              do j=1,kinfo%nk
                hc=hc+h(3+j)*kinfo%kv(:,j)
              end do
              hc=ha+hc
              if(sum(abs(hc-hr)) < epss) then !The good q_coeff are found!
                found=.true.
                exit do_ii
              end if
          end do !ia
        end do do_ii

      end if

    End Subroutine get_integer_indices

    Subroutine Write_Reflections(R,kinfo,lun)
      class(RefList_Type),  intent(in) :: R
      type(kvect_Info_Type),intent(in) :: kinfo
      integer, optional,    intent(in) :: lun

      integer :: iou,i,j,n,mr    !123456789012345678901234567890123456789
      character(len=39) :: fm="(i7, i4,a,3f8.4,2f12.3,f14.6,3f8.4,2i8)"
      character(len=:),allocatable :: line
      real(kind=cp), dimension(3)  :: hr

      iou=6
      if(present(lun)) iou=lun
      Select Case(kinfo%nk)
        case(0)
          line="List of read reflections"
        case default
          line="List of read reflections with extended integer indices"
      End Select
      write(unit=iou,fmt="(/,t5,a)") line
      line=repeat("=",len_trim(line))
      write(unit=iou,fmt="(t5,a)") line
      Select Type( ob => R%Ref)
        Type is (ObsRef)
          n=size(ob(1)%h)
          write(unit=fm(5:5),fmt="(i1)") n
          mr=0
          Select Case(kinfo%nk)
            case(0)         !123456789012345678901234567890123456789
              line=" NumRef   h   k   l      Hr      Kr      Lr       Intens        Sigma   SinTheta/Lambda"
            case(1)
              line=" NumRef   h   k   l   m      Hr      Kr      Lr       Intens        Sigma   SinTheta/Lambda"
            case(2)
              line=" NumRef   h   k   l   m   n      Hr      Kr      Lr       Intens        Sigma   SinTheta/Lambda"
            case(3)
              line=" NumRef   h   k   l   m   n   p      Hr      Kr      Lr       Intens        Sigma   SinTheta/Lambda"
          End Select
          write(unit=iou,fmt="(a)") line
          do i=1,R%nref
            hr=Ob(i)%h(1:3)
            do j=1,kinfo%nk
              hr=hr+Ob(i)%h(3+j)*kinfo%kv(:,j)
            end do
            if(Ob(i)%pfn == 1) then
              mr=mr+1
              write(unit=iou,fmt=fm) i,Ob(i)%h,"  ",Ob(i)%hr,Ob(i)%intens,Ob(i)%sigma,Ob(i)%s,hr,Ob(i)%pfn,Ob(i)%numor
            else
              write(unit=iou,fmt=fm) i,Ob(i)%h,"  ",Ob(i)%hr,Ob(i)%intens,Ob(i)%sigma,Ob(i)%s
            end if
          end do
      End Select
    End Subroutine Write_Reflections

    Subroutine Treat_Reflections(R,cond,cell,SpG,kinfo,tw,lun)
      Class(RefList_Type),  intent(in out) :: R
      type(Conditions_Type),intent(in)     :: cond
      type(Cell_G_type),    intent(in)     :: Cell
      class(SpG_type),      intent(in)     :: SpG
      type(kvect_Info_Type),intent(in)     :: kinfo
      type(Twin_Type),      intent(in)     :: tw
      integer,optional,     intent(in)     :: lun

    End Subroutine Treat_Reflections

  End Module DataRed_Reflections
