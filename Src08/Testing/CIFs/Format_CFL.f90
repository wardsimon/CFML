!!----
!!----
!!----
SubModule (CFML_IOForm) IOF_CFL
   Contains

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Class (SPG_Type),             intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++    Type(Job_Info_type), optional,intent(out)  :: Job_Info
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CFL File
    !!--++
    !!--++ Update: February - 2020
    !!
    Module Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,Type_Atm,CFrame,NPhase,Job_Info,xyz_type)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       class (Cell_Type),            intent(out)  :: Cell
       class (SPG_Type),             intent(out)  :: SpG
       Type (AtList_Type),           intent(out)  :: A
       character(len=*),             intent(in)   :: Type_Atm
       character(len=*),    optional,intent(in)   :: CFrame
       Integer,             optional,intent(in)   :: Nphase
       Type(Job_Info_type), optional,intent(out)  :: Job_Info
       character(len=*),    optional,intent(in)   :: xyz_type

       !---- Local variables ----!
       character(len=132)               :: line
       integer                          :: i, nauas, ndata, iph, n_ini,n_end,k
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip

       real(kind=cp),dimension(3):: vet

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!
       iph=1
       if (present(nphase)) iph=nphase
       if (present(Job_Info)) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Get_Job_Info(file_dat,n_ini,n_end,Job_info)
       end if

       !---- Reading Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if(present(CFrame)) then
         call Read_CFL_Cell(file_dat,n_ini,n_end,Cell,CFrame) !Read and construct Cell
       else
         call Read_CFL_Cell(file_dat,n_ini,n_end,Cell) !Read and construct Cell
       end if
       if (err_CFML%Ierr /= 0) return

       !---- Reading Space Group Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       call Read_CFL_SpG(file_dat,n_ini,n_end,SpG)
       if (err_CFML%Ierr /= 0) return

       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       !---- Calculating number of Atoms in the Phase ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (l_case(line(1:4)) == "atom")  nauas=nauas+1
       end do

       if (nauas > 0) then
          call Allocate_atom_list(nauas,A,Type_Atm)  !allocation space for Atom list
          call Read_CFL_Atoms(file_dat,n_ini,n_end,A,Type_Atm)
          if (err_CFML%Ierr /= 0) return

          do i=1,A%natoms
             vet=A%atom(i)%x
             A%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
             if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)
          end do

          Select Type (Cell)
            Type is(Cell_G_Type)
               do i=1,A%natoms
                  if (A%atom(i)%thtype == "aniso") then
                     select case (A%atom(i)%Utype)
                        case ("u_ij")
                           A%atom(i)%u(1:6) =  Get_U_from_Betas(A%atom(i)%u(1:6),Cell)
                        case ("b_ij")
                           A%atom(i)%u(1:6) =  Get_Betas_from_U(A%atom(i)%u(1:6),Cell)
                     end select
                     A%atom(i)%Utype="beta"
                  end if
               end do
          End Select
       end if
    End Subroutine Readn_Set_XTal_CFL
   !!----
   !!---- READ_CFL_ATOM
   !!----    Subroutine to read atoms parameters
   !!----
   !!----    Format: ATOM   Label  ChemSymb   x y z B Occ Us or Moment
   !!----
   !!----     For charge obtained from Label: Label[+/-][number]
   !!----
   !!----
   !!---- 25/06/2019
   !!
   Module Subroutine Read_CFL_Atoms(lines,n_ini, n_end, At_List,Type_Atm)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in)     :: lines
      integer,                        intent(in out) :: n_ini
      integer,                        intent(in)     :: n_end
      Type (AtList_Type),             intent(out)    :: At_List
      character(len=*),               intent(in)     :: Type_Atm

      !---- Local variables -----!
      character(len=12), parameter      :: DIGPM="0123456789+-"
      character(len=:), allocatable     :: line
      character(len=180), dimension(1)  :: extline
      character(len=10)                 :: dire
      character(len=20)                 :: label
      character(len=2)                  :: Chemsymb
      integer                           :: i, n, na, npos, nlong, iv
      real(kind=cp), dimension (20)     :: vet1, vet2


      !> Init
      call clear_error()
      call Allocate_Atom_List(0, At_list, Type_Atm)

      !> Calculate number of Atoms
      na=0
      do i=n_ini, n_end
         !> Checks
         if (len_trim(lines(i)) <=0) cycle                 ! Empty line
         line=adjustl(lines(i))
         if (line(1:1) =='!' .or. line(1:1) =='#') cycle   ! Comment line
         npos=index(u_case(line),'ATOM')
         if (npos > 0) na=na+1
      end do
      if (na == 0) return             ! No atoms in the lines

      !> Allocate List
      call Allocate_Atom_List(na, At_list,Type_Atm)

      !> Read Information
      na=0

      do i=n_ini, n_end
         !> Checks
         if (len_trim(lines(i)) <=0) cycle                 ! Empty line
         line=adjustl(lines(i))
         if (line(1:1) =='!' .or. line(1:1) =='#') cycle   ! Comment line

         !> Truncate line from symbols: # and !
         npos=index(line,'!')                             ! Cut line according to !
         if (npos > 0) line=line(:npos-1)
         npos=index(line,'#')
         if (npos > 0) line=line(:npos-1)

         !> Tabs
         do
            npos=index(line,TAB)
            if (npos == 0) exit
            line(npos:npos)=' '
         end do

         !> ATOM Directive
         call cut_string(line,nlong,dire)
         if (u_case(trim(dire)) /= 'ATOM') cycle

         na=na+1

         !> Label
         call cut_string(line,nlong,label)
         At_list%atom(na)%Lab=trim(label)

         !> Charge / Oxidation state
         !> Anions
         npos=index(label,'-')
         if (npos > 0) then
            iv=index(DIGPM,label(npos+1:npos+1))
            if (iv == 0 .or. iv > 10) then
               err_CFML%Ierr=1
               err_CFML%Msg='Read_CFL_Atoms@CFML_IOFORM: Error reading Charge information: '//trim(line)
               At_list%Natoms=Max(0,na-1)
               return
            end if
            At_list%atom(na)%charge=-(iv-1)
         end if

         !> Cations
         npos=index(label,'+')
         if (npos > 0) then
            iv=index(DIGPM,label(npos+1:npos+1))
            if (iv == 0 .or. iv > 10) then
               err_CFML%Ierr=1
               err_CFML%Msg='Read_CFL_Atoms@CFML_IOFORM: Error reading Charge information: '//trim(line)
               At_list%Natoms=Max(0,na-1)
               return
            end if
            At_list%atom(na)%charge=iv-1
         end if


         !> Chemical symbol / Scattering Factor
         call cut_string(line,nlong,dire)
         npos=index(DIGPM,dire(2:2))
         if (npos > 0) then
            Chemsymb=u_case(dire(1:1))
         else
            Chemsymb=u_case(dire(1:1))//l_case(dire(2:2))
         end if
         At_list%atom(na)%ChemSymb=ChemSymb
         At_list%atom(na)%SFacSymb=dire

         !> Parameters
         extline='atm '//trim(line)
         n=1
         call Read_Key_ValueSTD(extline,n,n,"atm",vet1,vet2,iv)
         if (iv <= 0) then
            err_CFML%Ierr=1
            err_CFML%Msg='Read_CFL_Atoms@CFML_IOFORM: Error reading Atom information: '//trim(line)
            At_list%Natoms=Max(0,na-1)
            return
         end if

         if (iv < 3) then
            err_CFML%Ierr=1
            err_CFML%Msg='Read_CFL_Atoms@CFML_IOFORM: Error reading Atom coordinates: '//trim(line)
            At_list%Natoms=Max(0,na-1)
            return
         end if

         !> Coordinates
         At_List%atom(na)%x=vet1(1:3)
         select type (atm => at_list%atom)
            class is (Atm_Std_Type)
               atm(na)%x_std=vet2(1:3)
         end select

         !> Biso
         if (iv > 3) then
            At_List%atom(na)%U_iso=vet1(4)
            select type (atm => at_list%atom)
               class is (Atm_Std_Type)
                  atm(na)%u_iso_std=vet2(4)
            end select
         end if

         !> Occ
         if (iv > 4) then
            At_List%atom(na)%Occ=vet1(5)
            select type (atm => at_list%atom)
               class is (Atm_Std_Type)
                  atm(na)%Occ_std=vet2(5)
            end select
         end if

         select case (iv)
            case (8)  ! Only moments
               At_List%atom(na)%Moment=vet1(6:8)
               select type (atm => at_list%atom)
                  class is (Atm_Std_Type)
                     atm(na)%Moment_std=vet2(6:8)
               end select
               At_List%atom(na)%Magnetic=.true.

            case (11) ! Only anisotropic parameters
               At_List%atom(na)%U=vet1(6:11)
               select type (atm => at_list%atom)
                  class is (Atm_Std_Type)
                     atm(na)%U_std=vet2(6:11)
               end select

            case (14) ! Anisotropic + Moments
               At_List%atom(na)%U=vet1(6:11)
               At_List%atom(na)%Moment=vet1(12:14)
               select type (atm => at_list%atom)
                  class is (Atm_Std_Type)
                     atm(na)%U_std=vet2(6:11)
                     atm(na)%Moment_std=vet2(12:14)
               end select

            case default
               err_CFML%Ierr=1
               err_CFML%Msg='Read_CFL_Atoms@CFML_IOFORM: Error reading Atom information: '//trim(line)
               At_list%Natoms=Max(0,na-1)
               return
         end select
         !Now verify if magnetic moments or anisotropic thermal parameters are provided in the next line
         line=adjustl(lines(i+1))
         call cut_string(line,nlong,dire)
         if (u_case(trim(dire)) == 'MOMENT') then !

         else

         if (u_case(lines(i+1)(1:4)) /= 'ATOM')  cycle
         end if
      end do
   End Subroutine Read_CFL_Atoms

   !!----
   !!---- READ_CFL_CELL
   !!----    Obtaining Cell Parameter from CFL Format
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_CFL_Cell(lines, n_ini, n_end, Cell, CFrame)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      class(Cell_Type),                intent(out)    :: Cell    ! Cell object
      Character(len=*), optional,      intent( in)    :: CFrame
      !---- Local variables -----!
      integer                              :: nlong, iv
      real(kind=cp), dimension (6)         :: vcell, std

      character(len=132)                   :: line
      character(len=4)                     :: dire

      !> Init
      call clear_error()

      call Read_Key_ValueSTD(lines,n_ini,n_end,"CELL",vcell,std,iv)
      if (iv /= 6) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: Problems reading cell parameters!"
         return
      end if
      if(present(CFrame)) then
         call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, CarType=CFrame,Vscell=std(1:3), Vsang=std(4:6))
      else
         call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, Vscell=std(1:3), Vsang=std(4:6))
      end if
   End Subroutine Read_CFL_Cell

    !!----
    !!---- Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
    !!----   character(len=*), dimension(:), intent( in) :: file_dat     !Lines of text (content of a file)
    !!----   integer,                        intent( in) :: i_ini,i_end  !Lines to explore
    !!----   type(job_info_type),            intent(out) :: Job_info     !Object to be constructed here
    !!----
    !!----
    !!----    Constructor of the object Job_info. The arrary of strings file_dat
    !!----    have to be provided as input. It contains lines corresponding to the
    !!----    input control file. The analysis of the command lines is not given here.
    !!----
    !!---- Update: February - 2005, February -2020
    !!
    Module Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
       !---- Arguments ----!
       character(len=*), dimension(:), intent( in) :: file_dat
       integer,                        intent( in) :: i_ini,i_end
       type(job_info_type),            intent(out) :: Job_info

       !---- Local Variables ----!
       integer                           :: i,nphas, ncmd,n_pat,ier, j
       integer, dimension(i_end-i_ini+1) :: ip,ic,ipt
       real(kind=cp)                     :: a1,a2,a3,a4,a5
       character(len=120)                :: line, fmtfields, fmtformat

       !--- Initialize FindFMT
       call Init_FindFMT(i_ini)
       nphas=0
       ncmd=0
       n_pat=0
       ip=i_end
       ic=0
       ipt=0
       Job_info%title=" General Job: CrysFML"
       Job_info%Num_Patterns=1

       do i=i_ini,i_end
          line=u_case(adjustl(file_dat(i)))
          if (line(1:5) == "TITLE") Job_info%title=line(7:)
          if (line(1:5) == "NPATT") then
             read(unit=line(7:), fmt=*,iostat=ier) Job_info%Num_Patterns
             if (ier /= 0) Job_info%Num_Patterns=1
          end if
          if (line(1:6) == "PHASE_") then
             nphas=nphas+1
             ip(nphas)=i
          end if
          if (line(1:4) == "CMDL") then
             ncmd=ncmd+1
             ic(ncmd)=i
          end if
          if (line(1:5) == "PATT_") then
             n_pat=n_pat+1
             ipt(n_pat)=i
          end if
       end do

       if (nphas == 0) then
          nphas=1
          ip(nphas)=0
       end if
       if (n_pat == 0) then
          n_pat=1
          ipt(n_pat) = 0
       end if

       if (Job_info%Num_Patterns /= n_pat) Job_info%Num_Patterns = n_pat
       Job_info%Num_Phases=nphas
       Job_info%Num_Cmd=ncmd

       if (allocated(Job_Info%Patt_typ)) deallocate(Job_Info%Patt_typ)
       allocate(Job_Info%Patt_typ(n_pat))

       if (allocated(Job_Info%Phas_nam)) deallocate(Job_Info%Phas_nam)
       allocate(Job_Info%Phas_nam(nphas))

       if (allocated(Job_Info%range_stl)) deallocate(Job_Info%range_stl)
       allocate(Job_Info%range_stl(n_pat))

       if (allocated(Job_Info%range_q)) deallocate(Job_Info%range_q)
       allocate(Job_Info%range_q(n_pat))

       if (allocated(Job_Info%range_d)) deallocate(Job_Info%range_d)
       allocate(Job_Info%range_d(n_pat))

       if (allocated(Job_Info%range_2theta)) deallocate(Job_Info%range_2theta)
       allocate(Job_Info%range_2theta(n_pat))

       if (allocated(Job_Info%range_energy)) deallocate(Job_Info%range_energy)
       allocate(Job_Info%range_energy(n_pat))

       if (allocated(Job_Info%range_tof)) deallocate(Job_Info%range_tof)
       allocate(Job_Info%range_tof(n_pat))

       if (allocated(Job_Info%lambda)) deallocate(Job_Info%lambda)
       allocate(Job_Info%lambda(n_pat))

       if (allocated(Job_Info%ratio)) deallocate(Job_Info%ratio)
       allocate(Job_Info%ratio(n_pat))

       if (allocated(Job_Info%dtt1)) deallocate(Job_Info%dtt1)
       allocate(Job_Info%dtt1(n_pat))

       if (allocated(Job_Info%dtt2)) deallocate(Job_Info%dtt2)
       allocate(Job_Info%dtt2(n_pat))

       !---- Initialize all variables
       Job_Info%Patt_typ    =" "
       Job_Info%Phas_nam    =" "
       Job_Info%range_stl%mina=0.0
       Job_Info%range_stl%maxb=0.0
       Job_Info%range_q%mina=0.0
       Job_Info%range_q%maxb=0.0
       Job_Info%range_d%mina=0.0
       Job_Info%range_d%maxb=0.0
       Job_Info%range_2theta%mina=0.0
       Job_Info%range_2theta%maxb=0.0
       Job_Info%range_Energy%mina=0.0
       Job_Info%range_Energy%maxb=0.0
       Job_Info%range_tof%mina=0.0
       Job_Info%range_tof%maxb=0.0
       Job_Info%Lambda%mina=0.0
       Job_Info%Lambda%maxb=0.0
       Job_Info%ratio = 0.0
       Job_Info%dtt1 = 0.0
       Job_Info%dtt2 = 0.0
       if (ncmd > 0) then
          if (allocated(Job_Info%cmd)) deallocate(Job_Info%cmd)
          allocate(Job_Info%cmd(ncmd))
          Job_Info%cmd=" "
       end if

       !---- Fill the different fields of Job_Info
       !---- Start with patterns
       fmtfields = "9fffff"

       !---- First asks if there is a PATT_ card, if not a standard is taken
       if (ipt(1) /= 0) then
          do n_pat=1, Job_info%Num_Patterns
             i=ipt(n_pat)
             line=u_case(adjustl(file_dat(i)))
             line=line(8:)
             call findfmt(0,line,fmtfields,fmtformat)
             read(unit=line,fmt=fmtformat) Job_Info%Patt_typ(n_pat), a1,a2,a3,a4,a5
             if (ierr_fmt /= 0) return
             line=u_case(Job_Info%Patt_typ(n_pat))

             select case(line(1:9))
                case("XRAY_2THE","NEUT_2THE","XRAY_SXTA","NEUT_SXTA")
                   if ( a1 <= 0.000001) a1=1.5405
                   if ( a2 <= 0.000001) then
                      a2=a1
                      a3=0.0
                   end if
                   if (a5 <= a4) a5=120.0
                   Job_Info%Lambda(n_pat)%mina=a1
                   Job_Info%Lambda(n_pat)%maxb=a2
                   Job_Info%ratio(n_pat)=a3
                   Job_Info%range_2theta(n_pat)%mina=a4
                   Job_Info%range_2theta(n_pat)%maxb=a5
                   a4=sind(0.5*a4)/a1
                   a5=sind(0.5*a5)/a2
                   Job_Info%range_stl(n_pat)%mina=a4
                   Job_Info%range_stl(n_pat)%maxb=a5
                   Job_Info%range_q(n_pat)%mina=a4*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
                   Job_Info%range_d(n_pat)%mina=0.5/a5
                   Job_Info%range_d(n_pat)%maxb=0.5/a4

                case("NEUT_TOF ")
                   if (a1 <= 0.000001) a1=1000.0
                   if (a4 <= a3) a4=2.0*abs(a3)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=a2
                   Job_Info%range_tof(n_pat)%mina=a3
                   Job_Info%range_tof(n_pat)%maxb=a4
                   Job_Info%range_d(n_pat)%mina=0.5*(-1.0+sqrt(1.0+4.0*a2*a3/a1/a1))
                   Job_Info%range_d(n_pat)%maxb=0.5*(-1.0+sqrt(1.0+4.0*a2*a4/a1/a1))
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

                case("XRAY_ENER")
                   if (a1 <= 0.000001) a1=12.4 !(=hc(keV.Angstr.)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=0.0
                   Job_Info%range_energy(n_pat)%mina=a3
                   Job_Info%range_energy(n_pat)%maxb=a4
                   if (a3 <= 0.00001) a3=0.01
                   if (a4 <= 0.00001) a4=2.00
                   Job_Info%range_d(n_pat)%mina=a1/a4
                   Job_Info%range_d(n_pat)%maxb=a1/a3
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

             end select
          end do

       else
          n_pat=1
          a1=1.5405
          a2=a1
          a3=0.0
          a4=0.0
          a5=120.0
          Job_Info%Patt_typ(n_pat)="XRAY_2THE"
          Job_Info%Lambda(n_pat)%mina=a1
          Job_Info%Lambda(n_pat)%maxb=a2
          Job_Info%ratio(n_pat)=a3
          Job_Info%range_2theta(n_pat)%mina=a4
          Job_Info%range_2theta(n_pat)%maxb=a5
          a4=sind(0.5*a4)/a1
          a5=sind(0.5*a5)/a2
          Job_Info%range_stl(n_pat)%mina=a4
          Job_Info%range_stl(n_pat)%maxb=a5
          Job_Info%range_q(n_pat)%mina=a4*4.0*pi
          Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
          Job_Info%range_d(n_pat)%mina=0.5/a5
          Job_Info%range_d(n_pat)%maxb=0.5/a4
       end if

       !---- Phase names
       if (ip(1) /= 0) then
          do i=1,nphas
             j=ip(i)
             line=adjustl(file_dat(j))
             Job_Info%Phas_nam(i)=line(8:)
          end do
       else
          Job_Info%Phas_nam(1)= Job_info%title
       end if

       !---- Command Lines, stored but not analysed here
       do i=1,ncmd
          j=ic(i)
          line=adjustl(file_dat(j))
          Job_Info%cmd(i)=line(8:)
       end do

       return
    End Subroutine Get_Job_Info

    Module Subroutine Read_CFL_SpG(lines, n_ini, n_end, SpG, xyz_type)
       !---- Arguments ----!
       character(len=*), dimension(:),  intent(in)     :: lines
       integer,                         intent(in out) :: n_ini
       integer,                         intent(in)     :: n_end
       class(SpG_Type),                 intent(out)    :: SpG
       character(len=*), optional,      intent(in)     :: xyz_type

       !--- Local Variables ---!
       integer :: i,j,ngen
       character(len=:),     allocatable :: line,uline,setting,strcode
       character(len=40), dimension(192) :: gen
       logical :: change_setting
       !Look for the appropriate keywords to construct the space group: Crystallographic, Shubnikov, or superspace

       ngen=0
       setting=" "
       change_setting=.false.
       strcode="xyz"
       if(present(xyz_type)) strcode=xyz_type
       !write(*,"(a,2i5)") " n_ini,n_end:",n_ini,n_end
       do i=n_ini,n_end
         line=trim(adjustl(lines(i)))
         if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
         j=index(line,"!")
         if(j /= 0) line=line(1:j-1)
         j=index(line,"::")
         if(j /= 0) then
           setting=trim(adjustl(line(j+2:)))
           if(len_trim(setting) /= 0) change_setting=.true.
           line=line(1:j-1)
         end if
        ! write(*,"(a)") "line: "//line
         j=index(line," ")
         uline=u_case(line(1:j-1))
         !write(*,"(a)") "uline: "//uline
         line=adjustl(line(j+1:))
         !write(*,"(a/)") "line: "//line
         Select Case(uline)
           Case("HALL","SPGR","SPACEG")
              call Set_SpaceGroup(line,SpG)
              exit
           Case("SHUB")
              call Set_SpaceGroup(line,"SHUBN",SpG)
              exit
           Case("SSG","SUPER","SSPG")
              call Set_SpaceGroup(line,"SUPER",SpG,strcode)
              exit
           Case("GENLIST","GENERATORS","LIST")
              call Set_SpaceGroup(line,SpG)
              exit
           Case("GEN","SYMM")
              ngen=ngen+1
              gen(ngen)=line
         End Select
       end do
       if(ngen > 0) then
           call Set_SpaceGroup("  ",SpG,ngen,gen)
       end if
       if(change_setting) then
         if(strcode == "xyz")  then
           call Change_Setting_SpaceG(setting, SpG)
         else
           call Change_Setting_SpaceG(setting, SpG,strcode)
         end if
       end if

    End Subroutine Read_CFL_SpG

End SubModule IOF_CFL