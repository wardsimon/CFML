!!----
!!----
!!----
SubModule (CFML_IOForm) IOF_001
   Contains
   
   !!----
   !!---- READ_SHX_ATOM
   !!----    Obtaining Atoms parameters from Shelx file (.ins or .res)
   !!----
   !!---- 26/06/2019 
   !!
   Module Subroutine Read_Shx_Atom(lines, n_ini, n_end, n_fvar, fvar, elem_type, cell, At_List)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in)      :: lines
      integer,                        intent(in out)  :: n_ini
      integer,                        intent(in)      :: n_end
      integer,                        intent(in)      :: n_fvar
      real(kind=cp), dimension(:),    intent(in)      :: fvar
      character(len=*), dimension(:), intent(in)      :: elem_type
      class(Cell_G_Type),             intent(in)      :: Cell
      type (AtList_type),             intent(out)     :: At_List

      !---- Local Variables ----!
      character(len=30),dimension(15) :: label
      character(len=80)               :: line
      character(len=2)                :: el
      integer                         :: i, j, n, nc, iv
      
      real(kind=cp)                   :: x, p, u
      real(kind=cp), dimension(15)    :: vet
      integer,       dimension(15)    :: ivet
      
      type(atlist_type)               :: Atm

      !> Init
      call clear_error()
      call Allocate_Atom_List(0, At_list)
      
      n=0
      call allocate_atom_list(n_end-n_ini+1,Atm)

      do i=n_ini,n_end
         line=lines(i)
         if (len_trim(line) == 0) cycle

         call get_words(line,label,nc)
         select case (nc)
            case (5) ! Atomname Sfac X Y Z
               call get_num(label(2),vet,ivet,iv)   ! Is Sfac integer?
               if (iv /= 1) cycle
               call get_num(label(3),vet,ivet,iv)   ! Is X real?
               if (iv /= 1) cycle
               call get_num(label(4),vet,ivet,iv)   ! Is Y real?
               if (iv /= 1) cycle
               call get_num(label(5),vet,ivet,iv)   ! Is Z real?
               if (iv /= 1) cycle

               !> Label
               n=n+1
               atm%atom(n)%lab=label(1)(1:4)
               
               !> Chemical Symbol
               call get_num(label(2),vet,ivet,iv)
               el=elem_type(ivet(1))
               atm%atom(n)%chemSymb=u_case(el(1:1))//l_case(el(2:2))
               
               !> Coordinates
               call get_num(label(3),vet,ivet,iv)
               atm%atom(n)%x(1)=vet(1)
               call get_num(label(4),vet,ivet,iv)
               atm%atom(n)%x(2)=vet(1)
               call get_num(label(5),vet,ivet,iv)
               atm%atom(n)%x(3)=vet(1)
               
               !> Thermal
               atm%atom(n)%utype="U"
               atm%atom(n)%thtype="iso"

            case (6) ! Atomname Sfac X Y Z Occ
               call get_num(label(2),vet,ivet,iv)   ! Is Sfac integer?
               if (iv /= 1) cycle
               call get_num(label(3),vet,ivet,iv)   ! Is X real?
               if (iv /= 1) cycle
               call get_num(label(4),vet,ivet,iv)   ! Is Y real?
               if (iv /= 1) cycle
               call get_num(label(5),vet,ivet,iv)   ! Is Z real?
               if (iv /= 1) cycle
               call get_num(label(6),vet,ivet,iv)   ! Is Occ real?
               if (iv /= 1) cycle

               !> Label
               n=n+1
               atm%atom(n)%lab=label(1)(1:4)
               
               !> Chemical Symbol
               call get_num(label(2),vet,ivet,iv)
               el=elem_type(ivet(1))
               atm%atom(n)%chemSymb=u_case(el(1:1))//l_case(el(2:2))
               
               !> Coordinates
               call get_num(label(3),vet,ivet,iv)
               atm%atom(n)%x(1)=vet(1)
               call get_num(label(4),vet,ivet,iv)
               atm%atom(n)%x(2)=vet(1)
               call get_num(label(5),vet,ivet,iv)
               atm%atom(n)%x(3)=vet(1)

               !> Occupancy                
               call get_num(label(6),vet,ivet,iv)
               atm%atom(n)%occ=vet(1)
               
               !> Thermal
               atm%atom(n)%utype="U"
               atm%atom(n)%thtype="iso"

            case (7,8) ! Atomname Sfac X Y Z Occ Uiso   
               !> TR: item 8 can be electronic density created by SHELXS
               call get_num(label(2),vet,ivet,iv)   ! Is Sfac integer?
               if (iv /= 1) cycle
               call get_num(label(3),vet,ivet,iv)   ! Is X real?
               if (iv /= 1) cycle
               call get_num(label(4),vet,ivet,iv)   ! Is Y real?
               if (iv /= 1) cycle
               call get_num(label(5),vet,ivet,iv)   ! Is Z real?
               if (iv /= 1) cycle
               call get_num(label(6),vet,ivet,iv)   ! Is Occ real?
               if (iv /= 1) cycle
               call get_num(label(7),vet,ivet,iv)   ! Is Uiso real?
               if (iv /= 1) cycle

               !> Label
               n=n+1
               atm%atom(n)%lab=label(1)(1:4)
               
               !> Chemical Symbol
               call get_num(label(2),vet,ivet,iv)
               el=elem_type(ivet(1))
               atm%atom(n)%chemSymb=u_case(el(1:1))//l_case(el(2:2))
               
               !> Coordinates
               call get_num(label(3),vet,ivet,iv)
               atm%atom(n)%x(1)=vet(1)
               call get_num(label(4),vet,ivet,iv)
               atm%atom(n)%x(2)=vet(1)
               call get_num(label(5),vet,ivet,iv)
               atm%atom(n)%x(3)=vet(1)

               !> Occupancy                
               call get_num(label(6),vet,ivet,iv)
               atm%atom(n)%occ=vet(1)
               
               !> U_iso
               call get_num(label(7),vet,ivet,iv)
               atm%atom(n)%U_iso=vet(1)
               
               !> Thermal
               atm%atom(n)%utype="U"
               atm%atom(n)%thtype="iso"
               
         case (9) ! Atomname Sfac X Y Z Occ U11 U22 = U33 U23 U13 U12
               call get_num(label(2),vet,ivet,iv)   ! Is Sfac integer?
               if (iv /= 1) cycle
               call get_num(label(3),vet,ivet,iv)   ! Is X real?
               if (iv /= 1) cycle
               call get_num(label(4),vet,ivet,iv)   ! Is Y real?
               if (iv /= 1) cycle
               call get_num(label(5),vet,ivet,iv)   ! Is Z real?
               if (iv /= 1) cycle
               call get_num(label(6),vet,ivet,iv)   ! Is Occ real?
               if (iv /= 1) cycle
               call get_num(label(7),vet,ivet,iv)   ! Is U11 real?
               if (iv /= 1) cycle
               call get_num(label(8),vet,ivet,iv)   ! Is U22 real?
               if (iv /= 1) cycle
               call get_num(lines(i+1),vet,ivet,iv) ! Are U33 U23 U13 U12?
               if (iv /= 4) cycle

               !> Label
               n=n+1
               atm%atom(n)%lab=label(1)(1:4)
               
               !> Chemical Symbol
               call get_num(label(2),vet,ivet,iv)
               el=elem_type(ivet(1))
               atm%atom(n)%chemSymb=u_case(el(1:1))//l_case(el(2:2))
               
               !> Coordinates
               call get_num(label(3),vet,ivet,iv)
               atm%atom(n)%x(1)=vet(1)
               call get_num(label(4),vet,ivet,iv)
               atm%atom(n)%x(2)=vet(1)
               call get_num(label(5),vet,ivet,iv)
               atm%atom(n)%x(3)=vet(1)

               !> Occupancy                
               call get_num(label(6),vet,ivet,iv)
               atm%atom(n)%occ=vet(1)
               
               !> U's (U11 U22 U33 U12 U13 U23 ) 
               call get_num(label(7),vet,ivet,iv)
               atm%atom(n)%U(1)=vet(1)
               call getnum(label(8),vet,ivet,iv)
               atm%atom(n)%U(2)=vet(1)
               
               call getnum(lines(i+1),vet,ivet,iv)
               atm%atom(n)%U(3)=vet(1)
               atm%atom(n)%U(4)=vet(4)
               atm%atom(n)%U(5)=vet(3)
               atm%atom(n)%U(6)=vet(2)
               
               !> Thermal
               atm%atom(n)%utype="U"
               atm%atom(n)%thtype="ani"
               
            case default
               cycle
         end select
      end do
      if (n <= 0) return
      
      !> Adjusting ...
      call Allocate_Atom_List(n, At_list) 
      At_List%atom=Atm%atom(1:n)
      call Allocate_Atom_List(0, Atm)

      !> Post calculations
      do i=1,n
         !> coordinates
         if (at_list%atom(i)%x(1) >= 10.0) at_list%atom(i)%x(1)=at_list%atom(i)%x(1)-10.0
         if (at_list%atom(i)%x(2) >= 10.0) at_list%atom(i)%x(2)=at_list%atom(i)%x(2)-10.0
         if (at_list%atom(i)%x(3) >= 10.0) at_list%atom(i)%x(3)=at_list%atom(i)%x(3)-10.0

         !> Occ
         if (abs(at_list%atom(i)%occ)  > 10.0) then
            x=at_list%atom(i)%occ
            if (x > 10.0) then
               at_list%atom(i)%occ=x-10.0
            else
               x=abs(at_list%atom(i)%occ)
               do j=2,n_fvar
                  if (x > 10.0*real(j) .and. x < 10.0*real(j+1)) then
                     p=x-10.0*real(j)
                     if (at_list%atom(i)%occ > 0.0) then
                        at_list%atom(i)%occ=p*fvar(j)
                     else
                        at_list%atom(i)%occ=p*(fvar(j)-1.0)
                     end if
                  end if
               end do
            end if
         end if

         !> Thermal factors
         if (at_list%atom(i)%thtype == "ani") then
            at_list%atom(i)%u_iso=U_Equiv(cell,at_list%atom(i)%u(1:6))  ! Uequi
         else
            if (at_list%atom(i)%u_iso < 0.0) then
               u=-at_list%atom(i)%u_iso
               if (u <= 5.0 .and. u >= 0.5) then
                  do j=i-1,1,-1
                     if (u_case(at_list%atom(j)%ChemSymb) == "H ") cycle
                     at_list%atom(i)%u_iso=u*U_Equiv(cell,at_list%atom(j)%u(1:6))  ! Uequi
                  end do
               end if
            end if
         end if

      end do

   End Subroutine Read_Shx_Atom
   
   !!----
   !!---- READ_SHX_CELL
   !!----    Obtaining Cell Parameter from Shelx format
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_Shx_Cell(lines, n_ini, n_end, Cell)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      class(Cell_Type),                intent(out)    :: Cell    ! Cell object

      !---- Local Variables ----!
      integer                      :: iv,z_shx
      integer,       dimension(10) :: ivet
      real(kind=cp), dimension(10) :: vet
      real(kind=cp)                :: lambda_shx
      real(kind=cp), dimension(6)  :: vcell
      real(kind=cp), dimension(6)  :: std

      !> Init
      call clear_error()

      !> Wave, CELL
      call read_key_value(lines,n_ini,n_end,"CELL",vet,ivet,iv)
      if (iv == 7) then
         lambda_shx = vet(1)
         vcell      = vet(2:7)
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_SHX_Cell@CFML_IOForm: Problems reading cell parameters!"  
         return 
      end if

      !> Z, STD 
      call read_key_value(lines,n_ini,n_end,"ZERR",vet,ivet,iv)
      if (iv == 7) then
         z_shx= ivet(1)
         std  = vet(2:7)
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_SHX_Cell@CFML_IOForm: Problems reading sigma(cell) parameters!"  
         return   
      end if

      call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, Vscell=std(1:3), Vsang=std(4:6))
      
   End Subroutine Read_Shx_Cell
   
   !!----
   !!---- READ_SHX_FVAR
   !!----    Obtaining Fvar parameters from Shelx file (.ins or .res)
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_Shx_Fvar(Lines,n_ini,n_end, n_fvar, fvar)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in)    :: lines
      integer,                        intent(in out):: n_ini
      integer,                        intent(in)    :: n_end
      integer,                        intent(out)   :: n_fvar
      real(kind=cp), dimension(:),    intent(out)   :: fvar

      !---- Local  variables ----!
      integer                      :: iv
      integer,       dimension(15) :: ivet
      real(kind=cp), dimension(15) :: vet

      !> Init
      call clear_error()
      n_fvar = 1
      fvar   = 1.0_cp

      call read_key_value(lines,n_ini,n_end,"FVAR",vet,ivet,iv)
      if (iv /= 0) then
         n_fvar=iv
         fvar=vet(1:iv)
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_SHX_FVar@CFML_IOForm: Problems reading the FVAR parameters!"   
      end if
   End Subroutine Read_Shx_Fvar
   
   !!----
   !!---- READ_SHX_WAVE
   !!----    Obtaining wavelength Parameter from Shelx format
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_Shx_Wave(lines, n_ini, n_end, Wave)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      real(kind=cp),                   intent(out)    :: Wave    ! Wavelength
      
      !---- Local Variables ----!
      integer                      :: iv
      integer,       dimension(10) :: ivet
      real(kind=cp), dimension(10) :: vet

      !> Init
      call clear_error()

      !> Wave, CELL
      call read_key_value(lines,n_ini,n_end,"CELL",vet,ivet,iv)
      if (iv == 7) then
         wave = vet(1)
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_SHX_Wave@CFML_IOForm: Problems reading wavelength parameter!"  
         return 
      end if
   End Subroutine Read_Shx_Wave
   
   !!----
   !!---- READ_SHX_Z
   !!----    Obtaining Cell Parameter from Shelx format
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_Shx_Z(lines, n_ini, n_end, Z)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      real(kind=cp),                   intent(out)    :: Z       ! Z number

      !---- Local Variables ----!
      integer                      :: iv
      integer,       dimension(10) :: ivet
      real(kind=cp), dimension(10) :: vet

      !> Init
      Z=0.0_cp
      call clear_error()

      !> Z 
      call read_key_value(lines,n_ini,n_end,"ZERR",vet,ivet,iv)
      if (iv == 7) then
         z= vet(1)
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_SHX_Z@CFML_IOForm: Problems reading Z value!"  
      end if
   End Subroutine Read_Shx_Z
   
   !!----
   !!---- READ_SHX_TITL
   !!----    Obtaining Title from Shelx file
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_Shx_Titl(lines,n_ini,n_end,Title)
      !---- Arguments ----!
      character(len=*),dimension(:), intent(in)     :: lines
      integer,                       intent(in out) :: n_ini
      integer,                       intent(in)     :: n_end
      character(len=*),              intent(out)    :: title

      !> Init
      Title=" "
      call Read_Key_StrVal(lines,n_ini,n_end,"TITL",title)
   End Subroutine Read_Shx_Titl
   
   !!----
   !!---- READ_SHX_LATT
   !!----    Obtaining lattice from Shelx file (.ins or .res)
   !!----
   !!----   Lattice
   !!----       1=P                     4=F              7=C
   !!----       2=I                     5=A
   !!----       3=rhombohedral obverse  6=B
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_Shx_Latt(lines,n_ini,n_end,latt)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in) :: lines
      integer,           intent(in out)          :: n_ini
      integer,           intent(in)              :: n_end
      integer,           intent(out)             :: latt

      !---- Local Variables ----!
      integer                     :: iv
      integer,       dimension(2) :: ivet
      real(kind=cp), dimension(2) :: vet

      !> Init
      latt=1
      
      call read_key_value(lines,n_ini,n_end,"LATT",vet,ivet,iv)
      if (iv == 1) latt = abs(ivet(1))
      
   End Subroutine Read_Shx_Latt
   
   !!----
   !!---- READ_SHX_CONT
   !!----    Obtaining Chemical contents from Shelx file (.ins or .res)
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_Shx_Cont(lines,n_ini,n_end, n_elem_type, elem_type, n_elem)
      !---- Arguments ----!
      character(len=*), dimension(:),           intent(in)      :: lines
      integer,                                  intent(in out)  :: n_ini
      integer,                                  intent(in)      :: n_end
      integer,                                  intent(out)     :: n_elem_type
      character(len=*), dimension(:),           intent(out)     :: elem_type
      real(kind=cp),    dimension(:), optional, intent(out)     :: n_elem

      !---- Local  variables ----!
      character(len=80)             :: line
      integer                       :: iv
      integer,       dimension(15)  :: ivet
      real(kind=cp), dimension(15)  :: vet

      !> Init
      call clear_error()
      n_elem_type = 0
      elem_type   = " "
      if (present(n_elem)) n_elem = 0.0_cp

      !> Define the elements types: C H N O
      call Read_Key_StrVal(lines,n_ini,n_end,"SFAC",line)
      if (len_trim(line) /=0) then
         call get_words(line,elem_type,n_elem_type)
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_SHX_CONT@CFML_IOForm: Problems reading the SFAC information!" 
         return   
      end if

      if (present(n_elem)) then
         call read_key_value(lines,n_ini,n_end,"UNIT",vet,ivet,iv)
         if (iv > 0) then
            n_elem=vet
         else
            err_CFML%Ierr=1
         err_CFML%Msg="Read_SHX_CONT@CFML_IOForm: Problems reading the UNITS!" 
         return
         end if   
      end if

   End Subroutine Read_Shx_Cont
   
   !!----
   !!---- Read_Shx_Symm
   !!----    Obtaining Symmetry Operators from Shelx file (.ins or .res)
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_Shx_Symm(lines,n_ini,n_end,n_oper,oper_symm)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in) :: lines
      integer,          intent(in out)           :: n_ini
      integer,          intent(in)               :: n_end
      integer,          intent(out)              :: n_oper
      character(len=*), dimension(:),intent(out) :: oper_symm

      !---- Local variables ----!
      character(len=80) :: line
      integer           :: nline

      !> Init
      call clear_error()
      n_oper=0
      oper_symm=" "

      do
         call Read_Key_StrVal(lines,n_ini,n_end,"SYMM",line)
         if (len_trim(line) /=0) then
            n_oper=n_oper+1
            oper_symm(n_oper)=trim(line)
            n_ini=n_ini+1
            nline=n_ini
         else
            exit
         end if
      end do
      n_ini=nline
   End Subroutine Read_Shx_Symm
    
   !!----
   !!---- WRITE_SHX_TEMPLATE
   !!----    Write a Shelx File
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Write_Shx_Template(filename,code,title,lambda,z,cell,spg,At_List)
      !---- Arguments ----!
      character(len=*),        intent(in) :: filename
      integer,                 intent(in) :: code        ! 0: Patterson, 1: DM   2: Refinement
      character(len=*),        intent(in) :: title
      real(kind=cp),           intent(in) :: lambda
      integer,                 intent(in) :: z
      class(cell_Type),        intent(in) :: cell
      class(SpG_Type),         intent(in) :: SpG
      type(atlist_type),       intent(in) :: at_List

      !---- Local Variables ----!
      logical                :: info

      integer                :: i,j,k,nc,iunit
      integer                :: nlat
      integer, dimension(15) :: z_cont

      !> Init
      info=.false.
      z_cont=0
      nc=0  !this depends on scattering factor?

      !> Ask if the file was oppened
      inquire(file=trim(filename),opened=info)
      if (info) then
         inquire(file=trim(filename),number=iunit)
         close(unit=iunit)
      end if

      !> Open
      open(newunit=iunit,file=trim(filename),status="unknown",action="write")
      rewind(unit=iunit)

      !> Title
      write(unit=iunit,fmt="(a)") "TITL "//trim(title)

      !> Lambda, Cell
      write(unit=iunit,fmt="(a,f8.5,3f9.4,3f7.3)") "CELL ",lambda,cell%cell,cell%ang

      !> Z, Std
      write(unit=iunit,fmt="(a,i3,a,3f8.4,3f7.3)") "ZERR ",z,"     ",cell%scell,cell%sang

      !> Latt
      nlat=1
      select case (spg%centred)
         case (0) ! Centric

         case (1) ! Acentric
            nlat=-1

         case (2) ! Not used in Shelx
            err_CFML%IErr=1
            err_CFML%Msg="Write_SHX_Template@CFML_IOFORM: Origin not at -1!"
            close(unit=iunit)
            return
      end select

      select case (spg%spg_lat)
         case ("I")
            nlat=2*nlat
         case ("R")
            nlat=3*nlat
         case ("F")
            nlat=4*nlat
         case ("A")
            nlat=5*nlat
         case ("B")
            nlat=6*nlat
         case ("C")
            nlat=7*nlat
      end select
      write(unit=iunit,fmt="(a,i2)") "LATT ",nlat

      !> Symm
      do i=2,spg%numops
         write(unit=iunit,fmt="(a)") "SYMM "//u_case(trim(spg%Symb_Op(i)))
      end do

      !> Sfac
      j=0
      do i=1,At_List%natoms
         if (j == 0) then
            j=1
            z_cont(j)=at_list%atom(i)%z
         else
            do k=1,j
               if (z_cont(k) == at_list%atom(i)%z) exit
            end do
            if (z_cont(k) /= at_list%atom(i)%z) then
               j=j+1
               z_cont(j)=at_List%atom(i)%z
            end if
         end if
      end do
      write(unit=iunit,fmt="(a)") "SFAC "

      !> Unit
      write(unit=iunit,fmt="(a)") "UNIT "

      select case (code)
         case (0) ! Shelxs - Patterson
            write(unit=iunit,fmt="(a)") "PATT "

         case (1) ! Shelxs - Direct Methods
            write(unit=iunit,fmt="(a)") "TREF "

         case (2) ! Shelxl - Refinement
            !> L.S.
            write(unit=iunit,fmt="(a)") "L.S. 10"

            !> Fvar
            write(unit=iunit,fmt="(a)") "FVAR 1.0"

            !> Weight
            write(unit=iunit,fmt="(a)") "WGHT 0.2"

            !> Fmap
            write(unit=iunit,fmt="(a)") "FMAP 2"

            !> Atoms
            do i=1,at_list%natoms
               write(unit=iunit,fmt="(a4,i3,4f11.5)") &
                    at_list%atom(i)%lab, nc, at_list%atom(i)%x, at_list%atom(i)%occ+10.0_cp
            end do
      end select

      !> Format
      write(unit=iunit,fmt="(a)") "HKLF 4"

      !>- End
      write(unit=iunit,fmt="(a)") "END "

   End Subroutine Write_Shx_Template
   
End SubModule IOF_001   