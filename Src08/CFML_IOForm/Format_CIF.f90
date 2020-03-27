!!----
!!----
!!----
SubModule (CFML_IOForm) IOF_003
   Contains


   !!----
   !!---- Read_Cif_Atom
   !!----    Obtaining Atoms parameters from Cif file. A control error is present.
   !!----
   !!---- 26/06/2019
   !!
   Module Subroutine Read_Cif_Atom(lines,n_ini,n_end, At_List)
      !---- Arguments ----!
      character(len=*), dimension(:),   intent(in)      :: lines
      integer,                          intent(in out)  :: n_ini
      integer,                          intent(in)      :: n_end
      type (atList_type),               intent(out)     :: At_List

      !---- Local Variables ----!
      character(len=12), parameter        :: DIGPM="0123456789+-"
      character(len=132), allocatable     :: line
      character(len=20),dimension(15)     :: label
      integer                             :: i, j, n, nc, nct, nline, iv, First, nline_big,num_ini,mm
      integer, dimension( 8)              :: lugar   !   1 -> label
                                                     !   2 -> Symbol
                                                     ! 3-5 -> coordinates
                                                     !   6 -> occupancy
                                                     !   7 -> Uequi
                                                     !   8 -> Biso
      real(kind=cp), dimension(1)     :: vet1,vet2
      integer,       dimension(1)     :: ivet
      type(atlist_type)               :: Atm

      type (atm_type)                 :: atm1
      type (atm_std_type)             :: atm2
      type (matm_std_type)            :: atm3
      type (atm_ref_type)             :: atm4

      class(atm_type), allocatable    :: atm5

      !> Init
      call clear_error()
      call allocate_atom_list(0,At_List)
      call allocate_atom_list(n_end-n_ini+1,Atm)

      lugar=0
      num_ini=n_ini

      !> Change of _atom_site by _atom_site_label in order to be able the reading
      !> the atoms positions even when the anisotropic parameters are given before
      call Read_Key_StrVal(lines,n_ini,n_end,"_atom_site_label",line)

      !> Look for the possibility that _atom_site_label is not the first item in the loop
      do i=n_ini,num_ini,-1
         line=adjustl(lines(i))
         if (line(1:) == "loop_") then
            n_ini=i+1
            exit
         end if
      end do
      j=0
      do i=n_ini,n_end
         line=adjustl(lines(i))
         if ("_atom_site_label" == line(1:16)) then
            j=j+1
            lugar(1)=j
            cycle
         end if
         if ("_atom_site_type_symbol" == line(1:22)) then
            j=j+1
            lugar(2)=j
            cycle
         end if
         if ("_atom_site_fract_x" == line(1:18)) then
            j=j+1
            lugar(3)=j
            cycle
         end if
         if ("_atom_site_fract_y" == line(1:18)) then
            j=j+1
            lugar(4)=j
            cycle
         end if
         if ("_atom_site_fract_z" == line(1:18)) then
            j=j+1
            lugar(5)=j
            cycle
         end if
         if ("_atom_site_occupancy" == line(1:20)) then
            j=j+1
            lugar(6)=j
            cycle
         end if
         if ("_atom_site_U_iso_or_equiv" == line(1:25)) then
            j=j+1
            lugar(7)=j
            cycle
         end if
         if ("_atom_site_B_iso_or_equiv" == line(1:25)) then
            j=j+1
            lugar(8)=j
            cycle
         end if
         if ("_atom_site_" == line(1:11)) then
            j=j+1
            cycle
         end if

         if ("_oxford_atom_site_" == line(1:18)) then
            j=j+1
            cycle
         end if

         nline=i
         exit
      end do

      if (any(lugar(3:5) == 0)) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_Cif_Atom@CFML_IOFORM: Error reading atoms in CIF format!"
         return
      end if

      nct=count(lugar > 0)
      nline_big=nline
      nline_ini=nline
      n=0
      do i=n_ini,n_end
         line=adjustl(lines(i))
         if (line(1:1) == "#" .or. line(1:1) == "?") cycle
         if (len_trim(line) == 0) exit
         if (line(1:1) == "_" .or. line(1:5) == "loop_") exit
         call get_words(line, label, nc)
         if (nc < nct) then
            nline=i
            exit
         end if

         n=n+1

         !> _atom_site_label
         atm%atom(n)%lab=label(lugar(1))

         !> _atom_site_type_symbol
         if (lugar(2) /= 0) then
            atm%atom(n)%SfacSymb=label(lugar(2))(1:4)
            if (index(DIGPM, label(lugar(2))(2:2)) /= 0 ) then
               atm%atom(n)%chemSymb=u_case(label(lugar(2))(1:1))
            else
               atm%atom(n)%chemSymb=u_case(label(lugar(2))(1:1))//l_case(label(lugar(2))(2:2))
            end if
         else
            if(index(DIGPM,label(lugar(1))(2:2)) /= 0 ) then
               atm%atom(n)%chemSymb=u_case(label(lugar(1))(1:1))
            else
               atm%atom(n)%chemSymb=u_case(label(lugar(1))(1:1))//l_case(label(lugar(1))(2:2))
            end if
            atm%atom(n)%SfacSymb=atm%atom(n)%chemSymb
         end if

         select type (at => atm%atom)
            type is (atm_type)
               call get_num(label(lugar(3)),vet1,ivet,iv)    ! _atom_site_fract_x
               at(n)%x(1)=vet1(1)
               call get_num(label(lugar(4)),vet1,ivet,iv)    ! _atom_site_fract_y
               at(n)%x(2)=vet1(1)
               call get_num(label(lugar(5)),vet1,ivet,iv)    ! _atom_site_fract_z
               at(n)%x(3)=vet1(1)

            class is (atm_std_type)
               call get_numstd(label(lugar(3)),vet1,vet2,iv)    ! _atom_site_fract_x
               at(n)%x(1)=vet1(1)
               at(n)%x_std(1)=vet2(1)
               call get_numstd(label(lugar(4)),vet1,vet2,iv)    ! _atom_site_fract_y
               at(n)%x(2)=vet1(1)
               at(n)%x_std(2)=vet2(1)
               call get_numstd(label(lugar(5)),vet1,vet2,iv)    ! _atom_site_fract_z
               at(n)%x(3)=vet1(1)
               at(n)%x_std(3)=vet2(1)
         end select

         !> _atom_site_occupancy
         if (lugar(6) /= 0) then
            call get_numstd(label(lugar(6)),vet1,vet2,iv)
         else
            vet1=1.0
            vet2=0.0_cp
         end if
         select type (at => atm%atom)
            type is (atm_type)
               at(n)%occ=vet1(1)

            class is (atm_std_type)
               at(n)%occ=vet1(1)
               at(n)%occ_std=vet2(1)
         end select

         if (lugar(7) /= 0) then
            call get_numstd(label(lugar(7)),vet1,vet2,iv)    ! _atom_site_U_iso_or_equiv
            select type (at => atm%atom)
               type is (atm_type)
                  at(n)%U_iso=vet1(1)
               class is (atm_std_type)
                  at(n)%U_iso=vet1(1)
                  at(n)%U_iso_std=vet2(1)
            end select
            atm%atom(n)%utype="U"

         else if (lugar(8) /= 0) then
            call get_numstd(label(lugar(8)),vet1,vet2,iv)    ! _atom_site_B_iso_or_equiv
            select type (at => atm%atom)
               type is (atm_type)
                  at(n)%U_iso=vet1(1)
               class is (atm_std_type)
                  at(n)%U_iso=vet1(1)
                  at(n)%U_iso_std=vet2(1)
            end select
            atm%atom(n)%utype="B"
         else
            select type (at => atm%atom)
               type is (atm_type)
                  at(n)%U_iso=0.0_cp
               class is (atm_std_type)
                  at(n)%U_iso=0.0_cp
                  at(n)%U_iso_std=0.0_cp
            end select
            atm%atom(n)%utype="U"
         end if
      end do

      if (nline >= nline_big) nline_big=nline

      !> Anisotropic
      nline_ini=num_ini !Changed to be able the reading of anisotropic parameters
                        !even if given before the coordinates
      lugar=0
      call Read_Key_StrVal(lines,n_ini,n_end,"_atom_site_aniso_", line)

      j=0
      do i=n_ini,n_end
         line=adjustl(lines(i))
         if ("_atom_site_aniso_label" == line(1:22)) then
            j=j+1
            lugar(1)=j
            cycle
         end if
         if ("_atom_site_aniso_type_symbol" == line(1:28)) then
            j=j+1
            lugar(8)=j
            cycle
         end if
         if ("_atom_site_aniso_U_11" == line(1:21)) then
            j=j+1
            lugar(2)=j
            cycle
         end if
         if ("_atom_site_aniso_U_22" == line(1:21)) then
            j=j+1
            lugar(3)=j
            cycle
         end if
         if ("_atom_site_aniso_U_33" == line(1:21)) then
            j=j+1
            lugar(4)=j
            cycle
         end if
         if ("_atom_site_aniso_U_12" == line(1:21)) then
            j=j+1
            lugar(5)=j
            cycle
         end if
         if ("_atom_site_aniso_U_13" == line(1:21)) then
            j=j+1
            lugar(6)=j
            cycle
         end if
         if ("_atom_site_aniso_U_23" == line(1:21)) then
            j=j+1
            lugar(7)=j
            cycle
         end if

         if ("_atom_site_aniso" == line(1:16) ) then
            j=j+1
            cycle
         end if

         nline=i
         exit
      end do
      if(nline >= nline_big) nline_big=nline
      if (all(lugar(1:7) > 0)) then        ! T.R. June 2017
         nct=count(lugar > 0)
         nline_ini=nline
         mm=0
         do i=n_ini,n_end
            line=adjustl(lines(i))
            if (line(1:1) == "#" .or. line(1:1) =="?") cycle
            if (len_trim(line) == 0) exit
            call getword(line,label,nc)
            if (nc < nct) then
               nline=i
               exit
            end if
            do j=1,n
               if (atm%atom(j)%thtype == "ani") cycle ! already assigned
               if (trim(atm%atom(j)%lab) /= trim(label(lugar(1))) ) cycle

               call get_numstd(label(lugar(2)),vet1,vet2,iv)    ! _atom_site_aniso_U_11
               select type (at => atm%atom)
                  type is (atm_type)
                     at(j)%u(1)    =vet1(1)
                  class is (atm_std_type)
                     at(j)%u(1)    =vet1(1)
                     at(j)%u_std(1)=vet2(1)
               end select

               call get_numstd(label(lugar(3)),vet1,vet2,iv)    ! _atom_site_aniso_U_22
               select type (at => atm%atom)
                  type is (atm_type)
                     at(j)%u(2)    =vet1(1)
                  class is (atm_std_type)
                     at(j)%u(2)    =vet1(1)
                     at(j)%u_std(2)=vet2(1)
               end select

               call get_numstd(label(lugar(4)),vet1,vet2,iv)    ! _atom_site_aniso_U_33
               select type (at => atm%atom)
                  type is (atm_type)
                     at(j)%u(3)    =vet1(1)
                  class is (atm_std_type)
                     at(j)%u(3)    =vet1(1)
                     at(j)%u_std(3)=vet2(1)
               end select

               call get_numstd(label(lugar(5)),vet1,vet2,iv)    ! _atom_site_aniso_U_12
               select type (at => atm%atom)
                  type is (atm_type)
                     at(j)%u(4)    =vet1(1)
                  class is (atm_std_type)
                     at(j)%u(4)    =vet1(1)
                     at(j)%u_std(4)=vet2(1)
               end select

               call get_numstd(label(lugar(6)),vet1,vet2,iv)    ! _atom_site_aniso_U_13
               select type (at => atm%atom)
                  type is (atm_type)
                     at(j)%u(5)    =vet1(1)
                  class is (atm_std_type)
                     at(j)%u(5)    =vet1(1)
                     at(j)%u_std(5)=vet2(1)
               end select

               call get_numstd(label(lugar(7)),vet1,vet2,iv)    ! _atom_site_aniso_U_23
               select type (at => atm%atom)
                  type is (atm_type)
                     at(j)%u(6)    =vet1(1)
                  class is (atm_std_type)
                     at(j)%u(6)    =vet1(1)
                     at(j)%u_std(6)=vet2(1)
               end select

               atm%atom(j)%thtype="ani"
               mm=mm+1
               exit
            end do
         end do

      end if
      if (nline >= nline_big) nline_big=nline
      n_ini=nline_big

      !> Look for the first atoms fully occupying the site and put it in first position
      !> This is needed for properly calculating the occupation factors
      !> after normalization in subroutine Readn_Set_XTal_CIF
      vet1(1)=maxval(atm%atom(1:n)%occ)  !Normalize occupancies
      atm%atom%occ=atm%atom%occ/vet1(1)
      First=1
      do i=1,n
         if (abs(atm%atom(i)%occ-1.0_cp) < EPS) then
            First=i
            exit
         end if
      end do

      !> Swapping the orinal atom at the first position with the first having full occupation
      if (First /= 1) Then
         select type (at => atm%atom)
            type is (atm_type)
               atm1=at(1)
               at(1)=at(first)
               at(first)=atm1

            type is (atm_std_type)
            type is (matm_std_type)
            type is (atm_ref_type)
         end select
      end if

      !> Put the first atom the first having a full occupation factor 1.0
      if (n > 0) then
         call allocate_atom_list(n,At_list)
         at_list%atom=atm%atom(1:n)
      end if

      call allocate_atom_list(0,Atm)

   End Subroutine Read_Cif_Atom

   !!----
   !!---- Read_Cif_Cell
   !!----    Read Cell Parameters from CIF format
   !!----
   !!---- Update: February - 2005
   !!
   Module Subroutine Read_Cif_Cell(lines,N_Ini,N_End,Cell)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      class(Cell_Type),                intent(out)    :: Cell    ! Cell object

      !---- Local Variables ----!
      integer                     :: iv,initl
      real(kind=cp), dimension(1) :: vet1,vet2
      real(kind=cp), dimension(6) :: vcell, std
      logical                     :: ierror

      !> Init
      call clear_error()
      vcell=0.0_cp; std=0.0_cp
      ierror=.false.

      !> Celda
      initl=n_ini  ! Preserve initial line => some CIF files have random order for cell parameters
      call read_key_valueSTD(lines,n_ini,n_end,"_cell_length_a",vet1,vet2,iv)
      if (iv == 1) then
         vcell(1)=vet1(1)
         std(1)=vet2(1)
      else
         ierror=.true.
      end if

      nline_ini=initl
      call read_key_valueSTD(lines,n_ini,n_end,"_cell_length_b",vet1,vet2,iv)
      if (iv == 1) then
         vcell(2)=vet1(1)
         std(2)=vet2(1)
      else
         ierror=.true.
      end if

      nline_ini=initl
      call read_key_valueSTD(lines,n_ini,n_end,"_cell_length_c",vet1,vet2,iv)
      if (iv == 1) then
         vcell(3)=vet1(1)
         std(3)=vet2(1)
      else
         ierror=.true.
      end if

      nline_ini=initl
      call read_key_valueSTD(lines,n_ini,n_end,"_cell_angle_alpha",vet1,vet2,iv)
      if (iv == 1) then
         vcell(4)=vet1(1)
         std(4)=vet2(1)
      else
         ierror=.true.
      end if

      nline_ini=initl
      call read_key_valueSTD(lines,n_ini,n_end,"_cell_angle_beta",vet1,vet2,iv)
      if (iv == 1) then
         vcell(5)=vet1(1)
         std(5)=vet2(1)
      else
         ierror=.true.
      end if

      nline_ini=initl
      call read_key_valueSTD(lines,n_ini,n_end,"_cell_angle_gamma",vet1,vet2,iv)
      if (iv == 1) then
         vcell(6)=vet1(1)
         std(6)=vet2(1)
      else
         ierror=.true.
      end if

      if (ierror) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CIF_Cell@CFML_IOForm: Problems reading cell parameters!"
         return
      end if
      call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, Vscell=std(1:3), Vsang=std(4:6))
   End Subroutine Read_Cif_Cell

   !!----
   !!---- READ_CIF_WAVE
   !!----    Read Wavelength in CIF Format
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_Wave(lines, n_ini, n_end, Wave)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      real(kind=cp),                   intent(out)    :: Wave    ! Wavelength

      !---- Local Variables ----!
      integer                    :: iv
      integer,dimension(1)       :: ivet
      real(kind=cp), dimension(1):: vet

      !> Init
      wave=0.0_cp
      call clear_error()

      call read_key_value(lines,n_ini,n_end, "_diffrn_radiation_wavelength",vet,ivet,iv)
      if (iv == 1) then
         wave=vet(1)
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CIF_Wave@CFML_IOForm: Problems reading wavelenth value!"
      end if
   End Subroutine Read_Cif_Wave

   !!----
   !!---- READ_CIF_Z
   !!----    Unit formula from Cif file
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_Z(lines, n_ini, n_end, Z)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      real(kind=cp),                   intent(out)    :: Z       ! Z number

      !---- Local Variables ----!
      integer                     :: iv
      integer,dimension(1)        :: ivet
      real(kind=cp), dimension(1) :: vet

      !> Init
      call clear_error()
      z=0

      call read_key_value(lines,n_ini,n_end, "_cell_formula_units_Z",vet,ivet,iv)
      if (iv == 1) then
         z=vet(1)
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CIF_Z@CFML_IOForm: Problems reading Z value!"
      end if
   End Subroutine Read_Cif_Z

   !!----
   !!---- Read_Cif_ChemicalName
   !!----    Obtaining Chemical Name from Cif file
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_ChemName(lines,N_ini,N_End,ChemName)
      !---- Arguments ----!
      character(len=*),  dimension(:), intent(in) :: lines
      integer,           intent(in out)           :: n_ini
      integer,           intent(in)               :: n_end
      character(len=*),  intent(out)              :: ChemName

      !---- Local variables ----!
      integer :: np1, np2

      !> Init
      ChemName=" "

      call Read_Key_StrVal(lines,n_ini,n_end, "_chemical_name_common",ChemName)

      if (len_trim(chemname) == 0) then
         call Read_Key_StrVal(lines,n_ini,n_end, "_chemical_name_systematic",ChemName)
      end if

      if (len_trim(chemname) > 0) then
         if (trim(chemname) =="; ?" .or. trim(chemname)=="#") chemname=" "
         np1=index(chemname,"'")
         np2=index(chemname,"'",back=.true.)
         if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
            chemname=chemname(np1+1:np2-1)
         else
            np1=index(chemname,'"')
            np2=index(chemname,'"',back=.true.)
            if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
               chemname=chemname(np1+1:np2-1)
            end if
         end if
      end if

   End Subroutine Read_Cif_ChemName

   !!----
   !!---- READ_CIF_CONT
   !!----    Obtaining the chemical contents from Cif file
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_Cont(lines,N_Ini,N_End,N_Elem_Type,Elem_Type,N_Elem)
      !---- Arguments ----!
      character(len=*), dimension(:),      intent(in)      :: lines
      integer,                             intent(in out)  :: n_ini
      integer,                             intent(in)      :: n_end
      integer,                             intent(out)     :: n_elem_type
      character(len=*), dimension(:),      intent(out)     :: elem_type
      real(kind=cp), dimension(:),optional,intent(out)     :: n_elem

      !---- Local  variables ----!
      character(len=132)                   :: line
      character(len=10), dimension(15)     :: label

      integer                     :: iv
      integer                     :: i,np1,np2,nlabel,nlong
      integer,       dimension(1) :: ivet
      real(kind=cp), dimension(1) :: vet

      !> Init
      call clear_error()

      n_elem_type = 0
      elem_type   = " "
      if (present(n_elem)) n_elem = 0.0_cp

      call Read_Key_StrVal(lines,n_ini,n_end, "_chemical_formula_sum",line)
      if (len_trim(line) ==0) line=lines(n_ini+1)
      line=adjustl(line)
      if (line(1:1) == "?") return

      np1=index(line,"'")
      np2=index(line,"'",back=.true.)
      nlabel=0
      if (np1 /= 0 .and. np2 /= 0 .and. np2 > np1) then
         call get_words(line(np1+1:np2-1), label, nlabel)
      end if
      if (nlabel /=0) then
         n_elem_type = nlabel
         do i=1,nlabel
            nlong=len_trim(label(i))
            select case (nlong)
                case (1)
                   elem_type(i)=label(i)(1:1)
                   if (present(n_elem)) n_elem(i)   = 1.0_cp

                case (2)
                   call get_num(label(i)(2:),vet,ivet,iv)
                   if (iv == 1) then
                      elem_type(i)=label(i)(1:1)
                      if (present(n_elem)) n_elem(i)   =vet(1)
                   else
                      elem_type(i)=label(i)(1:2)
                      if (present(n_elem)) n_elem(i)   = 1.0_cp
                   end if

                case (3:)
                   call get_num(label(i)(2:),vet,ivet,iv)
                   if (iv == 1) then
                      elem_type(i)=label(i)(1:1)
                      if (present(n_elem)) n_elem(i)   =vet(1)
                   else
                      call get_num(label(i)(3:),vet,ivet,iv)
                      if (iv == 1) then
                         elem_type(i)=label(i)(1:2)
                         if (present(n_elem)) n_elem(i)   =vet(1)
                      else
                         elem_type(i)=label(i)(1:2)
                         if (present(n_elem)) n_elem(i)   = 1.0
                      end if
                   end if
            end select
         end do
      end if

   End Subroutine Read_Cif_Cont

   !!----
   !!---- READ_CIF_PRESSURE
   !!----    Pressure and Sigma
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_Pressure(lines,N_ini,N_End, P, SigP)
      !---- Arguments ----!
      character(len=*),  dimension(:), intent(in) :: lines
      integer,           intent(in out)           :: n_ini
      integer,           intent(in)               :: n_end
      real(kind=cp),     intent(out)              :: p
      real(kind=cp),     intent(out)              :: sigp

      !---- Local Variables ----!
      integer                    :: iv
      real(kind=cp),dimension(1) :: vet1,vet2

      !> Init
      p=0.0_cp
      sigp=1.0e-5

      call read_key_valuestd(lines,n_ini,n_end, "_diffrn_ambient_pressure",vet1,vet2,iv)
      if (iv == 1) then
         p=vet1(1)*1.0e-6
         sigp=vet2(1)*1.0e-6
      end if

   End Subroutine Read_Cif_Pressure

   !!----
   !!---- Read_Cif_Title
   !!----    Obtaining Title from Cif file
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_Title(lines,N_Ini,N_End,Title)
      !---- Arguments ----!
      character(len=*),  dimension(:), intent(in) :: lines
      integer,           intent(in out)           :: n_ini
      integer,           intent(in)               :: n_end
      character(len=*),  intent(out)              :: title

      !---- Local variables ----!
      integer :: np, np1, np2

      !> Init
      title=" "
      call Read_Key_StrVal(lines,n_ini,n_end, "_publ_section_title",title)

      if (len_trim(title) ==0 ) title=adjustl(lines(n_ini+1))
      if (title =="; ?" .or. title=="#") then
         title=" "
      else
         np=len_trim(title)
         if (np <= 3) title=adjustl(lines(n_ini+2))
         np1=index(title,"'")
         np2=index(title,"'",back=.true.)
         if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
            title=title(np1+1:np2-1)
         else
            np1=index(title,'"')
            np2=index(title,'"',back=.true.)
            if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
               title=title(np1+1:np2-1)
            end if
         end if
      end if

   End Subroutine Read_Cif_Title

   !!----
   !!---- READ_CIF_TEMP
   !!----    Temperature and Sigma in Kelvin
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_Temp(lines,N_Ini,N_End,T,SigT)
      !---- Arguments ----!
      character(len=*),  dimension(:), intent(in) :: lines
      integer,           intent(in out)           :: n_ini
      integer,           intent(in)               :: n_end
      real(kind=cp),     intent(out)              :: T
      real(kind=cp),     intent(out)              :: sigT

      !---- Local Variables ----!
      integer                    :: iv
      real(kind=cp),dimension(1) :: vet1,vet2

      !> Init
      T=298.0_cp
      sigt=1.0_cp

      call read_key_valuestd(lines,n_ini,n_end, "_diffrn_ambient_temperature",vet1,vet2,iv)
      if (iv == 1) then
         t=vet1(1)
         sigt=vet2(1)
      end if

   End Subroutine Read_Cif_Temp

   !!----
   !!---- Read_Cif_Hall
   !!----    Obtaining the Hall symbol of the Space Group
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_Hall(lines, N_Ini, N_End, Hall)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in) :: lines
      integer,          intent(in out)           :: n_ini
      integer,          intent(in)               :: n_end
      character(len=*), intent(out)              :: Hall

      !---- Local variables ----!
      integer :: np1, np2

      !> Init
      Hall=" "
      call Read_Key_StrVal(lines,n_ini,n_end, "_symmetry_space_group_name_Hall",hall)
      if (len_trim(Hall)==0) Hall=adjustl(lines(n_ini+1))

      !> TR  feb. 2015 .(re-reading the same item with another name)
      if (len_trim(Hall) == 0) then
         call Read_Key_StrVal(lines,n_ini,n_end, "_space_group_name_Hall",hall)
         if (len_trim(Hall)==0) Hall=adjustl(lines(n_ini+1))
      end if

      if (trim(Hall) =="?" .or. trim(Hall)=="#") then
         Hall=" "
      else
         np1=index(Hall,"'")
         np2=index(Hall,"'",back=.true.)
         if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
            Hall=Hall(np1+1:np2-1)
         else
            np1=index(Hall,'"')
            np2=index(Hall,'"',back=.true.)
            if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
               Hall=Hall(np1+1:np2-1)
            else
               Hall=" "
            end if
         end if
      end if

   End Subroutine Read_Cif_Hall

   !!----
   !!---- READ_CIF_HM
   !!----    Obtaining the Herman-Mauguin symbol of Space Group
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_HM(lines, N_Ini, N_End, Spgr_Hm)
      !---- Arguments ----!
      character(len=*),  dimension(:), intent(in) :: lines
      integer,           intent(in out)           :: n_ini
      integer,           intent(in)               :: n_end
      character(len=*),  intent(out)              :: spgr_hm

      !---- Local variables ----!
      character(len=1) :: csym, csym2
      integer          :: np1, np2

      !> Init
      spgr_hm=" "
      np1=n_ini
      call Read_Key_Str(lines,n_ini,n_end, "_symmetry_space_group_name_H-M",spgr_hm)

      !> TR  feb. 2015 .(re-reading the same item with another name)
      if (len_trim(spgr_hm) == 0) then
         n_ini=np1
         spgr_hm = " "
         call Read_Key_Str(lines,n_ini,n_end, "_space_group_name_H-M_alt",spgr_hm)
         if (len_trim(spgr_hm) ==0 ) spgr_hm=adjustl(lines(n_ini+1))
      end if

      if (trim(spgr_hm) =="?" .or. trim(spgr_hm)=="#") then
         spgr_hm=" "
      else
         np1=index(spgr_hm,"'")
         np2=index(spgr_hm,"'",back=.true.)
         if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
            spgr_hm=spgr_hm(np1+1:np2-1)
         else
            np1=index(spgr_hm,'"')
            np2=index(spgr_hm,'"',back=.true.)
            if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
               spgr_hm=spgr_hm(np1+1:np2-1)
            else
               spgr_hm=" "
            end if
         end if
      end if

      !> Adapting Nomenclature from ICSD to our model
      np1=len_trim(spgr_hm)
      if (np1 > 0) then
         csym=u_case(spgr_hm(np1:np1))
         select case (csym)
            case("1")
               csym2=u_case(spgr_hm(np1-1:np1-1))
               if (csym2 == "Z" .or. csym2 =="S") then
                  spgr_hm=spgr_hm(:np1-2)//":1"
               end if

            case("S","Z")
               csym2=u_case(spgr_hm(np1-1:np1-1))
               select case (csym2)
                  case ("H")
                     spgr_hm=spgr_hm(:np1-2)
                  case ("R")
                     spgr_hm=spgr_hm(:np1-2)//":R"
                  case default
                     spgr_hm=spgr_hm(:np1-1)
               end select

            case("R")
               csym2=u_case(spgr_hm(np1-1:np1-1))
               if (csym2 == "H" ) then
                  spgr_hm=spgr_hm(:np1-2)
               else
                  spgr_hm=spgr_hm(:np1-1)//":R"
               end if
            case("H")
               spgr_hm=spgr_hm(:np1-1)
               csym2=u_case(spgr_hm(np1-1:np1-1))
               if(csym2 == ":") spgr_hm=spgr_hm(:np1-2)
         end select
      end if

   End Subroutine Read_Cif_HM

   !!----
   !!---- Read_Cif_Symm
   !!----    Obtaining Symmetry Operators from Cif file
   !!----
   !!---- 27/06/2019
   !!
   Module Subroutine Read_Cif_Symm(lines,N_Ini,N_End, N_Oper, Oper_Symm)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in)     :: lines
      integer,                        intent(in out) :: n_ini
      integer,                        intent(in)     :: n_end
      integer,                        intent(out)    :: n_oper
      character(len=*), dimension(:), intent(out)    :: oper_symm

      !---- Local variables ----!
      character(len=132) :: line
      integer            :: i,np1,np2

      !> Init
      n_oper=0
      oper_symm=" "

      np1=n_ini
      call Read_Key_StrVal(lines,n_ini,n_end, "_symmetry_equiv_pos_as_xyz",line)

      !> TR  feb. 2015 .(re-reading the same item with another name)
      if (n_ini == 1) then   ! TR june 2016
         n_ini=np1
         call Read_Key_StrVal(lines,n_ini,n_end, "_space_group_symop_operation_xyz",line)
      end if

      if (len_trim(line) /=0) then
         line=adjustl(line)

         if (line(1:1) /="#" .and. line(1:1) /= "?") then
            np1=index(line,"'")
            np2=index(line,"'",back=.true.)
            if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
               n_oper=n_oper+1
               oper_symm(n_oper)=line(np1+1:np2-1)
            else
               np1=index(line,'"')
               np2=index(line,'"',back=.true.)
               if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                  n_oper=n_oper+1
                  oper_symm(n_oper)=line(np1+1:np2-1)
               end if
            end if
         end if
      end if

      do i=n_ini+1,n_end
         line=adjustl(lines(i))
         if (len_trim(line) /=0) then
            if (line(1:1) /="#" .and. line(1:1) /= "?") then
               np1=index(line,"'")
               np2=index(line,"'",back=.true.)
               if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                  n_oper=n_oper+1
                  oper_symm(n_oper)=line(np1+1:np2-1)
               else
                  np1=index(line,'"')
                  np2=index(line,'"',back=.true.)
                  if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                     n_oper=n_oper+1
                     oper_symm(n_oper)=line(np1+1:np2-1)
                  end if
               end if
            end if
         else
            n_ini=i+1
            exit
         end if
      end do

   End Subroutine Read_Cif_Symm

   !!----
   !!---- WRITE_CIF_POWDER_PROFILE
   !!----    Write a Cif Powder Profile file
   !!----
   !!---- 28/06/2019
   !!
   Module Subroutine Write_Cif_Powder_Profile(filename)
      !---- Arguments ----!
      character(len=*), intent(in) :: filename

      !---- Local Variables ----!
      logical :: info
      integer :: iunit

      !> Init
      info=.false.
      iunit=0

      !> Is open the file?
      inquire(file=trim(filename),opened=info)
      if (info) then
         inquire(file=trim(filename),number=iunit)
         close(unit=iunit)
      end if

      !> Writting
      open(newunit=iunit,file=trim(filename),status="unknown",action="write")
      rewind(unit=iunit)

      !> Head
      write(unit=iunit,fmt="(a)") "data_profile"
      write(unit=iunit,fmt="(a)") " "
      write(unit=iunit,fmt="(a)")     "_pd_block_id      ?"

      !> Profile
      write(unit=iunit,fmt="(a)") " "

      write(unit=iunit,fmt="(a)") "loop_"
      write(unit=iunit,fmt="(a)") "_pd_proc_point_id"
      write(unit=iunit,fmt="(a)") "_pd_proc_2theta_corrected             # one of "
      write(unit=iunit,fmt="(a)") "_pd_proc_energy_incident              # these "
      write(unit=iunit,fmt="(a)") "_pd_proc_d_spacing                    # three"
      write(unit=iunit,fmt="(a)") "_pd_proc_intensity_net"
      write(unit=iunit,fmt="(a)") "_pd_calc_intensity_net "
      write(unit=iunit,fmt="(a)") "_pd_proc_ls_weight      "
      write(unit=iunit,fmt="(a)") "?     ?     ?     ?     ?     ?     ?"

      write(unit=iunit,fmt="(a)") " "
      write(unit=iunit,fmt="(a)") "# The following lines are used to test the character set of files sent by     "
      write(unit=iunit,fmt="(a)") "# network email or other means. They are not part of the CIF data set.        "
      write(unit=iunit,fmt="(a)") "# abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789              "
      write(unit=iunit,fmt="(a)") "# !@#$%^&*()_+{}:"//""""//"~<>?|\-=[];'`,./ "

      close (unit=iunit)
   End Subroutine Write_Cif_Powder_Profile

   !!----
   !!---- Write_Cif_Template
   !!----    Write a Cif File
   !!----
   !!---- 28/06/2019
   !!
   Module Subroutine Write_Cif_Template(filename, Cell, SpG, At_list, Type_data, Code)
      !---- Arguments ----!
      character(len=*),        intent(in) :: filename     ! Filename
      class(Cell_Type),        intent(in) :: Cell         ! Cell parameters
      class(SpG_Type),         intent(in) :: SpG          ! Space group information
      Type (AtList_Type),      intent(in) :: At_List      ! Atoms
      integer,                 intent(in) :: Type_data    ! 0,2:Single crystal diffraction; 1:Powder
      character(len=*),        intent(in) :: Code         ! Code or name of the structure

      !---- Local Variables ----!
      logical                                 :: info, aniso
      character(len=1), parameter             :: QMARK='?'
      character(len=132)                      :: line
      character(len=30)                       :: comm,adptyp
      character(len=30),dimension(6)          :: text
      real(kind=cp)                           :: u,su, ocf
      real(kind=cp), dimension(6)             :: Ua,sua,aux
      real(kind=cp), dimension(At_List%Natoms):: occup,soccup
      integer                                 :: iunit,i, j

      !> Init
      info=.false.
      iunit=0

      !> Is this file opened?
      inquire(file=trim(filename),opened=info)
      if (info) then
         inquire(file=trim(filename),number=iunit)
         close(unit=iunit)
      end if

      !> Writing
      open(newunit=iunit, file=trim(filename),status="unknown",action="write")
      rewind(unit=iunit)

      !> Head Information
      select case (type_data)
         case (0:1)
            write(unit=iunit,fmt="(a)") "##############################################################################"
            write(unit=iunit,fmt="(a)") "###    CIF submission form for molecular structure report (Acta Cryst. C)  ###"
            write(unit=iunit,fmt="(a)") "##############################################################################"
            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "#============================================================================="
            write(unit=iunit,fmt="(a)") "data_global"
            write(unit=iunit,fmt="(a)") "#============================================================================="
            write(unit=iunit,fmt="(a)") " "

         case (2:)
            write(unit=iunit,fmt="(a)") "##################################################################"
            write(unit=iunit,fmt="(a)") "###    CIF file from CrysFML, contains only structural data    ###"
            write(unit=iunit,fmt="(a)") "##################################################################"
      end select

      !> Processing Summary
      if (type_data < 2) then
         write(unit=iunit,fmt="(a)") "# PROCESSING SUMMARY (IUCr Office Use Only)"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_journal_data_validation_number      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_journal_date_recd_electronic        ?"
         write(unit=iunit,fmt="(a)") "_journal_date_to_coeditor            ?"
         write(unit=iunit,fmt="(a)") "_journal_date_from_coeditor          ?"
         write(unit=iunit,fmt="(a)") "_journal_date_accepted               ?"
         write(unit=iunit,fmt="(a)") "_journal_date_printers_first         ?"
         write(unit=iunit,fmt="(a)") "_journal_date_printers_final         ?"
         write(unit=iunit,fmt="(a)") "_journal_date_proofs_out             ?"
         write(unit=iunit,fmt="(a)") "_journal_date_proofs_in              ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_name               ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_code               ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_notes"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_journal_techeditor_code             ?"
         write(unit=iunit,fmt="(a)") "_journal_techeditor_notes"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_journal_coden_ASTM                  ?"
         write(unit=iunit,fmt="(a)") "_journal_name_full                   ?"
         write(unit=iunit,fmt="(a)") "_journal_year                        ?"
         write(unit=iunit,fmt="(a)") "_journal_volume                      ?"
         write(unit=iunit,fmt="(a)") "_journal_issue                       ?"
         write(unit=iunit,fmt="(a)") "_journal_page_first                  ?"
         write(unit=iunit,fmt="(a)") "_journal_page_last                   ?"
         write(unit=iunit,fmt="(a)") "_journal_paper_category              ?"
         write(unit=iunit,fmt="(a)") "_journal_suppl_publ_number           ?"
         write(unit=iunit,fmt="(a)") "_journal_suppl_publ_pages            ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Submission details
         write(unit=iunit,fmt="(a)") "# 1. SUBMISSION DETAILS"
         write(unit=iunit,fmt="(a)") " "

         write(unit=iunit,fmt="(a)") "_publ_contact_author_name            ?   # Name of author for correspondence"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_address             # Address of author for correspondence"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_email           ?"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_fax             ?"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_phone           ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_contact_letter"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_requested_journal              ?"
         write(unit=iunit,fmt="(a)") "_publ_requested_coeditor_name        ?"
         write(unit=iunit,fmt="(a)") "_publ_requested_category             ?   # Acta C: one of CI/CM/CO/FI/FM/FO"

         write(unit=iunit,fmt="(a)") "#=============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Title  and Author List
         write(unit=iunit,fmt="(a)") "# 3. TITLE AND AUTHOR LIST"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_section_title"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_title_footnote"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The loop structure below should contain the names and addresses of all "
         write(unit=iunit,fmt="(a)") "# authors, in the required order of publication. Repeat as necessary."

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _publ_author_name"
         write(unit=iunit,fmt="(a)") "    _publ_author_footnote"
         write(unit=iunit,fmt="(a)") "    _publ_author_address"
         write(unit=iunit,fmt="(a)") "?                                   #<--'Last name, first name' "
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Text
         write(unit=iunit,fmt="(a)") "# 4. TEXT"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_section_synopsis"
         write(unit=iunit,fmt="(a)") ";  ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_abstract"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";          "
         write(unit=iunit,fmt="(a)") "_publ_section_comment"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_exptl_prep      # Details of the preparation of the sample(s)"
         write(unit=iunit,fmt="(a)") "                              # should be given here. "
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_exptl_refinement"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_references"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_figure_captions"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_acknowledgements"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Identifier
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") "# If more than one structure is reported, the remaining sections should be "
         write(unit=iunit,fmt="(a)") "# completed per structure. For each data set, replace the '?' in the"
         write(unit=iunit,fmt="(a)") "# data_? line below by a unique identifier."
      end if !type_data < 2
      write(unit=iunit,fmt="(a)") " "

      if (len_trim(code) == 0) then
         write(unit=iunit,fmt="(a)") "data_?"
      else
         write(unit=iunit,fmt="(a)") "data_"//code(1:len_trim(code))
      end if
      write(unit=iunit,fmt="(a)") " "

      if (type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Chemical Data
         write(unit=iunit,fmt="(a)") "# 5. CHEMICAL DATA"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_chemical_name_systematic"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_chemical_name_common             ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_moiety          ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_structural      ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_analytical      ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_iupac           ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_sum             ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_weight          ?"
         write(unit=iunit,fmt="(a)") "_chemical_melting_point           ?"
         write(unit=iunit,fmt="(a)") "_chemical_compound_source         ?       # for minerals and "
         write(unit=iunit,fmt="(a)") "                                          # natural products"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _atom_type_symbol               "
         write(unit=iunit,fmt="(a)") "    _atom_type_description          "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_dispersion_real "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_dispersion_imag "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_source          "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_length_neutron       # include if applicable"
         write(unit=iunit,fmt="(a)") "    ?    ?    ?    ?    ?      ?    "
      end if !type_data < 2
      write(unit=iunit,fmt="(a)") " "
      write(unit=iunit,fmt="(a)") "#============================================================================="
      write(unit=iunit,fmt="(a)") " "

      !> Crystal Data
      select case (type_data)
         case (0,2) ! Single Crystal or structural data only
            write(unit=iunit,fmt="(a)") "# 6. CRYSTAL DATA"
         case (1) ! Powder Data + Crystal Data
            write(unit=iunit,fmt="(a)") "# 6. POWDER SPECIMEN AND CRYSTAL DATA"
      end select
      write(unit=iunit,fmt="(a)") " "

      write(unit=iunit,fmt="(a)") "_symmetry_cell_setting               ?"
      line=SpG%SPG_Symb
      write(unit=iunit,fmt="(a)") "_symmetry_space_group_name_H-M       '"//trim(line)//"'"
      line=SpG%Hall
      write(unit=iunit,fmt="(a)") "_symmetry_space_group_name_Hall      '"//trim(line)//"'"

      write(unit=iunit,fmt="(a)") " "
      write(unit=iunit,fmt="(a)") "loop_"
      write(unit=iunit,fmt="(a)") "    _symmetry_equiv_pos_as_xyz"
      do i=1,SpG%multip
         line="'"//trim(SpG%Sym_Op(i))//"'"
         write(iunit,'(a)') trim(line)
      end do
      write(unit=iunit,fmt="(a)") " "

      do i=1,3
         text(i)=string_numstd(cell%cell(i),cell%scell(i))
         text(i+3)=string_numstd(cell%ang(i),cell%sang(i))
      end do
      write(unit=iunit,fmt='(a)')       "_cell_length_a                       "//trim(adjustl(text(1)))
      write(unit=iunit,fmt='(a)')       "_cell_length_b                       "//trim(adjustl(text(2)))
      write(unit=iunit,fmt='(a)')       "_cell_length_c                       "//trim(adjustl(text(3)))
      write(unit=iunit,fmt='(a)')       "_cell_angle_alpha                    "//trim(adjustl(text(4)))
      write(unit=iunit,fmt='(a)')       "_cell_angle_beta                     "//trim(adjustl(text(5)))
      write(unit=iunit,fmt='(a)')       "_cell_angle_gamma                    "//trim(adjustl(text(6)))
      write(unit=iunit,fmt="(a,f14.4)") "_cell_volume                   ",Cell%Vol
      if (type_data < 2) then
         write(unit=iunit,fmt="(a)") "_cell_formula_units_Z                ?"
         write(unit=iunit,fmt="(a)") "_cell_measurement_temperature        ?"
         write(unit=iunit,fmt="(a)") "_cell_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
      end if

      select case (type_data)
         case (0) ! Single Crystal
            write(unit=iunit,fmt="(a)") "_cell_measurement_reflns_used        ?"
            write(unit=iunit,fmt="(a)") "_cell_measurement_theta_min          ?"
            write(unit=iunit,fmt="(a)") "_cell_measurement_theta_max          ?"

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_exptl_crystal_description           ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_colour                ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_size_max              ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_size_mid              ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_size_min              ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_size_rad              ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_density_diffrn        ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_density_meas          ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_density_method        ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_F_000                 ?"

         case (1) ! Powder Data
            write(unit=iunit,fmt="(a)") "# The next three fields give the specimen dimensions in mm.  The equatorial"
            write(unit=iunit,fmt="(a)") "# plane contains the incident and diffracted beam."

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_spec_size_axial               ?       # perpendicular to "
            write(unit=iunit,fmt="(a)") "                                          # equatorial plane"

            write(unit=iunit,fmt="(a)") "_pd_spec_size_equat               ?       # parallel to "
            write(unit=iunit,fmt="(a)") "                                          # scattering vector"
            write(unit=iunit,fmt="(a)") "                                          # in transmission"
            write(unit=iunit,fmt="(a)") "_pd_spec_size_thick               ?       # parallel to "
            write(unit=iunit,fmt="(a)") "                                          # scattering vector"
            write(unit=iunit,fmt="(a)") "                                          # in reflection"

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# The next five fields are character fields that describe the specimen."

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_spec_mounting                         # This field should be"
            write(unit=iunit,fmt="(a)") "                                          # used to give details of the "
            write(unit=iunit,fmt="(a)") "                                          # container."
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"
            write(unit=iunit,fmt="(a)") "_pd_spec_mount_mode               ?       # options are 'reflection'"
            write(unit=iunit,fmt="(a)") "                                          # or 'transmission'"
            write(unit=iunit,fmt="(a)") "_pd_spec_shape                    ?       # options are 'cylinder' "
            write(unit=iunit,fmt="(a)") "                                          # 'flat_sheet' or 'irregular'"
            write(unit=iunit,fmt="(a)") "_pd_char_particle_morphology      ?"
            write(unit=iunit,fmt="(a)") "_pd_char_colour                   ?       # use ICDD colour descriptions"

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# The following three fields describe the preparation of the specimen."
            write(unit=iunit,fmt="(a)") "# The cooling rate is in K/min.  The pressure at which the sample was "
            write(unit=iunit,fmt="(a)") "# prepared is in kPa.  The temperature of preparation is in K.        "

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_prep_cool_rate                ?"
            write(unit=iunit,fmt="(a)") "_pd_prep_pressure                 ?"
            write(unit=iunit,fmt="(a)") "_pd_prep_temperature              ?"
      end select
      if (type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The next four fields are normally only needed for transmission experiments."
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_coefficient_mu        ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_type       ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_process_details       ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_T_min      ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_T_max      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Experimental Data
         write(unit=iunit,fmt="(a)") "# 7. EXPERIMENTAL DATA"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_exptl_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

        if (type_data == 1) then
           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "# The following item is used to identify the equipment used to record "
           write(unit=iunit,fmt="(a)") "# the powder pattern when the diffractogram was measured at a laboratory "
           write(unit=iunit,fmt="(a)") "# other than the authors' home institution, e.g. when neutron or synchrotron"
           write(unit=iunit,fmt="(a)") "# radiation is used."

           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "_pd_instr_location"
           write(unit=iunit,fmt="(a)") "; ?"
           write(unit=iunit,fmt="(a)") ";"
           write(unit=iunit,fmt="(a)") "_pd_calibration_special_details           # description of the method used"
           write(unit=iunit,fmt="(a)") "                                          # to calibrate the instrument"
           write(unit=iunit,fmt="(a)") "; ?"
           write(unit=iunit,fmt="(a)") ";"
        end if

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "_diffrn_ambient_temperature          ?"
        write(unit=iunit,fmt="(a)") "_diffrn_radiation_type               ?"
        write(unit=iunit,fmt="(a)") "_diffrn_radiation_wavelength         ?"
        write(unit=iunit,fmt="(a)") "_diffrn_radiation_source             ?"
        write(unit=iunit,fmt="(a)") "_diffrn_source                       ?"
        write(unit=iunit,fmt="(a)") "_diffrn_source_target                ?"
        write(unit=iunit,fmt="(a)") "_diffrn_source_type                  ?"

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "_diffrn_radiation_monochromator      ?"
        write(unit=iunit,fmt="(a)") "_diffrn_measurement_device_type      ?"
        write(unit=iunit,fmt="(a)") "_diffrn_measurement_method           ?"
        write(unit=iunit,fmt="(a)") "_diffrn_detector_area_resol_mean     ?   # Not in version 2.0.1"
        write(unit=iunit,fmt="(a)") "_diffrn_detector                     ?"
        write(unit=iunit,fmt="(a)") "_diffrn_detector_type                ?   # make or model of detector"
        if (type_data == 1) then
           write(unit=iunit,fmt="(a)") "_pd_meas_scan_method                 ?   # options are 'step', 'cont',"
           write(unit=iunit,fmt="(a)") "                                         # 'tof', 'fixed' or"
           write(unit=iunit,fmt="(a)") "                                         # 'disp' (= dispersive)"
           write(unit=iunit,fmt="(a)") "_pd_meas_special_details"
           write(unit=iunit,fmt="(a)") ";  ?"
           write(unit=iunit,fmt="(a)") ";"
        end if

        select case (type_data)
           case (0)
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_number                ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_av_R_equivalents      ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_av_sigmaI/netI        ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_min             ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_max             ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_full            ?"
              write(unit=iunit,fmt="(a)") "_diffrn_measured_fraction_theta_max  ?"
              write(unit=iunit,fmt="(a)") "_diffrn_measured_fraction_theta_full ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_h_min           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_h_max           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_k_min           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_k_max           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_l_min           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_l_max           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_reduction_process     ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_diffrn_standards_number             ?"
              write(unit=iunit,fmt="(a)") "_diffrn_standards_interval_count     ?"
              write(unit=iunit,fmt="(a)") "_diffrn_standards_interval_time      ?"
              write(unit=iunit,fmt="(a)") "_diffrn_standards_decay_%            ?"
              write(unit=iunit,fmt="(a)") "loop_"
              write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_h"
              write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_k"
              write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_l"
              write(unit=iunit,fmt="(a)") "?   ?   ?"

           case (1)
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "#  The following four items give details of the measured (not processed)"
              write(unit=iunit,fmt="(a)") "#  powder pattern.  Angles are in degrees."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_meas_number_of_points         ?"
              write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_min         ?"
              write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_max         ?"
              write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_inc         ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# The following three items are used for time-of-flight measurements only."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_instr_dist_src/spec           ?"
              write(unit=iunit,fmt="(a)") "_pd_instr_dist_spec/detc          ?"
              write(unit=iunit,fmt="(a)") "_pd_meas_2theta_fixed             ?"

        end select

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "#============================================================================="
        write(unit=iunit,fmt="(a)") " "

        !> Refinement Data
        write(unit=iunit,fmt="(a)") "# 8. REFINEMENT DATA"

        write(unit=iunit,fmt="(a)") " "

        write(unit=iunit,fmt="(a)") "_refine_special_details"
        write(unit=iunit,fmt="(a)") "; ?"
        write(unit=iunit,fmt="(a)") ";"

        if (type_data == 1) then
           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "# Use the next field to give any special details about the fitting of the"
           write(unit=iunit,fmt="(a)") "# powder pattern."

           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "_pd_proc_ls_special_details"
           write(unit=iunit,fmt="(a)") "; ?"
           write(unit=iunit,fmt="(a)") ";"

           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "# The next three items are given as text."
           write(unit=iunit,fmt="(a)") " "

           write(unit=iunit,fmt="(a)") "_pd_proc_ls_profile_function      ?"
           write(unit=iunit,fmt="(a)") "_pd_proc_ls_background_function   ?"
           write(unit=iunit,fmt="(a)") "_pd_proc_ls_pref_orient_corr"
           write(unit=iunit,fmt="(a)") "; ?"
           write(unit=iunit,fmt="(a)") ";"
        end if

        select case (type_data)
           case (0)
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_reflns_number_total                 ?"
              write(unit=iunit,fmt="(a)") "_reflns_number_gt                    ?"
              write(unit=iunit,fmt="(a)") "_reflns_threshold_expression         ?"

           case (1)
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_R_factor         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_wR_factor        ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_wR_expected      ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# The following four items apply to angular dispersive measurements."
              write(unit=iunit,fmt="(a)") "# 2theta minimum, maximum and increment (in degrees) are for the "
              write(unit=iunit,fmt="(a)") "# intensities used in the refinement."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_min         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_max         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_inc         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_wavelength               ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_block_diffractogram_id        ?  # The id used for the block containing"
              write(unit=iunit,fmt="(a)") "                                     # the powder pattern profile (section 11)."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# Give appropriate details in the next two text fields."
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_info_excluded_regions    ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_info_data_reduction      ?"
        end select

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "_refine_ls_structure_factor_coef     ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_matrix_type               ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_R_I_factor                ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_R_Fsqd_factor             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_R_factor_all              ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_R_factor_gt               ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_wR_factor_all             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_wR_factor_ref             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_goodness_of_fit_all       ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_goodness_of_fit_ref       ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_restrained_S_all          ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_restrained_S_obs          ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_number_reflns             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_number_parameters         ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_number_restraints         ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_number_constraints        ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_hydrogen_treatment        ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_weighting_scheme          ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_weighting_details         ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_shift/su_max              ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_shift/su_mean             ?"
        write(unit=iunit,fmt="(a)") "_refine_diff_density_max             ?"
        write(unit=iunit,fmt="(a)") "_refine_diff_density_min             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_extinction_method         ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_extinction_coef           ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_details     ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_Flack       ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_Rogers      ?"

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "# The following items are used to identify the programs used."
        write(unit=iunit,fmt="(a)") " "

        write(unit=iunit,fmt="(a)") "_computing_data_collection           ?"
        write(unit=iunit,fmt="(a)") "_computing_cell_refinement           ?"
        write(unit=iunit,fmt="(a)") "_computing_data_reduction            ?"
        write(unit=iunit,fmt="(a)") "_computing_structure_solution        ?"
        write(unit=iunit,fmt="(a)") "_computing_structure_refinement      ?"
        write(unit=iunit,fmt="(a)") "_computing_molecular_graphics        ?"
        write(unit=iunit,fmt="(a)") "_computing_publication_material      ?"
      end if  !(type_data < 2) then
      write(unit=iunit,fmt="(a)") " "
      write(unit=iunit,fmt="(a)") "#============================================================================="
      write(unit=iunit,fmt="(a)") " "

      !> Atomic Coordinates and Displacement Parameters
      write(unit=iunit,fmt="(a)") "# 9. ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS"

      write(unit=iunit,fmt="(a)") " "

      write(unit=iunit,fmt="(a)") "loop_"
      write(unit=iunit,fmt='(a)') "    _atom_site_label"
      write(unit=iunit,fmt='(a)') "    _atom_site_type_symbol"
      write(unit=iunit,fmt='(a)') "    _atom_site_fract_x"
      write(unit=iunit,fmt='(a)') "    _atom_site_fract_y"
      write(unit=iunit,fmt='(a)') "    _atom_site_fract_z"
      write(unit=iunit,fmt='(a)') "    _atom_site_U_iso_or_equiv"
      write(unit=iunit,fmt='(a)') "    _atom_site_occupancy"
      write(unit=iunit,fmt='(a)') "    _atom_site_adp_type"
      write(unit=iunit,fmt='(a)') "    _atom_site_type_symbol"

      !> Calculation of the factor corresponding to the occupation factor provided in A
      do i=1,At_List%natoms
         occup(i)=At_list%Atom(i)%occ/(real(At_List%Atom(i)%mult)/real(SpG%multip))
        soccup(i)=At_list%Atom(i)%occ_std/(real(At_List%Atom(i)%mult)/real(SpG%multip))
      end do
      ocf=sum(abs(At_list%atom(1)%x-At_list%atom(2)%x))
      if ( ocf < 0.001) then
         ocf=occup(1)+occup(2)
      else
         ocf=occup(1)
      end if
      occup=occup/ocf; soccup=soccup/ocf
      aniso=.false.
      do i=1,At_List%natoms
         line=" "
         line(2:)= At_list%Atom(i)%Lab//"  "//At_list%Atom(i)%SfacSymb

         do j=1,3
            comm=string_numstd(At_list%Atom(i)%x(j),At_list%Atom(i)%x_std(j))
            line=trim(line)//" "//trim(comm)
         end do

         comm=" "
         select case (At_List%Atom(i)%Thtype)
            case ('iso')
               adptyp='Uiso'
               select case (trim(At_List%Atom(i)%UType))
                  case ("U")
                     u=At_list%Atom(i)%U_iso
                     su=At_list%Atom(i)%U_iso_std

                  case ("B")
                     u=At_list%Atom(i)%U_iso/(8.0*pi*pi)
                     su=At_list%Atom(i)%U_iso_std/(8.0*pi*pi)

                  case ("beta")
                     u=At_list%Atom(i)%U_iso
                     su=At_list%Atom(i)%U_iso_std
               end select
               comm=string_numstd(u,su)

            case ('ani')
               aniso=.true.
               adptyp='Uani'
               select case (trim(At_List%Atom(i)%UType))
                  case ("U")
                     ua=At_List%atom(i)%u
                     sua=At_List%atom(i)%u_std

                  case ("B")
                     ua=At_List%atom(i)%u/(8.0*pi*pi)
                     sua=At_List%atom(i)%u_std/(8.0*pi*pi)

                  case ("beta")
                     aux=At_list%atom(i)%u
                     ua=get_U_from_Betas(aux,cell)
                     aux=At_list%atom(i)%u_std
                     sua=get_U_from_Betas(aux,cell)
               end select
               u=(ua(1)+ua(2)+ua(3))/3.0
               su=(ua(1)+ua(2)+ua(3))/3.0
               com=string_numstd(u,su)

            case default
               adptyp='.'
         end select
         line=trim(line)//" "//trim(comm)

         comm=string_numstd(occup(i),soccup(i))
         line=trim(line)//" "//trim(comm)
         write(unit=iunit,fmt="(a)") trim(line)//" "//trim(adptyp)//" "//At_list%atom(i)%SfacSymb
      end do

      if (aniso) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_label "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_11  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_22  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_33  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_12  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_13  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_23  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_type_symbol"

         do i=1,At_List%natoms
            if (At_List%Atom(i)%thtype /= "ani") cycle

            line=" "
            line(2:)= At_list%Atom(i)%Lab

            select case (trim(At_List%Atom(i)%UType))
               case ("U")
                  ua=At_List%atom(i)%u
                  sua=At_List%atom(i)%u_std

               case ("B")
                  ua=At_List%atom(i)%u/(8.0*pi*pi)
                  sua=At_List%atom(i)%u_std/(8.0*pi*pi)

               case ("beta")
                  aux=At_list%atom(i)%u
                  ua=get_U_from_Betas(aux,cell)
                  aux=At_list%atom(i)%u_std
                  sua=get_U_from_Betas(aux,cell)
            end select

            do j=1,6
              comm=" "
              call setnum_std(ua(j),sua(j),comm)
              line=trim(line)//" "//trim(comm)
            end do
            write(iunit,"(a)") trim(line)//"  "//A%atom(i)%SfacSymb
         end do
      end if

      if(type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# Note: if the displacement parameters were refined anisotropically"
         write(unit=iunit,fmt="(a)") "# the U matrices should be given as for single-crystal studies."

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Molecular Geometry ----!
         write(unit=iunit,fmt="(a)") "# 10. MOLECULAR GEOMETRY"

         write(unit=iunit,fmt="(a)") " "


         write(unit=iunit,fmt="(a)") "_geom_special_details                ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _geom_bond_atom_site_label_1  "
         write(unit=iunit,fmt="(a)") "    _geom_bond_atom_site_label_2  "
         write(unit=iunit,fmt="(a)") "    _geom_bond_site_symmetry_1    "
         write(unit=iunit,fmt="(a)") "    _geom_bond_site_symmetry_2    "
         write(unit=iunit,fmt="(a)") "    _geom_bond_distance           "
         write(unit=iunit,fmt="(a)") "    _geom_bond_publ_flag          "
         write(unit=iunit,fmt="(a)") "    ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _geom_contact_atom_site_label_1 "
         write(unit=iunit,fmt="(a)") "    _geom_contact_atom_site_label_2 "
         write(unit=iunit,fmt="(a)") "    _geom_contact_distance          "
         write(unit=iunit,fmt="(a)") "    _geom_contact_site_symmetry_1   "
         write(unit=iunit,fmt="(a)") "    _geom_contact_site_symmetry_2   "
         write(unit=iunit,fmt="(a)") "    _geom_contact_publ_flag         "
         write(unit=iunit,fmt="(a)") "    ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_1 "
         write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_2 "
         write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_3 "
         write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_1   "
         write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_2   "
         write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_3   "
         write(unit=iunit,fmt="(a)") "_geom_angle                   "
         write(unit=iunit,fmt="(a)") "_geom_angle_publ_flag         "
         write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_1 "
         write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_2 "
         write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_3 "
         write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_4 "
         write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_1   "
         write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_2   "
         write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_3   "
         write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_4   "
         write(unit=iunit,fmt="(a)") "_geom_torsion                   "
         write(unit=iunit,fmt="(a)") "_geom_torsion_publ_flag         "
         write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_D "
         write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_H "
         write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_A "
         write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_D   "
         write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_H   "
         write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_A   "
         write(unit=iunit,fmt="(a)") "_geom_hbond_distance_DH       "
         write(unit=iunit,fmt="(a)") "_geom_hbond_distance_HA       "
         write(unit=iunit,fmt="(a)") "_geom_hbond_distance_DA       "
         write(unit=iunit,fmt="(a)") "_geom_hbond_angle_DHA         "
         write(unit=iunit,fmt="(a)") "_geom_hbond_publ_flag         "
         write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "


         !---- Final Informations ----!
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") "# Additional structures (last six sections and associated data_? identifiers) "
         write(unit=iunit,fmt="(a)") "# may be added at this point.                                                 "
         write(unit=iunit,fmt="(a)") "#============================================================================="

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The following lines are used to test the character set of files sent by     "
         write(unit=iunit,fmt="(a)") "# network email or other means. They are not part of the CIF data set.        "
         write(unit=iunit,fmt="(a)") "# abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789              "
         write(unit=iunit,fmt="(a)") "# !@#$%^&*()_+{}:"//""""//"~<>?|\-=[];'`,./ "
      end if

      close(unit=iunit)

      return
   End Subroutine Write_Cif_Template



End SubModule IOF_003