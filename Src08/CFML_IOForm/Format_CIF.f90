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
   Module Subroutine Read_Cif_ChemicalName(lines,N_ini,N_End,ChemName)
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

   End Subroutine Read_Cif_ChemicalName
   
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
    
   
   
End SubModule IOF_003   