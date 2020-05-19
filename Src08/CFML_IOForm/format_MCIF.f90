!!----
!!----
!!----
SubModule (CFML_IOForm) IO_MCIF

   !---- Local Variables ----!
   character(len=132)            :: line
   character(len=:), allocatable :: str
   integer                       :: j_ini, j_end

   Contains

   !!--++
   !!--++ Read_XTal_MCIF
   !!--++
   !!--++ Read a MCIF File
   !!--++
   !!--++ 17/05/2020
   !!
   Module Subroutine Read_XTal_MCIF(cif, Cell, Spg, AtmList, Nphase)
      !---- Arguments ----!
      type(File_Type),               intent(in)  :: cif
      class(Cell_Type),              intent(out) :: Cell
      class(SpG_Type),               intent(out) :: SpG
      Type(AtList_Type),             intent(out) :: Atmlist
      Integer,             optional, intent(in)  :: Nphase   ! Select the Phase to read

      !---- Local Variables ----!
      integer                        :: i, iph, nt_phases, it, n_ini,n_end
      integer, dimension(MAX_PHASES) :: ip

      real(kind=cp),dimension(6):: vet,vet2
      real(kind=cp)             :: val

      Type(Kvect_Info_Type)     :: Kvec

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_XTal_MCIF: No lines in the file!"
         return
      end if

      !> Calculating number of Phases
      nt_phases=0; ip=cif%nlines; ip(1)=1
      do i=1,cif%nlines
         line=adjustl(cif%line(i)%str)
         if (l_case(line(1:5)) == "data_" .and. l_case(line(1:11)) /= "data_global" )  then
            nt_phases=nt_phases+1
            ip(nt_phases)=i
         end if
      end do

      !> Read the Phase information
      iph=1
      if (present(nphase)) then
         iph=min(nphase, nt_phases)
         iph=max(1,iph)
      end if

      n_ini=ip(iph)
      n_end=ip(iph+1)

      !> Reading Cell Parameters
      call Read_Cif_Cell(cif,Cell,n_ini,n_end)
      if (Err_CFML%IErr==1) return

      !> Parent Propagation vector
      call Read_MCIF_Parent_Propagation_Vector(cif, Kvec,n_ini,n_end)

      !> Parent_Space_Group
      call Read_MCIF_Parent_SpaceG(cif,SpG,n_ini,n_end)

      !> Space_Group_Magn
      call Read_MCIF_SpaceG_Magn(cif,SpG,n_ini,n_end)
      line=SpG%BNS_symb
      line=pack_string(line)
      call Set_SpaceGroup(trim(line),"SHUBN",SpG)

      !> Atomos
      call Read_CIF_Atoms(cif, AtmList, n_ini, n_end)

      !> Moment
      call Read_MCIF_AtomSite_Moment(cif, AtmList, n_ini, n_end)


   End Subroutine Read_XTal_MCIF



   !!----
   !!---- WRITE_MCIF
   !!----
   !!----
   !!----
   !!---- 14/05/2020
   !!
   Module Subroutine Write_MCIF_Template(filename,Cell,SpG,AtmList)
      !---- Arguments ----!
      character(len=*),        intent(in) :: filename     ! Filename
      class(Cell_G_Type),      intent(in) :: Cell         ! Cell parameters
      class(SpG_Type),         intent(in) :: SpG          ! Space group information
      Type(AtList_Type),       intent(in) :: AtmList      ! Atoms

      !---- Local Variables ----!
      logical                        :: info
      integer                        :: i,j,ipr,L

      !> Init
      ipr=0

      !> Is this file opened?
      inquire(file=trim(filename),opened=info)
      if (info) then
         inquire(file=trim(filename),number=ipr)
         close(unit=ipr)
      end if

      !> Writing
      open(newunit=ipr, file=trim(filename),status="unknown",action="write")
      rewind(unit=ipr)

      !> Heading Information
      call write_cif_header(ipr)

      !> Extra items
      write(unit=ipr,fmt="(a)") " "
      write(unit=ipr,fmt="(a)") "_Neel_temperature  ?"
      write(unit=ipr,fmt="(a)") "_magn_diffrn_temperature  ?"
      write(unit=ipr,fmt="(a)") "_exptl_crystal_magnetic_properties_details"
      write(unit=ipr,fmt="(a)") ";"
      write(unit=ipr,fmt="(a)") ";"
      write(unit=ipr,fmt="(a)") "_active_magnetic_irreps_details"
      write(unit=ipr,fmt="(a)") ";"
      write(unit=ipr,fmt="(a)") ";"
      write(unit=ipr,fmt="(a)") " "

      !> Spg
      call Write_MCIF_Spg(ipr,Spg)

      !> Irrep
      write(unit=Ipr,fmt="(a)") "loop_"
      write(unit=Ipr,fmt="(a)") "_irrep_id"
      write(unit=Ipr,fmt="(a)") "_irrep_dimension"
      write(unit=Ipr,fmt="(a)") "_small_irrep_dimension"
      write(unit=Ipr,fmt="(a)") "_irrep_direction_type"
      write(unit=Ipr,fmt="(a)") "_irrep_action"
      write(unit=Ipr,fmt="(a)") "_irrep_modes_number"
      write(unit=Ipr,fmt="(a)") " ?  ?  ?  ?  ?  ?"
      write(unit=Ipr,fmt="(a)") " "

      !> Cell
      call write_cif_cell(Ipr,cell)

      !> Propag. Vectors
      select type(SpG)
         type is (SuperSpaceGroup_Type)
            if (SpG%nk > 0) then
               write(unit=Ipr,fmt="(a)") "_parent_propagation_vector.id"
               write(unit=Ipr,fmt="(a)") "_parent_propagation_vector.kxkykz"
               ! Incompleto
               !do i=1,SpG%nk
               !   line=Frac_Trans_2Dig(SpG%kv(:,i))
               !   line=adjustl(line(2:len_trim(line)-1))
               !   write(unit=Ipr,fmt="(a)") trim(SpG%kv_label(i))//"  '"//trim(line)//"'"
               !end do
               write(unit=Ipr,fmt="(a)") " "
            end if
      end select

      !> Atoms
      select case (l_case(Atmlist%atom(1)%UType))
         case ('beta')
            call write_cif_atoms(ipr,atmlist,spG,cell)
         case default
            call write_cif_atoms(ipr,atmlist,SpG)
      end select

      !> Moment
      call write_mcif_atomsite_moment(Ipr, Atmlist)

      !> close
      call write_cif_end(ipr)

   End Subroutine Write_MCIF_Template



   !!----
   !!---- WRITE_MCIF_PARENT_SPACEG
   !!----
   !!----
   !!---- 15/05/2020
   !!
   Module Subroutine Write_MCIF_Parent_SpaceG(Ipr, Spg)
      !---- Arguments ----!
      integer,          intent(in) :: Ipr
      class (Spg_type), intent(in) :: Spg

      !---- Local Variables ----!
      character(len=4) :: car

      !> Child_Transform_Pp_abc
      line=adjustl(SpG%tfrom_parent)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a,t50,a)") "_parent_space_group.child_transform_Pp_abc",trim(line)

      !> IT_number
      car="?"
      if (SpG%Parent_num > 0) write(unit=car,fmt='(i4)') SpG%Parent_num
      car=adjustl(car)
      write(unit=ipr,fmt="(a,t50,a)") "_parent_space_group.IT_number",trim(car)

      !> H-M_alt
      line=adjustl(SpG%Parent_spg)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a,t50,a)") "_parent_space_group.name_H-M_alt", trim(line)

      !> Setting
      line=adjustl(SpG%setting)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a,t50,a)") "_parent_space_group.reference_setting", trim(line)

      !> Transform_Pp_abc
      line=adjustl(SpG%mat2std)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a,t50,a)") "_parent_space_group.transform_Pp_abc", trim(line)

   End Subroutine Write_MCIF_Parent_SpaceG

   !!----
   !!---- WRITE_MCIF_SPG
   !!----
   !!----
   !!---- 15/05/2020
   !!
   Module Subroutine Write_MCIF_Spg(Ipr, Spg)
      !---- Arguments ----!
      integer,          intent(in) :: Ipr
      class (Spg_type), intent(in) :: Spg

      !---- Local Variables ----!
      integer                        :: i,j,L
      type(rational), dimension(3,3) :: unidad


      write(unit=ipr,fmt="(a)") " "

      !> setting
      if (SpG%standard_setting) then
         write(unit=Ipr,fmt="(a)") "_magnetic_space_group_standard_setting  'yes'"
      else
         write(unit=Ipr,fmt="(a)") "_magnetic_space_group_standard_setting  'no'"
      end if

      !> Parent Space group
      line=adjustl(SpG%tfrom_parent)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a)") "_parent_space_group.child_transform_Pp_abc  "//trim(line)  ! Real(4,4)
      !write(unit=ipr,fmt="(a)") "_parent_space_group.transform_Pp_abc    ?"

      line=adjustl(SpG%Parent_spg)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a)")    "_parent_space_group.name_H-M             "//trim(line)

      if (SpG%Parent_num > 0) then
         write(unit=ipr,fmt="(a,i6)")    "_parent_space_group.IT_number            ",Spg%Parent_num
      end if

      !> Space_group_Magn
      line=adjustl(SpG%BNS_symb)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a)") "_space_group_magn.name_BNS         "//trim(line)

      line=adjustl(SpG%BNS_num)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a)") "_space_group_magn.number_BNS       "//trim(line)

      line=adjustl(SpG%OG_symb)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a)") "_space_group_magn.name_OG          "//trim(line)

      line=adjustl(SpG%OG_num)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a)") "_space_group_magn.number_OG        "//trim(line)

      !> Space group SymOp Operation
      write(unit=ipr,fmt="(a)")  "loop_"
      write(unit=ipr,fmt="(a)")  "    _space_group_symop_magn_operation.id"
      write(unit=ipr,fmt="(a)")  "    _space_group_symop_magn_operation.xyz"
      L=SpG%Numops
      if (SpG%centred == 2) L=L*2
      do i=1,L
         line=trim(l_case(SpG%Symb_Op(i)))
         line="'"//trim(line)//"'"
         write(unit=ipr,fmt="(i4,a)") i,"  "//trim(line)
      end do
      write(unit=Ipr,fmt="(a)") " "

      !> Centering
      if (SpG%Num_Lat > 1 .or. SpG%Num_aLat > 0) then
         call Rational_Identity_Matrix(unidad)
         write(unit=ipr,fmt="(a)")  "loop_"
         write(unit=ipr,fmt="(a)")  "    _space_group_symop_magn_centering.id"
         write(unit=ipr,fmt="(a)")  "    _space_group_symop_magn_centering.xyz"

         j=1
         line=trim(l_case(SpG%Symb_Op(j)))
         line="'"//trim(line)//"'"
         write(unit=ipr,fmt="(i4,a)") j,"  "//trim(line)
         do i=L+1,SpG%Multip
            if (rational_equal(SpG%Op(i)%Mat,unidad)) then
               j=j+1
               line=trim(l_case(SpG%Symb_Op(j)))
               line="'"//trim(line)//"'"
               write(unit=ipr,fmt="(i4,a)") j,"  "//trim(line)
            end if
         end do
         write(unit=Ipr,fmt="(a)") " "
      end if

   End Subroutine Write_MCIF_Spg

   !!----
   !!---- WRITE_MCIF_SPACEG_SYMOP_MAGN_CENTERING
   !!----
   !!----
   !!---- 18/05/2020
   !!
   Module Subroutine Write_MCIF_SpaceG_SymOP_Magn_Centering(Ipr, Spg)
      !---- Arguments ----!
      integer,          intent(in) :: Ipr
      class (Spg_type), intent(in) :: Spg

      !---- Local Variables ----!
      integer                        :: i,j,L
      type(rational), dimension(3,3) :: unidad

      !> Init
      call Rational_Identity_Matrix(unidad)

      !> Centering
      if (SpG%Num_Lat > 1 .or. SpG%Num_aLat > 0) then
         write(unit=ipr,fmt="(a)")  "loop_"
         write(unit=ipr,fmt="(a)")  "    _space_group_symop_magn_centering.id"
         write(unit=ipr,fmt="(a)")  "    _space_group_symop_magn_centering.xyz"

         do i=1,SpG%Multip
            if (rational_equal(SpG%Op(i)%Mat(1:3,1:3),unidad)) then
               j=j+1
               line=trim(l_case(SpG%Symb_Op(i)))
               line="'"//trim(line)//"'"
               write(unit=ipr,fmt="(i4,5x,a)") j,trim(line)
            end if
         end do
         write(unit=Ipr,fmt="(a)") " "
      end if

   End Subroutine Write_MCIF_SpaceG_SymOP_Magn_Centering

   !!----
   !!---- WRITE_MCIF_SPACEG_SYMOP_MAGN_OPERATION
   !!----
   !!----
   !!---- 19/05/2020
   !!
   Module Subroutine Write_MCIF_SpaceG_SymOP_Magn(Ipr, Spg)
      !---- Arguments ----!
      integer,          intent(in) :: Ipr
      class (Spg_type), intent(in) :: Spg

      !---- Local Variables ----!
      integer                        :: i,j,L
      type(rational), dimension(3,3) :: unidad

      !> Init
      call Rational_Identity_Matrix(unidad)

      !> Operations
      write(unit=ipr,fmt="(a)")  "loop_"
      write(unit=ipr,fmt="(a)")  "    _space_group_symop_magn_operation.id"
      write(unit=ipr,fmt="(a)")  "    _space_group_symop_magn_operation.xyz"

      do i=1,SpG%Multip
         line=trim(l_case(SpG%Symb_Op(i)))
         line="'"//trim(line)//"'"
         write(unit=ipr,fmt="(i4,5x,a)") i,trim(line)
      end do
      write(unit=Ipr,fmt="(a)") " "

   End Subroutine Write_MCIF_SpaceG_SymOP_Magn

   !!----
   !!---- WRITE_MCIF_SPACEG_MAGN
   !!----
   !!----
   !!---- 18/05/2020
   !!
   Module Subroutine Write_MCIF_SpaceG_Magn(Ipr, Spg)
      !---- Arguments ----!
      integer,          intent(in) :: Ipr
      class (Spg_type), intent(in) :: Spg

      !---- Local Variables ----!
      integer :: np1,np2

      !> Name_BNS
      line=adjustl(SpG%BNS_symb)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a,t50,a)") "_space_group_magn.name_BNS",trim(line)

      !> Name_OG
      line=adjustl(SpG%OG_symb)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a,t50,a)") "_space_group_magn.name_OG",trim(line)

      !> number_BNS
      line=adjustl(SpG%BNS_num)
      if (len_trim(line) ==0) line="?"
      write(unit=ipr,fmt="(a,t50,a)") "_space_group_magn.number_BNS",trim(line)

      !> number_OG
      line=adjustl(SpG%OG_num)
      if (len_trim(line) ==0) line="?"
      write(unit=ipr,fmt="(a,t50,a)") "_space_group_magn.number_OG",trim(line)

      !> OG_wavevector_kxkykz (not implemented)

      !> point_group_name
      line=adjustl(SpG%mag_pg)
      if (len_trim(line) ==0) then
         line="?"
      else
         np1=index(line,'"')
         np2=index(line,"'")
         if (np1 > 0 .and. np2==0) then
            line="'"//trim(line)//"'"
         else if (np1 == 0 .and. np2 > 0) then
            line='"'//trim(line)//'"'
         end if
      end if
      write(unit=ipr,fmt="(a,t50,a)") "_space_group_magn.point_group_name",trim(line)

      !> point_group_number (not implemented)

      !> transform_BNS_Pp (not implemented)

      !> transform_BNS_Pp_abc

      !> transform_OG_Pp (not implemented)

      !> transform_OG_Pp_abc

      !> --------
      !> SGG Zone
      !> --------

      !> SSG_Name
      line=adjustl(SpG%SSG_symb)
      if (len_trim(line) ==0) then
         line="?"
      else
         line="'"//trim(line)//"'"
      end if
      write(unit=ipr,fmt="(a,t50,a)") "_space_group_magn.ssg_name",trim(line)

      !> SSG_Number
      line=adjustl(SpG%SSG_Nlabel)
      if (len_trim(line) ==0) line="?"
      write(unit=ipr,fmt="(a,t50,a)") "_space_group_magn.ssg_number",trim(line)


      write(unit=Ipr,fmt="(a)") " "

   End Subroutine Write_MCIF_SpaceG_Magn

   !!----
   !!---- WRITE_MCIF_SPACEG_MAGN_SSG_TRANSF
   !!----
   !!----
   !!---- 18/05/2020
   !!
   Module Subroutine Write_MCIF_SpaceG_Magn_SSG_Transf(Ipr, Spg)
      !---- Arguments ----!
      integer,          intent(in) :: Ipr
      class (Spg_type), intent(in) :: Spg

      !---- Local Variables ----!
      integer :: np1,np2



   End Subroutine Write_MCIF_SpaceG_Magn_SSG_Transf

   !!----
   !!---- READ_MCIF_ATOMSITE_MOMENT
   !!----
   !!----
   !!---- 19/05/2020
   !!
   Module Subroutine Read_MCIF_AtomSite_Moment(cif, AtmList,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),   intent(in)    :: cif
      Type(AtList_Type), intent(inout) :: AtmList
      integer, optional, intent(in)  :: i_ini,i_end   ! Index to Finish

      !---- Local Variables ----!
      logical                                      :: found, is_new
      character(len=40), dimension(:), allocatable :: dire
      integer, dimension(15)                       :: lugar
      integer, dimension(3)                        :: ivet
      integer                                      :: i, j, k,nl, np, nt,ic, iv
      integer                                      :: i1,i2,i3
      real(kind=cp), dimension(3)                  :: vet1,vet2

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_AtomSite_Moment: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop
      found=.false.
      str="_atom_site_moment."
      nl=len_trim(str)
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         !> search the loop
         do j=i-1,j_ini,-1
            line=adjustl(cif%line(j)%str)
            if (len_trim(line) <=0) cycle
            if (line(1:1) == '#') cycle

            npos=index(line,'loop_')
            if (npos ==0) cycle
            found=.true.
            j_ini=j+1
            exit
         end do
         exit
      end do
      if (.not. found) return

      lugar=0
      j=0
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle
         if (line(1:nl) /= trim(str)) exit

         select case (trim(line))
            case ('_atom_site_moment.label')
               j=j+1
               lugar(1)=j
            case ('_atom_site_moment.Cartn_x')
               j=j+1
               lugar(2)=j
            case ('_atom_site_moment.Cartn_y')
               j=j+1
               lugar(3)=j
            case ('_atom_site_moment.Cartn_z')
               j=j+1
               lugar(4)=j
            case ('_atom_site_moment.crystalaxis_x')
               j=j+1
               lugar(5)=j
            case ('_atom_site_moment.crystalaxis_y')
               j=j+1
               lugar(6)=j
            case ('_atom_site_moment.crystalaxis_z')
               j=j+1
               lugar(7)=j
            case ('_atom_site_moment.spherical_azimuthal')
               j=j+1
               lugar(8)=j
            case ('_atom_site_moment.spherical_modulus')
               j=j+1
               lugar(9)=j
            case ('_atom_site_moment.spherical_polar')
               j=j+1
               lugar(10)=j
            case ('_atom_site_moment.Cartn')
               j=j+1
               lugar(11)=j
            case ('_atom_site_moment.crystalaxis')
               j=j+1
               lugar(12)=j
            case ('_atom_site_moment.symmform')
               j=j+1
               lugar(13)=j
            case ('_atom_site_moment.modulation_flag')
               j=j+1
               lugar(14)=j
            case ('_atom_site_moment.refinement_flags_magnetic')
               j=j+1
               lugar(15)=j
         end select
      end do

      np=count(lugar > 0)
      if (np == 0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_AtomSite_Moment: Check the loop! "
         return
      end if

      if (allocated(dire)) deallocate(dire)
      allocate(dire(np))
      dire=" "

      !> Read vales
      is_new=.true.
      do i=j_ini, j_end
         line=adjustl(cif%line(i)%str)

         if (is_new) then
            if (line(1:1) == '#') cycle
            if (line(1:1) == "_" ) cycle
            if (len_trim(line) <=0) exit
            j=1
            nt=0
         end if

         call get_words(line,dire(j:),ic)
         if (ic ==0) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_AtomSite_Moment: Revise the line: "//trim(line)
            return
         end if
         nt=nt+ic
         if (nt < np) then
            j=ic+1
            is_new=.false.
            cycle
         end if

         !> Which atom
         do k=1,Atmlist%natoms
            if (l_case(trim(dire(lugar(1)))) /= l_case(trim(atmlist%atom(k)%lab))) cycle

            !> Cartn[xyz]
            if (all(lugar(2:4) > 0)) then
               i1=index(dire(lugar(2)),'(')
               i2=index(dire(lugar(3)),'(')
               i3=index(dire(lugar(4)),'(')
               if (i1==0 .and. i2==0 .and. i3==0) then
                  line=trim(dire(lugar(2)))//"  "//trim(dire(lugar(3)))//"  "//trim(dire(lugar(3)))
                  call get_num(line,vet1,ivet,iv)
                  if (iv == 3) Atmlist%atom(k)%moment=vet1
               else
                  select type(at=>Atmlist%atom)
                     class is (Atm_std_type)
                        line=trim(dire(lugar(2)))//"  "//trim(dire(lugar(3)))//"  "//trim(dire(lugar(3)))
                        call get_numstd(line,vet1,vet2,iv)
                        if (iv == 3) then
                           At(k)%moment=vet1
                           At(k)%moment_std=vet2
                        end if
                  end select
               end if

            !> Crystalaxis[xyz]
            else if (all(lugar(5:7) > 0)) then
               i1=index(dire(lugar(5)),'(')
               i2=index(dire(lugar(6)),'(')
               i3=index(dire(lugar(7)),'(')
               if (i1==0 .and. i2==0 .and. i3==0) then
                  line=trim(dire(lugar(5)))//"  "//trim(dire(lugar(6)))//"  "//trim(dire(lugar(7)))
                  call get_num(line,vet1,ivet,iv)
                  if (iv == 3) Atmlist%atom(k)%moment=vet1
               else
                  select type(at=>Atmlist%atom)
                     class is (Atm_std_type)
                        line=trim(dire(lugar(5)))//"  "//trim(dire(lugar(6)))//"  "//trim(dire(lugar(7)))
                        call get_numstd(line,vet1,vet2,iv)
                        if (iv == 3) then
                           At(k)%moment=vet1
                           At(k)%moment_std=vet2
                        end if
                  end select
               end if

            !> Spherical
            else if (all(lugar(8:10) > 0)) then
               i1=index(dire(lugar(8)) ,'(')
               i2=index(dire(lugar(9)) ,'(')
               i3=index(dire(lugar(10)),'(')
               if (i1==0 .and. i2==0 .and. i3==0) then
                  line=trim(dire(lugar(9)))//"  "//trim(dire(lugar(10)))//"  "//trim(dire(lugar(8)))
                  call get_num(line,vet1,ivet,iv)
                  if (iv == 3)  Atmlist%atom(k)%moment=vet1
               else
                  select type(at=>Atmlist%atom)
                     class is (Atm_std_type)
                        line=trim(dire(lugar(9)))//"  "//trim(dire(lugar(10)))//"  "//trim(dire(lugar(8)))
                        call get_numstd(line,vet1,vet2,iv)
                        if (iv == 3) then
                           At(k)%moment=vet1
                           At(k)%moment_std=vet2
                        end if
                  end select
               end if

            !> Cartn
            else if (lugar(11) > 0) then
               i1=index(dire(lugar(11)),'(')
               i2=index(dire(lugar(11)+1),'(')
               i3=index(dire(lugar(11)+2),'(')
               if (i1==0 .and. i2==0 .and. i3==0) then
                  line=trim(dire(lugar(11)))//"  "//trim(dire(lugar(11)+1))//"  "//trim(dire(lugar(11)+2))
                  call get_num(line, vet1,ivet,iv)
                  if (iv == 3)  Atmlist%atom(k)%moment=vet1
               else
                  select type(at=>Atmlist%atom)
                     class is (Atm_std_type)
                        line=trim(dire(lugar(11)))//"  "//trim(dire(lugar(11)+1))//"  "//trim(dire(lugar(11)+2))
                        call get_numstd(line, vet1,vet2,iv)
                        if (iv == 3) then
                           At(k)%moment=vet1
                           At(k)%moment_std=vet2
                        end if
                  end select
               end if

            !> crystalaxis
            else if (lugar(12) > 0) then
               i1=index(dire(lugar(12)),'(')
               i2=index(dire(lugar(12)+1),'(')
               i3=index(dire(lugar(12)+2),'(')
               if (i1==0 .and. i2==0 .and. i3==0) then
                  line=trim(dire(lugar(12)))//"  "//trim(dire(lugar(12)+1))//"  "//trim(dire(lugar(12)+2))
                  call get_num(line, vet1,ivet,iv)
                  if (iv == 3)  Atmlist%atom(k)%moment=vet1
               else
                  select type(at=>Atmlist%atom)
                     class is (Atm_std_type)
                        line=trim(dire(lugar(12)))//"  "//trim(dire(lugar(12)+1))//"  "//trim(dire(lugar(12)+2))
                        call get_numstd(line, vet1,vet2,iv)
                        if (iv == 3) then
                           At(k)%moment=vet1
                           At(k)%moment_std=vet2
                        end if
                  end select
               end if
            end if
            atmlist%atom(k)%mom=maxval(abs(atmlist%atom(k)%moment))

            !> Symmform
            if (lugar(13) > 0) then
            end if

            !> Modulation
            if (lugar(14) > 0) then
            end if

            !> Refinement
            if (lugar(15) > 0) then
            end if
         end do

         is_new=.true.
      end do

      write(unit=Ipr,fmt="(a)") " "

   End Subroutine Read_MCIF_AtomSite_Moment

   !!----
   !!---- WRITE_MCIF_ATOMSITE_MOMENT
   !!----
   !!----
   !!---- 15/05/2020
   !!
   Module Subroutine Write_MCIF_AtomSite_Moment(Ipr, AtmList)
      !---- Arguments ----!
      integer,           intent(in) :: Ipr
      Type(AtList_Type), intent(in) :: AtmList

      !---- Local Variables ----!
      integer                         :: i, j
      character(len=40), dimension(3) :: text

      !> Moment
      write(unit=Ipr,fmt="(a)") "loop_"
      write(unit=Ipr,fmt="(a)") "_atom_site_moment.label"
      write(unit=Ipr,fmt="(a)") "_atom_site_moment.crystalaxis_x"
      write(unit=Ipr,fmt="(a)") "_atom_site_moment.crystalaxis_y"
      write(unit=Ipr,fmt="(a)") "_atom_site_moment.crystalaxis_z"

      select type(at=>atmlist%atom)
         type is (Atm_Type)
            do i=1,Atmlist%natoms
               if (At(i)%mom < 0.001) cycle
               do j=1,3
                  text(j)=String_NumStd(At(i)%moment(j),0.0)
               end do
               write(unit=Ipr,fmt="(a8,3a12)") At(i)%lab, (text(j),j=1,3)
            end do

         class is (Atm_std_type)
            do i=1,Atmlist%natoms
               if (At(i)%mom < 0.001) cycle
               do j=1,3
                  text(j)=String_NumStd(At(i)%moment(j),At(i)%moment_std(j))
               end do
               write(unit=Ipr,fmt="(a8,3a12)") At(i)%lab,(text(j),j=1,3)
            end do
      end select
      write(unit=Ipr,fmt="(a)") " "

   End Subroutine Write_MCIF_AtomSite_Moment

   !!----
   !!---- READ_MCIF_PARENT_PROPAGATION_VECTOR
   !!----
   !!---- 17/05/2020 11:45:38
   !!
   Module Subroutine Read_MCIF_Parent_Propagation_Vector(cif, Kvec,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),       intent(in)  :: cif
      Type(Kvect_Info_Type), intent(out) :: Kvec
      integer, optional,     intent(in)  :: i_ini,i_end   ! Index to Finish

      !---- Local Variables ----!
      logical                :: found
      integer, dimension( 2) :: lugar   !   1:id, 2: kxkykz
      integer, dimension(3)  :: ivet
      integer                :: i,j,nl,iv
      integer                :: np1,np2

      real(kind=cp),     dimension(3) :: vet
      character(len=20), dimension(3) :: cvet

      !> Init
      call Allocate_KVector(0, 0, Kvec)

      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Parent_Propagation_Vector: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop for Propagation vectors
      found=.false.
      str="_parent_propagation_vector.kxkykz"
      nl=len_trim(str)
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         !> search the loop
         do j=i-1,j_ini,-1
            line=adjustl(cif%line(j)%str)
            if (len_trim(line) <=0) cycle
            if (line(1:1) == '#') cycle

            npos=index(line,'loop_')
            if (npos ==0) cycle
            j_ini=j+1
            found=.true.
            exit
         end do
         exit
      end do
      if (.not. found) return

      lugar=0
      j=0
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (index(line,'_parent')==0) exit
         if (len_trim(line) <=0) exit

         select case (trim(line))
            case ('_parent_propagation_vector.id')
               j=j+1
               lugar(1)=j
            case ('_parent_propagation_vector.kxkykz')
               j=j+1
               lugar(2)=j
         end select
      end do

      !> Number of k-vectors
      j=0
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (index(line,'_parent') > 0) cycle
         if (len_trim(line) <=0) exit
         j=j+1
      end do
      if (j ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Parent_propagation_Vector: 0 K-vectors"
         return
      end if

      call allocate_Kvector(j,0,kvec)

      !> Read k-vectors
      j=0
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (index(line,'_parent') > 0) cycle
         if (len_trim(line) <=0) exit

         !> eliminar tabs
         do
            iv=index(line,TAB)
            if (iv == 0) exit
            line(iv:iv)=' '
         end do

         !> Id
         !if (lugar(1) == 1) call Cut_String(line)

         !> vector
         np1=index(line,'[')
         np2=index(line,']')

         call get_num(line(np1+1:np2-1),vet,ivet,iv)
         if (iv == 3) then
            j=j+1
            Kvec%kv(:,j)=vet
            cycle
         end if
         call clear_error()

         call get_words(line(np1+1:np2-1),cvet,iv)
         call clear_error()
         if (iv == 3) then
            vet=0.0_cp

            vet(1)=read_fract(cvet(1))
            if (err_CFML%Ierr==1) return

            vet(2)=read_fract(cvet(2))
            if (err_CFML%Ierr==1) return

            vet(3)=read_fract(cvet(3))
            if (err_CFML%Ierr==1) return

            j=j+1
            Kvec%kv(:,j)=vet
         end if

      end do

   End Subroutine Read_MCIF_Parent_Propagation_Vector

   !!----
   !!---- WRITE_MCIF_PARENT_PROPAGATION_VECTOR
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Write_MCIF_Parent_Propagation_Vector(Ipr,KVec)
      !---- Arguments ----!
      integer, intent(in)               :: Ipr
      Type(Kvect_Info_Type), intent(in) :: Kvec

      !---- Local Variables ----!
      character(len=3) :: car
      integer          :: i,np1,np2

      !> K-vectors
      if (Kvec%nk <=0) return

      write(unit=Ipr,fmt="(a)") "loop_"
      write(unit=Ipr,fmt="(a)") "    _parent_propagation_vector.id"
      write(unit=Ipr,fmt="(a)") "    _parent_propagation_vector.kxkykz"
      do i=1,kvec%nk
         write(unit=car,fmt='(i3)') i
         car=adjustl(car)

         line=Frac_Trans_2Dig(kvec%kv(:,i))
         np1=index(line,'(')
         np2=index(line,')')
         if (np1 > 0) line(np1:np1)=" "
         if (np2 > 0) line(np2:np2)=" "
         line=adjustl(line)
         do
            np1=index(line,',')
            if (np1 ==0) exit
            line(np1:np1)=' '
         end do

         write(unit=ipr,fmt='(a, t10,a)') "   k"//trim(car), "["//trim(line)//"]"
      end do
      write(unit=ipr,fmt="(a)") " "

   End Subroutine Write_MCIF_Parent_Propagation_Vector

   !!----
   !!---- READ_MCIF_SPACEG_MAGN
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Read_MCIF_SpaceG_Magn(cif,Spg,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),       intent(in)    :: cif
      class(SpG_Type),       intent(inout) :: SpG
      integer, optional,     intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      character(len=1)            :: car
      integer                     :: i,iv,np
      integer,       dimension(1) :: ivet
      real(kind=cp), dimension(1) :: vet

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Parent_SpaceG: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Name_BNS
      str="_space_group_magn.name_BNS"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%BNS_symb=trim(adjustl(line))

         exit
      end do

      !> Name_OG
      str="_space_group_magn.name_OG"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%OG_symb=trim(adjustl(line))

         exit
      end do

      !> Number_BNS
      str="_space_group_magn.number_BNS"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%BNS_num=trim(adjustl(line))

         exit
      end do

      !> Number_OG
      str="_space_group_magn.number_OG"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%OG_num=trim(adjustl(line))

         exit
      end do

      !> OG_wavevector_kxkykz (Not implemented)

      !> point_group_name
      str="_space_group_magn.point_group_name"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         car=line(1:1)
         np=len_trim(line)

         if (car == line(np:np) .and. (car =='"' .or. car =="'")) then
            line(1:1)=" "
            line(np:np)=" "
         end if
         line=adjustl(line)
         SpG%mag_pg=trim(line)

         exit
      end do

      !> point_group_number (Not implemented)

      !> transform_BNS_Pp (Not implemented)
      !> transform_OG_Pp (Not implemented)

      !> transform_BNS_Pp_abc
      str="_space_group_magn.transform_BNS_Pp_abc"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         !SpG%XXXX=trim(adjustl(line))

         exit
      end do

      !> transform_OG_Pp_abc
      str="_space_group_magn.transform_OG_Pp_abc"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         !SpG%XXXX=trim(adjustl(line))

         exit
      end do

      !> --------
      !> SSG Zone
      !> --------

      !> SSG_Name
      str="_space_group_magn.ssg_name"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%SSG_Symb=trim(adjustl(line))

         exit
      end do

      !> SSG_number
      str="_space_group_magn.ssg_number"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%SSG_NLabel=trim(adjustl(line))

         exit
      end do

   End Subroutine Read_MCIF_SpaceG_Magn

   !!----
   !!---- READ_MCIF_SPACEG_MAGN_TRANSFORMS
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Read_MCIF_SpaceG_Magn_Transf(cif, Spg,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),       intent(in)    :: cif
      class(SpG_Type),       intent(inout) :: SpG
      integer, optional,     intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      logical                     :: found
      integer                     :: i,j,iv,np
      integer,       dimension(5) :: lugar
      integer,       dimension(1) :: ivet
      real(kind=cp), dimension(1) :: vet

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_Magn_Transf: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop
      found=.false.
      str="_space_group_magn_transforms."
      nl=len_trim(str)
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         !> eliminar tabs
         do
            iv=index(line,TAB)
            if (iv == 0) exit
            line(iv:iv)=' '
         end do

         npos=index(line,str)
         if (npos ==0) cycle

         !> search the loop
         do j=i-1,j_ini,-1
            line=adjustl(cif%line(j)%str)
            if (len_trim(line) <=0) cycle
            if (line(1:1) == '#') cycle

            npos=index(line,'loop_')
            if (npos ==0) cycle
            found=.true.
            j_ini=j+1
            exit
         end do
         exit
      end do
      if (.not. found) return

      lugar=0
      j=0
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle
         if (line(1:nl) /='_space_group_magn_transforms.') exit

         select case (trim(line))
            case ('_space_group_magn_transforms.id')
               j=j+1
               lugar(1)=j
            case ('_space_group_magn_transforms.Pp_abc')
               j=j+1
               lugar(2)=j
            case ('_space_group_magn_transforms.Pp')
               j=j+1
               lugar(3)=j
            case ('_space_group_magn_transforms.description')
               j=j+1
               lugar(4)=j
            case ('_space_group_magn_transforms.source')
               j=j+1
               lugar(5)=j
         end select
      end do

      !> reading transformations
      np=0
      do i=j_ini, j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit
         if (line(1:1) == "_" .or. line(1:5) == "loop_") exit

         !!!! Completar
      end do

   End Subroutine Read_MCIF_SpaceG_Magn_Transf

   !!----
   !!---- READ_MCIF_SPACEG_SYMOP_MAGN_CENTERING
   !!----      ???????
   !!---- 17/05/2020
   !!
   Module Subroutine Read_MCIF_SpaceG_SymOP_Magn_Centering(cif, Spg,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),       intent(in)    :: cif
      class(SpG_Type),       intent(inout) :: SpG
      integer, optional,     intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      logical                     :: found
      character(len=1)            :: sc
      character(len=40)           :: line2
      integer                     :: i,j,iv,np,nl,nl1
      integer,       dimension(3) :: lugar

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Centering: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop
      found=.false.
      str="_space_group_symop_magn_centering."
      nl=len_trim(str)
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         !> search the loop
         do j=i-1,j_ini,-1
            line=adjustl(cif%line(j)%str)
            if (len_trim(line) <=0) cycle
            if (line(1:1) == '#') cycle

            npos=index(line,'loop_')
            if (npos ==0) cycle
            found=.true.
            j_ini=j+1
            exit
         end do
         exit
      end do
      if (.not. found) return

      lugar=0
      j=0
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle
         if (line(1:nl) /='_space_group_symop_magn_centering.') exit

         select case (trim(line))
            case ('_space_group_symop_magn_centering.id')
               j=j+1
               lugar(1)=j
            case ('_space_group_symop_magn_centering.xyz')
               j=j+1
               lugar(2)=j
            case ('_space_group_symop_magn_centering.description')
               j=j+1
               lugar(3)=j
         end select
      end do

      !> reading symop for centering or anticentering traslations
      np=count(lugar > 0 )
      j=0
      do i=j_ini, j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit
         if (line(1:1) == "_") cycle

         !> first word is the ID?
         if (lugar(1) == 1) call Cut_String(line,nl1,line2)

         line=adjustl(line)

         !> comilla simple como marcador de caracter
         select case (line(1:1))
            case ("'")
               sc="'"

            case ('"')
               sc='"'

            case default
               sc=" "
         end select

         !!!! COMPLETAR
      end do

   End Subroutine Read_MCIF_SpaceG_SymOP_Magn_Centering

   !!----
   !!---- READ_MCIF_SPACEG_MAGN_TRANSFORMS
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Read_MCIF_SpaceG_Magn_SSG_Transf(cif,Spg,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),       intent(in)    :: cif
      class(SpG_Type),       intent(inout) :: SpG
      integer, optional,     intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      character(len=1)            :: car
      integer                     :: i,iv,np
      integer,       dimension(1) :: ivet
      real(kind=cp), dimension(1) :: vet

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_Magn_SSG_Transf: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

   End Subroutine Read_MCIF_SpaceG_Magn_SSG_Transf

   !!----
   !!---- READ_MCIF_PARENT_SPACE_GROUP
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Read_MCIF_Parent_SpaceG(cif,Spg,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),       intent(in)    :: cif
      class(SpG_Type),       intent(inout) :: SpG
      integer, optional,     intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      integer                     :: i,iv
      integer,       dimension(1) :: ivet
      real(kind=cp), dimension(1) :: vet

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Parent_SpaceG: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Child_Transform_Pp_abc
      str="_parent_space_group.child_transform_Pp_abc"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%tfrom_parent=trim(adjustl(line))

         exit
      end do

      !> IT_number
      str="_parent_space_group.IT_number"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         call Get_num(line,vet,ivet,iv)
         if (iv /=1) then
            err_CFML%Ierr=1
            err_CFML%Msg="_parent_space_group.IT_number is a number. Please, check it!"
            return
         end if

         Spg%Parent_num=ivet(1)
         exit
      end do

      !> H-M
      str="_parent_space_group.name_H-M_alt"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%Parent_spg=trim(adjustl(line))

         exit
      end do

      !> Reference_setting
      str="_parent_space_group.reference_setting"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%setting=trim(adjustl(line))

         exit
      end do

      !> Transorm_Pp_abc
      str="_parent_space_group.transform_Pp_abc"
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         call cut_string(line)
         line=adjustl(line)

         do
            iv=index(line,"'")
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if

            iv=index(line,'"')
            if (iv > 0) then
               line(iv:iv)=" "
               cycle
            end if
            exit
         end do
         SpG%mat2std=trim(adjustl(line))

         exit
      end do

   End Subroutine Read_MCIF_Parent_SpaceG

   !!----
   !!---- WRITE_MCIF_SPACEG_MAGN_TRANSF
   !!----
   !!----
   !!---- 18/05/2020
   !!
   Module Subroutine Write_MCIF_SpaceG_Magn_Transf(Ipr, SpG)
      !---- Arguments ----!
      integer,          intent(in) :: Ipr
      class(SpG_Type),  intent(in) :: Spg

      !---- Local Variables ----!
      integer                         :: i, j
      character(len=40), dimension(3) :: text

      !> Moment
      write(unit=Ipr,fmt="(a)") "loop_"
      write(unit=Ipr,fmt="(a)") "   _space_group_magn_transforms.id"
      write(unit=Ipr,fmt="(a)") "   _space_group_magn_transforms.Pp_abc"
      write(unit=Ipr,fmt="(a)") "   _space_group_magn_transforms.description"
      write(unit=Ipr,fmt="(a)") "   _space_group_magn_transforms.source"

      write(unit=Ipr,fmt="(a)") "   ?  ?   ?   ?"

      write(unit=Ipr,fmt="(a)") " "

   End Subroutine Write_MCIF_SpaceG_Magn_Transf


End SubModule IO_MCIF

