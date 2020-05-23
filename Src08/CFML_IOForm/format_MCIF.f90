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
   !!--++ IS_SSG_STRUCT
   !!--++
   !!--++ Get the number of elements there are in a loop_
   !!--++
   !!--++ 21/05/2020
   !!
   Module Function Is_SSG_Struct(cif, i_ini,i_end) Result(ok)
      !---- Arguments ----!
      type(File_Type),   intent(in) :: cif
      integer, optional, intent(in) :: i_ini, i_end
      logical                       :: ok

      !---- local Variables ----!
      integer :: i

      !> Init
      ok=.false.
      if (cif%nlines <=0) return

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop
      do i=j_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,'_ssg_')
         if (npos ==0) cycle

         ok=.true.
         exit
      end do

   End Function Is_SSG_Struct

   !!--++
   !!--++ GET_NELEM_LOOP
   !!--++
   !!--++ Get the number of elements there are in a loop_
   !!--++
   !!--++ 17/05/2020
   !!
   Module Function Get_NElem_Loop(cif, keyword, i_ini,i_end) Result(N)
      !---- Arguments ----!
      type(File_Type),   intent(in) :: cif
      character(len=*),  intent(in) :: keyword
      integer, optional, intent(in) :: i_ini, i_end
      integer                       :: n

      !---- local Variables ----!
      logical :: found
      integer :: i,j,k_ini,k_end

      !> Init
      N=0
      if (cif%nlines <=0) return
      if (len_trim(keyword) <=0) return

      k_ini=1; k_end=cif%nlines
      if (present(i_ini)) k_ini=i_ini
      if (present(i_end)) k_end=i_end

      !> Search loop
      found=.false.
      str=trim(keyword)
      nl=len_trim(str)
      do i=k_ini,k_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         !> search the loop
         do j=i-1,k_ini,-1
            line=adjustl(cif%line(j)%str)
            if (len_trim(line) <=0) cycle
            if (line(1:1) == '#') cycle

            npos=index(line,'loop_')
            if (npos ==0) cycle
            k_ini=j+1
            found=.true.
            exit
         end do
         exit
      end do
      if (.not. found) return

      j=0
      do i=k_ini,k_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle
         if (line(1:nl) /=str) exit

         j=j+1
      end do

      k_ini=k_ini+j
      do i=k_ini, k_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit
         n=n+1
      end do

   End Function Get_NElem_Loop

   !!--++
   !!--++ Read_XTal_MCIF
   !!--++
   !!--++ Read a MCIF File
   !!--++
   !!--++ 17/05/2020
   !!
   Module Subroutine Read_XTal_MCIF(cif, Cell, Spg, AtmList, Kvec, Nphase)
      !---- Arguments ----!
      type(File_Type),                 intent(in)  :: cif
      class(Cell_Type),                intent(out) :: Cell
      class(SpG_Type),                 intent(out) :: SpG
      Type(AtList_Type),               intent(out) :: Atmlist
      Type(Kvect_Info_Type), optional, intent(out) :: Kvec
      Integer,               optional, intent(in)  :: Nphase   ! Select the Phase to read

      !---- Local Variables ----!
      logical                                      :: SSG
      character(len=60), dimension(:), allocatable :: symop
      integer                                      :: nsym,ncen
      integer                                      :: i, iph, nt_phases, it, n_ini,n_end
      integer, dimension(MAX_PHASES)               :: ip


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
      if (Err_CFML%IErr==1) then
         call error_message(err_CFML%Msg)
         return
      end if

      !> SSG Space group?
      ssg=Is_SSG_struct(cif,n_ini_n_end)

      !> Space group determination
      if (.not. ssg) Then
         ncen=get_Nelem_Loop(cif,'_space_group_symop_magn_centering.')
         nsym=get_Nelem_Loop(cif,'_space_group_symop_magn_operation.')
      else
         ncen=get_Nelem_Loop(cif,'_space_group_symop_magn_ssg_centering.')
         nsym=get_Nelem_Loop(cif,'_space_group_symop_magn_ssg_operation.')
      end if

      if (ncen+nsym > 0) then
         !> Allocating generators
         if (allocated(symop)) deallocate(symop)
         allocate(symop(ncen+nsym))
         symop=" "

         if (.not. ssg) then
            call Read_MCIF_SpaceG_SymOP_Magn_Centering(cif, ncen, symop,n_ini,n_end)
            call Read_MCIF_SpaceG_SymOP_Magn_Operation(cif, nsym, symop(ncen+1:),n_ini,n_end)
         else
            call Read_MCIF_SpaceG_SymOP_Magn_Ssg_Centering(cif, ncen, symop,n_ini,n_end)
            call Read_MCIF_SpaceG_SymOP_Magn_Ssg_Operation(cif, nsym, symop(ncen+1:),n_ini,n_end)
         end if
         call Set_SpaceGroup("  ",SpG,ncen+nsym,symop)
      end if
      if (Err_CFML%IErr==1) then
         call error_message(err_CFML%Msg)
         return
      end if

      !> Parent Propagation vector / Cell_wave_vector
      if (.not. ssg) then
         if (present(Kvec)) call Read_MCIF_Parent_Propagation_Vector(cif, Kvec,n_ini,n_end)
      else
         select type (Spg)
            type is (SuperSpaceGroup_Type)
               call Read_MCIF_Cell_Wave_Vector(cif, SpG,i_ini=n_ini, i_end=n_end)
               if (Err_CFML%IErr==1) return

               call Read_MCIF_AtomSite_Fourier_Wave_Vector(cif, SpG,i_ini=n_ini, i_end=n_end)
               if (Err_CFML%IErr==1) return

         end select
         if (present(Kvec)) then
            call Read_MCIF_Cell_Wave_Vector(cif, Kvec=Kvec,i_ini=n_ini, i_end=n_end)
            if (Err_CFML%IErr==1) return
            call Read_MCIF_AtomSite_Fourier_Wave_Vector(cif,  Kvec=Kvec,i_ini=n_ini, i_end=n_end)
            if (Err_CFML%IErr==1) return
         end if
      end if

      !> Atomos
      call Read_CIF_Atoms(cif, AtmList, n_ini, n_end)
      do i=1,Atmlist%natoms
         if (Atmlist%atom(i)%mult ==0) Atmlist%atom(i)%mult=Get_Multip_Pos(Atmlist%atom(i)%x,SpG)
      end do

      !> Moment
      call Read_MCIF_AtomSite_Moment(cif, AtmList, n_ini, n_end)
      if (Err_CFML%IErr==1) return

      !> Moment-Fourier
      call Read_MCIF_AtomSite_Moment_Fourier(cif, AtmList,n_ini,n_end)
      if (Err_CFML%IErr==1) then
         call error_message(err_CFML%Msg)
         return
      end if


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

      !> Spg
      call Write_MCIF_Spg(ipr,Spg)

      !> Propag. Vectors
      !call Write_mcif_parent_propagation_Vector(Ipr,Kvect)


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
      call write_mcif_spaceg_symop_magn_operation(Ipr,SpG)

      !> Centering
      call write_mcif_spaceg_symop_magn_centering(Ipr,SpG)

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
   Module Subroutine Write_MCIF_SpaceG_SymOP_Magn_Operation(Ipr, Spg)
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

      do i=1,SpG%Numops
         line=trim(l_case(SpG%Symb_Op(i)))
         line="'"//trim(line)//"'"
         write(unit=ipr,fmt="(i4,5x,a)") i,trim(line)
      end do
      write(unit=Ipr,fmt="(a)") " "

   End Subroutine Write_MCIF_SpaceG_SymOP_Magn_Operation

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
         if (line(1:nl) /= str) exit

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

      if (all(lugar(2:4) > 0)) then
         atmlist%mcomp="Cartesian"
      else if (all(lugar(5:7) > 0)) then
         atmlist%mcomp="Crystal"
      else if (all(lugar(8:10) > 0)) then
         atmlist%mcomp="Spherical"
      else if (lugar(11) > 0) then
         atmlist%mcomp="Cartesian"
      else if (lugar(12) > 0) then
         atmlist%mcomp="Crystal"
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
            atmlist%atom(k)%magnetic=.true.

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
   !!---- READ_MCIF_ATOMSITE_MOMENT_FOURIER
   !!----
   !!----
   !!---- 23/05/2020
   !!
   Module Subroutine Read_MCIF_AtomSite_Moment_Fourier(cif, AtmList,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),   intent(in)    :: cif
      Type(AtList_Type), intent(inout) :: AtmList
      integer, optional, intent(in)    :: i_ini,i_end   ! Index to Finish

      !---- Local Variables ----!
      logical                                      :: found
      character(len=40), dimension(:), allocatable :: dire
      integer, dimension(15)                       :: lugar
      integer, dimension(3)                        :: ivet
      integer                                      :: i, j, k,nl, np, nq, ic, iv, k_ini
      integer                                      :: i1,i2,i3
      real(kind=cp), dimension(3)                  :: vet1,vet2
      real(kind=cp), dimension(2)                  :: xv,xv_std

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop
      found=.false.
      str="_atom_site_moment_Fourier"
      nl=len_trim(str)

      k_ini=j_ini
      do i=k_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         !> search the loop
         do j=i-1,k_ini,-1
            line=adjustl(cif%line(j)%str)
            if (len_trim(line) <=0) cycle
            if (line(1:1) == '#') cycle

            npos=index(line,'loop_')
            if (npos ==0) cycle
            found=.true.
            k_ini=j+1
            exit
         end do
         exit
      end do
      if (.not. found) return

      lugar=0
      j=0
      do i=k_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle
         if (line(1:nl) /= str) exit

         select case (trim(line))
            case ('_atom_site_moment_Fourier.atom_site_label')
               j=j+1
               lugar(1)=j
            case ('_atom_site_moment_Fourier.axis')
               j=j+1
               lugar(2)=j
            case ('_atom_site_moment_Fourier.id')
               j=j+1
               lugar(3)=j
            case ('_atom_site_moment_Fourier.wave_vector_seq_id')
               j=j+1
               lugar(4)=j
            case ('_atom_site_moment_Fourier_param.id')
               j=j+1
               lugar(5)=j
            case ('_atom_site_moment_Fourier_param.cos')
               j=j+1
               lugar(6)=j
            case ('_atom_site_moment_Fourier_param.sin')
               j=j+1
               lugar(7)=j
            case ('_atom_site_moment_Fourier_param.modulus')
               j=j+1
               lugar(8)=j
            case ('_atom_site_moment_Fourier_param.phase')
               j=j+1
               lugar(9)=j
            case ('_atom_site_moment_Fourier_param.cos_symmform')
               j=j+1
               lugar(10)=j
            case ('_atom_site_moment_Fourier_param.sin_symmform')
               j=j+1
               lugar(11)=j
            case ('_atom_site_moment_Fourier_param.modulus_symmform')
               j=j+1
               lugar(12)=j
            case ('_atom_site_moment_Fourier_param.phase_symmform')
               j=j+1
               lugar(13)=j
         end select
      end do

      np=count(lugar > 0)
      if (np == 0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Check the loop! "
         return
      end if

      !> Check
      if (lugar(1) ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: No atom label was identified! "
         return
      end if

      if (lugar(2) ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: No specified the coordinate system!"
         return
      end if

      if (lugar(4) ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: No specified the identification for wave vector!"
         return
      end if

      if (lugar(6) > 0) then
         if (lugar(7) ==0 .or. lugar(8)> 0 .or. lugar(9)> 0) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Using only cosine/sine component of the magnetic modulation!"
            return
         end if
      end if
      if (lugar(7) > 0) then
         if (lugar(6) ==0 .or. lugar(8)> 0 .or. lugar(9)> 0) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Using only cosine/sine component of the magnetic modulation!"
            return
         end if
      end if
      if (lugar(8) > 0) then
         if (lugar(9) ==0 .or. lugar(6)> 0 .or. lugar(7)> 0) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Using only modulus/phase component of the magnetic modulation!"
            return
         end if
      end if
      if (lugar(9) > 0) then
         if (lugar(8) ==0 .or. lugar(6)> 0 .or. lugar(7)> 0) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Using only modulus/phase component of the magnetic modulation!"
            return
         end if
      end if

      if (allocated(dire)) deallocate(dire)
      allocate(dire(np))
      dire=" "

      !> Read vales
      k_ini=k_ini+np
      do i=k_ini, j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit

         call get_words(line,dire,ic)
         if (ic ==0) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Revise the line: "//trim(line)
            return
         end if

         !> Which atom
         do k=1,Atmlist%natoms
            if (l_case(trim(dire(lugar(1)))) /= l_case(trim(atmlist%atom(k)%lab))) cycle

            !> wave vector id
            call get_num(dire(lugar(4)),vet1,ivet,iv)
            if (iv /=1) then
               err_CFML%Ierr=1
               err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: ID for Wave vector identification!"
               return
            end if
            nq=ivet(1)

            if (lugar(6) > 0) then
               !> Cos/Sin
               call get_numstd(dire(lugar(6)),vet1,vet2,iv)
               if (iv /=1) then
                  err_CFML%Ierr=1
                  err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Error in cosine component!"
                  return
               end if
               xv(1)=vet1(1)
               xv_std(1)=vet2(1)

               call get_numstd(dire(lugar(7)),vet1,vet2,iv)
               if (iv /=1) then
                  err_CFML%Ierr=1
                  err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Error in sine component!"
                  return
               end if
               xv(2)=vet1(1)
               xv_std(2)=vet2(1)

               !> axis system
               select case (trim(dire(lugar(2))))
                  case ('Cx','x')
                     i1=1
                     i2=4
                  case ('Cy','y')
                     i1=2
                     i2=5
                  case ('Cz','z')
                     i1=3
                     i2=6
                  case ('mod')
                  case ('pol')
                  case ('azi')
                  case ('a1')
                  case ('a2')
                  case ('a3')
                  case default
                     err_CFML%Ierr=1
                     err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Axis system is not defined!"
                     return
               end select

            else if (lugar(8) > 0) then
               !>Modulus/Phase
               call get_numstd(dire(lugar(8)),vet1,vet2,iv)
               if (iv /=1) then
                  err_CFML%Ierr=1
                  err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Error in modulus component!"
                  return
               end if
               xv(1)=vet1(1)
               xv_std(1)=vet2(1)

               call get_numstd(dire(lugar(9)),vet1,vet2,iv)
               if (iv /=1) then
                  err_CFML%Ierr=1
                  err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Error in phase component!"
                  return
               end if
               xv(2)=vet1(1)
               xv_std(2)=vet2(1)

               !> axis system
               select case (trim(dire(lugar(2))))
                  case ('Cx','x')
                  case ('Cy','y')
                  case ('Cz','z')
                  case ('mod')
                  case ('pol')
                  case ('azi')
                  case ('a1')
                  case ('a2')
                  case ('a3')
                  case default
                     err_CFML%Ierr=1
                     err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Axis system is not defined!"
                     return
               end select

            else
               err_CFML%Ierr=1
               err_CFML%Msg="Read_MCIF_AtomSite_Moment_Fourier: Error in components of Modulation vector!"
               return
            end if

            select type(At => Atmlist%atom)
               class is (MAtm_Std_Type)
                  At(k)%Mcs(i1,nq)=xv(1)
                  At(k)%Mcs(i2,nq)=xv(2)
                  At(k)%Mcs_std(i1,nq)=xv_std(1)
                  At(k)%Mcs_std(i2,nq)=xv_std(2)
                  At(k)%Magnetic=.true.

                  At(k)%n_mc=max(At(k)%n_mc,nq)   !!!! Es correcto?????
                  At(k)%pmc_q=1     !!!!! CAMBIAR
            end select
            exit
         end do   ! Over atomList

      end do

      write(unit=Ipr,fmt="(a)") " "

   End Subroutine Read_MCIF_AtomSite_Moment_Fourier

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
      write(unit=Ipr,fmt="(a)") "    _atom_site_moment.label"
      select case (trim(l_case(atmlist%mcomp)))
         case ("cartesian")
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.Cartn_x"
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.Cartn_y"
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.Cartn_z"

         case ("crystal")
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.crystalaxis_x"
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.crystalaxis_y"
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.crystalaxis_z"

         case ("spherical")
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.spherical_modulus"
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.spherical_polar"
            write(unit=Ipr,fmt="(a)") "    _atom_site_moment.spherical_azimuthal"
      end select

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
   !!---- READ_MCIF_CELL_WAVE_VECTOR
   !!----
   !!---- Load Wave vector / Modulation vector (nk)
   !!----
   !!---- 22/05/2020
   !!
   Module Subroutine Read_MCIF_Cell_Wave_Vector(cif, SpG, Kvec, i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),                 intent(in)    :: cif
      class(SpG_Type),       optional, intent(inout) :: SpG
      Type(Kvect_Info_Type), optional, intent(out)   :: Kvec
      integer,               optional, intent(in)    :: i_ini,i_end   ! Index to Finish

      !---- Local Variables ----!
      logical                         :: found
      integer, dimension(4)           :: lugar   !   1:id, 2: x, 3: y, 4:z
      character(len=40), dimension(4) :: dire
      integer, dimension(3)           :: ivet
      integer                         :: i,j,k,nl,iv,np,k_ini
      real(kind=cp),     dimension(3) :: vet,vet2,xv
      Type(Kvect_Info_Type)           :: Kv

      !> Init
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop for Wave vectors
      found=.false.
      str="_cell_wave_vector_"
      nl=len_trim(str)
      k_ini=j_ini
      do i=k_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         !> search the loop
         do j=i-1,k_ini,-1
            line=adjustl(cif%line(j)%str)
            if (len_trim(line) <=0) cycle
            if (line(1:1) == '#') cycle

            npos=index(line,'loop_')
            if (npos ==0) cycle
            k_ini=j+1
            found=.true.
            exit
         end do
         exit
      end do
      if (.not. found) return

      lugar=0
      j=0
      do i=k_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (index(line,str)==0) exit

         select case (trim(line))
            case ('_cell_wave_vector_seq_id')
               j=j+1
               lugar(1)=j
            case ('_cell_wave_vector_x')
               j=j+1
               lugar(2)=j
            case ('_cell_wave_vector_y')
               j=j+1
               lugar(3)=j
            case ('_cell_wave_vector_z')
               j=j+1
               lugar(4)=j
         end select
      end do

      np=count(lugar > 0)
      if (np ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: Error in loop"
         return
      end if

      !> Number of modulation vectors
      k=get_NElem_Loop(cif,str,j_ini,j_end)
      if (k ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: No modulation vectors was determined!"
         return
      end if

      call allocate_Kvector(k,0,kv)

      !> Read k-vectors
      k_ini=k_ini+np

      do i=k_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit

         call get_words(line, dire,ic)

         if (ic /= np) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: Error reading the values in loop"
            return
         end if

         !> Id
         call get_num(dire(lugar(1)),vet,ivet,iv)
         if (iv /=1) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: Error reading the ID value"
            return
         end if
         if (ivet(1) == 0 .or. ivet(1) > k) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: Error reading the ID value"
            return
         end if
         j=ivet(1)

         !> x
         call get_numstd(dire(lugar(2)),vet,vet2,iv)
         if (iv /=1) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: Error reading the X value"
            return
         end if
         xv(1)=vet(1)

         !> y
         call get_numstd(dire(lugar(3)),vet,vet2,iv)
         if (iv /=1) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: Error reading the Y value"
            return
         end if
         xv(2)=vet(1)

         !> z
         call get_numstd(dire(lugar(4)),vet,vet2,iv)
         if (iv /=1) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Cell_Wave_Vector: Error reading the Z value"
            return
         end if
         xv(3)=vet(1)

         Kv%kv(:,j)=xv
      end do

      if (present(SpG)) then
         select type (SpG)
            type is (SuperSpaceGroup_Type)
               Spg%nk=Kv%nk
               if (allocated(Spg%kv)) deallocate(Spg%kv)
               if (allocated(Spg%sintlim)) deallocate(Spg%sintlim)
               if (allocated(Spg%nharm)) deallocate(Spg%nharm)
               allocate(Spg%kv(3,Spg%nk))
               allocate(Spg%sintlim(Spg%nk))
               allocate(Spg%nharm(Spg%nk))
               Spg%kv=Kv%kv
               Spg%sintlim=1.0_cp
               Spg%nharm=0
         end select
      end if

      if (present(Kvec)) then
         call allocate_Kvector(Kv%nk,0,kvec)
         Kvec%kv=Kv%kv
      end if

      call allocate_Kvector(0,0,kv)

   End Subroutine Read_MCIF_Cell_Wave_Vector

   !!----
   !!---- READ_MCIF_ATOM_SITE_FOURIER_WAVE_VECTOR
   !!----
   !!---- Load Q_coeff. Nk vectors should be ready
   !!----
   !!---- 23/05/2020
   !!
   Module Subroutine Read_MCIF_AtomSite_Fourier_Wave_Vector(cif, SpG, Kvec, i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),                 intent(in)    :: cif
      class(SpG_Type),       optional, intent(inout) :: SpG
      Type(Kvect_Info_Type), optional, intent(inout) :: Kvec
      integer,               optional, intent(in)    :: i_ini,i_end   ! Index to Finish

      !---- Local Variables ----!
      logical                          :: found
      integer, dimension(8)            :: lugar
      character(len=40), dimension(10) :: dire
      integer, dimension(10)           :: ivet
      integer                          :: i,j,k,nq,nl,nk,iv,np, k_ini
      integer                          :: k1,k2
      real(kind=cp), dimension(10)     :: vet

      Type(Kvect_Info_Type)            :: Kv

      !> Init
      nk=0;nq=0

      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop for Wave vectors
      found=.false.
      str="_atom_site_Fourier_wave_vector"
      nl=len_trim(str)
      k_ini=j_ini
      do i=k_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (len_trim(line) <=0) cycle
         if (line(1:1) == '#') cycle

         npos=index(line,str)
         if (npos ==0) cycle

         !> search the loop
         do j=i-1,k_ini,-1
            line=adjustl(cif%line(j)%str)
            if (len_trim(line) <=0) cycle
            if (line(1:1) == '#') cycle

            npos=index(line,'loop_')
            if (npos ==0) cycle
            k_ini=j+1
            found=.true.
            exit
         end do
         exit
      end do
      if (.not. found) return

      lugar=0
      j=0
      do i=k_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (index(line,str)==0) exit

         select case (trim(line))
            case ('_atom_site_Fourier_wave_vector.seq_id','_atom_site_Fourier_wave_vector_seq_id')
               j=j+1
               lugar(1)=j
            case ('_atom_site_Fourier_wave_vector.x','_atom_site_Fourier_wave_vector_x')
               j=j+1
               lugar(2)=j
            case ('_atom_site_Fourier_wave_vector.y','_atom_site_Fourier_wave_vector_y')
               j=j+1
               lugar(3)=j
            case ('_atom_site_Fourier_wave_vector.z','_atom_site_Fourier_wave_vector_z')
               j=j+1
               lugar(4)=j
            case ('_atom_site_Fourier_wave_vector.q1_coeff','_atom_site_Fourier_wave_vector_q1_coeff')
               j=j+1
               lugar(5)=j
            case ('_atom_site_Fourier_wave_vector.q2_coeff','_atom_site_Fourier_wave_vector_q2_coeff')
               j=j+1
               lugar(6)=j
            case ('_atom_site_Fourier_wave_vector.q3_coeff','_atom_site_Fourier_wave_vector_q3_coeff')
               j=j+1
               lugar(7)=j
            case ('_atom_site_Fourier_wave_vector.q_coeff','_atom_site_Fourier_wave_vector_q_coeff')
               j=j+1
               lugar(8)=j
         end select
      end do

      np=count(lugar > 0)
      if (np ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error in loop"
         return
      end if

      !> Number of Q Coeff
      Nq=get_NElem_Loop(cif,str,j_ini,j_end)
      if (Nq ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error in loop"
         return
      end if

      !> Copy info
      if (present(SpG) .and. (.not. present(Kvec))) then
         select type(Spg)
            type is (SuperSpaceGroup_Type)
               nk=Spg%nk
               if (nk ==0) then
                  err_CFML%Ierr=1
                  err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error in number of modulation vectors"
                  return
               end if
               call allocate_Kvector(nk,Nq,kv)
               Kv%kv=Spg%kv
               kv%sintlim=Spg%sintlim
               kv%nharm=Spg%nharm
         end select

      else if (present(kvec) .and. (.not. present(SpG))) then
         nk=Kvec%nk
         if (nk ==0) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error in number of modulation vectors"
            return
         end if
         call allocate_Kvector(nk,Nq,kv)
         Kv%kv=Kvec%kv
         kv%sintlim=Kvec%sintlim
         kv%nharm=Kvec%nharm
      else
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error Copying modulation vectors"
         return
      end if

      !> Checks
      if (lugar(1) ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: No ID"
         return
      end if

      !> Read Q_coeff
      k_ini=k_ini+np

      do i=k_ini,j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit

         call get_words(line, dire,ic)

         !> Id
         call get_num(dire(lugar(1)),vet,ivet,iv)
         if (iv /=1) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error reading the ID value"
            return
         end if
         if (ivet(1) == 0 ) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error reading the ID value"
            return
         end if
         j=ivet(1)

         if (lugar(8) > 0) then
            k1=index(line,'[')
            k2=index(line,']')
            if (k1 ==0 .or. k2 ==0) then
               err_CFML%Ierr=1
               err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error reading the Q_Coeff"
               return
            end if
            call get_num(line(k1+1:k2-1),vet,ivet,iv)
            if (iv > 0 .and. iv == nk) Kv%q_coeff(:,j)=ivet(1:iv)
         else
            do k=1,nk
               if (lugar(4+k)==0) then
                  err_CFML%Ierr=1
                  err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error reading the Qi_Coeff"
                  return
               end if

               call get_num(dire(lugar(4+k)),vet,ivet,iv)
               if (iv /= 1) then
                  err_CFML%Ierr=1
                  err_CFML%Msg="Read_MCIF_Atom_Site_Fourier_Wave_Vector: Error reading the Qi_Coeff"
                  return
               end if
               Kv%q_coeff(k,j)=ivet(1)
            end do
         end if

      end do

      if (present(SpG)) then
         select type (SpG)
            type is (SuperSpaceGroup_Type)
                  if (allocated(Spg%q_coeff)) deallocate(Spg%q_coeff)
                  allocate(Spg%q_coeff(nk,nq))
                  Spg%nq=nq
                  SpG%q_coeff=kv%q_coeff

         end select
      end if

      if (present(Kvec)) then
         if (allocated(Kvec%q_coeff)) deallocate(Kvec%q_coeff)
         allocate(Kvec%q_coeff(nk,nq))
         Kvec%nq=nq
         Kvec%q_coeff=kv%q_coeff
      end if

      call allocate_Kvector(0,0,kv)

   End Subroutine Read_MCIF_AtomSite_Fourier_Wave_Vector

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
   !!---- WRITE_MCIF_PARENT_PROPAGATION_VECTOR
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Write_MCIF_Cell_Wave_Vector(Ipr,SpG,KVec)
      !---- Arguments ----!
      integer, intent(in)                         :: Ipr
      class(SpG_Type),       optional, intent(in) :: Spg
      Type(Kvect_Info_Type), optional, intent(in) :: Kvec

      !---- Local Variables ----!
      integer               :: i
      Type(Kvect_Info_Type) :: K

      !> K-vectors
      if (present(Spg).and. (.not. present(Kvec))) then
         select type (SpG)
            type is (SuperSpaceGroup_Type)
               call allocate_Kvector(SpG%nk,0,k)
               K%kv=SpG%kv
         end select

      else if (present(Kvec) .and. (.not. present(SpG))) then
         call allocate_Kvector(kvec%nk,0,k)
         K%kv=Kvec%kv
      else
         return
      end if
      if (K%nk <=0) return

      write(unit=Ipr,fmt="(a)") "loop_"
      write(unit=Ipr,fmt="(a)") "    _cell_wave_vector_seq_id"
      write(unit=Ipr,fmt="(a)") "    _cell_wave_vector_x"
      write(unit=Ipr,fmt="(a)") "    _cell_wave_vector_y"
      write(unit=Ipr,fmt="(a)") "    _cell_wave_vector_z"
      do i=1,k%nk
         write(unit=Ipr,fmt='(5x,i4, 3f12.6)') i, K%kv(:,i)
      end do
      write(unit=ipr,fmt="(a)") " "

      call allocate_Kvector(0,0,k)

   End Subroutine Write_MCIF_Cell_Wave_Vector

   !!----
   !!---- WRITE_MCIF_ATOM_SITE_FOURIER_WAVE_VECTOR
   !!----
   !!---- 23/05/2020
   !!
   Module Subroutine Write_MCIF_AtomSite_Fourier_Wave_Vector(Ipr,SpG,KVec)
      !---- Arguments ----!
      integer, intent(in)                         :: Ipr
      class(SpG_Type),       optional, intent(in) :: Spg
      Type(Kvect_Info_Type), optional, intent(in) :: Kvec

      !---- Local Variables ----!
      integer                     :: i,j
      real(kind=cp), dimension(3) :: xv
      Type(Kvect_Info_Type)       :: K

      !> K-vectors
      if (present(Spg).and. (.not. present(Kvec))) then
         select type (SpG)
            type is (SuperSpaceGroup_Type)
               call allocate_Kvector(SpG%nk,Spg%nq,k)
               K%kv=SpG%kv
               K%sintlim=Spg%sintlim
               K%nharm=Spg%nharm
               K%q_coeff=Spg%q_coeff
         end select

      else if (present(Kvec) .and. (.not. present(SpG))) then
         call allocate_Kvector(kvec%nk,Kvec%nq,k)
         K%kv=Kvec%kv
         K%sintlim=Kvec%sintlim
         K%nharm=Kvec%nharm
         K%q_coeff=Kvec%q_coeff
      else
         return
      end if

      !> Check
      if (K%nk ==0 .or. K%nq ==0) return

      write(unit=Ipr,fmt="(a)") "loop_"
      write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.seq_id"
      write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.x"
      write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.y"
      write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.z"
      select case (K%nk)
         case (1)
            write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.q1_coeff"
         case (2)
            write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.q1_coeff"
            write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.q2_coeff"
         case (3)
            write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.q1_coeff"
            write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.q2_coeff"
            write(unit=Ipr,fmt="(a)") "    _atom_site_Fourier_wave_vector.q3_coeff"
      end select
      do i=1,k%nq
         xv=0.0_cp
         do j=1,k%nk
            xv=xv+ k%q_coeff(j,i)*k%kv(:,j)
         end do
         write(unit=line, fmt='(3i4)') k%q_coeff(:,i)
         write(unit=Ipr,fmt='(10x,i4, 3f12.6,3x,a)') i, xv, trim(line)
      end do
      write(unit=ipr,fmt="(a)") " "

      call allocate_Kvector(0,0,k)

   End Subroutine Write_MCIF_AtomSite_Fourier_Wave_Vector

   !!----
   !!---- WRITE_MCIF_PARENT_PROPAGATION_VECTOR
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Write_MCIF_Cell_Modulation_Dimension(Ipr,SpG)
      !---- Arguments ----!
      integer,         intent(in) :: Ipr
      class(SpG_Type), intent(in) :: Spg

      !---- Local Variables ----!
      integer               :: i
      Type(Kvect_Info_Type) :: K

      !> d
      write(unit=Ipr,fmt="(a,5x,i3)") "_cell_modulation_dimension", SpG%D-4
      write(unit=ipr,fmt="(a)") " "

   End Subroutine Write_MCIF_Cell_Modulation_Dimension

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
         SpG%mat2std_shu=trim(adjustl(line))

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
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Read_MCIF_SpaceG_SymOP_Magn_Centering(cif, nsym, symop,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),                intent(in)    :: cif
      integer,                        intent(out)   :: nsym
      character(len=*), dimension(:), intent(out)   :: symop
      integer, optional,              intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      logical                                      :: found
      character(len=80), dimension(3)              :: dire
      integer                                      :: i,j,iv,np,nl,nl1,ic
      integer,       dimension(3)                  :: lugar

      !> Init
      nsym=0
      symop=" "
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
         if (line(1:nl) /=str) exit

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
      if (np ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Centering: No items for traslations"
         return
      end if

      !> How many operators
      j_ini=j_ini+j

      do i=j_ini, j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit

         call get_words(line,dire,ic)
         if (ic /= np) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Centering: Problems reading symmetry operators"
            return
         end if

         line=adjustl(dire(lugar(2)))
         if (index(line,"'") > 0) then
            do
               iv=index(line,"'")
               if (iv ==0) exit
               line(iv:iv)=" "
            end do
         end if
         if (index(line,'"') > 0) then
            do
               iv=index(line,'"')
               if (iv ==0) exit
               line(iv:iv)=" "
            end do
         end if

         nsym=nsym+1
         symop(nsym)=trim(adjustl(line))
      end do

   End Subroutine Read_MCIF_SpaceG_SymOP_Magn_Centering

   !!----
   !!---- READ_MCIF_SPACEG_SYMOP_MAGN_SSG_CENTERING
   !!----
   !!---- 17/05/2020
   !!
   Module Subroutine Read_MCIF_SpaceG_SymOP_Magn_Ssg_Centering(cif, nsym, symop,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),                intent(in)    :: cif
      integer,                        intent(out)   :: nsym
      character(len=*), dimension(:), intent(out)   :: symop
      integer, optional,              intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      logical                                      :: found
      character(len=80), dimension(3)              :: dire
      integer                                      :: i,j,iv,np,nl,nl1,ic
      integer,       dimension(2)                  :: lugar

      !> Init
      nsym=0
      symop=" "

      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Ssg_Centering: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop
      found=.false.
      str="_space_group_symop_magn_ssg_centering."
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
         if (line(1:nl) /=str) exit

         select case (trim(line))
            case ('_space_group_symop_magn_ssg_centering.id')
               j=j+1
               lugar(1)=j
            case ('_space_group_symop_magn_ssg_centering.algebraic')
               j=j+1
               lugar(2)=j
         end select
      end do

      !> reading symop for centering or anticentering traslations
      np=count(lugar > 0 )
      if (np ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Ssg_Centering: No items for traslations"
         return
      end if

      j_ini=j_ini+j
      do i=j_ini, j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit

         call get_words(line,dire,ic)
         if (ic /= np) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Ssg_Centering: Problems reading symmetry operators"
            return
         end if

         line=adjustl(dire(lugar(2)))
         if (index(line,"'") > 0) then
            do
               iv=index(line,"'")
               if (iv ==0) exit
               line(iv:iv)=" "
            end do
         end if
         if (index(line,'"') > 0) then
            do
               iv=index(line,'"')
               if (iv ==0) exit
               line(iv:iv)=" "
            end do
         end if

         nsym=nsym+1
         symop(nsym)=trim(adjustl(line))
      end do

   End Subroutine Read_MCIF_SpaceG_SymOP_Magn_Ssg_Centering

   !!----
   !!---- READ_MCIF_SPACEG_SYMOP_MAGN_OPERATION
   !!----
   !!---- 21/05/2020
   !!
   Module Subroutine Read_MCIF_SpaceG_SymOP_Magn_Operation(cif, nsym, symop,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),                intent(in)    :: cif
      integer,                        intent(out)   :: nsym
      character(len=*), dimension(:), intent(out)   :: symop
      integer, optional,              intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      logical                                      :: found
      character(len=80), dimension(3)              :: dire
      integer                                      :: i,j,iv,np,nl,nl1,ic
      integer,       dimension(3)                  :: lugar

      !> Init
      nsym=0
      symop=" "
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Operation: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop
      found=.false.
      str="_space_group_symop_magn_operation."
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
         if (line(1:nl) /=str) exit

         select case (trim(line))
            case ('_space_group_symop_magn_operation.id')
               j=j+1
               lugar(1)=j
            case ('_space_group_symop_magn_operation.xyz')
               j=j+1
               lugar(2)=j
            case ('_space_group_symop_magn_operation.description')
               j=j+1
               lugar(3)=j
         end select
      end do

      !> reading symop for centering or anticentering traslations
      np=count(lugar > 0 )
      if (np ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Operation: No items for traslations"
         return
      end if

      j_ini=j_ini+j
      do i=j_ini, j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit

         call get_words(line,dire,ic)
         if (ic /= np) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Operations: Problems reading symmetry operators"
            return
         end if

         line=adjustl(dire(lugar(2)))
         if (index(line,"'") > 0) then
            do
               iv=index(line,"'")
               if (iv ==0) exit
               line(iv:iv)=" "
            end do
         end if
         if (index(line,'"') > 0) then
            do
               iv=index(line,'"')
               if (iv ==0) exit
               line(iv:iv)=" "
            end do
         end if

         nsym=nsym+1
         symop(nsym)=trim(adjustl(line))
      end do

   End Subroutine Read_MCIF_SpaceG_SymOP_Magn_Operation

   !!----
   !!---- READ_MCIF_SPACEG_SYMOP_MAGN_SSG_OPERATION
   !!----
   !!---- 21/05/2020
   !!
   Module Subroutine Read_MCIF_SpaceG_SymOP_Magn_Ssg_Operation(cif, nsym, symop,i_ini,i_end)
      !---- Arguments ----!
      Type(File_Type),                intent(in)    :: cif
      integer,                        intent(out)   :: nsym
      character(len=*), dimension(:), intent(out)   :: symop
      integer, optional,              intent(in)    :: i_ini,i_end

      !---- Local Variables ----!
      logical                                      :: found
      character(len=80), dimension(3)              :: dire
      integer                                      :: i,j,iv,np,nl,nl1,ic
      integer,       dimension(3)                  :: lugar

      !> Init
      nsym=0
      symop=" "
      call clear_error()
      if (cif%nlines <=0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Ssg_Operation: 0 lines "
         return
      end if

      j_ini=1; j_end=cif%nlines
      if (present(i_ini)) j_ini=i_ini
      if (present(i_end)) j_end=i_end

      !> Search loop
      found=.false.
      str="_space_group_symop_magn_ssg_operation."
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
         if (line(1:nl) /=str) exit

         select case (trim(line))
            case ('_space_group_symop_magn_ssg_operation.id')
               j=j+1
               lugar(1)=j
            case ('_space_group_symop_magn_ssg_operation.algebraic')
               j=j+1
               lugar(2)=j
         end select
      end do

      !> reading symop for centering or anticentering traslations
      np=count(lugar > 0 )
      if (np ==0) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Ssg_Operation: No items for traslations"
         return
      end if

      j_ini=j_ini+j
      do i=j_ini, j_end
         line=adjustl(cif%line(i)%str)

         if (line(1:1) == '#') cycle
         if (len_trim(line) <=0) exit

         call get_words(line,dire,ic)
         if (ic /= np) then
            err_CFML%Ierr=1
            err_CFML%Msg="Read_MCIF_SpaceG_SymOP_Magn_Ssg_Operations: Problems reading symmetry operators"
            return
         end if

         line=adjustl(dire(lugar(2)))
         if (index(line,"'") > 0) then
            do
               iv=index(line,"'")
               if (iv ==0) exit
               line(iv:iv)=" "
            end do
         end if
         if (index(line,'"') > 0) then
            do
               iv=index(line,'"')
               if (iv ==0) exit
               line(iv:iv)=" "
            end do
         end if

         nsym=nsym+1
         symop(nsym)=trim(adjustl(line))
      end do

   End Subroutine Read_MCIF_SpaceG_SymOP_Magn_Ssg_Operation
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

