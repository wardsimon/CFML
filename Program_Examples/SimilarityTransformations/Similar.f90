   Program Similar
!!---------------------------------------------------------------------------------
!!--- Program SIMILAR
!!--- Purpose: Generate files with geometric transformations of atom positions.
!!---          Change of Space Group Settings
!!---          Relations group-subgroup
!!--- (C) Created by JRC at Laboratorire Leon Brillouin, January 2002. Based in an old
!!---     Fortran 77 program written in 1990 at ILL. Is was re-written based in CrysFML
!!---     Updated in January 2014 (JRC)
!!--- Author: Juan Rodriguez-Carvajal (Institut Laue-Langevin)
!!---------------------------------------------------------------------------------
    use CFML_GlobalDeps
    Use CFML_Crystallographic_Symmetry, only: Set_SpaceGroup, Write_SpaceGroup, Space_Group_Type, &
                                              get_stabilizer, Get_multip_pos, Get_SymSymb, &
                                              get_orbit,applyso, get_T_SubGroups  , &
                                              Lattice_trans, similar_transf_SG
    Use CFML_String_Utilities, only: l_case, number_lines, reading_lines, pack_string, u_case, get_transf, &
                                     ERR_String_Mess, ERR_String
    Use CFML_Crystal_Metrics
    Use CFML_Atom_TypeDef,     only:  Atom_list_Type, Allocate_Atom_list, Write_Atom_List, &
                                      Set_Atom_Equiv_List,Atom_Equiv_List_Type

    Use CFML_Geometry_Calc,    only: point_list_type, get_transf_list, allocate_point_list, &
                                     deallocate_point_list,set_orbits_inlist,Set_New_AsymUnit, &
                                     err_geom_Mess,err_geom
    Use CFML_IO_formats

    Use CFML_Math_3D, only: set_eps, determ_a

    Implicit None

    type(file_list_type) :: file_dat
    integer                                 :: nlines, nlong, n_ini,n_end

    type (Crystal_Cell_Type)                :: Cell, Cell_n
    type (Space_Group_Type)                 :: SpaceGroup,SpaceGroup_n, SubG
    type (Space_Group_Type),dimension (512) :: SubGroup
    type (Atom_list_Type)                   :: A, A_n   !List of atoms in the asymmetric unit

    type (Atom_Equiv_List_Type)             :: Ate,Ate_n  !List of all atoms in the cell

    type (Point_list_type)                  :: pl, pl_n

    character(len=1)         :: ans
    character(len=5)         :: so_ord
    character(len=16)        :: nam
    character(len=20)        :: spp, sppg,symb !symbol of space group
    character(len=80)        :: line , title, cmdline
    character(len=256)       :: filcod,outfil,texto
    integer, parameter       :: lun1=1,lun2=6,lun=2
    integer                  :: i, j, numops, ier, ln, nauas, natc, iid,  nmag, len_cmdline, &
                                lenf, lr, nsg, nat, i1,i2
    integer                  :: l,ng, indx, k, order,mulg, norbi, n, m, ifail, indice,i_cfl,mult
    real                     :: seconds, End_time, start_time, rminutes, hours, det,occ
    integer, dimension(3,3)  :: Mat        !Auxiliary matrix
    real, dimension(3,3)     :: trans      !matrix transforming the cells
    real, dimension(3  )     :: orig       !origin of the transformed cell in old cell
    real, dimension(3  )     :: xp         !auxiliary 3D-vector
    real, dimension(6  )     :: cel        !cell parameters
    logical                  :: iprin, trans_given, subg_given, index_given, fix_given, full_given, esta
    integer, dimension (48)  :: spg_ptr, ptr, norb, sbg_ptr
    real,   dimension(3,48)  :: xo !orbits matrix
    character(len=512)       :: cosets, stabil
    character(len=1)         :: fix_lat
    integer                  :: narg

      narg=COMMAND_ARGUMENT_COUNT()
      len_cmdline=0
      if(narg > 0) then
              call GET_COMMAND_ARGUMENT(1,cmdline)
              len_cmdline=len_trim(cmdline)
              outfil=" "
              if(narg > 1) then
                call GET_COMMAND_ARGUMENT(1,outfil)
                i=index(outfil,".")
                if(i /= 0) outfil=outfil(1:i-1)
              end if
      end if

      call set_eps(0.001)   !Needed for atom position comparisons

      write(unit=*,fmt="(/,/,6(a,/))")                                                 &
      "                    ------ PROGRAM SIMILAR ------",                             &
      "                    ---- Version 2.0 Jan-2002----",                             &
      "  **************************************************************************",  &
      "  * Cell and coordinate transformations, Space Group Settings, coordinates *",  &
      "  **************************************************************************",  &
      "                     (LLB-SPSMS JRC- January 2002 )"

    DO
      if(len_cmdline /=0) then
        lenf=index(cmdline,".")-1
        if(lenf <= 0) then
          filcod=cmdline
        else
          filcod=cmdline(1:lenf)
        end if
        if(len_trim(outfil) == 0) outfil=filcod
        ln=len_trim(filcod)
        lr=len_trim(outfil)
      else
        write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.cfl/cif (give xx): "
        read(unit=*,fmt="(a)") filcod
        ln=len_trim(filcod)
        write(unit=*,fmt="(3a)",advance="no") " => Code of the output files (.coo,.sym) ( <cr>= ",filcod(1:ln),") :"
        read(unit=*,fmt="(a)") outfil
        lr=len_trim(outfil)
        IF(lr == 0) THEN
          outfil=filcod
          lr=len_trim(outfil)
        END IF
      end if

     inquire(file=trim(filcod)//".cif",exist=esta)
     if(esta) then
       call Readn_set_Xtal_Structure(trim(filcod)//".cif",Cell,SpaceGroup,A,Mode="CIF")
     else
       inquire(file=trim(filcod)//".cfl",exist=esta)
       if( .not. esta) then
         write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cif (or .cfl) doesn't exist!"
         stop
       end if
       call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpaceGroup,A,Mode="CFL",file_list=file_dat)
     end if
     If(err_form) then
       write(unit=*,fmt="(a)") trim(err_form_mess)
     else
       exit
     end if
    End do

    open(unit=lun,file=outfil(1:lr)//".coo",status="replace",action="write")

    !Writing titles and content of the input file
      write(unit=lun,fmt="(/,/,6(a,/))")                                               &
      "                    ------ PROGRAM SIMILAR ------",                             &
      "                    ---- Version 2.0 Jan-2002----",                             &
      "  **************************************************************************",  &
      "  * Cell and coordinate transformations, Space Group Settings, coordinate  *",  &
      "  **************************************************************************",  &
      "                     (LLB-SPSMS JRC- January 2002 )"
    write(unit=lun,fmt="(/,a,a,/)")" => Content of the input file: ",filcod(1:ln)//".cfl"

 !  Second reading and start calculations
 !  -------------------------------------
    CALL CPU_TIME(start_time)

    n_ini=1
    n_end=file_dat%nlines
    trans_given=.false.

     !call Read_File_transf(file_dat%line,n_ini,n_end,trans,orig)
     !if(.not. err_form) then
     !  trans_given=.true.
     !else
     !  write(unit=*,fmt="(a)") " => Error: "//trim(Err_Form_Mess)
     !end if

    n_ini=1
    subg_given = .false.
    call Read_File_Spg (file_dat%line,n_ini,n_end,sppg,"sub")
     if(.not. err_form) then
       subg_given=.true.
     else
       write(unit=*,fmt="(a)") " => No explicit subgroup given "
     end if
     index_given=.false.
     indice=48
     fix_given=.false.
     full_given=.false.
     ERR_String=.false.
    do i=4,n_end
      line=l_case(adjustl(file_dat%line(i)))
      j=index(line," ")
      if(line(1:5) == "trans") then
        line=adjustl(line(j:))
        call  Get_Transf(line,trans,orig)
        if(ERR_String) then
          write(unit=*,fmt="(a)") " => Error: "//trim(Err_String_Mess)
        else
          trans_given=.true.
        end if
        cycle
      end if
      if(index(line,"index") /= 0) then
        read(unit=file_dat%line(i)(6:),fmt=*,iostat=ier) indice
        if(ier == 0) index_given=.true.
        cycle
      end if
      if(index(line,"fix_lat") /= 0) then
        fix_lat=adjustl(file_dat%line(i)(8:))
        if(index("PABCIRF",fix_lat) /= 0) fix_given=.true.
        cycle
      end if
      if(index(line,"full") /= 0) full_given=.true.
    end do

 ! Print the whole information on cell, space group and atoms

    call Write_SpaceGroup(SpaceGroup,lun,full=.true.)
    call Write_Crystal_Cell(Cell,lun)
    call Write_Atom_List(A,lun=lun)

    call Set_Atom_Equiv_List(SpaceGroup,cell,A,Ate,lun)

    do
      if(.not. trans_given) then
        write(unit=*,fmt="(/,a)") " => Give the transformation Matrix (S1...S9 numbers)"
        write(unit=*,fmt="(a)")"     "
        write(unit=*,fmt="(a)")"       A = S1 a + S2 b + S3 c"
        write(unit=*,fmt="(a)")"       B = S4 a + S5 b + S6 c"
        write(unit=*,fmt="(a)")"       C = S7 a + S8 b + S9 c"
        write(unit=*,fmt="(a)")"     "
        write(unit=*,fmt="(a)",advance="no")" => :"
        read(unit=*,fmt=*) (trans(i,:),i=1,3)
        write(unit=*,fmt="(a)",advance="no")" => Give the new origin (in old basis): "
        read(unit=*,fmt=*) orig
      end if
      write(unit=lun,fmt="(/,a,/,a,/)")" => Change of Space group setting according to the transformation:", &
                          "    --------------------------------------------------------------"
         write(unit=lun,fmt="(a)") "                Matrix M (A'= M A)              Origin"
      write(unit=*,fmt="(/,a,/,a,/)")" => Change of Space group setting according to the transformation:", &
                          "    --------------------------------------------------------------"
         write(unit=*,fmt="(a)") "                  Matrix M (A'= M A)              Origin"
      do i=1,3
         write(unit=lun,fmt="(a,(3f8.4,a,f10.5))") "            ",trans(i,1:3), "          ", orig(i)
         write(unit=*,  fmt="(a,(3f8.4,a,f10.5))") "            ",trans(i,1:3), "          ", orig(i)
      end do
      det=determ_a(trans)
      det=abs(det)
     !call change_setting_SG(SpaceGroup,trans,orig,SpaceGroup_n) !This is for non-conventional
      write(unit=lun,fmt="(/,a,/)")" => Maximal subgroup of the original space group according to the cell transformation:"
      if(fix_given) then
         call Similar_Transf_SG(trans,orig,SpaceGroup,SpaceGroup_n,fix_lat=fix_lat)
      else
         call Similar_Transf_SG(trans,orig,SpaceGroup,SpaceGroup_n)
      end if
      call Write_SpaceGroup(SpaceGroup_n,lun, .true.)
      call Change_Setting_Cell(Cell,trans,Cell_n)

      write(unit=lun,fmt="(/,a,/,a,/)")" => Change of Unit cell according to the transformation:", &
                          "    ----------------------------------------------------"
      write(unit=lun,fmt="(a)") "                Matrix M (A'= M A)               Origin"

      do i=1,3
         write(unit=lun,fmt="(a,(3f8.4,a,f10.5))") "            ",trans(i,1:3), "          ", orig(i)
      end do
      call Write_Crystal_Cell(Cell_n,lun)

      ! It is needed to increase the asymmetric unit because some atoms may be lost
      ! after decreasind the number of symmetry operators in the new subgroup.
      ! call change_setting_atoms(Cell_n,A,trans,orig,A_n)  <= not adequate in this context

      write(unit=*,fmt="(a)")" => Creating The new Asymmetric unit"
      !call Allocate_Atom_List(100,A_n)
      !call Set_new_AsymUnit(cell_n,SpaceGroup_n,Ate,trans,orig,A_n,debug="D")
      call Set_new_AsymUnit(SpaceGroup_n,Ate,trans,orig,A_n,debug="D")
      if(err_geom) then
        write(unit=*,fmt="(a)") trim(err_geom_Mess)
      end if
      ! Write the atoms the in new asymmetric unit
      call Write_Atom_List(A_n,lun=lun)

      write(unit=*,fmt="(a)")"    ..............done!"
      write(unit=lun,fmt="(/,a,/,a,/)")" => List of all atoms in the new cell", &
                          "    ---------------------------------"
      call Set_Atom_Equiv_List(SpaceGroup_n,cell_n,A_n,Ate_n,lun)


      call allocate_point_list(SpaceGroup%multip,pl,ier)  !Allocate the point list for original cell
      if(ier /= 0) then
         write(unit=*,fmt="(a)")  " => Error allocating  ", SpaceGroup%multip," points for PL"
      end if
      if(det >= 1.0) then
        nat=SpaceGroup%multip*det+1
      else
        nat=SpaceGroup%multip+1
      end if
      call Allocate_Point_List(nat,pl_n,ier)   !allocate the new point list in the new cell
      if(ier /= 0) then
         write(unit=*,fmt="(a)")  " => Error allocating  ", nat," points for PL_n"
      end if


     !Determine all translationengleiche subgroups of the new space group
      CALL CPU_TIME(seconds)
      write(unit=*,fmt="(a)")  " => Start calculation of subgroups: "
          call get_T_SubGroups(SpaceGroup_n,SubGroup,nsg)
      CALL CPU_TIME(end_time)
      write(unit=*,fmt="(a,i3,a,f9.3,a)") " =>", nsg, " Subgroups calculated !, CPU-time: ",end_time-seconds, " seconds"

      write(unit=lun,fmt="(/,/,a,a,a,/)") " => LIST of Translationengleiche Subgroups for: ",&
                                               SpaceGroup_n%Spg_Symb,SpaceGroup_n%hall
      do i=1,nsg
        indx=SpaceGroup_n%Multip/SubGroup(i)%multip
        if(index_given .and. indx > indice) cycle
        ng=SubGroup(i)%numops
        write(unit=lun,fmt="(4a,i2,30a)") " => ", SubGroup(i)%Spg_Symb, SubGroup(i)%hall,&
          " Index: [",indx,"]   ->  { ", ( trim(SubGroup(i)%SymopSymb(l))//" ; ",l=1,ng-1),&
          trim(SubGroup(i)%SymopSymb(ng))," }    ", trim(SubGroup(i)%centre)
      end do


      write(unit=lun,fmt="(///,a)") "    ---------------------------------------------------------  "
      write(unit=lun,fmt="(a  )")   "    Decomposition of orbits in Translationengleiche Subgroups  "
      write(unit=lun,fmt="(a,/)")   "    ---------------------------------------------------------  "
      !
      !  Writing the notation of the orbits with respect to the original group
      !
       write(unit=lun,fmt="(a,a)") "    Below, the atoms of the asymmetric unit of the original space group : ", &
                  trim(SpaceGroup%Spg_Symb)//"   "//trim(SpaceGroup%hall)
       write(unit=lun,fmt="(a)")   "    Are labelled as follows : "

       do n=1,A%natoms  !loop over atoms in the asymmetric unit of the original group

         Call get_orbit(A%atom(n)%x,SpaceGroup,pl%np,pl%x)
         do j=1,pl%np
           write(unit=pl%nam(j),fmt="(a,i3,a)") trim(A%atom(n)%lab)//"(",j,")"
           pl%nam(j)=pack_string(pl%nam(j))
         end do
         write(unit=lun,fmt="(a,48a)") "     Atom: "//trim(A%atom(n)%lab)//" -> ", ( trim(pl%nam(j))//"   ",j=1,pl%np)
       end do

      write(unit=lun,fmt="(/,a)")   "    The representative atoms of the different sub-orbits are obtained by appending"
      write(unit=lun,fmt="(  a)")   "    an underscore '_' followed by a figure numbering the different orbits  "
      write(unit=lun,fmt="(  a)")   "    Example: an atom like Fe1(3)_4 corresponds to the orbit 4 that is obtained from the "
      write(unit=lun,fmt="(  a)")   "             sublattice (3) of the atom Fe1 in the asymmetric unit of the original group "

      if(index_given) write(unit=lun,fmt="(/,a,i2,a,/)") "    Only sub-groups of index below or equal to [",indice,"] are output"

      CALL CPU_TIME(seconds)
      write(unit=*,fmt="(a)")  " => Start Decomposition of orbits for subgroup: "

      Mat=0
      do i=1,3
         Mat(i,i)=1 !identity matrix
      end do
      do i=1,nsg
        indx=SpaceGroup_n%Multip/SubGroup(i)%multip
        if(index_given .and. indx > indice .or. indx == 1) cycle
        write(*,"(a,i4,a,i2,a)") " => SubGroup #",i," of index [",indx,"]"
        write(unit=lun,fmt="(/,a)") "    --------------------------------------------------------------------------------------"
        write(unit=lun,fmt="(4a,i2,a)")  " => Decomposition of orbits for subgroup: ", &
                                          SubGroup(i)%Spg_Symb, SubGroup(i)%hall,"  Index [",indx,"]"
        write(unit=texto,fmt="(a,i2)") trim(SubGroup(i)%Spg_Symb)//"_"//trim(SubGroup(i)%hall)//"_ind",indx
        texto=Pack_String(texto)
        do j=1,len_trim(texto)
          if(texto(j:j) == "'" .or. texto(j:j) == '"') texto(j:j)="p"
          if(texto(j:j) == "/" ) texto(j:j)="_"
        end do
        write(unit=lun,fmt="(a,/)") "    --------------------------------------------------------------------------------------"

        open(newunit=i_cfl,file=trim(texto)//".cfl",status="replace",action="write")
        write(unit=i_cfl,fmt="(a)") "TITLE   CFL-file generated by SIMILAR: "//trim(texto)
        write(unit=i_cfl,fmt="(a)") "! The symmetry operators and the atom positions come from the transformation"
        write(unit=i_cfl,fmt="(a)") "! "//trim(SpaceGroup_n%SG_setting)//"  Applied to space group: "//trim(SpaceGroup%Spg_Symb)
        write(unit=i_cfl,fmt="(a)") "! The maximal subgroup with this setting is: "//trim(SpaceGroup_n%Spg_Symb)
        write(unit=i_cfl,fmt="(a)") "! This gives rise to a series of tranlationengleiche subgroups, the current subgroup is: "
        write(unit=i_cfl,fmt="(a,i3,a)") "! "//SubGroup(i)%Spg_Symb//"  "//SubGroup(i)%hall//"  Index [",indx,"]"
        write(unit=i_cfl,fmt="(a)") "!            a           b           c       alpha    beta   gamma"
        write(unit=i_cfl,fmt="(a,3f12.5,3f8.3)") "Cell  ",cell_n%cell,cell_n%ang
        write(unit=i_cfl,fmt="(a)") "!"
        if(index(SubGroup(i)%Spg_Symb,"Unknown") /= 0) then
           do j=2,SubGroup(i)%NumOps
             write(unit=i_cfl,fmt="(a)") "GENR  "//trim(SubGroup(i)%SymOpSymb(j))
           end do
           NumOps=SubGroup(i)%NumOps
           if(SubGroup(i)%centred /= 1) write(unit=i_cfl,fmt="(a)") "GENR  "//trim(SubGroup(i)%SymOpSymb(NumOps+1))
           if(SubGroup(i)%NumLat > 1) then
             do j=2,SubGroup(i)%NumLat
               xp=SubGroup(i)%latt_trans(:,j)
               call Get_SymSymb(Mat,xp,symb)
               write(unit=i_cfl,fmt="(a)") "GENR  "//trim(symb)
             end do
           end if
        else
           write(unit=i_cfl,fmt="(a,3f12.5,3f8.3)") "SpGR  "//trim(SubGroup(i)%Spg_Symb)
        end if
        write(unit=i_cfl,fmt="(a)") "!"
        write(unit=i_cfl,fmt="(a)") "!    Label     Sfac      x         y         z       Biso       occ    2*Spin     Q"

        !Construct for each atom in the asymmetric unit of the initial space group the whole list
        !of equivalent points
        do n=1,A%natoms  !loop over atoms in the asymmetric unit of the original group

          Call get_orbit(A%atom(n)%x,SpaceGroup,pl%np,pl%x)
          do j=1,pl%np
            write(unit=pl%nam(j),fmt="(a,i3,a)") trim(A%atom(n)%lab)//"(",j,")"
            pl%nam(j)=pack_string(pl%nam(j))
          end do

          call get_transf_list(trans,orig,pl,pl_n,ifail)  !Transform the list in SpaceGroup to the list in SPGn

          if(ifail == 0) then
             call set_orbits_inlist(SubGroup(i),pl_n)  !Determines the decomposition of orbits w.r.t. SubGroup(i)
             norbi=maxval(pl_n%p)

             if(full_given) then

               write(unit=lun,fmt="(a,i3,a)")     " => Atom number: ",n," of label: "//trim(A%atom(n)%lab)
               write(unit=lun,fmt="(a,2(i3,a))")  "    List of equivalent atoms in the subgroup: ",pl_n%np,&
                                                  " atoms in " ,norbi, " orbits"
               k=0
               do j=1,pl_n%np
                 l=index(pl_n%nam(j),"_")
                 nam=pl_n%nam(j)(1:l)
                 write(unit=nam(l+1:),fmt="(i5)") pl_n%p(j)
                 nam=pack_string(nam)
                 write(unit=lun,fmt="(a,3f9.4,i6)")  "    "//nam, pl_n%x(:,j), pl_n%p(j)
                 if(pl_n%p(j) > k ) then
                   k=k+1
                   if(k <= norbi ) then
                     i1=index(nam,"("); i2=index(nam,")")
                     if(i1 /= 0 .and. i2 /=0 ) then
                        nam(i1:i2) =" "
                     end if
                     nam=pack_string(nam)
                     Call get_orbit(pl_n%x(:,j),SubGroup(i),mult,xo)
                     occ=real(mult)/real(SubGroup(i)%Multip)

                     write(unit=i_cfl,fmt="(a15,a4,5f10.5,2f8.3)")  "Atom "//nam,A%atom(n)%SfacSymb, pl_n%x(:,j),&
                                                              A%atom(n)%Biso,occ,A%atom(n)%Moment,A%atom(n)%Charge
                   end if
                 end if
               end do

             else

               k=0
               do j=1,pl_n%np
                 if(pl_n%p(j) > k ) then
                   l=index(pl_n%nam(j),"_")
                   nam=pl_n%nam(j)(1:l)
                   write(unit=nam(l+1:),fmt="(i5)") pl_n%p(j)
                   nam=pack_string(nam)
                   write(unit=lun,fmt="(a,3f9.4,i6)")  "    "//nam, pl_n%x(:,j), pl_n%p(j)
                   i1=index(nam,"("); i2=index(nam,")")
                   if(i1 /= 0 .and. i2 /=0 ) then
                      nam(i1:i2) =" "
                   end if
                   nam=pack_string(nam)
                   Call get_orbit(pl_n%x(:,j),SubGroup(i),mult,xo)
                   occ=real(mult)/real(SubGroup(i)%Multip)
                   write(unit=i_cfl,fmt="(a15,a4,5f10.5,2f8.3)")  "Atom "//nam,A%atom(n)%SfacSymb, pl_n%x(:,j),&
                                                      A%atom(n)%Biso,occ,A%atom(n)%Moment,A%atom(n)%Charge
                   k=k+1
                   if(k > norbi) exit
                 end if
               end do
             end if
          end if

        end do !n
        write(unit=i_cfl,fmt="(a)") "!"
        write(unit=i_cfl,fmt="(a)") "transf    a,b,c;0,0,0"

        close(unit=i_cfl)
      end do
      CALL CPU_TIME(end_time)
      write(unit=*,fmt="(a,f9.3,a)") " => Decompositon calculated !, CPU-time: ",end_time-seconds, " seconds"

      if(.not. trans_given) then
        write(unit=*,fmt="(a)",advance="no") " => Another setting change?: "
        read(unit=*,fmt="(a)") ans
        if(ans == "y" .or. ans == "Y") then
          trans_given=.false.
          cycle
        end if
      end if
      exit

    end do

! Orbits with respect to a subgroup of the original space group in the same cell

    if(subg_given) then

      write(unit=lun,fmt="(/,/,a)") " ---------------------------------------------------------- "
      write(unit=lun,fmt="(a  )")   " => Decomposition of atom in orbits w.r.t. a given subgroup"
      write(unit=lun,fmt="(a,/,/)") " ---------------------------------------------------------- "
      call Set_SpaceGroup(Sppg,subG,force_hall="f_hall")

      write(unit=lun,fmt="(a,a)")   " => Information on the given subgroup: ",Sppg
      call Write_SpaceGroup(subG,lun)

       do n=1,A%natoms  !loop over atoms in the asymmetric unit

         Call get_orbit(A%atom(n)%x,SpaceGroup,pl%np,pl%x)
         do j=1,pl%np
           write(unit=pl%nam(j),fmt="(a,i3,a)") trim(A%atom(n)%lab)//"(",j,")_"
           pl%nam(j)=pack_string(pl%nam(j))
         end do

         call set_orbits_inlist(subG,pl)  !Determines the decomposition of orbits w.r.t. SubG
         norbi=maxval(pl%p)
         write(unit=lun,fmt="(a,i3,a)")     " => Atom number: ",n," of label: "//trim(A%atom(n)%lab)

         write(unit=lun,fmt="(a,2(i3,a))")  "    List of equivalent atoms in the subgroup: ",pl%np,&
                                            " atoms in " ,norbi, " orbits"
         do j=1,pl%np
           l=index(pl%nam(j),"_")
           nam=pl%nam(j)(1:l)
           write(unit=nam(l+1:),fmt="(i5)") pl%p(j)
           nam=pack_string(nam)
           write(unit=lun,fmt="(a,3f9.4,i6)")  "    "//nam, pl%x(:,j), pl%p(j)
         end do

       end do !n


    end if
    CALL CPU_TIME(end_time)
    write(unit=*,fmt="(/a,f9.3,a)")  " => Program finished normally, CPU-time: ",end_time-start_time, " seconds"
    write(unit=lun,fmt="(/a,f9.3,a)")  " => Program finished normally, CPU-time: ",end_time-start_time, " seconds"
    write(unit=*,fmt="(a,a)") " => Results in file: ",  trim(outfil)//".coo"
    stop
   End Program Similar
