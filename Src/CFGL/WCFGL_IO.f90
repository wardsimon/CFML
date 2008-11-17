module WCFGL_IO
!------------------------------------------------
! Written by Laurent C.Chapon
! July 2004.
! Updated : October 2004.
! IO for Winteracter
!------------------------------------------------
  use WCFGL_metrix
  use WCFGL_atomic_table
  use WCFGL_atom_tree
  use WCFGL_matom_tree
  use WCFGL_bond_tree
  use WCFGL_poly_tree
  use string_utilities, only: u_case
  use Math_gen, only: sind,cosd
  use WCFGL_display
  use WCFGL_quaternion
  use WCFGL_objects_definition, only: arrow_pos

  implicit none

     Logical, public          :: error_IO = .false.
     character(len=80),public :: mess_error_IO=" "

  contains
!------------------------------------------------------------------------------
  subroutine read_fst_file(filename,qview)
    character(len=*), intent(in)         :: filename
    type(quat), optional,intent(out)     :: qview
    ! Local variables
    type(quat)                           :: rotx, roty, rotz
    integer                              :: ier, is_box, is_cell, is_space, is_atom, is_matom,&
                                            is_scale, is_k, is_lattice, is_view, is_spher, is_rotxyz,  &
                                            is_rotax, is_xsym, is_msym, is_color, is_radius, no_line,&
                                            is_bond, is_multiple, is_nodisplay,is_poly,&
                                            is_genr,is_conn, is_skp, is_bkg, num_g,is_group,is_edge_color,&
                                            num_k,num_xsym,num_msym,kchoice,msymchoice,num_skp,mstart, &
                                            is_edges, is_envelop, is_envelop_color,is_molecule,&
                                            i

    logical                              :: ierror, dead, multiple,g_begin,mphase_begin,&
                                            k_begin,msym_begin,matom_begin,group,edges,matom_envelop,matom_edges,matom_envcolor
    character(len=256)                   :: line,upline
    character(len=80)                    :: spacegr
    character(len=2)                     :: symbol, symbol1, symbol2
    character(len=10)                    :: label, label1, label2
    character(len=1)                     :: lattice
    real, dimension(3)                   :: pos, a, angle, axis, edgecolor,pos1,pos2
    real, dimension(3,24)                :: all_k=0.0
    real,dimension(4)                    :: color, background_color, envelop_color
    real                                 :: xmin, xmax, norm, theta, phi, angx,angy,angz, &
                                            ymin, ymax, zmin, zmax, scal,  &
                                            radius,minvalue, maxvalue,skreal_a, skreal_b,&
                                            skreal_c, skim_a, skim_b, skim_c, phase
    complex,dimension(3,24)              :: my_skj
    real   ,dimension(24)                :: my_phikj
    character(len=80), dimension(192)    :: xsym , genr
    character(len=80), dimension(192,24) :: all_msym
    logical,           dimension(24)     :: envelop=.false.
    real, dimension(4,24)                :: my_color
    integer, dimension(24)               :: my_point2k, my_point2msym
    !-------Initialise parameters
    no_line=0
    num_g=0
    num_k       =0
    num_xsym    =0
    num_msym    =0
    mstart=0
    scal  =1.0
    arrow_pos  = -0.5
    ierror        = .false.
    dead          = .false.
    g_begin       = .true.
    mphase_begin  = .false.
    k_begin       = .true.
    msym_begin    = .true.
    group         = .false.
    envelop       = .false.
    matom_envcolor= .false.
    molecule      = .false.
    maxval_poly = 1.9
    pos1=2.0
    pos2=-2.0
    call define_spacegroup("P 1")
    !-------End initialise parameters


    open(unit=1,file=trim(filename),status="old",iostat=ier)

    if (ier/=0) then
      ierror =.true.
    else
      call init_box()
      call empty_atom_list()
      call empty_bond_list()
      call empty_matom_list()
      call empty_poly_list()
      call erase_magnetic_phase()
      if(present(qview))  qview=quat(1.0,(/0.0,0.0,0.0/))
      do
        read (unit=1,fmt="(a)", iostat=ier) line
        ! READING ERROR
        if (ier /= 0) then
          close(unit=1)
          exit
        ! BEGIN
        else
          ! READ COMMENT
          if (index(line(1:1),"!")/=0 .or. index(line(1:1),"#")/=0) cycle
          upline=u_case(line)
          no_line=no_line+1

          ! DETECT KEYWORDS
          is_box       = index(upline,"BOX")
          is_cell      = index(upline,"CELL")
          is_space     = index(upline,"SPACEG")
          is_atom      = index(upline,"ATOM")
          is_matom     = index(upline,"MATOM")
          is_scale     = index(upline,"SCALE")
          is_radius    = index(upline,"RADIUS")
          is_color     = index(upline,"COLOR")
          is_K         = index(upline,"K ")
          is_lattice   = index(upline,"LATTICE")
          is_xsym      = index(upline,"SYMM")
          is_msym      = index(upline,"MSYM")
          is_nodisplay = index(upline,"NODISPLAY")
          is_bond      = index(upline,"BOND")
          is_conn      = index(upline,"CONN")
          is_multiple  = index(upline,"MULTIPLE")
          is_genr      = index(upline,"GENR")
          is_skp       = index(upline,"SKP")
          is_bkg       = index(upline,"BKG")
          is_group     = index(upline,"GROUP")
          is_view      = index(upline,"VIEW")
          is_spher     = index(upline,"SPHER")
          is_rotxyz    = index(upline,"ROTXYZ")
          is_rotax     = index(upline,"ROTAX")
          is_poly      = index(upline,"POLY")
          is_edges     = index(upline,"EDGE")
          is_edge_color= index(upline,"EDGECOL")
          is_envelop   = index(upline,"ENVELOP")
          is_molecule  = index(upline,"MOLECULE")
          if(is_molecule /= 0) molecule=.true.
          is_envelop_color = index(upline,"ENVELOPCOL")

          if(index(upline,"ARROW_DISP") /= 0) arrow_pos=0.0
          ! READ OPTIONAL KEYWORDS

          if (is_group /= 0) group=.true.

          dead=(is_nodisplay/=0)
          multiple=(is_multiple/=0)

          if (is_color /= 0) then
            read(unit=line(is_color+5:),fmt=*,iostat=ier) color
            if (ier /= 0) then
              ierror=.true.
           end if
          end if

          if (is_radius /= 0) then
            read(unit=line(is_radius+6:),fmt=*,iostat=ier) radius
            if (ier /= 0) then
              ierror=.true.
            end if
          end if

          if (is_scale /= 0 .and. mphase_begin) then
            read(unit=line(is_scale+5:),fmt=*,iostat=ier) scal
            if(ier /= 0) then
              ierror=.true.
            end if
          end if




          if(present(qview)) then

              if (is_view /= 0 ) then
                read(unit=line(is_view+4:),fmt=*,iostat=ier) axis
                if (ier /= 0) then
                  ierror=.true.
                else
                  norm=.norm. axis
                  if(norm > 0.00001) then
                   axis=axis/norm
                   !Rotation angle
                   phi=acosd(axis(3))
                   axis=(/-axis(2),axis(1),0.0/)
                   norm=.norm. axis
                   if(norm > 0.00001) then
                     axis=axis/norm
                     qview=quat(cosd(0.5*phi),axis*sind(0.5*phi))
                   end if
                  end if
                end if
              end if

              if (is_rotxyz /= 0) then
                read(unit=line(is_rotxyz+6:),fmt=*,iostat=ier) angx,angy,angz
                if (ier /= 0) then
                  ierror=.true.
                else
                  rotx=quat(cosd(angx*0.5),(/sind(angx*0.5),0.0,0.0/))
                  roty=quat(cosd(angy*0.5),(/0.0,sind(angy*0.5),0.0/))
                  rotz=quat(cosd(angz*0.5),(/0.0,0.0,sind(angz*0.5)/))
                  qview=rotz*roty*rotx
                end if
              end if

              if (is_rotax /= 0) then
                read(unit=line(is_rotax+5:),fmt=*,iostat=ier) phi,axis
                if (ier /= 0) then
                  ierror=.true.
                else
                  qview=q_from_rot(phi,axis)
                end if
              end if

             if (is_spher /= 0) then
                read(unit=line(is_spher+5:),fmt=*,iostat=ier) theta,phi
                if (ier /= 0) then
                  ierror=.true.
                else
                  axis=(/-sind(phi)*sind(theta),cosd(phi)*sind(theta),0.0/)
                  norm=.norm. axis
                  if(norm > 0.00001) then
                   axis=axis/norm
                   qview=quat(cosd(0.5*theta),axis*sind(0.5*theta))
                  end if
                end if
              end if

          end if !present(view)

          ! END OF OPTIONAL KEYWORDS

          ! READ SPACE GROUP
          !-----If space group symbol is given
          if (is_space /= 0) then
           read(unit=line(is_space+6:),fmt="(a)") spacegr
           call define_spacegroup(spacegr)
          end if
          !-----If generators are given
          if (is_genr /= 0 .and. g_begin) then
            num_g=num_g+1
            read(unit=line(is_genr+4:),fmt="(a)") genr(num_g)
            do  !repeat reading until continuous GENR lines are exhausted
              read (unit=1,fmt="(a)", iostat=ier) line
              upline=u_case(line)
              no_line=no_line+1
              is_genr = index(upline,"GENR")
              if( is_genr /= 0) then
                 num_g=num_g+1
                 read(unit=line(is_genr+4:),fmt="(a)") genr(num_g)
              else
                no_line=no_line-1
                backspace(unit=1)
                g_begin=.false.
                exit
              end if
            end do
            call define_spacegroup("  ",genr,num_g)
            cycle
          end if
          ! END OF READ SPACE GROUP

          ! BEGIN MAGNETIC PHASE
          if (index(line,"{")/=0) then
            mphase_begin=.true.
          end if
          ! END MAGNETIC PHASE
          if (index(line,"}")/=0) then
            mphase_begin=.false.
            k_begin     =.true.
            num_xsym    =0
            num_k       =0
            num_msym    =0
          end if
          ! READ BACKGROUND COLOR
          if (is_bkg /=0) then
            read(unit=line(is_bkg+3:),fmt=*,iostat=ier) background_color
            if (ier ==0) then
              background_color=background_color
            end if
          end if
          ! READ BOX DIMENSION
          if (is_box /= 0) then
            read(unit=line(is_box+3:),fmt=*,iostat=ier) xmin,xmax,ymin,ymax,zmin,zmax
            if (ier /= 0) then
              call init_box()
            else
              call define_box(xmin,xmax,ymin,ymax,zmin,zmax)
            end if
          end if
          ! READ CELL PARAMETERS
          if (is_cell /= 0) then
            read(unit=line(is_cell+4:),fmt=*,iostat=ier) a, angle
            if (ier /= 0) then
              call define_cell((/10.0,10.0,10.0/),(/90.0,90.0,90.0/))
            else
              call define_cell(a,angle)
              current_cell%dead=dead
              current_cell%multiple=multiple
            end if
          end if
          ! READ ATOMS
          if (is_atom /= 0 .and. (.not.(mphase_begin))) then
            read(unit=line(is_atom+4:),fmt=*,iostat=ier) label, symbol, pos
            if (ier /= 0) then
              ierror=.true.
            else
              do i=1,3
                if(pos(i) < pos1(i)) pos1(i) = pos(i)
                if(pos(i) > pos2(i)) pos2(i) = pos(i)
              end do
              if (is_radius == 0) radius=get_radius_from_symbol(symbol)
              if (is_color  == 0) then
                call push_atom(label=label, symbol=symbol, xf=pos, Biso=0.5, radius=radius,dead=dead)
              else
                call push_atom(label, symbol, pos, 0.5, radius,dead,color)
              end if
            end if
          end if
          ! READ BONDS BETWEEN LABELS
          if (is_bond /= 0) then
            read(unit=line(is_bond+4:),fmt=*,iostat=ier) label1, label2, minvalue, maxvalue
            if (ier /= 0) then
              ierror=.true.
            else
              if (is_radius == 0) radius=1.0
              if (is_color == 0) then
                call push_bond_by_label(label1,label2,minvalue,maxvalue,radius,dead)
              else
                call push_bond_by_label(label1,label2,minvalue,maxvalue,radius,dead,color)
              end if
            end if
          end if
          if (is_poly /= 0) then
          read(unit=line(is_poly+4:),fmt=*,iostat=ier) label1
            if (ier /= 0) then
              ierror=.true.
            else
            if (is_radius /= 0)  read(unit=line(is_radius+5:),fmt=*,iostat=ier) radius
              if (is_edges == 0) then
                edges=.false.
              else
                edges=.true.
              end if
                if (is_edge_color /=0) then
                  read(unit=line(is_edge_color+7:),fmt=*,iostat=ier) edgecolor
                else
                  edgecolor=(/0.0,0.0,0.0/)
              end if
            if (is_color ==0) then
              call push_poly(label1,dead=dead,show_edges=edges,edge_color=edgecolor,edge_radius=radius)
            else
              call push_poly(label1,dead=dead,color=color,show_edges=edges,edge_color=edgecolor,edge_radius=radius)
            end if
          end if
          end if

          ! READ BONDS BETWEEN SYMBOL
          if (is_conn /= 0) then
            read(unit=line(is_conn+4:),fmt=*,iostat=ier) symbol1, symbol2, minvalue, maxvalue
            if(maxvalue > maxval_poly) maxval_poly = maxvalue
            if (ier /= 0) then
              ierror=.true.
            else
              if (is_radius == 0) radius=1.0
              if (is_color == 0) then
                call push_bond_by_symbol(symbol1,symbol2,minvalue,maxvalue,radius,dead)
              else
                call push_bond_by_symbol(symbol1,symbol2,minvalue,maxvalue,radius,dead, color)
              end if
            end if
         end if

          ! BEGIN MAGNETIC PHASE
          if (mphase_begin) then
            ! READ MAGNETIC LATTICE
            if (is_lattice /= 0) then
              read(unit=line(is_lattice+7:),fmt=*,iostat=ier) lattice
              if (ier /= 0) then
                ierror=.true.
              end if
            end if
            ! READ PROPAGATION VECTOR
            if (is_k /= 0) then
              num_k=num_k+1
             read(unit=line(is_k+1:),fmt=*) all_k(:,num_k)
              do  !repeat reading until continuous K lines are exhausted
                read (unit=1,fmt="(a)", iostat=ier) line
                upline=u_case(line)
                no_line=no_line+1
                is_k = index(upline,"K")
                if( is_k /= 0) then
                   num_k=num_k+1
                   read(unit=line(is_k+1:),fmt=*) all_k(:,num_k)
                else
                  no_line=no_line-1
                  backspace(unit=1)
                  k_begin=.false.
                  exit
                end if
              end do
              cycle
            end if
            ! READ XSYM OPERATORS
            if (is_xsym /= 0) then
              num_xsym=num_xsym+1
              num_msym=0
              read(unit=line(is_xsym+4:),fmt="(a)") xsym(num_xsym)
              msym_begin=.true.
            end if
            ! READ MSYM OPERATORS
            if (is_msym /= 0 .and. msym_begin) then
              num_msym=num_msym+1
              read(unit=line(is_msym:),fmt="(a)") all_msym(num_xsym,num_msym)
              do  !repeat reading until continuous MSYM lines are exhausted
                read (unit=1,fmt="(a)", iostat=ier) line
                upline=u_case(line)
                no_line=no_line+1
                is_msym = index(upline,"MSYM")
                if( is_msym /= 0) then
                   num_msym=num_msym+1
                   read(unit=line(is_msym:),fmt="(a)") all_msym(num_xsym,num_msym)
                else
                  no_line=no_line-1
                  backspace(unit=1)
                  msym_begin=.false.
                  exit
                end if
              end do
              cycle
            end if
            ! READ MATOM
            if (is_matom /=0) then
              if (mstart==0) call define_magnetic_phase(lattice,all_k,xsym,all_msym,num_k,num_xsym,num_msym)
              mstart=mstart+1
              read(unit=line(is_matom+5:),fmt=*, iostat=ier) label, symbol, pos
              if (ier /= 0) then
                ierror=.true.
              else
                matom_begin=.true.
                num_skp=0
              end if
              if (is_color  == 0) color =get_color_from_symbol(symbol)
              if (is_envelop /=0 ) then
              matom_envelop=.true.
              else
              matom_envelop=.false.
              end if
              if (is_edges /=0 ) then
              matom_edges=.true.
              else
              matom_edges=.false.
              end if
              if (is_radius /= 0)  read(unit=line(is_radius+5:),fmt=*,iostat=ier) radius
              if (is_edge_color /=0) then
                  read(unit=line(is_edge_color+7:),fmt=*,iostat=ier) edgecolor
                else
                  edgecolor=(/0.0,0.0,0.0/)
              end if
              if (is_envelop_color /=0) then
              read(unit=line(is_envelop_color+10:),fmt=*,iostat=ier) envelop_color
              matom_envcolor=.true.
              if (ier /= 0) matom_envcolor=.false.
              end if
              end if
            ! READ FOURIER COEFFICIENTS AND PHASE
            if (is_skp /= 0 .and. matom_begin) then
              num_skp=num_skp+1
              read(unit=line(is_skp+3:),fmt=*) kchoice, msymchoice, skreal_a, skreal_b,&
                                                   skreal_c, skim_a, skim_b, skim_c, phase

              my_skj(:,num_skp)=(/CMPLX(skreal_a,skim_a),CMPLX(skreal_b,skim_b),CMPLX(skreal_c,skim_c)/)
              my_phikj(num_skp)=phase
              my_color(:,num_skp)=color
              my_point2k(num_skp)=kchoice
              my_point2msym(num_skp)=msymchoice
              is_color=0
              do  !repeat reading until continuous SKP lines are exhausted
                read (unit=1,fmt="(a)", iostat=ier) line
                upline=u_case(line)
                is_color = index(upline,"COLOR")
                if (is_color/=0) then
                  read(unit=line(is_color+5:),fmt=*) color
                end if
                upline=u_case(line)
                no_line=no_line+1
                is_skp = index(upline,"SKP")
                if (is_skp /= 0) then
                   num_skp=num_skp+1
                   read(unit=line(is_skp+3:),fmt=*) kchoice, msymchoice, skreal_a, skreal_b,&
                                                        skreal_c, skim_a, skim_b, skim_c, phase

                   my_skj(:,num_skp)=(/CMPLX(skreal_a,skim_a),CMPLX(skreal_b,skim_b),CMPLX(skreal_c,skim_c)/)
                   my_phikj(num_skp)=phase
                   my_color(:,num_skp)=color
                   my_point2k(num_skp)=kchoice
                   my_point2msym(num_skp)=msymchoice
                else
                  no_line=no_line-1
                  backspace(unit=1)
                  matom_begin=.false.
                  if (matom_envcolor) then
                  call push_matom(label, symbol,pos, my_point2k, my_point2msym, my_Skj,&
                                   my_phikj,my_color,scal,.false.,group,num_skp,matom_envelop,envelop_color=envelop_color,edges=matom_edges,edges_color=(/edgecolor,1.0/),edges_radius=radius)
                  else
				  call push_matom(label, symbol,pos, my_point2k, my_point2msym, my_Skj,&
                                   my_phikj,my_color,scal,.false.,group,num_skp,matom_envelop,edges=matom_edges,edges_color=(/edgecolor,1.0/),edges_radius=radius)
                  end if
                  matom_envcolor=.false.
                  group=.false.
                  envelop=.false.
                  edges=.false.
                  scal=1.0
                  exit
                end if
              end do
              cycle
            end if
            ! END OF READ COEFFICIENTS AND PHASE
          end if
          ! END OF MAGNETIC PHASE
        end if
      end do
      close(unit=1)
      if(molecule) then
        call define_box(pos1(1),pos2(1),pos1(2),pos2(2),pos1(3),pos2(3))
      end if
    end if
    return

  end subroutine read_fst_file
!------------------------------------------------------------------------------
  subroutine write_fst_file(filename,qview)
    character(len=*), intent(in) :: filename
    type(quat),       intent(in) :: qview
    character(len=9)  :: multiple_string
    character(len=10) :: dead_string
    real              :: phi
    real, dimension(3):: axis
    integer           :: ier,i
    logical           :: op

    inquire(file=filename,opened=op)
    if(op) then
      error_io=.true.
      mess_error_io= " The file is being edited! Please close the file before saving!"
      return
    end if
    open(unit=1,file=filename,action="write",status="replace", iostat=ier)
    if(ier /= 0) then
      error_io=.true.
      mess_error_io= " Error opening the file "//trim(filename)
      return
    end if
     ! BACKGROUND COLOR
        write(unit=1,fmt='(a,4F7.3)') "BKG", background_color
     ! CELL
      if (associated(current_cell)) then

        if (current_cell%multiple) then
          multiple_string=" MULTIPLE"
        else
          multiple_string=" "
        end if

        if (current_cell%dead) then
          dead_string=" NODISPLAY"
        else
          dead_string=" "
        end if

          write(unit=1,fmt='(a,6F10.4,a,4f7.3,2a)') "CELL ",&
          current_cell%a, current_cell%angle, " COLOR", current_cell%color, &
          multiple_string, dead_string
      end if
     ! BOX SIZE
      if (associated(current_box)) &
        write(unit=1,fmt='(a,6F8.3)') "BOX", current_box

     ! CURRENT ORIENTATION
        call rot_from_q(qview,phi,axis)
        write(unit=1,fmt='(a,f9.3,3F9.5)') "ROTAX ", phi,axis

     ! SPACE GROUP
      if (associated(current_space_group)) then
      	 if(current_space_group%SPG_Symb == "unknown") then
      	  do i=1,current_space_group%multip
      	    write(unit=1,fmt='(a)')"GENR "//trim(current_space_group%symopsymb(i))
      	  end do
      	 else
           write(unit=1,fmt='(2a)') "SPACEG ", current_space_group%SPG_Symb
         end if
      end if
     ! ATOMS
      call write_atom_list(1)
     ! BONDS
      call write_bond_list(1)
     ! MAGNETIC STRUCTURE
      call write_magnetic_structure(1)
     ! POLYHEDRAS
      call write_poly_list(1)
     ! CLOSE
      close(unit=1)

    return

  end subroutine write_fst_file
!------------------------------------------------------------------------------
end module WCFGL_IO