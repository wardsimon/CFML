!     Last change:  TR   11 Oct 2007    6:53 pm

 ! lecture fichier .CIF
 !  source: EDPCR (JGP)

! lecture d'un fichier CIF

subroutine read_CIF_input_file(input_file, input_string)
 USE macros_module
 USE IO_module,                        ONLY : write_info
 USE cryscalc_module
 USE CFML_String_Utilities,            ONLY : number_lines,   reading_lines
 USE CFML_IO_Formats,                  ONLY : Read_CIF_title, read_CIF_cell, read_cif_hm, read_cif_hall,  &
                                              read_cif_symm,  Read_Cif_Cont, Read_Cif_atom
 USE CFML_Crystal_Metrics,             ONLY : Set_Crystal_Cell, convert_u_betas, U_Equiv
 USE CFML_Crystallographic_Symmetry,   ONLY : set_spacegroup,   Space_Group_Type, get_hallsymb_from_gener,  &
                                              Get_Multip_Pos, Wyckoff_type, Get_orbit, Get_stabilizer
 USE CFML_Atom_TypeDef,                ONLY : Atom_list_Type,  Deallocate_atom_list
 USE CFML_GlobalDeps,                  ONLY : sp, cp
 USE CFML_Math_General,                ONLY : set_epsg, set_epsg_default




 implicit none
  CHARACTER (LEN=*), INTENT(IN)            :: input_file
  CHARACTER (len=*), INTENT(IN)            :: input_string
  INTEGER                                  :: i, ier, long_input_string

  ! local variable for .CIF file
  integer                                     :: nb_lines, npos
  character(len=256),dimension(:),allocatable :: file_cif
  character(len=60),dimension(48)             :: car_symop

  integer                                     :: n_elem_atm    ! N. of different species
  !real            ,  dimension(15)           :: n_elem        ! Number of elementos into the same species
  character(len=2),  dimension(15)            :: elem_atm      ! Character to identify the specie
  character(len=40)                           :: symb_spgr
  logical                                     :: input_out
  logical                                     :: tmp_on_screen, input_string_cell

 ! sept 2014 -------------------------------------------------
  REAL (kind=cp), DIMENSION(3)               :: r
  !type (Wyckoff_type)                        :: wyckoff      !
  !real (kind=cp) , dimension(3,192)          :: orbit
  !integer, dimension(192)                    :: effsym  ! Pointer to effective symops
  !integer, dimension(192)                    :: ss_ptr  ! Pointer to stabilizer    >> oct. 2011
  !integer                                    :: order
  !real(kind=cp),     dimension(3,48)         :: atr
 ! -----------------------------------------------------------


  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_CIF_INPUT_FILE ("//trim(input_string)//")")

  IF (INI_create%PRF .and. .not. ON_SCREEN_PRF) then
   tmp_on_screen = on_screen
   on_screen = .false.
  else
   tmp_on_screen = .true.
  end if

  long_input_string = len_trim(input_string)
  input_out           = .true.
  input_string_cell   = .false.
  if(long_input_string == 6) then
   if(input_string(1:6) == 'NO_OUT') input_out = .false.
  elseif(long_input_string == 4) then
   if(input_string(1:4) == 'CELL')   input_string_CELL = .true.
  endif



  close(unit=CIF_read_unit)  ! indispensable pour G95

  call number_lines(trim(input_file),nb_lines)
  if (nb_lines ==0) then
   call write_info(' > No lines can be read. Program will be stopped')
   stop
  end if

  if (allocated(file_cif)) deallocate(file_cif)

  allocate(file_cif(nb_lines), stat=ier)
  if(ier/=0) call write_alloc_error("file_cif")


  file_cif=' '

 !---- Cargando Fichero en variable ----!
  call reading_lines(trim(input_file), nb_lines, file_cif)

  !---- TITL ----!
  npos=1
  call Read_Cif_Title(file_cif, npos, nb_lines, main_title)
  if(input_out .and. ON_SCREEN .and. write_details) then
   call write_info(' ' )
   call write_info('  . TITL: '//trim(main_title))
  endif

  !---- CELL ----!
  npos=1
  call Read_Cif_Cell(file_cif, npos, nb_lines, unit_cell%param, unit_cell%param_esd)
  known_cell_esd = .true.
  !call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)
  call create_CELL_object()

  if(.not. input_out) return
  IF(unit_cell%volume < 0.1) call volume_calculation('no_out')  ! << oct. 2011

  keyword_CELL  = .true.
  if(ON_SCREEN .and. write_details) then
   call write_info(' ')
   write(message_text,'(a,6F10.4)') '  . CELL: ', (unit_cell%param(i), i=1,6)
   call write_info(TRIM(message_text))
   write(message_text,'(a,6F10.4)') '          ', (unit_cell%param_esd(i), i=1,6)
   call write_info(TRIM(message_text))
  endif

  IF(input_string_CELL) return
  if(FCF_plot_stl)      return


  !---- OBTAIN SPACE GROUP ----!
  symb_spgr=' '
  npos=1
  call read_cif_hm(file_cif ,npos, nb_lines, symb_spgr)

  if(index(symb_spgr, '::R') /=0) symb_spgr = replace_car(symb_spgr, '::R', ':R')
  if(index(symb_spgr, ':H')  /=0) symb_spgr = symb_spgr(1:index(symb_spgr, ':H')-2)
  !if(index(symb_spgr, ' H')  /=0) symb_spgr = symb_spgr(1:index(symb_spgr, ' H')-1)
  if(index(symb_spgr, ':')   /=0) symb_spgr = symb_spgr(1:index(symb_spgr, ':') -1)

  if(symb_spgr(2:2) /= ' ') symb_spgr = symb_spgr(1:1) // ' ' // symb_spgr(2:)

  if (len_trim(symb_spgr) > 0) then
   call Set_SpaceGroup(symb_spgr, SPG)
   if(ON_SCREEN .and. write_details) then
    call write_info(' ')
    write(message_text, '(a)') '     >> SPACE GROUP extracted from "_symmetry_space_group_name_H-M" '
    call write_info(trim(message_text))
    write(message_text, '(a)') '                               or  "_space_group_name_H-M_alt" CIF field'
    call write_info(trim(message_text))
   endif

  else
   if(ON_SCREEN .and. write_details .and. npos /=1) then
    call write_info(' ')
    call write_info('     >> Space group can not be extracted from "_symmetry_space_group_name_H-M CIF" field.')
    call write_info('     >> Probable causes: . "_symmetry_space_group_name_H-M CIF" field is missing.')
    call write_info('                         . quotes are missing for the space group symbol.')
    call write_info(' ')
   end if
   npos=1
   call read_cif_hall(file_cif ,npos, nb_lines, symb_spgr)

   if (len_trim(symb_spgr) > 0) then
    call Set_SpaceGroup(symb_spgr, SPG)
    if(ON_SCREEN .and. write_details) then
     call write_info(' ')
     call write_info('     >> SPACE GROUP extracted from "_symmetry_space_group_name_Hall" CIF field')
    endif

   else
    npos=1
    call read_cif_symm(file_cif, npos ,nb_lines, nb_symm_op, car_symop)
    if(ON_SCREEN .and. write_details .and. nb_symm_op == 0) then
     call write_info('')
     call write_info('     >> No symmetry operator: Space group can not be deduced !!')
     call write_info('')
     return
    end if
    symb_spgr=' '
    call set_spacegroup("  ", SPG, car_symop, nb_symm_op,'GEN')

    !call get_hallsymb_from_gener(SPG)
    if(on_screen .and. write_details) then
     if(SPG%NumSpg == 0) then
      call write_info('')
      call write_info('     >> Space group can not be deduced from symmetry operators list !!')
     else
      call write_info('')
      call write_info('     >> Space group deduced from symmetry operators list !!')
     end if
    end if
   end if
  end if

  space_group_symbol = SPG%Spg_Symb
  !space_group_multip = SPG%multip
  !IF(LEN_TRIM(space_group_symbol) /=0 .and. SPG%NumSPG /=0) then
  IF(SPG%NumSPG /=0) then
   keyword_SPGR   = .true.
   WRITE_SPG_info = .true.
   if(ON_SCREEN .and. write_details) then
    call write_info(' ')
    call write_info('  . SPACE GROUP: '//TRIM(space_group_symbol))
   endif
  else
   keyword_SPGR   = .false.
   WRITE_SPG_info = .false.
   if(ON_SCREEN .and. write_details) then
    call write_info(' ')
    call write_info('     >> Unknown space group !!')
   endif
  endif


  !---- ATOMS ----!
  npos=1
  call Read_Cif_Cont(file_cif, npos, nb_lines, n_elem_atm, elem_atm)

  npos=1
  call Read_Cif_Atom(file_cif, npos, nb_lines, nb_atom, Atm_list)

  IF(nb_atom /=0) then
   if(ON_SCREEN .and. write_details) then
    call write_info(' ')
    call write_info('  . ATOM LIST: ')
    call write_info(' ')
   endif


   do i=1,nb_atom
    atom_label(i)     = Atm_list%atom(i)%lab
    atom_typ(i)       = Atm_list%atom(i)%ChemSymb
    atom_coord(1:3,i) = Atm_list%atom(i)%x
    !atom_mult(i)      = Atm_list%atom(i)%mult

    ! sept. 2014 -------------------
    r(1:3)= atom_coord(1:3,i)  ! coord. atomiques
    call get_site_multiplicity(r, atom_mult(i))
    !if(atom_mult(i) == 0) then
    !r(1:3)= atom_coord(1:3,i)  ! coord. atomiques
    !call Get_Orbit(r, SPG, wyckoff%num_orbit, orbit, effsym)
    !call Get_stabilizer(r, SPG, order, ss_ptr, atr)
    !atom_mult(i)      = Get_Multip_Pos(r, SPG)
    !end if
    !-------------------------------
    atom_occ(i)       = Atm_list%atom(i)%occ * real (Get_Multip_Pos(atom_coord(1:3,i) ,SPG)) / real(SPG%multip)
    atom_occ_perc(i)  = Atm_list%atom(i)%occ

    if (Atm_list%atom(i)%thtype=='aniso') then
 !    n_t(i)=2
     atom_adp_aniso(1:6,i) = Atm_list%atom(i)%u(1:6)
     !atom_adp_aniso(1:6,i) = convert_u_betas(atom_adp_aniso(1:6,i), crystal_cell)
     atom_adp_equiv(i)     = 8.0*pi*pi*U_Equiv(crystal_cell, Atm_list%atom(i)%u(1:6))
    else
 !    n_t(i)=0
      !atom_adp_equiv(i)    = Atm_list%atom(i)%u(7)*8.0*pi*pi
      atom_adp_equiv(i)    = Atm_list%atom(i)%Ueq
      atom_adp_equiv(i)    = 8.0*pi*pi*Atm_list%atom(i)%Ueq
    end if

    atom_Biso(i)      = atom_adp_equiv(i)

     if(ON_SCREEN .and. write_details) then
      !write(message_text,'(2x,2(1x,a),4F10.5,I4)')  atom_label(i)(1:6), atom_typ (i)(1:6) , &
  !                                              atom_coord(1:3,i), atom_occ_perc(i), atom_mult(i)
      !call write_info(TRIM(message_text))
      if(Atm_list%atom(i)%thtype=='aniso') then
       write(message_text,'(2x,2(1x,a6),2x,4F10.5,I4,6F10.5)')  atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
                                                                atom_occ_perc(i), atom_mult(i), Atm_list%atom(i)%u(1:6)
      else
       write(message_text,'(2x,2(1x,a6),2x,4F10.5,I4, F10.5)')  atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
                                                                atom_occ_perc(i), atom_mult(i), Atm_list%atom(i)%Ueq
      endif
      call write_info(TRIM(message_text))
     endif

   end do
  else
   if(on_screen .and. write_details) then
    call write_info('')
    call write_info('  > No atoms founded in the CIF file !')
    call write_info('')
   end if
  endif

  call Deallocate_atom_list(Atm_list)
  if (allocated(file_cif)) deallocate(file_cif)

  ! pour compatibilite avec CRYSCALC
   nb_atoms_type                = n_elem_atm
   SFAC%type(1:nb_atoms_type)   = elem_atm(1:n_elem_atm)
   !SFAC%number(1:nb_atoms_type) = n_elem(1:nb_atoms_type)

   !call set_epsg_default()

   IF (INI_create%PRF .and. .not. ON_SCREEN_PRF) on_screen = tmp_on_screen

return
end subroutine read_CIF_input_file
!------------------------------------------------------------------------------

subroutine read_CIF_input_file_TR(input_unit, input_CIF_file)
 USE cryscalc_module, ONLY : keyword_WAVE, wavelength, keyword_SIZE, Z_unit, keyword_ZUNIT,    &
                             crystal, unit_cell, nb_atoms_type, molecule,                      &
                             SADABS, SFAC, sto, Z_unit,                                        &
                             keyword_create_ARCHIVE, on_screen, on_screen_CIF_item, debug_proc
 USE CIF_module
 USE macros_module,   ONLY : nombre_de_colonnes, remove_car
 USE IO_module,       ONLY : write_info
 implicit none
  INTEGER, INTENT(IN)            :: input_unit
  CHARACTER (LEN=*), INTENT(IN)  :: input_CIF_file
  CHARACTER (LEN=256)            :: string_value, read_line
  CHARACTER (LEN=256)            :: CIF_string
  INTEGER                        :: nb_arg, i1, i2, long
  LOGICAL                        :: ok

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_CIF_INPUT_FILE_TR")

  open (unit= input_unit, file=trim(input_CIF_file)) ! le fichier etait ferme ?
  ! sample ID
   call get_sample_ID(input_unit)


  ! wave:
  CIF_string = '_diffrn_radiation_wavelength'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok) then
   keyword_WAVE = .true.
   READ(string_value, *) wavelength
   read(string_value, *) CIF_cell_measurement%wavelength
   call incident_beam
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  endif


  ! crystal description
  rewind(unit=input_unit)
  CIF_string = '_exptl_crystal_description'
  call get_champ_value(input_unit,  trim(CIF_string), string_value, ok)

  if(ok .and. crystal%morph(1:1) == '?') then
   crystal%morph = trim(string_value)
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  ! crystal colour
  CIF_string = '_exptl_crystal_colour'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. crystal%color(1:1) == '?') then
   crystal%color = trim(string_value)
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  ! SIZE
  CIF_string = '_exptl_crystal_size_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok) then
   IF(string_value(1:1) /= '?')  then
    READ(string_value, *)  crystal%size(1)
    READ(string_value, *)  crystal%size_max
    if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
   endif
   CIF_string = '_exptl_crystal_size_mid'
   call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
   IF(ok) then
    IF(string_value(1:1) /= '?') then
     READ(string_value, *) crystal%size(2)
     READ(string_value, *) crystal%size_mid
     if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
    endif
    CIF_string = '_exptl_crystal_size_min'
    call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
    IF(ok) then
     IF(string_value(1:1) /= '?') then
      READ(string_value, *) crystal%size(3)
      READ(string_value, *) crystal%size_min
      if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
     endif
     keyword_SIZE = .true.
    endif
   endif
   if(keyword_size) call crystal_volume_calculation('no_out')
  endif

  IF(crystal%size(1) < 0.000001 .AND.  crystal%size(2) < 0.000001  .AND. crystal%size(3) < 0.000001) keyword_size = .false.


 ! ZUNIT
  CIF_string = '_cell_formula_units_Z'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok .and. string_value(1:1) /= '1' ) then
   keyword_ZUNIT = .true.
   READ(string_value, *) Z_unit
   READ(string_value, '(a)') CIF_parameter%cell_formula_units_Z
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  endif
  molecule%Z_unit = int(Z_unit)

 ! chemical formula
  CIF_string = '_chemical_formula_sum'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  string_value = remove_car(string_value, "'")
  string_value = adjustl(string_value)
  long = len_trim(string_value)
  IF(ok .and. long/=0 .and. string_value(1:1) /= '?') then
   !call remove_car(string_value, "'")
   string_value = remove_car(string_value, "'")
   call nombre_de_colonnes(string_value, nb_arg)
   call decode_CHEM_string(string_value, nb_arg)
   molecule%formula = string_value
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  endif

  !call Get_CIF_value(input_unit, string_value, molecule%formula)

   ! chemical formula weight
  CIF_string = '_chemical_formula_weight'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok .and. CIF_parameter%formula_weight(1:1) == '?') then
   read(string_value, '(a)') CIF_parameter%formula_weight
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  endif


  !! july 2015
  ! ! space group
  ! CIF_string = '_symmetry_space_group_name_H-M'
  ! call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  ! IF(ok) then
  !  READ(string_value, *) CIF_parameter%symmetry_space_group
  !  if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  ! end if
  !! -------------------------------------------------------------------

    ! wave:
  CIF_string = '_cell_measurement_wavelength'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok) then
   READ(string_value, *) CIF_cell_measurement%wavelength
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if



! cell parameters and volume
  CIF_string =  '_cell_length_a'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok .and. CIF_parameter%CELL_length_a(1:1) == '?') then
   READ(string_value, '(a)') CIF_parameter%CELL_length_a
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string =  '_cell_length_b'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok .and. CIF_parameter%CELL_length_b(1:1) == '?') then
   READ(string_value, '(a)') CIF_parameter%CELL_length_b
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_cell_length_c'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok .and. CIF_parameter%CELL_length_c(1:1) == '?') then
   READ(string_value, '(a)') CIF_parameter%CELL_length_c
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  endif

  CIF_string = '_cell_angle_alpha'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok .and. CIF_parameter%CELL_angle_alpha(1:1) == '?') then
   READ(string_value, '(a)') CIF_parameter%CELL_angle_alpha
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_cell_angle_beta'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok .and. CIF_parameter%CELL_angle_beta(1:1) == '?') then
   READ(string_value, '(a)') CIF_parameter%CELL_angle_beta
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_cell_angle_gamma'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  IF(ok .and. CIF_parameter%CELL_angle_gamma(1:1) == '?') then
   READ(string_value, '(a)') CIF_parameter%CELL_angle_gamma
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_cell_volume'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%CELL_volume(1:1) == '?') then
   read(string_value, '(a)') CIF_parameter%Cell_volume
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  if(CIF_parameter%Cell_volume(1:1) /='?') then
   i1 = INDEX(string_value, '(')
   i2 = INDEX(string_value, ')')
   IF(i2 > i1) then
    READ(string_value(1:i1-1),   *) unit_cell%volume
    READ(string_value(i1+1:i2-1),*) unit_cell%volume_ESD
   endif
  endif


  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%CELL_volume(1:1) == '?') then
   read(string_value, '(a)') CIF_parameter%Cell_volume
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  if(.not. keyword_create_archive) return

 ! density
  CIF_string = '_exptl_crystal_density_meas'
  call get_champ_value(input_unit, trim(CIF_string),   string_value, ok)
  if(ok .and. crystal%density_meas(1:1) =='?') then
   crystal%density_meas = trim(string_value)
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_exptl_crystal_density_diffrn'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. crystal%density_diffrn(1:1) =='?') then
   crystal%density_diffrn = trim(string_value)
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_exptl_crystal_density_method'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. crystal%density_method(1:1) =='?') then
   crystal%density_method = trim(string_value)
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  ! F000
  CIF_string = '_exptl_crystal_F_000'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. crystal%F000(1:1) == '?' ) then
   crystal%F000 = trim(string_value)
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  ! THmin, THmax
  CIF_string = '_cell_measurement_reflns_used'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if (ok  .and. CIF_cell_measurement%reflns_used < 0.) then
   READ(string_value, *) CIF_cell_measurement%reflns_used
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_cell_measurement_theta_min'
  call get_champ_value(input_unit, trim(CIF_string),    string_value, ok)
  if (ok  .and. CIF_cell_measurement%theta_min < 0.) then
   READ(string_value, *) CIF_cell_measurement%theta_min
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_cell_measurement_theta_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if (ok .and. CIF_cell_measurement%theta_max < 0.) then
   READ(string_value, *) CIF_cell_measurement%theta_max
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_theta_min'
  call get_champ_value(input_unit, trim(CIF_string),    string_value, ok)
  if (ok .and. CIF_diffrn_reflns%theta_min < 0.) then
   READ(string_value, *) CIF_diffrn_reflns%theta_min
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_diffrn_reflns_theta_max'
  call get_champ_value(input_unit, trim(CIF_string),    string_value, ok)
  if (ok .and. CIF_diffrn_reflns%theta_max < 0.) then
   READ(string_value, *) CIF_diffrn_reflns%theta_max
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if


  ! temperature
  CIF_string = '_cell_measurement_temperature'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  !if (ok) READ(string_value, *) CIF_cell_measurement%temperature
  if (ok .and. CIF_cell_measurement%temperature(1:1) == '?') then
   READ(string_value, '(a)') CIF_cell_measurement%temperature
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  ! absorption
  CIF_string = '_exptl_absorpt_coefficient_mu'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)

  if(ok .and. SADABS%absorption_coeff(1:1) == '?') then
   read(string_value,'(a)') SADABS%absorption_coeff
   read(string_value,'(a)') CIF_parameter%absorption_coefficient_mu
   if(ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  !
  CIF_string = '_diffrn_reflns_av_R_equivalents'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%diffrn_reflns_av_R_equivalents == '?') then
   read(string_value,'(a)') CIF_parameter%diffrn_reflns_av_R_equivalents
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string =  '_diffrn_reflns_av_sigmai/neti'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%diffrn_reflns_av_R_sigma == '?')  then
   read(string_value,'(a)') CIF_parameter%diffrn_reflns_av_R_sigma
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_av_uneti/neti'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%diffrn_reflns_av_R_sigma == '?') then
   read(string_value,'(a)') CIF_parameter%diffrn_reflns_av_R_sigma
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_number'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%diffrn_reflns_number == '?') then
   read(string_value, '(a)') CIF_parameter%diffrn_reflns_number
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_limit_h_min'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%h_min == '?') then
   read(string_value, *) CIF_parameter%h_min
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string =  '_diffrn_reflns_limit_h_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%h_max == '?') then
   read(string_value, *) CIF_parameter%h_max
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_limit_k_min'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%k_min == '?') then
   read(string_value, *) CIF_parameter%k_min
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_limit_k_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%k_max == '?') then
   read(string_value, *) CIF_parameter%k_max
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_limit_l_min'
  call get_champ_value(input_unit, trim(CIF_string) , string_value, ok)
  if(ok .and. CIF_parameter%l_min == '?') then
   read(string_value, *) CIF_parameter%l_min
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_limit_l_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%l_max == '?') then
   read(string_value, *) CIF_parameter%l_max
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_theta_min'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%theta_min == '?') then
   read(string_value, *) CIF_parameter%theta_min
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_theta_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%theta_max == '?') then
   read(string_value, *) CIF_parameter%theta_max
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_reflns_theta_full'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%theta_full == '?') then
   read(string_value, *) CIF_parameter%theta_full
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_measured_fraction_theta_full'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  !if(ok .and. CIF_parameter_DEVICE%diffrn_theta_full == '?') then
  if(ok) then
   read(string_value, '(a)') CIF_parameter_DEVICE%diffrn_theta_full
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_diffrn_measured_fraction_theta_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  !if(ok .and. CIF_parameter_DEVICE%diffrn_theta_max == '?') then
  if(ok)  then
   read(string_value, '(a)') CIF_parameter_DEVICE%diffrn_theta_max
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_reflns_number_total'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%reflns_number_total == '?') then
   read(string_value, '(a)') CIF_parameter%reflns_number_total
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_reflns_number_gt'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%reflns_number_gt == '?') then
   read(string_value, '(a)') CIF_parameter%reflns_number_gt
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_refine_ls_weighting_details'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_weighting_details == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_weighting_details
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_reflns_Friedel_coverage'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%Friedel_coverage == '?') then
   read(string_value, '(a)') CIF_parameter%Friedel_coverage
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_reflns_Friedel_fraction_full'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%Friedel_fraction_full == '?') then
   read(string_value, '(a)') CIF_parameter%Friedel_fraction_full
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_reflns_Friedel_fraction_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%Friedel_fraction_max == '?') then
   read(string_value, '(a)') CIF_parameter%Friedel_fraction_max
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if


  CIF_string = '_computing_structure_refinement'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. string_value(1:1) /= '?') then
   read(string_value, '(a)') CIF_parameter_DEVICE%computing_structure_refinement
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if


  CIF_string = '_atom_sites_solution_primary'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%atom_sites_solution_1 == '?') then
   read(string_value, '(a)') CIF_parameter%atom_sites_solution_1
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_atom_sites_solution_secondary'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%atom_sites_solution_2 == '?') then
   read(string_value, '(a)') CIF_parameter%atom_sites_solution_2
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_atom_sites_solution_hydrogens'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%atom_sites_solution_H == '?') then
   read(string_value, '(a)') CIF_parameter%atom_sites_solution_H
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  ! H treatment
  call determine_H_treatment(input_unit, CIF_parameter%refine_ls_H_treatment)


  CIF_string = '_refine_ls_extinction_method'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_extinction_method == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_extinction_method
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_extinction_coef'
  call get_champ_value(input_unit, trim(CIF_string) , string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_extinction_coef == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_extinction_coef
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_refine_ls_number_reflns'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_number_reflns == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_number_reflns
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_number_parameters'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_number_parameters == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_number_parameters
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_number_restraints'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_number_restraints == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_number_restraints
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_R_factor_all'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_R_factor_all == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_R_factor_all
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_R_factor_gt'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_R_factor_gt == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_R_factor_gt
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_wR_factor_ref'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_wR_factor_ref == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_wR_factor_ref
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_wR_factor_gt'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_wR_factor_gt == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_wR_factor_gt
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_refine_ls_goodness_of_fit_ref'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_goodness_of_fit_ref == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_goodness_of_fit_ref
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_restrained_S_all'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_restrained_S_all == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_restrained_S_all
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_shift/su_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_shift_su_max == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_shift_su_max
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_shift/su_mean'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_ls_shift_su_mean == '?') then
   read(string_value, '(a)') CIF_parameter%refine_ls_shift_su_mean
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_refine_ls_abs_structure_flack'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%abs_structure_Flack == '?') then
   read(string_value, '(a)') CIF_parameter%abs_structure_Flack
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_ls_abs_structure_details'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%abs_structure_details(1) == '?') then
   n_details = 1
   if(string_value(1:1) == ';') then
    CIF_parameter%abs_structure_details(n_details) = ";"
    do
    read(unit = input_unit , fmt='(a)')  read_line
    read_line = adjustl(read_line)
    n_details = n_details + 1
    read(read_line, '(a)') CIF_parameter%abs_structure_details(n_details)
    if(read_line(1:1) == ';') exit
   end do
   else
    read(string_value, '(a)') CIF_parameter%abs_structure_details(n_details)
   endif
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if


  CIF_string = '_refine_diff_density_max'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_diff_density_max == '?') then
   read(string_value, '(a)') CIF_parameter%refine_diff_density_max
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_diff_density_min'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_diff_density_min == '?') then
   read(string_value, '(a)') CIF_parameter%refine_diff_density_min
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if
  CIF_string = '_refine_diff_density_rms'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%refine_diff_density_rms == '?') then
   read(string_value, '(a)') CIF_parameter%refine_diff_density_rms
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if


  CIF_string = '_shelx_res_checksum'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%shelx_res_checksum == '?') then
   read(string_value, '(a)') CIF_parameter%shelx_res_checksum
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_shelx_hkl_checksum'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%shelx_hkl_checksum == '?') then
   read(string_value, '(a)') CIF_parameter%shelx_hkl_checksum
   !if(ON_screen_CIF_item) call write_info('  > '//CIF_string(1:48)//'extracted from '//trim(input_CIF_file))
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  CIF_string = '_shelx_fab_checksum'
  call get_champ_value(input_unit, trim(CIF_string), string_value, ok)
  if(ok .and. CIF_parameter%shelx_fab_checksum == '?') then
   read(string_value, '(a)') CIF_parameter%shelx_fab_checksum
   !if(ON_screen_CIF_item) call write_info('  > '//CIF_string(1:48)//'extracted from '//trim(input_CIF_file))
   if (ON_screen_CIF_item) call Write_extracted_cif_parameter(CIF_string(1:48), trim(input_CIF_file), trim(string_value))
  end if

  return
end subroutine read_CIF_input_file_TR

!-------------------------------------------------------------------------------
 subroutine Write_extracted_cif_parameter(input_string_1, input_CIF_file, input_string_2)
  USE IO_module
  implicit none
  character (len=*), intent(in) :: input_string_1, input_CIF_file, input_string_2

  call write_info('  > '//input_string_1//'extracted from '//trim(input_CIF_file)//' ['//trim(input_string_2)//']')

  return
 end subroutine Write_extracted_cif_parameter

!-------------------------------------------------------------------------------
subroutine get_cell_parameters_from_cif_file(file_unit, cell_param, cell_param_esd)
 USE IO_module
 USE cryscalc_module, only : debug_proc
 USE CIF_module,      only : CIF_parameter
 implicit none
  INTEGER, INTENT(IN)                  :: file_unit
  REAL,    INTENT(OUT),   dimension(*) :: cell_param, cell_param_esd
  character (len=256)                  :: line
  integer                              :: ier
  character (len=64)                   :: cell_a_string, cell_b_string, cell_c_string
  character (len=64)                   :: cell_alfa_string, cell_beta_string, cell_gamma_string
  LOGICAL , DIMENSION(6)               :: logical_cell

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_CELL_PARAMETERS_FROM_CIF_FILE")

  logical_cell(1:6) = .false.
  REWIND(UNIT=file_unit)
  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info(' ... No cell parameters in the import.cif file !!')
    call write_info('')
    stop
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading import.cif file !!')
    call write_info('')
    stop
   endif

   if(trim(line(1:14))=='_cell_length_a') then
    read(line(16:),'(a)') cell_a_string
    call get_cell_param(cell_a_string, cell_param(1), cell_param_esd(1), 'a')
    logical_cell(1) = .true.
    CIF_parameter%cell_length_a = adjustl(cell_a_string)


   elseif(trim(line(1:14))=='_cell_length_b') then
    read(line(16:),'(a)') cell_b_string
    call get_cell_param(cell_b_string, cell_param(2), cell_param_esd(2), 'b')
    logical_cell(2) = .true.
    CIF_parameter%cell_length_b = adjustl(cell_b_string)


   elseif(trim(line(1:14))=='_cell_length_c') then
    read(line(16:),'(a)') cell_c_string
    call get_cell_param(cell_c_string, cell_param(3), cell_param_esd(3), 'c')
    logical_cell(3) = .true.
    CIF_parameter%cell_length_c = adjustl(cell_c_string)


   elseif(trim(line(1:17))=='_cell_angle_alpha') then
    read(line(19:),'(a)') cell_alfa_string
    call get_cell_param(cell_alfa_string, cell_param(4), cell_param_esd(4), 'alpha')
    logical_cell(4) = .true.
    CIF_parameter%cell_angle_alpha = adjustl(cell_alfa_string)


   elseif(trim(line(1:16))=='_cell_angle_beta') then
    read(line(18:),'(a)') cell_beta_string
    call get_cell_param(cell_beta_string, cell_param(5), cell_param_esd(5), 'beta')
    logical_cell(5) = .true.
    CIF_parameter%cell_angle_beta = adjustl(cell_beta_string)

   elseif(trim(line(1:17))=='_cell_angle_gamma') then
    read(line(19:),'(a)') cell_gamma_string
    call get_cell_param(cell_gamma_string, cell_param(6), cell_param_esd(6), 'gamma')
    logical_cell(6) = .true.
    CIF_parameter%cell_angle_gamma = adjustl(cell_gamma_string)
   endif

   IF(logical_cell(1) .AND. logical_cell(2) .AND. LOGICAL_cell(3)) then
    if(logical_cell(4) .AND. logical_cell(5) .AND. LOGICAL_cell(6)) then
     !call write_cell_parameters(cell_param)
     exit
    endif
   endif
  end do



 RETURN
end subroutine get_cell_parameters_from_cif_file

!----------------------------------------------------------------------
subroutine get_CIF_parameters_from_import_CIF(file_unit)
 use CRYSCALC_module, only : debug_proc, crystal, molecule, keyword_SIZE, keyword_CHEM, &
                             SFAC, nb_atoms_type, sto,Z_unit, ON_screen, WRITE_details
 USE CIF_module,      only : CIF_parameter, CIF_parameter_DEVICE, CIF_CELL_measurement
 USE IO_module,       only : write_info
 USE macros_module,   only : nombre_de_colonnes

 implicit none
  integer, intent (in) :: file_unit
  character (len=256)                  :: line
  integer                              :: ier, nb_col, long
  character (len=256)                  :: cif_string

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_CIF_PARAMETERS_FROM_CIF_FILE")


  ON_screen = WRITE_details

  REWIND(UNIT=file_unit)
  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    exit
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading import.cif file !!')
    call write_info('')
    stop
   endif


   if(trim(line(1:21))=='_chemical_formula_sum') then
    read(line(22:),'(a)') molecule%formula
    molecule%formula = adjustl(molecule%formula)
    long = len_trim(molecule%formula)
    !if(molecule%formula(1:1) == "'" .and. molecule%formula(long:long) == "'") molecule%formula = molecule%formula(2:long-1)
    if(molecule%formula(1:1) == "'") molecule%formula = molecule%formula(2:)
    long = len_trim(molecule%formula)
    if(molecule%formula(long:long) == "'") molecule%formula = molecule%formula(:long-1)
    if(molecule%formula(1:1) /= '?' .and. len_trim(molecule%formula)/=0) then
     call nombre_de_colonnes(molecule%formula, nb_col)
     call decode_CHEM_string(molecule%formula, nb_col)
     call atomic_identification()
     call molecular_weight
     SFAC%number(1:nb_atoms_type) =  sto(1:nb_atoms_type) * Z_unit
    end if


   elseif(trim(line(1:29))=='_cell_measurement_temperature') then
    read(line(30:),'(a)') CIF_cell_measurement%temperature
    CIF_cell_measurement%temperature = adjustl(CIF_cell_measurement%temperature)

   elseif(trim(line(1:29))=='_cell_measurement_reflns_used') then
    !read(line(30:), *) CIF_cell_measurement%reflns_used
    read(line(30:), '(a)') CIF_string
    CIF_string = adjustl(CIF_string)
    if(CIF_string(1:1) /= '?') then
     read(CIF_string, *,iostat=ier) CIF_cell_measurement%reflns_used
      if(ier /=0) then
      call write_info(" Error reading '_cell_measurement_reflns_used' string value !")
      stop
     end if
    end if


   elseif(trim(line(1:27))=='_cell_measurement_theta_min') then
    read(line(28:), '(a)') CIF_string
    CIF_string = adjustl(CIF_string)
    if(CIF_string(1:1) /= "?") then
     read(CIF_string, *, iostat=ier) CIF_cell_measurement%theta_min
     if(ier /=0) then
      call write_info(" Error reading '_cell_measurement_theta_min' string value !")
      stop
     end if
    end if

   elseif(trim(line(1:27))=='_cell_measurement_theta_max') then
    read(line(28:), '(a)') CIF_string
    CIF_string = adjustl(CIF_string)
    if(CIF_string(1:1) /= "?") then
     read(CIF_string, *, iostat=ier) CIF_cell_measurement%theta_max
     if(ier /=0) then
      call write_info(" Error reading '_cell_measurement_theta_max' string value !")
      stop
      end if
     end if

   elseif(trim(line(1:23))=='_diffrn_radiation_probe') then
    read(line(24:), '(a)') CIF_parameter_DEVICE%diffrn_radiation_probe
    CIF_parameter_DEVICE%diffrn_radiation_probe = adjustl(CIF_parameter_DEVICE%diffrn_radiation_probe)
   elseif(trim(line(1:31))=='_diffrn_radiation_monochromator') then
    read(line(32:), '(a)') CIF_parameter_DEVICE%diffrn_radiation_monochromator
    CIF_parameter_DEVICE%diffrn_radiation_monochromator = adjustl(CIF_parameter_DEVICE%diffrn_radiation_monochromator)
   elseif(trim(line(1:24))=='_diffrn_radiation_source') then
    read(line(25:), '(a)') CIF_parameter_DEVICE%diffrn_radiation_source
    CIF_parameter_DEVICE%diffrn_radiation_source = adjustl(CIF_parameter_DEVICE%diffrn_radiation_source)

   elseif(trim(line(1:26))=='_exptl_crystal_description') then
    read(line(27:), '(a)') crystal%morph
    crystal%morph = adjustl(crystal%morph)

   elseif(trim(line(1:21))=='_exptl_crystal_colour') then
    read(line(22:), '(a)') crystal%color
    crystal%color = adjustl(crystal%color)

   elseif(trim(line(1:23))=='_exptl_crystal_size_max') then
    read(line(24:), '(a)') CIF_string
    CIF_string = adjustl(CIF_string)
    if(CIF_string(1:1) /= "?") then
     read(CIF_string, *, iostat=ier) crystal%size(3)
     if(ier /=0) then
      call write_info(" Error reading '_exptl_crystal_size_max' string value !")
      stop
     end if
     crystal%size_max = crystal%size(3)
    end if


   elseif(trim(line(1:23))=='_exptl_crystal_size_mid') then
    read(line(24:), '(a)') CIF_string
    CIF_string = adjustl(CIF_string)
    if(CIF_string(1:1) /= "?") then
     read(CIF_string, *, iostat=ier) crystal%size(2)
     if(ier /=0) then
      call write_info(" Error reading '_exptl_crystal_size_mid' string value !")
      stop
     end if
     crystal%size_mid = crystal%size(2)
    end if

   elseif(trim(line(1:23))=='_exptl_crystal_size_min') then
    read(line(24:), '(a)') CIF_string
    CIF_string = adjustl(CIF_string)
    if(CIF_string(1:1) /= "?") then
     read(CIF_string, *, iostat=ier) crystal%size(1)
     if(ier /=0) then
      call write_info(" Error reading '_exptl_crystal_size_min' string value !")
      stop
     end if
    end if
    crystal%size_min = crystal%size(1)
    keyword_SIZE = .true.

   endif
  end do

 return
end subroutine get_CIF_parameters_from_import_CIF

!----------------------------------------------------------------------

subroutine Get_UB_matrix_from_CIF_file(file_unit, UB_matrix)
 use cryscalc_module, only : UB_mat_log, debug_proc
 implicit none
  INTEGER, INTENT(IN)                    :: file_unit
  REAL,    INTENT(OUT),   dimension(3,3) :: UB_matrix
  character (len=256)                    :: line, UB_string
  integer                                :: ier

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_UB_MATRIX_FROM_CIF_FILE")

  UB_mat_log = .false.
  REWIND(UNIT=file_unit)
  do
   read(unit=file_unit, fmt= '(a)', iostat=ier) line
   if(ier /=0) return

   if(trim(line(1:27))=='_diffrn_orient_matrix_UB_11') then
    read(line(28:),fmt='(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string, *) UB_matrix(1,1)
   elseif(trim(line(1:27))=='_diffrn_orient_matrix_UB_21') then
    read(line(28:),'(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string,*) UB_matrix(2,1)
   elseif(trim(line(1:27))=='_diffrn_orient_matrix_UB_31') then
    read(line(28:),'(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string, *) UB_matrix(3,1)
   elseif(trim(line(1:27))=='_diffrn_orient_matrix_UB_12') then
    read(line(28:),'(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string, *) UB_matrix(1,2)
   elseif(trim(line(1:27))=='_diffrn_orient_matrix_UB_22') then
    read(line(28:),'(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string, *) UB_matrix(2,2)
   elseif(trim(line(1:27))=='_diffrn_orient_matrix_UB_32') then
    read(line(28:),'(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string, *) UB_matrix(3,2)
   elseif(trim(line(1:27))=='_diffrn_orient_matrix_UB_13') then
    read(line(28:),'(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string, *) UB_matrix(1,3)
   elseif(trim(line(1:27))=='_diffrn_orient_matrix_UB_23') then
    read(line(28:),'(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string, *) UB_matrix(2,3)
   elseif(trim(line(1:27))=='_diffrn_orient_matrix_UB_33') then
    read(line(28:),'(a)', iostat=ier) UB_string
    if(ier/=0) return
    read(UB_string, *) UB_matrix(3,3)
    exit
   endif
  end do

   UB_mat_log = .true.

 return
end subroutine Get_UB_matrix_from_CIF_file
!----------------------------------------------------------------------

subroutine get_cell_param(cell_string, cell_value, cell_value_esd, input_string)
 use cryscalc_module, only : debug_proc
 implicit none
  CHARACTER (LEN=*) , INTENT(IN)        :: cell_string
  REAL              , INTENT(OUT)       :: cell_value
  REAL              , INTENT(OUT)       :: cell_value_esd
  character (len=*)                     :: input_string
  INTEGER                               :: i1, i2, i3

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_CELL_PARAM ("//trim(input_string)//")")

  ! lecture de la valeur avant la parenthese
  i1 = INDEX (cell_string, '.')
  i2 = INDEX (cell_string, '(')
  if (i2 ==0) then
   READ(cell_string, *) cell_value
  else
   READ(cell_string(1:i2-1), *) cell_value
   i3 = INDEX(cell_string, ')')
   READ(cell_string(i2+1:i3-1),*) cell_value_esd
   if(i1 /=0) then
    cell_value_esd = cell_value_esd * 10.**(i1-i2+1)
   endif
  end if



 RETURN
end subroutine get_cell_param

!---------------------------------------------------------------------------

subroutine get_wavelength_from_cif_file(file_unit, wave)
 use cryscalc_module, only : debug_proc
 use IO_module
 implicit none
  INTEGER, INTENT(IN)                  :: file_unit
  REAL,    INTENT(OUT)                 :: wave
  character (len=256)                  :: line
  CHARACTER (LEN=256)                  :: champ_text
  REAL                                 :: champ_value
  integer                              :: ier
  LOGICAL                              :: logical_wave

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_WAVELENGTH_FROM_CIF_FILE")

  logical_wave = .false.

  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info('  ... No wavelength in the import.cif file !!')
    call write_info('')
    return
   ELSEIF(ier >0) then
    call write_info(' ')
    call write_info('  ... Problem reading wavelength value in import.cif file !!')
    call write_info('')
    return
   endif

   line = ADJUSTL(line)

   ! cas d'un fichier import.CIF cree par Nonius
   !if(trim(line(1:31))=='_diffrn_radiation_wavelength_id') then
   ! READ(UNIT=file_unit, '(a)', IOSTAT=ier)line
   ! if(ier/=0) then
   !  call write_info('')
   !  call write_info(' ... Problem reading import.cif file !!')
   !  call write_info('')
   !  stop
   ! endif
   ! read(line,*) wave
   ! exit
   !endif

   if(trim(line(1:28))=='_diffrn_radiation_wavelength' .or. trim(line(1:28)) == '_cell_measurement_wavelength') then
    READ(line, *, IOSTAT=ier ) champ_text, champ_value
    !if(ier/=0) then
    ! call write_info('')
    ! call write_info(' ... Problem reading wavelength in import.cif file !!')
    ! call write_info('')
    ! stop
    !endif
    IF(ier/=0) then
     ! cas d'un fichier import.cif cree par Nonius
     READ(file_unit, '(a)', IOSTAT=ier) line
     IF(ier/=0) then
      call write_info('')
      call write_info(' ... Problem reading wavelength !!')
      call write_info('')
      stop
     endif
     if(trim(line(1:31))=='_diffrn_radiation_wavelength_id') then
      READ(file_unit, '(a)', IOSTAT=ier)line
      if(ier/=0) then
       call write_info('')
       call write_info(' ... Problem reading import.cif file !!')
       call write_info('')
       stop
      END if
     endif
     read(line,*) wave
     exit
    endif

    IF(champ_value > 0.01) then
     wave = champ_value
     exit
    endif
   endif


  END do

 return
end subroutine get_wavelength_from_cif_file
!--------------------------------------------------------------------------

subroutine get_crystal_system_from_CIF_file(file_unit, crystal_system)
 USE macros_module,   ONLY : l_case
 USE cryscalc_module, only  : debug_proc
 USE IO_module
 implicit none
  INTEGER,            INTENT(IN)       :: file_unit
  character (len=16), INTENT(OUT)      :: crystal_system
  character (len=256)                  :: line
  CHARACTER (LEN=256)                  :: champ_text, champ_value
  integer                              :: ier

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_CRYSTAL_SYSTEM_FROM_CIF_FILE")

  rewind(file_unit)
  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info(' ... No crystal system in the import.cif file !!')
    call write_info('')
    return
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading crystal system in import.cif file !!')
    call write_info('')
    return
   endif

   line = ADJUSTL(line)


   if(trim(line(1:22))=='_symmetry_cell_setting') then
    READ(line, *, IOSTAT=ier ) champ_text, champ_value
    IF(ier/=0) then
      call write_info('')
      call write_info(' ... Problem reading crystal_system !!')
      call write_info('')
      return
    endif

    exit
   endif
  END do
  champ_value = l_case(champ_value)

  if(champ_value(1:9) == 'triclinic')         then
   crystal_system = "TRIC"
  elseif(champ_value(1:10) == 'monoclinic')   then
   crystal_system = "MONO"
  elseif(champ_value(1:12) == 'orthorhombic') then
   crystal_system = "ORTHO"
  elseif(champ_value(1:10) == 'tetragonal')   then
   crystal_system = "TETRA"
  elseif(champ_value(1:8)  == 'trigonal')     then
   crystal_system = "TRIG"
  elseif(champ_value(1:9)  == 'hexagonal')    then
   crystal_system = "HEXA"
  elseif(champ_value(1:5)  == 'cubic')        then
   crystal_system = "CUB"
  endif


 return
end subroutine get_crystal_system_from_CIF_file


!--------------------------------------------------------------------------

subroutine get_H_M_from_CIF_file(file_unit, H_M)
 USE macros_module,   ONLY : l_case
 USE cryscalc_module, ONLY : write_details, debug_proc
 USE IO_module
 implicit none
  INTEGER,           INTENT(IN)        :: file_unit
  character (len=32), INTENT(OUT)      :: H_M
  character (len=256)                  :: line
  CHARACTER (LEN=256)                  :: champ_text, champ_value
  integer                              :: ier

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_H_M_FROM_CIF_FILE")

  rewind(file_unit)
  H_M = '?'     ! new 19.06.2015
  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    if(write_details) then
     call write_info('')
     call write_info(' ... No Hermann Mauguin space group in the import.cif file !!')
     call write_info('')
    end if
    return
   ELSEIF(ier >0) then
    if(write_details) then
     call write_info('')
     call write_info(' ... Problem reading crystal system in import.cif file !!')
     call write_info('')
    end if
    return
   endif

   line = ADJUSTL(line)


   if(trim(line(1:30))=='_symmetry_space_group_name_H-M') then
    READ(line, *, IOSTAT=ier ) champ_text, champ_value
    IF(ier/=0) then
      call write_info('')
      call write_info(' ... Problem reading crystal_system !!')
      call write_info('')
      return
    endif

    exit
   endif
  END do

  H_M = adjustl(champ_value(1:32))


 return
end subroutine get_H_M_from_CIF_file

!--------------------------------------------------------------------------

subroutine get_champ_value(CIF_unit, string_text, string_value, ok)
 use cryscalc_module, only : debug_proc
 use macros_module,   only : l_case
 implicit none
 INTEGER          , INTENT(IN)    :: CIF_unit
 CHARACTER (LEN=*), INTENT(IN)    :: string_text
 CHARACTER (LEN=*), INTENT(OUT)   :: string_value
 logical          , INTENT(INOUT) :: ok
 CHARACTER (LEN=256)              :: READ_line
 integer                          :: ier, long_string, long_line


 if(debug_proc%level_3)  call write_debug_proc_level(3, "GET_CHAMP_VALUE ("//trim(string_text)//")")



 REWIND(UNIT=CIF_unit)


 ok = .false.
 do
  READ(CIF_unit, '(a)', IOSTAT=ier) read_line
  if(ier/=0) return
  if(ier <0) return

  !if(index(read_line, '_shelx_res_file') /=0) return      ! evite de lire les fichiers .RES et .HKL (gain de temps)
  !if(index(read_line, '_shelx_hkl_file') /=0) return      ! evite de lire les fichiers .RES et .HKL (gain de temps)

  read_line   = ADJUSTL(read_line)
  long_line   = len_trim(read_line)
  !read_line   = l_case(read_line)
  if(long_line == 0) cycle
  long_string = LEN_TRIM(string_text)


  IF(l_case(READ_line(1:long_string)) == TRIM(l_case(string_text))) then
   if(long_line == long_string) then  ! la valeur du champ n'est pas sur la même ligne
    read(CIF_unit, '(a)', iostat=ier) read_line
    if(ier/=0) return
    read(read_line, '(a)') string_value
   else
     READ(read_line(long_string+1:), '(a)') string_value
   endif
   string_value = ADJUSTL(string_value)
   IF(LEN_TRIM(string_value) == 0)   return
   IF(string_value(1:1)      == '?') return
   ok = .true.
   return
  endif
 ENDDO


 return
end subroutine get_champ_value

!!-------------------------------------------------------------------!!

subroutine read_lines_before_hkl(file_unit, cos_exist)
 use cryscalc_module, only : debug_proc
 implicit none
  INTEGER, INTENT(IN)             :: file_unit
  LOGICAL, INTENT(OUT)            :: cos_exist
  ! local variables
  INTEGER                         :: i, num_ligne, ier
  CHARACTER(LEN=128)              :: line

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_LINES_BEFORE_HKL")

  num_ligne = 0
  rewind(UNIT=file_unit)

!  do
!   read(UNIT=file_unit, '(a)', iostat=ier) line
!   if(ier /=0) exit
!   num_ligne = num_ligne +1
!   if (line(1:22) =='_refln_F_squared_sigma' ) then
!    READ(UNIT=file_unit,'(a)',IOSTAT=Ier) (line, i=1,7)
!    if (line(1:34) == '_refln_nonius_diffracted_cos_cstar' .or. line(1:34) == '_nonius_refln_diffracted_cos_cstar') then
!     cos_exist = .true.
!     num_ligne = num_ligne + 7
!     exit
!    else
!     cos_exist = .false.
!     exit
!    end if
!   endif
!  end do

  cos_exist = .false.
loop_1:  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier /=0) exit
   num_ligne = num_ligne +1
   if (line(1:22) =='_refln_F_squared_sigma' ) then
    do
     READ(file_unit, '(a)', IOSTAT=ier) line
     IF(line(1:1)/='_') exit loop_1
     num_ligne = num_ligne + 1
     if (line == '_refln_nonius_diffracted_cos_cstar' .or. line == '_nonius_refln_diffracted_cos_cstar') cos_exist = .true.

    end do

   endif
  end do loop_1



  rewind(unit=file_unit)
  ! lecture des lignes d'entete
  read(file_unit, '(a)') (line, i=1,num_ligne)

 RETURN
END subroutine read_lines_before_hkl

!----------------------------------------------------------------------------
subroutine read_SQUEEZE_file(ID_sample)
 ! lecture du fichier au format CIF cree par la procedure SQUEEZE
 ! ecriture du contenu dans le fichier archive.CIF
 use cryscalc_module, only : SQUEEZE, include_SQUEEZE, debug_proc
 implicit none
  character (len=*) , intent(in) :: ID_sample
  logical               :: file_exist
  integer               :: i, i_error
  character (len=256)   :: read_line
  character (len=256), dimension(3) :: tmp_sqz_file

  tmp_sqz_file(1) = SQUEEZE%file
  tmp_sqz_file(2) = 'platon_sqr.sqf'
  tmp_sqz_file(3) = trim(ID_sample)//'.sqz'


 if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_SQUEEZE_FILE")

 !inquire(file= trim(SQUEEZE%file), exist= file_exist)
 SQUEEZE%procedure = .false.
 !if(.not. file_exist) then
 ! inquire(file = 'platon_sqr.sqf', exist=file_exist)
 ! if(.not. file_exist)  then
 !  inquire(file = trim(ID_sample)//'.sqz', exist=file_exist)
 !  if(.not. file_exist) return
 !  squeeze%file = trim(ID_sample)//'.sqz'
 ! else
 !  squeeze%file = 'platon_sqr.sqf'
 ! end if
 !end if

 file_exist = .false.
 do i=1, 3
  inquire(file= trim(tmp_sqz_file(i)), exist= file_exist)
  if(file_exist) exit
 end do

 if(.not. file_exist) return
 if(.not. include_SQUEEZE) return

 squeeze%file = trim(tmp_sqz_file(i))
 open(unit=SQUEEZE%unit, file=trim(SQUEEZE%file))
  do
   read(SQUEEZE%unit, '(a)', iostat= i_error) read_line
   if(i_error /=0) exit

   SQUEEZE%nb_lines = SQUEEZE%nb_lines + 1
   SQUEEZE%read_line(SQUEEZE%nb_lines) = read_line
  end do

 close (unit=SQUEEZE%unit)
 SQUEEZE%procedure = .true.

 call write_CIF_file("SQUEEZE")

 return
end subroutine read_SQUEEZE_file

! ----------------------------------------------------------------------------------
 subroutine read_and_modif_ARCHIVE(input_unit)
  use cryscalc_module, only : CIF_unit, CIF_archive_unit, AUTHOR, DEVICE, SQUEEZE, IMPORT_CIF_unit,    &
                              SADABS, CIFdep, ACTA, CIF_format80, include_RES_file, SPG,               &
                              include_HKL_file, tmp_unit, cryscalc, archive_cif, molecule, debug_proc
  use CIF_module
  use macros_module,   only : u_case, test_file_exist, Get_Wingx_job, Get_current_folder_name, Get_sample_ID_from_folder_name
  use IO_module,       only : write_info
  implicit none
   integer, intent(in)               :: input_unit
   character (len=256)               :: read_line
   integer                           :: i_error, i, j, i1, i2, long
   integer                           :: n_line
   logical                           :: file_exist
   logical                           :: diffracto
   logical                           :: correction_absorption
   logical                           :: squeeze_results
   logical                           :: author_data

   character (len=256)               :: CIF_string, tmp_string
   character (len=16), dimension(16) :: CIF_str
   character (len=64)                :: fmt_, fmt_final
   character (len=256)               :: job, sample_ID, diffraction_date
   character (len=256)               :: wingx_structure_dir, current_folder
   character (len=256)               :: wingx_res_file, wingx_hkl_file
   character (len=256)               :: CIF_sep_line
   logical                           :: contain_H
   logical                           :: sym_op_included

   if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_AND_MODIF_ARCHIVE")


   fmt_      = ''
   fmt_final = ''
   job       = "?"
   sample_ID = "?"
   write(CIF_sep_line, '(a,76a1,a)') '#', ('-',i=1,76), '#'

  rewind(unit=input_unit)
  close(unit = CIF_unit)
  open(unit = CIF_unit, file = 'cryscalc.cif')

  close(unit=CIF_archive_unit)
  if(include_HKL_file) then
   open(unit=CIF_archive_unit, file="archive_cryscalc_hkl.cif")
  else
   open(unit=CIF_archive_unit, file="archive_cryscalc.cif")
  end if

  ! premiere lecture du fichier archive.cif
  open (unit=input_unit, file=trim(archive_cif))  !! le fichier a été fermé je ne sais ou !
  n_line = 0
  do
   read(input_unit, '(a)', iostat= i_error) read_line
   if(i_error < 0) then
    n_line = 0
    exit
   endif
   if(i_error /=0) exit
   n_line  = n_line + 1
   i1 = index(read_line , "SECTION 2. COMPOUND(S) DETAILS")
   if(i1 /=0) exit
  end do


 ! >>> new sept. 09

  ! determination du nombre de champ dans la boucle _atom_site_
  call Get_CIF_champ_nb(input_unit, '_atom_site_label', CIF_atom_site_item_nb, CIF_atom_site_loop_numor)

  ! determination des positions des differents items dans la boucle
  call Get_CIF_numors(input_unit)

  ! determination du nombre de champ dans la boucle _atom_site_aniso
  call Get_CIF_champ_nb(input_unit, '_atom_site_aniso_label', CIF_atom_site_aniso_item_nb, CIF_atom_site_aniso_loop_numor)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! pour compatibilite de la routine avec le fichier archive.cif cree par WinGX
  ! a partir d'un fichier.CIF cree par SHELXL-2014 (les operateurs de symetrie ne sont
  ! pas inclus dans l'archive.cif)

  rewind(unit = input_unit)
  ! recherche de "_space_group_symop_operation_xyz
  sym_op_included = .false.
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit ! fin du fichier
   if(index(read_line, "_space_group_symop_operation_xyz") /=0 .or. &
      index(read_line, "_symmetry_equiv_pos_as_xyz") /=0) then
    sym_op_included = .true.
    exit
   end if
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! molecular moiety
  rewind(unit = input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   long = len_trim('_chemical_formula_moiety')
   if(index(read_line, '_chemical_formula_moiety') /=0) then
    read(read_line(long+1:),'(a)') CIF_parameter%formula_moiety
    CIF_parameter%formula_moiety = adjustl(CIF_parameter%formula_moiety)
    if(CIF_parameter%formula_moiety(1:1) == "'") CIF_parameter%formula_moiety = CIF_parameter%formula_moiety(2:)
    long = len_trim(CIF_parameter%formula_moiety)
    if(CIF_parameter%formula_moiety(long:long) == "'") CIF_parameter%formula_moiety = CIF_parameter%formula_moiety(1:long-1)
    exit
   end if
  end do

  ! molecular sum
  rewind(unit = input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   long = len_trim('_chemical_formula_sum')
   if(index(read_line, '_chemical_formula_sum') /=0) then
    read(read_line(long+1:),'(a)') CIF_parameter%formula_sum
    CIF_parameter%formula_sum = adjustl(CIF_parameter%formula_sum)
    if(CIF_parameter%formula_sum(1:1) == "'") CIF_parameter%formula_sum = CIF_parameter%formula_sum(2:)
    long = len_trim(CIF_parameter%formula_sum)
    if(CIF_parameter%formula_sum(long:long) == "'") CIF_parameter%formula_sum = CIF_parameter%formula_sum(1:long-1)
    exit
   end if

  end do



! >>> new sept. 09
!
! determination du _atom_site_refinement_flags pour les atomes d'hydrogene
  rewind (unit=input_unit)
  ! lecture des lignes avant la liste des atomes
  do i = 1, CIF_atom_site_loop_numor + CIF_atom_site_item_nb
   read(input_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit
  end do

  contain_H = .false.
  do
   read(input_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit
   if(len_trim(read_line) == 0) exit
   read(read_line, *) CIF_str(1:CIF_atom_site_item_nb)

   if(len_trim(CIF_str(CIF_atom_site_type_symbol_numor)) == 1 .and. &
               CIF_str(CIF_atom_site_type_symbol_numor)(1:1) == 'H') then
      contain_H = .true.
   else
   !  contain_H = .false.
      cycle
   endif

   if(len_trim(CIF_str(CIF_atom_site_type_symbol_numor)) == 1 .and. &
      CIF_str(CIF_atom_site_type_symbol_numor)(1:1) == 'H'    .and. &
      CIF_str(CIF_atom_site_calc_flag_numor)(1:1)   /= 'c')  then
      CIF_parameter%H_treatment = 'mixed'
   endif
  end do


!-------------------------------------------------------------------------------

  ! diffracto ?
  diffracto = .false.
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   if(index(read_line, '_diffrn_measurement_device_type')/=0) then
    diffracto = .true.
    exit
   endif
  end do

  ! correction d'absorption ?
  correction_absorption = .false.
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   if(index(read_line, '_exptl_absorpt_process_details')/=0) then
    correction_absorption = .true.
    exit
   endif
  end do


  ! squeeze  ?
  squeeze_results = .false.
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   if(index(read_line, 'SQUEEZE RESULTS')/=0) then
    squeeze_results = .true.
    exit
   endif
  end do


  ! author data  ?
  author_data = .false.
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   if(index(read_line, 'AUTHORS DETAILS')/=0) then
    author_data = .true.
    exit
   endif
  end do


  ! lecture à partir de la ligne contenant "SECTION 2. COMPOUND(S) DETAILS"
  rewind(unit=input_unit)
  if(n_line /=0) then
   do i=1, n_line
    read(input_unit, '(a)', iostat = i_error) read_line
   end do
  endif

 rewind(unit = input_unit)

  do
   read(input_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit

   if(.not. author_data) then
    if(index(read_line, '#  CIF files read  :' ) /=0) then
     call write_CIF(CIF_archive_unit, trim(read_line))
     if(CIFDEP) then
      call write_CIF_AUTHOR(CIF_archive_unit)
      if (ACTA) call write_CIF_text(CIF_archive_unit)
     endif
     cycle
    endif
   endif

  ! if(index(read_line, 'data_') /=0 .and. job=='?') then
  !  call Get_current_folder_name(current_folder)
  !  i1 = index(current_folder, '_')
  !  if(i1 /=0) then
  !   write(job, '(a)') current_folder(i1+1:)
  !  else
  !   call Get_Wingx_job(job, wingx_structure_dir)
  !  endif
  !  write(CIF_archive_unit, '(2a)') 'data_', trim(job)
 !   cycle
 !  endif


   if(index(read_line, 'data_') /=0 .and. sample_ID=='?') then
    call Get_current_folder_name(current_folder)
    call Get_sample_ID_from_folder_name(current_folder, sample_ID)

    i1 = index(current_folder, '\', back=.true.)
    if(i1 /=0) then
     i2 = index(current_folder(i1+1:),"_")
     if(i2 /=0) then
     write(sample_ID, '(a)') current_folder(i1+1:i1+i2-1)
    endif
   endif
   if(sample_ID == '?') call Get_Wingx_job(sample_ID, wingx_structure_dir)
   if(sample_ID == '?') sample_ID = trim(read_line(6:))

   diffraction_date = ''
   i1 = index(current_folder, '_', back=.true.)
   if (i1 /=0) then
    diffraction_date = current_folder(i1+1:)
    end if

    write(CIF_archive_unit, '(a)') trim(CIF_sep_line)
    long = len_trim(CRYSCALC%version)
    write(fmt_ , '(a,i2,a)') '(3a,', 27-long, 'x,a)'
    write(CIF_archive_unit, fmt=trim(fmt_))     '#  CIF file completed and formatted by CRYSCALC (', trim(CRYSCALC%version),')','#'
   !write(CIF_archive_unit, '(3a,19x,a)')      '#  CIF file completed and formatted by CRYSCALC (', trim(CRYSCALC%version),')','#'

    if(len_trim(diffraction_date) /= 0) then
     long = len_trim(diffraction_date)
     write(fmt_ , '(a,i2,a)') '(2a,', 41-long, 'x,a)'
     write(CIF_archive_unit, fmt = trim(fmt_))  '#  Date of diffraction experiment : ', trim(diffraction_date),"#"
     !write(CIF_archive_unit, '(2a,34x,a)')     '#  Date of diffraction experiment : ', trim(diffraction_date),"#"
    end if
    write(CIF_archive_unit, '(a)') trim(CIF_sep_line)
    write(CIF_archive_unit, '(a)')        ''

    write(CIF_archive_unit, '(a)') trim(CIF_sep_line)
    write(CIF_archive_unit, '(a)') trim(CIF_sep_line)
    write(CIF_archive_unit, '(a)')        ''
    if( sample_ID /='?') then
     write(CIF_archive_unit, '(2a)') 'data_', u_case(trim(sample_ID))
    else
     write(CIF_archive_unit, '(a)') trim(read_line)
    endif
    CIF_parameter%sample_ID = sample_ID
    cycle
   endif

   if(index(read_line, '_chemical_formula_moiety') /=0) then
    write(CIF_archive_unit, '(4a)') '_chemical_formula_moiety                ', "'",trim(CIF_parameter%formula_sum), "'"
    cycle
   end if
   if(index(read_line, '_chemical_formula_sum') /=0) then
    write(CIF_archive_unit, '(4a)') '_chemical_formula_sum                   ', "'",trim(CIF_parameter%formula_sum), "'"
    cycle
   end if

   ! ------------- oct. 2015. ------------------------------------------------
   if(index(read_line, '_cell_length_a') /=0 .and. .not. sym_op_included) then
    write(CIF_archive_unit, '(a)') 'loop_'
    write(CIF_archive_unit, '(a)') '_space_group_symop_operation_xyz'
    do j=1, SPG%Multip   !
     write(CIF_archive_unit, '(3a)') "'", trim(SPG%SymopSymb(j)), "'"
    end do

    write(CIF_archive_unit, '(a)') ' '
   end if
   ! -----------------------------------------------------------------------------

   if(index(read_line, '_exptl_absorpt_coefficient_mu') /=0 .and.  &
      .not. correction_absorption) then
      ! les info. sur les corrections d'absorption sont recuperees dans le fichier CRYSCALC.cif
    do
     read(CIF_unit, '(a)', iostat=i_error) read_line

     if(i_error /=0) exit
     if(index(read_line, '_exptl_absorpt_coefficient_mu') /=0) then

      write(CIF_archive_unit, '(a)') trim(read_line)
      !do i = 1,7
      do
       read(cif_unit, '(a)', iostat=i_error) read_line
       if(i_error /=0) exit
       if(read_line(1:11) == '#----------') exit
       write(CIF_archive_unit, '(a)') trim(read_line)
      end do
      exit
     endif
    end do

    cycle
   end if

   if(index(read_line, '_exptl_absorpt_correction_type') /=0) then
    if(.not. SQUEEZE_results .and. SQUEEZE%procedure) then
     write(CIF_archive_unit, '(a)') ""
     WRITE(CIF_archive_unit, '(a)') "#----------------------------------------------------------------------------#"
     WRITE(CIF_archive_unit, '(a)') "#                   SQUEEZE RESULTS                                          #"
     WRITE(CIF_archive_unit, '(a)') "#----------------------------------------------------------------------------#"
     write(CIF_archive_unit, '(a)') ""
     do i=2, SQUEEZE%nb_lines
      write(CIF_archive_unit, '(a)') trim(SQUEEZE%read_line(i))
     end do
     cycle
    else
     cycle
    endif
   endif


   !if(.not. squeeze_results) then
   ! if(index(read_line, '_exptl_absorpt_correction_type') /=0) then
   !  if(SQUEEZE%procedure)  then
   !   write(unit=CIF_archive_unit, '(a)') ""
   !   WRITE(UNIT=CIF_archive_unit, '(a)') "#----------------------------------------------------------------------------#"
   !   WRITE(unit=CIF_archive_unit, '(a)') "#                   SQUEEZE RESULTS                                          #"
   !   WRITE(UNIT=CIF_archive_unit, '(a)') "#----------------------------------------------------------------------------#"
   !   write(unit=CIF_archive_unit, '(a)') ""
   !   do i=2, SQUEEZE%nb_lines
   !    write(unit=CIF_archive_unit, '(a)') trim(SQUEEZE%read_line(i))
   !   end do
   !  endif
   !  cycle
   ! endif
   !else
   ! if(index(read_line, '_exptl_absorpt_correction_type') /=0) cycle
   !endif


   if(index(read_line, '_atom_sites_solution_hydrogens') /=0 .and. .not. contain_H) cycle

   if(index(read_line, '_refine_ls_hydrogen_treatment') /=0) then
    if(contain_H) then
     write(CIF_archive_unit, '(2a)') "_refine_ls_hydrogen_treatment           ",trim(CIF_parameter%H_treatment)
    endif
    cycle
   end if
   !if(index(read_line, 'mix') /=0) then
   ! write(CIF_archive_unit, '(a)') "_refine_ls_hydrogen_treatment           constr"
   ! cycle
   !endif


   if(index(read_line, '_diffrn_ambient_temperature')   /= 0    .and.  &
      .not. diffracto) then
    !if(u_case(DEVICE%diffracto(1:4)) == 'APEX' .or. u_case(DEVICE%diffracto(1:3)) == 'X2S') then
    if(DEVICE%APEX .or. DEVICE%X2S) then
     write(CIF_archive_unit, '(2a)') "_diffrn_measurement_device_type         ", &
                                      trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
     write(CIF_archive_unit, '(2a)') "_diffrn_measurement_method              ", &
                                      trim(CIF_parameter_DEVICE%diffrn_measurement_method)
     write(CIF_archive_unit, '(2a)') "_diffrn_detector                        ", &
                                      trim(CIF_parameter_DEVICE%diffrn_detector)
     write(CIF_archive_unit, '(a)') trim(read_line)
    !elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
    elseif(DEVICE%KCCD) then
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_source                          ", trim(CIF_parameter_DEVICE%diffrn_source)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_detector                        ", trim(CIF_parameter_DEVICE%diffrn_detector)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_detector_area_resol_mean        ", &
                                     trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_measurement_device              ", &
                                     trim(CIF_parameter_DEVICE%diffrn_measurement_device)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_measurement_device_type         ", &
                                     trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_measurement_method              ", trim(CIF_parameter_DEVICE%diffrn_measurement_method)
     WRITE(CIF_archive_unit, '(a)') trim(read_line)
    !elseif(u_case(DEVICE%diffracto(1:10)) == 'D8 VENTURE' .or. u_case(DEVICE%diffracto(1:10)) == 'D8_VENTURE' .or. &
    !       u_case(DEVICE%diffracto(1:3))  == 'D8V') then
    elseif(DEVICE%D8V) then
     write(CIF_archive_unit, '(2a)') "_diffrn_measurement_device_type         ", &
                                      trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
     write(CIF_archive_unit, '(2a)') "_diffrn_measurement_method              ", &
                                      trim(CIF_parameter_DEVICE%diffrn_measurement_method)
     write(CIF_archive_unit, '(2a)') "_diffrn_detector                        ", &
                                      trim(CIF_parameter_DEVICE%diffrn_detector)
     WRITE(CIF_archive_unit, '(2a)') "_diffrn_radiation_source                ", &
                                      trim(CIF_parameter_DEVICE%diffrn_radiation_source)
     WRITE(CIF_archive_unit, '(2a)') "_diffrn_radiation_monochromator         ", &
                                      trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator)
     WRITE(CIF_archive_unit, '(a)') trim(read_line)
    endif
    cycle
   endif

   !if(index(read_line, '_diffrn_radiation_monochromator')   /= 0  .and. u_case(DEVICE%diffracto(1:10)) == 'D8 VENTURE') then
   if(index(read_line, '_diffrn_radiation_monochromator')   /= 0  .and. DEVICE%D8V) then
    cycle
   endif


   if(index(read_line, '_computing_structure_solution') /=0) then
    write(CIF_archive_unit, '(2a)') "_computing_structure_solution           ", &
                                     trim(CIF_parameter_DEVICE%computing_structure_solution)
    cycle
   endif

   if(index(read_line, '_computing_publication_material') /=0) then
    write(CIF_archive_unit, '(a)')     "_computing_publication_material          "
    WRITE(CIF_archive_unit, '(a)')     ";"
    WRITE(CIF_archive_unit, '(10x,a)') trim(CIF_parameter_DEVICE%computing_publication_material_1)
    WRITE(CIF_archive_unit, '(10x,a)') trim(CIF_parameter_DEVICE%computing_publication_material_2)
    WRITE(CIF_archive_unit, '(a)')     ";"
    if(len_trim(adjustl(read_line)) == len_trim('_computing_publication_material')) then
     read(input_unit, '(a)', iostat = i_error) read_line
     if(i_error /=0) exit
    endif

    cycle
   endif

   ! -------  juillet 2011 : add jobname.res
   if(index(read_line, '_refine_diff_density_rms') /= 0) then
    write(CIF_archive_unit, '(a)')   trim(read_line)
    write(CIF_archive_unit, '(a)')   ''
    if(include_RES_file) then
     call Get_Wingx_job(job, wingx_structure_dir)
     if(job(1:1) /= '?') then
      write(wingx_res_file, '(4a)') trim(wingx_structure_dir), '\', trim(job), '.res'
      call test_file_exist(trim(wingx_res_file), file_exist, 'out')
      if(file_exist) then
       call write_info('')
       call write_info('   ... Include '//trim(wingx_res_file)//' SHELXL file.')
       call write_info('')
       write(CIF_archive_unit, '(a)')   '_iucr_refine_instructions_details'
       write(CIF_archive_unit, '(a)')   ';'
       open(unit=tmp_unit, file=trim(wingx_res_file))
       write(CIF_archive_unit, '(2a)') ' .res SHELXL file: ', trim(job)//'.res'
       write(CIF_archive_unit, '(a)')  '.................................................................'
       do
        read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
        if(i_error /=0) exit
        write(CIF_archive_unit, '(a)') trim(read_line)
       end do
       close(unit=tmp_unit)
       write(CIF_archive_unit, '(a)')   ';'
       !write(CIF_archive_unit, '(a)')   '_shelx_res_checksum              ?'
       write(CIF_archive_unit, '(a)')   '# End of res file'
        ! -------  juillet 2014 : add jobname.hkl
       if(include_HKL_file) then
        ! include .hkl file
         write(wingx_hkl_file, '(4a)') trim(wingx_structure_dir),'\', trim(job), '.hkl'
        call test_file_exist(trim(wingx_hkl_file), file_exist, 'out')
        if(file_exist) then
         call write_info('')
         call write_info('   ... Include '//trim(wingx_hkl_file)//' SHELXL file.')
         call write_info('')
         write(CIF_archive_unit, '(a)') ''
         !write(CIF_archive_unit, '(a)') '_shelx_hkl_file'
         write(CIF_archive_unit, '(a)') '_iucr_refine_reflections_details'
         write(CIF_archive_unit, '(a)') ';'
         open(unit=tmp_unit, file=trim(wingx_hkl_file))
          write(CIF_archive_unit, '(2a)') ' .hkl SHELXL file: ', trim(job)//'.hkl'
          write(CIF_archive_unit, '(a)')  '.................................................................'
         do
          read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
          if(i_error /=0) exit
          write(CIF_archive_unit, '(a)') trim(read_line)
         end do
         close(unit=tmp_unit)
         write(CIF_archive_unit, '(a)')   ';'
         !write(CIF_archive_unit, '(a)')   '_shelx_hkl_checksum              ?'
         write(CIF_archive_unit, '(a)')   '# End of hkl file'
        endif

       end if
       write(CIF_archive_unit, '(a)') trim(CIF_sep_line)
      endif
     endif
    end if
    cycle
   endif
   ! -----------------------------------------------------------

   if(index(read_line, '_atom_type_scat_source') /=0)  then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:4)
     i1 = index(read_line, "'")
     read(read_line(i1:), '(a)') CIF_string
     if(CIF_str(3)(1:1) /= '-') then
      write(CIF_archive_unit, '(2a4, 1x,2a8, a)') CIF_str(1:4), trim(CIF_string)
     else
      write(CIF_archive_unit, '(2a4, a8,1x,a8, a)') CIF_str(1:4), trim(CIF_string)
     endif
    end do
    write(CIF_archive_unit, '(a)') ''

   endif


  ! if(index(read_line, '_atom_type_scat_source') /=0)  then
  !  write(CIF_archive_unit, '(a)') trim(read_line)
!    do
!     read(input_unit, '(a)', iostat=i_error) read_line
!     if(i_error /=0) exit
!     if(len_trim(read_line) == 0) exit
!     write(CIF_archive_unit, '(a)') trim(read_line)
!    end do
!    write(CIF_archive_unit, '(a)') ''
!
!   endif




!  >>> new sept. 09

   if(index(read_line, "_atom_site_label") /=0  .and. &
     len_trim(adjustl(read_line)) ==len_trim("_atom_site_label")) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do i=1, CIF_atom_site_item_nb - 1
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     write(CIF_archive_unit, '(a)') trim(read_line)
    end do

    ! boucle sur les atomes et ecriture (formattee si demandé) des positions atomiques
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:CIF_atom_site_item_nb)

     if(CIF_atom_site_refinement_flags_numor /=0) then
      long = len_trim(CIF_str(CIF_atom_site_refinement_flags_numor))
      if(index(CIF_str(CIF_atom_site_refinement_flags_numor), 'P') /=0) then ! occ. partielle
       call get_fmt(4, CIF_str(CIF_atom_site_fract_x_numor : CIF_atom_site_fract_x_numor +3), fmt_, 0)
       if(CIF_atom_site_type_symbol_numor == 2) then  ! symbol en position 2
        if(CIF_atom_site_item_nb == 12) then
         write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a10,a5,a4,2a2)'
        else
         if(len_trim(CIF_str(CIF_atom_site_symm_multiplicity_numor)) == 1) then
          write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a10,a2,a5,a4,2a2)'
         else
          write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a10,a3,a5,a4,2a2)'
         end if
        endif
       elseif(CIF_atom_site_type_symbol_numor == CIF_atom_site_item_nb) then ! _symbol en derniere position
        if(CIF_atom_site_item_nb == 12) then
         write(fmt_final, '(3a)') '(a5,', trim(fmt_), ',a5,a10,a5,a4,2a2,a3)'
        else
         if(len_trim(CIF_str(CIF_atom_site_symm_multiplicity_numor)) == 1) then
          write(fmt_final, '(3a)') '(a5,', trim(fmt_), ',a5,a10,a2,a5,a4,2a2,a3)'
         else
          write(fmt_final, '(3a)') '(a5,', trim(fmt_), ',a5,a10,a3,a5,a4,2a2,a3)'
         endif
        endif
       endif

      else
       call get_fmt(4, CIF_str(CIF_atom_site_fract_x_numor : CIF_atom_site_fract_x_numor +3), fmt_, 0)
       if(CIF_atom_site_type_symbol_numor == 2) then  ! symbol en position 2
        if(CIF_atom_site_item_nb==12) then
         write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2,a5,a3,2a2)'
        else
         if(len_trim(CIF_str(CIF_atom_site_symm_multiplicity_numor)) == 1) then
          write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2,a2,a5,a3,2a2)'
         else
          write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2,a3,a5,a3,2a2)'
         end if
        endif
       elseif(CIF_atom_site_type_symbol_numor == CIF_atom_site_item_nb) then ! _symbol en derniere position

       if(CIF_atom_site_item_nb==12) then
        write(fmt_final, '(3a)') '(a5,', trim(fmt_), 'a5,a2,a5,a3,2a2,a3)'
       else
        if(len_trim(CIF_str(CIF_atom_site_symm_multiplicity_numor)) == 1) then
         write(fmt_final, '(3a)')  '(a5,', trim(fmt_), 'a5,a2,a2,a5,a3,2a2,a3)'
        else
         write(fmt_final, '(3a)')  '(a5,', trim(fmt_), 'a5,a2,a3,a5,a3,2a2,a3)'
        endif
       endif
      endif
     end if

    else  ! CIF_atom_site_refinement_flags_numor /= 0  (SHELX-2014)
      call get_fmt(4, CIF_str(CIF_atom_site_fract_x_numor : CIF_atom_site_fract_x_numor +3), fmt_, 0)
      if(CIF_atom_site_type_symbol_numor == 2) then  ! symbol en position 2
       if(CIF_atom_site_item_nb==11) then
        if(len_trim(CIF_str(CIF_atom_site_occupancy_numor)) == 1) then     ! occupation = 1
         write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2,a5,a3,a2)'
        else
         write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a10,a5,a3,a2)'
        endif
       endif

      elseif(CIF_atom_site_type_symbol_numor == CIF_atom_site_item_nb) then ! _symbol en derniere position

       if(CIF_atom_site_item_nb==11) then
        if(len_trim(CIF_str(CIF_atom_site_occupancy_numor)) == 1) then     ! occupation = 1
         write(fmt_final, '(3a)') '(a5,', trim(fmt_), 'a5,a2,a5,a3,a2,a3)'
        else
         write(fmt_final, '(3a)') '(a5,', trim(fmt_), 'a5,a7,a5,a3,a2,a3)'
        endif
       endif
      endif
     end if   ! CIF_atom_site_refinement_flags_numor /= 0

     write(tmp_string, fmt_final) CIF_str(1:CIF_atom_site_item_nb)
     if(len_trim(tmp_string) > 80 .and. .not. CIF_format80) then
      write(CIF_archive_unit, '(a)') trim(read_line)
     else
      write(CIF_archive_unit, fmt_final) CIF_str(1:CIF_atom_site_item_nb)
     endif


     !if(CIF_str(12)(1:1) == '.') then
     ! call get_fmt(4, CIF_str(3:6), fmt_, 0)
     ! write(fmt_, '(3a)') '(a5,a3', trim(fmt_), ',a5,a2,a5,3a3)'
     ! write(CIF_archive_unit, fmt_) CIF_str(1:12)
     !else ! cas d'un desordre traite avec PART 1 / PÄRT 2
     ! write(CIF_archive_unit, '(a)') trim(read_line)
     !endif
    end do
    write(CIF_archive_unit, '(a)') ''

   endif


!  >> new sept. 09

    if(index(read_line, "_atom_site_aniso_label") /=0) then
     write(CIF_archive_unit, '(a)') trim(read_line)
     do i=1, CIF_atom_site_aniso_item_nb - 1
      read(input_unit, '(a)', iostat=i_error) read_line
      if(i_error /=0) exit
      write(CIF_archive_unit, '(a)') trim(read_line)
     end do

     do
      read(input_unit, '(a)', iostat=i_error) read_line
      if(i_error /=0) exit
      if(len_trim(read_line) == 0) exit
      read_line = adjustl(read_line)
      read(read_line, *) CIF_str(1:7)
      call get_fmt(6, CIF_str(2:7), fmt_, 0)

      write(fmt_final, '(3a)') '(a4,', trim(fmt_),')'
      write(tmp_string, fmt_final) CIF_str(1:7)
      if(len_trim(tmp_string) > 80 .and. .not. CIF_format80) then
       write(CIF_archive_unit, '(a)') trim(read_line)
      else
       write(CIF_archive_unit, fmt_final) CIF_str(1:7)
      endif
     end do
     write(CIF_archive_unit, '(a)') ''

    endif

!-----------------------------------------------------------------------------------------------


   if(index(read_line, '_geom_bond_publ_flag') /=0) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:5)
     write(CIF_archive_unit, '(2a5, a12,2a12)') CIF_str(1:5)
    end do
    write(CIF_archive_unit, '(a)') ''
   endif

   if(index(read_line, '_geom_angle_publ_flag') /=0) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:7)
     call get_fmt(1, CIF_str(4:4), fmt_, 100)
     write(fmt_final, '(3a)') '(3a5,', trim(fmt_), ',3a8)'
     write(CIF_archive_unit, fmt_final) CIF_str(1:7)
    end do
    write(CIF_archive_unit, '(a)') ''
   endif


   if(index(read_line, '_geom_torsion_publ_flag') /=0) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:10)

     !call verif_linearity(CIF_str(1:4), linear)  ! calculs.F90
     !if(linear) cycle

     call get_fmt(1, CIF_str(5:5), fmt_, -100)
     write(fmt_final, '(3a)') '(4a5', trim(fmt_), ',x,5a8)'
     write(CIF_archive_unit, fmt_final) CIF_str(1:10)
    end do
    write(CIF_archive_unit, '(a)') ''
   endif

   if(index(read_line, '_geom_hbond_site_symmetry_A') /=0) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:8)
     !call get_fmt(1, CIF_str(5:5), fmt_, -100)
     !write(fmt_final, '(3a)') '(4a5', trim(fmt_), ',5a6)'
     !write(CIF_archive_unit, fmt_) CIF_str(1:8)
     write(CIF_archive_unit, '(3a5, 4a12,a6)') CIF_str(1:8)
    end do
    write(CIF_archive_unit, '(a)') ''

   endif


   write(CIF_archive_unit, '(a)') trim(read_line)
  end do

   if(sample_ID /= '?') then
    write(CIF_archive_unit, '(a)') trim(CIF_sep_line)
    long = len_trim(sample_ID)
    write(fmt_, '(a,i2,a)') '(2a,', 53-long, 'a1, a)'

    write(CIF_archive_unit, trim(fmt_)) '#===END OF CIF FILE FOR ', u_case(trim(sample_ID)), (' ',i=1,53-long), '#'
    !write(CIF_archive_unit, '(2a)') '# END OF CIF FILE FOR ', u_case(trim(sample_ID))
   endif
   write(CIF_archive_unit, '(a)') trim(CIF_sep_line)
  close(unit=CIF_unit)

  return
 end subroutine read_and_modif_ARCHIVE


!---------------------------------------------------------------------------
 subroutine Get_fmt(n, string, fmt_, level)
  implicit none
  integer,            intent(in)               :: n
  character (len=*),  dimension(n), intent(in) :: string
  character (len=64), intent(out)              :: fmt_
  integer,            intent(in)               :: level
  real                                         :: value
  integer                                      :: i, i1, len_max, long_max

  fmt_ = ''
  if(level == 0) then
   long_max = 11
  else
   long_max = 11
  endif

  !len_max = len_trim(string(1))
  !if(n > 2) then
  !do i=2, n
  ! if(len_trim(string(i)) > len_max) len_max = len_trim(string(i))
  !end do
  !if(len_max < long_max) len_max = long_max
  !len_max = len_max + 1
  !end if
  len_max = long_max + 1


  do i = 1, n
   i1 = index(string(i), '(')
   if(i1 ==0) then
    read(string(i), *) value
   else
    read(string(i)(1:i1-1), *) value
   endif

   if(level == 0) then
    if(string(i)(1:1) == "-") then
     !write(fmt_, '(2a)') trim(fmt_), ',a13'
     write(fmt_, '(2a,i2)') trim(fmt_), ',a', len_max+1
    else
     if(i==1) then
      !write(fmt_, '(2a)') trim(fmt_), '1x,a12'
      write(fmt_, '(2a,i2)') trim(fmt_), '1x,a', len_max
     else
      !write(fmt_, '(2a)') trim(fmt_), ',1x,a12'
      write(fmt_, '(2a,i2)') trim(fmt_), ',1x,a', len_max
     endif
    endif


   elseif(level == 100) then
    if(value > 99.99) then
     !write(fmt_, '(2a)') trim(fmt_), ',a13'
     write(fmt_, '(2a,i2)') trim(fmt_), ',a',len_max+1
    else
     if(i==1) then
      !write(fmt_, '(2a)') trim(fmt_), '1x,a12'
      write(fmt_, '(2a,i2)') trim(fmt_), '1x,a',len_max
     else
      !write(fmt_, '(2a)') trim(fmt_), ',1x,a12'
      write(fmt_, '(2a,i2)') trim(fmt_), ',1x,a', len_max
     endif
    endif

   elseif(level == -100) then
    if(value < -99.99) then
     !write(fmt_, '(2a)') trim(fmt_), ',a13'
     write(fmt_, '(2a,i2)') trim(fmt_), ',a', len_max+1
    elseif(value < -9.99) then
     !write(fmt_, '(2a)') trim(fmt_), ',1x,a12'
     write(fmt_, '(2a,i2)') trim(fmt_), ',1x,a', len_max
    elseif(string(i)(1:1) == '-') then
     !write(fmt_, '(2a)') trim(fmt_), ',2x,a11'
     write(fmt_, '(2a,i2)') trim(fmt_), ',2x,a', len_max-1
    elseif(value < 9.99) then
     !write(fmt_, '(2a)') trim(fmt_), ',3x,a10'
     write(fmt_, '(2a,i2)') trim(fmt_), ',3x,a', len_max-2
    elseif(value < 99.99) then
     !write(fmt_, '(2a)') trim(fmt_), ',2x,a11'
     write(fmt_, '(2a,i2)') trim(fmt_), ',2x,a', len_max-1
    else
     !write(fmt_, '(2a)') trim(fmt_), ',1x,a12'
     write(fmt_, '(2a,i2)') trim(fmt_), ',1x,a', len_max
    endif
   endif
  end do

  if(fmt_(1:1) == ',') fmt_ = trim(fmt_(2:))

  return
 end subroutine Get_fmt

!---------------------------------------------------------------------------
 subroutine Get_SADABS_ratio
  use cryscalc_module, only : SADABS, IMPORT_cif_unit, debug_proc
  use macros_module,   only : test_file_exist
  implicit none
   integer                  :: i_error,i1
   character (len=256)      :: read_line
   logical                  :: file_exist

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_SADABS_RATIO")

  ! lecture du fichier import.cif pour recuperer la valeur de Tmin/Tmax
  close(unit=IMPORT_CIF_unit)
  open(unit = IMPORT_CIF_unit, file = 'import.cif', iostat=i_error)
  call test_file_exist("import.cif", file_exist, 'out')
  do
   read(IMPORT_CIF_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit
   read_line = adjustl(read_line)
   i1 = index(read_line, '# Ratio of minimum to maximum apparent transmission:')
   if (i1 /=0) then
    i1 = index(read_line, ':')
    read(read_line(i1+1:), *) SADABS%ratio
    if(SADABS%ratio < 0. .or. SADABS%ratio > 1.) SADABS%ratio = -999.
   endif
  end do
  close(unit=IMPORT_CIF_unit)

  if(SADABS%ratio < 0.) then
   open(unit = IMPORT_CIF_unit, file = 'import.cif', iostat=i_error)
   call test_file_exist("import.cif", file_exist, 'out')
   do
    read(IMPORT_CIF_unit, '(a)', iostat = i_error) read_line
    if(i_error /=0) exit
    read_line = adjustl(read_line)
    i1 = index(read_line, '# Estimated minimum and maximum transmission:')
    if (i1 /=0) then
     i1 = index(read_line, ':')
     read(read_line(i1+1:), *) SADABS%Tmin, SADABS%Tmax
     SADABS%ratio = SADABS%Tmin / SADABS%Tmax
    endif
   end do
   close(unit=IMPORT_CIF_unit)
  endif



 return
end subroutine Get_SADABS_ratio

!----------------------------------------------------------------------------------
! recheche de la position d'un champ particulier dans une boucle loop_ du fichier.CIF

subroutine  Get_CIF_champ_numor(input_unit, input_string, output_numor)
 use cryscalc_module, only : debug_proc
 implicit none
  integer,           intent(in)    :: input_unit
  character (len=*), intent(in)    :: input_string
  integer,           intent(out)   :: output_numor
  CHARACTER (LEN=256)              :: READ_line
  integer                          :: i_error

  if(debug_proc%level_3)  call write_debug_proc_level(3, "GET_CIF_CHAMP_NUMOR ("//trim(input_string)//")")

  output_numor = 0
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error /=0) return
   if(index(read_line, trim(input_string)) /=0) then

    do
     backspace(unit=input_unit)
     backspace(unit=input_unit)
     read(input_unit, '(a)', iostat=i_error) read_line

     if(i_error /=0) return
     read_line = adjustl(read_line)
     if(read_line(1:5) == 'loop_') exit
    end do
    do
     output_numor = output_numor + 1
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) return
     read_line = adjustl(read_line)
     if(index(read_line, trim(input_string)) /=0) return  ! exit
     !output_numor = output_numor + 1
    end do

   endif
  end do

   return
  end subroutine Get_CIF_champ_numor

!--------------------------------------------------------------------
subroutine Get_CIF_champ_nb(input_unit, input_string, output_nb, output_nb2)
  use cryscalc_module, only : debug_proc
   implicit none
   integer          , intent(in)    :: input_unit
   character (len=*), intent(in)    :: input_string
   integer,           intent(out)   :: output_nb, output_nb2
   CHARACTER (LEN=256)              :: READ_line
   integer                          :: i_error

   if(debug_proc%level_3)  call write_debug_proc_level(3, "GET_CIF_CHAMP_NB ("//trim(input_string)//")")

   output_nb  = 0
   output_nb2 = 0
   rewind(unit = input_unit)
   do

    read(input_unit, '(a)', iostat=i_error) read_line
    if(i_error /=0) return
    output_nb2 = output_nb2 + 1
    read_line = adjustl(read_line)
    if(index(read_line, trim(input_string)) == 0) cycle

    do
     backspace(unit=input_unit)
     backspace(unit=input_unit)
     output_nb2 = output_nb2 - 1
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) return
     read_line = adjustl(read_line)
     if(read_line(1:5) == 'loop_') exit
    end do

    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) return
     read_line = adjustl(read_line)
     if(read_line(1:11) /= "_atom_site_") return
     output_nb = output_nb + 1
    end do

   end do

   return
end subroutine Get_CIF_champ_nb

!-------------------------------------------------------------------------------------------------------
 subroutine determine_H_treatment(input_unit, ls_H_treatment)
  USE CIF_module
  implicit none
  integer            , intent(in)    :: input_unit
  logical                            :: contain_H
  character (len=256), intent(out)   :: ls_H_treatment
  !integer                            :: CIF_atom_site_item_nb                ! nombre de champs dans la boucle _atom_site_
  !integer                            :: CIF_atom_site_label_numor
  !integer                            :: CIF_atom_site_type_symbol_numor      ! numero du champ
  !integer                            :: CIF_atom_site_loop_numor             ! numero de la ligne du champ _loop
  !integer                            :: CIF_atom_site_calc_flag_numor        ! numero du champ
  character (len=256)                :: read_line
  character (len=256), dimension(15) :: CIF_str
  integer                            :: i, i_error

! determination du nombre de champ dans la boucle _atom_site_
  call Get_CIF_champ_nb(input_unit, '_atom_site_label', CIF_atom_site_item_nb, CIF_atom_site_loop_numor)
   if(CIF_atom_site_item_nb == 0) return

! champ _atom_site_type_symbol
  call Get_CIF_champ_numor(input_unit, '_atom_site_type_symbol', CIF_atom_site_type_symbol_numor)

! champ _atom_site_calc_flag
  call Get_CIF_champ_numor(input_unit, '_atom_site_calc_flag', CIF_atom_site_calc_flag_numor)


! determination du _atom_site_refinement_flags pour les atomes d'hydrogene
  rewind (unit=input_unit)
  ! lecture des lignes avant la liste des atomes
  do i = 1, CIF_atom_site_loop_numor + CIF_atom_site_item_nb
   read(input_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit
  end do

  contain_H = .false.
  ls_H_treatment = 'constr'
  do
   read(input_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit
   if(len_trim(read_line) == 0) exit
   read_line = adjustl(read_line)
   if(read_line(1:1) == "#") cycle
   read(read_line, *, iostat=i_error) CIF_str(1:CIF_atom_site_item_nb)
   if(i_error /=0) exit

   if(len_trim(CIF_str(CIF_atom_site_type_symbol_numor)) == 1 .and. &
               CIF_str(CIF_atom_site_type_symbol_numor)(1:1) == 'H') then
     contain_H = .true.
     if (CIF_str(CIF_atom_site_calc_flag_numor)(1:1)   /= 'c')  then
      ls_H_treatment = 'mixed'
      exit
     end if
   else
   !  contain_H = .false.
      cycle
   endif

   !if(len_trim(CIF_str(CIF_atom_site_type_symbol_numor)) == 1 .and. &
   !   CIF_str(CIF_atom_site_type_symbol_numor)(1:1) == 'H'    .and. &
   !   CIF_str(CIF_atom_site_calc_flag_numor)(1:1)   /= 'c')  then
   !   ls_H_treatment = 'mixed'
!      exit
!   endif
  end do

  return
 end subroutine determine_H_treatment

!-------------------------------------------------------------------------------------------------------

subroutine Create_CEL_from_CIF
! creation d'un fichier .CEL à partir d'un fichier.CIF (ou .INS)

 USE cryscalc_module, only : keyword_read_CIF, keyword_read_PCR, keyword_read_INS, message_text,    &
                             CEL_unit, INS_file_name, CIF_file_name, CEL_file_name, PCR_file_name,   &
                             unit_cell, SPG, nb_atom, atom_coord, atom_typ, atom_occ, atom_occ_perc, &
                             debug_proc
 USE IO_module,                       ONLY : write_info
 USE CFML_Scattering_Chemical_Tables, ONLY : chem_info, num_chem_info,   set_chem_info
 USE macros_module,                   ONLY : u_case

 implicit none
  integer                 :: i, j

 if(debug_proc%level_2)  write(debug_proc%unit, '(a)') " . routine : CREATE_CEL_FROM_CIF"

  if (.not. keyword_read_CIF .and. .not. keyword_read_INS .and. .not. keyword_read_PCR) then
   call write_info('')
   if(.not. keyword_read_CIF) then
    call write_info(' ! A CIF file has to be read first !')
   elseif(.not. keyword_read_INS) then
    call write_info(' ! A INS file has to be read first !')
   elseif(.not. keyword_read_PCR) then
    call write_info(' ! A PCR file has to be read first !')
   endif
   call write_info('')
   return
  end if

  if(keyword_read_CIF) then
   i = index(CIF_file_name, '.', back=.true.)
   if (i/=0) then
    write(CEL_file_name, '(2a)') CIF_file_name(1:i-1), '_cif.CEL'
   else
    write(CEL_file_name, '(2a)') trim(CIF_file_name), '_cif.CEL'
   endif
  elseif(keyword_read_INS) then
   i = index(INS_file_name, '.', back=.true.)
   if (i/=0) then
    write(CEL_file_name, '(2a)') INS_file_name(1:i-1), '_ins.CEL'
   else
    write(CEL_file_name, '(2a)') trim(INS_file_name), '_ins.CEL'
   endif
  elseif(keyword_read_PCR) then
   i = index(PCR_file_name, '.', back=.true.)
   if (i/=0) then
    write(CEL_file_name, '(2a)') PCR_file_name(1:i-1), '_pcr.CEL'
   else
    write(CEL_file_name, '(2a)') trim(PCR_file_name), '_pcr.CEL'
   endif
  endif

  CALL set_chem_info

  open(unit = CEL_unit, file = trim(CEL_file_name))
    write(message_text, '(a,6F10.4)')   'CELL ', (unit_cell%param(i), i=1,6)
    write(CEL_unit, '(a)') adjustl(trim(message_text))
    if(nb_atom > 1) write(CEL_unit, '(a,i4)')     'natom ', nb_atom
    do i=1, nb_atom
     do j=1, Num_Chem_Info
       if (u_case(atom_typ(i)(1:))== u_case(chem_info(j)%symb(1:))) exit
     end do
     write(message_text,  '(a2,2x, i3,4F10.5)')  adjustl(atom_typ(i)), chem_info(j)%Z , atom_coord(1:3,i), atom_occ_perc(i)
     write(CEL_unit, '(a)') adjustl(trim(message_text))
    end do
    write(message_text, '(a,i3)')     'rgnr ',  SPG%NumSpg
    write(CEL_unit, '(a)') adjustl(trim(message_text))
  close(unit = CEL_unit)

  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(CEL_file_name), ' file has been created.'
  call write_info(trim(message_text))
  call write_info('')

  return
end subroutine create_CEL_from_CIF

!-------------------------------------------------------------------------------------------------------

subroutine Create_ACE_from_CIF
! creation d'un fichier .ACE pour CaRIne à partir d'un fichier.CIF (ou .INS)

 USE cryscalc_module, only : keyword_read_CIF, keyword_read_PCR, keyword_read_INS, message_text,   &
                             ACE_unit, INS_file_name, CIF_file_name, ACE_file_name, PCR_file_name,  &
                             unit_cell, SPG, nb_atom, atom_coord, atom_typ, atom_occ, atom_occ_perc, &
                             debug_proc
 USE IO_module,       ONLY : write_info
 USE CFML_Scattering_Chemical_Tables, ONLY : chem_info, num_chem_info,   set_chem_info
 USE macros_module,                   ONLY : u_case

  implicit none
  integer                           :: i
  character (len=32), dimension(6)  :: cell_string
  character (len=64)                :: ACE_string
  character (len=32), dimension(3)  :: coord_string
  character (len=1)                 :: tab
  real, parameter                   :: eps = 0.00001

  if(debug_proc%level_2)  call write_debug_proc_level(2, "CREATE_ACE_FROM_CIF")


  tab = char(9)  ! tabulation

  if (.not. keyword_read_CIF .and. .not. keyword_read_INS .and. .not. keyword_read_PCR) then
   call write_info('')
   if(.not. keyword_read_CIF) then
    call write_info(' ! A CIF file has to be read first !')
   elseif(.not. keyword_read_INS) then
    call write_info(' ! A INS file has to be read first !')
   elseif(.not. keyword_read_PCR) then
    call write_info(' ! A PCR file has to be read first !')
   endif
   call write_info('')
   return
  end if

  if(keyword_read_CIF) then
   i = index(CIF_file_name, '.', back=.true.)
   if (i/=0) then
    write(ACE_file_name, '(2a)') CIF_file_name(1:i-1), '_cif.ACE'
   else
    write(ACE_file_name, '(2a)') trim(CIF_file_name), '_cif.ACE'
   endif
  elseif(keyword_read_INS) then
   i = index(INS_file_name, '.', back=.true.)
   if (i/=0) then
    write(ACE_file_name, '(2a)') INS_file_name(1:i-1), '_ins.ACE'
   else
    write(ACE_file_name, '(2a)') trim(INS_file_name), '_ins.ACE'
   endif
  elseif(keyword_read_PCR) then
   i = index(PCR_file_name, '.', back=.true.)
   if (i/=0) then
    write(ACE_file_name, '(2a)') PCR_file_name(1:i-1), '_pcr.ACE'
   else
    write(ACE_file_name, '(2a)') trim(PCR_file_name), '_pcr.ACE'
   endif
  endif

  CALL set_chem_info



  open(unit = ACE_unit, file = trim(ACE_file_name))
   write(ACE_unit, '(a)')  'CaRIne Crystallography 3.0'
   write(ACE_unit, '(2a)') 'Cell Standard ASCII File=', trim(ACE_file_name)
   write(ACE_unit, '(a)')  'File Version Number=1'
   write(ACE_unit, '(a)')  ''
   write(ACE_unit, '(a)')  'Colors Refered to: Mendeleïev Table'
   write(ACE_unit, '(a)')  ''
   write(ACE_unit, '(a)')  '----------------------'
   write(ACE_unit, '(a)')  'Cell Definition'
   write(ACE_unit, '(a)')  '----------------------'
   write(ACE_unit, '(a)')  ''
   write(ACE_unit, '(a)')  'Cell Parameters(Å and °)'
   do i=1, 6
    write(cell_string(i), '(F15.5)') unit_cell%param(i)
   end do
   write(ACE_unit, '(2a,a1,2a,a1,2a)')  'a=',trim(adjustl(cell_string(1))), tab,   &
                                        'b=',trim(adjustl(cell_string(2))), tab,   &
                                        'c=',trim(adjustl(cell_string(3)))

   write(ACE_unit, '(2a,a1,2a,a1,2a)')  'alpha=',trim(adjustl(cell_string(4))), tab,   &
                                        'beta=', trim(adjustl(cell_string(5))), tab,   &
                                        'gamma=',trim(adjustl(cell_string(6)))

   write(ACE_unit, '(a)')          ''
   write(ACE_unit, '(a)')          'System=CEL'
   write(ACE_string, *) SPG%NumSpg
   write(ACE_unit, '(2a)')         'Space Group Number=',trim(adjustl(ACE_string))
   if(SPG%NumSpg == 146 .or. SPG%NumSpg == 148 .or. SPG%NumSpg == 155 .or. SPG%NumSpg == 160 .or. &
      SPG%NumSpg == 161 .or. SPG%NumSpg == 166 .or. SPG%NumSpg == 167) then   ! groupes trigonaux R
      if((abs(unit_cell%param(4) -  90.) < eps) .and. (abs(unit_cell%param(4) - 90.) < eps) .and.   &
         (abs(unit_cell%param(6) - 120.) < eps)) then
         write(ACE_unit, '(a)')          'Hexagonal description'
      else
         write(ACE_unit, '(a)')          'Trigonal description'
      endif

   endif
   write(ACE_unit, '(a)')          'Only Non Equivalent Positions are listed'
   write(ACE_string, *) nb_atom
   write(ACE_unit, '(2a)')         'Number of positions in Cell=',trim(adjustl(ACE_string))
   write(ACE_unit, '(a)')          ''
   write(ACE_unit, '(6(a,a1),a)')  'Atom',tab,'Oxi.',tab,'X',tab,'Y',tab,'Z',tab,'R(Å)',tab,'Occ.'
   do i=1,nb_atom
    write(ACE_string, '(F5.2)') atom_occ_perc(i)
    write(coord_string(1), '(F8.5)')  atom_coord(1:1,i)
    write(coord_string(2), '(F8.5)')  atom_coord(2:2,i)
    write(coord_string(3), '(F8.5)')  atom_coord(3:3,i)
    write(ACE_unit, '(6(a,a1),a)')  trim(adjustl(atom_typ(i))), tab, '0+', tab,          &
                                    trim(adjustl(coord_string(1))), tab,                 &
                                    trim(adjustl(coord_string(2))), tab,                 &
                                    trim(adjustl(coord_string(3))), tab,                 &
                                    '1',tab, trim(adjustl(ACE_string))
   end do

  close(unit = ACE_unit)

  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(ACE_file_name), ' file has been created.'
  call write_info(trim(message_text))
  call write_info('')

  return
end subroutine create_ACE_from_CIF

!---------------------------------------------------------------------------------------------------------
subroutine create_PAT_PRF
 use cryscalc_module,   only : keyword_read_CIF, keyword_read_INS, CIF_file_name, INS_file_name, HKL_2theta, &
                               X_min, X_max, create_PAT, write_HKL, keyword_GENHKL, wavelength, beam_type, pdp_simu
 USE Pattern_profile_module, only : X_pattern, N_pattern
 USE wavelength_module,      ONLY : X_target


 implicit none


  HKL_2THETA     = .true.
  create_PAT     = .true.
  write_HKL      = .true.
  keyword_GENHKL = .true.
  wavelength  = pdp_simu%wave
  beam_type   = pdp_simu%beam

  !if(beam_type(1:8) == 'neutrons') then
  if(pdp_simu%beam(1:8) == 'neutrons') then
   X_min = N_pattern%xmin
   X_max = N_pattern%xmax
  else
   X_min = X_pattern%xmin
   X_max = X_pattern%xmax
  end if

  call generate_HKL()

  return



end subroutine create_PAT_PRF


!----------------------------------------------------------------------------------------------------
subroutine Get_CIF_numors(input_unit)
 use CIF_module

 implicit none
  integer, intent(in) :: input_unit

 ! champ _atom_site_type_symbol
  call Get_CIF_champ_numor(input_unit, '_atom_site_type_symbol', CIF_atom_site_type_symbol_numor)

  ! champ _atom_site_fract_x,y,z
  call Get_CIF_champ_numor(input_unit, '_atom_site_fract_x', CIF_atom_site_fract_x_numor)
  call Get_CIF_champ_numor(input_unit, '_atom_site_fract_y', CIF_atom_site_fract_y_numor)
  call Get_CIF_champ_numor(input_unit, '_atom_site_fract_z', CIF_atom_site_fract_z_numor)

 ! champ _atom_site_U
  call Get_CIF_champ_numor(input_unit, '_atom_site_U_iso_or_equiv', CIF_atom_site_U_numor)

 ! champ _atom_site_adp_type
  call Get_CIF_champ_numor(input_unit, '_atom_site_adp_type', CIF_atom_adp_type_numor)

 ! champ _atom_site_calc_flag
  call Get_CIF_champ_numor(input_unit, '_atom_site_calc_flag', CIF_atom_site_calc_flag_numor)

 ! champ _atom_site_refinement_flags
  call Get_CIF_champ_numor(input_unit, '_atom_site_refinement_flags', CIF_atom_site_refinement_flags_numor)

 ! champ _atom_site_refinement_flags_occupancy
  call Get_CIF_champ_numor(input_unit, '_atom_site_refinement_flags_posn', CIF_atom_site_refinement_flags_posn_numor)

 ! champ _atom_site_refinement_flags_occupancy
  call Get_CIF_champ_numor(input_unit, '_atom_site_refinement_flags_adp', CIF_atom_site_refinement_flags_adp_numor)

 ! champ _atom_site_refinement_flags_occupancy
  call Get_CIF_champ_numor(input_unit, '_atom_site_refinement_flags_occupancy', CIF_atom_site_refinement_flags_occ_numor)

 ! champ _atom_site_occupancy
  call Get_CIF_champ_numor(input_unit, '_atom_site_occupancy', CIF_atom_site_occupancy_numor)

 ! champ _atom_site_symmetry_multiplicity / _atom_site_site_symmetry_order
  call Get_CIF_champ_numor(input_unit, '_atom_site_symmetry_multiplicity', CIF_atom_site_symm_multiplicity_numor)
  if(CIF_atom_site_symm_multiplicity_numor == 0) then
  call Get_CIF_champ_numor(input_unit, '_atom_site_site_symmetry_order', CIF_atom_site_symm_multiplicity_numor)
  end if


 ! champ _atom_site_disorder_assembly
  call Get_CIF_champ_numor(input_unit, '_atom_site_disorder_assembly', CIF_atom_site_disorder_assembly_numor)

 ! champ _atom_site_disorder_group
  call Get_CIF_champ_numor(input_unit, '_atom_site_disorder_group', CIF_atom_site_disorder_group_numor)


 return
end subroutine Get_CIF_numors
!----------------------------------------------------------------------------------------------------------

subroutine Extract_from_CIF
 ! extract .res and .hkl file from a CIF file

  use cryscalc_module, only : archive_cif, CIF_unit, INS_unit, HKL_unit
  use macros_module,   only : test_file_exist
  use IO_module

  implicit none
   logical                     :: file_exist
   character (len=256)         :: read_line
   character (len=256)         :: RES_file, HKL_file, FAB_file
   integer                     :: i1, long, i_error, pv

  i1 = index(archive_CIF, '.')
  if(i1 == 0) then
   archive_cif = trim(archive_cif)//'.cif'
   i1 = index(archive_CIF, '.')
  end if

  call test_file_exist(trim(archive_cif), file_exist, 'out')
  if(.not. file_exist) call end_of_program

  RES_file = archive_CIF(1:i1-1)//'_.res'
  HKL_file = archive_CIF(1:i1-1)//'_.hkl'
  FAB_file = archive_CIF(1:i1-1)//'_.fab'


  open(unit=CIF_unit, file = trim(archive_cif), iostat=i_error)
   do
    read(unit=CIF_unit, fmt='(a)', iostat=i_error) read_line
    if(i_error /=0) exit
    !read_line = adjustl(read_line)
    if(read_line(1:1) == '#') cycle

    long = len_trim(read_line)
    if(long == 0) cycle
    if(long == 15) then
     if(read_line(1:15) == '_shelx_res_file') then
      open(unit= INS_unit, file=trim(RES_file))
      call write_info('    >> Extract .RES file : '//trim(RES_file))
      do
       read(unit=CIF_unit, fmt='(a)', iostat=i_error) read_line
       if(i_error /=0) exit
       !read_line = adjustl(read_line)
       long = len_trim(read_line)
       if(long == 1) then
        if(read_line(1:1) == ';') cycle
       end if
       if(long >= 19) then
        if(read_line(1:19) == '_shelx_res_checksum') exit
       end if
       write(INS_unit, '(a)') trim(read_line)
      end do
      close (unit = INS_unit)
     end if

     if(read_line(1:15) == '_shelx_hkl_file') then
      open(unit= HKL_unit, file=trim(HKL_file))
      call write_info('    >> Extract .HKL file : '//trim(HKL_file))
      do
       read(unit=CIF_unit, fmt='(a)', iostat=i_error) read_line
       if(i_error /=0) exit
       !read_line = adjustl(read_line)
       long = len_trim(read_line)
       if(long == 1) then
        if(read_line(1:1) == ';') cycle
       end if
       if(long >= 19) then
        if(read_line(1:19) == '_shelx_hkl_checksum') exit
       end if
       write(HKL_unit, '(a)') trim(read_line)
      end do
      close (unit = HKL_unit)
     end if

     if(read_line(1:15) == '_shelx_fab_file') then
      pv = 0
      open(unit= HKL_unit, file=trim(FAB_file))
      call write_info('    >> Extract .FAB file : '//trim(FAB_file))
      do
       read(unit=CIF_unit, fmt='(a)', iostat=i_error) read_line
       if(i_error /=0) exit
       !read_line = adjustl(read_line)
       long = len_trim(read_line)
       if(long == 1) then
        if(read_line(1:1) == ';') then
         pv = pv + 1
         if(pv ==1) cycle
         if(pv ==2) exit
        end if
       end if
       !if(long >= 33) then
       ! if(read_line(1:33) == '_platon_squeeze_void_probe_radius') exit
       !end if
       write(HKL_unit, '(a)') trim(read_line)
      end do
      close (unit = HKL_unit)
     end if


    end if

   end do
  close(unit=CIF_unit)

 return
end subroutine Extract_from_CIF

!!------------------------------------------------------------------------------------------------------------!!
 subroutine write_CIF_output(input_string)
  use CRYSCALC_module, only : CIF_read_unit, CIF_file_name
  USE IO_module
  implicit none
  character (len=*), intent(in)     :: input_string
  character (len=256)               :: read_line, output_string
  character (len=16), dimension(16) :: CIF_str
  character (len=64)                :: fmt_, fmt_final
  integer                           :: long, i_error
  logical                           :: ask_bonds, ask_angles, ask_torsion, ask_htab
  logical                           :: ans_bonds, ans_angles, ans_torsion, ans_htab

  ask_bonds    = .false.
  ask_angles   = .false.
  ask_torsion  = .false.
  ask_htab     = .false.
  ans_bonds    = .false.
  ans_angles   = .false.
  ans_torsion  = .false.
  ans_htab     = .false.


  long = len_trim(input_string)
  if(long == 5) then
   if(input_string(1:5) == "BONDS") ask_bonds = .true.
  elseif(long == 6) then
   if(input_string(1:6) == "ANGLES") ask_angles = .true.
  elseif(long == 7) then
   if(input_string(1:7) == "TORSION") ask_torsion = .true.
  elseif(long == 4) then
   if(input_string(1:4) == "HTAB") ask_htab = .true.
  end if

  OPEN(UNIT=CIF_read_unit, FILE=TRIM(CIF_file_name), ACTION="read")

    do
      read(CIF_read_unit, '(a)', iostat=i_error) read_line
      if(i_error /=0) exit
      !if(len_trim(read_line) == 0) exit
      read_line = adjustl(read_line)

      if(ask_BONDS) then
       if(index(read_line, '_geom_bond_publ_flag') /=0) then
        do
         read(CIF_read_unit, '(a)', iostat=i_error) read_line
         if(i_error /=0) exit
          if(len_trim(read_line) == 0) exit
          !read(read_line, *) CIF_str(1:5)
         !write(output_string, '(2a5, a12,2a12)') CIF_str(1:5)
         read(read_line, *) CIF_str(1:3)
         write(output_string, '(2a5, a12)') CIF_str(1:3)
         call write_info(trim(output_string))
        end do
        ans_BONDS = .true.
       end if

      elseif(ask_ANGLES) then
       if(index(read_line, '_geom_angle_publ_flag') /=0) then
        do
         read(CIF_read_unit, '(a)', iostat=i_error) read_line
         if(i_error /=0) exit
          if(len_trim(read_line) == 0) exit
          !read(read_line, *) CIF_str(1:7)
         !call get_fmt(1, CIF_str(4:4), fmt_, 100)
         !write(fmt_final, '(3a)') '(3a5,', trim(fmt_), ',3a8)'
         !write(output_string, fmt_final) CIF_str(1:7)
          read(read_line, *) CIF_str(1:4)
         call get_fmt(1, CIF_str(4:4), fmt_, 100)
         write(fmt_final, '(3a)') '(3a5,', trim(fmt_), ')'
         write(output_string, fmt_final) CIF_str(1:4)
         call write_info(trim(output_string))
        end do
        ans_ANGLES = .true.
       end if

         elseif(ask_TORSION) then
       if(index(read_line, '_geom_torsion_publ_flag') /=0) then
        do
         read(CIF_read_unit, '(a)', iostat=i_error) read_line
         if(i_error /=0) exit
          if(len_trim(read_line) == 0) exit
          !read(read_line, *) CIF_str(1:10)
         !call get_fmt(1, CIF_str(5:5), fmt_, -100)
         !write(fmt_final, '(3a)') '(4a5,', trim(fmt_), ',x,5a8)'
         !write(output_string, fmt_final) CIF_str(1:10)
          read(read_line, *) CIF_str(1:5)
         call get_fmt(1, CIF_str(5:5), fmt_, -100)
         write(fmt_final, '(3a)') '(4a5,', trim(fmt_), ')'
         write(output_string, fmt_final) CIF_str(1:5)
         call write_info(trim(output_string))
        end do
        ans_TORSION = .true.
       end if

        elseif(ask_HTAB) then
       if(index(read_line, '_geom_hbond_site_symmetry_A') /=0) then
        do
         read(CIF_read_unit, '(a)', iostat=i_error) read_line
         if(i_error /=0) exit
          if(len_trim(read_line) == 0) exit
         read(read_line, *) CIF_str(1:8)
         write(output_string, '(3a5, 4a12,a6)') CIF_str(1:8)
         call write_info(trim(output_string))
        end do
        ans_HTAB = .true.
       end if
      end if
     end do


  CLOSE(UNIT=CIF_read_unit)

  if(ask_BONDS .and. .not. ans_BONDS) then
   call write_info('')
   call write_info('  No bonds founded in the CIF file !!')
   call write_info('')
   return
  end if

  if(ask_BONDS .and. .not. ans_BONDS) then
   call write_info('')
   call write_info('  No bonds founded in the CIF file !!')
   call write_info('')
   return
  end if

  if(ask_TORSION .and. .not. ans_TORSION) then
   call write_info('')
   call write_info('  No torsion angles founded in the CIF file !!')
   call write_info('')
   return
  end if

  if(ask_HTAB .and. .not. ans_HTAB) then
   call write_info('')
   call write_info('  No hydrogen bonds founded in the CIF file !!')
   call write_info('')
   return
  end if





   return
 end subroutine write_CIF_output
!!------------------------------------------------------------------------------------------------------------!!
