subroutine read_cif_and_create_ARCHIVE
! read CIF files and create ARCHIVE.cif
 use cryscalc_module, only : nb_input_CIF_files, input_CIF_file, CIF_file_exist, CIF_unit, CRYSCALC, CIF_read_unit, &
                             molecule, absorption, SADABS, SPG, device, experiment, Proposer, SFAC, nb_atoms_type,  &
                             include_RES_file, include_HKL_file, SQUEEZE, CIF_format80, on_screen_CIF_item, CIFDEP, &
                             Author, include_experimenter_name, input_line, write_details, message_text
 use CIF_module
 use macros_module,   only : Get_current_folder_name, Get_sample_ID_from_folder_name, u_case, l_case, test_file_exist, &
                             Nombre_de_colonnes, Replace_car, Get_Proposer_from_current_folder_name
 use IO_module,       only : read_input_line, write_info

 implicit none
  integer                           :: i, i1, i2, i3, long, i_error
  character (len=64)                :: fmt_, fmt_final
  character (len=256)               :: current_folder
  character (len=256)               :: sample_ID, diffraction_date, experimenter
  character (len=256)               :: CIF_sep_line
  character (len=256)               :: read_line
  character (len=256)               :: CIF_string
  character (len=16), dimension(16) :: CIF_str
  CHARACTER (LEN=256)               :: string_value
  logical                           :: on_one_line
  logical                           :: file_exist, import_cif_exist, ok, part_occ
  integer                           :: nb_files_exist, nb_col
  character (len=8)                 :: ans



  ! test si les fichiers .CIF existent
  nb_files_exist = 0
  import_cif_exist = .false.
  do i=1, nb_input_CIF_files
   call test_file_exist(input_CIF_file(i), file_exist, 'out')
   if(file_exist) then
    CIF_file_exist(i) = .true.
    nb_files_exist =  nb_files_exist + 1
    if(i==2) import_cif_exist = .true.
   end if
  end do
  if(nb_files_exist == 0) then
   call write_info('')
   call write_info('CIF files do not exist. Program will be stopped !')
   call write_info('')
   stop
  elseif(nb_files_exist /= nb_input_CIF_files) then
   call write_info('')
   call write_info(' >>> Final archive CIF file will not be completed !!')
   call write_info('     Do you want to continue (y/n[=def]) ?')
   do
    call read_input_line(input_line)
    read(input_line, *) ans
    !read(*,*) ans
    if(ans(1:1) == 'Y' .or. ans(1:1) == 'y') then
     exit
    else
     stop
    end if
   end do


   call write_info('')

  end if


  write(CIF_sep_line, '(a,76a1,a)') '#', ('-',i=1,76), '#'

  close(unit=CIF_unit)
  if(include_HKL_file) then
   open(unit=CIF_unit, file = "cryscalc_archive_hkl.cif")
  else
   open(unit=CIF_unit, file = "cryscalc_archive.cif")
  end if

  write(CIF_unit, '(a)') '#'
  write(CIF_unit, '(a)') '#'

  write(CIF_unit, '(a)') trim(CIF_sep_line)
  write(CIF_string, '(3a)') "CIF file created and formatted by CRYSCALC (",  trim(CRYSCALC%version), ")"
  call center_CIF_line(trim(CIF_string))
  write(CIF_unit, '(a)') trim(CIF_sep_line)
  WRITE(CIF_unit, '(a)') trim(CIF_sep_line_empty)
  call write_cif_line("Input CIF files:")

  !long = len_trim(CRYSCALC%version)
  !write(fmt_ , '(a,i2,a)') '(3a,', 30-long, 'x,a)'
  !write(CIF_unit, fmt=trim(fmt_))        '# CIF file created and formatted by CRYSCALC (', trim(CRYSCALC%version),')','#'
  !write(CIF_unit, '(a,59x,a)')           '# Input CIF files:','#'

  call write_info('   Input CIF files used to create final archive:')
  !do i=1, nb_input_CIF_files
  do i=1, nb_files_exist
   !long=len_trim(input_cif_file(i))
   !write(fmt_ , '(a,i2,a)') '(a,i2,2a,', 71-long, 'x,a)'
   !write(CIF_unit, fmt=trim(fmt_))        '# ', i, '. ',trim(Input_CIF_file(i)),'#'

   call test_file_exist(input_CIF_file(i), file_exist, 'no_out')
   !call test_file_exist(input_CIF_file(i), file_exist, 'out')
   if(.not. file_exist) cycle

   write(CIF_string, '(i2,2a)') i, '. ',trim(Input_CIF_file(i))
   call write_cif_line(trim(CIF_string))
   write(message_text, fmt='(5x,i2,2a)')    i, '. ' , trim(Input_CIF_file(i))
   call write_info(trim(message_text))
  end do
  call write_info('')
  WRITE(CIF_unit, '(a)') trim(CIF_sep_line_empty)
  write(CIF_unit, '(a)') trim(CIF_sep_line)


  if(on_screen_CIF_item) call write_info('')

  call Get_current_folder_name(current_folder)
  !diffraction_date = ''
  !i1 = index(current_folder, '_', back=.true.)
  !if (i1 /=0) then
  ! diffraction_date = current_folder(i1+1:)
  ! i1 = index(diffraction_date, '\')
  ! if (i1 /=0) diffraction_date = diffraction_date(1:i1-1)
  ! if(len_trim(diffraction_date) /= 0) then
  !  call Date_format(diffraction_date)
  !  write(CIF_string, '(2a)') 'Date of experiment: ', trim(diffraction_date)
  !  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
  !  call write_cif_line(trim(CIF_string))
  ! end if
  !end if

  sample_ID        = ''
  diffraction_date = ''
  experimenter     = ''

  if(import_cif_exist) then
   call Get_items_from_import_cif_header(diffraction_date, experimenter, Proposer%name, sample_ID)
   if(len_trim(diffraction_date) /=0 .and. diffraction_date(1:1) /= '?') then
    call Date_format_2(diffraction_date)
   else
    diffraction_date(1:1) = '?'
   end if
   write(CIF_string, '(2a)') 'Date of experiment: ', trim(diffraction_date)
   WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
   call write_cif_line(trim(CIF_string))
  end if



  if(include_experimenter_name) then
   if(len_trim(experimenter) /= 0) then
    write(CIF_string, fmt = '(2a)') 'Experimenter      : ', trim(experimenter)
    call write_cif_line(trim(CIF_string))
	WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
   else
    if(AUTHOR%first_name(1:1) /= "?" .and.  AUTHOR%name(1:1) /= "?" .and. AUTHOR%TEAM(1:1) /= '?') then
     !WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
     write(CIF_string, fmt = '(4a)') 'Experimenter      : ', &
                                      trim(AUTHOR%first_name)//' '//trim(AUTHOR%name)//" ["//trim(AUTHOR%team)//"]"
     call write_cif_line(trim(CIF_string))
	 WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
    end if
   end if

   long = len_trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
   write(CIF_string, fmt = '(2a)')  'Instrument        : ', trim(CIF_parameter_DEVICE%diffrn_measurement_device_type(2:long-1))
   call write_cif_line(trim(CIF_string))
   write(CIF_string, fmt = '(2a)')  'Laboratory        : ', trim(DEVICE%lab)
   call write_cif_line(trim(CIF_string))
   !WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
  end if

  !call Get_experiment_features
  !if(len_trim(experiment%date) /=0) then
  ! write(CIF_string, '(2a)')   "Date of experiment: ", trim(Experiment%date)
  ! call write_cif_line(trim(CIF_string))
  !end if

  if(len_trim(sample_ID) == 0) then
   ! call Get_sample_ID_from_folder_name(current_folder, sample_ID) ! comment\'e en janvier 2016
   call Get_sample_ID_from_CIF(input_CIF_file(1), sample_ID)
   !call Get_sample_ID_from_import_cif(sample_ID)   ! oct. 2016
  end if

  write(CIF_string, '(2a)') "Sample ID         : ", trim(sample_ID)
  call write_cif_line(trim(CIF_string))

  !write(CIF_string, '(2a)') "Proposer          : ", trim(sample_ID)

  !call Get_current_folder_name(current_folder)
  !call Get_Proposer_from_current_folder_name(current_folder, proposer%name, proposer%team)
  !Proposer%name = replace_car(Proposer%name, '_', '. ')

  !long = len_trim(Proposer%team)
  !team_date = .false.
  !if(long == 4) then
  ! if(Proposer%team(1:4) == '2015' .or. Proposer%team(1:4) == '2016') then
  !  team_date = .true.
  ! end if
  !end if

  !if(team_date) call Get_team_ISCR(trim(Proposer%name), Proposer%team)
  !call test_ISCR(trim(proposer%team), Proposer%ISCR)

  !if(proposer%name(1:1) /= '?' .and. proposer%team(1:1) /= '?') then
  ! if(Proposer%ISCR) then
  !  write(CIF_string, '(5a)')  "Proposer          : ", trim(proposer%name), ' [',trim(proposer%team),'/ISCR]'
  ! else
  !  write(CIF_string, '(5a)')  "Proposer          : ", trim(proposer%name), ' [',trim(proposer%team),']'
  ! end if
  ! call write_cif_line(trim(CIF_string))
  !end if

  if(proposer%name(1:1) /= '?' ) then
   write(CIF_string, '(2a)')  "Proposer          : ", trim(proposer%name)
   call write_cif_line(trim(CIF_string))
  end if



  WRITE(CIF_unit, '(a)') trim(CIF_sep_line_empty)
  write(CIF_unit, '(a)') trim(CIF_sep_line)
  write(CIF_unit, '(a)') '#'
  write(CIF_unit, '(a)') ''

  if(CIFDEP) then
   call write_CIF_AUTHOR(CIF_unit)
   call write_CIF_text(CIF_unit)
  endif


  !call Get_sample_ID_from_folder_name(current_folder, sample_ID)

  !if(sample_ID == '?') call Get_Wingx_job(sample_ID, wingx_structure_dir)
  !if(sample_ID == '?') sample_ID = trim(read_line(6:))
  write(CIF_unit, '(2a)') 'data_', u_case(trim(sample_ID))
  write(CIF_unit, '(a)') ''
   call write_CIF_file('CHEMICAL_INFORMATION')
  write(CIF_unit, '(a)') ''
  write(CIF_unit, '(a)') '_chemical_name_systematic'
  write(CIF_unit, '(a)') ';'
  write(CIF_unit, '(a)') ' ?'
  write(CIF_unit, '(a)') ';'

  do i=1, nb_input_CIF_files
  ! lecture fichier.CIF principal  (#1)
   if(CIF_file_exist(i)) then
    close(unit=CIF_read_unit)
    open(unit = CIF_read_unit, file = trim(input_CIF_file(i)))
    call read_CIF_input_file(trim(input_CIF_file(i)), '?')
    call read_CIF_input_file_TR(CIF_read_unit, trim(input_CIF_file(i)))
   endif
  end do

  !! lecture fichier.CIF
  !close (unit=CIF_read_unit)
  !open(unit = CIF_read_unit, file = trim(input_CIF_file(2)))
  !call read_CIF_input_file(trim(input_CIF_file(2)), '?')
  !call read_CIF_input_file_TR(CIF_read_unit)

  ! >> creation de l'objet MOLECULE
  ! call atomic_identification()
  ! call get_content
  ! call molecular_weight
  ! call density_calculation

  ! le bon Z est dans le fichier.CIF principal
  close(unit=CIF_read_unit)

  if(.not. CIF_file_exist(1)) then
   call write_info('')
   call write_info(' '//trim(input_CIF_file(1))//' file does not exist. Program will be stopped!')
   call write_info('')
   stop
  end if

  open(unit=CIF_read_unit, file = trim(input_CIF_file(1)))

  CIF_string = '_cell_formula_units_Z'
  ok = .false.
  call get_champ_value(CIF_read_unit, trim(CIF_string), string_value, ok)
  IF(ok) READ(string_value, '(a)') CIF_parameter%cell_formula_units_Z

  ! new mai 2015 ---------------------------------------------------------
  ! la bonne formule est dans le fichier.CIF principal
  CIF_string = '_chemical_formula_sum'
  ok = .false.
  call get_champ_value(CIF_read_unit, trim(CIF_string), string_value, ok)
  long = len_trim(string_value)
  if(ok) read(string_value(2:long-1), '(a)') molecule%formula
  call check_formula(molecule%formula)   ! verif. si stoechiometrie nulle
  ! -----------------------------------------------------------------------

  ! determination du nombre de champ dans la boucle _atom_site_
  call Get_CIF_champ_nb(CIF_read_unit, '_atom_site_label', CIF_atom_site_item_nb, CIF_atom_site_loop_numor)


  ! determination des positions des differents items dans la boucle
  call Get_CIF_numors(CIF_read_unit)

  molecule%formula = u_case(molecule%formula)
  call Nombre_de_colonnes(molecule%formula, nb_col)
  read(molecule%formula, *)  CIF_str(1:nb_col)

  CIF_string = ''
  do i=1, nb_col
   long = len_trim(CIF_str(i))
   if(long == 1) then
    write(CIF_string, '(3a)') trim(CIF_string), ' ', u_case((CIF_str(i)(1:1)))
   else
    write(CIF_string, '(4a)') trim(CIF_string), ' ', u_case((CIF_str(i)(1:1))), l_case((CIF_str(i)(2:)))
   end if
  end do
  write(CIF_unit, '(3a)')     "_chemical_formula_moiety              '", trim(adjustl(CIF_string)),"'"
  write(CIF_unit, '(3a)')     "_chemical_formula_sum                 '", trim(adjustl(CIF_string)),"'"

  !write(CIF_unit, '(3a)')     "_chemical_formula_moiety         '", trim(adjustl(molecule%formula)),"'"
  !write(CIF_unit, '(3a)')     "_chemical_formula_sum            '", trim(adjustl(molecule%formula)),"'"
   !write(CIF_unit, '(a,F6.2)') "_chemical_formula_weight          ", molecule%weight
  write(CIF_unit, '(2a)')     "_chemical_formula_weight               ", trim(CIF_parameter%formula_weight)
  write(CIF_unit, '(a)')      "_chemical_compound_source             'synthesis as described'"
  if(CIF_parameter%abs_structure_flack /= '?') then
  write(CIF_unit, '(a)')      "_chemical_absolute_configuration      'rm/ad/rmad/syn/unk/.'"
  end if

  call write_CIF_file('UNIT_CELL_INFO')

  call GET_SG_from_CIF(trim(input_CIF_file(1)))    ! le bon groupe est dans le fichier .CIF principal

  ! ------------- july 2015   --------------------------------------------
  ! recup. des op. de sym. (l'ordre peut etre different de celui de CFML !
  !call get_SG_op(trim(input_CIF_file(1)))
  CIF_string = '_space_group_symop_operation_xyz'
  ok = .false.

  !inquire(unit=CIF_read_unit, opened=file_opened)
  !if(.not. file_opened) open(unit=CIF_read_unit, file = trim(input_CIF_file(1)))

  call get_champ_value(CIF_read_unit, trim(CIF_string), string_value, ok)
  if(ok) then
   backspace(unit=CIF_read_unit)
   do i=1, SPG%Multip
    read(CIF_read_unit, '(a)', iostat=i_error) read_line
    if(i_error /=0) exit
    if(len_trim(read_line) == 0) exit
    CIF_parameter%sym_op(i) = read_line
   end do
  end if

  ! ----------------------------------------------------------------------


  call write_CIF_file('SPACE_GROUP')
  call write_CIF_file('CELL_PARAM_ESD')
  call write_CIF_file('CIF_THMIN,THMAX')

  !call write_CIF_file('CELL_WAVE')
  write(CIF_unit, '(a,F10.5)')      "_cell_measurement_wavelength        ", CIF_cell_measurement%wavelength

  call write_CIF_file('CRYSTAL_INFORMATION')
  call write_CIF_file('CRYSTAL_DENSITY')
  call write_CIF_file('ABSORPTION')
  read(SADABS%absorption_coeff, *) ABSORPTION%mu
  ABSORPTION%mu = 10. * ABSORPTION%mu
  call write_mu('no')
  call write_CIF_file('SADABS')
  call write_CIF_file('TMIN')
  call write_CIF_file('TMAX')


  call read_SQUEEZE_file(trim(sample_ID))
  !if(SQUEEZE%procedure) then
  ! write(CIF_unit, '(a)') ""
  ! WRITE(CIF_unit, '(a)') "#----------------------------------------------------------------------------#"
  ! WRITE(CIF_unit, '(a)') "#                   SQUEEZE RESULTS                                          #"
  ! WRITE(CIF_unit, '(a)') "#----------------------------------------------------------------------------#"
  ! write(CIF_unit, '(a)') ""
  ! do i=2, SQUEEZE%nb_lines
  !  write(CIF_unit, '(a)') trim(SQUEEZE%read_line(i))
  ! end do
  !endif

  call write_CIF_file('DATACOL')
  WRITE(CIF_unit ,'(2a)') "_diffrn_ambient_temperature             ", TRIM(CIF_cell_measurement%temperature)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_wavelength            ", trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_type                  ", trim(CIF_parameter_DEVICE%diffrn_radiation_type)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_probe                 ", trim(CIF_parameter_DEVICE%diffrn_radiation_probe)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_source                ", trim(CIF_parameter_DEVICE%diffrn_radiation_source)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_monochromator         ", trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_av_r_equivalents         ", trim(CIF_parameter%diffrn_reflns_av_R_equivalents)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_av_uneti/neti            ", trim(CIF_parameter%diffrn_reflns_av_R_sigma)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_number                   ", trim(CIF_parameter%diffrn_reflns_number)

  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_limit_h_min              ", trim(CIF_parameter%h_min)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_limit_h_max              ", trim(CIF_parameter%h_max)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_limit_k_min              ", trim(CIF_parameter%k_min)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_limit_k_max              ", trim(CIF_parameter%k_max)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_limit_l_min              ", trim(CIF_parameter%l_min)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_limit_l_max              ", trim(CIF_parameter%l_max)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_theta_min                ", trim(CIF_parameter%theta_min)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_theta_max                ", trim(CIF_parameter%theta_max)
  WRITE(CIF_unit, '(2a)') "_diffrn_reflns_theta_full               ", trim(CIF_parameter%theta_full)
  WRITE(CIF_unit, '(2a)') "_diffrn_measured_fraction_theta_full    ", trim(CIF_parameter_DEVICE%diffrn_theta_full)
  WRITE(CIF_unit, '(2a)') "_diffrn_measured_fraction_theta_max     ", trim(CIF_parameter_DEVICE%diffrn_theta_max)
  WRITE(CIF_unit, '(2a)') "_reflns_number_total                    ", trim(CIF_parameter%reflns_number_total)
  WRITE(CIF_unit, '(2a)') "_reflns_number_gt                       ", trim(CIF_parameter%reflns_number_gt)
  WRITE(CIF_unit, '(2a)') "_reflns_threshold_expression            ", ">2sigma(I)"

  if(CIF_parameter%Friedel_coverage(1:1) /="?") then
   WRITE(CIF_unit, '(2a)') "_reflns_Friedel_coverage                ", trim(CIF_parameter%Friedel_coverage)
   WRITE(CIF_unit, '(2a)') "_reflns_Friedel_fraction_max            ", trim(CIF_parameter%Friedel_fraction_max)
   WRITE(CIF_unit, '(2a)') "_reflns_Friedel_fraction_full           ", trim(CIF_parameter%Friedel_fraction_full)
   WRITE(CIF_unit, '(a)')  ""
   WRITE(CIF_unit, '(a)')  "_reflns_special_details"
   WRITE(CIF_unit, '(a)')  ";"
   WRITE(CIF_unit, '(a)')  "Reflections were merged by SHELXL according to the crystal"
   WRITE(CIF_unit, '(a)')  "class for the calculation of statistics and refinement."
   WRITE(CIF_unit, '(a)')  ""
   WRITE(CIF_unit, '(a)')  "_reflns_Friedel_fraction is defined as the number of unique"
   WRITE(CIF_unit, '(a)')  "Friedel pairs measured divided by the number that would be"
   WRITE(CIF_unit, '(a)')  "possible theoretically, ignoring centric projections and"
   WRITE(CIF_unit, '(a)')  "systematic absences."
   WRITE(CIF_unit, '(a)')  ";"
  endif



  call write_CIF_file('COMPUTER_PROGRAMS')

  call write_CIF_file('REFINEMENT_INFO')
  WRITE(CIF_unit, '(a)')  "_refine_ls_weighting_details            "
  WRITE(CIF_unit, '(2a)')  "        ", trim(CIF_parameter%refine_ls_weighting_details)

  WRITE(CIF_unit, '(2a)')  "_atom_sites_solution_primary            ", trim(CIF_parameter%atom_sites_solution_1)
  WRITE(CIF_unit, '(2a)')  "_atom_sites_solution_secondary          ", trim(CIF_parameter%atom_sites_solution_2)
  WRITE(CIF_unit, '(2a)')  "_atom_sites_solution_hydrogens          ", trim(CIF_parameter%atom_sites_solution_H)


  WRITE(CIF_unit, '(2a)')  "_refine_ls_hydrogen_treatment           ", trim(CIF_parameter%refine_ls_H_treatment)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_extinction_method            ", trim(CIF_parameter%refine_ls_extinction_method)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_extinction_coef              ", trim(CIF_parameter%refine_ls_extinction_coef)

  WRITE(CIF_unit, '(2a)')  "_refine_ls_number_reflns                ", trim(CIF_parameter%refine_ls_number_reflns)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_number_parameters            ", trim(CIF_parameter%refine_ls_number_parameters)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_number_restraints            ", trim(CIF_parameter%refine_ls_number_restraints)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_R_factor_all                 ", trim(CIF_parameter%refine_ls_R_factor_all)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_R_factor_gt                  ", trim(CIF_parameter%refine_ls_R_factor_gt)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_wR_factor_ref                ", trim(CIF_parameter%refine_ls_wR_factor_ref)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_wR_factor_gt                 ", trim(CIF_parameter%refine_ls_wR_factor_gt)

  WRITE(CIF_unit, '(2a)')  "_refine_ls_goodness_of_fit_ref          ", trim(CIF_parameter%refine_ls_goodness_of_fit_ref)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_restrained_S_all             ", trim(CIF_parameter%refine_ls_restrained_S_all)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_shift/su_max                 ", trim(CIF_parameter%refine_ls_shift_su_max)
  WRITE(CIF_unit, '(2a)')  "_refine_ls_shift/su_mean                ", trim(CIF_parameter%refine_ls_shift_su_mean)
  WRITE(CIF_unit, '(2a)')  "_refine_diff_density_max                ", trim(CIF_parameter%refine_diff_density_max)
  WRITE(CIF_unit, '(2a)')  "_refine_diff_density_min                ", trim(CIF_parameter%refine_diff_density_min)
  WRITE(CIF_unit, '(2a)')  "_refine_diff_density_rms                ", trim(CIF_parameter%refine_diff_density_rms)
  if(CIF_parameter%abs_structure_flack /= '?') then
  !WRITE(CIF_unit, '(2a)')  "_refine_ls_abs_structure_details        ", trim(CIF_parameter%abs_structure_details)
  WRITE(CIF_unit, '(a)')   "_refine_ls_abs_structure_details        "
  do i=1, n_details
   WRITE(CIF_unit, '(a)')   trim(CIF_parameter%abs_structure_details(i))
  end do
  WRITE(CIF_unit, '(2a)')  "_refine_ls_abs_structure_Flack          ", trim(CIF_parameter%abs_structure_Flack)
  end if


  if(include_RES_file) then
   WRITE(CIF_unit, '(a)')  ""
   call include_RES_file_into_CIF
   write(CIF_unit, '(a)') trim(CIF_sep_line)
  end if

  if(ON_screen_CIF_item) then
   call write_info('  > _chemical_formula_sum                           extracted from '//trim(input_CIF_file(1)) // &
                   ' [' // trim(molecule%formula) // ']')
   call write_info('  > _space group                           features extracted from '//trim(input_CIF_file(1)) // &
                   ' [' // trim(SPG%SPG_Symb) // ']')
  end if

  call write_CIF_file('ATOMS_HEADER')


  ! lecture fichier .CIF
  close(unit=CIF_read_unit)
  open(unit = CIF_read_unit, file = trim(input_CIF_file(1)))
  do
   read(unit = CIF_read_unit, fmt='(a)', iostat=i_error) read_line
   if(i_error /=0) exit

   if(index(read_line, '_atom_type_scat_source') /=0)  then
    if(ON_screen_CIF_item) then
     call write_info('  > _atom_type_                            features extracted from '// trim(input_CIF_file(1)))
    end if

    write(CIF_unit, '(a)') 'loop_'
    write(CIF_unit, '(a)') '    _atom_type_symbol'
    write(CIF_unit, '(a)') '    _atom_type_description'
    write(CIF_unit, '(a)') '    _atom_type_scat_dispersion_real'
    write(CIF_unit, '(a)') '    _atom_type_scat_dispersion_imag'
    write(CIF_unit, '(a)') '    _atom_type_scat_source'

    do
     read(CIF_read_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     on_one_line = .false.
     i3 = index(read_line, "'International Tables")  ! test si le champ est sur 1 ou 2 lignes
     if(i3 /=0) on_one_line = .true.
     read(read_line, *) CIF_str(1:4)
     do i=1, 2
      i1 = index(CIF_str(i), "'")
      i2 = index(CIF_str(i), "'", back=.true.)
      if(i1 /=0 .and. i2 /=0 .and. i2>i1) CIF_str(i) = CIF_str(i)(i1+1:i2-1)
     end do
     if(.not. on_one_line) then
      read(CIF_read_unit, '(a)', iostat=i_error) read_line
      if(i_error /=0) exit
      read_line = adjustl(read_line)
     end if

     if(CIF_str(3)(1:1) /= '-') then
      if(.not. on_one_line) then
       write(CIF_unit, '(2a4, 1x,2a8, a)') CIF_str(1:4), trim(read_line)
      else
       write(CIF_unit, '(2a4, 1x,2a8, a)') CIF_str(1:4), trim(read_line(i3:))
      end if
     else
      if(.not. on_one_line) then
       write(CIF_unit, '(2a4, a8,1x,a8, a)') CIF_str(1:4), trim(read_line)
      else
       write(CIF_unit, '(2a4, a8,1x,a8, a)') CIF_str(1:4), trim(read_line(i3:))
      end if

     endif
    end do
    write(CIF_unit, '(a)') ''

   elseif(index(read_line, '_atom_site_disorder_group') /=0)  then
    if(ON_screen_CIF_item) then
     call write_info('  > _atom_site_                            features extracted from '// trim(input_CIF_file(1)))
    end if
    write(CIF_unit, '(a)') 'loop_'
    if(CIF_atom_site_item_nb == 13) then
     write(CIF_unit, '(a)') '    _atom_site_label'
     write(CIF_unit, '(a)') '    _atom_site_type_symbol'
     write(CIF_unit, '(a)') '    _atom_site_fract_x'
     write(CIF_unit, '(a)') '    _atom_site_fract_y'
     write(CIF_unit, '(a)') '    _atom_site_fract_z'
     write(CIF_unit, '(a)') '    _atom_site_U_iso_or_equiv'
     write(CIF_unit, '(a)') '    _atom_site_adp_type'
     write(CIF_unit, '(a)') '    _atom_site_occupancy'
     write(CIF_unit, '(a)') '    _atom_site_symmetry_multiplicity'
     write(CIF_unit, '(a)') '    _atom_site_calc_flag'
     write(CIF_unit, '(a)') '    _atom_site_refinement_flags'
     write(CIF_unit, '(a)') '    _atom_site_disorder_assembly'
     write(CIF_unit, '(a)') '    _atom_site_disorder_group'
    elseif(CIF_atom_site_item_nb == 15) then
     write(CIF_unit, '(a)') '    _atom_site_label'
     write(CIF_unit, '(a)') '    _atom_site_type_symbol'
     write(CIF_unit, '(a)') '    _atom_site_fract_x'
     write(CIF_unit, '(a)') '    _atom_site_fract_y'
     write(CIF_unit, '(a)') '    _atom_site_fract_z'
     write(CIF_unit, '(a)') '    _atom_site_U_iso_or_equiv'
     write(CIF_unit, '(a)') '    _atom_site_adp_type'
     write(CIF_unit, '(a)') '    _atom_site_occupancy'
     write(CIF_unit, '(a)') '    _atom_site_site_symmetry_order'
     write(CIF_unit, '(a)') '    _atom_site_calc_flag'
     write(CIF_unit, '(a)') '    _atom_site_refinement_flags_posn'
     write(CIF_unit, '(a)') '    _atom_site_refinement_flags_adp'
     write(CIF_unit, '(a)') '    _atom_site_refinement_flags_occupancy'
     write(CIF_unit, '(a)') '    _atom_site_disorder_assembly'
     write(CIF_unit, '(a)') '    _atom_site_disorder_group'
    end if

    do
     read(CIF_read_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     !read(read_line, *) CIF_str(1:13)
     read(read_line, *) CIF_str(1:CIF_atom_site_item_nb)
     if(CIF_format80) then
      call get_fmt(4, CIF_str(3:6), fmt_, 0)

     ! dec. 2015 : format non optimise pour
      !     _atom_site_refinement_flags_posn
      !     _atom_site_refinement_flags_adp

      part_occ=.false.
      if(CIF_atom_site_refinement_flags_numor /=0) then
       if(index(CIF_str(CIF_atom_site_refinement_flags_numor), 'P') /=0) part_occ = .true.
      end if
      if(.not. part_occ) then
       if(CIF_atom_site_refinement_flags_occ_numor /=0) then
        if(index(CIF_str(CIF_atom_site_refinement_flags_occ_numor), 'P') /=0 .and. &
           len_trim(CIF_str(CIF_atom_site_occupancy_numor)) /=1) then
            part_occ = .true.
        end if
       end if
      end if

      !if(index(CIF_str(CIF_atom_site_refinement_flags_numor),     'P') /=0 .or. &
      !   index(CIF_str(CIF_atom_site_refinement_flags_occ_numor), 'P') /=0 .and. &
      !  len_trim(CIF_str(CIF_atom_site_occupancy_numor)) /=1 ) then ! occ. partielle

      if(part_occ) then
       if(len_trim(CIF_str(CIF_atom_site_symm_multiplicity_numor)) == 1) then
        write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a10,1x,a2,a5,a3,a2,a4,2a2)'
       else
        write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a12,1x,a3,a5,a3,a2,a4,2a2)'
       endif

      else
       if(CIF_atom_site_symm_multiplicity_numor/=0) then
        if(len_trim(CIF_str(CIF_atom_site_symm_multiplicity_numor)) == 1) then
         !write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2, 1x,a2,a5,a3,a2,a3,2a2)'
         write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2, a2,a5,a2,a2,a2,2a2)'
        else
         !write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a12,1x,a3,a5,a3,a2,a3,2a2)'
         write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a12,1x,a3,a5,a2,a2,a2,2a2)'
        endif
       else
        !write(fmt_final, '(3a)')  '(a5,a3,', trim(fmt_), ',a5,a2, 1x,a3,a5,a3,a2,a3,2a2)'
        write(fmt_final, '(3a)')  '(a5,a3,', trim(fmt_), ',a5,a2, a3,a5,a2,a2,a2,2a2)'
       endif
      end if
      !write(CIF_unit, fmt_final) CIF_str(1:13)

      write(CIF_unit, fmt_final) CIF_str(1:CIF_atom_site_item_nb)
     else
      write(CIF_unit, '(a)') trim(read_line)
     end if !  fin de la condition 'if (CIF_format80
    end do
    write(CIF_unit, '(a)') ''

   elseif(index(read_line, '_atom_site_aniso_U_12') /=0)  then
    if(ON_screen_CIF_item) then
     call write_info('  > _atom_type_aniso_                      features extracted from '// trim(input_CIF_file(1)))
    end if

    write(CIF_unit, '(a)') 'loop_'
    write(CIF_unit, '(a)') '    _atom_site_aniso_label'
    write(CIF_unit, '(a)') '    _atom_site_aniso_U_11'
    write(CIF_unit, '(a)') '    _atom_site_aniso_U_22'
    write(CIF_unit, '(a)') '    _atom_site_aniso_U_33'
    write(CIF_unit, '(a)') '    _atom_site_aniso_U_23'
    write(CIF_unit, '(a)') '    _atom_site_aniso_U_13'
    write(CIF_unit, '(a)') '    _atom_site_aniso_U_12'
    do
     read(CIF_read_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:7)
     if(.not. CIF_format80) then
      write(CIF_unit, '(a)') trim(read_line)
     else
      call get_fmt(6, CIF_str(2:7), fmt_, 0)
      write(fmt_final, '(3a)') '(a4,', trim(fmt_),')'
      write(CIF_unit, fmt_final) CIF_str(1:7)
     end if
    end do
    write(CIF_unit, '(a)') ''

   elseif(index(read_line, '_geom_bond_publ_flag') /=0) then
    if(ON_screen_CIF_item) then
     call write_info('  > _geom_bond_                            features extracted from '//trim(input_CIF_file(1)))
    end if

    call write_CIF_FILE('MOLECULAR_GEOMETRY')
    write(CIF_unit, '(a)') 'loop_'
    write(CIF_unit, '(a)') '    _geom_bond_atom_site_label_1'
    write(CIF_unit, '(a)') '    _geom_bond_atom_site_label_2'
    write(CIF_unit, '(a)') '    _geom_bond_distance'
    write(CIF_unit, '(a)') '    _geom_bond_site_symmetry_2'
    write(CIF_unit, '(a)') '    _geom_bond_publ_flag'
    do
     read(CIF_read_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:5)
     write(CIF_unit, '(2a5, a12,2a12)') CIF_str(1:5)
    end do
    write(CIF_unit, '(a)') ''

   elseif(index(read_line, '_geom_angle_publ_flag') /=0) then
    if(ON_screen_CIF_item) then
     call write_info('  > _geom_angle_                           features extracted from ' //trim(input_CIF_file(1)))
    end if

    write(CIF_unit, '(a)') 'loop_'
    write(CIF_unit, '(a)') '    _geom_angle_atom_site_label_1'
    write(CIF_unit, '(a)') '    _geom_angle_atom_site_label_2'
    write(CIF_unit, '(a)') '    _geom_angle_atom_site_label_3'
    write(CIF_unit, '(a)') '    _geom_angle'
    write(CIF_unit, '(a)') '    _geom_angle_site_symmetry_1'
    write(CIF_unit, '(a)') '    _geom_angle_site_symmetry_3'
    write(CIF_unit, '(a)') '    _geom_angle_publ_flag'
    do
     read(CIF_read_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:7)
     call get_fmt(1, CIF_str(4:4), fmt_, 100)
     write(fmt_final, '(3a)') '(3a5,', trim(fmt_), ',3a8)'
     write(CIF_unit, fmt_final) CIF_str(1:7)
    end do
    write(CIF_unit, '(a)') ''

   elseif(index(read_line, '_geom_torsion_publ_flag') /=0) then
    if(ON_screen_CIF_item) then
     call write_info('  > _geom_torsion_                         features extracted from '//trim(input_CIF_file(1)))
    end if

    write(CIF_unit, '(a)') 'loop_'
    write(CIF_unit, '(a)') '    _geom_torsion_atom_site_label_1'
    write(CIF_unit, '(a)') '    _geom_torsion_atom_site_label_2'
    write(CIF_unit, '(a)') '    _geom_torsion_atom_site_label_3'
    write(CIF_unit, '(a)') '    _geom_torsion_atom_site_label_4'
    write(CIF_unit, '(a)') '    _geom_torsion'
    write(CIF_unit, '(a)') '    _geom_torsion_site_symmetry_1'
    write(CIF_unit, '(a)') '    _geom_torsion_site_symmetry_2'
    write(CIF_unit, '(a)') '    _geom_torsion_site_symmetry_3'
    write(CIF_unit, '(a)') '    _geom_torsion_site_symmetry_4'
    write(CIF_unit, '(a)') '    _geom_torsion_publ_flag'
    do
     read(CIF_read_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:10)
     call get_fmt(1, CIF_str(5:5), fmt_, -100)
     write(fmt_final, '(3a)') '(4a5,', trim(fmt_), ',x,5a8)'

     write(CIF_unit, fmt_final) CIF_str(1:10)
    end do
    write(CIF_unit, '(a)') ''

   elseif(index(read_line, '_geom_hbond_site_symmetry_A') /=0) then
    if(ON_screen_CIF_item) then
     call write_info('  > _geom_hbond_                           features extracted from '// trim(input_CIF_file(1)))
    end if

    write(CIF_unit, '(a)') 'loop_'
    write(CIF_unit, '(a)') '    _geom_hbond_atom_site_label_D'
    write(CIF_unit, '(a)') '    _geom_hbond_atom_site_label_H'
    write(CIF_unit, '(a)') '    _geom_hbond_atom_site_label_A'
    write(CIF_unit, '(a)') '    _geom_hbond_distance_DH'
    write(CIF_unit, '(a)') '    _geom_hbond_distance_HA'
    write(CIF_unit, '(a)') '    _geom_hbond_distance_DA'
    write(CIF_unit, '(a)') '    _geom_hbond_angle_DHA'
    write(CIF_unit, '(a)') '    _geom_hbond_site_symmetry_A'
    do
     read(CIF_read_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:8)
     write(CIF_unit, '(3a5, 4a12,a6)') CIF_str(1:8)
    end do
    write(CIF_unit, '(a)') ''

   endif

  end do


  write(CIF_unit, '(a)') trim(CIF_sep_line)
  long = len_trim(sample_ID)
  !write(fmt_ , '(a,i2,a)') '(2a,', 53-long, 'x,2a)'
  !!write(CIF_unit, fmt=trim(fmt_))        '#===END OF CIF FILE FOR ', trim(u_case(sample_ID)),'#'
  !write(CIF_unit, trim(fmt_)) '#===END OF CIF FILE FOR ', u_case(trim(sample_ID)), (' ',i=1,53-long), '#'

  write(fmt_, '(a,i2,a)') '(2a,', 53-long, 'a1, a)'
  write(CIF_unit, trim(fmt_)) '#===END OF CIF FILE FOR ', u_case(trim(sample_ID)), (' ',i=1,53-long), '#'
  write(CIF_unit, '(a)') trim(CIF_sep_line)

  if(ON_screen_CIF_item) call write_info('')

  if(.not. write_details) call write_info('')
  if(include_HKL_file) then
  call write_info('   >>> CRYSCALC_ARCHIVE_hkl.cif file has been created.')
  else
  call write_info('   >>> CRYSCALC_ARCHIVE.cif file has been created.')
  end if
  call write_info('')

  if(CIF_parameter%abs_structure_flack /= '?') then
   call write_info('')
   call write_info("  !!! Please adjust the _chemical_absolute_configuration      'rm/ad/rmad/syn/unk/.' line.")
   call write_info('')
  end if

  stop

 return
end subroutine read_cif_and_create_ARCHIVE

!!----------------------------------------------------------------------------------------------------------!!
subroutine Include_RES_file_into_CIF
 use CRYSCALC_module, only  : input_CIF_file, CIF_unit, tmp_unit, include_HKL_file
 use CIF_module
 use macros_module,   only  : test_file_exist
 use IO_module,       only  : write_info
 implicit none
 integer                         :: i, i1, i_error, long
 integer                         :: pv
 CHARACTER (LEN=256)             :: job_RES_file, job_HKL_file, job_FAB_file
 LOGICAL                         :: file_exist
 CHARACTER (LEN=256)             :: read_line


   i1 = index(input_CIF_file(1), '.')
   if(i1 == 0) return

   job_RES_file = input_CIF_file(1)(1:i1)//'res'
   call test_file_exist(trim(job_RES_file), file_exist, 'out')

   job_HKL_file = input_CIF_file(1)(1:i1)//'hkl'

   if(file_exist) then
    call write_info('')
    call write_info('   ... Include '//trim(job_RES_file)//' SHELXL file.')
    !call write_info('')
    if(include_HKL_file) then
     write(CIF_unit, '(a,79a1)')  '#', ('.',i=1,79)
     write(CIF_unit, '(2a)') '# .RES SHELXL file: ', trim(job_RES_file)
     write(CIF_unit, '(2a)') '# .HKL SHELXL file: ', trim(job_HKL_file)
     job_FAB_file = input_CIF_file(1)(1:i1)//'fab'
     call test_file_exist(trim(job_FAB_file), file_exist, 'no_out')
     if(file_exist) then
     write(CIF_unit, '(2a)') '# .FAB SHELXL file: ', trim(job_FAB_file)
     end if
     write(CIF_unit, '(a,79a1)')  '#', ('.',i=1,79)
     write(CIF_unit, '(a)')  ''
     write(CIF_unit, '(a)') '_shelx_res_file'
    else
     call write_info('')
     !write(CIF_unit, '(a)') '_iucr_refine_instructions_details'
     write(CIF_unit, '(a)') '_shelx_res_file'
    end if
    write(CIF_unit, '(a)')   ';'
    open(unit=tmp_unit, file=trim(job_res_file))
    if(.not. include_HKL_file) then
     write(CIF_unit, '(2a)') ' .res SHELXL file: ', trim(job_RES_file)
     write(CIF_unit, '(a,68a1)')  '#', ('-',i=1,68)
    end if
    do
     read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     write(CIF_unit, '(a)') trim(read_line)
    end do
    close(unit=tmp_unit)
    write(CIF_unit, '(a)')   ';'
    if(CIF_parameter%shelx_res_checksum(1:1) /= '?') then
     write(CIF_unit, '(a)')   ''
     write(CIF_unit, '(2a)')   '_shelx_res_checksum              ', trim(CIF_parameter%shelx_res_checksum)
     write(CIF_unit, '(a)')   ''
    end if
    write(CIF_unit, '(2a,45a1)')   '# ............... End of .RES file ', ('.',i=1,45)

    !if(include_HKL_file) then
     job_HKL_file = input_CIF_file(1)(1:i1)//'hkl'
     if(include_HKL_file) then
      call test_file_exist(trim(job_HKL_file), file_exist, 'out')
     else
      call test_file_exist(trim(job_HKL_file), file_exist, 'no_out')
     end if
     if(file_exist) then
      if(include_HKL_file) then
       !call write_info('')
       call write_info('   ... Include '//trim(job_HKL_file)//' SHELXL file.')
       !call write_info('')
       write(CIF_unit, '(a)') ''
       write(CIF_unit, '(a)') '_shelx_hkl_file'
       !write(CIF_unit, '(a)') '_iucr_refine_reflections_details'
       write(CIF_unit, '(a)') ';'
       open(unit=tmp_unit, file=trim(job_HKL_file))
       !if(.not. include_HKL_file) then
       !write(CIF_unit, '(2a)') ' .HKL SHELXL file: ', trim(job_HKL_file)
       !write(CIF_unit, '(a,68a1)')  '#', ('-',i=1,68)
       !end if
       do
        read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
        if(i_error /=0) exit
        write(CIF_unit, '(a)') trim(read_line)
       end do
       close(unit=tmp_unit)
       write(CIF_unit, '(a)')   ';'
       if(CIF_parameter%shelx_hkl_checksum(1:1) /= '?') then
        write(CIF_unit, '(a)')   ''
        write(CIF_unit, '(2a)')   '_shelx_hkl_checksum              ', trim(CIF_parameter%shelx_hkl_checksum)
        write(CIF_unit, '(a)')   ''
       end if
       write(CIF_unit, '(2a,45a1)')   '# ............... End of .HKL file ', ('.',i=1,45)
      else
       write(CIF_unit, '(a)')   ''
       write(CIF_unit, '(a)')   '#'
       write(CIF_unit, '(2a)')  '# .HKL SHELXL file: ', trim(job_HKL_file)
       write(CIF_unit, '(a)')   '#'
      endif
     endif

     job_FAB_file = input_CIF_file(1)(1:i1)//'fab'
     !call test_file_exist(trim(job_FAB_file), file_exist, 'out')
     call test_file_exist(trim(job_FAB_file), file_exist, 'no_out')
     if(file_exist) then
      if(include_HKL_file) then
       !call write_info('')
       call write_info('   ... Include '//trim(job_FAB_file)//' SHELXL file.')
       call write_info('')
       write(CIF_unit, '(a)') ''
       write(CIF_unit, '(a)') '_shelx_fab_file'
       write(CIF_unit, '(a)') ';'
       open(unit=tmp_unit, file=trim(job_FAB_file))
       !if(.not. include_HKL_file) then
       !write(CIF_unit, '(2a)') ' .FAB SHELXL file: ', trim(job_FAB_file)
       !write(CIF_unit, '(a,68a1)')  '#', ('-',i=1,68)
       !end if
       pv = 0
       do
        read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
        if(i_error /=0) exit
        long = len_trim(read_line)
        !if(read_line(1:23) == '_platon_squeeze_details') then
        ! read_line = trim(read_line)//'                           '//'?'
        !end if
        !if(read_line(1:1) == ";") exit
        if(read_line(1:1) == ";") then
         pv = pv + 1
         read_line =')'
        end if
        write(CIF_unit, '(a)') trim(read_line)
        if(pv == 2) write(CIF_unit, '(a)') ';'
       end do
       close(unit=tmp_unit)
       if(pv==0) write(CIF_unit, '(a)')   ';'   ! necessaire pour l'ancien format du .fab qui ne contient pas de 'squeeze_details
       if(CIF_parameter%shelx_fab_checksum(1:1) /= '?') then
        write(CIF_unit, '(a)')   ''
        write(CIF_unit, '(2a)')   '_shelx_fab_checksum              ', trim(CIF_parameter%shelx_fab_checksum)
        write(CIF_unit, '(a)')   ''
       end if
       write(CIF_unit, '(2a,43a1)')   '# ............... End of .FAB file ', ('.',i=1,42),"#"
      else
       write(CIF_unit, '(a)')   ''
       write(CIF_unit, '(a)')   '#'
       write(CIF_unit, '(2a)')  '# .FAB SHELXL file: ', trim(job_FAB_file)
       write(CIF_unit, '(a)')   '#'
      endif
     else
      call write_info('')
     endif

    !end if
   endif

  return
end subroutine Include_RES_file_into_CIF


!!--------------------------------------------------------------------------------------------!!
subroutine GET_SG_from_CIF(input_file)
 USE cryscalc_module,                  ONLY : SPG, nb_symm_op, on_screen, CIF_read_unit
 USE macros_module,                    ONLY : replace_car
 USE CFML_String_Utilities,            ONLY : number_lines, reading_lines
 USE CFML_IO_Formats,                  ONLY : read_cif_hm, read_cif_symm, read_cif_hall
 USE CFML_Crystallographic_Symmetry,   ONLY : set_spacegroup, get_hallsymb_from_gener
 USE IO_module,                        ONLY : write_info

 implicit none
  CHARACTER (LEN=*), INTENT(IN)               :: input_file
  character(len=40)                           :: symb_spgr
  integer                                     :: npos, nb_lines, ier
  character(len=256),dimension(:),allocatable :: file_cif
  character(len=64), dimension(48)            :: car_symop

  close(unit=CIF_read_unit)
  call number_lines(trim(input_file),nb_lines)
  if (nb_lines ==0) then
   call write_info(' > No lines can be read. Program will be stopped')
   stop
  end if
  if (allocated(file_cif)) deallocate(file_cif)
  allocate(file_cif(nb_lines), stat=ier)
   if(ier/=0) call write_alloc_error("file_cif")
  file_cif=' '


  call reading_lines(trim(input_file), nb_lines, file_cif)

!---- OBTAIN SPACE GROUP ----!
  symb_spgr=' '
  npos=1
  call read_cif_hm(file_cif ,npos, nb_lines, symb_spgr)

  if(index(symb_spgr, '::R') /=0) symb_spgr = replace_car(symb_spgr, '::R', ':R')
  if(index(symb_spgr, ':H')  /=0) symb_spgr = symb_spgr(1:index(symb_spgr, ':H')-1)

  if (len_trim(symb_spgr) > 0) then
   !temp_string = SPG%Hall
   call Set_SpaceGroup(symb_spgr, SPG)
   !SPG%Hall = temp_string
  else
   npos=1
   call read_cif_hall(file_cif ,npos, nb_lines, symb_spgr)
   if (len_trim(symb_spgr) > 0) then
    call Set_SpaceGroup(symb_spgr, SPG)
   else
    npos=1
    call read_cif_symm(file_cif, npos ,nb_lines, nb_symm_op, car_symop)
    symb_spgr=' '
    !call set_spacegroup(symb_spgr, SPG, car_symop, nb_symm_op,'gen')
    !call get_hallsymb_from_gener(SPG)
    call set_spacegroup("  ", SPG, car_symop, nb_symm_op,'GEN')

    if(on_screen) then
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


  open(unit=CIF_read_unit, file= trim(input_file))

  return
end subroutine GET_SG_from_CIF

!!--------------------------------------------------------------------------------------------!!
subroutine GET_SG_op(input_file)
 USE cryscalc_module,                  ONLY : SPG, nb_symm_op, on_screen
 USE macros_module,                    ONLY : replace_car
 USE CFML_String_Utilities,            ONLY : number_lines, reading_lines
 USE CFML_IO_Formats,                  ONLY : read_cif_hm, read_cif_symm, read_cif_hall
 USE CFML_Crystallographic_Symmetry,   ONLY : set_spacegroup, get_hallsymb_from_gener
 USE IO_module,                        ONLY : write_info

 implicit none
  CHARACTER (LEN=*), INTENT(IN)               :: input_file
  character(len=40)                           :: symb_spgr
  integer                                     :: npos, nb_lines
  character(len=256),dimension(:),allocatable :: file_cif
  character(len=64), dimension(48)            :: car_symop

  call number_lines(trim(input_file),nb_lines)
  if (nb_lines ==0) then
   call write_info(' > No lines can be read. Program will be stopped')
   stop
  end if
  if (allocated(file_cif)) deallocate(file_cif)
  allocate(file_cif(nb_lines))
  file_cif=' '


  call reading_lines(trim(input_file), nb_lines, file_cif)

!---- OBTAIN SPACE GROUP ----!
  symb_spgr=' '
  npos=1
  call read_cif_hm(file_cif ,npos, nb_lines, symb_spgr)

  if(index(symb_spgr, '::R') /=0) symb_spgr = replace_car(symb_spgr, '::R', ':R')
  if(index(symb_spgr, ':H')  /=0) symb_spgr = symb_spgr(1:index(symb_spgr, ':H')-1)

  if (len_trim(symb_spgr) > 0) then
   call Set_SpaceGroup(symb_spgr, SPG)
  else
   npos=1
   call read_cif_hall(file_cif ,npos, nb_lines, symb_spgr)
   if (len_trim(symb_spgr) > 0) then
    call Set_SpaceGroup(symb_spgr, SPG)
   else
    npos=1
    call read_cif_symm(file_cif, npos ,nb_lines, nb_symm_op, car_symop)
    symb_spgr=' '
    !call set_spacegroup(symb_spgr, SPG, car_symop, nb_symm_op,'gen')
    !call get_hallsymb_from_gener(SPG)
    call set_spacegroup("  ", SPG, car_symop, nb_symm_op,'GEN')

    if(on_screen) then
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

  return
end subroutine GET_SG_op

!!--------------------------------------------------------------------------------------------!!
subroutine Date_format(date)
 use macros_module, only : check_character, replace_car
 implicit none
  character (len=*), intent(inout) :: date
  integer                          :: i, i1, i2, long
  logical                          :: alpha_char, numeric_char, previous_numeric
  logical, dimension(32)           :: change
  character (len=16)               :: mm, aa


  long = len_trim(date)

  previous_numeric = .true.
  i1 = 0
  i2 = 0
  do i=1, long
   change(i) = .false.
   call check_character(date(i:i), alpha_char, numeric_char)
    if((alpha_char    .and.       previous_numeric) .or.   &
      (numeric_char  .and. .not. previous_numeric)) then
      change(i) = .true.
      if(i1 ==0) then
       i1 = i
      elseif(i2 == 0) then
       i2 = i
       exit
      end if
   end if
   previous_numeric = .false.
   if(numeric_char)  previous_numeric = .true.
  end do

  !date = replace_car(date, date(i1:i1), ' '//date(i1:i1))
  !date = replace_car(date, date(i2+1:i2+1), ' '//date(i2+1:i2+1))

  if(i2 == 0) return    ! le repertoire ne contient pas la date

  mm = date(i1:i2-1)
  call Check_date_mois(mm)
  aa = date(i2:)
  call Check_date_an(aa)
  date = date(1:i1-1) // ' ' //trim(mm)//' '// trim(aa)


 return
end subroutine Date_format

!!--------------------------------------------------------------------------------------------!!
subroutine Date_format_2(date)
! transformation date format "jj/mm/aa" en "jj mois annee"
 use macros_module, only : check_character, replace_car
 implicit none
  character (len=*), intent(inout) :: date
  integer                          :: i1, i2, long
  character (len=16)               :: jj, mm, aa


  long = len_trim(date)

  i1 = index(date, "/")
  i2 = index(date, "/", back=.true.)

  jj = date(1:i1-1)
  mm = date(i1+1:i2-1)
  aa = date(i2+1:)

  if(mm(1:2) == '01') then
   mm = 'janvier'
  elseif(mm(1:2) == '02') then
   mm = 'fevrier'
  elseif(mm(1:2) == '03') then
   mm = 'mars'
  elseif(mm(1:2) == '04') then
   mm = 'avril'
  elseif(mm(1:2) == '05') then
   mm = 'mai'
  elseif(mm(1:2) == '06') then
   mm = 'juin'
  elseif(mm(1:2) == '07') then
   mm = 'juillet'
  elseif(mm(1:2) == '08') then
   mm = 'aout'
  elseif(mm(1:2) == '09') then
   mm = 'septembre'
  elseif(mm(1:2) == '10') then
   mm = 'octobre'
  elseif(mm(1:2) == '11') then
   mm = 'novembre'
  elseif(mm(1:2) == '12') then
   mm = 'decembre'
  end if

  if(aa(1:2) == "20") then
   date = date(1:i1-1) // ' ' //trim(mm)//' '// trim(aa)
  else
   date = date(1:i1-1) // ' ' //trim(mm)//' '// "20"//trim(aa)
  end if

 return
end subroutine Date_format_2

!!--------------------------------------------------------------------------------------------!!
subroutine Check_date_mois(mm)
 implicit none
   character(len=*), intent(inout) :: mm
  integer                          :: long


  long = len_trim(mm)
  if(long ==3) then
   if(mm(1:3) == 'jan') then
    mm = 'janvier'
   elseif(mm(1:3) == 'fev') then
    mm = 'fevrier'
   elseif(mm(1:3) == 'sep') then
    mm = 'septembre'
   elseif(mm(1:3) == 'oct') then
    mm = 'octobre'
   elseif(mm(1:3) == 'nov') then
    mm = 'novembre'
   elseif(mm(1:3) == 'dec') then
    mm = 'decembre'
   end if
  elseif(long == 4) then
   if(mm(1:4) == 'janv') then
    mm = 'janvier'
   elseif(mm(1:4) == 'juil') then
    mm = 'juillet'
   elseif(mm(1:4) == 'sept') then
    mm = 'septembre'
   end if
  end if

 return
end subroutine Check_date_mois

!!--------------------------------------------------------------------------------------------!!
subroutine Check_date_an(aa)
 implicit none
  character (len=*), intent(inout) :: aa
  integer                          :: long

  long = len_trim(aa)
  if (long == 2) aa = '20'//aa(1:2)

  return
end subroutine Check_date_an

!!--------------------------------------------------------------------------------------------!!
subroutine Get_sample_ID_from_import_cif(sample_ID)
 use CRYSCALC_module, only : CIF_read_unit
 implicit none
 character (len=256), intent(inout)    :: sample_ID
 character (len=256)                   :: read_line
 integer                               :: i_error
 integer                               :: i1, i2


  open(unit=CIF_read_unit,  file='import.cif')

  do
   read(CIF_read_unit, '(a)', iostat=i_error) read_line
   if(i_error /=0) exit
   if(len_trim(read_line) == 0) cycle
   i1 = index(read_line, "Sample ID")
   if( i1/=0) then
    i1 = index(read_line, ":")
    i2 = index(read_line, "#", back=.true.)
    if(i1 /=0 .and. i2/=0 .and. i2 > i1 +1) then
     read(read_line(i1+1:i2-1), '(a)') sample_ID
     sample_ID = trim(sample_ID)
     sample_ID = adjustl(sample_ID)
     exit
    end if
   end if

  end do

  close(unit=CIF_read_unit)

 return
end subroutine Get_sample_ID_from_import_cif

!!--------------------------------------------------------------------------------------------!!
subroutine Get_items_from_import_cif_header(diffraction_date, experimenter, proposer_name, sample_ID)
 use CRYSCALC_module, only : CIF_read_unit
 implicit none
 character (len=256), intent(inout)    :: diffraction_date, experimenter, Proposer_name, sample_ID
 character (len=256)                   :: read_line
 integer                               :: i_error
 integer                               :: i1, i2

  open(unit=CIF_read_unit,  file='import.cif')

  do
   read(CIF_read_unit, '(a)', iostat=i_error) read_line
   if(i_error /=0) exit
   if(len_trim(read_line) == 0) cycle


   i1 = index(read_line, "Date of experiment")
   if( i1/=0) then
    i1 = index(read_line, ":")
    i2 = index(read_line, "#", back=.true.)
    if(i1 /=0 .and. i2/=0 .and. i2 > i1 +1) then
     read(read_line(i1+1:i2-1), '(a)') diffraction_date
     diffraction_date = trim(diffraction_date)
     diffraction_date = adjustl(diffraction_date)
    end if
   end if

   i1 = index(read_line, "Experimenter")
   if( i1/=0) then
    i1 = index(read_line, ":")
    i2 = index(read_line, "#", back=.true.)
    if(i1 /=0 .and. i2/=0 .and. i2 > i1 +1) then
     read(read_line(i1+1:i2-1), '(a)') experimenter
     experimenter = trim(experimenter)
     experimenter = adjustl(experimenter)
    end if
   end if


   i1 = index(read_line, "Proposer")
   if( i1/=0) then
    i1 = index(read_line, ":")
    i2 = index(read_line, "#", back=.true.)
    if(i1 /=0 .and. i2/=0 .and. i2 > i1 +1) then
     read(read_line(i1+1:i2-1), '(a)') Proposer_name
     Proposer_name = trim(Proposer_name)
     Proposer_name = adjustl(Proposer_name)
    end if
   end if


   i1 = index(read_line, "Sample ID")
   if( i1/=0) then
    i1 = index(read_line, ":")
    i2 = index(read_line, "#", back=.true.)
    if(i1 /=0 .and. i2/=0 .and. i2 > i1 +1) then
     read(read_line(i1+1:i2-1), '(a)') sample_ID
     sample_ID = trim(sample_ID)
     sample_ID = adjustl(sample_ID)
     exit
    end if
   end if

  end do

  close(unit=CIF_read_unit)

 return
end subroutine Get_items_from_import_cif_header

!!--------------------------------------------------------------------------------------------!!
subroutine Get_sample_ID_from_CIF(input_file, sample_ID)
 use CRYSCALC_module, only : CIF_read_unit
 implicit none
 character (len=256), intent(inout)    :: input_file, sample_ID
 character (len=256)                   :: read_line
 integer                               :: i_error
 integer                               :: i1

  open(unit=CIF_read_unit,  file=trim(input_file))

  do
   read(CIF_read_unit, '(a)', iostat=i_error) read_line
   if(i_error /=0) exit
   if(len_trim(read_line) == 0) cycle

   i1 = index(read_line, "data_")
   if( i1/=0) then
    read(read_line(i1+5:), '(a)') sample_ID
    sample_ID = trim(sample_ID)
    sample_ID = adjustl(sample_ID)
    exit
   end if
  end do

  close(unit=CIF_read_unit)

 return
end subroutine Get_sample_ID_from_CIF

!!--------------------------------------------------------------------------------------------!!
 subroutine Check_formula(formula)
! verif. si stoechiometrie nulle
!   ex : transforme C12 H12 Cl0 N4 O2 en C12 H12 N4 O2

  use macros_module, only : nombre_de_colonnes, check_character, replace_car
  implicit none
  character (len=64), intent(inout) :: formula
  integer                           :: i, j, nb_col, i_error, long
  character (len=8), dimension(50)  :: string
  LOGICAL                           :: alpha_char, numeric_char
  real                              :: stoechio
  real, parameter                   :: eps = 0.0001


  call nombre_de_colonnes(formula, nb_col)
  READ(formula, *, IOSTAT=i_error) string(1:nb_col)
  do i=1, nb_col
   long = len_trim(string(i))
   do j=1, long
    numeric_char = .false.
    alpha_char   = .false.
    call check_character(string(i)(j:j), alpha_char, numeric_char)

    if(numeric_char) then
     read(string(i)(j:long), *) stoechio
     if (stoechio < eps) string(i) = ""
     exit
    end if
   end do
  end do

  formula = ''
  do i=1, nb_col
   formula = trim(formula)//' '//trim(string(i))
  end do
  formula = adjustl(formula)

  return
 end subroutine Check_formula
!--------------------------------------------------------------------------------------------
