!     Last change:  TR   12 Oct 2007    6:03 pm
!subroutine write_CIF_file(input_string, n_value, input_value)
subroutine write_CIF_file(input_string)
 USE cryscalc_module, ONLY : CIF_unit, nb_symm_op, symm_op_string, crystal, F000,             &
                             keyword_CELL, keyword_WAVE, keyword_SIZE, wavelength, Z_unit,    &
                             keyword_create_CIF, keyword_create_archive, experiment,          &
                             CRYSCALC, SAINT, SADABS, TWINABS_file, SHELX,                    &
                             ABSORPTION, ABS_CRYSALIS, SPG, unit_cell, molecule,              &
                             nb_atom, atom_label, atom_typ, atom_coord, atom_Ueq, atom_occ,   &
                             atom_adp_type, atom_occ_perc, atom_adp_aniso,                    &
                             device, AUTHOR, UB_matrix, P4P_file_name, SQUEEZE,               &
                             keyword_modif_ARCHIVE, ABS_file_name, RAW_file_name, debug_proc, &
                             EXP_SCAN
 USE CIF_module
 USE text_module, only     : CIF_lines_nb, CIF_title_line
 USE hkl_module
 USE macros_module, only   : test_file_exist, l_case, u_case, Get_current_folder_name, Get_Proposer_from_current_folder_name, &
                             replace_car

 implicit none
  CHARACTER (LEN=*),                     INTENT(IN) :: input_string
  INTEGER                                           :: i
  CHARACTER (LEN=36), DIMENSION(7)                  :: CELL_param_CIF_string
  CHARACTER (LEN=6)                                 :: adp_type



 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_CIF_file ("//trim(input_string)//")")

  cell_param_CIF_string(1) =  "_cell_length_a                      "
  cell_param_CIF_string(2) =  "_cell_length_b                      "
  cell_param_CIF_string(3) =  "_cell_length_c                      "
  cell_param_CIF_string(4) =  "_cell_angle_alpha                   "
  cell_param_CIF_string(5) =  "_cell_angle_beta                    "
  cell_param_CIF_string(6) =  "_cell_angle_gamma                   "
  cell_param_CIF_string(7) =  "_cell_volume                        "

  select case (input_string)

      case ('ABSORPTION')
  CALL WRITE_CIF_section("ABSORPTION CORRECTION")


  !if(len_trim(SADABS%absorption_coeff) /=0) then
  if(SADABS%absorption_coeff(1:1) /= '?') then
   WRITE(CIF_unit, '(2a)')      "_exptl_absorpt_coefficient_mu                         ", trim(SADABS%absorption_coeff)
   !WRITE(CIF_unit, '(2a)')      "_exptl_absorpt_coefficient_mu                           ", trim(CIF_parameter%absorption_coefficient_mu)

  else
   IF(ABSORPTION%mu >  0.) then
   WRITE(CIF_unit, '(a,F7.3)')   "_exptl_absorpt_coefficient_mu                         ", ABSORPTION%mu * 0.1
   else
   WRITE(CIF_unit, '(a)')        "_exptl_absorpt_coefficient_mu                        ?"
   endif
  endif


      case ('TMIN')
  IF(ABSORPTION%Tmin > 0.) then
   if(keyword_modif_archive) then
    if(ABSORPTION%Tmax > 0. .and. SADABS%ratio > 0.) then
     WRITE(CIF_unit, '(a,F7.3)') "_exptl_absorpt_correction_T_min                       ", ABSORPTION%Tmax*SADABS%ratio
    else
     WRITE(CIF_unit, '(a,F7.3)') "_exptl_absorpt_correction_T_min                       ", ABSORPTION%Tmin*0.985
    endif
   else
    WRITE(CIF_unit, '(a,F7.3)')  "_exptl_absorpt_correction_T_min                       ", ABSORPTION%Tmin
   endif
  else
   WRITE(CIF_unit, '(a)')        "_exptl_absorpt_correction_T_min                      ?"
  endif

      case ('TMAX')
  IF(ABSORPTION%Tmax > 0.) then
  WRITE(CIF_unit, '(a,F7.3)')    "_exptl_absorpt_correction_T_max                       ", ABSORPTION%Tmax
  else
  WRITE(CIF_unit, '(a)')         "_exptl_absorpt_correction_T_max                      ?"
  endif
  IF(.not. keyword_create_ARCHIVE) then
  !IF(LEN_TRIM(SADABS%line_ratio) /=0) then
  IF(SADABS%line_ratio(1:1) /= '?') then
   WRITE(CIF_unit, '(a)')  ""
   WRITE(CIF_unit, '(a)')  "#----------------------------- Remark ---------------------------------#"
   WRITE(CIF_unit, '(a)')  "# Tmax and Tmin values correspond to EXPECTED values calculated        #"
   WRITE(CIF_unit, '(a)')  "# from crystal size. In the case of absorption correction performed    #"
   WRITE(CIF_unit, '(a)')  "# with SADABS program, Tmax should be given as Tmax_expected and       #"
   WRITE(CIF_unit, '(a)')  "# and Tmin = Tmax * 'relative-correction-factor'.                      #"
   WRITE(CIF_unit, '(a)')  "# SADABS output:                                                       #"
   !WRITE(CIF_unit, '(3a)') "# ", TRIM(SADABS%line_ratio), "         #"
   call write_cif_line(trim(SADABS%line_ratio)) 
   WRITE(CIF_unit, '(a)')  "#----------------------------------------------------------------------#"
  ELSEIF(LEN_TRIM(SADABS%line_estimated_Tmin_Tmax) /=0) then
   if(SADABS%line_estimated_Tmin_Tmax(1:1) /= '?') then
   WRITE(CIF_unit, '(a)')  ""
   WRITE(CIF_unit, '(a)')  "#------------------------------- Remark -------------------------------------#"
   if(TWINABS_file) then
   WRITE(CIF_unit, '(a)')  "# TWINABS output:                                                            #"
   else
   WRITE(CIF_unit, '(a)')  "# SADABS output:                                                             #"
   endif
   !WRITE(CIF_unit, '(3a)') "# ", TRIM(SADABS%line_estimated_Tmin_Tmax), "                #"
   call write_cif_line(trim(SADABS%line_estimated_Tmin_Tmax))   
   WRITE(CIF_unit, '(a)')  "# The ratio of these values is more reliable than their absolute values!     #"
   WRITE(CIF_unit, '(a)')  "#----------------------------------------------------------------------------#"
   end if
  endif
  end if

      case ('SADABS')
  WRITE(CIF_unit, '(a)') trim(SADABS%type)
  if(.not. TWINABS_file) then
   do i=1, 4
    WRITE(CIF_unit, '(a)') trim(SADABS%details(i))
   end do
  else
   do i=1, 4
    WRITE(CIF_unit, '(a)') trim(SADABS%tw_details(i))
   end do
  end if

       case ('SHELX')
  do i=1, 6
   WRITE(CIF_unit, '(a)') trim(SHELX%details(i))
  end do


      case ('ABS_CRYSALIS')
  WRITE(CIF_unit, '(a)') trim(ABS_CRYSALIS%type)
  do i=1, 6
   WRITE(CIF_unit, '(a)') trim(ABS_CRYSALIS%details(i))
  end do

      case ('SQUEEZE')
  call write_CIF_section("SQUEEZE RESULTS")
  do i=2, SQUEEZE%nb_lines
   write(CIF_unit, '(a)') trim(SQUEEZE%read_line(i))
  end do

          case ('DATACOL')
  CALL write_CIF_section("DATA COLLECTION")
  WRITE(CIF_unit, '(2a)') "_diffrn_measurement_device_type        ", trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
  WRITE(CIF_unit, '(2a)') "_diffrn_measurement_method             ", trim(CIF_parameter_DEVICE%diffrn_measurement_method)
  WRITE(CIF_unit, '(2a)') "_diffrn_detector                       ", trim(CIF_parameter_DEVICE%diffrn_detector)

          case ('COMPUTER_PROGRAMS')
  call write_CIF_section("COMPUTER PROGRAMS USED")
  if(keyword_create_archive) then
   WRITE(CIF_unit, '(2a)') "_computing_data_collection       ", trim(CIF_parameter_DEVICE%computing_data_reduction)
   WRITE(CIF_unit, '(2a)') "_computing_cell_refinement       ", trim(CIF_parameter_DEVICE%computing_cell_refinement)
   WRITE(CIF_unit, '(2a)') "_computing_data_reduction        ", trim(CIF_parameter_DEVICE%computing_data_reduction)
  end if
  WRITE(CIF_unit, '(3a)') "_computing_structure_solution    '",trim(CIF_parameter_DEVICE%computing_structure_solution), "'"
  WRITE(CIF_unit, '(2a)') "_computing_structure_refinement  ", trim(CIF_parameter_DEVICE%computing_structure_refinement)
  WRITE(CIF_unit, '(2a)') "_computing_molecular_graphics    ", trim(CIF_parameter_DEVICE%computing_molecular_graphics)
  WRITE(CIF_unit, '(3a)') "_computing_publication_material  '", trim(CIF_parameter_DEVICE%computing_publication_material_2),"'"

      case ('X2S')
  call write_CIF_section("DATA COLLECTION")
  WRITE(CIF_unit, '(2a)') "_diffrn_measurement_device_type    ", trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
  WRITE(CIF_unit, '(2a)') "_diffrn_measurement_method         ", trim(CIF_parameter_DEVICE%diffrn_measurement_method)
  WRITE(CIF_unit, '(2a)') "_diffrn_detector                   ", trim(CIF_parameter_DEVICE%diffrn_detector)
  WRITE(CIF_unit, '(2a)') "_diffrn_detector_area_resol_mean   ", trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean)
!  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_wavelength       ", trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength)
!  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_type             ", trim(CIF_parameter_DEVICE%diffrn_radiation_type)
  !WRITE(CIF_unit, '(2a)') "_diffrn_radiation_probe            ", trim(CIF_parameter_DEVICE%diffrn_radiation_probe)
  !WRITE(CIF_unit, '(2a)') "_diffrn_radiation_source           ", trim(CIF_parameter_DEVICE%diffrn_radiation_source)
  !WRITE(CIF_unit, '(2a)') "_diffrn_radiation_monochromator    ", trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator)

        case ('APEX', 'D8V', 'D8V_MO', 'D8V_CU')
  call write_CIF_section("DATA COLLECTION")
  WRITE(CIF_unit, '(2a)') "_diffrn_measurement_device_type    ", trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
  WRITE(CIF_unit, '(2a)') "_diffrn_measurement_method         ", trim(CIF_parameter_DEVICE%diffrn_measurement_method)
  WRITE(CIF_unit, '(2a)') "_diffrn_detector                   ", trim(CIF_parameter_DEVICE%diffrn_detector)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_wavelength       ", trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_type             ", trim(CIF_parameter_DEVICE%diffrn_radiation_type)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_probe            ", trim(CIF_parameter_DEVICE%diffrn_radiation_probe)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_source           ", trim(CIF_parameter_DEVICE%diffrn_radiation_source)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_monochromator    ", trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator)

  call write_CIF_section("COMPUTER PROGRAMS USED")
  WRITE(CIF_unit, '(2a)') "_computing_data_collection       ", trim(CIF_parameter_DEVICE%computing_data_reduction)
  WRITE(CIF_unit, '(2a)') "_computing_cell_refinement       ", trim(CIF_parameter_DEVICE%computing_cell_refinement)
  WRITE(CIF_unit, '(2a)') "_computing_data_reduction        ", trim(CIF_parameter_DEVICE%computing_data_reduction)
  !WRITE(CIF_unit, '(2a)') "_computing_structure_solution    ", trim(CIF_parameter_DEVICE%computing_structure_solution)
  !WRITE(CIF_unit, '(2a)') "_computing_structure_refinement  ", trim(CIF_parameter_DEVICE%computing_structure_refinement)
  !WRITE(CIF_unit, '(2a)') "_computing_molecular_graphics    ", trim(CIF_parameter_DEVICE%computing_molecular_graphics)
  !WRITE(CIF_unit, '(2a)') "_computing_publication_material  ", trim(CIF_parameter_DEVICE%computing_publication_material)


      case ('KCCD')
  call write_CIF_section("DATA COLLECTION")
  WRITE(CIF_unit, '(2a)')"_diffrn_source                         ", trim(CIF_parameter_DEVICE%diffrn_source)
  WRITE(CIF_unit, '(2a)')"_diffrn_radiation_wavelength           ", trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength)
  WRITE(CIF_unit, '(2a)')"_diffrn_radiation_type                 ", trim(CIF_parameter_DEVICE%diffrn_radiation_type)
  WRITE(CIF_unit, '(2a)')"_diffrn_radiation_probe                ", trim(CIF_parameter_DEVICE%diffrn_radiation_probe)
  WRITE(CIF_unit, '(2a)')"_diffrn_radiation_source               ", trim(CIF_parameter_DEVICE%diffrn_radiation_source)
  WRITE(CIF_unit, '(2a)')"_diffrn_radiation_monochromator        ", trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator)
  WRITE(CIF_unit, '(2a)')"_diffrn_detector                       ", trim(CIF_parameter_DEVICE%diffrn_detector)
  WRITE(CIF_unit, '(2a)')"_diffrn_detector_area_resol_mean       ", trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean)
  WRITE(CIF_unit, '(2a)')"_diffrn_measurement_device             ", trim(CIF_parameter_DEVICE%diffrn_measurement_device)
  WRITE(CIF_unit, '(2a)')"_diffrn_measurement_device_type        ", trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
  WRITE(CIF_unit, '(2a)')"_diffrn_measurement_method             ", trim(CIF_parameter_DEVICE%diffrn_measurement_method)

        case ('SAPHIRE', 'XCALIBUR', 'SUPERNOVA')
  call write_CIF_section("DATA COLLECTION")
  WRITE(CIF_unit, '(a)') "_diffrn_measurement_device_type    ", trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
  WRITE(CIF_unit, '(a)') "_diffrn_measurement_method         ", trim(CIF_parameter_DEVICE%diffrn_measurement_method)
  WRITE(CIF_unit, '(a)') "_diffrn_detector                   ", trim(CIF_parameter_DEVICE%diffrn_detector)
  WRITE(CIF_unit, '(a)') "_diffrn_detector_area_resol_mean   ", trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean)

  call write_CIF_section("COMPUTER PROGRAMS USED")
  WRITE(CIF_unit, '(a)') "_computing_data_collection      "
  WRITE(CIF_unit, '(10x,a)')                      trim(CIF_parameter_DEVICE%computing_data_collection)
  WRITE(CIF_unit, '(a)') "_computing_cell_refinement      "
  WRITE(CIF_unit, '(10x,a)')                      trim(CIF_parameter_DEVICE%computing_cell_refinement)
  WRITE(CIF_unit, '(a)') "_computing_data_reduction       "
  WRITE(CIF_unit, '(10x,a)')                      trim(CIF_parameter_DEVICE%computing_data_reduction)

  !WRITE(CIF_unit, '(2a)') "_computing_structure_solution   ", trim(CIF_parameter_DEVICE%computing_structure_solution)
  !WRITE(CIF_unit, '(2a)') "_computing_structure_refinement ", trim(CIF_parameter_DEVICE%computing_structure_refinement)
  !WRITE(CIF_unit, '(2a)') "_computing_molecular_graphics   ", trim(CIF_parameter_DEVICE%computing_molecular_graphics)
  !WRITE(CIF_unit, '(a)')  "_computing_publication_material "
  !WRITE(CIF_unit, '(a)')  ";"
  !WRITE(CIF_unit, '(10x,a)')  trim(CIF_parameter_DEVICE%computing_publication_material_1)
  !WRITE(CIF_unit, '(10x,a)')  trim(CIF_parameter_DEVICE%computing_publication_material_2)
  !WRITE(CIF_unit, '(a)')  ";"




    case ('DATA_LIMITS')
  WRITE(CIF_unit, '(a)')  ""
  WRITE(CIF_unit, '(a,I6)') "_diffrn_reflns_number                   ", n_ref
  WRITE(CIF_unit, '(a,I6)') "_diffrn_reflns_limit_h_min              ", MINVAL(h(1:n_ref))
  WRITE(CIF_unit, '(a,I6)') "_diffrn_reflns_limit_h_max              ", MAXVAL(h(1:n_ref))
  WRITE(CIF_unit, '(a,I6)') "_diffrn_reflns_limit_k_min              ", MINVAL(k(1:n_ref))
  WRITE(CIF_unit, '(a,I6)') "_diffrn_reflns_limit_k_max              ", MAXVAL(k(1:n_ref))
  WRITE(CIF_unit, '(a,I6)') "_diffrn_reflns_limit_l_min              ", MINVAL(l(1:n_ref))
  WRITE(CIF_unit, '(a,I6)') "_diffrn_reflns_limit_l_max              ", MAXVAL(l(1:n_ref))
  IF(keyword_CELL .AND. keyword_WAVE) then
   WRITE(CIF_unit, '(a,F6.2)') "_diffrn_reflns_theta_min                ", MINVAL(theta_hkl(1:n_ref))
   WRITE(CIF_unit, '(a,F6.2)') "_diffrn_reflns_theta_max                ", MAXVAL(theta_hkl(1:n_ref))
  endif

    case ('HKL_DATA')
  call write_CIF_section("HKL DATA")
  WRITE(CIF_unit, '(a)')  "loop_"
  WRITE(CIF_unit, '(a)')  "_refln_index_h"
  WRITE(CIF_unit, '(a)')  "_refln_index_k"
  WRITE(CIF_unit, '(a)')  "_refln_index_l"
  WRITE(CIF_unit, '(a)')  "_refln_F_squared_meas"
  WRITE(CIF_unit, '(a)')  "_refln_F_squared_sigma"
  WRITE(CIF_unit, '(a)')  "_refln_scale_group_code"
  do i= 1, n_ref
   WRITE(CIF_unit, '(i4,2i5,2F9.2,I5)')      h(i), k(i), l(i), F2(i), sig_F2(i), code(i)   ! sans les cos. dir.
  end do




      case ('EVAL_PROGRAMS', 'DENZO_PROGRAMS', 'CRYSALIS_PROGRAMS')
  call write_CIF_section("COMPUTER PROGRAMS USED")
  WRITE(CIF_unit, '(2a)') "_computing_data_collection     ", trim(CIF_parameter_DEVICE%computing_data_collection)
  WRITE(CIF_unit, '(2a)') "_computing_cell_refinement     ", trim(CIF_parameter_DEVICE%computing_cell_refinement)
  WRITE(CIF_unit, '(2a)') "_computing_data_reduction      ", trim(CIF_parameter_DEVICE%computing_data_reduction)





      case ('CRYSTAL_INFORMATION')
  call write_CIF_section("CRYSTAL INFORMATION")
  WRITE(CIF_unit, '(2a)')       "_exptl_crystal_description             ", TRIM(crystal%morph)
  WRITE(CIF_unit, '(2a)')       "_exptl_crystal_colour                  ", TRIM(crystal%color)
  IF(keyword_SIZE) then
   WRITE(CIF_unit, '(a,F6.3)')  "_exptl_crystal_size_max                ", crystal%size_max
   WRITE(CIF_unit, '(a,F6.3)')  "_exptl_crystal_size_mid                ", crystal%size_mid
   WRITE(CIF_unit, '(a,F6.3)')  "_exptl_crystal_size_min                ", crystal%size_min
  else
   WRITE(CIF_unit, '(a)')       "_exptl_crystal_size_max                ?"
   WRITE(CIF_unit, '(a)')       "_exptl_crystal_size_mid                ?"
   WRITE(CIF_unit, '(a)')       "_exptl_crystal_size_min                ?"
  endif
  if(crystal%faces_nb /=0) then
   WRITE(CIF_unit, '(a)')       ""
   WRITE(CIF_unit, '(a)')       "loop_"
   WRITE(CIF_unit, '(a)')       "_exptl_crystal_face_index_h"
   WRITE(CIF_unit, '(a)')       "_exptl_crystal_face_index_k"
   WRITE(CIF_unit, '(a)')       "_exptl_crystal_face_index_l"
   WRITE(CIF_unit, '(a)')       "_exptl_crystal_face_perp_dist"
   do i=1, crystal%faces_nb
    WRITE(CIF_unit, '(a)') trim(crystal%face_line(i))
   end do
  endif
  if(.not. keyword_create_archive) WRITE(CIF_unit, '(a)')        ""

      case ('F000')
  WRITE(CIF_unit, '(a, F8.1)')  "_exptl_crystal_F_000                          ", F000

      case ('CHEMICAL_INFORMATION')
  call write_CIF_section("CHEMICAL INFORMATION")


      case('CHEMICAL_2')
  IF(molecule%common_name(1:1) /= '?') then
   WRITE(CIF_unit, '(a)')       ""
   WRITE(CIF_unit, '(3a)')      "_chemical_name_common                  '", TRIM(molecule%common_name),"'"
  endif
  WRITE(CIF_unit, '(3a)')       "_chemical_formula_sum                  '",  trim(adjustl(molecule%formula)),"'"
  IF(molecule%formula(1:1) /= '?') then
   !WRITE(CIF_unit, '(3a)')      "_chemical_formula_sum                  '", adjustl(trim(molecule%formula)),"'"
   WRITE(CIF_unit, '(a, F8.2)') "_chemical_formula_weight                ",  Molecule%weight
   !WRITE(CIF_unit, '(a)')       " "
  else
   WRITE(CIF_unit, '(a)') "_chemical_formula_weight                ?"
  endif
  WRITE(CIF_unit, '(a)')       " "

      case ('CRYSTAL_DENSITY')
  WRITE(CIF_unit, '(2a)') "_exptl_crystal_density_meas           ",  crystal%density_meas
  WRITE(CIF_unit, '(2a)') "_exptl_crystal_density_diffrn         ",  crystal%density_diffrn
  WRITE(CIF_unit, '(2a)') "_exptl_crystal_density_method         ",  crystal%density_method
  WRITE(CIF_unit, '(2a)') "_exptl_crystal_F_000                  ",  crystal%F000
  WRITE(CIF_unit, '(a)')  "_exptl_special_details"
  WRITE(CIF_unit, '(a)')  ";"
  WRITE(CIF_unit, '(a)')  " ?"
  WRITE(CIF_unit, '(a)')  ";"

      case ('WAVE')
  WRITE(CIF_unit, '(a, F10.5)') "_diffrn_radiation_wavelength            ", wavelength

      case ('CELL_WAVE')
  WRITE(CIF_unit, '(a, F10.5)') "_cell_measurement_wavelength            ", CIF_cell_measurement%wavelength

      case ('ZUNIT')
  WRITE(CIF_unit, '(a, F4.1)') "_cell_formula_units_Z                    ", Z_unit


      case ('UNIT_CELL_INFO')
 call write_CIF_section("UNIT CELL INFORMATION")


      case ('SPACE_GROUP')

  WRITE(CIF_unit, '(a)')     ""
  if(SPG%Numspg ==0) then
   !WRITE(CIF_unit, '(a)')    "_symmetry_cell_setting                   '?'"
   !WRITE(CIF_unit, '(a)')    "_symmetry_space_group_name_H-M           '?'"
   !WRITE(CIF_unit, '(a)')    "_symmetry_space_group_name_Hall          '?'"
   !WRITE(CIF_unit, '(a)')    "_symmetry_Int_Tables_number              '?'"
   WRITE(CIF_unit, '(a)')    "_space_group_crystal_system              '?'"
   WRITE(CIF_unit, '(a)')    "_space_group_name_H-M_alt                '?'"
   WRITE(CIF_unit, '(a)')    "_space_group_name_Hall                   '?'"
   WRITE(CIF_unit, '(a)')    "_space_group_IT_number                   '?'"
  else
   !WRITE(CIF_unit, '(2a)')    "_symmetry_cell_setting                   ", TRIM(SPG%CrystalSys)
   !WRITE(CIF_unit, '(3a)')    "_symmetry_space_group_name_H-M          '", TRIM(SPG%SPG_Symb),  "'"
   !WRITE(CIF_unit, '(3a)')    "_symmetry_space_group_name_Hall         '", TRIM(SPG%Hall),      "'"
   !WRITE(CIF_unit, '(a,I3)')  "_symmetry_Int_Tables_number              ", SPG%NumSpg
   WRITE(CIF_unit, '(2a)')    "_space_group_crystal_system              ", l_case(TRIM(SPG%CrystalSys))
   WRITE(CIF_unit, '(3a)')    "_space_group_name_H-M_alt                '", TRIM(SPG%SPG_Symb),  "'"
   !WRITE(CIF_unit, '(3a)')    "_space_group_name_H-M_alt                '", TRIM(CIF_parameter%symmetry_space_group),  "'"

   WRITE(CIF_unit, '(3a)')    "_space_group_name_Hall                   '", TRIM(adjustl(SPG%Hall)),      "'"
   WRITE(CIF_unit, '(a,I3)')  "_space_group_IT_number                   ", SPG%NumSpg
  endif
  WRITE(CIF_unit, '(a)')     "loop_"
  !WRITE(CIF_unit, '(a)')     "  _symmetry_equiv_pos_as_xyz"
  WRITE(CIF_unit, '(a)')     "  _space_group_symop_operation_xyz"
  if(keyword_create_ARCHIVE) then
   do i=1, SPG%Multip   !
    !write(CIF_unit, '(3a)') "'", trim(SPG%SymopSymb(i)), "'"
    if(adjustl(CIF_parameter%sym_op(i)(1:1)) /= '?') then
     write(CIF_unit, '(a)') trim(CIF_parameter%sym_op(i))
    else
     write(CIF_unit, '(3a)') "'", trim(SPG%SymopSymb(i)), "'"
    endif
   end do

  else
   if(SPG%NumSpg /=0) then
    do i=1, SPG%Multip   !
     write(CIF_unit, '(3a)') "'", trim(SPG%SymopSymb(i)), "'"
    end do
   else
    do i=1, nb_symm_op
     WRITE(CIF_unit, '(7a)')  "'", TRIM(symm_op_string(1,i)), ", ", TRIM(symm_op_string(2,i)), ", ", TRIM(symm_op_string(3,i)), "'"
    end do
   end if
  end if


  WRITE(CIF_unit, '(a)')     ""


 case ('SPACE_GROUP_PYMOL')
  if(SPG%Numspg ==0) then
   WRITE(CIF_unit, '(a)')    "_symmetry_cell_setting                   '?'"
   WRITE(CIF_unit, '(a)')    "_symmetry_space_group_name_H-M           '?'"
   WRITE(CIF_unit, '(a)')    "_symmetry_space_group_name_Hall          '?'"
   WRITE(CIF_unit, '(a)')    "_symmetry_Int_Tables_number              '?'"
  else
   WRITE(CIF_unit, '(2a)')    "_symmetry_cell_setting                   ", TRIM(SPG%CrystalSys)
   WRITE(CIF_unit, '(3a)')    "_symmetry_space_group_name_H-M          '", TRIM(SPG%SPG_Symb),  "'"
   WRITE(CIF_unit, '(3a)')    "_symmetry_space_group_name_Hall         '", TRIM(SPG%Hall),      "'"
   WRITE(CIF_unit, '(a,I3)')  "_symmetry_Int_Tables_number              ", SPG%NumSpg
  endif
  WRITE(CIF_unit, '(a)')     "loop_"
  WRITE(CIF_unit, '(a)')     "_symmetry_equiv_pos_site_id"
  WRITE(CIF_unit, '(a)')     "_symmetry_equiv_pos_as_xyz"
  if(nb_symm_op == 0) nb_symm_op = SPG%Multip

  do i=1, nb_symm_op
   if(i< 10) then
    WRITE(CIF_unit, '(i1,1x,5a)')  i, TRIM(symm_op_string(1,i)), ", ", TRIM(symm_op_string(2,i)), ", ", TRIM(symm_op_string(3,i))
   elseif(i < 100) then
    WRITE(CIF_unit, '(i2,1x,5a)')  i, TRIM(symm_op_string(1,i)), ", ", TRIM(symm_op_string(2,i)), ", ", TRIM(symm_op_string(3,i))
   else
    WRITE(CIF_unit, '(i3,1x,5a)')  i, TRIM(symm_op_string(1,i)), ", ", TRIM(symm_op_string(2,i)), ", ", TRIM(symm_op_string(3,i))
   endif
  end do


      case ('CRYSTAL_SYSTEM')
  IF(unit_cell%crystal_system(1:1) /= '?') then
   WRITE(CIF_unit, '(a)')     ""
   WRITE(CIF_unit, '(3a)')    "_symmetry_cell_setting                  '", TRIM(unit_cell%crystal_system),"'"
  endif
  if(unit_cell%H_M(1:1) /= '?') then
   WRITE(CIF_unit, '(3a)')    "_symmetry_space_group_name_H-M          '", TRIM(unit_cell%H_M),"'"
  END if
  WRITE(CIF_unit, '(a)')     ""

      case ('CELL_PARAM')
  WRITE(CIF_unit, '(a)')     ""
  do i=1, 6
   WRITE(CIF_unit, '(a,F10.5)') CELL_param_CIF_string(i), unit_cell%param(i)
  end do
  WRITE(CIF_unit, '(a,F14.5)')  CELL_param_CIF_string(7), unit_cell%volume
  WRITE(CIF_unit ,'(2a)')     "_cell_measurement_temperature           ", TRIM(CIF_cell_measurement%temperature)


      case ('CELL_PARAM_ESD')
  IF(CIF_parameter%cell_length_a(1:1) /= '?') then
   WRITE(CIF_unit, '(2a)')   "_cell_length_a                         ", TRIM(CIF_parameter%cell_length_a)
  else
   call write_CIF_cell_parameter(1)
  endif

  IF(CIF_parameter%cell_length_b(1:1) /= '?') then
   WRITE(CIF_unit, '(2a)')   "_cell_length_b                         ", TRIM(CIF_parameter%cell_length_b)
  else
   call write_CIF_cell_parameter(2)
  endif

  IF(CIF_parameter%cell_length_c(1:1) /= '?') then
   WRITE(CIF_unit, '(2a)')   "_cell_length_c                         ", TRIM(CIF_parameter%cell_length_c)
  else
   call write_CIF_cell_parameter(3)
  endif

  IF(CIF_parameter%cell_angle_alpha(1:1) /= '?') then
   WRITE(CIF_unit, '(2a)')   "_cell_angle_alpha                      ", TRIM(CIF_parameter%cell_angle_alpha)
  else
   call write_CIF_cell_parameter(4)
  endif

  IF(CIF_parameter%cell_angle_beta(1:1) /= '?') then
   WRITE(CIF_unit, '(2a)')   "_cell_angle_beta                       ", TRIM(CIF_parameter%cell_angle_beta)
  else
   call write_CIF_cell_parameter(5)
  endif

  IF(CIF_parameter%cell_angle_gamma(1:1) /= '?') then
   WRITE(CIF_unit, '(2a)')   "_cell_angle_gamma                      ", TRIM(CIF_parameter%cell_angle_gamma)
  else
   call write_CIF_cell_parameter(6)
  endif


  IF(CIF_parameter%cell_volume(1:1) /= '?') then
   WRITE(CIF_unit, '(2a)')   "_cell_volume                           ", TRIM(CIF_parameter%cell_volume)
  else
   call write_CIF_cell_parameter(7)
  endif

  WRITE(CIF_unit ,'(2a)')      "_cell_formula_units_Z                 ", TRIM(CIF_parameter%cell_formula_units_Z)
  WRITE(CIF_unit ,'(2a)')      "_cell_measurement_temperature         ", TRIM(CIF_cell_measurement%temperature)


      case ('CIF_THMIN,THMAX')

  WRITE( CIF_unit, '(a,I6)')    "_cell_measurement_reflns_used          ", CIF_cell_measurement%reflns_used
  WRITE( CIF_unit, '(a,F6.2)')  "_cell_measurement_theta_min            ", CIF_cell_measurement%theta_min
  WRITE( CIF_unit, '(a,F6.2)')  "_cell_measurement_theta_max            ", CIF_cell_measurement%theta_max

      case ('AUTHOR')
  call write_CIF_section("AUTHORS DETAILS")
  write( CIF_unit, '(a)')  "_publ_contact_author               'ROISNEL, Thierry'                     #"
  write( CIF_unit, '(a)')  "_publ_contact_author_email          thierry.roisnel@univ-rennes1.fr       #"

  !write( CIF_unit, '(a)')  ""
  !write( CIF_unit, '(a)')  "#-------------------------------------------------------------------------#"
  !write( CIF_unit, '(a)')  "#                       AUTHORS DETAILS                                   #"
  !write( CIF_unit, '(a)')  "#"
  !write( CIF_unit, '(a)')  "_publ_contact_author               'ROISNEL, Thierry'                     #"
  !write( CIF_unit, '(a)')  "_publ_contact_author_email          thierry.roisnel@univ-rennes1.fr       #"
  !write( CIF_unit, '(a)')  "#                                                                         #"

  call write_CIF_section("SUBMISSION DETAILS")
  !write( CIF_unit, '(a)')  "#------------------ SUBMISSION DETAILS -----------------------------------#"
  !write( CIF_unit, '(a)')  ""
  write( CIF_unit, '(a)')  "# Name and address of author for correspondence"
  write( CIF_unit, '(a)')  ""
  write( CIF_unit, '(a)')  "_publ_contact_author_name      'Roisnel, Thierry"
  write( CIF_unit, '(a)')  "_publ_contact_author_address"
  write( CIF_unit, '(a)')  ";"
  write( CIF_unit, '(a)')  "Centre de Diffractom\'etrie X"
  write( CIF_unit, '(a)')  "Institut des Sciences Chimiques de Rennes"
  write( CIF_unit, '(a)')  "UMR6226 CNRS - Universit\'e de Rennes 1"
  write( CIF_unit, '(a)')  "B\^at. 10B, Campus de Beaulieu "
  write( CIF_unit, '(a)')  "Avenue du G\'en\'eral Leclerc"
  write( CIF_unit, '(a)')  "35042 Rennes, France "
  write( CIF_unit, '(a)')  ";"
  write( CIF_unit, '(a)')  "_publ_contact_author_email       thierry.roisnel@univ-rennes1.fr"
  write( CIF_unit, '(a)')  "_publ_contact_author_phone         '33(2)23235902'"
  write( CIF_unit, '(a)')  "_publ_contact_author_fax           '33(2)99383487'"
  write( CIF_unit, '(a)')  ""

      case ('ATOMS_HEADER')
  call write_CIF_section("ATOMIC TYPES, COORDINATES AND THERMAL PARAMETERS")

      case('ATOMS_HEADER_LOOP')
  write( CIF_unit, '(a)') "loop_"
  write( CIF_unit, '(a)') "   _atom_site_label"
  write( CIF_unit, '(a)') "   _atom_site_type_symbol"
  write( CIF_unit, '(a)') "   _atom_site_fract_x"
  write( CIF_unit, '(a)') "   _atom_site_fract_y"
  write( CIF_unit, '(a)') "   _atom_site_fract_z"
  write( CIF_unit, '(a)') "   _atom_site_U_iso_or_equiv"
  write( CIF_unit, '(a)') "   _atom_site_adp_type"
  write( CIF_unit, '(a)') "   _atom_site_occupancy"

      case ('ATOM')
  do i=1, nb_atom
   if(atom_adp_type(i)(1:5) == 'aniso') then
    adp_type = ' Uani '
   else
    adp_type = ' Uiso '
   end if
   WRITE( CIF_unit, '(4a,4F9.5,a,F8.5)')  atom_label(i)(1:4), ' ', atom_typ(i)(1:4), ' ', &
                                          atom_coord(1,i), atom_coord(2,i), atom_coord(3,i), atom_Ueq(i), &
                                          adp_type, atom_occ_perc(i)
   !WRITE( CIF_unit, '(4a,4F9.5,a,F8.5)')  atom_label(i)(1:4), ' ', atom_typ(i)(1:4), ' ', &
   !                                       atom_coord(1,i), atom_coord(2,i), atom_coord(3,i), atom_Ueq(i), &
   !                                       ' Uiso ', atom_occ_perc(i)
   !WRITE( CIF_unit, '(4a,4F9.5,a,F8.5)')  atom_label(i)(1:4), ' ', atom_typ(i)(1:4), ' ', &
   !                                       atom_coord(1,i), atom_coord(2,i), atom_coord(3,i), atom_Ueq(i), &
   !                                       ' Uiso ', atom_occ(i)
   END do

      case('ATOMS_ANISO_HEADER')
  write(CIF_unit, '(a)') " "
  write(CIF_unit, '(a)') "loop_"
  write(CIF_unit, '(a)') "     _atom_site_aniso_label"
  write(CIF_unit, '(a)') "     _atom_site_aniso_U_11"
  write(CIF_unit, '(a)') "     _atom_site_aniso_U_22"
  write(CIF_unit, '(a)') "     _atom_site_aniso_U_33"
  write(CIF_unit, '(a)') "     _atom_site_aniso_U_23"
  write(CIF_unit, '(a)') "     _atom_site_aniso_U_13"
  write(CIF_unit, '(a)') "     _atom_site_aniso_U_12"

     case('ATOMS_ANISO')
  do i=1, nb_atom
   write(CIF_unit, '(a4,6F10.5)') atom_label(i), atom_adp_aniso(1:3,i), atom_adp_aniso(6,i), atom_adp_aniso(5,i), &
                                  atom_adp_aniso(4,i)
  end do

      case('ATOMS_PYMOL')
  write( CIF_unit, '(a)') "loop_"
  write( CIF_unit, '(a)') "_atom_site_label"
  write( CIF_unit, '(a)') "_atom_site_type_symbol"
  write( CIF_unit, '(a)') "_atom_site_fract_x"
  write( CIF_unit, '(a)') "_atom_site_fract_y"
  write( CIF_unit, '(a)') "_atom_site_fract_z"
  do i=1, nb_atom
   WRITE( CIF_unit, '(4a,3F9.5)')  atom_label(i)(1:4), ' ', atom_typ(i)(1:4), ' ', &
                                          atom_coord(1,i), atom_coord(2,i), atom_coord(3,i)
  END do


      case ('P4P')
  open (CIF_unit, FILE='P4P.cif')
  do i=1, CIF_lines_nb
   WRITE(CIF_unit, '(a)') trim(CIF_title_line(i))
  end do
  WRITE(CIF_unit, '(a)') ""
  WRITE(CIF_unit, '(a)') "data_P4P"
  WRITE(CIF_unit, '(a)') ""

      case ('IMPORT')
  call create_import_cif

      case ('MATRIX_SAINT')
  call write_CIF_section("ORIENTATION MATRIX")
  WRITE(CIF_unit, '(a)')  "_diffrn_orient_matrix_type   'by SAINT Bruker-AXS'"
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_11  ", UB_matrix(1,1)
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_21  ", UB_matrix(2,1)
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_31  ", UB_matrix(3,1)
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_12  ", UB_matrix(1,2)
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_22  ", UB_matrix(2,2)
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_32  ", UB_matrix(3,2)
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_13  ", UB_matrix(1,3)
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_23  ", UB_matrix(2,3)
  WRITE(CIF_unit, '(a,F10.7)') "_diffrn_orient_matrix_UB_33  ", UB_matrix(3,3)
  WRITE(CIF_unit, '(a)')  ""


      case ('DIST')
  write( CIF_unit, '(a)') ''
  write( CIF_unit, '(a)') 'loop_'
  write( CIF_unit, '(a)') '_geom_bond_atom_site_label_1'
  write( CIF_unit, '(a)') '_geom_bond_atom_site_label_2'
  write( CIF_unit, '(a)') '_geom_bond_distance'
  write( CIF_unit, '(a)') '_geom_bond_site_symmetry_2'
  write( CIF_unit, '(a)') '_geom_bond_publ_flag'

      case('DIST_values')
  do i=1, CIF_dist%n_text
   write( CIF_unit, '(a)') trim(CIF_dist%text(i))
  end do

      case ('ANG')
  write( CIF_unit, '(a)') ''
  write( CIF_unit, '(a)') 'loop_'
  write( CIF_unit, '(a)') '_geom_angle_atom_site_label_1'
  write( CIF_unit, '(a)') '_geom_angle_atom_site_label_2'
  write( CIF_unit, '(a)') '_geom_angle_atom_site_label_3'
  write( CIF_unit, '(a)') '_geom_angle'
  write( CIF_unit, '(a)') '_geom_angle_site_symmetry_1'
  write( CIF_unit, '(a)') '_geom_angle_site_symmetry_3'
  write( CIF_unit, '(a)') '_geom_angle_publ_flag'


      case ('REFINEMENT_INFO')
  call write_CIF_section("REFINEMENT INFORMATION")
  write( CIF_unit, '(a)') '_refine_special_details'
  write( CIF_unit, '(a)') ';'
  write( CIF_unit, '(a)') 'Refinement of F^2^ against ALL reflections.  The weighted R-factor wR and'
  write( CIF_unit, '(a)') 'goodness of fit S are based on F^2^, conventional R-factors R are based'
  write( CIF_unit, '(a)') 'on F, with F set to zero for negative F^2^. The threshold expression of'
  write( CIF_unit, '(a)') 'F^2^ > 2sigma(F^2^) is used only for calculating R-factors(gt) etc. and is'
  write( CIF_unit, '(a)') 'not relevant to the choice of reflections for refinement.  R-factors based'
  write( CIF_unit, '(a)') 'on F^2^ are statistically about twice as large as those based on F, and R-'
  write( CIF_unit, '(a)') 'factors based on ALL data will be even larger.'
  write( CIF_unit, '(a)') ';'
  write( CIF_unit, '(a)') '_refine_ls_structure_factor_coef        Fsqd'
  write( CIF_unit, '(a)') '_refine_ls_matrix_type                  full'
  write( CIF_unit, '(a)') '_refine_ls_weighting_scheme             calc'

      case ('MOLECULAR_GEOMETRY')
  call write_CIF_section("MOLECULAR GEOMETRY")
  write( CIF_unit, '(a)') '_geom_special_details'
  write( CIF_unit, '(a)') ';'
  write( CIF_unit, '(a)') ' All esds (except the esd in the dihedral angle between two l.s. planes)'
  write( CIF_unit, '(a)') ' are estimated using the full covariance matrix.  The cell esds are taken'
  write( CIF_unit, '(a)') ' into account individually in the estimation of esds in distances, angles'
  write( CIF_unit, '(a)') ' and torsion angles; correlations between esds in cell parameters are only'
  write( CIF_unit, '(a)') ' used when they are defined by crystal symmetry.  An approximate (isotropic)'
  write( CIF_unit, '(a)') ' treatment of cell esds is used for estimating esds involving l.s. planes.'
  write( CIF_unit, '(a)') ';'



      case default
  WRITE(CIF_unit, '(a)')  ''


  end select

 return

end subroutine write_CIF_file

!------------------------------------------------------------------------------
subroutine create_CIF_P4P_import(input_string)
! creation du fichier "import.CIF" apres lecture fichier .P4P et .HKL
 USE cryscalc_module, ONLY : known_cell_esd, wavelength, CIF_unit, crystal, UB_matrix,                 &
                             DEVICE, keyword_SIZE, keyword_WAVE, molecule, debug_proc
 USE CIF_module,      ONLY : CIF_CELL_measurement
 USE IO_module
 implicit none
  CHARACTER (LEN=*),   INTENT(IN)      :: input_string

  if(debug_proc%level_2)  call write_debug_proc_level(2, "create_CIF_P4P_import ("//trim(input_string)//")")



  select case (input_string)
      case ('P4P')
   call write_CIF_file('P4P')

      case ('IMPORT')
   call write_CIF_file('IMPORT')
   call write_CIF_file('MATRIX_SAINT')
   call write_CIF_file('CHEMICAL_INFORMATION')
   call write_CIF_file('CHEMICAL_2')

      case default
  end select

  call write_CIF_file('UNIT_CELL_INFO')
  call write_CIF_file('CRYSTAL_SYSTEM')


  IF(known_cell_ESD) then
   call write_CIF_file('CELL_PARAM_ESD')
  else
   call write_CIF_file('CELL_PARAM')
  endif

  call write_CIF_file('CIF_THMIN,THMAX')

  !IF(keyword_WAVE) call write_CIF_file('WAVE')
  IF(keyword_WAVE) call write_CIF_file('CELL_WAVE')

  call write_CIF_file('CRYSTAL_INFORMATION')
  call write_CIF_file('ABSORPTION')
  call write_CIF_file('SADABS')
  !call Get_SADABS_ratio()
  call write_CIF_file('TMIN')
  call write_CIF_file('TMAX')
  if(DEVICE%APEX .or. DEVICE%X2S) then
   call write_CIF_file('APEX')
  elseif(DEVICE%D8V) then
   call write_CIF_file('D8V')
  end if
 close (unit=CIF_unit)

 if(len_trim(input_string) == 3) then
  if(input_string(1:3) == 'P4P') then
   call write_info('')
   call write_info(' > P4P.cif file has been created.')
   call write_info('')
  end if
 end if

 return
end subroutine create_CIF_P4P_import

!------------------------------------------------------------
subroutine verif_CIF_character(input_string)
 USE macros_module,   ONLY : replace_car
 USE accents_module
 USE cryscalc_module, ONLY : debug_proc
 implicit none
  CHARACTER (LEN=256), intent(INOUT) :: input_string
  INTEGER                            :: i, j

  if(debug_proc%level_3)  call write_debug_proc_level(3, "verif_CIF_character ("//trim(input_string)//")")

  do i=1, LEN_TRIM(input_string)

   do j=1, nb_char
    if(input_string(i:i) == accents(1, j)(1:1)) then
     input_string = replace_car(input_string, accents(1,j)(1:1), trim(accents(2,j)))!
     exit
    endif
   end do

  end do

 RETURN
end subroutine verif_CIF_character

!---------------------------------------------------------------
subroutine write_CIF_cell_parameter(i)
 USE cryscalc_module, ONLY : unit_cell , CIF_unit, P4P_file_name, debug_proc
 implicit none
 INTEGER, INTENT(IN)                               :: i
 CHARACTER (LEN=16)                                :: ESD_string
 CHARACTER (LEN=36), DIMENSION(7)                  :: CELL_param_CIF_string
 REAL                                              :: uc_param, uc_param_esd
 CHARACTER (LEN=16)                                :: fmt_1, fmt_2

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_CIF_cell_parameter")

  cell_param_CIF_string(1) =  "_cell_length_a                      "
  cell_param_CIF_string(2) =  "_cell_length_b                      "
  cell_param_CIF_string(3) =  "_cell_length_c                      "
  cell_param_CIF_string(4) =  "_cell_angle_alpha                   "
  cell_param_CIF_string(5) =  "_cell_angle_beta                    "
  cell_param_CIF_string(6) =  "_cell_angle_gamma                   "
  cell_param_CIF_string(7) =  "_cell_volume                        "

  ESD_string = "?"

  if(i /=7) then
   uc_param     = unit_cell%param(i)
   uc_param_ESD = unit_cell%param_ESD(i)
  else
   uc_param     = unit_cell%volume
   uc_param_ESD = unit_cell%volume_ESD
  endif

 ! fmt_1 : format ecriture avec esd
 ! fmt_2 : format ecriture sans esd
 call create_ESD_string(i, ESD_string, fmt_1, fmt_2)


 !IF(i<=3) then                  ! parametres de maille a, b, c
 ! if(P4P_file_name /= '')then   ! lecture d'un fichier.P4P: les CELL et CELLSD sont en F10.4
 !  m_expo = 4
 ! else
 !  m_expo = 5
 ! endif

 !ELSEIF(i<=6) then          ! parametres de maille alpha, beta, gamma
 ! m_expo = 4

 !ELSEIF(i==7) then          ! volume
 ! m_expo = 3

 !endif

 !write(fmt_1, '(a,i1,a)') '(a,F10.', m_expo, ',3a)'
 !write(fmt_2, '(a,i1,a)') '(a,F10.', m_expo, ')'
 !if(uc_param_ESD < 10.) then
 ! write(fmt_3, '(a,i1,a,i1,a)') '(F', m_expo+2, '.', m_expo, ')'
 !elseif(uc_param_ESD < 100.) then
 ! write(fmt_3, '(a,i1,a,i1,a)') '(F', m_expo+3, '.', m_expo, ')'
 !endif


 !IF(INT(uc_param_ESD  * 10.**m_expo +0.5) /=0) then
!
 ! if(uc_param_ESD < 1.) then
 !  write(ESD_string, *) int(uc_param_ESD*10**m_expo + 0.5)
 ! elseif(uc_param_ESD < 10.) then
 !  write(ESD_string, fmt_3) uc_param_ESD
 ! elseif(uc_param_ESD < 100.) then
 !  write(ESD_string, fmt_3) uc_param_ESD
 ! endif
 ! ESD_string = adjustl(ESD_string)
 ! write( CIF_unit, fmt_1) CELL_param_CIF_string(i), uc_param, '(', trim(ESD_string), ')'

 !else
 ! WRITE( CIF_unit, fmt_2) CELL_param_CIF_string(i), uc_param
 !endif

  if(ESD_string(1:1) /= "?") then
   write( CIF_unit, fmt_1) CELL_param_CIF_string(i), uc_param, '(', trim(ESD_string), ')'
  else
   WRITE( CIF_unit, fmt_2) CELL_param_CIF_string(i), uc_param
  end if

 return
end subroutine write_CIF_cell_parameter

!-------------------------------------------------------------------------
subroutine Create_ESD_string(i, ESD_string, fmt_1, fmt_2)
 use Cryscalc_module, only : unit_cell, input_CELL, P4P_file_name, keyword_read_ins
 implicit none
  integer, intent(in)               :: i
  character (len=16), intent(inout) :: ESD_string
  character (len=16), intent(inout) :: fmt_1, fmt_2
  REAL                              :: uc_param, uc_param_esd
  integer                           :: m_expo
  CHARACTER (LEN=16)                :: fmt_3


  ! fmt_1 : format ecriture avec esd
  ! fmt_2 : format ecriture sans esd
  ! fmt_3 : format ecriture esd

  ! i=1,6 : cell_param
  ! i=7   : cell_volume

  if(i /=7) then
   uc_param     = unit_cell%param(i)
   uc_param_ESD = unit_cell%param_ESD(i)
  else
   uc_param     = unit_cell%volume
   uc_param_ESD = unit_cell%volume_ESD
  endif

  IF(i<=3) then                  ! parametres de maille a, b, c
   !if(P4P_file_name /= '')then   ! lecture d'un fichier.P4P: les CELL et CELLSD sont en F10.4
    m_expo = 4
   !else
   ! m_expo = 5
   !endif
  ELSEIF(i<=6) then              ! parametres de maille alpha, beta, gamma
   m_expo = 4
   if(input_CELL%INS .or. keyword_read_ins) m_expo = 3

  ELSEIF(i==7) then          ! volume
   m_expo = 3
  endif

  write(fmt_1, '(a,i1,a)') '(a,F10.', m_expo, ',3a)'
  !write(fmt_2, '(a,i1,a)') '(a,F10.', m_expo, ')'
  write(fmt_2, '(a,i1,a)') '(a,F10.', m_expo, ',a)'
  if(uc_param_ESD < 10.) then
   write(fmt_3, '(a,i1,a,i1,a)') '(F', m_expo+2, '.', m_expo, ')'
  elseif(uc_param_ESD < 100.) then
   write(fmt_3, '(a,i1,a,i1,a)') '(F', m_expo+3, '.', m_expo, ')'
  endif

  IF(INT(uc_param_ESD  * 10.**m_expo +0.5) /=0) then
   if(uc_param_ESD < 1.) then
    write(ESD_string, *) int(uc_param_ESD*10**m_expo + 0.5)
   elseif(uc_param_ESD < 10.) then
    write(ESD_string, fmt_3) uc_param_ESD
   elseif(uc_param_ESD < 100.) then
    write(ESD_string, fmt_3) uc_param_ESD
   endif
   ESD_string = adjustl(ESD_string)
  endif



 return
end subroutine Create_ESD_string
!--------------------------------------------------------------------

subroutine write_CIF(input_unit, input_string)
 use cryscalc_module, only : CIF_unit, CIF_archive_unit

 implicit none
  integer,           intent(in)        :: input_unit
  character (len=*), intent(in)        :: input_string

    WRITE(input_unit,         '(a)') trim(input_string)

 return
end subroutine write_CIF
!--------------------------------------------------------------------
subroutine write_CIF_author(input_unit)
 use cryscalc_module, only : AUTHOR, keyword_modif_archive, debug_proc
 use Accents_module,  ONLY : def_accents

 implicit none
  integer, intent(in)          :: input_unit
  character (len=256)          :: CIF_string

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_CIF_author")


   !call write_CIF(input_unit, "#----------------------------------------------------------------------------#")
   call write_CIF(input_unit, "")
   !if(keyword_modif_archive) then
    call write_CIF(input_unit, "data_global")
   !else
   ! call write_CIF(input_unit, "data_cryscalc")
   !endif
   call write_CIF_section("AUTHORS DETAILS")
   call write_CIF(input_unit, "#   >> Person making the deposition :")
   call write_CIF(input_unit, "")
   call write_CIF(input_unit, "")


   call def_accents
   call verif_CIF_character(AUTHOR%NAME)
   call verif_CIF_character(AUTHOR%first_name)
   call verif_CIF_character(AUTHOR%address)
   call verif_CIF_character(AUTHOR%email)
   call verif_CIF_character(AUTHOR%web)

   WRITE(CIF_string, '(5a)') "_publ_contact_author_name          '", TRIM(AUTHOR%name),", ", TRIM(AUTHOR%first_name), "'"
    call write_CIF(input_unit, trim(CIF_string))
   WRITE(CIF_string, '(2a)') "_publ_contact_author_email          ", trim(AUTHOR%email)
    call write_CIF(input_unit, trim(CIF_string))
   WRITE(CIF_string, '(a)')  "_publ_contact_author_address        "
    call write_CIF(input_unit, trim(CIF_string))
    call write_CIF_address(input_unit, AUTHOR%address)

   call write_CIF_section("PUBLICATION DETAILS")
   call write_CIF(input_unit, "#     Provide these details if the structure has been published,")
   call write_CIF(input_unit, "#     accepted or submitted for publication")
   call write_CIF(input_unit, "#")
   call write_CIF(input_unit, "#     The CCDC journal deposition number, eg. 182/357,")
   call write_CIF(input_unit, "#     should be included only if it has been assigned by the journal")
   call write_CIF(input_unit, "")
   call write_CIF(input_unit, "loop_")
   call write_CIF(input_unit, "_publ_author_name")
   call write_CIF(input_unit, "_publ_author_address")
   call write_CIF(input_unit, "")


   WRITE(CIF_string, '(5a)') "     '",TRIM(AUTHOR%name),", ", TRIM(AUTHOR%first_name), "'"
    call write_CIF(input_unit, trim(CIF_string))
    call write_CIF_address(input_unit, AUTHOR%address)
    call write_CIF(input_unit, "")
    call write_CIF(input_unit, "_journal_name_full                 'private communication'")
    call write_CIF(input_unit, "_journal_volume                    ?")
    call write_CIF(input_unit, "_journal_page_first                ?")
    call write_CIF(input_unit, "_journal_page_last                 ?")
    call write_CIF(input_unit, "_journal_year                      ?")
    call write_CIF(input_unit, "_ccdc_journal_depnumber            ?")
    call write_CIF(input_unit, "")
    !call write_CIF(input_unit, "#----------------------------------------------------------------------------#")


 return
end subroutine write_CIF_author

!---------------------------------------------------------------------
subroutine write_CIF_text(input_unit)
 use cryscalc_module, only : SQUEEZE, structure_refinement, debug_proc
 use CIF_module,      only : CIF_parameter, CIF_sep_line_full
 implicit none  
  integer,   intent(in) :: input_unit
  logical               :: SHELXL_2014
  character (len=512)   :: CIF_string
  integer               :: long


  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_CIF_text")

   call write_CIF_section("TEXT")
   call write_CIF(input_unit, "_publ_section_abstract")
   call write_CIF(input_unit, ";")
   call transf_moiety_string("CIF", CIF_parameter%formula_moiety, CIF_string)
   call write_CIF(input_unit, "The molecule of the title compound "//trim(CIF_string))
   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, "")
   call write_CIF(input_unit, "_publ_section_comment")
   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, "")
   call write_CIF(input_unit, "_publ_section_related_literature")
   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, "")
   call write_CIF(input_unit, "_publ_section_exptl_prep")
   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, "")
   call write_CIF(input_unit, "_publ_section_exptl_refinement")
   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, "All non-H atoms were refined with anisotropic atomic displacement parameters.")
   if(CIF_parameter%H_treatment(1:5) == 'mixed') then
    call write_CIF(input_unit, "Except XXX linked atoms there were introduced in the structural model")
    call write_CIF(input_unit, "through Fourier difference maps analysis,")
   endif
   call write_CIF(input_unit, "H atoms were finally included in their calculated positions and treated as riding")
   call write_CIF(input_unit, "on their parent atom.")

   SHELXL_2014 = .false.
   long = len_trim(structure_refinement%name)
   if(long == 11) then
   if(structure_refinement%name(1:11) == 'SHELXL-2014' .or. &
      structure_refinement%name(1:11) == 'SHELXL-2016') SHELXL_2014 = .true.
   endif
   
   if(SQUEEZE%procedure) then
    if(SHELXL_2014) then
     call write_CIF(input_unit, "The contribution of the disordered solvents to the structure factors was ")
     call write_CIF(input_unit, "calculated by the PLATON/SQUEEZE procedure and then taken into account ")
     call write_CIF(input_unit, "in the final SHELXL-2014 least-square refinement.")
    else
     call write_CIF(input_unit, "The contribution of the disordered solvents to the calculated structure factors")
     call write_CIF(input_unit, "was estimated following the BYPASS algorithm, implemented as the SQUEEZE option")
     call write_CIF(input_unit, "in PLATON. A new data set, free of solvent contribution, was then used in the")
     call write_CIF(input_unit, "final refinement.")
	end if 
   endif

   call write_CIF(input_unit, ";")
   call write_CIF(input_unit, "")
   call write_CIF(input_unit, "#")
   WRITE(input_unit, '(a)')  trim(CIF_sep_line_full)
   call write_CIF(input_unit, "")





 return
end subroutine write_CIF_text


!---------------------------------------------------------------------
subroutine write_CIF_address(input_unit, address)
 USE cryscalc_module, only : debug_proc
 implicit none
  integer,           intent(in)         :: input_unit
  character (len=*), intent(in)         :: address
  character (len=256)                   :: CIF_string, address_part
  integer                               :: i, i1

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_CIF_address")


  ! ----------- new : sept. 2011 ---------------------
  write(address_part, '(a)') trim(address)

  call write_CIF(input_unit, ";")
  do
   i1 = index(trim(address_part), ",")
   if(i1 /=0) then
    write(CIF_string, '(2a)') "     ", trim(adjustl(address_part(1:i1-1)))
    call write_CIF(input_unit, trim(CIF_string))
    address_part = address_part(i1+1:)
   else
    call write_CIF(input_unit, "     "//trim(adjustl(address_part)))
    exit
   endif

  end do
  call write_CIF(input_unit, ";")

  return
  ! ---------------------------------------------------------------

   IF(LEN_TRIM(address) > 70) then
    do i=1, INT(LEN_TRIM(address)/70) + 1
     IF(i==1) then
      WRITE(CIF_string, '(2a)') ";     ", address(1:70)
       call write_CIF(input_unit, trim(CIF_string))
     else
      WRITE(CIF_string, '(2a)') "      ", trim(address(1+70*(i-1):1+70*(i-1)+69))
       call write_CIF(input_unit, trim(CIF_string))
     endif
    end do
    WRITE(CIF_string, '(a)') ";"
     call write_CIF(input_unit, trim(CIF_string))
   else
    WRITE(CIF_string, '(2a)') "      ",TRIM(address)
     call write_CIF(input_unit, trim(CIF_string))
   endif

 return
end subroutine write_CIF_address

!----------------------------------------------------------------------------------------
subroutine create_CIF_PYMOL_file(input_file, extension)
 use cryscalc_module, only   : CIF_PYMOL_file_name, debug_proc
 implicit none
 CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
 integer                              :: i

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_CIF_PYMOL_file")

 i= INDEX(input_file, '.')
 WRITE(CIF_PYMOL_file_name, '(a)') input_file(1:i-1)//'_'//TRIM(extension)//'_pml.CIF'
 call Create_CIF_file_for_pymol

 return
end subroutine create_CIF_PYMOL_file

!----------------------------------------------------------------------------------------
subroutine Create_CIF_file_for_pymol
 use cryscalc_module, only : CIF_pymol_file_name, CIF_unit, debug_proc
 use CIF_module,      only : CIF_parameter
 use IO_module,       only : write_info

 implicit none

 if(debug_proc%level_2)  call write_debug_proc_level(2, "create_CIF_file_for_PYMOL")


 OPEN(UNIT=CIF_unit, FILE=TRIM(CIF_PYMOL_file_name), ACTION="write")
  write(UNIT=CIF_unit, fmt='(a)') '#######################################################'

  write(unit=CIF_unit, fmt='(2a)') 'data_', trim(CIF_parameter%sample_ID)
  call write_CIF_file("SPACE_GROUP_PYMOL")
  call write_CIF_file("CELL_PARAM_ESD")
  call write_CIF_file("ATOMS_PYMOL")

  write(UNIT=CIF_unit, fmt='(a)') ''
  write(UNIT=CIF_unit, fmt='(a)') '#END'
  write(UNIT=CIF_unit, fmt='(a)') '#######################################################'

 CLOSE(UNIT=CIF_unit)

 call write_info('')
 call write_info('   >> '//trim(CIF_PYMOL_file_name)//' CIF file for Pymol has been created.')
 call write_info('')

 return
end subroutine Create_CIF_file_for_pymol


!-----------------------------------------------------------------------------------------------

subroutine Get_sample_ID(CIF_unit)
 use cryscalc_module, only : debug_proc
 use CIF_module,      only : CIF_parameter
  implicit none
  INTEGER, INTENT(IN)      :: CIF_unit
  integer                  :: i_error
  character (len=256)      :: read_line

  if(debug_proc%level_2)  call write_debug_proc_level(2, "get_sample_ID")


 REWIND(UNIT=CIF_unit)

 do
   read(CIF_unit, '(a)', iostat = i_error) read_line

   if(i_error /=0) exit
   if(index(read_line, 'data_') /=0 ) then
    read(read_line(6:), '(a)') CIF_parameter%sample_ID
    exit
   endif
  enddo

 return
end subroutine Get_sample_ID

!-------------------------------------------------------------------------------------------------
subroutine create_Struct_CIF
 use Cryscalc_module,   only : CIF_unit, nb_atoms_type, SFAC, keyword_SIZE, crystal, SPG, Z_unit, molecule
 use CIF_module,        ONLY : CIF_parameter, CIF_parameter_DEVICE, CIF_cell_measurement
 USE HKL_module,        ONLY : HKL_file
 USE IO_module,         ONLY : write_info

 implicit none
 integer              :: i, i1
 CHARACTER (LEN=10)   :: date, time


 i1 = index(HKL_file%name, '.')

 open(unit=CIF_unit, file="struct.cif")
  write(CIF_unit, '(a)') 'data_'//HKL_file%name(1:i1-1)
  write(CIF_unit, '(a)') "_audit_creation_method                 'CRYSCALC routine'"

  call date_and_time(date, time)
  WRITE(CIF_unit, '(a,2x,13a)')"_audit_creation_date" , "'",date(7:8),'-',date(5:6),'-',date(1:4), '  at ', &
                                 time(1:2),':',time(3:4),':',time(5:6),"'"
  WRITE(CIF_unit, '(a)')  ""
  WRITE(CIF_unit, '(a)')  "loop_"
  WRITE(CIF_unit, '(a)')  "  _atom_type_symbol"
  WRITE(CIF_unit, '(a)')  "  _atom_type_number_in_cell"
  do i = 1, nb_atoms_type
   write(CIF_unit, '(a,1x,I5)')  SFAC%type(i), INT(SFAC%number(i))
  end do
  WRITE(CIF_unit, '(2a)')     "_cell_length_a                         ", TRIM(CIF_parameter%cell_length_a)
  WRITE(CIF_unit, '(2a)')     "_cell_length_b                         ", TRIM(CIF_parameter%cell_length_b)
  WRITE(CIF_unit, '(2a)')     "_cell_length_c                         ", TRIM(CIF_parameter%cell_length_c)
  WRITE(CIF_unit, '(2a)')     "_cell_angle_alpha                      ", TRIM(CIF_parameter%cell_angle_alpha)
  WRITE(CIF_unit, '(2a)')     "_cell_angle_beta                       ", TRIM(CIF_parameter%cell_angle_beta)
  WRITE(CIF_unit, '(2a)')     "_cell_angle_gamma                      ", TRIM(CIF_parameter%cell_angle_gamma)
  WRITE(CIF_unit, '(a, I4)')  "_cell_formula_units_Z                  ", int(Z_unit)
  WRITE(CIF_unit ,'(2a)')     "_cell_measurement_temperature          ", TRIM(CIF_cell_measurement%temperature)
  WRITE(CIF_unit, '(2a)')     "_cell_measurement_wavelength           ", trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength)
  WRITE(CIF_unit, '(a,I6)')   "_cell_measurement_reflns_used          ", CIF_cell_measurement%reflns_used
  WRITE(CIF_unit, '(a,F6.2)') "_cell_measurement_theta_min            ", CIF_cell_measurement%theta_min
  WRITE(CIF_unit, '(a,F6.2)') "_cell_measurement_theta_max            ", CIF_cell_measurement%theta_max
  WRITE(CIF_unit ,'(2a)')     "_diffrn_ambient_temperature            ", TRIM(CIF_cell_measurement%temperature)
  WRITE(CIF_unit, '(2a)')     "_diffrn_radiation_probe                ", trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator)
  WRITE(CIF_unit, '(2a)')     "_diffrn_radiation_source               ", trim(CIF_parameter_DEVICE%diffrn_radiation_source)
  WRITE(CIF_unit, '(2a)')     "_diffrn_radiation_monochromator        ", trim(CIF_parameter_DEVICE%diffrn_radiation_probe)

  WRITE(CIF_unit, '(2a)')           "_exptl_crystal_description             ", TRIM(crystal%morph)
  WRITE(CIF_unit, '(2a)')           "_exptl_crystal_colour                  ", TRIM(crystal%color)

  IF(keyword_SIZE) then
   WRITE(CIF_unit, '(a,F6.3)')      "_exptl_crystal_size_max                ", crystal%size_max
   WRITE(CIF_unit, '(a,F6.3)')      "_exptl_crystal_size_mid                ", crystal%size_mid
   WRITE(CIF_unit, '(a,F6.3)')      "_exptl_crystal_size_min                ", crystal%size_min
  else
   WRITE(CIF_unit, '(a)')           "_exptl_crystal_size_max                ?"
   WRITE(CIF_unit, '(a)')           "_exptl_crystal_size_mid                ?"
   WRITE(CIF_unit, '(a)')           "_exptl_crystal_size_min                ?"
  endif
  WRITE(CIF_unit, '(a, F9.3)')      "_exptl_crystal_density_method          'not measured'"

  WRITE(CIF_unit, '(3a)')           "_chemical_formula_moiety               '", adjustl(trim(molecule%formula)), "'"
  WRITE(CIF_unit, '(3a)')           "_chemical_formula_sum                  '", adjustl(trim(molecule%formula)), "'"
  WRITE(CIF_unit, '(a, F9.3)')      "_chemical_compound_source              'synthesis as described'"

  WRITE(CIF_unit, '(2a)')           "_symmetry_cell_setting                 ", TRIM(SPG%CrystalSys)
  WRITE(CIF_unit, '(3a)')           "_symmetry_space_group_name_H-M        '", TRIM(SPG%SPG_Symb),  "'"
  WRITE(CIF_unit, '(3a)')           "_symmetry_space_group_name_Hall       '", TRIM(SPG%Hall),      "'"
  WRITE(CIF_unit, '(a,I3)')         "_symmetry_Int_Tables_number            ", SPG%NumSpg

  WRITE(CIF_unit, '(2a)') "_computing_data_collection             ", trim(CIF_parameter_DEVICE%computing_data_reduction)
  WRITE(CIF_unit, '(2a)') "_computing_cell_refinement             ", trim(CIF_parameter_DEVICE%computing_cell_refinement)
  WRITE(CIF_unit, '(2a)') "_computing_data_reduction              ", trim(CIF_parameter_DEVICE%computing_data_reduction)
  WRITE(CIF_unit, '(2a)') "_computing_structure_solution          ", trim(CIF_parameter_DEVICE%computing_structure_solution)
  WRITE(CIF_unit, '(2a)') "_computing_structure_refinement        ", trim(CIF_parameter_DEVICE%computing_structure_refinement)
  WRITE(CIF_unit, '(2a)') "_computing_molecular_graphics          ", trim(CIF_parameter_DEVICE%computing_molecular_graphics)
  WRITE(CIF_unit, '(3a)') "_computing_publication_material        '",trim(CIF_parameter_DEVICE%computing_publication_material_2),"'"

 close(unit=CIF_unit)

  call write_info('')
  call write_info('   >> struct.cif file has been created')
  call write_info('')

 return
end subroutine create_Struct_CIF

!--------------------------------------------------------------------------------------------------------
 subroutine write_cif_line(import_string)
  use cryscalc_module, only : cif_unit
  use CIF_module,      only : CIF_line_width
  implicit none
  character (len=*), intent(in)     :: import_string
  integer                           :: long
  character (len=32)                :: fmt_


    long = len_trim(import_string)
    write(fmt_, '(a,i2,a)') '(a,1x,a,', CIF_line_width - long -3, 'x,a)'
	!write(fmt_, '(a,i2,a)') '(a,a,', CIF_line_width - long -2, 'x,a)'
    WRITE(CIF_unit, fmt = trim(fmt_)) "#", trim(import_string), "#"

  return

 end subroutine write_cif_line

!--------------------------------------------------------------------------------------------------------
subroutine center_CIF_line(import_string)
  use CIF_module,      only : CIF_line_width
  use cryscalc_module, only : cif_unit
  implicit none
  character (len=*), intent(in)     :: import_string
  character (len=80)                :: new_string
  character (len=32)                :: fmt_
  integer                           :: i, long_1, long

  long = len_trim(import_string)
  !long_1 = CIF_line_width/2 - int(long/2)  
  long_1 = int((CIF_line_width - long)/2)  -1 
  

  write(fmt_, '(a,i2,a)') '(', long_1, 'a1,a)'
  write(new_string, trim(fmt_)) (' ', i=1, long_1), trim(import_string)
 !call write_cif_line(trim(new_string))

  long = len_trim(new_string) +2 
  write(fmt_, '(a,i2,a)') '(a,a,', CIF_line_width - long, 'x,a)'
  WRITE(CIF_unit, fmt = trim(fmt_)) "#", trim(new_string), "#"

  return
 end subroutine center_CIF_line


!----------------------------------------------------------------------------------------------------------
 subroutine write_CIF_section(section_title)
  use cryscalc_module, only       : CIF_unit
  use CIF_module,      only       : CIF_sep_line_full, CIF_sep_line_empty

  implicit none
  character (len=*), intent(in)   :: section_title


  WRITE(CIF_unit, '(a)')        ""
  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_full)
  CALL center_CIF_line(trim(section_title))
  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_full)
  WRITE(CIF_unit, '(a)')        " "

  return
 end subroutine write_CIF_section

!---------------------------------------------------------------------------------------------------------
 subroutine Create_import_cif
  USE Cryscalc_module, only : CIF_unit, AUTHOR, DEVICE, Crystal, Experiment, Proposer,                                          &
                              CRYSCALC, SAINT, SADABS, TWINABS_file, TWINABS_4, P4P_file_name, ABS_file_name, RAW_file_name,    &
                              HKL_file_name, EXP_scan

  USE CIF_module
  USE hkl_module
  USE IO_module
  USE macros_module, only   : test_file_exist, l_case, u_case, Get_current_folder_name, Get_Proposer_from_current_folder_name, &
                              replace_car


  implicit none
   character (len=80)                                :: import_string, string
   CHARACTER (LEN=10)                                :: date, time
   character (len=256)                               :: current_folder
   CHARACTER (LEN=256)                               :: SHELX_HKL_file
   LOGICAL                                           :: file_exist
   integer                                           :: i, i1, long

  open (CIF_unit, FILE='import.cif')
  WRITE(CIF_unit, '(a)')  ""
  call date_and_time(date, time)

  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_full)
  write(import_string, '(3a)') "CIF file created by CRYSCALC (",  trim(CRYSCALC%version), ")"
  call center_CIF_line(trim(import_string))
  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_full)
  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)

 
  long = len_trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
  write(import_string, '(2a)') "Instrument             : ", trim(CIF_parameter_DEVICE%diffrn_measurement_device_type(2:long-1))
  call write_cif_line(trim(import_string))


  write(import_string, '(2a)') "Laboratory             : ", trim(DEVICE%lab)
  call write_cif_line(trim(import_string))

  call Get_experiment_features
  if(len_trim(experiment%date) /=0) then
   write(import_string, '(2a)')   "Date of experiment     : ", trim(Experiment%date)
   call write_cif_line(trim(import_string))
  else
   write(import_string, '(2a)')   "Date of data reduction : ", trim(SAINT%date)
   call write_cif_line(trim(import_string))
  end if

  call Get_current_folder_name(current_folder)
  call Get_Proposer_from_current_folder_name(current_folder, proposer%name, proposer%team)
  Proposer%name = replace_car(Proposer%name, '_', '. ')
  long = len_trim(Proposer%name)
  if(long == 5) then
   if(Proposer%name(1:5) == 'guest') then
    Proposer%name = 'Unknown'
    proposer%team = 'Unknown'
   end if
  end if

  call test_ISCR(trim(proposer%team), Proposer%ISCR)


  if(proposer%name(1:1) /= '?' .and. proposer%team(1:1) /= '?') then
   if(Proposer%ISCR) then
    write(import_string, '(5a)')  "Proposer               : ", trim(proposer%name), ' [',trim(proposer%team),'/ISCR]'
   else
    write(import_string, '(5a)')  "Proposer               : ", trim(proposer%name), ' [',trim(proposer%team),']'
   end if
  call write_cif_line(trim(import_string))
  end if

  write(import_string, '(5a)') "Experimenter           : ", trim(AUTHOR%first_name), " ", u_case(AUTHOR%name(1:1)), &
                                                                                     l_case(trim(AUTHOR%name(2:)))
  call write_cif_line(trim(import_string))

  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
  write(import_string, '(2a)')    "Sample ID              : ", trim(CIF_parameter%sample_ID)
  call write_cif_line(trim(import_string))

  write(import_string, '(a,2(1x,F5.2))') "Sample mosaic          : ", crystal%mosaic(1:2)
  call write_cif_line(trim(import_string))

  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
  !WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_full)
  write(import_string, '(2a)')           "X-ray target           : ", trim(experiment%target)
  call write_cif_line(trim(import_string))
  write(import_string, '(2a)')           "X-ray voltage          : ", trim(experiment%voltage)
  call write_cif_line(trim(import_string))
  write(import_string, '(2a)')           "X-ray current          : ", trim(experiment%current)
  call write_cif_line(trim(import_string))
  write(import_string, '(2a)')           "Detector               : ", trim(experiment%detector)
  call write_cif_line(trim(import_string))
  write(import_string, '(2a)')           "Detector temperature   : ", trim(experiment%detector_TEMP)
  call write_cif_line(trim(import_string))

  write(import_string, '(2a)')           "Experiment temperature : ", trim(experiment%TEMP)
  call write_cif_line(trim(import_string))

  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)

  if(EXP_SCAN%nb /=0) then
   WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_full)
   WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
   call center_CIF_line('Data collection')
   WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)

   write(import_string, '(2a)') "Software        : ", Experiment%ACQ_soft
   call write_cif_line(trim(import_string))

   write(import_string, '(a,i2)')    "Number of scans : ",EXP_SCAN%nb
   call write_cif_line(trim(import_string))
   if(EXP_SCAN%STL) then
    write(import_string, '(a,i2,a)') "Option          : shutterless mode"
    call write_cif_line(trim(import_string))
   end if

   if(EXP_scan%sfrm_exist(1)) then
    call write_info('')
    call write_cif_line('')
    write(import_string, '(a)') "    Scan  Time(s)  Width  DX (mm)  Frames   Theta   Omega     Phi     Chi"
    call write_cif_line(trim(import_string))
    call write_info(trim(import_string))
   end if

   do i=1, EXP_SCAN%nb
    if(EXP_scan%sfrm_exist(i)) then
     if(DEVICE%D8V .and. EXP_SCAN%nb_frames(i) == 180 .and. EXP_scan%theta(i) < 0.01 .and. EXP_scan%phi(i) < 0.01) &
        EXP_scan%axis(i) = 'Fast'
     if(i<10) then
      write(import_string, '(1x,i1,1x,a,2x,F6.1,1x,F6.2, 5x,F4.1, 4x,i4, 4F8.2)')  &
            i, EXP_scan%axis(i)(1:6), EXP_SCAN%exposition_time(i), EXP_SCAN%frame_width(i), EXP_SCAN%DX(i), EXP_SCAN%nb_frames(i), &
               EXP_scan%theta(i), EXP_scan%omega(i), EXP_scan%phi(i), EXP_scan%chi(i)
     else
      write(import_string, '(1x,i2,1x,a,2x, F6.1,1x,F6.2, 5x,F4.1, 4x,i4, 4F8.2)')  &
            i, EXP_scan%axis(i)(1:6), EXP_SCAN%exposition_time(i), EXP_SCAN%frame_width(i), EXP_SCAN%DX(i), EXP_SCAN%nb_frames(i), &
               EXP_scan%theta(i), EXP_scan%omega(i), EXP_scan%phi(i), EXP_scan%chi(i), EXP_scan%phi(i), EXP_scan%chi(i)
     endif
     call write_cif_line(trim(import_string))
     call write_info(trim(import_string))

    end if

   end do
  else
   write(import_string, '(a)') "Number of scans        : ?"
   call write_cif_line(trim(import_string))
   !call write_info(trim(import_string))
  end if
  call write_info("")

  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_full)
  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
  call center_CIF_line('Data reduction')
  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)

  call write_cif_line('Input files:')

  write(string, '(3a)') " . P4P file (SAINT ",trim(SAINT%version), ")"
  write(import_string, '(3a)') string(1:38), " : ", TRIM(P4P_file_name)
  !write(import_string, '(4a)') " . P4P file (SAINT ",trim(SAINT%version), ") : ", TRIM(P4P_file_name)
  call write_cif_line(trim(import_string))

  if(len_trim(HKL_file_name) == 0) then
   i1 = INDEX(P4P_file_name, '.')
   SHELX_HKL_file = P4P_file_name(1:i1) //'HKL'
  else
   SHELX_HKL_file = trim(HKL_file_name)
  end if

  if(len_trim(RAW_file_name) /=0) then
   write(string, '(3a)') " . RAW file (SAINT ",trim(SAINT%version), ")"
   write(import_string, '(3a)') string(1:38), " : ", TRIM(RAW_file_name)
   !write(import_string, '(4a)') " . RAW file (SAINT ",trim(SAINT%version), ") : ", TRIM(RAW_file_name)
   call write_cif_line(trim(import_string))
  else
   if(cos_exist) then
    write(string, '(3a)') " . RAW file (SAINT ",trim(SAINT%version), ")"
    write(import_string, '(3a)') string(1:38), " : ", TRIM(SHELX_HKL_file)
    !write(import_string, '(4a)') " . RAW file (SAINT ",trim(SAINT%version), ") : ", TRIM(SHELX_HKL_file)
    call write_cif_line(trim(import_string))
   else
    if(len_trim(SADABS%version) == 0) then
     write(string, '(3a)') " . HKL file (",trim(SADABS%name), ")"
     write(import_string, '(3a)') string(1:38), " : ", TRIM(SHELX_HKL_file)
     !write(import_string, '(4a)') " . HKL file (", trim(SADABS%name), ")       : ", TRIM(SHELX_HKL_file)
    else
     write(string, '(5a)') " . HKL file (", trim(SADABS%name), " ", trim(SADABS%version), ")"
     write(import_string, '(3a)') string(1:38), " : ", TRIM(SHELX_HKL_file)
     !write(import_string, '(6a)') " . HKL file (", trim(SADABS%name), " ", trim(SADABS%version), "): ", TRIM(SHELX_HKL_file)
    endif
    call write_cif_line(trim(import_string))
   end if
  endif

  !if(.not. TWINABS_file) then
   if(len_trim(ABS_file_name) /=0) then
    call test_file_exist(trim(ABS_file_name), file_exist, 'no_out')
    if(file_exist) then
     write(string, '(5a)') " . ABS file (", trim(SADABS%name), " ", trim(SADABS%version), ")"
     write(import_string, '(3a)') string(1:38), " : ", TRIM(ABS_file_name)
     !write(import_string, '(6a)') " . ABS file (", trim(SADABS%name), " ", trim(SADABS%version), "): ", TRIM(ABS_file_name)
     call write_cif_line(trim(import_string))
    end if
   end if
   WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)

   if(SADABS%point_group(1:1) /= '?') then
    if(.not. TWINABS_file) then
     write(import_string, '(2a)')   "SADABS point group for merging  :", trim(SADABS%point_group)
    else
     write(import_string, '(2a)')   "Equivalent reflections defined according to point group  :", trim(SADABS%point_group)
    end if
    call write_cif_line(trim(import_string))
    WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
   end if
  !end if

  !if(.not. twinabs_4) then
   if(SADABS%nb_ref /= 0) then
    write(import_string, '(a,i8)') "Number of corrected reflections : ", SADABS%nb_ref
    call write_cif_line(trim(import_string))
    WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
   end if

  !else
   if(twinabs_4 .and. SADABS%nb_obs /=0 .and. SADABS%nb_obs_3s /=0) then
    call write_cif_line('')
    write(import_string, '(a,F7.4,a,I8,a)') "Rint = ", SADABS%Rint, ' for all ', SADABS%nb_obs, ' observations and'
    call write_cif_line(trim(import_string))
    write(import_string, '(a,F7.4,a,I8,a)') "Rint = ", SADABS%Rint_3s, ' for all ', SADABS%nb_obs_3s, &
                                            ' observations with I > 3sigma(I)'
    call write_cif_line(trim(import_string))
    WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
   end if
  !end if



  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_empty)
  WRITE(CIF_unit, '(a)')  trim(CIF_sep_line_full)

  WRITE(CIF_unit, '(a)')  ""
  !if(CIF_parameter%sample_ID(1:1) == "?") then
   WRITE(CIF_unit, '(a)')  "data_collect"
  !else
  ! WRITE(CIF_unit, '(2a)')  "data_", trim(CIF_parameter%sample_ID)
  !endif
  WRITE(CIF_unit, '(a)')  ""


  return
 end subroutine Create_import_cif

!--------------------------------------------------------------------------------------------------------
subroutine Get_experiment_features
  use Cryscalc_module, only : SAINT, tmp_unit, DEVICE, experiment, EXP_scan
  USE CIF_module
  use macros_module,   only : l_case, replace_car
  use IO_module,       only : write_info

  character (len=256)               :: sfrm_file
  character (len=16)                :: new_string
  logical                           :: file_exist
  integer                           :: i, i1, i_error
  character (len=8192)              :: read_line
  real                              :: x1, x2, x3, x4


  if(DEVICE%D8V) then
   write(sfrm_file, '(4a)') trim(SAINT%source), '_', trim(CIF_parameter%sample_ID), '_01_0001.sfrm'
  else
   write(sfrm_file, '(3a)') trim(CIF_parameter%sample_ID), '_01_0001.sfrm'
  endif

  sfrm_file = '..\'//trim(sfrm_file)

  inquire(file = trim(sfrm_file), exist=file_exist)
  if(file_exist) then
   ! date of file
   !call system("dir /TW "//trim(sfrm_file)// " > temp")
   !open(unit=tmp_unit, file='temp')
   !do i =1, 6
   ! read(tmp_unit, '(a)') read_line
   !end do
   !Experiment_date = read_line(1:10)
   !close(unit=tmp_unit)
   !call system("del temp")

   ! lecture entete fichier .SFRM
   ! recup date, temperature, dx
   open(unit=tmp_unit, file=trim(sfrm_file))
    do
     read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(index(read_line, 'CREATED:') /=0)  then
     i1 = index(read_line, 'CREATED:') + 8
     read(read_line(i1:), '(a)') Experiment%date
     if(DEVICE%D8V) then
      if(index(Experiment%date, '-') /=0) then
       Experiment%date = replace_car(Experiment%date, '-', '/')
       if(index(Experiment%date, 'Jan') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Jan', '01')
       elseif(index(Experiment%date, 'Feb') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Feb', '02')
       elseif(index(Experiment%date, 'Mar') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Mar', '03')
       elseif(index(Experiment%date, 'Apr') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Apr', '04')
       elseif(index(Experiment%date, 'May') /=0) then
        Experiment%date = replace_car(Experiment%date, 'May', '05')
       elseif(index(Experiment%date, 'Jun') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Jun', '06')
       elseif(index(Experiment%date, 'Jul') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Jul', '07')
       elseif(index(Experiment%date, 'Aug') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Aug', '08')
       elseif(index(Experiment%date, 'Sep') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Sep', '09')
       elseif(index(Experiment%date, 'Oct') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Oct', '10')
       elseif(index(Experiment%date, 'Nov') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Nov', '11')
       elseif(index(Experiment%date, 'Dec') /=0) then
        Experiment%date = replace_car(Experiment%date, 'Dec', '12')
       end if
      end if
     else
      new_string=Experiment%date(4:5)//'/'//Experiment%date(1:2)//"/20"//Experiment%date(7:8)
      Experiment%date = new_string
     end if
    end if

    !temperature
    !if(index(read_line, 'TEMP:', back=.true.) /=0)  then
    ! i1 = index(read_line, 'TEMP:', back=.true.) + 5
    ! read(read_line(i1:), *) x1
    ! if(x1>0.1) then
    !  write(experiment%TEMP, '(F8.2)') x1
    !  CIF_cell_measurement%temperature = Experiment%TEMP
    !  experiment%TEMP = trim(experiment%TEMP)//' K'
    ! endif
    !end if
    if(index(read_line, 'LOWTEMP:') /=0)  then
     i1 = index(read_line, 'LOWTEMP:') + 8
     read(read_line(i1:), *) x1, x2, x3
     x2 = x2/100. + 273.15
     write(experiment%TEMP, '(F8.2)') x2
     CIF_cell_measurement%temperature = Experiment%TEMP
     Experiment%TEMP = adjustl(Experiment%TEMP)
     if(experiment%TEMP(1:4) == "0.00") then
      experiment%TEMP = "room temp."
      experiment%TEMP = "295.00"
     else
      experiment%TEMP = trim(experiment%TEMP)//' K'
     end if
     x3 = x3/100.
     write(experiment%detector_TEMP, '(F8.2)') x3
     experiment%detector_TEMP = trim(experiment%detector_TEMP) // ' C'
     experiment%detector_TEMP = adjustl(experiment%detector_TEMP)
    end if


    !if(index(read_line, 'DISTANC:') /=0)  then
    !i1 = index(read_line, 'DISTANC:') + 8
    ! read(read_line(i1:), *) x1, x2
    ! x2 = 10.*x2
    ! write(experiment%DX, '(F5.2)') X2
    !end if

    if(index(read_line, 'TARGET :') /=0)  then
     i1 = index(read_line, 'TARGET :') + 8
     read(read_line(i1:), *) experiment%target
    end if

   if(index(read_line, 'SOURCEK:') /=0)  then
    i1 = index(read_line, 'SOURCEK:') + 8
    read(read_line(i1:), *) x1
    write(experiment%voltage, '(F5.2,a)') x1, ' KV'
   end if

   if(index(read_line, 'SOURCEM:') /=0)  then
    i1 = index(read_line, 'SOURCEM:') + 8
    read(read_line(i1:), *) x1
    write(experiment%current, '(F5.2,a)') x1, ' mA'
   end if


   if(index(read_line, 'PROGRAM:') /=0)  then
    i1 = index(read_line, 'PROGRAM:') + 8
    read(read_line(i1:i1+48), '(a)') experiment%ACQ_soft
   end if

   if(index(read_line, 'DETTYPE:') /=0)  then
    i1 = index(read_line, 'DETTYPE:') + 8
    read(read_line(i1:), '(a)') experiment%detector
    exit
   end if

  end do
  close(unit=tmp_unit)
  else
   call write_info('')
   call write_info(' > ' // trim(SFRM_file) // ' does not exist! Scan features can not be extracted.')
   return
  end if

 ! recuperation de divers parametres dans .sfrm :
 !    . type de scan
 !    . temps d'exposition / frame
 !    . largeur angulaire / frame
 !    . valeurs des angles de depart de theta, omega, phi et chi
 !    . distance cristal-detecteur

  do i=1, EXP_SCAN%nb
   if(DEVICE%D8V) then
    if(i < 10) then
     write(sfrm_file, '(4a,I1,a)') trim(SAINT%source), '_', trim(CIF_parameter%sample_ID), '_0', i, '_0001.sfrm'
    else
    end if
   else
    write(sfrm_file, '(2a,I1,a)') trim(CIF_parameter%sample_ID), '_0', i, '_0001.sfrm'
   endif
   sfrm_file = '..\'//trim(sfrm_file)
   inquire(file = trim(sfrm_file), exist=file_exist)
   if(file_exist) then
    EXP_scan%sfrm_exist(i) = .true.
    open(unit=tmp_unit, file=trim(sfrm_file))
     do
      read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
      if(i_error /=0) exit

      !if(.not. DEVICE%D8V) then
      !if(index(read_line, 'TYPE   :') /=0)  then
      ! i1 = index(read_line, 'TYPE   :') + 8
      ! !read(read_line(i1:), '(a)') EXP_scan%axis(i)
      ! read(read_line(i1:), *) EXP_scan%axis(i)
      !end if
      !end if

      if(index(read_line, 'AXIS   :') /=0)  then
       i1 = index(read_line, 'AXIS   :') + 8
       read(read_line(i1:), *) EXP_scan%axis(i)
       if(EXP_scan%axis(i)(1:1) == "1") then
        EXP_scan%axis(i) = "Theta"
       elseif(EXP_scan%axis(i)(1:1) == "2") then
        EXP_scan%axis(i) = "Omega"
       elseif(EXP_scan%axis(i)(1:1) == "3") then
        EXP_scan%axis(i) = "Phi"
       elseif(EXP_scan%axis(i)(1:1) == "4") then
        EXP_scan%axis(i) = "Chi"
       end if
      end if


      if(index(read_line, 'CUMULAT:') /=0)  then
       i1 = index(read_line, 'CUMULAT:') + 8
       read(read_line(i1:), *) x1
       EXP_SCAN%exposition_time(i) = x1
      end if

       if(index(read_line, 'RANGE  :') /=0)  then
       i1 = index(read_line, 'RANGE  :') + 8
       read(read_line(i1:), *) x1
       EXP_SCAN%frame_width(i) = x1
      end if


      if(index(read_line, 'ANGLES :') /=0)  then
       i1 = index(read_line, 'ANGLES :') + 8
       read(read_line(i1:), *) x1,x2,x3,x4
       EXP_SCAN%theta(i) = x1
       EXP_SCAN%omega(i) = x2
       EXP_SCAN%phi(i) = x3
       EXP_SCAN%chi(i) = x4
      end if

      if(index(read_line, 'DISTANC:') /=0)  then
       i1 = index(read_line, 'DISTANC:') + 8
       read(read_line(i1:), *) x1, x2
       x2 = 10.*x2
       EXP_SCAN%DX(i) = X2
      end if

     end do
    close(unit=tmp_unit)
   else
    EXP_scan%sfrm_exist(i) = .false.
    call write_info(' > ' // trim(SFRM_file) // ' does not exist! Scan features can not be extracted.')
    return
   end if


  end do

  return
 end subroutine Get_experiment_features

!---------------------------------------------------------------------------------------
 subroutine test_ISCR(proposer_team, ISCR)
  USE Macros_module, only : u_case
  implicit none
  character (len=*), intent(in)     :: proposer_team
  logical,           intent(inout)  :: ISCR
  integer                           :: long

  ISCR = .false.
  long = len_trim(proposer_team)


  if(long == 3) then
   if(u_case(proposer_team(1:3)) == 'COS'  .or. &
      u_case(proposer_team(1:3)) == 'CPM'  .or. &
      u_case(proposer_team(1:3)) == 'CSM'  .or. &
      u_case(proposer_team(1:3)) == 'CTI'  .or. &
      u_case(proposer_team(1:3)) == 'OMC')   then
    ISCR = .true.
    return
   endif

  elseif(long == 4) then
   if(u_case(proposer_team(1:4)) == 'ICMV')   then
    ISCR = .true.
    return
   endif

  elseif(long == 5) then
   if(u_case(proposer_team(1:5)) == 'CDIFX'  .or. &
      u_case(proposer_team(1:5)) == 'MACSE'  .or. &
      u_case(proposer_team(1:5)) == 'PNSCM')   then
     ISCR = .true.
     return
   endif

  elseif(long == 6) then
   if(u_case(proposer_team(1:6)) == 'CORINT')   then
     ISCR = .true.
     return
   endif

  end if



  return
 end subroutine test_ISCR


!---------------------------------------------------------------------------------------
 subroutine Get_team_ISCR(proposer_name, proposer_team)
  USE Cryscalc_module, only : ISCR
  USE Macros_module,   only : u_case
  implicit none
  character (len=*), intent(in)     :: proposer_name
  character (len=*), intent(inout)  :: proposer_team
  integer                           :: i, long

  call ISCR_effectifs

  do i=1, ISCR%nb
   long = len_trim(proposer_name)
   if(u_case(proposer_name(1:long)) == u_case(ISCR%name(i))) then
    Proposer_team = ISCR%team(i)
    return
   end if
  end do


  !long = len_trim(proposer_name)
  !if(long == 10) then
  ! if(u_case(proposer_name(1:10)) == 'J CRASSOUS') then
  !  Proposer_team = 'OMC'
    !return
   !end if
  !end if



  return
 end subroutine Get_team_ISCR


!--------------------------------------------------------------------------------------------
subroutine ISCR_effectifs
  ! lecture du fichier ISCR_effectif.txt (effectif ISCR)
  use Cryscalc_module, only : CRYSCALC, ISCR, tmp_unit
  implicit none
   integer             :: i, long, i_error
   logical             :: file_exist
   character (len=256) :: read_line


  call getenv('CRYSCALC', cryscalc%path_name)

  long = len_trim(cryscalc%path_name)
  if (long ==0) return

  if(cryscalc%path_name(1:1) == '"' .and. cryscalc%path_name(long:long) == '"') then
   cryscalc%path_name = cryscalc%path_name(2:long-1)
   long = len_trim(cryscalc%path_name)
  end if

  if(cryscalc%path_name(long:long) == '\') then
   cryscalc%path_name = cryscalc%path_name(1:long-1)
  endif
  ISCR%effectifs  = trim(cryscalc%path_name)//'\ISCR_effectif.txt'

  inquire (FILE=TRIM(ISCR%effectifs), EXIST=file_exist)
  if(.not. file_exist) return

  open (unit = tmp_unit, file=trim(ISCR%effectifs))
   i = 0
   do
    read(unit = tmp_unit, fmt='(a)', iostat = i_error) read_line
    if(i_error /=0) exit
    read_line = adjustl(read_line)
    if(read_line(1:1) == '!' .or. read_line(1:1) == '#') cycle
    if(len_trim(read_line) == 0) cycle
    i = i + 1
    read_line = adjustl(read_line)
    read(read_line(1:16),  '(a)') ISCR%name(i)
    read(read_line(40:50), '(a)') ISCR%team(i)
    ISCR%name(i) = trim(ISCR%name(i))
    ISCR%team(i) = trim(ISCR%team(i))

   end do
  close (unit = tmp_unit)
  ISCR%nb = i

  return
 end subroutine ISCR_effectifs
