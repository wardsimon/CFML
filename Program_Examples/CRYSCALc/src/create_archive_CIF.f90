subroutine read_cif_and_create_ARCHIVE
! read CIF files and create ARCHIVE.cif 
 use cryscalc_module, only : nb_input_CIF_files, input_CIF_file, CIF_file_exist, CIF_unit, CRYSCALC, CIF_read_unit, &
                             molecule, absorption, SADABS_absorption_coeff, SPG, device, SFAC, nb_atoms_type,       &
							 include_RES_file, SQUEEZE, CIF_format80, message_text							 
 use CIF_module 
 use macros_module,   only : Get_current_folder_name, Get_sample_ID_from_folder_name, u_case, test_file_exist
 use IO_module,       only : write_info
 implicit none
  integer                           :: i, i1, i2, long, i_error
  character (len=64)                :: fmt_, fmt_final
  character (len=256)               :: current_folder
  character (len=256)               :: job, sample_ID, diffraction_date  
  character (len=256)               :: CIF_sep_line
  character (len=256)               :: read_line
  character (len=256)               :: CIF_string
  character (len=16), dimension(16) :: CIF_str
  CHARACTER (LEN=256)               :: string_value
  logical                           :: file_exist, ok
  integer                           :: nb_files_exist
  
  ! test si les fichiers .CIF existent
  nb_files_exist = 0
  do i=1, nb_input_CIF_files
   call test_file_exist(input_CIF_file(i), file_exist, 'out')
   if(file_exist) then
    CIF_file_exist(i) = .true.
	nb_files_exist =  nb_files_exist + 1
   end if	
  end do
  if(nb_files_exist == 0) then
   call write_info('')
   call write_info('CIF files does not exist. Program will be stopped !')
   call write_info('')
   stop
  end if
  
  
  write(CIF_sep_line, '(a,76a1,a)') '#', ('-',i=1,76), '#'
  
  
  
  close(unit=CIF_unit)
  open(unit=CIF_unit, file="cryscalc_archive.cif")

  write(CIF_unit, '(a)') '#'
  write(CIF_unit, '(a)') '#'

  write(CIF_unit, '(a)') trim(CIF_sep_line)
  long = len_trim(CRYSCALC%version)
  write(fmt_ , '(a,i2,a)') '(3a,', 30-long, 'x,a)'
  write(CIF_unit, fmt=trim(fmt_))        '# CIF file created and formatted by CRYSCALC (', trim(CRYSCALC%version),')','#'
  write(CIF_unit, '(a,59x,a)')           '# Input CIF files:','#'
  call write_info('   Input CIF files:')
  do i=1, nb_input_CIF_files
   long=len_trim(input_cif_file(i))
   write(fmt_ , '(a,i2,a)') '(a,i2,1x,a,', 72-long, 'x,a)'
   write(CIF_unit, fmt=trim(fmt_))        '# ', i, trim(Input_CIF_file(i)),'#'
   write(message_text, fmt='(5x,i2,2a)')    i, '. ' , trim(Input_CIF_file(i))  
   call write_info(trim(message_text))
  end do
	
	
  call Get_current_folder_name(current_folder)

  diffraction_date = ''
  i1 = index(current_folder, '_', back=.true.)
  if (i1 /=0) diffraction_date = current_folder(i1+1:)
  i1 = index(diffraction_date, '\')
  if (i1 /=0) diffraction_date = diffraction_date(1:i1-1)

  if(len_trim(diffraction_date) /= 0) then
	long = len_trim(diffraction_date)
	write(fmt_ , '(a,i2,a)') '(2a,', 41-long, 'x,a)'
    write(CIF_unit, fmt = trim(fmt_)) 	   '#  Date of diffraction experiment : ', trim(diffraction_date),"#"	 	
  end if 
  write(CIF_unit, '(a)') trim(CIF_sep_line)
  write(CIF_unit, '(a)') '#'
  write(CIF_unit, '(a)') ''
 
  call Get_sample_ID_from_folder_name(current_folder, sample_ID)
  !if(sample_ID == '?') call Get_Wingx_job(sample_ID, wingx_structure_dir)
  if(sample_ID == '?') sample_ID = trim(read_line(6:))
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
  open(unit=CIF_read_unit, file = trim(input_CIF_file(1)))

  CIF_string = '_cell_formula_units_Z'
  ok = .false.
  call get_champ_value(CIF_read_unit, trim(CIF_string), string_value, ok)
  IF(ok) READ(string_value, '(a)') CIF_parameter%cell_formula_units_Z  
  
  ! determination du nombre de champ dans la boucle _atom_site_
  call Get_CIF_champ_nb(CIF_read_unit, '_atom_site_label', CIF_atom_site_item_nb, CIF_atom_site_loop_numor)

  ! determination des positions des differents items dans la boucle
  call Get_CIF_numors(CIF_read_unit)
  
  write(CIF_unit, '(3a)')     "_chemical_formula_moiety         '", trim(adjustl(molecule%formula)),"'"
  write(CIF_unit, '(3a)')     "_chemical_formula_sum            '", trim(adjustl(molecule%formula)),"'"
  !write(CIF_unit, '(a,F6.2)') "_chemical_formula_weight          ", molecule%weight
  write(CIF_unit, '(2a)')     "_chemical_formula_weight          ", trim(CIF_parameter%formula_weight)
  write(CIF_unit, '(a)')      "_chemical_compound_source        'synthesis as described'"
  
  call write_CIF_file('UNIT_CELL_INFO')
  
  call GET_SG_from_CIF(trim(input_CIF_file(1)))    ! le bon groupe est dans le fichier .CIF principal
  
  call write_CIF_file('SPACE_GROUP')
  call write_CIF_file('CELL_PARAM_ESD')
  call write_CIF_file('CIF_THMIN,THMAX')
  
  !call write_CIF_file('CELL_WAVE')
  write(CIF_unit, '(a,F10.5)')      "_cell_measurement_wavelength        ", CIF_cell_measurement%wavelength
  
  call write_CIF_file('CRYSTAL_INFORMATION')
  call write_CIF_file('CRYSTAL_DENSITY')
  call write_CIF_file('ABSORPTION')
  read(SADABS_absorption_coeff, *) absorption%mu 
  absorption%mu = 10. * absorption%mu 
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
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_wavelength            ", trim(CIF_parameter%diffrn_radiation_wavelength)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_type                  ", trim(CIF_parameter%diffrn_radiation_type)
  WRITE(CIF_unit, '(3a)') "_diffrn_radiation_source                '", trim(CIF_parameter%diffrn_radiation_source),"'"
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_monochromator         ", trim(CIF_parameter%diffrn_radiation_monochromator)
  WRITE(CIF_unit, '(2a)') "_diffrn_radiation_probe                 ", trim(CIF_parameter%diffrn_radiation_probe)
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
  WRITE(CIF_unit, '(2a)') "_diffrn_measured_fraction_theta_full    ", trim(CIF_parameter%diffrn_theta_full)
  WRITE(CIF_unit, '(2a)') "_diffrn_measured_fraction_theta_max     ", trim(CIF_parameter%diffrn_theta_max)
  WRITE(CIF_unit, '(2a)') "_reflns_number_total                    ", trim(CIF_parameter%reflns_number_total)
  WRITE(CIF_unit, '(2a)') "_reflns_number_gt                       ", trim(CIF_parameter%reflns_number_gt)
  WRITE(CIF_unit, '(2a)') "_reflns_threshold_expression            ", ">2sigma(I)"
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
  
  if(include_RES_file) then
   WRITE(CIF_unit, '(a)')  ""
   call include_RES_file_into_CIF
   write(CIF_unit, '(a)') trim(CIF_sep_line)
  end if
  call write_CIF_file('ATOMS_HEADER')
  
  
  ! lecture fichier .CIF 
  close(unit=CIF_read_unit)
  open(unit = CIF_read_unit, file = trim(input_CIF_file(1)))
  do   
   read(unit = CIF_read_unit, fmt='(a)', iostat=i_error) read_line
   if(i_error /=0) exit
   
   if(index(read_line, '_atom_type_scat_source') /=0)  then
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
	 read(read_line, *) CIF_str(1:4)
	 do i=1, 2
	  i1 = index(CIF_str(i), "'")
	  i2 = index(CIF_str(i), "'", back=.true.)
	  if(i1 /=0 .and. i2 /=0 .and. i2>i1) CIF_str(i) = CIF_str(i)(i1+1:i2-1)
	 end do 
	 read(CIF_read_unit, '(a)', iostat=i_error) read_line
	 if(i_error /=0) exit
     
	 if(CIF_str(3)(1:1) /= '-') then
	 write(CIF_unit, '(2a4, 1x,2a8, a)') CIF_str(1:4), trim(read_line)
	 else
     write(CIF_unit, '(2a4, a8,1x,a8, a)') CIF_str(1:4), trim(read_line)
	 endif
    end do 
    write(CIF_unit, '(a)') ''   
    
   elseif(index(read_line, '_atom_site_disorder_group') /=0)  then
    write(CIF_unit, '(a)') 'loop_'
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
    do
     read(CIF_read_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
	 if(len_trim(read_line) == 0) exit
	 read(read_line, *) CIF_str(1:13)
	 if(CIF_format80) then 
	  call get_fmt(4, CIF_str(3:6), fmt_, 0)
	  if(index(CIF_str(CIF_atom_site_refinement_flags_numor), 'P') /=0) then ! occ. partielle
       if(len_trim(CIF_str(CIF_atom_site_symm_multiplicity_numor)) == 1) then		
	    write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a6,a2,a5,a4,2a2)' 
	   else
        write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a6,a3,a5,a4,2a2)' 	   
       endif

	  else
	   if(len_trim(CIF_str(CIF_atom_site_symm_multiplicity_numor)) == 1) then
	    write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2,a2,a5,a3,2a2)'            
	   else
        write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2,a3,a5,a3,2a2)'            
	   end if 
	  end if 
      write(CIF_unit, fmt_final) CIF_str(1:13)
     else	  
	  write(CIF_unit, '(a)') trim(read_line)
	 end if	 
	end do
	write(CIF_unit, '(a)') ''
 
   elseif(index(read_line, '_atom_site_aniso_U_12') /=0)  then
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
     write(fmt_final, '(3a)') '(4a5', trim(fmt_), ',x,5a8)'	 
     write(CIF_unit, fmt_final) CIF_str(1:10)     
    end do
	write(CIF_unit, '(a)') ''

   elseif(index(read_line, '_geom_hbond_site_symmetry_A') /=0) then   
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
  write(fmt_ , '(a,i2,a)') '(2a,', 53-long, 'x,a)'
  write(CIF_unit, fmt=trim(fmt_))        '#===END OF CIF FILE FOR ', trim(u_case(sample_ID)),'#'
  write(CIF_unit, '(a)') trim(CIF_sep_line)    
  
  call write_info('   >>> CRYSCALC_ARCHIVE.cif file has been created.')
  call write_info('')
  stop
	
 return
end subroutine read_cif_and_create_ARCHIVE
!----------------------------------------------------------------------------------------------------------


subroutine Include_RES_file_into_CIF
 use CRYSCALC_module, only  : input_CIF_file, CIF_unit, tmp_unit, include_HKL_file
 use macros_module,   only  : test_file_exist
 use IO_module,       only  : write_info
 implicit none
 integer                         :: i1, i_error
 CHARACTER (LEN=256)             :: job_RES_file, job_HKL_file
 LOGICAL                         :: file_exist
 CHARACTER (LEN=256)             :: read_line


   i1 = index(input_CIF_file(1), '.')
   if(i1 == 0) return
   
   job_RES_file = input_CIF_file(1)(1:i1)//'res'   
   call test_file_exist(trim(job_RES_file), file_exist, 'out')
   if(file_exist) then
	call write_info('')
	call write_info('  ... Include '//trim(job_RES_file)//' SHELXL file.')
	call write_info('')
	write(CIF_unit, '(a)')   '_iucr_refine_instructions_details'
	write(CIF_unit, '(a)')   ';'
	open(unit=tmp_unit, file=trim(job_res_file))
	write(CIF_unit, '(2a)') ' .res SHELXL file: ', trim(job_RES_file)
	write(CIF_unit, '(a)')  '.................................................................'
	do
	 read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
	 if(i_error /=0) exit
	 write(CIF_unit, '(a)') trim(read_line)
    end do 
	close(unit=tmp_unit)	
	write(CIF_unit, '(a)')   ';'
	write(CIF_unit, '(a)')   '# End of res file'
	
    if(include_HKL_file) then
	 job_HKL_file = input_CIF_file(1)(1:i1)//'hkl'   
     call test_file_exist(trim(job_HKL_file), file_exist, 'out')
     if(file_exist) then
	  call write_info('')
	  call write_info('  ... Include '//trim(job_HKL_file)//' SHELXL file.')
	  call write_info('')
	  write(CIF_unit, '(a)') ''
	  !write(CIF_unit, '(a)') '_shelx_hkl_file'
	  write(CIF_unit, '(a)') '_iucr_refine_reflections_details'
	  write(CIF_unit, '(a)') ';'
	  open(unit=tmp_unit, file=trim(job_HKL_file))
	   write(CIF_unit, '(2a)') ' .hkl SHELXL file: ', trim(job_HKL_file)
	   write(CIF_unit, '(a)')  '.................................................................'
	   do
	    read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
	    if(i_error /=0) exit
	    write(CIF_unit, '(a)') trim(read_line)
	   end do 		
	  close(unit=tmp_unit)
	  write(CIF_unit, '(a)')   ';'
	  write(CIF_unit, '(a)')   '# End of hkl file'
	 endif
	end if
   endif	

  return	  
end subroutine Include_RES_file_into_CIF  	


!------------------------------------------------------------------
subroutine GET_SG_from_CIF(input_file)
 USE cryscalc_module,                  ONLY : SPG, nb_symm_op, on_screen
 USE macros_module,                    ONLY : replace_car
 USE CFML_String_Utilities,            ONLY : number_lines, reading_lines
 USE CFML_IO_Formats,                  ONLY : read_cif_hm, read_cif_symm, read_cif_hall
 USE CFML_Crystallographic_Symmetry,   ONLY : set_spacegroup, get_hallsymb_from_gener
 USE IO_module,                        ONLY : write_info
 
 implicit none
  CHARACTER (LEN=*), INTENT(IN)               :: input_file
  character(len=40)                           :: symb_spgr  
  integer                                     :: npos, nb_lines, i
  character(len=256),dimension(:),allocatable :: file_cif
  character(len=40), dimension(48)            :: car_symop
  
  call number_lines(trim(input_file),nb_lines)
  if (nb_lines ==0) then
   call write_info(' > No lines could be read. Program will be stopped')
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
	
    if(SPG%NumSpg == 0 .and. on_screen) then
	 call write_info('')
	 call write_info('     >> Space group can not be deduced from symmetry operators list !!')   
	!else
	! call write_info('')
	! call write_info('     >> Space group deduced from symmetry operators list !!')   
	end if
   end if
   
  end if

   return
end subroutine GET_SG_from_CIF