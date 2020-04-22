
module structural_report_module
  implicit none
  CHARACTER (LEN=256)                 :: HTML_structural_report_file
  CHARACTER (LEN=256)                 :: TEXT_structural_report_file
  CHARACTER (LEN=256)                 :: LATEX_structural_report_file
  CHARACTER (LEN=256)                 :: PDF_structural_report_file
  CHARACTER (LEN=32)                  :: job
  character (len=256), dimension(192) :: op_string

end module structural_report_module

!---------------------------------------------------------------------------------------------------------------
module LATEX_module
  implicit none

 character (len=256), dimension(2)  :: logo


end module LATEX_module

!---------------------------------------------------------------------------------------------------------------
module special_characters_module

  character (len=32), parameter         :: HTML_elect_density_unit      = " e<sup>-</sup>.Å<sup>-3</sup>"
  character (len=32), parameter         :: LATEX_elect_density_unit     = " $e^-.\AA^{-3}$"
  character (len=32), parameter         :: LATEX_alt_elect_density_unit = "  \(\overline{e}\).\AA\(\sp{-3}\)"


end module special_characters_module

!---------------------------------------------------------------------------------------------------------------
subroutine create_structural_report
 USE cryscalc_module, only :  debug_proc

 USE CIF_module

 implicit none

 if(debug_proc%level_2)  call write_debug_proc_level(2, "create_structural_report")

  call Get_CIF_parameters
  call create_report_header
  call create_report


  !call create_structural_study


 return
end subroutine create_structural_report

!---------------------------------------------------------------------------------------------------------------
subroutine Get_CIF_parameters()
 USE structural_report_module,        ONLY : job, op_string
 USE MACROS_module,                   ONLY : test_file_exist,  l_case
 USE cryscalc_module,                 ONLY : CIF_unit, archive_CIF, debug_proc, include_squeeze, SPG, debug_proc
 use CIF_module,                      ONLY : CIF_parameter, CIF_parameter_DEVICE
 USE CFML_Crystallographic_Symmetry,  ONLY : set_spacegroup
 USE IO_module

 implicit none
 LOGICAL                             :: file_exist
 CHARACTER (LEN=256)                 :: CIF_input_line, CIF_field, CIF_field_value
 INTEGER                             :: i_error, long_field, long, i1
 integer                             :: n_op

 CHARACTER (LEN=64)                  :: device_mark, device_type
 LOGICAL                             :: symm_id

 if(debug_proc%level_2)  call write_debug_proc_level(2, "get_CIF_parameter")


 symm_id = .false.
 include_squeeze = .false.

 file_exist = .false.
 call test_file_exist(TRIM(archive_CIF), file_exist,'out')
 IF(.NOT. file_exist) then
  !call write_info('')
  !call write_info('   archive.cif file is missing !!')
  !call write_info('')
  if(debug_proc%write) then
   call write_debug_proc(' '//trim(archive_cif)//' file is missing !', '')
   call write_debug_proc(' ', '')
   call write_debug_proc(' CRYSCALC will be stopped.', '')
   call write_debug_proc('', '')
  end if
  stop
 endif



 open (unit = CIF_unit, file = TRIM(archive_CIF))
  do
   READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
   IF(i_error < 0) exit
   IF(i_error /=0) then
    call write_info('')
    call write_info(' !! Error reading archive.CIF file !!')
    call write_info('')
    return
   endif
   IF(LEN_TRIM(CIF_input_line) == 0) cycle
   if(CIF_input_line(1:32) == '#  CIF produced by WinGX routine') CIF_parameter%WinGX_used = .true.
   IF(CIF_input_line(1:1) == "#")    cycle


   !read(CIF_input_line, '(a)', IOSTAT=i_error) CIF_field
   !IF(i_error /=0) then
   ! call write_info('')
   ! call write_info(' !! Error reading archive.CIF file !!')
   ! call write_info('')
   ! return
   !endif
   CIF_field = adjustl(CIF_input_line)
   i1 =index(CIF_field, " ")
   if(i1 > 1) then
    CIF_field =  CIF_field(1:i1-1)
   else
    call write_info('')
    call write_info(' !! Error reading archive.CIF file !!')
    call write_info('')
    return
   endif

   ! june 2012 : get job name
   if(CIF_input_line(1:5) == 'data_') then
    write(job, '(a)') trim(CIF_input_line(6:))
    cycle
   endif
   ! ---------------------------

   CIF_field = l_case(CIF_field)
   long_field = LEN_TRIM(CIF_field)
   CIF_field_value  = CIF_input_line(long_field+1:)
   CIF_field_value  = ADJUSTL(CIF_field_value)
   long = LEN_TRIM(CIF_field_value)


   IF(CIF_field_value(1:1)== "'") then
    if(CIF_field_value(long:long) == "'" ) CIF_field_value = CIF_field_value(2:long-1)
   ENDIF


   select case (CIF_field(1:long_field))

       case ('_chemical_formula_moiety')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%formula_moiety)

       case ('_chemical_formula_sum')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%formula_sum)

       case ('_chemical_formula_weight')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%formula_weight)

       case ('_symmetry_cell_setting', '_space_group_crystal_system')
        if (CIF_parameter%symmetry_cell_setting == '?') then
         call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%symmetry_cell_setting)
        end if

       case ('_symmetry_space_group_name_h-m', '_space_group_name_h-m', '_space_group_name_h-m_alt')
        if(CIF_parameter%symmetry_space_group =='?' ) then
         call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%symmetry_space_group)
        end if

       case ('_symmetry_int_tables_number', '_space_group_it_number')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%symmetry_IT_number)

       case ('_symmetry_equiv_pos_site_id', '_space_group_symop_id')
        symm_id = .true.

       case ('_symmetry_equiv_pos_as_xyz', '_space_group_symop_operation_xyz')
        n_op = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         if(i_error /=0) exit
         if(len_trim(CIF_input_line)==0) exit
         if(CIF_input_line(1:1) == '#' .or. CIF_input_line(1:1) == '!' ) exit
         if(CIF_input_line(1:1) == '_') then
          backspace(unit = CIF_unit)
          exit
         endif
         n_op = n_op + 1
         if(.not. symm_id) then
          read(CIF_input_line, '(a)') op_string(n_op)
         else
          CIF_input_line = adjustl(CIF_input_line)
          op_string(n_op) = CIF_input_line(index(CIF_input_line, ' '):)
         endif
         op_string(n_op) = adjustl(op_string(n_op))
         op_string(n_op) = trim(op_string(n_op))
         if(op_string(n_op)(1:1) == "'") then
          op_string(n_op) = op_string(n_op)(2:len_trim(op_string(n_op))-1) ! enleve les caracteres "'"
         endif
        end do

        if(CIF_parameter%symmetry_space_group == "?") then
         ! creation du groupe d'espace a partir des operateurs de symetrie
         call set_spacegroup("  ", SPG, op_string, n_op,'GEN')
         CIF_parameter%symmetry_space_group  = SPG%spg_symb
         CIF_parameter%symmetry_cell_setting = SPG%CrystalSys
        end if

       case ('_cell_length_a')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_length_a)

       case ('_cell_length_b')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_length_b)

       case ('_cell_length_c')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_length_c)

       case ('_cell_angle_alpha')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_angle_alpha)

       case ('_cell_angle_beta')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_angle_beta)

       case ('_cell_angle_gamma')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_angle_gamma)

       case ('_cell_volume')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_volume)

       case ('_cell_formula_units_z')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%formula_units_Z)

       case ('_exptl_crystal_density_diffrn')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%exptl_density)

       case ('_cell_measurement_reflns_used')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_reflns_used)

       case ('_cell_measurement_theta_min')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_theta_min)

       case ('_cell_measurement_theta_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%cell_theta_max)


       CASE ('_diffrn_reflns_number')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%diffrn_reflns_number)

       case ('_reflns_number_total')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%reflns_number_total)

       case ('_diffrn_reflns_av_r_equivalents')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%diffrn_reflns_av_R_equivalents)

       case ('_diffrn_reflns_av_sigmai/neti', '_diffrn_reflns_av_uneti/neti')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%diffrn_reflns_av_R_sigma)


       case ('_refine_ls_number_parameters')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%refine_ls_number_parameters)

       case ('_refine_ls_number_restraints')
        !READ(CIF_field_value, '(a)') CIF_parameter%restraints_number
         call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%restraints_number)

       case ('_refine_ls_abs_structure_flack')
         call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%abs_structure_Flack)

       case ('_refine_ls_wr_factor_gt')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%refine_ls_wR_factor_gt)

       CASE ('_refine_ls_r_factor_gt')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%refine_ls_R_factor_gt)

       CASE ('_refine_ls_r_factor_all')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%refine_ls_R_factor_all)

       case ('_refine_ls_wr_factor_ref')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%refine_ls_wR_factor_ref)

       CASE ('_refine_ls_goodness_of_fit_ref')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%chi2)

       CASE ('_reflns_number_gt')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%reflns_number_gt)

       case ('_refine_diff_density_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%refine_diff_density_max)

       case ('_refine_diff_density_min')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%refine_diff_density_min)

       case ('_refine_diff_density_rms')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%refine_diff_density_rms)

       case ('_diffrn_measurement_device_type')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter_DEVICE%diffracto_device)
        IF(CIF_parameter_DEVICE%diffracto_device(1:6) == 'APEXII') then
         device_mark =  'Bruker-AXS'
         device_type =  'APEXII Kappa-CCD diffractometer'
        ELSEIF(CIF_parameter_DEVICE%diffracto_device(1:3) == 'D8V') then
         device_mark =  'Bruker-AXS'
         device_type =  'D8 Venture Kappa-CCD diffractometer'
        elseif(CIF_parameter_DEVICE%diffracto_device(1:3) == 'X2S') then
         device_mark = 'Bruker-AXS'
         device_type = 'SMART X2S benchtop diffractometer'
        elseif(CIF_parameter_DEVICE%diffracto_device(1:8) == 'KappaCCD') then
         device_mark = 'Nonius'
         device_type = 'KappaCCD'
        ELSEIF(CIF_parameter_DEVICE%diffracto_device(1:11) == 'CCD Saphire') then
         device_mark = 'Oxford Diffraction'
         device_type = 'Xcalibur Saphire 3'
        endif


       case ('_diffrn_ambient_temperature')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter_DEVICE%diffracto_temperature)

       case ('_diffrn_radiation_type')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter_DEVICE%diffracto_radiation_type)

       case ('_diffrn_radiation_probe')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter_DEVICE%diffracto_radiation_probe)
       case ('_diffrn_radiation_source')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter_DEVICE%diffracto_radiation_source)
       case ('_diffrn_radiation_monochromator')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter_DEVICE%diffracto_radiation_monochromator)


       case ('_diffrn_radiation_wavelength')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter_DEVICE%diffracto_radiation_wavelength)

       CASE ('_refine_ls_hydrogen_treatment')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%H_treatment)

       case ('_diffrn_reflns_theta_min')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%theta_min)

       case ('_diffrn_reflns_theta_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%theta_max)

       case ('_diffrn_reflns_limit_h_min')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%h_min)
       case ('_diffrn_reflns_limit_h_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%h_max)
       case ('_diffrn_reflns_limit_k_min')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%k_min)
       case ('_diffrn_reflns_limit_k_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%k_max)
       case ('_diffrn_reflns_limit_l_min')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%l_min)
       case ('_diffrn_reflns_limit_l_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%l_max)

       case ('_exptl_crystal_size_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%crystal_size_max)
       case ('_exptl_crystal_size_mid')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%crystal_size_mid)
       case ('_exptl_crystal_size_min')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%crystal_size_min)

       case ('_exptl_crystal_colour')
        READ(CIF_field_value, '(a)') CIF_parameter%crystal_colour

       case ('_exptl_crystal_f_000')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%F000)

       case ('_diffrn_measured_fraction_theta_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%completeness)

       case ('_exptl_absorpt_coefficient_mu')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%exptl_mu)

       case ('_exptl_absorpt_correction_type')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%absorption_correction_type)

       case ('_exptl_absorpt_correction_t_min')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%T_min)

       case ('_exptl_absorpt_correction_t_max')
        call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%T_max)

       case ('_platon_squeeze_void_nr')
        include_squeeze = .true.

       case default
   end select


  END do
 close (UNIT = CIF_unit)

  return

 end subroutine Get_CIF_parameters

!-------------------------------------------------------------------------
subroutine  create_report_header
  USE structural_report_module
  USE MACROS_module,                ONLY : u_case, Get_wingx_job,  Get_current_folder_name
  USE cryscalc_module,              ONLY : HTML_unit, text_unit, LATEX_unit, archive_CIF,  text_report, HTML_report, &
                                           LATEX_report, cryscalc, author, report_header, debug_proc
  USE IO_module


  implicit none
  CHARACTER (LEN=256)                 :: LATEX_string
  INTEGER                             :: i1
  CHARACTER (LEN=10)                  :: date, time
  character (len=256)                 :: wingx_structure_dir
  character (len=10)                  :: AUTHOR_initiales

 if(debug_proc%level_2)  call write_debug_proc_level(2, "create_report_header")

 i1 = INDEX(archive_cif, '.', back=.true.)
 if(HTML_report) then
  WRITE(HTML_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.html'
  open (UNIT=HTML_unit, FILE=TRIM(HTML_structural_report_file))
 endif

 if(text_report) then
  WRITE(TEXT_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.txt'
  open (UNIT=text_unit, FILE=TRIM(TEXT_structural_report_file))
 endif

 if(latex_report) then
  WRITE(LATEX_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.ltx'
  open (UNIT=latex_unit, FILE=TRIM(LATEX_structural_report_file))
  WRITE(PDF_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.pdf'
 endif


 if(HTML_report) then
  WRITE(HTML_unit, '(a)')  "<!-- Structural report (HTML file) created by CRYSCALC.exe"
  WRITE(HTML_unit, '(2a)')  "      Starting CIF file : ", trim(archive_CIF)
  !WRITE(HTML_unit, '(2a)') "      CRYSCALC setting file : ", trim(cryscalc%ini)
  WRITE(HTML_unit, '(a)')  ""
  WRITE(HTML_unit, '(a)')  "     Web site : http://www.cdifx.univ-rennes/cryscalc.htm"
  WRITE(HTML_unit, '(5a)') "     Version : ", trim(cryscalc%version), " [", trim(cryscalc%author), "]"
  WRITE(HTML_unit, '(a)') ""
  if(len_trim(cryscalc%report_css) ==0) then
   WRITE(HTML_unit, '(a)')  "     Default css style"
  else
   WRITE(HTML_unit, '(2a)') "     css file : ", trim(cryscalc%report_css)
  endif
  WRITE(HTML_unit, '(a)') "!-->"

  WRITE(HTML_unit, '(a)') "<HTML>"
  WRITE(HTML_unit, '(a)') "<HEAD>"
  call write_HTML_css(HTML_unit, 'report')
  WRITE(HTML_unit, '(3a)') "<TITLE>Structural report [", trim(AUTHOR%team),"]</TITLE>"
  WRITE(HTML_unit, '(a)') "</HEAD>"
  WRITE(HTML_unit, '(a)') "<BODY BGCOLOR='#FFFFFF' style='font-family:Times new roman; font-size:14; line-height:150%'>"
 endif

 if(LATEX_report) then
  call Latex_preambule
 endif


  call date_and_time(date, time)
  !call Get_Wingx_job(job, wingx_structure_dir)

 if(report_header) then

   if(AUTHOR%name /= '?' .and. AUTHOR%first_name /= '?') then
    write(AUTHOR_initiales, '(4a)') u_case(AUTHOR%first_name(1:1)), '.', u_case(AUTHOR%name(1:1)), '.'
    if(HTML_report) then
      WRITE(HTML_unit, '(a)')  "<hr>"
     WRITE(HTML_unit, '(18a)') '<p class="header">HTML report created by ', trim(AUTHOR_initiales),       &
                                            " [",trim(AUTHOR%team), "]  -  ",                         &
                                            " Date: " , date(7:8),'-',date(5:6),'-',date(1:4), '  at ', time(1:2),':',  &
                                                        time(3:4),':',time(5:6), "<br>"
    endif
    if(LATEX_report) then
     write(LATEX_unit, '(18a)') '\header{\LaTeX report created by ', trim(AUTHOR_initiales),     &
                                " [",trim(AUTHOR%team), "]  -  ",                              &
                                " Date: " , date(7:8),'-',date(5:6),'-',date(1:4), '  at ', time(1:2),':',  &
                                            time(3:4),':',time(5:6), "} \\"
    endif

   ELSE
    if(HTML_report) then
     WRITE(HTML_unit, '(a)')  "<hr>"
     WRITE(HTML_unit, '(13a)')' <p class="header"> Date: ' , date(7:8),'-',date(5:6),'-',date(1:4), '  at ', time(1:2),':', &
                              time(3:4),':',time(5:6), "<br>"
    endif
    if(LATEX_report) then
     WRITE(LATEX_unit, '(13a)')'\header{Date: ' , date(7:8),'-',date(5:6),'-',date(1:4), '  at ', time(1:2),':', &
                                                  time(3:4),':',time(5:6), "} \\"
    endif
   ENDIF

   !if(wingx_structure_dir =='?') call Get_current_folder_name(wingx_structure_dir)
   call Get_current_folder_name(wingx_structure_dir)

   !i1 = index(wingx_structure_dir, '_')
   !if(i1 /=0 .and. i1 > 1) then
   ! write(job, '(a)') wingx_structure_dir(1:i1-1)
   !endif
  if(HTML_report) then
   WRITE(HTML_unit, '(3a)') 'Working directory: <font face="courier">', trim(wingx_structure_dir), '</font><br>'
   WRITE(HTML_unit, '(3a)') "Input CIF file : ", trim(archive_CIF), '<br>'
   WRITE(HTML_unit, '(5a)') "CRYSCALC version : ", trim(cryscalc%version), " [", trim(cryscalc%author), "]"
   WRITE(HTML_unit, '(a)')  "</p class='header'>"
   WRITE(HTML_unit, '(a)')  "<hr>"
  endif
  if(LATEX_report) then
   LATEX_string = wingx_structure_dir
   call check_LATEX_file_name(LATEX_string)
   WRITE(LATEX_unit, '(3a)') '\header{Working directory: ', trim(LATEX_string), '} \\'
   LATEX_string = archive_cif
   call check_LATEX_file_name(LATEX_string)
   WRITE(LATEX_unit, '(3a)') "\header{Input CIF file : ", trim(LATEX_string), '} \\'
   WRITE(LATEX_unit, '(5a)') "\header{CRYSCALC version : ", trim(cryscalc%version), " [", trim(cryscalc%author), "]} \\"
   write(LATEX_unit, '(a)')  "\rule{\linewidth}{0.5pt}"
   write(LATEX_unit, '(a)')  "\vspace{1.0cm}"
   WRITE(LATEX_unit, '(a)')  ""
  endif
 endif

 return
end subroutine create_report_header

!----------------------------------------------------------------------------------------------------------------------------------
subroutine  create_report
 USE structural_report_module
 USE MACROS_module,                ONLY : test_file_exist, remove_car, nombre_de_colonnes, l_case, u_case,             &
                                          check_character, Get_wingx_job,  Get_current_folder_name, replace_car, multiple, &
                                          replace_car2
 USE cryscalc_module,              ONLY : CIF_unit, HTML_unit, text_unit, LATEX_unit, archive_CIF, long_report,        &
                                          CIF_torsion_limit, text_report, HTML_report, LATEX_report, report_logo,      &
                                          SQUEEZE, my_browser, my_pdflatex, cryscalc, message_text, debug_proc,        &
                                          cryscalc, structure_solution, structure_refinement, author,                  &
                                          report_header, include_SQUEEZE, debug_proc
 USE CIF_module,                   ONLY : CIF_parameter, CIF_parameter_DEVICE, CIF_CELL_measurement
 USE external_applications_module, ONLY : launch_browser, launch_word, launch_pdflatex
 USE IO_module
 USE LATEX_module,                 ONLY : logo
 USE Special_characters_module


 implicit none
 LOGICAL                             :: file_exist
 CHARACTER (LEN=256)                 :: CIF_input_line, CIF_field, output_line
 CHARACTER (LEN=512)                 :: HTML_string
 CHARACTER (LEN=256), dimension(4)   :: report_string
 CHARACTER (LEN=512)                 :: LATEX_string
 CHARACTER (LEN=256)                 :: DOS_command
 CHARACTER (LEN=256)                 :: solving_method
 CHARACTER (LEN=32)                  :: GIF_file, PNG_file, JPG_file
 CHARACTER (LEN=32)                  :: IMG_file
 CHARACTER (LEN=32)                  :: job_latex
 LOGICAL                             :: PNG_type, JPG_type
 INTEGER                             :: i_error, long_field, long, i, i1, i2
 integer                             :: n, op_numor

 integer,   dimension(500)           :: op_n       ! numero de l'op. de sym.
 integer,   dimension(500,3)         :: op_t       ! partie translation
 integer,   dimension(3)             :: t
 integer                             :: n_sym_htab, n_sym_dist, n_sym_ang, n_sym_tor
 integer                             :: nb_col

 CHARACTER (LEN=12), dimension(500) :: site_sym
 CHARACTER (len=12), dimension(500) :: dico
 INTEGER           , dimension(500) :: num_site_sym

 CHARACTER (LEN=12), dimension(2)    :: dist_atom
 CHARACTER (LEN=12)                  :: dist_value
 CHARACTER (LEN=12)                  :: dist_sym
 CHARACTER (LEN=12), dimension(3)    :: ang_atom
 CHARACTER (LEN=12)                  :: ang_value
 real                                :: ang_value_real
 CHARACTER (LEN=12)                  :: ang_esd_string
 CHARACTER (LEN=12), dimension(2)    :: ang_sym
 CHARACTER (LEN=12), dimension(4)    :: torsion_ang_atom
 CHARACTER (LEN=12)                  :: torsion_ang_value
 CHARACTER (LEN=12), dimension(4)    :: torsion_sym
 real                                :: torsion_value_real
 CHARACTER (LEN=12)                  :: torsion_esd_string
 INTEGER                             :: torsion_esd


 CHARACTER (LEN=12)                  :: atom_label, atom_typ, atom_x, atom_y, atom_z, atom_Ueq, atom_adp_type, atom_occ
 CHARACTER (LEN=12)                  :: atom_U11, atom_U22, atom_U33, atom_U23, atom_U13, atom_U12
 CHARACTER (len=12)                  :: site_label_D, site_label_H, site_label_A
 CHARACTER (len=12)                  :: dist_DH, dist_HA, dist_DA, angle_DHA, site_sym_A
 LOGICAL                             :: alpha_car, num_car, feder
 LOGICAL                             :: SHELXL_2014, created_by_CRYSCALC


 if(debug_proc%level_2)  call write_debug_proc_level(2, "create_report")

 n_sym_htab = 0
 n_sym_dist = 0
 n_sym_ang  = 0
 n_sym_tor  = 0



 if(HTML_report) then
  WRITE(HTML_unit, '(a)')  "<br>"
  WRITE(HTML_unit, '(a)')  "<div class='cadre'>"
  !WRITE(HTML_unit, '(3a)') "<p class='title_main'><i>", trim(job),"</i>: Crystal structure report</p><br>"
  WRITE(HTML_unit, '(3a)') "<p class='title_main'><font face='courier'>", trim(job),"</font>: Crystal structure report</p><br>"
  WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;X-ray crystallographic study</p>"
 endif

 if(LATEX_report) then
  write(LATEX_unit, '(a)') '\footnotesize'
  WRITE(LATEX_unit, '(a)')  ""
  !job = replace_car(job, '_', '/')
  !job = replace_car(job, "/", "\_")
  job_latex = replace_car(job, '_', '\_')
  WRITE(LATEX_unit, '(4a)')  '\titre{', trim(job_latex), '}{: Crystal structure report','}'
  WRITE(LATEX_unit, '(a)')   '\SousTitre{X-ray crystallographic study}'
 endif


  if(CIF_parameter%formula_moiety == '?' .or. len_trim(CIF_parameter%formula_moiety) ==0)  &
     CIF_parameter%formula_moiety = CIF_parameter%formula_sum

  !----------
  long = len_trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
  if(CIF_parameter_DEVICE%diffrn_measurement_device_type(1:1)== "'"           .and. &
     CIF_parameter_DEVICE%diffrn_measurement_device_type(long:long)== "'") then
    if(CIF_parameter_DEVICE%diffrn_measurement_device_type(2:11) == "D8 VENTURE") then
     feder = .true.
    else
     feder = .false.
    end if
   else
    WRITE(HTML_unit, '(4a)') TRIM(CIF_parameter_DEVICE%diffrn_measurement_device_type), " diffractometer,"
    IF(CIF_parameter_DEVICE%diffrn_measurement_device_type(1:10) == "D8 VENTURE")  then
     feder = .true.
    else
     feder = .false.
    end if
   endif
  !--------


  ! oct. 2016
  ! ----------- insertion de l'image du cristal dans .HTML et .Latex ----------------
  ! recherche de la presence de fichier.PNG/.JPG
  if(.not. TEXT_report) then
   file_exist = .false.
   JPG_type = .false.
   PNG_type = .false.

   do
    write(PNG_file, '(2a)') trim(job), ".png"
    call test_file_exist(trim(PNG_file), file_exist, 'no_out')
    if(file_exist) then
    ! janvier 2017 : conversion fichier.png en fichier.jpg (plus petit) avec Convert (ImageMagick)
     write(JPG_file, '(2a)') trim(job), ".jpg"
     write(DOS_command, '(4a)') 'convert -bordercolor rgb(255,255,255) -border 1 ', trim(PNG_file), ' ', trim(JPG_file)
     call system(trim(DOS_command))
     call test_file_exist(trim(JPG_file), file_exist, 'no_out')
     if(file_exist) then
      if(.not. LATEX_report) then
       write(IMG_file, '(2a)') trim(JPG_file)
       JPG_type = .true.
      else      ! le .JPG cree par CONVERT ne convient pas pour PDFLATEX
       write(IMG_file, '(2a)') trim(PNG_file)
       PNG_type = .true.
      end if
      !call system("del "//trim(PNG_file))
     else
      write(IMG_file, '(2a)') trim(PNG_file)
      PNG_type = .true.
     end if
     exit
    end if

    write(IMG_file, '(2a)') trim(job), ".jpg"
    call test_file_exist(trim(IMG_file), file_exist, 'no_out')
    if(file_exist) then
     JPG_type = .true.
     exit
    end if

    exit
   end do



   if(file_exist) then
    if(HTML_report) then
     WRITE(HTML_unit, '(a)') ""
     WRITE(HTML_unit, '(a)') "<br>"
     WRITE(HTML_unit, '(3a)') '<center><img width=300 src="', trim(IMG_file),'"></center>'
    endif

    if(LATEX_report) then
     WRITE(LATEX_unit, '(a)') ""
     write(LATEX_unit, '(a)')  '\begin{figure}[h]'
     write(LATEX_unit, '(a)')  ' \centering'
     write(LATEX_unit, '(3a)') ' \includegraphics[width=200.pt]{',trim(IMG_file),'}'
     write(LATEX_unit, '(a)')  '\end{figure}'
    endif
   endif
  end if ! fin de la condition if (.not. TXT_report)

  ! ---------------------------------------------


  if(index(CIF_parameter_DEVICE%diffrn_measurement_device_type, 'X2S') /=0) then
   if(HTML_report .or. LATEX_report) then
    call report_crystal_study("X2S")
   endif
  else


  if(HTML_report) then
   CALL transf_moiety_string("HTML", CIF_parameter%formula_moiety , HTML_string)
   WRITE(HTML_unit, '(a)') "<p class='retrait_1'>"
   WRITE(HTML_unit, '(a)')  ""
   !WRITE(HTML_unit, '(5a)', advance='NO') "(",TRIM(CIF_parameter%formula_moiety),"); M = ",TRIM(CIF_parameter%formula_weight),"."
   WRITE(HTML_unit, '(5a)', advance='NO') "(",TRIM(HTML_string),"); <i>M</i> = ",TRIM(CIF_parameter%formula_weight),".&nbsp;"

   long = len_trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
   if(CIF_parameter_DEVICE%diffrn_measurement_device_type(1:1)== "'"           .and. &
      CIF_parameter_DEVICE%diffrn_measurement_device_type(long:long)== "'") then
    if(feder) then
     WRITE(HTML_unit, '(4a)') CIF_parameter_DEVICE%diffrn_measurement_device_type(2:long-1), " diffractometer  [*],&nbsp;"
    else
     WRITE(HTML_unit, '(4a)') CIF_parameter_DEVICE%diffrn_measurement_device_type(2:long-1), " diffractometer,&nbsp;"
    end if
   else
    IF(feder)  then
     WRITE(HTML_unit, '(4a)') TRIM(CIF_parameter_DEVICE%diffrn_measurement_device_type), " diffractometer [*],&nbsp;"
    else
     WRITE(HTML_unit, '(4a)') TRIM(CIF_parameter_DEVICE%diffrn_measurement_device_type), " diffractometer,&nbsp;"
    end if
   endif

   !long = len_trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
   !if(CIF_parameter_DEVICE%diffrn_measurement_device_type(1:1)== "'"           .and. &
   !   CIF_parameter_DEVICE%diffrn_measurement_device_type(long:long)== "'") then
   ! WRITE(HTML_unit, '(4a)') CIF_parameter_DEVICE%diffrn_measurement_device_type(2:long-1), " diffractometer"
   ! if(CIF_parameter_DEVICE%diffrn_measurement_device_type(2:11) == "D8 VENTURE") then
   !  write(HTML_unit, '(a)', advance='NO') " [*],&nbsp;"
   !  feder = .true.
   ! else
   !  write(HTML_unit, '(a)', advance='NO') ",&nbsp;"
   !  feder = .false.
   ! end if
   !else
   ! WRITE(HTML_unit, '(4a)') TRIM(CIF_parameter_DEVICE%diffrn_measurement_device_type), " diffractometer,"
   ! IF(CIF_parameter_DEVICE%diffrn_measurement_device_type(1:10) == "D8 VENTURE")  then
   !  write(HTML_unit, '(a)', advance='NO') " [*], "
   !  feder = .true.
   ! else
   !  write(HTML_unit, '(a)', advance='NO') " ,&nbsp;"
   !  feder = .false.
   ! end if
   !endif
  endif


  if(LATEX_report) then
   CALL transf_moiety_string("LATEX", CIF_parameter%formula_moiety , LATEX_string)
   !LATEX_string = CIF_parameter%formula_moiety
   WRITE(LATEX_unit, '(a)')  ""
   WRITE(LATEX_unit, '(a)')  "{\setlength{\baselineskip}{1.2\baselineskip}"
   WRITE(LATEX_unit, '(5a)', advance='NO') "\noindent ($",TRIM(LATEX_string),"$); $M = ",TRIM(CIF_parameter%formula_weight),"$. "

   long = len_trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
   if(CIF_parameter_DEVICE%diffrn_measurement_device_type(1:1)== "'"           .and. &
      CIF_parameter_DEVICE%diffrn_measurement_device_type(long:long)== "'") then
    WRITE(LATEX_unit, '(4a)') CIF_parameter_DEVICE%diffrn_measurement_device_type(2:long-1), " diffractometer"
    if(feder) then
     write(LATEX_unit, '(a)') " [*],~"
    else
     write(LATEX_unit, '(a)') " ,~"
    end if
   else
    WRITE(LATEX_unit, '(4a)') TRIM(CIF_parameter_DEVICE%diffrn_measurement_device_type), " diffractometer"
    IF(feder)  then
     write(LATEX_unit, '(a)', advance='NO') " [*], "
    else
     write(LATEX_unit, '(a)', advance='NO') " ,;"
    end if
   endif


   !long = len_trim(CIF_parameter_DEVICE%diffrn_measurement_device_type)
   !if(CIF_parameter_DEVICE%diffrn_measurement_device_type(1:1)== "'"           .and. &
   !   CIF_parameter_DEVICE%diffrn_measurement_device_type(long:long)== "'") then
   ! WRITE(LATEX_unit, '(4a)') CIF_parameter_DEVICE%diffrn_measurement_device_type(2:long-1), " diffractometer"
   ! if(CIF_parameter_DEVICE%diffrn_measurement_device_type(2:11) == "D8 VENTURE") then
   !  write(LATEX_unit, '(a)') " [*],~"
   !  feder = .true.
   ! else
   !  write(LATEX_unit, '(a)') " ,~"
   !  feder = .false.
   ! end if
   !else
   ! WRITE(LATEX_unit, '(4a)') TRIM(CIF_parameter_DEVICE%diffrn_measurement_device_type), " diffractometer"
   ! IF(CIF_parameter_DEVICE%diffrn_measurement_device_type(1:10) == "D8 VENTURE")  then
   !  write(LATEX_unit, '(a)', advance='NO') " [*], "
   !  feder = .true.
   ! else
   !  write(LATEX_unit, '(a)', advance='NO') " ,;"
   !  feder = .false.
   ! end if
   !endif
  endif


  IF(CIF_parameter_DEVICE%diffracto_radiation_type(1:5) == 'MoK\a') then
   if(HTML_report) then
    WRITE(HTML_unit, '(a)', advance='NO') "Mo-K&alpha; radiation "
   endif
   if(LATEX_report) then
    WRITE(LATEX_unit, '(a)', advance='NO') "$Mo-K_{\alpha}$ radiation "
   endif
  ELSEIF(CIF_parameter_DEVICE%diffracto_radiation_type(1:5) == 'CuK\a') then
   if(HTML_report) then
    WRITE(HTML_unit, '(a)', advance='NO') "Cu-K&alpha; radiation "
   endif
   if(LATEX_report) then
    WRITE(LATEX_unit, '(a)', advance='NO') "$Cu-K_{\alpha}$ radiation "
   endif

  else
   if(HTML_report) then
    WRITE(HTML_unit, '(a)', advance='NO') TRIM(CIF_parameter_DEVICE%diffracto_radiation_type), " radiation"
   endif
   if(LATEX_report) then
    WRITE(LATEX_unit, '(a)', advance='NO') TRIM(CIF_parameter_DEVICE%diffracto_radiation_type), " radiation"
   endif
  endif
  if(HTML_report) then
   WRITE(HTML_unit, '(3a)') "(&lambda; = ", TRIM(CIF_parameter_DEVICE%diffracto_radiation_wavelength)," Å),"
   WRITE(HTML_unit, '(3a)', advance='NO') "<i>T</i> = ", TRIM(CIF_parameter_DEVICE%diffracto_temperature), " K; "
  endif
  if(LATEX_report) then
   WRITE(LATEX_unit, '(3a)') "($\lambda$ = ", TRIM(CIF_parameter_DEVICE%diffracto_radiation_wavelength)," \AA),"
   WRITE(LATEX_unit, '(3a)', advance='NO') "$T = ", TRIM(CIF_parameter_DEVICE%diffracto_temperature), " K$; "
  endif

  alpha_car    = .false.
  num_car      = .false.
  HTML_string  = ''
  LATEX_string = ''
  long = len_trim(HTML_string)

  if(HTML_report) then
   HTML_string = CIF_parameter%symmetry_space_group
   call check_letter_string(HTML_string, .false.)
   call check_indice_string(HTML_string)
   !WRITE(HTML_unit, '(6a)') TRIM(CIF_parameter%symmetry_cell_setting), " <i>", trim(HTML_string),  &
   !                        "</i> (I.T.#", TRIM(CIF_parameter%symmetry_IT_number), "), "
   WRITE(HTML_unit, '(6a)') TRIM(CIF_parameter%symmetry_cell_setting), " ", trim(HTML_string),  &
                           " (I.T.#", TRIM(CIF_parameter%symmetry_IT_number), "), "
  endif

  if(LATEX_report) then
   LATEX_string = CIF_parameter%symmetry_space_group
   call check_letter_string(LATEX_string, .false.)
   call check_indice_string(LATEX_string)
   WRITE(LATEX_unit, '(6a)') TRIM(CIF_parameter%symmetry_cell_setting), " $", trim(LATEX_string),  &
                           "$ (I.T.\#", TRIM(CIF_parameter%symmetry_IT_number), "), "
  endif


  if(HTML_report) then
   call report_cell_parameters("HTML")
   WRITE(HTML_unit, '(3a)', advance='NO') "&nbsp;<i>Z</i> = ",TRIM(CIF_parameter%formula_units_Z), ", "
   WRITE(HTML_unit, '(3a)') "<i>d</i> = ",TRIM(CIF_parameter%exptl_density), " g.cm<sup>-3</sup>, "
   WRITE(HTML_unit, '(3a)') "&mu; = ",TRIM(CIF_parameter%exptl_mu), " mm<sup>-1</sup>. "
  endif

  if(LATEX_report) then
   call report_cell_parameters("LATEX")
   WRITE(LATEX_unit, '(3a)', advance='NO') "~$Z = ",TRIM(CIF_parameter%formula_units_Z), "$, "
   WRITE(LATEX_unit, '(3a)') "$d = ",TRIM(CIF_parameter%exptl_density), "~g.cm^{-3}$, "
   WRITE(LATEX_unit, '(3a)') "$\mu = ",TRIM(CIF_parameter%exptl_mu), "~mm^{-1}$. "
  endif


  long = len_trim(structure_solution%name)
  if(long == 6) then
   if(u_case(structure_solution%name) == 'SHELXT') then
    solving_method = 'dual-space algorithm'
   endif
  else
   solving_method = 'direct methods'
  end if

  SHELXL_2014 = .false.
  long = len_trim(structure_refinement%name)
  if(long == 11) then
   if(structure_refinement%name(1:11) == 'SHELXL-2014' .or. &
      structure_refinement%name(1:11) == 'SHELXL-2016') SHELXL_2014 = .true.
  endif


  if(HTML_report) then
   !WRITE(HTML_unit, '(a)') "The structure was solved by direct methods using the SIR97 program [1], "
   WRITE(HTML_unit, '(4a)') "The structure was solved by ", trim(solving_method), " using the <i>",  &
                             u_case(trim(structure_solution%name)), "</i> program [1], "
   WRITE(HTML_unit, '(a)')  "and then refined with full-matrix least-square methods based on <i>F</i><sup>2</sup>"
   if(CIF_parameter%WinGX_used) then
    WRITE(HTML_unit, '(4a)') "(<i>", u_case(trim(structure_refinement%name)), "</i>) [2] ", &
                             "with the aid of the <i>WINGX</i> [3] program."
    if(include_SQUEEZE) then
      if(SHELXL_2014) then
      WRITE(HTML_unit, '(3a)')  "The contribution of the disordered solvents to the structure factors was ", &
                                "calculated by the <i>PLATON</i> SQUEEZE procedure [4] and then taken into account ", &
                                "in the final <i>SHELXL-2014</i> least-square refinement."
     else
      WRITE(HTML_unit, '(4a)') "The contribution of the disordered solvents to the calculated structure factors was ",    &
                               "estimated following the <i>BYPASS</i> algorithm [4], implemented as the <i>SQUEEZE</i> ", &
                               "option in <i>PLATON</i> [5]. A new data set, free of solvent contribution, was then ",    &
                               "used in the final refinement."
     end if
    end if
   else
    WRITE(HTML_unit, '(3a)') "(<i>", u_case(trim(structure_refinement%name)), "</i>) [2]."
    if(include_SQUEEZE) then
      if(SHELXL_2014) then
      WRITE(HTML_unit, '(3a)')  "The contribution of the disordered solvents to the structure factors was ", &
                                "calculated by the <i>PLATON</i> SQUEEZE procedure [3] and then taken into account ", &
                                "in the final <I>SHELXL-2014</i> least-square refinement."
     else
      WRITE(HTML_unit, '(4a)') "The contribution of the disordered solvents to the calculated structure factors was ",    &
                               "estimated following the <i>BYPASS</i> algorithm [3], implemented as the <i>SQUEEZE</i> ", &
                               "option in <i>PLATON</i> [4]. A new data set, free of solvent contribution, was then ",    &
                               "used in the final refinement."
     end if
    end if

   end if


   !if(SQUEEZE%procedure) then

   WRITE(HTML_unit, '(a)') "All non-hydrogen atoms were refined with anisotropic atomic displacement parameters. "
   IF(CIF_parameter%H_treatment(1:5) == 'mixed') then
    WRITE(HTML_unit, '(2a)') "Except <span style='color:red'><b>XXX</b></span> linked hydrogen atoms that were introduced ", &
                             "in the structural model through Fourier difference maps analysis,"
   endif
   WRITE(HTML_unit, '(2a)') "H atoms were finally included in their calculated positions. A final refinement on ", &
                            "<i>F</i><sup>2</sup> with "
   WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%reflns_number_total), " unique intensities and "
   WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%refine_ls_number_parameters), &
                            " parameters converged at &omega;<i>R</i>(<i>F</i><sup>2</sup>) = "
   WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%refine_ls_wR_factor_gt), " (<i>R</i>(<i>F</i>) = "
   WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%refine_ls_R_factor_gt), ") for "
   WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%reflns_number_gt), " observed reflections with <i>I</i> > 2&sigma;(<i>I</i>)."
   WRITE(HTML_unit, '(a)') "</p>"
  endif

   if(LATEX_report) then
    WRITE(LATEX_unit, '(5a)') "The structure was solved by ", trim(solving_method), " using the ",  &
                              u_case(trim(structure_solution%name)), " program [1], "
    if(CIF_parameter%WinGX_used) then
     WRITE(LATEX_unit, '(3a)') "and then refined with full-matrix least-square methods based on $F^2$ (", &
                               u_case(trim(structure_refinement%name)), ") [2] with the aid of the WINGX [3] program."
     if(include_SQUEEZE) then
       if(SHELXL_2014) then
       WRITE(LATEX_unit, '(3a)') "The contribution of the disordered solvents to the structure factors was ", &
                                 "calculated by the \textit{PLATON} SQUEEZE procedure [4] and then taken into account", &
                                 " in the final \textit{SHELXL-2014} least-square refinement."
      else
       WRITE(LATEX_unit, '(4a)') "The contribution of the disordered solvents to the calculated structure factors was ", &
                                 "estimated following the \textit{BYPASS} algorithm [4], implemented as the \textit{SQUEEZE} ", &
                                 "option in \textit{PLATON} [5]. A new data set, free of solvent contribution, was then ",    &
                                 "used in the final refinement."
      end if
     end if


    else
     WRITE(LATEX_unit, '(3a)') "and then refined with full-matrix least-square methods based on $F^2$ (", &
                               u_case(trim(structure_refinement%name)), " program [2])."
     if(include_SQUEEZE) then
      if(SHELXL_2014) then
       WRITE(LATEX_unit, '(3a)') "The contribution of the disordered solvents to the structure factors was ", &
                                 "calculated by the \textit{PLATON} SQUEEZE procedure [3] and then taken into account", &
                                 " in the final \textit{SHELXL-2014} least-square refinement."
      else
       WRITE(LATEX_unit, '(4a)') "The contribution of the disordered solvents to the calculated structure factors was ", &
                                 "estimated following the \textit{BYPASS} algorithm [3], implemented as the \textit{SQUEEZE} ", &
                                 "option in \textit{PLATON} [4]. A new data set, free of solvent contribution, was then ",    &
                                 "used in the final refinement."
      end if
     endif

    endif


    WRITE(LATEX_unit, '(a)') "All non-hydrogen atoms were refined with anisotropic atomic displacement parameters. "
    IF(CIF_parameter%H_treatment(1:5) == 'mixed') then
     WRITE(LATEX_unit, '(2a)') "Except \textcolor{red_25b}{XXX} linked hydrogen atoms that were introduced ", &
                               "in the structural model through Fourier difference maps analysis,"
    endif
    WRITE(LATEX_unit, '(4a)') "H atoms were finally included in their calculated positions. A final refinement on ", &
                              "$F^2$ with ", TRIM(CIF_parameter%reflns_number_total), " unique intensities and "
    WRITE(LATEX_unit, '(2a)') TRIM(CIF_parameter%refine_ls_number_parameters), &
                            " parameters converged at $\omega R(F^2) =~"
    WRITE(LATEX_unit, '(3a)') TRIM(CIF_parameter%refine_ls_wR_factor_gt), "~(RF =~", TRIM(CIF_parameter%refine_ls_R_factor_gt)
    WRITE(LATEX_unit, '(3a)') ")~$for ", TRIM(CIF_parameter%reflns_number_gt), " observed reflections with ($I > 2\sigma $)."
    WRITE(LATEX_unit, '(a)') ""
    WRITE(LATEX_unit, '(a)') "\par}"
   endif



  if(HTML_report) then
   WRITE(HTML_unit, '(a)') "<p class='retrait_1'>"
   WRITE(HTML_unit, '(3a)') "[1] &nbsp; ", trim(structure_solution%reference),   "<br>"
   WRITE(HTML_unit, '(3a)') "[2] &nbsp; ", trim(structure_refinement%reference), "<br>"
   if(CIF_parameter%WinGX_used) then
    WRITE(HTML_unit, '(a)') "[3] &nbsp; L. J. Farrugia, J. Appl. Cryst., 2012, 45, 849-854<br>"
    if(include_SQUEEZE) then
     if(SHELXL_2014) then
      WRITE(HTML_unit, '(a)') "[4] &nbsp; A.L. Spek, Acta Cryst. (2015) C71, 9-18<br>"
     else
      WRITE(HTML_unit, '(a)') "[4] &nbsp; P. v.d. Sluis and A.L. Spek, Acta Cryst. (1990) A46, 194-201<br>"
      WRITE(HTML_unit, '(a)') "[5] &nbsp; A. L. Spek, J. Appl. Cryst. (2003), 36, 7-13<br>"
     end if
    end if
   else
    if(include_SQUEEZE) then
     if(SHELXL_2014) then
      WRITE(HTML_unit, '(a)') "[3] &nbsp; A. L. Spek, Acta Cryst. (2015), C71, 9-18<br>"
     else
      WRITE(HTML_unit, '(a)') "[3] &nbsp; P. v.d. Sluis and A.L. Spek, Acta Cryst. (1990) A46, 194-201<br>"
      WRITE(HTML_unit, '(a)') "[4] &nbsp; A. L. Spek, J. Appl. Cryst. (2003), 36, 7-13<br>"
     end if
    end if
   end if

   IF(feder) WRITE(HTML_unit, '(a)') '<br><span class="founds">[*] Thanks to FEDER founds</span>'
   WRITE(HTML_unit, '(a)') "</p>"
   WRITE(HTML_unit, '(a)') "<br>"
  end if

  if(LATEX_report) then
   WRITE(LATEX_unit, '(a)') ""
   WRITE(LATEX_unit, '(a)') "\begin{enumerate}"
   WRITE(LATEX_unit, '(2a)') "\item ", trim(structure_solution%reference)
   WRITE(LATEX_unit, '(2a)') "\item ", trim(structure_refinement%reference)
   !WRITE(LATEX_unit, '(a)')  "\item  L. J. Farrugia, J. Appl. Cryst., 1999, 32, 837-838"
   if(CIF_parameter%WinGX_used) then
   WRITE(LATEX_unit, '(a)')  "\item   L. J. Farrugia, J. Appl. Cryst., 2012, 45, 849-854"
   end if
   if(SQUEEZE%procedure) then
    if(SHELXL_2014) then
     WRITE(LATEX_unit, '(a)') "\item A.L. Spek, Acta Cryst. (2015) C71, 9-18"
    else
     WRITE(LATEX_unit, '(a)') "\item P. v.d. Sluis and A.L. Spek, Acta Cryst. (1990) A46, 194-201"
     WRITE(LATEX_unit, '(a)') "\item A. L. Spek, J. Appl. Cryst. (2003), 36, 7-13"!
    end if
   end if
   WRITE(LATEX_unit, '(a)') "\end{enumerate}"

   IF(feder) then
    WRITE(LATEX_unit, '(a)') "\vspace{0.5cm}"
    WRITE(LATEX_unit, '(a)') "\header{[*] Thanks to the FEDER founds}"
    write(LATEX_unit, '(a)') '\footnotesize'
   end if
   WRITE(LATEX_unit, '(a)') ""

  end if
  end if   ! fin de la condition sur le type de diffracto. utilise

  long = len_trim(CIF_parameter%crystal_colour)
  if(CIF_parameter%crystal_colour(1:1) == '"' .and. CIF_parameter%crystal_colour(long:long) == '"') then
   CIF_parameter%crystal_colour = CIF_parameter%crystal_colour(2:long-1)
  end if
  if(CIF_parameter%crystal_colour(1:1) == "'" .and. CIF_parameter%crystal_colour(long:long) == "'") then
   CIF_parameter%crystal_colour = CIF_parameter%crystal_colour(2:long-1)
  end if


  if(HTML_report) then

   WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Structural data</p>"

   WRITE(HTML_unit, '(a)')  "<pre style='font-size:14' class='color_grey1'>"
   CALL transf_moiety_string("HTML", CIF_parameter%formula_sum , HTML_string)
   WRITE(HTML_unit, '(10x,2a)') "Empirical formula                      ", TRIM(adjustl(HTML_string))
   if(CIF_parameter%formula_sum(1:len_trim(CIF_parameter%formula_sum)) /=   &
      CIF_parameter%formula_moiety(1:len_trim(CIF_parameter%formula_moiety)))  then
    CALL transf_moiety_string("HTML", CIF_parameter%formula_moiety , HTML_string)
    WRITE(HTML_unit, '(10x,2a)') "Extended  formula                      ", TRIM(adjustl(HTML_string))
   end if
   WRITE(HTML_unit, '(10x,2a)') "Formula weight                         ", TRIM(CIF_parameter%formula_weight)
   WRITE(HTML_unit, '(10x,3a)') "Temperature                            ", &
                                 TRIM(CIF_parameter_DEVICE%diffracto_temperature)," <i>K</i>"
   WRITE(HTML_unit, '(10x,3a)') "Wavelength                             ", &
                                 TRIM(CIF_parameter_DEVICE%diffracto_radiation_wavelength), " Å "
   HTML_string = CIF_parameter%symmetry_space_group
   call check_letter_string(HTML_string, .false.)
   call check_indice_string(HTML_string)

   !WRITE(HTML_unit, '(10x,6a)') "Crystal system, space group            ", TRIM(CIF_parameter%symmetry_cell_setting), ", ", &
   !                             "<i>",TRIM(HTML_string),"</i>"
   WRITE(HTML_unit, '(10x,4a)') "Crystal system, space group            ", TRIM(CIF_parameter%symmetry_cell_setting), ", ", &
                                TRIM(HTML_string)
   WRITE(HTML_unit, '(10x,6a)') "Unit cell dimensions                   ", "a = ", TRIM(CIF_parameter%cell_length_a), &
                                " Å, &alpha; = ", TRIM(CIF_parameter%cell_angle_alpha), " °"
   WRITE(HTML_unit, '(10x,6a)') "                                       ", "b = ", TRIM(CIF_parameter%cell_length_b), &
                                " Å, &beta; = ", TRIM(CIF_parameter%cell_angle_beta),  " °"
   WRITE(HTML_unit, '(10x,6a)') "                                       ", "c = ", TRIM(CIF_parameter%cell_length_c), &
                                " Å, &gamma; = ", TRIM(CIF_parameter%cell_angle_gamma), " °"
   WRITE(HTML_unit, '(10x,3a)') "Volume                                 ", TRIM(CIF_parameter%cell_volume), " Å<sup>3</sup>"
   WRITE(HTML_unit, '(10x,5a)') "<i>Z</i>, Calculated density                  ", TRIM(CIF_parameter%formula_units_Z), ', ', &
                                TRIM(CIF_parameter%exptl_density), " (g.cm<sup>-3</sup>)"
   WRITE(HTML_unit, '(10x,3a)') "Absorption coefficient                 ", TRIM(CIF_parameter%exptl_mu), " mm<sup>-1</sup>"
   IF(long_report) then
    WRITE(HTML_unit, '(10x,2a)') "<i>F</i>(000)                                 ", TRIM(CIF_parameter%F000)
    WRITE(HTML_unit, '(10x,7a)') "Crystal size                           ", TRIM(CIF_parameter%crystal_size_max), " x ", &
                               TRIM(CIF_parameter%crystal_size_mid), " x ", TRIM(CIF_parameter%crystal_size_min), " mm"
    WRITE(HTML_unit, '(10x,2a)') "Crystal color                          ", TRIM(CIF_parameter%crystal_colour)
    WRITE(HTML_unit, '(10x,5a)') "Theta range for data collection        ", TRIM(CIF_parameter%theta_min), " to ", &
                                TRIM(CIF_parameter%theta_max), " °"
    WRITE(HTML_unit, '(10x,4a)') "h_min, h_max                           ", TRIM(CIF_parameter%h_min), ", ", &
                                TRIM(CIF_parameter%h_max)
    WRITE(HTML_unit, '(10x,4a)') "k_min, k_max                           ", TRIM(CIF_parameter%k_min), ", ", &
                                TRIM(CIF_parameter%k_max)
    WRITE(HTML_unit, '(10x,4a)') "l_min, l_max                           ", TRIM(CIF_parameter%l_min), ", ", &
                                TRIM(CIF_parameter%l_max)
    WRITE(HTML_unit, '(10x,7a)') "Reflections collected / unique         ", TRIM(CIF_parameter%diffrn_reflns_number), &
                               " / ", TRIM(CIF_parameter%reflns_number_total), " [R(int)<sup>a</sup> = ",             &
                               TRIM(CIF_parameter%diffrn_reflns_av_R_equivalents), "]"
    WRITE(HTML_unit, '(10x,2a)') "Reflections [I>2&sigma;]                     ", CIF_parameter%reflns_number_gt
    WRITE(HTML_unit, '(10x,2a)') "Completeness to theta_max              ", TRIM(CIF_parameter%completeness)
    WRITE(HTML_unit, '(10x,2a)') "Absorption correction type             ", TRIM(CIF_parameter%absorption_correction_type)
    if(CIF_parameter%absorption_correction_type(1:4) /= 'none') then
     WRITE(HTML_unit, '(10x,4a)') "Max. and min. transmission             ", TRIM(CIF_parameter%T_max), " , ", &
                               TRIM(CIF_parameter%T_min)
    endif
    WRITE(HTML_unit, '(10x,3a)') "Refinement method                      ", "Full-matrix least-squares on <i>F</i><sup>2</sup>"
   endif
   IF(long_report) then
    WRITE(HTML_unit, '(10x,6a)') "Data / restraints / parameters         ", TRIM(CIF_parameter%reflns_number_total), " / ", &
                               TRIM(CIF_parameter%restraints_number), " / ", TRIM(CIF_parameter%refine_ls_number_parameters)
    if(CIF_parameter%abs_structure_Flack(1:1) /= '?') then
    WRITE(HTML_unit, '(10x,2a)') "Flack parameter                        ", TRIM(CIF_parameter%abs_structure_Flack)
    end if
    WRITE(HTML_unit, '(10x,2a)') "<sup>b</sup>S (Goodness-of-fit)                    ", TRIM(CIF_parameter%Chi2)
   else
    WRITE(HTML_unit, '(10x,2a)') "Unique reflections                     ", TRIM(CIF_parameter%reflns_number_total)
    WRITE(HTML_unit, '(10x,2a)') "Refined parameters                     ", TRIM(CIF_parameter%refine_ls_number_parameters)
   endif
   WRITE(HTML_unit, '(10x,5a)') "Final <i>R</i> indices [I>2&sigma;]                 ", "<i>R</i>1<sup><i>c</i></sup> = ",  &
                                TRIM(CIF_parameter%refine_ls_R_factor_gt),  ", <i>wR</i>2<sup><i>d</i></sup> = ",          &
                                TRIM(CIF_parameter%refine_ls_wR_factor_gt)
   WRITE(HTML_unit, '(10x,5a)') "<i>R</i> indices (all data)                   ", "<i>R</i>1<sup><i>c</i></sup> = ",        &
                               TRIM(CIF_parameter%refine_ls_R_factor_all), ", <i>wR</i>2<sup><i>d</i></sup> = ",           &
                               TRIM(CIF_parameter%refine_ls_wR_factor_ref)
   WRITE(HTML_unit, '(10x,5a)') "Largest diff. peak and hole            ", TRIM(CIF_parameter%refine_diff_density_max),     &
                               " and ", TRIM(CIF_parameter%refine_diff_density_min), trim(HTML_elect_density_unit)
   WRITE(HTML_unit, '(a)')      ""
   WRITE(HTML_unit, '(12x,2a)') "<sup><i>a</i></sup><i>R<sub>int</sub></i>&nbsp;=&nbsp;&#8721; &#124;F<sub>o</sub><sup>2</sup> ", &
                                "&nbsp;-&nbsp;< F<sub>o</sub><sup>2</sup>>&#124; / &#8721;[F<sub>o</sub><sup>2</sup>]<br>"

   WRITE(HTML_unit, '(12x,2a)') "<sup><i>b</i></sup><i>S</i> = {&#8721; [<i>w</i>(F<sub>o</sub><sup>2</sup> ", &
                                " -  F<sub>c</sub><sup>2</sup>)<sup>2</sup>] / (<i>n</i> - <i>p</i>)}<sup>1/2</sup> <br>"

   WRITE(HTML_unit, '(12x,2a)')      "<sup><i>c</i></sup><i>R</i>1 = &#8721; &#124; &#124;<i>F<sub>o</sub></i>&#124; - &#124;", &
                                    "<i>F<sub>c</sub></i>&#124; &#124; / &#8721; &#124;<i>F<sub>o</sub></i>&#124;<br>"
   WRITE(HTML_unit, '(12x,3a)')      "<sup><i>d</i></sup><i>wR</i>2 = {&#8721; [<i>w</i>(<i>F<sub>o</sub></i><sup>2</sup> - ",  &
                                    " <i>F<sub>c</sub></i><sup>2</sup>)<sup>2</sup>] / &#8721; ",                              &
                                    "[<i>w</i>(<i>F<sub>o</sub></i><sup>2</sup>)<sup>2</sup>]}<sup>1/2</sup> <br>"
   WRITE(HTML_unit, '(12x,3a)')      "&nbsp;<i>w</i> = 1 / [&sigma;(<i>F</i><sub>o</sub><sup>2</sup>) + a<i>P</i><sup>2</sup> + ", &
                                     "b<i>P</i>] where <i>P</i> = [2<i>F</i><sub>c</sub><sup>2</sup> + ", &
                                     "MAX(<i>F</i><sub>o</sub><sup>2</sup>, 0)] /3"


   WRITE(HTML_unit, '(a)')  "</pre>"
  endif


  if(LATEX_report) then
   WRITE(LATEX_unit, '(a)')  "\newpage"
   WRITE(LATEX_unit, '(a)') "\SousTitre{Structural data}"
   WRITE(LATEX_unit, '(a)') "\begin{alltt}"

   CALL transf_moiety_string("LATEX_alltt", CIF_parameter%formula_sum , LATEX_string)
   WRITE(LATEX_unit, '(10x,4a)') "Empirical formula                      ", "\(",TRIM(adjustl(LATEX_string)),"\)"
   if(CIF_parameter%formula_sum(1:len_trim(CIF_parameter%formula_sum)) /=   &
      CIF_parameter%formula_moiety(1:len_trim(CIF_parameter%formula_moiety)))  then
    CALL transf_moiety_string("LATEX_alltt", CIF_parameter%formula_moiety , LATEX_string)
    WRITE(LATEX_unit, '(10x,4a)') "Extended  formula                      ", "\(",TRIM(adjustl(LATEX_string)),"\)"
   end if
   WRITE(LATEX_unit, '(10x,2a)') "Formula weight                         ", TRIM(CIF_parameter%formula_weight)
   WRITE(LATEX_unit, '(10x,3a)') "Temperature                            ", TRIM(CIF_parameter_DEVICE%diffracto_temperature)," K"
   WRITE(LATEX_unit, '(10x,3a)') "Wavelength                             ", &
                                  TRIM(CIF_parameter_DEVICE%diffracto_radiation_wavelength), " \AA "
   LATEX_string = CIF_parameter%symmetry_space_group
   !call check_letter_string(LATEX_string, .true.)
   call check_indice_string(LATEX_string)

   WRITE(LATEX_unit, '(10x,6a)') "Crystal system, space group            ", TRIM(CIF_parameter%symmetry_cell_setting), ", ", &
                                "\(",TRIM(LATEX_string),"\)"
   WRITE(LATEX_unit, '(10x,6a)') "Unit cell dimensions                   ", "a = ", TRIM(CIF_parameter%cell_length_a), &
                                " \AA, \(\alpha\) = ", TRIM(CIF_parameter%cell_angle_alpha), " \degre"
   WRITE(LATEX_unit, '(10x,6a)') "                                       ", "b = ", TRIM(CIF_parameter%cell_length_b), &
                                " \AA, \(\beta\) = ", TRIM(CIF_parameter%cell_angle_beta),  " \degre"
   WRITE(LATEX_unit, '(10x,6a)') "                                       ", "c = ", TRIM(CIF_parameter%cell_length_c), &
                                " \AA, \(\gamma\) = ", TRIM(CIF_parameter%cell_angle_gamma), " \degre"
   WRITE(LATEX_unit, '(10x,3a)') "Volume                                 ", TRIM(CIF_parameter%cell_volume), " \AA\(\sp{3}\)"
   WRITE(LATEX_unit, '(10x,5a)') "Z, Calculated density                  ", TRIM(CIF_parameter%formula_units_Z), ', ', &
                                TRIM(CIF_parameter%exptl_density), " (\(g.cm\sp{-1}\))"
   WRITE(LATEX_unit, '(10x,3a)') "Absorption coefficient                 ", TRIM(CIF_parameter%exptl_mu), " \(mm\sp{-1}\)"

   WRITE(LATEX_unit, '(10x,2a)') "F(000)                                 ", TRIM(CIF_parameter%F000)
   WRITE(LATEX_unit, '(10x,7a)') "Crystal size                           ", TRIM(CIF_parameter%crystal_size_max), " x ", &
                               TRIM(CIF_parameter%crystal_size_mid), " x ", TRIM(CIF_parameter%crystal_size_min), " mm"
   WRITE(LATEX_unit, '(10x,2a)') "Crystal color                          ", TRIM(CIF_parameter%crystal_colour)
   WRITE(LATEX_unit, '(10x,5a)') "Theta range for data collection        ", TRIM(CIF_parameter%theta_min), " to ", &
                                TRIM(CIF_parameter%theta_max), " \degre"
   WRITE(LATEX_unit, '(10x,4a)') "h_min, h_max                           ", TRIM(CIF_parameter%h_min), ", ", &
                                TRIM(CIF_parameter%h_max)
   WRITE(LATEX_unit, '(10x,4a)') "k_min, k_max                           ", TRIM(CIF_parameter%k_min), ", ", &
                                TRIM(CIF_parameter%k_max)
   WRITE(LATEX_unit, '(10x,4a)') "l_min, l_max                           ", TRIM(CIF_parameter%l_min), ", ", &
                                TRIM(CIF_parameter%l_max)
   WRITE(LATEX_unit, '(10x,7a)') "Reflections collected / unique         ", TRIM(CIF_parameter%diffrn_reflns_number), &
                               " / ", TRIM(CIF_parameter%reflns_number_total), " [\(\sp{a}\)R(int) = ",                &
                               TRIM(CIF_parameter%diffrn_reflns_av_R_equivalents), "]"
   WRITE(LATEX_unit, '(10x,2a)') "Reflections [\(I>2\sigma\)]                   ",    CIF_parameter%reflns_number_gt
   WRITE(LATEX_unit, '(10x,2a)') "Completeness to theta_max              ", TRIM(CIF_parameter%completeness)
   WRITE(LATEX_unit, '(10x,2a)') "Absorption correction type             ", TRIM(CIF_parameter%absorption_correction_type)
   if(CIF_parameter%absorption_correction_type(1:4) /= 'none') then
     WRITE(LATEX_unit, '(10x,4a)') "Max. and min. transmission             ", TRIM(CIF_parameter%T_max), " , ", &
                               TRIM(CIF_parameter%T_min)
   endif
   WRITE(LATEX_unit, '(10x,3a)') "Refinement method                      ", "Full-matrix least-squares on \(F\sp{2}\)"
   WRITE(LATEX_unit, '(10x,6a)') "Data / restraints / parameters         ", TRIM(CIF_parameter%reflns_number_total), " / ", &
                               TRIM(CIF_parameter%restraints_number), " / ", TRIM(CIF_parameter%refine_ls_number_parameters)
   WRITE(LATEX_unit, '(10x,2a)') "\(\sp{b}\)Goodness-of-fit                       ", TRIM(CIF_parameter%Chi2)


   WRITE(LATEX_unit, '(10x,5a)') "Final R indices [\(I>2\sigma\)]               ", "\(\sp{c}R\sb{1}\) = ",  &
                                TRIM(CIF_parameter%refine_ls_R_factor_gt),  ", \(\sp{d}wR\sb{2}\) = ",      &
                                TRIM(CIF_parameter%refine_ls_wR_factor_gt)
   WRITE(LATEX_unit, '(10x,5a)') "R indices (all data)                   ", "\(\sp{c}R\sb{1}\) = ",        &
                               TRIM(CIF_parameter%refine_ls_R_factor_all), ", \(\sp{d}wR\sb{2}\) = ",      &
                               TRIM(CIF_parameter%refine_ls_wR_factor_ref)
   WRITE(LATEX_unit, '(10x,5a)') "Largest diff. peak and hole            ", TRIM(CIF_parameter%refine_diff_density_max),     &
                               " and ", TRIM(CIF_parameter%refine_diff_density_min), "  \(\overline{e}\).\AA\(\sp{-3}\)"
   WRITE(LATEX_unit, '(a)')  "\end{alltt}"
   WRITE(LATEX_unit, '(a)')  ""
   WRITE(LATEX_unit, '(a)')  "\vspace{0.5cm}"
   WRITE(LATEX_unit, '(a)')  "\indent"
   WRITE(LATEX_unit, '(12x,a)')  "$^a R_{int}=\frac{\sum \mid {F_o^2 - \langle F_o^2 \rangle} \mid}{\sum F_o^2}$"
   WRITE(LATEX_unit, '(a)')  ""
   WRITE(LATEX_unit, '(a)')  "\vspace{0.25cm}"
   WRITE(LATEX_unit, '(a)')  "\indent"
   WRITE(LATEX_unit, '(12x,a)')  "$^b S = \left\{ \frac{\sum \lbrack w(F_o^2 - F_c^2)^2 \rbrack}{n-p} \right\}^{1/2}$"
   WRITE(LATEX_unit, '(a)')  ""


   WRITE(LATEX_unit, '(a)')  "\vspace{0.25cm}"
   WRITE(LATEX_unit, '(a)')  "\indent"
   WRITE(LATEX_unit, '(12x,a)')  "$^c R_1 = \frac{\sum{\mid F_o\mid - \mid F_c\mid}}{\sum{\mid F_o\mid}}$"
   WRITE(LATEX_unit, '(a)')  ""
   WRITE(LATEX_unit, '(a)')  "\vspace{0.25cm}"
   WRITE(LATEX_unit, '(a)')  "\indent"
   WRITE(LATEX_unit, '(12x,2a)') "$^d wR_2 = \left\{ \frac{\sum{\lbrack w(F_o^2 -  F_c^2)^2\rbrack}}",  &
                                 "{\sum{\lbrack w(F_o^2)^2 \rbrack}} \right\} ^{1/2}$"
   WRITE(LATEX_unit, '(a)')  ""
   WRITE(LATEX_unit, '(a)')  "\vspace{0.25cm}"
   WRITE(LATEX_unit, '(a)')  "\indent"
   WRITE(LATEX_unit, '(12x,2a)') "$ w = 1. / \lbrack\sigma(F_o^2) + (aP)^2 + bP \rbrack$ with ",   &
                                 "$P = \lbrack 2F_c^2 + Max(F_o^2, 0)\rbrack / 3$"

   WRITE(LATEX_unit, '(a)')  ""
  endif




  if(text_report) then
   write(text_unit, '(a)') ''
   write(text_unit, '(2x,3a)') 'Table 1: Crystal data and structure refinement for ', trim(job),'.'
   write(text_unit, '(a)') ''
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,2a)') "Empirical formula                      ", TRIM(CIF_parameter%formula_sum)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,2a)') "Formula weight                         ", TRIM(CIF_parameter%formula_weight)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,3a)') "Temperature                            ", TRIM(CIF_parameter_DEVICE%diffracto_temperature)," K"
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,3a)') "Wavelength                             ", &
                                 TRIM(CIF_parameter_DEVICE%diffracto_radiation_wavelength), " Å "
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,4a)') "Crystal system, space group            ", TRIM(CIF_parameter%symmetry_cell_setting), ", ", &
                                TRIM(CIF_parameter%symmetry_space_group)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,6a)') "Unit cell dimensions                   ", "a = ", TRIM(CIF_parameter%cell_length_a), &
                               " Å, alpha = ", TRIM(CIF_parameter%cell_angle_alpha), " °"
   WRITE(text_unit, '(10x,6a)') "                                       ", "b = ", TRIM(CIF_parameter%cell_length_b), &
                               " Å, beta = ", TRIM(CIF_parameter%cell_angle_beta),  " °"
   WRITE(text_unit, '(10x,6a)') "                                       ", "c = ", TRIM(CIF_parameter%cell_length_c), &
                               " Å, gamma = ", TRIM(CIF_parameter%cell_angle_gamma), " °"
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,3a)') "Volume                                 ", TRIM(CIF_parameter%cell_volume), " Å3"

   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,5a)') "Z, Calculated density                  ", TRIM(CIF_parameter%formula_units_Z), ', ', &
                               TRIM(CIF_parameter%exptl_density), " (g.cm-3)"
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,3a)') "Absorption coefficient                 ", TRIM(CIF_parameter%exptl_mu), " mm-1"

   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,2a)') "F(000)                                 ", TRIM(CIF_parameter%F000)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,7a)') "Crystal size                           ", TRIM(CIF_parameter%crystal_size_max), " x ", &
                               TRIM(CIF_parameter%crystal_size_mid), " x ", TRIM(CIF_parameter%crystal_size_min), " mm"
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,2a)') "Crystal color                          ", TRIM(CIF_parameter%crystal_colour)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,5a)') "Theta range for data collection        ", TRIM(CIF_parameter%theta_min), " to ", &
                                TRIM(CIF_parameter%theta_max), " °"
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,4a)') "h_min, h_max                           ", TRIM(CIF_parameter%h_min), ", ", &
                                TRIM(CIF_parameter%h_max)
   WRITE(text_unit, '(10x,4a)') "k_min, k_max                           ", TRIM(CIF_parameter%k_min), ", ", &
                                TRIM(CIF_parameter%k_max)
   WRITE(text_unit, '(10x,4a)') "l_min, l_max                           ", TRIM(CIF_parameter%l_min), ", ", &
                                TRIM(CIF_parameter%l_max)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,7a)') "Reflections collected / unique         ", TRIM(CIF_parameter%diffrn_reflns_number), &
                               " / ", TRIM(CIF_parameter%reflns_number_total), " [R(int) = ",                       &
                               TRIM(CIF_parameter%diffrn_reflns_av_R_equivalents), "]"
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,2a)') "Reflections [I>2sigma(I)]              ", CIF_parameter%reflns_number_gt
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,2a)') "Completeness to theta_max              ", TRIM(CIF_parameter%completeness)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,2a)') "Absorption correction type             ", TRIM(CIF_parameter%absorption_correction_type)
   if(CIF_parameter%absorption_correction_type(1:4) /= 'none') then
    write(text_unit, '(a)') ''
    WRITE(text_unit, '(10x,4a)') "Max. and min. transmission             ", TRIM(CIF_parameter%T_max), " , ", &
                                 TRIM(CIF_parameter%T_min)
   endif
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,3a)') "Refinement method                      ", "Full-matrix least-squares on F^2"

   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,6a)') "Data / restraints / parameters         ", TRIM(CIF_parameter%reflns_number_total), " / ", &
                               TRIM(CIF_parameter%restraints_number), " / ", TRIM(CIF_parameter%refine_ls_number_parameters)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,2a)') "Goodness-of-fit                        ", TRIM(CIF_parameter%Chi2)

   if(CIF_parameter%abs_structure_Flack(1:1) /= '?') then
    write(text_unit, '(a)') ''
    WRITE(text_unit, '(10x,2a)')"Flack parameter                        ", TRIM(CIF_parameter%abs_structure_Flack)
    end if

   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,5a)') "Final R indices [I>2sigma(I)]          ", "R1  = ",  &
                                TRIM(CIF_parameter%refine_ls_R_factor_gt),  ", wR2  = ",          &
                                TRIM(CIF_parameter%refine_ls_wR_factor_gt)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,5a)') "R indices (all data)                   ", "R1  = ",        &
                               TRIM(CIF_parameter%refine_ls_R_factor_all), ", wR2  = ",           &
                               TRIM(CIF_parameter%refine_ls_wR_factor_ref)
   write(text_unit, '(a)') ''
   WRITE(text_unit, '(10x,5a)') "Largest diff. peak and hole            ", TRIM(CIF_parameter%refine_diff_density_max),     &
                               " and ", TRIM(CIF_parameter%refine_diff_density_min), "  e.Å-3"
   WRITE(text_unit, '(a)')      ""



  endif




!--------------------------------------------------------

  if(HTML_report) then
   WRITE(HTML_unit, '(a)') "<br>"
   ! positions atomiques
   WRITE(HTML_unit, '(a)')  "<p class='title_2'>"
   WRITE(HTML_unit, '(2a)') "&nbsp;&nbsp;Atomic coordinates, site occupancy (%) and equivalent isotropic displacement ", &
                            " parameters (Å<sup>2</sup>). "
   WRITE(HTML_unit, '(a)')  "U(eq) is defined as one third of the trace of the orthogonalized U<sub>ij</sub> tensor."
   WRITE(HTML_unit, '(a)')  "</p>"
   WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
  endif

  if(text_report) then
   WRITE(text_unit, '(a)')  ""
   WRITE(text_unit, '(a)')  ""
   WRITE(text_unit, '(2x,a)')  "Table 2: Atomic coordinates, site occupancy (%) and equivalent isotropic displacement "
   WRITE(text_unit, '(2x,3a)') "         parameters (Å^2) for ", trim(job), "."
   WRITE(text_unit, '(2x,a)')  "         U(eq) is defined as one third of the trace of the orthogonalized Uij tensor."
   WRITE(text_unit, '(a)')  ""
  endif

  if(LATEX_report) then   ! positions atomiques
   WRITE(LATEX_unit, '(a)')   "\newpage"
   WRITE(LATEX_unit, '(3a)')  "\SousTitre{Atomic coordinates, site occupancy (\%) and equivalent isotropic displacement ", &
                              "parameters (\AA$^2$). U(eq) is defined as one third of the trace of the ",           &
                              "orthogonalized $U_{ij}$ tensor.}"
   WRITE(LATEX_unit, '(a)')  ""
   WRITE(LATEX_unit, '(a)') "\begin{verbatim}"
  endif


!----------------------------------------------------------------------------------
   ! verif. si le fichier.CIF a ete cree par CRYSCALC et est formatte
  created_by_cryscalc = .false.
  OPEN(UNIT = CIF_unit, FILE=TRIM(archive_cif))
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0)            exit
     IF(LEN_TRIM(CIF_input_line)==0) cycle
     if(index(CIF_input_line, 'CIF file created and formatted by CRYSCALC') /=0) then
      created_by_cryscalc = .true.
      exit
     end if
     if(index(CIF_input_line, 'data_') /=0) exit
    end do
  CLOSE(unit=CIF_unit)
!----------------------------------------------------------------------------------

   OPEN(UNIT = CIF_unit, FILE=TRIM(archive_cif))
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0)            exit
     IF(LEN_TRIM(CIF_input_line)==0) cycle
     READ(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
         case ('_atom_site_label')
          if(HTML_report) then
           WRITE(HTML_unit, '(10x,a)') "Atom     x            y            z               occ.        U(eq)"
           WRITE(HTML_unit, '(a)')     ""
          end if
          if(text_report) then
           WRITE(text_unit, '(5x,82a1)')  ("-",i=1,82)
           WRITE(text_unit, '(10x,a)') "Atom     x            y            z               occ.        U(eq)"
           WRITE(text_unit, '(5x,82a1)')  ("-",i=1,82)
           WRITE(text_unit, '(a)')     ""
          endif
          if(LATEX_report) then
           WRITE(LATEX_unit, '(5x,a)') "Atom     x            y            z               occ.        U(eq)"
           WRITE(LATEX_unit, '(a)')    ""
          endif
          do
           READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
           IF(i_error < 0) exit
           CIF_input_line = ADJUSTL(CIF_input_line)
           IF(CIF_input_line(1:1) == '_')     cycle
           IF(LEN_TRIM(CIF_input_line) == 0)  exit
           IF(CIF_input_line(1:5) == 'loop_') exit
           IF(CIF_input_line(1:1) == ';')     exit
           IF(CIF_input_line(1:1) == '#')     exit
           do
            CIF_input_line = TRIM(CIF_input_line)
            long= LEN_TRIM(CIF_input_line)
            IF(CIF_input_line(long:long) == '?') then
             CIF_input_line = CIF_input_line(1:long-1)
            ELSEIF(CIF_input_line(long:long) == '.') then
             CIF_input_line = CIF_input_line(1:long-1)
            else
             exit
            endif
           enddo
           CIF_parameter%atom = CIF_input_line
           call nombre_de_colonnes(CIF_parameter%atom, nb_col)
            if(nb_col ==5) then
            READ(CIF_parameter%atom, *) atom_label, atom_typ, atom_x, atom_y, atom_z
            atom_Ueq = '?'
            atom_occ = '?'
           else
            READ(CIF_parameter%atom, *) atom_label, atom_typ, atom_x, atom_y, atom_z, atom_Ueq, atom_adp_type, atom_occ
           endif
           if(created_by_cryscalc) then   ! le fichier .CIF n'est pas formate
            write(output_line, '(5x,2(a,4x),2a)') CIF_parameter%atom(1:4), CIF_parameter%atom(9:47), atom_occ(1:12), atom_Ueq(1:12)
           else                           ! le fichier .CIF n'est pas formate
            WRITE(output_line, '(5x,8a)') atom_label(1:4), atom_x(1:12), atom_y(1:12), atom_z(1:12), &
                                           atom_Ueq(1:12), atom_adp_type(1:12), atom_occ(1:12)
           end if
           if(HTML_report) then
            write(HTML_unit, '(5x,a)') trim(output_line)
            !WRITE(HTML_unit, '(10x,a,4x,a,4x,2a)') CIF_parameter%atom(1:4), CIF_parameter%atom(9:47), atom_occ(1:12), atom_Ueq(1:12)
           endif
           if(text_report) then
            write(text_unit, '(5x,a)') trim(output_line)
            !WRITE(text_unit, '(10x,a,4x,a,4x,2a)') CIF_parameter%atom(1:4), CIF_parameter%atom(9:47), atom_occ(1:12), atom_Ueq(1:12)
            !WRITE(text_unit, '(10x,8a)') atom_label(1:4), atom_typ(1:4), atom_x(1:12), atom_y(1:12), atom_z(1:12), &
            !                                       atom_Ueq(1:12), atom_adp_type(1:12), atom_occ(1:12)
           endif
           if(LATEX_report) then
            !WRITE(LATEX_unit, '(5x,6(a,3x))') atom_label(1:12), atom_x(1:12), atom_y(1:12), atom_z(1:12), atom_occ(1:12), &
            !                                  atom_Ueq(1:12)
            !WRITE(LATEX_unit, '(5x,a,4x,a,4x,2a)') CIF_parameter%atom(1:4), CIF_parameter%atom(9:47), atom_occ(1:12), atom_Ueq(1:12)
            write(LATEX_unit, '(a)') trim(output_line)
           endif
          END do


         case default
     end select
    end do
    close (UNIT = CIF_unit)
    if(HTML_report) then
     WRITE(HTML_unit, '(a)')  "</pre>"
     WRITE(HTML_unit, '(a)')  "<br>"
    end if
    if(LATEX_report) then
     WRITE(LATEX_unit, '(a)') '\end{verbatim}'
     WRITE(LATEX_unit, '(a)') ''
    end if
    if(TEXT_report) WRITE(text_unit, '(5x,80a1)')  ("-",i=1,80)


!-----------------------------------------------------------------------------
 ! Anisotropic displacement parameters
   if(HTML_report) then
    WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
   endif
   if(LATEX_report) then
    WRITE(LATEX_unit, '(a)') "\newpage"
    WRITE(LATEX_unit, '(a)') "\noindent"
    WRITE(LATEX_unit, '(a)') "\SousTitre{Anisotropic displacement parameters (\AA$^2$)}"
    WRITE(LATEX_unit, '(a)') ""
   endif

   OPEN(UNIT = CIF_unit, FILE=TRIM(archive_cif))
    n = 0
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0)            exit
     IF(LEN_TRIM(CIF_input_line)==0) cycle
     READ(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
         case ('_atom_site_aniso_label')


          do
           READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
           IF(i_error < 0) exit
           CIF_input_line = ADJUSTL(CIF_input_line)
           IF(CIF_input_line(1:1) == '_')     cycle
           IF(LEN_TRIM(CIF_input_line) == 0)  exit
           IF(CIF_input_line(1:5) == 'loop_') exit
           IF(CIF_input_line(1:1) == ';')     exit
           IF(CIF_input_line(1:1) == '#')     exit
           n = n + 1
           if (n == 1) then
            if(HTML_report) then
             WRITE(HTML_unit, '(a)')  "<p class='title_2'>&nbsp;&nbsp;Anisotropic displacement parameters (Å<sup>2</sup>)</p>"
             WRITE(HTML_unit, '(3a)') "<p class='retrait_1'>The anisotropic displacement factor exponent takes the form: ", &
                                      "-2&pi;<sup>2</sup> [ h<sup>2</sup> a*<sup>2</sup> U<sub>11</sub> + ... ", &
                                      "+ 2 h k a* b* U<sub>12</sub> ]. </p>"
!             WRITE(HTML_unit, '(10x,2a)') "Atom           U11            U22            U33            ",   &
!                                          "U23            U13            U12"
             WRITE(HTML_unit, '(10x,a)') "Atom     U11          U22          U33          U23          U13          U12"
             WRITE(HTML_unit, '(a)')  ""
            endif
            if(text_report) then
             WRITE(text_unit, '(a)') ""
             WRITE(text_unit, '(a)') ""
             WRITE(text_unit, '(2x,3a)') "Table 3: Anisotropic displacement parameters (Å^2) for ",trim(job),"."
             WRITE(text_unit, '(2x,a)')  "         The anisotropic displacement factor exponent takes the form:"
             WRITE(text_unit, '(2x,a)')  "         -2pi^2 [ h^2 a*^2 U11> + ...  + 2 h k a* b* U12 ]."
             WRITE(text_unit, '(a)')  ""
             WRITE(text_unit, '(5x,92a1)')  ("-",i=1,92)
             !WRITE(text_unit, '(10x,2a)') "Atom           U11            U22            U33            ",   &
             !                             "U23            U13            U12"
             WRITE(text_unit, '(10x,a)') "Atom     U11          U22          U33          U23          U13          U12"

             WRITE(text_unit, '(5x,92a1)')  ("-",i=1,92)
             WRITE(text_unit, '(a)')  ""
            endif
            if(LATEX_report) then
             WRITE(LATEX_unit, '(2a)') "The anisotropic displacement factor exponent takes the form:", &
                                       "$-2\pi [h^2 a^{*2} U_{11} + ... + 2hka^*b^*U_{12}]$"
             WRITE(LATEX_unit, '(a)')  ""
             WRITE(LATEX_unit, '(a)')  "\begin{verbatim}"
             !WRITE(LATEX_unit, '(5x,2a)') "Atom           U11            U22            U33            ",   &
             !                             "U23            U13            U12"
             WRITE(LATEX_unit, '(10x,a)') "Atom     U11          U22          U33          U23          U13          U12"

             WRITE(LATEX_unit, '(a)')  ""
            endif
           endif
           CIF_parameter%atom = CIF_input_line
           READ(CIF_parameter%atom, *) atom_label, atom_U11, Atom_U22, atom_U33, atom_U23, atom_U13, atom_U12
           if(created_by_cryscalc) then
            WRITE(output_line, '(10x,a,4x,a)')  CIF_parameter%atom(1:4), trim(CIF_parameter%atom(5:))
           else
             WRITE(output_line, '(10x,7(a,3x))') atom_label(1:4), atom_U11(1:12), atom_U22(1:12), atom_U33(1:12), &
                                                atom_U23(1:12), atom_U13(1:12), atom_U12(1:12)
           end if
           if(HTML_report) then
!            WRITE(HTML_unit, '(10x,7(a,3x))')  atom_label(1:12), atom_U11(1:12), atom_U22(1:12), atom_U33(1:12), &
!                                               atom_U23(1:12), atom_U13(1:12), atom_U12(1:12)
!            WRITE(HTML_unit, '(10x,a,4x,a)')  CIF_parameter%atom(1:4), CIF_parameter%atom(5:)
            WRITE(HTML_unit, '(a)') trim(output_line)
           endif
           if(text_report) then
!            WRITE(text_unit, '(10x,7(a,3x))')  atom_label(1:12), atom_U11(1:12), atom_U22(1:12), atom_U33(1:12), &
!                                               atom_U23(1:12), atom_U13(1:12), atom_U12(1:12)
!            WRITE(text_unit, '(10x,a,4x,a)')  CIF_parameter%atom(1:4), CIF_parameter%atom(5:)
             WRITE(text_unit, '(a)') trim(output_line)
           endif

           if(LATEX_report) then
            !WRITE(LATEX_unit, '(5x,7(a,3x))')  atom_label(1:5), atom_U11(1:12), atom_U22(1:12), atom_U33(1:12), &
            !                                   atom_U23(1:12), atom_U13(1:12), atom_U12(1:12)
            !WRITE(LATEX_unit, '(10x,a,4x,a)')  CIF_parameter%atom(1:4), CIF_parameter%atom(5:)
            WRITE(LATEX_unit, '(a)') trim(output_line)
           endif

          END do


         case default
     end select
    end do
    close (UNIT = CIF_unit)
    if(HTML_report   .and. n /=0)  WRITE(HTML_unit, '(a)')  "</pre>"
    if(LATEX_report  .and. n /=0)  WRITE(LATEX_unit, '(a)')  "\end{verbatim}"
    if(TEXT_report   .and. n /=0)  WRITE(text_unit, '(5x,110a1)')  ("-",i=1,110)


!-----------------------------------------------------------------------------




   ! distances interatomiques

   open (unit = CIF_unit, file = TRIM(archive_CIF))
   n = 0

    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0) exit
     IF(LEN_TRIM(CIF_input_line) == 0) cycle
     read(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
       case ('_geom_bond_atom_site_label_1')
        dico = ''
        n_sym_dist   = 0
        num_site_sym = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         CIF_input_line = ADJUSTL(CIF_input_line)
         IF(CIF_input_line(1:1) == '_') cycle
         IF(LEN_TRIM(CIF_input_line) == 0)  exit
         IF(CIF_input_line(1:5) == 'loop_') exit
         IF(CIF_input_line(1:1) == ';')     exit
         IF(CIF_input_line(1:1) == '#')     exit
         n = n + 1
         if(n==1) then
          if(HTML_report) then
           WRITE(HTML_unit, '(a)') "<br><p class='title_2'>&nbsp;&nbsp;Bond lengths [Å]</p>"
           WRITE(HTML_unit, '(a)')  "<pre style='font-size:14' class='color_grey1'>"
          endif
          if(text_report) then
           WRITE(text_unit, '(a)') ""
           WRITE(text_unit, '(a)') ""
           WRITE(text_unit, '(2x,3a)') "Table 4: Bond lengths [Å] for ", trim(job), '.'
           WRITE(text_unit, '(a)')  ""
          endif
          if(LATEX_report) then
           WRITE(LATEX_unit, '(a)') ""
           WRITE(LATEX_unit, '(a)') "\newpage"
           WRITE(LATEX_unit, '(a)') "\SousTitre{Bond length [\AA]}"
           WRITE(LATEX_unit, '(a)') ""
           WRITE(LATEX_unit, '(a)') "\begin{verbatim}"
          endif
         end if
         CIF_parameter%distance = CIF_input_line
         READ(CIF_parameter%distance, *) dist_atom(1), dist_atom(2), dist_value, dist_sym
         dist_sym = adjustl(dist_sym)


         report_string(1)  = ''

         if(dist_sym(1:1) /= '.') then
          n_sym_dist = n_sym_dist + 1
          site_sym(n_sym_dist) = dist_sym
          call Get_num_site_sym(n_sym_dist, site_sym, dist_sym, num_site_sym)
          call get_sym_op(dist_sym, op_numor, t)
          op_n(num_site_sym(n_sym_dist))   = op_numor
          op_t(num_site_sym(n_sym_dist),:) = t(:)

          call build_string_n(num_site_sym(n_sym_dist), report_string(1))
         endif

         write(report_string(1), '(2a)') trim(dist_atom(2)), trim(report_string(1))
         if(HTML_report) then
          write(HTML_string, '(10x,5a)') dist_atom(1)(1:4), ' - ' , report_string(1)(1:8), ' = ', TRIM(dist_value)
          !if(n == 2*int(n/2)) then
          if(multiple(n,2)) then
           write(HTML_unit, '(3a)') '<span class="ligne_1">',HTML_string(1:70),'</span>'
          else
           write(HTML_unit, '(3a)') '<span class="ligne_2">',HTML_string(1:70),'</span>'
          end if
          !write(HTML_unit, '(10x,5a)') dist_atom(1)(1:4), ' - ' , report_string(1)(1:8), ' = ', TRIM(dist_value)
         endif
         if(text_report) write(text_unit, '(10x,5a)') dist_atom(1)(1:4), ' - ' , report_string(1)(1:8), ' = ', TRIM(dist_value)
         if(LATEX_report) then
          write(LATEX_unit, '(5x,5a)') dist_atom(1)(1:4), ' - ' , report_string(1)(1:8), ' = ', TRIM(dist_value)
         endif
        end do


       case default
     end select
    END do
   close (UNIT = CIF_unit)
   if(HTML_report .and. n /=0) then
    WRITE(HTML_unit, '(a)')  "</pre>"
    WRITE(HTML_unit, '(a)')  "<br>"
   endif
   if(LATEX_report .and. n /=0)   WRITE(LATEX_unit, '(a)') "\end{verbatim} "


   if(n_sym_dist /=0) then
    if(HTML_report) then
     write(HTML_unit, '(a)')  ""
     write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
     WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey'>"
    endif
    if(text_report) then
     write(text_unit, '(a)')  ""
     write(text_unit, '(a)')  "   Symmetry transformations used to generate equivalent atoms:"
     WRITE(text_unit, '(a)') ""
    endif
    if(LATEX_report) then
     write(LATEX_unit, '(a)')  ""
     write(LATEX_unit, '(a)') "Symmetry transformations used to generate equivalent atoms:"
     write(LATEX_unit, '(a)') ""
    endif
    if(HTML_report)   call write_sym_transf(maxval(num_site_sym), op_n, op_t)
    if(LATEX_report)  call write_sym_transf(maxval(num_site_sym), op_n, op_t)
    if(TEXT_report)   call write_sym_transf(maxval(num_site_sym), op_n, op_t)
    if(HTML_report)  WRITE(HTML_unit, '(a)')  "</pre></ul>"

   endif





!-------------------------------------------------------------------------
   ! angles interatomiques


   open (unit = CIF_unit, file = TRIM(archive_CIF))
   n = 0
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0) exit
     IF(LEN_TRIM(CIF_input_line) == 0) cycle
     read(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
       case ('_geom_angle_atom_site_label_1')
        n_sym_ang    = 0
        num_site_sym = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         CIF_input_line = ADJUSTL(CIF_input_line)
         IF(CIF_input_line(1:1) == '_')     cycle
         IF(LEN_TRIM(CIF_input_line) == 0)  exit
         IF(CIF_input_line(1:5) == 'loop_') exit
         IF(CIF_input_line(1:1) == ';')     exit
         IF(CIF_input_line(1:1) == '#')     exit
         n = n + 1
         CIF_parameter%angle = CIF_input_line
         READ(CIF_parameter%angle, *) ang_atom(1), ang_atom(2), ang_atom(3), ang_value, ang_sym(1), ang_sym(2)
         ang_sym(1:2) = adjustl(ang_sym(1:2))

         i1 = index(ang_value, '(')
         i2 = index(ang_value, ')')
         if(i1 /=0 .and. i2 /=0 .and. i2 > i1) then
          read(ang_value(1:i1-1),    *) ang_value_real
          read(ang_value(i1:i2), *) ang_esd_string
         else
          read(ang_value, *) ang_value_real
          ang_esd_string = ''
         end if

         if(n==1) then
          if(HTML_report) then
           WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Angles [°]</p>"
           WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
          endif
          if(text_report) then
           WRITE(text_unit, '(a)') ""
           WRITE(text_unit, '(a)') ""
           WRITE(text_unit, '(2x,3a)') "Table 5: Angles [deg] for ", trim(job), "."
           WRITE(text_unit, '(a)')  ""
          endif
          if(LATEX_report) then
           WRITE(LATEX_unit, '(a)') ""
           WRITE(LATEX_unit, '(a)') "\newpage"
           WRITE(LATEX_unit, '(a)') "\SousTitre{Angles [\degre]}"
           WRITE(LATEX_unit, '(a)') ""
           WRITE(LATEX_unit, '(a)') "\begin{verbatim}"
          endif
         end if

         do i=1, 2
          report_string(i) = ''
          if(ang_sym(i)(1:1) /= '.') then
           n_sym_ang = n_sym_ang + 1
           site_sym(n_sym_ang) = ang_sym(i)
           call Get_num_site_sym(n_sym_ang, site_sym, ang_sym(i), num_site_sym)
           call get_sym_op(ang_sym(i), op_numor, t)
           op_n(num_site_sym(n_sym_ang))   = op_numor
           op_t(num_site_sym(n_sym_ang),:) = t(:)

           call build_string_n(num_site_sym(n_sym_ang), report_string(i))
          endif
         end do

         write(report_string(1), '(2a)') trim(ang_atom(1)), trim(report_string(1))
         write(report_string(2), '(2a)') trim(ang_atom(3)), trim(report_string(2))
         if(HTML_report) then
          !write(HTML_string, '(10x,7a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
          !                               report_string(2)(1:8), ' = ', trim(ang_value)
          write(HTML_string, '(10x,6a,F7.2,a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
                                                report_string(2)(1:8), ' = ', ang_value_real, ang_esd_string
          !if(n == 2*int(n/2)) then
          if(multiple(n,2)) then
           write(HTML_unit, '(3a)') '<span class="ligne_1">',HTML_string(1:70),'</span>'
          else
           write(HTML_unit, '(3a)') '<span class="ligne_2">',HTML_string(1:70),'</span>'
          end if
         endif
         if(text_report) then
          !write(text_unit, '(10x,7a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
          !                             report_string(2)(1:8), ' = ', trim(ang_value)
          write(text_unit, '(10x,6a,F7.2,a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
                                              report_string(2)(1:8), ' = ', ang_value_real, ang_esd_string
         endif
         if(LATEX_report) then
          !write(LATEX_unit, '(5x,7a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
          !                             report_string(2)(1:8), ' = ', trim(ang_value)
          write(LATEX_unit, '(10x,6a,F7.2,a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
                                               report_string(2)(1:8), ' = ', ang_value_real, ang_esd_string
         endif

        end do
       case default
     end select
    END do
   close (UNIT = CIF_unit)
   if(HTML_report .and. n /=0) then
    WRITE(HTML_unit, '(a)')  "</pre>"
    WRITE(HTML_unit, '(a)')  "<br>"
   endif
   if(LATEX_report .and. n /=0) WRITE(LATEX_unit, '(a)') "\end{verbatim}"


  if(n_sym_ang /=0) then
   if(HTML_report) then
    WRITE(HTML_unit, '(a)')  "</pre>"
    write(HTML_unit, '(a)')  ""
    write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
    WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
   end if
   if(text_report) then
    write(text_unit, '(a)')  ""
    write(text_unit, '(a)')  "   Symmetry transformations used to generate equivalent atoms:"
    WRITE(text_unit, '(a)') ""
   endif
   if(LATEX_report) then
    write(LATEX_unit, '(a)')  ""
    write(LATEX_unit, '(a)')  " Symmetry transformations used to generate equivalent atoms:"
    WRITE(LATEX_unit, '(a)')  ""
   endif
   if(HTML_report)  call write_sym_transf(maxval(num_site_sym), op_n, op_t)
   if(LATEX_report) call write_sym_transf(maxval(num_site_sym), op_n, op_t)
   if(TEXT_report)  call write_sym_transf(maxval(num_site_sym), op_n, op_t)
  endif
  if(HTML_report)    WRITE(HTML_unit, '(a)')  "</pre></ul>"


!-------------------------------------------------------------------------
   ! angles de torsion interatomiques
   open (unit = CIF_unit, file = TRIM(archive_CIF))
   n = 0
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0) exit
     IF(LEN_TRIM(CIF_input_line) == 0) cycle
     read(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
       case ('_geom_torsion_atom_site_label_1')
        n_sym_tor = 0
        num_site_sym = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         CIF_input_line = ADJUSTL(CIF_input_line)
         IF(CIF_input_line(1:1) == '_') cycle
         IF(LEN_TRIM(CIF_input_line) == 0)  exit
         IF(CIF_input_line(1:5) == 'loop_') exit
         IF(CIF_input_line(1:1) == ';')     exit
         IF(CIF_input_line(1:1) == '#')     exit

         n = n + 1
         if(n==1) then
          if(HTML_report) then
           WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Torsion angles [°]</p>"
           WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
          endif
          if(text_report) then
           WRITE(text_unit, '(a)') ""
           WRITE(text_unit, '(a)') ""
           WRITE(text_unit, '(2x,3a)') "Table 6: Torsion angles [deg] for ", trim(job), '.'
           WRITE(text_unit, '(a)')  ""
          endif
          if(LATEX_report) then
           WRITE(LATEX_unit, '(a)') ""
           WRITE(LATEX_unit, '(a)') "\newpage"
           WRITE(LATEX_unit, '(a)') "\SousTitre{Torsion angles [\degre]}"
           WRITE(LATEX_unit, '(a)') ""
           WRITE(LATEX_unit, '(a)') "\begin{verbatim}"
          endif
         endif
         CIF_parameter%torsion_angle = CIF_input_line
         READ(CIF_parameter%torsion_angle, *) (torsion_ang_atom(i), i=1,4), torsion_ang_value, (torsion_sym(i), i=1,4)
         torsion_sym(1:4) = adjustl(torsion_sym(1:4))
         i1 = index(torsion_ang_value, '(')
         i2 = index(torsion_ang_value, ')')
         if(i1 /=0 .and. i2 /=0 .and. i2 > i1) then
          read(torsion_ang_value(1:i1-1), *) torsion_value_real
          !read(torsion_ang_value(i1+1:i2-1), *) torsion_value_real
          read(torsion_ang_value(i1:i2),*)      torsion_esd_string
          read(torsion_ang_value(i1+1:i2-1), *) torsion_esd
         else
          read(torsion_ang_value, *) torsion_value_real
         endif
         !if(abs(torsion_value_real) > CIF_torsion_limit) then
         ! n = n -1
         ! cycle
         !endif

         do i=1, 4
          report_string(i) = ''
          if(torsion_sym(i) /= '.') then
           n_sym_tor = n_sym_tor + 1
           site_sym(n_sym_tor) = torsion_sym(i)
           call Get_num_site_sym(n_sym_tor, site_sym, torsion_sym(i), num_site_sym)
           call get_sym_op(torsion_sym(i), op_numor, t)
           op_n(num_site_sym(n_sym_tor))   = op_numor
           op_t(num_site_sym(n_sym_tor),:) = t(:)
           call build_string_n(num_site_sym(n_sym_tor), report_string(i))
          endif
         end do


       write(report_string(1), '(2a)') trim(torsion_ang_atom(1)), trim(report_string(1))
       write(report_string(2), '(2a)') trim(torsion_ang_atom(2)), trim(report_string(2))
       write(report_string(3), '(2a)') trim(torsion_ang_atom(3)), trim(report_string(3))
       write(report_string(4), '(2a)') trim(torsion_ang_atom(4)), trim(report_string(4))
       if(HTML_report) then
        !WRITE(HTML_string, '(10x, 8a,1a12)') report_string(1)(1:8), ' - ', &
        !                                report_string(2)(1:8), ' - ', &
        !                                report_string(3)(1:8), ' - ', &
        !                                report_string(4)(1:8), ' = ', trim(torsion_ang_value)

        WRITE(HTML_string, '(10x, 8a,F7.2,a)') report_string(1)(1:8), ' - ', &
                                        report_string(2)(1:8), ' - ', &
                                        report_string(3)(1:8), ' - ', &
                                        report_string(4)(1:8), ' = ', torsion_value_real, trim(torsion_esd_string)

        !if(n == 2*int(n/2)) then
        if(multiple(n,2)) then
         write(HTML_unit, '(3a)') '<span class="ligne_1">',HTML_string(1:70),'</span>'
        else
         write(HTML_unit, '(3a)') '<span class="ligne_2">',HTML_string(1:70),'</span>'
        end if
       endif
       if(text_report) then
 !       WRITE(text_unit, '(10x, 9a)') report_string(1)(1:8), ' - ', &
 !                                     report_string(2)(1:8), ' - ', &
 !                                     report_string(3)(1:8), ' - ', &
 !                                     report_string(4)(1:8), ' = ', trim(torsion_ang_value)
         WRITE(text_unit, '(10x, 8a,F7.2,a)') report_string(1)(1:8), ' - ', &
                                             report_string(2)(1:8), ' - ', &
                                             report_string(3)(1:8), ' - ', &
                                             report_string(4)(1:8), ' = ', torsion_value_real, trim(torsion_esd_string)

       endif
       if(LATEX_report) then
        !WRITE(LATEX_unit, '(5x, 9a)') report_string(1)(1:8), ' - ', &
        !                              report_string(2)(1:8), ' - ', &
        !                              report_string(3)(1:8), ' - ', &
        !                              report_string(4)(1:8), ' = ', trim(torsion_ang_value)
        WRITE(LATEX_unit, '(5x, 8a,F7.2,a)') report_string(1)(1:8), ' - ', &
                                      report_string(2)(1:8), ' - ', &
                                      report_string(3)(1:8), ' - ', &
                                      report_string(4)(1:8), ' = ', torsion_value_real, trim(torsion_esd_string)
       endif
      end do


       case default
     end select
    END do
   close (UNIT = CIF_unit)
   if(HTML_report  .and. n /=0)   WRITE(HTML_unit, '(a)')   "</pre><br>"
   if(LATEX_report .and. n /=0)   WRITE(LATEX_unit, '(a)')  "\end{verbatim}"


  if(n_sym_tor /=0) then
    if(HTML_report) then
     write(HTML_unit, '(a)')  ""
     write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
     WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
    endif
    if(text_report) then
     write(text_unit, '(a)')  ""
     write(text_unit, '(a)')  "   Symmetry transformations used to generate equivalent atoms:"
     WRITE(text_unit, '(a)') ""
    endif
    if(LATEX_report) then
     write(LATEX_unit, '(a)')  ""
     write(LATEX_unit, '(a)')  " Symmetry transformations used to generate equivalent atoms:"
     WRITE(LATEX_unit, '(a)')  ""
    endif
    if(HTML_report)  call write_sym_transf(maxval(num_site_sym), op_n, op_t)
    if(LATEX_report) call write_sym_transf(maxval(num_site_sym), op_n, op_t)
    if(TEXT_report)  call write_sym_transf(maxval(num_site_sym), op_n, op_t)
   endif

   if(HTML_report) WRITE(HTML_unit, '(a)')  "</pre></ul>"


!-------------------------------------------------------------------------
   ! liaisons H


   n_sym_htab = 0
   open (unit = CIF_unit, file = TRIM(archive_CIF))
   n = 0
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0) exit
     IF(LEN_TRIM(CIF_input_line) == 0) cycle
     read(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
       case ('_geom_hbond_atom_site_label_D')
        n_sym_htab = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         CIF_input_line = ADJUSTL(CIF_input_line)
         IF(CIF_input_line(1:1) == '_') cycle
         IF(LEN_TRIM(CIF_input_line) == 0)  exit
         IF(CIF_input_line(1:5) == 'loop_') exit
         IF(CIF_input_line(1:1) == ';')     exit
         IF(CIF_input_line(1:1) == '#')     exit
         n = n + 1
         CIF_parameter%Hbond = CIF_input_line
         READ(CIF_parameter%Hbond, *) site_label_D, site_label_H, site_label_A, dist_DH, dist_HA, dist_DA, angle_DHA, &
                                    site_sym_A
         if( n == 1) then
          if(HTML_report) then
           WRITE(HTML_unit, '(a)') "<br>"
           WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Hydrogen bonds [Å and deg.] </p>"
           WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
          endif
          if(text_report) then
           WRITE(text_unit, '(a)') ""
           WRITE(text_unit, '(a)') ""
           WRITE(text_unit, '(2x,3a)') "Table 7: Hydrogen bonds [A and deg.] for ", trim(job), '.'
           WRITE(text_unit, '(a)')  ""
          endif
          if(LATEX_report) then
           WRITE(LATEX_unit, '(a)') "\newpage"
           WRITE(LATEX_unit, '(a)') "\SousTitre{Hydrogen bonds [\AA~ and \degre]}"
           WRITE(LATEX_unit, '(a)') ""
           WRITE(LATEX_unit, '(a)') "\begin{verbatim}"
          endif
         end if

         report_string(1) = ''
         if(site_sym_A(1:1) /= '.') then
          n_sym_htab = n_sym_htab + 1
          site_sym(n_sym_htab) = site_sym_A
          call Get_num_site_sym(n_sym_htab, site_sym, site_sym_A, num_site_sym)
          call Get_sym_op(site_sym_A, op_numor, t)
          op_n(num_site_sym(n_sym_htab))    = op_numor
          op_t(num_site_sym(n_sym_htab), :) = t(:)
          call build_string_n(num_site_sym(n_sym_htab), report_string(1))
         endif

         if(HTML_report) then
           WRITE(HTML_string, '(10x,11a)')  site_label_D(1:4), ' - ', site_label_H(1:4), '.... ', site_label_A(1:4), &
                                            report_string(1)(1:8),  dist_DH(1:12), dist_HA(1:12), dist_DA(1:12), &
                                            angle_DHA(1:12)
          !if(n == 2*int(n/2)) then
           if(multiple(n,2)) then
            write(HTML_unit, '(3a)') '<span class="ligne_1">',HTML_string(1:92),'</span>'
           else
            write(HTML_unit, '(3a)') '<span class="ligne_2">',HTML_string(1:92),'</span>'
           end if
          endif
          if(text_report) then
           WRITE(text_unit, '(10x,11a)')  site_label_D(1:4), ' - ', site_label_H(1:4), '.... ', site_label_A(1:4), &
                                          report_string(1)(1:8),  dist_DH(1:12), dist_HA(1:12), dist_DA(1:12), &
                                          angle_DHA(1:12)
          endif
          if(LATEX_report) then
           WRITE(LATEX_unit, '(5x,11a)')  site_label_D(1:4), ' - ', site_label_H(1:4), '.... ', site_label_A(1:4), &
                                          report_string(1)(1:8),  dist_DH(1:12), dist_HA(1:12), dist_DA(1:12), &
                                          angle_DHA(1:12)
          endif

        end do


       case default
     end select
    END do
   close (UNIT = CIF_unit)

   if(HTML_report  .and. n/=0)  write(HTML_unit, '(a)')  "</pre>"
   if(LATEX_report .and. n/=0)  write(LATEX_unit, '(a)')  "\end{verbatim}"

   if(n_sym_htab /=0) then
    if(HTML_report) then
     WRITE(HTML_unit, '(a)')  "</pre>"
     write(HTML_unit, '(a)')  ""
     write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
     WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
    endif
    if(text_report) then
     write(text_unit, '(a)')  ""
     write(text_unit, '(a)')  "   Symmetry transformations used to generate equivalent atoms:"
     WRITE(text_unit, '(a)') ""
    endif
    if(LATEX_report) then
     WRITE(LATEX_unit, '(a)')  ""
     write(LATEX_unit, '(a)')  " Symmetry transformations used to generate equivalent atoms:"
     WRITE(LATEX_unit, '(a)')  ""
    endif
    if(HTML_report)  call write_sym_transf(maxval(num_site_sym), op_n, op_t)
    if(LATEX_report) call write_sym_transf(maxval(num_site_sym), op_n, op_t)
    if(TEXT_report)  call write_sym_transf(maxval(num_site_sym), op_n, op_t)
   endif
   if(HTML_report) WRITE(HTML_unit, '(a)')  "</pre></ul>"



!-------------------------------------------------------------------------


  ! recherche de la presence de fichier.GIF
  if(.not. TEXT_report) then

  file_exist = .false.
  do
   i1 = index(CIF_parameter_DEVICE%diffracto_temperature, "(")
   if(i1 /=0) then
    write(GIF_file, '(4a)') trim(job), "_", CIF_parameter_DEVICE%diffracto_temperature(1:i1-1), "K_ortep.gif"
   else
    write(GIF_file, '(4a)') trim(job), "_", TRIM(CIF_parameter_DEVICE%diffracto_temperature), "K_ortep.gif"
   endif
   call test_file_exist(trim(GIF_file), file_exist, 'out')
   if(file_exist) exit

   GIF_file = "platon_ortep.gif"
   call test_file_exist(trim(GIF_file), file_exist, 'out')
   if(file_exist) exit

   GIF_file= "ortep.gif"
   call test_file_exist(trim(GIF_file), file_exist, 'out')
   if(file_exist) exit

   GIF_file = "platon_jobte.gif"
   call test_file_exist(trim(GIF_file), file_exist, 'out')
   if(file_exist) exit

   GIF_file = "platon_Ite.gif"
   call test_file_exist(trim(GIF_file), file_exist, 'out')
   exit
  end do

  if(file_exist) then
   if(HTML_report) then
    WRITE(HTML_unit, '(a)') ""
    WRITE(HTML_unit, '(a)') "<br>"
    WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Structure visualisation</p>"
    WRITE(HTML_unit, '(3a)') '<center><img width=600 src="', trim(GIF_file),'"></center>'
   endif
   if(LATEX_report) then   ! .gif file can not be included in PDFLATEX
   ! convertion du fichier .gif en .png avec Convert (ImageMagick)
    i1 = index(GIF_file, '.', back=.true.)
    PNG_file = GIF_file(1:i1)//'png'

    write(DOS_command, '(4a)') 'convert -bordercolor rgb(255,255,255) -border 1 ', trim(GIF_file), ' ', trim(PNG_file)
    call system(trim(DOS_command))

    WRITE(LATEX_unit, '(a)') ""
    WRITE(LATEX_unit, '(a)') "\newpage"
    WRITE(LATEX_unit, '(a)') "\SousTitre{Structure visualisation}"
    WRITE(LATEX_unit, '(a)') ""
    write(LATEX_unit, '(a)')  '\begin{figure}[h]'
    write(LATEX_unit, '(a)')  ' \centering'
    write(LATEX_unit, '(3a)') ' \includegraphics[width=350.pt]{',trim(PNG_file),'}'
    write(LATEX_unit, '(a)')  '\end{figure}'
   endif
  endif
  end if ! fin de la condition if (.not. TXT_report)



  if(HTML_report) then
   WRITE(HTML_unit, '(a)')  "<HR>"
   WRITE(HTML_unit, '(a)')  "</div class='cadre'>"
   WRITE(HTML_unit, '(a)')  "<p class='retrait_1'>"
   WRITE(HTML_unit, '(6a)') "<i>This HTML report has been created through ",                        &
                            "<A HREF='http://", trim(CRYSCALC%url), "'>CRYSCALC</A> (", &
                            trim(cryscalc%version), ")."
   WRITE(HTML_unit, '(5a)') "Please report bugs and problems to <A HREF='mailto:", &
                             trim(CRYSCALC%mail), "'>", trim(CRYSCALC%mail), "</A></p>"

   WRITE(HTML_unit, '(a)')  "</BODY>"
   WRITE(HTML_unit, '(a)')  "</HTML>"
   CLOSE(UNIT=HTML_unit)
  endif

  if(LATEX_report) then
   call LATEX_End
   close(unit= LATEX_unit)
  endif

  if(text_report) then
   close(unit=text_unit)
   call write_info('')
   call write_info('    >>> '//TRIM(TEXT_structural_report_file)//' file has been created.')
   call write_info('')
  endif

  if(HTML_report) then
   call write_info('')
   call write_info('    >>> '//TRIM(HTML_structural_report_file)//' file has been created.')
   call write_info('')
  endif

  if(LATEX_report) then
   call write_info('')
   call write_info('    >>> '//TRIM(LATEX_structural_report_file)//' file has been created.')
   call write_info('')
  endif


! call launch_browser('\structural_report.HTML', 'internal')
 if(HTML_report) then
  IF(my_browser%exist) then
   call launch_browser(TRIM(HTML_structural_report_file))
   call write_info('')
   call write_info('  Please wait. Your browser will be launched to display the HTML report.')
   call write_info('')
  else
   call write_info('')
   write(message_text, '(3a)') '  No browser defined in the ', trim(cryscalc%ini), ' setting file !!'
   call write_info(trim(message_text))
   call write_info('')
  endif
 endif

 if(LATEX_report) then
  IF(my_pdfLATEX%exist) then
   call launch_pdfLATEX(TRIM(LATEX_structural_report_file))
   call write_info('')
   call write_info('  Please wait. Your LATEX to PDF convertor will be launched.')
   call write_info('')
   call write_info(' LATEX pdf : '//trim(PDF_structural_report_file))
   call system(TRIM(PDF_structural_report_file))

   call system('del '//trim(report_logo(1)))
   call system('del '//trim(report_logo(2)))

  else
   call write_info('')
   write(message_text, '(3a)') '  No LATEX to PDF convertor defined in the ', trim(cryscalc%ini), ' setting file !!'
   call write_info(trim(message_text))
   call write_info('')
  endif
 endif


  return

end subroutine create_report



!---------------------------------------------------------------------
subroutine get_sym_op(input_string, numor_op, t)
 use cryscalc_module, only : message_text
 use io_module
 implicit none
  character (len=*), intent(in)       :: input_string
  integer              , intent(out)  :: numor_op
  integer, dimension(3), intent(out)  :: t
  integer                             :: i1, i2


  i1 = index(input_string, '_')
  i2 = index(input_string, '.')

  if(i1==0) then
   if(i2/=0) then
    numor_op = 1
   else
    read(input_string, *) numor_op   ! decommenté oct 2011
   endif
   t(1) = 5
   t(2) = 5
   t(3) = 5
  else
   read(input_string(1:i1-1), *) numor_op
   read(input_string(i1+1: i1+1), *) t(1)
   read(input_string(i1+2: i1+2), *) t(2)
   read(input_string(i1+3: i1+3), *) t(3)
  endif

  t(:) = t(:) - 5

 return
end subroutine get_sym_op


!------------------------------------------------------------------------
subroutine TRANSF_moiety_string(code, input_string, output_string)
 use macros_module,   only : nombre_de_colonnes, check_character
 use cryscalc_module, only : debug_proc
 implicit none
  character (len=*),   intent(in)     :: code             !(HTML/CIF/HTML)
  character (len=*),   intent(inout)  :: input_string
  character (len=512), intent(out)    :: output_string
  character (len=256)                 :: tmp_string
  character (len=8), dimension(2)     :: balise
  integer                             :: i, j, i1, nb_col, num, num_alpha
  character (len=8),  dimension(32)   :: part_string
  LOGICAL                             :: alpha_char, numeric_char, espace_vide
  character (len=1)                   :: espace
  character (len=8)                   :: espace_2

 if(debug_proc%level_3)  call write_debug_proc_level(3, "transf_moiety_string ("//trim(input_string)//")")


  output_string = ''
  espace_vide   = .false.
  if(code == 'HTML') then
   balise(1) = '<sub>'
   balise(2) = '</sub>'
   espace    = ' '
   espace_2  = '&nbsp;'   ! espace insecable
  elseif(code == 'CIF') then
   balise(1) = '~'
   balise(2) = '~'
   espace    = ' '
   espace_2  = ' '
   espace_vide = .true.
  elseif(code == 'LATEX') then
   balise(1) = '_{'
   balise(2) = '}'
   espace    = ' '
   espace_2    = '~'   ! espace insecable
  elseif(code == 'LATEX_alltt') then
   balise(1) = '\sb{'
   balise(2) = '}'
   espace    = ' '
   espace_2  = ' '
   espace_vide = .true.
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nb_col = 0
  tmp_string = input_string
  do
   i1 = index(tmp_string, ' ')
   if(i1 ==0) exit
   if(len_trim(tmp_string) ==0) exit
   nb_col = nb_col + 1
   read(tmp_string(1:i1-1), '(a)') part_string(nb_col)
   tmp_string = tmp_string(i1+1:)
  end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! pb si la chaine contient le caractere virgule ','
  !input_string = adjustl(input_string)
  !call nombre_de_colonnes(input_string, nb_col)
  !if(nb_col == 0) then
  ! output_string = input_string
  ! return
  !endif
  !read(input_string, *) part_string(1:nb_col)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1, nb_col

   num       = 0
   num_alpha = 0

   do j=1, len_trim(part_string(i))
    !num       = 0
    !num_alpha = 0

    call check_character (part_string(i)(j:j), alpha_char, numeric_char)
    if(alpha_char) then
     num_alpha = num_alpha + 1
     if(num ==0) then
      if(num_alpha == 1 .and. i/=1) then
       !write(OUTPUT_string, '(3a)') trim(OUTPUT_string), espace(1:1) , part_string(i)(j:j)  ! espaces
       if(.not. espace_vide) then
        write(OUTPUT_string, '(3a)') trim(OUTPUT_string), trim(espace_2) , part_string(i)(j:j)  ! espaces
       else
        write(OUTPUT_string, '(3a)') trim(OUTPUT_string), espace_2(1:1) , part_string(i)(j:j)  ! espaces
       end if
       !write(OUTPUT_string, '(2a)') trim(OUTPUT_string),  part_string(i)(j:j)        ! pas d'espaces dans la chaine
      else
       write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)
      endif
     else
      if(num_alpha == 1) then
       write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)
      else
       write(OUTPUT_string, '(3a)') trim(OUTPUT_string), trim(balise(2)), part_string(i)(j:j)
      endif
      !write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)

     end if
     num = 0

    else
     num       = num+1
     if(num_alpha /=0) then
      if(num == 1 .and. j/=1) then
       write(OUTPUT_string, '(3a)') trim(OUTPUT_string), trim(balise(1)), part_string(i)(j:j)
      else
       write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)
      endif
     else
      !write(OUTPUT_string, '(3a)') trim(OUTPUT_string), espace(1:1), part_string(i)(j:j)
      if(num == 1 .and. i/=1) then
       if(.not. espace_vide) then
        write(OUTPUT_string, '(3a)') trim(OUTPUT_string), trim(espace_2), part_string(i)(j:j)
       else
        write(OUTPUT_string, '(3a)') trim(OUTPUT_string), espace_2(1:1), part_string(i)(j:j)
       end if
      else
       write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)
      endif
      num_alpha = 0
     endif
     !num_alpha = 0
    end if
   end do


   if(num/=0) then
    write(OUTPUT_string, '(2a)') trim(OUTPUT_string), trim(balise(2))
   !else
    !write(OUTPUT_string, '(2a)') trim(OUTPUT_string), espace(1:1)
   endif

  end do



 return
end subroutine TRANSF_moiety_string


!---------------------------------------------------------------------

subroutine Get_num_site_sym(n_sym, site_sym, current_sym, num_site_sym)
 implicit none
  integer,                             intent(in)    :: n_sym
  CHARACTER (LEN=12), dimension(500),  intent(in)    :: site_sym
  CHARACTER (LEN=12),                  intent(in)    :: current_sym
  INTEGER           , dimension(500),  intent(inout) :: num_site_sym
  INTEGER                                            :: i


  do i = 1, n_sym-1
   if(current_sym == site_sym(i)) then
    num_site_sym(n_sym) = num_site_sym(i)
    return
   end if
  end do

  num_site_sym(n_sym) = MAXVAL(num_site_sym) + 1
  return

end subroutine Get_num_site_sym

!---------------------------------------------------------------------
subroutine  Get_num_site_sym_new(n_sym, current_sym, num_site_sym, dico)
 implicit none
  integer,                             intent(in)    :: n_sym
  CHARACTER (LEN=12),                  intent(in)    :: current_sym
  INTEGER           , dimension(500),  intent(inout) :: num_site_sym
  character (len=12), dimension(500),  intent(inout) :: dico

  INTEGER                                            :: i

  num_site_sym(1) = 0
  do i = 1, n_sym

   if(dico(i) == '') then
    dico(i) = current_sym
    num_site_sym(n_sym) = num_site_sym(i)+1
    return
   else
    if(dico(i) == current_sym) then
     num_site_sym(n_sym) = i
     return
    end if
   end if

  end do

 return
end subroutine Get_num_site_sym_new

!---------------------------------------------------------------------
subroutine build_string_n(n, report_string)

 implicit none
  integer,            intent(in)          :: n
  character(len=256), intent(inout)       :: report_string

  if(n < 10) then
   write(report_string, '(a,i1)') '_#', n
  elseif(n < 100) then
   write(report_string, '(a,i2)') '_#', n
  else
   write(report_string, '(a,i3)') '_#', n
  endif

 return
end subroutine build_string_n

!---------------------------------------------------------------------
subroutine write_sym_transf(n_sym,  op_n, op_t)
 use structural_report_module, only : op_string
 use cryscalc_module,          only : HTML_unit, text_unit, LATEX_unit, HTML_report, text_report, LATEX_report, message_text
 use IO_module
 implicit none

  integer,                             intent(in) :: n_sym

  integer, dimension(500),             intent(in) :: op_n       ! numero de l'op. de sym. associe
  integer, dimension(500,3),           intent(in) :: op_t       ! partie translation
  integer                                         :: i




    do i=1, n_sym
     !if(op_n(i) < 10) then
     if(i < 10) then
      if(HTML_report) then
       write(HTML_unit, '(10x,a,i1,3x,2a,3(i2,a))') " #",i,  op_string(op_n(i))(1:24), &
                                                 "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
      endif
      if(text_report) then
       write(text_unit, '(10x,a,i1,3x,2a,3(i2,a))') " #",i, op_string(op_n(i))(1:24), &
                                                    "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
      endif
      if(LATEX_report) then
       if(i==1) write(LATEX_unit, '(a)') '\begin{verbatim}'
        write(LATEX_unit, '(5x,a,i1,3x,2a,3(i2,a))') " #",i,  op_string(op_n(i))(1:24), &
                                                     "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
        if(i==n_sym) write(LATEX_unit, '(a)') '\end{verbatim}'
       endif
     !elseif(op_n(i) < 100) then
     elseif(i < 100) then
      if(HTML_report) then
       write(HTML_unit, '(10x,a,i2,2x,2a,3(i2,a))') " #",i, op_string(op_n(i))(1:24), &
                                                    "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
      endif
      if(text_report) then
       write(text_unit, '(10x,a,i2,2x,2a,3(i2,a))') " #",i, op_string(op_n(i))(1:24), &
                                                    "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
      endif
      if(LATEX_report) then
       if(i==1) write(LATEX_unit, '(a)') '\begin{verbatim}'
        write(LATEX_unit, '(5x,a,i2,2x,2a,3(i2,a))') " #",i, op_string(op_n(i))(1:24), &
                                                     "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
       if(i==n_sym) write(LATEX_unit, '(a)') '\end{verbatim}'
      endif

     else
      if(HTML_report) then
       write(HTML_unit, '(10x,a,i3,2x,2a,3(i2,a))') " #",i, op_string(op_n(i))(1:24), &
                                                    "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
      endif
      if(text_report) then
       write(text_unit, '(10x,a,i3,2x,2a,3(i2,a))') " #",i, op_string(op_n(i))(1:24), &
                                                    "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
      endif
      if(LATEX_report) then
       if(i==1) write(LATEX_unit, '(a)') '\begin{verbatim}'
       write(LATEX_unit, '(5x,a,i3,2x,2a,3(i2,a))') " #",i, op_string(op_n(i))(1:24), &
                                                   "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
       if(i==n_sym) write(LATEX_unit, '(a)') '\end{verbatim}'
      endif
     endif
    end do

 return
end subroutine write_sym_transf

!------------------------------------------------------------------------------------------------------------
subroutine Get_CIF_value(CIF_unit, CIF_value_IN, CIF_value_OUT)
 USE IO_module, only  : write_info

 implicit none
  integer,            intent(in)   :: CIF_unit
  character(len=256), intent(IN)    :: CIF_value_IN
  character(len=256), intent(OUT)   :: CIF_value_OUT
  character(len=256)               :: read_line
  integer                          :: i_error, long


        READ(CIF_value_IN, '(a)') CIF_value_OUT
        IF(LEN_TRIM(CIF_value_OUT) == 0) then
         READ(CIF_unit, '(a)', IOSTAT=i_error) read_line
         IF(i_error < 0) return
         IF(i_error /=0) then
          call write_info('')
          call write_info(' !! Error reading archive.CIF file !!')
         call write_info('')
         return
         endif
         READ(read_line, '(a)') CIF_value_OUT
         CIF_value_OUT = ADJUSTL(CIF_value_OUT)
        endif


        long = LEN_TRIM(CIF_value_OUT)
        IF(CIF_value_OUT(1:1)== "'" .AND. CIF_value_OUT(long:long) == "'") CIF_value_OUT = CIF_value_OUT(2:long-1)


  return
end subroutine Get_CIF_value



!!-------------------------------------------------------------------------------------!!
subroutine check_indice_string(string)
 ! met un caractere en indice si necessaire (chaine = groupe d'espace)
 use cryscalc_module, only : HTML_report, LATEX_report
 use macros_module,  only : replace_car2
 implicit none
 character (len=*), intent(inout)  :: string


 if(HTML_report) then
  if(index(string, "21") /= 0)   string = replace_car2(string, "21", "2<sub>1</sub>")
  if(index(string, "31") /= 0)   string = replace_car2(string, "31", "3<sub>1</sub>")
  if(index(string, "32") /= 0)   string = replace_car2(string, "32", "3<sub>2</sub>")
  if(index(string, "41") /= 0)   string = replace_car2(string, "41", "4<sub>1</sub>")
  if(index(string, "42") /= 0)   string = replace_car2(string, "42", "4<sub>2</sub>")
  if(index(string, "43") /= 0)   string = replace_car2(string, "43", "4<sub>3</sub>")
  if(index(string, "61") /= 0)   string = replace_car2(string, "61", "6<sub>1</sub>")
  if(index(string, "62") /= 0)   string = replace_car2(string, "62", "6<sub>2</sub>")
  if(index(string, "63") /= 0)   string = replace_car2(string, "63", "6<sub>3</sub>")
  if(index(string, "64") /= 0)   string = replace_car2(string, "64", "6<sub>4</sub>")
  if(index(string, "65") /= 0)   string = replace_car2(string, "65", "6<sub>5</sub>")
 endif

 if(LATEX_report) then
  if(index(string, "21") /= 0)   string = replace_car2(string, "21", "2\sb{1}")
  if(index(string, "31") /= 0)   string = replace_car2(string, "31", "3\sb{1}")
  if(index(string, "32") /= 0)   string = replace_car2(string, "32", "3\sb{2}")
  if(index(string, "41") /= 0)   string = replace_car2(string, "41", "4\sb{1}")
  if(index(string, "42") /= 0)   string = replace_car2(string, "42", "4\sb{2}")
  if(index(string, "43") /= 0)   string = replace_car2(string, "43", "4\sb{3}")
  if(index(string, "61") /= 0)   string = replace_car2(string, "61", "6\sb{1}")
  if(index(string, "62") /= 0)   string = replace_car2(string, "62", "6\sb{2}")
  if(index(string, "63") /= 0)   string = replace_car2(string, "63", "6\sb{3}")
  if(index(string, "64") /= 0)   string = replace_car2(string, "64", "6\sb{4}")
  if(index(string, "65") /= 0)   string = replace_car2(string, "65", "6\sb{5}")
  if(index(string, "-1") /= 0)   string = replace_car2(string, "-1", "\overline{1}")
  if(index(string, "-3") /= 0)   string = replace_car2(string, "-3", "\overline{3}")
  if(index(string, "-4") /= 0)   string = replace_car2(string, "-4", "\overline{4}")

 endif

 return
 end subroutine Check_indice_string


!!--------------------------------------------------------------------------!!
subroutine Check_letter_string(string, latexTT)
 ! met un caractere en indice si necessaire (chaine = groupe d'espace)
 use cryscalc_module, only : HTML_report, LATEX_report
 use macros_module,   only : remove_car, check_character
 implicit none
 character (len=*), intent(inout)  :: string
 logical, intent(in)               :: latexTT
 character (len=256)               :: SG_string
 integer                           :: i
 logical                           :: alpha_car, numeric_car


 !string = remove_car(string, ' ')
 !SG_string = string
 SG_string = ''

 do i=1 , len_trim(string)
  if(string(i:i) == '/') then
   SG_string = trim(SG_string)//string(i:i)
   cycle
  elseif(string(i:i) == ' ') then
   if(HTML_report) then
    SG_string = trim(SG_string)//"&nbsp;"
    cycle
   elseif(LATEX_report) then
    if(latexTT) then
     SG_string = trim(SG_string)//' '
    else
     SG_string = trim(SG_string)//'~'
    end if
    cycle
   endif
  end if

  alpha_car   = .false.
  numeric_car = .false.
  call check_character(string(i:i), alpha_car, numeric_car)
  if(alpha_car) then
   if(HTML_report) then
    SG_string = trim(SG_string)//"<i>"//string(i:i)//"</i>"
   elseif(LATEX_report) then
    SG_string = trim(SG_string)//"\textit{"//string(i:i)//"}"
   endif
   !if(i==1) then
   ! if(HTML_report) then
   !  SG_string = trim(SG_string)//"&nbsp;"
   ! elseif(LATEX_report) then
   !  if(latexTT) then
   !   SG_string = trim(SG_string)//' '
   !  else
   !   SG_string = trim(SG_string)//'~'
   !  endif
   ! end if
   !end if
  else
   if(i==2 .and. LATEX_report .and. latexTT) then
    SG_string = trim(SG_string)//' '//string(i:i)
   else
    SG_string = trim(SG_string)//string(i:i)
   end if
  end if
 end do

 string = SG_string

  return
end subroutine Check_letter_string

!------------------------------------------------------------------------------------------------------
subroutine Latex_preambule
 use LATEX_module,    only : logo
 use cryscalc_module, only : LATEX_unit, archive_cif, cryscalc, report_logo, debug_proc

 implicit none
  character (len=256), dimension(2)   :: logo_full
  logical,             dimension(2)   :: file_exist
  character (len=256)                 :: DOS_command
  integer                             :: i, nb_files_exist, num_file_exist

  if(debug_proc%level_2)  call write_debug_proc_level(2, "latex_preambule")

  logo_full = ''
  !logo(1) = '\img\CDIFX_logo.jpg'
  !logo(2) = '\img\ISCR_logo.jpg'

 WRITE(LATEX_unit, '(a)')  "% Structural report (LATEX file) created by CRYSCALC.exe"
 WRITE(LATEX_unit, '(2a)') "%      Starting CIF file : ", trim(archive_CIF)
 WRITE(LATEX_unit, '(a)')  "%"
 WRITE(LATEX_unit, '(a)')  "%     Web site : http://www.cdifx.univ-rennes/cryscalc.htm"
 WRITE(LATEX_unit, '(5a)') "%     Version : ", trim(cryscalc%version), " [", trim(cryscalc%author), "]"
 WRITE(LATEX_unit, '(a)')  "%"

 write(LATEX_unit, '(a)') '\documentclass[11pt , a4paper]{article}'
 write(LATEX_unit, '(a)') '%% ///////   PREAMBULE  ////// %%'
 write(LATEX_unit, '(a)') '\usepackage[francais,english]{babel} % adaptation de LaTeX à la langue française'
 write(LATEX_unit, '(a)') '% geometrie de la page'
 write(LATEX_unit, '(a)') '\usepackage[dvips,lmargin=2.5cm,rmargin=2.5cm,tmargin=2.5cm,bmargin=2.5cm]{geometry}'
 write(LATEX_unit, '(a)') '\usepackage{color}                 % utilisation des couleurs'
 write(LATEX_unit, '(a)') "\usepackage{graphicx}              % insertion d'images"
 write(LATEX_unit, '(a)') '\usepackage{fancybox}              % fonctions encadrement'
 write(LATEX_unit, '(a)') '\renewcommand{\thefootnote}{\*}    % supprime le numero de la footnote'
 write(LATEX_unit, '(a)') "\usepackage[flushmargin]{footmisc} % supprime l'indentation de la footnote"
 write(LATEX_unit, '(a)') "\usepackage[pdftex,colorlinks=true,urlcolor=gris50,pdfstartview=FitH]{hyperref}   % hyperliens"
 write(LATEX_unit, '(a)') ''
 write(LATEX_unit, '(a)') "\usepackage{alltt}  % permet l'ecriture mathematique a l'interieur de l'environnement verbatim"
 write(LATEX_unit, '(a)') "\usepackage[T1]{fontenc}"
 write(LATEX_unit, '(a)') '%%%%%%%%%% defintion de couleurs %%%%'
 write(LATEX_unit, '(a)') '\definecolor{gris50}{gray}{0.50}'
 write(LATEX_unit, '(a)') '\definecolor{gris75}{gray}{0.75}'
 write(LATEX_unit, '(a)') '\definecolor{gris_clair} {rgb}{0.96078, 0.96078, 0.96078}    % rgb(245,245,245)  #f5f5f5'
 write(LATEX_unit, '(a)') '\definecolor{violet}     {rgb}{0.5,     0,       0.5}        % rgb{128,0,128)    #800080'
 write(LATEX_unit, '(a)') '\definecolor{red_25b}    {rgb}{0.5,     0,       0.1}        % rgb{127,0,25)     #7F0019'
 write(LATEX_unit, '(a)') '\definecolor{vert_titre} {rgb}{0.95294, 0.98039, 0.90196}    % rgb(243,250,230)  #F3FAE6'
 write(LATEX_unit, '(a)') '\definecolor{bleu_titre} {rgb}{0.843,   0.890,   0.929}      % rgb(215,227,237)  #D7E3ED'
 write(LATEX_unit, '(a)') '\definecolor{fonce_titre}{rgb}{0.28125, 0.22266, 0.38672}    % rgb(72, 57,99)    #483963'
 write(LATEX_unit, '(a)') '%%%%%%%%%%%%% macros %%%%%%%%%%%%%%%%'
 write(LATEX_unit, '(a)') '\newcommand{\bs}   {$\backslash$}'
 write(LATEX_unit, '(a)') '\newcommand{\header} [1]{'
 write(LATEX_unit, '(a)') '\tiny'
 write(LATEX_unit, '(a)') '\noindent'
 write(LATEX_unit, '(a)') '\textcolor{gris75}{#1}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') ''
 write(LATEX_unit, '(a)') '% titre'
 write(LATEX_unit, '(a)') '\newcommand{\titre}  [2] {'
 write(LATEX_unit, '(a)') ' \begin{center}'
 write(LATEX_unit, '(a)') '\fcolorbox{black}{bleu_titre}'
 write(LATEX_unit, '(a)') '{'
 write(LATEX_unit, '(a)') '\parbox{8cm}{'
 write(LATEX_unit, '(a)') '\begin{center}'
 write(LATEX_unit, '(a)') '\textcolor{fonce_titre}{\texttt{\textbf{#1}}}'
 write(LATEX_unit, '(a)') '\textcolor{fonce_titre}{\textit{#2}}'
 write(LATEX_unit, '(a)') '\end{center}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '\vspace{0.5cm}'
 write(LATEX_unit, '(a)') '\end{center}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '%'
 write(LATEX_unit, '(a)') '% sous-titre'
 write(LATEX_unit, '(a)') '\newcommand{\SousTitre}  [1] {'
 write(LATEX_unit, '(a)') '\vspace{0.5cm}'
 write(LATEX_unit, '(a)') '\noindent'
 write(LATEX_unit, '(a)') '\fcolorbox{gris_clair}{gris_clair}'
 write(LATEX_unit, '(a)') '{'
 write(LATEX_unit, '(a)') '\parbox{14.5cm}{'
 write(LATEX_unit, '(a)') '\textcolor{violet}{\hspace{0.0cm}#1}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '\vspace{0.5cm}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '%'
 write(LATEX_unit, '(a)') '%% ///////   FIN DU PREAMBULE  ////// %%'
 write(LATEX_unit, '(a)') ''
 write(LATEX_unit, '(a)') '\begin{document}'
 write(LATEX_unit, '(a)') '\renewcommand*\familydefault{\sfdefault} %% Only if the base font of the document is to be sans serif'
 write(LATEX_unit, '(a)') '\normalfont'
 write(LATEX_unit, '(a)') ''

 if(len_trim(cryscalc%path_name) /= 0) then
  !file_exist = .false.
  !write(logo_full(1), '(3a)') trim(cryscalc%path_name), '\img\', trim(logo(1))
  !inquire(file = trim(logo_full(1)), exist = file_exist(1))
  !write(logo_full(2), '(3a)') trim(cryscalc%path_name), '\img\', trim(logo(2))
  !inquire(file = trim(logo_full(2)), exist = file_exist(2))

  !if(file_exist(1) .and. file_exist(2)) then
  ! ! les fichiers doivent etre dans le repertoire de travail
  ! write(DOS_command, '(3a)') 'copy ', trim(logo_full(1)), ' .'
  ! call system(trim(DOS_command))
!
!   write(DOS_command, '(3a)') 'copy ', trim(logo_full(2)), ' .'
!   call system(trim(DOS_command))
!   write(LATEX_unit, '(a)') ''
!   write(LATEX_unit, '(a)') '\begin{figure}[h]'
!   !write(LATEX_unit, '(a)') '% \includegraphics[width=50.pt]{img/ISCR_logo.png} \includegraphics[width=40.pt]{img/cdifx_logo.jpg}'
!   write(LATEX_unit, '(5a)') ' \includegraphics[height=30.pt]{', trim(logo(1)), '}  \includegraphics[height=30.pt]{', &
!                             trim(logo(2)), '}'
!   write(LATEX_unit, '(a)') '\end{figure}'
!   write(LATEX_unit, '(a)') ''
!  endif

  nb_files_exist = 0
  do i=1, 2
   file_exist = .false.
   write(logo_full(i), '(3a)') trim(cryscalc%path_name), '\img\', trim(report_logo(i))
   inquire(file = trim(logo_full(i)), exist = file_exist(i))
   if(file_exist(i)) then
    nb_files_exist = nb_files_exist + 1
    num_file_exist = i
    ! les fichiers doivent etre dans le repertoire de travail
    write(DOS_command, '(3a)') 'copy ', trim(logo_full(i)), ' .'
    call system(trim(DOS_command))
   end if
  end do

  if(nb_files_exist/=0) then
   write(LATEX_unit, '(a)') ''
   write(LATEX_unit, '(a)') '\begin{figure}[h]'
   if(nb_files_exist == 2) then
    write(LATEX_unit, '(5a)') ' \includegraphics[height=30.pt]{', trim(report_logo(1)), '}  \includegraphics[height=30.pt]{', &
                              trim(report_logo(2)), '}'
   else
    write(LATEX_unit, '(3a)') ' \includegraphics[height=30.pt]{', trim(logo(num_file_exist)), '}'
   end if
   write(LATEX_unit, '(a)') '\end{figure}'
   write(LATEX_unit, '(a)') ''
  endif

 endif

 return
end subroutine Latex_preambule

!----------------------------------------------------------------------------------------------------------------------------
subroutine Latex_End
 use LATEX_module,    only : logo
 use cryscalc_module, only : LATEX_unit, cryscalc, report_logo, debug_proc

 implicit none

 if(debug_proc%level_2)  call write_debug_proc_level(2, "latex_end")

 write(LATEX_unit, '(a)')''
 write(LATEX_unit, '(6a)') '\footnote{\noindent \textcolor{gris50}{This structural report has been created through ', &
                           '\href{http://', trim(CRYSCALC%url), '}{CRYSCALC} (', trim(cryscalc%version), ').'

 write(LATEX_unit, '(5a)') 'Please report bugs and problems to \href{mailto:', &
                           trim(cryscalc%mail), "}{", trim(CRYSCALC%mail), "}.}}"

 write(LATEX_unit, '(a)') '\end{document}'

 !call system('del '//report_trim(logo(1)))
 !call system('del '//report_trim(logo(2)))

 return
end subroutine Latex_End

!------------------------------------------------------------------------------------
subroutine check_LATEX_file_name(LATEX_string)
 use macros_module, only : replace_car
 implicit none
 character (len=*), intent(inout)  :: LATEX_string

  LATEX_string = replace_car(LATEX_string, '\',   '\bs ')
  LATEX_string = replace_car(LATEX_string, '_',   '\_')

 return
end subroutine check_LATEX_file_name


!------------------------------------------------------------------------------------
subroutine report_crystal_study(diffracto)
 use structural_report_module, only : job
 use cryscalc_module,          only : HTML_unit, LATEX_unit, HTML_report, LATEX_report, LATEX_unit, debug_proc
 use CIF_module,               only : CIF_parameter
 use special_characters_module
 implicit none
  character (len=*), intent(in)    :: diffracto    ! APEX, KCCD, X2S
  CHARACTER (LEN=512)              :: LATEX_string
  real                             :: T_min, T_max

  if(debug_proc%level_2)  call write_debug_proc_level(2, "report_crystal_study")

  if(len_trim(diffracto) == 3 .and. diffracto(1:3) == "X2S") then
   !if(len_trim(doc_type) == 4 .and. doc_type(1:4) == "HTML") then
   if(HTML_report) then
    WRITE(HTML_unit, '(a)')   "<p class='retrait_1'>"
    WRITE(HTML_unit, '(a)')   ""
    write(HTML_unit, '(5a)')  "A ", trim(CIF_parameter%crystal_colour), " crystal of ", trim(job), " with the dimensions of "
    write(HTML_unit, '(10a)') TRIM(CIF_parameter%crystal_size_max), " mm x ", TRIM(CIF_parameter%crystal_size_mid), " mm x ", &
         TRIM(CIF_parameter%crystal_size_min), " mm mounted on a Mitegen Micromount was automatically centered on a Bruker ", &
         "SMART X2S benchtop crystallographic system. Intensity measurements were performed using monochromated (doubly ",    &
         "curved silicon crystal) Mo-K&alpha; radiation (&lambda;=0.71073 Å) from a sealed microfocus tube. Generator ",      &
         "settings were 50 kV, 1 mA. Data collection temperature was 21°C. Data were acquired using three sets of &Omega; ",  &
         "scans at different &phi; settings."
         ! The frame width was 0.5° with an exposure time of 5.0 s."
    write(HTML_unit, '(a)')  ""
    write(HTML_unit, '(a)')  "The detailed data collection stategy was as follows:<br>"
    write(HTML_unit, '(a)')  "&nbsp;&nbsp;  Detector distance : 40 mm<br>"
    write(HTML_unit, '(a)')  "&nbsp;&nbsp;  Detector swing angle (fixed 2&theta;) : -20°"
    write(HTML_unit, '(a)')  "<pre>"
    write(HTML_unit, '(10x,a)')  "  Run       &Omega;(start)     &Omega;(end)         &Phi;     Frames"
    write(HTML_unit, '(10x,a)')  "    1           -20.       160.       0.0        360"
    write(HTML_unit, '(10x,a)')  "    2           -20.       100.     120.0        240"
    write(HTML_unit, '(10x,a)')  "    3           -20.        40.     240.0        120"
    write(HTML_unit, '(10x,a)')  "</pre>"
    write(HTML_unit, '(a)')  ""
    WRITE(HTML_unit, '(a)')   "<p class='retrait_1'>"
    write(HTML_unit, '(7a)') "APEX2 [1] software was used for preliminary determination of the unit cell. Determination of ",    &
                             "integrated intensities and unit cell refinement were performed using SAINT [2]. The integration ", &
                             "of the data yielded a total of ", TRIM(CIF_parameter%diffrn_reflns_number), " reflections to a ",  &
                             "maximum &theta; angle of ", TRIM(CIF_parameter%theta_max), " °."
    write(HTML_unit, '(3a)') "The constants for the ", TRIM(CIF_parameter%symmetry_cell_setting), " unit cell are "

    call report_cell_parameters("HTML")

    WRITE(HTML_unit, '(2a)') "They are based upon the refinement of the XYZ-centroid of ", trim(CIF_parameter%cell_reflns_used)
    WRITE(HTML_unit, '(5a)') " reflections above 20.0 I/&sigma;(I) with ",trim(CIF_parameter%cell_theta_min), " <= ", &
                            trim(CIF_parameter%cell_theta_max), '.'
    if(CIF_parameter%T_min /= '?' .and. CIF_parameter%T_max /= '?') then
     read(CIF_parameter%T_min, *) T_min
     read(CIF_parameter%T_max, *) T_max
    else
     T_min = 1.
     T_max = 1.
    endif
    WRITE(HTML_unit, '(2a,F5.1,6a)') "Data were corrected for absorption effects with SADABS [3] using the multiscan technique. ", &
                                     "The ratio of minimum apparent transmission is ", 100.*T_min/T_max, " :100. The average ",    &
                                     "residual for symmetry equivalent reflections is R<sub>int</sub> = ",                         &
                                     TRIM(CIF_parameter%diffrn_reflns_av_R_equivalents), " and R&sigma; = " ,                      &
                                     TRIM(CIF_parameter%diffrn_reflns_av_R_sigma), "."
    WRITE(HTML_unit, '(7a)') "XPREP [4] determined the space goup to be ", trim(CIF_parameter%symmetry_space_group), " with Z = ", &
                             trim(CIF_parameter%formula_units_Z), " for the fomula unit ", trim(CIF_parameter%formula_moiety), "."

    WRITE(HTML_unit, '(21a)') "The structure was solved with XS [5] and subsequent structure refinements were performed with ",  &
                              "XL [6]. The final anisotropic full-matrix least-squares refinements on F<sup>2</sup> with ",      &
                              trim(CIF_parameter%refine_ls_number_parameters), " variables converged at R<sub>1</sub> = ",       &
                              TRIM(CIF_parameter%refine_ls_R_factor_gt),  " for the observed data and wR<sub>2</sub> = ",        &
                              TRIM(CIF_parameter%refine_ls_wR_factor_gt), " for all data. The goodness-of-fit was ",             &
                              TRIM(CIF_parameter%Chi2), ". The largest peak on the final difference electron density synthesis ",&
                              "was ", TRIM(CIF_parameter%refine_diff_density_max),  trim(HTML_elect_density_unit), " and the ",  &
                              "deepest hole was ", TRIM(CIF_parameter%refine_diff_density_min), trim(HTML_elect_density_unit),   &
                              " with an RMS deviation of ", trim(CIF_parameter%refine_diff_density_rms),                         &
                              trim(HTML_elect_density_unit), "."
    WRITE(HTML_unit, '(a)')   "<p class='retrait_1'>"
    WRITE(HTML_unit, '(a)')   "[1] APEX2 version 2009.9 (Bruker AXS Inc.)<br>"
    WRITE(HTML_unit, '(a)')   "[2] SAINT version 7.68 A (Bruker AXS Inc., 2009)<br>"
    WRITE(HTML_unit, '(a)')   "[3] SADABS version 2008/1 (Sheldrick, Bruker AXS Inc.)<br>"
    WRITE(HTML_unit, '(a)')   "[4] XPREP version 2008/2 (Sheldrick, Bruker AXS Inc.)<br>"
    WRITE(HTML_unit, '(a)')   "[5] XS version 2008/1 (George M. Sheldrick, <i>Acta Cryst.</i> (2008), A<b>64</b>,112-122)<br>"
    WRITE(HTML_unit, '(a)')   "[6] XL version 2008/4 (George M. Sheldrick, <i>Acta Cryst.</i> (2008), A<b>64</b>,112-122)<br>"
    WRITE(HTML_unit, '(a)')   ""
   endif


   if(LATEX_report) then

    WRITE(LATEX_unit, '(a)')   ""
    write(LATEX_unit, '(5a)')  "A ", trim(CIF_parameter%crystal_colour), " crystal of ", trim(job), " with the dimensions of "
    write(LATEX_unit, '(10a)') TRIM(CIF_parameter%crystal_size_max), " mm x ", TRIM(CIF_parameter%crystal_size_mid), " mm x ",    &
              TRIM(CIF_parameter%crystal_size_min), "mm mounted on a Mitegen Micromount was automatically centered on a Bruker ", &
              "SMART X2S benchtop crystallographic system. Intensity measurements were performed using monochromated (doubly ",   &
          "curved silicon crystal) $Mo-K_{\alpha}$ radiation ($\lambda$ =0.71073 \AA) from a sealed microfocus tube. Generator ", &
       "settings were 50 kV, 1 mA. Data collection temperature was 21\degre C. Data were acquired using three sets of $\Omega$ ", &
       "scans at different $\phi$ settings."
       ! The frame width was 0.5\degre with an exposure time of 5.0 s."
    write(LATEX_unit, '(a)')  ""
    write(LATEX_unit, '(a)')  "The detailed data collection stategy was as follows:\\"
    write(LATEX_unit, '(a)')  "\hspace*{1.cm} Detector distance : 40 mm\\"
    write(LATEX_unit, '(a)')  "\hspace*{1.cm} Detector swing angle (fixed $2\theta$) : -20\degre"
    write(LATEX_unit, '(a)')  "\begin{alltt}"
    write(LATEX_unit, '(10x,a)')  "  Run      \(\Omega\)(start)     \(\Omega\)(end)         \(\phi\)     Frames"
    write(LATEX_unit, '(10x,a)')  "    1           -20.       160.       0.0        360"
    write(LATEX_unit, '(10x,a)')  "    2           -20.       100.     120.0        240"
    write(LATEX_unit, '(10x,a)')  "    3           -20.        40.     240.0        120"
    write(LATEX_unit, '(10x,a)')  "\end{alltt}"
    write(LATEX_unit, '(7a)') "APEX2 [1] software was used for preliminary determination of the unit cell. Determination of ",    &
                              "integrated intensities and unit cell refinement were performed using SAINT [2]. The integration ", &
                              "of the data yielded a total of ", TRIM(CIF_parameter%diffrn_reflns_number), " reflections to a ",  &
                              "maximum $\theta$ angle of ", TRIM(CIF_parameter%theta_max), " °.\\"
    write(LATEX_unit, '(3a)') "The constants for the ", TRIM(CIF_parameter%symmetry_cell_setting), " unit cell are "

    call report_cell_parameters("LATEX")

    WRITE(LATEX_unit, '(2a)') "They are based upon the refinement of the XYZ-centroid of ", trim(CIF_parameter%cell_reflns_used)
    WRITE(LATEX_unit, '(5a)') " reflections above 20.0 $I/\sigma(I)$ with ",trim(CIF_parameter%cell_theta_min), " $\le$ ", &
                            trim(CIF_parameter%cell_theta_max), '.\\'
    if(CIF_parameter%T_min /= '?' .and. CIF_parameter%T_max /= '?') then
     read(CIF_parameter%T_min, *) T_min
     read(CIF_parameter%T_max, *) T_max
    else
     T_min = 1.
     T_max = 1.
    endif

    WRITE(LATEX_unit,'(2a,F4.1,6a)') "Data were corrected for absorption effects with SADABS [3] using the multiscan technique. ", &
                                     "The ratio of minimum apparent transmission is ", 100.*T_min/T_max, " :100. The average ",    &
                                     "residual for symmetry equivalent reflections is $R_{int}$ = ",                         &
                                     TRIM(CIF_parameter%diffrn_reflns_av_R_equivalents), " and $R_{\sigma}$ = " ,            &
                                     TRIM(CIF_parameter%diffrn_reflns_av_R_sigma), "."

    LATEX_string = CIF_parameter%symmetry_space_group
    call check_indice_string(LATEX_string)
    WRITE(LATEX_unit, '(3a)') "XPREP [4] determined the space goup to be $", trim(LATEX_string), "$ with $Z = "
    CALL transf_moiety_string("LATEX", CIF_parameter%formula_moiety , LATEX_string)
    WRITE(LATEX_unit, '(4a)')  trim(CIF_parameter%formula_units_Z), "$ for the fomula unit $", trim(LATEX_string), "$.\\"
    WRITE(LATEX_unit, '(21a)') "The structure was solved with XS [5] and subsequent structure refinements were performed with ",  &
                               "XL [6]. The final anisotropic full-matrix least-squares refinements on $F^2$ with ",      &
                               trim(CIF_parameter%refine_ls_number_parameters), " variables converged at $R_1$ = ",       &
                                 TRIM(CIF_parameter%refine_ls_R_factor_gt),  " for the observed data and $\omega R_2$ = ",        &
                               TRIM(CIF_parameter%refine_ls_wR_factor_gt), " for all data. The goodness-of-fit was ",             &
                               TRIM(CIF_parameter%Chi2), ". The largest peak on the final difference electron density synthesis ",&
                               "was ", TRIM(CIF_parameter%refine_diff_density_max),  trim(LATEX_elect_density_unit), " and the ", &
                                "deepest hole was ", TRIM(CIF_parameter%refine_diff_density_min), trim(LATEX_elect_density_unit), &
                               " with an RMS deviation of ", trim(CIF_parameter%refine_diff_density_rms),                         &
                               trim(LATEX_elect_density_unit), "."
    WRITE(LATEX_unit, '(a)')   "\vspace{0.5cm}"
    WRITE(LATEX_unit, '(2a)')  "\begin{enumerate}"
    WRITE(LATEX_unit, '(2a)')  "\item ",  "APEX2 version 2009.9 (Bruker AXS Inc.)"
    WRITE(LATEX_unit, '(2a)')  "\item ",  "SAINT version 7.68 A (Bruker AXS Inc., 2009)"
    WRITE(LATEX_unit, '(2a)')  "\item ",  "SADABS version 2008/1 (Sheldrick, Bruker AXS Inc.)"
    WRITE(LATEX_unit, '(2a)')  "\item ",  "XPREP version 2008/2 (Sheldrick, Bruker AXS Inc.)"
    WRITE(LATEX_unit, '(2a)')  "\item ",  "XS version 2008/1 (George M. Sheldrick, Acta Cryst. (2008), A64,112-122)"
    WRITE(LATEX_unit, '(2a)')  "\item ",  "XL version 2008/4 (George M. Sheldrick, Acta Cryst. (2008), A64,112-122)"
    WRITE(LATEX_unit, '(a)')   "\end{enumerate}"
   endif
  end if


 return
end subroutine report_crystal_study

!------------------------------------------------------------------------------------
subroutine report_cell_parameters(doc_type)
 use CRYSCALC_module, only  : HTML_unit, LATEX_unit, debug_proc
 use CIF_module,      only  : CIF_parameter
 use macros_module,   only  : l_case
 implicit none
  character (len=*), intent(in)     :: doc_type

  if(debug_proc%level_3)  call write_debug_proc_level(3, "report_cell_parameters")


  if(len_trim(doc_type) == 4 .and. doc_type(1:4) == "HTML") then
   IF(l_case(CIF_parameter%symmetry_cell_setting(1:9)) == 'triclinic') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
     WRITE(HTML_unit, '(3a)') "&alpha; = ",TRIM(CIF_parameter%cell_angle_alpha), ", "
     WRITE(HTML_unit, '(3a)') "&beta; = ",TRIM(CIF_parameter%cell_angle_beta),  ", "
     WRITE(HTML_unit, '(3a)') "&gamma; = ",TRIM(CIF_parameter%cell_angle_gamma), " °, "
   ELSEIF(l_case(CIF_parameter%symmetry_cell_setting(1:10)) == 'monoclinic') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
     WRITE(HTML_unit, '(3a)') "&beta; = ",TRIM(CIF_parameter%cell_angle_beta), " °, "
   ELSEIF(l_case(CIF_parameter%symmetry_cell_setting(1:12)) == 'orthorhombic') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
   ELSEIF(l_case(CIF_parameter%symmetry_cell_setting(1:10)) == 'tetragonal') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
   ELSEIF(l_case(CIF_parameter%symmetry_cell_setting(1:9)) == 'hexagonal' .or.    &
          l_case(CIF_parameter%symmetry_cell_setting(1:8)) == 'trigonal') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
   ELSEIF(l_case(CIF_parameter%symmetry_cell_setting(1:11)) == 'rhomboedral') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), " Å, "
     WRITE(HTML_unit, '(3a)') "&alpha; = ",TRIM(CIF_parameter%cell_angle_alpha), " °, "
   ELSEIF(l_case(CIF_parameter%symmetry_cell_setting(1:5)) == 'cubic') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_c), " Å, "
   ENDIF
   WRITE(HTML_unit, '(3a)', advance='NO') "<i>V</i> = ",TRIM(CIF_parameter%cell_volume), " Å<sup>3</sup>."

   return
  ENDIF

  if(len_trim(doc_type) == 5 .and. doc_type(1:5) == "LATEX") then
   IF(CIF_parameter%symmetry_cell_setting(1:9) == 'triclinic') then
    WRITE(LATEX_unit, '(3a)') "$a = ",TRIM(CIF_parameter%cell_length_a), ", "
    WRITE(LATEX_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
    WRITE(LATEX_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), "$ \AA, "
    WRITE(LATEX_unit, '(3a)') "$\alpha; = ",TRIM(CIF_parameter%cell_angle_alpha), ", "
    WRITE(LATEX_unit, '(3a)') "\beta; = ",TRIM(CIF_parameter%cell_angle_beta),  ", "
    WRITE(LATEX_unit, '(3a)') "\gamma; = ",TRIM(CIF_parameter%cell_angle_gamma), "$ \degre, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:10) == 'monoclinic') then
    WRITE(LATEX_unit, '(3a)') "$a = ",TRIM(CIF_parameter%cell_length_a), ", "
    WRITE(LATEX_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
    WRITE(LATEX_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), "$ \AA, "
    WRITE(LATEX_unit, '(3a)') "$\beta; = ",TRIM(CIF_parameter%cell_angle_beta), " $\degre, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:12) == 'orthorhombic') then
    WRITE(LATEX_unit, '(3a)') "$a = ",TRIM(CIF_parameter%cell_length_a), ", "
    WRITE(LATEX_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
    WRITE(LATEX_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), "$ \AA, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:10) == 'tetragonal') then
    WRITE(LATEX_unit, '(3a)') "$a = ",TRIM(CIF_parameter%cell_length_a), ", "
    WRITE(LATEX_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), "$ \AA, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:9) == 'hexagonal' .or.    &
          CIF_parameter%symmetry_cell_setting(1:8) == 'trigonal') then
    WRITE(LATEX_unit, '(3a)') "$a = ",TRIM(CIF_parameter%cell_length_a), ", "
    WRITE(LATEX_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), "$ \AA, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:11) == 'rhomboedral') then
    WRITE(LATEX_unit, '(3a)') "$a = ",TRIM(CIF_parameter%cell_length_a), "$ \AA, "
    WRITE(LATEX_unit, '(3a)') "$\alpha; = ",TRIM(CIF_parameter%cell_angle_alpha), "$ \degre, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:5) == 'cubic') then
    WRITE(LATEX_unit, '(3a)') "$a = ",TRIM(CIF_parameter%cell_length_c), "$ \AA, "
   ENDIF
   WRITE(LATEX_unit, '(3a)', advance='NO') "$V = ",TRIM(CIF_parameter%cell_volume), "$ \AA$^3$, "

   return
  ENDIF


 return
end subroutine report_cell_parameters








