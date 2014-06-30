
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

  character (len=32), parameter         :: HTML_elect_density_unit      = " e<sup>-</sup>.�<sup>-3</sup>"
  character (len=32), parameter         :: LATEX_elect_density_unit     = " $e^-.\AA^{-3}$"
  character (len=32), parameter         :: LATEX_alt_elect_density_unit = "  \(\overline{e}\).\AA\(\sp{-3}\)"
  

end module special_characters_module
!---------------------------------------------------------------------------------------------------------------

subroutine create_structural_report
 USE cryscal_module, only : debug_proc
 
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
 USE structural_report_module,     ONLY : job, op_string
 USE MACROS_module,                ONLY : test_file_exist,  l_case                                            
 USE cryscal_module,               ONLY : CIF_unit, CIF_parameter, archive_CIF, debug_proc, include_squeeze, debug_proc
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
   call write_debug_proc(' CRYSCAL will be stopped.', '')
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
		call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%symmetry_cell_setting)		

       case ('_symmetry_space_group_name_h-m', '_space_group_name_h-m', '_space_group_name_h-m_alt')
		call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%symmetry_space_group)		

       case ('_symmetry_int_tables_number', '_space_group_it_number')
		call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%symmetry_IT_number)		

       case ('_symmetry_equiv_pos_site_id')
	    symm_id = .true.
		
       case ('_symmetry_equiv_pos_as_xyz')
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
		  op_string(n_op) = CIF_input_line(index(CIF_input_line, ' '):)
		 endif
		 op_string(n_op) = adjustl(op_string(n_op))
         op_string(n_op) = trim(op_string(n_op))
         if(op_string(n_op)(1:1) == "'") then
		  op_string(n_op) = op_string(n_op)(2:len_trim(op_string(n_op))-1) ! enleve les caracteres "'"
		 endif
        end do
		

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
 		call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%diffracto_device)		
        IF(CIF_parameter%diffracto_device(1:6) == 'APEXII') then
         device_mark =  'Bruker-AXS'
         device_type =  'APEXII Kappa-CCD diffractometer'
        elseif(CIF_parameter%diffracto_device(1:3) == 'X2S') then
         device_mark = 'Bruker-AXS'
         device_type = 'SMART X2S benchtop diffractometer'        
		elseif(CIF_parameter%diffracto_device(1:8) == 'KappaCCD') then
         device_mark = 'Nonius'
         device_type = 'KappaCCD'
        ELSEIF(CIF_parameter%diffracto_device(1:11) == 'CCD Saphire') then
         device_mark = 'Oxford Diffraction'
         device_type = 'Xcalibur Saphire 3'
        endif


       case ('_diffrn_ambient_temperature')
 		call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%diffracto_temperature)		

       case ('_diffrn_radiation_type')
 		call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%diffracto_radiation_type)		

       case ('_diffrn_radiation_source')
 		call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%diffracto_radiation_source)		

       case ('_diffrn_radiation_wavelength')
 		call Get_CIF_value(CIF_unit, CIF_field_value, CIF_parameter%diffracto_radiation_wavelength)		

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
  USE cryscal_module,               ONLY : HTML_unit, text_unit, LATEX_unit, archive_CIF,  text_report, HTML_report, & 
                                           LATEX_report, cryscal_ini, cryscal_version, cryscal_author, author, report_header, &
                                           cryscal_report_css, debug_proc
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
  WRITE(HTML_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.HTML' 
  open (UNIT=HTML_unit, FILE=TRIM(HTML_structural_report_file))
 endif
 
 if(text_report) then
  WRITE(TEXT_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.TXT' 
  open (UNIT=text_unit, FILE=TRIM(TEXT_structural_report_file))
 endif
 
 if(latex_report) then
  WRITE(LATEX_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.LTX' 
  open (UNIT=latex_unit, FILE=TRIM(LATEX_structural_report_file))
  WRITE(PDF_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.pdf'   
 endif
 
 
 if(HTML_report) then
  WRITE(HTML_unit, '(a)')  "<!-- Structural report (HTML file) created by CRYSCAL.exe"
  WRITE(HTML_unit, '(2a)')  "      Starting CIF file : ", trim(archive_CIF)
  !WRITE(HTML_unit, '(2a)') "      CRYSCAL setting file : ", trim(cryscal_ini)
  WRITE(HTML_unit, '(a)')  ""
  WRITE(HTML_unit, '(a)')  "     Web site : http://www.cdifx.univ-rennes/cryscal.htm"
  WRITE(HTML_unit, '(5a)') "     Version : ", trim(cryscal_version), " [", trim(cryscal_author), "]"
  WRITE(HTML_unit, '(a)') ""
  if(len_trim(cryscal_report_css) ==0) then
   WRITE(HTML_unit, '(a)')  "     Default css style"
  else
   WRITE(HTML_unit, '(2a)') "     css file : ", trim(cryscal_report_css)
  endif
  WRITE(HTML_unit, '(a)') "!-->"
  
  WRITE(HTML_unit, '(a)') "<HTML>"
  WRITE(HTML_unit, '(a)') "<HEAD>"
  call write_HTML_css(HTML_unit, 'report')
  WRITE(HTML_unit, '(a)') "<TITLE>Structural report</TITLE>"
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
   WRITE(HTML_unit, '(5a)') "CRYSCAL version : ", trim(cryscal_version), " [", trim(cryscal_author), "]</p>"
   WRITE(HTML_unit, '(a)')  "<hr>" 
  endif 
  if(LATEX_report) then
   LATEX_string = wingx_structure_dir
   call check_LATEX_file_name(LATEX_string)   
   WRITE(LATEX_unit, '(3a)') '\header{Working directory: ', trim(LATEX_string), '} \\'
   LATEX_string = archive_cif
   call check_LATEX_file_name(LATEX_string)   
   WRITE(LATEX_unit, '(3a)') "\header{Input CIF file : ", trim(LATEX_string), '} \\'
   WRITE(LATEX_unit, '(5a)') "\header{CRYSCAL version : ", trim(cryscal_version), " [", trim(cryscal_author), "]} \\"
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
                                          check_character, Get_wingx_job,  Get_current_folder_name, replace_car, multiple
 USE cryscal_module,               ONLY : CIF_unit, HTML_unit, text_unit, LATEX_unit, CIF_parameter, archive_CIF, long_report, &
                                          CIF_CELL_measurement, CIF_torsion_limit, text_report, HTML_report, LATEX_report,     &
                                          SQUEEZE, my_browser, my_pdflatex, cryscal_ini, message_text, debug_proc,             &
										  cryscal_version, cryscal_author, structure_solution, structure_refinement, author,   &
										  report_header, include_SQUEEZE, debug_proc
 USE external_applications_module, ONLY : launch_browser, launch_word, launch_pdflatex
 USE IO_module
 USE LATEX_module,                 ONLY : logo
 USE Special_characters_module


 implicit none 
 LOGICAL                             :: file_exist
 CHARACTER (LEN=256)                 :: CIF_input_line, CIF_field, CIF_field_value
 CHARACTER (LEN=256)                 :: HTML_string, new_HTML_string
 CHARACTER (LEN=256), dimension(4)   :: report_string
 CHARACTER (LEN=256)                 :: LATEX_string
 CHARACTER (LEN=256)                 :: DOS_command
 CHARACTER (LEN=32)                  :: GIF_file, PNG_file
 INTEGER                             :: i_error, long_field, long, i, i1, i2 
 integer                             :: n, n_op, op_numor
 
 integer,   dimension(500)           :: op_n       ! numero de l'op. de sym.
 integer,   dimension(500,3)         :: op_t       ! partie translation
 integer,   dimension(3)             :: t
 integer                             :: n_sym_htab, n_sym_dist, n_sym_ang, n_sym_tor
 integer                             :: nb_col
 CHARACTER (LEN=64)                  :: device_mark, device_type

 CHARACTER (LEN=12), dimension(500) :: site_sym
 CHARACTER (len=12), dimension(500) :: dico
 INTEGER           , dimension(500) :: num_site_sym

 CHARACTER (LEN=12), dimension(2)    :: dist_atom
 CHARACTER (LEN=12)                  :: dist_value
 CHARACTER (LEN=12)                  :: dist_sym

 CHARACTER (LEN=12), dimension(3)    :: ang_atom
 CHARACTER (LEN=12)                  :: ang_value
 CHARACTER (LEN=12), dimension(2)    :: ang_sym

 CHARACTER (LEN=12), dimension(4)    :: torsion_ang_atom
 CHARACTER (LEN=12)                  :: torsion_ang_value 
 CHARACTER (LEN=12), dimension(4)    :: torsion_sym
 real                                :: torsion_value_real


 CHARACTER (LEN=12)                  :: atom_label, atom_typ, atom_x, atom_y, atom_z, atom_Ueq, atom_adp_type, atom_occ
 CHARACTER (LEN=12)                  :: atom_U11, atom_U22, atom_U33, atom_U23, atom_U13, atom_U12
 CHARACTER (len=12)                  :: site_label_D, site_label_H, site_label_A
 CHARACTER (len=12)                  :: dist_DH, dist_HA, dist_DA, angle_DHA, site_sym_A

 LOGICAL                             :: alpha_car, num_car
 CHARACTER (LEN=10)                  :: date, time  
 character (len=256)                 :: wingx_structure_dir
 character (len=10)                  :: AUTHOR_initiales

 real                                :: T_min, T_max


 if(debug_proc%level_2)  call write_debug_proc_level(2, "create_report")

 n_sym_htab = 0
 n_sym_dist = 0
 n_sym_ang  = 0
 n_sym_tor  = 0


 
 if(HTML_report) then
  WRITE(HTML_unit, '(a)')  "<br>"    
  WRITE(HTML_unit, '(3a)') "<p class='title_main'><i>", trim(job),"</i>: Crystal structure report</p><br>"
  WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;X-ray crystallographic study</p>"
 endif
 
 if(LATEX_report) then 
  write(LATEX_unit, '(a)') '\footnotesize'
  WRITE(LATEX_unit, '(a)')  ""
  !job = replace_car(job, '_', '/')
  !job = replace_car(job, "/", "\_")
  job = replace_car(job, '_', '\_')
  WRITE(LATEX_unit, '(4a)')  '\titre{', trim(job), ': Crystal structure report','}'  
  WRITE(LATEX_unit, '(a)')   '\SousTitre{X-ray crystallographic study}'  
 endif
  
  
  if(CIF_parameter%formula_moiety == '?' .or. len_trim(CIF_parameter%formula_moiety) ==0)  &
     CIF_parameter%formula_moiety = CIF_parameter%formula_sum

  if(index(CIF_parameter%diffrn_measurement_device_type, 'X2S') /=0) then
   if(HTML_report .or. LATEX_report) then
 	call report_crystal_study("X2S")   
   endif	
  else
  
	 
  if(HTML_report) then
   CALL transf_moiety_string("HTML", CIF_parameter%formula_moiety , HTML_string)
   WRITE(HTML_unit, '(a)') "<p class='retrait_1'>"
   WRITE(HTML_unit, '(a)')  ""
   !WRITE(HTML_unit, '(5a)', advance='NO') "(",TRIM(CIF_parameter%formula_moiety),"); M = ",TRIM(CIF_parameter%formula_weight),"."
   WRITE(HTML_unit, '(5a)', advance='NO') "(",TRIM(HTML_string),"); <i>M</i> = ",TRIM(CIF_parameter%formula_weight),". "
   long = len_trim(CIF_parameter%diffrn_measurement_device_type)
   if(CIF_parameter%diffrn_measurement_device_type(1:1)== "'"           .and. &
      CIF_parameter%diffrn_measurement_device_type(long:long)== "'") then
    WRITE(HTML_unit, '(4a)') CIF_parameter%diffrn_measurement_device_type(2:long-1), " diffractometer,"
   else
    WRITE(HTML_unit, '(4a)') TRIM(CIF_parameter%diffrn_measurement_device_type), " diffractometer,"
   endif
  endif

  if(LATEX_report) then  
   CALL transf_moiety_string("LATEX", CIF_parameter%formula_moiety , LATEX_string)
   !LATEX_string = CIF_parameter%formula_moiety
   WRITE(LATEX_unit, '(a)')  ""
   WRITE(LATEX_unit, '(a)')  "{\setlength{\baselineskip}{1.2\baselineskip}"   
   WRITE(LATEX_unit, '(5a)', advance='NO') "($",TRIM(LATEX_string),"$); $M = ",TRIM(CIF_parameter%formula_weight),"$. "
   long = len_trim(CIF_parameter%diffrn_measurement_device_type)
   if(CIF_parameter%diffrn_measurement_device_type(1:1)== "'"           .and. &
      CIF_parameter%diffrn_measurement_device_type(long:long)== "'") then
    WRITE(LATEX_unit, '(4a)') CIF_parameter%diffrn_measurement_device_type(2:long-1), " diffractometer,"
   else
    WRITE(LATEX_unit, '(4a)') TRIM(CIF_parameter%diffrn_measurement_device_type), " diffractometer,"
   endif
  endif

  
  IF(CIF_parameter%diffracto_radiation_type(1:5) == 'MoK\a') then
   if(HTML_report) then  
    WRITE(HTML_unit, '(a)', advance='NO') "Mo-K&alpha; radiation "
   endif
   if(LATEX_report) then  
    WRITE(LATEX_unit, '(a)', advance='NO') "$Mo-K_{\alpha}$ radiation "
   endif   
  else
   if(HTML_report) then
    WRITE(HTML_unit, '(a)', advance='NO') TRIM(CIF_parameter%diffracto_radiation_type), " radiation"
   endif
   if(LATEX_report) then
    WRITE(LATEX_unit, '(a)', advance='NO') TRIM(CIF_parameter%diffracto_radiation_type), " radiation"
   endif	
  endif
  if(HTML_report) then
   WRITE(HTML_unit, '(3a)') "(&lambda; = ", TRIM(CIF_parameter%diffracto_radiation_wavelength)," �),"
   WRITE(HTML_unit, '(3a)', advance='NO') "<i>T</i> = ", TRIM(CIF_parameter%diffracto_temperature), " K; "
  endif 
  if(LATEX_report) then
   WRITE(LATEX_unit, '(3a)') "($\lambda$ = ", TRIM(CIF_parameter%diffracto_radiation_wavelength)," \AA),"
   WRITE(LATEX_unit, '(3a)', advance='NO') "$T = ", TRIM(CIF_parameter%diffracto_temperature), " K$; "
  endif 

  alpha_car    = .false.
  num_car      = .false.
  HTML_string  = ''
  LATEX_string = ''
  long = len_trim(HTML_string)
  
  if(HTML_report) then   
   HTML_string = CIF_parameter%symmetry_space_group   
   call check_indice_string(HTML_string) 
   WRITE(HTML_unit, '(6a)') TRIM(CIF_parameter%symmetry_cell_setting), " <i>", trim(HTML_string),  &
                           "</i> (I.T.#", TRIM(CIF_parameter%symmetry_IT_number), "), "
  endif
  
  if(LATEX_report) then
   LATEX_string = CIF_parameter%symmetry_space_group
   call check_indice_string(LATEX_string)
   WRITE(LATEX_unit, '(6a)') TRIM(CIF_parameter%symmetry_cell_setting), " $", trim(LATEX_string),  &
                           "$ (I.T.\#", TRIM(CIF_parameter%symmetry_IT_number), "), "
  endif
  
						   
  if(HTML_report) then
   call report_cell_parameters("HTML")
   WRITE(HTML_unit, '(3a)', advance='NO') "<i>Z</i> = ",TRIM(CIF_parameter%formula_units_Z), ", "
   WRITE(HTML_unit, '(3a)') "<i>d</i> = ",TRIM(CIF_parameter%exptl_density), " g.cm<sup>-3</sup>, "
   WRITE(HTML_unit, '(3a)') "&mu; = ",TRIM(CIF_parameter%exptl_mu), " mm<sup>-1</sup>. "
  endif 

  if(LATEX_report) then
   call report_cell_parameters("LATEX")
   
   WRITE(LATEX_unit, '(3a)', advance='NO') "$Z = ",TRIM(CIF_parameter%formula_units_Z), "$, "
   WRITE(LATEX_unit, '(3a)') "$d = ",TRIM(CIF_parameter%exptl_density), "~g.cm^{-3}$, "
   WRITE(LATEX_unit, '(3a)') "$\mu = ",TRIM(CIF_parameter%exptl_mu), "~mm^{-1}$. "
  endif  

  
  !SQUEEZE%procedure = .false.
  !if(include_SQUEEZE) then
  !inquire(file= trim(SQUEEZE%file), exist= file_exist)
  !if(file_exist) then
  ! SQUEEZE%procedure = .true.
  !else   
  ! inquire(file= "sqz.sqf", exist= file_exist)
  ! if(file_exist)  then
  !  SQUEEZE%procedure = .true. 
  !  SQUEEZE%file = "sqz.sqf"
  ! else
  !  inquire(file="platon_sqr.sqf", exist = file_exist)  ! cree par PLATON version jan. 2013
  !	if(file_exist) then
  !	 SQUEEZE%procedure != .true. 
  !   SQUEEZE%file = "platon_sqr.sqf"
  !	endif
  ! endif	
  !endif
  !end if
  
  if(HTML_report) then
   !WRITE(HTML_unit, '(a)') "The structure was solved by direct methods using the SIR97 program [1], "
   WRITE(HTML_unit, '(3a)') "The structure was solved by direct methods using the <i>",  &
                             u_case(trim(structure_solution%name)), "</i> program [1], "
   WRITE(HTML_unit, '(3a)') "and then refined with full-matrix least-square methods based on <i>F</i><sup>2</sup> (<i>", &
                             u_case(trim(structure_refinement%name)), "</i>) [2] "
   WRITE(HTML_unit, '(a)')  "with the aid of the <i>WINGX</i> [3] program."

   !if(SQUEEZE%procedure) then
   if(include_SQUEEZE) then
    WRITE(HTML_unit, '(a)') "The contribution of the disordered solvents to the calculated structure factors was "
    WRITE(HTML_unit, '(3a)') "estimated following the <i>BYPASS</i> algorithm [4], implemented as the <i>SQUEEZE</i> ", &
                             "option in <i>PLATON</i> [5]. A new data set, free of solvent contribution, was then ",    &
                            "used in the final refinement."
   endif

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
    WRITE(LATEX_unit, '(3a)') "The structure was solved by direct methods using the ",  &
                              u_case(trim(structure_solution%name)), " program [1], "
    WRITE(LATEX_unit, '(3a)') "and then refined with full-matrix least-square methods based on $F^2$ (", &
                              u_case(trim(structure_refinement%name)), ") [2] with the aid of the WINGX [3] program."
    !if(SQUEEZE%procedure) then
	if(include_SQUEEZE) then
     WRITE(LATEX_unit, '(a)') "The contribution of the disordered solvents to the calculated structure factors was "
     WRITE(LATEX_unit, '(3a)') "estimated following the BYPASS algorithm [4], implemented as the SQUEEZE ", &
                               "option in PLATON [5]. A new data set, free of solvent contribution, was then ",    &
                               "used in the final refinement."
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
   !WRITE(HTML_unit, '(a)') "[3] &nbsp; L. J. Farrugia, J. Appl. Cryst., 1999, 32, 837-838<br>"
   WRITE(HTML_unit, '(a)') "[3] &nbsp; L. J. Farrugia, J. Appl. Cryst., 2012, 45, 849-854<br>"
   !if(SQUEEZE%procedure) then
   if(include_SQUEEZE) then
    WRITE(HTML_unit, '(a)') "[4] &nbsp; P. v.d. Sluis and A.L. Spek, Acta Cryst. (1990) A46, 194-201<br>"
    WRITE(HTML_unit, '(a)') "[5] &nbsp; A. L. Spek, J. Appl. Cryst. (2003), 36, 7-13<br>"
   end if
   WRITE(HTML_unit, '(a)') "</p>"
   WRITE(HTML_unit, '(a)') "<br>"
  endif
  if(LATEX_report) then
   WRITE(LATEX_unit, '(a)') ""
   WRITE(LATEX_unit, '(a)') "\begin{enumerate}"
   WRITE(LATEX_unit, '(2a)') "\item ", trim(structure_solution%reference)
   WRITE(LATEX_unit, '(2a)') "\item ", trim(structure_refinement%reference)
   !WRITE(LATEX_unit, '(a)')  "\item  L. J. Farrugia, J. Appl. Cryst., 1999, 32, 837-838"
   WRITE(LATEX_unit, '(a)')  "\item   L. J. Farrugia, J. Appl. Cryst., 2012, 45, 849-854"
   if(SQUEEZE%procedure) then
    WRITE(LATEX_unit, '(a)') "\item P. v.d. Sluis and A.L. Spek, Acta Cryst. (1990) A46, 194-201"
    WRITE(LATEX_unit, '(a)') "\item A. L. Spek, J. Appl. Cryst. (2003), 36, 7-13"
   end if
   WRITE(LATEX_unit, '(a)') "\end{enumerate}"
   WRITE(LATEX_unit, '(a)') ""
  endif
  end if   ! fin de la condition sur le type de diffracto. utilise

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
   WRITE(HTML_unit, '(10x,3a)') "Temperature                            ", TRIM(CIF_parameter%diffracto_temperature)," <i>K</i>"
   WRITE(HTML_unit, '(10x,3a)') "Wavelength                             ", TRIM(CIF_parameter%diffracto_radiation_wavelength), " � "
   HTML_string = CIF_parameter%symmetry_space_group
   call check_indice_string(HTML_string)
   
   WRITE(HTML_unit, '(10x,6a)') "Crystal system, space group            ", TRIM(CIF_parameter%symmetry_cell_setting), ", ", &
                                "<i>",TRIM(HTML_string),"</i>"							   
   WRITE(HTML_unit, '(10x,6a)') "Unit cell dimensions                   ", "a = ", TRIM(CIF_parameter%cell_length_a), &
                                " �, &alpha; = ", TRIM(CIF_parameter%cell_angle_alpha), " �"
   WRITE(HTML_unit, '(10x,6a)') "                                       ", "b = ", TRIM(CIF_parameter%cell_length_b), &
                                " �, &beta; = ", TRIM(CIF_parameter%cell_angle_beta),  " �"
   WRITE(HTML_unit, '(10x,6a)') "                                       ", "c = ", TRIM(CIF_parameter%cell_length_c), &
                                " �, &gamma; = ", TRIM(CIF_parameter%cell_angle_gamma), " �"
   WRITE(HTML_unit, '(10x,3a)') "Volume                                 ", TRIM(CIF_parameter%cell_volume), " �<sup>3</sup>"
   WRITE(HTML_unit, '(10x,5a)') "<i>Z</i>, Calculated density                  ", TRIM(CIF_parameter%formula_units_Z), ', ', &
                                TRIM(CIF_parameter%exptl_density), " (g.cm<sup>-3</sup>)"
   WRITE(HTML_unit, '(10x,3a)') "Absorption coefficient                 ", TRIM(CIF_parameter%exptl_mu), " mm<sup>-1</sup>"
   IF(long_report) then
    WRITE(HTML_unit, '(10x,2a)') "<i>F</i>(000)                                 ", TRIM(CIF_parameter%F000)
    WRITE(HTML_unit, '(10x,7a)') "Crystal size                           ", TRIM(CIF_parameter%crystal_size_max), " x ", &
                               TRIM(CIF_parameter%crystal_size_mid), " x ", TRIM(CIF_parameter%crystal_size_min), " mm"
    WRITE(HTML_unit, '(10x,2a)') "Crystal color                          ", TRIM(CIF_parameter%crystal_colour)
    WRITE(HTML_unit, '(10x,5a)') "Theta range for data collection        ", TRIM(CIF_parameter%theta_min), " to ", &
                                TRIM(CIF_parameter%theta_max), " �"
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
    WRITE(HTML_unit, '(10x,2a)') "<sup>b</sup>Goodness-of-fit                       ", TRIM(CIF_parameter%Chi2)
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
   WRITE(LATEX_unit, '(10x,3a)') "Temperature                            ", TRIM(CIF_parameter%diffracto_temperature)," K"
   WRITE(LATEX_unit, '(10x,3a)') "Wavelength                             ", TRIM(CIF_parameter%diffracto_radiation_wavelength), &
                                 " \AA "
   LATEX_string = CIF_parameter%symmetry_space_group
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
   WRITE(LATEX_unit, '(10x,5a)') "Z, Calculated density                  ", TRIM(CIF_parameter%formula_units_Z), ' , ', &
                                TRIM(CIF_parameter%exptl_density), " (\(g.cm\sp{-1}\))"
   WRITE(LATEX_unit, '(10x,3a)') "Absorption coefficient                 ", TRIM(CIF_parameter%exptl_mu), " \(mm\sp{-1}\)"
   
   WRITE(LATEX_unit, '(10x,2a)') "F(000)                                 ", TRIM(CIF_parameter%F000)
   WRITE(LATEX_unit, '(10x,7a)') "Crystal size                           ", TRIM(CIF_parameter%crystal_size_max), " x ", &
                               TRIM(CIF_parameter%crystal_size_mid), " x ", TRIM(CIF_parameter%crystal_size_min), " mm"
   WRITE(LATEX_unit, '(10x,2a)') "Crystal color                          ", TRIM(CIF_parameter%crystal_colour)
   WRITE(LATEX_unit, '(10x,5a)') "Theta range for data collection        ", TRIM(CIF_parameter%theta_min), " to ", &
                                TRIM(CIF_parameter%theta_max), " \degre"
   WRITE(LATEX_unit, '(10x,4a)') "h_min, h_max                           ", TRIM(CIF_parameter%h_min), " , ", &
                                TRIM(CIF_parameter%h_max)
   WRITE(LATEX_unit, '(10x,4a)') "k_min, k_max                           ", TRIM(CIF_parameter%k_min), " , ", &
                                TRIM(CIF_parameter%k_max)
   WRITE(LATEX_unit, '(10x,4a)') "l_min, l_max                           ", TRIM(CIF_parameter%l_min), " , ", &
                                TRIM(CIF_parameter%l_max)
   WRITE(LATEX_unit, '(10x,7a)') "Reflections collected / unique         ", TRIM(CIF_parameter%diffrn_reflns_number), &
                               " / ", TRIM(CIF_parameter%reflns_number_total), " [\(\sp{a}\)R(int) = ",                &
                               TRIM(CIF_parameter%diffrn_reflns_av_R_equivalents), "]"
   WRITE(LATEX_unit, '(10x,2a)') "Reflections [\(I>2\sigma\)]                   ",	CIF_parameter%reflns_number_gt
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
   write(text_unit, '(2x,3a)') 'Table 1. Crystal data and structure refinement for ', trim(job),'.'
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,2a)') "Empirical formula                      ", TRIM(CIF_parameter%formula_sum)
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,2a)') "Formula weight                         ", TRIM(CIF_parameter%formula_weight)
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,3a)') "Temperature                            ", TRIM(CIF_parameter%diffracto_temperature)," K"
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,3a)') "Wavelength                             ", TRIM(CIF_parameter%diffracto_radiation_wavelength), " � "
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,4a)') "Crystal system, space group            ", TRIM(CIF_parameter%symmetry_cell_setting), ", ", &
                                TRIM(CIF_parameter%symmetry_space_group)
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,6a)') "Unit cell dimensions                   ", "a = ", TRIM(CIF_parameter%cell_length_a), &
                               " �, alpha = ", TRIM(CIF_parameter%cell_angle_alpha), " �"
   WRITE(text_unit, '(10x,6a)') "                                       ", "b = ", TRIM(CIF_parameter%cell_length_b), &
                               " �, beta = ", TRIM(CIF_parameter%cell_angle_beta),  " �"
   WRITE(text_unit, '(10x,6a)') "                                       ", "c = ", TRIM(CIF_parameter%cell_length_c), &
                               " �, gamma = ", TRIM(CIF_parameter%cell_angle_gamma), " �"
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,3a)') "Volume                                 ", TRIM(CIF_parameter%cell_volume), " �3"
  
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,5a)') "Z, Calculated density                  ", TRIM(CIF_parameter%formula_units_Z), ' , ', &
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
                                TRIM(CIF_parameter%theta_max), " �"
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,4a)') "h_min, h_max                           ", TRIM(CIF_parameter%h_min), " , ", &
                                TRIM(CIF_parameter%h_max)
   WRITE(text_unit, '(10x,4a)') "k_min, k_max                           ", TRIM(CIF_parameter%k_min), " , ", &
                                TRIM(CIF_parameter%k_max)
   WRITE(text_unit, '(10x,4a)') "l_min, l_max                           ", TRIM(CIF_parameter%l_min), " , ", &
                                TRIM(CIF_parameter%l_max)
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,7a)') "Reflections collected / unique         ", TRIM(CIF_parameter%diffrn_reflns_number), &
                               " / ", TRIM(CIF_parameter%reflns_number_total), " [R(int) = ",                       &
                               TRIM(CIF_parameter%diffrn_reflns_av_R_equivalents), "]"
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,2a)') "Reflections [I>2sigma(I)]              ",	CIF_parameter%reflns_number_gt
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
  
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,5a)') "Final R indices [I>2sigma(I)]          ", "R1a = ",  &
                                TRIM(CIF_parameter%refine_ls_R_factor_gt),  ", wR2b = ",          &
                                TRIM(CIF_parameter%refine_ls_wR_factor_gt)
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,5a)') "R indices (all data)                   ", "R1a = ",        &
                               TRIM(CIF_parameter%refine_ls_R_factor_all), ", wR2b = ",           &
                               TRIM(CIF_parameter%refine_ls_wR_factor_ref)
   write(text_unit, '(a)') ''	
   WRITE(text_unit, '(10x,5a)') "Largest diff. peak and hole            ", TRIM(CIF_parameter%refine_diff_density_max),     &
                               " and ", TRIM(CIF_parameter%refine_diff_density_min), "  e.�-3"
   WRITE(text_unit, '(a)')      ""   
  endif
  
  
  

!--------------------------------------------------------

  if(HTML_report) then
   WRITE(HTML_unit, '(a)') "<br>"
   ! positions atomiques
   WRITE(HTML_unit, '(a)')  "<p class='title_2'>"
   WRITE(HTML_unit, '(2a)') "&nbsp;&nbsp;Atomic coordinates, site occupancy (%) and equivalent isotropic displacement ", &
                            " parameters (�<sup>2</sup> x 10<sup>3</sup>). "
   WRITE(HTML_unit, '(a)')  "U(eq) is defined as one third of the trace of the orthogonalized U<sub>ij</sub> tensor."
   WRITE(HTML_unit, '(a)')  "</p>"
   WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
  endif
  
  if(text_report) then
   WRITE(text_unit, '(a)')  ""
   WRITE(text_unit, '(a)')  ""
   WRITE(text_unit, '(2x,a)')  "Table 2. Atomic coordinates, site occupancy (%) and equivalent isotropic displacement " 
   WRITE(text_unit, '(2x,3a)') "         parameters (�^2 x 10^3) for ", trim(job), "."
   WRITE(text_unit, '(2x,a)')  "         U(eq) is defined as one third of the trace of the orthogonalized Uij tensor."
   WRITE(text_unit, '(a)')  ""
  endif
  
  if(LATEX_report) then   ! positions atomiques
   WRITE(LATEX_unit, '(a)')   "\newpage"
   WRITE(LATEX_unit, '(3a)')  "\SousTitre{Atomic coordinates, site occupancy (\%) and equivalent isotropic displacement ", &
                              "parameters ($A^2$ x $10^3$). U(eq) is defined as one third of the trace of the ",           &
						 	  "orthogonalized $U_{ij}$ tensor.}"
   WRITE(LATEX_unit, '(a)')  ""
   WRITE(LATEX_unit, '(a)') "\begin{verbatim}"
  endif


!---------------------------------------------------------------------
   


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
           WRITE(HTML_unit, '(10x,a)') "Atom           x              y              z              occ.               U(eq)"
           WRITE(HTML_unit, '(a)')     ""
		  end if 
		  if(text_report) then
 		   WRITE(text_unit, '(5x,100a1)')  ("-",i=1,100)
           WRITE(text_unit, '(10x,a)') "Atom           x              y              z              occ.               U(eq)"
		   WRITE(text_unit, '(5x,100a1)')  ("-",i=1,100)
           WRITE(text_unit, '(a)')     ""
		  endif
		  if(LATEX_report) then
		   WRITE(LATEX_unit, '(5x,a)') "Atom           x              y              z              occ.               U(eq)"
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
		   if(HTML_report) then
            WRITE(HTML_unit, '(10x,6(a,3x))') atom_label(1:12), atom_x(1:12), atom_y(1:12), atom_z(1:12), atom_occ(1:12), &
		                                      atom_Ueq(1:12)
           endif
           if(text_report) then
		    WRITE(text_unit, '(10x,6(a,3x))') atom_label(1:12), atom_x(1:12), atom_y(1:12), atom_z(1:12), atom_occ(1:12), &
		                                     atom_Ueq(1:12)
           endif
           if(LATEX_report) then
		    WRITE(LATEX_unit, '(5x,6(a,3x))') atom_label(1:12), atom_x(1:12), atom_y(1:12), atom_z(1:12), atom_occ(1:12), &
		                                      atom_Ueq(1:12)
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
             WRITE(HTML_unit, '(a)')  "<p class='title_2'>&nbsp;&nbsp;Anisotropic displacement parameters (�<sup>2</sup>)</p>"
             WRITE(HTML_unit, '(3a)') "<p class='retrait_1'>The anisotropic displacement factor exponent takes the form: ", &
                                      "-2&pi;<sup>2</sup> [ h<sup>2</sup> a*<sup>2</sup> U<sub>11</sub> + ... ", &
                                      "+ 2 h k a* b* U<sub>12</sub> ]. </p>"             
             WRITE(HTML_unit, '(10x,2a)') "Atom           U11            U22            U33            ",   &
                                          "U23            U13            U12"
             WRITE(HTML_unit, '(a)')  ""			 
			endif 
			if(text_report) then
			 WRITE(text_unit, '(a)') ""
			 WRITE(text_unit, '(a)') ""
			 WRITE(text_unit, '(2x,3a)') "Table 3. Anisotropic displacement parameters (�^2) for ",trim(job),"."
             WRITE(text_unit, '(2x,a)')  "         The anisotropic displacement factor exponent takes the form:"
             WRITE(text_unit, '(2x,a)')  "         -2pi^2 [ h^2 a*^2 U11> + ...  + 2 h k a* b* U12 ]."
             WRITE(text_unit, '(a)')  ""
			 WRITE(text_unit, '(5x,110a1)')  ("-",i=1,110)
             WRITE(text_unit, '(10x,2a)') "Atom           U11            U22            U33            ",   &
                                          "U23            U13            U12"
			 WRITE(text_unit, '(5x,110a1)')  ("-",i=1,110)
             WRITE(text_unit, '(a)')  ""
			endif
			if(LATEX_report) then
			 WRITE(LATEX_unit, '(2a)') "The anisotropic displacement factor exponent takes the form:", &
                                       "$-2\pi [h^2 a^{*2} U_{11} + ... + 2hka^*b^*U_{12}]$"
             WRITE(LATEX_unit, '(a)')  ""
             WRITE(LATEX_unit, '(a)')  "\begin{verbatim}"		
             WRITE(LATEX_unit, '(5x,2a)') "Atom           U11            U22            U33            ",   &
                                          "U23            U13            U12"
             WRITE(LATEX_unit, '(a)')  ""			 
			endif
           endif
           CIF_parameter%atom = CIF_input_line
           READ(CIF_parameter%atom, *) atom_label, atom_U11, Atom_U22, atom_U33, atom_U23, atom_U13, atom_U12
           if(HTML_report) then
		    WRITE(HTML_unit, '(10x,7(a,3x))')  atom_label(1:12), atom_U11(1:12), atom_U22(1:12), atom_U33(1:12), &
                                               atom_U23(1:12), atom_U13(1:12), atom_U12(1:12)
           endif											   
		   if(text_report) then
            WRITE(text_unit, '(10x,7(a,3x))')  atom_label(1:12), atom_U11(1:12), atom_U22(1:12), atom_U33(1:12), &
                                               atom_U23(1:12), atom_U13(1:12), atom_U12(1:12)
	       endif				   

		   if(LATEX_report) then
            WRITE(LATEX_unit, '(5x,7(a,3x))')  atom_label(1:5), atom_U11(1:12), atom_U22(1:12), atom_U33(1:12), &
                                               atom_U23(1:12), atom_U13(1:12), atom_U12(1:12)
	       endif				   

          END do


         case default
     end select
    end do
    close (UNIT = CIF_unit)
	if(HTML_report   .and. n /=0)   WRITE(HTML_unit, '(a)')  "</pre>"
    if(LATEX_report  .and. n /=0)  WRITE(LATEX_unit, '(a)')  "\end{verbatim}"

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
           WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Bond lengths [�]</p>"
           WRITE(HTML_unit, '(a)')  "<pre style='font-size:14' class='color_grey1'>"
		  endif 
		  if(text_report) then
		   WRITE(text_unit, '(a)') ""
		   WRITE(text_unit, '(a)') ""
		   WRITE(text_unit, '(2x,3a)') "Table 4. Bond lengths [�] for ", trim(job), '.'
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
   if(LATEX_report .and. n /=0) then
    WRITE(LATEX_unit, '(a)') "\end{verbatim} "
   endif	

   if(n_sym_dist /=0) then
    if(HTML_report) then
     write(HTML_unit, '(a)')  ""
     write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
     WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey'>"
	endif 
    if(text_report) then
	 write(text_unit, '(a)')  ""
     write(text_unit, '(a)')  " Symmetry transformations used to generate equivalent atoms:"
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
         if(n==1) then
		  if(HTML_report) then
           WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Angles [�]</p>"
           WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
		  endif 
          if(text_report) then
           WRITE(text_unit, '(a)') ""
	       WRITE(text_unit, '(2x,3a)') "Table 5. Angles [deg] for ", trim(job), "."
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
		  write(HTML_string, '(10x,7a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
		                                 report_string(2)(1:8), ' = ', trim(ang_value)
          !if(n == 2*int(n/2)) then
		  if(multiple(n,2)) then
		   write(HTML_unit, '(3a)') '<span class="ligne_1">',HTML_string(1:70),'</span>'
		  else	
		   write(HTML_unit, '(3a)') '<span class="ligne_2">',HTML_string(1:70),'</span>'
		  end if	
         endif
         if(text_report) then
          write(text_unit, '(10x,7a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
		                               report_string(2)(1:8), ' = ', trim(ang_value)
         endif		 
		 if(LATEX_report) then
		  write(LATEX_unit, '(5x,7a)') report_string(1)(1:8), ' - ', ang_atom(2)(1:8), ' - ', &
		                               report_string(2)(1:8), ' = ', trim(ang_value)
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
   if(LATEX_report .and. n /=0) then
    WRITE(LATEX_unit, '(a)') "\end{verbatim}"
   endif

  if(n_sym_ang /=0) then
    WRITE(HTML_unit, '(a)')  "</pre>"
    write(HTML_unit, '(a)')  ""
    write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
    WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
    if(text_report) then
	 write(text_unit, '(a)')  ""
     write(text_unit, '(a)')  " Symmetry transformations used to generate equivalent atoms:"
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
           WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Torsion angles [�]</p>"         
           WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
		  endif 
		  if(text_report) then
           WRITE(text_unit, '(a)') ""
	       WRITE(text_unit, '(2x,3a)') "Table 6. Torsion angles [deg] for ", trim(job), '.'
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
		  read(torsion_ang_value(i1+1:i2-1), *) torsion_value_real
		 else
		  read(torsion_ang_value, *) torsion_value_real
		 endif 
		 if(torsion_value_real > CIF_torsion_limit) then
		  n = n -1
		  cycle
		 endif
		 
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
		  WRITE(HTML_string, '(10x, 9a)') report_string(1)(1:8), ' - ', &
		                                  report_string(2)(1:8), ' - ', &
		 			                      report_string(3)(1:8), ' - ', &
	                 		              report_string(4)(1:8), ' = ', trim(torsion_ang_value)
		  !if(n == 2*int(n/2)) then
		  if(multiple(n,2)) then
		   write(HTML_unit, '(3a)') '<span class="ligne_1">',HTML_string(1:70),'</span>'
		  else	
		   write(HTML_unit, '(3a)') '<span class="ligne_2">',HTML_string(1:70),'</span>'
		  end if	
         endif										
		 if(text_report) then
          WRITE(text_unit, '(10x, 9a)') report_string(1)(1:8), ' - ', &
		                                report_string(2)(1:8), ' - ', &
								        report_string(3)(1:8), ' - ', &
									    report_string(4)(1:8), ' = ', trim(torsion_ang_value)
		 endif		 
		 if(LATEX_report) then
          WRITE(LATEX_unit, '(5x, 9a)') report_string(1)(1:8), ' - ', &
		                                report_string(2)(1:8), ' - ', &
								        report_string(3)(1:8), ' - ', &
							            report_string(4)(1:8), ' = ', trim(torsion_ang_value)
		 endif		 
		end do


       case default
     end select
    END do
   close (UNIT = CIF_unit)
   if(HTML_report  .and. n /=0)    WRITE(HTML_unit, '(a)')   "</pre><br>"
   if(LATEX_report .and. n /=0)    WRITE(LATEX_unit, '(a)')  "\end{verbatim}"

  if(n_sym_tor /=0) then
    if(HTML_report) then
     write(HTML_unit, '(a)')  ""
     write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
     WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
	endif 
    if(text_report) then
	 write(text_unit, '(a)')  ""
     write(text_unit, '(a)')  " Symmetry transformations used to generate equivalent atoms:"
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
           WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Hydrogen bonds [� and deg.] </p>"
           WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
		  endif 
		  if(text_report) then
           WRITE(text_unit, '(a)') ""
		   WRITE(text_unit, '(a)') ""
	       WRITE(text_unit, '(2x,3a)') "Table 7. Hydrogen bonds [A and deg.] for ", trim(job), '.'
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
     write(text_unit, '(a)')  " Symmetry transformations used to generate equivalent atoms:"
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
 
   file_exist = .false.
  do  
   i1 = index(CIF_parameter%diffracto_temperature, "(")  
   if(i1 /=0) then
    write(GIF_file, '(4a)') trim(job), "_", CIF_parameter%diffracto_temperature(1:i1-1), "K_ortep.gif"  
   else	
    write(GIF_file, '(4a)') trim(job), "_", TRIM(CIF_parameter%diffracto_temperature), "K_ortep.gif"  
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
  
  


  if(HTML_report) then
   WRITE(HTML_unit, '(a)')  "<HR>"
   WRITE(HTML_unit, '(a)')  "<p class='retrait_1'>"
   WRITE(HTML_unit, '(4a)') "<i>This HTML report has been created through ",                        &
                            "<A HREF='http://www.cdifx.univ-rennes1.fr/cryscal.htm'>CRYSCAL</A> (", &
                            trim(cryscal_version), ")."
   WRITE(HTML_unit, '(2a)') "Please report bugs and problems to <A HREF='mailto:cdifx@univ-rennes1.fr'>", &
                            "cdifx@univ-rennes1.fr</A></p>"
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
   write(message_text, '(3a)') '  No browser defined in the ', trim(cryscal_ini), ' setting file !!'
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
   
   call system('del '//trim(logo(1)))
   call system('del '//trim(logo(2)))

  else
   call write_info('')
   write(message_text, '(3a)') '  No LATEX to PDF convertor defined in the ', trim(cryscal_ini), ' setting file !!'
   call write_info(trim(message_text))
   call write_info('')
  endif
 endif
 

  return

end subroutine create_report



!---------------------------------------------------------------------
subroutine get_sym_op(input_string, numor_op, t)
 use cryscal_module, only : message_text
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
    read(input_string, *) numor_op   ! decomment� oct 2011
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
 use macros_module,  only : nombre_de_colonnes, check_character
 use cryscal_module, only : debug_proc
 implicit none
  character (len=*),   intent(in)     :: code             !(HTML/CIF/HTML)
  character (len=*),   intent(inout)  :: input_string
  character (len=256), intent(out)    :: output_string
  character (len=256)                 :: tmp_string
  character (len=8), dimension(2)     :: balise
  integer                             :: i, j, i1, long, nb_col, num, num_alpha 
  character (len=8),  dimension(32)   :: part_string
  LOGICAL                             :: alpha_char, numeric_char
  character (len=1)                   :: espace
  character (len=8)                   :: espace_2

 if(debug_proc%level_3)  call write_debug_proc_level(3, "transf_moiety_string ("//trim(input_string)//")")


  output_string = ''
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
	   write(OUTPUT_string, '(3a)') trim(OUTPUT_string), trim(espace_2) , part_string(i)(j:j)  ! espaces	   
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
	   write(OUTPUT_string, '(3a)') trim(OUTPUT_string), trim(espace_2), part_string(i)(j:j)	  	   
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
!
subroutine  Get_num_site_sym_new(n_sym, site_sym, current_sym, num_site_sym, dico)
 implicit none
  integer,                             intent(in)    :: n_sym
  CHARACTER (LEN=12), dimension(500),  intent(in)    :: site_sym
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
 use cryscal_module, only : HTML_unit, text_unit, LATEX_unit, HTML_report, text_report, LATEX_report, message_text
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



!--------------------------------------------------------------------------
subroutine check_indice_string(string)
 ! met un caractere en indice si necessaire (chaine = groupe d'espace)
 use cryscal_module, only : HTML_report, LATEX_report
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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine Latex_preambule
 use LATEX_module,   only : logo
 use cryscal_module, only : LATEX_unit, archive_cif, cryscal_version, cryscal_author, cryscal_path_name, debug_proc
 
 implicit none
  character (len=256), dimension(2)   :: logo_full
  logical,             dimension(2)   :: file_exist
  character (len=256)                 :: DOS_command
  
  if(debug_proc%level_2)  call write_debug_proc_level(2, "latex_preambule")

  logo_full = ''
  logo(1) = 'CDIFX_logo.jpg'
  logo(2) = 'ISCR_logo.jpg'
 
 WRITE(LATEX_unit, '(a)')  "% Structural report (LATEX file) created by CRYSCAL.exe"
 WRITE(LATEX_unit, '(2a)') "%      Starting CIF file : ", trim(archive_CIF)
 WRITE(LATEX_unit, '(a)')  "%"
 WRITE(LATEX_unit, '(a)')  "%     Web site : http://www.cdifx.univ-rennes/cryscal.htm"
 WRITE(LATEX_unit, '(5a)') "%     Version : ", trim(cryscal_version), " [", trim(cryscal_author), "]"
 WRITE(LATEX_unit, '(a)')  "%"

 write(LATEX_unit, '(a)') '\documentclass[11pt , a4paper]{article}'
 write(LATEX_unit, '(a)') '%% ///////   PREAMBULE  ////// %%'
 write(LATEX_unit, '(a)') '\usepackage[francais,english]{babel} % adaptation de LaTeX � la langue fran�aise'
 write(LATEX_unit, '(a)') '% geometrie de la page'
 write(LATEX_unit, '(a)') '\usepackage[dvips,lmargin=2.5cm,rmargin=2.5cm,tmargin=2.5cm,bmargin=2.5cm]{geometry}'
 write(LATEX_unit, '(a)') '\usepackage{color}               % utilisation des couleurs'
 write(LATEX_unit, '(a)') "\usepackage{graphicx}            % insertion d'images"
 write(LATEX_unit, '(a)') '\usepackage{fancybox}            % fonctions encadrement'
 write(LATEX_unit, '(a)') '\renewcommand{\thefootnote}{\*}  % supprime le numero de la footnote'
 write(LATEX_unit, '(a)') "\usepackage[pdftex,colorlinks=true,urlcolor=gris50,pdfstartview=FitH]{hyperref}   % hyperliens"
 write(LATEX_unit, '(a)') ''
 write(LATEX_unit, '(a)') "\usepackage{alltt}  % permet l'ecriture mathematique a l'interieur de l'environnement verbatim"
 write(LATEX_unit, '(a)') '%%%%%%%%%% defintion de couleurs %%%%'
 write(LATEX_unit, '(a)') '\definecolor{gris50}{gray}{0.50}'
 write(LATEX_unit, '(a)') '\definecolor{gris75}{gray}{0.75}'
 write(LATEX_unit, '(a)') '\definecolor{gris_clair} {rgb}{0.96078, 0.96078, 0.96078}    % rgb(245,245,245)  #f5f5f5'
 write(LATEX_unit, '(a)') '\definecolor{violet}     {rgb}{0.5,     0,       0.5}        % rgb{128,0,128)    #800080'
 write(LATEX_unit, '(a)') '\definecolor{red_25b}    {rgb}{0.5,     0,       0.1}        % rgb{127,0,25)     #7F0019'
 write(LATEX_unit, '(a)') '\definecolor{vert_titre} {rgb}{0.95294, 0.98039, 0.90196}    % rgb(243,250,230)  #F3FAE6'
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
 write(LATEX_unit, '(a)') '\newcommand{\titre}  [1] {'
 write(LATEX_unit, '(a)') ' \begin{center}'
 write(LATEX_unit, '(a)') '\fcolorbox{black}{vert_titre}'
 write(LATEX_unit, '(a)') '{'
 write(LATEX_unit, '(a)') '\parbox{8cm}{'
 write(LATEX_unit, '(a)') '\begin{center}'
 write(LATEX_unit, '(a)') '\textcolor{fonce_titre}{\textit{#1}}'
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

 if(len_trim(cryscal_path_name) /= 0) then
  file_exist = .false.
  write(logo_full(1), '(3a)') trim(cryscal_path_name), '\img\', trim(logo(1))
  inquire(file = trim(logo_full(1)), exist = file_exist(1))
  write(logo_full(2), '(3a)') trim(cryscal_path_name), '\img\', trim(logo(2))
  inquire(file = trim(logo_full(2)), exist = file_exist(2))
  
  if(file_exist(1) .and. file_exist(2)) then
   ! les fichiers doivent etre dans le repertoire de travail   
   write(DOS_command, '(3a)') 'copy ', trim(logo_full(1)), ' .'     
   call system(trim(DOS_command))
   
   write(DOS_command, '(3a)') 'copy ', trim(logo_full(2)), ' .'
   call system(trim(DOS_command))   
   
  
   write(LATEX_unit, '(a)') ''
   write(LATEX_unit, '(a)') '\begin{figure}[h]'
   !write(LATEX_unit, '(a)') '% \includegraphics[width=50.pt]{img/ISCR_logo.png} \includegraphics[width=40.pt]{img/cdifx_logo.jpg}'
   write(LATEX_unit, '(5a)') ' \includegraphics[width=50.pt]{', trim(logo(1)), '}  \includegraphics[width=50.pt]{', &
                             trim(logo(2)), '}'
   write(LATEX_unit, '(a)') '\end{figure}'
   write(LATEX_unit, '(a)') ''
   
   
  endif
 endif

 return
end subroutine Latex_preambule

!----------------------------------------------------------------------------------------------------------------------------


subroutine Latex_End
 use LATEX_module,   only : logo
 use cryscal_module, only : LATEX_unit, cryscal_version, debug_proc
 
 implicit none
 
 if(debug_proc%level_2)  call write_debug_proc_level(2, "latex_end")

 write(LATEX_unit, '(a)')''
 write(LATEX_unit, '(4a)') '\footnote{\textcolor{gris50}{This structural report has been created through ', &
                           '\href{http://www.cdifx.univ-rennes1.fr/CRYSCAL}{CRYSCAL} (', trim(cryscal_version), ').'
 write(LATEX_unit, '(a)') 'Please report bugs and problems to \href{mailto:cdifx@univ-rennes1.fr}{cdifx@univ-rennes1.fr}}}'

 write(LATEX_unit, '(a)') '\end{document}'

 !call system('del '//trim(logo(1)))
 !call system('del '//trim(logo(2)))
   
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
 use cryscal_module,           only : HTML_unit, LATEX_unit, CIF_parameter, HTML_report, LATEX_report, LATEX_unit, debug_proc
 use special_characters_module
 implicit none
  character (len=*), intent(in)    :: diffracto    ! APEX, KCCD, X2S
  CHARACTER (LEN=256)              :: HTML_string, LATEX_string
  real                             :: T_min, T_max
  
  if(debug_proc%level_2)  call write_debug_proc_level(2, "report_crystal_study")

  if(len_trim(diffracto) == 3 .and. diffracto(1:3) == "X2S") then
   !if(len_trim(doc_type) == 4 .and. doc_type(1:4) == "HTML") then
   if(HTML_report) then
    WRITE(HTML_unit, '(a)')   "<p class='retrait_1'>"
    WRITE(HTML_unit, '(a)')   ""
	write(HTML_unit, '(5a)')  "A ", trim(CIF_parameter%crystal_colour), " crystal of ", trim(job), " with the dimensions of " 
	write(HTML_unit, '(10a)') TRIM(CIF_parameter%crystal_size_max), " mm x ", TRIM(CIF_parameter%crystal_size_mid), " mm x ",     &
	         TRIM(CIF_parameter%crystal_size_min), " mm mounted on a Mitegen Micromount was automatically centered on a Bruker ", &
			 "SMART X2S benchtop crystallographic system. Intensity measurements were performed using monochromated (doubly ",    &
			 "curved silicon crystal) Mo-K&alpha; radiation (&lambda;=0.71073 �) from a sealed microfocus tube. Generator ",      &
			 "settings were 50 kV, 1 mA. Data collection temperature was 21�C. Data were acquired using three sets of &Omega; ",  &
			 "scans at different &phi; settings."
	!		 The frame width was 0.5� with an exposure time of 5.0 s."
	write(HTML_unit, '(a)')  ""
    write(HTML_unit, '(a)')  "The detailed data collection stategy was as follows:<br>"
    write(HTML_unit, '(a)')  "&nbsp;&nbsp;  Detector distance : 40 mm<br>"
	write(HTML_unit, '(a)')  "&nbsp;&nbsp;  Detector swing angle (fixed 2&theta;) : -20�"	
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
							 "maximum &theta; angle of ", TRIM(CIF_parameter%theta_max), " �."
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
	!		 The frame width was 0.5\degre with an exposure time of 5.0 s."
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
							 "maximum $\theta$ angle of ", TRIM(CIF_parameter%theta_max), " �.\\"
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
   
	WRITE(LATEX_unit, '(2a,F4.1,6a)') "Data were corrected for absorption effects with SADABS [3] using the multiscan technique. ", &
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
 use cryscal_module, only  : HTML_unit, LATEX_unit, CIF_parameter, debug_proc
 implicit none
  character (len=*), intent(in)     :: doc_type
  
  if(debug_proc%level_3)  call write_debug_proc_level(3, "report_cell_parameters")

  
  if(len_trim(doc_type) == 4 .and. doc_type(1:4) == "HTML") then
   IF(CIF_parameter%symmetry_cell_setting(1:9) == 'triclinic') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " �, "
     WRITE(HTML_unit, '(3a)') "&alpha; = ",TRIM(CIF_parameter%cell_angle_alpha), ", "
     WRITE(HTML_unit, '(3a)') "&beta; = ",TRIM(CIF_parameter%cell_angle_beta),  ", "
     WRITE(HTML_unit, '(3a)') "&gamma; = ",TRIM(CIF_parameter%cell_angle_gamma), " �, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:10) == 'monoclinic') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " �, "
     WRITE(HTML_unit, '(3a)') "&beta; = ",TRIM(CIF_parameter%cell_angle_beta), " �, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:12) == 'orthorhombic') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " �, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:10) == 'tetragonal') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " �, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:9) == 'hexagonal' .or.    &
           CIF_parameter%symmetry_cell_setting(1:8) == 'trigonal') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
     WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " �, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:11) == 'rhomboedral') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), " �, "
     WRITE(HTML_unit, '(3a)') "&alpha; = ",TRIM(CIF_parameter%cell_angle_alpha), " �, "
   ELSEIF(CIF_parameter%symmetry_cell_setting(1:5) == 'cubic') then
     WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_c), " �, "
   ENDIF
   WRITE(HTML_unit, '(3a)', advance='NO') "<i>V</i> = ",TRIM(CIF_parameter%cell_volume), " �<sup>3</sup>."
   
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







