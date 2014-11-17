!     Last change:  TR   11 Dec 2007   12:21 pm
!Last change:  TR   22 Sep 2005   10:34 am

subroutine write_CRYSCALC_title()
 USE IO_module
 USE text_module,   ONLY : title_lines_nb, title_line
 implicit none
  integer :: i

 call def_title_lines()
 call def_CIF_file_title()
 do i= 1, title_lines_nb
  call write_info(TRIM(title_line(i)))
 end do


end subroutine write_cryscalc_title

!----------------------------------------------------
subroutine write_header 
 USE IO_module,      ONLY : write_info
 USE text_module,    ONLY : header_lines_nb, header_line
 implicit none
  INTEGER :: i

 call def_header_lines()
 do i= 1,  header_lines_nb
  call write_info(TRIM(header_line(i)))
 end do


 RETURN
end subroutine write_header

!--------------------------------------------------------------
subroutine write_cryscalc_CLA
 USE IO_module,      only : write_info
 USE text_module,    only : CLA_nb, CLA_lines_nb, CLA_line
 implicit none
  integer                 :: i,k
 
 call Def_command_line_arguments
 do i=1, CLA_nb
  do k=1, CLA_lines_nb(i)
  call write_info(trim(CLA_line(i,k)))
  end do
  call write_info('')
 end do
 
end subroutine write_cryscalc_CLA
 

!-------------------------------------------------------------------------------------------------------

subroutine HELP_on_line
 use cryscalc_module
 use macros_module
 USE IO_module,    only : write_info
 USE text_module
  implicit none
   integer                       :: i, k, help_numor
   LOGICAL                       :: keyword_ok


  IF(nb_help == 0) return


  IF(nb_help == nb_help_max) then
   call write_info('')
   call write_info('  >>> LIST OF CRYSCALC KEYWORDS:')
   call write_info('      -------------------------')
   call write_info('')
  endif

  keyword_ok = .false.

  call def_KEYWORDS_lines()

  do i=1, nb_help
   !select case(trim(help_string(i)))
   select case(trim(help_arg(i)))

    case ("ABSENT_HKL", "HKL_ABSENT")
     call write_help_lines(HELP_ABSENT_HKL_numor)
     keyword_ok = .true.
	 
	case ("ABSORPTION", "ABSORPTION_CALC", "CALC_ABSORPTION", "MU", "MU_CALC", "CALC_MU")	
     call write_help_lines(HELP_ABSORPTION_numor)
     keyword_ok = .true.

    case ("ACTA", "CIF")
     call write_help_lines(HELP_ACTA_numor)
     keyword_ok = .true.

    case ("ANG", "ANGLE")
     call write_help_lines(HELP_ANG_numor)
     keyword_ok = .true.

    case ("APPLY_OP", "APPLY_SYMMETRY_OPERATOR")
     call write_help_lines(HELP_APPLY_OP_numor)
     keyword_ok = .true.

    case ("ATOM", "ATM")
     call write_help_lines(HELP_ATOM_numor)
     keyword_ok = .true.

    case ("ATOM_LIST", "ATOM_LST", "LIST_ATOM_LIST", "LIST_ATOMS", "LST_ATOMS", "WRITE_ATOMS", "WRITE_ATMS")
     call write_help_lines(HELP_ATOM_LIST_numor)
     keyword_ok = .true.
	 
	case("WRITE_ADP", "WRITE_UIJ")
	 call write_help_lines(HELP_WRITE_ADP_numor)
	 keyword_ok = .true.

    case ("BARY", "CENTROID")
     call write_help_lines(HELP_BARY_numor)
     keyword_ok = .true.

    case ("BEAM")
     call write_help_lines(HELP_BEAM_numor)
     keyword_ok = .true.

    case ("CELL", "CELL_PARAMETERS")
     call write_help_lines(HELP_CELL_numor)
     keyword_ok = .true.

    CASE ("CHEM", "CHEM_FORM", "CHEMICAL_FORMULA")
     call write_help_lines(HELP_CHEM_numor)
     keyword_ok = .true.

    CASE ("CONN", "CONNECT")
     call write_help_lines(HELP_CONN_numor)
     keyword_ok = .true.
     
    CASE ("CONT")
     call write_help_lines(HELP_CONT_numor)
     keyword_ok = .true.

    CASE ("CREATE_ACE")
     call write_help_lines(HELP_CREATE_ACE_numor)
     keyword_ok = .true.

    CASE ("CREATE_CEL")
     call write_help_lines(HELP_CREATE_CEL_numor)
     keyword_ok = .true.

    CASE ("CREATE_CFL")
     call write_help_lines(HELP_CREATE_CFL_numor)
     keyword_ok = .true.

    CASE ("CREATE_FST")
     call write_help_lines(HELP_CREATE_FST_numor)
     keyword_ok = .true.

    CASE ("CREATE_INS")
     call write_help_lines(HELP_CREATE_INS_numor)
     keyword_ok = .true.

    case ("REPORT", "CREATE_REPORT")
     call write_help_lines(HELP_CREATE_REPORT_numor)
     keyword_ok = .true.
	 
    case ('CREATE_SOLVE', 'CREATE_TO_SOLVE', 'CREATE_FILES_TO_SOLVE', 'SOLVE')
      call write_help_lines(HELP_CREATE_SOLVE_numor)
	  keyword_ok = .true.
	  
    case ("CREATE_TIDY", "CREATE_TIDY_FILE", "CREATE_TIDY_INPUT_FILE") 
	 call write_help_lines(HELP_CREATE_TIDY_numor)
     keyword_ok = .true.

    case ("D_HKL", "DHKL")
     call write_help_lines(HELP_D_HKL_numor)
     keyword_ok = .true.

    case ("D_STAR", "D_STAR_HKL", "DSTAR", "DSTARHKL", "DSTAR_HKL")
     call write_help_lines(HELP_D_STAR_numor)
     keyword_ok = .true.

    case ("DATA_DENSITY", "DENSITY_DATA", "DATA_ATOMIC_DENSITY", "ATOMIC_DENSITY")
     call write_help_lines(HELP_DATA_ATOMIC_DENSITY_numor)
     keyword_ok = .true.

    case ("DATA_RADIUS", "RADIUS_DATA", "DATA_ATOMIC_RADIUS", "ATOMIC_RADIUS_DATA")
     call write_help_lines(HELP_DATA_ATOMIC_RADIUS_numor)
     keyword_ok = .true.

    case ("DATA_WEIGHT", "WEIGHT_DATA", "DATA_ATOMIC_WEIGHT", "ATOMIC_WEIGHT")
     call write_help_lines(HELP_DATA_ATOMIC_WEIGHT_numor)
     keyword_ok = .true.

   case ("DATA_NEUTRONS", "NEUTRONS_DATA", "DATA_NEUTRON", "NEUTRON_DATA")
     call write_help_lines(HELP_DATA_NEUTRONS_numor)
     keyword_ok = .true.

    case ("DATA_XRAYS", "XRAYS_DATA", "DATA_XRAY", "XRAY_DATA")
     call write_help_lines(HELP_DATA_XRAYS_numor)
     keyword_ok = .true.


    CASE ('DIAG', 'DIAG_MAT', 'DIAG_MATR', 'DIAG_MATRIX')
     call write_help_lines(HELP_DIAG_MAT_numor)
     keyword_ok = .true.
	 
    case ("DIR_ANG", "DIRANG", "DIRECT_ANGLE")
     call write_help_lines(HELP_DIR_ANG_numor)
     keyword_ok = .true.


    case ("DIST", "DISTANCE", "ATOMIC_DISTANCE")
     call write_help_lines(HELP_DIST_numor)
     keyword_ok = .true.
	 
	CASE ('DIST_DHA', 'DHA', 'POS_H', 'CALC_POS_H')	
	 call write_help_lines(HELP_DIST_DHA_numor)
     keyword_ok = .true.

    case ("EDIT")
     call write_help_lines(HELP_EDIT_numor)
     keyword_ok = .true.

    case ("EQUIV", "EQUIV_HKL", "SEARCH_EQUIV", "SEARCH_EQUIV_HKL", "FIND_EQUIV", "FIND_EQUIV_HKL")
     call write_help_lines(HELP_EQUIV_numor)
     keyword_ok = .true.

    case ("EXIT", "X", "QUIT", "END", "STOP")
     call write_help_lines(HELP_EXIT_numor)
     keyword_ok = .true.

    case ("FILE")
     call write_help_lines(HELP_FILE_numor)
     keyword_ok = .true.

    case ("FIND_HKL", "FINDHKL", "SEARCH_HKL")
     call write_help_lines(HELP_FIND_HKL_numor)
     keyword_ok = .true.

	case ("FRIEDEL", "FRIEDEL_PAIRS")
     call write_help_lines(HELP_Friedel_pairs_numor)
     keyword_ok = .true.
    
    case ("GEN_HKL", "GENERATE_HKL", "GENERATE_HKL_LIST")
     call write_help_lines(HELP_GEN_HKL_numor)
     keyword_ok = .true.

    case ("HEADER", "HEAD")
     call write_help_lines(HELP_HEADER_numor)
     keyword_ok = .true.

    case ("HEX_RHOMB", "HEXA_RHOMB", "HEX_TO_RHOMB", "HEXA_TO_RHOMB")
     call write_help_lines(HELP_HEX_RHOMB_numor)
     keyword_ok = .true.

    case ("HKL")
     call write_help_lines(HELP_HKL_numor)
     keyword_ok = .true.

    case ("HKL_NEG", "HKL_NEGATIVE", "NEG_HKL", "NEGATIVE_HKL")
     call write_help_lines(HELP_HKL_NEG_numor)
     keyword_ok = .true.

    case ("HKL_POS", "HKL_POSITIVE", "POS_HKL", "POSITIVE_HKL")
     call write_help_lines(HELP_HKL_POS_numor)
     keyword_ok = .true.

    case ("INSIDE")
     call write_help_lines(HELP_INSIDE_numor)
     keyword_ok = .true.

    case ("LIST_EXTI", "LIST_EXTI_RULE", "LST_EXTI", "LST_EXTI_RULE")
     call write_help_lines(HELP_LIST_EXTI_numor)
     keyword_ok = .true.

 
    case ("KEY", "KEYS", "LST_KEYS", "LIST_KEYS", "LST_KEYWORDS", "LIST_KEYWORDS")
     call write_help_lines(HELP_LIST_KEYS_numor)
     keyword_ok = .true.

    case ("LST_LAUE", "LIST_LAUE", "LST_LAUE_CLASS", "LIST_LAUE_CLASS")
     call write_help_lines(HELP_LIST_LAUE_numor)
     keyword_ok = .true.

    CASE ("LST_MAT", "LST_MATR", "LST_MATRIX",  "LIST_MAT",  "LIST_MATR",  "LIST_MATRIX", "LIST_TRANSFORMATION_MATRIX")
     call write_help_lines(HELP_LIST_MATR_numor)
     keyword_ok = .true.

    CASE ("LST_SG", "LIST_SG", "LIST_SPACE_GROUPS")
     call write_help_lines(HELP_LIST_SG_numor)
     keyword_ok = .true.

    case ("MAG", "MAGNETIC")
     call write_help_lines(HELP_MAG_numor)
     keyword_ok = .true.

    case ("MAN", "HELP")
     call write_help_lines(HELP_MAN_numor)
     keyword_ok = .true.
	 	
    case ("MAN_HTML", "HTML_MAN", "HTML")
     call write_help_lines(HELP_MAN_HTML_numor)
     keyword_ok = .true.

    case ("MAT", "MATR", "MATRIX")
     call write_help_lines(HELP_MATR_numor)
     keyword_ok = .true.

    case ("MATMUL")
     call write_help_lines(HELP_MATMUL_numor)
     keyword_ok = .true.
     
    case ("MENDEL")
     call write_help_lines(HELP_MENDEL_numor)
     keyword_ok = .true.

    case ("TRICLINIC", "TRICL")
     call write_help_lines(HELP_TRICLINIC_numor)
     keyword_ok = .true.

    case ("MONOCLINIC", "MONOC", "MONOCL")
     call write_help_lines(HELP_MONOCLINIC_numor)
     keyword_ok = .true.

    case ("NEWS")
     call write_help_lines(HELP_NEWS_numor)
     keyword_ok = .true.

    case ("NIGGLI", "NIGGLI_CELL")
     call write_help_lines(HELP_NIGGLI_cell_numor)
     keyword_ok = .true.


    case ("OBV_REV", "OBVERSE_REVERSE")
     call write_help_lines(HELP_OBV_REV_numor)
     keyword_ok = .true.

    case ("P4P")
     call write_help_lines(HELP_P4P_numor)
     keyword_ok = .true.

    case ("PAUSE")
     call write_help_lines(HELP_PAUSE_numor)
     keyword_OK = .true. 

    case ("PERMUT")
     call write_help_lines(HELP_PERMUT_numor)
     keyword_ok = .true.

    case ("QHKL")
     call write_help_lines(HELP_Q_HKL_numor)
     keyword_ok = .true.

    case ("QVEC", "Q_VEC", "Q_VECTOR", "KVEC", "K_VEC", "K_VECTOR", "MOD_VEC", "MODULATION_VECTOR")
     call write_help_lines(HELP_QVEC_numor)
     keyword_ok = .true.

    case ("READ_CEL", "READ_CEL_FILE", "READ_POWDERCELL")
     call write_help_lines(HELP_READ_CEL_numor)
     keyword_ok = .true.

    case ("READ_FACES", "FACES")
	 call write_help_lines(HELP_READ_FACES_numor)
	 keyword_ok = .true.
	 
	case ("READ_INS", "READ_INS_FILE", "INS_FILE")
     call write_help_lines(HELP_READ_INS_numor)
     keyword_ok = .true.

    case ("READ_CIF", "READ_CIF_FILE", "CIF_FILE")
     call write_help_lines(HELP_READ_CIF_numor)
     keyword_ok = .true.
    
    case ("READ_NREPORT", "READ_NREPORT_HTML", "READ_HTMLREPORT")
     call write_help_lines(HELP_READ_NREPORT_numor)
     keyword_ok = .true.

	case ("READ_PCR", "READ_PCR_FILE", "PCR_FILE")
     call write_help_lines(HELP_READ_PCR_numor)
     keyword_ok = .true.

    case ("READ_TIDY_OUT")
     call write_help_lines(HELP_READ_TIDY_out_numor)
     keyword_ok = .true.

    case ("REC_ANG", "RECANG", "RECIPROCAL_ANGLE")
     call write_help_lines(HELP_REC_ANG_numor)
     keyword_ok = .true.

    case ("REF_APEX")
     call write_help_lines(HELP_REF_APEX_numor)
     keyword_ok = .true.
	 
	case ("REF_X2S")
      call write_help_lines(HELP_REF_X2S_numor)
      keyword_ok = .true.	  

    case ("REF_DENZO")
     call write_help_lines(HELP_REF_DENZO_numor)
     keyword_ok = .true.

    case ("REF_EVAL")
     call write_help_lines(HELP_REF_EVAL_numor)
     keyword_ok = .true.

    case ("REF_KCCD")
     call write_help_lines(HELP_REF_KCCD_numor)
     keyword_ok = .true.
     
    case ("REF_SADABS")
     call write_help_lines(HELP_REF_SADABS_numor)
     keyword_ok = .true.  

    case ("RESET", "RAZ", "INIT", "INITIALIZATION")
     call write_help_lines(HELP_RESET_numor)
     keyword_ok = .true.

    case ("RINT", "R_INT")
     call write_help_lines(HELP_RINT_numor)
     keyword_ok = .true.

    case ("MERGE")
     call write_help_lines(HELP_MERGE_numor)
     keyword_ok = .true.

    case ("RHOMB_HEX", "RHOMB_HEXA", "RHOMB_TO_HEX", "RHOMB_TO_HEXA")
     call write_help_lines(HELP_RHOMB_HEX_numor)
     keyword_ok = .true.

    case ("SEARCH_EXTI", "FIND_EXTI")
     call write_help_lines(HELP_SEARCH_EXTI_numor)
     keyword_ok = .true.

    case ("SEARCH_SPGR", "SEARCH_SPACE_GROUP", "SEARCH_GROUP",   &
          "CHECK_SPGR",  "CHECK_SPACE_GROUP",  "CHECK_GROUP")
     call write_help_lines(HELP_SEARCH_SPGR_numor)
     keyword_ok = .true.     

    case ("SET")
     call write_help_lines(HELP_SET_numor)
     keyword_ok = .true.

    case ("SETTING")
     call write_help_lines(HELP_SETTING_numor)
     keyword_ok = .true.

    CASE("SFAC")
     call write_help_lines(HELP_SFAC_numor)
     keyword_ok = .true.

    case ("SF_HKL", "SFAC_HKL")
     call write_help_lines(HELP_SFHKL_numor)
     keyword_ok = .true.

    case ("SG", "SPGR", "SPACE_GROUP")
     call write_help_lines(HELP_SG_numor)
     keyword_ok = .true.

    case ("SG_ALL",  "SP_ALL",  "SG_INFO_EXTI", "SP_INFO_EXTI", "SG_EXTI_INFO", "SP_EXTI_INFO")
     call write_help_lines(HELP_SG_ALL_numor)
     keyword_ok = .true.


    case ("SG_EXTI", "SP_EXTI", "SG_EXTINCTIONS", "SPACE_GROUP_EXTI", "SPACE_GROUP_EXTINCTIONS", &
	      "LIST_EXTINCTIONS, LIST_SPACE_GROUP_EXINCTIONS")
     call write_help_lines(HELP_SG_EXTI_numor)
     keyword_ok = .true.

    case ("SG_INFO", "SP_INFO", "SPACE_GROUP_INFO", "LIST_SPACE_GROUP_INFO")
     call write_help_lines(HELP_SG_INFO_numor)
     keyword_ok = .true.

    case ("SHAN", "SHANNON")
     call write_help_lines(HELP_SHANNON_numor)
     keyword_ok = .true.

    CASE("SHELL")
     call write_help_lines(HELP_SHELL_numor)
     keyword_ok = .true.

    case ("SHIFT_2TH", "SHIFT_2THETA", "2TH_SHIFT", "2THETA_SHIFT")
     call write_help_lines(HELP_SHIFT_2TH_numor)
     keyword_ok = .true.

    case ("SITE_INFO", "LIST_SITE_INFO")
     call write_help_lines(HELP_SITE_INFO_numor)
     keyword_ok = .true.

    CASE("SIZE", "CRYSTAL_SIZE")
     call write_help_lines(HELP_SIZE_numor)
     keyword_ok = .true.

    CASE("SORT")
     call write_help_lines(HELP_SORT_numor)
     keyword_ok = .true.

    case ("STL", "SINTHETA/WAVE", "SINTHETA/LAMBDA")
     call write_help_lines(HELP_STL_numor)
     keyword_ok = .true.

	case ('STAR_K')
     call write_help_lines(HELP_STAR_K_numor)
     keyword_ok = .true.
	 
    case ("SYMM", "SYM", "SYMMETRY_OPERATOR")
     call write_help_lines(HELP_SYMM_numor)
     keyword_ok = .true.

    case ("SYST", "CMD", "COMMAND", "DOS", "DOS_COMMAND")
     call write_help_lines(HELP_SYST_numor)
     keyword_ok = .true.

    case ("THERM", "THERMAL", "ADP")
     call write_help_lines(HELP_THERM_numor)
     keyword_ok = .true.
	
	case ("THERM_SHELX", "THERMAL_SHELX", "ADP_SHELX")
     call write_help_lines(HELP_THERM_SHELX_numor)
     keyword_ok = .true. 

    case ('THETA', 'TH', 'TH_HKL', 'THETAHKL', 'THETA_HKL')
     call write_help_lines(HELP_THETA_numor)
     keyword_ok = .true.

    case ("TITL", "TITLE")
     call write_help_lines(HELP_TITL_numor)
     keyword_ok = .true.

    case ("TRANSLATION", "TRANS", "TRANSLATE", "MOVE" )
     call write_help_lines(HELP_TRANSLATION_numor)
     keyword_ok = .true.

    case ("TRANSMISSION")
     call write_help_lines(HELP_TRANSMISSION_numor)
     keyword_ok = .true.

    case ('2THETA', '2TH', '2TH_HKL', '2THETA_HKL', 'TWO_THETA', 'TWO_THETA_HKL')
     call write_help_lines(HELP_TWO_THETA_numor)
     keyword_ok = .true.

    case ('HEXA_TWIN', 'HEXA_TWINNING',  'HEXAGONAL_TWIN', 'HEXAGONAL_TWINNING', &
          'TWIN_HEXA', 'TWIN_HEXAGONAL', 'TWINNING_HEXA',  'TWINNING_HEXAGONAL')
     call write_help_lines(HELP_TWIN_HEXA_numor)
     keyword_ok = .true.
 
    case ('TWIN_PSEUDO_HEXA')
     call write_help_lines(HELP_TWIN_PSEUDO_HEXA_numor)
     keyword_ok = .true.
     
	CASE ('UB_MAT', 'UB_MATRIX')
	 call  write_help_lines(HELP_UB_matrix_numor)
     keyword_ok = .true.
	 
    CASE("UNIT")
     call write_help_lines(HELP_UNIT_numor)
     keyword_ok = .true.

    case ("WAVE", "WAVELENGTH", "WL")
     call write_help_lines(HELP_WAVE_numor)
     keyword_ok = .true.

    case ("WEB", "INTERNET")
     call write_help_lines(HELP_WEB_numor)
     keyword_ok = .true.


    case ("WRITE_CELL", "OUTPUT_CELL")
     call write_help_lines(HELP_WRITE_CELL_numor)
     keyword_ok = .true.

    case ("WRITE_CHEM", "OUTPUT_CHEM")
     call write_help_lines(HELP_WRITE_CHEM_numor)
     keyword_ok = .true.

    case ("WRITE_QVEC", "OUTPUT_QVEC")
     call write_help_lines(HELP_WRITE_QVEC_numor)
     keyword_ok = .true.



   !CASE ("LST_HKL", "LIST_HKL", "LST_HKL", "LIST_HKL_EXTI", "LST_HKL_EXTI")
   !CASE ("WRITE_HKL", "WRITE_HKL_LIST", "WRITE_HKL_LST")
    CASE ("FIND_HKL_LIST", "FIND_HKL_LST", "EXTRACT_HKL_LIST", "EXTRACT_HKL_LST") 
     call write_help_lines(HELP_FIND_HKL_LIST_numor)
     keyword_ok = .true.

    case ("WRITE_SG", 'WRITE_SPACE_GROUP')
     call write_help_lines(HELP_WRITE_SG_numor)
     keyword_ok = .true.

    CASE ("WRITE_SYM_OP", "WRITE_SYMM_OP", "WRITE_SYM", "WRITE_SYMM", "WRITE_SYMMETRY_OPERATORS", "WRITE_SYMOP")
     call write_help_lines(HELP_WRITE_SYM_OP_numor)
     keyword_ok = .true.

    case ("WRITE_WAVE", "OUTPUT_WAVE")
     call write_help_lines(HELP_WRITE_WAVE_numor)
     keyword_ok = .true.
     
    case ("WRITE_DEVICE", "OUTPUT_DEVICE")
     call write_help_lines(HELP_WRITE_DEVICE_numor)
     keyword_ok = .true. 

    case ("WRITE_BEAM", "WRITE_INCIDENT_BEAM", "OUTPUT_BEAM", "OUTPUT_INCIDENT_BEAM")
     call write_help_lines(HELP_WRITE_BEAM_numor)
     keyword_ok = .true.

    case ("XRAYS_WAVELENGTH", "X_WAVE")
     call write_help_lines(HELP_X_WAVE_numor)
     keyword_ok = .true.

    case ("ZUNIT", "Z", "Z_UNIT")
     call write_help_lines(HELP_ZUNIT_numor)
     keyword_ok = .true.

    case default
     call write_info('')
     call write_info(' ... unknown '//trim(help_arg(i))//' keyword ...')
     call write_info('')
     
   end select
  end do

  !IF(.not.  keyword_ok) then
  !  call write_info('')
  !  call write_info(' ... unknown MAN argument ...')
  !  call write_info('')
  !endif




 return
end subroutine HELP_on_line

!------------------------------------------------------------------------------------
subroutine write_help_lines(help_numor)
 USE text_module,   ONLY : HELP_line, HELP_lines_nb
 USE IO_module,     ONLY : write_info
 INTEGER, INTENT(IN) :: help_numor
 INTEGER             :: k

 do k=1, HELP_lines_nb(help_numor)
  call write_info(TRIM(HELP_line(help_numor,k)))
 end do

end subroutine write_help_lines
!------------------------------------------------------------------------------------

subroutine expert_help
 use IO_module,   ONLY : write_info
 
 call write_info('')
 call write_info(' > List of specific keywords (for experts only) :')
 call write_info('')
 call write_info('  . DEBUG         :               access to debug mode (cryscalc_debug.txt)')
 call write_info('  . DEBUG_2       :               access to debug mode (level #2)')
 call write_info('  . DEBUG_3       :               access to debug mode (level #3)')
 call write_info('  . FIC           =               FILE import.cif')
 call write_info('  . FIT           =               FILE import_trans.hkl')
 call write_info('  . FSD           =               FILE import_shell_d.hkl')
 call write_info('  . FST           =               FILE import_shell_theta.hkl')
 call write_info('  . PDP           =               GEN_HKL theta_min=0. theta_max=120 pat')
 call write_info('  . PDP_Cu        =               GEN_HKL theta_min=0. theta_max=120 pat (Ka Cu radiation)')
 call write_info('  . SD7           =               SHELL d 0.77 7')
 call write_info('  . ST25          =               SHELL theta 2.5 25.')
 call write_info('  . S_G, GSG, SSG =               SEARCH_GROUP get')
 call write_info('  . WC            =               WRITE_CELL')
 call write_info('')
 call write_info(' > List of expert options in cryscalc.ini') 
 call write_info('')
 call write_info('  . bond_str_out  =               1/0') 
 call write_info('  . debug_level_2 =               1/0')
 call write_info('  . debug_level_3 =               1/0')
 
 call write_info('')


 return
end subroutine expert_help