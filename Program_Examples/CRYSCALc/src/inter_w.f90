!     Last change:  TR   17 Jul 2007    4:07 pm
subroutine interactive_mode(input_string)
 USE macros_module, ONLY : u_case 
 USE IO_module
 USE cryscalc_module
 USE HKL_module,    ONLY : HKL_list, search_EQUIV, search_FRIEDEL, requested_H
 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: input_string
  CHARACTER (LEN=256)           :: read_line
  INTEGER                       :: i1, i_error
  CHARACTER (LEN=32)            :: current_keyword

  if(debug_proc%level_2)  call write_debug_proc_level(2, "interactive_mode ("//trim(input_string)//")")
  
 IF(input_string(1:8) == 'keyboard') then
  close(unit=CFL_read_unit)
  open (UNIT=CFL_read_unit, file = 'cryscalc.cfl')
 endif
 !open (UNIT=11, FILE="cryscalc.CFL")

 do
  current_keyword = ''

  IF(input_string(1:8) == 'keyboard') then
   call write_info(' ')
   call write_info(' > Enter keyword: ')
   call read_input_line(input_line)
   READ(input_line,'(a)', IOSTAT=i_error) read_line
   IF(i_error /=0) cycle

  ELSEIF(input_string(1:4) == 'file') then
   READ(UNIT=CFL_read_unit, fmt='(a)', IOSTAT=i_error) read_line
   !READ(UNIT=input_unit, '(a)', IOSTAT=i_error) read_line
   IF(i_error < 0) EXIT   ! fin du fichier
   ! stocke la commande dans le fichier cryscalc.cfl
   ! write(unit = CFL_unit, fmt='(a)') trim(read_line)
  endif

  IF(LEN_TRIM(read_line) == 0) cycle
  IF(read_line(1:1) == '! ' .or. read_line(1:1) == '#' ) CYCLE  ! commentaire
  i1 = index(read_line, '!')
  if(i1 /=0) read_line = trim(read_line(1:i1-1))  
 

  read_line = ADJUSTL(read_line)
  read_line = u_case(read_line)

  !IF(read_line(1:3) == 'END'  .OR. read_line(1:4) == 'QUIT') call end_of_program
  !IF(read_line(1:4) == 'EXIT' .OR. read_line(1:4) == 'STOP') call end_of_program
  !IF(LEN_TRIM(read_line) == 1) then
  ! IF(read_line(1:1) =='X' .OR. read_line(1:1) == 'Q' .or. read_line(1:1) == '0') call end_of_program
  !endif
  IF(read_line(1:3) == 'END'  .OR. read_line(1:4) == 'QUIT') exit
  IF(read_line(1:4) == 'EXIT' .OR. read_line(1:4) == 'STOP') exit
  IF(LEN_TRIM(read_line) == 1) then
   IF(read_line(1:1) =='X'  .OR. read_line(1:1) == 'Q'  .or. read_line(1:1) == '0')  exit
  endif
  IF(read_line(1:2) =='XX' .OR. read_line(1:2) == 'QQ' .or. read_line(1:2) == '00') call end_of_program

  READ(read_line, *) current_keyword
  current_keyword = ADJUSTL(current_keyword)
  current_keyword = u_case(current_keyword)

 ! open (UNIT=CFL_unit, file = 'cryscalc.cfl')
  unknown_keyword     = .false.
  unknown_CFL_keyword = .false.

  call write_info('')
  call write_info('    >> Input command: '//trim(read_line))
  call write_info('')

  call identification_CFL_keywords(read_line)       ! read_CFL_.F90
  call identification_keywords(read_line)           ! read_KEYW.F90

  IF(unknown_keyword .and. unknown_CFL_keyword) then
   call write_info('')
   call write_info('  > Unknown '//TRIM(current_keyword)//' keyword !')
   call write_info('')
   cycle
  endif

  call run_keyword_interactive(current_keyword)

 END do

 !CLOSE(UNIT=CFL_unit)
  close(unit=input_unit)

END subroutine interactive_mode

!------------------------------------------------------------------------------------------

 subroutine run_keyword_interactive(current_keyword)
 USE cryscalc_module
 USE HKL_module
 USE wavelength_module
 USE IO_module
 use external_applications_module, ONLY : launch_browser
 USE realwin, ONLY : SPAWN

 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: current_keyword
  LOGICAL                       :: create_file

if(debug_proc%level_2)  call write_debug_proc_level(2, "run_keyword_interactive ("//trim(current_keyword)//")")
  
  select case (TRIM(current_keyword))
      CASE ('EXPERT_MODE', 'EXPERT')
	   if(expert_mode) then
	    call write_info('')
	    call write_info(' > Expert mode activated.')
	    call write_info('')
	   else
	    call write_info('')
	    call write_info(' > Expert mode ON.')
	    call write_info('')	   
	   end if
    
      CASE ('USER_MODE', 'USER')
	   if(expert_mode) then
	    call write_info('')
	    call write_info(' > Expert mode OFF.')
	    call write_info('')
	   else
	    call write_info('')
	    call write_info(' > User mode ON.')
	    call write_info('')	   
	   end if
	   
     case ("DEBUG", "DEBUG_ON", "DEBUG ON")
	   if(expert_mode) then
	    call write_info('')
	    call write_info(' > Debug mode activated.')
	    call write_info('')
		call write_debug_level_1
	   end if	
      case ("&DEBUG", "&DEBUG_ON", "&DEBUG ON")	   
	    call write_info('')
	    call write_info(' > Debug mode activated.')
	    call write_info('')
		call write_debug_level_1	   
	
      case ("DEBUG_2", "DEBUG_2_ON", "DEBUG_2 ON", "DEBUG 2 ON")
	   if(expert_mode) then
	    call write_info('')
	    call write_info(' > Debug mode (level 2) activated.')
	    call write_info('')
	   end if	
      case ("&DEBUG_2", "&DEBUG_2_ON", "&DEBUG_2 ON", "&DEBUG 2 ON")	   
	    call write_info('')
	    call write_info(' > Debug mode (level 2) activated.')
	    call write_info('')	   

      case ("DEBUG_3", "DEBUG_3_ON", "DEBUG_3 ON", "DEBUG 3 ON")
	   if(expert_mode) then
	    call write_info('')
	    call write_info(' > Debug mode (level 3) activated.')
	    call write_info('')
	   end if	

      case ("&DEBUG_3", "&DEBUG_3_ON", "&DEBUG_3 ON", "&DEBUG 3 ON")	   
	    call write_info('')
	    call write_info(' > Debug mode (level 3) activated.')
	    call write_info('')	   
		
      case ("DEBUG_OFF", "DEBUG OFF")
	   if(expert_mode) then
   	    call write_info('')
	    call write_info(' > Debug mode OFF.')
	    call write_info('')
	   end if	
	  
	  case ("&DEBUG_OFF", "&DEBUG OFF")	   
   	    call write_info('')
	    call write_info(' > Debug mode OFF.')
	    call write_info('')
	   	


      CASE ('APPLY_OP', 'APPLY_SYMMETRY_OPERATOR', 'APPLY_SYM_OP', 'APPLY_SYMOP')        
       IF(nb_atom /=0)  call apply_symm

   	  CASE ('STAR_K')
        IF(keyword_QVEC .and. keyword_SPGR) call star_K()	  

      case ('CELL', 'CELL_PARAMETERS')
        IF(keyword_CELL) call volume_calculation('out')

      case ('CREATE_ACE')
        call create_ACE_from_CIF
        
      case ('CREATE_CEL')
        call create_CEL_from_CIF
        
      case ('CREATE_CFL')
        if(LEN_TRIM(CIF_file_name) /= 0) then
		 call create_CFL_file(TRIM(CIF_file_name), 'CIF')
        elseif(LEN_TRIM(INS_file_name) /= 0) then
		 call create_CFL_file(TRIM(INS_file_name), 'INS')
        elseif(LEN_TRIM(PCR_file_name) /= 0) then
		 call create_CFL_file(TRIM(PCR_file_name), 'PCR')
		else
     	 !call write_info('')
		 !call write_info('  !! A .CIF, .INS or .PCR file has to be read !!')
		 !call write_info('')
		 call create_CFL_file('', '')
        endif 		
        
      case ('CREATE_FST')
        if(LEN_TRIM(CIF_file_name) /= 0) then
		 call create_FST_file(TRIM(CIF_file_name), 'CIF')
        elseif(LEN_TRIM(INS_file_name) /= 0) then
		 call create_CFL_file(TRIM(INS_file_name), 'INS')
        elseif(LEN_TRIM(PCR_file_name) /= 0) then
		 call create_FST_file(TRIM(PCR_file_name), 'PCR')       
		else
      	 !call write_info('')
		 !call write_info('  !! A .CIF, .INS or .PCR file has to be read !!')
		 !call write_info('')
		 call create_FST_file('', '')       
        endif		
        
      case ('CREATE_INS')     
        create_file = .false.
        if(LEN_TRIM(CIF_file_name) /= 0) then
		 call create_INS_file(TRIM(CIF_file_name), 'CIF')
		 create_file = .true.
		endif
        if(LEN_TRIM(PCR_file_name) /= 0) then
		 call create_INS_file(TRIM(PCR_file_name), 'PCR')
		 create_file = .true.
		endif
        if(.not. create_file) call create_INS_file('cryscalc.ins', '')

      case ('CREATE_TIDY', 'CREATE_TIDY_FILE', 'CREATE_TIDY_INPUT_FILE')
	    if(LEN_TRIM(CIF_FILE_name) /=0)  then
		 call create_TIDY_file(TRIM(CIF_file_name), 'CIF')
		elseif(LEN_TRIM(INS_FILE_name) /=0)  then
		 call create_TIDY_file(TRIM(INS_file_name), 'INS')
        elseif(LEN_TRIM(PCR_FILE_name) /=0)  then
		 call create_TIDY_file(TRIM(PCR_file_name), 'PCR')
		else
		 !call write_info('')
		 !call write_info('  !! A .CIF, .INS or .PCR file has to be read !!')
		 !call write_info('')
		 call create_TIDY_file('', '')
        endif		

      case ('NIGGLI', 'NIGGLI_CELL')
        IF(keyword_NIGGLI) call Niggli_cell_CFML()

      case ('EDIT')
        IF(keyword_EDIT) call edit_a_file

      case ('INSIDE')
        if(nb_atom /=0)  call inside_unit_cell()

      case ('HEADER', 'HEAD')
       call write_header 

      case ('HELP', 'MAN')
        call help_on_line()
		
      case ('HELP_EXPERT', 'MAN_EXPERT')
        if(expert_mode) call expert_help	  

      case ('&HELP_EXPERT', '&MAN_EXPERT')
        call expert_help	
		

      case ('KEY', 'KEYS', 'LST_KEYS', 'LIST_KEYS', 'LST_KEYWORDS', 'LIST_KEYWORDS')
        call keys_on_line()

      CASE ('LST_ATOMS', 'ATOM_LST', 'ATOM_LIST', 'LIST_ATOM_LIST', 'LIST_ATOMS', 'WRITE_ATOMS', 'WRITE_ATMS', 'WRITE_ADP', 'WRITE_UIJ')
        IF(nb_atom /=0)  then
		if(.not. write_atoms_cartesian) then
		 call write_atom_list   !
		else
		 call write_atom_list_cart
        endif		
	    else
		 call write_info('')
		 call write_info(' > No input atoms !')
		 call write_info('')
        endif		

      case ('LST_LAUE', 'LIST_LAUE', 'LST_LAUE_CLASS', 'LIST_LAUE_CLASS')
        call list_Laue_class()

      CASE ('LST_MAT', 'LST_MATR', 'LST_MATRIX',  'LIST_MAT',  'LIST_MATR',  'LIST_MATRIX', 'LIST_TRANSFORMATION_MATRIX')
       call write_list_matrice()

      case ('LST_SG', 'LIST_SG', 'LSPGR', 'LIST_SPACE_GROUPS' )
        call list_space_groups()

      CASE ("WRITE_SYM_OP", "WRITE_SYMM_OP", "WRITE_SYM", "WRITE_SYMM", "WRITE_SYMMETRY_OPERATORS", "WRITE_SYMOP")
       IF(WRITE_SHELX_symm_op) then
		 call write_SHELX_symm_card
        ELSEIF(WRITE_symm_op) then
		 call write_symm_op_mat()
		ENDIF 


      case ('MAN_HTML', 'HTML_MAN', 'HTML')
        call create_CRYSCALC_HTML()

      case ('MATMUL')
        call multi_matrix()

      CASE ('MAT', 'MATR', 'MATRIX', 'USER_MAT', 'USER_MATR', 'USER_MATRIX')
        IF (keyword_MAT) then
         call WRITE_matrice()
         if (keyword_CELL)    call transf_cell_parameters
         if (nb_hkl  /=0)     call transf_HKL
         !if (nb_atom ==0) then
         ! call get_matrice_inverse_transposee
         !else
         ! call transf_coord()
         !endif
         call transf_coord()
         if (keyword_FILE)  then
          call transf_HKL_file()
          known_theta = .false.
          ordered_hkl = .false.
          HKL_data_known = .false.
          HKL_file%name = HKL_file%transf
          HKL_file%SHELX = .true.
          call read_and_sort_HKL("sort")
         endif

         IF(ABS(ABS(Mat_det) -1.) < 0.001) then
          !IF (input_INS) then
          ! call create_NEW_INS()
          !else
          ! IF (keyword_SPGR )   call find_new_group()
          !endif
          IF (keyword_SPGR )  call find_new_group()
          !IF (keyword_read_INS)      call create_NEW_INS()
          if(keyword_read_INS)        call create_TRANSF_INS()

         endif
        endif

      CASE ('UB_MAT', 'UB_MATRIX')
       if(keyword_UB_MAT) call write_UB_matrix

      CASE ('DIAG', 'DIAG_MAT', 'DIAG_MATR', 'DIAG_MATRIX')
        IF (keyword_DIAG) then
         call WRITE_matrice()
		 call DIAG_matrice()
        endif

		
      case ('NEWS')
        call write_cryscalc_NEWS('screen')

      case ('P4P', 'READ_P4P')
        if (keyword_P4P) then
         !call read_P4P('')
         call read_P4P(trim(P4P_file_name))
         call create_CIF_P4P_import('P4P')
         !call read_SHELX_HKL
         !stop
        endif

      case ('REPORT', 'CREATE_REPORT')
        call create_structural_report   ! create_HTML.F90

      case ('SETTING', 'SETTINGS')
        call output_setting

      case ('SFAC', 'CONT', 'CHEM', 'CHEM_FORM', 'CHEMICAL_FORMULA')
       IF(keyword_SFAC_UNIT .or. keyword_CONT .or. keyword_CHEM)  then
        call atomic_identification()
        IF(keyword_CELL)                      call atomic_density_calculation()
        IF(keyword_ZUNIT)	                  call molecular_weight()
        IF(keyword_CELL .and. keyword_ZUNIT)  call density_calculation()
        IF(keyword_CELL)                      call absorption_calculation()
       endif

     case ('MU', 'CALC_MU', 'MU_CALC', 'ABSORPTION', 'ABSORPTION_CALC', 'CALC_ABSORPTION')
       if(keyword_MU) call Absorption_routine	  
	   
     CASE ('SITE_INFO', 'LIST_SITE_INFO')
        !site_info_all_atoms = .true.
        IF(nb_atom /=0 .AND. keyword_SPGR) then
		 if(write_PCR_site_info) then
		  call get_site_info('pcr')
		 elseif(write_PCR_mag_site_info) then
		  call get_site_info('pcr_mag')
		 else
      	  call get_SITE_info('')
		 end if
		end if 

      case ('SIZE', 'CRYSTAL_SIZE')
        if (keyword_SIZE) then
         call get_crystal_SIZE_min_max
         call crystal_volume_calculation
        endif


      CASE ('SPGR',                    'SG',               'SPACE_GROUP',                                          &
            'SG_INFO',                 'SP_INFO',          'SPACE_GROUP_INFO',            'LIST_SPACE_GROUP_INFO', &
            'SG_EXTI',                 'SP_EXTI',          'SG_EXTINCTIONS',              'SPACE_GROUP_EXTI',      &
            'SPACE_GROUP_EXTINCTIONS', 'LIST_EXTINCTIONS', 'LIST_SPACE_GROUP_EXINCTIONS',                          &
            'SG_ALL',                  'SP_ALL',           'SG_INFO_EXTI',                'SP_INFO_EXTI',          &
            'SG_EXTI_INFO',            'SP_EXTI_INFO',                                                             &
            'SG_SUB',                  'SG_SUBGROUP')
        call space_group_info
        if (keyword_FILE) call get_numref_SPG

      CASE ('SYMM', 'SYM', 'SYMMETRY_OPERATOR')                   ! entree d'un operateur de symetrie
        call decode_sym_op()        ! cryscalc_symm.F90

      case ('TRANSLATION', 'TRANS', 'TRANSLATE', 'MOVE')
        IF(nb_atom /=0)  call transl_coord()

      CASE ('THERM', 'THERMAL', 'ADP', 'THERM_SHELX', 'THERMAL_SHELX', 'ADP_SHELX')     ! conversion des parametres d'agitation thermique
       IF(keyword_THERM) then
        IF(THERM_aniso) then
         call calc_therm_aniso
        else
         call calc_therm_iso
        endif
       endif

      CASE ('GEN_HKL', 'GENERATE_HKL', 'GENERATE_HKL_LIST')   ! generation des reflections
	   if(keyword_GENHKL) call generate_HKL()
	   
	  CASE ('PDP', 'PATTERN', '&PDP', '&PATTERN')   ! generation des reflections
        if(keyword_GENHKL) then
		 tmp_real = wavelength
		 wavelength = X_target(3)%wave(1)   ! Ka1_Cu
 		 call generate_HKL()
		 wavelength = tmp_real
		endif 

	  CASE ('PDP_CU', 'PATTERN_CU', '&PDP_CU', '&PATTERN_CU')   ! generation des reflections	    
        if(keyword_GENHKL) then
		 tmp_real = wavelength
		 wavelength = X_target(3)%wave(1)   ! Ka1_Cu
		 call generate_HKL()
		 wavelength = tmp_real
		endif 

      CASE ('ANG', 'ANGLE')
        IF (nb_ang_calc /=0  .and. keyword_CELL)  call calcul_angles             ! angles

      CASE ('DIST', 'DISTANCE', 'ATOMIC_DISTANCE')            ! distances
        IF (nb_dist_calc /=0 .and. keyword_CELL)  call calcul_distances

      CASE ('DIST_', 'DIST_X', 'DIST_MULT')   ! distances
        IF (keyword_DIST_X .and. nb_dist_calc /=0 .and. keyword_CELL)  call calcul_distances
		
      CASE ('DIST_PLUS', 'DIST_+')            ! distances
        IF (keyword_DIST_plus .and. nb_dist_calc /=0 .and. keyword_CELL)  call calcul_distances
		
      CASE ('DIST_DHA', 'DHA', 'POS_H', 'CALC_POS_H')		
        IF (keyword_DHA .and. nb_dist_calc /=0 .and. keyword_CELL)  call calcul_distances

      case ('CONN', 'CONNECT', 'CONNECTIVITY')  ! calcul connectivite
        if (keyword_CONN .and. keyword_CELL .and. keyword_SPGR) call calcul_connectivity_atom

      CASE ('BARY', 'CENTROID')
        IF (nb_bary_calc /=0)  call calcul_barycentre         ! barycentres

      case ('DIR_ANG', 'DIRANG', 'DIRECT_ANGLE')
        if (nb_da /=0)  call calcul_vector_angle('direct')

      case ('REC_ANG', 'RECANG', 'RECIPROCAL_ANGLE')
        if (nb_ra /=0)  call calcul_vector_angle('reciprocal')

      case ('STL', 'STL_HKL', 'SINTHETA/WAVE', 'SINTHETA/LAMBDA')
       if (nb_stl_value /=0)  call X_space_calculation('STL')

      case ('D_HKL', 'DHKL')
       if (nb_dhkl_value /=0) call X_space_calculation('DHKL')

      CASE ('HKL')
       if(nb_hkl /=0 .and. keyword_CELL)    call calcul_dhkl

      CASE ('SF_HKL', 'SFAC_HKL')
       if(nb_hkl /=0 .and. keyword_CELL .and. keyword_SPGR ) call Calcul_SFAC_HKL

      case ('D_STAR', 'D_STAR_HKL', 'DSTAR', 'DSTARHKL', 'DSTAR_HKL')
       if (nb_dstar_value /=0) call X_space_calculation('DSTAR')

      case ('Q_HKL', 'QHKL')
       if (nb_Qhkl_value /=0) call X_space_calculation('QHKL')

      case ('THETA', 'TH', 'TH_HKL', 'THETAHKL', 'THETA_HKL')
       if (nb_theta_value /=0) call X_space_calculation('THETA')

      case ('2THETA', '2TH', '2TH_HKL', '2THETA_HKL', 'TWO_THETA', 'TWO_THETA_HKL')
       if (nb_2theta_value /=0) call X_space_calculation('2THETA')


      CASE ('FILE', 'FIC', '&FIC', 'FSD', '&FSD', 'FST', '&FST', 'FIT', '&FIT')                              ! lecture fichier.HKL ou .CIF
       IF(keyword_FILE) then
	    call allocate_HKL_arrays
        call def_HKL_rule()
        call read_and_sort_hkl('sort')
        ! cas d'un fichier import.cif : calcul systematique du Rint -----
        ! si presence du champ '_symmetry_cell_setting'
        if (crystal_system(1:1) /= '?') then
         call get_SPG(crystal_system)
         !call calcul_global_Rint()
        endif
        ! ---------------------------------------------------------------
       endif
       IF(keyword_SPGR) call get_numref_SPG()

      case ('SEARCH_EXTI', 'FIND_EXTI')                       ! recherche des extinctions
        call search_exti()
        
      case ("SEARCH_SPGR", "SEARCH_SPACE_GROUP", "SEARCH_GROUP",   &
            "CHECK_SPGR",  "CHECK_SPACE_GROUP",  "CHECK_GROUP",    &
			"S_G", "&S_G", "GSG", "&GSG")
        call search_SPGR()


      case ('FIND_HKL', 'SEARCH_HKL')
        call search_hkl()                                     ! recherche reflection particuliere

      case ('EQUIV', 'EQUIV_HKL')
       call search_hkl_EQUIV()

      !case ('LIST_HKL', 'LST_HKL', 'LIST_HKL_EXTI', 'LST_HKL_EXTI')
      !CASE ("WRITE_HKL", "WRITE_HKL_LIST", "WRITE_HKL_LST")
       CASE ("FIND_HKL_LIST", "FIND_HKL_LST", "EXTRACT_HKL_LIST", "EXTRACT_HKL_LST") 
        if (HKL_list%EXTI_number/=0) call search_HKL_list()     ! recherche reflections obeissant à une règle d'extinctions

      case ('LIST_EXTI', 'LIST_EXTI_RULE', 'LST_EXTI', 'LST_EXTI_RULE')
        call list_EXTI_RULE()

      case ('HKL_POS', 'POS_HKL', 'HKL_POSITIVE', 'POSITIVE_HKL')
        if (keyword_FIND_HKL_POS) call search_HKL_F2('POS')

      case ('HKL_NEG', 'NEG_HKL', 'HKL_NEGATIVE', 'NEGATIVE_HKL')
       if (keyword_FIND_HKL_NEG)  call search_HKL_F2('NEG')

      case ('ABSENT_HKL', 'HKL_ABSENT')
        IF(keyword_FIND_HKL_ABSENT) call search_HKL_F2('ABSENT')

      case ('RINT', 'R_INT')
        call calcul_global_Rint("")

      case ('MERGE')
        call calcul_global_Rint('merge')
        !call calcul_merge_HKL
		
      case ('FRIEDEL', 'FRIEDEL_PAIRS')
        call get_Friedel_pairs_number()	  
		
      case ('SHELL', 'SD7', 'SD7', 'ST25', '&ST25')
       IF(nb_shell /=0) call read_and_sort_hkl('shell')

      case ('SORT')
       IF(nb_sort /=0)  call read_and_sort_hkl('sort')

      CASE ('OBV_REV', 'OBVERSE_REVERSE', 'TWIN_OBVERSE_REVERSE', 'TWIN_OBV_REV', &
            'TWINNING_OBVERSE_REVERSE', 'TWINNING_OBV_REV')
        if (keyword_FILE) call twin_obverse_reverse

      CASE ('TRICLINIC', 'TRICL')
       call transf_triclinic  


      CASE ('MONOCLINIC', 'MONOC', 'MONOCL')
        call transf_monoclinic

      CASE ('HEXA_TWIN', 'HEXA_TWINNING',  'HEXAGONAL_TWIN', 'HEXAGONAL_TWINNING', &
            'TWIN_HEXA', 'TWIN_HEXAGONAL', 'TWINNING_HEXA',  'TWINNING_HEXAGONAL')
        call twin_hexa
        
      CASE ('TWIN_PSEUDO_HEXA')
        call twin_pseudo_hexa   

      case ('PAUSE')
        pause

      case ('PERMUT', 'PERMUTATION', 'PERMUT_ABC', 'PERMUTATION_ABC')
        call permutation_abc

      CASE ('RHOMB_HEX', 'RHOMB_HEXA', 'RHOMB_TO_HEX', 'RHOMB_TO_HEXA')
        call transf_rhomb_to_hex

      CASE ('HEX_RHOMB', 'HEXA_RHOMB', 'HEX_TO_RHOMB', 'HEXA_TO_RHOMB')
        call transf_hex_to_rhomb

      CASE ('SYST', 'CMD', 'COMMAND', 'DOS', 'DOS_COMMAND', 'DIR')
        !if (LEN_TRIM(SYST_command)/=0)     call system(TRIM(SYST_command))
        if (LEN_TRIM(SYST_command)/=0)     CALL spawn (file_name=TRIM(SYST_command))  ! << NOK

      case ('RESET', 'RAZ', 'INITIALIZE', 'INIT')
        if (keyword_RESET) then
         tmp_logical = mode_interactif
         call Deallocate_HKL_arrays
		 call deallocate_PGF_data_arrays
		 call Def_transformation_matrix
         call cryscalc_INIT         
         call read_cryscalc_ini()
		 mode_interactif = tmp_logical
        endif 


      case ('MENDEL')
        IF(mendel_atom_nb /=0) call write_atomic_features()            ! mendel.F90

      case ('SHAN', 'SHANNON')
        call write_shannon_lines()

      case ('MAG', 'MAGNETIC', 'MAGNETISM')
        call write_mag_lines()

      case ('DATA_NEUTRONS', 'NEUTRONS_DATA', 'DATA_NEUTRON', 'NEUTRON_DATA')
        call write_data('neutrons')      ! mendel

      case ('DATA_XRAYS', 'XRAYS_DATA', 'DATA_XRAY', 'XRAY_DATA')
        call write_data('xrays')         ! mendel

      case ('DATA_DENSITY', 'DENSITY_DATA', 'DATA_ATOMIC_DENSITY', 'ATOMIC_DENSITY')
        call write_data('density')

      case ('DATA_WEIGHT',  'WEIGHT_DATA',  'DATA_ATOMIC_WEIGHT',  'ATOMIC_WEIGHT')
        call write_data('weight')

      case ('DATA_RADIUS',  'RADIUS_DATA',  'DATA_ATOMIC_RADIUS',  'ATOMIC_RADIUS')
        call write_data('radius')

      case ('WEB', 'INTERNET')
        if (keyword_WEB) call launch_browser(TRIM(URL_address))


      case ('WRITE_CELL', 'OUTPUT_CELL' , 'WC', '&WC')	    ! oct. 2011 
	    IF(.NOT. keyword_CELL) then
         call write_info('')
         call write_info('  Cell parameters has to be known for WRITE_CELL keyword')
         call write_info('')
         return
        endif
        !IF (keyword_CELL) call volume_calculation('out')   ! commenté oct. 2011
		if(.not. write_cell_cart) then
		 call volume_calculation('out')
		else
		 call write_crystal_cell_Cart
        endif		


      case ('WRITE_CHEM', 'WRITE_CHEMICAL_FORMULA', 'WRITE_MOLECULE', 'OUTPUT_CHEM', 'OUTPUT_CHEMICAL_FORMULA',  'OUTPUT_MOLECULE')
        call write_molecular_features

      case ('WRITE_QVEC', 'OUTPUT_QVEC')
        if (keyword_WRITE_QVEC) call write_QVEC
        
      case ('WRITE_SG', 'WRITE_SPACE_GROUP')
        !if (keyword_SPGR) call write_space_group(SPG%NumSPG,SPG%NumSPG )
        !if (keyword_SPGR) call write_current_space_group(SPG%NumSPG)
        if (keyword_SPGR) call write_current_space_group(space_group_symbol)

      case ('WRITE_WAVE', 'OUTPUT_WAVE')
        if (keyword_WAVE) call write_wave_features()
        
      case ('WRITE_DEVICE', 'OUTPUT_DEVICE')
        call write_DEVICE_features
  
      case ('WRITE_BEAM', 'WRITE_INCIDENT_BEAM', 'OUTPUT_BEAM', 'OUTPUT_INCIDENT_BEAM')
        call write_beam_features

! references --------------------------!		
      case ('REF_KCCD', 'KCCD')
        call write_REF('KCCD')

      case ('REF_APEX', 'REF_APEXII' , 'WRITE_APEX', 'WRITE_APEXII',  'APEX', 'APEXII')
        call write_REF('APEX')

      case ('REF_EVAL', 'REF_EVALCCD', 'WRITE_EVAL', 'WRITE_EVALCCD', 'EVAL', 'EVALCCD')
        call write_REF('EVAL')

      case ('REF_DENZO', 'WRITE_DENZO', 'DENZO')
        call write_REF('DENZO')
        
      case ('REF_SADABS', 'REF_SAD', 'SADABS')
        call write_REF('SADABS')
	  
	  case ('REF_ABS_CRYSALIS')
        call write_REF('ABS_CRYSALIS')
      
	  case ('REF_X2S', 'REF_SMART_X2S')
        call write_REF('X2S')	  
		
      case ('REF_SUPERNOVA')
        call write_REF('SUPERNOVA')	  
		
	  case ('REF_XCALIBUR')
        call write_REF('XCALIBUR')	  
		
!------------------------------------------------		


      case ('READ_NREPORT', 'READ_NREPORT_HTML', 'READ_HTMLREPORT')
        IF (keyword_read_NREPORT) call read_nreport_HTML() 

      case ('READ_CEL', 'READ_CEL_FILE', 'READ_POWDERCELL')
        IF (keyword_read_CEL) then
         open(UNIT=CEL_read_unit, FILE=trim(CEL_file_name), ACTION='read')
         call read_CEL_input_file(TRIM(CEL_file_name))
         close(UNIT=CEL_read_unit)
		 !if (create_CIF_PYMOL) call Create_CIF_file_for_pymol   ! creation fichier CIF compatible pymol
		 if (create_CIF_PYMOL) call create_CIF_PYMOL_file(TRIM(CEL_file_name), 'cel')
        endif

      case ('READ_CIF', 'READ_CIF_FILE', 'CIF_FILE')
        if (keyword_read_CIF) THEN
         OPEN(UNIT=CIF_read_unit, FILE=TRIM(CIF_file_name), ACTION="read")
         call check_CIF_input_file(TRIM(CIF_file_name))
         call read_CIF_input_file(TRIM(CIF_file_name), '?')
         call read_CIF_input_file_TR(CIF_read_unit)
         CLOSE(UNIT=CIF_read_unit)
		 !if (create_CIF_PYMOL) call Create_CIF_file_for_pymol   ! creation fichier CIF compatible pymol
		 if (create_CIF_PYMOL) call create_CIF_PYMOL_file(TRIM(CIF_file_name), 'cif')
         if (keyword_SPGR )    call space_group_info
         if (keyword_SYMM)     CALL decode_sym_op
        end if
		
      case ('READ_FACES', 'FACES')
	   if(keyword_read_FACES) then
	    open(unit=INS_read_unit, file=trim(FACES_file_name), action='read')
		call read_FACES_file(trim(FACES_file_name))
		close(unit=INS_read_unit)
	   end if
		
      case ('READ_INS', 'READ_INS_FILE', 'INS_FILE')
        if (keyword_read_INS) THEN
         OPEN(UNIT=INS_read_unit, FILE=TRIM(INS_file_name), ACTION="read")
         call read_INS_input_file(TRIM(INS_file_name), "ALL")
         call read_INS_SHELX_lines
         !call read_input_file_KEYWORDS(TRIM(INS_file_name))
         CLOSE(UNIT=INS_read_unit)
		 !if (create_CIF_PYMOL) call Create_CIF_file_for_pymol ! creation fichier CIF compatible pymol
		 if (create_CIF_PYMOL) call create_CIF_PYMOL_file(TRIM(INS_file_name), 'ins')
         if (keyword_SPGR )    call space_group_info
         if (keyword_SYMM)     CALL decode_sym_op
        end if

      case ('READ_PCR', 'READ_PCR_FILE', 'PCR_FILE')
        if (keyword_read_PCR) THEN
         OPEN(UNIT=PCR_read_unit, FILE=TRIM(PCR_file_name), ACTION="read")
         call read_PCR_input_file()
         CLOSE(UNIT=PCR_read_unit)
		 !if (create_CIF_PYMOL) call Create_CIF_file_for_pymol ! creation fichier CIF compatible pymol
		 if (create_CIF_PYMOL) call create_CIF_PYMOL_file(TRIM(PCR_file_name), 'pcr')
         if (keyword_SPGR )    call space_group_info
         if (keyword_SYMM)     CALL decode_sym_op
        end if

      case ('READ_TIDY_OUT')
       if(keyword_read_TIDY_out) then
	    OPEN(unit=TIDY_read_unit, file=trim(TIDY_out_file_name), ACTION="read")
		 call read_TIDY_out_input_file(trim(TIDY_out_file_name))
		close (unit=TIDY_read_unit)
		if (keyword_SPGR )    call space_group_info        
       end if
	   
      case ('SIR', 'SIR97')
        if (keyword_SIR) call create_SIR_file

		
      case('VER', 'VERSION')
        if(keyword_VERSION) call write_version

      case ('XRAYS_WAVELENGTH', 'X_WAVE')
        if (keyword_X_WAVE) CALL write_Xrays_wavelength

      case default

  end select
 !END do

 !close (UNIT=CFL_unit)
 !close (UNIT=11)

END subroutine run_keyword_interactive


!------------------------------------------------------------------------------

