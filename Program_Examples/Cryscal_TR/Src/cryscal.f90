!     Last change:  TR   10 Oct 2007    6:54 pm

!------------------------------------------------------------------------
subroutine read_keywords_from_file(arg_file)
 USE macros_module
 USE cryscal_module, ONLY : input_unit, nb_help, nb_help_max, HELP_string, message_text,                      &
                            input_INS, input_CFL, input_CIF, input_PCR, P4P_file_name, input_line,            &
                            keyword_modif_archive, keyword_WRITE_REF_DENZO, nb_atoms_type,                    &
                            keyword_create_CEL, keyword_create_CFL, keyword_create_INS, keyword_create_CIF,   &
                            keyword_create_ACE, keyword_create_FST,                                           &
                            keyword_read_CIF, keyword_read_PCR, keyword_read_INS,                             &
                            ACE_file_name, CIF_file_name, PCR_file_name, INS_file_name,                       &
                            PCR_read_unit, INS_read_unit, CFL_read_unit

 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)            :: arg_file
  CHARACTER (LEN=256)                      :: input_file
  CHARACTER (LEN=256)                      :: extension, arg_string
  LOGICAL                                  :: CFL_file
  LOGICAL                                  :: file_exist
  INTEGER                                  :: long
  INTEGER                                  :: i, nb_col
  LOGICAL                                  :: step_by_step
  LOGICAL                                  :: lecture_ok

  step_by_step = .false.
  if (LEN_TRIM(arg_file) ==0) then    ! option #1 du menu principal
   do
    input_file = ''
    call write_info(' ')
    call write_info(' >> Enter input file: ')
    call read_input_line(input_line)
    READ(input_line,'(a)') input_file
    input_file = u_case(input_file)
    input_file = ADJUSTL(input_file)
 
    IF(input_file(1:4) == 'HELP' .or. input_file(1:3) == 'MAN') then
     read(input_file(5:), '(a)') arg_string
     call nombre_de_colonnes(arg_string, nb_col)
     IF(nb_col /=0) then
      nb_help = nb_col
      read(input_file(5:), *) HELP_string(1:nb_help)
     END if
     !call write_HELP()
     IF(nb_help == nb_help_max) call write_header
     call HELP_on_line
     return
    ELSEIF(input_file(1:3) == 'KEY') then
     call write_KEYWORD('screen')
    elseif(LEN_TRIM(input_file) == 0) then
     !call write_header()
     !stop
     return
    endif

    i = index(input_file,'.')
    if (i==0) input_file = trim(input_file)//'.cfl'
    call test_file_exist(input_file, file_exist)
    IF(.NOT. file_exist) then    
     input_file = ''
     cycle
    endif
    exit
   END do

  else
   input_file   = arg_file
   step_by_step = .true.
  endif


 ! type of input file
 input_file = l_case(input_file)
 i = INDEX(input_file, '.')
 long = LEN_TRIM(input_file)
 extension = input_file(i+1:long)

 

 call write_info(' ')
 call write_info('  Input CRYSCAL file: '//trim(input_file))
 call write_info(' ')
 
 OPEN(UNIT=input_unit, FILE=TRIM(input_file), ACTION="read")


 select case (extension(1:3))

     case ("ins", "res")
      input_INS = .true.
      open(UNIT=INS_read_unit, FILE=TRIM(input_file), ACTION='read')
      call read_INS_input_file(TRIM(input_file), "ALL")
      call read_input_file_KEYWORDS(TRIM(input_file))
      if (keyword_create_CFL) call create_CFL_file(TRIM(input_file), extension)
      if (keyword_create_FST) call create_FST_file(TRIM(input_file), extension)
      if(keyword_create_CEL) then
       keyword_read_INS = .true.
       INS_file_name = input_file
       call create_CEL_from_CIF
      endif
      if(keyword_create_ACE) then
       keyword_read_INS = .true.
       INS_file_name = input_file
       call create_ACE_from_CIF
      endif


     case ("pcr")
      input_PCR = .true.
      OPEN(UNIT=PCR_read_unit, FILE=TRIM(input_file), ACTION="read")
      call read_PCR_input_file(TRIM(input_file))
      if (keyword_create_CFL)  call create_CFL_file(TRIM(input_file), extension)
      if (keyword_create_FST)  call create_FST_file(TRIM(input_file), extension)
      if(keyword_create_CEL) then
       keyword_read_PCR = .true.
       PCR_file_name = input_file
       call create_CEL_from_CIF
      endif
      if(keyword_create_ACE) then
       keyword_read_PCR = .true.
       PCR_file_name = input_file
       call create_ACE_from_CIF
      endif

     case ("cfl")
      input_CFL = .true.
      IF(step_by_step) then             ! a partir de la ligne de commande
       call interactive_mode('file')
      else                              ! option #1 du menu principal
        !call read_CFL_input_file(TRIM(input_file))              ! read_CFL.F90
       call read_CFL_input_file()              ! read_CFL.F90
       call read_input_file_KEYWORDS(TRIM(input_file))         ! read_KEYW.F90
      endif

     case ("cif")
      input_CIF = .true.
      call check_CIF_input_file(TRIM(input_file))
      IF(input_file(:) /= 'archive.cif') then
       call get_P4P_file_name(P4P_file_name)
       IF(LEN_TRIM(P4P_file_name) /=0) call read_P4P_file(P4P_file_name, lecture_OK)
      endif


      call read_CIF_input_file(TRIM(input_file), '?')          ! read_CIF_file.F90
      call read_CIF_input_file_TR(input_unit)                  ! read_CIF_file.F90
      !call read_SQUEEZE_file
      !call read_input_file_KEYWORDS(TRIM(input_file))
      
      !if(nb_atoms_type /=0) then
       if(keyword_create_CFL) call create_CFL_file(TRIM(input_file), extension)
       if(keyword_create_INS) call create_INS_file(TRIM(input_file), extension)  
       if(keyword_create_FST) call create_FST_file(TRIM(input_file), extension)     
      !endif
     
     
      if(keyword_create_ACE) then
       keyword_read_CIF = .true.
       CIF_file_name = input_file
       call create_ACE_from_CIF
      endif

      if(keyword_create_CEL) then
       keyword_read_CIF = .true.
       CIF_file_name = input_file
       call create_CEL_from_CIF
      endif

      !if(keyword_modif_ARCHIVE) call read_and_modif_archive(input_unit)

     case default
      !call read_CFL_input_file(TRIM(input_file))
      call read_CFL_input_file(input_unit)
      call read_input_file_KEYWORDS(TRIM(input_file))

 end select

end subroutine read_keywords_from_file

!-------------------------------------------------------------------
subroutine run_keywords()
 USE cryscal_module
 USE HKL_module
 USE IO_module

 if(ON_SCREEN) then
  call write_info('')
  call write_info('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  endif

 IF(keyword_setting)                       call output_setting()              ! cryscal.init.F90
 if(keyword_HELP)                          call HELP_on_line         !         help.F90
 IF(keyword_HEADER)                        call write_header         !         help.F90
 IF(keyword_KEY)                           call KEYS_on_line         !      cryscal.F90
 IF(keyword_create_CRYSCAL_HTML)           call create_CRYSCAL_HTML           !  create_HTML.F90

 IF(keyword_CELL)                          call volume_calculation('out')   !       calculs.F90
 IF(keyword_NIGGLI)                        call Niggli_cell_CFML()          !        niggli.F90
 IF(keyword_WRITE_CELL)                    call volume_calculation('out')   !       calculs.F90
 IF(keyword_WRITE_QVEC)                    call write_QVEC()                ! cryscal_write.F90
 if(keyword_WRITE_WAVE .and. keyword_WAVE) call write_WAVE_features()       ! cryscal_write.F90
 if(keyword_WRITE_BEAM .and. keyword_BEAM) call write_BEAM_features()       ! crsycal_write.F90     
 if(keyword_WRITE_DEVICE)                  call write_DEVICE_features()     ! cryscal_write.F90   
 IF(keyword_WRITE_SG   .AND. keyword_SPGR) call space_group_info()          !   space_group.F90

 IF(keyword_SIZE)                          call crystal_volume_calculation('out') ! calculs.F90

 IF(keyword_SFAC_UNIT .or. keyword_CONT .or. keyword_CHEM)  then
  call atomic_identification()                                              ! calculs.F90
  IF(keyword_CELL)                      call atomic_density_calculation()   ! calculs.F90
  IF(keyword_ZUNIT)	                call molecular_weight()             ! calculs.F90
  IF(keyword_CELL .and. keyword_ZUNIT)  call density_calculation()          ! calculs.F90
  IF(keyword_CELL)                      call absorption_calculation()       ! mu_calc.F90
 END if


 ! liste groupes d'espace
 IF(keyword_LSPGR)  call list_space_groups()              ! cryscal_lsg.F90


 ! calcul des d_ 
 if(nb_hkl /=0           .and. keyword_CELL)        call calcul_dhkl
 if(nb_hkl_SFAC_calc /=0 .and. keyword_SFAC_HKL)    call calcul_SFAC_hkl


 ! transformation de la maille
 if (keyword_MATR) then
  call WRITE_matrice()
  if (keyword_CELL)  call transf_cell_parameters
  if (nb_hkl  /=0)   call transf_HKL
  if (keyword_FILE)  call transf_HKL_file
  call transf_coord()
  IF(input_INS .and. ABS(Mat_det) -1. < 0.001 ) call create_NEW_INS()
 end if

 IF(keyword_LST_MAT) call write_list_matrice()

 IF(keyword_TRANSL    .and. nb_atom /=0)  call transl_coord

 IF(keyword_INSIDE    .AND. nb_atom /=0)  call inside_unit_cell

 IF(keyword_ATOM_list .AND. nb_atom /=0)  call write_atom_list

 if (keyword_SPGR )    call space_group_info                 ! space_group.F90

 ! generation des equivalents
 if (keyword_SYMM ) call decode_sym_op()

 ! conversion des parametres d'agitation thermique
 if (keyword_THERM)   then
  IF(THERM_aniso) then
   call calc_THERM_aniso()
  else
   call calc_THERM_iso
  endif
 endif

 ! generation des reflections
 IF(keyword_GENHKL)   call generate_HKL()

 ! distance calculation between atom_1 and atom_2 ---------------------
 IF(nb_dist_calc /=0 .and. keyword_CELL)  call calcul_distances

 ! angles
 if (nb_ang_calc /=0 .and. keyword_CELL)  call calcul_angles

 ! barycentres
 if (nb_bary_calc /=0) call calcul_barycentre

 if(keyword_CONN .and. keyword_CELL .and. keyword_SPGR) call calcul_connectivity_atom

 ! calcul des X_space
 if (keyword_STL    .AND. nb_stl_value    /=0)  call X_space_calculation('STL')
 if (keyword_dhkl   .AND. nb_dhkl_value   /=0)  call X_space_calculation('DHKL')
 if (keyword_dstar  .AND. nb_dstar_value  /=0)  call X_space_calculation('DSTAR')
 if (keyword_Qhkl   .AND. nb_Qhkl_value   /=0)  call X_space_calculation('QHKL')
 if (keyword_theta  .AND. nb_theta_value  /=0)  call X_space_calculation('THETA')
 if (keyword_2theta .AND. nb_2theta_value /=0)  call X_space_calculation('2THETA')

 ! calcul des angles entre 2 vecteurs
 if (nb_da /=0)  call calcul_vector_angle('direct')
 if (nb_ra /=0)  call calcul_vector_angle('reciprocal')

 ! tri fichier.HKL
 !if (keyword_FILE  .and. nb_sort /=0) call read_and_sort_hkl()
 if (keyword_FILE) then
  call def_HKL_rule()
  if (nb_shell /=0) then
   call read_and_sort_hkl('shell')
  else
   call read_and_sort_hkl('sort')
  endif
!  if (nb_shell /=0)       call read_and_sort_hkl('shell')
!  if (nb_sort  /=0)       call read_and_sort_hkl('sort')
  IF(keyword_RINT)        call calcul_global_Rint()

  IF(keyword_search_exti) call search_exti()     ! recherche extinctions systematiques
  IF(keyword_search_SPGR) call search_SPGR()     ! recherche d'un groupe d'espace
  IF(keyword_FIND_HKL)    call search_hkl()      ! recherche reflection hkl particuliere
  IF(keyword_FIND_HKL_list .and. HKL_list%EXTI_number/=0)  call search_HKL_list() ! recherche reflections obeissant à une règle d'extinctions
  IF(keyword_FIND_HKL_pos)     call search_HKL_F2('POS')
  IF(keyword_FIND_HKL_neg)     call search_HKL_F2('NEG')
  IF(keyword_FIND_HKL_ABSENT)  call search_HKL_F2('ABSENT')

  if(keyword_OBV_REV)     call twin_obverse_reverse    ! obv_rev.F90

  IF(keyword_SPGR)        call get_numref_SPG()
 end if

 if (WRITE_triclinic_transf)  call transf_triclinic      ! transf.F90
 if (WRITE_monoclinic_transf) call transf_monoclinic     ! transf.F90
 if (WRITE_rhomb_hex_transf)  call transf_rhomb_to_hex   ! transf.F90
 if (WRITE_hex_rhomb_transf)  call transf_hex_to_rhomb   ! transf.F90
 if (WRITE_permutation_abc)   CALL permutation_abc       ! transf.F90
 if (WRITE_twin_hexa)         CALL twin_hexa             ! transf.F90
 if (WRITE_twin_pseudo_hexa)  CALL twin_pseudo_hexa      ! transf.F90


 ! keyword_SYST
 if (keyword_SYST)            call system(TRIM(SYST_command))

 IF(keyword_MENDEL .AND. mendel_atom_nb/=0) call write_atomic_features  ! mendel.F90


 ! DATA
 if (keyword_DATA_neutrons)         call write_data('neutrons')      ! mendel
 if (keyword_DATA_xrays)            call write_data('xrays')         ! mendel
 if (keyword_DATA_atomic_density)   call write_data('density')
 if (keyword_DATA_atomic_radius)    call write_data('radius')
 if (keyword_DATA_atomic_weight)    call write_data('weight')

 if (keyword_X_WAVE)         CALL write_Xrays_wavelength

 call read_SQUEEZE_file


 ! CIF file for APEX
 IF(keyword_create_CIF) then
  IF(keyword_WRITE_REF_APEX)   call write_REF('APEX')
  IF(keyword_WRITE_REF_KCCD)   call write_REF('KCCD')
  IF(keyword_WRITE_REF_EVAL)   call write_REF('EVAL')
  IF(keyword_WRITE_REF_DENZO)  call write_REF('DENZO')
  IF(keyword_WRITE_REF_SADABS) call write_REF('SADABS')
 endif
 if(keyword_PAUSE) pause

end subroutine run_keywords


!------------------------------------------------------------------------

subroutine KEYS_on_line()
 USE cryscal_module, ONLY : nb_help_max, help_string, write_keys, message_text, KEYS_unit
 USE IO_module

 implicit none
  INTEGER                       :: i , nb_keys_to_write

  nb_keys_to_write = 0
  do i=1, nb_help_max
   IF(write_keys(i)) then
    nb_keys_to_write = nb_keys_to_write + 1
    !if(input_string(1:6) == 'screen') then
     IF(nb_keys_to_write == 1) then
      call write_info('')
      call write_info('   > CRYSCAL keywords list:')
      call write_info('')
     endif
     WRITE(message_text, '(5x,2a)') '. ', TRIM(help_string(i))
     call write_info(TRIM(message_text))
    !elseif(input_string(1:4) == 'file') then
    ! if(nb_keys_to_write == 1) then
    !  write(unit = KEYS_unit, '(a)') ''
    !  write(unit = KEYS_unit, '(a)') '   > CRYSCAL keywords list::'
    !  write(unit = KEYS_unit, '(a)') ''
    ! endif
    ! WRITE(message_text, '(5x,2a)') '. ', TRIM(help_string(i))
    ! write(unit = KEYS_unit, '(a)') trim(message_text)
    !endif
   endif
  end do

  IF(nb_keys_to_write == 0) then
   call write_info('')
   call write_info('  ... unknown argument for keyword KEYS ...')
   call write_info('')
  endif

 RETURN
end subroutine KEYS_on_line

!--------------------------------------------------------------------

subroutine check_CIF_input_file(input_file)
 USE cryscal_module, ONLY : keyword_create_CIF, CIF_unit, CIF_archive_unit, &
                            keyword_WRITE_REF_APEX, keyword_WRITE_REF_KCCD, &
                            AUTHOR, DEVICE
 USE text_module,    ONLY : CIF_lines_nb, CIF_title_line
 use macros_module,  only : u_case
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)       :: input_file
  character (len=256)                 :: CIF_string
  INTEGER                             :: i, i1, i2, i_error

  IF(input_file(:) == 'archive.cif') then
   keyword_create_CIF = .true.

  else
   i1 = INDEX(input_file, ':')
   IF(i1/=0) then
    i2 = INDEX(input_file, '\', back=.TRUE.)
    IF(i2 ==0) then
     IF(input_file(i1+1:) == 'archive.cif') keyword_create_CIF = .true.
    else
     IF(input_file(i2+1:) == 'archive.cif') keyword_create_CIF = .true.
    endif
   endif
  endif

  IF(keyword_create_CIF) then
   if(u_case(DEVICE%diffracto(1:4)) == 'APEX') then
    keyword_WRITE_REF_APEX =.true.
   elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
    keyword_WRITE_REF_KCCD = .true.
   endif
   OPEN(UNIT=CIF_ARCHIVE_unit, FILE = 'archive_cryscal.cif', iostat=i_error)
   if(i_error /=0) then
    call write_info('')
    call write_info(' Error opening archive_cryscal.cif file. Program will be stopped.')
    call write_info('')
    stop
   endif
   OPEN(UNIT=CIF_unit,         FILE = 'cryscal.cif',         iostat=i_error)
   if(i_error /=0) then
    call write_info('')
    call write_info(' Error opening cryscal.cif file. Program will be stopped.')
    call write_info('')
    stop
   endif

   call verif_CIF_character(AUTHOR%NAME)
   call verif_CIF_character(AUTHOR%first_name)
   call verif_CIF_character(AUTHOR%address)
   call verif_CIF_character(AUTHOR%email)
   call verif_CIF_character(AUTHOR%web)

   do i=1, CIF_lines_nb
    call write_CIF(CIF_unit, trim(CIF_title_line(i)))
   end do

   call write_CIF_AUTHOR(CIF_unit)

  endif

 return
END subroutine check_CIF_input_file

!--------------------------------------------------------------------

subroutine get_P4P_file_name()
 USE cryscal_module, ONLY : P4P_file_name
 implicit none
  INTEGER              :: i, ier, n_P4P
  CHARACTER (LEN=256)  :: read_line

  call system('dir *.p4p /B > p4p.lst')      !! message 'fichier introuvable' si pas de fichier.P4P !
  open (UNIT=9, FILE='p4p.lst')
  n_P4P = 0
  do
   READ(UNIT=9, '(a)', IOSTAT=ier) read_line
   IF(ier <0) EXIT ! fin du fichier
   n_p4p = n_p4p + 1
  END do
  CLOSE(UNIT=9)
  call system('del p4p.lst')

  IF(n_p4p == 1) then
   P4P_file_name = read_line
  else
   P4P_file_name = ''
  endif


 RETURN
end subroutine get_P4P_file_name

!--------------------------------------------------------------------
