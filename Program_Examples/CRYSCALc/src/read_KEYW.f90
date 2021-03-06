!     Last change:  TR    5 Sep 2007    6:21 pm
subroutine read_input_file_KEYWORDS(input_unit)
 USE cryscalc_module, ONLY : CFL_read_unit, keyword_BEAM, keyword_SFAC_UNIT, keyword_CONT, keyword_CHEM, keyword_ZUNIT, &
                             debug_proc
 USE macros_module,   ONLY : u_case
 USE IO_module,       ONLY : write_info
 implicit none
  integer, intent(in)                      :: input_unit
  CHARACTER (LEN=256)                      :: read_line
  INTEGER                                  :: i_error


  if(debug_proc%level_2)  call write_debug_proc_level(2, "read_input_file_KEYWORDS")

  REWIND(UNIT=input_unit)

 do        ! lecture du fichier d'entree
  READ(input_unit, '(a)', IOSTAT=i_error) read_line

  IF(i_error < 0) EXIT   ! fin du fichier
  if(i_error > 0) then
   call write_info('')
   call write_info('Error reading file')
   call write_info('')
   return
  end if
  read_line = ADJUSTL(read_line)
  read_line = u_case(TRIM(read_line))

  IF (LEN_TRIM(read_line) == 0)                           CYCLE  ! ligne vide
  IF (read_line(1:1) == '! ' .or. read_line(1:1) == '#' ) CYCLE  ! commentaire
  !IF (read_line(1:1) == '_')                              CYCLE  ! cas d'un fichier.CIF

  call identification_keywords(read_line)
 ! call run_keyword_interactive(current_keyword)

 END do


 IF(.NOT. keyword_BEAM .AND. ( keyword_SFAC_UNIT .or. keyword_CONT .or. keyword_CHEM)) then
  call incident_beam()
 endif

 IF(keyword_CHEM) then
  IF(.NOT. keyword_ZUNIT) then
   call write_info('')
   call write_info('  >> ZUNIT keyword is missing ! ')
   call write_info('')
   stop
  endif
 endif


  RETURN
end subroutine read_input_file_KEYWORDS


!----------------------------------------------------------------------
subroutine identification_keywords(read_line)
 USE cryscalc_module
 USE CIF_module
 USE Pattern_profile_module,    ONLY : size_broadening, particle_size, X_pattern, N_pattern
 USE hkl_module,                ONLY : n_sig, threshold, ratio_criteria, requested_H, requested_H_string, search_H_string,    &
                                       HKL_list, HKL_list_NEG, HKL_list_POS, HKL_file, ordered_HKL,                           &
                                       HKL_data_known, HKL_list_ABSENT, search_equiv, search_friedel, HKL_rule_nb

 use CFML_crystallographic_symmetry, only : set_spacegroup
 USE MATRIX_list_module
 USE macros_module
 USE IO_module,                 ONLY : write_info
 USE matrix_module,             ONLY : get_mat_from_setting
 USE wavelength_module
 USE atome_module
 USE text_module ,              ONLY : CIF_lines_nb, CIF_title_line

 implicit none
  CHARACTER (LEN=*), INTENT(IN)              :: read_line
  CHARACTER (LEN=32)                         :: input_keyword
  INTEGER                                    :: long_kw, long, long1, long2
  CHARACTER (LEN=256)                        :: arg_line, new_line, arg1
  CHARACTER (LEN=12), DIMENSION(nb_help_max) :: temp_string
  INTEGER                                    :: i, j, i1, i2, num, i_error, nb_arg
  INTEGER                                    :: nb_arg_sep, nb_sup
  INTEGER                                    :: current_i1, current_i2, i_ok
  INTEGER                                    :: i_CONN_dmin, i_CONN_dmax
  REAL,               DIMENSION(10)          :: var
  CHARACTER (LEN=64), DIMENSION(20)          :: arg_string, sup_string
  CHARACTER (LEN=64)                         :: arg
  LOGICAL                                    :: file_exist, string_ok
  LOGICAL                                    :: UB_mat_1i2i3i
  LOGICAL                                    :: bary_ALL, OK


  if(debug_proc%level_2)  call write_debug_proc_level(2, "identification_keywords ("//trim(read_line)//")")

  arg_string    = ''
  var           = 0.
  bary_ALL      = .false.
  !write_details = .true.

  READ(read_line, *) input_keyword
  input_keyword = ADJUSTL(input_keyword)
  long_kw = LEN_TRIM(input_keyword)
  if (input_keyword(long_kw:long_kw) == '=') then
   input_keyword = input_keyword(1:long_kw)
   long_kw = LEN_TRIM(input_keyword)
  end if

  READ(read_line(long_kw+1:),'(a)') arg_line
  arg_line = ADJUSTL(arg_line)

  long=LEN_TRIM(arg_line)
  !IF(arg_line(1:1) == "'" .AND. arg_line(long:long)=="'") arg_line = arg_line(2:long-1) ! chaine entre ''
  !call remove_car(arg_line, "'")    ! caractere "'" dans la chaine
  arg_line = remove_car(arg_line, "'")
  IF(arg_line(1:1) == '=') then
   arg_line = arg_line(2:)
   arg_line = ADJUSTL(arg_line)
  endif


  call nombre_de_colonnes(arg_line, nb_arg)
  IF(nb_arg/=0)  then
   !READ(arg_line, *, IOSTAT=i_error) arg_string(1:nb_arg)    !!  << pb si presence d'une , dans l'argument (ex: FMT=(3i4,2f8.2)
   !IF(i_error /=0) then
   ! call write_info('')
   ! WRITE(message_text, '(2a)') '  input line: ', TRIM(read_line)
   ! call write_info(TRIM(message_text))
   ! WRITE(message_text, '(3a)') '  ... Error reading ', TRIM(input_keyword), ' arguments ...'
   ! call write_info(TRIM(message_text))
   ! return
   !endif

   call get_arg_string(arg_line, nb_arg, arg_string)          !!  << new routine de lecture des arguments (fev. 2013)
   !do i=1, nb_arg
   ! arg_string(i) = u_case(arg_string(i))
   ! long = len_trim(arg_string(i))
   ! if(long >= 6) then
   ! if(arg_string(i)(1:6) == 'NO_OUT' .or. arg_string(i)(1:6) == 'NO_DET') then
     ! WRITE_details = .false.
     !end if
    !end if
    !if(long >=10) then
    ! if(arg_string(i)(1:10) == 'NO_DETAILS') then
    !  WRITE_details = .false.
    ! end if
    !end if
   !end do
  endif

  !if(.not. WRITE_details) then
  ! nb_arg = nb_arg - 1
  ! arg_line = remove_car(arg_line, 'NO_DETAILS')
  ! arg_line = remove_car(arg_line, 'NO_DET')
  ! arg_line = remove_car(arg_line, 'NO_OUT')
  ! call get_arg_string(arg_line, nb_arg, arg_string)
  ! do i=1, nb_arg
  !  arg_string(i) = u_case(arg_string(i))
  ! end do
  !end if


  select case (TRIM(input_keyword))

   case ('EXPERT_MODE', 'EXPERT', 'EXPERT_ON')
    expert_mode = .true.
	if (nb_arg /=0) then
	 long = len_trim(arg_string(1))
	 if(long == 2) then
	  if (arg_string(i)(1:2) == 'ON') expert_mode = .true.
	 elseif(long == 3) then
	  if (arg_string(i)(1:3) == 'OFF') expert_mode = .false.
	 end if
	end if

   case ('USER_MODE', 'USER', 'EXPERT_OFF')
    if(expert_mode) expert_mode = .false.

   case ("DEBUG", "DEBUG_ON", "DEBUG ON")
    if(expert_mode) then
     if(.not. debug_proc%write) call open_debug_file
     debug_proc%write   = .true.
     debug_proc%level_1 = .true.
     debug_proc%level_3 = .false.
    end if

   case ("&DEBUG", "&DEBUG_ON", "&DEBUG ON")
    if(.not. debug_proc%write) call open_debug_file
    debug_proc%write   = .true.
    debug_proc%level_1 = .true.
    debug_proc%level_3 = .false.

   case ("DEBUG_2", "DEBUG_2_ON", "DEBUG_2 ON", "DEBUG 2 ON")
    if(expert_mode) then
     if(.not. debug_proc%write) call open_debug_file
     debug_proc%write   = .true.
     debug_proc%level_2 = .true.
     debug_proc%level_3 = .false.
    end if

   case ("&DEBUG_2", "&DEBUG_2_ON", "&DEBUG_2 ON", "&DEBUG 2 ON")
    if(.not. debug_proc%write) call open_debug_file
    debug_proc%write   = .true.
    debug_proc%level_2 = .true.
    debug_proc%level_3 = .false.

   case ("DEBUG_3", "DEBUG_3_ON", "DEBUG_3 ON", "DEBUG 3 ON")
    if(expert_mode) then
     if(.not. debug_proc%write) call open_debug_file
     debug_proc%write   = .true.
     debug_proc%level_2 = .true.
     debug_proc%level_3 = .true.
    end if

   case ("&DEBUG_3", "&DEBUG_3_ON", "&DEBUG_3 ON", "&DEBUG 3 ON")
    if(.not. debug_proc%write) call open_debug_file
    debug_proc%write   = .true.
    debug_proc%level_2 = .true.
    debug_proc%level_3 = .true.


   case ("DEBUG_OFF", "DEBUG OFF")
    if(expert_mode) then
     debug_proc%write   = .false.
     debug_proc%level_1 = .false.
     debug_proc%level_2 = .false.
     debug_proc%level_3 = .false.
     close(unit = debug_proc%unit)
    end if

   case ("&DEBUG_OFF", "&DEBUG OFF")
     debug_proc%write   = .false.
     debug_proc%level_1 = .false.
     debug_proc%level_2 = .false.
     debug_proc%level_3 = .false.
     close(unit = debug_proc%unit)

   case ('NO_DETAILS', 'NO_DET', 'NO_OUT')
    if(expert_mode) then
     keyword_NO_details = .true.
     write_details      = .false.
    end if

   case ('&NO_DETAILS', '&NO_DET', '&NO_OUT')
    keyword_NO_details = .true.
    write_details      = .false.


   case ('ACTA', 'CIF', 'CREATE_CIF')
    keyword_create_CIF = .true.
    !OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
    !do i=1, CIF_lines_nb
    ! WRITE(CIF_unit, '(a)') trim(CIF_title_line(i))
    !end do
    !WRITE(CIF_unit, '(a)') ''
    !WRITE(CIF_unit, '(a)') 'data_cryscalc'
    !WRITE(CIF_unit, '(a)') ''

   case ('ACTA_OFF', 'CIF_OFF')
    keyword_create_CIF = .false.


   case ('CREATE_ACE')
    keyword_create_ACE = .true.

   case ('CREATE_CEL')
    keyword_create_CEL = .true.

   case ('CREATE_CFL')
    keyword_create_CFL = .true.

   case ('CREATE_FST')
    keyword_create_FST = .true.
    create_FST_poly    = .false.
    create_FST_mole    = .false.
    FST_no_H           = .false.
    launch_FP_STUDIO   = .false.
    if(nb_arg /=0) then
     do i=1, nb_arg
      if(len_trim(arg_string(i))== 4 .and. arg_string(i)(1:4) == "POLY") create_FST_POLY  = .true.
      if(len_trim(arg_string(i))== 4 .and. arg_string(i)(1:4) == "MOLE") create_FST_MOLE  = .true.
      if(len_trim(arg_string(i))== 4 .and. arg_string(i)(1:4) == "NO_H") FST_no_H         = .true.
      if(len_trim(arg_string(i))== 3 .and. arg_string(i)(1:3) == "RUN")  launch_FP_STUDIO = .true.
     end do
    end if

   !case ('CREATE_CIF_PYMOL')
   ! create_CIF_PYMOL = .true.

   case ('CREATE_INS')
    keyword_create_INS = .true.
    if(nb_arg > 0) then
     do i=1, nb_arg
      if(len_trim(arg_string(i)) == 4 .and. arg_string(i)(1:4) == 'NO_H') no_H = .true.
     end do
    end if

   case ('CREATE_PCR')
    keyword_create_PCR = .true.

   case ('CREATE_SOLVE', 'CREATE_TO_SOLVE', 'CREATE_FILES_TO_SOLVE', 'SOLVE')
    keyword_create_SOLVE = .true.

   case ('CREATE_TIDY', 'CREATE_TIDY_FILE', 'CREATE_TIDY_INPUT_FILE')
    keyword_create_TIDY = .true.

   CASE ('LST_ATOMS', 'ATOM_LST', 'ATOM_LIST', 'LIST_ATOM_LIST', 'LIST_ATOMS', 'WRITE_ATOMS', 'WRITE_ATMS')
    keyword_atom_list = .true.
    keyword_ADP_list  = .false.
    write_atoms_cartesian = .false.
    if(nb_arg /=0) then
     do i=1, nb_arg
      if(arg_string(i)(1:4) == 'CART') then
       write_atoms_cartesian = .true.
       if(len_trim(arg_string(i)) > 5 .and. arg_string(i)(5:5) == '_') then
        if (arg_string(i)(6:6) == "A") then
         cartesian_frame%type   = "A"
         cartesian_frame%string = "x // a"
        elseif (arg_string(i)(6:6) == "C") then
         cartesian_frame%type = "C"
         cartesian_frame%string = "x // c"
        end if
       end if
      !elseif(arg_string(i)(1:3) == 'SHP' .or. arg_string(i)(1:5) == 'SHAPE') then
      ! create_SHAPE_file = .true.
      elseif(arg_string(i)(1:4) == 'IN_A') then
       write_atoms_in_A = .true.
      end if
     end do
    end if

   CASE ('WRITE_ADP')
    keyword_ADP_list  = .true.
    keyword_atom_list = .false.
    ADP_details       = .false.
    if(nb_arg /= 0 .and. arg_string(1)(1:7) == 'DETAILS') ADP_details = .true.


   case ('EDIT')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'EDIT')
     return
    endif
    file_to_edit = arg_string(1)
    call test_file_exist(file_to_edit, file_exist, 'out')
    IF(.not. file_exist) then
     file_to_edit = ''
     return
    endif
    keyword_EDIT = .true.

   CASE ('FILE', 'READ_HKL')
    HKL_file%NAME      = ''
    HKL_file%plot      = .false.
    HKL_file%read_NEG  = .true.
    HKL_data_known     = .false.
    HKL_file%cif       = .false.
    HKL_file%SHELX     = .false.
    HKL_file%HKLF5     = .false.
    HKL_file%final_y   = .false.
    HKL_file%RAW       = .false.
    HKL_file%M91       = .false.
    HKL_file%M95       = .false.
    HKL_file%INT       = .false.
    HKL_file%COL       = .false.
    file_out           = .false.
    ordered_hkl        = .false.
    keyword_read_HKLF5 = .false.
    !WRITE_details      = .true.
    nb_shell = 0
    nb_sort  = 0

    IF(nb_arg == 0) then
     call check_arg_nb('0', 'FILE')
     return
    endif
    READ(arg_string(1), *) HKL_file%NAME
    !IF(nb_arg >1) then
    ! READ(arg_string(2), *) arg1
    ! arg1 = ADJUSTL(arg1)
    !endif

    i = INDEX(HKL_file%NAME, '.')
    if (i==0) then
     HKL_file%NAME = TRIM(HKL_file%NAME)//'.hkl'
     HKL_file%SHELX = .true.
    else
     long = LEN_TRIM(HKL_file%NAME)
     IF(HKL_file%NAME(long-2:long) == 'HKL')     HKL_file%SHELX   = .true.
     IF(HKL_file%NAME(long-2:long) == 'CIF')     HKL_file%CIF     = .true.
     IF(HKL_file%NAME(long-2:long) == 'RAW')     HKL_file%RAW     = .true.
     IF(HKL_file%NAME(long-2:long) == 'M91')     HKL_file%M91     = .true.
     IF(HKL_file%NAME(long-2:long) == 'M95')     HKL_file%M95     = .true.
     ! IF(HKL_file%NAME(1:long)      == 'FINAL.Y') HKL_file%final_y = .true.  ! feb. 2013 : final_y ne contient pas les F2 !!
     IF(HKL_file%NAME(long-2:long) == 'INT')     HKL_file%INT     = .true.
     IF(HKL_file%NAME(long-2:long) == 'COL')     HKL_file%COL     = .true.
    endif
    call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
    IF(.NOT. file_exist) then
     HKL_file%name = ''
     return
    endif
    keyword_FILE = .true.

    hkl_format_SHELX = .true.
    hkl_format = hkl_SHELX_format

    !if(HKL_file%SHELX) hkl_format = hkl_SHELX_format

    do i=2, nb_arg
     IF(arg_string(i)(1:4) == 'PLOT')     then
      HKL_file%plot     = .true.
     elseIF(arg_string(i)(1:3) == 'NEG')  then
      HKL_file%read_NEG = .false.
     elseif(arg_string(i)(1:4) == 'OUT_') then
      read(arg_string(i)(5:),*, iostat=i_error) num
      if(i_error == 0) then
       file_out_n = num
       file_out = .true.
      end if
     elseif(arg_string(i)(1:4) == 'STAT') then
      hkl_statistics = .true.
     elseif(arg_string(i)(1:7) == 'NO_STAT') then
      hkl_statistics = .false.
     elseif(arg_string(i)(1:4) == 'FREE') then
      hkl_format_free = .true.
     elseif(arg_string(i)(1:5) == 'SHELX') then
      hkl_format_SHELX = .true.
      hkl_format = hkl_SHELX_format
     elseif(arg_string(i)(1:4) == 'FMT=') then
      i1 = index(read_line, 'FMT=')
      if(i1/=0) then
       !read(read_line(i1+4:), '(a)') hkl_format
       read(arg_string(i)(5:), '(a)') hkl_format
       call verif_format(hkl_format)
       if(hkl_format(1:4) == 'FREE') then
        hkl_format_free = .true.
       else
        hkl_format = verif_hkl_format(hkl_format)
       endif
      endif
     !elseif(arg_string(i)(1:6) == 'NO_OUT' .or. arg_string(1:6)(i) == 'NO_DET') then
     ! WRITE_details = .false.
     endif
    end do


   case ('FCF_FILE', 'FILE_FCF', 'READ_FCF')
    FCF_file_name    = ''
    FCF_plot     = .false.
    FCF_plot_stl = .false.

    IF(nb_arg == 0) then
     call check_arg_nb('0', 'FILE_FCF')
     return
    endif
    READ(arg_string(1), *) FCF_file_name
    !IF(nb_arg >1) then
    ! READ(arg_string(2), *) arg1
    ! arg1 = ADJUSTL(arg1)
    !endif

    i = INDEX(FCF_file_name, '.')
    if (i==0) then
     FCF_file_name = TRIM(FCF_file_name)//'.fcf'
    endif
    call test_file_exist(trim(FCF_file_name), file_exist, 'out')
    IF(.NOT. file_exist) then
     FCF_file_name = ''
     return
    endif
    keyword_FILE_FCF = .true.

    do i=2, nb_arg
    if(len_trim(arg_string(i)) >= 4 ) then
     IF(arg_string(i)(1:4) == 'PLOT')     then
      FCF_plot = .true.
      if(len_trim(arg_string(i)) == 8) then
       if(arg_string(i)(1:8) == 'PLOT_STL') FCF_plot_STL = .true.
      end if
     elseif(arg_string(i)(1:4) == 'OUT_') then
      file_out = .true.
      read(arg_string(i)(5:),*) file_out_n
     endif
    end if
    end do

    case ('HKL_DIFF', 'DIFF_HKL')
     IF(nb_arg /= 2) then
      call check_arg_nb('0', 'HKL_DIFF')
      return
     endif
     HKL_file_diff(1) = arg_string(1)
     HKL_file_diff(2) = arg_string(2)
     do i=1, 2
      call test_file_exist(HKL_file_diff(i), file_exist, 'out')
      if(.not. file_exist) return
     end do
     keyword_HKL_diff = .true.

! mode expert ----------------
   case('FIC')   ! 'FILE import.cif
    if(expert_mode) then
     HKL_file%NAME    = ''
     HKL_file%plot     = .false.
     HKL_file%read_NEG = .true.
     HKL_data_known    = .false.
     HKL_file%cif      = .false.
     HKL_file%SHELX    = .false.
     HKL_file%final_y  = .false.
     HKL_file%RAW      = .false.
     HKL_file%M91      = .false.
     HKL_file%M95      = .false.
     HKL_file%INT      = .false.
     HKL_file%COL      = .false.
     file_out          = .false.
     nb_shell = 0
     nb_sort  = 0
     HKL_file%NAME = 'import.cif'
     HKL_file%CIF  = .true.
     call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
     IF(.NOT. file_exist) then
      HKL_file%name = ''
      return
     endif
     keyword_FILE = .true.
    else
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    end if

   case('&FIC')   ! 'FILE import.cif
     HKL_file%NAME    = ''
     HKL_file%plot     = .false.
     HKL_file%read_NEG = .true.
     HKL_data_known    = .false.
     HKL_file%cif      = .false.
     HKL_file%SHELX    = .false.
     HKL_file%final_y  = .false.
     HKL_file%RAW      = .false.
     HKL_file%M91      = .false.
     HKL_file%M95      = .false.
     HKL_file%INT      = .false.
     HKL_file%COL      = .false.
     file_out          = .false.
     nb_shell = 0
     nb_sort  = 0
     HKL_file%NAME = 'import.cif'
     HKL_file%CIF  = .true.
     call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
     IF(.NOT. file_exist) then
      HKL_file%name = ''
      return
     endif
     keyword_FILE = .true.


   case('FST')   ! FILE import_shell_theta.hkl
    if(expert_mode) then
     HKL_file%NAME    = ''
     HKL_file%plot     = .false.
     HKL_file%read_NEG = .true.
     HKL_data_known    = .false.
     HKL_file%cif      = .false.
     HKL_file%SHELX    = .false.
     HKL_file%final_y  = .false.
     HKL_file%RAW      = .false.
     HKL_file%M91      = .false.
     HKL_file%M95      = .false.
     HKL_file%INT      = .false.
     HKL_file%COL      = .false.
     file_out          = .false.
     nb_shell = 0
     nb_sort  = 0
     HKL_file%NAME = 'import_shell_theta.hkl'
     HKL_file%SHELX = .true.
     call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
     IF(.NOT. file_exist) then
      HKL_file%name = ''
      return
     endif
     keyword_FILE = .true.
    else
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    end if

    case('&FST')   ! FILE import_shell_theta.hkl
     HKL_file%NAME    = ''
     HKL_file%plot     = .false.
     HKL_file%read_NEG = .true.
     HKL_data_known    = .false.
     HKL_file%cif      = .false.
     HKL_file%SHELX    = .false.
     HKL_file%final_y  = .false.
     HKL_file%RAW      = .false.
     HKL_file%M91      = .false.
     HKL_file%M95      = .false.
     HKL_file%INT      = .false.
     HKL_file%COL      = .false.
     file_out          = .false.
     nb_shell = 0
     nb_sort  = 0
     HKL_file%NAME = 'import_shell_theta.hkl'
     HKL_file%SHELX = .true.
     call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
     IF(.NOT. file_exist) then
      HKL_file%name = ''
      return
     endif
     keyword_FILE = .true.

   case('FSD')   ! FILE import_shell_d.hkl
    if(expert_mode) then
     HKL_file%NAME    = ''
     HKL_file%plot     = .false.
     HKL_file%read_NEG = .true.
     HKL_data_known    = .false.
     HKL_file%cif      = .false.
     HKL_file%SHELX    = .false.
     HKL_file%final_y  = .false.
     HKL_file%RAW      = .false.
     HKL_file%M91      = .false.
     HKL_file%M95      = .false.
     HKL_file%INT      = .false.
     HKL_file%COL      = .false.
     file_out          = .false.
     nb_shell = 0
     nb_sort  = 0
     HKL_file%NAME = 'import_shell_d.hkl'
     HKL_file%SHELX = .true.
     call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
     IF(.NOT. file_exist) then
      HKL_file%name = ''
      return
     endif
     keyword_FILE = .true.
    else
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    end if

    case('&FSD')   ! FILE import_shell_d.hkl
     HKL_file%NAME    = ''
     HKL_file%plot     = .false.
     HKL_file%read_NEG = .true.
     HKL_data_known    = .false.
     HKL_file%cif      = .false.
     HKL_file%SHELX    = .false.
     HKL_file%final_y  = .false.
     HKL_file%RAW      = .false.
     HKL_file%M91      = .false.
     HKL_file%M95      = .false.
     HKL_file%INT      = .false.
     HKL_file%COL      = .false.
     file_out          = .false.
     nb_shell = 0
     nb_sort  = 0
     HKL_file%NAME = 'import_shell_d.hkl'
     HKL_file%SHELX = .true.
     call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
     IF(.NOT. file_exist) then
      HKL_file%name = ''
      return
     endif
     keyword_FILE = .true.

   case('FIT')   ! FILE import_trans.hkl
    if(expert_mode) then
     HKL_file%NAME    = ''
     HKL_file%plot     = .false.
     HKL_file%read_NEG = .true.
     HKL_data_known    = .false.
     HKL_file%cif      = .false.
     HKL_file%SHELX    = .false.
     HKL_file%final_y  = .false.
     HKL_file%RAW      = .false.
     HKL_file%M91      = .false.
     HKL_file%M95      = .false.
     HKL_file%INT      = .false.
     HKL_file%COL      = .false.
     file_out          = .false.
     nb_shell = 0
     nb_sort  = 0
     HKL_file%NAME = 'import_trans.hkl'
     HKL_file%SHELX = .true.
     call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
     IF(.NOT. file_exist) then
      HKL_file%name = ''
      return
     endif
     keyword_FILE = .true.
    else
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    end if

    case('&FIT')   ! FILE import_trans.hkl
     HKL_file%NAME    = ''
     HKL_file%plot     = .false.
     HKL_file%read_NEG = .true.
     HKL_data_known    = .false.
     HKL_file%cif      = .false.
     HKL_file%SHELX    = .false.
     HKL_file%final_y  = .false.
     HKL_file%RAW      = .false.
     HKL_file%M91      = .false.
     HKL_file%M95      = .false.
     HKL_file%INT      = .false.
     HKL_file%COL      = .false.
     file_out          = .false.
     nb_shell = 0
     nb_sort  = 0
     HKL_file%NAME = 'import_trans.hkl'
     HKL_file%SHELX = .true.
     call test_file_exist(trim(HKL_file%NAME), file_exist, 'out')
     IF(.NOT. file_exist) then
      HKL_file%name = ''
      return
     endif
     keyword_FILE = .true.


!-----------------------------------

   case ('CELL_ESD', 'ESD_CELL')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'CELL_ESD')
     return
    endif
    call Get_cell_param_from_line(arg_line, "ESD_", unit_cell%param_ESD, ok)
    keyword_CELL_ESD = .true.
    known_cell_ESD   = .true.


   case ('P4P', 'READ_P4P')
    keyword_P4P = .true.
    if(nb_arg/=0) then
     i1 = INDEX(arg_string(1), '.')
     if (i1==0) then
      P4P_file_name = TRIM(arg_string(1))//'.P4P'
     else
      READ(arg_string(1), *, iostat = i_error) P4P_file_name
     end if
    else
     P4P_file_name = ""
    endif
    keyword_P4P = .true.

   case ('SEARCH_EXTI', 'FIND_EXTI')
    !n_sig = 0.
    ratio_criteria = 0.03
    HKL_list%ALL = .false.
    if(nb_arg /=0) then
     do i=1, nb_arg
      if(arg_string(i)(1:3) == 'ALL') then
       HKL_list%ALL = .true.
      end if
     end do
    end if

    if(.not. HKL_list%ALL) then
     if(nb_arg/=0) then
      READ(arg_string(1),*, IOSTAT=i_error) var(1)
      IF(i_error ==0) n_sig = var(1)
     end if
     IF(nb_arg >1) then
      READ(arg_string(2),*, IOSTAT=i_error) var(2)
      IF(i_error ==0) ratio_criteria = var(2)
     endif
    end if
    keyword_search_exti = .true.

    case ("SEARCH_SPGR", "SEARCH_SPACE_GROUP", "SEARCH_GROUP",    &
          "CHECK_SPGR",  "CHECK_SPACE_GROUP",  "CHECK_GROUP")

     check_cent            = .true.
     check_all_sg          = .false.
     check_all_sg          = .true.
     ONLY_monoclinic_angle = ONLY_monoclinic_angle_INI
     if(unit_cell%H_M(1:1) == '?') then   ! new 19.06.2015
      check_only_C = .false.              ! new 19.06.2015
      check_cent   = .true.
     else                                 ! new 19.06.2015
      if(unit_cell%H_M(1:1) /= 'P') then  ! new 19.06.2015
       check_only_C = .true.              ! new 19.06.2015
      end if                              ! new 19.06.2015
     end if                               ! new 19.06.2015
     get_SPGR     = .false.
     if (nb_arg /=0) then
      do i=1, nb_arg
       if(len_trim(arg_string(i)) > 2 .and. arg_string(i)(1:3) == 'GET') then
        get_SPGR = .true.
        exit
       end if
      end do

      current_i1 = 0
      current_i2 = 0

      do i=1, nb_arg
       IF(arg_string(i)(1:4) == 'TRIC') then
        crystal_system = "TRIC"
        current_i1 = i
        exit
       elseIF(arg_string(i)(1:4) == 'MONO'  .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'M')) then
        crystal_system = "MONO"
        current_i1 = i
        exit
       elseIF(arg_string(i)(1:5) == 'ORTHO' .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'O')) then
        crystal_system = "ORTHO"
        current_i1 = i
        exit
       elseIF(arg_string(i)(1:5) == 'TETRA' .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'T')) then
        crystal_system = "TETRA"
        current_i1 = i
        exit
       elseIF(arg_string(i)(1:4) == 'TRIG'  .or. arg_string(i)(1:5) == 'RHOMB' &
                                            .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'R')) then
        crystal_system = "TRIG"
        current_i1 = i
        exit
       elseIF(arg_string(i)(1:4) == 'HEXA'  .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'H')) then
        crystal_system = "HEXA"
        current_i1 = i
        exit
       elseIF(arg_string(i)(1:3) == 'CUB'   .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'C')) then
        crystal_system = "CUB"
        current_i1 = i
        exit
       endif
      end do

      do i=1, nb_arg
       if(len_trim(arg_string(i))==3 .and. arg_string(i)(1:3) == 'ALL') then
        check_all_sg = .true.
        check_only_C = .false.      ! new 19.06.2015
        current_i2 = i
        ONLY_monoclinic_angle = .false.
        exit
       elseif(len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'P') then
        check_cent   = .false.
        check_only_C = .false.      ! new 19.06.2015
        current_i2 = i
        exit
       elseif((len_trim(arg_string(i)) == 8 .and. arg_string(i)(1:8) == 'CENTERED') .or. &
              (len_trim(arg_string(i)) == 5 .and. arg_string(i)(1:5) == 'NOT_P')    .or. &
              (len_trim(arg_string(i)) == 4 .and. arg_string(i)(1:4) == 'NO_P') ) then
        check_only_C = .true.
        check_all_sg = .false.
        current_i2 = i
        exit
       elseif((len_trim(arg_string(i)) == 3 .and. arg_string(i)(1:4) == 'OMA') .or. &
              (len_trim(arg_string(i)) == 4 .and. arg_string(i)(1:4) == 'ONLY')) then
        ONLY_monoclinic_angle = .true.
       endif
      end do

      if(current_i1 /=0 .or. current_i2 /=0) then
       if(nb_arg > 1) then
        do i=1, nb_arg
         if(i==current_i1 .or. i==current_i2) cycle
         read(arg_string(i), *, iostat=i_error) var(1)
         if(i_error ==0) then
          n_sig = var(1)
          exit
         endif
        end do
       endif

       if(nb_arg > 2) then
        i_ok = 0
        do i=1, nb_arg
         if(i==current_i1 .or. i==current_i2) cycle
         i_ok = i_ok + 1
         read(arg_string(i), *, iostat=i_error) var(i)
         if(i_error /=0) exit
         if(i_ok == 1) then
          n_sig = var(i)
         elseif(i_ok == 2) then
          threshold = var(i)
          exit
         endif
        end do
       endif

      else ! le systeme n'est pas precise
       if(nb_arg == 1) then
        read(arg_string(1), *, iostat=i_error) var(1)
        if (i_error==0) n_sig = var(1)
       else
        read(arg_string(1), *, iostat=i_error) var(1)
        if (i_error==0) n_sig = var(1)
        read(arg_string(2), *, iostat=i_error) var(2)
        if (i_error==0) threshold = var(2)
       endif
      endif

     end if  ! fin de la condition if(n_arg /=0
     keyword_search_spgr = .true.



   case ('GSG', 'SSG', 'SGG')    ! S_G = SEARCH_GROUP GET
    if (expert_mode) then
     get_SPGR             = .true.
     keyword_search_spgr  = .true.
    else
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    end if

   case ('&GSG', '&SSG', '&SGG')   ! S_G = SEARCH_GROUP GET
    get_SPGR             = .true.
    keyword_search_spgr  = .true.

   case ('SEARCH_MONO', 'SEARCH_MONOCLINIC', 'SEARCH_MONO_ANGLE', 'SEARCH_MONOCLINIC_ANGLE', &
         'GET_MONO',    'GET_MONOCLINIC',    'GET_MONO_ANGLE',    'GET_MONOCLINIC_ANGLE')
    keyword_search_mono  = .true.
    if(nb_arg /=0) then
     read(arg_string(1), *, iostat=i_error) var(1)
     if(i_error == 0) search_mono_criteria(1) = var(1)
     if(nb_arg > 1) then
      read(arg_string(2), *, iostat=i_error) var(2)
      if(i_error == 0) search_mono_criteria(2) = var(2)
     end if
    end if

   case ('SEARCH_TETRA', 'SEARCH_TETRAGONAL', 'GET_TETRA', 'GET_TETRAGONAL')
    keyword_search_tetra  = .true.
    if(nb_arg /=0) then
     read(arg_string(1), *, iostat=i_error) var(1)
     if(i_error == 0) search_tetra_criteria(1) = var(1)
     if(nb_arg > 1) then
      read(arg_string(2), *, iostat=i_error) var(2)
      if(i_error == 0) search_tetra_criteria(2) = var(2)
      if(nb_arg > 2) then
       read(arg_string(3), *, iostat=i_error) var(3)
       if(i_error == 0) search_tetra_criteria(3) = var(3)
       if (search_tetra_criteria(3) > 1.) search_tetra_criteria(3) = search_tetra_criteria(3)/100.
      end if
     end if
    end if

   case ('FIND_HKL', 'SEARCH_HKL')
    IF(nb_arg < 3) then
     call check_arg_nb('3', 'FIND_HKL')
     keyword_FIND_HKL = .false.
     return
    endif
    !READ(arg_string(1), * ) requested_H(1)
    !READ(arg_string(2), * ) requested_H(2)
    !READ(arg_string(3), * ) requested_H(3)

    requested_H_string(1:3) = arg_string(1:3)
    if(requested_H_string(1)(1:1) == 'H'  .OR.  &
       requested_H_string(2)(1:1) == 'K'  .OR.  &
       requested_H_string(3)(1:1) == 'L')  then
     search_H_string = .true.
    else
     READ(arg_string(1), * ) requested_H(1)
     READ(arg_string(2), * ) requested_H(2)
     READ(arg_string(3), * ) requested_H(3)
     search_H_string = .false.
    endif

    search_equiv   = .false.
    search_friedel = .false.
    keyword_FIND_HKL = .true.
    IF(nb_arg == 3) return
    do i=1,nb_arg-3
     IF(arg_string(i+3)(1:5) == 'EQUIV'   ) search_equiv   = .true.
     IF(arg_string(i+3)(1:7) == 'FRIEDEL' ) search_friedel = .true.
    end do

   case ('EQUIV', 'EQUIV_HKL', 'SEARCH_EQUIV', 'SEARCH_EQUIV_HKL', 'FIND_EQUIV', 'FIND_EQUIV_HKL')
    IF(nb_arg < 3) then
     call check_arg_nb('3', 'EQUIV')
     return
    endif
    READ(arg_string(1), * ) requested_H(1)
    READ(arg_string(2), * ) requested_H(2)
    READ(arg_string(3), * ) requested_H(3)
    search_friedel = .false.
    keyword_FIND_HKL_EQUIV = .true.
    IF(nb_arg == 3) return
    do i=1,nb_arg-3
     IF(arg_string(i+3)(1:7) == 'FRIEDEL' ) search_friedel = .true.
    end do


   case ('ABSENT_HKL', 'HKL_ABSENT')
    HKL_list_ABSENT%OUT   = .false.
    HKL_list_ABSENT%ALL   = .false.
    HKL_list_ABSENT%WRITE = .false.
    !n_sig = 0.
    !IF(nb_arg ==0) then
    ! call check_arg_nb('0', 'HKL_ABSENT')
    ! return
    !endif
    IF(nb_arg /=0) then
     do i=1, nb_arg
      IF(arg_string(i)(1:3) == 'ALL')  then
       HKL_list_ABSENT%ALL   = .true.
      elseIF(arg_string(i)(1:3) == 'OUT') then
       HKL_list_ABSENT%OUT   = .true.
      elseIF(arg_string(i)(1:5) == 'WRITE')   then
       HKL_list_ABSENT%WRITE = .true.
      else
       read(arg_string(i),*, iostat=i_error) var(1)
       if(i_error /=0) cycle
       if(i_error == 0) n_sig = var(1)
       !READ(arg_string(i),*) n_sig
       HKL_list_ABSENT%ALL = .false.
      endif
     end do
    END if
    CLOSE(UNIT=31)
    IF(HKL_list_ABSENT%WRITE) then
     i1 = INDEX(HKL_file%NAME,'.')
     WRITE(HKL_file%OUTPUT, '(a,a)') HKL_file%NAME(1:i1-1),'_ABSENT.hkl'
     OPEN(UNIT=31, FILE=TRIM(HKL_file%OUTPUT))
    endif
    keyword_FIND_HKL_ABSENT    = .true.


   CASE ('HKL_NEG', 'HKL_NEGATIVE', 'NEG_HKL', 'NEGATIVE_HKL')
    HKL_list_NEG%OUT   = .false.
    HKL_list_NEG%ALL   = .false.
    HKL_list_NEG%WRITE = .false.
    !n_sig = 0.
    !IF(nb_arg ==0) then
    ! call check_arg_nb('0', 'HKL_NEG')
    ! return
    !endif
    IF(nb_arg /=0) then
     do i=1, nb_arg
      IF(arg_string(i)(1:3) == 'ALL')  then
       HKL_list_NEG%ALL   = .true.
      elseIF(arg_string(i)(1:3) == 'OUT') then
       HKL_list_NEG%OUT   = .true.
      elseIF(arg_string(i)(1:5) == 'WRITE')   then
       HKL_list_NEG%WRITE = .true.
      !else
      ! READ(arg_string(i),*) n_sig
      else
       read(arg_string(i), *, iostat=i_error) var(i)
       if(i_error /= 0) cycle
       if(i_error == 0) n_sig = var(i)
       HKL_list_NEG%ALL   = .false.
      endif
     end do
    endif
    CLOSE(UNIT=31)
    IF(HKL_list_NEG%WRITE) then
     i1 = INDEX(HKL_file%NAME,'.')
     WRITE(HKL_file%OUTPUT, '(a,a)') HKL_file%NAME(1:i1-1),'_NEG.hkl'
     OPEN(UNIT=31, FILE=TRIM(HKL_file%OUTPUT))
    endif
    keyword_FIND_HKL_NEG = .true.


   CASE ('HKL_POS', 'HKL_POSITIVE', 'POS_HKL', 'POSITIVE_HKL')
    HKL_list_POS%OUT   = .false.
    HKL_list_POS%ALL   = .false.
    HKL_list_POS%WRITE = .false.
    !n_sig = 0.
    !IF(nb_arg ==0) then
    ! call check_arg_nb('0', 'HKL_POS')
    ! return
    !endif
    if(nb_arg /=0) then
     do i=1, nb_arg
      IF(arg_string(i)(1:3) == 'ALL')  then
       HKL_list_POS%ALL   = .true.
      elseIF(arg_string(i)(1:3) == 'OUT') then
       HKL_list_POS%OUT   = .true.
      elseIF(arg_string(i)(1:5) == 'WRITE')   then
       HKL_list_POS%WRITE = .true.
      else
       !READ(arg_string(i),*) n_sig
       read(arg_string(i), *, iostat=i_error) var(i)
       if(i_error /=0) cycle
       if(i_error == 0) n_sig = var(i)
       HKL_list_POS%ALL   = .true.
      endif
     end do
    END if
    CLOSE(UNIT=31)
    IF(HKL_list_POS%WRITE) then
     i1 = INDEX(HKL_file%NAME,'.')
     WRITE(HKL_file%OUTPUT, '(a,a)') HKL_file%NAME(1:i1-1),'_POS.hkl'
     OPEN(UNIT=31, FILE=TRIM(HKL_file%OUTPUT))
    endif
    keyword_FIND_HKL_POS = .true.


   case ('STRUCT', 'STRUCTURE_FACTOR', 'FSTR')
    keyword_STRUCT_factor = .true.
    if (nb_arg ==0) then
     call check_arg_nb('0', 'STRUCTURE_FACTOR')
     return
    end if


   case ('RINT', 'R_INT')
    keyword_RINT = .true.
    if (nb_arg /=0) then
     call get_SPG(trim(arg_string(1)))
    end if

   case ("MERGE")
    keyword_MERGE_hkl = .true.
    if (nb_arg /=0) then
     IF(arg_string(1)(1:4) == 'TRIC') then
      space_group_symbol = 'P 1'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:4) == 'MONO')  then
      space_group_symbol = 'P 2'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:5) == 'ORTHO') then
      space_group_symbol = 'P 2 2 2'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:5) == 'TETRA') then
      space_group_symbol = 'P 4'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:4) == 'TRIG')  then
      space_group_symbol = 'P 3'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:4) == 'HEXA')  then
      space_group_symbol = 'P 6'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:3) == 'CUB')   then
      space_group_symbol = 'P 2 3'
      call set_spacegroup(space_group_symbol, SPG)
     endif
    end if

    case ('FRIEDEL', 'FRIEDEL_PAIRS')
     keyword_get_Friedel_pairs_number = .true.
     if (nb_arg /=0) then
      call get_SPG(trim(arg_string(1)))
     end if

   !case ('LIST_HKL', 'LST_HKL', 'LIST_HKL_EXTI', 'LST_HKL_EXTI')
   !CASE ("WRITE_HKL", "WRITE_HKL_LIST", "WRITE_HKL_LST")
    CASE ("FIND_HKL_LIST", "FIND_HKL_LST", "EXTRACT_HKL_LIST", "EXTRACT_HKL_LST")
    if (nb_arg == 0) then
     call check_arg_nb('0', 'FIND_HKL_LIST')
     HKL_list%EXTI_number = 0
     return
    endif

    if(arg_string(1)(1:1)     == 'A') then
     HKL_list%EXTI_number = 13
    elseif(arg_string(1)(1:1) == 'B') then
     HKL_list%EXTI_number = 14
    elseif(arg_string(1)(1:1) == 'C') then
     HKL_list%EXTI_number = 15
    elseif(arg_string(1)(1:1) == 'F') then
     HKL_list%EXTI_number = 16
    elseif(arg_string(1)(1:1) == 'I') then
     HKL_list%EXTI_number = 17
    else
     if(arg_string(1)(1:1) == '#') arg_string(1) = arg_string(1)(2:)
     READ(arg_string(1),*) var(1)
     IF(abs(var(1)) < 1. .or. abs(var(1)) > HKL_rule_nb) then
      call write_info('... Wrong selection rule number ...')
      return
     endif
     HKL_list%EXTI_number = INT(var(1))
    endif

    HKL_list%OUT      = .false.
    HKL_list%ALL      = .false.
    HKL_list%WRITE    = .false.
    HKL_list%SUPPRESS = .false.
    IF(nb_arg > 1) then
     do i=2, nb_arg
      IF(arg_string(i)(1:3) == 'OUT') then
       HKL_list%OUT = .TRUE.
       IF(len_trim(arg_string(i))>3 .and. arg_string(i)(1:4) == 'OUT_') read(arg_string(i)(5:),*) file_out_n
      ELSEIF(arg_string(i)(1:3) == 'ALL') then
       HKL_list%ALL = .TRUE.
      ELSEIF(arg_string(i)(1:5) == 'WRITE') then
       HKL_list%WRITE = .true.
      ELSEIF(arg_string(i)(1:8) == 'SUPPRESS' .or. arg_string(i)(1:6) == 'REMOVE') then
       HKL_list%SUPPRESS = .true.
      else
       READ(arg_string(i),*, iostat=i_error) n_sig
      endif
     end do
    endif
    keyword_FIND_HKL_list = .true.
    IF(HKL_list%WRITE) then
     CLOSE(UNIT=HKL_list_out1_unit)
     i1 = INDEX(HKL_file%NAME,'.')
     IF(i1/=0) then
      if(HKL_list%EXTI_number > 0) then
       IF(HKL_list%EXTI_number < 10) then
        WRITE(HKL_file%OUTPUT, '(a,a,i1,a)') HKL_file%NAME(1:i1-1),'_exti_',HKL_list%EXTI_number,'.hkl'
       else
        WRITE(HKL_file%OUTPUT, '(a,a,i2,a)') HKL_file%NAME(1:i1-1),'_exti_',HKL_list%EXTI_number,'.hkl'
       endif
      else
       IF(abs(HKL_list%EXTI_number) < 10) then
        WRITE(HKL_file%OUTPUT, '(a,a,i1,a)') HKL_file%NAME(1:i1-1),'_exti_m',ABS(HKL_list%EXTI_number),'.hkl'
       else
        WRITE(HKL_file%OUTPUT, '(a,a,i2,a)') HKL_file%NAME(1:i1-1),'_exti_m',ABS(HKL_list%EXTI_number),'.hkl'
       endif
      endif

      OPEN(UNIT=HKL_list_out1_unit, FILE=TRIM(HKL_file%OUTPUT))
     endif
    endif

    IF(HKL_list%SUPPRESS) then
     CLOSE(UNIT=HKL_list_out2_unit)
     i1 = INDEX(HKL_file%NAME,'.')
     IF(i1/=0) then
      if(HKL_list%EXTI_number > 0) then
       IF(HKL_list%EXTI_number < 10) then
        WRITE(HKL_file%OUTPUT2, '(a,a,i1,a)') HKL_file%NAME(1:i1-1),'_suppress_exti_',HKL_list%EXTI_number,'.hkl'
       else
        WRITE(HKL_file%OUTPUT2, '(a,a,i2,a)') HKL_file%NAME(1:i1-1),'_suppress_exti_',HKL_list%EXTI_number,'.hkl'
       endif
      else
       IF(abs(HKL_list%EXTI_number) < 10) then
        WRITE(HKL_file%OUTPUT2, '(a,a,i1,a)') HKL_file%NAME(1:i1-1),'_suppress_exti_',HKL_list%EXTI_number,'.hkl'
       else
        WRITE(HKL_file%OUTPUT2, '(a,a,i2,a)') HKL_file%NAME(1:i1-1),'_suppress_exti_',HKL_list%EXTI_number,'.hkl'
       endif
      endif
      OPEN(UNIT=HKL_list_out2_unit, FILE=TRIM(HKL_file%OUTPUT2))
     endif
    endif


   case ('LIST_EXTI', 'LIST_EXTI_RULE', 'LST_EXTI', 'LST_EXTI_RULE')
    keyword_LIST_EXTI_RULE = .true.

   case ('MATMUL')
    keyword_MATMUL = .true.

   CASE ('MAT', 'MATR', 'MATRIX', 'TRANS', 'TRANSF')
    ! matrice de transformation
    keyword_MAT = .false.
    matrix_num = 0
    IF(nb_arg == 1) then
     IF(arg_string(1)(1:1) == 'I' .OR. arg_string(1)(1:2) == '+I') then
      matrix_num = 1
     ELSEIF(arg_string(1)(1:2) == '-I') then
      matrix_num = 2
     ELSEIF(arg_string(1)(1:1) == '#') then
      READ(arg_string(1)(2:),*,  iostat=i_error) num
      if(i_error == 0) then
       matrix_num = num
      else
       call write_info('')
       WRITE(message_text, '(a)') ' ... Incorrect way to input the matric. Please enter "MAN MAT" keyword for more details.'
       call write_info(trim(message_text))
       call write_info('')
       return
      endif
     else
      call check_arg_nb('9', 'MATRIX')
      return
     endif

    ELSEIF(nb_arg ==3)  then   ! entree de la matrice sous la forme d'un repere abc
     do i=1,3
      call get_mat_from_setting(arg_STRING(i), i, 'A',  Mat)
      call get_mat_from_setting(arg_STRING(i), i, 'B',  Mat)
      call get_mat_from_setting(arg_STRING(i), i, 'C',  Mat)
     end do

    ELSEIF(nb_arg == 2) then
     IF    (arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'R') then
      matrix_num = 8
     elseif(arg_string(1)(1:1) == 'R' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 9
     elseif(arg_string(1)(1:5) == 'R_REV' .AND. arg_string(2)(1:5) == 'R_OBV') then
      matrix_num = 10
     elseif(arg_string(1)(1:5) == 'R_OBV' .AND. arg_string(2)(1:5) == 'R_REV') then
      matrix_num = 11
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'A') then
      matrix_num = 12
     elseif(arg_string(1)(1:1) == 'A' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 13
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'B') then
      matrix_num = 14
     elseif(arg_string(1)(1:1) == 'B' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 15
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'C') then
      matrix_num = 16
     elseif(arg_string(1)(1:1) == 'C' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 17
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'F') then
      matrix_num = 18
     elseif(arg_string(1)(1:1) == 'F' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 19
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'I') then
      matrix_num = 20
     elseif(arg_string(1)(1:1) == 'I' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 21
     else
      matrix_num = 0
      WRITE(message_text, '(a)') '  > Unknown matrix !! '
      call write_info(TRIM(message_text))
      return
     endif

    !elseif(nb_arg /=9) then
    elseif(nb_arg <9) then
     call check_arg_nb('9', 'MATRIX')
     return
    endif

    IF(matrix_num /=0 .and. (matrix_num <1 .OR. matrix_num > max_mat_nb+max_user_mat_nb)) then
     call write_info('')
     WRITE(message_text, '(a,i2,a)') ' ... Numor of matrix has to be in the 1-', max_mat_nb+max_user_mat_nb, ' range...'
     call write_info(trim(message_text))
     call write_info('')
     return
    endif
    IF(matrix_num /=0 ) then
     !if(matrix_num <=max_mat_nb+max_user_mat_nb) then
     ! call def_transformation_matrix()
     !endif
     Mat(1:3,1:3) = transf_mat(1:3,1:3,matrix_num)

    ELSEIF(nb_arg /=3) then
     i1 = index(arg_line, '/')
     if(i1 /=0) then ! elements fractionnaires
      call Get_matrix_coord(arg_line)
     else
      READ(arg_string(1), *) Mat(1,1)
      READ(arg_string(2), *) Mat(1,2)
      READ(arg_string(3), *) Mat(1,3)
      READ(arg_string(4), *) Mat(2,1)
      READ(arg_string(5), *) Mat(2,2)
      READ(arg_string(6), *) Mat(2,3)
      READ(arg_string(7), *) Mat(3,1)
      READ(arg_string(8), *) Mat(3,2)
      READ(arg_string(9), *) Mat(3,3)
     endif
    endif
    keyword_MAT = .true.
    call check_matrice(Mat)
    !WRITE(message_text, '(a,F8.3)') '  > Matrix determinant: ', Mat_det
    !call write_info(TRIM(message_text))
    !call write_info('')

    do i=1, nb_arg
     if(len_trim(arg_string(i)) >= 6 .and. arg_string(i)(1:6) == 'UPDATE') then
      update_parameters = .true.
      exit
     elseif(len_trim(arg_string(i)) >= 9 .and. arg_string(i)(1:9) == 'NO_UPDATE') then
      update_parameters = .false.
      exit
     end if
    end do


   CASE ('USER_MAT', 'USER_MATR', 'USER_MATRIX')
    if(user_mat_nb == 0) then
     call write_info('')
     call write_info(' ... No matrix provided by user in the setting.ini file !')
     call write_info('')
     return
    endif

    ! matrice de transformation
    keyword_MAT = .false.
    matrix_num  = 0
    matrix_text = '?'
    IF(nb_arg /= 0) then
     IF(arg_string(1)(1:1) == '#') then
      READ(arg_string(1)(2:),*) matrix_num
     ELSEIF(arg_string(1)(1:1) == '$') then
      read(arg_string(1)(2:),*) matrix_text
      matrix_text = adjustl(matrix_text)
      do i=1, user_mat_nb
       if(l_case(matrix_text(1:len_trim(matrix_text))) == l_case(user_mat_text(i)(1:len_trim(user_mat_text(i))))) then
        matrix_num = i
        exit
       else
        call write_info('')
        WRITE(message_text, '(a)') ' ... Matrix not defined in the setting file !'
        call write_info(trim(message_text))
        call write_info('')
        return
       endif
      end do

     else
      read(arg_string(1), *, iostat=i_error) matrix_num
      if(i_error/=0) then
       matrix_num = 1
       return
      endif
     endif
    endif

    IF(matrix_num /=0 .and. (matrix_num <1 .OR. matrix_num > user_mat_nb)) then
     call write_info('')
     WRITE(message_text, '(a,i2,a)') ' ... Numor of user matrix has to be in the 1-', user_mat_nb, ' range...'
     call write_info(trim(message_text))
     call write_info('')
     return
    endif
    IF(matrix_num /=0 ) then
     matrix_num = matrix_num + max_mat_nb
     Mat(1:3,1:3) = transf_mat(1:3,1:3,matrix_num)
    endif
    keyword_MAT = .true.
    call check_matrice(Mat)
    WRITE(message_text, '(a,F8.3)') '  > Matrix determinant: ', Mat_det
    call write_info(TRIM(message_text))
    call write_info('')

   CASE ('GET_TRANSF_MAT', 'GET_TRANS_MAT', 'GTM')
    keyword_Get_transf_mat = .true.
    if(nb_arg > 0) then
     if(arg_string(1)(1:1) == "P" .or. arg_string(1)(1:1) == "A" .or. arg_string(1)(1:1) == "B" .or. &
        arg_string(1)(1:1) == "C" .or. arg_string(1)(1:1) == "I" .or. arg_string(1)(1:1) == "R" .or. &
        arg_string(1)(1:1) == "F") reduce_BL = arg_string(1)(1:1)
    end if


   CASE ('REDUCE', 'REDUCE_CELL', 'REDUCED', 'REDUCED_CELL')
    keyword_reduce = .true.
    reduce_BL = 'P'
    if(nb_arg > 0) then
     if(arg_string(1)(1:1) == "P" .or. arg_string(1)(1:1) == "A" .or. arg_string(1)(1:1) == "B" .or. &
        arg_string(1)(1:1) == "C" .or. arg_string(1)(1:1) == "I" .or. arg_string(1)(1:1) == "R" .or. &
        arg_string(1)(1:1) == "F") reduce_BL = arg_string(1)(1:1)
    end if

   CASE ('UB_MAT', 'UB_MATRIX')
    ! UB matrix
    keyword_UB_MAT = .false.
    if(nb_arg <9) then
     call check_arg_nb('9', 'UB_MATRIX')
     return
    endif

    if(nb_arg == 9) then
     UB_mat_1i2i3i = .true.
    elseif(nb_arg == 10) then
     if(arg_string(10) == '112131122232132333') then
      UB_mat_1i2i3i = .true.
     elseif(arg_string(10) == '111213212223313233' .or. arg_string(10)=='P4P') then
      UB_mat_1i2i3i = .false.
     else
      return
     endif
    else
     return
    endif

    if(UB_mat_1i2i3i) then
     READ(arg_string(1), *) UB_matrix(1,1)
     READ(arg_string(2), *) UB_matrix(2,1)
     READ(arg_string(3), *) UB_matrix(3,1)
     READ(arg_string(4), *) UB_matrix(1,2)
     READ(arg_string(5), *) UB_matrix(2,2)
     READ(arg_string(6), *) UB_matrix(3,2)
     READ(arg_string(7), *) UB_matrix(1,3)
     READ(arg_string(8), *) UB_matrix(2,3)
     READ(arg_string(9), *) UB_matrix(3,3)
    else
     READ(arg_string(1), *) UB_matrix(1,1)
     READ(arg_string(2), *) UB_matrix(1,2)
     READ(arg_string(3), *) UB_matrix(1,3)
     READ(arg_string(4), *) UB_matrix(2,1)
     READ(arg_string(5), *) UB_matrix(2,2)
     READ(arg_string(6), *) UB_matrix(2,3)
     READ(arg_string(7), *) UB_matrix(3,1)
     READ(arg_string(8), *) UB_matrix(3,2)
     READ(arg_string(9), *) UB_matrix(3,3)
    endif
    keyword_UB_MAT = .true.



   CASE ('LST_MAT', 'LST_MATR', 'LST_MATRIX',  'LIST_MAT',  'LIST_MATR',  'LIST_MATRIX', 'LIST_TRANSFORMATION_MATRIX', &
         'MAT_LST', 'MAT_LIST', 'MATR_LST', 'MATR_LIST')
    keyword_LST_MAT  = .true.

   CASE ('HKLF5', 'CREATE_HKLF5')
    if(nb_arg <9) then
     call check_arg_nb('9', 'HKLF5')
     return
    endif

    i1 = index(arg_line, '/')
    if(i1 /=0) then    ! elements fractionnaires
     call Get_matrix_coord(arg_line)
    else
     READ(arg_string(1), *) Mat(1,1)
     READ(arg_string(2), *) Mat(1,2)
     READ(arg_string(3), *) Mat(1,3)
     READ(arg_string(4), *) Mat(2,1)
     READ(arg_string(5), *) Mat(2,2)
     READ(arg_string(6), *) Mat(2,3)
     READ(arg_string(7), *) Mat(3,1)
     READ(arg_string(8), *) Mat(3,2)
     READ(arg_string(9), *) Mat(3,3)
    end if
    call check_matrice(Mat)
    WRITE(message_text, '(a,F8.3)') '  > Matrix determinant: ', Mat_det
    call write_info(TRIM(message_text))
    call write_info('')
    keyword_MAT   = .true.
    keyword_HKLF5 = .true.

   case ('FILE_HKLF5', 'READ_HKLF5')
    HKL_file%HKLF5    = .false.
    clusters_out = .false.
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'READ_HKLF5')
     return
    endif
    READ(arg_string(1), *) HKL_file%NAME
    i = INDEX(HKL_file%NAME, '.')
    if (i==0) then
     HKL_file%NAME  = TRIM(HKL_file%NAME)//'.hkl'
     HKL_file%HKLF5 = .true.
    endif
    if(nb_arg ==2) then
     long = len_trim(arg_string(2))
     if(long == 3) then
      if(arg_string(2)(1:3) == 'OUT') Clusters_out = .true.
     end if
    end if
    keyword_read_HKLF5 = .true.

   CASE ('DIAG', 'DIAG_MAT', 'DIAG_MATR', 'DIAG_MATRIX')
    ! matrice 3*3 a diagonaliser
    keyword_DIAG = .false.

    if(nb_arg /=9) then
     call check_arg_nb('9', 'DIAG_MATRIX')
     return
    endif
    i1 = index(arg_line, '/')
    if(i1 /=0) then ! elements fractionnaires
     call Get_matrix_coord(arg_line)
    else
     !READ(arg_string(1:3), *) Mat(1,1:3)  >> warning 402 with IFORT : in call to I/O read routine,
     !READ(arg_string(4:6), *) Mat(2,1:3)  >> an array temporary was created for
     !READ(arg_string(7:9), *) Mat(3,1:3)  >> argument #1
     READ(arg_string(1), *) Mat(1,1)
     READ(arg_string(2), *) Mat(1,2)
     READ(arg_string(3), *) Mat(1,3)
     READ(arg_string(4), *) Mat(2,1)
     READ(arg_string(5), *) Mat(2,2)
     READ(arg_string(6), *) Mat(2,3)
     READ(arg_string(7), *) Mat(3,1)
     READ(arg_string(8), *) Mat(3,2)
     READ(arg_string(9), *) Mat(3,3)

    endif


    keyword_DIAG = .true.
    call check_matrice(Mat)
    WRITE(message_text, '(a,F8.3)') '  > Matrix determinant: ', Mat_det
    call write_info(TRIM(message_text))
    call write_info('')


   CASE ('TRANSLATION', 'TRANSLATE', 'MOVE')
    ! translation
    if (nb_arg <3) then
     call check_arg_nb('3', 'TRANSLATION')
     return
    end if
    if(nb_arg == 3)  then
     READ(arg_string(1:3), *) translat(1:3)
     translat(4) = 1.
    else
     READ(arg_string(1:4), *) translat(1:4)
     if(translat(4) > 0. ) then
      translat(4) = 1.
     else
      translat(4) = -1.
     endif
    endif
    keyword_TRANSL = .true.


   case ('SET')
    IF(nb_arg/=2) then
     call check_arg_nb('2', 'SET')
     return
    endif
    IF(arg_string(1)(1:7) == 'BROWSER') then
     my_browser%name = arg_string(2)
    elseIF(arg_string(1)(1:6) == 'EDITOR')  then
     my_editor%name  = arg_string(2)
    elseIF(arg_string(1)(1:8) == 'PDFLATEX')  then
     my_pdflatex%name  = arg_string(2)
    else
     call check_arg_nb('unknown', 'SET')
     return
    endif
    IF(LEN_TRIM(my_browser%name) == 0 .AND. LEN_TRIM(my_editor%name) == 0 .and. LEN_TRIM(my_pdflatex%name) == 0) return
    keyword_set    = .true.

   case ('REPORT', 'CREATE_REPORT')
    HTML_report  = .true.
    TEXT_report  = .false.
    LATEX_report = .false.
    long_report = .false.
    if (nb_arg ==0) then
     archive_CIF = 'archive.cif'
    ELSEIF(nb_arg==1) then
     IF(arg_string(1)(1:4)=='LONG' .OR. arg_string(1)(1:3)=='EXT') then
      long_report = .TRUE.
      archive_CIF = 'archive.cif'
     elseif(arg_string(1)(1:3) == 'TXT') then
      text_report  = .true.
      long_report = .true.
      archive_CIF = 'archive.cif'
     elseif(arg_string(1)(1:5) == 'LATEX') then
      latex_report = .true.
      long_report  = .true.
      archive_cif  = 'archive.cif'
     endif
    ELSEIF(nb_arg==2) then
     IF(arg_string(1)(1:4)=='LONG' .OR. arg_string(1)(1:3)=='EXT') then
      long_report = .true.
      archive_CIF = arg_string(2)
    ELSEIF(arg_string(1)(1:3) == 'TXT') then
      long_report  = .true.
      text_report  = .true.
      archive_cif  = arg_string(2)
     ELSEIF(arg_string(1)(1:5) == 'LATEX') then
      long_report  = .true.
      latex_report = .true.
      archive_cif = arg_string(2)
     ELSEIF(arg_string(2)(1:4)=='LONG' .OR. arg_string(2)(1:3)=='EXT') then
      long_report = .true.
      archive_CIF = arg_string(1)
     ELSEIF(arg_string(2)(1:3) == 'TXT') then
      long_report  = .true.
      text_report  = .true.
      archive_cif  = arg_string(1)
     ELSEIF(arg_string(2)(1:5) == 'LATEX') then
      long_report  = .true.
      latex_report = .true.
      archive_cif = arg_string(1)
     endif
    end if
    keyword_create_REPORT = .true.

   CASE ('SETTING', 'SETTINGS')
    keyword_setting = .true.


   CASE ('SAVE_SETTING', 'SAVE_SETTINGS')
    keyword_save_setting = .true.

   CASE ('SHELL')
    shell_arg_min_max = .false.
    shell_out         = .false.

    IF(nb_arg ==0) then
     call check_arg_nb('0', 'SHELL')
     return
    endif
    IF(mode_interactif) then
     nb_shell = 1
    else
     nb_shell = nb_shell + 1
    endif
    shell_plot(nb_shell) = .false.

    IF(arg_string(1)(1:1) == 'D')    then
     shell_type(nb_shell) = 'd'
    elseif(arg_string(1)(1:3) == 'STL')   then
     shell_type(nb_shell) = 'stl'
    elseif(arg_string(1)(1:5) == 'THETA') then
     shell_type(nb_shell) = 'theta'
    elseif(arg_string(1)(1:3) == 'INT')   then
     shell_type(nb_shell) = 'int'
    elseif(arg_string(1)(1:4) == 'ISIG')  then
     shell_type(nb_shell) = 'isig'
    else
     call check_arg_nb('unknown', 'SHELL')
     return
    endif

    IF(nb_arg >=3) then
     READ(arg_string(2), *) shell_min(nb_shell)
     READ(arg_string(3), *) shell_max(nb_shell)
     shell_arg_min_max = .true.
    endif

    !IF(nb_arg > 1 .AND. arg_string(2)(1:4) == 'PLOT') shell_plot(nb_shell) = .true.
    if(nb_arg > 1) then
     do i=1, nb_arg
      if(arg_string(i)(1:4) == 'PLOT') shell_plot(nb_shell) = .true.
      if(arg_string(i)(1:3) == 'OUT')  shell_out = .true.
     end do
    end if

!-----------------------------------------------------
   CASE ('SD7')   ! mode expert : SD7 = SHELL d 0.77 7.
    if(expert_mode) then
     nb_shell = 1
     shell_type(nb_shell) = 'd'
     shell_min(1) = 0.77
     shell_max(1) = 7.
     shell_arg_min_max = .true.
    else
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    endif
   CASE ('&SD7')   ! mode expert : SD7 = SHELL d 0.77 7.
    nb_shell = 1
    shell_type(nb_shell) = 'd'
    shell_min(1) = 0.77
    shell_max(1) = 7.
    shell_arg_min_max = .true.

   CASE ('ST25')   ! mode expert : ST25 : SHELL THETA 2.5 25
    if(expert_mode) then
     nb_shell = 1
     shell_type(nb_shell) = 'theta'
     shell_min(1) = 0.25
     shell_max(1) = 25.
     shell_arg_min_max = .true.
    else
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    endif

    CASE ('&ST25')   ! mode expert : ST25 : SHELL THETA 2.5 25
     nb_shell = 1
     shell_type(nb_shell) = 'theta'
     shell_min(1) = 0.25
     shell_max(1) = 25.
     shell_arg_min_max = .true.
 !------------------------------------------------------------------


   CASE ('SORT')
    IF(nb_arg ==0) then
     call check_arg_nb('0', 'SORT')
     return
    endif
    IF(mode_interactif) then
     nb_sort = 1
    else
     nb_sort = nb_sort + 1
    endif

    !IF(nb_arg > 1 .AND. arg_string(2)(1:4) == 'PLOT') sort_plot(nb_sort) = .true.
    sort_plot = .false.
    sort_out  = .false.
    IF(nb_arg > 1) then
     do i=2, nb_arg
      IF(arg_string(i)(1:4) == 'PLOT') then
       sort_plot(nb_sort) = .TRUE.
      elseIF(arg_string(i)(1:3) == 'OUT') then
       sort_out(nb_sort) = .true.
       IF(arg_string(i)(4:4) == '_')  READ(arg_string(i)(5:),*) sort_out_n(nb_sort)
      endif
     END do
    endif

    IF(arg_string(1)(1:1) == 'D')     then
     sort_type(nb_sort) = 'd'
    ELSEIF(arg_string(1)(1:3) == 'STL')   then
     sort_type(nb_sort) = 'stl'
    ELSEIF(arg_string(1)(1:5) == 'THETA') then
     sort_type(nb_sort) = 'theta'
    ELSEIF(arg_string(1)(1:3) == 'INT')   then
     sort_type(nb_sort) = 'int'
    ELSEIF(arg_string(1)(1:4) == 'ISIG')  then
     sort_type(nb_sort) = 'isig'
    ELSE
     call check_arg_nb('unknown', 'SORT')
     return
    ENDIF


   CASE ('GEN_HKL', 'GENERATE_HKL', 'GENERATE_HKL_LIST')
    i1 = INDEX(read_line, '=')
    if (i1==0) then
     call check_arg_nb('', 'GEN_HKL')
     return
    endif

    if (arg_string(1)(1:3) == '2TH') then
     HKL_2THETA = .true.
    ELSEIF(arg_string(1)(1:2) == 'TH') then
     HKL_THETA = .true.
    ELSEIF(arg_string(1)(1:3) == 'STL') then
     HKl_STL = .true.
    ELSEIF(arg_string(1)(1:1) == 'D') then
     HKL_D = .true.
    ELSEIF(arg_string(1)(1:1) == 'Q') then
     HKL_Q = .true.
    end if

    READ(read_line(i1+1:), *) X_min

    i2 = INDEX(read_line, 'MAX=', back=.TRUE.)
    IF(i2==0) return
    new_line = read_line(i2+4:)
    call nombre_de_colonnes(new_line, nb_arg)
    !if(nb_arg == 1) then
    ! read(new_line, *) X_max
    ! arg_line = ''
    !else
    ! read(new_line, *) X_max, arg_line
    !endif

    !!  - fev. 2011 !!
    write_HKL       = .false.
    on_screen_PRF   = .false.
    create_PAT      = .false.
    size_broadening = .false.
    INI_create%PRF  = .false.
    read(new_line, *) arg_string(1:nb_arg)
    read(arg_string(1), *) X_max
    if(nb_arg /=1) then
     do i=2, nb_arg
      if(arg_string(i)(1:3) == 'OUT') write_HKL  = .true.
      if(arg_string(i)(1:7) == 'PATTERN' .or. arg_string(i)(1:3) == 'PAT'   .or.    &
         arg_string(i)(1:7) == 'PROFILE' .or. arg_string(i)(1:4) == 'PROF'  .or.    &
         arg_string(i)(1:3) == 'PRF')       create_PAT = .true.
      if(arg_string(i)(1:5) == 'SIZE=') then
       read(arg_string(i)(6:), *) particle_size
       if (create_pat) size_broadening = .true.
      endif
      if(arg_string(i)(1:4) == 'PM2K') pm2k_out = .true.
     end do
    end if
    !if(create_PAT) write_HKL=.true.
    !!-----------------!!


    !READ(read_line(i2+1:), *) X_max, arg_line
    !if(len_trim(arg_line) /=0 ) then
    ! write_HKL = .false.
    ! arg_line = adjustl(arg_line)
    ! arg_line = u_case(arg_line)
    ! if(arg_line(1:3) =='OUT') write_HKL = .true.
    !else
    ! WRITE_HKL  = .false.
    !endif

    keyword_GENHKL = .true.

    !call write_info('')
    !call write_info('  > GEN_HKL ')
    !call write_info('')

!--- expert mode ---------------
   CASE ('PDP', 'PATTERN', 'PAT', 'PDP_CU', 'PATTERN_Cu')
    if(expert_mode) then
     HKL_2THETA     = .true.
     !X_min          = 0.
     !X_max          = 120.
     X_min          = X_pattern%xmin
     X_max          = X_pattern%xmax
     create_PAT     = .true.
     write_HKL      = .true.
     keyword_GENHKL = .true.
    end if

   CASE ('&PDP', '&PATTERN', '&PAT', '&PDP_CU', '&PATTERN_Cu')
     HKL_2THETA     = .true.
     !X_min          = 0.
     !X_max          = 120.
     X_min          = X_pattern%xmin
     X_max          = X_pattern%xmax
     create_PAT     = .true.
     write_HKL      = .true.
     keyword_GENHKL = .true.
!------------------------------

   CASE ('GEN_SAT')
    keyword_GENSAT = .true.

   CASE ('INSIDE')
    keyword_INSIDE = .true.

   CASE ('DIFF', 'DIFFERENCE', 'DIFF_CALC', 'DIFFERENCE_CALCULATION')
    IF(nb_arg < 2) then
     call check_arg_nb('2', 'DIST')
     keyword_DIFF = .false.
     return
    endif
    if(step_by_step) then
     nb_diff_calc = 1
    else
     nb_diff_calc = nb_diff_calc + 1
    end if

    READ(arg_string(1:2), *) atom1_diff(nb_diff_calc), atom2_diff(nb_diff_calc)
    atom1_diff(nb_diff_calc) = ADJUSTL(atom1_diff(nb_diff_calc))
    atom2_diff(nb_diff_calc) = ADJUSTL(atom2_diff(nb_diff_calc))
    keyword_DIFF = .true.

   CASE ('DIST', 'DISTANCE', 'ATOMIC_DISTANCE')
    keyword_DIST_X    = .false.
    keyword_DIST_plus = .false.
    keyword_DHA       = .false.
    IF(nb_arg < 2) then
     call check_arg_nb('2', 'DIST')
     return
    endif
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  CELL keyword is mandatory for interatomic distance calculations')
     call write_info('')
     return
    endif
    if(step_by_step) then
     nb_dist_calc = 1
    else
     nb_dist_calc = nb_dist_calc + 1
    end if

    READ(arg_string(1:2), *) atom1_dist(nb_dist_calc), atom2_dist(nb_dist_calc)
    atom1_dist(nb_dist_calc) = ADJUSTL(atom1_dist(nb_dist_calc))
    atom2_dist(nb_dist_calc) = ADJUSTL(atom2_dist(nb_dist_calc))

   CASE ('CONN', 'CONNECT', 'CONNECTIVITY')
    !if(nb_arg < 1) then
    ! call check_arg_nb('1', 'CONN')
    ! return
    !end if
    CONN_all           = .false.
    CONN_all_X         = .false.
    CONN_self          = .false.
    CONN_ang           = .false.
    CONN_out_condensed = .false.
    CONN_excluded      = .false.
    calcul_BVS         = .false.
    create_SHAPE_file  = .false.
    poly_vol_calc      = .false.
    CIF_dist%n_text = 0

    if(nb_arg /=0) then
     string_ok = .false.
     call test_arg_string('BVS',arg_string, nb_arg, string_ok)
     if(string_ok) calcul_BVS = .true.
    end if

    if(nb_arg /=0) then
     string_ok = .false.
     call test_arg_string('SHAPE',arg_string, nb_arg, string_ok)
     if(string_ok) create_SHAPE_file = .true.
    end if

    if(nb_arg /=0) then
     string_ok = .false.
     call test_arg_string('VOL',arg_string, nb_arg, string_ok)
     if(string_ok) poly_vol_calc = .true.
    end if

    if(nb_arg /=0) then
     string_ok = .false.
     call test_arg_string('ANG',arg_string, nb_arg, string_ok)
     if(string_ok) CONN_ang = .true.
    end if

    if(nb_arg /=0) then
     string_ok = .false.
     call test_arg_string('CONDENSED',arg_string, nb_arg, string_ok)
     if(string_ok) then
      CONN_out_condensed = .true.
     else
      call test_arg_string('SHORT',arg_string, nb_arg, string_ok)
      if(string_ok) CONN_out_condensed = .true.
     end if
    end if


    if(nb_arg /=0) then
     string_ok = .false.
     call test_arg_string('SELF',arg_string, nb_arg, string_ok)
     if(.not. string_OK) call test_arg_string('AUTO',arg_string, nb_arg, string_ok)
     if(string_ok) then
      CONN_self          = .true.
      CONN_all           = .false.
      CONN_all_X         = .false.
      CONN_ONLY_X        = .false.
      CONN_ang           = .false.
      !CONN_out_condensed = .false.
      CONN_excluded      = .false.
      calcul_BVS         = .false.
      create_SHAPE_file  = .false.
      poly_vol_calc      = .false.
     end if
    end if


    ! test if CONN_dmin= is present
    i_CONN_dmin = 0
    if(nb_arg /=0) then
     do i=1, nb_arg
      if(len_trim(arg_string(i)) > 4 .and. arg_string(i)(1:4) == 'MIN=') then
       read(arg_string(i)(5:), *, iostat=i_error) CONN_dmin
       if(i_error /=0)      CONN_dmin = CONN_dmin_ini
       if(CONN_dmin < 0.05) CONN_dmin = CONN_dmin_ini
       i_CONN_dmin = i
      end if
     end do
     if(i_CONN_dmin /=0) then
      if(i_CONN_dmin /= nb_arg) then
       do i=i_CONN_dmin, nb_arg
        arg_string(i) = arg_string(i+1)
       end do
      end if
      nb_arg = nb_arg - 1
     end if
    endif

    ! test if CONN_dmax= is present
    i_CONN_dmax = 0
    if(nb_arg /=0) then
     do i=1, nb_arg
      if(len_trim(arg_string(i)) > 4 .and.arg_string(i)(1:4) == 'MAX=') then
       read(arg_string(i)(5:), *, iostat=i_error) CONN_dmax
       if(i_error /=0) CONN_dmax = CONN_dmax_ini
       if(CONN_dmax > 25.) CONN_dmax = CONN_dmax_ini
       i_CONN_dmax = i
      end if
     end do
      if(i_CONN_dmax /=0) then
      if(i_CONN_dmax /= nb_arg) then
       do i=i_CONN_dmax, nb_arg
        arg_string(i) = arg_string(i+1)
       end do
      end if
      nb_arg = nb_arg - 1
     end if
    endif

    ! test if NO_ is present
    if(nb_arg /=0) then
     do i=1, nb_arg
      if(len_trim(arg_string(i)) > 3 .and. arg_string(i)(1:3) == 'NO_')then
       CONN_excluded = .true.
       read(arg_string(i)(4:), *) ATOM_CONN%excluded
      end if
     end do
    end if

    if(nb_arg == 0) CONN_all = .true.

    if (nb_arg /=0) then
     if(len_trim(arg_string(1)) == 3 .and. arg_string(1)(1:3) == 'ALL') then
      CONN_all = .true.
     elseif(len_trim(arg_string(1)) > 4 .and. arg_string(1)(1:4) == 'ALL_') then
      CONN_all_X = .true.
      read(arg_string(1)(5:), *) CONN_species
     elseif(len_trim(arg_string(1)) > 5 .and. arg_string(1)(1:5) == 'ONLY_') then
      CONN_ONLY_X = .true.
      read(arg_string(1)(6:),*) CONN_species
     else
      read(arg_string(1), *) atom_CONN%label
      atom_CONN%label = adjustl(atom_CONN%label)
     endif
    endif


    if(nb_arg > 1 .and. .not. CONN_excluded) then
     read(arg_string(nb_arg), *, iostat=i_error) CONN_dmax
     if(i_error /=0) CONN_dmax = CONN_dmax_ini
     if(CONN_dmax < 0.1 .or. CONN_dmax > 25.) CONN_dmax = CONN_dmax_ini
    endif
    if(CONN_dmin > CONN_dmax) then
     CONN_dmin = CONN_dmin_ini
     CONN_dmax = CONN_dmax_ini
    end if
    keyword_CONN = .true.



   CASE ('DIST_', 'DIST_X', 'DIST_MULT')
    keyword_DIST_X    = .false.
    keyword_DIST_plus = .false.
    keyword_DHA       = .false.
    IF(nb_arg < 3) then
     call check_arg_nb('3', 'DIST_X')
     return
    endif
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  CELL keyword is mandatory for interatomic distance calculations')
     call write_info('')
     return
    endif
    !nb_dist_calc = nb_dist_calc + 1
    nb_dist_calc = 1
    call get_arg_ABx(arg_string, DIST_coef, atom1_dist(nb_dist_calc), atom2_dist(nb_dist_calc), lecture_ok)
    !READ(arg_string(1)  , *) DIST_coef
    !READ(arg_string(2:3), *) atom1_dist(nb_dist_calc), atom2_dist(nb_dist_calc)
    atom1_dist(nb_dist_calc) = ADJUSTL(atom1_dist(nb_dist_calc))
    atom2_dist(nb_dist_calc) = ADJUSTL(atom2_dist(nb_dist_calc))
    keyword_DIST_X = .true.

   CASE ('DIST_PLUS', 'DIST_+')
    keyword_DIST_X    = .false.
    keyword_DIST_plus = .false.
    keyword_DHA       = .false.
    IF(nb_arg < 3) then
     call check_arg_nb('3', 'DIST_+')
     return
    endif
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  CELL keyword is mandatory for interatomic distance calculations')
     call write_info('')
     return
    endif
    !nb_dist_calc = nb_dist_calc + 1
    nb_dist_calc = 1
    lecture_ok = .false.
    call get_arg_ABx(arg_string, DIST_plus, atom1_dist(nb_dist_calc), atom2_dist(nb_dist_calc), lecture_ok)
    if(.not. lecture_ok) then
     call write_info('')
     call write_info(' ! Wrong arguments for DIST_plus keyword !')
     call write_info('')
     return
    end if

    !READ(arg_string(1)  , *) DIST_plus
    !READ(arg_string(2:3), *) atom1_dist(nb_dist_calc), atom2_dist(nb_dist_calc)
    atom1_dist(nb_dist_calc) = ADJUSTL(atom1_dist(nb_dist_calc))
    atom2_dist(nb_dist_calc) = ADJUSTL(atom2_dist(nb_dist_calc))
    keyword_DIST_plus = .true.

   CASE ('DIST_DHA', 'DHA', 'POS_H', 'CALC_POS_H')
    keyword_DIST_X    = .false.
    keyword_DIST_plus = .false.
    keyword_DHA       = .false.
    if(nb_arg < 2) then
     call check_arg_nb('2', 'DHA')
     return
    end if
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  CELL keyword is mandatory for interatomic distance calculations')
     call write_info('')
     return
    endif
    nb_dist_calc = 1
    READ(arg_string(1:2), *) atom1_dist(nb_dist_calc), atom2_dist(nb_dist_calc)
    atom1_dist(nb_dist_calc) = ADJUSTL(atom1_dist(nb_dist_calc))
    atom2_dist(nb_dist_calc) = ADJUSTL(atom2_dist(nb_dist_calc))
    if(nb_arg > 3)  read(arg_string(3), *) DIST_AH
    keyword_DHA = .true.

   CASE ('ANG', 'ANGLE')
    IF(nb_arg < 3) then
     call check_arg_nb('', 'ANG')
     return
    endif
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  CELL keyword is mandatory for interatomic angle calculation')
     call write_info('')
     return
    endif

    !nb_ang_calc = nb_ang_calc + 1
    nb_ang_calc = 1

    if (nb_arg == 3) then
     READ(arg_string(1:3),*) atom1_ang(nb_ang_calc),  atom2_ang(nb_ang_calc), atom3_ang(nb_ang_calc)
     atom1_ang(nb_ang_calc) = ADJUSTL(atom1_ang(nb_ang_calc))
     atom2_ang(nb_ang_calc) = ADJUSTL(atom2_ang(nb_ang_calc))
     atom3_ang(nb_ang_calc) = ADJUSTL(atom3_ang(nb_ang_calc))
     atom4_ang(nb_ang_calc) = 'ZZZ'

    ELSEIF (nb_arg == 4) then
     READ(arg_string(1:4),*) atom1_ang(nb_ang_calc),  atom2_ang(nb_ang_calc), atom3_ang(nb_ang_calc), atom4_ang(nb_ang_calc)
     atom1_ang(nb_ang_calc) = ADJUSTL(atom1_ang(nb_ang_calc))
     atom2_ang(nb_ang_calc) = ADJUSTL(atom2_ang(nb_ang_calc))
     atom3_ang(nb_ang_calc) = ADJUSTL(atom3_ang(nb_ang_calc))
     atom4_ang(nb_ang_calc) = ADJUSTL(atom4_ang(nb_ang_calc))

    else
     call write_info('')
     call write_info('    ... uncorrect list of atoms ...')
     call write_info('')
     nb_ang_calc = nb_ang_calc - 1
    end if


   case ('DIR_ANG', 'DIRANG', 'DIRECT_ANGLE')
    IF(nb_arg <6 ) then
     call check_arg_nb('6', 'DIR_ANG')
    endif
    READ(arg_string(1:6), *, IOSTAT=i_error) var(1:6)
    IF(i_error /=0) then
     call check_arg_nb('uncorrect', 'DIR_ANG')
    endif
    IF(mode_interactif) then
     nb_da = 1
    else
     nb_da = nb_da + 1
    endif
    U1_da(1:3, nb_da) = var(1:3)
    U2_da(1:3, nb_da) = var(4:6)
    keyword_DIR_ANG = .true.



   case ('REC_ANG', 'RECANG', 'RECIPROCAL_ANGLE')
    IF(nb_arg <6 ) then
     call check_arg_nb('6', 'REC_ANG')
    endif
    READ(arg_string(1:6), *, IOSTAT=i_error) var(1:6)
    IF(i_error /=0) then
     call check_arg_nb('uncorrect', 'REC_ANG')
    endif
    IF(mode_interactif) then
     nb_ra = 1
    else
     nb_ra = nb_ra + 1
    endif
    U1_ra(1:3, nb_ra) = var(1:3)
    U2_ra(1:3, nb_ra) = var(4:6)
    keyword_REC_ANG = .true.

   CASE ('PLANE', 'PLAN')
    !IF(nb_arg < 3 ) then
    ! call check_arg_nb('3', 'PLANE')
    ! return
    !endif
    IF(mode_interactif) then
     nb_plane = 1
    else
     nb_plane = nb_plane + 1
    endif
    if(nb_arg > 20) nb_arg = 20   ! 20 atomes max. pour definir le plan moyen
    atom_plane_nb(nb_plane) = nb_arg
    do i=1, nb_arg
     read(arg_string(i), *) atom_plane(i, nb_plane)
    end do
    !read(arg_string(1:nb_arg), *) atom_plane(1:nb_arg, nb_plane)
    !READ(arg_string(1:3),*) atom1_plane(nb_plane),  atom2_plane(nb_plane), atom3_plane(nb_plane)
    !atom1_plane(nb_plane) = ADJUSTL(atom1_plane(nb_plane))
    !atom2_plane(nb_plane) = ADJUSTL(atom2_plane(nb_plane))
    !atom3_plane(nb_plane) = ADJUSTL(atom3_plane(nb_plane))

    keyword_PLANE = .true.


   CASE ('BARY', 'CENTROID')
    keyword_BARY = .false.
    if(nb_arg == 1 .and. arg_string(1) == 'ALL') then
     if(nb_atom /=0) then
      nb_arg = nb_atom
      bary_all = .true.
     else
      call write_info('! No input atoms !')
      return
     endif
    else
     IF(nb_arg < 2) then
      call check_arg_nb('at least 2', 'BARY')
      return
     endif

     if(nb_arg >= 3 .and. arg_string(2)(1:1) == ">") then
      nb_sup = nb_arg - 3
      do i=1, nb_sup
       sup_string(i) = arg_string(3+i)
      end do
      nb_arg_sep = 0
      ok = .false.
      call get_bary_atoms(arg_string, nb_arg_sep, ok)
      if(.not. ok) return
      do i=1, nb_sup
       arg_string(nb_arg_sep+i) = sup_string(i)
      end do
      nb_arg = nb_arg_sep + nb_sup
     end if
    end if
    nb_bary_calc = nb_bary_calc + 1
    keyword_BARY = .true.

    nb_atom_bary(nb_bary_calc) = nb_arg
    if(bary_all) then
     atom_bary(nb_bary_calc, 1:nb_arg) = atom_label(1:nb_arg)
    else
     !read(arg_string(1:nb_arg), *) atom_bary(nb_bary_calc,1:nb_arg)  ! <<< warning with ifc
     do i=1, nb_arg
      read(arg_string(i), *) atom_bary(nb_bary_calc,i)
     end do
    endif

   case ('TOLMAN_CONE_ANGLE', 'TOLMAN_CONE', 'TOLMAN_ANGLE', 'CONE_ANGLE', 'CONE', 'TOLMAN')
    keyword_TOLMAN = .false.
    modif_rvdw     = .false.
    if(nb_arg <5) then
     call check_arg_nb('5', 'TOLMAN_ANGLE')
     return
    end if
    !nb_tolman_calc = nb_tolman_calc + 1
    nb_tolman_calc = 1
    do i=1, 5
     read(arg_string(i), *) atom_tolman(nb_tolman_calc,i)
    end do
    if(nb_arg >5) then
     if(nb_arg > 8) nb_arg = 8
     do i=6, nb_arg
      read(arg_string(i), *) r_vdw(nb_tolman_calc, i-5)
     end do
     modif_rvdw = .true.
    end if
    keyword_TOLMAN = .true.


   case ('STL', 'STL_HKL', 'SINTHETA/WAVE', 'SINTHETA/LAMBDA')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'STL')
     return
    endif
    nb_STL_value = nb_arg
    READ(arg_string(1:nb_arg), *) STL_value(1:nb_arg)
    keyword_STL = .true.

   case ('D_HKL', 'DHKL')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'D_HKL')
     return
    endif
    nb_dhkl_value = nb_arg
    READ(read_line(long_kw+1:), *) (dhkl_value(i), i=1, nb_dhkl_value)
    keyword_DHKL = .true.

   case ('D_STAR', 'D_STAR_HKL', 'DSTAR', 'DSTARHKL', 'DSTAR_HKL')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'D_STAR')
     return
    endif
    nb_dstar_value = nb_arg
    READ(arg_string(1:nb_arg), *) dstar_value(1:nb_arg)
    keyword_DSTAR = .true.

   case ('Q_HKL', 'QHKL')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'Q_HKL')
     return
    endif
    nb_Qhkl_value = nb_arg
    READ(arg_string(1:nb_arg), *) Qhkl_value(1:nb_arg)
    keyword_QHKL = .true.

   case ('THETA', 'TH', 'TH_HKL', 'THETAHKL', 'THETA_HKL')
    IF(nb_arg== 0) then
     call check_arg_nb('0', 'THETA_HKL')
     return
    endif
    nb_theta_value = nb_arg
    READ(arg_string(1:nb_arg), *) theta_value(1:nb_arg)
    keyword_theta = .true.

   case ('2THETA', '2TH', '2TH_HKL', '2THETA_HKL', 'TWO_THETA', 'TWO_THETA_HKL')
    IF(nb_arg == 0) then
     call check_arg_nb('0', '2THETA_HKL')
     return
    endif
    nb_2theta_value = nb_arg
    READ(arg_string(1:nb_arg), *) TWO_theta_value(1:nb_arg)
    keyword_2theta = .true.


   case ('RESET', 'RAZ', 'INITIALIZE', 'INIT')
    keyword_RESET = .true.

   CASE ('SYST', 'CMD', 'COMMAND', 'DOS', 'DOS_COMMAND')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'SYST')
     return
    endif
    READ(arg_line,'(a)') SYST_command
    SYST_command = ADJUSTL(SYST_command)
    IF(LEN_TRIM(SYST_command) /=0) keyword_SYST =.true.

   case ('DIR')
    IF(nb_arg == 0) then
     SYST_command = 'dir'
    else
     WRITE(SYST_command, '(2a)') 'dir ', TRIM(arg_line)
    endif
    keyword_SYST = .true.



   CASE ('THERM', 'THERMAL', 'ADP', 'THERM_SHELX', 'THERMAL_SHELX', 'ADP_SHELX', 'THERM_SHELXL', 'THERMAL_SHELXL', 'ADP_SHELXL')
    THERM%Uiso              = .false.
    THERM%Biso              = .false.
    THERM%Uij               = .false.
    THERM%Bij               = .false.
    THERM%Beta              = .false.
    THERM%aniso             = .false.
    keyword_THERM           = .false.

    if(input_keyword(1:11) == 'THERM_SHELX'     .or.   &
       input_keyword(1:13) == 'THERMAL_SHELX'   .or.   &
       input_keyword(1:9)  == 'ADP_SHELX')  then
     keyword_THERM_SHELX = .true.
    endif
    IF(nb_arg <2) then
     call check_arg_nb('at least 2', 'THERM')
     return
    endif

    IF(arg_string(1)(1:4)     == 'BISO'    .or. arg_string(1)(1:5) == 'B_ISO') then
     THERM%Biso = .true.
    ELSEIF(arg_string(1)(1:4) == 'UISO'    .or. arg_string(1)(1:5) == 'U_ISO') then
     THERM%Uiso = .true.
    ELSEIF(arg_string(1)(1:4) == 'U_IJ'    .OR. arg_string(1)(1:3) == 'UIJ') then
     THERM%Uij   = .true.
     THERM%aniso = .true.
    ELSEIF(arg_string(1)(1:4) == 'B_IJ'    .OR. arg_string(1)(1:3) == 'BIJ') then
     THERM%Bij   = .true.
     THERM%aniso = .true.
    ELSEIF(arg_string(1)(1:7) == 'BETA_IJ' .OR. arg_string(1)(1:6) == 'BETAIJ' .or. arg_string(1)(1:4) == 'BETA') then
     THERM%Beta  = .true.
     THERM%aniso = .true.
    endif
    THERM%nb_values = nb_arg - 1
    IF(THERM%aniso) then
     IF(.NOT. keyword_CELL) then
      call write_info('')
      call write_info('  Cell parameters have to be known for anisotropic ADP calculations')
      call write_info('')
      return
     endif
     if (THERM%nb_values <6) then
      call check_arg_nb('6', 'THERM')
      return
     endif
    !else
    ! THERM%nb_values = 6
    endif
    READ(arg_string(2:nb_arg),*) therm%values(1:THERM%nb_values)

    keyword_THERM = .true.

   CASE ('TRANSMISSION')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'TRANSMISSION')
     return
    endif
    keyword_transmission = .true.
    nb_dim_transmission = nb_arg
    READ(arg_string(1:nb_arg),*)  dim_transmission(1:nb_arg)


   case ('SG_ALL',  'SP_ALL')
    write_SPG_info      = .true.
    write_SPG_info_all  = .true.
    write_SPG_exti      = .true.
    write_SPG_subgroups = .true.

   CASE ('SG_INFO', 'SP_INFO', 'SPACE_GROUP_INFO', 'LIST_SPACE_GROUP_INFO')
    WRITE_SPG_info      = .true.
    write_SPG_info_all  = .false.
    write_SPG_exti      = .false.
    write_SPG_subgroups = .false.
    IF(nb_arg /=0) then
     do i=1, nb_arg
      long = len_trim(arg_string(i))
      if(long >=3 .and. arg_string(i)(1:3)=='ALL') then
       write_SPG_info_all = .true.
       write_SPG_exti     = .true.
      !elseif(long==4 .and. arg_string(i)(1:4) == 'EXTI') then
      ! write_SPG_exti     = .true.
      end if
     end do
    end if

   CASE ('SG_EXTI', 'SP_EXTI', 'SG_EXTINCTIONS', 'SPACE_GROUP_EXTI', 'SPACE_GROUP_EXTINCTIONS')
    WRITE_SPG_exti = .true.
    WRITE_SPG_info = .false.

   CASE ('SG_SUB', 'SG_SUBGROUP')
    write_SPG_subgroups = .true.
    write_SPG_info      = .false.
    write_SPG_info_all  = .false.


   !CASE ('SYM_OP', 'SYMM_OP', 'SYM_OPERATOR', 'SYMM_OPERATOR')
   CASE ("WRITE_SYM_OP", "WRITE_SYMM_OP", "WRITE_SYM", "WRITE_SYMM", "WRITE_SYMMETRY_OPERATORS", "WRITE_SYMOP")
    WRITE_symm_op             = .true.
    WRITE_symm_op_condensed   = .false.
    WRITE_SYMM_op_on_one_line = .false.
    WRITE_SHELX_symm_op       = .false.
    if (nb_arg /=0) then
     do i=1, nb_arg
      if(arg_string(i)(1:5) == 'SHELX') then
       WRITE_SHELX_symm_op  = .true.
       exit
      elseif(arg_string(i)(1:8) == 'ONE_LINE') then
       WRITE_symm_op_on_one_line = .true.
       exit
      elseif(arg_string(i)(1:9) == 'CONDENSED') then
       WRITE_symm_op_condensed = .true.
       exit
      end if
     end do
    endif


   CASE ('APPLY_OP', 'APPLY_SYMMETRY_OPERATOR', 'APPLY_SYM_OP', 'APPLY_SYMOP')
    WRITE_APPLY_symm = .true.


   CASE ('SITE_INFO', 'LIST_SITE_INFO', 'GEN_EQUIV_ATOMS', 'WYCKOFF')
    WRITE_site_info         = .true.
    WRITE_PCR_site_info     = .false.
    WRITE_PCR_mag_site_info = .false.

    if(nb_arg /=0) then
     do i=1, nb_arg
      if(len_trim(arg_string(i)) == 3 .and. arg_string(i)(1:3) == 'PCR')     WRITE_PCR_site_info     = .true.
      if(len_trim(arg_string(i)) == 3 .and. arg_string(i)(1:3) == 'MAG')     WRITE_PCR_mag_site_info = .true.
      if(len_trim(arg_string(i)) == 7 .and. arg_string(i)(1:7) == 'PCR_MAG') WRITE_PCR_mag_site_info = .true.
      if(arg_string(i)(1:3) == 'ALL')     site_info_all_atoms = .true.
     end do
    end if

    IF(nb_arg==0 .or. arg_string(1)(1:3) == 'ALL') then
     site_info_all_atoms = .true.
    ELSEIF(nb_arg == 1 .and. len_trim(arg_string(1))==3 .and. arg_string(1)(1:3) == 'PCR' ) then
     WRITE_PCR_mag_site_info = .false.
     WRITE_PCR_site_info     = .true.
     site_info_all_atoms     = .true.
    ELSEIF(nb_arg == 1 .and. len_trim(arg_string(1))==7 .and. arg_string(1)(1:7) == 'PCR_MAG' ) then
     WRITE_PCR_mag_site_info = .true.
     WRITE_PCR_site_info     = .false.
     site_info_all_atoms     = .true.
    ELSEIF(nb_arg == 1 .and. len_trim(arg_string(1))==3 .and. arg_string(1)(1:3) == 'MAG' ) then
     WRITE_PCR_mag_site_info = .true.
     WRITE_PCR_site_info     = .false.
     site_info_all_atoms     = .true.
    else
     site_info_all_atoms = .false.
     nb_atom_site_info = nb_arg
     READ(arg_string(1:nb_arg),*) site_info_label(1:nb_arg)
     site_info_label(1:nb_arg) = ADJUSTL(site_info_label(1:nb_arg))
    endif

   CASE ('SHIFT_2TH', 'SHIFT_2THETA', '2TH_SHIFT', '2THETA_SHIFT')
    IF(nb_arg == 0) then
     call check_arg_nb('0', 'SHIFT_2TH')
     return
    endif

    keyword_SH_2th = .true.
    READ(arg_string(1),*, iostat=i_error) var(1)
    IF(i_error /=0) then
     call error_message('SHIFT_2TH')
     return
    endif

    if(nb_arg > 1) then
     READ(arg_string(2),*, iostat=i_error) var(2)
     IF(i_error /=0) then
      call error_message('SHIFT_2TH')
      return
     endif
    endif

    if(nb_arg > 2) then
     READ(arg_string(3),*, iostat=i_error) var(3)
     IF(i_error /=0) then
      call error_message('SHIFT_2TH')
      return
     endif
    endif

    shift_2theta(1) = var(1)
    shift_2theta(2) = var(2)
    shift_2theta(3) = var(3)

    call write_info('')
    WRITE(message_text,'( a,F10.5)') '  > Shift 2theta (deg):', shift_2theta(1)
    call write_info(TRIM(message_text))
    WRITE(message_text,'( a,F10.5)') '     . cosTheta shift :', shift_2theta(2)
    call write_info(TRIM(message_text))
    WRITE(message_text,'( a,F10.5)') '     . sinTheta shift :', shift_2theta(3)
    call write_info(TRIM(message_text))
    keyword_SH_2th = .true.


   CASE ('HELP', 'MAN')
    keyword_HELP = .true.
    nb_help = nb_help_max
    !IF(nb_col /=0) then
    ! nb_help = nb_col
    ! read(read_line(long_kw+1:), *) HELP_string(1:nb_help)
    !endif
    IF(nb_arg >1) then
     nb_help = nb_arg
     READ(arg_string(1:nb_arg),*) HELP_arg(1:nb_arg)
    ELSEIF(nb_arg == 1) then
     i1 = INDEX(arg_string(1), '*')
     IF(i1==0) then
      nb_help = 1
      READ(arg_string(1), *) HELP_arg(1)
      return
     else
      IF(i1==1) then
       READ(arg_string(1)(i1+1:),*) arg1
      else
       READ(arg_string(1)(1:i1-1),*) arg1
      endif
      nb_help = 0
      temp_string = help_string
      do i=1, nb_help_max
       i1=INDEX(temp_string(i),TRIM(arg1))
       if (i1/=0) then
        nb_help= nb_help+1
        help_arg(nb_help) = temp_string(i)
       end if
      end do
     endif
    END if

   case ('MAN_HTML', 'HTML_MAN', 'HTML')
    keyword_create_CRYSCALC_HTML = .true.
    browse_cryscalc_HTML         = .false.
    IF(arg_string(1)(1:6) == 'BROWSE') browse_cryscalc_HTML = .true.

   case ('MAN_EXPERT', 'HELP_EXPERT')
    if(.not. expert_mode) then
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    endif

   case ('&MAN_EXPERT', '&HELP_EXPERT')

   case ('NEWS')
    keyword_NEWS     = .true.
	news_only_expert = .false.
    !if(nb_arg /=0) then
    ! read(arg_string(1), *) news_year
    !else
    ! news_year = 'all'
    !end if
	if(nb_arg == 0) then
	 news_year = 'all'
	else
	 do i=1, nb_arg
	  if(arg_string(i) == 'EXPERT_ONLY'      .or. arg_string(i) == 'ONLY_EXPERT'      .or. &
	     arg_string(i) == 'EXPERT_MODE'      .or. arg_string(i) == 'EXPERT_MODE_ONLY' .or. &
		 arg_string(i) == 'ONLY_EXPERT_MODE' .or. arg_string(i) == 'ONLY_MODE_EXPERT' .or. &
		 arg_string(i) == 'EXPERT') then
	   news_only_expert = .true.
	  else
       read(arg_string(i), *) news_year
	  end if
	 end do
    end if


   case ('HEADER', 'HEAD')
    keyword_HEADER = .true.

   CASE ('KEY', 'KEYS', 'LST_KEYS', 'LIST_KEYS', 'LST_KEYWORDS', 'LIST_KEYWORDS')
    keyword_KEY = .true.
    IF(nb_arg == 0) then
     write_keys(1: nb_help_max) = .true.
    ELSE
     READ(read_line(long_kw+1:), *) arg_line
     arg_line = adjustl(arg_line)
     i1=INDEX(arg_string(1), '*')
     if (i1==0) then
      write_keys(1:nb_help_max) = .false.
      do i=1, nb_help_max
       do j=1, nb_arg
        if (arg_string(j) == help_string(i)) write_keys(i) = .true.
       END do
      end do
      return
     elseif(i1 == 1) then
      READ(arg_string(1)(2:),*) arg1
      temp_string = help_string
      do i=1, nb_help_max
       i1 = INDEX(temp_string(i), TRIM(arg1))
       if (i1/=0) then
        write_keys(i) = .true.
       else
        write_keys(i) = .false.
       end if
      end do

     elseif(i1 > 1) then
      READ(arg_string(1)(1:i1-1),*) arg1
      temp_string = help_string
      do i=1, nb_help_max
       i1 = INDEX(temp_string(i), TRIM(arg1))
       if (i1/=0) then
        write_keys(i) = .true.
       else
        write_keys(i) = .false.
       end if
      end do
     end if
    endif


   CASE ('LIST_SG', 'LST_SG', 'LIST_SPACE_GROUPS')
    list_sg(1:7)         = .false.
    list_sg_Bravais(1:7) = .true.
    list_sg_centric(1:2) = .false.
    list_sg_laue(1:14)   = .true.
    list_sg_multip       = .false.
    list_sg_enantio      = .false.
    list_sg_chiral       = .false.
    list_sg_polar        = .false.

    IF(nb_arg == 0) then
     list_sg(1:7)  = .true.
     keyword_LSPGR = .true.
     return
    endif


    do i=1, nb_arg
     arg = l_case(arg_string(i))

     select case (arg)
       case ('all')
         list_sg(1:7) = .true.

       case ('tri', 'tric', 'tricl', 'tricli', 'triclini', 'triclinic')
         list_sg(1) = .true.

       case ('mono', 'monoc', 'monocl', 'monocli', 'monoclini', 'monoclinic')
         list_sg(2) = .true.

       case ('ortho', 'orthor', 'orthorh', 'orthorho', 'orthorhom', 'orthorhomb', 'orthorhombi', 'orthorhombic')
         list_sg(3) = .true.

       case ('tetra',  'tetrag',  'tetrago',  'tetragon',  'tetragona',  'tetragonal', &
             'quadra', 'quadrat', 'quadrati', 'quadratiq', 'quadratiqu', 'quadratique')
         list_sg(4) =  .true.

       case ('trig', 'trigo', 'trigon', 'trigona', 'trigonal')
         list_sg(5) = .true.

       case ('hex', 'hexa', 'hexag', 'hexago', 'hexagon', 'hexagona', 'hexagonal')
         list_sg(6) = .true.

       case ('cub', 'cubi', 'cubic')
         list_sg(7) = .true.

       case ('centric', 'centro')
         list_sg_centric(1:2) = .false.
         list_sg_centric(1) = .true.
         IF(nb_arg == 1) list_sg(1:7) = .true.

       case ('acentric', 'non-centro')
        list_sg_centric(1:2) = .false.
        list_sg_centric(2) = .true.
        IF(nb_arg == 1) list_sg(1:7) = .true.

       case ('enantio', 'enantiomorphic')
        list_sg_enantio = .true.
        IF(nb_arg == 1) list_sg(1:7) = .true.

       case('chiral')
        list_sg_chiral = .true.
        if(nb_arg == 1) list_sg(1:7) = .true.

       case('polar')
        list_sg_polar = .true.
        if(nb_arg == 1) list_sg(1:7) = .true.

       case ('mult', 'multip', 'multipl', 'multiplicity')
        list_sg_multip = .true.
        IF(nb_arg == 1) list_sg(1:7) = .true.

       case ('p')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(1)   = .true.
       case ('a')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(2)   = .true.
       case ('b')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(3)   = .true.
       case ('c')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(4)   = .true.
       case ('i')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(5)   = .true.
       case ('f')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(6)   = .true.
       case ('r')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(7)   = .true.

       case ('laue_1')    ! -1
        list_sg_laue(1:14) = .false.
        list_sg_laue(1)    = .true.
        list_sg(1)         = .true.
       case ('laue_2')    ! 2/m
        list_sg_laue(1:14) = .false.
        list_sg_laue(2)    = .true.
        list_sg(2)         = .true.
       case ('laue_3')    ! mmm
        list_sg_laue(1:14) = .false.
        list_sg_laue(3)    = .true.
        list_sg(3)         = .true.
       case ('laue_4')    ! 4/m
        list_sg_laue(1:14) = .false.
        list_sg_laue(4)    = .true.
        list_sg(4)         = .true.
       case ('laue_5')   ! 4/mmm
        list_sg_laue(1:14) = .false.
        list_sg_laue(5)    = .true.
        list_sg(4)         = .true.
       case ('laue_6')   ! -3 (rhomb. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(6)    = .true.
        list_sg(5)         = .true.
       case ('laue_7')   ! -3 (hex. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(7)    = .true.
        list_sg(5)         = .true.
       case ('laue_8')  ! -3m (rhomb. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(8)    = .true.
        list_sg(5)         = .true.
       case ('laue_9')  ! -31m (hex. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(9)    = .true.
        list_sg(5)         = .true.
       case ('laue_10') ! -3m1 (hex. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(10)   = .true.
        list_sg(5)         = .true.
       case ('laue_11') ! 6/m
        list_sg_laue(1:14) = .false.
        list_sg_laue(11)    = .true.
        list_sg(6)         = .true.
       case ('laue_12') ! 6/mmm
        list_sg_laue(1:14) = .false.
        list_sg_laue(12)    = .true.
        list_sg(6)         = .true.
       case ('laue_13') ! m3
        list_sg_laue(1:14) = .false.
        list_sg_laue(13)    = .true.
        list_sg(7)         = .true.
       case ('laue_14') ! m3m
        list_sg_laue(1:14) = .false.
        list_sg_laue(14)    = .true.
        list_sg(7)         = .true.

     end select

    end do
    keyword_LSPGR = .true.

   case ('LST_LAUE', 'LIST_LAUE', 'LST_LAUE_CLASS', 'LIST_LAUE_CLASS')
    keyword_laue = .true.

   case ('PAUSE')
    keyword_PAUSE = .true.

   case ('REM', '!', '#')
    keyword_REM = .true.


   case ('PERMUT', 'PERMUTATION', 'PERMUT_ABC', 'PERMUTATION_ABC')
    write_permutation_abc= .true.

   CASE ('TRICLINIC', 'TRICL')
    write_triclinic_transf = .true.

   CASE ('MONOCLINIC', 'MONOC', 'MONOCL')
    WRITE_monoclinic_transf = .true.
    !WRITE_monoclinic_out    = .true.
    !if(nb_arg /=0) then
    ! if(len_trim(arg_string(1)) >= 6) then
    !  if(u_case(arg_string(1)(1:6)) == 'NO_OUT' .or. u_case(arg_string(1)(1:6)) == 'NO_DET') write_monoclinic_out = .false.
    ! end if
    !endif

   CASE ('HEXA_TWIN', 'HEXA_TWINNING',  'HEXAGONAL_TWIN', 'HEXAGONAL_TWINNING', &
         'TWIN_HEXA', 'TWIN_HEXAGONAL', 'TWINNING_HEXA',  'TWINNING_HEXAGONAL')
    WRITE_twin_hexa = .true.

   CASE ('TWIN_PSEUDO_HEXA')
    WRITE_twin_pseudo_hexa = .true.


   CASE ('OBV_REV', 'OBVERSE_REVERSE', 'TWIN_OBVERSE_REVERSE', 'TWIN_OBV_REV', &
            'TWINNING_OBVERSE_REVERSE', 'TWINNING_OBV_REV')
    keyword_OBV_REV = .true.
    if (nb_arg/=0) then
     READ(arg_string(1), *) OBV_REV_twin_matrix
     IF(OBV_REV_twin_matrix /=2)  OBV_REV_twin_matrix = 1
    end if

   CASE ('RHOMB_HEX', 'RHOMB_HEXA', 'RHOMB_TO_HEX', 'RHOMB_TO_HEXA')
    WRITE_rhomb_hex_transf = .true.

   CASE ('HEX_RHOMB', 'HEXA_RHOMB', 'HEX_TO_RHOMB', 'HEXA_TO_RHOMB')
    WRITE_hex_rhomb_transf = .true.

   case ('SHAN', 'SHANNON')
    shannon_atom_label = ''
    !IF(nb_arg == 0) then
    ! call check_arg_nb('0', 'SHANNON')
    ! return
    !endif
    IF(nb_arg /=0 ) READ(arg_string(1),*) shannon_atom_label
    keyword_SHANNON = .true.

   case ('MAG', 'MAGNETISM', 'MAGNETIC')
    mag_atom_label = ''
    IF(nb_arg /=0) READ(arg_string(1), *) mag_atom_label
    keyword_MAG = .true.

   case ('MENDEL')
    mendel_atom_nb   = 0
    mendel_plot      = .false.
    IF(nb_arg == 0 .or. arg_string(1) == 'PLOT')  then
     call check_arg_nb('0', 'MENDEL')
     return
    endif

    if (arg_string(nb_arg) == 'PLOT') then
     mendel_plot = .true.
     mendel_atom_nb = nb_arg - 1
    else
     mendel_atom_nb = nb_arg
    endif
    READ(arg_string(1:nb_arg),*) mendel_atom_label(1:nb_arg)
    keyword_MENDEL = .true.

   case ('DATA_NEUTRONS', 'DATA_NEUTRON', 'NEUTRON_DATA', 'NEUTRONS_DATA')
    data_neutrons_PLOT         = .false.
    DATA_neutrons_RE_PLOT_ALL  = .false.
    data_n_RE                  = .false.
    data_neutrons_RE(1:14)     = .false.
    if(nb_arg/=0) then
     !data_neutrons_RE(1:10) = .false.
     !DATA_n_RE        = .false.
     do i=1, nb_arg
       long = len_trim(arg_string(i))
       if(long >=4 .and. arg_string(i)(1:4) == 'PLOT')      data_neutrons_PLOT        = .true.
       if(long >=8 .and. arg_string(i)(1:8) == 'PLOT_ALL')  DATA_neutrons_RE_PLOT_ALL = .true.
       !if(arg_string(i) == 'PLOT')   data_neutrons_PLOT = .true.

       if(arg_string(i) == 'SM_NAT') then
        data_neutrons_RE(1)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'SM_149') then
        data_neutrons_RE(2)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'EU_NAT') then
        data_neutrons_RE(3)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'EU_151') then
        data_neutrons_RE(4)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'GD_NAT') then
        data_neutrons_RE(5)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'GD_155') then
        data_neutrons_RE(6)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'GD_157') then
        data_neutrons_RE(7)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'DY_164') then
        data_neutrons_RE(8)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'ER_NAT') then
        data_neutrons_RE(9)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'ER_167') then
        data_neutrons_RE(10)  = .true.
        DATA_n_RE = .true.
       end if
       if(arg_string(i) == 'YB_NAT') then
        data_neutrons_RE(11) = .true.
        DATA_n_RE = .true.
       endif
       if(arg_string(i) == 'YB_168') then
        data_neutrons_RE(12) = .true.
        DATA_n_RE = .true.
       endif
       if(arg_string(i) == 'YB_174') then
        data_neutrons_RE(13) = .true.
        DATA_n_RE = .true.
       endif
       if(arg_string(i) == 'LU_176') then
        data_neutrons_RE(14) = .true.
        DATA_n_RE = .true.
       endif
     end do

     do i=1, 14
      if(data_neutrons_RE(i)) then
       current_RE   = i
       current_RE_n = RE_n(i)
       current_RE_label = RE_label(i)
       exit
      end if
     end do
    end if
    keyword_DATA_NEUTRONS = .true.

   case ('DATA_XRAYS', 'DATA_XRAY', 'XRAYS_DATA', 'XRAY_DATA')
    data_xrays_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_Xrays_PLOT = .true.
    keyword_DATA_XRAYS = .true.

   case ('DATA_DENSITY', 'DENSITY_DATA', 'DATA_ATOMIC_DENSITY', 'ATOMIC_DENSITY')
    data_atomic_density_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_atomic_density_PLOT = .true.
    keyword_DATA_ATOMIC_DENSITY = .true.

   case ('DATA_RADIUS', 'RADIUS_DATA', 'DATA_ATOMIC_RADIUS', 'ATOMIC_RADIUS')
    data_atomic_radius_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_atomic_radius_PLOT = .true.
    keyword_DATA_ATOMIC_RADIUS  = .true.

   case ('DATA_WEIGHT', 'WEIGHT_DATA', 'DATA_ATOMIC_WEIGHT', 'ATOMIC_WEIGHT')
    data_atomic_weight_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_atomic_weight_PLOT = .true.
    keyword_DATA_ATOMIC_WEIGHT = .true.


   case ('UPDATE', 'DOWNLOAD')
    keyword_download_EXE = .true.

   case ('WEB', 'INTERNET')
    IF(nb_arg == 0)  then
     call check_arg_nb('0', 'WEB')
     return
    endif
    URL_address = ''
    if (arg_string(1)(1:5) == 'CDIFX') then
     URL_address = 'http://www.cdifx.univ-rennes1.fr'
    elseif(arg_string(1)(1:8) == 'CRYSCALC') then
     URL_address = "http://"//trim(CRYSCALC%url)
    else
     long1 = LEN_TRIM(arg_string(1))
     do i = 1, WEB%num_site
      long2 = LEN_TRIM(WEB%NAME(i))
      if (arg_string(1)(1:long1) == WEB%NAME(i)(1:long2)) then
       URL_address = WEB%address(i)
       exit
      end if
     end do
     IF(LEN_TRIM(URL_address) == 0) URL_address = arg_string(1)
    end if
    keyword_WEB = .true.

   case ('WRITE_CELL', 'OUTPUT_CELL', '&WC')
    !IF(.NOT. keyword_CELL) then
    ! call write_info('')
    ! call write_info('  Cell parameters has to be known for WRITE_CELL keyword')
    ! call write_info('')
    ! return
    !endif
    if(nb_arg /=0) then
     write_cell_cart = .false.
     do i=1, nb_arg
      if(arg_string(i)(1:4) == 'CART') then
       write_cell_cart = .true.
       if (arg_string(i)(6:6) == "A") then
        cartesian_frame%type = "A"
        cartesian_frame%string = "a // x"
       elseif (arg_string(i)(6:6) == "C") then
        cartesian_frame%type = "C"
        cartesian_frame%string = "c // x"
       end if
      end if
     end do
    end if
    keyword_WRITE_CELL = .true.


   case ('WC')
    if(expert_mode)  then
     if(nb_arg /=0) then
      write_cell_cart = .false.
      do i=1, nb_arg
       if(arg_string(i)(1:4) == 'CART') then
        write_cell_cart = .true.
        if(len_trim(arg_string(i)) > 5 .and. arg_string(i)(5:5) == '_') then
         if (arg_string(i)(6:6) == "A") then
          cartesian_frame%type = "A"
          cartesian_frame%string = "a // x"
         elseif (arg_string(i)(6:6) == "C") then
          cartesian_frame%type = "C"
          cartesian_frame%string = "c // x"
         end if
        end if
       end if
      end do
     end if
     keyword_WRITE_CELL = .true.
    else
     call write_info('')
     call write_info(' > Unknown '//TRIM(input_keyword)//' keyword !')
     call write_info('')
    endif


   case ('WRITE_CHEM', 'WRITE_CHEMICAL_FORMULA', 'WRITE_MOLECULE', 'OUTPUT_CHEM', 'OUTPUT_CHEMICAL_FORMULA',  'OUTPUT_MOLECULE')
    IF(len_trim(molecule%formula) == 0) then
     call write_info('')
     call write_info('  Chemical formula not known !')
     call write_info('')
     return
    endif
    keyword_WRITE_CHEM = .true.

   case ('WRITE_QVEC', 'OUTPUT_QVEC')
    if(.not. keyword_QVEC) then
     call write_info('')
     call write_info('  QVEC has to be known for WRITE_QVEC keyword')
     call write_info('')
     return
    endif
    keyword_WRITE_QVEC = .true.

   case ('WRITE_SUPERCELL', 'OUTPUT_SUPERCELL')
    if(.not. keyword_SUPERCELL) then
     call write_info('')
     call write_info('  SUPERCELL has to be known for WRITE_SUPERCELL keyword')
     call write_info('')
     return
    endif
    keyword_WRITE_SUPERCELL = .true.

   case ('WRITE_SG', 'WRITE_SPACE_GROUP')
    IF(.NOT. keyword_SPGR) then
     call write_info('')
     call write_info('  SPGR keyword is mandatory for WRITE_SG keyword')
     call write_info('')
     return
    endif
    list_sg(1:7)  = .true.
    keyword_WRITE_SG = .true.


   case ('WRITE_WAVE', 'OUTPUT_WAVE')
    IF(.NOT. keyword_wave) then
     call write_info('')
     call write_info('  Wave has to be known for WRITE_WAVE keyword')
     call write_info('')
     return
    endif
    keyword_WRITE_WAVE = .true.

   case ('WRITE_Z', 'WRITE_ZUNIT')
    if(.not. keyword_Zunit) then
     call write_info('')
     call write_info('  ZUNIT has to be known for WRITE_ZUNIT keyword')
     call write_info('')
     return
    endif
    keyword_write_zunit = .true.

   case ('WRITE_DEVICE',         'OUTPUT_DEVICE',         'DEVICE',     &
         'WRITE_DIFFRACTO',      'OUTPUT_DIFFRACTO',      'DIFFRACTO',  &
         'WRITE_DIFFRACTOMETER', 'OUTPUT_DIFFRACTOMETER', 'DIFFRACTOMETER')
    keyword_WRITE_DEVICE =.true.

   case ('WRITE_BEAM', 'WRITE_INCIDENT_BEAM')
    keyword_WRITE_BEAM = .true.


   case ('NIGGLI', 'NIGGLI_CELL')
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info(' CELL paramters has to be known')
     call write_info('')
     return
    endif
    keyword_NIGGLI = .true.


   case ('KAPPA', 'KAPPA_TO_EULER')
    motor = 0.
    if(nb_arg ==3) then
     do i=1, nb_arg
      read(arg_string(i), *, iostat=i_error) var(i)
      if(i_error ==0) motor(i) = var(i)
     end do
    end if
    keyword_kappa_to_euler = .true.

   case ('EULER', 'EULER_TO_KAPPA')
    motor = 0.
    if(nb_arg ==3) then
     do i=1, nb_arg
      read(arg_string(i), *, iostat=i_error) var(i)
      if(i_error == 0) motor(i) = var(i)
     end do
    end if
    keyword_Euler_to_kappa = .true.

! references --------- !
   case ('REF_KCCD', 'KCCD')
    keyword_WRITE_REF_KCCD = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_APEX', 'REF_APEXII' , 'WRITE_APEX', 'WRITE_APEXII',  'APEX', 'APEXII')
    keyword_WRITE_REF_APEX = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_D8V_CU')
    keyword_WRITE_REF_D8V_Cu = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_D8V_MO')
    keyword_WRITE_REF_D8V_Mo = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_X2S', 'REF_SMART_X2S')
    keyword_WRITE_REF_X2S = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_SUPERNOVA')
    keyword_WRITE_REF_SUPERNOVA = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_XCALIBUR')
    keyword_WRITE_REF_XCALIBUR = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_EVAL', 'REF_EVALCCD', 'WRITE_EVAL', 'WRITE_EVALCCD', 'EVAL', 'EVALCCD')
    keyword_WRITE_REF_EVAL = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_DENZO', 'WRITE_DENZO', 'DENZO')
    keyword_WRITE_REF_DENZO = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_SHELX')
    keyword_WRITE_REF_SHELX = .true.
	if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

    case ('REF_SIR')
     keyword_WRITE_REF_SIR  = .true.

    case ("REF_SPF", "REF_SUPERFLIP")
     keyword_WRITE_REF_SPF  = .true.

   case ('REF_SADABS', 'REF_SAD', 'SADABS')
    keyword_WRITE_REF_SADABS = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)

   case ('REF_ABS_CRYSALIS')
    keyword_WRITE_REF_ABS_CRYSALIS = .true.
    if(nb_arg /=0) call check_write_CIF(nb_arg, arg_string)
! ------------------- !


   case ('READ_NREPORT', 'READ_NREPORT_HTML', 'READ_HTMLREPORT')
    keyword_read_NREPORT = .true.
    browse_nreport       = .false.
    if(nb_arg /=0 .and. arg_string(1) (1:6) == 'BROWSE') browse_nreport = .true.

   case ('READ_CEL', 'READ_CEL_FILE', 'READ_POWDERCELL')
    if(nb_arg == 0) then
     call check_arg_nb('0', 'READ_CEL')
     return
    endif
    i1 = index(arg_string(1), '.')
    if (i1 == 0) then
     CEL_file_name = trim(arg_string(1))//'.CEL'
    else
     READ(arg_string(1), *) CEL_file_name
    endif
    call test_file_exist(trim(CEL_file_name), file_exist, 'out')
    if(.not. file_exist) then
     CEL_file_name = ''
     return
    endif
    create_CIF_pymol = .false.
    if (nb_arg > 1) then
     do i=1, nb_arg
      if(arg_string(i)(1:5) == 'PYMOL') then
       create_CIF_PYMOL = .true.
      end if
     end do
     i1 = index(CEL_file_name, '.', back=.true.)
     write(CIF_pymol_file_name, '(a)') CEL_file_name(1:i1-1)//'_'//trim(CEL_file_name(i1+1:))//'_pml.CIF'
    end if

    keyword_read_CEL  = .true.

   case ('READ_CIF', 'READ_CIF_FILE', 'CIF_FILE')
    IF(nb_arg == 0)  then
     call check_arg_nb('0', 'READ_CIF')
     return
    endif
    i1 = INDEX(arg_string(1), '.')
    if (i1==0) then
     CIF_file_name = TRIM(arg_string(1))//'.CIF'
    else
     READ(arg_string(1), *) CIF_file_name
    end if
    call test_file_exist(TRIM(CIF_file_name), file_exist, 'out')
    if(.not. file_exist) then
     CIF_file_name = ''
     keyword_read_CIF = .false.
     return
    endif
    create_CIF_pymol = .false.
    if (nb_arg > 1) then
     do i=1, nb_arg
      if(arg_string(i)(1:5) == 'PYMOL') then
       create_CIF_PYMOL = .true.
      end if
     end do
     i1 = index(CIF_file_name, '.', back=.true.)
     write(CIF_pymol_file_name, '(a)') CIF_file_name(1:i1-1)//'_pml'//trim(CIF_file_name(i1:))
    end if

    keyword_read_CIF = .true.
    keyword_CELL     = .true.

   case ('WRITE_BONDS', 'WRITE_ANGLES', 'WRITE_TORSION', 'WRITE_TORSION_ANGLES', 'WRITE_HTAB')
    if(expert_mode) then
     if(.not. keyword_read_cif) then
      call write_info('')
      call write_info('  No input CIF file !!')
      call write_info('')
     end if
    end if

   case ('READ_FACES', 'FACES')
    !if(nb_arg == 0) then
    ! call check_arg_nb('0', 'READ_FACES')
    ! return
    !end if
    if(nb_arg == 0) then
     FACES_file_name = "absorb.ins"
    else
     i1 = index(arg_string(1), '.')
     if(i1 == 0) then
      FACES_file_name = trim(arg_string(1))//'.INS'
     else
      read(arg_string(1), *, iostat=i_error) FACES_file_name
     end if
    end if
    call test_file_exist(trim(FACES_file_name), file_exist, 'out')
    if(.not. file_exist) then
     FACES_file_name = ''
     return
    end if
    keyword_READ_FACES = .true.


   case ('READ_INS', 'READ_INS_FILE', 'INS_FILE')
    IF(nb_arg == 0)  then
     call check_arg_nb('0', 'READ_INS')
     return
    endif
    i1 = INDEX(arg_string(1), '.')
    if (i1==0) then
     INS_file_name = TRIM(arg_string(1))//'.INS'
    else
     READ(arg_string(1), *, iostat=i_error) INS_file_name
    end if
    create_CIF_pymol = .false.
    read_Q_peaks     = .false.
    if (nb_arg > 1) then
     do i=1, nb_arg
      long = len_trim(arg_string(i))
      if(long==5 .and. arg_string(i)(1:5) == 'PYMOL') then
       create_CIF_PYMOL = .true.
      elseif(long==7 .and. arg_string(i)(1:7) == 'Q_PEAKS') then
       read_Q_peaks = .true.
      end if
     end do
     i1 = index(INS_file_name, '.', back=.true.)
     write(CIF_pymol_file_name, '(a)') INS_file_name(1:i1-1)//'_'//trim(INS_file_name(i1+1:))//'_pml.CIF'
    end if

    call test_file_exist(trim(INS_file_NAME), file_exist, 'out')
    IF(.NOT. file_exist) then
     INS_file_name = ''
     return
    endif

    keyword_read_INS = .true.
    keyword_CELL     = .true.

   case ('READ_PCR', 'READ_PCR_FILE', 'PCR_FILE')
    IF(nb_arg == 0)  then
     call check_arg_nb('0', 'READ_PCR')
     return
    endif
    i1 = INDEX(arg_string(1), '.PCR')
    if (i1==0) then
     PCR_file_name = TRIM(arg_string(1))//'.PCR'
    else
     READ(arg_string(1), *, iostat=i_error) PCR_file_name
    end if
    call test_file_exist(trim(PCR_file_name), file_exist, 'out')
    IF(.NOT. file_exist) then
     PCR_file_name = ''
     return
    endif
    keyword_read_PCR = .true.


  case ('ABIN')
   if(expert_mode) then
    IF(nb_arg == 0)  then
     call check_arg_nb('0', 'ABIN')
     return
    endif
    i1 = INDEX(arg_string(1), '.')
    if (i1==0) then
     INS_file_name = TRIM(arg_string(1))//'.INS'
    else
     READ(arg_string(1), *, iostat=i_error) INS_file_name
    end if
    call test_file_exist(trim(INS_file_NAME), file_exist, 'out')
    IF(.NOT. file_exist) then
     INS_file_name = ''
     return
    endif
    keyword_ABIN = .true.
   end if

  case ('READ_TIDY_OUT')
    IF(nb_arg == 0)  then
     !call check_arg_nb('0', 'READ_TIDY_OUT')
     !return
      TIDY_out_file_name = 'stidy.out'
    else
     READ(arg_string(1), *, iostat=i_error) TIDY_out_file_name
    endif

    call test_file_exist(trim(TIDY_out_file_name), file_exist, 'out')
    IF(.NOT. file_exist) then
     TIDY_out_file_name = ''
     return
    endif
    keyword_read_TIDY_out = .true.


   case ('SIR', 'SIR97')
    keyword_SIR = .true.

   case ('XRAYS_WAVELENGTH', 'X_WAVE')
    keyword_X_WAVE = .true.
    IF(nb_arg /= 0) then
     X_target(1: tabulated_target_nb)%write = .false.
     do i=1, tabulated_target_nb
      do j=1, nb_arg
       if (arg_string(j)(1:2) == u_case(X_target(i)%label)) X_target(i)%write = .true.
      END do
     end do
    else
     X_target(1: tabulated_target_nb)%write = .true.
    endif

   case('VER', 'VERSION')
    keyword_version  = .true.

   CASE ('END', 'EXIT', 'QUIT')
    return

  ! CASE default
  !  call write_info('')
  !  call write_info('  > Unknown '//TRIM(input_keyword)//' keyword !')
  !  call write_info('')

   case default

    unknown_keyword = .true.

  end select


 return
END subroutine identification_keywords

!------------------------------------------------------------------------------


subroutine check_arg_nb(input_string, input_keyword)
 USE IO_module, ONLY : write_info
 implicit none
  CHARACTER (LEN=*),    INTENT(IN) :: input_string
  CHARACTER (LEN=*),    INTENT(IN) :: input_keyword

  select case (input_string)
      case ('0')
     call write_info('')
     call write_info('  Warning: Argument(s) is(are) mandatory for '//TRIM(input_keyword)// ' keyword ...')
     call write_info('')

      case ('-1')
     call write_info('')
     call write_info('  Warning: At least one argument is mandatory for '//TRIM(input_keyword)// ' keyword ...')
     call write_info('')

      case ('9')
     call write_info('')
     call write_info('  Warning: Arguments (components of the 3*3 matrix) are mandatory for MATRIX keyword ...')
     call write_info('')

      case ('unknown')
     call write_info('')
     call write_info(' Warning: Unknown argument for '//TRIM(input_keyword)// ' ...')
     call write_info('')

      case default
     call write_info('')
     call write_info('  Warning: '//TRIM(input_string)//' arguments are mandatory for '//TRIM(input_keyword)// ' keyword ...')
     call write_info('')

  end select

end subroutine check_arg_nb

!------------------------------------------------------------------------------------------------

subroutine get_SPG(input_string)
 use cryscalc_module,                only  : space_group_symbol, SPG, unit_cell
 use macros_module,                  only  : l_case
 use CFML_crystallographic_symmetry, only  : set_spacegroup

 implicit none
  character (len=*), intent(in)  :: input_string
  integer                        :: long_input_string
  LOGICAL                        :: input_tricl, input_mono, input_ortho, input_tetra, input_trig, input_hexa, input_cub

  input_tricl = .false.
  input_mono  = .false.
  input_ortho = .false.
  input_tetra = .false.
  input_trig  = .false.
  input_hexa  = .false.
  input_cub   = .false.

  long_input_string = len_trim(input_string)


  select case(l_case(trim(input_string)))

    case('tricl', 'triclinic')
     input_tricl = .true.

    case('mono', 'monocl', 'monoclinic', 'monoclinique')
     input_mono = .true.

    case('trig', 'trigonal')
     input_trig = .true.

    case('hexa', 'hexagonal')
     input_hexa = .true.

    case('ortho', 'orthorhombic', 'orthorhombique')
     input_ortho = .true.

    case('tetra', 'tetrag', 'tetragonal')
     input_tetra = .true.

    case('cub', 'cubic', 'cubique')
     input_cub = .true.

    case default
     input_tricl = .true.

  end select

  !if(long_input_string == 4) then
  ! if(input_string(1:4) == 'TRIC') input_tricl = .true.
  ! if(input_string(1:4) == 'MONO') input_mono  = .true.
  ! if(input_string(1:4) == 'TRIG') input_trig  = .true.
  ! if(input_string(1:4) == 'HEXA') input_hexa  = .true.
  !elseif(long_input_string == 5) then
  ! if(input_string(1:5) == 'ORTHO') input_ortho = .true.
  ! if(input_string(1:5) == 'TETRA') input_tetra = .true.
  !elseif(long_input_string == 3) then
  ! if(input_string(1:3) == 'CUB')   input_cub = .true.
  !endif


    IF(input_tricl) then
       space_group_symbol = 'P 1'
       call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_mono)  then
      space_group_symbol = 'P 2'
      if(unit_cell%H_M(1:1) == 'C') space_group_symbol = 'C 2'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_ortho) then
      space_group_symbol = 'P 2 2 2'
      if(unit_cell%H_M(1:1) == 'F') space_group_symbol = 'F 2 2 2'
      if(unit_cell%H_M(1:1) == 'I') space_group_symbol = 'I 2 2 2'
      if(unit_cell%H_M(1:1) == 'C') space_group_symbol = 'C 2 2 2'
      if(unit_cell%H_M(1:1) == 'B') space_group_symbol = 'B 2 2 2'
      if(unit_cell%H_M(1:1) == 'A') space_group_symbol = 'A 2 2 2'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_tetra) then
      space_group_symbol = 'P 4'
      if(unit_cell%H_M(1:1) == 'I') space_group_symbol = 'I 4'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_trig)  then
      space_group_symbol = 'P 3'
      if(unit_cell%H_M(1:1) == 'R') space_group_symbol = 'R 3'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_hexa)  then
      space_group_symbol = 'P 6'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_cub)   then
      space_group_symbol = 'P 2 3'
      if(unit_cell%H_M(1:1) == 'I') space_group_symbol = 'I 2 3'
      if(unit_cell%H_M(1:1) == 'F') space_group_symbol = 'F 2 3'

      call set_spacegroup(space_group_symbol, SPG)
     endif

 return
end subroutine get_SPG

!---------------------------------------------------------------------------------------------------------

subroutine get_matrix_coord(input_line)
 USE cryscalc_module, ONLY : mat
 use macros_module,   only : nombre_de_colonnes
 USE Definition_fractions

 implicit none
  CHARACTER (LEN=*), INTENT(INOUT)            :: input_line
  CHARACTER (LEN=32)                          :: string
  INTEGER                                     :: i1, i,j, k

  call def_fractions


  do i=1, 3
   do j=1, 3
    input_line = ADJUSTL(input_line)
    i1 = INDEX(input_line, ' ')

    string = input_line(1:i1-1)
    string = ADJUSTL(string)
    if(string(1:1) == '+') string = string(2:)
    do k = 1, nb_fraction
     IF(string(1:3) == ratio_string_pos(k)(1:3)) then
      Mat(i,j) = ratio_real_pos(k)
      exit

     ELSEIF(string(1:4) == ratio_string_neg(k)(1:4)) then
      Mat(i,j) = ratio_real_neg(k)
      exit

     else
      READ(string, *) Mat(i,j)
     endif
    end do
    input_line = input_line(i1+1:)
   end do   ! loop_i
  end do    ! loop_i

 return
end subroutine get_matrix_coord

!---------------------------------------------------------------------------------
subroutine verif_format(hkl_format)
 use macros_module, only : remove_car
 implicit none
  character (len=*) , intent(inout)  :: hkl_format
  character (len=32) :: new_format
  integer                            :: i1, long


  long = len_trim(hkl_format)

  i1 = index(hkl_format, "'")
  if(i1 /=0)   hkl_format = remove_car(hkl_format, "'")

  i1 = index(hkl_format, '"')
  if(i1 /=0)   hkl_format = remove_car(hkl_format, '"')

  i1 = index(hkl_format, '(')
  if(i1 /=0)   hkl_format = remove_car(hkl_format, '(')

  i1 = index(hkl_format, ')')
  if(i1 /=0)   hkl_format = remove_car(hkl_format, ')')

  write(new_format , '(3a)') '(', trim(hkl_format), ')'
  hkl_format = trim(new_format)


 return
end subroutine verif_format

subroutine get_arg_string(arg_line, nb_arg, arg_string)
 implicit none
  character (len=*),                 intent(in)    :: arg_line
  integer,                           intent(in)    :: nb_arg
  character (len=64), dimension(20), intent(inout) :: arg_string
  character (len=256)                              :: tmp_line
  integer                                          :: i, i1

  tmp_line = arg_line
  do i = 1, nb_arg
   tmp_line = adjustl(tmp_line)
   i1 = index(tmp_line, ' ')
   arg_string(i) = tmp_line(1:i1-1)
   tmp_line = tmp_line(i1:)
  end do

 return
end subroutine Get_arg_string

!----------------------------------------------------------------------
subroutine test_arg_string(input_string,arg_string, nb_arg, string_ok)
 use cryscalc_module, only  : cartesian_frame
 implicit none
  character (len=*)                     , intent(in)    :: input_string
  integer                               , intent(inout) :: nb_arg
  character (len=*), dimension(nb_arg+1), intent(inout) :: arg_string
  logical                               , intent(inout) :: string_ok
  integer                                               :: i, i1, i_string, l1, l2



    ! test if input_string is present
    i_string  = 0


    l1 = len_trim(input_string)
    do i = 1, nb_arg
     i1 = index(arg_string(i), '_')
     l2 = len_trim(arg_string(i))

     if(i1==0) then
      if(arg_string(i)(1:l2) == input_string(1:l1)) then
       string_ok = .true.
       i_string    = i
       exit
      end if
     else
      if(arg_string(i)(1:i1-1) == input_string(1:l1)) then
       string_ok = .true.
       i_string  = i
       if(arg_string(i)(1:7) == "SHAPE_A") then
        cartesian_frame%type = "A"
        cartesian_frame%string = "a // x"
       elseif (arg_string(i)(1:7) == "SHAPE_C") then
        cartesian_frame%type = "C"
        cartesian_frame%string = "c // x"
       end if
      end if
     endif
    end do

    if(string_ok) then
     if(i_string /= nb_arg) then
      do i=i_string, nb_arg
       arg_string(i) = arg_string(i+1)
      end do
     end if
     nb_arg = nb_arg - 1
    end if


 return
end subroutine test_arg_string
!-------------------------------------------------------------------------------------------
subroutine check_write_CIF(nb_arg, arg_string)
 use cryscalc_module, only : WRITE_ref_CIF
 implicit none
  integer, intent(in)                               :: nb_arg
  character (len=64), dimension(nb_arg), intent(in) :: arg_string
  integer                                           :: i, long



    WRITE_ref_CIF = .false.
    do i=1, nb_arg
     long = len_trim(arg_string(i))
     if(long >=3) then
      if(arg_string(i)(1:3) == 'CIF') then
       WRITE_ref_CIF = .true.
       exit
      end if
     end if
     if(long >=4) then
      if(arg_string(i)(1:4) == 'ACTA') then
       WRITE_ref_CIF = .true.
       exit
      end if
     end if
    end do

   return
end subroutine check_write_CIF
!-------------------------------------------------------------------------------------------
subroutine get_arg_ABx(arg_string, DIST_real, atm_1, atm_2, lecture_OK)
 implicit none
  character(len=*), dimension(3), intent(in)      :: arg_string
  real,                           intent(inout)   :: DIST_real
  character(len=8),               intent(inout)   :: atm_1, atm_2
  logical,                        intent(out)     :: lecture_OK
  integer                                         :: i_error
  character(len=8)                                :: loc_str
  real                                            :: loc_real
  real, parameter                                 :: eps = 0.0000001

  lecture_ok = .false.

  ! teste A B x
  read(arg_string(1), *, iostat = i_error) loc_str
  if(i_error /=0) return
  atm_1 = loc_str
  read(arg_string(2), *, iostat = i_error) loc_str
  if(i_error /=0) return
  atm_2 = loc_str
  read(arg_string(3), *, iostat = i_error) loc_real
  if(i_error /=0) then ! teste x A B
   read(arg_string(1), *, iostat = i_error) loc_real
   if(i_error /=0) return
   if(ABS(loc_real) > eps)  DIST_real = loc_real
   read(arg_string(2), *, iostat = i_error) loc_str
   if(i_error /=0) return
   atm_1 = loc_str
   read(arg_string(3), *, iostat = i_error) loc_str
   atm_2 = loc_str
  else
   if(ABS(loc_real) > eps)  DIST_real = loc_real
  end if

  lecture_ok = .true.

 return
end subroutine get_arg_ABx

!--------------------------------------------------------------------
subroutine Get_bary_atoms(arg_string, nb_arg, ok)
 use macros_module, only : check_character, remove_car
 USE IO_module,     only : write_info

 implicit none
 character(len=64), dimension(20), intent(inout) :: arg_string
 integer,                          intent(inout) :: nb_arg
 LOGICAL,                          intent(inout) :: ok
 integer                                         :: i, pos_num, i1, i2
 CHARACTER(len=8)                                :: bary_label
 LOGICAL                                         :: alpha_char, numeric_char

 ! decodage d'une commande du type C10 > C15

 alpha_char    = .false.
 numeric_char  = .false.
 do i=1, len_trim(arg_string(1))
  call check_character(arg_string(1)(i:i), alpha_char, numeric_char)
  if(numeric_char) then
   pos_num = i
   exit
  end if
 end do

 if(arg_string(1)(1:pos_num-1) /= arg_string(3)(1:pos_num-1)) then
  call write_info('! Wrong atomic labels !')
  ok = .false.
  return
 end if

 read(arg_string(1)(1:pos_num-1), *) bary_label

 read(arg_string(1)(pos_num:), *) i1
 read(arg_string(3)(pos_num:), *) i2
 if(i2>i1) then
  nb_arg = i2 - i1 +1
  do i=1, nb_arg
   write(arg_string(i), '(a,i4)') trim(bary_label),i1 + (i-1)
   arg_string(i) = remove_car(arg_string(i), ' ')
  end do
 else
  call write_info('! Wrong atomic labels !')
  ok = .false.
  return
 endif

 ok = .true.

 return
end subroutine Get_bary_atoms

!-------------------------------------------------------------------------------------------
