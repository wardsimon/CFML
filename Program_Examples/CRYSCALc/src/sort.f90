!     Last change:  TR   23 Apr 2007    7:00 pm
!     Last change:  TR   23 Apr 2007    7:00 pm
!
subroutine read_and_sort_HKL(input_string, input_string_2)
 USE cryscalc_module, ONLY     : nb_sort,  sort_type,  cut_off,             &
                                 nb_shell, shell_type, shell_plot,          &
                                 keyword_CELL, keyword_WAVE, keyword_FILE,  &
                                 HKL_unit, message_text, lecture_OK, WRITE_details, debug_proc
 USE macros_module,  ONLY      : test_file_exist
 USE IO_module,      ONLY      : write_info
 USE hkl_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_string
  CHARACTER (LEN=*), INTENT(INOUT)     :: input_string_2
  INTEGER                              :: i
  LOGICAL                              :: file_exist
  !LOGICAL                              :: lecture_ok
  CHARACTER (LEN=4)                    :: ext_string
  INTEGER                              :: long_input_string
  LOGICAL                              :: input_sort, input_shell

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_AND_SORT_HKL ("//trim(input_string)//")")


  !lecture_ok = .false.
  input_sort  = .false.
  input_shell = .false.

  long_input_string = len_trim(input_string)
  if(long_input_string == 4) then
   IF(input_string(1:4) == 'sort') input_sort = .true.
  elseif(long_input_string == 5) then
   IF(input_string(1:5) == 'shell') input_shell = .true.
  endif


  IF(.NOT. keyword_FILE) then
   call write_info('')
   call write_info('  >>> FILE is a mandatory keyword !')
   call write_info('')
   return
  endif
  IF(HKL_file%INT) THEN
   ext_string = '.int'
  ELSEIF(HKL_file%COL) then
   ext_string = '.col'
  ELSEIF(HKL_file%M91) THEN
   ext_string = '.m91'
  ELSEIF(HKL_file%M95) then
   ext_string = '.m95'
  else
   ext_string = '.hkl'
  endif

  call test_file_exist(HKL_file%NAME, file_exist, 'out')
  IF(.not. file_exist) then
   HKL_file%NAME = ''
   return
  endif

  if(WRITE_details) then
  call write_info(' ')
  call write_info('  >> HKL file name: '// TRIM(HKL_file%NAME))
  call write_info(' ')
  end if

  i = INDEX(HKL_file%NAME, ".")
  if (i/=0) then
   IF(input_sort) then
    HKL_file%d     = TRIM(HKL_file%NAME(1:i-1))//'_sort_d'//TRIM(ext_string)
    HKL_file%stl   = TRIM(HKL_file%NAME(1:i-1))//'_sort_stl'//TRIM(ext_string)
    HKL_file%theta = TRIM(HKL_file%NAME(1:i-1))//'_sort_theta'//TRIM(ext_string)
    HKL_file%I     = TRIM(HKL_file%NAME(1:i-1))//'_sort_i'//TRIM(ext_string)
    HKL_file%Isig  = TRIM(HKL_file%NAME(1:i-1))//'_sort_isig'//TRIM(ext_string)
   ELSEIF(input_shell) then
    HKL_file%d     = TRIM(HKL_file%NAME(1:i-1))//'_shell_d'//TRIM(ext_string)
    HKL_file%stl   = TRIM(HKL_file%NAME(1:i-1))//'_shell_stl'//TRIM(ext_string)
    HKL_file%theta = TRIM(HKL_file%NAME(1:i-1))//'_shell_theta'//TRIM(ext_string)
    HKL_file%I     = TRIM(HKL_file%NAME(1:i-1))//'_shell_i'//TRIM(ext_string)
    HKL_file%Isig  = TRIM(HKL_file%NAME(1:i-1))//'_shell_isig'//TRIM(ext_string)
   endif
  END if


  IF(nb_shell /=0) then
   OPEN(UNIT=HKL_unit, FILE=TRIM(HKL_file%name))

   do i=1, nb_shell
    select case (shell_type(i))
     case ('d')
      if(write_details) then
      call write_info(' ')
      call write_info('  >>> Create '// TRIM(HKL_file%d)// ' file ...')
      call write_info(' ')
      end if

      OPEN(UNIT=12, FILE=TRIM(HKL_file%d))
      IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
      IF(.NOT. lecture_ok) return
      call classement(i,'shell_d')
      CLOSE(UNIT=12)

     case ('stl')
      if (.not. keyword_CELL) then
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> cell parameters are mandatory to calculate SinTheta/lambda <<<')
       call write_info(' ')
       end if
       cycle
      ELSEIF(.NOT. keyword_WAVE) then
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> wavelength is a mandatory parameter to calculate SinTheta/lambda <<<')
       call write_info(' ')
       end if
       cycle
      end if

      if(WRITE_details) then
      call write_info(' ')
      call write_info('  >>> Create '// TRIM(HKL_file%stl)// ' file ...')
      call write_info(' ')
      end if

      OPEN(UNIT=12, FILE=TRIM(HKL_file%stl))
      IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
      IF(.NOT. lecture_ok) return
      call classement(i,'shell_stl')
      CLOSE(UNIT=12)

     case ('theta')
       if (.not. keyword_CELL) then
        if(WRITE_details) then
        call write_info(' ')
        call write_info('  >>> cell parameters are mandatory to calculate Theta angles <<<')
        call write_info(' ')
        end if
        cycle
       ELSEIF(.NOT. keyword_WAVE) then
        if(WRITE_details) then
        call write_info(' ')
        call write_info('  >>> wavelength is a mandatory parameter to calculate Theta angles <<<')
        call write_info(' ')
        end if
        cycle
       end if
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> Create '// TRIM(HKL_file%theta)// ' file ...')
       call write_info(' ')
       end if

       OPEN(UNIT=12, FILE=TRIM(HKL_file%theta))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
       IF(.NOT. lecture_ok) return
       call classement(i,'shell_theta')
       close (UNIT=12)

     case ('int')
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> Create '// TRIM(HKL_file%I)// ' file ...')
       call write_info(' ')
       end if


       OPEN(UNIT=12, FILE=TRIM(HKL_file%i))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
       IF(.NOT. lecture_ok) return
       call classement(i,'shell_int')
       CLOSE(UNIT=12)

     case ('isig')
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> Create '// TRIM(HKL_file%Isig)// ' file ...')
       call write_info(' ')
       end if

       OPEN(UNIT=12, FILE=TRIM(HKL_file%isig))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
       IF(.NOT. lecture_ok) return
       call classement(i,'shell_isig')
       CLOSE(UNIT=12)

     case default
    end select
   end do

  !endif
  !
  !if(nb_sort == 0) then
  elseif(nb_sort == 0) then
   OPEN(UNIT=HKL_unit, FILE=TRIM(HKL_file%NAME))
   IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
   IF(.NOT. lecture_ok) return
   call classement(0, 'no_sort')

  else
   do i=1, nb_sort
    OPEN(UNIT=HKL_unit, FILE=TRIM(HKL_file%NAME))

    select case (sort_type(i))
     case ('d')
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> SORT HKL file in increasing d(A) and create '// TRIM(HKL_file%d)// ' file ...')
       call write_info(' ')
       end if

       OPEN(UNIT=12, FILE=TRIM(HKL_file%d))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
       IF(.NOT. lecture_ok) return
       call classement(i,'sort_d')
       CLOSE(UNIT=12)

     case ('stl')
       if (.not. keyword_CELL) then
        if(WRITE_details) then
        call write_info(' ')
        call write_info('  >>> cell parameters are mandatory to calculate SinTheta/lambda <<<')
        call write_info(' ')
        end if
        cycle
       ELSEIF(.NOT. keyword_WAVE) then
        if(WRITE_details) then
        call write_info(' ')
        call write_info('  >>> wavelength is a mandatory parameter to calculate SinTheta/lambda <<<')
        call write_info(' ')
        end if
        cycle
       end if
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> SORT HKL file in increasing SinTheta/lambda and create '//TRIM(HKL_file%stl)// ' file ...')
       call write_info(' ')
       end if

       OPEN(UNIT=12, FILE=TRIM(HKL_file%stl))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
       IF(.NOT. lecture_ok) return
       call classement(i,'sort_stl')
       close (UNIT=12)

     case ('theta')
       if (.not. keyword_CELL) then
        call write_info(' ')
        call write_info('  >>> cell parameters are mandatory to calculate Theta angles <<<')
        call write_info(' ')
        cycle
       ELSEIF(.NOT. keyword_WAVE) then
        call write_info(' ')
        call write_info('  >>> wavelength is a mandatory parameter to calculate Theta angles <<<')
        call write_info(' ')
        cycle
       end if
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> SORT HKL file in increasing Theta and create '// TRIM(HKL_file%theta)// ' file ...')
       call write_info(' ')
       end if

       OPEN(UNIT=12, FILE=TRIM(HKL_file%theta))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
       IF(.NOT. lecture_ok) return
       call classement(i,'sort_theta')
       close (UNIT=12)

     case ('int')
       if(write_details) then
       call write_info(' ')
       call write_info('  >>> SORT HKL file in decreasing intensities and create '// TRIM(HKL_file%I)// ' file ...')
       call write_info(' ')
       end if

       OPEN(UNIT=12, FILE=TRIM(HKL_file%i))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
       IF(.NOT. lecture_ok) return
       call classement(i,'sort_int')
       CLOSE(UNIT=12)

     case ('isig')
       if(WRITE_details) then
       call write_info(' ')
       call write_info('  >>> SORT HKL file in decreasing I/sig and create '// TRIM(HKL_file%Isig)// ' file ...')
       call write_info(' ')
       end if

       OPEN(UNIT=12, FILE=TRIM(HKL_file%isig))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok, input_string_2)
       IF(.NOT. lecture_ok) return
       call classement(i,'sort_isig')
       CLOSE(UNIT=12)

     case default
    end select

   end do
  endif

  CLOSE(UNIT=HKL_unit)

 return
end subroutine read_and_sort_HKL

!----------------------------------------------------------------------------------------------------------------

subroutine read_data_file(ok, input_string)
 USE cryscalc_module,   ONLY         : cut_off, unit_cell, wave => wavelength,              &
                                       keyword_CELL, keyword_WAVE, message_text, HKL_unit,  &
                                       known_theta, keyword_create_CIF, hkl_format_free,    &
                                       hkl_format_SHELX, hkl_SHELX_format, hkl_format, crystal_system , UB_matrix, &
                                       WRITE_details, debug_proc
 USE HKL_module
 USE IO_module,         ONLY         : write_info
 USE macros_module,     ONLY         : nombre_de_colonnes, u_case
 USE CFML_Math_General, ONLY         : sort

 implicit none
 LOGICAL, INTENT(OUT)              :: ok
 character (len=*), intent(inout)  :: input_string
 integer                           :: i, i1, ier

 integer                           :: h_, k_, l_, m_, code_, numor
 integer                           :: nvk
 real                              :: F2_, sig_F2_
 REAL, DIMENSION(6)                :: cosdir
 LOGICAL                           :: HKL_ok
 real, parameter                   :: eps = 0.000001
 CHARACTER (LEN=256)               :: read_line
 CHARACTER (LEN=256)               :: M91_string
 REAL,    DIMENSION(5)             :: M95_r
 INTEGER, DIMENSION(3)             :: M95_i
 CHARACTER (LEN=256)               :: INT_string, format_
 real, dimension(4)                :: INT_real
 INTEGER                           :: div_10, nb_col
 REAL                              :: max_F2
 logical                           :: QVEC_type
 character (len=8)                 :: nref_string, out_string


 if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_DATA_FILE")

 div_10 = 0


  !lecture fichier .HKL format SHELX ou fichier .CIF
  ! + creation fichier avec sintheta/l dans la 1ere colonne

  if(WRITE_details) then
   call write_info('')
   if(.not. HKL_file%read_NEG) then
    call write_info('  ... Reading data (All negative intensities reflections will be replaced by 0.0001) ...')
   else
    call write_info('  ... Reading data (All negative intensities reflections will be treated as they are) ...')
   end if
   call write_info('')
  end if

  ok = .false.
  code       = 0
  n_ref      = 0
  n_ref_2    = 0
  n_ref_3    = 0
  n_ref_neg  = 0
  cosdir     = 0.
  HKL_file%HKL = ''

  IF(HKL_file%SHELX) THEN     ! lecture d'un fichier.HKL
   !if(hkl_format_free) then ! verifie qu'il y a au moins 5 colonnes (format libre)
    read(HKL_unit, '(a)', iostat=ier) read_line
    if(ier /=0) then
     call write_info('')
     call write_info(' > Problem reading HKL file in free format !')
     call write_info('')
     return
    endif
    call nombre_de_colonnes(read_line, nb_col)
    if(hkl_format_free .and. nb_col < 5) then
     call write_info('')
     call write_info(' > Problem reading HKL file in free format : at least 5 columns have to be present !')
     call write_info('')
     return
    else
     if(u_case(hkl_format(1:11)) == '(3I4,2F8.2)') then
      if(nb_col == 6) then
       hkl_format = '(3I4,2F8.4,I4)'
      elseif(nb_col == 12) then
       hkl_format = hkl_SHELX_format
      end if
     end if
    endif
    rewind(unit=HKL_unit)



   if(debug_proc%level_2) then
    if(hkl_format_free) then
     call write_debug_proc_level(2, "READ (in free format) and CHECK HKL")
    else
     !call write_debug_proc_level(2, "READ (format="//trim(hkl_format)//") and CHECK HKL !!")
    endif
   end if

   do
    if(hkl_format_free) then
     read(HKL_unit, *, iostat=ier) h_, k_, l_, F2_, sig_F2_
     if(ier >0) then
      write(nref_string, '(I8)') n_ref+1
      write(message_text, '(5x,3a)') ' Error reading line ', trim(adjustl(nref_string)), ' ! Program will stop reading file.'
      call write_info(trim(message_text))
      return
     end if
     if(ier < 0) exit  ! fin du fichier
      code_ = 0
     !write(31,'(3I4,2F8.2,I4,6F8.5)', iostat=ier) h_,k_,l_, F2_, sig_F2_
    else
     !READ(HKL_unit,'(3I4,2F8.2,I4,6F8.5)', iostat=ier) h_,k_,l_, F2_, sig_F2_,code_,cosdir(1:6)
     if(nb_col == 5) then
      READ(HKL_unit, fmt=trim(hkl_format), iostat=ier) h_,k_,l_, F2_, sig_F2_
      if(ier >0) then
       write(nref_string, '(I8)') n_ref+1
       write(message_text, '(5x,3a)') ' Error reading line ', trim(adjustl(nref_string)), ' ! Program will stop reading file.'
       call write_info(trim(message_text))
       return
      end if
      if (ier < 0) exit ! fin du fichier atteinte
      code_ = 0
     elseif(nb_col == 6) then
      READ(HKL_unit, fmt=trim(hkl_format), iostat=ier) h_,k_,l_, F2_, sig_F2_, code_
      if(ier >0) then
       write(nref_string, '(I8)') n_ref+1
       write(message_text, '(5x,3a)') ' Error reading line ', trim(adjustl(nref_string)), ' ! Program will stop reading file.'
       call write_info(trim(message_text))
       return
      end if
      if (ier < 0) exit ! fin du fichier atteinte
     else
      READ(HKL_unit, fmt=trim(hkl_format), iostat=ier) h_,k_,l_, F2_, sig_F2_,code_,cosdir(1:6)
      if (ier < 0) exit ! fin du fichier atteinte
     end if
    endif
    IF(h_==0 .AND. k_==0 .AND. l_==0) exit
    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)

    if (n_ref ==1) then
     if (cosdir(1) < eps .and. cosdir(2) < eps .and. cosdir(3) < eps .and. cosdir(4) < eps .and. &
         cosdir(5) < eps .and. cosdir(6) < eps) then
      cos_exist = .false.
     else
      cos_exist = .true.
     end if
    endif
   end do
   ok = .true.
   HKL_data_known = .true.
   !HKL_file%HKL   = HKL_file%NAME

   div_10 = 0
   max_F2 = MAXVAL(F2(1:n_ref))
   do
    if (max_F2 > 99999.99) then
      max_F2 = max_F2 / 10.
     div_10 = div_10 + 1
    else
     exit
    end if
   END do

   if(div_10 /=0) then
    i1 = index(HKL_file%NAME, '.')
    if(i1 /=0) then
     HKL_file%HKL =  HKL_file%NAME(1:i1-1)//'_SX.hkl'
    else
     HKL_file%HKL =  '_SX.hkl'
    endif
    open (UNIT=31, FILE=TRIM(HKL_file%HKL))


    do i=1, n_ref
     F2(i) = F2(i) / 10.**div_10
     sig_F2(i) = sig_F2(i) / 10.**div_10
     !write(31, '(3I4, 2F8.2,I4, 6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i), code(i), cos_dir(i, 1:6)
     write(31, '(3I4, 2F8.2,I4, 6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i), code(i), cos_dir(i, 1), &
                                         cos_dir(i, 2),cos_dir(i, 3),cos_dir(i, 4),cos_dir(i, 5),cos_dir(i, 6)
    END do
   else
    HKL_file%HKL   = HKL_file%NAME
   end if

  ELSEIF(HKL_file%CIF) then ! lecture fichier import.CIF

   call get_wavelength_from_CIF_file(HKL_unit, wave)
   keyword_WAVE = .true.
   call get_cell_parameters_from_CIF_file(HKL_unit, unit_cell%param, unit_cell%param_esd )
   call get_UB_matrix_from_CIF_file(HKL_unit, UB_matrix)
   call volume_calculation('no_out')
   !call volume_calculation('out')
   call get_crystal_system_from_CIF_file(HKL_unit, crystal_system)
   keyword_CELL = .true.
   call get_CIF_parameters_from_import_CIF(HKL_unit)
   call get_H_M_from_CIF_file(HKL_unit, unit_cell%H_M)
   call read_lines_before_hkl(HKL_unit, cos_exist)
   IF(.not. cos_exist) then
    HKL_file%HKL = 'import.hkl'
   else
    HKL_file%HKL = 'import_cos.hkl'
   endif
   open (UNIT=31, FILE=TRIM(HKL_file%HKL))

   if(debug_proc%level_2) then
    if(cos_exist) then
     call write_debug_proc_level(2, "READ (with dir. cos.) and CHECK HKL")
    else
     call write_debug_proc_level(2, "READ (without dir. cos.) and CHECK HKL")
    endif
   endif
   do
    if (.not. cos_exist) then
     READ(HKL_unit, '(i4,2i5,2F9.2,I5)', IOSTAT=ier) h_,k_,l_,F2_, sig_F2_, code_  ! sans les cos. dir.
     if (ier < 0 ) exit
     IF(h_==0 .AND. k_==0 .AND. l_==0) exit
     WRITE(31, '(3i4,2F8.2,I4)')h_,k_,l_,F2_, sig_F2_, code_
    else
     READ(HKL_unit,'(i4,2i5,2F9.2,I5)', IOSTAT=ier) h_,k_,l_,F2_, sig_F2_, code_
     !write(*,'(i4,2i5,2F9.2,I5)', IOSTAT=ier) h_,k_,l_,F2_, sig_F2_, code_

      if (ier /=0 ) exit                                                         ! avec les cos. dir.
      IF(h_==0 .AND. k_==0 .AND. l_==0) exit
     READ(HKL_unit,'(F8.5,5F9.5)',IOSTAT=ier) cosdir(1:6)
    !READ(HKL_unit,'(6F9.5)',IOSTAT=ier) cosdir(1:6)

    if (ier /=0 ) exit
     write(31,'(3I4,2F8.2,I4,6F8.5)', iostat=ier) h_,k_,l_, F2_, sig_F2_,code_, cosdir(1:6)
    endif

    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)
   END do

   ok = .true.
   HKL_data_known = .true.


  ELSEIF(HKL_file%FINAL_y) then
   call get_wavelength_from_FINAL_Y_file(HKL_unit, wave)
   keyword_WAVE = .true.
   call get_cell_parameters_from_FINAL_Y_file(HKL_unit, unit_cell%param)
   call volume_calculation('no_out')
   keyword_CELL = .true.

   call get_modulation_vect(HKL_unit, qvec_type, nvk)

   if(.not. QVEC_type ) then
    HKL_file%HKL = 'final_y.HKL'
   else
    HKL_file%HKL = 'final_y_mod.HKL'
   endif
   OPEN (UNIT=31, FILE=TRIM(HKL_file%HKL))

   if(.not. QVEC_type) then
    do
     call read_reflexion_from_FINAL_Y(HKL_unit, h_, k_, l_, F2_, sig_F2_, HKL_ok)
     IF(.NOT. HKL_ok) exit
     IF(h_==0 .AND. k_==0 .AND. l_==0) exit

     call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)
     WRITE(31, '(3i4,2F8.2)') h_,k_,l_,F2_, sig_F2_
    END do
   else
    do
     call read_mod_reflexion_from_FINAL_Y(HKL_unit, h_, k_, l_, m_, F2_, sig_F2_, nvk, HKL_ok)
     IF(.NOT. HKL_ok) exit
     !IF(h_==0 .AND. k_==0 .AND. l_==0) exit

     call check_HKLM(h_, k_, l_, m_, F2_, sig_F2_, code_, cosdir)
     WRITE(31, '(4i4,2F8.2)') h_,k_,l_,m_, F2_, sig_F2_
    END do

   endif
   ok = .true.
   HKL_data_known = .true.


  ELSEIF(HKL_file%RAW) then

   i1 = INDEX(HKL_file%NAME, ".")
   if (i1/=0) then
    HKL_file%HKL = TRIM(HKL_file%NAME(1:i1-1))//'_raw.hkl'
   else
    HKL_file%HKL = TRIM(HKL_file%NAME(1:))//'_raw.hkl'
   endif
   OPEN (UNIT=31, FILE=TRIM(HKL_file%HKL))

   if(debug_proc%level_2) call write_debug_proc_level(2, "READ .RAW file (with dir. cos.) and CHECK HKL")
   do
    READ(HKL_unit, '(a)', IOSTAT=ier) read_line
    IF(ier/=0) exit
    READ(read_line,      '(3I4)')       h_, k_, l_
    READ(read_line(13:), *)             F2_, sig_F2_
    READ(read_line(29:), '(I4,6F8.5)')  code_, cosdir(1:6)

    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)
   END do

   ok = .true.
   HKL_data_known = .true.

   div_10 = 0
   max_F2 = MAXVAL(F2(1:n_ref))
   do
    if (max_F2 > 99999.99) then
     div_10 = div_10 + 1
    else
     exit
    end if
   END do

   do i=1, n_ref
    if(div_10 /=0) then
     F2(i)     = F2(i)     / 10**div_10
     sig_F2(i) = sig_F2(i) / 10**div_10
    end if

    !write(31, '(3I4, 2F8.2,I4, 6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i), code(i), cos_dir(i, 1:6)
     write(31, '(3I4, 2F8.2,I4, 6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i), code(i), cos_dir(i, 1), &
                                         cos_dir(i, 2),cos_dir(i, 3),cos_dir(i, 4),cos_dir(i, 5),cos_dir(i, 6)

   END do


  ELSEIF(HKL_file%M91) then
   i1 = INDEX(HKL_file%NAME, ".")
   if (i1/=0) then
    HKL_file%HKL = TRIM(HKL_file%NAME(1:i1-1))//'_91.hkl'
   else
    HKL_file%HKL = TRIM(HKL_file%NAME(1:))//'_m91.hkl'
   endif
   OPEN (UNIT=31, FILE=TRIM(HKL_file%HKL))

   if(debug_proc%level_2) call write_debug_proc_level(2, "READ .M91 and CHECK HKL")
   do
    READ(HKL_unit, '(a)', IOSTAT=ier) read_line
    IF(ier/=0) exit
    READ(read_line,      '(3I4,2F9.1,a)') h_, k_, l_, F2_, sig_F2_, M91_string
    IF(h_ == 999) exit
    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)
    HKL_line(n_ref) = M91_string
   END do
   ok = .true.
   HKL_data_known = .true.

   do i=1, n_ref
    write( 31, '(3I4, 2F8.2)') h(i),k(i),l(i), F2(i)/10., sig_F2(i)/10.
   end do

  ELSEIF(HKL_file%M95) then
   i1 = INDEX(HKL_file%NAME, ".")
   if (i1/=0) then
    HKL_file%HKL = TRIM(HKL_file%NAME(1:i1-1))//'_m95.hkl'
   else
    HKL_file%HKL = TRIM(HKL_file%NAME(1:))//'_m95.hkl'
   endif
   OPEN (UNIT=31, FILE=TRIM(HKL_file%HKL))


   READ(HKL_unit, '(a)', IOSTAT=ier) read_line
   IF(ier/=0) return
   M95_comment_line_nb = 0
   IF(read_line(1:19) == 'import_report_begin') then
    do i=1, 8
     READ(HKL_unit, '(a)', IOSTAT=ier) read_line
     IF(ier/=0) return
     M95_comment_line(i) = read_line
    end do
    M95_comment_line_nb = 9
   else
    BACKSPACE(UNIT=HKL_unit)
   endif

   if(debug_proc%level_2) call write_debug_proc_level(2, "READ .M95 and CHECK HKL")
   do
    READ(HKL_unit, '(a)', IOSTAT=ier) read_line
    IF(ier/=0) exit
    READ(read_line,      '(I6, 3I4, 4F7.2, 2E15.6, F10.3, 2I2)') M95_i(1), h_, k_, l_, M95_r(1:4), &
                                                                 F2_, sig_F2_, M95_r(5), M95_i(2:3)
    READ(HKL_unit, '(a)', IOSTAT=ier) read_line
    IF(ier/=0) exit

    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)
    HKL_line(n_ref) = read_line
    cos_dir(n_ref,1:4)  = M95_r(1:4)
    cos_dir(n_ref,5)    = M95_r(5)
    HKL_flag(n_ref,1:3) = M95_i(1:3)
   END do
   ok = .true.
   HKL_data_known = .true.

   do i=1, n_ref
    write(31, '(3I4, 2F8.2)') h(i),k(i),l(i), F2(i)/10., sig_F2(i)/10.
   end do

  ELSEIF(HKL_file%INT) then



   i1 = INDEX(HKL_file%NAME, ".")
   if (i1/=0) then
    HKL_file%HKL = TRIM(HKL_file%NAME(1:i1-1))//'_int.hkl'
   else
    HKL_file%HKL = TRIM(HKL_file%NAME(1:))//'_int.hkl'
   endif
   OPEN (UNIT=31, FILE=TRIM(HKL_file%HKL))



   READ(HKL_unit, '(a)') INT_string
   READ(HKl_unit, '(a)') format_
   READ(HKL_unit, *) WAVE
   if(debug_proc%level_2) call write_debug_proc_level(2, "READ and CHECK HKL")
   do
    READ(HKL_unit, '(a)', IOSTAT=ier) read_line
    IF(ier/=0) exit
    READ(read_line,    fmt = trim(format_)) h_, k_, l_, F2_, sig_F2_, code_, INT_real(1), INT_real(2), INT_real(3), INT_real(4)

    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)

   END do
   ok = .true.
   HKL_data_known = .true.

   do i=1, n_ref
    write( 31, '(3I4, 2F8.2)') h(i),k(i),l(i), F2(i), sig_F2(i)
   end do

  ELSEIF(HKL_file%COL) then

   i1 = INDEX(HKL_file%NAME, ".")
   if (i1/=0) then
    HKL_file%HKL = TRIM(HKL_file%NAME(1:i1-1))//'_col.hkl'
   else
    HKL_file%HKL = TRIM(HKL_file%NAME(1:))//'_col.hkl'
   endif
   OPEN (UNIT=31, FILE=TRIM(HKL_file%HKL))

   if(debug_proc%level_2) call write_debug_proc_level(2, "READ .COL and CHECK HKL")
   do
    READ(HKL_unit, '(a)', IOSTAT=ier) read_line
    IF(ier/=0) exit
    READ(read_line, '(i6,3i4,2F10.2,4F8.2)') numor, h_, k_, l_, F2_, sig_F2_, INT_real(1), INT_real(2), INT_real(3), INT_real(4)
    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)
   END do
   ok = .true.
   HKL_data_known = .true.

   do i=1, n_ref
    write( 31, '(3I4, 2F8.2)') h(i),k(i),l(i), F2(i), sig_F2(i)
   end do

  else
   call write_info('')
   call write_info('  >>>> Unknown extension for data file <<<<')
   call write_info('')
   CLOSE(UNIT = HKL_unit)
   CLOSE(unit=31)
   ok = .false.
   return
  endif
  CLOSE(UNIT=HKL_unit)
  CLOSE(unit=31)


   if(len_trim(input_string) /=6) then
   if(WRITE_details) then
   call write_info('')
   write(message_text,'(a,I8)')  '  >> Total number of reflections:          ', n_ref
   call write_info(TRIM(message_text))
   end if
   IF(n_ref == 0) then
    ok = .false.
    return
   endif

   if(WRITE_details)  then
   write(message_text,'(a,I8,2x,a,F6.2,a)')  '  >> Number of reflections with I > 2sig.: ', n_ref_2,  &
                                             '(', (real(n_ref_2)/real(n_ref))*100., '%)'
   call write_info(TRIM(message_text))
   write(message_text,'(a,I8,2x,a,F6.2,a)')  '  >> Number of reflections with I > 3sig.: ', n_ref_3,  &
                                             '(', (real(n_ref_3)/real(n_ref))*100., '%)'
   call write_info(TRIM(message_text))
   write(message_text,'(a,I8,2x,a,F6.2,a)')  '  >> Number of reflections with I < 0.   : ', n_ref_neg, &
                                             '(', (real(n_ref_neg)/real(n_ref))*100., '%)'
   call write_info(TRIM(message_text))
   call write_info('')

   WRITE(message_text,'(a,I4,2x,I4)')       '                           . h_min, h_max: ', MINVAL(h(1:n_ref)), &
                                                                          MAXVAL(h(1:n_ref))
   call write_info(TRIM(message_text))
   WRITE(message_text,'(a,I4,2x,I4)')       '                           . k_min, k_max: ', MINVAL(k(1:n_ref)), &
                                                                          MAXVAL(k(1:n_ref))
   call write_info(TRIM(message_text))
   WRITE(message_text,'(a,I4,2x,I4)')       '                           . l_min, l_max: ', MINVAL(l(1:n_ref)), &
                                                                          MAXVAL(l(1:n_ref))
   call write_info(TRIM(message_text))
   call write_info('')

   IF(MAXVAL(code(1:n_ref)) /=0) then
    WRITE(message_text,'(a,I3)')            '                           . number of scans: ', MAXVAL(code(1:n_ref))
    call write_info(TRIM(message_text))
    call write_info('')
   endif

   ! <i>, <sig>
   write(message_text, '(a,F15.2)')   '                           . Max. of intensity: ', MAXVAL(F2(1:n_ref))
   call write_info(trim(message_text))
   write(message_text, '(a,F15.2)')   '                           . Min. of intensity: ', MINVAL(F2(1:n_ref))
   call write_info(trim(message_text))
   write(message_text, '(a,F15.2)')   '                           . <F2>             : ', SUM(F2(1:n_ref))/n_ref
   call write_info(trim(message_text))
   write(message_text, '(a,F15.2)')   '                           . <sig>            : ', SUM(SIG_F2(1:n_ref))/n_ref
   call write_info(trim(message_text))
   write(message_text, '(a,F15.2)')   '                           . <F2/sig>         : ', SUM(F2(1:n_ref)/sig_F2(1:n_ref))/n_ref
   call write_info(trim(message_text))
   write(message_text, '(a,F15.2)')   '                           . <F2>/<sig>       : ', SUM(F2(1:n_ref))/SUM(sig_F2(1:n_ref))
   call write_info(trim(message_text))
   call write_info('')

   !if(MAXVAL(F2(1:n_ref)) > 99999.99) then
   ! hkl_format = '(3I4,2F10.2,I4,6F8.5)'
   !elseif(MAXVAL(F2(1:n_ref)) > 9999999.99) then
   ! hkl_format = '(3I4,2F12.2,I4,6F8.5)'
   !elseif(MAXVAL(F2(1:n_ref)) > 999999999.99) then
   ! hkl_format = '(3I4,2F12.2,I4,6F8.5)'
   !end if


    IF(keyword_CELL) then
    WRITE(message_text,'(a,F6.3,2x,F6.3)')  '                           . STL_min, STL_max (A-1)     :  ',  &
                                                                          MINVAL(sintheta_lambda(1:n_ref)), &
                                                                          MAXVAL(sintheta_lambda(1:n_ref))
    call write_info(TRIM(message_text))
    WRITE(message_text,'(a,F6.3,2x,F6.3)')  '                           . d_min, d_max (A)           :  ', &
                                                                          MINVAL(d_hkl(1:n_ref)),          &
                                                                          MAXVAL(d_hkl(1:n_ref))
    call write_info(TRIM(message_text))
    IF(keyword_WAVE) then
     WRITE(message_text,'(a,F6.2,2x,F6.2)') '                           . theta_min, theta_max (deg.):  ', &
                                                                          MINVAL(theta_hkl(1:n_ref)),      &
                                                                          MAXVAL(theta_hkl(1:n_ref))
     call write_info(TRIM(message_text))
     call write_info('')
    endif
   endif
   end if
   end if ! fin de la condition "if len_trim(input_string) /=6"


   !IF(LEN_TRIM(HKL_file%HKL) /=0 .and. .not. hkl_format_SHELX) then
   IF(.not. HKL_file%SHELX .or. div_10 /=0) then
    call write_info(' ')
    write (message_text, '(3a)') '   ==> ', trim(HKL_file%HKL), ' FILE has been created (SHELX format)'
    call write_info(TRIM(message_text))
    !IF(HKL_file%RAW .AND. div_10 /=0) then
    !IF(div_10 /=0) then
    ! write (message_text, '(a,I2)') ' >> F2 and sigmas have been divided by 10**',div_10
    ! call write_info(TRIM(message_text))
    !endif
    call write_info(' ')
   endif

   IF(div_10 /=0) then

    write(out_string, '(I6)') 10*div_10
    write(message_text, '(2x,3a)') ' >> F2 and sigmas have been divided by ', trim(adjustl(out_string)),'.'
    call write_info(TRIM(message_text))
    close(unit=31)
   endif
   call write_info(' ')



   IF(keyword_create_CIF) then
    call write_CIF_file('DATA_LIMITS')
    call write_CIF_file('HKL_DATA')
   endif


 return
end subroutine read_data_file

!----------------------------------------------------------------------------------------------

 subroutine READ_FCF_file()
  USE cryscalc_module, only : FCF_file_name, FCF_plot, FCF_plot_stl, HKL_unit, file_out, file_out_n,    &
                              cell_star, cos_angle_star,                                                &
                              PGF_file, PGF_data, allocate_PGF_data_arrays, deallocate_PGF_data_arrays, &
                              winplotr_exe, message_text, lecture_ok, debug_proc
  USE hkl_module,      only : h,k,l, F2, F2c, sig_F2, n_ref, n_ref_eff, sinTheta_lambda
  USE macros_module,   only : multiple
  USE IO_module,       ONLY : write_info
  implicit none
   character (len=256)           :: read_line
   character (len=2)             :: hkl_status
   integer                       :: i, i1, i_error, long
   real                          :: Fo, Fc, num_R1, den_R1, R1, weight, num_wR2, den_wR2, wR2
   real                          :: Q_hkl, d_hkl, stl_hkl
   real, dimension(3)            :: ref_H

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_FCF_FILE")

  call write_info(' ')
  call write_info('  >> FCF file name: '// TRIM(FCF_file_name))
  call write_info(' ')

  if(FCF_plot_stl) then
   call READ_CIF_input_file(TRIM(FCF_file_name), 'out')
  end if



  close (unit=HKL_unit)
  open(unit=HKL_unit, file=trim(FCF_file_name))
   do
    read(unit=HKL_unit, fmt='(a)', iostat=i_error) read_line
    if(i_error /=0) then
     lecture_ok = .false.
     exit
    end if
    read_line = adjustl(read_line)
    long = len_trim('_refln_observed_status')
    if(read_line(1:long) /= '_refln_observed_status') cycle
    lecture_ok = .true.
    exit
   end do

   if(lecture_ok) then
    n_ref     = 0
    n_ref_eff = 0
    R1        = 0.
    wR2       = 0.
    num_R1    = 0.
    den_R1    = 0.
    num_wR2   = 0.
    den_wR2   = 0.


    do
     read(unit=HKL_unit, fmt='(a)', iostat=i_error) read_line
     if(i_error /=0 .or. len_trim(read_line) == 0) exit
     n_ref = n_ref + 1
     read(read_line,*, iostat=i_error) h(n_ref), k(n_ref), l(n_ref), F2c(n_ref), F2(n_ref), sig_F2(n_ref), hkl_status
     if(file_out .and. multiple(n_ref,file_out_n)) then
      write(message_text,'(2x,I8,2x,3I4, 3F15.2)') n_ref, h(n_ref), k(n_ref), l(n_ref), F2c(n_ref), F2(n_ref), sig_F2(n_ref)
      call write_info(trim(message_text))
     end if
     if(F2(n_ref) > 0.) then
      Fo = sqrt(F2(n_ref))
     else
      Fo = 0.
     end if
     Fc = sqrt(F2c(n_ref))
     num_R1 = num_R1 + abs(Fo - Fc)
     den_R1 = den_R1 + Fo
     !if(F2(n_ref) > 0.) then

     weight = 1./sig_F2(n_ref)**2.
     num_wR2 = num_wR2 + weight*(F2(n_ref) - F2c(n_ref))**2.
     den_wR2 = den_wR2 + weight*F2(n_ref)**2.
     !end if

     if(FCF_plot_stl) then
      if(n_ref == 1)  call calcul_cell_star()
      ref_H(1) = h(n_ref)
      ref_H(2) = k(n_ref)
      ref_H(3) = l(n_ref)
      call Calcul_Q_d_stl(ref_H, Q_hkl, d_hkl, stl_hkl)
      sinTheta_lambda(n_ref) = stl_hkl

      !Q_hkl = h(n_ref)**2 * a_star**2 + k(n_ref)**2 * b_star**2 + l(n_ref)**2 * c_star**2            &
      !       + 2*h(n_ref)*k(n_ref) * a_star*b_star*cos_gama_star                                      &
      !       + 2*k(n_ref)*l(n_ref) * b_star*c_star*cos_alfa_star                                      &
      !       + 2*l(n_ref)*h(n_ref) * c_star*a_star*cos_beta_star

      !d_hkl = 1/sqrt(Q_hkl)
      !stl_hkl = 1/(2*d_hkl)

     end if
    end do
   end if

   R1 = num_R1 / den_R1
   wR2 = sqrt(num_wR2 / den_wR2)

   if(file_out) call write_info('')
   write(message_text, '(5x,a,F10.4)') ' R1 = ', R1
   call write_info(trim(message_text))
   write(message_text, '(5x,a,F10.4)') 'wR2 = ', wR2
   call write_info(trim(message_text))
  !write(*,*) '  R1 = ', R1
  !write(*,*) ' wR2 = ', wR2


  close(unit=HKL_unit)
  if(FCF_plot) then
   i1 =index(FCF_file_name, '.')
   call Allocate_PGF_data_arrays(n_ref)
   if(.not. FCF_plot_stl) then
    do i=1, n_ref
     PGF_data%X(i) = F2(i)
     PGF_data%Y(i) = F2c(i)
     WRITE(pgf_data%string(i), '(a,3I4,a)') '(', h(i), k(i), l(i), ')'
    end do
    pgf_file%name = FCF_file_name(1:i1-1)//'_FCF.pgf'
    call create_PGF_file(TRIM(PGF_file%name), PGF_data%X, PGF_data%Y, PGF_data%string, n_ref, "fcf")
   else
    do i=1, n_ref
     !PGF_data%X(i) = stl(i)
     PGF_data%X(i) = sinTheta_lambda(i)
     PGF_data%Y(i) = F2c(i) - F2(i)
     WRITE(pgf_data%string(i), '(a,3I4,a)') '(', h(i), k(i), l(i), ')'
    end do

    pgf_file%name = FCF_file_name(1:i1-1)//'_FCF_stl.pgf'
    call create_PGF_file(TRIM(PGF_file%name), PGF_data%X, PGF_data%Y, PGF_data%string, n_ref, "fcf_stl")
   end if
   call deallocate_PGF_data_arrays

   IF(LEN_TRIM(winplotr_exe) == 0) then
    call write_info('')
    call write_info(' Warning: WinPLOTR is not installed ')
    call write_info('')
   else
    !call system('winplotr '//trim(pgf_file%name), .true. )    ! lf95
    call system('winplotr '//trim(pgf_file%name))              ! g95
   endif
  end if


  return
 end subroutine READ_FCF_file

!----------------------------------------------------------------------------------------------
subroutine check_HKL(H_h, H_k, H_l, H_F2, H_sig, code_, H_cos)
 USE HKL_module
 USE cryscalc_module, ONLY : keyword_cell, keyword_WAVE, known_theta,  unit_cell, wave => wavelength, &
                             file_out, file_out_n, hkl_format_free, hkl_format, message_text, debug_proc
 USE macros_module,   ONLY : multiple
 USE IO_module,       ONLY : write_info

 Implicit none
 INTEGER, INTENT(IN)            :: H_h, H_k, H_l
 REAL,    INTENT(INOUT)         :: H_F2, H_sig
 INTEGER, INTENT(IN)            :: code_
 REAL, DIMENSION(6), INTENT(IN) :: H_cos
 REAL                           :: Q, Z
 real, parameter                :: eps = 0.000001
 real, parameter                :: pi=3.1415926535897932
 integer                        :: ier

 !if(debug_proc%level_3)  call write_debug_proc_level(3, "CHECK_HKL")


    n_ref = n_ref + 1
    call check_nref(n_ref)
    !IF(H_F2 < 0.00001) n_ref_neg = n_ref_neg + 1
    !IF(.not. HKL_file%read_NEG) call check_F2(H_F2, H_sig)
    if(H_F2 < 0.00001) then
     IF(.not. HKL_file%read_NEG) then
      call check_F2(H_F2, H_sig)
     else
      n_ref_neg = n_ref_neg + 1
     end if
    end if

    if(H_F2 >= 3*H_sig) n_ref_3 = n_ref_3 + 1
    if(H_F2 >= 2*H_sig) n_ref_2 = n_ref_2 + 1

    IF(keyword_CELL) then
     call calcul_Q(real(H_h), real(H_k), real(H_l), unit_cell%param, Q)

     sintheta_lambda(n_ref) = SQRT(Q)/2.
     d_hkl(n_ref)           = 1/(2*  sintheta_lambda(n_ref))

     IF(keyword_WAVE) then
      IF (d_hkl(n_ref) > eps) then
       Z =wave / (2. * d_hkl(n_ref))
      else
       return
      endif

      IF (Z**2 < 1.) then
       theta_hkl(n_ref) = ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
      else
       theta_hkl(n_ref) = 90.
      endif
      known_theta = .true.


     endif
    endif

    h(n_ref)               = H_h
    k(n_ref)               = H_k
    l(n_ref)               = H_l
    F2(n_ref)              = H_F2
    if(H_sig < 0.01) then
     write(message_text, '(5x,a,3I4,a,F10.2,a)') '!! (',H_h, H_k, H_l, ') : strange sigma (', H_sig, ' fixed to 0.01) !!'
     call write_info(trim(message_text))
     sig_F2(n_ref) = 0.01
    else
     sig_F2(n_ref)          = H_sig
    end if
    I_sigma(n_ref)         = H_F2/ H_sig
    I_(n_ref)              = H_F2
    code(n_ref)            = code_
    cos_dir(n_ref,1:6)     = H_cos(1:6)
    !som_F2                 = som_F2 + F2(n_ref)
    !som_sig                = som_sig + sig_F2(n_ref)
    !som_F2_sig             = som_F2_sig + I_sigma(n_ref)

    WRITE(HKL_string(n_ref), '(a,3I4,a,F15.5)') '(', H_h, H_k, H_l, ')   F2= ', F2(n_ref)

    !if(file_out .and. file_out_n*int(n_ref/file_out_n) == n_ref ) then
    if(file_out .and. multiple(n_ref,file_out_n)) then
     if(hkl_format_free) then
      write(message_text,'(I8,2x,3I4,2F15.2)', iostat=ier) n_ref, H_h,H_k,H_l, H_F2, H_sig
     else
      write(message_text,'(I8,2x,3I4,2F15.5,I4)', iostat=ier) n_ref, H_h,H_k,H_l, H_F2, H_sig,code_
      !long = len_trim(hkl_format)
      !write(fmt_, '(23a)') '(I8,2x', trim(hkl_format(2:long)),')'
      !write(message_text, trim(fmt_)) n_ref, H_h,H_k,H_l, H_F2, H_sig
     endif
     call write_info(trim(message_text))
    endif

  return
end subroutine check_HKL

!----------------------------------------------------------------------------------------------
subroutine check_HKLM(H_h, H_k, H_l, H_m, H_F2, H_sig, code_, H_cos)
 USE HKL_module
 USE cryscalc_module, ONLY : keyword_cell, keyword_WAVE, known_theta,  unit_cell, wave => wavelength, Qvec, &
                             file_out, file_out_n, message_text, debug_proc
 USE IO_module,       ONLY : write_info
 USE macros_module,   ONLY : multiple

 Implicit none
 INTEGER, INTENT(IN)            :: H_h, H_k, H_l, H_m
 REAL,    INTENT(INOUT)         :: H_F2, H_sig
 INTEGER, INTENT(IN)            :: code_
 REAL, DIMENSION(6), INTENT(IN) :: H_cos
 REAL                           :: Q, Z
 real, parameter                :: eps = 0.000001
 real, parameter                :: pi=3.1415926535897932
 integer                        :: ier

 if(debug_proc%level_2)  call write_debug_proc_level(2, "CHECK_HKLM")

    n_ref = n_ref + 1
    call check_nref(n_ref)
    IF(H_F2 < 0.00001) n_ref_neg = n_ref_neg + 1
    IF(.not. HKL_file%read_NEG) call check_F2(H_F2, H_sig)

    if(H_F2 >= 3*H_sig) n_ref_3 = n_ref_3 + 1
    if(H_F2 >= 2*H_sig) n_ref_2 = n_ref_2 + 1

    IF(keyword_CELL) then
     call calcul_Q(real(H_h)+H_m*Qvec(1),real(H_k)+H_m*Qvec(2),real(H_l)+H_m*QVEC(3), unit_cell%param, Q)

     sintheta_lambda(n_ref) = SQRT(Q)/2.
     d_hkl(n_ref)           = 1/(2*  sintheta_lambda(n_ref))

     IF(keyword_WAVE) then
      IF (d_hkl(n_ref) > eps) then
       Z =wave / (2. * d_hkl(n_ref))
      else
       return
      endif

      IF (Z**2 < 1.) then
       theta_hkl(n_ref) = ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
      else
       theta_hkl(n_ref) = 90.
      endif
      known_theta = .true.
     endif
    endif

    h(n_ref)               = H_h
    k(n_ref)               = H_k
    l(n_ref)               = H_l
    F2(n_ref)              = H_F2
    if(H_sig < 0.01) then
     write(message_text, '(5x,a,3I4,a,F10.2,a)') '!! (',H_h, H_k, H_l, ') : strange sigma (', H_sig, ' fixed to 0.01) !!'
     call write_info(trim(message_text))
     sig_F2(n_ref) = 0.01
    else
     sig_F2(n_ref)          = H_sig
    endif
    I_sigma(n_ref)         = H_F2/ H_sig
    I_(n_ref)              = H_F2
    code(n_ref)            = code_
    cos_dir(n_ref,1:6)     = H_cos(1:6)

    WRITE(HKL_string(n_ref), '(a,4I4.0,a)') '(', H_h, H_k, H_l, H_m, ')'

   if(file_out .and. multiple(n_ref,file_out_n)) then
    write(message_text,'(I8,2x,4I4,2F15.2)', iostat=ier) n_ref, H_h,H_k,H_l, H_m, H_F2, H_sig
    call write_info(trim(message_text))
   endif

  return
end subroutine check_HKLM

!----------------------------------------------------------------------------------------------------------------
subroutine check_F2(F2_, sig_F2_)
 implicit none
 REAL, INTENT(INOUT) :: F2_, sig_F2_


  if(F2_     < 0.0001)  F2_= 0.0001
  if(sig_F2_ < 0.00001) sig_F2_=sqrt(abs(F2_))

 return
END subroutine check_F2

!----------------------------------------------------------------------------------------------------------------
subroutine check_nref(n_ref)
 !USE HKL_module, ONLY : max_ref
 USE cryscalc_module, ONLY : max_ref
 USE IO_module,       ONLY : write_info
 Implicit none
 INTEGER, INTENT(IN)  :: n_ref

    IF(n_ref > max_ref) then
     call write_info('')
     call write_info('  >> Number of reflections larger than arrays dimensions.')
     call write_info('     Change the hkl array in the cryscalc.ini setting file')
     call write_info('     if lower than 500000 or contact the author of CRYSCALC')
     call write_info('     program (T.R./ CDIFX Rennes):')
     call write_info('      thierry.roisnel@univ-rennes1.fr')
     call write_info('     Program will be stopped.')
     call write_info(' ')
     stop
    endif
  return
end subroutine check_nref

!----------------------------------------------------------------------------------------------------------------
subroutine classement(n, input_string)
 USE cryscalc_module,   ONLY       : cut_off, unit_cell, wave => wavelength, sort_plot, shell_plot,    &
                                     sort_out, sort_out_n, write_details,                              &
                                     keyword_CELL, keyword_WAVE, message_text,                         &
                                     nb_shell, shell_arg_min_max, shell_min, shell_max, shell_out,     &
                                     PGF_data, PGF_file, known_theta, winplotr_exe, input_line,        &
                                     allocate_PGF_data_arrays, hkl_statistics, debug_proc
 USE CFML_Math_General, ONLY       : sort
 USE hkl_module
 USE IO_module

 implicit none
 INTEGER,           INTENT(IN)     :: n            ! numero du passage dans la routine
 character (len=*), intent(in)     :: input_string
 INTEGER                           :: long_input_string
 integer                           :: i, ier
 real                              :: max_stl,   expected_max_stl
 real                              :: min_stl,   expected_min_stl
 REAL                              :: max_d_hkl, expected_max_d_hkl
 REAL                              :: min_d_hkl, expected_min_d_hkl
 REAL                              :: max_theta, expected_max_theta
 REAL                              :: min_theta, expected_min_theta
 real                              :: min_cut_off, expected_min_cut_off
 real                              :: max_cut_off, expected_max_cut_off
 real                              :: intensity_min, expected_min_intensity
 real                              :: intensity_max, expected_max_intensity
 real                              :: val_1, val_2
 real, parameter                   :: eps = 0.000001
 INTEGER,  ALLOCATABLE, DIMENSION(:) :: ordered_array

 LOGICAL                           :: sort_X
 LOGICAL                           :: X_stl, X_d, X_theta, X_int, X_isig
 LOGICAL                           :: shell_X, shell_stl, shell_d, shell_theta, shell_int, shell_isig
 LOGICAL                           :: input_sort, input_no_sort, ok

 if(debug_proc%level_2)  call write_debug_proc_level(2, "CLASSEMENT ("//trim(input_string)//")")

  long_input_string = LEN_TRIM(input_string)
  input_sort    = .false.
  input_no_sort = .false.
  if(long_input_string == 5) then
   if(input_string(1:5) == 'sort_') input_sort = .true.
  elseif(long_input_string == 7) then
   if(input_string(1:7) == 'no_sort') input_no_sort = .true.
  endif


  IF(input_sort .or. input_no_sort) then
   if (n==0 .or. n==1) then
    if(hkl_statistics .and. write_details) call test_ratio_I_sig()
    !call statistics_on_E2()

    IF(keyword_CELL) then
     if(hkl_statistics .and. write_details) call stat_repartition_sintheta_lambda()    ! ok

     !call stat_repartition_d_hkl(F2, sig_F2, d_hkl, n_ref)                       ! nok
     IF(keyword_WAVE) then
      if(hkl_statistics .and. write_details) call stat_repartition_theta()                    ! ok
     endif

     IF(HKL_file%plot) then
      !call create_PGF_file(sintheta_lambda,F2, h,k,l, n_ref, "sinTheta")
      Pgf_file%name = 'F2_versus_sintheta_lambda.pgf'
      call create_PGF_file(TRIM(Pgf_file%name), sintheta_lambda,F2, HKL_string, n_ref, "sinTheta")
      !call system('winplotr F2_versus_sintheta_lambda.pgf', .true. )
      IF(LEN_TRIM(winplotr_exe) == 0) then
       call write_info('')
       call write_info(' Warning: WinPLOTR is not installed ')
       call write_info('')
      else
       !call system('winplotr '//trim(pgf_file%name), .true. )    ! lf95
       call system('winplotr '//trim(pgf_file%name))              ! g95
      endif
     END if

    endif
   endif
  endif


   IF(ALLOCATED(ordered_array)) DEALLOCATE(ordered_array)
   ALLOCATE(ordered_array(n_ref), stat=ier)
    if(ier/=0) call write_alloc_error("ordered_arrays")


   select case (input_string)
     case ('sort_stl')
       call write_info('')
       call write_info('  ... Sorting hkl in increasing SinTheta...')
       call write_info('')

       ! tri par ordre croissant des X
       call sort(sintheta_lambda, n_ref, ordered_array )     ! math_gen module
       call sort_arrays(input_string, '+', ordered_array, sintheta_lambda, 'out')

     case ('sort_d')
       call write_info('')
       call write_info(' ... Sorting hkl in increasing d_hkl...')
       call write_info('')

       call sort(d_hkl, n_ref, ordered_array)
       call sort_arrays(input_string,'+', ordered_array, d_hkl, 'out')
       call stat_repartition_d_hkl(F2, sig_F2, d_hkl, n_ref)

     case ('sort_theta')
       call write_info('')
       call write_info(' ... Sorting hkl in increasing theta...')
       call sort(theta_hkl, n_ref, ordered_array)
       call sort_arrays(input_string,'+', ordered_array, theta_hkl, 'out')


     case ('sort_int')
       call write_info('')
       call write_info('  ... Sorting hkl in decreasing intensities ...')
       call write_info('')

      call sort(I_, n_ref, ordered_array)
      call sort_arrays(input_string,'-', ordered_array, I_, 'out')

      F2(1:n_ref) = I_(1:n_ref)

     CASE ('sort_isig')
       call write_info('')
       call write_info('  ... Sorting hkl in decreasing intensities/sigma ...')
       call write_info('')

      call sort(I_sigma, n_ref, ordered_array)
      call sort_arrays(input_string,'-', ordered_array, I_sigma, 'out')

      F2(1:n_ref) = I_(1:n_ref)

     case default
   end select



   IF(ALLOCATED(ordered_array)) DEALLOCATE(ordered_array)

   IF(keyword_CELL) then
    min_stl = minval(sintheta_lambda(1:n_ref))   ! min_res not correct !!!!????
    max_stl = maxval(sintheta_lambda(1:n_ref))

    min_d_hkl = MINVAL(d_hkl(1:n_ref))
    max_d_hkl = MAXVAL(d_hkl(1:n_ref))

    min_theta = MINVAL(theta_hkl(1:n_ref))
    max_theta = MAXVAL(theta_hkl(1:n_ref))
   endif

   min_cut_off = minval(I_sigma(1:n_ref))
   max_cut_off = maxval(I_sigma(1:n_ref))

   intensity_min = minval(F2(1:n_ref))
   intensity_max = maxval(F2(1:n_ref))

   ! tres utile si on ne fait pas de cut_off
   if (keyword_CELL) then
    expected_min_stl = min_stl
    expected_max_stl = max_stl

    expected_min_d_hkl = min_d_hkl
    expected_max_d_hkl = max_d_hkl

    expected_min_theta = min_theta
    expected_max_theta = max_theta
   endif

   expected_min_cut_off = min_cut_off
   expected_max_cut_off = max_cut_off

   expected_min_intensity = intensity_min
   expected_max_intensity = intensity_max

   sig_coef = 1.
   ! -----------------------------------------------

   sort_X      = .false.
   X_stl       = .false.
   X_d         = .false.
   X_theta     = .false.
   X_int       = .false.
   X_isig      = .false.
   shell_X     = .false.
   shell_stl   = .false.
   shell_theta = .false.
   shell_d     = .false.
   shell_int   = .false.
   shell_isig  = .false.

   IF(keyword_CELL) then
    !sort_X      = .false.
    !X_stl       = .false.
    !X_d         = .false.
    !X_theta     = .false.
    !X_int       = .false.
    !X_isig      = .false.
    !shell_X     = .false.
    !shell_stl   = .false.
    !shell_theta = .false.
    !shell_d     = .false.
    !shell_int   = .false.
    !shell_isig  = .false.

    if(long_input_string >= 5) then
     if(input_string(1:5) == 'sort_') sort_X = .true.
    endif
    if(long_input_string == 7) then
     if(input_string(1:7) == 'shell_d') shell_d = .true.
    endif
    if(long_input_string == 9) then
     if(input_string(1:9) == 'shell_stl') shell_stl = .true.
     if(input_string(1:9) == 'shell_int') shell_int = .true.
    endif
    if(long_input_string == 10) then
     if(input_string(1:10) == 'shell_isig') shell_isig = .true.
    endif
    if(long_input_string == 11) then
     if(input_string(1:11) == 'shell_theta') shell_theta = .true.
    endif
    if (input_string(long_input_string-2:long_input_string) == 'stl')   X_stl   = .true.
    if (input_string(long_input_string  :long_input_string) == 'd')     X_d     = .true.
    if (input_string(long_input_string-4:long_input_string) == 'theta') X_theta = .true.
    if (input_string(long_input_string-2:long_input_string) == 'int')   X_int   = .true.
    if (input_string(long_input_string-3:long_input_string) == 'isig')  X_isig  = .true.


    IF(sort_X .or. shell_stl) then
     call write_info('')
     write(message_text,'(a,2(2x,F8.5))') '  Experimental (sinTheta/lambda) shell: ',   min_stl, max_stl
     call write_info(TRIM(message_text))
     call write_info('')

     if (X_stl) then
      IF(.NOT. shell_arg_min_max) then
       call write_info(' > Enter expected (sinTheta/lambda) min. and max. values [0. 0.: exp. values]: ')
       call read_input_line(input_line)
       read(input_line,*, iostat=ier)  val_1, val_2
       if(ier ==0) then
        expected_min_stl = val_1
        expected_max_stl = val_2
       else
        return
       endif

       if (expected_min_stl < eps  .and. expected_max_stl < eps) then      ! 0. 0.
        expected_min_stl = min_stl
        expected_max_stl = max_stl
       end if
       if (expected_min_stl > expected_max_stl) then
        call write_info('')
        call write_info(' !!! Min_value > Max_value !!!')
        stop
       endif
      else
       expected_min_stl = shell_min(nb_shell)
       expected_max_stl = shell_max(nb_shell)
       WRITE(message_text,'(a,2F10.5)') '  >> Expected (sinTheta/lambda) min. and max. values: ', expected_min_stl, expected_max_stl
       call write_info(trim(message_text))
      endif
     endif
    endif

    IF(sort_X .or. shell_d) then
     write(message_text,'(a,2(2x,F8.5))') '  Experimental d_hkl(A) shell:          ',   min_d_hkl, max_d_hkl
     call write_info(TRIM(message_text))
     call write_info('')
     if (X_d) then
      IF(.NOT. shell_arg_min_max) then
       call write_info('  > Enter expected d_hkl min. and max. values [0. 0.: exp. values]: ')
       call read_input_line(input_line)
       read(input_line,*, iostat=ier)  val_1, val_2
       if(ier == 0) then
        expected_min_d_hkl = val_1
        expected_max_d_hkl = val_2
       else
        return
       end if

       if (expected_min_d_hkl < eps  .and. expected_max_d_hkl < eps) then      ! 0. 0.
        expected_min_d_hkl = min_d_hkl
        expected_max_d_hkl = max_d_hkl
       end if
       if (expected_min_d_hkl > expected_max_d_hkl) then
        call write_info('')
        call write_info(' !!! Min_value > Max_value !!!')
        stop
       endif
      else
       expected_min_d_hkl = shell_min(nb_shell)
       expected_max_d_hkl = shell_max(nb_shell)
       WRITE(message_text,'(a,2F10.5)') '  >> Expected d_hkl min. and max. values: ',  expected_min_d_hkl, expected_max_d_hkl
       call write_info(trim(message_text))
      endif
     endif
    endif

    IF(sort_X .or. shell_theta) then
     write(message_text,'(a,2(2x,F8.3))') '  Experimental theta(deg) shell:        ',   min_theta, max_theta
     call write_info(TRIM(message_text))
     call write_info('')
     if (X_theta) then
      IF(.NOT. shell_arg_min_max) then
       call write_info(' > Enter expected theta min. and max. values [0. 0.: exp. values]: ')
       call read_input_line(input_line)
       read(input_line,*, iostat=ier)  val_1, val_2
       if(ier==0) then
        expected_min_theta = val_1
        expected_max_theta = val_2
       else
        return
       end if

       if (expected_min_theta < eps  .and. expected_max_theta < eps) then      ! 0. 0.
        expected_min_theta = min_theta
        expected_max_theta = max_theta
       end if
       if (expected_min_theta > expected_max_theta) then
        call write_info('')
        call write_info(' !!! Min_value > Max_value !!!')
        stop
       endif
      else
       expected_min_theta = shell_min(nb_shell)
       expected_max_theta = shell_max(nb_shell)
       WRITE(message_text, '(a, 2F10.5)') '  >>  Expected theta min. and max. values: ', expected_min_theta, expected_max_theta
       call write_info(trim(message_text))
      endif
     endif
    endif
   ENDIF ! fin de la condition if (keyword_CELL


   IF(sort_X .or. shell_int) then
    IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
     write(message_text,'(a,2F8.2)')      '  Experimental Intensities shell:        ', intensity_min, intensity_max
    else
     write(message_text,'(a,2F9.1)')      '  Experimental Intensities shell:        ', intensity_min, intensity_max
    endif
    call write_info(TRIM(message_text))
    call write_info('')

    if (X_int) then
     IF(.NOT. shell_arg_min_max) then
      call write_info(' > Enter expected Intensity cut-off min. and max. values [0. 0.: exp. values]: ')
      call read_input_line(input_line)
      read(input_line,*, iostat=ier)  val_1, val_2
      if(ier == 0) then
       expected_min_intensity = val_1
       expected_max_intensity = val_2
      else
       return
      end if
      if (expected_min_intensity < eps .and. expected_max_intensity < eps) then
       expected_min_intensity = intensity_min
       expected_max_intensity = intensity_max
      endif
      if (expected_min_intensity > expected_max_intensity) then
       call write_info('')
       call write_info(' !!! Min_value > Max_value !!!')
      endif
     else
      expected_min_intensity = shell_min(nb_shell)
      expected_max_intensity = shell_max(nb_shell)
      WRITE(message_text, '(a,2F10.5)') '  >>  Expected Intensity cut-off min. and max. values: ', &
                                         expected_min_intensity, expected_max_intensity
      call write_info(trim(message_text))
     endif
    endif
   endif

   IF(sort_X .or. shell_isig) then
    call write_info('')
    IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
     write(message_text,'(a,2F8.2)')       '  Experimental (I/sig) shell:           ', min_cut_off, max_cut_off
    else
     write(message_text,'(a,2F9.1)')       '  Experimental (I/sig) shell:           ', min_cut_off, max_cut_off
    endif
    call write_info(TRIM(message_text))
    call write_info('')

    IF(X_isig) then
     IF(.NOT. shell_arg_min_max) then
      call write_info(' > Enter expected I/sig cut-off min. and max. values [0. 0.: exp. values]: ')
      call read_input_line(input_line)
      read(input_line,*, iostat=ier)  val_1, val_2
      if(ier == 0) then
       expected_min_cut_off = val_1
       expected_max_cut_off = val_2
      else
       return
      end if
      if (expected_min_cut_off < eps .and. expected_max_cut_off < eps) then
       expected_min_cut_off = min_cut_off
       expected_max_cut_off = max_cut_off
      endif
      if (expected_min_cut_off > expected_max_cut_off) then
       call write_info('')
       call write_info('!!! Min_value > Max_value !!!')
      endif
     else
      expected_min_cut_off = shell_min(nb_shell)
      expected_max_cut_off = shell_max(nb_shell)
      WRITE(message_text, '(a,2F10.5)') '  >> Expected I/sig cut-off min. and max. values: ', &
                                        expected_min_cut_off, expected_max_cut_off
      call write_info(trim(message_text))
     endif
    endif
   endif

   call write_info('')
   IF(cut_off) then
    call write_info(' > Enter divider sigma coefficient: ')
    call read_input_line(input_line)
    read(input_line,*) sig_coef
    if (sig_coef < eps)   sig_coef = 1.

    write(message_text,*, iostat=ier) ' > Divider sigma coefficient: ' , val_1
    if(ier == 0) then
     sig_coef = val_1
    else
     return
    end if
    call write_info(TRIM(message_text))
   else
    sig_coef = 1.
   endif


   ! creation fichier .HKL / .CIF
   IF(M95_comment_line_nb/=0) then
    WRITE(12, '(a)') 'import_report_begin'
    do i=1,8
     WRITE(12, '(a)') TRIM(M95_comment_line(i))
    end do
   endif

   call Allocate_PGF_data_arrays(n_ref)

    !call write_info('')
    select case (input_string)
     case ('sort_d', 'shell_d')
      n_ref_eff      = 0
      n_ref_excluded = 0
      do i=1, n_ref
       if(d_hkl(i) >= expected_min_d_hkl .and. d_hkl(i) <= expected_max_d_hkl)  then
        call write_sorted_file(i, ok)
        if(.not. ok) return
        n_ref_eff = n_ref_eff + 1
        pgf_data%X(n_ref_eff) = d_hkl(i)
        pgf_data%Y(n_ref_eff) = F2(i)
        pgf_data%h(n_ref_eff) = h(i)
        pgf_data%k(n_ref_eff) = k(i)
        pgf_data%l(n_ref_eff) = l(i)
        WRITE(pgf_data%string(n_ref_eff), '(a,3I4,a)') '(', h(i), k(i), l(i), ')'

        IF(sort_out(n) .AND. n_ref_eff <=sort_out_n(n)) then
         if(i==1) then
          call write_info('    h   k   l      F2     sig          d (A)')
          call write_info('')
         endif
         IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
          WRITE(message_text, '(1x,3I4,2F8.2,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), d_hkl(i)
         else
          WRITE(message_text, '(1x,3I4,2F9.1,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), d_hkl(i)
         endif
         call write_info(TRIM(message_text))
        endif
       else
        n_ref_excluded = n_ref_excluded + 1
        if(shell_out) then
         IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
          WRITE(message_text, '(1x,3I4,2F8.2,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), d_hkl(i)
         else
          WRITE(message_text, '(1x,3I4,2F9.1,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), d_hkl(i)
         endif
         if(n_ref_excluded == 1) then
          call write_info("  List of excluded reflections:")
          call write_info("    h   k   l      F2  sig_F2          d_hkl")
          call write_info("")
         end if
         call write_info(TRIM(message_text))
        end if
       endif
      end do

      IF(sort_out(n) .or. shell_out) call write_info(' ')
      call write_info('  >>> '//TRIM(HKL_file%d)// ' file has been created.')
      WRITE(message_text, '(a,I8)')          '     . number of reflections          : ', n_ref_eff
      call write_info(TRIM(message_text))
      write(message_text, '(a,I8,a,F7.3,a)') '     . number of excluded reflections : ', n_ref - n_ref_eff, &
                                             ' (', 100.*(n_ref - n_ref_eff)/n_ref, '%)'
      call write_info(trim(message_text))
      call write_info(' ')
      IF(sort_plot(n) .or. shell_plot(n)) then
       !call create_PGF_file(pgf_data%X, pgf_data%Y, pgf_data%h, pgf_data%k, pgf_data%l, n_ref_eff, "d_hkl")
       pgf_file%name = 'F2_versus_d_hkl.pgf'
       call create_PGF_file(TRIM(Pgf_file%name), pgf_data%X, pgf_data%Y, pgf_data%string, n_ref_eff, "d_hkl")
       IF(LEN_TRIM(winplotr_exe) == 0) then
        call write_info('')
        call write_info(' Warning: WinPLOTR is not installed ')
        call write_info('')
       else
        !call system('winplotr '//trim(pgf_file%name), .true. )  ! lf95
        call system('winplotr '//trim(pgf_file%name))            ! g95
       endif
      END if


     case ('sort_stl', 'shell_stl')
      n_ref_eff      = 0
      n_ref_excluded = 0
      do i=1, n_ref
       if(sinTheta_lambda(i) >= expected_min_stl .and. sinTheta_lambda(i) <= expected_max_stl) then
        call write_sorted_file(i, ok)
        if(.not. ok) return
        n_ref_eff = n_ref_eff + 1
        pgf_data%X(n_ref_eff) = sinTheta_lambda(i)
        pgf_data%Y(n_ref_eff) = F2(i)
        pgf_data%h(n_ref_eff) = h(i)
        pgf_data%k(n_ref_eff) = k(i)
        pgf_data%l(n_ref_eff) = l(i)
        WRITE(pgf_data%string(n_ref_eff), '(a,3I4,a)') '(', h(i), k(i), l(i), ')'
        IF(sort_out(n) .AND. n_ref_eff<=sort_out_n(n)) then
         if(i==1) then
          call write_info('    h   k   l      F2     sig      stl (A-1)')
          call write_info('')
         endif
         IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
          WRITE(message_text, '(1x,3I4,2F8.2,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), sinTheta_lambda(i)
         else
          WRITE(message_text, '(1x,3I4,2F9.1,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), sinTheta_lambda(i)
         endif
         call write_info(TRIM(message_text))
        endif
       else
        n_ref_excluded = n_ref_excluded + 1
        if(shell_out) then
         IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
          WRITE(message_text, '(1x,3I4,2F8.2,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), sinTheta_lambda(i)
         else
          WRITE(message_text, '(1x,3I4,2F9.1,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), sinTheta_lambda(i)
         endif
         if(n_ref_excluded == 1) then
          call write_info("  List of excluded reflections:")
          call write_info("    h   k   l      F2  sig_F2        stl_hkl")
          call write_info("")
         end if
         call write_info(TRIM(message_text))
        end if
       endif
      end do

      IF(sort_out(n) .or. shell_out) call write_info(' ')
      call write_info('  >>> '//TRIM(HKL_file%stl)// ' file has been created.')
      WRITE(message_text, '(a,I8)')          '     . number of reflections          : ', n_ref_eff
      call write_info(TRIM(message_text))
      write(message_text, '(a,I8,a,F7.3,a)') '     . number of excluded reflections : ', n_ref - n_ref_eff, &
                                             ' (', 100.*(n_ref - n_ref_eff)/n_ref, '%)'
      call write_info(trim(message_text))

      call write_info(' ')
      IF(sort_plot(n) .or. shell_plot(n)) then
       !call create_PGF_file(pgf_data%X, pgf_data%Y, pgf_data%h, pgf_data%k, pgf_data%l, n_ref_eff, "sinTheta")
       pgf_file%name = 'F2_versus_sintheta_lambda.pgf'
       call create_PGF_file(TRIM(Pgf_file%name), pgf_data%X, pgf_data%Y, pgf_data%string, n_ref_eff, "sinTheta")
       IF(LEN_TRIM(winplotr_exe) == 0) then
        call write_info('')
        call write_info(' Warning: WinPLOTR is not installed ')
        call write_info('')
       else
        !call system('winplotr '//trim(pgf_file%name), .true. )   ! lf95
        call system('winplotr '//trim(pgf_file%name))             ! g95
       endif
      END if


     case ('sort_theta', 'shell_theta')
      n_ref_eff      = 0
      n_ref_excluded = 0
      do i=1, n_ref
       if(Theta_hkl(i) >= expected_min_theta .and. theta_hkl(i) <= expected_max_theta) then
        call write_sorted_file(i, ok)
        if(.not. ok) return
        n_ref_eff = n_ref_eff + 1
        pgf_data%X(n_ref_eff)= Theta_hkl(i)
        pgf_data%Y(n_ref_eff) = F2(i)
        pgf_data%h(n_ref_eff) = h(i)
        pgf_data%k(n_ref_eff) = k(i)
        pgf_data%l(n_ref_eff) = l(i)
        WRITE(pgf_data%string(n_ref_eff), '(a,3I4,a)') '(', h(i), k(i), l(i), ')'
        IF(sort_out(n) .AND. n_ref_eff <=sort_out_n(n)) then
         if(i==1) then
          call write_info('    h   k   l      F2     sig    theta (deg)')
          call write_info('')
         endif
         IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
          WRITE(message_text, '(1x,3I4,2F8.2,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), theta_hkl(i)
         else
          WRITE(message_text, '(1x,3I4,2F9.1,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), theta_hkl(i)
         endif
         call write_info(TRIM(message_text))
        endif
       else
        n_ref_excluded = n_ref_excluded + 1
        if(shell_out) then
          IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
          WRITE(message_text, '(1x,3I4,2F8.2,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), theta_hkl(i)
         else
          WRITE(message_text, '(1x,3I4,2F9.1,F15.5)')  h(i), k(i), l(i), F2(i), sig_F2(i), theta_hkl(i)
         endif
         if(n_ref_excluded == 1) then
          call write_info("  List of excluded reflections:")
          call write_info("    h   k   l      F2  sig_F2      theta_hkl")
          call write_info("")
         end if
         call write_info(TRIM(message_text))
        end if
       END if
      end do
      IF(sort_out(n) .or. shell_out) call write_info(' ')
      call write_info('  >>> '//TRIM(HKL_file%theta)// ' file has been created.')
      WRITE(message_text, '(a,I8)')          '     . number of reflections          : ', n_ref_eff
      call write_info(TRIM(message_text))
      write(message_text, '(a,I8,a,F7.3,a)') '     . number of excluded reflections : ', n_ref - n_ref_eff, &
                                             ' (', 100.*(n_ref - n_ref_eff)/n_ref, '%)'
      call write_info(trim(message_text))
      call write_info(' ')
      IF(sort_plot(n) .or. shell_plot(n)) then
      ! call create_PGF_file(pgf_data%X, pgf_data%Y, pgf_data%h, pgf_data%k, pgf_data%l, n_ref_eff, "sinTheta")
       pgf_file%name = 'F2_versus_Theta.pgf'
       call create_PGF_file(TRIM(Pgf_file%name), pgf_data%X, pgf_data%Y, pgf_data%string, n_ref_eff, "theta")
       IF(LEN_TRIM(winplotr_exe) == 0) then
        call write_info('')
        call write_info(' Warning: WinPLOTR is not installed ')
        call write_info('')
       else
        !call system('winplotr '//trim(pgf_file%name),  .true. )   ! lf95
        call system('winplotr '//trim(pgf_file%name))              ! g95
       endif
      END if

     case ('sort_int', 'shell_int')
      if(debug_proc%level_2)  call write_debug_proc_level(2, "WRITE_SORTED_FILE")
      n_ref_eff      = 0
      n_ref_excluded = 0
      do i=1, n_ref
       if(F2(i) >= expected_min_intensity .and. F2(i)  <= expected_max_intensity) then
        call write_sorted_file(i, ok)
        if(.not. ok) return
        n_ref_eff = n_ref_eff + 1
       else
        n_ref_excluded = n_ref_excluded + 1
        if(shell_out) then
          IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
          WRITE(message_text, '(1x,3I4,2F8.2)')  h(i), k(i), l(i), F2(i), sig_F2(i)
         else
          WRITE(message_text, '(1x,3I4,2F9.1)')  h(i), k(i), l(i), F2(i), sig_F2(i)
         endif
         if(n_ref_excluded == 1) then
          call write_info("  List of excluded reflections:")
          call write_info("    h   k   l      F2  sig_F2")
          call write_info("")
         end if
         call write_info(TRIM(message_text))
        end if
       end if
      end do
      IF(sort_out(n) .or. shell_out) call write_info(' ')
      call write_info('  >>> '//TRIM(HKL_file%I)// ' file has been created.')
      WRITE(message_text, '(a,I8)')          '     . number of reflections          : ', n_ref_eff
      call write_info(TRIM(message_text))
      write(message_text, '(a,I8,a,F7.3,a)') '     . number of excluded reflections : ', n_ref - n_ref_eff, &
                                             ' (', 100.*(n_ref - n_ref_eff)/n_ref, '%)'
      call write_info(trim(message_text))
      call write_info(' ')


     case ('sort_isig', 'shell_isig')
      n_ref_eff      = 0
      n_ref_excluded = 0
      do i=1, n_ref
       if(F2(i)/sig_F2(i) >= expected_min_cut_off .and. F2(i)/sig_F2(i)  <= expected_max_cut_off) then
        call write_sorted_file(i, ok)
        if(.not. ok) return
        n_ref_eff = n_ref_eff + 1
       else
        n_ref_excluded = n_ref_excluded + 1
        if(shell_out) then
          IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
           WRITE(message_text, '(1x,3I4,3F8.2)')  h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i)
          else
          WRITE(message_text, '(1x,3I4,3F9.1)')   h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i)
          endif
          if(n_ref_excluded == 1) then
           call write_info("  List of excluded reflections:")
           call write_info("    h   k   l      F2  sig_F2  F2/sig")
           call write_info("")
          end if
         call write_info(TRIM(message_text))
        end if
       end if
      end do
      IF(sort_out(n) .or. shell_out) call write_info(' ')
      call write_info('  >>> '//TRIM(HKL_file%Isig)// ' file has been created.')
      WRITE(message_text, '(a,I8)')          '     . number of reflections          : ', n_ref_eff
      call write_info(TRIM(message_text))
      write(message_text, '(a,I8,a,F7.3,a)') '     . number of excluded reflections : ', n_ref - n_ref_eff, &
                                             ' (', 100.*(n_ref - n_ref_eff)/n_ref, '%)'
      call write_info(trim(message_text))
      call write_info(' ')

    end select



 RETURN
end subroutine classement


!---------------------------------------------------------------------------------------------

subroutine calcul_Q(h,k,l,cell,Q)
 use cryscalc_module, only : debug_proc
 implicit none
  real,               intent(in)       :: h,k,l
  real, dimension(6), intent(in)       :: cell
  real,               intent(out)      :: Q
  real                                 :: cos_alfa, cos_beta, cos_gama
  real                                 :: cos_a2, cos_b2, cos_g2
  real                                 :: volume
  REAL                                 :: pi
  real, dimension(3)                   :: cell_star      ! parametre de la maille reciproque
  real, dimension(3)                   :: cell_rad, cos_cell_star
  real                                 :: Q1, Q2, Q3, Q4, Q5, Q6

  !if(debug_proc%level_3)  call write_debug_proc_level(3, "CALCUL_Q")


  pi = 4 * ATAN(1.)

 ! Qhkl = (h.as + k.bs + l.cs)**2
 !      = h**2.as**2 + k**2.bs**2 + l**2.cs**2
 !        + 2hk.as.bs.cosgs
 !        + 2kl.bs.cs.cosas
 !        + 2lh.cs.as.cosbs

 ! calcul du volume de la maille
  Cos_Alfa = COS(cell(4) * pi / 180.)
  Cos_beta = COS(cell(5) * pi / 180.)
  Cos_gama = COS(cell(6) * pi / 180.)

  cos_a2 = Cos_Alfa ** 2
  cos_b2 = cos_Beta ** 2
  cos_g2 = cos_Gama ** 2

  volume = SQRT(1 + 2 * Cos_Alfa * cos_beta * cos_gama - cos_a2 - cos_b2 - cos_g2)
  volume = volume * cell(1) * cell(2) * cell(3)

 ! calcul des parametres de la maille reciproque



  cell_rad(1)  = cell(4) * pi / 180.
  cell_rad(2)  = cell(5) * pi / 180.
  cell_rad(3)  = cell(6) * pi / 180.

  cell_star(1) = cell(2) * cell(3) * SIN(cell_rad(1)) / Volume
  cell_star(2) = cell(3) * cell(1) * SIN(cell_rad(2)) / Volume
  cell_star(3) = cell(1) * cell(2) * SIN(cell_rad(3)) / Volume

  cos_cell_star(1) = (cos(cell_rad(2)) * cos(cell_rad(3)) - cos(cell_rad(1))) / (sin(cell_rad(2)) * sin(cell_rad(3)))
  cos_cell_star(2) = (cos(cell_rad(3)) * cos(cell_rad(1)) - cos(cell_rad(2))) / (sin(cell_rad(3)) * sin(cell_rad(1)))
  cos_cell_star(3) = (cos(cell_rad(1)) * cos(cell_rad(2)) - cos(cell_rad(3))) / (sin(cell_rad(1)) * sin(cell_rad(2)))

  Q1 = (h**2)*(cell_star(1)**2)
  Q2 = (k**2)*(cell_star(2)**2)
  Q3 = (l**2)*(cell_star(3)**2)
  Q4 = 2*h*k*cell_star(1)*cell_star(2)*cos_cell_star(3)
  Q5 = 2*k*l*cell_star(2)*cell_star(3)*cos_cell_star(1)
  Q6 = 2*l*h*cell_star(3)*cell_star(1)*cos_cell_star(2)
  Q  = Q1 + Q2 + Q3 + Q4 + Q5 + Q6


 return
end subroutine calcul_Q

!-----------------------------------------------------------------------------------------------------------------------

 subroutine classement_par_ordre_croissant(sort_type, N, Z, indice_h, indice_k, indice_l, col4, col5, i_code,  &
                                           cos1, cos2, cos3, cos4, cos5, cos6)


 ! classement d'un tableau X,Y... par ordre croissant des X
  USE cryscalc_module, ONLY : message_text, debug_proc
  USE macros_module,   ONLY : multiple
  USE IO_module,       ONLY : write_info

  implicit none
   character  (len=*),              intent(in)              :: sort_type
   integer,                         INTENT(IN)              :: N
   REAL,              DIMENSION(N), INTENT(INOUT)           :: Z, col4, col5
   integer,           dimension(N), intent(inout)           :: indice_h, indice_k, indice_l, i_code
   REAL,              DIMENSION(N), INTENT(INOUT), OPTIONAL :: cos1, cos2, cos3, cos4, cos5, cos6

   !local variables
   character (len=32)             :: comment
   integer                        :: i, j
   real                           :: real_temp
   integer                        :: i_temp

   if(debug_proc%level_2)  call write_debug_proc_level(2, "CLASSEMENT_PAR_ORDRE_CROISSANT")

   if (sort_type(1:3) == 'stl') then
    comment = 'Sintheta/lambda'
   ELSEIF(sort_type(1:1)== 'd') then
    comment = 'd_hkl'
   ELSEIF(sort_type(1:5) == 'theta') then
    comment = 'theta'
   ELSEIF(sort_type(1:1) == 'i') then
    comment = 'intensity'
   ELSEIF(sort_type(1:4) == 'isig') then
    comment = 'I/sig'
   endif


   do i=1 , N -1
    do j=i+1 , N
     if (Z(i) > Z(j)) then
        real_temp = Z(i)
        Z(i) = Z(j)
        Z(j) = real_temp

        i_temp = indice_h(i)
        indice_h(i) = indice_h(j)
        indice_h(j) = i_temp

        i_temp = indice_k(i)
        indice_k(i) = indice_k(j)
        indice_k(j) = i_temp

        i_temp = indice_l(i)
        indice_l(i) = indice_l(j)
        indice_l(j) = i_temp

        real_temp = col4(i)
        col4(i) = col4(j)
        col4(j) = real_temp

        real_temp = col5(i)
        col5(i) = col5(j)
        col5(j) = real_temp

        i_temp = i_code(i)
        i_code(i) = i_code(j)
        i_code(j) = i_temp

        if (PRESENT(cos1)) then
         real_temp = cos1(i)
         cos1(i) = cos1(j)
         cos1(j) = real_temp
        endif

        if (PRESENT(cos2)) then
         real_temp = cos2(i)
         cos2(i) = cos2(j)
         cos2(j) = real_temp
        endif

        if (PRESENT(cos3)) then
         real_temp = cos3(i)
         cos3(i) = cos3(j)
         cos3(j) = real_temp
        endif

        if (PRESENT(cos4)) then
         real_temp = cos4(i)
         cos4(i) = cos4(j)
         cos4(j) = real_temp
        endif

        if (PRESENT(cos5)) then
         real_temp = cos5(i)
         cos5(i) = cos5(j)
         cos5(j) = real_temp
        endif

        if (PRESENT(cos6)) then
         real_temp = cos6(i)
         cos6(i) = cos6(j)
         cos6(j) = real_temp
        endif

     END if
    END do

    !if (1000*int(i/1000) == i ) then
    if(multiple(i, 1000)) then
     write(message_text,'(a,i6,a,a,a,F10.5,a,3I4,a)') '   ...',i, ' ...   ',trim(comment), ' = ', Z(i), &
                                                     '  (',indice_h(i), indice_k(i),indice_l(i),')'
     call write_info(TRIM(message_text))
    endif
   END do

   write(message_text,'(a,i6,a,a,a,F10.5,a,3I4,a)')   '   ...',n, ' ...   ',trim(comment), ' = ', Z(n), &
                                                     '  (',indice_h(n), indice_k(n),indice_l(n),')'
   call write_info(TRIM(message_text))

   return
 end subroutine classement_par_ordre_croissant


!-------------------------------------------------------------------

 subroutine classement_par_ordre_decroissant(sort_type, N, Z, indice_h, indice_k, indice_l, col4, col5, &
                                             i_code,cos1, cos2, cos3, cos4, cos5, cos6)

   USE cryscalc_module, ONLY : message_text, debug_proc
   USE macros_module,   ONLY : multiple
   USE IO_module,       ONLY : write_info

 ! classement d'un tableau X,Y... par ordre decroissant des X

  implicit none
   character  (len=*),              intent(in)              :: sort_type
   integer,                         INTENT(IN)              :: N
   REAL,              DIMENSION(N), INTENT(INOUT)           :: Z, col4, col5
   integer,           dimension(N), intent(inout)           :: indice_h, indice_k, indice_l, i_code
   REAL,              DIMENSION(N), INTENT(INOUT), OPTIONAL :: cos1, cos2, cos3, cos4, cos5, cos6

   ! local variables
   character (len=32)             :: comment
   integer                        :: i, j
   real                           :: real_temp
   integer                        :: i_temp

   if(debug_proc%level_2)  call write_debug_proc_level(2, "CLASSEMENT_PAR_ORDRE_DECROISSANT")

   if (sort_type(1:3) == 'int' ) then
    comment = 'Intensity'
   elseif(sort_type(1:5) == 'isig') then
    comment = 'I/sig'
   end if

   do i=1 , N -1

    do j=i+1 , N
     if (Z(i) < Z(j)) then
        real_temp = Z(i)
        Z(i) = Z(j)
        Z(j) = real_temp

        i_temp = indice_h(i)
        indice_h(i) = indice_h(j)
        indice_h(j) = i_temp

        i_temp = indice_k(i)
        indice_k(i) = indice_k(j)
        indice_k(j) = i_temp

        i_temp = indice_l(i)
        indice_l(i) = indice_l(j)
        indice_l(j) = i_temp

        real_temp = col4(i)
        col4(i) = col4(j)
        col4(j) = real_temp

        real_temp = col5(i)
        col5(i) = col5(j)
        col5(j) = real_temp

        i_temp = i_code(i)
        i_code(i) = i_code(j)
        i_code(j) = i_temp

        if (PRESENT(cos1)) then
         real_temp = cos1(i)
         cos1(i) = cos1(j)
         cos1(j) = real_temp
        endif

        if (PRESENT(cos2)) then
         real_temp = cos2(i)
         cos2(i) = cos2(j)
         cos2(j) = real_temp
        endif

        if (PRESENT(cos3)) then
         real_temp = cos3(i)
         cos3(i) = cos3(j)
         cos3(j) = real_temp
        endif

        if (PRESENT(cos4)) then
         real_temp = cos4(i)
         cos4(i) = cos4(j)
         cos4(j) = real_temp
        endif

        if (PRESENT(cos5)) then
         real_temp = cos5(i)
         cos5(i) = cos5(j)
         cos5(j) = real_temp
        endif

        if (PRESENT(cos6)) then
         real_temp = cos6(i)
         cos6(i) = cos6(j)
         cos6(j) = real_temp
        endif

     END if
    END do

    !if (1000*int(i/1000) == i) then
    if(multiple(i, 1000)) then
     write(message_text,'(a,i6,a,a,a,F10.5)') '   ...',i, ' ...   ',trim(comment), ' = ', Z(i)
     call write_info(TRIM(message_text))
    endif

   END do

   write(message_text,'(a,i6,a,a,a,F10.5)') '   ...',n, ' ...   ',trim(comment), ' = ', Z(n)
   call write_info(TRIM(message_text))

   return
 end subroutine classement_par_ordre_decroissant

!--------------------------------------------------------------------------
subroutine sort_arrays(sort_type, sort_sign, ordered, Z, input_string)

 ! classement d'un tableau X,Y... par ordre croissant des X
   USE cryscalc_module, ONLY : message_text, max_ref, WRITE_details, debug_proc
   USE IO_module,       ONLY : write_info
   USE macros_module,   ONLY : multiple
   USE HKL_module,      ONLY : n_ref, h, k,l, F2, sig_F2, code, cos_dir, HKL_flag

  implicit none
   character  (len=*),                intent(in)              :: sort_type
   CHARACTER  (LEN=*),                INTENT(IN)              :: sort_sign
   INTEGER,       DIMENSION(N_ref),   INTENT(IN)              :: ordered
   REAL,          DIMENSION(N_ref),   INTENT(INOUT)           :: Z
   CHARACTER(len=*),                  intent(in)              :: input_string

   !local variables
   character (len=32)             :: comment
   integer                        :: i, rang
   REAL,    allocatable, DIMENSION(:)      :: temp_Z, temp_F2, temp_sig_F2
   REAL,    allocatable, DIMENSION(:,:)    :: temp_cos
   integer, allocatable, DIMENSION(:,:)    :: temp_H
   INTEGER, allocatable, DIMENSION(:)      :: temp_code
   integer                                 :: ier

   if(debug_proc%level_2)  call write_debug_proc_level(2, "SORT_ARRAYS ("//trim(input_string)//")")

   if (ALLOCATED(temp_Z))      DEALLOCATE(temp_Z)
   if (ALLOCATED(temp_F2))     DEALLOCATE(temp_F2)
   if (ALLOCATED(temp_sig_F2)) DEALLOCATE(temp_sig_F2)
   if (ALLOCATED(temp_cos))    DEALLOCATE(temp_cos)
   if (ALLOCATED(temp_H))      DEALLOCATE(temp_H)
   if (ALLOCATED(temp_code))   DEALLOCATE(temp_code)

   ALLOCATE(temp_Z(Max_ref), stat=ier)
    if(ier/=0) call write_alloc_error("temp_Z")
   ALLOCATE(temp_F2(Max_ref), stat=ier)
    if(ier/=0) call write_alloc_error("temp_F2")
   ALLOCATE(temp_sig_F2(Max_ref), stat=ier)
    if(ier/=0) call write_alloc_error("temp_CIF")
   ALLOCATE(temp_cos(Max_ref,6), stat=ier)
    if(ier/=0) call write_alloc_error("temp_COS")
   ALLOCATE(temp_code(Max_ref), stat=ier)
    if(ier/=0) call write_alloc_error("temp_code")
   ALLOCATE(temp_H(Max_ref,3), stat=ier)
    if(ier/=0) call write_alloc_error("temp_H")


   select case (sort_type)
     case ('sort_stl')
      comment = 'Sintheta/lambda'

     case ('sort_d')
      comment = 'd_hkl'

     case ('sort_theta')
      comment = 'theta(deg)'

     case ('sort_int')
      comment = 'Intensity'

     case ('sort_isig')
      comment = 'I/sig'

     case default
      comment = ''
   end select

   select case (sort_sign)
       case ("+")
        do i=1, n_ref
         rang = ordered(i)
         temp_Z(i)        = Z(rang)
         temp_H(i,1)      = h(rang)
         temp_H(i,2)      = k(rang)
         temp_H(i,3)      = l(rang)
         temp_F2(i)       = F2(rang)
         temp_sig_F2(i)   = sig_F2(rang)
         temp_code(i)     = code(rang)
         temp_cos(i,1:6)  = cos_dir(rang,1:6)

         !if (input_string(1:3) =='out' .and. 1000*int(i/1000) == i ) then
         if(input_string(1:3) == 'out' .and. multiple(i, 1000)) then
          write(message_text,'(a,i6,a,a,a,F10.5,a,3I4,a)')   '   ...',i, ' ...   ',trim(comment), ' = ', temp_Z(i), &
                                                             '  (',temp_H(i,1), temp_H(i,2),temp_H(i,3),')'
          call write_info(TRIM(message_text))
         endif
        end do
        !Z       = temp_Z
        !h       = temp_H(:,1)
        !k       = temp_H(:,2)
        !l       = temp_H(:,3)
        !F2      = temp_F2
        !sig_F2  = temp_sig_F2
        !code(1:n_ref)    = temp_code(1:n_ref)
        !cos_dir = temp_cos


       case ("-")
        do i=1, n_ref
         rang = ordered(n_ref-i+1)
         temp_Z(i)        = Z(rang)
         temp_H(i,1)      = h(rang)
         temp_H(i,2)      = k(rang)
         temp_H(i,3)      = l(rang)
         temp_F2(i)       = F2(rang)
         temp_sig_F2(i)   = sig_F2(rang)
         temp_code(i)     = code(rang)
         temp_cos(i,1:6)  = cos_dir(rang,1:6)
         !if (input_string(1:3) =='out' .and. 1000*int(i/1000) == i ) then
         if(input_string(1:3) == 'out' .and. multiple(i, 1000)) then
          write(message_text,'(a,i6,a,a,a,F10.5,a,3I4,a)')   '   ...',i, ' ...   ',trim(comment), ' = ', temp_Z(i), &
                                                             '  (',temp_H(i,1), temp_H(i,2),temp_H(i,3),')'
          call write_info(TRIM(message_text))
         endif
        end do
        !Z       = temp_Z
        !h       = temp_h(:,1)
        !k       = temp_h(:,2)
        !l       = temp_h(:,3)
        !F2      = temp_F2
        !sig_F2  = temp_sig_F2
        !code    = temp_code
        !cos_dir = temp_cos

   end select

   Z       = temp_Z
   h       = temp_H(:,1)
   k       = temp_H(:,2)
   l       = temp_H(:,3)
   F2      = temp_F2
   sig_F2  = temp_sig_F2
   code    = temp_code
   cos_dir = temp_cos
   !HKL_flag(1:n_ref, 1) = ordered(1:n_ref)


   if(input_string(1:3) =='out') then
   write(message_text,'(a,i6,a,a,a,F10.5,a,3I4,a)')   '   ...',n_ref, ' ...   ',trim(comment), ' = ', Z(n_ref), &
                                                     '  (',h(n_ref), k(n_ref),l(n_ref),')'
   call write_info(TRIM(message_text))
   endif

   if (ALLOCATED(temp_Z))      DEALLOCATE(temp_Z)
   if (ALLOCATED(temp_F2))     DEALLOCATE(temp_F2)
   if (ALLOCATED(temp_sig_F2)) DEALLOCATE(temp_sig_F2)
   if (ALLOCATED(temp_cos))    DEALLOCATE(temp_cos)
   if (ALLOCATED(temp_H))      DEALLOCATE(temp_H)
   if (ALLOCATED(temp_code))   DEALLOCATE(temp_code)

   return
 end subroutine sort_arrays


!-------------------------------------------------------------------
 subroutine ordered_real_array(n, array, ordered, sort_type)
  implicit none
   INTEGER,               INTENT(IN)     :: n
   REAL,    DIMENSION(N), INTENT(INOUT)  :: array
   INTEGER, DIMENSION(N), INTENT(IN)     :: ordered
   CHARACTER(LEN=*),      INTENT(IN)     :: sort_type
   REAL,    DIMENSION(N)                 :: real_array
   INTEGER                               :: i

  IF(sort_type == '+') then
   do i = 1, n
    real_array(i) = array(ordered(i))
   end do
  ELSEIF(sort_type == '-') then
   do i=1, n
    real_array(i) = array(ordered(n-i))
   end do
  endif

  array = real_array

 end subroutine ordered_real_array

!-------------------------------------------------------------------
 subroutine ordered_int_array(n, array, ordered, sort_type)
  implicit none
   INTEGER,               INTENT(IN)     :: n
   integer, DIMENSION(N), INTENT(INOUT)  :: array
   INTEGER, DIMENSION(N), INTENT(IN)     :: ordered
   CHARACTER(LEN=*),      INTENT(IN)     :: sort_type
   integer, DIMENSION(N)                 :: int_array
   INTEGER                               :: i

  IF(sort_type == '+') then
   do i = 1, n
    int_array(i) = array(ordered(i))
   end do
  ELSEIF(sort_type=='-') then
   do i = 1, n
    int_array(i) = array(ordered(n-i))
   end do
  endif
  array = int_array

 end subroutine ordered_int_array


!-----------------------------------------------------------------
subroutine write_sorted_file(i, ok)
 use hkl_module
 use io_module
 use cryscalc_module, only : hkl_format, message_text, debug_proc
 implicit none
  integer, intent(in)    :: i
  logical, intent(inout) :: ok
  INTEGER                :: j, i_error

  !if(debug_proc%level_2)  call write_debug_proc_level(2, "WRITE_SORTED_FILE")

   hkl_format             = '(3I4,2F8.2,I4,6F8.5)'

   ok = .false.
   if (cos_exist) then
     !write(12, '(3I4,2F8.2,I4,6F8.5)', iostat=i_error)h(i), k(i), l(i), F2(i), sig_f2(i)/sig_coef , code(i), &
     !                                (cos_dir(i,j),j=1,6)
     write(12, trim(hkl_format), iostat=i_error)h(i), k(i), l(i), F2(i), sig_f2(i)/sig_coef , code(i), &
                                     (cos_dir(i,j),j=1,6)

     if(i_error /=0) then
      write(message_text, '(I8,a,3I4,a)') i, ' (',h(i), k(i), l(i),')'
      call write_info('  >> Wrong writing format at line '// trim(adjustl(message_text)) // '. Writing file routine is stopped.')
      return
     end if
    !IF(cif_file) then
    ! write(13, '(i4,2i5,2F9.2,I5)')   h(i), k(i), l(i), F2(i), sig_F2(i)/sig_coef, code(i)
    ! WRITE(13, '(F8.5,5F9.5)')        (cos_dir(j),j=1,6)
    !endif

   else

    IF(HKL_file%M91) then
     write(12, '(3I4,2F9.1,I4)')      h(i), k(i), l(i), F2(i), sig_f2(i)/sig_coef, TRIM(HKL_line(i))
    ELSEIF(HKL_file%M95) then
     WRITE(12, '(I6,3I4,4F7.2,2E15.6, F10.3, 2i2)', iostat=i_error) HKL_flag(i,1), h(i), k(i), l(i), cos_dir(i,1:4), &
                                                    F2(i), sig_F2(i)/sig_coef, cos_dir(i,5), HKL_flag(i,2:3)
     WRITE(12, '(a)')                               TRIM(HKL_line(i))
    else
     write(12, '(3I4,2F8.2,i4)', iostat=i_error)      h(i), k(i), l(i), F2(i), sig_f2(i)/sig_coef, code(i)
     if(i_error >0) then
      call write_info(' >> Wrong writing format. Writing file routine is stopped.')
      return
     end if
    endif
    !IF(cif_file) then
     ! WRITE(13, '(i4,2i5,2F9.2,I5)')      h(i), k(i), l(i), F2(i), sig_F2(i)/sig_coef, code(i)   ! sans les cos. dir.
    !endif
   end if

   ok = .true.
 return
end subroutine write_sorted_file

!--------------------------------------------------------------------------------------------------------------
subroutine HKL_diff_routine
 ! creation d'un fichier .HKL contenant I_HKL(1) - I_HKL(2) pour les reflexionx communes
 !
 use CRYSCALC_module, only : HKL_file_diff, HKL_SHELX_format
 use macros_module,   only : l_case
 USE IO_module,       ONLY : write_info
 implicit none
  character (len=256)              :: read_line
  integer                          :: i, i1, i_error
  character (len=32), dimension(2) :: fmt_
  integer                          :: h_1, k_1, l_1, h_2, k_2, l_2
  real                             :: F2_1, sig_1, F2_2, sig_2, F2_diff, sig_diff
  logical                          :: INT_file


  open(unit = 11, file=trim(HKL_file_diff(1)))
  open(unit = 12, file=trim(HKL_file_diff(2)))
  open(unit = 13, file="HKL_diff.hkl")

  do i=1, 2
   i1 = index(l_case(HKL_file_diff(i)), '.int')
   if(i1 /=0) then ! format .int pour FP
    read(unit = 10+i, fmt='(a)', iostat = i_error) read_line   ! TITLE
    read(unit = 10+i, fmt='(a)', iostat = i_error) read_line   ! format
    if(len_trim(read_line) /=0) read(read_line, '(a)') fmt_(i)
    read(unit = 10+i, fmt='(a)', iostat = i_error) read_line   ! WL, item1, item2
    INT_file = .true.
   else
    fmt_(i) = trim(hkl_SHELX_format)
    INT_file = .false.
   end if
  end do


  do
   read(unit = 11, fmt='(a)', iostat = i_error) read_line
   if(i_error < 0) exit
   if(len_trim(read_line) == 0) exit
   read(read_line, fmt=trim(fmt_(1))) h_1, k_1, l_1, F2_1, sig_1
   if(h_1 == 0 .and. k_1 == 0 .and. l_1 == 0) exit

   ! recherche de la reflexion dans le second fichier
   rewind(unit=12)
   if(INT_file) then
    do i = 1, 3
     read(unit = 12, fmt='(a)', iostat = i_error) read_line
     if(i_error < 0) exit
     if(len_trim(read_line) == 0) exit
    end do
   end if
   do
    read(unit = 12, fmt='(a)', iostat = i_error) read_line
    if(i_error < 0) exit
    if(len_trim(read_line) == 0) exit
    !write(*,*) ' long : ', trim(read_line), '  ', len_trim(read_line)
    read(read_line, fmt=trim(fmt_(2))) h_2, k_2, l_2, F2_2, sig_2
    if(h_2 == 0 .and. k_2 == 0 .and. l_2 == 0) exit

    if(h_1 == h_2 .and. k_1 == k_2 .and. l_1 == l_2) then
     F2_diff = F2_1 - F2_2
     sig_diff = sqrt(sig_1**2 + sig_2**2)
     write(unit=13 , fmt='(3I4, 2F8.2)') h_1, k_1, l_1, F2_diff, sig_diff
     exit
    end if
   end do

  end do

  close(unit = 13)
  close(unit = 12)
  close(unit = 11)

  call write_info('')
  call write_info('  >> HKL_diff.hkl file has been created.')
  call write_info('')

 return
end subroutine HKL_diff_routine


!----------------------------------------------------------------------
subroutine read_HKLF5_file
 use HKL_module,      only : HKL_file
 use CRYSCALC_module, only : HKL_unit, lecture_OK, message_text, clusters_out
 use MACROS_module,   only : test_file_exist
 USE IO_module,       only : write_info
 implicit none
  logical                :: file_exist
  character (len=256)    :: read_line
  integer                :: i_error
  integer                :: i, h_, k_, l_, m_, nb_domains
  integer                :: current_domain, nr_common, nr_total, nr_eff, n_n
  integer, dimension(5)  :: nr
  real                   :: F2_, sig_F2_
  REAL, DIMENSION(6)     :: cosdir
  integer                :: n_cl                     ! nombre de clusters de reflexions
  integer, dimension(10) :: cl_h, cl_k, cl_l, cl_m   ! indices associes
  real,    dimension(10) :: cl_F2, cl_sig

  file_exist = .false.

  call test_file_exist(HKL_file%NAME, file_exist, 'out')
  IF(.not. file_exist) then
   HKL_file%NAME = ''
   return
  endif

  call write_info(' ')
  call write_info('  >> HKLF5 file name: '// TRIM(HKL_file%NAME))
  call write_info(' ')

  OPEN(UNIT=HKL_unit, FILE=TRIM(HKL_file%name))
   nb_domains     = 1
   current_domain = 1
   nr        = 0
   nr_eff    = 0
   nr_common = 0
   nr_total  = 0
   cosdir = 0.
   n_n    = 0
   n_cl   = 0

   !do
   ! read(unit=HKL_unit, fmt='(a)', iostat=i_error) read_line
   ! if(i_error /=0) then
   !  lecture_ok = .false.
   !  exit
   ! end if
    !read(read_line, *) h_,k_,l_, F2_, sig_F2_, m_
    !if(h_ == 0 .and. k_ == 0 .and. l_ == 0) exit
    !nr_total = nr_total + 1

    !if(abs(m_) > nb_domains) nb_domains = abs(m_)
    !if(m_ < 0) then
    ! current_domain = -1
    ! n_n = n_n + 1
    ! if(n_n == 1) then
    !  nr_common = nr_common +1
    ! else
    !  write(*,*) " cluster  ", h_, k_, l_
    ! end if
    !end if

    !if(m_ > 0) then
    ! n_n = 0
    ! if(current_domain > 0) then
    !  nr(m_) = nr(m_) + 1
    !  current_domain = m_
    ! else
     !  current_domain = 1
    ! end if
    !end if
   !end do

loop_1: do
    n_n = 0
    read(unit=HKL_unit, fmt='(a)', iostat=i_error) read_line
    if(i_error /=0) then
     lecture_ok = .false.
     exit
    endif

    read(read_line, *) h_, k_, l_, F2_, sig_F2_, m_
    if(h_ == 0 .and. k_ == 0 .and. l_ == 0) exit
    nr_total = nr_total + 1

    if(abs(m_) > nb_domains) nb_domains = abs(m_)

    if(m_ > 0) then
     nr(m_) = nr(m_) + 1
     nr_eff = nr_eff + 1
     cycle
    end if

    n_n = n_n + 1
    cl_h(n_n)   = h_
    cl_k(n_n)   = k_
    cl_l(n_n)   = l_
    cl_m(n_n)   = m_
    cl_F2(n_n)  = F2_
    cl_sig(n_n) = sig_F2_

    do
     read(unit=HKL_unit, fmt='(a)', iostat=i_error) read_line
     if(i_error /=0) then
      lecture_ok = .false.
      exit loop_1
     endif
     read(read_line, *) h_, k_, l_, F2_, sig_F2_, m_
     if(h_ == 0 .and. k_ == 0 .and. l_ == 0) exit
     nr_total = nr_total + 1

     n_n = n_n + 1
     cl_h(n_n) = h_
     cl_k(n_n) = k_
     cl_l(n_n) = l_
     cl_m(n_n) = m_
     cl_F2(n_n)  = F2_
     cl_sig(n_n) = sig_F2_


     if(m_ > 0) then
      nr_common = nr_common +1
      nr_eff = nr_eff + 1
      exit
     end if
    end do

    if(n_n > 2) then
     n_cl = n_cl + 1
     if (clusters_out) then
      do i=1, n_n
       if(i==1) then
        write(message_text,"(a,i4,5x,4I4,2F8.2)") "        cluster #", n_cl, cl_h(i), cl_k(i), cl_l(i), cl_m(i), cl_F2(i), cl_sig(i)
       else
        write(message_text,"(a,9x,4I4,2F8.2)")    "                 ", cl_h(i), cl_k(i), cl_l(i), cl_m(i), cl_F2(i), cl_sig(i)
       end if
       call write_info(trim(message_text))
      end do
      call write_info('')
     end if
    end if

   end do loop_1



  close(unit= HKL_unit)

  write(message_text, '(a, i6)')       '   . Total number of reflexions                   = ', nr_total
  call write_info(trim(message_text))
  write(message_text, '(a, i6)')       '   . Effective number of observed reflexions      = ', nr_eff
  call write_info(trim(message_text))
  write(message_text, '(a, i6)')       '   . Number of domains                            = ', nb_domains
  call write_info(trim(message_text))

  do i=1, nb_domains
   write(message_text, '(a, i2,a,i6)') '   . Number of individual reflexions in domain ', i, ' = ', nr(i)
   call write_info(trim(message_text))
  end do
  write(message_text, '(a, i6)')       '   . Number of common reflexions                  = ', nr_common
  call write_info(trim(message_text))
  write(message_text, '(a, i6)')       '   . Number of clusters of reflexions             = ', n_cl
  call write_info(trim(message_text))

 return
end subroutine read_HKLF5_file
