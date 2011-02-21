!     Last change:  TR   23 Apr 2007    7:00 pm
!
subroutine read_and_sort_HKL(input_string)
 USE cryscal_module, ONLY      : nb_sort,  sort_type,  cut_off,             &
                                 nb_shell, shell_type, shell_plot,          &
                                 keyword_CELL, keyword_WAVE, keyword_FILE,  &
                                 HKL_unit, message_text
 USE HKL_module,     ONLY      : HKL_file
 USE macros_module,  ONLY      : test_file_exist
 USE IO_module,      ONLY      : write_info
 USE hkl_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_string
  INTEGER                              :: i
  LOGICAL                              :: file_exist
  LOGICAL                              :: lecture_ok
  CHARACTER (LEN=4)                    :: ext_string


  IF(.NOT. keyword_FILE) then
   call write_info('')
   call write_info('  >>> FILE is a mandatory keyword !')
   call write_info('')
   return
  endif
  IF(HKL_file%M91) THEN
   ext_string = '.m91'
  ELSEIF(HKL_file%M95) then
   ext_string = '.m95'
  else
   ext_string = '.hkl'
  endif

  call test_file_exist(HKL_file%NAME, file_exist)
  IF(.not. file_exist) then
   HKL_file%NAME = ''
   return
  endif

  call write_info(' ')
  call write_info('   >> HKL file name: '// TRIM(HKL_file%NAME))
  call write_info(' ')
  i = INDEX(HKL_file%NAME, ".")
  if (i/=0) then
   IF(input_string(1:4) == 'sort') then

    HKL_file%d     = TRIM(HKL_file%NAME(1:i-1))//'_sort_d'//TRIM(ext_string)
    HKL_file%stl   = TRIM(HKL_file%NAME(1:i-1))//'_sort_stl'//TRIM(ext_string)
    HKL_file%theta = TRIM(HKL_file%NAME(1:i-1))//'_sort_theta'//TRIM(ext_string)
    HKL_file%I     = TRIM(HKL_file%NAME(1:i-1))//'_sort_i'//TRIM(ext_string)
    HKL_file%Isig  = TRIM(HKL_file%NAME(1:i-1))//'_sort_isig'//TRIM(ext_string)
   ELSEIF(input_string(1:5) == 'shell') then
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
      call write_info(' ')
      call write_info('  >>> Create '// TRIM(HKL_file%d)// ' file ...')
      call write_info(' ')

      OPEN(UNIT=12, FILE=TRIM(HKL_file%d))
      IF(.not. HKL_data_known) call read_data_file(lecture_ok)
      IF(.NOT. lecture_ok) return
      call classement(i,'shell_d')
      CLOSE(UNIT=12)

     case ('stl')
      if (.not. keyword_CELL) then
       call write_info(' ')
       call write_info('  >>> cell parameters are mandatory to calculate SinTheta/lambda <<<')
       call write_info(' ')
       cycle
      ELSEIF(.NOT. keyword_WAVE) then
       call write_info(' ')
       call write_info('  >>> wavelength is a mandatory parameter to calculate SinTheta/lambda <<<')
       call write_info(' ')
       cycle
      end if

      call write_info(' ')
      call write_info('  >>> Create '// TRIM(HKL_file%stl)// ' file ...')
      call write_info(' ')

      OPEN(UNIT=12, FILE=TRIM(HKL_file%stl))
      IF(.not. HKL_data_known) call read_data_file(lecture_ok)
      IF(.NOT. lecture_ok) return
      call classement(i,'shell_stl')
      CLOSE(UNIT=12)

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
       call write_info(' ')
       call write_info('  >>> Create '// TRIM(HKL_file%theta)// ' file ...')
       call write_info(' ')

       OPEN(UNIT=12, FILE=TRIM(HKL_file%theta))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok)
       IF(.NOT. lecture_ok) return
       call classement(i,'shell_theta')
       close (UNIT=12)

     case ('int')
       call write_info(' ')
       call write_info('  >>> Create '// TRIM(HKL_file%I)// ' file ...')

       OPEN(UNIT=12, FILE=TRIM(HKL_file%i))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok)
       IF(.NOT. lecture_ok) return
       call classement(i,'shell_int')
       CLOSE(UNIT=12)

     case ('isig')
       call write_info(' ')
       call write_info('  >>>Create '// TRIM(HKL_file%Isig)// ' file ...')

       OPEN(UNIT=12, FILE=TRIM(HKL_file%isig))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok)
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
   IF(.not. HKL_data_known) call read_data_file(lecture_ok)
   IF(.NOT. lecture_ok) return
   call classement(0, 'no_sort')

  else
   do i=1, nb_sort
    OPEN(UNIT=HKL_unit, FILE=TRIM(HKL_file%NAME))

    select case (sort_type(i))
     case ('d')
       call write_info(' ')
       call write_info('  >>> SORT HKL file in increasing d(A) and create '// TRIM(HKL_file%d)// ' file ...')
       call write_info(' ')

       OPEN(UNIT=12, FILE=TRIM(HKL_file%d))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok)
       IF(.NOT. lecture_ok) return
       call classement(i,'sort_d')
       CLOSE(UNIT=12)

     case ('stl')
       if (.not. keyword_CELL) then
        call write_info(' ')
        call write_info('  >>> cell parameters are mandatory to calculate SinTheta/lambda <<<')
        call write_info(' ')
        cycle
       ELSEIF(.NOT. keyword_WAVE) then
        call write_info(' ')
        call write_info('  >>> wavelength is a mandatory parameter to calculate SinTheta/lambda <<<')
        call write_info(' ')
        cycle
       end if
       call write_info(' ')
       call write_info('  >>> SORT HKL file in increasing SinTheta/lambda and create '//TRIM(HKL_file%stl)// ' file ...')
       call write_info(' ')

       OPEN(UNIT=12, FILE=TRIM(HKL_file%stl))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok)
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
       call write_info(' ')
       call write_info('  >>> SORT HKL file in increasing Theta and create '// TRIM(HKL_file%theta)// ' file ...')
       call write_info(' ')

       OPEN(UNIT=12, FILE=TRIM(HKL_file%theta))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok)
       IF(.NOT. lecture_ok) return
       call classement(i,'sort_theta')
       close (UNIT=12)

     case ('int')
       call write_info(' ')
       call write_info('  >>> SORT HKL file in decreasing intensities and create '// TRIM(HKL_file%I)// ' file ...')

       OPEN(UNIT=12, FILE=TRIM(HKL_file%i))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok)
       IF(.NOT. lecture_ok) return
       call classement(i,'sort_int')
       CLOSE(UNIT=12)

     case ('isig')
       call write_info(' ')
       call write_info('  >>> SORT HKL file in decreasing I/sig and create '// TRIM(HKL_file%Isig)// ' file ...')

       OPEN(UNIT=12, FILE=TRIM(HKL_file%isig))
       IF(.not. HKL_data_known) call read_data_file(lecture_ok)
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

subroutine read_data_file(ok)
 USE cryscal_module,    ONLY         : cut_off, unit_cell, wave => wavelength,              &
                                       keyword_CELL, keyword_WAVE, message_text, HKL_unit,  &
                                       PGF_data, PGF_file, known_theta, keyword_create_CIF, &
                                       crystal_system 
 USE HKL_module
 USE IO_module,         ONLY          : write_info
 USE CFML_Math_General, ONLY          : sort

 implicit none
 LOGICAL, INTENT(OUT)              :: ok
 integer                           :: i, i1, ier

 integer                           :: h_, k_, l_, m_, code_
 real                              :: F2_, sig_F2_
 REAL, DIMENSION(6)                :: cosdir
 LOGICAL                           :: HKL_ok
 real, parameter                   :: eps = 0.000001
 CHARACTER (LEN=256)               :: read_line
 CHARACTER (LEN=256)               :: M91_string
 REAL,    DIMENSION(5)             :: M95_r
 INTEGER, DIMENSION(3)             :: M95_i

 INTEGER                           :: div_10
 REAL                              :: max_F2
 
 logical                           :: QVEC_type


  !lecture fichier .HKL format SHELX ou fichier .CIF
  ! + creation fichier avec sintheta/l dans la 1ere colonne


  call write_info(' ')
  call write_info('  ... Reading data ...')
  call write_info(' ')
  ok = .false.

  n_ref     = 0
  n_ref_2   = 0
  n_ref_3   = 0
  n_ref_neg = 0
  cosdir    = 0.
  HKL_file%HKL = ''

  IF(HKL_file%SHELX) THEN     ! lecture d'un fichier.HKl format SHELX
   do
    READ(HKL_unit,'(3I4,2F8.2,I4,6F8.5)', iostat=ier) h_,k_,l_, F2_, sig_F2_,code_,cosdir(1:6)
    if (ier < 0) exit ! fin du fichier atteinte
    IF(h_==0 .AND. k_==0 .AND. l_==0) exit

    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)

    if (n_ref ==1) then
     if (cosdir(1) < eps .and.  cosdir(2) < eps .and. cosdir(3) < eps .and. cosdir(4) < eps .and. &
	     cosdir(5) < eps .and. cosdir(6) < eps) then
      cos_exist = .false.
     else
      cos_exist = .true.
     end if
    endif
   end do
   ok = .true.
   HKL_data_known = .true.
   HKL_file%HKL   = HKL_file%NAME

  ELSEIF(HKL_file%CIF) then ! lecture fichier import.CIF

   call get_wavelength_from_CIF_file(HKL_unit, wave)
   keyword_WAVE = .true.
   call get_cell_parameters_from_CIF_file(HKL_unit, unit_cell%param, unit_cell%param_esd )
   call volume_calculation('no_out')
   !call volume_calculation('out')
   call get_crystal_system_from_CIF_file(HKL_unit, crystal_system)
   keyword_CELL = .true.
   call read_lines_before_hkl(HKL_unit, cos_exist)
   IF(.not. cos_exist) then
    HKL_file%HKL = 'import.HKL'
   else
    HKL_file%HKL = 'import_cos.HKL'
   endif
   open (UNIT=31, FILE=TRIM(HKL_file%HKL))

   do
    if (.not. cos_exist) then
     READ(HKL_unit, '(i4,2i5,2F9.2,I5)', IOSTAT=ier) h_,k_,l_,F2_, sig_F2_, code_  ! sans les cos. dir.
      if (ier < 0 ) exit
      IF(h_==0 .AND. k_==0 .AND. l_==0) exit
     WRITE(31, '(3i4,2F8.2,I4)')h_,k_,l_,F2_, sig_F2_, code_
    else
     READ(HKL_unit,'(i4,2i5,2F9.2,I5)', IOSTAT=ier) h_,k_,l_,F2_, sig_F2_, code_
      if (ier /=0 ) exit                                                         ! avec les cos. dir.
      IF(h_==0 .AND. k_==0 .AND. l_==0) exit
     READ(HKL_unit,'(F8.5,5F9.5)',IOSTAT=ier) cosdir(1:6)
      if (ier /=0 ) exit
     write(31,'(3I4,2F8.2,I4,6F8.5)', iostat=ier) h_,k_,l_, F2_, sig_F2_,code_, cosdir(1:6)
    endif

    call check_HKL(h_, k_, l_, F2_, sig_F2_, code_, cosdir)

   END do
   ok = .true.
   HKL_data_known = .true.
   CLOSE(UNIT=31)

  ELSEIF(HKL_file%FINAL_y) then
   call get_wavelength_from_FINAL_Y_file(HKL_unit, wave)
   keyword_WAVE = .true.
   call get_cell_parameters_from_FINAL_Y_file(HKL_unit, unit_cell%param)
   call volume_calculation('no_out')
   keyword_CELL = .true.
   
   call get_modulation_vect(HKL_unit, qvec_type)

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
     call read_mod_reflexion_from_FINAL_Y(HKL_unit, h_, k_, l_, m_, F2_, sig_F2_, HKL_ok)
     IF(.NOT. HKL_ok) exit
     IF(h_==0 .AND. k_==0 .AND. l_==0) exit

     call check_HKLM(h_, k_, l_, m_, F2_, sig_F2_, code_, cosdir)
     WRITE(31, '(4i4,2F8.2)') h_,k_,l_,m_, F2_, sig_F2_
    END do

   endif 
   ok = .true.
   HKL_data_known = .true.
   CLOSE(UNIT=31)

  ELSEIF(HKL_file%RAW) then

   i1 = INDEX(HKL_file%NAME, ".")
   if (i1/=0) then
    HKL_file%HKL = TRIM(HKL_file%NAME(1:i1-1))//'_raw.hkl'
   else
    HKL_file%HKL = TRIM(HKL_file%NAME(1:))//'_raw.hkl'
   endif
   OPEN (UNIT=31, FILE=TRIM(HKL_file%HKL))

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
   do
    max_F2 = MAXVAL(F2(1:n_ref))
    if (max_F2 > 99999.99) then
     F2(1:n_ref)     = F2(1:n_ref)     / 10
     sig_F2(1:n_ref) = sig_F2(1:n_ref) / 10
     div_10 = div_10 + 1
    else
     exit
    end if
   END do

   do i=1, n_ref
    write(31, '(3I4, 2F8.2,I4, 6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i), code(i), cos_dir(i, 1:6)
   END do
   CLOSE(UNIT=31)

  ELSEIF(HKL_file%M91) then
   i1 = INDEX(HKL_file%NAME, ".")
   if (i1/=0) then
    HKL_file%HKL = TRIM(HKL_file%NAME(1:i1-1))//'_91.hkl'
   else
    HKL_file%HKL = TRIM(HKL_file%NAME(1:))//'_m91.hkl'
   endif
   OPEN (UNIT=31, FILE=TRIM(HKL_file%HKL))

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
    write( 1, '(3I4, 2F8.2)') h(i),k(i),l(i), F2(i)/10., sig_F2(i)/10.
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

  else
   call write_info('')
   call write_info('  >>>> Unknown data file format <<<<')
   call write_info('')
   CLOSE(UNIT = HKL_unit)
   return
  endif
  CLOSE(UNIT=HKL_unit)




   call write_info('')
   write(message_text,'(a,I8)')  '  >> Total number of reflections:          ', n_ref
   call write_info(TRIM(message_text))
   IF(n_ref == 0) then
    ok = .false.
    return
   endif

   write(message_text,'(a,I8)')  '  >> Number of reflections with I > 2sig.: ', n_ref_2
   call write_info(TRIM(message_text))
   write(message_text,'(a,I8)')  '  >> Number of reflections with I > 3sig.: ', n_ref_3
   call write_info(TRIM(message_text))
   write(message_text,'(a,I8)')  '  >> Number of reflections with I < 0.   : ', n_ref_neg
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
    WRITE(message_text,'(a,I2)')            '                           . number of scans: ', MAXVAL(code(1:n_ref))
    call write_info(TRIM(message_text))
    call write_info('')
   endif

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

   IF(LEN_TRIM(HKL_file%HKL) /=0) then
    call write_info(' ')
    write (message_text, '(3a)') '   ==> ', trim(HKL_file%HKL), ' FILE has been created (SHELX format)'
    call write_info(TRIM(message_text))
    IF(HKL_file%RAW .AND. div_10 /=0) then
     write (message_text, '(a,I2)') ' >> F2 and sigmas have been divided by 10**',div_10
     call write_info(TRIM(message_text))
    endif
    call write_info(' ')
   endif


   IF(keyword_create_CIF) then
    call write_CIF_file('DATA_LIMITS')
    call write_CIF_file('HKL_DATA')
   endif
   

 return
end subroutine read_data_file

!----------------------------------------------------------------------------------------------
subroutine check_HKL(H_h, H_k, H_l, H_F2, H_sig, code_, H_cos)
 USE HKL_module
 USE cryscal_module, ONLY : keyword_cell, keyword_WAVE, known_theta,  unit_cell, wave => wavelength

 INTEGER, INTENT(IN)            :: H_h, H_k, H_l
 REAL,    INTENT(INOUT)         :: H_F2, H_sig
 INTEGER, INTENT(IN)            :: code_
 REAL, DIMENSION(6), INTENT(IN) :: H_cos
 REAL                           :: Q, Z
 real, parameter                :: eps = 0.000001
 real, parameter                :: pi=3.1415926535897932

    n_ref = n_ref + 1
    call check_nref(n_ref)
    IF(H_F2 < 0.00001) n_ref_neg = n_ref_neg + 1
    IF(.not. HKL_file%read_NEG) call check_F2(H_F2, H_sig)

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
    sig_F2(n_ref)          = H_sig
    I_sigma(n_ref)         = H_F2/ H_sig
    I_(n_ref)              = H_F2
    code(n_ref)            = code_
    cos_dir(n_ref,1:6)     = H_cos(1:6)

    WRITE(HKL_string(n_ref), '(a,3I4,a)') '(', H_h, H_k, H_l, ')'

  return
end subroutine check_HKL
!----------------------------------------------------------------------------------------------
subroutine check_HKLM(H_h, H_k, H_l, H_m, H_F2, H_sig, code_, H_cos)
 USE HKL_module
 USE cryscal_module, ONLY : keyword_cell, keyword_WAVE, known_theta,  unit_cell, wave => wavelength, Qvec

 INTEGER, INTENT(IN)            :: H_h, H_k, H_l, H_m
 REAL,    INTENT(INOUT)         :: H_F2, H_sig
 INTEGER, INTENT(IN)            :: code_
 REAL, DIMENSION(6), INTENT(IN) :: H_cos
 REAL                           :: Q, Z
 real, parameter                :: eps = 0.000001
 real, parameter                :: pi=3.1415926535897932

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
    sig_F2(n_ref)          = H_sig
    I_sigma(n_ref)         = H_F2/ H_sig
    I_(n_ref)              = H_F2
    code(n_ref)            = code_
    cos_dir(n_ref,1:6)     = H_cos(1:6)

    WRITE(HKL_string(n_ref), '(a,4I4.0,a)') '(', H_h, H_k, H_l, H_m, ')'

  return
end subroutine check_HKLM

!----------------------------------------------------------------------------------------------------------------
subroutine check_F2(F2_, sig_F2_)
 REAL, INTENT(INOUT) :: F2_, sig_F2_

  if(F2_     < 0.0001)  F2_= 0.0001
  if(sig_F2_ < 0.00001) sig_F2_=sqrt(abs(F2_))

 return
END subroutine check_F2

!----------------------------------------------------------------------------------------------------------------
subroutine check_nref(n_ref)
 USE HKL_module, ONLY : max_ref
 USE IO_module,  ONLY : write_info
 INTEGER, INTENT(IN)  :: n_ref

    IF(n_ref > max_ref) then
     call write_info('')
     call write_info('  >> Number of reflections larger than arrays dimensions.')
     call write_info('     Please inform the author of CRYSCAL program (T.R./ CDIFX Rennes):')
     call write_info('      thierry.roisnel@univ-rennes1.fr')
     call write_info('     Program will be stopped.')
     call write_info(' ')
     stop
    endif
  return
end subroutine check_nref
!----------------------------------------------------------------------------------------------------------------
subroutine classement(n, input_string)
 USE cryscal_module, ONLY          : cut_off, unit_cell, wave => wavelength, sort_plot, shell_plot,    &
                                     sort_out, sort_out_n,                                             &
                                     keyword_CELL, keyword_WAVE, message_text,                         &
                                     nb_shell, shell_arg_min_max, shell_min, shell_max,                &
                                     PGF_data, PGF_file, known_theta, winplotr_exe, input_line
 USE CFML_Math_General,       ONLY          : sort
 USE hkl_module
 USE IO_module
 
 implicit none
 INTEGER,           INTENT(IN)     :: n            ! numero du passage dans la routine
 character (len=*), intent(in)     :: input_string
 INTEGER                           :: long_input_string
 integer                           :: i
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
 real, parameter                   :: eps = 0.000001
 INTEGER,  ALLOCATABLE, DIMENSION(:) :: ordered_array


  long_input_string = LEN_TRIM(input_string)


  IF(input_string(1:5) == 'sort_' .or. input_string(1:7) == 'no_sort') then
   if (n==0 .or. n==1) then
    call test_ratio_I_sig()
    !call statistics_on_E2()

    IF(keyword_CELL) then
     call stat_repartition_sintheta_lambda()    ! ok

     !call stat_repartition_d_hkl(F2, sig_F2, d_hkl, n_ref)                       ! nok
     IF(keyword_WAVE) then
      call stat_repartition_theta()                    ! ok
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
   ALLOCATE(ordered_array(n_ref))

   select case (input_string)
     case ('sort_stl')
       call write_info('')
       call write_info('  ... Sorting hkl in increasing SinTheta...')
       call write_info('')

       ! tri par ordre croissant des X
       call sort(sintheta_lambda, n_ref, ordered_array )     ! math_gen module
       call sort_arrays(input_string, '+', ordered_array, sintheta_lambda)

     case ('sort_d')
       call write_info('')
       call write_info(' ... Sorting hkl in increasing d_hkl...')
       call write_info('')

       call sort(d_hkl, n_ref, ordered_array)
       call sort_arrays(input_string,'+', ordered_array, d_hkl)
       call stat_repartition_d_hkl(F2, sig_F2, d_hkl, n_ref)

     case ('sort_theta')
       call write_info('')
       call write_info(' ... Sorting hkl in increasing theta...')
       call sort(theta_hkl, n_ref, ordered_array)
       call sort_arrays(input_string,'+', ordered_array, theta_hkl)


     case ('sort_int')
       call write_info('')
       call write_info('  ... Sorting hkl in decreasing intensities ...')
       call write_info('')

      call sort(I_, n_ref, ordered_array)
      call sort_arrays(input_string,'-', ordered_array, I_)

      F2(1:n_ref) = I_(1:n_ref)

     CASE ('sort_isig')
       call write_info('')
       call write_info('  ... Sorting hkl in decreasing intensities/sigma ...')
       call write_info('')

      call sort(I_sigma, n_ref, ordered_array)
      call sort_arrays(input_string,'-', ordered_array, I_sigma)

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


   IF(keyword_CELL) then
    IF(input_string(1:5) == 'sort_' .or. input_string(1:9) == 'shell_stl') then
     call write_info('')
     write(message_text,'(a,2(2x,F8.5))') '  Experimental (sinTheta/lambda) shell: ',   min_stl, max_stl
     call write_info(TRIM(message_text))
     call write_info('')

     if (input_string(long_input_string-2:long_input_string) == 'stl' ) then
      IF(.NOT. shell_arg_min_max) then
       call write_info(' > Enter expected (sinTheta/lambda) min. and max. values [0. 0.: exp. values]: ')
       call read_input_line(input_line)
       read(input_line,*)  expected_min_stl, expected_max_stl

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

    IF(input_string(1:5) == 'sort_' .or. input_string(1:7) == 'shell_d') then
     write(message_text,'(a,2(2x,F8.5))') '  Experimental d_hkl(A) shell:          ',   min_d_hkl, max_d_hkl
     call write_info(TRIM(message_text))
     call write_info('')
     if (input_string(long_input_string:long_input_string) == 'd' ) then
      IF(.NOT. shell_arg_min_max) then
       call write_info('  > Enter expected d_hkl min. and max. values [0. 0.: exp. values]: ')
       call read_input_line(input_line)
       read(input_line,*)  expected_min_d_hkl, expected_max_d_hkl

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

    IF(input_string(1:5) == 'sort_' .or. input_string(1:11) == 'shell_theta') then
     write(message_text,'(a,2(2x,F8.3))') '  Experimental theta(deg) shell:        ',   min_theta, max_theta
     call write_info(TRIM(message_text))
     call write_info('')
     if (input_string(long_input_string-4:long_input_string) == 'theta' ) then
      IF(.NOT. shell_arg_min_max) then
       call write_info(' > Enter expected theta min. and max. values [0. 0.: exp. values]: ')
       call read_input_line(input_line)
       read(input_line,*)  expected_min_theta, expected_max_theta

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


   IF(input_string(1:5) == 'sort_' .or. input_string(1:9) == 'shell_int') then
    IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
     write(message_text,'(a,2F8.2)')      '  Experimental Intensities shell:        ', intensity_min, intensity_max
    else
     write(message_text,'(a,2F9.1)')      '  Experimental Intensities shell:        ', intensity_min, intensity_max
    endif
    call write_info(TRIM(message_text))
    call write_info('')

    if (input_string(long_input_string-2:long_input_string) == 'int' ) then
     IF(.NOT. shell_arg_min_max) then
      call write_info(' > Enter expected Intensity cut-off min. and max. values [0. 0.: exp. values]: ')
      call read_input_line(input_line)
      read(input_line,*)  expected_min_intensity, expected_max_intensity
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

   IF(input_string(1:5) == 'sort_' .or. input_string(1:10) == 'shell_isig') then
    IF(.not. HKL_file%M91 .and. .not. HKL_file%M95) then
     write(message_text,'(a,2F8.2)')       '  Experimental (I/sig) shell:           ', min_cut_off, max_cut_off
    else
     write(message_text,'(a,2F9.1)')       '  Experimental (I/sig) shell:           ', min_cut_off, max_cut_off
    endif
    call write_info(TRIM(message_text))
    call write_info('')

    IF(input_string(long_input_string-3:long_input_string) == 'isig' ) then
     IF(.NOT. shell_arg_min_max) then
      call write_info(' > Enter expected I/sig cut-off min. and max. values [0. 0.: exp. values]: ')
      call read_input_line(input_line)
      read(input_line,*)  expected_min_cut_off, expected_max_cut_off
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

    write(message_text,*) ' > Divider sigma coefficient: ' , sig_coef
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

    call write_info('')
    select case (input_string)
     case ('sort_d', 'shell_d')
      n_ref_eff      = 0
      do i=1, n_ref
       if(d_hkl(i) >= expected_min_d_hkl .and. d_hkl(i) <= expected_max_d_hkl)  then
        call write_sorted_file(i)
        n_ref_eff = n_ref_eff + 1        
       	pgf_data(n_ref_eff)%X = d_hkl(i)
       	pgf_data(n_ref_eff)%Y = F2(i)
       	pgf_data(n_ref_eff)%h = h(i)
       	pgf_data(n_ref_eff)%k = k(i)
       	pgf_data(n_ref_eff)%l = l(i)
        WRITE(pgf_data(n_ref_eff)%string, '(a,3I4,a)') '(', h(i), k(i), l(i), ')'
        
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
       endif
      end do
      IF(sort_out(n)) call write_info(' ')
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
      n_ref_eff = 0
      do i=1, n_ref
       if(sinTheta_lambda(i) >= expected_min_stl .and. sinTheta_lambda(i) <= expected_max_stl) then
        call write_sorted_file(i)
        n_ref_eff = n_ref_eff + 1
       	pgf_data(n_ref_eff)%X = sinTheta_lambda(i)
       	pgf_data(n_ref_eff)%Y = F2(i)
       	pgf_data(n_ref_eff)%h = h(i)
       	pgf_data(n_ref_eff)%k = k(i)
       	pgf_data(n_ref_eff)%l = l(i)
        WRITE(pgf_data(n_ref_eff)%string, '(a,3I4,a)') '(', h(i), k(i), l(i), ')'
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
       endif
      end do
      IF(sort_out(n)) call write_info(' ')
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
      n_ref_eff = 0
      do i=1, n_ref
       if(Theta_hkl(i) >= expected_min_theta .and. theta_hkl(i) <= expected_max_theta) then
        call write_sorted_file(i)
        n_ref_eff = n_ref_eff + 1
       	pgf_data(n_ref_eff)%X = Theta_hkl(i)
       	pgf_data(n_ref_eff)%Y = F2(i)
       	pgf_data(n_ref_eff)%h = h(i)
       	pgf_data(n_ref_eff)%k = k(i)
       	pgf_data(n_ref_eff)%l = l(i)
        WRITE(pgf_data(n_ref_eff)%string, '(a,3I4,a)') '(', h(i), k(i), l(i), ')'
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
       END if
      end do
      IF(sort_out(n)) call write_info(' ')
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
      do i=1, n_ref
       if(F2(i) >= expected_min_intensity .and. F2(i)  <= expected_max_intensity)  call write_sorted_file(i)       	
      end do
      call write_info(' ')
      call write_info('  >>> '//TRIM(HKL_file%I)// ' file has been created.')

     case ('sort_isig', 'shell_isig')
       do i=1, n_ref
       if(F2(i)/sig_F2(i) >= expected_min_cut_off .and. F2(i)/sig_F2(i)  <= expected_max_cut_off) &
	      call write_sorted_file(i)       	
      end do
      call write_info(' ')
      call write_info('  >>> '//TRIM(HKL_file%Isig)// ' file has been created.')

    end select






 RETURN
end subroutine classement


!---------------------------------------------------------------------------------------------

subroutine calcul_Q(h,k,l,cell,Q)
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
  USE cryscal_module, ONLY : message_text
  USE IO_module,      ONLY : write_info

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

    if (1000*int(i/1000) == i ) then
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

   USE cryscal_module, ONLY : message_text
 USE IO_module,        ONLY : write_info

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

    if (1000*int(i/1000) == i) then
     write(message_text,'(a,i6,a,a,a,F10.5)') '   ...',i, ' ...   ',trim(comment), ' = ', Z(i)
     call write_info(TRIM(message_text))
    endif

   END do

   write(message_text,'(a,i6,a,a,a,F10.5)') '   ...',n, ' ...   ',trim(comment), ' = ', Z(n)
   call write_info(TRIM(message_text))

   return
 end subroutine classement_par_ordre_decroissant

!--------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

 subroutine sort_arrays(sort_type, sort_sign, ordered, Z)

 ! classement d'un tableau X,Y... par ordre croissant des X
   USE cryscal_module, ONLY : message_text
 USE IO_module,        ONLY : write_info
   USE HKL_module,     ONLY : n_ref, h, k,l, F2, sig_F2, code, cos_dir, HKL_flag

  implicit none
   character  (len=*),                intent(in)              :: sort_type
   CHARACTER  (LEN=*),                INTENT(IN)              :: sort_sign
   INTEGER,       DIMENSION(N_ref),   INTENT(IN)              :: ordered
   REAL,          DIMENSION(N_ref),   INTENT(INOUT)           :: Z

   !local variables
   character (len=32)             :: comment
   integer                        :: i, rang
   REAL,    allocatable, DIMENSION(:)      :: temp_Z, temp_F2, temp_sig_F2
   REAL,    allocatable, DIMENSION(:,:)    :: temp_cos
   integer, allocatable, DIMENSION(:,:)    :: temp_H
   INTEGER, allocatable, DIMENSION(:)      :: temp_code

   if (ALLOCATED(temp_Z))      DEALLOCATE(temp_Z)
   if (ALLOCATED(temp_F2))     DEALLOCATE(temp_F2)
   if (ALLOCATED(temp_sig_F2)) DEALLOCATE(temp_sig_F2)
   if (ALLOCATED(temp_cos))    DEALLOCATE(temp_cos)
   if (ALLOCATED(temp_H))      DEALLOCATE(temp_H)
   if (ALLOCATED(temp_code))   DEALLOCATE(temp_code)

   ALLOCATE(temp_Z(500000))
   ALLOCATE(temp_F2(500000))
   ALLOCATE(temp_sig_F2(500000))
   ALLOCATE(temp_cos(500000,6))
   ALLOCATE(temp_code(500000))
   ALLOCATE(temp_H(500000,3))




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

         if (1000*int(i/1000) == i ) then
           write(message_text,'(a,i6,a,a,a,F10.5,a,3I4,a)') '   ...',i, ' ...   ',trim(comment), ' = ', temp_Z(i), &
                                                           '  (',temp_H(i,1:3),')'
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
         if (1000*int(i/1000) == i ) then
           write(message_text,'(a,i6,a,a,a,F10.5,a,3I4,a)') '   ...',i, ' ...   ',trim(comment), ' = ', temp_Z(i), &
                                                           '  (',temp_H(i,1:3),')'
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



   write(message_text,'(a,i6,a,a,a,F10.5,a,3I4,a)')   '   ...',n_ref, ' ...   ',trim(comment), ' = ', Z(n_ref), &
                                                     '  (',h(n_ref), k(n_ref),l(n_ref),')'
   call write_info(TRIM(message_text))

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
subroutine write_sorted_file(i)
 use hkl_module
 implicit none
  integer, intent(in) :: i
  INTEGER             :: j


        if (cos_exist) then
         write(12, '(3I4,2F8.2,I4,6F8.5)')h(i), k(i), l(i), F2(i), sig_f2(i)/sig_coef , code(i), &
                                         (cos_dir(i,j),j=1,6)
         !IF(cif_file) then
         ! write(13, '(i4,2i5,2F9.2,I5)')   h(i), k(i), l(i), F2(i), sig_F2(i)/sig_coef, code(i)
         ! WRITE(13, '(F8.5,5F9.5)')        (cos_dir(j),j=1,6)
         !endif

        else

         IF(HKL_file%M91) then
          write(12, '(3I4,2F9.1,I4)')      h(i), k(i), l(i), F2(i), sig_f2(i)/sig_coef, TRIM(HKL_line(i))
         ELSEIF(HKL_file%M95) then
          WRITE(12, '(I6,3I4,4F7.2,2E15.6, F10.3, 2i2)') HKL_flag(i,1), h(i), k(i), l(i), cos_dir(i,1:4), &
		                                                 F2(i), sig_F2(i)/sig_coef, cos_dir(i,5), HKL_flag(i,2:3)
          WRITE(12, '(a)')                               TRIM(HKL_line(i))
         else
          write(12, '(3I4,2F8.2,i4)')      h(i), k(i), l(i), F2(i), sig_f2(i)/sig_coef, code(i)
         endif
         !IF(cif_file) then
         ! WRITE(13, '(i4,2i5,2F9.2,I5)')      h(i), k(i), l(i), F2(i), sig_F2(i)/sig_coef, code(i)   ! sans les cos. dir.
         !endif
        end if

 return
end subroutine write_sorted_file
