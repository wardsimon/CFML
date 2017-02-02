!     Last change:  TR    6 Nov 2006    3:12 pm

subroutine read_SHELX_HKL(input_string)
 ! lecture fichier.HKL format SHELX cree par SAINT/SADABS
 USE IO_module,       ONLY  : write_info
 USE cryscalc_module, ONLY  : P4P_file_name, HKL_file_name, keyword_P4P, message_text, hkl_format, HKL_SHELX_format, CIF_unit, &
                              TWINABS_file, debug_proc
 USE HKL_module,      ONLY  : cos_exist, cos_dir
 USE macros_module,   ONLY  : u_case

 implicit none
  character (len=*), intent(in) :: input_string
  !INTEGER              :: nb_HKL_file
  INTEGER              :: i_iostat, i1, i2
  CHARACTER (LEN=256)  :: read_line
  CHARACTER (LEN=256)  :: SHELX_HKL_file
  INTEGER              :: h, k, l, code
  REAL                 :: F2, sig
  LOGICAL              :: file_exist
  real, parameter      :: eps = 0.000001
  REAL, DIMENSION(6)   :: cosdir

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_SHELX_HKL ("//trim(input_string)//")")



  ! fichier .HKL ?
  ! 1.   xx_Om.P4P
  !   1.1  xx_0m.HKL
  !   1.2  xx.HKL
  ! 2.   xx.P4P
  !   2.1  xx.HKL
  !   2.2  xx_0m.HKL

  if(len_trim(HKL_file_name) ==0) then
   i1 = INDEX(u_case(P4P_file_name), '_0M.P4P')
   if(i1 /=0) then                                  !  1.
    i2 = INDEX(u_case(P4P_file_name), '.')
    SHELX_HKL_file = P4P_file_name(1:i2)//'HKL'     !   1.1
    inquire(file = trim(SHELX_HKL_file), exist=file_exist)
    if(.not. file_exist) then
     SHELX_HKL_file = P4P_file_name(1:i1-1)//'.HKL' !   1.2
     inquire(file = trim(SHELX_HKL_file), exist=file_exist)
     if(.not. file_exist) call Stop_program(TRIM(SHELX_HKL_file))
    end if
   else
    i1 = index(u_case(P4P_file_name), '.PAP')              ! 2.
    SHELX_HKL_file = P4P_file_name(1:i1)//'HKL'            !  2.1
    inquire(file = trim(SHELX_HKL_file), exist = file_exist)
    if(.not. file_exist) then
     SHELX_HKL_file = P4P_file_name(1:i1-1)//'_0M.HKL'     !  2.2
     inquire (file = trim(SHELX_HKL_file), exist = file_exist)
     if(.not. file_exist) call Stop_program(TRIM(SHELX_HKL_file))
    end if
   endif
  else
   SHELX_HKl_file = HKL_file_name
  end if




  if(input_string(1:3) == 'cos') then
   call write_info('')
   WRITE(message_text, '(2a)') '  . HKL file: ', TRIM(SHELX_HKL_file)
   call write_info(TRIM(message_text))
  end if

  OPEN(UNIT=66, FILE=TRIM(SHELX_HKL_file))

   if(input_string(1:3) == 'cos') then
    ! cos exist ?
    cos_exist = .false.
    if(HKL_format(1:11) == '(3I4,2F8.2)') HKL_format = HKL_SHELX_format
    read(66, fmt=trim(HKL_format), iostat=i_iostat)  h,k,l, F2, sig, code, cosdir(1:6)
    if(i_iostat /=0) then
     call write_info('')
     write(message_text, '(2x,3a)') '! Error reading .hkl file. Wrong format ', trim(HKL_format), ' ?'
     call write_info(trim(message_text))
     call write_info('')
     close(unit=66)
     stop
    endif
    if (cosdir(1) < eps .and. cosdir(2) < eps .and. cosdir(3) < eps .and. cosdir(4) < eps .and. &
        cosdir(5) < eps .and. cosdir(6) < eps) then
      cos_exist = .false.
    else
      cos_exist = .true.
    end if
    close(unit=66)
    return
   end if
  !---------------------------------------

  OPEN(UNIT=CIF_unit, FILE='import.cif', POSITION = 'append')

  if(.not. cos_exist) then
   if(.not. TWINABS_file) then
    call write_CIF_section("HKL DATA from SADABS")
   else
    call write_CIF_section("HKL DATA from TWINABS")
   end if   
  else
   call write_CIF_section("HKL DATA from SAINT")
  end if
  WRITE(CIF_unit, '(a)') 'loop_'
   WRITE(CIF_unit, '(a)') '_refln_index_h'
   WRITE(CIF_unit, '(a)') '_refln_index_k'
   WRITE(CIF_unit, '(a)') '_refln_index_l'
   WRITE(CIF_unit, '(a)') '_refln_F_squared_meas'
   WRITE(CIF_unit, '(a)') '_refln_F_squared_sigma'
   WRITE(CIF_unit, '(a)') '_refln_scale_group_code'
   if(cos_exist) then
   WRITE(CIF_unit, '(a)') '_refln_nonius_incident_cos_astar'
   WRITE(CIF_unit, '(a)') '_refln_nonius_incident_cos_bstar'
   WRITE(CIF_unit, '(a)') '_refln_nonius_incident_cos_cstar'
   WRITE(CIF_unit, '(a)') '_refln_nonius_diffracted_cos_astar'
   WRITE(CIF_unit, '(a)') '_refln_nonius_diffracted_cos_bstar'
   WRITE(CIF_unit, '(a)') '_refln_nonius_diffracted_cos_cstar'
   endif
   do
    READ(66, '(a)', IOSTAT=i_iostat) read_line
    IF(i_iostat <0) exit
    if(len_trim(read_line) == 0) exit
    READ(read_line, fmt=trim(HKL_format))   h, k, l, F2, sig, code, cosdir(1:6)
    WRITE(CIF_unit, '(I4,2I5,2F9.2,I5)') h, k, l, F2, sig, code
    if(cos_exist) WRITE(CIF_unit,'(F8.5, 5(1x,F8.5))') cosdir(1:6)
   END do
  call write_CIF_section("END OF CIF FILE")
  CLOSE(UNIT=CIF_unit)
  close(UNIT=66)

 return
end subroutine read_SHELX_HKL


!-------------------------------------------------------------------------

subroutine read_RAW_file(input_string)
! lecture fichier.RAW cree par SAINT
 USE IO_module,       ONLY  : write_info
 USE HKL_module
 USE cryscalc_module, ONLY  : RAW_file_name, message_text, debug_proc
 use macros_module,   ONLY  : test_file_exist

 implicit none
  CHARACTER(len=*), INTENT(IN) :: input_string
  !INTEGER              :: nb_HKL_file
  INTEGER              :: i, i_iostat, i1
  CHARACTER (LEN=256)  :: read_line
  CHARACTER (LEN=256)  :: SHELX_HKL_file
  INTEGER              :: h_, k_, l_, code_
  REAL                 :: F2_, sig_
  REAL, dimension(6)   :: cosdir
  INTEGER              :: div_10
  real                 :: max_F2
  LOGICAL              :: file_exist


  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_RAW_FILE ("//trim(input_string)//")")

  call test_file_exist(trim(RAW_file_name), file_exist, 'out')
  if(.not. file_exist) return

  if(input_string(1:1) == '1') then
   i1 = INDEX(RAW_file_name, '.')
   SHELX_HKL_file = RAW_file_name(1:i1-1) //'_raw.HKL'
   OPEN(UNIT=67, FILE=trim(SHELX_HKL_file))
  elseif(input_string(1:1) == '2') then
   OPEN(UNIT=67, FILE='import.cif', POSITION = 'append')
  end if

  OPEN(UNIT=66, FILE=TRIM(RAW_file_name))
   n_ref = 0
   do
    READ(66, '(a)', IOSTAT=i_iostat) read_line
    IF(i_iostat <0) exit

    READ(read_line,      '(3I4)')       h_, k_, l_
    READ(read_line(13:), *)             F2_, sig_
    READ(read_line(29:), '(I4,6F8.5)')  code_, cosdir(1:6)
    !if(input_string(1:1) == '1')  write(67, '(3I4, 2F8.2,I4, 6F8.5)') h_,k_,l_, F2_, sig_, code_, cosdir(1:6)

    n_ref = n_ref+1
    h(n_ref)               = h_
    k(n_ref)               = k_
    l(n_ref)               = l_
    F2(n_ref)              = F2_
    sig_F2(n_ref)          = sig_
    code(n_ref)            = code_
    cos_dir(n_ref,1:6)     = cosdir(1:6)

   END do
  !if(input_string(1:1) == '1') CLOSE(UNIT=67)
  close(UNIT=66)

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

  if(input_string(1:1) == '1') then
   do i=1, n_ref
   ! write(67, '(3I4, 2F8.2,I4, 6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i), code(i), cos_dir(i,1:6) ! pb avec i_fort !
    write(67, '(3I4, 2F8.2,I4, 6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i), code(i), &
                                        cos_dir(i,1), cos_dir(i,2),cos_dir(i,3), cos_dir(i,4), cos_dir(i,5), cos_dir(i,6)
   end do
   call write_info('')
   call write_info(' >> '//trim(SHELX_HKl_file) // ' file (SHELX format) has been created.')
   IF(div_10 /=0) then
    write (message_text, '(a,I2)') ' >> F2 and sigmas have been divided by 10**',div_10
    call write_info(TRIM(message_text))
   endif
   call write_info('')
  elseif(input_string(1:1) == '2') then
   WRITE(67, '(a)') ''
   WRITE(67, '(a)') '#----------------------------------------------------------------------------#'
   WRITE(67, '(a)') '#                   HKL DATA from SAINT                                      #'
   WRITE(67, '(a)') '#----------------------------------------------------------------------------#'
   WRITE(67, '(a)') ''
   WRITE(67, '(a)') 'loop_'
   WRITE(67, '(a)') '_refln_index_h'
   WRITE(67, '(a)') '_refln_index_k'
   WRITE(67, '(a)') '_refln_index_l'
   WRITE(67, '(a)') '_refln_F_squared_meas'
   WRITE(67, '(a)') '_refln_F_squared_sigma'
   WRITE(67, '(a)') '_refln_scale_group_code'
   WRITE(67, '(a)') '_refln_nonius_incident_cos_astar'
   WRITE(67, '(a)') '_refln_nonius_incident_cos_bstar'
   WRITE(67, '(a)') '_refln_nonius_incident_cos_cstar'
   WRITE(67, '(a)') '_refln_nonius_diffracted_cos_astar'
   WRITE(67, '(a)') '_refln_nonius_diffracted_cos_bstar'
   WRITE(67, '(a)') '_refln_nonius_diffracted_cos_cstar'

   do i=1, n_ref
    WRITE(67, '(I4,2I5,2F9.2,I5)') h(i), k(i), l(i), F2(i), sig_F2(i), code(i)
    !WRITE(67, '(F8.5, 5(1x,F8.5))') cos_dir(i,1:6)  ! >> pb avec ifort !
    WRITE(67, '(F8.5, 5(1x,F8.5))') cos_dir(i,1), cos_dir(i,2),cos_dir(i,3), cos_dir(i,4), cos_dir(i,5), cos_dir(i,6)
   end do
  end if
  call write_CIF_section("END OF CIF FILE")
  close(unit=67)
 return
end subroutine read_RAW_file


!----------------------------------------------------------------------------------
subroutine Stop_program(SHELX_HKL_file)
 USE IO_module,       ONLY  : write_info
 USE cryscalc_module, ONLY  : message_text

 implicit none
  character (len=*), intent(in)      :: SHELX_HKL_file

    call write_info('')
    WRITE(message_text, '(2x,3a)') '!! ',TRIM(SHELX_HKL_file), ' does not exist. '
    call write_info(TRIM(message_text)// ' Program is stopped !!')
    call write_info('')
    stop

 return
end subroutine Stop_program
