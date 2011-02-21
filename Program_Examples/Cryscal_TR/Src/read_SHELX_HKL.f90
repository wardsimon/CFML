!     Last change:  TR    6 Nov 2006    3:12 pm
subroutine read_SHELX_HKL()
 ! lecture fichier.HKL format SHELX cree par SADABS
 USE IO_module,      ONLY  : write_info
 USE cryscal_module, ONLY  : P4P_file_name, keyword_P4P, message_text

 implicit none
  !INTEGER              :: nb_HKL_file
  INTEGER              :: i_iostat, i1
  CHARACTER (LEN=256)  :: read_line
  CHARACTER (LEN=256)  :: SHELX_HKL_file
  INTEGER              :: h, k, l, code
  REAL                 :: F2, sig
  LOGICAL              :: file_exist


  i1 = INDEX(P4P_file_name, '.')
  SHELX_HKL_file = P4P_file_name(1:i1) //'HKL'
  inquire (FILE=TRIM(SHELX_HKL_file), EXIST= file_exist)
  IF(.NOT. file_exist) then
   call write_info('')
   WRITE(message_text, '(2x,3a)') '!! ',TRIM(SHELX_HKL_file), ' does not exist  !! '
   call write_info(TRIM(message_text))
   if(keyword_P4P) pause
   return
  endif


  !call system('dir *.hkl /B > hkl.lst')
  !OPEN(UNIT=66, FILE='hkl.lst')
  ! nb_HKL_file = 0
  ! do
  !  READ(UNIT=66, '(a)', IOSTAT=i_iostat) read_line
  !  IF(i_iostat <0) EXIT ! fin du fichier
  !  nb_HKL_file = nb_HKL_file + 1
  !  IF(nb_HKL_file  > 1) exit
  !  SHELX_HKL_file = ADJUSTL(read_line)
  ! END do
  !CLOSE(UNIT=66)
  !call system('del hkl.lst')

  !IF(nb_HKL_file ==0) then
  ! call write_info('')
  ! WRITE(message_text, '(a)') '  No HKL file in the current folder !! '
  ! call write_info(TRIM(message_text))
  ! return
  !endif



  call write_info('')
  WRITE(message_text, '(2a)') '  . HKL file: ', TRIM(SHELX_HKL_file)
  call write_info(TRIM(message_text))

  OPEN(UNIT=66, FILE=TRIM(SHELX_HKL_file))
  OPEN(UNIT=67, FILE='import.cif', POSITION = 'append')
   WRITE(67, '(a)') ''
   WRITE(67, '(a)') '#----------------------------------------------------------------------------#'
   WRITE(67, '(a)') '#                   HKL DATA from SADABS                                     #'
   WRITE(67, '(a)') '#----------------------------------------------------------------------------#'
   WRITE(67, '(a)') ''
   WRITE(67, '(a)') 'loop_'
   WRITE(67, '(a)') '_refln_index_h'
   WRITE(67, '(a)') '_refln_index_k'
   WRITE(67, '(a)') '_refln_index_l'
   WRITE(67, '(a)') '_refln_F_squared_meas'
   WRITE(67, '(a)') '_refln_F_squared_sigma'
   WRITE(67, '(a)') '_refln_scale_group_code'

   do
    READ(66, '(a)', IOSTAT=i_iostat) read_line
    IF(i_iostat <0) exit
    READ(read_line, '(3I4,2F8.2,I4)')   h, k, l, F2, sig, code
    WRITE(67, '(I4,2I5,2F9.2,I5)') h, k, l, F2, sig, code
   END do
  CLOSE(UNIT=67)
  close(UNIT=66)


end subroutine read_SHELX_HKL


!-------------------------------------------------------------------------

subroutine read_RAW_file
! lecture fichier.RAW cree par SAINT
 USE IO_module,      ONLY  : write_info
 USE HKL_module
 USE cryscal_module, ONLY  : RAW_file_name, message_text

 implicit none
  !INTEGER              :: nb_HKL_file
  INTEGER              :: i_iostat, i1
  CHARACTER (LEN=256)  :: read_line
  CHARACTER (LEN=256)  :: SHELX_HKL_file
  INTEGER              :: h_, k_, l_, code_
  REAL                 :: F2_, sig_
  REAL, dimension(6)   :: cosdir
  INTEGER              :: div_10
  real                 :: max_F2


  i1 = INDEX(RAW_file_name, '.')
  SHELX_HKL_file = RAW_file_name(1:i1-1) //'_raw.HKL'


  OPEN(UNIT=66, FILE=TRIM(RAW_file_name))
  OPEN(UNIT=67, FILE=trim(SHELX_HKL_file))
   n_ref = 0
   do
    READ(66, '(a)', IOSTAT=i_iostat) read_line
    IF(i_iostat <0) exit
    
    READ(read_line,      '(3I4)')       h_, k_, l_
    READ(read_line(13:), *)             F2_, sig_
    READ(read_line(29:), '(I4,6F8.5)')  code_, cosdir(1:6)

    write(67, '(3I4, 2F8.2,I4, 6F8.5)') h_,k_,l_, F2_, sig_, code_, cosdir(1:6)
    
    n_ref = n_ref+1
    h(n_ref)               = h_
    k(n_ref)               = k_
    l(n_ref)               = l_
    F2(n_ref)              = F2_
    sig_F2(n_ref)          = sig_ 
    code(n_ref)            = code_
    cos_dir(n_ref,1:6)     = cosdir(1:6)

   END do
  CLOSE(UNIT=67)
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

  
  call write_info('')
  call write_info(' >> '//trim(SHELX_HKl_file) // ' file (SHELX format) has been created.')
  IF(div_10 /=0) then
   write (message_text, '(a,I2)') ' >> F2 and sigmas have been divided by 10**',div_10
   call write_info(TRIM(message_text))
  endif
  call write_info('')


 return
end subroutine read_RAW_file