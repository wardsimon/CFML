!     Last change:  TR   11 Dec 2007   11:13 am
! calculette cristallographique:

program crystallographic_calculations
 USE IO_module
 USE macros_module
 use cryscal_module
 USE hkl_module,                       ONLY: HKL_list

 implicit none
  CHARACTER (len=256)                      :: cmd_line
  CHARACTER (LEN=256)                      :: input_file, input_arg, arg_string, input_string
  CHARACTER (LEN=4)                        :: choice
  LOGICAL                                  :: file_exist
  LOGICAL                                  :: arg_keyword
  INTEGER                                  :: len_cmd_line, nb_arg
  INTEGER                                  :: i, i1, i2, nb_col
  character (len=256), dimension(10)       :: cmd_arg
  INTEGER                                  :: i_arg_DEBUG, i_arg_NOOUT, i_arg_CIFDEP, i_arg_ACTA


 ! >>>>  initialisation
  call Def_transformation_matrix
  call cryscal_init()
  tmp_logical = .false.
 

 ! >>>>  analyse de la ligne de commande
  input_file = ''
  input_arg  = ''
  cmd_line   = ''
  cmd_arg    = ''
  arg_keyword = .false.
  len_cmd_line = 0

  ! ligne de commande ?
  !CALL GETCL(cmd_line)  ! not valid in G95
  !len_cmd_line = len_trim(cmd_line)
  ! oct. 08
  !call nombre_de_colonnes(cmd_line, nb_arg)
  
  ! g95 :
  nb_arg = IARGC()   
  !

  if(nb_arg > 9) then
   write(*,*) ' '
   write(*,*) ' Number of arguments in the command line too large (max. = 10) !'
   write(*,*) ' Programm will be stopped.'
   write(*,*) ' '
   pause
   stop
  endif
  if(nb_arg /=0) then
   !read(cmd_line, *) (cmd_arg(i), i=1,nb_arg)     
   do i=1, nb_arg
    call GetArg(i, cmd_arg(i))   ! g95
	len_cmd_line = len_cmd_line + len_trim(cmd_arg(i)) + 1
   end do
   i_arg_DEBUG  = 0
   i_arg_NOOUT  = 0
   i_arg_CIFDEP = 0
   i_arg_ACTA   = 0
   do i = 1, nb_arg
    if(u_case(cmd_arg(i)(1:5)) == 'DEBUG' .or. u_case(cmd_arg(i)(1:3)) == 'LOG')  then
     DEBUG_file%write = .true.
     i_arg_DEBUG      = i     
    end if
    if(u_case(cmd_arg(i)(1:6)) == 'NO_OUT') then
     ON_SCREEN      = .false.
     i_arg_NOOUT    = i
    endif 
    if(u_case(cmd_arg(i)(1:6)) == 'CIFDEP') then
     CIFDEP         = .true.
     i_arg_CIFDEP   = i
    endif
    if(u_case(cmd_arg(i)(1:4)) == 'ACTA')   then
     ACTA           = .true.
     CIFDEP         = .true.
     i_arg_ACTA     = i
    endif
   end do 
   
   if(DEBUG_file%write) then
    open(unit = DEBUG_file%unit, file = 'cryscal_debug.txt')
    write(DEBUG_file%unit, '(2a)') ' COMMAND LINE : ', trim(cmd_line)
    if(nb_arg /=0) then
     do i=1, nb_arg
      write(DEBUG_file%unit, '(2x,a,i1,2a)') ' . arg_', i, ' : ', trim(cmd_arg(i))
     end do
	end if 
    write(DEBUG_file%unit,*) '' 
   end if
    
    ! creation de la ligne de commande exempte des arguments DEBUG, NOOUT, CIFDEP, ACTA
   if(i_arg_DEBUG /=0 .or. i_arg_NOOUT /=0 .or. i_arg_CIFDEP /=0 .or. i_arg_ACTA /=0) then
    cmd_line = ''
    do i=1, nb_arg 
     if (i/=i_arg_DEBUG .and. i/=i_arg_NOOUT .and. i/=i_arg_CIFDEP .and. i/=i_arg_ACTA) then
      cmd_line = trim(cmd_line)// ' ' //trim(cmd_arg(i))
     endif
    end do
    len_cmd_line = len_trim(cmd_line)
    call nombre_de_colonnes(cmd_line, nb_arg)
   end if
  end if 

 ! >>>>  lecture du fichier de configuration CRYSCAL.ini 
  call read_cryscal_ini()
  
 ! >>>> interpretation de l'argument #1 de la ligne de commande
  if(len_cmd_line /=0) then  
   !input_arg = u_case(cmd_line)
   input_arg = u_case(cmd_arg(1))
   i1 = INDEX(input_arg, '.', back=.true.)
   
   IF(input_arg(1:4) == 'HELP' .or. input_arg(1:3) == 'MAN') then
    keyword_HELP = .true.
    arg_keyword  = .true.
    read(cmd_line(5:), '(a)') arg_string
    call nombre_de_colonnes(arg_string, nb_col)
    IF(nb_col /=0) then
     nb_help = nb_col
     read(cmd_line(5:), *) HELP_arg(1:nb_help)
     do i=1, nb_help
      HELP_arg(i) = u_case(HELP_arg(i))
     end do
    END if
    !OPEN(UNIT = 2, FILE='CRYSCAL_manual.txt', ACTION='write', STATUS='replace')
	OPEN(UNIT = 2, FILE='CRYSCAL_manual.txt')
            

   ELSEIF(input_arg(1:4) == 'HTML') then
    keyword_create_CRYSCAL_HTML = .true.
    arg_keyword                 = .true.
    if(input_arg(1:11) == 'HTML_BROWSE' .or. input_arg(1:10) == 'HTM_BROWSE' ) browse_cryscal_HTML = .true.
    
   ELSEIF(input_arg(1:4) == 'NEWS') then
    keyword_create_CRYSCAL_news = .true.
    arg_keyword                 = .true.
    !OPEN(UNIT = NEWS_unit, FILE = 'cryscal_news.txt', ACTION='write', STATUS = 'replace')
    !OPEN(UNIT = 2, FILE = 'CRYSCAL_news.txt', ACTION='write', STATUS = 'replace')
	OPEN(UNIT = 2, FILE = 'CRYSCAL_news.txt')
    
   ELSEIF(input_arg(1:3) == 'KEY' .or. input_arg(1:4) == 'KEYS') then
    keyword_KEY = .true.
    arg_keyword = .true.
    !OPEN(UNIT=KEYS_unit, FILE='CRYSCAL_keys.txt',   ACTION='write', STATUS='replace')
    !OPEN(UNIT = 2, FILE = 'CRYSCAL_keys.txt',   ACTION='write', STATUS='replace')
	OPEN(UNIT = 2, FILE = 'CRYSCAL_keys.txt')

    !call write_KEYWORD()
   
   ELSEIF(input_arg(1:3) == 'CLA' ) then
    keyword_CLA = .true.
    arg_keyword = .true.
    !open(unit = 2, FILE = 'CRYSCAL_cla.txt', ACTION = 'write', STATUS = 'replace')
	open(unit = 2, FILE = 'CRYSCAL_cla.txt')
     
   ELSEIF(input_arg(1:3) == 'P4P' .and. len_trim(input_arg)==3)  then
    keyword_P4P = .true.
    arg_keyword = .true.
    P4P_file_name = ''
    
   ELSEIF(input_arg(i1+1:i1+3) == 'P4P') then
    keyword_P4P = .true.
    arg_keyword = .true.
    i1 = INDEX(input_arg, '"')
    i2 = INDEX(input_arg, '"', back=.true.)
    IF(i1 /=0 .AND. i2 /=0) then
     P4P_file_name = input_arg(i1+1:i2-1)
    else
     P4P_file_name = input_arg
    endif
   
   elseif(input_arg(i1+1:i1+3) == "RAW") then
    keyword_RAW = .true.
    arg_keyword = .true.
    i1 = INDEX(input_arg, '"')
    i2 = INDEX(input_arg, '"', back=.true.)
    IF(i1 /=0 .AND. i2 /=0) then
     RAW_file_name = input_arg(i1+1:i2-1)
    else
     RAW_file_name = input_arg
    endif

    
   ELSEIF(input_arg(1:6) == 'REPORT' .or. input_arg(1:11) == 'REPORT_LONG') then
    keyword_create_report = .true.
    arg_keyword = .true.
    long_report = .false.
    if(input_arg(1:11) == 'REPORT_LONG') long_report = .true.
    !if(nb_arg > 1  .and. u_case(cmd_arg(2)(1:5)) /= 'DEBUG'   &
	!               .and. u_case(cmd_arg(2)(1:3)) /= 'LOG'     &
    !               .and. u_case(cmd_arg(2)(1:6)) /= 'NO_OUT'  &
    !               .and. u_case(cmd_arg(2)(1:6)) /= 'CIFDEP'  &
    !               .and. u_case(cmd_arg(2)(1:4)) == 'ACTA') then
    if(nb_arg > 1) then
     read(cmd_arg(2), *) archive_CIF
    else
     archive_CIF = 'archive.cif'
    endif    
        
    
   ELSEIF(input_arg(1:11) == 'ARCHIVE.CIF') then
    keyword_modif_archive = .true.
    input_file = input_arg(1:11)
    arg_keyword = .true.
   
   elseif(input_arg(1:12) == 'SOLVE_TO_INS' .or. input_arg(1:10) == 'CREATE_INS') then
    keyword_SOLVE_to_INS = .true.
    arg_keyword =  .true.
   endif

   if(.not. keyword_HELP .and. .not. keyword_KEY .and. .not. keyword_create_CRYSCAL_HTML &
                                                 .and. .not. keyword_create_CRYSCAL_NEWS) then
    !open(UNIT=2, FILE='cryscal.log',        ACTION='write', status='replace')   
	open(UNIT=2, FILE='cryscal.log')   
   endif 
   
 
   IF(.NOT. arg_keyword) then
    input_file = input_arg
    i = index(input_file, '.', back=.true.)
    if(i==0) then
     input_file = trim(input_file) //'.cfl'
     input_CFL = .true.
    else
     i1 = INDEX(input_file, '"')
     i2 = INDEX(input_file, '"', back=.TRUE.)
     IF(i1 /=0 .AND. i2/=0) input_file = input_file(i1+1:i2-1)
     IF(input_file(i+1:) =='CFL') input_CFL = .true.
    endif
    call test_file_exist(input_file, file_exist)
    if(.not. file_exist) input_file = ''
   ENDIF
  
  else
   !open(UNIT=2, FILE='CRYSCAL.log',        ACTION='write', status='replace')    
   open(UNIT=2, FILE='CRYSCAL.log')    
  endif   ! fin de la condition if (len_cmd_line /=0


 if(nb_arg /=0) then
  !open(UNIT=2, FILE='CRYSCAL.log',        ACTION='write', status='replace')   
  open(UNIT=2, FILE='CRYSCAL.log')   
  call write_info('')
  call write_info('  . CRYSCAL arguments : ' //trim(cmd_line))
  call write_info('')
  !close(unit = 2)
 endif

  
 !call write_CRYSCAL_title()   ! help.F90

 if (keyword_HELP) then
  close(unit=2)
  !OPEN(UNIT=2, FILE='CRYSCAL_manual.txt', ACTION='write', STATUS='replace')
  OPEN(UNIT=2, FILE='CRYSCAL_manual.txt')
  call write_CRYSCAL_title
  call Write_HELP()
 end if
 
 if (keyword_KEY) then
  close(unit=2)
  !OPEN(UNIT=2, FILE='CRYSCAL_keys.txt',   ACTION='write', STATUS='replace')
  OPEN(UNIT=2, FILE='CRYSCAL_keys.txt')
  call write_CRYSCAL_title
  call Write_KEYWORD() 
 end if
 
 IF (keyword_create_CRYSCAL_HTML) then
  call write_CRYSCAL_title
  call Write_cryscal_HTML 
 endif
 
 if (keyword_CLA) then
  close(unit=2)
  !OPEN(UNIT=2, FILE='CRYSCAL_cla.txt', ACTION='write', STATUS='replace')
  OPEN(UNIT=2, FILE='CRYSCAL_cla.txt')
   call write_CRYSCAL_title
   call write_info('')
   call write_info(' >>> List of CRYSCAL commands line arguments :')
   call write_info('     ----------------------------------------')
   call write_info('')
   call write_cryscal_CLA
   call write_info('')
   call write_info('CRYSCAL news are stored in the CRYSCAL_cla.txt file.')
   call write_info('')
   stop
 endif
 
 if (keyword_create_CRYSCAL_NEWS) then
  close(unit=2)
  !OPEN(UNIT=2, FILE='cryscal_news.txt',   ACTION='write', STATUS='replace')
  OPEN(UNIT=2, FILE='cryscal_news.txt')
  call write_CRYSCAL_title

  call write_cryscal_NEWS('text')
  call write_info('')
  call write_info('CRYSCAL news are stored in the CRYSCAL_news.txt file.')
  call write_info('')
  stop
 endif
 
 
  call write_CRYSCAL_title()   ! help.F90

 !if (keyword_create_CEL)          call create_CEL_from_CIF

 IF (keyword_P4P) THEN
  call read_P4P(P4P_file_name)
  IF(molecule%formula(1:1) /= '?' ) then
   call nombre_de_colonnes(molecule%formula, nb_col)
   call decode_CHEM_string(molecule%formula, nb_col)
   call atomic_identification()
   call molecular_weight
  ENDIF

  call read_SADABS_file()
  call create_CIF_P4P_import('IMPORT')
  call read_SHELX_HKL
  call write_info('')
  call write_info('    >>> import.cif file has been created.')
  call write_info('')
  close(unit=2)
  call system('del CRYSCAL.log')
  stop
 
 elseif(keyword_RAW) then
  call read_RAW_file
  stop
  
 ELSEIF(keyword_create_REPORT) then
  call create_structural_report()

  stop
  
 !elseif(keyword_create_CEL) then
 ! OPEN(UNIT=CIF_read_unit, FILE=TRIM(CIF_file_name), ACTION="read")
 !  call check_CIF_input_file(TRIM(CIF_file_name))
 !  call read_CIF_input_file(TRIM(CIF_file_name), '?')
 !  call read_CIF_input_file_TR(CIF_read_unit)
 ! CLOSE(UNIT=CIF_read_unit)
 ! keyword_read_CIF = .true.
 ! call create_CEL_from_CIF
 ! stop
   
 elseif(keyword_SOLVE_to_INS) then
  call create_INS_from_SOLVE()
  stop 
 END IF


 if(len_trim(input_arg) /=0) then
  IF(input_CFL) then
   call read_keywords_from_file(TRIM(input_file))
  else
   call read_keywords_from_file(TRIM(input_file))                      ! cryscal.F90
   call run_keywords()                                                 ! cryscal.F90 
   if(keyword_modif_ARCHIVE) call read_and_modif_archive(input_unit)   ! read_cif_file.F90
   
  endif

 else
  do
   call WRITE_info("  ")
   call WRITE_info("  ****           CRYSCAL main menu          ****")
   call write_info("  ")
   call WRITE_info("     1. Interactive mode")
   call WRITE_info("     2. Read keywords list from input file")
   call WRITE_info("     3. User's guide")
   call WRITE_info("     4. List of keywords")
   call write_info("     5. List of command line arguments")
   call WRITE_info("     6. What's new in Cryscal ?")
   call WRITE_info("     0. Stop")
   call WRITE_info("  ")
   call WRITE_info("    Enter your choice [def=1] : ")

   !call read_input_line(input_string)
   !READ(input_string, '(a)') choice
   
   read(*,'(a)') choice
   
   choice = ADJUSTL(choice)
   if(len_trim(choice) == 0) choice = '1'
   
   select case (choice)
       case ('0', 'x', 'X', 'xx', 'XX', 'xX', 'Xx')
        stop

       case ('1')
        mode_interactif = .true.
        call interactive_mode('keyboard')
        !exit

       case ('2')
        mode_interactif = .false.
        call read_keywords_from_file('')
        call run_keywords()
        !exit

       case ('3')
        CLOSE(UNIT=2)
        !OPEN(UNIT=2, FILE='cryscal_manual.txt', ACTION='write', STATUS='replace')
		OPEN(UNIT=2, FILE='cryscal_manual.txt')
         call write_CRYSCAL_title
         IF(nb_help == nb_help_max) call write_header 
         call HELP_on_line
         call write_info('')
         call write_info(' >>> List of CRYSCAL commands line arguments :')
         call write_info('     ----------------------------------------')
         call write_info('')
         call write_cryscal_CLA
        CLOSE(UNIT=2)
        open(UNIT=2, FILE='CRYSCAL.log',        ACTION='write', position='append')
        call WRITE_info(' ')
        call WRITE_info('   CRYSCAL manual stored in CRYSCAL_manual.txt file. ')
        call WRITE_info(' ')
        !open(UNIT=2, FILE='CRYSCAL.log',        ACTION='write', position='append')
         call write_cryscal_title()

       CASE('4')
        CLOSE(UNIT=2)
        !OPEN(UNIT=2, FILE='CRYSCAL_keys.txt', ACTION='write', status='replace')
		OPEN(UNIT=2, FILE='CRYSCAL_keys.txt')
         call write_CRYSCAL_title         
         call keys_on_line()
        CLOSE(UNIT=2)
        open(UNIT=2, FILE='CRYSCAL.log',        ACTION='write', position='append')
        call WRITE_info(' ')
        call WRITE_info(' ')
        call WRITE_info('   CRYSCAL keys stored in CRYSCAL_keys.txt file. ')
        call WRITE_info(' ')
        !open(UNIT=2, FILE='CRYSCAL.log',        ACTION='write', position='append')
        call write_cryscal_title()

       CASE('5')
        CLOSE(UNIT=2)
        OPEN(UNIT=2, FILE='CRYSCAL_cla.txt', ACTION='write', status='replace')
         call write_CRYSCAL_title         
         call write_info('')
         call write_info(' >>> List of CRYSCAL commands line arguments :')
         call write_info('     ----------------------------------------')
         call write_info('')
         call write_cryscal_CLA
        CLOSE(UNIT=2)
        open(UNIT=2, FILE='CRYSCAL.log',        ACTION='write', position='append')
        call WRITE_info(' ')
        call WRITE_info('   CRYSCAL commands line arguments stored in CRYSCAL_cla.txt file. ')
        call WRITE_info(' ')
        call write_cryscal_title()
        
       CASE('6')
        CLOSE(UNIT=2)
        OPEN(UNIT=2, FILE='CRYSCAL_news.txt', ACTION='write', status='replace')
         call write_CRYSCAL_title         
         call write_info('')
         call write_cryscal_NEWS('screen') 
        close(unit = 2)
        open(UNIT=2, FILE='CRYSCAL.log',        ACTION='write', position='append')
        call WRITE_info(' ')
        call WRITE_info('   CRYSCAL news in CRYSCAL_news.txt file. ')
        call WRITE_info(' ')
        call write_cryscal_title()
        
       case default
   end select
  END do
 endif

 call end_of_program()

end program crystallographic_calculations


!--------------------------------------------------------------------
subroutine end_of_program
 USE IO_module
 USE cryscal_module, ONLY : keyword_create_CIF, keyword_modif_ARCHIVE


 call WRITE_info(' ')
 call WRITE_info(' ')
 call WRITE_info('   All CRYSCAL results are stored in CRYSCAL.log file. ')
 call WRITE_info(' ')

 IF(keyword_create_CIF) then
  call WRITE_info(' ')
  call WRITE_info(' ')
  if (keyword_modif_ARCHIVE) then
   call WRITE_info('   Results in CIF format are stored in ARCHIVE_CRYSCAL.CIF file. ')
  else 
   call WRITE_info('   Results in CIF format are stored in CRYSCAL.CIF file. ')
  endif 
  call WRITE_info(' ')
 ENDIF

 stop

end subroutine end_of_program
