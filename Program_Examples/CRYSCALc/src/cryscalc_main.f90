!     Last change:  TR   11 Dec 2007   11:13 am
! calculette cristallographique:
! calcul de distances interatomiques et d'angles entre 2 vecteurs
! note: d:\mes_documents\notes\math\produit_scalaire.doc

program crystallographic_calculations
 USE IO_module
 USE macros_module
 use cryscalc_module
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
  INTEGER                                  :: i_arg_LOG, i_arg_NOOUT, i_arg_CIFDEP, i_arg_ACTA


 ! >>>>  initialisation
  call cryscalc_init()
  call Def_transformation_matrix
   

 ! >>>>  analyse de la ligne de commande
  input_file = ''
  input_arg  = ''
  arg_keyword = .false.

  ! ligne de commande ?
  CALL GETCL(cmd_line)
  len_cmd_line = len_trim(cmd_line)
  
  ! oct. 08
  call nombre_de_colonnes(cmd_line, nb_arg)
  if(nb_arg > 9) then
   write(*,*) ' '
   write(*,*) ' Number of arguments in the command line too large (max. = 10) !'
   write(*,*) ' Programm will be stopped.'
   write(*,*) ' '
   pause
   stop
  endif
  if(nb_arg /=0) then
   read(cmd_line, *) (cmd_arg(i), i=1,nb_arg)   
   i_arg_LOG    = 0
   i_arg_NOOUT  = 0
   i_arg_CIFDEP = 0
   i_arg_ACTA   = 0
   do i = 1, nb_arg
    if(u_case(cmd_arg(i)(1:3)) == 'LOG')  then
     LOG_file%write = .true.
     i_arg_LOG      = i     
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
   
   if(LOG_file%write) then
    open(unit = LOG_file%unit, file = 'cryscalc.log')
    write(unit=LOG_file%unit, '(2a)') ' COMMAND LINE : ', trim(cmd_line)
    do i=1, nb_arg
     write(unit=LOG_file%unit, '(2x,a,i1,2a)') ' . arg_', i, ' : ', trim(cmd_arg(i))
    end do
    write(unit=LOG_file%unit,*) '' 
   end if
    
    ! creation de la ligne de commande exempte des arguments LOG, NOOUT, CIFDEP, ACTA
   if(i_arg_LOG /=0 .or. i_arg_NOOUT /=0 .or. i_arg_CIFDEP /=0 .or. i_arg_ACTA /=0) then
    cmd_line = ''
    do i=1, nb_arg 
     if (i/=i_arg_LOG .and. i/=i_arg_NOOUT .and. i/=i_arg_CIFDEP .and. i/=i_arg_ACTA) then
      cmd_line = trim(cmd_line)// ' ' //trim(cmd_arg(i))
     endif
    end do
    len_cmd_line = len_trim(cmd_line)
    call nombre_de_colonnes(cmd_line, nb_arg)
   end if
  end if 

 ! >>>>  lecture du fichier de configuration CRYSCALC.ini 
  call read_cryscalc_ini()
  
 ! >>>> interpretation de l'argument #1 de la ligne de commande
  if(len_cmd_line /=0) then  
   !input_arg = u_case(cmd_line)
   input_arg = u_case(cmd_arg(1))
   i1 = INDEX(input_arg, '.', back=.true.)
   
   IF((len_trim(input_arg) == 4 .and. input_arg(1:4) == 'HELP') .or. &
      (len_trim(input_arg) == 3 .and. input_arg(1:3) == 'MAN')) then
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
    OPEN(UNIT = 2, FILE='CRYSCALC_manual.txt', ACTION='write', STATUS='replace')
            

   ELSEIF(len_trim(input_arg) == 4 .and. input_arg(1:4) == 'HTML') then
    keyword_create_CRYSCALC_HTML = .true.
    arg_keyword                 = .true.
    if(input_arg(1:11) == 'HTML_BROWSE' .or. input_arg(1:10) == 'HTM_BROWSE' ) browse_cryscalc_HTML = .true.
    
   ELSEIF(len_trim(input_arg) == 4 .and. input_arg(1:4) == 'NEWS') then
    keyword_create_CRYSCALc_news = .true.
    arg_keyword                  = .true.
    !OPEN(UNIT = NEWS_unit, FILE = 'cryscalc_news.txt', ACTION='write', STATUS = 'replace')
    OPEN(UNIT = 2, FILE = 'cryscalc_news.txt', ACTION='write', STATUS = 'replace')
    
   ELSEIF((len_trim(input_arg)==3 .and. input_arg(1:3) == 'KEY') .or. &
          (len_trim(input_arg)==4 .and. input_arg(1:4) == 'KEYS')) then
    keyword_KEY = .true.
    arg_keyword = .true.
    !OPEN(UNIT=KEYS_unit, FILE='CRYSCALC_keys.txt',   ACTION='write', STATUS='replace')
    OPEN(UNIT = 2, FILE = 'CRYSCALC_keys.txt',   ACTION='write', STATUS='replace')

    !call write_KEYWORD()
   
   ELSEIF(len_trim(input_arg) == 3 .and. input_arg(1:3) == 'CLA' ) then
    keyword_CLA = .true.
    arg_keyword = .true.
    open(unit = 2, FILE = 'CRYSCALC_cla.txt', ACTION = 'write', STATUS = 'replace')
     
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
    !if(nb_arg > 1  .and. u_case(cmd_arg(2)(1:3)) /= 'LOG'     &
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
   
   elseif(input_arg(1:12) == 'SOLVE_TO_INS' .or. input_arg(1:10) == 'CREATE_JOB_INS') then
    keyword_SOLVE_to_INS = .true.
    arg_keyword =  .true.
    create_INS%ANIS = .false.
    if(nb_arg > 1 .and. u_case(cmd_arg(2)(1:4)) == 'ANIS') create_INS%ANIS = .true.
	
   endif

   if(.not. keyword_HELP .and. .not. keyword_KEY .and. .not. keyword_create_CRYSCALC_HTML &
                                                 .and. .not. keyword_create_CRYSCALC_NEWS) then
    open(UNIT=2, FILE='cryscalc.out',        ACTION='write', status='replace')   
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
   open(UNIT=2, FILE='CRYSCALC.out',        ACTION='write', status='replace')    
  endif   ! fin de la condition if (len_cmd_line /=0


 if(nb_arg /=0) then
  open(UNIT=2, FILE='CRYSCALC.out',        ACTION='write', status='replace')   
  call write_info('')
  call write_info('  . CRYSCALC arguments : ' //trim(cmd_line))
  call write_info('')
  !close(unit = 2)
 endif

  
 !call write_CRYSCALC_title()   ! help.F90

 if (keyword_HELP) then
  close(unit=2)
  OPEN(UNIT=2, FILE='CRYSCALC_manual.txt', ACTION='write', STATUS='replace')
  call write_CRYSCALC_title
  call Write_HELP()
 end if
 
 if (keyword_KEY) then
  close(unit=2)
  OPEN(UNIT=2, FILE='CRYSCALC_keys.txt',   ACTION='write', STATUS='replace')
  call write_CRYSCALC_title
  call Write_KEYWORD() 
 end if
 
 IF (keyword_create_CRYSCALC_HTML) then
  call write_CRYSCALC_title
  call Write_cryscalc_HTML 
 endif
 
 if (keyword_CLA) then
  close(unit=2)
  OPEN(UNIT=2, FILE='CRYSCALC_cla.txt', ACTION='write', STATUS='replace')
   call write_CRYSCALC_title
   call write_info('')
   call write_info(' >>> List of CRYSCALC commands line arguments :')
   call write_info('     -----------------------------------------')
   call write_info('')
   call write_cryscalc_CLA
   call write_info('')
   call write_info('CRYSCALC commands line arguments are stored in the CRYSCALC_cla.txt file.')
   call write_info('')
   stop
 endif
 
 if (keyword_create_CRYSCALC_NEWS) then
  close(unit=2)
  OPEN(UNIT=2, FILE='cryscalc_news.txt',   ACTION='write', STATUS='replace')
  call write_CRYSCALC_title

  call write_cryscalc_NEWS('text')
  call write_info('')
  call write_info('CRYSCALC news are stored in the CRYSCALC_news.txt file.')
  call write_info('')
  stop
 endif
 
 
  call write_CRYSCALC_title()   ! help.F90

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
  call system('del CRYSCALC.out')
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
   call read_keywords_from_file(TRIM(input_file))                      ! cryscalc.F90
   call run_keywords()                                                 ! cryscalc.F90 
   if(keyword_modif_ARCHIVE) call read_and_modif_archive(input_unit)   ! read_cif_file.F90
   
  endif

 else
  do
   call WRITE_info("  ")
   call WRITE_info("  ****           CRYSCALC main menu          ****")
   call write_info("  ")
   call WRITE_info("     1. Interactive mode")
   call WRITE_info("     2. Read keywords list from input file")
   call WRITE_info("     3. User's guide")
   call WRITE_info("     4. List of keywords")
   call write_info("     5. List of command line arguments")
   call WRITE_info("     6. What's new in Cryscalc ?")
   call WRITE_info("     0. Stop")
   call WRITE_info("  ")
   call WRITE_info("    Enter your choice [def=1] : ")

   !call read_input_line(input_string)
   !READ(input_string, '(a)') choice
   
   read(*,'(a)') choice
   
   choice = ADJUSTL(choice)
   if(len_trim(choice) == 0) choice = '1'
   
   select case (choice)
       case ('0', 'x', 'X')
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
        OPEN(UNIT=2, FILE='cryscalc_manual.txt', ACTION='write', STATUS='replace')
         call write_CRYSCALC_title
         IF(nb_help == nb_help_max) call write_header 
         call HELP_on_line
         call write_info('')
         call write_info(' >>> List of CRYSCALC commands line arguments :')
         call write_info('     -----------------------------------------')
         call write_info('')
         call write_cryscalc_CLA
        CLOSE(UNIT=2)
        open(UNIT=2, FILE='CRYSCALC.out',        ACTION='write', position='append')
        call WRITE_info(' ')
        call WRITE_info('   CRYSCALC manual stored in CRYSCALC_manual.txt file. ')
        call WRITE_info(' ')
        !open(UNIT=2, FILE='CRYSCALC.out',        ACTION='write', position='append')
         call write_cryscalc_title()

       CASE('4')
        CLOSE(UNIT=2)
        OPEN(UNIT=2, FILE='CRYSCALC_keys.txt', ACTION='write', status='replace')
         call write_CRYSCALC_title         
         call keys_on_line()
        CLOSE(UNIT=2)
        open(UNIT=2, FILE='CRYSCALC.out',        ACTION='write', position='append')
        call WRITE_info(' ')
        call WRITE_info(' ')
        call WRITE_info('   CRYSCALC keys stored in CRYSCALC_keys.txt file. ')
        call WRITE_info(' ')
        !open(UNIT=2, FILE='CRYSCALC.out',        ACTION='write', position='append')
        call write_cryscalc_title()

       CASE('5')
        CLOSE(UNIT=2)
        OPEN(UNIT=2, FILE='CRYSCALC_cla.txt', ACTION='write', status='replace')
         call write_CRYSCALC_title         
         call write_info('')
         call write_info(' >>> List of CRYSCALC commands line arguments :')
         call write_info('     -----------------------------------------')
         call write_info('')
         call write_cryscalc_CLA
        CLOSE(UNIT=2)
        open(UNIT=2, FILE='CRYSCALC.out',        ACTION='write', position='append')
        call WRITE_info(' ')
        call WRITE_info('   CRYSCALC commands line arguments stored in CRYSCALC_cla.txt file. ')
        call WRITE_info(' ')
        call write_cryscalc_title()
        
       CASE('6')
        CLOSE(UNIT=2)
        OPEN(UNIT=2, FILE='cryscalc_news.txt', ACTION='write', status='replace')
         call write_CRYSCALC_title         
         call write_info('')
         call write_cryscalc_NEWS('screen') 
        close(unit = 2)
        open(UNIT=2, FILE='CRYSCALC.out',        ACTION='write', position='append')
        call WRITE_info(' ')
        call WRITE_info('   CRYSCALC news in cryscalc_news.txt file. ')
        call WRITE_info(' ')
        call write_cryscalc_title()
        
       case default
   end select
  END do
 endif

 call end_of_program()

end program crystallographic_calculations


!--------------------------------------------------------------------
subroutine end_of_program
 USE IO_module
 USE cryscalc_module, ONLY : keyword_create_CIF, keyword_modif_ARCHIVE


 call WRITE_info(' ')
 call WRITE_info(' ')
 call WRITE_info('   All CRYSCALC results are stored in CRYSCALC.out file. ')
 call WRITE_info(' ')

 IF(keyword_create_CIF) then
  call WRITE_info(' ')
  call WRITE_info(' ')
  if (keyword_modif_ARCHIVE) then
   call WRITE_info('   Results in CIF format are stored in ARCHIVE_CRYSCALC.CIF file. ')
  else 
   call WRITE_info('   Results in CIF format are stored in CRYSCALC.CIF file. ')
  endif 
  call WRITE_info(' ')
 ENDIF

 stop

end subroutine end_of_program
