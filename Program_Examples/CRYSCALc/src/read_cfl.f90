!     Last change:  TR   15 Jun 2007   11:00 am



subroutine read_CFL_input_file(inp_unit)
 USE macros_module
 use cryscalc_module

 implicit none
  integer, intent(in)                      :: inp_unit
  CHARACTER (LEN=256)                      :: read_line
  integer                                  :: i_error, i1

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_CFL_INPUT_FILE")

! ------------------------------
 !open(unit = inp_unit, file = input_file)
 do        ! lecture du fichier d'entree
  !READ(unit=input_unit, fmt='(a)', IOSTAT=i_error) read_line
  !READ(unit=CFL_read_unit, fmt='(a)', IOSTAT=i_error) read_line
  READ(unit=inp_unit, fmt='(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) EXIT   ! fin du fichier
  read_line = ADJUSTL(read_line)
  read_line = u_case(TRIM(read_line))

  IF (LEN_TRIM(read_line) == 0) cycle
  IF (read_line(1:1) == '! ' .or. read_line(1:1) == '#') cycle
  
  i1 = index(read_line, '!')
  if(i1 /=0) read_line = trim(read_line(1:i1-1))
  write(*,*) ' >> ', trim(read_line)
  pause

  call identification_CFL_keywords(read_line)
  call identification_keywords(read_line)
  
 ! call run_keyword_interactive(current_keyword)

 END do
 close(unit=inp_unit)

! !IF(.NOT. keyword_BEAM .AND. keyword_WAVE .and. keyword_SFAC_UNIT ) then
! IF(.NOT. keyword_BEAM .AND. ( keyword_SFAC_UNIT .or. keyword_CONT .or. keyword_CHEM)) then
!  call incident_beam()
! endif
!
! IF(keyword_CHEM) then
!  IF(.NOT. keyword_ZUNIT) then
!   call write_info('')
!   call write_info('  >> ZUNIT keyword is missing ! ')
!   call write_info('')
!   stop
!  !else
!  !do i = 1, nb_atoms_type
!  !  SFAC_number(i) = SFAC_number(i) * Z_unit
!  ! END do
!  endif
! endif

  return
 end subroutine read_CFL_input_file
!---------------------------------------------------------------------------------

subroutine incident_beam()
 USE cryscalc_module, ONLY : X_rays, neutrons, wavelength, LOCK_wave_value, debug_proc
 USE wavelength_module
 implicit none
  integer             :: i

  if(debug_proc%level_2)  call write_debug_proc_level(2, "INCIDENT_BEAM")
  
  X_rays   = .true.
  neutrons = .false.
  
  X_target(1:tabulated_target_nb)%logic = .false.
  
  do i=1, tabulated_target_nb
   if(ABS(wavelength - X_target(i)%wave(1)) < LOCK_wave_value) then
    wavelength = X_target(i)%wave(1)
    X_target(i)%logic = .true.
    anti_cathode    = .true.
    exit
   end if
  end do
  
  !if(.not. anti_cathode) then
  ! X_target(2)%logic = .true.  ! Mo
  ! anti_cathode = .true.
  !endif

  return
end subroutine incident_beam
!-------------------------------------------------------------------------
 subroutine identification_CFL_keywords(read_line)
  USE cryscalc_module
  USE wavelength_module
  USE macros_module
  USE text_module ,              ONLY : CIF_lines_nb, CIF_title_line

  USE IO_module
  implicit none
  CHARACTER (LEN=*), INTENT(IN)            :: read_line
  CHARACTER (LEN=32)                       :: input_keyword, input_arg
  INTEGER                                  :: long_kw
  REAL,               DIMENSION(10)        :: var
  CHARACTER (LEN=256)                      :: new_line, arg_line
  CHARACTER (LEN=64), DIMENSION(20)        :: arg_string
  INTEGER                                  :: i
  integer                                  :: i1, i2, i_pos, i_error, nb_arg
  !LOGICAL                                  :: lecture_ok

  if(debug_proc%level_2)  call write_debug_proc_level(2, "IDENTIFICATION_CFL_KEYWORDS")
  
  READ(read_line, *) input_keyword
  input_keyword = ADJUSTL(input_keyword)
  long_kw = LEN_TRIM(input_keyword)
  if (input_keyword(long_kw:long_kw) == '=') then
   input_keyword = input_keyword(1:long_kw)
   long_kw = LEN_TRIM(input_keyword)
  end if

  READ(read_line(long_kw+1:),'(a)') arg_line
  arg_line = ADJUSTL(arg_line)
  IF(arg_line(1:1) == '=') then
   arg_line = arg_line(2:)
   arg_line = ADJUSTL(arg_line)
  endif
  call nombre_de_colonnes(arg_line, nb_arg)
  IF(nb_arg/=0)  then
   READ(arg_line, *, IOSTAT=i_error) arg_string(1:nb_arg)
   IF(i_error /=0) then
    call write_info('')
    WRITE(message_text, '(3a)') ' ... Error reading ', TRIM(input_keyword), ' arguments ...'
    call write_info(TRIM(message_text))
    return
   endif
   do i=1, nb_arg
    arg_string(i) = u_case(arg_string(i))
   end do
  endif

  IF(mode_interactif)   WRITE(CFL_unit, '(a)') TRIM(read_line)

 
  select case (TRIM(input_keyword))

   case ('ACTA', 'CIF', 'CREATE_CIF')
    keyword_create_CIF = .true.
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
    do i=1, CIF_lines_nb
     WRITE(CIF_unit, '(a)') trim(CIF_title_line(i))
    end do
    WRITE(CIF_unit, '(a)') ''
    WRITE(CIF_unit, '(a)') 'data_cryscalc'
    WRITE(CIF_unit, '(a)') ''
    
   case ('CREATE_ACE')
    keyword_create_ACE = .true.

   case ('CREATE_CEL')
    keyword_create_CEL = .true.
    
   case ('CREATE_CFL')
    keyword_create_CFL = .true.
    
   case ('CREATE_FST')
    keyword_create_FST = .true.
	
   !case ('CREATE_CIF_PYMOL')
   ! create_CIF_PYMOL = .true.   
     
   case ('CREATE_INS')
    keyword_create_INS = .true.  
      
 
   case ('TITL', 'TITLE')

    READ(arg_line, '(a)') main_title
    call write_info(' ' )
    call write_info('  . TITL: '//trim(main_title))

   CASE ('CELL', 'CELL_PARAMETERS')
    input_CELL_P4P  = .false.
    input_CELL_M50  = .false.
    input_CELL_INS  = .false.
	input_CELL_CIF  = .false.
	input_CELL_RED  = .false.
	
    input_CELL = .false.
    IF(nb_arg /=0) then
     i1 = INDEX(arg_string(1), '.')
     IF(i1/=0) then
      IF(arg_string(1)(i1:) == '.P4P') then
       input_CELL_P4P  = .true.
       input_CELL      = .true.
       WRITE(P4P_file_name, '(a)') arg_string(1)
      ELSEIF(arg_string(1)(i1:) == '.M50') then
       input_CELL_M50  = .true.
       input_CELL      = .true.
       WRITE(M50_file_name, '(a)') arg_string(1)
	  ELSEIF(arg_string(1)(i1:) == '.CIF' .OR.  arg_string(1)(i1:) == '.CIF_OD') then
       input_CELL_CIF  = .true.
       input_CELL      = .true.	
       WRITE(CIF_file_name, '(a)') arg_string(1)	   
      ELSEIF(arg_string(1)(i1:) == '.INS' .OR.  arg_string(1)(i1:) == '.RES') then
       input_CELL_INS = .true.
       input_CELL     = .true.
       WRITE(INS_file_name, '(a)') arg_string(1)
      ELSEIF(arg_string(1)(i1:) == '.X') then
       input_CELL_X = .true.
       input_CELL   = .true.
       WRITE(X_file_name, '(a)') arg_string(1)
      ELSEIF(arg_string(1)(i1:) == '.RMAT') then
       input_CELL_RMAT = .true.
       input_CELL      = .true.
       WRITE(RMAT_file_name, '(a)') arg_string(1)
	  ELSEIF(arg_string(1)(i1:) == '.RED') then
       input_CELL_RED  = .true.
       input_CELL	   = .true.
	   WRITE(RED_file_name, '(a)') arg_string(1)
      endif
     endif
    endif

    IF(.not. input_CELL) then
     IF(nb_arg == 1) then
      arg_string(2:3) = arg_string(1)
      arg_string(4:6) = '90.0'
     elseif(nb_arg == 2) then      ! systeme tetragonal
      arg_string(3) = arg_string(2)  ! parametre C
      arg_string(2) = arg_string(1)      
      arg_string(4:6) = '90.0' 
     elseIF(nb_arg == 3) then
      arg_string(4:6) = '90.0'
     ELSEIF(nb_arg == 4) then
      arg_string(5) = arg_string(4)
      arg_string(4) = '90.0'
      arg_string(6) = '90.0'
     elseif (nb_arg < 6) then
      call check_arg('6' , 'CELL')
      return
     endif
    endif

    var(1:10) = 0.
    IF(.NOT. input_CELL) then
     do i=1, 6
      i1 = INDEX(arg_string(i), '(')
      i2 = INDEX(arg_string(i), ')')
      IF(i1/=0 .AND. i2/=0) then
       READ(arg_string(i)(1:i1-1),*, IOSTAT=i_error) var(i)
       IF(i_error /=0) then
        call error_message('CELL')
        return
       endif
       !READ(arg_string(i)(i1+1:i2-1),*, IOSTAT=i_error) cell_param_esd(i)
       !IF(i_error /=0) call error_message('CELL')
       known_cell_esd = .true.
      else
       READ(arg_string(i), *, IOSTAT=i_error) var(i)
       IF(i_error /=0) then
        call error_message('CELL')
        return
       endif
      endif
     end do

     !READ(arg_string(1:6), *, IOSTAT=i_error) (var(i), i=1,6)
     !IF(i_error /=0) call error_message('CELL')

     !if(var(1) < 2.) then ! for compatibility with INS SHELX input file (CELL keyword)
     IF(nb_arg == 7) then
      READ(arg_string(1:7), *, IOSTAT=i_error) (var(i), i=1,7)
      IF(i_error /=0) then
       call error_message('CELL')
       return
      endif
      wavelength = var(1)
      call write_info(' ')
      WRITE(message_text,'( a,F10.5)') '  > WAVELENGTH (A):', wavelength
       call write_info(TRIM(message_text))
      IF(keyword_create_CIF)  call write_CIF_file('WAVE')
      keyword_WAVE = .true.
      unit_cell%param(1:6) = var(2:7)
     else
      unit_cell%param(1:6) = var(1:6)
     endif

    else
     IF(input_CELL_P4P) then
      call read_P4P_file(TRIM(P4P_file_name))
      IF(.NOT. lecture_ok) then
       call write_info('')
       call write_info('   >>> Problem reading P4P file <<<')
       call write_info('')
       return
      endif
     ELSEIF(input_CELL_M50) then
      call read_M50_file(TRIM(M50_file_name), lecture_ok)
      IF(.NOT. lecture_ok) then
       call write_info('')
       call write_info('   >>> Problem reading M50 file <<<')
       call write_info('')
       return
      endif
	 ELSEIF(input_CELL_CIF) then
      call READ_CIF_input_file(TRIM(CIF_file_name), 'out')
      
	  
     ELSEIF(input_CELL_INS) then
      call read_INS_input_file(TRIM(INS_file_name), "CELL")
       !IF(.NOT. lecture_ok) then
       ! call write_info('')
       ! call write_info('   >>> Problem reading INS file <<<')
       ! call write_info('')
       ! return
       !endif
     ELSEIF(input_CELL_X) then
      call READ_X_input_file(TRIM(X_file_name), lecture_OK)
      IF(.NOT. lecture_OK) then
       call write_info('')
       call write_info('   >>> Problem reading DENZO .x file')
       call write_info('')
      endif
     ELSEIF(input_CELL_RMAT) then
      call READ_RMAT_input_file(TRIM(RMAT_file_name), lecture_OK)
      IF(.NOT. lecture_OK) then
       call write_info('')
       call write_info('   >>> Problem reading DIRAX .rmat file')
       call write_info('')
      endif
	 ELSEIF(input_CELL_RED) then
      call READ_RED_input_file(TRIM(RED_file_name), lecture_OK)
      IF(.NOT. lecture_OK) then
       call write_info('')
       call write_info('   >>> Problem reading DATARED .red file')
       call write_info('')
      endif
     endif
    endif

    call write_info('')
    WRITE(message_text,'( a,6F10.4)') '  > CELL PARAMETERS: ', (unit_cell%param(i),i=1,6)
    call write_info(TRIM(message_text))
    keyword_CELL = .true.

    !call volume_calculation('out')
    !ESD volume
    !call setnum_std(volum(1,k,n_pat),SQRT(volum(2,k,n_pat)),line)
    IF(keyword_FILE .AND. .NOT. known_theta) call calcul_theta()


   CASE ('SIZE', 'CRYSTAL_SIZE')
    IF(nb_arg < 3) then
     call check_arg('3', 'SIZE')
     return
    endif

    var(1:10) = 0.
    read(arg_string(1:3), *, IOSTAT=i_error) (var(i), i=1,3)
    IF(i_error /=0) then
     call error_message('SIZE')
     return
    endif

    crystal%size(1:3) = var(1:3)
    keyword_SIZE = .true.

   CASE ('WAVE', 'WAVELENGTH', 'WL')
    IF(nb_arg ==0) then
     call check_arg('0', 'WAVE')
     return
    endif
    var(1:10) = 0.
	if(beam_type(1:7) == 'x_rays') then
    call get_X_radiation(arg_string(1))
    !if(.not. keyword_WAVE) then
    if(.not. anti_cathode) then
     read(arg_string(1), *,  IOSTAT= i_error) var(1)
     IF(i_error /=0) then
      call error_message('WAVE')
      return
     endif     
     wavelength = var(1)
    endif
    endif    

    call write_wave_features()
    keyword_WAVE = .true.
    IF(keyword_create_CIF)  call write_CIF_file('WAVE')

    !anti_cathode = .false.
    !do i=1, tabulated_target_nb 
    ! IF(X_target(i)%logic) anti_cathode = .true.
    !ENDDO
    
    IF(.not. anti_cathode) call incident_beam       !

    IF(keyword_FILE .AND. .NOT. known_theta) call calcul_theta()


   CASE ('BEAM')
    IF(nb_arg ==0) then
     call check_arg('0', 'BEAM')
     return
    endif
    read(arg_string(1),*, iostat=i_error) input_arg
    IF(i_error /=0) then
     call error_message('BEAM')
     return
    endif
    input_arg = u_case(TRIM(input_arg))
    IF(input_arg(1:1) == 'X') then
     X_rays    = .true.
     neutrons  = .false.
     electrons = .false.
     beam_type = 'x_rays'
     !call get_X_radiation(beam_type)
	 call get_X_radiation(trim(input_arg))
     keyword_wave = .true.
    ELSEIF(input_arg(1:4) == 'NEUT') then
     neutrons  = .true.
     X_rays    = .false.
     electrons = .false.
     beam_type = 'neutrons'
     if(nb_arg==2) then
      read(arg_string(2), *, iostat=i_error) var(2)
      if(i_error==0) then
       read(arg_string(2),*) wavelength
       keyword_wave = .true.
      endif
     endif
    ELSEIF(input_arg(1:4) == 'ELEC') then
     electrons = .true.
     neutrons  = .false.
     X_rays    = .false.     
     beam_type = 'electrons'
     if(nb_arg==2) then
      read(arg_string(2), *, iostat=i_error) var(2)
      if(i_error==0) then
       read(arg_string(2),*) wavelength
       keyword_wave = .true.
      endif
     endif
    endif
    keyword_beam = .true.
	call write_info('    . beam_type: '// trim(beam_type))


   CASE ('SPGR', 'SG', 'SPACE_GROUP')
    IF(nb_arg ==0) then
     call check_arg('0', 'SPGR')
     return
    endif
    read(arg_line, '(a)') space_group_symbol
    call write_info('')
    WRITE(message_text,'(2a)')   '  > SPACE GROUP : ', TRIM(space_group_symbol)
    call write_info(TRIM(message_text))
    keyword_SPGR = .true.


   CASE ('Z', 'ZUNIT', 'Z_UNIT' )
    IF(nb_arg ==0) then
     call check_arg('0', 'ZUNIT')
     return
    endif

    var(1:10) = 0.
    read(arg_string(1), *,  IOSTAT=i_error) var(1)
    IF(i_error /=0) then
     call error_message('Zunit')
     return
    endif

    Z_unit = var(1)
    call write_info('')
    write(message_text, '(a,F5.1)')  '  > Z (num. of molecular unit): ',Z_unit
    call write_info(TRIM(message_text))
    call write_info('')
    keyword_ZUNIT = .true.

    IF(keyword_CHEM) sfac_number(1:nb_atoms_type) = sto(1:nb_atoms_type) * Z_unit

    IF(keyword_create_CIF)  call write_CIF_file('ZUNIT')
    molecule%Z_unit = int(Z_unit)

   CASE ('SFAC')
    if(nb_arg == 0) then
     call check_arg('0', 'SFAC')
     return
    endif
    nb_atoms_type = nb_arg
    READ(arg_string(1:nb_arg),*) SFAC_type(1:nb_arg)
    call write_info('')
    WRITE(message_text,'(2a)') '  > ',TRIM(read_line)
    call write_info(TRIM(message_text))
    do i = 1, nb_atoms_type
     SFAC_type(i) = u_case(SFAC_type(i))
    end do


    IF(mode_interactif) then
     call write_info('')
     call write_info('  > Enter UNIT values:  ')
     call read_input_line(input_line)
     READ(input_line,*) (SFAC_number(i),i=1,nb_atoms_type)

    else
     READ(1, '(a)') new_line
     if (new_line(1:4)=='UNIT') then
      READ(new_line(5:), *) (SFAC_number(i),i=1,nb_atoms_type)

      WRITE(message_text,'(2a)') '  > ',TRIM(read_line)
      call write_info(TRIM(message_text))
      call write_info('')

     else
      call write_info('  ! UNIT string is missing ! ')
      stop
     end if
    endif

    keyword_SFAC_UNIT = .true.
    keyword_CHEM      = .false.
    keyword_CONT      = .false.

    call get_content()


   CASE ('CONT')
    IF(nb_arg ==0) then
     call check_arg('0', 'CONT')
     return
    endif

    nb_atoms_type = INT(nb_arg / 2.)
    READ(arg_line, *)( SFAC_type(i), SFAC_number(i), i = 1, nb_atoms_type)
    call write_info('')
    WRITE(message_text,'(2a)') '  > ',TRIM(read_line)
    call write_info(TRIM(message_text))

    do i = 1, nb_atoms_type
     SFAC_type(i) = u_case(SFAC_type(i))
    end do
    keyword_CONT      = .true.
    keyword_CHEM      = .false.
    keyword_SFAC_UNIT = .false.

    molecule%content = TRIM(read_line)

   CASE ('CHEM', 'CHEM_FORM', 'CHEMICAL_FORMULA')
    IF(nb_arg ==0) then
     call check_arg('0', 'CHEM')
     return
    endif
    call decode_CHEM_string(arg_line, nb_arg)

   case ('MU', 'CALC_MU', 'MU_CALC', 'ABSORPTION', 'ABSORPTION_CALC', 'CALC_ABSORPTION')
    keyword_MU = .true.
    
   CASE ('ATOM', 'ATM')
    IF(nb_arg <5) then
     call check_arg('5', 'ATOM')
     return
    endif

    nb_atom = nb_atom + 1

    READ(read_line(long_kw+1:), '(a)', IOSTAT=i_error) input_line
    IF(i_error /=0) then
     call error_message('ATOM')
     return
    endif
    input_line = adjustl(input_line)
    i = INDEX(input_line, ' ')
    if (i/=0) input_line = input_line(i+1:)
    input_line = ADJUSTL(input_line)

    i = INDEX(input_line, ' ')
    if (i/=0) input_line = input_line(i+1:)
    input_line = ADJUSTL(input_line)

    READ(arg_string(1),*) atom_label(nb_atom)
    READ(arg_string(2),*) atom_typ(nb_atom)

	
    do i=3, nb_arg
     i1 = INDEX(arg_string(i), '(')
     i2 = INDEX(arg_string(i), ')')
     IF(i1 /=0 .AND. i2/=0) then
      READ(arg_string(i)(1:i1-1), *, IOSTAT=i_error) arg_string(i)
      IF(i_error /=0) then
       call error_message('ATOM')
       return
      endif
     endif
    end do

    ! recherche de coordonnes sous forme fractionnaire
    atom_Biso(nb_atom)     = 0.
    atom_occ(nb_atom)      = 1.
    atom_occ_perc(nb_atom) = 1.
    i = INDEX(input_line, '/')
    IF (i/=0) THEN
     call get_atom_coord(input_line)
    else
     call nombre_de_colonnes(input_line, nb_arg)
      READ(arg_string(3), *) atom_coord(1, nb_atom)  ! x
      READ(arg_string(4), *) atom_coord(2, nb_atom)  ! y
      READ(arg_string(5), *) atom_coord(3, nb_atom)  ! z
      IF(nb_arg > 3) READ(arg_string(6), *) atom_Biso(nb_atom)      ! Biso
      IF(nb_arg > 4) READ(arg_string(7), *) atom_occ_perc(nb_atom)       ! occ

     !if(nb_arg ==3)      then    !  lecture x,y,z
     ! READ(input_line, *, IOSTAT=i_error) (atom_coord(i, nb_atom) ,i=1,3)
     !elseif(nb_arg == 4) then    !  lecture x,y,z,Biso
     ! READ(input_line, *, IOSTAT=i_error) (atom_coord(i, nb_atom) ,i=1,3), atom_Biso(nb_atom)
     !elseif(nb_arg == 5) then    !  lecture x,y,z,Biso, occ
     ! READ(input_line, *, IOSTAT=i_error) (atom_coord(i, nb_atom) ,i=1,3), atom_Biso(nb_atom), atom_occ_perc(nb_atom)
     !endif
    END IF


    call write_info(' ')
    if(i /=0)  then
      WRITE(message_text,'(a,2a6,5(1x,F9.6))') '  > ATOM: ', trim(atom_label(nb_atom)),trim(atom_typ(nb_atom)),  &
	                                           (atom_coord(i,nb_atom),i=1,3), atom_Biso(nb_atom), atom_occ_perc(nb_atom)
    else
     IF(nb_arg == 3) then
      WRITE(message_text,'(a,2a6,3(1x,F9.6))') '  > ATOM: ', trim(atom_label(nb_atom)),trim(atom_typ(nb_atom)),  &
	                                           (atom_coord(i,nb_atom),i=1,3)
     ELSEIF(nb_arg == 4) then
      WRITE(message_text,'(a,2a6,4(1x,F9.6))') '  > ATOM: ', trim(atom_label(nb_atom)),trim(atom_typ(nb_atom)),  &
	                                           (atom_coord(i,nb_atom),i=1,3), atom_Biso(nb_atom)
     ELSEIF(nb_arg == 5) then
      WRITE(message_text,'(a,2a6,5(1x,F9.6))') '  > ATOM: ', trim(atom_label(nb_atom)),trim(atom_typ(nb_atom)),  &
	                                           (atom_coord(i,nb_atom),i=1,3), atom_Biso(nb_atom), atom_occ_perc(nb_atom)
     endif
    endif
    call write_info(TRIM(message_text))

   ! WRITE(message_text,*) ' 2: ', atom_label(nb_atom),atom_typ(nb_atom)
   ! call write_info(trim(message_text))
   ! WRITE(message_text,*) ' 3: ', TRIM(input_line)
   ! call write_info(trim(message_text))
   ! pause
   !stop


   CASE ('HKL' )
    IF(nb_arg <3) then
     call check_arg('3', 'HKL')
     return
    endif
    IF(mode_interactif) then
     nb_hkl = 1
    else
     nb_hkl = nb_hkl + 1
    endif
    var(1:10) = 0.
    read(arg_string(1:3), *, IOSTAT=i_error)  var(1:3)
    IF(i_error /=0) then
     call error_message('HKL')
     return
    endif


    H(1,nb_hkl) = var(1)
    H(2,nb_hkl) = var(2)
    H(3,nb_hkl) = var(3)
    WRITE(message_text,'(a,3F6.2)') '  > HKL: ', H(1,nb_hkl), H(2,nb_hkl), H(3,nb_hkl)
    call write_info(TRIM(message_text))

    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('   CELL parameters have to be known for d_hkl calculation ...')
     call write_info('')
    endif

   CASE ('SF_HKL', 'SFAC_HKL')
    IF(nb_arg <3) then
     call check_arg('3', 'SFAC_HKL')
     return
    endif
    IF(mode_interactif) then
     nb_hkl_SFAC_calc = 1
    else
     nb_hkl_SFAC_calc = nb_hkl_SFAC_calc + 1
    endif
    var(1:10) = 0.
    
    read(arg_string(1:3), *, IOSTAT=i_error)  var(1:3)
    IF(i_error /=0) then
     call error_message('SFAC_HKL')
     return
    endif


    H(1,nb_hkl_SFAC_calc) = var(1)
    H(2,nb_hkl_SFAC_calc) = var(2)
    H(3,nb_hkl_SFAC_calc) = var(3)
    WRITE(message_text,'(a,3F6.2)') '  > Structure factor calculation for : ', &
	                                H(1,nb_hkl_SFAC_calc), H(2,nb_hkl_SFAC_calc), H(3,nb_hkl_SFAC_calc)
    call write_info(TRIM(message_text))

    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('   CELL parameters have to be known for structure factor calculation ...')
     call write_info('')
    endif
    keyword_SFAC_HKL = .true.


   CASE ('QVEC', 'Q_VEC', 'Q_VECTOR', 'KVEC', 'K_VEC', 'K_VECTOR', 'MOD_VEC', 'MODULATION_VECTOR')
    IF(nb_arg <3) then
     call check_arg('3', 'QVEC')
     return
    endif
    
    ! recherche forme fractionnaire
    i = index(arg_line, '/')
    if(i/=0) then
     call GET_QVEC_coord(arg_line)
    else
     var(1:10) = 0.
     READ(arg_string(1:3), *, IOSTAT=i_error) var(1:3)
     IF(i_error /=0) then
      call error_message('QVEC')
      return
     endif
     qvec(1) = var(1)
     qvec(2) = var(2)
     qvec(3) = var(3)
    endif 
    keyword_QVEC = .true.
    WRITE(message_text, '(a,3F6.2)') '  > QVEC: ', qvec(1:3)
    call write_info(TRIM(message_text))



   CASE ('SYMM', 'SYM', 'SYMMETRY_OPERATOR')
    IF(nb_arg ==0) then
     call check_arg('0', 'SYMM')
     return
    endif

    ! operateurs de symétrie sous la forme matricielle
    nb_symm_op = nb_symm_op + 1

    arg_line = l_case(arg_line)
    i_pos = INDEX(arg_line, 'x')
    arg_line = ADJUSTL(arg_line)

    if (i_pos==0) THEN  ! operateurs de symétrie sous la forme matricielle
     READ(arg_line, *, IOSTAT=i_error) (symm_op_rot (1,i,nb_symm_op),i=1,3), symm_op_trans(1,nb_symm_op),   &
                                       (symm_op_rot (2,i,nb_symm_op),i=1,3), symm_op_trans(2,nb_symm_op),   &
                                       (symm_op_rot (3,i,nb_symm_op),i=1,3), symm_op_trans(3,nb_symm_op)
     IF(i_error/=0)  then
      call error_message('SYMM')
      return
     endif

     keyword_SYMM   = .true.
     symm_op_xyz    = .false.
     symm_op_mat    = .true.

    else            ! operateurs de symétrie sous la forme x,y,z
     !call remove_car(arg_line, ' ')   ! remove blank in string
     arg_line = remove_car(arg_line, ' ')
     !call replace_car(arg_line, ',', ' ')
     arg_line = replace_car(arg_line,  ',', ' ')

     i1=INDEX(arg_line, ' ')
     symm_op_string(1,nb_symm_op) = arg_line(1:i1-1)
     i2=INDEX(TRIM(arg_line), ' ' , back=.TRUE.)
     symm_op_string(2,nb_symm_op) = arg_line(i1+1:i2-1)
     symm_op_string(3,nb_symm_op) = arg_line(i2+1:)
     keyword_SYMM   = .true.
     symm_op_xyz    = .true.
     symm_op_mat    = .false.
    end if

   case default

    unknown_CFL_keyword = .true.


  end select

 return
end subroutine identification_CFL_keywords

!----------------------------------------------------------------------------------------------






!----------------------------------------------------------------
subroutine check_matrice(M)
 USE cryscalc_module, ONLY : Mat_integer, Mat_det, debug_proc
 USE matrix_module,   ONLY : matrix_determinant_33
 implicit none
  REAL, DIMENSION(3,3), INTENT(IN)      :: M
  REAL, parameter                       :: EPS=0.001
  INTEGER                               :: i, j, n

  if(debug_proc%level_2)  write(debug_proc%unit, '(a)') " . routine : CHECK_MATRICE"
  
  ! verifie si la matrice (3,3) est entiere ou reelle

  Mat_integer = .false.
  n=0

  do i=1, 3
   do j=1,3
    if (ABS(M(i,j) - floor(M(i,j)+.5)) < eps ) n=n+1
   end do
  end do

  if (n==9) Mat_integer = .true.

  Mat_det = matrix_determinant_33(M)


 RETURN
end subroutine check_matrice

!--------------------------------------------------------------------

subroutine get_atom_coord(input_line)
 USE cryscalc_module, ONLY : atom_coord, atom_Biso, atom_occ, atom_occ_perc, nb_atom, debug_proc
 use macros_module,   only : nombre_de_colonnes
 implicit none
  CHARACTER (LEN=*), INTENT(INOUT)            :: input_line
  CHARACTER (LEN=32)                          :: string

  integer, parameter                          :: nb_fraction = 11
  CHARACTER (LEN=12), DIMENSION(nb_fraction)  :: ratio_string_neg, ratio_string_pos
  real,               DIMENSION(nb_fraction)  :: ratio_real_neg,   ratio_real_pos
  INTEGER                                     :: i, nb_col

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_ATOM_COORD")
  
  ratio_string_pos(1:nb_fraction) = (/'1/2',                            &
                                      '1/3', '2/3',                     &
                                      '1/4', '3/4',                     &
                                      '1/6', '5/6',                     &
                                      '1/8', '3/8', '5/8', '7/8'/)

  ratio_string_neg(1:nb_fraction) = (/'-1/2',                           &
                                      '-1/3', '-2/3',                   &
                                      '-1/4', '-3/4',                   &
                                      '-1/6', '-5/6',                   &
                                      '-1/8', '-3/8', '-5/8', '-7/8'/)

  ratio_real_pos(1:nb_fraction) = (/ 1./2.,                     &
                                     1./3., 2./3.,              &
                                     1./4., 3./4.,              &
                                     1./6., 5./6.,              &
                                     1./8., 3./8., 5./.8,  7./8. /)


  ratio_real_neg(1:nb_fraction) = -ratio_real_pos(1:nb_fraction)


  ! coord. x
  input_line = ADJUSTL(input_line)
  i = INDEX(input_line, ' ')

  string = input_line(1:i-1)
  string = ADJUSTL(string)
  input_line = input_line(i+1:)

  do i = 1, nb_fraction
   IF(string(1:3) == ratio_string_pos(i)(1:3)) then
    atom_coord(1, nb_atom) = ratio_real_pos(i)
    exit

   ELSEIF(string(1:4) == ratio_string_neg(i)(1:4)) then
    atom_coord(1, nb_atom) = ratio_real_neg(i)
    exit

   else
    READ(string, *) atom_coord(1, nb_atom)
   endif
  end do

  ! coord y
  input_line = ADJUSTL(input_line)
  i = INDEX(input_line, ' ')
  string = input_line(1:i-1)
  string = ADJUSTL(string)
  input_line = input_line(i+1:)

  do i = 1, nb_fraction
   IF(string(1:3) == ratio_string_pos(i)(1:3)) then
    atom_coord(2, nb_atom) = ratio_real_pos(i)
    exit

   ELSEIF(string(1:4) == ratio_string_neg(i)(1:4)) then
    atom_coord(2, nb_atom) = ratio_real_neg(i)
    exit

   else
    READ(string, *) atom_coord(2, nb_atom)
   endif
  end do



  ! coord z
  input_line = ADJUSTL(input_line)
  i = INDEX(input_line, ' ')
  string = input_line(1:i-1)
  string = ADJUSTL(string)
  do i = 1, nb_fraction
   IF(string(1:3) == ratio_string_pos(i)(1:3)) then
    atom_coord(3, nb_atom) = ratio_real_pos(i)
    exit

   ELSEIF(string(1:4) == ratio_string_neg(i)(1:4)) then
    atom_coord(3, nb_atom) = ratio_real_neg(i)
    exit

   else
    READ(string, *) atom_coord(3, nb_atom)
   endif
  end do

  ! Biso? occ ?
  input_line = ADJUSTL(input_line)
  i = INDEX(input_line, ' ')
  string = input_line(i:)
  string = ADJUSTL(string)
  call nombre_de_colonnes(string, nb_col)
  if(nb_col == 1) then
   read(string, *) atom_Biso(nb_atom)
   atom_occ_perc(nb_atom) = 1.
  elseif(nb_col == 2) then
   read(string, *) atom_Biso(nb_atom), atom_occ_perc(nb_atom)
  endif

 return
end subroutine get_atom_coord
!--------------------------------------------------------------------

subroutine get_qvec_coord(input_line)
 USE cryscalc_module, ONLY : qvec, debug_proc
 use macros_module,   only : nombre_de_colonnes
 USE Definition_fractions
 implicit none
  CHARACTER (LEN=*), INTENT(INOUT)            :: input_line
  CHARACTER (LEN=32)                          :: string
  INTEGER                                      :: i, k, nb_col

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_QVEC_COORD")
  
  call Def_fractions
 
  
  do k=1, 3
   input_line = ADJUSTL(input_line)
   i = INDEX(input_line, ' ')

   string = input_line(1:i-1)
   string = ADJUSTL(string)
   if(string(1:1) == '+') string = string(2:)
   
   do i = 1, nb_fraction
    IF(string(1:3) == ratio_string_pos(i)(1:3)) then
     Qvec(k) = ratio_real_pos(i)
     exit

    ELSEIF(string(1:4) == ratio_string_neg(i)(1:4)) then
     Qvec(k) = ratio_real_neg(i)
     exit

    else
     READ(string, *) Qvec(k)
    endif
   end do
   input_line = input_line(i+1:)
  end do
  
 return
end subroutine get_qvec_coord

!-----------------------------------------------------------------
subroutine decode_CHEM_string(arg_line, nb_arg)
 USE cryscalc_module, ONLY : nb_atoms_type, SFAC_type, SFAC_number, STO, keyword_ZUNIT, Z_UNIT, &
                             keyword_CHEM, keyword_CONT, keyword_SFAC_UNIT, debug_proc
 USE MACROS_module,   ONLY : check_character, nombre_de_colonnes
 USE IO_module,       ONLY : write_info

 implicit none
  CHARACTER(LEN=*),   INTENT(INOUT)  :: arg_line
  INTEGER,            INTENT(INOUT)  :: nb_arg
  CHARACTER(len=256)                 :: new_line, new_arg_line
  !local variables
  CHARACTER (LEN=64), DIMENSION(20)  :: arg_string
  INTEGER                            :: i, j
  CHARACTER (LEN=1)                  :: car
  LOGICAL                            :: alpha_char, numeric_char
  !INTEGER                            :: n_alpha, n_num, n_blank

  if(debug_proc%level_2)  call write_debug_proc_level(2, "DECODE_CHEM_STRING")
  
!! Insere des blancs si necessaire dans la chaine arg_line
     if(nb_arg == 1) call Insert_blank(arg_line, nb_arg)
	  
!    n_alpha = 0
!    n_num   = 0
!    n_blank = 0
!    new_arg_line = arg_line
!    if(nb_arg == 1) then
!     do i = 1, len_trim(arg_line)
!      car = arg_line(i:i)
!      call check_character(car, alpha_char, numeric_char)
!      if(alpha_char)   n_alpha = n_alpha + 1
!      if(numeric_char) n_num   = n_num + 1
!      if(i==1 .and. numeric_char) then
!       call write_info('   ... WRONG CHEMICAL FORMULA ...')
!       call write_info('')
!       return       
!      endif
!      if(alpha_char .and.  n_num /=0) then
!      write(*,*) ' i  = ', i
!       write(*,*) ' n_num  = ', n_num
!       write(*,*) ' n_alpha = ', n_alpha
!       write(*,*) ' >> ' , new_arg_line(1:i-1)
!       write(*,*) '    ' , trim(arg_line(i:))
!       write(new_line, '(3a)') new_arg_line(1:i+n_blank-1), ' ', trim(arg_line(i:))
!         write(*,*) ' >> ', trim(new_line)
!         new_arg_line = new_line
!         write(*,*) ' >> ', trim(new_arg_line)
!   
!       n_num   = 0
!       n_alpha = 0
!       n_blank = n_blank + 1
!      endif
!     end do
!     arg_line = new_arg_line
!     call nombre_de_colonnes(arg_line, nb_arg)
!     
!     write(*,*) ' >> ', trim(arg_line)
!    end if
     
!! ----------------------------------------------------------

   
    nb_atoms_type = nb_arg
    READ(arg_line,*) arg_string(1:nb_arg)

na: do i=1, nb_atoms_type
     do j= 1,LEN_TRIM(arg_string(i))
      car = arg_string(i)(j:j)
      alpha_char   = .false.
      numeric_char = .false.
      call check_character(car, alpha_char, numeric_char)
      IF(j==1 .AND. numeric_char) then
       call write_info('   ... WRONG CHEMICAL FORMULA ...')
       call write_info('')
       return
      endif

      IF(numeric_char) then
       READ(arg_string(i)(1:j-1),*) SFAC_type(i)        
       READ(arg_string(i)(j:),   *) STO(i)
       CYCLE na
      endif

      IF(j==LEN_TRIM(arg_string(i)) .AND. .NOT. numeric_char) then
       READ(arg_string(i)(1:j),*) SFAC_type(i)
       STO(i) = 1
       CYCLE na
      endif

     end do
    end do na

    if (.not. keyword_ZUNIT) then
     Z_unit = 1
     keyword_ZUNIT = .true.
    end if
    SFAC_number(1:nb_atoms_type) =  sto(1:nb_atoms_type) * Z_unit

    call get_Content()

    keyword_CHEM      = .true.
    keyword_CONT      = .false.
    keyword_SFAC_UNIT = .false.


 return
end subroutine decode_CHEM_string
!-----------------------------------------------------------------
subroutine get_content()
 USE cryscalc_module, ONLY : SFAC_type, SFAC_number, nb_atoms_type, molecule
 implicit none
  INTEGER             :: i

   molecule%content = ''
   do i=1, nb_atoms_type
    WRITE(molecule%content , '(4a,F6.2)') TRIM(molecule%content), ' ', TRIM(SFAC_type(i)), ' ',  SFAC_number(i)
   end do

 return
end subroutine get_content

!--------------------------------------------------------------------------------------------------
Subroutine Insert_blank(arg_line, nb_arg)
 USE cryscalc_module, only : debug_proc
 USE macros_module,   only : check_character, nombre_de_colonnes
 USE IO_module,       ONLY : write_info

 
  implicit none
   character (len=*), intent(inout) :: arg_line
   integer          , intent(inout) :: nb_arg
   integer                          :: n_alpha
   integer                          :: n_num
   integer                          :: n_blank
   character (len=256)              :: new_line, new_arg_line
   integer                          :: i
   character (len=1)                :: car
   logical                          :: alpha_char, numeric_char
   
   if(debug_proc%level_2)  call write_debug_proc_level(2, "INSERT_BLANK")
   
!! Insere des blancs si necessaire dans la chaine arg_line
!! ex: AA11BB22C33  ==>> AA11 BB22 C33

    n_alpha = 0
    n_num   = 0
    n_blank = 0
    new_arg_line = arg_line

	do i = 1, len_trim(arg_line)
     car = arg_line(i:i)
     call check_character(car, alpha_char, numeric_char)
     if(alpha_char)   n_alpha = n_alpha + 1
     if(numeric_char) n_num   = n_num + 1
     if(i==1 .and. numeric_char) then
      call write_info('   ... WRONG CHEMICAL FORMULA ...')
      call write_info('')
      return       
     endif
     if(alpha_char .and.  n_num /=0) then
      write(new_line, '(3a)') new_arg_line(1:i+n_blank-1), ' ', trim(arg_line(i:))
      new_arg_line = new_line
       n_num   = 0
       n_alpha = 0
       n_blank = n_blank + 1
      endif
     end do
     arg_line = new_arg_line
     nb_arg   = n_blank + 1 
     
     
    return
end subroutine Insert_blank     
!! ----------------------------------------------------------
