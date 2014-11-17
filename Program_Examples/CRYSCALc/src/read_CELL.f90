!     Last change:  TR   12 Oct 2007    9:24 am

subroutine read_P4P(P4P_input_file)
! . selection du fichier.P4P
! . lecture fichier.P4P
! . creation P4P.cif
! . creation apex.cif

 USE IO_module
 USE cryscalc_module, ONLY : P4P_file_name, message_text, keyword_CELL, input_line, lecture_OK, debug_proc 
 implicit none
  CHARACTER (LEN=*), intent(IN)        :: P4P_input_file
  INTEGER                              :: i, ier, nb_P4P, selected_P4P
  CHARACTER (LEN=256)                  :: P4P_line
  CHARACTER (LEN=256), DIMENSION(100)  :: P4P_file
  REAL                                 :: esd

  if(debug_proc%level_2)  call write_debug_proc_level(2, "read_P4P")
  
  IF(LEN_TRIM(P4P_input_file) ==0) then
   call system('dir *.P4P /B > P4P.lst')
   OPEN(UNIT=66, FILE='P4P.lst')
    nb_P4P = 0
    do
     READ(66, '(a)', IOSTAT=ier) P4P_line
     IF(ier <0) EXIT ! fin du fichier
     nb_P4P = nb_P4P + 1
     IF(nb_P4P > 100) exit
     P4P_file(nb_P4P) = ADJUSTL(P4P_line)
    END do
   CLOSE(UNIT=66)

   IF(nb_P4P ==0) then
    call write_info('')
    call write_info('  !! No .P4P file in the current folder !!')
    call write_info('')
    return
   ELSEIF(nb_P4P > 100) then
    call write_info('')
    call write_info('  !! Number of .P4P file in the current folder is too large !!')
    call write_info('')
    return

   ELSEIF(nb_P4P == 1) then
    selected_P4P = 1
    call write_info('')
    WRITE(message_text, '(2a)') '  . P4P file: ', TRIM(P4P_file(selected_P4P))
    call write_info(TRIM(message_text))
   else
    call WRITE_info('')
    do i=1, nb_P4P
     WRITE(message_text,'(2x,I3,2a)') i, ': ', TRIM(P4P_file(i))
     call write_info(trim(message_text))
    END do
    call WRITE_info('')
    do
     call WRITE_info(' >> Select .P4P file: ')
     call read_input_line(input_line)
     READ(input_line,*) selected_P4P
     IF(selected_P4P > 0 .and. selected_P4P <= nb_P4P)  exit
    END do
   endif

   call system('del P4P.lst')

   P4P_file_name = P4P_file(selected_P4P)
  else
   P4P_file_name = P4P_input_file
  endif
  lecture_ok = .false.
  
  call read_P4P_file(TRIM(P4P_file_name))

  IF(.NOT. lecture_ok) then
   call write_info('')
   call write_info(' !! Problem reading .P4P file !!')
   call write_info('')
  endif
  
  !call volume_calculation('no_out')
 
  keyword_CELL = .true.


 return
end subroutine read_P4P
!-------------------------------------------------------------------------------------------------------------------

subroutine read_P4P_file(P4P_file)
 USE cryscalc_module,     ONLY : P4P_read_unit, known_cell_esd, wavelength, keyword_WAVE, unit_cell,       &
                                 message_text, keyword_SIZE, crystal, UB_matrix, UB_mat_log, molecule,     &
 	 						     CIF_string, lecture_OK, debug_proc							 
 USE CIF_module,          ONLY : CIF_cell_measurement, CIF_parameter								 
 USE macros_module,       ONLY : test_file_exist, nombre_de_colonnes, l_case, u_case
 USE IO_module,           ONLY : write_info
 USE CFML_Geometry_SXTAL, ONLY : CELL_Fr_UB


 implicit none
  CHARACTER (LEN=*), INTENT(IN)  :: P4P_file
  ! local variables
  CHARACTER (LEN=256)              :: P4P_line 
  INTEGER                          :: ier, i, i1, i2, long
  INTEGER                          :: nb_arg
  CHARACTER(LEN=16), DIMENSION(12) :: arg_string
  LOGICAL                          :: file_exist
  LOGICAL                          :: UB_mat_1, UB_mat_2, UB_mat_3

  
  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_P4P_FILE")
  
  arg_string = '?'
  lecture_ok = .false.
  file_exist = .false.
  P4P_line   = ''
  CIF_string = ''
  UB_mat_1   = .false.
  UB_mat_2   = .false.
  UB_mat_3   = .false.
  UB_mat_log = .false.

  call test_file_exist(TRIM(P4P_file),file_exist, 'out' )
  IF(.NOT. file_exist) return

  OPEN(UNIT= P4P_read_unit, FILE=TRIM(P4P_file))
   do
    READ(P4P_read_unit, '(a)', IOSTAT=ier) P4P_line
    IF(ier <0) EXIT ! fin du fichier
	long = len_trim(P4P_line)
	!if(P4P_line(1:5) == 'REF05') cycle
	if(P4P_line(1:5) == 'REF05') exit
	
	call nombre_de_colonnes(P4P_line, nb_arg)
	IF(nb_arg/=0)  then
     READ(P4P_line, *) arg_string(1:nb_arg)
	 arg_string(1:nb_arg) = adjustl(arg_string(1:nb_arg))
    else
	 cycle
	end if	
 
    select case (arg_string(1)) 
	    case ('FILEID')
		 i1 = index(P4P_line, ":", back=.true.)
		 if(i1 /=0)  read(P4P_line(i1+3:), *) CIF_parameter%sample_ID
		 CIF_parameter%sample_ID = adjustl(CIF_parameter%sample_ID)
	 
	   
        case ('CELL')
         READ(P4P_line(5:),*) unit_cell%param(1:6), unit_cell%volume
         lecture_ok = .true.
 
        case ('CELLSD')
         READ(P4P_line(7:), *) unit_cell%param_esd(1:6), unit_cell%volume_esd
         known_cell_esd = .true.

        case ('CSIZE')
         IF(TRIM(arg_string(2)) /= '?') then
          READ(P4P_line(6:), *) crystal%size(1:3)
          keyword_SIZE = .true.
          call get_crystal_SIZE_min_max
         endif

		CASE ('FACE')
         crystal%faces_nb = crystal%faces_nb + 1
         read(P4P_line(6:), '(a)') crystal%face_line(crystal%faces_nb)
		 
		CASE ('SAINGL')
         READ(P4P_line(7:), *) arg_string(1:3)
         READ(arg_string(1), *) CIF_cell_measurement%reflns_used
         READ(arg_string(2), *) CIF_cell_measurement%theta_min
         READ(arg_string(3), *) CIF_cell_measurement%theta_max
         !exit ! sinon la lecture du facies ne se fait pas 
		 
        
         
        case ('SOURCE')
         READ(P4P_line(7:), *) arg_string(1:2)
         READ(arg_string(2), *) wavelength
         keyword_WAVE = .true.

       

        case ('ORT1')
         READ(P4P_line(5:), *, iostat=ier) UB_matrix(1,1), UB_matrix(1,2), UB_matrix(1,3)
		 if(ier == 0) UB_mat_1 = .true.
		
        case ('ORT2')
         READ(P4P_line(5:), *, iostat=ier) UB_matrix(2,1), UB_matrix(2,2), UB_matrix(2,3)
		 if(ier == 0) UB_mat_2 = .true.

        case ('ORT3')
         READ(P4P_line(5:), *, iostat=ier) UB_matrix(3,1), UB_matrix(3,2), UB_matrix(3,3)
		 if(ier == 0) UB_mat_3 = .true.
		 if(UB_mat_1 .and. UB_mat_2 .and. UB_mat_3) UB_mat_log = .true.

		 
        case ('MORPH')
         READ(P4P_line(6:), *) crystal%morph

        case ('CCOLOR')
         READ(P4P_line(7:), *) crystal%color

        case ('CHEM')
         READ(P4P_line(5:), '(a)') molecule%formula
         molecule%formula = ADJUSTL(molecule%formula)

        case ('BRAVAIS')
         !!READ(P4P_line(8:), '(a)') unit_cell%Bravais
         !!unit_cell%Bravais = ADJUSTL(unit_cell%Bravais)
         !READ(P4P_line(8:), '(a)') CIF_string
         !CIF_string = ADJUSTL(CIF_string)
         !i1 = INDEX(CIF_string, '(')
         !i2 = INDEX(CIF_string, ')')
         !unit_cell%crystal_system = '?'
         !unit_cell%Bravais        = '?'
         !if (i1 > 1) then
         ! READ(CIF_string(1:i1-1), *) unit_cell%crystal_system
         ! IF(i2 /=0) READ(CIF_string(i2+1:), *) unit_cell%Bravais
         !else
         ! READ(CIF_string, *) unit_cell%crystal_system, unit_cell%Bravais
         !endif
		 
		 i1 = index(P4P_line, '(')
		 if(i1 /=0) then
		  read(P4P_line(8:i1-1), *)   unit_cell%crystal_system
		 else
          read(P4P_line(8:),     *)   unit_cell%crystal_system
         endif
         read(P4P_line(long:long), *) unit_cell%Bravais		 
		
		 !read(P4P_line(8:), *)        unit_cell%crystal_system
		 !read(P4P_line(long:long), *) unit_cell%Bravais
         unit_cell%crystal_system = l_case(unit_cell%crystal_system)
         unit_cell%crystal_system = ADJUSTL(unit_cell%crystal_system)
		 unit_cell%Bravais        = u_case(unit_cell%Bravais)
         unit_cell%Bravais        = ADJUSTL(unit_cell%Bravais)
		 
		 
		 ! creation de la chaine unit_cell%H_M
		 select case(unit_cell%crystal_system)
		  case ("cubic")
		   write(unit_cell%H_M, '(2a)') trim(unit_cell%Bravais),"23"  
		   
		  case ("tetragonal")		  
           write(unit_cell%H_M, '(2a)') trim(unit_cell%Bravais),"4"  
			
	      case ("trigonal", "rhomboedral")
           write(unit_cell%H_M, '(2a)') trim(unit_cell%Bravais),"3"  
   		   
	      case ("hexagonal")
		   if(arg_string(3) == 'RO') then
		    write(unit_cell%H_M, '(a)') "R3"  
		   else
            write(unit_cell%H_M, '(2a)') trim(unit_cell%Bravais),"6"  
		   endif	
			
		  case ("orthorhombic")
           write(unit_cell%H_M, '(2a)') trim(unit_cell%Bravais),"222"  
		  
		  case ("monoclinic")
           write(unit_cell%H_M, '(2a)') trim(unit_cell%Bravais),"2"  
		   
		  case("triclinic")
           write(unit_cell%H_M, '(2a)') trim(unit_cell%Bravais),"1"  
		 		  		
		 end select 
 
       case ('TITLE')
        READ(P4P_line(6:), '(a)') CIF_string
        CIF_string = ADJUSTL(CIF_string)
        IF(CIF_string(1:14) == 'Integration of') then
         READ(CIF_string(15:), *) molecule%common_name
         molecule%common_name = ADJUSTL(molecule%common_name)
		 !CIF_parameter%sample_ID = molecule%common_name
        endif
		
	   case ('REF05', 'DATA')
        exit	   

       case default
    end select
   END do



  close (UNIT=P4P_read_unit)
  lecture_OK = .true.
  
 RETURN
end subroutine read_P4P_file

!----------------------------------------------------------------------
subroutine read_M50_file(M50_file, lecture_ok)
 USE cryscalc_module, ONLY : M50_read_unit, unit_cell, known_cell_esd, wavelength, keyword_WAVE, debug_proc
 USE macros_module,   ONLY : test_file_exist

 implicit none
  CHARACTER (LEN=*), INTENT(IN)  :: M50_file
  LOGICAL,           INTENT(OUT) :: lecture_ok
  ! local variables
  CHARACTER (LEN=256)              :: M50_line
  INTEGER                          :: ier
  CHARACTER(LEN=16), DIMENSION(10) :: arg_string
  LOGICAL                          :: file_exist

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_M50_FILE")

  lecture_ok = .false.
  file_exist = .false.

  call test_file_exist(TRIM(M50_file),file_exist, 'out' )
  IF(.NOT. file_exist) return

  OPEN(UNIT= M50_read_unit, FILE=TRIM(M50_file))
   do
    READ(M50_read_unit, '(a)', IOSTAT=ier) M50_line
    IF(ier <0) EXIT ! fin du fichier
    READ(M50_line, *) arg_string(1)
    arg_string(1) = ADJUSTL(arg_string(1))

    select case (arg_string(1))
        case ('cell')
         READ(M50_line(5:),*) unit_cell%param(1:6)
         lecture_ok = .true.

        case ('esdcell')
         READ(M50_line(8:), *) unit_cell%param_esd(1:6)
         known_cell_esd = .true.
         exit

        case ('lambda')
         READ(M50_line(7:), *) wavelength
         keyword_WAVE = .true.

        case default
    end select
   END do



  close (UNIT=M50_read_unit)
 RETURN
end subroutine read_M50_file

!----------------------------------------------------------------------
subroutine read_RED_input_file(RED_file, lecture_ok)
 USE cryscalc_module,                  ONLY : RED_read_unit, unit_cell, known_cell_esd, wavelength, keyword_WAVE, SPG, &
                                              space_group_symbol, keyword_SPGR, debug_proc 
 USE macros_module,                    ONLY : test_file_exist, u_case
 USE CFML_Crystallographic_Symmetry,   ONLY : set_spacegroup

 implicit none
  CHARACTER (LEN=*), INTENT(IN)  :: RED_file
  LOGICAL,           INTENT(OUT) :: lecture_ok
  ! local variables
  CHARACTER (LEN=256)              :: RED_line
  INTEGER                          :: ier
  CHARACTER(LEN=16), DIMENSION(10) :: arg_string
  LOGICAL                          :: file_exist
  character(len=40)                :: symb_spgr
  
  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_RED_INPUT_FILE")

  lecture_ok = .false.
  file_exist = .false.

  call test_file_exist(TRIM(RED_file),file_exist, 'out' )
  IF(.NOT. file_exist) return

  OPEN(UNIT= RED_read_unit, FILE=TRIM(RED_file))
   do
    READ(RED_read_unit, '(a)', IOSTAT=ier) RED_line
    IF(ier <0) EXIT ! fin du fichier
    READ(RED_line, *) arg_string(1)
    arg_string(1) = ADJUSTL(arg_string(1))

    select case (u_case(arg_string(1)))
        case ('CELL')
         READ(RED_line(5:),*) unit_cell%param(1:6)
         lecture_ok = .true.        

        case ('WAVE')
         READ(RED_line(7:), *) wavelength
         keyword_WAVE = .true.

		case ('SPGR')
         read(RED_line(5:), *) symb_spgr		
		 call Set_SpaceGroup(symb_spgr, SPG)
		 keyword_SPGR   = .true.  
		 
        case default
    end select
   END do



  close (UNIT=RED_read_unit)
 RETURN
end subroutine read_RED_input_file


!----------------------------------------------------------------------
subroutine read_X_input_file(X_file, lecture_ok)
 USE cryscalc_module, ONLY : X_read_unit, unit_cell, wavelength, keyword_WAVE, debug_proc
 USE macros_module,   ONLY : test_file_exist

 implicit none
  CHARACTER (LEN=*), INTENT(IN)  :: X_file
  LOGICAL,           INTENT(OUT) :: lecture_ok
  ! local variables
  CHARACTER (LEN=256)              :: X_line
  INTEGER                          :: ier
  CHARACTER(LEN=16), DIMENSION(10) :: arg_string
  LOGICAL                          :: file_exist

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_X_INPUT_FILE")
  
  lecture_ok = .false.
  file_exist = .false.

  call test_file_exist(TRIM(X_file),file_exist, 'out' )
  IF(.NOT. file_exist) return

  OPEN(UNIT= X_read_unit, FILE=TRIM(X_file))
   do
    READ(X_read_unit, '(a)', IOSTAT=ier) X_line
    IF(ier <0) EXIT ! fin du fichier
    READ(X_line, *) arg_string(1)
    arg_string(1) = ADJUSTL(arg_string(1))

    select case (arg_string(1))
        case ('unit')
         READ(X_line(10:),*) unit_cell%param(1:6)
         lecture_ok = .true.

        case ('wavelength')
         READ(X_line(11:), *) wavelength
         keyword_WAVE = .true.

        case default
    end select
   END do



  close (UNIT=X_read_unit)
 RETURN
end subroutine read_X_input_file

!----------------------------------------------------------------------
subroutine read_RMAT_input_file(RMAT_file, lecture_ok)
 USE cryscalc_module, ONLY : RMAT_read_unit, unit_cell, debug_proc
 USE macros_module,   ONLY : test_file_exist

 implicit none
  CHARACTER (LEN=*), INTENT(IN)  :: RMAT_file
  LOGICAL,           INTENT(OUT) :: lecture_ok
  ! local variables
  CHARACTER (LEN=256)              :: RMAT_line
  INTEGER                          :: ier
  CHARACTER(LEN=16), DIMENSION(10) :: arg_string
  LOGICAL                          :: file_exist

  if(debug_proc%level_2)  write(debug_proc%unit, '(a)') " . routine : READ_RMAT_INPUT_FILE"
  
  lecture_ok = .false.
  file_exist = .false.

  call test_file_exist(TRIM(RMAT_file),file_exist, 'out' )
  IF(.NOT. file_exist) return

  OPEN(UNIT= RMAT_read_unit, FILE=TRIM(RMAT_file))
   do
    READ(RMAT_read_unit, '(a)', IOSTAT=ier) RMAT_line
    IF(ier <0) EXIT ! fin du fichier
    READ(RMAT_line, *) arg_string(1)
    arg_string(1) = ADJUSTL(arg_string(1))

    select case (arg_string(1))
        case ('CELL')
         READ(RMAT_line(5:),*) unit_cell%param(1:6)
         lecture_ok = .true.

        case default
    end select
   END do



  close (UNIT=RMAT_read_unit)
 RETURN
end subroutine read_RMAT_input_file

!------------------------------------------------------------

subroutine read_SADABS_file()
 USE Macros_module,   ONLY : test_file_exist
 USE cryscalc_module, ONLY : P4P_file_name, ABS_file_name, ABS_read_unit, SADABS_line_ratio, SADABS_line_estimated_Tmin_Tmax, &
                             SADABS_absorption_coeff, debug_proc
 USE IO_module

 implicit none
  INTEGER                       :: i1, i_error, long
  CHARACTER (LEN=256)           :: read_line
  LOGICAL                       :: file_exist

  if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_SADABS_FILE")
  
  SADABS_line_ratio               = '?'
  SADABS_line_estimated_Tmin_Tmax = '?'
  SADABS_absorption_coeff         = '?'
  i1 = INDEX(P4P_file_name, '_0M.P4P')
  IF(i1 == 0)   return
  

  ABS_file_name = P4P_file_name(1:i1-1)//'m.abs'
  call test_file_exist(trim(ABS_file_name), file_exist, 'no_out')
  
  if(file_exist) then
   OPEN(UNIT=ABS_read_unit, FILE=TRIM(ABS_file_name))
   do
    READ(ABS_read_unit, '(a)', IOSTAT=i_error) read_line
    IF(i_error /=0) exit
    read_line = ADJUSTL(read_line)
	
    long = LEN_TRIM('Ratio of minimum to maximum apparent transmission:')
    IF(read_line(1:long) == 'Ratio of minimum to maximum apparent transmission:') then
     SADABS_line_ratio = read_line
     exit
    endif
   end do
   return
   
  else
   ! mars 2011 : new version of SADABS
   ABS_file_name = P4P_file_name(1:i1-1)//'.abs'
   call test_file_exist(trim(ABS_file_name), file_exist, 'out')
   if(.not. file_exist)  return
   
  
   OPEN(UNIT=ABS_read_unit, FILE=TRIM(ABS_file_name))
   do
    READ(ABS_read_unit, '(a)', IOSTAT=i_error) read_line
    IF(i_error /=0) exit
    read_line = ADJUSTL(read_line)
	long = LEN_TRIM('Linear absorption coefficient set to')
    IF(read_line(1:long) == 'Linear absorption coefficient set to') then
	 i1 = index(read_line, 'mm-1')
	 if(i1 /=0) then
	  read(read_line(long+1:i1-1), '(a)') SADABS_absorption_coeff
	 endif
     !exit
    endif
	
	
    long = LEN_TRIM('Estimated minimum and maximum transmission:')
    IF(read_line(1:long) == 'Estimated minimum and maximum transmission:') then
     SADABS_line_estimated_Tmin_Tmax = read_line
     exit
    endif
   end do
  endif 
  
  call write_info('')
  call write_info('  . SADABS output file: '//trim(ABS_file_name))
  
 RETURN
end subroutine read_SADABS_file
