!     Last change:  TR   11 Oct 2007    6:53 pm

 ! lecture fichier .CIF
 !  source: EDPCR (JGP)

! lecture d'un fichier CIF

subroutine read_CIF_input_file(input_file, input_string)
 USE macros_module
 USE IO_module,                        ONLY : write_info
 USE cryscal_module
 USE CFML_String_Utilities,            ONLY : number_lines,   reading_lines
 USE CFML_IO_Formats,                  ONLY : Read_CIF_title, read_CIF_cell, read_cif_hm, read_cif_hall,  &
                                              read_cif_symm,  Read_Cif_Cont, Read_Cif_atom
 USE CFML_Crystal_Metrics,             ONLY : Set_Crystal_Cell, convert_u_betas, U_Equiv
 USE CFML_Crystallographic_Symmetry,   ONLY : set_spacegroup,   Space_Group_Type, get_hallsymb_from_gener,  &
                                              Get_Multip_Pos
 USE CFML_Atom_TypeDef,                ONLY : Atom_list_Type,  Deallocate_atom_list

 implicit none
  CHARACTER (LEN=*), INTENT(IN)            :: input_file
  CHARACTER (len=*), INTENT(IN)            :: input_string
  INTEGER                                  :: i, long_input_string

  ! local variable for .CIF file
  integer                                    :: nb_lines, npos
  character(len=80),dimension(:),allocatable :: file_cif
  character(len=40),dimension(48)            :: car_symop

  integer                                    :: n_elem_atm    ! N. of different species
  !real            ,  dimension(15)           :: n_elem        ! Number of elementos into the same species
  character(len=2),  dimension(15)           :: elem_atm      ! Character to identify the specie
  character(len=40)                          :: symb_spgr
  
  logical                                    :: input_out

  long_input_string = len_trim(input_string)
  input_out = .true.
  if(long_input_string == 6) then
   if(input_string(1:6) == 'NO_OUT') input_out = .false.
  endif
   
  call number_lines(trim(input_file),nb_lines)
  if (nb_lines ==0) then
   call write_info(' > No lines could be read. Program will be stopped')
   stop
  end if

  if (allocated(file_cif)) deallocate(file_cif)
  allocate(file_cif(nb_lines))
  file_cif=' '

 !---- Cargando Fichero en variable ----!
  call reading_lines(trim(input_file), nb_lines, file_cif)

  !---- TITL ----!
  npos=1
  call Read_Cif_Title(file_cif, npos, nb_lines, main_title)
  if(input_out .and. ON_SCREEN) then
   call write_info(' ' )
   call write_info('  . TITL: '//trim(main_title))
  endif 

  !---- CELL ----!
  npos=1
  call Read_Cif_Cell(file_cif, npos, nb_lines, unit_cell%param, unit_cell%param_esd)
  known_cell_esd = .true.
  call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)
  if(.not. input_out) return
  IF(unit_cell%volume < 0.1) call volume_calculation('no_out')  ! << oct. 2011
  
  keyword_CELL  = .true.
  if(ON_SCREEN) then
   call write_info(' ')
   write(message_text,'(a,6F10.4)') '  . CELL: ', (unit_cell%param(i), i=1,6)
   call write_info(TRIM(message_text))
   write(message_text,'(a,6F10.4)') '          ', (unit_cell%param_esd(i), i=1,6)
   call write_info(TRIM(message_text))
  endif 


  !---- OBTAIN SPACE GROUP ----!
  symb_spgr=' '
  npos=1
  call read_cif_hm(file_cif ,npos, nb_lines, symb_spgr)
  if (len_trim(symb_spgr) > 0) then
   call Set_SpaceGroup(symb_spgr, SPG)
  else
   npos=1
   call read_cif_hall(file_cif ,npos, nb_lines, symb_spgr)
   if (len_trim(symb_spgr) > 0) then
    call Set_SpaceGroup(symb_spgr, SPG)
   else
    npos=1
    call read_cif_symm(file_cif, npos ,nb_lines, nb_symm_op, car_symop)
    symb_spgr=' '
    call set_spacegroup(symb_spgr, SPG, car_symop, nb_symm_op,'gen')
    call get_hallsymb_from_gener(SPG)
   end if
  end if

  space_group_symbol = SPG%Spg_Symb
  !space_group_multip = SPG%multip
  IF(LEN_TRIM(space_group_symbol) /=0) then
   keyword_SPGR   = .true.
   WRITE_SPG_info = .true.
   if(ON_SCREEN) then
    call write_info(' ')
    call write_info('  . SPACE GROUP: '//TRIM(space_group_symbol))
   endif 
  else
   keyword_SPGR   = .false.
   WRITE_SPG_info = .false.
   if(ON_SCREEN) then
    call write_info(' ')
    call write_info('     >> Unknown space group !!')
   endif
  endif


  !---- ATOMS ----!
  npos=1
  call Read_Cif_Cont(file_cif, npos, nb_lines, n_elem_atm, elem_atm)

  npos=1
  call Read_Cif_Atom(file_cif, npos, nb_lines, nb_atom, Atm_list)

  IF(nb_atom /=0) then
   if(ON_SCREEN) then
    call write_info(' ')
    call write_info('  . ATOM LIST: ')
    call write_info(' ')
   endif

   do i=1,nb_atom
    atom_label(i)     = Atm_list%atom(i)%lab
    atom_typ(i)       = Atm_list%atom(i)%ChemSymb
    atom_coord(1:3,i) = Atm_list%atom(i)%x
    !atom_mult(i)      = Atm_list%atom(i)%mult
    atom_mult(i)      = Get_Multip_Pos(atom_coord(1:3,i), SPG)
    atom_occ(i)       = Atm_list%atom(i)%occ * real (Get_Multip_Pos(atom_coord(1:3,i) ,SPG)) / real(SPG%multip)
    atom_occ_perc(i)  = Atm_list%atom(i)%occ 

    if (Atm_list%atom(i)%thtype=='aniso') then
 !    n_t(i)=2
     atom_adp_aniso(1:6,i) = Atm_list%atom(i)%u(1:6)
     !atom_adp_aniso(1:6,i) = convert_u_betas(atom_adp_aniso(1:6,i), crystal_cell)
     atom_adp_equiv(i)     = 8.0*pi*pi*U_Equiv(crystal_cell, Atm_list%atom(i)%u(1:6))
    else
 !    n_t(i)=0
      !atom_adp_equiv(i)    = Atm_list%atom(i)%u(7)*8.0*pi*pi
      atom_adp_equiv(i)    = Atm_list%atom(i)%Ueq
	  atom_adp_equiv(i)    = 8.0*pi*pi*Atm_list%atom(i)%Ueq
    end if

!     atom_label(i)     = Atm_list%atom(i)%lab
!     atom_typ (i)      = Atm_list%atom(i)%ChemSymb
!     atom_coord(1:3,i) = Atm_list%atom(i)%x
!     atom_occ(i)       = Atm_list%atom(i)%occ
!     atom_mult(i)      = Atm_list%atom(i)%mult
!     atom_mult(i)      = Get_Multip_Pos(atom_coord(1:3,i), SPG)
!
!     atom_U(i)         = Atm_list%atom(i)%u(7)
     !atom_Biso(i)      = Atm_list%atom(i)%u(7)*8.0*pi*pi
	 atom_Biso(i)      = atom_adp_equiv(i) 

     if(ON_SCREEN) then
      !write(message_text,'(2x,2(1x,a),4F10.5,I4)')  atom_label(i)(1:6), atom_typ (i)(1:6) , &
	  !                                              atom_coord(1:3,i), atom_occ_perc(i), atom_mult(i)
      !call write_info(TRIM(message_text))
	  if(Atm_list%atom(i)%thtype=='aniso') then
       write(message_text,'(2x,2(1x,a6),2x,4F10.5,I4,6F10.5)')  atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
	                                                            atom_occ_perc(i), atom_mult(i), Atm_list%atom(i)%u(1:6)
      else
	  write(message_text,'(2x,2(1x,a6),2x,4F10.5,I4, F10.5)')  atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
	                                                            atom_occ_perc(i), atom_mult(i), Atm_list%atom(i)%Ueq
      endif	  
      call write_info(TRIM(message_text))
     endif 

   end do
  endif

  call Deallocate_atom_list(Atm_list)
  if (allocated(file_cif)) deallocate(file_cif)

  ! pour compatibilite avec CRYSCAL
   nb_atoms_type                = n_elem_atm
   SFAC_type(1:nb_atoms_type)   = elem_atm(1:n_elem_atm)
   !SFAC_number(1:nb_atoms_type) = n_elem(1:nb_atoms_type)
   

return
end subroutine read_CIF_input_file
!------------------------------------------------------------------------------

subroutine read_CIF_input_file_TR(input_unit)
 USE cryscal_module, ONLY : keyword_WAVE, wavelength, keyword_SIZE, Z_unit,          &
                            keyword_ZUNIT, CIF_diffrn_reflns, CIF_cell_measurement,  &
                            crystal, CIF_parameter, unit_cell, nb_atoms_type, molecule
 USE macros_module,  ONLY : nombre_de_colonnes, remove_car
 implicit none
  INTEGER, INTENT(IN)       :: input_unit
  CHARACTER (LEN=256)       :: string_value
  INTEGER                   :: nb_arg, i1, i2
  LOGICAL                   :: ok

  ! wave:
  call get_champ_value(input_unit, '_diffrn_radiation_wavelength', string_value, ok)
  IF(ok) then
   keyword_WAVE = .true.
   READ(string_value, *) wavelength
   call incident_beam
  endif

  ! SIZE
  call get_champ_value(input_unit, '_exptl_crystal_size_max', string_value, ok)
  IF(ok) then
   IF(string_value(1:1) /= '?')  then
    READ(string_value, *)  crystal%size(1)
    READ(string_value, *)  crystal%size_max
   endif
   call get_champ_value(input_unit, '_exptl_crystal_size_mid', string_value, ok)
   IF(ok) then
    IF(string_value(1:1)/= '?') then
     READ(string_value, *) crystal%size(2)
     READ(string_value, *) crystal%size_mid
    endif
    call get_champ_value(input_unit, '_exptl_crystal_size_min', string_value, ok)
    IF(ok) then
     IF(string_value(1:1) /= '?') then
      READ(string_value, *) crystal%size(3)
      READ(string_value, *) crystal%size_min
     endif
     keyword_SIZE = .true.
    endif
   endif
  endif

  IF(crystal%size(1) < 0.000001 .AND.  crystal%size(2) < 0.000001  .AND. crystal%size(3) < 0.000001) keyword_size = .false.



 ! ZUNIT
  call get_champ_value(input_unit,'_cell_formula_units_Z',        string_value, ok)
  IF(ok) then
   keyword_ZUNIT = .true.
   READ(string_value, *) Z_unit
  endif
  molecule%Z_unit = int(Z_unit)

 ! chemical formula
  call get_champ_value(input_unit,'_chemical_formula_sum',        string_value, ok)
  IF(ok) then
   !call remove_car(string_value, "'")
   string_value = remove_car(string_value, "'")
   call nombre_de_colonnes(string_value, nb_arg)
   call decode_CHEM_string(string_value, nb_arg)
   molecule%formula = string_value
  endif
  
  !call Get_CIF_value(input_unit, string_value, molecule%formula)


 ! THmin, THmax


  call get_champ_value(input_unit, '_cell_measurement_reflns_used',    string_value, ok)
  if (ok) READ(string_value, *) CIF_cell_measurement%reflns_used

  call get_champ_value(input_unit, '_cell_measurement_theta_min',    string_value, ok)
  if (ok) READ(string_value, *) CIF_cell_measurement%theta_min
  call get_champ_value(input_unit, '_cell_measurement_theta_max',    string_value, ok)
  if (ok) READ(string_value, *) CIF_cell_measurement%theta_max

  call get_champ_value(input_unit, '_diffrn_reflns_theta_min',    string_value, ok)
  if (ok) READ(string_value, *) CIF_diffrn_reflns%theta_min
  call get_champ_value(input_unit, '_diffrn_reflns_theta_max',    string_value, ok)
  if (ok) READ(string_value, *) CIF_diffrn_reflns%theta_max



! cell parameters and volume

  call get_champ_value(input_unit, '_cell_length_a', string_value, ok)
  IF(ok) READ(string_value, *) CIF_parameter%CELL_length_a

  call get_champ_value(input_unit, '_cell_length_b', string_value, ok)
  IF(ok) READ(string_value, *) CIF_parameter%CELL_length_b

  call get_champ_value(input_unit, '_cell_length_c', string_value, ok)
  IF(ok) READ(string_value, *) CIF_parameter%CELL_length_c

  call get_champ_value(input_unit, '_cell_angle_alpha', string_value, ok)
  IF(ok) READ(string_value, *) CIF_parameter%CELL_angle_alpha

  call get_champ_value(input_unit, '_cell_angle_beta', string_value, ok)
  IF(ok) READ(string_value, *) CIF_parameter%CELL_angle_beta

  call get_champ_value(input_unit, '_cell_angle_gamma', string_value, ok)
  IF(ok) READ(string_value, *) CIF_parameter%CELL_angle_gamma

  call get_champ_value(input_unit, '_cell_volume', string_value, ok)               
  if(ok) read(string_value, *) CIF_parameter%Cell_volume                           
  if(CIF_parameter%Cell_volume(1:1) /='?') then
   i1 = INDEX(string_value, '(')
   i2 = INDEX(string_value, ')')
   IF(i2 > i1) then
    READ(string_value(1:i1-1),   *) unit_cell%volume
    READ(string_value(i1+1:i2-1),*) unit_cell%volume_ESD
   endif
  endif


  ! temperature
  call get_champ_value(input_unit, '_cell_measurement_temperature', string_value, ok)
  if (ok) READ(string_value, *) CIF_cell_measurement%temperature

  return
end subroutine read_CIF_input_file_TR

!------------------------------------------------------------------------------
subroutine get_cell_parameters_from_cif_file(file_unit, cell_param, cell_param_esd)
 USE IO_module
 implicit none
  INTEGER, INTENT(IN)                  :: file_unit
  REAL,    INTENT(OUT),   dimension(*) :: cell_param, cell_param_esd
  character (len=256)                  :: line
  integer                              :: ier
  character (len=64)                   :: cell_a_string, cell_b_string, cell_c_string
  character (len=64)                   :: cell_alfa_string, cell_beta_string, cell_gamma_string
  LOGICAL , DIMENSION(6)               :: logical_cell


  logical_cell(1:6) = .false.
  REWIND(UNIT=file_unit)
  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info(' ... No cell parameters in the import.cif file !!')
    call write_info('')
    stop
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading import.cif file !!')
    call write_info('')
    stop
   endif

   if(trim(line(1:14))=='_cell_length_a') then
    read(line(16:),'(a)') cell_a_string
    call get_cell_param(cell_a_string, cell_param(1), cell_param_esd(1))
    logical_cell(1) = .true.

   elseif(trim(line(1:14))=='_cell_length_b') then
    read(line(16:),'(a)') cell_b_string
    call get_cell_param(cell_b_string, cell_param(2), cell_param_esd(2))
    logical_cell(2) = .true.

   elseif(trim(line(1:14))=='_cell_length_c') then
    read(line(16:),'(a)') cell_c_string
    call get_cell_param(cell_c_string, cell_param(3), cell_param_esd(3))
    logical_cell(3) = .true.

   elseif(trim(line(1:17))=='_cell_angle_alpha') then
    read(line(19:),'(a)') cell_alfa_string
    call get_cell_param(cell_alfa_string, cell_param(4), cell_param_esd(4))
    logical_cell(4) = .true.

   elseif(trim(line(1:16))=='_cell_angle_beta') then
    read(line(18:),'(a)') cell_beta_string
    call get_cell_param(cell_beta_string, cell_param(5), cell_param_esd(5))
    logical_cell(5) = .true.

   elseif(trim(line(1:17))=='_cell_angle_gamma') then
    read(line(19:),'(a)') cell_gamma_string
    call get_cell_param(cell_gamma_string, cell_param(6), cell_param_esd(6))
    logical_cell(6) = .true.
   endif

   IF(logical_cell(1) .AND. logical_cell(2) .AND. LOGICAL_cell(3)) then
    if(logical_cell(4) .AND. logical_cell(5) .AND. LOGICAL_cell(6)) then
     !call write_cell_parameters(cell_param)
     exit
    endif
   endif
  end do



 RETURN
end subroutine get_cell_parameters_from_cif_file


!----------------------------------------------------------------------

subroutine get_cell_param(cell_string, cell_value, cell_value_esd)
 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: cell_string
  REAL              , INTENT(OUT)       :: cell_value
  REAL              , INTENT(OUT)       :: cell_value_esd
  INTEGER                               :: i1, i2



  ! lecture de la valeur avant la parenthese
  i1 = INDEX (cell_string, '(')
  if (i1 ==0) then
   READ(cell_string, *) cell_value
  else
   READ(cell_string(1:i1-1), *) cell_value
   i2 = INDEX(cell_string, ')')
   READ(cell_string(i1+1:i2-1),*) cell_value_esd
  end if



 RETURN
end subroutine get_cell_param

!---------------------------------------------------------------------------

subroutine get_wavelength_from_cif_file(file_unit, wave)
 use IO_module
 implicit none
  INTEGER, INTENT(IN)                  :: file_unit
  REAL,    INTENT(OUT)                 :: wave
  character (len=256)                  :: line
  CHARACTER (LEN=256)                  :: champ_text
  REAL                                 :: champ_value
  integer                              :: ier
  LOGICAL                              :: logical_wave


  logical_wave = .false.

  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info(' ... No wavelength in the import.cif file !!')
    call write_info('')
    return
   ELSEIF(ier >0) then
    call write_info(' ')
    call write_info(' ... Problem reading wavelength value in import.cif file !!')
    call write_info('')
    return
   endif

   line = ADJUSTL(line)

   ! cas d'un fichier import.CIF cree par Nonius
   !if(trim(line(1:31))=='_diffrn_radiation_wavelength_id') then
   ! READ(UNIT=file_unit, '(a)', IOSTAT=ier)line
   ! if(ier/=0) then
   !  call write_info('')
   !  call write_info(' ... Problem reading import.cif file !!')
   !  call write_info('')
   !  stop
   ! endif
   ! read(line,*) wave
   ! exit
   !endif

   if(trim(line(1:28))=='_diffrn_radiation_wavelength') then
    READ(line, *, IOSTAT=ier ) champ_text, champ_value
    !if(ier/=0) then
    ! call write_info('')
    ! call write_info(' ... Problem reading wavelength in import.cif file !!')
    ! call write_info('')
    ! stop
    !endif
    IF(ier/=0) then
     ! cas d'un fichier import.cif cree par Nonius
     READ(file_unit, '(a)', IOSTAT=ier) line
     IF(ier/=0) then
      call write_info('')
      call write_info(' ... Problem reading wavelength !!')
      call write_info('')
      stop
     endif
     if(trim(line(1:31))=='_diffrn_radiation_wavelength_id') then
      READ(file_unit, '(a)', IOSTAT=ier)line
      if(ier/=0) then
       call write_info('')
       call write_info(' ... Problem reading import.cif file !!')
       call write_info('')
       stop
      END if
     endif
     read(line,*) wave
     exit
    endif

    IF(champ_value > 0.01) then
     wave = champ_value
     exit
    endif
   endif


  END do

 return
end subroutine get_wavelength_from_cif_file
!--------------------------------------------------------------------------

subroutine get_crystal_system_from_CIF_file(file_unit, crystal_system)
 USE macros_module, ONLY : l_case
 USE IO_module
 implicit none
  INTEGER,           INTENT(IN)        :: file_unit
  character (len=16), INTENT(OUT)       :: crystal_system
  character (len=256)                  :: line
  CHARACTER (LEN=256)                  :: champ_text, champ_value
  integer                              :: ier
  LOGICAL                              :: logical_wave

  rewind(file_unit)
  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info(' ... No crystal system in the import.cif file !!')
    call write_info('')
    return
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading crystal system in import.cif file !!')
    call write_info('')
    return
   endif

   line = ADJUSTL(line)

   
   if(trim(line(1:22))=='_symmetry_cell_setting') then
    READ(line, *, IOSTAT=ier ) champ_text, champ_value
    IF(ier/=0) then
      call write_info('')
      call write_info(' ... Problem reading crystal_system !!')
      call write_info('')
      return
    endif  
    
    exit
   endif 
  END do
  champ_value = l_case(champ_value)

  if(champ_value(1:9) == 'triclinic')         then
   crystal_system = "TRIC"
  elseif(champ_value(1:10) == 'monoclinic')   then
   crystal_system = "MONO"
  elseif(champ_value(1:12) == 'orthorhombic') then
   crystal_system = "ORTHO"
  elseif(champ_value(1:10) == 'tetragonal')   then
   crystal_system = "TETRA"
  elseif(champ_value(1:8)  == 'trigonal')     then
   crystal_system = "TRIG"
  elseif(champ_value(1:9)  == 'hexagonal')    then
   crystal_system = "HEXA"
  elseif(champ_value(1:5)  == 'cubic')        then
   crystal_system = "CUB"
  endif



end subroutine get_crystal_system_from_CIF_file


!--------------------------------------------------------------------------

subroutine get_H_M_from_CIF_file(file_unit, H_M)
 USE macros_module,  ONLY : l_case
 USE IO_module
 implicit none
  INTEGER,           INTENT(IN)        :: file_unit
  character (len=32), INTENT(OUT)      :: H_M
  character (len=256)                  :: line
  CHARACTER (LEN=256)                  :: champ_text, champ_value
  integer                              :: ier
  LOGICAL                              :: logical_wave

  rewind(file_unit)
  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info(' ... No Hermann Mauguin space group in the import.cif file !!')
    call write_info('')
    return
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading crystal system in import.cif file !!')
    call write_info('')
    return
   endif

   line = ADJUSTL(line)

   
   if(trim(line(1:30))=='_symmetry_space_group_name_H-M') then
    READ(line, *, IOSTAT=ier ) champ_text, champ_value
    IF(ier/=0) then
      call write_info('')
      call write_info(' ... Problem reading crystal_system !!')
      call write_info('')
      return
    endif  
    
    exit
   endif 
  END do
  
  H_M = adjustl(champ_value(1:32))



end subroutine get_H_M_from_CIF_file

!--------------------------------------------------------------------------

subroutine get_champ_value(CIF_unit, string_text, string_value, ok)
 implicit none
 INTEGER          , INTENT(IN)    :: CIF_unit
 CHARACTER (LEN=*), INTENT(IN)    :: string_text
 CHARACTER (LEN=*), INTENT(OUT)   :: string_value
 logical          , INTENT(INOUT) :: ok
 CHARACTER (LEN=256)              :: READ_line
 integer                          :: ier, long_string, long_line


 REWIND(UNIT=CIF_unit)

 ok = .false.
 do
  READ(CIF_unit, '(a)', IOSTAT=ier) read_line
  if(ier/=0) return

  read_line = ADJUSTL(read_line)
  long_string = LEN_TRIM(string_text)
  long_line   = len_trim(read_line)
  
  IF(READ_line(1:long_string) == TRIM(string_text)) then
   if(long_line == long_string) then  ! la valeur du champ n'est pas sur la même ligne
    read(CIF_unit, '(a)', iostat=ier) read_line
	if(ier/=0) return
	read(read_line, '(a)') string_value
   else
     READ(read_line(long_string+1:), '(a)') string_value
   endif
   string_value = ADJUSTL(string_value)
   IF(LEN_TRIM(string_value) == 0)   return
   IF(string_value(1:1)      == '?') return
   ok = .true.
   return
  endif
 ENDDO
              

 return
end subroutine get_champ_value
!--------------------------------------------------------------------------

subroutine read_lines_before_hkl(file_unit, cos_exist)
 implicit none
  INTEGER, INTENT(IN)             :: file_unit
  LOGICAL, INTENT(OUT)            :: cos_exist
  ! local variables
  INTEGER                         :: i, num_ligne, ier
  CHARACTER(LEN=128)              :: line


  num_ligne = 0
  rewind(UNIT=file_unit)

!  do
!   read(UNIT=file_unit, '(a)', iostat=ier) line
!   if(ier /=0) exit
!   num_ligne = num_ligne +1
!   if (line(1:22) =='_refln_F_squared_sigma' ) then
!    READ(UNIT=file_unit,'(a)',IOSTAT=Ier) (line, i=1,7)
!    if (line(1:34) == '_refln_nonius_diffracted_cos_cstar' .or. line(1:34) == '_nonius_refln_diffracted_cos_cstar') then
!     cos_exist = .true.
!     num_ligne = num_ligne + 7
!     exit
!    else
!     cos_exist = .false.
!     exit
!    end if
!   endif
!  end do

  cos_exist = .false.
loop_1:  do
   read(file_unit, '(a)', iostat=ier) line
   if(ier /=0) exit
   num_ligne = num_ligne +1
   if (line(1:22) =='_refln_F_squared_sigma' ) then
    do
     READ(file_unit, '(a)', IOSTAT=ier) line
     IF(line(1:1)/='_') exit loop_1
     num_ligne = num_ligne + 1
     if (line == '_refln_nonius_diffracted_cos_cstar' .or. line == '_nonius_refln_diffracted_cos_cstar') cos_exist = .true.

    end do

   endif
  end do loop_1



  rewind(unit=file_unit)
  ! lecture des lignes d'entete
  read(file_unit, '(a)') (line, i=1,num_ligne)

 RETURN
END subroutine read_lines_before_hkl

!----------------------------------------------------------------------------
subroutine read_SQUEEZE_file
 use cryscal_module, only : SQUEEZE
 implicit none
  logical               :: file_exist
  integer               :: i_error
  character (len=256)   :: read_line

 inquire(file= trim(SQUEEZE%file), exist= file_exist)
 SQUEEZE%procedure = .false.
 if(.not. file_exist) return
 
 
 open(unit=SQUEEZE%unit, file=trim(SQUEEZE%file))
  do
   read(SQUEEZE%unit, '(a)', iostat= i_error) read_line
   if(i_error /=0) exit
   SQUEEZE%nb_lines = SQUEEZE%nb_lines + 1
   SQUEEZE%read_line(SQUEEZE%nb_lines) = read_line
  end do
  
 close (unit=SQUEEZE%unit)
 SQUEEZE%procedure = .true.

 call write_CIF_file("SQUEEZE")
 
 return
end subroutine read_SQUEEZE_file

! ----------------------------------------------------------------------------------
 subroutine read_and_modif_ARCHIVE(input_unit)
  use cryscal_module, only : CIF_unit, CIF_archive_unit, AUTHOR, DEVICE, SQUEEZE, IMPORT_CIF_unit, &
                             CIF_parameter, SADABS_ratio, CIFdep, ACTA, CIF_format80, include_RES_file, tmp_unit
  use macros_module,  only : u_case, test_file_exist, Get_Wingx_job
  use IO_module,      only : write_info
  implicit none
   integer, intent(in)               :: input_unit
   character (len=256)               :: read_line
   integer                           :: i_error, i, i1, long
   integer                           :: n_line 
   logical                           :: file_exist
   logical                           :: diffracto
   logical                           :: correction_absorption
   logical                           :: squeeze_results
   logical                           :: author_data
   
   character (len=256)               :: CIF_string, tmp_string
   character (len=16), dimension(16) :: CIF_str
   character (len=16), dimension(8)  :: CIF_str_fmt
   character (len=64)                :: fmt_, fmt_final
   character (len=32)                :: job
   character (len=256)               :: wingx_structure_dir
   character (len=256)               :: wingx_res_file
   
   
   ! caracteristiques de la boucle _atom_site
   integer                           :: CIF_atom_site_item_nb                ! nombre de champs dans la boucle _atom_site_ 
   integer                           :: CIF_atom_site_label_numor
   integer                           :: CIF_atom_site_type_symbol_numor      ! numero du champ
   integer                           :: CIF_atom_site_fract_x_numor
   integer                           :: CIF_atom_site_fract_y_numor
   integer                           :: CIF_atom_site_fract_z_numor
   integer                           :: CIF_atom_site_U_numor
   integer                           :: CIF_atom_adp_type_numor
   integer                           :: CIF_atom_site_calc_flag_numor        ! numero du champ
   integer                           :: CIF_atom_site_refinement_flags_numor ! numero du champ
   integer                           :: CIF_atom_site_loop_numor             ! numero de la ligne du champ _loop   
   integer                           :: CIF_atom_site_occupancy_numor
   integer                           :: CIF_atom_site_disorder_assembly_numor
   integer                           :: CIF_atom_site_disorder_group_numor
   
   ! caracteristiques de la boucle _atom_type_
   integer                           :: CIF_atom_type_loop_numor
   integer                           :: CIF_atom_type_item_nb                   ! nombre de champs dans la boucle
   integer                           :: CIF_atom_type_symbol_numor              ! numero du champ _symbol dans la boucle

   ! caracteristiques de la boucle _atom_site_aniso
   integer                           :: CIF_atom_site_aniso_item_nb
   integer                           :: CIF_atom_site_aniso_label_numor
   integer                           :: CIF_atom_site_aniso_loop_numor             ! numero de la ligne du champ _loop   

   logical                           :: linear
   logical                           :: contain_H
  
  

   fmt_      = ''
   fmt_final = ''
  
  rewind(unit=input_unit)
  close(unit = CIF_unit)
  open(unit = CIF_unit, file = 'cryscal.cif')
  
  close(unit=CIF_archive_unit)
  open(unit=CIF_archive_unit, file="archive_cryscal.cif")
  
  ! premiere lecture du fichier archive.cif
  n_line = 0
  do 
   read(input_unit, '(a)', iostat= i_error) read_line
   if(i_error < 0) then
    n_line = 0
    exit
   endif
   if(i_error /=0) exit
   n_line  = n_line + 1
   i1 = index(read_line , "SECTION 2. COMPOUND(S) DETAILS")
   if(i1 /=0) exit 
  end do


!------------------------------------------------------------------------------------
! sept. 09

  ! champ "_atom_site_symmetry_multiplicity" present ?
  ! erreur dans SHELXL: le champ  est '_atom_site_symetry_multiplicity' (1 seul m)
  ! au lieu de "_atom_site_symmetry_multiplicity"
  ! TR nov. 08 : erreur corrigee dans SHELXL
  ! par consequent, l'archive.CIF cree par WinGX contiendra ou pas ce champ
  ! et le champ _atom_site_refinement_flags sera en position 12 ou 13.
  !CIF_atom_site_item_nb                = 12
  !CIF_atom_site_calc_flag_numor        = 9
  !CIF_atom_site_refinement_flags_numor = 10
  !CIF_atom_site_type_symbol_numor      = 2
  
  
  !rewind(unit=input_unit)
  !do
  ! read(unit=input_unit, '(a)', iostat=i_error) read_line
  ! if(i_error<0) exit
  ! if(index(read_line, '_atom_site_symmetry_multiplicity')/=0) then    
  !  CIF_atom_site_item_nb                = 13
  !  CIF_atom_site_calc_flag_numor        = 10    
  !  CIF_atom_site_refinement_flags_numor = 11
  !  exit
  ! end if
  !end do
 
 ! >>> new sept. 09
  
  ! determination du nombre de champ dans la boucle _atom_site_
  call Get_CIF_champ_nb('_atom_site_label', CIF_atom_site_item_nb, CIF_atom_site_loop_numor)
     
  ! champ _atom_site_type_symbol
  call Get_CIF_champ_numor('_atom_site_type_symbol', CIF_atom_site_type_symbol_numor)

  ! champ _atom_site_fract_x,y,z
  call Get_CIF_champ_numor('_atom_site_fract_x', CIF_atom_site_fract_x_numor)
  call Get_CIF_champ_numor('_atom_site_fract_y', CIF_atom_site_fract_y_numor)
  call Get_CIF_champ_numor('_atom_site_fract_z', CIF_atom_site_fract_z_numor)
  
 ! champ _atom_site_U
  call Get_CIF_champ_numor('_atom_site_U_iso_or_equiv', CIF_atom_site_U_numor)
  
 ! champ _atom_site_adp_type
  call Get_CIF_champ_numor('_atom_site_adp_type', CIF_atom_adp_type_numor)

 ! champ _atom_site_calc_flag
  call Get_CIF_champ_numor('_atom_site_calc_flag', CIF_atom_site_calc_flag_numor)
  
 ! champ _atom_site_refinement_flags
  call Get_CIF_champ_numor('_atom_site_refinement_flag', CIF_atom_site_refinement_flags_numor)
  
 ! champ _atom_site_occupancy
  call Get_CIF_champ_numor('_atom_site_occupancy', CIF_atom_site_occupancy_numor)
  
 ! champ _atom_site_disorder_assembly
  call Get_CIF_champ_numor('_atom_site_disorder_assembly', CIF_atom_site_disorder_assembly_numor)
  
 ! champ _atom_site_disorder_group
  call Get_CIF_champ_numor('_atom_site_disorder_group', CIF_atom_site_disorder_group_numor)
  
  ! determination du nombre de champ dans la boucle _atom_site_aniso
  call Get_CIF_champ_nb('_atom_site_aniso_label', CIF_atom_site_aniso_item_nb, CIF_atom_site_aniso_loop_numor)
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  ! molecular moiety
  rewind(unit = input_unit) 
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   long = len_trim('_chemical_formula_moiety') 
   if(index(read_line, '_chemical_formula_moiety') /=0) then
    read(read_line(long+1:),'(a)') CIF_parameter%formula_moiety
    CIF_parameter%formula_moiety = adjustl(CIF_parameter%formula_moiety)    
    if(CIF_parameter%formula_moiety(1:1) == "'") CIF_parameter%formula_moiety = CIF_parameter%formula_moiety(2:)
    long = len_trim(CIF_parameter%formula_moiety)
    if(CIF_parameter%formula_moiety(long:long) == "'") CIF_parameter%formula_moiety = CIF_parameter%formula_moiety(1:long-1)
    exit
   end if
  end do
  
  
!-----------------------------------------------------------------------------------------------  
! sept 09
!

!  rewind(unit=input_unit)
!  ! determination du _atom_site_refinement_flags pour les atomes d'hydrogen
!  do 
!   read(unit=input_unit, '(a)', iostat=i_error) read_line
!   if(i_error < 0) exit
!   if(index(read_line, '_atom_site_disorder_group') /=0) then
!    do 
!     read(unit=input_unit, '(a)', iostat=i_error) read_line
!     if(i_error < 0) exit
!     if(len_trim(read_line) == 0) exit
!     write(*,*) ' CIF_atom_site_item_numor : ', CIF_atom_site_item_nb
!     write(*,*) ' read_line  : ', trim(read_line)
!     read(read_line, *) CIF_str(1:CIF_atom_site_item_nb)
!     ! CIF_str(CIF_atom_site_type_symbol_numor)  : _atom_site_type_symbol
!     ! CIF_str(9) : _atom_site_calc_flag
!     if(CIF_str(CIF_atom_site_type_symbol_numor)(1:1) == 'H' .and. CIF_str(CIF_atom_site_calc_flag_numor)(1:1) /= 'c')  then
!      CIF_parameter%H_treatment = 'mixed'
!      exit
!     endif
!    end do
!    exit 
!   endif  
!  end do 

! >>> new sept. 09
!
! determination du _atom_site_refinement_flags pour les atomes d'hydrogene
  rewind (unit=input_unit)
  ! lecture des lignes avant la liste des atomes
  do i = 1, CIF_atom_site_loop_numor + CIF_atom_site_item_nb
   read(input_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit   
  end do
  
  contain_H = .false.
  do 
   read(input_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit 
   if(len_trim(read_line) == 0) exit
   read(read_line, *) CIF_str(1:CIF_atom_site_item_nb)
 
   if(len_trim(CIF_str(CIF_atom_site_type_symbol_numor)) == 1 .and. &
               CIF_str(CIF_atom_site_type_symbol_numor)(1:1) == 'H') then
      contain_H = .true.
   else
   !  contain_H = .false.
      cycle
   endif         

   if(len_trim(CIF_str(CIF_atom_site_type_symbol_numor)) == 1 .and. &
      CIF_str(CIF_atom_site_type_symbol_numor)(1:1) == 'H'    .and. &
      CIF_str(CIF_atom_site_calc_flag_numor)(1:1)   /= 'c')  then
      CIF_parameter%H_treatment = 'mixed'      
   endif     
  end do

  
!-------------------------------------------------------------------------------  
  
  ! diffracto ?
  diffracto = .false.
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   if(index(read_line, '_diffrn_measurement_device_type')/=0) then
    diffracto = .true.
    exit
   endif
  end do
  
  ! correction d'absorption ?
  correction_absorption = .false.
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   if(index(read_line, '_exptl_absorpt_process_details')/=0) then
    correction_absorption = .true.
    exit
   endif
  end do
  
  
  ! squeeze  ?
  squeeze_results = .false.
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   if(index(read_line, 'SQUEEZE RESULTS')/=0) then
    squeeze_results = .true.
    exit
   endif
  end do

  ! author data  ?
  author_data = .false.
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error < 0) exit
   if(index(read_line, 'AUTHORS DETAILS')/=0) then
    author_data = .true.
    exit
   endif
  end do


  ! lecture à partir de la ligne contenant "SECTION 2. COMPOUND(S) DETAILS"
  rewind(unit=input_unit)
  if(n_line /=0) then
   do i=1, n_line
    read(input_unit, '(a)', iostat = i_error) read_line    
   end do
  endif 
 
 
 
  do
   read(input_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit
   if(.not. author_data) then
    if(index(read_line, '#  CIF files read  :' ) /=0) then
     call write_CIF(CIF_archive_unit, trim(read_line))
     if(CIFDEP) then
      call write_CIF_AUTHOR(CIF_archive_unit)
      if (ACTA) call write_CIF_text(CIF_archive_unit)
     endif
     cycle
    endif
   endif
   
   
   
   if(index(read_line, '_exptl_absorpt_coefficient_mu') /=0 .and.  &
      .not. correction_absorption) then
      ! les info. sur les corrections d'absorption sont recuperees dans le fichier CRYSCAL.cif
    do
     read(CIF_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
    
     if(index(read_line, '_exptl_absorpt_coefficient_mu') /=0) then
      write(CIF_archive_unit, '(a)') trim(read_line)
      !do i = 1,7
      do
       read(cif_unit, '(a)', iostat=i_error) read_line
       if(i_error /=0) exit
        if(read_line(1:11) == '#----------') exit
       write(CIF_archive_unit, '(a)') trim(read_line)
      end do
     endif
    end do
 
    cycle   
   end if

   if(index(read_line, '_exptl_absorpt_correction_type') /=0) then
    if(.not. SQUEEZE_results .and. SQUEEZE%procedure) then
     write(CIF_archive_unit, '(a)') ""
     WRITE(CIF_archive_unit, '(a)') "#----------------------------------------------------------------------------#"
     WRITE(CIF_archive_unit, '(a)') "#                   SQUEEZE RESULTS                                          #"
     WRITE(CIF_archive_unit, '(a)') "#----------------------------------------------------------------------------#"
     write(CIF_archive_unit, '(a)') ""
     do i=2, SQUEEZE%nb_lines
      write(CIF_archive_unit, '(a)') trim(SQUEEZE%read_line(i))
     end do
     cycle
    else
     cycle
    endif 
   endif


   !if(.not. squeeze_results) then
   ! if(index(read_line, '_exptl_absorpt_correction_type') /=0) then
   !  if(SQUEEZE%procedure)  then
   !   write(unit=CIF_archive_unit, '(a)') ""
   !   WRITE(UNIT=CIF_archive_unit, '(a)') "#----------------------------------------------------------------------------#"
   !   WRITE(unit=CIF_archive_unit, '(a)') "#                   SQUEEZE RESULTS                                          #"
   !   WRITE(UNIT=CIF_archive_unit, '(a)') "#----------------------------------------------------------------------------#"
   !   write(unit=CIF_archive_unit, '(a)') ""
   !   do i=2, SQUEEZE%nb_lines
   !    write(unit=CIF_archive_unit, '(a)') trim(SQUEEZE%read_line(i))
   !   end do
   !  endif
   !  cycle
   ! endif
   !else
   ! if(index(read_line, '_exptl_absorpt_correction_type') /=0) cycle
   !endif


   if(index(read_line, '_atom_sites_solution_hydrogens') /=0 .and. .not. contain_H) cycle
   
   if(index(read_line, '_refine_ls_hydrogen_treatment') /=0) then
    if(contain_H) then
     write(CIF_archive_unit, '(2a)') "_refine_ls_hydrogen_treatment           ",trim(CIF_parameter%H_treatment)
    endif
    cycle
   end if 
   !if(index(read_line, 'mix') /=0) then
   ! write(CIF_archive_unit, '(a)') "_refine_ls_hydrogen_treatment           constr"
   ! cycle
   !endif

   if(index(read_line, '_diffrn_ambient_temperature')   /= 0    .and.  &
      .not. diffracto) then
    if(u_case(DEVICE%diffracto(1:4)) == 'APEX') then
     write(CIF_archive_unit, '(2a)') "_diffrn_measurement_device_type         ", trim(CIF_parameter%diffrn_measurement_device_type)
     write(CIF_archive_unit, '(2a)') "_diffrn_measurement_method              ", trim(CIF_parameter%diffrn_measurement_method)
     write(CIF_archive_unit, '(2a)') "_diffrn_detector                        ", trim(CIF_parameter%diffrn_detector)
     write(CIF_archive_unit, '(a)') trim(read_line)
    elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_source                          ", trim(CIF_parameter%diffrn_source)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_detector                        ", trim(CIF_parameter%diffrn_detector)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_detector_area_resol_mean        ", trim(CIF_parameter%diffrn_detector_area_resol_mean)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_measurement_device              ", trim(CIF_parameter%diffrn_measurement_device)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_measurement_device_type         ", trim(CIF_parameter%diffrn_measurement_device_type)
     WRITE(CIF_archive_unit, '(2a)')"_diffrn_measurement_method              ", trim(CIF_parameter%diffrn_measurement_method)
     WRITE(CIF_archive_unit, '(a)') trim(read_line)
    endif 
    cycle
   endif
   
   if(index(read_line, '_computing_structure_solution') /=0) then
    write(CIF_archive_unit, '(2a)') "_computing_structure_solution           ",trim(CIF_parameter%computing_structure_solution)
    cycle
   endif
   
   if(index(read_line, '_computing_publication_material') /=0) then
    write(CIF_archive_unit, '(a)')     "_computing_publication_material          "
    WRITE(CIF_archive_unit, '(a)')     ";"
    WRITE(CIF_archive_unit, '(10x,a)') trim(CIF_parameter%computing_publication_material_1)
    WRITE(CIF_archive_unit, '(10x,a)') trim(CIF_parameter%computing_publication_material_2)
    WRITE(CIF_archive_unit, '(a)')     ";" 
    if(len_trim(adjustl(read_line)) == len_trim('_computing_publication_material')) then   
     read(input_unit, '(a)', iostat = i_error) read_line
     if(i_error /=0) exit
    endif
   
    cycle
   endif
  
   ! -------  juillet 2011 : add jobname.res
   if(index(read_line, '_refine_diff_density_rms') /= 0) then
    write(CIF_archive_unit, '(a)')   trim(read_line)
	write(CIF_archive_unit, '(a)')   ''
    if(include_RES_file) then
	 call Get_Wingx_job(job, wingx_structure_dir)
	 if(job(1:1) /= '?') then
	  write(wingx_res_file, '(4a)') trim(wingx_structure_dir), '\', trim(job), '.res'
	  call test_file_exist(trim(wingx_res_file), file_exist)
	  if(file_exist) then
	   call write_info('')
	   call write_info('  ... Include '//trim(wingx_res_file)//' SHELXL file.')
	   call write_info('')
	   write(CIF_archive_unit, '(a)')   '_iucr_refine_instructions_details'
	   write(CIF_archive_unit, '(a)')   ';'
	   open(unit=tmp_unit, file=trim(wingx_res_file))
	    write(CIF_archive_unit, '(2a)') ' .res file for SHELXL : ', trim(job)//'.res'
		write(CIF_archive_unit, '(a)')  '.................................................................'
	    do
	     read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
	     if(i_error /=0) exit
	     write(CIF_archive_unit, '(a)') trim(read_line)
		end do 
	   close(unit=tmp_unit)	
	   write(CIF_archive_unit, '(a)')   ';'
	  endif	
	 endif 
	end if 
	cycle
   endif
   ! -----------------------------------------------------------
   
   if(index(read_line, '_atom_type_scat_source') /=0)  then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:4)
     i1 = index(read_line, "'")
     read(read_line(i1:), '(a)') CIF_string
     write(CIF_archive_unit, '(2a4, 2a8, a)') CIF_str(1:4), trim(CIF_string)
    end do 
    write(CIF_archive_unit, '(a)') ''   
    cycle
   endif
   
   
 
  
!  >>> new sept. 09

   if(index(read_line, "_atom_site_label") /=0  .and. &
     len_trim(adjustl(read_line)) ==len_trim("_atom_site_label")) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do i=1, CIF_atom_site_item_nb - 1
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     write(CIF_archive_unit, '(a)') trim(read_line)
    end do
  
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit     
     read(read_line, *) CIF_str(1:CIF_atom_site_item_nb)
    
     long = len_trim(CIF_str(CIF_atom_site_refinement_flags_numor))     
     !if(CIF_str(CIF_atom_site_refinement_flags_numor)(long:long) == 'P') then
     if(index(CIF_str(CIF_atom_site_refinement_flags_numor), 'P') /=0) then
     ! write(CIF_archive_unit, '(a)') trim(read_line)
      
      !call get_fmt(4, CIF_str(3:6), fmt_, 0)      
      call get_fmt(4, CIF_str(CIF_atom_site_fract_x_numor : CIF_atom_site_fract_x_numor +3), fmt_, 0)
      if(CIF_atom_site_type_symbol_numor == 2) then  ! symbol en position 2
       if(CIF_atom_site_item_nb == 12) then
        write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a10,a5,a4,2a2)' 
       else
        write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a10,a3,a5,a4,2a2)' 
       endif
      elseif(CIF_atom_site_type_symbol_numor == CIF_atom_site_item_nb) then ! _symbol en derniere position
       if(CIF_atom_site_item_nb == 12) then
        write(fmt_final, '(3a)') '(a5,', trim(fmt_), ',a5,a10,a5,a4,2a2,a3)'  
       else
        write(fmt_final, '(3a)') '(a5,', trim(fmt_), ',a5,a10,a3,a5,a4,2a2,a3)' 
       endif
      endif
      !write(CIF_archive_unit, fmt_final) CIF_str(1:CIF_atom_site_item_nb)
   
     else
      !call get_fmt(4, CIF_str(3:6), fmt_, 0)      
      
      call get_fmt(4, CIF_str(CIF_atom_site_fract_x_numor : CIF_atom_site_fract_x_numor +3), fmt_, 0)
      if(CIF_atom_site_type_symbol_numor == 2) then  ! symbol en position 2
       if(CIF_atom_site_item_nb==12) then
        write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,a2,a5,a3,2a2)'        
       else
        write(fmt_final, '(3a)') '(a5,a3,', trim(fmt_), ',a5,2a2,a5,a3,2a2)'            
       endif
      elseif(CIF_atom_site_type_symbol_numor == CIF_atom_site_item_nb) then ! _symbol en derniere position
      
       if(CIF_atom_site_item_nb==12) then
        write(fmt_final, '(3a)') '(a5,', trim(fmt_), 'a5,a2,a5,a3,2a2,a3)'        
       else
        write(fmt_final, '(3a)')  '(a5,', trim(fmt_), 'a5,2a2,a5,a3,2a2,a3)'            
       endif      
      endif
            	 
     endif
	 
	 write(tmp_string, fmt_final) CIF_str(1:CIF_atom_site_item_nb)
     if(len_trim(tmp_string) > 80 .and. .not. CIF_format80) then
      write(CIF_archive_unit, '(a)') trim(read_line)
     else
      write(CIF_archive_unit, fmt_final) CIF_str(1:CIF_atom_site_item_nb)  
     endif          
   
     
     !if(CIF_str(12)(1:1) == '.') then
     ! call get_fmt(4, CIF_str(3:6), fmt_, 0)
     ! write(fmt_, '(3a)') '(a5,a3', trim(fmt_), ',a5,a2,a5,3a3)'        
     ! write(CIF_archive_unit, fmt_) CIF_str(1:12)
     !else ! cas d'un desordre traite avec PART 1 / PÄRT 2
     ! write(CIF_archive_unit, '(a)') trim(read_line)
     !endif
    end do
    write(CIF_archive_unit, '(a)') ''
   
   endif 
   

!  >> new sept. 09

    if(index(read_line, "_atom_site_aniso_label") /=0) then
     write(CIF_archive_unit, '(a)') trim(read_line)
     do i=1, CIF_atom_site_aniso_item_nb - 1
      read(input_unit, '(a)', iostat=i_error) read_line
      if(i_error /=0) exit
      write(CIF_archive_unit, '(a)') trim(read_line)
     end do

     do
      read(input_unit, '(a)', iostat=i_error) read_line
      if(i_error /=0) exit
      if(len_trim(read_line) == 0) exit
      read_line = adjustl(read_line)
      read(read_line, *) CIF_str(1:7)
      call get_fmt(6, CIF_str(2:7), fmt_, 0)
            
      write(fmt_final, '(3a)') '(a4,', trim(fmt_),')'    
      write(tmp_string, fmt_final) CIF_str(1:7)
      if(len_trim(tmp_string) > 80 .and. .not. CIF_format80) then
	   write(CIF_archive_unit, '(a)') trim(read_line)
	  else
       write(CIF_archive_unit, fmt_final) CIF_str(1:7)
      endif
     end do
     write(CIF_archive_unit, '(a)') ''
    
	endif
    
!-----------------------------------------------------------------------------------------------   
   
   
   if(index(read_line, '_geom_bond_publ_flag') /=0) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:5)
     write(CIF_archive_unit, '(2a5, a12,2a12)') CIF_str(1:5)
    end do
    write(CIF_archive_unit, '(a)') ''
   endif
 
   if(index(read_line, '_geom_angle_publ_flag') /=0) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:7)
     call get_fmt(1, CIF_str(4:4), fmt_, 100)
     write(fmt_final, '(3a)') '(3a5,', trim(fmt_), ',3a8)'
     write(CIF_archive_unit, fmt_final) CIF_str(1:7)
    end do
    write(CIF_archive_unit, '(a)') ''
   endif

 
   if(index(read_line, '_geom_torsion_publ_flag') /=0) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:10)
     
     !call verif_linearity(CIF_str(1:4), linear)  ! calculs.F90
     !if(linear) cycle
     
     call get_fmt(1, CIF_str(5:5), fmt_, -100)
     write(fmt_final, '(3a)') '(4a5', trim(fmt_), ',x,5a8)'
     write(CIF_archive_unit, fmt_final) CIF_str(1:10)     
    end do
    write(CIF_archive_unit, '(a)') ''
   endif

   if(index(read_line, '_geom_hbond_site_symmetry_A') /=0) then
    write(CIF_archive_unit, '(a)') trim(read_line)
    do
     read(input_unit, '(a)', iostat=i_error) read_line
     if(i_error /=0) exit
     if(len_trim(read_line) == 0) exit
     read(read_line, *) CIF_str(1:8)
     !call get_fmt(1, CIF_str(5:5), fmt_, -100)
     !write(fmt_final, '(3a)') '(4a5', trim(fmt_), ',5a6)'
     !write(CIF_archive_unit, fmt_) CIF_str(1:8)
     write(CIF_archive_unit, '(3a5, 4a12,a6)') CIF_str(1:8)
    end do
    write(CIF_archive_unit, '(a)') ''

   endif
   
    
   write(CIF_archive_unit, '(a)') trim(read_line)
  end do
 
  close(unit=CIF_unit)
  return
 end subroutine read_and_modif_ARCHIVE
 
 
!---------------------------------------------------------------------------
 subroutine Get_fmt(n, string, fmt_, level)
  implicit none
  integer,            intent(in)               :: n
  character (len=*),  dimension(n), intent(in) :: string
  character (len=64), intent(out)              :: fmt_
  integer,            intent(in)               :: level
  real                                         :: value
  integer                                      :: i, i1
  
  fmt_ = ''
  do i = 1, n
   i1 = index(string(i), '(')
   if(i1 ==0) then 
    read(string(i), *) value
   else
    read(string(i)(1:i1-1), *) value
   endif
   
   if(level == 0) then
    if(string(i)(1:1) == "-") then
     write(fmt_, '(2a)') trim(fmt_), ',a13'
    else
     if(i==1) then
      write(fmt_, '(2a)') trim(fmt_), '1x,a12'
     else
      write(fmt_, '(2a)') trim(fmt_), ',1x,a12'
     endif
    endif
    !if(value > 0.) then
    ! write(fmt_, '(2a)') trim(fmt_), ',1x,a13'
    !else
    ! write(fmt_, '(2a)') trim(fmt_), ',a14'
    !end if
    
   elseif(level == 100) then     
    if(value > 99.99) then
     write(fmt_, '(2a)') trim(fmt_), ',a13'
    else
     if(i==1) then
      write(fmt_, '(2a)') trim(fmt_), '1x,a12'
     else
      write(fmt_, '(2a)') trim(fmt_), ',1x,a12'
     endif
    endif 
    
   elseif(level == -100) then
    if(value < -99.99) then
     write(fmt_, '(2a)') trim(fmt_), ',a13'
    elseif(value < -9.99) then
     write(fmt_, '(2a)') trim(fmt_), ',1x,a12'
    elseif(string(i)(1:1) == '-') then
     write(fmt_, '(2a)') trim(fmt_), ',2x,a11'
    elseif(value < 9.99) then
     write(fmt_, '(2a)') trim(fmt_), ',3x,a10'
    elseif(value < 99.99) then
     write(fmt_, '(2a)') trim(fmt_), ',2x,a11' 
    else
     write(fmt_, '(2a)') trim(fmt_), ',1x,a12'
    endif 
   endif
  end do
  
  return
 end subroutine Get_fmt
 
!---------------------------------------------------------------------------
 subroutine Get_SADABS_ratio
  use  cryscal_module, only : SADABS_ratio, SADABS_Tmin, SADABS_Tmax, IMPORT_cif_unit
  use  macros_module,  only : test_file_exist
  implicit none
   integer                  :: i_error,i1
   character (len=256)      :: read_line
   logical                  :: file_exist
  
  ! lecture du fichier import.cif pour recuperer la valeur de Tmin/Tmax  
  close(unit=IMPORT_CIF_unit)
  open(unit = IMPORT_CIF_unit, file = 'import.cif', iostat=i_error)
  call test_file_exist("import.cif", file_exist)
  do
   read(IMPORT_CIF_unit, '(a)', iostat = i_error) read_line
   if(i_error /=0) exit
   read_line = adjustl(read_line)
   i1 = index(read_line, '# Ratio of minimum to maximum apparent transmission:')
   if (i1 /=0) then
    i1 = index(read_line, ':')
    read(read_line(i1+1:), *) SADABS_ratio
    if(SADABS_ratio < 0. .or. SADABS_ratio > 1.) SADABS_ratio = -999.
   endif 
  end do
  close(unit=IMPORT_CIF_unit)
  
  if(SADABS_ratio < 0.) then
   open(unit = IMPORT_CIF_unit, file = 'import.cif', iostat=i_error)
   call test_file_exist("import.cif", file_exist)
   do
    read(IMPORT_CIF_unit, '(a)', iostat = i_error) read_line
    if(i_error /=0) exit
    read_line = adjustl(read_line)
    i1 = index(read_line, '# Estimated minimum and maximum transmission:')
    if (i1 /=0) then
     i1 = index(read_line, ':')
     read(read_line(i1+1:), *) SADABS_Tmin, SADABS_Tmax
     SADABS_ratio = SADABS_Tmin / SADABS_Tmax
    endif 
   end do
   close(unit=IMPORT_CIF_unit)  
  endif

 

 return
end subroutine Get_SADABS_ratio

!----------------------------------------------------------------------------------
! recheche de la position d'un champ particulier dans une boucle loop_ du fichier.CIF

subroutine  Get_CIF_champ_numor(input_string, output_numor)
 use cryscal_module, only : input_unit
 implicit none
  character (len=*), intent(inout) :: input_string  
  integer,           intent(out)   :: output_numor
  CHARACTER (LEN=256)              :: READ_line
  integer                          :: i_error

  output_numor = 1
  rewind(unit=input_unit)
  do
   read(input_unit, '(a)', iostat=i_error) read_line
   if(i_error /=0) return
   if(index(read_line, trim(input_string)) /=0) then

    do 
     backspace(unit=input_unit)
     backspace(unit=input_unit)
     read(input_unit, '(a)', iostat=i_error) read_line
     
     if(i_error /=0) return
     read_line = adjustl(read_line)
     if(read_line(1:5) == 'loop_') exit
    end do
    do 
      read(input_unit, '(a)', iostat=i_error) read_line
      if(i_error /=0) return
      read_line = adjustl(read_line)
      if(index(read_line, trim(input_string)) /=0) exit
      output_numor = output_numor + 1
    end do     
    
   endif
  end do
  
   return
  end subroutine Get_CIF_champ_numor
  
!--------------------------------------------------------------------
subroutine Get_CIF_champ_nb(input_string, output_nb, output_nb2)
  use cryscal_module, only : input_unit
   implicit none
   character (len=*), intent(inout) :: input_string  
   integer,           intent(out)   :: output_nb, output_nb2
   CHARACTER (LEN=256)              :: READ_line
   integer                          :: i_error

   output_nb  = 0
   output_nb2 = 0
   rewind(unit = input_unit)
   do
    read(input_unit, '(a)', iostat=i_error) read_line
    if(i_error /=0) return
    output_nb2 = output_nb2 + 1
    read_line = adjustl(read_line)
    if(index(read_line, trim(input_string)) == 0) cycle
    
    do 
     backspace(unit=input_unit)
     backspace(unit=input_unit)
     output_nb2 = output_nb2 - 1
     read(input_unit, '(a)', iostat=i_error) read_line    
     if(i_error /=0) return
     read_line = adjustl(read_line)
     if(read_line(1:5) == 'loop_') exit
    end do
    
    do 
     read(input_unit, '(a)', iostat=i_error) read_line    
     if(i_error /=0) return
     read_line = adjustl(read_line)
     if(read_line(1:11) /= "_atom_site_") return
     output_nb = output_nb + 1     
    end do
    
   end do

   return
end subroutine Get_CIF_champ_nb


!-------------------------------------------------------------------------------------------------------

subroutine Create_CEL_from_CIF
! creation d'un fichier .CEL à partir d'un fichier.CIF

 USE cryscal_module, only : keyword_read_CIF, keyword_read_PCR, keyword_read_INS, message_text,  &
                           CEL_unit, INS_file_name, CIF_file_name, CEL_file_name, PCR_file_name, &
                           unit_cell, SPG, nb_atom, atom_coord, atom_typ, atom_occ, atom_occ_perc
 USE IO_module,      ONLY : write_info
 USE CFML_Scattering_Chemical_Tables, ONLY : chem_info, num_chem_info,   set_chem_info

 USE macros_module,                   ONLY : u_case

  
  implicit none
  integer                 :: i, j


  if (.not. keyword_read_CIF .and. .not. keyword_read_INS .and. .not. keyword_read_PCR) then  
   call write_info('')
   if(.not. keyword_read_CIF) then
    call write_info(' ! A CIF file has to be read first !')
   elseif(.not. keyword_read_INS) then
    call write_info(' ! A INS file has to be read first !')
   elseif(.not. keyword_read_PCR) then
    call write_info(' ! A PCR file has to be read first !')
   endif 
   call write_info('')
   return
  end if

  if(keyword_read_CIF) then
   i = index(CIF_file_name, '.', back=.true.)
   if (i/=0) then
    write(CEL_file_name, '(2a)') CIF_file_name(1:i-1), '_cif.CEL'
   else
    write(CEL_file_name, '(2a)') trim(CIF_file_name), '_cif.CEL'
   endif 
  elseif(keyword_read_INS) then
   i = index(INS_file_name, '.', back=.true.)
   if (i/=0) then
    write(CEL_file_name, '(2a)') INS_file_name(1:i-1), '_ins.CEL'
   else
    write(CEL_file_name, '(2a)') trim(INS_file_name), '_ins.CEL'
   endif 
  elseif(keyword_read_PCR) then  
   i = index(PCR_file_name, '.', back=.true.)
   if (i/=0) then
    write(CEL_file_name, '(2a)') PCR_file_name(1:i-1), '_pcr.CEL'
   else
    write(CEL_file_name, '(2a)') trim(PCR_file_name), '_pcr.CEL'
   endif 
  endif
 
  CALL set_chem_info

  open(unit = CEL_unit, file = trim(CEL_file_name))
    write(message_text, '(a,6F10.4)')   'CELL ', (unit_cell%param(i), i=1,6)
    write(CEL_unit, '(a)') adjustl(trim(message_text))
    if(nb_atom > 1) write(CEL_unit, '(a,i4)')     'natom ', nb_atom    
    do i=1, nb_atom 
     do j=1, Num_Chem_Info
       if (u_case(atom_typ(i)(1:))== u_case(chem_info(j)%symb(1:))) exit
     end do
     write(message_text,  '(a2,2x, i3,4F10.5)')  adjustl(atom_typ(i)), chem_info(j)%Z , atom_coord(1:3,i), atom_occ_perc(i)
     write(CEL_unit, '(a)') adjustl(trim(message_text))         
    end do 
    write(message_text, '(a,i3)')     'rgnr ',  SPG%NumSpg
    write(CEL_unit, '(a)') adjustl(trim(message_text))
   
     

  close(unit = CEL_unit)
  
  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(CEL_file_name), ' file has been created.'  
  call write_info(trim(message_text))
  call write_info('')
  return
  
 



end subroutine create_CEL_from_CIF  

!-------------------------------------------------------------------------------------------------------

subroutine Create_ACE_from_CIF
! creation d'un fichier .ACE pour CaRIne à partir d'un fichier.CIF

 USE cryscal_module, only : keyword_read_CIF, keyword_read_PCR, keyword_read_INS, message_text,  &
                           ACE_unit, INS_file_name, CIF_file_name, ACE_file_name, PCR_file_name, &
                           unit_cell, SPG, nb_atom, atom_coord, atom_typ, atom_occ, atom_occ_perc
 USE IO_module,      ONLY : write_info
 USE CFML_Scattering_Chemical_Tables, ONLY : chem_info, num_chem_info,   set_chem_info

 USE macros_module,                   ONLY : u_case

  
  implicit none
  integer                           :: i, j
  character (len=32), dimension(6)  :: cell_string
  character (len=64)                :: ACE_string
  character (len=32), dimension(3)  :: coord_string
  character (len=1)                 :: tab    
  real, parameter                   :: eps = 0.00001

  tab = char(9)  ! tabulation
  
  if (.not. keyword_read_CIF .and. .not. keyword_read_INS .and. .not. keyword_read_PCR) then  
   call write_info('')
   if(.not. keyword_read_CIF) then
    call write_info(' ! A CIF file has to be read first !')
   elseif(.not. keyword_read_INS) then
    call write_info(' ! A INS file has to be read first !')
   elseif(.not. keyword_read_PCR) then
    call write_info(' ! A PCR file has to be read first !')
   endif 
   call write_info('')
   return
  end if

  if(keyword_read_CIF) then
   i = index(CIF_file_name, '.', back=.true.)
   if (i/=0) then
    write(ACE_file_name, '(2a)') CIF_file_name(1:i-1), '_cif.ACE'
   else
    write(ACE_file_name, '(2a)') trim(CIF_file_name), '_cif.ACE'
   endif 
  elseif(keyword_read_INS) then
   i = index(INS_file_name, '.', back=.true.)
   if (i/=0) then
    write(ACE_file_name, '(2a)') INS_file_name(1:i-1), '_ins.ACE'
   else
    write(ACE_file_name, '(2a)') trim(INS_file_name), '_ins.ACE'
   endif 
  elseif(keyword_read_PCR) then  
   i = index(PCR_file_name, '.', back=.true.)
   if (i/=0) then
    write(ACE_file_name, '(2a)') PCR_file_name(1:i-1), '_pcr.ACE'
   else
    write(ACE_file_name, '(2a)') trim(PCR_file_name), '_pcr.ACE'
   endif 
  endif
 
  CALL set_chem_info
  
  

  open(unit = ACE_unit, file = trim(ACE_file_name))
   write(ACE_unit, '(a)')  'CaRIne Crystallography 3.0'
   write(ACE_unit, '(2a)') 'Cell Standard ASCII File=', trim(ACE_file_name) 
   write(ACE_unit, '(a)')  'File Version Number=1'
   write(ACE_unit, '(a)')  ''
   write(ACE_unit, '(a)')  'Colors Refered to: Mendeleïev Table'
   write(ACE_unit, '(a)')  ''
   write(ACE_unit, '(a)')  '----------------------'
   write(ACE_unit, '(a)')  'Cell Definition'
   write(ACE_unit, '(a)')  '----------------------'
   write(ACE_unit, '(a)')  ''
   write(ACE_unit, '(a)')  'Cell Parameters(Å and °)'
   do i=1, 6
    write(cell_string(i), '(F15.5)') unit_cell%param(i)
   end do
   write(ACE_unit, '(2a,a1,2a,a1,2a)')  'a=',trim(adjustl(cell_string(1))), tab,   &
                                        'b=',trim(adjustl(cell_string(2))), tab,   &
                                        'c=',trim(adjustl(cell_string(3)))

   write(ACE_unit, '(2a,a1,2a,a1,2a)')  'alpha=',trim(adjustl(cell_string(4))), tab,   &
                                        'beta=', trim(adjustl(cell_string(5))), tab,   &
                                        'gamma=',trim(adjustl(cell_string(6)))

   write(ACE_unit, '(a)')          '' 
   write( ACE_unit, '(a)')          'System=CEL'
   write(ACE_string, *) SPG%NumSpg
   write(ACE_unit, '(2a)')         'Space Group Number=',trim(adjustl(ACE_string))
   if(SPG%NumSpg == 146 .or. SPG%NumSpg == 148 .or. SPG%NumSpg == 155 .or. SPG%NumSpg == 160 .or. &
      SPG%NumSpg == 161 .or. SPG%NumSpg == 166 .or. SPG%NumSpg == 167) then   ! groupes trigonaux R
      if((abs(unit_cell%param(4) -  90.) < eps) .and. (abs(unit_cell%param(4) - 90.) < eps) .and.   &
         (abs(unit_cell%param(6) - 120.) < eps)) then
         write(ACE_unit, '(a)')          'Hexagonal description'
      else
         write( ACE_unit, '(a)')          'Trigonal description'
      endif   
      
   endif
   write(ACE_unit, '(a)')          'Only Non Equivalent Positions are listed'
   write(ACE_string, *) nb_atom
   write(ACE_unit, '(2a)')         'Number of positions in Cell=',trim(adjustl(ACE_string))
   write(ACE_unit, '(a)')          '' 
   write(ACE_unit, '(6(a,a1),a)')  'Atom',tab,'Oxi.',tab,'X',tab,'Y',tab,'Z',tab,'R(Å)',tab,'Occ.'
   do i=1,nb_atom
    write(ACE_string, '(F5.2)') atom_occ_perc(i)
    write(coord_string(1), '(F8.5)')  atom_coord(1:1,i)
    write(coord_string(2), '(F8.5)')  atom_coord(2:2,i)
    write(coord_string(3), '(F8.5)')  atom_coord(3:3,i)
    write(ACE_unit, '(6(a,a1),a)')  trim(adjustl(atom_typ(i))), tab, '0+', tab,         &
                                    trim(adjustl(coord_string(1))), tab,                 &
                                    trim(adjustl(coord_string(2))), tab,                 &
                                    trim(adjustl(coord_string(3))), tab,                 &
                                    '1',tab, trim(adjustl(ACE_string))
   end do

  close(unit = ACE_unit)

!    write(unit = message_text, '(a,6F10.4)')   'CELL ', (unit_cell%param(i), i=1,6)
!    write(unit = CEL_unit, '(a)') adjustl(trim(message_text))
!    if(nb_atom > 1) write(unit = CEL_unit, '(a,i4)')     'natom ', nb_atom    
!    do i=1, nb_atom 
!     do j=1, Num_Chem_Info
!       if (u_case(atom_typ(i)(1:))== u_case(chem_info(j)%symb(1:))) exit
!     end do
!     write(message_text,  '(a2,2x, i3,4F10.5)')  adjustl(atom_typ(i)), chem_info(j)%Z , atom_coord(1:3,i), atom_occ_perc(i)
!     write(unit = CEL_unit, '(a)') adjustl(trim(message_text))         
!    end do 
!    write(unit = message_text, '(a,i3)')     'rgnr ',  SPG%NumSpg
!    write(unit = CEL_unit, '(a)') adjustl(trim(message_text))
!   
!      

! close(unit = CEL_unit)
  
  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(ACE_file_name), ' file has been created.'  
  call write_info(trim(message_text))
  call write_info('')
 
!  call write_info('')
!  call write_info(' This routine is not yet implement in CRYSCAL. Sorry !')
!  call write_info('')
 
  return
  
 

end subroutine create_ACE_from_CIF  