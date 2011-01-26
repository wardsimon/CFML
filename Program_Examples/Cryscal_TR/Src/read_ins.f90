!     Last change:  TR    2 Mar 2007    2:56 pm



subroutine read_INS_input_file(input_file, input_string)
 USE macros_module
 USE IO_module,                        ONLY : write_info
 USE cryscal_module
 USE wavelength_module
 USE CFML_String_Utilities,            ONLY : number_lines,   reading_lines
 USE CFML_IO_Formats,                  ONLY : Read_Shx_Titl, Read_Shx_Cell, Read_Shx_Latt, Read_Shx_Symm,   &
                                              Read_Shx_Cont, Read_Shx_Fvar, Read_Shx_Atom
 USE CFML_Crystal_Metrics,             ONLY : Set_Crystal_Cell
 USE CFML_Crystallographic_Symmetry,   ONLY : set_spacegroup,   Space_Group_Type, get_hallsymb_from_gener,  &
                                              Get_Multip_Pos
 USE CFML_Atom_TypeDef,                ONLY : Atom_list_Type,  Deallocate_atom_list

 USE SHELX_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)            :: input_file
  CHARACTER (LEN=*), INTENT(IN)            :: input_string 
  CHARACTER (LEN=256)                      :: read_line
  INTEGER                                  :: i, n, i1, i2, i_error
  integer                                  :: j, k, long
  INTEGER, DIMENSION(500)                  :: n_sl            ! nombre d'atomes avec le meme label

  REAL, DIMENSION(10)                      :: var
  character (len=10)                       :: tmp_string

  ! local variable for SHELX.INS file
  integer                                    :: nb_lines, npos
  character(len=80),dimension(:),allocatable :: fileshx
  character(len=40)                          :: car_gsp

  CHARACTER(LEN=48)                          :: op_string
  !type (Crystal_Cell_Type)                   :: crystal_cell


 ! lecture fichier SHELXL.INS
 ! read .INS input file for SHELX
 ! extract space group from the lattice and symmetry operators
 !  source: EDPCR (JGP)

  call number_lines(trim(input_file),nb_lines)
  
  if (nb_lines ==0) then
   call write_info(' > No lines could be read. Program will be stopped')
   stop
  end if

  if (allocated(fileshx)) deallocate(fileshx)
  allocate(fileshx(nb_lines))
  fileshx=' '

 !---- Cargando Fichero en variable ----!
  call reading_lines(trim(input_file), nb_lines, fileshx)

 !---- TITL ----!
  npos=1
  call Read_Shx_Titl(fileshx, npos, nb_lines, main_Title)
  if(input_string(1:6) /= 'NO_OUT') then 
   call write_info(' ' )
   call write_info('  . TITL: '//trim(main_title))
  endif 

 !---- CELL / ZERR ----!
  call Read_Shx_Cell(fileshx, npos, nb_lines, unit_cell%param, unit_cell%param_esd, wavelength, Z_unit_INS)
  known_cell_esd = .true.
  call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)
  keyword_CELL  = .true.
  keyword_WAVE  = .true.
  keyword_ZUNIT = .true.
 
  if(input_string(1:6) /= 'NO_OUT') then
   call write_info(' ')
   write(message_text,'(a,6F10.4)') '  . CELL: ', (unit_cell%param(i), i=1,6)
   call write_info(TRIM(message_text))
   write(message_text,'(a,6F10.4)') '          ', (unit_cell%param_esd(i), i=1,6)
   call write_info(TRIM(message_text))
   write(message_text,'(a,F10.5)')  '  . WAVE: ', wavelength
   call write_info(TRIM(message_text))
   write(message_text,'(a,I3)')     '  . Z:    ', Z_unit_INS
   call write_info(TRIM(message_text))
  endif 


  call volume_calculation('')
  !call incident_beam()
  IF(input_string(1:4) == "CELL") return


 !---- OBTAIN SPACE GROUP (LATT / SYMM) ----!
 ! >> LATT
 call Read_Shx_Latt(fileshx, npos, nb_lines, n_latt)

 ! >> SYMM
 call Read_Shx_Symm(fileshx, npos, nb_lines, nb_symm_op, car_symop)
 symm_nb = nb_symm_op

 if (n_latt > 0) then
  nb_symm_op=nb_symm_op+1
  car_symop(nb_symm_op)='-X,-Y,-Z'
  centro_string = 'centric'
 else
  centro_string = 'acentric'
 end if

 keyword_SYMM = .true.

 select case (abs(n_latt))
    case (2) ! I
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y+1/2,Z+1/2'
       nlatt_string = 'I centred)'
    case (3) ! Rom, Hex
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+2/3,Y+1/3,Z+1/3'
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/3,Y+2/3,Z+2/3'
       nlatt_string = 'Rhomb. Hex)'
    case (4) ! F
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X,Y+1/2,Z+1/2'
       nlatt_string = 'F centred)'
    case (5) ! A
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X,Y+1/2,Z+1/2'
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y,Z+1/2'
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y+1/2,Z'
       nlatt_string = 'A centred)'
    case (6) ! B
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y,Z+1/2'
       nlatt_string = 'B centred)'
    case (7) ! C
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y+1/2,Z'
       nlatt_string = 'C centred)'
 end select ! nl
 symm_nb = nb_symm_op

 if(input_string(1:6) /= 'NO_OUT') then
  call write_info(' ')
  IF(ABS(n_latt) /= 1) then
   write(message_text,'(a,i2,4a)') '  . LATT: ',n_latt, ' (',TRIM(centro_string), ' / ',trim(nlatt_string)
  else
   write(message_text,'(a,i2,3a)') '  . LATT: ',n_latt, ' (',TRIM(centro_string), ' )'
  endif
  call write_info(TRIM(message_text))
 endif 


 do i=1, nb_symm_op
   op_string = car_symop(i)
   !call remove_car(op_string, ' ')   ! remove blank in string
   op_string = remove_car(op_string, ' ')
   !call replace_car(op_string, ',', ' ')
   op_string = replace_car(op_string,',', ' ') 
   i1=INDEX(op_string, ' ')
   symm_op_string(1,i) = op_string(1:i1-1)
   i2=INDEX(TRIM(op_string), ' ' , back=.TRUE.)
   symm_op_string(2,i) = op_string(i1+1:i2-1)
   symm_op_string(3,i) = op_string(i2+1:)
 end do

 if(input_string(1:6) /= 'NO_OUT') then
  if(nb_symm_op /=0) then
   call write_info(' ')
   write(message_text,*) ' . Number of symmetry operations: ',nb_symm_op
   call write_info(TRIM(message_text))
   do i=1, nb_symm_op
    write(message_text,'(5x,a,i3,4a)') '#',i ,': ', (symm_op_string(n, i), n=1,3)
    call write_info(TRIM(message_text))
   end do
  end if
 endif
 
 ! >>> DEDUCE SPACE GROUP
  car_gsp=' '
  call set_spacegroup(car_gsp, SPG, car_symop, nb_symm_op, 'gen')
  call get_hallsymb_from_gener(SPG)
  space_group_symbol = SPG%Spg_Symb
  !space_group_multip = SPG%multip
  if(input_string(1:6) /= 'NO_OUT') then 
   IF(LEN_TRIM(space_group_symbol) /=0) then
    keyword_SPGR = .true.
    call write_info(' ')
    call write_info('     >> SPACE GROUP : '//TRIM(space_group_symbol))
   else
    keyword_SPGR = .false.
    call write_info(' ')
    call write_info('     >> Unknown space group !!')
   endif
  endif 



 !---- ATOMS ----!
  call Read_Shx_Cont(fileshx, npos, nb_lines, n_elem_atm, elem_atm, n_elem)
  keyword_SFAC_UNIT = .true.
  call Read_Shx_Fvar(fileshx, npos, nb_lines, n_fvar, fvar)

  if(input_string(1:6) /= 'NO_OUT') then
   if(n_elem_atm < 10) then
    write(fmt_sfac, '(a,i1,a)') '(a,',n_elem_atm,'(2x,a))'
    write(fmt_unit, '(a,i1,a)') '(a,',n_elem_atm,'F6.1)'
   elseif(n_elem_atm < 100) then
    write(fmt_sfac, '(a,i2,a)') '(a,',n_elem_atm,'(2x,a))'
    write(fmt_unit, '(a,i2,a)') '(a,',n_elem_atm,'F6.1)'
   endif

   call write_info('')
   write(message_text,fmt_sfac) '  . SFAC: ', elem_atm(1:n_elem_atm)
   call write_info(TRIM(message_text))
   write(message_text,fmt_unit) '  . UNIT: ', n_elem(1:n_elem_atm)
   call write_info(TRIM(message_text))

   call write_info('')
   write(message_text,*) ' . FVAR: ',fvar(1:n_fvar)
   call write_info(TRIM(message_text))
  end if
  
  IF(n_fvar > 1) npos = npos + 1
  call Read_Shx_Atom(fileshx, npos, nb_lines, n_fvar, fvar, elem_atm, crystal_cell, Atm_list)
  !! pb si pas d'atome dans le fichier.ins: cas d'un fichier pour SHELXS !!

   nb_atom = Atm_list%natoms
   call write_info('')
   do i=1,nb_atom

     atom_label(i)     = Atm_list%atom(i)%lab
     atom_type (i)     = Atm_list%atom(i)%ChemSymb
     atom_coord(1:3,i) = Atm_list%atom(i)%x
     atom_occ(i)       = Atm_list%atom(i)%occ
     atom_mult(i)      = Atm_list%atom(i)%mult
     atom_mult(i)      = Get_Multip_Pos(atom_coord(1:3,i), SPG)

     !atom_Ueq(i)       = Atm_list%atom(i)%u(7)
     !atom_Biso(i)      = Atm_list%atom(i)%u(7)*8.0*pi*pi

     atom_Ueq(i)       = Atm_list%atom(i)%Ueq
     atom_Biso(i)      = Atm_list%atom(i)%Ueq*8.0*pi*pi
     if (atom_Biso(i) < 0.0) atom_Biso(i)=1.0

     if(input_string(1:6) /= 'NO_OUT') then
      write(message_text,'(2x,i3, 2(1x,a4),5F10.5,I4)') i,  atom_label(i), atom_type (i) , atom_coord(1:3,i), atom_occ(i) ,  atom_Ueq(i), atom_mult(i)
      call write_info(TRIM(message_text))
     endif 
   end do

   if(input_string(1:6) /= 'NO_OUT')  call Deallocate_atom_list(Atm_list)
   if (allocated(fileshx)) deallocate(fileshx)
   
   if(input_string(1:6) == 'NO_OUT') then  ! CRYSCAL create_ins
   ! modif. des labels pour éviter les doublons
    do i=1, Atm_list%natoms
     n_sl(i) = 1
     k       = 0
     do j=i+1, Atm_list%natoms
      if(Atm_list%atom(j)%lab == Atm_list%atom(i)%lab) then
       n_sl(i) = n_sl(i) + 1
       k = k + 1
       long = len_trim(Atm_list%atom(j)%lab)
       write(tmp_string, '(a,i1)') Atm_list%atom(j)%lab(1:long),k
       Atm_list%atom(j)%lab = tmp_string       
       !if(k==1) then
       ! long = len_trim(Atm_list%atom(i)%lab)
       ! write(tmp_string, '(a,a1)') Atm_list%atom(i)%lab(1:long),"0" 
       ! Atm_list%atom(i)%lab = tmp_string 
       !endif 
       !write(*,*) Atm_list%atom(j)%lab ,'  ',  Atm_list%atom(i)%lab

      end if
     end do 
    end do
   endif

  ! pour compatibilite avec CRYSCAL
   nb_atoms_type                = n_elem_atm
   SFAC_type(1:nb_atoms_type)   = elem_atm(1:n_elem_atm)
   SFAC_number(1:nb_atoms_type) = n_elem(1:nb_atoms_type)
   !atom_label(1:nb_atom)        = atom_type(1:nb_atom)
   Z_unit                       = REAL(Z_unit_INS)

! ------------------------------

 do        ! lecture du fichier d'entree
  READ(UNIT=1, '(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) EXIT   ! fin du fichier
  read_line = ADJUSTL(read_line)
  read_line = u_case(TRIM(read_line))

  IF (LEN_TRIM(read_line) == 0) cycle
  IF (read_line(1:1) == '! ' .or. read_line(1:1) == '#') cycle

  if(read_line(1:4) == 'SIZE') then
   read(read_line(5:), *, IOSTAT=i_error) (var(i), i=1,3)
   IF(i_error /=0) then
    call error_message('SIZE')
     return
    endif

   crystal%size(1:3) = var(1:3)
   keyword_SIZE = .true.
  end if
 END do

 return
end subroutine read_INS_input_file

!--------------------------------------------------------------------------------------------------
subroutine read_INS_SHELX_lines() 
 USE SHELX_module,   only : SHELX_line_nb, SHELX_line
 use cryscal_module, only : INS_read_unit
 use macros_module,  only : u_case
  
 implicit none
  character (len=256)                      :: read_line
  integer                                  :: i_error
  
  read_line = ''  
  SHELX_line(:) = ''
  SHELX_line_nb = 0
  
  rewind(unit = INS_read_unit)
  do
   read(unit = INS_read_unit , '(a)', iostat=i_error) read_line
   if(i_error /=0) exit
   read_line = adjustl(read_line)
   read_line = u_case(read_line)
   
   select case(read_line(1:4))
    case ('TITL', 'CELL', 'LATT', 'SYMM', 'SFAC', 'UNIT')
     cycle
    
    case ('FVAR') 
     exit
     
    case default
     SHELX_line_nb = SHELX_line_nb + 1
     read(read_line, '(a)', iostat=i_error) SHELX_line(SHELX_line_nb)
      
   end select  
   
   
  end do
  
end subroutine read_INS_SHELX_lines
!----------------------------------------------------------------
subroutine create_TRANSF_ins
 use cryscal_module, only : INS_unit, unit_cell, SPG, wavelength, Z_unit, Mat, nb_atom, &
                            ATOM_type, Atom_Ueq, Atom_occ, ATOM_label, new_atom_coord
 use SHELX_module 
 use macros_module,  only : u_case
 use IO_module

 implicit none
  integer           :: i, atom_order
 
 open (unit = INS_unit, file = 'CRYSCAL_transf.ins')
  !WRITE(UNIT = INS_unit, '(3(a, 3(1x,F7.4)))')  'TITL  after transformation with matrix: ', Mat(1,:), '  ', Mat(2,:), '  ', Mat(3,:)
  WRITE(UNIT = INS_unit, '(2a)')            'TITL new space group: ', TRIM(SPG%SPG_Symb)
  WRITE(UNIT = INS_unit, '(a,F8.5,6F9.4)')  'CELL ', wavelength, unit_cell%new_param(1:6)
  WRITE(UNIT = INS_unit, '(a,F8.2,6F9.4)')  'ZERR ', Z_unit,     unit_cell%new_param_ESD(1:6)
  
  call get_SHELX_nlatt
  WRITE(UNIT = INS_unit, '(a,I4)')          'LATT ', n_latt
  
 !---- Symm ----!
  do i=2,SPG%numops
   write(UNIT = INS_unit, fmt="(a)") "SYMM "//u_case(SPG%symopsymb(i))
  end do

 !---- SFAC, UNIT
  if(n_elem_atm < 10) then
   write(fmt_sfac, '(a,i1,a)') '(a,',n_elem_atm,'(2x,a))'
   write(fmt_unit, '(a,i1,a)') '(a,',n_elem_atm,'F6.1)'
  elseif(n_elem_atm < 100) then
   write(fmt_sfac, '(a,i2,a)') '(a,',n_elem_atm,'(2x,a))'
   write(fmt_unit, '(a,i2,a)') '(a,',n_elem_atm,'F6.1)'
  endif
  WRITE(fmt_fvar, '(a,i1,a)') '(a,',n_fvar,'F8.4)'
  write(unit = INS_unit, fmt_sfac) 'SFAC ', elem_atm(1:n_elem_atm)
  write(unit = INS_unit, fmt_unit) 'UNIT ', n_elem(1:n_elem_atm)
  !--- non interpreted lines -----------
  if(SHELX_line_nb /=0) then
   do i=1, SHELX_line_nb
    write(unit = INS_unit, '(a)') trim(SHELX_line(i))
   end do
  endif 

  write(unit = INS_unit, fmt_fvar) 'FVAR ', fvar(1:n_fvar)


 !--- Atoms -----!
  do i=1,nb_atom
    call get_atom_order(atom_type(i), atom_order)
    write(UNIT=ins_unit,'(a,1x,I4,5F10.5)')  atom_label(i), atom_order , new_atom_coord(1:3,i), 10.+atom_occ(i) ,  atom_Ueq(i)
  end do

  WRITE(UNIT = ins_unit, '(a)') 'HKLF  4'
  WRITE(UNIT = ins_unit, '(a)') 'END'

  call write_info('')
  call write_info('   >> The CRYSCAL_transf.ins file has been created.')
  call write_info('')
 
 close (unit = INS_unit)
end subroutine create_TRANSF_ins 

!----------------------------------------------------------------

subroutine create_NEW_ins
 USE cryscal_module
 USE macros_module,                    ONLY : remove_car, replace_car, l_case
 use CFML_Crystallographic_Symmetry,        only : set_spacegroup,   get_hallsymb_from_gener
 USE SHELX_module
 USE IO_module

 implicit none
  INTEGER                          :: i, i1, i2
  INTEGER                          :: op_num, atom_order
  CHARACTER (LEN=12), DIMENSION(3) :: op_STRING
  REAL, DIMENSION(3,3)             :: op_ROT
  REAL, DIMENSION(3)               :: op_TRANS
  REAL, DIMENSION(3,3)             :: new_op_ROT
  REAL, DIMENSION(3)               :: new_op_TRANS



 OPEN(UNIT = INS_unit, FILE='CRYSCAL_new.ins')
  WRITE(UNIT = INS_unit, '(3(a, 3(1x,F7.4)))')   'TITL  after transformation with matrix: ', Mat(1,:), '  ', Mat(2,:), '  ', Mat(3,:)
  WRITE(UNIT = INS_unit, '(a,F8.5,6F9.4)')  'CELL ', wavelength, unit_cell%new_param(1:6)
  WRITE(UNIT = INS_unit, '(a,F8.2,6F9.4)')  'ZERR ', Z_unit,     unit_cell%new_param_ESD(1:6)
  WRITE(UNIT = INS_unit, '(a,I4)')          'LATT ', n_latt

 do op_num = 1, symm_nb
   input_line = TRIM(car_symop(op_num))
   !call remove_car(input_line, ' ')
   input_line = remove_car(input_line, ' ')
   !call replace_car(input_line, ',', ' ')
   input_line = replace_car(input_line, ',', ' ')
   i1=INDEX(input_line, ' ')
   i2=INDEX(TRIM(input_line), ' ', back=.TRUE.)
   op_STRING(1) = input_line(1    : i1-1)
   op_STRING(2) = input_line(i1+1 : i2-1)
   op_STRING(3) = input_line(i2+1 : )
   op_STRING = ADJUSTL(op_STRING)
   call get_op_ROT_TRANS(op_STRING, op_ROT, op_TRANS)

   new_op_ROT   = MATMUL(MAT, op_ROT)
   new_op_TRANS = MATMUL(MAT, op_TRANS)

   call get_op_STRING(new_op_ROT, new_op_TRANS, op_string)
   op_string = ADJUSTL(op_string)
  
  ! ------------------------------------------------
   !call replace_car(op_string(1), "y", "x")
   !call replace_car(op_string(1), "z", "x")
   !call replace_car(op_string(2), "x", "y")
   !call replace_car(op_string(2), "z", "y")
   !call replace_car(op_string(3), "x", "z")
   !call replace_car(op_string(3), "y", "z")

   op_string(1) = replace_car(op_string(1), "y", "x")
   op_string(1) = replace_car(op_string(1), "z", "x")
   op_string(1) = replace_car(op_string(2), "x", "y")
   op_string(1) = replace_car(op_string(2), "z", "y")
   op_string(1) = replace_car(op_string(3), "x", "z")
   op_string(1) = replace_car(op_string(3), "y", "z")
   
  ! -----------------------------------------------------

   if (op_num == 1) then
    IF(ABS(n_latt) /= 1) then
     write(message_text,'(a,i2,4a)') '  LATT: ',n_latt, ' (',TRIM(centro_string), ' / ',trim(nlatt_string)
    else
     write(message_text,'(a,i2,3a)') '  LATT: ',n_latt, ' (',TRIM(centro_string), ' )'
    endif
    call write_info(TRIM(message_text))
   endif

   WRITE(car_symop(op_num), '(5a)')  op_string(1),  ',',  op_string(2), ',', op_string(3)
   !call remove_car(car_symop(op_num), ' ')
   car_symop(op_num) = remove_car(car_symop(op_num), ' ')

   WRITE(message_text, '(a,i2,1x,6a)') '  SYMM# ', op_num, ': ', TRIM(op_string(1)),  ',  ',  TRIM(op_string(2)), ',  ', TRIM(op_string(3))
   call write_info(TRIM(message_text))
   WRITE(UNIT = INS_UNIT, '(2a)') 'SYMM ',TRIM(car_symop(op_num))

 END do
  ! deduction du nouveau symbol du groupe
do op_num = 1, symm_nb
   WRITE(message_text, *) op_num, car_symop(op_num)
   call write_info(trim(message_text))
end do

!  call set_spacegroup(' ', SPG, car_symop, nb_symm_op, 'gen')
  call set_spacegroup(' ', SPG, car_symop, symm_nb, 'gen')
  call get_hallsymb_from_gener(SPG)

  WRITE(message_text, '(2a)') '  > New Space group: ', TRIM(SPG%SPG_Symb)
  call write_info('')
  call write_info(TRIM(message_text))


  if(n_elem_atm < 10) then
   write(fmt_sfac, '(a,i1,a)') '(a,',n_elem_atm,'(2x,a))'
   write(fmt_unit, '(a,i1,a)') '(a,',n_elem_atm,'F6.1)'
  elseif(n_elem_atm < 100) then
   write(fmt_sfac, '(a,i2,a)') '(a,',n_elem_atm,'(2x,a))'
   write(fmt_unit, '(a,i2,a)') '(a,',n_elem_atm,'F6.1)'
  endif
  WRITE(fmt_fvar, '(a,i1,a)') '(a,',n_fvar,'F8.4)'
  write(unit = INS_unit, fmt_sfac) 'SFAC ', elem_atm(1:n_elem_atm)
  write(unit = INS_unit, fmt_unit) 'UNIT ', n_elem(1:n_elem_atm)
  write(unit = INS_unit, fmt_fvar) 'FVAR ', fvar(1:n_fvar)

  do i=1,nb_atom
    !write(UNIT=ins_unit,'(a,1x,a4,5F10.5)')  atom_label(i), atom_type (i) , new_atom_coord(1:3,i), 10.+atom_occ(i) ,  atom_Ueq(i)
    call get_atom_order(atom_type(i), atom_order)
    write(UNIT=ins_unit,'(a,1x,I4,5F10.5)')  atom_label(i), atom_order , new_atom_coord(1:3,i), 10.+atom_occ(i) ,  atom_Ueq(i)
  end do

  WRITE(UNIT = ins_unit, '(a)') 'HKLF  4'
  WRITE(UNIT = ins_unit, '(a)') 'END'


 close (unit = INS_unit)

  call write_info('')
  call write_info('  CRYSCAL_new.INS file has been created ')
  call write_info('')

 RETURN
END subroutine create_NEW_ins

!--------------------------------------------------------------
subroutine get_atom_order(atom_string, atom_order)
 USE cryscal_module, ONLY : nb_atoms_type , SFAC_type
 USE SHELX_module
 implicit none
  CHARACTER (LEN=*), INTENT(IN)         :: atom_string
  INTEGER,           INTENT(OUT)        :: atom_order
  ! local variables
  INTEGER                               :: i, j

   do j=1, nb_atoms_type
    IF(atom_string(1:) == elem_atm(j)(1:) .or. atom_string(1:) == SFAC_type(j)(1:)) then
     atom_order = j
     return
    endif
   end do

 RETURN
end subroutine get_atom_order


