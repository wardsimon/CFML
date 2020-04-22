!     Last change:  TR    1 Feb 2007    5:58 pm



subroutine read_PCR_input_file
 USE macros_module
 use cryscalc_module
 USE IO_module
 use CFML_Crystallographic_Symmetry,        only : Get_Multip_Pos

 implicit none
  CHARACTER (LEN=256)                      :: read_line
  INTEGER                                  :: i, n_atom, n_lines, i_error
  integer                                  :: so_1, so_2, N_t
  real                                     :: f, f1

  if(debug_proc%level_2)  call write_debug_proc_level(2, "read_PCR_input_file")

  !REWIND(UNIT=PCR_read_unit)


! TITL: keyword COMM
  READ(UNIT=PCR_read_unit, fmt='(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) return   ! fin du fichier
  read_line = ADJUSTL(read_line)
  read_line = u_case(TRIM(read_line))

  ! TITL
  if (read_line(1:4) == 'COMM') then
   READ(read_line(5:), '(a)')  main_title
   if(write_details) then
   call write_info(' ' )
   call write_info('  . TITL: '//trim(main_title))
   end if
  else
   main_title = ''
  end if

 !"! lambda1"
 do
  READ(UNIT=PCR_read_unit, fmt='(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) EXIT   ! fin du fichier
  read_line = l_case(read_line)

  i = INDEX(read_line,'! lambda1')
  if (i/=0) then
   READ(UNIT=PCR_read_unit, fmt=*, IOSTAT=i_error) wavelength
   keyword_WAVE = .true.

   if(write_details) then
   write(message_text,fmt='(a,F10.5)')  '  > WAVE: ', wavelength
   call write_info(TRIM(message_text))
   end if
   exit
  end if
 end do



 !"!Nat"
 do
  READ(UNIT=PCR_read_unit, fmt='(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) EXIT   ! fin du fichier

  i = INDEX(read_line,'!Nat')
  if (i/=0) then
   READ(UNIT=PCR_read_unit, fmt=*, IOSTAT=i_error) nb_atom
   exit
  end if
 END do



 !"!Atom Typ"
 do
  READ(UNIT=PCR_read_unit, fmt='(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) EXIT   ! fin du fichier
  read_line = l_case(read_line)

  i = INDEX(read_line,'!atom')

  if (i/=0) then
   do
    READ(UNIT=PCR_read_unit, fmt='(a)', IOSTAT=i_error) read_line
    IF(i_error /=0) cycle
    read_line=ADJUSTL(read_line)
    IF(read_line(1:1) == '!') cycle
    BACKSPACE(UNIT=PCR_read_unit)
    exit
   end do

   do n_atom=1, nb_atom
    READ(UNIT=PCR_read_unit, fmt=*, IOSTAT=i_error) atom_label(n_atom),atom_typ(n_atom), (atom_coord(i, n_atom) ,i=1,3), &
                                             atom_Biso(n_atom), atom_Occ(n_atom), so_1, so_2, N_t
    IF(N_t == 2) then
     n_lines = 3
    else
     n_lines = 1
    endif

    do i=1, n_lines
     READ(UNIT=PCR_read_unit, fmt='(a)', IOSTAT=i_error) read_line
     IF(i_error < 0) EXIT   ! fin du fichier
    end do
    if(write_details) then
    WRITE(message_text,fmt='(a,2a6,5(1x,F8.5))') '  > ATOM: ', trim(atom_label(n_atom)),trim(atom_typ(n_atom)),  &
                                                 (atom_coord(i,n_atom),i=1,3), atom_Biso(n_atom), atom_occ(n_atom)
    call write_info(TRIM(message_text))
    end if

   end do

   exit
  end if
 END do

 REWIND (UNIT=PCR_read_unit)
 ! "<--Space group symbol"
 do
  READ(UNIT=PCR_read_unit, fmt='(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) EXIT   ! fin du fichier

  i = INDEX(read_line,'<--Space group symbol')
  if (i/=0) then
   READ(read_line(1:i-1), fmt='(a)', IOSTAT=i_error) space_group_symbol
   IF(i_error /= 0) EXIT
   space_group_symbol = ADJUSTL(space_group_symbol)
   if (write_details) then
   WRITE(message_text, fmt='(2a)') '  > Space group : ', space_group_symbol
   call write_info(TRIM(message_text))
   end if
   keyword_SPGR = .true.
   call space_group_info()
   exit
  end if
 END do

 if(keyword_SPGR) then
  atom_mult(1) = Get_Multip_Pos(atom_coord(1:3,1), SPG)
  f1 = atom_occ(1) * REAL(SPG%Multip) / REAL(atom_mult(1))
  f1 = 1./f1
  do i = 1, nb_atom
   atom_mult(i) = Get_Multip_Pos(atom_coord(1:3,i), SPG)
   f =  SPG%Multip / atom_mult(i)
   atom_occ_perc(i) = atom_occ(i) *  f * f1
   end do
 endif


 !"!     a          b         c        alpha      beta       gamma"
 do
  READ(UNIT=PCR_read_unit, fmt='(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) EXIT   ! fin du fichier
  read_line = l_case(read_line)
  i = INDEX(read_line,'!     a          b         c        alpha      beta       gamma')

  if (i/=0) then
   READ(UNIT=PCR_read_unit, fmt=*, IOSTAT=i_error) (unit_cell%param(i),i=1,6)
   !IF(i_error /= 0) EXIT
   if(write_details) then
   WRITE(message_text,fmt='( a,6F10.4)') '  > CELL PARAMETERS: ', (unit_cell%param(i),i=1,6)
   call write_info(TRIM(message_text))
   end if
   keyword_CELL = .true.
   exit
  endif

 END do

 return
end subroutine read_PCR_input_file

!----------------------------------------------------------------

subroutine read_CEL_input_file(input_file)
 USE cryscalc_module, only         : CEL_read_unit, unit_cell, SPG, space_group_symbol, nb_atom, &
                                     atom_label, atom_typ, atom_coord, atom_occ_perc, atom_occ, atom_Biso, &
                                     keyword_CELL, keyword_SPGR, debug_proc
 USE macros_module,   only         : u_case, nombre_de_colonnes
 USE IO_module
 !use CFML_crystallographic_symmetry, ONLY : set_spacegroup

 implicit none
  character (len=*), intent(in)    :: input_file
  character (len=256)              :: read_line, adjusted_line
  integer                          :: i_error, nb_col
  integer                          :: i, n_atom
  integer                          :: atomic_number

  if(debug_proc%level_2)  call write_debug_proc_level(2, "read_CEL_input_file")


  ! line 1 : CELL
  read(unit = CEL_read_unit, fmt='(a)', iostat=i_error) read_line
  if(i_error /=0) then
   call write_info('')
   call write_info(' Error reading '//trim(input_file)// ' at line 1 (CELL).')
   call write_info('')
   close(unit=CEL_read_unit)
   return
  end if
  adjusted_line = adjustl(read_line)
  if(u_case(adjusted_line(1:4)) == 'CELL') then
   read(adjusted_line(6:), fmt=*) (unit_cell%param(i),i=1,6)
  else
   call write_info('')
   call write_info(' Error reading ' //trim(input_file)// ' at line 1 (CELL).')
   call write_info('')
   close(unit=CEL_read_unit)
   return
  endif
  keyword_CELL = .true.

  ! line 2 : natoms
   read(unit = CEL_read_unit, fmt='(a)', iostat=i_error) read_line
  if(i_error /=0) then!
   call write_info('')
   call write_info(' Error reading ' //trim(input_file) // ' at line 2 (natom).')
   call write_info('')
   close(CEL_read_unit)
   return
  end if
  adjusted_line = adjustl(read_line)
  if(u_case(adjusted_line(1:5)) == 'NATOM') then
   read(adjusted_line(7:), fmt=*) nb_atom
  else
   nb_atom = 1
   backspace(unit = CEL_read_unit)
  endif
  do n_atom = 1, nb_atom
   read(unit = CEL_read_unit, fmt='(a)', iostat = i_error) read_line
   if(i_error /=0) exit
   adjusted_line = adjustl(read_line)
   call nombre_de_colonnes(adjusted_line, nb_col)
   if(nb_col < 5) then
    call write_info('')
    call write_info('  > Problem reading '//trim(input_file)//' file. Wrong format ?')
    call write_info('')
    close(unit = CEL_read_unit)
    return
   elseif(nb_col == 7) then
    read(adjusted_line, fmt=*) atom_label(n_atom), atomic_number, (atom_coord(i, n_atom) ,i=1,3), &
                              atom_occ_perc(n_atom), atom_Biso(n_atom)
   elseif(nb_col == 6) then
    read(adjusted_line, fmt=*) atom_label(n_atom), atomic_number, (atom_coord(i, n_atom) ,i=1,3), &
                               atom_occ_perc(n_atom)
   else
    read(adjusted_line, fmt=*) atom_label(n_atom), atomic_number, (atom_coord(i, n_atom) ,i=1,3)
    atom_occ_perc(n_atom) = 1.
   endif

   call get_atom_typ(atomic_number, atom_typ(n_atom))
  end do



  ! ligne suivante : groupe d'espace
  read(unit = CEL_read_unit, fmt='(a)', iostat=i_error) read_line
  if(i_error /=0) then!
   call write_info('')
   call write_info(' Error reading ' //trim(input_file) // ' (rgnr).')
   call write_info('')
   close(unit = CEL_read_unit)
   return
  end if
  adjusted_line = adjustl(read_line)
  if(u_case(adjusted_line(1:4)) == 'RGNR') then
   read(adjusted_line(6:), '(a)') space_group_symbol
   call space_group_info()
   if(SPG%NumSpg /=0) keyword_SPGR = .true.
  end if

  close(CEL_read_unit)



 return
end subroutine read_CEL_input_file

!----------------------------------------------------------------------------------------------

subroutine Get_atom_typ(atomic_number, atomic_symbol)
 USE cryscalc_module,                  ONLY : known_atomic_label
 USE CFML_Scattering_Chemical_Tables, ONLY : set_chem_info, chem_info, Num_Chem_Info

 implicit none
  integer, intent(in)                    :: atomic_number
  character (len=2), intent(out)         :: atomic_symbol
  integer                                :: i

  CALL set_chem_info

  do i=1 , Num_Chem_Info
   if(atomic_number == chem_info(i)%Z) then
    atomic_symbol = chem_info(i)%symb
    exit
   end if
  end do


 return
end subroutine Get_atom_typ
