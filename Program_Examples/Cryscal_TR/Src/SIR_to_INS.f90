subroutine create_INS_from_SOLVE

! . lecture d'un fichier de sortie de SIR(92,97,2002,2004) ou SHELXS
! . lecture du fichiers struct.cif
! . creation d'un fichier.INS pour SHELXL:
!     . wavelength
!     . unit cell and esd's
!     . atoms renumbering
!     . cosmectic keywords: ACTA, BOND$H ...

 USE IO_module
 USE macros_module,  only : test_file_exist, u_case
 use cryscal_module, only : tmp_unit, message_text, Atm_list, unit_cell, Z_unit_ins, Create_INS_temperature, Create_INS_U_threshold
 implicit none
  logical                            :: file_exist
  character (len=256)                :: struct_CIF, job_RES, job_INS
  integer                            :: i, i_error, nb_res, selected_res
  integer                            :: long
  integer                            :: Atom_order
  logical                            :: skip
  real                               :: skip_level
  integer                            :: skip_at_nb
  character (len=256)                :: read_line
  character (len=256), dimension(10) :: res_file   ! nom des fichiers .RES dans le repertoire
  real                               :: temperature
  logical                            :: shelxs_file

 shelxs_file = .false.
 
 skip_level = Create_INS_U_threshold
 
 call write_info(" ")
 call write_info("*****************************************************************")
 call write_info(" ")
 call write_info("    . read STRUCT.CIF file and get cell parameters with esd's    ")
 call write_info("    . read import.RES file created by SIRxx or SHELXS            ")
 call write_info("")
 call write_info("    >> create job.INS file for SHELXL with correct esd's         ")
 call write_info("       and different useful SHELXL keywords (ACTA, BOND$H...)    ")
 call write_info(" ")
 call write_info("*****************************************************************")
 call write_info(" ")

 struct_CIF = 'struct.cif'
 call test_file_exist("struct.cif",  file_exist)
 if(.not. file_exist) return
 call write_info('')
 call write_info('  . CIF file : '//trim(struct_cif))
 call write_info('')
 
 job_res = 'import.res'
 inquire(file = trim(job_res), exist = file_exist)
 if(.not. file_exist) then
  call system('dir *.res /B /ON > res_files.lst')
  open(unit = tmp_unit, file = 'res_files.lst')
  nb_res = 0
  do
   read(unit = tmp_unit, fmt='(a)', iostat = i_error) read_line
   if(i_error /=0) exit
   nb_res = nb_res + 1
   read(read_line, '(a)') res_file(nb_res)
   if(nb_res == 10) exit
  end do
  close(unit = tmp_unit)
  
  call system('del res_files.lst')
  if(nb_res > 1) then
   do
    write(*,*) ' '
    do i = 1, nb_res
     WRITE(*,'(5x,i1,2a)')  i, '. ', TRIM(res_file(i))
    end do
    WRITE(*,*) ''
    WRITE(*,*) '  > Select the .RES file to extract atomic positions:'
    READ(*,*) selected_res
    IF(selected_res >= 1 .and. selected_res <= nb_res) exit
   END do
   job_RES = res_file(selected_res)
  else
   job_RES = res_file(1)
  endif
 end if  ! fin de la condition if(import.res not exist)
 call write_info('')
 call write_info('  . RES file : '//trim(job_res))
 call write_info('')

 if(Create_INS_temperature < -900.) then
  write(*,*) 
  write(*,*) ' > Enter temperature (in Celsius) : '
  read(*, '(a)') read_line
  read_line = u_case(read_line)
  read_line = adjustl(read_line)
  if(len_trim(read_line) == 0) then
   temperature = 20.
  elseif(read_line(1:3) == 'AMB' .or. read_line(1:2) == 'RT') then
   temperature = 20.
  elseif(read_line(1:2) == 'BT'  .or. read_line(1:2) == 'LT') then
   temperature = -173.
  else
   long = len_trim(read_line)
   if(read_line(long:long) == 'K') then
    read(read_line(1:long-1), *) temperature
    temperature = temperature - 273.
   else
    read(read_line, *) temperature
   endif
  end if 
  if(abs(temperature) < .1) temperature = 20.
 else
  temperature = Create_INS_temperature 
 endif
 Create_INS_temperature = temperature
 
 call write_info('')
 write(message_text, '(a,F8.2)') '   . Temperature (C): ', temperature 
 call write_info(trim(message_text))
 
 if(skip_level < -900.) then
  skip = .false.
  write(*,*) ' '
  write(*,*) ' > Skip atoms with large Uiso (y/n) ? '
  read(*, '(a)') read_line
  read_line = adjustl(read_line)
  read_line = u_case(read_line)
  if(read_line(1:1) == 'Y') then
   skip = .true.
   write(*,*) ' > Enter skip level (def.=0.1) : '
   read(*,'(a)') read_line
   read_line = adjustl(read_line)
   if(len_trim(read_line) == 0) then
    skip_level = 0.1
   else
    read(read_line, *) skip_level
   endif
  else
   skip = .false.
  endif
 else
  skip  = .true.  
 endif
 
 if(skip) then
  if(skip_level < 0.) skip_level = 0.1
  Create_INS_U_threshold = skip_level
  call write_info('')
  write(message_text, '(a,F5.2)') '   . U_threshold: ', Create_INS_U_threshold
  call write_info(trim(message_text))
 endif 

 
 ! lecture fichier.INS pour recuperer Z_unit et atom_list
 ! create_ins: call read_INS_file_CFML(job_res)
 call read_INS_input_file(trim(job_res), 'NO_OUT')

 
 ! lecture fichier STRUCT.CIF pour recuperer cell et cell_esd
 open(unit = tmp_unit, file = trim(struct_CIF))
 call read_CIF_input_file(trim(struct_CIF), 'NO_OUT')
 
 ! lecture fichier.RES
 open(unit = tmp_unit+1, file=trim(job_RES))
 
 job_INS = 'job.ins'
 open(unit = tmp_unit+2, file=trim(job_INS))
  do 
   read(unit = tmp_unit+1, fmt='(a)', iostat = i_error) read_line
   read_line = adjustl(read_line)
   if(len_trim(read_line) == 0) shelxs_file = .true.
   
   select case (read_line(1:4))
     case ('CELL')
       write(unit=tmp_unit+2, fmt='(a,6F10.4)') 'CELL    0.71073', (unit_cell%param(i), i=1,6)
       
     case ('ZERR')
       write(unit=tmp_unit+2, fmt='(a6,I8,6F10.4)') 'ZERR  ',Z_unit_INS, (unit_cell%param_ESD(i),i=1,6) 
       
     case ('OMIT')   ! l'instruction OMIT n'existe pas dans le fichier cree par SIR2004
     
     case ('LIST')
       write(unit=tmp_unit+2, fmt='(a)') 'LIST 4'
       
     case ('PLAN')
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  PLAN -20  : Fourier peak list'  
       WRITE(unit=tmp_unit+2, fmt='(a)')         'PLAN -20'  
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  ACTA      : create .CIF file'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'ACTA'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  BOND $H   : include H in bond lengths / angles table'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'BOND $H'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  HTAB      : analyse all hydrogen bonds'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'HTAB'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  CONF      : all torsion angles except involving hydrogen'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'CONF'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'WPDB -2     !'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  TEMP      : temperature in celcius'
       WRITE(unit=tmp_unit+2, fmt='(a6,F8.0,a)') 'TEMP  ', temperature, '    ! temperature in celcius'
       if(shelxs_file) then
        write(tmp_unit+2, fmt='(a)')             'FVAR 1.0'  ! le fichier a ete cree par SHELXS et ne contient pas la ligne FVAR
       endif
        
     case ('FVAR', 'MOLE')
       !if(read_line(1:4) == 'MOLE') then
       ! write(tmp_unit+2, '(a)')             'FVAR 1.0'  ! le fichier a ete cree par SHELXS et ne contient pas la ligne FVAR
       !endif
       WRITE(unit=tmp_unit+2, fmt='(a)') trim(read_line)
       
       skip_at_nb = 0
       do i=1, Atm_list%natoms
        if(skip) then
         if(Atm_list%atom(i)%Ueq > skip_level .or. Atm_list%atom(i)%Ueq < 0.00002) then
          skip_at_nb = skip_at_nb + 1
          cycle
         endif 
        endif
        call get_atom_order(Atm_list%atom(i)%ChemSymb, atom_order)
        write(tmp_unit+2,fmt='(a4,1x,I4,5F10.5)') Atm_list%atom(i)%lab,   atom_order , Atm_list%atom(i)%x, &
                                              10.+Atm_list%atom(i)%occ ,  Atm_list%atom(i)%Ueq
                                              
                                             
       end do
       WRITE(tmp_unit+2, fmt='(a)') ' '
       WRITE(tmp_unit+2, fmt='(a)') 'HKLF    4'
       WRITE(tmp_unit+2, fmt='(a)') 'END'
       exit
 
     case('END') 
       WRITE(unit=tmp_unit+2, fmt='(a)') trim(read_line)       
       exit
         
     case default
       WRITE(unit=tmp_unit+2, fmt='(a)') trim(read_line)
       
      end select
     end do
     close(unit=tmp_unit+2)
     
    
    call write_info('')
    write(message_text, fmt='(3a)') '   ', trim(job_INS), ' has been created.'
    call write_info(trim(message_text))
    call write_info('')
    if(skip_at_nb /=0) then
     write(message_text, fmt='(a,i3,a)') '   ', skip_at_nb, ' atoms have been skipped (too large Uiso).'
     call write_info(trim(message_text))     
    else
     call write_info('  >> No skipped atoms due to large Uiso.')
    endif
    call write_info('') 
    
    call write_info('')
    call write_info('   >>> job.ins file for SHELXL has been created.')
    call write_info('')
     
  
end subroutine create_INS_from_SOLVE


!------------------------------------------------------------------
!subroutine get_atom_order(atom_string, atom_order)
!
! use SHELX_module
! use cryscal_module, only   : Atm_list
! implicit none
!  character (len=*), intent(in)      :: atom_string
!  integer,           intent(out)     :: atom_order
!  ! local variables
!  integer                            :: i, j
!  
!  do i=1, Atm_list%natoms
!   do j = 1, n_elem_atm
!    if(atom_string(1:) == elem_atm(j)(1:)) then
!     atom_order = j
!     return
!    endif
!   end do
!  end do
!  
!  
!  return   
! 
! end subroutine get_atom_order
