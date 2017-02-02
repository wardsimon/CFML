subroutine create_INS_from_SOLVE

! . lecture d'un fichier de sortie de SIR(92,97,2002,2004) ou SHELXS
! . lecture du fichiers struct.cif
! . creation d'un fichier.INS pour SHELXL:
!     . wavelength
!     . unit cell and esd's
!     . atoms renumbering
!     . cosmectic keywords: ACTA, BOND$H ...

 USE IO_module
 USE macros_module,   only : test_file_exist, u_case, Get_current_folder_name, Get_sample_ID_from_folder_name, nombre_de_colonnes
 use cryscalc_module, only : tmp_unit, tmp_2_unit, message_text, Atm_list, unit_cell, Z_unit_ins, Create_INS,   &
                             get_sample_ID, sample_job, input_line, debug_proc
 implicit none
  logical                            :: file_exist
  character (len=256)                :: struct_CIF, job_RES, job_INS, job_HKL
  integer                            :: i, i1, i_error, nb_res, selected_res
  integer                            :: nb_arg, long
  integer                            :: Atom_order
  logical                            :: skip
  real                               :: skip_level
  integer                            :: skip_at_nb
  character (len=256)                :: read_line
  character (len=256)                :: folder_name
  character (len=256), dimension(26) :: res_file   ! nom des fichiers .RES dans le repertoire
  character (len=256), dimension(26) :: hkl_file   ! nom des fichiers .HKL dans le repertoire
  character (len=256), dimension(26) :: solution_spg        ! space group
  character (len=256), dimension(26) :: solution_R1         ! R1
  real                               :: temperature
  logical                            :: shelxs_file, shelx_solve, H_present
  character (len=256)                :: SHELX_title, SHELX_line
  character (len=4), dimension(20)   :: SFAC_type

  if(debug_proc%level_2)  call write_debug_proc_level(2, "create_INS_file_from_solve")

 shelxs_file = .false.
 shelx_solve = .false.
 H_present   = .false.
 skip_at_nb  = 0


 skip_level = Create_INS%U_threshold
 if(.not. get_sample_ID) then
  job_ins = 'job.ins'
 else
  call Get_current_folder_name(folder_name)
  call get_sample_ID_from_folder_name(folder_name, sample_job)
  if(sample_job(1:3) == 'job') then
   job_INS = 'job.ins'
  else
   job_INS = trim(sample_job)//'.ins'
  endif
 end if

 call write_info(" ")
 call write_info("*****************************************************************")
 call write_info(" ")
 call write_info("    . read STRUCT.CIF file and get cell parameters with esd's    ")
 call write_info("    . read import.RES file created by SIRxx or SHELXS            ")
 call write_info("")
 call write_info("    >> create .INS file for SHELXL with correct esd's            ")
 call write_info("       and different useful SHELXL keywords (ACTA, BOND$H...)    ")
 call write_info(" ")
 call write_info("*****************************************************************")
 call write_info(" ")

 struct_CIF = 'struct.cif'
 call test_file_exist("struct.cif",  file_exist, 'out')
 if(.not. file_exist) return
 call write_info('')
 call write_info('  . CIF file : '//trim(struct_cif))

 job_res = 'import.res'
 inquire(file = trim(job_res), exist = file_exist)
 if(.not. file_exist) then
  call system('dir import*.res /B /ON > res_files.lst')
  open(unit = tmp_unit, file = 'res_files.lst')
  nb_res = 0
  solution_spg = "?"
  solution_R1  = "?"
  do
   read(unit = tmp_unit, fmt='(a)', iostat = i_error) read_line
   if(i_error /=0) exit
   nb_res = nb_res + 1
   read(read_line, '(a)') res_file(nb_res)
   i1 = index(res_file(nb_res), ".")
   hkl_file(nb_res) = res_file(nb_res)(1:i1)//'hkl'
    if(nb_res == 26) exit
   ! lecture du fichier .RES et recup. du groupe d'espace
   open(unit=tmp_2_unit, file=trim(res_file(nb_res)))
    do
     read(unit=tmp_2_unit, fmt='(a)', iostat=i_error) read_line
     if(i_error /= 0) exit
     if(read_line(1:22) == 'REM SHELXT solution in') write(solution_spg(nb_res), '(a)') read_line(24:)
     if(read_line(1:6)  == 'REM R1') then
      i1 = index (read_line, ',')
      if(i1 /=0) solution_R1(nb_res) = read_line(7:i1-1)
      exit
     end if
    end do
   close(unit=tmp_2_unit)
  end do
  close(unit = tmp_unit)

  call system('del res_files.lst')
  if(nb_res > 1) then
   do
    call write_info(' ')
    do i = 1, nb_res
     WRITE(message_text,'(5x,i2,7a)')  i, '. ', TRIM(res_file(i)), ' [',trim(solution_spg(i)),' / R1 = ',trim(solution_R1(i)), ']'
     call write_info(trim(message_text))
    end do
    call write_info(' ')
    call write_info('  > Select the .RES file to extract atomic positions: ')
    call read_input_line(input_line)
    read(input_line, *) selected_res
    !READ(*,*) selected_res
    IF(selected_res >= 1 .and. selected_res <= nb_res) exit
   END do
   job_RES = res_file(selected_res)
   job_HKL = hkl_file(selected_res)
  else
   job_RES = res_file(1)
   job_HKL = hkl_file(1)
  endif
 else
  job_HKL = 'import.hkl'
 end if  ! fin de la condition if(import.res not exist)
 call write_info('  . RES file : '//trim(job_res))
 call write_info('  . HKL file : '//trim(job_hkl))

 if(Create_INS%temperature < -900.) then
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
  temperature = Create_INS%temperature
 endif
 Create_INS%temperature = temperature

 call write_info('')
 write(message_text, '(a,F8.2)') '  . Temperature (C): ', temperature
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
  Create_INS%U_threshold = skip_level
  write(message_text, '(a,F5.2)') '  . U_threshold: ', Create_INS%U_threshold
 endif

 ! creation du fichier job.hkl
 call system("copy "//trim(job_hkl)//" job.hkl")


 ! lecture fichier.INS pour recuperer Z_unit et atom_list
 ! create_ins: call read_INS_file_CFML(job_res)
 call read_INS_input_file(trim(job_res), 'NO_OUT')

 ! lecture fichier STRUCT.CIF pour recuperer cell et cell_esd
 open(unit = tmp_unit, file = trim(struct_CIF))
 call read_CIF_input_file(trim(struct_CIF), 'NO_OUT')

 ! lecture fichier.RES
 open(unit = tmp_unit+1, file=trim(job_RES))

 open(unit = tmp_unit+2, file=trim(job_INS))
  do
   read(unit = tmp_unit+1, fmt='(a)', iostat = i_error) read_line
   read_line = adjustl(read_line)
   if(len_trim(read_line) == 0) shelxs_file = .true.

   select case (read_line(1:4))
     case ('TITL')
      read(read_line(5:), '(a)') SHELX_title
      SHELX_title = adjustl(SHELX_title)
      if(SHELX_title(1:12) == 'import_a.res') shelx_solve = .true.
      write(unit=tmp_unit+2, fmt='(a)') trim(read_line)

     case ('CELL')
      if(shelx_solve) then
       write(unit=tmp_unit+2, fmt='(a)') trim(read_line)
      else
       write(unit=tmp_unit+2, fmt='(a,6F10.4)') 'CELL    0.71073', (unit_cell%param(i), i=1,6)
      end if

     case ('ZERR')
      if(shelx_solve) then
       write(unit=tmp_unit+2, fmt='(a)') trim(read_line)
      else
       write(unit=tmp_unit+2, fmt='(a6,I8,6F10.4)') 'ZERR  ',Z_unit_INS, (unit_cell%param_ESD(i),i=1,6)
      end if

     case ('OMIT')   ! l'instruction OMIT n'existe pas dans le fichier cree par SIR2004

     case ('SFAC')
      read(read_line(5:), '(a)') SHELX_line
      call nombre_de_colonnes(SHELX_line, nb_arg)
      read(SHELX_line, *)  SFAC_type(1:nb_arg)
      do i=1, nb_arg
       long = len_trim(SFAC_type(i))
       if(long == 1) then
        if (SFAC_type(i)(1:1) == "H") then
         H_present = .true.
         exit
        end if
       end if
      end do
      write(unit=tmp_unit+2, fmt='(a)') trim(read_line)

     case ('LIST')
       write(unit=tmp_unit+2, fmt='(a)') 'LIST 4'

     case ('L.S.')
       write(unit=tmp_unit+2, fmt='(a)') 'L.S. 4'

     case ('BOND')
      !read(read_line(5:), '(a)') SHELX_line
      !read(SHELX_line, *)  SFAC_type(1)
      !long = len_trim(SFAC_type(1))
      !write(*,*) 'SFAC_type(1)(1:2): ', SFAC_type(1)(1:2)
      !write(*,*) 'H_present: ', H_present
      !pause
      !if(long == 2) then
      ! if(SFAC_type(1)(1:2) == "$H" .and. H_present) write(unit=tmp_unit+2, fmt='(a)') 'BOND $H'
      !end if


     case ('PLAN')
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  PLAN -20  : Fourier peak list'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'PLAN -20'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  ACTA      : create .CIF file'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'ACTA'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  BOND $H   : include H in bond lengths / angles table'

       if(H_present) then
       WRITE(unit=tmp_unit+2, fmt='(a)')         'BOND $H'
       end if
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  HTAB      : analyse all hydrogen bonds'
       if(H_present) then
       WRITE(unit=tmp_unit+2, fmt='(a)')         'HTAB'
       end if
       WRITE(unit=tmp_unit+2, fmt='(a)')         'REM  CONF      : all torsion angles except involving hydrogen'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'CONF'
       WRITE(unit=tmp_unit+2, fmt='(a)')         'WPDB -2     !'
       if(create_INS%ANIS) then
        WRITE(unit=tmp_unit+2, fmt='(a)')         'ANIS        !'
       end if
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
     close(unit=tmp_unit+2)


    if(skip_at_nb /=0) then
     write(message_text, fmt='(a,i3,a)') '   ', skip_at_nb, ' atoms have been skipped (too large Uiso).'
     call write_info(trim(message_text))
    else
     call write_info('  >> No skipped atoms due to large Uiso.')
    endif
    call write_info('')

    call write_info('')
    call write_info('   >>> '//trim(job_ins)//' file for SHELXL has been created.')
    call write_info('')

 return
end subroutine create_INS_from_SOLVE


!------------------------------------------------------------------
!subroutine get_atom_order(atom_string, atom_order)
!
! use SHELX_module
! use cryscalc_module, only   : Atm_list
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
