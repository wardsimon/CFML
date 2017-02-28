!     Last change:  TR   27 Nov 2006    2:44 pm
!     Last change:  TR   27 Nov 2006    2:44 pm
!------------------------------------------------------------------------

subroutine Write_HELP
 USE cryscalc_module, ONLY : nb_help, nb_help_max, debug_proc
 USE IO_module,       ONLY : write_info

 implicit none

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_help")

  IF(nb_help == nb_help_max) call write_header
  call HELP_on_line

  IF(nb_help == nb_help_max) then
   call write_info(' >>> List of CRYSCALC commands line arguments :')
   call write_info('     ----------------------------------------')
   call write_info('')
   call write_cryscalc_CLA
   call write_info(' ')
   call write_info(' ')
   call write_info('   CRYSCALC manual stored in CRYSCALC_manual.txt file. ')
   call write_info(' ')
  endif
  stop

  return
end subroutine Write_HELP

!------------------------------------------------------------------------
 subroutine Write_cryscalc_HTML
  use cryscalc_module,              ONLY : my_browser, browse_cryscalc_HTML, debug_proc
  use IO_module,                    ONLY : write_info
  USE external_applications_module, ONLY : launch_browser

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_cryscalc_HTML")

  !browse_cryscalc_HTML = .false.
  call create_CRYSCALC_HTML
  call write_info(' ')
  call write_info(' ')

  stop

  return
 end subroutine Write_cryscalc_html
!------------------------------------------------------------------------

subroutine Write_version
 use cryscalc_module,      ONLY : cryscalc, CFML_version, Architecture
 use IO_module,            ONLY : write_info

 call Compiler_type
 call Cryscalc_CFML_Versions
 call write_info('')
 call write_info('   > CRYSCALC version      : '//trim(cryscalc%version))
 call write_info('   > Author                : '//trim(cryscalc%author))
 if(cryscalc%date(1:1) /= "?") then
 call write_info('   > Compilation date      : '//trim(cryscalc%date))
 end if
 call write_info('   > Fortran compiler      : '//trim(cryscalc%compiler))
 if(CRYSCALC%option(1:1) /= '?') then
 call write_info('   > Option                : '//trim(cryscalc%option))
 end if
 call write_info('   > CFML version          : '//trim(CFML_version))
 call write_info('   > Architecture          : '//trim(Architecture))

 call write_info('')


 return
end subroutine Write_version

!------------------------------------------------------------------------
subroutine Write_KEYWORD()
  USE IO_module, only : write_info

  implicit none

  call KEYS_on_line()

  call write_info(' ')
  call write_info(' ')
  call write_info('   CRYSCALC keys stored in CRYSCALC_KEYS.txt file. ')
  call write_info(' ')
  stop

  return
end subroutine Write_KEYWORD

!------------------------------------------------------------------------
subroutine write_wave_features()
 USE cryscalc_module, ONLY: ON_SCREEN, write_details, wavelength, keyword_beam, message_text, wavelength, pi, &
                            neutrons, DEVICE, debug_proc
 USE IO_module,       ONLY: write_info
  implicit none
   REAL                  :: Ki
   REAL                  :: neutron_velocity
   REAL                  :: neutron_energy
   REAL                  :: neutron_temperature
   REAL                  :: X_energy

   if(debug_proc%level_2)  call write_debug_proc_level(2, "write_wave_features")

   Ki = 2.*pi/wavelength
   neutron_velocity    =  3.956 / wavelength              ! vitesse des neutrons
   neutron_energy      = (9.045 / wavelength) ** 2        ! energie des neutrons
   neutron_temperature = (30.81 / wavelength) ** 2        ! temperature des neutrons
   X_energy =  12.398 / wavelength                        ! energie des rayons-X

   if(.not. ON_SCREEN .or. .not. write_details) return

   !if(DEVICE%diffracto(1:1) /= "?") then
   ! call write_info('')
   ! call write_info('  > DIFFRACTOMETER: '// trim(DEVICE%diffracto))
   !endif
   !if(DEVICE%lab(1:1) /= "?") then
   ! call write_info('')
   ! call write_info('  > Laboratory: '// trim(DEVICE%lab))
   !endif

   call write_info('')
   WRITE(message_text,'( a,F10.5)') '  > WAVELENGTH (A):', wavelength
   call write_info(TRIM(message_text))
   if (keyword_BEAM) then
    IF(neutrons) then
     WRITE(message_text, '(10x,a,F10.5,a)') '. neutron velocity    [3.956/lambda]      : ', neutron_velocity, ' Km/s'
     call write_info(TRIM(message_text))
     WRITE(message_text, '(10x,a,F10.5,a)') '. neutron energy      [(9.045/lambda)**2] : ', neutron_energy,   ' meV'
     call write_info(TRIM(message_text))
     WRITE(message_text, '(10x,a,F10.5,a)') '. neutron temperature [(30.81/lambda)**2] : ', neutron_temperature, ' K'
     call write_info(TRIM(message_text))
     WRITE(message_text, '(10x,a,F10.5,a)') '. neutron wave vector [2pi/lambda]        : ', Ki, ' A-1'
     call write_info(TRIM(message_text))
    else
     WRITE(message_text, '(10x,a,F10.5,a)') '. X-rays energy      [12.398/lambda]      : ', X_energy,   ' KeV'
     call write_info(TRIM(message_text))
     WRITE(message_text, '(10x,a,F10.5,a)') '. X-rays wave vector [2pi/lambda]         : ', Ki, ' A-1'
     call write_info(TRIM(message_text))
    endif
   else
     WRITE(message_text, '(10x,a,F10.5,a)') '. wave vector         [2pi/lambda]        : ', Ki, ' A-1'
     call write_info(TRIM(message_text))
     WRITE(message_text, '(10x,a,F10.5,a)') '. neutron velocity    [3.956/lambda]      : ', neutron_velocity, ' Km/s'
     call write_info(TRIM(message_text))
     WRITE(message_text, '(10x,a,F10.5,a)') '. neutron energy      [(9.045/lambda)**2] : ', neutron_energy,   ' meV'
     call write_info(TRIM(message_text))
     WRITE(message_text, '(10x,a,F10.5,a)') '. neutron temperature [(30.81/lambda)**2] : ', neutron_temperature, ' K'
     call write_info(TRIM(message_text))
     WRITE(message_text, '(10x,a,F10.5,a)') '. X-rays energy       [12.398/lambda]     : ', X_energy,   ' KeV'
     call write_info(TRIM(message_text))
   end if
   call write_info('')


 RETURN
end subroutine write_wave_features

!-----------------------------------------------------------------------------------------------------------------------
subroutine write_device_features()
 USE cryscalc_module, ONLY: ON_SCREEN, wavelength, message_text, DEVICE, debug_proc
 USE IO_module,       ONLY: write_info
  implicit none

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_device_features")

   if(.not. ON_SCREEN) return

   if(DEVICE%diffracto(1:1) /= "?") then
    call write_info('')
    call write_info('  > DIFFRACTOMETER: '// trim(DEVICE%diffracto))
   endif
   if(DEVICE%lab(1:1) /= "?") then
    call write_info('')
    call write_info('  > Laboratory    : '// trim(DEVICE%lab))
   endif

   call write_info('')
   WRITE(message_text,'( a,F10.5)') '  > WAVELENGTH (A):', wavelength
   call write_info(TRIM(message_text))
   call write_info('')


 RETURN
end subroutine write_device_features


!-------------------------------------------------------------------------
subroutine write_BEAM_features
 USE cryscalc_module, ONLY: ON_SCREEN, beam_type, debug_proc
 USE IO_module,       ONLY: write_info
  implicit none

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_beam_features")

  if(.not. ON_SCREEN) return
  call write_info('')
  if(beam_type(1:8) == 'neutrons') then
   call write_info(' > Incident radiation beam : neutrons')
  elseif(beam_type(1:9) == 'electrons') then
   call write_info(' > Incident radiation beam : electrons')
  else
   call write_info(' > Incident radiation beam : X-rays')
  endif

  call write_wave_features
 return
end subroutine write_BEAM_features
!-------------------------------------------------------------------------

subroutine write_Xrays_wavelength()
 USE cryscalc_module,           ONLY : message_text, debug_proc
 USE IO_module,                 ONLY : write_info
 USE CFML_Scattering_Chemical_Tables
 USE wavelength_module
 implicit none
  INTEGER           :: i
  REAL              :: mean_Ka

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_xrays_wavelength")

  call write_info('')
  call write_info('  >> Main Xrays wavelength (A):')
  call write_info('')
  call write_info('             Ka1         Ka2        <Ka>       <Kb1>')
  call write_info('')



  do i=1, tabulated_target_nb
   IF(.NOT. X_target(i)%WRITE) cycle
   mean_Ka = (2* X_target(i)%wave(1) + X_target(i)%wave(2))/3
   WRITE(message_text,'(2x,a2,4(2x,F10.5))') X_target(i)%label, X_target(i)%wave(1:2), mean_Ka, X_target(i)%wave(3)
   call write_info(TRIM(message_text))
 end do

   call write_info('')

 return
end subroutine write_Xrays_wavelength

!-------------------------------------------------------------------------
 subroutine Write_QVEC
  USE IO_module,               ONLY : write_info
  USE cryscalc_module,         ONLY : message_text, qvec, debug_proc
  implicit none

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_QVEC")

   call write_info('')
   write(message_text, '(a,3F6.2)') '  > Modulation wave vector Qvec = ', qvec(1:3)
   call write_info(trim(message_text))
   call write_info('')

  return
 end subroutine Write_QVEC


!-------------------------------------------------------------------------
 subroutine Write_SUPERCELL
  USE IO_module,               ONLY : write_info
  USE cryscalc_module,         ONLY : message_text, keyword_CELL, sup_cell, unit_cell, crystal_cell,write_details, debug_proc
  USE CFML_Crystal_Metrics,    ONLY : Set_Crystal_Cell
  implicit none
  integer           :: coef_vol

  coef_vol = sup_cell(1) * sup_cell(2) * sup_cell(3)

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_SUPERCELL")

  if(write_details) then
  call write_info('')
  WRITE(message_text, '(a,3(I2,a))') '  > SUPERCELL: ', sup_cell(1), 'a * ', sup_cell(2), 'b * ', sup_cell(3), 'c'
  call write_info(trim(message_text))
  call write_info('')
  end if

  if(keyword_CELL) then
   unit_cell%param(1) = sup_cell(1)*unit_cell%param(1)
   unit_cell%param(2) = sup_cell(2)*unit_cell%param(2)
   unit_cell%param(3) = sup_cell(3)*unit_cell%param(3)
   unit_cell%volume   = coef_vol*unit_cell%volume

   if(write_details) then
   WRITE(message_text, '(a,F15.5, a)') '           a = ', unit_cell%param(1), ' A'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F15.5, a)') '           b = ', unit_cell%param(2), ' A'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F15.5, a)') '           c = ', unit_cell%param(3), ' A'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F15.5, a)') '        alfa = ', unit_cell%param(4), ' deg.'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F15.5, a)') '        beta = ', unit_cell%param(5), ' deg.'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F15.5, a)') '       gamma = ', unit_cell%param(6), ' deg.'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F15.5, a)') '      Volume = ', unit_cell%volume, ' A3 '
   call write_info(TRIM(message_text))
   end if

   call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)

  end if

  return
 end subroutine Write_SUPERCELL


!-------------------------------------------------------------------------
subroutine write_ZUNIT_features
 USE cryscalc_module, ONLY: ON_SCREEN, molecule, message_text, debug_proc
 USE IO_module,       ONLY: write_info
  implicit none

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_Zunit_features")

  if(.not. ON_SCREEN) return
  call write_info('')
   write(message_text, '(a,F3.0)') '  > Z unit number of molecular unit) = ', molecule%Z_unit
  call write_info(trim(message_text))
  call write_info('')

 return
end subroutine  write_ZUNIT_features

!-------------------------------------------------------------------------
subroutine write_space_group(n_out, i1, i2, n_enantio, n_chiral, n_polar)
 use CFML_crystallographic_symmetry, only : set_spacegroup
 USE cryscalc_module,           ONLY : list_sg, list_sg_centric, list_sg_multip, list_sg_enantio, list_sg_chiral, &
                                       list_sg_polar, message_text, SPG, debug_proc
 USE IO_module,                 ONLY : write_info
 implicit none
  INTEGER, INTENT(INOUT)    :: n_out
  INTEGER, INTENT(IN)       :: i1, i2    ! numero du groupe d'espace
  INTEGER, INTENT(INOUT)    :: n_enantio ! nombre de groupes d'espace enantiomorphes repondant au critere recherche
  INTEGER, INTENT(INOUT)    :: n_chiral  ! nombre de groupes d'espace chiraux repondant au critere recherche
  INTEGER, INTENT(INOUT)    :: n_polar   ! nombre de groupes d'espace polaires repondant au critere recherche
  CHARACTER (LEN=3)         :: string_i
  CHARACTER (LEN=12)        :: PG_string, Laue_string
  CHARACTER (LEN=1)         :: enantio_string, chiral_string, polar_string
  INTEGER                   :: i, n
  LOGICAL                   :: ok, enantio, chiral, polar

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_space_group")

  ! numero du groupe d'espace:          SPG%NumSpg
  ! symbole:                            SPG%SPG_Symb
  ! centro:                             SPG%Centred  =0 Centric(-1 no at origin)
  !                                                  =1 Acentric
  !                                                  =2 Centric(-1 at origin)
  n=0
  !n_enantio = 0   ! initialise dans la routine list_space_groups (cryscalc_lsg_cfml.F90)
  !n_chiral  = 0
  !n_polar   = 0

  do i = i1, i2
   WRITE(string_i, '(i3)') i
   call set_spacegroup(string_i, SPG)

   IF(list_sg_centric(1) .AND. SPG%centred==1) cycle
   IF(list_sg_centric(2) .AND. SPG%centred/=1) cycle

   call test_Bravais(SPG%SPG_symb, ok)
   IF(.NOT. ok) cycle

   call test_laue(SPG%Laue, ok)
   IF(.NOT. ok) cycle

   call test_enantio(i, enantio)
   if(enantio) n_enantio = n_enantio + 1
   if(.not. enantio .and. list_sg_enantio) cycle

   call test_chiral(i, chiral)
   if(chiral) n_chiral = n_chiral + 1
   if(.not. chiral .and. list_sg_chiral) cycle

   !call test_polar(i, polar)
   call test_polar_2(SPG%PG, polar)
   if(polar) n_polar = n_polar + 1
   if(.not. polar .and. list_sg_polar) cycle


   if(.not. enantio .or. list_sg_enantio) then
    enantio_string = ' '
   else
    enantio_string = '*'
   endif

   if(.not. chiral .or. list_sg_chiral) then
    chiral_string = ' '
   else
    chiral_string = '+'
   endif

   if(.not. polar .or. list_sg_polar) then
    polar_string = ' '
   else
    polar_string ='§'
   endif

   write(Laue_string, '(3a)') '(', trim(SPG%laue), ')'
   write(PG_string,   '(3a)') '[', trim(SPG%PG), ']'

   n_out = n_out + 1
   IF(list_sg_multip) then
    WRITE(message_text, '(5x,I4,5x,a,I3,a,5x,4a,10x,a,5x,2a,I3)') n_out, 'IT# ', i,'.', chiral_string(1:1),enantio_string(1:1), &
                                                             polar_string(1:1), SPG%SPG_symb(1:10) , Laue_string(1:10),         &
                                                             PG_string(1:8), '    mult.=', SPG%multip


   else
    WRITE(message_text, '(5x,I4,5x,a,I3,a,5x,4a,10x,a,5x,a)') n_out, 'IT# ', i,'.', chiral_string(1:1),enantio_string(1:1),  &
                                                             polar_string(1:1), SPG%SPG_symb(1:10) , Laue_string(1:10),      &
                                                            PG_string(1:8)

   endif
   call write_info(TRIM(message_text))
   n=n+1
  END do

  IF(n==0) then
   call write_info('  >> No space group founded with the required conditions.')
   call write_info('')
  endif



 RETURN

end subroutine write_space_group

!-------------------------------------------------------------
!subroutine write_current_space_group(SG_numor)
subroutine write_current_space_group(SG_symbol)
 USE CRYSCALC_module,                ONLY : keyword_create_CIF
 USE IO_module,                      ONLY : write_info
 USE CFML_crystallographic_symmetry, ONLY : set_spacegroup, searchop, write_sym
 USE cryscalc_module,                ONLY : SPG, message_text, debug_proc
 USE CFML_symmetry_tables,           ONLY : intsymoh,x_oh, intsymd6h,x_d6h

 implicit none
  !INTEGER, INTENT(IN)             :: SG_numor
  CHARACTER (LEN=*), INTENT(IN)   :: SG_symbol

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_current_space_group")

  !WRITE(string_numor, '(i3)') SG_numor
  !call set_spacegroup(string_numor, SPG)
  if(SPG%numSPG == 0) then
   call write_info("  !! Unknown space group !!")
   call write_info("")
   return
  end if

  call set_spacegroup(SG_symbol, SPG)

  call write_info('')
  call write_info('  > Space group:')
  call write_info('')
  WRITE(message_text, '(a,i3)') '   .  Number of space group: ', SPG%NumSPG
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,a)')  '   . Hermann-Mauguin symbol: ', SPG%SPG_symb
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,a)')  '   .       Crystal symmetry: ', SPG%CrystalSys
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,a)')  '   .             Laue class: ', SPG%Laue
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,a)')  '   .           Lattice type: ', SPG%Bravais
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,I3)') '   .   General multiplicity: ', SPG%Multip
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,a)')  '   .         Centrosymmetry: ', SPG%centre
  call write_info(TRIM(message_text))

  if(keyword_create_CIF) then
   call write_CIF_file('SPACE_GROUP')
  end if

 return
end subroutine write_current_space_group

!---------------------------------------------------------------------------------------------
subroutine write_atom_list()
 USE cryscalc_module, ONLY : message_text, nb_atom, keyword_SPGR, keyword_ADP_list,         &
                             atom_label, atom_typ, atom_coord, atom_occ, atom_occ_perc,     &
                             atom_mult, atom_Biso,  atom_occ, atom_Ueq, atom_adp_type,      &
                             SPG, input_PCR, keyword_create_CIF, write_atoms_in_A,          &
                             keyword_read_PCR, therm, ADP_details, ADP_comment,             &
                             atom_adp_aniso, unit_cell, keyword_CELL, CIF_unit, debug_proc
 use CFML_Crystallographic_Symmetry,    ONLY : Get_Multip_Pos
 USE IO_module,                         ONLY : write_info
 USE CFML_GlobalDeps,       ONLY : sp


 implicit none
  INTEGER              :: i, j
  REAL                 :: occ_perc  ! % occ. site
  REAL                 :: f1, f
  REAL(kind=sp), DIMENSION(3)      :: rms
  CHARACTER(len=128)   :: title_Uij, title_Uij_rms



 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_atoms_list")

 call write_info('')
 call write_info('  > ATOMS LIST:')
 call write_info('')

 if(keyword_ADP_list) then
  write(title_Uij, '(a)')      '   NAME     LABEL    U_11      U_22      U_33      U_23      U_13      U_12'
  write(title_Uij_rms, '(2a)') '   NAME     LABEL    U_11      U_22      U_33      U_23      U_13      U_12', &
                               '     rms_U1   rms_U2   rms_U3'
  if(ADP_details) then
   call write_info(trim(title_Uij))
  else
   call write_info(trim(title_Uij_rms))
  end if
 else
  IF(keyword_SPGR) then
   if(input_PCR .or. keyword_read_PCR) then
    if(.not. write_atoms_in_A) then
     call write_info('         NAME     LABEL          x         y         z      Biso       Occ    Occ(%)   Mult')
    else
     call write_info('         NAME     LABEL          x         y         z      Biso       Occ    Occ(%)   Mult' // &
                     '         x(A)      y(A)      z(A)')
    endif
   else
    if(.not. write_atoms_in_A) then
     call write_info('         NAME     LABEL          x         y         z       Beq      Occ     Occ(%)   Mult')
    else
     call write_info('         NAME     LABEL          x         y         z       Beq       Occ    Occ(%)   Mult' // &
                     '         x(A)      y(A)      z(A)')
    endif
   endif
  else
   if(.not. write_atoms_in_A) then
    call write_info( '   NAME     LABEL          x         y         z       Beq       Occ')
   else
    call write_info( '   NAME     LABEL          x         y         z       Beq       Occ         x(A)      y(A)      z(A)')
   endif
  endif
 endif
 call write_info('')



 IF(keyword_SPGR .and. (input_PCR .or. keyword_read_PCR)) then
  atom_mult(1) = Get_Multip_Pos(atom_coord(1:3,1), SPG)
  f1 = atom_occ(1) * REAL(SPG%Multip) / REAL(atom_mult(1))
  f1 = 1./f1
 endif

 !if(keyword_create_CIF) call write_CIF_file('ATOMS_ANISO_HEADER')

 Therm%Uij = .true.
 !ADP_details   = .false.
 do i= 1, nb_atom


 if(keyword_ADP_list ) then
  if(atom_adp_type(i)(1:5) == 'isotr') cycle

  if(ADP_details) then

   !call write_info('')
   write(message_text, '(2x,80a1)')  ('-', j=1,80)
   call write_info(TRIM(message_text))
   write(message_text,'(3x,a4,5x,a4,2x,6F10.5)')    atom_label(i), atom_typ (i) , atom_adp_aniso(1:3,i), &
                                                    atom_adp_aniso(6,i), atom_adp_aniso(5,i), atom_adp_aniso(4,i)
   call write_info(TRIM(message_text))
   write(message_text, '(2x,80a1)')  ('-', j=1,80)
   call write_info(TRIM(message_text))
   call write_info('')
  end if
  call check_adp(atom_adp_aniso(1:6,i), rms, adp_comment, ADP_details)
  if(ADP_details) then
   call write_info(trim(title_Uij_rms))
   call write_info('')
  end if
  write(message_text,'(3x,a4,5x,a4,2x,6F10.5,3(1x,F8.5), 3x, a)')    atom_label(i), atom_typ (i) , atom_adp_aniso(1:3,i), &
                                                     atom_adp_aniso(6,i), atom_adp_aniso(5,i), atom_adp_aniso(4,i), &
                                                     rms(1:3), trim(ADP_comment)
  !if(keyword_create_CIF) then
  ! write(CIF_unit, '(a4,6F10.5)') atom_label(i), atom_adp_aniso(1:3,i), atom_adp_aniso(6,i), atom_adp_aniso(5,i), &
  !                                atom_adp_aniso(4,i)
  !end if
 else


  IF(keyword_SPGR) then
   atom_mult(i) = Get_Multip_Pos(atom_coord(1:3,i), SPG)

   if(input_PCR .or. keyword_read_PCR) then
    f =  SPG%Multip / atom_mult(i)
    occ_perc = atom_occ(i) *  f * f1
    atom_occ_perc(i) =  occ_perc

    if(write_atoms_in_A .and. keyword_CELL) then
     write(message_text,'(3x,i3,3x,a4,5x,a4,2x,6F10.5,3x,I4,3x,3F10.5)')  i, atom_label(i), atom_typ (i) , atom_coord(1:3,i), &
                                                                       atom_Biso(i),  atom_occ(i),    occ_perc, atom_mult(i), &
                                                                       atom_coord(1,i)*unit_cell%param(1), &
                                                                       atom_coord(2,i)*unit_cell%param(2), &
                                                                       atom_coord(3,i)*unit_cell%param(3)
    else
     write(message_text,'(3x,i3,3x,a4,5x,a4,2x,6F10.5,3x,I4)')  i, atom_label(i), atom_typ (i) , atom_coord(1:3,i), &
                                                                atom_Biso(i),  atom_occ(i),    occ_perc, atom_mult(i)
    endif

   else
    !write(message_text,'(3x,i3,3x,a4,5x,a4,2x,5F10.5,3x,I4)')  i, atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
    !                                                           atom_Biso(i),  atom_occ_perc(i),    atom_mult(i)


    atom_occ(i) = real(atom_mult(i))/SPG%multip
    if(write_atoms_in_A .and. keyword_CELL) then
     write(message_text,'(3x,i3,3x,a4,5x,a4,2x,6F10.5,3x,I4,3x,3F10.5)')  i, atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
                                                               atom_Biso(i),  atom_occ(i),  atom_occ_perc(i),  atom_mult(i), &
                                                               atom_coord(1,i)*unit_cell%param(1), &
                                                               atom_coord(2,i)*unit_cell%param(2), &
                                                               atom_coord(3,i)*unit_cell%param(3)
    else
     write(message_text,'(3x,i3,3x,a4,5x,a4,2x,6F10.5,3x,I4)')  i, atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
                                                               atom_Biso(i),  atom_occ(i),  atom_occ_perc(i),  atom_mult(i)
    endif

   endif
  else
   if(write_atoms_in_A .and. keyword_CELL) then
   write(message_text,'(3x,a4,5x,a4,2x,5F10.5,3x,3F10.5)')  atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
                                                             atom_Biso(i),  atom_occ_perc(i),    &
                                                             atom_coord(1,i)*unit_cell%param(1), &
                                                             atom_coord(2,i)*unit_cell%param(2), &
                                                             atom_coord(3,i)*unit_cell%param(3)
   else
    write(message_text,'(3x,a4,5x,a4,2x,5F10.5)')     atom_label(i), atom_typ(i) , atom_coord(1:3,i), &
                                                      atom_Biso(i),  atom_occ_perc(i)
   endif
  endif
 end if
 call write_info(TRIM(message_text))
 if(ADP_details) then
  call write_info('')
  call write_info('')
 end if
 !if(write_atoms_in_A .and. keyword_CELL .and. keyword_read_PCR) then  ! ???? !
 ! write(message_text, '(97x,a)') '0.00000   0.00000   0.00000'
 ! call write_info(message_text)
 !end if

 end do

 IF(keyword_create_CIF) then
  if(keyword_ADP_list) then
   call write_CIF_file('ATOMS_ANISO_HEADER')
   call write_CIF_file('ATOMS_ANISO')
  else
   call write_CIF_file('ATOMS_HEADER')
   call write_CIF_file('ATOMS_HEADER_LOOP')
   call write_CIF_file('ATOM')
  end if
 endif

 return
end subroutine write_atom_list

!-----------------------------------------------------------------------------------------------------
subroutine write_atom_list_cart
 use io_module
 USE cryscalc_module,       ONLY : message_text, nb_atom, atom_label, atom_typ, atom_coord, crystal_cell, &
                                   write_details, debug_proc,                                             &
                                   cartesian_frame, create_SHAPE_file, tmp_unit, tmp_2_unit, main_title
 USE CFML_crystal_metrics, only : cart_vector
 USE CFML_GlobalDeps,      ONLY : sp
 implicit none
  integer         :: i
  real (kind=sp), dimension(3)   :: cart_vect
  character (len=12)             :: string
  character (len=256)            :: file_xyz

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_atoms_list_cart")


 call  create_CELL_object

 call write_info('')
 !call write_info(' Cartesian frame type : '//cartesian_frame%type//' ' //trim(cartesian_frame%string))
 call write_info(' Cartesian frame type : ' //trim(cartesian_frame%string))
 call write_info('')

 if(create_SHAPE_file) then
  open(unit = tmp_unit, file="cryscalc_shape.dat")
  write(unit = tmp_unit, fmt='(a)')    '! Input file for SHAPE, created by CRYSCALC (TR/CDIFX-ISCR Rennes)'
  write(unit = tmp_unit, fmt='(2a)') '$ ', trim(main_title)
  write(unit = tmp_unit, fmt='(a)')  '! Ligands    central atom (line 1 in the atoms list)'
  write(string, fmt='(2(1x,i3))') nb_atom-1, 1
  write(unit = tmp_unit, fmt='(a)') adjustl(string)
  write(unit = tmp_unit, fmt='(a)')  '! Polyedron will be compared to any kind of polyedra'
  write(unit = tmp_unit, fmt='(a)') '1 2 3 4 5 6 7 8 9 10 11 12 13'
  write(unit = tmp_unit, fmt='(a)')  '! Label for the structure'
  if(len_trim(main_title) < 16) then
   write(unit = tmp_unit, fmt='(a)') trim(main_title)
  else
   write(unit = tmp_unit, fmt='(a)') main_title(1:15)
  end if
  write(unit = tmp_unit, fmt='(a)')  '! List of cartesian atomic coordinates'
 end if

 if(cartesian_frame%type == 'A') then
  file_xyz = "cryscalc_cart_A.xyz"
 else
  file_xyz = "cryscalc_cart_C.xyz"
 end if
 open(unit = tmp_2_unit, file=trim(file_xyz))
 write(unit = tmp_2_unit, fmt='(I6)') nb_atom
 write(unit = tmp_2_unit, fmt='(a)') main_title


 do i= 1, nb_atom
  cart_vect(1:3) = cart_vector('D', atom_coord(1:3, i), crystal_cell)

  if(write_details) then
  write(message_text,'(3x,i3,3x,a4,5x,a4,2x,3F12.6)')     i, atom_label(i), atom_typ(i) , cart_vect(1:3)
  call write_info(TRIM(message_text))
  end if

  write(unit=tmp_2_unit, fmt='(a4,3F12.6)') atom_typ(i), cart_vect(1:3)

  if(create_SHAPE_file) write(unit = tmp_unit, fmt='(a4,2x,3F12.6)')     atom_label(i) , cart_vect(1:3)
 end do
 close(unit = tmp_2_unit)

 call write_info('')
 call write_info('   > '//trim(file_xyz)//' file has been created.')
 call write_info('')


  if(create_SHAPE_file) then
   close(unit=tmp_unit)
   call write_info('')
   call write_info('   > cryscalc_shape.dat file for SHAPE has been created.')
   call write_info('')
  end if



 return
end subroutine write_atom_list_cart

!-----------------------------------------------------------------------------------------------------
subroutine write_molecular_features
 use cryscalc_module, only : molecule, message_text, debug_proc
 use IO_module,       only : write_info
 implicit none

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_molecular_features")

 call write_info('')
 call write_info('  >> Molecular features:')
 call write_info('')
 write(message_text, '(2a)')        '    . formula                                     : ', trim(molecule%formula)
 call write_info(trim(message_text))

 if(molecule%Z_unit /=0) then
  write(message_text, '(a, I3)')    '    . number of formula units                     : ', molecule%Z_unit
  call write_info(trim(message_text))
 endif

 IF(molecule%content(1:1) /= '?') then
  write(message_text, '(2a)')       '    . content                                     : ', trim(molecule%content)
  call write_info(trim(message_text))
 endif

 write(message_text, '(a, F8.2)')   '    . weight                                      : ', molecule%weight
 call write_info(trim(message_text))
 if(molecule%density > 0.01) then
  write(message_text, '(a, F9.3)')  '    . density                                     : ', molecule%density
  call write_info(trim(message_text))
 endif

 if(molecule%Z /=0) then
  WRITE(message_text, '(a,I6)')     '    . total number of electrons in the molecule   : ', molecule%Z
  call write_info(TRIM(message_text))
 endif

 call write_info('')

 return
end subroutine write_molecular_features


!-----------------------------------------------------------------------------------------------------
subroutine write_REF(input_string)
 USE cryscalc_module, ONLY : keyword_create_CIF, SW_EVAL, SADABS, SHELX, SIR, SPF, ABS_CRYSALIS, CIF_unit, WRITE_ref_CIF, &
                             write_details, debug_proc
 USE CIF_module,      ONLY : CIF_parameter_DEVICE_KCCD,     CIF_parameter_DEVICE_APEX,      CIF_parameter_DEVICE_X2S,    &
                             CIF_parameter_DEVICE_XCALIBUR, CIF_parameter_DEVICE_SUPERNOVA, CIF_parameter_DEVICE_D8V_Cu, &
                             CIF_parameter_DEVICE_D8V_Mo
 USE IO_module,       ONLY : write_info
 implicit none
 CHARACTER(LEN=*), INTENT(IN) :: input_string
 INTEGER                      :: long_input_string, i, k
 LOGICAL                      :: input_KCCD, input_APEX, input_X2S, input_XCALIBUR, input_SUPERNOVA
 LOGICAL                      :: input_EVAL, input_DENZO, input_SADABS, input_SHELX, input_ABS_CRYSALIS
 LOGICAL                      :: input_SIR, input_SPF
 LOGICAL                      :: input_D8V_Cu, input_D8V_Mo
 LOGICAL                      :: file_opened

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_ref ("//trim(input_string)//")")

 long_input_string = len_trim(input_string)
 input_KCCD          = .false.
 input_APEX          = .false.
 input_X2S           = .false.
 input_SIR           = .false.
 input_SPF           = .false.
 input_XCALIBUR      = .false.
 input_SUPERNOVA     = .false.
 input_EVAL          = .false.
 input_DENZO         = .false.
 input_SADABS        = .false.
 input_SHELX         = .false.
 input_ABS_CRYSALIS  = .false.
 input_D8V_Cu        = .false.
 input_D8V_Mo        = .false.

 if(long_input_string == 3) then
  if(input_string(1:3)  == 'X2S')           input_X2S = .true.
  if(input_string(1:3)  == 'SIR')           input_SIR = .true.
  if(input_string(1:3)  == 'SPF')           input_SPF = .true.
 elseif(long_input_string == 4) then
  if(input_string(1:4)  == 'KCCD')          input_KCCD = .true.
  if(input_string(1:4)  == 'APEX')          input_APEX = .true.
  if(input_string(1:4)  == 'EVAL')          input_EVAL = .true.
 elseif(long_input_string == 5) then
  if(input_string(1:5)  == 'DENZO')         input_DENZO = .true.
  if(input_string(1:5)  == 'SHELX')         input_SHELX = .true.
 elseif(long_input_string == 6) then
  if(input_string(1:6)  == 'SADABS')        input_SADABS = .true.
  if(input_string(1:6)  == 'D8V_CU')        input_D8V_Cu = .true.
  if(input_string(1:6)  == 'D8V_MO')        input_D8V_Mo = .true.
 elseif(long_input_string == 8)  then
  if(input_string(1:8)  == 'XCALIBUR')      input_XCALIBUR = .true.
 elseif(long_input_string == 9)  then
  if(input_string(1:9)  == 'SUPERNOVA')     input_SUPERNOVA    = .true.
 elseif(long_input_string == 12)  then
  if(input_string(1:12) == 'ABS_CRYSALIS')  input_ABS_CRYSALIS = .true.
 endif


 IF(input_KCCD) then
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   DATA COLLECTION                                          #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  call write_info("_diffrn_source                          "//trim(CIF_parameter_DEVICE_KCCD%diffrn_source))
  call write_info("_diffrn_radiation_wavelength            "//trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_wavelength))
  call write_info("_diffrn_radiation_type                  "//trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_type))
  call write_info("_diffrn_radiation_source                "//trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_source))
  call write_info("_diffrn_radiation_monochromator         "//trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_monochromator))
  call write_info("_diffrn_radiation_probe                 "//trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_probe))
  call write_info("_diffrn_detector                        "//trim(CIF_parameter_DEVICE_KCCD%diffrn_detector))
  call write_info("_diffrn_detector_area_resol_mean        "//trim(CIF_parameter_DEVICE_KCCD%diffrn_detector_area_resol_mean))
  call write_info("_diffrn_measurement_device              "//trim(CIF_parameter_DEVICE_KCCD%diffrn_measurement_device))
  call write_info("_diffrn_measurement_device_type         "//trim(CIF_parameter_DEVICE_KCCD%diffrn_measurement_device_type))
  call write_info("_diffrn_measurement_method              "//trim(CIF_parameter_DEVICE_KCCD%diffrn_measurement_method))

  IF(keyword_create_CIF .or. WRITE_ref_CIF) then
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   !call write_CIF_file("KCCD")
   call write_CIF_section("DATA COLLECTION")
   !call test_file_open(unit=CIF_unit)
   WRITE(CIF_unit, '(2a)')"_diffrn_source                         ", trim(CIF_parameter_DEVICE_KCCD%diffrn_source)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_wavelength           ", trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_wavelength)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_type                 ", trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_type)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_probe                ", trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_probe)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_source               ", trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_source)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_monochromator        ", trim(CIF_parameter_DEVICE_KCCD%diffrn_radiation_monochromator)
   WRITE(CIF_unit, '(2a)')"_diffrn_detector                       ", trim(CIF_parameter_DEVICE_KCCD%diffrn_detector)
   WRITE(CIF_unit, '(2a)')"_diffrn_detector_area_resol_mean       ", trim(CIF_parameter_DEVICE_KCCD%diffrn_detector_area_resol_mean)
   WRITE(CIF_unit, '(2a)')"_diffrn_measurement_device             ", trim(CIF_parameter_DEVICE_KCCD%diffrn_measurement_device)
   WRITE(CIF_unit, '(2a)')"_diffrn_measurement_device_type        ", trim(CIF_parameter_DEVICE_KCCD%diffrn_measurement_device_type)
   WRITE(CIF_unit, '(2a)')"_diffrn_measurement_method             ", trim(CIF_parameter_DEVICE_KCCD%diffrn_measurement_method)
   write(CIF_unit, '(a)') ''
  endif


 ELSEIF(input_APEX) then
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_DEVICE_APEX%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_DEVICE_APEX%diffrn_measurement_method))
   call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_probe))
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_DEVICE_APEX%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_DEVICE_APEX%diffrn_detector_area_resol_mean))

   IF(keyword_create_CIF .or. WRITE_ref_CIF)   then
   !call write_CIF_file('APEX')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   call write_CIF_section("DATA COLLECTION")
   WRITE(CIF_unit, '(2a)') "_diffrn_measurement_device_type    ", trim(CIF_parameter_DEVICE_APEX%diffrn_measurement_device_type)
   WRITE(CIF_unit, '(2a)') "_diffrn_measurement_method         ", trim(CIF_parameter_DEVICE_APEX%diffrn_measurement_method)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_wavelength       ", trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_wavelength)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_type             ", trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_type)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_source           ", trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_source)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_monochromator    ", trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_monochromator)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_probe            ", trim(CIF_parameter_DEVICE_APEX%diffrn_radiation_probe)
   WRITE(CIF_unit, '(2a)') "_diffrn_detector                   ", trim(CIF_parameter_DEVICE_APEX%diffrn_detector)
   WRITE(CIF_unit, '(2a)') "_diffrn_detector_area_resol_mean   ", trim(CIF_parameter_DEVICE_APEX%diffrn_detector_area_resol_mean)
   write(CIF_unit, '(a)') ''
  endif

 ELSEIF(input_X2S) then
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_DEVICE_X2S%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_DEVICE_X2S%diffrn_measurement_method))
   call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_probe))
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_DEVICE_X2S%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_DEVICE_X2S%diffrn_detector_area_resol_mean))

  IF(keyword_create_CIF .or. WRITE_ref_CIF)   then
   !call write_CIF_file('X2S')
    inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   call write_CIF_section("DATA COLLECTION")
   WRITE(CIF_unit, '(2a)') "_diffrn_measurement_device_type    ", trim(CIF_parameter_DEVICE_X2S%diffrn_measurement_device_type)
   WRITE(CIF_unit, '(2a)') "_diffrn_measurement_method         ", trim(CIF_parameter_DEVICE_X2S%diffrn_measurement_method)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_wavelength       ", trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_wavelength)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_type             ", trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_type)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_source           ", trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_source)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_monochromator    ", trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_monochromator)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_probe            ", trim(CIF_parameter_DEVICE_X2S%diffrn_radiation_probe)
   WRITE(CIF_unit, '(2a)') "_diffrn_detector                   ", trim(CIF_parameter_DEVICE_X2S%diffrn_detector)
   WRITE(CIF_unit, '(2a)') "_diffrn_detector_area_resol_mean   ", trim(CIF_parameter_DEVICE_X2S%diffrn_detector_area_resol_mean)
   write(CIF_unit, '(a)') ''
  endif

 ELSEIF(input_D8V_CU) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   DATA COLLECTION                                          #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_measurement_device_type))
  call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_measurement_method))
  call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_wavelength))
  call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_type))
  call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_source))
  call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_monochromator))
  call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_probe))
  call write_info("_diffrn_detector                  "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_detector))
  !call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_DEVICE_D8V_CU%diffrn_detector_area_resol_mean))
  IF(keyword_create_CIF .or. WRITE_ref_CIF)   then
   !call write_CIF_file('D8V_CU')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   call write_CIF_section("DATA COLLECTION")
   WRITE(CIF_unit, '(2a)') "_diffrn_measurement_device_type    ", trim(CIF_parameter_DEVICE_D8V_CU%diffrn_measurement_device_type)
   WRITE(CIF_unit, '(2a)') "_diffrn_measurement_method         ", trim(CIF_parameter_DEVICE_D8V_CU%diffrn_measurement_method)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_wavelength       ", trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_wavelength)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_type             ", trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_type)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_source           ", trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_source)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_monochromator    ", trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_monochromator)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_probe            ", trim(CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_probe)
   WRITE(CIF_unit, '(2a)') "_diffrn_detector                   ", trim(CIF_parameter_DEVICE_D8V_CU%diffrn_detector)
   write(CIF_unit, '(a)') ''
  endif

 ELSEIF(input_D8V_MO) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   DATA COLLECTION                                          #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_measurement_device_type))
  call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_measurement_method))
  call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_wavelength))
  call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_type))
  call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_source))
  call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_monochromator))
  call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_probe))
  call write_info("_diffrn_detector                  "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_detector))
  !call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_DEVICE_D8V_MO%diffrn_detector_area_resol_mean))

  IF(keyword_create_CIF .or. WRITE_ref_CIF)   then
   !call write_CIF_file('D8V_MO')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   call write_CIF_section("DATA COLLECTION")
   WRITE(CIF_unit, '(2a)') "_diffrn_measurement_device_type    ", trim(CIF_parameter_DEVICE_D8V_Mo%diffrn_measurement_device_type)
   WRITE(CIF_unit, '(2a)') "_diffrn_measurement_method         ", trim(CIF_parameter_DEVICE_D8V_Mo%diffrn_measurement_method)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_wavelength       ", trim(CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_wavelength)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_type             ", trim(CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_type)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_source           ", trim(CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_source)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_monochromator    ", trim(CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_monochromator)
   WRITE(CIF_unit, '(2a)') "_diffrn_radiation_probe            ", trim(CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_probe)
   WRITE(CIF_unit, '(2a)') "_diffrn_detector                   ", trim(CIF_parameter_DEVICE_D8V_Mo%diffrn_detector)
   write(CIF_unit, '(a)') ''
  endif


 ELSEIF(input_XCALIBUR) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   DATA COLLECTION                                          #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_measurement_device_type))
  call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_measurement_method))
  call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_wavelength))
  !call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_type))
  call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_source))
  call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_monochromator))
  call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_probe))
  call write_info("_diffrn_detector                  "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_detector))
  call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_detector_area_resol_mean))

  IF(keyword_create_CIF .or. WRITE_ref_CIF)   then
   !call write_CIF_file('XCALIBUR')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   call write_CIF_section("DATA COLLECTION")
   WRITE(CIF_unit, '(2a)')"_diffrn_measurement_device_type   ", trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_measurement_device_type)
   WRITE(CIF_unit, '(2a)')"_diffrn_measurement_method        ", trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_measurement_method)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_wavelength      ", trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_wavelength)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_source          ", trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_source)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_monochromator   ", trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_monochromator)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_probe           ", trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_probe)
   WRITE(CIF_unit, '(2a)')"_diffrn_detector                  ", trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_detector)
   WRITE(CIF_unit, '(2a)')"_diffrn_detector_area_resol_mean  ", trim(CIF_parameter_DEVICE_XCALIBUR%diffrn_detector_area_resol_mean)
   write(CIF_unit, '(a)') ''
  endif

 ELSEIF(input_SUPERNOVA) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   DATA COLLECTION                                          #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_measurement_device_type))
  call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_measurement_method))
  call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_wavelength))
  call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_type))
  call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_source))
  call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_monochromator))
  call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_probe))
  call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_source))
  call write_info("_diffrn_detector                  "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_detector))
  call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_DEVICE_SUPERNOVA%diffrn_detector_area_resol_mean))

  IF(keyword_create_CIF .or. WRITE_ref_CIF)   then
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   call write_CIF_section("DATA COLLECTION")
   WRITE(CIF_unit, '(2a)')"_diffrn_measurement_device_type   ", trim(CIF_parameter_DEVICE_Supernova%diffrn_measurement_device_type)
   WRITE(CIF_unit, '(2a)')"_diffrn_measurement_method        ", trim(CIF_parameter_DEVICE_Supernova%diffrn_measurement_method)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_wavelength      ", trim(CIF_parameter_DEVICE_Supernova%diffrn_radiation_wavelength)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_source          ", trim(CIF_parameter_DEVICE_Supernova%diffrn_radiation_source)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_monochromator   ", trim(CIF_parameter_DEVICE_Supernova%diffrn_radiation_monochromator)
   WRITE(CIF_unit, '(2a)')"_diffrn_radiation_probe           ", trim(CIF_parameter_DEVICE_Supernova%diffrn_radiation_probe)
   WRITE(CIF_unit, '(2a)')"_diffrn_detector                  ", trim(CIF_parameter_DEVICE_Supernova%diffrn_detector)
   WRITE(CIF_unit, '(2a)')"_diffrn_detector_area_resol_mean  ", trim(CIF_parameter_DEVICE_Supernova%diffrn_detector_area_resol_mean)
   write(CIF_unit, '(a)') ''
  endif


 ELSEIF(input_EVAL) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   COMPUTER PROGRAMS USED                                   #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  !call write_info("_computing_data_collection      "//trim(CIF_parameter_DEVICE_KCCD%computing_data_collection))
  !call write_info("_computing_cell_refinement      "//trim(CIF_parameter_DEVICE_KCCD%computing_cell_refinement))
  !call write_info("_computing_data_reduction       "//trim(CIF_parameter_DEVICE_KCCD%computing_data_reduction))
  call write_info("_computing_data_collection      "//trim(SW_EVAL%data_collection))
  call write_info("_computing_cell_refinement      "//trim(SW_EVAL%cell_refinement))
  call write_info("_computing_data_reduction       "//trim(SW_EVAL%data_collection))
  call write_info("")

  IF(keyword_create_CIF .or. WRITE_ref_CIF) then
   !call write_CIF_file('EVAL_PROGRAMS')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   call write_CIF_section("COMPUTER PROGRAMS USED")
   WRITE(CIF_unit, '(2a)') "_computing_data_collection     ", trim(SW_EVAL%data_collection)
   WRITE(CIF_unit, '(2a)') "_computing_cell_refinement     ", trim(SW_EVAL%cell_refinement)
   WRITE(CIF_unit, '(2a)') "_computing_data_reduction      ", trim(SW_EVAL%data_collection)
   write(CIF_unit, '(a)') ''
  endif

 ELSEIF(input_DENZO) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   COMPUTER PROGRAMS USED                                   #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  call write_info("_computing_data_collection      "//trim(CIF_parameter_DEVICE_KCCD%computing_data_collection))
  call write_info("_computing_cell_refinement      "//trim(CIF_parameter_DEVICE_KCCD%computing_cell_refinement))
  call write_info("_computing_data_reduction       "//trim(CIF_parameter_DEVICE_KCCD%computing_data_reduction))
  call write_info("")

  IF(keyword_create_CIF .or. WRITE_ref_CIF)  then
   !call write_CIF_file('DENZO_PROGRAMS')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   call write_CIF_section("COMPUTER PROGRAMS USED")
   WRITE(CIF_unit, '(2a)') "_computing_data_collection     ", trim(CIF_parameter_DEVICE_KCCD%computing_data_collection)
   WRITE(CIF_unit, '(2a)') "_computing_cell_refinement     ", trim(CIF_parameter_DEVICE_KCCD%computing_cell_refinement)
   WRITE(CIF_unit, '(2a)') "_computing_data_reduction      ", trim(CIF_parameter_DEVICE_KCCD%computing_data_reduction)
   write(CIF_unit, '(a)') ''
  endif

 ELSEIF(input_SADABS) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   ABSORPTION CORRECTION                                    #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  call write_info(trim(SADABS%type))
  do i=1, 4
   call write_info(trim(SADABS%details(i)))
  end do

  if(keyword_create_CIF .or. WRITE_ref_CIF) then
   !call write_CIF_file('SADABS')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   WRITE(CIF_unit, '(a)') trim(SADABS%type)
   do i=1, 4
    WRITE(CIF_unit, '(a)') trim(SADABS%details(i))
   end do
   write(CIF_unit, '(a)') ''
  endif

 ELSEIF(input_SHELX) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                   COMPUTER PROGRAMS USED                                   #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  do i=1, 6
   call write_info(trim(SHELX%details(i)))
  end do
  if(write_details) then
   do k=1, 3
    call write_info('') 
    do i=10+4*(k-1), 13+4*(k-1)
     call write_info(trim(SHELX%details(i)))
    end do
   end do	
  end if


  if(keyword_create_CIF .or. WRITE_ref_CIF) then
   !call write_CIF_file('SHELX')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
   do i=1, 6
    WRITE(CIF_unit, '(a)') trim(SHELX%details(i))
   end do
   write(CIF_unit, '(a)') ''
  endif
  
 ELSEIF(input_SIR) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                        SIR TEAM PROGRAMS                                   #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  if(write_details) then
   do k=1, 3
    call write_info('') 
    do i=1+4*(k-1), 4+4*(k-1)
     call write_info(trim(SIR%details(i)))
    end do
   end do	
  else
   call write_info('') 
   do k=1, 3
     call write_info(trim(SIR%details(1+4*(k-1))))
   end do	
 
  end if

 ELSEIF(input_SPF) then
  call write_info("")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("#                        SUPERFLIP PROGRAM                                   #")
  call write_info("#----------------------------------------------------------------------------#")
  call write_info("")
  if(write_details) then
   do k=1, 1
    call write_info('') 
    do i=1+4*(k-1), 4+4*(k-1)
     call write_info(trim(SPF%details(i)))
    end do
   end do	
  else
   do k=1, 1
     call write_info(trim(SPF%details(1+4*(k-1))))
   end do	
  end if


  ELSEIF(input_ABS_CRYSALIS) then
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   ABSORPTION CORRECTION                                    #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info(trim(ABS_CRYSALIS%type))
   do i=1, 6
    call write_info(trim(ABS_CRYSALIS%details(i)))
   end do

   if(keyword_create_CIF .or. WRITE_ref_CIF) then
    !call write_CIF_file('ABS_CRYSALIS')
   inquire(file="crycalc.cif", opened=file_opened)
   if(file_opened) then
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif', position="append")
   else
    OPEN(UNIT=CIF_unit, FILE='cryscalc.cif')
   end if
    WRITE(CIF_unit, '(a)') trim(ABS_CRYSALIS%type)
    do i=1, 6
     WRITE(CIF_unit, '(a)') trim(ABS_CRYSALIS%details(i))
    end do
    write(CIF_unit, '(a)') ''
   endif
  end if


 return
END subroutine write_REF

!------------------------------------------------------------------------------------------------------------------------------------------
subroutine write_crystal_cell_Cart
 use cryscalc_module,      only : tmp_unit, crystal_cell
 use CFML_crystal_metrics, only : write_crystal_cell
 USE IO_module,            ONLY: write_info

 implicit none
 integer                        :: i_error
 character (len=256)            :: read_line

 call create_cell_object
 open(unit=tmp_unit, file='cfml_cell.out')
   call write_crystal_cell(crystal_cell, tmp_unit)
 close(unit=tmp_unit)

 open(unit=tmp_unit, file='cfml_cell.out')
  do
   read(unit=tmp_unit, fmt='(a)', iostat=i_error) read_line
   if(i_error /=0) exit
   call write_info(' '//trim(read_line))
  end do
 close(unit=tmp_unit)

 call system("del cfml_cell.out")

 return
end subroutine write_crystal_cell_Cart

!------------------------------------------------------------------------------------------------------------------------------
subroutine Generate_atoms_in_supercell(n_at, at_coord, at_label, at_typ, at_biso, at_occ_perc)
 use Cryscalc_module, only : sup_cell, SUPERCELL_pcr, keyword_WRITE_SUPERCELL, nb_atom, nb_atom_max, atom_label, &
                             atom_typ, atom_coord, atom_orbit, atom_mult, atom_occ_perc, atom_Biso, SPG,         &
                             message_text, write_details, debug_proc
 Use CFML_Crystallographic_Symmetry,  only: Set_SpaceGroup
 USE IO_module,       ONLY : write_info
 implicit none
  integer,                              intent(in)    :: n_at
  real,              dimension(3,n_at), intent(in)    :: at_coord
  character (len=*), dimension(n_at),   intent(inout) :: at_label
  character (len=*), dimension(n_at),   intent(in)    :: at_typ
  real,              dimension(n_at),   intent(in)    :: at_biso
  real,              dimension(n_at),   intent(in)    :: at_occ_perc
  integer                 :: i, n, nt, n_total, m1, m2, m3
  real, dimension(3)      :: sup_coord, new_coord


  logical, dimension(3) :: sup

  if(debug_proc%level_2)  call write_debug_proc_level(2, "Generate_atoms_in_SUPERCELL")

  sup = .false.
  do i=1, 3
   if(sup_cell(i) > 1) sup(i) = .true.
  end do


  !if(.not. sup(1) .and. .not. sup(2) .and. .not. sup(3)) then
  ! call write_info('')
  ! call write_info('    !! No SUPERCELL !!')
  ! call write_info('')
  ! return
  !end if

  n_total = n_at * sup_cell(1) * sup_cell(2) * sup_cell(3)

  if(write_details) then
   call write_info('')
   call write_info('  > Atoms in SUPERCELL: ')
   call write_info('')
   write(message_text, '(a,i6)')  ' Number of atoms in the supercell: ', n_total
   call write_info(trim(message_text))
   call write_info('')
  end if

  nt = 0
  do i=1, n_at
   n = 0
   new_coord(1:3) =  at_coord(1:3,i)/sup_cell(1:3)

   m1 = 0
   if(sup_cell(1) > 1) then
    do
     m1 = m1 + 1
     if(m1 > sup_cell(1)) exit
     n  = n  + 1
     nt = nt + 1
     if(nt > nb_atom_max) then
      if(write_details) then
      call write_info('')
      write(message_text, '(a,i4,a)') ' > Number of atoms in supercell is greater than allowed limit (', nb_atom_max,') !!'
      call write_info(trim(message_text))
      call write_info('')
      end if
      return
     end if

     sup_coord(1) = new_coord(1) + (m1-1) * 1./sup_cell(1)
     sup_coord(2) = at_coord(2,i)/sup_cell(2)
     sup_coord(3) = at_coord(3,i)/sup_cell(3)
     call write_SUPERCELL_coord(n, nt, at_label(i), at_typ(i), at_occ_perc(i), sup_coord, '1')
     call Get_new_label(at_label(i), n, atom_label(nt))
     atom_typ(nt)        = at_typ(i)
     atom_coord(1:3, nt) = sup_coord(1:3)
     atom_Biso(nt)       = at_biso(i)
     atom_occ_perc(nt)   = at_occ_perc(i)


     m2 = 1
     if(sup_cell(2) > 1) then
      do
       m2 = m2 + 1
       if(m2 > sup_cell(2)) exit
        n  = n  + 1
        nt = nt + 1
        if(nt > nb_atom_max) then
         if(write_details) then
         call write_info('')
         write(message_text, '(a,i4,a)') ' > Number of atoms in supercell is greater than allowed limit (', nb_atom_max,') !!'
         call write_info(trim(message_text))
         call write_info('')
         end if
         return
        end if
        sup_coord(1) = new_coord(1) + (m1-1) * 1./sup_cell(1)
        sup_coord(2) = new_coord(2) + (m2-1) * 1./sup_cell(2)
        sup_coord(3) = at_coord(3,i)/sup_cell(3)
        call write_SUPERCELL_coord(n, nt, at_label(i), at_typ(i), at_occ_perc(i), sup_coord, '1+2')
        call Get_new_label(at_label(i), n, atom_label(nt))
        !atom_label(nt) = at_label(i)
        atom_typ(nt)   = at_typ(i)
        atom_coord(1:3, nt) = sup_coord(1:3)
        atom_Biso(nt)  = at_biso(i)
        atom_occ_perc(nt)   = at_occ_perc(i)

        m3 = 1
        if(sup_cell(3) > 1) then
         do
          m3 = m3 + 1
          if(m3 > sup_cell(3)) exit
          n  = n  + 1
          nt = nt + 1
          if(nt > nb_atom_max) then
           if(write_details) then
           call write_info('')
           write(message_text, '(a,i4,a)') ' > Number of atoms in supercell is greater than allowed limit (', nb_atom_max,') !!'
           call write_info(trim(message_text))
           call write_info('')
           end if
           return
          end if
          sup_coord(1) = new_coord(1) + (m1-1) * 1./sup_cell(1)
          sup_coord(2) = new_coord(2) + (m2-1) * 1./sup_cell(2)
          sup_coord(3) = new_coord(3) + (m3-1) * 1./sup_cell(3)
          !call Get_format(n, fmt_, '')
          call write_SUPERCELL_coord(n, nt, at_label(i), at_typ(i), at_occ_perc(i), sup_coord, '1+2+3')
          call Get_new_label(at_label(i), n, atom_label(nt))
          !atom_label(nt) = at_label(i)
          atom_typ(nt)   = at_typ(i)
          atom_coord(1:3, nt) = sup_coord(1:3)
          atom_Biso(nt)  = at_biso(i)
          atom_occ_perc(nt)   = at_occ_perc(i)
         end do
        end if
      end do
     end if ! sup_cell(2) > 1

     m3 = 1
     if(sup_cell(3) > 1) then
      do
       m3 = m3 + 1
       if(m3 > sup_cell(3)) exit
       n  = n  + 1
       nt = nt + 1
       if(nt > nb_atom_max) then
        if(write_details) then
        call write_info('')
        write(message_text, '(a,i4,a)') ' > Number of atoms in supercell is greater than allowed limit (', nb_atom_max,') !!'
        call write_info(trim(message_text))
        call write_info('')
        end if
        return
       end if
       sup_coord(1) = new_coord(1) + (m1-1) * 1./sup_cell(1)
       sup_coord(2) = new_coord(2)
       sup_coord(3) = new_coord(3) + (m3-1) * 1./sup_cell(3)
       !call Get_format(n, fmt_, '')
       call write_SUPERCELL_coord(n, nt, at_label(i), at_typ(i), at_occ_perc(i), sup_coord, '1+3')
       call Get_new_label(at_label(i), n, atom_label(nt))
       !atom_label(nt) = at_label(i)
       atom_typ(nt)   = at_typ(i)
       atom_coord(1:3, nt) = sup_coord(1:3)
       atom_Biso(nt)  = at_biso(i)
       atom_occ_perc(nt)   = at_occ_perc(i)
      end do
     end if
     end do


   else     ! sup_cell(1) > 1
    if(sup_cell(2) > 1) then
     m2 = 0
     do
      m2 = m2 + 1
      if(m2 > sup_cell(2)) exit
      n  = n  + 1
      nt = nt + 1
      if(nt > nb_atom_max) then
       if(write_details) then
       call write_info('')
       write(message_text, '(a,i4,a)') ' > Number of atoms in supercell is greater than allowed limit (', nb_atom_max,') !!'
       call write_info(trim(message_text))
       call write_info('')
       end if
       return
      end if
      sup_coord(1) = new_coord(1)
      sup_coord(2) = new_coord(2) + (m2-1) * 1./sup_cell(2)
      sup_coord(3) = new_coord(3)
      !call Get_format(n, fmt_, '')
      call write_SUPERCELL_coord(n, nt, at_label(i), at_typ(i), at_occ_perc(i), sup_coord, '2')
      call Get_new_label(at_label(i), n, atom_label(nt))
      !atom_label(nt) = at_label(i)
      atom_typ(nt)   = at_typ(i)
      atom_coord(1:3, nt) = sup_coord(1:3)
      atom_Biso(nt)  = at_biso(i)
      atom_occ_perc(nt)   = at_occ_perc(i)

      m3 = 1
      do
       m3 = m3 + 1
       if(m3 > sup_cell(3)) exit
       n  = n  + 1
       nt = nt + 1
       if(nt > nb_atom_max) then
        if(write_details) then
        call write_info('')
        write(message_text, '(a,i4,a)') ' > Number of atoms in supercell is greater than allowed limit (', nb_atom_max,') !!'
        call write_info(trim(message_text))
        call write_info('')
        end if
        return
       end if
       sup_coord(1) = new_coord(1)
       sup_coord(2) = new_coord(2) + (m2-1) * 1./sup_cell(2)
       sup_coord(3) = new_coord(3) + (m3-1) * 1./sup_cell(3)
       !call Get_format(n, fmt_, '')
       call write_SUPERCELL_coord(n, nt, atom_orbit%label(i), atom_orbit%typ(i), at_occ_perc(i), sup_coord, '2+3')
       call Get_new_label(at_label(i), n, atom_label(nt))
       !atom_label(nt) = at_label(i)
       atom_typ(nt)   = at_typ(i)
       atom_coord(1:3, nt) = sup_coord(1:3)
       atom_Biso(nt)  = at_biso(i)
       atom_occ_perc(nt)   = at_occ_perc(i)
      end do

     end do



    else   ! sup_cell(3) > 1
     m3 = 0
     do
      m3 = m3 + 1
      if(m3 > sup_cell(3)) exit
      n  = n  + 1
      nt = nt + 1
      if(nt > nb_atom_max) then
       if(write_details) then
       call write_info('')
       write(message_text, '(a,i4,a)') ' > Number of atoms in the the supercell is greater than allowed limit (', nb_atom_max,') !!'
       call write_info(trim(message_text))
       call write_info('')
       end if
       return
      end if
      sup_coord(1) = new_coord(1)
      sup_coord(2) = new_coord(2)
      sup_coord(3) = new_coord(3) + (m3-1) * 1./sup_cell(3)
      !call Get_format(n, fmt_, '')
      call write_SUPERCELL_coord(n, nt, at_label(i), at_typ(i), at_occ_perc(i), sup_coord, '3')
      call Get_new_label(at_label(i), n, atom_label(nt))
      !atom_label(nt) = at_label(i)
      atom_typ(nt)   = at_typ(i)
      atom_coord(1:3, nt) = sup_coord(1:3)
      atom_Biso(nt)  = at_biso(i)
      atom_occ_perc(nt)   = at_occ_perc(i)
     end do


    end if


   end if   ! sup_cell(1) > 1
  end do    ! i=1, n_at


   !new space_group
  call set_SpaceGroup("P 1",SPG)    !Space group is set to P 1  before structure factor calculations

   !Atom list objet
  nb_atom      = n_total
  atom_mult(1:nb_atom)     = 1.0
  !atom_occ_perc(1:nb_atom) = 1.0



 return
end subroutine Generate_atoms_in_supercell

!--------------------------------------------------------------------------------------------------------
subroutine write_SUPERCELL_coord(n, nt, atom_label, atom_type, atom_occ_perc, sup_coord, input_string)
 use CRYSCALC_module, only : SUPERCELL_pcr, nb_atom_orbit, message_text, write_details
 USE IO_module,       ONLY : write_info
 implicit none
  integer,            intent(in)   :: n, nt
  character (len=*),  intent(inout):: atom_label
  character (len=*),  intent(in)   :: atom_type
  real,               intent(in)   :: atom_occ_perc
  real, dimension(3), intent(in)   :: sup_coord
  character (len=*),  intent(in)   :: input_string
  character (len=16)               :: new_label
  character (len=48)               :: fmt_

  if(.not. write_details) return

  new_label = ''

  if(SUPERCELL_pcr) then
   !call Get_format(n, fmt_, 'pcr')
   !write(message_text, fmt=trim(fmt_))  trim(atom_label), '_', n, atom_type, sup_coord(1), &
   !                                                                          sup_coord(2), &
   !                                                                          sup_coord(3), '  0.25000   1.00000'
   !call write_info(trim(message_text))
   !if(nb_atom_orbit == 0) then
   !call write_info('                    0.00      0.00      0.00     0.00      0.00')
   !else
   !call write_info('                      0.00      0.00      0.00     0.00      0.00')
   !endif

   call Get_format_(n, atom_label, new_label, fmt_, 'pcr')
   !write(message_text, fmt=trim(fmt_))  new_label(1:8), atom_type, sup_coord(1), &
   !                                                                sup_coord(2), &
   !                                                                sup_coord(3), '  0.25000   1.00000'
   write(message_text, fmt=trim(fmt_))  new_label(1:12), atom_type, sup_coord(1), &
                                                                   sup_coord(2), &
                                                                   sup_coord(3), '  0.25000', atom_occ_perc
   call write_info(trim(message_text))
   call write_info('                           0.00      0.00      0.00     0.00      0.00')

  else

   call Get_format_(n, atom_label, new_label, fmt_, '')
   write(message_text, fmt=trim(fmt_)) nt, new_label(1:12), atom_type, sup_coord(1), &
                                                                      sup_coord(2), &
                                                                      sup_coord(3), '   ['//input_string//']'
   call write_info(trim(message_text))

   !call Get_format(n, fmt_, '')
   !write(message_text, fmt=trim(fmt_))  nt, trim(atom_label), '_', n, atom_type, sup_coord(1), &
   !                                                                              sup_coord(2), &
   !                                                                              sup_coord(3), '   ['//input_string//']'

   !call write_info(trim(message_text))

  end if


 return
end subroutine write_SUPERCELL_coord

!---------------------------------------------------------------------------

subroutine Get_format_(n, atom_label, new_label, fmt_, input_string)
 implicit none
 integer,            intent(in)    :: n
 character (len=*),  intent(in)    :: atom_label
 character (len=16), intent(inout) :: new_label
 character (len=48), intent(inout) :: fmt_
 character (len=*),  intent(in)    :: input_string


 call Get_new_label(atom_label, n, new_label)

 if(len_trim(input_string) == 3) then
  if(input_string(1:3) == 'pcr') then
   !fmt_ = '(a, 5x,a, 3F10.5,a)'
   fmt_ = '(a, 5x,a, 3F10.5,a, F10.5)'
  end if
 else
  fmt_ = '(5x,i4,2x,a, 5x,a, 3F10.5,a)'
 end if

 return
end subroutine Get_format_

!---------------------------------------------------------------------------

subroutine Get_format(n, fmt_, input_string)
 implicit none
 integer,            intent(in)    :: n
 character (len=48), intent(inout) :: fmt_
 character (len=*),  intent(in)    :: input_string

  if(len_trim(input_string) == 3) then
   if(input_string(1:3) == 'pcr') then
   if(n < 10) then
    fmt_ = '(2a,i1, 5x,a, 3F10.5,a)'
   elseif(n < 100) then
    fmt_ = '(2a,i2, 4x,a, 3F10.5,a)'
   elseif(n < 1000) then
    fmt_ = '(2a,i3, 3x,a, 3F10.5,a)'
   elseif(n < 10000) then
    fmt_ = '(2a,i4, 2x,a, 3F10.5,a)'
   endif
   end if
  else
   if(n < 10) then
    fmt_ = '(5x,i4,2x,2a,i1, 5x,a, 3F10.5,a)'
   elseif(n < 100) then
    fmt_ = '(5x,i4,2x,2a,i2, 4x,a, 3F10.5,a)'
   elseif(n < 1000) then
    fmt_ = '(5x,i4,2x,2a,i3, 3x,a, 3F10.5,a)'
   elseif(n < 10000) then
    fmt_ = '(5x,i4,2x,2a,i4, 2x,a, 3F10.5,a)'
   endif
  end if

 return
end subroutine    Get_format


!-----------------------------------------------------------------------------------
subroutine Get_new_label(atom_label, n, new_label)
 implicit none
  character (len=16), intent(in)     :: atom_label
  integer,           intent(in)     :: n
  character (len=16),intent(inout) :: new_label

  if(n < 10) then
   write(new_label, '(2a,i1)') trim(atom_label),'_',n
  elseif(n < 100) then
   write(new_label, '(2a,i2)') trim(atom_label),'_',n
  elseif(n < 1000) then
   write(new_label, '(2a,i3)') trim(atom_label),'_',n
  elseif(n < 10000) then
   write(new_label, '(2a,i4)') trim(atom_label),'_',n
  end if


   return
 end subroutine Get_new_label

!------------------------------------------------------------------------------------

subroutine write_cryscalc_cmd_log(cmd_line)
 use CRYSCALC_module, only : tmp_unit
 implicit none
  CHARACTER (LEN=*), intent(in)         :: cmd_line
  CHARACTER (LEN=10)                    :: date, time


  open(unit=tmp_unit, file="cryscalc_cmd.log", access="append")
   call date_and_time(date, time)
   !write(tmp_unit, '(2a)') "CRYSCALC ", trim(cmd_line)
   WRITE(tmp_unit, '(11a,10x,2a)') date(7:8),'-', date(5:6),'-', date(1:4), '  at ', &
                                   time(1:2),':', time(3:4),':', time(5:6), "CRYSCALC ", trim(cmd_line)
  close (unit=tmp_unit)


  return
end subroutine write_cryscalc_cmd_log
