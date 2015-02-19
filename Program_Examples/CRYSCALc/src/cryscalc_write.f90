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
  
 end subroutine Write_cryscalc_html
!------------------------------------------------------------------------

subroutine Write_version
 use cryscalc_module,      ONLY : cryscalc, CFML_version
 use IO_module,            ONLY : write_info
 
 call write_info('')
 call write_info('   > CRYSCALC version: '//trim(cryscalc%version))
 call write_info('   > Author:           '//trim(cryscalc%author))
 call write_info('   > CFML version:     '//trim(CFML_version))

 call write_info('')

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

end subroutine Write_KEYWORD

!------------------------------------------------------------------------


subroutine write_wave_features()
 USE cryscalc_module, ONLY: ON_SCREEN, wavelength, keyword_beam, message_text, wavelength, pi, neutrons, DEVICE, debug_proc
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

   if(.not. ON_SCREEN) return
   
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
    call write_info('  > Laboratory: '// trim(DEVICE%lab))
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
  INTEGER           :: i, n
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

subroutine write_space_group(i1, i2, n_enantio, n_chiral, n_polar)
 use CFML_crystallographic_symmetry, only : set_spacegroup
 USE cryscalc_module,           ONLY : list_sg, list_sg_centric, list_sg_multip, list_sg_enantio, list_sg_chiral, &
                                       list_sg_polar, message_text, SPG, debug_proc
 USE IO_module,                 ONLY : write_info
 implicit none
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
  
  ! numero du groupe d'espace: 		    SPG%NumSpg
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

   call test_polar(i, polar)
   if(polar) n_polar = n_polar + 1
   if(.not. polar .and. list_sg_polar) cycle

   write(Laue_string, '(3a)') '(', trim(SPG%laue), ')'	   
   write(PG_string,   '(3a)') '[', trim(SPG%PG), ']'	
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

   IF(list_sg_multip) then
	WRITE(message_text, '(10x,a,I3,a,5x,4a,10x,a,5x,2a,I3)') 'IT# ', i,'.', chiral_string(1:1),enantio_string(1:1),     &
	                                                         polar_string(1:1), SPG%SPG_symb(1:10) , Laue_string(1:10), &
															 PG_string(1:8), '    mult.=', SPG%multip 	
	
    
   else
    WRITE(message_text, '(10x,a,I3,a,5x,4a,10x,a,5x,a)')     'IT# ', i,'.', chiral_string(1:1),enantio_string(1:1),     &
	                                                         polar_string(1:1), SPG%SPG_symb(1:10) , Laue_string(1:10), &
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
 USE IO_module,                      ONLY : write_info
 USE CFML_crystallographic_symmetry, ONLY : set_spacegroup, searchop, write_sym
 USE cryscalc_module,                ONLY : SPG, message_text, debug_proc
 USE CFML_symmetry_tables,           ONLY : intsymoh,x_oh, intsymd6h,x_d6h

 implicit none
  !INTEGER, INTENT(IN)             :: SG_numor
  CHARACTER (LEN=*), INTENT(IN)   :: SG_symbol
  CHARACTER(LEN=3)                :: string_numor
  INTEGER                         :: i, j, i1, i2
  INTEGER, DIMENSION(192)         :: indx_op

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


 return
end subroutine write_current_space_group

!---------------------------------------------------------------------------------------------

subroutine write_atom_list()
 USE cryscalc_module, ONLY : message_text, nb_atom, keyword_SPGR, keyword_ADP_list,         &
                             atom_label, atom_typ, atom_coord, atom_occ, atom_occ_perc,     &
                             atom_mult, atom_Biso,  atom_occ, atom_Ueq, SPG, input_PCR,     &
                             keyword_create_CIF, write_atoms_in_A, keyword_read_PCR,        &
							 atom_adp_aniso, unit_cell, keyword_CELL, debug_proc
 use CFML_Crystallographic_Symmetry,    ONLY : Get_Multip_Pos
  USE IO_module,                         ONLY : write_info

 implicit none
  INTEGER              :: i
  REAL                 :: occ_perc  ! % occ. site
  REAL                 :: f1, f
 
 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_atoms_list")
 
 call write_info('')
 call write_info('  > ATOMS LIST:')
 call write_info('')

 if(keyword_ADP_list) then
    call write_info('   NAME     LABEL    U_11      U_22      U_33      U_23      U_13      U_12')
 else
  IF(keyword_SPGR) then
   if(input_PCR .or. keyword_read_PCR) then
    if(.not. write_atoms_in_A) then
     call write_info('         NAME     LABEL    x         y         z         Biso      Occ       Occ(%)    Mult')
	else
	 call write_info('         NAME     LABEL    x         y         z         Biso      Occ       Occ(%)    Mult' // &
                     '         x(A)      y(A)      z(A)')
    endif	
   else
    if(.not. write_atoms_in_A) then   
     call write_info('         NAME     LABEL    x         y         z         Biso      Occ       Occ(%)    Mult')
	else
	 call write_info('         NAME     LABEL    x         y         z         Biso      Occ       Occ(%)    Mult' // &
                     '         x(A)      y(A)      z(A)')	
    endif	
   endif
  else
   if(.not. write_atoms_in_A) then
    call write_info( '         NAME     LABEL    x         y         z         Beq       Occ')
   else
    call write_info( '         NAME     LABEL    x         y         z         Beq       Occ         x(A)      y(A)      z(A)')
   endif   
  endif
 endif 
 call write_info('')

 
 IF(keyword_SPGR .and. (input_PCR .or. keyword_read_PCR)) then
  atom_mult(1) = Get_Multip_Pos(atom_coord(1:3,1), SPG)
  f1 = atom_occ(1) * REAL(SPG%Multip) / REAL(atom_mult(1))
  f1 = 1./f1
 endif

 do i= 1, nb_atom
  
 if(keyword_ADP_list) then
   write(message_text,'(3x,a4,5x,a4,2x,6F10.5)')     atom_label(i), atom_typ (i) , atom_adp_aniso(1:3,i), &
                                                     atom_adp_aniso(6,i), atom_adp_aniso(5,i), atom_adp_aniso(4,i)

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
 if(write_atoms_in_A .and. keyword_CELL .and. keyword_read_PCR) then
  write(message_text, '(97x,a)') '0.00000   0.00000   0.00000'
  call write_info(message_text)
 end if 
 
 end do
  

 IF(keyword_create_CIF) then
   call write_CIF_file('ATOMS_HEADER')
   call write_CIF_file('ATOMS_HEADER_LOOP')   
   call write_CIF_file('ATOM')
 endif

end subroutine write_atom_list

!-----------------------------------------------------------------------------------------------------
subroutine write_atom_list_cart
 use io_module
 USE cryscalc_module,       ONLY : message_text, nb_atom, atom_label, atom_typ, atom_coord, crystal_cell, &
                                   debug_proc,                                                            &       
                                   cartesian_frame, create_SHAPE_file, tmp_unit, main_title
 USE CFML_crystal_metrics, only : cart_vector
 USE CFML_GlobalDeps,      ONLY : sp
 implicit none
  integer         :: i
  real (kind=sp), dimension(3)   :: cart_vect
  character (len=12)             :: string
 
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
 
  
 do i= 1, nb_atom
  cart_vect(1:3) = cart_vector('D', atom_coord(1:3, i), crystal_cell)
  
  write(message_text,'(3x,i3,3x,a4,5x,a4,2x,3F12.6)')     i, atom_label(i), atom_typ(i) , cart_vect(1:3)
  call write_info(TRIM(message_text))
  
  if(create_SHAPE_file) write(unit = tmp_unit, fmt='(a4,2x,3F12.6)')     atom_label(i) , cart_vect(1:3)
 end do					

  if(create_SHAPE_file) then
   close(unit=tmp_unit)
   call write_info('')
   call write_info('   > cryscalc_shape.dat for SHAPE has been created.')
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
 USE cryscalc_module, ONLY : keyword_create_CIF, EVAL, SADABS, ABS_CRYSALIS, debug_proc
 USE CIF_module,      ONLY : CIF_parameter_KCCD, CIF_parameter_APEX, CIF_parameter_X2S, CIF_parameter_XCALIBUR, &
                             CIF_parameter_SUPERNOVA, CIF_parameter_D8_VENTURE_Cu, CIF_parameter_D8_VENTURE_Mo
 USE IO_module,       ONLY: write_info
 implicit none
 CHARACTER(LEN=*), INTENT(IN) :: input_string
 INTEGER                      :: long_input_string, i
 LOGICAL                      :: input_KCCD, input_APEX, input_X2S, input_XCALIBUR, input_SUPERNOVA
 LOGICAL                      :: input_EVAL, input_DENZO, input_SADABS, input_ABS_CRYSALIS
 LOGICAL                      :: input_D8_VENTURE_Cu, input_D8_VENTURE_Mo

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_ref ("//trim(input_string)//")")
 
 long_input_string = len_trim(input_string)
 input_KCCD          = .false.
 input_APEX          = .false.
 input_X2S           = .false.
 input_XCALIBUR      = .false.
 input_SUPERNOVA     = .false. 
 input_EVAL          = .false.
 input_DENZO         = .false.
 input_SADABS        = .false.
 input_ABS_CRYSALIS  = .false.
 input_D8_VENTURE_Cu = .false.
 input_D8_VENTURE_Mo = .false.

 if(long_input_string == 3) then 
  if(input_string(1:3)  == 'X2S')           input_X2S = .true.
 elseif(long_input_string == 4) then
  if(input_string(1:4)  == 'KCCD')          input_KCCD = .true.
  if(input_string(1:4)  == 'APEX')          input_APEX = .true.
  if(input_string(1:4)  == 'EVAL')          input_EVAL = .true.
 elseif(long_input_string == 5) then 
  if(input_string(1:5)  == 'DENZO')         input_DENZO = .true.
 elseif(long_input_string == 6) then
  if(input_string(1:6)  == 'SADABS')        input_SADABS = .true.
 elseif(long_input_string == 8)  then
  if(input_string(1:8)  == 'XCALIBUR')      input_XCALIBUR = .true.
 elseif(long_input_string == 9)  then
  if(input_string(1:9)  == 'SUPERNOVA')     input_SUPERNOVA    = .true.
 elseif(long_input_string == 12)  then
  if(input_string(1:12) == 'ABS_CRYSALIS')  input_ABS_CRYSALIS = .true.
 elseif(long_input_string == 13) then
  if(input_string(1:15) == 'D8_VENTURE_CU') input_D8_VENTURE_Cu = .true.
  if(input_string(1:15) == 'D8_VENTURE_MO') input_D8_VENTURE_Mo = .true.
 endif 
  
  
  
 IF(input_KCCD) then

  IF(keyword_create_CIF) then
   call write_CIF_file("KCCD")
  else
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_source                          "//trim(CIF_parameter_KCCD%diffrn_source))
   call write_info("_diffrn_radiation_wavelength            "//trim(CIF_parameter_KCCD%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type                  "//trim(CIF_parameter_KCCD%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source                "//trim(CIF_parameter_KCCD%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator         "//trim(CIF_parameter_KCCD%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe                 "//trim(CIF_parameter_KCCD%diffrn_radiation_probe))
   call write_info("_diffrn_radiation_source                "//trim(CIF_parameter_KCCD%diffrn_radiation_source))
   call write_info("_diffrn_detector                        "//trim(CIF_parameter_KCCD%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean        "//trim(CIF_parameter_KCCD%diffrn_detector_area_resol_mean))
   call write_info("_diffrn_measurement_device              "//trim(CIF_parameter_KCCD%diffrn_measurement_device))
   call write_info("_diffrn_measurement_device_type         "//trim(CIF_parameter_KCCD%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method              "//trim(CIF_parameter_KCCD%diffrn_measurement_method))
  endif


 ELSEIF(input_APEX) then
  IF(keyword_create_CIF)   then
   call write_CIF_file('APEX')
  else
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_APEX%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_APEX%diffrn_measurement_method))
   call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_APEX%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_APEX%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_APEX%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_APEX%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_APEX%diffrn_radiation_probe))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_APEX%diffrn_radiation_source))
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_APEX%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_APEX%diffrn_detector_area_resol_mean))
  endif

 ELSEIF(input_X2S) then
  IF(keyword_create_CIF)   then
   call write_CIF_file('X2S')
  else
   call write_info("")   
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_X2S%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_X2S%diffrn_measurement_method))
   call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_X2S%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_X2S%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_X2S%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_X2S%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_X2S%diffrn_radiation_probe))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_X2S%diffrn_radiation_source))
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_X2S%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_X2S%diffrn_detector_area_resol_mean))
  endif

 ELSEIF(input_D8_VENTURE_CU) then
  IF(keyword_create_CIF)   then
   call write_CIF_file('D8_VENTURE_CU')
  else
   call write_info("")   
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_measurement_method))
   call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_radiation_probe))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_radiation_source))
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_D8_VENTURE_CU%diffrn_detector_area_resol_mean))
  endif

 ELSEIF(input_D8_VENTURE_MO) then
  IF(keyword_create_CIF)   then
   call write_CIF_file('D8_VENTURE_MO')
  else
   call write_info("")   
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_measurement_method))
   call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_radiation_probe))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_radiation_source))
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_D8_VENTURE_MO%diffrn_detector_area_resol_mean))
  endif

  
 ELSEIF(input_XCALIBUR) then
  IF(keyword_create_CIF)   then
   call write_CIF_file('XCALIBUR')
  else
   call write_info("")   
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_XCALIBUR%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_XCALIBUR%diffrn_measurement_method))
   call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_XCALIBUR%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_XCALIBUR%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_XCALIBUR%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_XCALIBUR%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_XCALIBUR%diffrn_radiation_probe))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_XCALIBUR%diffrn_radiation_source))
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_XCALIBUR%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_XCALIBUR%diffrn_detector_area_resol_mean))
  endif
 
 ELSEIF(input_SUPERNOVA) then
  IF(keyword_create_CIF)   then
   call write_CIF_file('SUPERNOVA')
  else
   call write_info("")   
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   DATA COLLECTION                                          #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_diffrn_measurement_device_type   "//trim(CIF_parameter_SUPERNOVA%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method        "//trim(CIF_parameter_SUPERNOVA%diffrn_measurement_method))
   call write_info("_diffrn_radiation_wavelength      "//trim(CIF_parameter_SUPERNOVA%diffrn_radiation_wavelength))
   call write_info("_diffrn_radiation_type            "//trim(CIF_parameter_SUPERNOVA%diffrn_radiation_type))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_SUPERNOVA%diffrn_radiation_source))
   call write_info("_diffrn_radiation_monochromator   "//trim(CIF_parameter_SUPERNOVA%diffrn_radiation_monochromator))
   call write_info("_diffrn_radiation_probe           "//trim(CIF_parameter_SUPERNOVA%diffrn_radiation_probe))
   call write_info("_diffrn_radiation_source          "//trim(CIF_parameter_SUPERNOVA%diffrn_radiation_source))
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_SUPERNOVA%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_SUPERNOVA%diffrn_detector_area_resol_mean))
  endif

  
 ELSEIF(input_EVAL) then

  IF(keyword_create_CIF) then
   call write_CIF_file('EVAL_PROGRAMS')
  else
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   COMPUTER PROGRAMS USED                                   #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   !call write_info("_computing_data_collection      "//trim(CIF_parameter_KCCD%computing_data_collection))
   !call write_info("_computing_cell_refinement      "//trim(CIF_parameter_KCCD%computing_cell_refinement))
   !call write_info("_computing_data_reduction       "//trim(CIF_parameter_KCCD%computing_data_reduction))
   
   call write_info("_computing_data_collection      "//trim(EVAL%data_collection))
   call write_info("_computing_cell_refinement      "//trim(EVAL%cell_refinement))
   call write_info("_computing_data_reduction       "//trim(EVAL%data_collection))
   
   
   call write_info("")
  endif

 ELSEIF(input_DENZO) then

  IF(keyword_create_CIF)  then
   call write_CIF_file('DENZO_PROGRAMS')
  else
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   COMPUTER PROGRAMS USED                                   #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info("_computing_data_collection      "//trim(CIF_parameter_KCCD%computing_data_collection))
   call write_info("_computing_cell_refinement      "//trim(CIF_parameter_KCCD%computing_cell_refinement))
   call write_info("_computing_data_reduction       "//trim(CIF_parameter_KCCD%computing_data_reduction))
   call write_info("")
  endif

 ELSEIF(input_SADABS) then
  if(keyword_create_CIF) then
   call write_CIF_file('SADABS')
  else
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   ABSORPTION CORRECTION                                    #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info(trim(SADABS%type))
   do i=1, 4
    call write_info(trim(SADABS%details(i)))
   end do	
   
  endif
  
  ELSEIF(input_ABS_CRYSALIS) then
  if(keyword_create_CIF) then
   call write_CIF_file('ABS_CRYSALIS')
  else
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   ABSORPTION CORRECTION                                    #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info(trim(ABS_CRYSALIS%type))
   do i=1, 6
    call write_info(trim(ABS_CRYSALIS%details(i)))
   end do	   
  endif
 endif


 return
END subroutine write_REF

!------------------------------------------------------------------------------------------------------------------------------------------
subroutine write_crystal_cell_Cart
 use cryscalc_module,      only : tmp_unit, crystal_cell
 use CFML_crystal_metrics, only : write_crystal_cell
 USE IO_module,            ONLY: write_info

 implicit none
 integer                        :: i_error, i
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