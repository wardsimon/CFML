!     Last change:  TR   27 Nov 2006    2:44 pm
!------------------------------------------------------------------------

subroutine Write_HELP 
 USE cryscal_module, ONLY : nb_help, nb_help_max
 USE IO_module,      ONLY : write_info

 implicit none

  IF(nb_help == nb_help_max) call write_header 
  call HELP_on_line

  IF(nb_help == nb_help_max) then 
   call write_info(' >>> List of CRYSCAL commands line arguments :')
   call write_info('     ----------------------------------------')
   call write_info('')
   call write_cryscal_CLA 
   call write_info(' ')
   call write_info(' ')
   call write_info('   CRYSCAL manual stored in CRYSCAL_manual.txt file. ')
   call write_info(' ')
  endif
  stop

end subroutine Write_HELP

!------------------------------------------------------------------------
 subroutine Write_cryscal_HTML
  use cryscal_module,               ONLY : my_browser, browse_cryscal_HTML
  use IO_module,                    ONLY : write_info
  USE external_applications_module, ONLY : launch_browser

  !browse_cryscal_HTML = .false. 
  call create_CRYSCAL_HTML
  call write_info(' ')
  call write_info(' ')
  
  
 

  stop
  return
  
 end subroutine Write_cryscal_html
!------------------------------------------------------------------------

subroutine Write_KEYWORD()
  USE IO_module, only : write_info

  implicit none

  call KEYS_on_line()

  call write_info(' ')
  call write_info(' ')
  call write_info('   CRYSCAL keys stored in CRYSCAL_KEYS.txt file. ')
  call write_info(' ')
  stop

end subroutine Write_KEYWORD

!------------------------------------------------------------------------


subroutine write_wave_features()
 USE cryscal_module, ONLY: ON_SCREEN, wavelength, keyword_beam, message_text, wavelength, pi, neutrons, DEVICE
 USE IO_module,      ONLY: write_info
  implicit none
   REAL                  :: Ki
   REAL                  :: neutron_velocity
   REAL                  :: neutron_energy
   REAL                  :: neutron_temperature
   REAL                  :: X_energy

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
 USE cryscal_module, ONLY: ON_SCREEN, wavelength, message_text, DEVICE
 USE IO_module,      ONLY: write_info
  implicit none

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
 USE cryscal_module, ONLY: ON_SCREEN, beam_type
 USE IO_module,      ONLY: write_info
  implicit none
 
  if(.not. ON_SCREEN) return
  call write_info('')
  if(beam_type(1:8) == 'neutrons') then
   call write_info(' > Incident radiation beam : neutrons')
  else
   call write_info(' > Incident radiation beam : X-rays')
  endif
  
  call write_wave_features
 return
end subroutine write_BEAM_features
!-------------------------------------------------------------------------

subroutine write_Xrays_wavelength()
 USE cryscal_module,            ONLY : message_text
 USE IO_module,                 ONLY : write_info
 USE CFML_Scattering_Chemical_Tables
 USE wavelength_module
 implicit none
  INTEGER           :: i, n
  REAL              :: mean_Ka

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
  USE cryscal_module,          ONLY : message_text, qvec
  implicit none
    
   call write_info('')
   write(message_text, '(a,3F6.2)') '  > Modulation wave vector Qvec = ', qvec(1:3)
   call write_info(trim(message_text))
   call write_info('')
   
  return
 end subroutine Write_QVEC
!-------------------------------------------------------------------------

subroutine write_space_group(i1, i2)
 use CFML_crystallographic_symmetry, only : set_spacegroup
 USE cryscal_module,            ONLY : list_sg, list_sg_centric, list_sg_multip, message_text, SPG
 USE IO_module,                 ONLY : write_info
 implicit none
  INTEGER, INTENT(IN)       :: i1, i2   ! numero du groupe d'espace
  CHARACTER (LEN=3)         :: string_i
  INTEGER                   :: i, n
  LOGICAL                   :: ok

  ! numero du groupe d'espace: 		SPG%NumSpg
  ! symbole:                            SPG%SPG_Symb
  ! centro:                             SPG%Centred  =0 Centric(-1 no at origin)
  !                                                  =1 Acentric
  !                                                  =2 Centric(-1 at origin)
  n=0
  do i = i1, i2
   WRITE(string_i, '(i3)') i
   call set_spacegroup(string_i, SPG)

   IF(list_sg_centric(1) .AND. SPG%centred==1) cycle
   IF(list_sg_centric(2) .AND. SPG%centred/=1) cycle

   call test_Bravais(SPG%SPG_symb, ok)
   IF(.NOT. ok) cycle

   call test_laue(SPG%Laue, ok)
   IF(.NOT. ok) cycle

   IF(list_sg_multip) then
    WRITE(message_text, '(10x,a,I3,a,5x,a,10x,3a,I3)') 'IT# ', i,'.', SPG%SPG_symb(1:10) , '(',TRIM(SPG%laue),')    mult.=', &
	                                                    SPG%multip
   else
    WRITE(message_text, '(10x,a,I3,a,5x,a,10x,3a)')    'IT# ', i,'.', SPG%SPG_symb(1:10) , '(',TRIM(SPG%laue),')'
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
 USE cryscal_module,                 ONLY : SPG, message_text
 USE CFML_symmetry_tables,           ONLY : intsymoh,x_oh, intsymd6h,x_d6h

 implicit none
  !INTEGER, INTENT(IN)             :: SG_numor
  CHARACTER (LEN=*), INTENT(IN)   :: SG_symbol
  CHARACTER(LEN=3)                :: string_numor
  INTEGER                         :: i, j, i1, i2
  INTEGER, DIMENSION(192)         :: indx_op

  !WRITE(string_numor, '(i3)') SG_numor
  !call set_spacegroup(string_numor, SPG)

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
 USE cryscal_module, ONLY : message_text, nb_atom, keyword_SPGR,                           &
                            atom_label, atom_type, atom_coord, atom_occ, atom_occ_perc,    &
                            atom_mult, atom_Biso,  atom_occ, atom_Ueq, SPG, input_PCR,     &
                            keyword_create_CIF, keyword_read_PCR
 use CFML_Crystallographic_Symmetry,        only : Get_Multip_Pos

 USE IO_module,      ONLY : write_info

 implicit none
  INTEGER              :: i
  REAL                 :: occ_perc  ! % occ. site
  REAL                 :: f1, f

 call write_info('')
 call write_info('  > ATOMS LIST:')
 call write_info('')


 IF(keyword_SPGR) then
  if(input_PCR .or. keyword_read_PCR) then
   call write_info('   NAME     LABEL    x         y         z         Biso      Occ       Occ(%)    Mult')
  else
   call write_info('   NAME     LABEL    x         y         z         Biso      Occ(%)    Mult')
  endif
 else
  call write_info( '   NAME     LABEL    x         y         z         Biso      Occ')
 endif
 call write_info('')

 IF(keyword_SPGR .and. (input_PCR .or. keyword_read_PCR)) then
  atom_mult(1) = Get_Multip_Pos(atom_coord(1:3,1), SPG)
  f1 = atom_occ(1) * REAL(SPG%Multip) / REAL(atom_mult(1))
  f1 = 1./f1
 endif

 do i= 1, nb_atom
  IF(keyword_SPGR) then
   atom_mult(i) = Get_Multip_Pos(atom_coord(1:3,i), SPG)
   if(input_PCR .or. keyword_read_PCR) then
    f =  SPG%Multip / atom_mult(i)
    occ_perc = atom_occ(i) *  f * f1
    atom_occ_perc(i) =  occ_perc
    
    write(message_text,'(3x,a4,5x,a4,2x,6F10.5,I4)')  atom_label(i), atom_type (i) , atom_coord(1:3,i), &
                                                      atom_Biso(i),  atom_occ(i),    occ_perc, atom_mult(i)
   else
    write(message_text,'(3x,a4,5x,a4,2x,5F10.5,I4)')  atom_label(i), atom_type (i) , atom_coord(1:3,i), &
                                                      atom_Biso(i),  atom_occ_perc(i),    atom_mult(i)

   endif
  else
   
   write(message_text,'(3x,a4,5x,a4,2x,5F10.5)')     atom_label(i), atom_type (i) , atom_coord(1:3,i), &
                                                     atom_Biso(i),  atom_occ_perc(i)

  endif
  call write_info(TRIM(message_text))
 end do

 IF(keyword_create_CIF) then
   call write_CIF_file('ATOMS_HEADER')
   call write_CIF_file('ATOM')
 endif

end subroutine write_atom_list


!-----------------------------------------------------------------------------------------------------

subroutine write_molecular_features
 use cryscal_module, only : molecule, message_text
 use IO_module,      only : write_info
 implicit none

 call write_info('')
 call write_info('  >> Molecular features:')
 call write_info('')
 write(message_text, '(2a)')       '    . formula : ', trim(molecule%formula)
 call write_info(trim(message_text))
 write(message_text, '(a, F8.2)')  '    . weight  : ', molecule%weight
 call write_info(trim(message_text))
 if(molecule%density > 0.01) then
  write(message_text, '(a, F9.3)') '    . density : ', molecule%density
  call write_info(trim(message_text))
 endif
 IF(LEN_TRIM(molecule%content) /=0) then
  write(message_text, '(2a)')      '    . content : ', trim(molecule%content)
  call write_info(trim(message_text))
 endif
 call write_info('')

 return
end subroutine write_molecular_features


!-----------------------------------------------------------------------------------------------------
subroutine write_REF(input_string)
 USE cryscal_module, ONLY: keyword_create_CIF, CIF_parameter_KCCD, CIF_parameter_APEX, CIF_parameter_XCALIBUR, EVAL, SADABS
 USE IO_module,      ONLY: write_info
 implicit none
 CHARACTER(LEN=*), INTENT(IN) :: input_string


 IF(input_string(1:4) == 'KCCD') then

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
   call write_info("_diffrn_detector                        "//trim(CIF_parameter_KCCD%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean        "//trim(CIF_parameter_KCCD%diffrn_detector_area_resol_mean))
   call write_info("_diffrn_measurement_device              "//trim(CIF_parameter_KCCD%diffrn_measurement_device))
   call write_info("_diffrn_measurement_device_type         "//trim(CIF_parameter_KCCD%diffrn_measurement_device_type))
   call write_info("_diffrn_measurement_method              "//trim(CIF_parameter_KCCD%diffrn_measurement_method))
  endif


 elseIF(input_string(1:4) == 'APEX') then
  IF(keyword_create_CIF)   then
   call write_CIF_file('APEX')
  else
   call write_info("")
   !call write_info("#----------------------------------------------------------------------------#")
   !call write_info("#                   COMPUTER PROGRAMS USED                                   #")
   !call write_info("#----------------------------------------------------------------------------#")
   !call write_info("")
   !call write_info("_computing_data_collection       "//trim(CIF_parameter_APEX%computing_data_collection))
   !call write_info("_computing_cell_refinement       "//trim(CIF_parameter_APEX%computing_cell_refinement))
   !call write_info("_computing_data_reduction        "//trim(CIF_parameter_APEX%computing_data_reduction))
   !call write_info("_computing_structure_solution    "//trim(CIF_parameter_APEX%computing_structure_solution))
   !call write_info("_computing_structure_refinement  "//trim(CIF_parameter_APEX%computing_structure_refinement))
   !call write_info("_computing_molecular_graphics    "//trim(CIF_parameter_APEX%computing_molecular_graphics))
   !call write_info("_computing_publication_material  "//trim(CIF_parameter_APEX%computing_publication_material))
   !call write_info("")

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
   call write_info("_diffrn_detector                  "//trim(CIF_parameter_APEX%diffrn_detector))
   call write_info("_diffrn_detector_area_resol_mean  "//trim(CIF_parameter_APEX%diffrn_detector_area_resol_mean))
  endif


 ELSEIF(input_string(1:4) == 'EVAL') then

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

 ELSEIF(input_string(1:5) == 'DENZO') then

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

 ELSEIF(input_string(1:6) == 'SADABS') then
  if(keyword_create_CIF) then
   call write_CIF_file('SADABS')
  else
   call write_info("")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("#                   ABSORPTION CORRECTION                                    #")
   call write_info("#----------------------------------------------------------------------------#")
   call write_info("")
   call write_info(trim(SADABS%type))
   call write_info(trim(SADABS%details(1)))
   call write_info(trim(SADABS%details(2)))
   call write_info(trim(SADABS%details(3)))
   call write_info(trim(SADABS%details(4)))
  endif
 endif


 return
END subroutine write_REF

