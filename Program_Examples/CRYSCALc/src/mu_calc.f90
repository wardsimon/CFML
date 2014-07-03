!     Last change:  TR   15 Jun 2007   11:00 am
subroutine absorption_calculation
 USE Neutrons_data
 USE Xrays_data
 USE cryscalc_module, ONLY : ON_SCREEN, neutrons, X_rays, wavelength,  message_text, keyword_BEAM,   &
                             known_data_neutrons, debug_proc
 USE wavelength_module
 USE IO_module ,      ONLY : write_info

 implicit none
  INTEGER             :: k

  if(debug_proc%level_2)  call write_debug_proc_level(2, "absorption_calculation")
  
 if(ON_SCREEN) then
 call write_info( '                                                    ')
 call write_info( '         Absorption coefficient calculation  (mu)   ')
 call write_info( '         ----------------------------------------   ')
 call write_info( '              mu = Si Ni.siT                        ')
 call write_info( '         with:   Ni: number of atoms i/cm3          ')
 call write_info( '                siT: Total interaction cross section')
 call write_info( '                                                    ')
 endif


 IF(keyword_BEAM) then
  if (neutrons) then
   if(ON_SCREEN) then
   call write_info('')
   call write_info('  > Incident beam radiation: neutrons')
   WRITE(message_text,'(1x,a,F10.6)')        '            Wavelength (A): ', wavelength
   call write_info(TRIM(message_text))
   endif

   !call F000_calculation("neutrons")
   if(.not. known_data_neutrons) then
    call definition_data_neutrons()     ! neutrons_mod.F90
    known_data_neutrons = .true.
   endif
   call neutron_cross_section()        !
   call neutron_absorption()
   call write_mu('yes')

  ELSEIF(X_rays) then
   call F000_calculation("X")

   if(ON_SCREEN) then
   call write_info('')
   call write_info('  > Incident beam radiation: X_rays')
   endif

   IF(.not. anti_cathode) then
    write(message_text, '(8x, a,F10.5,a)') '!! Xrays cross section are not tabulated for the current wavelength (', &
	                                       wavelength,  ' A)'
    call write_info(trim(message_text))
    write(message_text, '(8x,a)')          '   Molybdenum wavelength values will be used !!'
    call write_info(trim(message_text))
    call write_info('')
    !if(.not. anti_cathode) then
     X_target(2)%logic = .true.  ! Mo
     anti_cathode = .true.
    !endif
   endif

   do k=1, tabulated_target_nb
    if(anti_cathode) then
     if(X_target(k)%logic) then
      call Calcul_X_mu(X_target(k)%label)
      exit
     endif
    else
     call Calcul_X_mu(X_target(k)%label)
     if(ON_SCREEN) call write_info('')
    endif
   end do

  end if

 else
   call F000_calculation("X")
   if(ON_SCREEN) then
   call write_info('')
   call write_info('  > Incident beam radiation: X_rays')
   endif

   do k=1, tabulated_target_nb
    if(anti_cathode) then
     if(X_target(k)%logic) then
      call Calcul_X_mu(X_target(k)%label)
      exit
     endif
    else
     call Calcul_X_mu(X_target(k)%label)
     if(ON_SCREEN) call write_info('')
    endif
   end do


 ENDIF ! keyword_BEAM


 RETURN
end subroutine absorption_calculation
!-------------------------------------------------------------------

subroutine calcul_X_mu(input_string)
 USE Neutrons_data
 USE Xrays_data
 USE cryscalc_module, ONLY : ON_SCREEN, neutrons, X_rays, wavelength,  message_text, keyword_BEAM, &
                             known_data_x, debug_proc
 USE wavelength_module
 USE IO_module ,      ONLY : write_info
 implicit none
  CHARACTER (LEN=*), INTENT(IN)   :: input_string

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_X_mu ("//trim(input_string)//")")
  
  if(ON_SCREEN) then
   IF(input_string(1:2) == 'Ag') then
    WRITE(message_text, '(a,F10.5,a)')  '      ----  Wavelength Ka_Ag = ', X_target(1)%wave(1), ' A '
   elseIF(input_string(1:2) == 'Mo') then
    WRITE(message_text, '(a,F10.5,a)')  '      ----  Wavelength Ka_Mo = ', X_target(2)%wave(1), ' A '
   elseIF(input_string(1:2) == 'Cu') then
    WRITE(message_text, '(a,F10.5,a)')  '      ----  Wavelength Ka_Cu = ', X_target(3)%wave(1), ' A '
   elseIF(input_string(1:2) == 'Ni') then
    WRITE(message_text, '(a,F10.5,a)')  '      ----  Wavelength Ka_Ni = ', X_target(4)%wave(1), ' A '
   elseIF(input_string(1:2) == 'Co') then
    WRITE(message_text, '(a,F10.5,a)')  '      ----  Wavelength Ka_Co = ', X_target(5)%wave(1), ' A '
   elseIF(input_string(1:2) == 'Fe') then
    WRITE(message_text, '(a,F10.5,a)')  '      ----  Wavelength Ka_Fe = ', X_target(6)%wave(1), ' A '
   elseIF(input_string(1:2) == 'Cr') then
    WRITE(message_text, '(a,F10.5,a)')  '      ----  Wavelength Ka_Cr = ', X_target(7)%wave(1), ' A '
   endif
   call write_info(TRIM(message_text))
  endif

  !call F000_calculation("X")
  if (.not. known_data_x) then
   call definition_data_X() ! xrays_mod.F90
   known_data_X = .true.
  endif
  call X_attenuation_calculation()
  call write_mu('yes')


 RETURN
end subroutine calcul_X_mu
!-------------------------------------------------------------------
subroutine write_mu(input_string)
 USE cryscalc_module, ONLY : ON_SCREEN, message_text, keyword_SIZE, crystal, X_rays, neutrons,      &
                             keyword_TRANSMISSION, nb_dim_transmission, dim_transmission, keyword_create_CIF, &
                             keyword_WRITE_REF_APEX, keyword_WRITE_REF_X2S, absorption, SADABS_ratio,         &
					 		SADABS_Tmin, SADABS_Tmax, debug_proc
 USE IO_module,       ONLY : write_info
 USE math_module,     ONLY : transmission

 implicit none
  CHARACTER (LEN=*), INTENT(IN)   :: input_string
  REAL                            :: d_min, d_max, T_min, T_max, T
  INTEGER                         :: i
  REAL, parameter                 :: EPS=0.001
  
 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_mu ("//trim(input_string)//")")
 
 if(ON_SCREEN) then
  call write_info(' ')
  WRITE(message_text,'(a,F10.5,a)') '   > mu:   ',absorption%mu,     ' cm-1'
   call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F10.5,a)') '           ',0.1*absorption%mu, ' mm-1'
   call write_info(TRIM(message_text))
 endif

 IF(keyword_create_CIF)  call write_CIF_file('ABSORPTION')

 if (keyword_SIZE) then
  if(ON_SCREEN .and. crystal%radius < eps) then
   call write_info(' ')
   WRITE(message_text,'(a,F10.5,a)')  '   > R   : ', crystal%radius, ' mm'
   WRITE(message_text,'(a,F10.5)')    '   > mu*R: ',0.1*absorption%mu*  crystal%radius
   call write_info(TRIM(message_text))
   call write_info(' ')
  endif

  d_min = MINVAL(crystal%size)
  d_max = MAXVAL(crystal%size)

  absorption%Tmin = transmission(0.1*absorption%mu, d_max)
  absorption%Tmax = transmission(0.1*absorption%mu, d_min)

  if(ON_SCREEN) then
   call write_info('')
   WRITE(message_text, '(a,F7.4,a,F10.4)') '    . Min of transmission (for x =  ',d_max, ' mn) = ', absorption%Tmin
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F7.4,a,F10.4)') '    . Max of transmission (for x =  ',d_min, ' mn) = ', absorption%Tmax
   call write_info(TRIM(message_text))
  end if

  IF(keyword_create_CIF) then
   IF(keyword_WRITE_REF_APEX)  call write_CIF_file('SADABS')
   IF(keyword_WRITE_REF_X2S)   call write_CIF_file('SADABS')
   call Get_SADABS_ratio()

   call write_CIF_file('TMIN')
   call write_CIF_file('TMAX')
   if(ON_SCREEN) then
    call write_info("")
    call write_info("   Tmax and Tmin values correspond to EXPECTED values calculated from crystal size." )
    if(SADABS_ratio > 0.) then
     call write_info("   If SADABS program has been used for absorption correction, it provided one ")
     call write_info("   'relative-correction-factor'. In such a case, Tmax should be given as")
     call write_info("   Tmax_expected and Tmin = Tmax * 'relative-correction-factor'")
     write(message_text, '(3x,a,F8.4,a, F8.4)') ' Tmin = Tmax * ', sadabs_ratio, ' = ', absorption%Tmax * sadabs_ratio
     call write_info(trim(message_text))
    elseif(SADABS_Tmin > 0. .and. SADABS_Tmax > 0.) then
     write(message_text, '(3x,a,F8.4,a, F8.4)') ' Tmin = Tmax * ', sadabs_ratio, ' = ', absorption%Tmax * sadabs_ratio
     call write_info(trim(message_text))
    endif
    call write_info("")
   endif

  endif

 end if

 IF(keyword_TRANSMISSION .AND. nb_dim_transmission /=0) then
  if(ON_SCREEN) then
   call write_info('')
   call write_info('      dimension (mm)             transmission coefficient')
   do i=1, nb_dim_transmission
    T = transmission(0.1*absorption%mu, dim_transmission(i))
    WRITE(message_text, '(10x,F7.4,20x,F10.4)') dim_transmission(i), T
   call write_info(TRIM(message_text))
   end do
   call write_info('')
  endif
 endif

 IF(input_string(1:3) == 'yes' .and. ON_SCREEN) then
  if (X_rays) then
   call write_info(' ')
   call write_info('  Xrays total interaction cross sections values (in barns) are extracted from: ')
   call write_info('    Crystallographic Internationales Tables vol.C  1995, p.193-199 ')
   call write_info(' ')
  elseif(neutrons) then
   call write_info(' ')
   call write_info('  Neutron cross-sections values (in barns) are extracted from:')
   call write_info('    V.F. Sears, Neutron News, vol.3-3, 1992, 26-37           ')
   call write_info(' ')
  endif
 endif

 RETURN
end subroutine write_mu
!-------------------------------------------------------------------

subroutine neutron_cross_section()
 USE atome_module
 USE cryscalc_module, ONLY         : ON_SCREEN, nb_atoms_type, wavelength, atom_typ, SFAC_type, sfac_number, num_atom, &
                                     message_text, debug_proc
 USE IO_module,  ONLY              : write_info

  implicit none
  ! local variables
  INTEGER                          :: i
  REAL, parameter                  :: pi = 3.14159
  REAL                             :: neutron_velocity

  if(debug_proc%level_2)  call write_debug_proc_level(2, "neutron_cross_section")
  
  neutron_velocity = 3952. / wavelength
  if(ON_SCREEN) then
   WRITE(message_text,'(1x,a,F7.2)')        '    neutron velocity (m/s): ', neutron_velocity
   call write_info(TRIM(message_text))

   call write_info(' ')
   call write_info(' Atom   num. at.   nb of atoms     coh_scat.     inc_scat.      abs. (2200m/s)              abs. ')
   call write_info(' ')
  endif

  do i=1, nb_atoms_type
   atom(num_atom(i))%N_SED_coh = atom(num_atom(i))%bcoh **2. *4. * pi
   atom(num_atom(i))%N_SED_inc = atom(num_atom(i))%SEDinc
   atom(num_atom(i))%N_SE_absorption = atom(num_atom(i))%SEA*2200./neutron_velocity

   if(ON_SCREEN) then
    WRITE(message_text,'(2x,a,i8,F8.2,4(10x,F8.3))') SFAC_type(i), num_atom(i), sfac_number(i),                 &
                                                     atom(num_atom(i))%N_SED_coh,  atom(num_atom(i))%N_SED_inc, &
													 atom(num_atom(i))%SEA,                                     &
                                                     atom(num_atom(i))%N_SE_absorption
    call write_info(TRIM(message_text))
   endif
  end do

 RETURN
end subroutine neutron_cross_section

!--------------------------------------------------------------
subroutine neutron_absorption()
 USE atome_module
 USE cryscalc_module, ONLY          : absorption, nb_atoms_type, nb_at, num_atom

 implicit none
   INTEGER                    :: i
   REAL                       :: N_absorption, N_diff_coh, N_diff_incoh

   N_absorption  = 0.
   N_diff_coh    = 0.
   N_diff_incoh = 0.
   absorption%mu = 0.

   do i=1, nb_atoms_type
    N_absorption = N_absorption  + nb_at(i)*  atom(num_atom(i))%N_SE_absorption * 1.E-24
    N_diff_coh   = N_diff_coh    + nb_at(i)*  atom(num_atom(i))%N_SED_coh * 1.E-24
    N_diff_incoh = N_diff_incoh  + nb_at(i)*  atom(num_atom(i))%N_SED_inc * 1.E-24
   end do
   absorption%mu =  N_absorption  + N_diff_coh + N_diff_incoh

 RETURN
end subroutine    neutron_absorption

!-------------------------------------------------

subroutine X_attenuation_calculation()
 USE atome_module,    ONLY : atom
 USE cryscalc_module,  ONLY : ON_SCREEN, absorption, nb_atoms_type, SFAC_type, SFAC_number, num_atom, nb_at, message_text
 USE wavelength_module
 USE IO_module,       ONLY : write_info

 implicit none
  INTEGER                    :: i, k


 IF(X_target(1)%logic) THEN   ! Ag
   do i=1, nb_atoms_type
    atom(num_atom(i))%X_attenuation =  atom(num_atom(i))%tics(1)
   end do
 ELSEIF(X_target(2)%logic) then          ! Mo
   do i=1, nb_atoms_type
    atom(num_atom(i))%X_attenuation =  atom(num_atom(i))%tics(2)
   end do
 elseif (X_target(3)%logic) then         ! Cu
   do i=1, nb_atoms_type
    atom(num_atom(i))%X_attenuation =  atom(num_atom(i))%tics(3)
   end do
 elseif (X_target(4)%logic) then         ! Ni
   do i=1, nb_atoms_type
    atom(num_atom(i))%X_attenuation =  atom(num_atom(i))%tics(4)
   end do
 elseif (X_target(5)%logic) then         ! Co
   do i=1, nb_atoms_type
    atom(num_atom(i))%X_attenuation =  atom(num_atom(i))%tics(5)
   end do
 elseif (X_target(6)%logic) then         ! Fe
   do i=1, nb_atoms_type
    atom(num_atom(i))%X_attenuation =  atom(num_atom(i))%tics(6)
   end do
 elseif (X_target(7)%logic) then         ! Cr
   do i=1, nb_atoms_type
    atom(num_atom(i))%X_attenuation =  atom(num_atom(i))%tics(7)
   end do
 end if

 if(ON_SCREEN) then
  call write_info('')
  call write_info(' Atom      Z       nb of at.   Total cross section (barns)')
  call write_info(' ')
 endif

 absorption%mu = 0.
 do i=1, nb_atoms_type
  if(ON_SCREEN) then
   !WRITE(message_text,'(2x,a,i8,F8.2,1(10x,F12.3))') SFAC_type(i), num_atom(i), SFAC_number(i), atom(num_atom(i))%X_attenuation
   WRITE(message_text,'(2x,a,i8,F8.2,1(10x,F12.3))') SFAC_type(i), atom(num_atom(i))%Z, SFAC_number(i),  &
                                                     atom(num_atom(i))%X_attenuation
   call write_info(trim(message_text))
  endif
  absorption%mu = absorption%mu +  nb_at(i)*atom(num_atom(i))%X_attenuation*1.E-24
 end do


 return
end subroutine X_attenuation_calculation


!--------------------------------------------------------------
subroutine F000_calculation(beam)
 USE cryscalc_module
 USE atomic_data 
 USE IO_module, ONLY : write_info

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: beam
  INTEGER                              :: i

  F000 = 0.    ! pouvoir diffusant dans la maille
  if (beam(1:1)=='X') then
   ! F000: nombre total d'electrons dans la maille
   do i=1, nb_atoms_type
    !F000 = F000 + SFAC_number(i)*num_atom(i)
	F000 = F000 + SFAC_number(i)*atom(Num_atom(i))%Z
   end do
  else
   ! F000: somme des
  end if

  if(ON_SCREEN) then
   call write_info('')
   WRITE(message_text, '(a,F8.1)')    '   >> F000 = ', real(F000)
   call write_info(TRIM(message_text))
   call write_info('')
  endif

  IF(keyword_create_CIF)  call write_CIF_file('F000')



 return
end subroutine F000_calculation

!--------------------------------------------------------------


subroutine Absorption_routine
 use cryscalc_module, only : keyword_CELL, keyword_WAVE, keyword_CONT, keyword_CHEM, message_text, debug_proc
 USE IO_module ,     ONLY : write_info

 implicit none
 
 if(debug_proc%level_2)  call write_debug_proc_level(2, "absorption_routine")
 
 if(.not. keyword_CELL) then
  write(unit=message_text, fmt='(2a)') ' Cell parameters are mandatory for a absorption coefficient calculation.', & 
                                       ' Please enter CELL keyword to input cell parameters.'
  call write_info(trim(message_text))
  call write_info('')
  return
 endif

 if(.not. keyword_WAVE) then
  write(unit=message_text, fmt='(2a)') ' Wavelength value of incident beam is mandatory for a absorption coefficient', &
                                       ' calculation. Please enter WAVE keyword to input incident wavelength.'
  call write_info(trim(message_text))
  call write_info('')
  return
 endif
 
! if(.not. keyword_CONT .and. .not. keyword_CHEM) then
!  write(unit=message_text, fmt='(2a)') ' Cell content is mandatory for a absorption coefficient calculation.', &
!                                       ' Please enter CONT keyword to input cell content.'
!  call write_info(trim(message_text))
!  call write_info('')
!  return
! endif
 
 call atomic_identification()
 call atomic_density_calculation()
 call absorption_calculation
 
 return
end subroutine Absorption_routine
