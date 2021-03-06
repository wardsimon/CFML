!     Last change:  TR   27 Nov 2006    4:26 pm
subroutine write_mag_lines ()
 USE cryscalc_module, ONLY : mag_atom_label, message_text, known_mag_lines
 use magnetic_module
 USE macros_module,   ONLY : u_case
 USE IO_module,       ONLY : write_info
  implicit none
   INTEGER           :: i, i1, long


  IF(.NOT. known_mag_lines) then
   call set_mag_3d
   call set_mag_4f()
   known_mag_lines = .true.
  endif

  i1=0
  long = LEN_TRIM(mag_atom_label)
  ! if (long==0) return
  do i=1, nb_mag_3d_lines
   IF(u_case(mag_atom_label(1:long)) == u_case(mag_3d(i)%ion(1:long))) then
    i1=i1+1
    IF(i1==1) then
     call write_info('')
     call write_info('   ion    ec    level   g       p_calc_g     S     p_calc_S   p_exp ')
     call write_info('')
    endif
    WRITE(message_text, '(4(3x,a),4x,2(3x,a),4x,2(5x,a))') mag_3d(i)%ion,      mag_3d(i)%ec, mag_3d(i)%level,    mag_3d(i)%g,    &
                                                           mag_3d(i)%p_calc_g, mag_3d(i)%S,  mag_3d(i)%p_calc_S, mag_3d(i)%p_exp
    call write_info(TRIM(message_text))
   endif
  end do
  !IF(i1 /=0) return

  !i1=0
  long = LEN_TRIM(mag_atom_label)
  ! if (long==0) return
  do i=1, nb_mag_4f_lines
   IF(u_case(mag_atom_label(1:long)) == u_case(mag_4f(i)%ion(1:long))) then
    i1=i1+1
    IF(i1==1) then
     call write_info('')
     call write_info('   ion    ec     level    L    S    J      g        p       gJ ')
     call write_info('')
    endif
    WRITE(message_text, '(9(3x,a))') mag_4f(i)%ion, mag_4f(i)%ec, mag_4f(i)%level, mag_4f(i)%L,    &
                                     mag_4f(i)%S,   mag_4f(i)%J,  mag_4f(i)%g,     mag_4f(i)%p, mag_4f(i)%gJ
    call write_info(TRIM(message_text))
   endif
  end do

   if(i1==0) then
    call write_info('')
    call write_info(' Please check the input magnetic atom (ex: Mn3+ and not Mn+3)')
    call write_info('')
   end if


  return
end subroutine write_mag_lines

!------------------------------------------------------------------
subroutine write_shannon_lines ()
 USE cryscalc_module, ONLY : shannon_atom_label, message_text, known_shannon_lines
 use shannon_module
 USE macros_module,   ONLY : u_case
 USE IO_module,       ONLY : write_info
  implicit none
   INTEGER           :: i, i1, long


  IF(.NOT. known_shannon_lines) then
   call set_shannon()
   known_shannon_lines = .true.
  endif

  i1=0
  long = LEN_TRIM(shannon_atom_label)
  ! if (long==0) return
  do i=1, nb_shannon_lines
   IF(u_case(shannon_atom_label(1:long)) == SHANNON(i)%ion(1:long)) then
    i1=i1+1
    IF(i1==1) then
     call write_info('    ion   ec     CN      SP    CR     IR ')
     call write_info('')
    endif
    WRITE(message_text, '(4(3x,a), 2x, 2F7.3)') SHANNON(i)%ion, SHANNON(i)%ec, SHANNON(i)%CN, SHANNON(i)%SP,    &
                                                SHANNON(i)%CR,  SHANNON(i)%IR
    call write_info(TRIM(message_text))
   endif
  end do

  return
end subroutine write_shannon_lines

!---------------------------------------------------------------
subroutine write_atomic_features()
 USE cryscalc_module,                 ONLY : mendel_atom_label, mendel_atom_nb, message_text, pi, &
                                             known_atomic_label, known_atomic_features, known_data_neutrons, known_data_x, &
                                             debug_proc
 USE CFML_Scattering_Chemical_Tables, ONLY : set_chem_info, set_xray_form, set_delta_fp_fpp
 USE IO_module,                       ONLY : write_info
 use atomic_data
 USE Neutrons_data
 USE Xrays_data

 implicit none
  INTEGER                  :: i
  CHARACTER (LEN=4)        :: symb_car
  CHARACTER (LEN=6)        :: symbol
  LOGICAL                  :: ok

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_atomic_features")

  if(.not. known_atomic_label)  then
   call definition_atomic_label    ! %symbol, %name
   known_atomic_label = .true.
  endif
  if(.not. known_atomic_features) then
   call def_atomic_features()      ! %weight, %config_electr, %density, %radius
   known_atomic_features = .true.
  endif
  if(.not. known_data_neutrons)  then
   call definition_data_neutrons()
   call definition_data_neutrons_rare_earth()
   known_data_neutrons = .true.
  endif
  if(.not. known_data_X) then
   call definition_data_X
   known_data_X = .true.
  endif

  ! definition des caracteristiques atomiques, coef. de diffusion=f(s), diff. anomale =f(wave)
  ! dans CFML
  CALL set_chem_info
  call Set_Xray_Form    ! definition des facteurs de diffusion
  call Set_Delta_Fp_Fpp ! definition des facteurs de diffusion anomale fp et fpp



  do i=1, mendel_atom_nb
   call verif2_atomic_symbol(mendel_atom_label(i), symb_car, ok)
   if(.not. ok) then
    call verif_atomic_number(mendel_atom_label(i), symb_car, ok)
    if(ok) call verif_atomic_symbol(mendel_atom_label(i), symb_car, ok)
    IF(.not. ok) return
   end if

   call write_info('')
   call symbol_to_write(mendel_atom_label(i),  symbol)
   WRITE(message_text, '(3a)')        '  >> ', TRIM(symbol), ':'
   call write_info(TRIM(message_text))
   call write_info('')



   call write_atomic_features_from_CFML(symb_car)
   call write_atomic_features_from_ATOME_module(symb_car)
   call write_X_rays_scatt_factors(mendel_atom_label(i))
   call write_X_rays_anomalous_factors(symb_car)

  end do

 RETURN
end subroutine  write_atomic_features

!-----------------------------------------------------------------------------
subroutine verif_atomic_symbol(input_symb_car, symb_car, ok)
 USE atome_module,    ONLY : atom
 USE macros_module,   ONLY : u_case
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : message_text

 implicit none
  CHARACTER (LEN=*), INTENT(IN)   :: input_symb_car
  CHARACTER (LEN=*), INTENT(OUT)  :: symb_car
  LOGICAL,           INTENT(OUT)  :: ok
  INTEGER                         :: i1, j

  i1 = INDEX(input_symb_car, '+')
  IF(i1>1) then
   symb_car = input_symb_car(1:i1-1)
  else
   i1 = INDEX(input_symb_car, '-')
   IF(i1>1) then
    symb_car = input_symb_car(1:i1-1)
   else
    symb_car = input_symb_car(1:)
   endif
  endif

  j = 0
  do
   j=j+1
   if (j>201) then
    !call write_info('')
    !call write_info('   !!! '//TRIM(input_symb_car)// ': incorrect atomic symbol !!')
    ok = .false.
    return
   end if
   if (u_case(symb_car(1:))== atom(j)%symbol(1:)) exit    ! atome_module

  end do

  ok = .true.
  !call write_info('')

  !call symbol_to_write(input_symb_car,  symbol)
  !WRITE(message_text, '(3a)')        '  >> ', TRIM(symbol), ':'
  !call write_info(TRIM(message_text))
  !call write_info('')

 return
end subroutine verif_atomic_symbol

!!-----------------------------------------------------------------------------!!
! le symbole ne peut pas etre sous la forme M+x ou Mx+
subroutine verif2_atomic_symbol(input_symb_car, symb_car, ok)

 USE atome_module,    ONLY : atom
 USE macros_module,   ONLY : u_case, remove_car
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : message_text

 implicit none
  CHARACTER (LEN=*), INTENT(IN)   :: input_symb_car
  CHARACTER (LEN=*), INTENT(OUT)  :: symb_car
  LOGICAL,           INTENT(OUT)  :: ok
  INTEGER                         :: i
  character (len=1)               :: i_string

  symb_car = input_symb_car
  symb_car = remove_car(symb_car, '+')
  symb_car = remove_car(symb_car, '-')
  do i=1, 9
   write(i_string, '(i1)') i
  symb_car = remove_car(symb_car, i_string)
  end do

  ok = .true.
  return

end subroutine verif2_atomic_symbol


!-----------------------------------------------------------------------------
subroutine verif_atomic_number(input_symb_car, symb_car, ok)
 USE atome_module,    ONLY : atom
 USE macros_module,   ONLY : u_case
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : message_text

 implicit none
  CHARACTER (LEN=*), INTENT(INOUT) :: input_symb_car
  CHARACTER (LEN=*), INTENT(OUT)   :: symb_car
  LOGICAL,           INTENT(OUT)   :: ok
  INTEGER                          :: atom_number

  ok = .false.
  read(input_symb_car, *) atom_number
  if(atom_number == 0) return

  if(atom_number == 201) then
   input_symb_car = "D"
   symb_car       = "D"
  elseif(atom_number > 98) then
   call write_info('')
   call write_info('   !!! '//TRIM(input_symb_car)// ': incorrect atomic number !!')
   return
  else
   input_symb_car = atom(atom_number)%symbol
  endif

  ok = .true.

  !call write_info('')
  !WRITE(message_text, '(3a)')        '  >> ', TRIM(input_symb_car), ':'
  !call write_info(TRIM(message_text))
  !call write_info('')

 return
end subroutine verif_atomic_number

!-----------------------------------------------------------------------------
subroutine write_atomic_features_from_ATOME_module(symb_car)
 USE atome_module,    ONLY : atom
 USE macros_module,   ONLY : u_case
 USE Io_module,       ONLY : write_info
 USE cryscalc_module, ONLY : message_text, pi, debug_proc
 USE wavelength_module
 implicit none
  CHARACTER (LEN=*), INTENT(IN)  :: symb_car
  INTEGER                        :: j, k

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_atomic_features_from_atome_module")

   j=0
   do
    j=j+1
    if (u_case(symb_car(1:))== atom(j)%symbol(1:)) exit    ! atome_module
   end do

   WRITE(message_text, '(a,F8.4)')    '    . Atomic density           : ', atom(j)%density
   call write_info(TRIM(message_text))
   WRITE(message_text, '(2a)')        '    . Electronic configuration : ', atom(j)%config_electr
   call write_info(TRIM(message_text))

   call write_info('')
   WRITE(message_text, '(a)')         '   . Neutron data: '
   call write_info(TRIM(message_text))
   IF(ABS(atom(j)%bcoh) > 0.00001) then
    WRITE(message_text, '(a,F10.4,a)') '    . coherent scattering length:           ', atom(j)%bcoh,                ' 10-12 cm'
    call write_info(TRIM(message_text))
    WRITE(message_text, '(a,F10.4,a)') '    . coherent scattering cross section:    ', 4 * pi * atom(j)%bcoh ** 2,  ' barns'
    call write_info(TRIM(message_text))
   else
    WRITE(message_text, '(a)') '    . coherent scattering length:           unknown'
    WRITE(message_text, '(a)') '    . coherent scattering cross section:    unknown'
   endif
   IF(atom(j)%SEDinc > 0.0001) then
    WRITE(message_text, '(a,F10.4,a)') '    . incoherent scattering cross section:  ', atom(j)%SEDinc,              ' barns'
   else
    WRITE(message_text, '(a)')         '    . incoherent scattering cross section:  unknown or zero'
   endif
   call write_info(TRIM(message_text))
   IF(atom(j)%SEA > 0.0001) then
    WRITE(message_text, '(a,F10.4,a)') '    . absorption cross section (v=2200m/s): ', atom(j)%SEA,                 ' barns'
   else
    WRITE(message_text, '(a)')         '    . absorption cross section (v=2200m/s): unknown or zero'
   endif
   call write_info(TRIM(message_text))
   call write_info('')
   IF(j>200) return

   WRITE(message_text, '(a)')         '   . X-rays data: '
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a)')         '      - Mass Attenuation Coefficient:'
   call write_info(TRIM(message_text))

   do k=1, tabulated_target_nb
    IF(atom(j)%cam(k) > 0.) then
     WRITE(message_text, '(3a,F8.3,a)')  '        . for ', X_target(k)%label, ' Ka:     ',atom(j)%cam(k) ,' cm2/g'
    else
     WRITE(message_text, '(3a)')         '        . for ', X_target(k)%label, ' Ka:     unknown'
    endif
    call write_info(TRIM(message_text))
   end do


   WRITE(message_text, '(a)')         '      - Total interaction cross sections:'
   call write_info(TRIM(message_text))
   do k=1, tabulated_target_nb
    if (atom(j)%tics(k) > 0.) then
     WRITE(message_text, '(3a,F10.3,a)') '        . for ', X_target(k)%label, ' Ka:     ',atom(j)%tics(k) ,' barns'
    else
     WRITE(message_text, '(3a)')         '        . for ', X_target(k)%label, ' Ka:     unknown'
    endif
    call write_info(TRIM(message_text))
   END do



 return
end subroutine write_atomic_features_from_ATOME_module

!-------------------------------------------------------------------
subroutine write_atomic_features_from_CFML(symb_car)
 USE CFML_Scattering_Chemical_Tables, ONLY : chem_info, num_chem_info
 USE macros_module,                   ONLY : u_case, l_case
 USE Io_module,                       ONLY : write_info
 USE cryscalc_module,                 ONLY : message_text, mendel_atom_label, debug_proc
 implicit none
  character (LEN=*), INTENT(IN) :: symb_car
  !local variables:
  INTEGER                       :: i, j, n_ox
  CHARACTER (LEN=16)            :: fmt_, symbol

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_atomic_features_from__CFML")

  do j=1, Num_Chem_Info
    if (u_case(symb_car(1:))== u_case(chem_info(j)%symb(1:))) exit
  end do


  call write_info('   . Atomic features:')
  call symbol_to_write(chem_info(j)%symb, symbol)
  !if(len_trim(chem_info(j)%symb) == 1) then
  ! WRITE(message_text,'(2a)')      '    . Symbole                  : ', u_case(chem_info(j)%symb(1:1))
  !else
  ! WRITE(message_text,'(2a)')      '    . Symbole                  : ', trim(u_case(chem_info(j)%symb(1:1))//l_case(chem_info(j)%symb(2:)))
  !end if
  WRITE(message_text,'(2a)')      '    . Symbole                  : ', trim(symbol)

  call write_info(TRIM(message_text))
  WRITE(message_text,'(2a)')      '    . Name                     : ', TRIM(chem_info(j)%NAME)
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,I3)')    '    . Atomic number            : ', chem_info(j)%Z
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.4)')  '    . Atomic weight            : ', chem_info(j)%AtWe
  call write_info(TRIM(message_text))
  IF(chem_info(j)%RCov > 0.01) then
   WRITE(message_text,'(a,F8.4)') '    . Covalent radius (A)      : ', chem_info(j)%RCov
   call write_info(TRIM(message_text))
  endif
  IF(chem_info(j)%RWaals > 0.01) then
   WRITE(message_text,'(a,F8.4)') '    . van der Waals radius (A) : ', chem_info(j)%RWaals
   call write_info(TRIM(message_text))
  endif
  IF(chem_info(j)%VAtm > 0.01) then
   WRITE(message_text,'(a,F8.4)') '    . Atomic volume (A3)       : ', chem_info(j)%VAtm
   call write_info(TRIM(message_text))
  endif

  do i=5,1,-1
   IF( chem_info(j)%Oxid(i) /= 0) exit
  end do
  n_ox = i
  IF(n_ox /=0) then
   WRITE(fmt_, '(a,i1,a)') '(a,', n_ox,'i3)'
   WRITE(message_text, fmt_)      '    . Oxydation state          : ', chem_info(j)%Oxid(1:n_ox)
   call write_info(TRIM(message_text))

   WRITE(fmt_, '(a,i1,a)') '(a,', n_ox,'F8.4)'
   WRITE(message_text, fmt_)      '    . Ionic radius (A)         : ', chem_info(j)%Rion(1:n_ox)
   call write_info(TRIM(message_text))
  endif

end subroutine write_atomic_features_from_CFML
!-------------------------------------------------------------------

subroutine write_X_rays_scatt_factors(symb_car)
 USE CFML_Scattering_Chemical_Tables
 USE cryscalc_module,                  ONLY : message_text, mendel_atom_nb,  mendel_atom_label, mendel_plot, PGF_file, PGF_data, &
                                              winplotr_exe, allocate_PGF_data_arrays, debug_proc
 USE macros_module,                    ONLY : u_case
 USE IO_module,                        ONLY : write_info


    !!----    Coefficients for calculating the X-ray scattering factors
    !!----        f(s) = Sum_{i=1,4} { a(i) exp(-b(i)*s^2) } + c
    !!----
    !!----    where s=sinTheta/Lambda

 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: symb_car
  CHARACTER (len=6)             :: symbol
  INTEGER                       :: i,j, k
  REAL, DIMENSION(100)          :: stl, f0
  LOGICAL                       :: ok


  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_X_rays_scatt_factors")

  !call Set_Xray_Form ! definition des facteurs de diffusion

  call Allocate_PGF_data_arrays(76)

  !do i=1, mendel_atom_nb
   ok = .false.
   do j=1, Num_Xray_Form
    if (u_case(symb_car(1:))== u_case(Xray_form(j)%symb(1:))) then
     ok = .true.
     exit
    endif
   end do

   IF(.not. ok) then
    call write_info('')
    WRITE(message_text, '(2a)') '  . X-ray scattering factors not tabulated for: ', TRIM(symb_car)
    call write_info(TRIM(message_text))
    call write_info('')
    return
   endif

   call write_info('')
   call write_info('  . Coefficients for calculating the X-ray scattering factors')
   call write_info( '       f(s) = Sum_{i=1,4} { a(i) exp(-b(i)*s^2) } + c')
   call write_info( '')
   call write_info( '    where s=sinTheta/Lambda')
   call write_info( '')
   !call write_info( '   => X-ray scattering coeff. (A1,A2,A3,A4,B1,B2,B3,B4,C):')
   call write_info( '   => X-ray scattering coeff. (A1,B1,A2,B2,A3,B3,A4,B4,C):')
   call write_info( '')

   call symbol_to_write(Xray_form(j)%symb,  symbol)
   write(message_text, '(2a)')      '   Element: ', trim(symbol)
   !write(message_text, '(2a)')      '   Element: ', u_case(Xray_form(j)%symb)
   call write_info(TRIM(message_text))
   !WRITE(message_text, '(a,9F9.4)') '   coef:    ', Xray_form(j)%A(1:4), Xray_form(j)%B(1:4), Xray_form(j)%C
   WRITE(message_text, '(a,9F9.4)') '   coef:    ', Xray_form(j)%A(1), Xray_form(j)%B(1), &
                                                    Xray_form(j)%A(2), Xray_form(j)%B(2), &
                                                    Xray_form(j)%A(3), Xray_form(j)%B(3), &
                                                    Xray_form(j)%A(4), Xray_form(j)%B(4), Xray_form(j)%C
   call write_info(TRIM(message_text))
   call write_info('')

   write(message_text, '(5x,a,F8.3)')     'at s = 0 : f(O) = ', sum(Xray_form(j)%A(1:4)) + Xray_form(j)%C
   call write_info(trim(message_text))

   IF(mendel_plot) then
    do i=1,76
     stl(i) = 0.02*(i-1)
     f0(i) = 0
     do k=1,4
      f0(i) = f0(i) + Xray_form(j)%A(k)*EXP(-Xray_form(j)%B(k)*stl(i)**2)
     end do
     f0(i) = f0(i) + Xray_form(j)%C
     WRITE(pgf_data%string(i), '(a,F6.2,a,F6.3)') 'stl (A-1) = ', stl(i), '  f0 = ', f0(i)
    end do
    WRITE(Pgf_file%NAME, '(a)') TRIM(u_case(Xray_form(j)%symb))//'_X_scat_fact_versus_stl.pgf'
    call create_PGF_file(TRIM(Pgf_file%name), stl, f0, pgf_data%string, 76, "f0 :"//u_case(Xray_form(j)%symb))
    IF(LEN_TRIM(winplotr_exe) == 0) then
     call write_info('')
     call write_info(' Warning: WinPLOTR is not installed ')
     call write_info('')
    else
     !call system('winplotr '//trim(pgf_file%name), .true. )   ! lf95
     call system('winplotr '//trim(pgf_file%name))             ! g95
    endif
   endif

  !END do
 RETURN
end subroutine write_X_rays_scatt_factors

!-------------------------------------------------------------------
subroutine write_X_rays_anomalous_factors(symb_car)
 USE CFML_Scattering_Chemical_Tables
 USE cryscalc_module,            ONLY : message_text, mendel_atom_nb,  mendel_atom_label, debug_proc
 USE macros_module,              ONLY : u_case
 USE IO_module,                  ONLY : write_info

 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: symb_car
  CHARACTER (len=6)             :: symbol
  INTEGER                       :: j
  LOGICAL                       :: ok

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_X_rays_anomalous_factors")

  !call Set_Delta_Fp_Fpp ! definition des facteurs de diffusion anomale fp et fpp

  !do i=1, mendel_atom_nb
  ! i1 = INDEX(mendel_atom_label(i), '+')
  ! IF(i1>1) then
  !  symb_car = mendel_atom_label(i)(1:i1-1)
  ! else
  !  i1 = INDEX(mendel_atom_label(i), '-')
  !  IF(i1>1) then
  !   symb_car = mendel_atom_label(i)(1:i1-1)
  !  else
  !   symb_car = mendel_atom_label(i)(1:)
  !  endif
  ! endif

   ok = .false.
   do j=1, Num_Delta_Fp
    if (u_case(symb_car(1:))== u_case(Anomalous_ScFac(j)%symb(1:))) then
     ok = .true.
     exit
    endif
   end do

   IF(.not. ok) then
    call write_info('')
    WRITE(message_text, '(2a)') '  . Values of Delta_Fp and Delta_Fpp not tabulated for: ', TRIM(symb_car)
    call write_info(TRIM(message_text))
    call write_info('')
    return
   endif

   call write_info('')
   call write_info('  . Values of Delta_Fp and Delta_Fpp anomalous dispersion coefficients')
   call write_info('    for 5 common radiations')
   call write_info('')
   call write_info('    Wavelenghts:       Cr        Fe        Cu        Mo        Ag')
   call write_info('         Lambda   2.28962   1.93597   1.54051   0.70926  0.556363')
   call write_info('')
   call symbol_to_write(Anomalous_ScFac(j)%symb, symbol)
   write(message_text, '(2a)')      '   Element: ', trim(Symbol)
   !write(message_text, '(2a)')      '   Element: ', u_case(Anomalous_ScFac(j)%symb)
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,5(2x,F8.3))') '   Fp coef:    ', Anomalous_ScFac(j)%fp(1:5)
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,5(2x,F8.3))') '   Fpp coef:   ', Anomalous_ScFac(j)%fpp(1:5)
   call write_info(TRIM(message_text))

   call write_info('')

  !END do

END subroutine write_X_rays_anomalous_factors

!-------------------------------------------------------------------
subroutine WRITE_data(input_string)
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : message_text, pgf_data, pgf_file,                           &
                             known_atomic_label, known_atomic_features, known_data_neutrons, known_data_x,  &
                             data_neutrons_PLOT, data_Xrays_PLOT, data_atomic_density_PLOT, data_atomic_weight_PLOT,      &
                             data_atomic_radius_PLOT,   winplotr_exe, allocate_PGF_data_arrays, debug_proc
 USE wavelength_module
 use atomic_data
 use atome_module
 USE Neutrons_data
 USE Xrays_data

 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: input_string
  INTEGER                       :: i, k
  REAL                          :: Xmin, Xmax, Ymin, Ymax
  REAL                          :: lambda_E


  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_data ("//trim(input_string)//")")

  !call Allocate_PGF_data_arrays(96)
  call Allocate_PGF_data_arrays(100)  ! nov. 2014 : permet de visualiser Gd_157 (100 valeurs)

  if(.not. known_atomic_label)  then
   call definition_atomic_label    ! %symbol, %name
   known_atomic_label = .true.
  endif
  if(.not. known_atomic_features) then
   call def_atomic_features()      ! %weight, %config_electr, %density, %radius
   known_atomic_features = .true.
  endif
  if(.not. known_data_neutrons)  then
   call definition_data_neutrons()
   call definition_data_neutrons_Rare_earth()
   known_data_neutrons = .true.
  endif
  if(.not. known_data_X) then
   call definition_data_X
   known_data_X = .true.
  endif

  select case (input_string)
      case ('neutrons')
       call write_info('     Element        coherent scattering length (10-12 cm)')
       call write_info('                    incoherent scattering cross section (in barns)')
       call write_info('                    absorption cross section (in barns) for v=2200m/s')
       call write_info('')
       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,10x,3F15.5)') i,  atom(i)%symbol, atom(i)%bcoh, atom(i)%sedinc, atom(i)%sea
        call write_info(TRIM(message_text))
       END do

       IF(data_neutrons_PLOT) then
        pgf_file%name = 'neutrons_bcoh.pgf'
        do i=1,96
         pgf_data%X(i) = REAL(i)
         WRITE(pgf_data%string(i), '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%bcoh
        END do
        call create_PGF_file(TRIM(pgf_file%name), pgf_data%X, atom(1:96)%bcoh, pgf_data%string, 96, "bcoh")
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         !call system('winplotr '//trim(pgf_file%name), .true. )
         call system('winplotr '//trim(pgf_file%name))
        endif

        pgf_file%name = 'neutrons_sedinc.pgf'
        do i=1,96
         WRITE(pgf_data%string(i), '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%sedinc
        end do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom(1:96)%sedinc, pgf_data%string, 96, "sedinc")
        !IF(LEN_TRIM(winplotr_exe) /=0) call system('winplotr '//trim(pgf_file%name), .true. )    ! lf95
        IF(LEN_TRIM(winplotr_exe) /=0) call system('winplotr '//trim(pgf_file%name))              ! g95

        pgf_file%name = 'neutrons_sea.pgf'
        do i=1,96
         WRITE(pgf_data%string(i), '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%sea
        END do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom(1:96)%sea, pgf_data%string, 96, "sea")
        !IF(LEN_TRIM(winplotr_exe) /=0) call system('winplotr '//trim(pgf_file%name), .true. )
        IF(LEN_TRIM(winplotr_exe) /=0) call system('winplotr '//trim(pgf_file%name))
       END if

      case ('neutrons_RE_ALL')
       call write_info('     Element : '//trim(current_RE_label))
       call write_info('')
       call write_info('       b : neutrons scattering length (10^-12cm)')
       call write_info('')
       call write_info('            E(ev)     l(A)     Re(b)     Im(b)       |b|')
       call write_info('')

       if(current_RE == 1) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Sm_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 2) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Sm_149(1:current_RE_n, 1:4)
       elseif(current_RE == 3) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Eu_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 4) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Eu_151(1:current_RE_n, 1:4)
       elseif(current_RE == 5) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Gd_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 6) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Gd_155(1:current_RE_n, 1:4)
       elseif(current_RE == 7) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Gd_157(1:current_RE_n, 1:4)
       elseif(current_RE == 8) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Dy_164(1:current_RE_n, 1:4)
       elseif(current_RE == 9) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Er_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 10) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Er_167(1:current_RE_n, 1:4)
       elseif(current_RE == 11) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Yb_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 12) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Yb_168(1:current_RE_n, 1:4)
       elseif(current_RE == 13) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Yb_174(1:current_RE_n, 1:4)
       elseif(current_RE == 14) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Lu_176(1:current_RE_n, 1:4)
      end if
       do i=1, current_RE_n
        lambda_E = 9.045/SQRT(current_RE_data(i, 1)*1000.)
        WRITE(message_text, '(2x,I4,5(2X,F8.3))') i, current_RE_data(i, 1), lambda_E,0.10*current_RE_data(i, 2:4)
        call write_info(TRIM(message_text))
       end do

       IF(data_neutrons_PLOT) then
        pgf_file%name = 'neutrons_'//trim(current_RE_label)//".pgf"
        open(unit=11, file=trim(PGF_file%name))
         WRITE(11,'(a)')      '# .PGF (winPLOTR Graphics file) created by CRYSCALC:'
         WRITE(11,'(a)')      '#'
         WRITE(11,'(a)')      '# DATA FROM:  Atomic data and nuclear data tables 44, 191-207 (1990)'
         WRITE(11,'(a)')      '#             Resonance effects in neutron scattering lengths or rare-earth nuclides'
         WRITE(11,'(a)')      '#             J.E. Lynn and P.A. Seeger, Los Alamos National Laboratory'
         WRITE(11,'(2a)')     '# MAIN LEGEND TEXT: Neutron coherent scattering length for ', trim(current_RE_label)
         WRITE(11,'(a)')      '# X LEGEND TEXT   : Energy (ev)'
         WRITE(11,'(a)')      '# Y LEGEND TEXT   : b (10^-12 cm) [Re, Im, Mod]'
         Xmin = minval(current_RE_data(1:current_RE_n, 1))
         Xmax = maxval(current_RE_data(1:current_RE_n, 1))
         Ymin = 0.1*minval(current_RE_data(1:current_RE_n, 2:4))
         Ymax = 0.1*maxval(current_RE_data(1:current_RE_n, 2:4))
         write(11, '(a,2(1x,F15.5))')  '#   YMIN YMAX     : ', Ymin, Ymax
         WRITE(11,'(a)')               '# NUMBER OF PATTERNS: 3'

        ! data 1 : Re
         WRITE(11,'(a1,80a1)')  '#',('-',i=1,80)
         WRITE(11,'(a,i6)')     '# >>>>>>>> PATTERN #: ', 1
         write(11,'(a)')        '#             TITLE : Re_b = f(energy)'
         write(11,'(a,i8)')     '#  NUMBER OF POINTS : ', current_RE_n
         write(11,'(a)')        '#            MARKER : 4'
         write(11,'(a)')        '#              SIZE : 1.5'
         write(11,'(a)')        '#             STYLE : 1'
         write(11,'(a)')        '#   DATA: X Y COMM'
        do i=1, current_RE_n
         pgf_data%X(i) = current_RE_data(i, 1)
         lambda_E = 9.045/SQRT(current_RE_data(i, 1)*1000.)
         WRITE(pgf_data%string(i), '(3(a,F8.3),a)') 'E(ev)= ', current_RE_data(i, 1), &
                                                    ', l(A)= ', lambda_E, &
                                                    ', Re(b)= ', 0.1*current_RE_data(i,2), ' 10^-12cm'
        end do
        call create_PGF_file_multi(11,pgf_data%X, 0.1*current_RE_data(1:current_RE_n,2), pgf_data%string, current_RE_n)

        ! data 2 : Im
         WRITE(11,'(a1,80a1)')  '#',('-',i=1,80)
         WRITE(11,'(a,i6)')     '# >>>>>>>> PATTERN #: ', 2
         write(11,'(a)')        '#             TITLE : Im_b = f(energy)'
         write(11,'(a,i8)')     '#  NUMBER OF POINTS : ', current_RE_n
         write(11,'(a)')        '#            MARKER : 4'
         write(11,'(a)')        '#              SIZE : 1.5'
         write(11,'(a)')        '#             STYLE : 1'
         write(11,'(a)')        '#   DATA: X Y COMM'
        do i=1, current_RE_n
         pgf_data%X(i) = current_RE_data(i, 1)
         lambda_E = 9.045/SQRT(current_RE_data(i, 1)*1000.)
         WRITE(pgf_data%string(i), '(3(a,F8.3),a)') 'E(ev)= ', current_RE_data(i, 1), &
                                            ', l(A)= ', lambda_E, &
                                            ', Im(b)= ', 0.1*current_RE_data(i,3), ' 10^-12cm'
        end do
        call create_PGF_file_multi(11,pgf_data%X, 0.1*current_RE_data(1:current_RE_n,3), pgf_data%string, current_RE_n)

        ! data 3 : module
         WRITE(11,'(a1,80a1)')  '#',('-',i=1,80)
         WRITE(11,'(a,i6)')     '# >>>>>>>> PATTERN #: ', 3
         write(11,'(a)')        '#             TITLE : Mod(b) = f(energy)'
         write(11,'(a,i8)')     '#  NUMBER OF POINTS : ', current_RE_n
         write(11,'(a)')        '#            MARKER : 4'
         write(11,'(a)')        '#              SIZE : 1.5'
         write(11,'(a)')        '#             STYLE : 1'
         write(11,'(a)')        '#   DATA: X Y COMM'
         do i=1, current_RE_n
          pgf_data%X(i) = current_RE_data(i, 1)
          lambda_E = 9.045/SQRT(current_RE_data(i, 1)*1000.)
          WRITE(pgf_data%string(i), '(3(a,F8.3),a)') 'E(ev)= ', current_RE_data(i, 1), &
                                                     ', l(A)= ', lambda_E, &
                                                     ', |b|= ', 0.1*current_RE_data(i,4), ' 10^-12cm'
         end do
         call create_PGF_file_multi(11,pgf_data%X, 0.1*current_RE_data(1:current_RE_n,4), pgf_data%string, current_RE_n)
         WRITE(11,'(a)') '# END OF FILE'
         close(unit=11)
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         call system('winplotr '//trim(pgf_file%name))
        endif
       end if


      case ('neutrons_RE')
       !current_RE_label = RE_label(current_RE)
       call write_info('     Element : '//trim(current_RE_label))
       call write_info('')
       call write_info('       b : neutrons scattering length (10^-12cm)')
       call write_info('')
       call write_info('            E(ev)     l(A)     Re(b)     Im(b)       |b|')
       call write_info('')

       if(current_RE == 1) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Sm_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 2) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Sm_149(1:current_RE_n, 1:4)
       elseif(current_RE == 3) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Eu_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 4) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Eu_151(1:current_RE_n, 1:4)
       elseif(current_RE == 5) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Gd_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 6) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Gd_155(1:current_RE_n, 1:4)
       elseif(current_RE == 7) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Gd_157(1:current_RE_n, 1:4)
       elseif(current_RE == 8) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Dy_164(1:current_RE_n, 1:4)
       elseif(current_RE == 9) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Er_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 10) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Er_167(1:current_RE_n, 1:4)
       elseif(current_RE == 11) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Yb_nat(1:current_RE_n, 1:4)
       elseif(current_RE == 12) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Yb_168(1:current_RE_n, 1:4)
       elseif(current_RE == 13) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Yb_174(1:current_RE_n, 1:4)
       elseif(current_RE == 14) then
        current_RE_data(1:current_RE_n, 1:4) = RE_Lu_176(1:current_RE_n, 1:4)
       end if


       do i=1, current_RE_n
        lambda_E = 9.045/SQRT(current_RE_data(i, 1)*1000.)
        WRITE(message_text, '(2x,I4,5(2X,F8.3))') i, current_RE_data(i, 1), lambda_E,0.10*current_RE_data(i, 2:4)
        call write_info(TRIM(message_text))
       end do

      IF(data_neutrons_PLOT) then
        pgf_file%name = 'neutrons_'//trim(current_RE_label)//"_mod.pgf"
        do i=1, current_RE_n
         pgf_data%X(i) = current_RE_data(i, 1)
         lambda_E = 9.045/SQRT(current_RE_data(i, 1)*1000.)
         WRITE(pgf_data%string(i), '(3(a,F8.3),a)') 'E(ev)= ', current_RE_data(i, 1), &
                                            ', l(A)= ', lambda_E, &
                                            ', b= ', 0.1*current_RE_data(i,4), ' 10^-12cm'
        end do
        call create_PGF_file(TRIM(pgf_file%name), pgf_data%X, 0.1*current_RE_data(1:current_RE_n,4), pgf_data%string, &
                             current_RE_n, 'RE_data_mod')

        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         call system('winplotr '//trim(pgf_file%name))
        endif

        pgf_file%name = 'neutrons_'//trim(current_RE_label)//"_Re.pgf"
        do i=1, current_RE_n
         pgf_data%X(i) = current_RE_data(i, 1)
         lambda_E = 9.045/SQRT(current_RE_data(i, 1)*1000.)
         WRITE(pgf_data%string(i), '(3(a,F8.3),a)') 'E(ev)= ', current_RE_data(i, 1), &
                                                    ', l(A)= ', lambda_E, &
                                                    ', b= ', 0.1*current_RE_data(i,2), ' 10^-12cm'
        end do
        call create_PGF_file(TRIM(pgf_file%name), pgf_data%X, 0.1*current_RE_data(1:current_RE_n,2), pgf_data%string, &
                             current_RE_n, 'RE_data_Re')
        IF(LEN_TRIM(winplotr_exe) /=0) call system('winplotr '//trim(pgf_file%name))


        pgf_file%name = 'neutrons_'//trim(current_RE_label)//"_Im.pgf"
        do i=1, current_RE_n
         pgf_data%X(i) = current_RE_data(i, 1)
         lambda_E = 9.045/SQRT(current_RE_data(i, 1)*1000.)
         WRITE(pgf_data%string(i), '(3(a,F8.3),a)') 'E(ev)= ', current_RE_data(i, 1), &
                                                    ', l(A)= ', lambda_E, &
                                                    ', b= ', 0.1*current_RE_data(i,3), ' 10^-12cm'
        end do
        call create_PGF_file(TRIM(pgf_file%name), pgf_data%X, 0.1*current_RE_data(1:current_RE_n,3), pgf_data%string, &
                             current_RE_n, 'RE_data_Im')
        IF(LEN_TRIM(winplotr_exe) /=0) call system('winplotr '//trim(pgf_file%name))



      END if



      case ('xrays')
       call write_info('                    total interaction cross section (in barns)')
       write(message_text, '(2a)')  '     Element               Ag             Mo             Cu             ', &
                                    'Ni             Co             Fe             Cr'

       call write_info(trim(message_text))
       call write_info('')

       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,7F15.5)') i,  atom(i)%symbol,  atom(i)%tics(1), atom(i)%tics(2), atom(i)%tics(3), &
                                                       atom(i)%tics(4), atom(i)%tics(5), atom(i)%tics(6), atom(i)%tics(7)
        call write_info(TRIM(message_text))
       end do
       !IF(data_xrays_PLOT) then
       ! do k = 1, tabulated_target_nb
       !  pgf_file%NAME = X_target(k)%tics_file
       !  do i = 1, 96
       !   pgf_data(i)%X = REAL(i)
       !   WRITE(pgf_data(i)%string, '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%tics(k)
       !  end do
       !  call create_PGF_file(TRIM(PGF_file%NAME), pgf_data%X, atom%tics(k), pgf_data%string, 96, X_target(k)%label//'_tics')
       !  IF(LEN_TRIM(winplotr_exe) == 0) then
       !   call write_info('')
       !   call write_info(' Warning: WinPLOTR is not installed ')
       !   call write_info('')
       !  else
       !   call system('winplotr '//TRIM(PGF_file%NAME), .TRUE.)
       !  endif
       ! end do
       !END if

       IF(data_Xrays_plot) then
        pgf_file%NAME = "X_rays_tics.pgf"
        open(unit=11, file=trim(pgf_file%name))
         WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Total X-rays interaction cross section'
         WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
         WRITE(11,'(a)')      '# Y LEGEND TEXT   : barns'
         Xmin = 1.
         Xmax = 96.
         Ymin = minval(atom%tics(tabulated_target_nb))
         Ymax = maxval(atom%tics(tabulated_target_nb))
         write(11, '(a,2(1x,F6.2),a)') '#   XMIN XMAX     : ', 0., real(int((10*Xmax+1)))/10
         write(11, '(a,2(1x,F15.5))')  '#   YMIN YMAX     : ', Ymin, Ymax
         WRITE(11,'(a,i6)')     '# NUMBER OF PATTERNS: ',  tabulated_target_nb
         do k = 1, tabulated_target_nb
          WRITE(11,'(a1,70a1)')  '#',('-',i=1,70)
          WRITE(11,'(a,i6)')     '# >>>>>>>> PATTERN #: ', k
          write(11,'(2a)')       '#             TITLE : ', X_target(k)%label
          write(11,'(a,i8)')     '#  NUMBER OF POINTS : ', 96
          write(11,'(a)')        '#            MARKER : 4'
          write(11,'(a)')        '#              SIZE : 1.5'
          write(11,'(a)')        '#             STYLE : 1'
          write(11,'(a)')        '#   DATA: X Y COMM'
          do i = 1, 96
           pgf_data%X(i) = REAL(i)
           WRITE(pgf_data%string(i), '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%tics(k)
          end do
          call create_PGF_file_multi(11,pgf_data%X, atom%tics(k), pgf_data%string, 96)
         END do
         WRITE(11,'(a)') '# END OF FILE'
         CLOSE(11)

         IF(LEN_TRIM(winplotr_exe) == 0) then
          call write_info('')
          call write_info(' Warning: WinPLOTR is not installed ')
          call write_info('')
         else
          !call system('winplotr '//TRIM(PGF_file%NAME), .TRUE.) ! lf95
          call system('winplotr '//TRIM(PGF_file%NAME))          ! g95
         endif
       endif


       call write_info('')
       call write_info('                    Mass attenuation coefficient (cm2/g)')
       write(message_text, '(2a)') '     Element               Ag             Mo             Cu             ', &
                                   'Ni             Co             Fe             Cr'
       call write_info(trim(message_text))
       call write_info('')

       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,7F15.5)') i,  atom(i)%symbol, atom(i)%cam(1), atom(i)%cam(2), atom(i)%cam(3), &
                                                       atom(i)%cam(4), atom(i)%cam(5), atom(i)%cam(6), atom(i)%cam(7)
        call write_info(TRIM(message_text))
       end do
       !IF(data_xrays_PLOT) then
       ! do k = 1, tabulated_target_nb
       !  pgf_file%NAME = X_target(k)%cam_file
       !  do i = 1, 96
       !   pgf_data(i)%X = REAL(i)
       !   WRITE(pgf_data(i)%string, '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%cam(k)
       !  end do
       !  call create_PGF_file(TRIM(PGF_file%NAME), pgf_data%X, atom%cam(k), pgf_data%string, 96, X_target(k)%label//'_cam')
       !  IF(LEN_TRIM(winplotr_exe) == 0) then
       !   call write_info('')
       !   call write_info(' Warning: WinPLOTR is not installed ')
       !   call write_info('')
       !  else
       !   call system('winplotr '//TRIM(PGF_file%NAME), .TRUE.)
       !  endif
       ! end do
       !END if
       IF(data_Xrays_plot) then
        pgf_file%NAME = "X_rays_mac.pgf"
        open(unit=11, file=trim(pgf_file%name))
         WRITE(11,'(a)')      '# MAIN LEGEND TEXT: Mass attenuation coefficient'
         WRITE(11,'(a)')      '# X LEGEND TEXT   : Atomic number'
         WRITE(11,'(a)')      '# Y LEGEND TEXT   : cm2/g'
         Xmin = 1.
         Xmax = 96.
         Ymin = minval(atom%cam(tabulated_target_nb))
         Ymax = maxval(atom%cam(tabulated_target_nb))
         write(11, '(a,2(1x,F6.2),a)') '#   XMIN XMAX     : ', 0., real(int((10*Xmax+1)))/10
         write(11, '(a,2(1x,F15.5))')  '#   YMIN YMAX     : ', Ymin, Ymax
         WRITE(11,'(a,i6)')     '# NUMBER OF PATTERNS: ',  tabulated_target_nb
         do k = 1, tabulated_target_nb
          WRITE(11,'(a1,70a1)')  '#',('-',i=1,70)
          WRITE(11,'(a,i6)')     '# >>>>>>>> PATTERN #: ', k
          write(11,'(2a)')       '#             TITLE : ', X_target(k)%label
          write(11,'(a,i8)')     '#  NUMBER OF POINTS : ', 96
          write(11,'(a)')        '#            MARKER : 4'
          write(11,'(a)')        '#              SIZE : 1.5'
          write(11,'(a)')        '#             STYLE : 1'
          write(11,'(a)')        '#   DATA: X Y COMM'
          do i = 1, 96
           pgf_data%X(i) = REAL(i)
           WRITE(pgf_data%string(i), '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%cam(k)
          end do
          call create_PGF_file_multi(11,pgf_data%X, atom%cam(k), pgf_data%string, 96)
         END do
         WRITE(11,'(a)') '# END OF FILE'
         CLOSE(UNIT=11)

         IF(LEN_TRIM(winplotr_exe) == 0) then
          call write_info('')
          call write_info(' Warning: WinPLOTR is not installed ')
          call write_info('')
         else
          !call system('winplotr '//TRIM(PGF_file%NAME), .TRUE.)  ! lf95
          call system('winplotr '//TRIM(PGF_file%NAME))           ! g95
         endif
       endif


      case('density')
       call write_info('     Element        density')
       call write_info('')

       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,10x,F15.5)') i,  atom(i)%symbol, atom(i)%density
        call write_info(TRIM(message_text))
       end do

       IF(data_atomic_density_PLOT) then
        pgf_file%name = 'atomic_density.pgf'
        do i=1, 96
         pgf_data%X(i) = REAL(i)
         WRITE(pgf_data%string(i), '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%density
        end do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%density, pgf_data%string, 96, "density")
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         !call system('winplotr '//trim(pgf_file%name), .true. )  ! lf95
         call system('winplotr '//trim(pgf_file%name))            ! g95
        endif
       ENDIF

      case ('radius')
       call write_info('     Element        atomic radius')
       call write_info('')

       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,10x,F15.5)') i,  atom(i)%symbol, atom(i)%radius
        call write_info(TRIM(message_text))
       end do

       IF(data_atomic_radius_PLOT) then
        pgf_file%name = 'atomic_radius.pgf'
        do i=1, 96
         pgf_data%X(i) = REAL(i)
         WRITE(pgf_data%string(i), '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%radius
        end do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%radius, pgf_data%string, 96, "radius")
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         !call system('winplotr '//trim(pgf_file%name), .true. )  ! lf95
         call system('winplotr '//trim(pgf_file%name))            ! g95
        endif
       ENDIF


      case('weight')
       call write_info('     Element        weight')
       call write_info('')

       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,10x,F15.5)') i,  atom(i)%symbol, atom(i)%weight
        call write_info(TRIM(message_text))
       end do
       if (data_atomic_weight_PLOT) then
        pgf_file%name = 'atomic_weight.pgf'
        do i=1,96
         pgf_data%X(i) = REAL(i)
         WRITE(pgf_data%string(i), '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%weight
        end do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%weight, pgf_data%string, 96, "weight")
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         !call system('winplotr '//trim(pgf_file%name), .true. )    ! lf95
         call system('winplotr '//trim(pgf_file%name))              ! g95
        endif
       end if

      case default
  end select
 RETURN
end subroutine write_data

!!------------------------------------------------------------------------------------------
subroutine symbol_to_write(symbol, ws)
 use macros_module, only : u_case, l_case
 implicit none
  character (len=*), intent(in)  :: symbol
  character (len=4), intent(out) :: ws


  if(len_trim(symbol) == 1) then
   ws = u_case(symbol(1:1))
  else
   ws = u_case(symbol(1:1))//l_case(symbol(2:))
  end if

   return
end subroutine symbol_to_write
