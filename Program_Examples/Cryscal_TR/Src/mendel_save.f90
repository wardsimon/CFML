!     Last change:  TR   30 Aug 2006    6:08 pm
subroutine write_atomic_features()
 USE cryscal_module,                  ONLY : mendel_atom_label, mendel_atom_nb, message_text, pi, &
                                             known_atomic_label, known_atomic_features, known_data_neutrons, known_data_x
 USE CFML_Scattering_Chemical_Tables, ONLY : set_chem_info
 USE macros_module,                   ONLY : u_case
 USE IO_module,                       ONLY : write_info
 USE atome_module,                    ONLY : atom
 use atomic_data
 USE Neutrons_data
 USE Xrays_data

 implicit none
  INTEGER                  :: i, j, i1
  CHARACTER (LEN=4)        :: symb_car

  if(.not. known_atomic_label)  then
   call definition_atomic_label    ! %symbol, %name
   known_atomic_label = .true.
  endif
  if(.not. known_atomic_features) then
   call def_atomic_features()      ! %weight, %config_electr, %densite, %rayon
   known_atomic_features = .true.
  endif
  if(.not. known_data_neutrons)  then
   call definition_data_neutrons()
   known_data_neutrons = .true.
  endif
  if(.not. known_data_X) then
   call definition_data_X
   known_data_X = .true.
  endif

  CALL set_chem_info()

  do i=1, mendel_atom_nb
   i1 = INDEX(mendel_atom_label(i), '+')
   IF(i1>1) then
    symb_car = mendel_atom_label(i)(1:i1-1)
   else
    i1 = INDEX(mendel_atom_label(i), '-')
    IF(i1>1) then
     symb_car = mendel_atom_label(i)(1:i1-1)
    else
     symb_car = mendel_atom_label(i)(1:)
    endif
   endif

   j=0
   do
    j=j+1
    if (j>201) then
     call write_info('')
     call write_info('   !!! '//TRIM(mendel_atom_label(i))// ': incorrect atomic symbol !!')
     return
    end if

    !if (u_case(mendel_atom_label(i)(1:))== atom(j)%symbol(1:)) exit
    if (u_case(symb_car(1:))== atom(j)%symbol(1:)) exit

   end do


   ! j !
   call write_info('')
   WRITE(message_text, '(3a)')        '  >> ', TRIM(mendel_atom_label(i)), ':'
   call write_info(TRIM(message_text))
   call write_info('')
   IF(j>200) then
    call write_atomic_features_from_CFML(108)  ! D
   else
    call write_atomic_features_from_CFML(j)
   endif

   !IF(j>200) then
   !WRITE(message_text, '(a,I3)')      '   . atomic number:                        ', j-200
   !else
   !WRITE(message_text, '(a,I3)')      '   . atomic number:                        ', j
   !endif
   !call write_info(TRIM(message_text))
   !WRITE(message_text, '(2a)')        '   . atomic symbol:                        ', TRIM(atom(j)%symbol)
   !call write_info(TRIM(message_text))
   !WRITE(message_text, '(2a)')        '   . atomic name:                          ', TRIM(atom(j)%NAME)
   !call write_info(TRIM(message_text))
   !WRITE(message_text, '(a,F8.4)')    '   . atomic weight:                        ', atom(j)%weight
   !call write_info(TRIM(message_text))

   WRITE(message_text, '(a,F8.4)')    '    . Atomic density:           ', atom(j)%densite
   call write_info(TRIM(message_text))
   WRITE(message_text, '(2a)')        '    . Electronic configuration: ', atom(j)%config_electr
   call write_info(TRIM(message_text))

   call write_info('')
   WRITE(message_text, '(a)')         '   . Neutron data: '
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F10.4,a)') '    . coherent scattering length:           ', atom(j)%bcoh,                ' 10-12 cm'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F10.4,a)') '    . coherent scattering cross section:    ', 4 * pi * atom(j)%bcoh ** 2,  ' barns'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F10.4,a)') '    . incoherent scattering cross section:  ', atom(j)%SEDinc,              ' barns'
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F10.4,a)') '    . absorption cross section (v=2200m/s): ', atom(j)%SEA,                 ' barns'
   call write_info(TRIM(message_text))
   call write_info('')
   IF(j>200) cycle

   WRITE(message_text, '(a)')         '   . X-rays data: '
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a)')         '      - Mass Attenuation Coefficient:'
   call write_info(TRIM(message_text))

   if (atom(j)%cam_Mo_ka > 0.) then
   WRITE(message_text, '(a,F8.3,a)')  '        . for Mo Ka:     ',atom(j)%cam_Mo_ka ,' cm2/g'
   else
   WRITE(message_text, '(a)')         '        . for Mo Ka:     unknown'
   endif
   call write_info(TRIM(message_text))

   if (atom(j)%cam_Cu_ka > 0.) then
   WRITE(message_text, '(a,F8.3,a)')  '        . for Cu Ka:     ',atom(j)%cam_Cu_ka ,' cm2/g'
   else
   WRITE(message_text, '(a)')         '        . for Cu Ka:     unknown'
   endif
   call write_info(TRIM(message_text))

   if (atom(j)%cam_Co_ka > 0.) then
   WRITE(message_text, '(a,F8.3,a)')  '        . for Co Ka:     ',atom(j)%cam_Co_ka ,' cm2/g'
   else
   WRITE(message_text, '(a)')         '        . for Co Ka:     unknown'
   endif
   call write_info(TRIM(message_text))

   WRITE(message_text, '(a)')         '      - Total interaction cross sections:'
   call write_info(TRIM(message_text))
   if (atom(j)%tics_Mo_ka > 0.) then
   WRITE(message_text, '(a,F10.3,a)') '        . for Mo Ka:     ',atom(j)%tics_mo_ka ,' barns'
   else
   WRITE(message_text, '(a)')         '        . for Mo Ka:     unknown'
   endif
   call write_info(TRIM(message_text))

   if (atom(j)%tics_Cu_ka > 0.) then
   WRITE(message_text, '(a,F10.3,a)') '        . for Cu Ka:     ',atom(j)%tics_cu_ka ,' barns'
   else
   WRITE(message_text, '(a)')         '        . for Cu Ka:     unknown'
   endif
   call write_info(TRIM(message_text))

   if (atom(j)%tics_Co_ka > 0.) then
   WRITE(message_text, '(a,F10.3,a)') '        . for Co Ka:     ',atom(j)%tics_Co_ka ,' barns'
   else
   WRITE(message_text, '(a)')         '        . for Co Ka:     unknown'
   endif
   call write_info(TRIM(message_text))



  end do

  call write_X_rays_scatt_factors()       ! mendel.F90
  call write_X_rays_anomalous_factors()

 RETURN
end subroutine  write_atomic_features

!-------------------------------------------------------------------
subroutine write_atomic_features_from_CFML(input_i)
 USE Scattering_Chemical_Tables, ONLY: chem_info
 USE IO_module,                  ONLY: write_info
 USE cryscal_module,             ONLY : message_text
 implicit none
  INTEGER, INTENT(IN) :: input_i
  !local variables:
  INTEGER             :: i, n_ox
  CHARACTER (LEN=16)  :: fmt_

  call write_info('   . Atomic features:')
  WRITE(message_text,'(2a)')     '    . Symbole:                  ', TRIM(chem_info(input_i)%symb)
  call write_info(TRIM(message_text))
  WRITE(message_text,'(2a)')     '    . Name:                     ', TRIM(chem_info(input_i)%NAME)
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,I3)')   '    . Atomic number:            ', chem_info(input_i)%Z
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.4)') '    . Atomic weight:            ', chem_info(input_i)%AtWe
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.4)') '    . Covalent radius (A):      ', chem_info(input_i)%RCov
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.4)') '    . van der Waals radius (A): ', chem_info(input_i)%RWaals
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.4)') '    . Atomic volume (A3)        ', chem_info(input_i)%VAtm
  call write_info(TRIM(message_text))

  do i=5,1,-1
   IF( chem_info(input_i)%Oxid(i) /= 0) exit
  end do
  n_ox = i
  WRITE(fmt_, '(a,i1,a)') '(a,', n_ox,'i2)'
  WRITE(message_text, fmt_) '    . Oxydation state:          ', chem_info(input_i)%Oxid(1:n_ox)
  call write_info(TRIM(message_text))

  WRITE(fmt_, '(a,i1,a)') '(a,', n_ox,'F8.4)'
  WRITE(message_text, fmt_) '    . Ionic radius (A):         ', chem_info(input_i)%Rion(1:n_ox)
  call write_info(TRIM(message_text))


end subroutine write_atomic_features_from_CFML
!-------------------------------------------------------------------

subroutine write_X_rays_scatt_factors
 USE Scattering_Chemical_Tables
 USE cryscal_module,             ONLY : message_text, mendel_atom_nb,  mendel_atom_label
 USE macros_module,              ONLY : u_case
 USE IO_module,                  ONLY : write_info



    !!----    Coefficients for calculating the X-ray scattering factors
    !!----        f(s) = Sum_{i=1,4} { a(i) exp(-b(i)*s^2) } + c
    !!----
    !!----    where s=sinTheta/Lambda

 implicit none
  INTEGER                  :: i,j

  call Set_Xray_Form ! definition des facteurs de diffusion

  do i=1, mendel_atom_nb
   do j=1, Num_Xray_Form
    if (u_case(mendel_atom_label(i)(1:))== u_case(Xray_form(j)%symb(1:))) exit

   end do

   call write_info('')
   call write_info('  Coefficients for calculating the X-ray scattering factors')
   call write_info( '       f(s) = Sum_{i=1,4} { a(i) exp(-b(i)*s^2) } + c')
   call write_info( '')
   call write_info( '    where s=sinTheta/Lambda')
   call write_info( '')
   call write_info( '   => X-ray scattering coeff. (A1,A2,A3,A4,B1,B2,B3,B4,C):')
   call write_info( '')
   write(message_text, '(2a)')      '   Element: ', u_case(Xray_form(j)%symb)
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,9F9.4)') '   coef:    ', Xray_form(j)%A(1:4), Xray_form(j)%B(1:4), Xray_form(j)%C
   call write_info(TRIM(message_text))
   call write_info('')

  END do
 RETURN
end subroutine write_X_rays_scatt_factors
!
!-------------------------------------------------------------------
subroutine write_X_rays_anomalous_factors
 USE Scattering_Chemical_Tables
 USE cryscal_module,             ONLY : message_text, mendel_atom_nb,  mendel_atom_label
 USE macros_module,              ONLY : u_case
 USE IO_module,                  ONLY : write_info
 
 implicit none
  INTEGER                  :: i,j, i1
  CHARACTER (LEN=4)        :: symb_car

  call Set_Delta_Fp_Fpp ! definition des facteurs de diffusion anomale fp et fpp

  do i=1, mendel_atom_nb
   i1 = INDEX(mendel_atom_label(i), '+')
   IF(i1>1) then
    symb_car = mendel_atom_label(i)(1:i1-1)
   else
    i1 = INDEX(mendel_atom_label(i), '-')
    IF(i1>1) then
     symb_car = mendel_atom_label(i)(1:i1-1)
    else
     symb_car = mendel_atom_label(i)(1:)
    endif
   endif

   do j=1, Num_Xray_Form
    if (u_case(symb_car(1:))== u_case(Anomalous_ScFac(j)%symb(1:))) exit

   end do

   call write_info('')
   call write_info('  . Values of Delta_Fp and Delta_Fpp anomalous dispersion coefficients')
   call write_info('    for 5 common radiation')
   call write_info('')
   call write_info('    Wavelenghts:       Cr        Fe        Cu        Mo        Ag')
   call write_info('         Lambda   2.28962   1.93597   1.54051   0.70926  0.556363')
   call write_info('')
   write(message_text, '(2a)')      '   Element: ', u_case(Anomalous_ScFac(j)%symb)
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,5(2x,F8.3))') '   Fp coef:    ', Anomalous_ScFac(j)%fp(1:5)
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,5(2x,F8.3))') '   Fpp coef:   ', Anomalous_ScFac(j)%fpp(1:5)
   call write_info(TRIM(message_text))

   call write_info('')

  END do

END subroutine write_X_rays_anomalous_factors
!-------------------------------------------------------------------

subroutine WRITE_data(input_string)
 USE IO_module,      ONLY : write_info
 USE cryscal_module, ONLY : message_text, pgf_data, pgf_file,                                              &
                            known_atomic_label, known_atomic_features, known_data_neutrons, known_data_x,  &
                            data_neutrons_PLOT, data_Xrays_PLOT, data_density_PLOT, data_weight_PLOT,      &
                            winplotr_exe
 USE atome_module,   ONLY : atom
 use atomic_data
 USE Neutrons_data
 USE Xrays_data

 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: input_string
  INTEGER                       :: i


  if(.not. known_atomic_label)  then
   call definition_atomic_label    ! %symbol, %name
   known_atomic_label = .true.
  endif
  if(.not. known_atomic_features) then
   call def_atomic_features()      ! %weight, %config_electr, %densite, %rayon
   known_atomic_features = .true.
  endif
  if(.not. known_data_neutrons)  then
   call definition_data_neutrons()
   known_data_neutrons = .true.
  endif
  if(.not. known_data_X) then
   call definition_data_X
   known_data_X = .true.
  endif

  select case (input_string)
      case ('neutrons')
       call write_info('     Element        coherent scattering length (in 10-12 cm)')
       call write_info('')
       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,10x,3F15.5)') i,  atom(i)%symbol, atom(i)%bcoh, atom(i)%sedinc, atom(i)%sea
        call write_info(TRIM(message_text))
       END do

       IF(data_neutrons_PLOT) then
        pgf_file%name = 'neutrons_bcoh.pgf'
        do i=1,96
         pgf_data(i)%X = REAL(i)
         WRITE(pgf_data(i)%string, '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%bcoh
        END do
        call create_PGF_file(TRIM(pgf_file%name), pgf_data%X, atom%bcoh, pgf_data%string, 96, "bcoh")
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         call system('winplotr '//trim(pgf_file%name), .true. )
        endif

        pgf_file%name = 'neutrons_sedinc.pgf'
        do i=1,96
         WRITE(pgf_data(i)%string, '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%sedinc
        end do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%sedinc, pgf_data%string, 96, "sedinc")
        IF(LEN_TRIM(winplotr_exe) /=0) call system('winplotr '//trim(pgf_file%name), .true. )

        pgf_file%name = 'neutrons_sea.pgf'
        do i=1,96
         WRITE(pgf_data(i)%string, '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%sea
        END do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%sea, pgf_data%string, 96, "sea")
        IF(LEN_TRIM(winplotr_exe) /=0) call system('winplotr '//trim(pgf_file%name), .true. )
       END if

      case ('xrays')
       call write_info('                    total interaction cross section (in barns)')
       call write_info('     Element        Mo              Cu            Co')
       call write_info('')

       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,10x,3F15.5)') i,  atom(i)%symbol, atom(i)%tics_mo_ka, atom(i)%tics_cu_ka, atom(i)%tics_co_ka
        call write_info(TRIM(message_text))
       end do

       IF(data_xrays_PLOT) then
        pgf_file%name = 'xrays_mo.pgf'
        do i=1,96
         pgf_data(i)%X = REAL(i)
         WRITE(pgf_data(i)%string, '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%tics_mo_ka
        end do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%tics_mo_ka, pgf_data%string, 96, "x_mo")
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         call system('winplotr '//trim(pgf_file%name), .true. )
        endif

        pgf_file%name = 'xrays_cu.pgf'
        do i=1,96
         WRITE(pgf_data(i)%string, '(2a,F10.5)') atom(i)%symbol, ': ', atom(i)%tics_cu_ka
        END do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%tics_cu_ka, pgf_data%string, 96, "x_cu")
        IF(LEN_TRIM(winplotr_exe) /= 0) call system('winplotr '//trim(pgf_file%name), .true. )

        pgf_file%name = 'xrays_co.pgf'
        do i=1,96
         WRITE(pgf_data(i)%string, '(2a,F10.5)') atom(i)%symbol, ': ', atom(i)%tics_co_ka
        END do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%tics_co_ka, pgf_data%string, 96, "x_co")
        IF(LEN_TRIM(winplotr_exe) /= 0)call system('winplotr '//trim(pgf_file%name), .true. )
       END if

      case('density')
       call write_info('     Element        density')
       call write_info('')

       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,10x,F15.5)') i,  atom(i)%symbol, atom(i)%densite
        call write_info(TRIM(message_text))
       end do

       IF(data_density_PLOT) then
        pgf_file%name = 'atomic_density.pgf'
        do i=1, 96
         pgf_data(i)%X = REAL(i)
         WRITE(pgf_data(i)%string, '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%densite
        end do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%densite, pgf_data%string, 96, "density")
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         call system('winplotr '//trim(pgf_file%name), .true. )
        endif
       ENDIF

      case('weight')
       call write_info('     Element        weight')
       call write_info('')

       do i=1, 96
        WRITE(message_text, '(2x,i3,5X,a,10x,F15.5)') i,  atom(i)%symbol, atom(i)%weight
        call write_info(TRIM(message_text))
       end do
       if (data_weight_PLOT) then
        pgf_file%name = 'atomic_weight.pgf'
        do i=1,96
         pgf_data(i)%X = REAL(i)
         WRITE(pgf_data(i)%string, '(2a,F15.5)') atom(i)%symbol, ': ', atom(i)%weight
        end do
        call create_PGF_file(TRIM(pgf_file%name),pgf_data%X, atom%weight, pgf_data%string, 96, "weight")
        IF(LEN_TRIM(winplotr_exe) == 0) then
         call write_info('')
         call write_info(' Warning: WinPLOTR is not installed ')
         call write_info('')
        else
         call system('winplotr '//trim(pgf_file%name), .true. )
        endif
       end if

      case default
  end select
 RETURN
end subroutine write_data
