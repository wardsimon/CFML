!     Last change:  TR   30 Jan 2007    4:38 pm

subroutine write_matrice()
 USE cryscalc_module,     ONLY : M => Mat, Mat_det, message_text, debug_proc
 USE IO_module,           ONLY : write_info
 USE MATRIX_list_module,  ONLY : matrix_num, transf_mat_text
 use matrix_module,       only : matrix_determinant_33

if(debug_proc%level_2)  call write_debug_proc_level(2, "write_matrice")

 IF(matrix_num /=0) then
  call write_info('')
  WRITE(message_text, '(2x,a)') TRIM(transf_mat_text(matrix_num))
  call write_info(TRIM(message_text))
 endif

 call write_info('')
 call write_info('  Transformation matrix for metrics:')
 write(message_text,'(10x,a,3F9.4,a)') '|', M(1,1), M(1,2), M(1,3), '|'
  call write_info(TRIM(message_text))
 write(message_text,'(10x,a,3F9.4,a)') '|', M(2,1), M(2,2), M(2,3), '|'
  call write_info(TRIM(message_text))
 write(message_text,'(10x,a,3F9.4,a)') '|', M(3,1), M(3,2), M(3,3), '|'
  call write_info(TRIM(message_text))


 Mat_det = matrix_determinant_33(M)
 WRITE(message_text, '(a,F8.3)') '  > Matrix determinant: ', Mat_det
 call write_info(TRIM(message_text))
 call write_info('')

 return
end subroutine write_matrice
!!--------------------------------------------------------------------------------------------!!

subroutine  transf_cell_parameters
 USE cryscalc_module, ONLY  : M => mat, Mat_det, unit_cell, pi, message_text,   &
                              WRITE_triclinic_transf, WRITE_monoclinic_transf,  &
                              WRITE_details, crystal_cell,                      &
                              WRITE_twin_hexa, WRITE_twin_pseudo_hexa, update_parameters, debug_proc
 USE IO_module,      ONLY   : write_info
 USE CFML_Crystal_Metrics, ONLY : Crystal_cell_type, Change_Setting_Cell


 implicit none
  REAL, DIMENSION(6)          :: tmp_cell
  type (Crystal_Cell_Type)    :: New_cell

  if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_cell_parameters")


  call Change_Setting_Cell(crystal_cell, M, New_cell)   ! le calcul des new esd n'est pas effectue dans la routine

  unit_cell%new_param(1:3) = New_cell%cell(1:3)
  unit_cell%new_param(4:6) = New_cell%ang(1:3)



 if(WRITE_details) then
  call write_info('')
  call write_info('  New cell parameters:')
  write(message_text,'(10x,6F9.4)') unit_cell%new_param(1:6)
  call write_info(TRIM(message_text))
  !write(message_text,'(10x,6F9.4)') unit_cell%new_param_ESD(1:6)
  !call write_info(TRIM(message_text))
  call write_info('')
 end if

 IF(unit_cell%volume < 0.1) call volume_calculation('')
 if(WRITE_details) then
  WRITE(message_text,'(a,F12.2,a)') '  New cell volume : ', abs(Mat_det) * unit_cell%volume, ' A3'
  call write_info(TRIM(message_text))
  call write_info('')
 end if

 tmp_cell = unit_cell%param

 if (      .not. WRITE_triclinic_transf  .and.  .not. WRITE_monoclinic_transf &
     .and. .not. WRITE_twin_hexa         .and.  .not. WRITE_twin_pseudo_hexa) then

    if(update_parameters)  then
     unit_cell%param = unit_cell%new_param
     unit_cell%param_ESD = unit_cell%new_param_ESD
    end if
 end if
 !unit_cell%param = tmp_cell
 !call volume_calculation('out')



 return
end subroutine  transf_cell_parameters

!!--------------------------------------------------------------------------------------------!!

subroutine  transf_cell_parameters_old
 USE cryscalc_module, ONLY  : M => mat, Mat_det, unit_cell, pi, message_text,   &
                              WRITE_triclinic_transf, WRITE_monoclinic_transf,  &
                              WRITE_details,                                    &
                              WRITE_twin_hexa, WRITE_twin_pseudo_hexa, update_parameters, debug_proc
 USE IO_module,      ONLY   : write_info
 USE matrix_module,  ONLY   : MV_product
 implicit none
  real                               :: cosAlfa1, cosBeta1, cosGamma1
  real                               :: cosAlfa2, cosBeta2, cosGamma2
  real                               :: a2b2, a2c2, b2c2
  REAL, DIMENSION(6)                 :: tmp_cell

  if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_cell_parameters")



  CosAlfa1  = COS(unit_cell%param(4)  * pi / 180.)
  Cosbeta1  = COS(unit_cell%param(5)  * pi / 180.)
  Cosgamma1 = COS(unit_cell%param(6)  * pi / 180.)



! parametres a, b, c
  unit_cell%new_param(1) =   M(1,1)**2. * unit_cell%param(1)**2. +  M(1,2) ** 2. * unit_cell%param(2) **2. &
                           + M(1,3)**2. * unit_cell%param(3)**2.  +                                        &
                2*M(1,1) * M(1,2) * unit_cell%param(1) * unit_cell%param(2) *  cosGamma1 +                 &
                2*M(1,1) * M(1,3) * unit_cell%param(1) * unit_cell%param(3) *  cosBeta1  +                 &
                2*M(1,2) * M(1,3) * unit_cell%param(2) * unit_cell%param(3) *  cosAlfa1


  unit_cell%new_param(2) =   M(2,1)**2. * unit_cell%param(1)**2. +  M(2,2) ** 2. * unit_cell%param(2) **2. &
                           + M(2,3)**2. * unit_cell%param(3)**2.  +                                        &
                2*M(2,1) * M(2,2) * unit_cell%param(1) * unit_cell%param(2) *  cosGamma1 +                 &
                2*M(2,1) * M(2,3) * unit_cell%param(1) * unit_cell%param(3) *  cosBeta1  +                 &
                2*M(2,2) * M(2,3) * unit_cell%param(2) * unit_cell%param(3) *  cosAlfa1

  unit_cell%new_param(3) =   M(3,1)**2. * unit_cell%param(1)**2. +  M(3,2) ** 2. * unit_cell%param(2) **2. &
                           + M(3,3)**2. * unit_cell%param(3)**2.  +                                        &
                2*M(3,1) * M(3,2) * unit_cell%param(1) * unit_cell%param(2) *  cosGamma1 +                 &
                2*M(3,1) * M(3,3) * unit_cell%param(1) * unit_cell%param(3) *  cosBeta1  +                 &
                2*M(3,2) * M(3,3) * unit_cell%param(2) * unit_cell%param(3) *  cosAlfa1

  unit_cell%new_param(1:3) = sqrt(unit_cell%new_param(1:3))

 ! esd des parametres
  unit_cell%new_param_ESD(1) =   M(1,1)**2. * unit_cell%param_ESD(1)**2. +  M(1,2) ** 2. * unit_cell%param_ESD(2) **2. &
                               + M(1,3)**2. * unit_cell%param_ESD(3)**2.  +                                            &
                2*M(1,1) * M(1,2) * unit_cell%param_ESD(1) * unit_cell%param_ESD(2) *  cosGamma1 +                     &
                2*M(1,1) * M(1,3) * unit_cell%param_ESD(1) * unit_cell%param_ESD(3) *  cosBeta1  +                     &
                2*M(1,2) * M(1,3) * unit_cell%param_ESD(2) * unit_cell%param_ESD(3) *  cosAlfa1

  unit_cell%new_param_ESD(2) =   M(2,1)**2. * unit_cell%param_ESD(1)**2. +  M(2,2) ** 2. * unit_cell%param_ESD(2) **2. &
                               + M(2,3)**2. * unit_cell%param_ESD(3)**2.  +                                            &
                2*M(2,1) * M(2,2) * unit_cell%param_ESD(1) * unit_cell%param_ESD(2) *  cosGamma1 +                     &
                2*M(2,1) * M(2,3) * unit_cell%param_ESD(1) * unit_cell%param_ESD(3) *  cosBeta1  +                     &
                2*M(2,2) * M(2,3) * unit_cell%param_ESD(2) * unit_cell%param_ESD(3) *  cosAlfa1

  unit_cell%new_param_ESD(3) =   M(3,1)**2. * unit_cell%param_ESD(1)**2. +  M(3,2) ** 2. * unit_cell%param_ESD(2) **2. &
                               + M(3,3)**2. * unit_cell%param_ESD(3)**2.  +                                            &
                2*M(3,1) * M(3,2) * unit_cell%param_ESD(1) * unit_cell%param_ESD(2) *  cosGamma1 +                     &
                2*M(3,1) * M(3,3) * unit_cell%param_ESD(1) * unit_cell%param_ESD(3) *  cosBeta1  +                     &
                2*M(3,2) * M(3,3) * unit_cell%param_ESD(2) * unit_cell%param_ESD(3) *  cosAlfa1

  unit_cell%new_param_ESD(1:3) = sqrt(ABS(unit_cell%new_param_ESD(1:3)))


! angles

  a2b2  =    M(1,1) * M(2,1) * unit_cell%param(1) ** 2. +  M(1,2) * M(2,2) * unit_cell%param(2) ** 2.  +   &
             M(1,3) * M(2,3) * unit_cell%param(3) ** 2.  +                                                 &
             unit_cell%param(1) * unit_cell%param(2) * cosGamma1 * (M(1,1) * M(2,2) + M(1,2) * M(2,1))  +  &
             unit_cell%param(1) * unit_cell%param(3) * cosbeta1  * (M(1,1) * M(2,3) + M(1,3) * M(2,1))  +  &
             unit_cell%param(2) * unit_cell%param(3) * cosAlfa1  * (M(1,2) * M(2,3) + M(1,3) * M(2,2))


  a2c2  =    M(1,1) * M(3,1) * unit_cell%param(1) ** 2. +  M(1,2) * M(3,2) * unit_cell%param(2) ** 2.  +   &
             M(1,3) * M(3,3) * unit_cell%param(3) ** 2.  +                                                 &
             unit_cell%param(1) * unit_cell%param(2) * cosGamma1 * (M(1,1) * M(3,2) + M(1,2) * M(3,1))  +  &
             unit_cell%param(1) * unit_cell%param(3) * cosbeta1  * (M(1,1) * M(3,3) + M(1,3) * M(3,1))  +  &
             unit_cell%param(2) * unit_cell%param(3) * cosAlfa1  * (M(1,2) * M(3,3) + M(1,3) * M(3,2))


  b2c2  =    M(2,1) * M(3,1) * unit_cell%param(1) ** 2. +  M(2,2) * M(3,2) * unit_cell%param(2) ** 2.  +   &
             M(2,3) * M(3,3) * unit_cell%param(3) ** 2.  +                                                 &
             unit_cell%param(1) * unit_cell%param(2) * cosGamma1 * (M(2,1) * M(3,2) + M(2,2) * M(3,1))  +  &
             unit_cell%param(1) * unit_cell%param(3) * cosbeta1  * (M(2,1) * M(3,3) + M(2,3) * M(3,1))  +  &
             unit_cell%param(2) * unit_cell%param(3) * cosAlfa1  * (M(2,2) * M(3,3) + M(2,3) * M(3,2))

  cosGamma2 = a2b2 / (unit_cell%new_param(1)*unit_cell%new_param(2))
  cosBeta2  = a2c2 / (unit_cell%new_param(1)*unit_cell%new_param(3))
  cosAlfa2  = b2c2 / (unit_cell%new_param(2)*unit_cell%new_param(3))
  unit_cell%new_param(4) = ACOS(cosAlfa2)  * 180./pi
  unit_cell%new_param(5) = ACOS(cosBeta2)  * 180./pi
  unit_cell%new_param(6) = ACOS(cosGamma2) * 180./pi

 ! do i= 4,6
 !  if(unit_cell%new_param(i) < 90. )  unit_cell%new_param(i) = 180. -  unit_cell%new_param(i)
 ! end do

! esd des angles

  CosAlfa1  = COS(unit_cell%param_ESD(4)  * pi / 180.)
  Cosbeta1  = COS(unit_cell%param_ESD(5)  * pi / 180.)
  Cosgamma1 = COS(unit_cell%param_ESD(6)  * pi / 180.)


  a2b2  =    M(1,1) * M(2,1) * unit_cell%param_ESD(1) ** 2. +  M(1,2) * M(2,2) * unit_cell%param_ESD(2) ** 2.  +   &
             M(1,3) * M(2,3) * unit_cell%param_ESD(3) ** 2. +                                          &
             unit_cell%param_ESD(1) * unit_cell%param_ESD(2) * cosGamma1 * (M(1,1) * M(2,2) + M(1,2) * M(2,1))  +  &
             unit_cell%param_ESD(1) * unit_cell%param_ESD(3) * cosbeta1  * (M(1,1) * M(2,3) + M(1,3) * M(2,1))  +  &
             unit_cell%param_ESD(2) * unit_cell%param_ESD(3) * cosAlfa1  * (M(1,2) * M(2,3) + M(1,3) * M(2,2))


  a2c2  =    M(1,1) * M(3,1) * unit_cell%param_ESD(1) ** 2. +  M(1,2) * M(3,2) * unit_cell%param_ESD(2) ** 2.  +   &
             M(1,3) * M(3,3) * unit_cell%param_ESD(3) ** 2.  +                                          &
             unit_cell%param_ESD(1) * unit_cell%param_ESD(2) * cosGamma1 * (M(1,1) * M(3,2) + M(1,2) * M(3,1))  +  &
             unit_cell%param_ESD(1) * unit_cell%param_ESD(3) * cosbeta1  * (M(1,1) * M(3,3) + M(1,3) * M(3,1))  +  &
             unit_cell%param_ESD(2) * unit_cell%param_ESD(3) * cosAlfa1  * (M(1,3) * M(3,2) + M(1,3) * M(3,2))


  b2c2  =    M(2,1) * M(3,1) * unit_cell%param_ESD(1) ** 2. +  M(2,2) * M(3,2) * unit_cell%param_ESD(2) ** 2.  +   &
             M(2,3) * M(3,3) * unit_cell%param_ESD(3) ** 2.  +                                          &
             unit_cell%param_ESD(1) * unit_cell%param_ESD(2) * cosGamma1 * (M(2,1) * M(3,2) + M(2,2) * M(3,1))  +  &
             unit_cell%param_ESD(1) * unit_cell%param_ESD(3) * cosbeta1  * (M(2,1) * M(3,3) + M(2,3) * M(3,1))  +  &
             unit_cell%param_ESD(2) * unit_cell%param_ESD(3) * cosAlfa1  * (M(2,3) * M(3,2) + M(2,3) * M(3,2))

  cosGamma2 = a2b2 / (unit_cell%param_ESD(1)*unit_cell%param_ESD(2))
  cosBeta2  = a2c2 / (unit_cell%param_ESD(1)*unit_cell%param_ESD(3))
  cosAlfa2  = b2c2 / (unit_cell%param_ESD(2)*unit_cell%param_ESD(3))

  unit_cell%new_param_ESD(4) = ACOS(cosAlfa2)  * 180./pi
  unit_cell%new_param_ESD(5) = ACOS(cosBeta2)  * 180./pi
  unit_cell%new_param_ESD(6) = ACOS(cosGamma2) * 180./pi

  !if(abs(Mat_det) == 1) then
  ! unit_cell%new_param_ESD(4:6) = MV_product(unit_cell%param_ESD(4:6),  M)
  ! unit_cell%new_param_ESD(4:6) = ABS(unit_cell%new_param_ESD(4:6))
  !else
  ! unit_cell%new_param_ESD(4:6) = 0.
  !end if


 if(WRITE_details) then
  call write_info('')
  call write_info('  New cell parameters:')
  write(message_text,'(10x,6F9.4)') unit_cell%new_param(1:6)
  call write_info(TRIM(message_text))
  write(message_text,'(10x,6F9.4)') unit_cell%new_param_ESD(1:6)
  call write_info(TRIM(message_text))
  call write_info('')
 end if

 IF(unit_cell%volume < 0.1) call volume_calculation('')
 if(WRITE_details) then
  WRITE(message_text,'(a,F12.2,a)') '  New cell volume : ', abs(Mat_det) * unit_cell%volume, ' A3'
  call write_info(TRIM(message_text))
  call write_info('')
 end if

 tmp_cell = unit_cell%param

 if (      .not. WRITE_triclinic_transf  .and.  .not. WRITE_monoclinic_transf &
     .and. .not. WRITE_twin_hexa         .and.  .not. WRITE_twin_pseudo_hexa) then

    if(update_parameters)  then
     unit_cell%param = unit_cell%new_param
     unit_cell%param_ESD = unit_cell%new_param_ESD
    end if
 end if
 !unit_cell%param = tmp_cell
 !call volume_calculation('out')



 return
end subroutine  transf_cell_parameters_old


!----------------------------------------------------------------------

subroutine transf_HKL()
 USE cryscalc_module, ONLY      : M => Mat, nb_hkl, H, message_text, debug_proc
 USE IO_module,       ONLY      : write_info
 USE matrix_module,   ONLY      : MV_product
 implicit none
  !real                         :: new_h, new_k, new_l
  REAL, DIMENSION(3)           :: new_H, Hr
  INTEGER                      :: i

 if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_hkl")

 do i=1, nb_hkl

  !new_H = MV_product(real(H(:,i)),  M)  ! >> warning with IFC
  Hr = real(H(:,i))
  new_H = MV_product(Hr,  M)

  WRITE(message_text,'(a,3F6.2,10x,a,3F6.2)') '  h1 k1 l1 : ', H(1,i), H(2,i), H(3,i) , &
                                              '  h2 k2 l2 : ', new_H(1), new_H(2), new_H(3)
  call write_info(TRIM(message_text))

 END do

 call write_info('')
 return
end subroutine transf_HKL

!----------------------------------------------------------------------

subroutine transf_coord()
 USE cryscalc_module, ONLY      : Mit, nb_atom, atom_coord, new_atom_coord, atom_label, atom_typ, message_text, update_parameters, &
                                  debug_proc
 USE IO_module,       ONLY      : write_info
 use matrix_module,   ONLY      : MV_product
 implicit none
  !real                         :: new_x, new_y, new_z
  !REAL, DIMENSION(3)           :: new_coord
  INTEGER                      :: i

 if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_coord")

 !call inversion_matrice()
 call get_matrice_inverse_transposee()
 call write_atomic_matrix()
 IF(nb_atom /=0) then
  call write_info('  New atomic coordinates: ')
  call write_info('')
  do i=1, nb_atom
   new_atom_coord(:,i) = MV_product(atom_coord(:, i), Mit)


   WRITE(message_text,'(a,2a6,3(1x,F10.6))') '  ATOM ', trim(atom_label(i)),trim(atom_typ(i)),  new_atom_coord(:,i)
   call write_info(TRIM(message_text))
   if(update_parameters) atom_coord(:,i) = new_atom_coord(:,i)
  END do
  call write_info('')
 endif

 return
end subroutine transf_coord
!----------------------------------------------------------------------
subroutine write_atomic_matrix()
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : Mit, message_text, debug_proc

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_atomic_matrix")

 call write_info('')
 call write_info('  Transformation matrix for atomic coordinates:')
 write(message_text,'(10x,a,3F9.4,a)') '|',Mit(1,1), Mit(1,2), Mit(1,3),'|'
  call write_info(TRIM(message_text))
 write(message_text,'(10x,a,3F9.4,a)') '|',Mit(2,1), Mit(2,2), Mit(2,3), '|'
  call write_info(TRIM(message_text))
 write(message_text,'(10x,a,3F9.4,a)') '|',Mit(3,1), Mit(3,2), Mit(3,3), '|'
  call write_info(TRIM(message_text))
 call write_info('')
 WRITE(message_text,'(5x,a,3(F9.4,a))') ' x2 = ',Mit(1,1),' * x1  + ',Mit(1,2), ' * y1  + ',Mit(1,3), ' * z1   '
  call write_info(TRIM(message_text))
 WRITE(message_text,'(5x,a,3(F9.4,a))') ' y2 = ',Mit(2,1),' * x1  + ',Mit(2,2), ' * y1  + ',Mit(2,3), ' * z1   '
  call write_info(TRIM(message_text))
 WRITE(message_text,'(5x,a,3(F9.4,a))') ' z2 = ',Mit(3,1),' * x1  + ',Mit(3,2), ' * y1  + ',Mit(3,3), ' * z1   '
  call write_info(TRIM(message_text))
 call write_info('')


end subroutine write_atomic_matrix
!----------------------------------------------------------------------

subroutine transl_coord()
 USE cryscalc_module, ONLY : translat, nb_atom, atom_coord, atom_label, atom_typ, message_text, update_parameters, debug_proc
 USE IO_module,       ONLY : write_info
 implicit none

  REAL, DIMENSION(3)           :: new_coord
  INTEGER                      :: i

  if(debug_proc%level_2)  call write_debug_proc_level(2, "transl_coord")

  call write_info('')
  WRITE(message_text, '(a,3(1x,F10.6),5x,I3)') '  Unit cell translation:', translat(1:3), INT(translat(4))
  call write_info(TRIM(message_text))
  call write_info(' ')
  call write_info('  ... new coordinates :')
  call write_info(' ')

  do i=1, nb_atom
   new_coord(1:3) = translat(4)*atom_coord(1:3, i) + translat(1:3)
   WRITE(message_text,'(a,2a6,3(1x,F10.6))') '  ATOM ', trim(atom_label(i)),trim(atom_typ(i)),  new_coord
   call write_info(TRIM(message_text))
   if(update_parameters) atom_coord(:, i) = new_coord
  END do

  call write_info('')

  return
end subroutine transl_coord

!----------------------------------------------------------------------

subroutine inside_unit_cell
 USE cryscalc_module, ONLY : nb_atom, atom_coord, atom_label, atom_typ, message_text, debug_proc
 USE IO_module,       ONLY : write_info
 use matrix_module,   ONLY : MV_product
 implicit none
  INTEGER                       :: i, j

  if(debug_proc%level_2)  call write_debug_proc_level(2, "inside_unit_cell")

  call write_info('')
  call write_info(' > List of atoms with atomic coordinates inside the unit cell')
  call write_info(' -------------------------------------------------------------')
  call write_info('')

  do i=1, nb_atom
   do j=1,3
    if(atom_coord(j,i) < 0.) then
     atom_coord(j,i) = atom_coord(j,i) + 1.
    elseif(atom_coord(j,i) > 1.) then
     atom_coord(j,i) = atom_coord(j,i) - 1.
    end if
   END do
   WRITE(message_text,'(a,2a6,3(1x,F9.6))') '  ATOM ', trim(atom_label(i)),trim(atom_typ(i)),  atom_coord(:,i)
   call write_info(TRIM(message_text))
  END do

  return
end subroutine inside_unit_cell

!!-------------------------------------------------------------------------------------
! put atomic coordinated inside the current unit cell (between 0. and 1.)
subroutine inside_coord(coord)
 implicit none
  real, dimension(3), intent(inout)  :: coord
  integer                            :: i

  do i=1, 3
   if(coord(i) < 0.) then
    coord(i) = coord(i) + 1
   elseif(coord(i) > 1.) then
    coord(i) = coord(i) - 1.
   end if
  end do

 return
end subroutine inside_coord



!!----------------------------------------------------------------------------------------------

subroutine get_matrice_inverse_transposee()
 USE cryscalc_module, ONLY : M => Mat, Mi => Mat_inv, message_text, Mit, Mat_det, debug_proc
 USE matrix_module
 USE IO_module,      ONLY : write_info
 implicit none
  real, dimension(3,3)              :: Mp   ! matrice des complements algebriques
  REAL, DIMENSION(3,3)              :: Mpt  ! matrice transposee

  if(debug_proc%level_2)  call write_debug_proc_level(2, "get_matrice_inverse_transposee")

! l'inversion d'une matrice 3*3 est une matrice 3*3
! M-1 * M = 1
! M-1 = 1/Det(M) * M't
!   det(M): determinant de la matrice M
!   M'    : matrice des compléments algébriques dont les éléments sont affectés
!           du signe (-1)**i+j

 ! calcul du determinant de la matrice suivant la premiere ligne
   Mat_det = matrix_determinant_33(M)
   if (abs(Mat_det) < 0.0001) then
    call write_info ('')
    call write_info ('  !!! matrix determinant equal to zero !!!')
    return
   endif

 ! calcul de M':
   call get_matrix_complements_algebriques_33(M, Mp)

 ! M' transposée: Mpt(i,j) = Mp(j,i)
   call get_matrix_transposee_33(Mp, Mpt)

 ! matrice inverse: M-1(i,j) = M'T(i,j)/det
   Mi(:,:) =  Mpt(:,:)/Mat_det

 ! transposee de la matrice inverse
   call get_matrix_transposee_33(Mi, Mpt)


 !call write_info('')
 !call write_info('  Transformation matrix for atomic coordinates:')
 !write(message_text,'(10x,a,3F9.4,a)') '|',Mpt(1,1), Mpt(1,2), Mpt(1,3),'|'
 ! call write_info(TRIM(message_text))
 !write(message_text,'(10x,a,3F9.4,a)') '|',Mpt(2,1), Mpt(2,2), Mpt(2,3), '|'
 ! call write_info(TRIM(message_text))
 !write(message_text,'(10x,a,3F9.4,a)') '|',Mpt(3,1), Mpt(3,2), Mpt(3,3), '|'
 ! call write_info(TRIM(message_text))
 !call write_info('')
 !WRITE(message_text,'(5x,a,3(F9.4,a))') ' x2 = ',Mpt(1,1),' * x1  + ',Mpt(1,2), ' * y1  + ',Mpt(1,3), ' * z1   '
 ! call write_info(TRIM(message_text))
 !WRITE(message_text,'(5x,a,3(F9.4,a))') ' y2 = ',Mpt(2,1),' * x1  + ',Mpt(2,2), ' * y1  + ',Mpt(2,3), ' * z1   '
 ! call write_info(TRIM(message_text))
 !WRITE(message_text,'(5x,a,3(F9.4,a))') ' z2 = ',Mpt(3,1),' * x1  + ',Mpt(3,2), ' * y1  + ',Mpt(3,3), ' * z1   '
 ! call write_info(TRIM(message_text))
 !call write_info('')

  Mit = Mpt

 RETURN
end subroutine get_matrice_inverse_transposee


!!--------------------------------------------------------------------------------------
subroutine  transf_HKL_file()
 use cryscalc_module
 USE HKL_module,    ONLY  : HKL_file
 USE IO_module,     ONLY  : write_info
 USE matrix_module, ONLY  : MV_product
  implicit none
   integer                           :: i, i_error
   INTEGER, DIMENSION(3)             :: h_
   real,    DIMENSION(3)             :: hr_
   integer                           :: code_
   real                              :: F2_, sig_F2_
   real                              :: cos1,  cos2, cos3, cos4, cos5, cos6
   REAL, DIMENSION(3)                :: new_H
   logical                           :: new_indices_entiers
   REAL, parameter                   :: EPS=0.001
   integer                           :: n

  if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_HKl_file")

   !lecture fichier .HKL format SHELX
   ! + creation fichier avec sintheta/l dans la 1ere colonne

   open (unit=HKL_unit, file=trim(HKL_file%HKL), iostat=i_error)
   if(i_error/=0) then
    call write_info('')
    call write_info('   >>> Error opening the '// trim(HKL_file%HKL)// ' file <<<')
    call write_info('')
    return
   endif
   i=index(HKL_file%NAME, '.', back=.true.)
   open(unit=23, file=HKL_file%NAME(1:i-1)//'_trans.HKL', iostat=i_error)
   if(i_error/=0) then
    call write_info('')
    call write_info('   >>> Error opening the '// HKL_file%NAME(1:i-1)//'_trans.HKL'// ' file <<<')
    call write_info('')
    return
   endif

   if(.not. Mat_integer) then
    open(unit=24, file=HKL_file%NAME(1:i-1)//'_ini_1.HKL', iostat=i_error)
    open(unit=25, file=HKL_file%NAME(1:i-1)//'_new_2.HKL', iostat=i_error)
    open(unit=26, file=HKL_file%NAME(1:i-1)//'_ini_2.HKL', iostat=i_error)
   endif


   !call write_info('')
   !call write_info('  ... Reading data in '// trim(HKL_file%NAME)// ' ...')
   !call write_info('')


   do
    READ(HKL_unit,'(3I4,2F8.2,I4,6F8.5)', iostat=i_error) h_(1), h_(2), h_(3), F2_, sig_F2_, code_,  &
                                                          cos1,  cos2, cos3, cos4, cos5, cos6
    if (i_error < 0) exit ! fin du fichier atteinte
    !new_H = MV_product(real(h_), Mat)  ! >> warning with IFC
    hr_ = real(h_)
    new_H = MV_product(hr_, Mat)

    new_indices_entiers = .false.
    n=0
    do i=1, 3
     !if (ABS(H_(i) - floor(H_(i)+.5)) < eps ) n=n+1
     if (ABS(new_H(i) - floor(new_H(i)+.5)) < eps ) n=n+1
    end do
    if(n==3) new_indices_entiers = .true.


!    IF(Mat_integer) then
!     write(3, '(3I4,2F8.2,I4)' ) floor(new_H(1)+.5),floor(new_H(2)+.5),floor(new_H(3)+.5), F2_, sig_F2_, code_
!     write(4, '(3I4,2F8.2,I4)' ) H_(1), H_(2), H_(3), F2_, sig_F2_, code_
!    else
!     write(3, '(3(1x,F7.3),1x,2F8.2,I4)') new_H(1), new_H(2), new_H(3), F2_, sig_F2_, code_
!     write(4, '(3I4,2F8.2,I4)' ) H_(1), H_(2), H_(3), F2_, sig_F2_, code_
!    endif

     if(new_indices_entiers) then
      write(23, '(3I4,2F8.2,I4)' ) floor(new_H(1)+.5),floor(new_H(2)+.5),floor(new_H(3)+.5), F2_, sig_F2_, code_
      if(.not. Mat_integer)  write(24, '(3I4,2F8.2,I4)' ) H_(1), H_(2), H_(3), F2_, sig_F2_, code_
     else
      write(25, '(3(1x,F7.3),1x,2F8.2,I4)') new_H(1), new_H(2), new_H(3), F2_, sig_F2_, code_
      write(26, '(3I4,2F8.2,I4)' ) H_(1), H_(2), H_(3), F2_, sig_F2_, code_
     endif

   end do

   i=index(HKL_file%NAME, '.', back=.true.)
   write(HKL_file%transf, '(2a)') HKL_file%NAME(1:i-1), '_trans.HKL'
   call write_info(' ')
   call write_info('  >> HKL file after matrix transformation: '// trim(HKL_file%transf))
   call write_info(' ')
   if(.not. Mat_integer) then
    call write_info('  >> HKL file after matrix transformation (reals indices): '// HKL_file%NAME(1:i-1)//'_new_2.HKL')
    call write_info(' ')
    call write_info('  >> Initial HKL file with "even" indices: '// HKL_file%NAME(1:i-1)//'_ini_1.HKL')
    call write_info(' ')
    call write_info('  >> Initial HKL file with "odd" indices : '// HKL_file%NAME(1:i-1)//'_ini_2.HKL')
    call write_info(' ')
   endif

   CLOSE(UNIT=23)
   close(unit=24)
   close(unit=25)
   close(unit=26)

  return
end subroutine transf_HKL_file


!!----------------------------------------------------------------------------------------------
subroutine   transf_triclinic
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : keyword_CELL, Mat, Mat_det, debug_proc
 USE Matrix_module,   ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none
  integer                :: mat_numor

 if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_triclinic")

  call write_info(' ')
  call write_info(' --------- Triclinic  transformations -----------')


  do mat_numor = 3, 7
   matrix_num = mat_numor
   call write_info('')
   call write_info('  *** M '//trim(transf_mat_text(mat_numor))//' ***')
   Mat(:,:) = transf_mat(:,:, mat_numor)
   call write_matrice()
   Mat_det = matrix_determinant_33(Mat)
   if (keyword_CELL) call transf_cell_parameters
  end do



 return
end subroutine   transf_triclinic


!!-----------------------------------------------------------------
subroutine   transf_monoclinic
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : keyword_CELL, Mat, Mat_det, unit_cell, message_text, WRITE_details, debug_proc
 USE Matrix_module,   ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none
  integer                        :: i, mat_numor
  real, dimension(6,MAT_mono_nb) :: tcp

 if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_monoclinic")

  if(WRITE_details) then
   call write_info(' ')
   call write_info(' --------- Monoclinic transformations -----------')
  end if

  !do mat_numor = 27, 28
  do mat_numor = MAT_mono_ini, MAT_mono_fin
   matrix_num = mat_numor
   if(WRITE_details) call write_info('')
   !call write_info('  *** M '//trim(transf_mat_text(mat_numor))//' ***')
   Mat(:,:) = transf_mat(:,:, mat_numor)
   if(WRITE_details) call write_matrice()
   Mat_det = matrix_determinant_33(Mat)
   if (keyword_CELL) call transf_cell_parameters
   do i=1, 6
    tcp(i, mat_numor-MAT_mono_ini+1) = unit_cell%new_param(i)
   end do
  end do

  if(WRITE_details) then
  call write_info('')
  call write_info('   . examples:')
  call write_info('')
  call write_info('')
  call write_info('      P 21/n  ==> P 21/c :   matrix M #29, #33')
  call write_info('      P 21/a  ==> P 21/c :   matrix M #30, #32')
  call write_info('')
  call write_info('      I 2/a   ==> C 2/c  :   matrix M #27')
  call write_info('      C 2/c   ==> I 2/a  :   matrix M #28')
  call write_info('')
  end if


  call write_info('     Summary of monoclinic settings:')
  write(message_text, '(8x,a, 6F10.4,5x,a)') '. cell #1: ', unit_cell%param(:), '[Mat. #1,  #30, #32]'
  call write_info(trim(message_text))
  write(message_text, '(8x,a, 6F10.4,5x,a)') '. cell #2: ', tcp(:,1), '[Mat. #27, #29, #33]'
  call write_info(trim(message_text))
  write(message_text, '(8x,a, 6F10.4,5x,a)') '. cell #3: ', tcp(:,2), '[Mat. #28, #31, #34]'
  call write_info(trim(message_text))


 return
end subroutine   transf_monoclinic

!!-------------------------------------------------------------------------------------
subroutine   twin_hexa
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : keyword_CELL, Mat, Mat_det, debug_proc
 USE Matrix_module,   ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none
  integer                :: mat_numor

  if(debug_proc%level_2)  call write_debug_proc_level(2, "twin_hexa")

  call write_info(' ')
  call write_info(' --------- Hexagonal transformations -----------')


  do mat_numor = MAT_hexa_ini, MAT_hexa_fin
   matrix_num = mat_numor
   call write_info('')
   call write_info('  *** M '//trim(transf_mat_text(mat_numor))//' ***')
   Mat(:,:) = transf_mat(:,:, mat_numor)
   call write_matrice()
   Mat_det = matrix_determinant_33(Mat)
   if (keyword_CELL) call transf_cell_parameters
  end do



 return
end subroutine   twin_hexa


!!-------------------------------------------------------------------------------------------------
subroutine   twin_pseudo_hexa
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : keyword_CELL, Mat, Mat_det, debug_proc
 USE Matrix_module,   ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none
  integer       :: mat_numor

 if(debug_proc%level_2)  call write_debug_proc_level(2, "twin_pseudo_hexa")

  call write_info(' ')
  call write_info(' --------- Pseudo-hexagonal twinning in a monoclinic cell -----------')

  do mat_numor = MAT_mono_ini, MAT_mono_fin
   matrix_num = mat_numor
   call write_info('')
   call write_info('  *** M '//trim(transf_mat_text(mat_numor))//' ***')
   Mat(:,:) = transf_mat(:,:, mat_numor)
   call write_matrice()
   Mat_det = matrix_determinant_33(Mat)
   if (keyword_CELL) call transf_cell_parameters
  end do

 return
end subroutine   twin_pseudo_hexa


!!----------------------------------------------------------------------------------------
subroutine   permutation_abc
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : keyword_CELL, Mat, Mat_det, Mat_det, debug_proc
 use matrix_module,   only : matrix_determinant_33
 USE Matrix_list_module

 implicit none
  integer               :: mat_numor

 if(debug_proc%level_2)  call write_debug_proc_level(2, "permutation_abc")

  call write_info(' ')
  call write_info(' --------- Transformations -----------')

  do mat_numor = MAT_ortho_ini, MAT_ortho_fin
   matrix_num = mat_numor
   call write_info('')
   call write_info('  *** M '//trim(transf_mat_text(mat_numor))//' ***')
   Mat(:,:) = transf_mat(:,:, mat_numor)
   call write_matrice()
   Mat_det = matrix_determinant_33(Mat)
   if (keyword_CELL) call transf_cell_parameters
  end do

 return
end subroutine   permutation_abc



!-----------------------------------------------------------------
subroutine   transf_rhomb_to_hex
 USE IO_module,              ONLY : write_info
 USE cryscalc_module,        ONLY : keyword_CELL, Mat, Mat_det, debug_proc
 USE Matrix_module,          ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none

  if(debug_proc%level_2)  call write_debug_proc_level(2, "trasnf_rhomb_to_hex")

  call write_info(' ')
  call write_info('  ---------Special transformations -----------')
  call write_info(' ')
  call write_info('                     Tr_h')
  call write_info('  Rhomboedral cell  =====> Hexagonal cell')
  call write_info('')
  call write_info(' ')
  call write_info(' ')
  call write_info('               (  1 -1  0)')
  call write_info('  with:  Tr_h =(  0  1 -1)')
  call write_info('               (  1  1  1)')
  call write_info(' ')
  call write_info('         (equiv keyword : MAT#8)')
  call write_info(' ')

  Mat(:,:) = transf_mat(:,:, 8)
  Mat_det = matrix_determinant_33(Mat)
  IF(keyword_CELL) call transf_cell_parameters


 return
end subroutine   transf_rhomb_to_hex
!-----------------------------------------------------------------
subroutine   transf_hex_to_rhomb
 USE cryscalc_module,        ONLY : keyword_CELL, Mat, Mat_det, debug_proc
 USE IO_module,              ONLY : write_info
 USE Matrix_module,          ONLY : matrix_determinant_33
 USE MATRIX_list_module

 !USE MATRIX_list_module
 implicit none

  if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_hex_to_rhomb")

  call write_info(' ')
  call write_info('  ---------Special transformations -----------')
  call write_info(' ')
  call write_info(' ')
  call write_info('                   Th_r')
  call write_info('  Hexagonal cell  =====>  Rhomboedral cell')
  call write_info('')
  call write_info(' ')
  call write_info(' ')
  call write_info('               ( 2/3  1/3  1/3)')
  call write_info('  with:  Th_r =(-1/3  1/3  1/3)')
  call write_info('               (-1/3 -2/3  1/3)')
  call write_info(' ')
  call write_info('         (equiv keyword : MAT#9)')
  call write_info(' ')

  Mat(:,:) = transf_mat(:,:, 9)
  Mat_det = matrix_determinant_33(Mat)
  IF(keyword_CELL) call transf_cell_parameters


 return
end subroutine   transf_hex_to_rhomb
!-----------------------------------------------------------------

subroutine write_list_matrice()
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY : message_text, Mat, Mat_det, Mit, debug_proc
 USE matrix_module,   ONLY : matrix_determinant_33
 USE MATRIX_list_module
 implicit none
  INTEGER :: i
  character (len=24) :: fmt_

  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_list_matrice")

  call write_info('')
  call write_info('  > List of transformation matrices implemented in CRYSCALC: ')
  call write_info('')

  do i=1, Max_mat_nb + user_mat_nb
   Mat(:,:) =  transf_mat(:,:,i)
   call get_matrice_inverse_transposee()

   call write_info('')
   WRITE(message_text, '(2x,a)') TRIM(transf_mat_text(i))
   call write_info(TRIM(message_text))
   call write_info('')
   call write_info('        Direct matrix:                       Inverse matrix:')

   fmt_= '(2x,3F10.5,5x,3F10.5)'
   WRITE(message_text, fmt=trim(fmt_)) transf_mat(1,1,i), transf_mat(1,2,i), transf_mat(1,3,i), Mit(1,1), Mit(2,1), Mit(3,1)
   call write_info(TRIM(message_text))
   WRITE(message_text, fmt=trim(fmt_)) transf_mat(2,1,i), transf_mat(2,2,i), transf_mat(2,3,i),  Mit(1,2), Mit(2,2), Mit(3,2)
   call write_info(TRIM(message_text))
   WRITE(message_text, fmt=trim(fmt_)) transf_mat(3,1,i), transf_mat(3,2,i), transf_mat(3,3,i),  Mit(1,3), Mit(2,3), Mit(3,3)
   call write_info(TRIM(message_text))


   Mat_det  = matrix_determinant_33(Mat)
   WRITE(message_text, '(10x,a,F5.2)') '> det = ', Mat_det
   call write_info(TRIM(message_text))

  end do

 return
end subroutine write_list_matrice
!-----------------------------------------------------------------

subroutine find_new_group()
 USE cryscalc_module,           ONLY : SPG, Mat, message_text, nb_symm_op, Mit, Space_group_symbol, debug_proc
 USE macros_module,             ONLY : remove_car, replace_car, replace_car2
 USE IO_module,                 ONLY : write_info
 use CFML_Crystallographic_Symmetry, ONLY : set_spacegroup,   get_hallsymb_from_gener, space_group_type

 USE SHELX_module,              ONLY : car_symop
 USE matrix_module,             ONLY : MATMUL_33, MV_product

 implicit none
  INTEGER                          :: op_num
  INTEGER                          :: i, i1, i2
  CHARACTER (LEN=40)               :: input_line
  CHARACTER (LEN=12), DIMENSION(3) :: op_STRING
  REAL, DIMENSION(3,3)             :: op_ROT
  REAL, DIMENSION(3)               :: op_TRANS
  REAL, DIMENSION(3,3)             :: new_op_ROT
  REAL, DIMENSION(3)               :: new_op_TRANS
  LOGICAL                          :: trigonal, hexagonal, tetragonal, cubic, permut

  TYPE (space_group_type)                      :: new_SPG


  if(debug_proc%level_2)  call write_debug_proc_level(2, "find_new_group")

 SPG%Bravais = ADJUSTL(SPG%Bravais)
 new_SPG%NumSPG = 0

 ! creation de la carte des operateurs de symmetrie

 select case (SPG%Bravais(1:1))
    case ('I') ! I
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y+1/2,Z+1/2'

    case ('R') ! Rom, Hex
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+2/3,Y+1/3,Z+1/3'
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/3,Y+2/3,Z+2/3'

    case ('F') ! F
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X,Y+1/2,Z+1/2'

    case ('A') ! A
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X,Y+1/2,Z+1/2'
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y,Z+1/2'
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y+1/2,Z'

    case ('B') ! B
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y,Z+1/2'

    case ('C') ! C
       nb_symm_op=nb_symm_op+1
       car_symop(nb_symm_op)='X+1/2,Y+1/2,Z'
  END select

  SPG%hexa = .false.
  permut   = .true.
  trigonal = .false.
  if(len_trim(SpG%CrystalSys) == 8 .and. SpG%CrystalSys(1:8) == "Trigonal")  then
   trigonal = .true.
   permut   = .false.
  end if
  hexagonal = .false.
  if(len_trim(SpG%CrystalSys) == 9 .and. SpG%CrystalSys(1:9) == "Hexagonal") then
   hexagonal = .true.
   permut    = .false.
  end if
  if(trigonal .or. hexagonal) SPG%hexa = .true.
  tetragonal = .false.
  if(len_trim(SpG%CrystalSys) == 10 .and. SpG%CrystalSys(1:10) == "Tetragonal") then
   tetragonal = .true.
   permut     = .false.
  end if
  cubic = .false.
  if(len_trim(SpG%CrystalSys) == 5 .and. SpG%CrystalSys(1:5) == "Cubic") then
   cubic  = .true.
   permut = .false.
  end if


  !do op_num =2, SPG%NumOps
  do op_num=1, SPG%multip

   !write(*,*) 'sym op. : ', TRIM(SPG%SymopSymb(op_num))
   !write(*,*) '  ', TRIM(SPG%SymopSymb(op_num))
   input_line = TRIM(SPG%SymopSymb(op_num))
   input_line = remove_car(input_line, ' ')
   input_line = replace_car(input_line, ',', ' ')
   i1=INDEX(input_line, ' ')
   i2=INDEX(TRIM(input_line), ' ', back=.TRUE.)
   op_STRING(1) = input_line(1    : i1-1)
   op_STRING(2) = input_line(i1+1 : i2-1)
   op_STRING(3) = input_line(i2+1 : )
   op_STRING = ADJUSTL(op_STRING)
   call get_op_ROT_TRANS(op_STRING, op_ROT, op_TRANS)

   call get_matrice_inverse_transposee   !!! very important !!
   new_op_ROT   = MATMUL(ABS(Mit), op_ROT)
   new_op_TRANS = MATMUL(Mit, op_TRANS)

   !new_op_ROT   = MATMUL(ABS(MAT), op_ROT)
   !new_op_TRANS = MATMUL(MAT, op_TRANS)
   !new_op_trans = op_trans

   call get_op_STRING(new_op_ROT, new_op_TRANS, op_string)
   op_string = ADJUSTL(op_string)


   if(permut) then
    op_string(1) = replace_car(op_string(1), "y", "x")
    op_string(1) = replace_car(op_string(1), "z", "x")
    op_string(2) = replace_car(op_string(2), "x", "y")
    op_string(2) = replace_car(op_string(2), "z", "y")
    op_string(3) = replace_car(op_string(3), "x", "z")
    op_string(3) = replace_car(op_string(3), "y", "z")
   end if

   if(trigonal) then
    do i=1, 3
     i1 = index(op_string(i), "-1/3")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "-1/3", "+1/3")
     IF(op_string(i)(1:1) == '+') op_string(i) = op_string(i)(2:)
     i1 = index(op_string(i), "-2/3")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "-2/3", "+2/3")
     IF(op_string(i)(1:1) == '+') op_string(i) = op_string(i)(2:)
    end do
   else
    do i=1, 3
     i1 = index(op_string(i), "-1/2")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "-1/2", "+1/2")
     i1 = index(op_string(i), "-1/4")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "-1/4", "+3/4")
     i1 = index(op_string(i), "-3/4")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "-3/4", "+1/4")


     IF(op_string(i)(1:1) == '+') op_string(i) = op_string(i)(2:)

     i1 = index(op_string(i), "x+x")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "x+x", "x")
     i1 = index(op_string(i), "-x-x")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "-x-x", "-x")
     i1 = index(op_string(i), "y+y")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "y+y", "y")
     i1 = index(op_string(i), "-y-y")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "-y-y", "-y")
     i1 = index(op_string(i), "z+z")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "z+z", "z")
     i1 = index(op_string(i), "-z-z")
     if(i1 /=0) op_string(i) = replace_car2(op_string(i), "-z-z", "-z")

    end do

   endif


   WRITE(car_symop(op_num), '(5a)')  op_string(1),  ',',  op_string(2), ',', op_string(3)
   car_symop(op_num) = remove_car(car_symop(op_num), ' ')
  end do


  ! deduction du nouveau symbol du groupe
 !WRITE(SPG_num_string, '(i3)')SPG%NumSpg

 nb_symm_op = SPG%multip
 call write_info('  > New symmetry operators list: ')
 do op_num=1, nb_symm_op
  write(message_text,'(5x,i3,2x,a)') op_num, car_symop(op_num)
  call write_info(TRIM(message_text))
 end do
 !call set_spacegroup(TRIM(SPG_num_string), SPG, car_symop, nb_symm_op, 'GEN')
 call set_spacegroup(" ", new_SPG, car_symop, nb_symm_op, 'GEN')


  ! call get_hallsymb_from_gener(new_SPG) ! pas necessaire car routine deja appelee par "set_spacegroup" avec option "GEN"

  if(new_SPG%NumSPG == 0) then
   WRITE(message_text, '(a)') '  > New Space group: can not be deduced from symmetry operators list'
  else
   WRITE(message_text, '(2a)') '  > New Space group: ', TRIM(new_SPG%SPG_Symb)
  end if
  call write_info('')
  call write_info(TRIM(message_text))


  IF(new_SPG%NumSPG /=0)  then
   SPG = new_SPG
   space_group_symbol = new_SPG%SPG_Symb
   call set_spacegroup(new_SPG%SPG_Symb, SPG)
   call write_symm_op_reduced_set()
  end if


 RETURN
end subroutine find_new_group


!----------------------------------------------------------------------

subroutine transf_FACES()
 USE cryscalc_module, ONLY      : FACES_file_name, M => Mat, unit_cell, crystal, message_text, debug_proc, tmp_unit, &
                                  keyword_CELL
 USE IO_module,       ONLY      : write_info
 USE matrix_module,   ONLY      : MV_product
 implicit none
  !real                         :: new_h, new_k, new_l
  REAL, DIMENSION(3)           :: new_H, Hr
  INTEGER                      :: i, i1
  character (len=256)          :: new_faces_file_name

 if(debug_proc%level_2)  call write_debug_proc_level(2, "transf_faces")

 i1 = index(FACES_file_name, '.', back=.true.)
 if(i1 /=0) then
  new_faces_file_name = FACES_file_name(1:i1-1)//'_new'//FACES_file_name(i1:)
 else
  new_faces_file_name = 'absorb_new.ins'
 end if

 open (unit = tmp_unit, file = trim(new_faces_file_name), action='write')
 if (keyword_CELL) then
  write(tmp_unit, '(a)')         '# The unit cell upon which faces are based:'
  write(tmp_unit, '(a, 6F10.4)') '# CELL ', unit_cell%new_param(1:6)
 end if

 do i=1, crystal%faces_nb

  Hr = real(crystal%face_index(:,i))
  new_H = MV_product(Hr,  M)

  WRITE(message_text,'(a,3F6.2,10x,a,3F6.2)') '  h1 k1 l1 : ', Hr(1), Hr(2), Hr(3), &
                                              '  h2 k2 l2 : ', new_H(1), new_H(2), new_H(3)
  call write_info(TRIM(message_text))

  !write(tmp_unit, '(a,3(1x,F6.2),1x,F9.4)') 'FACE ', new_H(1), new_H(2), new_H(3), crystal%face_dim(i)
  write(tmp_unit, '(a,3(1x,I5),1x,F9.4)') 'FACE ', int(new_H(1)), int(new_H(2)), int(new_H(3)), crystal%face_dim(i)
 END do

 close(unit=tmp_unit)

 call write_info('')
 call write_info(trim(new_faces_file_name)// ' has been created.')
 call write_info('')

 return
end subroutine transf_FACES

!----------------------------------------------------------------------

subroutine Create_HKLF5_file
 use CRYSCALC_module, only : keyword_FILE, keyword_MAT, Mat, HKL_unit, message_text, debug_proc
 use HKL_module
 USE matrix_module,   ONLY : MV_product
 use IO_module,       only : write_info

 implicit none
  integer                 :: i, j, i1, n, m
  integer                 :: n1, n2
  character (len=256)     :: HKLF5_file
  REAL, DIMENSION(3)      :: new_H, Hr
  logical                 :: new_indices_entiers
  REAL, parameter         :: EPS=0.15

 if(debug_proc%level_2)  call write_debug_proc_level(2, "Create_HKLF5")


 if(.not. keyword_MAT) then
  call write_info('')
  call write_info('  >> Transformation matrix not known !')
  call write_info('')
  return
 end if

 if(.not. keyword_FILE) then
  call write_info('')
  call write_info('  >> hkl file not known !')
  call write_info('')
  return
 end if


 i1 = index(HKL_file%name, '.')
 if(i1/=0) then
  write(HKLF5_file, '(2a)') HKL_file%name(1:i1-1), '_hklf5.hkl'
  open(unit=HKL_unit, file=trim(HKLF5_file))

 else
  call write_info('')
  call write_info('  >> Erroneous hkl file name !')
  call write_info('')
  return
 end if

 n1 = 0
 n2 = 0
 do i=1, n_ref
  Hr(1) = real(h(i))
  Hr(2) = real(k(i))
  Hr(3) = real(l(i))
  new_H = MV_product(Hr,  Mat)

  new_indices_entiers = .false.
  n=0
  do j=1, 3
   if (ABS(new_H(j) - floor(new_H(j)+.5)) < overlapping_criteria ) n=n+1
  end do
  if(n==3) new_indices_entiers = .true.

  if(new_indices_entiers) then
   n2 = n2 + 1
   m = -2
   write(unit=HKL_unit, fmt='(3I4,2F8.2,I4)') int(new_H(1)), int(new_H(2)), int(new_H(3)), F2(i), sig_F2(i), m
   m = 1
   write(unit=HKL_unit, fmt='(3I4,2F8.2,I4)') h(i), k(i), l(i), F2(i), sig_F2(i), m
  else
   n1 = n1 + 1
   m = 1
   write(unit=HKL_unit, fmt='(3I4,2F8.2,I4)') h(i), k(i), l(i), F2(i), sig_F2(i), m
  end if
 END do

 close(unit=HKL_unit)

 call write_info('')
 write(message_text, '(3a)') '  >> ', trim(HKLF5_file) , ' has been created (hklf5 format):'
 call write_info(trim(message_text))
 write(message_text, '(a,F6.2)') '      . reflections overlapping criteria= ', overlapping_criteria
 call write_info(trim(message_text))
 write(message_text, '(a,i6)')   '      . number of overlapped reflections= ', n2
 call write_info(trim(message_text))
 write(message_text, '(a,i6)')   '      . number of single reflections=     ', n1
 call write_info(trim(message_text))

 call write_info('')

 return
end subroutine Create_HKLF5_file

 !----------------------------------------------------------------------------------------------------------------------

 subroutine diffractometer_geometry(input_string, motor)
  use cryscalc_module, only : pi, message_text, input_line
  use IO_module,       only : read_input_line, write_info
  implicit none
   character (len=*), intent(in) :: input_string
   real, dimension(3)            :: motor
   real                          :: Phi_K, Phi_E, Om_K, Om_E, Kappa, Chi
   real                          :: alpha, sin_chi, sin_delta, cos_delta, delta, sin_kappa
   real, parameter               :: eps = 0.001
   LOGICAL                       :: input_values

   alpha = 50. * pi / 180.      ! en rad.

   input_values = .true.
   if(abs(motor(1)) < eps .and. abs(motor(2)) < eps .and. abs(motor(3)) < eps) input_values = .false.

   if(input_string(1:5) == 'KAPPA') then
    call write_info('')
    call write_info(' Diffractometer geometry : KAPPA')
    call write_info('')

    if(.not. input_values) then
     call write_info('   > Enter Phi, Omega and Kappa (in deg.) : ')
     call read_input_line(input_line)
     read(input_line, *) Phi_K, Om_K, Kappa
    else
     Phi_K = motor(1)
     Om_K  = motor(2)
     Kappa = motor(3)
    end if


    kappa = kappa * pi / 180.    ! en rad
    ! calcul de Chi_euler: sin(Chi_euler/2) = sin(alfa).sin(kappa/2)
    sin_chi = SIN(alpha) * SIN(kappa/2)
    chi     = (2*ASIN(sin_chi) * 180. / pi) * pi/180. ! en rad.

    ! calcul de delta:     cos(delta) = cos(kappa/2) / cos(chi_euler/2)
    cos_delta = COS(kappa/2) / COS(chi/2)
    delta = ACOS(cos_delta) * 180. / pi     ! en deg.

    ! calcul de omega_euler et phi_euler: om_euler  = om_kappa + delta
    !                                     phi_euler = phi_kappa + delta
    Om_E  = Om_K  + delta
    Phi_E = Phi_K + delta

    call write_info('')
    call write_info('     Kappa geometry:      Phi           Om        Kappa        Delta')
    call write_info('     ----------------- ')
    WRITE(message_text,'(16x,4(x,F12.3))' ) Phi_K, Om_K, Kappa*180./pi, delta
    call write_info(TRIM(message_text))

    call write_info('')
    call write_info(' =>> Eulerian geometry:   Phi           Om          Chi')
    call write_info('     ----------------- ')
    WRITE(message_text,'(16x,3(x,F12.3))' ) Phi_E, Om_E, Chi*180./pi
    call write_info(TRIM(message_text))
    call write_info('')


   elseif(input_string(1:5) == 'EULER') then
    call write_info('')
    call write_info(' Diffractometer geometry : EULER')
    call write_info('')
    if(.not. input_values) then
     call write_info('   > Enter Phi, Omega and Chi (in deg.) : ')
     call read_input_line(input_line)
     read(input_line, *) Phi_E, Om_E, Chi
    else
     Phi_E = motor(1)
     Om_E  = motor(2)
     Chi   = motor(3)
    end if

    if(Chi > 180.) Chi = 360.-Chi        ! new janvier 2017
    Chi = Chi * pi / 180.                ! en rad.

    ! calcul de delta: sin(delta) = cotg(alfa) * tg(chi_euler/2)
    sin_delta = TAN(Chi/2) / TAN(alpha)
    delta   = ASIN(sin_delta) *180./pi  ! en deg.

    ! calcul de omega_kappa et phi_kappa: om_kappa  = om_euler  - delta
    !                                     phi_kappa = phi_euler - delta
    Om_K  = Om_E  - delta
    Phi_K = Phi_E - delta

    ! calcul de kappa: sin(kappa/2) =  sin(chi_euler/2) / sin(alfa)
    sin_kappa = SIN(Chi/2) / SIN(alpha)
    Kappa = 2* ASIN(sin_kappa) * 180./pi   ! en deg.

    call write_info('')
    call write_info('     Eulerian geometry:   Phi           Om          Chi        Delta')
    call write_info('     ----------------- ')
    WRITE(message_text,'(16x,4(x,F12.3))' ) Phi_E, Om_E, Chi*180./pi, delta
    call write_info(TRIM(message_text))
     call write_info('')

    call write_info(' =>> Kappa geometry:      Phi           Om        Kappa')
    call write_info('     ----------------- ')
    WRITE(message_text,'(16x,3(x,F12.3))' ) Phi_K, Om_K, Kappa
    call write_info(TRIM(message_text))
    call write_info('')

   end if

  return
 end subroutine diffractometer_geometry

 !!--------------------------------------------------------------------------------------!!
 subroutine reduce_cell
  ! provides conventional unit cell parameters and transformation matrix from
  ! input cell paramters
  ! based on Get_Conv_CELL.F90 (CRYSFML\program_examples\metrics)

  use CRYSCALC_module,       only : unit_cell, crystal_cell, reduce_BL, unit_cell, write_details, message_text
  use CFML_Crystal_Metrics,  only : set_crystal_Cell, Crystal_Cell_Type, Get_Primitive_Cell, Niggli_Cell, Get_twofold_axes, &
                                    Twofold_Axes_Type, Get_Conventional_Cell
  use CFML_Math_3D,          only : determ_A, determ_V, invert_A
  use CFML_Math_General,     only : sort
  USE CFML_GlobalDeps,       ONLY : cp
  USE IO_module,             ONLY : write_info


  implicit none
   real(kind=cp), dimension(3,3)    :: transfm, trans, prod, finalm
   integer,       dimension(3,3)    :: tr
   real(kind=cp), dimension(12)     :: del
   type(Crystal_Cell_Type)          :: cell_P, celln, cell
   Type(Twofold_Axes_Type)          :: twofold, twf
   type(Twofold_Axes_Type), dimension(12) :: otwf   !Individual two-fold axes for searching
                                                     !monoclinic cells
   real(kind=cp)                    :: rmi, rma, tol, told, det
   integer                          :: i, j, n, p, nold, ntwot
   character(len=11)                :: metr
   integer,       dimension(12)     :: ind
   character(len=80)                :: message
   logical                          :: ok

   call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)

   call Get_Primitive_Cell(reduce_BL, crystal_cell, cell_P, transfm)
   call Niggli_Cell(cell_P,celln=celln,trans=trans)

   tol  = 3.0  ! angular tolerance
   told = 0.2  ! distance tolerance
   call Get_twofold_axes(celln, tol, twofold)

   write(message_text,fmt="(a,6f10.4,a)")  "             Input Cell (Aic): ",unit_cell%param,"  Centring: "//reduce_BL
   call write_info(trim(message_text))
   write(message_text,fmt="(a,6f10.4)")    "   Primitive input Cell (Api): ",Cell_p%cell,Cell_p%ang
   call write_info(trim(message_text))
   write(message_text,fmt="(a,6f10.4)")    "            Niggli Cell  (AN): ",celln%cell, celln%ang
   call write_info(trim(message_text))
   call write_info('')
   write(message_text,fmt="(a)")    "             (Aic) : formal column matrix containing the input cell vectors"
   call write_info(trim(message_text))
   write(message_text,fmt="(a)")    "             (Api) : formal column matrix containing the primitive cell vectors"
   call write_info(trim(message_text))
   write(message_text,fmt="(a)")    "              (AN) : formal column matrix containing the Niggli cell vectors"
   call write_info(trim(message_text))
   write(message_text,fmt="(a)")    "             (Acc) : formal column matrix containing the conventional cell vectors"
   call write_info(trim(message_text))
   call write_info('')
   write(message_text,fmt="(a)")   "                            (Api) = M (Aic)                          (AN) = N (Api)"
   call write_info(trim(message_text))
   do i=1,3
     write(message_text,fmt="(a,3f12.6,tr5,3f12.6)") "     TransF:    ",transfm(i,1),transfm(i,2), transfm(i,3), &
                                                                        trans(i,1), trans(i,2), trans(i,3)
     call write_info(trim(message_text))
   end do
   prod = matmul(trans, transfm)
   metr= "    Pseudo-"
   ntwot=0



   if( twofold%ntwo > 0) then
    call write_info('')
    write(message_text,fmt="(a)")        " => Two-fold axes (indices in the Niggli cell)"
    call write_info(trim(message_text))
    write(message_text,fmt="(a)")        " =>       Direct       Reciprocal    Dot    Cross      Length "
    call write_info(trim(message_text))
    rma=-100.0
    do i=1,twofold%ntwo
     write(message_text,fmt="(a,3i4,tr4,3i4,i6,f10.3,f12.5)")  "     ", &
           twofold%dtwofold(:,i),twofold%rtwofold(:,i),twofold%dot(i),twofold%cross(i),twofold%maxes(i)
     call write_info(trim(message_text))
     if( twofold%cross(i) > rma) rma= twofold%cross(i)
     del(i)=twofold%cross(i)
    end do
    call sort(del,twofold%ntwo,ind)

    ntwot = twofold%ntwo
    do i=1,twofold%ntwo
     j=ind(i)
     del(twofold%ntwo-i+1)=twofold%cross(j)  !reorder by decreasing errors
     otwf(twofold%ntwo-i+1)%ntwo=1
     otwf(twofold%ntwo-i+1)%tol=tol
     otwf(twofold%ntwo-i+1)%caxes(:,1)=twofold%caxes(:,j)
     otwf(twofold%ntwo-i+1)%dtwofold(:,1)=twofold%dtwofold(:,j)
     otwf(twofold%ntwo-i+1)%rtwofold(:,1)=twofold%rtwofold(:,j)
     otwf(twofold%ntwo-i+1)%dot(1)=twofold%dot(j)
     otwf(twofold%ntwo-i+1)%cross(1)=twofold%cross(j)
     otwf(twofold%ntwo-i+1)%maxes(1)=twofold%maxes(j)
     otwf(twofold%ntwo-i+1)%a(:)=twofold%a(:)
     otwf(twofold%ntwo-i+1)%b(:)=twofold%b(:)
     otwf(twofold%ntwo-i+1)%c(:)=twofold%c(:)
    end do

    !!!------ FIRST: test with all twofold axes found
    call Get_Conventional_Cell(twofold,Cell,tr,message,ok,told)
    det=determ_A(tr)
    !Here Tr is the matrix transforming the Niggli cell to the conventional cell
    if(ok) then
     if(rma < 0.1)  metr="Metrically "
     if(abs(det) > 0) then
      call Write_New_Cell(metr, message, tr, prod, finalm, cell)
      call Write_New_Monoc_Cell(finalm, message, cell)
     end if
     !call Write_Crystal_Cell(Cell)

     if(.not. write_details) return

    !!!------ SECOND: select smaller blocks of twofold axes
    if( rma >= 0.1) then
     p=twofold%ntwo+1
     del(1)=del(1)-0.06
     nold=0
     do j=1,p-1       !Loop decreasing the tolerance for selecting better 2-fold axes
      rmi=del(j)+0.05
      if(rmi < 0.1) then
        rmi=0.1
        metr="Metrically "
      else
        metr="    Pseudo-"
      end if

      n=0
      do i=1,twofold%ntwo
       if(twofold%cross(i) > rmi) cycle
       n=n+1
       twf%caxes(:,n)    = twofold%caxes(:,i)
       twf%dtwofold(:,n) = twofold%dtwofold(:,i)
       twf%rtwofold(:,n) = twofold%rtwofold(:,i)
       twf%dot(n)        = twofold%dot(i)
       twf%cross(n)      = twofold%cross(i)
       twf%maxes(n)      = twofold%maxes(i)
      end do

      if( n == nold ) then
       cycle
      else
       nold=n
      end if
      if( n >= 1 ) then
       twf%ntwo=n
       twf%tol=twofold%tol
       twf%a=twofold%a
       twf%b=twofold%b
       twf%c=twofold%c
       call Get_Conventional_Cell(twf,Cell,tr,message,ok)
       det=determ_A(tr)
       if(ok) then
        call Write_New_Cell(metr, message, tr, prod, finalm, cell)
        call Write_New_Monoc_Cell(finalm, message, cell)
        call write_info("")
        write(message_text, fmt="(a,i3,a,f10.4)")  " => Two-fold axes:  ",n, "  Angular Discrepancy (deg)", rmi
        call write_info(trim(message_text))
        write(message_text,fmt="(a)")              " =>       Direct       Reciprocal    Dot    Cross      Length "
        call write_info(trim(message_text))
        do i=1,n
         write(message_text,fmt="(a,3i4,tr4,3i4,i6,f10.3,f12.5)")  "     ",&
                     twf%dtwofold(:,i),twf%rtwofold(:,i), twf%dot(i),twf%cross(i),twf%maxes(i)
         call write_info(trim(message_text))
        end do
        !else
        !  write(unit=io, fmt="(  a)")   " => No lattice results: "//trim(message)
       end if
      end if
      if(n <= 1) exit
      twofold=twf   !copy back on Twofold
     end do
    end if !rma >0.1
   else
    write(message_text,fmt="(a)") " => An unexpected error has occurred: "//trim(message)
    call write_info(trim(message_text))
    write(message_text,fmt="(a)") " => Change the angular/distance tolerance to obtain proper two-fold axes."
    call write_info(trim(message_text))
   end if !ok

  else
   !The Niggli cell is accepted as triclinic cell
   call write_info('')
   write(message_text,fmt="(a,3f10.5,3f9.3)") "  Cell (Triclinic, Niggli Cell): ",Celln%cell,Celln%ang
   call write_info(trim(message_text))
   call write_info('')
   write(message_text,fmt="(a,3f12.6,a)")     "                                   /",prod(1,1),prod(1,2),prod(1,3),"\"
   call write_info(trim(message_text))
   write(message_text,fmt="(a,3f12.6,a)")     "     Final Tranformation Matrix:  | ",prod(2,1),prod(2,2),prod(2,3)," |"
   call write_info(trim(message_text))
   write(message_text,fmt="(a,3f12.6,a)")     "                                   \",prod(3,1),prod(3,2),prod(3,3),"/"
   call write_info(trim(message_text))
  end if
  if(ntwot == 0) return


    !!!------ THIRD: Determine all monoclinic cells considering one-by-one all the two-fold axes
    write(message_text,fmt="(a,i3,a)")       " => There are ",ntwot," two-fold axes"
    call write_info('')
    call write_info(trim(message_text))
    call write_info('')
    write(message_text, fmt="(a,120a1)") " ",("*",i=1,120)
    call write_info(trim(message_text))
    write(message_text,fmt="(a)")   &
       " => MONOCLINIC CELLS: In this section each two-fold axis is considered as a unique b-axis of M-cells"
    call write_info(trim(message_text))
    write(message_text,fmt="(a)")   &
       "    The metrics of the monoclinic cells may be of higher symmetry, e.g. beta=90.0"
    call write_info(trim(message_text))
     write(message_text, fmt="(a,120a1)") " ",("*",i=1,120)
    call write_info(trim(message_text))

    do j=1,ntwot
     call Get_Conventional_Cell(otwf(j),Cell,tr,message,ok)
     det=determ_A(tr)
     rma=otwf(j)%cross(1)
     if(rma < 0.1) then
      metr="Metrically "
     else
      metr="    Pseudo-"
     end if
     if(abs(det) > 0) then
      call Write_New_Cell(metr, message, tr, prod, finalm, cell)
      call Write_New_Monoc_Cell(finalm, message, cell)
     end if
      call write_info('')
     write(message_text, fmt="(a,f10.4)")  " => Two-fold axis along monoclinic b-axis,  Angular Discrepancy (deg)", rma
     call write_info(trim(message_text))
     write(message_text,fmt="(a)")           " =>       Direct       Reciprocal    Dot    Cross      Length "
     call write_info(trim(message_text))
     write(message_text,fmt="(a,3i4,tr4,3i4,i6,f10.3,f12.5)")  "     ",&
         otwf(j)%dtwofold(:,1),otwf(j)%rtwofold(:,1), otwf(j)%dot(1),otwf(j)%cross(1),otwf(j)%maxes(1)
     call write_info(trim(message_text))
    end do


  return
 end subroutine reduce_cell

!!-------------------------------------------------------------------------------------------------------------------
  Subroutine Write_New_Cell(metr, message, tr, prod, finalm, cell)
   use cryscalc_module,       only : write_details, message_text
   USE CFML_GlobalDeps,       ONLY : cp
   use CFML_Math_3D,          only : determ_A, invert_A
   use CFML_Crystal_Metrics,  only : Crystal_Cell_Type
   USE IO_module,             ONLY : write_info
   implicit none
   character(len=11),             intent(in)    :: metr
   character(len=80),             intent(in)    :: message
   integer,       dimension(3,3), intent(in)    :: tr
   real(kind=cp), dimension(3,3), intent(in)    :: prod
   real(kind=cp), dimension(3,3), intent(out)   :: finalm
   type(Crystal_Cell_Type),       intent(in)    :: cell
   real(kind=cp), dimension(3,3)                :: invm

   integer        :: i
   real(kind=cp)  :: det

   call write_info("")
   write(message_text, fmt="(a, 120a1)") " ",("-",i=1,120)
   call write_info(trim(message_text))
   write(message_text,fmt="(  a)")   &
   " => The new Cell is "//metr//trim(message)//" and the transformation matrix from then Niggli cell is:"
   call write_info(trim(message_text))
   finalm=matmul(real(tr,kind=cp),prod)
   det=determ_A(tr)
   write(message_text, fmt="(a, 120a1)") " ",("-",i=1,120)
   call write_info(trim(message_text))
   call write_info("")
   write(message_text,fmt="(a,i3,2i4,a)")     "                         /",tr(1,1), tr(1,2),tr(1,3)," \"
   call write_info(trim(message_text))
   write(message_text,fmt="(a,i3,2i4,a,i4)")  "  (Acc) = Tr (AN);  Tr: | ",tr(2,1), tr(2,2),tr(2,3),"  |   Determinant: ",nint(det)
   call write_info(trim(message_text))
   write(message_text,fmt="(a,i3,2i4,a)")     "                         \",tr(3,1), tr(3,2),tr(3,3), " /"
   call write_info(trim(message_text))
   call write_info("")
   write(message_text,fmt="(a,3f10.5,3f9.3)") "  Conventional Cell: ",Cell%cell,Cell%ang
   call write_info(trim(message_text))

   call write_info("")
   write(message_text,fmt="(a,3f12.6,a)")     "                                   /",finalm(1,1),finalm(1,2),finalm(1,3), "\"
   call write_info(trim(message_text))
   write(message_text,fmt="(a,3f12.6,2a)")    "     Final Tranformation Matrix:  | ",finalm(2,1),finalm(2,2),finalm(2,3)," |", &
                                              "   (Acc) = Ftr (Aic)"
   call write_info(trim(message_text))
   write(message_text,fmt="(a,3f12.6,a)")     "                                   \",finalm(3,1),finalm(3,2),finalm(3,3),"/"
   call write_info(trim(message_text))
   call write_info("")
   invm=Invert_A(finalm)
   det=determ_A(finalm)
   write(message_text,fmt="(a,f12.6)")        "     Determinant: ",det
   call write_info(trim(message_text))
   call write_info('')
   write(message_text,fmt="(a,3f12.6,a)")     "                                   /",invm(1,1),invm(1,2),invm(1,3),  "\"
   call write_info(trim(message_text))
   write(message_text,fmt="(a,3f12.6,2a)")    "   Inverse Tranformation Matrix:  | ",invm(2,1),invm(2,2),invm(2,3),  " |", &
                                              "   (Aic) = (Ftr)^(-1) (Acc)"
   call write_info(trim(message_text))
   write(message_text,fmt="(a,3f12.6,a)")     "                                   \",invm(3,1),invm(3,2),invm(3,3),  "/"
   call write_info(trim(message_text))
   call write_info("")
   return
 End Subroutine Write_New_Cell

!!------------------------------------------------------------------------------------------------------------------
 Subroutine Write_New_Monoc_Cell(finalm, message, cell)
  USE CFML_GlobalDeps,       ONLY : cp
  use CFML_Crystal_Metrics,  only : Change_Setting_Cell, Crystal_Cell_Type
  use CFML_Math_3D,          only : determ_A, invert_A
  use cryscalc_module,       only : message_text
  USE IO_module,             ONLY : write_info
  implicit none
   real(kind=cp), dimension(3,3),             intent(inout)    :: finalm
   character(len=80),             intent(in)  :: message
   type(Crystal_Cell_Type),       intent(in)  :: cell
   type(Crystal_Cell_Type)          :: cellt
   real(kind=cp)                    :: det
   real(kind=cp), dimension(3,3)    :: invm
   logical                          :: cell_trans
   real(kind=cp), dimension(3,3)    :: mat


   cell_trans=.false.
   if(trim(message) == "Monoclinic, A-centred cell") then
    mat=reshape( (/0.0,0.0,1.0, 0.0,-1.0,0.0, 1.0,0.0,0.0/),(/3,3/))
    call Change_Setting_Cell(Cell,Mat,Cellt)
    cell_trans=.true.
   else if(trim(message) == "Monoclinic, I-centred cell") then
    mat=reshape( (/1.0,0.0,1.0, 0.0,1.0,0.0, 0.0,0.0,1.0/),(/3,3/))
    call Change_Setting_Cell(Cell,Mat,Cellt)
    if(Cellt%ang(2) < 80.0) then
     mat=reshape( (/1.0,0.0,-1.0, 0.0, 1.0,0.0, 0.0,0.0, 1.0/),(/3,3/))
     call Change_Setting_Cell(Cell,Mat,Cellt)
    end if
    cell_trans=.true.
   end if
   if(cell_trans) then
    finalm=matmul(Mat,finalm)
    write(message_text,fmt="(a,3f10.5,3f9.3)") "  Equivalent C-centred Cell: ",Cellt%cell,Cellt%ang
    call write_info(trim(message_text))
    call write_info('')
    write(message_text,fmt="(a,3f12.6,a)")     "                                   /",finalm(1,1),finalm(1,2),finalm(1,3),"\"
    call write_info(trim(message_text))
    write(message_text,fmt="(a,3f12.6,2a)")    "     Final Tranformation Matrix:  | ",finalm(2,1),finalm(2,2),finalm(2,3)," |", &
                                               "     (Acc) = Ftr (Aic)"
    call write_info(trim(message_text))
    write(message_text,fmt="(a,3f12.6,a)")     "                                   \",finalm(3,1),finalm(3,2),finalm(3,3),"/"
    call write_info(trim(message_text))
    call write_info('')

    det=determ_A(finalm)
    invm=Invert_A(finalm)
    write(message_text,fmt="(a,f12.6)")      "     Determinant: ",det
    call write_info(trim(message_text))
    write(message_text,fmt="(a,3f12.6,a)")   "                                   /",invm(1,1),invm(1,2),invm(1,3),  "\"
    call write_info(trim(message_text))
    write(message_text,fmt="(a,3f12.6,2a)")  "   Inverse Tranformation Matrix:  | ",invm(2,1),invm(2,2),invm(2,3),  " |", &
                                             "    (Aic) = (Ftr)^(-1) (Acc)"
    call write_info(trim(message_text))
    write(message_text,fmt="(a,3f12.6,a)")   "                                   \",invm(3,1),invm(3,2),invm(3,3),  "/"
    call write_info(trim(message_text))
   end if

   return
 End Subroutine Write_New_Monoc_Cell



 !!------------------------------------------------------------------------------------------------------!!

 subroutine Get_Transformation_matrix
  USE CRYSCALC_module, only  : keyword_cell, unit_cell, message_text
  USE MACROS_module,   only  : nombre_de_colonnes
  USE IO_module
  USE CFML_Crystal_Metrics,  ONLY  : Crystal_cell_type, set_crystal_Cell, Get_transfm_Matrix, Change_Setting_Cell
  USE CFML_GlobalDeps,       ONLY : cp
  implicit none
   character (len=256)             :: input_line
   real, dimension(6)              :: cell_param
   TYPE (crystal_cell_type)        :: cell_1, cell_2
   logical                         :: ok
   integer                         :: i
   real(kind=cp), dimension(3,3)   :: trm



  if(.not. keyword_CELL) then
   call write_info('')
   call write_info(' >> Input cell parameters for cell #1: ')
   call read_input_line(input_line)
   call Get_cell_param_from_line(input_line, "CELL", cell_param, ok)
   if(.not. ok) return
   call set_crystal_Cell(cell_param(1:3), cell_param(4:6), cell_1)
  else
   !cell_1%cell(1:3) = unit_cell%param(1:3)
   !cell_1%ang(1:3)  = unit_cell%param(4:6)
   call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), cell_1)
  end if


  call write_info('')
  call write_info(' >> Input cell parameters for cell #2: ')
  call read_input_line(input_line)
  call Get_cell_param_from_line(input_line, "CELL", cell_param, ok)
  if(.not. ok) return
  call set_crystal_Cell(cell_param(1:3), cell_param(4:6), cell_2)

  call write_info('')
  WRITE(message_text,'( a,6F12.5)') '  > Parameters for CELL #1: ', (cell_1%cell(i),i=1,3), (cell_1%ang(i),i=1,3)
  call write_info(TRIM(message_text))
  WRITE(message_text,'( a,6F12.5)') '  > Parameters for CELL #2: ', (cell_2%cell(i),i=1,3), (cell_2%ang(i),i=1,3)
  call write_info(TRIM(message_text))

  call Get_transfm_Matrix(cell_1, cell_2, trm, ok)

  if(ok)  then
   call write_info('')
   write(message_text, '(a)')         '   Transformation matrix from cell #1 to cell #2:  '
   call write_info(trim(message_text))
   write(message_text, '(25x,3F4.0)')  trm(1, 1), trm(1, 2) , trm(1, 3)
   call write_info(trim(message_text))
   write(message_text, '(25x,3F4.0)')  trm(2, 1), trm(2, 2) , trm(2, 3)
   call write_info(trim(message_text))
   write(message_text, '(25x,3F4.0)')  trm(3, 1), trm(3, 2) , trm(3, 3)
   call write_info(trim(message_text))

  else
   call write_info('')
   write(message_text, '(a)')         '   No transformation matrix (with det=1) from cell #1 to cell #2 has been founded.'
   call write_info(trim(message_text))
   write(message_text, '(a)')         '   Check input cell paramters !'
   call write_info(trim(message_text))
  end if

  return
 end subroutine Get_Transformation_matrix

 !!-----------------------------------------------------------------------------------------------------

 subroutine Get_cell_param_from_line(input_line, input_string, param, ok)
  USE MACROS_module,   only  : nombre_de_colonnes
  USE IO_module
  implicit none
   character (len=*),  intent(in)  :: input_line
   character (len=*),  intent(in)  :: input_string
   real, dimension(6), intent(out) :: param
   logical,            intent(out) :: ok
   integer                         :: nb_arg, ier


   call nombre_de_colonnes(input_line, nb_arg)
   if(nb_arg == 1) then
    read(input_line, *, iostat=ier) param(1)
    param(2) = param(1)
    param(3) = param(1)
    if(input_string(1:4) == 'CELL') then
     param(4) = 90.
     param(5) = 90.
     param(6) = 90.
    else
     param(4) = 0.
     param(5) = 0.
     param(6) = 0.
    end if
   elseif(nb_arg == 2) then
    read(input_line, *, iostat=ier) param(1), param(3)
    param(2) = param(1)
    if(input_string(1:4) == 'CELL') then
     param(4) = 90.
     param(5) = 90.
     param(6) = 90.
    else
     param(4) = 0.
     param(5) = 0.
     param(6) = 0.
    end if
   elseif(nb_arg == 3) then
    read(input_line, *, iostat=ier) param(1), param(2), param(3)
    if(input_string(1:4) == 'CELL') then
     param(4) = 90.
     param(5) = 90.
     param(6) = 90.
    else
     param(4) = 0.
     param(5) = 0.
     param(6) = 0.
    endif
   elseif(nb_arg == 6) then
    read(input_line, *, iostat=ier) param(1), param(2), param(3), &
                                    param(4), param(5), param(6)
   else
    if(input_string(1:4) == 'CELL') then
     call write_info('  !! Wrong input cell parameters !!')
    else
     call write_info('  !! Wrong input ESD cell parameters !!')
    end if
    call write_info('')
    ok = .false.
    return
   end if

   ok = .true.
  return
 end subroutine Get_cell_param_from_line

