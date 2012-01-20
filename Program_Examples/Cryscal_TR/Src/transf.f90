!     Last change:  TR   30 Jan 2007    4:38 pm

subroutine write_matrice()
 USE cryscal_module,      ONLY : M => Mat, Mat_det, message_text
 USE IO_module,           ONLY : write_info
 USE MATRIX_list_module,  ONLY : matrix_num, transf_mat_text
 use matrix_module,       only : matrix_determinant_33
 


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
 WRITE(message_text, '(a,F6.3)') '  > Matrix determinant: ', Mat_det 
 call write_info(TRIM(message_text))
 call write_info('')


end subroutine write_matrice
!--------------------------------------------------------
subroutine  transf_cell_parameters
 USE cryscal_module, ONLY   : M => mat, Mat_det, unit_cell, pi, message_text,   &
                              WRITE_triclinic_transf, WRITE_monoclinic_transf,  &
                              WRITE_twin_hexa, WRITE_twin_pseudo_hexa, update_parameters
 USE IO_module,      ONLY   : write_info
 implicit none
  real                               :: cosAlfa1, cosBeta1, cosGamma1
  real                               :: cosAlfa2, cosBeta2, cosGamma2
  real                               :: a2b2, a2c2, b2c2
  INTEGER                            :: i
  REAL, DIMENSION(6)                 :: tmp_cell


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
!  a2b2  =    M(1,1) * M(2,1) * ESD_cell(1) ** 2. +  M(1,2) * M(2,2) * ESD_cell(2) ** 2.  +   &
!             M(1,3) * M(2,3) * ESD_cell(3) ** 2.  +                                          &
!             ESD_cell(1) * ESD_cell(2) * cosGamma1 * (M(1,1) * M(2,2) + M(1,2) * M(2,1))  +  &
!             ESD_cell(1) * ESD_cell(3) * cosbeta1  * (M(1,1) * M(2,3) + M(1,3) * M(2,1))  +  &
!             ESD_cell(2) * ESD_cell(3) * cosAlfa1  * (M(1,3) * M(2,3) + M(1,3) * M(2,2))
!
!
!  a2c2  =    M(1,1) * M(3,1) * ESD_cell(1) ** 2. +  M(1,2) * M(3,2) * ESD_cell(2) ** 2.  +   &
!             M(1,3) * M(3,3) * ESD_cell(3) ** 2.  +                                          &
!             ESD_cell(1) * ESD_cell(2) * cosGamma1 * (M(1,1) * M(3,2) + M(1,2) * M(3,1))  +  &
!             ESD_cell(1) * ESD_cell(3) * cosbeta1  * (M(1,1) * M(3,3) + M(1,3) * M(3,1))  +  &
!             ESD_cell(2) * ESD_cell(3) * cosAlfa1  * (M(1,3) * M(3,2) + M(1,3) * M(3,2))
!
!
!  b2c2  =    M(2,1) * M(3,1) * ESD_cell(1) ** 2. +  M(2,2) * M(3,2) * ESD_cell(2) ** 2.  +   &
!             M(2,3) * M(3,3) * ESD_cell(3) ** 2.  +                                          &
!             ESD_cell(1) * ESD_cell(2) * cosGamma1 * (M(2,1) * M(3,2) + M(2,2) * M(3,1))  +  &
!             ESD_cell(1) * ESD_cell(3) * cosbeta1  * (M(2,1) * M(3,3) + M(2,3) * M(3,1))  +  &
!             ESD_cell(2) * ESD_cell(3) * cosAlfa1  * (M(2,3) * M(3,2) + M(2,3) * M(3,2))
!
!  cosGamma2 = a2b2 / (ESD_new_cell(1)*ESD_new_cell(2))
!  cosBeta2  = a2c2 / (ESD_new_cell(1)*ESD_new_cell(3))
!  cosAlfa2  = b2c2 / (ESD_new_cell(2)*ESD_new_cell(3))
!
!
!  ESD_new_cell(4) = ACOS(cosAlfa2)  * 180./pi
!  ESD_new_cell(5) = ACOS(cosBeta2)  * 180./pi
!  ESD_new_cell(6) = ACOS(cosGamma2) * 180./pi




 call write_info('')
 call write_info('  New cell parameters:')
 write(message_text,'(10x,6F9.4)') unit_cell%new_param(1:6)
  call write_info(TRIM(message_text))
 !write(message_text,'(10x,3F9.4)') unit_cell%new_param_ESD(1:3)
 ! call write_info(TRIM(message_text))

 call write_info('')
 IF(unit_cell%volume < 0.1) call volume_calculation('')

 WRITE(message_text,'(a,F12.2,a)') '  New cell volume : ', abs(Mat_det) * unit_cell%volume, ' A3'
 call write_info(TRIM(message_text))
 call write_info('')

 tmp_cell = unit_cell%param
 
 if (      .not. WRITE_triclinic_transf  .and.  .not. WRITE_monoclinic_transf &
     .and. .not. WRITE_twin_hexa         .and.  .not. WRITE_twin_pseudo_hexa) then
	 
    if(update_parameters)  unit_cell%param = unit_cell%new_param
 end if	  
 !unit_cell%param = tmp_cell
 !call volume_calculation('out')




end subroutine  transf_cell_parameters

!----------------------------------------------------------------------

subroutine transf_HKL()
 USE cryscal_module, ONLY      : M => Mat, nb_hkl, H, message_text
 USE IO_module,      ONLY      : write_info
 USE matrix_module,  ONLY      : MV_product
 implicit none
  !real                         :: new_h, new_k, new_l
  REAL, DIMENSION(3)           :: new_H
  INTEGER                      :: i



 do i=1, nb_hkl

  new_H = MV_product(real(H(:,i)),  M)

  WRITE(message_text,'(a,3F6.2,10x,a,3F6.2)') '  h1 k1 l1 : ', H(1,i), H(2,i), H(3,i) , &
                                               '  h2 k2 l2 : ', new_H(1), new_H(2), new_H(3)
  call write_info(TRIM(message_text))

 END do

call write_info('')
end subroutine transf_HKL

!----------------------------------------------------------------------

subroutine transf_coord()
 USE cryscal_module, ONLY      : Mit, nb_atom, atom_coord, new_atom_coord, atom_label, atom_type, message_text, update_parameters
 USE IO_module,      ONLY      : write_info
 use matrix_module,  ONLY      : MV_product
 implicit none
  !real                         :: new_x, new_y, new_z
  !REAL, DIMENSION(3)           :: new_coord
  INTEGER                      :: i



 !call inversion_matrice()
 call get_matrice_inverse_transposee()
 call write_atomic_matrix()
 IF(nb_atom /=0) then
  call write_info('  New atomic coordinates: ')
  call write_info('')
  do i=1, nb_atom
   new_atom_coord(:,i) = MV_product(atom_coord(:, i), Mit)
   WRITE(message_text,'(a,2a6,3(1x,F10.6))') '  ATOM: ', trim(atom_label(i)),trim(atom_type(i)),  new_atom_coord(:,i)
   call write_info(TRIM(message_text))
   if(update_parameters) atom_coord(:,i) = new_atom_coord(:,i)
  END do
  call write_info('')
 endif


end subroutine transf_coord
!----------------------------------------------------------------------
subroutine write_atomic_matrix()
 USE IO_module,      ONLY : write_info
 USE cryscal_module, ONLY : Mit, message_text

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
 USE cryscal_module, ONLY : translat, nb_atom, atom_coord, atom_label, atom_type, message_text, update_parameters
 USE IO_module,      ONLY : write_info
 implicit none

  REAL, DIMENSION(3)           :: new_coord
  INTEGER                      :: i


 call write_info('')
 WRITE(message_text, '(a,3(1x,F10.6))') '  Unit cell translation:', translat(1:3)
 call write_info(TRIM(message_text))
 call write_info(' ')
 call write_info('  ... new coordinates :')
 call write_info(' ')

 do i=1, nb_atom

  new_coord = atom_coord(:, i) + translat

  WRITE(message_text,'(a,2a6,3(1x,F10.6))') '  ATOM: ', trim(atom_label(i)),trim(atom_type(i)),  new_coord
  call write_info(TRIM(message_text))

  if(update_parameters) atom_coord(:, i) = new_coord
 END do

 call write_info('')

end subroutine transl_coord

!----------------------------------------------------------------------

subroutine inside_unit_cell()
 USE cryscal_module, ONLY : nb_atom, atom_coord, atom_label, atom_type, message_text
 USE IO_module,      ONLY : write_info
 use matrix_module,  ONLY : MV_product
 implicit none
  !real                   :: new_x, new_y, new_z
  REAL, DIMENSION(3)      :: new_coord
  INTEGER                 :: i, j


 call write_info('')
 call write_info(' > List of atoms with atomic coordinates inside the unit cell')
 call write_info(' -------------------------------------------------------------')
 call write_info('')
 do i=1, nb_atom

  do j=1,3
   do
    IF(atom_coord(j,i) < 0.) then
     atom_coord(j,i) = atom_coord(j,i) + 1
    else
     exit
    endif
   END do

   do
    IF(atom_coord(j,i) > 1.) then
     atom_coord(j,i) = atom_coord(j,i) - 1
    else
     exit
    endif
   END do

  END do

  WRITE(message_text,'(a,2a6,3(1x,F9.6))') '  ATOM: ', trim(atom_label(i)),trim(atom_type(i)),  atom_coord(:,i)
  call write_info(TRIM(message_text))

 END do

 call write_info('')

end subroutine inside_unit_cell



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_matrice_inverse_transposee()
 USE cryscal_module, ONLY : M => Mat, Mi => Mat_inv, message_text, Mit, Mat_det
 USE matrix_module
 USE IO_module,      ONLY : write_info
 implicit none
  real, dimension(3,3)              :: Mp   ! matrice des complements algebriques
  REAL, DIMENSION(3,3)              :: Mpt  ! matrice transposee

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


!------------------------------------------------
subroutine  transf_HKL_file()
 use cryscal_module
 USE HKL_module,    ONLY  : HKL_file
 USE IO_module,     ONLY  : write_info
 USE matrix_module, ONLY  : MV_product
  implicit none
   integer                           :: i, i_error
   INTEGER, DIMENSION(3)             :: h_
   integer                           :: code_
   real                              :: F2_, sig_F2_
   real                              :: cos1,  cos2, cos3, cos4, cos5, cos6
   REAL, DIMENSION(3)                :: new_H
   logical                           :: new_indices_entiers
   REAL, parameter                   :: EPS=0.001
   integer                           :: n



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
   

   call write_info('')
   call write_info('  ... Reading data in '// trim(HKL_file%NAME)// ' ...')
   call write_info('')


   do
    READ(HKL_unit,'(3I4,2F8.2,I4,6F8.5)', iostat=i_error) h_(1), h_(2), h_(3), F2_, sig_F2_, code_,  &
	                                                      cos1,  cos2, cos3, cos4, cos5, cos6
    if (i_error < 0) exit ! fin du fichier atteinte
    
    new_H = MV_product(REAL(h_), Mat)
    
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

end subroutine transf_HKL_file


!-----------------------------------------------------------------
subroutine   transf_triclinic
 USE IO_module,      ONLY : write_info
 USE cryscal_module, ONLY : keyword_CELL, Mat, Mat_det
 USE Matrix_module,  ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none
  integer                :: mat_numor



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


!-----------------------------------------------------------------
subroutine   transf_monoclinic
 USE IO_module,      ONLY : write_info
 USE cryscal_module, ONLY : keyword_CELL, Mat, Mat_det
 USE Matrix_module,  ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none
  integer          :: mat_numor


  call write_info(' ')
  call write_info(' --------- Monoclinic transformations -----------')


  !do mat_numor = 27, 28
  do mat_numor = 27, 31
   matrix_num = mat_numor
   call write_info('')
   !call write_info('  *** M '//trim(transf_mat_text(mat_numor))//' ***')
   Mat(:,:) = transf_mat(:,:, mat_numor)
   call write_matrice()
   Mat_det = matrix_determinant_33(Mat)
   if (keyword_CELL) call transf_cell_parameters
  end do 

  call write_info('')
  call write_info('   . examples:')
  call write_info('')
  !call write_info('      P 21/c (a<c) ==> P 21/n (a<c):   matrix M #29')
  !call write_info('      P 21/c (a<c) ==> P 21/a (a<c):   matrix M #28')
  call write_info('')
  !call write_info('      P 21/n (a<c) ==> P 21/a (a>c):   matrix M #27')
  call write_info('      P 21/n (a<c) ==> P 21/c (a<c):   matrix M #28')
  call write_info('')
  !call write_info('      P 21/a (a<c) ==> P 21/n (a<c):   matrix M #28')
  call write_info('      P 21/a (a<c) ==> P 21/c (a>c):   matrix M #30')
  call write_info('')




 return
end subroutine   transf_monoclinic

!-----------------------------------------------------------------
subroutine   twin_hexa
 USE IO_module,      ONLY : write_info
 USE cryscal_module, ONLY : keyword_CELL, Mat, Mat_det
 USE Matrix_module,  ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none
  integer                :: mat_numor



  call write_info(' ')
  call write_info(' --------- Hexagonal transformations -----------')


  do mat_numor = 22, 26
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


!-----------------------------------------------------------------
subroutine   twin_pseudo_hexa
 USE IO_module,      ONLY : write_info
 USE cryscal_module, ONLY : keyword_CELL, Mat, Mat_det
 USE Matrix_module,  ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none
  integer       :: mat_numor



  call write_info(' ')
  call write_info(' --------- Pseudo-hexagonal twinning in a monoclinic cell -----------')


  do mat_numor = 27, 31
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


!-----------------------------------------------------------------
subroutine   permutation_abc
 USE IO_module,      ONLY : write_info
 USE cryscal_module, ONLY : keyword_CELL, Mat, Mat_det, Mat_det
 use matrix_module,  only : matrix_determinant_33
 USE Matrix_list_module
 
 implicit none
  integer               :: mat_numor



  call write_info(' ')
  call write_info(' --------- Transformations -----------')


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
end subroutine   permutation_abc



!-----------------------------------------------------------------
subroutine   transf_rhomb_to_hex
 USE IO_module,              ONLY : write_info
 USE cryscal_module,         ONLY : keyword_CELL, Mat, Mat_det
 USE Matrix_module,          ONLY : matrix_determinant_33
 USE MATRIX_list_module

 implicit none


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
  call write_info(' ')
  
  Mat(:,:) = transf_mat(:,:, 8)
  Mat_det = matrix_determinant_33(Mat)
  IF(keyword_CELL) call transf_cell_parameters


 return
end subroutine   transf_rhomb_to_hex
!-----------------------------------------------------------------
subroutine   transf_hex_to_rhomb
 USE cryscal_module,         ONLY : keyword_CELL, Mat, Mat_det
 USE IO_module,              ONLY : write_info
 USE Matrix_module,          ONLY : matrix_determinant_33
 USE MATRIX_list_module
 
 !USE MATRIX_list_module
 implicit none
   
 
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
  call write_info(' ')
  
  Mat(:,:) = transf_mat(:,:, 9)
  Mat_det = matrix_determinant_33(Mat)
  IF(keyword_CELL) call transf_cell_parameters


 return
end subroutine   transf_hex_to_rhomb
!-----------------------------------------------------------------

subroutine write_list_matrice()
 USE IO_module,      ONLY : write_info
 USE cryscal_module, ONLY : message_text, Mat, Mat_det, Mit
 USE matrix_module,  ONLY : matrix_determinant_33
 USE MATRIX_list_module
 implicit none
  INTEGER :: i

  call write_info('')
  call write_info('  > List of transformation matrices implemented in CRYSCAL: ')
  call write_info('')

  do i=1, Max_mat_nb + user_mat_nb
   Mat(:,:) =  transf_mat(:,:,i)
   call get_matrice_inverse_transposee()

!   call write_info('')
!  ! WRITE(message_text, '(a)', ADVANCE='no') TRIM(transf_mat_text(i))
!   WRITE(message_text, '(2x,a)') TRIM(transf_mat_text(i))
!   call write_info(TRIM(message_text))
!   call write_info('')
!   call write_info('        Direct matrix:')
!
!   WRITE(message_text, '(2x,3F10.5)')     transf_mat(1,:,i)
!   call write_info(TRIM(message_text))
!   WRITE(message_text, '(2x,3F10.5)')     transf_mat(2,:,i)
!   call write_info(TRIM(message_text))
!   WRITE(message_text, '(2x,3F10.5)')     transf_mat(3,:,i)
!   call write_info(TRIM(message_text))
!   Mat_det  = matrix_determinant_33(Mat)
!
!   WRITE(message_text, '(a,F5.2)') '    > det = ', Mat_det
!   call write_info(TRIM(message_text))
!
!   call write_info('')
!   call write_info('        Inverse matrix:')
!   WRITE(message_text, '(2x,3F10.5)')     Mit(:,1)
!   call write_info(TRIM(message_text))
!   WRITE(message_text, '(2x,3F10.5)')     Mit(:,2)
!   call write_info(TRIM(message_text))
!   WRITE(message_text, '(2x,3F10.5)')     Mit(:,3)
!   call write_info(TRIM(message_text))


   call write_info('')
   WRITE(message_text, '(2x,a)') TRIM(transf_mat_text(i))
   call write_info(TRIM(message_text))
   call write_info('')
   call write_info('        Direct matrix:                       Inverse matrix:')

   WRITE(message_text, '(2x,3F10.5,5x,3F10.5)')     transf_mat(1,:,i),  Mit(:,1)
   call write_info(TRIM(message_text))
   WRITE(message_text, '(2x,3F10.5,5x,3F10.5)')     transf_mat(2,:,i),  Mit(:,2)
   call write_info(TRIM(message_text))
   WRITE(message_text, '(2x,3F10.5,5x,3F10.5)')     transf_mat(3,:,i),  Mit(:,3)
   call write_info(TRIM(message_text))

   Mat_det  = matrix_determinant_33(Mat)
   WRITE(message_text, '(10x,a,F5.2)') '> det = ', Mat_det
   call write_info(TRIM(message_text))




  end do

end subroutine write_list_matrice
!-----------------------------------------------------------------

subroutine find_new_group()
 USE cryscal_module,            ONLY : SPG, Mat, message_text, nb_symm_op, Mit
 USE macros_module,             ONLY : remove_car, replace_car
 USE IO_module,                 ONLY : write_info
 use CFML_Crystallographic_Symmetry, ONLY : set_spacegroup,   get_hallsymb_from_gener
 USE SHELX_module,              ONLY : car_symop
 USE matrix_module,             ONLY : MATMUL_33, MV_product

 implicit none
  INTEGER                          :: op_num
  INTEGER                          :: i1, i2
  CHARACTER (LEN=40)               :: input_line
  CHARACTER (LEN=12), DIMENSION(3) :: op_STRING
  REAL, DIMENSION(3,3)             :: op_ROT
  REAL, DIMENSION(3)               :: op_TRANS
  REAL, DIMENSION(3,3)             :: new_op_ROT
  REAL, DIMENSION(3)               :: new_op_TRANS
  CHARACTER (LEN=12)               :: SPG_num_string

 SPG%Bravais = ADJUSTL(SPG%Bravais)

 ! creation de la carte des operateurs de symmetrie
 
 !write(*,*) ' num op : ', SPG%numops, SPG%multip
 !write(*,*) ' centro : ', SPG%centred
 
 ! SPG%NumOps : nombre d'operateurs de symetrie reduit, y compris l'identite
 !              
 !nb_symm_op = SPG%NumOps  - 1
 !do op_num = 1, nb_symm_op
 ! car_symop(op_num) = SPG%SymopSymb(op_num+1)
 ! write(*,*) op_num, car_symop(op_num)
 !end do
 
 !do op_num= 1, SPG%Multip
 ! WRITE(*,*)   op_num, SPG%SymopSymb(op_num)
 !end do

 
! IF(SPG%Centred  /=1)  then
!  nb_symm_op = nb_symm_op + 1
!  car_symop(nb_symm_op)='-X,-Y,-Z'
! endif

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


  !do op_num =2, SPG%NumOps
   do op_num=1, SPG%multip

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

   call get_op_STRING(new_op_ROT, new_op_TRANS, op_string)
   op_string = ADJUSTL(op_string)

   op_string(1) = replace_car(op_string(1), "y", "x")
   op_string(1) = replace_car(op_string(1), "z", "x")
   op_string(2) = replace_car(op_string(2), "x", "y")
   op_string(2) = replace_car(op_string(2), "z", "y")
   op_string(3) = replace_car(op_string(3), "x", "z")
   op_string(3) = replace_car(op_string(3), "y", "z")

   IF(op_string(1)(1:1) == '+') op_string(1) = op_string(1)(2:)
   IF(op_string(2)(1:1) == '+') op_string(2) = op_string(2)(2:)
   IF(op_string(3)(1:1) == '+') op_string(3) = op_string(3)(2:)

   WRITE(car_symop(op_num), '(5a)')  op_string(1),  ',',  op_string(2), ',', op_string(3)
   car_symop(op_num) = remove_car(car_symop(op_num), ' ')
   
  end do


  ! deduction du nouveau symbol du groupe  
 WRITE(SPG_num_string, '(i3)')SPG%NumSpg
 !write(*,*) SPG%numspg, SPG_num_string
 !nb_symm_op = nb_symm_op + 1
 !car_symop(1)="x,y,z"
 
 nb_symm_op = SPG%multip
 !write(*,*) ' nb_sym_op :  ', nb_symm_op
 !do op_num=1, nb_symm_op
 ! write(*,*) op_num, car_symop(op_num)
 !end do
 call set_spacegroup(TRIM(SPG_num_string), SPG, car_symop, nb_symm_op, 'GEN') 
 
 !call set_spacegroup(TRIM(SPG_num_string), SPG)
  
  call get_hallsymb_from_gener(SPG)

  WRITE(message_text, '(2a)') '  > New Space group: ', TRIM(SPG%SPG_Symb)
  call write_info('')
  call write_info(TRIM(message_text))

  IF(SPG%NumSPG /=0)  call write_symm_op_reduced_set()


 RETURN
end subroutine find_new_group


