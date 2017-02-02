!     Last change:  TR   17 Jul 2007    4:50 pm


!------------------------------------------------------------------------------
!
subroutine get_label_atom_coord(input_string, i, n, atom_index,ok)
! determination des coord. atomiques a partir du label
! rq: le label peut etre de la forme AA_$n, ou:
!                                     AA: label de l'atome dans la liste des atomes
!                                     $n: numero de l'operateur de symetrie (dans la liste) a appliquer a l'atome
 use cryscalc_module, ONLY      : R => symm_op_rot,  T => symm_op_trans, nb_dist_calc, nb_diff_calc, nb_symm_op, &
                                  nb_atom, atom1_dist, atom2_dist, atom1_ang, atom2_ang, atom3_ang, atom4_ang,   &
                                  atom1_plane, atom2_plane, atom3_plane, atom1_diff, atom2_diff,                 &
                                  atom_plane, atom_plane_nb, &
                                  atom_label, atom_coord, new_coord, message_text, debug_proc
 USE IO_module
 USE macros_module, only        : u_case

 implicit none
 integer, intent(inout)       :: atom_index
 CHARACTER(LEN=*), INTENT(IN) :: input_string
 INTEGER, INTENT(IN)          :: i,n
 LOGICAL, INTENT(INOUT)       :: ok
 INTEGER                      :: i1, j,long, num_sym_op
 CHARACTER(LEN=8)             :: label
 logical                      :: input_dist, input_ang, input_plane, input_diff

 if(debug_proc%level_2)  call write_debug_proc_level(2, "get_label_atom_coord ("//trim(input_string)//")")

 input_dist  = .false.
 input_ang   = .false.
 input_plane = .false.
 input_diff  = .false.


 if(len_trim(input_string) == 4) then
  if(input_string(1:4) == 'dist') input_dist = .true.
  if(input_string(1:4) == 'diff') input_diff = .true.
 elseif(len_trim(input_string) == 3) then
  if(input_string(1:3) == 'ang') input_ang = .true.
 elseif(len_trim(input_string) == 5) then
  if(input_string(1:5) == 'plane') input_plane = .true.
 endif


 IF(input_dist) then
  IF(n==1) then
   label = atom1_dist(i)
  ELSEIF(n==2) then
   label = atom2_dist(i)
  endif
 endif

 IF(input_diff) then
  IF(n==1) then
   label = atom1_diff(i)
  ELSEIF(n==2) then
   label = atom2_diff(i)
  endif
 endif


 IF(input_ang) then
  IF(n==1) then
   label = atom1_ang(i)
  ELSEif(n==2) then
   label = atom2_ang(i)
  ELSEif(n==3) then
   label = atom3_ang(i)
  ELSEif(n==4) then
   label = atom4_ang(i)
  endif
 endif

 !IF(input_plane) then
 ! IF(n==1) then
 !  label = atom1_plane(i)
 ! ELSEif(n==2) then
 !  label = atom2_plane(i)
 ! ELSEif(n==3) then
 !  label = atom3_plane(i)
 ! endif
 !endif

 if(input_plane) then
  label = atom_plane(n, i)
 end if



 i1 = index(label, '$*')
 if(i1 /=0) then
  call write_info('')
  call write_info(' ... Wrong input label ...')
  return
 end if

 i1 = INDEX(label, '_$')
 if (i1 /=0) then
  long = LEN_TRIM(label)
  READ(label(i1+2:long),*) num_sym_op        ! numero de l'operateur de symetrie
  IF(num_sym_op >  nb_symm_op) then
   call write_info('')
   WRITE(message_text, '(2x,a)') ' ... Symmetry operator not available ...'
   call write_info(TRIM(message_text))
   call write_info('')
   return
  endif

  label = label(1:i1-1)     ! label de l'atome
  atom_index = 0
  do j=1, nb_atom
   if (u_case(label) == u_case(atom_label(j))) atom_index = j
  end do

  IF(atom_index == 0) then
   call write_info('')
   WRITE(message_text, '(2x,3a)') ' ... ', TRIM(label), ' atom not available ...'
   call write_info(TRIM(message_text))
   call write_info('')
   return
  endif

 ! apply symm. op
   new_coord(1) = R(1,1, num_sym_op) * atom_coord(1, atom_index)  + &
                  R(1,2, num_sym_op) * atom_coord(2, atom_index)  + &
                  R(1,3, num_sym_op) * atom_coord(3, atom_index)  + T(1, num_sym_op)
   new_coord(2) = R(2,1, num_sym_op) * atom_coord(1, atom_index)  + &
                  R(2,2, num_sym_op) * atom_coord(2, atom_index)  + &
                  R(2,3, num_sym_op) * atom_coord(3, atom_index)  + T(2, num_sym_op)
   new_coord(3) = R(3,1, num_sym_op) * atom_coord(1, atom_index)  + &
                  R(3,2, num_sym_op) * atom_coord(2, atom_index)  + &
                  R(3,3, num_sym_op) * atom_coord(3, atom_index)  + T(3, num_sym_op)
 else
  atom_index = 0

  do j= 1, nb_atom
    if (u_case(label) == u_case(atom_label(j)) ) atom_index = j
  end do
  IF(atom_index == 0) then
   call write_info('')
   WRITE(message_text, '(2x,3a)') ' ... ', TRIM(label), ' atom not available ...'
   call write_info(TRIM(message_text))
   call write_info('')
   return
  endif
  new_coord(1:3) = atom_coord(1:3,atom_index)

 end if

  ok= .true.

 return
END subroutine get_label_atom_coord


!-------------------------------------------------------------------------------
subroutine  calcul_angles
 use cryscalc_module,         ONLY : nb_ang_calc, atom1_ang, atom2_ang, atom3_ang, atom4_ang, SP_value, new_coord, message_text, &
                                     debug_proc
 USE CFML_Math_General,       ONLY : acosd
 !USE CFML_Constants,          ONLY : sp
 USE CFML_GlobalDeps,                 ONLY : sp
 USE IO_module

 implicit none
  integer                 :: i
  REAL, DIMENSION(3)      :: atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord
  REAL(KIND=sp)           :: angle
  CHARACTER(LEN=16)       :: label_1, label_2, label_3, label_4
  LOGICAL                 :: ok
  integer                 :: atom_index

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_angles")

  call write_info('')
  call write_info('    . ANGLES CALCULATIONS')
  call write_info('    ---------------------')
  call write_info('')

  do i = 1, nb_ang_calc
   label_1 = atom1_ang(i)
   label_2 = atom2_ang(i)
   label_3 = atom3_ang(i)
   label_4 = atom4_ang(i)

   ok = .false.
   call get_label_atom_coord('ang', i,1, atom_index, ok)
   IF(.NOT. ok) cycle
   atom_1_coord(1:3) = new_coord(1:3)

   ok = .false.
   call get_label_atom_coord('ang', i,2, atom_index, ok)
   IF(.NOT. ok) cycle
   atom_2_coord(1:3) = new_coord(1:3)

   ok = .false.
   call get_label_atom_coord('ang', i,3, atom_index, ok)
   IF(.NOT. ok) cycle
   atom_3_coord(1:3) = new_coord(1:3)

!   ok = .false.
!   call get_label_atom_coord('ang', i,4, ok)
!   IF(.NOT. ok) cycle
!   atom_4_coord(1:3) = new_coord(1:3)

   IF(atom4_ang(i)(1:3) == 'ZZZ') then
    atom_4_coord(1) = -99.
    atom_4_coord(2) = -99.
    atom_4_coord(3) = -99.
    call angle_calculation(atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord, angle)

    WRITE(message_text,'(10x,7a,F7.2, a)') '    . ', label_1, ' - ', label_2, &
                                                              ' - ', label_3, &
                                                              ' = ', angle, ' deg.'
    call write_info(TRIM(message_text))

   else
    ok = .false.
    call get_label_atom_coord('ang', i,4, atom_index, ok)
    IF(.NOT. ok) cycle
    atom_4_coord(1:3) = new_coord(1:3)

    call angle_calculation(atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord, angle)

    WRITE(message_text,'(10x,9a,F7.2, a)') '    . ', label_2, ' - ', label_1, ' ^ ',               &
                                                     label_3, ' - ', label_4, ' = ',angle, ' deg.'
    call write_info(TRIM(message_text))
   endif


  END do

  return
end subroutine  calcul_angles

!-------------------------------------------------------------------------------
subroutine angle_calculation(coord_1, coord_2, coord_3, coord_4, angle)
 USE cryscalc_module,         ONLY : SP_value, debug_proc
 USE CFML_Math_General,       ONLY : acosd
 !USE CFML_Constants,          ONLY : sp
 USE CFML_GlobalDeps,         ONLY : sp

 implicit none
  real, dimension(3), intent(in)    :: coord_1, coord_2, coord_3, coord_4
  real,               intent(inout) :: angle
  real                              :: dist_12, dist_34
  real                              :: cos_angle

  if(debug_proc%level_3)  call write_debug_proc_level(3, "angle_calculation")


  if(coord_4(1) < -98. .and. coord_4(2) < -98. .and. coord_4(3) < -98.) then
    ! distance entre 1 et 2
    call distance_calculation(coord_1(:), coord_2(:), dist_12)

    ! distance entre 2 et 3
    call distance_calculation(coord_2(:), coord_3(:), dist_34)

    ! produit scalaire 21.23
    call scalar_product(coord_2(:), coord_1(:), coord_2(:), coord_3(:))
    cos_angle = SP_value / (dist_12 * dist_34)
    IF((1.-ABS(cos_angle))  < 0.00001) cos_angle=SIGN(1.,cos_angle)
    angle = ACOSd(cos_angle)

  else
    ! distance entre 1 et 2
    call distance_calculation(coord_1(:), coord_2(:), dist_12)

    ! distance entre 3 et 4
    call distance_calculation(coord_3(:), coord_4(:), dist_34)

    ! produit scalaire 21.34
    call scalar_product(coord_2(:), coord_1(:), coord_3(:), coord_4(:))
    cos_angle = SP_value / (dist_12 * dist_34)
    IF((1.-ABS(cos_angle))  < 0.00001) cos_angle=SIGN(1.,cos_angle)
    angle = ACOSd(cos_angle)

  endif

 return
end subroutine angle_calculation

!-------------------------------------------------------------------------------
subroutine volume_calculation(input_string)
 use cryscalc_module,  only       :  ON_SCREEN, pi, unit_cell, keyword_create_CIF, CIF_unit, &
                                     DC_ort, ort_DC, RC_ort, ort_RC, GMD, GMR, keyword_CELL, &
                                     known_cell_ESD, input_CELL, keyword_read_ins, message_text, &
                                     UB_mat_log, debug_proc, crystal_cell, write_details

 use CIF_module,       only       :  CIF_parameter

 USE CFML_Math_General      ,  ONLY   :  cosd, sind, acosd
 USE math_module,              ONLY   :  orthog, metric, matinv
 USE CFML_crystal_metrics,     only   :  Set_Crystal_Cell, write_crystal_cell, Volume_sigma_from_cell
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: input_string
  REAL                          :: CosAlfa, CosBeta, CosGamma
  REAL                          :: cosa2, cosb2, cosg2
  INTEGER                       :: i
  INTEGER                       :: ifail
  character (len=64)            :: uc_string
  character (len=16)            :: ESD_string
  real, parameter               :: eps = 0.00001
  character (len=16)            :: fmt_1, fmt_2
  !type (crystal_cell_type)      :: celda

  if(debug_proc%level_2)  call write_debug_proc_level(2, "volume_calculation ("//trim(input_string)//")")

  !call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)
  call create_CELL_object()


  if(ON_SCREEN .and. write_details) then
   IF(input_string=='out') then
    call write_info(' ')
    call write_info('          Crystallographic unit cell ')
    call write_info('          --------------------------')
    call write_info(' ')
   endif
  endif

  ! calcul du volume
  CosAlfa  = COSd(unit_cell%param(4))
  Cosbeta  = COSd(unit_cell%param(5))
  Cosgamma = COSd(unit_cell%param(6))

  cosa2 = CosAlfa  ** 2
  cosb2 = cosBeta  ** 2
  cosg2 = cosGamma ** 2

  IF(2 * CosAlfa * cosbeta * cosgamma - cosa2 - cosb2 - cosg2 < -1.) then
   unit_cell%param(4) = 180. - unit_cell%param(4)
   unit_cell%param(5) = 180. - unit_cell%param(5)
   unit_cell%param(6) = 180. - unit_cell%param(6)
   CosAlfa  = COSd(unit_cell%param(4))
   Cosbeta  = COSd(unit_cell%param(5))
   Cosgamma = COSd(unit_cell%param(6))
   cosa2 = CosAlfa  ** 2
   cosb2 = cosBeta  ** 2
   cosg2 = cosGamma ** 2
  endif

  if(debug_proc%level_3) call write_debug_proc_level(3, "volume (direct space)")

  unit_cell%volume = SQRT(1 + 2 * CosAlfa * cosbeta * cosgamma - cosa2 - cosb2 - cosg2)
  unit_cell%volume = unit_cell%volume * unit_cell%param(1) * unit_cell%param(2) * unit_cell%param(3)

  ! reciprocal parameters
  if(debug_proc%level_3) call write_debug_proc_level(3, "reciprocal parameters")

  unit_cell%rec_param(1) = unit_cell%param(2)*unit_cell%param(3)*sind(unit_cell%param(4))/unit_cell%volume
  unit_cell%rec_param(2) = unit_cell%param(3)*unit_cell%param(1)*sind(unit_cell%param(5))/unit_cell%volume
  unit_cell%rec_param(3) = unit_cell%param(1)*unit_cell%param(2)*sind(unit_cell%param(6))/unit_cell%volume


  unit_cell%rec_param(4) = (cosd(unit_cell%param(5))*cosd(unit_cell%param(6)) - cosd(unit_cell%param(4))) &
                         / (sind(unit_cell%param(5))*sind(unit_cell%param(6)))
  unit_cell%rec_param(5) = (cosd(unit_cell%param(4))*cosd(unit_cell%param(6)) - cosd(unit_cell%param(5))) &
                         / (sind(unit_cell%param(4))*sind(unit_cell%param(6)))
  unit_cell%rec_param(6) = (cosd(unit_cell%param(5))*cosd(unit_cell%param(4)) - cosd(unit_cell%param(6))) &
                         / (sind(unit_cell%param(5))*sind(unit_cell%param(4)))

  do i=4,6
   unit_cell%rec_param(i) = acosd(unit_cell%rec_param(i))
  end do
  if(debug_proc%level_3) call write_debug_proc_level(3, "volume (reciprocal space)")

  unit_cell%rec_volume = 1./unit_cell%volume

  !-------------------------------------------------
  Call ORTHOG(unit_cell%param(1:3),      unit_cell%param(4:6),      DC_ORT)
  Call ORTHOG(unit_cell%rec_param(1:3),  unit_cell%rec_param(4:6),  RC_ORT)
  Call METRIC(unit_cell%param(1:3),      unit_cell%param(4:6),      GMD)
  Call METRIC(unit_cell%rec_param(1:3),  unit_cell%rec_param(4:6),  GMR)
  Call MATINV(DC_ORT, ORT_DC, IFAIL)
  call MATINV(RC_ORT, ORT_DC, IFAIL)
  IF(IFAIL.EQ.1) THEN
   call write_info(' => Bad cell parameters!')
   keyword_CELL = .false.
   return
  ENDIF
  !-------------------------------------------------


  if(ON_SCREEN .and. write_details) then
  IF(input_string=='out') then
   call write_info('            >>> Direct cell parameters:  ')

   if(known_cell_ESD) then
    ESD_string = '?'
    call Create_ESD_string(1, ESD_string, fmt_1, fmt_2)
    if(ESD_string(1:1) /= '?') then
     write(uc_string, fmt_1) '                 .     a = ', unit_cell%param(1), '(', trim(ESD_string),') A'
    else
     write(uc_string, fmt_2) '                 .     a = ', unit_cell%param(1), ' A'
    end if
    call write_info(TRIM(uc_string))

    ESD_string = '?'
    call Create_ESD_string(2, ESD_string, fmt_1, fmt_2)
    if(ESD_string(1:1) /= '?') then
     write(uc_string, fmt_1) '                 .     b = ', unit_cell%param(2), '(', trim(ESD_string),') A'
    else
     write(uc_string, fmt_2) '                 .     b = ', unit_cell%param(2), ' A'
    end if
    call write_info(TRIM(uc_string))

    ESD_string = '?'
    call Create_ESD_string(3, ESD_string, fmt_1, fmt_2)
    if(ESD_string(1:1) /= '?') then
     write(uc_string, fmt_1) '                 .     c = ', unit_cell%param(3), '(', trim(ESD_string),') A'
    else
     write(uc_string, fmt_2) '                 .     c = ', unit_cell%param(3), ' A'
    end if
    call write_info(TRIM(uc_string))

    ESD_string = '?'
    call Create_ESD_string(4, ESD_string, fmt_1, fmt_2)
    if(ESD_string(1:1) /= '?') then
     write(uc_string, fmt_1) '                 .  alfa = ', unit_cell%param(4), '(', trim(ESD_string),') deg.'
    else
     write(uc_string, fmt_2) '                 .  alfa = ', unit_cell%param(4), ' deg'
    end if
    call write_info(TRIM(uc_string))

    ESD_string = '?'
    call Create_ESD_string(5, ESD_string, fmt_1, fmt_2)
    if(ESD_string(1:1) /= '?') then
     write(uc_string, fmt_1) '                 .  beta = ', unit_cell%param(5), '(', trim(ESD_string),') deg'
    else
     write(uc_string, fmt_2) '                 .  beta = ', unit_cell%param(5), ' deg'
    end if
    call write_info(TRIM(uc_string))

    ESD_string = '?'
    call Create_ESD_string(6, ESD_string, fmt_1, fmt_2)
    if(ESD_string(1:1) /= '?') then
     write(uc_string, fmt_1) '                 . gamma = ', unit_cell%param(6), '(', trim(ESD_string),') deg'
    else
     write(uc_string, fmt_2) '                 . gamma = ', unit_cell%param(6), ' deg'
    end if
    call write_info(TRIM(uc_string))

    ESD_string = '?'
    call Create_ESD_string(7, ESD_string, fmt_1, fmt_2)
    if(ESD_string(1:1) /= '?') then
     write(uc_string, fmt_1) '      . unit cell volume = ', unit_cell%volume, '(', trim(ESD_string),') A3'
    else
     call Volume_sigma_from_cell(unit_cell%param(1:3), unit_cell%param(4:6), unit_cell%param_ESD(1:3), unit_cell%param_ESD(4:6), &
                                 unit_cell%volume, unit_cell%volume_ESD)
     call Create_ESD_string(7, ESD_string, fmt_1, fmt_2)
     write(uc_string, fmt_1) '      . unit cell volume = ', unit_cell%volume, '(', trim(ESD_string),') A3'
     !write(uc_string, fmt_2) '      . unit cell volume = ', unit_cell%volume, ' A3'
    end if
    call write_info('')
    call write_info(TRIM(uc_string))

   else
    WRITE(message_text,'(a,F12.5,a)') '                 .     a = ',unit_cell%param(1),' A'
    call write_info(TRIM(message_text))
    WRITE(message_text,'(a,F12.5,a)') '                 .     b = ',unit_cell%param(2),' A'
    call write_info(TRIM(message_text))
    WRITE(message_text,'(a,F12.5,a)') '                 .     c = ',unit_cell%param(3),' A'
    call write_info(TRIM(message_text))
    WRITE(message_text,'(a,F12.5,a)') '                 .  alfa = ',unit_cell%param(4),' deg.'
    call write_info(TRIM(message_text))
    WRITE(message_text,'(a,F12.5,a)') '                 .  beta = ',unit_cell%param(5),' deg.'
    call write_info(TRIM(message_text))
    WRITE(message_text,'(a,F12.5,a)') '                 . gamma = ',unit_cell%param(6),' deg.'
    call write_info(TRIM(message_text))

    call write_info('')
    WRITE(message_text,'(a,F15.2,a)') '                 . unit cell volume : ', unit_cell%volume, ' A3'
    call write_info(TRIM(message_text))

   end if

   call write_info('')
   call write_info('            >>> Reciprocal cell parameters:  ')

   WRITE(message_text,'(a,F10.5,a)') '                 .     a* = ',unit_cell%rec_param(1),' (A-1)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .     b* = ',unit_cell%rec_param(2),' (A-1)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .     c* = ',unit_cell%rec_param(3),' (A-1)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .  alfa* = ',unit_cell%rec_param(4),' (deg-1)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .  beta* = ',unit_cell%rec_param(5),' (deg-1)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 . gamma* = ',unit_cell%rec_param(6),' (deg-1)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a)')         ' '
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 . reciprocal volume = ', unit_cell%rec_volume, ' A-3'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a)')         ' '
    call write_info(TRIM(message_text))
  end if
  end if

  !if(keyword_create_CIF .and. input_string== "out") then
  if(keyword_create_CIF) then
   write(CIF_parameter%cell_length_a, *) unit_cell%param(1)
   write(CIF_parameter%cell_length_b, *) unit_cell%param(2)
   write(CIF_parameter%cell_length_c, *) unit_cell%param(3)
   write(CIF_parameter%cell_angle_alpha, *)  unit_cell%param(4)
   write(CIF_parameter%cell_angle_beta,  *)  unit_cell%param(5)
   write(CIF_parameter%cell_angle_gamma, *)  unit_cell%param(6)
   call write_CIF_file('UNIT_CELL_INFO')
   IF(known_cell_ESD) then
    call write_CIF_file('CELL_PARAM_ESD')
   else
    call write_CIF_file('CELL_PARAM')
   endif
  endif


  if(UB_mat_log .And. input_string=='out')  call write_UB_matrix

  !if(input_string =='out') call write_crystal_cell(crystal_cell)

 return
end subroutine volume_calculation

!---------------------------------------------------------------------------------------------
 subroutine Get_esd_string(esd_real, esd_string)
 ! transforme esd type 0.000x en (x)

  implicit none
   real, intent(inout)             :: esd_real
   character (len=64), intent(out) :: esd_string
   character (len=32)              :: string
   integer                         :: esd_int
   real, parameter                 :: eps = 0.00001
   integer                         :: i1


   !do
   !write(*,*) 'esd_real : ', esd_real, '   ', int(esd_real),  '   ', nint(esd_real)
   ! if(esd_real - int(esd_real) < eps) exit
!    esd_real = 10.*esd_real
!   end do
!   write(esd_string, '(I12)') INT(esd_real)

    write(string, '(F10.4)') esd_real
    i1 = index(string, '.')
    string = string(i1+1:)
    read(string, *) esd_int
    write(esd_string, '(I12)') esd_int

 return
end subroutine Get_esd_string

!------------------------------------------------------------------------------
subroutine write_UB_matrix
 USE cryscalc_module,          only  : tmp_unit, UB_matrix, unit_cell, keyword_cell, message_text, debug_proc
 USE CFML_Geometry_SXTAL,      ONLY  : CELL_fr_UB
 USE IO_module


 implicit none
  integer               :: i, i1, ier
  character (len=256)   :: read_line

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_UB_matrix")


  write (message_text, '(a,3(1x,F10.5))') '  UB_matrix 11 12 13: ', UB_matrix(1,1),   UB_matrix(1,2),  UB_matrix(1,3)
  call write_info(trim(message_text))
  write (message_text, '(a,3(1x,F10.5))') '  UB_matrix 21 22 23: ', UB_matrix(2,1),   UB_matrix(2,2),  UB_matrix(2,3)
  call write_info(trim(message_text))
  write (message_text, '(a,3(1x,F10.5))') '  UB_matrix 31 32 33: ', UB_matrix(3,1),   UB_matrix(3,2),  UB_matrix(3,3)
  call write_info(trim(message_text))

  open(unit=tmp_unit, file='tmp_ub.txt')
   call CELL_fr_UB(UB_matrix, tmp_unit)
  close(unit=tmp_unit)

  open(unit=tmp_unit, file='tmp_ub.txt')
  do i=1,5
   read(unit=tmp_unit, fmt='(a)', iostat=ier) read_line
   if(ier ==0) then
    call write_info(trim(read_line))
   end if
   if(i==3) then    ! lecture des parametres de maille (directs) sur la ligne 3
    i1 = index(read_line, ':')
    if(i1 /=0) then
     read(read_line(31+1:), *, iostat=ier) unit_cell%param(1:6)
     if(ier == 0)  keyword_cell = .true.
    end if
   end if
  end do
  close(unit=tmp_unit)
  call system('del tmp_ub.txt')



 return
end subroutine write_UB_matrix

!------------------------------------------------------------------------------
subroutine get_crystal_SIZE_min_max()
 USE cryscalc_module, ONLY : crystal, debug_proc
  implicit none
   INTEGER               :: i, i_mid
   INTEGER, DIMENSION(1) :: i_min, i_max

   if(debug_proc%level_2)  call write_debug_proc_level(2, "get_crystal_size_min_max")


   i_max = MAXLOC(crystal%size)
   i_min = MINLOC(crystal%size)
   do i = 1, 3
    if (i_max(1) == i) cycle
    if (i_min(1) == i) cycle
    i_mid = i
   end do

   crystal%size_min = MINVAL(crystal%SIZE)
   crystal%size_max = MAXVAL(crystal%SIZE)
   crystal%size_mid = crystal%SIZE(i_mid)

  RETURN
end subroutine get_crystal_SIZE_min_max

!-------------------------------------------------------------------------------
subroutine crystal_volume_calculation(input_string)
 USE cryscalc_module,   ONLY : ON_SCREEN, crystal, pi, keyword_create_CIF, CIF_unit, keyword_SIZE, message_text, &
                               debug_proc
 USE macros_module,     ONLY : l_case
 USE IO_module

 implicit none
   CHARACTER (len=*), intent(in) :: input_string

   if(debug_proc%level_2)  call write_debug_proc_level(2, "crystal_volume_calculation")

   if(l_case(input_string(1:3)) == 'out') then
    call write_info(' ')
    call write_info('          Crystal dimensions ')
    call write_info('          ------------------')
    call write_info(' ')
   endif

   crystal%volume =  crystal%size(1)* crystal%size(2)* crystal%size(3)
   crystal%radius =  ((3.*crystal%volume) / (4.*pi) )  ** (1./3.)

   if(l_case(input_string(1:3)) == 'out') then
    call write_info('')
    write(message_text,'(10x,a,3(1x,F6.3))')'  > Crystal size (mm)   : ',  crystal%size(1:3)
    call write_info(TRIM(message_text))
    call write_info('')
    write(message_text,'(a,F8.5)')      '      . Crystal volume (mm3)             : ',  &
                                       crystal%size(1)* crystal%size(2)* crystal%size(3)
    call write_info(TRIM(message_text))
    crystal%radius = ((3.*crystal%volume) / (4.*pi) )  ** (1./3.)
    write(message_text,'(a,F8.5)')      '      . Crystal equiv. radius <Req> (mm) : ', crystal%radius
    call write_info(TRIM(message_text))
    call write_info('')
   end if

   IF(keyword_create_CIF)  call write_CIF_file('CRYSTAL_INFORMATION')

  return
END subroutine crystal_volume_calculation

!-------------------------------------------------------------------------------
subroutine calcul_dhkl
 use cryscalc_module, only         : pi, unit_cell, wavelength, keyword_WAVE,  keyword_SPGR,          &
                                     cell_star, cos_angle_star,                                       &
                                     nb_hkl, H, shift_2theta, keyword_QVEC, Qvec, message_text, SPG,  &
                                     debug_proc

 USE CFML_Math_General,              ONLY : sind
 USE CFML_Reflections_Utilities,     ONLY : HKL_Absent
 USE CFML_crystallographic_symmetry, ONLY : Space_Group_Type, set_spacegroup

 USE IO_module
 implicit none

 integer                              :: i
 real                                 :: Q_hkl, d_hkl, stl_hkl
 real                                 :: Z, angle_2theta
 INTEGER                              :: nb_ref, n
 REAL, DIMENSION(3)                   :: HQ
 CHARACTER (len=1)                    :: ind_H
 logical                              :: absent, angle_OK

 if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_dhkl")

 IF(unit_cell%volume < 0.1) call volume_calculation('out')

 call write_info('')
 call write_info('    . d_hkl CALCULATIONS')
 call write_info('      ------------------')
 call write_info('')


 IF(keyword_WAVE) then
  write(message_text, '(5x,a,F10.5,a)') 'Wave = ', wavelength, ' A'
  call write_info(trim(message_text))
  call write_info('')
  call write_info('        h      k      l               d_hkl(A)   stl_hkl(A-1)     Q_hkl(A-1)   theta(deg.)    2theta(deg.)')
  call write_info(' ')
 else
  call write_info('        h      k      l               d_hkl(A)   stl_hkl(A-1)     Q_hkl(A-1)')
  call write_info(' ')
 endif

 call calcul_cell_star

 !a_star = unit_cell%param(2)*unit_cell%param(3) * sind(unit_cell%param(4)) / unit_cell%Volume
 !b_star = unit_cell%param(3)*unit_cell%param(1) * sind(unit_cell%param(5)) / unit_cell%Volume
 !c_star = unit_cell%param(1)*unit_cell%param(2) * sind(unit_cell%param(6)) / unit_cell%Volume

 !alfa_rad = unit_cell%param(4) * pi/180.
 !beta_rad = unit_cell%param(5) * pi/180.
 !gama_rad = unit_cell%param(6) * pi/180.

 !cos_alfa_star = (cos(beta_rad) * cos(gama_rad) - cos(alfa_rad) ) / (sin(beta_rad) * sin(gama_rad))
 !cos_beta_star = (cos(gama_rad) * cos(alfa_rad) - cos(beta_rad) ) / (sin(gama_rad) * sin(alfa_rad))
 !cos_gama_star = (cos(alfa_rad) * cos(beta_rad) - cos(gama_rad) ) / (sin(alfa_rad) * sin(beta_rad))

 IF(keyword_QVEC) then
  nb_ref   = 3*nb_hkl
 else
  nb_ref   = nb_hkl
 END if


 n = 1
 do i = 1, nb_ref

  if(keyword_SPGR) then
   HQ(:) = H(: ,i)
   absent=HKL_Absent(HQ(:), SPG)
   if(absent) then
    write(message_text, '(5x,   3(1x,F6.2),5xa)')   HQ(1:3), ' <<< not allowed !'
    call write_info(trim(message_text))
    cycle
   endif
  end if

  IF(keyword_QVEC) then
   IF(i == 1+3*(n-1)) then
    IF(n/=1) call write_info(' ')
    HQ(:) = H(: ,n)                  ! reflection principale
    ind_H = ' '
   ELSEIF(i == 2 + 3*(n-1)) then
    HQ(:) = H(: ,n) - qvec(:)        ! sat -
    ind_H = '-'
   ELSEIF(i == 3 + 3*(n-1)) then
    HQ(:) = H(: ,n) + qvec(:)        ! sat +
    n = n + 1
    ind_H = '+'
   endif
  else
   HQ(:) = H(: ,i)
  END if

  call Calcul_Q_d_stl(HQ, Q_hkl, d_hkl, stl_hkl)

  !Q_hkl = HQ(1)**2 * a_star**2 + HQ(2)**2 * b_star**2 + HQ(3)**2 * c_star**2              &
  !      + 2*HQ(1)*HQ(2) * a_star*b_star*cos_gama_star                                      &
  !      + 2*HQ(2)*HQ(3) * b_star*c_star*cos_alfa_star                                      &
  !      + 2*HQ(3)*HQ(1) * c_star*a_star*cos_beta_star

  !d_hkl = 1/sqrt(Q_hkl)
  !stl_hkl = 1/(2*d_hkl)
  !Q_hkl = 4*pi*stl_hkl

  IF(keyword_WAVE) then
   Z = wavelength / (2 * d_hkl)
   IF (Z**2 < 1.) then
    angle_2theta = 2 * ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
    angle_ok = .true.
   else
    angle_2theta = 180.
    angle_OK = .false.
   endif
   if(angle_OK) then
   angle_2theta = angle_2theta + shift_2theta(1) + shift_2theta(2)*COS(Angle_2theta*pi/180) +  &
                                                   shift_2theta(3)*SIN(Angle_2theta*pi/180)
   end if
   if(keyword_QVEC) then
    if(angle_OK) then
    write(message_text, '(4x,a1,3(1x,F6.2),5x,5(5x,F10.4))')      ind_H(:), HQ(1:3), d_hkl, stl_hkl, Q_hkl, &
                                                                  angle_2theta/2., angle_2theta
    else
    write(message_text, '(4x,a1,3(1x,F6.2),5x,3(5x,F10.4),5x,a)') ind_H(:), HQ(1:3), d_hkl, stl_hkl, Q_hkl, &
                                                                  'out of angular range !!'
    end if
   else
    if(angle_OK) then
    write(message_text, '(5x,   3(1x,F6.2),5x,5(5x,F10.4))')         HQ(1:3), d_hkl, stl_hkl, Q_hkl, angle_2theta/2., angle_2theta
    else
    write(message_text, '(5x,   3(1x,F6.2),5x,3(5x,F10.4),5x,a)')    HQ(1:3), d_hkl, stl_hkl, Q_hkl,  ' out of angular range !!'
    end if
   end if
   call write_info(TRIM(message_text))
  else
   if(keyword_QVEC) then
    write(message_text, '(4x,a1,3(1x,F6.2),5x,3(5x,F10.4))') ind_H(:), HQ(1:3), d_hkl, stl_hkl, Q_hkl
   else
    write(message_text, '(5x,   3(1x,F6.2),5x,3(5x,F10.4))')           HQ(1:3), d_hkl, stl_hkl, Q_hkl
   endif
    call write_info(TRIM(message_text))
  endif


 end do
 call write_info('')

 return
end subroutine calcul_dhkl

!----------------------------------------------------------
subroutine Calcul_cell_star
 USE Cryscalc_module,           only : unit_cell, cell_star, cos_angle_star, pi
 USE CFML_Math_General,         ONLY : sind

 implicit none
  real :: alfa_rad, beta_rad, gamma_rad


  cell_star(1) = unit_cell%param(2)*unit_cell%param(3) * sind(unit_cell%param(4)) / unit_cell%Volume
  cell_star(2) = unit_cell%param(3)*unit_cell%param(1) * sind(unit_cell%param(5)) / unit_cell%Volume
  cell_star(3) = unit_cell%param(1)*unit_cell%param(2) * sind(unit_cell%param(6)) / unit_cell%Volume

  alfa_rad  = unit_cell%param(4) * pi/180.
  beta_rad  = unit_cell%param(5) * pi/180.
  gamma_rad = unit_cell%param(6) * pi/180.

  cos_angle_star(1) = (cos(beta_rad)  * cos(gamma_rad) - cos(alfa_rad) )   / (sin(beta_rad)  * sin(gamma_rad))
  cos_angle_star(2) = (cos(gamma_rad) * cos(alfa_rad)  - cos(beta_rad) )  / (sin(gamma_rad) * sin(alfa_rad))
  cos_angle_star(3) = (cos(alfa_rad)  * cos(beta_rad)  - cos(gamma_rad) ) / (sin(alfa_rad)  * sin(beta_rad))

  return
end subroutine calcul_cell_star

!------------------------------------------------------------

subroutine Calcul_Q_d_stl(ref_H, Q_hkl, d_hkl, stl_hkl)
 USE Cryscalc_module, only : cell_star, cos_angle_star, pi
 implicit none
  real, dimension(3), intent(in)    :: ref_H
  real,               intent(inout) :: Q_hkl, d_hkl, stl_hkl

  Q_hkl = ref_H(1)**2 * cell_star(1)**2 + ref_H(2)**2 * cell_star(2)**2 + ref_H(3)**2 * cell_star(3)**2       &
        + 2*ref_H(1)*ref_H(2) * cell_star(1)*cell_star(2)*cos_angle_star(3)                        &
        + 2*ref_H(2)*ref_H(3) * cell_star(2)*cell_star(3)*cos_angle_star(1)                        &
        + 2*ref_H(3)*ref_H(1) * cell_star(3)*cell_star(1)*cos_angle_star(2)

  d_hkl = 1/sqrt(Q_hkl)
  stl_hkl = 1/(2*d_hkl)
  Q_hkl = 4*pi*stl_hkl

  return
end subroutine Calcul_Q_d_stl

!----------------------------------------------------------
subroutine atomic_density_calculation()
 USE cryscalc_module, ONLY : ON_SCREEN, nb_atoms_type, SFAC, nb_at, unit_cell, message_text, &
                             debug_proc
 USE IO_module

 implicit none
 !local variables
 REAL                                      :: density
 INTEGER                                   :: i

 if(debug_proc%level_2)  call write_debug_proc_level(2, "atomic_density_calculation")

 if(ON_SCREEN ) then
  call write_info(' ')
  call write_info('          Atomic density ')
  call write_info('          --------------')
  call write_info(' ')

  write(message_text,'(a)')      '    Atom     nb atoms/cm3 E22'
   call write_info(TRIM(message_text))
  write(message_text,'(a)')      ' '
   call write_info(TRIM(message_text))
 endif


! calcul du nombre d'atomes i/cm3
 density = 0.
 do i = 1, nb_atoms_type
  nb_at(i) = SFAC%number(i) / unit_cell%volume * 1.E24

  if(ON_SCREEN) then
   write(message_text, '(a8,10x,F9.2)') TRIM(SFAC%type(i)),  nb_at(i)/1.e22
   call write_info(TRIM(message_text))
  endif

  density  = density +  nb_at(i)
 end do

 if(ON_SCREEN) then
  call write_info('')
  WRITE(message_text,'(5x,a,F9.2,a)') "           => Atomic density = ", density/1.e22,' E22 atoms/cm3'
  call write_info(TRIM(message_text))
  call write_info('')
 endif


 return
END subroutine atomic_density_calculation

!------------------------------------------------------------------------------
subroutine atomic_identification()
 USE atomic_data
 USE cryscalc_module, only          : nb_atoms_type, SFAC, num_atom, known_atomic_label, known_atomic_features, &
                                      debug_proc
 USE macros_module, ONLY            : u_case
 USE IO_module

 implicit none
  INTEGER                          :: i, j

  if(debug_proc%level_3)  call write_debug_proc_level(3, "atomic_identification")

  if(.not. known_atomic_label)    then
   call definition_atomic_label    ! %symbol, %name
   known_atomic_label = .true.
  endif
  if(.not. known_atomic_features) then
   call def_atomic_features        ! %weight, ...
   known_atomic_features = .true.
  endif


! reconnaissance des atomes
 do i=1, nb_atoms_type
  SFAC%type(i) = ADJUSTL(SFAC%type(i))
  !
  j=0
  do
    j=j+1
    if (j>201) then
     call write_info('')
     call write_info('   !!! '//TRIM(SFAC%type(i))// ': incorrect atomic symbol !!')
     return
    end if

    if (u_case(SFAC%type(i)(1:))==atom(j)%symbol(1:)) exit

  end do
  Num_atom(i) = j

 end do


 return
end subroutine atomic_identification


!------------------------------------------------------------------------------
subroutine molecular_weight()
 USE atomic_data
 USE cryscalc_module, only          : ON_SCREEN, write_details, nb_atoms_type, SFAC, sto, num_atom, Z_unit,              &
                                      keyword_CHEM, neutrons, keyword_create_CIF, CIF_unit, message_text, &
                                      molecule, nb_atom_no_H, debug_proc
 USE IO_module,      only           : write_info
 USE macros_module,  only           : l_case

 implicit none
  INTEGER                                     :: i
  REAL,              DIMENSION(nb_atoms_type) :: Weight_percent, atomic_percent
  CHARACTER(LEN=16), DIMENSION(nb_atoms_type) :: labl
  REAL                                        :: mol_weight, Total_sto
  INTEGER                                     :: n_electrons
  CHARACTER (LEN=16)                          :: fmt_

  if(debug_proc%level_2)  call write_debug_proc_level(2, "molecular_weight")

 if(ON_SCREEN .and. write_details) then
 call write_info(' ')
 call write_info('          Molecular features ')
 call write_info('          ------------------')
 endif

 labl = '?'

 do i=1, nb_atoms_type

  ! molecular formula
  IF(LEN_TRIM(SFAC%type(i))==2) SFAC%type(i)(2:2) = l_case(SFAC%type(i)(2:2))

  IF(.NOT. keyword_CHEM)  then
   sto(i) = SFAC%number(i) / Z_unit
  !else
  ! sto(i) = SFAC%number(i)
  endif

  IF(INT(sto(i)) < 10) then
   if (ABS(INT(sto(i))-sto(i)) < 0.01) then
    WRITE(labl(i), '(a,I1)') trim(SFAC%type(i)), INT(sto(i))
   else
    WRITE(labl(i), '(a,F5.2)') trim(SFAC%type(i)), sto(i)
   endif
  ELSEIF(INT(sto(i)) < 100) THEN
   if (ABS(INT(sto(i))-sto(i)) < 0.01) then
    WRITE(labl(i), '(a,I2)') TRIM(SFAC%type(i)), INT(sto(i))
   else
    WRITE(labl(i), '(a,F6.2)') trim(SFAC%type(i)), sto(i)
   endif
  ELSEIF(INT(sto(i)) < 1000) THEN
   if (ABS(INT(sto(i))-sto(i)) < 0.01) then
    WRITE(labl(i), '(a,I3)') TRIM(SFAC%type(i)), INT(sto(i))
   else
    WRITE(labl(i), '(a,F7.2)') trim(SFAC%type(i)), sto(i)
   endif
  endif

 end do


 ! molecular weight et nombre d'electrons
 Mol_Weight  = 0.
 Total_sto   = 0.
 n_electrons = 0
 do  i=1, nb_atoms_type
  Mol_Weight = Mol_Weight + atom(Num_atom(i))%weight * sto(i)
  total_sto = total_sto + sto(i)
  n_electrons = n_electrons + atom(Num_atom(i))%Z * sto(i)

  if(len_trim(SFAC%type(i)) == 1 .and. SFAC%type(i)(1:1) == 'H') cycle
  nb_atom_no_H = nb_atom_no_H + sto(i)*molecule%Z_unit

 end do
 molecule%weight = mol_weight
 molecule%Z      = n_electrons




 ! weight percentage
 do i=1, nb_atoms_type
  Weight_percent(i) = 100* (atom(Num_atom(i))%weight * sto(i))/Mol_weight
  Atomic_percent(i) = 100* (sto(i)/total_sto)
 end do

 if(nb_atoms_type < 10) then
  write(fmt_, '(a,i1,a)') '(', nb_atoms_type, '(1x,a))'
 else
  write(fmt_, '(a,i2,a)') '(', nb_atoms_type, '(1x,a))'
 endif
 if(nb_atoms_type /=0)  write(molecule%formula, fmt_) (TRIM(labl(i)),i=1,nb_atoms_type)

 if(ON_SCREEN .and. write_details) then
 call write_info('')
 WRITE(message_text, '(2a)')       '   >> Molecular formula: ',trim(molecule%formula)
 call write_info(TRIM(message_text))

 call write_info('')
 WRITE(message_text,'(a,F8.2)') '   >> Molecular weight: ',Molecule%weight
  call write_info(TRIM(message_text))
 call write_info('')
 WRITE(message_text, '(a,I6)')       '   >> Total number of electrons in the molecule: ', molecule%Z
  call write_info(TRIM(message_text))

 write(message_text,'(a)')      ' '
  call write_info(TRIM(message_text))
 write(message_text,'(a)')      '    Atom       Z      Atomic weight    Stoechiometry        % atomic        % weight'
  call write_info(TRIM(message_text))
 write(message_text,'(a)')      ' '
  call write_info(TRIM(message_text))
 do i=1, nb_atoms_type
  write(message_text, '(a8,4x,I4,11x,F8.4,11x,F6.2, 2(10x,F6.2))') TRIM(SFAC%type(i)),  atom(Num_atom(i))%Z, &
                                                                   atom(Num_atom(i))%weight, &
                                                                   sto(i), Atomic_percent(i), Weight_percent(i)
  call write_info(TRIM(message_text))
 end do
 endif

 IF(keyword_create_CIF)  then
  call write_CIF_file('CHEMICAL_INFORMATION')
  call write_CIF_file('CHEMICAL_2')
 END IF

 ! calcul du F000
 !if (.not. neutrons) call F000_calculation('X')


 return
end subroutine molecular_weight

!------------------------------------------------------------------------------
subroutine density_calculation()
 USE cryscalc_module, ONLY        : ON_SCREEN, write_details, molecule, Z_unit, unit_cell, nb_atom_no_H, molecule, &
                                    keyword_create_CIF, CIF_unit, message_text, debug_proc
 USE IO_module

 implicit none
  REAL, parameter                 :: avogadro = 6.022

  if(debug_proc%level_2)  call write_debug_proc_level(2, "density_calculation")

  !molecule%density = (Z_unit * Molecule%weight) / (0.1* avogadro * unit_cell%volume  )
  molecule%density = (molecule%Z_unit * Molecule%weight) / (0.1* avogadro * unit_cell%volume  )

  if(ON_SCREEN .and. write_details) then
  call write_info('')
  WRITE(message_text,'(1x,a,F9.3)')  "   >> Density (g/cm3)                      = ", molecule%density
  call write_info(TRIM(message_text))
  call write_info('')

  if(nb_atom_no_H /=0) then
   WRITE(message_text,'(1x,a,I6)')   "   >> Number of no H atom in the unit cell = ", nb_atom_no_H
   call write_info(TRIM(message_text))
   WRITE(message_text,'(1x,a,F9.3)') "   >> Mean no H atomic volume (A3)         = ", unit_cell%volume / nb_atom_no_H
   call write_info(TRIM(message_text))
   call write_info('')

  end if

  end if

  IF(keyword_create_CIF)  call write_CIF_file('CRYSTAL_DENSITY')


 RETURN
end subroutine density_calculation

!------------------------------------------------------------------------------
subroutine calcul_barycentre
 use cryscalc_module
 USE IO_module

  implicit none
   integer                      :: i , j , f, t
   real, dimension(3)           :: bary_coord, inside_coord
   integer, dimension(500)      :: bary_atom    ! numero des atomes dans le calcul du barycentre
   integer, dimension(3)        :: move


 if(debug_proc%level_2)  call write_debug_proc_level(2, "barycentre")

 do i = 1, nb_bary_calc


  ! quels atomes
  do j= 1, nb_atom
    do f=1,nb_atom_bary(i)
     if (atom_bary(i,f) == atom_label(j)) then
      bary_atom(f) = j
     end if
    end do
  end do

  do f=1, nb_atom_bary(i)
   if(bary_atom(f) == 0) then
    call write_info('')
    write(message_text, '(10x, 3a)') '>> ', trim(atom_bary(i,f)), ' : wrong input label !!'
    call write_info(trim(message_text))
    call write_info('')
    nb_bary_calc = nb_bary_calc - 1
    return
   end if
  end do



  bary_coord(1) = sum(atom_coord(1, bary_atom(1:nb_atom_bary(i)))) /   nb_atom_bary(i)
  bary_coord(2) = sum(atom_coord(2, bary_atom(1:nb_atom_bary(i)))) /   nb_atom_bary(i)
  bary_coord(3) = sum(atom_coord(3, bary_atom(1:nb_atom_bary(i)))) /   nb_atom_bary(i)


  call write_info('')
  write(message_text,'(a,I3)')      '  ==> Centroid #: ', i
   call write_info(TRIM(message_text))
  call write_info('')
  do j=1, nb_atom_bary(i)
   write(message_text,'(5x,2a,3F10.5, 5x, a,3F10.5)')   atom_label(bary_atom(j))(1:4), ': '  ,  &
                                                        atom_coord(1,bary_atom(j)) , &
                                                        atom_coord(2,bary_atom(j)) , &
                                                        atom_coord(3,bary_atom(j)) , 'v = ',  &
                                                        bary_coord(1)-atom_coord(1,bary_atom(j)), &
                                                        bary_coord(2)-atom_coord(2,bary_atom(j)), &
                                                        bary_coord(3)-atom_coord(3,bary_atom(j))

   call write_info(TRIM(message_text))
  end do
  call write_info('')
  write(message_text,'(5x,a,3F10.5)')   'centroid coordinates: ', bary_coord(1:3)
  call write_info(TRIM(message_text))

  inside_coord = bary_coord
  move = 0
  t    = 0
  do j=1, 3
   do while (inside_coord(j) < 0.0)
    move(j) = move(j) + 1
    inside_coord(j) = inside_coord(j)  + 1
    t = t + 1
   end do
   do while (inside_coord(j) > 1.0)
    move(j) = move(j) - 1
    inside_coord(j) = inside_coord(j)  - 1
    t = t + 1
   end do
  end do

  if(t /=0) then
   call write_info('')
   write(message_text, '(a,3(2x,i2))') '  SHELX instruction : MOVE ', move(1:3)
   call write_info(trim(message_text))
  endif
 end do

 return
end subroutine calcul_barycentre

!------------------------------------------------------------------------------
subroutine calcul_Tolman_cone_angle
 use cryscalc_module,  only : nb_tolman_calc, nb_atom, atom_label, atom_typ, atom_coord, atom_tolman, r_vdw, &
                              modif_rvdw, message_text, debug_proc, mendel_atom_nb, known_atomic_label, pi, write_details
 USE macros_module,    only : u_case
 USE CFML_Scattering_Chemical_Tables, ONLY : set_chem_info, chem_info, num_chem_info
 use atomic_data
 USE IO_module

  implicit none
   integer                      :: i , j , f, k, l, m, long
   real, dimension(3)           :: virtual_coord, theta
   integer, dimension(500)      :: tolman_atom    ! numero des atomes dans le calcul du barycentre
   real, dimension(3)           :: distance, rvdw, angle_rvdw
   real                         :: angle, sum_angles, Tolman_angle
   CHARACTER (LEN=4)            :: symb_car
   LOGICAL                      :: ok



 if(debug_proc%level_2)  call write_debug_proc_level(2, "Tolman angle")



 do i = 1, nb_tolman_calc

  ! quels atomes
  do j= 1, nb_atom
    do f=1,5
     if (u_case(atom_tolman(i,f)) == u_case(atom_label(j))) then
      tolman_atom(f) = j
     end if
    end do
  end do

  if(.not. known_atomic_label)  then
   call definition_atomic_label    ! %symbol, %name
   known_atomic_label = .true.
  endif
  CALL set_chem_info

  if(write_details) then
  !write(message_text, '(10x,3a,i2,3F10.5)') ' Metal atom    : ', atom_tolman(nb_tolman_calc, 1), atom_typ(Tolman_atom(1)), Tolman_atom(1), atom_coord(1:3,Tolman_atom(1))
  write(message_text, '(10x,3a,3F10.5)') ' Metal atom    : ', atom_tolman(nb_tolman_calc, 1), atom_typ(Tolman_atom(1)),  &
                                                              atom_coord(1:3,Tolman_atom(1))
  call write_info(TRIM(message_text))
  !write(message_text, '(10x,3a,i2,3F10.5)') ' Centered atom : ', atom_tolman(nb_tolman_calc, 2), atom_typ(Tolman_atom(2)), Tolman_atom(2), atom_coord(1:3,Tolman_atom(2))
  write(message_text, '(10x,3a,3F10.5)') ' Centered atom : ', atom_tolman(nb_tolman_calc, 2), atom_typ(Tolman_atom(2)), &
                                                              atom_coord(1:3,Tolman_atom(2))
  call write_info(TRIM(message_text))
  end if

  do k=1, 3
   m = 2 + k
   call verif_atomic_symbol(atom_typ(Tolman_atom(m)), symb_car, ok)
   if(.not. ok) then
    call verif_atomic_number(atom_typ(Tolman_atom(m)), symb_car, ok)
    if(ok) call verif_atomic_symbol(atom_typ(Tolman_atom(m)), symb_car, ok)
    IF(.not. ok) return
   end if
   do l=1, Num_Chem_Info
    if (u_case(symb_car(1:))== u_case(chem_info(l)%symb(1:))) exit
   end do
   !write(message_text, '(10x,a,I2,3a, i2,3F10.5,5x,a, F8.4)')  ' Ligand ', k, '     : ', atom_tolman(nb_tolman_calc, m), &
   !                                                                                      atom_typ(Tolman_atom(m)),       &
   !                                                                                       Tolman_atom(m),                 &
   !                                                                                      atom_coord(1:3,Tolman_atom(m)), &
   !                                                                               '        R_vdw  = ', chem_info(l)%RWaals
   if(r_vdw(i,k) < 0.) then
    rvdw(k) = chem_info(l)%RWaals
   else
    rvdw(k) = r_vdw(i,k)
   end if

   if(.not. modif_rvdw) then
    long = len_trim(atom_typ(Tolman_atom(m)))
    if(long == 1) then
     if(atom_typ(Tolman_atom(m))(1:1) == 'H') rvdw(k) = 1.2
    end if
   end if

   if(write_details) then
   write(message_text, '(10x,a,I2,3a, 3F10.5,5x,a, F6.2)')  ' Ligand ', k, '     : ', atom_tolman(nb_tolman_calc, m), &
                                                                                      atom_typ(Tolman_atom(m)),       &
                                                                                      atom_coord(1:3,Tolman_atom(m)), &
                                                                                  '   R_vdw  = ', rvdw(k)
   call write_info(TRIM(message_text))
   end if
   !write(message_text, '(10x,a,F8.4)')       '        R_vdw  = ', chem_info(l)%RWaals
   !call write_info(TRIM(message_text))
   !rvdw(k) = chem_info(l)%RWaals
  end do




  do f=1, 5
   if(tolman_atom(f) == 0) then
    call write_info('')
    write(message_text, '(10x, 3a)') '>> ', trim(atom_tolman(i,f)), ' : wrong input label !!'
    call write_info(trim(message_text))
    call write_info('')
    nb_tolman_calc = nb_tolman_calc - 1
    return
   end if
  end do

  !call write_info('')
  !call write_info('      ... Distance M - X  ...')
  !call distance_calculation(atom_coord(1:3,Tolman_atom(1)) , atom_coord(1:3,Tolman_atom(2)) , distance(1))
  !write(message_text, '(10x, a, F10.5,a)') " distance M - X = ", distance(1), " A"
  !call write_info(trim(message_text))

  if(write_details) then
  call write_info('')
  call write_info('      ... Distances metal - ligand  ...')
  end if
  do j=1,3
   call distance_calculation(atom_coord(1:3,Tolman_atom(1)) , atom_coord(1:3,Tolman_atom(2+j)) , distance(j))
   if(write_details) then
   write(message_text, '(10x, a, I2, a, F10.5,a)') " distance Metal - ligand ", j, " = ", distance(j), " A"
   call write_info(trim(message_text))
   end if
  end do

  virtual_coord(1) = -99.
  virtual_coord(2) = -99.
  virtual_coord(3) = -99.
  sum_angles = 0.
  if(write_details) then
  call write_info('')
  call write_info('      ... Angles ligand - metal - X ...')
  end if
  do j=1, 3
   call angle_calculation(atom_coord(1:3,Tolman_atom(2+j)), atom_coord(1:3,Tolman_atom(1)), atom_coord(1:3, Tolman_atom(2)), &
                          virtual_coord, angle)
   if(write_details) then
   write(message_text, '(10x, a, i1,a, F8.2,a)')      " Angle (  X - M - L_", j, ")       = ", angle, " deg."
   call write_info(trim(message_text))
   end if
   theta(j) = angle
   sum_angles = sum_angles + angle
   angle_rvdw(j) = 180./pi * asin(rvdw(j)/distance(j))
   theta(j) = theta(j) + angle_rvdw(j)

   if(write_details) then
   write(message_text, '(10x, a, i1,a, i1,a,F8.2,a)') " Angle (L_", j," - M - L_", j, "_VdW )  = ", angle_rvdw(j), " deg."
   call write_info(trim(message_text))

   write(message_text, '(10x, a, i1,a, F8.2,a)')      " Angle (  X - M - L_", j, "_VdW )  = ", theta(j), " deg."
   call write_info(trim(message_text))
   call write_info('')
   end if
  end do



  !call write_info('')
  !write(message_text, '(10x, a, F8.2,a)') " >>> Sum of angles     = ", sum_angles, " deg."
  !call write_info(trim(message_text))
  if(write_details) then
  call write_info('')
  write(message_text, '(10x, a, F8.2,a)') " >>> Cone angle = ", 2.*sum_angles/3., " deg."
  call write_info(trim(message_text))
  write(message_text, '(10x, a)') "     [2 * <angle X-M-L>]"
  call write_info(trim(message_text))
  end if



  Tolman_angle = sum(theta(1:3)) / 3
  Tolman_angle = 2. * Tolman_angle
  call write_info('')
  write(message_text, '(10x, a, F8.2,a)') " >>> Tolman cone angle = ", Tolman_angle, " deg."
  call write_info(trim(message_text))
  write(message_text, '(10x, a)') "     [2 * <angle X-M-L_vdw>]"
  call write_info(trim(message_text))

 end do

 return
end subroutine calcul_Tolman_cone_angle

!----------------------------------------------------------------
subroutine Calcul_plane
  use cryscalc_module, only : atom1_plane,  atom2_plane,  atom3_plane, nb_plane, &
                              atom_plane, atom_plane_nb, new_coord, message_text
  use CFML_math_3D,    only : Get_plane_from_points, Get_plane_from_3points, ERR_Math3D, ERR_Math3D_Mess
  USE CFML_GlobalDeps, ONLY : cp
  USE io_module
  implicit none
   integer                 :: i, j
   real, dimension(3,20)   :: atom_coord
   real (kind=cp), dimension(4)      :: plane, SPlane
   CHARACTER(len=16), dimension(20)  :: label
   integer, dimension(20)  :: atom_plane_index
   logical                 :: OK
   real                    :: A, B, C, D

   do i=1, nb_plane

    do j=1, atom_plane_nb(i)
     label(j) = atom_plane(j, i)
     ok = .false.
     call get_label_atom_coord('plane', i, j, atom_plane_index(j), ok)
     if(.not. ok) exit
     atom_coord(1:3, j) = new_coord
    end do
    if (.not. ok) cycle

    do j=1, atom_plane_nb(i)
     write(message_text,'(8x,a4,a,3F10.5)') trim(label(j)), ': ', atom_coord(1,j), atom_coord(2,j), atom_coord(3,j)
     call write_info(trim(message_text))
    end do

    if(atom_plane_nb(i) == 3) then
     call Get_plane_from_3points(atom_coord(i,1), atom_coord(i,2), atom_coord(i,3), A, B, C, D)
     call write_info('')
     write(message_text,'(5x,a)') ' >> Plane equation : Ax + By + Cz + D = 0'
     call write_info(trim(message_text))
     write(message_text,'(10x,a,F10.5)')  'A = ', A
     call write_info(trim(message_text))
     write(message_text,'(10x,a,F10.5)')  'B = ', B
     call write_info(trim(message_text))
     write(message_text,'(10x,a,F10.5)')  'C = ', C
     call write_info(trim(message_text))
     write(message_text,'(10x,a,F10.5)')  'D = ', D
     call write_info(trim(message_text))

    else
     call Get_plane_from_points(atom_plane_nb(i), atom_coord, Plane, SPlane)
     if(ERR_Math3D) then
      call write_info('')
      call write_info(trim(ERR_Math3D_Mess))
      call write_info('')
      return
     end if
     call write_info('')
     write(message_text,'(5x,a)') ' >> Plane equation : Ax + By + Cz + D = 0'
     call write_info(trim(message_text))
     !write(message_text,'(10x,a,2F10.5)')  'A = ', plane(1), splane(1)
     write(message_text,'(10x,a,F10.5)')  'A = ', plane(1)
     call write_info(trim(message_text))
     !write(message_text,'(10x,a,2F10.5)')  'B = ', plane(2), splane(2)
     write(message_text,'(10x,a,F10.5)')  'B = ', plane(2)
     call write_info(trim(message_text))
     !write(message_text,'(10x,a,2F10.5)')  'C = ', plane(3), splane(3)
     write(message_text,'(10x,a,F10.5)')  'C = ', plane(3)
     call write_info(trim(message_text))
     !write(message_text,'(10x,a,2F10.5)')  'D = ', plane(4), splane(4)
     write(message_text,'(10x,a,F10.5)')  'D = ', plane(4)
     call write_info(trim(message_text))
    end if


    ! new routine in dec. 2016
    !call Get_plane_from_points(3, atom_coord, Plane, SPlane)

    !call write_info('')
    !write(message_text,'(5x,a)') ' >> Plane equation : Ax + By + Cz + D = 0'
    !call write_info(trim(message_text))
    !write(message_text,'(10x,a,2F10.5)')  'A = ', plane(1), splane(1)
    !call write_info(trim(message_text))
    !write(message_text,'(10x,a,2F10.5)')  'B = ', plane(2), splane(2)
    !call write_info(trim(message_text))
    !write(message_text,'(10x,a,2F10.5)')  'C = ', plane(3), splane(3)
    !call write_info(trim(message_text))
    !write(message_text,'(10x,a,2F10.5)')  'D = ', plane(4), splane(4)
    !call write_info(trim(message_text))

   end do


  return
 end subroutine Calcul_plane

!----------------------------------------------------------------
subroutine Calcul_plane_old
  use cryscalc_module, only : atom1_plane,  atom2_plane,  atom3_plane, nb_plane, new_coord, message_text
  use CFML_math_3D,    only : Get_plane_from_points, Get_plane_from_3points
  USE CFML_GlobalDeps, ONLY : cp
  USE io_module
  implicit none
   integer                 :: i
   REAL, DIMENSION(3)      :: atom_1_coord, atom_2_coord, atom_3_coord
   real, dimension(3,3)    :: atom_coord
   CHARACTER(LEN=16)       :: label_1, label_2, label_3
   integer                 :: atom_index
   logical                 :: OK
   real                    :: A, B, C, D


   do i=1, nb_plane
    label_1 = atom1_plane(i)
    label_2 = atom2_plane(i)
    label_3 = atom3_plane(i)

    ok = .false.
    call get_label_atom_coord('plane', i, 1, atom_index, ok)
    IF(.NOT. ok) cycle
    atom_1_coord(1:3) = new_coord(1:3)
    atom_coord(1, 1:3) = new_coord(1:3)

    ok = .false.
    call get_label_atom_coord('plane', i, 2, atom_index, ok)
    IF(.NOT. ok) cycle
    atom_2_coord(1:3) = new_coord(1:3)
    atom_coord(2, 1:3) = new_coord(1:3)

    ok = .false.
    call get_label_atom_coord('plane', i, 3, atom_index, ok)
    IF(.NOT. ok) cycle
    atom_3_coord(1:3) = new_coord(1:3)
    atom_coord(3, 1:3) = new_coord(1:3)


    !write(message_text,'(8x,2a,3F10.5)') trim(label_1), ': ', atom_1_coord
    ! call write_info(trim(message_text))
    !write(message_text,'(8x,2a,3F10.5)') trim(label_2), ': ', atom_2_coord
    ! call write_info(trim(message_text))
    !write(message_text,'(8x,2a,3F10.5)') trim(label_3), ': ', atom_3_coord
    ! call write_info(trim(message_text))

    write(message_text,'(8x,2a,3F10.5)') trim(label_1), ': ', atom_coord(1,1), atom_coord(1,2), atom_coord(1,3)
     call write_info(trim(message_text))
    write(message_text,'(8x,2a,3F10.5)') trim(label_2), ': ', atom_coord(2,1), atom_coord(2,2), atom_coord(2,3)
     call write_info(trim(message_text))
    write(message_text,'(8x,2a,3F10.5)') trim(label_3), ': ', atom_coord(3,1), atom_coord(3,2), atom_coord(3,3)
     call write_info(trim(message_text))


    !call Get_plane_from_points(atom_1_coord, atom_2_coord, atom_3_coord, A, B, C, D)
    call Get_plane_from_3points(atom_1_coord, atom_2_coord, atom_3_coord, A, B, C, D)
    call write_info('')
    write(message_text,'(5x,a)') ' >> Plane equation : Ax + By + Cz + D = 0'
    call write_info(trim(message_text))
    write(message_text,'(10x,a,F10.5)')  'A = ', A
    call write_info(trim(message_text))
    write(message_text,'(10x,a,F10.5)')  'B = ', B
    call write_info(trim(message_text))
    write(message_text,'(10x,a,F10.5)')  'C = ', C
    call write_info(trim(message_text))
    write(message_text,'(10x,a,F10.5)')  'D = ', D
    call write_info(trim(message_text))

    ! new routine in dec. 2016
    !call Get_plane_from_points(3, atom_coord, Plane, SPlane)

    !call write_info('')
    !write(message_text,'(5x,a)') ' >> Plane equation : Ax + By + Cz + D = 0'
    !call write_info(trim(message_text))
    !write(message_text,'(10x,a,2F10.5)')  'A = ', plane(1), splane(1)
    !call write_info(trim(message_text))
    !write(message_text,'(10x,a,2F10.5)')  'B = ', plane(2), splane(2)
    !call write_info(trim(message_text))
    !write(message_text,'(10x,a,2F10.5)')  'C = ', plane(3), splane(3)
    !call write_info(trim(message_text))
    !write(message_text,'(10x,a,2F10.5)')  'D = ', plane(4), splane(4)
    !call write_info(trim(message_text))

   end do


  return
 end subroutine Calcul_plane_old

!----------------------------------------------------------------
subroutine calcul_2theta(stl, angle_2theta)
 use cryscalc_module, only     : pi, wavelength, debug_proc
 implicit none
  real, intent(in)             :: stl
  real, intent(out)            :: angle_2theta

  real                         :: Z

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_2theta")

   Z = stl * wavelength
   if (Z**2 <1) then
    angle_2theta = 2. * ATAN(Z / (SQRT(1. - Z** 2))) * 180. / PI
   ELSE
    angle_2theta = 180.
   END IF

 return
end subroutine calcul_2theta

!-----------------------------------------------------------------
!subroutine calc_THERM_ANISO()   sept. 2011  >>>> therm.F90
!
! USE cryscalc_module,        ONLY : THERM_uij, THERM_bij, THERM_beta, therm_values,  unit_cell, pi, message_text, crystal_cell
! USE CFML_Crystal_Metrics,  only : U_equiv, Convert_U_B, convert_U_betas, convert_B_U, convert_B_Betas, &
!                                   convert_Betas_U, convert_Betas_B
! !USE CFML_Math_General,    ONLY : sp
! !USE  CFML_Constants,       ONLY : sp
! USE CFML_math_3D,          ONLY : matrix_diageigen
! USE CFML_GlobalDeps,       ONLY : sp
! use CFML_Crystal_Metrics,  only : Set_Crystal_Cell  !, Crystal_cell_type
! USE IO_module

! implicit none
!  integer                          :: i
!  REAL(kind=sp), DIMENSION(6)      :: Xij
!  REAL(kind=sp), DIMENSION(3,3)    :: M_U, eigen
!  REAL(kind=sp), DIMENSION(3)      :: rms
!  REAL, DIMENSION(6)               :: new_ADP
!  REAL(kind=sp)                    :: Ueq
!  !type (Crystal_Cell_Type)         :: crystal_cell

!  call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)

! call write_info('')
! call write_info('   >> ANISOTROPIC THERMAL PARAMETERS CONVERSION:')
! call write_info('')

! Xij(1:6) = therm_values(1:6)

! M_U = reshape((/Xij(1),Xij(4),Xij(5), Xij(4),Xij(2),Xij(6), Xij(5),Xij(6),Xij(3) /),(/3,3/))
! call matrix_diageigen(M_U, rms,eigen)
! do i=1, 3
!  if(rms(i) < 0.0) then
!   write(message_text, '(a)') "   -> Matrix U non-positive definite!"
!   call write_info(trim(message_text))
!   call write_info('')
!  endif
! end do

! IF(THERM_Uij) then
!  Ueq = U_equiv(Crystal_cell, Xij)
!  WRITE(message_text,'(a,6F8.5)') '   >> U_ij (A2):   ', Xij
!  call write_info(TRIM(message_text))

!  WRITE(message_text,'(a,F8.5)')  '   >> Ueq  (A2):   ', Ueq
!  call write_info(TRIM(message_text))
!  new_ADP = convert_U_B(Xij)
!  WRITE(message_text,'(a,6F8.5)') '   >> B_ij (A2):   ', new_ADP
!  call write_info(TRIM(message_text))

!  new_ADP = convert_U_betas(Xij, crystal_cell)
!  WRITE(message_text,'(a,6F8.5)') '   >> Beta_ij:     ', new_ADP
!  call write_info(TRIM(message_text))

! ELSEIF(therm_Bij) then
!  WRITE(message_text,'(a,6F8.5)') '   >> B_ij (A2):   ', Xij
!  call write_info(TRIM(message_text))

!  new_ADP = convert_B_U(Xij)
!  WRITE(message_text,'(a,6F8.5)') '   >> U_ij (A2):   ', new_ADP
!  call write_info(TRIM(message_text))
!  Ueq = U_equiv(Crystal_cell, new_ADP)
!  WRITE(message_text,'(a,F8.5)')  '   >> Ueq  (A2):   ', Ueq
!  call write_info(TRIM(message_text))

!  new_ADP = convert_B_betas(Xij, crystal_cell)
!  WRITE(message_text,'(a,6F8.5)') '   >> Beta_ij:     ', new_ADP
!  call write_info(TRIM(message_text))

! ELSEIF(therm_BETA) then
!  WRITE(message_text,'(a,6F8.5)') '   >> Beta_ij:     ', Xij
!  call write_info(TRIM(message_text))

!  new_ADP = convert_Betas_B(Xij, crystal_cell)
!  WRITE(message_text,'(a,6F8.5)') '   >> B_ij (A2):   ', new_ADP
!  call write_info(TRIM(message_text))

!  new_ADP = convert_Betas_U(Xij, crystal_cell)
!  WRITE(message_text,'(a,6F8.5)') '   >> U_ij (A2):   ', new_ADP
!  call write_info(TRIM(message_text))

!  Ueq = U_equiv(Crystal_cell, new_ADP)
!  WRITE(message_text,'(a,F8.5)')  '   >> Ueq  (A2):   ', Ueq
!  call write_info(TRIM(message_text))

! endif

!end subroutine calc_THERM_ANISO

!--------------------------------------------------------------------------
! calcul de l'angle entre 2 vecteurs
subroutine  calcul_vector_angle(input_string)
 USE cryscalc_module, ONLY : DC_ort, RC_ort, nb_da, nb_ra, U1_da, U2_da, U1_ra, U2_ra, gmd, gmr, message_text, &
                             debug_proc
 USE math_module,     ONLY : length, matvec, scalpr
 USE CFML_Math_General,       ONLY : acosd      ! crysfml
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: input_string
  INTEGER                       :: i, nb_ang
  REAL, DIMENSION(3)            :: V1, V2
  REAL                          :: V1_long, V2_long
  REAL, DIMENSION(3,3)          :: ort, gm
  REAL, DIMENSION(3)            :: aux1
  real                          :: aux
  REAL                          :: angle
  CHARACTER (LEN=12)            :: space_string
  logical                       :: input_direct, input_reciprocal

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_vector_angle ("//trim(input_string)//")")

  input_direct     = .false.
  input_reciprocal = .false.


  if(len_trim(input_string) == 6) then
   if(input_string(1:6) == 'direct') input_direct = .true.
  elseif(len_trim(input_string) ==  10) then
   if(input_string(1:10) == 'reciprocal') input_reciprocal = .true.
  endif

  call write_info('')
  IF(input_direct) then
   nb_ang = nb_da
   V1(1:3) = U1_da(1:3, nb_da)
   V2(1:3) = U2_da(1:3, nb_da)
   ort = DC_ort
   gm = gmd
   space_string = 'A'
   call write_info('   >> Calculation of angle between two vectors of the direct space:')
  endif

  IF(input_reciprocal) then
   nb_ang = nb_ra
   V1(1:3) = U1_ra(1:3, nb_ra)
   V2(1:3) = U2_ra(1:3, nb_ra)
   ort = RC_ort
   gm = gmr
   space_string = 'A-1'
   call write_info('   >> Calculation of angle between two vectors of the reciprocal space:')
  endif

  do i=1, nb_ang
    call length(ORT, V1, V1_long)
    call length(ORT, V2, V2_long)
    call matvec(gm, V1, aux1 )
    call scalpr(V2, aux1, aux)
    angle = acosd(aux / V1_long / V2_long)

    call write_info('')
    write(message_text, '(a,3f8.4,a,f8.4,1x,a)') ' -> Vector 1: ', V1,' Length: ',V1_long, TRIM(space_string)
    call write_info(TRIM(message_text))
    write(message_text, '(a,3f8.4,a,f8.4,1x,a)') ' -> Vector 2: ', V2,' Length: ',V2_long, TRIM(space_string)
    call write_info(TRIM(message_text))
    write(message_text, '(a,f8.4,a)')         ' -> Angle between V1 and V2: ',angle,' degrees'
    call write_info(TRIM(message_text))
  end do

 RETURN
end subroutine   calcul_vector_angle

!--------------------------------------------------------------------------------
subroutine calcul_theta()
 USE HKL_module
 USE cryscalc_module, ONLY : unit_cell, wavelength, known_theta, keyword_CELL, keyword_WAVE, debug_proc
 USE IO_module
 implicit none
  INTEGER                        :: i
  REAL                           :: Q, Z
  real, parameter                :: pi=3.1415926535897932
  real, parameter                :: eps = 0.000001

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_theta")


  IF(.NOT. keyword_CELL) then
   call write_info('')
   call write_info('  !! CELL keyword mandatory for the Rint calculation procedure !!')
   call write_info('')
   return
  endif

  IF(.NOT. keyword_WAVE) then
   call write_info('')
   call write_info('  !! WAVE keyword mandatory for the Rint calculation procedure !!')
   call write_info('')
   return
  END if


  if(debug_proc%level_3)  call write_debug_proc_level(3, "CALCUL_Q for all reflections")

  do i=1, n_ref
   call calcul_Q(real(h(i)), real(k(i)), real(l(i)),unit_cell%param, Q)

   sintheta_lambda(i) = SQRT(Q)/2.
   d_hkl(i)           = 1/(2*  sintheta_lambda(i))

   IF (d_hkl(i) > eps) then
    Z =wavelength / (2. * d_hkl(i))
   else
    cycle
   endif

   IF (Z**2 < 1.) then
    theta_hkl(i) = ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
   else
    theta_hkl(i) = 90.
   endif
  end do
  known_theta = .true.


 return
end subroutine calcul_theta

!----------------------------------------------------------
subroutine multi_matrix
 USE cryscalc_module, ONLY : message_text, input_line, debug_proc
 USE IO_module

 implicit none
  INTEGER                     :: i_error, i, j, k
  REAL, DIMENSION(3,3)        :: M1, M2, M

  if(debug_proc%level_2)  call write_debug_proc_level(2, "multi_matrix")

 do
  call write_info('')
  call write_info('  > Enter matrix #1 components (m11 m12 m13 ... m33):  ')
  call read_input_line(input_line)
  READ(input_line,*, IOSTAT=i_error)  M1(1,1), M1(1,2), M1(1,3), M1(2,1), M1(2,2), M1(2,3), M1(3,1), M1(3,2), M1(3,3)
  IF(i_error ==0) exit
  call write_info('')
  call write_info('  ... bad input for matrix #1 !!! ')
 END do

 do
  call write_info('')
  call write_info('  > Enter matrix #2 components (m11 m12 m13 ... m33):  ')
  call read_input_line(input_line)
  READ(input_line,*, IOSTAT=i_error)  M2(1,1), M2(1,2), M2(1,3), M2(2,1), M2(2,2), M2(2,3), M2(3,1), M2(3,2), M2(3,3)
  IF(i_error ==0) exit
  call write_info('')
  call write_info('  ... bad input for matrix #2 !!! ')
 END do

 call write_info('')

 ! la multiplication de 2 matrices 3*3 est une matrice 3*3
! cij = ai1.b1j +  ai2.b2j + ai3.b3j

 do i=1,3
  do j=1,3
   M(i,j)=0.
   do k=1,3
      M(i,j) = M(i,j) + M1(i,k)*M2(k,j)
   end do
  end do
 end do

 call write_info('  . matrix 1:')
 WRITE(message_text,'(10x,3F9.4)') M1(1,1), M1(1,2), M1(1,3)
 call write_info(TRIM(message_text))
 WRITE(message_text,'(10x,3F9.4)') M1(2,1), M1(2,2), M1(2,3)
 call write_info(TRIM(message_text))
 WRITE(message_text,'(10x,3F9.4)') M1(3,1), M1(3,2), M1(3,3)
 call write_info(TRIM(message_text))

 call write_info('')
 call write_info('  . matrix 2:')
 WRITE(message_text,'(10x,3F9.4)') M2(1,1), M2(1,2), M2(1,3)
 call write_info(TRIM(message_text))
 WRITE(message_text,'(10x,3F9.4)') M2(2,1), M2(2,2), M2(2,3)
 call write_info(TRIM(message_text))
 WRITE(message_text,'(10x,3F9.4)') M2(3,1), M2(3,2), M2(3,3)
 call write_info(TRIM(message_text))


 call write_info(' ')
 call write_info('  . matrix 1 * matrix 2:')
 WRITE(message_text,'(10x,3F9.4)') M(1,1), M(1,2), M(1,3)
 call write_info(TRIM(message_text))
 WRITE(message_text,'(10x,3F9.4)') M(2,1), M(2,2), M(2,3)
 call write_info(TRIM(message_text))
 WRITE(message_text,'(10x,3F9.4)') M(3,1), M(3,2), M(3,3)
 call write_info(TRIM(message_text))
 call write_info('')

 RETURN
END subroutine multi_matrix

!---------------------------------------------------------------------------

subroutine verif_linearity(str, linear)
! verification de la linearite de str1-str2-str3 et str2-str3-str4
 use cryscalc_module, only : atom1_ang, atom2_ang, atom3_ang, atom4_ang,  new_coord, debug_proc

 implicit none
  character (len=*), dimension(4), intent(in)    :: str
  logical,                         intent(inout) :: linear
  REAL, DIMENSION(3)                             :: atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord
  real                                           :: angle
  logical                                        :: ok
  integer                                        :: atom_index

  if(debug_proc%level_3)  call write_debug_proc_level(3, "verif_linearity")

  linear = .false.


  ! A1 - A2 - A3
  atom1_ang(1) = str(1)
  atom2_ang(1) = str(2)
  atom3_ang(1) = str(3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 1, atom_index, ok)
  IF(.NOT. ok) return
  atom_1_coord(1:3) = new_coord(1:3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 2, atom_index, ok)
  IF(.NOT. ok) return
  atom_2_coord(1:3) = new_coord(1:3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 3, atom_index, ok)
  IF(.NOT. ok) return
  atom_3_coord(1:3) = new_coord(1:3)

  atom_4_coord(1:3) = -99.

  call angle_calculation(atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord, angle)
  if(angle > 179.) then
   linear = .true.
   return
  endif

  ! A2 - A3 - A4
  atom1_ang(1) = str(2)
  atom2_ang(1) = str(3)
  atom3_ang(1) = str(4)

  ok = .false.
  call get_label_atom_coord('ang', 1, 1, atom_index, ok)
  IF(.NOT. ok) return
  atom_1_coord(1:3) = new_coord(1:3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 2, atom_index, ok)
  IF(.NOT. ok) return
  atom_2_coord(1:3) = new_coord(1:3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 3, atom_index, ok)
  IF(.NOT. ok) return
  atom_3_coord(1:3) = new_coord(1:3)

  call angle_calculation(atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord, angle)
  if(angle > 179.) then
   linear = .true.
   return
  endif



 return
end subroutine verif_linearity

!--------------------------------------------------------------------
subroutine estimate_Zunit
 use Cryscalc_module, only : molecule, unit_cell, SFAC, sto, nb_atoms_type, Z_unit, keyword_Zunit
 implicit none
  real                 :: estimated_Z

  estimated_Z = 1.5 * 0.6023 * unit_cell%volume / molecule%weight
  Z_unit = nint(estimated_Z + 0.5)

  SFAC%number(1:nb_atoms_type) =  sto(1:nb_atoms_type) * Z_unit
  call get_content               ! molecule%content
  keyword_Zunit = .true.

end subroutine estimate_Zunit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      function FETA(X)
      USE pattern_profile_module, ONLY  : eta
      implicit none
       REAL, INTENT(IN)    :: x
       REAL                :: feta

       feta= eta - (1.36603 - 0.47719*X + 0.11116*X*X)*X
       return
      end function
!
!------------------------------------------------------------------
      function Fgau(X)
      USE pattern_profile_module, ONLY : HG, HL, H
      implicit none
       REAL, INTENT(IN)    :: x
       REAL                :: fgau

        fgau = H - (X**5. +2.69269*X**4.*HL + 2.42843*X**3.*HL**2. + 4.47163*X**2.*HL**3. +   &
                   0.07842*X*HL**4. + HL**5.)**0.2
      return
      end function

!------------------------------------------------------------------
      FUNCTION zbrent(func,x1,x2,tol)
       USE IO_module
       USE cryscalc_module, only : lecture_OK

       implicit none
       REAL               :: zbrent
       REAL               :: func
       REAL, INTENT(IN)   :: x1, x2, tol
       INTEGER, PARAMETER :: ITMAX = 100
       REAL, parameter    :: eps = 3.e-08
       integer            :: iter
       REAL               :: a,b,c,d,e, fa,fb, fc, xm,p,q,r,s
       REAL               :: tol1

      a=x1
      b=x2
      fa=func(a)
      fb=func(b)

      if((fa > 0. .and. fb > 0.) .or. (fa < 0. .and. fb < 0.)) then
       !call write_info('')
       !call write_info(' > Wrong values to apply T.C.H. formulae !!')
       !call write_info('')
       lecture_ok = .false.
       !stop     ! ' => root must be bracketed for zbrent'
      else
       lecture_ok = .true.
      endif


      c=b
      fc=fb
      do iter=1,ITMAX
        if((fb > 0. .and. fc > 0.) .or. (fb < 0. .and. fc < 0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        !if(abs(xm) < tol1 .or. fb == 0.)then
        if(abs(xm) < tol1 .or. ABS(fb) < eps)then
          zbrent=b
          return
        endif
        if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa
          !if(a ==c) then
          IF(ABS(a-c) < eps) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p > 0.) q=-q
          p=abs(p)
          if(2.*p < min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) > tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
     END do

      zbrent=b
      return
      END function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
