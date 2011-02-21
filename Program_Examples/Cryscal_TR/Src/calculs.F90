!     Last change:  TR   17 Jul 2007    4:50 pm


!------------------------------------------------------------------------------
!
subroutine get_label_atom_coord(input_string, i,n, ok)
! determination des coord. atomiques a partir du label
! rq: le label peut etre de la forme AA_$n, ou:
!                                     AA: label de l'atome dans la liste des atomes
!                                     $n: numero de l'operateur de symetrie (dans la liste) a appliquer a l'atome
 use cryscal_module, ONLY       : R => symm_op_rot,  T => symm_op_trans, nb_dist_calc, nb_symm_op,              &
                                  nb_atom, atom1_dist, atom2_dist, atom1_ang, atom2_ang, atom3_ang, atom4_ang,  &
                                  atom_label, atom_coord, new_coord, message_text
 USE IO_module
 USE macros_module, only        : u_case

 implicit none
 CHARACTER(LEN=*), INTENT(IN) :: input_string
 INTEGER, INTENT(IN)          :: i,n
 LOGICAL, INTENT(INOUT)       :: ok
 INTEGER                      :: i1, j, atom_index, long, num_sym_op
 CHARACTER(LEN=6)             :: label

 IF(input_string(1:4) == 'dist') then
  IF(n==1) then
   label = atom1_dist(i)
  ELSEif(n==2) then
   label = atom2_dist(i)
  endif
 ELSEIF(input_string(1:3) == 'ang') then
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
 use cryscal_module,          ONLY : nb_ang_calc, atom1_ang, atom2_ang, atom3_ang, atom4_ang, SP_value, new_coord, message_text
 USE CFML_Math_General,       ONLY : acosd
 !USE CFML_Constants,          ONLY : sp
 USE CFML_GlobalDeps,                 ONLY : sp
 USE IO_module

 implicit none
  integer                 :: i
  REAL, DIMENSION(3)      :: atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord
  REAL(KIND=sp)           :: dist_12, dist_34
  REAL(KIND=sp)           :: cos_ang, angle
  CHARACTER(LEN=16)       :: label_1, label_2, label_3, label_4
  LOGICAL                 :: ok


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
   call get_label_atom_coord('ang', i,1, ok)
   IF(.NOT. ok) cycle
   atom_1_coord(1:3) = new_coord(1:3)

   ok = .false.
   call get_label_atom_coord('ang', i,2, ok)
   IF(.NOT. ok) cycle
   atom_2_coord(1:3) = new_coord(1:3)

   ok = .false.
   call get_label_atom_coord('ang', i,3, ok)
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
    call get_label_atom_coord('ang', i,4, ok)
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
 USE cryscal_module,          ONLY : SP_value
 USE CFML_Math_General,       ONLY : acosd
 !USE CFML_Constants,          ONLY : sp
 USE CFML_GlobalDeps,         ONLY : sp

 implicit none
  real, dimension(3), intent(in)    :: coord_1, coord_2, coord_3, coord_4
  real,               intent(inout) :: angle
  real                              :: dist_12, dist_34
  real                              :: cos_angle

  if(coord_4(1) < -98. .and. coord_4(2) < -98. .and. coord_4(3) < -98.) then
    ! distance entre 1 et 2
    call distance_calculation(coord_1(:), coord_2(:), dist_12)

    ! distance entre 2 et 3
    call distance_calculation(coord_2(:), coord_3(:), dist_34)

    ! produit scalaire 21.23
    call scalar_product(coord_2(:), coord_1(:), coord_2(:), coord_3(:))
    cos_angle = SP_value / (dist_12 * dist_34)
    IF((1-ABS(cos_angle))  < 0.0001) cos_angle=SIGN(1.,cos_angle)

    angle = ACOSd(cos_angle)

  else
    ! distance entre 1 et 2
    call distance_calculation(coord_1(:), coord_2(:), dist_12)

    ! distance entre 3 et 4
    call distance_calculation(coord_3(:), coord_4(:), dist_34)

    ! produit scalaire 21.34
    call scalar_product(coord_2(:), coord_1(:), coord_3(:), coord_4(:))
    cos_angle = SP_value / (dist_12 * dist_34)
    IF((1-ABS(cos_angle))  < 0.0001) cos_angle=SIGN(1.,cos_angle)

    angle = ACOSd(cos_angle)

  endif

 return
end subroutine angle_calculation

!-------------------------------------------------------------------------------
subroutine volume_calculation(input_string)
 use cryscal_module,  only        :  ON_SCREEN, pi, unit_cell, keyword_create_CIF, CIF_unit, &
                                     DC_ort, ort_DC, RC_ort, ort_RC, GMD, GMR, keyword_CELL, message_text

 USE CFML_Math_General      ,  ONLY        :  cosd, sind, acosd
 USE math_module,     ONLY        :  orthog, metric, matinv
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: input_string
  REAL                          :: CosAlfa, CosBeta, CosGamma
  REAL                          :: cosa2, cosb2, cosg2
  INTEGER                       :: i
  INTEGER                       :: ifail

  if(ON_SCREEN) then
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

  unit_cell%volume = SQRT(1 + 2 * CosAlfa * cosbeta * cosgamma - cosa2 - cosb2 - cosg2)
  unit_cell%volume = unit_cell%volume * unit_cell%param(1) * unit_cell%param(2) * unit_cell%param(3)

  ! reciprocal parameters
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


  if(ON_SCREEN) then
  IF(input_string=='out') then
   call write_info('            >>> Direct cell parameters:  ')

   WRITE(message_text,'(a,F10.5,a)') '                 .     a = ',unit_cell%param(1),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .     b = ',unit_cell%param(2),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .     c = ',unit_cell%param(3),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .  alfa = ',unit_cell%param(4),' (deg)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .  beta = ',unit_cell%param(5),' (deg)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 . gamma = ',unit_cell%param(6),' (deg)'
    call write_info(TRIM(message_text))

   call write_info('')
   WRITE(message_text,'(a,F15.2,a)') '                 . unit cell volume : ', unit_cell%volume, ' A3'
   call write_info(TRIM(message_text))

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

  if(keyword_create_CIF) then
   call write_CIF_file('UNIT_CELL_INFO')
   call write_CIF_file('CELL_PARAM')
  endif

 return
end subroutine volume_calculation

!------------------------------------------------------------------------------
subroutine get_crystal_SIZE_min_max()
 USE cryscal_module, ONLY : crystal
  implicit none
   INTEGER               :: i, i_mid
   INTEGER, DIMENSION(1) :: i_min, i_max

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
subroutine crystal_volume_calculation()
 USE cryscal_module,    ONLY : ON_SCREEN, crystal, pi, keyword_create_CIF, CIF_unit, keyword_SIZE, message_text
 USE IO_module

 implicit none
   INTEGER                  :: i, i_mid
   INTEGER, DIMENSION(1)    :: i_min, i_max

   if(ON_SCREEN) then
   call write_info(' ')
   call write_info('          Crystal dimensions ')
   call write_info('          ------------------')
   call write_info(' ')
   endif

   crystal%volume =  crystal%size(1)* crystal%size(2)* crystal%size(3)
   crystal%radius =  ((3.*crystal%volume) / (4.*pi) )  ** (1./3.)

   if(ON_SCREEN) then
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

!-------------------------------------------------------------------------------



subroutine calcul_dhkl
 use cryscal_module, only          : pi, unit_cell, wavelength, keyword_WAVE,  keyword_SPGR, &
                                     nb_hkl, H, shift_2theta, keyword_QVEC, Qvec, message_text, SPG
 USE CFML_Math_General,              ONLY : sind
 USE CFML_Reflections_Utilities,     ONLY : HKL_Absent
 USE CFML_crystallographic_symmetry, ONLY : Space_Group_Type, set_spacegroup

 USE IO_module
 implicit none

 integer                              :: i
 real                                 :: Q_hkl, d_hkl, stl_hkl
 real                                 :: a_star, b_star, c_star
 real                                 :: alfa_rad, beta_rad, gama_rad
 real                                 :: cos_alfa_star, cos_beta_star, cos_gama_star
 real                                 :: Z, angle_2theta
 INTEGER                              :: nb_ref, n
 REAL, DIMENSION(3)                   :: HQ
 CHARACTER (len=1)                    :: ind_H
 logical                              :: absent


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


 a_star = unit_cell%param(2)*unit_cell%param(3) * sind(unit_cell%param(4)) / unit_cell%Volume
 b_star = unit_cell%param(3)*unit_cell%param(1) * sind(unit_cell%param(5)) / unit_cell%Volume
 c_star = unit_cell%param(1)*unit_cell%param(2) * sind(unit_cell%param(6)) / unit_cell%Volume

 alfa_rad = unit_cell%param(4) * pi/180.
 beta_rad = unit_cell%param(5) * pi/180.
 gama_rad = unit_cell%param(6) * pi/180.

 cos_alfa_star = (cos(beta_rad) * cos(gama_rad) - cos(alfa_rad) ) / (sin(beta_rad) * sin(gama_rad))
 cos_beta_star = (cos(gama_rad) * cos(alfa_rad) - cos(beta_rad) ) / (sin(gama_rad) * sin(alfa_rad))
 cos_gama_star = (cos(alfa_rad) * cos(beta_rad) - cos(gama_rad) ) / (sin(alfa_rad) * sin(beta_rad))

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

  Q_hkl = HQ(1)**2 * a_star**2 + HQ(2)**2 * b_star**2 + HQ(3)**2 * c_star**2              &
        + 2*HQ(1)*HQ(2) * a_star*b_star*cos_gama_star                                      &
        + 2*HQ(2)*HQ(3) * b_star*c_star*cos_alfa_star                                      &
        + 2*HQ(3)*HQ(1) * c_star*a_star*cos_beta_star

  d_hkl = 1/sqrt(Q_hkl)
  stl_hkl = 1/(2*d_hkl)
  Q_hkl = 4*pi*stl_hkl

  IF(keyword_WAVE) then
   Z = wavelength / (2 * d_hkl)
   IF (Z**2 < 1.) then
    angle_2theta = 2 * ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
   else
    angle_2theta = 180.
   endif
   angle_2theta = angle_2theta + shift_2theta

   if(keyword_QVEC) then
    write(message_text, '(4x,a1,3(1x,F6.2),5x,5(5x,F10.4))') ind_H(:), HQ(1:3), d_hkl, stl_hkl, Q_hkl, angle_2theta/2., angle_2theta
   else
    write(message_text, '(5x,   3(1x,F6.2),5x,5(5x,F10.4))')           HQ(1:3), d_hkl, stl_hkl, Q_hkl, angle_2theta/2., angle_2theta
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

end subroutine calcul_dhkl


!----------------------------------------------------------
subroutine atomic_density_calculation( )
 USE cryscal_module, ONLY : ON_SCREEN, nb_atoms_type, SFAC_number, SFAC_type, nb_at, unit_cell, message_text
 USE IO_module

 implicit none
 !local variables
 REAL                                      :: density
 INTEGER                                   :: i

 if(ON_SCREEN) then
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
  nb_at(i) = SFAC_number(i) / unit_cell%volume * 1.E24

  if(ON_SCREEN) then
  write(message_text, '(a8,10x,F9.2)') TRIM(SFAC_type(i)),  nb_at(i)/1.e22
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
 USE cryscal_module, only           : nb_atoms_type, SFAC_type, num_atom, known_atomic_label, known_atomic_features
 USE macros_module, ONLY            : u_case
 USE IO_module

 implicit none
  INTEGER                          :: i, j

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
  SFAC_type(i) = ADJUSTL(SFAC_type(i))
  !
  j=0
  do
    j=j+1
    if (j>201) then
     call write_info('')
     call write_info('   !!! '//TRIM(SFAC_type(i))// ': incorrect atomic symbol !!')
     return
    end if

    if (u_case(SFAC_type(i)(1:))==atom(j)%symbol(1:)) exit

  end do
  Num_atom(i) = j

 end do


 return
end subroutine atomic_identification


!------------------------------------------------------------------------------

subroutine molecular_weight()
 USE atomic_data
 USE cryscal_module, only           : ON_SCREEN, nb_atoms_type, SFAC_number, SFAC_type, sto, num_atom, Z_unit,       &
                                      keyword_CHEM, neutrons, keyword_create_CIF, CIF_unit, message_text, &
                                      molecule
 USE IO_module,      only           : write_info
 USE macros_module,  only           : l_case

 implicit none
  INTEGER                                     :: i, long
  REAL,              DIMENSION(nb_atoms_type) :: Weight_percent, atomic_percent
  CHARACTER(LEN=16), DIMENSION(nb_atoms_type) :: labl
  REAL                                        :: mol_weight, Total_sto
  CHARACTER (LEN=16)                          :: fmt_

 if(ON_SCREEN) then
 call write_info(' ')
 call write_info('          Molecular features ')
 call write_info('          ------------------')
 endif

 do i=1, nb_atoms_type

  ! molecular formula
  IF(LEN_TRIM(SFAC_type(i))==2) SFAC_type(i)(2:2) = l_case(SFAC_type(i)(2:2))

  IF(.NOT. keyword_CHEM)  then
   sto(i) = SFAC_number(i) / Z_unit
  !else
  ! sto(i) = SFAC_number(i)
  endif

  IF(INT(sto(i)) < 10) then
   if (ABS(INT(sto(i))-sto(i)) < 0.01) then
    WRITE(labl(i), '(a,I1)') trim(SFAC_type(i)), INT(sto(i))
   else
    WRITE(labl(i), '(a,F5.2)') trim(SFAC_type(i)), sto(i)
   endif
  ELSEIF(INT(sto(i)) < 100) THEN
   if (ABS(INT(sto(i))-sto(i)) < 0.01) then
    WRITE(labl(i), '(a,I2)') TRIM(SFAC_type(i)), INT(sto(i))
   else
    WRITE(labl(i), '(a,F6.2)') trim(SFAC_type(i)), sto(i)
   endif
  ELSEIF(INT(sto(i)) < 1000) THEN
   if (ABS(INT(sto(i))-sto(i)) < 0.01) then
    WRITE(labl(i), '(a,I3)') TRIM(SFAC_type(i)), INT(sto(i))
   else
    WRITE(labl(i), '(a,F7.2)') trim(SFAC_type(i)), sto(i)
   endif
  endif

 end do


 ! molecular weight
 Mol_Weight = 0.
 Total_sto  = 0.
 do  i=1, nb_atoms_type
  Mol_Weight = Mol_Weight + atom(Num_atom(i))%weight * sto(i)
  total_sto = total_sto + sto(i)
 end do
 molecule%weight = mol_weight

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
 write(molecule%formula, fmt_) (TRIM(labl(i)),i=1,nb_atoms_type)

 if(ON_SCREEN) then
 call write_info('')
 WRITE(message_text, '(2a)')       '   >> Molecular formula: ',trim(molecule%formula)
 call write_info(TRIM(message_text))

 call write_info('')
 WRITE(message_text,'(a,F8.2)') '   >> Molecular weight: ',Molecule%weight
  call write_info(TRIM(message_text))
 write(message_text,'(a)')      ' '
  call write_info(TRIM(message_text))
 write(message_text,'(a)')      '    Atom     Atomic weight Stoechiometry        % atomic        % weight'
  call write_info(TRIM(message_text))
 write(message_text,'(a)')      ' '
  call write_info(TRIM(message_text))
 do i=1, nb_atoms_type
  write(message_text, '(a8,10x,F8.4,8x,F6.2, 2(10x,F6.2))') TRIM(SFAC_type(i)),  atom(Num_atom(i))%weight, &
                                                            sto(i), Atomic_percent(i), Weight_percent(i)
  call write_info(TRIM(message_text))
 end do
 endif

 IF(keyword_create_CIF)  call write_CIF_file('CHEMICAL_INFORMATION')

 ! calcul du F000
 !if (.not. neutrons) call F000_calculation('X')


 return
end subroutine molecular_weight
!------------------------------------------------------------------------------

subroutine density_calculation()
 USE cryscal_module, ONLY         : ON_SCREEN, molecule, Z_unit, unit_cell, keyword_create_CIF, CIF_unit, message_text
 USE IO_module

 implicit none
  REAL, parameter                 :: avogadro = 6.022
  REAL                            :: density

  molecule%density = (Z_unit * Molecule%weight) / (0.1* avogadro * unit_cell%volume  )

  if(ON_SCREEN) then
  call write_info('')
  WRITE(message_text,'(5x,a,F9.3)') "   >> Density (g/cm3)    = ", molecule%density
  call write_info(TRIM(message_text))
  call write_info('')
  end if

  IF(keyword_create_CIF)  call write_CIF_file('CRYSTAL_DENSITY')


 RETURN
end subroutine density_calculation
!------------------------------------------------------------------------------

subroutine calcul_barycentre
 use cryscal_module
 USE IO_module

  implicit none
   integer                      :: i , j , f
   real, dimension(3)           :: bary_coord
   integer, dimension(500)      :: bary_atom    ! numero des atomes dans le calcul du barycentre


 do i = 1, nb_bary_calc
  ! quels atomes
  do j= 1, nb_atom
    do f=1,nb_atom_bary(i)
     if (atom_bary(i,f) == atom_label(j)) bary_atom(f) = j
    end do
  end do

  bary_coord(1) = sum(atom_coord(1, bary_atom(1:nb_atom_bary(i)))) /   nb_atom_bary(i)
  bary_coord(2) = sum(atom_coord(2, bary_atom(1:nb_atom_bary(i)))) /   nb_atom_bary(i)
  bary_coord(3) = sum(atom_coord(3, bary_atom(1:nb_atom_bary(i)))) /   nb_atom_bary(i)


  call write_info('')
  write(message_text,'(a,I3)')      '  ==> Centroid #: ', i
   call write_info(TRIM(message_text))
  call write_info('')
  do j=1, nb_atom_bary(i)
   write(message_text,'(5x,2a,3F10.5)')       atom_label(j)(1:4), ': '  ,  atom_coord(1,bary_atom(j)) , &
                                                                           atom_coord(2,bary_atom(j)) , &
                                                                           atom_coord(3,bary_atom(j))
   call write_info(TRIM(message_text))
  end do
  call write_info('')
  write(message_text,'(5x,a,3F10.5)')   'centroid coordinates: ', bary_coord(1:3)
   call write_info(TRIM(message_text))


 end do

 return
end subroutine calcul_barycentre
!----------------------------------------------------------------

subroutine calcul_2theta(stl, angle_2theta)
 use cryscal_module, only       : pi, wavelength
 implicit none
  real, intent(in)             :: stl
  real, intent(out)            :: angle_2theta

  real                         :: Z

   Z = stl * wavelength
   if (Z**2 <1) then
    angle_2theta = 2. * ATAN(Z / (SQRT(1. - Z** 2))) * 180. / PI
   ELSE
    angle_2theta = 180.
   END IF

 return
end subroutine calcul_2theta

!-----------------------------------------------------------------

subroutine calc_therm_iso()
 USE cryscal_module, ONLY : THERM_Biso, THERM_Uiso, nb_therm_values, therm_values, pi, message_text
 USE IO_module

 implicit none
  integer                  :: i
  real                     :: new_therm_value

 call write_info('')
 call write_info('   >> ISOTROPIC THERMAL PARAMETERS CONVERSION:')
 call write_info('')

 if (THERM_Biso) then      ! Biso --> Uiso
  call write_info('       Biso       Uiso')

  do i=1, nb_therm_values
   new_therm_value = therm_values(i) / (8*pi**2)
   WRITE(message_text, '(2(5x,F8.4))') therm_values(i), new_therm_value
   call write_info(TRIM(message_text))
  end do
  return

 ELSEIF(THERM_Uiso) then   ! Uiso --> Biso
  call write_info('       Uiso       Biso')

  do i=1, nb_therm_values
   new_therm_value = therm_values(i) * (8*pi**2)
   WRITE(message_text,'(2(5x,F8.4))') therm_values(i), new_therm_value
   call write_info(TRIM(message_text))
  end do
  return
 end if



 RETURN
end subroutine calc_therm_iso

!--------------------------------------------------------------------------

subroutine calc_THERM_ANISO()
 USE cryscal_module,        ONLY : THERM_uij, THERM_bij, THERM_beta, therm_values,  unit_cell, pi, message_text, crystal_cell
 USE CFML_Crystal_Metrics,  only : U_equiv, Convert_U_B, convert_U_betas, convert_B_U, convert_B_Betas, &
                                   convert_Betas_U, convert_Betas_B
 !USE CFML_Math_General,    ONLY : sp
 !USE  CFML_Constants,       ONLY : sp
 USE CFML_math_3D,          ONLY : matrix_diageigen
 USE CFML_GlobalDeps,       ONLY : sp
 use CFML_Crystal_Metrics,  only : Set_Crystal_Cell  !, Crystal_cell_type
 USE IO_module

 implicit none
  integer                          :: i
  REAL(kind=sp), DIMENSION(6)      :: Xij
  REAL(kind=sp), DIMENSION(3,3)    :: M_U, eigen
  REAL(kind=sp), DIMENSION(3)      :: rms
  REAL, DIMENSION(6)               :: new_ADP
  REAL(kind=sp)                    :: Ueq
  !type (Crystal_Cell_Type)         :: crystal_cell

  call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)

 call write_info('')
 call write_info('   >> ANISOTROPIC THERMAL PARAMETERS CONVERSION:')
 call write_info('')

 Xij(1:6) = therm_values(1:6)

 M_U = reshape((/Xij(1),Xij(4),Xij(5), Xij(4),Xij(2),Xij(6), Xij(5),Xij(6),Xij(3) /),(/3,3/))
 call matrix_diageigen(M_U, rms,eigen)
 do i=1, 3
  if(rms(i) < 0.0) then
   write(message_text, '(a)') "   -> Matrix U non-positive definite!"
   call write_info(trim(message_text))
   call write_info('')
  endif
 end do

 IF(THERM_Uij) then
  Ueq = U_equiv(Crystal_cell, Xij)
  WRITE(message_text,'(a,6F8.5)') '   >> U_ij (A2):   ', Xij
  call write_info(TRIM(message_text))

  WRITE(message_text,'(a,F8.5)')  '   >> Ueq  (A2):   ', Ueq
  call write_info(TRIM(message_text))
  new_ADP = convert_U_B(Xij)
  WRITE(message_text,'(a,6F8.5)') '   >> B_ij (A2):   ', new_ADP
  call write_info(TRIM(message_text))

  new_ADP = convert_U_betas(Xij, crystal_cell)
  WRITE(message_text,'(a,6F8.5)') '   >> Beta_ij:     ', new_ADP
  call write_info(TRIM(message_text))

 ELSEIF(therm_Bij) then
  WRITE(message_text,'(a,6F8.5)') '   >> B_ij (A2):   ', Xij
  call write_info(TRIM(message_text))

  new_ADP = convert_B_U(Xij)
  WRITE(message_text,'(a,6F8.5)') '   >> U_ij (A2):   ', new_ADP
  call write_info(TRIM(message_text))
  Ueq = U_equiv(Crystal_cell, new_ADP)
  WRITE(message_text,'(a,F8.5)')  '   >> Ueq  (A2):   ', Ueq
  call write_info(TRIM(message_text))

  new_ADP = convert_B_betas(Xij, crystal_cell)
  WRITE(message_text,'(a,6F8.5)') '   >> Beta_ij:     ', new_ADP
  call write_info(TRIM(message_text))

 ELSEIF(therm_BETA) then
  WRITE(message_text,'(a,6F8.5)') '   >> Beta_ij:     ', Xij
  call write_info(TRIM(message_text))

  new_ADP = convert_Betas_B(Xij, crystal_cell)
  WRITE(message_text,'(a,6F8.5)') '   >> B_ij (A2):   ', new_ADP
  call write_info(TRIM(message_text))

  new_ADP = convert_Betas_U(Xij, crystal_cell)
  WRITE(message_text,'(a,6F8.5)') '   >> U_ij (A2):   ', new_ADP
  call write_info(TRIM(message_text))

  Ueq = U_equiv(Crystal_cell, new_ADP)
  WRITE(message_text,'(a,F8.5)')  '   >> Ueq  (A2):   ', Ueq
  call write_info(TRIM(message_text))

 endif

end subroutine calc_THERM_ANISO

!--------------------------------------------------------------------------
! calcul de l'angle entre 2 vecteurs
subroutine  calcul_vector_angle(input_string)
 USE cryscal_module, ONLY : DC_ort, RC_ort, nb_da, nb_ra, U1_da, U2_da, U1_ra, U2_ra, gmd, gmr, message_text
 USE math_module,    ONLY : length, matvec, scalpr
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

  call write_info('')
  IF(input_string(1:6) == 'direct') then
   nb_ang = nb_da
   V1(1:3) = U1_da(1:3, nb_da)
   V2(1:3) = U2_da(1:3, nb_da)
   ort = DC_ort
   gm = gmd
   space_string = 'A'
   call write_info('   >> Calculation of angle between two vectors of the direct space:')

  ELSEIF(input_string(1:10) == 'reciprocal') then
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
    write(message_text, '(a,f8.4,a)')         ' -> Angle between them: ',angle,' degrees'
    call write_info(TRIM(message_text))
  end do

 RETURN
end subroutine   calcul_vector_angle

!--------------------------------------------------------------------------------

subroutine calcul_theta()
 USE HKL_module
 USE cryscal_module, ONLY : unit_cell, wavelength, known_theta, keyword_CELL, keyword_WAVE
 USE IO_module
 implicit none
  INTEGER                        :: i
  REAL                           :: Q, Z
  real, parameter                :: pi=3.1415926535897932
  real, parameter                :: eps = 0.000001

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
 USE cryscal_module, ONLY : message_text, input_line
 USE IO_module

 implicit none
  INTEGER                     :: i_error, i, j, k
  REAL, DIMENSION(3,3)        :: M1, M2, M

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
 use cryscal_module, only : atom1_ang, atom2_ang, atom3_ang, atom4_ang,  new_coord

 implicit none
  character (len=*), dimension(4), intent(in)    :: str
  logical,                         intent(inout) :: linear
  REAL, DIMENSION(3)                             :: atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord
  real                                           :: angle
  logical                                        :: ok

  linear = .false.


  ! A1 - A2 - A3
  atom1_ang(1) = str(1)
  atom2_ang(1) = str(2)
  atom3_ang(1) = str(3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 1, ok)
  IF(.NOT. ok) return
  atom_1_coord(1:3) = new_coord(1:3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 2, ok)
  IF(.NOT. ok) return
  atom_2_coord(1:3) = new_coord(1:3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 3, ok)
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
  call get_label_atom_coord('ang', 1, 1, ok)
  IF(.NOT. ok) return
  atom_1_coord(1:3) = new_coord(1:3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 2, ok)
  IF(.NOT. ok) return
  atom_2_coord(1:3) = new_coord(1:3)

  ok = .false.
  call get_label_atom_coord('ang', 1, 3, ok)
  IF(.NOT. ok) return
  atom_3_coord(1:3) = new_coord(1:3)

  call angle_calculation(atom_1_coord, atom_2_coord, atom_3_coord, atom_4_coord, angle)
  if(angle > 179.) then
   linear = .true.
   return
  endif



 return
end subroutine verif_linearity
