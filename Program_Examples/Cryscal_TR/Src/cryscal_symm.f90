!     Last change:  TR   27 Oct 2006    3:19 pm

subroutine decode_sym_op
 USE cryscal_module, ONLY : symm_op_mat, write_symm_op, nb_atom, keyword_SPGR, write_site_info, write_APPLY_symm
 implicit none

   IF(symm_op_mat) THEN    ! forme matricielle
    call get_op_xyz()
   ELSE                    ! forme x,y,z
    !call write_symm_op_xyz( )
    call get_op_mat()
   endif

   IF(WRITE_symm_op) call write_symm_op_mat()
   IF(nb_atom /=0) then
    if(keyword_SPGR) then
     IF(WRITE_SITE_info ) call get_SITE_info()
    else
     IF(WRITE_APPLY_symm) call apply_symm()
    endif
   endif

 RETURN
end subroutine decode_sym_op
!-------------------------------------------------------------------------------

subroutine write_symm_op_xyz()
 USE cryscal_module, ONLY : symm_op_string, nb_symm_op, message_text
 USE IO_module,      ONLY : write_info
 implicit none
  integer       :: i, num

 do num = 1, nb_symm_op
  call write_info('')
  WRITE(message_text,'(a,i3,4a)') '   SYMM #',num,': ', (symm_op_string(i, num), i=1,3)
  call write_info(TRIM(message_text))
  call write_info('')
 end do

end subroutine write_symm_op_xyz



!--------------------------------------------------------------
! determination de la matrice 3*3 et de la translation associées à un operateur de symetrie
subroutine get_op_mat
 USE cryscal_module, ONLY : symm_op_string , R => symm_op_rot,  T => symm_op_trans, nb_symm_op
 implicit none
  INTEGER                             :: num
  CHARACTER (LEN=12), DIMENSION(3)    :: op_STRING
  REAL,               DIMENSION(3,3)  :: op_ROT
  REAL,               DIMENSION(3)    :: op_TRANS



 do num= 1, nb_symm_op
  op_STRING(:) = symm_op_string(:,num)

  call get_op_ROT_TRANS(op_STRING, op_ROT, op_TRANS)

  R(:,:,num)= op_ROT(:,:)
  T(:,num)  = op_TRANS(:)


 end do


END subroutine get_op_mat

!-----------------------------------------------------------------------------------------
subroutine get_op_ROT_TRANS(op_STRING, op_ROT, op_TRANS)
 USE macros_module, ONLY    : remove_car, l_case
 implicit none
  CHARACTER (LEN=*), DIMENSION(3),    INTENT(INOUT) :: op_STRING
  REAL,              DIMENSION(3,3),  INTENT(INOUT) :: op_ROT
  REAL,              DIMENSION(3),    INTENT(INOUT) :: op_TRANS

  INTEGER              :: i


  do i=1,3
   ! remove blank in string
   !call remove_car(op_STRING(i), ' ')
   op_string(i) = remove_car(op_string(i), ' ')
   op_STRING(i) = l_case(op_STRING(i))

   ! partie rotationnelle
   call get_op_ROT(op_STRING(i), i,'x', op_ROT)
   call get_op_ROT(op_STRING(i), i,'y', op_ROT)
   call get_op_ROT(op_STRING(i), i,'z', op_ROT)

   ! partie translationnelle
   call get_op_TRANS(op_STRING(i), i, op_TRANS)
  end do

 RETURN
end subroutine get_op_ROT_TRANS

!-----------------------------------------------------------------

subroutine get_op_xyz()
 USE cryscal_module, ONLY : nb_symm_op, symm_op_ROT, symm_op_TRANS, symm_op_STRING
 implicit none
  INTEGER                         :: op_num
  REAL, DIMENSION(3,3)            :: op_ROT
  REAL, DIMENSION(3)              :: op_TRANS
  CHARACTER(LEN=12), DIMENSION(3) :: op_string

  do op_num = 1, nb_symm_op

   op_ROT(:,:)  = symm_op_ROT(:,:, op_num)
   op_TRANS(:)  = symm_op_TRANS(:, op_num)

   call get_op_string(op_ROT, op_TRANS, op_STRING)

   symm_op_STRING(:,op_num) = op_STRING(:)

  end do

 return
end subroutine get_op_xyz
!-----------------------------------------------------------------
subroutine get_op_STRING(op_ROT, op_TRANS, op_STRING)
! creation operateur chaine de caractere a partir de la matrice de rotaion et de la translation
 USE MACROS_module,  ONLY : remove_car,  l_case

 implicit none
  REAL,              DIMENSION(3,3), INTENT(INOUT) :: op_ROT
  REAL,              DIMENSION(3),   INTENT(INOUT) :: op_TRANS
  CHARACTER(LEN=12), DIMENSION(3),   INTENT(INOUT) :: op_STRING
  ! local variables
  INTEGER                            :: i, k
  REAL, parameter                    :: eps = 0.001
  CHARACTER(LEN=2), DIMENSION(6)     :: string_1
  CHARACTER(LEN=3), DIMENSION(6)     :: string_2
  CHARACTER(LEN=32)                  :: trans_string

  character(len=6), dimension(7)    :: ratio_string_pos, ratio_string_neg
  real            , dimension(7)    :: ratio_real_pos,   ratio_real_neg
  LOGICAL                           :: ratio

  CHARACTER (LEN=4)                 :: string

  ratio = .false.

  string_1(1:6) = (/'+x',  '+y',  '+z',     &
                    '-x',  '-y',  '-z'   /)
  string_2(1:6) = (/'+2x', '+2y', '+2z',    &
                    '-2x', '-2y', '-2z'   /)

  ratio_string_pos(1:7)  = (/'+1/2',             &
                              '+1/3', '+2/3',    &
                              '+1/4', '+3/4',    &
                              '+1/6', '+5/6'/)

  ratio_string_neg(1:7) = (/'-1/2',            &
                            '-1/3',  '-2/3',   &
                            '-1/4',  '-3/4',   &
                            '-1/6',  '-5/6' /)

  ratio_real_pos(1:7)  = (/0.5,                    &
                           0.3333333, 0.6666667,   &
                           0.25,      0.75,        &
                           0.1666667, 0.8333333/)

  ratio_real_neg(1:7) = -ratio_real_pos(1:7)



  op_STRING = ''

  do k=1, 3
   trans_string        = ''

   do i=1,3
    WRITE(string, '(F4.1)') op_ROT(k,i)
    string = ADJUSTL(string)
    IF(string(1:2) == '1.') then
     op_STRING(k) = TRIM(op_STRING(k))//string_1(i)
    elseIF(string(1:3) == '-1.') then
     op_STRING(k) = TRIM(op_STRING(k))//string_1(i+3)
    elseIF(string(1:2) == '2.') then
     op_STRING(k) = TRIM(op_STRING(k))//string_2(i)
    elseIF(string(1:3) == '-2.') then
     op_STRING(k) = TRIM(op_STRING(k))//string_2(i+3)
    endif

    !IF(op_ROT(k,i) == 1.) then
    ! op_STRING(k) = TRIM(op_STRING(k))//string_1(i)
    !elseIF(op_ROT(k,i) == -1.) then
    ! op_STRING(k) = TRIM(op_STRING(k))//string_1(i+3)
    !elseIF(op_ROT(k,i) == 2.) then
    ! op_STRING(k) = TRIM(op_STRING(k))//string_2(i)
    !elseIF(op_ROT(k,i) == -2.) then
    ! op_STRING(k) = TRIM(op_STRING(k))//string_2(i+3)
    !endif
   end do   ! end loop on i

   IF(ABS(op_TRANS(k))> eps) then
    IF(op_TRANS(k) > eps) then
     DO i=1,7
      IF(ABS(op_TRANS(k)-ratio_real_pos(i)) < eps) then
       trans_string = ratio_string_pos(i)
       ratio = .true.
      endif
     END DO
     IF(.not. ratio) WRITE(trans_string, '(a,F6.3)') "+",op_TRANS(k)

    else
     DO i=1,7
      IF(ABS(op_TRANS(k)-ratio_real_neg(i)) < eps) then
       trans_string = ratio_string_neg(i)
       ratio = .true.
      endif
     END DO
     IF(.not. ratio) WRITE(trans_string, '(F6.3)') op_TRANS(k)
    endif

   endif

   IF(LEN_TRIM(trans_string)/=0) then
    op_STRING(k) =  TRIM(op_STRING(k))//ADJUSTL(TRIM(trans_string))
   endif

  END do ! end loop on k






END subroutine get_op_STRING

!------------------------------------------------------------

!subroutine get_op_mat_ROT(ligne, input_string, num)
! USE cryscal_module, ONLY : symm_op_string, R => symm_op_rot,  T => symm_op_trans
!
! implicit none
! INTEGER, INTENT(IN)               :: ligne
! CHARACTER(LEN=1), INTENT(IN)      :: input_string
! integer, intent(in)               :: num
! integer                           :: col
! CHARACTER(LEN=1)                  :: string
! INTEGER                           :: k
! REAL                              :: sign_R
!
!
! k= INDEX(symm_op_string(ligne, num),input_string)
! IF    (input_string =='x') then
!  col = 1
! ELSEIF(input_string =='y') then
!  col = 2
! ELSEIF(input_string =='z') then
!  col = 3
! endif
!
!
! IF(k==0) then
!  R(ligne, col, num) = 0.
!  sign_R            = 1.
! else
!
!  IF(k==1) then
!   R(ligne, col, num) = 1
!   sign_R             = 1.
!
!  ELSE
!
!   READ(symm_op_string(ligne, num)(k-1:k-1),*) string
!   IF(string == "+") then
!    sign_R             = 1.
!    R(ligne, col, num) = 1.
!   ELSEIF(string == '-') then
!    sign_R             = -1.
!    R(ligne, col, num) = 1
!   else
!    READ(string, *) R(ligne, col, num)
!    READ(symm_op_string(ligne, num)(k-2:k-2),*) string
!    IF(string == "+") then
!     sign_R        = 1.
!    ELSEIF(string == '-') then
!     sign_R        = -1.
!    endif
!   endif
!
!
!  endif
! endif
!
! R(ligne, col, num) = sign_R * R(ligne,col, num)
!
!
!
! return
!end subroutine  get_op_mat_ROT

!-------------------------------------------------------------------------------------
! determination de la partie rotationnelle de l'operateur de symmetrie sous la forme matricielle
subroutine get_op_ROT(op_string, ligne, input_string, rot)
! USE cryscal_module, ONLY : symm_op_string, R => symm_op_rot,  T => symm_op_trans

 implicit none
 CHARACTER (LEN=*),    INTENT(IN)    :: op_string
 INTEGER,              INTENT(IN)    :: ligne
 CHARACTER(LEN=1),     INTENT(IN)    :: input_string
 REAL, DIMENSION(3,3), INTENT(INOUT) :: rot
 !integer, intent(in)               :: num
 integer                           :: col
 CHARACTER(LEN=1)                  :: string
 INTEGER                           :: k
 REAL                              :: sign_R


 k= INDEX(op_string,input_string)
 IF    (input_string =='x') then
  col = 1
 ELSEIF(input_string =='y') then
  col = 2
 ELSEIF(input_string =='z') then
  col = 3
 endif


 IF(k==0) then
  Rot(ligne, col) = 0.
  sign_R          = 1.
 else

  IF(k==1) then
   Rot(ligne, col ) = 1
   sign_R           = 1.

  ELSE

   READ(op_string(k-1:k-1),*) string
   IF(string == "+") then
    sign_R          = 1.
    Rot(ligne, col) = 1.
   ELSEIF(string == '-') then
    sign_R          = -1.
    Rot(ligne, col) = 1
   else
    READ(string, *) Rot(ligne, col)
    READ(op_string(k-2:k-2),*) string
    IF(string == "+") then
     sign_R        = 1.
    ELSEIF(string == '-') then
     sign_R        = -1.
    endif
   endif


  endif
 endif

 Rot(ligne, col) = sign_R * Rot(ligne,col)



 return
end subroutine  get_op_ROT

!--------------------------------------------------------------
!
!subroutine get_op_mat_TRANS(ligne, num)
! USE cryscal_module, ONLY : symm_op_string, T => symm_op_trans
!
! implicit none
! INTEGER, INTENT(IN)               :: ligne
! integer, intent(in)               :: num
! CHARACTER(LEN=12)                 :: new_string
!
! INTEGER                           :: i, k, long, long_ratio
! character(len=6), dimension(8)    :: ratio_string_pos, ratio_string_plus, ratio_string_neg
! real            , dimension(8)    :: ratio_real_pos,   ratio_real_neg
! logical                           :: ratio
!
! ratio = .false.
!
! ratio_string_pos(1:8)  = (/'1/2', '3/2',    &
!                            '1/3', '2/3',    &
!                            '1/4', '3/4',    &
!                            '1/6', '5/6'/)
!
! ratio_string_plus(1:8)  = (/'+1/2',  '+3/2',     &
!                             '+1/3',  '+2/3',     &
!                             '+1/4',  '+3/4',     &
!                             '+1/6',  '+5/6'/)
!
! ratio_string_neg(1:8)  = (/'-1/2',  '-3/2',     &
!                            '-1/3',  '-2/3',     &
!                            '-1/4',  '-3/4',     &
!                            '-1/6',  '-5/6'/)
!
! ratio_real_pos(1:8)  = (/0.5,       1.5,        &
!                          0.3333333, 0.6666667,  &
!                          0.25,      0.75,       &
!                          0.1666667, 0.8333333/)
!
! ratio_real_neg(1:8) = -ratio_real_pos(1:8)
!
! new_string = symm_op_string(ligne, num)
! k= INDEX(symm_op_string(ligne, num),'x')
! call get_new_string(new_string, k)
!
!
! k= INDEX(new_string,'y')
! call get_new_string(new_string, k)
!
!
! k= INDEX(new_string,'z')
! call get_new_string(new_string, k)
!
!
! IF(LEN_TRIM(new_string) /=0) then
!  ! permet de reconnaitre les fractions 1/2, 1/3 ... ----------------!
!  long = len_trim(new_string)                                        !
!  do i=1,8                                                           !
!   long_ratio = len_trim(ratio_string_pos(i))                        !
!   if(new_string(1:long)==ratio_string_pos(i)(1:long_ratio)) then    !
!    T(ligne, num) = ratio_real_pos(i)                                !
!    ratio = .true.                                                   !
!    exit                                                             !
!   endif                                                             !
!  end do                                                             !
!                                                                     !
!  do i=1,8                                                           !
!   long_ratio = len_trim(ratio_string_plus(i))                       !
!   if(new_string(1:long)==ratio_string_plus(i)(1:long_ratio)) then   !
!    T(ligne, num) = ratio_real_pos(i)                                !
!    ratio = .true.                                                   !
!    exit                                                             !
!   endif                                                             !
!  end do                                                             !
!                                                                     !
!  if(.not. ratio) then                                               !
!   do i=1,8                                                         !
!    long_ratio = len_trim(ratio_string_neg(i))                       !
!    if(new_string(1:long)==ratio_string_neg(i)(1:long_ratio)) then   !
!     T(ligne, num) = ratio_real_neg(i)                               !
!     ratio = .true.                                                  !
!     exit                                                            !
!    endif                                                            !
!   end do                                                            !
!  endif                                                              !
! !-------------------------------------------------------------------!
!
!  if(.not. ratio) READ(new_string,*) T(ligne, num)
! else
!  T(ligne, num) = 0.
! endif
!
!
! return
!end subroutine get_op_mat_TRANS
!
!--------------------------------------------------------------
! determination de la partie translationnelle de l'operateur de symmetrie sous forme vectorielle
subroutine get_op_TRANS(op_string, ligne, op_trans)

 implicit none
 CHARACTER (LEN=*),    INTENT(IN)    :: op_string
 INTEGER,              INTENT(IN)    :: ligne
 REAL, DIMENSION(3),   INTENT(INOUT) :: op_trans

 CHARACTER(LEN=12)                   :: new_string
 INTEGER                             :: i, k, long, long_ratio
 character(len=6), dimension(8)      :: ratio_string_pos, ratio_string_plus, ratio_string_neg
 real            , dimension(8)      :: ratio_real_pos,   ratio_real_neg
 logical                             :: ratio

 ratio = .false.
 ratio_string_pos(1:8)  = (/'1/2', '3/2',         &
                            '1/3', '2/3',         &
                            '1/4', '3/4',         &
                            '1/6', '5/6'/)

 ratio_string_plus(1:8)  = (/'+1/2',  '+3/2',     &
                             '+1/3',  '+2/3',     &
                             '+1/4',  '+3/4',     &
                             '+1/6',  '+5/6'/)

 ratio_string_neg(1:8)  = (/'-1/2',  '-3/2',      &
                            '-1/3',  '-2/3',      &
                            '-1/4',  '-3/4',      &
                            '-1/6',  '-5/6'/)

 ratio_real_pos(1:8)  = (/0.5,       1.5,         &
                          0.3333333, 0.6666667,   &
                          0.25,      0.75,        &
                          0.1666667, 0.8333333/)

 ratio_real_neg(1:8) = -ratio_real_pos(1:8)

 new_string = op_string
 k= INDEX(op_string,'x')
 call get_new_string(new_string, k)

 k= INDEX(new_string,'y')
 call get_new_string(new_string, k)

 k= INDEX(new_string,'z')
 call get_new_string(new_string, k)


 IF(LEN_TRIM(new_string) /=0) then
  ! permet de reconnaitre les fractions 1/2, 1/3 ... ----------------!
  long = len_trim(new_string)                                        !
  do i=1,8                                                           !
   long_ratio = len_trim(ratio_string_pos(i))                        !
   if(new_string(1:long)==ratio_string_pos(i)(1:long_ratio)) then    !
    op_trans(ligne) = ratio_real_pos(i)                                !
    ratio = .true.                                                   !
    exit                                                             !
   endif                                                             !
  end do                                                             !
                                                                     !
  do i=1,8                                                           !
   long_ratio = len_trim(ratio_string_plus(i))                       !
   if(new_string(1:long)==ratio_string_plus(i)(1:long_ratio)) then   !
    op_trans(ligne) = ratio_real_pos(i)                                !
    ratio = .true.                                                   !
    exit                                                             !
   endif                                                             !
  end do                                                             !
                                                                     !
  if(.not. ratio) then                                               !
   do i=1,8                                                         !
    long_ratio = len_trim(ratio_string_neg(i))                       !
    if(new_string(1:long)==ratio_string_neg(i)(1:long_ratio)) then   !
     op_trans(ligne) = ratio_real_neg(i)                               !
     ratio = .true.                                                  !
     exit                                                            !
    endif                                                            !
   end do                                                            !
  endif                                                              !
 !-------------------------------------------------------------------!

  if(.not. ratio) READ(new_string,*) op_trans(ligne)
 else
  op_trans(ligne) = 0.
 endif


 return
end subroutine get_op_TRANS

!-----------------------------------------------------------------------------------------------

subroutine get_new_string(new_string, k)
 implicit none
  CHARACTER (LEN=*), INTENT(INOUT)   :: new_string
  INTEGER          , INTENT(IN)      :: k
  CHARACTER(LEN=12)                  :: reduced_string
  CHARACTER(LEN=1)                   :: string



 IF(k/=0) then
  if (k == 1) then
   reduced_string = new_string(2:)
  else
   READ(new_string(k-1:k-1),*) string
   IF(string=='+' .OR. string=='-') then
    IF(k-1 ==1) then
     reduced_string = new_string(k+1:)
    else
     reduced_string = new_string(1:k-2) // new_string(k+1:)
    endif
   else
    READ(new_string(k-2:k-2),*) string
    IF(string=='+' .OR. string=='-') then
     IF(k-2 ==1) then
      reduced_string = new_string(k+1:)
     else
      reduced_string = new_string(1:k-3) // new_string(k+1:)
     endif
    endif
   endif
  end if
  new_string =  reduced_string
 endif


 RETURN
end subroutine get_new_string
!--------------------------------------------------------------

subroutine write_symm_op_mat()
 USE  cryscal_module, ONLY :    keyword_SPGR, R => symm_op_rot,  T => symm_op_trans, nb_symm_op, symm_op_mat,   &
                                symm_op_string, symm_op_xyz, space_group_symbol, message_text, Mat_det,         &
								keyword_create_CIF
 USE  IO_module,      ONLY :    write_info
 implicit none
  integer                :: i, num
  CHARACTER (LEN=64)     :: operation_string


 call write_info('')
 call write_info('')
 if(keyword_SPGR) then
  call write_info('  >>> SYMMETRY OPERATORS FOR: '// TRIM(space_group_symbol)//' <<<')
  call write_info('      --------------------------------------------')
 else
  call write_info('  >>> SYMMETRY OPERATORS: <<<')
  call write_info('      ------------------')
 endif

 call write_info('')

 do num = 1, nb_symm_op
  call get_operation(operation_string, R(:,:,num), T(:,num))

  call write_info('')
  WRITE(message_text,'(a,i3,5a)') '   . SYMM #',num,':    ', (symm_op_string(i, num), i=1,3), '       ==> '//TRIM(operation_string)
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,F5.2)') '     _det = ', Mat_det
  call write_info(TRIM(message_text))
  call write_info('')

  !call get_operation(operation_string, R(:,:,num), T(:,num))
  !call write_info('       ==> '//TRIM(operation_string))
  !call write_info(' ')


  IF(symm_op_xyz) then
   write(message_text, '(5x,a,3(2x,3F4.0,2x,F6.3))')  ' * Matrix (r+t): ', (R(1,i,num),i=1,3), T(1,num),   &
                                                                           (R(2,i,num),i=1,3), T(2,num),   &
                                                                           (R(3,i,num),i=1,3), T(3,num)
   call write_info(TRIM(message_text))
  endif

  call write_info (' ')
  call write_info ('       > rotational part:')
  write(message_text,'(10x,a,3F9.4,a)') "|",R(1,1,num), R(1,2,num), R(1,3,num)," |"
   call write_info(TRIM(message_text))
  write(message_text,'(10x,a,3F9.4,a)') "|",R(2,1,num), R(2,2,num), R(2,3,num)," |"
   call write_info(TRIM(message_text))
  write(message_text,'(10x,a,3F9.4,a)') "|",R(3,1,num), R(3,2,num), R(3,3,num)," |"
   call write_info(TRIM(message_text))
  WRITE(message_text,'(a)') ' '
   call write_info(TRIM(message_text))
  WRITE(message_text,'(a)') '       > translational part:'
   call write_info(TRIM(message_text))
  write(message_text,'(10x,3F9.4)') (T(i, num), i=1,3)
   call write_info(TRIM(message_text))

  WRITE(message_text,'(a)') ' '
   call write_info(TRIM(message_text))
  WRITE(message_text,'(5x,4(a,F9.4))') ' new_x = ',R(1,1,num),' * x  + ',R(1,2,num), ' * y  + ',R(1,3,num), ' * z  + ', T(1,num)
   call write_info(TRIM(message_text))
  WRITE(message_text,'(5x,4(a,F9.4))') ' new_y = ',R(2,1,num),' * y  + ',R(2,2,num), ' * y  + ',R(2,3,num), ' * z  + ', T(2,num)
   call write_info(TRIM(message_text))
  WRITE(message_text,'(5x,4(a,F9.4))') ' new_z = ',R(3,1,num),' * z  + ',R(3,2,num), ' * y  + ',R(3,3,num), ' * z  + ', T(3,num)
   call write_info(TRIM(message_text))
  WRITE(message_text,'(a)') ' '
   call write_info(TRIM(message_text))

 end do

 
 IF(keyword_create_CIF) then
   call write_CIF_file('SPACE_GROUP')  
 endif

end subroutine write_symm_op_mat
!------------------------------------------------------------------------------------

subroutine apply_symm()
 USE cryscal_module, ONLY      : R => symm_op_rot,  T => symm_op_trans, nb_atom, atom_coord, atom_label, atom_typ,         &
                                 nb_symm_op, message_text,                                                                  &
                                 keyword_create_CIF, CIF_unit
 USE IO_module,      ONLY      : write_info
 implicit none
  !real                         :: new_x, new_y, new_z
  real, dimension(3)            :: new_coord
  INTEGER                      :: i, num


 !call write_info('')
 !call write_info('  > New coordinates:')
 !call write_info('')

 IF(keyword_create_CIF) then
  WRITE(CIF_unit, '(a)') 'loop_'
  WRITE(CIF_unit, '(a)') '_atom_site_label'
  WRITE(CIF_unit, '(a)') '_atom_site_type_symbol'
  WRITE(CIF_unit, '(a)') '_atom_site_fract_x'
  WRITE(CIF_unit, '(a)') '_atom_site_fract_y'
  WRITE(CIF_unit, '(a)') '_atom_site_fract_z'
  WRITE(CIF_unit, '(a)') '_atom_site_U_iso_or_equiv'
  WRITE(CIF_unit, '(a)') '_atom_site_adp_type'
  WRITE(CIF_unit, '(a)') '_atom_site_occupancy'
  WRITE(CIF_unit, '(a)') '_atom_site_symmetry_multiplicity'
 endif



 do i=1, nb_atom

  call write_info('')
  write(message_text,'(a,I3,a,a6,a,3(2x,F8.5))') ' * ATOM #',i, ': ', atom_label(i),': ', atom_coord(1,i), &
                                                 atom_coord(2,i), atom_coord(3,i)
  call write_info(trim(message_text))

  call write_info('')
  call write_info('  > New coordinates:')
  call write_info('')

  do num = 1, nb_symm_op
   call Get_new_coord(i, num, new_coord)
   
   !new_coord(1) = R(1,1, num) * atom_coord(1, i)  +  R(1,2, num) * atom_coord(2, i) +  R(1,3, num) * atom_coord(3,i) &
   !             + T(1, num)
   !new_coord(2) = R(2,1, num) * atom_coord(1, i)  +  R(2,2, num) * atom_coord(2, i) +  R(2,3, num) * atom_coord(3,i) &
   !             + T(2, num)
   !new_coord(3) = R(3,1, num) * atom_coord(1, i)  +  R(3,2, num) * atom_coord(2, i) +  R(3,3, num) * atom_coord(3,i) &
   !             + T(3, num)

   if(num<10) then
    write(message_text,'(a,I3,a,a6,a,i1,a6,3(1x,F8.5))') '  Symmetry operator # ',num, ':  ATOM: ', &
	                                                   trim(atom_label(i)),'_',num,                 &
                                                       trim(atom_typ(i)),  new_coord(1:3)
   elseif(num < 100) then
    write(message_text,'(a,I3,a,a6,a,i2,a6,3(1x,F8.5))') '  Symmetry operator # ',num, ':  ATOM: ', &
	                                                   trim(atom_label(i)),'_',num,                 &
                                                       trim(atom_typ(i)),  new_coord(1:3)
   else
    write(message_text,'(a,I3,a,a6,a,i3,a6,3(1x,F8.5))') '  Symmetry operator # ',num, ':  ATOM: ', &
	                                                   trim(atom_label(i)),'_',num,                 &
                                                       trim(atom_typ(i)),  new_coord(1:3)
   endif

    call write_info(TRIM(message_text))
    call write_info('')
   IF(keyword_create_CIF) then
    IF(num<10) then
     WRITE(CIF_unit, '(a,a,i1,a6, 3(1x,F8.5),a)') trim(atom_label(i)),'_',num, trim(atom_typ(i)),  & 
	                                              new_coord(1:3), '  0.05   Uiso  1  1'
    ELSEif(num<100) then
     WRITE(CIF_unit, '(a,a,i2,a6, 3(1x,F8.5),a)') trim(atom_label(i)),'_',num, trim(atom_typ(i)),  & 
	                                              new_coord(1:3), '  0.05   Uiso  1  1'
    else
     WRITE(CIF_unit, '(a,a,i3,a6, 3(1x,F8.5),a)') trim(atom_label(i)),'_',num, trim(atom_typ(i)),  & 
	                                              new_coord(1:3), '  0.05   Uiso  1  1'
    endif
   endif
  END do
 end do


 call write_info('')


end subroutine apply_symm

!--------------------------------------------------------------------------------------
subroutine get_operation(operation_string, mat_33, translation)
 USE cryscal_module, ONLY :  T => symm_op_trans, Mat_det, message_text
 use matrix_module,  only :  matrix_determinant_33
 USE IO_module,      ONLY :  write_info

 implicit none
  CHARACTER(LEN=64),    INTENT(INOUT)   :: operation_string
  REAL, DIMENSION(3,3), INTENT(IN)      :: mat_33
  REAL, DIMENSION(3),   INTENT(IN)      :: translation
  REAL                                  :: sg_translation
  LOGICAL                               :: identity, inversion
  LOGICAL                               :: t_ok

 ! calcul du determinant de la partie rotationnelle de la matrice
 ! calcul du determinant de la matrice suivant la premiere ligne

 Mat_det = matrix_determinant_33(Mat_33)


 if (abs(Mat_det -1) > 0.0001 .and. ABS(Mat_det+1) > 0.0001) then
  WRITE(message_text, '(a,F5.2)') '  Mat_det: ', Mat_det
  call write_info(TRIM(message_text))
  operation_string = 'No valid symmetry operation'
  return
 endif

 sg_translation = translation(1) + translation(2) + translation(3)
 !WRITE(message_text,*)  ' translation: ', sg_translation
 !call write_info(trim(message_text))

 call matrice_identity(Mat_33, identity, inversion)

 t_ok = .false.
 IF((mat_33(1,1) > 0. .AND. translation(1) > 0.) .or.   &
    (mat_33(2,2) > 0. .AND. translation(2) > 0.) .or.   &
    (mat_33(3,3) > 0. .AND. translation(3) > 0.)) then
     t_ok = .true.
 endif

 !IF(ABS(sg_translation) < 0.0001) then ! pas de translation
 IF(.NOT. t_ok) then
  IF(identity) then
   operation_string = 'Identity'
  ELSEIF(inversion .and. sg_translation < 0.0001) then
   operation_string = 'Inversion'
  elseif(ABS(Mat_det - 1) < 0.0001) THEN        ! det =  1
   operation_string = 'Pure rotation'
  ELSEIF(ABS(Mat_det+1) < 0.0001) THEN          ! det = -1
   operation_string = 'Rotoinversion'
  endif

 else
  call matrice_identity(Mat_33, identity, inversion)
  IF(identity) then
    operation_string = 'Pure translation'
  else
   IF(ABS(Mat_det - 1) < 0.0001) THEN            ! det =  1
    operation_string = 'Screw axis'
   ELSEIF(ABS(Mat_det+1) < 0.0001) THEN          ! det = -1
    operation_string = 'Glide plane'
   endif
  end if
 endif                              ! translation

 RETURN
end subroutine get_operation

!---------------------------------------------------------------------------

subroutine matrice_identity(Mat_33, identity, inversion)
 implicit none
  REAL, DIMENSION(3,3), INTENT(IN)     :: Mat_33
  LOGICAL,              INTENT(INOUT)  :: identity, inversion

  identity  = .false.
  inversion = .false.
  if (ABS(Mat_33(1,1) - 1.) < 0.0001 .and. ABS(Mat_33(2,2) - 1.) < 0.0001 .and. ABS(Mat_33(3,3) - 1.) < 0.0001) then
   if (ABS(Mat_33(1,2)) < 0.0001 .AND. ABS(Mat_33(1,3)) < 0.0001) then
    if (ABS(Mat_33(2,1)) < 0.0001 .AND. ABS(Mat_33(2,3)) < 0.0001) then
     if (ABS(Mat_33(3,1)) < 0.0001 .AND. ABS(Mat_33(3,2)) < 0.0001) then
      identity = .true.
      return
     endif
    endif
   endif
  endif

  if (ABS(Mat_33(1,1) + 1.) < 0.0001 .and. ABS(Mat_33(2,2) + 1.) < 0.0001 .and. ABS(Mat_33(3,3) + 1.) < 0.0001) then
   if (ABS(Mat_33(1,2)) < 0.0001 .AND. ABS(Mat_33(1,3)) < 0.0001) then
    if (ABS(Mat_33(2,1)) < 0.0001 .AND. ABS(Mat_33(2,3)) < 0.0001) then
     if (ABS(Mat_33(3,1)) < 0.0001 .AND. ABS(Mat_33(3,2)) < 0.0001) then
      inversion = .true.
      return
     endif
    endif
   endif
  endif


 RETURN
end subroutine matrice_identity

!----------------------------------------------------------------------------
subroutine Get_new_coord(i, num, new_coord)
 USE cryscal_module, ONLY      : R => symm_op_rot,  T => symm_op_trans, atom_coord
 
 implicit none
  integer,            intent(in)      :: i                ! numero de l'atome
  integer,            intent(in)      :: num              ! numero de l'operateur de sym.
  real, dimension(3), intent(out)     :: new_coord
  
  
  
  new_coord(1) = R(1,1, num) * atom_coord(1, i)  +  R(1,2, num) * atom_coord(2, i) +  R(1,3, num) * atom_coord(3,i)  + T(1, num)
  new_coord(2) = R(2,1, num) * atom_coord(1, i)  +  R(2,2, num) * atom_coord(2, i) +  R(2,3, num) * atom_coord(3,i)  + T(2, num)
  new_coord(3) = R(3,1, num) * atom_coord(1, i)  +  R(3,2, num) * atom_coord(2, i) +  R(3,3, num) * atom_coord(3,i)  + T(3, num)


end subroutine Get_new_coord