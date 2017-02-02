
!-----------------------------------------------------------------

subroutine calc_therm_iso()
 USE cryscalc_module, ONLY : THERM, pi, message_text, debug_proc
 USE IO_module

 implicit none
  integer                  :: i
  real                     :: new_therm_value

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calc_therm_iso")

 call write_info('')
 call write_info('   >> ISOTROPIC THERMAL PARAMETERS CONVERSION:')
 call write_info('')

 if (THERM%Biso) then      ! Biso --> Uiso
  call write_info('         Biso         Uiso')

  do i=1, THERM%nb_values
   new_therm_value = therm%values(i) / (8*pi**2)
   WRITE(message_text, '(2(5x,F8.4))') therm%values(i), new_therm_value
   call write_info(TRIM(message_text))
  end do
  return

 ELSEIF(THERM%Uiso) then   ! Uiso --> Biso
  call write_info('         Uiso         Biso')

  do i=1, THERM%nb_values
   new_therm_value = therm%values(i) * (8*pi**2)
   WRITE(message_text,'(2(5x,F8.4))') therm%values(i), new_therm_value
   call write_info(TRIM(message_text))
  end do
  return
 end if



 RETURN
end subroutine calc_therm_iso

!--------------------------------------------------------------------------

subroutine calc_THERM_ANISO()
 USE cryscalc_module,       ONLY : keyword_THERM_SHELX, THERM,  unit_cell,   &
                                   pi, message_text, crystal_cell, ADP_details, ADP_comment, debug_proc
 USE CFML_Crystal_Metrics,  only : U_equiv, Convert_U_B, convert_U_betas, convert_B_U, convert_B_Betas, &
                                   convert_Betas_U, convert_Betas_B
 USE CFML_Math_General,     ONLY : acosd
 USE CFML_math_3D,          ONLY : matrix_diageigen
 USE CFML_GlobalDeps,       ONLY : sp
 use CFML_Crystal_Metrics,  only : Set_Crystal_Cell   ! Crystal_cell_type
 USE IO_module

 implicit none
  REAL(kind=sp), DIMENSION(3)      :: rms
  
  if(debug_proc%level_2)  call write_debug_proc_level(2, "calc_therm_aniso")

  call write_info('')
  call write_info('   >> ANISOTROPIC THERMAL PARAMETERS CONVERSION:')
  call write_info('')

  ADP_details = .true.
  call check_ADP(therm%values, rms, ADP_comment, ADP_details)



 return
end subroutine calc_THERM_ANISO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diag_matrice
! diagonalisation d'une matric 3*3
 USE cryscalc_module,       ONLY : Mat, message_text, debug_proc
 USE CFML_math_3D,          ONLY : matrix_diageigen
 USE IO_module
 USE CFML_GlobalDeps,       ONLY : sp

 implicit none
  REAL(kind=sp), DIMENSION(3,3)    :: M_U, eigen
  REAL(kind=sp), DIMENSION(3)      :: rms
  integer                          :: i, j

  if(debug_proc%level_2)  call write_debug_proc_level(2, "diag_matrice")

  M_U = Mat
  call matrix_diageigen(M_U, rms,eigen)

  write(message_text,'(a)') '      Eigen Value        ----      Eigen vector             '
  call write_info(trim(message_text))
  do i=1, 3
  rms(i) = rms(i)
  WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)
  call write_info(trim(message_text))
  enddo
  call write_info('')

 return
end subroutine diag_matrice


!----------------------------------------------------------------------------
subroutine CHECK_ADP(ADP, rms, ADP_comment, details)
 USE CFML_math_3D,          ONLY : matrix_diageigen
 USE CFML_GlobalDeps,       ONLY : sp
 USE cryscalc_module,       ONLY : crystal_cell, pi, keyword_THERM_SHELX, THERM, message_text
 USE CFML_Crystal_Metrics,  only : U_equiv, convert_U_betas, convert_B_betas, convert_Betas_U, convert_Betas_B
 USE CFML_Math_General,     ONLY : acosd
 USE IO_module

 implicit none
  REAL(kind=sp), DIMENSION(6), intent(inout)      :: ADP
  REAL(kind=sp), DIMENSION(3), intent(out)        :: rms
  LOGICAL,                     intent(in)         :: details
  REAL(kind=sp), DIMENSION(6)      :: Xij, tmp_Xij
  CHARACTER(LEN=48), intent(INOUT) :: ADP_comment
  integer                          :: i, j, ineg
  integer, dimension(1)            :: i1, i2, i3

  REAL(kind=sp), DIMENSION(3,3)    :: M_U, eigen, beta, aux
  REAL(kind=sp), DIMENSION(3,3)    :: L, LT
  !REAL(kind=sp), DIMENSION(3)      :: rms
  REAL, DIMENSION(6)               :: new_ADP
  REAL(kind=sp)                    :: Ueq, Uiso


  call create_CELL_object
  L=crystal_cell%Cr_Orth_cel   !Cell transforming  column vectors in direct lattice to Cartesian basis
  LT=Transpose(L)

 Xij = ADP
 if(keyword_THERM_SHELX) then ! 11 22 33 23 13 12
  tmp_Xij = Xij
  Xij(1:3) = tmp_Xij(1:3)
  Xij(4) = tmp_Xij(6)
  Xij(5) = tmp_Xij(5)
  Xij(6) = tmp_Xij(4)
 endif
 !write(*,*) 'Xij : ', ADP


  !Xij = convert_U_betas(Xij, crystal_cell)
 if(therm%Uij) then
  Xij = convert_U_betas(Xij, crystal_cell)
  THERM%beta = .true.
  THERM%Uij  = .false.
  THERM%Bij  = .false.
 elseif(THERM%Bij) then
  Xij= convert_B_betas(Xij, crystal_cell)
  THERM%beta = .true.
  THERM%Bij  = .false.
  THERM%Uij  = .false.
 endif



  DO j=1,3                  !Transformation to 3x3 symmetric matrix
   beta(j,j)=Xij(j)
  END DO
  beta(1,2)=Xij(4)
  beta(1,3)=Xij(5)
  beta(2,3)=Xij(6)
  beta(2,1)=beta(1,2)
  beta(3,1)=beta(1,3)
  beta(3,2)=beta(2,3)

  beta=beta*1E4
  !write(*, '(a, 3(x, 3F10.5))')  ' beta :  ', beta(1:3,1),  beta(1:3, 2),  beta(1:3,3)


 DO i=1,3
  DO j=1,3
   M_U(j,i)=L(i,j)/Crystal_cell%cell(j)
   beta(i,j)=beta(i,j)/2.0/pi/pi
  END DO
 END DO

 aux=matmul(L,beta)
 beta=matmul(aux,LT)
 CALL matrix_diageigen(beta,rms,eigen)




 !ADP_comment = ""
 !do i=1,3
 !rms(i) = rms(i)*1.E-04
  !WRITE(*,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)

  !if(rms(i) < 0.) then
  ! ADP_comment = '-> Matrix U non-positive definite! '
  !end if
 !end do

 !if(.not. details) return

  if(details) then
  write(message_text,'(a)') '     U-Eigen Value (A**2) ----      Eigen vector(Orth. syst.)'
  call write_info(trim(message_text))
  end if
  ADP_comment = ''
  ineg = 0
  do i=1, 3
   rms(i) = rms(i)*1.E-04
   if(rms(i) < 0.) then
    ADP_comment = '-> Matrix U non-positive definite!'
    ineg = 1
   end if
   if(details) then
    WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)
    call write_info(trim(message_text))
    !WRITE(message_text,'(2a)')'   ', trim(ADP_comment)
    !call write_info(trim(message_text))
    !call write_info('')
   end if
  enddo
  if(.not. details) return


  call write_info('')

  ineg = 0
  Uiso = 0
  !ADP_comment = ''
  do i=1,3
   if(rms(i) < 0.) then
    ineg =1
    WRITE(message_text,'(2a)')'   ', trim(ADP_comment)
    call write_info(trim(message_text))
    call write_info('')
    return
   else
    Uiso= Uiso + rms(i)
    rms(i) = sqrt(rms(i))
   endif
  end do
  Uiso = Uiso/3.


  write(message_text, '(a)') '     Principal mean square atomic displacements U'
  call write_info(trim(message_text))

  i1(1) = 0
  i2(1) = 0
  i3(1) = 0
  i1 = maxloc(rms)
  i3 = minloc(rms)
  if(i1(1) == 1) then
   if(i3(1) == 2) then
    i2(1) = 3
   else
    i2(1) = 2
   endif
  elseif(i1(1)==2) then
   if(i3(1)==1) then
    i2(1)=3
   else
    i2(1)=1
   endif
  elseif(i1(1)==3) then
   if(i3(1)==2) then
    i2(1)=1
   else
    i2(1)=2
   endif
  endif

  write(message_text, '(5x,3F10.5)') rms(i1(1))**2, rms(i2(1))**2, rms(i3(1))**2
  call write_info(trim(message_text))
  call write_info('')



  if(ineg == 0) then
   write(message_text, '(a)') '     R.M.S. in ANGSTROMS ------       Angles with A, B, C'
   call write_info(trim(message_text))
   !LT=acosd(LT)
   LT=matmul(M_U, eigen)
   DO i=1,3
    DO j=1,3
     IF(ABS(LT(i,j)) > 1.0) LT(i,j)=sign(1.0,LT(i,j))
     LT(i,j)=acosd(LT(i,j))
    END DO
   END DO

   DO i=1,3
   WRITE(message_text,'((5x,F10.5,a,3(1x,F10.3)))')    rms(i),'          ----',(LT(j,i),j=1,3)
   call write_info(trim(message_text))
   END DO
   call write_info('')
  end if


  write(message_text, '(a)')      '                        11        22        33        12        13        23'
  call write_info(trim(message_text))
  call write_info('')



  WRITE(message_text,'(a,6F10.5)') '   >> Beta_ij:     ', Xij
  call write_info(TRIM(message_text))

  new_ADP = convert_Betas_B(Xij, crystal_cell)
  WRITE(message_text,'(a,6F10.5)') '   >> B_ij (A2):   ', new_ADP
  call write_info(TRIM(message_text))

  new_ADP = convert_Betas_U(Xij, crystal_cell)
  WRITE(message_text,'(a,6F10.5)') '   >> U_ij (A2):   ', new_ADP
  call write_info(TRIM(message_text))

  Ueq = U_equiv(Crystal_cell, new_ADP)
  WRITE(message_text,'(a,F10.5)')  '   >> Ueq  (A2):   ', Ueq
  call write_info(TRIM(message_text))
  !WRITE(message_text,'(a,F10.5)')  '   >> Uiso (A2):   ', Uiso
  !call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F10.5)')  '   >> Beq  (A2):   ', Ueq*8.*pi**2.
  call write_info(TRIM(message_text))
  call write_info('')
  !WRITE(message_text,'(a,F10.5)')  '   >> Biso (A2):   ', Uiso*8.*pi**2.
  !call write_info(TRIM(message_text))


 !write(*,*) ' ADP comment : ', trim(ADP_comment)

 return
end subroutine CHECK_ADP
