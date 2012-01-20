
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
 USE cryscal_module,        ONLY : keyword_THERM_SHELX, THERM_uij, THERM_bij, THERM_beta, therm_values,  unit_cell,   &
                                   pi, message_text, crystal_cell
 USE CFML_Crystal_Metrics,  only : U_equiv, Convert_U_B, convert_U_betas, convert_B_U, convert_B_Betas, &
                                   convert_Betas_U, convert_Betas_B
 USE CFML_Math_General,     ONLY : acosd 
 USE CFML_math_3D,          ONLY : matrix_diageigen
 USE CFML_GlobalDeps,       ONLY : sp
 use CFML_Crystal_Metrics,  only : Set_Crystal_Cell   ! Crystal_cell_type

 
 USE IO_module

 implicit none
  integer                          :: i, j, ineg
  REAL(kind=sp), DIMENSION(6)      :: Xij, tmp_Xij
  REAL(kind=sp), DIMENSION(3,3)    :: M_U, eigen, beta, aux
  REAL(kind=sp), DIMENSION(3,3)    :: L, LT
  REAL(kind=sp), DIMENSION(3)      :: rms
  REAL, DIMENSION(6)               :: new_ADP, new_ADP2
  REAL(kind=sp)                    :: Ueq, Uiso

  call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell, 'A')

  !write(unit=message_text, fmt="(a, 1x, 3F10.5)")  ' CELL :  ', crystal_cell%cell(1:3)
  !call write_info(trim(message_text))
  !write(unit=message_text, fmt="(a, 1x, 3F10.5)")  ' ANG  :  ', crystal_cell%ang(1:3)
  !call write_info(trim(message_text))
 
 
 
  L=crystal_cell%Cr_Orth_cel   !Cell transforming  column vectors in direct lattice to Cartesian basis
  LT=Transpose(L)
  
  !write(message_text, '(a, 3(x, 3F10.5))')  ' L :  ', L(1:3,1),  L(1:3, 2),  L(1:3,3)
  !call write_info(trim(message_text))
  !write(message_text, '(a, 3(x, 3F10.5))')  ' LT:  ', LT(1:3,1), LT(1:3, 2), LT(1:3,3)
  !call write_info(trim(message_text))

 call write_info('')
 call write_info('   >> ANISOTROPIC THERMAL PARAMETERS CONVERSION:')
 call write_info('')

 Xij(1:6) = therm_values(1:6)
 
 if(keyword_THERM_SHELX) then ! 11 22 33 23 13 12
  tmp_Xij = Xij
  Xij(1:3) = tmp_Xij(1:3)
  Xij(4) = tmp_Xij(6)
  Xij(5) = tmp_Xij(5)
  Xij(6) = tmp_Xij(4)
 endif
 
 call write_info(trim(message_text))
 if(therm_Uij) then
  Xij = convert_U_betas(Xij, crystal_cell) 
  THERM_beta = .true.
  THERM_Uij  = .false.
  THERM_Bij  = .false.
 elseif(THERM_Bij) then
  Xij= convert_B_betas(Xij, crystal_cell) 
  THERM_beta = .true.
  THERM_Bij  = .false.
  THERM_Uij  = .false.  
 endif 

 IF(THERM_Uij) then
  
  
  
  M_U = reshape((/Xij(1),Xij(4),Xij(5), Xij(4),Xij(2),Xij(6), Xij(5),Xij(6),Xij(3) /),(/3,3/))
  M_U = M_U*1.E04
!     write(message_text, '(a, 1(x, 3F10.5))')  ' M_U :  ', M_U(1:3,1)
!	 call write_info(trim(message_text))
!     write(message_text, '(a, 1(x, 3F10.5))')  ' M_U :  ', M_U(1:3,2)
!	 call write_info(trim(message_text))
!     write(message_text, '(a, 1(x, 3F10.5))')  ' M_U :  ', M_U(1:3,3)
!	 call write_info(trim(message_text))
!     call write_info('')
	 
 call matrix_diageigen(M_U, rms,eigen)
 
 ! do i=1, 3  
 !  if(rms(i) < 0.0) then
 !   write(message_text, '(a)') "   -> Matrix U non-positive definite!"
 !   call write_info(trim(message_text))
 !   call write_info('')
 !   ineg = 1
 !  else
 !   ineg = 0 
 !   !rms(i) = sqrt(rms(i))
 !   !write(message_text, '(1(2x,F8.5))') rms(i) 
 !   !call write_info(trim(message_text))
 !  endif
 ! end do
  
 ! write(message_text,'(a)') '    U-Eigen Value (A**2) ----      Eigen vector(Orth. syst.)'
 ! call write_info(trim(message_text))
 ! do i=1, 3  
 ! WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)
 ! call write_info(trim(message_text))
 ! enddo
  
 ! write(message_text, '(3(2x,F8.5))') eigen(1:3, 1) 
 ! call write_info(trim(message_text))
 ! write(message_text, '(3(2x,F8.5))') eigen(1:3, 2) 
 ! call write_info(trim(message_text))
 ! write(message_text, '(3(2x,F8.5))') eigen(1:3, 3) 
 ! call write_info(trim(message_text))
  
  !rms(1:3) = sqrt(rms(1:3))
  !write(message_text, '(5x, a, 3F8.5)') " RMS (A)  :   " , rms(1:3) 
  !call write_info(trim(message_text))
  !call write_info('')
  
  write(message_text,'(a)') '    U-Eigen Value (A**2) ----      Eigen vector(Orth. syst.)'
  call write_info(trim(message_text))
  do i=1, 3  
  rms(i) = rms(i)*1E-04
  WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)
  call write_info(trim(message_text))
  enddo
  call write_info('')
  
  ineg = 0
  Uiso = 0
  do i=1,3
   if(rms(i) < 0.) then
    ineg = 1
	WRITE(message_text,'(a)')'   -> Matrix U non-positive definite! '
	call write_info(trim(message_text)) 
	call write_info('')
   else
    Uiso= Uiso + rms(i)
    rms(i) = sqrt(rms(i))
   endif
  end do   
  Uiso = Uiso / 3.
  
  DO i=1,3
   DO j=1,3
    M_U(j,i)=L(i,j)/Crystal_cell%cell(j)
   END DO
  END DO
	
	
  if(ineg == 0) then
   write(message_text, '(a)') '     R.M.S. in ANGSTROMS ------       Angles with A, B, C'
   call write_info(trim(message_text))
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
  

  !Ueq = U_equiv(Crystal_cell, Xij)
  WRITE(message_text,'(a,6F10.5)') '   >> U_ij (A2):   ', Xij
  call write_info(TRIM(message_text))
  
  new_ADP = convert_U_B(Xij)
  WRITE(message_text,'(a,6F10.5)') '   >> B_ij (A2):   ', new_ADP
  call write_info(TRIM(message_text))

  new_ADP = convert_U_betas(Xij, crystal_cell)
  WRITE(message_text,'(a,6F10.5)') '   >> Beta_ij:     ', new_ADP
  call write_info(TRIM(message_text))
  
  Ueq = U_equiv(Crystal_cell, Xij)
  WRITE(message_text,'(a,F10.5)')  '   >> Ueq  (A2):   ', Ueq
  call write_info(TRIM(message_text)) 
  WRITE(message_text,'(a,F10.5)')  '   >> Uiso (A2):   ', Uiso
  call write_info(TRIM(message_text)) 
  WRITE(message_text,'(a,F10.5)')  '   >> Beq  (A2):   ', Ueq*8.*pi**2.
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F10.5)')  '   >> Biso (A2):   ', Uiso*8.*pi**2.
  call write_info(TRIM(message_text))


 ELSEIF(therm_Bij) then
  new_ADP = convert_B_U(Xij)
  M_U = reshape((/new_ADP(1),new_ADP(4),new_ADP(5), new_ADP(4),new_ADP(2),new_ADP(6),new_ADP(5),new_ADP(6),new_ADP(3) /),(/3,3/))
  M_U = M_U*1.e04
  
   !write(message_text, '(a, 1(x, 3F10.5))')  ' M_U :  ', M_U(1:3,1)
   !call write_info(trim(message_text))
   !write(message_text, '(a, 1(x, 3F10.5))')  ' M_U :  ', M_U(1:3,2)
   !call write_info(trim(message_text))
   !write(message_text, '(a, 1(x, 3F10.5))')  ' M_U :  ', M_U(1:3,3)
   !call write_info(trim(message_text))

  call matrix_diageigen(M_U, rms,eigen)
   
  
!  write(message_text,'(a)') '    U-Eigen Value (A**2) ----      Eigen vector(Orth. syst.)'
!  call write_info(trim(message_text))
!  do i=1, 3  
!  WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)
!  call write_info(trim(message_text))
!  enddo

  write(message_text,'(a)') '    U-Eigen Value (A**2) ----      Eigen vector(Orth. syst.)'
  call write_info(trim(message_text))
  do i=1, 3 
  rms(i) = rms(i) * 1.E-04  
  WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)
  call write_info(trim(message_text))
  enddo
  call write_info('')
 
  ineg = 0 
  Uiso = 0
  do i=1,3
   if(rms(i) < 0.) then
    ineg = 1
	WRITE(message_text,'(a)')'   -> Matrix U non-positive definite! '
	call write_info(trim(message_text)) 
	call write_info('')
   else    
    Uiso = Uiso + rms(i)
    rms(i) = sqrt(rms(i))
   endif
  end do   
  Uiso = Uiso / 3.
  
  

  DO i=1,3
   DO j=1,3
    M_U(j,i)=L(i,j)/Crystal_cell%cell(j)
   END DO
  END DO
 
  
  if(ineg == 0) then
   write(message_text, '(a)') '     R.M.S. in ANGSTROMS ------       Angles with A, B, C'
   call write_info(trim(message_text)) 
   LT = MATMUL(M_U, eigen)
   LT=acosd(LT)
   DO i=1,3
   WRITE(message_text,'((5x,F10.5,a,3(1x,F10.3)))')    rms(i),'          ----',(LT(j,i),j=1,3)
   call write_info(trim(message_text))
   END DO
   call write_info('')
  end if
  

 

  write(message_text, '(a)')      '                        11        22        33        12        13        23'
  call write_info(trim(message_text))
  call write_info('')
  
  WRITE(message_text,'(a,6F10.5)') '   >> B_ij (A2):   ', Xij
  call write_info(TRIM(message_text))

  new_ADP = convert_B_U(Xij)
  WRITE(message_text,'(a,6F10.5)') '   >> U_ij (A2):   ', new_ADP
  call write_info(TRIM(message_text))
  Ueq = U_equiv(Crystal_cell, new_ADP)
  
  new_ADP = convert_B_betas(Xij, crystal_cell)
  WRITE(message_text,'(a,6F10.5)') '   >> Beta_ij  :   ', new_ADP
  call write_info(TRIM(message_text))
  
  WRITE(message_text,'(a,F10.5)')  '   >> Ueq  (A2):   ', Ueq
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F10.5)')  '   >> Uiso (A2):   ', Uiso
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F10.5)')  '   >> Beq  (A2):   ', Ueq*8.*pi**2.
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F10.5)')  '   >> Biso (A2):   ', Uiso*8.*pi**2.
  call write_info(TRIM(message_text))


  

 ELSEIF(therm_BETA) then
 
 !beta=reshape((/Xij(1),Xij(4),Xij(5), Xij(4),Xij(2),Xij(6), Xij(5),Xij(6), Xij(3) /),(/3,3/))
 !beta=beta*0.5/pi/pi
 !beta=matmul(matmul(crystal_cell%Cr_Orth_cel,beta),transpose(crystal_cell%Cr_Orth_cel))
 !call matrix_diageigen(beta,rms,eigen)
 !write(message_text,fmt="(a)") "               U-Eigen Value(A**2) ----       Eigen vector(Orth. syst.)     R.M.S (Angstroms)"
 !call write_info(trim(message_text))
 !do j =1,3
 ! if (rms(j) < 0.0)  then
 !  write(message_text,fmt="((t16,f10.5,a,3(tr1,f10.5),a))")       rms(j), "          --- ", eigen(:,j),"   -> Matrix U non-positive definite!"
 ! else
 !  write(message_text,fmt="((t16,f10.5,a,3(tr1,f10.5),a,f14.5))") rms(j), "          ---(", eigen(:,j),")", sqrt(rms(j))
 ! end if
 ! call write_info(trim(message_text))
 !end do
 
  
 
  DO j=1,3                  !Transformation to 3x3 symmetric matrix
   beta(j,j)=Xij(j)
  END DO
  beta(1,2)=Xij(4)
  beta(1,3)=Xij(5)
  beta(2,3)=Xij(6)  
  beta(2,1)=beta(1,2)
  beta(3,1)=beta(1,3)
  beta(3,2)=beta(2,3)
  
  
 ! 
 !DO i=1,3
 !       DO j=1,3
 !         !u(j,i)=l(i,j)/crystal_cell%cell(j)
 !         beta(i,j)=beta(i,j)/2.0/pi/pi
 !       END DO
 !     END DO 
  
  !write(message_text, '(3(2x,3F8.5))') beta(1:3, 1), beta(1:3,2), beta(1:3,3)
  !call write_info(trim(message_text))
  !aux=matmul(L,beta)
  !beta=matmul(aux,LT)
  ! write(message_text, '(3(2x,3F8.5))') beta(1:3, 1), beta(1:3,2), beta(1:3,3)
  !call write_info(trim(message_text))
  !CALL matrix_diageigen(beta,rms,aux)

  !  do i=1, 3  
  !WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(aux(j,i),j=1,3)
  !call write_info(trim(message_text))
  !enddo

 
 
 
  !new_ADP2 = convert_Betas_B(Xij, crystal_cell)
  !new_ADP = convert_B_U(new_ADP2)
  !M_U = reshape((/new_ADP(1),new_ADP(4),new_ADP(5), new_ADP(4),new_ADP(2),new_ADP(6),new_ADP(5),new_ADP(6),new_ADP(3) /),(/3,3/))
  !call matrix_diageigen(M_U, rms, eigen)
 
 
 !write(message_text, '(3(1x,3F8.5))') M_U(1:3, 1), M_U(1:3, 2), M_U(1:3, 3)
 !call write_info(trim(message_text))
 !   do i=1, 3  
 ! WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)
 ! call write_info(trim(message_text))
 ! enddo

 beta=beta*1E4
 !write(message_text, '(a, 3(x, 3F10.5))')  ' beta :  ', beta(1:3,1),  beta(1:3, 2),  beta(1:3,3)
 !call write_info(trim(message_text))

 DO i=1,3
  DO j=1,3
   M_U(j,i)=L(i,j)/Crystal_cell%cell(j)
   beta(i,j)=beta(i,j)/2.0/pi/pi
  END DO
 END DO
 
 aux=matmul(L,beta)
 beta=matmul(aux,LT)
 CALL matrix_diageigen(beta,rms,eigen)

 !write(message_text, '(a, 1(x, 3F10.5))')  ' beta :  ', beta(1:3,1)
 !call write_info(trim(message_text))
 !write(message_text, '(a, 1(x, 3F10.5))')  ' beta :  ', beta(1:3,2)
 !call write_info(trim(message_text))
 !write(message_text, '(a, 1(x, 3F10.5))')  ' beta :  ', beta(1:3,3)
 !call write_info(trim(message_text))

  write(message_text,'(a)') '    U-Eigen Value (A**2) ----      Eigen vector(Orth. syst.)'
  call write_info(trim(message_text))
  do i=1, 3  
  rms(i) = rms(i)*1.E-04
  WRITE(message_text,'((5x,F10.5,a,3(1x,F10.5)))')    rms(i),'          ----',(eigen(j,i),j=1,3)
  call write_info(trim(message_text))
  enddo
  call write_info('')
  
  ineg = 0
  Uiso = 0
  do i=1,3
   if(rms(i) < 0.) then
    ineg = 1
	WRITE(message_text,'(a)')'   -> Matrix U non-positive definite! '
	call write_info(trim(message_text)) 
	call write_info('')
   else
    Uiso= Uiso + rms(i)
    rms(i) = sqrt(rms(i))	
   endif
  end do  
  Uiso = Uiso/3.  
  
  

  
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
  
  
  !write(message_text, '(5x, a, 3F10.5)') " RMS (A)  :   " , rms(1:3) 
  !call write_info(trim(message_text))
  !call write_info('')

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
  !WRITE(message_text,'(a,F10.5)')  '   >> Biso (A2):   ', Uiso*8.*pi**2.
  !call write_info(TRIM(message_text))

 endif

end subroutine calc_THERM_ANISO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diag_matrice
! diagonalisation d'une matric 3*3
 USE cryscal_module,        ONLY : Mat, message_text
 USE CFML_math_3D,          ONLY : matrix_diageigen
 USE IO_module
 USE CFML_GlobalDeps,       ONLY : sp
 
 implicit none
  REAL(kind=sp), DIMENSION(3,3)    :: M_U, eigen
  REAL(kind=sp), DIMENSION(3)      :: rms
  integer                          :: i, j
 
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