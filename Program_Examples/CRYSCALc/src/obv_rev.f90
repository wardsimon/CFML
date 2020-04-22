!     Last change:  TR   13 Jul 2006    4:22 pm

subroutine twin_obverse_reverse

 USE cryscalc_module, ONLY : keyword_file, HKL_unit, HKLF5_unit, twin_matrix => OBV_REV_twin_matrix, message_text, debug_proc
 USE macros_module,   ONLY : multiple
 USE HKL_module,      ONLY : n_ref, h,k,l, F2, sig_F2, cos_dir, HKL_file
 USE IO_module,       ONLY : write_info

!--------------------------------------------------------------
!
! lecture d'un fichier au format SHELX (HKLF4):
!  analyse des reflections:
!   . -h +k +l = 3n  >>> refl. type obverse
!   . +h -k +l = 3n  >>> refl. type reverse
!   . -h +k +l <> 3n et
!   . +h -k +l <> 3n    >>> refl. ni obverse ni reverse
!

!  >> calcul des intensités moyennes <I_obv>, <I_rev>, <I_NoNr>
!  >> calcul des intensités/sigmas moyennes <I_obv/sig>, <I_rev>/sig, <I_NoNr/sig>
!  >> %obv = <I_obv> / (<I_obv> + <I_rev>)
!
! le fichier '_hklf5.hkl' est genere en tenant compte des instructions contenues dans :
! "Refinement of obverse/reverse twins"
! R. Herbst-Irmer and G. Sheldrick
! Acta Cryst. (2002) B58, 477-481
!-----------------------------------------------------------
!
 implicit none
  integer                           :: i, j, i1, code
  character (len=256)               :: HKLF5_file
  integer                           :: hkl_obv, hkl_rev
  integer                           :: n, n_obv, n_rev, n_nonr, n_oar
  real                              :: Imoy_obv,      Imoy_rev,      Imoy_nonr,      Imoy_oar
  real                              :: I_sig_moy_obv, I_sig_moy_rev, I_sig_moy_nonr, I_sig_moy_oar
  real                              :: fraction_rev, fraction_rev2


  if(debug_proc%level_2)  call write_debug_proc_level(2, "twin_obverse_reverse")


  if (.not. keyword_FILE) then
   call write_info('')
   call write_info(' Please input hkl file !')
   call write_info('')
   return
  end if

  n      = 0
  n_obv  = 0
  n_rev  = 0
  n_nonr = 0
  Imoy_obv  = 0.
  Imoy_rev  = 0.
  Imoy_nonr = 0.
  I_sig_moy_obv  = 0.
  I_sig_moy_rev  = 0.
  I_sig_moy_nonr = 0.


  i1 = INDEX(HKL_file%NAME, '.')
  HKLF5_file     = TRIM(HKL_file%NAME(1:i1-1))//'_hklf5.hkl'

  open (unit=HKLF5_unit, file=trim(HKLF5_file))

  ! lecture du fichier
  do i=1, n_ref

    !READ(unit=HKl_unit,'(3I4,2F8.2,I4,6F8.5)', iostat=ier) h,k,l, F2, sig_F2,code,(cos_dir(i),i=1,6)
    !if (ier < 0) THEN ! fin du fichier atteinte
    ! exit
    !end if
    !if (h==0 .and. k==0 .and. l==0)  exit

    ! analyse des reflections
    hkl_obv = -h(i)+k(i)+l(i)
    hkl_rev =  h(i)-k(i)+l(i)

    !if (3*int(hkl_obv/3)==hkl_obv)     then
    if(multiple(hkl_obv,3)) then
     n_obv = n_obv+1
     Imoy_obv = Imoy_obv + F2(i)
     I_sig_moy_obv =  I_sig_moy_obv + F2(i)/sig_F2(i)
     code = 1   ! domaine 1 exclusivement
    end if

    !if (3*int(hkl_rev/3)==hkl_rev) then
    if(multiple(hkl_rev, 3)) then
     n_rev = n_rev+1
     Imoy_rev = Imoy_rev + F2(i)
     I_sig_moy_rev =  I_sig_moy_rev + F2(i)/sig_F2(i)
     code = 2   ! domaine 2 uniquement (à exclure)
    end if

    !if (3*int(hkl_obv/3)==hkl_obv .and. 3*int(hkl_rev/3)==hkl_rev) then
    if(multiple(hkl_obv,3) .and. multiple(hkl_rev,3)) then
     n_oar = n_oar + 1
     Imoy_oar = Imoy_oar + F2(i)
     I_sig_moy_oar =  I_sig_moy_oar + F2(i)/sig_F2(i)
     code = -2   ! domaine 1 + domaine 2
    endif

    !if (3*int(hkl_obv/3)/=hkl_obv .and. 3*int(hkl_rev/3)/=hkl_rev) then
    if(.not. multiple(hkl_obv, 3) .and. .not. multiple(hkl_rev,3)) then
     n_nonr = n_nonr + 1
     Imoy_nonr = Imoy_nonr + F2(i)
     I_sig_moy_nonr =  I_sig_moy_nonr + F2(i)/sig_F2(i)
     code = 0   ! reflection absente: a exclure
    endif


    n = n + 1

    if (code /=0) then
     if (code ==-2) then   ! contributions des 2 domaines
      if (twin_matrix ==1) then
       WRITE(HKLF5_unit,'(3I4,2F8.2,I4,6F8.5)') -h(i),-k(i), l(i), F2(i), sig_F2(i), -2,(cos_dir(i,j),j=1,6)
      else
       WRITE(HKLF5_unit,'(3I4,2F8.2,I4,6F8.5)') -k(i),-h(i),-l(i), F2(i), sig_F2(i), -2,(cos_dir(i,j),j=j,6)
      endif
      WRITE(HKLF5_unit,'(3I4,2F8.2,I4,6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i),1,(cos_dir(i,j),j=1,6)

     else
      WRITE(HKLF5_unit,'(3I4,2F8.2,I4,6F8.5)') h(i),k(i),l(i), F2(i), sig_F2(i),code,(cos_dir(i,j),j=1,6)
     end if
    end if
  end do

  !close (unit=1)
  close (unit=HKLF5_unit)

  if (n_obv  /=0)  Imoy_obv  = Imoy_obv  / n_obv
  if (n_rev  /=0)  Imoy_rev  = Imoy_rev  / n_rev
  if (n_nonr /=0)  Imoy_nonr = Imoy_nonr / n_nonr
  if (n_oar  /=0)  Imoy_oar  = Imoy_oar  / n_oar

  if (n_obv  /=0)  I_sig_moy_obv  = I_sig_moy_obv  / n_obv
  if (n_rev  /=0)  I_sig_moy_rev  = I_sig_moy_rev  / n_rev
  if (n_nonr /=0)  I_sig_moy_nonr = I_sig_moy_nonr / n_nonr
  if (n_oar  /=0)  I_sig_moy_oar  = I_sig_moy_oar  / n_oar

  fraction_rev   =  Imoy_rev / (Imoy_obv+Imoy_rev)
  fraction_rev2  =  I_sig_moy_rev / (I_sig_moy_obv+I_sig_moy_rev)


  call write_info(' ')
  WRITE(message_text,'(a,I8,a)')  '   ANALYSIS OF THE ',n_ref, ' reflections:'
  call write_info(TRIM(message_text))
  call write_info(' ')
  call write_info(' ')
  call write_info('  >> Obverse reflections (-h +k +l = 3n) :  ')
  WRITE(message_text,'(a,I8)')    '       . number of reflections : ', n_obv
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')  '       . mean intensity       : ', Imoy_obv
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')  '       . mean intensity/sigma : ', I_sig_moy_obv
  call write_info(TRIM(message_text))

  call write_info(' ')
  call write_info(' ')
  call write_info('  >> Reverse reflections (+h -k +l = 3n) :  ')
  WRITE(message_text,'(a,I8)')    '       . number of reflections : ', n_rev
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')  '       . mean intensity       : ', Imoy_rev
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')  '       . mean intensity/sigma : ', I_sig_moy_rev
  call write_info(TRIM(message_text))

  call write_info(' ')
  call write_info(' ')
  call write_info('  >> Obverse and reverse reflections :  ')
  WRITE(message_text,'(a,I8)')    '       . number of reflections : ', n_oar
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')  '       . mean intensity       : ', Imoy_oar
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')  '       . mean intensity/sigma : ', I_sig_moy_oar
  call write_info(TRIM(message_text))

  call write_info(' ')
  call write_info(' ')
  call write_info('  >> Neither obverse nor reverse reflections :  ')
  WRITE(message_text,'(a,I8)')    '       . number of reflections : ', n_nonr
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')  '       . mean intensity       : ', Imoy_nonr
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')  '       . mean intensity/sigma : ', I_sig_moy_nonr
  call write_info(TRIM(message_text))


  call write_info(' ')
  call write_info(' ')
  call write_info('  >> Estimated fractionnal contribution of reverse component: ')
  WRITE(message_text,'(a,F6.3)')  '       . from mean intensities:            ',fraction_rev
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F6.3)')  '       . from mean intensity/sigma values: ',fraction_rev2
  call write_info(TRIM(message_text))

  call write_info(' ')
  call write_info(' ')
  WRITE(message_text,'(2a)')      '   >>>>> ',trim(HKLF5_file)//' file has been created (HKLF5 format)'
  call write_info(TRIM(message_text))


 return
end subroutine twin_obverse_reverse

