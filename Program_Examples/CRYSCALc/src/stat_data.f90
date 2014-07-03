!     Last change:  TR   19 Sep 2006    4:01 pm
!--------------------------------------------------------------------------
 subroutine test_ratio_I_sig()
 USE cryscalc_module, ONLY : message_text, debug_proc 
 use hkl_module,      only : F2, sig_F2, n_ref
 USE IO_module,       ONLY : write_info
  implicit none
   ! local variables
   integer                                     :: nb_I_neg, nb_sig0, nb_sig1, nb_sig2, nb_sig3, nb_sig4, nb_sig5, nb_sig10
   real                                        :: F2_max, ratio
   integer, dimension(11)                      :: nb_ratio
   integer                                     :: i, n

   if(debug_proc%level_2)  call write_debug_proc_level(2, "test_ratio_I_sig")
  
   ! statistique sur les I/sig
   nb_I_neg = 0
   nb_sig0  = 0
   nb_sig1  = 0
   nb_sig2  = 0
   nb_sig3  = 0
   nb_sig4  = 0
   nb_sig5  = 0
   nb_sig10 = 0

   do i =1, n_ref

    if (F2(i) < 0.) then
      nb_I_neg = nb_I_neg + 1

    elseif(F2(i)/sig_F2(i) >= 10.) then
      nb_sig10 = nb_sig10 + 1

    elseif(F2(i)/sig_F2(i) >=  5.) then
      nb_sig5 = nb_sig5 + 1

    elseif(F2(i)/sig_F2(i) >=  4.) then
      nb_sig4 = nb_sig4 + 1

    elseif(F2(i)/sig_F2(i) >=  3.) then
      nb_sig3 = nb_sig3 + 1

    elseif(F2(i)/sig_F2(i) >=  2.) then
      nb_sig2 = nb_sig2 + 1

    elseif(F2(i)/sig_F2(i) >=  1.) then
      nb_sig1 = nb_sig1 + 1

    else
      nb_sig0 = nb_sig0 + 1
    endif

   end do

   write(message_text,'(a)')              '   ---- Statistics on F2/sig: '
    call write_info(TRIM(message_text))
   write(message_text,'(a,I8,a,F7.3,a)')  '         . F2/sig >= 10.: ',nb_sig10     ,  '  (', 100*real(nb_sig10)/REAL(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,'(a,I8,a,F7.3,a)')  '         . F2/sig >=  5.: ',nb_sig5      ,  '  (', 100*real(nb_sig5)/REAL(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,'(a,I8,a,F7.3,a)')  '         . F2/sig >=  4.: ',nb_sig4      ,  '  (', 100*real(nb_sig4)/REAL(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,'(a,I8,a,F7.3,a)')  '         . F2/sig >=  3.: ',nb_sig3      ,  '  (', 100*real(nb_sig3)/REAL(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,'(a,I8,a,F7.3,a)')  '         . F2/sig >=  2.: ',nb_sig2      ,  '  (', 100*real(nb_sig2)/REAL(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,'(a,I8,a,F7.3,a)')  '         . F2/sig >=  1.: ',nb_sig1      ,  '  (', 100*real(nb_sig1)/REAL(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,'(a,I8,a,F7.3,a)')  '         . F2/sig  <  1.: ',nb_sig0      ,  '  (', 100*real(nb_sig0)/REAL(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,'(a,I8,a,F7.3,a)')  '         .     F2  <  0.: ',nb_I_neg     ,  '  (', 100*real(nb_I_neg)/REAL(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,'(a)')              ' '
    call write_info(TRIM(message_text))


   ! statistique sur les F2/F2_max
   nb_ratio(:)=0
   F2_max = maxval(F2(:))

   do i=1, n_ref
    ratio = F2(i) / F2_max

    if (ratio > 0.9) then
     nb_ratio(11) = nb_ratio(11) + 1

    elseif(ratio > 0.8) then
     nb_ratio(10) = nb_ratio(10) + 1

    elseif(ratio > 0.7) then
     nb_ratio(9) = nb_ratio(9) + 1

    elseif(ratio > 0.6) then
     nb_ratio(8) = nb_ratio(8) + 1

    elseif(ratio > 0.5) then
     nb_ratio(7) = nb_ratio(7) + 1

    elseif(ratio > 0.4) then
     nb_ratio(6) = nb_ratio(6) + 1

    elseif(ratio > 0.3) then
     nb_ratio(5) = nb_ratio(5) + 1

    elseif(ratio > 0.2) then
     nb_ratio(4) = nb_ratio(4) + 1

    elseif(ratio > 0.1) then
     nb_ratio(3) = nb_ratio(3) + 1

    elseif(ratio > 0.0) then
     nb_ratio(2) = nb_ratio(2) + 1

    else
     nb_ratio(1) = nb_ratio(1) + 1

    endif

   end do

   write(message_text,'(a)')         '   ---- Statistics on F2/F2_max: '
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 1.0 >=  F2/F2_max > 0.9: ',nb_ratio(11) , ' (', &
                                                              100*real(nb_ratio(11))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.9 >=  F2/F2_max > 0.8: ',nb_ratio(10) , ' (', &
                                                              100*real(nb_ratio(10))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.8 >=  F2/F2_max > 0.7: ',nb_ratio(9)  , ' (', &
                                                              100*real(nb_ratio(9))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.7 >=  F2/F2_max > 0.6: ',nb_ratio(8)  , ' (', &
                                                              100*real(nb_ratio(8))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.6 >=  F2/F2_max > 0.5: ',nb_ratio(7)  , ' (', &
                                                              100*real(nb_ratio(7))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.5 >=  F2/F2_max > 0.4: ',nb_ratio(6)  , ' (', &
                                                              100*real(nb_ratio(6))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.4 >=  F2/F2_max > 0.3: ',nb_ratio(5)  , ' (', &
                                                              100*real(nb_ratio(5))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.3 >=  F2/F2_max > 0.2: ',nb_ratio(4)  , ' (', &
                                                              100*real(nb_ratio(4))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.2 >=  F2/F2_max > 0.1: ',nb_ratio(3)  , ' (', &
                                                              100*real(nb_ratio(3))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.1 >=  F2/F2_max > 0.0: ',nb_ratio(2)  , ' (', &
                                                              100*real(nb_ratio(2))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a,I8,a,F7.3,a)')      '         . 0.0 <=  F2/F2_max      : ',nb_ratio(1)  , ' (', &
                                                              100*real(nb_ratio(1))/real(n_ref), ' %)'
    call write_info(TRIM(message_text))
   write(message_text,*)             ' '
    call write_info(TRIM(message_text))


   n = nb_ratio(11)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.9: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(10)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.8: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(9)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.7: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(8)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.6: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(7)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.5: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(6)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.4: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(5)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.3: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(4)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.2: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(3)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.1: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = n + nb_ratio(2)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max > 0.0: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   n = nb_ratio(1)
   write(message_text,fmt='(a,F7.3,a)')  '         . F2/F2_max < 0.0: ',100*real(n)/real(n_ref)  ,' %'
    call write_info(TRIM(message_text))
   write(message_text,fmt='(a)')         ' '
    call write_info(TRIM(message_text))




  return
 end subroutine test_ratio_I_sig
!--------------------------------------------------------------------------


   ! statistiques sur la repartition en fonction de sintheta/lambda

subroutine  stat_repartition_sintheta_lambda()
 USE cryscalc_module, ONLY   : message_text, debug_proc
 use hkl_module,      ONLY   : F2, sintheta_lambda, n_ref
 USE IO_module,       ONLY   : write_info
 implicit none
  ! local variables
  real                                  :: sintheta_lambda_max
  integer, dimension(15)                :: nb_shell
  real,    dimension(15)                :: F2_mean_shell, F2_min_shell, F2_max_shell
  integer                               :: i
  character(len=128)                    :: fmt_, range_

  if(debug_proc%level_2)  call write_debug_proc_level(2, "stat_repartition_sintheta_lambda")
  
  sintheta_lambda_max = maxval(sintheta_lambda(1:n_ref))

  nb_shell(1:15)=0
  F2_mean_shell(1:15) = 0.
  F2_min_shell(1:15)  = 1.E09
  F2_max_shell(1:15)  = -1.E09

  if(debug_proc%level_2)  call write_debug_proc_level(2, "statistic_shell_stl")
  do i=1, n_ref
   if (sintheta_lambda(i) < 0.1) then
    call statistic_shell(1, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 0.2) then
    call statistic_shell(2, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 0.3) then
    call statistic_shell(3, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 0.4) then
    call statistic_shell(4, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 0.5) then
    call statistic_shell(5, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 0.6) then
    call statistic_shell(6, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 0.7) then
    call statistic_shell(7, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 0.8) then
    call statistic_shell(8, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 0.9) then
    call statistic_shell(9, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 1.0) then
    call statistic_shell(10, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 1.1) then
    call statistic_shell(11, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 1.2) then
    call statistic_shell(12, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 1.3) then
    call statistic_shell(13, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 1.4) then
    call statistic_shell(14, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)

   elseif(sintheta_lambda(i) < 1.5) then
    call statistic_shell(15, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell, 15)
   endif
  end do

  do i=1, 15
   F2_mean_shell(i) =   F2_mean_shell(i)/nb_shell(i)
  end do
   call write_info(' ')
   call write_info('   ---- Number of reflections in sintheta/lambda ranges:')
   call write_info(' ')   
   call write_info('  sinTheta/lambda range       reflections            <F2>        F2_min        F2_max')

   fmt_ = '(a,3x,I8,a,F7.3,a,3(5x,F9.2))'
  do i = 1, 15
   if (nb_shell(i) ==0) cycle
   WRITE(range_,'(8x,F4.1,a,F4.1)') REAL(i-1)/10, ' -- ',REAL(i)/10
   WRITE(message_text, fmt_)  TRIM(range_),nb_shell(i),  ' (',100*real(nb_shell(i))/real(n_ref),' %)',   &
                             F2_mean_shell(i) , F2_min_shell(i) , F2_max_shell(i)
   call write_info(TRIM(message_text))
  end do
  call write_info(' ')

 return

end subroutine  stat_repartition_sintheta_lambda

!-----------------------------------------------------------------------------------------
   ! statistiques sur la repartition en fonction de d_hkl

subroutine  stat_repartition_d_hkl()
 USE cryscalc_module, ONLY   : message_text, debug_proc
 use hkl_module,      ONLY   :  F2, sig_F2, d_hkl, n_ref
 USE IO_module,       ONLY   : write_info
 implicit none
  ! local variables
  real                                  :: d_hkl_max, d_hkl_min
  integer, dimension(10)                :: nb_shell
  REAL,    DIMENSION(10)                :: d_shell_min, d_shell_max
  real,    dimension(10)                :: F2_mean_shell, F2_min_shell, F2_max_shell
  integer                               :: i
  character(len=128)                    :: fmt_, range_

  if(debug_proc%level_2)  call write_debug_proc_level(2, "stat_repartition_d_hkl")
  
  d_hkl_max = maxval(d_hkl(1:n_ref))
  d_hkl_min = maxval(d_hkl(1:n_ref))

  nb_shell(1:10)=0
  F2_mean_shell(1:10) = 0.
  F2_min_shell(1:10)  = 1.E09
  F2_max_shell(1:10)  = -1.E09



  do i=1, n_ref
   if (i <= 0.1*n_ref) then
    call statistic_d_shell(1, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   elseif (i <= 0.2*n_ref) then
    call statistic_d_shell(2, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   elseif (i <= 0.3*n_ref) then
    call statistic_d_shell(3, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   elseif (i <= 0.4*n_ref) then
    call statistic_d_shell(4, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   elseif (i <= 0.5*n_ref) then
    call statistic_d_shell(5, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   elseif (i <= 0.6*n_ref) then
    call statistic_d_shell(6, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   elseif (i <= 0.7*n_ref) then
    call statistic_d_shell(7, nb_shell,  d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   elseif (i <= 0.8*n_ref) then
    call statistic_d_shell(8, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   elseif (i <= 0.9*n_ref) then
    call statistic_d_shell(9, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
   else
    call statistic_d_shell(10, nb_shell, d_hkl(i),  F2(i),F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)

   end if

  end do

  do i=1, 10
   F2_mean_shell(i) =   F2_mean_shell(i)/nb_shell(i)
  end do

   call write_info(' ')
   call write_info('   ---- Number of reflections in d_hkl ranges:')
   call write_info(' ')
   call write_info('        d_hkl range        reflections           <F2>        F2_min        F2_max')

   fmt_ = '(a,3x,I8,6x,3(5x,F9.2))'
  do i = 1, 10
   if (nb_shell(i) ==0) cycle

   WRITE(range_,'(8x,F5.2,a,F5.2)') d_shell_min(i), ' -- ',d_shell_max(i)
   WRITE(message_text, fmt_)  TRIM(range_),nb_shell(i),   F2_mean_shell(i) , F2_min_shell(i) , F2_max_shell(i)
   call write_info(TRIM(message_text))
  end do
  call write_info('')

 return

end subroutine  stat_repartition_d_hkl
!-------------------------------------------------------------------------------------

   ! statistiques sur la repartition en fonction de theta

subroutine  stat_repartition_theta ()
 USE cryscalc_module, ONLY   : message_text, debug_proc
 USE hkl_module,      ONLY   : F2, sig_F2, theta_hkl, n_ref
 USE IO_module,       ONLY   : write_info
 implicit none
  ! local variables
  real                                  :: theta_max, theta_min
  integer, dimension(12)                :: nb_shell
  REAL,    DIMENSION(12)                :: theta_shell
  real,    dimension(12)                :: F2_mean_shell, F2_min_shell, F2_max_shell
  integer                               :: i
  character(len=128)                    :: fmt_, range_

  if(debug_proc%level_2)  call write_debug_proc_level(2, "stat_repartition_theta")
  
  theta_min = minval(theta_hkl(1:n_ref))
  theta_max = maxval(theta_hkl(1:n_ref))

  nb_shell(1:12)=0
  F2_mean_shell(1:12) = 0.
  F2_min_shell(1:12)  = 1.E09
  F2_max_shell(1:12)  = -1.E09

  theta_shell(1) =   5.
  theta_shell(2) =  10.
  theta_shell(3) =  15.
  theta_shell(4) =  20.
  theta_shell(5) =  25.
  theta_shell(6) =  30.
  theta_shell(7) =  35.
  theta_shell(8) =  40.
  theta_shell(9) =  45.
  theta_shell(10) = 50.
  theta_shell(11) = 60.

  if(debug_proc%level_2)  call write_debug_proc_level(2, "statistic_shell_theta")
  do i=1, n_ref
   if (theta_hkl(i) < theta_shell(1)) then
    call statistic_shell(1, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(2)) then
    call statistic_shell(2, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(3)) then
    call statistic_shell(3, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(4)) then
    call statistic_shell(4, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(5)) then
    call statistic_shell(5, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(6)) then
    call statistic_shell(6, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(7)) then
    call statistic_shell(7, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(8)) then
    call statistic_shell(8, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(9)) then
    call statistic_shell(9, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(10)) then
    call statistic_shell(10, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   elseif(theta_hkl(i) < theta_shell(11)) then
    call statistic_shell(11, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   else
    call statistic_shell(12, F2(i), nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,12)

   endif
  end do





  do i=1, 12
   F2_mean_shell(i) =   F2_mean_shell(i)/nb_shell(i)
  end do

   call write_info('')
   call write_info('   ---- Number of reflections in theta ranges:')
   call write_info('')
   call write_info('        theta range        reflections           <F2>        F2_min        F2_max')

   fmt_ = '(a,3x,I8,6x,3(5x,F9.2))'
  do i = 1, 12
   if (nb_shell(i) ==0) cycle
   if (i==1) then
    WRITE(range_, '(8x,F5.2,a,F5.2)')  0.00, ' -- ',theta_shell(i)
   ELSEIF(i==12) then
    WRITE(range_, '(8x,F5.2,a )')      theta_shell(i-1), ' -- '
   else
    WRITE(range_,'(8x,F5.2,a,F5.2)')   theta_shell(i-1), ' -- ',theta_shell(i)
   endif
   WRITE(message_text, fmt_)  TRIM(range_),nb_shell(i),   F2_mean_shell(i) , F2_min_shell(i) , F2_max_shell(i)
   call write_info(TRIM(message_text))
  end do
  call write_info('')

 return

end subroutine  stat_repartition_theta



!-----------------------------------------------------------------------------------------
subroutine  statistic_d_shell(num_shell, nb_shell, d_hkl_, F2_,F2_mean_shell, F2_min_shell, F2_max_shell, d_shell_min, d_shell_max)
 use cryscalc_module, only : debug_proc
 implicit none
  INTEGER, INTENT(IN)                    :: num_shell
  REAL,    INTENT(IN)                    :: d_hkl_
  REAL,    INTENT(IN)                    :: F2_
  integer, intent(inout), dimension(10)  :: nb_shell
  REAL,    INTENT(INOUT), DIMENSION(10)  :: F2_mean_shell, F2_min_shell, F2_max_shell
  REAL,    INTENT(INOUT), DIMENSION(10)  :: d_shell_min, d_shell_max

  if(debug_proc%level_2)  call write_debug_proc_level(2, "statistic_d_shell")

  nb_shell(num_shell) = nb_shell(num_shell) + 1
  F2_mean_shell(num_shell) = F2_mean_shell(num_shell) + F2_
  if (nb_shell(num_shell) == 1) then
   d_shell_min(num_shell) = d_hkl_
  else
   d_shell_max(num_shell) = d_hkl_
  endif

  if (F2_ < F2_min_shell(num_shell))    F2_min_shell(num_shell) = F2_
  if (F2_ > F2_max_shell(num_shell))    F2_max_shell(num_shell) = F2_

 RETURN
end subroutine statistic_d_shell


!-------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------



subroutine statistic_shell(num_shell, F2_, nb_shell, F2_mean_shell, F2_min_shell, F2_max_shell,dim_shell)
 use cryscalc_module, only : debug_proc
 implicit none
  integer, intent(in)                          :: num_shell
  real,    intent(in)                          :: F2_
  INTEGER, INTENT(IN)                          :: dim_shell
  integer, intent(inout), dimension(dim_shell) :: nb_shell
  real,    intent(inout), dimension(dim_shell) :: F2_mean_shell, F2_min_shell, F2_max_shell

  !if(debug_proc%level_2)  call write_debug_proc_level(2, "statistic_shell")

   nb_shell(num_shell) = nb_shell(num_shell) + 1
   F2_mean_shell(num_shell) = F2_mean_shell(num_shell) + F2_

   if (F2_ < F2_min_shell(num_shell))    F2_min_shell(num_shell) = F2_

   if (F2_ > F2_max_shell(num_shell))    F2_max_shell(num_shell) = F2_


 return
end subroutine statistic_shell

!--------------------------------------------------------------------------------------------------

subroutine statistics_on_E2()
 USE IO_module,       ONLY : write_info
 USE cryscalc_module, ONLY: message_text, debug_proc
 use  HKL_module
 implicit none
  INTEGER   :: i
  real      :: n_1, n_2, n_3

  if(debug_proc%level_2)  call write_debug_proc_level(2, "statistics_on_E2")

 write(message_text,'(a)')              '   ---- Statistics on E and E2: '
  call write_info(TRIM(message_text))

  F2_mean = SUM(F2(1:n_ref)) / n_ref
  E2(1:n_ref) = F2(1:n_ref) / F2_mean
  E2_mean = SUM(E2(1:n_ref)) / n_ref
  E2m1_mean = SUM(ABS(E2(1:n_ref)-1))/n_ref

  E(1:n_ref)  = SQRT(E2(1:n_ref))
  E_mean = SUM(E(1:n_ref)) / n_ref

  E3(1:n_ref) = E(1:n_ref)**3
  E4(1:n_ref) = E(1:n_ref)**4
  E5(1:n_ref) = E(1:n_ref)**5
  E6(1:n_ref) = E(1:n_ref)**6

  E3_mean = SUM(E3(1:n_ref))/ n_ref
  E4_mean = SUM(E4(1:n_ref))/ n_ref
  E5_mean = SUM(E5(1:n_ref))/ n_ref
  E6_mean = SUM(E6(1:n_ref))/ n_ref

  n_1 = 0.
  n_2 = 0.
  n_3 = 0.
  do i=1, n_ref
   IF(ABS(E(i)) > 3.0) then
    n_3 = n_3 +1
    cycle
   endif
   IF(ABS(E(i)) > 2.0) then
    n_2 = n_2 +1
    cycle
   endif
   IF(ABS(E(i)) > 1.0) then
    n_1 = n_1 +1
    cycle
   endif
  end do
  n_1 = (n_1 /n_ref) * 100
  n_2 = (n_2 /n_ref) * 100
  n_3 = (n_3 /n_ref) * 100

  WRITE(message_text,*) '              centro     non-centro           exp.'
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F6.3,10x,F6.3,10xF6.3)') '  <ABS(E)**2>:  ', 1.000, 1.000, E2_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F6.3,10x,F6.3,10xF6.3)') '  <ABS(E)**3>:  ', 1.000, 1.000, E3_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F6.3,10x,F6.3,10xF6.3)') '  <ABS(E)**4>:  ', 1.000, 1.000, E4_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F6.3,10x,F6.3,10xF6.3)') '  <ABS(E)**5>:  ', 1.000, 1.000, E5_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F6.3,10x,F6.3,10xF6.3)') '  <ABS(E)**6>:  ', 1.000, 1.000, E6_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F6.3,10x,F6.3,10xF6.3)') '  <ABS(E2-1)>:  ', 0.968, 0.736, E2m1_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F6.3,10x,F6.3,10xF6.3)') '     <ABS(E)>:  ', 0.798, 0.886, E_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F5.2,10x,F5.2,10xF5.2)') '  %ABS(E)>1.0:  ', 31.7, 36.8 , n_1
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F5.2,10x,F5.2,10xF5.2)') '  %ABS(E)>2.0:  ',  4.6,  1.8 , n_2
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F5.2,10x,F5.2,10xF5.2)') '  %ABS(E)>3.0:  ',  0.3,  0.02, n_3
  call write_info(TRIM(message_text))

 RETURN
end subroutine statistics_on_E2
