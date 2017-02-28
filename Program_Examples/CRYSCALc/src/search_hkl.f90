!     Last change:  TR   26 Apr 2007    2:19 pm
subroutine search_exti()
 USE IO_module,       ONLY:  write_info
 USE cryscalc_module, ONLY : keyword_file, message_text, debug_proc
 USE hkl_module
 implicit none
  INTEGER               :: i

 if(debug_proc%level_2)  call write_debug_proc_level(2, "search_exti")

 IF(.NOT. keyword_file) then
  call write_info('')
  call write_info('  !! FILE keyword mandatory for the search extinctions procedure !!')
  call write_info('')
  return
 endif

 call write_info('')
 call write_info('  ... Search extinction rules: ')
 if(.not.HKL_list%ALL) then
  WRITE(message_text, '(a,F4.0,a)') '       . criteria for I/sig  :          ', n_sig
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,F7.3)')   '       . criteria for <I_odd>/<I_even>: ', ratio_criteria
  call write_info(TRIM(message_text))
 end if
 call write_info(' ')

 call write_info('')
 WRITE(message_text, '(38a1,a)') (' ',i=1,38), 'all   odd     <I_odd>    <I_even>  <I_odd>/<I_even>'
 call write_info(TRIM(message_text))
 call write_info('')

 do i=1, HKL_rule_nb
  call search_HKL_exti(i)
  call write_hkl_exti(i, HKL_rule(i))
  if(i==3 .or. i==6  .or. i==9  .or. i==12 .or. i==15  .or. i==20 .or. i==23   &
          .or. i==26 .or. i==27 .or. i==30) call write_info('')
 end do



 RETURN
end subroutine search_exti
!------------------------------------------------------

subroutine search_HKL_exti(rule_number)
 USE macros_module,   ONLY : multiple
 USE HKL_module
 USE cryscalc_module, only : debug_proc
 implicit none
  INTEGER, INTENT(IN)          :: rule_number
  INTEGER                      :: i

  if(debug_proc%level_2)  call write_debug_proc_level(2, "search_HKl_exti")

  n_pair    = 0
  n_impair  = 0
  F2_pair   = 0.
  F2_impair = 0.
  HKL_search_ok = .false.

 select case (ABS(rule_number))
 case (1)      ! 1. h00   h = 2n+1 ----------------------
     do i=1, n_ref
      HKL_search_ok(i)=.false.
      IF(.NOT. HKL_list%ALL) then
       if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
       IF(F2(i) < 0.01) cycle
      endif
      IF(h(i)/=0 .AND. k(i)==0 .AND. l(i)==0) then
       !IF(h(i)==2*INT(h(i)/2)) then
       if(multiple(h(i),2)) then
        n_pair = n_pair + 1
        F2_pair = F2_pair + F2(i)
        if(rule_number < 0) HKL_search_ok(i) = .true.
       else
        n_impair = n_impair + 1
        F2_impair = F2_impair + F2(i)
        if(rule_number > 0) HKL_search_ok(i) = .true.
       endif
      endif
     end do


 CASE (2)   ! 2. 0k0 k=2n+1 ------------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==0 .AND. k(i)/=0 .AND. l(i)==0) then
   !IF(k(i)==2*INT(k(i)/2)) then
   if(multiple(k(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE (3)   ! 3. 00l  l=2n+1-------------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==0 .AND. k(i)==0 .AND. l(i)/=0) then
   !IF(l(i)==2*INT(l(i)/2)) then
   if(multiple(l(i),2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do

  CASE(4)    ! 4. 0kl k=2n+1  ----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==0 .AND. k(i)/=0 .AND. l(i)/=0) then
   !IF(k(i)==2*INT(k(i)/2)) then
   if(multiple(k(i),2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(5)   ! 5. 0kl l=2n +1 -----------------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==0 .AND. k(i)/=0 .AND. l(i)/=0) then
   !IF(l(i)==2*INT(l(i)/2)) then
   if(multiple(l(i),2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(6)   ! 6. 0kl k+l=2n+1  -----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==0 .AND. k(i)/=0 .AND. l(i)/=0) then
   !IF(k(i)+l(i)==2*INT((k(i)+l(i))/2)) then
   if(multiple(k(i)+l(i),2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(7)     ! 7. h0l h=2n+1  ----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)==0 .AND. l(i)/=0) then
   !IF(h(i)==2*INT(h(i)/2)) then
   if(multiple(h(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(8)   ! 8. h0l l=2n+1 -----------------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)==0 .AND. l(i)/=0) then
   !IF(l(i)==2*INT(l(i)/2)) then
   if(multiple(l(i),2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(9)  ! 9. h0l h+l=2n+1  -----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)==0 .AND. l(i)/=0) then
   !IF(h(i)+l(i)==2*INT((h(i)+l(i))/2)) then
   if(multiple(h(i)+l(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do



 CASE(10)  ! 10. hk0 h=2n+1  ----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)/=0 .AND. l(i)==0) then
   !IF(h(i)==2*INT(h(i)/2)) then
   if(multiple(h(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do

 CASE(11)   ! 11. hk0 k=2n+1 -----------------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)/=0 .AND. l(i)==0) then
   !IF(k(i)==2*INT(k(i)/2)) then
   if(multiple(k(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(12)   ! 12. hk0 h+k=2n+1  -----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)/=0 .AND. l(i)==0) then
   !IF(h(i)+k(i)==2*INT((h(i)+k(i))/2)) then
   if(multiple(h(i)+k(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 !------------------ new . dec 10 : hhl, hkk, hkh
 CASE(13)   ! 13. hhl h+l=2n+1  -----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==k(i) .AND. l(i)/=0) then
   !IF(h(i)+l(i)==2*INT((h(i)+l(i))/2)) then
   if (multiple(h(i)+l(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(14)   ! 14. hkk k+h=2n+1  -----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(k(i)==l(i) .AND. h(i)/=0) then
   !IF(k(i)+h(i)==2*INT((k(i)+h(i))/2)) then
   if(multiple(k(i)+h(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(15)   ! 15. hkh h+k=2n+1  -----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==l(i) .AND. k(i)/=0) then
   !IF(h(i)+k(i)==2*INT((h(i)+k(i))/2)) then
   if(multiple(h(i)+k(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do

 !------------------


 CASE(16)  ! 16. hkl k+l=2n+1      A   -----------------------
 do i=1, n_ref
  HKL_search_ok = .false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  if(k(i)/=0 .and. l(i)/=0) then
   !IF(k(i)+l(i)==2*INT((k(i)+l(i))/2)) then
   if(multiple(k(i)+l(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  end if
 end do


 CASE(17)  ! 17. hkl h+l=2n+1      B   -----------------------
 do i=1, n_ref
  HKL_search_ok = .false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  if(h(i)/=0 .and. l(i)/=0) then
   !IF(h(i)+l(i)==2*INT((h(i)+l(i))/2)) then
   if(multiple(h(i)+l(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  end if
 end do


 CASE(18)  ! 18. hkl h+k=2n+1      C   -----------------------
 do i=1, n_ref
  HKL_search_ok = .false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  if(h(i)/=0 .and. k(i)/=0) then
   !IF(h(i)+k(i)==2*INT((h(i)+k(i))/2)) then
   if(multiple(h(i)+k(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  end if
 end do



 CASE(19)  ! 19. hkl h+k+l=2n+1  I   -----------------------
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)/=0 .AND. l(i)/=0) then
   !IF(h(i)+k(i)+l(i)==2*INT((h(i)+k(i)+l(i))/2)) then
   if(multiple(h(i)+k(i)+l(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(20)  ! 20. hkl not all odd/even F
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)/=0 .AND. l(i)/=0) then
   !IF(h(i)==2*INT(h(i)/2) .and. k(i)==2*INT(k(i)/2) .and. l(i)==2*INT(l(i)/2)) then
   if(multiple(h(i), 2) .and. multiple(k(i),2) .and. multiple(l(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(21)   ! 21. h00     h=4n+1 41 .  .
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)==0 .AND. l(i)==0) then
   !IF(h(i)==4*INT(h(i)/4) ) then
   if(multiple(h(i), 4)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(22)   ! 22. 0k0     k=4n+1 .  41 .
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==0 .AND. k(i)/=0 .AND. l(i)==0) then
   !IF(k(i)==4*INT(k(i)/4) ) then
   if(multiple(k(i), 4)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do

 CASE(23)  ! 23. 00l     00l     l=4n+1 .  . 41
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==0 .AND. k(i)==0 .AND. l(i)/=0) then
   !IF(l(i)==4*INT(l(i)/4) ) then
   if(multiple(l(i), 4)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 case(24)    ! 24. 0kl   k+l/=4n+1 d  .  .
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)==0 .AND. k(i)/=0 .AND. l(i)/=0) then
   !IF(k(i)+l(i)==4*INT((k(i)+l(i))/4) ) then
   if(multiple(k(i)+l(i), 4)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(25)  ! 25. h0l   h+l=4n+1 .  d  .
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)==0 .AND. l(i)/=0) then
   !IF(h(i)+l(i)==4*INT((h(i)+l(i))/4) ) then
   if(multiple(h(i)+l(i), 4)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(26)  ! 26. hk0   h+k=4n+1 .  .  d
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif
  IF(h(i)/=0 .AND. k(i)/=0 .AND. l(i)==0) then
   !IF(h(i)+k(i)==4*INT((h(i)+k(i))/4) ) then
   if(multiple(h(i)+k(i), 4)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(27)  !27.  h-hl, h0l, 0kl    l=2n+1       '
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif

  IF(h(i)==0 .or. k(i)==0 .or. h(i)+k(i) == 0 ) then
   !IF(l(i)==2*INT(l(i)/2) ) then
   if(multiple(l(i),2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
  endif
 end do


 CASE(28)  !28.  hkl h=2n+1       '
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif

   !IF(h(i)==2*INT(h(i)/2) ) then
   if(multiple(h(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
 end do


 CASE(29)  !29.  hkl k=2n+1       '
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif

   !IF(k(i)==2*INT(k(i)/2) ) then
   if(multiple(k(i), 2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
 end do


 CASE(30)  !30.  hkl l=2n+1       '
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif

   !IF(l(i)==2*INT(l(i)/2) ) then
   if(multiple(l(i),2)) then
    n_pair = n_pair + 1
    F2_pair = F2_pair + F2(i)
    if(rule_number < 0) HKL_search_ok(i) = .true.
   else
    n_impair = n_impair + 1
    F2_impair = F2_impair + F2(i)
    if(rule_number > 0) HKL_search_ok(i) = .true.
   endif
 end do

 CASE(31)  !31.  hkl h=2n+1, k=2n+1       '
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif

  if(rule_number > 0) then ! recherche de h impair et k impair
   !if(h(i)/=2*int(h(i)/2) .and. k(i)/=2*int(k(i)/2)) then
   if(.not. multiple(h(i), 2) .and. .not. multiple(k(i),2)) then
    n_pair  = n_pair + 1
    F2_pair = F2_pair + F2(i)
    HKL_search_ok(i) = .true.
   else
    n_impair  = n_impair + 1
    F2_impair =  F2_impair + F2(i)
   endif

  else                     ! recherche de h pair et k pair
   !if(h(i)==2*int(h(i)/2) .and. k(i)==2*int(k(i)/2)) then
   if(multiple(h(i),2) .and. multiple(k(i),2)) then
    n_pair  = n_pair + 1
    F2_pair = F2_pair + F2(i)
    HKL_search_ok(i) = .true.
   else
    n_impair  = n_impair + 1
    F2_impair =  F2_impair + F2(i)
   endif

  endif

 end do

 CASE(32)  !32.  hkl h=2n+1, l=2n+1       '
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif

  if(rule_number > 0) then ! recherche de h impair et l impair
   !if(h(i)/=2*int(h(i)/2) .and. l(i)/=2*int(l(i)/2)) then
   if(.not. multiple(h(i), 2) .and. .not. multiple(l(i),2)) then
    n_pair  = n_pair + 1
    F2_pair = F2_pair + F2(i)
    HKL_search_ok(i) = .true.
   else
    n_impair  = n_impair + 1
    F2_impair =  F2_impair + F2(i)
   endif

  else                     ! recherche de h pair et l pair
   !if(h(i)==2*int(h(i)/2) .and. l(i)==2*int(l(i)/2)) then
   if(multiple(h(i),2) .and. multiple(l(i), 2)) then
    n_pair  = n_pair + 1
    F2_pair = F2_pair + F2(i)
    HKL_search_ok(i) = .true.
   else
    n_impair  = n_impair + 1
    F2_impair =  F2_impair + F2(i)
   endif

  endif
 end do

 CASE(33)  !33.  hkl k=2n+1, l=2n+1       '
 do i=1, n_ref
  HKL_search_ok(i)=.false.
  IF(.NOT. HKL_list%ALL) then
    if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
    IF(F2(i) < 0.01) cycle
  endif

  if(rule_number > 0) then ! recherche de k impair et l impair
   !if(k(i)/=2*int(k(i)/2) .and. l(i)/=2*int(l(i)/2)) then
   if(.not. multiple(k(i), 2) .and. .not. multiple(l(i), 2)) then
    n_pair  = n_pair + 1
    F2_pair = F2_pair + F2(i)
    HKL_search_ok(i) = .true.
   else
    n_impair  = n_impair + 1
    F2_impair =  F2_impair + F2(i)
   endif

  else                     ! recherche de k pair et l pair
   !if(k(i)==2*int(k(i)/2) .and. l(i)==2*int(l(i)/2)) then
    if(multiple(k(i),2) .and. multiple(l(i),2)) then
    n_pair  = n_pair + 1
    F2_pair = F2_pair + F2(i)
    HKL_search_ok(i) = .true.
   else
    n_impair  = n_impair + 1
    F2_impair =  F2_impair + F2(i)
   endif

  endif
 end do


  end select



 return
end subroutine search_HKL_exti

!------------------------------------------------------------------------------------------------
subroutine write_hkl_exti(exti_num, exti_string)
 USE cryscalc_module, ONLY : message_text, debug_proc
 USE IO_module,       ONLY : write_info
 USE hkl_module,      ONLY : n_pair, n_impair, F2_pair, F2_impair, mean_F2_pair, mean_F2_impair, I_ratio, ratio_criteria
 implicit none
  INTEGER,           INTENT(IN)   :: exti_num
  CHARACTER (LEN=*), INTENT(IN)   :: exti_string
 ! REAL                            :: I_critere

 !I_critere  = 0.03

 if(debug_proc%level_2)  call write_debug_proc_level(2, "write_HKL_exti")

 !IF(n_pair==0 .AND. n_impair==0)       return
 IF(n_pair   /=0)                      mean_F2_pair   = F2_pair/n_pair
 IF(n_impair /=0)                      mean_F2_impair = F2_impair/n_impair
 IF(n_pair   /= 0 .and. n_impair /=0)  I_ratio        = mean_F2_impair / mean_F2_pair

   !'all       odd      <I_odd>    <I_even>  <I_odd>/<I_even>'

 IF(n_pair /=0 .and. n_impair==0) then
  WRITE(message_text, '(1x,I2,a,2I8,12x,1F12.2,18x,a)') exti_num, exti_string, n_pair   , 0       , mean_F2_pair,  '   ???'
  call write_info(TRIM(message_text))
  return
 elseif(n_pair ==0 .and. n_impair/=0) then
  WRITE(message_text, '(1x,I2,a,2I8,1F12.2)') exti_num, exti_string, n_impair , n_impair, mean_F2_impair
  call write_info(TRIM(message_text))
  return
 ELSEIF(n_pair == 0 .AND. n_impair ==0) then
  WRITE(message_text, '(1x,I2,a,I8)') exti_num, exti_string, n_pair+n_impair
  call write_info(TRIM(message_text))
  return
 endif

 IF(I_ratio < ratio_criteria) then
  WRITE(message_text, '(1x,I2,a,2I8,2F12.2,8x,F10.3,a)') exti_num, exti_string, n_pair + n_impair, n_impair, &
                                                         mean_F2_impair, mean_F2_pair, I_ratio, '   <<<'
 else
  WRITE(message_text, '(1x,I2,a,2I8,2F12.2,8x,F10.3,a)') exti_num, exti_string, n_pair + n_impair, n_impair, &
                                                         mean_F2_impair, mean_F2_pair, I_ratio, ''
 endif
 call write_info(TRIM(message_text))

 return

END subroutine write_hkl_exti


!--------------------------------------------------------------------------------------------------------
subroutine search_hkl()
 USE hkl_module
 USE IO_module,                   ONLY : write_info
 USE cryscalc_module,             ONLY : message_text, keyword_FILE, SPG, hkl_format_free, debug_proc
 USE CFML_Reflections_Utilities,  ONLY : HKL_absent,  HKL_equiv

 implicit none
  INTEGER                         :: i, n_req
  REAL                            :: F2_req, sig_F2_req, ratio_req
  REAL                            :: wF2_req, sum_weight
  !REAL, DIMENSION(3)             :: ref_H
  INTEGER,           DIMENSION(3) :: ref_H
  LOGICAL                         :: equiv_H
  CHARACTER (len=256)             :: fmt_1, fmt_2

  if(debug_proc%level_2)  call write_debug_proc_level(2, "search_hkl")

  fmt_1 = '(1x,2I8, 2x, 3I4,2x, 2F15.2,F10.3,2x,I6)'
  fmt_2 = '(1x,2I8, 2x, 3I4,2x, 2F8.2, F10.3,2x,I6)'


  IF(.NOT. keyword_file) then
   call write_info('')
   call write_info('  !! FILE keyword mandatory for the search reflections procedure !!')
   call write_info('')
   return
  endif

  IF(search_equiv .AND. SPG%NumSpg == 0) then
   call write_info('')
   call write_info('  !! SPGR keyword mandatory for the search equivalent reflections procedure !!')
   call write_info('')
   return
  endif


  call write_info('')
  if(.not. search_H_string) then
   IF(hkl_absent(requested_H, SPG)) then
    WRITE(message_text, '(a,3I4,a)') '  ... Search for reflections: ', requested_H(1:3), ' (forbidden reflection)'
   else
    WRITE(message_text, '(a,3I4)')   '  ... Search for reflections: ', requested_H(1:3)
   endif
  else
   write(message_text, '(a,3(1x,a4))') ' ... Search for reflections: ', requested_H_string(1:3)
  endif
  call write_info(TRIM(message_text))
  call write_info('')
  call write_info('                      h   k   l        F2     sig    F2/sig    code')
  call write_info('')

  n_req = 0
  F2_req     = 0.
  wF2_req    = 0.
  sum_weight = 0.
  sig_F2_req = 0.
  ratio_req  = 0.
  HKL_search_ok = .FALSE.


  if(.not. search_H_string) then
   do i=1, n_ref
    ref_H(1) = INT(h(i))
    ref_H(2) = INT(k(i))
    ref_H(3) = INT(l(i))
    equiv_H = .false.
    IF(search_EQUIV) then
     equiv_H = HKL_equiv(INT(ref_H), INT(requested_H), SPG, search_friedel)
     IF((ref_H(1) == requested_H(1) .and. ref_H(2) == requested_H(2) .and. ref_H(3) == requested_H(3)) .or. equiv_H) then
      n_req = n_req + 1
      if(hkl_format_free) then
       WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
      else
       WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
      endif
     else
      cycle
     endif
    else
     IF(search_friedel) then
      IF((ref_H(1) ==  requested_H(1) .And. ref_H(2) ==  requested_H(2) .and. ref_H(3) ==  requested_H(3)) .OR.     &
         (ref_H(1) == -requested_H(1) .And. ref_H(2) == -requested_H(2) .and. ref_H(3) == -requested_H(3))) then
       n_req = n_req + 1
       if(hkl_format_free) then
        WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       else
        WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       endif
      else
       cycle
      endif
     else
      IF(ref_H(1) == requested_H(1) .And. ref_H(2) == requested_H(2) .and. ref_H(3) == requested_H(3)) then
       n_req = n_req + 1
       if(hkl_format_free) then
        WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       else
        WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       endif
      else
       cycle
      endif
     endif
    endif

    call write_info(TRIM(message_text))
    !n_req = n_req + 1
    F2_req = F2_req + F2(i)
    wF2_req = wF2_req + F2(i)*1/sig_F2(i)**2
    sig_F2_req = sig_F2_req + sig_F2(i)
    sum_weight = sum_weight + 1/sig_F2(i)**2
    ratio_req = ratio_req + F2(i)/sig_F2(i)
    HKL_search_ok(i) = .true.

   end do  ! fin de la boucle i=1, n_ref

  else
   ! recherche particuliere d'une serie de reflexion (ex: h 0 0)
   do i=1, n_ref
    ref_H(1) = INT(h(i))
    ref_H(2) = INT(k(i))
    ref_H(3) = INT(l(i))

    if(requested_H_string(1) == 'H') then
     if(requested_H_string(2) == 'K') then        ! recherche reflexion avec h et k quelconques
      read(requested_H_string(3), *) requested_H(3)
      IF(ref_H(3) == requested_H(3)) then
       n_req = n_req + 1
       if(hkl_format_free) then
        WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       else
        WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       endif
      else
       cycle
      endif

     elseif(requested_H_string(3) == 'L') then   ! recherche reflexion avec h et l quelconques
      read(requested_H_string(2), *) requested_H(2)
      IF(ref_H(2) == requested_H(2)) then
       n_req = n_req + 1
       if(hkl_format_free) then
        WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       else
        WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       endif
      else
       cycle
      endif

     else                                         ! recherche reflexion avec h quelconque
      read(requested_H_string(2), *) requested_H(2)
      read(requested_H_string(3), *) requested_H(3)
      IF(ref_H(2) == requested_H(2) .and. ref_H(3) == requested_H(3)) then
       n_req = n_req + 1
       if(hkl_format_free) then
        WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       else
        WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       endif
      else
       cycle
      endif
     endif

    elseif(requested_H_string(2) == 'K') then
     if(requested_H_string(3) == 'L') then     ! recherche reflexion avec k et l quelconques
      IF(ref_H(1) == requested_H(1)) then
       n_req = n_req + 1
       if(hkl_format_free) then
        WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       else
        WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       endif
      else
       cycle
      endif

     else                                          ! recherche reflexion avec k quelconque
      read(requested_H_string(1), *) requested_H(1)
      read(requested_H_string(3), *) requested_H(3)
      IF(ref_H(1) == requested_H(1) .and. ref_H(3) == requested_H(3)) then
       n_req = n_req + 1
       if(hkl_format_free) then
        WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       else
        WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
       endif
      else
       cycle
      endif
     endif

    elseif(requested_H_string(3) == 'L') then   ! recherche reflexion avec l quelconque
     read(requested_H_string(1), *) requested_H(1)
     read(requested_H_string(2), *) requested_H(2)
     IF(ref_H(1) == requested_H(1) .and. ref_H(2) == requested_H(2)) then
      n_req = n_req + 1
      if(hkl_format_free) then
       WRITE(message_text, fmt=trim(fmt_1)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
      else
       WRITE(message_text, fmt=trim(fmt_2)) n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i), F2(i)/sig_F2(i), code(i)
      endif
     else
      cycle
     endif
    endif

    call write_info(TRIM(message_text))
    F2_req = F2_req + F2(i)
    wF2_req = wF2_req + F2(i)*1/sig_F2(i)**2
    sig_F2_req = sig_F2_req + sig_F2(i)
    sum_weight = sum_weight + 1/sig_F2(i)**2
    ratio_req = ratio_req + F2(i)/sig_F2(i)
    HKL_search_ok(i) = .true.

   end do
  endif

  IF(n_req == 0) then
   call write_info('')
   call write_info('   >> No reflection founded !')
   call write_info('')
   return
  endif
  F2_mean  = F2_req/n_req
  !wF2_mean    = wF2_req/n_req
  !sig_F2_mean = sig_F2_req/n_req
  wF2_mean    = wF2_req/sum_weight
  sig_F2_mean = SQRT(1/sum_weight)
  ratio_mean  = ratio_req/n_req

  call write_info('')
  IF(n_req >1) then
   call calcul_Rint()
   WRITE(message_text,'(a,F8.2)') '    . Rint     = ', Rint
   call write_info(TRIM(message_text))
  endif
  WRITE(message_text,'(a,F8.2)') '    . <F2>     = ', wF2_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)') '    . <sig>    = ', sig_F2_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)') '    . <F2/sig> = ', ratio_mean
  call write_info(TRIM(message_text))



 RETURN
end subroutine search_hkl
!--------------------------------------------------------------------------------------------------------
subroutine search_hkl_EQUIV()
 USE hkl_module
 USE IO_module,                   ONLY : write_info
 USE cryscalc_module,             ONLY : message_text, keyword_FILE, SPG, debug_proc
 USE CFML_Reflections_Utilities,  ONLY : HKL_equiv, HKL_equiv_list, HKL_absent

 implicit none
  INTEGER                              :: i, ier, n_req
  REAL                                 :: F2_req, sig_F2_req, ratio_req
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: list_equiv_H
  INTEGER                              :: mult, n_equiv_H

  if(debug_proc%level_2)  call write_debug_proc_level(2, "search_HKL_equiv")

  IF(SPG%NumSpg == 0) then
   call write_info('')
   call write_info('  !! SPGR keyword mandatory for the search equivalent reflections procedure !!')
   call write_info('')
   return
  endif

  call allocate_HKL_arrays

  n_equiv_H = SPG%numops*2
  if (ALLOCATED(list_equiv_H)) DEALLOCATE(list_equiv_H)
  ALLOCATE(list_equiv_H(3,n_equiv_H), stat=ier)
   if(ier/=0) call write_alloc_error("list_equiv_H")

  ! ---------- feb. 2017 : check if systematic absence ----------------------
  if(hkl_absent(requested_H, SPG)) then
   call write_info("   Warning: Requested reflection IS a SYSTEMATIC absence in the current space group !")
   return   
  end if
  ! ------------------------------------------------------------------------- 
   
  call write_info('')
  WRITE(message_text, '(a,3I4)') '  ... Search for equivalent reflections of: ', requested_H(1), requested_H(2), requested_H(3)
  call write_info(TRIM(message_text))
  !if(search_Friedel) then
  ! call write_info(              '       Friedel reflections : yes')
  !else
  ! call write_info(              '       Friedel reflections : no')
  !end if
  call write_info('')
  n_req = 0
  F2_req     = 0.
  sig_F2_req = 0.
  ratio_req  = 0.
  HKL_search_ok = .FALSE.

  call HKL_Equiv_list(requested_H, SPG, search_Friedel, mult, list_equiv_H)
  do i=1, mult
    WRITE(message_text, '(1x,I4, 2x, 3I4)') i, list_equiv_H(1:3, i)
   call write_info(TRIM(message_text))
  END do

  IF(mult == 0) then
   call write_info('')
   call write_info('   >> No reflection founded !')
   call write_info('')
   if (ALLOCATED(list_equiv_H)) DEALLOCATE(list_equiv_H)
   return
  endif

  if (ALLOCATED(list_equiv_H)) DEALLOCATE(list_equiv_H)

 RETURN
end subroutine search_hkl_EQUIV

!--------------------------------------------------------------------------------------------------------
subroutine search_hkl_F2(input_string)
 USE hkl_module
 USE IO_module,                    ONLY : write_info
 USE cryscalc_module,              ONLY : message_text, keyword_FILE, keyword_find_HKL_ABSENT, SPG, debug_proc
 Use CFML_Reflections_Utilities,   ONLY : HKL_absent,  HKL_equiv

 implicit none
  CHARACTER(LEN=*), INTENT(IN) :: input_string
  INTEGER                      :: i, n_req, long_input_string
  REAL                         :: F2_req, sig_F2_req, ratio_req
  REAL                         :: wF2_req, sum_weight
  INTEGER, DIMENSION(3)        :: ref_H
  LOGICAL                      :: input_NEG, input_POS, input_ABSENT


  if(debug_proc%level_2)  call write_debug_proc_level(2, "search_hkl_F2 ("//trim(input_string)//")")

  long_input_string = len_trim(input_string)
  input_POS    = .false.
  input_NEG    = .false.
  input_ABSENT = .false.

  if(long_input_string == 3) then
   if(input_string(1:3) == 'POS') input_POS = .true.
   if(input_string(1:3) == 'NEG') input_NEG = .true.
  elseif(long_input_string == 6) then
   if(input_string(1:6) == 'ABSENT') input_ABSENT = .true.
  endif

  IF(.NOT. keyword_file) then
   call write_info('')
   call write_info('  !! FILE keyword mandatory for the search extinctions procedure !!')
   call write_info('')
   return
  endif

  IF(keyword_find_HKL_ABSENT .AND. SPG%NumSpg == 0) then
   call write_info('')
   call write_info('  !! SPGR keyword mandatory for the search absent reflections procedure !!')
   call write_info('')
   return
  endif

  call write_info('')
  IF(input_POS) then
   call write_info( '  ... Search for positive intensity reflections: ')
  ELSEIF(input_NEG) then
   call write_info( '  ... Search for negative intensity reflections: ')
  ELSEIF(input_ABSENT) then
   call write_info( '  ... Search for observed absent reflections: ')
  endif
  call write_info('')


  select case (input_string)
      case ('POS')
  n_req = 0
  F2_req     = 0.
  wF2_req    = 0.
  sig_F2_req = 0.
  sum_weight = 0.
  ratio_req  = 0.
  do i=1, n_ref
   IF(F2(i) < 0.01) cycle
   if(F2(i) < n_sig*sig_F2(i)) cycle
   n_req = n_req + 1

   IF(HKL_list_POS%OUT) then
    WRITE(message_text, '(1x,2I8,2x, 3I4,2F8.2)') n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i)
    call write_info(TRIM(message_text))
   endif
   IF(HKL_list_POS%WRITE) then
    if(cos_exist) then
     WRITE(message_text, '(3I4,2F8.2,I4, 6F8.5)') h(i), k(i), l(i), F2(i), sig_F2(i), code(i), cos_dir(i,1:6)
    else
     WRITE(message_text, '(3I4,2F8.2,I4)') h(i), k(i), l(i), F2(i), sig_F2(i), code(i)
    end if
    WRITE(31, '(a)') TRIM(message_text)
   endif
   F2_req     = F2_req     + F2(i)
   wF2_req    = wF2_req    + F2(i)*1/sig_F2(i)**2
   sum_weight = sum_weight + 1/sig_F2(i)**2
   sig_F2_req = sig_F2_req + sig_F2(i)
   ratio_req  = ratio_req  + F2(i)/sig_F2(i)
  end do

      case ('NEG')
  n_req = 0
  F2_req     = 0.
  wF2_req    = 0.
  sig_F2_req = 0.
  sum_weight = 0.
  ratio_req  = 0.
  do i=1, n_ref
   if (F2(i) > 0.) cycle
   n_req = n_req + 1

   IF(HKL_list_NEG%OUT) then
    WRITE(message_text, '(1x,2I8,2x, 3I4,2F8.2)') n_req, i, h(i), k(i), l(i), F2(i), sig_F2(i)
    call write_info(TRIM(message_text))
   endif
   IF(HKL_list_NEG%WRITE) then
    if(cos_exist) then
     WRITE(message_text, '(3I4,2F8.2,I4, 6F8.5)') h(i), k(i), l(i), F2(i), sig_F2(i), code(i), cos_dir(i,1:6)
    else
     WRITE(message_text, '(3I4,2F8.2,I4)') h(i), k(i), l(i), F2(i), sig_F2(i), code(i)
    end if
    WRITE(31, '(a)') TRIM(message_text)
   endif
   F2_req     = F2_req     + F2(i)
   wF2_req    = wF2_req    + F2(i)*1/sig_F2(i)**2
   sum_weight = sum_weight + 1/sig_F2(i)**2
   sig_F2_req = sig_F2_req + sig_F2(i)
   ratio_req  = ratio_req  + F2(i)/sig_F2(i)
  end do

      case ('ABSENT')
   n_req      = 0
   F2_req     = 0.
   wF2_req    = 0.
   sum_weight = 0.
   sig_F2_req = 0.
   ratio_req  = 0.
   do i=1, n_ref
    ref_H(1) = h(i)
    ref_H(2) = k(i)
    ref_H(3) = l(i)
    IF(.NOT. HKL_list_ABSENT%ALL) then
     if(F2(i) > 0. .and. F2(i)< n_sig*sig_F2(i)) cycle
     IF(F2(i) < 0.01) cycle
    endif

    IF(.not. HKL_absent(ref_H, SPG)) cycle
    n_req      = n_req + 1
    F2_req     = F2_req     + F2(i)
    wF2_req    = wF2_req    + F2(i)*1/sig_F2(i)**2
    sum_weight = sum_weight + 1/sig_F2(i)**2
    sig_F2_req = sig_F2_req + sig_F2(i)
    ratio_req = ratio_req   + F2(i)/sig_F2(i)
    IF(HKL_list_ABSENT%OUT) then
     WRITE(message_text, '(1x, 2I8,2x, 3I4, 2F8.2)') n_req, i, ref_H(1:3), F2(i), sig_F2(i)
     call write_info(TRIM(message_text))
    endif
    IF(HKL_list_ABSENT%WRITE) then
     if(cos_exist) then
      WRITE(message_text, '(3I4,2F8.2,I4,6F8.5)') ref_H(1:3), sig_F2(i), code(i), cos_dir(i,1:6)
     else
      WRITE(message_text, '(3I4,2F8.2,I4)') ref_H(1:3), sig_F2(i), code(i)
     end if
     call write_info(TRIM(message_text))
    endif
   end do

  end select


  IF(n_req == 0) then
   call write_info('')
   call write_info('  >> No reflection founded !')
   call write_info('')
   return
  endif

  F2_mean     = F2_req/n_req
  !wF2_mean    = wF2_req/n_req
  !sig_F2_mean = sig_F2_req/n_req

  wF2_mean    = wF2_req/sum_weight
  sig_F2_mean = SQRT(1/sum_weight)

  ratio_mean  = ratio_req/n_req

  call write_info('')
  IF(input_POS) then
   WRITE(message_text,'(a,I8)')    '   . Number of positive intensity reflections: ', n_req
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F4.0)') '     Criteria: F2/sig > ',  n_sig
   call write_info(TRIM(message_text))

  ELSEIF(input_NEG) then
   WRITE(message_text,'(a,I8)')    '   . Number of negative intensity reflections: ', n_req
   call write_info(TRIM(message_text))

  ELSEIF(input_ABSENT) then
   WRITE(message_text, '(a,I8)')   '   . Total number of reflections:             ', n_ref
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,I8)')   '   . Number of authorized reflections:        ', n_ref - n_req
   call write_info(TRIM(message_text))

   WRITE(message_text,'(a,I8)')    '   . Number of systematic absence violations: ', n_req
   call write_info(TRIM(message_text))
   IF(.NOT. HKL_list_ABSENT%ALL) then
    call write_info(           '     Criterias: . F2 > 0.')
    WRITE(message_text, '(a,F4.0)')'                . F2/sig > ',n_sig
    call write_info(TRIM(message_text))
   else
    call write_info(           '     Criteria: no')
   endif
  endif

  !WRITE(message_text,'(a,F8.2)')   '    . <F2>     = ', F2_mean
  !call write_info(TRIM(message_text))
  !WRITE(message_text,'(a,F8.2)')   '    . <wF2>    = ', wF2_mean
  !call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')   '    . <F2>     = ', wF2_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')   '    . <sig>    = ', sig_F2_mean
  call write_info(TRIM(message_text))
  WRITE(message_text,'(a,F8.2)')   '    . <F2/sig> = ', ratio_mean
  call write_info(TRIM(message_text))

  IF(HKL_list_POS%WRITE .or. HKL_list_NEG%write .or. HKl_list_ABSENT%write) then
   CLOSE(UNIT=31)
   call write_info('')
   WRITE(message_text, '(3a)') '  > The ', TRIM(HKL_file%OUTPUT), ' file has been created.'
   call write_info(trim(message_text))
   call write_info('')
  endif


 RETURN
end subroutine search_hkl_F2

!------------------------------------------------------------
subroutine search_HKL_list()
 USE hkl_module
 USE cryscalc_module, ONLY : message_text, keyword_file, HKL_list_out1_unit, HKL_list_out2_unit, file_out_n, debug_proc
 USE IO_module,       ONLY : write_info
 USE macros_module,   ONLY : multiple
  implicit none
  INTEGER           :: i
  REAL              :: sig_impair, ratio_impair, sum_weight

  if(debug_proc%level_2)  call write_debug_proc_level(2, "search_HKL_list")


  IF(.NOT. keyword_file) then
   call write_info('')
   call write_info('  !! FILE keyword mandatory for the search extinctions procedure !!')
   call write_info('')
   return
  endif


 call write_info('')
 if(HKL_list%EXTI_number > 0) then
  WRITE(message_text, '(2a)') '  ... Search for reflections: ', TRIM(HKL_rule(HKL_list%EXTI_number))
 else
  WRITE(message_text, '(2a)') '  ... Search for reflections: ', TRIM(HKL_rule(abs(HKL_list%EXTI_number)+HKL_rule_nb))
 endif
 call write_info(TRIM(message_text))
 call write_info('')

 call search_HKL_exti(HKL_list%EXTI_number)

 n_impair     = 0
 F2_impair    = 0.
 sig_impair   = 0.
 ratio_impair = 0.
 sum_weight   = 0.

 do i=1, n_ref
  if (.not. HKL_search_ok(i)) then
   IF(HKL_list%SUPPRESS) then
    if(cos_exist) then
     WRITE(message_text, '(3I4, 2F8.2,I4,6F8.5)') h(i), k(i), l(i), F2(i), sig_F2(i), code(i), cos_dir(i,1:6)
    else
     WRITE(message_text, '(3I4, 2F8.2,I4)')       h(i), k(i), l(i), F2(i), sig_F2(i), code(i)
    endif
    write(HKL_list_out2_unit,'(a)') TRIM(message_text)
   endif
   cycle
  endif

  n_impair     = n_impair + 1

  IF(HKL_list%OUT .and. multiple(n_impair,file_out_n)) then
   !WRITE(message_text, '(I6,x,I6,2x, 3I4, 2F8.2)') n_impair, i, h(i), k(i), l(i), F2(i), sig_F2(i)
   WRITE(message_text, '(I8, 1x, I8,2x, 3I4, 2F8.2)') n_impair, i, h(i), k(i), l(i), F2(i), sig_F2(i)
   call write_info(TRIM(message_text))
  endif

  IF(HKL_list%WRITE) then
   if(cos_exist) then
    WRITE(message_text, '(3I4, 2F8.2,I4,6F8.5)') h(i), k(i), l(i), F2(i), sig_F2(i), code(i), cos_dir(i,1:6)
   else
    WRITE(message_text, '(3I4, 2F8.2,I4)')       h(i), k(i), l(i), F2(i), sig_F2(i), code(i)
   endif
   write(HKL_list_out1_unit,'(a)') TRIM(message_text)
  endif


  F2_impair    = F2_impair  + F2(i)
  sum_weight   = sum_weight + 1/sig_F2(i)**2

  sig_impair   = sig_impair + sig_F2(i)
  ratio_impair = ratio_impair + F2(i)/sig_F2(i)
 end do

 F2_mean       = F2_impair/n_impair
 !sig_F2_mean   = sig_impair/n_impair
 sig_F2_mean   = SQRT(1/sum_weight)
 ratio_mean    = ratio_impair/n_impair

 IF(n_impair > 1) call calcul_Rint()

 call write_info('')
 if(HKL_list%EXTI_number > 0) then
  WRITE(message_text, '(2a)')  '  >>>> reflections : ', TRIM(HKL_rule(HKL_list%EXTI_number))
 else
  WRITE(message_text, '(2a)')  '  >>>> reflections : ', TRIM(HKL_rule(abs(HKL_list%EXTI_number) + HKL_rule_nb) )
 endif
 call write_info(TRIM(message_text))
 call write_info('')
 IF(.NOT. HKL_list%ALL) then
  call write_info(            '    Criterias: . F2 > 0.')
  WRITE(message_text, '(a,F4.0)')'               . F2/sig > ',n_sig
  call write_info(TRIM(message_text))
 else
  call write_info(           '    Criteria: no')
 endif
 WRITE(message_text, '(a,I8)')   '         . number of reflections: ', n_impair
 call write_info(TRIM(message_text))
 IF(n_impair ==0) return

 IF(n_impair > 1) then
  call write_info('')
  WRITE(message_text, '(a,F8.2)') '         . Rint     = ', Rint
 endif
 !call write_info(TRIM(message_text))
 WRITE(message_text, '(a,F8.2)') '         . <F2>     = ', F2_mean
 call write_info(TRIM(message_text))
 WRITE(message_text, '(a,F8.2)') '         . <sig>    = ', sig_F2_mean
 call write_info(TRIM(message_text))
 WRITE(message_text, '(a,F8.2)') '         . <F2/sig> = ', ratio_mean
 call write_info(TRIM(message_text))
 call write_info('')

 IF(HKL_list%WRITE) then
  CLOSE(UNIT=HKL_list_out1_unit)
  call write_info('')
  WRITE(message_text, '(3a)') '  > The ', TRIM(HKL_file%OUTPUT), ' file has been created.'
  call write_info(trim(message_text))
  call write_info('')
 endif

 IF(HKL_list%SUPPRESS) then
  CLOSE(UNIT=HKL_list_out2_unit)
  call write_info('')
  WRITE(message_text, '(3a)') '  > The ', TRIM(HKL_file%OUTPUT2), ' file has been created.'
  call write_info(trim(message_text))
  call write_info('')
 endif

 RETURN
end subroutine search_HKL_list
!-------------------------------------------------------
subroutine list_EXTI_RULE
USE HKL_module
USE IO_module,       ONLY : write_info
USE cryscalc_module, ONLY : message_text, debug_proc
 implicit none
  INTEGER             :: i

  if(debug_proc%level_2)  call write_debug_proc_level(2, "list_exti_rule")

  call def_HKL_rule

  do i=1, HKL_rule_nb
   IF(i<10) then
    WRITE(message_text, '(a,i1,a)') '  #',i, TRIM(HKL_rule(i))
   else
    WRITE(message_text, '(a,i2,a)')  ' #',i, TRIM(HKL_rule(i))
   endif
   call write_info(TRIM(message_text))
  end do

 return
END subroutine list_EXTI_RULE
!-------------------------------------------------------
subroutine calcul_Rint()
 USE cryscalc_module, only : debug_proc
 USE HKL_module
 implicit none
  INTEGER             :: i
  REAL                :: F2_sum

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_Rint")

  Rint = 0.
  n_ref_eff = 0
  F2_sum = 0
  do i=1, n_ref
   if (.not. HKL_search_ok(i)) cycle
   Rint =  Rint + ABS(F2(i) - F2_mean)
   F2_sum = F2_sum + F2(i)
   n_ref_eff = n_ref_eff + 1
  end do
 ! Rint = 100*(Rint / MAX(1.,F2_sum))
  Rint = 100*Rint / F2_sum

 RETURN
end subroutine calcul_Rint
!-------------------------------------------------------

 subroutine search_monoclinic_angle
  use CRYSCALC_module,                only : keyword_FILE, space_group_symbol, SPG, message_text, unit_cell, &
                                             search_mono_criteria, write_details
  use CFML_crystallographic_symmetry, ONLY : set_spacegroup
  USE HKL_module,                     ONLY : Rint
  USE IO_module,                      ONLY : write_info
   implicit none
    real, dimension(3)               :: R_int
    integer, dimension(1)            :: i_mono
    integer                          :: i, j, n_too_far, ok
    character (len=16), dimension(3) :: tested_SG
    real, dimension(3)               :: diff, diff_ratio
    character (len=8), dimension(3)  :: angle_string, cell_string
    logical, dimension(3)            :: too_far
    character (len=8), dimension(3)  :: mat_string


 ! recherche de l'angle monoclinic dans une maille pseudo-orthorhombique (3 angles à 90 deg.)
   cell_string(1) = 'cell_a'
   cell_string(2) = 'cell_b'
   cell_string(3) = 'cell_c'
   angle_string(1) = 'alfa'
   angle_string(2) = 'beta'
   angle_string(3) = 'gamma'
   mat_string(1)   = 'c a b'
   mat_string(2)   = 'a b c'
   mat_string(3)   = 'b c a'
   tested_SG(1)    = 'P 2 1 1'
   tested_SG(2)    = 'P 1 2 1'
   tested_SG(3)    = 'P 1 1 2'

   if(write_details) then
   do i=1,3
    write(message_text, '(5x,3a,F5.2,a)')  '. ', cell_string(i)(1:6), ' = ', unit_cell%param(i), ' A'
    call write_info(trim(message_text))
   end do
   do i=1,3
    write(message_text, '(5x,3a,F5.2,a)')  '. ', angle_string(i)(1:6), ' = ', unit_cell%param(i+3), ' deg.'
    call write_info(trim(message_text))
   end do
   write(message_text, '(5x,a,F5.2,a)')  '. delta_max  = ', search_mono_criteria(1), ' deg.'
   call write_info(trim(message_text))
   write(message_text, '(5x,a,F6.2,a)')  '. ratio_Rint = ', 100*search_mono_criteria(2), ' %'
   call write_info(trim(message_text))
   call write_info('')
   end if

   too_far = .false.
   n_too_far = 0
   do i=1, 3
    diff(i) = abs(90. - unit_cell%param(i+3))
    if(diff(i) > search_mono_criteria(1)) then
     too_far(i) = .true.
     n_too_far = n_too_far + 1
    end if
   end do
   if(write_details) then
    if(n_too_far > 1 ) then
     do i=1, 3
      if(too_far(i)) then
       write(message_text, '(3a,F10.5,a,F5.2,a)') '    => Angle ', angle_string(i)(1:5), ' (',  unit_cell%param(i+3), &
                                                  ') too far from 90.0 deg. !! [delta_max = ',search_mono_criteria(1), ']'
       call write_info(trim(message_text))
      end if
     end do
     return
    end if
   end if


   IF(.NOT. keyword_file) then
   call write_info('')
   call write_info('  !! FILE keyword mandatory for the Rint calculation procedure !!')
   call write_info('')
   return
  endif

 ! calcul du Rint pour les 3 groupes "P211", "P121" et "P112"
  do i=1, 3
   R_int(i) = 100.
   ! test si les 2 angles non monocliniques sont proches de 90.  --------
   if(write_details) then
    write(message_text, '(2a)') '    ..... Check if monoclinic angle is ', angle_string(i)(1:5)
    call write_info(trim(message_text))
   end if

   too_far = .false.
   n_too_far = 0
   do j=1, 3
    if (j/=i) then
     diff(j) = abs(90. - unit_cell%param(j+3))
     if(diff(j) > search_mono_criteria(1)) then
      too_far(j) = .true.
      n_too_far = n_too_far + 1
      if(write_details) then
       write(message_text, '(5x,3a,F10.5,a,F5.2,a)') '=> Angle ', angle_string(j)(1:6), ' (',  unit_cell%param(j+3), &
                                                     ') too far from 90.0 deg. !! [delta_max = ',search_mono_criteria(1), ']'
       call write_info(trim(message_text))
       call write_info('')
      end if
     !cycle
     end if
    end if
   end do
   if(n_too_far /=0) cycle
   !----------------------------------------------------------------------

   call set_spacegroup(trim(tested_SG(i)), SPG)
   if(write_details) then
    call calcul_global_Rint("mono", "out")
   else
    call calcul_global_Rint("mono", "no_out")
   endif
   R_int(i) = Rint
   if(write_details) then
    write(message_text, '(8x,a,F6.2)') ' Rint : ', Rint
    call write_info(trim(message_text))
    call write_info('')
   end if
  end do

  ! determination du Rint le plus faible
  i_mono = MINLOC(R_int)
  ok = 0
  do i=1, 3
   if(i == i_mono(1)) cycle
   diff(i) = abs(R_int(i) - MINVAL(R_int))
   diff_ratio(i) = diff(i)/MINVAL(R_int)
   !write (*,*) ' ABS(Rint - Rint_min) : ', diff(i) , diff_perc(i)
   if(diff_ratio(i) > search_mono_criteria(2)) then
    ok = ok + 1
   end if
  end do

  if(ok == 2 .or. n_too_far /=0) then
   if(i_mono(1) /=2) then
    write(message_text, '(3a,F10.5,a)')      '   >> Probable monoclinic angle : ', angle_string(i_mono(1)), &
                                             '=',unit_cell%param(3+i_mono(1)),' deg.'
    call write_info(trim(message_text))
    write(message_text, '(a,F6.2)')  '      Rint : ', MINVAL(R_int)
    call write_info(trim(message_text))
    write(message_text, '(2a)')     '      > Transformation matrix to get beta monoclinic angle : ', trim(mat_string(i_mono(1)))
    call write_info(trim(message_text))
    if(i_mono(1) == 1) then
     write(message_text, '(8x,3a,F5.2,a)')  '. ', cell_string(1)(1:6), ' = ', unit_cell%param(3), ' A'
     call write_info(trim(message_text))
     write(message_text, '(8x,3a,F5.2,a)')  '. ', cell_string(2)(1:6), ' = ', unit_cell%param(1), ' A'
     call write_info(trim(message_text))
     write(message_text, '(8x,3a,F5.2,a)')  '. ', cell_string(3)(1:6), ' = ', unit_cell%param(2), ' A'
     call write_info(trim(message_text))

     write(message_text, '(8x,3a,F5.2,a)')  '. ', angle_string(1)(1:6), ' = ', unit_cell%param(6), ' deg.'
     call write_info(trim(message_text))
     write(message_text, '(8x,3a,F5.2,a)')  '. ', angle_string(2)(1:6), ' = ', unit_cell%param(4), ' deg.'
     call write_info(trim(message_text))
     write(message_text, '(8x,3a,F5.2,a)')  '. ', angle_string(3)(1:6), ' = ', unit_cell%param(5), ' deg.'
     call write_info(trim(message_text))
    elseif(i_mono(1) == 3) then
     write(message_text, '(8x,3a,F5.2,a)')  '. ', cell_string(1)(1:6), ' = ', unit_cell%param(2), ' A'
     call write_info(trim(message_text))
     write(message_text, '(8x,3a,F5.2,a)')  '. ', cell_string(2)(1:6), ' = ', unit_cell%param(3), ' A'
     call write_info(trim(message_text))
     write(message_text, '(8x,3a,F5.2,a)')  '. ', cell_string(3)(1:6), ' = ', unit_cell%param(1), ' A'
     call write_info(trim(message_text))

     write(message_text, '(8x,3a,F5.2,a)')  '. ', angle_string(1)(1:6), ' = ', unit_cell%param(5), ' deg.'
     call write_info(trim(message_text))
     write(message_text, '(8x,3a,F5.2,a)')  '. ', angle_string(2)(1:6), ' = ', unit_cell%param(6), ' deg.'
     call write_info(trim(message_text))
     write(message_text, '(8x,3a,F5.2,a)')  '. ', angle_string(3)(1:6), ' = ', unit_cell%param(4), ' deg.'
     call write_info(trim(message_text))

    end if

   else
    write(message_text, '(a)')      '   >> Probable monoclinic angle : beta (no transformation)'
    call write_info(trim(message_text))
    write(message_text, '(a,F6.2)')  '      Rint : ', MINVAL(R_int)
    call write_info(trim(message_text))
   end if
   !call write_info(trim(message_text))
   call write_info('')
  else
   write(message_text, '(a,F5.2,a)')      '   >> No obvious monoclinic angle !! [criteria : ', search_mono_criteria(2), ']'
   call write_info(trim(message_text))

  end if


  return
 end subroutine search_monoclinic_angle

 !-------------------------------------------------------

 subroutine search_tetragonal_axis
  use CRYSCALC_module,                only : keyword_FILE, space_group_symbol, SPG, message_text, unit_cell, &
                                             search_tetra_criteria, Mat, mat_det, Mat_integer, get_SPGR, crystal_system, &
                                             write_details
  use CFML_crystallographic_symmetry, ONLY : set_spacegroup
  USE HKL_module,                     ONLY : HKL_file, Rint, ordered_hkl, HKL_data_known
  USE Matrix_list_module

  USE IO_module,                      ONLY : write_info
   implicit none
    real, dimension(3)                :: R_int
    integer, dimension(1)             :: i_tetra
    integer                           :: i, i1, n_too_far, ok
    character (len=16), dimension(3)  :: tetra_SPG
    character (len=8),  dimension(3)  :: tetra_setting
    real                              :: mean_cell_param
    real, dimension(3)                :: diff, diff_ratio
    character (len=1), dimension(3)   :: param_string
    character (len=8), dimension(3)   :: angle_string
    logical, dimension(3)             :: too_far
    character (len=8), dimension(3)   :: mat_string
    character (len=256), dimension(3) :: ini_file, transf_file
    character (len=256)               :: command_DOS, file_prefixe


 ! recherche de l'axe tetrafgonal dans une maille pseudo-cubic
   param_string(1) = 'a'
   param_string(2) = 'b'
   param_string(3) = 'c'
   angle_string(1) = 'alfa'
   angle_string(2) = 'beta'
   angle_string(3) = 'gamma'
   mat_string(1)   = 'b c a'
   mat_string(2)   = 'a b c'
   mat_string(3)   = 'b c a'
   tetra_setting(1) = 'abc'
   tetra_setting(2) = 'cab'
   tetra_setting(3) = 'bca'

   i1 = index(HKL_file%hkl, '.hkl')
   write(file_prefixe,   '(a)')  HKL_file%hkl(1:i1-1)
   write(ini_file(1),    '(2a)') trim(file_prefixe), '.hkl'
   write(ini_file(2),    '(2a)') trim(file_prefixe), '_trans.hkl'
   write(ini_file(3),    '(2a)') trim(file_prefixe), '_trans_trans.hkl'
   write(transf_file(1), '(2a)') trim(file_prefixe), '_abc.hkl'
   write(transf_file(2), '(2a)') trim(file_prefixe), '_cab.hkl'
   write(transf_file(3), '(2a)') trim(file_prefixe), '_bca.hkl'



   do i=1,3
    write(message_text, '(5x,3a,F10.5,a)')  '. ', param_string(i), ' = ', unit_cell%param(i), ' A'
    call write_info(trim(message_text))
   end do
   do i=1,3
    write(message_text, '(5x,3a,F6.2,a)')  '. ', angle_string(i)(1:5), ' = ', unit_cell%param(i+3), ' deg.'
    call write_info(trim(message_text))
   end do
   write(message_text, '(5x,a,F5.2,a)')  '. delta_max_parameters  = ', search_tetra_criteria(1), ' A'
   call write_info(trim(message_text))
   write(message_text, '(5x,a,F5.2,a)')  '. delta_max_angles      = ', search_tetra_criteria(2), ' deg.'
   call write_info(trim(message_text))
   write(message_text, '(5x,a,F6.2,a)')  '. ratio_Rint            = ', 100*search_tetra_criteria(3), ' %'
   call write_info(trim(message_text))
   call write_info('')


   too_far = .false.
   n_too_far = 0
   mean_cell_param = SUM(unit_cell%param(1:3)) / 3
   do i=1, 3
    diff(i) = abs(mean_cell_param - unit_cell%param(i))
    if(diff(i) > search_tetra_criteria(1)) then
     too_far(i) = .true.
     n_too_far = n_too_far + 1
    end if
   end do
   if(n_too_far > 0 ) then
    do i=1, 3
     if(too_far(i)) then
      write(message_text, '(3a,F10.5,a,F5.2,a)') '    => Cell parameter ', param_string(i)(1:1), ' (',  unit_cell%param(i), &
                                                 ') too far from the others !! [delta_max = ',search_tetra_criteria(1), ' A]'
      call write_info(trim(message_text))
     end if
    end do
    return
   end if

   too_far = .false.
   n_too_far = 0
   do i=1, 3
    diff(i) = abs(90. - unit_cell%param(i+3))
    if(diff(i) > search_tetra_criteria(1)) then
     too_far(i) = .true.
     n_too_far = n_too_far + 1
    end if
   end do
   if(n_too_far > 0 ) then
    do i=1, 3
     if(too_far(i)) then
      write(message_text, '(3a,F10.5,a,F5.2,a)') '    => Angle ', angle_string(i)(1:5), ' (',  unit_cell%param(i+3), &
                                                 ') too far from 90.0 deg. !! [delta_max = ',search_tetra_criteria(1), ']'
      call write_info(trim(message_text))
     end if
    end do
    return
   end if


   IF(.NOT. keyword_file) then
   call write_info('')
   call write_info('  !! FILE keyword mandatory for the Rint calculation procedure !!')
   call write_info('')
   return
  endif

  ! setting a b c
  get_SPGR     = .true.
  crystal_system = "TETRA"
  call search_SPGR('no_out')
  tetra_SPG(1) = trim(SPG%SPG_Symb)
  if(write_details) then
   call calcul_global_Rint("", "out")
  else
   call calcul_global_Rint("", "no_out")
  end if
  R_int(1) = Rint

  ! setting c a b
  Mat(:,:) = transf_mat(:,:, 7)    ! c a b
  Mat_det = 1.
  Mat_integer = .true.
  call write_info('    >>> Matrix transformation : c a b')
  call transf_cell_parameters
  call transf_HKL_file
  HKL_file%NAME = trim(ini_file(2))
  HKL_data_known = .false.
  HKL_file%SHELX = .true.
  HKL_file%CIF   = .false.
  ordered_hkl    = .false.
  if(write_details) then
   call read_and_sort_hkl('sort', 'out')
  else
   call read_and_sort_hkl('sort', 'no_out')
  endif
  get_SPGR     = .true.
  if(write_details) then
   call search_SPGR('out')
  else
   call search_SPGR('no_out')
  endif
  tetra_SPG(2) = trim(SPG%SPG_Symb)
  if(write_details) then
  call calcul_global_Rint("", "out")
  else
  call calcul_global_Rint("", "no_out")
  endif
  R_int(2) = Rint
  !write(*,*) ' Setting cab : Rint = ', Rint, '  ', trim(SPG%SPG_Symb)

  !pause
  ! setting b c a
  Mat(:,:) = transf_mat(:,:, 7)    ! c a b
  Mat_det = 1.
  Mat_integer = .true.
  call write_info('    >>> Matrix transformation : c a b')
  call transf_cell_parameters
  call transf_HKL_file
  HKL_file%NAME = trim(ini_file(3))
  HKL_data_known = .false.
  HKL_file%SHELX = .true.
  HKL_file%CIF   = .false.
  ordered_hkl    = .false.
  if(write_details) then
  call read_and_sort_hkl('sort', 'out')
  else
  call read_and_sort_hkl('sort', 'no_out')
  end if
  get_SPGR     = .true.
  if(write_details) then
  call search_SPGR('out')
  else
  call search_SPGR('no_out')
  end if
  tetra_SPG(3) = trim(SPG%SPG_Symb)
  if(write_details) then
  call calcul_global_Rint("", "out")
  else
  call calcul_global_Rint("", "no_out")
  end if
  R_int(3) = Rint
   !write(*,*) ' Setting bca : Rint = ', Rint, '  ', trim(SPG%SPG_Symb)


 ! calcul du Rint pour les 3 settings abc, cab et bca
  do i=1, 3
   write(message_text, '(8x,2a,F6.2,3a)') trim(tetra_setting(i)),' setting : Rint = ', R_int(i) , '   [', trim(tetra_SPG(i)),']'
   call write_info(trim(message_text))
  end do
  call write_info('')

  ! determination du Rint le plus faible
  i_tetra = MINLOC(R_int)
  ok = 0
  do i=1, 3
   if(i == i_tetra(1)) cycle
   diff(i) = abs(R_int(i) - MINVAL(R_int))
   diff_ratio(i) = diff(i)/MINVAL(R_int)
   if(diff_ratio(i) > search_tetra_criteria(3)) then
    ok = ok + 1
   end if
  end do

  if(ok /=0) then
   write(message_text, '(5a)')  '     >> Probable tetragonal axis : ', param_string(i_tetra(1)), ' (setting = ', &
                                trim(tetra_setting(i_tetra(1)))   ,')'
   call write_info(trim(message_text))
   write(message_text, '(2a)')  '     >> HKL file : ', trim(transf_file(i_tetra(1)))
   call write_info(trim(message_text))
   call write_info("")
   if(i_tetra(1) == 1) then
    write(command_dos, '(a)') 'copy import.hkl ', trim(transf_file(1))
   elseif(i_tetra(1) == 2) then
    write(command_dos, '(2a)') 'copy import_trans.hkl ', trim(transf_file(2))
   elseif(i_tetra(1) == 3) then
    write(command_dos, '(2a)') 'copy import_trans_trans.hkl ', trim(transf_file(3))
   end if
   call system(trim(command_dos))
   do i=1, 3
    write(command_dos, '(2a)') 'del ', trim(ini_file(i))
   end do
  else
   write(message_text, '(a,F5.2,a)') '   >> No obvious preferential tetragonal axis !! [criteria : ', search_tetra_criteria(3), ']'
   call write_info(trim(message_text))
  end if


  !if(ok == 2) then
  ! if(i_mono(1) /=2) then
  !  write(message_text, '(3a,F10.5,a)')      '   >> Probable monoclinic angle : ', angle_string(i_mono(1)), &
  !                                         '=',unit_cell%param(3+i_mono(1)),' deg.'
    !call write_info(trim(message_text))
    !write(message_text, '(a,F6.2)')  '      Rint : ', MINVAL(R_int)
    !call write_info(trim(message_text))
    !write(message_text, '(2a)')     '      > Transformation matrix to get beta monoclinic angle : ', trim(mat_string(i_mono(1)))
   !else
    !write(message_text, '(a)')      '   >> Probable monoclinic angle : beta (no transformation)'
    !call write_info(trim(message_text))
    !write(message_text, '(a,F6.2)')  '      Rint : ', MINVAL(R_int)
   !end if
   !call write_info(trim(message_text))
   !call write_info('')
  !else
   !write(message_text, '(a,F5.2,a)')      '   >> No obvious monoclinic angle !! [criteria : ', search_mono_criteria(2), ']'
   !call write_info(trim(message_text))
  !end if


  return
 end subroutine search_tetragonal_axis


!-------------------------------------------------------
subroutine calcul_global_Rint(input_string, input_string_2)
! extrait de 'datared.F90' JRC
 USE cryscalc_module,             ONLY : SPG, keyword_FILE, keyword_CELL, keyword_WAVE, message_text, known_theta, &
                                         WRITE_details, debug_proc

 USE HKL_module
 Use CFML_Reflections_Utilities,  ONLY : HKL_absent,  HKL_equiv
 USE CFML_Math_General,           ONLY : sort
 USE IO_module,                   ONLY : write_info
 implicit none
  CHARACTER (LEN=*), INTENT(IN)    :: input_string
  CHARACTER (LEN=*), INTENT(IN)    :: input_string_2
  CHARACTER  (LEN=8)               :: out_string
  INTEGER                          :: i,i1, j,j2
  REAL                             :: F2_sum

  INTEGER, DIMENSION(n_ref) :: itreat !
  INTEGER, DIMENSION(n_ref) :: n_equiv, ini, fin
  REAL,    DIMENSION(n_ref) :: weight, F2_av, sig_av
  INTEGER, DIMENSION(n_ref) :: redundancy
  INTEGER, DIMENSION(3)     :: H1, H2
  LOGICAL                   :: absent
  INTEGER                   :: n_rejected ! nombre de reflections rejetees: eteintes par le groupe d'esapce
  INTEGER                   :: n_ok       ! nombre de reflections conservees
  INTEGER                   :: n_ok_2     ! nombre de reflections conservees avec I>2sig
  INTEGER                   :: n_ok_3     ! nombre de reflections conservees avec I>3sig
  INTEGER                   :: ns
  LOGICAL                   :: merge, friedel
  REAL                      :: total, sig, sigg
  REAL                      :: sum_an, sum_aw, sum_anw
  REAL                      :: aver_sig, Rwint
  REAL                      :: sum_weight
  INTEGER                   :: long_input_string
  INTEGER                   :: ier

  INTEGER,  ALLOCATABLE, DIMENSION(:) :: ordered_array

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_global_Rint ("//trim(input_string)//")")

  long_input_string = len_trim(input_string)
  merge = .false.
  if(input_string == "merge") merge = .true.


  IF(SPG%NumSpg == 0) then
   call write_info('')
   call write_info('  !! SPGR keyword mandatory for the Rint calculation procedure !!')
   call write_info('')
   return
  endif

  IF(.NOT. keyword_file) then
   call write_info('')
   call write_info('  !! FILE keyword mandatory for the Rint calculation procedure !!')
   call write_info('')
   return
  endif

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

  IF(merge) then
   i1 = INDEX(HKL_file%NAME, '.')
   IF(i1 /=0) then
    WRITE(HKL_file%MERGE , '(2a)') HKL_file%NAME(1:i1-1), '_merge.hkl'
   else
    WRITE(HKL_file%MERGE, '(a)') 'merge.hkl'
   endif
   OPEN(UNIT = 31, FILE =TRIM(HKL_file%MERGE))
  endif



  !
  IF(.NOT. known_theta) call calcul_theta()

  if(.not. WRITE_details) then
   out_string = 'no_out'
  else
   out_string = input_string_2
  end if

  ! Order the reflections by ascending Theta
  if (.not. ordered_hkl) then
   IF(ALLOCATED(ordered_array)) DEALLOCATE(ordered_array)
   ALLOCATE(ordered_array(n_ref), stat=ier)
    if(ier/=0) call write_alloc_error("ordered_array")
   if(WRITE_details) then
   call write_info('')
   call write_info('  ... Sort reflections in ascending theta values ...')
   end if
   call sort(theta_hkl, n_ref, ordered_array)
   if(len_trim(out_string) == 6) then
    if(out_string(1:6) == 'no_out') then
     call write_info('')
     call sort_arrays('theta','+', ordered_array, theta_hkl,'no_out')
    end if
   else
    call sort_arrays('theta','+', ordered_array, theta_hkl,'out')
   end if
   call write_info('')
   ordered_hkl =.true.
  END if
  !


  write(message_text, '(3a)') '   ... R_internal calculation in ', trim(SPG%SPG_Symb), ' ...'
  call write_info(trim(message_text))
  if(.not. input_string == 'mono') call write_info('')

  friedel = .true.
  IF(SPG%centred /= 2) friedel=.false.

  n_ok       = 0
  n_ok_2     = 0
  n_ok_3     = 0
  itreat(:)  = 0
  n_equiv(:) = 0
  F2_av(:)   = 0.
  ini(:)     = 0
  fin(:)     = 0
  total      = 0
  n_rejected = 0
  F2_sum     = 0.
  absent     = .true.
  redundancy = 0.


  do i=1, n_ref    ! boucle sur les n_ref reflections
   IF(itreat(i) /=0) cycle
   H1(1) = h(i)
   H1(2) = k(i)
   H1(3) = l(i)
   absent = hkl_absent(H1, SPG)
   IF(absent) then
    n_rejected = n_rejected + 1
    cycle
   endif

   n_ok = n_ok + 1          ! update the number of independent reflections
   !IF(MOD(i,1000) == 0) then
   ! WRITE(message_text, '(a,3i10)') '  ==> reflection: ', i, n_ok, n_rejected
   ! call write_info(trim(message_text))
   !endif

   IF(F2(i) > 2*sig_F2(i)) n_ok_2 = n_ok_2 + 1
   IF(F2(i) > 3*sig_F2(i)) n_ok_3 = n_ok_3 + 1
   itreat(i) = i            ! Make this reflection treated
   sig = 1.0/sig_F2(i)**2
   ini(n_ok) = i
   fin(n_ok) = i
   n_equiv(n_ok) = 1

   do j = i+1, n_ref  ! boucles sur les reflections suivantes: recherche des reflections  equivalentes
    IF(known_theta .and. (ABS(theta_hkl(i) - theta_HKL(j)) > 0.2)) exit
    H2(1) = h(j)
    H2(2) = k(j)
    H2(3) = l(j)
    if (HKL_equiv(H1, H2, SPG, friedel)) then

     itreat(j)     = i
     n_equiv(n_ok) = n_equiv(n_ok) + 1
     fin(n_ok)     = j
     sig           = sig + 1.0/sig_F2(j)**2  ! Sw: Somme des 1/sigmas**
    end if
   end do

   ns = 0
   do j = ini(n_ok), fin(n_ok)
    if (itreat(j) == i) then
     ns = ns + 1
     weight(ns) = (1.0/sig_F2(j)**2)/sig     ! w = wj / Sw
    end if
   end do

   F2_sum = 0
   ns     = 0
   do j=ini(n_ok), fin(n_ok)
    IF(itreat(j) == i) then
     ns = ns + 1
     F2_sum = F2_sum + F2(j)*weight(ns)      ! Somme des F2
    endif
   end do

   F2_av(n_ok) = F2_sum
   F2_sum = 0.
   ns     = 0
   do j=ini(n_ok), fin(n_ok)
    IF(itreat(j) == i) then
     ns=ns+1
     F2_sum = F2_sum + weight(ns)*(F2_av(n_ok) - F2(j))**2    !
    endif
   end do
   sig_av(n_ok) = SQRT(F2_sum)
   sigg = SUM(sig_F2(ini(n_ok):fin(n_ok))) / MAX(1.0, REAL(fin(n_ok) - ini(n_ok)))
   IF(sig_av(n_ok) < sigg) sig_av(n_ok) = sigg
   total = total + sig_av(n_ok)

   redundancy(n_ok) = real(ns)


   !write(message_text,fmt="(a,3i4,3i6,2f12.3)")"   =>",h1,n_equiv(n_ok),ini(n_ok),fin(n_ok),F2_av(n_ok),sig_av(n_ok)
   !call write_info(TRIM(message_text))
  end do  ! fin de la boucle sur les n_ref reflections


  ! seconde boucle sur les reflections pour calculer le R-int
  ns = 0
  F2_sum  = 0.0
  sum_an  = 0.0
  sum_aw  = 0.0
  sum_anw = 0.0

  do i=1, n_ok
   F2_mean = 0.0
   sum_weight = 0.0
   j2 = ini(i)
   IF(n_equiv(i) < 2) cycle
   sig = 0.0
   do j=ini(i), fin(i)
    IF(itreat(j) == j2) then
     sig = sig + 1.0/sig_F2(j)**2     ! Somme des 1/sig**2
    END if
   END do
   sig = 1.0/sig
   do j=ini(i), fin(i)
    if (itreat(j) == j2) then
     ns=ns+1
     F2_sum  = F2_sum  + ABS(F2_av(i)-F2(j))
     sum_an  = sum_an  + F2(j)
     sum_aw  = sum_aw  + sig*((F2_av(i) - F2(j))**2)/sig_F2(j)**2
     sum_anw = sum_anw + sig*(F2(j)/sig_F2(j))**2

     if (merge)  then
      F2_mean    = F2_mean     + F2(j)/1/sig_F2(j)**2
      sum_weight = sum_weight + 1/sig_F2(j)**2
     endif
    end if
   end do
   IF(merge) then
    !F2_mean = F2_mean / n_equiv(i)
    F2_mean = F2_mean /sum_weight

    WRITE(31, '(3I4, 2F8.2)') h(ini(i)), k(ini(i)), l(ini(i)), F2_mean, SQRT(sig)
   endif

  end do
  Rint  = 100.0*F2_sum/MAX(1.0, sum_an)
  Rwint = 100.0*SQRT(sum_aw/MAX(1.0,sum_anw))
  aver_sig = total/REAL(n_ok)
  !n_equiv = ns

  if(len_trim(input_string) == 4) then
   if (input_string(1:4) == 'mono') return
  end if

  !if(.not. input_string(1:4) == 'mono') then
  if(len_trim(input_string_2) /= 6) then
  if(WRITE_details) then
  call write_info('')
  WRITE(message_text, '(a,I8)')   '  . Total number of reflections                        : ', n_ref
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,I8)')   '  . Number of valid independent reflections            : ', n_ok
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,I8)')   '  . Number of valid independent reflections with I>2sig: ', n_ok_2
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,I8)')   '  . Number of valid independent reflections with I>3sig: ', n_ok_3
  call write_info(TRIM(message_text))

  WRITE(message_text, '(a,I8)')   '  . Number of equivalent          reflections          : ', ns
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,I8)')   '  . Number of rejected (absences) reflections          : ', n_rejected
  call write_info(TRIM(message_text))

  WRITE(message_text, '(a,F10.1)') '  . Mean redundancy                                    : ', real(SUM(redundancy(1:n_ok))) / n_ok
  call write_info(TRIM(message_text))
  end if

  IF(Rint < 0.00001) then
   WRITE(message_text, '(a,F6.2,a)') '  . R_int (%)     = ', Rint, '  (HKL data have probably been merged)'
  else
   WRITE(message_text, '(a,F6.2)')   '  . R_int (%)     = ', Rint
  endif
  call write_info(TRIM(message_text))

  if(WRITE_details) then
  WRITE(message_text, '(a,F6.2)')    '  . wR_int (%)    = ', Rwint
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,F8.2)')    '  . <sig> for equiv. reflections   = ', aver_sig
  call write_info(TRIM(message_text))
  end if
  end if

  !call calcul_completude(n_ref, n_ok, input_string_2)
  call calcul_completude(n_ref, n_ok, out_string)
  !WRITE(message_text, '(a,F6.2)') '  . Rw_int = ', Rwint
  !call write_info(TRIM(message_text))
  call write_info('')


  IF(merge) then
   WRITE(message_text, '(2a)') '  >>> merged HKL file: ', TRIM(HKL_file%MERGE)
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,I8)')     '  . Number of valid independent reflections           : ', n_ok
   call write_info(TRIM(message_text))
   CLOSE(UNIT=31)
  endif

  !end if ! fin de la condition if .not. monoclinic


 RETURN
end subroutine calcul_global_Rint


!-------------------------------------------------------
subroutine Get_Friedel_pairs_number

 USE cryscalc_module,             ONLY : SPG, unit_cell, crystal_cell, wavelength, pi, &
                                         keyword_FILE, keyword_CELL, keyword_WAVE, message_text, known_theta, debug_proc
 USE HKL_module
 Use CFML_Reflections_Utilities,  ONLY : HKL_absent,  HKL_equiv, HKL_gen, Hkl_GEN_SXTAL, Get_MaxNumRef, reflect_type
 USE CFML_Crystal_Metrics,        ONLY : Set_Crystal_Cell
 

 USE CFML_Math_General,           ONLY : sort
 USE IO_module,                   ONLY : write_info
 implicit none
  INTEGER                   :: n_ok       ! nombre de reflections conservees
  INTEGER                   :: n_ok_2     ! nombre de reflections conservees avec I>2sig
  INTEGER                   :: n_ok_c, n_ok_nc
  INTEGER                   :: n_ok_2_c, n_ok_2_nc
  INTEGER                   :: nb_Friedel, nb_Friedel_2
  INTEGER                   :: ier
  INTEGER                   :: num_ref
  REAL                      :: theta_min, theta_max
  REAL                      :: STL_min, STL_max
  TYPE(Reflect_type),         ALLOCATABLE, DIMENSION(:) :: reflex_HKL



  INTEGER,  ALLOCATABLE, DIMENSION(:) :: ordered_array

  if(debug_proc%level_2)  call write_debug_proc_level(2, "get_Friedel_pairs_number")

  IF(SPG%NumSpg == 0) then
   call write_info('')
   call write_info('  !! SPGR keyword mandatory for the Rint calculation procedure !!')
   call write_info('')
   return
  endif

  IF(.NOT. keyword_file) then
   call write_info('')
   call write_info('  !! FILE keyword mandatory for the Rint calculation procedure !!')
   call write_info('')
   return
  endif

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

  IF(SPG%centred == 2) then
   call write_info('')
   call write_info('  !! Number of Friedel pairs can be calculated only for acentric space groups !!')
   call write_info('')
   return
  end if

  !
  IF(.NOT. known_theta) call calcul_theta()


  ! Order the reflections by ascending Theta
  if (.not. ordered_hkl) then
   IF(ALLOCATED(ordered_array)) DEALLOCATE(ordered_array)
   ALLOCATE(ordered_array(n_ref), stat=ier)
    if(ier/=0) call write_alloc_error("ordered_array")
   !call write_info('')
   !call write_info('  ... Sort reflections in ascending theta values ...')
   call sort(theta_hkl, n_ref, ordered_array)
   call sort_arrays('theta','+', ordered_array, theta_hkl, 'no_out')
   call write_info('')
   ordered_hkl =.true.
  END if
  !

  write(message_text, '(3a)') '   ... Get Friedel pairs number in ', trim(SPG%SPG_Symb), ' ...'
  call write_info(trim(message_text))
  call write_info('')

  call get_n_ok(n_ok, n_ok_2, .true.)
  n_ok_c   = n_ok
  n_ok_2_c = n_ok_2

  call get_n_ok(n_ok, n_ok_2, .false.)
  n_ok_nc       = n_ok
  n_ok_2_nc     = n_ok_2
  nb_Friedel    = n_ok_nc - n_ok_c
  nb_Friedel_2  = n_ok_2_nc - n_ok_2_c

  
  ! calcul du nombre de Friedel attendu
  theta_min = MINVAL(theta_hkl(1:n_ref))
  theta_max = MAXVAL(theta_hkl(1:n_ref))
  STL_min  = SIN(theta_min  * PI / 180) / wavelength
  STL_max  = SIN(theta_max  * PI / 180) / wavelength
  Num_ref  =  Get_MaxNumRef(STL_max, Unit_cell%volume, STL_min ) 
  IF(ALLOCATED(reflex_HKL))        DEALLOCATE(reflex_HKL)
  ALLOCATE(reflex_HKL(Num_ref), stat=ier)
   if(ier/=0) call write_alloc_error("reflex_HKL")
     
  call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell) 
  call HKL_gen(crystal_cell, SPG, .true., STL_min, STL_max, Num_ref, reflex_HKL )
  
  WRITE(message_text, '(a,I8)')            "  . Number of expected Friedel pairs         : ", Num_ref 
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,I8,a,F6.2,a)')   "  . Number of valid Friedel pairs            : ", nb_Friedel, &
                                           "  (", 100.*nb_Friedel/Num_ref, "%)"
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,I8)')            "  . Number of valid Friedel pairs with I>2sig: ", nb_Friedel_2
  call write_info(TRIM(message_text))

  IF(ALLOCATED(reflex_HKL))        DEALLOCATE(reflex_HKL)
  ALLOCATE(reflex_HKL(Num_ref), stat=ier)


 RETURN
end subroutine Get_Friedel_pairs_number
!-----------------------------------------------------------------------

subroutine Get_n_ok(n_ok, n_ok_2, friedel)

 USE cryscalc_module,             ONLY : SPG, keyword_FILE, keyword_CELL, keyword_WAVE, message_text, known_theta
 USE HKL_module
 Use CFML_Reflections_Utilities,  ONLY : HKL_absent,  HKL_equiv
 USE CFML_Math_General,           ONLY : sort
 USE IO_module,                   ONLY : write_info
 implicit none
  integer, intent(out)          :: n_ok, n_ok_2
  LOGICAL, intent(in)           :: friedel
  INTEGER                       :: i, j

  INTEGER, DIMENSION(n_ref) :: itreat !
  INTEGER, DIMENSION(n_ref) :: ini, fin
  INTEGER, DIMENSION(3)     :: H1, H2
  LOGICAL                   :: absent
  INTEGER                   :: n_rejected ! nombre de reflections rejetees: eteintes par le groupe d'esapce



  n_ok       = 0
  n_ok_2     = 0
  itreat(:)  = 0
  ini(:)     = 0
  fin(:)     = 0

  n_rejected = 0


  do i=1, n_ref    ! boucle sur les n_ref reflections
   IF(itreat(i) /=0) cycle
   H1(1) = h(i)
   H1(2) = k(i)
   H1(3) = l(i)
   absent = hkl_absent(H1, SPG)
   IF(absent) then
    n_rejected = n_rejected + 1
    cycle
   endif

   n_ok = n_ok + 1          ! update the number of independent reflections


   IF(F2(i) > 2*sig_F2(i)) n_ok_2 = n_ok_2 + 1
   itreat(i) = i            ! Make this reflection treated
   ini(n_ok) = i
   fin(n_ok) = i


   do j = i+1, n_ref  ! boucles sur les reflections suivantes: recherche des reflections  equivalentes
    IF(known_theta .and. (ABS(theta_hkl(i) - theta_HKL(j)) > 0.2)) exit
    H2(1) = h(j)
    H2(2) = k(j)
    H2(3) = l(j)
    if (HKL_equiv(H1, H2, SPG, friedel)) then
     itreat(j)     = i
     fin(n_ok)     = j
    end if
   end do
  end do

 return
end subroutine Get_n_ok

!-----------------------------------------------------------------------

subroutine calcul_completude(n_ref, n_ok, input_string)
 USE cryscalc_module,             ONLY : unit_cell, message_text, pi, wavelength, crystal_cell, debug_proc
 Use CFML_Reflections_Utilities,  ONLY : Get_MaxNumRef, HKL_gen, reflect_type
 USE cryscalc_module,             ONLY : SPG
 USE CFML_Crystal_Metrics,        ONLY : Set_Crystal_Cell
 USE HKL_module,                  ONLY : theta_HKL
 USE IO_module,                   ONLY : write_info
 implicit none
  INTEGER,          INTENT(IN)    :: n_ref
  INTEGER,          INTENT(IN)    :: n_ok
  CHARACTER(len=*), intent(in)    :: input_string
  character (len=256)             :: message_text2

  ! variables pour le calcul de la completude
  TYPE(Reflect_type),         ALLOCATABLE, DIMENSION(:) :: reflex_HKL
  REAL                      :: theta_min, theta_max
  REAL                      :: STL_min, STL_max
  INTEGER                   :: num_ref, exp_num_ref
  LOGICAL                   :: Friedel
  INTEGER                   :: ier


  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_completude")

   ! calcul de la completude
  theta_min = MINVAL(theta_hkl(1:n_ref))
  theta_max = MAXVAL(theta_hkl(1:n_ref))
  STL_min  = SIN(theta_min  * PI / 180) / wavelength
  STL_max  = SIN(theta_max  * PI / 180) / wavelength

  Num_ref  =  Get_MaxNumRef(STL_max, Unit_cell%volume, STL_min )
  IF(ALLOCATED(reflex_HKL))        DEALLOCATE(reflex_HKL)
  ALLOCATE(reflex_HKL(Num_ref), stat=ier)
   if(ier/=0) call write_alloc_error("reflex_HKL")

  friedel = .true.
  IF(SPG%centred /= 2)   friedel=.false.

  call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)
  call HKL_gen(crystal_cell, SPG, Friedel, STL_min, STL_max, Num_ref, reflex_HKL )

  !IF(SPG%centred /= 2) then
  ! exp_Num_ref = Num_ref * 2
  !else
  ! exp_Num_ref = num_ref
  !endif
  exp_num_ref = num_ref

  if(len_trim(input_string) /=6) then
   call write_info('')
   WRITE(message_text, '(a,F6.2)')   '  . Theta_min:       ', theta_min
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,F6.2)')   '  . Theta_max:       ', theta_max
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,I8)')     '  . Number of expected independent reflections         : ', exp_Num_ref
   call write_info(TRIM(message_text))
  end if

  if(n_ok < exp_num_ref) then
   write(message_text2, '(I8,a)') exp_num_ref - n_ok, ' missing reflections)'
   write(message_text,  '(a,F6.2,2a)')      '  . Completeness  = ', 100.*(REAL(n_ok)/REAL(exp_Num_ref)),  &
                                            ' (', adjustl(trim(message_text2))
   !WRITE(message_text, '(a,F6.2,a,I8,a)')   '  . Completeness: ', 100.*(REAL(n_ok)/REAL(exp_Num_ref)),  &
   !                                         ' (', exp_num_ref - n_ok, ' missing reflections)'
  else
   WRITE(message_text, '(a,F6.2)')   '  . Completeness  = ', 100.*(REAL(n_ok)/REAL(exp_Num_ref))
  endif
  call write_info(TRIM(message_text))
  
  IF(ALLOCATED(reflex_HKL))        DEALLOCATE(reflex_HKL)
  ALLOCATE(reflex_HKL(Num_ref), stat=ier)


 return
end subroutine calcul_completude

!-----------------------------------------------------------------------
subroutine get_numref_SPG()
 USE hkl_module
 USE IO_module,                   ONLY : write_info
 USE cryscalc_module,             ONLY : message_text, SPG, debug_proc
 Use CFML_Reflections_Utilities,  ONLY : HKL_absent

 implicit none
  INTEGER                     :: i, n_req
  INTEGER, DIMENSION(3)       :: ref_H


  if(debug_proc%level_2)  call write_debug_proc_level(2, "get_numref_SPG")

 n_ref_2 = 0
 n_ref_3 = 0
 n_req   = 0

 do i=1, n_ref
  ref_H(1) = h(i)
  ref_H(2) = k(i)
  ref_H(3) = l(i)

  IF(HKL_absent(ref_H, SPG)) cycle
  n_req = n_req + 1

  if(F2(i) >= 2*sig_F2(i)) n_ref_2 = n_ref_2 + 1
  if(F2(i) >= 3*sig_F2(i)) n_ref_3 = n_ref_3 + 1
 end do


 call write_info('')
 WRITE(message_text, '(3a)') '  Reflections in agreement with ', TRIM(SPG%SPG_Symb), ' space group:'
 call write_info(TRIM(message_text))
 call write_info('')

 write(message_text,'(a,I8)')  '  >> Systematic absence violations:        ', n_ref - n_req
 call write_info(TRIM(message_text))
 write(message_text,'(a,I8)')  '  >> Total number of reflections:          ', n_req
 call write_info(TRIM(message_text))

 IF(n_req == 0) return
 write(message_text,'(a,I8)')  '  >> Number of reflections with I > 2sig.: ', n_ref_2
 call write_info(TRIM(message_text))
 write(message_text,'(a,I8)')  '  >> Number of reflections with I > 3sig.: ', n_ref_3
 call write_info(TRIM(message_text))
 call write_info('')


 return
end subroutine get_numref_SPG

!-------------------------------------------------------
subroutine calcul_merge_HKL()
! extrait de 'datared.F90' JRC
 USE cryscalc_module,             ONLY : SPG, keyword_FILE, keyword_CELL, keyword_WAVE, message_text, known_theta, debug_proc
 USE HKL_module
 Use CFML_Reflections_Utilities,  ONLY : HKL_absent,  HKL_equiv
 USE CFML_Math_General,           ONLY : sort
 USE IO_module,                   ONLY : write_info
 implicit none
  INTEGER                   :: i,j,j2
  INTEGER, DIMENSION(n_ref) :: itreat !
  INTEGER, DIMENSION(n_ref) :: n_equiv, ini, fin
  INTEGER, DIMENSION(3)     :: H1, H2
  LOGICAL                   :: absent
  LOGICAL                   :: friedel
  INTEGER                   :: i1, ier, n_ok
  real                      :: sum_weight
  INTEGER,  ALLOCATABLE, DIMENSION(:) :: ordered_array

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_MERGE_hkl")

  IF(SPG%NumSpg == 0) then
   call write_info('')
   call write_info('  !! SPGR keyword mandatory for the HKL merging procedure !!')
   call write_info('')
   return
  endif

  IF(.NOT. keyword_file) then
   call write_info('')
   call write_info('  !! FILE keyword mandatory for the HKL merging procedure !!')
   call write_info('')
   return
  endif

!  IF(.NOT. keyword_CELL) then
!   call write_info('')
!   call write_info('  !! CELL keyword mandatory for the HKL merging procedure !!')
!   call write_info('')
!   return
!  endif
!
!  IF(.NOT. keyword_WAVE) then
!   call write_info('')
!   call write_info('  !! WAVE keyword mandatory for the HKL merging procedure !!')
!   call write_info('')
!   return
!  END if

  i1 = INDEX(HKL_file%NAME, '.')
  IF(i1 /=0) then
   WRITE(HKL_file%MERGE , '(2a)') HKL_file%NAME(1:i1-1), '_merge.hkl'
  else
   WRITE(HKL_file%MERGE, '(a)') 'merge.hkl'
  endif
  OPEN(UNIT = 31, FILE =TRIM(HKL_file%MERGE))

  ! Order the reflections by ascending Theta
  if (.not. ordered_hkl) then
   IF(ALLOCATED(ordered_array)) DEALLOCATE(ordered_array)
   ALLOCATE(ordered_array(n_ref), stat=ier)
    if(ier/=0) call write_alloc_error("ordered_array")
   call write_info('')
   call write_info('  ... Sort reflections in ascending theta values ...')
   call sort(theta_hkl, n_ref, ordered_array)
   call sort_arrays('theta','+', ordered_array, theta_hkl, 'out')
   call write_info('')
   ordered_hkl =.true.
  END if
  !

  call write_info('   ... Merging the data ...')
  call write_info('')

  friedel = .true.
  IF(SPG%centred /= 2) friedel=.false.

  n_ok       = 0
  itreat(:)  = 0

  !----------------
  do i=1, n_ref
   IF(itreat(i) /=0) cycle
   H1(1) = h(i)
   H1(2) = k(i)
   H1(3) = l(i)
   absent = hkl_absent(H1, SPG)
   IF(absent) cycle

   n_ok = n_ok + 1
   itreat(i) = i
   n_equiv(n_ok) = 1

   ini(n_ok) = i
   fin(n_ok) = i

   !boucles sur les reflections suivantes: . recherche des reflections equivalentes
   !                                       . calcul de la somme des weight
   do j=i+1, n_ref
    IF(known_theta .and. (ABS(theta_hkl(i) - theta_HKL(j)) > 0.2)) exit
    H2(1) = h(j)
    H2(2) = k(j)
    H2(3) = l(j)
    IF(HKL_equiv(H1, H2, SPG, Friedel)) then
     itreat(j)     = i
     fin(n_ok)     = j
     n_equiv(n_ok) = n_equiv(n_ok) + 1
    endif
   end do
  end do


  ! seconde boucle sur les reflections
  do i=1, n_ok
   F2_mean = 0.0
   sum_weight = 0.0
   j2 = ini(i)
   IF(n_equiv(i) < 2) cycle
   do j=ini(i), fin(i)
    if (itreat(j) == j2) then
     F2_mean    = F2_mean    + F2(j)/1/sig_F2(j)**2
     sum_weight = sum_weight + 1/sig_F2(j)**2
    end if
   end do
   F2_mean = F2_mean /sum_weight
   WRITE(31, '(3I4, 2F8.2)') h(ini(i)), k(ini(i)), l(ini(i)), F2_mean, SQRT(1/sum_weight)
  end do

  WRITE(message_text, '(2a)') '  >>> merged HKL file: ', TRIM(HKL_file%MERGE)
  call write_info(TRIM(message_text))
  WRITE(message_text, '(a,I8)')     '  . Number of valid independent reflections            : ', n_ok
  call write_info(TRIM(message_text))

  CLOSE(UNIT=31)
RETURN
end subroutine calcul_merge_HKL

