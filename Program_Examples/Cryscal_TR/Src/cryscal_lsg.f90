!     Last change:  TR   27 Jun 2006    6:28 pm
subroutine list_space_groups()
 USE list_space_groups_module, ONLY : sg, definition_space_groups
 USE cryscal_module,           ONLY : list_sg, list_sg_centric, message_text, known_space_groups, debug_proc
 USE IO_module,                ONLY : write_info


  implicit none
   character (len=16), dimension(15) :: sg2
   integer                           :: i, i1, i2, j, n_ok
   LOGICAL                           :: ok

   if(debug_proc%level_2)  call write_debug_proc_level(2, "list_space_groups")
   
   IF(.NOT. known_space_groups) then   
    call definition_space_groups()
    known_space_groups = .true.
   endif


      call write_info(' ')
      call write_info('  >> LIST OF SPACE GROUPS:')
      call write_info(' -------------------------')
      call write_info(' ')

      n_ok = 0
      if(list_sg(1)) then          ! triclinique
       i1 = 1
       i2 = 2
       call write_info('  . triclinic:')
       call write_info('')
       do i = i1, i2
        j = INDEX(sg(i),'(')
        IF(list_sg_centric(1) .AND. sg(i)(j+1:j+1)=='n') cycle
        IF(list_sg_centric(2) .AND. sg(i)(j+1:j+1)=='c') cycle

        call test_Bravais(sg(i), ok)
        IF(.NOT. ok) cycle
        
        n_ok = n_ok + 1 
        WRITE(message_text, '(10x,I3,5x, a,I1,7x,a)') n_ok, 'IT# ', i,TRIM(sg(i))
        call write_info(TRIM(message_text))
       end do
       call write_info('')
      endif

      if(list_sg(2)) then      ! monoclinique
       i1 = 3
       i2 = 15
       call write_info('  . monoclinic:')
       call write_info('')
       do i = i1, i2
        j = INDEX(sg(i),'(')
        IF(list_sg_centric(1) .AND. sg(i)(j+1:j+1)=='n') cycle
        IF(list_sg_centric(2) .AND. sg(i)(j+1:j+1)=='c') cycle

        call test_Bravais(sg(i), ok)
        IF(.NOT. ok) cycle
 
        n_ok = n_ok + 1
        if(i<10) then
         WRITE(message_text, '(10x,I3,5x,a,I1,7x,a)') n_ok, 'IT# ',i,TRIM(sg(i))
        else
         WRITE(message_text, '(10x,I3,5x,a,I2,6x,a)') n_ok, 'IT# ',i,TRIM(sg(i))
        endif
        call write_info(TRIM(message_text))
       end do
       call write_info('')
      endif

      if(list_sg(3)) then      ! orthorhombique
       i1 = 16
       i2 = 74
       call write_info('  . orthorhombic:')
       call write_info('')
       do i = i1, i2
        j = INDEX(sg(i),'(')
        IF(list_sg_centric(1) .AND. sg(i)(j+1:j+1)=='n') cycle
        IF(list_sg_centric(2) .AND. sg(i)(j+1:j+1)=='c') cycle

        call test_Bravais(sg(i), ok)
        IF(.NOT. ok) cycle

        n_ok = n_ok + 1
        WRITE(message_text, '(10x,I3,5x,I2,6x,a)') n_ok, 'IT# ', i,TRIM(sg(i))
        call write_info(TRIM(message_text))
       end do
       call write_info('')
      endif

      if(list_sg(4)) then      ! tetragonal
       i1 = 75
       i2 = 142
       call write_info('  . tetragonal:')
       call write_info('')
       do i = i1, i2
        j = INDEX(sg(i),'(')
        IF(list_sg_centric(1) .AND. sg(i)(j+1:j+1)=='n') cycle
        IF(list_sg_centric(2) .AND. sg(i)(j+1:j+1)=='c') cycle
        call test_Bravais(sg(i), ok)
        IF(.NOT. ok) cycle
        
        n_ok = n_ok + 1
        if(i < 100) then
         WRITE(message_text, '(10x,I3,5x,I2,6x,a)') n_ok, 'IT# ', i,TRIM(sg(i))
        else
         WRITE(message_text, '(10x,I3,5x,I3,5x,a)') n_ok, 'IT# ', i,TRIM(sg(i))
        endif
        call write_info(TRIM(message_text))
       end do
       call write_info('')
      endif

      if(list_sg(5)) then      ! trigonal
       i1 = 143
       i2 = 167
       call write_info('  . trigonal:')
       call write_info('')
       do i = i1, i2
        j = INDEX(sg(i),'(')
        IF(list_sg_centric(1) .AND. sg(i)(j+1:j+1)=='n') cycle
        IF(list_sg_centric(2) .AND. sg(i)(j+1:j+1)=='c') cycle
        call test_Bravais(sg(i), ok)
        IF(.NOT. ok) cycle

        n_ok = n_ok + 1
        WRITE(message_text, '(10x,I3,5x,a,I3,5x,a)') n_ok, 'IT# ', i,TRIM(sg(i))
        call write_info(TRIM(message_text))
       end do
       call write_info('')
      endif

      if(list_sg(6)) then      ! hexagonal
       i1 = 168
       i2 = 194
       call write_info('  . hexagonal:')
       call write_info('')
       do i = i1, i2
        j = INDEX(sg(i),'(')
        IF(list_sg_centric(1) .AND. sg(i)(j+1:j+1)=='n') cycle
        IF(list_sg_centric(2) .AND. sg(i)(j+1:j+1)=='c') cycle
        call test_Bravais(sg(i), ok)
        IF(.NOT. ok) cycle

        n_ok = n_ok + 1
        WRITE(message_text, '(10x,I3,5x,a,I3,5x,a)') n_ok, 'IT# ', i,TRIM(sg(i))
        call write_info(TRIM(message_text))
       end do
       call write_info('')
      endif

      if(list_sg(7)) then      ! cubique
       i1 = 195
       i2 = 230
       call write_info('  . cubic:')
       call write_info('')
       do i = i1, i2
        j = INDEX(sg(i),'(')
        IF(list_sg_centric(1) .AND. sg(i)(j+1:j+1)=='n') cycle
        IF(list_sg_centric(2) .AND. sg(i)(j+1:j+1)=='c') cycle
        call test_Bravais(sg(i), ok)
        IF(.NOT. ok) cycle

        n_ok = n_ok + 1
        WRITE(message_text, '(10x,I3,5x,a,I3,5x,a)') n_ok, 'IT# ', i,TRIM(sg(i))
        call write_info(TRIM(message_text))
       end do
       call write_info('')
      endif


 return
end subroutine list_space_groups

!-------------------------------------------------------------
subroutine test_Bravais(space_group, ok)
 USE cryscal_module, ONLY : list_sg_bravais, debug_proc

 implicit none
  CHARACTER (LEN=*), INTENT(IN)    :: space_group
  LOGICAL,           INTENT(INOUT) :: ok
  INTEGER                          :: i

  if(debug_proc%level_2)  call write_debug_proc_level(2, "test_Bravais")
  ok=.false.

  IF(list_sg_bravais(1) .and. space_group(1:1) == "P") ok=.true.
  IF(list_sg_bravais(2) .and. space_group(1:1) == "A") ok=.true.
  IF(list_sg_bravais(3) .and. space_group(1:1) == "B") ok=.true.
  IF(list_sg_bravais(4) .and. space_group(1:1) == "C") ok=.true.
  IF(list_sg_bravais(5) .and. space_group(1:1) == "I") ok=.true.
  IF(list_sg_bravais(6) .and. space_group(1:1) == "F") ok=.true.
  IF(list_sg_bravais(7) .and. space_group(1:1) == "R") ok=.true.

  RETURN
end subroutine test_Bravais

