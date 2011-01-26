!     Last change:  TR   30 Jun 2006    9:15 am
subroutine list_space_groups()
 USE cryscal_module,            ONLY : list_sg
 USE IO_module,                 ONLY : write_info

  implicit none
   integer                           :: i1, i2


      call write_info(' ')
      call write_info('  >> LIST OF SPACE GROUPS:')
      call write_info(' -------------------------')
      call write_info(' ')


      if(list_sg(1)) then          ! triclinique
       i1 = 1
       i2 = 2
       call write_info('  . triclinic:')
       call write_info('')
       call write_space_group(i1, i2)
       call write_info('')
      endif

      if(list_sg(2)) then      ! monoclinique
       i1 = 3
       i2 = 15
       call write_info('  . monoclinic:')
       call write_info('')
       call write_space_group(i1, i2)
       call write_info('')
      endif

      if(list_sg(3)) then      ! orthorhombique
       i1 = 16
       i2 = 74
       call write_info('  . orthorhombic:')
       call write_info('')
       call write_space_group(i1, i2)
       call write_info('')
      endif

      if(list_sg(4)) then      ! tetragonal
       i1 = 75
       i2 = 142
       call write_info('  . tetragonal:')
       call write_info('')
       call write_space_group(i1, i2)
       call write_info('')
      endif

      if(list_sg(5)) then      ! trigonal
       i1 = 143
       i2 = 167
       call write_info('  . trigonal:')
       call write_info('')
       call write_space_group(i1, i2)
       call write_info('')
      endif

      if(list_sg(6)) then      ! hexagonal
       i1 = 168
       i2 = 194
       call write_info('  . hexagonal:')
       call write_info('')
       call write_space_group(i1, i2)
       call write_info('')
      endif

      if(list_sg(7)) then      ! cubique
       i1 = 195
       i2 = 230
       call write_info('  . cubic:')
       call write_info('')
       call write_space_group(i1, i2)
       call write_info('')
      endif



end subroutine list_space_groups
!-------------------------------------------------------------
subroutine test_Bravais(space_group, ok)
 USE cryscal_module, ONLY : list_sg_bravais

 implicit none
  CHARACTER (LEN=*), INTENT(IN)    :: space_group
  LOGICAL,           INTENT(INOUT) :: ok
  INTEGER                          :: i

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
!-------------------------------------------------------------
subroutine test_laue(laue, ok)
 USE cryscal_module, ONLY : list_sg_laue
  implicit none
   CHARACTER (LEN=*), INTENT(IN)    :: laue
   LOGICAL,           INTENT(INOUT) :: ok
   INTEGER                          :: i, long

   long = LEN_TRIM(laue)
   ok = .false.
   IF(list_sg_laue(1)  .AND. laue(1:long) == "-1")     ok = .true.
   IF(list_sg_laue(2)  .AND. laue(1:long) == "2/m")    ok = .true.
   IF(list_sg_laue(3)  .AND. laue(1:long) == "mmm")    ok = .true.
   IF(list_sg_laue(4)  .AND. laue(1:long) == "4/m")    ok = .true.
   IF(list_sg_laue(5)  .AND. laue(1:long) == "4/mmm")  ok = .true.
   IF(list_sg_laue(6)  .AND. laue(1:long) == "-3")     ok = .true.
   IF(list_sg_laue(7)  .AND. laue(1:long) == "-3")     ok = .true.
   IF(list_sg_laue(8)  .AND. (laue(1:long) == "-31m" .or. laue(1:long) == "-3m1"))    ok = .true.
   IF(list_sg_laue(9)  .AND. laue(1:long) == "-31m")   ok = .true.
   IF(list_sg_laue(10) .AND. laue(1:long) == "-3m1")   ok = .true.
   IF(list_sg_laue(11) .AND. laue(1:long) == "6/m")    ok = .true.
   IF(list_sg_laue(12) .AND. laue(1:long) == "6/mmm")  ok = .true.
   IF(list_sg_laue(13) .AND. laue(1:long) == "m-3")    ok = .true.
   IF(list_sg_laue(14) .AND. laue(1:long) == "m-3m")   ok = .true.

  RETURN
end subroutine test_laue
!-------------------------------------------------------------

subroutine list_laue_class()
 USE cryscal_module, ONLY : laue_class, message_text
 USE IO_module,      ONLY : write_info
 implicit none
  INTEGER :: i

 call def_laue_class()
 
 call write_info('          Crystal system           Laue class')
 call write_info('')               
 do i=1, 14
  WRITE(message_text, '(3x,a,i2,a,3x,a)') '[',i,']', TRIM(laue_class(i))
  call write_info(TRIM(message_text))
 end do

 return

end subroutine list_laue_class

!-------------------------------------------------------------
