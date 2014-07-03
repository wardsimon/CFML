!     Last change:  TR   30 Jun 2006    9:15 am
subroutine list_space_groups()
 USE cryscalc_module,           ONLY : list_sg, list_sg_enantio, list_sg_chiral, list_sg_polar,  message_text, debug_proc
 USE IO_module,                 ONLY : write_info

  implicit none
   integer                           :: i1, i2, n_enantio, n_chiral, n_polar

   if(debug_proc%level_2)  call write_debug_proc_level(2, "list_space_groups")
   
      n_enantio = 0
	  n_chiral  = 0
	  n_polar   = 0

      call write_info(' ')
      call write_info('  >> LIST OF SPACE GROUPS:')
      call write_info(' -------------------------')
      call write_info(' ')


      if(list_sg(1)) then          ! triclinique
       i1 = 1
       i2 = 2
       call write_info('  . triclinic:')
       call write_info('')
       call write_space_group(i1, i2, n_enantio, n_chiral, n_polar)
       call write_info('')
      endif

      if(list_sg(2)) then      ! monoclinique
       i1 = 3
       i2 = 15
       call write_info('  . monoclinic:')
       call write_info('')
       call write_space_group(i1, i2, n_enantio, n_chiral, n_polar)
       call write_info('')
      endif

      if(list_sg(3)) then      ! orthorhombique
       i1 = 16
       i2 = 74
       call write_info('  . orthorhombic:')
       call write_info('')
       call write_space_group(i1, i2, n_enantio, n_chiral, n_polar)
       call write_info('')
      endif

      if(list_sg(4)) then      ! tetragonal
       i1 = 75
       i2 = 142
       call write_info('  . tetragonal:')
       call write_info('')
       call write_space_group(i1, i2, n_enantio, n_chiral, n_polar)
       call write_info('')
      endif

      if(list_sg(5)) then      ! trigonal
       i1 = 143
       i2 = 167
       call write_info('  . trigonal:')
       call write_info('')
       call write_space_group(i1, i2, n_enantio, n_chiral, n_polar)
       call write_info('')
      endif

      if(list_sg(6)) then      ! hexagonal
       i1 = 168
       i2 = 194
       call write_info('  . hexagonal:')
       call write_info('')
       call write_space_group(i1, i2, n_enantio, n_chiral, n_polar)
       call write_info('')
      endif

      if(list_sg(7)) then      ! cubique
       i1 = 195
       i2 = 230
       call write_info('  . cubic:')
       call write_info('')
       call write_space_group(i1, i2, n_enantio, n_chiral, n_polar)
       call write_info('')
      endif

	  
	  if(n_enantio /=0) then
	   if(.not. list_sg_enantio) then
        call write_info(' ')
        call write_info('     (* enantiomorphic space group)')	    
       endif
	  end if
	  
	  if(n_chiral/=0) then
	   if(.not. list_sg_chiral) then
	    call write_info('     (+ chiral space group)')
	   end if
	  end if 
	  
	  if(n_polar /=0) then
	   if(.not. list_sg_polar) then
	    call write_info('     (§ polar space group)')
	    call write_info(' ') 
	   end if
	  end if 
	  
	 return 
end subroutine list_space_groups
!-------------------------------------------------------------
subroutine test_Bravais(space_group, ok)
 USE cryscalc_module, ONLY : list_sg_bravais

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
 USE cryscalc_module, ONLY : list_sg_laue
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
 USE cryscalc_module, ONLY : laue_class, message_text
 USE IO_module,       ONLY : write_info
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
subroutine test_enantio(i, ok)
 implicit none
  integer, intent(in)    :: i
  logical, intent(inout) :: ok
  
  ok = .false.
  if(i ==  76 .or. i ==  78 .or. &       ! P 41         P 42
	 i ==  91 .or. i ==  95 .or. &       ! P 41 2  2    P 43 2  2
	 i ==  92 .or. i ==  96 .or. &       ! P 41 21 2    P 43 21 2
	 i == 144 .or. i == 145 .or. &       ! P 31         P 32
     i == 151 .or. i == 153 .or. &       ! P 31 1  2    P 32 1  2
	 i == 152 .or. i == 154 .or. &       ! P 31 2  1    P 32 2  1
	 i == 169 .or. i == 170 .or. &       ! P 61         P 65
	 i == 171 .or. i == 172 .or. &       ! P 62         P 64
	 i == 178 .or. i == 179 .or. &       ! P 61 2  2    P 65 2  2
	 i == 180 .or. i == 181 .or. &	     ! P 62 2  2    P 64 2  2
	 i == 213 .or. i == 212 )    &       ! P 41 3  2    P 43 3  2
	 then
	  ok = .true.
  end if
	
 return
end subroutine test_enantio
!-------------------------------------------------------------

subroutine test_chiral(i, ok)
! www.ruppweb.org/xray/tutorial/enantion.htm

 implicit none
  integer, intent(in)    :: i
  logical, intent(inout) :: ok
  
  ok = .false.
  if(i ==   1 .or.                                                          & ! P 1
     i ==   2 .or. i ==   4 .or. i ==   5 .or.                              & ! P 2       P 21      C 2
	 i ==  16 .or. i ==  17 .or. i ==  18 .or. i ==  19 .or. i == 20  .or.  & ! P 2 2 2   P 2 2 21  P 21 21 2   P 21 21 21  C 2 2 21
	 i ==  21 .or. i ==  22 .or. i ==  23 .or. i ==  24 .or.                & ! C 2 2 2   F 2 2 2   I 2  2  2   I 21 21 21
	 i ==  75 .or. i ==  76 .or. i ==  77 .or. i ==  78 .or. i == 79  .or.  & ! P 4       P 41      P 42        P 43        I 4
	 i ==  80 .or. i ==  89 .or.                                            & ! I 41      P 4 2 2
	 i ==  90 .or. i ==  91 .or. i ==  92 .or. i ==  93 .or. i == 94  .or.  & ! P 4 21 2  P 41 2 2  P 41 21 2   P 42 2 2    P 42 21 2
	 i ==  95 .or. i ==  96 .or. i ==  97 .or. i ==  98 .or.                & ! P 43 2 2  P 43 21 2 I 4 2 2     I 41 2 2  
	 i == 143 .or. i == 144 .or. i == 145 .or. i == 146 .or.                & ! P 3       P 31      P 32        R 3
	 i == 149 .or. i == 150 .or. i == 151 .or. i == 152 .or.                & ! P 3 1 2   P 3 2 1   P 31 1 2    P 31 2 1
	 i == 153 .or. i == 154 .or. i == 155 .or.                              & ! P 32 1 2  P 32 2 1  R 3 2
     i == 168 .or. i == 169 .or. i == 170 .or. i == 171 .or. i == 172 .or.  & ! P 6       P 61      P 65        P 62        P 64	 
	 i == 173 .or. i == 177 .or.                                            & ! P 63      P 6 2 2
	 i == 178 .or. i == 179 .or. i == 180 .or. i == 181 .or. i == 182 .or.  & ! P 61 2 2  P 65 2 2  P 62 2 2    P 64 2 2    P 63 2 2
	 i == 195 .or. i == 196 .or. i == 197 .or. i == 198 .or. i == 199 .or.  & ! P 2 3     F 2 3     I 2 3       P 21 3      I 21 3
	 i == 207 .or. i == 208 .or. i == 209 .or.                              & ! P 4 3 2   I 4 3 2   F 4 3 2
	 i == 210 .or. i == 211 .or. i == 212 .or. i == 213 .or. i == 214)      & ! F 41 3 2  I 4 3 2   P 43 3 2    P 41 3 2    I 41 3 2
	then
	  ok = .true.
  end if
	

 return
end subroutine test_chiral

!-------------------------------------------------------------

subroutine test_polar(i, ok)
! www.ruppweb.org/xray/tutorial/enantion.htm

 implicit none
  integer, intent(in)    :: i
  logical, intent(inout) :: ok
  
  ok = .false.
  if(i ==   1 .or.                                                          & ! P 1
     i ==   2 .or. i ==   4 .or. i ==   5 .or.                              & ! P 2       P 21      C 2
	 
	 i ==  75 .or. i ==  76 .or. i ==  77 .or. i ==  78 .or. i == 79  .or.  & ! P 4       P 41      P 42        P 43        I 4
	 i ==  80 .or.                                                          & ! I 41      
	 
	 i == 143 .or. i == 144 .or. i == 145 .or. i == 146 .or.                & ! P 3       P 31      P 32        R 3
	 
	 i == 168 .or. i == 169 .or. i == 170 .or. i == 171 .or. i == 172 .or.  & ! P 6       P 61      P 65        P 62        P 64	 
	 i == 173 )                                                             & ! P 63      
	 
	
	then
	  ok = .true.
  end if
	

 return
end subroutine test_polar

!-------------------------------------------------------------
subroutine test_hkl_equiv_condition(i)
 USE IO_module,                 ONLY : write_info
 implicit none
 integer,  intent(in)            :: i
 character (len=256)             :: hkl_condition_string
  
 hkl_condition_string = ''
 
 if(i >= 195 .and. i <= 206) then
  hkl_condition_string = '   > h,k,l cyclically permutable !!'
 elseif(i>=207) then 
  hkl_condition_string = '   > h,k,l permutable.'
 else
  return
 end if
 
 call write_info('')
 call write_info(trim(hkl_condition_string))
 call write_info('')
 
 return
end subroutine test_hkl_equiv_condition
!-------------------------------------------------------------


subroutine Get_MOVE_string(MOVE_string)
use cryscalc_module,                 ONLY :  SPG, debug_proc
 implicit none
   
  character (len=32), intent(inout) :: MOVE_string
  character (len=6)                 :: m_a, m_b, m_c
  
  if(debug_proc%level_2)  call write_debug_proc_level(2, "get_MOVE_string")
   
   
  if(SPG%NumSPG == 43) then                                 ! F d d 2
   m_a = '  0.25'
   m_b = '  0.25'
   m_c = '  1.00'
  elseif(SPG%NumSPG == 80) then                             ! I 41 
   m_a = '  1.00'
   m_b = '  0.50'
   m_c = '  1.00'
  elseif(SPG%NumSPG == 98) then                             ! I 41 2 2  
   m_a = '  1.00'
   m_b = '  0.50'
   m_c = '  0.25'
  elseif(SPG%NumSPG == 109) then                            ! I 41 m d
   m_a = '  1.00'
   m_b = '  0.50'
   m_c = '  1.00'
  elseif(SPG%NumSPG == 110) then                            ! I 41 c d  
   m_a = '  1.00'
   m_b = '  0.50'
   m_c = '  1.00'
  elseif(SPG%NumSPG == 122) then                            ! I -4 2 d  
   m_a = '  1.00'
   m_b = '  0.50'
   m_c = '  0.25'
  elseif(SPG%NumSPG == 210) then                            ! F 41 3 2  
   m_a = '  0.25'
   m_b = '  0.25'
   m_c = '  0.25'
  else
   m_a = '  1.00'
   m_b = '  1.00'
   m_c = '  1.00'
  endif
     
  if(SPG%info(1:3) == "bca") then
   write(MOVE_string, '(5a)') 'MOVE   ', m_b, m_c, m_a, '  -1'
  elseif(SPG%info(1:3) == "cab") then
   write(MOVE_string, '(5a)') 'MOVE   ', m_c, m_a, m_b, '  -1'
  else	
   write(MOVE_string, '(5a)') 'MOVE   ', m_a, m_b, m_c, '  -1'
  end if
   
  
return
end subroutine Get_MOVE_string
