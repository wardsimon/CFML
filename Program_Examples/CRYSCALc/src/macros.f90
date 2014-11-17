!     Last change:  TR   29 Mar 2007   12:31 pm
!  definition de fonctions generales

MODULE MACROS_module
 implicit none

! oct. 09:  . routine remove_car  remplacee par function remove_car
!           . routine remove_line remplacee par function clear_string
!           . routine replace_car remplacee par function replace_car


 PRIVATE
 PUBLIC    :: Nombre_de_colonnes,        &
              !lower_case,               &
              !upper_case,               &
              l_case,                    &
              u_case,                    &
              test_file_exist,           &
              error_message,             &
              signe,                     &
              remove_car,                &
              clear_string,              &
              replace_car,               &
              replace_car1,              &
              replace_car2,              &
			  verif_hkl_format,          &
			  multiple,                  &
			  answer_yes,                &
			  answer_no,                 &
              write_message,             &
              check_character,           &
			  Get_current_folder_name,   &
			  Get_sample_ID_from_folder_name, &
              Get_Wingx_job


  interface signe
   module PROCEDURE signe_real
   module PROCEDURE signe_int
  end interface signe

 contains

!--------------------------------------------------------------------
! determination du nombre de colonnes dans une ligne

 subroutine Nombre_de_colonnes(line, nb_col)
  CHARACTER (LEN=*),   INTENT(IN)    :: line
  INTEGER,             INTENT(OUT)   :: nb_col
  CHARACTER (LEN=256000)             :: adjusted
  INTEGER                            :: i, ini

  IF(LEN_TRIM(line) == 0) then
   nb_col = 0
   return
  endif

  adjusted = adjustl(line)       !ajustement de la chaine a gauche
  nb_col = 1

  ini = INDEX(adjusted, '"')
  IF(ini==0) then
   ini=1
  else
   adjusted = adjusted(ini+1:)
   ini = INDEX(adjusted, '"') + 1
  endif

  do i = ini, LEN_TRIM(adjusted)
   if (adjusted(i:i) /= ' ') then
    if (i /= 1) then
     if (adjusted(i-1:i-1) == ' ') then
      nb_col = nb_col + 1
     end if
    end if
   end if
  end do

 end subroutine Nombre_de_colonnes

!------------------------------------------------
!    Subroutine Lower_Case(String)
!    !---- Convert string to Lowercase
!      Character (len=*), intent(in out) :: String
!      Integer                           :: i,Ln,iCar
!
!      Ln = Len(String)
!      do i=1,Ln
!        iCar = iChar(String(i:i))
!        if (iCar >= 65 .and. iCar <= 90 ) String(i:i) = Char(iCar+32)
!      end do
!      return
!    end Subroutine Lower_Case
!
!!------------------------------------------------
!   Subroutine Upper_Case(String)
!     !---- Convert string to uppercase
!      Character (len=*), intent(in out) :: String
!      Integer                           ::   i,Ln,iCar
!!
!      Ln = Len(String)
!      do i=1,Ln
!        iCar = iChar(String(i:i))
!        if (iCar >= 97 .and. iCar <= 122) String(i:i) = Char(iCar-32)
!      end do
!!
!      return
!      end Subroutine Upper_Case


!!----
    !!---- Character Function L_Case(Text) Result (Mtext)
    !!----    character (len=*), intent(in) :: text   !  In -> String: "InPUT Line"
    !!----    character (len=len(text))     :: mtex   ! Out -> String: "input line"
    !!----
    !!----    Conversion to lower case, text is not modified
    !!----
    !!---- Update: February - 2003
    !!
    Function L_Case(Text) Result (Mtext)
       !---- Argument ----!
       character (len=*), intent(in) :: text
       character (len=len(text))     :: mtext

       !---- Local variables ----!
       integer, parameter :: inc = ICHAR("A") - ICHAR("a")
       integer            :: leng, pos

       mtext=text
       leng=len_trim(mtext)
       do pos=1,leng
          if (mtext(pos:pos) >= "A" .and. mtext(pos:pos) <= "Z")           &
              mtext(pos:pos) = CHAR ( ICHAR(mtext(pos:pos)) - inc )
       end do

       return
    End Function L_Case

    !!----
    !!---- Character Function U_Case(Text) Result (Mtext)
    !!----    character (len=*), intent(in) :: text   !  In -> String:"Input Line"
    !!----    character (len=len(text))     :: mtext  ! Out -> String:"INPUT LINE"
    !!----
    !!----    Conversion to upper case, text is not modified
    !!----
    !!---- Update: February - 2003
    !!
    Function U_Case(Text) Result (Mtext)
       !---- Argument ----!
       character (len=*), intent(in) :: text
       character (len=len(text))     :: mtext

       !---- Local variables ----!
       integer, parameter :: inc = ICHAR("A") - ICHAR("a")
       integer            :: leng, pos

       mtext=text
       leng=len_trim(mtext)
       do pos=1,leng
          if (mtext(pos:pos) >= "a" .and. mtext(pos:pos) <= "z")           &
              mtext(pos:pos) = CHAR ( ICHAR(mtext(pos:pos)) + inc )
       end do

       return
    End Function U_Case


!-------------------------------------------------------------------

 subroutine test_file_exist(file_name, file_exist, input_string)
  use IO_module
  CHARACTER(LEN=*), INTENT(IN)          :: file_name
  LOGICAL,          INTENT(INOUT)       :: file_exist
  CHARACTER(LEN=*), INTENT(IN)          :: input_string
  character(len=256)                    :: text_info

  inquire (FILE=TRIM(file_name), EXIST=file_exist)
 ! if (.not. file_exist) then
 !  write(*,'(a)') ' '
 !  write(*,'(a)') '  >>> '// trim(file_name) //' does not exist !'
 !  WRITE(*,'(a)') ' '
 !  write(2,'(a)') ' '
 !  write(2,'(a)') '  >>> '// trim(file_name) //' does not exist !'
 !  WRITE(2,'(a)') ' '
!
!  ! stop
!  end if


  !if(len_trim(file_name) == 0) then
  ! call write_info('')
  ! call write_info(' archive.cif file not defined !')
  ! call write_info('')
  ! return
  !endif
  
  if(.not. file_exist .and. input_string(1:3) == 'out') then
   call write_info('')
   write(text_info, '(3a)') ' >>> ', trim(file_name), ' does not exist !'
   call write_info(trim(text_info))
   call write_info('')
  endif

  return
 end subroutine test_file_exist
!-------------------------------------------------------------------
!
!
subroutine error_message(error_string)
 CHARACTER (LEN=*), INTENT(IN)    :: error_string

 WRITE(*,'(a)' ) ' '
 WRITE(*,'(3a)') '   >> Error reading file at line: ',TRIM(error_string), ' !'
 WRITE(*,'(a)' ) ' '
 WRITE(2,'(a)' ) ' '
 WRITE(2,'(3a)') '   >> Error reading file at line: ',TRIM(error_string), ' !'
 WRITE(2,'(a)' ) ' '
 !stop
 return

end subroutine error_message



!-------------------------------------------------------------------

 elemental function signe_real(a) RESULT (signe_a)
  REAL, INTENT(IN) :: a
  REAL             :: signe_a

  IF(a < 0.) then
   signe_a = -1
  else
   signe_a = 1
  endif

 end function signe_real

 elemental function signe_int(a) RESULT (signe_a)
  integer, INTENT(IN) :: a
  REAL                :: signe_a

  IF(a < 0.) then
   signe_a = -1
  else
   signe_a = 1
  endif

 end function signe_int


!-------------------------------------------------------------------
!subroutine remove_car(chaine, caractere)
! ! enleve un caractere (ou chaine de caractere) d'une chaine
!
! implicit none
!  CHARACTER (LEN=*) ,   INTENT(IN OUT):: chaine
!  CHARACTER (LEN=*),    INTENT(IN)    :: caractere
!  INTEGER                             :: i, len_car
!
!  len_car = LEN_TRIM(caractere)
!  IF(len_car==0) len_car = 1
!
!  do
!   i = INDEX(chaine,caractere)
!   if (i==0) exit
!   IF (i> LEN_TRIM(chaine)) exit
!
!   if (i==1) then
!    chaine = chaine(len_car+1:LEN_TRIM(chaine))
!   elseif (i==LEN_TRIM(chaine)) then
!    chaine = chaine(1:i-1)
!   else
!    chaine = chaine(1:i-1) // chaine(i+len_car:)
!   endif
!  END do
!
! return
!end subroutine remove_car


!-------------------------------------------------------------------
 function remove_car(chaine, caractere) result (new_chaine)
  ! enleve un caractere (ou chaine de caractere) d'une chaine

  implicit none
   CHARACTER (LEN=*) ,   INTENT(IN)    :: chaine
   CHARACTER (LEN=*),    INTENT(IN)    :: caractere
   CHARACTER (len=len(chaine))         :: new_chaine

   INTEGER                             :: i, len_car

   len_car = LEN_TRIM(caractere)
   IF(len_car==0) len_car = 1

   new_chaine = chaine
   do
    i = INDEX(new_chaine,caractere)
    if (i==0) exit
    IF (i> LEN_TRIM(new_chaine)) exit

    if (i==1) then
     new_chaine = new_chaine(len_car+1:LEN_TRIM(new_chaine))
    elseif (i==LEN_TRIM(new_chaine)) then
     new_chaine = new_chaine(1:i-1)
    else
     new_chaine = new_chaine(1:i-1) // new_chaine(i+len_car:)
    endif
   END do

  return
 end function remove_car

!-------------------------------------------------------------------
! subroutine clear_string(chaine, caractere)
!  ! la chaine de caractere 'chaine' est remplacee par une chaine vide si elle contient
!  ! la chaine de caractere 'caractere'
!  ! ancien nom de la routine : remove_line
!
!  implicit none
!   CHARACTER (LEN=*) ,   INTENT(IN OUT):: chaine
!   CHARACTER (LEN=*),    INTENT(IN)    :: caractere
!   INTEGER                             :: i


!   i = INDEX(chaine,caractere)
!   if (i/=0) chaine = ''

!  return
! end subroutine clear_string

!-------------------------------------------------------------------
 function clear_string(chaine, caractere) result(new_chaine)
  ! la chaine de caractere 'chaine' est remplacee par une chaine vide si elle contient
  ! la chaine de caractere 'caractere'

  implicit none
   CHARACTER (LEN=*) ,   INTENT(IN)    :: chaine
   CHARACTER (LEN=*),    INTENT(IN)    :: caractere
   CHARACTER (LEN=len(chaine))         :: new_chaine
   INTEGER                             :: i

   new_chaine = chaine
   i = INDEX(new_chaine,caractere)
   if (i/=0) new_chaine = ''

  return
 end function clear_string

!-------------------------------------------------------------------
! subroutine replace_car3(chaine, car1, car2)
!  ! remplace un caractere (ou une chaine de caractere) par un autre (ou par une autre chaine de caractere)
!  ! dans une chaine

!  implicit none
!   CHARACTER (LEN=*),    INTENT(IN OUT):: chaine
!   CHARACTER (LEN=*),    INTENT(IN)    :: car1, car2
!   INTEGER                             :: i
!   integer                             :: len_car1


!   len_car1 = len_trim(car1)
!
!  do
!    i = INDEX(chaine,car1)
!
!    if (i==0) exit
!    IF (i> LEN_TRIM(chaine)) exit
!    chaine = chaine(1:i-1)//car2//chaine(i+len_car1:)
!   END do


!  return
! end subroutine replace_car3
 !-------------------------------------------------------------------
 function replace_car(chaine, car1, car2) result(new_chaine)
  ! remplace un caractere (ou une chaine de caractere) par un autre (ou par une autre chaine de caractere)
  ! dans une chaine

  implicit none
   CHARACTER (LEN=*),    INTENT(IN)    :: chaine
   CHARACTER (LEN=*),    INTENT(IN)    :: car1, car2
   CHARACTER (LEN=len(chaine))         :: new_chaine
   character (len=1)                   :: new_car
   INTEGER                             :: i
   integer                             :: len_car1
   logical                             :: modif_car2


   new_chaine = chaine
   !len_car1 = len_trim(car1)

   !do
   ! i = INDEX(new_chaine,car1)

   ! if (i==0) exit
   ! IF (i> LEN_TRIM(new_chaine)) exit
   ! new_chaine = new_chaine(1:i-1)//car2//new_chaine(i+len_car1:)
   !END do

   modif_car2 = .false.
   do i = 1, len_trim(car2)
    if(index(car1, car2(i:i)) /=0) then	  
     if(index(chaine, '/') /=0) then
	  new_car = "/"
	 else
      if(index(chaine, '/') /=0) then
	   new_car ="!"
      endif	  
     endif	 
	 new_chaine = replace_car2(chaine, car1(1:1), new_car)
	 new_chaine = replace_car2(new_chaine, new_car, car2)
	 
	 modif_car2 = .true.
	 exit
    end if
   end do	
   
   if(.not. modif_car2)  new_chaine =  replace_car2(chaine, car1, car2)

   return
 end function replace_car
 !-------------------------------------------------------------------
 function replace_car1(chaine, car1, car2) result(new_chaine)
  ! remplace un caractere (ou une chaine de caractere) par un autre (ou par une autre chaine de caractere)
  ! dans une chaine

  implicit none
   CHARACTER (LEN=*),    INTENT(IN)    :: chaine
   CHARACTER (LEN=*),    INTENT(IN)    :: car1, car2
   CHARACTER (LEN=len(chaine))         :: new_chaine
   character (len=1)                   :: new_car
   INTEGER                             :: i
   integer                             :: len_car1
   logical                             :: modif_car2


   new_chaine = chaine
   modif_car2 = .false.
   do i = 1, len_trim(car2)
    if(index(car1, car2(i:i)) /=0) then	  
     !if(index(chaine, '/') /=0) then
	 ! new_car = "/"
	 !else
     ! if(index(chaine, '/') /=0) then
	 !  new_car ="!"
     ! endif	  
     !endif	 
	 new_chaine = replace_car2(chaine, car1(1:1), new_car)
	 new_chaine = replace_car2(new_chaine, new_car, car2)
	 
	 modif_car2 = .true.
	 exit
    end if
   end do
   
   
   if(.not. modif_car2)  new_chaine =  replace_car2(chaine, car1, car2)

   return
 end function replace_car1
 
  !-------------------------------------------------------------------
 function replace_car2(chaine, car1, car2) result(new_chaine)
  ! remplace un caractere (ou une chaine de caractere) par un autre (ou par une autre chaine de caractere)
  ! dans une chaine

  implicit none
   CHARACTER (LEN=*),    INTENT(IN)    :: chaine
   CHARACTER (LEN=*),    INTENT(IN)    :: car1, car2
   CHARACTER (LEN=len(chaine))         :: new_chaine
   INTEGER                             :: i
   integer                             :: len_car1

   new_chaine = chaine
   len_car1 = len_trim(car1)

   do
    i = INDEX(new_chaine,car1)

    if (i==0) exit
    IF (i> LEN_TRIM(new_chaine)) exit
    new_chaine = new_chaine(1:i-1)//car2//new_chaine(i+len_car1:)
   END do


  return
 end function replace_car2
 


!-------------------------------------------------------------------
function verif_hkl_format(hkl_format) result(new_format)

 implicit none
  character (len=*) , intent(in)  :: hkl_format
  CHARACTER (LEN=32)              :: new_format
 
  integer                            :: i1, i2, long
  
  new_format = hkl_format
  
  long = len_trim(new_format)
  
  i1 = index(new_format, "'")
  if(i1 /=0)   new_format = remove_car(new_format, "'")
    
  i1 = index(new_format, '"')
  if(i1 /=0)   new_format = remove_car(new_format, '"')
  
  i1 = index(new_format, '(')
  if(i1 /=0)   new_format = remove_car(new_format, '(')
  
  i1 = index(new_format, ')')
  if(i1 /=0)   new_format = remove_car(new_format, ')')
  new_format = '('//trim(new_format)//')'
   
 
 return
end function verif_hkl_format
!-------------------------------------------------------------------

function multiple(n, m) result(multipl)   ! n multiple de m ?

 implicit none
  integer, intent(in)   :: n, m
  logical               :: multipl
  
  
  multipl = .false.
  if(m*int(n/m) == n) multipl = .true.

end function multiple

!-------------------------------------------------------------------
function answer_yes(input_string) result(answer_y)
 
 implicit none
  character (*), intent(in)  :: input_string
  logical                    :: answer_y
 
  answer_y = .false.
  if(input_string(1:1) == '1' .or. Input_string(1:1) == 'y'  .or. input_string(1:1) == 'Y') answer_y = .true.
 
end function answer_yes

!-------------------------------------------------------------------
function answer_no(input_string) result(answer_n)
 
 implicit none
  character (*), intent(in)  :: input_string
  logical                    :: answer_n
 
  answer_n = .false.
  if(input_string(1:1) == '0' .or. Input_string(1:1) == 'n'  .or. input_string(1:1) == 'n') answer_n = .true.
 
end function answer_no


!-------------------------------------------------------------------

 subroutine WRITE_message(string_text )
  ! ecrit un message à l'écran et dans le fichier CRYSCALC.log

  implicit none
   CHARACTER (LEN=*) ,   INTENT(IN)    :: string_text

   IF(LEN_TRIM(string_text) == 0) then
    WRITE(*,*) ' '
   else
    WRITE(*,'(a)') TRIM(string_text)
   endif
   WRITE(2,'(a)') TRIM(string_text)

  return
 end subroutine write_message


!---------------------------------------------------------------------
! verif. si un caractere est une lettre, symbole, chiffre

subroutine  check_character(input_character, alpha_char, numeric_char)
 implicit none
  CHARACTER (LEN=1), INTENT(IN)           :: input_character
  LOGICAL,           INTENT(INOUT)        :: alpha_char, numeric_char

  CHARACTER(LEN=1), DIMENSION(26)         :: alphabet_min, alphabet_maj
  CHARACTER(LEN=1), DIMENSION(10)         :: alpha
  CHARACTER(LEN=1), DIMENSION(12)         :: numeric
  INTEGER                                 :: i

  alphabet_min(1:26) = (/'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',  &
                         'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z' /)

  alphabet_maj(1:26) = (/'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',  &
                         'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' /)

  alpha(1:10)        = (/',', '!', '?', '#', ';', '_', '(', ')', '/', '\'/)

  numeric(1:12)      = (/'.', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)

  alpha_char   = .false.
  numeric_char = .false.


 ! test des lettres minuscules
  do i=1, 26
   if (input_character(1:1) == alphabet_min(i)(1:1)) then
    alpha_char = .true.
    return
   end if
  end do

 ! test des lettres majuscules
  do i=1, 26
   if (input_character(1:1) == alphabet_maj(i)(1:1)) then
    alpha_char = .true.
    return
   end if
  end do

 ! test des caracteres particuliers (;!?...)
 do i=1, 10
  if(input_character(1:1) == alpha(i)(1:1)) then
   alpha_char = .true.
   return
  endif
 end do

 ! test chiffres
  do i=1, 12
   if (input_character(1:1) == numeric(i)(1:1)) then
    numeric_char = .true.
    return
   end if
  end do



 return
end subroutine check_character
!------------------------------------------------------------------------------------

subroutine Get_current_folder_name(folder_name)

 character(len=256), intent(out) :: folder_name
 integer                         :: i, i1, i_error
 character(len=256)              :: read_line

 folder_name = '?'
 call system("dir   > temp_file")
  open (unit=22, file="temp_file")
   do i=1, 4
    read(22,'(a)', iostat=i_error) read_line
	if(i_error /=0) exit
	if(i==4 .and. index(read_line,':') /=0) then
	 folder_name = trim(read_line(index(read_line,':')-1:))
	endif
   end do

  close(unit=22)
  call system("del temp_file")

end subroutine Get_current_folder_name

!-------------------------------------------------------------------------
subroutine Get_sample_ID_from_folder_name(current_folder, sample_ID)
 character (len=*),   intent(in)    :: current_folder
 character (len=256), intent(out)   :: sample_ID
 integer                            :: i1, i2
 
   sample_ID = 'job'
   i1 = index(current_folder, '\', back=.true.) 
   if(i1 /=0) then
	i2 = index(current_folder(i1+1:),"_")
    if(i2 /=0) then	 
	  write(sample_ID, '(a)') current_folder(i1+1:i1+i2-1)
	 endif 
	endif
	
 return
end subroutine Get_sample_ID_from_folder_name 
	 
!------------------------------------------------------------------------------------
subroutine Get_WinGX_job(job, wingx_structure_directory)
 use IO_module


 implicit none
  character (len=32), intent(out) :: job
  character (len=256)             :: wingx_structure_directory
  character (len=256)             :: wingx_path_name
  character (len=256)             :: wingx_ini
  character (len=256)             :: read_line
  integer                         :: i_error, long, i1

  job                       = '?'
  wingx_structure_directory = '?'

  ! recherche de la variable d'environnement wingxdir
  call getenv('WINGXDIR', wingx_path_name)

  long = len_trim(wingx_path_name)
  if(long /=0) then
   if(wingx_path_name(long:long) == '\') then
    wingx_path_name = wingx_path_name(1:long-1)
   endif
   wingx_ini= trim(wingx_path_name) // '\wingx.ini'
  else
   call write_info('')
   call write_info('!! WinGX not installed !!')
   call write_info('')
   return
  endif

  open(unit=22, file = trim(wingx_ini), iostat=i_error)
   if(i_error /=0) then
    call write_info('')
    call write_info('Error opening the '//trim(wingx_ini)//' WinGX setting file !!')
    call write_info('')
    return
   end if

   do
    read(22, '(a)', iostat=i_error) read_line
    if(i_error /=0) exit
    read_line = adjustl(read_line)
	if(index(read_line, 'StructureDirectory=') /=0) then
     i1 = index(read_line, '=')
     read(read_line(i1+1 :), '(a)') Wingx_Structure_directory
	elseif(index(read_line, 'StructureName=') /=0) then
     i1 = index(read_line, '=')
     read(read_line(i1+1 :), '(a)') job
     exit
    endif
   end do
  close(unit=22)



 return
end subroutine Get_WinGX_job

!------------------------------------------------------------------------------------
END module macros_module


