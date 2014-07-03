!     Last change:  TR    2 May 2007    5:19 pm
!------------------------------------------------------------------------------
module external_applications_module

 PRIVATE
 PUBLIC    :: launch_browser,    &
              launch_word,       &
			  launch_pdflatex

 contains

subroutine launch_browser(HTML_file)
 USE cryscalc_module,   ONLY : my_browser, debug_proc
 USE macros_module,     ONLY : Get_current_folder_name
 implicit none
  CHARACTER (LEN=*),   INTENT(IN)      :: HTML_file
  CHARACTER (LEN=256)                  :: current_directory
  CHARACTER (LEN=256)                  :: DOS_command

  if(debug_proc%level_2)  call write_debug_proc_level(2, "launch_browser")
  
   WRITE(DOS_command, '(3a)') TRIM(my_browser%name),' ', TRIM(HTML_file)
   call system(trim(DOS_command))
   
  return
end subroutine launch_browser

subroutine launch_WORD(word_file)
 !USE RealWin,           ONLY : get_current_directory, spawn
 USE cryscalc_module,   ONLY : my_word, debug_proc
 USE macros_module,     ONLY : Get_current_folder_name
 implicit none
  CHARACTER (LEN=*),   INTENT(IN)      :: word_file
  CHARACTER (LEN=256)                  :: current_directory
  CHARACTER (LEN=256)                  :: DOS_command

  if(debug_proc%level_2)  call write_debug_proc_level(2, "launch_WORD")
  
  IF(LEN_TRIM(my_word%name) == 0) return
  ! repertoire de travail ?
  ! utilisation de la fonction 'get_current_directory' de RealWin
  !call get_current_directory(directory=current_directory)

   ! jan. 11 : fonction "get_current_directory" de RealWin remplacee par routine Get_current_folder_name
  call Get_current_folder_name(current_directory)

  IF(LEN_TRIM(word_file) /=0) then
   WRITE(DOS_command, '(4a)') TRIM(my_word%name),' ', TRIM(current_directory), TRIM(word_file)
  else
   WRITE(DOS_command, '(a)') TRIM(my_word%name)
  endif
   !call spawn(file_name=TRIM(DOS_command))
   call system(trim(DOS_command))
 
end subroutine launch_word

subroutine launch_pdflatex(pdflatex_file)
 USE cryscalc_module,   ONLY : my_pdflatex, debug_proc
 USE macros_module,     ONLY : Get_current_folder_name
 implicit none
  CHARACTER (LEN=*),   INTENT(IN)      :: pdflatex_file
  CHARACTER (LEN=256)                  :: current_directory
  CHARACTER (LEN=256)                  :: DOS_command

  if(debug_proc%level_2)  call write_debug_proc_level(2, "laucnh_pdflatex")
 
   
   WRITE(DOS_command, '(3a)') TRIM(my_pdflatex%name),' ', TRIM(pdflatex_file)
   call system(trim(DOS_command))
   
  return
end subroutine launch_pdflatex


END module external_applications_module


