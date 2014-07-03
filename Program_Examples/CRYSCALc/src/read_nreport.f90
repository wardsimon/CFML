!     Last change:  TR    2 May 2007    3:56 pm
! lecture du fichier nreport.HTML

subroutine read_nreport_html()
 USE cryscalc_module,              ONLY : message_text, my_browser, browse_nreport, debug_proc 
 USE macros_module,                ONLY : test_file_exist
 USE IO_module,                    ONLY : write_info
 use external_applications_module, ONLY : launch_browser
! USE RealWin,        ONLY : get_current_directory, spawn

 implicit none
   LOGICAL                                    :: file_exist
   CHARACTER (LEN=265)                        :: read_line
   INTEGER                                    :: i_error
   INTEGER                                    :: i1, i2
   CHARACTER (LEN=64)                         :: HTML_string
   INTEGER                                    :: HTML_string_long
   CHARACTER (LEN=256)                        :: current_directory
   CHARACTER (LEN=256)                        :: DOS_command


   if(debug_proc%level_2)  call write_debug_proc_level(2, "READ_NREPORT_HTML")
   
 call test_file_exist('nreport.html', file_exist, 'out')
 IF(.NOT. file_exist) return


 if (browse_nreport .and. my_browser%exist) then
  call launch_browser('nreport.html')
 else
  OPEN(UNIT=51, FILE='nreport.html')
  do
   READ(UNIT=51, fmt='(a)', IOSTAT=i_error) read_line
   IF(i_error < 0) EXIT   ! fin du fichier
   call HTML_to_text(read_line)
   call write_info('  '//TRIM(read_line))   
  END do
  close (UNIT=51)
 end if

 return
end subroutine read_nreport_html

!---------------------------------------------------------------------------------

subroutine HTML_to_text(string)
 use cryscalc_module, only : debug_proc
 USE macros_module,   ONLY: replace_car, replace_car2, remove_car, clear_string
 implicit none
 CHARACTER (LEN=*), INTENT(INOUT) :: string
 
 if(debug_proc%level_2)  call write_debug_proc_level(2, "HTML_TO_TEXT")

! call replace_car(string, "&deg;"             , " deg.")
! call replace_car(string, "&lt;"              , " < ")
! call replace_car(string, "<TD>"              , "  ")
! call replace_car(string, "</TH>"             , "  ")
! call replace_car(string, "&nbsp;"            , "   ")
! call replace_car(string, "<TD ALIGN=RIGHT>"  , "   ")


! call remove_car(string, "<P>" )
! call remove_car(string, "<BR>")
! call remove_car(string, "<HR>")
! call remove_car(string, "</TD>")
! call remove_car(string, "<TR>")
! call remove_car(string, "</TR>")
! call remove_car(string, "<TH ALIGN='LEFT'>" )
! call remove_car(string, "<TH>" )
! call remove_car(string, "<TT>" )
! call remove_car(string, "</TT>" )

! call remove_car(string, "<HTML>")
! call remove_car(string, "</HTML>")
! call remove_car(string, "<HEAD>")
! call remove_car(string, "</HEAD>")
! call remove_car(string, "<TITLE>")
! call remove_car(string, "</TITLE>")
! call remove_car(string, "<BODY>")
! call remove_car(string, "</BODY>")
! call remove_car(string, "<TABLE BORDER>")
! call remove_car(string, "</TABLE>")


! call remove_car(string, "<H1>")
! call remove_car(string, "</H1>")
! call remove_car(string, "<H2>")
! call remove_car(string, "</H2>")
! call remove_car(string, "<H3>")
! call remove_car(string, "</H3>")

! call remove_line(string, "<ADDRESS>")


 string = replace_car2(string, "&deg;"             , " deg.")
 string = replace_car(string, "&lt;"              , " < ")
 string = replace_car(string, "<TD>"              , "  ")
 string = replace_car(string, "</TH>"             , "  ")
 string = replace_car(string, "&nbsp;"            , "  ")
 string = replace_car(string, "<TD ALIGN=RIGHT>"  , "  ")

 string = remove_car(string, "<P>" )
 string = replace_car(string, "<BR>", " ")
 string = remove_car(string, "<HR>")
 string = remove_car(string, "</TD>")
 string = remove_car(string, "<TR>")
 string = remove_car(string, "</TR>")
 string = remove_car(string, "<TH ALIGN='LEFT'>" )
 string = remove_car(string, "<TH>" )
 string = remove_car(string, "<TT>" )
 string = remove_car(string, "</TT>" )
 
 string = remove_car(string, "<HTML>")
 string = remove_car(string, "</HTML>")
 string = remove_car(string, "<HEAD>")
 string = remove_car(string, "</HEAD>")
 string = remove_car(string, "<TITLE>")
 string = remove_car(string, "</TITLE>")
 string = remove_car(string, "<BODY>")
 string = remove_car(string, "</BODY>")
 string = remove_car(string, "<TABLE BORDER>")
 string = remove_car(string, "</TABLE>")

 string = remove_car(string, "<H1>")
 string = remove_car(string, "</H1>")
 string = replace_car(string, "<H2>" , "*** ")
 string = replace_car(string, "</H2>", " ***" )
 string = replace_car2(string, "<H3>" , " >> ")
 string = remove_car(string, "</H3>")


 string = clear_string(string, "<ADDRESS>")



 return
end subroutine HTML_to_text
