 subroutine  create_header
  USE MACROS_module,                ONLY : u_case, Get_wingx_job,  Get_current_folder_name
  USE cryscalc_module,              ONLY : HTML_unit, text_unit, LATEX_unit, archive_CIF,  text_report, HTML_report, &
                                           LATEX_report, CRYSCALC, author, report_header, debug_proc
  USE IO_module


  implicit none
  CHARACTER (LEN=256)                 :: HTML_structural_report_file
  CHARACTER (LEN=256)                 :: TEXT_structural_report_file
  CHARACTER (LEN=256)                 :: LATEX_structural_report_file
  CHARACTER (LEN=256)                 :: LATEX_string
  INTEGER                             :: i1
  CHARACTER (LEN=10)                  :: date, time
  character (len=256)                 :: wingx_structure_dir
  character (len=10)                  :: AUTHOR_initiales

   if(debug_proc%level_2)  call write_debug_proc_level(2, "create_header")


 i1 = INDEX(archive_cif, '.', back=.true.)
 if(HTML_report) then
  WRITE(HTML_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.HTML'
  open (UNIT=HTML_unit, FILE=TRIM(HTML_structural_report_file))
 endif

 if(text_report) then
  WRITE(TEXT_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.TXT'
  open (UNIT=text_unit, FILE=TRIM(TEXT_structural_report_file))
 endif

 if(latex_report) then
  WRITE(LATEX_structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.LTX'
  open (UNIT=latex_unit, FILE=TRIM(LATEX_structural_report_file))
 endif


 if(HTML_report) then
  WRITE(HTML_unit, '(a)')  "<!-- Structural report (HTML file) created by CRYSCALC.exe"
  WRITE(HTML_unit, '(2a)')  "      Starting CIF file : ", trim(archive_CIF)
  !WRITE(HTML_unit, '(2a)') "      CRYSCALC setting file : ", trim(cryscalc%ini)
  WRITE(HTML_unit, '(a)')  ""
  WRITE(HTML_unit, '(a)')  "     Web site : http://www.cdifx.univ-rennes/cryscalc.htm"
  WRITE(HTML_unit, '(5a)') "     Version : ", trim(cryscalc%version), " [", trim(cryscalc%author), "]"
  WRITE(HTML_unit, '(a)') "!-->"

  WRITE(HTML_unit, '(a)') "<HTML>"
  WRITE(HTML_unit, '(a)') "<HEAD>"
  call write_HTML_css(HTML_unit, 'report')
  WRITE(HTML_unit, '(a)') "<TITLE>Structural report</TITLE>"
  WRITE(HTML_unit, '(a)') "</HEAD>"
  WRITE(HTML_unit, '(a)') "<BODY BGCOLOR='#FFFFFF' style='font-family:Times new roman; font-size:14; line-height:150%'>"
 endif

 if(LATEX_report) then
  call Latex_preambule
 endif


  call date_and_time(date, time)
  !call Get_Wingx_job(job, wingx_structure_dir)

 if(report_header) then

   if(AUTHOR%name /= '?' .and. AUTHOR%first_name /= '?') then
    write(AUTHOR_initiales, '(4a)') u_case(AUTHOR%first_name(1:1)), '.', u_case(AUTHOR%name(1:1)), '.'
    if(HTML_report) then
     WRITE(HTML_unit, '(a)')  "<hr>"
     WRITE(HTML_unit, '(18a)') '<p class="header">HTML report created by ', trim(AUTHOR_initiales),       &
                                            " [",trim(AUTHOR%team), "]  -  ",                         &
                                            " Date: " , date(7:8),'-',date(5:6),'-',date(1:4), '  at ', time(1:2),':',  &
                                                        time(3:4),':',time(5:6), "<br>"
    endif
    if(LATEX_report) then
     write(LATEX_unit, '(18a)') '\header{\LaTeX report created by ', trim(AUTHOR_initiales),     &
                                " [",trim(AUTHOR%team), "]  -  ",                              &
                                " Date: " , date(7:8),'-',date(5:6),'-',date(1:4), '  at ', time(1:2),':',  &
                                            time(3:4),':',time(5:6), "} \\"
    endif

  ELSE
    if(HTML_report) then
     WRITE(HTML_unit, '(a)')  "<hr>"
     WRITE(HTML_unit, '(13a)')' <p class="header"> Date: ' , date(7:8),'-',date(5:6),'-',date(1:4), '  at ', time(1:2),':', &
                              time(3:4),':',time(5:6), "<br>"
    endif
    if(LATEX_report) then
     WRITE(LATEX_unit, '(13a)')'\header{Date: ' , date(7:8),'-',date(5:6),'-',date(1:4), '  at ', time(1:2),':', &
                                                  time(3:4),':',time(5:6), "} \\"
    endif
   ENDIF

   !if(wingx_structure_dir =='?') call Get_current_folder_name(wingx_structure_dir)
   call Get_current_folder_name(wingx_structure_dir)

   !i1 = index(wingx_structure_dir, '_')
   !if(i1 /=0 .and. i1 > 1) then
   ! write(job, '(a)') wingx_structure_dir(1:i1-1)
   !endif
  if(HTML_report) then
   WRITE(HTML_unit, '(3a)') 'Working directory: <font face="courier">', trim(wingx_structure_dir), '</font><br>'
   WRITE(HTML_unit, '(3a)') "Input CIF file : ", trim(archive_CIF), '<br>'
   WRITE(HTML_unit, '(5a)') "CRYSCALC version : ", trim(cryscalc%version), " [", trim(cryscalc_author), "]</p>"
   WRITE(HTML_unit, '(a)')  "<hr>"
  endif
  if(LATEX_report) then
   LATEX_string = wingx_structure_dir
   call check_LATEX_file_name(LATEX_string)
   WRITE(LATEX_unit, '(3a)') '\header{Working directory: ', trim(LATEX_string), '} \\'
   LATEX_string = archive_cif
   call check_LATEX_file_name(LATEX_string)
   WRITE(LATEX_unit, '(3a)') "\header{Input CIF file : ", trim(LATEX_string), '} \\'
   WRITE(LATEX_unit, '(5a)') "\header{CRYSCALC version : ", trim(cryscalc%version), " [", trim(cryscalc_author), "]} \\"
   write(LATEX_unit, '(a)')  "\rule{\linewidth}{0.5pt}"
   write(LATEX_unit, '(a)')  "\vspace{1.0cm}"
   WRITE(LATEX_unit, '(a)')  ""
  endif
 endif

 return
end subroutine create_header
