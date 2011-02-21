!     Last change:  TR   29 Nov 2007    6:18 pm
subroutine create_CRYSCAL_HTML()
! creation du manuel au format HTML
USE cryscal_module,               ONLY : nb_help_max, help_string, HTML_unit, browse_cryscal_HTML, my_browser, &
                                         cryscal_version, cryscal_author
USE text_module
USE external_applications_module, ONLY : launch_browser
USE IO_module




implicit none
 INTEGER             :: i, k
 CHARACTER (LEN=512) :: HTML_text
 CHARACTER (LEN=256) :: current_directory
 CHARACTER (LEN=256) :: DOS_command


 OPEN(UNIT=HTML_unit, FILE='cryscal.html')
  WRITE(HTML_unit, '(a)')  "<!-- CRYSCAL user's guide (HTML file) created by CRYSCAL.exe"
  WRITE(HTML_unit, '(a)')  "     Web site : http://www.cdifx.univ-rennes/cryscal.htm"
  WRITE(HTML_unit, '(2a)') "     Version : ", trim(cryscal_version)
  WRITE(HTML_unit, '(2a)') "     Author  : ", trim(cryscal_author)
  WRITE(HTML_unit, '(a)')  "!-->"

  WRITE(HTML_unit, '(a)') "<HTML>"
  WRITE(HTML_unit, '(a)') "<HEAD>"

  call write_HTML_css(HTML_unit)

  WRITE(HTML_unit, '(a)') '<A NAME="cryscal_main"></A>'
  WRITE(HTML_unit, '(a)') "<TITLE>CRYSCAL user's guide</TITLE>"
  WRITE(HTML_unit, '(a)') "</HEAD>"
  WRITE(HTML_unit, '(a)') '<BODY BGCOLOR="#FFFFFF">'
  WRITE(HTML_unit, '(a)') "<br><br>"
  WRITE(HTML_unit, '(a)') "<p class='title_main'>CRYSCAL</p>"
  !WRITE(HTML_unit, '(a)') "<center><h1><FONT COLOR="#FF0000">CRYSCAL</FONT></h1></center>"
  WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"


  ! TITLE
  WRITE(HTML_unit, '(a)') "<center>"
  do i=2, title_lines_nb
   WRITE(HTML_unit, '(a)') TRIM(title_line(i))
  end do
  WRITE(HTML_unit, '(a)') "</center>"
  WRITE(HTML_unit, '(a)') ""


  !WRITE(HTML_unit, '(a)') '<p class="title_4">'
  WRITE(HTML_unit, '(a)') '<ul>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#intro">&nbsp;&nbsp;Introduction</A>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#kw_list">&nbsp;&nbsp;List of CRYSCAL keywords</A>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#kw_details">&nbsp;&nbsp;Details of CRYSCAL keywords</A>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cla_section">&nbsp;&nbsp;List of CRYSCAL command line arguments</A>'
  WRITE(HTML_unit, '(a)') "<li class='item_2'><A HREF='#cryscal_news'>&nbsp;&nbsp;What's new ?</A>"
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cryscal_ini">&nbsp;&nbsp;CRYSCAL.ini setting file</A>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cryscal_cfl">&nbsp;&nbsp;Examples of .CFL input files</A>'

  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cryscal_links">&nbsp;&nbsp;CRYSCAL download and links</a>'
  WRITE(HTML_unit, '(a)') '</ul>'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '</p>'


  WRITE(HTML_unit, '(a)') '<A NAME="intro"></A><p class="title_1"> Introduction :<br></p>'

  ! HEADER
  call def_header_lines()
  do i=1, header_lines_nb
   WRITE(HTML_unit, '(a)') TRIM(header_line(i))
  end do
  WRITE(HTML_unit, '(a)') ""
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscal_main">[return]</a>'
  WRITE(HTML_unit, '(a)') "</pre>"

  WRITE(HTML_unit, '(a)') '<A NAME="kw_list"></A><p class="title_1"> List of CRYSCAL keywords :<br></p>'
  WRITE(HTML_unit, '(a)') '<ul><ul>'
  do i=1, nb_help_max
   WRITE(HTML_text, '(5a)') '   <li class="item"><A HREF="#', TRIM(help_string(i)), '">', TRIM(help_string(i)), '</A>'
   WRITE(HTML_unit, '(a)') TRIM(HTML_text)
  end do
  WRITE(HTML_unit, '(a)') '</ul></ul><br>'
  WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"


  call def_keywords_lines()
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscal_main">[return]</a>'

  WRITE(HTML_unit, '(a)') '<A NAME="kw_details"></A><p class="title_1"> Details of CRYSCAL keywords :<br></p>'
  do i=1, nb_help_max
   WRITE(HTML_unit, '(3a)') '<A NAME="', TRIM(help_string(i)), '"></A>'
   WRITE(HTML_unit, '(3a)') "<p class='title_3'>", trim(HELP_line(i,2)), "</p> <pre style='font-size:14'>"
   do k=4, HELP_lines_nb(i)
    WRITE(HTML_unit, '(a)') TRIM(HELP_line(i,k))
   end do
   WRITE(HTML_unit, '(a)')'<A HREF="#kw_list">[return to keywords list]</a>'
   WRITE(HTML_unit, '(a)')'<hr></pre>'
  end do



  call Def_command_line_arguments
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscal_main">[return]</a>'
  WRITE(HTML_unit, '(a)') '<A NAME="cla_section"></A><p class="title_1"> List of CRYSCAL command line arguments :<br></p>'
  !WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"
  !do i=1, CLA_lines_nb
  ! WRITE(HTML_unit, '(a)') TRIM(cla_line(i))
  !end do
  !WRITE(HTML_unit, '(a)') ""
  !WRITE(HTML_unit, '(a)') "</pre>"

  !WRITE(HTML_unit, '(a)') '<ul>'
  do i=1, CLA_nb
   WRITE(HTML_unit, '(3a)') "<p class='title_3'>", trim(CLA_line(i,1)), "</p><pre style='font-size:14'>"
   do k=2, CLA_lines_nb(i)
    WRITE(HTML_unit, '(a)') TRIM(CLA_line(i,k))
   end do
   WRITE(HTML_unit, '(a)') '</pre>'
   !WRITE(HTML_text, '(2a)') '   <li class="item">', TRIM(CLA_line(i))
   !WRITE(HTML_unit, '(a)') TRIM(HTML_text)
  end do
  !WRITE(HTML_unit, '(a)') '</ul><br>'
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscal_main">[return]</a>'
  WRITE(HTML_unit, '(a)') "</pre>"


! --------------- what's new in CRYSCAL ? ------------------------

  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') "<A NAME='cryscal_news'></A><p class='title_1'> What's new in CRYSCAL ?<br></p>"
  WRITE(HTML_unit, '(a)') '<pre>'
  call WRITE_CRYSCAL_NEWS("html")
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscal_main">[return]</a>'
  WRITE(HTML_unit, '(a)') '</pre>'


! --------------- setting file -----------------------------------
  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') '<A NAME="cryscal_ini"></A><p class="title_1"> CRYSCAL.ini setting file :<br></p>'
  WRITE(HTML_unit, '(a)') '<p class="retrait_1">A setting file can be used by <b>CRYSCAL</b>, containing '
  WRITE(HTML_unit, '(a)') 'the definition of different parameters such as external applications that can be'
  WRITE(HTML_unit, '(a)') 'launched from <b>CRYSCAL</b> (browser, editor ...) or defaut values about'
  WRITE(HTML_unit, '(a)') 'diffractometer, author, structure solution and refiement programs  ... '
  WRITE(HTML_unit, '(a)') 'This setting file, called <font face="courier">cryscal.ini</font> has to be located'
  WRITE(HTML_unit, '(a)') 'in the folder related to <b>CRYSCAL</b> through the <font face="courier">CRYSCAL</font>'
  WRITE(HTML_unit, '(a)') 'environment variable.</p>'
  WRITE(HTML_unit, '(a)') '<p class="retrait_1">Example of setting file:</p>'
  WRITE(HTML_unit, '(a)') '<pre>'
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[EXTERNAL APPLICATIONS]</font>'
  WRITE(HTML_unit, '(a)') '    browser = "C:\Program Files\Mozilla Firefox\firefox.exe"'
  WRITE(HTML_unit, '(a)') '    editor  = "C:\Program Files\Keditw\KEDITW32.EXE"'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[WEB ADDRESS]</font>'
  WRITE(HTML_unit, '(a)') '    fps        = www.ill.fr/dif/Soft/fp/'
  WRITE(HTML_unit, '(a)') '    cri        = www.cri.univ-rennes1.fr/'
  WRITE(HTML_unit, '(a)') '    cdifx      = www.cdifx.univ-rennes1.fr/'
  WRITE(HTML_unit, '(a)') '    cryscal    = www.cdifx.univ-rennes1.fr/cryscal'
  WRITE(HTML_unit, '(a)') '    reciprocs  = www.cdifx.univ-rennes1.fr/reciprocs'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[DEVICE]</font>'
  WRITE(HTML_unit, '(a)') '    diffractometer = APEXII AXS Bruker'
  WRITE(HTML_unit, '(a)') '    !diffractometer = KCCD Nonius'
  WRITE(HTML_unit, '(a)') '    laboratory     = CDIFX Rennes'
  WRITE(HTML_unit, '(a)') '    radiation      = X_Mo'
  WRITE(HTML_unit, '(a)') '    wave_A         = 0.71073'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[AUTHOR]</font>'
  WRITE(HTML_unit, '(a)') '    name        = ROISNEL'
  WRITE(HTML_unit, '(a)') '    first_name  = Thierry'
  WRITE(HTML_unit, '(2a)')'    address     = Centre de Diffractométrie X, UMR6511 CNRS Université de Rennes 1,', &
                                            'Sciences Chimiques de Rennes, 35042 RENNES Cedex France'
  WRITE(HTML_unit, '(a)') '    email       = thierry.roisnel@univ-rennes1.fr'
  WRITE(HTML_unit, '(a)') '    web         = www.cdifx.univ-rennes1.fr'

  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[CREATE INS]</font>'
  WRITE(HTML_unit, '(a)') '    temperature = 100K       ! experimental temperature value'
  WRITE(HTML_unit, '(a)') '    u_threshold = 0.1        ! atoms with U_iso > U_threshold will be excluded'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[PARAMETERS]</font>'
  WRITE(HTML_unit, '(a)') '    i_sig       = 3.         ! used in the SEARCH_GROUP procedure'
  WRITE(HTML_unit, '(a)') '    threshold   = 0.03       ! used in the SEARCH_GROUP procedure '
  WRITE(HTML_unit, '(a)') '    d_max_A     = 3.5        ! used with the CONN keyword (connectivity calculation)'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[COMMAND LINE ARGUMENTS]</font>'
  WRITE(HTML_unit, '(a)') '    create_ACE  = 1          ! .ACE file for Carine'
  WRITE(HTML_unit, '(a)') '    create_CEL  = 1          ! .CEL file for PowderCELL'
  WRITE(HTML_unit, '(a)') '    create_CFL  = 1          ! .CFL file for CRYSCAL'
  WRITE(HTML_unit, '(a)') '    create_INS  = 1          ! .INS file for SHELXL'
  WRITE(HTML_unit, '(a)') '    create_FST  = 1          ! .FST file for FP Studio'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[PROGRAMS]</font>'
  WRITE(HTML_unit, '(a)') '    structure_solution_name         = SIR97'
  WRITE(HTML_unit, '(3a)')'    structure_solution_reference    = A. Altomare, M. C. Burla, M. Camalli, G. Cascarano, ', &
                                                                'C. Giacovazzo, A. Guagliardi, A. G. G. Moliterni, ',   &
																'G. Polidori, R. Spagna, J. Appl. Cryst. (1999) 32, 115-119'
  WRITE(HTML_unit, '(a)') '    structure_solution_cif_ref      = SIR97 (Altomare et al., 1999)'
  WRITE(HTML_unit, '(a)') '    structure_refinement_name       = SHELXL-97'
  WRITE(HTML_unit, '(a)') '    structure_refinement_reference  = Sheldrick G.M., Acta Cryst. A64 (2008), 112-122'
  WRITE(HTML_unit, '(a)') '    structure_refinement_cif_ref    = SHELXL-97 (Sheldrick, 2008)'

  WRITE(HTML_unit, '(a)') '    absorption_correction_name      = SADABS'
  WRITE(HTML_unit, '(2a)')'    absorption_correction_reference = Sheldrick G.M. (2002), SADABS Bruker AXS Inc., Madison, ', &
                                                                'Wisconsin, USA'
  WRITE(HTML_unit, '(2a)')'    absorption_correction_cif_ref   = Sheldrick G.M. (2002), SADABS Bruker AXS Inc., Madison, ', &
                                                                'Wisconsin, USA'

  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '<A HREF="#cryscal_main">[return]</a>'
  WRITE(HTML_unit, '(a)') '</pre>'

  
! ----------------------------------------------------------------
  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') '<A NAME="cryscal_cfl"></A><p class="title_1"> Examples of .CFL input files :<br></p>'
  WRITE(HTML_unit, '(a)') '<ul>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscal/cryscal_pnma.cfl">', &
                           '  Pnma space group information</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscal/cryscal_pbnm_pnma.cfl">', &
                           '  Transformation from Pbnm to Pnma space group</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscal/cryscal_si_x.cfl">', &
                           '  Simulation of a X-ray diffraction pattern for Si</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscal/cryscal_si_n.cfl">', &
                           '  Simulation of a neutron diffraction pattern for Si</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscal/cryscal_y2o3.cfl">', &
                           '  Atom connectivity in Y2O3</a>'
  WRITE(HTML_unit, '(a)') '<p>'
  WRITE(HTML_unit, '(a)') '</ul>'
 
  
! ----------------------------------------------------------------
  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') '<A NAME="cryscal_links"></A><p class="title_1"> CRYSCAL links :<br></p>'
  WRITE(HTML_unit, '(a)') '<ul>'
  WRITE(HTML_unit, '(a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscal/cryscal.exe">CRYSCAL.exe</a>'
  WRITE(HTML_unit, '(a)') '<p>'

  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscal/cryscal.ini">', &
                          'Example of CRYSCAL setting file</a>'
  WRITE(HTML_unit, '(a)') '<p>'

  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.iucr.org/resources/commissions/crystallographic-computing/', &
                          'newsletters/1/crysfml" target="_blank">CrysFML:</A>'
  WRITE(HTML_unit, '(2a)')' <font color="black">Crystallographic Fortran Modules Library by J. Rodriguez-Carvajal and J. ',   &
                          'González-Platas</font>'
  WRITE(HTML_unit, '(a)') '</ul>'
  WRITE(HTML_unit, '(a)') ''

  !WRITE(HTML_unit, '(a)') '<p class="retrait_1"><b>Author:</b><br>'
  !WRITE(HTML_unit, '(a)') 'Thierry Roisnel <A HREF="mailto:cdifx@univ-rennes1.fr">cdifx.univ-rennes1.fr</A>'
  WRITE(HTMl_unit, '(a)')  "<br><br><hr>"
  WRITE(HTML_unit, '(3a)') "<p class='retrait_2'>last updated: <i>CDIFX - ", trim(cryscal_version), "</i></P>"
  !WRITE(HTML_unit, '(a)') 'Web site: <a href="www.cdifx.univ-rennes1.fr/CRYSCAL" target="_blank"> www.cdifx.univ-rennes1.fr/CRYSCAL</a>'

  WRITE(HTML_unit, '(a)') "</BODY>"
  WRITE(HTML_unit, '(a)') "</HTML>"
CLOSE(UNIT=HTML_unit)


  call write_info("")
  call write_info("   The CRYSCAL user's guide (cryscal.html) in HTML format has been created." )
  call write_info("")


  IF(browse_cryscal_HTML .and. my_browser%exist) then
   call launch_browser('cryscal.html')
  ENDIF


end subroutine create_CRYSCAL_HTML

!--------------------------------------------------------------------------------

subroutine write_HTML_css(HTML_unit)

 implicit none
  integer, intent(in)   :: HTML_unit
  WRITE(HTML_unit, '(a)') "<style TYPE='text/css'>"
  WRITE(HTML_unit, '(a)') "  h2          { font-family:'Trebuchet MS',Arial; color:blue;margin-left:50px; }"
  WRITE(HTML_unit, '(2a)') "  .title_main { font-family:'Trebuchet MS',Arial; font-size:24px; color:#AA0000;", &
                           " text-align:center; font-weight:bold; }"
  WRITE(HTML_unit, '(2a)') "  .title_1    { font-family:'Trebuchet MS',Arial; font-size:20px; color:#880000;",  &
                           " margin-left:20px; background-color:#f5f5f5; }"
  WRITE(HTML_unit, '(2a)') "  .title_2    { font-family:'Trebuchet MS',Arial; font-size:20px; color:#880000;",  &
                           " margin-left:50px; background-color:#f5f5f5; }"
  WRITE(HTML_unit, '(2a)') "  .title_3    { font-family:'Courier new';        font-size:14px; color:#660000;", &
                           " font-weight:bold; margin-left:30px;}"
  WRITE(HTML_unit, '(2a)') "  .title_4    { font-family:'Trebuchet MS',Arial; font-size:18px; color:#880000;",  &
                           " margin-left:20px; background-color:#f5f5f5; }"
  WRITE(HTML_unit, '(2a)') "  .retrait_1  { font-family:'Trebuchet MS',Arial; font-size:14px; margin-left:50px;", &
                           " margin-right:50px; text-align:justify; line-height:18.0pt;  }"
  WRITE(HTML_unit, '(2a)') "  .retrait_2  { font-family:'Trebuchet MS',Arial; font-size:14px; margin-left:50px;",  &
                           " margin-right:50px; text-align:justify; line-height:18.0pt; color:#AAAAAA; }"
  WRITE(HTML_unit, '(2a)') "  .retrait_3  { font-family:'Trebuchet MS',Arial; font-size:14px; margin-left:100px;", &
                           " margin-right:50px; text-align:justify; line-height:18.0pt;  }"
  WRITE(HTML_unit, '(a)') "  .item       { font-family:'Trebuchet MS',Arial; font-size:14px; color:#aa0000;}"
  WRITE(HTML_unit, '(a)') "  .item_2     { font-family:'Trebuchet MS',Arial; font-size:16px; color:#aa0000;}"
  WRITE(HTML_unit, '(a)') ""
  WRITE(HTML_unit, '(a)') " a:link,"
  WRITE(HTML_unit, '(a)') "  a:active"
  WRITE(HTML_unit, '(a)') "    { color:#666666;"
  WRITE(HTML_unit, '(a)') "      text-decoration:none;"
  WRITE(HTML_unit, '(a)') "      background-color:transparent;}"
  WRITE(HTML_unit, '(a)') "  a:visited"
  WRITE(HTML_unit, '(a)') "    { color:#666666;"
  WRITE(HTML_unit, '(a)') "      text-decoration:none;"
  WRITE(HTML_unit, '(a)') "      background-color:transparent;}"
  WRITE(HTML_unit, '(a)') "  a:hover"
  WRITE(HTML_unit, '(a)') "    { color:rgb(92,92,92);"
  WRITE(HTML_unit, '(a)') "      text-decoration:none;"
  WRITE(HTMl_unit, '(a)') "      font-weight:bold;"
  WRITE(HTML_unit, '(a)') "      background-color:transparent;}"

  WRITE(HTML_unit, '(a)') "</style>"

  return
end subroutine write_HTML_css
!--------------------------------------------------------------------------------


subroutine create_structural_report()
 USE MACROS_module,                ONLY : test_file_exist, remove_car, replace_car, nombre_de_colonnes, l_case, u_case, &
                                          check_character
 USE cryscal_module,               ONLY : CIF_unit, HTML_unit, CIF_parameter, archive_CIF, long_report, SQUEEZE, &
                                          my_browser, cryscal_ini, message_text, DEBUG_file, cryscal_version, cryscal_author, &
                                          structure_solution, structure_refinement
 USE external_applications_module, ONLY : launch_browser, launch_word
 USE IO_module


 implicit none
 LOGICAL                             :: file_exist
 CHARACTER (LEN=256)                 :: CIF_input_line, CIF_field, CIF_field_value
 CHARACTER (LEN=256)                 :: structural_report_file
 CHARACTER (LEN=256)                 :: HTML_string, new_HTML_string
 CHARACTER (LEN=256), dimension(4)   :: HTML_str
 CHARACTER (LEN=32)                  :: GIF_file
 INTEGER                             :: i_error, long_field, long, i, i1
 integer                             :: n, n_op, op_numor
 character (len=256), dimension(192) :: op_string
 integer,   dimension(500)           :: op_n       ! numero de l'op. de sym.
 integer,   dimension(500,3)         :: op_t       ! partie translation
 integer,   dimension(3)             :: t
 integer                             :: n_sym_htab, n_sym_dist, n_sym_ang, n_sym_tor
 integer                             :: nb_col
 CHARACTER (LEN=64)                  :: device_mark, device_type

 CHARACTER (LEN=12), dimension(1000) :: site_sym
 CHARACTER (len=12), dimension(1000) :: dico
 INTEGER           , dimension(1000) :: num_site_sym

 CHARACTER (LEN=12), dimension(2)    :: dist_atom
 CHARACTER (LEN=12)                  :: dist_value
 CHARACTER (LEN=12)                  :: dist_sym

 CHARACTER (LEN=12), dimension(3)    :: ang_atom
 CHARACTER (LEN=12)                  :: ang_value
 CHARACTER (LEN=12), dimension(2)    :: ang_sym

 CHARACTER (LEN=12), dimension(4)    :: torsion_ang_atom
 CHARACTER (LEN=12)                  :: torsion_ang_value
 CHARACTER (LEN=12), dimension(4)    :: torsion_sym


 CHARACTER (LEN=12)                  :: atom_label, atom_symbol, atom_x, atom_y, atom_z, atom_Ueq
 CHARACTER (LEN=12)                  :: atom_U11, atom_U22, atom_U33, atom_U23, atom_U13, atom_U12
 CHARACTER (len=12)                  :: site_label_D, site_label_H, site_label_A
 CHARACTER (len=12)                  :: dist_DH, dist_HA, dist_DA, angle_DHA, site_sym_A

 LOGICAL                             :: alpha_car, num_car



 file_exist = .false.
 call test_file_exist(TRIM(archive_CIF), file_exist)
 IF(.NOT. file_exist) then
  !call write_info('')
  !call write_info('   archive.cif file is missing !!')
  !call write_info('')
  if(DEBUG_file%write) then
   call write_DEBUG_file(' '//trim(archive_cif)//' file is missing !', '')
   call write_DEBUG_file(' ', '')
   call write_DEBUG_file(' CRYSCAL will be stopped.', '')
   call write_DEBUG_file('', '')
  end if
  stop
 endif



 open (unit = CIF_unit, file = TRIM(archive_CIF))
  do
   READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
   IF(i_error < 0) exit
   IF(i_error /=0) then
    call write_info('')
    call write_info(' !! Error reading archive.CIF file !!')
    call write_info('')
    return
   endif
   IF(LEN_TRIM(CIF_input_line) == 0) cycle
   read(CIF_input_line, *, IOSTAT=i_error) CIF_field
   IF(i_error /=0) then
    call write_info('')
    call write_info(' !! Error reading archive.CIF file !!')
    call write_info('')
    return
   endif


   CIF_field = l_case(CIF_field)
   long_field = LEN_TRIM(CIF_field)
   CIF_field_value  = CIF_input_line(long_field+1:)
   CIF_field_value  = ADJUSTL(CIF_field_value)
   long = LEN_TRIM(CIF_field_value)


   IF(CIF_field_value(1:1)== "'" .AND. CIF_field_value(long:long) == "'") CIF_field_value = CIF_field_value(2:long-1)


   select case (CIF_field(1:long_field))
       case ('_chemical_formula_moiety')
        READ(CIF_field_value, '(a)') CIF_parameter%formula_moiety
        !call remove_car(CIF_parameter%formula_moiety, ' ')
        !CIF_parameter%formula_moiety = remove_car(trim(CIF_parameter%formula_moiety), ' ')

       case ('_chemical_formula_sum')
        READ(CIF_field_value, '(a)') CIF_parameter%formula_sum


       case ('_chemical_formula_weight')
        READ(CIF_field_value, '(a)') CIF_parameter%formula_weight


       case ('_symmetry_cell_setting', '_space_group_crystal_system')
        READ(CIF_field_value, '(a)') CIF_parameter%symmetry_cell_setting


       case ('_symmetry_space_group_name_h-m', '_space_group_name_h-m')
        READ(CIF_field_value, '(a)') CIF_parameter%symmetry_space_group
        !READ(read_line(2:long-1), '(a)') CIF_parameter%symmetry_space_group

       case ('_symmetry_int_tables_number')
        READ(CIF_field_value, '(a)') CIF_parameter%symmetry_IT_number

       case ('_symmetry_equiv_pos_as_xyz')
        n_op = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         if(i_error /=0) exit
         if(len_trim(CIF_input_line)==0) exit
         if(CIF_input_line(1:1) == '#' .or. CIF_input_line(1:1) == '!' ) exit
         if(CIF_input_line(1:1) == '_') then
          backspace(unit = CIF_unit)
          exit
         endif
         n_op = n_op + 1
         read(CIF_input_line, '(a)') op_string(n_op)
         op_string(n_op) = trim(op_string(n_op))
         op_string(n_op) = op_string(n_op)(2:len_trim(op_string(n_op))-1) ! enleve les caracteres "'"
        end do

       case ('_cell_length_a')
        READ(CIF_field_value, '(a)') CIF_parameter%cell_length_a


       case ('_cell_length_b')
        READ(CIF_field_value, '(a)') CIF_parameter%cell_length_b


       case ('_cell_length_c')
        READ(CIF_field_value, '(a)') CIF_parameter%cell_length_c


       case ('_cell_angle_alpha')
        READ(CIF_field_value, '(a)') CIF_parameter%cell_angle_alpha


       case ('_cell_angle_beta')
        READ(CIF_field_value, '(a)') CIF_parameter%cell_angle_beta


       case ('_cell_angle_gamma')
        READ(CIF_field_value, '(a)') CIF_parameter%cell_angle_gamma


       case ('_cell_volume')
        READ(CIF_field_value, '(a)') CIF_parameter%cell_volume


       case ('_cell_formula_units_z')
        READ(CIF_field_value, '(a)') CIF_parameter%formula_units_Z


       case ('_exptl_crystal_density_diffrn')
        READ(CIF_field_value, '(a)') CIF_parameter%exptl_density



       CASE ('_diffrn_reflns_number')
        READ(CIF_field_value, '(a)') CIF_parameter%diffrn_reflns_number


       case ('_reflns_number_total')
        READ(CIF_field_value, '(a)') CIF_parameter%reflns_number_total


       case ('_diffrn_reflns_av_r_equivalents')
        READ(CIF_field_value, '(a)') CIF_parameter%diffrn_reflns_av_R_equivalents


       case ('_refine_ls_number_parameters')
        READ(CIF_field_value, '(a)') CIF_parameter%refine_ls_number_parameters


       case ('_refine_ls_number_restraints')
        READ(CIF_field_value, '(a)') CIF_parameter%restraints_number


       case ('_refine_ls_wr_factor_gt')
        READ(CIF_field_value, '(a)') CIF_parameter%refine_ls_wR_factor_gt


       CASE ('_refine_ls_r_factor_gt')
        READ(CIF_field_value, '(a)') CIF_parameter%refine_ls_R_factor_gt


       CASE ('_refine_ls_r_factor_all')
        READ(CIF_field_value, '(a)') CIF_parameter%refine_ls_R_factor_all


       case ('_refine_ls_wr_factor_ref')
        READ(CIF_field_value, '(a)') CIF_parameter%refine_ls_wR_factor_ref


       CASE ('_refine_ls_goodness_of_fit_ref')
        READ(CIF_field_value, '(a)') CIF_parameter%chi2


       CASE ('_reflns_number_gt')
        READ(CIF_field_value, '(a)') CIF_parameter%reflns_number_gt


       case ('_refine_diff_density_max')
        READ(CIF_field_value, '(a)') CIF_parameter%refine_diff_density_max


       case ('_refine_diff_density_min')
        READ(CIF_field_value, '(a)') CIF_parameter%refine_diff_density_min


       case ('_diffrn_measurement_device_type')
        READ(CIF_field_value, '(a)') CIF_parameter%diffracto_device
        IF(CIF_parameter%diffracto_device(1:6) == 'APEXII') then
         device_mark =  'Bruker-AXS'
         device_type =  'APEXII Kappa-CCD diffractometer'
        elseif(CIF_parameter%diffracto_device(1:8) == 'KappaCCD') then
         device_mark = 'Nonius'
         device_type = 'KappaCCD'
        ELSEIF(CIF_parameter%diffracto_device(1:11) == 'CCD Saphire') then
         device_mark = 'Oxford Diffraction'
         device_type = 'Xcalibur Saphire 3'
        endif


       case ('_diffrn_ambient_temperature')
        READ(CIF_field_value, '(a)') CIF_parameter%diffracto_temperature


       case ('_diffrn_radiation_type')
        READ(CIF_field_value, '(a)') CIF_parameter%diffracto_radiation_type

       case ('_diffrn_radiation_source')
        READ(CIF_field_value, '(a)') CIF_parameter%diffracto_radiation_source

       case ('_diffrn_radiation_wavelength')
        READ(CIF_field_value, '(a)') CIF_parameter%diffracto_radiation_wavelength


       CASE ('_refine_ls_hydrogen_treatment')
        READ(CIF_field_value, '(a)') CIF_parameter%H_treatment


       case ('_diffrn_reflns_theta_min')
        READ(CIF_field_value, '(a)') CIF_parameter%theta_min


       case ('_diffrn_reflns_theta_max')
        READ(CIF_field_value, '(a)') CIF_parameter%theta_max


       case ('_diffrn_reflns_limit_h_min')
        READ(CIF_field_value, '(a)') CIF_parameter%h_min
       case ('_diffrn_reflns_limit_h_max')
        READ(CIF_field_value, '(a)') CIF_parameter%h_max
       case ('_diffrn_reflns_limit_k_min')
        READ(CIF_field_value, '(a)') CIF_parameter%k_min
       case ('_diffrn_reflns_limit_k_max')
        READ(CIF_field_value, '(a)') CIF_parameter%k_max
       case ('_diffrn_reflns_limit_l_min')
        READ(CIF_field_value, '(a)') CIF_parameter%l_min
       case ('_diffrn_reflns_limit_l_max')
        READ(CIF_field_value, '(a)') CIF_parameter%l_max


       case ('_exptl_crystal_size_max')
        READ(CIF_field_value, '(a)') CIF_parameter%crystal_size_max
       case ('_exptl_crystal_size_mid')
        READ(CIF_field_value, '(a)') CIF_parameter%crystal_size_mid
       case ('_exptl_crystal_size_min')
        READ(CIF_field_value, '(a)') CIF_parameter%crystal_size_min

       case ('_exptl_crystal_colour')
        READ(CIF_field_value, '(a)') CIF_parameter%crystal_colour

       case ('_exptl_crystal_f_000')
        READ(CIF_field_value, '(a)') CIF_parameter%F000


       case ('_diffrn_measured_fraction_theta_max')
        READ(CIF_field_value, '(a)') CIF_parameter%completeness
        IF(LEN_TRIM(CIF_parameter%completeness) == 0) then
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         IF(i_error /=0) then
          call write_info('')
          call write_info(' !! Error reading archive.CIF file !!')
         call write_info('')
         return
         endif
         READ(CIF_input_line, '(a)') CIF_parameter%completeness
         CIF_parameter%completeness = ADJUSTL(CIF_parameter%completeness)
        endif


       case ('_exptl_absorpt_coefficient_mu')
        READ(CIF_field_value, '(a)') CIF_parameter%exptl_mu


       case ('_exptl_absorpt_correction_type')
        READ(CIF_field_value, '(a)') CIF_parameter%absorption_correction_type


       case ('_exptl_absorpt_correction_t_min')
        READ(CIF_field_value, '(a)') CIF_parameter%T_min

       case ('_exptl_absorpt_correction_t_max')
        READ(CIF_field_value, '(a)') CIF_parameter%T_max


       case default
   end select


  END do
 close (UNIT = CIF_unit)


 i1 = INDEX(archive_cif, '.', back=.true.)
 WRITE(structural_report_file, '(2a)') archive_CIF(1:i1-1), '_structural_report.HTML'
 open (UNIT=HTML_unit, FILE=TRIM(structural_report_file))
 !OPEN(UNIT=HTML_unit, FILE='structural_report.HTML')
  WRITE(HTML_unit, '(a)') "<!-- Structural report (HTML file) created by CRYSCAL.exe from 'archive.cif' file."
  WRITE(HTML_unit, '(a)')  "     Web site : http://www.cdifx.univ-rennes/cryscal.htm"
  WRITE(HTML_unit, '(5a)') "     Version : ", trim(cryscal_version), "(", trim(cryscal_author), ")"
  WRITE(HTML_unit, '(a)') "!-->"

  WRITE(HTML_unit, '(a)') "<HTML>"
  WRITE(HTML_unit, '(a)') "<HEAD>"
  call write_HTML_css(HTML_unit)
  !WRITE(HTML_unit, '(a)') "<style TYPE='text/css'>"
  !WRITE(HTML_unit, '(a)') "    h2         { font-family:'Trebuchet MS',Arial; color:blue; margin-left:50px; }"
  !WRITE(HTML_unit, '(a)') "    .title_0   { font-family:'Trebuchet MS',Arial; color:#AA0000; text-align:center; margin-left:50px; font-size:24px; font-weight:bold; }"
  !WRITE(HTML_unit, '(a)') "    .title_1   { font-family:'Trebuchet MS',Arial; color:#880000; margin-left:50px; font-size:20px; background-color:#f5f5f5; }"
  !WRITE(HTML_unit, '(a)') "    .retrait_1 { font-family:'Trebuchet MS',Arial; margin-left:50px; text-align:justify; line-height:18.0pt;  }"
  !WRITE(HTML_unit, '(a)') "    .retrait_2 { margin-left:80px; }"
  !WRITE(HTML_unit, '(a)') "</style>"
  WRITE(HTML_unit, '(a)') "<TITLE>Structural report</TITLE>"
  WRITE(HTML_unit, '(a)') "</HEAD>"
  WRITE(HTML_unit, '(a)') "<BODY BGCOLOR='#FFFFFF' style='font-family:Times new roman; font-size:14; line-height:150%'>"

  WRITE(HTML_unit, '(a)') "<br>"
  WRITE(HTML_unit, '(a)') "<p class='title_main'>Crystal structure report</p><br>"
  WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;X-ray crystallographic study</p>"


  WRITE(HTMl_unit, '(a)') "<p class='retrait_1'>"
  WRITE(HTML_unit, '(a)')  ""
  if(CIF_parameter%formula_moiety == '?' .or. len_trim(CIF_parameter%formula_moiety) ==0)  &
     CIF_parameter%formula_moiety = CIF_parameter%formula_sum

  CALL transf_moiety_string("HTML", CIF_parameter%formula_moiety , HTML_string)

  !WRITE(HTML_unit, '(5a)', advance='NO') "(",TRIM(CIF_parameter%formula_moiety),"); M = ",TRIM(CIF_parameter%formula_weight),"."
  WRITE(HTML_unit, '(5a)', advance='NO') "(",TRIM(HTML_string),"); <i>M</i> = ",TRIM(CIF_parameter%formula_weight),". "
  long = len_trim(CIF_parameter%diffrn_measurement_device_type)
  if(CIF_parameter%diffrn_measurement_device_type(1:1)== "'"           .and. &
     CIF_parameter%diffrn_measurement_device_type(long:long)== "'") then
   WRITE(HTML_unit, '(4a)') CIF_parameter%diffrn_measurement_device_type(2:long-1), " diffractometer,"
  else
   WRITE(HTML_unit, '(4a)') TRIM(CIF_parameter%diffrn_measurement_device_type), " diffractometer,"
  endif


  IF(CIF_parameter%diffracto_radiation_type(1:5) == 'MoK\a') then
  WRITE(HTML_unit, '(a)', advance='NO') "Mo-K&alpha; radiation "
  else
  WRITE(HTML_unit, '(a)', advance='NO') TRIM(CIF_parameter%diffracto_radiation_type), " radiation"
  endif
  WRITE(HTML_unit, '(3a)') "(&lambda; = ", TRIM(CIF_parameter%diffracto_radiation_wavelength)," Å),"
  WRITE(HTML_unit, '(3a)', advance='NO') "<i>T</i> = ", TRIM(CIF_parameter%diffracto_temperature), " K; "

  alpha_car   = .false.
  num_car     = .false.
  HTML_string = ''
  long = len_trim(HTML_string)
  do i = 1, len_trim(CIF_parameter%symmetry_space_group)
    call check_character (CIF_parameter%symmetry_space_group(i:i), alpha_car, num_car)
    !if(alpha_car) then
    ! write(HTML_string, '(4a)') trim(HTML_string), '<i>', CIF_parameter%symmetry_space_group(i:i), '</i>'
    !else
    ! write(HTML_string, '(2a)') trim(HTML_string),  CIF_parameter%symmetry_space_group(i:i)
    !endif

    if(len_trim(CIF_parameter%symmetry_space_group(i:i)) == 0) then
     long = len_trim(HTML_string)
     HTML_string =  HTML_string(1:long) // ' '
    else
     long = len_trim(HTML_string)
     if(alpha_car) then
      HTML_string = HTML_string(1:long) // '<i>' // CIF_parameter%symmetry_space_group(i:i) // '</i>'
     else
      HTML_string = HTML_string(1:long) // CIF_parameter%symmetry_space_group(i:i)
     endif
    endif
  end do
  !write(new_HTML_string, '(3a)') HTML_string(1:8) , ' ' , trim(HTML_string(9:))
  !HTML_string = replace_car(HTML_string, "_", " ")
  if(index(HTML_string, "21") /= 0)   HTML_string = replace_car(HTML_string, "21", "2<sub>1</sub>")
  if(index(HTML_string, "31") /= 0)   HTML_string = replace_car(HTML_string, "31", "3<sub>1</sub>")
  if(index(HTML_string, "32") /= 0)   HTML_string = replace_car(HTML_string, "32", "3<sub>2</sub>")
  if(index(HTML_string, "41") /= 0)   HTML_string = replace_car(HTML_string, "41", "4<sub>1</sub>")
  if(index(HTML_string, "42") /= 0)   HTML_string = replace_car(HTML_string, "42", "4<sub>2</sub>")
  if(index(HTML_string, "43") /= 0)   HTML_string = replace_car(HTML_string, "43", "4<sub>3</sub>")
  if(index(HTML_string, "61") /= 0)   HTML_string = replace_car(HTML_string, "61", "6<sub>1</sub>")
  if(index(HTML_string, "62") /= 0)   HTML_string = replace_car(HTML_string, "62", "6<sub>2</sub>")
  if(index(HTML_string, "63") /= 0)   HTML_string = replace_car(HTML_string, "63", "6<sub>3</sub>")
  if(index(HTML_string, "64") /= 0)   HTML_string = replace_car(HTML_string, "64", "6<sub>4</sub>")
  if(index(HTML_string, "65") /= 0)   HTML_string = replace_car(HTML_string, "65", "6<sub>5</sub>")

  WRITE(HTML_unit, '(6a)') TRIM(CIF_parameter%symmetry_cell_setting), " ", trim(HTML_string),  &
                           " (I.T.#", TRIM(CIF_parameter%symmetry_IT_number), "), "

  IF(CIF_parameter%symmetry_cell_setting(1:9) == 'triclinic') then
   WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
   WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
   WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
   WRITE(HTML_unit, '(3a)') "&alpha; = ",TRIM(CIF_parameter%cell_angle_alpha), ", "
   WRITE(HTML_unit, '(3a)') "&beta; = ",TRIM(CIF_parameter%cell_angle_beta),  ", "
   WRITE(HTML_unit, '(3a)') "&gamma = ",TRIM(CIF_parameter%cell_angle_gamma), " °, "
  ELSEIF(CIF_parameter%symmetry_cell_setting(1:10) == 'monoclinic') then
   WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
   WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
   WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
   WRITE(HTML_unit, '(3a)') "&beta; = ",TRIM(CIF_parameter%cell_angle_beta), " °, "
  ELSEIF(CIF_parameter%symmetry_cell_setting(1:12) == 'orthorhombic') then
   WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
   WRITE(HTML_unit, '(3a)') "b = ",TRIM(CIF_parameter%cell_length_b), ", "
   WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
  ELSEIF(CIF_parameter%symmetry_cell_setting(1:10) == 'tetragonal') then
   WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
   WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
  ELSEIF(CIF_parameter%symmetry_cell_setting(1:9) == 'hexagonal' .or.    &
         CIF_parameter%symmetry_cell_setting(1:8) == 'trigonal') then
   WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), ", "
   WRITE(HTML_unit, '(3a)') "c = ",TRIM(CIF_parameter%cell_length_c), " Å, "
  ELSEIF(CIF_parameter%symmetry_cell_setting(1:11) == 'rhomboedral') then
   WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_a), " Å, "
   WRITE(HTML_unit, '(3a)') "&alpha; = ",TRIM(CIF_parameter%cell_angle_alpha), " °, "
  ELSEIF(CIF_parameter%symmetry_cell_setting(1:5) == 'cubic') then
   WRITE(HTML_unit, '(3a)') "a = ",TRIM(CIF_parameter%cell_length_c), " Å, "
  endif

  WRITE(HTML_unit, '(3a)', advance='NO') "<i>V</i> = ",TRIM(CIF_parameter%cell_volume), " Å<sup>3</sup>, "
  WRITE(HTML_unit, '(3a)', advance='NO') "<i>Z</i> = ",TRIM(CIF_parameter%formula_units_Z), ", "
  WRITE(HTML_unit, '(3a)') "<i>d</i> = ",TRIM(CIF_parameter%exptl_density), " g.cm<sup>-3</sup>, "
  WRITE(HTML_unit, '(3a)') "&mu; = ",TRIM(CIF_parameter%exptl_mu), " mm<sup>-1</sup>. "


  !WRITE(HTML_unit, '(a)') "The structure was solved by direct methods using the SIR97 program [1], "
  WRITE(HTML_unit, '(3a)') "The structure was solved by direct methods using the <i>",  &
                            u_case(trim(structure_solution%name)), "</i> program [1], "
  WRITE(HTML_unit, '(3a)') "and then refined with full-matrix least-square methods based on <i>F</i><sup>2</sup> (<i>", &
                            u_case(trim(structure_refinement%name)), "</i>) [2] "
  WRITE(HTML_unit, '(a)')  "with the aid of the <i>WINGX</i> [3] program."


  SQUEEZE%procedure = .false.
  inquire(file= trim(SQUEEZE%file), exist= file_exist)
  if(file_exist) SQUEEZE%procedure = .true.
  if(SQUEEZE%procedure) then
   WRITE(HTML_unit, '(a)') "The contribution of the disordered solvents to the calculated structure factors was "
   WRITE(HTML_unit, '(3a)') "estimated following the <i>BYPASS</i> algorithm [4], implemented as the <i>SQUEEZE</i> ", &
                           "option in <i>PLATON</i> [5]. A new data set, free of solvent contribution, was then ",    &
                           "used in the final refinement."
  endif

  WRITE(HTML_unit, '(a)') "All non-hydrogen atoms were refined with anisotropic atomic displacement parameters. "
  IF(CIF_parameter%H_treatment(1:5) == 'mixed') then
   WRITE(HTML_unit, '(2a)') "Except <span style='color:red'><b>XXX</b></span> linked hydrogen atoms that were introduced ", &
                           "in the structural model through Fourier difference maps analysis,"
  endif
  WRITE(HTML_unit, '(2a)') "H atoms were finally included in their calculated positions. A final refinement on ", &
                           "<i>F</i><sup>2</sup> with "
  WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%reflns_number_total), " unique intensities and "
  WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%refine_ls_number_parameters), &
                           " parameters converged at &omega;<i>R</i>(<i>F</i><sup>2</sup>) = "
  WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%refine_ls_wR_factor_gt), " (<i>R</i>(<i>F</i>) = "
  WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%refine_ls_R_factor_gt), ") for "
  WRITE(HTML_unit, '(2a)') TRIM(CIF_parameter%reflns_number_gt), " observed reflections with <i>I</i> > 2&sigma;(<i>I</i>)."
  WRITE(HTML_unit, '(a)') "</p>"


  WRITE(HTMl_unit, '(a)') "<p class='retrait_1'>"
  !WRITE(HTML_unit, '(a)') "[1] &nbsp; A. Altomare, M. C. Burla, M. Camalli, G. Cascarano, C. Giacovazzo, A. Guagliardi, "
  !WRITE(HTML_unit, '(a)') "A. G. G. Moliterni, G. Polidori, R. Spagna, J. Appl. Cryst. (1999) 32, 115-119. <br>"
  !WRITE(HTML_unit, '(a)') "[2] &nbsp; SHELX97 - Programs for Crystal Structure Analysis (Release 97-2). G. M. Sheldrick, "
  !WRITE(HTMl_unit, '(a)') "Institüt für Anorganische Chemie der Universität, Tammanstrasse 4, D-3400 Göttingen, Germany, 1998<br>"
  !WRITE(HTML_unit, '(a)') "[2] &nbsp; Sheldrick, G. M. (2008) Acta Cryst. A64 112-122<br>"


  WRITE(HTML_unit, '(3a)') "[1] &nbsp; ", trim(structure_solution%reference),   "<br>"
  WRITE(HTML_unit, '(3a)') "[2] &nbsp; ", trim(structure_refinement%reference), "<br>"

  WRITE(HTML_unit, '(a)') "[3] &nbsp; L. J. Farrugia, J. Appl. Cryst., 1999, 32, 837-838<br>"
  if(SQUEEZE%procedure) then
  WRITE(HTML_unit, '(a)') "[4] &nbsp; P. v.d. Sluis and A.L. Spek, Acta Cryst. (1990) A46, 194-201<br>"
  WRITE(HTML_unit, '(a)') "[5] &nbsp; A. L. Spek, J. Appl. Cryst. (2003), 36, 7-13<br>"
  end if

  WRITE(HTMl_unit, '(a)') "</p>"
  WRITE(HTML_unit, '(a)') "<br>"



  WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Structural data</p>"


  WRITE(HTML_unit, '(a)')  "<pre style='font-size:14'>"
  CALL transf_moiety_string("HTML", CIF_parameter%formula_moiety , HTML_string)
  WRITE(HTML_unit, '(10x,2a)') "Empirical formula                      ", TRIM(HTML_string)
  !WRITE(HTML_unit, '(10x,2a)') "Empirical formula                      ", TRIM(CIF_parameter%formula_moiety)
  WRITE(HTML_unit, '(10x,2a)') "Formula weight                         ", TRIM(CIF_parameter%formula_weight)
  WRITE(HTML_unit, '(10x,3a)') "Temperature                            ", TRIM(CIF_parameter%diffracto_temperature)," <i>K</i>"
  WRITE(HTML_unit, '(10x,3a)') "Wavelength                             ", TRIM(CIF_parameter%diffracto_radiation_wavelength), " Å "
  HTML_string = CIF_parameter%symmetry_space_group
  if(index(HTML_string, "21") /= 0)   HTML_string = replace_car(HTML_string, "21", "2<sub>1</sub>")
  if(index(HTML_string, "31") /= 0)   HTML_string = replace_car(HTML_string, "31", "3<sub>1</sub>")
  if(index(HTML_string, "32") /= 0)   HTML_string = replace_car(HTML_string, "32", "3<sub>2</sub>")
  if(index(HTML_string, "41") /= 0)   HTML_string = replace_car(HTML_string, "41", "4<sub>1</sub>")
  if(index(HTML_string, "42") /= 0)   HTML_string = replace_car(HTML_string, "42", "4<sub>2</sub>")
  if(index(HTML_string, "43") /= 0)   HTML_string = replace_car(HTML_string, "43", "4<sub>3</sub>")
  if(index(HTML_string, "61") /= 0)   HTML_string = replace_car(HTML_string, "61", "6<sub>1</sub>")
  if(index(HTML_string, "62") /= 0)   HTML_string = replace_car(HTML_string, "62", "6<sub>2</sub>")
  if(index(HTML_string, "63") /= 0)   HTML_string = replace_car(HTML_string, "63", "6<sub>3</sub>")
  if(index(HTML_string, "64") /= 0)   HTML_string = replace_car(HTML_string, "64", "6<sub>4</sub>")
  if(index(HTML_string, "65") /= 0)   HTML_string = replace_car(HTML_string, "65", "6<sub>5</sub>")
  WRITE(HTML_unit, '(10x,6a)') "Crystal system, space group            ", TRIM(CIF_parameter%symmetry_cell_setting), ", ", &
                               "<i>",TRIM(HTML_string),"</i>"
  !WRITE(HTML_unit, '(10x,4a)') "Crystal system, space group            ", TRIM(CIF_parameter%symmetry_cell_setting), ", ", &
  !                              TRIM(CIF_parameter%symmetry_space_group)
  WRITE(HTML_unit, '(10x,6a)') "Unit cell dimensions                   ", "a = ", TRIM(CIF_parameter%cell_length_a), &
                               " Å, alpha = ", TRIM(CIF_parameter%cell_angle_alpha), " °"
  WRITE(HTML_unit, '(10x,6a)') "                                       ", "b = ", TRIM(CIF_parameter%cell_length_b), &
                               " Å, beta = ", TRIM(CIF_parameter%cell_angle_beta),  " °"
  WRITE(HTML_unit, '(10x,6a)') "                                       ", "c = ", TRIM(CIF_parameter%cell_length_c), &
                               " Å, gamma = ", TRIM(CIF_parameter%cell_angle_gamma), " °"
  WRITE(HTML_unit, '(10x,3a)') "Volume                                 ", TRIM(CIF_parameter%cell_volume), " Å<sup>3</sup>"
  WRITE(HTML_unit, '(10x,5a)') "<i>Z</i>, Calculated density                  ", TRIM(CIF_parameter%formula_units_Z), ' , ', &
                               TRIM(CIF_parameter%exptl_density), " (g.cm<sup>-3</sup>)"
  WRITE(HTML_unit, '(10x,3a)') "Absorption coefficient                 ", TRIM(CIF_parameter%exptl_mu), " mm<sup>-1</sup>"
  IF(long_report) then
  WRITE(HTML_unit, '(10x,2a)') "<i>F</i>(000)                                 ", TRIM(CIF_parameter%F000)
  WRITE(HTML_unit, '(10x,7a)') "Crystal size                           ", TRIM(CIF_parameter%crystal_size_max), " x ", &
                               TRIM(CIF_parameter%crystal_size_mid), " x ", TRIM(CIF_parameter%crystal_size_min), " mm"
  WRITE(HTML_unit, '(10x,2a)') "Crystal color                          ", TRIM(CIF_parameter%crystal_colour)
  WRITE(HTML_unit, '(10x,5a)') "Theta range for data collection        ", TRIM(CIF_parameter%theta_min), " to ", &
                                TRIM(CIF_parameter%theta_max), " °"
  WRITE(HTML_unit, '(10x,4a)') "h_min, h_max                           ", TRIM(CIF_parameter%h_min), " , ", &
                                TRIM(CIF_parameter%h_max)
  WRITE(HTML_unit, '(10x,4a)') "k_min, k_max                           ", TRIM(CIF_parameter%k_min), " , ", &
                                TRIM(CIF_parameter%k_max)
  WRITE(HTML_unit, '(10x,4a)') "l_min, l_max                           ", TRIM(CIF_parameter%l_min), " , ", &
                                TRIM(CIF_parameter%l_max)
  WRITE(HTMl_unit, '(10x,7a)') "Reflections collected / unique         ", TRIM(CIF_parameter%diffrn_reflns_number), &
                               " / ", TRIM(CIF_parameter%reflns_number_total), " [R(int) = ",                       &
                               TRIM(CIF_parameter%diffrn_reflns_av_R_equivalents), "]"
  WRITE(HTML_unit, '(10x,2a)') "Completeness to theta_max              ", TRIM(CIF_parameter%completeness)
  WRITE(HTML_unit, '(10x,2a)') "Absorption correction type             ", TRIM(CIF_parameter%absorption_correction_type)
  if(CIF_parameter%absorption_correction_type(1:4) /= 'none') then
  WRITE(HTML_unit, '(10x,4a)') "Max. and min. transmission             ", TRIM(CIF_parameter%T_max), " , ", &
                               TRIM(CIF_parameter%T_min)
  endif
  WRITE(HTML_unit, '(10x,3a)') "Refinement method                      ", "Full-matrix least-squares on <i>F</i><sup>2</sup>"
  endif
  IF(long_report) then
  WRITE(HTML_unit, '(10x,6a)') "Data / restraints / parameters         ", TRIM(CIF_parameter%reflns_number_total), " / ", &
                               TRIM(CIF_parameter%restraints_number), " / ", TRIM(CIF_parameter%refine_ls_number_parameters)
  WRITE(HTML_unit, '(10x,2a)') "Goodness-of-fit                        ", TRIM(CIF_parameter%Chi2)
  else
  WRITE(HTML_unit, '(10x,2a)') "Unique reflections                     ", TRIM(CIF_parameter%reflns_number_total)
  WRITE(HTML_unit, '(10x,2a)') "Refined parameters                     ", TRIM(CIF_parameter%refine_ls_number_parameters)
  endif
  !WRITE(HTML_unit, '(10x,5a)') "Final R indices [I>2sigma(I)]          ", "<i>R</i>1<sup><i>a</i></sup> = ",     &
  !                              TRIM(CIF_parameter%refine_ls_R_factor_gt),  ", <i>wR</i>2<sup><i>b</i></sup> = ", &
  !                              TRIM(CIF_parameter%refine_ls_wR_factor_gt)
  WRITE(HTML_unit, '(10x,5a)') "Final <i>R</i> indices [I>2&sigma;]                 ", "<i>R</i>1<sup><i>a</i></sup> = ",  &
                                TRIM(CIF_parameter%refine_ls_R_factor_gt),  ", <i>wR</i>2<sup><i>b</i></sup> = ",          &
                                TRIM(CIF_parameter%refine_ls_wR_factor_gt)
  WRITE(HTML_unit, '(10x,5a)') "<i>R</i> indices (all data)                   ", "<i>R</i>1<sup><i>a</i></sup> = ",        &
                               TRIM(CIF_parameter%refine_ls_R_factor_all), ", <i>wR</i>2<sup><i>b</i></sup> = ",           &
                               TRIM(CIF_parameter%refine_ls_wR_factor_ref)
  WRITE(HTML_unit, '(10x,5a)') "Largest diff. peak and hole            ", TRIM(CIF_parameter%refine_diff_density_max),     &
                               " and ", TRIM(CIF_parameter%refine_diff_density_min), "  e.Å<sup>-3</sup>"
  WRITE(HTML_unit, '(a)')      ""
  WRITE(HTML_unit, '(12x,2a)')      "<sup><i>a</i></sup><i>R</i>1 = &#8721; &#124; &#124;<i>F<sub>o</sub></i>&#124; - &#124;", &
                                    "<i>F<sub>c</sub></i>&#124; &#124; / &#8721; &#124;<i>F<sub>o</sub></i>&#124;<br>"
  WRITE(HTML_unit, '(12x,3a)')      "<sup><i>b</i></sup><i>wR</i>2 = {&#8721; [<i>w</i>(<i>F<sub>o</sub></i><sup>2</sup> - ",  &
                                    " <i>F<sub>c</sub></i><sup>2</sup>)<sup>2</sup>] / &#8721; ",                              &
                                    "[<i>w</i>(<i>F<sub>o</sub></i><sup>2</sup>)<sup>2</sup>]}<sup>1/2</sup>"
  WRITE(HTML_unit, '(a)')  "</pre>"
  !WRITE(HTML_unit, '(a)') "<p class='retrait_3'>"
  !WRITE(HTML_unit, '(12x,2a)')     "<sup><i>a</i></sup><i>R</i>1 = &#8721; &#124; &#124;<i>F<sub>o</sub></i>&#124; - &#124;", &
  !                                 "<i>F<sub>c</sub></i>&#124; &#124; / &#8721; &#124;<i>F<sub>o</sub></i>&#124;<br>"
  !WRITE(HTML_unit, '(12x,3a)')     "<sup><i>b</i></sup><i>wR</i>2 = {&#8721; [<i>w</i>(<i>F<sub>o</sub></i><sup>2</sup> -  ", &
  !                                 "<i>F<sub>c</sub></i><sup>2</sup>)<sup>2</sup>] / &#8721; ",                               &
  !                                 "[<i>w</i>(<i>F<sub>o</sub></i><sup>2</sup>)<sup>2</sup>]}<sup>1/2</sup>"
  !write(HTML_unit, '(a)')  "</p>"


!--------------------------------------------------------


  WRITE(HTML_unit, '(a)') "<br>"


  IF(long_report) then
   ! positions atomiques
   WRITE(HTML_unit, '(a)')  "<p class='title_2'>"
   WRITE(HTML_unit, '(2a)') "&nbsp;&nbsp;Atomic coordinates and equivalent isotropic displacement parameters ", &
                            "(Å<sup>2</sup> x 10<sup>3</sup>). "
   WRITE(HTML_unit, '(a)')  "U(eq) is defined as one third of the trace of the orthogonalized U<sub>ij</sub> tensor."
   WRITE(HTML_unit, '(a)')  "</p>"


!---------------------------------------------------------------------
   WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"


   OPEN(UNIT = CIF_unit, FILE=TRIM(archive_cif))
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0)            exit
     IF(LEN_TRIM(CIF_input_line)==0) cycle
     READ(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
         case ('_atom_site_label')
          WRITE(HTML_unit, '(10x,2a)') "Atom           x              y              z              U(eq)"
          WRITE(HTMl_unit, '(10x,a)')  ""
          do
           READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
           IF(i_error < 0) exit
           CIF_input_line = ADJUSTL(CIF_input_line)
           IF(CIF_input_line(1:1) == '_')     cycle
           IF(LEN_TRIM(CIF_input_line) == 0)  exit
           IF(CIF_input_line(1:5) == 'loop_') exit
           IF(CIF_input_line(1:1) == ';')     exit
           IF(CIF_input_line(1:1) == '#')     exit
           do
            CIF_input_line = TRIM(CIF_input_line)
            long= LEN_TRIM(CIF_input_line)
            IF(CIF_input_line(long:long) == '?') then
             CIF_input_line = CIF_input_line(1:long-1)
            ELSEIF(CIF_input_line(long:long) == '.') then
             CIF_input_line = CIF_input_line(1:long-1)
            else
             exit
            endif
           enddo
           CIF_parameter%atom = CIF_input_line
           call nombre_de_colonnes(CIF_parameter%atom, nb_col)
           if(nb_col ==5) then
            READ(CIF_parameter%atom, *) atom_label, atom_symbol, atom_x, atom_y, atom_z
            atom_Ueq = '?'
           else
            READ(CIF_parameter%atom, *) atom_label, atom_symbol, atom_x, atom_y, atom_z, atom_Ueq
           endif
           WRITE(HTML_unit, '(10x,5(a,3x))') atom_label(1:12), atom_x(1:12), atom_y(1:12), atom_z(1:12), atom_Ueq(1:12)
          END do


         case default
     end select
    end do
    close (UNIT = CIF_unit)
    WRITE(HTML_unit, '(a)')  "</pre>"
    WRITE(HTML_unit, '(a)')  "<br>"




  ! Anisotropic displacement parameters
  ! WRITE(HTML_unit, '(a)') "<p class='title_2'>Anisotropic displacement parameters (Å<sup>2</sup>)</p>"
  ! WRITE(HTML_unit, '(a)') "<p class='retrait_1'>The anisotropic displacement factor exponent takes the form:"
  ! WRITE(HTML_unit, '(2a)') "-2&pi;<sup>2</sup> [ h<sup>2</sup> a*<sup>2</sup> U<sub>11</sub> + ... + ", &
  !                          "2 h k a* b* U<sub>12</sub> ]."
  ! WRITE(HTML_unit, '(a)') "</p>"


!-----------------------------------------------------------------------------
   WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"
   OPEN(UNIT = CIF_unit, FILE=TRIM(archive_cif))
    n = 0
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0)            exit
     IF(LEN_TRIM(CIF_input_line)==0) cycle
     READ(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
         case ('_atom_site_aniso_label')
          !WRITE(HTML_unit, '(10x,3a)') "Atom           U11            U22            U33            U23            U13            U12"
          !WRITE(HTML_unit, '(a)')  ""


          do
           READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
           IF(i_error < 0) exit
           CIF_input_line = ADJUSTL(CIF_input_line)
           IF(CIF_input_line(1:1) == '_')     cycle
           IF(LEN_TRIM(CIF_input_line) == 0)  exit
           IF(CIF_input_line(1:5) == 'loop_') exit
           IF(CIF_input_line(1:1) == ';')     exit
           IF(CIF_input_line(1:1) == '#')     exit
           n = n + 1
           if (n == 1) then
            WRITE(HTML_unit, '(a)')  "<p class='title_2'>&nbsp;&nbsp;Anisotropic displacement parameters (Å<sup>2</sup>)</p>"
            WRITE(HTML_unit, '(a)')  "<p class='retrait_1'>The anisotropic displacement factor exponent takes the form:"
            WRITE(HTML_unit, '(2a)') "-2&pi;<sup>2</sup> [ h<sup>2</sup> a*<sup>2</sup> U<sub>11</sub> + ... ", &
                                     "+ 2 h k a* b* U<sub>12</sub> ]."
            WRITE(HTML_unit, '(a)')  "</p>"
            WRITE(HTML_unit, '(10x,2a)') "Atom           U11            U22            U33            ",   &
                                         "U23            U13            U12"
            WRITE(HTML_unit, '(a)')  ""
           endif
           CIF_parameter%atom = CIF_input_line
           READ(CIF_parameter%atom, *) atom_label, atom_U11, Atom_U22, atom_U33, atom_U23, atom_U13, atom_U12
           WRITE(HTML_unit, '(10x,7(a,3x))')  atom_label(1:12), atom_U11(1:12), atom_U22(1:12), atom_U33(1:12), &
                                              atom_U23(1:12), atom_U13(1:12), atom_U12(1:12)


          END do


         case default
     end select
    end do
    close (UNIT = CIF_unit)
    WRITE(HTML_unit, '(a)')  "</pre>"


!-----------------------------------------------------------------------------




   ! distances interatomiques
   !WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Bond lengths [Å]</p>"
   !!WRITE(HTML_unit, '(a)')  "<font face='courier new'>"
   !WRITE(HTML_unit, '(a)')  "<pre style='font-size:14'>"


   open (unit = CIF_unit, file = TRIM(archive_CIF))
   n = 0


    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0) exit
     IF(LEN_TRIM(CIF_input_line) == 0) cycle
     read(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
       case ('_geom_bond_atom_site_label_1')
        dico = ''
        n_sym_dist = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         CIF_input_line = ADJUSTL(CIF_input_line)
         IF(CIF_input_line(1:1) == '_') cycle
         IF(LEN_TRIM(CIF_input_line) == 0)  exit
         IF(CIF_input_line(1:5) == 'loop_') exit
         IF(CIF_input_line(1:1) == ';')     exit
         IF(CIF_input_line(1:1) == '#')     exit
         n = n + 1
         if(n==1) then
          WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Bond lengths [Å]</p>"
          WRITE(HTML_unit, '(a)')  "<pre style='font-size:14'>"
         end if
         CIF_parameter%distance = CIF_input_line
         READ(CIF_parameter%distance, *) dist_atom(1), dist_atom(2), dist_value, dist_sym
         dist_sym = adjustl(dist_sym)
         !WRITE(HTML_unit, '(6a)') dist_atom(1)(1:4), ' - ' , dist_atom(2)(1:4), ' = ', TRIM(dist_value), '<br>'

         HTML_str(1) = ''
         if(dist_sym(1:1) /= '.') then
          n_sym_dist = n_sym_dist + 1
          site_sym(n_sym_dist) = dist_sym
          call Get_num_site_sym(n_sym_dist, site_sym, dist_sym, num_site_sym)
          !call Get_num_site_sym_new(n_sym_dist, site_sym, dist_sym, num_site_sym, dico)


          call get_sym_op(dist_sym, op_numor, t)
          op_n(num_site_sym(n_sym_dist))   = op_numor
          op_t(num_site_sym(n_sym_dist),:) = t(:)

          call build_HTML_string(num_site_sym(n_sym_dist), HTML_str(1))
         endif
         write(HTML_unit, '(10x,6a)') dist_atom(1)(1:4), ' - ' , dist_atom(2)(1:4), trim(HTML_str(1)), ' = ', TRIM(dist_value)
        end do


       case default
     end select
    END do
   close (UNIT = CIF_unit)
   !WRITE(HTML_unit, '(a)')  "</font>"
   WRITE(HTML_unit, '(a)')  "</pre>"
   WRITE(HTML_unit, '(a)')  "<br>"

   if(n_sym_dist /=0) then
    WRITE(HTML_unit, '(a)')  "</pre>"
    write(HTML_unit, '(a)')  ""
    write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
    WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"
    call write_HTML_sym(maxval(num_site_sym),op_string, op_n, op_t)
   endif
   WRITE(HTML_unit, '(a)')  "</pre></ul>"




!-------------------------------------------------------------------------
   ! angles interatomiques
   !WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Angles [°]</p>"
   !!WRITE(HTML_unit, '(a)') "<font face='courier new'>"
   !WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"



   open (unit = CIF_unit, file = TRIM(archive_CIF))
   n = 0
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0) exit
     IF(LEN_TRIM(CIF_input_line) == 0) cycle
     read(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
       case ('_geom_angle_atom_site_label_1')
        n_sym_ang    = 0
        num_site_sym = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         CIF_input_line = ADJUSTL(CIF_input_line)
         IF(CIF_input_line(1:1) == '_')     cycle
         IF(LEN_TRIM(CIF_input_line) == 0)  exit
         IF(CIF_input_line(1:5) == 'loop_') exit
         IF(CIF_input_line(1:1) == ';')     exit
         IF(CIF_input_line(1:1) == '#')     exit
         n = n + 1
         CIF_parameter%angle = CIF_input_line
         READ(CIF_parameter%angle, *) ang_atom(1), ang_atom(2), ang_atom(3), ang_value, ang_sym(1), ang_sym(2)
         ang_sym(1:2) = adjustl(ang_sym(1:2))
         if(n==1) then
          WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Angles [°]</p>"
          !WRITE(HTML_unit, '(a)') "<font face='courier new'>"
          WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"
         end if

         do i=1, 2
          HTML_str(i) = ''
          if(ang_sym(i)(1:1) /= '.') then
           n_sym_ang = n_sym_ang + 1
           site_sym(n_sym_ang) = ang_sym(i)
           call Get_num_site_sym(n_sym_ang, site_sym, ang_sym(i), num_site_sym)
           call get_sym_op(ang_sym(i), op_numor, t)
           op_n(num_site_sym(n_sym_ang))   = op_numor
           op_t(num_site_sym(n_sym_ang),:) = t(:)

           call build_HTML_string(num_site_sym(n_sym_ang), HTML_str(i))
          endif
         end do

         write(HTML_unit, '(10x,9a)') ang_atom(1)(1:4), trim(HTML_str(1)), ' - ', ang_atom(2)(1:4), ' - ', &
                                      ang_atom(3)(1:4), trim(HTML_str(2)), ' = ', trim(ang_value)

        end do
       case default
     end select
    END do
   close (UNIT = CIF_unit)
   !WRITE(HTML_unit, '(a)')  "</font>"
   WRITE(HTML_unit, '(a)')  "</pre>"
   WRITE(HTML_unit, '(a)')  "<br>"

  if(n_sym_ang /=0) then
    WRITE(HTML_unit, '(a)')  "</pre>"
    write(HTML_unit, '(a)')  ""
    write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
    WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"

    call write_HTML_sym(maxval(num_site_sym), op_string, op_n, op_t)

   endif
   WRITE(HTML_unit, '(a)')  "</pre></ul>"


!-------------------------------------------------------------------------
   ! angles de torsion interatomiques
   !WRITE(HTML_unit, '(a)') "<p class='title_2">&nbsp;&nbsp;Torsion angles [°]</p>"
   !!WRITE(HTML_unit, '(a)') "<font face='courier new'>"
   !WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"



   open (unit = CIF_unit, file = TRIM(archive_CIF))
   n = 0
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0) exit
     IF(LEN_TRIM(CIF_input_line) == 0) cycle
     read(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
       case ('_geom_torsion_atom_site_label_1')
        n_sym_tor = 0
        num_site_sym = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         CIF_input_line = ADJUSTL(CIF_input_line)
         IF(CIF_input_line(1:1) == '_') cycle
         IF(LEN_TRIM(CIF_input_line) == 0)  exit
         IF(CIF_input_line(1:5) == 'loop_') exit
         IF(CIF_input_line(1:1) == ';')     exit
         IF(CIF_input_line(1:1) == '#')     exit

         n = n + 1
         if(n==1) then
          WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Torsion angles [°]</p>"
          !WRITE(HTML_unit, '(a)') "<font face='courier new'>"
          WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"
         endif
         CIF_parameter%torsion_angle = CIF_input_line
         READ(CIF_parameter%torsion_angle, *) (torsion_ang_atom(i), i=1,4), torsion_ang_value, (torsion_sym(i), i=1,4)
         torsion_sym(1:4) = adjustl(torsion_sym(1:4))

         do i=1, 4
          HTML_str(i) = ''
          if(torsion_sym(i) /= '.') then
           n_sym_tor = n_sym_tor + 1
           site_sym(n_sym_tor) = torsion_sym(i)
           call Get_num_site_sym(n_sym_tor, site_sym, torsion_sym(i), num_site_sym)
           call get_sym_op(torsion_sym(i), op_numor, t)
           op_n(num_site_sym(n_sym_tor))   = op_numor
           op_t(num_site_sym(n_sym_tor),:) = t(:)
           call build_HTML_string(num_site_sym(n_sym_tor), HTML_str(i))
          endif
         end do

         WRITE(HTML_unit, '(10x,13a)') torsion_ang_atom(1)(1:4), trim(HTML_str(1)), ' - ',  &
                                       torsion_ang_atom(2)(1:4), trim(HTML_str(2)), ' - ',  &
                                       torsion_ang_atom(3)(1:4), trim(HTML_str(3)), ' - ',  &
                                       torsion_ang_atom(4)(1:4), trim(HTML_str(4)), ' = ',  TRIM(torsion_ang_value)
        end do


       case default
     end select
    END do
   close (UNIT = CIF_unit)
   !WRITE(HTML_unit, '(a)')  "</font>"
   WRITE(HTML_unit, '(a)')  "</pre>"

  if(n_sym_tor /=0) then
    WRITE(HTML_unit, '(a)')  "</pre>"
    write(HTML_unit, '(a)')  ""
    write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
    WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"
    call write_HTML_sym(maxval(num_site_sym), op_string, op_n, op_t)
   endif
   WRITE(HTML_unit, '(a)')  "</pre></ul>"


!-------------------------------------------------------------------------
   ! liaisons H



   open (unit = CIF_unit, file = TRIM(archive_CIF))
    do
     READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
     IF(i_error < 0) exit
     IF(LEN_TRIM(CIF_input_line) == 0) cycle
     read(CIF_input_line, *, IOSTAT=i_error) CIF_field
     IF(i_error /=0) then
      call write_info('')
      call write_info(' !! Error reading archive.CIF file !!')
      call write_info('')
      return
     endif


     long_field = LEN_TRIM(CIF_field)
     select case (CIF_field(1:long_field))
       case ('_geom_hbond_atom_site_label_D')
        n_sym_htab = 0
        do
         READ(CIF_unit, '(a)', IOSTAT=i_error) CIF_input_line
         IF(i_error < 0) exit
         CIF_input_line = ADJUSTL(CIF_input_line)
         IF(CIF_input_line(1:1) == '_') cycle
         IF(LEN_TRIM(CIF_input_line) == 0)  exit
         IF(CIF_input_line(1:5) == 'loop_') exit
         IF(CIF_input_line(1:1) == ';')     exit
         IF(CIF_input_line(1:1) == '#')     exit
         CIF_parameter%Hbond = CIF_input_line
         READ(CIF_parameter%Hbond, *) site_label_D, site_label_H, site_label_A, dist_DH, dist_HA, dist_DA, angle_DHA, &
                                      site_sym_A
         call get_sym_op(site_sym_A, op_numor, t)
         n_sym_htab = n_sym_htab + 1
         op_n(n_sym_htab) = op_numor
         op_t(n_sym_htab,:) = t(:)


         if(n_sym_htab ==1) then
          WRITE(HTML_unit, '(a)') "<br>"
          WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Hydrogen bonds [A and deg.] </p>"
          WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"
         endif
         if(op_numor < 10) then
         WRITE(HTML_unit, '(10x,6a,i1,5a)')  site_label_D(1:4), ' - ', site_label_H(1:4), '.... ', site_label_A(1:4), "#", &
                                             n_sym_htab, '   ', dist_DH(1:12), dist_HA(1:12), dist_DA(1:12), angle_DHA(1:12)
         elseif(op_numor < 100) then
         WRITE(HTML_unit, '(10x,6a,i2,5a)')  site_label_D(1:4), ' - ', site_label_H(1:4), '.... ', site_label_A(1:4), "#", &
                                             n_sym_htab, '   ', dist_DH(1:12), dist_HA(1:12), dist_DA(1:12), angle_DHA(1:12)
         else
         WRITE(HTML_unit, '(10x,6a,i3,5a)')  site_label_D(1:4), ' - ', site_label_H(1:4), '.... ', site_label_A(1:4), "#", &
                                             n_sym_htab, '   ', dist_DH(1:12), dist_HA(1:12), dist_DA(1:12), angle_DHA(1:12)
         endif



        end do


       case default
     end select
    END do
   close (UNIT = CIF_unit)

   if(n_sym_htab /=0) then
    WRITE(HTML_unit, '(a)')  "</pre>"
    write(HTML_unit, '(a)')  ""
    write(HTML_unit, '(a)')  "<p class='retrait_1'>Symmetry transformations used to generate equivalent atoms:</p>"
    WRITE(HTML_unit, '(a)') "<pre style='font-size:14'>"
    call write_HTML_sym(n_sym_htab, op_string, op_n, op_t)
   endif
   WRITE(HTML_unit, '(a)')  "</pre></ul>"



!-------------------------------------------------------------------------
  endif

  ! recherche de la presence de fichier.BMP;JPG;PNG
  !WRITE(HTML_unit, '(a)') '<h2><blockquote><br>Structure visualisation</h2>'
  !WRITE(HTML_unit, '(a)') '<blockquote><img src="job.png">'
  !WRITE(HTML_unit, '(a)') "</blockquote></blockquote>"


  file_exist = .false.
  GIF_file = "platon_ortep.gif"
  call test_file_exist(trim(GIF_file), file_exist)
  if(.not. file_exist) then
   GIF_file   = "ortep.gif"
   call test_file_exist(trim(GIF_file), file_exist)
   if(.not. file_exist) then
    GIF_file = "platon_jobte.gif"
    call test_file_exist(trim(GIF_file), file_exist)
   endif
  endif

  if(file_exist) then
    WRITE(HTML_unit, '(a)') ""
    WRITE(HTML_unit, '(a)') "<br>"
    WRITE(HTML_unit, '(a)') "<p class='title_2'>&nbsp;&nbsp;Structure visualisation</p>"
    WRITE(HTML_unit, '(3a)') '<center><img width=600 src="', trim(GIF_file),'"></center>'
  endif




  WRITE(HTML_unit, '(a)')  "<HR>"
  WRITE(HTML_unit, '(a)')  "<p class='retrait_1'>"
  WRITE(HTML_unit, '(4a)') "<i>This HTML report has been created through ",                        &
                           "<A HREF='http://www.cdifx.univ-rennes1.fr/cryscal.htm'>CRYSCAL</A> (", &
                           trim(cryscal_version), ")."
  WRITE(HTML_unit, '(2a)') "Please report bugs and problems to <A HREF='mailto:cdifx@univ-rennes1.fr'>", &
                           "cdifx.univ-rennes1.fr</A></p>"
  WRITE(HTML_unit, '(a)')  "</BODY>"
  WRITE(HTML_unit, '(a)')  "</HTML>"
CLOSE(UNIT=HTML_unit)

  call write_info('')
  call write_info('    >>> '//TRIM(structural_report_file)//' file has been created.')
  call write_info('')


! call launch_browser('\structural_report.HTML', 'internal')
 IF(my_browser%exist) then
  call launch_browser(TRIM(structural_report_file))
  call write_info('')
  call write_info('  Please wait. Your browser will be launched to display the HTML report.')
  call write_info('')
 else
  call write_info('')
  write(message_text, '(3a)') '  No browser defined in the ', trim(cryscal_ini), ' setting file !!'
  call write_info(trim(message_text))
  call write_info('')
 endif
 !call launch_word('')
 return




end subroutine create_structural_report



!---------------------------------------------------------------------
subroutine get_sym_op(input_string, numor_op, t)
 implicit none
  character (len=*), intent(in)       :: input_string
  integer              , intent(out)  :: numor_op
  integer, dimension(3), intent(out)  :: t
  integer                             :: i1

  i1 = index(input_string, '_')
  if(i1==0) then
    read(input_string, *) numor_op
    t(1) = 5
    t(2) = 5
    t(3) = 5
  else
   read(input_string(1:i1-1), *) numor_op
   read(input_string(i1+1: i1+1), *) t(1)
   read(input_string(i1+2: i1+2), *) t(2)
   read(input_string(i1+3: i1+3), *) t(3)
  endif

  t(1:3) = t(1:3) - 5

 return
end subroutine get_sym_op


!------------------------------------------------------------------------


subroutine TRANSF_moiety_string(code, input_string, output_string)
 use macros_module, only : nombre_de_colonnes, check_character
 implicit none
  character (len=*),   intent(in)     :: code             !(HTML/CIF)
  character (len=*),   intent(inout)  :: input_string
  character (len=256), intent(out)    :: output_string
  character (len=256)                 :: tmp_string
  character (len=8), dimension(2)     :: balise
  integer                             :: i, j, i1, long, nb_col, num, num_alpha
  character (len=8),  dimension(32)   :: part_string
  LOGICAL                             :: alpha_char, numeric_char



  output_string = ''
  if(code(1:4) == 'HTML') then
   balise(1) = '<sub>'
   balise(2) = '</sub>'
  elseif(code(1:3) == 'CIF') then
   balise(1) = '~'
   balise(2) = '~'
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nb_col = 0
  tmp_string = input_string
  do
   i1 = index(tmp_string, ' ')
   if(i1 ==0) exit
   if(len_trim(tmp_string) ==0) exit
   nb_col = nb_col + 1
   read(tmp_string(1:i1-1), '(a)') part_string(nb_col)
   tmp_string = tmp_string(i1+1:)
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! pb si la chaine contient le caractere virgule ','
  !input_string = adjustl(input_string)
  !call nombre_de_colonnes(input_string, nb_col)
  !if(nb_col == 0) then
  ! output_string = input_string
  ! return
  !endif
  !read(input_string, *) part_string(1:nb_col)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1, nb_col

   num       = 0
   num_alpha = 0
   do j=1, len_trim(part_string(i))
    !num       = 0
    !num_alpha = 0

    call check_character (part_string(i)(j:j), alpha_char, numeric_char)
    if(alpha_char) then
     num_alpha = num_alpha + 1
     if(num ==0) then
      if(num_alpha == 1 .and. i/=1) then
       write(OUTPUT_string, '(3a)') trim(OUTPUT_string), ' ' , part_string(i)(j:j)
      else
       write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)
      endif
     else
      if(num_alpha == 1) then
       write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)
      else
       write(OUTPUT_string, '(3a)') trim(OUTPUT_string), trim(balise(2)), part_string(i)(j:j)
      endif
      !write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)

     end if
     num = 0

    else
     num       = num+1
     if(num_alpha /=0) then
      if(num == 1 .and. j/=1) then
       write(OUTPUT_string, '(3a)') trim(OUTPUT_string), trim(balise(1)), part_string(i)(j:j)
      else
       write(OUTPUT_string, '(2a)') trim(OUTPUT_string), part_string(i)(j:j)
      endif
     else
      write(OUTPUT_string, '(3a)') trim(OUTPUT_string), ' ', part_string(i)(j:j)
      num_alpha = 0
     endif
     !num_alpha = 0
    end if
   end do

   if(num/=0) then
    write(OUTPUT_string, '(2a)') trim(OUTPUT_string), trim(balise(2))
   else
    write(OUTPUT_string, '(2a)') trim(OUTPUT_string), ' '
   endif

  end do


 return
end subroutine TRANSF_moiety_string


!---------------------------------------------------------------------

subroutine  Get_num_site_sym(n_sym, site_sym, current_sym, num_site_sym)
 implicit none
  integer,                             intent(in)    :: n_sym
  CHARACTER (LEN=12), dimension(1000), intent(in)    :: site_sym
  CHARACTER (LEN=12),                  intent(in)    :: current_sym
  INTEGER           , dimension(1000), intent(inout) :: num_site_sym
  INTEGER                                            :: i


  do i = 1, n_sym -1
   if(current_sym == site_sym(i)) then
    num_site_sym(n_sym) = num_site_sym(i)
    return
    !exit
   endif
  end do

  !do i=1, n_sym -1
  ! if(current_sym /= site_sym(i)) then
  !  exit
  ! endif
  !end do

  num_site_sym(n_sym) = MAXVAL(num_site_sym) +1


 return
end subroutine Get_num_site_sym

!---------------------------------------------------------------------

subroutine  Get_num_site_sym_new(n_sym, site_sym, current_sym, num_site_sym, dico)
 implicit none
  integer,                             intent(in)    :: n_sym
  CHARACTER (LEN=12), dimension(1000), intent(in)    :: site_sym
  CHARACTER (LEN=12),                  intent(in)    :: current_sym
  INTEGER           , dimension(1000), intent(inout) :: num_site_sym
  character (len=12), dimension(1000), intent(inout) :: dico

  INTEGER                                            :: i

  num_site_sym(1) = 0
  do i = 1, n_sym

   if(dico(i) == '') then
    dico(i) = current_sym
    num_site_sym(n_sym) = num_site_sym(i)+1
    return
   else
    if(dico(i) == current_sym) then
     num_site_sym(n_sym) = i
     return
    end if
   end if

  end do




 return
end subroutine Get_num_site_sym_new

!---------------------------------------------------------------------
subroutine build_HTML_string(n, HTML_string)

 implicit none
  integer,            intent(in)          :: n
  character(len=256), intent(inout)       :: HTML_string

  if(n < 10) then
   write(HTML_string, '(a,i1)') '#', n
  elseif(n < 100) then
   write(HTML_string, '(a,i2)') '#', n
  else
   write(HTML_string, '(a,i3)') '#', n
  endif

 return
end subroutine build_HTML_string

!---------------------------------------------------------------------
subroutine write_HTML_sym(n_sym, op_string, op_n, op_t)
 use cryscal_module, only : HTML_unit
 implicit none

  integer,                             intent(in) :: n_sym
  character (len=256), dimension(192), intent(in) :: op_string
  integer, dimension(500),             intent(in) :: op_n       ! numero de l'op. de sym. associe
  integer, dimension(500,3),           intent(in) :: op_t       ! partie translation
  integer                                         :: i


    do i=1, n_sym
     if(op_n(i) < 10) then
      write(HTML_unit, '(10x,a,i1,3a,3(i2,a))') " #",i, "  ", op_string(op_n(i))(1:24), &
                                                "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
     elseif(op_n(i) < 100) then
      write(HTML_unit, '(10x,a,i2,3a,3(i2,a))') " #",i, "  ", op_string(op_n(i))(1:24), &
                                                "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
     else
      write(HTML_unit, '(10x,a,i3,3a,3(i2,a))') " #",i, "  ", op_string(op_n(i))(1:24), &
                                                "  T = [", op_t(i,1), ", ", op_t(i,2), ", ", op_t(i,3), "]"
     endif
    end do

 return
end subroutine write_HTML_sym


