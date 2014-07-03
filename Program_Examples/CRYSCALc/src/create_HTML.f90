!     Last change:  TR   29 Nov 2007    6:18 pm
subroutine create_CRYSCALC_HTML()
! creation du manuel au format HTML
USE cryscalc_module,               ONLY : nb_help_max, help_string, HTML_unit, browse_cryscalc_HTML, my_browser, &
                                          cryscalc, max_ref, debug_proc
USE text_module
USE external_applications_module, ONLY : launch_browser
USE MATRIX_list_module
USE IO_module

implicit none
 INTEGER             :: i, k
 CHARACTER (LEN=512) :: HTML_text
 CHARACTER (LEN=256) :: current_directory
 CHARACTER (LEN=256) :: DOS_command

 if(debug_proc%level_2)  call write_debug_proc_level(2, "create_CRYSCALC_HTML")
 
 
 OPEN(UNIT=HTML_unit, FILE='cryscalc.HTML')
  WRITE(HTML_unit, '(a)')  "<!-- CRYSCALC user's guide (HTML file) created by CRYSCALC.exe"
  WRITE(HTML_unit, '(a)')  "     Web site : http://www.cdifx.univ-rennes/cryscalc.htm"
  WRITE(HTML_unit, '(2a)') "     Version : ", trim(cryscalc%version)
  WRITE(HTML_unit, '(2a)') "     Author  : ", trim(cryscalc%author)
  WRITE(HTML_unit, '(a)')  "!-->"

  WRITE(HTML_unit, '(a)') "<HTML>"
  WRITE(HTML_unit, '(a)') "<HEAD>"

  call write_HTML_css(HTML_unit, 'cryscalc')

  WRITE(HTML_unit, '(a)') '<A NAME="cryscalc_main"></A>'
  WRITE(HTML_unit, '(a)') "<TITLE>CRYSCALC user's guide</TITLE>"
  WRITE(HTML_unit, '(a)') "</HEAD>"
  WRITE(HTML_unit, '(a)') '<BODY BGCOLOR="#FFFFFF">'
  WRITE(HTML_unit, '(a)') "<br><br>"
  WRITE(HTML_unit, '(a)') "<p class='title_main'>CRYSCALC</p>"
  !WRITE(HTML_unit, '(a)') "<center><h1><FONT COLOR="#FF0000">CRYSCALC</FONT></h1></center>"
  WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"


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
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#kw_list">&nbsp;&nbsp;List of CRYSCALC keywords</A>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#kw_details">&nbsp;&nbsp;Details of CRYSCALC keywords</A>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cla_section">&nbsp;&nbsp;List of CRYSCALC command line arguments</A>'
  WRITE(HTML_unit, '(a)') "<li class='item_2'><A HREF='#cryscalc_news'>&nbsp;&nbsp;What's new ?</A>"
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cryscalc_ini">&nbsp;&nbsp;CRYSCALC.ini setting file</A>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cryscalc_cfl">&nbsp;&nbsp;Examples of .CFL input files</A>'

  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cryscalc_links">&nbsp;&nbsp;CRYSCALC download and links</a>'
  WRITE(HTML_unit, '(a)') '<li class="item_2"><A HREF="#cryscalc_linux">&nbsp;&nbsp;Linux version</a>'
  WRITE(HTML_unit, '(a)') '</ul>'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '</p>'


  WRITE(HTML_unit, '(a)') '<A NAME="intro"></A><p class="title_1">&nbsp;&nbsp;Introduction :<br></p>'

  ! HEADER
  call def_header_lines()
  do i=1, header_lines_nb
   WRITE(HTML_unit, '(a)') TRIM(header_line(i))
  end do
  WRITE(HTML_unit, '(a)') ""
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscalc_main">[return]</a>'
  WRITE(HTML_unit, '(a)') "</pre>"

  WRITE(HTML_unit, '(a)') '<A NAME="kw_list"></A><p class="title_1">&nbsp;&nbsp;List of CRYSCALC keywords :<br></p>'
  WRITE(HTML_unit, '(a)') '<ul><ul>'
  do i=1, nb_help_max
   WRITE(HTML_text, '(5a)') '   <li class="item"><A HREF="#', TRIM(help_string(i)), '">', TRIM(help_string(i)), '</A>'
   WRITE(HTML_unit, '(a)') TRIM(HTML_text)
  end do
  WRITE(HTML_unit, '(a)') '</ul></ul><br>'
  WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'> "


  call def_keywords_lines()
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscalc_main">[return]</a>'

  WRITE(HTML_unit, '(a)') '<A NAME="kw_details"></A><p class="title_1">&nbsp;&nbsp;Details of CRYSCALC keywords :<br></p>'
  do i=1, nb_help_max
   WRITE(HTML_unit, '(3a)') '<A NAME="', TRIM(help_string(i)), '"></A>'
   WRITE(HTML_unit, '(3a)') "<p class='title_3'>", trim(HELP_line(i,2)), "</p> <pre style='font-size:14' class='color_grey1'>"
   do k=4, HELP_lines_nb(i)
    WRITE(HTML_unit, '(a)') TRIM(HELP_line(i,k))
   end do
   WRITE(HTML_unit, '(a)')'<A HREF="#kw_list">[return to keywords list]</a>'
   WRITE(HTML_unit, '(a)')'<hr></pre>'
  end do



  call Def_command_line_arguments
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscalc_main">[return]</a>'
  WRITE(HTML_unit, '(a)') '<A NAME="cla_section"></A><p class="title_1"> List of CRYSCALC command line arguments :<br></p>'
  !WRITE(HTML_unit, '(a)') "<pre style='font-size:14' class='color_grey1'>"
  !do i=1, CLA_lines_nb
  ! WRITE(HTML_unit, '(a)') TRIM(cla_line(i))
  !end do
  !WRITE(HTML_unit, '(a)') ""
  !WRITE(HTML_unit, '(a)') "</pre>"

  !WRITE(HTML_unit, '(a)') '<ul>'
  do i=1, CLA_nb
   WRITE(HTML_unit, '(3a)') "<p class='title_3'>", trim(CLA_line(i,1)), "</p><pre style='font-size:14' class='color_grey1'>"
   do k=2, CLA_lines_nb(i)
    WRITE(HTML_unit, '(a)') TRIM(CLA_line(i,k))
   end do
   WRITE(HTML_unit, '(a)') '</pre>'
   !WRITE(HTML_text, '(2a)') '   <li class="item">', TRIM(CLA_line(i))
   !WRITE(HTML_unit, '(a)') TRIM(HTML_text)
  end do
  !WRITE(HTML_unit, '(a)') '</ul><br>'
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscalc_main">[return]</a>'
  WRITE(HTML_unit, '(a)') "</pre>"


! --------------- what's new in CRYSCALC ? ------------------------

  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') "<A NAME='cryscalc_news'></A><p class='title_1'>&nbsp;&nbsp;What's new in CRYSCALC ?<br></p>"
  WRITE(HTML_unit, '(a)') '<pre>'
  call WRITE_CRYSCALC_NEWS("HTML")
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)')'<A HREF="#cryscalc_main">[return]</a>'
  WRITE(HTML_unit, '(a)') '</pre>'


! --------------- setting file -----------------------------------
  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') '<A NAME="cryscalc_ini"></A><p class="title_1">&nbsp;&nbsp;CRYSCALC.ini setting file :<br></p>'
  WRITE(HTML_unit, '(a)') '<p class="retrait_1">A setting file can be used by <b>CRYSCALC</b>, containing '
  WRITE(HTML_unit, '(a)') 'the definition of different parameters such as external applications that can be'
  WRITE(HTML_unit, '(a)') 'launched from <b>CRYSCALC</b> (browser, editor ...) or defaut values about'
  WRITE(HTML_unit, '(a)') 'diffractometer, author, structure solution and refinement programs  ... '
  WRITE(HTML_unit, '(a)') 'This setting file, called <font face="courier">cryscalc.ini</font> has to be located'
  WRITE(HTML_unit, '(a)') 'in the folder related to <b>CRYSCALC</b> through the <font face="courier">CRYSCALC</font>'
  WRITE(HTML_unit, '(a)') 'environment variable.</p>'
  WRITE(HTML_unit, '(a)') '<p class="retrait_1">Example of setting file:</p>'
  WRITE(HTML_unit, '(a)') '<pre>'
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[EXTERNAL APPLICATIONS]</font>'
  WRITE(HTML_unit, '(a)') '    browser = "C:\Program Files\Mozilla Firefox\firefox.exe"'
  WRITE(HTML_unit, '(a)') '    editor  = "C:\Program Files\Keditw\KEDITW32.EXE"'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[WEB ADDRESS]</font>'
  WRITE(HTML_unit, '(a)') '    fps        = www.ill.fr/dif/Soft/fp/'
  WRITE(HTML_unit, '(a)') '    cdifx      = www.cdifx.univ-rennes1.fr/'
  WRITE(HTML_unit, '(a)') '    cryscalc   = www.cdifx.univ-rennes1.fr/cryscalc'
  WRITE(HTML_unit, '(a)') '    reciprocs  = www.cdifx.univ-rennes1.fr/reciprocs'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[DEVICE]</font>'
  WRITE(HTML_unit, '(a)') '    diffractometer = APEXII AXS Bruker'
  WRITE(HTML_unit, '(a)') '    !diffractometer = KCCD Nonius'
  WRITE(HTML_unit, '(a)') '    !diffractometer = Bruker SMART X2S benchtop'
  WRITE(HTML_unit, '(a)') '    laboratory     = CDIFX Rennes'
  WRITE(HTML_unit, '(a)') '    radiation      = X_Mo'
  WRITE(HTML_unit, '(a)') '    wave_A         = 0.71073'
  WRITE(HTML_unit, '(a)') '    temperature    = 150(2)'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[USER]</font>'
  WRITE(HTML_unit, '(a)') '    name        = ROISNEL'
  WRITE(HTML_unit, '(a)') '    first_name  = Thierry'
  WRITE(HTML_unit, '(2a)')'    address     = Centre de Diffractométrie X, UMR6511 CNRS Université de Rennes 1,', &
                                            'Sciences Chimiques de Rennes, 35042 RENNES Cedex France'
  WRITE(HTML_unit, '(a)') '    email       = thierry.roisnel@univ-rennes1.fr'
  WRITE(HTML_unit, '(a)') '    web         = www.cdifx.univ-rennes1.fr'
  WRITE(HTML_unit, '(a)') '    team        = CDIFX'
  

  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[ARRAYS DIMENSIONS]</font>'
  WRITE(HTML_unit, '(a)') '    hkl_reflections = 200000       ! max. number of hkl reflections in a file'

  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[CREATE INS]</font>'
  WRITE(HTML_unit, '(a)') '    get_sample_ID = 0          ! get sample ID (default=job)'
  WRITE(HTML_unit, '(a)') '    temperature   = 100K       ! experimental temperature value'
  WRITE(HTML_unit, '(a)') '    u_threshold   = 0.1        ! atoms with U_iso > U_threshold will be excluded'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[PARAMETERS]</font>'
  WRITE(HTML_unit, '(a)') '    i_sig       = 3.         ! used in the SEARCH_GROUP procedure'
  WRITE(HTML_unit, '(a)') '    threshold   = 0.03       ! used in the SEARCH_GROUP procedure '
  WRITE(HTML_unit, '(a)') '    d_max_A     = 3.5        ! used with the CONN keyword (connectivity calculation)'
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[COMMAND LINE ARGUMENTS]</font>'
  WRITE(HTML_unit, '(a)') '    create_ACE       = 1          ! .ACE file for Carine'
  WRITE(HTML_unit, '(a)') '    create_CEL       = 1          ! .CEL file for PowderCELL'
  WRITE(HTML_unit, '(a)') '    create_CFL       = 1          ! .CFL file for CRYSCALC'
  WRITE(HTML_unit, '(a)') '    create_INS       = 1          ! .INS file for SHELXL'
  WRITE(HTML_unit, '(a)') '    create_FST       = 1          ! .FST file for FP Studio'
  WRITE(HTML_unit, '(a)') '    create_CIF_pymol = 0          ! X_pml.CIF for PYMOL'
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
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[OPTIONS]</font>'
  WRITE(HTML_unit, '(a)') '    LOCK_wave_value              = 0.02      ! lock current wavelength to anticathode value'
  WRITE(HTML_unit, '(2a)')'    CIF_format80                 = 0         ! formatted line, when creating a CIF file, ', &
                                                                        ' if more than 80 characters'
  WRITE(HTML_unit, '(2a)')'    CIF_torsion_limit            = 170.      ! exclude torsion angle if greater than this limit'
  WRITE(HTML_unit, '(a)') '    include_RES_file             = 1         ! include .RES file in the archive_cryscalc.cif file'
  WRITE(HTML_unit, '(2a)')'    update_parameters            = 1         ! update parameters after transformation', &
                                                                        '(cell parameters, atomic coordinates) '
  WRITE(HTML_unit, '(2a)')'    report_header                = 1         ! Write header in structural report'
  WRITE(HTML_unit, '(2a)')'    skip_start_menu              = 0         ! Skip start menu'
  WRITE(HTML_unit, '(2a)')'    hkl_statistics               = 1         ! Ouput statistics on hkl reflections'
  WRITE(HTML_unit, '(2a)')'    hkl_format                   = 3I4,2F8.2 ! format for .hkl file (h,k,l,F2,sig)' 
  WRITE(HTML_unit, '(2a)')'    cartesian_frame_type         = A         ! A: x//a ; C: x //c' 
  
  WRITE(HTML_unit, '(a)') ''
  WRITE(HTML_unit, '(a)') '    <font color="darkred">[PATTERN SIMULATION (Pseudo-Voigt profile)]  </font>'
  
  WRITE(HTML_unit, '(2a)') '    X_profile_U                  = 0.0055    ! U value of the Cagliotti formula : ', &
                                                                        'FWHM2 = U*TAN**2(theta) + V*TAN(theta) + W'
  WRITE(HTML_unit, '(a)')  '    X_profile_V                  = -0.0015   ! V value '
  WRITE(HTML_unit, '(a)')  '    X_profile_W                  = 0.0036    ! W value '
  WRITE(HTML_unit, '(a)')  '    X_profile_eta0               = 0.3       ! Lorentzian components : eta = eta° + 2theta * eta1'
  WRITE(HTML_unit, '(a)')  '    X_profile_eta1               = 0.         '
  WRITE(HTML_unit, '(a)')  '    X_pattern_step               = 0.01       '
  WRITE(HTML_unit, '(a)')  '    X_pattern_scale              = 1000.      '
  WRITE(HTML_unit, '(a)')  '    X_pattern_background         = 50.        '
  WRITE(HTML_unit, '(a)')  '    N_profile_U                  = 0.0146     '
  WRITE(HTML_unit, '(a)')  '    N_profile_V                  = -0.0375    '
  WRITE(HTML_unit, '(a)')  '    N_profile_W                  = 0.0475     '
  WRITE(HTML_unit, '(a)')  '    N_profile_eta0               = 0.01       '
  WRITE(HTML_unit, '(a)')  '    N_profile_eta1               = 0.         '
  WRITE(HTML_unit, '(a)')  '    N_pattern_step               = 0.025      '
  WRITE(HTML_unit, '(a)')  '    N_pattern_scale              = 1000.      '
  WRITE(HTML_unit, '(a)')  '    N_pattern_background         = 50.        '

  
  WRITE(HTML_unit, '(a)') ''
  if(user_mat_nb /=0) then
    WRITE(HTML_unit, '(a)') '    <font color="darkred">[USER TRANSFORMATION MATRICES]</font>'    
	do i=1, user_mat_nb
	 WRITE(HTML_unit, '(a,i1,a,3(2x,3F4.0),2a)') '    MAT_',i,'                        = ', transf_mat(:,:, max_mat_nb+i), &
	                                             '   !  ', trim(user_mat_text(i))
	end do 
    WRITE(HTML_unit, '(a)') ''	 
  endif
  WRITE(HTML_unit, '(a)') '<A HREF="#cryscalc_main">[return]</a>'
  WRITE(HTML_unit, '(a)') '</pre>'

  
! ----------------------------------------------------------------
  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') '<A NAME="cryscalc_cfl"></A><p class="title_1">&nbsp;&nbsp;Examples of .CFL input files :<br></p>'
  WRITE(HTML_unit, '(a)') '<ul>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/rint.cfl">', &
                           '  Calculation of internal R factor from hkl data included in a import.cif file</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_pnma.cfl">', &
                           '  Pnma space group information</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_pbnm_pnma.cfl">', &
                           '  Transformation from Pbnm to Pnma space group</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_si_x.cfl">', &
                           '  Simulation of a X-ray diffraction pattern for Si</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_si_x_100.cfl">', &
                           '  Simulation of a X-ray diffraction pattern for Si (particles size=100A)</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_si_n.cfl">', &
                           '  Simulation of a neutron diffraction pattern for Si</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_y2o3.cfl">', &
                           '  Atom connectivity in Y2O3</a><p>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_y2o3_bvs.cfl">', &
                           '  Bond valence calculation in Y2O3</a>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_ambi_mu.cfl">', &
                           '  X-ray absorption coefficient calculation (case of ammonium bitartrate)</a>'
  WRITE(HTML_unit, '(a)') '<p>'
  WRITE(HTML_unit, '(a)') '</ul>'
 
  
! ----------------------------------------------------------------
  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') '<A NAME="cryscalc_links"></A><p class="title_1">&nbsp;&nbsp;CRYSCALC links :<br></p>'
  WRITE(HTML_unit, '(a)') '<ul>'
  WRITE(HTML_unit, '(2a)') ' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc.exe">', &
                          'CRYSCAL.exe (for Windows)</a>'
  WRITE(HTML_unit, '(a)') '<p>'

  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc.ini">', &
                          'Example of CRYSCAL setting file</a>'
  WRITE(HTML_unit, '(a)') '<p>'
  
  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc.css">', &
                          "Example of CSS file for HTML user's guide</a>"
  WRITE(HTML_unit, '(a)') '<p>'
  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc_report.css">', &
                          "Example of CSS file for HTML structural report</a>"
  WRITE(HTML_unit, '(a)') '<p>'
  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="https://forge.epn-campus.eu/projects/crysfml">', &
                          "CRYSFML repository</a>"
  
  WRITE(HTML_unit, '(a)') '<p>'
  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.iucr.org/resources/commissions/crystallographic-computing/', &
                          'newsletters/1/crysfml" target="_blank">CrysFML:</A>'
  WRITE(HTML_unit, '(2a)')' <font color="black">Crystallographic Fortran Modules Library by J. Rodriguez-Carvajal and J. ',   &
                          'González-Platas</font>'
  WRITE(HTML_unit, '(a)') '</ul>'
  WRITE(HTML_unit, '(a)') ''
  
  ! ----------------------------------------------------------------
  WRITE(HTML_unit, '(a)') '<br>'
  WRITE(HTML_unit, '(a)') '<A NAME="cryscalc_linux"></A><p class="title_1">&nbsp;&nbsp;CRYSCALC under Linux:<br></p>'
  WRITE(HTML_unit, '(a)') '<ul>'
  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/linux/readme.txt">', &
                          ' Readme.txt</a>'
  WRITE(HTML_unit, '(a)') '<p>'
  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/linux/cleaned.sh">', &
                          ' Script for compilation</a>'
  WRITE(HTML_unit, '(a)') '<p>'
  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/linux/cryscalc_init_linux.', &
                          'f90">cryscalc_init_linux.f90</a>'
  WRITE(HTML_unit, '(a)') '<p>'
  WRITE(HTML_unit, '(2a)')' <li class="item"><a href="http://www.cdifx.univ-rennes1.fr/progs/cryscalc/linux/cryscalc_linux">', &
                          'cryscalc_linux (executable)</a>'
  WRITE(HTML_unit, '(a)') '<p>'

  WRITE(HTML_unit, '(a)') '</ul>'
  WRITE(HTML_unit, '(a)') ''

  

  !WRITE(HTML_unit, '(a)') '<p class="retrait_1"><b>Author:</b><br>'
  !WRITE(HTML_unit, '(a)') 'Thierry Roisnel <A HREF="mailto:cdifx@univ-rennes1.fr">cdifx.univ-rennes1.fr</A>'
  WRITE(HTML_unit, '(a)')  "<br><br><hr>"
  WRITE(HTML_unit, '(3a)') "<p class='retrait_2'>last updated: <i>TR / CDIFX - ", trim(cryscalc%version), "</i></P>"
  !WRITE(HTML_unit, '(a)') 'Web site: <a href="www.cdifx.univ-rennes1.fr/CRYSCALC" target="_blank"> www.cdifx.univ-rennes1.fr/CRYSCALC</a>'

  WRITE(HTML_unit, '(a)') "</BODY>"
  WRITE(HTML_unit, '(a)') "</HTML>"
CLOSE(UNIT=HTML_unit)


  call write_info("")
  call write_info("   The CRYSCALC user's guide (cryscalc.HTML) in HTML format has been created." )
  call write_info("")

  IF(browse_cryscalc_HTML .and. my_browser%exist) then
   call launch_browser('cryscalc.HTML')
  ENDIF

 return
end subroutine create_CRYSCALC_HTML

!--------------------------------------------------------------------------------

subroutine write_HTML_css(HTML_unit, input_string)
 use cryscalc_module, only : tmp_unit, cryscalc, debug_proc
 implicit none
  integer, intent(in)           :: HTML_unit
  character (len=*), intent(in) :: input_string  !  input_string = cryscalc / report
  character (len=256)           :: css_file
  character (len=1024)          :: read_line, adjusted_line
  logical                       :: file_exist
  integer                       :: i_error
  
  if(debug_proc%level_2)  call write_debug_proc_level(2, "write_HTML_css ("//trim(input_string)//")")

  
  if(input_string == 'cryscalc') then
   css_file = cryscalc%css
  elseif(input_string == 'report') then
   css_file = cryscalc%report_css
  endif
 
  inquire(file = trim(css_file), exist=file_exist)

  if(.not. file_exist) then  
  WRITE(HTML_unit, '(a)') "<style TYPE='text/css'>"
  WRITE(HTML_unit, '(2a)') "  .title_main { font-family:'Trebuchet MS',Arial; font-size:24px; color:#AA0000;", &
                           " margin-left:50px; margin-right:50px; text-align:center; font-weight:bold; }"
  WRITE(HTML_unit, '(2a)') "  .title_1    { font-family:'Trebuchet MS',Arial; font-size:20px; color:#880000;",  &
                           " margin-left:20px; background-color:#f5f5f5; }"
  WRITE(HTML_unit, '(2a)') "  .title_2    { font-family:'Trebuchet MS',Arial; font-size:20px; color:#880000;",  &
                           " margin-left:50px; margin-right:50px; background-color:#f5f5f5; }"
  WRITE(HTML_unit, '(2a)') "  .title_3    { font-family:'Courier new';        font-size:14px; color:#660000;", &
                           " font-weight:bold; margin-left:30px;}"
  WRITE(HTML_unit, '(2a)') "  .header   { font-family:'Trebuchet MS',Arial; font-size:10px; color:#CCCCCC;",  &
                           " margin-left:50px; line-height: 10pt;}"
  WRITE(HTML_unit, '(2a)') "  .retrait_1  { font-family:'Trebuchet MS',Arial; font-size:14px; margin-left:50px;", &
                           " margin-right:50px; text-align:justify; line-height:18.0pt;  }"
  WRITE(HTML_unit, '(2a)') "  .retrait_2  { font-family:'Trebuchet MS',Arial; font-size:14px; margin-left:50px;",  &
                           " margin-right:50px; text-align:justify; line-height:18.0pt; color:#AAAAAA; }"
  WRITE(HTML_unit, '(a)') "  .item       { font-family:'Trebuchet MS',Arial; font-size:14px; color:#aa0000;}"
  WRITE(HTML_unit, '(a)') "  .item_2     { font-family:'Trebuchet MS',Arial; font-size:16px; color:#aa0000;}"
  WRITE(HTML_unit, '(a)') "  .ligne_1    { background-color : #FAFAFA;}"
  WRITE(HTML_unit, '(a)') "  .ligne_2    { background-color : #EEEEEE;}"
  WRITE(HTML_unit, '(a)') "  .color_grey1{color:#323232;}"
  WRITE(HTML_unit, '(2a)') "  hr {border-color: #CCC; background-color: #CCC; height:1px; padding: 0; border:0;", &
                           " margin-left:50px; margin-right:50px;}"
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
  WRITE(HTML_unit, '(a)') "      font-weight:bold;"
  WRITE(HTML_unit, '(a)') "      background-color:transparent;}"
  WRITE(HTML_unit, '(a)') "</style>"
  else
   open(unit = tmp_unit, file=trim(css_file))
    do
	 read(tmp_unit, '(a)', iostat=i_error) read_line 
	 if(i_error < 0) exit
	 if(len_trim(read_line) == 0) cycle
	 adjusted_line = adjustl(read_line)
	 if(adjusted_line(1:1) == '/' .or. adjusted_line(1:1) == "*") cycle
	 write(HTML_unit, '(a)') trim(read_line)
	end do
   close(unit=tmp_unit)
  endif
 

  return
end subroutine write_HTML_css
!--------------------------------------------------------------------------------


