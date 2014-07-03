
module LATEX_module

 character (len=256), dimension(2)  :: logo

end module LATEX_module


subroutine Latex_preambule
 use LATEX_module,    only : logo
 use cryscalc_module, only : LATEX_unit, archive_cif, cryscalc 
 implicit none
  character (len=256), dimension(2)   :: logo_full
  logical,             dimension(2)   :: file_exist
  character (len=256)                 :: DOS_command
  
  logo_full = ''
  logo(1) = 'CDIFX_logo.jpg'
  logo(2) = 'ISCR_logo.png'
 
 WRITE(LATEX_unit, '(a)')  "% Structural report (LATEX file) created by CRYSCALC.exe"
 WRITE(LATEX_unit, '(2a)') "%      Starting CIF file : ", trim(archive_CIF)
 WRITE(LATEX_unit, '(a)')  "%"
 WRITE(LATEX_unit, '(a)')  "%     Web site : http://www.cdifx.univ-rennes/cryscalc.htm"
 WRITE(LATEX_unit, '(5a)') "%     Version : ", trim(cryscalc%version), " [", trim(cryscalc%author), "]"
 WRITE(LATEX_unit, '(a)')  "%"

 write(LATEX_unit, '(a)') '\documentclass[11pt , a4paper]{article}'
 write(LATEX_unit, '(a)') '%% ///////   PREAMBULE  ////// %%'
 write(LATEX_unit, '(a)') '\usepackage[francais,english]{babel} % adaptation de LaTeX à la langue française'
 write(LATEX_unit, '(a)') '% geometrie de la page'
 write(LATEX_unit, '(a)') '\usepackage[dvips,lmargin=2.5cm,rmargin=2.5cm,tmargin=2.5cm,bmargin=2.5cm]{geometry}'
 write(LATEX_unit, '(a)') '\usepackage{color}               % utilisation des couleurs'
 write(LATEX_unit, '(a)') "\usepackage{graphicx}            % insertion d'images"
 write(LATEX_unit, '(a)') '\usepackage{fancybox}            % fonctions encadrement'
 write(LATEX_unit, '(a)') '\renewcommand{\thefootnote}{\*}  % supprime le numero de la footnote'
 write(LATEX_unit, '(a)') "\usepackage[pdftex,colorlinks=true,urlcolor=gris50,pdfstartview=FitH]{hyperref}   % creation d'hyperliens"
 write(LATEX_unit, '(a)') ''
 write(LATEX_unit, '(a)') "\usepackage{alltt}  % permet l'ecriture mathematique a l'interieur de l'environnement verbatim"
 write(LATEX_unit, '(a)') '%%%%%%%%%% defintion de couleurs %%%%'
 write(LATEX_unit, '(a)') '\definecolor{gris50}{gray}{0.50}'
 write(LATEX_unit, '(a)') '\definecolor{gris75}{gray}{0.75}'
 write(LATEX_unit, '(a)') '\definecolor{gris_clair} {rgb}{0.96078, 0.96078, 0.96078}    % rgb(245,245,245)  #f5f5f5'
 write(LATEX_unit, '(a)') '\definecolor{violet}     {rgb}{0.5,     0,       0.5}        % rgb{128,0,128)    #800080'
 write(LATEX_unit, '(a)') '\definecolor{red_25b}    {rgb}{0.5,     0,       0.1}        % rgb{127,0,25)     #7F0019'
 write(LATEX_unit, '(a)') '\definecolor{vert_titre} {rgb}{0.95294, 0.98039, 0.90196}    % rgb(243,250,230)  #F3FAE6'
 write(LATEX_unit, '(a)') '\definecolor{fonce_titre}{rgb}{0.28125, 0.22266, 0.38672}    % rgb(72, 57,99)    #483963'
 write(LATEX_unit, '(a)') '%%%%%%%%%%%%% macros %%%%%%%%%%%%%%%%'
 write(LATEX_unit, '(a)') '\newcommand{\bs}   {$\backslash$}'
 write(LATEX_unit, '(a)') '\newcommand{\header} [1]{'
 write(LATEX_unit, '(a)') '\tiny'
 write(LATEX_unit, '(a)') '\noindent'
 write(LATEX_unit, '(a)') '\textcolor{gris75}{#1}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') ''
 write(LATEX_unit, '(a)') '% titre'
 write(LATEX_unit, '(a)') '\newcommand{\titre}  [1] {'
 write(LATEX_unit, '(a)') ' \begin{center}'
 write(LATEX_unit, '(a)') '\fcolorbox{black}{vert_titre}'
 write(LATEX_unit, '(a)') '{'
 write(LATEX_unit, '(a)') '\parbox{8cm}{'
 write(LATEX_unit, '(a)') '\begin{center}'
 write(LATEX_unit, '(a)') '\textcolor{fonce_titre}{\textit{#1}}'
 write(LATEX_unit, '(a)') '\end{center}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '\vspace{0.5cm}'
 write(LATEX_unit, '(a)') '\end{center}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '%'
 write(LATEX_unit, '(a)') '% sous-titre'
 write(LATEX_unit, '(a)') '\newcommand{\SousTitre}  [1] {'
 write(LATEX_unit, '(a)') '\vspace{0.5cm}'
 write(LATEX_unit, '(a)') '\noindent'
 write(LATEX_unit, '(a)') '\fcolorbox{gris_clair}{gris_clair}'
 write(LATEX_unit, '(a)') '{'
 write(LATEX_unit, '(a)') '\parbox{14.5cm}{'
 write(LATEX_unit, '(a)') '\textcolor{violet}{\hspace{0.0cm}#1}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '\vspace{0.5cm}'
 write(LATEX_unit, '(a)') '}'
 write(LATEX_unit, '(a)') '%'
 write(LATEX_unit, '(a)') '%% ///////   FIN DU PREAMBULE  ////// %%'
 write(LATEX_unit, '(a)') ''
 write(LATEX_unit, '(a)') '\begin{document}'
 write(LATEX_unit, '(a)') '\renewcommand*\familydefault{\sfdefault} %% Only if the base font of the document is to be sans serif'
 write(LATEX_unit, '(a)') '\normalfont'
 write(LATEX_unit, '(a)') ''
 
 if(len_trim(cryscalc%path_name) /= 0) then
  file_exist = .false.
  write(logo_full(1), '(3a)') trim(cryscalc%path_name), '\img\', trim(logo(1))
  inquire(file = trim(logo_full(1)), exist = file_exist(1))
  write(logo_full(2), '(3a)') trim(cryscalc%path_name), '\img\', trim(logo(2))
  inquire(file = trim(logo_full(2)), exist = file_exist(2))
  
  if(file_exist(1) .and. file_exist(2)) then
   ! les fichiers doivent etre dans le repertoire de travail   
   write(DOS_command, '(3a)') 'copy ', trim(logo_full(1)), ' .'     
   call system(trim(DOS_command))
   
   write(DOS_command, '(3a)') 'copy ', trim(logo_full(2)), ' .'
   call system(trim(DOS_command))   
   
  
   write(LATEX_unit, '(a)') ''
   write(LATEX_unit, '(a)') '\begin{figure}[h]'
   !write(LATEX_unit, '(a)') '% \includegraphics[width=50.pt]{img/ISCR_logo.png} \includegraphics[width=40.pt]{img/cdifx_logo.jpg}'
   write(LATEX_unit, '(5a)') ' \includegraphics[width=50.pt]{', trim(logo(1)), '}  \includegraphics[width=40.pt]{', trim(logo(2)), '}'
   write(LATEX_unit, '(a)') '\end{figure}'
   write(LATEX_unit, '(a)') ''
   
   
  endif
 endif

  return
end subroutine Latex_preambule




subroutine Latex_End
 use LATEX_module,    only : logo
 use cryscalc_module, only : LATEX_unit
 
 implicit none
 
 write(LATEX_unit, '(a)')''
 write(LATEX_unit, '(4a)') '\footnote{\textcolor{gris50}{This structural report has been created through ', &
                           '\href{http://www.cdifx.univ-rennes1.fr/CRYSCALC}{CRYSCALC} (',  trim(cryscalc%version), ').'
 write(LATEX_unit, '(a)') 'Please report bugs and problems to \href{mailto:cdifx@univ-rennes1.fr}{cdifx@univ-rennes1.fr}}}'

 write(LATEX_unit, '(a)') '\end{document}'

 !call system('del '//trim(logo(1)))
 !call system('del '//trim(logo(2)))
   
 return
end subroutine Latex_End


subroutine check_LATEX_file_name(LATEX_string)
 use macros_module, only : replace_car
 implicit none
 character (len=*), intent(inout)  :: LATEX_string

  !LATEX_string = replace_car(LATEX_string, '\',   '/bs ')
  !LATEX_string = replace_car(LATEX_string, '/bs', '\bs ')
  LATEX_string = replace_car(LATEX_string, '\', '\bs ')
  
  !LATEX_string = replace_car(LATEX_string, '_',   '\/')
  !LATEX_string = replace_car(LATEX_string, '\/',  '\_')
  LATEX_string = replace_car(LATEX_string, '_',  '\_')


 return
end subroutine check_LATEX_file_name