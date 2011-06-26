!     Last change:  TR   17 Jul 2007    4:52 pm
subroutine write_cryscal_NEWS(input_string)
 use cryscal_module, only : news_year
 USE IO_module, ONLY : write_info

  implicit none
   character (len=*), intent(in) :: input_string

  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '  Main new features implemented in CRYSCAL:')
  call write_cryscal_news_line(input_string,  '  (for more details, see the CRYSCAL manual (MAN keyword))')
  call write_cryscal_news_line(input_string,  ' ')
  call write_cryscal_news_line(input_string,  '  T. Roisnel / CDIFX Rennes')
  call write_cryscal_news_line(input_string,  '')


  if(news_year(1:3) == 'all' .or. news_year(1:2) == '11' .or. news_year(1:4) == '2011') then
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . June 11 :')
  call write_cryscal_news_line(input_string,  '     # CRYSCAL has been compiled with new version of CFML (5.00)')

  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . May 11 :')
  call write_cryscal_news_line(input_string,  '     # Bug has been corrected in the structure factor calculation routine')
  call write_cryscal_news_line(input_string,  '     # Total number of electrons is output after CHEM keyword.')
  call write_cryscal_news_line(input_string,  '     # The matrix used for the hexagonal to rhomboedral system was corrupted.')

  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . April 11 :')
  call write_cryscal_news_line(input_string,  '     # Some changes in the routine to get new space group after matrix')
  call write_cryscal_news_line(input_string,  '       transformation')
  call write_cryscal_news_line(input_string,  '     # Some changes in the monoclinic transformation matrix list')
  call write_cryscal_news_line(input_string,  '     # A non systematic bug in the bonds distribution list (CONN keyword)')
  call write_cryscal_news_line(input_string,  '       has been corrected')


  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Feb. 11 :')
  call write_cryscal_news_line(input_string,  '     # Stucture factors calculation can be performed for electrons')
  call write_cryscal_news_line(input_string,  '       diffraction (GEN_HKL and SF_HKL routines). This has to be')
  call write_cryscal_news_line(input_string,  '       specified by the BEAM keyword and "ELECTRONS" argument')
  call write_cryscal_news_line(input_string,  '     # Some examples of CFL input files can be downloaded from the')
  call write_cryscal_news_line(input_string,  '       CRYSCAL web site (www.cdifx.univ-rennes1.fr/cryscal)')
  call write_cryscal_news_line(input_string,  '     # Some changes in output files :')
  call write_cryscal_news_line(input_string,  '        . cryscal.log is renamed cryscal_debug.txt')
  call write_cryscal_news_line(input_string,  '        . cryscal.out is renamed cryscal.log')
  call write_cryscal_news_line(input_string,  '       LOG argument has been replaced by DEBUG argument')
  call write_cryscal_news_line(input_string,  '     # New PAT argument for GEN_HKL keyword allows to generate a')
  call write_cryscal_news_line(input_string,  '       diffraction pattern. PRF file is automatically plotted')
  call write_cryscal_news_line(input_string,  '       with the WinPLOTR program if installed.')
  call write_cryscal_news_line(input_string,  '     # I/Imax is now calculated in the GEN_HKL routine when working')
  call write_cryscal_news_line(input_string,  '       in the 2theta space and neutron or Cu_K_alpha1 X-ray radiation.')
  call write_cryscal_news_line(input_string,  '       Intensity is calculated as follows:')
  call write_cryscal_news_line(input_string,  '                        I=mult * Lp * F^2')
  call write_cryscal_news_line(input_string,  '       where : mult is the multicicity of the reflexion')
  call write_cryscal_news_line(input_string,  '               Lp the Lorentz-polarization factor, calculated by:')
  call write_cryscal_news_line(input_string,  '               Lp=(1-K+K*CTHM*cos2(2theta)/(2sin2(theta)cos(theta)')
  call write_cryscal_news_line(input_string,  '                CTHM=cos2(2theta_monok) [CTHM=0.79]')
  call write_cryscal_news_line(input_string,  '                K=0.  for neutrons')
  call write_cryscal_news_line(input_string,  '                K=0.5 for unpolarized X-ray radiation')

  call write_cryscal_news_line(input_string,  '     # Cosmetic changes in the Fortran codes to allow the compilation')
  call write_cryscal_news_line(input_string,  '       with the free G95 Fortran compiler.')

  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Jan. 11 :')
  call write_cryscal_news_line(input_string,  '     # CREATE_FST keyword allows to create a .FST file for')
  call write_cryscal_news_line(input_string,  '       FullProf Studio after reading a CIF file.')
  call write_cryscal_news_line(input_string,  '     # "create_FST" field can be input in the "[COMMAND LINE ARGUMENTS]" part ')
  call write_cryscal_news_line(input_string,  '       in the setting file (cryscal.ini)')

  endif

  if(news_year(1:3) == 'all' .or. news_year(1:2) == '10' .or. news_year(1:4) == '2010') then
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Dec. 10 :')
  call write_cryscal_news_line(input_string,  '     # New extinction rules for FIND_HKL_LIST keyword:')
  call write_cryscal_news_line(input_string,  '       . hhl with h+l=2n')
  call write_cryscal_news_line(input_string,  '       . hkk with k+h=2n')
  call write_cryscal_news_line(input_string,  '       . hkh with h+k=2n')
  call write_cryscal_news_line(input_string,  '     # Cosmetic changes in archive_cryscal.cif file: in the case')
  call write_cryscal_news_line(input_string,  '       of samples without H atoms, "_atom_sites_solution_hydrogens" and ')
  call write_cryscal_news_line(input_string,  '       "_refine_ls_hydrogen_tretament" cif field lines are removed')
  call write_cryscal_news_line(input_string,  '')

  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Nov. 10 :')
  call write_cryscal_news_line(input_string,  '     # Connectivity calculations are followed by a bonds distribution list')
  call write_cryscal_news_line(input_string,  '     # Transformation matrix components can be input as fractional values')
  call write_cryscal_news_line(input_string,  '        Example : matrix 1/2 1/2 0.   -1/2 1/2 0. 0. 0. 1')
  call write_cryscal_news_line(input_string,  '     # Cosmetic changes in HTML structural report')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Oct. 10 :')
  call write_cryscal_news_line(input_string,  '     # A threshold value can be given as argument to the SEARCH_GROUP keyword')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Sept. 10 :')
  call write_cryscal_news_line(input_string,  '     # FIND_HKL_list argument value can be negative and allows to search')
  call write_cryscal_news_line(input_string,  '       reflections with opposite rule than for positive value.')
  call write_cryscal_news_line(input_string,  '     # SIR_TO_INS command line argument has been replaced by SOLVE_TO_INS')
  call write_cryscal_news_line(input_string,  '       and allows now to create INS file for SHELXL from SHELXS ouput file')
  call write_cryscal_news_line(input_string,  '       as well as SIRxx output file.)')
  call write_cryscal_news_line(input_string,  '     # CREATE_ACE keyword allows to create a .ACE file for CaRIne')
  call write_cryscal_news_line(input_string,  '       after reading a CIF file.')
  call write_cryscal_news_line(input_string,  '     # "create_CEL" field has been added in the "[COMMAND LINE ARGUMENTS]" part ')
  call write_cryscal_news_line(input_string,  '       in the setting file (cryscal.ini)')
  call write_cryscal_news_line(input_string,  '     # NIGGLI/NIGGLI_CELL keyword has been added and allows to determine')
  call write_cryscal_news_line(input_string,  '       the Niggli cell from any triclinic cell')
  call write_cryscal_news_line(input_string,  '     # [CREATE INS] section has been added in the crysal.ini setting file to')
  call write_cryscal_news_line(input_string,  '       define temperature and thermal parameter threshold to skip atoms')
  call write_cryscal_news_line(input_string,  '       This avoids to enter these values related to the CREATE_INS keyword')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . June 10 :')
  call write_cryscal_news_line(input_string,  '     # Bugs in the FIND_HKL_LIST routine has been corrected')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . April 10 :')
  call write_cryscal_news_line(input_string,  '     # Some items of the experimental part in the HTML structure report are now' &
                                              //' in italic,')
  call write_cryscal_news_line(input_string,  '       in agreement with published articles.')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . March 10 :')
  call write_cryscal_news_line(input_string,  '     # WRITE_DEVICE keyword')
  call write_cryscal_news_line(input_string,  '     # WRITE_HKL keyword has been replaced by FIND_HKL_LIST.')
  call write_cryscal_news_line(input_string,  '       The SUPPRESS/REMOVE argument for FIND_HKL_LIST keyword')
  call write_cryscal_news_line(input_string,  '       has been added, leading to the creation of a new file free of ')
  call write_cryscal_news_line(input_string,  '       data obeying the current selection rule.')
  call write_cryscal_news_line(input_string,  '     # QVEC components can be input as fractionnal values.')
  call write_cryscal_news_line(input_string,  '       The following fractionnal absolute values are : ')
  call write_cryscal_news_line(input_string,  '       1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 2/5, 3/5, 4/5, 1/6, 5/6,')
  call write_cryscal_news_line(input_string,  '       1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1/8, 3/8, 5/8, 7/8,')
  call write_cryscal_news_line(input_string,  '       1/9, 2/9, 4/9, 5/9, 7/9 and 8/9')
  call write_cryscal_news_line(input_string,  '     # GEN_HKL keyword leads to a structure factor calculation if atoms')
  call write_cryscal_news_line(input_string,  '       has been input.')
  call write_cryscal_news_line(input_string,  '     # SF_HKL keyword leads to a structure factor calculation for')
  call write_cryscal_news_line(input_string,  '       a given hkl reflexion')
  call write_cryscal_news_line(input_string,  '     # New WRITE_BEAM and WRITE_QVEC keywords')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Feb. 10:')
  call write_cryscal_news_line(input_string,  '     # CREATE_CEL keyword allows to create a .CEL file for PowderCELL')
  call write_cryscal_news_line(input_string,  '       after reading a CIF file.')
  call write_cryscal_news_line(input_string,  '     # CREATE_INS keyword allows to create a .INS file for SHELXL')
  call write_cryscal_news_line(input_string,  '       after reading a CIF file.')
  call write_cryscal_news_line(input_string,  '     # CREATE_CFL keyword allows to create a .CFL file for CRYSCAL')
  call write_cryscal_news_line(input_string,  '       after reading a CIF file.')
  call write_cryscal_news_line(input_string,  '     # New "[COMMAND LINE ARGUMENTS]" part in the setting file (cryscal.ini)')
  call write_cryscal_news_line(input_string,  '       with the following fields :')
  call write_cryscal_news_line(input_string,  '        . create_CEL')
  call write_cryscal_news_line(input_string,  '        . create_INS')
  call write_cryscal_news_line(input_string,  '        . create_CFL')
  call write_cryscal_news_line(input_string,  '       Putting the corresponding values to 1 will create .CEL, .INS and .CFL')
  call write_cryscal_news_line(input_string,  '       respectively.')
  call write_cryscal_news_line(input_string,  '     # cryscal file.p4p: if crystal faces are present in the .P4P file,')
  call write_cryscal_news_line(input_string,  '        they are extracted and saved in the import.cif.')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Jan. 10:')
  call write_cryscal_news_line(input_string,  '     # Bug in _trans.HKL file has been corrected.')
  call write_cryscal_news_line(input_string,  '')
  endif

  if(news_year(1:3) == 'all' .or. news_year(1:2) == '09' .or. news_year(1:4) == '2009') then
  call write_cryscal_news_line(input_string,  '   . Nov. 09:')
  call write_cryscal_news_line(input_string,  '     # If a "platon_ortep.gif" file is present in the current folder, it is')
  call write_cryscal_news_line(input_string,  '       automatically incorporated in the HTML report file')
  call write_cryscal_news_line(input_string,  '     # If final.y file contains QVEC field and related parameters, ')
  call write_cryscal_news_line(input_string,  '       a hklm file is created')

  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Oct. 09:')
  call write_cryscal_news_line(input_string,  '     # If a hkl file is input, the transformed hkl file is')
  call write_cryscal_news_line(input_string,  '       loaded automatically after MATRIX keyword')
  call write_cryscal_news_line(input_string,  '     # Shannon table value for magnesium is corrected (Mg+2)')
  call write_cryscal_news_line(input_string,  "     # Cosmetic changes in the HTML documents (report and user's guide) created")
  call write_cryscal_news_line(input_string,  '       by CRYSCAL.')
  call write_cryscal_news_line(input_string,  '     # Space group number has been added in the HTML report')
  call write_cryscal_news_line(input_string,  '     # SHELX reference has been changed to : ')
  call write_cryscal_news_line(input_string,  '       G. M. Sheldrick, Acta Cryst A, 2008, A64, 112-122')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Sept. 09:')
  call write_cryscal_news_line(input_string,  '     # TWIN_PSEUDO_HEXA keyword')
  call write_cryscal_news_line(input_string,  '     # TWIN_HEXA keyword')
  call write_cryscal_news_line(input_string,  '     # Bugs in the routine to create the archive_cryscal.cif file has ')
  call write_cryscal_news_line(input_string,  '       been corrected, specially when the archive.cif file is created')
  call write_cryscal_news_line(input_string,  '       from the "create CIF / ACTA-C" procedure in WinGX, where the ')
  call write_cryscal_news_line(input_string,  '       order of some CIF fields is changed.')
  call write_cryscal_news_line(input_string,  '     # CRYSCAL has been compiled with new version of CFML (4.00)')
  call write_cryscal_news_line(input_string,  '     # Bugs in the connectivity calculation routine has been corrected')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Feb. 09:')
  call write_cryscal_news_line(input_string,  '     # CHECK_GROUP can now determine trigonal space group')
  call write_cryscal_news_line(input_string,  '     # SG_ALL keyword output sub-groups of the current space group')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Jan. 09:')
  call write_cryscal_news_line(input_string,  '     # .RAW file can be given as argument in the command line')
  call write_cryscal_news_line(input_string,  '       i.e. d:\> CRYSCAL my_file.RAW')
  call write_cryscal_news_line(input_string,  '       CRYSCAL is then reading my_file.RAW and will create a')
  call write_cryscal_news_line(input_string,  '       my_file_RAW.HKL file in a SHELX format. This file also contains')
  call write_cryscal_news_line(input_string,  '       dir. cos. for further absorption correction.')
  call write_cryscal_news_line(input_string,  '       Same behavior can be optained by associating a selected RAW file to')
  call write_cryscal_news_line(input_string,  '       CRYSCAL program in the Windows file folder')
  end if
  if(news_year(1:3) == 'all' .or. news_year(1:2) == '08' .or. news_year(1:4) == '2008') then
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . Nov. 08:')
  call write_cryscal_news_line(input_string,  '     # GEN_HKL accepts now Q_min and Q_max arguments')
  call write_cryscal_news_line(input_string,  '     # SEARCH_HKL arguments can be h, k or l letters.')
  call write_cryscal_news_line(input_string,  '       ex: SEARCH_HKL h 0 0')
  call write_cryscal_news_line(input_string,  '           SEARCH_HKL 1 k 0')
  call write_cryscal_news_line(input_string,  '     # SITE_INFO keyword gets constraints on anisotropic ADP, using the')
  call write_cryscal_news_line(input_string,  '       routine implemented in FullProf')
  call write_cryscal_news_line(input_string,  '     # new available arguments for CRYSCAL in command line:')
  call write_cryscal_news_line(input_string,  '       - LOG    : create a cryscal.log file')
  call write_cryscal_news_line(input_string,  '       - NO_OUT : no output information are written on screen')
 !call write_cryscal_news_line(input_string,  '       - CIFDEP : combined with ARCHIVE.CIF, author name and address'
 !call write_cryscal_news_line(input_string,  '                  are stored in the cryscal_archive.cif file')
 !call write_cryscal_news_line(input_string,  '       - ACTA   : combined with ARCHIVE.CIF, several CIF fields related')
 !call write_cryscal_news_line(input_string,  '                  to CIF file deposition are stored in the')
 !call write_cryscal_news_line(input_string,  '                  cryscal_archive.cif file')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . oct. 08:')
  call write_cryscal_news_line(input_string,  '     # CRYSCAL archive.cif : for archive.cif files created with the new version of')
  call write_cryscal_news_line(input_string,  '       WinGX (sept. 2008), CRYSCAL is skipping the whole part of the cif file ')
  call write_cryscal_news_line(input_string,  '       containing programs references, to created the cryscal_archive.cif file')
  call write_cryscal_news_line(input_string,  ' ')
  call write_cryscal_news_line(input_string,  '     # new argument for CRYSCAL in command line : CRYSCAL CREATE_INS/SIR_TO_INS')
  call write_cryscal_news_line(input_string,  '       Create job.ins file from :')
  call write_cryscal_news_line(input_string,  '         . output.RES file created by SIR programs')
  call write_cryscal_news_line(input_string,  '         . struct.cif file created by WinGX')
  call write_cryscal_news_line(input_string,  '       This job.ins can be directly used for structure refinement with')
  call write_cryscal_news_line(input_string,  "       SHELXL, and contains right cell parameters and esd's, Mo wavelength")
  call write_cryscal_news_line(input_string,  '       and useful instructions such as ACTA, BOND$H, CONF, TEMP ...')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . sept. 08:')
  call write_cryscal_news_line(input_string,  '     # CRYSCAL archive.cif : CRYSCAL is reading archive.cif and cryscal.cif')
  call write_cryscal_news_line(input_string,  '       files and creating completed archive_cryscal.cif file')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '     # keyword PAUSE: pause in th execution of the requested commands')
  call write_cryscal_news_line(input_string,  '       This keyword can be useful when comands are executed from a CFL')
  call write_cryscal_news_line(input_string,  '       commands file')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . july 08:')
  call write_cryscal_news_line(input_string,  '     # REPORT/REPORT_long: CIF fields can be independly in lower and')
  call write_cryscal_news_line(input_string,  '       upper cases')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . april 08:')
  call write_cryscal_news_line(input_string,  '     # CONN keyword: calculation of connectivity around a selected atom')
  call write_cryscal_news_line(input_string,  '       given as argument.')
  end if
  if(news_year(1:3) == 'all' .or. news_year(1:2) == '07' .or. news_year(1:4) == '2007') then
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . sept. 07:')
  call write_cryscal_news_line(input_string,  '     # SEARCH_SPGR keyword: search for a space group, given a hkl list and')
  call write_cryscal_news_line(input_string,  '       a given crystal system (Thanks to JRC for the CHECK_GROUP routine)')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . july 07:')
  call write_cryscal_news_line(input_string,  '     # MATMUL keyword for 3*3 matrix multiplication')
  call write_cryscal_news_line(input_string,  '     # CIF file can be given as a second argument if "REPORT" or "REPORT_long"')
  call write_cryscal_news_line(input_string,  '       is given at first. A CIF_file_structural_report.HTML is then created')
  call write_cryscal_news_line(input_string,  '       ex: d:\cryscal report_long my_CIF_file.CIF')
  call write_cryscal_news_line(input_string,  '     # P4P file can be given as argument in the command line to run CRYSCAL,')
  call write_cryscal_news_line(input_string,  '       i.e. CRYSCAL my_file.P4P.')
  call write_cryscal_news_line(input_string,  '       CRYSCAL is then reading my_file.P4P and my_file.HKL to create')
  call write_cryscal_news_line(input_string,  '       import.CIF file for WinGX')
  call write_cryscal_news_line(input_string,  '       Same behavior can be optained by associating a selected P4P file to')
  call write_cryscal_news_line(input_string,  '       CRYSCAL program in the Windows file folder')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . april 07:')
  call write_cryscal_news_line(input_string,  '     # MERGE keyword: merge equivalent reflections of the current HKL file')
  call write_cryscal_news_line(input_string,  '       ')
  call write_cryscal_news_line(input_string,  '     # REPORT keyword can be interpreted in command line:')
  call write_cryscal_news_line(input_string,  '       \> cryscal report')
  call write_cryscal_news_line(input_string,  '       CRYSCAL is looking in the current folder for the presence of a')
  call write_cryscal_news_line(input_string,  '       "archive.cif" file: "structural_report.html" file is then created ')
  call write_cryscal_news_line(input_string,  '       and contains text about the crystallographic study. The browser')
  call write_cryscal_news_line(input_string,  '       defined in the "cryscal.ini" file is then launch')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . march 07:')
  call write_cryscal_news_line(input_string,  '     # .x file (created by DENZO) and .rmat file (created by DIRAX) can be')
  call write_cryscal_news_line(input_string,  '       passed as argument for CELL keyword')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '   . feb. 07:')
  call write_cryscal_news_line(input_string,  '     # RESET keyword for input parameters and arrays initialization')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '  . jan. 07:')
  call write_cryscal_news_line(input_string,  '    # THERM keyword can performed conversion of anisotropic displacement '  &
                                            //'parameters:')
  call write_cryscal_news_line(input_string,  '      new available arguments: U_ij, B_ij, Beta_ij')
  call write_cryscal_news_line(input_string,  '    # DIR keyword has been added and corresponds to the DIR DOS command.')
  call write_cryscal_news_line(input_string,  '      Arguments may follow this keyword.')
  call write_cryscal_news_line(input_string,  '    # Wcryscal for Windows has been created.')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '  . nov. 06:')
  call write_cryscal_news_line(input_string,  '    # X rays data for Ag, Fe and Cr have been tabulated')
  call write_cryscal_news_line(input_string,  '    # \> cryscal P4P:')
  call write_cryscal_news_line(input_string,  '       CRYSCAL is looking in the current folder, for a P4P file (created')
  call write_cryscal_news_line(input_string,  '       by SAINT) and a HKL file (created by SADABS). Import.cif file is then')
  call write_cryscal_news_line(input_string,  '       created and can be directly read by WinGX as a KappaCCD file')
  endif
  if(news_year(1:3) == 'all' .or. news_year(1:2) == '06' .or. news_year(1:4) == '2006') then
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '  . oct. 06:')
  call write_cryscal_news_line(input_string,  '    # MAG keyword:     output magnetic features for a 3d or 4f ion')
  call write_cryscal_news_line(input_string,  '    # SHANNON keyword: get effective ionic radii from Shannon article ' &
                                            //'(Acta Cryst. 1976, A32, 751)')
  call write_cryscal_news_line(input_string,  '    # P4P keyword:     read P4P file created by SAINT (Bruker-AXS)')
  call write_cryscal_news_line(input_string,  '')
  call write_cryscal_news_line(input_string,  '  . sept. 06:')
  call write_cryscal_news_line(input_string,  '    # FILE keyword: .m91 and .m95 files created by JANA can be read')
  call write_cryscal_news_line(input_string,  '    # CELL keyword: .m50 file can be read')
  call write_cryscal_news_line(input_string,  '')
  endif

 return
end subroutine write_cryscal_NEWS

!--------------------------------------------------------------------------------------------------
 subroutine write_cryscal_news_line(input_string, input_line)
  USE CRYSCAL_module, ONLY : HTML_unit, NEWS_unit
  USE IO_module,      ONLY : write_info

  implicit none
   character (len=*), intent(in)    :: input_string
   character (len=*), intent(inout) :: input_line

  if(input_string(1:6) == 'screen') then
   call write_info(trim(input_line))
  elseif(input_string(1:4) == 'html') then
   write(HTML_unit, '(a)') trim(input_line)
  elseif(input_string(1:4) == 'text') then
   write(2, '(a)') trim(input_line)
  end if

  return

 end subroutine write_cryscal_NEWS_line

!-----------------------------------------------------------------------------------------------------

!subroutine write_cryscal_CLA
! USE text_module
! USE IO_module,      ONLY : write_info
! implicit none
!  integer       :: i, k


!  call Def_command_line_arguments
!
!  call write_info('')
!  call write_info(' List of CRYSCAL command line arguments : ')
!
!  do i=1, CLA_nb
!   call write_info(trim(CLA_line(i,1)))
!   do k=2, CLA_lines_nb(i)
!    call write_info(TRIM(CLA_line(i,k)))
!   end do
!  end do
!  call write_info('')
!
! return
!end subroutine write_cryscal_CLA
