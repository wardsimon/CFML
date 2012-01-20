!     Last change:  TR    2 May 2008   12:34 pm
!     Last change:  TR    2 May 2008   12:34 pm
!--------------------------------------------------------------------

module text_module
 USE cryscal_module, ONLY: nb_help_max
 implicit none
 INTEGER, parameter                    :: title_lines_nb      = 19
 INTEGER, parameter                    :: header_lines_nb     = 58
 INTEGER, parameter                    :: CIF_lines_nb        = 10


 CHARACTER (LEN=80), DIMENSION(title_lines_nb)       :: title_line
 CHARACTER (LEN=80), DIMENSION(header_lines_nb)      :: header_line
 CHARACTER (LEN=80), DIMENSION(CIF_lines_nb)         :: CIF_title_line

 INTEGER,            DIMENSION(nb_help_max)          :: HELP_lines_nb
 CHARACTER (LEN=80), DIMENSION(nb_help_max,100)      :: HELP_line

 INTEGER, parameter                                  :: CLA_nb = 13
 INTEGER,            DIMENSION(CLA_nb)               :: CLA_lines_nb
 CHARACTER (LEN=80), DIMENSION(CLA_nb, 10)           :: CLA_line



end module text_module
!-------------------------------------------------------------------------------

subroutine def_title_lines
 USE CRYSCAL_module, ONLY : cryscal_version
 USE TEXT_module,    ONLY : title_line

 title_line(1) =' '
 title_line(2) ='  *************************************************************'
 title_line(3) ='  *                                                           *'
 title_line(4) ='  *                 CRYSTALLOGRAPHIC CALCULATIONS             *'
 title_line(5) ='  *                                                           *'
 title_line(6) ='  *                                                           *'
 title_line(7) ='  *                (T.R. / CDIFX - UMR6226 Rennes)            *'
 title_line(8) ='  *                           '// cryscal_version(1:16)// '                *'
 title_line(9) ='  *                                                           *'
 title_line(10)='  *            with courtesy of JRC and JGP for CFML          *'
 title_line(11)='  *                                                           *'
 title_line(12)='  *                                                           *'
 title_line(13)='  *              (beta version still under tests)             *'
 title_line(14)='  *               contact: cdifx@univ-rennes1.fr              *'
 title_line(15)='  *                                                           *'
 title_line(16)='  *      Web site : www.cdifx.univ-rennes1.fr/cryscal.htm     *'
 title_line(17)='  *                                                           *'
 title_line(18)='  *************************************************************'
 title_line(19)=' '

end subroutine  def_title_lines

!-----------------------------------------------------------------------------------------
subroutine def_CIF_file_title
 use text_module, only : CIF_title_line

 CIF_title_line    = ''

 CIF_title_line(1) = '###########################################################'
 CIF_title_line(2) = '#                   CIF file created by CRYSCAL           #'
 CIF_title_line(3) = '#                  T.R. / CDIFX Rennes / 2006-11          #'
 CIF_title_line(4) = '#                www.cdifx.univ-rennes1.fr/cryscal        #'
 CIF_title_line(5) = '###########################################################'

end subroutine def_CIF_file_title
!-------------------------------------------------------------------------------


subroutine def_header_lines()
 USE text_module, ONLY :    header_line
 integer               :: n

 n= 0
 n = n + 1; header_line(n) = ' '
 n = n + 1; header_line(n) = ' '
 n = n + 1; header_line(n) = '   CRYSCAL has been created to perform basic crystallographic calculations or'
 n = n + 1; header_line(n) = '   get crystallographic informations.  '
 n = n + 1; header_line(n) = '   CRYSCAL has been written in Fortran 95, and uses the crystallographic '
 n = n + 1; header_line(n) = '   calculations facilities of the Crystallographic Fortran Modules Librairies'
 n = n + 1; header_line(n) = '   written by J. Rodriguez-Carvajal (ILL-Grenoble, France) and J. Gonzalez'
 n = n + 1; header_line(n) = '   (Univ. La Laguna, Spain).'
 n = n + 1; header_line(n) = '    '
 n = n + 1; header_line(n) = ' '
 n = n + 1; header_line(n) = '    Facilities implemented in CRYSCAL:'
 n = n + 1; header_line(n) = '     . unit cell volume calculation'
 n = n + 1; header_line(n) = '     . space group informations: space group features, Wyckoff positions, '
 n = n + 1; header_line(n) = '       symmetry operators, extinctions ...'
 n = n + 1; header_line(n) = '     . calculation of d_hkl, Q_hkl, 2theta_hkl (including a propagation wave'
 n = n + 1; header_line(n) = '       vector)'
 n = n + 1; header_line(n) = '     . hkl generation for a given space group' 
 n = n + 1; header_line(n) = '     . structure factor calculation (Xrays, neutrons, electrons)'
 n = n + 1; header_line(n) = '     . simulation of pouwder diffraction pattern (X, neutrons)'
 n = n + 1; header_line(n) = '     . geometric calculations: interatomic distances, angles, connectivity, '
 n = n + 1; header_line(n) = '       centroid coordinates, angles between 2 vectors in direct and reciprocal '
 n = n + 1; header_line(n) = '       space, ...'
 n = n + 1; header_line(n) = '     . atomic features: weight, density, electronic configuration, ionic and'
 n = n + 1; header_line(n) = '       Shannon radii, neutron data, X-rays data ...'
 n = n + 1; header_line(n) = '     . molecular informations: molecular weight, density ...'
 n = n + 1; header_line(n) = '     . absorption coefficient calculation (X-rays, neutrons) and '
 n = n + 1; header_line(n) = '       transmission calculation'
 n = n + 1; header_line(n) = '     . transformation of unit cell, atomic coodinates and hkl files '
 n = n + 1; header_line(n) = '     . statistics on hkl file and sort of hkl data'
 n = n + 1; header_line(n) = '     . search for systematic extinctions and space group'
 n = n + 1; header_line(n) = '     ... '
 n = n + 1; header_line(n) = '      '
 n = n + 1; header_line(n) = '   CRYSCAL can be run through an input file containing a list of keywords,'
 n = n + 1; header_line(n) = '   defining the type of crystallographic calculations that will be performed,'
 n = n + 1; Header_line(n) = '   or in an interactive mode, by entering keywords at the CRYSCAL prompt. '
 n = n + 1; header_line(n) = '      '
 n = n + 1; header_line(n) = '   Crystallographic features can be read from different types of input files : '
 n = n + 1; header_line(n) = '                                               .CFL'
 n = n + 1; header_line(n) = '                                               .INS/.RES file (SHELXL)'
 n = n + 1; header_line(n) = '                                               .CIF file'
 n = n + 1; header_line(n) = '                                               .PCR file (FullProf)'
 n = n + 1; header_line(n) = '                                               .CEL file (PowderCELL)'
 n = n + 1; header_line(n) = '    '
 n = n + 1; header_line(n) = '   Alternatively, particular jobs can be performed by CRYSCAL when special'
 n = n + 1; header_line(n) = '   arguments are passed to CRYSCAL through the command line (see CRYSCAL'
 n = n + 1; header_line(n) = '   command line arguments section).'
 n = n + 1; header_line(n) = '    '
 n = n + 1; header_line(n) = '   Online help can be obtained by typing "MAN" or "HELP" at the "Enter input'
 n = n + 1; Header_line(n) = '   file" (menu option #1) or "Enter keyword:" (menu option #2) prompt:'
 n = n + 1; header_line(n) = '      d:> cryscal'
 n = n + 1; header_line(n) = '          > Enter keyword : man'
 n = n + 1; header_line(n) = '   or launching CRYSCAL program with "MAN" or "HELP" as argument:'
 n = n + 1; header_line(n) = '      d:> cryscal man'
 n = n + 1; header_line(n) = '   Details on the meaning of keywords can be obtained by typing the'
 n = n + 1; header_line(n) = '   corresponding keyword at the "Enter keyword :" prompt or launching CRYSCAL '
 n = n + 1; header_line(n) = '   program with the corresponding keyword(s) as argument(s):'
 n = n + 1; header_line(n) = '      d:> cryscal MAN CELL'
 n = n + 1; header_line(n) = ''

end subroutine def_header_lines
!-------------------------------------------------------------------------------


subroutine def_keywords_lines()
 USE text_module
 USE cryscal_module
  implicit none
  integer                    :: n, HELP_numor

 ! 1
  n = 0    ;  HELP_numor = HELP_ABSENT_HKL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > ABSENT_HKL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             list observed reflections (with F2>0.) that shoud'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           be absent for a given space group'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal arguments: . arg = "ALL": all the violations reflections are'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          output'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = "OUT": requested reflections are output on'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          the screen'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = "WRITE": requested reflextions are output'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                            in a HKL file'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = real_value (n_sig): only reflections with'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                                       I/sig > n_sig are output'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE, SPGR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  ABSENT_HKL, HKL_ABSENT'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  
  ! 2
  n = 0    ;  HELP_numor = HELP_ABSORPTION_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > ABSORPTION:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             absorption coefficient calculation'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keywords:  CELL, WAVE, CONT / CHEM'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  ABSORPTION, ABSORPTION_CALC, CALC_ABSORPTION'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MU, CALC_MU, MU_CALC' 
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n
  

 ! 3
  n = 0    ;  HELP_numor = HELP_ACTA_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > ACTA:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             create CRYSCAL.CIF file containing calculation'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           results in a CIF format'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  ACTA, CIF, CREATE_CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 4
  n = 0    ; HELP_numor = HELP_ANG_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > ANG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           3 character strings '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             calculation of the angle defined by three atoms'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           labelled by the atom labels (cf ATOM keyword)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             ANG C8 C10 C9'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  or '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           4 character strings '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             calculation of the angle defined by the segments '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           constituted by atoms couples labelled as S_1, S_2 '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           and S_3, S_4 respectively (cf ATOM keyword)'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . remark:              atom label can refer to equivalent position through'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           a particular symmetry operator (cf SYMM keyword). '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           This is specified by adding "_$n" to the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           atomic label, with n refering to the number of '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           symmetry operator in the list.'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             ANG O1 C8 C9 C10_$1'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  ANG, ANGLE'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 5
  n = 0    ; HELP_numor = HELP_APPLY_OP_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > APPLY_OP:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:              apply the symmetry operators on the atomic '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           positions'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SYMM, ATOM'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  APPLY_OP, APPLY_SYMMETRY_OPERATOR,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           APPLY_SYM_OP, APPLY_SYMOP'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 6
  n = 0    ;  HELP_numor = HELP_ATOM_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > ATOM:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           2 character strings and 5 reals '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             . string #1: atomic label'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . string #2: atom type (see SFAC)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . x,y,z atomic reduced coordinates'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . Biso, site occupancy (%)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . remarks:             . if Biso   is missing: Biso = 0.0'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if Occ(%) is missing: Occ(%) = 1.0 (site is fully '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             occupied)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             ATOM O1   O   0.04356  0.03008  0.39001   0.35  1.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           ATOM C8   C   0.02071 -0.12436  0.36957   0.30  1.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           ATOM C9   C  -0.27497 -0.07538  0.27585   0.30  1.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           ATOM C10  C  -0.16896 -0.18823  0.36382   0.32  1.'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  ATOM, ATM'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor ) = n

 ! 7
  n= 0     ;  HELP_numor = HELP_ATOM_LIST_numor
  n = n + 1;  HELP_line(HELP_numor, n)   = ''
  n = n + 1;  HELP_line(HELP_numor, n)   = ' > ATOM_LIST:'
  n = n + 1;  HELP_line(HELP_numor, n)   = ''
  n = n + 1;  HELP_line(HELP_numor, n)   = '   . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n)   = '   . arguments:           no'
  n = n + 1;  HELP_line(HELP_numor, n)   = '   . meaning:             list the atoms (type, labels, coordinates, ...)'
  n = n + 1;  HELP_line(HELP_numor, n)   = '   . optional keyword:    SPGR'
  n = n + 1;  HELP_line(HELP_numor, n)   = '   . identical keywords:  ATOM_LIST, LIST_ATOM_LIST, LIST_ATOMS, LST_ATOMS,'
  n = n + 1;  HELP_line(HELP_numor, n)   = '                          WRITE_ATOMS, WRITE_ATMS'
  n = n + 1;  HELP_line(HELP_numor, n)  = ''
  HELP_lines_nb(HELP_numor) = n

 ! 8
  n = 0    ;  HELP_numor = HELP_BARY_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > BARY:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           n character strings '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             calculation of the coordinates of the centroid of the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           n atoms known by the atom label S (cf ATOM keyword)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             BARY C8 C10 C9 C11'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  BARY, CENTROID'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 9
  n = 0    ;  HELP_numor = HELP_BEAM_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > BEAM:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           character string '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             type of the incident beam:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - BEAM NEUT  for neutrons'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - BEAM ELECTRONS  for electrons'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - X_Ag for X Rays (Silver K_alpha)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - X_Mo for X Rays (Molybdenum K_alpha)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - X_Cu for X Rays (Copper K_alpha)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - X_Ni for X Rays (Cobalt K_alpha)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - X_Co for X Rays (Cobalt K_alpha)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - X_Fe for X Rays (Iron K_alpha)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - X_Cr for X Rays (Chromium K_alpha)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:            BEAM X_Mo'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 10
  n = 0    ;  HELP_numor = HELP_CELL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CELL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           6 reals or 1 character string'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:              - reals: unit cell parameters in A (a, b, c) and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                     angle in deg. (alfa, beta, gamma)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            - character string: file name containing unit cell'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              parameters.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              Following files can be read:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . INS/RES file for SHELXL'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . P4P  file created by SAINT'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . M50  file created by JANA'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . X    file created by DENZO'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . RMAT file created by DIRAX'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              unit cell volume calculation'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           direct and reciprocal unit cell parameters'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             CELL 7.6520 7.8450 11.0760 90. 90. 90.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           CELL my_saint_data.P4P'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           CELL my_janat_data.m50'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords   CELL, CELL_PARAMETERS'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 11
  n = 0    ;  HELP_numor = HELP_CHEM_numor
  n = n + 1;  HELP_line(HELP_numor, n)  = ' '
  n = n + 1;  HELP_line(HELP_numor, n)  = '  > CHEM:'
  n = n + 1;  HELP_line(HELP_numor, n)  = ' '
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . arguments:           n "El_i_n_i" character strings (without blank'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                           character between label and number)'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . meaning:             Molecular chemical formula: El_i is the chemical'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                           symbol of the species i and n_i is the corresponding'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                           number of atoms in the formula unit'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . output:              molecular weight, total number of electrons,'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                           atomic and weight percentage'
  n = n + 1;  HELP_line(HELP_numor, n)  = ' '
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . mandatory keyword:   ZUNIT'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . example:             CHEM C4 O6 H9 N1'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . identical keywords:  CHEM, CHEM_FORM, CHEMICAL_FORMULA'
  n = n + 1;  HELP_line(HELP_numor, n)  = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 12
  n = 0    ;  HELP_numor = HELP_CONN_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CONN:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            atom_label + dist_max'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             determine the connectivity around the atom with'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           interatomic distances calculated until dist_max value.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           dist_max (< 15. A) is optionnal. Default value = 3.0 A.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           if atom_label=ALL, the program will calculate'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           connectivity around all atoms.'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              interatomic distances. IF CIF/ACTA keyword is'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           entered, the created CIF file will contains all the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculated distances in CIF format.'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keywords:  SPGR, ATOM, CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  CONN, CONNECT, CONNECTIVITY'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n



 ! 13
  n = 0    ;  HELP_numor = HELP_CONT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CONT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           n "El_i  n_i"(character string, real) couples'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             unit cell contents: El_i is the chemical symbol of'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           the species i and n_i is the corresponding number'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           of atoms in the unit cell'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional keyword:    ZUNIT'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             CONT C 16.  O 24.  H 36.  N  4.'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


  ! 14 
  n = 0    ;  HELP_numor = HELP_CREATE_ACE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CREATE_ACE :'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             create .ACE file for CaRIne from a CIF file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   READ_CIF file.cif'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . dependent parameter: CREATE_ACE parameter value in the CRYSCAL.INI'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           setting file. If equal to 1, a .ACE file will'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           be automatically created if a .CIF file is  '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           given as argument when CRYSCAL is launching '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           from a command line :'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           > cryscal file.cif '
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  ! 15
  n = 0    ;  HELP_numor = HELP_CREATE_CEL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CREATE_CEL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             create .CEL file for PowderCELL from a CIF file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   READ_CIF file.cif'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . dependent parameter: CREATE_CEL parameter value in the CRYSCAL.INI'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           setting file. If equal to 1, a .CEL file will'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           be automatically created if a .CIF file is  '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           given as argument when CRYSCAL is launching '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           from a command line :'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           > cryscal file.cif '
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  ! 16
  n = 0    ;  HELP_numor = HELP_CREATE_CFL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CREATE_CFL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             create .CFL file for CRYSCAL from a CIF file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   READ_CIF file.cif'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . dependent parameter: CREATE_CFL parameter value in the CRYSCAL.INI'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           setting file. If equal to 1, a .CFL file will'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           be automatically created if a .CIF file is  '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           given as argument when CRYSCAL is launching '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           from a command line :'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           > cryscal file.cif '
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  ! 17 
  n = 0    ;  HELP_numor = HELP_CREATE_FST_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CREATE_FST'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             create .FST file for FP Studio'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   READ_CIF file.cif'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . dependent parameter: CREATE_fst parameter value in the CRYSCAL.INI'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           setting file. If equal to 1, a .fst file will'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           be automatically created if a .CIF file is  '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           given as argument when CRYSCAL is launching '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           from a command line :'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           > cryscal file.cif '
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  ! 18 
  n = 0    ;  HELP_numor = HELP_CREATE_INS_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CREATE_INS'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             create .INS file for SHELXL from a CIF file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   READ_CIF file.cif'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . dependent parameter: CREATE_INS parameter value in the CRYSCAL.INI'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           setting file. If equal to 1, a .INS file will'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           be automatically created if a .CIF file is  '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           given as argument when CRYSCAL is launching '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           from a command line :'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           > cryscal file.cif '
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 19
  n = 0    ;  HELP_numor = HELP_CREATE_REPORT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > CREATE_REPORT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                output keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional argument:   .CIF file name'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             create a STRUCTURAL_REPORT.HTML file in HTML format'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           from the reading of the ARCHIVE.CIF file present in'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           the current folder, and launch the browser with this'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           HTML file. The .CIF file can be explicitely defined'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           with the argument.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           If "long" or "ext" is given as argument, a longer'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           report will be created, containing more informations,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           included distances and angles.'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  REPORT, CREATE_REPORT'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            report'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           report long my_struct.cif'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 20
  n = 0    ;  HELP_numor = HELP_D_HKL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > D_HKL: '
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            real values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             d_hkl(A) values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             Q(A-1), SinTheta/lambda(A-1)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           theta(deg) for known wavelength'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal keyword:   WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  D_HKL, DHKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             DHKL 0.77'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 21
  n = 0    ;  HELP_numor = HELP_D_STAR_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > D_STAR: '
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            real values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             1/d_hkl (A-1'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             d(A)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           theta(deg)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal keyword:   WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  D_STAR, D_STAR_HKL, DSTAR, DSTARHKL, DSTAR_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             D_STAR  0.5'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 22
  n = 0    ;  HELP_numor = HELP_DATA_ATOMIC_DENSITY_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > DATA_ATOMIC_DENSITY:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '     . type:               OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . optionnal argument: PLOT'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . meaning:            list atomic density data for all atoms '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           if arg=PLOT: create a PGF file and plot it with'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           WinPLOTR'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . identical keywords: DATA_DENSITY, DENSITY_DATA, DATA_ATOMIC_DENSITY,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           ATOMIC_DENSITY'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n


 ! 23
  n = 0    ;  HELP_numor = HELP_DATA_ATOMIC_RADIUS_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > DATA_ATOMIC_RADIUS:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '     . type:               OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . optionnal argument: PLOT'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . meaning:            list atomic radius data for all atoms '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           if arg=PLOT: create a PGF file and plot it with'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           WinPLOTR'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . identical keywords: DATA_RADIUS, RADIUS_DATA, DATA_ATOMIC_RADIUS,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           ATOMIC_RADIUS'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 24
  n = 0    ;  HELP_numor = HELP_DATA_ATOMIC_WEIGHT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > DATA_ATOMIC_WEIGHT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '     . type:               OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . optionnal argument: PLOT'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . meaning:            list atomic weight data for all atoms '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           if arg=PLOT: create a PGF file and plot it with'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           WinPLOTR'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . identical keywords: DATA_WEIGHT, WEIGHT_DATA, DATA_ATOMIC_WEIGHT,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           ATOMIC_WEIGHT'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n


 ! 25
  n = 0    ;  HELP_numor = HELP_DATA_NEUTRONS_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > DATA_NEUTRONS:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '     . type:               OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . optionnal argument: PLOT'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . meaning:            list neutrons data for all atoms (coherent'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           scattering length, incoherent scattering '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           cross-section, absorption cross-section)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           Neutron data are extracted from :'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           V.F. Sears'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           Neutron News, vol.3 n°3, 1992, 26-37'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           if arg=PLOT: create a PGF file and plot it with'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           WinPLOTR'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . see:                Xrays_DATA'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . identical keywords: DATA_NEUTRONS, NEUTRONS_DATA, DATA_NEUTRON,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           NEUTRON_DATA'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 26
  n = 0    ;  HELP_numor = HELP_DATA_XRAYS_numor
  n = n + 1;  HELP_line(HELP_numor, n)  = ''
  n = n + 1;  HELP_line(HELP_numor, n)  = '  > DATA_XRAYS:'
  n = n + 1;  HELP_line(HELP_numor, n)  = ''
  n = n + 1;  HELP_line(HELP_numor, n)  = '     . type:               OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n)  = '     . optionnal argument: PLOT'
  n = n + 1;  HELP_line(HELP_numor, n)  = '     . meaning:            list X-ray data for all atoms (total interaction'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                           cross section for Ag, Mo, Cu, Co, Fe and Cr radiations)'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                           X-ray data are extracted from :'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                            Tables Internationales vol.C  1995, p.200-206,'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                            Tables Internationales vol.C  1995, p. 193-199'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                           if arg=PLOT: create a PGF file and plot it with'
  n = n + 1;  HELP_line(HELP_numor, n)  = '                           WinPLOTR'
  n = n + 1;  HELP_line(HELP_numor, n)  = '     . see:                Xrays_DATA'
  n = n + 1;  HELP_line(HELP_numor, n)  = '     . identical keywords: DATA_XRAYS, XRAYS_DATA, DATA_XRAY, XRAY_DATA'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = 12

 
  ! 27
  n = 0    ;  HELP_numor = HELP_DIAG_MAT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > DIAG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           9 reals '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             transformation (3,3) matrix components'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            MATR   0  0  1    0  1  0    -1  0  -1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MATR   0.5 0.5 0   -0.5 0.5 0   0 0 1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MATR   1/2 1/2 0   -1/2 1/2 0   0 0 1'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              Diagonalization of the 3*3 matrix and output'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           the Eigen values and Eigen vectors'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           (dependent keyword = FILE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  DIAG, DIAG_MAT, DIAG_MATR, DIAG_MATRIX'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 
 
 ! 28
  n = 0    ;  HELP_numor = HELP_DIR_ANG_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > DIR_ANG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           2*3 reals'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             calculation of the angle between 2 vectors in the '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           direct space. The 3 first real values are related '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           to the coordinates of the first vector and the 3 '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           last real values to the coordinates of the second '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           vector'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  DIR_ANG, DIRANG, DIRECT_ANGLE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             DIR_ANG 1. 0. 0.   0. 1. 0.'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 29
  n = 0    ;  HELP_numor = HELP_DIST_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > DIST:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           2 character strings '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             calculation of the interatomic distance between 2'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           atoms labelled by their atom labels (cf ATOM keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . remark:              atom label can refer to equivalent position through'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           a particular symmetry operator (cf SYMM keyword)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           atomic label, with n refering to the number of '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           symmetry operator in the list.'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             DIST C8 C10'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           DIST C10 C9'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           DIST C10 C9_$1'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  DIST, DISTANCE, ATOMIC_DISTANCE'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 30
  n = 0    ;  HELP_numor = HELP_EDIT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > EDIT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                output keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            1 character string (filename)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             open the file given as argument with the editor defined'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           in the CRYSCAL.ini setting file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             edit my_cryscal.cfl'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . related keywords:    SET, SETTING'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 31
  n = 0    ;  HELP_numor = HELP_EQUIV_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > EQUIV:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           3 integer values: h, k, l'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal arguments: FRIEDEL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             search for equivalent reflections (Space group is'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           mandatory)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = FRIEDEL: search for Friedel reflections'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              list of required hkl reflections'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SPGR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  EQUIV, EQUIV_HKL, SEARCH_EQUIV, SEARCH_EQUIV_HKL, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           FIND_EQUIV, FIND_EQUIV_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             EQUIV 1 1 0 FRIEDEL'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 32 
  n = 0    ;  HELP_numor = HELP_EXIT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > EXIT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '     . type:               OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . argument:           no'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . meaning:            exit from the "enter keyword" procedure to come'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           back to the CRYSCAL main menu'
  n = n + 1;  HELP_line(HELP_numor, n) = '     . identical keywords: EXIT, X, QUIT, END, STOP'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) =  n

 ! 33
  n = 0    ;  HELP_numor = HELP_FILE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > FILE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           1 character string reals '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             read .HKL file name,  containing a list of '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           h, k, l, F2,  sig_F2'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . available formats for HKL file:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . SHELX type: 3I4, 2F8.2'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . CIF format'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . final.y (created by EVALCCD)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . RAW format (created by SAINT)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               . M91/M95 format (created by JANA)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              . analysis of the hkl reflections:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . number of hkl reflections'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . d(A) and sinTheta/lambda collected ranges '
  n = n + 1;  HELP_line(HELP_numor, n) = '                              (dependent keyword = CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . Theta collected range (dependent keyword = CELL, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                              WAVE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . statistic on collected data'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if space group is known (see SPGR keyword), the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             number of reflections in agreement is output'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal arguments: NEG: negative reflections are treated as follows:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             if (F2      < 0.0001)  F2 = 0.0001'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             if (sig_F2_ < 0.00001) sig_F2_=sqrt(abs(F2)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           PLOT (dependent keywords = CELL, WAVE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . action:              plot the F2=f(sinTheta/lambda) curve with WinPLOTR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            FILE my_HKL_FILE.HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           FILE import.CIF plot'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 34
  n = 0    ;  HELP_numor = HELP_FIND_HKL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > FIND_HKL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           3 integer values: h, k, l'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           or'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           characters with one or two index characters being h,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           k or l'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal arguments: EQUIV, FRIEDEL (only for integer values of h, k and l'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             search a particular hkl reflection in a '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           reflections list'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = EQUIV:   search for equivalent reflections '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (Space group is mandatory)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = FRIEDEL: search for Friedel reflections'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              list of required hkl reflections'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  FIND_HKL, FINDHKL, SEARCH_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            FIND_HKL 1 1 0 FRIEDEL'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           FIND_HKL h 0 0'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           FIND_HKL 2 k 0'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           FIND_HKL -1 k l'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n


 ! 35
  n = 0    ;  HELP_numor = HELP_FIND_HKL_LIST_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > FIND_HKL_LIST:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '   . type:                 OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '   . mandatory argument:   1 integer value corresponding to the numor in the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           following the extinction rules list:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               1. h00 h=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               2. 0k0 k=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               3. 00l l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               4. 0kl k=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               5. 0kl l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               6. 0kl k+l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               7. h0l h=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               8. h0l l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               9. h0l h+l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              10. hk0 h=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              11. kk0 k=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              12. hk0 h+k=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              13. hhl h+l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              14. hkk k+h=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              15. hkh h+k=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              16. hkl k+l=2n+1 (A)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              17. hkl h+l=2n+1 (B)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              18. hkl h+k+2n+1 (C)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              19. hkl h+k+l=2n+1 (I)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              20. hkl not all odd/even (F)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              21. h00 h=4n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              22. 0k0 k=4n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              23. 00l l=4n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              24. 0kl k+l=4n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              25. h0l h+l=4n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              26. hk0 h+k=4n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              27. h-hl, h0l, 0kl l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              28. hkl h=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              29. hkl k=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              30. hkl l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              31. hkl h=2n+1,k=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              32. hkl h=2n+1,l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              33. hkl k=2n+1,l=2n+1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            or a character string in the following list :'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            A, B, C, I, F corresponding to the #13, #14, #15,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            #16 and #17 in the previous list.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            If the integer value is negative, the opposite'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            condition will be searched in the HKL file.'
  n = n + 1;  HELP_line(HELP_numor, n) = '   . meaning:               search, in a reflections list, those obeying '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            to the selected extinction rule'
  n = n + 1;  HELP_line(HELP_numor, n) = '   . optionnal arguments:  . arg = "OUT": the requested list will be output on '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          the screen'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = "ALL": no condition about F2 and F2/sig is '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          required'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = "WRITE": requested list is stored in a HKL '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                            file'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = "SUPPRESS"/"REMOVE": requested list is'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                    removed from the initial file'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . arg = real value: F2/sig value'
  n = n + 1;  HELP_line(HELP_numor, n) = '   . mandatory keyword:   FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '   . examples:            FIND_HKL_LIST 5 OUT ALL'
  n = n + 1;  HELP_line(HELP_numor, n) = '                          FIND_HKL_LIST I WRITE SUPPRESS'
  n = n + 1;  HELP_line(HELP_numor, n) = '   . identical keywords:  FIND_HKL_LIST, FIND_HKL_LST, EXTRACT_HKL_LIST, EXTRACT_HKL_LST'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor ) = n



 ! 36
  n = 0    ;  HELP_numor = HELP_GEN_HKL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > GEN_HKL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                calculation keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           character strings '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             generate a hkl reflections list in a particular '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           scattering range:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . STL_min=xx  STL_max=xx for SinTheta/lambda '
  n = n + 1;  HELP_line(HELP_numor, n) = '                              range, in A-1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . d_min=xx    d_max=xx   for d range, in A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . Q_min=xx    Q_max=xx   for Q=4pi*SinTheta/lambda'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              range, in A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . Theta_min=xx  Theta_max=xx    for Theta range, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                              in deg.(dependent keyword = WAVE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . 2Theta_min=xx 2Theta_max=xx   for 2Theta range, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                              in deg. (dependent keyword = WAVE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           If atoms are input, a structure factor calculation is'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           done.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           In the case where scattering variable is 2theta, a'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            list of I/Imax is output.'  
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional arguments    "OUT": a list of generated reflections will be ouput'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            "PAT": In the case where scattering variable is '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            2theta, atoms are input and space group is know,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            a diffraction pattern is then created. Two different'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            files are ouput :' 
  n = n + 1;  HELP_line(HELP_numor, n) = '                                - cryscal_pat.xy (2 columns)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                - cryscal_pat.prf (FullProf format)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            If WinPLOTR is installed on the workstation, it is '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            automatically launched with the PRF file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              - calculation of interplanar distance d_hkl, in A '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keyword = CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - calculation of Bragg angle 2theta_hkl, in deg.  '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keyword = CELL, WAVE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - structure factor calculation'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywors = ATOM)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keywords:  CELL, SPGR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples             GEN_HKL theta_min=0. THETA_max=25. OUT'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           GEN_HKL d_min=0.5.   d_max=5.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           GEN_HKL 2THETA_min=0. 2THETA_max=125. OUT PAT'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  GEN_HKL, GENERATE_HKL, GENERATE_HKL_LIST'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 37
  n = 0    ;  HELP_numor = HELP_HEADER_numor
  n = n + 1;  HELP_line(HELP_HEADER_numor, n) = ''
  n = n + 1;  HELP_line(HELP_HEADER_numor, n) = ' > HEADER:'
  n = n + 1;  HELP_line(HELP_HEADER_numor, n) = ''
  n = n + 1;  HELP_line(HELP_HEADER_numor, n) = '   . type:                 OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_HEADER_numor, n) = '   . argument:             no'
  n = n + 1;  HELP_line(HELP_HEADER_numor, n) = '   . meaning:              write header text of CRYSCAL'
  n = n + 1;  HELP_line(HELP_HEADER_numor, n) = '   . identical keywords:   HEADER, HEAD'
  n = n + 1;  HELP_line(HELP_HEADER_numor, n) = ' '
  HELP_lines_nb(HELP_HEADER_numor) = n


 ! 38
  n = 0    ;  HELP_numor = HELP_HEX_RHOMB_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > HEX_RHOMB:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             give transformation matrix from hexagonal to '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           rhomboedral setting'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  HEX_RHOMB, HEXA_RHOMB, HEX_TO_RHOMB, HEXA_TO_RHOMB'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) =  n


 ! 39
  n = 0    ;  HELP_numor = HELP_HKL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > HKL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           3 reals '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             h,k,l, Miller indices'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              calculation of interplanar distance d_hkl, in A '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           (dependent keyword = CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculation of Bragg angle 2theta_hkl, in deg.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           (dependent keyword = CELL, WAVE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             HKL 1. 0. 1.'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 40
  n = 0    ;  HELP_numor = HELP_HKL_NEG_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > HKL_NEG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:               output keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:          no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:            search for negative intensity reflections in a HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '                          file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal keywords: . arg = "OUT": requested reflections are output on the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                         screen'
  n = n + 1;  HELP_line(HELP_numor, n) = '                          . arg = "WRITE": requested reflextions are output in a'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                           HKL file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:  FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords: HKL_NEG, HKL_NEGATIVE, NEG_HKL, NEGATIVE_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . relared keyword:    HKL_POS'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 41
  n = 0    ;  HELP_numor = HELP_HKL_POS_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > HKL_POS:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:               output keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:          no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:            search for positive intensity reflections in a HKL '
  n = n + 1;  HELP_line(HELP_numor, n) = '                          file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal keywords: . arg = "OUT": requested reflections are output on the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                         screen'
  n = n + 1;  HELP_line(HELP_numor, n) = '                          . arg = "WRITE": requested reflextions are output in a'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                           HKL file'
  n = n + 1;  HELP_line(HELP_numor, n) = '                          . arg = real_value (n_sig): only positive reflections'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                  with I > n_sig*sig(F2) are output'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:  FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords: HKL_POS, HKL_POSITIVE, POS_HKL, POSITIVE_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . relared keyword:    HKL_NEG'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 42
  n = 0    ;  HELP_numor = HELP_INSIDE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > INSIDE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                output keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           no '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             put the atoms of the atom list inside the unit cell'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              list of atoms with atomic coordinates in the 0.0 - '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           1.0 range'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   ATOM or atom list in a .INS/.RES file'
  n = n + 1;  HELP_line(HELP_numor, n ) = ' '
  HELP_lines_nb(HELP_numor) = n



 ! 43
  n = 0    ;  HELP_numor = HELP_LIST_EXTI_numor
  n = n + 1;  HELP_line(HELP_LIST_EXTI_numor, n) = ''
  n = n + 1;  HELP_line(HELP_LIST_EXTI_numor, n) = ' > LIST_EXTI_RULE:'
  n = n + 1;  HELP_line(HELP_LIST_EXTI_numor, n) = ''
  n = n + 1;  HELP_line(HELP_LIST_EXTI_numor, n) = '   . type:                 OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_LIST_EXTI_numor, n) = '   . argument:             no'
  n = n + 1;  HELP_line(HELP_LIST_EXTI_numor, n) = '   . meaning:              list the extinction rules implemented in CRYSCAL'
  n = n + 1;  HELP_line(HELP_LIST_EXTI_numor, n) = '   . identical keywords:   LIST_EXTI, LIST_EXTI_RULE, LST_EXTI, LST_EXTI_RULE'
  n = n + 1;  HELP_line(HELP_LIST_EXTI_numor, n) = ''
  HELP_lines_nb(HELP_LIST_EXTI_numor) = n

 ! 44
  n = 0    ;  HELP_numor = HELP_LIST_KEYS_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > LIST_KEYS:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  character string containing "*" character'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             list of keywords:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if no argument: all the keywords are listed'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if "*" character is present in the character '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             string: all the keywords containing the '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             the character string before or after the "*" '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             will be listed'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  KEY, KEYS, LST_KEYS, LIST_KEYS, LST_KEYWORDS, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           LIST_KEYWORDS'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            . KEYS'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . KEYS HKL*'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor ) = n

 ! 45
  n = 0    ;  HELP_numor = HELP_LIST_LAUE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > LIST_LAUE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             list of Laue classes'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  "LST_LAUE", "LIST_LAUE, "LST_LAUE_CLASS",'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           "LIST_LAUE_CLASS"'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            . LIST_LAUE'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor ) = n


 ! 46
  n = 0    ;  HELP_numor = HELP_LIST_MATR_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > LIST_MAT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             list of transformation matrices implemented in CRYSCAL'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           (direct and inverse matrices)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  LST_MAT, LST_MATR, LST_MATRIX, LIST_MAT, LIST_MATR,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           LIST_MATRIX, LIST_TRANSFORMATION_MATRIX'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) =  n

 ! 47
  n = 0    ;  HELP_numor = HELP_LIST_SG_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > LIST_SG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           character strings'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           available arguments (order is not important):'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . symmetry arguments: "all", "centric", "acentric", '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             "triclinic", "monoclinic", "orthorhombic", '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             "tetragonal","trigonal", "hexagonal", "cubic"'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . "chiral", "enantio", "polar"'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . Bravais arguments: "P", "A", "B", "C", "I", '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             "F", "R" '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . Laue class arguments: "laue_n", with n relative'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             to the numor of Laue class of the space group:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               n=1: -1                 n=8:  -3m (rhomb. axes)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               n=2: 2/m                n=9:  -31m (hex. axes)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               n=3: mmm                n=10: -3m1 (hex. axes)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               n=4: 4/m                n=11: 6/m '
  n = n + 1;  HELP_line(HELP_numor, n) = '                               n=5: 4/mmm              n=12: 6/mmm '
  n = n + 1;  HELP_line(HELP_numor, n) = '                               n=6: -3 (rhomb. axes)   n=13: m3'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               n=7: -3 (hex. axes)     n=14: m3m'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . "multi": output the general multiciplicity of the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                      space group'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             list space groups'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            LST_SG tetragonal centric:  list tetragonal and '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           centric space groups'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           LST_SG acentric monoclinic C: list monoclinic, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            acentric and C centred space groups '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           LST_SG laue_4: list tetragonal space groups with 4/m, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            Laue class'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           LST_SG trigonal enantio'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           LST_SG cubic chiral'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  LIST_SG, LST_SG, LIST_SPACE_GROUPS'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 48
  n = 0    ;  HELP_numor = HELP_MAG_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > MAG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          character string (atom or ion label)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             list magnetic features of the  current ion or atom: '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           electronic configuration, level, magnetic moment ...'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           for 3d and 4f elements'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  MAG, MAGNETIC'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            mag mn+3'
  n = n + 1;  HELP_line(HELP_numor, n) = '    .                      mag TB'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 49
  n = 0    ;  HELP_numor = HELP_MAN_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > MAN:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal arguments: character strings'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             get the CRYSCAL manual'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              list the meaning and use of CRYSCAL keywords'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if no arguments: all the keywords are listed'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . optionnal arguments are keywords name'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if "*" character in the character string: all the '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             keywords containing the characte string before'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             or after the "*" will be listed'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            MAN wave cell'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MAN ANG*'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  MAN, HELP'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 50
  n = 0    ;  HELP_numor = HELP_MAN_HTML_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > MAN_HTML:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = "    . optionnal argument:  'browse'"
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             get the CRYSCAL manual in HTML format'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           open the HTML file with the current browser'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              list the meaning and use of CRYSCAL keywords'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            MAN_HTML'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  MAN_HTML, HTML_MAN, HTML'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 51
  n = 0    ;  HELP_numor = HELP_MATMUL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > MATMUL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '      . type:                calculation keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '      . arguments:           no '
  n = n + 1;  HELP_line(HELP_numor, n) = '      . meaning:             input of components of 3*3 matrices as follows:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             M11, M12, M13, M21, M22, M23, M31, M32, M33'
  n = n + 1;  HELP_line(HELP_numor, n) = '      . ouput:               calculation of the M1xM2 matrix '
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n


 ! 52
  n = 0    ;  HELP_numor = HELP_MATR_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > MATR:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           9 reals '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             transformation (3,3) matrix components'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            MATR   0  0  1    0  1  0    -1  0  -1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MATR   0.5 0.5 0   -0.5 0.5 0   0 0 1'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MATR   1/2 1/2 0   -1/2 1/2 0   0 0 1'
  n = n + 1;  HELP_line(HELP_numor, n) = '  or'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           "#" 1 integer'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             number of matrix in the matrices list implemented in'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           CRYSCAL (see LIST_MATR keyword)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             MATR  #3'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . special arguments  :     "I"  : identity matrix'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              "-I" : inverse matrix'
  n = n + 1;  HELP_line(HELP_numor, n) = '  or'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           3 character strings: a, -a, b, -b, c, -c'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             setting of the unit cell'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             MATR  -c b a'
  n = n + 1;  HELP_line(HELP_numor, n) = '  or'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           2 character strings: P, A, B, C, I, F, R_rev, R_obv'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             Bravais lattices of the original and final cells.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           The transformation matrices are related to matrices'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           from #8 to #21 in the matrix list'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             MATR C P'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MATR R_obv R_rev'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              calculation of new cell parameters (dependent '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           keyword = CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculation of new atomic coordinates (dependent '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           keyword = ATOM)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculation of new hkl Miller indices (dependent'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           keyword = HKL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           create a new file containing the new hkl indices '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           (dependent keyword = FILE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  MAT, MATR, MATRIX'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 53
  n = 0    ;  HELP_numor = HELP_MENDEL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > MENDEL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTOUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           character strings'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             atom type'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              get atomic features: atomic number, weight, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            density, electronic configuration, oxydation'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            states, ionic radius, ...'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           get neutron data: scattering length, scattering'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            and absorption cross sections '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           get Xrays data: mass attenuation coefficient, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            total interaction cross sections, coefficient)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            for scattering factor calculation =f(stl),'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            values of D_fp and D_fpp anomalous dispersion'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            coefficient (for Cr, Fe, Cu, Mo and Ag radiations)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal keyword:   PLOT: create a .PGF file containing the scattering'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            factor values versus SinTheta/lambda and plot it'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            with WinPLOTR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            MENDEL TI C'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MENDEL Cu+2 PLOT'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 54
  n = 0    ;  HELP_numor = HELP_MERGE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > MERGE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  symmetry (character string): "TRIC", "MONO", "ORTHO", '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           "TETRA", "TRIG", "HEXA", "CUB"'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             merge the data in the current symmetry and create'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           a xxx_merge.HKL file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE, SPGR if no symmetry argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             merge monoclinic'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n


 ! 55 
  n = 0    ;  HELP_numor = HELP_MONOCLINIC_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > MONOCLINIC:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             . list the transformation matrices for equivalent '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             monoclinic settings'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if CELL exists: apply the transformation matrices'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             to give new monoclinic cell parameters'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional keyword:    CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  MONOCLINIC, MONOC, MONOCL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 56
  n = 0    ;  HELP_numor = HELP_NEWS_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > NEWS:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             list the last new facilities implemented in CRYSCAL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  specified year'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             NEWS 2010'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  NEWS'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 57
  n = 0    ;  HELP_numor = HELP_NIGGLI_CELL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > NIGGLI_CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              transform the current triclinic cell to the Niggli'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           cell.'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . related keyword:     CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  NIGGLi, NIGGLI_CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 58
  n = 0    ;  HELP_numor = HELP_OBV_REV_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > OBV_REV:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          1 integer'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             twin matrix number:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . 1: matrix law = (-1 0 0  0-1 0  0 0 1): twofold'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                 axis parallel to threefold '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . 2: matrix law = ( 0-1 0 -1 0 0  0 0-1): twofold'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                 axis parallel to (a-b)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             analysis of the reflections versus the parity'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . remark:              if no argument is given, first matrix law is taken'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  OBV_REV, OBVERSE_REVERSE, TWIN_OBVERSE_REVERSE'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           TWIN_OBV_REV, TWINNING_OBVERSE_REVERSE'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           TWINNING_OBV_REV'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 59
  n = 0    ;  HELP_numor = HELP_P4P_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ""
  n = n + 1;  HELP_line(HELP_numor, n) = "  > P4P:"
  n = n + 1;  HELP_line(HELP_numor, n) = ""
  n = n + 1;  HELP_line(HELP_numor, n) = "    . type:                OUTPUT keyword"
  n = n + 1;  HELP_line(HELP_numor, n) = "    . optional arguments : P4P file name"
  n = n + 1;  HELP_line(HELP_numor, n) = "    . meaning:             read the .P4P file (created by SAINT) present in the"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           current folder (if more than one .P4P file is present,"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           select a P4P file in a .P4P file list) and create"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           a P4P.cif file containing cell parameters (and esd's),"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           wavelength, number of reflexions used to refined the"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           cell parameters and related theta_min and theta_max."
  n = n + 1;  HELP_line(HELP_numor, n) = "    . identical keywords:  P4P, READ_P4P"
  n = n + 1;  HELP_line(HELP_numor, n) = ""
  HELP_lines_nb(HELP_numor) = n


 ! 60
  n = 0    ;  HELP_numor = HELP_PAUSE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ""
  n = n + 1;  HELP_line(HELP_numor, n) = "  > PAUSE:"
  n = n + 1;  HELP_line(HELP_numor, n) = ""
  n = n + 1;  HELP_line(HELP_numor, n) = "    . type:                "
  n = n + 1;  HELP_line(HELP_numor, n) = "    . arguments :          no"
  n = n + 1;  HELP_line(HELP_numor, n) = "    . meaning:             create a break in the execution of the requested"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           commands. This keyword can be useful when commands"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           are executed through a .CFL command file"
  n = n + 1;  HELP_line(HELP_numor, n) = ""
  HELP_lines_nb(HELP_numor) = n


 ! 61
  n = 0    ;  HELP_numor = HELP_PERMUT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > PERMUT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             apply matrix #3, #4, #5, #6 and #7 (see LIST_MATR '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           keyword) to the current cell parameters. This '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           corresponds to axes permutation '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             new cell parameters and new volume, obtained after'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           applying matrix #n'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           if cell parameters are unknown: list the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           transformation matrix'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 62 
  n = 0    ;  HELP_numor = HELP_Q_HKL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > QHKL: '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            real values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             Q=4pi*sinTheta/lambda values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             d_hkl(A), SinTheta/lambda(A-1)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           theta(deg) for known wavelength'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal keyword:   WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  QHKL, Q_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             QHKL 10.'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 63
  n = 0    ;  HELP_numor = HELP_QVEC_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > QVEC:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                input keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          3 reals'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             components of the propagation vector'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             apply the propagation vector on the hkl list and '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculate the corresponding d value'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  QVEC, Q_VEC, Q_VECTOR, KVEC, K_VEC, K_VECTOR,  '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           MOD_VEC, MODULATION_VECTOR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             QVEC 0.5 0.5 0.'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 64
  n = 0    ;  HELP_numor = HELP_READ_CEL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > READ_CEL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                input keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            1 character string corresponding to a CEL file name'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             read a CEL format file (PowderCELL format) to extract'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           structural features as: space group, cell parameters'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           and atomic positions.'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  READ_CEL, READ_CEL_FILE, READ_POWDERCELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             READ_CEL my_file.cel'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 65
  n = 0    ;  HELP_numor = HELP_READ_CIF_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > READ_CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                input keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            1 character string corresponding to a CIF file name'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             read a CIF format file to extract structural features:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           as: space group and/or symmetry operators, cell '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           parameters, atomic positions, wave ...'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  READ_CIF, READ_CIF_FILE, CIF_FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             READ_CIF my_file.cif'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 66
  n = 0    ;  HELP_numor = HELP_READ_INS_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > READ_INS'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                input keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            1 character string corresponding to a INS/RES file'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           name (SHELX)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             read a .INS/.RES SHELX format file to extract '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           structural features as: space group and/or symmetry'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           operators, cell parameters, atomic positions, ...'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  READ_INS, READ_INS_FILE, INS_FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             READ_INS my_file.ins'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 67
  n = 0    ;  HELP_numor = HELP_READ_PCR_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > READ_PCR'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                input keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            1 character string corresponding to a PCR file name'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             read a PCR FullProf format file to extract structural'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           features as: space group and/or symmetry operators,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           cell parameters, atomic positions, wave ...'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  READ_PCR, READ_PCR_FILE, PCR_FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             READ_PCR my_file.pcr'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 68
  n = 0    ;  HELP_numor = HELP_READ_NREPORT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = " "
  n = n + 1;  HELP_line(HELP_numor, n) = "  > READ_NREPORT:"
  n = n + 1;  HELP_line(HELP_numor, n) = " "
  n = n + 1;  HELP_line(HELP_numor, n) = "    . type:                input keyword"
  n = n + 1;  HELP_line(HELP_numor, n) = "    . argument:            no"
  n = n + 1;  HELP_line(HELP_numor, n) = "    . meaning:             read nreport.html file created by the Nonius software"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           suite for KCCD diffractometer and launch the browser defined"
  n = n + 1;  HELP_line(HELP_numor, n) = "                           in the 'cryscal.ini' setting file"
  n = n + 1;  HELP_line(HELP_numor, n) = "    . identical keywords:  READ_NREPORT, READ_HTML_REPORT, READ_NREPORTHTML"
  n = n + 1;  HELP_line(HELP_numor, n) = " "
  HELP_lines_nb(HELP_numor) = n


 ! 69
  n = 0    ;  HELP_numor = HELP_REC_ANG_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > REC_ANG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           2*3 reals'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             calculation of the angle between 2 vectors in the '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           reciprocal space. The 3 first real values are related'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           to the coordinates of the first vector and the 3 last'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           real values to the coordinates of the second vector'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  REC_ANG, RECANG, RECIPROCAL_ANGLE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             REC_ANG 1. 0. 0.   0. 1. 0.'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 70
  n = 0    ;  HELP_numor = HELP_REF_APEX_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = ' > REF_APEX'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write APEXII programs and device references'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  REF_APEX, REF_APEXII, WRITE_APEX, WRITE_APEXII,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           APEX, APEXII'
  HELP_lines_nb(HELP_numor) =  n

 ! 71
  n = 0    ;  HELP_numor = HELP_REF_EVAL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = ' > REF_DENZO'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write DENZO/SCALEPACK references'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  REF_DENZO, WRITE_DENZO, DENZO'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) =  n

 ! 72
  n = 0    ;  HELP_numor = HELP_REF_DENZO_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = ' > REF_EVAL'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write EVALCCD references'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  REF_EVAL, REF_EVALCCD, WRITE_EVAL, WRITE_EVALCCD,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           EVAL, EVALCCD'
  HELP_lines_nb(HELP_numor) =  n

 ! 73
  n = 0    ;  HELP_numor = HELP_REF_KCCD_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = ' > REF_KCCD'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write KCCD and device references'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  REF_KCCD, KCCD'
  HELP_lines_nb(HELP_numor) =  n

 ! 74
  n = 0    ;  HELP_numor = HELP_REF_SADABS_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = ' > REF_SADABS'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write SADABS references'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  REF_SAD, SADABS'
  HELP_lines_nb(HELP_numor) =  n




 ! 75
  n = 0    ;  HELP_numor = HELP_RESET_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = ' > RESET'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                input keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             initialization of all input parameters'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           and arrays'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  RESET, RAZ, INIT, INITIALIZATION'
  HELP_lines_nb(HELP_numor) =  n


 ! 76
  n = 0    ;  HELP_numor = HELP_RINT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > RINT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  symmetry (character string): "TRIC", "MONO", "ORTHO", '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           "TETRA", "TRIG", "HEXA", "CUB"'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             calculate the internal Rint value'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculate the completeness of the data collection'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE, SPGR if no symmetry argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  RINT, R_INT'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             RINT mono'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n


 ! 77
  n = 0    ;  HELP_numor = HELP_RHOMB_HEX_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > RHOMB_HEX:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             give transformation matrix from rhomboedral to'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           hexagonal setting'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  RHOMB_HEX, RHOMB_HEXA, RHOMB_TO_HEX, RHOMB_TO_HEXA'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) =  n

 ! 78
  n = 0    ;  HELP_numor = HELP_SEARCH_EXTI_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SEARCH_EXTI:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           2 optionnal real arguments: n_sig, ratio_criteria'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             analyse hkl reflections list and search systematic'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           extinctions'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             . n_sig: only reflections with I/sig > n_sig will'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               be analysed'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             . ratio_criteria: '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                if <F2_odd>/<F2_even> < ratio_criteria, the '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                current reflection type is considered as '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                extinction rule'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SEARCH_EXTI, FIND_EXTI'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             SEARCH_EXTI 2 0.02'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 79
  n = 0    ;  HELP_numor = HELP_SEARCH_SPGR_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SEARCH_SPACE_GROUP:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                calculation keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional arguments:  2 real : n_sig and threshold'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           1 character string, defining the crystal system:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            M/mono, O/ortho, T/tetra, R/trig, H/hexa, C/cub'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             search in a hkl file systematic extinctions and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           propose space groups in agreement with observed'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           extinctions. '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . n_sig: only reflections with I/sig > n_sig are '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                     considered. Default value is n_sig=3.'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . a threshold is also applied: only the '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            reflections with F2 > threshold * max(F2) are '
  n = n + 1;  HELP_line(HELP_numor, n) = '                            considered. Default value is threshold=0.03'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           Default value for crystal system is monoclinic'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SEARCH_SPGR, SEARCH_SPACE_GROUP, SEARCH_GROUP, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           CHECK_SPGR,  CHECK_SPACE_GROUP,  CHECK_GROUP'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n




 ! 80
  n = 0    ;  HELP_numor = HELP_SET_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SET:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:        INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:   2 string arguments'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:     definition of external applications:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                    . if arg(1) = "EDITOR":  arg(2) is the browser application'
  n = n + 1;  HELP_line(HELP_numor, n) = '                    . if arg(1) = "BROWSER": arg(2) is the editor application'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:     set browser "C:\Program Files\Mozilla Firefox\firefox.exe"'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 81
  n = 0    ;  HELP_numor = HELP_SETTING_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SETTING:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:        OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:   no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:     output the different parameters defined in the CRYSCAL.INI'
  n = n + 1;  HELP_line(HELP_numor, n) = '                   setting file'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 82 
  n = 0    ;  HELP_numor = HELP_SFAC_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SFAC:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           n character strings '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning.:            list of chemical elements in the molecule'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              - atomic density calculation, in at/cm3 '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keyword = CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - absorption coefficient calculation '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywords = CELL, BEAM)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - determination of the molecular formula '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywords = Z)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - calculation of the molecular weight '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywords = Z)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - calculation of the atomic and weight '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             percentage (dependent keywords = Z)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - determination of the density '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywords = CELL, Z)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   UNIT'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             SFAC C O H N'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  ! 83
  n = 0    ;  HELP_numor = HELP_SFHKL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SF_HKL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           3 reals '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             h,k,l, Miller indices'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              calculation of crystallographic structure'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           factor for the hkl reflection'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           (dependent keyword = CELL, WAVE, SPGR, ATOM)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             HKL 1 0 1'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   CELL, WAVE, SPGR, ATOM'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . related keyword :    BEAM'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SF_HKL, SFAC_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n



 ! 84
  n = 0    ;  HELP_numor = HELP_SG_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          character string'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             space group (Hall Mauguin symbol, number in IT'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             . informations on the space group: list of symmetry'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             operators and symmetry symbols, Wyckoff positions,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             extinctions  ...   '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . get informations on the Wyckoff site: multiplicity,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             list of symmetry operators and symmetry elements of'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             the site point group'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if HKL file has been previously read (see FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             keyword), the number of reflections in agreement'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             with the space group is output'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . see:                 SG_INFO, SG_EXTI, LST_SYM_OP, SITE_INFO'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             SPGR P 21 21 21 '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SPGR, SG, SPACE_GROUP'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 85
  n = 0    ;  HELP_numor = HELP_SG_ALL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SG_ALL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             get informations on the current space group:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            symmetry operator (complete and reduced set), Wyckoff'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            positions, extinction rules ...'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SPGR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SG_ALL, SP_ALL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) =  n

 ! 86
  n = 0    ;  HELP_numor = HELP_SG_EXTI_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SG_EXTI:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             list of extinctions rules of the current space'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           group'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SPGR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SP_EXTI, SP_EXTI, SG_EXTINCTIONS, SPACE_GROUP_EXTI,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           SPACE_GROUP_EXTINCTIONS'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 87
  n = 0    ;  HELP_numor = HELP_SG_INFO_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SG_INFO:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             informations on the space group: list of'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           symmetry operators and symmetry symbols, Wyckoff'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           positions'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SPGR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal argument:  ALL: list all symmetry operators (including '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                inversion and lattice centring translations)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SG_INFO, SP_INFO, SPACE_GROUP_INFO, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           LIST_SPACE_GROUP_INFO'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 88
  n = 0    ;  HELP_numor = HELP_SG_SUB_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SG_SUB:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             informations on the space group: list of'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           Translationengleische Subgroups for a given'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           space group'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SPGR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SG_SUB, SG_SUBGROUP'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 89
  n = 0    ;  HELP_numor = HELP_SHANNON_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SHANNON:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            character string'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             get effective ionic radii from the Shannon table'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           (Acta Cryst 1976, A32, 751)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . ec:     electronic configuration'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . CN:     coordinence'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . SP:     configuration de spin'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . CR:     crystal radius'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            . IR:     effective radius'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SHAN, SHANNON'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            SHANNON Pb'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           SHANNON CU+2'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n



 ! 90
  n = 0    ;  HELP_numor = HELP_SHELL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SHELL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           2 character strings (one mandatory argument) + 2 reals'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             arg_1 = d:     keeps reflections in the d_min and '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          d_max range'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_shell_d.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          dependent keyword: CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           arg_1 = stl:   keeps reflections in the stl_min and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          stl_max range'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_shell_stl.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          dependent keyword: CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           arg_1 = theta: keeps reflections in the theta_min and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          theta_max range'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_shell_theta.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          dependent keyword: CELL, WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           arg_1 = int:   keeps reflections in the int_min and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          int_max range'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_shell_i.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           arg_1 = isig:  keeps reflections in the i/sig_min and '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          i/sig_max range'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_shell_isig.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           arg_2 = plot:  create a .PGF file from the created '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          HKL file and plot it with WinPLOTR '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          (excepted for arg_1=int and '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          arg_1 = isig)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           real_1:        X_min value'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           real_2:        X_max value'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional argument:   arg_2 = plot'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            SHELL d plot'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           SHELL theta'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 91
  n = 0    ;  HELP_numor = HELP_SHIFT_2TH_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SHIFT_2TH:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          1 real'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             diffractometer 2theta shift'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             apply the diffractometer 2theta shift to the '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculated 2theta value:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            2theta = 2theta_calc + shift_2theta'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   HKL, GEN_HKL, GEN_SAT'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SHIFT_2TH, SHIFT_2THETA, 2TH_SHIFT, 2THETA_SHIFT'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 92
  n = 0    ;  HELP_numor = HELP_SITE_INFO_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SITE_INFO:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional arguments : list of atomic labels (character strings)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             give information on Wyckoff atomic positions and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           apply the symmetry operators of the current space'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           group on the atomic positions'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            if no argument:    all atoms of the atoms list '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                           (cf ATOM keyword) are considered'
  n = n + 1;  HELP_line(HELP_numor, n) = '                            if arg(1) = "ALL": all atoms of the atoms list'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          (cf ATOM keyword) are considered'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SPGR, ATOM'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             SITE_INFO O1 C8'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SITE_INFO, LIST_SITE_INFO'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 93
  n = 0    ;  HELP_numor = HELP_SIZE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SIZE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           3 reals '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             crystal dimensions in mm'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              crystal volume calculation'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculation of the radius of a sphere with identical'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           volume'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             SIZE 0.11 0.13 0.122'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SIZE, CRYSTAL_SIZE'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 94
  n = 0    ;  HELP_numor = HELP_SORT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SORT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           2 character strings '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             arg_1 = d:     sort HKL file in increasing d_hkl '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          (dependent keyword: CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_sort_d.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           arg_1 = stl:   sort HKL file in increasing '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          sinTheta/lambda'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          (dependent keyword: CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_sort_stl.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           arg_1 = theta: sort HKL file in increasing Theta '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          (dependent keyword: CELL, WAVE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_sort_theta.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '    .                      arg_1 = int:   sort HKL file in decreasing F2'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL_sort_i.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '    .                      arg_1 = isig:  sort HKL file in decreasing F2/sigma'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          created hkl file: HKL__sort_isig.hkl'
  n = n + 1;  HELP_line(HELP_numor, n) = '    .                      arg_2 = plot:  create a .PGF file from the created'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          HKL file and plot it with WinPLOTR '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                         (excepted for arg_1=int and '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                          arg_1 = isig)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           arg_2 = out_n: list the n first sorted reflections'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   FILE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional argument:   arg_2 = plot'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            SORT d plot'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           SORT stl'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           SORT theta OUT_10'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 95
  n = 0    ;  HELP_numor = HELP_STL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > STL: '
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            real values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             sinTheta / lambda values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             d_hkl(A), Q(A-1)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           theta(deg) for known wavelength'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optionnal keyword:   WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  STL, STL_HKL, SINTHETA/WAVE, SINTHETA/LAMBDA'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             STL 0.70'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor ) = n

 ! 96
  n = 0    ;  HELP_numor = HELP_SYMM_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SYMM:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          12 reals'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             symmetry operator in numeric form:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             S11 S12 S13  T1   S21 S22 S23 T2   S31 S32 S33 T3'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               Sij: components of the rotational part of the '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                    symmetry operator'
  n = n + 1;  HELP_line(HELP_numor, n) = '                               Ti:  components of the translationnal part of '
  n = n + 1;  HELP_line(HELP_numor, n) = '                                    the symmetry operator'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             SYMM   1  0  0   0.5    0 -1  0   0.5   0  0 -1  0. '
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '   or'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          string'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             symmetry operator in alphanumeric form'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             SYMM   x+1/2, -y+1/2, -z'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             apply the current symmetry operator to atomic '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           coordinates (dependent keyword = ATOM)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional keyword:    APPLY_OP, SYM_OP'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SYMM, SYM, SYMMETRY_OPERATOR'
  n = n + 1;  HELP_line(HELP_numor, n) = '  '
  HELP_lines_nb(HELP_numor) = n

 ! 97
  n = 0    ;  HELP_numor = HELP_SYST_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > SYST:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                EXTERNAL keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           1 character string'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             launch an external command or program'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             SYST dir *.CFL / P'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  SYST, CMD, COMMAND, DOS, DOS_COMMAND'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) =  n

 ! 98
  n = 0    ;  HELP_numor = HELP_THERM_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > THERM:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT/OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           1 character string + reals'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             conversion of atomic displacement parameters:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . available arguments: '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "Uiso", "U_iso": conversion from Uiso to Biso'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "Biso", "B_iso": conversion from Biso to Uiso'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "U_ij", "Uij"  : conversion from U_ij to B_ij and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                                beta_ij'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "U_ij", "Uij"  : conversion from U_ij to B_ij and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                                beta_ij'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "beta_ij", "betaij", "beta"": conversion from '
  n = n + 1;  HELP_line(HELP_numor, n) = '                               beta_ij to B_ij and U_ij'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . reals values correspond to isotropic thermal'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             parameters values (Biso or Uiso) or anisotropic'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             thermal parameters (Uij, Bij, Beta) in the following'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             order: 11, 22, 33, 12, 13, 23'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            THERM Biso 0.52 0.76 0.35'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           THERM B_ij 0.01 0.01 0.01 0.0 0.0 0.0'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword    CELL for anisotropic thermal parameters conversion'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  THERM, THERMAL, ADP'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  
   ! 99
  n = 0    ;  HELP_numor = HELP_THERM_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > THERM_SHELX:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT/OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           1 character string + reals'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             conversion of atomic displacement parameters:'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . available arguments: '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "Uiso", "U_iso": conversion from Uiso to Biso'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "Biso", "B_iso": conversion from Biso to Uiso'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "U_ij", "Uij"  : conversion from U_ij to B_ij and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                                beta_ij'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "U_ij", "Uij"  : conversion from U_ij to B_ij and'
  n = n + 1;  HELP_line(HELP_numor, n) = '                                                beta_ij'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             - "beta_ij", "betaij", "beta"": conversion from '
  n = n + 1;  HELP_line(HELP_numor, n) = '                               beta_ij to B_ij and U_ij'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . reals values correspond to isotropic thermal'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             parameters values (Biso or Uiso) or anisotropic'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             thermal parameters (Uij, Bij, Beta) in the following'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             order: 11, 22, 33, 23, 13, 12 (as in SHELXL program)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            THERM Biso 0.52 0.76 0.35'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           THERM B_ij 0.01 0.01 0.01 0.0 0.0 0.0'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword    CELL for anisotropic thermal parameters conversion'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  THERM_SHELX, THERMAL_SHELX, ADP_SHELX'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 100
  n = 0    ;  HELP_numor = HELP_THETA_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > THETA: '
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            real values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             Theta (deg) values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             d_hkl(A), Q(A-1), SinTheta/lambda'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  THETA, THETA_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             THETA 27.5'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 101
  n = 0    ;  HELP_numor = HELP_TITL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > TITL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           character strings'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             title of the job'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             Ammonium bitartrate    '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  TITL, TITLE'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 102
  n = 0    ;  HELP_numor = HELP_TRANSLATION_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > TRANSLATION:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          translation components (3 reals)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             apply the translation on atomic coordinates'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   ATOM'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             TRANSLATION 0.5 0.5 0.5'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  TRANSLATION, TRANS, TRANSLATE, MOVE'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 103
  n = 0    ;  HELP_numor = HELP_TRANSMISSION_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > TRANSMISSION:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT/OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          list of distances (in mm)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             calculation of transmission coefficient for'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           different distances given as arguments'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   CELL, UNIT, SFAC'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 104
  n = 0    ;  HELP_numor = HELP_TWIN_HEXA_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > TWIN_HEXA:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             . list the transformation matrices for hexagonal setting'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if CELL exists: apply the transformation matrices'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             to give new cell parameters'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional keyword:    CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  HEXA_TWIN, HEXA_TWINNING, HEXAGONAL_TWIN, HEXAGONAL_TWINNING'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           TWIN_HEXA, TWIN_HEXAGONAL, TWINNING_HEXA, TWINNING_HEXAGONAL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

! 105
  n = 0    ;  HELP_numor = HELP_TRICLINIC_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > TRICLINIC:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             . list the transformation matrices for equivalent '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             triclinic unit cells'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if CELL exists: apply the transformation matrices'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             to give new triclinic cell parameters'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional keyword:    CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  TRICLINIC, TRICL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n



 ! 106 
  n = 0    ;  HELP_numor = HELP_TWIN_PSEUDO_HEXA_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > TWIN_PSEUDO_HEXA:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             . list the transformation matrices for pseudo-hexagonal'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              setting in a monoclinic unit cell'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . if CELL exists: apply the transformation matrices'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             to give new cell parameters'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional keyword:    CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  TWIN_PSEUDO_HEXA'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 107
  n = 0    ;  HELP_numor = HELP_TWO_THETA_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > TWO_THETA: '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                CALCULATION keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            real values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             2Theta (deg) values'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             d_hkl(A), Q(A-1), SinTheta/lambda'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  TWO_THETA, 2THETA, 2THETA_HKL, TWO_THETA_HKL'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             TWO_THETA 50.'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 108
  n = 0    ;  HELP_numor = HELP_UNIT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > UNIT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           n character strings '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             list of chemical elements in the molecule'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              - atomic density calculation, in at/cm3 '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keyword = CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - absorption coefficient calculation'
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywords = CELL, BEAM)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - determination of the molecular formula '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywords = Z)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - calculation of the molecular weight '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywords = Z)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - calculation of the atomic and weight '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             percentage (dependent keywords = Z)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           - determination of the density '
  n = n + 1;  HELP_line(HELP_numor, n) = '                             (dependent keywords = CELL, Z)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SFAC'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             UNIT 16. 24. 36. 4.'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  ! 109
  n = 0    ;  HELP_numor = HELP_MATR_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > USER_MAT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           "#" + 1 integer'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             number of user matrix in the matrices list defined  by'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           the user in the cryscal.ini setting file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             USER_MAT  #3'
  n = n + 1;  HELP_line(HELP_numor, n) = '  or'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           "$" character string'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             this character string has to be one of the comment text'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           accompagnying the matrix defined by the user in the'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           cryscal.ini setting file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             MATR  $2a'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . output:              calculation of new cell parameters (dependent '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           keyword = CELL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculation of new atomic coordinates (dependent '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           keyword = ATOM)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           calculation of new hkl Miller indices (dependent'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           keyword = HKL)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           create a new file containing the new hkl indices '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           (dependent keyword = FILE)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  USER_MAT, USER_MATR, USER_MATRIX'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

  
  
 ! 110
  n = 0    ;  HELP_numor = HELP_WAVE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WAVE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           1 real or 1 character string'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             . real value: wavelength value in A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           . character string: '
  n = n + 1;  HELP_line(HELP_numor, n) = '                              - "X_Ag" or "XAG": wavelength = 0.556363 A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              - "X_Mo" or "XMO": wavelength = 0.71073 A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              - "X_Cu" or "XCU": wavelength = 1.5406 A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              - "X_Ni" or "XNI": wavelength = 1.65794 A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              - "X_Co" or "XCo": wavelength = 1.78892 A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              - "X_Fe" or "XFE": wavelength = 1.93597 A'
  n = n + 1;  HELP_line(HELP_numor, n) = '                              - "X_Cr" or "XCR": wavelength = 2.28962 A'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . examples:            WAVE 0.71073'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           WAVE XMO'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WAVE, WAVELENGTH, WL'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n


 ! 111
  n = 0    ;  HELP_numor = HELP_WEB_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WEB:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            web site name or address'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             web site name can be one of the following: CDIFX'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           or site name defined in the cryscal.ini setting file'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           in the [WEB] section'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory argument:  CELL or FILE .CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WEB, INTERNET'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) =  n

 ! 112
  n = 0    ;  HELP_numor = HELP_WRITE_BEAM_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WRITE_BEAM:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write incident beam type)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword :  WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WRITE_BEAM, WRITE_INCIDENT_BEAM, OUTPUT_BEAM)'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           OUTPUT_INCIDENT_BEAM)'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n


 ! 113
  n = 0    ;  HELP_numor = HELP_WRITE_CELL_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WRITE_CELL:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write cell parameters (direct and reciprocal space)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword :  CELL or FILE .CIF'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WRITE_CELL, OUTPUT_CELL'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 114
  n = 0    ;  HELP_numor = HELP_WRITE_CHEM_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n)  = '  > WRITE_CHEM:'
  n = n + 1;  HELP_line(HELP_numor, n)  = ''
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . meaning:             write molecular features (formula, weight )'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . mandatory keywords:  CHEM, CONT/ZUNIT, SFAC/UNIT/ZUNIT'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    . identical keywords:  WRITE_CHEM, WRITE_CHEMICAL_FORMULA,'
  n = n + 1;  HELP_line(HELP_numor, n)  = '    .                      OUTPUT_CHEM'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 115
  n = 0    ;  HELP_numor = HELP_WRITE_DEVICE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WRITE_DEVICE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write device features defined in the cryscal.ini setting file'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WRITE_DEVICE, OUTPUT_DEVICE'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) =  n


  ! 116
  n = 0    ;  HELP_numor = HELP_WRITE_QVEC_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WRITE_QVEC:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write modulation wave vector components'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword :  QVEC'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WRITE_QVEC, OUTPUT_QVEC'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

 ! 117
  n = 0    ;  HELP_numor = HELP_WRITE_SG_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WRITE_SG:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write space group features'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword :  SPGR'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WRITE_SG, WRITE_SPGR, WRITE_SPACE_GROUP'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  HELP_lines_nb(HELP_numor) = n

! 118
  n = 0    ;  HELP_numor = HELP_WRITE_SYM_OP_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WRITE_SYM_OP:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments :          no argument'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . outputs:             list of symmetry operators (alphanumeric form, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           rotational and translational parts with the SYMM'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword:   SPGR or SYMM'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WRITE_SYM_OP, WRITE_SYMM_OP, WRITE_SYMP, WRITE_SYMM,'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           WRITE_SYMMETRY_OPERATORS'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 119
  n = 0    ;  HELP_numor = HELP_WRITE_WAVE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '  > WRITE_WAVE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                OUTPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . argument:            no'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             write cell parameters (direct and reciprocal space)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . mandatory keyword :  WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  WRITE_WAVE, OUTPUT_WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = ''
  n = n + 1;  HELP_lines_nb(HELP_numor) =  9

 ! 120
  n = 0    ;  HELP_numor = HELP_X_WAVE_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > X_WAVE:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                output keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . optional argument:   target nature (Ag, Mo, Cu, Ni, Co, Fe, Cr)'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             if no optional argument, output Ka1 and Ka2 wavelength'
  n = n + 1;  HELP_line(HELP_numor, n) = '                           value (in A) for main radiations: Ag, Mo, Cu, Ni, Co, '
  n = n + 1;  HELP_line(HELP_numor, n) = '                           Fe and Cr'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  XRAYS_WAVELENGTH, X_WAVE'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

 ! 121
  n = 0    ;  HELP_numor = HELP_ZUNIT_numor
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '  > ZUNIT:'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  n = n + 1;  HELP_line(HELP_numor, n) = '    . type:                INPUT keyword'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . arguments:           1 real'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . meaning:             number of formula units'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . example:             4'
  n = n + 1;  HELP_line(HELP_numor, n) = '    . identical keywords:  ZUNIT, Z, Z_UNIT'
  n = n + 1;  HELP_line(HELP_numor, n) = ' '
  HELP_lines_nb(HELP_numor) = n

end subroutine def_keywords_lines

!-------------------------------------------------------------------------------------------------------------
subroutine Def_command_line_arguments
 use Text_module, only : CLA_line, CLA_lines_nb
 implicit none
  integer   :: CLA_numor, n

 CLA_numor = 0


 n = 0; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)  = "   > MAN  : "
 n = n + 1; CLA_line(CLA_numor, n)  = "       Create 'cryscal.txt' file user's guide (text format)."
 n = n + 1; CLA_line(CLA_numor, n)  = "       ex : > cryscal man"
 CLA_lines_nb(CLA_numor) = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)  = "   > HTML  : "
 n = n + 1; CLA_line(CLA_numor, n)  = "       Create 'cryscal.html' file user's guide (HTML format)."
 n = n + 1; CLA_line(CLA_numor, n)  = "       ex : > cryscal HTML"
 CLA_lines_nb(CLA_numor) = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)  = "   > HTML_BROWSE  : "
 n = n + 1; CLA_line(CLA_numor, n)  = "       Create 'cryscal.html' file user's guide (HTML format) and launch the browser"
 n = n + 1; CLA_line(CLA_numor, n)  = "       previously defined in the 'cryscal.ini' setting file'"
 n = n + 1; CLA_line(CLA_numor, n)  = "       ex : > cryscal HTML_BROWSE"
 CLA_lines_nb(CLA_numor) = n



 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > KEY/KEYS : "
 n = n + 1; CLA_line(CLA_numor, n)   = "       Create 'cryscal_keys.txt' file containing the keywords list"
 n = n + 1; CLA_line(CLA_numor, n)   = "        ex : > cryscal KEYS"
 CLA_lines_nb(CLA_numor)  = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > P4P  : "
 n = n + 1; CLA_line(CLA_numor, n)   = "       . List the .P4P files in the current folder."
 n = n + 1; CLA_line(CLA_numor, n)   = "       . Select a .P4P file"
 n = n + 1; CLA_line(CLA_numor, n)   = "       . Search for the corresponding .HKL file and output SADABS file"
 n = n + 1; CLA_line(CLA_numor, n)   = "         to create an 'import.cif' file."
 n = n + 1; CLA_line(CLA_numor, n)   = "        ex : > cryscal P4P"
 CLA_lines_nb(CLA_numor)  = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > filename.P4P : "
 n = n + 1; CLa_line(CLA_numor, n)   = "        Search for the corresponding .HKL file and output SADABS file"
 n = n + 1; CLA_line(CLA_numor, n)   = "        to create an 'import.cif' file."
 n = n + 1; CLA_line(CLA_numor, n)   = "         ex : > cryscal mycrystal.P4P"
 CLA_lines_nb(CLA_numor)  = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > filename.RAW : "
 n = n + 1; CLA_line(CLA_numor, n)   = "        Read the filename.RAW file (created by SAINT) and create"
 n = n + 1; CLA_line(CLA_numor, n)   = "        a HKL file with the SHELX format (3I4,2F8.2)."
 n = n + 1; CLA_line(CLA_numor, n)   = "         ex : > cryscal mycryscal.RAW"
 CLA_lines_nb(CLA_numor)  = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > REPORT/REPORT_LONG : "
 n = n + 1; CLA_line(CLA_numor, n)   = "        Create a HTML report from 'archive.cif' file or from the"
 n = n + 1; CLA_line(CLA_numor, n)   = "        .CIF file give as a second argument."
 n = n + 1; CLA_line(CLA_numor, n)   = "         ex : > cryscal report"
 CLA_lines_nb(CLA_numor)  = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > ARCHIVE.CIF : "
 n = n + 1; CLA_line(CLA_numor, n)   = "       . Read 'archive.cif' CIF file"
 n = n + 1; CLA_line(CLA_numor, n)   = "       . Complete this 'archive.cif' with additionnal CIF fields related to"
 n = n + 1; CLA_line(CLA_numor, n)   = "         absorption correction, squeeze procedure, hydrogen treatment, "
 n = n + 1; CLA_line(CLA_numor, n)   = "         diffractometer, structure solution program ... "
 n = n + 1; CLA_line(CLa_numor, n)   = "         A new 'cryscal_archive.cif' file is then created."
 n = n + 1; CLA_line(CLA_numor, n)   = "       Remark : diffractometer and structure solution program can be defined"
 n = n + 1; CLA_line(CLA_numor, n)   = "                in the 'cryscal.ini' setting file in the [DEVICE] and [PROGRAMS]"
 n = n + 1; CLA_line(CLA_numor, n)   = "                parts respectively."
 CLA_lines_nb(CLA_numor)  = n

 ! n = 0 ; CLA_numor = CLA_numor + 1
 ! n = n + 1; CLA_line(CLA_numor, n)   = "   > ACTA : "
 ! n = n + 1; CLA_line(CLA_numor, n)   = "      Combined with 'ARCHIVE.CIF' argument, this optional ACTA argument"
 ! n = n + 1; CLA_line(CLA_numor, n)   = "      completes the 'cryscal_archive.cif' file with several CIF fields related"
 ! n = n + 1; CLA_line(CLA_numor, n)   = "      to CIF file deposition."
 ! n = n + 1; CLA_line(CLA_numor, n)   = "      Can be combined with 'CIFDEP' argument."
 ! CLA_lines_nb(CLA_numor)  = n

  n = 0; CLA_numor = CLA_numor + 1
  n = n + 1; CLA_line(CLA_numor, n)   = "   > CIFDEP : "
  n = n + 1; CLA_line(CLA_numor, n)   = "       Combined with 'ARCHIVE.CIF' argument, this optional argument"
  n = n + 1; CLA_line(CLA_numor, n)   = "       completes the 'cryscal_archive.cif' file with CIF fields related"
  n = n + 1; CLA_line(CLA_numor, n)   = "       to the author of the deposition (name, address), extracted from"
  n = n + 1; CLA_line(CLA_numor, n)   = "       the 'cryscal.ini' setting file ([AUTHOR] part)"
 ! n = n + 1; CLa_line(CLA_numor, n)   = "       Can be combined with 'ACTA' argument."
  CLA_lines_nb(CLA_numor)  = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > SOLVE_TO_INS/CREATE_INS : "
 n = n + 1; CLA_line(CLA_numor, n)   = "       . read STRUCT.CIF file and get cell parameters with esd's"
 n = n + 1; CLA_line(CLA_numor, n)   = "       . read import.RES created by SIRxx or SHELXS"
 n = n + 1; CLA_line(CLA_numor, n)   = "       . create job.INS file for SHELXL with correct esd's and different"
 n = n + 1; CLa_line(CLA_numor, n)   = "         useful SHELXL keywords (ACTA, BOND$H ...)"
 CLA_lines_nb(CLA_numor)  = n




 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > DEBUG: "
 n = n + 1; CLA_line(CLA_numor, n)   = "       Combined with all previous arguments, this optional argument"
 n = n + 1; CLA_line(CLA_numor, n)   = "       will create a 'cryscal_debug.txt' file containing informations"
 n = n + 1; CLA_line(CLA_numor, n)   = "       about the values of some variables during the CRYSCAL run."
 n = n + 1; CLA_line(CLA_numor, n)   = "       This argument can be useful to detect the origin of the bug"
 n = n + 1; CLA_line(CLA_numor, n)   = "       in a CRYSCAL crash."
 CLA_lines_nb(CLA_numor)  = n

 n = 0 ; CLA_numor = CLA_numor + 1
 n = n + 1; CLA_line(CLA_numor, n)   = "   > NO_OUT : "
 n = n + 1; CLA_line(CLA_numor, n)   = "       Combined with all previous arguments, this optional argument"
 n = n + 1; CLA_line(CLA_numor, n)   = "       avoids to write CRYSCAL results lines in the screen and"
 n = n + 1; CLA_line(CLA_numor, n)   = "       in the 'cryscal.log' file."
 CLA_lines_nb(CLA_numor)  = n



end subroutine Def_command_line_arguments
!-------------------------------------------------------------------------------------------------------------
