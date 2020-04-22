!     Last change:  TR   11 Oct 2007    6:11 pm
!     Last change:  TR   11 Oct 2007    6:11 pm
subroutine cryscalc_init()
 USE cryscalc_module
 USE pattern_profile_module
 USE wavelength_module
 USE atome_module
 USE HKL_module
 USE CIF_module
 USE MATRIX_list_module, only : user_mat_nb


  implicit none
   INTEGER   :: numor, i

 ! initialisation ---------------------------------------
  winplotr_path_name  = ''
  winplotr_exe        = ''

  CRYSCALC%version      = 'March 2017'
  CRYSCALC%date         = '?'
  CRYSCALC%author       = 'Thierry Roisnel (CDIFX / ISCR UMR6226 - Rennes)'
  CRYSCALC%mail         = "thierry.roisnel@univ-rennes1.fr"
  CRYSCALC%url          = "www.cdifx.univ-rennes1.fr/cryscalc"
  CRYSCALC%ini          = ''
  CRYSCALC%css          = ''
  CRYSCALC%report_css   = ''
  CRYSCALC%path_name    = ''
  CRYSCALC%url_exe      = 'http://www.cdifx.univ-rennes1.fr/progs/cryscalc/cryscalc.exe'
  CRYSCALC%compiler     = "?"
  CRYSCALC%option       = "?"

  my_editor%name    = '?'
  my_browser%name   = '?'
  my_word%name      = '?'
  my_pdflatex%name  = '?'

  WEB%NAME          = '?'
  WEB%ADDRESS       = '?'


! initialisation des champs lus dans le fichier de config. CRYSCALC.INI
  my_editor%exist   = .false.
  my_browser%exist  = .false.
  my_word%exist     = .false.
  my_pdflatex%exist = .false.

  AUTHOR%string      = '?'
  AUTHOR%name        = '?'
  AUTHOR%first_name  = '?'
  AUTHOR%initiales   = '?'
  AUTHOR%ADDRESS     = '?'
  AUTHOR%email       = '?'
  AUTHOR%WEB         = '?'
  AUTHOR%TEAM        = '?'

  DEVICE%string    = "?"
  DEVICE%diffracto = "?"
  DEVICE%lab       = "?"
  DEVICE%radiation = "?"
  DEVICE%wave      = "?"
  DEVICE%APEX      = .false.
  DEVICE%D8V       = .false.
  DEVICE%D8VMo     = .false.
  DEVICE%D8VCu     = .false.
  DEVICE%KCCD      = .false.
  DEVICE%X2S       = .false.
  DEVICE%Xcalibur  = .false.
  DEVICE%SuperNova = .false.


  !Structure_solution%title       = "SIR97: a new tool for crystal structure determination and refinement"
  !Structure_solution%DOI         = "10.1107/S0021889898007717"
  Structure_solution%name         = 'SIR97'
  Structure_solution%reference    = 'A. Altomare, M. C. Burla, M. Camalli, G. Cascarano, '// &
                                    'C. Giacovazzo, A. Guagliardi, A. G. G. Moliterni, '//   &
                                    'G. Polidori, R. Spagna, J. Appl. Cryst. (1999) 32, 115-119'
  Structure_solution%CIF_ref      = 'SIR97 (Altomare et al., 1999)'


  !Structure_solution%title        = "SIR2004: an improved tool for crystal structure determination and refinement"
  !Structure_solution%DOI          = "10.1107/S002188980403225X"
  !Structure_solution%name         = 'SIR2004'
  !Structure_solution%reference    = 'M.C. Burla, R. Caliandro, M. Camalli, B. Carrozzini, '//  &
  !                                   'G. Cascarano, L. De Caro, C. Giacovazzo, G. Polidori, '// &
  !  								  'R. Spagna, J. Appl. Cryst. (2005) 38, 381-388'
  !Structure_solution%CIF_ref      = 'SIR2004 (Burla., 2005)  '

  !Structure_solution%title        = 'SHELXT - Integrated space-goup and crystal-structure determination'
  !Structure_solution%DOI          = '10.1107/S2053273314026370'
  !Structure_solution%name         = 'SHELXT'
  !Structure_solution%reference    = 'Sheldrick G.M., Acta Cryst. A71 (2015) 3-8'
  !Structure_solution%CIF_ref      = 'SHELXT (Sheldrick G.M., 2015)'

  !Structure_solution%title        = 'Superflip - a computer program for the solution of crystal structures by charge flipping in arbitrary dimensions'
  !Structure_solution%DOI          = '10.1107/S0021889807029238'
  !Structure_solution%name         = 'Superflip'
  !Structure_solution%reference    = 'Palatinus L., Chapuis G., J. Appl. Cryst. (2007, 40, 786-790'
  !Structure_solution%CIF_ref      = 'SUPERFLIP (Palatinus et al, 2007)'


  CIF_parameter%atom_sites_solution_1 = "direct"    ! "structure-invariant direct methods"
  CIF_parameter%atom_sites_solution_2 = "difmap"    ! "difference Fourier map"
  !Structure_refinement%title     = 'A short history of SHELX'
  !Structure_refinement%DOI       = '10.1107/S0108767307043930'
  Structure_refinement%name       = 'SHELXL-97'
  Structure_refinement%reference  = 'Sheldrick G.M., Acta Cryst. A64 (2008), 112-122'
  Structure_refinement%CIF_ref    = 'SHELXL-97 (Sheldrick, 2008)'

  !Structure_refinement%title     = 'Crystal structure refinement with SHELXL'
  !Structure_refinement%DOI       = '10.1107/S2053229614024218'
  Structure_refinement%name       = 'SHELXL-2014'
  Structure_refinement%reference  = 'SHELXL-2014/6: Sheldrick G.M., Acta Cryst. C71 (2015) 3-8'
  Structure_refinement%CIF_ref    = 'SHELXL-2014/6 (Sheldrick, 2015)'

  Absorption_correction%name      = '?'
  Absorption_correction%reference = '?'
  Absorption_correction%CIF_ref   = '?'
  !SHELXL_2014                     = .false.

!-------------------------------------------------------------------------------

  debug_proc%unit    =  10
  debug_proc%write   = .false.
  debug_proc%level_1 = .false.
  debug_proc%level_2 = .false.
  debug_proc%level_3 = .false.


  ON_SCREEN          = .true.
  ON_SCREEN_PRF      = .false.
  ON_screen_CIF_item = .false.
  CIFdep             = .false.
  ACTA               = .false.
  step_by_step       = .true.

  ACE_file_name            = ''
  CEL_file_name            = ''
  CIF_file_name            = ''
  CIF_pymol_file_name      = ''
  input_CIF_file           = ''
  CIF_file_exist           = .false.
  nb_input_CIF_files       = 0
  archive_CIF              = ''
  INS_file_name            = ''
  PCR_file_name            = ''
  P4P_file_name            = ''
  M50_file_name            = ''
  X_file_name              = ''
  RED_file_name            = ''
  TIDY_out_file_name       = ''
  RMAT_file_name           = ''
  RAW_file_name            = ''
  HKL_file_name            = ''
  FCF_file_name            = ''
  FCF_plot                 = .false.
  FCF_plot_stl             = .false.
  HKL_file_diff            = ''
  ABS_file_name            = ''
  FACES_file_name          = ''
  main_title               = ''
  SAINT%version            = ''
  SAINT%date               = ''
  SAINT%source             = ''
  Experiment%date          = '?'
  !Experiment%DX            = '?'
  Experiment%temp          = '?'
  Experiment%target        = '?'
  Experiment%voltage       = '?'
  Experiment%current       = '?'
  Experiment%detector      = '?'
  Experiment%detector_temp = '?'
  Experiment%ACQ_soft      = '?'
  Proposer%name            = '?'
  Proposer%team            = '?'
  Proposer%ISCR            = .false.
  ISCR%name                = '?'
  ISCR%team                = '?'
  ISCR%effectifs           = '?'
  ISCR%nb                  = 0


  TWINABS_file        = .false.
  TWINABS_4           = .false.
  SADABS%name         = 'SADABS'
  SADABS%version      = ''
  SADABS%line_ratio               = '?'
  SADABS%line_estimated_Tmin_Tmax = '?'
  SADABS%absorption_coeff         = '?'
  SADABS%ratio   = -999.
  SADABS%Tmin    = -999.
  SADABS%Tmax    = -999.
  SADABS%point_group  = '?'
  SADABS%nb_ref       = 0
  SADABS%nb_obs       = 0
  SADABS%nb_obs_3s    = 0
  SADABS%Rint         = 0.
  SADABS%Rint_3s      = 0.

  SHELX%name         = 'SHELX'
  SHELX%version      = ''
  SHELX%details      = '?'

  SIR%name           = '?'
  SIR%details        = '?'

  SPF%name           = '?'
  SPF%details        = '?'

  nb_da          = 0
  nb_ra          = 0
  nb_plane       = 0
  nb_ang_calc    = 0   ! nombre d'angles a calculer
  nb_bary_calc   = 0   ! nombre de barycentres a calculer
  nb_tolman_calc = 0   ! nombre d'angles de cone de Tolman a calculer

  atom1_dist    = ''
  atom2_dist    = ''
  atom1_diff    = ''
  atom2_diff    = ''
  atom1_ang     = ''
  atom2_ang     = ''
  atom3_ang     = ''
  atom4_ang     = ''
  atom_bary     = ''
  atom1_plane   = ''
  atom2_plane   = ''
  atom3_plane   = ''
  atom_plane    = ''
  atom_plane_nb = 0
  atom_tolman   = ''
  r_vdw         = -1.
  modif_rvdw    = .false.


  EXP_SCAN%nb               = 0
  EXP_SCAN%exposition_time  = 0.
  EXP_SCAN%frame_width      = 0.
  EXP_SCAN%DX               = 0.
  EXP_scan%theta            = 0.
  EXP_scan%omega            = 0.
  EXP_scan%phi              = 0.
  EXP_scan%chi              = 0.
  EXP_scan%kappa            = 0.
  EXP_SCAN%nb_frames        = 0
  EXP_SCAN%STL              = .false.
  EXP_SCAN%sfrm_exist       = .false.
  EXP_SCAN%axis             = "?"

  nb_atom           = 0
  nb_atom_no_H      = 0
  nb_atom_orbit     = 0
  nb_hkl            = 0
  nb_hkl_SFAC_calc  = 0
  nb_dist_calc      = 0
  nb_diff_calc      = 0

  CONN_dmin          = 0.5
  CONN_dmax          = 3.
  CONN_dmax_ini      = 3.
  CONN_dmin_ini      = 0.5
  CONN_all           = .false.
  CONN_all_X         = .false.
  CONN_ONLY_X        = .false.
  CONN_self          = .false.
  CONN_ang           = .false.
  CONN_out_condensed = .false.
  CONN_excluded      = .false.
  ATOM_CONN%label    = "?"
  ATOM_CONN%type     = "?"
  ATOM_CONN%excluded = "?"
  CONN_species       = "?"
  CONN_out           = .false.

  calcul_BVS   = .false.
  DIST_coef    = 0.
  DIST_plus    = 0.
  DIST_AH      = 0.9
  nb_sort      = 0
  nb_shell     = 0
  nb_symm_op   = 0
  wavelength   = 0.

  n_sig         = 3.
  threshold     = 0.03
  user_mat_nb   = 0

  sample_job             = "job"
  get_sample_ID          = .false.
  Create_INS%temperature = -999.
  Create_INS%U_threshold = -999.

  Max_ref                = Max_ref_allowed   ! Max. number of hkl reflections
  LOCK_wave_value        = 0.02
  CIF_format80           = .true.
  CIF_torsion_limit      = 170.
  include_RES_file       = .false.
  include_HKL_file       = .false.
  include_experimenter_name = .false.
  CMD_include_HKL_file   = 0
  update_parameters      = .true.
  report_header          = .true.
  report_logo(1)         = 'CDIFX_logo.jpg'
  report_logo(2)         = 'ISCR_logo.jpg'
  expert_mode            = .false.
  news_only_expert       = .false.
  skip_start_menu        = .true.
  hkl_statistics         = .true.
  hkl_format_free        = .false.
  hkl_format_SHELX       = .true.
  hkl_format             = '(3I4,2F8.2,I4,6F8.5)'
  hkl_SHELX_format       = hkl_format
  cartesian_frame%type   = "A"
  cartesian_frame%string = "x // a"

  pdp_simu%cu            = .true.
  pdp_simu%beam          = "x-rays"
  pdp_simu%wave          = 1.5406

  keep_bond_str_out      = .false.



  X_PV%U                  = 0.0055
  X_PV%V                  = -0.0015
  X_PV%W                  = 0.0036
  X_PV%eta0               = 0.3
  X_PV%eta1               = 0.
  X_pattern%WDT           = 20.
  X_pattern%background    = 50.
  X_pattern%scale         = 0.01
  X_pattern%step          = 0.01
  X_pattern%xmin          = 0
  X_pattern%xmax          = 120.
  X_pattern%job           = 2


  N_PV%U                  = 0.0146        ! D2B l=1.59 A a3=10 min
  N_PV%V                  = -0.0375
  N_PV%W                  = 0.0475
  N_PV%eta0               = 0.01
  N_PV%eta1               = 0.
  N_pattern%WDT           = 20.
  N_pattern%background    = 50.
  N_pattern%scale         = 100.
  N_pattern%step          = 0.025
  N_pattern%xmin          = 0
  N_pattern%xmax          = 140.
  N_pattern%job           = 3


  X_target(1:tabulated_target_nb)%logic  = .false.
  X_target(1:tabulated_target_nb)%write  = .true.
  X_target(1:tabulated_target_nb)%label   = (/'Ag',      'Mo',    'Cu',   'Ni', 'Co',    'Fe',    'Cr'/)
  ! Kalpha_1
  X_target(1:tabulated_target_nb)%wave(1) = (/ 0.55942,  0.70932, 1.5406, 1.65794, 1.78900, 1.93608, 2.28975 /)
  ! Kalpha_2
  X_target(1:tabulated_target_nb)%wave(2) = (/ 0.5638,   0.71361, 1.5444, 1.66178, 1.79289, 1.94002, 2.29365 /)
  ! Kbeta
  X_target(1:tabulated_target_nb)%wave(3) = (/ 0.4971,   0.63230, 1.3922, 1.50016, 1.62082, 1.75664, 2.08491 /)
  X_target(1:tabulated_target_nb)%tics_FILE = (/ 'xrays_tics_ag.pgf', 'xrays_tics_mo.pgf', 'xrays_tics_cu.pgf', &
                                                 'xrays_tics_ni.pgf', 'xrays_tics_co.pgf', 'xrays_tics_fe.pgf', &
                                                 'xrays_tics_cr.pgf' /)
  X_target(1:tabulated_target_nb)%cam_FILE  = (/ 'xrays_mac_ag.pgf',  'xrays_mac_mo.pgf',  'xrays_mac_cu.pgf',  &
                                                 'xrays_mac_cu.pgf',  'xrays_mac_co.pgf',  'xrays_mac_fe.pgf',  &
                                                 'xrays_mac_cr.pgf'  /)

  space_group_symbol  = 'ZZZ'

  input_INS            = .false.
  input_PCR            = .false.
  input_CIF            = .false.
  input_CFL            = .false.
  input_CELL%P4P       = .false.
  input_CELL%M50       = .false.
  input_CELL%INS       = .false.
  input_CELL%x         = .false.
  input_CELL%rmat      = .false.
  input_CELL%CIF       = .false.
  input_CELL%RED       = .false.
  input_CELL%file      = .false.

  mode_interactif         = .false.
  keyword_CELL            = .false.
  keyword_CELL_ESD        = .false.
  keyword_NIGGLI          = .false.
  keyword_WAVE            = .false.
  keyword_WRITE_CELL      = .false.
  write_cell_cart         = .false.
  keyword_WRITE_CHEM      = .false.
  keyword_WRITE_WAVE      = .false.
  keyword_WRITE_DEVICE    = .false.
  keyword_WRITE_SG        = .false.
  keyword_WRITE_BEAM      = .false.
  keyword_WRITE_QVEC      = .false.
  keyword_WRITE_SUPERCELL = .false.
  keyword_WRITE_ZUNIT     = .false.
  keyword_BEAM            = .false.
  keyword_SIZE            = .false.
  keyword_SPGR            = .false.
  get_SPGR                = .false.
  check_cent              = .true.
  check_all_SG            = .false.
  check_only_C            = .false.
  keyword_LSPGR           = .false.
   list_sg(1:7)          = .false.
   list_sg_Bravais(1:7)  = .false.
   list_sg_centric(1:2)  = .false.
   list_sg_laue(1:14)    = .false.
   list_sg_multip        = .false.
   list_sg_enantio       = .false.
   list_sg_chiral        = .false.
   list_sg_polar         = .false.
  keyword_LAUE           = .false.
  keyword_ATOM_list      = .false.
   write_atoms_cartesian = .false.
   write_atoms_in_A      = .false.
  keyword_ADP_list       = .false.
  ADP_details            = .false.
  ADP_comment            = ''
  keyword_SFAC_UNIT      = .false.
  keyword_CONT           = .false.
  keyword_CHEM           = .false.
  keyword_ZUNIT          = .false.
  keyword_MU             = .false.
  keyword_MAT            = .false.
  keyword_REDUCE         = .false.
  reduce_BL              = 'P'
  keyword_get_transf_mat = .false.
  keyword_HKLF5          = .false.
  keyword_LST_MAT        = .false.
  keyword_TRANSL         = .false.
  keyword_FILE           = .false.
  keyword_FILE_FCF       = .false.
  keyword_X_WAVE         = .false.
  keyword_P4P            = .false.
  keyword_RAW            = .false.
  keyword_HKL            = .false.
  keyword_HKL_diff       = .false.
  keyword_MATMUL         = .false.
  keyword_DIAG           = .false.
  keyword_VERSION        = .false.
  keyword_NO_details     = .false.
  keyword_Kappa_to_Euler = .false.
  keyword_Euler_to_Kappa = .false.
  motor                = 0.

  keyword_INSIDE       = .false.
  keyword_PAUSE        = .false.
  keyword_REM          = .false.

  HKL_file%name        = ''  ! nom du fichier.HKL
  HKL_file%output      = ''  ! fichier de sortie .HKL
  HKL_file%HKL         = ''  ! fichier .HKL
  HKL_file%d           = ''
  HKL_file%stl         = ''
  HKL_file%theta       = ''
  HKL_file%I           = ''
  HKL_file%Isig        = ''
  HKL_file%merge       = ''
  HKL_file%transf      = ''
  HKL_file%plot        = .false.
  HKL_file%read_NEG    = .false.
  HKL_file%CIF         = .false.
  HKL_file%SHELX       = .false.
  HKL_file%HKLF5       = .false.
  HKL_file%final_y     = .false.
  HKL_file%RAW         = .false.
  HKL_file%M91         = .false.
  HKL_file%M95         = .false.
  HKL_file%INT         = .false.
  HKL_file%COL         = .false.
  M95_comment_line     = ''
  M95_comment_line_nb  = 0

  keyword_STRUCT_factor = .false.

  keyword_read_CEL      = .false.
  keyword_read_CIF      = .false.
  keyword_read_INS      = .false.
  keyword_ABIN          = .false.
  keyword_read_HKLF5    = .false.
  clusters_out          = .false.
  keyword_read_FACES    = .false.
  read_Q_peaks          = .false.
  keyword_read_PCR      = .false.
  keyword_read_TIDY_out = .false.

  keyword_create_ACE       = .false.
  keyword_create_CEL       = .false.
  keyword_create_CFL       = .false.
  keyword_create_CIF       = .false.
  keyword_create_FST       = .false.
  keyword_create_PRF       = .false.
  keyword_create_PCR       = .false.
  create_FST_POLY          = .false.
  create_FST_MOLE          = .false.
  FST_no_H                 = .false.
  no_H                     = .false.
  launch_FP_Studio         = .false.
  keyword_create_INS       = .false.
  keyword_create_SOLVE     = .false.
  keyword_create_TIDY      = .false.
  create_CIF_PYMOL         = .false.
  create_SHAPE_file        = .false.
  poly_vol_calc            = .false.
  INI_create%ACE           = .false.
  INI_create%CEL           = .false.
  INI_create%CFL           = .false.
  INI_create%INS           = .false.
  INI_create%FST           = .false.
  INI_create%CIF_PYMOL     = .false.
  INI_create%PCR           = .false.
  INI_create%PRF           = .false.
  include_SQUEEZE          = .true.
  EXTRACT_hkl_ins_from_CIF = .false.

  keyword_create_REPORT  = .false.
  keyword_read_NREPORT   = .false.
  browse_nreport         = .false.
  keyword_modif_ARCHIVE  = .false.
  keyword_create_ARCHIVE = .false.
  keyword_SOLVE_to_INS   = .false.
  keyword_SIR            = .false.
  keyword_NEWS           = .false.
  news_year              = 'all'

  keyword_EDIT         = .false.
  file_to_edit         = ''
  keyword_SET          = .false.
  keyword_SETTING      = .false.
  keyword_save_setting = .false.

  keyword_QVEC         = .false.
  keyword_SUPERCELL    = .false.
  SUPERCELL_pcr        = .false.
  keyword_SYMM         = .false.
  keyword_BARY         = .false.      ! calcul du barycentre
  keyword_TOLMAN       = .false.
  keyword_DIFF         = .false.
  keyword_DIST         = .false.      ! calcul de distances interatomiques
  keyword_DIST_X       = .false.
  keyword_DIST_plus    = .false.
  keyword_DHA          = .false.
  keyword_CONN         = .false.
  keyword_ANG          = .false.      ! calcul d'angles
  keyword_PLANE        = .false.
  keyword_TRANSMISSION = .false.
  keyword_GENHKL       = .false.      ! generation de HKL
  keyword_GENSAT       = .false.      ! generation des satellites
  keyword_SFAC_HKL     = .false.      ! calcul de facteur de structure
  keyword_OBV_REV      = .FALSE.      ! analyse maclage obverse/reverse
  OBV_REV_twin_matrix  = 1            ! matrice de maclage:  (-1 0 0  0-1 0  0 0 1): twofold axis parallel to threefold '
 !OBV_REV_twin_matrix  = 2            ! matrice de maclage:  ( 0-1 0 -1 0 0  0 0-1): twofold axis parallel to (a-b)'
  write_HKL            = .false.      ! sortie des HKL
  create_PAT           = .false.      ! creation d'un diagramme
  PM2K_out             = .false.
  lecture_OK           = .false.

  keyword_SH_2th       = .false.
  SHIFT_2theta(1)      = 0.
  SHIFT_2theta(2)      = 0.
  SHIFT_2theta(3)      = 0.


  HKL_2THETA           = .false.
  HKL_THETA            = .false.
  HKL_STL              = .false.
  HKL_Q                = .false.
  HKL_D                = .false.

  X_min                = 0.
  X_max                = 0.

  WRITE(CIF_sep_line_empty, '(a,76a1,a)')  "#", (" ",i=1,CIF_line_width-2), "#"
  WRITE(CIF_sep_line_full,  '(a,76a1,a)')  "#", ("-",i=1,CIF_line_width-2), "#"

  CIF_diffrn_reflns%theta_min       = -1.
  CIF_diffrn_reflns%theta_max       = -1.
  CIF_cell_measurement%theta_min    = -1.
  CIF_cell_measurement%theta_max    = -1.
  CIF_cell_measurement%reflns_used  = -1.
  CIF_cell_measurement%temperature  = '?'
  CIF_cell_measurement%wavelength   = -1.

  CIF_dist%n_text          = 0
  CIF_dist%max_text        = 2000

  ! CIF parameters (extraits de archive.cif et utilises pour generer le fichier structural_REPORT.HTML
  CIF_parameter%sample_ID                        = '?'
  CIF_parameter%formula_moiety                   = '?'
  CIF_parameter%formula_sum                      = '?'
  CIF_parameter%formula_weight                   = '?'
  CIF_parameter%formula_units_Z                  = '?'

  CIF_parameter%symmetry_cell_setting            = '?'
  CIF_parameter%symmetry_space_group             = '?'
  CIF_parameter%symmetry_IT_number               = '?'
  CIF_parameter%sym_op                           = '?'
  CIF_parameter%cell_length_a                    = '?'
  CIF_parameter%cell_length_b                    = '?'
  CIF_parameter%cell_length_c                    = '?'
  CIF_parameter%cell_angle_alpha                 = '?'
  CIF_parameter%cell_angle_beta                  = '?'
  CIF_parameter%cell_angle_gamma                 = '?'
  CIF_parameter%cell_volume                      = '?'
  CIF_parameter%cell_formula_units_Z             = '?'
  CIF_parameter%exptl_density                    = '?'
  CIF_parameter%exptl_mu                         = '?'
  CIF_parameter%diffrn_reflns_number             = '?'
  CIF_parameter%diffrn_reflns_av_R_equivalents   = '?'
  CIF_parameter%diffrn_reflns_av_R_sigma         = '?'

  CIF_parameter%reflns_number_total              = '?'
  CIF_parameter%reflns_number_gt                 = '?'
  CIF_parameter%refine_ls_number_parameters      = '?'
  CIF_parameter%refine_ls_wR_factor_gt           = '?'
  CIF_parameter%refine_ls_R_factor_gt            = '?'
  CIF_parameter%refine_ls_R_factor_all           = '?'
  CIF_parameter%refine_ls_wR_factor_ref          = '?'
  CIF_parameter%refine_diff_density_max          = '?'
  CIF_parameter%refine_diff_density_min          = '?'
  CIF_parameter%refine_diff_density_rms          = '?'
  CIF_parameter%refine_ls_weighting_details      = '?'
  CIF_parameter%atom_sites_solution_1            = '?'
  CIF_parameter%atom_sites_solution_2            = '?'
  CIF_parameter%atom_sites_solution_H            = '?'
  CIF_parameter%refine_ls_H_treatment            = '?'
  CIF_parameter%refine_ls_extinction_method      = '?'
  CIF_parameter%refine_ls_extinction_coef        = '?'
  CIF_parameter%refine_ls_number_reflns          = '?'
  CIF_parameter%refine_ls_number_parameters      = '?'
  CIF_parameter%refine_ls_number_restraints      = '?'
  CIF_parameter%refine_ls_goodness_of_fit_ref    = '?'
  CIF_parameter%refine_ls_restrained_S_all       = '?'
  CIF_parameter%refine_ls_shift_su_max           = '?'
  CIF_parameter%refine_ls_shift_su_mean          = '?'

  CIF_parameter%H_treatment                      = 'constr'
  CIF_parameter%atom                             = '?'
  CIF_parameter%distance                         = '?'
  CIF_parameter%angle                            = '?'
  CIF_parameter%torsion_angle                    = '?'
  CIF_parameter%Hbond                            = '?'
  CIF_parameter%theta_min                        = '?'
  CIF_parameter%theta_max                        = '?'
  CIF_parameter%cell_theta_min                   = '?'
  CIF_parameter%cell_theta_max                   = '?'
  CIF_parameter%cell_reflns_used                 = '?'
  CIF_parameter%F000                             = '?'
  CIF_parameter%crystal_size_min                 = '?'
  CIF_parameter%crystal_size_mid                 = '?'
  CIF_parameter%crystal_size_max                 = '?'
  CIF_parameter%crystal_colour                   = '?'
  CIF_parameter%h_min                            = '?'
  CIF_parameter%h_max                            = '?'
  CIF_parameter%k_min                            = '?'
  CIF_parameter%k_max                            = '?'
  CIF_parameter%l_min                            = '?'
  CIF_parameter%l_max                            = '?'
  CIF_parameter%theta_full                       = '?'
  CIF_parameter%theta_max                        = '?'
  CIF_parameter%completeness                     = '?'
  CIF_parameter%absorption_correction_type       = '?'
  CIF_parameter%absorption_coefficient_mu        = '?'
  CIF_parameter%T_min                            = '?'
  CIF_parameter%T_max                            = '?'
  CIF_parameter%restraints_number                = '?'
  CIF_parameter%CHI2                             = '?'
  CIF_parameter%crystal_system                   = '?'
  CIF_parameter%Bravais                          = '?'
  CIF_parameter%abs_structure_Flack              = '?'
  CIF_parameter%abs_structure_details            = '?'
  CIF_parameter%shelx_res_checksum               = '?'
  CIF_parameter%shelx_hkl_checksum               = '?'
  CIF_parameter%shelx_fab_checksum               = '?'
  CIF_parameter%Friedel_coverage                 = "?"
  CIF_parameter%Friedel_fraction_full            = "?"
  CIF_parameter%Friedel_fraction_max             = "?"

  CIF_parameter%WinGX_used                           = .false.

  CIF_parameter_DEVICE%diffracto_device                  = '?'
  CIF_parameter_DEVICE%diffracto_radiation_type          = '?'
  CIF_parameter_DEVICE%diffracto_radiation_source        = '?'
  CIF_parameter_DEVICE%diffracto_radiation_wavelength    = '?'
  CIF_parameter_DEVICE%diffracto_radiation_probe         = '?'
  CIF_parameter_DEVICE%diffracto_radiation_monochromator = '?'
  CIF_parameter_DEVICE%diffracto_temperature             = '?'

  CIF_parameter_DEVICE%diffrn_source                     = '?'
  CIF_parameter_DEVICE%diffrn_radiation_wavelength       = '?'
  CIF_parameter_DEVICE%diffrn_radiation_type             = '?'
  CIF_parameter_DEVICE%diffrn_radiation_source           = '?'
  CIF_parameter_DEVICE%diffrn_radiation_monochromator    = '?'
  CIF_parameter_DEVICE%diffrn_radiation_probe            = '?'

  CIF_parameter_DEVICE%diffrn_measurement_device        = '?'
  CIF_parameter_DEVICE%diffrn_measurement_device_type   = '?'
  CIF_parameter_DEVICE%diffrn_measurement_method        = '?'
  CIF_parameter_DEVICE%diffrn_detector                  = '?'
  CIF_parameter_DEVICE%diffrn_detector_area_resol_mean  = '?'
  CIF_parameter_DEVICE%diffrn_theta_full                = '?'
  CIF_parameter_DEVICE%diffrn_theta_max                 = '?'

  CIF_parameter_DEVICE%computing_data_collection        = '?'
  CIF_parameter_DEVICE%computing_cell_refinement        = '?'
  CIF_parameter_DEVICE%computing_data_reduction         = '?'
  CIF_parameter_DEVICE%computing_structure_solution     = '?'
  CIF_parameter_DEVICE%computing_structure_refinement   = '?'
  CIF_parameter_DEVICE%computing_molecular_graphics     = '?'
  CIF_parameter_DEVICE%computing_publication_material_1 = '?'
  CIF_parameter_DEVICE%computing_publication_material_2 = '?'

  write(CIF_parameter_DEVICE%computing_structure_solution,   '(3a)') "'", trim(structure_solution%CIF_ref),   "'"
  write(CIF_parameter_DEVICE%computing_structure_refinement, '(3a)') "'", trim(structure_refinement%CIF_ref), "'"
  !CIF_parameter_DEVICE%computing_molecular_graphics         = "'Ortep-3 for Windows (Farrugia, 1997)'"
  CIF_parameter_DEVICE%computing_molecular_graphics         = "'SXGRAPH (Farrugia, 1999), Mercury (CSD, 2016)'"
  CIF_parameter_DEVICE%computing_publication_material_1     = "WinGX publication routines (Farrugia, 2012),"
  CIF_parameter_DEVICE%computing_publication_material_2     = "CRYSCALC (T. Roisnel, local program, 2016)"

  ! -------------- WinGX features -----------------------------------------
  ! ref   : Farrugia L. J., J. Appl. Crystallogr. 1999, 32, 837-838
  ! title : WinGX suite for small-molecule single-crystal crystallography
  ! DOI   : 10.1107/S0021889899006020
  ! -----------------------------------------------------------------------




  SW_DENZO%data_collection                              = "'Collect (Nonius BV, 1997-2000)'"
  SW_DENZO%cell_refinement                              = "'HKL Scalepack (Otwinowski & Minor 1997)'"
  SW_DENZO%data_reduction                               = "'HKL Denzo and Scalepack (Otwinowski & Minor 1997)'"
  SW_EVAL%data_collection                               = "'Collect (Nonius BV, 1997-2000)'"
  SW_EVAL%cell_refinement                               = "'Dirax/lsq (Duisenberg & Schreurs, 1989-2000)'"
  SW_EVAL%data_reduction                                = "'EvalCCD (Duisenberg & Schreurs 1990-2000)'"
  !SW_APEX%data_collection                               = "'Bruker SMART (Bruker, 2002)'"
  !SW_APEX%cell_refinement                               = "'Bruker SMART (Bruker, 2002)'"
  !SW_APEX%data_reduction                                = "'Bruker SAINT (Bruker, 2002)'"
  SW_APEX%data_collection                               = "'Bruker APEX2 (Bruker, 2014)'"
  SW_APEX%cell_refinement                               = "'Bruker APEX2 (Bruker, 2014)'"
  SW_APEX%data_reduction                                = "'Bruker APEX2 (Bruker, 2014)'"

  SW_X2S%data_collection                               = "'GIS (Bruker)'"
  SW_X2S%cell_refinement                               = "'APEX2 (Bruker, 2010); SAINT (Bruker, 2009)'"
  SW_X2S%data_reduction                                = "'SAINT (Bruker, 2009); XPREP (Sheldrick, 2008)'"

  SW_D8V_CU%data_collection                            = "'Bruker APEX3 (Bruker, 2015)'"
  SW_D8V_CU%cell_refinement                            = "'APEX3 (Bruker, 2015)'"
  SW_D8V_CU%data_reduction                             = "'SAINT (Bruker, 2014)'"
  SW_D8V_MO%data_collection                            = "'Bruker APEX3 (Bruker, 2015)'"
  SW_D8V_MO%cell_refinement                            = "'APEX3 (Bruker, 2015); SAINT (Bruker, 2014)'"
  SW_D8V_MO%data_reduction                             = "'SAINT (Bruker, 2014)'"

  SW_XCALIBUR%data_collection                          = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"
  SW_XCALIBUR%cell_refinement                          = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"
  SW_XCALIBUR%data_reduction                           = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"
  SW_SUPERNOVA%data_collection                         = "'CrysAlisPro, Agilent Technologies,Version 1.171.36'"
  SW_SUPERNOVA%cell_refinement                         = "'CrysAlisPro, Agilent Technologies,Version 1.171.36'"
  SW_SUPERNOVA%data_reduction                          = "'CrysAlisPro, Agilent Technologies,Version 1.171.36'"


  CIF_parameter_DEVICE_KCCD%diffrn_source                        = "'Enraf Nonius FR590'"
  CIF_parameter_DEVICE_KCCD%diffrn_radiation_wavelength          = '0.71073'
  CIF_parameter_DEVICE_KCCD%diffrn_radiation_type                = 'MoK\a'
  CIF_parameter_DEVICE_KCCD%diffrn_radiation_source              = "'fine-focus sealed tube'"
  CIF_parameter_DEVICE_KCCD%diffrn_radiation_monochromator       = 'graphite'
  CIF_parameter_DEVICE_KCCD%diffrn_radiation_probe               = 'x-ray'
  CIF_parameter_DEVICE_KCCD%diffrn_detector                      = "'CCD plate'"
  CIF_parameter_DEVICE_KCCD%diffrn_detector_area_resol_mean      = '9'
  CIF_parameter_DEVICE_KCCD%diffrn_measurement_device            = "'95mm CCD camera on \k-goniostat'"
  CIF_parameter_DEVICE_KCCD%diffrn_measurement_device_type       = "'KappaCCD, Nonius'"
  CIF_parameter_DEVICE_KCCD%diffrn_measurement_method            = "'CCD rotation images, thin slices'"
  CIF_parameter_DEVICE_KCCD%computing_data_collection            = SW_DENZO%data_collection
  CIF_parameter_DEVICE_KCCD%computing_cell_refinement            = SW_DENZO%cell_refinement
  CIF_parameter_DEVICE_KCCD%computing_data_reduction             = SW_DENZO%data_reduction
  CIF_parameter_DEVICE_KCCD%computing_structure_solution         = '?'
  CIF_parameter_DEVICE_KCCD%computing_structure_refinement       = '?'
  CIF_parameter_DEVICE_KCCD%computing_molecular_graphics         = '?'
  CIF_parameter_DEVICE_KCCD%computing_publication_material_1     = '?'
  CIF_parameter_DEVICE_KCCD%computing_publication_material_2     = '?'


  CIF_parameter_DEVICE_APEX%diffrn_measurement_device_type       = "'APEXII, Bruker-AXS'"
  CIF_parameter_DEVICE_APEX%diffrn_measurement_method            = "'CCD rotation images, thin slices'"
  CIF_parameter_DEVICE_APEX%diffrn_detector                      = "'CCD plate'"
  CIF_parameter_DEVICE_APEX%diffrn_radiation_wavelength          = '0.71073'
  CIF_parameter_DEVICE_APEX%diffrn_radiation_type                = 'MoK\a'
  CIF_parameter_DEVICE_APEX%diffrn_radiation_source              = "'fine-focus sealed tube'"
  CIF_parameter_DEVICE_APEX%diffrn_radiation_monochromator       = 'graphite'
  CIF_parameter_DEVICE_APEX%diffrn_radiation_probe               = 'x-ray'
  CIF_parameter_DEVICE_APEX%diffrn_measurement_device            = '?'
  CIF_parameter_DEVICE_APEX%diffrn_source                        = '?'
  CIF_parameter_DEVICE_APEX%diffrn_detector_area_resol_mean      = '?'
  CIF_parameter_DEVICE_APEX%computing_data_collection            = SW_APEX%data_collection
  CIF_parameter_DEVICE_APEX%computing_cell_refinement            = SW_APEX%cell_refinement
  CIF_parameter_DEVICE_APEX%computing_data_reduction             = SW_APEX%data_reduction

  CIF_parameter_DEVICE_APEX%computing_structure_solution         = CIF_parameter_DEVICE%computing_structure_solution
  CIF_parameter_DEVICE_APEX%computing_structure_refinement       = CIF_parameter_DEVICE%computing_structure_refinement
  CIF_parameter_DEVICE_APEX%computing_molecular_graphics         = CIF_parameter_DEVICE%computing_molecular_graphics
  CIF_parameter_DEVICE_APEX%computing_publication_material_1     = CIF_parameter_DEVICE%computing_publication_material_1
  CIF_parameter_DEVICE_APEX%computing_publication_material_2     = CIF_parameter_DEVICE%computing_publication_material_2


  CIF_parameter_DEVICE_D8V_CU%diffrn_measurement_device_type       = "'D8 VENTURE Bruker AXS'"
  CIF_parameter_DEVICE_D8V_CU%diffrn_measurement_method            = "'rotation images, thin slices'"
  CIF_parameter_DEVICE_D8V_CU%diffrn_detector                      = "'(CMOS) PHOTON 100'"
  CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_wavelength          = '1.54178'
  CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_type                = 'CuK\a'
  CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_probe               = 'x-ray'
  CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_source              = "'Incoatec microfocus sealed tube'"
  CIF_parameter_DEVICE_D8V_CU%diffrn_radiation_monochromator       = "'multilayer monochromator'"
  CIF_parameter_DEVICE_D8V_CU%diffrn_measurement_device            = '?'
  CIF_parameter_DEVICE_D8V_CU%diffrn_source                        = 'Incoatec microfocus sealed tube'
  CIF_parameter_DEVICE_D8V_CU%diffrn_detector_area_resol_mean      = '?'
  CIF_parameter_DEVICE_D8V_CU%computing_data_collection            = SW_D8V_CU%data_collection
  CIF_parameter_DEVICE_D8V_CU%computing_cell_refinement            = SW_D8V_CU%cell_refinement
  CIF_parameter_DEVICE_D8V_CU%computing_data_reduction             = SW_D8V_CU%data_reduction
  CIF_parameter_DEVICE_D8V_CU%computing_structure_solution         = CIF_parameter_DEVICE%computing_structure_solution
  CIF_parameter_DEVICE_D8V_CU%computing_structure_refinement       = CIF_parameter_DEVICE%computing_structure_refinement
  CIF_parameter_DEVICE_D8V_CU%computing_molecular_graphics         = CIF_parameter_DEVICE%computing_molecular_graphics
  CIF_parameter_DEVICE_D8V_CU%computing_publication_material_1     = CIF_parameter_DEVICE%computing_publication_material_1
  CIF_parameter_DEVICE_D8V_CU%computing_publication_material_2     = CIF_parameter_DEVICE%computing_publication_material_2


  CIF_parameter_DEVICE_D8V_MO%diffrn_measurement_device_type       = "'D8 VENTURE Bruker AXS'"
  CIF_parameter_DEVICE_D8V_MO%diffrn_measurement_method            = "'rotation images'"
  CIF_parameter_DEVICE_D8V_MO%diffrn_detector                      = "'(CMOS) PHOTON 100'"
  CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_wavelength          = '0.71073'
  CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_type                = 'MoK\a'
  CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_probe               = 'x-ray'
  CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_source              = "'Incoatec microfocus sealed tube'"
  CIF_parameter_DEVICE_D8V_MO%diffrn_radiation_monochromator       = "'multilayer monochromator'"
  CIF_parameter_DEVICE_D8V_MO%diffrn_measurement_device            = '?'
  CIF_parameter_DEVICE_D8V_MO%diffrn_source                        = 'Incoatec microfocus sealed tube'
  CIF_parameter_DEVICE_D8V_MO%diffrn_detector_area_resol_mean      = '?'
  CIF_parameter_DEVICE_D8V_MO%computing_data_collection            = SW_D8V_MO%data_collection
  CIF_parameter_DEVICE_D8V_MO%computing_cell_refinement            = SW_D8V_MO%cell_refinement
  CIF_parameter_DEVICE_D8V_MO%computing_data_reduction             = SW_D8V_MO%data_reduction
  CIF_parameter_DEVICE_D8V_MO%computing_structure_solution         = CIF_parameter_DEVICE%computing_structure_solution
  CIF_parameter_DEVICE_D8V_MO%computing_structure_refinement       = CIF_parameter_DEVICE%computing_structure_refinement
  CIF_parameter_DEVICE_D8V_MO%computing_molecular_graphics         = CIF_parameter_DEVICE%computing_molecular_graphics
  CIF_parameter_DEVICE_D8V_MO%computing_publication_material_1     = CIF_parameter_DEVICE%computing_publication_material_1
  CIF_parameter_DEVICE_D8V_MO%computing_publication_material_2     = CIF_parameter_DEVICE%computing_publication_material_2

  CIF_parameter_DEVICE_X2S%diffrn_measurement_device_type       = "'Bruker SMART X2S benchtop'"
  CIF_parameter_DEVICE_X2S%diffrn_measurement_method            = "'omega scans'"
  CIF_parameter_DEVICE_X2S%diffrn_detector                      = "'CCD plate'"
  CIF_parameter_DEVICE_X2S%diffrn_radiation_wavelength          = '0.71073'
  CIF_parameter_DEVICE_X2S%diffrn_radiation_type                = 'MoK\a'
  CIF_parameter_DEVICE_X2S%diffrn_radiation_source              = "'XOS X-beam microfocus source'"
  CIF_parameter_DEVICE_X2S%diffrn_radiation_monochromator       = "'doubly curved silicon crystal'"
  CIF_parameter_DEVICE_X2S%diffrn_radiation_probe               = 'x-ray'
  CIF_parameter_DEVICE_X2S%diffrn_measurement_device            = '?'
  CIF_parameter_DEVICE_X2S%diffrn_source                        = '?'
  CIF_parameter_DEVICE_X2S%diffrn_detector_area_resol_mean      = '?'
  CIF_parameter_DEVICE_X2S%computing_data_collection            = SW_X2S%data_collection
  CIF_parameter_DEVICE_X2S%computing_cell_refinement            = SW_X2S%cell_refinement
  CIF_parameter_DEVICE_X2S%computing_data_reduction             = SW_X2S%data_reduction

  CIF_parameter_DEVICE_X2S%computing_structure_solution         = CIF_parameter_DEVICE%computing_structure_solution
  CIF_parameter_DEVICE_X2S%computing_structure_refinement       = CIF_parameter_DEVICE%computing_structure_refinement
  CIF_parameter_DEVICE_X2S%computing_molecular_graphics         = CIF_parameter_DEVICE%computing_molecular_graphics
  CIF_parameter_DEVICE_X2S%computing_publication_material_1     = CIF_parameter_DEVICE%computing_publication_material_1
  CIF_parameter_DEVICE_X2S%computing_publication_material_2     = CIF_parameter_DEVICE%computing_publication_material_2




  CIF_parameter_DEVICE_XCALIBUR%diffrn_measurement_device_type   = "'CCD Saphire 3 Xcalibur'"
  CIF_parameter_DEVICE_XCALIBUR%diffrn_measurement_method        = "'omega scan'"
  CIF_parameter_DEVICE_XCALIBUR%diffrn_detector                  = "'CCD plate'"
  CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_wavelength      = '0.71073'
  CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_type            = 'MoK\a'
  CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_monochromator   = 'graphite'
  CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_probe           = 'x-ray'
  CIF_parameter_DEVICE_XCALIBUR%diffrn_radiation_source          = "'sealed X-ray tube'"
  CIF_parameter_DEVICE_XCALIBUR%diffrn_detector_area_resol_mean  = '19.64'
  CIF_parameter_DEVICE_XCALIBUR%computing_data_collection        = SW_XCALIBUR%data_collection
  CIF_parameter_DEVICE_XCALIBUR%computing_cell_refinement        = SW_XCALIBUR%cell_refinement
  CIF_parameter_DEVICE_XCALIBUR%computing_data_reduction         = SW_XCALIBUR%data_reduction
  CIF_parameter_DEVICE_XCALIBUR%computing_structure_solution     = CIF_parameter_DEVICE%computing_structure_solution
  CIF_parameter_DEVICE_XCALIBUR%computing_structure_refinement   = CIF_parameter_DEVICE%computing_structure_refinement
  CIF_parameter_DEVICE_XCALIBUR%computing_molecular_graphics     = CIF_parameter_DEVICE%computing_molecular_graphics
  CIF_parameter_DEVICE_XCALIBUR%computing_publication_material_1 = CIF_parameter_DEVICE%computing_publication_material_1
  CIF_parameter_DEVICE_XCALIBUR%computing_publication_material_2 = CIF_parameter_DEVICE%computing_publication_material_2


  CIF_parameter_DEVICE_SUPERNOVA%diffrn_measurement_device_type   = "'SuperNova, Dual, Cu at zero, Atlas'"
  CIF_parameter_DEVICE_SUPERNOVA%diffrn_measurement_method        = "'omega scans'"
  CIF_parameter_DEVICE_SUPERNOVA%diffrn_detector                  = "'CCD plate'"
  CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_wavelength      = ''
  CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_monochromator   = 'multiplayer optics'
  CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_probe           = 'x-ray'
  CIF_parameter_DEVICE_SUPERNOVA%diffrn_radiation_source          = 'X-ray microsource'

  CIF_parameter_DEVICE_SUPERNOVA%diffrn_detector_area_resol_mean  = '5.2474'
  CIF_parameter_DEVICE_SUPERNOVA%computing_data_collection        = SW_SUPERNOVA%data_collection
  CIF_parameter_DEVICE_SUPERNOVA%computing_cell_refinement        = SW_SUPERNOVA%cell_refinement
  CIF_parameter_DEVICE_SUPERNOVA%computing_data_reduction         = SW_SUPERNOVA%data_reduction
  CIF_parameter_DEVICE_SUPERNOVA%computing_structure_solution     = CIF_parameter_DEVICE%computing_structure_solution
  CIF_parameter_DEVICE_SUPERNOVA%computing_structure_refinement   = CIF_parameter_DEVICE%computing_structure_refinement
  CIF_parameter_DEVICE_SUPERNOVA%computing_molecular_graphics     = CIF_parameter_DEVICE%computing_molecular_graphics
  CIF_parameter_DEVICE_SUPERNOVA%computing_publication_material_1 = CIF_parameter_DEVICE%computing_publication_material_1
  CIF_parameter_DEVICE_SUPERNOVA%computing_publication_material_2 = CIF_parameter_DEVICE%computing_publication_material_2

  SADABS%type       = "_exptl_absorpt_correction_type                       multi-scan"
  SADABS%details(1) = "_exptl_absorpt_process_details"
  SADABS%details(2) = ";"
  SADABS%details(3) = "  [Sheldrick, G.M. (2014). SADABS Bruker AXS Inc., Madison, Wisconsin, USA]"
  SADABS%details(4) = ";"

  SADABS%tw_details(1) = "_exptl_absorpt_process_details"
  SADABS%tw_details(2) = ";"
  SADABS%tw_details(3) = "  [Sheldrick, G.M. (2012). TWINABS Bruker AXS Inc., Madison, Wisconsin, USA]"
  SADABS%tw_details(4) = ";"


  SHELX%type        = ""
  SHELX%details(1)  = ";"
  SHELX%details(2)  = "_computing_structure_solution    'SHELXS-97:     Sheldrick G.M., Acta Cryst. A64 (2008), 112-122'"
  SHELX%details(3)  = "_computing_structure_solution    'SHELXT:        Sheldrick G.M., Acta Cryst. A71 (2015) 3-8'"
  SHELX%details(4)  = "_computing_structure_refinement  'SHELXL-97:     Sheldrick G.M., Acta Cryst. A64 (2008), 112-122'"
  SHELX%details(5)  = "_computing_structure_refinement  'SHELXL-2014/6: Sheldrick G.M., Acta Cryst. C71 (2015) 3-8'"
  SHELX%details(6)  = ";"

  SHELX%details(10) = "SHELXL-97"
  SHELX%details(11) = "  Sheldrick G.M., Acta Cryst. A64 (2008), 112-122"
  SHELX%details(12) = "  A short history of SHELX"
  SHELX%details(13) = "  DOI: http://dx.doi.org/10.1107/S0108767307043930"
  SHELX%details(14) = "SHELXL-2014"
  SHELX%details(15) = "  Sheldrick G.M., Acta Cryst. C71 (2015) 3-8"
  SHELX%details(16) = "  Crystal structure refinement with SHELXL"
  SHELX%details(17) = "  DOI: http://dx.doi.org/10.1107/S2053229614024218"
  SHELX%details(18) = "SHELXT"
  SHELX%details(19) = "  Sheldrick G.M., Acta Cryst. A71 (2015) 3-8"
  SHELX%details(20) = "  SHELXT - Integrated space-group and crystal-structure determination"
  SHELX%details(21) = "  DOI: http://dx.doi.org/10.1107/S2053273314026370"

  SIR%details(1)    = "SIR97"
  SIR%details(2)    = "  A. Altomare, M. C. Burla, M. Camalli, G. Cascarano, C. Giacovazzo, A. Guagliardi, A. G. G. Moliterni, "// &
                      "G. Polidori, R. Spagna, J. Appl. Cryst. (1999) 32, 115-119"
  SIR%details(3)    = "  SIR97: a new tool for crystal structure determination and refinement"
  SIR%details(4)    = "  DOI: http://dx.doi.org/10.1107/S0021889898007717"
  SIR%details(5)    = "SIR2004"
  SIR%details(6)    = "  M.C. Burla, R. Caliandro, M. Camalli, B. Carrozzini, G. Cascarano, L. De Caro, C. Giacovazzo, " // &
                      "G. Polidori, R. Spagna, J. Appl. Cryst. (2005) 38, 381-388"
  SIR%details(7)    = "  SIR2004: an improved tool for crystal structure determination and refinement"
  SIR%details(8)    = "  DOI: http://dx.doi.org/10.1107/S002188980403225X"
  SIR%details(9)    = "SIR2014"
  SIR%details(10)   = "  M.C. Burla, R. Caliandro, B. Carrozzini, G. L. Cascarano, C. Cuocci, C. Giacovazzo, M. Mallamo, "// &
                      "A. Mazzoneb and G. Polidori, J. Appl. Cryst. (2015). 48, 306–309"
  SIR%details(11)   = "  SIR2004: an improved tool for crystal structure determination and refinement"
  SIR%details(12)   = "  DOI: http://dx.doi.org/10.1107/S1600576715001132"

  SPF%details(1)    = "SUPERFLIP"
  SPF%details(2)    = "  Palatinus L., Chapuis G., J. Appl. Cryst. (2007, 40, 786-790"
  SPF%details(3)    = "  Superflip - a computer program for the solution of crystal structures by charge flipping " // &
                      "in arbitrary dimensions"
  SPF%details(4)    = "  DOI: http://dx.doi.org/10.1107/S0021889807029238"


  ABS_CRYSALIS%type       = "_exptl_absorpt_correction_type                          multi-scan"
  ABS_CRYSALIS%details(1) = "_exptl_absorpt_process_details"
  ABS_CRYSALIS%details(2) = ";"
  ABS_CRYSALIS%details(3) = "  CrysAlisPro, Agilent Technologies, Version 1.171.36.28c"
  ABS_CRYSALIS%details(4) = "  Empirical absorption correction using spherical harmonics"
  ABS_CRYSALIS%details(5) = "  emplemented in SCALE3 ABSPACK scaling algorithm."
  ABS_CRYSALIS%details(6) = ";"

  ! -----------------------------------------------------------------

  SQUEEZE%procedure        = .false.
  SQUEEZE%read_line(:)     = '?'
  SQUEEZE%unit             = 61
  SQUEEZE%file             = 'platon.sqf'
  SQUEEZE%nb_lines         = 0

  unit_cell%param          = 0.
  unit_cell%param_ESD      = 0.
  unit_cell%rec_param      = 0.
  unit_cell%new_param      = 0.
  unit_cell%new_param_ESD  = 0.
  unit_cell%Niggli         = 0.
  unit_cell%volume         = 0.
  unit_cell%volume_ESD     = 0.
  unit_cell%rec_volume     = 0.
  unit_cell%crystal_system = '?'
  unit_cell%Bravais        = '?'
  unit_cell%lattice        = '?'
  unit_cell%H_M            = '?'

  unit_cell_2%param        = 0.
  unit_cell_2%param_ESD    = 0.


  UB_matrix      = 0.
  UB_mat_log     = .false.
  keyword_UB_mat = .false.


  crystal%size           = 0.
  crystal%size_min       = 0.
  crystal%size_max       = 0.
  crystal%size_mid       = 0.
  crystal%volume         = 0.
  crystal%radius         = 0.
  crystal%mosaic         = 0.
  crystal%morph          = '?'
  crystal%color          = '?'
  crystal%density_meas   = '?'
  crystal%density_diffrn = '?'
  crystal%density_method = '?'
  crystal%F000           = '?'

  crystal_system       = '?'
  crystal%face_line(:) = '?'
  crystal%faces_nb     = 0
  crystal%face_index(:,:)    = 0
  crystal%face_dim(:)  = 0.

  ABSORPTION%mu        = 0.
  ABSORPTION%Tmin      = 0.
  ABSORPTION%Tmax      = 0.


  molecule%common_name = '?'
  molecule%formula     = '?'
  molecule%content     = '?'
  molecule%weight      = 0.0
  molecule%density     = 0.0
  molecule%Z           = 0
  molecule%Z_unit      = 1
  Z_unit               = 1.
  Z_unit_INS           = 1




  F000                 = 0.

  beam_type               = 'x-rays'
  X_rays                  = .true.
  neutrons                = .false.
  electrons               = .false.

  Mat_integer             = .false.
  sort_plot               = .false.
  sort_out                = .false.
  file_out                = .false.
  sort_out_n              = 1
  file_out_n              = 1
  shell_plot              = .false.
  shell_arg_min_max       = .false.
  shell_out               = .false.
  symm_op_xyz             = .false.
  symm_op_mat             = .false.

  WRITE_SPG_info             = .false.
  write_SPG_info_all         = .false.
  WRITE_SPG_exti             = .false.
  write_SPG_all              = .false.
  write_SPG_subgroups        = .false.
  WRITE_symm_op              = .false.
  WRITE_symm_op_condensed    = .false.
  WRITE_SYMM_op_on_one_line  = .false.
  WRITE_SHELX_symm_op        = .false.
  WRITE_APPLY_symm           = .false.
  WRITE_STAR_K               = .false.
  routine_TR                 = .false.
  WRITE_SITE_info            = .false.
  WRITE_PCR_site_info        = .false.
  WRITE_PCR_mag_site_info    = .false.
  WRITE_triclinic_transf     = .false.
  WRITE_monoclinic_transf    = .false.
  !WRITE_monoclinic_out       = .true.
  WRITE_details              = .true.
  WRITE_rhomb_hex_transf     = .false.
  WRITE_hex_rhomb_transf     = .false.
  WRITE_permutation_abc      = .false.
  WRITE_twin_hexa            = .false.
  WRITE_twin_pseudo_hexa     = .false.
  ONLY_monoclinic_angle      = .false.
  ONLY_monoclinic_angle_ini  = .false.
  keyword_LST_MAT            = .false.

  keyword_THERM           = .false.
  keyword_THERM_SHELX     = .false.
  THERM%Uiso              = .false.
  THERM%Biso              = .false.
  THERM%Uij               = .false.
  THERM%Bij               = .false.
  THERM%Beta              = .false.
  THERM%aniso             = .false.


  keyword_KEY             = .false.
  keyword_HELP            = .false.
  keyword_CLA             = .false.
  keyword_HEADER          = .false.
  keyword_SYST            = .false.
  keyword_RESET           = .false.

  keyword_create_CRYSCALC_news   = .false.
  keyword_create_CRYSCALC_HTML   = .false.
  browse_cryscalc_HTML           = .false.

  keyword_WRITE_REF_APEX          = .false.
  keyword_WRITE_REF_DENZO         = .false.
  keyword_WRITE_REF_EVAL          = .false.
  keyword_write_REF_KCCD          = .false.
  keyword_WRITE_REF_SHELX         = .false.
  keyword_WRITE_REF_SADABS        = .false.
  keyword_WRITE_REF_ABS_CRYSALIS  = .false.
  keyword_write_REF_SUPERNOVA     = .false.
  keyword_WRITE_REF_X2S           = .false.
  keyword_WRITE_REF_D8V_Cu        = .false.
  keyword_WRITE_REF_D8V_Mo        = .false.
  keyword_write_REF_XCALIBUR      = .false.
  WRITE_ref_CIF                   = .false.


  unknown_keyword              = .false.
  unknown_CFL_keyword          = .false.
  keyword_SEARCH_EXTI          = .false.
  keyword_SEARCH_SPGR          = .false.
  keyword_SEARCH_MONO          = .false.
  keyword_SEARCH_TETRA         = .false.
  keyword_FIND_HKL             = .false.
  keyword_FIND_HKL_EQUIV       = .false.
  keyword_FIND_HKL_ABSENT      = .false.
  keyword_FIND_HKL_NEG         = .false.
  keyword_FIND_HKL_POS         = .false.
  keyword_FIND_HKL_LIST        = .false.

  keyword_LIST_EXTI_RULE  = .false.
  search_H_string         = .false.
  search_equiv            = .false.
  search_friedel          = .false.
  search_mono_criteria(1) = 2.5    ! ecart max. / 90.
  search_mono_criteria(2) = 0.2    ! ecart du Rint_mono_min /  Rint_mono(1:3)
  search_mono_criteria(2) = 0.2    ! ecart du Rint_mono_min /  Rint_mono(1:3)
  ONLY_monoclinic_angle     = .false.
  ONLY_monoclinic_angle_INI = .false.

  search_tetra_criteria(1) = 0.5   ! ecart max. / <cell_param>
  search_tetra_criteria(2) = 2.5   ! ecart max. / 90.
  search_tetra_criteria(3) = 0.2   ! ecart du Rint_tetra_min /  Rint_tetra(1:3)

  overlapping_criteria     = 0.15


  ordered_HKL                      = .false.
  keyword_Rint                     = .false.
  keyword_merge_HKL                = .false.
  keyword_get_Friedel_pairs_number = .false.

  keyword_STL         = .false.
  keyword_dhkl        = .false.
  keyword_theta       = .false.
  keyword_2theta      = .false.
  keyword_Qhkl        = .false.
  nb_STL_value    = 0
  nb_dhkl_value   = 0
  nb_Qhkl_value   = 0
  nb_theta_value  = 0
  nb_2theta_value = 0

  keyword_DIR_ANG = .false.
  keyword_REC_ANG = .false.

  keyword_mendel = .false.
  mendel_plot    = .false.
  mendel_atom_nb = 0
  mendel_atom_label = ''

  keyword_shannon    = .false.
  shannon_atom_label = ''

  keyword_mag        = .false.
  mag_atom_label     = ''

  keyword_DATA_neutrons     = .false.
  data_neutrons_PLOT        = .false.
  DATA_neutrons_RE_PLOT_ALL = .false.
  DATA_n_RE                 = .false.

  keyword_DATA_Xrays    = .false.
  data_Xrays_PLOT       = .false.

  keyword_DATA_atomic_density  = .false.
  data_atomic_density_PLOT     = .false.

  keyword_DATA_atomic_radius   = .false.
  data_atomic_radius_PLOT      = .false.

  keyword_DATA_atomic_weight   = .false.
  data_atomic_weight_PLOT      = .false.

  known_atomic_label    = .false.
  known_atomic_features = .false.
  known_data_neutrons   = .false.
  known_data_X          = .false.
  known_space_groups    = .false.
  known_theta           = .false.
  known_cell_esd        = .false.
  known_shannon_lines   = .false.
  known_mag_lines       = .false.

  keyword_WEB           = .false.
  keyword_download_EXE  = .false.
  URL_address           = ''

  translat          = 0.
  symm_op_rot       = 0.
  symm_op_trans     = 0.

  nb_help = nb_help_max
  write_keys = .true.
                                   !1                      2                      3                      4                    5
  HELP_string( 1: nb_help_max) = &
  (/'ABSENT_HKL         ', 'ABSORPTION         ', 'ACTA               ', 'ANG                ', 'APPLY_OP           ', &
    'ATOM               ', 'ATOM_LIST          ', 'BARY               ', 'BEAM               ', 'CELL               ', &
    'CELL_ESD           ', 'CHEM               ', 'CONN               ', 'CONT               ', 'CREATE_ACE         ', &
    'CREATE_CEL         ', 'CREATE_CFL         ', 'CREATE_FST         ', 'CREATE_INS         ', 'CREATE_PCR         ', &
    'CREATE_REPORT      ', 'CREATE_SOLVE       ', 'CREATE_TIDY        ', 'D_HKL              ', 'D_STAR             ', &
    'DATA_ATOMIC_DENSITY', 'DATA_ATOMIC_RADIUS ', 'DATA_ATOMIC_WEIGHT ', 'DATA_NEUTRONS      ', 'DATA_XRAYS         ', &
    'DIAG_MAT           ', 'DIFF               ', 'DIR_ANG            ', 'DIST               ', 'DIST_DHA           ', &
    'EDIT               ', 'EQUIV_HKL          ', 'EULER_TO_KAPPA     ', 'EXIT               ', 'FCF_FILE           ', &
    'FILE               ', 'FIND_HKL           ', 'FIND_HKL_LIST      ', 'FRIEDEL            ', 'GEN_HKL            ', &
    'GET_TRANS_MAT      ', 'HEADER             ', 'HEX_RHOMB          ', 'HKL                ', 'HKL_DIFF           ', &
    'HKL_NEG            ', 'HKL_POS            ', 'HKLF5              ', 'INSIDE             ', 'KAPPA_TO_EULER     ', &
    'LIST_EXTI/LST_EXTI ', 'LIST_KEYS/LST_KEYS ', 'LIST_LAUE/LST_LAUE ', 'LIST_MATR/LST_MATR ', 'LIST_SG/LST_SG     ', &
    'MAG                ', 'MAN                ', 'MAN_HTML           ', 'MATMUL             ', 'MATR               ', &
    'MENDEL             ', 'MERGE              ', 'MONOCLINIC         ', 'NEWS               ', 'NIGGLI             ', &
    'OBV_REV            ', 'P4P                ', 'PAUSE              ', 'PERMUT             ', 'PLANE              ', &
    'Q_HKL              ', 'QVEC               ', 'READ_CEL           ', 'READ_CIF           ', 'READ_FACES         ', &
    'READ_HKLF5         ', 'READ_INS           ', 'READ_NREPORT       ', 'READ_PCR           ', 'READ_TIDY_OUT      ', &
    'REC_ANG            ', 'REDUCE_CELL        ', 'REF_ABS_CRYSALIS   ', 'REF_D8V_CU         ', 'REF_D8V_MO         ', &
    'REF_APEX           ', 'REF_DENZO          ', 'REF_EVAL           ', 'REF_KCCD           ', 'REF_SADABS         ', &
    'REF_SHELX          ', 'REF_SIR            ', 'REF_SUPERFLIP      ', 'REF_SUPERNOVA      ', 'REF_X2S            ', &
	'REF_XCALIBUR       ', 'RESET              ', 'RINT               ', 'RHOMB_HEX          ', 'SAVE_SETTINGS      ', &
	'SEARCH_EXTI        ', 'SEARCH_MONO        ', 'SEARCH_SPGR        ', 'SEARCH_TETRA       ', 'SET                ', &
	'SETTING            ', 'SFAC               ', 'SF_HKL             ', 'SG                 ', 'SG_ALL             ', &
	'SG_EXTI            ', 'SG_INFO            ', 'SG_SUB             ', 'SHANNON            ', 'SHELL              ', &
	'SHIFT_2TH          ', 'SITE_INFO          ', 'SIZE               ', 'SORT               ', 'STAR_K             ', &
	'STL                ', 'SUPERCELL          ', 'SYMM               ', 'SYST               ', 'THERM              ', &
	'THERM_SHELX        ', 'THETA              ', 'TITL               ', 'TOLMAN_ANGLE       ', 'TRANSLATION        ', &
	'TRANSMISSION       ', 'UB_MATRIX          ', 'UPDATE             ', 'USER_MAT           ', 'TRICLINIC          ', &
	'TWIN_HEXA          ', 'TWIN_PSEUDO_HEXA   ', 'TWO_THETA          ', 'UNIT               ', 'WAVE               ', &
	'WEB                ', 'WRITE_ADP          ', 'WRITE_BEAM         ', 'WRITE_CELL         ', 'WRITE_CHEM         ', &
	'WRITE_DEVICE       ', 'WRITE_QVEC         ', 'WRITE_SG           ', 'WRITE_SUPERCELL    ', 'WRITE_SYM_OP       ', &
	'WRITE_WAVE         ', 'WRITE_ZUNIT        ', 'X_WAVE             ', 'ZUNIT              ' /)

  HELP_arg(1:nb_help_max) = HELP_string(1:nb_help_max)

  numor = 1;            HELP_ABSENT_HKL_numor          =  numor  !  1
  numor = numor + 1;    HELP_ABSORPTION_numor          =  numor  !  2
  numor = numor + 1;    HELP_ACTA_numor                =  numor  !  3
  numor = numor + 1;    HELP_ANG_numor                 =  numor  !  4
  numor = numor + 1;    HELP_APPLY_OP_numor            =  numor  !  5
  numor = numor + 1;    HELP_ATOM_numor                =  numor  !  6
  numor = numor + 1;    HELP_ATOM_LIST_numor           =  numor  !  7
  numor = numor + 1;    HELP_BARY_numor                =  numor  !  8
  numor = numor + 1;    HELP_BEAM_numor                =  numor  !  9
  numor = numor + 1;    HELP_CELL_numor                =  numor  ! 10
  numor = numor + 1;    HELP_CELL_ESD_numor            =  numor  ! 11
  numor = numor + 1;    HELP_CHEM_numor                =  numor  ! 12
  numor = numor + 1;    HELP_CONN_numor                =  numor  ! 13
  numor = numor + 1;    HELP_CONT_numor                =  numor  ! 14
  numor = numor + 1;    HELP_CREATE_ACE_numor          =  numor  ! 15
  numor = numor + 1;    HELP_CREATE_CEL_numor          =  numor  ! 16
  numor = numor + 1;    HELP_CREATE_CFL_numor          =  numor  ! 17
  numor = numor + 1;    HELP_CREATE_FST_numor          =  numor  ! 18
  numor = numor + 1;    HELP_CREATE_INS_numor          =  numor  ! 19
  numor = numor + 1;    HELP_CREATE_PCR_numor          =  numor  ! 20
  numor = numor + 1;    HELP_CREATE_REPORT_numor       =  numor  ! 21
  numor = numor + 1;    HELP_CREATE_SOLVE_numor        =  numor  ! 22
  numor = numor + 1;    HELP_CREATE_TIDY_numor         =  numor  ! 23
  numor = numor + 1;    HELP_D_HKL_numor               =  numor  ! 24
  numor = numor + 1;    HELP_D_STAR_numor              =  numor  ! 25
  numor = numor + 1;    HELP_DATA_ATOMIC_DENSITY_numor =  numor  ! 26
  numor = numor + 1;    HELP_DATA_ATOMIC_RADIUS_numor  =  numor  ! 27
  numor = numor + 1;    HELP_DATA_ATOMIC_WEIGHT_numor  =  numor  ! 28
  numor = numor + 1;    HELP_DATA_NEUTRONS_numor       =  numor  ! 29
  numor = numor + 1;    HELP_DATA_XRAYS_numor          =  numor  ! 30
  numor = numor + 1;    HELP_DIAG_MAT_numor            =  numor  ! 31
  numor = numor + 1;    HELP_DIFF_numor                =  numor  ! 32
  numor = numor + 1;    HELP_DIR_ANG_numor             =  numor  ! 33
  numor = numor + 1;    HELP_DIST_numor                =  numor  ! 34
  numor = numor + 1;    HELP_DIST_DHA_numor            =  numor  ! 35
  numor = numor + 1;    HELP_EDIT_numor                =  numor  ! 36
  numor = numor + 1;    HELP_EQUIV_numor               =  numor  ! 37
  numor = numor + 1;    HELP_EULER_TO_KAPPA_numor      =  numor  ! 38
  numor = numor + 1;    HELP_EXIT_numor                =  numor  ! 39
  numor = numor + 1;    HELP_FCF_FILE_numor            =  numor  ! 40
  numor = numor + 1;    HELP_FILE_numor                =  numor  ! 41
  numor = numor + 1;    HELP_FIND_HKL_numor            =  numor  ! 42
  numor = numor + 1;    HELP_FIND_HKL_LIST_numor       =  numor  ! 43
  numor = numor + 1;    HELP_FRIEDEL_pairs_numor       =  numor  ! 44
  numor = numor + 1;    HELP_GEN_HKL_numor             =  numor  ! 45
  numor = numor + 1;    HELP_GET_TRANS_MAT_numor       =  numor  ! 46
  numor = numor + 1;    HELP_HEADER_numor              =  numor  ! 47
  numor = numor + 1;    HELP_HEX_RHOMB_numor           =  numor  ! 48
  numor = numor + 1;    HELP_HKL_numor                 =  numor  ! 49
  numor = numor + 1;    HELP_HKL_diff_numor            =  numor  ! 50
  numor = numor + 1;    HELP_HKL_NEG_numor             =  numor  ! 51
  numor = numor + 1;    HELP_HKL_POS_numor             =  numor  ! 52
  numor = numor + 1;    HELP_HKLF5_numor               =  numor  ! 53
  numor = numor + 1;    HELP_INSIDE_numor              =  numor  ! 54
  numor = numor + 1;    HELP_KAPPA_TO_EULER_numor      =  numor  ! 55
  numor = numor + 1;    HELP_LIST_EXTI_numor           =  numor  ! 56
  numor = numor + 1;    HELP_LIST_KEYS_numor           =  numor  ! 57
  numor = numor + 1;    HELP_LIST_LAUE_numor           =  numor  ! 58
  numor = numor + 1;    HELP_LIST_MATR_numor           =  numor  ! 59
  numor = numor + 1;    HELP_LIST_SG_numor             =  numor  ! 60
  numor = numor + 1;    HELP_MAG_numor                 =  numor  ! 61
  numor = numor + 1;    HELP_MAN_numor                 =  numor  ! 62
  numor = numor + 1;    HELP_MAN_HTML_numor            =  numor  ! 63
  numor = numor + 1;    HELP_MATMUL_numor              =  numor  ! 64
  numor = numor + 1;    HELP_MATR_numor                =  numor  ! 65
  numor = numor + 1;    HELP_MENDEL_numor              =  numor  ! 66
  numor = numor + 1;    HELP_MERGE_numor               =  numor  ! 67
  numor = numor + 1;    HELP_MONOCLINIC_numor          =  numor  ! 68
  numor = numor + 1;    HELP_NEWS_numor                =  numor  ! 69
  numor = numor + 1;    HELP_NIGGLI_CELL_numor         =  numor  ! 70
  numor = numor + 1;    HELP_OBV_REV_numor             =  numor  ! 71
  numor = numor + 1;    HELP_P4P_numor                 =  numor  ! 72
  numor = numor + 1;    HELP_PAUSE_numor               =  numor  ! 73
  numor = numor + 1;    HELP_PERMUT_numor              =  numor  ! 74
  numor = numor + 1;    HELP_PLANE_numor               =  numor  ! 75
  numor = numor + 1;    HELP_Q_HKL_numor               =  numor  ! 76
  numor = numor + 1;    HELP_QVEC_numor                =  numor  ! 77
  numor = numor + 1;    HELP_READ_CEL_numor            =  numor  ! 78
  numor = numor + 1;    HELP_READ_CIF_numor            =  numor  ! 79
  numor = numor + 1;    HELP_READ_FACES_numor          =  numor  ! 80
  numor = numor + 1;    HELP_READ_HKLF5_file_numor     =  numor  ! 81
  numor = numor + 1;    HELP_READ_INS_numor            =  numor  ! 82
  numor = numor + 1;    HELP_READ_NREPORT_numor        =  numor  ! 83
  numor = numor + 1;    HELP_READ_PCR_numor            =  numor  ! 84
  numor = numor + 1;    HELP_READ_TIDY_out_numor       =  numor  ! 85
  numor = numor + 1;    HELP_REC_ANG_numor             =  numor  ! 86
  numor = numor + 1;    HELP_REDUCE_CELL_numor         =  numor  ! 87
  numor = numor + 1;    HELP_REF_ABS_CRYSALIS_numor    =  numor  ! 88
  numor = numor + 1;    HELP_REF_APEX_numor            =  numor  ! 89
  numor = numor + 1;    HELP_REF_D8V_Cu_numor          =  numor  ! 90
  numor = numor + 1;    HELP_REF_D8V_Mo_numor          =  numor  ! 91
  numor = numor + 1;    HELP_REF_DENZO_numor           =  numor  ! 92
  numor = numor + 1;    HELP_REF_EVAL_numor            =  numor  ! 93
  numor = numor + 1;    HELP_REF_KCCD_numor            =  numor  ! 94
  numor = numor + 1;    HELP_REF_SADABS_numor          =  numor  ! 95
  numor = numor + 1;    HELP_REF_SHELX_numor           =  numor  ! 96
  numor = numor + 1;    HELP_REF_SIR_numor             =  numor  ! 97
  numor = numor + 1;    HELP_REF_SPF_numor             =  numor  ! 98
  numor = numor + 1;    HELP_REF_SUPERNOVA_numor       =  numor  ! 99
  numor = numor + 1;    HELP_REF_X2S_numor             =  numor  !100
  numor = numor + 1;    HELP_REF_XCALIBUR_numor        =  numor  !101
  numor = numor + 1;    HELP_RESET_numor               =  numor  !102
  numor = numor + 1;    HELP_RINT_numor                =  numor  !103
  numor = numor + 1;    HELP_RHOMB_HEX_numor           =  numor  !104
  numor = numor + 1;    HELP_SAVE_SETTINGS_numor       =  numor  !105
  numor = numor + 1;    HELP_SEARCH_EXTI_numor         =  numor  !106
  numor = numor + 1;    HELP_SEARCH_MONO_numor         =  numor  !107
  numor = numor + 1;    HELP_SEARCH_SPGR_numor         =  numor  !108
  numor = numor + 1;    HELP_SEARCH_TETRA_numor        =  numor  !109
  numor = numor + 1;    HELP_SET_numor                 =  numor  !110
  numor = numor + 1;    HELP_SETTING_numor             =  numor  !111
  numor = numor + 1;    HELP_SFAC_numor                =  numor  !112
  numor = numor + 1;    HELP_SFHKL_numor               =  numor  !113
  numor = numor + 1;    HELP_SG_numor                  =  numor  !114
  numor = numor + 1;    HELP_SG_ALL_numor              =  numor  !115
  numor = numor + 1;    HELP_SG_EXTI_numor             =  numor  !116
  numor = numor + 1;    HELP_SG_INFO_numor             =  numor  !117
  numor = numor + 1;    HELP_SG_SUB_numor              =  numor  !118
  numor = numor + 1;    HELP_SHANNON_numor             =  numor  !119
  numor = numor + 1;    HELP_SHELL_numor               =  numor  !120
  numor = numor + 1;    HELP_SHIFT_2TH_numor           =  numor  !121
  numor = numor + 1;    HELP_SITE_INFO_numor           =  numor  !122
  numor = numor + 1;    HELP_SIZE_numor                =  numor  !123
  numor = numor + 1;    HELP_SORT_numor                =  numor  !124
  numor = numor + 1;    HELP_STAR_K_numor              =  numor  !125
  numor = numor + 1;    HELP_STL_numor                 =  numor  !126
  numor = numor + 1;    HELP_SUPERCELL_numor           =  numor  !127
  numor = numor + 1;    HELP_SYMM_numor                =  numor  !128
  numor = numor + 1;    HELP_SYST_numor                =  numor  !129
  numor = numor + 1;    HELP_THERM_numor               =  numor  !130
  numor = numor + 1;    HELP_THERM_SHELX_numor         =  numor  !121
  numor = numor + 1;    HELP_THETA_numor               =  numor  !132
  numor = numor + 1;    HELP_TITL_numor                =  numor  !133
  numor = numor + 1;    HELP_TOLMAN_ANGLE_numor        =  numor  !134
  numor = numor + 1;    HELP_TRANSLATION_numor         =  numor  !135
  numor = numor + 1;    HELP_TRANSMISSION_numor        =  numor  !136
  numor = numor + 1;    HELP_TRICLINIC_numor           =  numor  !137
  numor = numor + 1;    HELP_TWIN_HEXA_numor           =  numor  !138
  numor = numor + 1;    HELP_TWIN_PSEUDO_HEXA_numor    =  numor  !139
  numor = numor + 1;    HELP_TWO_THETA_numor           =  numor  !140
  numor = numor + 1;    HELP_UB_matrix_numor           =  numor  !141
  numor = numor + 1;    HELP_UNIT_numor                =  numor  !142
  numor = numor + 1;    HELP_UPDATE_numor              =  numor  !143
  numor = numor + 1;    HELP_USER_MAT_numor            =  numor  !144
  numor = numor + 1;    HELP_WAVE_numor                =  numor  !145
  numor = numor + 1;    HELP_WEB_numor                 =  numor  !146
  numor = numor + 1;    HELP_WRITE_ADP_numor           =  numor  !147
  numor = numor + 1;    HELP_WRITE_BEAM_numor          =  numor  !148
  numor = numor + 1;    HELP_WRITE_CELL_numor          =  numor  !149
  numor = numor + 1;    HELP_WRITE_CHEM_numor          =  numor  !150
  numor = numor + 1;    HELP_WRITE_DEVICE_numor        =  numor  !151
  numor = numor + 1;    HELP_WRITE_QVEC_numor          =  numor  !152
  numor = numor + 1;    HELP_WRITE_SG_numor            =  numor  !153
  numor = numor + 1;    HELP_WRITE_SUPERCELL_numor     =  numor  !154
  numor = numor + 1;    HELP_WRITE_SYM_OP_numor        =  numor  !155
  numor = numor + 1;    HELP_WRITE_WAVE_numor          =  numor  !156
  numor = numor + 1;    HELP_WRITE_ZUNIT_numor         =  numor  !157
  numor = numor + 1;    HELP_X_wave_numor              =  numor  !158
  numor = numor + 1;    HELP_ZUNIT_numor               =  numor  !159


 end subroutine cryscalc_init
 !------------------------------------------------------

 subroutine read_cryscalc_ini()
  USE cryscalc_module, ONLY : debug_proc, INI_unit, cryscalc, winplotr_exe, my_editor, my_browser, my_word, &
                              my_pdflatex, WEB, AUTHOR, DEVICE, cryscalc, gfortran,                             &
                              wavelength, keyword_beam, keyword_WAVE, neutrons, X_rays,                         &
                              SW_DENZO, SW_EVAL, SW_APEX, SW_X2S, SW_SUPERNOVA, SW_XCALIBUR,                    &
                              CONN_dmax_ini, CONN_dmax, CONN_dmin_ini, CONN_dmin,                               &
                              CIF_format80, include_RES_file, include_HKL_file, CMD_include_HKL_file,           &
                              include_experimenter_name, CIFDEP,                                                &
                              update_parameters, report_header, report_logo, LOCK_wave_value,  Max_ref,         &
                              expert_mode, hkl_statistics, hkl_format, cartesian_frame, keep_bond_str_out,      &
                              skip_start_menu, pdp_simu, CIF_torsion_limit,                                     &
                              keyword_create_ACE, keyword_create_CEL, keyword_create_CFL, keyword_create_FST,   &
                              keyword_create_INS, keyword_create_PRF, keyword_create_PCR,                       &
                              create_CIF_PYMOL, INI_create, search_mono_criteria, ONLY_monoclinic_angle,        &
                              ONLY_monoclinic_angle_ini,                                                        &
                              structure_solution, structure_refinement, absorption_correction,                  &
                              Create_INS, get_sample_ID,                                                        &
                              message_text, winplotr_path_name, on_screen, on_screen_prf
  USE HKL_module,     ONLY : n_sig, threshold, MAX_ref_allowed, overlapping_criteria
  USE Wavelength_module
  USE USER_module
  USE CIF_module
  USE pattern_profile_module
  USE macros_module,  ONLY : u_case, l_case, verif_hkl_format, answer_yes, answer_no
  USE MATRIX_list_module
  USE IO_module,                 ONLY : write_info


  implicit none
   INTEGER                :: long, iostat_err, i, i1, i2, n_s
   LOGICAL                :: file_exist
   CHARACTER (LEN=256)    :: read_line, message2_text
   CHARACTER (len=256)    :: INI_string
   real                   :: INI_real
   real, dimension(3,3)   :: arg_mat


   ! recherche du nom du repertoire associé à CRYSCALC
   call getenv('CRYSCALC', cryscalc%path_name)
   long = len_trim(cryscalc%path_name)
   if(long /=0) then
    if(cryscalc%path_name(1:1) == '"' .and. cryscalc%path_name(long:long) == '"') then
     cryscalc%path_name = cryscalc%path_name(2:long-1)
     long = len_trim(cryscalc%path_name)
    end if
    if(cryscalc%path_name(long:long) == '\') then
     cryscalc%path_name = cryscalc%path_name(1:long-1)
    endif
    cryscalc%css        = trim(cryscalc%path_name)//'\cryscalc.css'
    cryscalc%report_css = trim(cryscalc%path_name)//'\cryscalc_report.css'
   endif

  ! recherche du fichier cryscalc.ini
  !1. dans le repertoire courant
  file_exist = .false.
  inquire(file='cryscalc.ini', exist = file_exist)
  if(file_exist) then
   cryscalc%ini = 'cryscalc.ini'
   call getenv('WINPLOTR', winplotr_path_name)
   long = len_trim(winplotr_path_name)
   if(long /=0) then
    if(winplotr_path_name(long:long) == '\') then
     winplotr_path_name = winplotr_path_name(1:long-1)
    endif
    winplotr_exe = trim(winplotr_path_name) // '\winplotr.exe'
   endif
  else

   ! >>
   ! 2. dans le repertoire ou est installe CRYSCALC
   long = len_trim(cryscalc%path_name)
   if(long /=0) then
    if(cryscalc%path_name(1:1) == '"' .and. cryscalc%path_name(long:long) == '"') then
     cryscalc%ini        = cryscalc%path_name(1:long-1) // '\cryscalc.ini"'
     cryscalc%css        = cryscalc%path_name(1:long-1) // '\cryscalc.css"'
     cryscalc%report_css = cryscalc%path_name(1:long-1) // '\cryscalc_report.css"'
    elseif(cryscalc%path_name(1:1) == "'" .and. cryscalc%path_name(long:long) == '"') then
     cryscalc%ini        = cryscalc%path_name(1:long-1) // "\cryscalc.ini'"
     cryscalc%css        = cryscalc%path_name(1:long-1) // "\cryscalc.css'"
     cryscalc%report_css = cryscalc%path_name(1:long-1) // "\cryscalc_report.css'"
    else
     cryscalc%ini        = trim(cryscalc%path_name) // '\cryscalc.ini'
     cryscalc%css        = trim(cryscalc%path_name) // '\cryscalc.css'
     cryscalc%report_css = trim(cryscalc%path_name) // '\cryscalc_report.css'
    end if
   endif

   call getenv('WINPLOTR', winplotr_path_name)
   long = len_trim(winplotr_path_name)
   if(long /=0) then
    if(winplotr_path_name(long:long) == '\') then
     winplotr_path_name = winplotr_path_name(1:long-1)
    endif
    winplotr_exe = trim(winplotr_path_name) // '\winplotr.exe'
    if(len_trim(cryscalc%ini) == 0) then
    ! 3. dans le repertoire ou est installe WinPLOTR
      cryscalc%ini  = trim(winplotr_path_name)//'\cryscalc.ini'
     !cryscalc%css = trim(cryscal%path_name) // '\cryscalc.css'
     !cryscalc%report_css = trim(cryscal%path_name) // '\cryscalc_report.css'
    endif
   endif


   if(len_trim(cryscalc%ini) == 0 .and. len_trim(winplotr_exe) == 0) then
    if(debug_proc%write) then
     call write_debug_proc('', '')
     call write_debug_proc('No environment variable defined for CRYSCALC.', '')
     call write_debug_proc('', '')
    endif
    return
   else
    if(len_trim(winplotr_exe) == 0) then
     if(debug_proc%write) then
      call write_debug_proc('', '')
      call write_debug_proc('WinPLOTR not installed or wrong WinPLOTR environment variable.', '')
      call write_debug_proc('', '')
     endif
    end if
   endif
  endif
   ! <<


   if(debug_proc%write) then
    call write_debug_proc('','')
    call write_debug_proc('CRYSCALC_ini', trim(cryscalc%ini))
    !call write_debug_proc('','')
    call write_debug_proc('CRYSCALC_css', trim(cryscalc%css))
    !call write_debug_proc('','')
    call write_debug_proc('CRYSCALC_report_css', trim(cryscalc%report_css))
    call write_debug_proc('','')
    call write_debug_proc('WinPLOTR_path', trim(winplotr_path_name))
    call write_debug_proc('WinPLOTR_exe', trim(winplotr_exe))
    call write_debug_proc('','')
   endif

   long = len_trim(cryscalc%ini)
   if(cryscalc%ini(1:1) == '"' .and. cryscalc%ini(long:long) == '"') cryscalc%ini = cryscalc%ini(2:long-1)

   INQUIRE(FILE=TRIM(cryscalc%ini), EXIST=file_exist)
   IF(.NOT. file_exist) then
    call write_info('')
    call write_info('!! '//trim(CRYSCALC%ini)//' does not exist. Default values will be used.')
    call write_info('')
    call CIF_default_values
    if(debug_proc%write) then
     call write_debug_proc('','')
     call write_debug_proc('!! '//trim(CRYSCALC%ini)//' does not exist. Default values will be used.','')
     call write_debug_proc('','')

     call write_default_values
    endif
    cryscalc%ini = ""
    return
   endif



   OPEN(UNIT=INI_unit, FILE=TRIM(cryscalc%ini))

!---------------- EDITOR, BROWSER, WORD ---------------------------------------
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    IF(LEN_TRIM(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = u_case(ADJUSTL(read_line))
    i1 = INDEX(read_line,'BROWSER')
    IF(i1==0) cycle
    i2 = INDEX(read_line,'=')
    IF(i2==0) cycle
    READ(READ_line(i2+1:), '(a)') my_browser%name
    my_browser%name = ADJUSTL(my_browser%name)
   END do
   long = len_trim(my_browser%name)
   if(long /=0)  then
    if(my_browser%name(1:1) == '"' .and. my_browser%name(long:long) == '"') then
     inquire(file=my_browser%name(2:long-1), exist = my_browser%exist)
    else
     inquire(file=trim(my_browser%name), exist = my_browser%exist)
    end if
   endif
   if(debug_proc%write) then
    call write_debug_proc('BROWSER', trim(my_browser%name))
   endif



   REWIND(UNIT=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    IF(LEN_TRIM(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = u_case(ADJUSTL(read_line))
    i1 = INDEX(read_line,'EDITOR')
    IF(i1==0) cycle
    i2 = INDEX(read_line,'=')
    IF(i2==0) cycle
    READ(READ_line(i2+1:), '(a)') my_editor%name
    my_editor%name = ADJUSTL(my_editor%name)
   END do
   long = len_trim(my_editor%name)
   if(long /=0)  then
    if(my_editor%name(1:1) == '"' .and. my_editor%name(long:long) == '"') then
     inquire(file=my_editor%name(2:long-1), exist = my_editor%exist)
    else
     inquire(file=trim(my_editor%name), exist = my_editor%exist)
    end if
   endif
   if(debug_proc%write) then
    call write_debug_proc('EDITOR', trim(my_editor%name))
   endif



   REWIND(UNIT=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    IF(LEN_TRIM(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = u_case(ADJUSTL(read_line))
    i1 = INDEX(read_line,'WORD')
    IF(i1==0) cycle
    i2 = INDEX(read_line,'=')
    IF(i2==0) cycle
    READ(READ_line(i2+1:), '(a)') my_word%name
    my_word%name= ADJUSTL(my_word%name)
   END do
   long = len_trim(my_word%name)
   if(long /=0)  then
    if(my_word%name(1:1) == '"' .and. my_word%name(long:long) == '"') then
     inquire(file=my_word%name(2:long-1), exist = my_word%exist)
    else
     inquire(file=trim(my_word%name), exist = my_word%exist)
    end if
   endif
   if(debug_proc%write) then
    call write_debug_proc('WORD', trim(my_word%name))
   endif

   REWIND(UNIT=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    IF(LEN_TRIM(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = u_case(ADJUSTL(read_line))
    i1 = INDEX(read_line,'PDFLATEX')
    IF(i1==0) cycle
    i2 = INDEX(read_line,'=')
    IF(i2==0) cycle
    READ(READ_line(i2+1:), '(a)') my_pdflatex%name
    my_pdflatex%name= ADJUSTL(my_pdflatex%name)
   END do
   long = len_trim(my_pdflatex%name)
   if(long /=0)  then
    if(my_pdflatex%name(1:1) == '"' .and. my_pdflatex%name(long:long) == '"') then
     inquire(file=my_pdflatex%name(2:long-1), exist = my_pdflatex%exist)
    else
     inquire(file=trim(my_pdflatex%name), exist = my_pdflatex%exist)
    end if
   endif
   if(debug_proc%write) then
    call write_debug_proc('PDFLATEX', trim(my_pdflatex%name))
   endif



!-------- WEB --------------------------------------------------------------------
   REWIND(UNIT=INI_unit)
   WEB%num_site = 0
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    IF(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    IF(u_case(read_line(1:4)) == '[WEB') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)          exit
      IF(LEN_TRIM(read_line) ==0) exit
      IF(read_line(1:1) == '[')   exit  ! ??
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      WEB%num_site = WEB%num_site + 1
      IF(WEB%num_site > 10) exit
      i1 = INDEX(read_line, '=')
      IF(i1 ==0) cycle
      READ(read_line(1:i1-1), '(a)', IOSTAT=iostat_err) WEB%NAME(WEB%num_site)
      IF(iostat_err /=0) then
       WEB%num_site = WEB%num_site - 1
       exit
      endif
      READ(read_line(i1+1:) , '(a)') WEB%ADDRESS(WEB%num_site)
      IF(iostat_err /=0) then
       WEB%num_site = WEB%num_site - 1
       exit
      endif

      WEB%NAME(WEB%num_site)    = u_case(WEB%NAME(WEB%num_site))
      WEB%ADDRESS(WEB%num_site) = ADJUSTL(WEB%ADDRESS(WEB%num_site))
      if(debug_proc%write) then
       call write_debug_proc ('WEB_name',    trim(WEB%NAME(WEB%num_site)))
       call write_debug_proc ('WEB_address', trim(WEB%ADDRESS(WEB%num_site)))
      endif

     end do

    endif
   END do


!--------- AUTHOR -----------------------------------------------------
   REWIND(UNIT=INI_unit)


   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    IF(read_line(1:7) == '[AUTHOR' .or. read_line(1:6) == '[USER]' .or. read_line(1:14) == '[EXPERIMENTER]') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit
      READ(read_line(i2+1:), '(a)') AUTHOR%string
      IF(l_case(read_line(1:4)) == 'name') then
       AUTHOR%name = ADJUSTL(AUTHOR%string)
       cycle
      ELSEIF(l_case(read_line(1:10))== 'first_name') then
       AUTHOR%first_name = ADJUSTL(AUTHOR%string)
       if(AUTHOR%first_name(1:1) /= "?" .and. AUTHOR%name(1:1) /= "?") then
        write(AUTHOR%initiales, '(4a)') u_case(AUTHOR%first_name(1:1)), '.', u_case(AUTHOR%name(1:1)), '.'
       end if
       cycle
      ELSEIF(l_case(read_line(1:7)) == 'address') then
       AUTHOR%address = ADJUSTL(AUTHOR%string)
       cycle
      ELSEIF(l_case(read_line(1:5)) == 'email') then
       AUTHOR%email = ADJUSTL(AUTHOR%string)
       cycle
      ELSEIF(l_case(read_line(1:3)) == 'web') then
       AUTHOR%WEB = ADJUSTL(AUTHOR%string)
       cycle
      ELSEIF(l_case(read_line(1:4)) == 'team') then
       AUTHOR%team = ADJUSTL(AUTHOR%string)
       cycle
      endif
     END do
    endif
   END do

   if(debug_proc%write) then
    call write_debug_proc ('USER_name',          trim(AUTHOR%name))
    call write_debug_proc ('USER_first_name',    trim(AUTHOR%first_name))
    call write_debug_proc ('USER_address',       trim(AUTHOR%address))
    call write_debug_proc ('USER_web',           trim(AUTHOR%web))
    call write_debug_proc ('USER_team',          trim(AUTHOR%team))
   endif


!--------------- DEVICE -------------------------------------------
   REWIND(UNIT=INI_unit)

   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    IF(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    read_line = u_case(ADJUSTL(read_line))
    IF(read_line(1:7) == '[DEVICE') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit
      READ(read_line(i2+1:), '(a)') DEVICE%string
      IF(l_case(read_line(1:9)) == 'diffracto') then
       DEVICE%diffracto = ADJUSTL(DEVICE%string)
       cycle
      ELSEIF(l_case(read_line(1:3))== 'lab') then
       DEVICE%lab = ADJUSTL(DEVICE%string)
       cycle
      ELSEIF(l_case(read_line(1:9)) == 'radiation') then
       DEVICE%radiation = ADJUSTL(DEVICE%string)
       cycle
      ELSEIF(l_case(read_line(1:4)) == 'wave') then
       DEVICE%wave = ADJUSTL(DEVICE%string)
       read(DEVICE%wave, *) wavelength
       CIF_cell_measurement%wavelength = wavelength
       keyword_WAVE = .true.
       cycle
      ELSEIF(l_case(read_line(1:11)) == 'temperature') then
       DEVICE%string = ADJUSTL(DEVICE%string)
       read(DEVICE%string, *) CIF_cell_measurement%temperature
       cycle
      endif
     END do
    endif
   END do

   anti_cathode = .false.

   if(u_case(DEVICE%radiation(1:1)) == 'X') then
    neutrons = .false.
    X_rays   = .true.
    if (DEVICE%wave == '?') call get_X_radiation(u_case(DEVICE%radiation))
    keyword_beam = .true.
   ELSEIF(u_case(DEVICE%radiation(1:)) == 'NEUT') then
    neutrons = .true.
    X_rays   = .false.
    keyword_beam = .true.
   endif


   ! -------------------- sept 2010  ---------------------------------
   if(u_case(DEVICE%radiation(1:4)) == 'X_Mo') then
    CIF_parameter_DEVICE%diffracto_radiation_type = 'MoK\a'
   elseif(u_case(DEVICE%radiation(1:4)) == 'X_Cu') then
    CIF_parameter_DEVICE%diffracto_radiation_type = 'CuK\a'
   endif
   CIF_parameter_DEVICE%diffracto_radiation_wavelength = DEVICE%wave
   ! ------------------------------------------------------------------


   if(u_case(DEVICE%diffracto(1:4)) == 'APEX') then
    DEVICE%APEX = .true.
    CIF_parameter_DEVICE = CIF_parameter_DEVICE_APEX
    !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_APEX%diffrn_measurement_device_type
    !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_APEX%diffrn_measurement_method
    !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_APEX%diffrn_detector
    !CIF_parameter_DEVICE%diffrn_radiation_wavelength          = CIF_parameter_DEVICE_APEX%diffrn_radiation_wavelength
    !CIF_parameter_DEVICE%diffrn_radiation_type                = CIF_parameter_DEVICE_APEX%diffrn_radiation_type
    !CIF_parameter_DEVICE%diffrn_radiation_source              = CIF_parameter_DEVICE_APEX%diffrn_radiation_source
    !CIF_parameter_DEVICE%diffrn_radiation_monochromator       = CIF_parameter_DEVICE_APEX%diffrn_radiation_monochromator
    !CIF_parameter_DEVICE%diffrn_radiation_probe               = CIF_parameter_DEVICE_APEX%diffrn_radiation_probe
    !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_APEX%computing_data_collection
    !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_APEX%computing_cell_refinement
    !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_APEX%computing_data_reduction
    !!CIF_parameter_DEVICE%computing_structure_solution         = CIF_parameter_DEVICE_APEX%computing_structure_solution
    !!CIF_parameter_DEVICE%computing_structure_refinement       = CIF_parameter_DEVICE_APEX%computing_structure_refinement
    !!CIF_parameter_DEVICE%computing_molecular_graphics         = CIF_parameter_DEVICE_APEX%computing_molecular_graphics
    !!CIF_parameter_DEVICE%computing_publication_material       = CIF_parameter_DEVICE_APEX%computing_publication_material

   elseif(u_case(DEVICE%diffracto(1:13)) == 'D8_VENTURE_CU' .or. u_case(DEVICE%diffracto(1:6))  == 'D8V_CU' .or. &
          u_case(DEVICE%diffracto(1:13)) == 'D8 VENTURE CU' .or. u_case(DEVICE%diffracto(1:6))  == 'D8V CU') then
    DEVICE%D8V   = .true.
    DEVICE%D8VCu = .true.
    DEVICE%radiation  =  'X_Cu'
    if (DEVICE%wave == '?') call get_X_radiation(u_case(DEVICE%radiation))
    CIF_parameter_DEVICE = CIF_parameter_DEVICE_D8V_CU
    !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_D8V_Cu%diffrn_measurement_device_type
    !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_D8V_Cu%diffrn_measurement_method
    !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_D8V_Cu%diffrn_detector
    !CIF_parameter_DEVICE%diffrn_radiation_wavelength          = CIF_parameter_DEVICE_D8V_Cu%diffrn_radiation_wavelength
    !CIF_parameter_DEVICE%diffrn_radiation_type                = CIF_parameter_DEVICE_D8V_Cu%diffrn_radiation_type
    !CIF_parameter_DEVICE%diffrn_radiation_source              = CIF_parameter_DEVICE_D8V_Cu%diffrn_radiation_source
    !CIF_parameter_DEVICE%diffrn_radiation_monochromator       = CIF_parameter_DEVICE_D8V_Cu%diffrn_radiation_monochromator
    !CIF_parameter_DEVICE%diffrn_radiation_probe               = CIF_parameter_DEVICE_D8V_Cu%diffrn_radiation_probe
    !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_D8V_Cu%computing_data_collection
    !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_D8V_Cu%computing_cell_refinement
    !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_D8V_Cu%computing_data_reduction

   elseif(u_case(DEVICE%diffracto(1:13)) == 'D8_VENTURE_MO' .or. u_case(DEVICE%diffracto(1:6))  == 'D8V_MO' .or. &
          u_case(DEVICE%diffracto(1:13)) == 'D8 VENTURE MO' .or. u_case(DEVICE%diffracto(1:6))  == 'D8V MO') then
    DEVICE%D8V   = .true.
    DEVICE%D8VMo = .true.
    DEVICE%radiation =  'X_Mo'
    if (DEVICE%wave == '?') call get_X_radiation(u_case(DEVICE%radiation))
    CIF_parameter_DEVICE = CIF_parameter_DEVICE_D8V_MO
    !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_D8V_Mo%diffrn_measurement_device_type
    !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_D8V_Mo%diffrn_measurement_method
    !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_D8V_Mo%diffrn_detector
    !CIF_parameter_DEVICE%diffrn_radiation_wavelength          = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_wavelength
    !CIF_parameter_DEVICE%diffrn_radiation_type                = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_type
    !CIF_parameter_DEVICE%diffrn_radiation_source              = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_source
    !CIF_parameter_DEVICE%diffrn_radiation_monochromator       = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_monochromator
    !CIF_parameter_DEVICE%diffrn_radiation_probe               = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_probe
    !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_D8V_Mo%computing_data_collection
    !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_D8V_Mo%computing_cell_refinement
    !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_D8V_Mo%computing_data_reduction

   elseif((len_trim(DEVICE%diffracto) == 10 .and. u_case(DEVICE%diffracto(1:10)) == 'D8_VENTURE') .or. &
          (len_trim(DEVICE%diffracto) == 10 .and. u_case(DEVICE%diffracto(1:10)) == 'D8 VENTURE') .or. &
          (len_trim(DEVICE%diffracto) ==  3 .and. u_case(DEVICE%diffracto(1:3))  == 'D8V'))         then
    DEVICE%D8V   = .true.
    CIF_parameter_DEVICE = CIF_parameter_DEVICE_D8V
    !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_D8V_Mo%diffrn_measurement_device_type
    !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_D8V_Mo%diffrn_measurement_method
    !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_D8V_Mo%diffrn_detector
    !CIF_parameter_DEVICE%diffrn_radiation_type                = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_type
    !CIF_parameter_DEVICE%diffrn_radiation_source              = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_source
    !CIF_parameter_DEVICE%diffrn_radiation_monochromator       = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_monochromator
    !CIF_parameter_DEVICE%diffrn_radiation_probe               = CIF_parameter_DEVICE_D8V_Mo%diffrn_radiation_probe
    !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_D8V_Mo%computing_data_collection
    !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_D8V_Mo%computing_cell_refinement
    !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_D8V_Mo%computing_data_reduction

   elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
    CIF_parameter_DEVICE = CIF_parameter_DEVICE_KCCD
    DEVICE%KCCD = .true.
    !CIF_parameter_DEVICE%diffrn_source                        = CIF_parameter_DEVICE_KCCD%diffrn_source
    !CIF_parameter_DEVICE%diffrn_radiation_wavelength          = CIF_parameter_DEVICE_KCCD%diffrn_radiation_wavelength
    !CIF_parameter_DEVICE%diffrn_radiation_type                = CIF_parameter_DEVICE_KCCD%diffrn_radiation_type
    !CIf_parameter_DEVICE%diffrn_radiation_source              = CIF_parameter_DEVICE_KCCD%diffrn_radiation_source
    !CIF_parameter_DEVICE%diffrn_radiation_monochromator       = CIF_parameter_DEVICE_KCCD%diffrn_radiation_monochromator
    !CIF_parameter_DEVICE%diffrn_radiation_probe               = CIF_parameter_DEVICE_KCCD%diffrn_radiation_probe
    !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_KCCD%diffrn_detector
    !CIF_parameter_DEVICE%diffrn_detector_area_resol_mean      = CIF_parameter_DEVICE_KCCD%diffrn_detector_area_resol_mean
    !CIF_parameter_DEVICE%diffrn_measurement_device            = CIF_parameter_DEVICE_KCCD%diffrn_measurement_device
    !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_KCCD%diffrn_measurement_device_type
    !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_KCCD%diffrn_measurement_method
    !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_KCCD%computing_data_collection
    !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_KCCD%computing_cell_refinement
    !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_KCCD%computing_data_reduction

    !SW_DENZO%data_collection                              = CIF_parameter_DEVICE_KCCD%computing_data_collection
    !SW_DENZO%cell_refinement                              = CIF_parameter_DEVICE_KCCD%computing_cell_refinement
    !SW_DENZO%data_reduction                               = CIF_parameter_DEVICE_KCCD%computing_data_reduction

    !CIF_parameter_DEVICE%computing_data_collection            = "'Collect (Nonius BV, 1997-2000)'"
    !CIF_parameter_DEVICE%computing_cell_refinement            = "'Dirax/lsq (Duisenberg & Schreurs, 1989-2000)'"
    !CIF_parameter_DEVICE%computing_data_reduction             = "'EvalCCD (Duisenberg & Schreurs 1990-2000)'"
    !SW_EVAL%data_collection                               = CIF_parameter_DEVICE%computing_data_collection
    !SW_EVAL%cell_refinement                               = CIF_parameter_DEVICE%computing_cell_refinement
    !SW_EVAL%data_reduction                                = CIF_parameter_DEVICE%computing_data_reduction


   elseif(u_case(DEVICE%diffracto(1:8)) == 'SAPPHIRE' .or. u_case(DEVICE%diffracto(1:8)) == 'XCALIBUR') then
    DEVICE%Xcalibur = .true.
    CIF_parameter_DEVICE = CIF_parameter_DEVICE_XCALIBUR
    !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_XCALIBUR%diffrn_measurement_device_type
    !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_XCALIBUR%diffrn_measurement_method
    !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_XCALIBUR%diffrn_detector
    !CIF_parameter_DEVICE%diffrn_detector_area_resol_mean      = CIF_parameter_DEVICE_XCALIBUR%diffrn_detector_area_resol_mean
    !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_XCALIBUR%computing_data_collection
    !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_XCALIBUR%computing_cell_refinement
    !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_XCALIBUR%computing_data_reduction

    SW_XCALIBUR%data_collection                           = CIF_parameter_DEVICE%computing_data_collection
    SW_XCALIBUR%cell_refinement                           = CIF_parameter_DEVICE%computing_cell_refinement
    SW_XCALIBUR%data_reduction                            = CIF_parameter_DEVICE%computing_data_reduction

   elseif(u_case(DEVICE%diffracto(1:9)) == 'SUPERNOVA' ) then
    DEVICE%SuperNova = .true.
    CIF_parameter_DEVICE = CIF_parameter_DEVICE_XCALIBUR
    !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_SUPERNOVA%diffrn_measurement_device_type
    !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_SUPERNOVA%diffrn_measurement_method
    !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_SUPERNOVA%diffrn_detector
    !CIF_parameter_DEVICE%diffrn_detector_area_resol_mean      = CIF_parameter_DEVICE_SUPERNOVA%diffrn_detector_area_resol_mean
    !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_SUPERNOVA%computing_data_collection
    !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_SUPERNOVA%computing_cell_refinement
    !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_SUPERNOVA%computing_data_reduction

    SW_SUPERNOVA%data_collection                          = CIF_parameter_DEVICE%computing_data_collection
    SW_SUPERNOVA%cell_refinement                          = CIF_parameter_DEVICE%computing_cell_refinement
    SW_SUPERNOVA%data_reduction                           = CIF_parameter_DEVICE%computing_data_reduction


   elseif(u_case(DEVICE%diffracto(1:3)) == 'X2S') then
    DEVICE%X2S = .true.
    CIF_parameter_DEVICE = CIF_parameter_DEVICE_APEX
    !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_X2S%diffrn_measurement_device_type
    !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_X2S%diffrn_measurement_method
    !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_X2S%diffrn_detector
    !CIF_parameter_DEVICE%diffrn_radiation_wavelength          = CIF_parameter_DEVICE_X2S%diffrn_radiation_wavelength
    !CIF_parameter_DEVICE%diffrn_radiation_type                = CIF_parameter_DEVICE_X2S%diffrn_radiation_type
    !CIF_parameter_DEVICE%diffrn_radiation_source              = CIF_parameter_DEVICE_X2S%diffrn_radiation_source
    !CIF_parameter_DEVICE%diffrn_radiation_monochromator       = CIF_parameter_DEVICE_X2S%diffrn_radiation_monochromator
    !CIF_parameter_DEVICE%diffrn_radiation_probe               = CIF_parameter_DEVICE_X2S%diffrn_radiation_probe
    !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_X2S%computing_data_collection
    !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_X2S%computing_cell_refinement
    !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_X2S%computing_data_reduction

   endif




   if(debug_proc%write) then
    call write_debug_proc ('DEVICE_diffractometer', trim(DEVICE%diffracto))
    call write_debug_proc ('DEVICE_lab',            trim(DEVICE%lab))
    call write_debug_proc ('DEVICE_radiation',      trim(DEVICE%radiation))
    call write_debug_proc ('DEVICE_wave',           trim(DEVICE%wave))
    !if(u_case(DEVICE%diffracto(1:4))  == 'APEX'          .or.  &
    !   u_case(DEVICE%diffracto(1:13)) == 'D8_VENTURE_CU' .or.  &
    !   u_case(DEVICE%diffracto(1:13)) == 'D8 VENTURE CU' .or.  &
    !   u_case(DEVICE%diffracto(1:6))  == 'D8V_CU'        .or.  &
    !   u_case(DEVICE%diffracto(1:6))  == 'D8V CU'        .or.  &
    !   u_case(DEVICE%diffracto(1:13)) == 'D8_VENTURE_MO' .or.  &
    !   u_case(DEVICE%diffracto(1:13)) == 'D8 VENTURE MO' .or.  &
    !   u_case(DEVICE%diffracto(1:6))  == 'D8V_MO'        .or.  &
    !   u_case(DEVICE%diffracto(1:6))  == 'D8V MO') then
    if(DEVICE%APEX .or. DEVICE%D8V) then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call WRITE_debug_proc('MEASUREMENT_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter_DEVICE%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter_DEVICE%diffrn_radiation_probe))

    !elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
    elseif(DEVICE%KCCD) then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call WRITE_debug_proc('MEASUREMENT_device',        trim(CIF_parameter_DEVICE%diffrn_measurement_device))
     call WRITE_debug_proc('DIFFRACTION_source',        trim(CIF_parameter_DEVICE%diffrn_source))
     call WRITE_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call WRITE_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter_DEVICE%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter_DEVICE%diffrn_radiation_probe))

    !elseif(u_case(DEVICE%diffracto(1:7)) == 'SAPHIRE' .or. u_case(DEVICE%diffracto(1:9)) == 'EXCALIBUR') then
    elseif(DEVICE%Xcalibur) then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call WRITE_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call WRITE_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter_DEVICE%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter_DEVICE%diffrn_radiation_probe))

    !elseif(u_case(DEVICE%diffracto(1:9)) == 'SUPERNOVA') then
    elseif(DEVICE%SuperNOVA) then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call WRITE_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call WRITE_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter_DEVICE%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter_DEVICE%diffrn_radiation_probe))

    !elseif(u_case(DEVICE%diffracto(1:3)) == 'X2S') then
    elseif(DEVICE%X2S) then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call WRITE_debug_proc('MEASUREMENT_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter_DEVICE%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter_DEVICE%diffrn_radiation_probe))


    endif
    call WRITE_debug_proc('COMPUTING_data_collection',  trim(CIF_parameter_DEVICE%computing_data_collection))
    call WRITE_debug_proc('COMPUTING_cell_refinement',  trim(CIF_parameter_DEVICE%computing_cell_refinement))
    call WRITE_debug_proc('COMPUTING_data_reduction',   trim(CIF_parameter_DEVICE%computing_data_reduction))

   endif



 !------------ PARAMETERS -------------------------------------------
   rewind(unit=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    IF(LEN_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    IF(read_line(1:12) == '[PARAMETERS]') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit
      READ(read_line(i2+1:), '(a)') INI_string
      IF(l_case(read_line(1:5)) == 'i_sig') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       n_sig = INI_real
       cycle
      elseif(l_case(read_line(1:9)) == 'threshold') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       threshold = INI_real
       cycle
      elseif(l_case(read_line(1:7)) == 'd_min_a') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       CONN_dmin_ini = INI_real
       CONN_dmin     = CONN_dmin_ini
       cycle
      elseif(l_case(read_line(1:7)) == 'd_max_a') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       CONN_dmax_ini = INI_real
       CONN_dmax     = CONN_dmax_ini
       cycle
      endif
     end do
    endif
   end do
   if(debug_proc%write) then
    write(read_line, '(F5.2,a)') n_sig,          '    ! sigma_I (used in the SEARCH_GROUP routine)'
    call write_debug_proc ('SEARCH_I_sig',      trim(read_line))
    write(read_line, '(F6.3,a)') threshold,      '   ! I > threshold * I_max (used in the SEARCH_GROUP routine)'
    call write_debug_proc ('SEARCH_threshold',  trim(read_line))
    write(read_line, '(F6.2,a)') CONN_dmax_ini,  '    ! (in A) for connectivity calculations '
    call write_debug_proc ('CONNECTIVITY_dmax', trim(read_line))
   endif

   !-------- ARRAYS DIMENSIONS    ---------------------------------------------------------------------
  rewind(unit = INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT = iostat_err) read_line
    if(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle

    if(read_line(1:19) == '[ARRAYS DIMENSIONS]') then
     do
      read(INI_unit, '(a)', iostat = iostat_err) read_line
      if(iostat_err /=0)          exit
      if(len_trim(read_line)==0)  exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = adjustl(read_line)
      i2 = index(read_line, '=')
      if(i2==0) exit

      read(read_line(i2+1:), '(a)') INI_string
      INI_string = adjustl(INI_string)
      if(l_case(read_line(1:15)) == 'hkl_reflections') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       if(INI_real < MAX_ref_allowed) Max_ref = INI_real
       cycle
      end if
     end do
    endif
   end do

   if(debug_proc%write) then
    write(read_line, '(I10,a)') Max_ref,      '    ! Max. number of hkl reflections in a file'
    call write_debug_proc ('Max_ref',  trim(read_line))
   endif



!-------- CREATE INS parameters     ---------------------------------------------------------------------
  rewind(unit = INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT = iostat_err) read_line
    if(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle

    if(read_line(1:12) == '[CREATE INS]') then
     do
      read(INI_unit, '(a)', iostat = iostat_err) read_line
      if(iostat_err /=0)          exit
      if(len_trim(read_line)==0)  exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = adjustl(read_line)
      i2 = index(read_line, '=')
      if(i2==0) exit

      read(read_line(i2+1:), '(a)') INI_string
      INI_string = adjustl(INI_string)
      if(l_case(read_line(1:13)) == 'get_sample_id' .and.  answer_yes(trim(INI_string))) then
       get_sample_ID = .true.
      elseif(l_case(read_line(1:11)) == 'temperature') then
       long = len_trim(INI_string)
       if(u_case(INI_string(long:long)) == 'K') then
        read(INI_string(1:long-1), *, iostat=iostat_err) INI_real
        if(iostat_err /=0) exit
        Create_INS%temperature = INI_real - 273
       else
        read(INI_string, *, iostat = iostat_err) INI_real
        if(iostat_err /=0) exit
        Create_INS%temperature = INI_real
       endif
       cycle
      elseif(l_case(read_line(1:11)) == 'u_threshold') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       Create_INS%U_threshold = INI_real
       cycle

      end if
     end do
    endif
   end do

   if(debug_proc%write) then
    if(get_sample_ID) then
     call write_debug_proc ('get sample_ID',  '             1     ! get sample_ID')
    else
     call write_debug_proc ('get sample_ID',  '             0     ! sample_ID = job')
    endif
    write(read_line, '(F8.2,a)') Create_INS%temperature,      '  ! experimental temperature'
    call write_debug_proc ('temperature',      trim(read_line))
    write(read_line, '(F8.2,a)') Create_INS%U_threshold,      '    ! atoms with Uiso > u_threshold are excluded'
    call write_debug_proc ('U_threshold',  trim(read_line))
   endif

!-------- CIF COMMAND LINE ARGUMENTS --------------------------------------------------------------------
   rewind(unit=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    IF(read_line(1:24) == '[COMMAND LINE ARGUMENTS]') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit

      READ(read_line(i2+1:), '(a)') INI_string
      INI_string = adjustl(INI_string)
      IF(l_case(read_line(1:10)) == 'create_ace' .and.  answer_yes(trim(INI_string))) then
       INI_create%ACE     = .true.
       keyword_create_ACE = .true.
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_cel' .and. answer_yes(trim(INI_string))) then
       INI_create%CEL     = .true.
       keyword_create_CEL = .true.
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_cfl' .and.answer_yes(trim(INI_string))) then
       INI_create%CFL     = .true.
       keyword_create_CFL = .true.
       cycle
      ELSEIF(L_case(read_line(1:10)) == 'create_fst' .AND. answer_yes(trim(INI_string))) then
       INI_create%FST     = .true.
       keyword_create_FST = .true.
       cycle
      ELSEIF(L_case(read_line(1:16)) == 'create_cif_pymol' .AND. answer_yes(trim(INI_string))) then
       INI_create%CIF_PYMOL      = .true.
       create_CIF_PYMOL  = .true.
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_ins' .and. answer_yes(trim(INI_string))) then
       INI_create%INS     = .true.
       keyword_create_INS = .true.
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_pcr' .and. answer_yes(trim(INI_string))) then
       INI_create%PCR     = .true.
       keyword_create_PCR = .true.
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'create_pat_prf' .and. answer_yes(trim(INI_string))) then
       INI_create%PRF     = .true.
       keyword_create_PRF = .true.
       if(l_case(read_line(15:16)) == '_1' .or. l_case(read_line(15:18)) == '_out' ) on_screen_prf = .true.
       cycle
      endif
     end do
    endif
   end do

   if(debug_proc%write) then
    if (keyword_create_ACE) then
     call write_debug_proc ('create_ACE',      '1     ! .ACE file for Carine')
    else
     call write_debug_proc ('create_ACE',      '0     ! .ACE file for Carine')
    endif
    if (keyword_create_CEL) then
     call write_debug_proc ('create_CEL',      '1     ! .CEL file for PowderCELL')
    else
     call write_debug_proc ('create_CEL',      '0     ! .CEL file for PowderCELL')
    endif
    if (keyword_create_CFL) then
     call write_debug_proc ('create_CFL',      '1     ! .CFL file for CRYSCALC')
    else
     call write_debug_proc ('create_CFL',      '0     ! .CFL file for CRYSCALC')
    endif
    if (keyword_create_FST) then
     call write_debug_proc ('create_FST',      '1     ! .FST file for FP Studio')
    else
     call write_debug_proc ('create_FST',      '0     ! .FST file for FP Studio')
    endif
    if (create_CIF_PYMOL) then
     call write_debug_proc ('create_CIF_PYMOL',    '1     ! .CIF file for PYMOL')
    else
     call write_debug_proc ('create_CIF_PYMOL',    '0     ! .CIF file for PYMOL')
    endif
    if (keyword_create_INS) then
     call write_debug_proc ('create_INS',      '1     ! .INS file for SHELXL')
    else
     call write_debug_proc ('create_INS',      '0     ! .INS file for SHELXL')
    endif
    if (keyword_create_PRF) then
     call write_debug_proc ('create_PAT_PRF',  '1     ! .PRF file for FullProf')
    else
     call write_debug_proc ('create_PAT_PRF',  '0     ! .PRF file for FullProf')
    endif
    if (keyword_create_PCR) then
     call write_debug_proc ('create_PCR',  '1     ! .PCR file for FullProf')
    else
     call write_debug_proc ('create_PCR',  '0     ! .PCR file for FullProf')
    endif
   endif


   !-------- OPTIONS --------------------------------------------------------------------
   rewind(unit=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    IF(read_line(1:9) == '[OPTIONS]') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit

      READ(read_line(i2+1:), '(a)') INI_string
      INI_string = adjustl(INI_string)
      IF(l_case(read_line(1:15)) == 'lock_wave_value') then
       read(INI_string, *) LOCK_wave_value
       cycle
      ELSEIF(l_case(read_line(1:12)) == 'cif_format80' .and. answer_no(trim(INI_string))) then
       CIF_format80     = .false.
       cycle
      ELSEIF(l_case(read_line(1:16)) == 'include_res_file' .and. answer_yes(trim(INI_string))) then
       include_RES_file = .true.
       cycle
      ELSEIF(l_case(read_line(1:16)) == 'include_hkl_file' .and. answer_yes(trim(INI_string))) then
       if(CMD_include_HKL_file==0) include_HKL_file = .true.
       cycle
      ELSEIF(l_case(read_line(1:20)) == 'include_experimenter' .and. answer_yes(trim(INI_string))) then
       include_experimenter_name = .true.
       cycle
      ELSEIF(l_case(read_line(1:17)) == 'update_parameters' .and. answer_no(trim(INI_string))) then
       update_parameters = .false.
       cycle
      ELSEIF((l_case(read_line(1:18)) == 'html_report_header' .or. l_case(read_line(1:13)) == 'report_header') &
             .and. answer_no(trim(INI_string))) then
       report_header = .false.
       cycle
      ELSEIF(l_case(read_line(1:13)) == 'report_logo_1') then
       report_logo(1) = INI_string
       cycle
      ELSEIF(l_case(read_line(1:13)) == 'report_logo_2') then
       report_logo(2) = INI_string
       cycle
      elseif(l_case(read_line(1:11)) == 'expert_mode'    .and. answer_yes(trim(INI_string))) then
       expert_mode = .true.
       cycle
      elseif(l_case(read_line(1:12)) == 'bond_str.out'    .and. answer_yes(trim(INI_string))) then
       if (expert_mode) keep_bond_str_out = .true.
       cycle
      elseif(l_case(read_line(1:5)) == 'debug' .and. answer_yes(trim(ini_string))) then
        if(l_case(read_line(1:13)) == 'debug_level_3'      .and. answer_yes(trim(ini_string))) then
        debug_proc%level_3 = .true.
        debug_proc%level_2 = .true.
       elseif(l_case(read_line(1:13)) == 'debug_level_2'  .and. answer_yes(trim(ini_string))) then
        debug_proc%level_2 = .true.
       end if
       debug_proc%write = .true.
       !call open_debug_file
       cycle

      !elseif(l_case(read_line(1:13)) == 'debug_level_3'  .and. answer_yes(trim(ini_string))) then
      ! debug_proc%level_2 = .true.
      ! debug_proc%level_3 = .true.
      ! debug_proc%write   = .true.
      ! call open_debug_file
      ! cycle
      !elseif(l_case(read_line(1:13)) == 'debug_level_2'  .and. answer_yes(trim(ini_string))) then
       ! debug_proc%level_2 = .true.
      ! debug_proc%level_3 = .false.
      ! debug_proc%write   = .true.
      ! call open_debug_file
      ! cycle

      elseif(l_case(read_line(1:14)) == 'hkl_statistics' .and. answer_no(trim(INI_string))) then
       hkl_statistics = .false.
       cycle
      elseif(l_case(read_line(1:10)) == 'hkl_format') then
       !call verif_format(trim(INI_string))
       !hkl_format = trim(INI_string)
       hkl_format = verif_hkl_format(trim(INI_string))
       cycle

      elseif(l_case(read_line(1:15)) == 'skip_start_menu' .and. answer_no(trim(INI_string))) then
       skip_start_menu = .false.
       !if(gfortran) then
       ! call write_info(' Press <Enter> key to continue.')
       ! read(*,*)
       !else
       ! pause
       !end if
       cycle
      elseif(l_case(read_line(1:6)) == 'pdp_cu' .and. answer_no(trim(INI_string))) then
       pdp_simu%cu = .false.
       cycle
      elseif(l_case(read_line(1:8)) == 'pdp_beam') then
        if(ini_string(1:1) == "n" .or. ini_string(1:1) == "N") then
         pdp_simu%beam = "neutrons"
        end if
        cycle
      elseif(l_case(read_line(1:8)) == 'pdp_wave') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err /=0) exit
       pdp_simu%wave = INI_real
       cycle
      elseif(l_case(read_line(1:20)) == "cartesian_frame%type") then
       if(ini_string(1:1) == "c" .or. ini_string(1:1) == "C") then
        cartesian_frame%type   = "C"
        cartesian_frame%string = "x // c"
       end if
       cycle
      elseif(l_case(read_line(1:17)) == 'cif_torsion_limit') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       if (ABS(INI_real) < 180. ) CIF_torsion_limit =  INI_real
       cycle
      elseif(l_case(read_line(1:20)) == 'ref_overlap_criteria') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       if (ABS(INI_real) < 0.25 ) overlapping_criteria =  INI_real
       cycle
      elseif(l_case(read_line(1:20)) == 'search_mono_criteria') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       if (ABS(INI_real) < 2.5 ) search_mono_criteria(1) =  INI_real
       cycle
      elseif(l_case(read_line(1:14)) == 'search_mono_sg' .and. .not. answer_no(trim(INI_string))) then
       ONLY_monoclinic_angle_ini = .true.
       ONLY_monoclinic_angle     = .true.
       cycle
      elseif(l_case(read_line(1:19)) == 'search_sg_only_mono' .and. .not. answer_no(trim(INI_string))) then
       ONLY_monoclinic_angle_ini = .true.
       ONLY_monoclinic_angle     = .true.
       cycle

      endif


     end do
    endif
   end do

   if(debug_proc%write) then
    write(message2_text, '(F6.3,a)') LOCK_wave_value, '   ! '
    call write_debug_proc ('LOCK_wave_value', trim(message2_text))
    if (CIF_format80) then
     call write_debug_proc ('CIF_format80',         '1     ! ')
    else
     call write_debug_proc ('CIF_format80',         '0     ! ')
    endif
    write(message2_text, '(F6.2,a)') CIF_torsion_limit, '   ! '
    call write_debug_proc ('CIF_torsion_limit', trim(message2_text))
    if (include_res_file) then
     call write_debug_proc ('Include_RES_file',      '1     ! ')
    else
     call write_debug_proc ('Include_RES_file',      '0     ! ')
    endif
    if (include_hkl_file) then
     call write_debug_proc ('Include_HKL_file',      '1     ! ')
    else
     call write_debug_proc ('Include_HKL_file',      '0     ! ')
    endif
    if (update_parameters) then
     call write_debug_proc ('Update_parameters',     '1     ! ')
    else
     call write_debug_proc ('Update_parameters',     '0     ! ')
    endif
    if (report_header) then
     call write_debug_proc ('Report_header'    ,    '1     ! ')
    else
     call write_debug_proc ('Report_header'    ,    '0     ! ')
    endif
    write(message2_text, '(2a)') trim(report_logo(1)), '   ! '
    call write_debug_proc ('Report_logo_1', trim(message2_text))
    write(message2_text, '(2a)') trim(report_logo(2)), '   ! '
    call write_debug_proc ('Report_logo_2', trim(message2_text))

    if(expert_mode) then
     call write_debug_proc('Expert_mode'       ,    '1     ! ')
     if(keep_bond_str_out) then
     call write_debug_proc('Keep bond_str.out' ,    '1     ! ')
     else
     call write_debug_proc('Keep bond_str.out' ,    '0     ! ')
     endif
    else
     call write_debug_proc('Expert_mode'       ,    '0     ! ')
    endif
    if(debug_proc%level_3) then
     call write_debug_proc('Debug_mod_level_3' ,    '1     ! ')
    elseif(debug_proc%level_2) then
     call write_debug_proc('Debug_mod_level_2' ,    '1     ! ')
    endif
    if(skip_start_menu) then
     call write_debug_proc('Skip_start_menu'   ,     '1     !')
    else
     call write_debug_proc('Skip_start_menu'   ,     '0     !')
    end if
    if(hkl_statistics) then
     call write_debug_proc('hkl_statistics'     ,    '1     ! ')
    else
     call write_debug_proc('hkl_statistics'     ,    '0     ! ')
    endif
    call write_debug_proc('hkl_format'          ,     trim(hkl_format)//'  !')
    call write_debug_proc('cartesian_frame_type',     cartesian_frame%type// ' '//trim(cartesian_frame%string))
    if(expert_mode) then
     if(pdp_simu%cu) then
      call write_debug_proc('Ka1 Cu for pdp simulation'     ,    '1     ! ')
     else
      call write_debug_proc('Ka1 Cu for pdp simulation'     ,    '0     ! ')
     endif
     if(pdp_simu%beam == 'neutrons') then
      call write_debug_proc('Pattern simulation',  'N  ! Neutrons')
     else
      call write_debug_proc('Pattern simulation',  'X  ! X-rays')
     endif
     write(message2_text, '(F10.5,a)') pdp_simu%wave, '   ! wavelength'
     call write_debug_proc ('Pattern simulation', trim(message2_text))
    endif
    write(message2_text, '(F6.3,a)') overlapping_criteria, '   ! '
    call write_debug_proc ('ref_overlap_criteria', trim(message2_text))
    write(message2_text, '(F6.3,a)') search_mono_criteria(1), '   ! '
    call write_debug_proc ('search_mono_criteria', trim(message2_text))

    if(ONLY_monoclinic_angle) then
     call write_debug_proc('search_SG_only_mono'     ,    '1     ! ')
    else
     call write_debug_proc('search_SG_only_mono'     ,    '0     ! ')
    endif

   endif

      !-------- ARCHIVE_AND_REPORT --------------------------------------------------------------------
   rewind(unit=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    IF(read_line(1:20) == '[ARCHIVE_AND_REPORT]') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit

      READ(read_line(i2+1:), '(a)') INI_string
      INI_string = adjustl(INI_string)
      IF(l_case(read_line(1:12)) == 'cif_format80' .and. answer_no(trim(INI_string))) then
       CIF_format80     = .false.
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'cif_author' .and. answer_yes(trim(INI_string))) then
       CIFDEP     = .true.
       cycle
      ELSEIF(l_case(read_line(1:17)) == 'cif_torsion_limit') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       if (ABS(INI_real) < 180. ) CIF_torsion_limit =  INI_real
       cycle
      ELSEIF(l_case(read_line(1:16)) == 'include_res_file' .and. answer_yes(trim(INI_string))) then
       include_RES_file = .true.
       cycle
      ELSEIF(l_case(read_line(1:16)) == 'include_hkl_file' .and. answer_yes(trim(INI_string))) then
       if(CMD_include_HKL_file==0) include_HKL_file = .true.
       cycle
      ELSEIF(l_case(read_line(1:20)) == 'include_experimenter' .and. answer_yes(trim(INI_string))) then
       include_experimenter_name = .true.
       cycle
      ELSEIF((l_case(read_line(1:18)) == 'html_report_header' .or. l_case(read_line(1:13)) == 'report_header') &
             .and. answer_no(trim(INI_string))) then
       report_header = .false.
       cycle
      ELSEIF(l_case(read_line(1:13)) == 'report_logo_1') then
       report_logo(1) = INI_string
       cycle
      ELSEIF(l_case(read_line(1:13)) == 'report_logo_2') then
       report_logo(2) = INI_string
       cycle
      endif
     end do
    endif
   end do

   if(debug_proc%write) then
    if (CIF_format80) then
     call write_debug_proc ('CIF_format80',         '1     ! ')
    else
     call write_debug_proc ('CIF_format80',         '0     ! ')
    endif
    write(message2_text, '(F6.2,a)') CIF_torsion_limit, '   ! '
    call write_debug_proc ('CIF_torsion_limit', trim(message2_text))

    if (CIFDEP) then
     call write_debug_proc ('CIF_author',           '1     ! ')
    else
     call write_debug_proc ('CIF_author',           '0     ! ')
    endif
    write(message2_text, '(F6.2,a)') CIF_torsion_limit, '   ! '
    call write_debug_proc ('CIF_torsion_limit', trim(message2_text))
    if (include_res_file) then
     call write_debug_proc ('Include_RES_file',      '1     ! ')
    else
     call write_debug_proc ('Include_RES_file',      '0     ! ')
    endif
    if (include_hkl_file) then
     call write_debug_proc ('Include_HKL_file',      '1     ! ')
    else
     call write_debug_proc ('Include_HKL_file',      '0     ! ')
    endif
    if (include_experimenter_name) then
     call write_debug_proc ('Include_experimenter',  '1     ! ')
    else
     call write_debug_proc ('Include_experimenter',      '0     ! ')
    endif
    if (report_header) then
     call write_debug_proc ('Report_header'    ,    '1     ! ')
    else
     call write_debug_proc ('Report_header'    ,    '0     ! ')
    endif
    write(message2_text, '(2a)') trim(report_logo(1)), '   ! '
    call write_debug_proc ('Report_logo_1', trim(message2_text))
    write(message2_text, '(2a)') trim(report_logo(2)), '   ! '
    call write_debug_proc ('Report_logo_2', trim(message2_text))
   endif


   !-------- PATTERN SIMULATION --------------------------------------------------------------------
   rewind(unit=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    IF(read_line(1:19) == '[PATTERN SIMULATION') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit

      READ(read_line(i2+1:), '(a)') INI_string
      INI_string = adjustl(INI_string)
      IF(l_case(read_line(1:11)) == 'x_profile_u') then
       read(INI_string, *) X_PV%U
       cycle
      ELSEIF(l_case(read_line(1:11)) == 'x_profile_v') then
       read(INI_string, *) X_PV%V
       cycle
      ELSEIF(l_case(read_line(1:11)) == 'x_profile_w') then
       read(INI_string, *) X_PV%W
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'x_profile_eta0') then
       read(INI_string, *) X_PV%eta0
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'x_profile_eta1') then
       read(INI_string, *) X_PV%eta1
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'x_pattern_step') then
       read(INI_string, *) X_pattern%step
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'x_pattern_xmin') then
       read(INI_string, *) X_pattern%xmin
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'x_pattern_xmax') then
       read(INI_string, *) X_pattern%xmax
       cycle

      ELSEIF(l_case(read_line(1:20)) == 'x_pattern_background') then
       read(INI_string, *) X_pattern%background
       cycle
      ELSEIF(l_case(read_line(1:15)) == 'x_pattern_scale') then
       read(INI_string, *) X_pattern%scale
       cycle

      ELSEIF(l_case(read_line(1:11)) == 'n_profile_u') then
       read(INI_string, *) N_PV%U
       cycle
      ELSEIF(l_case(read_line(1:11)) == 'n_profile_v') then
       read(INI_string, *) N_PV%V
       cycle
      ELSEIF(l_case(read_line(1:11)) == 'n_profile_w') then
       read(INI_string, *) N_PV%W
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'n_profile_eta0') then
       read(INI_string, *) N_PV%eta0
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'n_profile_eta1') then
       read(INI_string, *) N_PV%eta1
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'n_pattern_step') then
       read(INI_string, *) N_pattern%step
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'n_pattern_xmin') then
       read(INI_string, *) N_pattern%xmin
       cycle
      ELSEIF(l_case(read_line(1:14)) == 'n_pattern_xmax') then
       read(INI_string, *) N_pattern%xmax
       cycle

      ELSEIF(l_case(read_line(1:20)) == 'n_pattern_background') then
       read(INI_string, *) N_pattern%background
       cycle
      ELSEIF(l_case(read_line(1:15)) == 'n_pattern_scale') then
       read(INI_string, *) N_pattern%scale
       cycle
      endif
     end do
    endif
   end do

   if(debug_proc%write) then
    write(message2_text, '(F10.5,a)') X_PV%U, '   ! '
    call write_debug_proc ('X_profile_U', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_PV%V, '   ! '
    call write_debug_proc ('X_profile_V', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_PV%W, '   ! '
    call write_debug_proc ('X_profile_W', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_PV%eta0,'   ! '
    call write_debug_proc ('X_profile_eta0', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_PV%eta1,'   ! '
    call write_debug_proc ('X_profile_eta1', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_pattern%step, '   ! '
    call write_debug_proc ('Pattern step', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_pattern%xmin, '   ! '
    call write_debug_proc ('Pattern Xmin', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_pattern%xmax, '   ! '
    call write_debug_proc ('Pattern Xmax', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_pattern%background, '   ! '
    call write_debug_proc ('X_pattern background', trim(message2_text))
    write(message2_text, '(F10.5,a)') X_pattern%scale, '   ! '
    call write_debug_proc ('X_pattern scale', trim(message2_text))

    write(message2_text, '(F10.5,a)') N_PV%U, '   ! '
    call write_debug_proc ('N_profile_U', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_PV%V, '   ! '
    call write_debug_proc ('N_profile_V', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_PV%W, '   ! '
    call write_debug_proc ('N_profile_W', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_PV%eta0,'   ! '
    call write_debug_proc ('N_profile_eta0', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_PV%eta1,'   ! '
    call write_debug_proc ('N_profile_eta1', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_pattern%step, '   ! '
    call write_debug_proc ('Pattern step', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_pattern%xmin, '   ! '
    call write_debug_proc ('Pattern Xmin', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_pattern%xmax, '   ! '
    call write_debug_proc ('Pattern Xmax', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_pattern%background, '   ! '
    call write_debug_proc ('Neutron pattern background', trim(message2_text))
    write(message2_text, '(F10.5,a)') N_pattern%scale, '   ! '
    call write_debug_proc ('Neutron pattern scale', trim(message2_text))

   endif


   !-------- USER's MATRICES --------------------------------------------------------------------
   rewind(unit=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)


    IF(read_line(1:30) == '[USER TRANSFORMATION MATRICES]') then
     user_mat_nb = 0
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      long = len_trim(read_line)

      !i1 = index(read_line, '!')
      !if(i1 /=0) then
      ! if(i1 == 1) cycle
      ! read_line = read_line(1:i1-1)
      !endif
      !read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit
      READ(read_line(i2+1:), '(a)') INI_string
      INI_string = adjustl(INI_string)
      IF(l_case(read_line(1:5)) == 'mat_1') then
       arg_mat(:,:) = 0.
       read(INI_string, *, iostat=iostat_err) arg_mat(:,1), arg_mat(:,2), arg_mat(:,3)
       if(iostat_err ==0)   transf_mat(:,:,max_mat_nb+1) = arg_mat(:,:)
       user_mat_nb = user_mat_nb  + 1
       i1 = index(read_line, '!')
       if(i1 /= long) then
        user_mat_text(user_mat_nb) = '?'
        read(read_line(i1+1:long), '(a)', iostat=iostat_err) INI_string
        if(iostat_err ==0) then
         user_mat_text(user_mat_nb) = INI_string
         user_mat_text(user_mat_nb) = adjustl(user_mat_text(user_mat_nb))
        endif
       end if
       cycle
      ELSEIF(l_case(read_line(1:5)) == 'mat_2') then
       arg_mat(:,:) = 0.
       read(INI_string, *, iostat=iostat_err) arg_mat(:,1), arg_mat(:,2), arg_mat(:,3)
       if(iostat_err ==0)   transf_mat(:,:,max_mat_nb+2) = arg_mat(:,:)
       user_mat_nb = user_mat_nb  + 1
       i1 = index(read_line, '!')
       if(i1 /= long) then
        user_mat_text(user_mat_nb) = '?'
        read(read_line(i1+1:long), '(a)', iostat=iostat_err) INI_string
        if(iostat_err ==0) then
         user_mat_text(user_mat_nb) = INI_string
         user_mat_text(user_mat_nb) = adjustl(user_mat_text(user_mat_nb))
        end if
       end if
       cycle
      ELSEIF(l_case(read_line(1:5)) == 'mat_3') then
       arg_mat(:,:) = 0.
       read(INI_string, *, iostat=iostat_err) arg_mat(:,1), arg_mat(:,2), arg_mat(:,3)
       if(iostat_err ==0)   transf_mat(:,:,max_mat_nb+3) = arg_mat(:,:)
       user_mat_nb = user_mat_nb  + 1
       i1 = index(read_line, '!')
       if(i1 /= long) then
        user_mat_text(user_mat_nb) = '?'
        read(read_line(i1+1:long), '(a)', iostat=iostat_err) INI_string
        if(iostat_err ==0) then
         user_mat_text(user_mat_nb) = INI_string
         user_mat_text(user_mat_nb) = adjustl(user_mat_text(user_mat_nb))
        endif
       end if
       cycle
      ELSEIF(l_case(read_line(1:5)) == 'mat_4') then
       arg_mat(:,:) = 0.
       read(INI_string, *, iostat=iostat_err) arg_mat(:,1), arg_mat(:,2), arg_mat(:,3)
       if(iostat_err ==0)   transf_mat(:,:,max_mat_nb+4) = arg_mat(:,:)
       user_mat_nb = user_mat_nb  + 1
       i1 = index(read_line, '!')
       if(i1 /= long) then
        user_mat_text(user_mat_nb) = '?'
        read(read_line(i1+1:long), '(a)', iostat=iostat_err) INI_string
        if(iostat_err ==0) then
         user_mat_text(user_mat_nb) = INI_string
         user_mat_text(user_mat_nb) = adjustl(user_mat_text(user_mat_nb))
        endif
       end if
       cycle
      ELSEIF(l_case(read_line(1:5)) == 'mat_5') then
       arg_mat(:,:) = 0.
       read(INI_string, *, iostat=iostat_err) arg_mat(:,1), arg_mat(:,2), arg_mat(:,3)
       if(iostat_err ==0)   transf_mat(:,:,max_mat_nb+5) = arg_mat(:,:)
       user_mat_nb = user_mat_nb  + 1
       i1 = index(read_line, '!')
       if(i1 /= long) then
        user_mat_text(user_mat_nb) = '?'
        read(read_line(i1+1:long), '(a)', iostat=iostat_err) INI_string
        if(iostat_err ==0) then
         user_mat_text(user_mat_nb) = INI_string
         user_mat_text(user_mat_nb) = adjustl(user_mat_text(user_mat_nb))
        endif
       end if
       cycle
      endif
     end do
    endif

   end do


   if(debug_proc%write) then
    do i=1, user_mat_nb
     write(message_text,  '(a,i2)')        'User matrix #', i
     write(message2_text, '(3(2x,3F4.0),2a)')  transf_mat(:,:, max_mat_nb+i), '   ! ' , trim(user_mat_text(i))
    call write_debug_proc (trim(message_text), trim(message2_text))
    end do
   endif

!-------- USER's SHORTCUTS --------------------------------------------------------------------
  rewind(unit=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    if(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)

    nb_shortcuts = 0
    IF(read_line(1:16) == '[USER SHORTCUTS]') then!

     n_s = 0
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = INDEX(read_line, '=')
      IF(i1 ==0) exit
      i2 = index(read_line, '!')
      if(i2 /=0 .and. i2 > i1) then
       i2 = i2 -1
      else
       i2 = len_trim(read_line)
      endif
      if(n_s == user_shortcut_max) exit
      n_s = n_s  + 1
      READ(read_line(1:i1-1), '(a)') shortcut_kw(n_s)
      shortcut_kw(n_s) = trim(shortcut_kw(n_s))
      READ(read_line(i1+1:i2),  '(a)') shortcut_details(n_s)
      shortcut_details(n_s) = adjustl(shortcut_details(n_s))
      shortcut_details(n_s) = trim(shortcut_details(n_s))
     end do
     nb_shortcuts = n_s
    end if
   end do

   if(debug_proc%write) then
    do i=1, nb_shortcuts
     write(message_text,  '(a,i2)')        'User shortcut #', i
     write(message2_text, '(3a)')          shortcut_kw(i)(1:20)," = " , trim(shortcut_details(i))
     call write_debug_proc (trim(message_text), trim(message2_text))
    end do
   endif


!-------- STRUCTURE PROGRAMS --------------------------------------------------------------------
   rewind(unit=INI_unit)
   do
    READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
    IF(iostat_err /=0)           exit
    IF(len_trim(read_line) == 0) cycle
    i1 = index(read_line, '!')
    if(i1 /=0) then
     if(i1 == 1) cycle
     read_line = read_line(1:i1-1)
    endif
    read_line = ADJUSTL(read_line)
    IF(read_line(1:20) == '[STRUCTURE PROGRAMS]' .or. read_line(1:10) == '[PROGRAMS]') then
     do
      READ(INI_unit, '(a)', IOSTAT=iostat_err) read_line
      IF(iostat_err /=0)           exit
      IF(LEN_TRIM(read_line) == 0) exit
      i1 = index(read_line, '!')
      if(i1 /=0) then
       if(i1 == 1) cycle
       read_line = read_line(1:i1-1)
      endif
      read_line = ADJUSTL(read_line)
      i2 = INDEX(read_line, '=')
      IF(i2 ==0) exit


      READ(read_line(i2+1:), '(a)') INI_string
      INI_string = adjustl(INI_string)
      IF(l_case(read_line(1:23)) == 'structure_solution_name') then
       structure_solution%name = INI_string
       long = len_trim(structure_solution%name)
       if(long == 3) then
        if(l_case(structure_solution%name(1:3)) == 'sir') then
         CIF_parameter%atom_sites_solution_1 = "direct"      ! "structure-invariant direct methods"
         cycle
        end if
       elseif(long == 5) then
        if(l_case(structure_solution%name(1:5)) == 'shelxt') then
         CIF_parameter%atom_sites_solution_1 = "dual"        ! "dual-space method"
         cycle
        end if
       elseif(long == 8) then
        if(l_case(structure_solution%name(1:8)) == 'superflip') then
         CIF_parameter%atom_sites_solution_1 = "iterative"   ! "iterative algorithm, e.g. charge flipping"
         cycle
        end if
       end if
      ELSEIF(l_case(read_line(1:28)) == 'structure_solution_reference') then
       structure_solution%reference = INI_string
       cycle
      ELSEIF(l_case(read_line(1:26)) == 'structure_solution_cif_ref') then
       structure_solution%cif_ref = INI_string
       CIF_parameter_DEVICE%computing_structure_solution = trim(structure_solution%cif_ref)
       cycle
      ELSEIF(l_case(read_line(1:25)) == 'structure_refinement_name') then
       structure_refinement%name = INI_string
       !if(len_trim(structure_refinement%name) == 11) then
       ! if(structure_refinement%name(1:11) == 'SHELXL-2014') SHELXL_2014 = .true.
       !endif
       cycle
      ELSEIF(l_case(read_line(1:30)) == 'structure_refinement_reference') then
       structure_refinement%reference = INI_string
       cycle
      ELSEIF(l_case(read_line(1:28)) == 'structure_refinement_cif_ref') then
       structure_refinement%cif_ref = INI_string
       CIF_parameter_DEVICE%computing_structure_refinement = structure_refinement%cif_ref
       cycle
      ELSEIF(l_case(read_line(1:26)) == 'absorption_correction_name') then
       absorption_correction%name = INI_string
       cycle
      ELSEIF(l_case(read_line(1:31)) == 'absorption_correction_reference') then
       absorption_correction%reference = INI_string
       cycle
      ELSEIF(l_case(read_line(1:29)) == 'absorption_correction_cif_ref') then
       absorption_correction%cif_ref = INI_string
       cycle
      endif
     end do
    endif
   end do

   if(debug_proc%write) then
    call write_debug_proc ('Structure_solution_name',         trim(structure_solution%name))
    call write_debug_proc ('Structure_solution_reference',    trim(structure_solution%reference))
    call write_debug_proc ('Structure_solution_cif_ref',      trim(structure_solution%cif_ref))
    call write_debug_proc ('Structure_refinement_name',       trim(structure_refinement%name))
    call write_debug_proc ('Structure_refinement_reference',  trim(structure_refinement%reference))
    call write_debug_proc ('Structure_refinement_cif_ref',    trim(structure_refinement%cif_ref))
    call write_debug_proc ('Absorption_correction_name',      trim(absorption_correction%name))
    call write_debug_proc ('Absorption_correction_reference', trim(absorption_correction%reference))
    call write_debug_proc ('Absorption_correction_cif_ref',   trim(absorption_correction%cif_ref))
    call write_debug_proc ('','')
    call write_debug_proc ('','')
   endif


   CLOSE(UNIT=INI_unit)

  RETURN
 end subroutine read_cryscalc_ini
!-----------------------------------------------------

 subroutine edit_a_file()
  USE cryscalc_module, ONLY : my_editor, file_to_edit
  USE IO_module,       ONLY : write_info

  IF(.not. my_editor%exist) then
   call write_info('')
   call write_info('  No editor defined in the CRYSCALC.ini setting file .')
   call write_info('')
   return
  endif

  call system(TRIM(my_editor%name)//' '//TRIM(file_to_edit))

  return
 end subroutine  edit_a_file

!------------------------------------------------------------------------------------
subroutine output_setting(input_string)
 USE cryscalc_module, ONLY :  INI_unit, cryscalc, my_browser, my_editor, my_pdflatex, WEB, DEVICE, AUTHOR, message_text,        &
                              CONN_dmax_ini, INI_create, get_sample_ID, INI_unit, search_mono_criteria,                         &
                              structure_solution, structure_refinement, absorption_correction, Create_INS, max_ref,             &
                              lock_wave_value, CIF_format80, include_RES_file, include_HKL_file, include_experimenter_name,     &
                              update_parameters , report_header, report_logo, expert_mode, ONLY_monoclinic_angle,               &
                              hkl_statistics, hkl_format, skip_start_menu, pdp_simu, CIF_torsion_limit, debug_proc,             &
                              cartesian_frame, CIFDEP, keep_bond_str_out
 USE pattern_profile_module
 USE CIF_module
 USE USER_module
 USE MATRIX_list_module, only : user_mat_text, user_mat_nb, transf_mat, max_mat_nb
 USE HKL_module,     ONLY :  n_sig, threshold, overlapping_criteria
 USE IO_module,      ONLY :  write_info
 implicit none
  character (len=*), intent(in)   :: input_string
  integer                 :: i, long
  logical                 :: save_settings
  character (len=256)     :: write_line

  save_settings = .false.
  long = len_trim(input_string)
  if(long == 12) then
   if(input_string(1:12) == 'cryscalc.ini') then
    save_settings = .true.
    open(unit=INI_unit, file='cryscalc.ini')
    write(INI_unit, '(a)') '***********************************************'
    write(INI_unit, '(a)') '*                                             *'
    write(INI_unit, '(a)') '*      C R Y S C A L C   s e t t i n g s      *'
    write(INI_unit, '(a)') '*                                             *'
    write(INI_unit, '(a)') '***********************************************'
    write(INI_unit, '(a)') ''
   end if
  end if

 IF(.not. save_settings .and. LEN_TRIM(cryscalc%ini) == 0) then
   call write_info('')
   call write_info('  No CRYSCALC.ini setting file has been defined.')
   call write_info('')
   return
 endif

 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[EXTERNAL APPLICATIONS]')
 call write_setting_line(save_settings, 'browser     = ' // TRIM(my_editor%name))
 call write_setting_line(save_settings, 'editor      = ' // TRIM(my_browser%name))
 call write_setting_line(save_settings, 'pdflatex    = ' // TRIM(my_pdflatex%name))
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[WEB ADDRESS]')
 if(WEB%num_site /= 0) then
  do i = 1, WEB%num_site
   write(write_line, '(a10,2a)') WEB%name(i), ' = ', trim(WEB%ADDRESS(i))
   call write_setting_line(save_settings, trim(write_line))
  end do
 end if
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[DEVICE]')
 call write_setting_line(save_settings, 'diffractometer = '// trim(DEVICE%diffracto))
 call write_setting_line(save_settings, 'lab            = '// trim(DEVICE%lab))
 call write_setting_line(save_settings, 'radiation      = '// trim(DEVICE%radiation))
 call write_setting_line(save_settings, 'wave_A         = '// trim(DEVICE%wave))
 call write_setting_line(save_settings, 'temperature    = '// trim(CIF_cell_measurement%temperature))
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[USER]')
 call write_setting_line(save_settings, 'name           = '//trim(AUTHOR%name))
 call write_setting_line(save_settings, 'first_name     = '//trim(AUTHOR%first_name))
 call write_setting_line(save_settings, 'address        = '//trim(AUTHOR%ADDRESS))
 call write_setting_line(save_settings, 'mail           = '//trim(AUTHOR%email))
 call write_setting_line(save_settings, 'web            = '//trim(AUTHOR%WEB))
 call write_setting_line(save_settings, 'team           = '//trim(AUTHOR%team))
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[ARRAYS DIMENSIONS]')
 write(write_line, '(a,I8)') 'hkl_reflections = ',   Max_ref
 call write_setting_line(save_settings, trim(write_line))
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[PARAMETERS]')
 write(write_line, '(a,F5.2,a)')  'I_sig          = ', n_sig,     '     ! used in the SEARCH_GROUP procedure'
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F5.2,a)')  'threshold      = ', threshold, '     ! used in the SEARCH_GROUP procedure'
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F5.2,2a)') 'd_max_A        = ', CONN_dmax_ini,   '     ! used with the CONN keyword', &
                                      ' (connectivity calculation)'
 call write_setting_line(save_settings, trim(write_line))
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[CREATE INS]')
 if(get_sample_ID) then
  write(write_line, '(a)') 'get_sample_ID     = 1          ! 0: sample_ID =job / 1: get sample ID from folder name'
 else
  write(write_line, '(a)') 'get_sample_ID     = 0          ! 0: sample_ID =job / 1: get sample ID from folder name'
 end if
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F8.2)')  'temperature       = ', Create_INS%temperature
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F4.1)')  'u_threshold       = ', Create_INS%U_threshold
 call write_setting_line(save_settings, trim(write_line))
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[COMMAND LINE ARGUMENT]')
 if(INI_create%ACE) then
  write(write_line, '(a)') 'create_ACE       = 1    ! .ACE file for Carine'
 else
  write(write_line, '(a)') 'create_ACE       = 0    ! .ACE file for Carine'
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(INI_create%CEL) then
  write(write_line, '(a)') 'create_CEL       = 1    ! .CEL file for PowderCELL'
 else
  write(write_line, '(a)') 'create_CEL       = 0    ! .CEL file for PowderCELL'
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(INI_create%CFL) then
  write(write_line, '(a)') 'create_CFL       = 1    ! .CFL file for CRYSCALC'
 else
  write(write_line, '(a)') 'create_CFL       = 0    ! .CFL file for CRYSCALC'
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(INI_create%FST) then
  write(write_line, '(a)') 'create_FST       = 1    ! .FST file for FPStudio'
 else
  write(write_line, '(a)') 'create_FST       = 0    ! .FST file for FPStudio'
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(INI_create%INS) then
  write(write_line, '(a)') 'create_INS       = 1    ! .INS file for SHELXL'
 else
  write(write_line, '(a)') 'create_INS       = 0    ! .INS file for SHELXL'
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(INI_create%CIF_PYMOL) then
  write(write_line, '(a)') 'create_CIF_PYMOL = 1    ! .CIF file for Pymol'
 else
  write(write_line, '(a)') 'create_CIF_PYMOL = 0    ! .CIF file for Pymol'
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(INI_create%PCR) then
  write(write_line, '(a)') 'create_PCR       = 1    ! .PCR file for FULLPROF'
 else
  write(write_line, '(a)') 'create_PCR       = 0    ! .PCR file for FULLPROF'
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(INI_create%PRF) then
  write(write_line, '(a)') 'create_PAT_PRF   = 1    ! .PRF file for FullProf '
 else
  write(write_line, '(a)') 'create_PAT_PRF   = 0    ! .PRF file for FullProf '
 endif
 call write_setting_line(save_settings, trim(write_line))
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[PROGRAMS]')
 call write_setting_line(save_settings,'structure_solution_name         = '// trim(structure_solution%name))
 call write_setting_line(save_settings,'structure_solution_reference    = '// trim(structure_solution%reference))
 call write_setting_line(save_settings,'structure_solution_cif_ref      = '// trim(structure_solution%cif_ref))
 call write_setting_line(save_settings,'structure_refinement_name       = '// trim(structure_refinement%name))
 call write_setting_line(save_settings,'structure_refinement_reference  = '// trim(structure_refinement%reference))
 call write_setting_line(save_settings,'structure_refinement_cif_ref    = '// trim(structure_refinement%cif_ref))
 call write_setting_line(save_settings,'absorption_correction_name      = '// trim(Absorption_correction%name))
 call write_setting_line(save_settings,'absorption_correction_reference = '// trim(Absorption_correction%reference))
 call write_setting_line(save_settings,'absorption_correction_cif_ref   = '// trim(Absorption_correction%cif_ref))
 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[OPTIONS]')
 write(write_line, '(a,F5.2,a)') 'LOCK_wave_value            = ', lock_wave_value,    &
                                     '       ! lock wavelength to anticathode value'
 call write_setting_line(save_settings, trim(write_line))
 if(update_parameters) then
  write(write_line, '(a)') 'update_parameters          = 1           ! update parameters after transformation'
 else
  write(write_line, '(a)') 'update_parameters          = 0           ! update parameters after transformation'
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(hkl_statistics) then
  write(write_line, '(a)') 'hkl_statistics             = 1           ! output statistics on hkl reflections '
 else
  write(write_line, '(a)') 'hkl_statistics             = 0           ! output statistics on hkl reflections '
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(skip_start_menu) then
  write(write_line, '(a)') 'skip_start_menu            = 1           ! Skip CRYSCALC start menu'
 else
  write(write_line, '(a)') 'skip_start_menu            = 0           ! Skip CRYSCALC start menu'
 end if
 call write_setting_line(save_settings, trim(write_line))
 write(write_line,  '(a)')  'hkl_format                 = '//trim(hkl_format) // '   ! format for .hkl file (h,k,lF2,sig)'
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a)') 'Cartesian frame type       = '//trim(cartesian_frame%type) // '             ! A: x//A ; C:x//c'
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F5.2,a)') 'ref_overlap_criteria       = ', overlapping_criteria,    &
                                     '       ! criteria on hkl index for overlapped reflections'
 call write_setting_line(save_settings, trim(write_line))
 write(write_line,'(a,F5.2,a)') 'search_mono_criteria       = ', search_mono_criteria(1),    &
                                     '       ! max. diff. between monoclinic angle and 90. '
 call write_setting_line(save_settings, trim(write_line))
 if(ONLY_monoclinic_angle) then
  write(write_line, '(2a)') 'search_SG_only_mono        = 1           ! output only monoclinic space groups compatible with ', &
                                                                      'unit cell metrics '
 else
  write(write_line, '(2a)') 'search_SG_only_mono        = 0           ! output not only monoclinic space groups compatible with ', &
                                                                      'unit cell metrics '
 end if
 call write_setting_line(save_settings, trim(write_line))
 if(hkl_statistics) then
  write(write_line, '(2a)') 'hkl_statistics             = 1           ! Output statistics on hkl reflections'
 else
  write(write_line, '(2a)') 'hkl_statistics             = 0           ! Output statistics on hkl reflections'
 endif
 call write_setting_line(save_settings, trim(write_line))
 call write_setting_line(save_settings, 'hkl_format                 = '//trim(hkl_format) &
                                                                       // '   ! format for .hkl file (h,k,l,F2,sig)')
 call write_setting_line(save_settings, 'Cartesian frame type       = '//trim(cartesian_frame%type) // '   ! A: x//A ; C:x//c')

 if(expert_mode) then
  call write_setting_line(save_settings,'expert_mode                = 1           ! Expert mode activated')
  if(keep_bond_str_out) then
   call write_setting_line(save_settings,'bond_str.out               = 1           ! keep bond_str.out')
  else
   call write_setting_line(save_settings,'bond_str.out               = 0           ! remove bond_str.out')
  endif
 else
  call write_setting_line(save_settings,'expert_mode                = 0           ! Expert mode OFF')
 end if
 if(debug_proc%level_2) then
  call write_setting_line(save_settings,'debug_mode level 2         = 1           ! Debug mode (level 2) activated')
 elseif(debug_proc%level_3) then
  call write_setting_line(save_settings,'debug_mode level 3         = 1           ! Debug mode (level 3) activated')
 endif
 if(pdp_simu%cu) then
  call write_setting_line(save_settings,'pdp_cu                     = 1           ' &
                                      //'! Ka Cu for powder diffraction pattern simulation')
 else
  call write_setting_line(save_settings,'pdp_cu                     = 0           ' &
                                      //'! Current wavelength for powder diffraction pattern simulation')
 endif
 if(pdp_simu%beam == 'neutrons') then
  call write_setting_line(save_settings,'pdp_beam                   = N           ! ' &
                //'Simulation of powder diffraction pattern : N for neutrons / X for X-rays')
 else
  call write_setting_line(save_settings,'pdp_beam                   = X           ! ' &
                //'Simulation of powder diffraction pattern : N for neutrons / X for X-rays')
 endif
 write(write_line, '(a,F10.5,a)') 'pdp_wave                   = ', pdp_simu%wave, &
                                  '  ! Wavelength used for powder diffraction pattern simulation'
 call write_setting_line(save_settings, trim(write_line))


 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[ARCHIVE_AND_REPORT]')
 if(CIF_format80) then
  call write_setting_line(save_settings, 'CIF_format80               = 1           ! formatted line, ' &
                               //'when creating a CIF file if more than 80 characters')
 else
  call write_setting_line(save_settings, 'CIF_format80               = 0           ! formatted line, ' &
                               //'when creating a CIF file if more than 80 characters')
 endif
 write(write_line, '(a,F6.2,a)') 'CIF_torsion_limit          = ', CIF_torsion_limit, &
                                     '      ! exclude torsion angle if greater'
 call write_setting_line(save_settings, trim(write_line))
 if(CIFDEP) then
  call write_setting_line(save_settings, 'CIF_author                 = 1           ! include author name and address in CIF file')
 else
  call write_setting_line(save_settings, 'CIF_author                 = 0           ! include author name and address in CIF file')
 end if
 if(include_RES_file) then
  call write_setting_line(save_settings, 'include_RES_file           = 1           ' &
                                       //'! include .RES file in the archive_cryscalc.cif file')
 else
  call write_setting_line(save_settings, 'include_RES_file           = 0           ' &
                                       //'! include .RES file in the archive_cryscalc.cif file')
 end if
 if(include_HKL_file) then
  call write_setting_line(save_settings, 'include_HKL_file           = 1           ' &
                                       //'! include .HKL file in the archive_cryscalc.cif file')
 else
  call write_setting_line(save_settings, 'include_HKL_file           = 0           ' &
                                       //'! include .HKL file in the archive_cryscalc.cif file')
 end if
 if(include_experimenter_name) then
  call write_setting_line(save_settings, 'include_experimenter       = 1           ' &
                                       //'! include experimenter name in the archive_cryscalc.cif file')
 else
  call write_setting_line(save_settings, 'include_experimenter       = 0           ' &
                                       //'! include experimenter name in the archive_cryscalc.cif file')
 end if
 if(report_header) then
  call write_setting_line(save_settings, 'report_header              = 1           ! Write header in structural report')
 else
  call write_setting_line(save_settings, 'report_header              = 0           ! Write header in structural report')
 end if
 call write_setting_line(save_settings, 'report_logo_1              = '//trim(report_logo(1))// ' !')
 call write_setting_line(save_settings, 'report_logo_2              = '//trim(report_logo(2))// ' !')

 call write_setting_line(save_settings, '')
 call write_setting_line(save_settings, '[PATTERN SIMULATION (Pseudo-Voigt function)]')
 write(write_line, '(a,F10.5)') 'X_profile_U            = ', X_PV%U
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_V            = ', X_PV%V
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_W            = ', X_PV%W
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_eta0         = ', X_PV%eta0
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_eta1         = ', X_PV%eta1
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_step         = ', X_pattern%step
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_xmin         = ', X_pattern%xmin
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_xmax         = ', X_pattern%xmax
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_background   = ', X_pattern%background
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'X_profile_scale        = ', X_pattern%scale
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_U            = ', N_PV%U
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_V            = ', N_PV%V
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_W            = ', N_PV%W
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_eta0         = ', N_PV%eta0
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_eta1         = ', N_PV%eta1
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_step         = ', N_pattern%step
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_xmin         = ', N_pattern%xmin
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_xmax         = ', N_pattern%xmax
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_background   = ', N_pattern%background
 call write_setting_line(save_settings, trim(write_line))
 write(write_line, '(a,F10.5)') 'N_profile_scale        = ', N_pattern%scale
 call write_setting_line(save_settings, trim(write_line))

 if(user_mat_nb /=0) then
  call write_setting_line(save_settings, '')
  call write_setting_line(save_settings, '[USER TRANSFORMATION MATRICES]')
  do i=1, user_mat_nb
    write(write_line,'(a,i1,18x,a,3(2x,3F3.0),2a)') 'MAT_',i, '=', transf_mat(:,:,max_mat_nb+i), '   !  ', &
                                                        trim(user_mat_text(i))
    call write_setting_line(save_settings, trim(write_line))
  end do
 end if

 if(nb_shortcuts /=0) then
  call write_setting_line(save_settings, '')
  call write_setting_line(save_settings, '[USER SHORTCUTS]       ! Only in expert mode')
  do i=1, nb_shortcuts
    write(write_line,'(3a)')          shortcut_kw(i)(1:20),"   = " , trim(shortcut_details(i))
    call write_setting_line(save_settings, trim(write_line))
  end do
 end if

 if(save_settings) then
  call write_info('')
  call write_info(' > CRYSCALC settings has been saved in cryscalc.ini file in the current folder.')
  call write_info('')
 end if


end subroutine output_setting

!--------------------------------------------------------------------------------
subroutine get_X_radiation(input_string)
 use wavelength_module
 use cryscalc_module, only : wavelength, keyword_WAVE, neutrons, X_rays

 implicit none
  character (len=*), intent(in)   :: input_string
  integer                         :: long_input_string
  logical                         :: target_Ag, target_Co, target_Cr, target_Cu, target_Fe, target_Mo, target_Ni

  target_Ag    = .false.
  target_Co    = .false.
  target_Cr    = .false.
  target_Cu    = .false.
  target_Fe    = .false.
  target_Mo    = .false.
  target_Ni    = .false.
  anti_cathode = .false.

  long_input_string = len_trim(input_string)

  if(long_input_string == 3) then
   if(input_string(1:3) == 'XAG') target_Ag = .true.
   if(input_string(1:3) == 'XCO') target_Co = .true.
   if(input_string(1:3) == 'XCR') target_Cr = .true.
   if(input_string(1:3) == 'XCU') target_Cu = .true.
   if(input_string(1:3) == 'XFE') target_Fe = .true.
   if(input_string(1:3) == 'XMO') target_Mo = .true.
   if(input_string(1:3) == 'XNi') target_Ni = .true.
  elseif(long_input_string == 4) then
   if(input_string(1:4) == 'X_AG') target_Ag = .true.
   if(input_string(1:4) == 'X_CO') target_Co = .true.
   if(input_string(1:4) == 'X_CR') target_Cr = .true.
   if(input_string(1:4) == 'X_CU') target_Cu = .true.
   if(input_string(1:4) == 'X_FE') target_Fe = .true.
   if(input_string(1:4) == 'X_MO') target_Mo = .true.
   if(input_string(1:4) == 'X_Ni') target_Ni = .true.
  endif


     IF(target_Ag) then
      wavelength = X_target(1)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(1)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Mo) then
      wavelength = X_target(2)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(2)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Cu) then
      wavelength = X_target(3)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(3)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Ni) then
      wavelength = X_target(4)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(4)%logic   = .true.
      anti_cathode = .true.
     ELSEIF(target_Co) then
      wavelength = X_target(5)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(5)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Fe) then
      wavelength = X_target(6)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(6)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Cr) then
      wavelength = X_target(7)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(7)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     else
      X_target(1:tabulated_target_nb)%logic    = .false.
      keyword_WAVE = .false.
      anti_cathode = .false.
     endif

 return
 end subroutine get_X_radiation


!-------------------------------------------------------------------
subroutine write_debug_proc(field, value)
 use cryscalc_module, only : debug_proc
 implicit none
  character (len=*), intent(in)           :: field
  !character (len=*), intent(in), optional :: value
  character (len=*), intent(in)           :: value
  integer                                 :: long, i
  character (len=24)                      :: fmt_


  long = len_trim(field)

!  if(.not. present(value)) then
!   write(debug_proc%unit, '(a)') trim(field)
!
!  else
   !if(value(1:1) == '?' .or. len_trim(value)==0) then
   ! write(debug_proc%unit, '(a)') ''
   !else
   ! write(fmt_, '(a,i2,a)') '(1x,a,', 32-long, 'x, 2a)'
   ! write(debug_proc%unit, trim(fmt_))  field,  ': ', adjustl(value)
   !endif
   if(len_trim(field) == 0) then
    write(debug_proc%unit, '(a)') ''
   elseif(field(1:1) == '!') then
    write(debug_proc%unit, '(1x,a)') trim(field)
   elseif(field(1:1) == '-') then
    write(debug_proc%unit, '(1x,60a1)')  ('-',i=1,60)
   else
    write(fmt_, '(a,i2,a)') '(1x,a,', 32-long, 'x, 2a)'
    write(debug_proc%unit, trim(fmt_))  field,  ': ', adjustl(value)
   endif
!  endif

   return
end subroutine write_debug_proc

!-------------------------------------------------------------------
subroutine write_debug_proc_level(level, field)
 use cryscalc_module, only : debug_proc
 use macros_module,   only : u_case, l_case
 implicit none
  integer,           intent(in)           :: level
  character (len=*), intent(in)           :: field


  select case (level)

   case(2)
    if(len_trim(field)/=0)  write(debug_proc%unit, '(1x,2a)') '. debug_level 2 :  routine = ',  trim(u_case(field))

   case(3)
    if(len_trim(field)/=0)  write(debug_proc%unit, '(1x,2a)') '        level 3 :  ...',  trim(l_case(field))

   end select

   return
end subroutine write_debug_proc_level

!-------------------------------------------------------------------------------
subroutine write_debug_level_1
 use cryscalc_module
 use CIF_module
 use pattern_profile_module
 use hkl_module
 use MATRIX_list_module
 use macros_module,  only : u_case

 implicit none
  character (len=256)             :: string, string_1
  integer                         :: i

    call write_debug_proc ('','')
    call write_debug_proc ('CRYSCALC_ini        :', trim(cryscalc%ini))
    call write_debug_proc ('CRYSCALC_css        :', trim(cryscalc%css))
    call write_debug_proc ('CRYSCALC_report_css :', trim(cryscalc%report_css))

    call write_debug_proc ('BROWSER',            trim(my_browser%name))
    call write_debug_proc ('EDITOR',             trim(my_editor%name))
    call write_debug_proc ('WORD',               trim(my_word%name))
    call write_debug_proc ('PDFLATEX',           trim(my_pdflatex%name))
    if(WEB%num_site /=0) then
     do i=1, WEB%num_site
      call write_debug_proc ('WEB_name',           trim(WEB%NAME(i)))
      call write_debug_proc ('WEB_address',        trim(WEB%ADDRESS(i)))
     end do
    else
     call write_debug_proc ('WEB_name',           '?')
     call write_debug_proc ('WEB_address',        '?')
    end if
    call write_debug_proc ('USER_name',          trim(AUTHOR%name))
    call write_debug_proc ('USER_first_name',    trim(AUTHOR%first_name))
    call write_debug_proc ('USER_address',       trim(AUTHOR%address))
    call write_debug_proc ('USER_web',           trim(AUTHOR%web))
    call write_debug_proc ('USER_team',          trim(AUTHOR%team))

    call write_debug_proc ('DEVICE_diffractometer', trim(DEVICE%diffracto))
    call write_debug_proc ('DEVICE_lab',            trim(DEVICE%lab))
    call write_debug_proc ('DEVICE_radiation',      trim(DEVICE%radiation))
    call write_debug_proc ('DEVICE_wave',           trim(DEVICE%wave))

    !if(u_case(DEVICE%diffracto(1:4)) == 'APEX') then
    if(DEVICE%APEX .or. DEVICE%D8V) then
     call write_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call write_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call write_debug_proc('MEASUREMENT_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call write_debug_proc('RADIATION_wavelength',      trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
     call write_debug_proc('RADIATION_type',            trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call write_debug_proc('RADIATION_source',          trim(CIF_parameter_DEVICE%diffrn_radiation_source))
     call write_debug_proc('RADIATION_monochromator',   trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
     call write_debug_proc('RADIATION_probe',           trim(CIF_parameter_DEVICE%diffrn_radiation_probe))

    !elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
    elseif(DEVICE%KCCD) then
     call write_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call write_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call write_debug_proc('MEASUREMENT_device',        trim(CIF_parameter_DEVICE%diffrn_measurement_device))
     call write_debug_proc('DIFFRACTION_source',        trim(CIF_parameter_DEVICE%diffrn_source))
     call write_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call write_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean))
     call write_debug_proc('RADIATION_wavelength',      trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
     call write_debug_proc('RADIATION_type',            trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call write_debug_proc('RADIATION_source',          trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call write_debug_proc('RADIATION_monochromator',   trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
     call write_debug_proc('RADIATION_probe',           trim(CIF_parameter_DEVICE%diffrn_radiation_probe))

    !elseif(u_case(DEVICE%diffracto(1:7)) == 'SAPHIRE' .or. u_case(DEVICE%diffracto(1:9)) == 'EXCALIBUR') then
    elseif(DEVICE%Xcalibur) then
     call write_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call write_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call write_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call write_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter_DEVICE%diffrn_detector_area_resol_mean))

    !elseif(u_case(DEVICE%diffracto(1:3)) == 'X2S') then
    elseif(DEVICE%X2S) then
     call write_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
     call write_debug_proc('MEASUREMENT_method',        trim(CIF_parameter_DEVICE%diffrn_measurement_method))
     call write_debug_proc('MEASUREMENT_detector',      trim(CIF_parameter_DEVICE%diffrn_detector))
     call write_debug_proc('RADIATION_wavelength',      trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
     call write_debug_proc('RADIATION_type',            trim(CIF_parameter_DEVICE%diffrn_radiation_type))
     call write_debug_proc('RADIATION_source',          trim(CIF_parameter_DEVICE%diffrn_radiation_source))
     call write_debug_proc('RADIATION_monochromator',   trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
     call write_debug_proc('RADIATION_probe',           trim(CIF_parameter_DEVICE%diffrn_radiation_probe))
    endif

    call write_debug_proc('COMPUTING_data_collection',  trim(CIF_parameter_DEVICE%computing_data_collection))
    call write_debug_proc('COMPUTING_cell_refinement',  trim(CIF_parameter_DEVICE%computing_cell_refinement))
    call write_debug_proc('COMPUTING_data_reduction',   trim(CIF_parameter_DEVICE%computing_data_reduction))

    write(string, '(F5.2,a)') n_sig,          '    ! sigma_I (used in the SEARCH_GROUP routine)'
    call write_debug_proc ('SEARCH_I_sig',      trim(string))
    write(string, '(F6.3,a)') threshold,      '   ! I > threshold * I_max (used in the SEARCH_GROUP routine)'
    call write_debug_proc ('SEARCH_threshold',  trim(string))
    write(string, '(F6.2,a)') CONN_dmax_ini,  '    ! (in A) for connectivity calculations '
    call write_debug_proc ('CONNECTIVITY_dmax', trim(string))

    write(string, '(I10,a)') Max_ref,      '    ! Max. number of hkl reflections in a file'
    call write_debug_proc ('Max_ref',  trim(string))
    write(string, '(F8.2,a)') Create_INS%temperature,      '  ! experimental temperature'
    call write_debug_proc ('temperature',      trim(string))
    write(string, '(F8.2,a)') Create_INS%U_threshold,      '    ! atoms with Uiso > u_threshold are excluded'
    call write_debug_proc ('U_threshold',  trim(string))

    if (keyword_create_ACE) then
     call write_debug_proc ('create_ACE',      '1     ! .ACE file for Carine')
    else
     call write_debug_proc ('create_ACE',      '0     ! .ACE file for Carine')
    endif
    if (keyword_create_CEL) then
     call write_debug_proc ('create_CEL',      '1     ! .CEL file for PowderCELL')
    else
     call write_debug_proc ('create_CEL',      '0     ! .CEL file for PowderCELL')
    endif
    if (keyword_create_CFL) then
     call write_debug_proc ('create_CFL',      '1     ! .CFL file for CRYSCALC')
    else
     call write_debug_proc ('create_CFL',      '0     ! .CFL file for CRYSCALC')
    endif
    if (keyword_create_FST) then
     call write_debug_proc ('create_FST',      '1     ! .FST file for FP Studio')
    else
     call write_debug_proc ('create_FST',      '0     ! .FST file for FP Studio')
    endif
    if (create_CIF_PYMOL) then
     call write_debug_proc ('create_CIF_PYMOL',    '1     ! .CIF file for PYMOL')
    else
     call write_debug_proc ('create_CIF_PYMOL',    '0     ! .CIF file for PYMOL')
    endif
    if (keyword_create_INS) then
     call write_debug_proc ('create_INS',      '1     ! .INS file for SHELXL')
    else
     call write_debug_proc ('create_INS',      '0     ! .INS file for SHELXL')
    endif
    if (keyword_create_PRF) then
     call write_debug_proc ('create_PAT_PRF',  '1     ! .PRF file for FullProf')
    else
     call write_debug_proc ('create_PAT_PRF',  '0     ! .PRF file for FullProf')
    endif
    if (keyword_create_PCR) then
     call write_debug_proc ('create_PCR',  '1     ! .PCR file for FullProf')
    else
     call write_debug_proc ('create_PCR',  '0     ! .PCR file for FullProf')
    endif

    write(string, '(F6.3,a)') LOCK_wave_value, '   ! '
    call write_debug_proc ('LOCK_wave_value', trim(string))
    if (CIF_format80) then
     call write_debug_proc ('CIF_format80',         '1     ! ')
    else
     call write_debug_proc ('CIF_format80',         '0     ! ')
    endif
    write(string, '(F6.2,a)') CIF_torsion_limit, '   ! '
    call write_debug_proc ('CIF_torsion_limit', trim(string))
    if (include_res_file) then
     call write_debug_proc ('Include_RES_file',      '1     ! ')
    else
     call write_debug_proc ('Include_RES_file',      '0     ! ')
    endif
    if (include_HKL_file) then
     call write_debug_proc ('Include_HKL_file',      '1     ! ')
    else
     call write_debug_proc ('Include_HKL_file',      '0     ! ')
    endif
    if (update_parameters) then
     call write_debug_proc ('Update_parameters',     '1     ! ')
    else
     call write_debug_proc ('Update_parameters',     '0     ! ')
    endif
    if (report_header) then
     call write_debug_proc ('Report_header'    ,    '1     ! ')
    else
     call write_debug_proc ('Report_header'    ,    '0     ! ')
    endif
    call write_debug_proc  ('Report_logo_1'    , trim(report_logo(1))// ' !')
    call write_debug_proc  ('Report_logo_2'    , trim(report_logo(2))// ' !')

    if(expert_mode) then
     call write_debug_proc('Expert_mode'       ,    '1     ! ')
    else
     call write_debug_proc('Expert_mode'       ,    '0     ! ')
    endif
    if(skip_start_menu) then
     call write_debug_proc('Skip_start_menu'   ,     '1     !')
    else
     call write_debug_proc('Skip_start_menu'   ,     '0     !')
    end if
    if(hkl_statistics) then
     call write_debug_proc('hkl_statistics'     ,    '1     ! ')
    else
     call write_debug_proc('hkl_statistics'     ,    '0     ! ')
    endif
    call write_debug_proc('hkl_format'          ,     trim(hkl_format)//'  !')
    call write_debug_proc('cartesian_frame_type',     trim(cartesian_frame%type)//'  !')
    if(expert_mode) then
    if(pdp_simu%cu) then
     call write_debug_proc('Ka1 Cu for pdp simulation'     ,    '1     ! ')
    else
     call write_debug_proc('Ka1 Cu for pdp simulation'     ,    '0     ! ')
    endif
   end if


    write(string, '(F8.5,a)') X_PV%U, '   ! '
    call write_debug_proc ('X_profile_U', trim(string))
    write(string, '(F8.5,a)') X_PV%V, '   ! '
    call write_debug_proc ('X_profile_V', trim(string))
    write(string, '(F8.5,a)') X_PV%W, '   ! '
    call write_debug_proc ('X_profile_W', trim(string))
    write(string, '(F8.5,a)') X_PV%eta0,'   ! '
    call write_debug_proc ('X_profile_eta0', trim(string))
    write(string, '(F8.5,a)') X_PV%eta1,'   ! '
    call write_debug_proc ('X_profile_eta1', trim(string))
    write(string, '(F8.5,a)') X_pattern%step, '   ! '
    call write_debug_proc ('Pattern step', trim(string))
    write(string, '(F10.5,a)') X_pattern%background, '   ! '
    call write_debug_proc ('X_pattern background', trim(string))
    write(string, '(F10.5,a)') X_pattern%scale, '   ! '
    call write_debug_proc ('X_pattern scale', trim(string))

    write(string, '(F8.5,a)') N_PV%U, '   ! '
    call write_debug_proc ('N_profile_U', trim(string))
    write(string, '(F8.5,a)') N_PV%V, '   ! '
    call write_debug_proc ('N_profile_V', trim(string))
    write(string, '(F8.5,a)') N_PV%W, '   ! '
    call write_debug_proc ('N_profile_W', trim(string))
    write(string, '(F8.5,a)') N_PV%eta0,'   ! '
    call write_debug_proc ('N_profile_eta0', trim(string))
    write(string, '(F8.5,a)') N_PV%eta1,'   ! '
    call write_debug_proc ('N_profile_eta1', trim(string))
    write(string, '(F8.5,a)') N_pattern%step, '   ! '
    call write_debug_proc ('Pattern step', trim(string))
    write(string, '(F10.5,a)') N_pattern%background, '   ! '
    call write_debug_proc ('Neutron pattern background', trim(string))
    write(string, '(F10.5,a)') N_pattern%scale, '   ! '
    call write_debug_proc ('Neutron pattern scale', trim(string))

    call write_debug_proc ('Structure_solution_name',         trim(structure_solution%name))
    call write_debug_proc ('Structure_solution_reference',    trim(structure_solution%reference))
    call write_debug_proc ('Structure_solution_cif_ref',      trim(structure_solution%cif_ref))
    call write_debug_proc ('Structure_refinement_name',       trim(structure_refinement%name))
    call write_debug_proc ('Structure_refinement_reference',  trim(structure_refinement%reference))
    call write_debug_proc ('Structure_refinement_cif_ref',    trim(structure_refinement%cif_ref))
    call write_debug_proc ('Absorption_correction_name',      trim(absorption_correction%name))
    call write_debug_proc ('Absorption_correction_reference', trim(absorption_correction%reference))
    call write_debug_proc ('Absorption_correction_cif_ref',   trim(absorption_correction%cif_ref))

    do i=1, user_mat_nb
     write(string ,  '(a,i2)')        'User matrix #', i
     write(string_1, '(3(2x,3F4.0),2a)')  transf_mat(:,:, max_mat_nb+i), '   ! ' , trim(user_mat_text(i))
     call write_DEBUG_proc (trim(string), trim(string_1))
    end do

 return
end subroutine write_debug_level_1

!-------------------------------------------------------------------------------
subroutine CIF_default_values
 use cryscalc_module, only : debug_proc, my_browser, my_editor, my_word, my_pdflatex, AUTHOR, DEVICE
 use CIF_module
 implicit none

 CIF_parameter_DEVICE = CIF_parameter_DEVICE_APEX
 !CIF_parameter_DEVICE%diffrn_measurement_device_type       = CIF_parameter_DEVICE_APEX%diffrn_measurement_device_type
 !CIF_parameter_DEVICE%diffrn_measurement_method            = CIF_parameter_DEVICE_APEX%diffrn_measurement_method
 !CIF_parameter_DEVICE%diffrn_detector                      = CIF_parameter_DEVICE_APEX%diffrn_detector
 !CIF_parameter_DEVICE%diffrn_radiation_wavelength          = CIF_parameter_DEVICE_APEX%diffrn_radiation_wavelength
 !CIF_parameter_DEVICE%diffrn_radiation_type                = CIF_parameter_DEVICE_APEX%diffrn_radiation_type
 !CIF_parameter_DEVICE%diffrn_radiation_source              = CIF_parameter_DEVICE_APEX%diffrn_radiation_type
 !CIF_parameter_DEVICE%diffrn_radiation_monochromator       = CIF_parameter_DEVICE_APEX%diffrn_radiation_monochromator
 !CIF_parameter_DEVICE%diffrn_radiation_probe               = CIF_parameter_DEVICE_APEX%diffrn_radiation_probe
 !CIF_parameter_DEVICE%computing_data_collection            = CIF_parameter_DEVICE_APEX%computing_data_collection
 !CIF_parameter_DEVICE%computing_cell_refinement            = CIF_parameter_DEVICE_APEX%computing_cell_refinement
 !CIF_parameter_DEVICE%computing_data_reduction             = CIF_parameter_DEVICE_APEX%computing_data_reduction
 !CIF_parameter_DEVICE%computing_structure_solution         = CIF_parameter_DEVICE_APEX%computing_structure_solution
 !CIF_parameter_DEVICE%computing_structure_refinement       = CIF_parameter_APEX%computing_structure_refinement
 !CIF_parameter_DEVICE%computing_molecular_graphics         = CIF_parameter_APEX%computing_molecular_graphics
 !!CIF_parameter_DEVICE%computing_publication_material       = CIF_parameter__DEVICEAPEX%computing_publication_material
 !CIF_parameter_DEVICE%computing_publication_material_1     = CIF_parameter_DEVICE_APEX%computing_publication_material_1
 !CIF_parameter_DEVICE%computing_publication_material_2     = CIF_parameter_DEVICE_APEX%computing_publication_material_2

 if(debug_proc%write) then
   call write_debug_proc ('-', '')
   call write_debug_proc ('BROWSER',                        trim(my_browser%name))
   call write_debug_proc ('EDITOR',                         trim(my_editor%name))
   call write_debug_proc ('WORD',                           trim(my_word%name))
   call write_debug_proc ('PDFLATEX',                       trim(my_pdflatex%name))

   call write_debug_proc ('USER_name',                    trim(AUTHOR%name))
   call write_debug_proc ('USER_first_name',              trim(AUTHOR%first_name))
   call write_debug_proc ('USER_first_address',           trim(AUTHOR%address))
   call write_debug_proc ('USER_first_web',               trim(AUTHOR%web))

   call write_debug_proc ('DEVICE_diffractometer',          trim(DEVICE%diffracto))
   call write_debug_proc ('DEVICE_lab',                     trim(DEVICE%lab))
   call write_debug_proc ('DEVICE_radiation',               trim(DEVICE%radiation))
   call write_debug_proc ('DEVICE_wave',                    trim(DEVICE%wave))
   call WRITE_debug_proc ('MEASUREMENT_device_type',        trim(CIF_parameter_DEVICE%diffrn_measurement_device_type))
   call WRITE_debug_proc ('MEASUREMENT_method',             trim(CIF_parameter_DEVICE%diffrn_measurement_method))
   call WRITE_debug_proc ('MEASUREMENT_detector',           trim(CIF_parameter_DEVICE%diffrn_detector))
   call WRITE_debug_proc ('RADIATION_wavelength',           trim(CIF_parameter_DEVICE%diffrn_radiation_wavelength))
   call WRITE_debug_proc ('RADIATION_type',                 trim(CIF_parameter_DEVICE%diffrn_radiation_type))
   call WRITE_debug_proc ('RADIATION_source',               trim(CIF_parameter_DEVICE%diffrn_radiation_source))
   call WRITE_debug_proc ('RADIATION_monochromator',        trim(CIF_parameter_DEVICE%diffrn_radiation_monochromator))
   call WRITE_debug_proc ('RADIATION_probe',                trim(CIF_parameter_DEVICE%diffrn_radiation_probe))
   call WRITE_debug_proc ('COMPUTING_data_collection',      trim(CIF_parameter_DEVICE%computing_data_collection))
   call WRITE_debug_proc ('COMPUTING_cell_refinement',      trim(CIF_parameter_DEVICE%computing_cell_refinement))
   call WRITE_debug_proc ('COMPUTING_data_reduction',       trim(CIF_parameter_DEVICE%computing_data_reduction))
   call WRITE_debug_proc ('COMPUTING_structure_solution',   trim(CIF_parameter_DEVICE%computing_structure_solution))
   call WRITE_debug_proc ('COMPUTING_structure_refinement', trim(CIF_parameter_DEVICE%computing_structure_refinement))
   call WRITE_debug_proc ('COMPUTING_molecular_graphics',   trim(CIF_parameter_DEVICE%computing_molecular_graphics))
   !call WRITE_debug_proc ('COMPUTING_publication_material', trim(CIF_parameter_DEVICE%computing_publication_material))
   call WRITE_debug_proc ('COMPUTING_publication_material', trim(CIF_parameter_DEVICE%computing_publication_material_1))
   call WRITE_debug_proc ('COMPUTING_publication_material', trim(CIF_parameter_DEVICE%computing_publication_material_2))



 endif

 return
end subroutine CIF_default_values

!!---------------------------------------------------------------------------------------
subroutine write_default_values
 use cryscalc_module
 use pattern_profile_module
 use hkl_module
 implicit none
   character (len=256)          :: string


   write(string, '(F5.2,a)') n_sig,          '    ! sigma_I (used in the SEARCH_GROUP routine)'
   call write_DEBUG_proc ('SEARCH_I_sig',      trim(string))
   write(string, '(F6.3,a)') threshold,      '   ! I > threshold * I_max (used in the SEARCH_GROUP routine)'
   call write_DEBUG_proc ('SEARCH_threshold',  trim(string))
   write(string, '(F6.2,a)') CONN_dmax_ini,  '    ! (in A) for connectivity calculations '
   call write_DEBUG_proc ('CONNECTIVITY_dmax', trim(string))
   write(string, '(I10,a)') Max_ref,      '    ! Max. number of hkl reflections in a file'
   call write_DEBUG_proc ('Max_ref',  trim(string))

   write(string, '(F8.2,a)') Create_INS%temperature,      '  ! experimental temperature'
   call write_DEBUG_proc ('temperature',      trim(string))
   write(string, '(F8.2,a)') Create_INS%U_threshold,      '  ! atoms with Uiso > u_threshold are excluded'
   call write_DEBUG_proc ('U_threshold',  trim(string))


   call write_DEBUG_proc ('create_ACE',          '0     ! .ACE file for Carine')
   call write_DEBUG_proc ('create_CEL',          '0     ! .CEL file for PowderCELL')
   call write_DEBUG_proc ('create_CFL',          '0     ! .CFL file for CRYSCALC')
   call write_DEBUG_proc ('create_FST',          '0     ! .FST file for FP Studio')
   call write_DEBUG_proc ('create_CIF_PYMOL',    '0     ! .CIF file for PYMOL')
   call write_DEBUG_proc ('create_INS',          '0     ! .INS file for SHELXL')
   call write_DEBUG_proc ('create_PAT_PRF',      '0     ! .PRF file for FullProf')



   write(string, '(F6.3,a)') LOCK_wave_value, '   ! '
   call write_DEBUG_proc ('LOCK_wave_value', trim(string))

   call write_DEBUG_proc ('CIF_format80'      ,    '1     ! ')
   write(string, '(F6.2,a)') CIF_torsion_limit,    '   ! '
   call write_DEBUG_proc ('CIF_torsion_limit' , trim(string))
   call write_DEBUG_proc ('Include_RES_file'  ,    '0     ! ')
   call write_DEBUG_proc ('Include_HKL_file'  ,    '0     ! ')
   call write_DEBUG_proc ('Include_experimenter' , '1     ! ')
   call write_DEBUG_proc ('Update_parameters' ,    '1     ! ')
   call write_DEBUG_proc ('Report_header'     ,    '1     ! ')
   call write_DEBUG_proc ('Report_logo(1)'    ,    trim(report_logo(1))//' !')
   call write_DEBUG_proc ('Report_logo(2)'    ,    trim(report_logo(2))//' !')
   call write_DEBUG_proc ('Expert_mode'       ,    '0     ! ')
   call write_DEBUG_proc ('Skip_start_menu'   ,    '0     !')
   call write_DEBUG_proc ('hkl_statistics'    ,    '1     ! ')

   call write_DEBUG_proc ('hkl_format'           ,     trim(hkl_format)//'  !')
   call write_DEBUG_proc ('cartesian_frame_type' ,     trim(cartesian_frame%string)//'  !')
   call write_DEBUG_proc ('Ka1 Cu for pdp simulation'     ,    '0     ! ')

   write(string, '(F6.3,a)') overlapping_criteria, '   ! '
   call write_DEBUG_proc ('ref_overlap_criteria', trim(string))
   write(string, '(F6.3,a)') search_mono_criteria(1), '   ! '
   call write_DEBUG_proc ('search_mono_criteria', trim(string))
   call write_DEBUG_proc ('search_SG_only_mono'    ,    '0     ! ')


    write(string, '(F8.5,a)') X_PV%U, '   ! '
    call write_DEBUG_proc ('X_profile_U', trim(string))
    write(string, '(F8.5,a)') X_PV%V, '   ! '
    call write_DEBUG_proc ('X_profile_V', trim(string))
    write(string, '(F8.5,a)') X_PV%W, '   ! '
    call write_DEBUG_proc ('X_profile_W', trim(string))
    write(string, '(F8.5,a)') X_PV%eta0,'   ! '
    call write_DEBUG_proc ('X_profile_eta0', trim(string))
    write(string, '(F8.5,a)') X_PV%eta1,'   ! '
    call write_DEBUG_proc ('X_profile_eta1', trim(string))
    write(string, '(F8.5,a)') X_pattern%step, '   ! '
    call write_DEBUG_proc ('X_Pattern step', trim(string))
    write(string, '(F10.5,a)') X_pattern%background, '   ! '
    call write_DEBUG_proc ('X_pattern background', trim(string))
    write(string, '(F10.5,a)') X_pattern%scale, '   ! '
    call write_DEBUG_proc ('X_pattern scale', trim(string))


    write(string, '(F8.5,a)') N_PV%U, '   ! '
    call write_DEBUG_proc ('N_profile_U', trim(string))
    write(string, '(F8.5,a)') N_PV%V, '   ! '
    call write_DEBUG_proc ('N_profile_V', trim(string))
    write(string, '(F8.5,a)') N_PV%W, '   ! '
    call write_DEBUG_proc ('N_profile_W', trim(string))
    write(string, '(F8.5,a)') N_PV%eta0,'   ! '
    call write_DEBUG_proc ('N_profile_eta0', trim(string))
    write(string, '(F8.5,a)') N_PV%eta1,'   ! '
    call write_DEBUG_proc ('N_profile_eta1', trim(string))
    write(string, '(F8.5,a)') N_pattern%step, '   ! '
    call write_DEBUG_proc ('X_Pattern step', trim(string))
    write(string, '(F10.5,a)') N_pattern%background, '   ! '
    call write_DEBUG_proc ('N_pattern background', trim(string))
    write(string, '(F10.5,a)') N_pattern%scale, '   ! '
    call write_DEBUG_proc ('N_pattern scale', trim(string))

 return

end subroutine write_default_values

!---------------------------------------------------------------------------------------------------------
subroutine open_debug_file
 use cryscalc_module, only : debug_proc, cryscalc
 use io_module
 implicit none
  integer             :: i_error


  open(unit = debug_proc%unit, file = 'cryscalc_debug.txt', status = "replace", iostat=i_error)
    if(i_error /=0) then
     call write_info('')
     call write_info('  !! Problem opening cryscalc_debug.txt file !!')
     call write_info('')
     debug_proc%write = .false.
     close(unit = debug_proc%unit)
     return
    end if

    write(debug_proc%unit, '(a)')  ' '
    write(debug_proc%unit, '(a)')  ' ***************************************************************************'
    write(debug_proc%unit, '(3a)') '      CRYSCALC_debug.txt file (CRYSCALC version : ',trim(cryscalc%version),')'
    if(debug_proc%level_3) then
    write(debug_proc%unit, '(3a)') '           . debug level : 3'
    elseif(debug_proc%level_2) then
    write(debug_proc%unit, '(3a)') '           . debug level : 2'
    elseif(debug_proc%level_1) then
    write(debug_proc%unit, '(3a)') '           . debug level : 1'
    end if
    write(debug_proc%unit, '(a)')  ' ***************************************************************************'
    write(debug_proc%unit, '(a)')  ' '


    return
 end subroutine open_debug_file

!!-------------------------------------------------------------------------------------------------------------------
subroutine write_setting_line(save_settings, input_string)
 use cryscalc_module,   only : ini_unit
 use io_module
 implicit none
  logical, intent(in)              :: save_settings ! T: in cryscalc.ini / F: on screen
  character (len=*), intent(in)    :: input_string
  character (len=256)              :: new_line
  integer                          :: long

  long = len_trim(input_string)

  if(save_settings) then
   write(ini_unit, '(a)') trim(input_string)
  else
   new_line = trim(input_string)
   if(long /=0) then
     if(input_string == '[USER SHORTCUTS]') then
     new_line = trim(new_line)//"   ! Only in expert mode"
    end if

    if(input_string(1:1) /= '[' .and. input_string(long:long) /= ']') then
     new_line = " >  "//trim(input_string)
    end if
   end if
   call write_info(trim(new_line))
  end if

  return
 end subroutine write_setting_line

!!-------------------------------------------------------------------------------------------------------------------

subroutine SPG_init
 use cryscalc_module,                ONLY : SPG
 use CFML_crystallographic_symmetry, ONLY : set_spacegroup, space_group_type

  SPG%NumSPG      = 0
  SPG%SPG_symb    = ""
  SPG%CrystalSys  = ""
  SPG%Laue        = ""
  SPG%Bravais     = ""
  SPG%Multip      = 0
  SPG%centre      = ""

 return
end subroutine SPG_init