!     Last change:  TR   11 Oct 2007    6:11 pm
subroutine cryscalc_init()
 USE cryscalc_module
 USE pattern_profile_module
 USE wavelength_module
 USE HKL_module
 USE CIF_module
 USE MATRIX_list_module, only : user_mat_nb


  implicit none
   INTEGER   :: numor, long

 ! initialisation ---------------------------------------
  winplotr_exe        = ''
  
  CRYSCALC%version      = 'Nov. 2014'          
  CRYSCALC%author       = 'Thierry Roisnel (CDIFX / ISCR UMR6226 - Rennes)'    
  CRYSCALC%mail         = "thierry.roisnel@univ-rennes1.fr"
  CRYSCALC%url          = "www.cdifx.univ-rennes1.fr/cryscalc"  
  CRYSCALC%ini          = ''
  CRYSCALC%css          = ''
  CRYSCALC%report_css   = ''
  CRYSCALC%path_name    = ''
     
  my_editor%name    = ''
  my_browser%name   = ''
  my_word%name      = ''
  my_pdflatex%name  = ''

  WEB%NAME          = ''
  WEB%ADDRESS       = ''


! initialisation des champs lus dans le fichier de config. CRYSCALC.INI
  my_editor%exist   = .false.
  my_browser%exist  = .false.
  my_word%exist     = .false.
  my_pdflatex%exist = .false.

  AUTHOR%string      = '?'
  AUTHOR%name        = '?'
  AUTHOR%first_name  = '?'
  AUTHOR%ADDRESS     = '?'
  AUTHOR%email       = '?'
  AUTHOR%WEB         = '?'
  AUTHOR%TEAM        = '?'

  DEVICE%string    = "?"
  DEVICE%diffracto = "?"
  DEVICE%lab       = "?"
  DEVICE%radiation = "?"
  DEVICE%wave      = "?"


  Structure_solution%name         = 'SIR97'
  Structure_solution%reference    = 'A. Altomare, M. C. Burla, M. Camalli, G. Cascarano, '// &
                                    'C. Giacovazzo, A. Guagliardi, A. G. G. Moliterni, '//   &
                                    'G. Polidori, R. Spagna, J. Appl. Cryst. (1999) 32, 115-119'
  Structure_solution%CIF_ref      = 'SIR97 (Altomare et al., 1999)'
  Structure_refinement%name       = 'SHELXL-97'
  Structure_refinement%reference  = 'Sheldrick G.M., Acta Cryst. A64 (2008), 112-122'
  Structure_refinement%CIF_ref    = 'SHELXL-97 (Sheldrick, 2008)'
  Absorption_correction%name      = '?'
  Absorption_correction%reference = '?'
  Absorption_correction%CIF_ref   = '?'

!-------------------------------------------------------------------------------

  debug_proc%unit    =  10
  debug_proc%write   = .false.
  debug_proc%level_1 = .false.
  debug_proc%level_2 = .false.
  debug_proc%level_3 = .false.
  

  ON_SCREEN          = .true.
  ON_SCREEN_PRF      = .false.
  CIFdep             = .false.
  ACTA               = .false.
  step_by_step       = .true.

  ACE_file_name       = ''
  CEL_file_name       = ''
  CIF_file_name       = ''
  CIF_pymol_file_name = ''
  input_CIF_file      = ''
  CIF_file_exist      = .false.  
  nb_input_CIF_files  = 0
  archive_CIF         = ''
  INS_file_name       = ''
  PCR_file_name       = ''
  P4P_file_name       = ''
  M50_file_name       = ''
  X_file_name         = ''
  RED_file_name       = ''
  TIDY_out_file_name  = ''
  RMAT_file_name      = ''
  RAW_file_name       = ''
  HKL_file_name       = ''
  ABS_file_name       = ''
  FACES_file_name     = ''
  main_title          = ''
  SADABS_line_ratio               = '?'
  SADABS_line_estimated_Tmin_Tmax = '?'
  SADABS_absorption_coeff         = '?'
  SADABS_ratio   = -999.
  SADABS_Tmin    = -999.
  SADABS_Tmax    = -999.

  nb_atom = 0
  nb_hkl  = 0
  nb_hkl_SFAC_calc = 0
  nb_dist_calc = 0

  CONN_dmin          = 0.5
  CONN_dmax          = 3.
  CONN_dmax_ini      = 3.
  CONN_dmin_ini      = 0.5
  CONN_all           = .false.
  CONN_all_X         = .false.  
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
  nb_ang_calc  = 0
  nb_bary_calc = 0
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

  Max_ref                = Max_allowed   ! Max. number of hkl reflections
  LOCK_wave_value        = 0.02
  CIF_format80           = .true.
  CIF_torsion_limit      = 170.
  include_RES_file       = .false.
  include_HKL_file       = .false.
  update_parameters      = .true.
  report_header          = .true.
  expert_mode            = .false.
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
  X_pattern%scale         = 100.
  X_pattern%step          = 0.01
  X_pattern%xmin          = 0
  X_pattern%xmax          = 120.


  N_PV%U                  = 0.0055
  N_PV%V                  = -0.0015
  N_PV%W                  = 0.0036
  N_PV%eta0               = 0.01
  N_PV%eta1               = 0.
  N_pattern%WDT           = 20.
  N_pattern%background    = 50.
  N_pattern%scale         = 100.
  N_pattern%step          = 0.025
  N_pattern%xmin          = 0
  N_pattern%xmax          = 140.



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
  input_CELL_P4P       = .false.
  input_CELL_M50       = .false.
  input_CELL_INS       = .false.
  input_CELL_x         = .false.
  input_CELL_rmat      = .false.
  input_CELL_CIF       = .false.
  input_CELL_RED       = .false.
  input_CELL           = .false.

  mode_interactif      = .false.
  keyword_CELL         = .false.
  keyword_NIGGLI       = .false.
  keyword_WAVE         = .false.
  keyword_WRITE_CELL   = .false.
  write_cell_cart      = .false.
  keyword_WRITE_CHEM   = .false.
  keyword_WRITE_WAVE   = .false.
  keyword_WRITE_DEVICE = .false.
  keyword_WRITE_SG     = .false.
  keyword_WRITE_BEAM   = .false.
  keyword_WRITE_QVEC   = .false.
  keyword_WRITE_ZUNIT  = .false.
  keyword_BEAM         = .false.
  keyword_SIZE         = .false.
  keyword_SPGR         = .false.
  get_SPGR             = .false.
  keyword_LSPGR        = .false.
   list_sg(1:7)        = .false.
   list_sg_Bravais(1:7)= .false.
   list_sg_centric(1:2)= .false.
   list_sg_laue(1:14)  = .false.
   list_sg_multip      = .false.
   list_sg_enantio     = .false.
   list_sg_chiral      = .false.
   list_sg_polar       = .false.
  keyword_LAUE         = .false.
  keyword_ATOM_list    = .false.
   write_atoms_cartesian = .false.   
   write_atoms_in_A      = .false.
  keyword_ADP_list     = .false.
  keyword_SFAC_UNIT    = .false.
  keyword_CONT         = .false.
  keyword_CHEM         = .false.
  keyword_ZUNIT        = .false.
  keyword_MU           = .false.
  keyword_MAT          = .false.
  keyword_LST_MAT      = .false.
  keyword_TRANSL       = .false.
  keyword_FILE         = .false.
  keyword_X_WAVE       = .false.
  keyword_P4P          = .false.
  keyword_RAW          = .false.
  keyword_HKL          = .false.
  keyword_MATMUL       = .false.
  keyword_DIAG         = .false.
  keyword_VERSION      = .false.

  keyword_INSIDE       = .false.
  keyword_PAUSE        = .false.

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
  create_FST_POLY          = .false.
  create_FST_MOLE          = .false.  
  FST_no_H                 = .false.
  launch_FP_Studio         = .false.
  keyword_create_INS       = .false.
  keyword_create_SOLVE     = .false.
  keyword_create_TIDY      = .false.  
  create_CIF_PYMOL         = .false.
  create_SHAPE_file        = .false.
  poly_vol_calc            = .false.
  INI_create_ACE           = .false.
  INI_create_CEL           = .false.
  INI_create_CFL           = .false.
  INI_create_INS           = .false.
  INI_create_FST           = .false.
  INI_create_CIF_PYMOL     = .false.
  INI_create_PRF           = .false.
  include_SQUEEZE          = .true.

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

  keyword_QVEC         = .false.
  keyword_SYMM         = .false.
  keyword_BARY         = .false.      ! calcul du barycentre
  keyword_DIST         = .false.      ! calcul de distances interatomiques
  keyword_DIST_X       = .false.
  keyword_DIST_plus    = .false.
  keyword_DHA          = .false.
  keyword_CONN         = .false.
  keyword_ANG          = .false.      ! calcul d'angles
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

  CIF_parameter%diffracto_device                 = '?'
  CIF_parameter%diffracto_radiation_type         = '?'
  CIF_parameter%diffracto_radiation_wavelength   = '?'
  CIF_parameter%diffracto_radiation_source       = '?'
  CIF_parameter%diffracto_temperature            = '?'

  CIF_parameter%diffrn_measurement_device        = '?'
  CIF_parameter%diffrn_measurement_device_type   = '?'
  CIF_parameter%diffrn_measurement_method        = '?'
  CIF_parameter%diffrn_detector                  = '?'
  CIF_parameter%diffrn_source                    = '?'
  CIF_parameter%diffrn_detector_area_resol_mean  = '?'

  CIF_parameter%symmetry_cell_setting            = '?'
  CIF_parameter%symmetry_space_group             = '?'
  CIF_parameter%symmetry_IT_number               = '?'
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
  CIF_parameter%diffrn_theta_full                = '?'
  CIF_parameter%diffrn_theta_max                 = '?'

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


  CIF_parameter%computing_data_collection        = '?'
  CIF_parameter%computing_cell_refinement        = '?'
  CIF_parameter%computing_data_reduction         = '?'
  CIF_parameter%computing_structure_solution     = '?'
  CIF_parameter%computing_structure_refinement   = '?'
  CIF_parameter%computing_molecular_graphics     = '?'
  CIF_parameter%computing_publication_material_1 = '?'
  CIF_parameter%computing_publication_material_2 = '?'


  write(CIF_parameter%computing_structure_solution,   '(3a)') "'", trim(structure_solution%CIF_ref),   "'"
  write(CIF_parameter%computing_structure_refinement, '(3a)') "'", trim(structure_refinement%CIF_ref), "'"

  CIF_parameter%computing_molecular_graphics         = "'Ortep-3 for Windows (Farrugia, 1997)'"
  !CIF_parameter%computing_publication_material_1     = "WinGX publication routines (Farrugia, 1999),"
  CIF_parameter%computing_publication_material_1     = "WinGX publication routines (Farrugia, 2012),"
  CIF_parameter%computing_publication_material_2     = "CRYSCALc (T. Roisnel, local program, 2014)"
  CIF_parameter%WinGX_used                           = .false.


  DENZO%data_collection                              = "'Collect (Nonius BV, 1997-2000)'"
  DENZO%cell_refinement                              = "'HKL Scalepack (Otwinowski & Minor 1997)'"
  DENZO%data_reduction                               = "'HKL Denzo and Scalepack (Otwinowski & Minor 1997)'"
  EVAL%data_collection                               = "'Collect (Nonius BV, 1997-2000)'"
  EVAL%cell_refinement                               = "'Dirax/lsq (Duisenberg & Schreurs, 1989-2000)'"
  EVAL%data_reduction                                = "'EvalCCD (Duisenberg & Schreurs 1990-2000)'"
  !APEX%data_collection                               = "'Bruker SMART (Bruker, 2002)'"
  !APEX%cell_refinement                               = "'Bruker SMART (Bruker, 2002)'"
  !APEX%data_reduction                                = "'Bruker SAINT (Bruker, 2002)'"
  APEX%data_collection                               = "'Bruker APEX2 (Bruker, 2006)'"
  APEX%cell_refinement                               = "'Bruker APEX2 (Bruker, 2006)'"
  APEX%data_reduction                                = "'Bruker APEX2 (Bruker, 2006)'"

  X2S%data_collection                               = "'GIS (Bruker)'"
  X2S%cell_refinement                               = "'APEX2 (Bruker, 2010); SAINT (Bruker, 2009)'"
  X2S%data_reduction                                = "'SAINT (Bruker, 2009); XPREP (Sheldrick, 2008)'"
  XCALIBUR%data_collection                          = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"
  XCALIBUR%cell_refinement                          = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"
  XCALIBUR%data_reduction                           = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"
  SUPERNOVA%data_collection                         = "'CrysAlisPro, Agilent Technologies,Version 1.171.36'"
  SUPERNOVA%cell_refinement                         = "'CrysAlisPro, Agilent Technologies,Version 1.171.36'"
  SUPERNOVA%data_reduction                          = "'CrysAlisPro, Agilent Technologies,Version 1.171.36'"


  CIF_parameter_KCCD%diffrn_source                        = "'Enraf Nonius FR590'"
  CIF_parameter_KCCD%diffrn_radiation_wavelength          = '0.71073'
  CIF_parameter_KCCD%diffrn_radiation_type                = 'MoK\a'
  CIF_parameter_KCCD%diffrn_radiation_source              = 'fine-focus sealed tube'
  CIF_parameter_KCCD%diffrn_radiation_monochromator       = 'graphite'
  CIF_parameter_KCCD%diffrn_radiation_probe               = 'x-ray'  
  CIF_parameter_KCCD%diffrn_detector                      = "'CCD plate'"
  CIF_parameter_KCCD%diffrn_detector_area_resol_mean      = '9'
  CIF_parameter_KCCD%diffrn_measurement_device            = "'95mm CCD camera on \k-goniostat'"
  CIF_parameter_KCCD%diffrn_measurement_device_type       = "'KappaCCD, Nonius'"
  CIF_parameter_KCCD%diffrn_measurement_method            = "'CCD rotation images, thin slices'"
  CIF_parameter_KCCD%computing_data_collection            = DENZO%data_collection
  CIF_parameter_KCCD%computing_cell_refinement            = DENZO%cell_refinement
  CIF_parameter_KCCD%computing_data_reduction             = DENZO%data_reduction


  CIF_parameter_APEX%diffrn_measurement_device_type       = "'APEXII, Bruker-AXS'"
  CIF_parameter_APEX%diffrn_measurement_method            = "'CCD rotation images, thin slices'"
  CIF_parameter_APEX%diffrn_detector                      = "'CCD plate'"
  CIF_parameter_APEX%diffrn_radiation_wavelength          = '0.71073'
  CIF_parameter_APEX%diffrn_radiation_type                = 'MoK\a'
  CIF_parameter_APEX%diffrn_radiation_source              = 'fine-focus sealed tube'
  CIF_parameter_APEX%diffrn_radiation_monochromator       = 'graphite'
  CIF_parameter_APEX%diffrn_radiation_probe               = 'x-ray'  
  CIF_parameter_APEX%diffrn_measurement_device            = '?'
  CIF_parameter_APEX%diffrn_source                        = '?'
  CIF_parameter_APEX%diffrn_detector_area_resol_mean      = '?'
  CIF_parameter_APEX%computing_data_collection            = APEX%data_collection
  CIF_parameter_APEX%computing_cell_refinement            = APEX%cell_refinement
  CIF_parameter_APEX%computing_data_reduction             = APEX%data_reduction

  CIF_parameter_APEX%computing_structure_solution         = CIF_parameter%computing_structure_solution
  CIF_parameter_APEX%computing_structure_refinement       = CIF_parameter%computing_structure_refinement
  CIF_parameter_APEX%computing_molecular_graphics         = CIF_parameter%computing_molecular_graphics
  CIF_parameter_APEX%computing_publication_material_1     = CIF_parameter%computing_publication_material_1
  CIF_parameter_APEX%computing_publication_material_2     = CIF_parameter%computing_publication_material_2

  CIF_parameter_X2S%diffrn_measurement_device_type       = "'Bruker SMART X2S benchtop'"
  CIF_parameter_X2S%diffrn_measurement_method            = "'omega scans'"
  CIF_parameter_X2S%diffrn_detector                      = "'CCD plate'"
  CIF_parameter_X2S%diffrn_radiation_wavelength          = '0.71073'
  CIF_parameter_X2S%diffrn_radiation_type                = 'MoK\a'
  CIF_parameter_X2S%diffrn_radiation_source              = "'XOS X-beam microfocus source'"
  CIF_parameter_X2S%diffrn_radiation_monochromator       = "'doubly curved silicon crystal'"
  CIF_parameter_X2S%diffrn_radiation_probe               = 'x-ray'  
  CIF_parameter_X2S%diffrn_measurement_device            = '?'
  CIF_parameter_X2S%diffrn_source                        = '?'
  CIF_parameter_X2S%diffrn_detector_area_resol_mean      = '?'
  CIF_parameter_X2S%computing_data_collection            = X2S%data_collection
  CIF_parameter_X2S%computing_cell_refinement            = X2S%cell_refinement
  CIF_parameter_X2S%computing_data_reduction             = X2S%data_reduction

  CIF_parameter_X2S%computing_structure_solution         = CIF_parameter%computing_structure_solution
  CIF_parameter_X2S%computing_structure_refinement       = CIF_parameter%computing_structure_refinement
  CIF_parameter_X2S%computing_molecular_graphics         = CIF_parameter%computing_molecular_graphics
  CIF_parameter_X2S%computing_publication_material_1     = CIF_parameter%computing_publication_material_1
  CIF_parameter_X2S%computing_publication_material_2     = CIF_parameter%computing_publication_material_2




  CIF_parameter_XCALIBUR%diffrn_measurement_device_type   = "'CCD Saphire 3 Xcalibur'"
  CIF_parameter_XCALIBUR%diffrn_measurement_method        = "'omega scan'"
  CIF_parameter_XCALIBUR%diffrn_detector                  = "'CCD plate'"
  CIF_parameter_XCALIBUR%diffrn_radiation_wavelength      = '0.71073'
  CIF_parameter_XCALIBUR%diffrn_radiation_monochromator   = 'graphite'
  CIF_parameter_XCALIBUR%diffrn_radiation_probe           = 'x-ray'
  CIF_parameter_XCALIBUR%diffrn_radiation_source          = 'sealed X-ray tube'
  CIF_parameter_XCALIBUR%diffrn_detector_area_resol_mean  = '19.64'
  CIF_parameter_XCALIBUR%computing_data_collection        = XCALIBUR%data_collection
  CIF_parameter_XCALIBUR%computing_cell_refinement        = XCALIBUR%cell_refinement
  CIF_parameter_XCALIBUR%computing_data_reduction         = XCALIBUR%data_reduction
  
  CIF_parameter_SUPERNOVA%diffrn_measurement_device_type   = "'SuperNova, Dual, Cu at zero, Atlas'"
  CIF_parameter_SUPERNOVA%diffrn_measurement_method        = "'omega scans'"
  CIF_parameter_SUPERNOVA%diffrn_detector                  = "'CCD plate'"
  CIF_parameter_SUPERNOVA%diffrn_radiation_wavelength      = ''
  CIF_parameter_SUPERNOVA%diffrn_radiation_monochromator   = 'multiplayer optics'
  CIF_parameter_SUPERNOVA%diffrn_radiation_probe           = 'x-ray'
  CIF_parameter_SUPERNOVA%diffrn_radiation_source          = 'X-ray microsource'

  CIF_parameter_SUPERNOVA%diffrn_detector_area_resol_mean  = '5.2474'
  CIF_parameter_SUPERNOVA%computing_data_collection        = SUPERNOVA%data_collection
  CIF_parameter_SUPERNOVA%computing_cell_refinement        = SUPERNOVA%cell_refinement
  CIF_parameter_SUPERNOVA%computing_data_reduction         = SUPERNOVA%data_reduction
  CIF_parameter_SUPERNOVA%computing_structure_solution     = CIF_parameter%computing_structure_solution
  CIF_parameter_SUPERNOVA%computing_structure_refinement   = CIF_parameter%computing_structure_refinement
  CIF_parameter_SUPERNOVA%computing_molecular_graphics     = CIF_parameter%computing_molecular_graphics
  CIF_parameter_SUPERNOVA%computing_publication_material_1 = CIF_parameter%computing_publication_material_1
  CIF_parameter_SUPERNOVA%computing_publication_material_2 = CIF_parameter%computing_publication_material_2


  SADABS%type       = "_exptl_absorpt_correction_type                       multi-scan"
  SADABS%details(1) = "_exptl_absorpt_process_details"
  SADABS%details(2) = ";"
  SADABS%details(3) = "  [Sheldrick, G.M. (2002). SADABS Bruker AXS Inc., Madison, Wisconsin, USA]"
  SADABS%details(4) = ";"
  
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


  UB_matrix     = 0.
  UB_mat_log     = .false.
  keyword_UB_mat = .false.


  crystal%size           = 0.
  crystal%size_min       = 0.
  crystal%size_max       = 0.
  crystal%size_mid       = 0.
  crystal%volume         = 0.
  crystal%radius         = 0.
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

  absorption%mu        = 0.
  absorption%Tmin      = 0.
  absorption%Tmax      = 0.


  molecule%common_name = '?'
  molecule%formula     = '?'
  molecule%content     = '?'
  molecule%weight      = 0.0
  molecule%density     = 0.0
  molecule%Z           = 0
  molecule%Z_unit      = 1
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
  symm_op_xyz             = .false.
  symm_op_mat             = .false.

  WRITE_SPG_info          = .false.
  write_SPG_info_all      = .true.
  WRITE_SPG_exti          = .false.
  write_SPG_all           = .false.
  write_SPG_subgroups     = .false.
  WRITE_symm_op           = .false.
  WRITE_SHELX_symm_op     = .false.
  WRITE_APPLY_symm        = .false.
  WRITE_STAR_K            = .false.
  WRITE_SITE_info         = .false.
  WRITE_PCR_site_info     = .false.
  WRITE_PCR_mag_site_info = .false.
  WRITE_triclinic_transf  = .false.
  WRITE_monoclinic_transf = .false.
  WRITE_rhomb_hex_transf  = .false.
  WRITE_hex_rhomb_transf  = .false.
  WRITE_permutation_abc   = .false.
  WRITE_twin_hexa         = .false.
  WRITE_twin_pseudo_hexa  = .false.
  keyword_LST_MAT         = .false.

  keyword_THERM           = .false.
  keyword_THERM_SHELX     = .false.
  THERM_Uiso              = .false.
  THERM_Biso              = .false.
  THERM_Uij               = .false.
  THERM_Bij               = .false.
  THERM_Beta              = .false.
  THERM_aniso             = .false.


  keyword_KEY             = .false.
  keyword_HELP            = .false.
  keyword_CLA             = .false.
  keyword_HEADER          = .false.
  keyword_SYST            = .false.
  keyword_RESET           = .false.

  keyword_create_CRYSCALC_news   = .false.
  keyword_create_CRYSCALC_HTML   = .false.
  browse_cryscalc_HTML           = .false.

  keyword_WRITE_REF_APEX         = .false.
  keyword_WRITE_REF_DENZO        = .false.
  keyword_WRITE_REF_EVAL         = .false.
  keyword_write_REF_KCCD         = .false.  
  keyword_WRITE_REF_SADABS       = .false.
  keyword_WRITE_REF_ABS_CRYSALIS = .false.
  keyword_write_REF_SUPERNOVA    = .false.
  keyword_WRITE_REF_X2S          = .false.
  keyword_write_REF_XCALIBUR     = .false.
  
  
  
  unknown_keyword              = .false.
  unknown_CFL_keyword          = .false.
  keyword_SEARCH_EXTI          = .false.
  keyword_SEARCH_SPGR          = .false.
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

  keyword_DATA_neutrons = .false.
  data_neutrons_PLOT    = .false.

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
    'CHEM               ', 'CONN               ', 'CONT               ', 'CREATE_ACE         ', 'CREATE_CEL         ', &
    'CREATE_CFL         ', 'CREATE_FST         ', 'CREATE_INS         ', 'CREATE_REPORT      ', 'CREATE_SOLVE       ', &
	'CREATE_TIDY        ', 'D_HKL              ', 'D_STAR             ', 'DATA_ATOMIC_DENSITY', 'DATA_ATOMIC_RADIUS ', &
	'DATA_ATOMIC_WEIGHT ', 'DATA_NEUTRONS      ', 'DATA_XRAYS         ', 'DIAG_MAT           ', 'DIR_ANG            ', &
	'DIST               ', 'DIST_DHA           ', 'EDIT               ', 'EQUIV_HKL          ', 'EXIT               ', &
	'FILE               ', 'FIND_HKL           ', 'FIND_HKL_LIST      ', 'FRIEDEL            ', 'GEN_HKL            ', &
	'HEADER             ', 'HEX_RHOMB          ', 'HKL                ', 'HKL_NEG            ', 'HKL_POS            ', &
	'INSIDE             ', 'LIST_EXTI          ', 'LIST_KEYS          ', 'LIST_LAUE          ', 'LIST_MATR          ', &
	'LIST_SG            ', 'MAG                ', 'MAN                ', 'MAN_HTML           ', 'MATMUL             ', &
	'MATR               ', 'OBV_REV            ', 'MENDEL             ', 'MERGE              ', 'MONOCLINIC         ', &
	'NEWS               ', 'NIGGLI             ', 'P4P                ', 'PAUSE              ', 'PERMUT             ', &
	'Q_HKL              ', 'QVEC               ', 'READ_CEL           ', 'READ_CIF           ', 'READ_FACES         ', &
	'READ_INS           ', 'READ_NREPORT       ', 'READ_PCR           ', 'READ_TIDY_OUT      ', 'REC_ANG            ', &
	'REF_ABS_CRYSALIS   ', 'REF_APEX           ', 'REF_DENZO          ', 'REF_EVAL           ', 'REF_KCCD           ', &
	'REF_SADABS         ', 'REF_SUPERNOVA      ', 'REF_X2S            ', 'REF_XCALIBUR       ', 'RESET              ', &
	'RINT               ', 'RHOMB_HEX          ', 'SEARCH_EXTI        ', 'SEARCH_SPGR        ', 'SET                ', &
	'SETTING            ', 'SFAC               ', 'SF_HKL             ', 'SG                 ', 'SG_ALL             ', &
	'SG_EXTI            ', 'SG_INFO            ', 'SG_SUB             ', 'SHANNON            ', 'SHELL              ', &
	'SHIFT_2TH          ', 'SITE_INFO          ', 'SIZE               ', 'SORT               ', 'STAR_K             ', &
	'STL                ', 'SYMM               ', 'SYST               ', 'THERM              ', 'THERM_SHELX        ', &
	'THETA              ', 'TITL               ', 'TRANSLATION        ', 'TRANSMISSION       ', 'UB_MATRIX          ', &
	'USER_MAT           ', 'TRICLINIC          ', 'TWIN_HEXA          ', 'TWIN_PSEUDO_HEXA   ', 'TWO_THETA          ', &
	'UNIT               ', 'WAVE               ', 'WEB                ', 'WRITE_ADP          ', 'WRITE_BEAM         ', &
	'WRITE_CELL         ', 'WRITE_CHEM         ', 'WRITE_DEVICE       ', 'WRITE_QVEC         ', 'WRITE_SYM_OP       ', &
	'WRITE_WAVE         ', 'WRITE_ZUNIT        ', 'X_WAVE             ', 'ZUNIT              ', 'WRITE_SG           '  /)

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
  numor = numor + 1;    HELP_CHEM_numor                =  numor  ! 11
  numor = numor + 1;    HELP_CONN_numor                =  numor  ! 12
  numor = numor + 1;    HELP_CONT_numor                =  numor  ! 13
  numor = numor + 1;    HELP_CREATE_ACE_numor          =  numor  ! 14
  numor = numor + 1;    HELP_CREATE_CEL_numor          =  numor  ! 15
  numor = numor + 1;    HELP_CREATE_CFL_numor          =  numor  ! 16
  numor = numor + 1;    HELP_CREATE_FST_numor          =  numor  ! 17
  numor = numor + 1;    HELP_CREATE_INS_numor          =  numor  ! 18
  numor = numor + 1;    HELP_CREATE_REPORT_numor       =  numor  ! 19
  numor = numor + 1;    HELP_CREATE_SOLVE_numor        =  numor  ! 20 
  numor = numor + 1;    HELP_CREATE_TIDY_numor         =  numor  ! 21
  numor = numor + 1;    HELP_D_HKL_numor               =  numor  ! 22
  numor = numor + 1;    HELP_D_STAR_numor              =  numor  ! 23
  numor = numor + 1;    HELP_DATA_ATOMIC_DENSITY_numor =  numor  ! 24
  numor = numor + 1;    HELP_DATA_ATOMIC_RADIUS_numor  =  numor  ! 25
  numor = numor + 1;    HELP_DATA_ATOMIC_WEIGHT_numor  =  numor  ! 26
  numor = numor + 1;    HELP_DATA_NEUTRONS_numor       =  numor  ! 27
  numor = numor + 1;    HELP_DATA_XRAYS_numor          =  numor  ! 28
  numor = numor + 1;    HELP_DIAG_MAT_numor            =  numor  ! 29
  numor = numor + 1;    HELP_DIR_ANG_numor             =  numor  ! 30
  numor = numor + 1;    HELP_DIST_numor                =  numor  ! 31
  numor = numor + 1;    HELP_DIST_DHA_numor            =  numor  ! 32
  numor = numor + 1;    HELP_EDIT_numor                =  numor  ! 33
  numor = numor + 1;    HELP_EQUIV_numor               =  numor  ! 34
  numor = numor + 1;    HELP_EXIT_numor                =  numor  ! 35
  numor = numor + 1;    HELP_FILE_numor                =  numor  ! 36
  numor = numor + 1;    HELP_FIND_HKL_numor            =  numor  ! 37
  numor = numor + 1;    HELP_FIND_HKL_LIST_numor       =  numor  ! 38
  numor = numor + 1;    HELP_FRIEDEL_pairs_numor       =  numor  ! 39
  numor = numor + 1;    HELP_GEN_HKL_numor             =  numor  ! 40
  numor = numor + 1;    HELP_HEADER_numor              =  numor  ! 41
  numor = numor + 1;    HELP_HEX_RHOMB_numor           =  numor  ! 42
  numor = numor + 1;    HELP_HKL_numor                 =  numor  ! 43
  numor = numor + 1;    HELP_HKL_NEG_numor             =  numor  ! 44
  numor = numor + 1;    HELP_HKL_POS_numor             =  numor  ! 45
  numor = numor + 1;    HELP_INSIDE_numor              =  numor  ! 46
  numor = numor + 1;    HELP_LIST_EXTI_numor           =  numor  ! 47
  numor = numor + 1;    HELP_LIST_KEYS_numor           =  numor  ! 48
  numor = numor + 1;    HELP_LIST_LAUE_numor           =  numor  ! 49
  numor = numor + 1;    HELP_LIST_MATR_numor           =  numor  ! 50
  numor = numor + 1;    HELP_LIST_SG_numor             =  numor  ! 51
  numor = numor + 1;    HELP_MAG_numor                 =  numor  ! 52
  numor = numor + 1;    HELP_MAN_numor                 =  numor  ! 53
  numor = numor + 1;    HELP_MAN_HTML_numor            =  numor  ! 54
  numor = numor + 1;    HELP_MATMUL_numor              =  numor  ! 55
  numor = numor + 1;    HELP_MATR_numor                =  numor  ! 56
  numor = numor + 1;    HELP_MENDEL_numor              =  numor  ! 57
  numor = numor + 1;    HELP_MERGE_numor               =  numor  ! 58
  numor = numor + 1;    HELP_MONOCLINIC_numor          =  numor  ! 59
  numor = numor + 1;    HELP_NEWS_numor                =  numor  ! 60
  numor = numor + 1;    HELP_NIGGLI_CELL_numor         =  numor  ! 61
  numor = numor + 1;    HELP_OBV_REV_numor             =  numor  ! 62
  numor = numor + 1;    HELP_P4P_numor                 =  numor  ! 63
  numor = numor + 1;    HELP_PAUSE_numor               =  numor  ! 64
  numor = numor + 1;    HELP_PERMUT_numor              =  numor  ! 65
  numor = numor + 1;    HELP_Q_HKL_numor               =  numor  ! 66
  numor = numor + 1;    HELP_QVEC_numor                =  numor  ! 67
  numor = numor + 1;    HELP_READ_CEL_numor            =  numor  ! 68
  numor = numor + 1;    HELP_READ_CIF_numor            =  numor  ! 69
  numor = numor + 1;    HELP_READ_FACES_numor          =  numor  ! 70
  numor = numor + 1;    HELP_READ_INS_numor            =  numor  ! 71
  numor = numor + 1;    HELP_READ_NREPORT_numor        =  numor  ! 72
  numor = numor + 1;    HELP_READ_PCR_numor            =  numor  ! 73
  numor = numor + 1;    HELP_READ_TIDY_out_numor       =  numor  ! 74
  numor = numor + 1;    HELP_REC_ANG_numor             =  numor  ! 75 
  numor = numor + 1;    HELP_REF_ABS_CRYSALIS_numor    =  numor  ! 76
  numor = numor + 1;    HELP_REF_APEX_numor            =  numor  ! 77
  numor = numor + 1;    HELP_REF_DENZO_numor           =  numor  ! 78
  numor = numor + 1;    HELP_REF_EVAL_numor            =  numor  ! 79
  numor = numor + 1;    HELP_REF_KCCD_numor            =  numor  ! 80
  numor = numor + 1;    HELP_REF_SADABS_numor          =  numor  ! 81
  numor = numor + 1;    HELP_REF_SUPERNOVA_numor       =  numor  ! 82
  numor = numor + 1;    HELP_REF_X2S_numor             =  numor  ! 83
  numor = numor + 1;    HELP_REF_XCALIBUR_numor        =  numor  ! 84
  numor = numor + 1;    HELP_RESET_numor               =  numor  ! 85
  numor = numor + 1;    HELP_RINT_numor                =  numor  ! 86
  numor = numor + 1;    HELP_RHOMB_HEX_numor           =  numor  ! 87
  numor = numor + 1;    HELP_SEARCH_EXTI_numor         =  numor  ! 88
  numor = numor + 1;    HELP_SEARCH_SPGR_numor         =  numor  ! 89
  numor = numor + 1;    HELP_SET_numor                 =  numor  ! 90
  numor = numor + 1;    HELP_SETTING_numor             =  numor  ! 91
  numor = numor + 1;    HELP_SFAC_numor                =  numor  ! 92
  numor = numor + 1;    HELP_SFHKL_numor               =  numor  ! 93
  numor = numor + 1;    HELP_SG_numor                  =  numor  ! 94
  numor = numor + 1;    HELP_SG_ALL_numor              =  numor  ! 95
  numor = numor + 1;    HELP_SG_EXTI_numor             =  numor  ! 96
  numor = numor + 1;    HELP_SG_INFO_numor             =  numor  ! 97
  numor = numor + 1;    HELP_SG_SUB_numor              =  numor  ! 98
  numor = numor + 1;    HELP_SHANNON_numor             =  numor  ! 99
  numor = numor + 1;    HELP_SHELL_numor               =  numor  !100
  numor = numor + 1;    HELP_SHIFT_2TH_numor           =  numor  !101
  numor = numor + 1;    HELP_SITE_INFO_numor           =  numor  !102
  numor = numor + 1;    HELP_SIZE_numor                =  numor  !103
  numor = numor + 1;    HELP_SORT_numor                =  numor  !104 
  numor = numor + 1;    HELP_STAR_K_numor              =  numor  !105
  numor = numor + 1;    HELP_STL_numor                 =  numor  !106
  numor = numor + 1;    HELP_SYMM_numor                =  numor  !107
  numor = numor + 1;    HELP_SYST_numor                =  numor  !108
  numor = numor + 1;    HELP_THERM_numor               =  numor  !109
  numor = numor + 1;    HELP_THERM_SHELX_numor         =  numor  !110
  numor = numor + 1;    HELP_THETA_numor               =  numor  !111
  numor = numor + 1;    HELP_TITL_numor                =  numor  !112
  numor = numor + 1;    HELP_TRANSLATION_numor         =  numor  !113
  numor = numor + 1;    HELP_TRANSMISSION_numor        =  numor  !114
  numor = numor + 1;    HELP_TRICLINIC_numor           =  numor  !115
  numor = numor + 1;    HELP_TWIN_HEXA_numor           =  numor  !116
  numor = numor + 1;    HELP_TWIN_PSEUDO_HEXA_numor    =  numor  !117
  numor = numor + 1;    HELP_TWO_THETA_numor           =  numor  !118
  numor = numor + 1;    HELP_UB_matrix_numor           =  numor  !119
  numor = numor + 1;    HELP_UNIT_numor                =  numor  !120
  numor = numor + 1;    HELP_USER_MAT_numor            =  numor  !121
  numor = numor + 1;    HELP_WAVE_numor                =  numor  !122
  numor = numor + 1;    HELP_WEB_numor                 =  numor  !123
  numor = numor + 1;    HELP_WRITE_ADP_numor           =  numor  !124
  numor = numor + 1;    HELP_WRITE_BEAM_numor          =  numor  !125
  numor = numor + 1;    HELP_WRITE_CELL_numor          =  numor  !126
  numor = numor + 1;    HELP_WRITE_CHEM_numor          =  numor  !127
  numor = numor + 1;    HELP_WRITE_DEVICE_numor        =  numor  !128
  numor = numor + 1;    HELP_WRITE_QVEC_numor          =  numor  !129
  numor = numor + 1;    HELP_WRITE_SG_numor            =  numor  !130
  numor = numor + 1;    HELP_WRITE_SYM_OP_numor        =  numor  !131
  numor = numor + 1;    HELP_WRITE_WAVE_numor          =  numor  !132
  numor = numor + 1;    HELP_WRITE_ZUNIT_numor         =  numor  !133
  numor = numor + 1;    HELP_X_wave_numor              =  numor  !134
  numor = numor + 1;    HELP_ZUNIT_numor               =  numor  !135



 end subroutine cryscalc_init
 !------------------------------------------------------

 subroutine read_cryscalc_ini()
  USE cryscalc_module, ONLY : debug_proc, INI_unit, cryscalc, winplotr_exe, my_editor, my_browser, my_word, &
                              my_pdflatex, WEB, AUTHOR, DEVICE, cryscalc,                                       &                              
                              wavelength, keyword_beam, keyword_WAVE, neutrons, X_rays,                         &
                              DENZO, EVAL, APEX, X2S, SUPERNOVA, XCALIBUR,                                      &
                              CONN_dmax_ini, CONN_dmax, CONN_dmin_ini, CONN_dmin,                               &
							  CIF_format80, include_RES_file, include_HKL_file, update_parameters,              &
                              report_header, LOCK_wave_value,  Max_ref,  expert_mode, hkl_statistics,           &
                              hkl_format, cartesian_frame, keep_bond_str_out,                                   &
                              skip_start_menu, pdp_simu, CIF_torsion_limit,                                       &
                              keyword_create_ACE, keyword_create_CEL, keyword_create_CFL, keyword_create_FST,   &
                              keyword_create_INS, keyword_create_PRF,                                           &
                              create_CIF_PYMOL, INI_create_CIF_PYMOL,                                           &
                              INI_create_ACE,     INI_create_CEL,     INI_create_CFL,     INI_create_FST,       &
                              INI_create_INS, INI_create_PRF,                                                   &
                              structure_solution, structure_refinement, absorption_correction,                  &
                              Create_INS, get_sample_ID,                    &
                              message_text, winplotr_path_name, on_screen, on_screen_prf
  USE HKL_module,     ONLY : n_sig, threshold, MAX_allowed
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
  inquire(file='cryscalc.ini', exist = file_exist)
  if(file_exist) then
   cryscalc%ini = 'cryscalc.ini'
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
    IF(read_line(1:7) == '[AUTHOR' .or. read_line(1:6) == '[USER]') then
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


   if(u_case(DEVICE%radiation(1:1)) == 'X') then
    neutrons = .false.
    X_rays   = .true.
    call get_X_radiation(u_case(DEVICE%radiation))
    keyword_beam = .true.
   ELSEIF(u_case(DEVICE%radiation(1:)) == 'NEUT') then
     neutrons = .true.
     X_rays   = .false.
     keyword_beam = .true.
   endif


   ! -------------------- sept 2010  ---------------------------------
   if(u_case(DEVICE%radiation(1:4)) == 'X_Mo') then
    CIF_parameter%diffracto_radiation_type = 'MoK\a'
   elseif(u_case(DEVICE%radiation(1:4)) == 'X_Cu') then
    CIF_parameter%diffracto_radiation_type = 'CuK\a'
   endif
   CIF_parameter%diffracto_radiation_wavelength = DEVICE%wave
   ! ------------------------------------------------------------------

   if(u_case(DEVICE%diffracto(1:4)) == 'APEX') then
    !CIF_parameter = CIF_parameter_APEX
    CIF_parameter%diffrn_measurement_device_type       = CIF_parameter_APEX%diffrn_measurement_device_type
    CIF_parameter%diffrn_measurement_method            = CIF_parameter_APEX%diffrn_measurement_method
    CIF_parameter%diffrn_detector                      = CIF_parameter_APEX%diffrn_detector
    CIF_parameter%diffrn_radiation_wavelength          = CIF_parameter_APEX%diffrn_radiation_wavelength
    CIF_parameter%diffrn_radiation_type                = CIF_parameter_APEX%diffrn_radiation_type
    CIF_parameter%diffrn_radiation_source              = CIF_parameter_APEX%diffrn_radiation_source
    CIF_parameter%diffrn_radiation_monochromator       = CIF_parameter_APEX%diffrn_radiation_monochromator
    CIF_parameter%diffrn_radiation_probe               = CIF_parameter_APEX%diffrn_radiation_probe
    CIF_parameter%computing_data_collection            = CIF_parameter_APEX%computing_data_collection
    CIF_parameter%computing_cell_refinement            = CIF_parameter_APEX%computing_cell_refinement
    CIF_parameter%computing_data_reduction             = CIF_parameter_APEX%computing_data_reduction
    !CIF_parameter%computing_structure_solution         = CIF_parameter_APEX%computing_structure_solution
    !CIF_parameter%computing_structure_refinement       = CIF_parameter_APEX%computing_structure_refinement
    !CIF_parameter%computing_molecular_graphics         = CIF_parameter_APEX%computing_molecular_graphics
    !CIF_parameter%computing_publication_material       = CIF_parameter_APEX%computing_publication_material


   elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
    !CIF_parameter = CIF_parameter_KCCD
    CIF_parameter%diffrn_source                        = CIF_parameter_KCCD%diffrn_source
    CIF_parameter%diffrn_radiation_wavelength          = CIF_parameter_KCCD%diffrn_radiation_wavelength
    CIF_parameter%diffrn_radiation_type                = CIF_parameter_KCCD%diffrn_radiation_type
    CIf_parameter%diffrn_radiation_source              = CIF_parameter_KCCD%diffrn_radiation_source
    CIF_parameter%diffrn_radiation_monochromator       = CIF_parameter_KCCD%diffrn_radiation_monochromator
    CIF_parameter%diffrn_radiation_probe               = CIF_parameter_KCCD%diffrn_radiation_probe
    CIF_parameter%diffrn_detector                      = CIF_parameter_KCCD%diffrn_detector
    CIF_parameter%diffrn_detector_area_resol_mean      = CIF_parameter_KCCD%diffrn_detector_area_resol_mean
    CIF_parameter%diffrn_measurement_device            = CIF_parameter_KCCD%diffrn_measurement_device
    CIF_parameter%diffrn_measurement_device_type       = CIF_parameter_KCCD%diffrn_measurement_device_type
    CIF_parameter%diffrn_measurement_method            = CIF_parameter_KCCD%diffrn_measurement_method
    CIF_parameter%computing_data_collection            = CIF_parameter_KCCD%computing_data_collection
    CIF_parameter%computing_cell_refinement            = CIF_parameter_KCCD%computing_cell_refinement
    CIF_parameter%computing_data_reduction             = CIF_parameter_KCCD%computing_data_reduction

    !DENZO%data_collection                              = CIF_parameter_KCCD%computing_data_collection
    !DENZO%cell_refinement                              = CIF_parameter_KCCD%computing_cell_refinement
    !DENZO%data_reduction                               = CIF_parameter_KCCD%computing_data_reduction

    !CIF_parameter%computing_data_collection            = "'Collect (Nonius BV, 1997-2000)'"
    !CIF_parameter%computing_cell_refinement            = "'Dirax/lsq (Duisenberg & Schreurs, 1989-2000)'"
    !CIF_parameter%computing_data_reduction             = "'EvalCCD (Duisenberg & Schreurs 1990-2000)'"
    !EVAL%data_collection                               = CIF_parameter%computing_data_collection
    !EVAL%cell_refinement                               = CIF_parameter%computing_cell_refinement
    !EVAL%data_reduction                                = CIF_parameter%computing_data_reduction


   elseif(u_case(DEVICE%diffracto(1:7)) == 'SAPHIRE' .or. u_case(DEVICE%diffracto(1:9)) == 'EXCALIBUR') then
    !CIF_parameter = CIF_parameter_XCALIBUR
    CIF_parameter%diffrn_measurement_device_type       = CIF_parameter_XCALIBUR%diffrn_measurement_device_type
    CIF_parameter%diffrn_measurement_method            = CIF_parameter_XCALIBUR%diffrn_measurement_method
    CIF_parameter%diffrn_detector                      = CIF_parameter_XCALIBUR%diffrn_detector
    CIF_parameter%diffrn_detector_area_resol_mean      = CIF_parameter_XCALIBUR%diffrn_detector_area_resol_mean
    CIF_parameter%computing_data_collection            = CIF_parameter_XCALIBUR%computing_data_collection
    CIF_parameter%computing_cell_refinement            = CIF_parameter_XCALIBUR%computing_cell_refinement
    CIF_parameter%computing_data_reduction             = CIF_parameter_XCALIBUR%computing_data_reduction

    XCALIBUR%data_collection                           = CIF_parameter%computing_data_collection
    XCALIBUR%cell_refinement                           = CIF_parameter%computing_cell_refinement
    XCALIBUR%data_reduction                            = CIF_parameter%computing_data_reduction

   elseif(u_case(DEVICE%diffracto(1:9)) == 'SUPERNOVA' ) then
    !CIF_parameter = CIF_parameter_XCALIBUR
    CIF_parameter%diffrn_measurement_device_type       = CIF_parameter_SUPERNOVA%diffrn_measurement_device_type
    CIF_parameter%diffrn_measurement_method            = CIF_parameter_SUPERNOVA%diffrn_measurement_method
    CIF_parameter%diffrn_detector                      = CIF_parameter_SUPERNOVA%diffrn_detector
    CIF_parameter%diffrn_detector_area_resol_mean      = CIF_parameter_SUPERNOVA%diffrn_detector_area_resol_mean
    CIF_parameter%computing_data_collection            = CIF_parameter_SUPERNOVA%computing_data_collection
    CIF_parameter%computing_cell_refinement            = CIF_parameter_SUPERNOVA%computing_cell_refinement
    CIF_parameter%computing_data_reduction             = CIF_parameter_SUPERNOVA%computing_data_reduction

    SUPERNOVA%data_collection                          = CIF_parameter%computing_data_collection
    SUPERNOVA%cell_refinement                          = CIF_parameter%computing_cell_refinement
    SUPERNOVA%data_reduction                           = CIF_parameter%computing_data_reduction

	
   elseif(u_case(DEVICE%diffracto(1:3)) == 'X2S') then
    !CIF_parameter = CIF_parameter_APEX
    CIF_parameter%diffrn_measurement_device_type       = CIF_parameter_X2S%diffrn_measurement_device_type
    CIF_parameter%diffrn_measurement_method            = CIF_parameter_X2S%diffrn_measurement_method
    CIF_parameter%diffrn_detector                      = CIF_parameter_X2S%diffrn_detector
    CIF_parameter%diffrn_radiation_wavelength          = CIF_parameter_X2S%diffrn_radiation_wavelength
    CIF_parameter%diffrn_radiation_type                = CIF_parameter_X2S%diffrn_radiation_type
    CIF_parameter%diffrn_radiation_source              = CIF_parameter_X2S%diffrn_radiation_source
    CIF_parameter%diffrn_radiation_monochromator       = CIF_parameter_X2S%diffrn_radiation_monochromator
    CIF_parameter%diffrn_radiation_probe               = CIF_parameter_X2S%diffrn_radiation_probe
    CIF_parameter%computing_data_collection            = CIF_parameter_X2S%computing_data_collection
    CIF_parameter%computing_cell_refinement            = CIF_parameter_X2S%computing_cell_refinement
    CIF_parameter%computing_data_reduction             = CIF_parameter_X2S%computing_data_reduction

   endif


   if(debug_proc%write) then
    call write_debug_proc ('DEVICE_diffractometer', trim(DEVICE%diffracto))
    call write_debug_proc ('DEVICE_lab',            trim(DEVICE%lab))
    call write_debug_proc ('DEVICE_radiation',      trim(DEVICE%radiation))
    call write_debug_proc ('DEVICE_wave',           trim(DEVICE%wave))
    if(u_case(DEVICE%diffracto(1:4)) == 'APEX') then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call WRITE_debug_proc('MEASUREMENT_detector',      trim(CIF_parameter%diffrn_detector))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))

    elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call WRITE_debug_proc('MEASUREMENT_device',        trim(CIF_parameter%diffrn_measurement_device))
     call WRITE_debug_proc('DIFFRACTION_source',        trim(CIF_parameter%diffrn_source))
     call WRITE_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter%diffrn_detector))
     call WRITE_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter%diffrn_detector_area_resol_mean))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))

    elseif(u_case(DEVICE%diffracto(1:7)) == 'SAPHIRE' .or. u_case(DEVICE%diffracto(1:9)) == 'EXCALIBUR') then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call WRITE_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter%diffrn_detector))
     call WRITE_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter%diffrn_detector_area_resol_mean))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))

    elseif(u_case(DEVICE%diffracto(1:9)) == 'SUPERNOVA') then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call WRITE_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter%diffrn_detector))
     call WRITE_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter%diffrn_detector_area_resol_mean))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))

    elseif(u_case(DEVICE%diffracto(1:3)) == 'X2S') then
     call WRITE_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call WRITE_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call WRITE_debug_proc('MEASUREMENT_detector',      trim(CIF_parameter%diffrn_detector))
     call WRITE_debug_proc('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call WRITE_debug_proc('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call WRITE_debug_proc('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_source))
     call WRITE_debug_proc('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call WRITE_debug_proc('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))


    endif
    call WRITE_debug_proc('COMPUTING_data_collection',  trim(CIF_parameter%computing_data_collection))
    call WRITE_debug_proc('COMPUTING_cell_refinement',  trim(CIF_parameter%computing_cell_refinement))
    call WRITE_debug_proc('COMPUTING_data_reduction',   trim(CIF_parameter%computing_data_reduction))

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
       if(INI_real < MAX_allowed) Max_ref = INI_real
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
       INI_create_ACE     = .true.
       keyword_create_ACE = .true.
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_cel' .and. answer_yes(trim(INI_string))) then
       INI_create_CEL     = .true.
       keyword_create_CEL = .true.
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_cfl' .and.answer_yes(trim(INI_string))) then
       INI_create_CFL     = .true.
       keyword_create_CFL = .true.
       cycle
      ELSEIF(L_case(read_line(1:10)) == 'create_fst' .AND. answer_yes(trim(INI_string))) then
       INI_create_FST     = .true.
       keyword_create_FST = .true.
       cycle
      ELSEIF(L_case(read_line(1:16)) == 'create_cif_pymol' .AND. answer_yes(trim(INI_string))) then
       INI_create_CIF_PYMOL      = .true.
       create_CIF_PYMOL  = .true.
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_ins' .and. answer_yes(trim(INI_string))) then
       INI_create_INS     = .true.
       keyword_create_INS = .true.
       cycle
	  ELSEIF(l_case(read_line(1:14)) == 'create_pat_prf' .and. answer_yes(trim(INI_string))) then
       INI_create_PRF     = .true.
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
       include_HKL_file = .true.
       cycle
      ELSEIF(l_case(read_line(1:17)) == 'update_parameters' .and. answer_no(trim(INI_string))) then
       update_parameters = .false.
       cycle
      ELSEIF((l_case(read_line(1:18)) == 'html_report_header' .or. l_case(read_line(1:13)) == 'report_header') &
             .and. answer_no(trim(INI_string))) then
       report_header = .false.
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
	   call open_debug_file
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
      elseif(l_case(read_line(1:17)) == 'CIF_torsion_limit') then
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
      if (ABS(INI_real) < 180. ) CIF_torsion_limit =  INI_real
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
       cycle
      ELSEIF(l_case(read_line(1:28)) == 'structure_solution_reference') then
       structure_solution%reference = INI_string
       cycle
      ELSEIF(l_case(read_line(1:26)) == 'structure_solution_cif_ref') then
       structure_solution%cif_ref = INI_string
       cycle
      ELSEIF(l_case(read_line(1:25)) == 'structure_refinement_name') then
       structure_refinement%name = INI_string
       cycle
      ELSEIF(l_case(read_line(1:30)) == 'structure_refinement_reference') then
       structure_refinement%reference = INI_string
       cycle
      ELSEIF(l_case(read_line(1:28)) == 'structure_refinement_cif_ref') then
       structure_refinement%cif_ref = INI_string
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

 end subroutine  edit_a_file

!-----------------------------------------------------
subroutine output_setting
 USE cryscalc_module, ONLY :  INI_unit, cryscalc, my_browser, my_editor, my_pdflatex, WEB, DEVICE, AUTHOR, message_text,        &
                              CONN_dmax_ini, INI_create_ACE, INI_create_CEL, INI_create_CFL, INI_create_FST, INI_create_INS,    &
							  INI_create_PRF, INI_create_CIF_PYMOL,                                                             & 
                              structure_solution, structure_refinement, absorption_correction, Create_INS, max_ref,             &
                              lock_wave_value, CIF_format80, include_RES_file, include_HKL_file,                                &
							  update_parameters , report_header, expert_mode,                                                   &
                              hkl_statistics, hkl_format, skip_start_menu, pdp_simu, CIF_torsion_limit, debug_proc,             &
							  cartesian_frame, keep_bond_str_out 
 USE pattern_profile_module
 USE CIF_module
 USE USER_module
 USE MATRIX_list_module, only : user_mat_text, user_mat_nb, transf_mat, max_mat_nb
 USE HKL_module,     ONLY :  n_sig, threshold
 USE IO_module,      ONLY :  write_info
 implicit none
  integer                 :: i

 IF(LEN_TRIM(cryscalc%ini) == 0) then
   call write_info('')
   call write_info('  No CRYSCALC.ini setting file has been defined.')
   call write_info('')
   return
 endif

 call write_info('')
 call write_info(' [EXTERNAL APPLICATIONS]')
 IF(my_editor%exist)    then
  call write_info('   > EDITOR   : '//TRIM(my_editor%name))
 else
  call write_info('   > EDITOR   :   not defined.')
 endif
 IF(my_browser%exist)   then
  call write_info('   > BROWSER  : '//TRIM(my_browser%name))
 else
  call write_info('   > BROWSER  :   not defined.')
 endif
 IF(my_pdflatex%exist)   then
  call write_info('   > PDFLATEX : '//TRIM(my_pdflatex%name))
 else
  call write_info('   > PDFLATEX :   not defined.')
 endif
 call write_info('')

 call write_info(' [WEB ADDRESS]')
 if(WEB%num_site /= 0) then
  do i = 1, WEB%num_site
   write(message_text, '(a,a10,2a)') '   > ',WEB%name(i), ' = ', trim(WEB%ADDRESS(i))
   call write_info(trim(message_text))
  end do
 endif
 call write_info('')

 call write_info(' [DEVICE]')
 call write_info('   > diffractometer = '// trim(DEVICE%diffracto))
 call write_info('   > lab            = '// trim(DEVICE%lab))
 call write_info('   > radiation      = '// trim(DEVICE%radiation))
 call write_info('   > wave_A         = '// trim(DEVICE%wave))
 call write_info('   > temperature    = '// trim(CIF_cell_measurement%temperature))
 call write_info('')

 call write_info(' [USER]')
 call write_info('   > name           = '//trim(AUTHOR%name))
 call write_info('   > first_name     = '//trim(AUTHOR%first_name))
 call write_info('   > address        = '//trim(AUTHOR%ADDRESS))
 call write_info('   > email          = '//trim(AUTHOR%email))
 call write_info('   > web            = '//trim(AUTHOR%WEB))
 call write_info('   > team           = '//trim(AUTHOR%team))
 call write_info('')

 call write_info(' [ARRAYS DIMENSIONS]')
 write(message_text, '(I8)') Max_ref
 call write_info('   > hkl_reflections = '//   trim(message_text))
 call write_info('')

 call write_info(' [PARAMETERS]')
 write(message_text, '(F5.2)') n_sig
 call write_info('   > I_sig          = '//   trim(message_text))
 write(message_text, '(F5.2)') threshold
 call write_info('   > threshold      ='//  trim(message_text))
 write(message_text, '(F5.2)') CONN_dmax_ini
 call write_info('   > d_max_A        = '//   trim(message_text))
 call write_info('')

 call write_info(' [CREATE INS]')
 write(message_text, '(F8.2)') Create_INS%temperature
 call write_info('   > Temperature    = '//   trim(message_text))
 write(message_text, '(F4.1)') Create_INS%U_threshold
 call write_info('   > U_threshold    ='//  trim(message_text))
 call write_info('')


 call write_info(' [COMMAND LINE ARGUMENT]')
 if(INI_create_ACE) then
  call write_info('   > create_ACE    = 1')
 else
  call write_info('   > create_ACE    = 0')
 endif
 if(INI_create_CEL) then
  call write_info('   > create_CEL       = 1')
 else
  call write_info('   > create_CEL       = 0')
 endif
 if(INI_create_CFL) then
  call write_info('   > create_CFL       = 1')
 else
  call write_info('   > create_CFL       = 0')
 endif
 if(INI_create_FST) then
  call write_info('   > create_FST       = 1')
 else
  call write_info('   > create_FST       = 0')
 endif
 if(INI_create_CIF_PYMOL) then
  call write_info('   > create_CIF_PYMOL = 1')
 else
  call write_info('   > create_CIF_PYMOL = 0')
 endif
 if(INI_create_INS) then
  call write_info('   > create_INS       = 1')
 else
  call write_info('   > create_INS       = 0')
 endif
 if(INI_create_PRF) then
  call write_info('   > create_PAT_PRF   = 1')
 else
  call write_info('   > create_PAT_PRF   = 0')
 endif
 call write_info('')

 call write_info(' [PROGRAMS]')
 call write_info('   > structure_solution_name         = '// trim(structure_solution%name))
 call write_info('   > structure_solution_reference    = '// trim(structure_solution%reference))
 call write_info('   > structure_solution_cif_ref      = '// trim(structure_solution%cif_ref))
 call write_info('   > structure_refinement_name       = '// trim(structure_refinement%name))
 call write_info('   > structure_refinement_reference  = '// trim(structure_refinement%reference))
 call write_info('   > structure_refinement_cif_ref    = '// trim(structure_refinement%cif_ref))
 call write_info('   > absorption_correction_name      = '// trim(Absorption_correction%name))
 call write_info('   > absorption_correction_reference = '// trim(Absorption_correction%reference))
 call write_info('   > absorption_correction_cif_ref   = '// trim(Absorption_correction%cif_ref))
 call write_info('')


 call write_info(' [OPTIONS]')
 write(message_text, '(a,F5.2,a)') '   > LOCK_wave_value            = ', lock_wave_value,    &
                                   '       ! lock wavelength to anticathode value'
 call write_info(trim(message_text))
 if(CIF_format80) then
  call write_info('   > CIF_format80               = 1           ! formatted line, when creating a CIF file, ')
 else
  call write_info('   > CIF_format80               = 0           ! formatted line, when creating a CIF file, ')
 endif
 call write_info('                                                if more than 80 characters')
 write(message_text, '(a,F6.2,a)') '   > CIF_torsion_limit          = ', CIF_torsion_limit, &
                                   '      ! exclude torsion angle if greater'
 call write_info(trim(message_text))
 if(include_RES_file) then
  call write_info('   > include_RES_file           = 1           ! include .RES file in the archive_cryscalc.cif file')
 else
  call write_info('   > include_RES_file           = 0           ! include .RES file in the archive_cryscalc.cif file')
 endif
 if(include_HKL_file) then
  call write_info('   > include_HKL_file           = 1           ! include .HKL file in the archive_cryscalc.cif file')
 else
  call write_info('   > include_HKL_file           = 0           ! include .HKL file in the archive_cryscalc.cif file')
 endif
 if(update_parameters) then
  call write_info('   > update_parameters          = 1           ! update parameters after transformation (cell parameters,' // &
                                                                ' atomic coordinates) ')
 else
  call write_info('   > update_parameters          = 0           ! update parameters after transformation (cell parameters,' // &
                                                                 ' atomic coordinates) ')
 endif

 if(report_header) then
  call write_info('   > report_header              = 1           ! Write header in structural report')
 else
  call write_info('   > report_header              = 0           ! Write header in structural report')
 endif

 if(expert_mode) then
  call write_info('   > expert_mode                = 1           ! Expert mode activated')
  if(keep_bond_str_out) then
  call write_info('   > bond_str.out               = 1           ! keep bond_str.out')
  else
  call write_info('   > bond_str.out               = 0           ! remove bond_str.out')
  endif
 else
  call write_info('   > expert_mode                = 0           ! Expert mode OFF')
 endif
 if(debug_proc%level_2) then
  call write_info('   > debug_mode level 2         = 1           ! Debug mode (level 2) activated')
 elseif(debug_proc%level_3) then
  call write_info('   > debug_mode level 3         = 1           ! Debug mode (level 3) activated')
 endif
 if(skip_start_menu) then
  call write_info('   > skip_start_menu            = 1           ! Skip start menu')
 else
  call write_info('   > skip_start_menu            = 0           ! Start menu ON')
 endif

 if(hkl_statistics) then
  call write_info('   > hkl_statistics             = 1           ! Output statistics on hkl reflections')
 else
  call write_info('   > hkl_statistics             = 0           ! Output statistics on hkl reflections')
 endif
  call write_info('   > hkl_format                 = '//trim(hkl_format) // '   ! format for .hkl file (h,k,lF2,sig)')
  call write_info('   > Cartesian frame type       = '//trim(cartesian_frame%type) // '   ! A: x//A ; C:x//c')  

 if(pdp_simu%cu) then
  call write_info('   > pdp_cu                     = 1           ! Ka Cu for powder diffraction pattern simulation')
 else
  call write_info('   > pdp_cu                     = 0           ! Current wavelength for powder diffraction pattern simulation')
 endif
 if(pdp_simu%beam == 'neutrons') then
  call write_info('   > pdp_beam                   = N           ! ' & 
                //'Simulation of powder diffraction pattern : N for neutrons / X for X-rays')
 else
  call write_info('   > pdp_beam                   = X           ! ' &
                //'Simulation of powder diffraction pattern : N for neutrons / X for X-rays')
 endif
 write(message_text, '(a,F10.5,a)') '   > pdp_wave                   = ', pdp_simu%wave, &
                                    '  ! Wavelength used for powder diffraction pattern simulation'
 call write_info(trim(message_text))
 
 call write_info('')
 call write_info('[PATTERN SIMULATION (Pseudo-Voigt function)]')
 call write_info(' X-ray pattern profile :')
 write(message_text, '(a,F10.5,a)') '   > X_profile_U            = ', X_PV%U,    &
                                   '       ! U value of the Cagliotti formula : FWHM2 = U*TAN**2(theta) + V*TAN(theta) + W'
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > X_profile_V            = ', X_PV%V
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > X_profile_W            = ', X_PV%W
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5,a)') '   > X_profile_eta0         = ', X_PV%eta0, &
                                   '       ! Lorentzian components : eta = eta° + 2theta * eta1'
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > X_profile_eta1         = ', X_PV%eta1
 call write_info(trim(message_text))
 call write_info(' X-ray Diffraction pattern :')
 write(message_text, '(a,F10.5)')   '   > X_Pattern_background     = ', X_pattern%background
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')  '   > X_Pattern_scale          = ', X_pattern%scale
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > X_Pattern_step           = ', X_pattern%step
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > X_Pattern_xmin           = ', X_pattern%xmin
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > X_Pattern_xmax           = ', X_pattern%xmax
 call write_info(trim(message_text))

 call write_info(' Neutron pattern profile :')
 write(message_text, '(a,F10.5,a)') '   > N_profile_U            = ', N_PV%U,    &
                                   '       ! U value of the Cagliotti formula : FWHM2 = U*TAN**2(theta) + V*TAN(theta) + W'
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > N_profile_V            = ', N_PV%V
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > N_profile_W            = ', N_PV%W
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5,a)') '   > N_profile_eta0         = ', N_PV%eta0, &
                                   '       ! Lorentzian components : eta = eta0 + 2theta * eta1'
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > N_profile_eta1         = ', N_PV%eta1
 call write_info(trim(message_text))
 call write_info(' Neutron diffraction pattern :')
 write(message_text, '(a,F10.5)')   '   > N_Pattern_background     = ', N_pattern%background
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')  '   > N_Pattern_scale          = ', N_pattern%scale
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > N_Pattern_step           = ', N_pattern%step
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > N_Pattern_xmin           = ', N_pattern%xmin
 call write_info(trim(message_text))
 write(message_text, '(a,F10.5)')   '   > N_Pattern_xmax           = ', N_pattern%xmax
 call write_info(trim(message_text))




 if(user_mat_nb /=0) then
  call write_info('')
  call write_info('[USER TRANSFORMATION MATRICES]')
  do i=1 , user_mat_nb
   write (message_text,'(3x,a,i1,22x,a,3(2x,3F3.0),2a)') '> MAT_',i, '=', transf_mat(:,:,max_mat_nb+i), '   !  ', &
                                                       trim(user_mat_text(i))
  call write_info(trim(message_text))
  enddo
 end if
 
 if(debug_proc%write) then
  call write_info('')
  call write_info('[USER SHORTCUTS]')
  do i=1, nb_shortcuts   
   write(message_text, '(3a)')          shortcut_kw(i)(1:20)," = " , trim(shortcut_details(i))	 
   call write_info(trim(message_text))
  end do
 endif



end subroutine output_setting

!--------------------------------------------------------------------------------


subroutine get_X_radiation(input_string)
 use wavelength_module
 use cryscalc_module, only : wavelength, keyword_WAVE, neutrons, X_rays

 implicit none
  character (len=*), intent(in)   :: input_string
  integer                         :: long_input_string
  logical                         :: target_Ag, target_Co, target_Cr, target_Cu, target_Fe, target_Mo, target_Ni

  target_Ag = .false.
  target_Co = .false.
  target_Cr = .false.
  target_Cu = .false.
  target_Fe = .false.
  target_Mo = .false.
  target_Ni = .false.

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
      wavelength =  X_target(1)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(1)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Mo) then
      wavelength =  X_target(2)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(2)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Cu) then
      wavelength =  X_target(3)%wave(1)
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
      wavelength =  X_target(5)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(5)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Fe) then
      wavelength =  X_target(6)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(6)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(target_Cr) then
      wavelength =  X_target(7)%wave(1)
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
  character (len=*), intent(in), optional :: value
  integer                                 :: long, i
  character (len=24)                      :: fmt_


  long = len_trim(field)

  if(.not. present(value)) then
   write(debug_proc%unit, '(a)') trim(field)

  else
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
  endif

   return
end subroutine write_debug_proc

!-------------------------------------------------------------------
subroutine write_debug_proc_level(level, field)
 use cryscalc_module, only : debug_proc
 use macros_module,   only : u_case, l_case
 implicit none
  integer,           intent(in)           :: level
  character (len=*), intent(in)           :: field
  integer                                 :: long
  character (len=24)                      :: fmt_


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

    if(u_case(DEVICE%diffracto(1:4)) == 'APEX') then
     call write_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call write_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call write_debug_proc('MEASUREMENT_detector',      trim(CIF_parameter%diffrn_detector))
     call write_debug_proc('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call write_debug_proc('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call write_debug_proc('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_source))
     call write_debug_proc('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call write_debug_proc('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))

    elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
     call write_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call write_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call write_debug_proc('MEASUREMENT_device',        trim(CIF_parameter%diffrn_measurement_device))
     call write_debug_proc('DIFFRACTION_source',        trim(CIF_parameter%diffrn_source))
     call write_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter%diffrn_detector))
     call write_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter%diffrn_detector_area_resol_mean))
     call write_debug_proc('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call write_debug_proc('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call write_debug_proc('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_type))
     call write_debug_proc('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call write_debug_proc('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))

    elseif(u_case(DEVICE%diffracto(1:7)) == 'SAPHIRE' .or. u_case(DEVICE%diffracto(1:9)) == 'EXCALIBUR') then
     call write_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call write_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call write_debug_proc('DIFFRACTION_detector',      trim(CIF_parameter%diffrn_detector))
     call write_debug_proc('DIFFRACTION_detector_area', trim(CIF_parameter%diffrn_detector_area_resol_mean))

    elseif(u_case(DEVICE%diffracto(1:3)) == 'X2S') then
     call write_debug_proc('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call write_debug_proc('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call write_debug_proc('MEASUREMENT_detector',      trim(CIF_parameter%diffrn_detector))
     call write_debug_proc('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call write_debug_proc('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call write_debug_proc('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_source))
     call write_debug_proc('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call write_debug_proc('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))
    endif

    call write_debug_proc('COMPUTING_data_collection',  trim(CIF_parameter%computing_data_collection))
    call write_debug_proc('COMPUTING_cell_refinement',  trim(CIF_parameter%computing_cell_refinement))
    call write_debug_proc('COMPUTING_data_reduction',   trim(CIF_parameter%computing_data_reduction))

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

 !CIF_parameter = CIF_parameter_APEX
 CIF_parameter%diffrn_measurement_device_type       = CIF_parameter_APEX%diffrn_measurement_device_type
 CIF_parameter%diffrn_measurement_method            = CIF_parameter_APEX%diffrn_measurement_method
 CIF_parameter%diffrn_detector                      = CIF_parameter_APEX%diffrn_detector
 CIF_parameter%diffrn_radiation_wavelength          = CIF_parameter_APEX%diffrn_radiation_wavelength
 CIF_parameter%diffrn_radiation_type                = CIF_parameter_APEX%diffrn_radiation_type
 CIF_parameter%diffrn_radiation_source              = CIF_parameter_APEX%diffrn_radiation_type
 CIF_parameter%diffrn_radiation_monochromator       = CIF_parameter_APEX%diffrn_radiation_monochromator
 CIF_parameter%diffrn_radiation_probe               = CIF_parameter_APEX%diffrn_radiation_probe
 CIF_parameter%computing_data_collection            = CIF_parameter_APEX%computing_data_collection
 CIF_parameter%computing_cell_refinement            = CIF_parameter_APEX%computing_cell_refinement
 CIF_parameter%computing_data_reduction             = CIF_parameter_APEX%computing_data_reduction
 CIF_parameter%computing_structure_solution         = CIF_parameter_APEX%computing_structure_solution
 CIF_parameter%computing_structure_refinement       = CIF_parameter_APEX%computing_structure_refinement
 CIF_parameter%computing_molecular_graphics         = CIF_parameter_APEX%computing_molecular_graphics
 !CIF_parameter%computing_publication_material       = CIF_parameter_APEX%computing_publication_material
 CIF_parameter%computing_publication_material_1     = CIF_parameter_APEX%computing_publication_material_1
 CIF_parameter%computing_publication_material_2     = CIF_parameter_APEX%computing_publication_material_2

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
   call WRITE_debug_proc ('MEASUREMENT_device_type',        trim(CIF_parameter%diffrn_measurement_device_type))
   call WRITE_debug_proc ('MEASUREMENT_method',             trim(CIF_parameter%diffrn_measurement_method))
   call WRITE_debug_proc ('MEASUREMENT_detector',           trim(CIF_parameter%diffrn_detector))
   call WRITE_debug_proc ('RADIATION_wavelength',           trim(CIF_parameter%diffrn_radiation_wavelength))
   call WRITE_debug_proc ('RADIATION_type',                 trim(CIF_parameter%diffrn_radiation_type))
   call WRITE_debug_proc ('RADIATION_source',               trim(CIF_parameter%diffrn_radiation_source))
   call WRITE_debug_proc ('RADIATION_monochromator',        trim(CIF_parameter%diffrn_radiation_monochromator))
   call WRITE_debug_proc ('RADIATION_probe',                trim(CIF_parameter%diffrn_radiation_probe))
   call WRITE_debug_proc ('COMPUTING_data_collection',      trim(CIF_parameter%computing_data_collection))
   call WRITE_debug_proc ('COMPUTING_cell_refinement',      trim(CIF_parameter%computing_cell_refinement))
   call WRITE_debug_proc ('COMPUTING_data_reduction',       trim(CIF_parameter%computing_data_reduction))
   call WRITE_debug_proc ('COMPUTING_structure_solution',   trim(CIF_parameter%computing_structure_solution))
   call WRITE_debug_proc ('COMPUTING_structure_refinement', trim(CIF_parameter%computing_structure_refinement))
   call WRITE_debug_proc ('COMPUTING_molecular_graphics',   trim(CIF_parameter%computing_molecular_graphics))
   !call WRITE_debug_proc ('COMPUTING_publication_material', trim(CIF_parameter%computing_publication_material))
   call WRITE_debug_proc ('COMPUTING_publication_material', trim(CIF_parameter%computing_publication_material_1))
   call WRITE_debug_proc ('COMPUTING_publication_material', trim(CIF_parameter%computing_publication_material_2))



 endif

 return
end subroutine CIF_default_values


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
    call write_DEBUG_proc ('Update_parameters' ,    '1     ! ')
    call write_DEBUG_proc ('Report_header'     ,    '1     ! ')
    call write_DEBUG_proc ('Expert_mode'       ,    '0     ! ')
    call write_DEBUG_proc ('Skip_start_menu'   ,    '0     !')
    call write_DEBUG_proc ('hkl_statistics'    ,    '1     ! ')

   call write_DEBUG_proc ('hkl_format'         ,     trim(hkl_format)//'  !')
   call write_DEBUG_proc ('Ka1 Cu for pdp simulation'     ,    '0     ! ')



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
 use cryscalc_module, only : debug_proc
 use io_module
 implicit none
  integer             :: i_error
 
 
  open(unit = debug_proc%unit, file = 'cryscalc_debug.txt', status = "replace", iostat=i_error)
	if(i_error /=0) then
     call write_info('')
     call write_info('  !! Problem opening cryscalc_debug.txt file !!')
     call write_info('')	
	 debug_proc%write = .false.
	end if

	return
 end subroutine open_debug_file
 
	