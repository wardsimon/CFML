!     Last change:  TR   11 Oct 2007    6:11 pm
subroutine cryscal_init()
 USE cryscal_module
 USE wavelength_module
 USE HKL_module


  implicit none
   INTEGER   :: numor, long

 ! initialisation ---------------------------------------
  winplotr_exe     = ''
  cryscal_ini      = ''

  my_editor%name    = ''
  my_browser%name   = ''
  my_word%name      = ''

! initialisation des champs lus dans le fichier de config. CRYSCAL.INI
  my_editor%exist   = .false.
  my_browser%exist  = .false.
  my_word%exist     = .false.

  AUTHOR%string      = '?'
  AUTHOR%name        = '?'
  AUTHOR%first_name  = '?'
  AUTHOR%ADDRESS     = '?'
  AUTHOR%email       = '?'
  AUTHOR%WEB         = '?'

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

  DEBUG_file%unit  =  10
  DEBUG_file%write = .false.
  ON_SCREEN        = .true.
  CIFdep           = .false.
  ACTA             = .false.

  ACE_file_name  = ''
  CEL_file_name  = ''
  CIF_file_name  = ''
  archive_CIF    = ''
  INS_file_name  = ''
  PCR_file_name  = ''
  P4P_file_name  = ''
  M50_file_name  = ''
  X_file_name    = ''
  RMAT_file_name = ''
  ABS_file_name  = ''
  SADABS_line    = ''
  SADABS_ratio   = -999.

  nb_atom = 0
  nb_hkl  = 0
  nb_hkl_SFAC_calc = 0
  nb_dist_calc = 0
  
  CONN_dmax    = 3.
  CONN_all     = .false.
  DIST_coef    = 0.
  nb_ang_calc  = 0
  nb_bary_calc = 0
  nb_sort      = 0
  nb_shell     = 0
  nb_symm_op   = 0
  wavelength   = 0.

  n_sig         = 3.
  threshold     = 0.03
  CONN_dmax_ini = 3.
  
  Create_INS_temperature = -999.
  Create_INS_u_threshold = -999.

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
  input_CELL           = .false.

  mode_interactif      = .false.
  keyword_CELL         = .false.
  keyword_NIGGLI       = .false.
  keyword_WAVE         = .false.
  keyword_WRITE_CELL   = .false.
  keyword_WRITE_CHEM   = .false.
  keyword_WRITE_WAVE   = .false. 
  keyword_WRITE_DEVICE = .false. 
  keyword_WRITE_SG     = .false.
  keyword_WRITE_BEAM   = .false.
  keyword_WRITE_QVEC   = .false.
  keyword_BEAM         = .false.
  keyword_SIZE         = .false.
  keyword_SPGR         = .false.
  keyword_LSPGR        = .false.
   list_sg(1:7)        = .false.
   list_sg_Bravais(1:7)= .false.
   list_sg_centric(1:2)= .false.
   list_sg_laue(1:14)  = .false.
   list_sg_multip      = .false.
  keyword_LAUE         = .false.
  keyword_ATOM_list    = .false.
  keyword_SFAC_UNIT    = .false.
  keyword_CONT         = .false.
  keyword_CHEM         = .false.
  keyword_ZUNIT        = .false.
  keyword_MATR         = .false.
  keyword_LST_MAT      = .false.
  keyword_TRANSL       = .false.
  keyword_FILE         = .false.
  keyword_X_WAVE       = .false.
  keyword_P4P          = .false.
  keyword_RAW          = .false.
  keyword_MATMUL       = .false.

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
  M95_comment_line     = ''
  M95_comment_line_nb  = 0

  keyword_STRUCT_factor= .false.

  keyword_read_CEL     = .false.
  keyword_read_CIF     = .false.
  keyword_read_INS     = .false.
  keyword_read_PCR     = .false.

  keyword_create_ACE   = .false.  
  keyword_create_CEL   = .false.
  keyword_create_CFL   = .false.
  keyword_create_CIF   = .false.
  keyword_create_FST   = .false.
  INI_create_ACE       = .false.
  INI_create_CEL       = .false.
  INI_create_CFL       = .false.
  INI_create_INS       = .false.
  INI_create_FST       = .false.

  keyword_create_REPORT = .false.
  keyword_read_NREPORT  = .false.
  keyword_modif_ARCHIVE = .false.
  keyword_SOLVE_to_INS  = .false.
  keyword_SIR           = .false.
  keyword_NEWS          = .false.
  news_year             = 'all'

  keyword_EDIT         = .false.
  file_to_edit         = ''
  keyword_SET          = .false.
  keyword_SETTING      = .false.

  keyword_QVEC         = .false.
  keyword_SYMM         = .false.
  keyword_BARY         = .false.      ! calcul du barycentre
  keyword_DIST         = .false.      ! calcul de distances interatomiques
  keyword_DIST_        = .false.
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

  keyword_SH_2th       = .false.
  SHIFT_2theta         = 0.

  
  HKL_2THETA           = .false.
  HKL_THETA            = .false.
  HKL_STL              = .false.
  HKL_Q                = .false.
  HKL_D                = .false.
  
  X_min                = 0.
  X_max                = 0.

  CIF_diffrn_reflns%theta_min       = 0.
  CIF_diffrn_reflns%theta_max       = 0.
  CIF_cell_measurement%theta_min    = 0.
  CIF_cell_measurement%theta_max    = 0.
  CIF_cell_measurement%reflns_used  = 0
  CIF_cell_measurement%temperature  = '?'

  CIF_dist%n_text          = 0
  CIF_dist%max_text        = 2000

  ! CIF parameters (extraits de archive.cif et utilises pour generer le fichier structural_REPORT.HTML
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
  CIF_parameter%exptl_density                    = '?'
  CIF_parameter%exptl_mu                         = '?'
  CIF_parameter%diffrn_reflns_number             = '?'
  CIF_parameter%diffrn_reflns_av_R_equivalents   = '?'
  CIF_parameter%reflns_number_total              = '?'
  CIF_parameter%reflns_number_gt                 = '?'
  CIF_parameter%refine_ls_number_parameters      = '?'
  CIF_parameter%refine_ls_wR_factor_gt           = '?'
  CIF_parameter%refine_ls_R_factor_gt            = '?'
  CIF_parameter%refine_ls_R_factor_all           = '?'
  CIF_parameter%refine_ls_wR_factor_ref          = '?'
  CIF_parameter%refine_diff_density_max          = '?'
  CIF_parameter%refine_diff_density_min          = '?'
  CIF_parameter%H_treatment                      = 'constr'
  CIF_parameter%atom                             = '?'
  CIF_parameter%distance                         = '?'
  CIF_parameter%angle                            = '?'
  CIF_parameter%torsion_angle                    = '?'
  CIF_parameter%Hbond                            = '?'
  CIF_parameter%theta_min                        = '?'
  CIF_parameter%theta_max                        = '?'
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
  CIF_parameter%completeness                     = '?'
  CIF_parameter%absorption_correction_type       = '?'
  CIF_parameter%T_min                            = '?'
  CIF_parameter%T_max                            = '?'
  CIF_parameter%restraints_number                = '?'
  CIF_parameter%CHI2                             = '?'
  CIF_parameter%crystal_system                   = '?'
  CIF_parameter%Bravais                          = '?'
  
  
  write(CIF_parameter%computing_structure_solution,   '(3a)') "'", trim(structure_solution%CIF_ref),   ")'"  
  write(CIF_parameter%computing_structure_refinement, '(3a)') "'", trim(structure_refinement%CIF_ref), ")'"
                                                           
  CIF_parameter%computing_molecular_graphics         = "'Ortep-3 for Windows (Farrugia, 1997)'"
  CIF_parameter%computing_publication_material_1     = " WinGX publication routines (Farrugia, 1999),"
  CIF_parameter%computing_publication_material_2     = " CRYSCAL (T. Roisnel, local program)"

  
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
  CRYSALIS%data_collection                           = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"
  CRYSALIS%cell_refinement                           = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"
  CRYSALIS%data_reduction                            = "'CrysAlis CCD, Oxford Diffraction Ltd.,Version 1.171.26'"

  
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

  
  CIF_parameter_XCALIBUR%diffrn_measurement_device_type   = "'CCD Saphire 3 Xcalibur'"
  CIF_parameter_XCALIBUR%diffrn_measurement_method        = "'omega scan'"
  CIF_parameter_XCALIBUR%diffrn_detector                  = "'CCD plate'"
  CIF_parameter_XCALIBUR%diffrn_radiation_wavelength      = '0.71073'
  CIF_parameter_XCALIBUR%diffrn_radiation_monochromator   = 'graphite'
  CIF_parameter_XCALIBUR%diffrn_radiation_probe           = 'x-ray'
  CIF_parameter_XCALIBUR%diffrn_detector_area_resol_mean  = '19.64'
  CIF_parameter_XCALIBUR%computing_data_collection        = CRYSALIS%data_collection
  CIF_parameter_XCALIBUR%computing_cell_refinement        = CRYSALIS%cell_refinement
  CIF_parameter_XCALIBUR%computing_data_reduction         = CRYSALIS%data_reduction

  SADABS%type       = "_exptl_absorpt_correction_type                              multi-scan"
  SADABS%details(1) = "_exptl_absorpt_process_details"
  SADABS%details(2) = ";"
  SADABS%details(3) = "  [Sheldrick, G.M. (2002). SADABS Bruker AXS Inc., Madison, Wisconsin, USA]"
  SADABS%details(4) = ";"
 
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


  UB_matrix = 0.


  crystal%size         = 0.
  crystal%size_min     = 0.
  crystal%size_max     = 0.
  crystal%size_mid     = 0.
  crystal%volume       = 0.
  crystal%radius       = 0.
  crystal%morph        = '?'
  crystal%color        = '?'
  crystal_system       = '?'
  crystal%face_line(:) = '?'
  crystal%faces_nb     = 0

  absorption%mu        = 0.
  absorption%Tmin      = 0.
  absorption%Tmax      = 0.


  molecule%common_name = '?'
  molecule%formula     = '?'
  molecule%content     = '?'
  molecule%weight      = 0.0
  molecule%density     = 0.0




  F000                 = 0.

  beam_type               = 'x_rays'
  X_rays                  = .true.
  neutrons                = .false.
  electrons               = .false.

  Mat_integer             = .false.
  sort_plot               = .false.
  sort_out                = .false.
  sort_out_n              = 1
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
  WRITE_APPLY_symm        = .false.
  WRITE_SITE_info         = .false.
  WRITE_triclinic_transf  = .false.
  WRITE_monoclinic_transf = .false.
  WRITE_rhomb_hex_transf  = .false.
  WRITE_hex_rhomb_transf  = .false.
  WRITE_permutation_abc   = .false.
  WRITE_twin_hexa         = .false.
  WRITE_twin_pseudo_hexa  = .false.
  keyword_LST_MAT         = .false.

  keyword_THERM           = .false.
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

  keyword_create_CRYSCAL_news  = .false.
  keyword_create_CRYSCAL_HTML  = .false.
  browse_cryscal_HTML          = .false.
  keyword_create_CIF           = .false.
  keyword_write_REF_KCCD       = .false.
  keyword_WRITE_REF_APEX       = .false.
  keyword_WRITE_REF_EVAL       = .false.
  keyword_WRITE_REF_DENZO      = .false.
  keyword_WRITE_REF_SADABS     = .false.
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
  
  ordered_HKL             = .false.
  keyword_Rint            = .false.
  keyword_merge_HKL       = .false.

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
  (/'ABSENT_HKL         ', 'ACTA               ', 'ANG                ', 'APPLY_OP           ', 'ATOM               ', &
    'ATOM_LIST          ', 'BARY               ', 'BEAM               ', 'CELL               ', 'CHEM               ', &
    'CONN               ', 'CONT               ', 'CREATE_ACE         ', 'CREATE_CEL         ', 'CREATE_CFL         ', &
    'CREATE_FST         ', 'CREATE_INS         ', 'CREATE_REPORT      ', 'D_HKL              ', 'D_STAR             ', &
    'DATA_ATOMIC_DENSITY', 'DATA_ATOMIC_RADIUS ', 'DATA_ATOMIC_WEIGHT ', 'DATA_NEUTRONS      ', 'DATA_XRAYS         ', &
    'DIR_ANG            ', 'DIST               ', 'EDIT               ', 'EQUIV_HKL          ', 'EXIT               ', &
    'FILE               ', 'FIND_HKL           ', 'FIND_HKL_LIST      ', 'GEN_HKL            ', 'HEADER             ', &
    'HEX_RHOMB          ', 'HKL                ', 'HKL_NEG            ', 'HKL_POS            ', 'INSIDE             ', &
    'LIST_EXTI          ', 'LIST_KEYS          ', 'LIST_LAUE          ', 'LIST_MATR          ', 'LIST_SG            ', &
    'MAG                ', 'MAN                ', 'MAN_HTML           ', 'MATMUL             ', 'MATR               ', &
    'OBV_REV            ', 'MENDEL             ', 'MERGE              ', 'MONOCLINIC         ', 'NEWS               ', &
    'NIGGLI             ', 'P4P                ', 'PAUSE              ', 'PERMUT             ', 'Q_HKL              ', &
    'QVEC               ', 'READ_CEL           ', 'READ_CIF           ', 'READ_INS           ', 'READ_PCR           ', &
    'READ_NREPORT       ', 'REC_ANG            ', 'REF_APEX           ', 'REF_DENZO          ', 'REF_EVAL           ', &
    'REF_KCCD           ', 'REF_SADABS         ', 'RESET              ', 'RINT               ', 'RHOMB_HEX          ', &
    'SEARCH_EXTI        ', 'SEARCH_SPGR        ', 'SET                ', 'SETTING            ', 'SFAC               ', &
    'SF_HKL             ', 'SG                 ', 'SG_ALL             ', 'SG_EXTI            ', 'SG_INFO            ', &
    'SG_SUB             ', 'SHANNON            ', 'SHELL              ', 'SHIFT_2TH          ', 'SITE_INFO          ', &
    'SIZE               ', 'SORT               ', 'STL                ', 'SYMM               ', 'SYST               ', &
    'THERM              ', 'THETA              ', 'TITL               ', 'TRANSLATION        ', 'TRANSMISSION       ', &
    'TRICLINIC          ', 'TWIN_HEXA          ', 'TWIN_PSEUDO_HEXA   ', 'TWO_THETA          ', 'UNIT               ', &
    'WEB                ', 'WRITE_BEAM         ', 'WRITE_CELL         ', 'WRITE_CHEM         ', 'WRITE_DEVICE       ', &
    'WRITE_SG           ', 'WRITE_SYM_OP       ', 'WRITE_WAVE         ', 'X_WAVE             ', 'ZUNIT              ', &
    'WAVE               ', 'WRITE_QVEC         '/)
						
  HELP_arg(1:nb_help_max) = HELP_string(1:nb_help_max)

  numor = 1;            HELP_ABSENT_HKL_numor          =  numor  !  1
  numor = numor + 1;    HELP_ACTA_numor                =  numor  !  2
  numor = numor + 1;    HELP_ANG_numor                 =  numor  !  3
  numor = numor + 1;    HELP_APPLY_OP_numor            =  numor  !  4
  numor = numor + 1;    HELP_ATOM_numor                =  numor  !  5
  numor = numor + 1;    HELP_ATOM_LIST_numor           =  numor  !  6
  numor = numor + 1;    HELP_BARY_numor                =  numor  !  7
  numor = numor + 1;    HELP_BEAM_numor                =  numor  !  8
  numor = numor + 1;    HELP_CELL_numor                =  numor  !  9
  numor = numor + 1;    HELP_CHEM_numor                =  numor  ! 10
  numor = numor + 1;    HELP_CONN_numor                =  numor  ! 11
  numor = numor + 1;    HELP_CONT_numor                =  numor  ! 12
  numor = numor + 1;    HELP_CREATE_ACE_numor          =  numor  ! 13  
  numor = numor + 1;    HELP_CREATE_CEL_numor          =  numor  ! 14
  numor = numor + 1;    HELP_CREATE_CFL_numor          =  numor  ! 15
  numor = numor + 1;    HELP_CREATE_FST_numor          =  numor  ! 16
  numor = numor + 1;    HELP_CREATE_INS_numor          =  numor  ! 17 
  numor = numor + 1;    HELP_CREATE_REPORT_numor       =  numor  ! 19
  numor = numor + 1;    HELP_D_HKL_numor               =  numor  ! 19
  numor = numor + 1;    HELP_D_STAR_numor              =  numor  ! 20
  numor = numor + 1;    HELP_DATA_ATOMIC_DENSITY_numor =  numor  ! 21
  numor = numor + 1;    HELP_DATA_ATOMIC_RADIUS_numor  =  numor  ! 22
  numor = numor + 1;    HELP_DATA_ATOMIC_WEIGHT_numor  =  numor  ! 23  
  numor = numor + 1;    HELP_DATA_NEUTRONS_numor       =  numor  ! 24
  numor = numor + 1;    HELP_DATA_XRAYS_numor          =  numor  ! 25
  numor = numor + 1;    HELP_DIR_ANG_numor             =  numor  ! 26
  numor = numor + 1;    HELP_DIST_numor                =  numor  ! 27
  numor = numor + 1;    HELP_EDIT_numor                =  numor  ! 28
  numor = numor + 1;    HELP_EQUIV_numor               =  numor  ! 29
  numor = numor + 1;    HELP_EXIT_numor                =  numor  ! 30
  numor = numor + 1;    HELP_FILE_numor                =  numor  ! 31
  numor = numor + 1;    HELP_FIND_HKL_numor            =  numor  ! 32
  numor = numor + 1;    HELP_FIND_HKL_LIST_numor       =  numor  ! 33  
  numor = numor + 1;    HELP_GEN_HKL_numor             =  numor  ! 34
  numor = numor + 1;    HELP_HEADER_numor              =  numor  ! 35
  numor = numor + 1;    HELP_HEX_RHOMB_numor           =  numor  ! 36
  numor = numor + 1;    HELP_HKL_numor                 =  numor  ! 37 
  numor = numor + 1;    HELP_HKL_NEG_numor             =  numor  ! 38
  numor = numor + 1;    HELP_HKL_POS_numor             =  numor  ! 39
  numor = numor + 1;    HELP_INSIDE_numor              =  numor  ! 40
  numor = numor + 1;    HELP_LIST_EXTI_numor           =  numor  ! 41
  numor = numor + 1;    HELP_LIST_KEYS_numor           =  numor  ! 42
  numor = numor + 1;    HELP_LIST_LAUE_numor           =  numor  ! 43
  numor = numor + 1;    HELP_LIST_MATR_numor           =  numor  ! 44
  numor = numor + 1;    HELP_LIST_SG_numor             =  numor  ! 45
  numor = numor + 1;    HELP_MAG_numor                 =  numor  ! 46
  numor = numor + 1;    HELP_MAN_numor                 =  numor  ! 47
  numor = numor + 1;    HELP_MAN_HTML_numor            =  numor  ! 48
  numor = numor + 1;    HELP_MATMUL_numor              =  numor  ! 49
  numor = numor + 1;    HELP_MATR_numor                =  numor  ! 50
  numor = numor + 1;    HELP_MENDEL_numor              =  numor  ! 51
  numor = numor + 1;    HELP_MERGE_numor               =  numor  ! 52
  numor = numor + 1;    HELP_MONOCLINIC_numor          =  numor  ! 53
  numor = numor + 1;    HELP_NEWS_numor                =  numor  ! 54
  numor = numor + 1;    HELP_NIGGLI_CELL_numor         =  numor  ! 55
  numor = numor + 1;    HELP_OBV_REV_numor             =  numor  ! 56
  numor = numor + 1;    HELP_P4P_numor                 =  numor  ! 57
  numor = numor + 1;    HELP_PAUSE_numor               =  numor  ! 58
  numor = numor + 1;    HELP_PERMUT_numor              =  numor  ! 59
  numor = numor + 1;    HELP_Q_HKL_numor               =  numor  ! 60
  numor = numor + 1;    HELP_QVEC_numor                =  numor  ! 61
  numor = numor + 1;    HELP_READ_CEL_numor            =  numor  ! 62
  numor = numor + 1;    HELP_READ_CIF_numor            =  numor  ! 63
  numor = numor + 1;    HELP_READ_INS_numor            =  numor  ! 64
  numor = numor + 1;    HELP_READ_PCR_numor            =  numor  ! 65
  numor = numor + 1;    HELP_READ_NREPORT_numor        =  numor  ! 66
  numor = numor + 1;    HELP_REC_ANG_numor             =  numor  ! 67
  numor = numor + 1;    HELP_REF_APEX_numor            =  numor  ! 68
  numor = numor + 1;    HELP_REF_DENZO_numor           =  numor  ! 69
  numor = numor + 1;    HELP_REF_EVAL_numor            =  numor  ! 70
  numor = numor + 1;    HELP_REF_KCCD_numor            =  numor  ! 71
  numor = numor + 1;    HELP_REF_SADABS_numor          =  numor  ! 72
  numor = numor + 1;    HELP_RESET_numor               =  numor  ! 73
  numor = numor + 1;    HELP_RINT_numor                =  numor  ! 74
  numor = numor + 1;    HELP_RHOMB_HEX_numor           =  numor  ! 75
  numor = numor + 1;    HELP_SEARCH_EXTI_numor         =  numor  ! 76
  numor = numor + 1;    HELP_SEARCH_SPGR_numor         =  numor  ! 77
  numor = numor + 1;    HELP_SET_numor                 =  numor  ! 78
  numor = numor + 1;    HELP_SETTING_numor             =  numor  ! 79
  numor = numor + 1;    HELP_SFAC_numor                =  numor  ! 80
  numor = numor + 1;    HELP_SFHKL_numor               =  numor  ! 81
  numor = numor + 1;    HELP_SG_numor                  =  numor  ! 82
  numor = numor + 1;    HELP_SG_ALL_numor              =  numor  ! 83
  numor = numor + 1;    HELP_SG_EXTI_numor             =  numor  ! 84
  numor = numor + 1;    HELP_SG_INFO_numor             =  numor  ! 85
  numor = numor + 1;    HELP_SG_SUB_numor              =  numor  ! 86 
  numor = numor + 1;    HELP_SHANNON_numor             =  numor  ! 87
  numor = numor + 1;    HELP_SHELL_numor               =  numor  ! 88
  numor = numor + 1;    HELP_SHIFT_2TH_numor           =  numor  ! 89
  numor = numor + 1;    HELP_SITE_INFO_numor           =  numor  ! 90
  numor = numor + 1;    HELP_SIZE_numor                =  numor  ! 91
  numor = numor + 1;    HELP_SORT_numor                =  numor  ! 92
  numor = numor + 1;    HELP_STL_numor                 =  numor  ! 93 
  numor = numor + 1;    HELP_SYMM_numor                =  numor  ! 94
  numor = numor + 1;    HELP_SYST_numor                =  numor  ! 95
  numor = numor + 1;    HELP_THERM_numor               =  numor  ! 96
  numor = numor + 1;    HELP_THETA_numor               =  numor  ! 97
  numor = numor + 1;    HELP_TITL_numor                =  numor  ! 98
  numor = numor + 1;    HELP_TRANSLATION_numor         =  numor  ! 99
  numor = numor + 1;    HELP_TRANSMISSION_numor        =  numor  !100
  numor = numor + 1;    HELP_TRICLINIC_numor           =  numor  !101
  numor = numor + 1;    HELP_TWIN_HEXA_numor           =  numor  !102
  numor = numor + 1;    HELP_TWIN_PSEUDO_HEXA_numor    =  numor  !103
  numor = numor + 1;    HELP_TWO_THETA_numor           =  numor  !104 
  numor = numor + 1;    HELP_UNIT_numor                =  numor  !105
  numor = numor + 1;    HELP_WAVE_numor                =  numor  !106
  numor = numor + 1;    HELP_WEB_numor                 =  numor  !107
  numor = numor + 1;    HELP_WRITE_BEAM_numor          =  numor  !108 
  numor = numor + 1;    HELP_WRITE_CELL_numor          =  numor  !109
  numor = numor + 1;    HELP_WRITE_CHEM_numor          =  numor  !110
  numor = numor + 1;    HELP_WRITE_DEVICE_numor        =  numor  !111 
  numor = numor + 1;    HELP_WRITE_QVEC_numor          =  numor  !112 
  numor = numor + 1;    HELP_WRITE_SG_numor            =  numor  !113
  numor = numor + 1;    HELP_WRITE_SYM_OP_numor        =  numor  !114
  numor = numor + 1;    HELP_WRITE_WAVE_numor          =  numor  !115
  numor = numor + 1;    HELP_X_wave_numor              =  numor  !116
  numor = numor + 1;    HELP_ZUNIT_numor               =  numor  !117  



 end subroutine cryscal_init
 !------------------------------------------------------

 subroutine read_cryscal_ini()
  USE cryscal_module, ONLY : DEBUG_file, INI_unit, cryscal_ini, winplotr_exe, my_editor, my_browser, my_word, &
                             WEB, AUTHOR, DEVICE,                                                           &
                             CIF_parameter, CIF_parameter_APEX, CIF_parameter_KCCD, CIF_parameter_XCALIBUR, &
                             wavelength, keyword_beam, keyword_WAVE, neutrons, X_rays,                      &
                             DENZO, EVAL, APEX, CRYSALIS,                                                   &
                             CONN_dmax_ini, CONN_dmax,                                                      &                             
                             keyword_create_ACE, keyword_create_CEL, keyword_create_CFL, keyword_create_FST, keyword_create_INS,&
                             INI_create_ACE,     INI_create_CEL,     INI_create_CFL,     INI_create_FST,     INI_create_INS,    &
                             structure_solution, structure_refinement, absorption_correction,               &
                             Create_INS_temperature, Create_INS_U_threshold
  USE HKL_module,     ONLY : n_sig, threshold
  USE macros_module,  ONLY : u_case, l_case

  implicit none
   CHARACTER (LEN=256)    :: cryscal_path_name, winplotr_path_name
   INTEGER                :: long, iostat_err, i1, i2
   LOGICAL                :: file_exist
   CHARACTER (LEN=256)    :: read_line
   CHARACTER (len=256)    :: INI_string
   real                   :: INI_real




   ! >>
   call getenv('CRYSCAL', cryscal_path_name)
   long = len_trim(cryscal_path_name)
   if(long /=0) then
    if(cryscal_path_name(long:long) == '\') then
     cryscal_path_name = cryscal_path_name(1:long-1)
    endif
    cryscal_ini= trim(cryscal_path_name) // '\cryscal.ini'
   endif
   
   call getenv('WINPLOTR', winplotr_path_name)
   long = len_trim(winplotr_path_name)
   if(long /=0) then
    if(winplotr_path_name(long:long) == '\') then
     winplotr_path_name = winplotr_path_name(1:long-1)
    endif
    winplotr_exe = trim(winplotr_path_name) // '\winplotr.exe'
    if(len_trim(cryscal_ini) == 0) cryscal_ini  = trim(winplotr_path_name)//'\cryscal.ini'
   endif
   
   if(len_trim(cryscal_ini) == 0 .and. len_trim(winplotr_exe) == 0) then
    if(DEBUG_file%write) then
     call write_DEBUG_file('', '')
     call write_DEBUG_file('No environment variable defined for CRYSCAL.', '')
     call write_DEBUG_file('', '')
    endif
    return
   endif
   
   ! <<

   !call getenv('CRYSCAL', cryscal_path_name)
   !long = len_trim(cryscal_path_name)
   !if(long /=0) then
   ! if(cryscal_path_name(long:long) == '\') then
   !  cryscal_path_name = cryscal_path_name(1:long-1)
   ! endif
   ! cryscal_ini  = trim(cryscal_path_name) // '\cryscal.ini'
   ! 
   ! call getenv('WINPLOTR', winplotr_path_name)
   ! long = LEN_TRIM(winplotr_path_name)
   ! if(long /=0) then
   !  if (winplotr_path_name(long:long) == '\') then
   !   winplotr_path_name = winplotr_path_name(1:long-1)
   !  end if
   !  winplotr_exe = trim(winplotr_path_name)//'\winplotr.exe'
   ! endif
   ! 
   !else
   ! call getenv('WINPLOTR',winplotr_path_name)
   ! long = LEN_TRIM(winplotr_path_name)
   ! if(long /=0) then
   !  if (winplotr_path_name(long:long) == '\') then
   !   winplotr_path_name = winplotr_path_name(1:long-1)
   !  end if
   !  winplotr_exe = trim(winplotr_path_name)//'\winplotr.exe'
   !  if(len_trim(cryscal_ini) == 0) cryscal_ini  = trim(winplotr_path_name)//'\cryscal.ini'
   ! else
   !  if(DEBUG_file%write) then
   !   call write_DEBUG_file('')
   !   call write_DEBUG_file('No environment variable defined for CRYSCAL.')
   !   call write_DEBUG_file('')
   !  endif
   !  return
   ! endif
   !endif

   if(DEBUG_file%write) then
    call write_DEBUG_file('', '')
    call write_DEBUG_file('CRYSCAL_ini', trim(cryscal_ini))
    call write_DEBUG_file('', '')
   endif


   INQUIRE(FILE=TRIM(cryscal_ini), EXIST=file_exist)
   IF(.NOT. file_exist) then
    if(DEBUG_file%write) then
     call write_DEBUG_file('', '')
     call write_DEBUG_file('!! CRYSCAL.ini does not exist.', '')
     call write_DEBUG_file('', '')

     call CIF_default_values     
    endif
    return
   endif


   OPEN(UNIT=INI_unit, FILE=TRIM(cryscal_ini))

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
   if(DEBUG_file%write) then
    call write_DEBUG_file('BROWSER', trim(my_browser%name))
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
   if(DEBUG_file%write) then
    call write_DEBUG_file('EDITOR', trim(my_editor%name))
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
   if(DEBUG_file%write) then
    call write_DEBUG_file('WORD', trim(my_word%name))
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
      if(DEBUG_file%write) then
       call write_DEBUG_file ('WEB_name',    trim(WEB%NAME(WEB%num_site)))
       call write_DEBUG_file ('WEB_address', trim(WEB%ADDRESS(WEB%num_site)))
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
    IF(read_line(1:7) == '[AUTHOR') then
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
      endif
     END do
    endif
   END do

   if(DEBUG_file%write) then
    call write_DEBUG_file ('AUTHOR_name',          trim(AUTHOR%name))
    call write_DEBUG_file ('AUTHOR_first_name',    trim(AUTHOR%first_name))
    call write_DEBUG_file ('AUTHOR_address',       trim(AUTHOR%address))
    call write_DEBUG_file ('AUTHOR_web',           trim(AUTHOR%web))
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

    CRYSALIS%data_collection                           = CIF_parameter%computing_data_collection
    CRYSALIS%cell_refinement                           = CIF_parameter%computing_cell_refinement
    CRYSALIS%data_reduction                            = CIF_parameter%computing_data_reduction

   endif


   if(DEBUG_file%write) then
    call write_DEBUG_file ('DEVICE_diffractometer', trim(DEVICE%diffracto))
    call write_DEBUG_file ('DEVICE_lab',            trim(DEVICE%lab))
    call write_DEBUG_file ('DEVICE_radiation',      trim(DEVICE%radiation))
    call write_DEBUG_file ('DEVICE_wave',           trim(DEVICE%wave))
    if(u_case(DEVICE%diffracto(1:4)) == 'APEX') then
     call WRITE_DEBUG_file('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call WRITE_DEBUG_file('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call WRITE_DEBUG_file('MEASUREMENT_detector',      trim(CIF_parameter%diffrn_detector))
     call WRITE_DEBUG_file('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call WRITE_DEBUG_file('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call WRITE_DEBUG_file('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_source))
     call WRITE_DEBUG_file('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call WRITE_DEBUG_file('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))

    elseif(u_case(DEVICE%diffracto(1:4)) == 'KCCD') then
     call WRITE_DEBUG_file('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call WRITE_DEBUG_file('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call WRITE_DEBUG_file('MEASUREMENT_device',        trim(CIF_parameter%diffrn_measurement_device))
     call WRITE_DEBUG_file('DIFFRACTION_source',        trim(CIF_parameter%diffrn_source))
     call WRITE_DEBUG_file('DIFFRACTION_detector',      trim(CIF_parameter%diffrn_detector))
     call WRITE_DEBUG_file('DIFFRACTION_detector_area', trim(CIF_parameter%diffrn_detector_area_resol_mean))
     call WRITE_DEBUG_file('RADIATION_wavelength',      trim(CIF_parameter%diffrn_radiation_wavelength))
     call WRITE_DEBUG_file('RADIATION_type',            trim(CIF_parameter%diffrn_radiation_type))
     call WRITE_DEBUG_file('RADIATION_source',          trim(CIF_parameter%diffrn_radiation_type))
     call WRITE_DEBUG_file('RADIATION_monochromator',   trim(CIF_parameter%diffrn_radiation_monochromator))
     call WRITE_DEBUG_file('RADIATION_probe',           trim(CIF_parameter%diffrn_radiation_probe))

    elseif(u_case(DEVICE%diffracto(1:7)) == 'SAPHIRE' .or. u_case(DEVICE%diffracto(1:9)) == 'EXCALIBUR') then
     call WRITE_DEBUG_file('MEASUREMENT_device_type',   trim(CIF_parameter%diffrn_measurement_device_type))
     call WRITE_DEBUG_file('MEASUREMENT_method',        trim(CIF_parameter%diffrn_measurement_method))
     call WRITE_DEBUG_file('DIFFRACTION_detector',      trim(CIF_parameter%diffrn_detector))
     call WRITE_DEBUG_file('DIFFRACTION_detector_area', trim(CIF_parameter%diffrn_detector_area_resol_mean))
    endif   
    call WRITE_DEBUG_file('COMPUTING_data_collection',  trim(CIF_parameter%computing_data_collection))
    call WRITE_DEBUG_file('COMPUTING_cell_refinement',  trim(CIF_parameter%computing_cell_refinement))
    call WRITE_DEBUG_file('COMPUTING_data_reduction',   trim(CIF_parameter%computing_data_reduction))

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
   if(DEBUG_file%write) then
    write(read_line, '(F5.2,a)') n_sig,          '    ! sigma_I (used in the SEARCH_GROUP routine)'
    call write_DEBUG_file ('SEARCH_I_sig',      trim(read_line))
    write(read_line, '(F6.3,a)') threshold,      '   ! I > threshold * I_max (used in the SEARCH_GROUP routine)'
    call write_DEBUG_file ('SEARCH_threshold',  trim(read_line))
    write(read_line, '(F6.2,a)') CONN_dmax_ini,  '    ! for connectivity calculations '
    call write_DEBUG_file ('CONNECTIVITY_dmax', trim(read_line))
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
      if(l_case(read_line(1:11)) == 'temperature') then
       long = len_trim(INI_string)
       if(u_case(INI_string(long:long)) == 'K') then
        read(INI_string(1:long-1), *, iostat=iostat_err) INI_real
        if(iostat_err /=0) exit        
        Create_INS_temperature = INI_real - 273
       else
        read(INI_string, *, iostat = iostat_err) INI_real
        if(iostat_err /=0) exit        
        Create_INS_temperature = INI_real
       endif
       cycle
      elseif(l_case(read_line(1:11)) == 'u_threshold') then 
       read(ini_string, *, iostat=iostat_err) INI_real
       if(iostat_err/=0) exit
       Create_INS_U_threshold = INI_real
       cycle

      end if
     end do
    endif
   end do
 
   if(DEBUG_file%write) then
    write(read_line, '(F8.2,a)') Create_INS_temperature,      '  ! experimental temperature'
    call write_DEBUG_file ('temperature',      trim(read_line))
    write(read_line, '(F6.3,a)') Create_INS_U_threshold,      '    ! atoms with Uiso > u_threshold are excluded'
    call write_DEBUG_file ('U_threshold',  trim(read_line))
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
      IF(l_case(read_line(1:10)) == 'create_ace' .and. INI_string(1:1) == '1') then
       INI_create_ACE     = .true.
       keyword_create_ACE = .true.      
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_cel' .and. INI_string(1:1) == '1') then
       INI_create_CEL     = .true.
       keyword_create_CEL = .true.      
       cycle
      ELSEIF(l_case(read_line(1:10)) == 'create_cfl' .and. INI_string(1:1) == '1') then
       INI_create_CFL     = .true.
       keyword_create_CFL = .true.      
       cycle 
      ELSEIF(L_case(read_line(1:10)) == 'create_fst' .AND. INI_string(1:1) == '1') then
       INI_create_FST     = .true.
       keyword_create_FST = .true.
       cycle 
      ELSEIF(l_case(read_line(1:10)) == 'create_ins' .and. INI_string(1:1) == '1') then
       INI_create_INS     = .true.
       keyword_create_INS = .true.      
       cycle
      endif
     end do
    endif 	
   end do
   
   if(DEBUG_file%write) then
    if (keyword_create_ACE) then
     call write_DEBUG_file ('create_ACE',      '1     ! .ACE file for Carine')
    else
     call write_DEBUG_file ('create_ACE',      '0     ! .ACE file for Carine')
    endif 
    if (keyword_create_CEL) then
     call write_DEBUG_file ('create_CEL',      '1     ! .CEL file for PowderCELL')
    else
     call write_DEBUG_file ('create_CEL',      '0     ! .CEL file for PowderCELL')
    endif 
    if (keyword_create_CFL) then
     call write_DEBUG_file ('create_CFL',      '1     ! .CFL file for CRYSCAL')
    else
     call write_DEBUG_file ('create_CFL',      '0     ! .CFL file for CRYSCAL')
    endif 
    if (keyword_create_FST) then
     call write_DEBUG_file ('create_FST',      '1     ! .FST file for FP Studio')
    else
     call write_DEBUG_file ('create_FST',      '0     ! .FST file for FP Studio')
    endif
    if (keyword_create_INS) then
     call write_DEBUG_file ('create_INS',      '1     ! .INS file for SHELXL')
    else
     call write_DEBUG_file ('create_INS',      '0     ! .INS file for SHELXL')
    endif 
    
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
    IF(read_line(1:20) == '[STRUCTURE PROGRAMS]') then
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
   
   if(DEBUG_file%write) then
    call write_DEBUG_file ('Structure_solution_name',         trim(structure_solution%name))
    call write_DEBUG_file ('Structure_solution_reference',    trim(structure_solution%reference))
    call write_DEBUG_file ('Structure_solution_cif_ref',      trim(structure_solution%cif_ref))
    call write_DEBUG_file ('Structure_refinement_name',       trim(structure_refinement%name))
    call write_DEBUG_file ('Structure_refinement_reference',  trim(structure_refinement%reference))   
    call write_DEBUG_file ('Structure_refinement_cif_ref',    trim(structure_refinement%cif_ref))   
    call write_DEBUG_file ('Absorption_correction_name',      trim(Absorption_correction%name))
    call write_DEBUG_file ('Absorption_correction_reference', trim(Absorption_correction%reference))   
    call write_DEBUG_file ('Absorption_correction_cif_ref',   trim(Absorption_correction%cif_ref))   
   endif





   CLOSE(UNIT=INI_unit)



  RETURN
 end subroutine read_cryscal_ini
!-----------------------------------------------------

 subroutine edit_a_file()
  USE cryscal_module, ONLY : my_editor, file_to_edit
  USE IO_module,      ONLY : write_info

  IF(.not. my_editor%exist) then
   call write_info('')
   call write_info('  No editor defined in the CRYSCAL.ini setting file .')
   call write_info('')
   return
  endif

  call system(TRIM(my_editor%name)//' '//TRIM(file_to_edit))

 end subroutine  edit_a_file

!-----------------------------------------------------
subroutine output_setting
 USE cryscal_module, ONLY :  INI_unit, cryscal_ini, my_browser, my_editor, WEB, DEVICE, AUTHOR, message_text,  &
                             CONN_dmax_ini, INI_create_ACE, INI_create_CEL, INI_create_CFL, INI_create_FST, INI_create_INS,    &
                             structure_solution, structure_refinement, absorption_correction,      &
                             Create_INS_temperature, Create_INS_U_threshold
 USE HKl_module,     ONLY :  n_sig, threshold
 USE IO_module,      ONLY :  write_info
 implicit none
  integer                 :: i

 IF(LEN_TRIM(cryscal_ini) == 0) then
   call write_info('')
   call write_info('  No CRYSCAL.ini setting file has been defined.')
   call write_info('')
   return
 endif

 call write_info('')
 call write_info(' [EXTERNAL APPLICATIONS]')
 IF(my_editor%exist)    then
  call write_info('  > EDITOR  : '//TRIM(my_editor%name))
 else
  call write_info('  > EDITOR  :   not defined.')
 endif
 IF(my_browser%exist)   then
  call write_info('  > BROWSER : '//TRIM(my_browser%name))
 else
  call write_info('  > BROWSER :   not defined.')
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
 call write_info('')

 call write_info(' [AUTHOR]')
 call write_info('   > name           = '//trim(AUTHOR%name))
 call write_info('   > first_name     = '//trim(AUTHOR%first_name))
 call write_info('   > address        = '//trim(AUTHOR%ADDRESS))
 call write_info('   > email          = '//trim(AUTHOR%email))
 call write_info('   > web            = '//trim(AUTHOR%WEB))
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
 write(message_text, '(F6.2)') Create_INS_temperature
 call write_info('   > Temperature    = '//   trim(message_text))
 write(message_text, '(F4.1)') Create_INS_U_threshold
 call write_info('   > U_threshold    ='//  trim(message_text))
 call write_info('')
 
 
 call write_info(' [COMMAND LINE ARGUMENT]')
 if(INI_create_ACE) then
  call write_info('   > create_ACE    = 1')
 else
  call write_info('   > create_ACE    = 0')  
 endif
 if(INI_create_CEL) then
  call write_info('   > create_CEL    = 1')
 else
  call write_info('   > create_CEL    = 0')  
 endif
 if(INI_create_CFL) then
  call write_info('   > create_CFL    = 1')
 else
  call write_info('   > create_CFL    = 0')  
 endif
 if(INI_create_FST) then
  call write_info('   > create_FST    = 1')
 else
  call write_info('   > create_FST    = 0')  
 endif
 
 if(INI_create_INS) then
  call write_info('   > create_INS    = 1')
 else
  call write_info('   > create_INS    = 0')  
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
 
  


 
 

end subroutine output_setting

!--------------------------------------------------------------------------------


subroutine get_X_radiation(input_string)
 use wavelength_module
 use cryscal_module, only : wavelength, keyword_WAVE, neutrons, X_rays 

 implicit none
  character (len=*), intent(in)   :: input_string


     IF(input_string(1:4) == 'X_AG' .or. input_string(1:3) == 'XAG') then
      wavelength =  X_target(1)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(1)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(input_string(1:4) == 'X_MO' .or. input_string(1:3) == 'XMO') then
      wavelength =  X_target(2)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(2)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(input_string(1:4) == 'X_CU' .or. input_string(1:3) == 'XCU') then
      wavelength =  X_target(3)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(3)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(input_string(1:4) == 'X_NI' .OR. input_string(1:3) == 'XNI') then
      wavelength = X_target(4)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(4)%logic   = .true.
      anti_cathode = .true.
     ELSEIF(input_string(1:4) == 'X_CO' .or. input_string(1:3) == 'XCO') then
      wavelength =  X_target(5)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(5)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(input_string(1:4) == 'X_FE' .or. input_string(1:3) == 'XFE') then
      wavelength =  X_target(6)%wave(1)
      X_target(1:tabulated_target_nb)%logic = .false.
      X_target(6)%logic   = .true.
      keyword_WAVE = .true.
      anti_cathode = .true.
     ELSEIF(input_string(1:4) == 'X_CR' .or. input_string(1:3) == 'XCR') then
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
subroutine write_DEBUG_file(field, value)
 use cryscal_module, only : DEBUG_file
 implicit none
  character (len=*), intent(in)           :: field
  character (len=*), intent(in), optional :: value
  integer                                 :: long
  character (len=24)                      :: fmt_

  long = len_trim(field)

  if(.not. present(value)) then
   write(DEBUG_file%unit, '(a)') trim(field)

  else
   !if(value(1:1) == '?' .or. len_trim(value)==0) then
   ! write(DEBUG_file%unit, '(a)') ''
   !else
    write(fmt_, '(a,i2,a)') '(1x,a,', 32-long, 'x, 2a)'
    write(DEBUG_file%unit, trim(fmt_))  field,  ': ', adjustl(value)
   !endif
  endif

end subroutine write_DEBUG_file


!-------------------------------------------------------------------------------
subroutine CIF_default_values
 use cryscal_module, only : CIF_parameter, CIF_parameter_KCCD, CIF_parameter_APEX, CIF_parameter_XCALIBUR, &
                            DEBUG_file, my_browser, my_editor, my_word, AUTHOR, DEVICE
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

 if(DEBUG_file%write) then
   call write_DEBUG_file ('BROWSER',                        trim(my_browser%name))
   call write_DEBUG_file ('EDITOR',                         trim(my_editor%name))
   call write_DEBUG_file ('WORD',                           trim(my_word%name))

   call write_DEBUG_file ('AUTHOR_name',                    trim(AUTHOR%name))
   call write_DEBUG_file ('AUTHOR_first_name',              trim(AUTHOR%first_name))
   call write_DEBUG_file ('AUTHOR_first_address',           trim(AUTHOR%address))
   call write_DEBUG_file ('AUTHOR_first_web',               trim(AUTHOR%web))

   call write_DEBUG_file ('DEVICE_diffractometer',          trim(DEVICE%diffracto))
   call write_DEBUG_file ('DEVICE_lab',                     trim(DEVICE%lab))
   call write_DEBUG_file ('DEVICE_radiation',               trim(DEVICE%radiation))
   call write_DEBUG_file ('DEVICE_wave',                    trim(DEVICE%wave))
   call WRITE_DEBUG_file ('MEASUREMENT_device_type',        trim(CIF_parameter%diffrn_measurement_device_type))
   call WRITE_DEBUG_file ('MEASUREMENT_method',             trim(CIF_parameter%diffrn_measurement_method))
   call WRITE_DEBUG_file ('MEASUREMENT_detector',           trim(CIF_parameter%diffrn_detector))
   call WRITE_DEBUG_file ('RADIATION_wavelength',           trim(CIF_parameter%diffrn_radiation_wavelength))
   call WRITE_DEBUG_file ('RADIATION_type',                 trim(CIF_parameter%diffrn_radiation_type))
   call WRITE_DEBUG_file ('RADIATION_source',               trim(CIF_parameter%diffrn_radiation_source))
   call WRITE_DEBUG_file ('RADIATION_monochromator',        trim(CIF_parameter%diffrn_radiation_monochromator))
   call WRITE_DEBUG_file ('RADIATION_probe',                trim(CIF_parameter%diffrn_radiation_probe))
   call WRITE_DEBUG_file ('COMPUTING_data_collection',      trim(CIF_parameter%computing_data_collection))
   call WRITE_DEBUG_file ('COMPUTING_cell_refinement',      trim(CIF_parameter%computing_cell_refinement))
   call WRITE_DEBUG_file ('COMPUTING_data_reduction',       trim(CIF_parameter%computing_data_reduction))
   call WRITE_DEBUG_file ('COMPUTING_structure_solution',   trim(CIF_parameter%computing_structure_solution))
   call WRITE_DEBUG_file ('COMPUTING_structure_refinement', trim(CIF_parameter%computing_structure_refinement))
   call WRITE_DEBUG_file ('COMPUTING_molecular_graphics',   trim(CIF_parameter%computing_molecular_graphics))
   !call WRITE_DEBUG_file ('COMPUTING_publication_material', trim(CIF_parameter%computing_publication_material))
   call WRITE_DEBUG_file ('COMPUTING_publication_material', trim(CIF_parameter%computing_publication_material_1))
   call WRITE_DEBUG_file ('COMPUTING_publication_material', trim(CIF_parameter%computing_publication_material_2))


 endif

 return
end subroutine CIF_default_values