! Last change:  TR   11 Dec 2007   12:19 pm
!--------------------------------------------------------------------


module cryscalc_module
 use CFML_crystallographic_symmetry, ONLY  : space_group_type
 use CFML_Atom_TypeDef,              only  : Atom_list_Type
 USE CFML_Crystal_Metrics,           ONLY  : Crystal_cell_type
 USE CFML_GlobalDeps,                ONLY  : sp, cp

 
  implicit none

  TYPE, PUBLIC :: CRYSCALC_type                 ! caracteristiques du programme CRYSCALC
   character (len=256)     :: version                ! version
   character (len=256)     :: author                 ! author
   character (len=256)     :: mail                   ! mail 
   character (len=256)     :: url                    ! cryscalc url
   character (len=256)     :: ini                    ! setting file
   character (len=256)     :: css                    ! css for HTML user's guide
   character (len=256)     :: report_css             ! css for HTML structural report
   character (len=256)     :: path_name              ! pathname 
  end type CRYSCALC_type
  type (CRYSCALC_type) :: CRYSCALC

  character (len=256), parameter               :: CFML_version = "version 5.0 (JRC, JGP)"
  character (LEN=256)                          :: winplotr_exe
  character (len=256)                          :: winplotr_path_name

  CHARACTER (LEN=512)                          :: message_text
  character (len=256)                          :: input_line
  character (len=256)                          :: CIF_string
  

  LOGICAL                                      :: input_INS               ! fichier d'entree INS
  LOGICAL                                      :: input_CFL               ! fichier d'entree CFL
  LOGICAL                                      :: input_PCR               ! fichier d'entree PCR
  LOGICAL                                      :: input_CIF               ! fichier d'entree CIF
  LOGICAL                                      :: input_CELL_P4P          ! fichier .P4P (lecture des parametres de maille)
  LOGICAL                                      :: input_CELL_m50          ! fichier .m50 Jana (lecture des parametres de maille)
  LOGICAL                                      :: input_CELL_INS          ! fichier .INS/.TES Shelx (lecture des parametres de maille)
  LOGICAL                                      :: input_CELL_x            ! fichier .x cree par DENZO
  LOGICAL                                      :: input_CELL_rmat         ! fichier .rmat cree par DIRAX
  LOGICAL                                      :: input_CELL_CIF          ! fichier .CIF
  LOGICAL                                      :: input_CELL_RED          ! ficheir .RED pour DATARED
  LOGICAL                                      :: input_CELL              ! fichier contenant les parametres de maille
  LOGICAL                                      :: step_by_step
  
  CHARACTER (LEN=80)                           :: main_title              ! Title in the INS.file
  character (len=16)                           :: keyword                 ! mot clé dans le fichier d'entree
  LOGICAL                                      :: unknown_keyword
  LOGICAL                                      :: unknown_CFL_keyword
  LOGICAL                                      :: mode_interactif         ! mode interactif
  LOGICAL                                      :: tmp_logical
  real                                         :: tmp_real

 ! real,               dimension(6)             :: cell_param              ! parametres de maille
 ! real,               dimension(6)             :: cell_param_ESD          ! ESD des parametres de maille
 ! REAL,               DIMENSION(6)             :: new_cell_param          ! parametres de maille apres transformation
 ! REAL,               DIMENSION(6)             :: new_cell_param_ESD      ! ESD des parametres de maille apres transformation
  REAL,               DIMENSION(3,3)           :: DC_ort, ort_DC
  REAL,               DIMENSION(3,3)           :: RC_ort, ort_RC
  REAL,               DIMENSION(3,3)           :: GMD
  REAL,               DIMENSION(3,3)           :: GMR

  REAL,               DIMENSION(3,3)           :: UB_matrix    ! matrice d'orientation
  LOGICAL                                      :: UB_mat_log
  LOGICAL                                      :: keyword_UB_mat

  real                                         :: wavelength              ! longueur d'onde
  real                                         :: Z_unit                  ! Z: nombre de molecules dans la maille
  integer                                      :: Z_unit_INS              ! Z: nombre de molecules dans la maille: entier dans fichier.INS
  character (len=20)                           :: space_group_symbol      ! groupe d'espace
  character (len=16)                           :: crystal_system          ! systeme cristallin

  integer, parameter                            :: nb_atom_max = 1000
  REAL(KIND=sp),      dimension(3,nb_atom_max)         :: atom_coord              ! coordonnees atomiques
  REAL,               DIMENSION(6,nb_atom_max)         :: atom_adp_aniso          ! parametres de deplacements atomiques
  REAL,               DIMENSION(nb_atom_max)           :: atom_adp_equiv          ! parametre de deplacements atomiques equivalent
  character (len=8),  dimension(nb_atom_max)           :: atom_label              ! label des atomes (ATOM)
  character (len=4),  dimension(nb_atom_max)           :: atom_typ                ! symbole des atomes (ATOM)
  logical,            dimension(nb_atom_max)           :: CONN_out                ! 
  real,               dimension(nb_atom_max)           :: atom_occ_perc           ! % occupation du site
  real,               dimension(nb_atom_max)           :: atom_occ                ! occupation du site = mult. site / mult. generale
  real,               dimension(nb_atom_max)           :: atom_Ueq                ! U
  real,               dimension(nb_atom_max)           :: atom_Biso               ! Biso
  integer,            dimension(nb_atom_max)           :: atom_mult               ! multiplicite du site
  REAL(KIND=sp),      dimension(3,nb_atom_max)         :: new_atom_coord          ! nouvelles coordonnees atomiques (apres MATR)
  REAL(KIND=sp),      DIMENSION(3)                     :: new_coord               ! nouvelles coordonnes atomiques

  integer                                              :: nb_atoms_type           ! nombre de type d'atomes de nature differente
  
  TYPE, public :: SFAC_type
   character (len=6),  dimension(nb_atom_max)          :: type               ! type des atomes (SFAC)
   real,               dimension(nb_atom_max)          :: number             ! nombre d'atomes de chaque espece (UNIT)  
  end type SFAC_type
  type (SFAC_type) :: SFAC
  !character (len=6),  dimension(nb_atom_max)           :: sfac_type               ! type des atomes (SFAC)
  !real,               dimension(nb_atom_max)           :: sfac_number             ! nombre d'atomes de chaque espece (UNIT)
  REAL,               DIMENSION(nb_atom_max)           :: sto                     ! stoechiometrie
  INTEGER,            DIMENSION(201)                   :: num_atom                ! numero atomique de l'atome
  real,               dimension(nb_atom_max)           :: nb_at                   ! nb_at= sfac_number/unit_cell_volume

  !real,               dimension(500)           :: h, k, l                 ! indices des reflections
  real,               dimension(3,500)         :: H                       ! indices des reflections entrees au clavier

  character (len=8),  dimension(500)           :: atom1_dist, atom2_dist  ! label de l'atome 1 et de l'atome 2 pour calculer la distance entre 1 et 2
  character (LEN=8),  DIMENSION(500)           :: atom1_ang,  atom2_ang,  atom3_ang, atom4_ang   ! labels des atomes 1,2,3 et 4 pour calculaer les angles
  character (len=8),  dimension(500, 100)      :: atom_bary               ! labels des atomes pour le calcul du barycentre

  TYPE, PUBLIC :: ATOM_CONN_type                         ! caracteristiques de l'atome pour le calcul de la connectivité
   character (len=8)     :: label                        ! label
   character (len=9)     :: type                         ! type
   character (len=8)     :: excluded                     ! atom to exclude from connectivity (ex: H)
  end type ATOM_CONN_type
  type (ATOM_CONN_TYPE) :: ATOM_CONN
  character (len=8)     :: CONN_species

  REAL,               dimension(3,3)           :: MAT                     ! composantes de la matrice de transformation 3*3
  REAL,               dimension(3,3)           :: MAT_inv                 ! composantes de la matrice inverse
  REAL,               DIMENSION(3,3)           :: Mit                     ! composantes de la matrice inverse transposee
  LOGICAL                                      :: Mat_integer             ! matrice composée d'entiers
  REAL                                         :: Mat_det                 ! determinant de la matrice


  CHARACTER (LEN=12), DIMENSION(3,200)         :: symm_op_STRING          ! operateurs de symétrie
  REAL,               DIMENSION(3,3,200)       :: symm_op_ROT             !  partie rotationnelle de la matrice
  REAL,               DIMENSION(3,200)         :: symm_op_TRANS           !    "    translationnelle
  integer                                      :: nb_symm_op              ! nombre d'operateurs de symetrie

  CHARACTER (LEN=256)                          :: ACE_file_name           ! nom du fichier.ACE
  CHARACTER (LEN=256)                          :: CEL_file_name           ! nom du fichier.CEL
  CHARACTER (LEN=256)                          :: CIF_file_name           ! nom du fichier.CIF
  CHARACTER (LEN=256)                          :: CIF_pymol_file_name     ! nom du fichier.CIF compatible Pymol
  CHARACTER (LEN=256), dimension(10)           :: input_CIF_file          ! fichiers .CIF utilises pour creer l'archive.CIF
  LOGICAL,             dimension(10)           :: CIF_file_exist          ! fichiers .CIF utilises pour creer l'archive.CIF
  INTEGER                                      :: nb_input_CIF_files      ! nombre de fichiers .CIF 
  
  CHARACTER (LEN=256)                          :: archive_CIF             ! nom de l'archive.CIF
  CHARACTER (LEN=256)                          :: INS_file_name           ! nom du fichier.INS
  CHARACTER (LEN=256)                          :: PCR_file_name           ! nom du fichier.PCR
  CHARACTER (LEN=256)                          :: P4P_file_name           ! nom du fichier.P4P (SAINT)
  CHARACTER (LEN=256)                          :: HKL_file_name           ! nom du fichier.HKL
  CHARACTER (LEN=256)                          :: RAW_file_name           ! nom du fichier.RAW
  CHARACTER (LEN=256)                          :: M50_file_name           ! nom du fichier.m50 (JANA)
  CHARACTER (LEN=256)                          :: X_file_name             ! nom du fichier.x (DENZO)
  CHARACTER (LEN=256)                          :: RMAT_file_name          ! nom du fichier.rmat (DIRAX)
  CHARACTER (LEN=256)                          :: ABS_file_name           ! nom du fichier.ABS (SADABS)
  CHARACTER (LEN=256)                          :: RED_file_name           ! nom du fichier.RED (DATARED)
  CHARACTER (LEN=256)                          :: TIDY_out_file_name      ! nom du fichier TIDY_out
  CHARACTER (LEN=256)                          :: FACES_file_name
  CHARACTER (LEN=256)                          :: SADABS_line_ratio                ! Ratio of minimum to maximum apparent transmission (SADABS output)
  CHARACTER (LEN=256)                          :: SADABS_line_estimated_Tmin_Tmax  ! Estimated Tmin and Tmax values (SADABS output)
  CHARACTER (LEN=256)                          :: SADABS_absorption_coeff    
  real                                         :: SADABS_ratio            ! value of the ratio Tmin/Tmax
  real                                         :: SADABS_Tmin
  real                                         :: SADABS_Tmax


  REAL                                         :: F000                    ! pouvoir diffusant dans la maille

  LOGICAL                                      :: keyword_CELL            ! connaissance des parametres de maille
  LOGICAL                                      :: keyword_NIGGLI
  logical                                      :: keyword_WAVE            ! connaissance de la longueur d'onde
  logical                                      :: keyword_BEAM            ! connaissance de la nature du faisceau incident

  character (len=16)                           :: beam_type
  LOGICAL                                      :: X_rays                  ! beam = X_rays
  !LOGICAL                                      :: X_Cu                    ! beam = X_rays Cu
  !LOGICAL                                      :: X_Mo                    ! beam = X_rays Mo
  !LOGICAL                                      :: X_Co                    ! beam = X_rays Co
  !LOGICAL                                      :: X_Cr                    ! beam = X_rays Cr
  !LOGICAL                                      :: X_Fe                    ! beam = X_rays Fe
  !LOGICAL                                      :: X_Ag                    ! beam = X_rays Ag

  LOGICAL                                      :: neutrons                ! beam = neutrons
  LOGICAL                                      :: electrons               ! bema = electrons
  logical                                      :: keyword_SFAC_UNIT       ! connaissance du contenu de la maille
  LOGICAL                                      :: keyword_CONT            ! connaissance du contenu de la maille
  LOGICAL                                      :: keyword_CHEM            ! connaissance du contenu de la maille
  logical                                      :: keyword_ZUNIT           ! connaissance du Z
  LOGICAL                                      :: keyword_MU              ! calcul du coef. d'absorption

  LOGICAL                                      :: keyword_WRITE_CELL      ! ecriture des parametres de maille
  LOGICAL                                      :: write_cell_cart
  LOGICAL                                      :: keyword_WRITE_CHEM      ! ecriture des caracteristiques de la molecule: formule, weight ...
  LOGICAL                                      :: keyword_WRITE_WAVE      ! ecriture de la longueur d'onde
  LOGICAL                                      :: keyword_WRITE_DEVICE    ! ecriture appareillage
  LOGICAL                                      :: keyword_WRITE_SG        ! ecriture du groupe d'espace
  LOGICAL                                      :: keyword_WRITE_BEAM      ! ecriture faisceau incident
  LOGICAL                                      :: keyword_WRITE_QVEC      ! ecriture vecteur de modulation
  LOGICAL                                      :: keyword_WRITE_ZUNIT     ! ecriture Zunit
  logical                                      :: keyword_SPGR            ! connaissance du groupe d'espace
  LOGICAL                                      :: get_SPGR                ! recuperation du automatique du G.E. apres SEARCH_GROUP
  LOGICAL                                      :: keyword_LSPGR           ! liste des groupes d'espace
  LOGICAL                                      :: keyword_LAUE            ! liste des classes de Laue
  CHARACTER (LEN=48), DIMENSION(14)            :: laue_class              !
  LOGICAL                                      :: keyword_ATOM_list       ! liste des atomes
  LOGICAL                                      :: write_atoms_cartesian  
  LOGICAL                                      :: write_atoms_in_A
  TYPE, PUBLIC :: cartesian_type
   CHARACTER (len=1)                            :: type
   CHARACTER (len=12)                           :: string  
  end type cartesian_type
  type (cartesian_type) :: cartesian_frame

  LOGICAL                                      :: keyword_ADP_list        ! liste des ADP des atomes
  CHARACTER (LEN=32), DIMENSION(250)           :: sg                      ! groupes d'espaces
  LOGICAL,            DIMENSION(7)             :: list_sg                 ! symmetry des groupes d'espace
  LOGICAL,            DIMENSION(7)             :: list_sg_Bravais         ! reseau de Bravais
  LOGICAL,            DIMENSION(2)             :: list_sg_centric         ! liste des groupes d'espace centric/acentric
  LOGICAL,            DIMENSION(14)            :: list_sg_laue            ! classe de Laue
  LOGICAL                                      :: list_sg_multip          ! multip. du groupe
  LOGICAL                                      :: list_sg_enantio         ! liste des groupes enantiomorphes
  LOGICAL                                      :: list_sg_chiral          ! liste des groupes chiraux 
  LOGICAL                                      :: list_sg_polar           ! liste des groupes polaires
  logical                                      :: keyword_SIZE            ! connaissance de la taille de l'echantillon
  LOGICAL                                      :: keyword_MAT             ! connaissance d'une matrice de transformation
  LOGICAL                                      :: keyword_LST_MAT         ! liste les principales matrices de transformation  
  LOGICAL                                      :: keyword_TRANSL          ! connaissance d'une translation
  REAL,               DIMENSION(4)             :: translat                ! translation
  LOGICAL                                      :: keyword_MATMUL          ! multiplication de 2 matrices 3*3
  LOGICAL                                      :: keyword_DIAG            ! diagonalisation d'une matrice 3*3

  LOGICAL                                      :: keyword_VERSION         ! version CRYSCALC
  
  LOGICAL                                      :: keyword_THERM           ! mot cle = THERM
  LOGICAL                                      :: keyword_THERM_SHELX
  LOGICAL                                      :: THERM_Uiso              !
  LOGICAL                                      :: THERM_Biso
  LOGICAL                                      :: THERM_Uij
  LOGICAL                                      :: THERM_Bij
  LOGICAL                                      :: THERM_Beta
  LOGICAL                                      :: THERM_aniso
  INTEGER                                      :: nb_therm_values         ! nombre de valeurs
  REAL,               DIMENSION(100)           :: therm_values            ! valeurs

  LOGICAL                                      :: keyword_FILE               ! existence d'un fichier.HKL (format SHELX ou import.cif)
  LOGICAL                                      :: keyword_FIND_HKL           ! recherche d'une reflection particuliere dans un fichier
  LOGICAL                                      :: keyword_find_HKL_EQUIV     ! recherche des reflections equivalentes
  LOGICAL                                      :: keyword_search_exti        ! recherche des extinctions
  LOGICAL                                      :: keyword_FIND_HKL_list      ! recherche d'une famille reflections (regle d'extinction paarticuliere)
  LOGICAL                                      :: keyword_LIST_EXTI_RULE     ! liste des regles d'extinctions possible
  logical                                      :: keyword_SEARCH_SPGR        ! recherche groupe d'espace
  LOGICAL                                      :: keyword_SYMM               ! existence d'opérateurs de symmétrie
  LOGICAL                                      :: symm_op_xyz                !  sous la forme x,y,z
  LOGICAL                                      :: symm_op_mat                !  sous la forme matricielle
  LOGICAL                                      :: keyword_find_HKL_NEG       ! recherche des reflections avec F2< 0. dans un fichier.HKL
  LOGICAL                                      :: keyword_find_HKL_POS       ! recherche des reflections avec F2> 0. dans un fichier.HKL
  LOGICAL                                      :: keyword_find_HKL_ABSENT    ! recherche des reflections qui devraient etre absentes
  LOGICAL                                      :: keyword_Rint               ! calcul du Rint
  LOGICAL                                      :: keyword_merge_HKL          !
  LOGICAL                                      :: keyword_get_Friedel_pairs_number ! determination du nombre de paires de Friedel

  LOGICAL                                      :: keyword_EDIT
  CHARACTER (LEN=256)                          :: file_to_edit
  LOGICAL                                      :: keyword_set
  LOGICAL                                      :: keyword_setting

  LOGICAL                                      :: keyword_STL             !
  LOGICAL                                      :: keyword_dhkl            !
  LOGICAL                                      :: keyword_dstar           !
  LOGICAL                                      :: keyword_theta           !
  LOGICAL                                      :: keyword_2theta          !
  LOGICAL                                      :: keyword_Qhkl            !
  LOGICAL                                      :: keyword_SFAC_hkl        !
  INTEGER                                      :: nb_STL_value
  INTEGER                                      :: nb_dhkl_value
  INTEGER                                      :: nb_dstar_value
  INTEGER                                      :: nb_Qhkl_value
  INTEGER                                      :: nb_theta_value
  INTEGER                                      :: nb_2theta_value
  REAL,               DIMENSION(100)           :: STL_value
  REAL,               DIMENSION(100)           :: dhkl_value
  REAL,               DIMENSION(100)           :: dstar_value
  REAL,               DIMENSION(100)           :: Qhkl_value
  REAL,               DIMENSION(100)           :: theta_value
  REAL,               DIMENSION(100)           :: Two_theta_value

  LOGICAL                                      :: keyword_QVEC            ! vecteur de modulation
  REAL,              DIMENSION(3)              :: qvec
  LOGICAL                                      :: keyword_BARY            ! calcul du barycentre
  LOGICAL                                      :: keyword_DIST            ! calcul de distances interatomiques
  LOGICAL                                      :: keyword_DIST_X
  LOGICAL                                      :: keyword_DIST_plus
  LOGICAL                                      :: keyword_DHA             ! calcul de la position d'H
  logical                                      :: keyword_CONN            ! connectivite autour d'un atome  
  LOGICAL                                      :: keyword_ANG             ! calcul d'angles
  LOGICAL                                      :: keyword_TRANSMISSION    ! calcul du coef. de transmission
  LOGICAL                                      :: keyword_GENHKL          ! generation de HKL
  LOGICAL                                      :: keyword_GENSAT          ! generation des satellites
  LOGICAL                                      :: keyword_OBV_REV         ! analyse maclage obverse/reverse
  INTEGER                                      :: OBV_REV_twin_matrix     ! matrice du maclage obverse/reverse

  LOGICAL                                      :: keyword_INSIDE          ! atomic coordinates in the 0.0 - 1.0 range
  logical                                      :: keyword_PAUSE

  LOGICAL                                      :: keyword_read_CEL        ! lecture fichier.CEL (PowderCELL)
  LOGICAL                                      :: keyword_read_CIF        ! lecture fichier.CIF
  LOGICAL                                      :: keyword_read_INS        ! lecture fichier.INS
  LOGICAL                                      :: keyword_read_FACES      ! lecture fichier absorb.ins
  LOGICAL                                      :: read_Q_peaks            ! lecture Q peaks dans fichier.INS/.RES
  LOGICAL                                      :: keyword_read_PCR        ! lecture fichier.PCR
  LOGICAL                                      :: keyword_read_NREPORT    ! lecture fichier nreport.HTML
  LOGICAL                                      :: keyword_read_TIDY_out   ! lecture fichier de sortie de TIDY
  LOGICAL                                      :: browse_nreport
  LOGICAL                                      :: keyword_SIR             ! creation d'un fichier pour SIR97

  LOGICAL                                      :: keyword_NEWS            ! nouveautes de CRYSCALC
  character (len=4)                            :: news_year               ! annee
  LOGICAL                                      :: keyword_SH_2th          ! decalage de 2theta
  REAL, dimension(3)                           :: shift_2theta            ! valeur du decalage en 2theta (constant, cosTheta, sinTheta)

  LOGICAL                                      :: lecture_ok
  LOGICAL                                      :: write_HKL               ! sortie des HKL
  LOGICAL                                      :: create_PAT              ! creation d'un diagramme
  LOGICAL                                      :: PM2K_out                ! creation d'un fichier reflexions pour PM2K (M. Leoni)
  LOGICAL                                      :: HKL_2THETA
  LOGICAL                                      :: HKL_THETA
  LOGICAL                                      :: HKL_STL
  LOGICAL                                      :: HKL_D
  LOGICAL                                      :: HKL_Q
  REAL                                         :: X_min, X_max


  LOGICAL                                      :: keyword_STRUCT_factor   ! calcul du facteur de structure

  LOGICAL                                      :: WRITE_SPG_info          ! liste informations sur le groupe d'espace
  LOGICAL                                      :: Write_SPG_info_all      !
  LOGICAL                                      :: WRITE_SPG_exti          ! liste les extinctions du groupe d'espace
  LOGICAL                                      :: write_SPG_all           ! liste des informations + extinctions du groupe d'espace
  LOGICAL                                      :: WRITE_SPG_subgroups     ! liste les subgroups

  LOGICAL                                      :: WRITE_SYMM_op           ! liste les operateurs de symetrie
  LOGICAL                                      :: WRITE_SHELX_SYMM_op     ! liste les op. de symetrie au format SHELX (LATT ; SYMM)
  LOGICAL                                      :: WRITE_APPLY_symm        ! applique les oper. de sym. à tous les atomes
  LOGICAL                                      :: WRITE_SITE_info         ! liste informations sur la symetrie du site des atomes
   INTEGER                                     :: nb_atom_site_info       ! nb d'atomes pour lesquels on ecrit les "site infos"
   CHARACTER (LEN=6), DIMENSION(500)           :: site_info_label         ! liste des atomes pour lesquels on ecrit les "site infos"
   LOGICAL                                     :: site_info_all_atoms     ! "site infos" pour tous les "nb_atom" atomes de la liste
  LOGICAL                                      :: WRITE_STAR_K            ! applique les parties rot. des op. sym. sur le vecteur de propagation
  LOGICAL                                      :: WRITE_PCR_SITE_info     ! liste des atomes equivalents par symetrie (format FullProf)
  LOGICAL                                      :: WRITE_PCR_mag_SITE_info ! idem pour atome magnetique
  LOGICAL                                      :: WRITE_triclinic_transf  ! matrices de transformation maille triclinique
  LOGICAL                                      :: WRITE_monoclinic_transf ! matrices de transformation systeme monoclinique
  LOGICAL                                      :: WRITE_rhomb_hex_transf  ! matrice  de transformation rhomboedrique ==> hexagonal
  LOGICAL                                      :: WRITE_hex_rhomb_transf  ! matrice  de transformation hexagonal     ==> rhomboedrique
  LOGICAL                                      :: WRITE_permutation_abc   ! matrices de transformation (systeme orthorhombique)
  LOGICAL                                      :: WRITE_twin_hexa         ! matrices de transformation dans un système hexagonal
  LOGICAL                                      :: WRITE_twin_pseudo_hexa  ! matrices de transformation dans un système monoclinique (derive d'une maille hexagonal)


  logical                                      :: keyword_HELP            ! aide en ligne
  INTEGER, PARAMETER                           :: nb_help_max = 135       ! nombre max. d'arguments de HELP
  character (len=19), dimension(nb_help_max)   :: HELP_string             ! liste des HELP
  character (len=19), dimension(nb_help_max)   :: HELP_arg                ! arguments de HELP


  integer                                      :: nb_help                 ! nombre d'arguments de HELP
  LOGICAL,            DIMENSION(nb_help_max)   :: write_keys
  LOGICAL                                      :: keyword_HEADER          ! entete de CRYSCALC
  LOGICAL                                      :: keyword_SYST            ! mot cle = SYST
  LOGICAL                                      :: keyword_RESET           ! mot cle = RESET/RAZ/INIT
  CHARACTER (LEN=256)                          :: SYST_command            !
  logical                                      :: keyword_KEY             ! liste des mots cles
  logical                                      :: keyword_CLA             ! liste des arguments en ligne de commande

  REAL, parameter                              :: pi=3.1415926535897932

  LOGICAL                                      :: keyword_X_WAVE          ! write X-rays Kalpha1, Kalpha2 wavelength
  REAL                                         :: LOCK_wave_value         ! comparaison with targets wavelength values
  LOGICAL                                      :: CIF_format80            ! 
  REAL                                         :: CIF_torsion_limit       ! 
  LOGICAL                                      :: include_RES_file        !
  CHARACTER (LEN=256)                          :: RES_file
  LOGICAL                                      :: include_HKL_file        !
  
  LOGICAL                                      :: include_SQUEEZE
  LOGICAL                                      :: update_parameters       ! update cell parameters, atomic coordinates ... after a matrix transdformation
  LOGICAL                                      :: report_header           ! write header in HTML report
  LOGICAL                                      :: expert_mode             ! allows specific instructions (ex: fic  = FILE import.cif)
  LOGICAL                                      :: skip_start_menu         ! skip start menu 
  LOGICAL                                      :: hkl_statistics          ! output statistics on hkl reflections
  LOGICAL                                      :: hkl_format_free         ! format libre pour un fichier .hkl
  LOGICAL                                      :: hkl_format_SHELX        ! format SHELX pour un fichier .hkl  
  CHARACTER (LEN=32)                           :: hkl_format              ! format de lecture d'un fichier .hkl (si non libre)
  CHARACTER (LEN=32)                           :: hkl_SHELX_format        ! format SHEXL de lecture d'un fichier .hkl (3I4,2F8.2)
  LOGICAL                                      :: keep_bond_str_out       ! keep bond_str.out file
  
  
  TYPE, PUBLIC :: pdp_simulation_type                 
   LOGICAL                                      :: cu                  ! use of Ka cu radiation for powder diffration pattern simulation
   character (len=16)                           :: beam
   REAL                                         :: wave  
  end type pdp_simulation_type
  type (pdp_simulation_type) :: pdp_simu

  
  integer                                      :: nb_atom                 ! nombre d'atomes dans la liste
  integer                                      :: nb_hkl                  ! nombre de reflections
  integer                                      :: nb_hkl_SFAC_calc        ! nombre de reflections pour le calcul de facteur de structure
  integer                                      :: nb_dist_calc            ! nombre de distances a calculer  
  real(kind=cp)                                :: CONN_dmax_ini           ! valeur initiale de CONN_dmax
  real(kind=cp)                                :: CONN_dmax               ! dist. max. pour calcul de la connectivité
  real(kind=cp)                                :: CONN_dmin               ! dist. min. pour calcul de la connectivité
  real(kind=cp)                                :: CONN_dmin_ini           ! valeur initiale de dist. min. pour calcul de la connectivité
  logical                                      :: CONN_all                ! calcul de connectivite pour tous les atomes
  logical                                      :: CONN_all_X              ! calcul de connectivite pour tous les atomes d'une meme espece
  logical                                      :: CONN_self               ! calcul de la connectivite par un atome identique
  LOGICAL                                      :: CONN_ang                ! calcul des angles
  LOGICAL                                      :: CONN_out_condensed      ! short output
  LOGICAL                                      :: CONN_excluded
  logical                                      :: calcul_BVS              ! calcul des BVS
  REAL                                         :: dist_coef
  REAL                                         :: DIST_plus
  REAL                                         :: DIST_AH
  integer                                      :: nb_ang_calc             ! nombre d'angles a calculer
  integer                                      :: nb_bary_calc            ! nombre de barycentres a calculer
  integer, dimension(100)                      :: nb_atom_bary            ! nombre d'atomes a considerer dans le calcul du barycentre

  INTEGER                                      :: nb_dim_transmission     ! nombre de dimensions pour le calcul de la transmission
  REAL,              DIMENSION(100)            :: dim_transmission        ! dimensions pour le calcul de la transmission


  INTEGER                                      :: nb_sort                 ! nombre de SORT
  CHARACTER (LEN=6), DIMENSION(100)            :: sort_type               ! type de tri (d, stl, I)
  LOGICAL,           DIMENSION(100)            :: sort_plot               ! trace du fichier apres SORT
  LOGICAL,           DIMENSION(100)            :: sort_out                !
  INTEGER,           DIMENSION(100)            :: sort_out_n
  LOGICAL                                      :: file_out
  INTEGER                                      :: file_out_n
  
  INTEGER                                      :: nb_shell                ! nombre de SHELL
  CHARACTER (LEN=6), DIMENSION(100)            :: shell_type              ! type de SHELL (d, stl, I)
  LOGICAL,           DIMENSION(100)            :: shell_plot              ! trace du fichier apres SHELL
  REAL,              DIMENSION(100)            :: shell_min, shell_max    !
  LOGICAL                                      :: shell_arg_min_max

  logical                                      :: cut_off                 ! coupure
  REAL                                         :: cut_off_min             ! valeur min. de la coupure
  REAL                                         :: cut_off_max             ! valeur max. de la coupure
  REAL                                         :: expected_cut_off_min    ! valeur min. attendue de la coupure
  REAL                                         :: expected_cut_off_max    ! valeur max. attendue de la coupure

  CHARACTER (len=256)                          :: sample_job
  LOGICAL                                      :: get_sample_ID

  TYPE, public :: CREATE_INS_type
   REAL                                         :: temperature
   REAL                                         :: U_threshold
   logical                                      :: ANIS
  end type CREATE_INS_type
  type (CREATE_INS_type) :: CREATE_INS

  !REAL                                         :: Create_INS%temperature, Create_INS%U_threshold

  real                                         :: max_res, expected_max_res
  real                                         :: min_res, expected_min_res
  REAL                                         :: max_d_hkl, expected_max_d_hkl
  REAL                                         :: min_d_hkl, expected_min_d_hkl
  REAL                                         :: max_theta, expected_max_theta
  REAL                                         :: min_theta, expected_min_theta

  real                                         :: intensity_min, intensity_max
  real                                         :: expected_intensity_min, expected_intensity_max
  real                                         :: sig_coef

  LOGICAL                                      :: keyword_DIR_ANG
  LOGICAL                                      :: keyword_REC_ANG
  REAL,  DIMENSION(3,50)                       :: U1_DA  ! vecteur U1 de l'espace direct
  REAL,  DIMENSION(3,50)                       :: U2_DA  ! vecteur U2 de l'espace direct
  REAL,  DIMENSION(3,50)                       :: U1_RA  ! vecteur U1 de l'espace reciproque
  REAL,  DIMENSION(3,50)                       :: U2_RA  ! vecteur U2 de l'espace reciproque
  INTEGER                                      :: nb_da  ! nombre de calcul d'angle dans l'espace direct
  INTEGER                                      :: nb_ra  ! nombre de calcul d'angle dans l'espace reciproque

  LOGICAL                                      :: keyword_create_REPORT
  LOGICAL                                      :: long_report
  LOGICAL                                      :: HTML_report
  LOGICAL                                      :: text_report
  LOGICAL                                      :: latex_report
  LOGICAL                                      :: keyword_create_CIF
  LOGICAL                                      :: keyword_create_CRYSCALC_HTML
  LOGICAL                                      :: keyword_create_CRYSCALC_NEWS
  LOGICAL                                      :: browse_cryscalc_HTML
  
  LOGICAL                                      :: keyword_WRITE_REF_APEX
  LOGICAL                                      :: keyword_WRITE_REF_DENZO
  LOGICAL                                      :: keyword_WRITE_REF_EVAL  
  LOGICAL                                      :: keyword_WRITE_REF_KCCD
  LOGICAL                                      :: keyword_WRITE_REF_SADABS
  LOGICAL                                      :: keyword_WRITE_REF_SUPERNOVA
  LOGICAL                                      :: keyword_WRITE_REF_ABS_CRYSALIS
  LOGICAL                                      :: keyword_WRITE_REF_X2S
  LOGICAL                                      :: keyword_WRITE_REF_XCALIBUR
  
  LOGICAL                                      :: keyword_modif_ARCHIVE
  LOGICAL                                      :: keyword_create_ARCHIVE 
  LOGICAL                                      :: keyword_SOLVE_to_INS

  LOGICAL                                      :: keyword_create_ACE    !  creation d'un fichier.ACE pour Carine
  LOGICAL                                      :: keyword_create_CEL    !  creation d'un fichier.CEL pour powdercell
  LOGICAL                                      :: keyword_create_INS    !  creation d'un fichier.INS pour SHELXL
  LOGICAL                                      :: keyword_create_CFL    !  creation d'un fichier.CFL pour CRYSCALC
  LOGICAL                                      :: keyword_create_FST    !  creation d'un fichier.FST pour FP Studio
  LOGICAL                                      :: keyword_create_PRF    !  creation d'un fichier.PRF pour FullProf
  LOGICAL                                      :: create_FST_POLY       !
  LOGICAL                                      :: create_FST_MOLE       !
  LOGICAL                                      :: FST_no_H              !
  LOGICAL                                      :: launch_FP_Studio      !  
  
  LOGICAL                                      :: keyword_create_TIDY   ! creation d'un fichier TIDY.dat pour STIDY  
  LOGICAL                                      :: keyword_create_SOLVE  ! creation de fichiers d'entree pour SHELXS, SIR, Superflip  
  LOGICAL                                      :: create_CIF_PYMOL      ! creation d'un fichier.CIF compatible Pymol
  LOGICAL                                      :: create_SHAPE_file     ! creation d'un fichier.SHP pour SHAPE
  LOGICAL                                      :: poly_vol_calc         ! calcul du volume d'un polyedre de coordination
  LOGICAL                                      :: INI_create_ACE
  LOGICAL                                      :: INI_create_CEL
  LOGICAL                                      :: INI_create_CFL
  LOGICAL                                      :: INI_create_FST
  LOGICAL                                      :: INI_create_CIF_PYMOL
  LOGICAL                                      :: INI_create_INS
  LOGICAL                                      :: INI_create_PRF

  LOGICAL                                      :: keyword_WEB
  CHARACTER (LEN=256)                          :: URL_address

  INTEGER, parameter                           :: Input_unit          = 11
  INTEGER, parameter                           :: HKL_unit            = 12
  INTEGER, parameter                           :: HKLF5_unit          = 13
  INTEGER, parameter                           :: HKL_list_out1_unit  = 31
  INTEGER, parameter                           :: HKL_list_out2_unit  = 32

  INTEGER, parameter          :: HELP_unit           = 38    ! unité logique attribue au fichier cryscalc_manual.txt
  INTEGER, parameter          :: KEYS_unit           = 39    ! unité logique attribue au fichier cryscalc_keys.txt
  INTEGER, parameter          :: INI_unit            = 40
  INTEGER, parameter          :: CIF_unit            = 41    ! unite logique attribuée au fichier CRYSCALC.CIF
  INTEGER, parameter          :: CIF_archive_unit    = 42
  INTEGER, parameter          :: IMPORT_CIF_unit     = 43
  INTEGER, parameter          :: CFL_unit            = 44    ! unite logique attribuee au fichier CRYSCALC.CFL
  INTEGER, parameter          :: INS_unit            = 45    ! unite logique attribuee au fichier CRYSCALC_new.INS
  INTEGER, parameter          :: HTML_unit           = 46    ! unite logique attribuee au fichier .HTML
  INTEGER, parameter          :: TEXT_unit           = 47    ! unite logique attribuee au fichier .TXT
  INTEGER, parameter          :: LATEX_unit          = 48    ! unite logique attribuee au fichier .LTX (LATEX)
  INTEGER, parameter          :: CEL_unit            = 49    ! unite logique attribuee au fichier .CEL
  INTEGER, parameter          :: ACE_unit            = 50    ! unite logique attribuee au fichier .ACE
  INTEGER, parameter          :: NEWS_unit           = 51    ! unité logique attribuee au fichier cryscalc_news.tx      

  INTEGER, parameter          :: CFL_read_unit       = 52
  INTEGER, parameter          :: CEL_read_unit       = 53      ! unite logique attribuee au fichier .CEL (en lecture)
  INTEGER, parameter          :: CIF_read_unit       = 54      ! unite logique attribuee au fichier .CIF (en lecture)
  INTEGER, parameter          :: CIF_pymol_unit      = 55      ! unite logique attribuee au fichier .CIF (en lecture)
  INTEGER, parameter          :: INS_read_unit       = 56      ! unite logique attribuee au fichier .INS (en lecture)
  INTEGER, parameter          :: PCR_read_unit       = 57      ! unite logique attribuee au fichier .PCR (en lecture)
  INTEGER, parameter          :: P4P_read_unit       = 58      ! unite logique attribuee au fichier .P4P (en lecture)
  INTEGER, parameter          :: M50_read_unit       = 59      ! unite logique attribuee au fichier .m50 de Jana (en lecture)
  INTEGER, parameter          :: X_read_unit         = 60      ! unite logique attribuee au fichier.x (DENZO)
  INTEGER, parameter          :: RMAT_read_unit      = 61      ! unite logique attribuee au fichier.RMAT (DIRAX)
  INTEGER, parameter          :: ABS_read_unit       = 62      ! unite logique attribuee au fichier.ABS (SADABS)
  INTEGER, parameter          :: RED_read_unit       = 63      ! unite logique attribuee au fichier.RED (DATARED)
  INTEGER, parameter          :: TIDY_read_unit      = 64      ! unite logique attribuee au fichier TIDY_out 
  
  INTEGER, parameter          :: PAT_unit            = 65      ! unite logique attribuee au fichier CRYSCALC_pat.xy
  INTEGER, parameter          :: PRF_unit            = 66      ! unite logique attribuee au fichier CRYSCALC_pat.PRF
  INTEGER, parameter          :: PM2K_unit           = 67      ! unite logique attribuee au fichier CRYSCALC_PM2K.inp



  INTEGER, parameter                           :: tmp_unit        = 22

  
  LOGICAL                                      :: keyword_P4P               ! lecture fichier.P4P: extraction cell_param, wave, ...
  logical                                      :: keyword_RAW               ! lecture fichier.RAW + creation fichier.HKL format SHELX
  logical                                      :: keyword_HKL               !

  !real                                         :: unit_cell_volume        ! volume de la maille
  !REAL                                         :: unit_cell_volume_esd
  !CHARACTER (LEN=64)                           :: unit_cell_Bravais

  REAL(KIND=sp)                                :: SP_value                ! produit scalaire

  TYPE (space_group_type)                      :: SPG
  TYPE (crystal_cell_type)                     :: crystal_cell
  TYPE (Atom_list_Type)                        :: Atm_list


  LOGICAL                                      :: keyword_mendel
  LOGICAL                                      :: mendel_plot      ! trace de f(s)
  CHARACTER (LEN=4) , DIMENSION(50)            :: mendel_atom_label
  INTEGER                                      :: mendel_atom_nb

  LOGICAL                                      :: keyword_shannon
  CHARACTER(LEN=4)                             :: shannon_atom_label

  LOGICAL                                      :: keyword_mag
  CHARACTER(LEN=4)                             :: mag_atom_label

  LOGICAL                                      :: keyword_DATA_neutrons
  logical                                      :: DATA_neutrons_PLOT

  LOGICAL                                      :: keyword_DATA_Xrays
  logical                                      :: data_Xrays_PLOT

  logical                                      :: keyword_DATA_atomic_density
  logical                                      :: data_atomic_density_PLOT
  logical                                      :: keyword_DATA_atomic_radius
  logical                                      :: data_atomic_radius_PLOT
  logical                                      :: keyword_DATA_atomic_weight
  logical                                      :: data_atomic_weight_PLOT

  logical                                      :: known_atomic_label
  logical                                      :: known_atomic_features
  logical                                      :: known_data_neutrons
  logical                                      :: known_data_X
  LOGICAL                                      :: known_space_groups
  LOGICAL                                      :: known_theta
  LOGICAL                                      :: known_cell_esd
  LOGICAL                                      :: known_shannon_lines
  LOGICAL                                      :: known_mag_lines


  INTEGER         :: HELP_ABSENT_HKL_numor
  INTEGER         :: HELP_ABSORPTION_numor
  INTEGER         :: HELP_ACTA_numor
  INTEGER         :: HELP_ANG_numor
  INTEGER         :: HELP_APPLY_OP_numor
  INTEGER         :: HELP_ATOM_numor
  INTEGER         :: HELP_ATOM_LIST_numor
  INTEGER         :: HELP_BARY_numor
  INTEGER         :: HELP_BEAM_numor
  INTEGER         :: HELP_CELL_numor
  INTEGER         :: HELP_CHEM_numor
  INTEGER         :: HELP_CONN_numor
  INTEGER         :: HELP_CONT_numor
  INTEGER         :: HELP_CREATE_ACE_numor
  INTEGER         :: HELP_CREATE_CEL_numor
  INTEGER         :: HELP_CREATE_CFL_numor
  INTEGER         :: HELP_CREATE_FST_numor
  INTEGER         :: HELP_CREATE_INS_numor
  INTEGER         :: HELP_CREATE_REPORT_numor
  INTEGER         :: HELP_CREATE_SOLVE_numor
  INTEGER         :: HELP_CREATE_TIDY_numor
  INTEGER         :: HELP_D_HKL_numor
  INTEGER         :: HELP_D_STAR_numor
  INTEGER         :: HELP_DATA_ATOMIC_DENSITY_numor
  INTEGER         :: HELP_DATA_ATOMIC_RADIUS_numor
  INTEGER         :: HELP_DATA_ATOMIC_WEIGHT_numor
  INTEGER         :: HELP_DATA_NEUTRONS_numor
  INTEGER         :: HELP_DATA_XRAYS_numor
  INTEGER         :: HELP_DIAG_MAT_numor
  INTEGER         :: HELP_DIR_ANG_numor
  INTEGER         :: HELP_DIST_numor
  INTEGER         :: HELP_DIST_DHA_numor
  INTEGER         :: HELP_EDIT_numor
  INTEGER         :: HELP_EQUIV_numor
  INTEGER         :: HELP_EXIT_numor
  INTEGER         :: HELP_FILE_numor
  INTEGER         :: HELP_FIND_HKL_numor
  INTEGER         :: HELP_FIND_HKL_LIST_numor
  INTEGER         :: HELP_FRIEDEL_PAIRS_numor
  INTEGER         :: HELP_GEN_HKL_numor
  INTEGER         :: HELP_HEADER_numor
  INTEGER         :: HELP_HKL_numor
  INTEGER         :: HELP_HKL_NEG_numor
  INTEGER         :: HELP_HKL_POS_numor
  INTEGER         :: HELP_HEX_RHOMB_numor
  INTEGER         :: HELP_INSIDE_numor
  INTEGER         :: HELP_LIST_EXTI_numor
  !INTEGER         :: HELP_WRITE_HKL_numor  << remplace par HELP_FIND_HKL_LIST_numor
  INTEGER         :: HELP_LIST_KEYS_numor
  INTEGER         :: HELP_LIST_LAUE_numor
  INTEGER         :: HELP_LIST_MATR_numor
  INTEGER         :: HELP_LIST_SG_numor
  INTEGER         :: HELP_WRITE_SYM_OP_numor
  INTEGER         :: HELP_MAG_numor
  INTEGER         :: HELP_MAN_numor
  INTEGER         :: HELP_MAN_HTML_numor
  INTEGER         :: HELP_MATMUL_numor
  INTEGER         :: HELP_MATR_numor
  INTEGER         :: HELP_MENDEL_numor
  INTEGER         :: HELP_MERGE_numor
  INTEGER         :: HELP_MONOCLINIC_numor
  INTEGER         :: HELP_NEWS_numor
  INTEGER         :: HELP_NIGGLI_CELL_numor
  INTEGER         :: HELP_OBV_REV_numor
  INTEGER         :: HELP_P4P_numor
  INTEGER         :: HELP_PAUSE_numor
  INTEGER         :: HELP_PERMUT_numor
  INTEGER         :: HELP_Q_HKL_numor
  INTEGER         :: HELP_QVEC_numor
  INTEGER         :: HELP_READ_CEL_numor
  INTEGER         :: HELP_READ_CIF_numor
  INTEGER         :: HELP_READ_FACES_numor
  INTEGER         :: HELP_READ_INS_numor
  INTEGER         :: HELP_READ_NREPORT_numor
  INTEGER         :: HELP_READ_PCR_numor
  INTEGER         :: HELP_READ_TIDY_out_numor
  INTEGER         :: HELP_REC_ANG_numor
  INTEGER         :: HELP_REF_ABS_CRYSALIS_numor
  INTEGER         :: HELP_REF_APEX_numor  
  INTEGER         :: HELP_REF_DENZO_numor
  INTEGER         :: HELP_REF_EVAL_numor
  INTEGER         :: HELP_REF_KCCD_numor
  INTEGER         :: HELP_REF_SADABS_numor
  INTEGER         :: HELP_REF_SUPERNOVA_numor
  INTEGER         :: HELP_REF_X2S_numor
  INTEGER         :: HELP_REF_XCALIBUR_numor
  INTEGER         :: HELP_RESET_numor
  INTEGER         :: HELP_RINT_numor
  INTEGER         :: HELP_RHOMB_HEX_numor
  INTEGER         :: HELP_SEARCH_EXTI_numor
  INTEGER         :: HELP_SEARCH_SPGR_numor
  INTEGER         :: HELP_SET_numor
  INTEGER         :: HELP_SETTING_numor
  INTEGER         :: HELP_SFHKL_numor
  INTEGER         :: HELP_SFAC_numor
  INTEGER         :: HELP_SG_numor
  INTEGER         :: HELP_SG_ALL_numor
  INTEGER         :: HELP_SG_EXTI_numor
  INTEGER         :: HELP_SG_INFO_numor
  INTEGER         :: HELP_SG_SUB_numor
  INTEGER         :: HELP_SHANNON_numor
  INTEGER         :: HELP_SHELL_numor
  INTEGER         :: HELP_SHIFT_2TH_numor
  INTEGER         :: HELP_SITE_INFO_numor
  INTEGER         :: HELP_SIZE_numor
  INTEGER         :: HELP_SORT_numor
  INTEGER         :: HELP_STAR_K_numor
  INTEGER         :: HELP_STL_numor
  INTEGER         :: HELP_SYMM_numor
  INTEGER         :: HELP_SYST_numor
  INTEGER         :: HELP_THERM_numor
  INTEGER         :: HELP_THERM_SHELX_numor
  INTEGER         :: HELP_THETA_numor
  INTEGER         :: HELP_TITL_numor
  INTEGER         :: HELP_TRICLINIC_numor
  INTEGER         :: HELP_TRANSLATION_numor
  INTEGER         :: HELP_TRANSMISSION_numor
  INTEGER         :: HELP_TWIN_HEXA_numor
  INTEGER         :: HELP_TWIN_PSEUDO_HEXA_numor
  INTEGER         :: HELP_TWO_THETA_numor
  INTEGER         :: HELP_UB_matrix_numor
  INTEGER         :: HELP_UNIT_numor
  INTEGER         :: HELP_USER_MAT_numor
  INTEGER         :: HELP_WAVE_numor
  INTEGER         :: HELP_WEB_numor
  INTEGER         :: HELP_WRITE_ADP_numor
  INTEGER         :: HELP_WRITE_BEAM_numor
  INTEGER         :: HELP_WRITE_CELL_numor
  INTEGER         :: HELP_WRITE_CHEM_numor
  INTEGER         :: HELP_WRITE_DEVICE_numor
  INTEGER         :: HELP_WRITE_QVEC_numor
  INTEGER         :: HELP_WRITE_SG_numor
  INTEGER         :: HELP_WRITE_WAVE_numor
  INTEGER         :: HELP_WRITE_ZUNIT_numor
  INTEGER         :: HELP_X_WAVE_numor
  INTEGER         :: HELP_ZUNIT_numor

  !------------------------------------------------------------------------------------------------
  ! variables dependent of command line arguments

  logical         :: ON_SCREEN                            ! affichage du texte à l'ecran
  logical         :: ON_screen_PRF
  logical         :: CIFdep
  logical         :: acta

  !---------------- type ---------------------------------------------------------------------------------

  TYPE, PUBLIC :: debug_proc_type
   INTEGER           :: unit
   LOGICAL           :: write
   LOGICAL           :: level_1
   LOGICAL           :: level_2
   LOGICAL           :: level_3
  end type debug_proc_type
  type (debug_proc_type) :: debug_proc


  TYPE, PUBLIC :: my_appli_type
   character (len=256)           :: name
   logical                       :: exist
  END TYPE my_appli_type
  type (my_appli_type) :: my_editor
  type (my_appli_type) :: my_browser
  type (my_appli_type) :: my_word
  type (my_appli_type) :: my_pdflatex

  TYPE, PUBLIC :: unit_cell_type
   real,               dimension(6)      :: param              ! parametres de maille
   real,               dimension(6)      :: param_ESD          ! ESD des parametres de maille
   real,               dimension(6)      :: rec_param          ! parametres de la maille reciproque
   REAL,               DIMENSION(6)      :: new_param          ! parametres de maille apres transformation
   REAL,               DIMENSION(6)      :: new_param_ESD      ! ESD des parametres de maille apres transformation
   REAL,               DIMENSION(6)      :: Niggli             ! parametres de la maille Niggli
   REAL,               DIMENSION(6)      :: tmp
   real                                  :: volume             ! volume de la maille
   REAL                                  :: volume_esd
   real                                  :: rec_volume         ! volume de la maille reciproque
   CHARACTER (LEN=64)                    :: crystal_system
   CHARACTER (LEN=1)                     :: Bravais
   CHARACTER (LEN=32)                    :: lattice
   CHARACTER (LEN=32)                    :: H_M
  END TYPE unit_cell_type
  TYPE (unit_cell_type) :: unit_cell

  TYPE, PUBLIC :: crystal_type
   REAL, DIMENSION(3)   :: size
   REAL                 :: size_min
   REAL                 :: size_max
   REAL                 :: size_mid
   REAL                 :: volume
   REAL                 :: radius
   CHARACTER (LEN=32)   :: morph
   CHARACTER (LEN=32)   :: color
   CHARACTER (LEN=32)   :: density_meas
   CHARACTER (LEN=32)   :: density_diffrn
   CHARACTER (LEN=32)   :: density_method
   CHARACTER (LEN=32)   :: F000
   CHARACTER (LEN=256), dimension(100)   :: face_line
   integer,             dimension(3,100) :: face_index
   real,                dimension(100)   :: face_dim
   INTEGER              :: faces_nb
  END TYPE crystal_type
  TYPE (crystal_type) :: crystal

  TYPE, PUBLIC :: absorption_type
   REAL                 :: mu
   REAL                 :: Tmin
   REAL                 :: Tmax
  END TYPE absorption_type
  TYPE (absorption_type) :: absorption


 
  TYPE, PUBLIC :: WEB_type
   INTEGER                               :: num_site
   CHARACTER (LEN=16),  DIMENSION(10)    :: name
   CHARACTER (LEN=256), DIMENSION(10)    :: address
  END TYPE WEB_type
  TYPE (WEB_type)  :: WEB

  TYPE, PUBLIC :: AUTHOR_type
   CHARACTER (LEN=256)              :: string
   CHARACTER (LEN=256)              :: name
   CHARACTER (LEN=256)              :: first_name
   CHARACTER (LEN=256)              :: address
   CHARACTER (LEN=256)              :: email
   CHARACTER (LEN=256)              :: web
   CHARACTER (LEN=256)              :: team
   CHARACTER (LEN=256)              :: init
  END TYPE AUTHOR_type
  TYPE (AUTHOR_type) :: AUTHOR

  type, public  :: DEVICE_type
   character (len=256)              :: string
   character (len=256)              :: diffracto ! nom du diffractometre
   character (len=256)              :: lab       ! labo., service
   character (len=256)              :: radiation ! type de radiation
   character (len=256)              :: wave      ! longueur d'onde
  end type DEVICE_type
  type (DEVICE_type)  :: DEVICE

  type, public  :: PROGRAM_type
   CHARACTER (len=256)               :: name      ! programme utilise
   CHARACTER (len=256)               :: reference ! reference associée
   CHARACTER (len=256)                :: CIF_ref   ! reference associée pour fichier.CIF (max = 80 car)
  END type PROGRAM_type
  type (PROGRAM_type) :: Structure_solution
  type (PROGRAM_type) :: Structure_refinement
  type (PROGRAM_type) :: Absorption_correction


  



!----------------------------------------------------------------------------------
  TYPE, public :: Diffracto_sw_type
   CHARACTER (len=80)    :: data_collection
   CHARACTER (len=80)    :: cell_refinement
   CHARACTER (len=80)    :: data_reduction
  END TYPE Diffracto_sw_type
  TYPE (Diffracto_sw_type) :: EVAL
  TYPE (Diffracto_sw_type) :: DENZO
  TYPE (Diffracto_sw_type) :: APEX
  TYPE (Diffracto_sw_type) :: X2S
  TYPE (Diffracto_sw_type) :: SUPERNOVA
  TYPE (Diffracto_sw_type) :: XCALIBUR

!---------------------------------------------------------------------------------

 TYPE, public :: Absorption_correction_features
  CHARACTER (len=80)                :: type
  CHARACTER (len=80), dimension(6)  :: details
 END TYPE Absorption_correction_features
 TYPE (Absorption_correction_features) :: SADABS
 TYPE (Absorption_correction_features) :: ABS_CRYSALIS


!----------------------------------------------------------------------------------
 type, public :: SQUEEZE_type
  logical                             :: procedure
  integer                             :: unit
  character (len=256), dimension(500) :: read_line
  character (len=256)                 :: file
  integer                             :: nb_lines
 end type SQUEEZE_type
 type (SQUEEZE_type) :: SQUEEZE

  !----------------------------------------------------------------------------------



  type, public :: molecule_features
   CHARACTER (LEN=64)   :: common_name
   character (len=64)   :: formula
   CHARACTER (LEN=256)  :: content
   real                 :: weight
   real                 :: density
   integer              :: Z        ! nombre total d'electrons dans la molecule
   integer              :: Z_unit   ! nombre de formula units
  end type molecule_features
  type(molecule_features)  :: molecule

  integer  :: max_ref 
  type, public :: pgf_data_features
   real,               dimension(:), allocatable  :: X
   real,               dimension(:), allocatable  :: Y
   integer,            dimension(:), allocatable  :: h, k, l
   CHARACTER (LEN=32), dimension(:), allocatable  :: string
  end type pgf_data_features
  !type(pgf_data_features), dimension(max_ref) :: pgf_data
  !type (PGF_data_features), dimension(:), allocatable : PGF_data
  type (PGF_data_features) :: PGF_data

  
  type, public :: pgf_file_features
   character (len=256)  :: name
   character (len=256)  :: X_legend
   character (len=256)  :: Y_legend
   character (len=256)  :: title
  end type pgf_file_features
  type (pgf_file_features) :: pgf_file
  ! ------------ type ---------------------------------------------------------------------------------


  
  contains
  
   subroutine allocate_PGF_data_arrays(dim)
    integer, intent(in) :: dim
	integer             :: ier
    
	if(allocated(PGF_data%X))   deallocate(PGF_data%X)
	 allocate(PGF_data%X(dim), stat=ier)
     if(ier/=0) call write_alloc_error('PGF_data_X')
	
	if(allocated(PGF_data%Y))   deallocate(PGF_data%Y)
	 allocate(PGF_data%Y(dim), stat=ier)
     if(ier/=0) call write_alloc_error('PGF_data_Y')
	 
	if(allocated(PGF_data%h))   deallocate(PGF_data%h)
	 allocate(PGF_data%h(dim), stat=ier)
     if(ier/=0) call write_alloc_error('PGF_data_h')
	 
	if(allocated(PGF_data%k))   deallocate(PGF_data%k)
	 allocate(PGF_data%k(dim), stat=ier)
     if(ier/=0) call write_alloc_error('PGF_data_k')
	 
	if(allocated(PGF_data%l))   deallocate(PGF_data%l)
	allocate(PGF_data%l(dim), stat=ier)
     if(ier/=0) call write_alloc_error('PGF_data_l')
	 
	if(allocated(PGF_data%string))   deallocate(PGF_data%string)
	 allocate(PGF_data%string(dim), stat=ier)
     if(ier/=0) call write_alloc_error('PGF_data_string')
	
	return
	
   end subroutine allocate_PGF_data_arrays
   
   
   subroutine deallocate_PGF_data_arrays()
    
	if(allocated(PGF_data%X))        deallocate(PGF_data%X)
	if(allocated(PGF_data%Y))        deallocate(PGF_data%Y)
	if(allocated(PGF_data%h))        deallocate(PGF_data%h)
	if(allocated(PGF_data%k))        deallocate(PGF_data%k)
	if(allocated(PGF_data%l))        deallocate(PGF_data%l)
	if(allocated(PGF_data%string))   deallocate(PGF_data%string)
	
	return
	
   end subroutine deallocate_PGF_data_arrays
   
   
   subroutine write_alloc_error(input_string)
    character (len=*), intent(in) :: input_string
    
    write(*,*) '  ! Problem to allocate memory for '//trim(input_string)//' array. Program will be stopped!'
    !call deallocate_HKL_arrays
    !call deallocate_PGF_data_arrays
    stop
	
   end subroutine write_alloc_error
   
end module cryscalc_module

!-----------------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------------------
module hkl_module
 use cryscalc_module, only : max_ref, deallocate_PGF_data_arrays, write_alloc_error
 
 implicit none

 integer                           :: n_ref          ! nombre de reflections dans le fichier
 integer                           :: n_ref_2        ! nombre de reflections avec I > 2sig
 integer                           :: n_ref_3        ! nombre de reflections avec I > 3sig
 INTEGER                           :: n_ref_eff      ! nombre de reflections dans le fichier final
 INTEGER                           :: n_ref_neg      ! nombre de reflections d'intensite negative
 

 !type, public :: reflexion_type
  integer,             dimension(:),   allocatable   :: h, k, l, code
  real,                dimension(:),   allocatable   :: F2, sig_F2, E, E2, E3, E4, E5, E6
  real,                dimension(:),   allocatable   :: sinTheta_lambda, d_hkl, theta_hkl
  LOGICAL,             dimension(:),   allocatable   :: HKL_impair
  CHARACTER (LEN=20),  DIMENSION(:),   allocatable   :: HKL_string
  CHARACTER (LEN=256), DIMENSION(:),   allocatable   :: HKL_line
  INTEGER,             DIMENSION(:,:), allocatable   :: HKL_flag
  INTEGER,             DIMENSION(:),   allocatable   :: HKL_good
  REAL,                DIMENSION(:,:), allocatable   :: cos_dir
  REAL,                DIMENSION(:),   allocatable   :: I_sigma,I_
  LOGICAL,             DIMENSION(:),   allocatable   :: HKL_search_ok
 !end type reflexion_type
 !type (reflexion_type) :: HKL

 integer, parameter                       :: MAX_allowed = 500000
 CHARACTER (LEN=256), DIMENSION(10)       :: M95_comment_line
 INTEGER                                  :: M95_comment_line_nb


 logical                                  :: cos_exist    ! fichier contenant les cos. dir.
 logical                                  :: cif_file     ! fichier *.CIF (ex: import.cif)
 real                                     :: sig_coef     ! coef. pour les sigmas

 LOGICAL                                  :: HKL_data_known


 INTEGER                           :: n_pair, n_impair
 REAL                              :: F2_pair, F2_impair, mean_F2_pair, mean_F2_impair
 REAL                              :: F2_mean, sig_F2_mean, ratio_mean
 REAL                              :: wF2_mean
 REAL                              :: E_mean, E2_mean, E2m1_mean
 REAL                              :: E3_mean, E4_mean, E5_mean, E6_mean
 real                              :: Rint
 REAL                              :: I_ratio, ratio_criteria
 REAL                              :: n_sig, threshold
 LOGICAL                           :: ordered_HKL

 character (len=4), dimension(3)           :: requested_H_string
 INTEGER,           DIMENSION(3)           :: requested_H
 integer                                   :: requested_HKL_list
 integer, parameter                        :: HKL_rule_nb =33
 CHARACTER(LEN=32), DIMENSION(2*HKL_rule_nb) :: HKL_rule
 LOGICAL                                   :: search_H_string
 LOGICAL                                   :: search_equiv
 LOGICAL                                   :: search_friedel



 type, public :: HKL_file_features
  CHARACTER (LEN=256)                          :: name           ! nom du fichier.HKL
  CHARACTER (LEN=256)                          :: output         ! fichier de sortie .HKL
  CHARACTER (LEN=256)                          :: output2        ! fichier de sortie .HKL #2
  CHARACTER (LEN=256)                          :: HKL            ! fichier .HKL
  CHARACTER (LEN=256)                          :: d
  CHARACTER (LEN=256)                          :: stl
  CHARACTER (LEN=256)                          :: theta
  CHARACTER (LEN=256)                          :: I
  CHARACTER (LEN=256)                          :: Isig
  CHARACTER (LEN=256)                          :: merge
  CHARACTER (LEN=256)                          :: transf
  LOGICAL                                      :: CIF            ! fichier *.cif (ex: import.cif)
  LOGICAL                                      :: SHELX          ! fichier *.HKL (SHELX format)
  LOGICAL                                      :: final_y        ! fichier final.y cree par EVALCCD
  LOGICAL                                      :: RAW            ! fichier .RAW cree par SAINT
  LOGICAL                                      :: M91            ! fichier .M91 cree par JANA
  LOGICAL                                      :: M95            ! fichier .M95 cree par JANA
  LOGICAL                                      :: INT            ! fichier .INT cree par DATARED
  LOGICAL                                      :: COL            ! fichier .COL cree par COLL5 (ILL, format 4)
  LOGICAL                                      :: plot           ! trace de F2=f(sinTheta/lambda)
  LOGICAL                                      :: read_NEG       ! lecture intensites negatives
 end type HKL_file_features
 type(HKL_file_features)  :: HKL_file

 TYPE, PUBLIC :: HKL_list_features
  INTEGER             :: EXTI_number        ! numero d'extinction systematique
  LOGICAL             :: ALL                ! aucun critere sur la selection
  LOGICAL             :: OUT                ! liste a l'ecran les reflexions compatibles avec EXTI_number
  LOGICAL             :: WRITE              ! ecrit dans un fichier les reflexions compatibles avec EXTI_number
  LOGICAL             :: SUPPRESS           ! supprime du fichier initial les reflexions compatibles avec EXTI_number
 END TYPE HKL_list_features
 TYPE (HKL_list_features) :: HKL_list
 TYPE (HKL_list_features) :: HKL_list_POS
 TYPE (HKL_list_features) :: HKL_list_NEG
 TYPE (HKL_list_features) :: HKL_list_ABSENT

 
 contains

 subroutine allocate_HKL_arrays
  integer        :: ier
  
  if (allocated(h))       deallocate(h)
   allocate(h(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('h')
   
  if (allocated(k))       deallocate(k)
   allocate(k(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('k')
   
  if (allocated(l))       deallocate(l)
   allocate(l(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('l')

  if (allocated(code))    deallocate(code)
   allocate(code(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('code')
   
  if (allocated(F2))       deallocate(F2)
   allocate(F2(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('F2')
  if (allocated(sig_F2))   deallocate(sig_F2)
   allocate(sig_F2(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('sig_F2')
  if (allocated(E))       deallocate(E)
   allocate(E(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('E')
  if (allocated(E2))      deallocate(E2)
   allocate(E2(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('E2')
  if (allocated(E3))      deallocate(E3)
   allocate(E3(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('E3')
  if (allocated(E4))      deallocate(E4)
   allocate(E4(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('E4')
  if (allocated(E5))      deallocate(E5)
   allocate(E5(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('E5')
  if (allocated(E6))      deallocate(E6)
   allocate(E6(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('E6')
   
  if (allocated(sinTheta_lambda)) deallocate(sinTheta_lambda)
   allocate(sinTheta_lambda(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('sinTheta_lambda')
  
  if (allocated(d_hkl)) deallocate(d_hkl)
   allocate(d_hkl(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('d_hkl')
   
  if (allocated(Theta_hkl)) deallocate(Theta_hkl)
   allocate(Theta_hkl(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('Theta_hkl')
   
  if (allocated(HKL_impair)) deallocate(HKL_impair)
   allocate(HKL_impair(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('HKL_impair')
   
  if (allocated(HKL_string)) deallocate(HKL_string)
   allocate(HKL_string(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('HKL_string')
   
  if (allocated(HKL_line)) deallocate(HKL_line)
   allocate(HKL_line(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('HKL_line')
   
  if (allocated(HKL_flag)) deallocate(HKL_flag)
   allocate(HKL_flag(Max_ref,3), stat=ier)
   if(ier/=0) call write_alloc_error('HKL_flag')
   
  if (allocated(HKL_good)) deallocate(HKL_good)
   allocate(HKL_good(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('HKL_good')
   
  if (allocated(cos_dir)) deallocate(cos_dir)
   allocate(cos_dir(Max_ref,6), stat=ier)
   if(ier/=0) call write_alloc_error('cos_dir')
   
  if (allocated(I_sigma)) deallocate(I_sigma)
   allocate(I_sigma(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('I_sigma')
   
  if (allocated(I_)) deallocate(I_)
   allocate(I_(Max_ref), stat=ier)
   if(ier/=0)  call write_alloc_error('I_')
   
  if (allocated(HKL_search_ok)) deallocate(HKL_search_ok)
   allocate(HKL_search_ok(Max_ref), stat=ier)
   if(ier/=0) call write_alloc_error('HKL_search_ok')

  return   
 end subroutine allocate_HKL_arrays
 
 
 
 subroutine deallocate_HKL_arrays
    
  if (allocated(h))               deallocate(h)
  if (allocated(k))               deallocate(k)  
  if (allocated(l))               deallocate(l)
  if (allocated(code))            deallocate(code)
  if (allocated(F2))              deallocate(F2)
  if (allocated(sig_F2))          deallocate(sig_F2)
  if (allocated(E))               deallocate(E)
  if (allocated(E2))              deallocate(E2)
  if (allocated(E3))              deallocate(E3)
  if (allocated(E4))              deallocate(E4)
  if (allocated(E5))              deallocate(E5)
  if (allocated(E6))              deallocate(E6)
  if (allocated(sinTheta_lambda)) deallocate(sinTheta_lambda)
  if (allocated(d_hkl))           deallocate(d_hkl)
  if (allocated(Theta_hkl))       deallocate(Theta_hkl)
  if (allocated(HKL_impair))      deallocate(HKL_impair)
  if (allocated(HKL_string))      deallocate(HKL_string)
  if (allocated(HKL_line))        deallocate(HKL_line)
  if (allocated(HKL_flag))        deallocate(HKL_flag)
  if (allocated(HKL_good))        deallocate(HKL_good)
  if (allocated(cos_dir) )        deallocate(cos_dir)
  if (allocated(I_sigma))         deallocate(I_sigma)
  if (allocated(I_))              deallocate(I_)
  if (allocated(HKL_search_ok))   deallocate(HKL_search_ok)
  
  return   
  
  
 end subroutine deallocate_HKL_arrays
 
 
 
END module hkl_module


!--------------------------------------------------------------------

module pattern_profile_module
 implicit none
 
  real              :: HG, HL, H
  real              :: FWHM, eta
  real              :: particle_size 
  logical           :: size_broadening

  type, public :: profile_param
   real                     :: U
   real                     :: V
   real                     :: W
   real                     :: X
   real                     :: Y
   real                     :: Z
   real                     :: eta0
   real                     :: eta1
  end type profile_param 
  type (profile_param) :: PV, X_PV, N_PV
  type (profile_param) :: TCH, X_TCH, N_TCH
  
  type, public :: pattern_param
   real                   :: step
   real                   :: wdt
   real                   :: background
   real                   :: scale    
   real                   :: Xmin
   real                   :: Xmax
  end type pattern_param
  type (pattern_param) :: pattern
  type (pattern_param) :: X_pattern
  type (pattern_param) :: N_pattern
  

end module pattern_profile_module


!--------------------------------------------------------------------
module atome_module
   implicit none

   TYPE, public :: atomic_features
    REAL                    :: weight
    CHARACTER(LEN=4)        :: symbol
    CHARACTER(LEN=13)       :: name
	INTEGER                 :: Z    ! nombre d'electrons
   ! neutrons data
    REAL                    :: bcoh             ! longueur de diffusion coherente neutronique
    REAL                    :: SEDinc           ! section efficace de diffusion incoherente
    REAL                    :: SEA              ! section efficice d'aborption
    REAL                    :: N_SED_coh        !  section efficace de diffusion coherente (=4pi. b**2)
    REAL                    :: N_SED_inc        !  section efficace de diffusion incoherente
    REAL                    :: N_SE_absorption  !  section efficace de diffusion d'absorption (a la longueur d'onde utilisee)
  
    REAL, DIMENSION(7)      :: cam        ! coef. absorption massique pour Ag, Mo, Cu, Co, Fe, Cr
    REAL, DIMENSION(7)      :: tics       ! total interaction cross section pour Ag, Mo, Cu, Co, Fe, Cr

    REAL                    :: X_attenuation

    character (len=22)      :: config_electr    ! configuration electronique
    REAL                    :: density
    REAL                    :: radius

   END TYPE atomic_features
   TYPE (atomic_features), DIMENSION(201) :: atom


end module atome_module

! --------------------------------------------------------------------
module wavelength_module

  INTEGER, parameter                 :: tabulated_target_nb = 7

  TYPE, PUBLIC :: X_target_type
   LOGICAL             :: logic
   LOGICAL             :: write
   REAL, DIMENSION(3)  :: wave    ! longueur d'onde associée (Ka1, Ka2, Kb1)
   CHARACTER (LEN=2)   :: label   ! nature anticathode
   CHARACTER (LEN=256) :: tics_file
   CHARACTER (LEN=256) :: cam_file
  END TYPE X_target_type

  TYPE (X_target_type), DIMENSION(tabulated_target_nb) :: X_target

  logical                                              :: anti_cathode


end module wavelength_module


!----------------------------------------------------------------------

module SHELX_module

  integer                                    :: n_elem_atm    ! N. of different species
  real            ,  dimension(15)           :: n_elem        ! Number of elementos into the same species
  character(len=2),  dimension(15)           :: elem_atm      ! Character to identify the specie
  integer                                    :: n_fvar        ! Out -> N. of parameters on FVAR
  real,dimension(10)                         :: fvar          ! Out -> values of FVAR
  character(len=24)                          :: fmt_sfac, fmt_unit, fmt_fvar
  integer                                    :: n_latt        !
  CHARACTER(LEN=12)                          :: nlatt_string  ! chaine de caractere associee
  character(len=40), dimension(48)           :: car_symop     ! carte des operateurs de symetrie
  CHARACTER(LEN=12)                          :: centro_string ! centric/acentric (chaine de caractere)
  INTEGER                                    :: symm_nb       ! nombre d'operateurs de symetrie
  integer                                    :: SHELX_line_nb ! nombre de lignes non interpretees
  character(len=80), dimension(50)           :: SHELX_line    ! ligne non interpretees

end module SHELX_module

!-------------------------------------------------------------------------------------------------------

module CIF_module
  ! parametres CIF (utiles pour generer le fichier structural_REPORT.html) -----------------------

  TYPE, PUBLIC  :: CIF_parameter_type 
   CHARACTER (LEN=256)    :: sample_ID  
   CHARACTER (LEN=256)    :: formula_moiety
   CHARACTER (LEN=256)    :: formula_sum
   CHARACTER (LEN=256)    :: formula_weight
   CHARACTER (LEN=256)    :: formula_units_Z
   CHARACTER (LEN=256)    :: diffracto_device
   CHARACTER (LEN=256)    :: diffracto_radiation_type
   CHARACTER (LEN=256)    :: diffracto_radiation_source
   CHARACTER (LEN=256)    :: diffracto_radiation_wavelength
   CHARACTER (LEN=256)    :: diffracto_temperature

   CHARACTER (LEN=256)    :: diffrn_source
   CHARACTER (LEN=256)    :: diffrn_radiation_wavelength
   CHARACTER (LEN=256)    :: diffrn_radiation_type
   CHARACTER (LEN=256)    :: diffrn_radiation_source
   CHARACTER (LEN=256)    :: diffrn_radiation_monochromator
   CHARACTER (LEN=256)    :: diffrn_radiation_probe   
   CHARACTER (LEN=256)    :: diffrn_measurement_device
   CHARACTER (LEN=256)    :: diffrn_measurement_device_type
   CHARACTER (LEN=256)    :: diffrn_measurement_method
   CHARACTER (LEN=256)    :: diffrn_detector
   CHARACTER (LEN=256)    :: diffrn_detector_area_resol_mean
   CHARACTER (LEN=256)    :: diffrn_theta_full
   CHARACTER (LEN=256)    :: diffrn_theta_max

   CHARACTER (LEN=256)    :: computing_data_collection
   CHARACTER (LEN=256)    :: computing_cell_refinement
   CHARACTER (LEN=256)    :: computing_data_reduction
   CHARACTER (LEN=256)    :: computing_structure_solution
   CHARACTER (LEN=256)    :: computing_structure_refinement
   CHARACTER (LEN=256)    :: computing_molecular_graphics
   !CHARACTER (LEN=80)    :: computing_publication_material
   CHARACTER (LEN=256)    :: computing_publication_material_1
   CHARACTER (LEN=256)    :: computing_publication_material_2

   CHARACTER (LEN=256)    :: symmetry_cell_setting
   CHARACTER (LEN=256)    :: symmetry_space_group
   CHARACTER (LEN=256)    :: symmetry_IT_number
   CHARACTER (LEN=256)    :: cell_length_a
   CHARACTER (LEN=256)    :: cell_length_b
   CHARACTER (LEN=256)    :: cell_length_c
   CHARACTER (LEN=256)    :: cell_angle_alpha
   CHARACTER (LEN=256)    :: cell_angle_beta
   CHARACTER (LEN=256)    :: cell_angle_gamma
   CHARACTER (LEN=256)    :: cell_volume
   CHARACTER (LEN=256)    :: cell_formula_units_Z
   CHARACTER (LEN=256)    :: exptl_density
   CHARACTER (LEN=256)    :: exptl_mu
   CHARACTER (LEN=256)    :: diffrn_reflns_number
   CHARACTER (LEN=256)    :: diffrn_reflns_av_R_equivalents
   CHARACTER (LEN=256)    :: diffrn_reflns_av_R_sigma
   CHARACTER (LEN=256)    :: reflns_number_total
   CHARACTER (LEN=256)    :: reflns_number_gt
   CHARACTER (LEN=256)    :: refine_ls_wR_factor_ref
   CHARACTER (LEN=256)    :: refine_ls_wR_factor_gt
   CHARACTER (LEN=256)    :: refine_ls_R_factor_gt
   CHARACTER (LEN=256)    :: refine_ls_R_factor_all
   CHARACTER (LEN=256)    :: refine_ls_goodness_of_fit_ref
   CHARACTER (len=256)    :: refine_ls_restrained_S_all
   CHARACTER (len=256)    :: refine_ls_shift_su_max
   CHARACTER (len=256)    :: refine_ls_shift_su_mean
   CHARACTER (LEN=256)    :: refine_diff_density_max
   CHARACTER (LEN=256)    :: refine_diff_density_min
   CHARACTER (LEN=256)    :: refine_diff_density_rms
   CHARACTER (LEN=256)    :: refine_ls_weighting_details
   CHARACTER (LEN=256)    :: atom_sites_solution_1
   CHARACTER (LEN=256)    :: atom_sites_solution_2
   CHARACTER (LEN=256)    :: atom_sites_solution_H
   CHARACTER (LEN=256)    :: refine_ls_H_treatment
   CHARACTER (LEN=256)    :: refine_ls_extinction_method
   CHARACTER (LEN=256)    :: refine_ls_extinction_coef
   CHARACTER (LEN=256)    :: refine_ls_number_reflns 
   CHARACTER (LEN=256)    :: refine_ls_number_parameters
   CHARACTER (LEN=256)    :: refine_ls_number_restraints
   CHARACTER (LEN=256)    :: H_treatment
   CHARACTER (LEN=256)    :: atom
   CHARACTER (LEN=256)    :: distance
   CHARACTER (LEN=256)    :: angle
   CHARACTER (LEN=256)    :: torsion_angle
   CHARACTER (LEN=256)    :: Hbond
   CHARACTER (LEN=256)    :: theta_min
   CHARACTER (LEN=256)    :: theta_max
   CHARACTER (LEN=256)    :: theta_full
   CHARACTER (LEN=256)    :: cell_theta_min
   CHARACTER (LEN=256)    :: cell_theta_max
   CHARACTER (LEN=256)    :: cell_reflns_used
   CHARACTER (LEN=256)    :: F000
   CHARACTER (LEN=256)    :: crystal_size_min, crystal_size_mid, crystal_size_max
   CHARACTER (LEN=256)    :: crystal_colour
   CHARACTER (LEN=256)    :: h_min, h_max, k_min, k_max, l_min, l_max
   CHARACTER (LEN=256)    :: completeness
   CHARACTER (LEN=256)    :: absorption_correction_type
   CHARACTER (LEN=256)    :: absorption_coefficient_mu
   CHARACTER (LEN=256)    :: T_min, T_max
   CHARACTER (LEN=256)    :: restraints_number
   CHARACTER (LEN=256)    :: Chi2
   CHARACTER (len=256)    :: crystal_system
   CHARACTER (LEN=256)    :: Bravais
   LOGICAL                :: WinGX_used
  END TYPE CIF_parameter_type
  TYPE (CIF_parameter_type) :: CIF_parameter
  TYPE (CIF_parameter_type) :: CIF_parameter_APEX
  TYPE (CIF_parameter_type) :: CIF_parameter_KCCD  
  TYPE (CIF_parameter_type) :: CIF_parameter_SUPERNOVA
  TYPE (CIF_parameter_type) :: CIF_parameter_X2S
  TYPE (CIF_parameter_type) :: CIF_parameter_XCALIBUR

  
  TYPE, PUBLIC:: CIF_reflns_type
   REAL                    :: theta_min
   REAL                    :: theta_max
   INTEGER                 :: reflns_used
   CHARACTER (LEN=64)      :: temperature
   REAL                    :: wavelength
  END TYPE cif_reflns_type
  TYPE (CIF_reflns_type) :: CIF_diffrn_reflns
  TYPE (CIF_reflns_type) :: CIF_cell_measurement

  TYPE, PUBLIC :: CIF_dist_type
   integer                                        :: max_text
   integer                                        :: n_text
   character (len=132), dimension(:), allocatable :: text
  END TYPE CIF_DIST_type
  type (CIF_DIST_type) :: CIF_DIST


  ! caracteristiques de la boucle _atom_site  
   integer                           :: CIF_atom_site_item_nb                ! nombre de champs dans la boucle _atom_site_ 
   integer                           :: CIF_atom_site_label_numor
   integer                           :: CIF_atom_site_type_symbol_numor      ! numero du champ
   integer                           :: CIF_atom_site_fract_x_numor
   integer                           :: CIF_atom_site_fract_y_numor
   integer                           :: CIF_atom_site_fract_z_numor
   integer                           :: CIF_atom_site_U_numor
   integer                           :: CIF_atom_adp_type_numor
   integer                           :: CIF_atom_site_calc_flag_numor        ! numero du champ
   integer                           :: CIF_atom_site_refinement_flags_numor ! numero du champ
   integer                           :: CIF_atom_site_loop_numor             ! numero de la ligne du champ _loop   
   integer                           :: CIF_atom_site_occupancy_numor
   integer                           :: CIF_atom_site_symm_multiplicity_numor 
   integer                           :: CIF_atom_site_disorder_assembly_numor
   integer                           :: CIF_atom_site_disorder_group_numor
   
   
   ! caracteristiques de la boucle _atom_type_
   integer                           :: CIF_atom_type_loop_numor
   integer                           :: CIF_atom_type_item_nb                   ! nombre de champs dans la boucle
   integer                           :: CIF_atom_type_symbol_numor              ! numero du champ _symbol dans la boucle

   ! caracteristiques de la boucle _atom_site_aniso
   integer                           :: CIF_atom_site_aniso_item_nb
   integer                           :: CIF_atom_site_aniso_label_numor
   integer                           :: CIF_atom_site_aniso_loop_numor          ! numero de la ligne du champ _loop   

end module CIF_module

!-------------------------------------------------------------------
subroutine def_laue_class()
 USE cryscalc_module, ONLY: laue_class
 implicit none

 laue_class(1)  = "triclinic                -1"
 laue_class(2)  = "monoclinic                2/m"
 laue_class(3)  = "orthorhombic              mmm"
 laue_class(4)  = "tetragonal                4/m"
 laue_class(5)  = "tetragonal                4/mmm"
 laue_class(6)  = "trigonal                  -3 (rhomb. axes)"
 laue_class(7)  = "trigonal                  -3 (hex. axes)"
 laue_class(8)  = "trigonal                  -3m (rhomb. axes)"
 laue_class(9)  = "trigonal                  -31m (hex. axes)"
 laue_class(10) = "trigonal                  -31m (hex. axes)"
 laue_class(11) = "hexagonal                 6/m (hex. axes)"
 laue_class(12) = "hexagonal                 6/mmm (hex. axes)"
 laue_class(13) = "cubic                     m3"
 laue_class(14) = "cubic                     m3m"

 return
end subroutine def_laue_class

!-------------------------------------------------------------------
subroutine def_HKL_rule()
 USE hkl_module, ONLY : HKL_rule, HKL_rule_nb
 implicit none

 HKL_rule(1)  = '  h00     h=2n+1  21 .  .  '
 HKL_rule(2)  = '  0k0     k=2n+1  . 21  .  '
 HKL_rule(3)  = '  00l     l=2n+1  .  . 21  '

 HKL_rule(4)  = '  0kl     k=2n+1  b  .  .  '
 HKL_rule(5)  = '  0kl     l=2n+1  c  .  .  '
 HKL_rule(6)  = '  0kl   k+l=2n+1  n  .  .  '

 HKL_rule(7)  = '  h0l     h=2n+1  .  a  .  '
 HKL_rule(8)  = '  h0l     l=2n+1  .  c  .  '
 HKL_rule(9)  = '  h0l   h+l=2n+1  .  n  .  '

 HKL_rule(10) = '  hk0     h=2n+1  .  .  a  '
 HKL_rule(11) = '  hk0     k=2n+1  .  .  b  '
 HKL_rule(12) = '  hk0   h+k=2n+1  .  .  n  '

 HKL_rule(13) = '  hhl     h+l=2n+1  .  .  '
 HKL_rule(14) = '  hkk     k+h=2n+1  .  .  '
 HKL_rule(15) = '  hkh     h+k=2n+1  .  .  '


 HKL_rule(16) = '  hkl   k+l=2n+1  A        '
 HKL_rule(17) = '  hkl   h+l=2n+1  B        '
 HKL_rule(18) = '  hkl   h+k=2n+1  C        '
 HKL_rule(19) = '  hkl not all odd/even F   '
 HKL_rule(20) = '  hkl h+k+l=2n+1  I        '

 HKL_rule(21) = '  h00     h=4n+1 41 .  .   '
 HKL_rule(22) = '  0k0     k=4n+1 .  41 .   '
 HKL_rule(23) = '  00l     l=4n+1 .  . 41   '

 HKL_rule(24) = '  0kl   k+l=4n+1 d  .  .   '
 HKL_rule(25) = '  h0l   h+l=4n+1 .  d  .   '
 HKL_rule(26) = '  hk0   h+k=4n+1 .  .  d   '

 HKL_rule(27) = '  h-hl, h0l, 0kl    l=2n+1 '

 HKL_rule(28) = '  hkl    h=2n+1 '
 HKL_rule(29) = '  hkl    k=2n+1 '
 HKL_rule(30) = '  hkl    l=2n+1 '

 HKL_rule(31) = '  hkl    h=2n+1 and k=2n+1'
 HKL_rule(32) = '  hkl    h=2n+1 and l=2n+1'
 HKL_rule(33) = '  hkl    k=2n+1 and l=2n+1'

! regles opposees
 HKL_rule(HKL_rule_nb+1)  = '  h00     h=2n'
 HKL_rule(HKL_rule_nb+2)  = '  0k0     k=2n'
 HKL_rule(HKL_rule_nb+3)  = '  00l     l=2n'

 HKL_rule(HKL_rule_nb+4)  = '  0kl     k=2n'
 HKL_rule(HKL_rule_nb+5)  = '  0kl     l=2n'
 HKL_rule(HKL_rule_nb+6)  = '  0kl   k+l=2n'

 HKL_rule(HKL_rule_nb+7)  = '  h0l     h=2n'
 HKL_rule(HKL_rule_nb+8)  = '  h0l     l=2n'
 HKL_rule(HKL_rule_nb+9)  = '  h0l   h+l=2n+'

 HKL_rule(HKL_rule_nb+10) = '  hk0     h=2n'
 HKL_rule(HKL_rule_nb+11) = '  hk0     k=2n'
 HKL_rule(HKL_rule_nb+12) = '  hk0   h+k=2n'

 HKL_rule(HKL_rule_nb+13) = '  hhl     h+l=2n'
 HKL_rule(HKL_rule_nb+14) = '  hkk     k+h=2n'
 HKL_rule(HKL_rule_nb+15) = '  hkh     k+k=2n'

 HKL_rule(HKL_rule_nb+16) = '  hkl   k+l=2n'
 HKL_rule(HKL_rule_nb+17) = '  hkl   h+l=2n'
 HKL_rule(HKL_rule_nb+18) = '  hkl   h+k=2n'
 HKL_rule(HKL_rule_nb+19) = '  hkl all odd/even'
 HKL_rule(HKL_rule_nb+20) = '  hkl h+k+l=2n'

 HKL_rule(HKL_rule_nb+21) = '  h00     h/=4n+1'
 HKL_rule(HKL_rule_nb+22) = '  0k0     k/=4n+1'
 HKL_rule(HKL_rule_nb+23) = '  00l     l/=4n+1'

 HKL_rule(HKL_rule_nb+24) = '  0kl   k+l/=4n+1'
 HKL_rule(HKL_rule_nb+25) = '  h0l   h+l/=4n+1'
 HKL_rule(HKL_rule_nb+26) = '  hk0   h+k/=4n+1'

 HKL_rule(HKL_rule_nb+27) = '  h-hl, h0l, 0kl    l=2n'

 HKL_rule(HKL_rule_nb+28) = '  hkl    h=2n'
 HKL_rule(HKL_rule_nb+29) = '  hkl    k=2n'
 HKL_rule(HKL_rule_nb+30) = '  hkl    l=2n'

 HKL_rule(HKL_rule_nb+31) = '  hkl    h=2n and k=2n'
 HKL_rule(HKL_rule_nb+32) = '  hkl    h=2n and l=2n'
 HKL_rule(HKL_rule_nb+33) = '  hkl    k=2n and l=2n'


end subroutine def_HKL_rule

!-----------------------------------------------------------------------------------------------
 module MATRIX_list_module
  implicit none
  integer, parameter                                               :: max_mat_nb      = 33
  integer, parameter                                               :: max_user_mat_nb = 5
  integer                                                          :: user_mat_nb
  REAL,                 DIMENSION(3,3,max_mat_nb+max_user_mat_nb)  :: transf_mat
  CHARACTER (LEN=256),  DIMENSION(max_mat_nb+max_user_mat_nb)      :: transf_mat_text
  CHARACTER (LEN=256),  DIMENSION(max_user_mat_nb)                 :: user_mat_text
  
  INTEGER                                                          :: matrix_num
  CHARACTER (len=256)                                              :: matrix_text

 end module MATRIX_list_module
!-----------------------------------------------------------------------------------------------

module USER_module
 implicit none
 integer, parameter                                 :: user_shortcut_max = 10
 integer                                            :: nb_shortcuts
 CHARACTER (LEN=32), dimension(user_shortcut_max)   :: shortcut_kw
 CHARACTER (LEN=32), dimension(user_shortcut_max)   :: shortcut_details
 
end module USER_module 
!-----------------------------------------------------------------------------------------------

subroutine def_transformation_matrix()
 USE MATRIX_list_module
 
  

 ! ' #1: a b c  ==>  a  b  c'
 ! ' #2: a b c  ==> -a -b -c'
 ! ' #3: a b c  ==>  a -c  b (orthorh. system)'
 ! ' #4: a b c  ==>  b  a -c (orthorh. system)'
 ! ' #5: a b c  ==> -c  b  a (orthorh. system)'
 ! ' #6: a b c  ==>  b  c  a (orthorh. system)'
 ! ' #7: a b c  ==>  c  a  b (orthorh. system)'
 ! ' #8: P      ==>  R      '
 ! ' #9: R      ==>  P      '
 ! '#10: R_reverse  ==>  O_obverse'
 ! '#11: R_obverse  ==>  R_reverse'
 ! '#12: P      ==>  A      '
 ! '#13: A      ==>  P      '
 ! '#14: P      ==>  B      '
 ! '#15: B      ==>  P      '
 ! '#16: P      ==>  C      '
 ! '#17: C      ==>  P      '
 ! '#18: P      ==>  F      '
 ! '#19: F      ==>  P      '
 ! '#20: P      ==>  I      '
 ! '#21: I      ==>  P      '
 ! '#22: a b c  ==> -a-b    a    c (hexa. setting)'
 ! '#23: a b c  ==>    b -a-b    c (hexa. setting)'
 ! '#24: a b c  ==>   -a  a+b   -c (hexa. setting)'
 ! '#25: a b c  ==>   -b   -a   -c (hexa. setting)'
 ! '#26: a b c  ==>  a+b   -b   -c (hexa. setting)'

 ! '#27: a b c  ==> -a-c    b    a (mono. setting)'
 ! '#28: a b c  ==>    c    b -a-c (mono. setting)'
 
 ! '#29: a b c  ==>   -a   -b  a+c (mono. setting)'
 ! '#30: a b c  ==>   -c   -b   -a (mono. setting)'
 ! '#31: a b c  ==>  a+c   -b   -c (mono. setting)'
 ! '#32: a b c  ==>    c   -b    a (mono. setting)'
 
 ! '#33: a b c  ==>    a   -b -a-c (mono. setting) : #33 = #32 x #27'
 
 ! + 5 matrices defined by the user in the setting file
 
 ! rem PLATON P21/n --> P21/c : 1  0  0    0 -1  0  -1  0 -1    #33      a  -b -a-c
 !                             -1  0  0    0 -1  0   1  0  1    #29     -a  -b  a+c
 !
 !            P21/a --> P21/c : 0  0 -1    0 -1  0  -1  0  0    #30     -c  -b  -a
 !                              0  0  1    0 -1  0   1  0  0    #32      c  -b   a 


  transf_mat_text(1)  = ' #1: a b c  ==>  a  b  c'
  transf_mat(1,:,1)  = (/ 1.,  0.,  0./)           !
  transf_mat(2,:,1)  = (/ 0.,  1.,  0./)           !
  transf_mat(3,:,1)  = (/ 0.,  0.,  1./)           !

  transf_mat_text(2)  = ' #2: a b c  ==> -a -b -c'
  transf_mat(1,:,2)  = (/-1.,  0.,  0./)           !
  transf_mat(2,:,2)  = (/ 0., -1.,  0./)           !
  transf_mat(3,:,2)  = (/ 0.,  0., -1./)           !

  !transf_mat_text(3)  = ' #3: a b c  ==>  a -c  b (orthorh. system)'
  transf_mat_text(3)  = ' #3: a b c  ==>  a -c  b'
  transf_mat(1,:,3)  = (/ 1.,  0.,  0./)           !
  transf_mat(2,:,3)  = (/ 0.,  0., -1./)           !
  transf_mat(3,:,3)  = (/ 0.,  1.,  0./)           !

  !transf_mat_text(4)  = ' #4: a b c  ==>  b  a -c (orthorh. system)'
  transf_mat_text(4)  = ' #4: a b c  ==>  b  a -c'
  transf_mat(1,:,4)  = (/ 0.,  1.,  0./)           !
  transf_mat(2,:,4)  = (/ 1.,  0.,  0./)           !
  transf_mat(3,:,4)  = (/ 0.,  0., -1./)           !

  !transf_mat_text(5)  = ' #5: a b c  ==> -c  b  a (orthorh. system)'
  transf_mat_text(5)  = ' #5: a b c  ==> -c  b  a'
  transf_mat(1,:,5)  = (/ 0.,  0., -1./)           !
  transf_mat(2,:,5)  = (/ 0.,  1.,  0./)           !
  transf_mat(3,:,5)  = (/ 1.,  0.,  0./)           !

  !transf_mat_text(6)  = ' #6: a b c  ==>  b  c  a (orthorh. system)'
  transf_mat_text(6)  = ' #6: a b c  ==>  b  c  a'
  transf_mat(1,:,6)  = (/ 0.,  1.,  0./)           !
  transf_mat(2,:,6)  = (/ 0.,  0.,  1./)           !
  transf_mat(3,:,6)  = (/ 1.,  0.,  0./)           !

  !transf_mat_text(7)  = ' #7: a b c  ==>  c  a  b (orthorh. system)'
  transf_mat_text(7)  = ' #7: a b c  ==>  c  a  b'
  transf_mat(1,:,7)  = (/ 0.,  0.,  1./)           !
  transf_mat(2,:,7)  = (/ 1.,  0.,  0./)           !
  transf_mat(3,:,7)  = (/ 0.,  1.,  0./)           !

  transf_mat_text(8)  = ' #8: P      ==>  R      '
  transf_mat(1,:,8)  = (/ 1., -1.,  0./)           !
  transf_mat(2,:,8)  = (/ 0.,  1., -1./)           !
  transf_mat(3,:,8)  = (/ 1.,  1.,  1./)           !

  transf_mat_text(9)  = ' #9: R      ==>  P      '
  transf_mat(1,:,9)  = (/  2./3.,  1./3.,  1./3./)    !
  transf_mat(2,:,9)  = (/ -1./3.,  1./3.,  1./3./)    !
  transf_mat(3,:,9)  = (/ -1./3., -2./3.,  1./3./)    !

  transf_mat_text(10) = '#12: R_reverse  ==>  O_obverse'
  transf_mat(1,:,10) = (/-1.,  0.,  0./)         !
  transf_mat(2,:,10) = (/ 0., -1.,  0./)         !
  transf_mat(3,:,10) = (/ 0.,  0.,  1./)         !

  transf_mat_text(11) = '#13: R_obverse  ==>  R_reverse'
  transf_mat(1,:,11) = (/-1.,  0.,  0./)         !
  transf_mat(2,:,11) = (/ 0., -1.,  0./)         !
  transf_mat(3,:,11) = (/ 0.,  0.,  1./)         !

  transf_mat_text(12) = '#10: P      ==>  A      '
  transf_mat(1,:,12) = (/-1.,  0.,  0./)         !
  transf_mat(2,:,12) = (/ 0., -1.,  1./)         !
  transf_mat(3,:,12) = (/ 0.,  1.,  1./)         !

  transf_mat_text(13) = '#11: A      ==>  P      '
  transf_mat(1,:,13) = (/-1.,  0.,  0./)         !
  transf_mat(2,:,13) = (/ 0., -0.5, 0.5/)        !
  transf_mat(3,:,13) = (/ 0.,  0.5, 0.5/)        !

  transf_mat_text(14) = '#14: P      ==>  B      '
  transf_mat(1,:,14) = (/-1.,  0.,  1./)         !
  transf_mat(2,:,14) = (/ 0., -1.,  0./)         !
  transf_mat(3,:,14) = (/ 1.,  0.,  1./)         !

  transf_mat_text(15) = '#15: B      ==>  P      '
  transf_mat(1,:,15) = (/-0.5,  0.,  0.5/)       !
  transf_mat(2,:,15) = (/ 0. , -1.,  0./)        !
  transf_mat(3,:,15) = (/ 0.5,  0.,  0.5/)       !

  transf_mat_text(16) = '#16: P      ==>  C      '
  transf_mat(1,:,16) = (/ 1.,  1.,  0./)         !
  transf_mat(2,:,16) = (/ 1., -1.,  0./)         !
  transf_mat(3,:,16) = (/ 0.,  0., -1./)         !

  transf_mat_text(17) = '#17: C      ==>  P      '
  transf_mat(1,:,17) = (/ 0.5,  0.5, 0./)        !
  transf_mat(2,:,17) = (/ 0.5, -0.5, 0./)        !
  transf_mat(3,:,17) = (/ 0.,   0., -1./)        !

  transf_mat_text(18) = '#18: P      ==>  F    (ex : rhomb. ==> cubic F)  '
  transf_mat(1,:,18) = (/-1.,  1.,  1./)         !
  transf_mat(2,:,18) = (/ 1., -1.,  1./)         !
  transf_mat(3,:,18) = (/ 1.,  1., -1./)         !

  transf_mat_text(19) = '#19: F      ==>  P    (ex : cubic F ==> rhomb.)  '
  transf_mat(1,:,19) = (/ 0.,  0.5, 0.5/)        !
  transf_mat(2,:,19) = (/ 0.5, 0. , 0.5/)        !
  transf_mat(3,:,19) = (/ 0.5, 0.5, 0./)         !

  transf_mat_text(20) = '#20: P      ==>  I    (ex : rhomb.  ==> cubic I)  '
  transf_mat(1,:,20) = (/ 0.,  1.,  1./)         !
  transf_mat(2,:,20) = (/ 1.,  0.,  1./)         !
  transf_mat(3,:,20) = (/ 1.,  1.,  0./)         !

  transf_mat_text(21) = '#21: I      ==>  P    (ex : cubic I ==> rhomb.)  '
  transf_mat(1,:,21) = (/-0.5,  0.5,  0.5/)      !
  transf_mat(2,:,21) = (/ 0.5, -0.5,  0.5/)      !
  transf_mat(3,:,21) = (/ 0.5,  0.5, -0.5/)      !


  transf_mat_text(22) = '#22: a b c  ==> -a-b    b  c (hexa. setting)'
  transf_mat(1,:,22) = (/-1.,  -1.,  0./)      !
  transf_mat(2,:,22) = (/ 1.,   0.,  0./)      !
  transf_mat(3,:,22) = (/ 0.,   0.,  1./)      !

  transf_mat_text(23) = '#23: a b c  ==>    b -a-b  c (hexa. setting)'
  transf_mat(1,:,23) = (/ 0.,   1.,  0./)      !
  transf_mat(2,:,23) = (/-1.,  -1.,  0./)      !
  transf_mat(3,:,23) = (/ 0.,   0.,  1./)      !

  transf_mat_text(24) = '#24: a b c  ==>   -a  a+b -c (hexa. setting)'
  transf_mat(1,:,24) = (/-1.,   0.,  0./)      !
  transf_mat(2,:,24) = (/ 1.,   1.,  0./)      !
  transf_mat(3,:,24) = (/ 0.,   0., -1./)      !

  transf_mat_text(25) = '#25: a b c  ==>   -b   -a -c (hexa. setting)'
  transf_mat(1,:,25) = (/ 0.,  -1.,  0./)      !
  transf_mat(2,:,25) = (/-1.,   0.,  0./)      !
  transf_mat(3,:,25) = (/ 0.,   0., -1./)      !

  transf_mat_text(26) = '#26: a b c  ==>  a+b   -b -c (hexa. setting)'
  transf_mat(1,:,26) = (/ 1.,   1.,  0./)      !
  transf_mat(2,:,26) = (/ 0.,  -1.,  0./)      !
  transf_mat(3,:,26) = (/ 0.,   0., -1./)      !

  transf_mat_text(27) = '#27: a b c  ==> -a-c    b    a (mono. setting)'
  transf_mat(1,:,27) = (/-1.,   0., -1./)      !
  transf_mat(2,:,27) = (/ 0.,   1.,  0./)      !
  transf_mat(3,:,27) = (/ 1.,   0.,  0./)      !

  transf_mat_text(28) = '#28: a b c  ==>    c    b -a-c (mono. setting)'
  transf_mat(1,:,28) = (/ 0.,   0.,  1./)      !
  transf_mat(2,:,28) = (/ 0.,   1.,  0./)      !
  transf_mat(3,:,28) = (/-1.,   0., -1./)      !

  transf_mat_text(29) = '#29: a b c  ==>   -a   -b  a+c (mono. setting)'
  transf_mat(1,:,29) = (/-1.,   0.,  0./)      !
  transf_mat(2,:,29) = (/ 0.,  -1.,  0./)      !
  transf_mat(3,:,29) = (/ 1.,   0.,  1./)      !

  transf_mat_text(30) = '#30: a b c  ==>   -c   -b   -a (mono. setting)'
  transf_mat(1,:,30) = (/ 0.,   0., -1./)      !
  transf_mat(2,:,30) = (/ 0.,  -1.,  0./)      !
  transf_mat(3,:,30) = (/-1.,   0.,  0./)      !

  transf_mat_text(31) = '#31: a b c  ==>  a+c   -b   -c (mono. setting)'
  transf_mat(1,:,31) = (/ 1.,   0.,  1./)      !
  transf_mat(2,:,31) = (/ 0.,  -1.,  0./)      !
  transf_mat(3,:,31) = (/ 0.,   0., -1./)      !

  transf_mat_text(32) = '#32: a b c  ==>    c   -b   a (mono. setting)'
  transf_mat(1,:,32) = (/ 0.,   0.,  1./)      !
  transf_mat(2,:,32) = (/ 0.,  -1.,  0./)      !
  transf_mat(3,:,32) = (/ 1.,   0.,  0./)      !

  transf_mat_text(33) = '#33: a b c  ==>    a   -b  -a-c (mono. setting)'
  transf_mat(1,:,33) = (/ 1.,   0.,  0./)      !
  transf_mat(2,:,33) = (/ 0.,  -1.,  0./)      !
  transf_mat(3,:,33) = (/-1.,   0., -1./)      !

  ! initialisation des matrices utilisateurs
  do i=1, 5
   write(unit=transf_mat_text(max_mat_nb+i), fmt='(a,i2,a,i1)')  '#', max_mat_nb+i, ': user matrix #', i
   write(unit=user_mat_text(i),              fmt='(a,i2,a,i1)')  '#u', i, ': user matrix #', i   
   transf_mat(1,:,max_mat_nb+i) = (/ 1.,   0.,  0./)      !
   transf_mat(2,:,max_mat_nb+i )= (/ 0.,   1.,  0./)      !
   transf_mat(3,:,max_mat_nb+i) = (/ 0.,   0.,  1./)      !
  end do
  
  
   
end subroutine def_transformation_matrix


!-------------------------------------------------------------------------------------
module Definition_fractions
  implicit none


  integer, parameter                          :: nb_fraction = 35
  CHARACTER (LEN=12), DIMENSION(nb_fraction)  :: ratio_string_neg, ratio_string_pos
  real,               DIMENSION(nb_fraction)  :: ratio_real_neg,   ratio_real_pos


  contains
  subroutine def_fractions

  ratio_string_pos(1:nb_fraction) = (/'1/2 ',                                                 &
                                      '1/3 ', '2/3 ', '4/3 ', '5/3 ',                         &
                                      '1/4 ', '3/4 ', '5/4 ',                                 &
                                      '1/5 ', '2/5 ', '3/5 ', '4/5 ', '6/5 ',                 &
                                      '1/6 ', '5/6 ', '7/6 ',                                 &
                                      '1/7 ', '2/7 ', '3/7 ', '4/7 ', '5/7 ', '6/7 ', '8/7 ', &
                                      '1/8 ', '3/8 ', '5/8 ', '7/8 ', '9/8 ',                 &
                                      '1/9 ', '2/9 ', '4/9 ', '5/9 ', '7/9 ', '8/9 ', '10/9' /)

  ratio_string_neg(1:nb_fraction) = (/'-1/2 ',                                                        &
                                      '-1/3 ', '-2/3 ', '-4/3 ', '-5/3 ',                             &
                                      '-1/4 ', '-3/4 ', '-5/4 ',                                      &
                                      '-1/5 ', '-2/5 ', '-3/5 ', '-4/5 ', '-6/5 ',                    &
                                      '-1/6 ', '-5/6 ', '-7/6 ',                                      &
                                      '-1/7 ', '-2/7 ', '-3/7 ', '-4/7 ', '-5/7 ', '-6/7 ', '-8/7 ',  &
                                      '-1/8 ', '-3/8 ', '-5/8 ', '-7/8 ', '-9/8 ',                    &
                                      '-1/9 ', '-2/9 ', '-4/9 ', '-5/9 ', '-7/9 ', '-8/9 ', '-10.9'  /)

  ratio_real_pos(1:nb_fraction) = (/  1./2.,                                           &
                                      1./3., 2./3., 4./3., 5./3.,                      &
                                      1./4., 3./4., 5./4.,                             &
                                      1./5., 2./5., 3./5., 4./5., 6./5.,               &
                                      1./6., 5./6., 7./6.,                             &
                                      1./7., 2./7., 3./7., 4./7., 5./7., 6./7., 8./7., &
                                      1./8., 3./8., 5./8., 7./8., 9./8.,               &
                                      1./9., 2./9., 4./9., 5./9., 7./9., 8./9., 10./9. /)
ratio_real_neg(1:nb_fraction) = -ratio_real_pos(1:nb_fraction)

  return
  end subroutine def_fractions

end Module Definition_Fractions


!-------------------------------------------------------------------------------------------------------------------

Module Accents_module
 implicit none
  integer, parameter                         :: nb_char = 14
  CHARACTER (LEN=8), DIMENSION(4, nb_char)   :: accents
  integer                                    :: n
  


  contains

  subroutine def_accents
    !                            caractere       format_CIF       HTML            HTML decimal      
	
    n=0
    n=n+1;  accents(1:4, n) = (/ "à       ",    "\`a     ",      "&agrave;",     "&#224;  "     /)
    n=n+1;  accents(1:4, n) = (/ "â       ",    "\^a     ",      "&acirc; ",     "&#226;  "     /)

    n=n+1;  accents(1:4, n) = (/ "é       ",    "\'e     ",      "&eacute;",     "&#233;  "     /)
    n=n+1;  accents(1:4, n) = (/ "è       ",    "\`e     ",      "&egrave;",     "&#232;  "     /)
    n=n+1;  accents(1:4, n) = (/ "ê       ",    "\^e     ",      "&ecirc; ",     "&#234;  "     /)
    n=n+1;  accents(1:4, n) = (/ "ë       ",    '\"e     ',      "&euml;  ",     "&#235;  "     /)

    n=n+1;  accents(1:4, n) = (/ "î       ",    "\^i     ",      "&icirc; ",     "&#238;  "     /)
    n=n+1;  accents(1:4, n) = (/ "ï       ",    '\"i     ',      "&iuml;  ",     "&#239;  "     /)

    n=n+1;  accents(1:4, n) = (/ "ô       ",    "\^o     ",      "&ocirc; ",     "&#244;  "     /)

    n=n+1;  accents(1:4, n) = (/ "ù       ",    "\`u     ",      "&ugrave;",     "&#249;  "     /)
    n=n+1;  accents(1:4, n) = (/ "û       ",    "\^u     ",      "&ucirc; ",     "&#251;  "     /)
    n=n+1;  accents(1:4, n) = (/ "ü       ",    '\"u     ',      "&uuml;  ",     "&#252;  "     /)

    n=n+1;  accents(1:4, n) = (/ "ç       ",    "\,c     ",      "&ccedil;",     "&#231;  "     /)
    n=n+1;  accents(1:4, n) = (/ "ñ       ",    "\~n     ",      "&ntilde;",     "&#241;  "     /)

	!exposant
    !n=n+1;  accents(1:4, n) = (/ "expo    ",    "^       ",      "<sup>   ",     "        "     /)
	
	

    return
  end subroutine def_accents
									  
End Module Accents_module


