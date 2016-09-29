!***********************************************************************
!                                                                      *
!     Copyright 1987-2002 Michael M. J. Treacy and Michael W. Deem     *
!                                                                      *
!***********************************************************************
!***********************************************************************
!*******************      Source file DIFFaX.f       *******************
!***********************************************************************
!***********************************************************************
!******************** version 1.813, 19th May, 2010 ********************
!***********************************************************************
!***********************************************************************
! This program calculates the powder diffraction spectrum of a crystal *
! formed from layers that stack coherently, but not deterministically. *
! The algorithm used is described in detail in "A General Recursion    *
! Method for Calculating Diffracted Intensities From Crystals          *
! Containing Planar Faults", by M.M.J. Treacy, M.W. Deem and           *
! J.M. Newsam, Proceedings of the Royal Society, London, A (1991) 433, *
! pp 499 - 520.                                                        *
!                                                                      *
! Source code written by Michael M. J. Treacy and Michael W. Deem.     *
!                                                                      *
! HISTORY:                                                             *
! 6/87-5/88; MMJT: Original versions were named 'betaPXD' and 'FauPXD'.*
! They were 'hardwired' for simulating PXD patterns from zeolite beta, *
! and the faujasite/'Breck's structure 6' zeolite family.              *
!                                                                      *
! 5-8/88; MWD: Completely rewritten in generalized form by Mike Deem   *
! for Cray and VAX, and named 'PDS'.                                   *
!                                                                      *
! 8/88-11/89; MMJT: Deem's code was extended and mostly rewritten.     *
! Control file option added. Improved handling of sharp peaks.         *
! Symmetry testing added. Layer 'stacking uncertainty' factors added.  *
! Selected area electron diffraction option added.                     *
! Explicit layer sequencing option. Optimization of layer form factor  *
! calculations. Renamed 'DIFFaX' for, (D)iffraction (I)ntensities      *
! (F)rom (Fa)ulted (X)tals.                                            *
!                                                                      *
! 12/89 - 3/90; MMJT: Finite crystal thickness now accepted under      *
!              'RECURSIVE' option. Self-consistency check of atomic    *
!               coordinates in data file. (v1.76)                      *
!                                                                      *
! 4/90 - 12/90; MMJT: Bug fixes. Streamlined 'data.sfc' file.          *
!               (v1.761 and v1.762)                                    *
!                                                                      *
! 1/91; MMJT: Eliminated the use of scratch file while reading data.   *
!             GETLNE now handles multiple, nested, comments (v1.763)   *
!                                                                      *
! 5/91; MMJT: Eliminated bug in default tolerance parameter. Added     *
!             average crystal composition printout to dump file.       *
!             (v1.764)                                                 *
!                                                                      *
! 8/91; MMJT: Replaced the LU decomposition routines CLUDCM and CLUBKS *
!             (from "Numerical Recipes") with the faster linpack       *
!             routines, CGEFA and CGESL (v1.765)                       *
!                                                                      *
! 8/91; MMJT: Improved sharp peak detection - use peak widths rather   *
!             than the more complicated phase coherence argument in    *
!             subroutine SHARP (v1.766)                                *
!                                                                      *
! 4/92; MMJT: Fixed bug in INTEN2 where last layer in an explicit      *
!             sequence was inadvertently assigned a scattering factor  *
!             of C_ONE. Improved error checking for explicit layers so *
!             that when alpha(j,i) = 0, an error is issued if j        *
!             follows i. XP_MAX increased to 5000. GETLNE now checks   *
!             that data lines do not exceed maximum length. (v1.767)   *
!                                                                      *
! 12/94; MMJT:Reinstated the use of the scratch file that had been     *
!             eliminated in v1.763. The Cray, and Microsoft fortran    *
!             compiler for PC, adhere to the FORTRAN77 standard and    *
!             do not allow unformatted reads/writes from/to strings.   *
!                                                           (v1.768)   *
!                                                                      *
! 1/95; MMJT: Finessed the diffraction symmetry detection routines.    *
!             Introduced the subroutine THRESH. (v1.769)               *
!                                                                      *
! 2/95; MMJT: Fixed glitches in THRESH, TST_MIR and TST_ROT that were  *
!             introduced in the 1/95 fix. (still v1.769)               *
!                                                                      *
! 3/95; MMJT: Implemented Debye-Scherrer type broadening due to finite *
!             lateral layer widths. Added CHWDTH, RDWDTH. Modified the *
!             way the powder pattern is written to the array spec, so  *
!             that the low angle intensity begins at spec(1). Atom     *
!             names are now case insensitive. The layer Bij factors    *
!             were reordered - B23, B31 become B13, B23.               *
!             GET_G modified to handle singularities better. (v1.80)   *
!                                                                      *
! 7/95; MMJT: Fixed rare zero integration range bug in GLQ16.          *
!             Fixed "fatsWalla" bug in GET_MAT.  (v1.801)              *
!                                                                      *
! 5/96; MMJT: Fixed a bug in the LL() function in INTEGR() which was   *
!             introduced by a "cosmetic" change in v1.801.  (v1.802)   *
!                                                                      *
! 10/96; MMJT:Changed eps3 to eps5 in CHWDTH function so that the      *
!             broadening tails extend further.  (v1.803)               *
!                                                                      *
! 7/97; MMJT: Added subroutines NXTARG and RDNMBR. These allow data    *
!             to be entered as fractions (ie 1/3). Improved robustness *
!             of the "fatswalla" interlayer uncertainty code. (v1.804) *
!                                                                      *
! 6/98; MMJT: Fixed bug in PV() that was introduced in v1.80. (v1.805) *
!                                                                      *
! 3/00; MMJT: Now allow 16-bit deep SADPs. (v1.806)                    *
!                                                                      *
! 8/00; MMJT: RDSTAK changed so that if a stacking probability is zero,*
!             the rest of the line is ignored.        (v1.807)         *
!                                                                      *
! 4/02; MMJT: Removed the non-standard calls to iidnnt in WRTSADP.     *
!                                                      (v1.808)        *
!                                                                      *
! 2/03; MMJT: Halved the value of ffhkcnst in GETSPC. The half-width   *
!             was being used instead of the FWHM, which made the shape *
!             broadening twice as large as it should be. (v1.809)      *
!                                                                      *
! 3/04; MMJT: Fixed a minor printing bug in TST_MIR.     (v1.810)      *
!                                                                      *
! 1/05; MMJT: Fixed some f77 compiler compatibility bugs.  (v1.811)    *
!                                                                      *
! 7/05; MMJT: Fixed a bug in EQUALB that caused DIFFaX to ignore the   *
!             sign of the Fats-Waller Bij terms.  (v1.812)             *
!                                                                      *
! 5/10; MMJT: Changed the way WRTSADP writes big-endian 16-bit data.   *
!             Improved scratch-file clean-up.                (v1.813)  *
!                                                                      *
!    WRTSADP has been removed in the present version for FAULTS (JRC)  *                                                       *
!                                                                      *
!***********************************************************************
!***************************** Legal Note ******************************
!***********************************************************************
!                                                                      *
! * * * * * * * * * *  DISCLAIMER OF WARRANTIES: * * * * * * * * * * * *
!                                                                      *
! The authors make no warranties whatsoever, express or implied, with  *
! respect to the DIFFaX software or any of its parts, nor do they      *
! warrant that the DIFFaX software, or any of its parts, will be       *
! error-free, will operate without interruption, or will be compatible *
! with any software or hardware possessed by the user.                 *
!                                                                      *
! * * * * * * * * * *  LIMITATION OF LIABILITY:  * * * * * * * * * * * *
!                                                                      *
! The authors will not be liable for any special, incidental, or       *
! consequential damages, even if informed of the possibility of such   *
! damages in advance.                                                  *
!                                                                      *
!***********************************************************************
!************************** DIFFaX file i/o. ***************************
!***********************************************************************
!                                                                      *
! * * * * OPTIONAL CONTROLFILE FOR AUTOMATIC RUNNING OF DIFFaX * * * * *
!                                                                      *
! DIFFaX first searches the current directory for a control file named *
! 'control.dif'. If it finds this file it opens it on unit 'cntrl'     *
! and this becomes the default input unit. Structure filenames,        *
! and the various parameters (which would normally be requested        *
! interactively) obviously must be in the correct sequence. The data   *
! read from 'control' is echoed on the default output device (ie. the  *
! screen, unit 'op') so the user can check that the responses          *
! are properly synchronized. If 'control.dif' does not exist, the      *
! default input device is the keyboard (unit number 'ip'), and the     *
! user is expected to answer the prompts. DIFFaX will loop through the *
! contents of 'control', and thus can be used to rerun DIFFaX on fresh *
! data files, without quitting. Under direction from a control file,   *
! normal termination of DIFFaX occurs when a filename 'END' is         *
! encountered. Interactively, DIFFaX will end normally when the user   *
! chooses not to return to the function menu.                          *
! The name of the control file is stored in the global character       *
! variable 'cfname', and is assigned in 'main'.                        *
!                                                                      *
!                                                                      *
! * * * * * * * * * * *  STRUCTURE INPUT FILE  * * * * * * * * * * * * *
!                                                                      *
! The structure input file is opened on unit 'df'. It can have any     *
! name except 'END' (case insensitive). For clarity it may be best to  *
! keep it short (less than 8 characters) and (optionally) with '.dat'  *
! appended. Output files use the input name up to the first blank      *
! (' ') or period ('.') as their root name. Thus, if 'beta.dat'        *
! (or 'beta') is the data input file name, then 'beta.spc' etc... will *
! be the form of the output file names.                                *
!                                                                      *
!                                                                      *
! * * * * * * * STRUCTURE FACTOR PARAMETER INPUT FILE  * * * * * * * * *
!                                                                      *
! The structure factor parameter file, 'data.sfc' is opened on unit    *
! 'sf'. If a file of name 'data.sfc' is not found, DIFFaX will abort.  *
! The name of the structure factor parameter file is stored in the     *
! global character variable 'sfname', and is assigned in 'main'.       *
!                                                                      *
!                                                                      *
! * * * * * * * * * * *  SPECTRUM OUTPUT FILE  * * * * * * * * * * * * *
!                                                                      *
! Spectra are output as text files on unit 'sp'. Each record contains  *
!        2theta     intensity     (instrumentally broadened intensity) *
! in tab-delimited format. 'Instrumentally broadened intensity' is     *
! output only if the pseudo-Voigt, Gaussian or Lorentzian options were *
! requested. Spectra output file names are in the form 'rootname.spc', *
! or alternatively, if that name is already taken, as 'rootname.spc#', *
! where #=1,2,3 etc...                                                 *
!                                                                      *
!                                                                      *
! * * * * * * * * * STREAK INTENSITIES OUTPUT FILE * * * * * * * * * * *
!                                                                      *
! Streak calculations are output on unit 'sk'. Streak output file      *
! names are in the form 'rootname.str', or alternatively, if that name *
! is already taken, as 'rootname.str#', where #=1,2,3 etc...           *
!                                                                      *
!                                                                      *
! * * *   SELECTED AREA DIFFRACTION PATTERN (SADP) OUTPUT FILE   * * * *
!                                                                      *
! Selected area diffraction pattern data is saved in binary format in  *
! a file named 'rootname.sadp' which is output on unit 'sad'. If that  *
! name is already taken, the alternative name 'rootname.sadp#' is      *
! used, where #=1,2,3 etc...                                           *
!                                                                      *
!                                                                      *
! * * * * * *  OPTIONAL DUMP FILE OF STRUCTURAL PARAMETERS * * * * * * *
!                                                                      *
! If the user requests a dump of the structure data file (as DIFFaX    *
! read it!) a dumpfile named 'rootname.dmp' is output on unit dmp. If  *
! that name is already taken, the alternative name 'rootname.dmp#' is  *
! used, where #=1,2,3 etc...This is valuable for debugging the input   *
! data file.                                                           *
!                                                                      *
!                                                                      *
! * * * * * * *  OPTIONAL DUMP FILE OF INTENSITIES FOUND   * * * * * * *
! * * * * * * * WHEN EVALUATING DIFFRACTION POINT SYMMETRY * * * * * * *
!                                                                      *
! The user may also output the history of the intensity values found   *
! when DIFFaX attempts to establish the point group symmetry of the    *
! diffraction output. This is useful when debugging the datafile. The  *
! intensity data is saved in a file named 'rootname.sym' which is      *
! output on unit 'sy'. If that name is already taken, the alternative  *
! name 'rootname.sym#' is used, where #=1,2,3 etc...                   *
!                                                                      *
!***********************************************************************
!******************************* DIFFaX ********************************
!***********************************************************************
!
! ______________________________________________________________________
! Title: DIFFaX
! Authors: MWD and MMJT
! Date: 23 Oct 1988
! Description: This is the main program. First, important global
! constants, such as PI, are defined. The name of the control file
! is assigned to cfname, and then FNDCTL searches for this file in the
! current directory. If found, the control file is opened and it
! becomes the default input device. If not found, then the keyboard is
! the standard input device. The user's data file, and the atomic
! scattering factor data file (whose name is contained in 'sfname')
! are then searched for in the current directory (GETFIL), and opened.
! The user's data file is then read (RDFILE). The standard scattering
! factor data file 'sfname' is then searched for data on the atom
! types specified by the user (SFC). The layer existence
! probabilities are calculated (GET_G). If the user data file
! requested EXPLICIT, RANDOM stacking, then DIFFaX computes a random
! layer sequence consistent with the stacking probabilities (GETLAY).
! Reciprocal lattice constants related to the unit cell are then
! calculated (SPHCST). If the user requested (either interactively,
! or through the control file) a dump of what DIFFaX read from the
! user's data file, then an annotated dump is generated (DUMP). DETUN
! then delicately adjusts the probability data so as to avoid zero
! determinants at the sharp peaks. The user is then asked if a dump
! of DIFFaX's symmetry evaluations is required, and then searches the
! data looking for simple opportunities to speed up the calculation
! (OPTIMZ). The user is then asked if he wants to calculate the
! intensity at a point (POINT), along a streak (GOSTRK), integrated
! within a defined interval (GOINTR), a powder pattern (GOSPEC) or
! a selected area diffraction pattern (GOSADP). If running
! interactively, the user can return to any of these menu options,
! except if GOSPEC was chosen, where DIFFaX will finish. If a control
! file is being used, then DIFFaX will return to the beginning if
! GOSPEC was chosen. If a new data file name is read then DIFFaX will
! run again. If the control file reads 'End' (case insensitive) as the
! new file name, then DIFFaX will finish.
! Note: The file names contained in 'cfname' and 'sfname', and the name
! 'End' are reserved names, and cannot be used by the user as data file
! names.
!
!      COMMON VARIABLES:
!            uses:  rndm, cntrl, CFile, SymGrpNo
!
!        modifies:  PI, PI2, DEG2RAD, RAD2DEG, DoDatdump,
!                   DoSymDump, cfname, sfname
! ______________________________________________________________________
!


 Module diffax_mod

  use CFML_GlobalDeps, only :  cp, sp, dp

!***********************************************************************
!*****************      Declaration of parameters      *****************
!***********************************************************************

      Implicit None
      Integer, Parameter ::   max_l=20, &    !MAX_L  -  the maximum number of layer types allowed
                              max_a=200,&    !MAX_A  -  the maximum number of atoms per layer
                              max_ta=20,&    !MAX_TA -  the maximum number of different atom types
                              max_sp=200001,& !MAX_SP -  the maximum number of points in the spectrum
                              max_cyc=20,&   !MAX_CYC-  the maximum number of iteration cycles
                              max_npar= 300, &      ! maximum number of parameters to refine
                              sadsize=256    !SADSIZE - the array size for the selected area diffraction pattern

      Integer, Parameter :: xp_max=5000, &   !XP_MAX   -  the maximum number of layers that can be
                                             !            explicitly sequenced non-recursively.
                              rcsv_max=1022,&!RCSV_MAX -  the maximum number of layers that can be
                                             !            explicitly sequenced recursively. RCSV_MAX should
                                             !            not exceed 2**MAX_BIN - 2
                              max_nam=31,&   !MAX_NAM   -  the maximum allowable number of characters in a filename
                              max_bin=10     !MAX_BIN   -  Maximum number of 'bits' to be used in binary
                                             !             representation of RCRSV_MAX+2

      Integer, Parameter :: ffact_size=201,&   !FFACT_SIZE-  Array size for pre-computed Lorentzian used for
                                               !             computing the lateral (a-b) size broadening.
                              n_sigmas=7       !N_SIGMAS  -  Number of half-widths to compute Lorentzian. The
                                               !             remainder of the array goes linearly to zero.

      Real(kind=dp),    Parameter :: inf_width=1.0D4  !inf_width -  Layer width in Angstroms that DIFFaX considers to
                                               !             be infinite, with no detectable size broadening
      Integer, Parameter ::   ip=5, &   !ip  -  standard input device number
                              op=6, &   !op  -  standard output device number
                              df=2, &   !df  -  unit that the structure data file will be read from
                              sf=4, &   !sf  -  unit that the standard scattering factor data will be read from
                              sy=11,&   !sy  -  unit that the symmetry evaluation data will be written to
                              unit_sp=12,&    !sk  -  unit that the streak data will be written to.
                              sk=13,&         !sp  -  unit that the spectrum data will be written to.
                              sa=14 ,&        !sa  -  unit that the 8-bit binary formatted selected area diffraction pattern data will be written to.
                              i_flts = 27, &  !Output flst file unit
                              i_fst = 28      !logical unit of FST file

      Integer, Parameter :: scrtch = 3      !scrtch  -  unit that the scratch file will be used on

      Integer, Parameter :: clip = 14       !CLIP    -  allowed length of filename appendages

      Integer, Parameter :: max_bckg_points=100               !Maximum of background points

      Integer, Parameter :: UNKNOWN = -1    !UNKNOWN -  flag indicating whether or not the symmetry
                                            !           point-group has been defined by the user
! define some useful numerical constants

      Complex(kind=dp), Parameter :: c_zero = (0.0D0,0.0D0), c_one = (1.0D0,0.0D0)
      Real(kind=dp),    Parameter :: zero = 0.0D0, quarter = 0.25D0, half = 0.5D0,     &
                               one = 1.0D0, two = 2.0D0, three = 3.0D0, four = 4.0D0,  &
                               five = 5.0D0, six = 6.0D0, eight = 8.0D0, ten = 10.0D0, &
                               twelve = 12.0D0, twenty = 20.0D0, fifty = 50.0D0,       &
                               hundred = 100.0D0, one_eighty = 180.0D0, two_fiftyfive = 255.0D0
      Real(kind=dp),    Parameter :: eps1 = 1.0D-1, eps2 = 1.0D-2, eps3 = 1.0D-3,  &
                               eps4 = 1.0D-4, eps5 = 1.0D-5, eps6 = 1.0D-6,  &
                               eps7 = 1.0D-7, eps8 = 1.0D-8, eps9 = 1.0D-9,  &
                               eps10 = 1.0D-10, eps14 = 1.0D-14

      Real(kind=dp),    Parameter :: EIGHTBITS= 256.0D0, FIFTEENBITS= 32768.0D0, SIXTEENBITS= 65536.0D0
!
!************   Description of variables
!
!
!   d-> indicates that the variable is read in from the user's data file.
!
!   i-> indicates that the variable is acquired interactively at run_time, either by
!       prompting the user for keyboard input, or (if present) by reading the control file.
!
!   s-> indicates that the variable is read in from the atomic scattering factor data file.
!       Normally, this data should never change.
!
!
!**********************************************************************
!**************      Declaration of COMMON variables      **************
!***********************************************************************

    ! save
!
!*********************     character variables
!

  Character (Len=4), dimension(max_ta)     :: atom_l  ! The name of each type of atom found in the structure data file
  Character (Len=4), dimension(max_a,max_l):: a_name  !d-> Name of each atom. DIFFaX expects 4 characters.
                                                      !    See file 'data.sfc' for allowed names
  Character (Len=12) :: pnt_grp     !d->  Symbolic name of the point group symmetry of
                                    !     the diffraction data.
  Character (Len=256) :: sfname, &  !     The name of the data file (including the path) containing the
                                    !     atomic scattering factor data (usually set to 'data.sfc')
                        cfname      !     The name of the control file containing the automated responses to
                                    !     DIFFaX's prompts (usually set to 'control.dif')

  Character (Len=31)    :: infile, outfile, outfile_notrepl
  Character (Len=31), dimension(:), allocatable    :: file_streak

  Character (Len=20), dimension(max_npar):: namepar
  Character (Len = 20)                   :: dfile, fmode  ! name of diffraction pattern file, file mode as in  winplotr
  Character (Len = 20)                   :: background_file  !name of backfround file
  Character (Len=20)                     :: mode             ! needed in read_background_file in case of interpolation
  Character (Len=20)                     :: lstype !type of explicit stacking
  Character (Len=132)                    :: ttl  !title
!
!**********************     logical variables
!
  Logical, dimension(max_l) :: one_b, &   !  TRUE if all Debye-Waller factors are the same
                               one_occup  !  TRUE if all occupancy factors are the same

  Logical, dimension(max_l, max_l) :: bs_zero,  & !TRUE if all layer stacking uncertainty factors in
                                                  !one transition are all zero.
                                      there       !TRUE if the transition j to i is non-zero



  Logical :: replace_files=.false., & !if the user makes this variable true with REPLACE_FILES command old files are overwritten
             only_real,   &   !TRUE if all layers are centrosymmetric
             same_bs,     &   !TRUE if all layer stacking uncertainty factors are identical
             all_bs_zero, &   !TRUE when all Bs_zero(i,j) = TRUE
             rot_only,    &   !TRUE if diffraction point group symmetry has no vertical mirrors
             cfile,       &   !TRUE if there is a file named 'control.dif' in the current directory.
                              !     This is the so-called 'control file' that automates the running of DIFFaX
             dodatdump,   &   !TRUE if the user wants a dump of the data file
             dosymdump,   &   !TRUE if the user wants to dump the output of the symmetry tests
             intp_f,      &   !TRUE if the form factors are to be interpolated for the purposes of the
                              !     Gauss-Legendre adaptive quadrature integration.
             trim_origin, &   !If TRUE, intensity near the origin will be ignored for the purposes of
                              !applying instrumental broadening
             recrsv,      &   !d-> TRUE if the recursive model is to be used to calculate diffraction
                              !    from a statistical ensemble of crystallites
             xplcit,      &   !d-> TRUE if the user specifies a layer sequence explicitly
                              !    'recrsv' and 'xplcit' cannot both be TRUE
             rndm,        &   !d-> TRUE if 'xplcit' is TRUE, and if the user wishes the computer to generate
                              !    a random sequence of layers biassed by the transition probabilities 'l_alpha'
             inf_thick,   &   !TRUE if crystal is to be treated as if it has an infinite number of layers
             has_l_mirror,&   !TRUE if the diffraction data has a mirror perpendicular to c* (along fault axis)
             h_mirror,    &   !TRUE if diffraction data has a mirror across the a*-c* plane
             k_mirror,    &   !TRUE if diffraction data has a mirror across the b*-c* plane
             hk_mirror        !TRUE if diffraction data has a mirror across the a*=b*, c* plane

  Logical :: check_sym,   &   !TRUE if user specified a point group symmetry to test against. If the user's
                              !     choice is incompatible with the input data, then this will be reset to FALSE
             same_rz,     &   !TRUE if all stacking vectors have the same z-component
             any_sharp,   &   !TRUE if DIFFaX suspects there are any sharp peaks
             same_layer,  &   !TRUE if all of the explicitly defined layers are identical
             finite_width     !d-> TRUE if layer widths are finite

  logical,  dimension(:), allocatable  :: fundamental  !array which defines if each layer is fundamental or not
  logical                              :: randm, semirandm, spcfc

  logical                              :: streakOrPowder = .false. !if TRUE then streak or powder fitting is handled
  logical                              :: unbroaden = .false. !if TRUE then streak calculates without broadening
  logical, dimension(10)               :: table = .false.  !if TRUE then a vector component of l_r and/or probabilities are written in the form of table
!
!*********************     integer*4 variables
!
  Integer, dimension(MAX_A,MAX_L) :: a_number   !d-> numeric label of atom in the layer
  Integer, dimension(MAX_A,MAX_L) :: a_type     !d-> type of each atom in each layer (a_type <= MAX_TA)
  Integer  :: blurring     !d->      -  Type of instrumental broadening to simulate.
                           !            Choices are; PS_VGT, GAUSS, LORENZ
  Integer  :: bitdepth = 16    !i->  The bit-depth of the selected area diffraction.
                           !     pattern (sadp) data. Can equal 8 or 16.
  Integer  :: bckg_points  ! number of background points

  Integer  :: CENTRO = 1   !  numeric constant (= 1) indicating layer has a center of symmetry.

  Integer  :: cntrl        !  -  The device number to read input from.
                           !  If the file 'control.dif' is present, this file
                           !  becomes the default input. Otherwise, it is the keyboard.
  Integer    :: conv_a = 0 , conv_b=0 , conv_c=0 , conv_d = 0 , conv_e = 0, conv_f=0 , conv_g=0 !labels of refinable parameters  used not to repeat calculations
  Integer, dimension(MAX_TA) :: e_sf    !s->   Electron scattering factors.


  Integer  :: full_brd     !i->  1 if full adaptive integration is to be carried out on the sharp spots.

  Integer  :: full_shrp    !i->  1 if full adaptive integration is to be carried out on the streaks.

  Integer  ::   h_bnd      !  Maximum value of h to explore.
  Integer  ::   k_bnd      !  Maximum value of k to explore.

  Integer, dimension(MAX_L) :: l_actual    ! Contains the layer number that layer i is structurally identical to.
                                           ! If all layers are unique, l_actual(i) = i; else, l_actual(i) <= i

  Real(kind=sp)  ::  l_cnt        ! Number of layers in explicit sequence. This is
                                  ! tallied by DIFFaX, by counting the layers.
  Integer, dimension(MAX_L)  :: l_n_atoms   ! number of atoms in each layer

  Integer, dimension(MAX_L)  :: l_symmetry  !d->  symmetry of layer (NONE or CENTRO)

  Integer, dimension(XP_MAX) :: l_seq , l_seqdef  !d-> array containing the explicitly defined sequence of layers.
                                        !    Used only if 'xplcit' = TRUE

  Integer  :: loglin      !i->  1 if logarithmic scaling of SADP data is required ,  0: - Logarithmic  1: - Linear

  Integer    :: numcal = 0

  Integer  :: GAUSS  = 1   ! numeric constant (= 1) indicating user wishes to simulate Gaussian
                           ! instrumental broadening, with a constant half width

  Integer  :: LORENZ = 2   ! numeric constant (= 2) indicating user wishes to simulate Lorentzian
                           ! instrumental broadening, with a constant half width

  Integer  :: PS_VGT = 3   !  numeric constant (= 3) indicating user wishes to simulate pseudo-Voigt
                           !  instrumental broadening.

  Integer  :: PV_GSS = 4   !  numeric constant (= 4) indicating user wishes to simulate Gaussian
                           !  instrumental broadening, with a variable half width

  Integer  :: PV_LRN = 5   !  numeric constant (= 5) indicating user wishes to simulate Lorentzian
                           !  instrumental broadening, with a variable half width
  Integer  ::  maxsad    !  Maximum intensity of sadp patterns.

  Integer  ::  n_actual  ! Number of unique layers ( <= n_layers).
  Integer  ::  n_atoms   ! Temporary variable holding the number of unique atoms in a given layer.

  Integer  ::  n_layers  !d-> Number of user-defined layers. n_layers <= MAX_L

  Integer  ::  no_trials !  Number of reciprocal space points to sample.
                         !  whilst evaluating diffraction symmetry.
  Integer  ::   NONE = 0 !  numeric constant (= 0) indicating layer has no center of symmetry.

  Integer  ::  rad_type  !d->  Type of radiation for which to calculate diffraction intensities.
                         !     Choices are: X_RAY, NEUTRN, ELECTN

  Integer  ::  punts     !     number of excluded points
  Integer  ::  sadblock  !     Length of a row of SADP data
  Integer  ::  SymGrpNo  !d->  Numeric descriptor of diffraction symmetry.

  Integer  ::  X_RAY  = 0   !     Numeric constant (= 0) indicating radiation type is X-rays.
  Integer  ::  NEUTRN = 1   !     Numeric constant (= 1) indicating radiation type is neutrons.
  Integer  ::  ELECTN = 2   !     Numeric constant (= 2) indicating radiation type is electrons.
!
  Integer, dimension(max_bin) ::  pow         !?
  Integer                     ::  max_pow     ! ?

  Integer                       :: n_high

  Integer                       :: opt
  Integer :: d_punt     !number of reflections
  Integer :: h_min, h_max, k_min, k_max ! lower and higher values of index h and k

  Integer,          dimension(:),     allocatable  :: original
  Integer                                          :: fls, lls, otls, stls

  Integer :: i_plane !Plane in reciprocal space: 1: k = 0.   2: h = 0.   3: h = k.   4: h = -k
  Integer :: i_adapt = 1  !Adaptive quadrature always applied
  Integer :: funct_num ! Function number: 1=streak, 3=powder pattern, 4=SADP

!
!**********************     REAL(kind=dp) variables
!
  Real(kind=dp)  ::  a0    !  One of seven reciprocal lattice constants
  Real(kind=dp)  ::  ab0   !  One of seven reciprocal lattice constants
  Real(kind=dp)  ::  b0    !  One of seven reciprocal lattice constants
  Real(kind=dp)  ::  bc0   !  One of seven reciprocal lattice constants
  Real(kind=dp)  ::  c0    !  One of seven reciprocal lattice constants
  Real(kind=dp)  ::  ca0   !  One of seven reciprocal lattice constants
  Real(kind=dp)  ::  d0    !  One of seven reciprocal lattice constants
  Real(kind=dp)  :: cell_a      !d->  -  Unit cell a axis dimension.
  Real(kind=dp)  :: cell_b      !d->  -  Unit cell b axis dimension.
  Real(kind=dp)  :: cell_c      !d->  -  Unit cell c axis dimension.
  Real(kind=dp)  :: cell_gamma  !d->  -  Angle between a and b axes. Angle between b and
                         !        c, and a and c axes is 90 degrees by default.
  Real(kind=dp)  :: d_ret
  Real(kind=dp)  :: a_B11,a_B22,a_B33,a_B12,a_B23,a_B31 ! The average values of the r_Bij arrays

  Real(kind=dp)  :: bnds_wt    !  Equals 1.0 if rot_only is TRUE, otherwise equals 0.5
                        !  if rot_only is FALSE (ie there is a vertical mirror)
  Real(kind=dp)  :: brightness ! For SADP simulations

  Integer                                   :: num_streak
  Integer, dimension(:), allocatable        :: h_streak     ! These variables are used for calculating the intensity along the streak
  Integer, dimension(:), allocatable        :: k_streak     !
  Real,    dimension(:), allocatable        :: l0_streak    !
  Real,    dimension(:), allocatable        :: l1_streak    !
  Real(kind=dp), dimension(:), allocatable  :: dl_streak    !
  Integer, dimension(:), allocatable        :: streak_flags
  Integer        :: adapt_quad   ! Variable for adaptive quadrature

  Real(kind=dp)  :: d_theta    !d->  Angular increment in PXD spectrum.

  Real(kind=dp)  :: DEG2RAD    !     Conversion factor for degrees to radians

  Real(kind=dp)  :: fatsWalla_hk ! temporary storage for Fats-Waller factor
  Real(kind=dp)  :: ffact_scale  ! Angular scale (radians) of array 'formfactor'.
  Real(kind=dp)  :: ffhkcnst     ! Constant associated with the form-factor half-width.
                          ! Depends on reflection indices h and k, as well as Wx and Wy.
  Real(kind=dp)  :: ffwdth       ! Form factor half width in reciprocal Angstroms.

  Real(kind=dp)  ::  h_end    !- (h_end,k_end) is a vector in the h-k plane of reciprocal space defining the upper boundary
                       !  of the wedge in reciprocal space to integrate within. This is defined by the
                       !  symmetry of the diffraction data.
  Real(kind=dp)  ::  h_start  !- (h_start,k_start) is a vector in the h-k plane of reciprocal space defining the lower boundary
                       !  of the wedge in reciprocal space to integrate within. This is defined by the
                       !  symmetry of the diffraction data.
  Real(kind=dp)  ::  k_end    !- (h_end,k_end) is a vector in the h-k plane of reciprocal space defining the upper boundary
                       !  of the wedge in reciprocal space to integrate within. This is defined by the
                       !  symmetry of the diffraction data.
  Real(kind=dp)  ::  k_start  !- (h_start,k_start) is a vector in the h-k plane of reciprocal space defining the lower boundary
                       !  of the wedge in reciprocal space to integrate within. This is defined by the
                       !  symmetry of the diffraction data.
  Real(kind=dp)  ::  FWHM     !d->    Full width half maximum of instrumental broadening.

  Real(kind=dp)  :: l_bnd       !    Maximum value of l to explore.
  Real(kind=dp)  :: l_rz        !    Value of Rz if same_rz = TRUE.
  Real(kind=dp)  :: lambda, lambda2, ratio      !d-> Radiation wavelength.
  Real(kind=dp)  :: max_angle   !    Maximum angle that intensity information is to
                         !    be taken from for the purposes of evaluating
                         !    diffraction symmetry.
  Real(kind=dp)  :: max_var     !    Maximum mean variation of intensities when a
                         !    given symmetry operator was applied.
  Real(kind=dp)  :: mltplcty    !    1/mltplcty is the fraction of reciprocal space
                         !    necessary to integrate over, as determined by
                         !    the diffraction point group symmetry.

  Real(kind=sp), dimension(max_npar)    :: mult,  gen
  Real(kind=dp)  :: PI          !    The value of pi, 3.141592653589793.....
  Real(kind=dp)  :: PI2         !    The value of 2*pi
  Real(kind=sp)  :: pv_gamma    !d-> Pseudo-Voigt gamma parameter.
  Real(kind=dp)  :: pv_u        !d-> Pseudo-Voigt u parameter.
  Real(kind=dp)  :: pv_v        !d-> Pseudo-Voigt v parameter.
  Real(kind=dp)  :: pv_w        !d-> Pseudo-Voigt w parameter.
  Real(kind=sp)  :: pv_x, pv_dg, pv_dl   ! Pseudo-voigt x parameter and gaussian and lorentzian average volumetric sizes
  !Real(kind=sp)  :: pv_hg, pv_hl         ! gaussian and lorentzian FWHM
  Real(kind=dp)  :: RAD2DEG     !    Conversion factor for radians to degrees
  Real(kind=dp)  :: scaleint    !    Intensity scaling factor used in calculating
                                !    the selected area diffraction patterns.
  Real(kind=dp)  :: th2_max     !d-> Upper bound of angle in PXD spectrum in radians.
  Real(kind=dp)  :: th2_min     !d-> Lower bound of angle in PXD spectrum in radians.

  Real(kind=cp)  :: thmax=0.0    !  Upper bound of angle in PXD spectrum in degrees.
  Real(kind=cp)  :: thmin=0.0    !  Lower bound of angle in PXD spectrum in degrees.
  Real(kind=cp)  :: step_2th=0.0 !  Step in 2theta (in degrees) for calculating the PDP
  Integer        :: i_geom=0     ! =0 for Bragg-Brentano, =1 for Debye-Scherrer

  Real(kind=dp)  :: theta1      !    angle relative to (1,0,0) of lower wedge bound
  Real(kind=dp)  :: theta2      !    angle relative to (1,0,0) of upper wedge bound
  Real(kind=dp)  :: tolerance   !d-> Maximum deviation that intensities can have
                         !    from symmetry related points if intensities are
                         !    to be considered equal.
  Real(kind=dp)  :: tiny_inty   !    a small intensity value used in the diffraction
                         !    symmetry checking routines. Intensities lower
                         !    than tiny_inty are treated as being close to zero.
  Real(kind=dp)  ::  Wa         !i-> In-plane width of crystal along a-direction.
  Real(kind=dp)  ::  Wb         !i-> In-plane width of crystal perpendicular to
                         !    a-direction. Wx and Wy in Angstroms.

  Real(kind=sp)  :: rpo         !lowest rp
  Real(kind=sp)  :: chi2o       !lowest chi2

  Real(kind=dp)  :: l_upper         !upper limit of l, used in SADP simulations

  Real(kind=dp) ,dimension(max_bckg_points)   ::  bckg_p ,bckg_v    ! background position and background value

  Real(kind=dp), dimension(MAX_A,MAX_L)   :: a_B       !d-> isotropic Debye-Waller factor for each atom in each layer

  Real(kind=dp), dimension(MAX_A,MAX_L)   :: a_occup   !d-> Occupancy of each atom in each layer
                                                !    (Normally this will lie between 0 and 1)

  Real(kind=dp), dimension(3,MAX_A,MAX_L) :: a_pos     !d->  x,y,z relative coordinates of each atom in each layer.

  Real(kind=sp), dimension (:), allocatable     :: dos_theta
  Integer, dimension(:,:), allocatable          :: hkl_list
  Integer                                       :: n_hkl

  Real(kind=dp), dimension(3,MAX_A,MAX_L) :: l_r       !d->   Array of layer stacking vectors. The order is (column, row).

  Real(kind=dp), dimension(MAX_SP)        :: brd_spc   ! Array holding the powder diffraction data
                                                       ! after instrumental broadening has been applied.
  Real(kind=dp), dimension(MAX_SP)        :: spec      ! Array holding the unbroadened powder diffraction data.
                                                       ! This array also holds the SADP image data.
  Real,          dimension(max_sp)        :: strkAngl  ! This array holds values if 2theta angles for pv_streak calculationns

  Real(kind=dp), dimension(MAX_TA)        :: n_sf      !s->  Neutron scattering factors.

  Real(kind=sp), dimension (max_npar)     :: vector    !vector containing all the parameters to optimize

  Real(kind=dp), dimension(MAX_A,MAX_L)   :: detune
    ! detune -  Array of small positive numbers whose purpose is to prevent the determinant of the
    !           recursion array 'mat' from becoming zero at the sharp peaks. This produces a singularity which
    !           is hard to integrate accurately over. In essence, the 'detune' parameters are small stacking
    !           uncertainty factors. The result is to reduce the value of l_alpha(j,i) by an amount
    !           detune(j,i), such that the sum over the alphas for stacking from a given layer do not quite
    !           add to unity.

  Real(kind=dp), dimension(FFACT_SIZE)  ::   formfactor
    !formfactor-  Array containing a normalized Lorentzian profile,the form factor due to in-plane size broadening.
    !             The profile is Lorentzian out to N_SIGMAS half-widths, and linear to zero from there. The linear
    !             portion has the same gradient as the last point of the Lorentzian portion, thus the sampling step
    !             is governed by N_SIGMAS as well as FFACT_SIZE.



  Real(kind=dp), dimension(MAX_L)      :: high_atom ! The highest atomic z-rel position in each layer
  Real(kind=dp), dimension(MAX_L)      :: low_atom  ! The lowest atomic z-rel position in each layer

  Real(kind=dp), dimension(MAX_L)      :: l_g ! Array of layer existence probabilities.
                                       ! These are determined by the transition
                                       ! probabilities 'l_alpha' entered by the user.

  Real(kind=dp), dimension(MAX_A,MAX_L):: hx_ky     !  Temporary storage of h*Rx + k*Ry, whilst
                                             !  l*Rz is being computed along the streaks.

  Real(kind=dp), dimension(MAX_A,MAX_L):: l_alpha   !d-> Array of layer transition probabilities. The order is (column,row)

  Real(kind=dp), dimension(MAX_A,MAX_L):: r_B11,r_B22,r_B33,r_B12,r_B23,r_B13
           ! d-> The 6 components of the anisotropic layer stacking uncertainties. These are
           !     equivalent to the atomic Debye-Waller factors, except they apply to the stacking
           !     vectors. These parameters allow for 'turbostratic' disorder such as is found
           !     in liquid crystals. These parameters are optional, and can be entered by the user
           !     in the PARAMETERS section enclosed in parentheses.

  Real(kind=dp), dimension(9,MAX_TA) :: x_sf !s-> X-ray scattering factors.

  Real(kind=sp),dimension (max_sp)  :: ycalcdef

 !********************     complex*16 variables

   Complex(kind=dp), dimension(MAX_L,MAX_L) :: l_phi ! Phases of components of 'mat'

   Complex(kind=dp), dimension(MAX_L,MAX_L) :: mat   ! Recursion matrix relating the scattering
                                                     ! from crystals centered on different layers
   Complex(kind=dp), dimension(MAX_L,MAX_L) :: mat1  ! Storage for intermediate 'mat' results.

   Complex(kind=dp) :: wavefn   !  Coherent wavefunction calculated for an
                                !  explicitly defined sequence of layers (if requested)

  Integer, parameter        :: max_excl=100 !Maximum number of excluded regions
  Integer                   :: nexcrg       !Number of excluded regions
  Real, dimension(max_excl) :: alow, ahigh  !Excluded regions


  integer, parameter        :: max_avercell=25 !max number of transition vectors needed for calculating the average cell

  Contains

    Subroutine Close_Faults()
      character(len=1) :: keyw
      write(unit=*,fmt="(/,a)") " => Press <Enter> to finish ..."
      read(unit=*,fmt="(a)") keyw
      stop
    End Subroutine Close_Faults

 End Module diffax_mod
!----------------------------------------------------------------------------------------------------
