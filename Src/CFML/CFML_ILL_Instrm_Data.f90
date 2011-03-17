!!----
!!---- Copyleft(C) 1999-2011,              Version: 5.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: CFML_ILL_Instrm_Data
!!----   INFO: Subroutines related to Instrument information from ILL
!!----
!!---- HISTORY
!!----    Update: 05/03/2011
!!----
!!--..    The default instrument cartesian frame is that defined by
!!--..       W.R.Busing & H.A.Levy (Acta Cryst. 22,457-464 (1967)
!!--..    The reciprocal Cartesian Frame used by Busing-Levy has "x" along a*
!!--..    "y" in the (a*,b*) plane and "z" perpendicular to that plane. The
!!--..    matrix passing from (h,k,l) in r.l.u. to the cartesian frames is
!!--..             | a*   b*.cos(gamma*)   c*.cos(beta*)            |
!!--..      [B] =  | 0    b*.sin(gamma*)  -c*.sin(beta*).cos(alpha) |
!!--..             ! 0          0          1/c                      |
!!--..
!!--..          [hc] = [B] [h]
!!--..
!!--..      Let U be the orthogonal matrix relating the phi-axis system
!!--..      attached to the phi-shaft of the instrument with the Cartesian crystal system
!!--..      then     [h-phi]= [U] [hc]      [h-phi] = [U] [B] [h] = [UB] [h]
!!--..
!!--..      Let us define three more cartesian systems attached to the Chi,Omega and
!!--..      Theta axes, coincident with phi-axis when all instrument angles are set
!!--..      to zero :
!!--..         [h-chi] = [PHI] [h-phi]
!!--..         [h-omg] = [CHI] [h-chi]
!!--..         [h-tet] = [OMG] [h-omg]
!!--..
!!--..               |  cos(phi)      sin(phi)      0  |             |zL-axis
!!--..       [PHI] = | -sin(phi)      cos(phi)      0  |             |
!!--..               |     0            0           1  |             |
!!--..                                             incoming beam ->__|_________yL-axis
!!--..               |  cos(chi)       0     sin(chi)  |            /
!!--..       [CHI] = |     0           1       0       |           /
!!--..               | -sin(chi)       0     cos(chi)  |          /xL-axis
!!--..
!!--..               |  cos(omg)      sin(omg)      0  |
!!--..       [OMG] = | -sin(omg)      cos(omg)      0  |
!!--..               |     0            0           1  |
!!--..
!!--..      Let us define 2 more Cartesian frames: The laboratory fixed system
!!--..      and the 2Theta-axis system attached to the 2theta shaft. All coincide
!!--..      for all angles equal to zero.
!!--..       [h-Lab] = [THE]  [h-tet] = [N] [h-omg]
!!--..       [h-tte] = [THE]t [h-tet] = [M] [h-omg]
!!--..     Where:
!!--..               |  cos(the)      sin(the)      0  |
!!--..       [THE] = | -sin(the)      cos(the)      0  |
!!--..               |     0            0           1  |
!!--..
!!--..                           |  cos(nu)      sin(nu)      0  |
!!--..       [N] = [THE] [OMG] = | -sin(nu)      cos(nu)      0  |  with nu=omg+tet
!!--..                           |     0            0         1  |
!!--..
!!--..                            |  cos(mu)      sin(mu)      0  |
!!--..       [M] = [THE]t [OMG] = | -sin(mu)      cos(mu)      0  |  with mu=omg-tet
!!--..                            |     0            0         1  |
!!--..    All angles, except Chi, are left-handed rotations about their respective axes.
!!--..
!!--..      Basic Diffractometer equations:
!!--..            q = |h| = sqrt(dot_product(hc,hc)), with [hc]=[B][h]
!!--..      Bragg equation: sin(theta) = Lambda.q/2
!!--..      The reflection [h] is in the diffraction position if
!!--..           [h-tet] = [OMG] [CHI] [PHI] [U] [B] [h]
!!--..      has the form [h-tet]t = (q, 0, 0)
!!--..
!!--..    The orientation matrix U can be obtained from two non-colinear reflections h1 & h2,
!!--..    provided the cell parameters are known
!!--..
!!--..    Two kind of Busing-Levy (BL) frames are possible:
!!--..    For both frames the origin is at sample position and the y-axis
!!--..    positive sense is along the secondary beam (from monochromator to sample)
!!--..    a)  z upward
!!--..    b)  z downward
!!--..    The x-axis makes a right-handed frame with the other axes. Positive 2theta
!!--..    angles are from y-axis towards x-axis
!!--..    A change of geometry is always possible by providing the components of the
!!--..    axes {e1,e2,e3} with respect to the standard {i,j,k} BL-frame
!!--..
!!--..    igeom=1: Bissectrice (PSI=0)
!!--..    igeom=2: Bissecting - HiCHI
!!--..    igeom=3: Normal beam
!!--..    igeom=4: Parallel (PSI=90)
!!----
!!---- DEPENDENCIES
!!--++   Use CFML_GlobalDeps,       only: sp, dp, cp, pi, to_deg, to_rad, eps, OPS, OPS_sep
!!--++   Use CFML_Math_General,     only: cosd,sind
!!--++   use CFML_String_Utilities, only: u_case, lcase, Get_LogUnit, Number_Lines
!!--++   use CFML_Math_3D,          only: err_math3d,err_math3d_mess, Cross_Product, Determ_A, Determ_V, &
!!--++                                    invert => Invert_A
!!----
!!----
!!---- VARIABLES
!!--..    Types
!!----    BASIC_NUMORC_TYPE
!!----    BASIC_NUMORI_TYPE
!!----    BASIC_NUMORR_TYPE
!!----    CALIBRATION_DETECTOR_TYPE
!!----    DIFFRACTOMETER_TYPE
!!----    GENERIC_NUMOR_TYPE
!!----    ILL_DATA_RECORD_TYPE
!!----    POWDER_NUMOR_TYPE
!!----    SXTAL_NUMOR_TYPE
!!----    SXTAL_ORIENT_TYPE
!!--..
!!----    CURRENT_INSTRM
!!--++    CURRENT_INSTRM_SET                [Private]
!!----    CURRENT_ORIENT
!!----    CYCLE_NUMBER
!!----    ERR_ILLDATA
!!----    ERR_ILLDATA_MESS
!!--++    GOT_ILL_DATA_DIRECTORY            [Private]
!!--++    GOT_ILL_TEMP_DIRECTORY            [Private]
!!----    ILL_DATA_DIRECTORY
!!----    ILL_TEMP_DIRECTORY
!!----    INSTRM_DIRECTORY
!!--++    INSTRM_DIRECTORY_SET              [Private]
!!--++    IVALUES                           [Private]
!!----    MACHINE_NAME
!!--++    N_KEYTYPES                        [Private]
!!--++    NL_KEYTYPES                       [Private]
!!--++    NTEXT                             [Private]
!!--++    NVAL_F                            [Private]
!!--++    NVAL_I                            [Private]
!!--++    RVALUES                           [Private]
!!--++    TEXT_ILL                          [Private]
!!--++    UNCOMPRESSCOMMAND                 [Private]
!!----    YEAR_ILLDATA
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ALLOCATE_NUMORS
!!--++       ALLOCATE_POWDER_NUMORS          [Overloaded]
!!--++       ALLOCATE_SXTAL_NUMORS           [Overloaded]
!!----       DEFINE_UNCOMPRESS_PROGRAM
!!----       GET_ABSOLUTE_DATA_PATH
!!----       GET_NEXT_YEARCYCLE
!!----       GET_SINGLE_FRAME_2D
!!----       INITIALIZE_DATA_DIRECTORY
!!--++       INITIALIZE_NUMORS_DIRECTORY     [Private]
!!--++       INITIALIZE_TEMP_DIRECTORY       [Private]
!!--++       NUMBER_KEYTYPES_ON_FILE         [Private]
!!----       POWDERNUMOR_TO_DIFFPATTERN
!!--++       READ_A_KEYTYPE                  [Private]
!!----       READ_CURRENT_INSTRM
!!--++       READ_F_KEYTYPE                  [Private]
!!--++       READ_I_KEYTYPE                  [Private]
!!--++       READ_J_KEYTYPE                  [Private]
!!----       READ_NUMOR
!!--++       READ_NUMOR_D19                  [Private]
!!--++       READ_NUMOR_D1B                  [Private]
!!--++       READ_NUMOR_D20                  [Private]
!!--++       READ_NUMOR_D2B                  [Private]
!!--++       READ_NUMOR_D9                   [Private]
!!--++       READ_POWDER_NUMOR               [Overloaded]
!!--++       READ_R_KEYTYPE                  [Private]
!!--++       READ_S_KEYTYPE                  [Private]
!!--++       READ_V_KEYTYPE                  [Private]
!!--++       READ_SXTAL_NUMOR                [Overloaded]
!!----       SET_CURRENT_ORIENT
!!----       SET_DEFAULT_INSTRUMENT
!!----       SET_ILL_DATA_DIRECTORY
!!----       SET_INSTRM_DIRECTORY
!!----       SET_KEYTYPES_ON_FILE            [Private]
!!----       UPDATE_CURRENT_INSTRM_UB
!!----       WRITE_CURRENT_INSTRM_DATA
!!----       WRITE_HEADERINFO_NUMOR
!!--++       WRITE_HEADERINFO_POWDER_NUMOR   [Overloaded]
!!--++       WRITE_HEADERINFO_SXTAL_NUMOR    [Overloaded]
!!----       WRITE_NUMOR_INFO
!!--++       WRITE_POWDER_NUMOR              [Overloaded]
!!--++       WRITE_SXTAL_NUMOR               [Overloaded]
!!----
!!
Module CFML_ILL_Instrm_Data
   !---- Use Modules ----!
   !use f2kcli !Comment for compliant F2003 compilers
   Use CFML_GlobalDeps
   Use CFML_Math_General,         only: cosd,sind, equal_vector
   use CFML_String_Utilities,     only: u_case, lcase, Get_LogUnit, Number_Lines, Reading_Lines
   use CFML_Math_3D,              only: err_math3d,err_math3d_mess, Cross_Product, Determ_A, Determ_V, &
                                        invert => Invert_A
   use CFML_Diffraction_Patterns, only: Diffraction_Pattern_Type, Allocate_Diffraction_Pattern

   !---- Variables ----!
   Implicit none

   private

   !---- Public Subroutines ----!
   public :: Set_Current_Orient, Read_Numor, Read_Current_Instrm, Write_Current_Instrm_data,     &
             Allocate_Numors, Set_ILL_data_directory, Set_Instrm_directory,                      &
             Update_Current_Instrm_UB, Set_Default_Instrument,Get_Single_Frame_2D,               &
             Initialize_Data_Directory, Get_Absolute_Data_Path, Get_Next_YearCycle,              &
             Write_Generic_Numor, Set_Instrm_Geometry_Directory, Write_Numor_Info,               &
             Define_Uncompress_Program, PowderNumor_To_DiffPattern, Write_HeaderInfo_Numor,      &
             Read_Calibration_File_D1A

   !---- Private Subroutines ----!
   private:: Initialize_Numors_Directory,Initialize_Temp_Directory,Number_Keytypes_On_File, &
             Read_A_Keytype,Read_F_Keytype,Read_I_Keytype,Read_J_Keytype,Read_R_Keytype,    &
             Read_S_Keytype,Read_V_Keytype,Set_Keytypes_On_File, Read_Powder_Numor,         &
             Read_SXTAL_Numor, Read_Numor_Generic, Read_Numor_D1B, Read_Numor_D20,          &
             Read_Numor_D9, Read_Numor_D19, Write_POWDER_Numor, Write_SXTAL_Numor,          &
             Write_HeaderInfo_POWDER_Numor, Write_HeaderInfo_SXTAL_Numor, Read_Numor_D2B,   &
             Allocate_SXTAL_numors, Allocate_Powder_Numors, Read_Numor_D1A


   !---- Definitions ----!

   !!----
   !!---- TYPE :: BASIC_NUMC_TYPE
   !!--..
   !!---- Type, public :: Basic_NumC_Type
   !!----    integer                                      :: N        ! Number of elements in this Type
   !!----    character(len=40), dimension(:), allocatable :: NameVar  ! Name of the diferents fields
   !!----    character(len=80), dimension(:), allocatable :: CValues  ! Character Values
   !!---- End Type Basic_NumC_Type
   !!----
   !!---- Update: April - 2009
   !!
   Type, public :: Basic_NumC_Type
      integer                                      :: N
      character(len=40), dimension(:), allocatable :: NameVar
      character(len=80), dimension(:), allocatable :: CValues
   End Type Basic_NumC_Type

   !!----
   !!---- TYPE :: BASIC_NUMI_TYPE
   !!--..
   !!---- Type, public :: Basic_NumI_Type
   !!----    integer                                      :: N        ! Number of elements in this Type
   !!----    character(len=40), dimension(:), allocatable :: NameVar  ! Name of the diferents fields
   !!----    integer,           dimension(:), allocatable :: IValues  ! Integer values
   !!---- End Type, private :: Basic_NumI_Type
   !!----
   !!---- Update: April - 2009
   !!
   Type, public :: Basic_NumI_Type
      integer                                      :: N
      character(len=40), dimension(:), allocatable :: NameVar
      integer,           dimension(:), allocatable :: IValues
   End Type Basic_NumI_Type

   !!----
   !!---- TYPE :: BASIC_NUMR_TYPE
   !!--..
   !!---- Type, public :: Basic_NumR_Type
   !!----    integer                                      :: N        ! Number of elements in this Type
   !!----    character(len=40), dimension(:), allocatable :: NameVar  ! Name of the diferents fields
   !!----    real(kind=cp),     dimension(:), allocatable :: RValues  ! Real Values
   !!---- End Type Basic_NumR_Type
   !!----
   !!---- Update: April - 2009
   !!
   Type, public :: Basic_NumR_Type
      integer                                      :: N
      character(len=40), dimension(:), allocatable :: NameVar
      real(kind=cp),     dimension(:), allocatable :: RValues
   End Type Basic_NumR_Type

   !!----
   !!---- TYPE :: CALIBRATION_DETECTOR_TYPE
   !!--..
   !!---- Type, public :: Calibration_Detector_Type
   !!----   character(len=4)                             :: Name_Instrm       ! Instrument Name
   !!----   integer                                      :: NDet              ! Number of Detectors
   !!----   integer                                      :: NPointsDet        ! Number of Points by Detector
   !!----   logical, dimension(:), allocatable           :: Active            ! Flag for active detector or not
   !!----   real(kind=cp), dimension(:), allocatable     :: PosX              ! Relative Positions of each Detector
   !!----   real(kind=cp), dimension(:,:), allocatable   :: Effic             ! Efficiency of each point detector (NpointsDetect,NDect)
   !!---- End Type Calibration_Detector_Type
   !!----
   !!---- Update: April - 2009
   !!
   Type, public :: Calibration_Detector_Type
      character(len=4)                             :: Name_Instrm       ! Instrument Name
      integer                                      :: NDet             ! Number of Detectors
      integer                                      :: NPointsDet       ! Number of Points by Detector
      logical, dimension(:), allocatable           :: Active            ! Flag for active detector or not
      real(kind=cp), dimension(:), allocatable     :: PosX               ! Relative Positions of each Detector
      real(kind=cp), dimension(:,:), allocatable   :: Effic             ! Efficiency of each point detector (NpointsDetect,NDect)
   End Type Calibration_Detector_Type

   !!----
   !!---- TYPE :: DIFFRACTOMETER_TYPE
   !!--..
   !!---- Type, public :: diffractometer_type
   !!----    character(len=80)                         :: info                 !information about the instrument
   !!----    character(len=10)                         :: name_inst            !Short name of the instrument
   !!----    character(len=15)                         :: geom                 !"Eulerian_4C","Kappa_4C","Lifting_arm","Powder","Laue"
   !!----    character(len=6)                          :: BL_frame             !Kind of BL-frame: "z-up" or "z-down"
   !!----    character(len=4)                          :: dist_units           !distance units: "mm  ","cm  ","inch"
   !!----    character(len=4)                          :: angl_units           !angle units: "rad","deg"
   !!----    character(len=30)                         :: detector_type        !"Point","Flat_rect","Cylin_ImPlate","Tube_PSD", Put ipsd=1,2,...
   !!----    real(kind=cp)                             :: dist_samp_detector   ! dist. to centre for: point, Flat_rect, Tube_PSD; radius for: Cylin_ImPlate
   !!----    real(kind=cp)                             :: wave_min,wave_max    !Minimum and maximum wavelengths (Laue diffractometers)
   !!----    real(kind=cp)                             :: vert                 !Vertical dimension
   !!----    real(kind=cp)                             :: horiz                !Horizontal dimension
   !!----    real(kind=cp)                             :: agap                 !gap between anodes
   !!----    real(kind=cp)                             :: cgap                 !gap between cathodes
   !!----    integer                                   :: np_vert              !number of pixels in vertical direction
   !!----    integer                                   :: np_horiz             !number of pixels in horizontal direction
   !!----    integer                                   :: igeom                !1: Bissectrice (PSI=0),2: Bissecting - HiCHI, 3: Normal beam, 4:Parallel (PSI=90)
   !!----    integer                                   :: ipsd                 !1: Flat,2: Vertically Curved detector (used in D19amd)
   !!----    real(kind=cp),dimension(3)                :: e1                   !Components of e1 in {i,j,k}
   !!----    real(kind=cp),dimension(3)                :: e2                   !Components of e2 in {i,j,k}
   !!----    real(kind=cp),dimension(3)                :: e3                   !Components of e3 in {i,j,k}
   !!----    integer                                   :: num_ang              !Number of angular motors
   !!----    character(len=12),dimension(15)           :: ang_names            !Name of angular motors
   !!----    real(kind=cp),dimension(15,2)             :: ang_limits           !Angular limits (up to 15 angular motors)
   !!----    real(kind=cp),dimension(15)               :: ang_offsets          !Angular offsets
   !!----    integer                                   :: num_disp             !Number of displacement motors
   !!----    character(len=12),dimension(10)           :: disp_names           !Name of displacement motors
   !!----    real(kind=cp),dimension(10,2)             :: disp_limits          !Displacement limits (up to 15 displacement motors)
   !!----    real(kind=cp),dimension(10)               :: disp_offsets         !Displacement offsets
   !!----    real(kind=cp),dimension(3 )               :: det_offsets          !Offsets X,Y,Z of the detector centre
   !!----    real(kind=cp),dimension(:,:), allocatable :: alphas               !Efficiency corrections for each pixel
   !!---- End Type diffractometer_type
   !!----
   !!----    Definition for Diffractometer type
   !!----
   !!---- Update: April - 2008
   !!
   Type, public :: diffractometer_type
      character(len=80)                         :: info                 !information about the instrument
      character(len=10)                         :: name_inst            !Short name of the instrument
      character(len=15)                         :: geom                 !"Eulerian_4C","Kappa_4C","Lifting_arm","Powder","Laue"
      character(len=6)                          :: BL_frame             !Kind of BL-frame: "z-up" or "z-down"
      character(len=4)                          :: dist_units           !distance units: "mm  ","cm  ","inch"
      character(len=4)                          :: angl_units           !angle units: "rad","deg"
      character(len=30)                         :: detector_type        !"Point","Flat_rect","Cylin_ImPlate","Tube_PSD", Put ipsd=1,2,...
      real(kind=cp)                             :: dist_samp_detector   ! dist. to centre for: point, Flat_rect, Tube_PSD; radius for: Cylin_ImPlate
      real(kind=cp)                             :: wave_min,wave_max    !Minimum and maximum wavelengths (Laue diffractometers)
      real(kind=cp)                             :: vert                 !Vertical dimension
      real(kind=cp)                             :: horiz                !Horizontal dimension
      real(kind=cp)                             :: agap                 !gap between anodes
      real(kind=cp)                             :: cgap                 !gap between cathodes
      integer                                   :: np_vert              !number of pixels in vertical direction
      integer                                   :: np_horiz             !number of pixels in horizontal direction
      integer                                   :: igeom                !1: Bissectrice (PSI=0),2: Bissecting - HiCHI, 3: Normal beam, 4:Parallel (PSI=90)
      integer                                   :: ipsd                 !1: Flat,2: Vertically Curved detector (used in D19amd)
      real(kind=cp),dimension(3)                :: e1                   !Components of e1 in {i,j,k}
      real(kind=cp),dimension(3)                :: e2                   !Components of e2 in {i,j,k}
      real(kind=cp),dimension(3)                :: e3                   !Components of e3 in {i,j,k}
      integer                                   :: num_ang              !Number of angular motors
      character(len=12),dimension(15)           :: ang_names            !Name of angular motors
      real(kind=cp),dimension(15,2)             :: ang_limits           !Angular limits (up to 15 angular motors)
      real(kind=cp),dimension(15)               :: ang_offsets          !Angular offsets
      integer                                   :: num_disp             !Number of displacement motors
      character(len=12),dimension(10)           :: disp_names           !Name of displacement motors
      real(kind=cp),dimension(10,2)             :: disp_limits          !Displacement limits (up to 15 displacement motors)
      real(kind=cp),dimension(10)               :: disp_offsets         !Displacement offsets
      real(kind=cp),dimension(3 )               :: det_offsets          !Offsets X,Y,Z of the detector centre
      real(kind=cp),dimension(:,:), allocatable :: alphas               !Efficiency corrections for each pixel
   End Type diffractometer_type

   !!----
   !!---- TYPE :: GENERIC_NUMOR_TYPE
   !!--..
   !!---- Type, public :: Generic_Numor_type
   !!----    integer                                    :: Numor       ! Numor
   !!----    character(len=4)                           :: Instr       ! Instrument on ILL
   !!----    character(len=10)                          :: ExpName     ! Experimental Name
   !!----    character(len=20)                          :: Date        ! Date
   !!----    character(len=80)                          :: Title       ! Title
   !!----    type(basic_numc_type)                      :: SampleID    ! Sample Identification
   !!----    type(basic_numr_type)                      :: DiffOpt     ! Diffractometer Optics and Reactor Parameters
   !!----    type(basic_numr_type)                      :: MonMPar     ! Monochromator Motor Parameters
   !!----    type(basic_numr_type)                      :: DiffMPar    ! Diffractometer Motor Parameters
   !!----    type(basic_numr_type)                      :: DetPar      ! Detector Parameters
   !!----    type(basic_numi_type)                      :: DACFlags    ! Data Acquisition Control
   !!----    type(basic_numr_type)                      :: DACParam    ! Data Acquisition Parameters
   !!----    type(basic_numr_type)                      :: SampleSt    ! Sample status
   !!----    type(basic_numi_type)                      :: ICounts     ! Counts as Integers
   !!----    type(basic_numr_type)                      :: RCounts     ! Counts as Reals
   !!---- End Type Generic_Numor_Type
   !!----
   !!----    Definition for Generic Numor type
   !!----
   !!---- Update: April - 2009
   !!
   Type, public :: Generic_Numor_type
      integer                                    :: Numor       ! Numor
      character(len=4)                           :: Instr       ! Instrument on ILL
      character(len=10)                          :: ExpName     ! Experimental Name
      character(len=20)                          :: Date        ! Date
      character(len=80)                          :: Title       ! Title
      type(basic_numc_type)                      :: SampleID    ! Sample Identification
      type(basic_numr_type)                      :: DiffOpt     ! Diffractometer Optics and Reactor Parameters
      type(basic_numr_type)                      :: MonMPar     ! Monochromator Motor Parameters
      type(basic_numr_type)                      :: DiffMPar    ! Diffractometer Motor Parameters
      type(basic_numr_type)                      :: DetPar      ! Detector Parameters
      type(basic_numi_type)                      :: DACFlags    ! Data Acquisition Control
      type(basic_numr_type)                      :: DACParam    ! Data Acquisition Parameters
      type(basic_numr_type)                      :: SampleSt    ! Sample status
      type(basic_numi_type)                      :: ICounts     ! Counts as Integers
      type(basic_numr_type)                      :: RCounts     ! Counts as Reals
   End Type Generic_Numor_Type

   !!----
   !!---- TYPE :: ILL_DATA_RECORD_TYPE
   !!--..
   !!---- Type, public :: ILL_data_record_type
   !!----    integer                     :: numor      ! data set numor.
   !!----    integer                     :: nset_prime ! set number (groups of 100000 numor).
   !!----    integer                     :: ntran      ! (key2) 0 or numcomp => data transferred?
   !!----    character(len=4)            :: inst_ch    ! instrument name (4 characters)
   !!----    character(len=22)           :: date_ch    ! measurement date (22 characters). !it was 18
   !!----    character(len=2)            :: fill_ch    ! 2 characters (key3) leader
   !!----    character(len=6)            :: user_ch    ! user name (6 characters)
   !!----    character(len=6)            :: lc_ch      ! local contact name (6 characters)
   !!----    character(len=72)           :: text_ch    ! commentary (72characters)
   !!----    character(len=8)            :: scan_motor ! principal scan motor name. (8 characters)
   !!----    integer                     :: nvers      !ival(1),  data version number
   !!----    integer                     :: ntype      !ival(2),  data type - single/multi/powder
   !!----    integer                     :: kctrl      !ival(3),  data function type
   !!----    integer                     :: manip      !ival(4),  principle scan angle
   !!----    integer                     :: nbang      !ival(5),  number of data saved
   !!----    integer                     :: nkmes      !ival(6),  pre-calculated number of points
   !!----    integer                     :: npdone     !ival(7),  actual number of points
   !!----    integer                     :: jcode      !ival(8),  count on monitor/time
   !!----    integer                     :: icalc      !ival(9),  angle calculation type
   !!----    integer                     :: ianal      !ival(10), analyser present (d10)
   !!----    integer                     :: imode      !ival(11), 2th motor sense (d10)
   !!----    integer                     :: itgv       !ival(12), d19/d9 fast measurement
   !!----    integer                     :: iregul     !ival(13), temperature monitor function
   !!----    integer                     :: ivolt      !ival(14), voltmeter function
   !!----    integer                     :: naxe       !ival(15), d10 (number of axes)
   !!----    integer                     :: npstart    !ival(16), point starting no frag. numor (d19/16)
   !!----    integer                     :: ilasti     !ival(17), elastic measurement (d10)
   !!----    integer                     :: isa        !ival(18), analyser motor sense (d10)
   !!----    integer                     :: flgkif     !ival(19), constant ki or kf (d10)
   !!----    integer                     :: ih_sqs     !ival(20), d10 sqs variation on h
   !!----    integer                     :: ik_sqs     !ival(21), d10 sqs variation on k
   !!----    integer                     :: nbsqs      !ival(22), d10 sqs slice number
   !!----    integer                     :: nb_cells   !ival(24), multi/powder data - number of detectors
   !!----    integer                     :: nfree1     !          data control (free).
   !!----    integer,dimension(11)       :: icdesc     !
   !!----    real(kind=cp), dimension(35):: valco      !rval( 1:35)
   !!----    real(kind=cp), dimension(10):: valdef     !rval(36:45)
   !!----    real(kind=cp), dimension(5) :: valenv     !rval(46:50)
   !!---- End Type ILL_data_record_type
   !!----
   !!----    Definition for Data Record type
   !!----
   !!---- Update: April - 2008
   !!
   Type, public :: ILL_data_record_type
      integer                     :: numor      ! data set numor.
      integer                     :: nset_prime ! set number (groups of 100000 numor).
      integer                     :: ntran      ! (key2) 0 or numcomp => data transferred?
      character(len=4)            :: inst_ch    ! instrument name (4 characters)
      character(len=22)           :: date_ch    ! measurement date (22 characters). !it was 18
      character(len=2)            :: fill_ch    ! 2 characters (key3) leader
      character(len=6)            :: user_ch    ! user name (6 characters)
      character(len=6)            :: lc_ch      ! local contact name (6 characters)
      character(len=72)           :: text_ch    ! commentary (72characters)
      character(len=8)            :: scan_motor ! principal scan motor name. (8 characters)
      integer                     :: nvers      !ival(1),  data version number
      integer                     :: ntype      !ival(2),  data type - single/multi/powder
      integer                     :: kctrl      !ival(3),  data function type
      integer                     :: manip      !ival(4),  principle scan angle
      integer                     :: nbang      !ival(5),  number of data saved
      integer                     :: nkmes      !ival(6),  pre-calculated number of points
      integer                     :: npdone     !ival(7),  actual number of points
      integer                     :: jcode      !ival(8),  count on monitor/time
      integer                     :: icalc      !ival(9),  angle calculation type
      integer                     :: ianal      !ival(10), analyser present (d10)
      integer                     :: imode      !ival(11), 2th motor sense (d10)
      integer                     :: itgv       !ival(12), d19/d9 fast measurement
      integer                     :: iregul     !ival(13), temperature monitor function
      integer                     :: ivolt      !ival(14), voltmeter function
      integer                     :: naxe       !ival(15), d10 (number of axes)
      integer                     :: npstart    !ival(16), point starting no frag. numor (d19/16)
      integer                     :: ilasti     !ival(17), elastic measurement (d10)
      integer                     :: isa        !ival(18), analyser motor sense (d10)
      integer                     :: flgkif     !ival(19), constant ki or kf (d10)
      integer                     :: ih_sqs     !ival(20), d10 sqs variation on h
      integer                     :: ik_sqs     !ival(21), d10 sqs variation on k
      integer                     :: nbsqs      !ival(22), d10 sqs slice number
      integer                     :: nb_cells   !ival(24), multi/powder data - number of detectors
      integer                     :: nfree1     !          data control (free).
      integer,dimension(11)       :: icdesc     !
      real(kind=cp), dimension(35):: valco      !rval( 1:35)
      real(kind=cp), dimension(10):: valdef     !rval(36:45)
      real(kind=cp), dimension(5) :: valenv     !rval(46:50)
   End Type ILL_data_record_type

   !!----
   !!---- TYPE :: POWDER_NUMOR_TYPE
   !!--..
   !!---- Type, public :: POWDER_Numor_type
   !!----    integer                                    :: numor       ! Numor
   !!----    integer                                    :: manip       ! principle scan angle
   !!----    integer                                    :: icalc       ! angle calculation type
   !!----    character(len=32)                          :: header      ! User, local contact, date
   !!----    character(len=4)                           :: Instrm      ! Instrument name
   !!----    character(len=32)                          :: title       !
   !!----    character(len=8)                           :: Scantype    ! omega, phi, etc...
   !!----    real(kind=cp), dimension(5)                :: angles      ! Angles: phi, chi, omega, 2theta(gamma), psi
   !!----    real(kind=cp), dimension(3)                :: scans       ! scan start, scan step, scan width
   !!----    real(kind=cp)                              :: monitor     !
   !!----    real(kind=cp)                              :: time        !
   !!----    real(kind=cp)                              :: wave        ! wavelength
   !!----    real(kind=cp), dimension(5)                :: conditions  ! Temp-s.pt,Temp-Regul,Temp-sample,Voltmeter,Mag.field
   !!----    integer                                    :: nbdata      ! Total number of pixels nx*ny = np_vert*np_horiz
   !!----    integer                                    :: nframes     ! Total number of frames
   !!----    integer                                    :: nbang       ! Total number of angles moved during scan
   !!----    integer, dimension(11)                     :: icdesc      ! Integer values
   !!----    real(kind=cp),allocatable,dimension(:,:)   :: tmc_ang     ! time,monitor,total counts, angles*1000
   !!----                                                              ! To be allocated as tmc_ang(nbang,nframes)
   !!----    real(kind=cp),allocatable,dimension(:,:)   :: counts      ! Counts array to be reshaped (np_vert,np_horiz,nframes) in case of 2D detectors
   !!----                                                              ! To be allocated as counts(nbdata,nframes)
   !!---- End Type POWDER_Numor_type
   !!----
   !!----    Definition for POWDER Numor type
   !!----
   !!---- Update: April - 2009
   !!
   Type, public :: POWDER_Numor_type
      integer                                    :: numor       ! Numor
      integer                                    :: manip       ! principle scan angle
      integer                                    :: icalc       ! angle calculation type
      character(len=32)                          :: header      ! User, local contact, date
      character(len=4)                           :: Instrm      ! Instrument name
      character(len=32)                          :: title       !
      character(len=8)                           :: Scantype    ! omega, phi, etc...
      real(kind=cp), dimension(5)                :: angles      ! Angles: phi, chi, omega, 2theta(gamma), psi
      real(kind=cp), dimension(3)                :: scans       ! scan start, scan step, scan width
      real(kind=cp)                              :: monitor     !
      real(kind=cp)                              :: time        !
      real(kind=cp)                              :: wave        ! wavelength
      real(kind=cp), dimension(5)                :: conditions  ! Temp-s.pt,Temp-Regul,Temp-sample,Voltmeter,Mag.field
      integer                                    :: nbdata      ! Total number of pixels nx*ny = np_vert*np_horiz
      integer                                    :: nframes     ! Total number of frames
      integer                                    :: nbang       ! Total number of angles moved during scan
      integer, dimension(11)                     :: icdesc      ! Integer values
      real(kind=cp),allocatable,dimension(:,:)   :: tmc_ang     ! time,monitor,total counts, angles*1000
                                                                ! To be allocated as tmc_ang(nbang,nframes)
      real(kind=cp),allocatable,dimension(:,:)   :: counts      ! Counts array to be reshaped (np_vert,np_horiz,nframes) in case of 2D detectors
                                                                ! To be allocated as counts(nbdata,nframes)
   End type POWDER_Numor_type

   !!----
   !!---- TYPE :: SXTAL_NUMOR_TYPE
   !!--..
   !!---- Type, public :: SXTAL_Numor_type
   !!----    integer                                    :: numor       ! Numor
   !!----    integer                                    :: manip       ! principle scan angle
   !!----    integer                                    :: icalc       ! angle calculation type
   !!----    character(len=32)                          :: header      ! User, local contact, date
   !!----    character(len=4)                           :: Instrm      ! Instrument name
   !!----    character(len=32)                          :: title       !
   !!----    character(len=8)                           :: Scantype    ! omega, phi, etc...
   !!----    real(kind=cp), dimension(3)                :: hmin        ! or h,k,l for omega-scans
   !!----    real(kind=cp), dimension(3)                :: hmax        !
   !!----    real(kind=cp), dimension(5)                :: angles      ! Angles: phi, chi, omega, 2theta(gamma), psi
   !!----    real(kind=cp), dimension(3,3)              :: UB          ! UB-matrix
   !!----    real(kind=cp), dimension(3)                :: dh          ! delta_h, delta_k, delta_l
   !!----    real(kind=cp), dimension(3)                :: scans       ! scan start, scan step, scan width
   !!----    real(kind=cp)                              :: preset      !
   !!----    real(kind=cp)                              :: wave        ! wavelength
   !!----    real(kind=cp)                              :: cpl_fact    ! Coupling Factor
   !!----    real(kind=cp), dimension(5)                :: conditions  ! Temp-s.pt,Temp-Regul,Temp-sample,Voltmeter,Mag.field
   !!----    integer                                    :: nbdata      ! Total number of pixels nx*ny = np_vert*np_horiz
   !!----    integer                                    :: nframes     ! Total number of frames
   !!----    integer                                    :: nbang       ! Total number of angles moved during scan
   !!----    integer, dimension(11)                     :: icdesc      ! Integer values
   !!----    real(kind=cp),allocatable,dimension(:,:)   :: tmc_ang     ! time,monitor,total counts, angles*1000
   !!----                                                              ! To be allocated as tmc_ang(nbang,nframes)
   !!----    real(kind=cp),allocatable,dimension(:,:)   :: counts      ! Counts array to be reshaped (np_vert,np_horiz,nframes) in case of 2D detectors
   !!----                                                              ! To be allocated as counts(nbdata,nframes)
   !!---- End Type SXTAL_Numor_type
   !!----
   !!----    Definition for XTAL Numor type
   !!----
   !!---- Update: April - 2008
   !!
   Type, public :: SXTAL_Numor_type
      integer                                    :: numor       ! Numor
      integer                                    :: manip       ! principle scan angle
      integer                                    :: icalc       ! angle calculation type
      character(len=32)                          :: header      ! User, local contact, date
      character(len=4)                           :: Instrm      ! Instrument name
      character(len=32)                          :: title       !
      character(len=8)                           :: Scantype    ! omega, phi, etc...
      real(kind=cp), dimension(3)                :: hmin        ! or h,k,l for omega-scans
      real(kind=cp), dimension(3)                :: hmax        !
      real(kind=cp), dimension(5)                :: angles      ! Angles: phi, chi, omega, 2theta(gamma), psi
      real(kind=cp), dimension(3,3)              :: UB          ! UB-matrix
      real(kind=cp), dimension(3)                :: dh          ! delta_h, delta_k, delta_l
      real(kind=cp), dimension(3)                :: scans       ! scan start, scan step, scan width
      real(kind=cp)                              :: preset      !
      real(kind=cp)                              :: wave        ! wavelength
      real(kind=cp)                              :: cpl_fact    ! Coupling Factor
      real(kind=cp), dimension(5)                :: conditions  ! Temp-s.pt,Temp-Regul,Temp-sample,Voltmeter,Mag.field
      integer                                    :: nbdata      ! Total number of pixels nx*ny = np_vert*np_horiz
      integer                                    :: nframes     ! Total number of frames
      integer                                    :: nbang       ! Total number of angles moved during scan
      integer, dimension(11)                     :: icdesc      ! Integer values
      real(kind=cp),allocatable,dimension(:,:)   :: tmc_ang     ! time,monitor,total counts, angles*1000
                                                                ! To be allocated as tmc_ang(nbang,nframes)
      real(kind=cp),allocatable,dimension(:,:)   :: counts      ! Counts array to be reshaped (np_vert,np_horiz,nframes) in case of 2D detectors
                                                                ! To be allocated as counts(nbdata,nframes)
   End type SXTAL_Numor_type

   !!----
   !!---- TYPE :: SXTAL_ORIENT_TYPE
   !!--..
   !!---- Type, public :: SXTAL_Orient_type
   !!----    real(kind=cp)                :: wave       ! Wavelength (in Laue machines any value between Lambda_min and Lambda_max)
   !!----    real(kind=cp),dimension(3,3) :: UB         ! UB matrix in Busing-Levy setting
   !!----    real(kind=cp),dimension(3,3) :: UBINV      ! Inverse of UB-matrix
   !!----    real(kind=cp),dimension(3,3) :: CONV       ! Conversion matrix to the local setting
   !!---- End Type SXTAL_Orient_type
   !!----
   !!----    Definition for XTAL Orientation Parameters
   !!----
   !!---- Update: April - 2008
   !!
   Type, public :: SXTAL_Orient_type
      real(kind=cp)                :: wave       !Wavelength (in Laue machines any value between Lambda_min and Lambda_max)
      real(kind=cp),dimension(3,3) :: UB         !UB matrix in Busing-Levy setting
      real(kind=cp),dimension(3,3) :: UBINV      !Inverse of UB-matrix
      real(kind=cp),dimension(3,3) :: CONV       !Conversion matrix to the local setting
   End type SXTAL_Orient_type

   !!----
   !!---- CURRENT_INSTRM
   !!----    type(diffractometer_type), public :: Current_Instrm
   !!----
   !!----    Define a Current Instrument varuable
   !!----
   !!---- Update: April - 2008
   !!
   type(diffractometer_type), public :: Current_Instrm

   !!--++
   !!--++ CURRENT_INSTRM_SET
   !!--++    logical, private :: Current_Instrm_set
   !!--++
   !!--++    set to .true. if the Current_Instrm has been set up by calling
   !!--++    to Read_Current_Instrm
   !!--++
   !!--++ Update: April - 2008
   !!
   logical, private :: Current_Instrm_set=.false.

   !!----
   !!---- CURRENT_ORIENT
   !!----    type(SXTAL_Orient_type), public :: Current_Orient
   !!----
   !!----    Define a Current variable
   !!----    instrument
   !!----
   !!---- Update: April - 2008
   !!
   type(SXTAL_Orient_type), public :: Current_Orient

   !!----
   !!---- CYCLE_NUMBER
   !!----    integer, public :: cycle_number
   !!----
   !!----    Value to give the cycle number of Reactor at ILL
   !!----
   !!---- Update: April - 2008
   !!
   integer, public  ::  cycle_number

   !!----
   !!---- ERR_ILLDATA
   !!----    logical, public:: ERR_ILLData
   !!----
   !!----    logical variable to taking the value .true. if an error in the module
   !!----    ILL_INSTRM_DATA occurs.
   !!----
   !!---- Update: April - 2008
   !!
   logical, public :: ERR_ILLData=.false.

   !!----
   !!---- ERR_ILLDATA_MESS
   !!----    character(len=150), public:: ERR_ILLData_Mess
   !!----
   !!----    String containing information about the last error
   !!----
   !!---- Update: April - 2008
   !!
   character(len=150), public :: ERR_ILLData_Mess=" "

   !!----
   !!---- ILL_DATA_DIRECTORY
   !!----    character(len=512), public :: ILL_data_directory
   !!----
   !!----    String containing information about the data directory for ILL
   !!----    Initialised (automatic save attribute) for Windows case but set
   !!----    in subroutine: Initialize_Data_Directory
   !!----
   !!---- Update: March - 2009
   !!
   character(len=512), public  ::  ILL_Data_Directory = "c:\diffraction_Windows\illdata\"

   !!----
   !!---- ILL_TEMP_DIRECTORY
   !!----    character(len=512), private :: ILL_temp_directory
   !!----
   !!----
   !!----    String containing information about the data directory for ILL
   !!----    Initialised (automatic save attribute) for Windows case but set
   !!----    In subroutine: Initialize_Data_Directory
   !!----
   !!---- Update: 10/03/2011
   !!
   character(len=512), public ::  ILL_Temp_directory = "c:\diffraction_Windows\illdata\"

   !!--++
   !!--++ IVALUES
   !!--++    integer, dimension(:), allocatable, private :: ivalues
   !!--++
   !!--++    (Private)
   !!--++    Integer vector to read integer values on keytype blocks
   !!--++
   !!--++ Update: April - 2008
   !!
   integer, dimension(:), allocatable, private :: ivalues

   !!----
   !!---- INSTRM_DIRECTORY
   !!----    character(len=512), public :: Instrm_directory
   !!----
   !!----    String containing information about the data directory for specific
   !!----    instrument
   !!----
   !!---- Update: April - 2008
   !!
   character(len=512), public  :: Instrm_Directory = " "

   !!--++
   !!--++ INSTRM_DIRECTORY_SET
   !!--++    logical, private:: Instrm_directory_set
   !!--++
   !!--++    logical variable taking the value .true. if set the instrument directory
   !!--++    in correct form
   !!--++
   !!--++ Update: April - 2008
   !!
   logical, private :: Instrm_Directory_Set=.false.
   !!----
   !!---- INSTRM_GEOMETRY_DIRECTORY
   !!----    character(len=512), public :: Instrm_Geometry_directory
   !!----
   !!----    String containing information about the directory containing .geom files
   !!----    for specific instrument
   !!----
   !!---- Update: July - 2010
   !!
   character(len=512), public  :: Instrm_Geometry_Directory = " "

   !!--++
   !!--++ INSTRM_GEOMETRY_DIRECTORY_SET
   !!--++    logical, private:: Instrm_Geometry_Directory_Set
   !!--++
   !!--++    logical variable taking the value .true. if set the instrument geometry directory
   !!--++    in correct form
   !!--++
   !!--++ Update: April - 2010
   !!
   logical, private :: Instrm_Geometry_Directory_Set=.false.

   !!--++
   !!--++ GOT_ILL_DATA_DIRECTORY
   !!--++    logical, private:: got_ILL_data_directory
   !!--++
   !!--++    Logical variable taking the value .true. if the instrument directory
   !!--++    has been obtained from subroutine Initialize_Data_Directory
   !!--++
   !!--++ Update: March - 2009
   !!
   logical, private :: got_ILL_data_directory = .false.

   !!--++
   !!--++ GOT_ILL_TEMP_DIRECTORY
   !!--++    logical, private:: got_ILL_temp_directory
   !!--++
   !!--++    Logical variable taking the value .true. if the instrument directory
   !!--++    for temporary uncompressed data has been obtained from subroutine
   !!--++    Initialize_Data_Directory
   !!--++
   !!--++ Update: March - 2009
   !!
   logical, private :: got_ILL_temp_directory = .false.

   !!----
   !!---- MACHINE_NAME
   !!----    character(len=8), public :: machine_name
   !!----
   !!----    String containing information about the Instrument name
   !!----
   !!---- Update: April - 2008
   !!
   character(len=8),   public  ::  machine_name

   !!--++
   !!--++ N_KEYTYPES
   !!--++    Integer, dimension(7), private :: N_KEYTYPES
   !!--++
   !!--++    Integer containing the number of keytypes into the numor
   !!--++    according to:
   !!--++       Index   1    2    3    4    5    6    7
   !!--++     KeyType   R    A    S    F    I    J    V
   !!--++
   !!--++ Update: March - 2009
   !!
   integer, dimension(7), private :: n_keytypes

   !!--++
   !!--++ NL_KEYTYPES
   !!--++    Integer, dimension(:,:,:), allocatable, private :: NL_KEYTYPES
   !!--++
   !!--++    Integer array where write the initial and final lines for each
   !!--++    keytype group.
   !!--++    Index 1: [1-7]
   !!--++    Index 2: Number of Blocks of this keytype
   !!--++    Index 3: [1,2], The first the beginning and the second the end
   !!--++
   !!--++ Update: March - 2009
   !!
   integer, dimension(:,:,:), allocatable, private :: nl_keytypes

   !!--++
   !!--++ NTEXT
   !!--++    Integer, private :: NTEXT
   !!--++
   !!--++    Integer containing the number of lines of Text for KeyType
   !!--++
   !!--++ Update: March - 2009
   !!
   integer, private :: ntext

   !!--++
   !!--++ NVAL_F
   !!--++    Integer, private :: NVAL_F
   !!--++
   !!--++    Integer containing the number of real values readed on keytype blocks
   !!--++
   !!--++ Update: March - 2009
   !!
   integer, private :: nval_f

   !!--++
   !!--++ NVAL_I
   !!--++    Integer, private :: NVAL_I
   !!--++
   !!--++    Integer containing the number of integer values readed on keytype blocks
   !!--++
   !!--++ Update: March - 2009
   !!
   integer, private :: nval_i

   !!--++
   !!--++ RVALUES
   !!--++    real(kind=cp), dimension(:), allocatable, private :: rvalues
   !!--++
   !!--++    (Private)
   !!--++    Vector to read real values on keytype blocks
   !!--++
   !!--++ Update: April - 2008
   !!
   real(kind=cp), dimension(:), allocatable, private :: rvalues

   !!--++
   !!--++ TEXT_ILL
   !!--++    character(len=80), dimension(:), allocatable, private :: text_ill
   !!--++
   !!--++    String containing the info from Numors blocks
   !!--++
   !!--++ Update: March - 2009
   !!
   character(len=80), dimension(:), allocatable, private :: text_ill

   !!--++
   !!--++ UNCOMPRESSCOMMAND
   !!--++    character(len=512), private :: uncompresscommand
   !!--++
   !!--++    String containing the command for uncompressing compressed data files
   !!--++
   !!--++ Update: March - 2009
   !!
   character(len=512), private  ::  uncompresscommand=' '

   !!----
   !!---- YEAR_ILLDATA
   !!----    Integer, public :: YEAR_ILLDATA
   !!----
   !!----    Integer containing the two last figures of the "year" needed to
   !!----    find datafiles
   !!----
   !!---- Update: March - 2009
   !!
   integer,            public  ::  year_illdata

   !!--++
   !!--++ SET_CALIBRATION_DETECTOR
   !!--++    logical, private:: Set_Calibration_Detector
   !!--++
   !!--++    Logical variable taking the value .true. if the Calibration Detector
   !!--++    is readed using the respective routine.
   !!--++
   !!--++ Update: 17/03/2011
   !!
   logical, private :: Set_Calibration_Detector = .false.

   !---- Interfaces - Overloaded ----!
   Interface  Allocate_Numors
      Module Procedure Allocate_Powder_Numors
      Module Procedure Allocate_SXTAL_Numors
   End Interface

   Interface  Read_Numor
      Module Procedure Read_Powder_Numor
      Module Procedure Read_SXTAL_Numor
   End Interface

   Interface  Write_Numor_Info
      Module Procedure Write_Powder_Numor
      Module Procedure Write_SXTAL_Numor
   End Interface

   Interface  Write_HeaderInfo_Numor
      Module Procedure Write_HeaderInfo_Powder_Numor
      Module Procedure Write_HeaderInfo_SXTAL_Numor
   End Interface

 Contains

    !!----
    !!---- Subroutine Allocate_POWDER_numors(num_max,ndata,num_ang,nframes,Num)
    !!----    integer,                                            intent(in)    :: num_max
    !!----    integer,                                            intent(in)    :: ndata
    !!----    integer,                                            intent(in)    :: num_ang
    !!----    integer,                                            intent(in)    :: nframes
    !!----    type(POWDER_Numor_type), allocatable, dimension(:), intent(in out):: Num
    !!----
    !!--<<    Subroutine allocating and initializing the array 'Num' of type
    !!----    Powder_Numor_type. The input arguments are:
    !!----        num_max : number of components of the array (dimension of Num)
    !!----        ndata   : number of pixels of a single frame
    !!----        num_ang : number of angles moved simultaneously during the scan
    !!----        nframes : number of frames (number of scan points)
    !!-->>
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_POWDER_Numors(Num_Max,Ndata,Num_Ang,Nframes,Num)
       !---- Arguments ----!
       integer, intent(in)                                                :: num_max
       integer, intent(in)                                                :: ndata
       integer, intent(in)                                                :: num_ang
       integer, intent(in)                                                :: nframes
       type(Powder_Numor_type), allocatable, dimension(:), intent(in out) :: Num

       !---- Local variables ----!
       integer :: i

       if (allocated(Num)) deallocate(Num)
       allocate(Num(num_max))

       do i=1, num_max
          Num(i)%numor=0
          Num(i)%manip=0
          Num(i)%icalc=0
          Num(i)%header=" "
          Num(i)%Instrm=" "
          Num(i)%title=" "
          Num(i)%scantype=" "
          Num(i)%angles=0.0
          Num(i)%scans=0.0
          Num(i)%monitor=0.0
          Num(i)%time=0.0
          Num(i)%wave=0.0
          Num(i)%conditions=0.0
          Num(i)%nbdata=0
          Num(i)%nbang=0
          Num(i)%nframes=0

          if(allocated(Num(i)%tmc_ang)) deallocate(Num(i)%tmc_ang)
          if (num_ang > 0 .and. nframes > 0) then
             allocate(Num(i)%tmc_ang(num_ang,nframes))
             Num(i)%tmc_ang=0.0
          end if

          if(allocated(Num(i)%counts)) deallocate(Num(i)%counts)
          if (ndata > 0 .and. nframes > 0) then
             allocate(Num(i)%counts(ndata,nframes))
             Num(i)%counts=0.0
          end if
       end do

       return
    End Subroutine Allocate_POWDER_Numors

    !!----
    !!---- Subroutine Allocate_SXTAL_Numors(Num_Max,Ndata,Num_Ang,Nframes,Num)
    !!----    integer,                                           intent(in)     :: num_max
    !!----    integer,                                           intent(in)     :: ndata
    !!----    integer,                                           intent(in)     :: num_ang
    !!----    integer,                                           intent(in)     :: nframes
    !!----    type(SXTAL_Numor_type), allocatable, dimension(:), intent(in out) :: Num
    !!----
    !!--<<    Subroutine allocating and initializing the array 'Num' of type SXTAL_Numor_type
    !!----    The input arguments are:
    !!----    num_max : number of components of the array (dimension of Num)
    !!----    ndata   : number of pixels of a single frame
    !!----    num_ang : number of angles moved simultaneously during the scan
    !!-->>    nframes : number of frames (number of scan points)
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_SXTAL_numors(num_max,ndata,num_ang,nframes,Num)
       !---- Arguments ----!
       integer,                                           intent(in)     :: num_max
       integer,                                           intent(in)     :: ndata
       integer,                                           intent(in)     :: num_ang
       integer,                                           intent(in)     :: nframes
       type(SXTAL_Numor_type), allocatable, dimension(:), intent(in out) :: Num

       !---- Local variables ----!
       integer :: i

       if(allocated(Num)) deallocate(Num)
       allocate(Num(num_max))

       do i=1, num_max
          Num(i)%numor=0
          Num(i)%manip=0
          Num(i)%icalc=0
          Num(i)%header=" "
          Num(i)%instrm=" "
          Num(i)%title=" "
          Num(i)%scantype=" "
          Num(i)%hmin=0.0
          Num(i)%hmax=0.0
          Num(i)%dh=0.0
          Num(i)%angles=0.0
          Num(i)%ub=0.0
          Num(i)%scans=0.0
          Num(i)%preset=0.0
          Num(i)%wave=0.0
          Num(i)%cpl_fact=0.0
          Num(i)%conditions=0.0
          Num(i)%nbdata=0
          Num(i)%nbang=0
          Num(i)%nframes=0

          if(allocated(Num(i)%tmc_ang)) deallocate(Num(i)%tmc_ang)
          if (num_ang > 0 .and. nframes > 0) then
             allocate(Num(i)%tmc_ang(num_ang,nframes))
             Num(i)%tmc_ang=0.0
          end if

          if(allocated(Num(i)%counts)) deallocate(Num(i)%counts)
          if (ndata > 0 .and. nframes > 0) then
             allocate(Num(i)%counts(ndata,nframes))
             Num(i)%counts=0.0
          end if
       end do

       return
    End Subroutine Allocate_SXTAL_Numors

    !!----
    !!---- Subroutine Define_Uncompress_Program(ProgName)
    !!----    character(len=*), intent(in) :: ProgName
    !!----
    !!---- Routine that define the uncompress program that you wants to use
    !!----
    !!---- Update:  April - 2009
    !!
    Subroutine Define_Uncompress_Program(ProgName)
       !---- Argument ----!
       character(len=*), intent(in) :: ProgName

       uncompresscommand=ProgName

       return
    End Subroutine Define_Uncompress_Program

    !!----
    !!---- Subroutine Get_Absolute_Data_Path(Numor,Instrm,Path,Iyear,Icycle, Actual_Path)
    !!----    integer,           intent(in)            :: numor
    !!----    character(len=*),  intent(in)            :: instrm
    !!----    character(len=*),  intent(out)           :: path
    !!----    integer, optional, intent(in)            :: iyear
    !!----    integer, optional, intent(in)            :: icycle
    !!----    character(len=*),  intent(out), optional :: Actual_Path
    !!----
    !!----    Finds the absolute path to any numor. The base directory
    !!----    is set by a call to 'initialize_data_directory'. The subroutine
    !!----    then searches for the numor in the following order:
    !!----    1. In the subdirectory defined by the year and cycle if passed
    !!----       as arguments to the subroutine (i.e args iyear, icycle).
    !!----    2. In the subdirectory defined by the year and cycle of the
    !!----       previous call to get_absolute_data_path (since numors are
    !!----       likely to be adjacent).
    !!----    3. In the 'data' subdirectory (since likely to process recent
    !!----       data).
    !!----    4. In the 'data-1' subdirectory (same logic as above)
    !!----    5. Working from the current year and most recent cycle and
    !!----       working back through cycles and year until the stopping
    !!----       at the first cycle of 1973, when the first data was archived
    !!----
    !!--..    Tries to find an uncompress numor first and then tries to find
    !!--..    a compressed numor (.Z extension). If found the numor is uncomp-
    !!--..    ressed in the a temporary directory if defined (see subroutine
    !!--..    'initialize_data_directory') or else into the current directory.
    !!--..
    !!--..    Nb At present no attempt is made to tidy up these uncompressed
    !!--..    numors, which could potentially litter the current directory.
    !!----
    !!---- Update: March - 2009
    !!
    Subroutine Get_Absolute_Data_Path(Numor, Instrm, Path, Iyear, Icycle, Actual_Path)

       !---- Arguments ----!
       integer,           intent(in)            :: numor
       character(len=*),  intent(in)            :: instrm
       character(len=*),  intent(out)           :: path
       integer, optional, intent(in)            :: iyear
       integer, optional, intent(in)            :: icycle
       character(len=*),  intent(out), optional :: Actual_Path

       !---- Local Variables ----!
       character(len=*),parameter :: Current_Data = 'data'
       character(len=*),parameter :: Previous_Data = 'data-1'
       character(len=*),parameter :: Extension = '.Z'

       character(len=6) :: numstr
       character(len=4) :: inst
       character(len=5) :: yearcycle
       character(len=7) :: subdir
       logical          :: exists

       ! Init value
       path=" "
       if (present(actual_path)) actual_path = path

       ! Some compilers are unable to test the existence of file: data_diretory//ops_sep//"." !!!!!
       ! so, got_ILL_data_directory may still be .false., ILL_data_directory has been set anyway
       ! in initialize_data_directory(), so let us continue and set got_ILL_data_directory=.true.
       ! if the associated numor file exist anywhere
       if (.not. got_ILL_data_directory) call initialize_data_directory()

       ! At this point the ILL Data and Temporal directories must be defined
       if (.not. got_ILL_Data_directory) return
       if (.not. got_ILL_Temp_Directory) return

       ! Uncompress program must be defined
       if (len_trim(uncompresscommand) == 0) then
          if (Ops == 1) then
             call define_uncompress_program('7z e')
          else
             call define_uncompress_program('gzip -q -d -c')
          end if
       end if

       ! Numor character
       write(unit=numstr,fmt='(i6.6)') numor

       ! Instrument
       inst=adjustl(instrm)
       call lcase(inst)

       ! Using Instrument Information
       select case(inst)
          case("d1b","d20","d9","d15","d19","d16")
             subdir = trim(inst)//"_"//numstr(2:2)//ops_sep

          case default ! d10, d3 etc don't divide numors up into sub directories
             subdir = ""
       end select

       if (present(iyear) .and. present(icycle)) then
          !> Using Year+Cycle
          write(unit=yearcycle,fmt="(i4.4,i1.1)") iyear,icycle
          yearcycle = yearcycle(len_trim(yearcycle)-2:len_trim(yearcycle))
          path = trim(ILL_data_directory)//trim(yearcycle)//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
          if (present(actual_path)) actual_path = path
          inquire(file=trim(path),exist=exists)
          if (exists) return ! found numor so return

          !> Using previous + compressed data
          path = trim(path)//Extension
          inquire(file=trim(path),exist=exists)
          if (exists) then ! uncompress into temp directory
             call system(trim(uncompresscommand)//' '//trim(path)//' > '//trim(ILL_temp_directory)//numstr)
             if (present(actual_path)) actual_path = path
             path = trim(ILL_temp_directory)//numstr
             return ! found numor so return
          end if

          ! The next, not haven't sense if you use the Year and cycle     ! by JGP

          !> Using Year + Cycle of a previous call to routine
          !write(unit=yearcycle,fmt="(i4.4,i1.1)") year_illdata, cycle_number
          !yearcycle = yearcycle(len_trim(yearcycle)-2:len_trim(yearcycle))
          !
          !path = trim(ILL_data_directory)//trim(yearcycle)//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
          !if (present(actual_path)) actual_path = path
          !inquire(file=trim(path),exist=exists)
          !if (exists) return ! found numor so return
          !
          !! Using previous + compressed data
          !path = trim(path)//Extension
          !if (present(actual_path)) actual_path = path
          !inquire(file=trim(path),exist=exists)
          !if (exists) then ! uncompress into temp directory
          !   call system(trim(uncompresscommand)//' '//trim(path)//' > '//trim(ILL_temp_directory)//numstr)
          !   if (present(actual_path)) actual_path = path
          !   path = trim(ILL_temp_directory)//numstr
          !   return ! found numor so return
          !end if
       end if

       !> Using Current_data
       path = trim(ILL_data_directory)//Current_Data//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
       if (present(actual_path)) actual_path = path
       inquire(file=trim(path),exist=exists)
       if (exists) return ! found numor so return

       ! Using Previous data
       path = trim(ILL_data_directory)//Previous_Data//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
       if (present(actual_path)) actual_path = path
       inquire(file=trim(path),exist=exists)
       if (exists) return ! found numor so return

       ! start from the most recent yearcycle and work search backwards
       call get_next_yearcycle(yearcycle,.true.)
       do
          if (yearcycle == "") exit

          path = trim(ILL_data_directory)//trim(yearcycle)//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
          if (present(actual_path)) actual_path = path
          inquire(file=trim(path),exist=exists)
          if (exists) return

          ! Using previous + compressed data
          path = trim(path)//Extension
          inquire(file=trim(path),exist=exists)
          if (exists) then ! uncompress into temp directory
             call system(trim(uncompresscommand)//' '//trim(path)//' > '//trim(ILL_temp_directory)//numstr)
             if (present(actual_path)) actual_path = path
             path = trim(ILL_temp_directory)//numstr
             return ! found numor so return
          end if

          call get_next_yearcycle(yearcycle)
       end do

       ! the numor wasn't found compressed or uncompressed anywhere
       path = " "
       if (present(actual_path)) actual_path = path

       return
    End Subroutine Get_Absolute_Data_Path

    !!----
    !!---- Subroutine Get_Next_YearCycle(YearCycle,Reset_To_Most_Recent)
    !!----    character(len=*), intent(out) :: yearcycle
    !!----    logical, optional, intent(in) :: reset_to_most_recent
    !!----
    !!----    Works back through the cycles and years, returning the previous
    !!----    yearcycle combination as a 3 character string i.e.
    !!----       if year_illdata = 6 and cycle_number = 5, returns '064'
    !!----       if year_illdata = 6 and cycle_number = 1, returns '057'
    !!----
    !!----    The reset_to_most_recent flag allows the year_illdata and cycle_number to
    !!----    be set to the most recent possible. If asked for a yearcycle before
    !!----    '731' (first cycle of 1973) then returns "", since no data was archived
    !!----    before this date.
    !!----    (Original code from Mike Turner)
    !!----
    !!---- Update: March - 2009
    !!
    Subroutine Get_Next_YearCycle(Yearcycle, reset_to_most_recent)
       !---- Argument ----!
       character(len=*), intent(out) :: yearcycle
       logical, optional, intent(in) :: reset_to_most_recent

       !---- Local Variables ----!
       integer, parameter    :: Maxcycle = 7
       integer, dimension(8) :: dt
       character(len=10)     :: date,time,zone

       if (present(reset_to_most_recent)) then
          if (reset_to_most_recent) then
             call date_and_time (date, time, zone, dt)
             year_illdata = mod(dt(1),100) ! last two digits of year
             cycle_number = Maxcycle
          end if
       else
          if (cycle_number /= 1) then
             cycle_number = cycle_number -1
          else
             if (year_illdata == 0) then
                year_illdata = 99
             else
                year_illdata = year_illdata - 1
             end if
             cycle_number = Maxcycle
          end if
       end if

       if (year_illdata == 72) then ! no more point searching
          yearcycle = ""
       else
          write(unit=yearcycle,fmt="(i4.4,i1.1)") year_illdata,cycle_number
          yearcycle = yearcycle(len_trim(yearcycle)-2:len_trim(yearcycle))

       end if

       return

    End Subroutine Get_Next_YearCycle

    !!----
    !!---- Subroutine Get_Single_Frame_2D(nfr,iord,snum,dat_2D,appl_alphas)
    !!----    integer,                       intent(in)  :: nfr,iord
    !!----    type(SXTAL_Numor_type),        intent(in)  :: snum
    !!----    real(kind=cp), dimension(:,:), intent(out) :: dat_2D
    !!----    logical, optional,             intent(in)  :: appl_alphas
    !!----
    !!--<<    Extracts into the real two-dimensional array dat_2D the counts
    !!----    of the frame number 'nfr' of the numor object 'snum', applying
    !!----    the efficiency corrections depending of the optional argument
    !!----    'appl_alphas'. The value of iord selects the type of ordering
    !!----    of counts in 'snum'
    !!----
    !!----    IORD = 1 for D19 Banana, IORD = 2 for D9/D10, IORD = 3 for D19 Flat
    !!-->>    IORD = 4 for ID20 trial.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Get_Single_Frame_2D(nfr,iord,snum,dat_2D,appl_alphas)
       !---- Arguments ----!
       integer,                       intent(in)  :: nfr,iord
       type(SXTAL_Numor_type),        intent(in)  :: snum
       real(kind=cp), dimension(:,:), intent(out) :: dat_2D
       logical, optional,             intent(in)  :: appl_alphas

       !---- IA = anodes (1-NANO), IC = cathodes (1-NCAT)
       !   IORD = 1 for D19 Banana, IORD = 2 for D9/D10, IORD = 3 for D19 Flat.
       !   IORD = 4 for ID20 trial.
       integer :: k,ia,ic,ncat,nano

       Select Case(iord)

         Case(1)   ! D19
           ! Read vertical-wire by vertical-wire, with zero at high nu, high gamma.
           ! D19 Banana
           ncat =  Current_instrm%np_vert
           nano =  Current_instrm%np_horiz
           k=0
           do ic=1,ncat
             do ia=1,nano
               k=k+1
               dat_2D(ia,ic)=snum%counts(k,nfr)
             end do
           end do

         Case(2)   ! D9/D10
           ! Read horizontal-wire by horizontal-wire, with zero at high nu, high gamma?
           ! D9, D10
           ncat =  Current_instrm%np_vert
           nano =  Current_instrm%np_horiz
           k=0
           do ic=1,ncat
             do ia=1,nano
               k=k+1
               dat_2D(ia,ic)=snum%counts(k,nfr)
             end do
           end do

         Case(3)   ! D19-Flat
           ! Read horizontal-wire by horizontal-wire, with zero at high nu, high gamma.
           ! D19 flat detector (but we reverse IA because of the funny different
           ! conversion for flat and curved detectors in subroutine d19amd.
           ncat =  Current_instrm%np_vert
           nano =  Current_instrm%np_horiz
           k=0
           do ia=nano,1,-1
             do ic=1,ncat
               k=k+1
               dat_2D(ia,ic)=snum%counts(k,nfr)
             end do
           end do

         Case(4)   ! ID20-trial
           ! Read vertical-wire by vertical-wire, with zero at low nu, high gamma.
           ! ID20 flat detector (but we reverse IA because of the funny different
           ! conversion for flat and curved detectors in subroutine d19amd).
           ncat =  Current_instrm%np_vert
           nano =  Current_instrm%np_horiz
           k=0
           do ic=ncat,1,-1
             do ia=nano,1,-1
               k=k+1
               dat_2D(ia,ic)=snum%counts(k,nfr)
             end do
           end do

       End Select

       if(present(appl_alphas)) then
         if(appl_alphas) dat_2D(:,:)= dat_2D(:,:)*Current_instrm%alphas(:,:)
       end if

       return
    End Subroutine Get_Single_Frame_2D

    !!----
    !!---- Subroutine Initialize_Data_Directory()
    !!----
    !!----    Call two subroutines: Initialize_Numors_Directory and Initialize_Temp_Directory .
    !!----    The first one is to set the ILL data directory and
    !!----    the second one to set the temporary directory.
    !!----    Both subroutines were originally coded by Mike Turner.
    !!----
    !!----
    !!---- Update: January - 2010
    !!
    Subroutine Initialize_Data_Directory()

       Call Initialize_Numors_Directory()
       Call Initialize_Temp_Directory()

       Return
    End Subroutine Initialize_Data_Directory

    !!----
    !!---- Subroutine Initialize_Numors_Directory()
    !!----
    !!....    Original code from Mike Turner (as well as the following comments)
    !!....    Depending on the operating system as reported by winteracter routine
    !!....    InfoOpSystem, assigns values to the following global public variables:
    !!....    sep, the path separator. "sep" has been replaced by "ops_sep" from
    !!....    CFML_GlobalDeps. InfoOpSystem is a function from Winteracter: replaced by
    !!....    OPS integer from CFML_GlobalDeps
    !!....
    !!----    Subroutine assigning values to the following global public variables:
    !!----    ILL_Data_Directory: Base directory for the ILL data. If found then
    !!----                        the flag got_ILL_data_directory is set to true.
    !!----
    !!--..    Original code from Mike Turner, changed to make the subroutine independent
    !!--..    of the Winteracter library.
    !!--..
    !!--..    The subroutine ask for enviroment variables TEMP and ILLDATA in first place and
    !!--..    and if it return a blank then the default is set as following:
    !!--..    ILL_Data_Directory -> \\Serdon\illdata
    !!--..    ILL_Data_Directory -> /net/serdon/illdata
    !!----
    !!----
    !!---- Update: January - 2010
    !!
    Subroutine Initialize_Numors_Directory()
       !---- Local Variables ----!
       Logical                     :: Exists
       Character(Len=*), Parameter :: Envvar = 'ILLDATA'
       Integer                     :: Nlong

       Exists = .False.

       !> Exist ILLDATA environment variable?
       Call Get_Environment_Variable(Envvar, Ill_Data_Directory)

       If (Len_Trim(Ill_Data_Directory) == 0) Then
          If (Ops == 1) Then
             Ill_Data_Directory = '\\Serdon\illdata'                     ! For Windows
          Else
             Ill_Data_Directory = '/net/serdon/Illdata'                  ! For Linux/MACoS
          End If
       End If

       nlong = Len_Trim(Ill_Data_Directory)
       if (nlong > 0 ) then
          If (Ill_Data_Directory(nlong:nlong) /= Ops_Sep) Ill_Data_Directory = trim(Ill_Data_Directory)//Ops_Sep
       end if
       Got_Ill_Data_Directory = Directory_Exists(Trim(Ill_Data_Directory))

       Return

    End Subroutine Initialize_Numors_Directory

    !!----
    !!---- Subroutine Initialize_Temp_Directory()
    !!----
    !!....    Original code from Mike Turner (as well as the following comments)
    !!....    Depending on the operating system as reported by winteracter routine
    !!....    InfoOpSystem, assigns values to the following global public variables:
    !!....    sep, the path separator. "sep" has been replaced by "ops_sep" from
    !!....    CFML_GlobalDeps. InfoOpSystem is a function from Winteracter: replaced by
    !!....    OPS integer from CFML_GlobalDeps
    !!....
    !!----    Subroutine assigning values to the following global public variables:
    !!----    ILL_Temp_Directory: A temporary directory (used for uncompressing). If
    !!----                        found the flag got_ILL_temp_directory is set to true.
    !!----
    !!--..    Original code from Mike Turner, changed to make the subroutine independent
    !!--..    of the Winteracter library.
    !!--..
    !!--..    The subroutine ask for enviroment variables TEMP and ILLDATA in first place and
    !!--..    and if it return a blank then the default is set as following:
    !!--..    Windows: ILL_Temp_Directory -> C:\Temp
    !!--..    Linux:   ILL_Temp_Directory -> $HOME/tmp
    !!----
    !!----
    !!---- Update: January - 2010
    !!
    Subroutine Initialize_Temp_Directory()

       !---- Parameters ----!
       Character (Len=*), Parameter :: Envvar1 = 'TEMP'
       Character (Len=*), Parameter :: Envvar2 = 'TMP'

       !---- Local Variables ----!
       Integer                      :: I

       !> Is defined TEMP environment variable?
       Call Get_Environment_Variable(Envvar1, Ill_Temp_Directory)

       If (Len_Trim(Ill_Temp_Directory) == 0) Then
          !> Is defined TMP environment variable?
          Call Get_Environment_Variable(Envvar2, Ill_Temp_Directory)
          If (Len_Trim(Ill_Temp_Directory) == 0) Then
             !> Default values set to Home user
             If (Ops == 1) Then
                Call Get_Environment_Variable('USERPROFILE', Ill_Temp_Directory)
             Else
                Call Get_Environment_Variable('HOME', Ill_Temp_Directory)
             End If
          End If
       End If

       I = Len_Trim(Ill_Temp_Directory)
       If (Ill_Temp_Directory(I:I) /= Ops_Sep) Ill_Temp_Directory = Trim(Ill_Temp_Directory)//Ops_Sep

       ! Temporal
       Got_Ill_Temp_Directory = Directory_Exists(Trim(Ill_Temp_Directory))

       Return
    End Subroutine Initialize_Temp_Directory


    !!--++
    !!--++ Subroutine Number_KeyTypes_on_File(filevar, nlines)
    !!--++    character(len=*),dimension(:), intent(in) :: filevar
    !!--++    integer,                       intent(in) :: nlines
    !!--++
    !!--++    1:R, 2:A, 3:S, 4:F, 5:I, 6:J, 7:V
    !!
    Subroutine Number_KeyTypes_on_File(filevar, nlines)
       !---- Arguments ----!
       character(len=*),dimension(:), intent(in) :: filevar
       integer,                       intent(in) :: nlines

       !---- Local Variables ----!
       integer :: i

       ! Init output variables
       n_keytypes=0

       do i=1,nlines
          select case (filevar(i)(1:10))
             case ('RRRRRRRRRR')
                n_keytypes(1)=n_keytypes(1)+1
             case ('AAAAAAAAAA')
                n_keytypes(2)=n_keytypes(2)+1
             case ('SSSSSSSSSS')
                n_keytypes(3)=n_keytypes(3)+1
             case ('FFFFFFFFFF')
                n_keytypes(4)=n_keytypes(4)+1
             case ('IIIIIIIIII')
                n_keytypes(5)=n_keytypes(5)+1
             case ('JJJJJJJJJJ')
                n_keytypes(6)=n_keytypes(6)+1
             case ('VVVVVVVVVV')
                n_keytypes(7)=n_keytypes(7)+1
          end select
       end do

       return
    End Subroutine Number_KeyTypes_on_File

    !!----
    !!---- Subroutine PowderNumor_To_DiffPattern(Numor,Pat)
    !!----    type(powder_numor_type),        intent(in)  :: Numor
    !!----    type(diffraction_pattern_type), intent(out) :: Pat
    !!----
    !!----
    !!---- Pass the information from Powder_Numor_Type to Diffraction_Pattern_type
    !!----
    !!----
    !!---- 10/03/2011
    !!
    Subroutine PowderNumor_To_DiffPattern(Numor,Pat)
       !---- Arguments ----!
       type(powder_numor_type),        intent(in)  :: Numor
       type(diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       integer :: i

       !> Initialize
       ERR_ILLData=.false.
       ERR_ILLData_Mess=' '

       if (numor%nbdata == 0) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess=' No points are defined in the present Numor'
          return
       end if

       ! Powder Instruments / Allocating points
       select case (trim(u_case(numor%instrm)))
          case ('D1A')
             call Allocate_Diffraction_Pattern(Pat,3000)

          case ('D1B','D20')
             call Allocate_Diffraction_Pattern(Pat,Numor%nbdata)

             ! Load values
             Pat%title=trim(Numor%title)

             pat%Scat_Var=trim(numor%ScanType)
             Pat%xmin=numor%scans(1)
             Pat%step=numor%scans(2)
             Pat%xmax=Pat%xmin+(Pat%npts-1)*Pat%step

             Pat%Monitor=Numor%monitor
             Pat%TSamp=Numor%conditions(3)
             Pat%TSet=Numor%conditions(1)

             Pat%Conv(1)=Numor%wave
             Pat%y=Numor%counts(:,1)
             do i=1,pat%npts
                if (pat%y(i) <= 0.00001) pat%y(i) = 1.0
                pat%sigma(i) = pat%y(i)
                pat%x(i)= pat%xmin+(i-1)*pat%step
             end do
             pat%ymin=minval(pat%y(1:pat%npts))
             pat%ymax=maxval(pat%y(1:pat%npts))
       end select

       return
    End Subroutine PowderNumor_To_DiffPattern

    !!--++
    !!--++ Subroutine Read_A_KeyType(filevar, n_ini, n_end, nchars, charline)
    !!--++    character(len=*), dimension(:), intent(in)  :: filevar   ! Input information
    !!--++    integer,                        intent(in)  :: n_ini     ! First line
    !!--++    integer,                        intent(in)  :: n_end     ! Last line to be read
    !!--++    integer,                        intent(out) :: nchars    ! Number of characters to be read from the next data
    !!--++    character(len=*),               intent(out) :: charline  ! String to be read
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Read_A_KeyType(filevar, n_ini, n_end, nchars, charline)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)  :: filevar
       integer,                        intent(in)  :: n_ini
       integer,                        intent(in)  :: n_end
       integer,                        intent(out) :: nchars
       character(len=*),               intent(out) :: charline

       !---- Local Variables ----!
       character(len=80)      :: line
       integer                :: i,j,n
       integer, dimension(10) :: ivet

       ! Init output values
       err_illdata=.false.

       nchars =0
       charline=' '
       ntext=0
       if (allocated(text_ill)) deallocate(text_ill)

       ! Check the correct KeyType A
       line=filevar(n_ini)
       if (line(1:10) /= 'AAAAAAAAAA') then
          err_illdata=.true.
          err_illdata_mess=' A bad Block Type has been found'
          return
       end if

       ! Getting information
       read(unit=filevar(n_ini+1),fmt='(10i8)') ivet
       nchars =ivet(1)
       ntext=ivet(2)

       if (ntext > 0) then
          allocate(text_ill(ntext))
          text_ill=' '
          do i=1,ntext
             j=(n_ini+2)+i
             if (j > n_end) then
                err_illdata=.true.
                err_illdata_mess=' Impossible to read a line for this block!'
                exit
             end if
             text_ill(i)=trim(filevar(j))
          end do
       end if

       if (nchars > 0) then
          charline=' '
          n=nchars/80
          if (mod(nchars,80) /= 0) n=n+1
          do i=1,n
             charline=trim(charline)//filevar(n_ini+ntext+1+i)//char(10)
          end do
       end if

       return
    End Subroutine Read_A_KeyType

    !!----
    !!---- Subroutine Read_Current_Instrm(filenam)
    !!----    character(len=*),  intent(in) :: filenam
    !!----
    !!----    Subroutine reading the file 'filenam' where the characteristics
    !!----    of the current instrument are written. The global Current_Instrm
    !!----    variable is filled after returning from this subroutine.
    !!----    In case of error the subroutine puts ERR_ILLData=.true.
    !!----    and fills the error message variable ERR_ILLData_Mess.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Read_Current_Instrm(filenam)
       !---- Argument ----!
       character(len=*),  intent(in) :: filenam

       !---- Local variables ----!
       character(len=120)             :: line
       character(len=10)              :: key
       integer                        :: i, j, lun, ier,npx,npz, n1,n2
       real(kind=cp), dimension(3,3)  :: set, ub
       real(kind=cp)                  :: wave

       logical                        :: read_set, read_alphas, read_ub

       npx=32  !default values
       npz=32
       call Get_LogUnit(lun)

       open(unit=lun,file=trim(filenam),status="old", action="read", position="rewind",iostat=ier)
       if (ier /= 0) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess="Error opening the file: "//trim(filenam)
          return
       end if

       read_set    = .false.
       read_alphas = .false.
       read_ub     = .false.

       Current_Instrm%det_offsets=0.0
       Current_Instrm%ang_Limits=0.0
       Current_Instrm%disp_Limits=0.0
       Current_Instrm%wave_min=0.5
       Current_Instrm%wave_max=0.5
       Current_Instrm%e1=(/1.0,0.0,0.0/)
       Current_Instrm%e2=(/0.0,1.0,0.0/)
       Current_Instrm%e3=(/0.0,0.0,1.0/)

       n1=0
       n2=0

       do
         read(unit=lun,fmt="(a)",iostat=ier) line
         if(ier /= 0) exit

         line=adjustl(line)
         if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
         i=index(line," ")
         key=u_case(line(1:i-1))

         Select Case(key)
            Case("INFO")
                Current_Instrm%info= adjustl(line(i+1:))

            Case("NAME")
                Current_Instrm%name_inst= adjustl(line(i+1:))

            Case("GEOM")
                Current_Instrm%geom= adjustl(line(i+1:))
                read(unit=Current_Instrm%geom,fmt=*) Current_Instrm%igeom
                j=index(Current_Instrm%geom," ")
                Current_Instrm%geom=adjustl(Current_Instrm%geom(j+1:))

            Case("BLFR")
                Current_Instrm%BL_frame = adjustl(line(i+1:))

            Case("DIST_UNITS")
                Current_Instrm%dist_units = adjustl(line(i+1:))

            Case("ANGL_UNITS")
                Current_Instrm%angl_units = adjustl(line(i+1:))

            Case("DET_TYPE")
                Current_Instrm%detector_type = adjustl(line(i+1:))
                Current_Instrm%ipsd=2 !Flat detector
                 j=index(Current_Instrm%detector_type,"ipsd")
                if( j /= 0)  &
                read(unit=Current_Instrm%detector_type(j+4:),fmt=*,iostat=ier) Current_Instrm%ipsd
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the type of detector (ipsd)"
                  return
                end if

            Case("DIST_DET")
                read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%dist_samp_detector
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the distance sample-detector"
                  return
                end if

            Case("WAVE")
                read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%wave_min
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the wavelength"
                  return
                else
                	Current_Instrm%wave_max=Current_Instrm%wave_min
                end if

            Case("WAVE_LIMITS")
                read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%wave_min, Current_Instrm%wave_max
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the wavelength limits"
                  return
                end if

            Case("DIM_XY")
                read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%horiz,    Current_Instrm%vert, &
                                                       Current_Instrm%np_horiz, Current_Instrm%np_vert
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the dimensions of the detector"
                  return
                end if
                npx = Current_Instrm%np_horiz
                npz = Current_Instrm%np_vert

            Case("GAPS_DET")
                read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%agap, Current_Instrm%cgap
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess=&
                  "Error in file: "//trim(filenam)//", reading the gaps between anodes and between cathodes"
                  return
                end if

            Case("UBMAT")
                do j=1,3
                  read(unit=lun,fmt=*,iostat=ier) ub(j,:)
                  if(ier /= 0) then
                    ERR_ILLData=.true.
                    ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the UB-matrix"
                    return
                  end if
                end do
                read_ub=.true.

            Case("SETTING")
                read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%e1, Current_Instrm%e2, Current_Instrm%e3
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the setting vectors"
                  return
                end if
                set(:,1)=Current_Instrm%e1
                set(:,2)=Current_Instrm%e2
                set(:,3)=Current_Instrm%e3
                read_set=.true.

            Case("NUM_ANG")
                 read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%num_ang
                 n1 = Current_Instrm%num_ang

            Case("NUM_DISP")
                 read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%num_disp
                 n2 = Current_Instrm%num_disp

            Case("ANG_LIMITS")
                if(n1 == 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", Number of angular motors missing!"
                  return
                end if
                do j=1,n1
                  read(unit=lun,fmt=*,iostat=ier) Current_Instrm%ang_names(j), Current_Instrm%ang_Limits(j,1:2), &
                                                  Current_Instrm%ang_offsets(j)
                  if(ier /= 0) then
                    ERR_ILLData=.true.
                    ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the angular limits"
                    return
                  end if
                end do

            Case("DISP_LIMITS")
                if(n2 == 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", Number of displacement motors missing!"
                  return
                end if
                do j=1,n2
                  read(unit=lun,fmt=*,iostat=ier) Current_Instrm%disp_names(j),Current_Instrm%disp_Limits(j,1:2), &
                                                  Current_Instrm%disp_offsets(1:n2)
                  if(ier /= 0) then
                    ERR_ILLData=.true.
                    ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the displacement limits"
                    return
                  end if
                end do

            Case("DET_OFF")
                read(unit=line(i+1:),fmt=*,iostat=ier) Current_Instrm%det_offsets
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the detector offsets"
                  return
                end if

            Case("DET_ALPHAS")
                if(allocated(Current_Instrm%alphas)) deallocate(Current_Instrm%alphas)
                allocate(Current_Instrm%alphas(npx,npz))

                do j=1,npz
                  read(unit=lun,fmt=*,iostat=ier) Current_Instrm%alphas(j,1:npx)
                  if(ier /= 0) then
                    ERR_ILLData=.true.
                    ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the efficiency of 2D-detector"
                   return
                  end if
                end do
                read_alphas= .true.

         End Select
       End do

       close(unit=lun)

       if(read_ub) then
         if (read_set) then
            call Set_Current_Orient(wave,ub,set)
         else
            call Set_Current_Orient(wave,ub)
         end if
       end if

       if (.not. read_alphas .and. index(Current_Instrm%geom,"Laue") == 0) then
          if (allocated(Current_Instrm%alphas)) deallocate(Current_Instrm%alphas)
          allocate(Current_Instrm%alphas(npx,npz))
          Current_Instrm%alphas(:,:)=1.0
       end if
       Current_Instrm_set=.true.

       return
    End Subroutine Read_Current_Instrm

    !!--++
    !!--++ Subroutine Read_F_KeyType(filevar, n_ini, n_end)
    !!--++    character(len=*), dimension(:), intent(in)  :: filevar   ! Input information
    !!--++    integer,                        intent(in)  :: n_ini     ! First line
    !!--++    integer,                        intent(in)  :: n_end     ! Last line to be read
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Read_F_KeyType(filevar, n_ini, n_end)
       !---- Arguments----!
       character(len=*), dimension(:), intent(in)  :: filevar
       integer,                        intent(in)  :: n_ini
       integer,                        intent(in)  :: n_end

       !---- Local Variables ----!
       character(len=80)      :: line
       integer                :: i,j,k,nl
       integer, dimension(10) :: ivet
       real, dimension(5)     :: vet

       ! Init output values
       err_illdata=.false.
       nval_f=0
       ntext=0

       if (allocated(text_ill)) deallocate(text_ill)
       if (allocated(rvalues)) deallocate(rvalues)

       ! Check the correct KeyType F
       line=filevar(n_ini)
       if (line(1:10) /= 'FFFFFFFFFF') then
          err_illdata=.true.
          err_illdata_mess=' A bad Block Type has been found'
          return
       end if

       ! Getting information
       read(unit=filevar(n_ini+1),fmt='(10i8)') ivet
       nval_f =ivet(1)
       ntext=ivet(2)

       if (ntext > 0) then
          allocate(text_ill(ntext))
          text_ill=' '
          do i=1,ntext
             j=(n_ini+1)+i
             text_ill(i)=filevar(j)
          end do
       end if

       if (nval_f > 0) then
          nl=nval_f/5
          if (mod(nval_f,5) /= 0) nl=nl+1
          allocate(rvalues(nval_f))
          rvalues=0.0

          j=1
          do i=n_ini+ntext+2, n_end
             read(unit=filevar(i),fmt='(5E16.8)') vet
             if (i /= n_end) then
                rvalues(j:j+4)=vet
                j=j+5
             else
                k=nval_f-(nl-1)*5
                rvalues(j:j+k-1)=vet(1:k)
             end if
          end do
       end if

       return
    End Subroutine Read_F_KeyType

    !!--++
    !!--++ Subroutine Read_I_KeyType(filevar, n_ini, n_end)
    !!--++    character(len=*), dimension(:), intent(in)  :: filevar   ! Input information
    !!--++    integer,                        intent(in)  :: n_ini     ! First line
    !!--++    integer,                        intent(in)  :: n_end     ! Last line to be read
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Read_I_KeyType(filevar, n_ini, n_end)
       !---- Arguments----!
       character(len=*), dimension(:), intent(in)  :: filevar
       integer,                        intent(in)  :: n_ini
       integer,                        intent(in)  :: n_end

       !---- Local Variables ----!
       character(len=80)      :: line
       integer                :: i,j,k, nl
       integer, dimension(10) :: ivet

       ! Init output values
       err_illdata=.false.

       nval_i=0
       ntext=0
       if (allocated(text_ill)) deallocate(text_ill)
       if (allocated(ivalues)) deallocate(ivalues)

       ! Check the correct KeyType I
       line=filevar(n_ini)
       if (line(1:10) /= 'IIIIIIIIII') then
          err_illdata=.true.
          err_illdata_mess=' A bad Block Type has been found'
          return
       end if

       ! Getting information
       read(unit=filevar(n_ini+1),fmt='(10i8)') ivet
       nval_i =ivet(1)
       ntext=ivet(2)

       if (ntext > 0) then
          allocate(text_ill(ntext))
          text_ill=' '
          do i=1,ntext
             j=(n_ini+1)+i
             text_ill(i)=filevar(j)
          end do
       end if

       if (nval_i > 0) then
          nl=nval_i/10
          if (mod(nval_i,10) /= 0) nl=nl+1
          allocate(ivalues(nval_i))
          ivalues=0

          j=1
          do i=n_ini+ntext+2, n_end
             read(unit=filevar(i),fmt='(10i8)') ivet
             if (i /= n_end) then
                ivalues(j:j+9)=ivet
                j=j+10
             else
                k=nval_i-(nl-1)*10
                ivalues(j:j+k-1)=ivet(1:k)
             end if
          end do
       end if

       return
    End Subroutine Read_I_KeyType

    !!--++
    !!--++ Subroutine Read_J_KeyType(filevar, n_ini, n_end)
    !!--++    character(len=*), dimension(:), intent(in)  :: filevar   ! Input information
    !!--++    integer,                        intent(in)  :: n_ini     ! First line
    !!--++    integer,                        intent(in)  :: n_end     ! Last line to be read
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Read_J_KeyType(filevar, n_ini, n_end)
       !---- Arguments----!
       character(len=*), dimension(:), intent(in)  :: filevar
       integer,                        intent(in)  :: n_ini
       integer,                        intent(in)  :: n_end

       !---- Local Variables ----!
       character(len=80)      :: line
       integer                :: i,j,k,nl
       integer, dimension(8)  :: ivet

       ! Init output values
       err_illdata=.false.

       nval_i=0
       ntext=0
       if (allocated(text_ill)) deallocate(text_ill)
       if (allocated(ivalues)) deallocate(ivalues)

       ! Check the correct KeyType J
       line=filevar(n_ini)
       if (line(1:10) /= 'JJJJJJJJJJ') then
          err_illdata=.true.
          err_illdata_mess=' A bad Block Type has been found'
          return
       end if

       ! Getting information
       read(unit=filevar(n_ini+1),fmt='(10i8)') ivet
       nval_i =ivet(1)
       ntext=ivet(2)

       if (ntext > 0) then
          allocate(text_ill(ntext))
          text_ill=' '
          do i=1,ntext
             j=(n_ini+1)+i
             text_ill(i)=filevar(j)
          end do
       end if

       if (nval_i > 0) then
          nl=nval_i/8
          if (mod(nval_i,8) /= 0) nl=nl+1
          allocate(ivalues(nval_i))
          ivalues=0

          j=1
          do i=n_ini+ntext+2, n_end
             read(unit=filevar(i),fmt='(8i10)') ivet
             if (i /= n_end) then
                ivalues(j:j+7)=ivet
                j=j+8
             else
                k=nval_i-(nl-1)*8
                ivalues(j:j+k-1)=ivet(1:k)
             end if
          end do

       end if

       return
    End Subroutine Read_J_KeyType

    !!----
    !!---- Subroutine Read_Numor_D1A(filevar,N)
    !!----    character(len=*),        intent(in)   :: fileinfo
    !!----    type(powder_numor_type), intent(out) :: n
    !!----
    !!---- Subroutine to read a Numor of D1A Instrument at ILL
    !!----
    !!---- Update: 15/03/2011
    !!
    Subroutine Read_Numor_D1A(fileinfo,N)
       !---- Arguments ----!
       character(len=*),         intent(in)   :: fileinfo
       type(powder_numor_type), intent(out)   :: n

       !---- Local Variables ----!
       character(len=80), dimension(:), allocatable :: filevar
       character(len=80)                            :: line
       integer                                      :: nlines
       integer                                      :: numor,idum
       integer                                      :: i,j

       err_illdata=.false.

       ! Detecting numor
       call Number_Lines(fileinfo,nlines)
       if (nlines <=0) then
          err_illdata=.true.
          err_illdata_mess=' Problems trying to read the Numor for D1A Instrument'
          return
       end if

       ! Allocating variables
       if (allocated(filevar)) deallocate(filevar)
       allocate(filevar(nlines))
       call Reading_Lines(fileinfo,nlines,filevar)

       ! Check format for D1B
       call Number_KeyTypes_on_File(filevar,nlines)
       if (.not. equal_vector(n_keytypes,(/1,2,1,1,2,0,0/),7)) then
          err_illdata=.true.
          err_illdata_mess='This numor does not correspond with D1A Format'
          return
       end if

       ! Defining the different blocks and load information on nl_keytypes
       call Set_KeyTypes_on_File(filevar,nlines)

       ! Numor
       call read_R_keyType(filevar,nl_keytypes(1,1,1),nl_keytypes(1,1,2),numor,idum)
       n%numor=numor

       ! Instr/Experimental Name/ Date
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (idum > 0) then
          n%instrm=line(1:4)
          n%header=line(5:32)
       end if
       if (index(u_case(n%instrm),'D1A') <=0) then
          err_illdata=.true.
          err_illdata_mess='This numor does not correspond with D1A Format'
          return
       end if

       ! Title/Sample
       call read_A_keyType(filevar,nl_keytypes(2,2,1),nl_keytypes(2,2,2),idum,line)
       if (idum > 0) then
          n%title=trim(line)
       end if

       ! Control Flags
       call read_I_keyType(filevar,nl_keytypes(5,1,1),nl_keytypes(5,1,2))
       if (nval_i > 0) then
          n%manip=ivalues(4)
          n%nbang=ivalues(5)
          n%nframes=ivalues(7)
          n%nbdata=ivalues(24)
       end if

       ! Real values
       call read_F_keyType(filevar,nl_keytypes(4,1,1),nl_keytypes(4,1,2))
       if (nval_f > 0) then
          n%scantype='2theta'
          n%wave=rvalues(18)
          n%scans(1)=rvalues(36)
          n%scans(2)=rvalues(37)
          n%conditions(1:5)=rvalues(46:50)  ! Temp-s, Temp-r, Temp-sample
       end if

       ! Allocating
       if (allocated(n%counts)) deallocate(n%counts)
       allocate(n%counts(n%nbdata,n%nframes))
       n%counts=0.0

       if (allocated(n%tmc_ang)) deallocate(n%tmc_ang)
       allocate(n%tmc_ang(n%nbang+2,n%nframes))
       n%tmc_ang=0.0

       ! Counts
       call read_I_keyType(filevar,nl_keytypes(5,2,1),nl_keytypes(5,2,2))
       if (nval_i > 0) then
          if (nval_i /= (n%nbdata+n%nbang+2)*n%nframes) then
             err_illdata=.true.
             err_illdata_mess='Counts problems in Numor format for D1A Instrument'
             return
          end if
          i=0
          do j=1,n%nframes
             n%tmc_ang(1,j)=real(ivalues(i+2))*0.001    ! Time (s)
             n%tmc_ang(2,j)=real(ivalues(i+1))          ! Monitor
             n%tmc_ang(3,j)=real(ivalues(i+3))*0.001    ! Angles
             n%counts(:,j)=ivalues(i+5:i+29)
             i=i+29
          end do
       end if

       return
    End Subroutine Read_Numor_D1A

    !!----
    !!---- Subroutine Read_Numor_D1B(filevar,N)
    !!----    character(len=*),        intent(in)   :: fileinfo
    !!----    type(generic_numor_type), intent(out) :: n
    !!----
    !!---- Subroutine to read a Numor of D1B Instrument at ILL
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Read_Numor_D1B(fileinfo,N)
       !---- Arguments ----!
       character(len=*),         intent(in)   :: fileinfo
       type(powder_numor_type), intent(out)   :: n

       !---- Local Variables ----!
       character(len=80), dimension(:), allocatable :: filevar
       character(len=80)                            :: line
       integer                                      :: nlines
       integer                                      :: numor,idum

       err_illdata=.false.

       ! Detecting numor
       call Number_Lines(fileinfo,nlines)
       if (nlines <=0) then
          err_illdata=.true.
          err_illdata_mess=' Problems trying to read the Numor for D1B Instrument'
          return
       end if

       ! Allocating variables
       if (allocated(filevar)) deallocate(filevar)
       allocate(filevar(nlines))
       call Reading_Lines(fileinfo,nlines,filevar)

       ! Check format for D1B
       call Number_KeyTypes_on_File(filevar,nlines)
       if (.not. equal_vector(n_keytypes,(/1,2,1,1,2,0,0/),7)) then
          err_illdata=.true.
          err_illdata_mess='This numor does not correspond with D1B Format'
          return
       end if

       ! Defining the different blocks and load information on nl_keytypes
       call Set_KeyTypes_on_File(filevar,nlines)

       ! Numor
       call read_R_keyType(filevar,nl_keytypes(1,1,1),nl_keytypes(1,1,2),numor,idum)
       n%numor=numor

       ! Instr/Experimental Name/ Date
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (idum > 0) then
          n%instrm=line(1:4)
          n%header=line(5:32)
       end if

       ! Title/Sample
       call read_A_keyType(filevar,nl_keytypes(2,2,1),nl_keytypes(2,2,2),idum,line)
       if (idum > 0) then
          n%title=trim(line)
       end if

       ! Control Flags
       call read_I_keyType(filevar,nl_keytypes(5,1,1),nl_keytypes(5,1,2))
       if (nval_i > 0) then
       end if

       ! Real values
       call read_F_keyType(filevar,nl_keytypes(4,1,1),nl_keytypes(4,1,2))
       if (nval_f > 0) then
          n%scantype='2theta'
          n%wave=rvalues(18)
          n%scans(1)=rvalues(7)
          n%scans(2)=0.2
          n%conditions(1:3)=rvalues(46:48)  ! Temp-s, Temp-r, Temp-sample
          n%time=rvalues(39)*0.001          ! T en seg
       end if

       ! Counts
       if (allocated(n%counts)) deallocate(n%counts)
       call read_I_keyType(filevar,nl_keytypes(5,2,1),nl_keytypes(5,2,2))
       if (nval_i > 0) then
          n%nbdata=nval_i-3
          n%monitor=real(ivalues(1))
          if (n%nbdata > 0) then
             allocate(n%counts(n%nbdata,1))
             n%counts(:,1)=real(ivalues(4:nval_i))
          end if
       end if

       return
    End Subroutine Read_Numor_D1B

    !!----
    !!---- Subroutine Read_Numor_D2B(filevar,N)
    !!----    character(len=*),        intent(in)   :: fileinfo
    !!----    type(Powder_numor_type), intent(out)  :: n
    !!----
    !!---- Subroutine to read a Numor of D2B Instrument at ILL
    !!----
    !!---- Counts: 128 x 128
    !!----
    !!---- Update: 15/03/2011
    !!
    Subroutine Read_Numor_D2B(fileinfo,N)
       !---- Arguments ----!
       character(len=*),         intent(in)   :: fileinfo
       type(Powder_numor_type),  intent(out)   :: n

       !---- Local Variables ----!
       character(len=80), dimension(:), allocatable :: filevar
       character(len=80)                            :: line
       character(len=5)                             :: car
       integer                                      :: i,nlines
       integer                                      :: numor,idum

       err_illdata=.false.

       ! Detecting numor
       call Number_Lines(fileinfo,nlines)
       if (nlines <=0) then
          err_illdata=.true.
          err_illdata_mess=' Problems trying to read the Numor for D2B Instrument'
          return
       end if

       ! Allocating variables
       if (allocated(filevar)) deallocate(filevar)
       allocate(filevar(nlines))
       call Reading_Lines(fileinfo,nlines,filevar)

       ! Defining the different blocks and load information on nl_keytypes
       call Number_KeyTypes_on_File(filevar,nlines)
       call Set_KeyTypes_on_File(filevar,nlines)

       ! Check format for D2B
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (index(u_case(line(1:4)),'D2B') <= 0) then
          err_illdata=.true.
          err_illdata_mess='This numor does not correspond with D2B Format'
          return
       end if

       ! Numor
       call read_R_keyType(filevar,nl_keytypes(1,1,1),nl_keytypes(1,1,2),numor,idum)
       n%numor=numor

       ! Instr/Experimental Name/ Date
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (idum > 0) then
          n%instrm=line(1:4)
          n%header=line(5:32)
       end if

       ! Title/Sample
       call read_A_keyType(filevar,nl_keytypes(2,2,1),nl_keytypes(2,2,2),idum,line)
       if (idum > 0) then
          n%title=trim(line(1:60))
          n%scantype=trim(line(73:))
       end if

       ! Control Flags
       call read_I_keyType(filevar,nl_keytypes(5,1,1),nl_keytypes(5,1,2))
       if (nval_i > 0) then
          n%manip=ivalues(4)              ! 1: 2Theta, 2: Omega, 3:Chi, 4: Phi
          n%nbang=ivalues(5)              ! Total number of angles moved during scan
          n%nframes=ivalues(7)            ! Frames medidos. En general igual que los solicitados
          n%icalc=ivalues(9)
          n%nbdata=ivalues(24)            ! Number of Points
          n%icdesc(1:7)=ivalues(25:31)
       end if

       ! Real values
       call read_F_keyType(filevar,nl_keytypes(4,1,1),nl_keytypes(4,1,2))
       if (nval_f > 0) then
          n%wave=rvalues(18)           ! Wavelength
          n%scans=rvalues(36:38)       ! Scan start, Scan step, Scan width
          n%conditions=rvalues(46:50)  ! Temp-s, Temp-r, Temp-sample. Voltmeter, Mag.Field
       end if

       ! Allocating
       if (allocated(n%counts)) deallocate(n%counts)
       allocate(n%counts(n%nbdata,n%nframes))
       n%counts=0.0

       if (allocated(n%tmc_ang)) deallocate(n%tmc_ang)
       allocate(n%tmc_ang(n%nbang+3,n%nframes))
       n%tmc_ang=0.0

       ! Loading Frames
       do i=1,n%nframes

          !Time/Monitor/Counts/Angles
          call read_F_keyType(filevar,nl_keytypes(4,i+1,1),nl_keytypes(4,i+1,2))
          if (nval_f > 0 .and. nval_f == (n%nbang+3)) then
             n%tmc_ang(1,i)=rvalues(1)*0.001 ! Time (s)
             n%tmc_ang(2:3,i)=rvalues(2:3)
             n%tmc_ang(4:nval_f,i)=rvalues(4:nval_f)*0.001  ! Angle
          else
             write(unit=car,fmt='(i5)') i
             car=adjustl(car)
             err_illdata=.true.
             err_illdata_mess='Problem reading Time, Monitor, Counts, Angles' &
                               //' parameters in the Frame: '//trim(car)
             return
          end if

          ! Counts
          call read_I_keyType(filevar,nl_keytypes(5,i+1,1),nl_keytypes(5,i+1,2))
          if (nval_i /= n%nbdata) then
             write(unit=car,fmt='(i5)') i
             car=adjustl(car)
             err_illdata=.true.
             err_illdata_mess='Problem reading Counts in the Frame: '//trim(car)
             return
          end if
          if (nval_i > 0) then
             n%counts(:,i)=ivalues(1:n%nbdata)
          end if
       end do

       return
    End Subroutine Read_Numor_D2B

    !!----
    !!---- Subroutine Read_Numor_D9(filevar,N)
    !!----    character(len=*),        intent(in)   :: fileinfo
    !!----    type(SXTAL_numor_type),  intent(out) :: n
    !!----
    !!---- Subroutine to read a Numor of D9 Instrument at ILL
    !!----
    !!---- Counts: 32 x 32 = 1024
    !!----
    !!---- Update: 14/03/2011
    !!
    Subroutine Read_Numor_D9(fileinfo,N)
       !---- Arguments ----!
       character(len=*),         intent(in)   :: fileinfo
       type(SXTAL_numor_type),   intent(out)   :: n

       !---- Local Variables ----!
       character(len=80), dimension(:), allocatable :: filevar
       character(len=80)                            :: line
       character(len=5)                             :: car
       integer                                      :: i,nlines
       integer                                      :: numor,idum

       err_illdata=.false.

       ! Detecting numor
       call Number_Lines(fileinfo,nlines)
       if (nlines <=0) then
          err_illdata=.true.
          err_illdata_mess=' Problems trying to read the Numor for D9 Instrument'
          return
       end if

       ! Allocating variables
       if (allocated(filevar)) deallocate(filevar)
       allocate(filevar(nlines))
       call Reading_Lines(fileinfo,nlines,filevar)

       ! Defining the different blocks and load information on nl_keytypes
       call Number_KeyTypes_on_File(filevar,nlines)
       call Set_KeyTypes_on_File(filevar,nlines)

       ! Check format for D9
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (index(line(1:4),'D9') <= 0) then
          err_illdata=.true.
          err_illdata_mess='This numor does not correspond with D9 Format'
          return
       end if

       ! Numor
       call read_R_keyType(filevar,nl_keytypes(1,1,1),nl_keytypes(1,1,2),numor,idum)
       n%numor=numor

       ! Instr/Experimental Name/ Date
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (idum > 0) then
          n%instrm=line(1:4)
          n%header=line(5:32)
       end if

       ! Title/Sample
       call read_A_keyType(filevar,nl_keytypes(2,2,1),nl_keytypes(2,2,2),idum,line)
       if (idum > 0) then
          n%title=trim(line(1:60))
          n%scantype=trim(line(73:))
       end if

       ! Control Flags
       call read_I_keyType(filevar,nl_keytypes(5,1,1),nl_keytypes(5,1,2))
       if (nval_i > 0) then
          n%manip=ivalues(4)              ! 1: 2Theta, 2: Omega, 3:Chi, 4: Phi
          n%nbang=ivalues(5)              ! Total number of angles moved during scan
          n%nframes=ivalues(7)            ! Frames medidos. En general igual que los solicitados
          n%icalc=ivalues(9)
          n%nbdata=ivalues(24)            ! Number of Points
          n%icdesc(1:7)=ivalues(25:31)
       end if

       ! Real values
       call read_F_keyType(filevar,nl_keytypes(4,1,1),nl_keytypes(4,1,2))
       if (nval_f > 0) then
          n%HMin=rvalues(1:3)          ! HKL min
          n%angles=rvalues(4:8)        ! Phi, Chi, Omega, 2Theta, Psi
          n%ub(1,:)=rvalues(9:11)      !
          n%ub(2,:)=rvalues(12:14)     ! UB Matrix
          n%ub(3,:)=rvalues(15:17)     !
          n%wave=rvalues(18)           ! Wavelength
          n%HMax=rvalues(22:24)        ! HKL max
          n%dh=rvalues(25:27)          ! Delta HKL
          n%scans=rvalues(36:38)       ! Scan start, Scan step, Scan width
          n%preset=rvalues(39)         ! Preset
          n%cpl_fact=rvalues(43)       ! Coupling factor
          n%conditions=rvalues(46:50)  ! Temp-s, Temp-r, Temp-sample. Voltmeter, Mag.Field
       end if

       ! Allocating
       if (allocated(n%counts)) deallocate(n%counts)
       allocate(n%counts(n%nbdata,n%nframes))
       n%counts=0.0

       if (allocated(n%tmc_ang)) deallocate(n%tmc_ang)
       allocate(n%tmc_ang(n%nbang+3,n%nframes))
       n%tmc_ang=0.0

       ! Loading Frames
       do i=1,n%nframes

          !Time/Monitor/Counts/Angles
          call read_F_keyType(filevar,nl_keytypes(4,i+1,1),nl_keytypes(4,i+1,2))
          if (nval_f > 0) then
             select case (nval_f)
                case (4)
                   n%tmc_ang(1,i)=rvalues(1)*0.001 ! Time (s)
                   n%tmc_ang(2:3,i)=rvalues(2:3)
                   n%tmc_ang(4,i)=rvalues(4)*0.001  ! Angle

                case (5)
                   n%tmc_ang(1,i)=rvalues(1)*0.001 ! Time (s)
                   n%tmc_ang(2:3,i)=rvalues(2:3)
                   n%tmc_ang(4:5,i)=rvalues(4:5)*0.001  ! Angle

                case default
                   write(unit=car,fmt='(i5)') i
                   car=adjustl(car)
                   err_illdata=.true.
                   err_illdata_mess='Problem reading Time, Monitor, Counts, Angles' &
                                    //' parameters in the Frame: '//trim(car)
                   return
             end select
          end if

          ! Counts
          call read_I_keyType(filevar,nl_keytypes(5,i+1,1),nl_keytypes(5,i+1,2))
          if (nval_i /= n%nbdata) then
             write(unit=car,fmt='(i5)') i
             car=adjustl(car)
             err_illdata=.true.
             err_illdata_mess='Problem reading Counts in the Frame: '//trim(car)
             return
          end if
          if (nval_i > 0) then
             n%counts(:,i)=ivalues(1:n%nbdata)
          end if
       end do

       return
    End Subroutine Read_Numor_D9

    !!----
    !!---- Subroutine Read_Numor_D19(filevar,N)
    !!----    character(len=*),        intent(in)   :: fileinfo
    !!----    type(SXTAL_numor_type),  intent(out) :: n
    !!----
    !!---- Subroutine to read a Numor of D19 Instrument at ILL
    !!----
    !!---- Counts: 640 x 256 = 1024
    !!----
    !!---- Update: 14/03/2011
    !!
    Subroutine Read_Numor_D19(fileinfo,N)
       !---- Arguments ----!
       character(len=*),         intent(in)   :: fileinfo
       type(SXTAL_numor_type),   intent(out)   :: n

       !---- Local Variables ----!
       character(len=80), dimension(:), allocatable :: filevar
       character(len=80)                            :: line
       character(len=5)                             :: car
       integer                                      :: i,nlines
       integer                                      :: numor,idum

       err_illdata=.false.

       ! Detecting numor
       call Number_Lines(fileinfo,nlines)
       if (nlines <=0) then
          err_illdata=.true.
          err_illdata_mess=' Problems trying to read the Numor for D19 Instrument'
          return
       end if

       ! Allocating variables
       if (allocated(filevar)) deallocate(filevar)
       allocate(filevar(nlines))
       call Reading_Lines(fileinfo,nlines,filevar)

       ! Defining the different blocks and load information on nl_keytypes
       call Number_KeyTypes_on_File(filevar,nlines)
       call Set_KeyTypes_on_File(filevar,nlines)

       ! Check format for D19
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (index(line(1:4),'D19') <= 0) then
          err_illdata=.true.
          err_illdata_mess='This numor does not correspond with D19 Format'
          return
       end if

       ! Numor
       call read_R_keyType(filevar,nl_keytypes(1,1,1),nl_keytypes(1,1,2),numor,idum)
       n%numor=numor

       ! Instr/Experimental Name/ Date
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (idum > 0) then
          n%instrm=line(1:4)
          n%header=line(5:32)
       end if

       ! Title/Sample
       call read_A_keyType(filevar,nl_keytypes(2,2,1),nl_keytypes(2,2,2),idum,line)
       if (idum > 0) then
          n%title=trim(line(1:60))
          n%scantype=trim(line(73:))
       end if

       ! Control Flags
       call read_I_keyType(filevar,nl_keytypes(5,1,1),nl_keytypes(5,1,2))
       if (nval_i > 0) then
          n%manip=ivalues(4)              ! 1: 2Theta, 2: Omega, 3:Chi, 4: Phi
          n%nbang=ivalues(5)              ! Total number of angles moved during scan
          n%nframes=ivalues(7)            ! Frames medidos. En general igual que los solicitados
          n%icalc=ivalues(9)
          n%nbdata=ivalues(24)            ! Number of Points
          n%icdesc(1:7)=ivalues(25:31)
       end if

       ! Real values
       call read_F_keyType(filevar,nl_keytypes(4,1,1),nl_keytypes(4,1,2))
       if (nval_f > 0) then
          n%HMin=rvalues(1:3)          ! HKL min
          n%angles=rvalues(4:8)        ! Phi, Chi, Omega, 2Theta, Psi
          n%ub(1,:)=rvalues(9:11)      !
          n%ub(2,:)=rvalues(12:14)     ! UB Matrix
          n%ub(3,:)=rvalues(15:17)     !
          n%wave=rvalues(18)           ! Wavelength
          n%HMax=rvalues(22:24)        ! HKL max
          n%dh=rvalues(25:27)          ! Delta HKL
          n%scans=rvalues(36:38)       ! Scan start, Scan step, Scan width
          n%preset=rvalues(39)         ! Preset
          n%cpl_fact=rvalues(43)       ! Coupling factor
          n%conditions=rvalues(46:50)  ! Temp-s, Temp-r, Temp-sample. Voltmeter, Mag.Field
       end if

       ! Allocating
       if (allocated(n%counts)) deallocate(n%counts)
       allocate(n%counts(n%nbdata,n%nframes))
       n%counts=0.0

       if (allocated(n%tmc_ang)) deallocate(n%tmc_ang)
       allocate(n%tmc_ang(n%nbang+3,n%nframes))
       n%tmc_ang=0.0

       ! Loading Frames
       do i=1,n%nframes

          !Time/Monitor/Counts/Angles
          call read_F_keyType(filevar,nl_keytypes(4,i+1,1),nl_keytypes(4,i+1,2))
          if (nval_f > 0 .and. nval_f == (n%nbang+3)) then
             n%tmc_ang(1,i)=rvalues(1)*0.001 ! Time (s)
             n%tmc_ang(2:3,i)=rvalues(2:3)
             n%tmc_ang(4:nval_f,i)=rvalues(4:nval_f)*0.001  ! Angle
          else
             write(unit=car,fmt='(i5)') i
             car=adjustl(car)
             err_illdata=.true.
             err_illdata_mess='Problem reading Time, Monitor, Counts, Angles' &
                               //' parameters in the Frame: '//trim(car)
             return
          end if

          ! Counts
          call read_I_keyType(filevar,nl_keytypes(5,i+1,1),nl_keytypes(5,i+1,2))
          if (nval_i /= n%nbdata) then
             write(unit=car,fmt='(i5)') i
             car=adjustl(car)
             err_illdata=.true.
             err_illdata_mess='Problem reading Counts in the Frame: '//trim(car)
             return
          end if
          if (nval_i > 0) then
             n%counts(:,i)=ivalues(1:n%nbdata)
          end if
       end do

       return
    End Subroutine Read_Numor_D19

    !!----
    !!---- Subroutine Read_Numor_Generic(filevar,N)
    !!----    character(len=*), intent(in) :: fileinfo
    !!----    type(generic_numor_type), intent(out) :: n
    !!----
    !!---- Read a Numor (Not yet operative)
    !!----
    !!---- Update: 11/03/2011
    !!
    Subroutine Read_Numor_Generic(fileinfo,N)
       !---- Arguments ----!
       character(len=*), intent(in) :: fileinfo
       type(generic_numor_type), intent(out) :: n

       !---- Local Variables ----!
       character(len=80), dimension(:), allocatable :: filevar
       character(len=40), dimension(5)              :: dire
       character(len=1024)                          :: line
       character(len=80)                            :: linec
       integer                                      :: nlines
       integer                                      :: i,j,numor,idum,nl

       err_illdata=.false.

       ! Detecting numor
       call Number_Lines(fileinfo,nlines)
       if (nlines <=0) then
          err_illdata=.true.
          err_illdata_mess=' Problems trying to read the numor for D20 Instrument'
          return
       end if

       ! Allocating variables
       if (allocated(filevar)) deallocate(filevar)
       allocate(filevar(nlines))
       call Reading_Lines(fileinfo,nlines,filevar)

       ! Check format for D20
       call Number_KeyTypes_on_File(filevar,nlines)
       if (.not. equal_vector(n_keytypes,(/1,2,1,6,0,1,0/),7)) then
          err_illdata=.true.
          err_illdata_mess='This numor does not correspond with D20 Format'
          return
       end if

       ! Defining the different blocks
       call Set_KeyTypes_on_File(filevar,nlines)

       ! Numor
       call read_R_keyType(filevar,nl_keytypes(1,1,1),nl_keytypes(1,1,2),numor,idum)
       n%numor=numor

       ! Instr/Experimental Name/ Date
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (idum > 0) then
          n%instr=line(1:4)
          n%expname=line(5:14)
          n%date=line(15:32)
       end if

       ! Title/Sample
       if (allocated(n%sampleid%namevar)) deallocate(n%sampleid%namevar)
       if (allocated(n%sampleid%cvalues)) deallocate(n%sampleid%cvalues)

       call read_A_keyType(filevar,nl_keytypes(2,2,1),nl_keytypes(2,2,2),idum,line)
       select case (idum)
          case (1:80)
             n%SampleID%n=4
             allocate(n%sampleid%namevar(4))
             allocate(n%sampleid%cvalues(4))

             ! Given names
             n%sampleid%namevar=' '
             n%sampleid%namevar(1)='Name'
             n%sampleid%namevar(2)='Local Contact'
             n%sampleid%namevar(3)='Experimentalist'
             n%sampleid%namevar(4)='Proposal Number'

             ! Load values
             n%sampleid%cvalues=' '
             n%title=trim(line)
             n%sampleid%cvalues(1)=trim(line)

          case (81:)
             nl=(idum/80)+mod(idum,80)
             n%SampleID%n=nl
             allocate(n%sampleid%namevar(nl))
             allocate(n%sampleid%cvalues(nl))

             ! Give names
             do i=1,nl
                j=index(line,char(9))
                linec=line(1:j-1)
                line=line(j+1:)

                j=index(linec,':')
                n%sampleid%namevar(i)=linec(1:j-1)
                n%sampleid%cvalues(i)=trim(linec(j+1:))
             end do
       end select

       ! Diffractometers Optics / Reactor parameters
       call read_F_keyType(filevar,nl_keytypes(4,1,1),nl_keytypes(4,1,2))
       if (nval_f > 0) then
          n%diffopt%n=nval_f
          if (allocated(n%diffopt%namevar)) deallocate(n%diffopt%namevar)
          if (allocated(n%diffopt%rvalues)) deallocate(n%diffopt%rvalues)
          allocate(n%diffopt%namevar(nval_f))
          allocate(n%diffopt%rvalues(nval_f))

          ! Given names
          n%diffopt%namevar=' '
          j=0
          do i=2,ntext
             read(unit=text_ill(i),fmt='(5a16)') dire
             do nl=1,5
                j=j+1
                n%diffopt%namevar(j)=dire(nl)
             end do
          end do

          ! Load values
          n%diffopt%rvalues=rvalues(1:nval_f)
       end if

       ! Monochromator Motor Parameters
       call read_F_keyType(filevar,nl_keytypes(4,2,1),nl_keytypes(4,2,2))
       if (nval_f > 0) then
          n%monmpar%n=nval_f
          if (allocated(n%monmpar%namevar)) deallocate(n%monmpar%namevar)
          if (allocated(n%monmpar%rvalues)) deallocate(n%monmpar%rvalues)
          allocate(n%monmpar%namevar(nval_f))
          allocate(n%monmpar%rvalues(nval_f))

          ! Given names
          n%monmpar%namevar=' '
          j=0
          do i=2,ntext
             read(unit=text_ill(i),fmt='(5a16)') dire
             do nl=1,5
                j=j+1
                n%monmpar%namevar(j)=dire(nl)
             end do
          end do

          ! Load values
          n%monmpar%rvalues=rvalues(1:nval_f)
       end if

       ! Diffractometer Motor Parameters
       call read_F_keyType(filevar,nl_keytypes(4,3,1),nl_keytypes(4,3,2))
       if (nval_f > 0) then
          n%diffmpar%n=nval_f
          if (allocated(n%diffmpar%namevar)) deallocate(n%diffmpar%namevar)
          if (allocated(n%diffmpar%rvalues)) deallocate(n%diffmpar%rvalues)
          allocate(n%diffmpar%namevar(nval_f))
          allocate(n%diffmpar%rvalues(nval_f))

          ! Given names
          n%diffmpar%namevar=' '
          j=0
          do i=2,ntext
             read(unit=text_ill(i),fmt='(5a16)') dire
             do nl=1,5
                j=j+1
                n%diffmpar%namevar(j)=dire(nl)
             end do
          end do

          ! Load values
          n%diffmpar%rvalues=rvalues(1:nval_f)
       end if

       ! Detector and DAS Parameters
       call read_F_keyType(filevar,nl_keytypes(4,4,1),nl_keytypes(4,4,2))
       if (nval_f > 0) then
          n%detpar%n=nval_f
          if (allocated(n%detpar%namevar)) deallocate(n%detpar%namevar)
          if (allocated(n%detpar%rvalues)) deallocate(n%detpar%rvalues)
          allocate(n%detpar%namevar(nval_f))
          allocate(n%detpar%rvalues(nval_f))

          ! Given names
          n%detpar%namevar=' '
          j=0
          do i=2,ntext
             read(unit=text_ill(i),fmt='(5a16)') dire
             do nl=1,5
                j=j+1
                n%detpar%namevar(j)=dire(nl)
             end do
          end do

          ! Load values
          n%detpar%rvalues=rvalues(1:nval_f)
       end if

       ! Data Acquisition Control Parameters
       call read_F_keyType(filevar,nl_keytypes(4,5,1),nl_keytypes(4,5,2))
       if (nval_f > 0) then
          n%dacparam%n=nval_f
          if (allocated(n%dacparam%namevar)) deallocate(n%dacparam%namevar)
          if (allocated(n%dacparam%rvalues)) deallocate(n%dacparam%rvalues)
          allocate(n%dacparam%namevar(nval_f))
          allocate(n%dacparam%rvalues(nval_f))

          ! Given names
          n%dacparam%namevar=' '
          j=0
          do i=2,ntext
             read(unit=text_ill(i),fmt='(5a16)') dire
             do nl=1,5
                j=j+1
                n%dacparam%namevar(j)=dire(nl)
             end do
          end do

          ! Load values
          n%dacparam%rvalues=rvalues(1:nval_f)
       end if

       ! Sample Status
       call read_F_keyType(filevar,nl_keytypes(4,6,1),nl_keytypes(4,6,2))
       if (nval_f > 0) then
          n%samplest%n=nval_f
          if (allocated(n%samplest%namevar)) deallocate(n%samplest%namevar)
          if (allocated(n%samplest%rvalues)) deallocate(n%samplest%rvalues)
          allocate(n%samplest%namevar(nval_f))
          allocate(n%samplest%rvalues(nval_f))

          ! Given names
          n%samplest%namevar=' '
          j=0
          do i=2,ntext
             read(unit=text_ill(i),fmt='(5a16)') dire
             do nl=1,5
                j=j+1
                n%samplest%namevar(j)=dire(nl)
             end do
          end do

          ! Load values
          n%samplest%rvalues=rvalues(1:nval_f)
       end if

       ! Counts Information
       call read_J_keyType(filevar,nl_keytypes(6,1,1),nl_keytypes(6,1,2))

       if (nval_i > 0) then
          n%icounts%n=nval_i
          if (allocated(n%icounts%namevar)) deallocate(n%icounts%namevar)
          if (allocated(n%icounts%ivalues)) deallocate(n%icounts%ivalues)
          allocate(n%icounts%namevar(1))
          allocate(n%icounts%ivalues(nval_i))

          ! Given names
          n%icounts%namevar(1)=' N. of Counts'

          ! Load values
          n%icounts%ivalues=ivalues(1:nval_i)
       end if

       return
    End Subroutine Read_Numor_Generic

    !!----
    !!---- Subroutine Read_Numor_D20(filevar,N)
    !!----    character(len=*), intent(in) :: fileinfo
    !!----    type(generic_numor_type), intent(out) :: n
    !!----
    !!---- Read a Numor for D20 Instrument at ILL
    !!----
    !!---- Update: 11/03/2011
    !!
    Subroutine Read_Numor_D20(fileinfo,N)
       !---- Arguments ----!
       character(len=*), intent(in) :: fileinfo
       type(powder_numor_type), intent(out) :: n

       !---- Local Variables ----!
       character(len=80), dimension(:), allocatable :: filevar
       character(len=1024)                          :: line
       integer                                      :: nlines
       integer                                      :: i,j,numor,idum

       err_illdata=.false.

       ! Detecting numor
       call Number_Lines(fileinfo,nlines)
       if (nlines <=0) then
          err_illdata=.true.
          err_illdata_mess=' Problems trying to read the numor for D20 Instrument'
          return
       end if

       ! Allocating variables
       if (allocated(filevar)) deallocate(filevar)
       allocate(filevar(nlines))
       call Reading_Lines(fileinfo,nlines,filevar)

       ! Check format for D20
       call Number_KeyTypes_on_File(filevar,nlines)
       if (.not. equal_vector(n_keytypes,(/1,2,1,6,0,1,0/),7)) then
          err_illdata=.true.
          err_illdata_mess='This numor does not correspond with D20 Format'
          return
       end if

       ! Defining the different blocks
       call Set_KeyTypes_on_File(filevar,nlines)

       ! Numor
       call read_R_keyType(filevar,nl_keytypes(1,1,1),nl_keytypes(1,1,2),numor,idum)
       n%numor=numor

       ! Instr/Experimental Name/ Date
       call read_A_keyType(filevar,nl_keytypes(2,1,1),nl_keytypes(2,1,2),idum,line)
       if (idum > 0) then
          n%instrm=line(1:4)
          n%header=line(5:32)
       end if

       ! Title
       call read_A_keyType(filevar,nl_keytypes(2,2,1),nl_keytypes(2,2,2),idum,line)
       if (idum > 0) then
          i=index(line,'Title :')
          j=index(line,'Local')
          if (i > 0) then
             n%title=line(i+7:j-2)
             n%title=trim(n%title)
          end if
       end if

       ! Wave
       call read_F_keyType(filevar,nl_keytypes(4,1,1),nl_keytypes(4,1,2))
       if (nval_f > 0) then
          n%wave=rvalues(4)
       end if

       ! Monochromator Motor Parameters
       !call read_F_keyType(filevar,nl_keytypes(4,2,1),nl_keytypes(4,2,2))
       !if (nval_f > 0) then
       !end if

       ! Diffractometer Motor Parameters
       !call read_F_keyType(filevar,nl_keytypes(4,3,1),nl_keytypes(4,3,2))
       !if (nval_f > 0) then
       !   ! Zero 2Theta?
       !end if

       ! Detector and DAS Parameters
       !call read_F_keyType(filevar,nl_keytypes(4,4,1),nl_keytypes(4,4,2))
       !if (nval_f > 0) then
       !end if

       ! Data Acquisition Control Parameters
       call read_F_keyType(filevar,nl_keytypes(4,5,1),nl_keytypes(4,5,2))
       if (nval_f > 0) then
          n%monitor=rvalues(51)
          n%time=rvalues(52)
       end if

       ! Sample Status
       call read_F_keyType(filevar,nl_keytypes(4,6,1),nl_keytypes(4,6,2))
       if (nval_f > 0) then
          n%conditions(1:3)=rvalues(2:4)
          n%ScanType='2theta'
          n%scans(1)=rvalues(5)
          n%scans(2)=0.1
       end if

       ! S Block

       ! Counts Information
       call read_J_keyType(filevar,nl_keytypes(6,1,1),nl_keytypes(6,1,2))
       if (nval_i > 0) then
          n%nbdata=nval_i
          if (n%nbdata > 0) then
             allocate(n%counts(n%nbdata,1))
             n%counts(:,1)=real(ivalues(1:nval_i))
          end if
       end if

       return
    End Subroutine Read_Numor_D20

    !!----
    !!---- Subroutine Read_Powder_Numor(PathNumor,Instrument,Num)
    !!----    character(len=*),            intent(in)    :: PathNumor
    !!----    character(len=*),            intent(in)    :: Instrument
    !!----    type(Powder_Numor_type),     intent(out)   :: Num
    !!----
    !!----    Read a Powder numor from the ILL database.
    !!----    In case of error the subroutine puts ERR_ILLData=.true.
    !!----    and fils the error message variable ERR_ILLData_Mess.
    !!----
    !!---- Update: 15/03/2011
    !!
    Subroutine Read_Powder_Numor(PathNumor,Instrument,Num)
       !---- Arguments ----!
       character(len=*),           intent(in)    :: PathNumor
       character(len=*),           intent(in)    :: Instrument
       type(Powder_Numor_Type),    intent(out)   :: Num

       !---- Local variables ----!
       character(len=512)     :: Path
       character(len=80)      :: Filename
       character(len=4)       :: Instr
       integer                :: n

       ! Initialize
       ERR_ILLData=.false.
       ERR_ILLData_Mess= ' '

       ! Check
       if (len_trim(PathNumor) <= 0 .or. len_trim(Instrument) <=0 ) return
       instr=u_case(adjustl(Instrument))

       ! Path + Numor
       n=index(pathnumor,ops_sep, back=.true.)
       if (n > 0) then
          path=pathnumor(:n)
          filename=pathnumor(n+1:)
       else
          path=' '
          filename=trim(pathnumor)
       end if

       ! Compressed Numor
       n=index(filename,'.Z')
       if (n > 0) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess= " Numor file is compressed. Please uncompress the numor before to use this routine"
          return
       end if

       ! Read Numor
       select case (trim(u_case(Instr)))
          case ('D1A')
             call Read_Numor_D1A(trim(path)//trim(filename),Num)

          case ('D1B')
             call Read_Numor_D1B(trim(path)//trim(filename),Num)

          case ('D2B')
             call Read_Numor_D2B(trim(path)//trim(filename),Num)

          case ('D20')
             call Read_Numor_D20(trim(path)//trim(filename),Num)

          case default
             ERR_ILLData=.true.
             ERR_ILLData_Mess= " Not Implemented for the Powder Instrument name: "//trim(instrument)
             return
       end select

       return
    End Subroutine Read_Powder_Numor

    !!----
    !!---- Subroutine Read_SXTAL_Numor(PathNumor,Instrument,Num)
    !!----    character(len=*),            intent(in)    :: PathNumor
    !!----    character(len=*),            intent(in)    :: Instrument
    !!----    type(SXTAL_Numor_type),      intent(out)   :: Num
    !!----
    !!----    Read a SXTAL numor from the ILL database.
    !!----    In case of error the subroutine puts ERR_ILLData=.true.
    !!----    and fils the error message variable ERR_ILLData_Mess.
    !!----
    !!---- Update: 15/03/2011 9:16:46
    !!
    Subroutine Read_SXTAL_Numor(PathNumor,Instrument,Num)
       !---- Arguments ----!
       character(len=*),           intent(in)    :: PathNumor
       character(len=*),           intent(in)    :: Instrument
       type(SXTAL_Numor_Type),     intent(out)   :: Num

       !---- Local variables ----!
       character(len=512)     :: Path
       character(len=80)      :: Filename
       character(len=4)       :: Instr
       integer                :: n

       ! Initialize
       ERR_ILLData=.false.
       ERR_ILLData_Mess= ' '

       ! Check
       if (len_trim(PathNumor) <= 0 .or. len_trim(Instrument) <=0 ) return
       instr=u_case(adjustl(Instrument))

       ! Path + Numor
       n=index(pathnumor,ops_sep, back=.true.)
       if (n > 0) then
          path=pathnumor(:n)
          filename=pathnumor(n+1:)
       else
          path=' '
          filename=trim(pathnumor)
       end if

       ! Compressed Numor
       n=index(filename,'.Z')
       if (n > 0) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess= " Numor file is compressed. Please uncompress the numor before to use this routine"
          return
       end if

       ! Read Numor
       select case (trim(Instr))
          case ('D9')
             call Read_Numor_D9(trim(path)//trim(filename),Num)

          case ('D19')
             call Read_Numor_D19(trim(path)//trim(filename),Num)

          case default
             ERR_ILLData=.true.
             ERR_ILLData_Mess= " Not Implemented for the SXTAL Instrument name: "//trim(instrument)
             return
       end select

       return
    End Subroutine Read_SXTAL_Numor

    !!--++
    !!--++ Subroutine Read_R_KeyType(filevar, n_ini, n_end, nrun, nvers)
    !!--++    character(len=*), dimension(:), intent(in)  :: filevar   ! Input information
    !!--++    integer,                        intent(in)  :: n_ini     ! First line
    !!--++    integer,                        intent(in)  :: n_end     ! Last line to be read
    !!--++    integer,                        intent(out) :: nrun      ! Run number of the data
    !!--++    integer,                        intent(out) :: nvers     ! Version data
    !!--++
    !!--++ (Private)
    !!--++ Read the Blocktype R on Numor at ILL.
    !!--++ If text is readen, then it is saved on Text_ILL variable and the number of lines is
    !!--++ saved on ntext
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Read_R_KeyType(filevar, n_ini, n_end, nrun, nvers)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)  :: filevar
       integer,                        intent(in)  :: n_ini
       integer,                        intent(in)  :: n_end
       integer,                        intent(out) :: nrun
       integer,                        intent(out) :: nvers

       !---- Local Variables ----!
       character(len=80)      :: line
       integer                :: i,j
       integer, dimension(10) :: ivet

       ! Init output values
       err_illdata=.false.

       nrun =0
       nvers=0
       ntext=0
       if (allocated(text_ill)) deallocate(text_ill)

       ! Check the correct KeyType R
       line=filevar(n_ini)
       if (line(1:10) /= 'RRRRRRRRRR') then
          err_illdata=.true.
          err_illdata_mess=' A bad Block Type has been found'
          return
       end if

       ! Getting information
       read(unit=filevar(n_ini+1),fmt='(10i8)') ivet
       nrun =ivet(1)   ! Numor
       ntext=ivet(2)   ! Number of lines of descriptive text
       nvers=ivet(3)   ! Version of the data

       if (ntext > 0) then
          allocate(text_ill(ntext))
          text_ill=' '
          do i=1,ntext
             j=(n_ini+1)+i
             if (j > n_end) then
                err_illdata=.true.
                err_illdata_mess='  Impossible to read a line for this block!'
                exit
             end if
             text_ill(i)=trim(filevar(j))
          end do
       end if

       return
    End Subroutine Read_R_KeyType

    !!--++
    !!--++ Subroutine Read_S_KeyType(filevar, n_ini, n_end, ispec, nrest, ntot, nrun, npars)
    !!--++    character(len=*), dimension(:), intent(in)  :: filevar   ! Input information
    !!--++    integer,                        intent(in)  :: n_ini     ! First line
    !!--++    integer,                        intent(in)  :: n_end     ! Last line to be read
    !!--++    integer,                        intent(out) :: ispec     ! Sub-spectrum number
    !!--++    integer,                        intent(out) :: nrest     ! Number of subspectra remaining after ispec
    !!--++    integer,                        intent(out) :: ntot      ! Total number of subspectra in the run
    !!--++    integer,                        intent(out) :: nrun      ! Current run number
    !!--++    integer,                        intent(out) :: npars     ! Number of parameters sections
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Read_S_KeyType(filevar, n_ini, n_end, ispec, nrest, ntot, nrun, npars)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)  :: filevar
       integer,                        intent(in)  :: n_ini
       integer,                        intent(in)  :: n_end
       integer,                        intent(out) :: ispec
       integer,                        intent(out) :: nrest
       integer,                        intent(out) :: ntot
       integer,                        intent(out) :: nrun
       integer,                        intent(out) :: npars

       !---- Local Variables ----!
       character(len=80)      :: line
       integer                :: i,j
       integer, dimension(10) :: ivet

       ! Init output variables
       err_illdata=.false.

       ispec=0
       nrest=0
       ntot=0
       nrun=0
       ntext=0
       npars=0
       if (allocated(text_ill)) deallocate(text_ill)

       ! Check the correct KeyType S
       line=filevar(n_ini)
       if (line(1:10) /= 'SSSSSSSSSS') then
          err_illdata=.true.
          err_illdata_mess=' A bad Block Type has been found'
          return
       end if

       ! Getting information
       read(unit=filevar(n_ini+1),fmt='(10i8)') ivet
       ispec=ivet(1)
       nrest=ivet(2)
       ntot =ivet(3)
       nrun =ivet(4)
       ntext=ivet(5)
       npars=ivet(6)

       if (ntext > 0) then
          allocate(text_ill(ntext))
          text_ill=' '
          do i=1,ntext
             j=(n_ini+1)+i
             if (j > n_end) then
                err_illdata=.true.
                err_illdata_mess=' Impossible to read a line for this block!'
                exit
             end if
             text_ill(i)=trim(filevar(j))
          end do
       end if

       return
    End Subroutine Read_S_KeyType

    !!--++
    !!--++ Subroutine Read_V_KeyType(filevar, n_ini, n_end)
    !!--++    character(len=*), dimension(:), intent(in)  :: filevar   ! Input information
    !!--++    integer,                        intent(in)  :: n_ini     ! First line
    !!--++    integer,                        intent(in)  :: n_end     ! Last line to be read
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Read_V_KeyType(filevar, n_ini, n_end)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)  :: filevar
       integer,                        intent(in)  :: n_ini
       integer,                        intent(in)  :: n_end

       !---- Local Variables ----!
       character(len=80)      :: line
       integer                :: i,j

       ! Init output values
       err_illdata=.false.
       ntext=0
       if (allocated(text_ill)) deallocate(text_ill)

       ! Check the correct KeyType V
       line=filevar(n_ini)
       if (line(1:10) /= 'VVVVVVVVVV') then
          err_illdata=.true.
          err_illdata_mess=' A bad Block Type has been found'
          return
       end if

       ! Getting information
       ntext=n_end-n_ini-1
       if (ntext > 0) then
          allocate(text_ill(ntext))
          text_ill=' '
          j=0
          do i=n_ini+1,n_end
             j=j+1
             text_ill(j)=filevar(i)
          end do
       end if

       return
    End Subroutine Read_V_KeyType

    !!----
    !!---- Subroutine Set_Current_Orient(wave,ub,setting)
    !!----   real(kind=cp),                         intent(in)   :: wave
    !!----   real(kind=cp), dimension(3,3),         intent(in)   :: ub
    !!----   real(kind=cp), dimension(3,3),optional,intent(in)   :: setting
    !!----
    !!----    Subroutine setting the Current_Orient global variable
    !!----    If the final UB matrix is singular an error is rised
    !!----    by putting ERR_ILLData=.true. and filling the
    !!----    error message variable ERR_ILLData_Mess.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Set_Current_Orient(wave,ub,setting)
       !---- Argument ----!
       real(kind=cp),                         intent(in)   :: wave
       real(kind=cp), dimension(3,3),         intent(in)   :: ub
       real(kind=cp), dimension(3,3),optional,intent(in)   :: setting

       !--- Local variables ---!
       real(kind=cp)                 :: det
       real(kind=cp), dimension(3,3) :: mat

       Current_Orient%conv=reshape( (/1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0/),(/3,3/))

       if (present(setting)) then
          mat=matmul(setting,ub)  !check that point
          Current_Orient%conv=setting
       else
          mat=ub
       end if

       det=determ_a(mat)
       if (abs(det) < eps) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess="Singular UB or Setting Matrix"
          return
       end if

       Current_Orient%wave=wave
       Current_Orient%ub=mat
       Current_Orient%ubinv=Invert(mat)

       return
    End Subroutine Set_Current_Orient

    !!----
    !!---- Subroutine Set_Default_Instrument(typ,wav)
    !!----    Character(len=*),   optional, intent(in) :: typ  !"Laue" for a Laue diffractometer
    !!----    real, dimension(2), optional, intent(in) :: wav  ! Lambda min and Lambda max
    !!----
    !!----    Construct the Current_Instrument as a default 4C diffractometer
    !!----    (or a Laue diffractometer if Typ and wav are provided)
    !!----    The UB matrix is set to a real matrix corresponding to a measurement
    !!----    done on D9. The characteristics of the diffractometer correspond to
    !!----    those of D9
    !!----
    !!---- Update: May - 2007
    !!
    Subroutine Set_Default_Instrument(typ,wav)
       Character(len=*),   optional, intent(in) :: typ
       real, dimension(2), optional, intent(in) :: wav
       !---- Local Variables ----!
       real(kind=cp)                 :: wave
       integer                       :: npx, npz
       real(kind=cp), dimension(3,3) :: ub

       Current_Instrm%det_offsets=0.0
       Current_Instrm%disp_offsets=0.0
       Current_Instrm%disp_Limits=0.0
       Current_Instrm%ang_Limits=0.0
       Current_Instrm%disp_names=" "
       Current_Instrm%ang_names=" "
       Current_Instrm%e1=(/1.0,0.0,0.0/)
       Current_Instrm%e2=(/0.0,1.0,0.0/)
       Current_Instrm%e3=(/0.0,0.0,1.0/)
       if(present(typ)) then
         Current_Instrm%info= "Default Laue diffractrometer"
         Current_Instrm%name_inst= "LaueDiff"
         Current_Instrm%geom="Laue"
         if(present(wav)) then
           Current_Instrm%wave_min=wav(1)
           Current_Instrm%wave_max=wav(2)
         else
           Current_Instrm%wave_min=0.7
           Current_Instrm%wave_max=5.5
         end if
         Current_Instrm%igeom=3
         wave=2.0
         Current_Instrm%dist_samp_detector=488.0
         Current_Instrm%np_horiz= 2560
         Current_Instrm%np_vert=  1980
         Current_Instrm%horiz= 256.0
         Current_Instrm%vert=  198.0
         npx = 2560
         npz = 1980
         Current_Instrm%agap=0
         Current_Instrm%cgap=0
         Current_Instrm%ang_names(1) ="Gamma"
         Current_Instrm%ang_Limits(1,1:2)=(/2.0,175.0/)
       else
         Current_Instrm%info= "Default 4-cercles diffractrometer"
         Current_Instrm%name_inst= "4C-Diff"
         Current_Instrm%geom="4C-Diff High-Chi, Eulerian cradle"
         Current_Instrm%igeom=2
         wave=0.71
         Current_Instrm%dist_samp_detector=488.0
         Current_Instrm%np_horiz= 32
         Current_Instrm%np_vert= 32
         Current_Instrm%horiz= 64.0
         Current_Instrm%vert=  64.0
         npx = 32
         npz = 32
         Current_Instrm%agap=2
         Current_Instrm%cgap=2
         Current_Instrm%ang_names(1) ="2Theta"
         Current_Instrm%ang_Limits(1,1:2)=(/2.0,130.0/)
       end if
       Current_Instrm%BL_frame="z-up"
       Current_Instrm%dist_units = "mm"
       Current_Instrm%angl_units = "deg"
       Current_Instrm%detector_type = "Flat_rect"
       Current_Instrm%ipsd=2          !Flat detector
       ub=reshape((/-0.0989455,   0.0671905,  -0.1005396, &
                     0.0045075,  -0.1487497,  -0.0642365, &
                    -0.1588914,  -0.0460609,   0.0607861/),(/3,3/))
       Current_Instrm%num_ang=4
       Current_Instrm%num_disp=0
       Current_Instrm%ang_names(2) ="Omega"
       Current_Instrm%ang_names(3) ="Chi"
       Current_Instrm%ang_names(4) ="Phi"
       Current_Instrm%ang_Limits(2,1:2)=(/1.0,49.0/)
       Current_Instrm%ang_Limits(3,1:2)=(/77.0,202.0/)
       Current_Instrm%ang_Limits(4,1:2)=(/-180.0,180.0/)
       if (allocated(Current_Instrm%alphas)) deallocate(Current_Instrm%alphas)
       allocate(Current_Instrm%alphas(npx,npz))
       Current_Instrm%alphas(:,:)=1.0
       call Set_Current_Orient(wave,ub)
       Current_Instrm_set=.true.

       return
    End Subroutine Set_Default_Instrument

    !!----
    !!---- Subroutine Set_ILL_Data_Directory(Filedir)
    !!----    character(len=*), intent(in) :: Filedir  !proposed location of ILL data
    !!----
    !!----    Assign the global public variable: ILL_data_directory
    !!----    If the directory doesn't exist the subroutine rises an error condition
    !!----    by putting ERR_ILLData=.true. and filling the error message
    !!----    variable ERR_ILLData_Mess.
    !!----
    !!---- Update: 08/03/2011
    !!
    Subroutine Set_ILL_Data_Directory(Filedir)
       !---- Arguments ----!
       character(len=*), intent(in) :: Filedir

       !---- Local Variables ----!
       integer            :: i

       !> Initialize
       Err_ILLData=.false.
       Err_ILLData_Mess=' '

       if (len_trim(filedir) == 0) then
          ILL_data_directory =" "      !data are in the current directory
          got_ILL_data_directory=.true.
          return
       end if

       !> Add separator if absent and ILL_data_directory is not the current directory
       i=len_trim(filedir)
       if (filedir(i:i) /= ops_sep) then
          ILL_data_directory=trim(filedir)//ops_sep
       else
          ILL_data_directory=trim(filedir)
       end if

       !> Check that the directory exist, otherwise rise an error condition
       if (.not. directory_exists(trim(ILL_data_directory))) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess="The ILL directory: '"//trim(ILL_data_directory)//"' doesn't exist"
          got_ILL_data_directory=.false.
       else
          got_ILL_data_directory=.true.
       end if

       return
    End Subroutine Set_ILL_Data_Directory

    !!----
    !!---- Subroutine Set_Instrm_directory(working_dir,instrm, iyear, icycle)
    !!----    character(len=*),  intent(in), optional :: instrm         ! name of the diffractometer
    !!----    character(len=*),  intent(in), optional :: working_dir    ! Data directory to search
    !!----    character(len=*),  intent(in), optional :: iyear          ! Year for search Numor
    !!----    character(len=*),  intent(in), optional :: icycle         ! Cycle for Search Numor
    !!----
    !!----    Assign the global public variable: Instrm_directory
    !!----    It is assumed that the subroutine Set_ILL_data_directory has already been called.
    !!----
    !!----    1� Look at the working_dir to define Instrm_directory
    !!----    2� Arguments instrm, iyear and icycle is used with ILL_Data_Directory
    !!----    3� If not year and cycle then use "data"
    !!----    If the directory doesn't exist the subroutine rises an error condition
    !!----    by putting ERR_ILLData=.true. and filling the error message
    !!----    variable ERR_ILLData_Mess.
    !!----
    !!---- Update: 08/03/2011
    !!
    Subroutine Set_Instrm_Directory(working_dir, instrm, iyear, icycle)
       !---- Argument ----!
       character(len=*),  intent(in), optional :: working_dir
       character(len=*),  intent(in), optional :: instrm  !Name of the instrument
       integer         ,  intent(in), optional :: iyear
       integer         ,  intent(in), optional :: icycle

       !---- Local Variables ----!
       integer          :: i
       character(len=5) :: yearcycle

       ! Initialize
       ERR_ILLData=.false.
       ERR_ILLData_Mess=""

       ! If a working directory is given as an argument, then use it directly as the instrument
       ! directory.
       if ( Present(working_dir) ) then
           Instrm_directory = trim(working_dir)
           i=len_trim(Instrm_directory)
           if ( i > 0) then
              if (Instrm_directory(i:i) /= ops_sep) Instrm_directory=trim(Instrm_directory)//ops_sep
           else
              ERR_ILLData=.true.
              ERR_ILLData_Mess="A working directory was empty"
              return
           end if
       else
           if ( Present(instrm) ) then
              ! If an instrument name is given as argument, then build the instrument directory
              ! using the path of the ILL internal database.
              ! If the year and the cycle are given as arguments, then construct the yearcycle directory.
              i=len_trim(instrm)
              if (i > 0) then
                 if ( Present(iyear) .and. Present(icycle) ) then
                    Write(Unit=Yearcycle, fmt='(i4.4,i1.1)') iyear, icycle
                    yearcycle = yearcycle(len_trim(yearcycle)-2:len_trim(yearcycle))
                    Instrm_directory = trim(ILL_data_directory)//trim(yearcycle)//ops_sep//trim(instrm)//ops_sep
                 else
                    ! Otherwise, use the current data directory as the base location for the instrument directory.
                    Instrm_directory = trim(ILL_data_directory)//'data'//ops_sep//trim(instrm)//ops_sep
                 end if
              else
                 ERR_ILLData=.true.
                 ERR_ILLData_Mess="An instrument name must be at least provided"
                 return
              end if
           else
              ERR_ILLData=.true.
              ERR_ILLData_Mess="A working directory or an instrument name must be at least provided."
              return
           end if
       end if

       !---- check that the directory exist, ----!
       !---- otherwise raise an error condition ----!
       Instrm_directory_set=directory_exists(trim(Instrm_directory))
       if (.not. Instrm_directory_set) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess="The INSTRM directory: '"//trim(Instrm_directory)//"' doesn't exist"
          Instrm_directory = " "
       end if

       return
    End Subroutine Set_Instrm_Directory
    !!----
    !!---- Subroutine Set_Instrm_Geometry_Directory(env_var)
    !!----    character(len=*),  intent(in), optional :: env_var
    !!----
    !!----    Assign the global public variable: Instrm_Geometry_Directory
    !!----    The optional intent 'in' variable env_var is the name of any
    !!----    directory set as the instrument geometry directory.
    !!----    If the directory doesn't exist the subroutine rises an error condition
    !!----    by putting ERR_ILLData=.true. and filling the error message
    !!----    variable ERR_ILLData_Mess. If the subroutine is invoked without argument
    !!----    the Instrm_Geometry_Directory="" corresponds to the current directory.
    !!----
    !!---- Update: July - 2010
    !!
    Subroutine Set_Instrm_Geometry_Directory(env_var)
       !---- Argument ----!
       character(len=*),  intent(in), optional :: env_var

       !---- Local Variables ----!
       integer:: i

       ! If an environment variable is given as an argument, then use it directly as the instrument
       ! geometry directory, otherwise take the current directory as the place where geometry instrument
       ! files are stored.
       if ( Present(env_var) ) then
           Instrm_Geometry_directory = trim(env_var)
           i=len_trim(Instrm_Geometry_Directory)
           if (Instrm_Geometry_Directory(i:i) /= ops_sep) Instrm_Geometry_Directory=trim(Instrm_Geometry_Directory)//ops_sep
           !---- check that the directory exist, ----!
           !---- otherwise raise an error condition ----!
           ERR_ILLData=.false.
           ERR_ILLData_Mess=" "
           Instrm_Geometry_directory_set=directory_exists(trim(Instrm_Geometry_Directory))
           if (.not. Instrm_Geometry_directory_set) then
              ERR_ILLData=.true.
              ERR_ILLData_Mess=&
              "The INSTRM Geometry directory: '"//trim(Instrm_Geometry_Directory)//"' doesn't exist, currect directory assumed"
              Instrm_Geometry_Directory = " "
           end if
       else
           Instrm_Geometry_Directory = " "       !Current directory for searching .geom files
           Instrm_Geometry_directory_set=.true.
       end if

       return
    End Subroutine Set_Instrm_Geometry_Directory

    !!--++
    !!--++ Subroutine Set_KeyTypes_on_File(filevar, nlines)
    !!--++    character(len=*),dimension(:), intent(in) :: filevar
    !!--++    integer,                       intent(in) :: nlines
    !!--++
    !---++ Set the information on NL_KEYTYPES
    !!--++
    !!--++ Update: April-2009
    !!
    Subroutine Set_KeyTypes_on_File(filevar, nlines)
       !---- Arguments ----!
       character(len=*),dimension(:), intent(in) :: filevar
       integer,                       intent(in) :: nlines

       !---- Local Variables ----!
       integer :: i,j,nl,ndim

       if (allocated(nl_keytypes)) deallocate(nl_keytypes)
       if (all(n_keytypes == 0)) return

       ndim=maxval(n_keytypes)
       allocate(nl_keytypes(7,ndim,2))
       nl_keytypes=0

       do i=1,7
          if (n_keytypes(i) == 0) cycle
          j=1
          do nl=1,nlines
             select case (i)
                case (1)
                   if (nl_keytypes(1,j,1) == 0) then
                      if (filevar(nl)(1:10) /= 'RRRRRRRRRR') cycle
                      nl_keytypes(1,j,1)=nl
                   else
                      if (filevar(nl)(1:10) == 'RRRRRRRRRR' .or. &
                          filevar(nl)(1:10) == 'AAAAAAAAAA' .or. &
                          filevar(nl)(1:10) == 'SSSSSSSSSS' .or. &
                          filevar(nl)(1:10) == 'FFFFFFFFFF' .or. &
                          filevar(nl)(1:10) == 'IIIIIIIIII' .or. &
                          filevar(nl)(1:10) == 'JJJJJJJJJJ' .or. &
                          filevar(nl)(1:10) == 'VVVVVVVVVV') then
                         nl_keytypes(1,j,2)=nl-1
                         j=j+1
                         if (j > n_keytypes(i)) exit
                      end if
                      if (filevar(nl)(1:10) == 'RRRRRRRRRR') then
                         nl_keytypes(1,j,1)=nl
                      end if
                   end if

                case (2)
                   if (nl_keytypes(2,j,1) == 0) then
                      if (filevar(nl)(1:10) /= 'AAAAAAAAAA') cycle
                      nl_keytypes(2,j,1)=nl
                   else
                      if (filevar(nl)(1:10) == 'RRRRRRRRRR' .or. &
                          filevar(nl)(1:10) == 'AAAAAAAAAA' .or. &
                          filevar(nl)(1:10) == 'SSSSSSSSSS' .or. &
                          filevar(nl)(1:10) == 'FFFFFFFFFF' .or. &
                          filevar(nl)(1:10) == 'IIIIIIIIII' .or. &
                          filevar(nl)(1:10) == 'JJJJJJJJJJ' .or. &
                          filevar(nl)(1:10) == 'VVVVVVVVVV') then
                         nl_keytypes(2,j,2)=nl-1
                         j=j+1
                         if (j > n_keytypes(i)) exit
                      end if
                      if (filevar(nl)(1:10) == 'AAAAAAAAAA') then
                         nl_keytypes(2,j,1)=nl
                      end if
                   end if

                case (3)
                   if (nl_keytypes(3,j,1) == 0) then
                      if (filevar(nl)(1:10) /= 'SSSSSSSSSS') cycle
                      nl_keytypes(3,j,1)=nl
                   else
                      if (filevar(nl)(1:10) == 'RRRRRRRRRR' .or. &
                          filevar(nl)(1:10) == 'AAAAAAAAAA' .or. &
                          filevar(nl)(1:10) == 'SSSSSSSSSS' .or. &
                          filevar(nl)(1:10) == 'FFFFFFFFFF' .or. &
                          filevar(nl)(1:10) == 'IIIIIIIIII' .or. &
                          filevar(nl)(1:10) == 'JJJJJJJJJJ' .or. &
                          filevar(nl)(1:10) == 'VVVVVVVVVV') then
                         nl_keytypes(3,j,2)=nl-1
                         j=j+1
                         if (j > n_keytypes(i)) exit
                      end if
                      if (filevar(nl)(1:10) == 'SSSSSSSSSS') then
                         nl_keytypes(3,j,1)=nl
                      end if
                   end if

                case (4)
                   if (nl_keytypes(4,j,1) == 0) then
                      if (filevar(nl)(1:10) /= 'FFFFFFFFFF') cycle
                      nl_keytypes(4,j,1)=nl
                   else
                      if (filevar(nl)(1:10) == 'RRRRRRRRRR' .or. &
                          filevar(nl)(1:10) == 'AAAAAAAAAA' .or. &
                          filevar(nl)(1:10) == 'SSSSSSSSSS' .or. &
                          filevar(nl)(1:10) == 'FFFFFFFFFF' .or. &
                          filevar(nl)(1:10) == 'IIIIIIIIII' .or. &
                          filevar(nl)(1:10) == 'JJJJJJJJJJ' .or. &
                          filevar(nl)(1:10) == 'VVVVVVVVVV') then
                         nl_keytypes(4,j,2)=nl-1
                         j=j+1
                         if (j > n_keytypes(i)) exit
                      end if
                      if (filevar(nl)(1:10) == 'FFFFFFFFFF') then
                         nl_keytypes(4,j,1)=nl
                      end if
                   end if

                case (5)
                   if (nl_keytypes(5,j,1) == 0) then
                      if (filevar(nl)(1:10) /= 'IIIIIIIIII') cycle
                      nl_keytypes(5,j,1)=nl
                   else
                      if (filevar(nl)(1:10) == 'RRRRRRRRRR' .or. &
                          filevar(nl)(1:10) == 'AAAAAAAAAA' .or. &
                          filevar(nl)(1:10) == 'SSSSSSSSSS' .or. &
                          filevar(nl)(1:10) == 'FFFFFFFFFF' .or. &
                          filevar(nl)(1:10) == 'IIIIIIIIII' .or. &
                          filevar(nl)(1:10) == 'JJJJJJJJJJ' .or. &
                          filevar(nl)(1:10) == 'VVVVVVVVVV') then
                         nl_keytypes(5,j,2)=nl-1
                         j=j+1
                         if (j > n_keytypes(i)) exit
                      end if
                      if (filevar(nl)(1:10) == 'IIIIIIIIII') then
                         nl_keytypes(5,j,1)=nl
                      end if
                   end if

                case (6)
                   if (nl_keytypes(6,j,1) == 0) then
                      if (filevar(nl)(1:10) /= 'JJJJJJJJJJ') cycle
                      nl_keytypes(6,j,1)=nl
                   else
                      if (filevar(nl)(1:10) == 'RRRRRRRRRR' .or. &
                          filevar(nl)(1:10) == 'AAAAAAAAAA' .or. &
                          filevar(nl)(1:10) == 'SSSSSSSSSS' .or. &
                          filevar(nl)(1:10) == 'FFFFFFFFFF' .or. &
                          filevar(nl)(1:10) == 'IIIIIIIIII' .or. &
                          filevar(nl)(1:10) == 'JJJJJJJJJJ' .or. &
                          filevar(nl)(1:10) == 'VVVVVVVVVV') then
                         nl_keytypes(6,j,2)=nl-1
                         j=j+1
                         if (j > n_keytypes(i)) exit
                      end if
                      if (filevar(nl)(1:10) == 'JJJJJJJJJJ') then
                         nl_keytypes(6,j,1)=nl
                      end if
                   end if

                case (7)
                   if (nl_keytypes(7,j,1) == 0) then
                      if (filevar(nl)(1:10) /= 'VVVVVVVVVV') cycle
                      nl_keytypes(7,j,1)=nl
                   else
                      if (filevar(nl)(1:10) == 'RRRRRRRRRR' .or. &
                          filevar(nl)(1:10) == 'AAAAAAAAAA' .or. &
                          filevar(nl)(1:10) == 'SSSSSSSSSS' .or. &
                          filevar(nl)(1:10) == 'FFFFFFFFFF' .or. &
                          filevar(nl)(1:10) == 'IIIIIIIIII' .or. &
                          filevar(nl)(1:10) == 'JJJJJJJJJJ' .or. &
                          filevar(nl)(1:10) == 'VVVVVVVVVV') then
                         nl_keytypes(7,j,2)=nl-1
                         j=j+1
                         if (j > n_keytypes(i)) exit
                      end if
                      if (filevar(nl)(1:10) == 'VVVVVVVVVV') then
                         nl_keytypes(7,j,1)=nl
                      end if
                   end if
             end select
          end do
       end do

       do i=1,7
          do j=1,n_keytypes(i)
             if (nl_keytypes(i,j,1) /= 0 .and. nl_keytypes(i,j,2)==0) then
                nl_keytypes(i,j,2)=nlines
                exit
             end if
          end do
       end do

       return
    End Subroutine Set_KeyTypes_on_File

    !!----
    !!---- Subroutine Update_Current_Instrm_UB(filenam,UB,wave)
    !!----    character(len=*),              intent(in) :: filenam
    !!----    real(kind=cp), dimension(3,3), intent(in) :: UB
    !!----    real(kind=cp),                 intent(in) :: wave
    !!----
    !!----    Subroutine updating the file 'filenam' where the characteristics
    !!----    of the current instrument are written. The global Current_Instrm
    !!----    variable is re-filled with new values of wavelength and UB-matrix.
    !!----    The file 'filenam' is re-written and the old version is saved with
    !!----    appended extension '.bak'.
    !!----    The Current_Orient global variable is also updated.
    !!----    In case of error the subroutine puts ERR_ILLData=.true.
    !!----    and fils the error message variable ERR_ILLData_Mess.
    !!----
    !!---- Update: December - 2005
    !!
    Subroutine Update_Current_Instrm_UB(filenam,UB,wave)
       !---- Arguments ----!
       character(len=*),              intent(in) :: filenam
       real(kind=cp), dimension(3,3), intent(in) :: UB
       real(kind=cp),                 intent(in) :: wave

       !---- Local variables ----!
       character(len=120), dimension(:), allocatable :: file_lines
       character(len=120)            :: line
       character(len=10)             :: key
       integer                       :: i, j, lun_out, lun, ier, nlines,jw,iw,iub
       real(kind=cp), dimension(3,3) :: set
       logical                       :: read_wave, read_UB

       if(.not. Current_Instrm_set) then
         ERR_ILLData=.true.
         ERR_ILLData_Mess=" Current Instrument not set! (call subroutine: Read_Current_Instrm) "
         return
       end if

       call Number_Lines(filenam,nlines)
       if(nlines == 0) then
         ERR_ILLData=.true.
         ERR_ILLData_Mess="Error opening the file: "//trim(filenam)//" => ZERO lines found!"
         return
       end if

       if(allocated(file_lines)) deallocate(file_lines)
       allocate(file_lines(nlines))

       call Get_LogUnit(lun)
       open(unit=lun,file=trim(filenam),status="old", action="read", position="rewind",iostat=ier)
       if(ier /= 0) then
         ERR_ILLData=.true.
         ERR_ILLData_Mess="Error opening the file: "//trim(filenam)
         return
       end if

       call Get_LogUnit(lun_out)
       open(unit=lun_out,file=trim(filenam)//".bak",status="replace", action="write",iostat=ier)
       if(ier /= 0) then
         ERR_ILLData=.true.
         ERR_ILLData_Mess="Error opening the backup file: "//trim(filenam)//".bak"
         return
       end if

       do i=1,nlines
          read(unit=lun,    fmt="(a)") line
         write(unit=lun_out,fmt="(a)") line
         file_lines(i)=line
       end do
       close(unit=lun)
       close(unit=lun_out)

       read_wave=.false.
       read_UB=.false.
       jw=0

       do i=1,nlines
          line=adjustl(file_lines(i))
          if(line(1:1) == "!" .or. line(1:1) == "#" .or. len_trim(line) == 0) cycle
          j=index(line," ")
          key=u_case(line(1:j-1))

          Select Case(key)
             Case("WAVE")
                read_wave=.true.
                jw = j
                iw = i
             Case("UBMAT")
                read_UB=.true.
                iub=i
          End Select

          if( read_wave .and. read_UB) exit
       End do

       !Update the corresponding lines
       write(unit=file_lines(iw)(jw+1:), fmt="(f10.5)") wave
       do i=1,3
          write(unit=file_lines(iub+i),   fmt="(3f12.7)") ub(i,:)
       end do

       set(:,1)=Current_Instrm%e1
       set(:,2)=Current_Instrm%e2
       set(:,3)=Current_Instrm%e3
       call Set_Current_Orient(wave,ub,set)

       call Get_LogUnit(lun)
       open(unit=lun,file=trim(filenam),status="replace", action="write",iostat=ier)
       if(ier /= 0) then
         ERR_ILLData=.true.
         ERR_ILLData_Mess="Error updating the file: "//trim(filenam)
         return
       end if

       do i=1,nlines
          write(unit=lun,fmt="(a)",iostat=ier) file_lines(i)
          if(ier /= 0) then
            ERR_ILLData=.true.
            write(unit=ERR_ILLData_Mess,fmt="(a,i4)")"Error updating the file: "//trim(filenam)//" at line:",i
            exit
          end if
       end do
       close(unit=lun)

       return
    End Subroutine Update_Current_Instrm_UB

    !!----
    !!---- Subroutine Write_Current_Instrm_data(lun)
    !!----    integer, optional, intent(in) :: lun
    !!----
    !!----    Writes the characteristics of the Current Instrument
    !!----    in the file of logical unit 'lun'
    !!----    If the subroutine is invoked without argument the subroutine
    !!----    outputs the information on the standard output (screen)
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Write_Current_Instrm_data(lun,fil)
       !---- Arguments ----!
       integer,         optional, intent(in) :: lun
       character(len=*),optional, intent(in) :: fil

       !--- Local variables ---!
       integer           :: ipr,i
       character(len=16) :: forma


       ipr=6
       if (present(lun)) ipr=lun
       if (.not. present(fil)) then
          write(unit=ipr,fmt="(a)") " ---------------------------------------------------------------------------"
          write(unit=ipr,fmt="(a)") " ---   INFORMATION ABOUT THE CURRENT INSTRUMENT AND ORIENTATION MATRIX   ---"
          write(unit=ipr,fmt="(a)") " ---------------------------------------------------------------------------"
          write(unit=ipr,fmt="(a)") " "

          write(unit=ipr,fmt="(a)")          "        INFO: "//Current_Instrm%info
          write(unit=ipr,fmt="(a)")          "        NAME: "//Current_Instrm%name_inst
          write(unit=ipr,fmt="(a)")          "        GEOM: "//Current_Instrm%geom
          write(unit=ipr,fmt="(a)")          "    BL_FRAME: "//Current_Instrm%BL_frame
          write(unit=ipr,fmt="(a)")          "  DIST_UNITS: "//Current_Instrm%dist_units
          write(unit=ipr,fmt="(a)")          "  ANGL_UNITS: "//Current_Instrm%angl_units
          write(unit=ipr,fmt="(a)")          "    DET_TYPE: "//Current_Instrm%detector_type
          write(unit=ipr,fmt="(a,f8.3,a)")   "  DIST_SDETR: ", Current_Instrm%dist_samp_detector," "//Current_Instrm%dist_units
          if(index(Current_Instrm%geom,"Laue") /= 0) then
            write(unit=ipr,fmt="(a,f8.3,a)") "  LAMBDA_MIN: ", Current_Instrm%wave_min," Angstroms"
            write(unit=ipr,fmt="(a,f8.3,a)") "  LAMBDA_MAX: ", Current_Instrm%wave_max," Angstroms"
          end if
          write(unit=ipr,fmt="(a,f8.3,a)")   "  HORIZ_SIZE: ", Current_Instrm%horiz," "//Current_Instrm%dist_units
          write(unit=ipr,fmt="(a,f8.3,a)")   "   VERT_SIZE: ", Current_Instrm%vert ," "//Current_Instrm%dist_units
          write(unit=ipr,fmt="(a,f8.3,a)")   "   ANODE_GAP: ", Current_Instrm%agap," "//Current_Instrm%dist_units
          write(unit=ipr,fmt="(a,f8.3,a)")   " CATHODE_GAP: ", Current_Instrm%cgap ," "//Current_Instrm%dist_units
          write(unit=ipr,fmt="(a,i5    )")   "  HORIZ_PXEL: ", Current_Instrm%np_horiz
          write(unit=ipr,fmt="(a,i5    )")   "   VERT_PXEL: ", Current_Instrm%np_vert

          if (Current_Instrm%num_ang > 0) then
             write(unit=ipr,fmt="(a,i5    )")   "  NUMBER_ANG: ", Current_Instrm%num_ang
             write(unit=ipr,fmt="(a    )")      "  Names of Angular Motors, Limits (Min-Max) and Offsets:   "
             do i=1,Current_Instrm%num_ang
                write(unit=ipr,fmt="(a,3f9.3)")  "         "//Current_Instrm%ang_names(i), &
                Current_Instrm%ang_limits(i,1:2), Current_Instrm%ang_offsets(i)
             end do
          end if

          if (Current_Instrm%num_disp > 0) then
             write(unit=ipr,fmt="(a,i5    )")   "  NUMBER_DISP: ", Current_Instrm%num_disp
             write(unit=ipr,fmt="(a    )")      "  Names of Displacement Motors, Limits (Min-Max) and Offsets:   "
             do i=1,Current_Instrm%num_disp
                write(unit=ipr,fmt="(a,3f9.3)")   "         "//Current_Instrm%disp_names(i), &
                Current_Instrm%ang_limits(i,1:2), Current_Instrm%disp_offsets(i)
             end do
          end if

          write(unit=ipr,fmt="(a,3f8.2 )")   "  DET_OFFSET: ", Current_Instrm%det_offsets(1:3)
          write(unit=ipr,fmt="(a,3f8.3 )")   "   e1_COMPNT: ", Current_Instrm%e1
          write(unit=ipr,fmt="(a,3f8.3 )")   "   e2_COMPNT: ", Current_Instrm%e2
          write(unit=ipr,fmt="(a,3f8.3 )")   "   e3_COMPNT: ", Current_Instrm%e3

          write(unit=ipr,fmt="(a)") " "
          write(unit=ipr,fmt="(a)") "       ORIENTATION MATRIX  UB                       INVERSE OF UB-MATRIX"
          write(unit=ipr,fmt="(a)") " ---------------------------------------------------------------------------"
          write(unit=ipr,fmt="(a)") " "
          write(unit=ipr,fmt="(tr2,3f11.6,tr8,3f11.5)")  Current_Orient%ub(1,:), Current_Orient%ubinv(1,:)
          write(unit=ipr,fmt="(tr2,3f11.6,tr8,3f11.5)")  Current_Orient%ub(2,:), Current_Orient%ubinv(2,:)
          write(unit=ipr,fmt="(tr2,3f11.6,tr8,3f11.5)")  Current_Orient%ub(3,:), Current_Orient%ubinv(3,:)
          write(unit=ipr,fmt="(a)") " "
          if(index(Current_Instrm%geom,"Laue") == 0) write(unit=ipr,fmt="(a,f8.4,a)")  &
                                         "  WAVELENGTH:",Current_Orient%wave," angstroms"

       else  !Write a *.geom file

          write(unit=ipr,fmt="(a)") "INFO "//Current_Instrm%info
          write(unit=ipr,fmt="(a)") "NAME "//Current_Instrm%name_inst
          write(unit=ipr,fmt="(a,i3,a)") "GEOM ",Current_Instrm%igeom,"  "//Current_Instrm%geom
          write(unit=ipr,fmt="(a)") "BLFR "//Current_Instrm%BL_frame
          write(unit=ipr,fmt="(a)") "DIST_UNITS "//Current_Instrm%dist_units
          write(unit=ipr,fmt="(a)") "ANGL_UNITS "//Current_Instrm%angl_units
          write(unit=ipr,fmt="(a,i3)") "DET_TYPE "//trim(Current_Instrm%detector_type)//" ipsd ", &
                                        Current_Instrm%ipsd
          write(unit=ipr,fmt="(a,f8.3)") "DIST_DET ",Current_Instrm%dist_samp_detector
          if(index(Current_Instrm%geom,"Laue") /= 0) then
            write(unit=ipr,fmt="(a,f8.3,a)") "LAMBDA_MIN: ", Current_Instrm%wave_min
            write(unit=ipr,fmt="(a,f8.3,a)") "LAMBDA_MAX: ", Current_Instrm%wave_max
          end if
          write(unit=ipr,fmt="(a,2f8.3,2i5)") "DIM_XY ", Current_Instrm%horiz,    Current_Instrm%vert, &
                                                    Current_Instrm%np_horiz, Current_Instrm%np_vert
          write(unit=ipr,fmt="(a,2f8.4)") "GAPS_DET ",Current_Instrm%agap, Current_Instrm%cgap
          write(unit=ipr,fmt="(a,f8.5)") "WAVE ",Current_Orient%wave
          if(index(Current_Instrm%geom,"Laue") == 0) then
            write(unit=ipr,fmt="(a)") "UBMAT"
            do i=1,3
               write(unit=ipr,fmt="(3f12.7)") Current_Orient%ub(i,:)
            end do
          End If
          write(unit=ipr,fmt="(a,9f6.1)") "SETTING ",Current_Instrm%e1, Current_Instrm%e2, Current_Instrm%e3

          write(unit=ipr,fmt="(a,i4)") "NUM_ANG  ",Current_Instrm%num_ang
          write(unit=ipr,fmt="(a,i4)") "NUM_DISP ",Current_Instrm%num_disp
          write(unit=ipr,fmt="(a)")    "ANG_LIMITS"
          do i=1,Current_Instrm%num_ang
             write(unit=ipr,fmt="(a,2f8.2,f10.4)") "      "//Current_Instrm%ang_names(i), &
                                                             Current_Instrm%ang_Limits(i,1:2), &
                                                             Current_Instrm%ang_offsets(i)
          end do
          if (Current_Instrm%num_disp > 0) then
             write(unit=ipr,fmt="(a)")    "DISP_LIMITS"
             do i=1,Current_Instrm%num_disp
                write(unit=ipr,fmt="(a,2f8.2,f10.4)") "      "//Current_Instrm%disp_names(i), &
                                                                Current_Instrm%disp_Limits(i,1:2), &
                                                                Current_Instrm%disp_offsets(i)
             end do
          end if
          write(unit=ipr,fmt="(a,3f10.4)") "DET_OFF ", Current_Instrm%det_offsets
          if(index(Current_Instrm%geom,"Laue") == 0 .and. Current_Instrm%np_vert < 513 .and. &
                   Current_Instrm%np_horiz < 513) then
            write(unit=ipr,fmt="(a)") "DET_ALPHAS"
            forma="(     f8.4)"
            write(unit=forma(2:6),fmt="(i5)") Current_Instrm%np_horiz
            do i=1,Current_Instrm%np_vert
               write(unit=ipr,fmt=forma) Current_Instrm%alphas(i,1:Current_Instrm%np_horiz)
            end do
          end if
       end if

       return
    End Subroutine Write_Current_Instrm_data

    !!----
    !!---- Subroutine Write_Generic_Numor(Num,lun)
    !!----    type(generic_numor_type), intent(in) :: N
    !!----    integer, optional,        intent(in) :: lun
    !!----
    !!---- Write the information of a Generic Numor
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Write_Generic_Numor(Num,lun)
       !---- Arguments ----!
       type(generic_numor_type), intent(in) :: Num
       integer, optional,        intent(in) :: lun

       !---- Local Variables ----!
       integer :: ilun,i

       ilun=6
       if (present(lun)) ilun=lun

       write(unit=ilun, fmt='(a)') '#### Numor Information ####'
       write(unit=ilun, fmt='(a,i6.6)') 'Numor: ',num%numor

       write(unit=ilun,fmt='(a,t50,a)') 'Instrument: '//trim(num%instr),'Date: '//trim(num%date)
       write(unit=ilun,fmt='(a)') 'Experimental Name: '//trim(num%expname)
       write(unit=ilun,fmt='(a)') ' '

       write(unit=ilun,fmt='(a)') '#---- Sample Information ----#'
       do i=1,num%sampleid%n
          write(unit=ilun,fmt='(a)') trim(num%sampleid%namevar(i))//': '//trim(num%sampleid%cvalues(i))
       end do
       write(unit=ilun,fmt='(a)') ' '

       select case (num%instr)
           case ('D20')
              write(unit=ilun,fmt='(a)') '#---- Diffractometer Optics and Reactor Parameters ----#'
              do i=1,num%diffopt%n
                 if (len_trim(num%diffopt%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(num%diffopt%namevar(i))//': ',num%diffopt%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Monochromator Motor Parameters ----#'
              do i=1,num%monmpar%n
                 if (len_trim(num%monmpar%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(num%monmpar%namevar(i))//': ',num%monmpar%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Diffractometer Motor Parameters ----#'
              do i=1,num%diffmpar%n
                 if (len_trim(num%diffmpar%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(num%diffmpar%namevar(i))//': ',num%diffmpar%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Detector and DAS Parameters ----#'
              do i=1,num%detpar%n
                 if (len_trim(num%detpar%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(num%detpar%namevar(i))//': ',num%detpar%rvalues(i)
              end do
             write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Data Acquisition Parameters ----#'
              do i=1,num%dacparam%n
                 if (len_trim(num%dacparam%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(num%dacparam%namevar(i))//': ',num%dacparam%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Sample Status ----#'
              do i=1,num%samplest%n
                 if (len_trim(num%samplest%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(num%samplest%namevar(i))//': ',num%samplest%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Counts Information ----#'
              write(unit=ilun,fmt='(a,i8)') trim(num%icounts%namevar(1))//': ',num%icounts%n
              write(unit=ilun,fmt='(a)') ' '

           case ('D1B')
              write(unit=ilun,fmt='(a)') '#---- Data Acquisition Flags ----#'
              do i=1,num%dacflags%n
                 if (len_trim(num%dacflags%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,i6)') trim(num%dacflags%namevar(i))//': ',num%dacflags%ivalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Data Acquisition Parameters ----#'
              do i=1,num%dacparam%n
                 if (len_trim(num%dacparam%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.4)') trim(num%dacparam%namevar(i))//': ',num%dacparam%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Counts Information ----#'
              write(unit=ilun,fmt='(a,g15.4)') trim(num%icounts%namevar(1))//': ',real(num%icounts%ivalues(1))
              write(unit=ilun,fmt='(a,g15.4)') trim(num%icounts%namevar(2))//': ',real(num%icounts%ivalues(2))/60000.0
              write(unit=ilun,fmt='(a,g15.4)') trim(num%icounts%namevar(3))//': ',real(num%icounts%ivalues(3))*0.001
              write(unit=ilun,fmt='(a)') ' '

       end select

       return
    End Subroutine Write_Generic_Numor

    !!----
    !!---- Subroutine Write_HeaderInfo_POWDER_Numor(Num,lun)
    !!----    type(POWDER_Numor_type), intent(in) :: Num
    !!----    integer,       optional, intent(in) :: lun
    !!----
    !!----    Writes the characteristics of the numor objet 'Num'
    !!----    of type POWDER_Numor_type in the file of logical unit 'lun'.
    !!----    If the subroutine is invoked without the 'lun' argument the subroutine
    !!----    outputs the information on the standard output (screen)
    !!----
    !!---- Update: 15/03/2011
    !!
    Subroutine Write_HeaderInfo_POWDER_Numor(Num,lun)
       !---- Arguments ----!
       type(POWDER_Numor_type), intent(in) :: Num
       integer,       optional, intent(in) :: lun

       !--- Local variables ---!
       integer       :: ipr

       ipr=6
       if (present(lun)) ipr=lun

       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a,i6,a)") " ---  Header Information about the NUMOR: ",Num%numor," of instrument "//trim(num%instrm)
       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)")        "         HEADER: "//trim(Num%header)
       write(unit=ipr,fmt="(a)")        "          TITLE: "//trim(Num%title )
       write(unit=ipr,fmt="(a)")        "      INTRUMENT: "//trim(Num%Instrm)
       write(unit=ipr,fmt="(a)")        "      SCAN_TYPE: "//trim(Num%Scantype)
       write(unit=ipr,fmt="(a,f14.3)") "     SCAN_START: ", Num%scans(1)
       write(unit=ipr,fmt="(a,f14.3)") "     SCAN_STEP : ", Num%scans(2)
       write(unit=ipr,fmt="(a,f14.3 )")  "      MONITOR: ", Num%monitor
       write(unit=ipr,fmt="(a,f14.3 )")  "    TIME(seg): ", Num%time
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-SETPNT: ", Num%conditions(1)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-REGUL.: ", Num%conditions(2)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-SAMPLE: ", Num%conditions(3)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3  )") " Magnetic Field: ", Num%conditions(5)
       write(unit=ipr,fmt="(a)") " "

       return
    End Subroutine Write_HeaderInfo_POWDER_Numor

    !!----
    !!---- Subroutine Write_POWDER_Numor(Num,lun)
    !!----    type(POWDER_Numor_type), intent(in) :: Num
    !!----    integer,      optional, intent(in) :: lun
    !!----
    !!----    Writes the characteristics of the numor objet 'Num'
    !!----    of type POWDER_Numor_type in the file of logical unit 'lun'.
    !!----    If the subroutine is invoked without the 'lun' argument the subroutine
    !!----    outputs the information on the standard output (screen)
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Write_POWDER_Numor(Num,lun)
       !---- Arguments ----!
       type(POWDER_Numor_type), intent(in) :: Num
       integer,       optional, intent(in) :: lun

       !--- Local variables ---!
       integer       :: ipr, i, j, nlines,n
       integer, dimension(10) :: cou

       ipr=6
       if (present(lun)) ipr=lun

       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a,i6,a)") " ---  Information about the NUMOR: ",Num%numor," of instrument "//trim(num%instrm)
       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)")        "         HEADER: "//Num%header
       write(unit=ipr,fmt="(a)")        "          TITLE: "//Num%title
       write(unit=ipr,fmt="(a)")        "      INTRUMENT: "//Num%Instrm
       write(unit=ipr,fmt="(a)")        "      SCAN_TYPE: "//Num%Scantype
       write(unit=ipr,fmt="(a,f14.3)") "      ANGLE_PHI: ", Num%angles(1)
       write(unit=ipr,fmt="(a,f14.3)") "      ANGLE_CHI: ", Num%angles(2)
       write(unit=ipr,fmt="(a,f14.3)") "    ANGLE_OMEGA: ", Num%angles(3)
       write(unit=ipr,fmt="(a,f14.3)") "    ANGLE_GAMMA: ", Num%angles(4)
       write(unit=ipr,fmt="(a,f14.3)") "      ANGLE_PSI: ", Num%angles(5)
       write(unit=ipr,fmt="(a,f14.3)") "     SCAN_START: ", Num%scans(1)
       write(unit=ipr,fmt="(a,f14.3)") "     SCAN_STEP : ", Num%scans(2)
       write(unit=ipr,fmt="(a,f14.3)") "     SCAN_WIDTH: ", Num%scans(3)
       write(unit=ipr,fmt="(a,f14.3 )")  "      MONITOR: ", Num%monitor
       write(unit=ipr,fmt="(a,f14.3 )")  "    TIME(min): ", Num%time
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-SETPNT: ", Num%conditions(1)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-REGUL.: ", Num%conditions(2)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-SAMPLE: ", Num%conditions(3)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3  )") "      VOLTMETER: ", Num%conditions(4)
       write(unit=ipr,fmt="(a,f14.3  )") " Magnetic Field: ", Num%conditions(5)
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a,i6)")     "  Calc. angl. type: ",Num%icalc
       write(unit=ipr,fmt="(a,i6)")     "  Principal  angle: ",Num%manip
       write(unit=ipr,fmt="(a,i6)")     "  Number of frames: ",Num%nframes
       write(unit=ipr,fmt="(a,i6)")     "  Number of angles: ",Num%nbang
       write(unit=ipr,fmt="(a,i6)")     "  Number of pixels: ",Num%nbdata
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)")        "  Profile with all counts in detector: "


       do i=1, Num%nframes
          nlines=Num%nbdata/10+1
          n=1
          do j=1,nlines
             cou=nint(Num%counts(n:n+9,i))
             write(unit=ipr,fmt="(10i8)") cou
             n=n+10
             if (n > Num%nbdata) exit
          end do
       end do

       return
    End Subroutine Write_POWDER_Numor

    !!----
    !!---- Subroutine Write_HeaderInfo_SXTAL_Numor(Num,lun)
    !!----    type(SXTAL_Numor_type), intent(in) :: Num
    !!----    integer,      optional, intent(in) :: lun
    !!----
    !!----    Writes the Header Information of the numor objet 'Num'
    !!----    of type SXTAL_Numor_type in the file of logical unit 'lun'.
    !!----    If the subroutine is invoked without the 'lun' argument the subroutine
    !!----    outputs the information on the standard output (screen)
    !!----
    !!---- Update: 15/03/2011
    !!
    Subroutine Write_HeaderInfo_SXTAL_Numor(Num,lun)
       !---- Arguments ----!
       type(SXTAL_Numor_type), intent(in) :: Num
       integer,      optional, intent(in) :: lun

       !--- Local variables ---!
       integer       :: ipr

       ipr=6
       if (present(lun)) ipr=lun

       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a,i6,a)") " ---  Information about the NUMOR: ",Num%numor," of instrument "//trim(Num%Instrm)
       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)")        "         HEADER: "//trim(Num%header)
       write(unit=ipr,fmt="(a)")        "          TITLE: "//trim(Num%title )
       write(unit=ipr,fmt="(a)")        "      INTRUMENT: "//trim(Num%Instrm)
       write(unit=ipr,fmt="(a)")        "      SCAN_TYPE: "//trim(Num%Scantype)
       write(unit=ipr,fmt="(a,3f8.3)")  "      HKL(Min.): ", Num%hmin
       write(unit=ipr,fmt="(a,3f8.3)")  "      HKL(Max.): ", Num%hmax
       if (Num%Scantype == "phi") then
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_START: ", Num%scans(1) ," r.l.u."
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_STEP : ", Num%scans(2) ," r.l.u."
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_WIDTH: ", Num%scans(3) ," r.l.u."
       else
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_START: ", Num%scans(1) ," "//Current_Instrm%angl_units
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_STEP : ", Num%scans(2) ," "//Current_Instrm%angl_units
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_WIDTH: ", Num%scans(3) ," "//Current_Instrm%angl_units
       end if
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-SETPNT: ", Num%conditions(1)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-REGUL.: ", Num%conditions(2)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-SAMPLE: ", Num%conditions(3)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3  )") "      VOLTMETER: ", Num%conditions(4)
       write(unit=ipr,fmt="(a,f14.3  )") " Magnetic Field: ", Num%conditions(5)
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)") "       ORIENTATION MATRIX  UB"
       write(unit=ipr,fmt="(a)") " --------------------------------------"
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(tr2,3f12.7)")  Num%ub(1,:)
       write(unit=ipr,fmt="(tr2,3f12.7)")  Num%ub(2,:)
       write(unit=ipr,fmt="(tr2,3f12.7)")  Num%ub(3,:)
       write(unit=ipr,fmt="(a)") " "

       return
    End Subroutine Write_HeaderInfo_SXTAL_Numor

    !!----
    !!---- Subroutine Write_SXTAL_Numor(Num,lun)
    !!----    type(SXTAL_Numor_type), intent(in) :: Num
    !!----    integer,      optional, intent(in) :: lun
    !!----
    !!----    Writes the characteristics of the numor objet 'Num'
    !!----    of type SXTAL_Numor_type in the file of logical unit 'lun'.
    !!----    If the subroutine is invoked without the 'lun' argument the subroutine
    !!----    outputs the information on the standard output (screen)
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Write_SXTAL_Numor(Num,lun)
       !---- Arguments ----!
       type(SXTAL_Numor_type), intent(in) :: Num
       integer,      optional, intent(in) :: lun

       !--- Local variables ---!
       integer       :: ipr, i, cou, ctot
       real(kind=cp) :: tim,mon,ang1,ang2

       ipr=6
       if (present(lun)) ipr=lun

       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a,i6,a)") " ---  Information about the NUMOR: ",Num%numor," of instrument "//trim(Num%Instrm)
       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)")        "         HEADER: "//Num%header
       write(unit=ipr,fmt="(a)")        "          TITLE: "//Num%title
       write(unit=ipr,fmt="(a)")        "      INTRUMENT: "//Num%Instrm
       write(unit=ipr,fmt="(a)")        "      SCAN_TYPE: "//Num%Scantype
       write(unit=ipr,fmt="(a,f8.3)")   "  COUPLING_FACT: ", Num%cpl_fact
       write(unit=ipr,fmt="(a,3f8.3)")  "      HKL(Min.): ", Num%hmin
       write(unit=ipr,fmt="(a,3f8.3)")  "      HKL(Max.): ", Num%hmax
       write(unit=ipr,fmt="(a,3f8.3)")  "      DELTA-hkl: ", Num%dh
       write(unit=ipr,fmt="(a,f14.3,a)") "      ANGLE_PHI: ", Num%angles(1)," "//Current_Instrm%angl_units
       write(unit=ipr,fmt="(a,f14.3,a)") "      ANGLE_CHI: ", Num%angles(2)," "//Current_Instrm%angl_units
       write(unit=ipr,fmt="(a,f14.3,a)") "    ANGLE_OMEGA: ", Num%angles(3)," "//Current_Instrm%angl_units
       write(unit=ipr,fmt="(a,f14.3,a)") "    ANGLE_GAMMA: ", Num%angles(4)," "//Current_Instrm%angl_units
       write(unit=ipr,fmt="(a,f14.3,a)") "      ANGLE_PSI: ", Num%angles(5)," "//Current_Instrm%angl_units
       if (Num%Scantype == "phi") then
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_START: ", Num%scans(1) ," r.l.u."
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_STEP : ", Num%scans(2) ," r.l.u."
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_WIDTH: ", Num%scans(3) ," r.l.u."
       else
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_START: ", Num%scans(1) ," "//Current_Instrm%angl_units
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_STEP : ", Num%scans(2) ," "//Current_Instrm%angl_units
          write(unit=ipr,fmt="(a,f14.3,a)") "     SCAN_WIDTH: ", Num%scans(3) ," "//Current_Instrm%angl_units
       end if
       write(unit=ipr,fmt="(a,f14.3 )") "         PRESET: ", Num%preset
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-SETPNT: ", Num%conditions(1)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-REGUL.: ", Num%conditions(2)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3,a)") "    Temp-SAMPLE: ", Num%conditions(3)," Kelvin"
       write(unit=ipr,fmt="(a,f14.3  )") "      VOLTMETER: ", Num%conditions(4)
       write(unit=ipr,fmt="(a,f14.3  )") " Magnetic Field: ", Num%conditions(5)
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)") "       ORIENTATION MATRIX  UB"
       write(unit=ipr,fmt="(a)") " --------------------------------------"
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(tr2,3f12.7)")  Num%ub(1,:)
       write(unit=ipr,fmt="(tr2,3f12.7)")  Num%ub(2,:)
       write(unit=ipr,fmt="(tr2,3f12.7)")  Num%ub(3,:)
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a,i6)")     "  Calc. angl. type: ",Num%icalc
       write(unit=ipr,fmt="(a,i6)")     "  Principal  angle: ",Num%manip
       write(unit=ipr,fmt="(a,i6)")     "  Number of frames: ",Num%nframes
       write(unit=ipr,fmt="(a,i6)")     "  Number of angles: ",Num%nbang
       write(unit=ipr,fmt="(a,i6)")     "  Number of pixels: ",Num%nbdata
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)")        "  Profile with all counts in detector: "

       do i=1, Num%nframes
          cou=nint(sum(Num%counts(:,i)))
          !ang= Num%scans(1)+ real(i-1)*Num%scans(2)
          tim = Num%tmc_ang(1,i)
          mon = Num%tmc_ang(2,i)
          ctot = nint(Num%tmc_ang(3,i))
          ang1 = Num%tmc_ang(4,i)
          ang2 = Num%tmc_ang(5,i)
          write(unit=ipr,fmt="(4f12.2,2i8)") tim,mon,ang1,ang2, ctot, cou
       end do

       return
    End Subroutine Write_SXTAL_Numor

    !!----
    !!---- Subroutine Read_Calibration_File_D1A(FileCal,Cal)
    !!----    character(len=*), intent(in)                 :: FileCal    ! Path+Filename of Calibration File
    !!----    type(calibration_detector_type), intent(out) :: Cal        ! Calibration Detector Object
    !!----
    !!---- Load the Calibration parameters for a D1A Instrument
    !!----
    !!---- 17/03/2011
    !!
    Subroutine Read_Calibration_File_D1A(FileCal,Cal)
        !---- Arguments ----!
        character(len=*), intent(in)                 :: FileCal
        type(calibration_detector_type), intent(out) :: Cal

        !---- Local Variables ----!
        character(len=512), dimension(:), allocatable :: filevar
        integer                                       :: iv,nlines
        integer, dimension(25)                        :: ivet
        real(kind=cp), dimension(25)                  :: vet

        ! Init
        err_illdata=.false.
        err_illdata_mess=' '

        set_calibration_detector=.false.

        ! Detecting numor
        call Number_Lines(FileCal,nlines)
        if (nlines <=0) then
           err_illdata=.true.
           err_illdata_mess=' Problems reading Calibration File for D1A Instrument'
           return
        end if

        ! Allocating variables
        if (allocated(filevar)) deallocate(filevar)
        allocate(filevar(nlines))
        call Reading_Lines(FileCal,nlines,filevar)

        ! Load Information on Cal Object
        Cal%Name_Instrm='D1A'
        Cal%NDet=25
        Cal%NPointsDet=1
        if (allocated(Cal%Active)) deallocate(Cal%Active)
        if (allocated(Cal%PosX))    deallocate(Cal%PosX)
        if (allocated(Cal%Effic))  deallocate(Cal%Effic)

        allocate(Cal%Active(Cal%NDet))
        allocate(Cal%PosX(Cal%NDet))
        allocate(Cal%Effic(Cal%NPointsDet,Cal%NDet))

        Cal%Active=.true.
        Cal%PosX=0.0
        Cal%Effic=1.0

        ! First line is not considered

        ! Second line is the relative position: Last to First
        call getnum(filevar(2), vet,ivet,iv)
        if (iv /= 25) then
           err_illdata=.true.
           err_illdata_mess=' Problems reading Positions of Detectors for D1A Instrument'
           deallocate(filevar)
           return
        end if
        Cal%PosX=vet

        ! Efficiency of each detector
        call getnum(filevar(2), vet,ivet,iv)
        if (iv /= 25) then
           err_illdata=.true.
           err_illdata_mess=' Problems reading Efficiency values for Detectors of D1A Instrument'
           deallocate(filevar)
           return
        end if
        Cal%Effic(:,1)=vet

        set_calibration_detector =.true.
        deallocate(filevar)

        return
    End Subroutine Read_Calibration_File_D1A


 End Module CFML_ILL_Instrm_Data