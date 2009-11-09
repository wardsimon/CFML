!!----
!!---- Copyleft(C) 1999-2009,              Version: 4.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: CFML_ILL_Instrm_Data
!!----   INFO: Subroutines related to Instrument information from ILL
!!----
!!---- HISTORY
!!----    Update: April - 2008
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
!!--++   Use CFML_GlobalDeps,        only: sp, dp, cp, pi, to_deg, to_rad, eps, OPS, OPS_sep
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
!!--++    ILL_TEMP_DIRECTORY                [Private]
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
!!----       ALLOCATE_POWDER_NUMORS
!!----       ALLOCATE_SXTAL_NUMORS
!!----       DEFINE_UNCOMPRESS_PROGRAM
!!----       GET_ABSOLUTE_DATA_PATH
!!----       GET_NEXT_YEARCYCLE
!!----       GET_SINGLE_FRAME_2D
!!----       INITIALIZE_DATA_DIRECTORY
!!--++       NUMBER_KEYTYPES_ON_FILE         [Private]
!!--++       READ_A_KEYTYPE                  [Private]
!!----       READ_CURRENT_INSTRM
!!--++       READ_F_KEYTYPE                  [Private]
!!--++       READ_I_KEYTYPE                  [Private]
!!--++       READ_J_KEYTYPE                  [Private]
!!----       READ_NUMOR_D1B
!!----       READ_NUMOR_D20
!!----       READ_POWDER_NUMOR
!!--++       READ_R_KEYTYPE                  [Private]
!!--++       READ_S_KEYTYPE                  [Private]
!!--++       READ_V_KEYTYPE                  [Private]
!!----       READ_SXTAL_NUMOR
!!----       SET_CURRENT_ORIENT
!!----       SET_DEFAULT_INSTRUMENT
!!----       SET_ILL_DATA_DIRECTORY
!!----       SET_INSTRM_DIRECTORY
!!----       SET_KEYTYPES_ON_FILE            [Private]
!!----       UPDATE_CURRENT_INSTRM_UB
!!----       WRITE_CURRENT_INSTRM_DATA
!!----       WRITE_GENERIC_NUMOR
!!----       WRITE_POWDER_NUMOR
!!----       WRITE_SXTAL_NUMOR
!!----
!!
Module CFML_ILL_Instrm_Data
   !---- Use Modules ----!
   !use f2kcli !Comment for compliant F2003 compilers
   Use CFML_GlobalDeps
   Use CFML_Math_General,     only: cosd,sind, equal_vector
   use CFML_String_Utilities, only: u_case, lcase, Get_LogUnit, Number_Lines, Reading_Lines
   use CFML_Math_3D,          only: err_math3d,err_math3d_mess, Cross_Product, Determ_A, Determ_V, &
                                    invert => Invert_A
   !---- Variables ----!
   Implicit none

   private

   !---- Public Subroutines ----!
   public :: Set_Current_Orient, Read_SXTAL_Numor, Read_Current_Instrm, Write_Current_Instrm_data,   &
             Allocate_SXTAL_numors, Write_SXTAL_Numor, Set_ILL_data_directory, Set_Instrm_directory, &
             Update_Current_Instrm_UB, Set_Default_Instrument,Get_Single_Frame_2D, &
             Initialize_Data_Directory, Get_Absolute_Data_Path, Get_Next_YearCycle, &
             Write_Generic_Numor, Read_Numor_D1B, Read_Numor_D20,                    &
             Allocate_Powder_Numors, Read_POWDER_Numor, Write_POWDER_Numor, Define_Uncompress_Program

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
   !!----    integer,dimension(7)        :: icdesc     !
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
      integer,dimension(7)        :: icdesc     !
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
   !!----    integer, dimension(7)                      :: icdesc      ! Integer values
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
      integer, dimension(7)                      :: icdesc      ! Integer values
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
   !!----    integer, dimension(7)                      :: icdesc      ! Integer values
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
      integer, dimension(7)                      :: icdesc      ! Integer values
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

   !!--++
   !!--++ ILL_TEMP_DIRECTORY
   !!--++    character(len=512), private :: ILL_temp_directory
   !!--++
   !!--++    (Private)
   !!--++    String containing information about the data directory for ILL
   !!--++    Initialised (automatic save attribute) for Windows case but set
   !!--++    In subroutine: Initialize_Data_Directory
   !!--++
   !!--++ Update: April - 2008
   !!
   character(len=512), private ::  ILL_temp_directory = "c:\diffraction_Windows\illdata\"

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
   character(len=512), public  :: Instrm_directory = " "

   !!--++
   !!--++ INSTRM_DIRECTORY_SET
   !!--++    logical, private:: Instrm_directory_set
   !!--++
   !!--++    logical variable taking the value .true. if set the instrument directory
   !!--++    in correct form
   !!--++
   !!--++ Update: April - 2008
   !!
   logical, private :: Instrm_directory_set=.false.

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
          allocate(Num(i)%tmc_ang(num_ang,nframes))

          if(allocated(Num(i)%counts)) deallocate(Num(i)%counts)
          allocate(Num(i)%counts(ndata,nframes))

          Num(i)%tmc_ang(:,:)=0.0
          Num(i)%counts(:,:)=0.0
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
          allocate(Num(i)%tmc_ang(num_ang,nframes))

          if(allocated(Num(i)%counts)) deallocate(Num(i)%counts)
          allocate(Num(i)%counts(ndata,nframes))

          Num(i)%tmc_ang(:,:)=0.0
          Num(i)%counts(:,:)=0.0
       end do

       return
    End Subroutine Allocate_SXTAL_Numors

    !!----
    !!---- Subroutine Define_Uncompress_Program(uncompress)
    !!----    character(len=*), intent(in) :: uncompress
    !!----
    !!---- Routine that define the uncompress program that you wants to use
    !!----
    !!---- Update:  April - 2009
    !!
    Subroutine Define_Uncompress_Program(uncompress)
       !---- Argument ----!
       character(len=*), intent(in) :: uncompress

       uncompresscommand=uncompress

       return
    End Subroutine Define_Uncompress_Program

    !!----
    !!---- Subroutine Get_Absolute_Data_Path(Numor,Instrm,Path,Iyear,Icycle)
    !!----    integer,           intent(in)  :: numor
    !!----    character(len=*),  intent(in)  :: instrm
    !!----    character(len=*),  intent(out) :: path
    !!----    integer, optional, intent(in)  :: iyear
    !!----    integer, optional, intent(in)  :: icycle
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
    Subroutine Get_Absolute_Data_Path(Numor, Instrm, Path, Iyear,Icycle)

       !---- Arguments ----!
       integer,           intent(in)  :: numor
       character(len=*),  intent(in)  :: instrm
       character(len=*),  intent(out) :: path
       integer, optional, intent(in)  :: iyear
       integer, optional, intent(in)  :: icycle

       !---- Local Variables ----!
       character(len=*),parameter :: Current_Data = 'data'
       character(len=*),parameter :: Previous_Data = 'data-1'
       character(len=*),parameter :: Extension = '.Z'

       character(len=6) :: numstr
       character(len=4) :: inst
       character(len=3) :: yearcycle
       character(len=7) :: subdir
       logical          :: exists

       integer          :: i

       ! Init value
       path=" "

       ! Some compilers are unable to test the existence of file: data_diretory//ops_sep//"." !!!!!
       ! so, got_ILL_data_directory may still be .false., ILL_data_directory has been set anyway
       ! in initialize_data_directory(), so let us continue and set got_ILL_data_directory=.true.
       ! if the associated numor file exist anywhere
       if (.not. got_ILL_data_directory) then
          call initialize_data_directory()
       else
          i=len_trim(ILL_data_directory)
          if (ILL_data_directory(i:i) /= ops_sep) then
             ILL_data_directory=trim(ILL_data_directory)//ops_sep
          end if
       end if

       ! At this point the ILL Data and Temporal directories must be defined
       if (.not. got_ILL_Data_directory) return
       if (.not. got_ILL_Temp_Directory) return

       ! Uncompress program must be defined
       if (len_trim(uncompresscommand) == 0) then
          select case (ops)
             case (1) ! Windows
                call define_uncompress_program('7z e -y -so')
             case (2:) ! Linux...
                call define_uncompress_program('gunzip -c')
          end select
       end if

       ! Numor character
       write(unit=numstr,fmt='(i6.6)') numor

       ! Instrument
       inst=adjustl(instrm)
       call lcase(inst)

       ! Using Instrument Information
       select case(inst)
          case("d1b","d20","d9","d15","d19")
             subdir = trim(inst)//"_"//numstr(2:2)//ops_sep
          case default ! d10, d3 etc don't divide numors up into sub directories
             subdir = ""
       end select

       if (present(iyear) .and. present(icycle)) then
          ! Using Year+Cycle
          write(unit=yearcycle,fmt="(i2.2,i1.1)") iyear,icycle
          path = trim(ILL_data_directory)//yearcycle//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
          inquire(file=trim(path),exist=exists)
          if (exists) return ! found numor so return

          ! Using previous + compressed data
          path = trim(path)//Extension
          inquire(file=trim(path),exist=exists)
          if (exists) then ! uncompress into temp directory
             call system(trim(uncompresscommand)//' '//trim(path)//' > '//trim(ILL_temp_directory)//numstr)
             path = trim(ILL_temp_directory)//numstr
             return ! found numor so return
          end if

          ! Using Year + Cycle of a previous call to routine
          write(unit=yearcycle,fmt="(i2.2,i1.1)") year_illdata,cycle_number
          path = trim(ILL_data_directory)//yearcycle//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
          inquire(file=trim(path),exist=exists)
          if (exists) return ! found numor so return

          ! Using previous + compressed data
          path = trim(path)//Extension
          inquire(file=trim(path),exist=exists)
          if (exists) then ! uncompress into temp directory
             call system(trim(uncompresscommand)//' '//trim(path)//' > '//trim(ILL_temp_directory)//numstr)
             path = trim(ILL_temp_directory)//numstr
             return ! found numor so return
          end if
       end if

       ! Using Current_data
       path = trim(ILL_data_directory)//Current_Data//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr

       inquire(file=trim(path),exist=exists)
       if (exists) return ! found numor so return

       ! Using Previous data
       path = trim(ILL_data_directory)//Previous_Data//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
       inquire(file=trim(path),exist=exists)
       if (exists) return ! found numor so return

       ! start from the most recent yearcycle and work search backwards
       call get_next_yearcycle(yearcycle,.true.)
       do
          if (yearcycle == "") exit

          path = trim(ILL_data_directory)//yearcycle//ops_sep//trim(inst)//ops_sep//trim(subdir)//numstr
          inquire(file=trim(path),exist=exists)
          if (exists) return

          path = trim(path)//Extension
          inquire(file=trim(path),exist=exists)
          if (exists) then ! uncompress into temp directory
             call system(trim(uncompresscommand)//' '//trim(path)//' > '//trim(ILL_temp_directory)//numstr)
             path = trim(ILL_temp_directory)//numstr
             return ! found numor so return
          end if
          call get_next_yearcycle(yearcycle)
       end do

       ! the numor was found compressed or uncompressed anywhere
       path = " "

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
          write(unit=yearcycle,fmt="(i2.2,i1.1)") year_illdata,cycle_number
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
    !!----    ILL_Temp_Directory: A temporary directory (used for uncompressing). If
    !!----                        found the flag got_ILL_temp_directory is set to true.
    !!----
    !!--..    Original code from Mike Turner, changed to make the subroutine independent
    !!--..    of the Winteracter library.
    !!--..
    !!--..    The subroutine ask for enviroment variables TEMP and ILLDATA in first place and
    !!--..    and if it return a blank then the default is set as following:
    !!--..    Windows: ILL_Temp_Directory -> C:\Temp     ILL_Data_Directory -> \\Serdon\illdata
    !!--..    Linux:   ILL_Temp_Directory -> $HOME/.tmp  ILL_Data_Directory -> /net/serdon/illdata
    !!----
    !!--..    NB: You need change this routine according to the compiler used !!!!!!
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Initialize_Data_Directory()
       !---- Local Variables ----!
       character(len=*), parameter :: Envvar1 = 'TEMP'
       character(len=*), parameter :: Envvar2 = 'ILLDATA'
       logical                     :: exists


       select case (OPS)
          case (1)  ! Windows
             ! Temporal
             call getenv(ENVVAR1, ILL_Temp_Directory)
             if (len_trim(ILL_Temp_directory) == 0) then
                ILL_Temp_directory='C:\Temp'
             end if

             ! ILL Data
             call getenv(ENVVAR2, ILL_Data_Directory)
             if (len_trim(ILL_Data_Directory) == 0) then
                ILL_Data_Directory = '\\Serdon\illdata'
             end if

          case (2:)  !Different variants of UNIX (Linux, MacOS, ...)
             ! Temporal
             call getenv(ENVVAR1, ILL_Temp_Directory)
             if (len_trim(ILL_Temp_directory) == 0) then
                call get_environment_variable('HOME', ILL_Temp_Directory)
                ILL_Temp_directory =trim(ILL_Temp_directory)//OPS_SEP//'tmp'
             end if

             ! ILL Data
             call get_environment_variable(ENVVAR2, ILL_Data_Directory)
             if (len_trim(ILL_Data_Directory) == 0) then
                ILL_Data_Directory = '/net/serdon/illdata'
             end if
       end select
       ! Temporal
       ILL_Temp_directory=trim(ILL_Temp_directory)//ops_sep
       exists=directory_exists(trim(ILL_temp_directory))
       if (exists) got_ILL_temp_directory = .true.

       ! ILL DATA
       ILL_Data_Directory=trim(ILL_Data_Directory)//ops_sep
       exists=directory_exists(trim(ILL_Data_Directory))
       if (exists) got_ILL_data_directory = .true.

       return
    End Subroutine Initialize_Data_Directory

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
          n=nchars/80 + mod(nchars,80)
          do i=1,n
             charline=trim(charline)//filevar(n_ini+1+i)//char(9)
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

       logical                        :: read_set, read_alphas

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
                read(unit=line(i+1:),fmt=*,iostat=ier) wave
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the wavelength"
                  return
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
                  ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading the distance sample-detector"
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

       if (read_set) then
          call Set_Current_Orient(wave,ub,set)
       else
          call Set_Current_Orient(wave,ub)
       end if

       if (.not. read_alphas) then
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
       integer                :: i,j,nl
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
          nl=5*(nval_f/5 + mod(nval_f,5))
          allocate(rvalues(nl))
          rvalues=0.0

          j=1
          do i=n_ini+ntext+2, n_end
             read(unit=filevar(i),fmt='(5E16.8)') vet
             rvalues(j:j+4)=vet
             j=j+5
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
       integer                :: i,j, nl
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
          nl=10*(nval_i/10 + mod(nval_i,10))
          allocate(ivalues(nl))
          ivalues=0

          j=1
          do i=n_ini+ntext+2, n_end
             read(unit=filevar(i),fmt='(10i8)') ivet
             ivalues(j:j+9)=ivet
             j=j+10
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
       integer                :: i,j,nl
       integer, dimension(10) :: ivet

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
          nl=8*(nval_i/8 + mod(nval_i,8))
          allocate(ivalues(nl))
          ivalues=0

          j=1
          do i=n_ini+ntext+2, n_end
             read(unit=filevar(i),fmt='(8i10)') ivet(1:8)
             ivalues(j:j+7)=ivet(1:8)
             j=j+8
          end do

       end if

       return
    End Subroutine Read_J_KeyType

    !!----
    !!---- Subroutine Read_Numor_D1B(filevar,N)
    !!----    character(len=*), intent(in)          :: fileinfo
    !!----    type(generic_numor_type), intent(out) :: n
    !!----
    !!---- Subroutine to read a Numor of D1B Instrument at ILL
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Read_Numor_D1B(fileinfo,N)
       !---- Arguments ----!
       character(len=*), intent(in)          :: fileinfo
       type(generic_numor_type), intent(out) :: n

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
          err_illdata_mess=' Problems trying to read the numor for D1B Instrument'
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
       call read_A_keyType(filevar,nl_keytypes(2,2,1),nl_keytypes(2,2,2),idum,line)
       if (idum > 0) then
          n%SampleID%n=4
          if (allocated(n%sampleid%namevar)) deallocate(n%sampleid%namevar)
          if (allocated(n%sampleid%cvalues)) deallocate(n%sampleid%cvalues)
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
       end if

       ! Control Flags
       call read_I_keyType(filevar,nl_keytypes(5,1,1),nl_keytypes(5,1,2))
       if (nval_i > 0) then
          n%dacflags%n=nval_i
          if (allocated(n%dacflags%namevar)) deallocate(n%dacflags%namevar)
          if (allocated(n%dacflags%ivalues)) deallocate(n%dacflags%ivalues)
          allocate(n%dacflags%namevar(nval_i))
          allocate(n%dacflags%ivalues(nval_i))

          ! Given names
          n%dacflags%namevar=' '
          n%dacflags%namevar( 5)='N. of additional parameters saved per data point'      ! Nbang (0-7)
          n%dacflags%namevar( 7)='Number of data points saved'                           ! Npdone (Nframes)
          n%dacflags%namevar( 8)='Count on Monitor or Time'                              ! Jcode (0/1)
          n%dacflags%namevar(24)='N. of count data'                                      ! Nbdata

          ! Load values
          n%dacflags%ivalues=ivalues(1:nval_i)
       end if

       ! Real values
       call read_F_keyType(filevar,nl_keytypes(4,1,1),nl_keytypes(4,1,2))
       if (nval_f > 0) then
          n%dacparam%n=nval_f
          if (allocated(n%dacparam%namevar)) deallocate(n%dacparam%namevar)
          if (allocated(n%dacparam%rvalues)) deallocate(n%dacparam%rvalues)
          allocate(n%dacparam%namevar(nval_f))
          allocate(n%dacparam%rvalues(nval_f))

          ! Given names
          n%dacparam%namevar=' '

          n%dacparam%namevar(18)='Wavelength'

          n%dacparam%namevar(36)='Scan Start'
          n%dacparam%namevar(37)='Scan Step'
          n%dacparam%namevar(38)='Scan Width'

          n%dacparam%namevar(39)='Preset monitor or Time'

          n%dacparam%namevar(46)='Temp-s.pt'
          n%dacparam%namevar(47)='Temp-Regul'
          n%dacparam%namevar(48)='Temp-Sample'

          ! Load values
          n%dacparam%rvalues=rvalues(1:nval_f)
       end if

       ! Counts
       call read_I_keyType(filevar,nl_keytypes(5,2,1),nl_keytypes(5,2,2))
       if (nval_i > 0) then
          n%icounts%n=nval_i
          if (allocated(n%icounts%namevar)) deallocate(n%icounts%namevar)
          if (allocated(n%icounts%ivalues)) deallocate(n%icounts%ivalues)
          allocate(n%icounts%namevar(3))
          allocate(n%icounts%ivalues(nval_i))

          ! Given Names
          n%icounts%namevar(1)='Monitor'
          n%icounts%namevar(2)='Time'
          n%icounts%namevar(3)='Initial Scan'

          ! Load values
          n%icounts%ivalues=ivalues(1:nval_i)
       end if

       return
    End Subroutine Read_Numor_D1B

    !!----
    !!---- Subroutine Read_Numor_D20(filevar,N)
    !!----    character(len=*), intent(in) :: fileinfo
    !!----    type(generic_numor_type), intent(out) :: n
    !!----
    !!---- Read a Numor for D20 Instrument at ILL
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Read_Numor_D20(fileinfo,N)
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
       print*, nval_i

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
    End Subroutine Read_Numor_D20

    !!----
    !!---- Subroutine Read_POWDER_Numor(Numor,Instrm,Pathdir,Snum,Info)
    !!----    integer,                 intent(in)    :: numor
    !!----    character(len=*),        intent(in)    :: instrm
    !!----    character(len=*),        intent(in)    :: pathdir
    !!----    type(POWDER_Numor_type), intent(in out):: snum
    !!----    logical, optional,       intent(in)    :: info
    !!----
    !!----    Read the numor "numor" from the ILL database and construct
    !!----    the object 'snum' of type POWDER_Numor_type.
    !!----    (not completely checked yet)
    !!----    In case of error the subroutine puts ERR_ILLData=.true.
    !!----    and fils the error message variable ERR_ILLData_Mess.
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Read_POWDER_Numor(numor,instrm,pathdir,snum,info)
       !---- Arguments ----!
       integer,                 intent(in)    :: numor
       character(len=*),        intent(in)    :: instrm
       character(len=*),        intent(in)    :: pathdir
       type(POWDER_Numor_type), intent(in out):: snum
       logical, optional,       intent(in)    :: info

       !---- Local variables ----!
       character(len=512)              :: filenam
       character(len=6)                :: inst
       character(len=8)                :: item
       character(len=80)               :: line
       integer                         :: i,j,n1,n2,n3, lun, ier
       integer, dimension(31)          :: ival
       real(kind=cp),    dimension(50) :: rval
       logical                         :: existe

       ! Instrument
       inst=adjustl(instrm)
       call lcase(inst)

       ! Construct the absolute path and filename to be read
       if (.not. got_ILL_data_directory  .or. len_trim(pathdir) == 0) then
          call Get_Absolute_Data_Path(numor,instrm,filenam)
       else
          write(unit=line,fmt="(i6.6)") numor
          i=len_trim(pathdir)
          if (pathdir(i:i) /= ops_sep) then
             filenam=trim(pathdir)//ops_sep//trim(line)
          else
             filenam=trim(pathdir)//trim(line)
          end if
       end if

       ! Check if the data exist uncompressed or compressed
       inquire(file=trim(filenam),exist=existe)
       if (.not. existe) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess=" Error: the file: "//trim(filenam)//", doesn't exist!"
          return
       end if

       call Get_LogUnit(lun)
       open(unit=lun,file=trim(filenam),status="old", action="read", position="rewind",iostat=ier)
       if (ier /= 0) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess=" Error opening the file: "//trim(filenam)//", for reading numors"
          return
       end if

       Select Case (inst)
          case("d1b")
             read(unit=lun,fmt="(a)",iostat=ier) line
             read(unit=lun,fmt=*,iostat=ier) snum%numor
             do i=1,3
                read(unit=lun,fmt="(a)",iostat=ier) line
             end do
             snum%header=line
             do i=1,3
                read(unit=lun,fmt="(a)",iostat=ier) line
             end do
             snum%title=line(1:32)
             do i=1,2
                read(unit=lun,fmt="(a)",iostat=ier) line
             end do
             read(unit=lun,fmt=*,iostat=ier) ival
             if (ier /= 0) then
                 ERR_ILLData=.true.
                 ERR_ILLData_Mess=" Error in file: "//trim(filenam)//", reading integer control values"
                 return
             end if
             do i=1,2
                read(unit=lun,fmt="(a)",iostat=ier) line
             end do
             read(unit=lun,fmt=*,iostat=ier) rval
             if (ier /= 0) then
                ERR_ILLData=.true.
                ERR_ILLData_Mess=" Error in file: "//trim(filenam)//", reading real control values"
                return
             end if
             snum%angles=rval(4:8)        !Phi, chi, omega,  2theta(gamma), psi      [valco(4:8)]
             snum%wave=rval(18)           ! wavelength             [valco(18)]
             snum%scans=rval(36:38)       ! start, step, end     [valdef(1:3)]
             snum%scans(2)=0.2
             snum%monitor=rval(39)        ! monitor               [valdef(4)]
             snum%conditions=rval(46:50)  ! Temp-s.pt, Temp-Regul, Temp-sample, Voltmeter, Mag.field [valenv(1:5))]
             snum%nbdata=ival(24)         ! nbdata
             snum%manip=ival(4)           ! manip
             snum%nbang=ival(5)           ! nbang
             snum%nframes=ival(7)         ! npdone
             snum%icalc=ival(9)           ! icalc
             snum%icdesc(1:7)=ival(25:31) ! icdesc
             if (present(info)) then
                if (info) then
                   do j=1,4
                      read(unit=lun,fmt="(a)",iostat=ier) line
                   end do
                   read(unit=lun,fmt="(a8,2i8)", iostat=ier) item,n2,n3
                   read(unit=item,fmt=*,iostat=ier) n1
                   if (ier /=0) n1=99999999
                   snum%scans(1)=real(n3)*0.001
                   snum%scans(3)= snum%scans(1) + snum%scans(2)*real(snum%nbdata-1)
                   snum%monitor=real(n1)
                   snum%time=real(n2)/60000.0
                   return
                end if
             end if
             if (allocated(snum%counts)) deallocate(snum%counts)
             allocate(snum%counts(ival(24),ival(7)))
             snum%counts=0.0
             if (allocated(snum%tmc_ang)) deallocate(snum%tmc_ang)
             allocate(snum%tmc_ang(ival(5)+3,ival(7)))   !normally five values, maximum nbang=2 angles
             snum%tmc_ang=0.0

             do i=1,snum%nframes
                do j=1,3
                   read(unit=lun,fmt="(a)",iostat=ier) line
                end do
                read(unit=lun,fmt=*,    iostat=ier) j      !read number of items to store in snum%counts
                if (j > snum%nbdata) then
                   read(unit=lun,fmt="(a8,2i8,7f8.0)", iostat=ier) item,n2,n3, snum%counts(1:7,i)
                   if (ier /= 0) then
                      ERR_ILLData=.true.
                      write(unit=ERR_ILLData_Mess,fmt="(a,i4)") &
                           " Error in file: "//trim(filenam)//", reading counts at frame #",i
                      return
                   end if
                   read(unit=lun,fmt=*, iostat=ier) snum%counts(8:snum%nbdata,i)
                   if (ier /= 0) then
                      ERR_ILLData=.true.
                      write(unit=ERR_ILLData_Mess,fmt="(a,i4)") &
                           " Error in file: "//trim(filenam)//", reading counts at frame #",i
                      return
                   end if

                   !Re-afect the values of n1,n2 and n3 to Monitor, Time and initial 2Theta
                   snum%scans(1)=real(n3)*0.001
                   snum%scans(3)= snum%scans(1) + snum%scans(2)*real(snum%nbdata-1)
                   read(unit=item,fmt=*,iostat=ier) n1
                   if (ier /=0) n1=99999999
                   snum%monitor=real(n1)
                   snum%time=real(n2)/60000.0
                else
                   read(unit=lun,fmt="(10f8.0)", iostat=ier) snum%counts(1:j,i)
                   if (ier /= 0) then
                      ERR_ILLData=.true.
                      write(unit=ERR_ILLData_Mess,fmt="(a,i4)") &
                           " Error in file: "//trim(filenam)//", reading counts at frame #",i
                      return
                   end if
                end if
             end do
             close(unit=lun)

          Case("d20")
             snum%scans(2)=0.1
             read(unit=lun,fmt="(a)",iostat=ier) line
             read(unit=lun,fmt=*,iostat=ier) snum%numor
             do i=1,4
                read(unit=lun,fmt="(a)",iostat=ier) line
             end do
             snum%header=line
             do i=1,3
                read(unit=lun,fmt="(a)",iostat=ier) line
             end do
             i=index(line,":")
             snum%header=trim(snum%header)//line(i+1:)
             do i=1,2
                read(unit=lun,fmt="(a)",iostat=ier) line
             end do
             read(unit=lun,fmt=*,iostat=ier) ival
             if (ier /= 0) then
                ERR_ILLData=.true.
                ERR_ILLData_Mess=" Error in file: "//trim(filenam)//", reading integer control values"
                return
             end if
             do i=1,2
                read(unit=lun,fmt="(a)",iostat=ier) line
             end do
             read(unit=lun,fmt=*,iostat=ier) rval
             if (ier /= 0) then
                ERR_ILLData=.true.
                ERR_ILLData_Mess=" Error in file: "//trim(filenam)//", reading real control values"
                return
             end if
             snum%angles=rval(4:8)        !Phi, chi, omega,  2theta(gamma), psi      [valco(4:8)]
             snum%wave=rval(18)           ! wavelength             [valco(18)]
             snum%scans=rval(36:38)       ! start, step, end     [valdef(1:3)]
             snum%scans(2)=0.2
             snum%monitor=rval(39)        ! monitor               [valdef(4)]
             snum%conditions=rval(46:50)  ! Temp-s.pt, Temp-Regul, Temp-sample, Voltmeter, Mag.field [valenv(1:5))]
             snum%nbdata=ival(24)         ! nbdata
             snum%manip=ival(4)           ! manip
             snum%nbang=ival(5)           ! nbang
             snum%nframes=ival(7)         ! npdone
             snum%icalc=ival(9)           ! icalc
             snum%icdesc(1:7)=ival(25:31) ! icdesc
             if (present(info)) then
                if (info) then
                   do j=1,4
                      read(unit=lun,fmt="(a)",iostat=ier) line
                   end do
                   read(unit=lun,fmt="(a8,2i8)", iostat=ier) item,n2,n3
                   read(unit=item,fmt=*,iostat=ier) n1
                   if(ier /=0) n1=99999999
                   snum%scans(1)=real(n3)*0.001
                   snum%scans(3)= snum%scans(1) + snum%scans(2)*real(snum%nbdata-1)
                   snum%monitor=real(n1)
                   snum%time=real(n2)/60000.0
                   return
                end if
             end if
             if (allocated(snum%counts)) deallocate(snum%counts)
             allocate(snum%counts(ival(24),ival(7)))
             snum%counts=0.0
             if (allocated(snum%tmc_ang)) deallocate(snum%tmc_ang)
             allocate(snum%tmc_ang(ival(5)+3,ival(7)))   !normally five values, maximum nbang=2 angles
             snum%tmc_ang=0.0

             do i=1,snum%nframes
                do j=1,3
                   read(unit=lun,fmt="(a)",iostat=ier) line
                end do
                read(unit=lun,fmt=*,    iostat=ier) j      !read number of items to store in snum%counts
                if (j > snum%nbdata) then
                   read(unit=lun,fmt="(a8,2i8,7f8.0)", iostat=ier) item,n2,n3, snum%counts(1:7,i)
                   if (ier /= 0) then
                      ERR_ILLData=.true.
                      write(unit=ERR_ILLData_Mess,fmt="(a,i4)") &
                           " Error in file: "//trim(filenam)//", reading counts at frame #",i
                      return
                   end if
                   read(unit=lun,fmt=*, iostat=ier) snum%counts(8:snum%nbdata,i)
                   if (ier /= 0) then
                      ERR_ILLData=.true.
                      write(unit=ERR_ILLData_Mess,fmt="(a,i4)") &
                           " Error in file: "//trim(filenam)//", reading counts at frame #",i
                      return
                   end if
                   ! Re-afect the values of n1,n2 and n3 to Monitor, Time and initial 2Theta
                   snum%scans(1)=real(n3)*0.001
                   snum%scans(3)= snum%scans(1) + snum%scans(2)*real(snum%nbdata-1)
                   read(unit=item,fmt=*,iostat=ier) n1
                   if(ier /=0) n1=99999999
                   snum%monitor=real(n1)
                   snum%time=real(n2)/60000.0
               else
                  read(unit=lun,fmt="(10f8.0)", iostat=ier) snum%counts(1:j,i)
                  if (ier /= 0) then
                     ERR_ILLData=.true.
                     write(unit=ERR_ILLData_Mess,fmt="(a,i4)") &
                          " Error in file: "//trim(filenam)//", reading counts at frame #",i
                     return
                  end if
               end if
            end do
            close(unit=lun)

         Case("d1a")
         Case("d2b")
         Case("d4c")

         Case default
            ERR_ILLData=.true.
            ERR_ILLData_Mess= " Error in file: "//trim(filenam)//", Incorrect Instrument name: "//inst
       End Select

       return
    End Subroutine Read_POWDER_Numor

    !!--++
    !!--++ Subroutine Read_R_KeyType(filevar, n_ini, n_end, nrun, nvers)
    !!--++    character(len=*), dimension(:), intent(in)  :: filevar   ! Input information
    !!--++    integer,                        intent(in)  :: n_ini     ! First line
    !!--++    integer,                        intent(in)  :: n_end     ! Last line to be read
    !!--++    integer,                        intent(out) :: nrun      ! Run number of the data
    !!--++    integer,                        intent(out) :: nvers     ! Version data
    !!--++
    !!--++ (Private)
    !!--++ Read the Blocktype R on Numor at ILL
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
       nrun =ivet(1)
       nvers=ivet(3)
       ntext=ivet(2)

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

    !!----
    !!---- Subroutine Read_SXTAL_Numor(Numor,Instr,Snum)
    !!----    integer,                intent(in)    :: numor
    !!----    character(len=*),       intent(in)    :: Instr
    !!----    type(SXTAL_Numor_type), intent(in out):: snum
    !!----
    !!----    Read the numor "numor" from the ILL database and construct
    !!----    the object 'snum' of type SXTAL_Numor_type.
    !!----    (not completely checked yet)
    !!----    In case of error the subroutine puts ERR_ILLData=.true.
    !!----    and fils the error message variable ERR_ILLData_Mess.
    !!----
    !!---- Update: December - 2005
    !!
    Subroutine Read_SXTAL_Numor(numor,Instrm,snum)
       !---- Arguments ----!
       integer,                intent(in)    :: numor
       character(len=*),       intent(in)    :: Instrm
       type(SXTAL_Numor_type), intent(in out):: snum

       !---- Local variables ----!
       character(len=280)              :: filenam
       character(len=6)                :: inst
       character(len=80)               :: line
       integer                         :: i,j, lun, long,ier
       integer, dimension(31)          :: ival
       real(kind=cp),    dimension(50) :: rval
       logical                         :: existe


       !---- Construct the absolute path and filename to be read ----!
       write(unit=inst,fmt='(i6.6)') numor
       line=trim(Instrm_directory)
       long=len_trim(line)
       if (line(long:long) /= ops_sep) line=trim(line)//ops_sep
       filenam=trim(line)//inst

       !---- Check if the data exist uncompressed or compressed
       inquire(file=trim(filenam),exist=existe)
       if (.not. existe) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess="The file: "//trim(filenam)//", doesn't exist!"
          return
       end if

       call Get_LogUnit(lun)
       open(unit=lun,file=trim(filenam),status="old", action="read", position="rewind",iostat=ier)
       if (ier /= 0) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess="Error opening the file: "//trim(filenam)//", for reading numors"
          return
       end if

       inst=instrm
       call lcase (inst)

       Select Case (trim(inst))
          Case("d9","d10","d19")
              read(unit=lun,fmt="(a)",iostat=ier) line
              read(unit=lun,fmt=*,iostat=ier) snum%numor
              do i=1,4
                read(unit=lun,fmt="(a)",iostat=ier) line
              end do
              snum%header=line
              do i=1,4
                read(unit=lun,fmt="(a)",iostat=ier) line
              end do
              snum%title=line(1:32)
              snum%Scantype=adjustl(line(71:80))
              do i=1,6
                read(unit=lun,fmt="(a)",iostat=ier) line
              end do
              read(unit=lun,fmt=*,iostat=ier) ival
              if(ier /= 0) then
                ERR_ILLData=.true.
                ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading integer control values"
                return
              end if
              do i=1,12
                read(unit=lun,fmt="(a)",iostat=ier) line
              end do
              read(unit=lun,fmt=*,iostat=ier) rval
              if(ier /= 0) then
                ERR_ILLData=.true.
                ERR_ILLData_Mess="Error in file: "//trim(filenam)//", reading real control values"
                return
              end if
              !valco      rval( 1:35)
              !valdef     rval(36:45)
              !valenv     rval(46:50)
              snum%hmin=rval(1:3)          !h,k,l          [valco( 1: 3)]
              snum%hmax=rval(22:24)        !h,k,l (max)    [valco(22:24)]
              snum%dh=rval(25:27)          !dh,dk,dl       [valco(25:27)]
              snum%angles=rval(4:8)        !Phi, chi, omega,  2theta(gamma), psi      [valco(4:8)]
              snum%ub(1,:)=rval(9:11)      !
              snum%ub(2,:)=rval(12:14)     ! UB matrix              [valco(9:17)]
              snum%ub(3,:)=rval(15:17)     !
              snum%wave=rval(18)           ! wavelength             [valco(18)]
              snum%cpl_fact=rval(43)       ! Coupling Factor        [valdef(8)]
              snum%scans=rval(36:38)       ! start, step, width     [valdef(1:3)]
              snum%preset=rval(39)         ! preset                 [valdef(4)]
              snum%conditions=rval(46:50)  ! Temp-s.pt, Temp-Regul, Temp-sample, Voltmeter, Mag.field [valenv(1:5))]
              snum%nbdata=ival(24)         ! nbdata
              snum%manip=ival(4)           ! manip
              snum%nbang=ival(5)           ! nbang
              snum%nframes=ival(7)         ! npdone
              snum%icalc=ival(9)           ! icalc
              snum%icdesc(1:7)=ival(25:31) ! icdesc
              if(allocated(snum%counts)) deallocate(snum%counts)
              allocate(snum%counts(ival(24),ival(7)))
              snum%counts=0.0
              if(allocated(snum%tmc_ang)) deallocate(snum%tmc_ang)
              allocate(snum%tmc_ang(ival(5)+3,ival(7)))   !normally five values, maximum nbang=2 angles
              snum%tmc_ang=0.0

              do i=1,snum%nframes
                do j=1,3
                   read(unit=lun,fmt="(a)",iostat=ier) line
                end do
                read(unit=lun,fmt=*,    iostat=ier) j      !read number of items to store in snum%tmc_ang
                read(unit=lun,fmt="(a)",iostat=ier) line
                read(unit=lun,fmt=*,    iostat=ier) snum%tmc_ang(1:j,i)
                snum%tmc_ang(4:j,i)= snum%tmc_ang(4:j,i)*0.001 !angles are stored here as milidegrees
                read(unit=lun,fmt="(a)",iostat=ier) line
                read(unit=lun,fmt=*,    iostat=ier) j      !read number of items to store in snum%counts
                read(unit=lun,fmt="(10f8.0)", iostat=ier) snum%counts(1:j,i)
                if(ier /= 0) then
                  ERR_ILLData=.true.
                  write(unit=ERR_ILLData_Mess,fmt="(a,i4)") &
                     "Error in file: "//trim(filenam)//", reading counts at frame #",i
                  return
                end if
              end do
              close(unit=lun)

          Case default
              ERR_ILLData=.true.
              ERR_ILLData_Mess= "Error in file: "//trim(filenam)//", Incorrect Instrument name: "//inst
       End Select

       return
    End Subroutine Read_SXTAL_Numor

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
    !!---- Subroutine Set_Default_Instrument()
    !!----
    !!----    Construct the Current_Instrument as a default 4C diffractometer
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
    !!---- Subroutine Set_ILL_Data_Directory(ILL_data_try)
    !!----    character(len=*), optional, intent(in) :: ILL_data_try  !proposed location of ILL data
    !!----
    !!----    Assign the global public variable: ILL_data_directory
    !!----    If ILL_data_try=' ' data are in the current directory
    !!----    If invoked without argument, the subroutine asks for the environment
    !!----    variable ILLDATA. If it is defined ILL_data_directory=ILLDATA,
    !!----    if not ILL_data_directory takes the default value defined in the
    !!----    declaration of the variable.
    !!----    If the directory doesn't exist the subroutine rises an error condition
    !!----    by putting ERR_ILLData=.true. and filling the error message
    !!----    variable ERR_ILLData_Mess.
    !!----
    !!---- Update: December - 2006
    !!
    Subroutine Set_ILL_Data_Directory(ILL_data_try)
       !---- Arguments ----!
       character(len=*), optional, intent(in) :: ILL_data_try

       !---- Local Variables ----!
       character(len=512) :: ILL_data
       integer            :: len_d
       logical            :: existe

       if (present(ILL_data_try)) then
          if (len_trim(ILL_data_try) == 0) then
             ILL_data_directory =" "      !data are in the current directory
             got_ILL_data_directory=.true.
             return
          else
             ILL_data_directory =trim(ILL_data_try)
          end if

       else
          call get_environment_variable("ILLDATA", ILL_data)
          if (len_trim(ILL_data) > 0) then
             ILL_data_directory = trim(ILL_data)
          else
             select case (ops)
                case (1) ! Windows
                   ILL_data_directory ='\\Serdon\illdata'
                case (2:)
                   ILL_data_directory ='/net/serdon/illdata'
             end select
          end if
       end if

       !---- Add separator if absent ----!
       len_d=len_trim(ILL_data_directory)
       if (ILL_data_directory(len_d:len_d) /= ops_sep) then
         ILL_data_directory=trim(ILL_data_directory)//ops_sep
       end if

       !---- Check that the directory exist, ----!
       !---- otherwise rise an error condition ----!
       err_ILLData=.false.
       existe=directory_exists(trim(ILL_data_directory))

       if (.not. existe) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess="The ILL directory: '"//trim(ILL_data_directory)//"' doesn't exist"
          got_ILL_data_directory=.false.
       else
          got_ILL_data_directory=.true.
       end if

       return
    End Subroutine Set_ILL_Data_Directory

    !!----
    !!---- Subroutine Set_Instrm_directory(instrm)
    !!----    character(len=*),  intent(in) :: instrm  !Name of the instrument
    !!----
    !!----    Assign the global public variable: Instrm_directory
    !!----    The intent 'in' variable instrm is the name of the diffractometer
    !!----    Instrm_directory=trim(ILL_data_directory)//trim(instrm)//ops_sep
    !!----    It is assumed that the subroutine Set_ILL_data_directory has
    !!----    already been called.
    !!----    If the directory doesn't exist the subroutine rises an error condition
    !!----    by putting ERR_ILLData=.true. and filling the error message
    !!----    variable ERR_ILLData_Mess.
    !!----
    !!---- Update: December - 2006
    !!
    Subroutine Set_Instrm_directory(instrm)
       !---- Argument ----!
       character(len=*),  intent(in) :: instrm  !Name of the instrument

       !---- Local Variables ----!
       logical :: existe

       Instrm_directory=trim(ILL_data_directory)//'data'//ops_sep//trim(instrm)//ops_sep

       !---- check that the directory exist, ----!
       !---- otherwise rise an error condition ----!
       ERR_ILLData=.false.
       existe=directory_exists(trim(Instrm_directory))
       if (.not. existe) then
          ERR_ILLData=.true.
          ERR_ILLData_Mess="The INSTRM directory: '"//trim(Instrm_directory)//"' doesn't exist"
          Instrm_directory_set=.false.
       end if

       Instrm_directory_set=.true.

       return
    End Subroutine Set_Instrm_directory

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
       character(len=10) :: forma


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
          write(unit=ipr,fmt="(a,f8.4,a)")   "  WAVELENGTH:",Current_Orient%wave," angstroms"

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
          write(unit=ipr,fmt="(a,2f8.3,2i5)") "DIM_XY ", Current_Instrm%horiz,    Current_Instrm%vert, &
                                                    Current_Instrm%np_horiz, Current_Instrm%np_vert
          write(unit=ipr,fmt="(a,2f8.4)") "GAPS_DET ",Current_Instrm%agap, Current_Instrm%cgap
          write(unit=ipr,fmt="(a,f8.5)") "WAVE ",Current_Orient%wave
          write(unit=ipr,fmt="(a)") "UBMAT"
          do i=1,3
             write(unit=ipr,fmt="(3f12.7)") Current_Orient%ub(i,:)
          end do
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

          write(unit=ipr,fmt="(a)") "DET_ALPHAS"
          forma="(   f8.4)"
          write(unit=forma(2:4),fmt="(i3)") Current_Instrm%np_horiz
          do i=1,Current_Instrm%np_vert
             write(unit=ipr,fmt=forma) Current_Instrm%alphas(i,1:Current_Instrm%np_horiz)
          end do
       end if

       return
    End Subroutine Write_Current_Instrm_data

    !!----
    !!---- Subroutine Write_Generic_Numor(N,lun)
    !!----    type(generic_numor_type), intent(in) :: N
    !!----    integer, optional,        intent(in) :: lun
    !!----
    !!---- Write the information of a Generic Numor
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Write_Generic_Numor(N,lun)
       !---- Arguments ----!
       type(generic_numor_type), intent(in) :: N
       integer, optional,        intent(in) :: lun

       !---- Local Variables ----!
       integer :: ilun,i

       ilun=6
       if (present(lun)) ilun=lun

       write(unit=ilun, fmt='(a)') '#### Numor Information ####'
       write(unit=ilun, fmt='(a,i6.6)') 'Numor: ',n%numor

       write(unit=ilun,fmt='(a,t50,a)') 'Instrument: '//trim(n%instr),'Date: '//trim(n%date)
       write(unit=ilun,fmt='(a)') 'Experimental Name: '//trim(n%expname)
       write(unit=ilun,fmt='(a)') ' '

       write(unit=ilun,fmt='(a)') '#---- Sample Information ----#'
       do i=1,n%sampleid%n
          write(unit=ilun,fmt='(a)') trim(n%sampleid%namevar(i))//': '//trim(n%sampleid%cvalues(i))
       end do
       write(unit=ilun,fmt='(a)') ' '

       select case (n%instr)
           case ('D20')
              write(unit=ilun,fmt='(a)') '#---- Diffractometer Optics and Reactor Parameters ----#'
              do i=1,n%diffopt%n
                 if (len_trim(n%diffopt%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(n%diffopt%namevar(i))//': ',n%diffopt%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Monochromator Motor Parameters ----#'
              do i=1,n%monmpar%n
                 if (len_trim(n%monmpar%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(n%monmpar%namevar(i))//': ',n%monmpar%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Diffractometer Motor Parameters ----#'
              do i=1,n%diffmpar%n
                 if (len_trim(n%diffmpar%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(n%diffmpar%namevar(i))//': ',n%diffmpar%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Detector and DAS Parameters ----#'
              do i=1,n%detpar%n
                 if (len_trim(n%detpar%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(n%detpar%namevar(i))//': ',n%detpar%rvalues(i)
              end do
             write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Data Acquisition Parameters ----#'
              do i=1,n%dacparam%n
                 if (len_trim(n%dacparam%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(n%dacparam%namevar(i))//': ',n%dacparam%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Sample Status ----#'
              do i=1,n%samplest%n
                 if (len_trim(n%samplest%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.8)') trim(n%samplest%namevar(i))//': ',n%samplest%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Counts Information ----#'
              write(unit=ilun,fmt='(a,i8)') trim(n%icounts%namevar(1))//': ',n%icounts%n
              write(unit=ilun,fmt='(a)') ' '

           case ('D1B')
              write(unit=ilun,fmt='(a)') '#---- Data Acquisition Flags ----#'
              do i=1,n%dacflags%n
                 if (len_trim(n%dacflags%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,i6)') trim(n%dacflags%namevar(i))//': ',n%dacflags%ivalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Data Acquisition Parameters ----#'
              do i=1,n%dacparam%n
                 if (len_trim(n%dacparam%namevar(i)) <= 0) cycle
                 write(unit=ilun,fmt='(a,g15.4)') trim(n%dacparam%namevar(i))//': ',n%dacparam%rvalues(i)
              end do
              write(unit=ilun,fmt='(a)') ' '

              write(unit=ilun,fmt='(a)') '#---- Counts Information ----#'
              write(unit=ilun,fmt='(a,g15.4)') trim(n%icounts%namevar(1))//': ',real(n%icounts%ivalues(1))
              write(unit=ilun,fmt='(a,g15.4)') trim(n%icounts%namevar(2))//': ',real(n%icounts%ivalues(2))/60000.0
              write(unit=ilun,fmt='(a,g15.4)') trim(n%icounts%namevar(3))//': ',real(n%icounts%ivalues(3))*0.001
              write(unit=ilun,fmt='(a)') ' '

       end select

       return
    End Subroutine Write_Generic_Numor

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
    Subroutine Write_POWDER_Numor(Num,inst,lun)
       !---- Arguments ----!
       type(POWDER_Numor_type), intent(in) :: Num
       character(len=*),        intent(in) :: inst
       integer,       optional, intent(in) :: lun

       !--- Local variables ---!
       integer       :: ipr, i, j, nlines,n
       integer, dimension(10) :: cou

       ipr=6
       if (present(lun)) ipr=lun

       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a,i6,a)") " ---  Information about the NUMOR: ",Num%numor," of instrument "//trim(inst)
       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)")        "         HEADER: "//Num%header
       write(unit=ipr,fmt="(a)")        "          TITLE: "//Num%title
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
       write(unit=ipr,fmt="(a,i6,a)") " ---  Information about the NUMOR: ",Num%numor," of instrument "//Current_Instrm%name_inst
       write(unit=ipr,fmt="(a)")      " ---------------------------------------------------------------------------"
       write(unit=ipr,fmt="(a)") " "
       write(unit=ipr,fmt="(a)")        "         HEADER: "//Num%header
       write(unit=ipr,fmt="(a)")        "          TITLE: "//Num%title
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

 End Module CFML_ILL_Instrm_Data
