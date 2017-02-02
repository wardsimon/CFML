
!     Last change:  TR   30 Aug 2006    2:20 pm
!
!----------------------------------------------------------------------
module atomic_data
 USE macros_module, ONLY       : u_case
 use atome_module,  only       : atom

 contains

 subroutine definition_atomic_label()
 implicit none
  INTEGER          :: i

 ! symboles
 atom(1:201)%symbol = ""
 atom( 1:10)%symbol   = (/"H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne"/)
 atom( 11:20)%symbol  = (/"Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca"/)
 atom( 21:30)%symbol  = (/"Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"/)
 atom( 31:40)%symbol  = (/"Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr"/)
 atom( 41:50)%symbol  = (/"Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn"/)
 atom( 51:60)%symbol  = (/"Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd"/)
 atom( 61:70)%symbol  = (/"Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"/)
 atom( 71:80)%symbol  = (/"Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg"/)
 atom( 81:90)%symbol  = (/"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th"/)
 atom( 91:100)%symbol = (/"Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm"/)
 atom(101:103)%symbol = (/"Md", "No", "Lw"/)
 atom(201:201)%symbol = (/"D "/)

 do i=1, 201
  atom(i)%symbol = u_case(atom(i)%symbol)
 end do

 ! nom des atomes           1234567890123     1234567890123    1234567890123     1234567890123    1234567890123     1234567890123    1234567890123    1234567890123    1234567890123    1234567890123
 atom(  1:201)%name    = ""
 atom(  1: 10)%name    = (/"Hydrogene    ",  "Helium       ", "Lithium      ",  "Beryllium    ", "Bore         ", &
                           "Carbone      ",  "Azote        ", "Oxygene      ",  "Fluor        ", "Neon         "/)
 atom( 11: 20)%name    = (/"Sodium       ",  "Magnesium    ", "Aluminium    ",  "Silicium     ", "Phosphore    ", &
                           "Soufre       ",  "Chlore       ", "Argon        ",  "Potassium    ", "Calcium      "/)
 atom( 21: 30)%name    = (/"Scandium     ",  "Titane       ", "Vanadium     ",  "Chrome       ", "Manganese    ", &
                           "Fer          ",  "Cobalt       ", "Nickel       ",  "Cuivre       ", "Zinc         "/)
 atom( 31: 40)%name    = (/"Gallium      ",  "Germanium    ", "Arsenic      ",  "Selenium     ", "Brome        ", &
                           "Krypton      ",  "Rubidium     ", "Strontium    ",  "Yttrium      ", "Zirconium    "/)
 atom( 41: 50)%name    = (/"Niobium      ",  "Molybdene    ", "Technetium   ",  "Ruthenium    ", "Rhodium      ", &
                           "Palladium    ",  "Argent       ", "Cadmium      ",  "Indium       ", "Etain        "/)
 atom( 51: 60)%name    = (/"Antimoine    ",  "Tellure      ", "Iode         ",  "Xenon        ", "Cesium       ", &
                           "Barium       ",  "Lanthane     ", "Cerium       ",  "Praseodyme   ", "Neodyme      "/)
 atom( 61: 70)%name    = (/"Promethium   ",  "Samarium     ", "Europium     ",  "Gadolinium   ", "Terbium      ", &
                           "Dysprosium   ",  "Holmium      ", "Erbium       ",  "Thullium     ", "Ytterbium    "/)
 atom( 71: 80)%name    = (/"Lutetium     ",  "Hafnium      ", "Tantale      ",  "Tungtene     ", "Rhenium      ", &
                           "Osmium       ",  "Iridium      ", "Platine      ",  "Or           ", "Mercure      "/)
 atom( 81: 90)%name    = (/"Thallium     ",  "Plomb        ", "Bismuth      ",  "Polonium     ", "Astate       ", &
                           "Radon        ",  "Francium     ", "Radium       ",  "Actinum      ", "Thorium      "/)
 atom( 91:100)%name    = (/"Protactinium ",  "Uranium      ", "Neptunium    ",  "Plutonium    ", "Americium    ", &
                           "Curium       ",  "Berkelium    ", "Californium  ",  "Einsteinium  ", "Fermium      "/)
 atom(101:110)%name    = (/"Mendelevium  ",  "Nobelium     ", "Lawrencium   ",  "Rutherfordium", "Dubnium      ", &
                           "Seaborgium   ",  "Bohrium      ", "Hassium      ",  "Meitneriume  ", "Ununnilium   "/)
 atom(201)%NAME        = "Deuterium"

 end subroutine definition_atomic_label

 !------------------------------------------
 subroutine def_atomic_features()

!-----------------------------------------------------------------
! masse molaire

 atom(  1:201)%weight = 0.
 atom(  1: 10)%weight = (/ 1.0079,    4.0026,  6.941,     9.01218, 10.81,    12.011, 14.0067,  15.9994, 18.998403, 20.179/)
 atom( 11: 20)%weight = (/22.98977,  24.305,  26.98154,  28.0855,  30.97376, 32.06,  35.453,   39.9484, 39.0983,   40.08/)
 atom( 21: 30)%weight = (/44.9559 ,  47.905,  50.9415,   51.996,   54.938,   55.847, 58.9332,  58.7,    63.546,    65.38/)
 atom( 31: 40)%weight = (/69.72,     72.59,   74.9216,   78.96,    79.904,   83.8,   85.4678,  87.62,   88.9059,   91.22/)
 atom( 41: 50)%weight = (/92.9064,   95.94,   98.,      101.07,   102.9055, 106.4,  107.868,  112.41,  114.82,    118.69/)
 atom( 51: 60)%weight = (/121.75,   127.6,   126.9045,  131.3,    132.9024, 137.33, 138.9055, 140.12,  140.9077,  144.24/)
 atom( 61: 70)%weight = (/145.,     150.4,   151.96,    157.25,   158.9254, 162.5,  164.9304, 167.26,  168.9342,  173.04/)
 atom( 71: 80)%weight = (/174.967,  178.49,  180.9479,  183.85,   186.207,  190.2,  192.22,   195.09,  196.9665,  200.59/)
 atom( 81: 90)%weight = (/204.37,   207.2,   208.9804,  209.,     210.,     222.,   223.,     226.,    227.0278,  232.0381/)
 atom( 91:100)%weight = (/231.0359, 238.029, 237.0482,  244.,     243.,     247.,   247.,     251.,    252.,      257./)
 atom(101:103)%weight = (/258.,     259.,    260./)
 atom(201:201)%weight = (/2.0158/)

 !
 atom(1:201)%Z = 0
 atom(  1: 10)%Z = (/  1,   2,   3,   4,   5,   6,   7,   8,   9,  10/)
 atom( 11: 20)%Z = (/ 11,  12,  13,  14,  15,  16,  17,  18,  19,  20/)
 atom( 21: 30)%Z = (/ 21,  22,  23,  24,  25,  26,  27,  28,  29,  30/)
 atom( 31: 40)%Z = (/ 31,  32,  33,  34,  35,  36,  37,  38,  39,  40/)
 atom( 41: 50)%Z = (/ 41,  42,  43,  44,  45,  46,  47,  48,  49,  50/)
 atom( 51: 60)%Z = (/ 51,  52,  53,  54,  55,  56,  57,  58,  59,  60/)
 atom( 61: 70)%Z = (/ 61,  62,  63,  64,  65,  66,  67,  68,  69,  70/)
 atom( 71: 80)%Z = (/ 71,  72,  73,  74,  75,  76,  77,  78,  79,  80/)
 atom( 81: 90)%Z = (/ 81,  82,  83,  84,  85,  86,  87,  88,  89,  90/)
 atom( 91:100)%Z = (/ 91,  92,  93,  94,  95,  96,  97,  98,  99, 100/)
 atom(101:103)%Z = (/101, 102, 103/)
 atom(201:201)%Z = (/  1/)

 ! configuration electronique
 atom(1:201)%config_electr   = ""
!                               1234567890123456789012    1234567890123456789012     1234567890123456789012     1234567890123456789012    1234567890123456789012
 atom(1:10)%config_electr    = (/"1s1                   ", "1s2                   ",  "1s2 2s1               ", &
                                 "1s2 2s2               ", "1s2 2s2 2p1           ",  "1s2 2s2 2p2           ", &
                                 "1s2 2s2 2p3           ", "1s2 2s2 2p4           ",  "1s2 2s2 2p5           ", &
                                 "1s2 2s2 2p6           "/)
 atom(11:20)%config_electr   = (/"[Ne] 3s1              ", "[Ne] 3s2              ",  "[Ne] 3s2 3p1          ", &
                                 "[Ne] 3s2 3p2          ", "[Ne] 3s2 3p3          ",  "[Ne] 3s2 3p4          ", &
                                 "[Ne] 3s2 3p5          ", "[Ne] 3s2 3p6          ",  "[Ar] 4s1              ", &
                                 "[Ar] 4s2              "/)
 atom(21:30)%config_electr   = (/"[Ar] 3d1 4s2          ", "[Ar] 3d2 4s2          ",  "[Ar] 3d3 4s2          ", &
                                 "[Ar] 3d5 4s1          ", "[Ar] 3d5 4s2          ",  "[Ar] 3d6 4s2          ", &
                                 "[Ar] 3d7 4s2          ", "[Ar] 3d8 4s2          ",  "[Ar] 3d10 4s1         ", &
                                 "[Ar] 3d10 4s2         "/)
 atom(31:40)%config_electr   = (/"[Ar] 3d10 4s2 4p1     ", "[Ar] 3d10 4s2 4p2     ",  "[Ar] 3d10 4s2 4p3     ", &
                                 "[Ar] 3d10 4s2 4p4     ", "[Ar] 3d10 4s2 4p5     ",  "[Ar] 3d10 4s2 4p6     ", &
                                 "[Kr] 5s1              ", "[Kr] 5s2              " , "[Kr] 4d1 5s2          ", &
                                 "[Kr] 4d2 5s2          "/)
 atom(41:50)%config_electr   = (/"[Kr] 4d4 5s1          ", "[Kr] 4d5 5s1          ",  "[Kr] 4d5 5s2          ", &
                                 "[Kr] 4d7 5s1          ", "[Kr] 4d8 5s1          ",  "[Kr] 4d10             ", &
                                 "[Kr] 4d10 5s1         ", "[Kr] 4d10 5s2         ",  "[Kr] 4d10 5s2 5p1     ", &
                                 "[Kr] 4d10 5s2 5p2     "/)
 atom(51:60)%config_electr   = (/"[Kr] 4d10 5s2 5p3     ", "[Kr] 4d10 5s2 5p4     ",  "[Kr] 4d10 5s2 5p5     ", &
                                 "[Kr] 4d10 5s2 5p6     ", "[Xe] 6s1              ",  "[Xe] 6s2              ", &
                                 "[Xe] 5d1 6s2          ", "[Xe] 4f1 5d1 6s2      ",  "[Xe] 4f3 5d0 6s2      ", &
                                 "[Xe] 4f4 5d0 6s2      "/)
 atom(61:70)%config_electr   = (/"[Xe] 4f5 5d0 6s2      ", "[Xe] 4f6 5d0 6s2      ",  "[Xe] 4f7 5d0 6s2      ", &
                                 "[Xe] 4f7 5d1 6s2      ", "[Xe] 4f9 5d0 6s2      ",  "[Xe] 4f10 5d0 6s2     ", &
                                 "[Xe] 4f11 5d0 6s2     ", "[Xe] 4f12 5d0 6s2     ",  "[Xe] 4f13 5d0 6s2     ", &
                                 "[Xe] 4f14 5d0 6s2     "/)
 atom(71:80)%config_electr   = (/"[Xe] 4f14 5d1 6s2     ", "[Xe] 4f14 5d2 6s2     ",  "[Xe] 4f14 5d3 6s2     ", &
                                 "[Xe] 4f14 5d4 6s2     ", "[Xe] 4f14 5d5 6s2     ",  "[Xe] 4f14 5d6 6s2     ", &
                                 "[Xe] 4f14 5d7 6s2     ", "[Xe] 4f14 5d9 6s1     ",  "[Xe] 4f14 5d10 6s1    ", &
                                 "[Xe] 4f14 5d10 6s2    "/)
 atom(81:90)%config_electr   = (/"[Xe] 4f14 5d10 6s2 6p1", "[Xe] 4f14 5d10 6s2 6p2",  "[Xe] 4f14 5d10 6s2 6p3", &
                                 "[Xe] 4f14 5d10 6s2 6p4", "[Xe] 4f14 5d10 6s2 6p5",  "[Xe] 4f14 5d10 6s2 6p6", &
                                 "[Rn] 7s1              ", "[Rn] 7s2              ",  "[Rn] 6d2 7s1          " ,&
                                 "[Rn] 5f0 6d2 7s2      "/)
 atom(91:100)%config_electr  = (/"[Rn] 5f2 6d1 7s2      ", "[Rn] 5f3 6d1 7s2      ",  "[Rn] 5f4 6d1 7s2      ", &
                                 "[Rn] 5f6 6d0 7s2      ", "[Rn] 5f7 6d0 7s2      ",  "[Rn] 5f7 6d1 7s2      ", &
                                 "[Rn] 5f8 6d1 7s2      ", "[Rn] 5f10 6d0 7s2     ",  "[Rn] 5f11 6d0 7s2     ", &
                                 "[Rn] 5f12 6d0 7s2     "/)
 atom(101:110)%config_electr = (/"[Rn] 5f13 + 6d0 7s2   ", "[Rn] 5f14 + 6d0 7s2   ",  "[Rn] 5f14 + 6d1 7s2   ", &
                                 "[Rn] 5f14 + 6d2 7s2   ", "[Rn] 5f14 + 6d3 7s2 ? ",  "[Rn] 5f14 + 6d4 7s2 ? ", &
                                 "[Rn] 5f14 + 6d5 7s2 ? ", "[Rn] 5f14 + 6d6 7s2 ? " , "[Rn] 5f14 + 6d7 7s2 ? ", &
                                  "[Rn] 5f14 + 6d8 7s2 ? "/)
 atom(201)%config_electr     = "1s1"


! desite atomique
 atom( 1:201)%density = 0.
 atom(  1:10)%density = (/   .0899,   .1787,    .53,   1.85,   2.34,   2.62 ,  1.251,   1.429,   1.696,   .901 /)
 atom( 11:20)%density = (/   .97  ,  1.74  ,   2.7 ,   2.33,   1.82,   2.07 ,  3.17 ,   1.784,    .86 ,  1.55  /)
 atom( 21:30)%density = (/  3.    ,  4.5   ,   5.8 ,   7.19,   7.43,   7.86 ,  8.9  ,   8.9  ,   8.96 ,  7.14  /)
 atom( 31:40)%density = (/  5.91  ,  5.32  ,   5.72,   4.8 ,   3.12,   3.74 ,  1.53 ,   2.6  ,   4.5  ,  6.49  /)
 atom( 41:50)%density = (/  8.55  , 10.2   ,  11.5 ,  12.2 ,  12.4 ,  12.   , 10.5  ,   8.65 ,   7.31 ,  7.3   /)
 atom( 51:60)%density = (/  6.68  ,  6.24  ,   4.92,   5.89,   1.87,   3.5  ,  6.7  ,   6.78 ,   6.77 ,  7.    /)
 atom( 61:70)%density = (/  6.475 ,  7.54  ,   5.26,   7.89,   8.27,   8.54 ,  8.8  ,   9.05 ,   9.33 ,  6.98  /)
 atom( 71:80)%density = (/  9.84  , 13.1   ,  16.6 ,  19.3 ,  21.  ,  22.4  , 22.5  ,  21.4  ,  19.3  , 13.53  /)
 atom( 81:90)%density = (/ 11.85  , 11.4   ,   9.8 ,   9.4 ,   0.  ,   9.91 ,  0.   ,   5.   ,  10.07 ,  11.7  /)
 atom( 91:96)%density = (/ 15.4   , 18.9   ,  20.4 ,  19.8 ,  13.6 ,  13.511/)


 ! rayon atomique:
 atom(1:201)%radius = 0.
 atom( 1:10)%radius =  (/0.78,  1.28,  1.52,  1.13,  0.83,  0.65,  0.71,  0.47,  0.71,  0.36 /)
 atom(11:20)%radius =  (/1.54,  1.60,  1.43,  1.17,  0.93,  1.04,  0.78,  1.74,  2.27,  1.97 /)
 atom(21:30)%radius =  (/1.61,  1.44,  1.32,  1.25,  1.24,  1.24,  1.25,  1.24,  1.28,  1.33 /)
 atom(31:40)%radius =  (/1.22,  1.23,  1.21,  2.15,  0.96,  0.88,  2.48,  2.15,  1.81,  1.60 /)
 atom(41:50)%radius =  (/1.43,  1.36,  1.36,  1.34,  1.34,  1.38,  1.44,  1.49,  1.63,  1.40 /)
 atom(51:60)%radius =  (/1.82,  1.43,  1.12,  2.18,  2.65,  2.17,  1.88,  1.83,  1.83,  1.82 /)
 atom(61:70)%radius =  (/1.81,  1.80,  2.04,  1.80,  1.78,  1.77,  1.77,  1.76,  1.75,  1.94 /)
 atom(71:80)%radius =  (/1.73,  1.56,  1.43,  1.37,  1.38,  1.35,  1.36,  1.38,  1.44,  1.60 /)
 atom(81:90)%radius =  (/1.70,  1.75,  1.55,  1.67,   0.0,  1.32,  2.70,  2.23,  1.88,  1.80 /)
 atom(91:95)%radius =  (/1.61,  1.39,  1.31,  1.51,   1.84 /)

  return
 end subroutine def_atomic_features


end module atomic_data

!!------------------------------------------------------------------------------------------------------------------------------!!
!     Last change:  TR   17 Oct 2006    3:34 pm

module shannon_module

 
!# Extrait de :
!#        R.D. Shannon
!#        Acta Cryst 1976, A32, 751
!#        Revised effectve ionic radii and systematic studies of interatomic distances
!#        in halides and chalogenides
!#
!# ec:     electronic configuration
!# CN:     coordinence
!# SP:     configuration de spin
!# CR:     crystal radius
!# IR:     effective radius
!#
!# ion   ec      CN      SP      CR      IR
!#

 implicit none

  private

  PUBLIC ::  set_SHANNON


  ! definitions
  TYPE, PUBLIC :: shannon_type
   CHARACTER (LEN=4)               :: ION        ! ion
   CHARACTER (LEN=4)               :: EC         ! electronic configuration
   CHARACTER (LEN=5)               :: CN         ! coordinence
   CHARACTER (LEN=2)               :: SP         ! spin
   REAL                            :: CR         ! crystal radius
   REAL                            :: IR         ! effective ionic radius
  END TYPE shannon_type
  TYPE (SHANNON_type), ALLOCATABLE, DIMENSION(:), PUBLIC  :: SHANNON

  INTEGER, PARAMETER, PUBLIC :: nb_shannon_lines = 495


contains
 subroutine set_shannon
  
  if (.not. ALLOCATED(SHANNON)) ALLOCATE(SHANNON(nb_shannon_lines))                

   SHANNON(  1) = SHANNON_type('AC+3',  '6p6 ',   'VI   ', '  ',     1.23  ,  1.12   )
   SHANNON(  2) = SHANNON_type('AG+1',  '4d10',   'II   ', '  ',     0.81  ,  0.67   )
   SHANNON(  3) = SHANNON_type('AG+1',  '4d10',   'IV   ', '  ',     1.14  ,  1.00   )
   SHANNON(  4) = SHANNON_type('AG+1',  '4d10',   'IVsq ', '  ',     1.16  ,  1.02   )
   SHANNON(  5) = SHANNON_type('AG+1',  '4d10',   'V    ', '  ',     1.23  ,  1.09   )
   SHANNON(  6) = SHANNON_type('AG+1',  '4d10',   'VI   ', '  ',     1.29  ,  1.15   )
   SHANNON(  7) = SHANNON_type('AG+1',  '4d10',   'VIII ', '  ',     1.42  ,  1.28   )
   SHANNON(  8) = SHANNON_type('AG+2',  '4d9 ',   'IVsq ', '  ',     0.93  ,  0.79   )
   SHANNON(  9) = SHANNON_type('AG+2',  '4d9 ',   'VI   ', '  ',     1.08  ,  0.94   )
   SHANNON( 10) = SHANNON_type('AG+3',  '4d8 ',   'IVsq ', '  ',     0.89  ,  0.75   )
   SHANNON( 11) = SHANNON_type('AL+3',  '2p6 ',   'IV   ', '  ',     0.53  ,  0.39   )
   SHANNON( 12) = SHANNON_type('AL+3',  '2p6 ',   'V    ', '  ',     0.62  ,  0.48   )
   SHANNON( 13) = SHANNON_type('AL+3',  '2p6 ',   'VI   ', '  ',     0.675 ,  0.535  )
   SHANNON( 14) = SHANNON_type('AM+2',  '5f7 ',   'VII  ', '  ',     1.35  ,  1.21   )
   SHANNON( 15) = SHANNON_type('AM+2',  '5f7 ',   'VIII ', '  ',     1.40  ,  1.26   )
   SHANNON( 16) = SHANNON_type('AM+2',  '5f7 ',   'IX   ', '  ',     1.45  ,  1.31   )
   SHANNON( 17) = SHANNON_type('AM+3',  '5f6 ',   'VI   ', '  ',     1.115 ,  0.975  )
   SHANNON( 18) = SHANNON_type('AM+3',  '5f6 ',   'VIII ', '  ',     1.23  ,  1.09   )
   SHANNON( 19) = SHANNON_type('AM+4',  '5f5 ',   'VI   ', '  ',     0.99  ,  0.85   )
   SHANNON( 20) = SHANNON_type('AM+4',  '5f5 ',   'VIII ', '  ',     1.09  ,  0.95   )
   SHANNON( 21) = SHANNON_type('AS+3',  '4s2 ',   'VI   ', '  ',     0.72  ,  0.58   )
   SHANNON( 22) = SHANNON_type('AS+5',  '3d10',   'IV   ', '  ',     0.475 ,  0.335  )
   SHANNON( 23) = SHANNON_type('AS+5',  '3d10',   'VI   ', '  ',     0.60  ,  0.46   )
   SHANNON( 24) = SHANNON_type('AT+7',  '5d10',   'VI   ', '  ',     0.76  ,  0.62   )
   SHANNON( 25) = SHANNON_type('AU+1',  '5d10',   'VI   ', '  ',     1.51  ,  1.37   )
   SHANNON( 26) = SHANNON_type('AU+3',  '5d8 ',   'IVsq ', '  ',     0.82  ,  0.68   )
   SHANNON( 27) = SHANNON_type('AU+3',  '5d8 ',   'VI   ', '  ',     0.99  ,  0.85   )
   SHANNON( 28) = SHANNON_type('AU+5',  '5d6 ',   'VI   ', '  ',     0.71  ,  0.57   )
   SHANNON( 29) = SHANNON_type('B+3 ',  '1s2 ',   'III  ', '  ',     0.15  ,  0.01   )
   SHANNON( 30) = SHANNON_type('B+3 ',  '1s2 ',   'IV   ', '  ',     0.25  ,  0.11   )
   SHANNON( 31) = SHANNON_type('B+3 ',  '1s2 ',   'VI   ', '  ',     0.41  ,  0.27   )
   SHANNON( 32) = SHANNON_type('BA+2',  '5p6 ',   'VI   ', '  ',     1.49  ,  1.35   )
   SHANNON( 33) = SHANNON_type('BA+2',  '5p6 ',   'VII  ', '  ',     1.52  ,  1.38   )
   SHANNON( 34) = SHANNON_type('BA+2',  '5p6 ',   'VIII ', '  ',     1.56  ,  1.42   )
   SHANNON( 35) = SHANNON_type('BA+2',  '5p6 ',   'IX   ', '  ',     1.61  ,  1.47   )
   SHANNON( 36) = SHANNON_type('BA+2',  '5p6 ',   'X    ', '  ',     1.66  ,  1.52   )
   SHANNON( 37) = SHANNON_type('BA+2',  '5p6 ',   'XI   ', '  ',     1.71  ,  1.57   )
   SHANNON( 38) = SHANNON_type('BA+2',  '5p6 ',   'XII  ', '  ',     1.75  ,  1.61   )
   SHANNON( 39) = SHANNON_type('BE+2',  '1s2 ',   'III  ', '  ',     0.30  ,  0.16   )
   SHANNON( 40) = SHANNON_type('BE+2',  '1s2 ',   'IV   ', '  ',     0.41  ,  0.27   )
   SHANNON( 41) = SHANNON_type('BE+2',  '1s2 ',   'VI   ', '  ',     0.59  ,  0.45   )
   SHANNON( 42) = SHANNON_type('BI+3',  '6s  ',   'V    ', '  ',     1.10  ,  0.96   )
   SHANNON( 43) = SHANNON_type('BI+3',  '6s  ',   'VI   ', '  ',     1.17  ,  1.03   )
   SHANNON( 44) = SHANNON_type('BI+3',  '6s  ',   'VIII ', '  ',     1.31  ,  1.17   )
   SHANNON( 45) = SHANNON_type('BI+5',  '5d10',   'VI   ', '  ',     0.90  ,  0.76   )
   SHANNON( 46) = SHANNON_type('BK+3',  '5f8 ',   'VI   ', '  ',     1.10  ,  0.96   )
   SHANNON( 47) = SHANNON_type('BK+4',  '5f7 ',   'VI   ', '  ',     0.97  ,  0.83   )
   SHANNON( 48) = SHANNON_type('BK+4',  '5f7 ',   'VIII ', '  ',     1.07  ,  0.93   )
   SHANNON( 49) = SHANNON_type('BR-1',  '4p6 ',   'VI   ', '  ',     1.82  ,  1.96   )
   SHANNON( 50) = SHANNON_type('BR+3',  '4p2 ',   'IVsq ', '  ',     0.73  ,  0.59   )
   SHANNON( 51) = SHANNON_type('BR+5',  '4s2 ',   'IIIpy', '  ',     0.45  ,  0.31   )
   SHANNON( 52) = SHANNON_type('BR+7',  '3d10',   'IV   ', '  ',     0.39  ,  0.25   )
   SHANNON( 53) = SHANNON_type('BR+7',  '3d10',   'VI   ', '  ',     0.53  ,  0.39   )
   SHANNON( 54) = SHANNON_type('C+4 ',  '1s2 ',   'III  ', '  ',     0.06  , -0.08   )
   SHANNON( 55) = SHANNON_type('C+4 ',  '1s2 ',   'IV   ', '  ',     0.29  ,  0.15   )
   SHANNON( 56) = SHANNON_type('C+4 ',  '1s2 ',   'VI   ', '  ',     0.30  ,  0.16   )
   SHANNON( 57) = SHANNON_type('CA+2',  '3p6 ',   'VI   ', '  ',     1.14  ,  1.00   )
   SHANNON( 58) = SHANNON_type('CA+2',  '3p6 ',   'VII  ', '  ',     1.20  ,  1.06   )
   SHANNON( 59) = SHANNON_type('CA+2',  '3p6 ',   'VIII ', '  ',     1.26  ,  1.12   )
   SHANNON( 60) = SHANNON_type('CA+2',  '3p6 ',   'IX   ', '  ',     1.32  ,  1.18   )
   SHANNON( 61) = SHANNON_type('CA+2',  '3p6 ',   'X    ', '  ',     1.37  ,  1.23   )
   SHANNON( 62) = SHANNON_type('CA+2',  '3p6 ',   'XII  ', '  ',     1.47  ,  1.34   )
   SHANNON( 63) = SHANNON_type('CD+2',  '4d10',   'IV   ', '  ',     0.92  ,  0.78   )
   SHANNON( 64) = SHANNON_type('CD+2',  '4d10',   'V    ', '  ',     1.01  ,  0.87   )
   SHANNON( 65) = SHANNON_type('CD+2',  '4d10',   'VI   ', '  ',     1.09  ,  0.95   )
   SHANNON( 66) = SHANNON_type('CD+2',  '4d10',   'VII  ', '  ',     1.17  ,  1.03   )
   SHANNON( 67) = SHANNON_type('CD+2',  '4d10',   'VIII ', '  ',     1.24  ,  1.10   )
   SHANNON( 68) = SHANNON_type('CD+2',  '4d10',   'XII  ', '  ',     1.45  ,  1.31   )
   SHANNON( 69) = SHANNON_type('CE+3',  '6s1 ',   'VI   ', '  ',     1.15  ,  1.01   )
   SHANNON( 70) = SHANNON_type('CE+3',  '6s1 ',   'VII  ', '  ',     1.21  ,  1.07   )
   SHANNON( 71) = SHANNON_type('CE+3',  '6s1 ',   'VIII ', '  ',     1.283 ,  1.143  )
   SHANNON( 72) = SHANNON_type('CE+3',  '6s1 ',   'IX   ', '  ',     1.336 ,  1.196  )
   SHANNON( 73) = SHANNON_type('CE+3',  '6s1 ',   'X    ', '  ',     1.39  ,  1.25   )
   SHANNON( 74) = SHANNON_type('CE+3',  '6s1 ',   'XII  ', '  ',     1.48  ,  1.34   )
   SHANNON( 75) = SHANNON_type('CE+4',  '5p6 ',   'VI   ', '  ',     1.01  ,  0.87   )
   SHANNON( 76) = SHANNON_type('CE+4',  '5p6 ',   'VIII ', '  ',     1.11  ,  0.97   )
   SHANNON( 77) = SHANNON_type('CE+4',  '5p6 ',   'X    ', '  ',     1.21  ,  1.07   )
   SHANNON( 78) = SHANNON_type('CE+4',  '5p6 ',   'XII  ', '  ',     1.28  ,  1.14   )
   SHANNON( 79) = SHANNON_type('CF+3',  '6d1 ',   'VI   ', '  ',     1.09  ,  0.95   )
   SHANNON( 80) = SHANNON_type('CF+4',  '5f8 ',   'VI   ', '  ',     0.961 ,  0.821  )
   SHANNON( 81) = SHANNON_type('CF+4',  '5f8 ',   'VIII ', '  ',     1.06  ,  0.92   )
   SHANNON( 82) = SHANNON_type('CL-1',  '3p6 ',   'VI   ', '  ',     1.67  ,  1.81   )
   SHANNON( 83) = SHANNON_type('CL+5',  '3s2 ',   'IIIpy', '  ',     0.26  ,  0.12   )
   SHANNON( 84) = SHANNON_type('CL+7',  '2p6 ',   'IV   ', '  ',     0.22  ,  0.08   )
   SHANNON( 85) = SHANNON_type('CL+7',  '2p6 ',   'VI   ', '  ',     0.41  ,  0.27   )
   SHANNON( 86) = SHANNON_type('CM+3',  '5f7 ',   'VI   ', '  ',     1.11  ,  0.97   )
   SHANNON( 87) = SHANNON_type('CM+4',  '5f6 ',   'VI   ', '  ',     0.99  ,  0.85   )
   SHANNON( 88) = SHANNON_type('CM+4',  '5f6 ',   'VIII ', '  ',     1.09  ,  0.95   )
   SHANNON( 89) = SHANNON_type('CO+2',  '3d7 ',   'IV   ', 'HS',     0.72  ,  0.58   )
   SHANNON( 90) = SHANNON_type('CO+2',  '3d7 ',   'V    ', '  ',     0.81  ,  0.67   )
   SHANNON( 91) = SHANNON_type('CO+2',  '3d7 ',   'VI   ', 'LS',     0.79  ,  0.65   )
   SHANNON( 92) = SHANNON_type('CO+2',  '3d7 ',   'VI   ', 'HS',     0.885 ,  0.745  )
   SHANNON( 93) = SHANNON_type('CO+2',  '3d7 ',   'VIII ', '  ',     1.04  ,  0.90   )
   SHANNON( 94) = SHANNON_type('CO+3',  '3d6 ',   'VI   ', 'LS',     0.685 ,  0.545  )
   SHANNON( 95) = SHANNON_type('CO+3',  '3d6 ',   'VI   ', 'HS',     0.75  ,  0.61   )
   SHANNON( 96) = SHANNON_type('CO+4',  '3d5 ',   'IV   ', '  ',     0.54  ,  0.40   )
   SHANNON( 97) = SHANNON_type('CO+4',  '3d5 ',   'VI   ', 'HS',     0.67  ,  0.53   )
   SHANNON( 98) = SHANNON_type('CR+2',  '3d4 ',   'VI   ', 'LS',     0.87  ,  0.73   )
   SHANNON( 99) = SHANNON_type('CR+2',  '3d4 ',   'VI   ', 'HS',     0.94  ,  0.80   )
   SHANNON(100) = SHANNON_type('CR+3',  '3d3 ',   'VI   ', '  ',     0.755 ,  0.615  )
   SHANNON(101) = SHANNON_type('CR+4',  '3d2 ',   'IV   ', '  ',     0.55  ,  0.41   )
   SHANNON(102) = SHANNON_type('CR+4',  '3d2 ',   'VI   ', '  ',     0.69  ,  0.55   )
   SHANNON(103) = SHANNON_type('CR+5',  '3d1 ',   'IV   ', '  ',     0.485 ,  0.345  )
   SHANNON(104) = SHANNON_type('CR+5',  '3d1 ',   'VI   ', '  ',     0.63  ,  0.49   )
   SHANNON(105) = SHANNON_type('CR+5',  '3d1 ',   'VIII ', '  ',     0.71  ,  0.57   )
   SHANNON(106) = SHANNON_type('CR+6',  '3p6 ',   'IV   ', '  ',     0.40  ,  0.26   )
   SHANNON(107) = SHANNON_type('CR+6',  '3p6 ',   'VI   ', '  ',     0.58  ,  0.44   )
   SHANNON(108) = SHANNON_type('CS+1',  '5p6 ',   'VI   ', '  ',     1.81  ,  1.67   )
   SHANNON(109) = SHANNON_type('CS+1',  '5p6 ',   'VIII ', '  ',     1.88  ,  1.74   )
   SHANNON(110) = SHANNON_type('CS+1',  '5p6 ',   'IX   ', '  ',     1.92  ,  1.78   )
   SHANNON(111) = SHANNON_type('CS+1',  '5p6 ',   'X    ', '  ',     1.95  ,  1.81   )
   SHANNON(112) = SHANNON_type('CS+1',  '5p6 ',   'XI   ', '  ',     1.99  ,  1.85   )
   SHANNON(113) = SHANNON_type('CS+1',  '5p6 ',   'XII  ', '  ',     2.02  ,  1.88   )
   SHANNON(114) = SHANNON_type('CU+1',  '3d10',   'II   ', '  ',     0.90  ,  0.46   )
   SHANNON(115) = SHANNON_type('CU+1',  '3d10',   'IV   ', '  ',     0.74  ,  0.60   )
   SHANNON(116) = SHANNON_type('CU+1',  '3d10',   'VI   ', '  ',     0.91  ,  0.77   )
   SHANNON(117) = SHANNON_type('CU+2',  '3d9 ',   'IV   ', '  ',     0.71  ,  0.57   )
   SHANNON(118) = SHANNON_type('CU+2',  '3d9 ',   'IVsq ', '  ',     0.71  ,  0.57   )
   SHANNON(119) = SHANNON_type('CU+2',  '3d9 ',   'V    ', '  ',     0.79  ,  0.65   )
   SHANNON(120) = SHANNON_type('CU+2',  '3d9 ',   'VI   ', '  ',     0.87  ,  0.73   )
   SHANNON(121) = SHANNON_type('CI+3',  '3d8 ',   'VI   ', 'LS',     0.68  ,  0.54   )
   SHANNON(122) = SHANNON_type('D+1 ',  '1s0 ',   'II   ', '  ',     0.04  , -0.10   )
   SHANNON(123) = SHANNON_type('DY+2',  '4f10',   'VI   ', '  ',     0.21  ,  1.07   )
   SHANNON(124) = SHANNON_type('DY+2',  '4f10',   'VII  ', '  ',     1.27  ,  1.13   )
   SHANNON(125) = SHANNON_type('DY+2',  '4f10',   'VIII ', '  ',     1.33  ,  1.19   )
   SHANNON(126) = SHANNON_type('DY+3',  '4f9 ',   'VI   ', '  ',     1.052 ,  0.912  )
   SHANNON(127) = SHANNON_type('DY+3',  '4f9 ',   'VII  ', '  ',     1.11  ,  0.97   )
   SHANNON(128) = SHANNON_type('DY+3',  '4f9 ',   'VIII ', '  ',     1.167 ,  1.027  )
   SHANNON(129) = SHANNON_type('DY+3',  '4f9 ',   'IX   ', '  ',     1.223 ,  1.083  )
   SHANNON(130) = SHANNON_type('ER+3',  '4f11',   'VI   ', '  ',     1.030 ,  0.890  )
   SHANNON(131) = SHANNON_type('ER+3',  '4f11',   'VII  ', '  ',     1.085 ,  0.945  )
   SHANNON(132) = SHANNON_type('ER+3',  '4f11',   'VIII ', '  ',     1.144 ,  1.004  )
   SHANNON(133) = SHANNON_type('ER+3',  '4f11',   'IX   ', '  ',     1.202 ,  1.062  )
   SHANNON(134) = SHANNON_type('EU+2',  '4f7 ',   'VI   ', '  ',     1.31  ,  1.17   )
   SHANNON(135) = SHANNON_type('EU+2',  '4f7 ',   'VII  ', '  ',     1.34  ,  1.20   )
   SHANNON(136) = SHANNON_type('EU+2',  '4f7 ',   'VIII ', '  ',     1.396 ,  1.25   )
   SHANNON(137) = SHANNON_type('EU+2',  '4f7 ',   'IX   ', '  ',     1.44  ,  1.30   )
   SHANNON(138) = SHANNON_type('EU+2',  '4f7 ',   'X    ', '  ',     1.49  ,  1.35   )
   SHANNON(139) = SHANNON_type('EU+3',  '4f6 ',   'VI   ', '  ',     1.087 ,  0.947  )
   SHANNON(140) = SHANNON_type('EU+3',  '4f6 ',   'VII  ', '  ',     1.15  ,  1.01   )
   SHANNON(141) = SHANNON_type('EU+3',  '4f6 ',   'VIII ', '  ',     1.206 ,  1.066  )
   SHANNON(142) = SHANNON_type('EU+3',  '4f6 ',   'IX   ', '  ',     1.260 ,  1.120  )
   SHANNON(143) = SHANNON_type('F-1 ',  '2p6 ',   'II   ', '  ',     1.145 ,  1.285  )
   SHANNON(144) = SHANNON_type('F-1 ',  '2p6 ',   'III  ', '  ',     1.16  ,  1.30   )
   SHANNON(145) = SHANNON_type('F-1 ',  '2p6 ',   'IV   ', '  ',     1.17  ,  1.31   )
   SHANNON(146) = SHANNON_type('F-1 ',  '2p6 ',   'VI   ', '  ',     1.19  ,  1.33   )
   SHANNON(147) = SHANNON_type('F+7 ',  '1s2 ',   'VI   ', '  ',     0.22  ,  0.08   )
   SHANNON(148) = SHANNON_type('FE+2',  '3d6 ',   'IV   ', 'HS',     0.77  ,  0.63   )
   SHANNON(149) = SHANNON_type('FE+2',  '3d6 ',   'IVsq ', 'HS',     0.78  ,  0.64   )
   SHANNON(150) = SHANNON_type('FE+2',  '3d6 ',   'VI   ', 'LS',     0.75  ,  0.61   )
   SHANNON(151) = SHANNON_type('FE+2',  '3d6 ',   'VI   ', 'HS',     0.920 ,  0.780  )
   SHANNON(152) = SHANNON_type('FE+2',  '3d6 ',   'VIII ', 'HS',     1.06  ,  0.92   )
   SHANNON(153) = SHANNON_type('FE+3',  '3d5 ',   'IV   ', 'HS',     0.63  ,  0.49   )
   SHANNON(154) = SHANNON_type('FE+3',  '3d5 ',   'V    ', '  ',     0.72  ,  0.58   )
   SHANNON(155) = SHANNON_type('FE+3',  '3d5 ',   'VI   ', 'LS',     0.69  ,  0.55   )
   SHANNON(156) = SHANNON_type('FE+3',  '3d5 ',   'VI   ', 'HS',     0.785 ,  0.645  )
   SHANNON(157) = SHANNON_type('FE+3',  '3d5 ',   'VIII ', 'HS',     0.92  ,  0.78   )
   SHANNON(158) = SHANNON_type('FE+4',  '3d4 ',   'VI   ', '  ',     0.725 ,  0.585  )
   SHANNON(159) = SHANNON_type('FE+6',  '3d2 ',   'IV   ', '  ',     0.39  ,  0.25   )
   SHANNON(160) = SHANNON_type('FR+1',  '6p6 ',   'VI   ', '  ',     1.94  ,  1.80   )
   SHANNON(161) = SHANNON_type('GA+3',  '3d10',   'IV   ', '  ',     0.61  ,  0.47   )
   SHANNON(162) = SHANNON_type('GA+3',  '3d10',   'V    ', '  ',     0.69  ,  0.55   )
   SHANNON(163) = SHANNON_type('GA+3',  '3d10',   'VI   ', '  ',     0.760 ,  0.620  )
   SHANNON(164) = SHANNON_type('GD+3',  '4f7 ',   'VI   ', '  ',     1.078 ,  0.938  )
   SHANNON(165) = SHANNON_type('GD+3',  '4f7 ',   'VII  ', '  ',     1.14  ,  1.00   )
   SHANNON(166) = SHANNON_type('GD+3',  '4f7 ',   'VIII ', '  ',     1.193 ,  1.053  )
   SHANNON(167) = SHANNON_type('GD+3',  '4f7 ',   'IX   ', '  ',     1.247 ,  1.107  )
   SHANNON(168) = SHANNON_type('GE+2',  '4s2 ',   'VI   ', '  ',     0.87  ,  0.73   )
   SHANNON(169) = SHANNON_type('GE+4',  '3d10',   'IV   ', '  ',     0.530 ,  0.390  )
   SHANNON(170) = SHANNON_type('GE+4',  '3d10',   'VI   ', '  ',     0.670 ,  0.530  )
   SHANNON(171) = SHANNON_type('H+1 ',  '1s0 ',   'I    ', '  ',    -0.24  , -0.38   )
   SHANNON(172) = SHANNON_type('H+1 ',  '1s0 ',   'II   ', '  ',    -0.04  ,  0.18   )
   SHANNON(173) = SHANNON_type('HF+4',  '4f14',   'IV   ', '  ',     0.72  ,  0.58   )
   SHANNON(174) = SHANNON_type('HF+4',  '4f14',   'VI   ', '  ',     0.85  ,  0.71   )
   SHANNON(175) = SHANNON_type('HF+4',  '4f14',   'VII  ', '  ',     0.90  ,  0.76   )
   SHANNON(176) = SHANNON_type('HF+4',  '4f14',   'VIII ', '  ',     0.97  ,  0.83   )
   SHANNON(177) = SHANNON_type('HG+1',  '6s1 ',   'III  ', '  ',     1.11  ,  0.97   )
   SHANNON(178) = SHANNON_type('HG+1',  '6s1 ',   'VI   ', '  ',     1.33  ,  1.19   )
   SHANNON(179) = SHANNON_type('HG+2',  '5d10',   'II   ', '  ',     0.83  ,  0.69   )
   SHANNON(180) = SHANNON_type('HG+2',  '5d10',   'IV   ', '  ',     1.10  ,  0.96   )
   SHANNON(181) = SHANNON_type('HG+2',  '5d10',   'VI   ', '  ',     1.16  ,  1.02   )
   SHANNON(182) = SHANNON_type('HG+2',  '5d10',   'VIII ', '  ',     1.28  ,  1.14   )
   SHANNON(183) = SHANNON_type('HO+3',  '4f10',   'VI   ', '  ',     1.041 ,  0.901  )
   SHANNON(184) = SHANNON_type('HO+3',  '4f10',   'VIII ', '  ',     1.155 ,  1.015  )
   SHANNON(185) = SHANNON_type('HO+3',  '4f10',   'IX   ', '  ',     1.212 ,  1.072  )
   SHANNON(186) = SHANNON_type('HO+3',  '4f10',   'X    ', '  ',     1.26  ,  1.12   )
   SHANNON(187) = SHANNON_type('I-1 ',  '5p6 ',   'VI   ', '  ',     2.06  ,  2.20   )
   SHANNON(188) = SHANNON_type('I+5 ',  '5s2 ',   'IIIpy', '  ',     0.58  ,  0.44   )
   SHANNON(189) = SHANNON_type('I+5 ',  '5s2 ',   'I    ', '  ',     1.09  ,  0.95   )
   SHANNON(190) = SHANNON_type('I+7 ',  '4d10',   'IV   ', '  ',     0.56  ,  0.42   )
   SHANNON(191) = SHANNON_type('I+7 ',  '4d10',   'VI   ', '  ',     0.67  ,  0.53   )
   SHANNON(192) = SHANNON_type('IN+3',  '4d10',   'IV   ', '  ',     0.76  ,  0.62   )
   SHANNON(193) = SHANNON_type('IN+3',  '4d10',   'VI   ', '  ',     0.910 ,  0.800  )
   SHANNON(194) = SHANNON_type('IN+3',  '4d10',   'VIII ', '  ',     1.06  ,  0.92   )
   SHANNON(195) = SHANNON_type('IR+3',  '5d6 ',   'VI   ', '  ',     0.82  ,  0.68   )
   SHANNON(196) = SHANNON_type('IR+4',  '5d5 ',   'VI   ', '  ',     0.765 ,  0.625  )
   SHANNON(197) = SHANNON_type('IR+5',  '5d4 ',   'VI   ', '  ',     0.71  ,  0.57   )
   SHANNON(198) = SHANNON_type('K+1 ',  '3p6 ',   'IV   ', '  ',     1.51  ,  1.37   )
   SHANNON(199) = SHANNON_type('K+1 ',  '3p6 ',   'VI   ', '  ',     1.52  ,  1.38   )
   SHANNON(200) = SHANNON_type('K+1 ',  '3p6 ',   'VII  ', '  ',     1.60  ,  1.46   )
   SHANNON(201) = SHANNON_type('K+1 ',  '3p6 ',   'VIII ', '  ',     1.65  ,  1.51   )
   SHANNON(202) = SHANNON_type('K+1 ',  '3p6 ',   'IX   ', '  ',     1.69  ,  1.55   )
   SHANNON(203) = SHANNON_type('K+1 ',  '3p6 ',   'X    ', '  ',     1.73  ,  1.59   )
   SHANNON(204) = SHANNON_type('K+1 ',  '3p6 ',   'XII  ', '  ',     1.78  ,  1.64   )
   SHANNON(205) = SHANNON_type('LA+3',  '4d10',   'VI   ', '  ',     1.172 ,  1.032  )
   SHANNON(206) = SHANNON_type('LA+3',  '4d10',   'VII  ', '  ',     1.24  ,  1.10   )
   SHANNON(207) = SHANNON_type('LA+3',  '4d10',   'VIII ', '  ',     1.300 ,  1.160  )
   SHANNON(208) = SHANNON_type('LA+3',  '4d10',   'IX   ', '  ',     1.356 ,  1.216  )
   SHANNON(209) = SHANNON_type('LA+3',  '4d10',   'X    ', '  ',     1.41  ,  1.27   )
   SHANNON(210) = SHANNON_type('LA+3',  '4d10',   'XII  ', '  ',     1.50  ,  1.36   )
   SHANNON(211) = SHANNON_type('LI+1',  '1s2 ',   'IV   ', '  ',     0.730 ,  0.590  )
   SHANNON(212) = SHANNON_type('LI+1',  '1s2 ',   'VI   ', '  ',     0.90  ,  0.76   )
   SHANNON(213) = SHANNON_type('LI+1',  '1s2 ',   'VIII ', '  ',     1.06  ,  0.92   )
   SHANNON(214) = SHANNON_type('LU+3',  '4f14',   'VI   ', '  ',     1.001 ,  0.861  )
   SHANNON(215) = SHANNON_type('LU+3',  '4f14',   'VIII ', '  ',     1.117 ,  0.977  )
   SHANNON(216) = SHANNON_type('LU+3',  '4f14',   'IX   ', '  ',     1.172 ,  1.032  )
   SHANNON(217) = SHANNON_type('MG+2',  '2p6 ',   'IV   ', '  ',     0.71  ,  0.57   )
   SHANNON(218) = SHANNON_type('MG+2',  '2p6 ',   'V    ', '  ',     0.80  ,  0.66   )
   SHANNON(219) = SHANNON_type('MG+2',  '2p6 ',   'VI   ', '  ',     0.860 ,  0.720  )
   SHANNON(220) = SHANNON_type('MG+2',  '2p6 ',   'VIII ', '  ',     1.03  ,  0.89   )
   SHANNON(221) = SHANNON_type('MN+2',  '3d5 ',   'IV   ', 'HS',     0.80  ,  0.66   )
   SHANNON(222) = SHANNON_type('MN+2',  '3d5 ',   'V    ', 'HS',     0.89  ,  0.75   )
   SHANNON(223) = SHANNON_type('MN+2',  '3d5 ',   'VI   ', 'LS',     0.81  ,  0.67   )
   SHANNON(224) = SHANNON_type('MN+2',  '3d5 ',   'VI   ', 'HS',     0.970 ,  0.830  )
   SHANNON(225) = SHANNON_type('MN+2',  '3d5 ',   'VII  ', 'HS',     1.04  ,  0.90   )
   SHANNON(226) = SHANNON_type('MN+2',  '3d5 ',   'VIII ', '  ',     1.10  ,  0.96   )
   SHANNON(227) = SHANNON_type('MN+3',  '3d4 ',   'V    ', '  ',     0.72  ,  0.58   )
   SHANNON(228) = SHANNON_type('MN+3',  '3d4 ',   'VI   ', 'LS',     0.72  ,  0.58   )
   SHANNON(229) = SHANNON_type('MN+3',  '3d4 ',   'VI   ', 'HS',     0.785 ,  0.645  )
   SHANNON(230) = SHANNON_type('MN+4',  '3d3 ',   'IV   ', '  ',     0.53  ,  0.39   )
   SHANNON(231) = SHANNON_type('MN+4',  '3d3 ',   'VI   ', '  ',     0.670 ,  0.530  )
   SHANNON(232) = SHANNON_type('MN+5',  '3d2 ',   'IV   ', '  ',     0.47  ,  0.33   )
   SHANNON(233) = SHANNON_type('MN+6',  '3d1 ',   'IV   ', '  ',     0.395 ,  0.255  )
   SHANNON(234) = SHANNON_type('MN+7',  '3p6 ',   'IV   ', '  ',     0.39  ,  0.25   )
   SHANNON(235) = SHANNON_type('MN+7',  '3p6 ',   'VI   ', '  ',     0.60  ,  0.46   )
   SHANNON(236) = SHANNON_type('MO+3',  '4d3 ',   'VI   ', '  ',     0.83  ,  0.69   )
   SHANNON(237) = SHANNON_type('MO+4',  '4d2 ',   'VI   ', '  ',     0.790 ,  0.650  )
   SHANNON(238) = SHANNON_type('MO+5',  '4d1 ',   'IV   ', '  ',     0.60  ,  0.46   )
   SHANNON(239) = SHANNON_type('MO+5',  '4d1 ',   'VI   ', '  ',     0.75  ,  0.61   )
   SHANNON(240) = SHANNON_type('MO+6',  '4p6 ',   'IV   ', '  ',     0.55  ,  0.41   )
   SHANNON(241) = SHANNON_type('MO+6',  '4p6 ',   'V    ', '  ',     0.64  ,  0.50   )
   SHANNON(242) = SHANNON_type('MO+6',  '4p6 ',   'VI   ', '  ',     0.73  ,  0.59   )
   SHANNON(243) = SHANNON_type('MO+6',  '4p6 ',   'VII  ', '  ',     0.87  ,  0.73   )
   SHANNON(244) = SHANNON_type('N-3 ',  '2p6 ',   'IV   ', '  ',     1.32  ,  1.46   )
   SHANNON(245) = SHANNON_type('N+3 ',  '2s2 ',   'VI   ', '  ',     0.30  ,  0.16   )
   SHANNON(246) = SHANNON_type('N+5 ',  '1s2 ',   'III  ', '  ',     0.044 , -0.104  )
   SHANNON(247) = SHANNON_type('N+5 ',  '1s2 ',   'VI   ', '  ',     0.27  ,  0.13   )
   SHANNON(248) = SHANNON_type('NA+1',  '2p6 ',   'IV   ', '  ',     1.13  ,  0.99   )
   SHANNON(249) = SHANNON_type('NA+1',  '2p6 ',   'V    ', '  ',     1.14  ,  1.00   )
   SHANNON(250) = SHANNON_type('NA+1',  '2p6 ',   'VI   ', '  ',     1.16  ,  1.02   )
   SHANNON(251) = SHANNON_type('NA+1',  '2p6 ',   'VII  ', '  ',     1.26  ,  1.12   )
   SHANNON(252) = SHANNON_type('NA+1',  '2p6 ',   'VIII ', '  ',     1.32  ,  1.18   )
   SHANNON(253) = SHANNON_type('NA+1',  '2p6 ',   'IX   ', '  ',     1.38  ,  1.24   )
   SHANNON(254) = SHANNON_type('NA+1',  '2p6 ',   'XII  ', '  ',     1.053 ,  1.39   )
   SHANNON(255) = SHANNON_type('NB+3',  '4d2 ',   'VI   ', '  ',     0.86  ,  0.72   )
   SHANNON(256) = SHANNON_type('NB+4',  '4d1 ',   'VI   ', '  ',     0.82  ,  0.68   )
   SHANNON(257) = SHANNON_type('NB+4',  '4d1 ',   'VIII ', '  ',     0.93  ,  0.79   )
   SHANNON(258) = SHANNON_type('NB+5',  '4p6 ',   'IV   ', '  ',     0.62  ,  0.48   )
   SHANNON(259) = SHANNON_type('NB+5',  '4p6 ',   'VI   ', '  ',     0.78  ,  0.64   )
   SHANNON(260) = SHANNON_type('NB+5',  '4p6 ',   'VII  ', '  ',     0.83  ,  0.69   )
   SHANNON(261) = SHANNON_type('NB+5',  '4p6 ',   'VIII ', '  ',     0.88  ,  0.74   )
   SHANNON(262) = SHANNON_type('ND+2',  '4f4 ',   'VIII ', '  ',     1.43  ,  1.29   )
   SHANNON(263) = SHANNON_type('ND+2',  '4f4 ',   'IX   ', '  ',     1.49  ,  1.35   )
   SHANNON(264) = SHANNON_type('ND+3',  '4f3 ',   'VI   ', '  ',     1.123 ,  0.983  )
   SHANNON(265) = SHANNON_type('ND+3',  '4f3 ',   'VIII ', '  ',     1.249 ,  1.109  )
   SHANNON(266) = SHANNON_type('ND+3',  '4f3 ',   'IX   ', '  ',     1.303 ,  1.163  )
   SHANNON(267) = SHANNON_type('ND+3',  '4f3 ',   'XII  ', '  ',     1.41  ,  1.27   )
   SHANNON(268) = SHANNON_type('NI+2',  '3d8 ',   'IV   ', '  ',     0.69  ,  0.55   )
   SHANNON(269) = SHANNON_type('NI+2',  '3d8 ',   'IVsq ', '  ',     0.63  ,  0.49   )
   SHANNON(270) = SHANNON_type('NI+2',  '3d8 ',   'V    ', '  ',     0.77  ,  0.63   )
   SHANNON(271) = SHANNON_type('NI+2',  '3d8 ',   'VI   ', '  ',     0.830 ,  0.690  )
   SHANNON(272) = SHANNON_type('NI+3',  '3d7 ',   'VI   ', 'LS',     0.70  ,  0.56   )
   SHANNON(273) = SHANNON_type('NI+3',  '3d7 ',   'VI   ', 'HS',     0.74  ,  0.60   )
   SHANNON(274) = SHANNON_type('NI+4',  '3d6 ',   'VI   ', 'LS',     0.62  ,  0.48   )
   SHANNON(275) = SHANNON_type('NO+2',  '5f14',   'VI   ', '  ',     1.24  ,  1.1    )
   SHANNON(276) = SHANNON_type('NP+2',  '5f5 ',   'VI   ', '  ',     1.24  ,  1.10   )
   SHANNON(277) = SHANNON_type('NP+3',  '5f4 ',   'VI   ', '  ',     1.15  ,  1.01   )
   SHANNON(278) = SHANNON_type('NP+4',  '5f3 ',   'VI   ', '  ',     1.01  ,  0.87   )
   SHANNON(279) = SHANNON_type('NP+4',  '5f3 ',   'VIII ', '  ',     1.12  ,  0.98   )
   SHANNON(280) = SHANNON_type('NP+5',  '5f2 ',   'VI   ', '  ',     0.89  ,  0.75   )
   SHANNON(281) = SHANNON_type('NP+6',  '5f1 ',   'VI   ', '  ',     0.86  ,  0.72   )
   SHANNON(282) = SHANNON_type('NP+7',  '6p6 ',   'VI   ', '  ',     0.85  ,  0.71   )
   SHANNON(283) = SHANNON_type('O-2 ',  '2p6 ',   'II   ', '  ',     1.21  ,  1.35   )
   SHANNON(284) = SHANNON_type('O-2 ',  '2p6 ',   'III  ', '  ',     1.22  ,  1.36   )
   SHANNON(285) = SHANNON_type('O-2 ',  '2p6 ',   'IV   ', '  ',     1.24  ,  1.38   )
   SHANNON(286) = SHANNON_type('O-2 ',  '2p6 ',   'VI   ', '  ',     1.26  ,  1.40   )
   SHANNON(287) = SHANNON_type('O-2 ',  '2p6 ',   'VIII ', '  ',     1.28  ,  1.32   )
   SHANNON(288) = SHANNON_type('OH-1',  '    ',   'II   ', '  ',     1.18  ,  1.32   )
   SHANNON(289) = SHANNON_type('OH-1',  '    ',   'III  ', '  ',     1.20  ,  1.34   )
   SHANNON(290) = SHANNON_type('OH-1',  '    ',   'IV   ', '  ',     1.21  ,  1.35   )
   SHANNON(291) = SHANNON_type('OH-1',  '    ',   'VI   ', '  ',     1.23  ,  1.37   )
   SHANNON(292) = SHANNON_type('OS+4',  '5d4 ',   'VI   ', '  ',     0.770 ,  0.630  )
   SHANNON(293) = SHANNON_type('OS+5',  '5d3 ',   'VI   ', '  ',     0.715 ,  0.575  )
   SHANNON(294) = SHANNON_type('OS+6',  '5d2 ',   'V    ', '  ',     0.63  ,  0.49   )
   SHANNON(295) = SHANNON_type('OS+6',  '5d2 ',   'VI   ', '  ',     0.685 ,  0.545  )
   SHANNON(296) = SHANNON_type('OS+7',  '5d1 ',   'VI   ', '  ',     0.665 ,  0.525  )
   SHANNON(297) = SHANNON_type('OS+8',  '5p6 ',   'IV   ', '  ',     0.53  ,  0.39   )
   SHANNON(298) = SHANNON_type('P+3 ',  '3s2 ',   'VI   ', '  ',     0.58  ,  0.44   )
   SHANNON(299) = SHANNON_type('P+5 ',  '2p6 ',   'IV   ', '  ',     0.31  ,  0.17   )
   SHANNON(300) = SHANNON_type('P+5 ',  '2p6 ',   'V    ', '  ',     0.43  ,  0.29   )
   SHANNON(301) = SHANNON_type('P+5 ',  '2p6 ',   'VI   ', '  ',     0.52  ,  0.38   )
   SHANNON(302) = SHANNON_type('PA+3',  '5f2 ',   'VI   ', '  ',     1.18  ,  1.04   )
   SHANNON(303) = SHANNON_type('PA+4',  '6d1 ',   'VI   ', '  ',     1.04  ,  0.90   )
   SHANNON(304) = SHANNON_type('PA+4',  '6d1 ',   'VIII ', '  ',     1.15  ,  1.01   )
   SHANNON(305) = SHANNON_type('PA+5',  '6p6 ',   'VI   ', '  ',     0.92  ,  0.78   )
   SHANNON(306) = SHANNON_type('PA+5',  '6p6 ',   'VIII ', '  ',     1.05  ,  0.91   )
   SHANNON(307) = SHANNON_type('PA+5',  '6p6 ',   'IX   ', '  ',     1.09  ,  0.95   )
   SHANNON(308) = SHANNON_type('PB+2',  '6s2 ',   'IVpy ', '  ',     1.12  ,  0.98   )
   SHANNON(309) = SHANNON_type('PB+2',  '6s2 ',   'VI   ', '  ',     1.33  ,  1.19   )
   SHANNON(310) = SHANNON_type('PB+2',  '6s2 ',   'VII  ', '  ',     1.37  ,  1.23   )
   SHANNON(311) = SHANNON_type('PB+2',  '6s2 ',   'VIII ', '  ',     1.43  ,  1.29   )
   SHANNON(312) = SHANNON_type('PB+2',  '6s2 ',   'IX   ', '  ',     1.49  ,  1.35   )
   SHANNON(313) = SHANNON_type('PB+2',  '6s2 ',   'X    ', '  ',     1.54  ,  1.40   )
   SHANNON(314) = SHANNON_type('PB+2',  '6s2 ',   'XI   ', '  ',     1.59  ,  1.45   )
   SHANNON(315) = SHANNON_type('PB+2',  '6s2 ',   'XII  ', '  ',     1.63  ,  1.49   )
   SHANNON(316) = SHANNON_type('PB+4',  '4d10',   'IV   ', '  ',     0.79  ,  0.65   )
   SHANNON(317) = SHANNON_type('PB+4',  '4d10',   'V    ', '  ',     0.87  ,  0.73   )
   SHANNON(318) = SHANNON_type('PB+4',  '4d10',   'VI   ', '  ',     0.915 ,  0.775  )
   SHANNON(319) = SHANNON_type('PB+4',  '4d10',   'VIII ', '  ',     1.08  ,  0.94   )
   SHANNON(320) = SHANNON_type('PD+1',  '4d9 ',   'II   ', '  ',     0.73  ,  0.59   )
   SHANNON(321) = SHANNON_type('PD+2',  '4d8 ',   'IVsq ', '  ',     0.78  ,  0.64   )
   SHANNON(322) = SHANNON_type('PD+2',  '4d8 ',   'VI   ', '  ',     1.00  ,  0.86   )
   SHANNON(323) = SHANNON_type('PD+3',  '4d7 ',   'VI   ', '  ',     0.90  ,  0.76   )
   SHANNON(324) = SHANNON_type('PD+4',  '4d6 ',   'VI   ', '  ',     0.755 ,  0.615  )
   SHANNON(325) = SHANNON_type('PM+3',  '4f4 ',   'VI   ', '  ',     1.11  ,  0.97   )
   SHANNON(326) = SHANNON_type('PM+3',  '4f4 ',   'VIII ', '  ',     1.233 ,  1.093  )
   SHANNON(327) = SHANNON_type('PM+3',  '4f4 ',   'IX   ', '  ',     1.284 ,  1.144  )
   SHANNON(328) = SHANNON_type('PO+4',  '6s2 ',   'VI   ', '  ',     1.08  ,  0.94   )
   SHANNON(329) = SHANNON_type('PO+4',  '6s2 ',   'VIII ', '  ',     1.22  ,  1.08   )
   SHANNON(330) = SHANNON_type('PO+6',  '5d10',   'VI   ', '  ',     0.81  ,  0.67   )
   SHANNON(331) = SHANNON_type('PR+3',  '4f2 ',   'VI   ', '  ',     1.13  ,  0.99   )
   SHANNON(332) = SHANNON_type('PR+3',  '4f2 ',   'VIII ', '  ',     1.266 ,  1.126  )
   SHANNON(333) = SHANNON_type('PR+3',  '4f2 ',   'IX   ', '  ',     1.319 ,  1.179  )
   SHANNON(334) = SHANNON_type('PR+4',  '4f1 ',   'VI   ', '  ',     0.99  ,  0.85   )
   SHANNON(335) = SHANNON_type('PR+4',  '4f1 ',   'VIII ', '  ',     1.10  ,    0.96 )
   SHANNON(336) = SHANNON_type('PT+2',  '5d8 ',   'IVsq ', '  ',     0.74  ,  0.60   )
   SHANNON(337) = SHANNON_type('PT+2',  '5d8 ',   'VI   ', '  ',     0.94  ,  0.80   )
   SHANNON(338) = SHANNON_type('PT+4',  '5d6 ',   'VI   ', '  ',     0.765 ,  0.625  )
   SHANNON(339) = SHANNON_type('PT+5',  '5d5 ',   'VI   ', '  ',     0.71  ,  0.57   )
   SHANNON(340) = SHANNON_type('PU+3',  '5f5 ',   'VI   ', '  ',     1.14  ,  1.00   )
   SHANNON(341) = SHANNON_type('PU+4',  '5f4 ',   'VI   ', '  ',     1.00  ,  0.86   )
   SHANNON(342) = SHANNON_type('PU+4',  '5f4 ',   'VIII ', '  ',     1.10  ,  0.96   )
   SHANNON(343) = SHANNON_type('PU+5',  '5f3 ',   'VI   ', '  ',     0.88  ,  0.74   )
   SHANNON(344) = SHANNON_type('PU+6',  '5f2 ',   'VI   ', '  ',     0.85  ,  0.71   )
   SHANNON(345) = SHANNON_type('RA+2',  '6p6 ',   'VIII ', '  ',     1.62  ,  1.48   )
   SHANNON(346) = SHANNON_type('RA+2',  '6p6 ',   'XII  ', '  ',     1.84  ,  1.70   )
   SHANNON(347) = SHANNON_type('RB+1',  '4p6 ',   'VI   ', '  ',     1.66  ,  1.52   )
   SHANNON(348) = SHANNON_type('RB+1',  '4p6 ',   'VII  ', '  ',     1.70  ,  1.56   )
   SHANNON(349) = SHANNON_type('RB+1',  '4p6 ',   'VIII ', '  ',     1.75  ,  1.61   )
   SHANNON(350) = SHANNON_type('RB+1',  '4p6 ',   'IX   ', '  ',     1.77  ,  1.63   )
   SHANNON(351) = SHANNON_type('RB+1',  '4p6 ',   'X    ', '  ',     1.80  ,  1.66   )
   SHANNON(352) = SHANNON_type('RB+1',  '4p6 ',   'XI   ', '  ',     1.83  ,  1.69   )
   SHANNON(353) = SHANNON_type('RB+1',  '4p6 ',   'XII  ', '  ',     1.86  ,  1.72   )
   SHANNON(354) = SHANNON_type('RB+1',  '4p6 ',   'XIV  ', '  ',     1.97  ,  1.83   )
   SHANNON(355) = SHANNON_type('RE+4',  '5d3 ',   'VI   ', '  ',     0.77  ,  0.63   )
   SHANNON(356) = SHANNON_type('RE+5',  '5d2 ',   'VI   ', '  ',     0.72  ,  0.58   )
   SHANNON(357) = SHANNON_type('RE+6',  '5d1 ',   'VI   ', '  ',     0.69  ,  0.55   )
   SHANNON(358) = SHANNON_type('RE+7',  '5p6 ',   'IV   ', '  ',     0.52  ,  0.38   )
   SHANNON(359) = SHANNON_type('RE+7',  '5p6 ',   'VI   ', '  ',     0.67  ,  0.53   )
   SHANNON(360) = SHANNON_type('RH+3',  '4d6 ',   'VI   ', '  ',     0.805 ,  0.60   )
   SHANNON(361) = SHANNON_type('RH+4',  '4d5 ',   'VI   ', '  ',     0.74  ,  0.60   )
   SHANNON(362) = SHANNON_type('RH+5',  '5d4 ',   'VI   ', '  ',     0.69  ,  0.55   )
   SHANNON(363) = SHANNON_type('RU+3',  '4d5 ',   'VI   ', '  ',     0.805 ,  0.665  )
   SHANNON(364) = SHANNON_type('RU+4',  '4d4 ',   'VI   ', '  ',     0.760 ,  0.620  )
   SHANNON(365) = SHANNON_type('RU+5',  '4d3 ',   'VI   ', '  ',     0.705 ,  0.565  )
   SHANNON(366) = SHANNON_type('RU+7',  '4d1 ',   'VI   ', '  ',     0.52  ,  0.38   )
   SHANNON(367) = SHANNON_type('RU+8',  '4p6 ',   'IV   ', '  ',     0.50  ,  0.36   )
   SHANNON(368) = SHANNON_type('S-2 ',  '3p6 ',   'VI   ', '  ',     1.70  ,  1.84   )
   SHANNON(369) = SHANNON_type('S+4 ',  '3s2 ',   'VI   ', '  ',     0.51  ,  0.37   )
   SHANNON(370) = SHANNON_type('S+6 ',  '2p6 ',   'IV   ', '  ',     0.26  ,  0.12   )
   SHANNON(371) = SHANNON_type('S+6 ',  '2p6 ',   'VI   ', '  ',     0.43  ,  0.29   )
   SHANNON(372) = SHANNON_type('SB+3',  '5s2 ',   'IVpy ', '  ',     0.90  ,  0.76   )
   SHANNON(373) = SHANNON_type('SB+3',  '5s2 ',   'V    ', '  ',     0.94  ,  0.80   )
   SHANNON(374) = SHANNON_type('SB+3',  '5s2 ',   'VI   ', '  ',     0.90  ,  0.76   )
   SHANNON(375) = SHANNON_type('SB+5',  '4d10',   'VI   ', '  ',     0.74  ,  0.60   )
   SHANNON(376) = SHANNON_type('SC+3',  '3p6 ',   'VI   ', '  ',     0.885 ,  0.745  )
   SHANNON(377) = SHANNON_type('SC+3',  '3p6 ',   'VIII ', '  ',     1.010 ,  0.870  )
   SHANNON(378) = SHANNON_type('SE-2',  '4p6 ',   'VI   ', '  ',     1.84  ,  1.98   )
   SHANNON(379) = SHANNON_type('SE+4',  '4s2 ',   'VI   ', '  ',     0.64  ,  0.50   )
   SHANNON(380) = SHANNON_type('SE+6',  '3d10',   'IV   ', '  ',     0.42  ,  0.28   )
   SHANNON(381) = SHANNON_type('SE+6',  '3d10',   'VI   ', '  ',     0.56  ,  0.42   )
   SHANNON(382) = SHANNON_type('SI+4',  '2p6 ',   'IV   ', '  ',     0.40  ,  0.26   )
   SHANNON(383) = SHANNON_type('SI+4',  '2p6 ',   'VI   ', '  ',     0.540 ,  0.400  )
   SHANNON(384) = SHANNON_type('SM+2',  '4f6 ',   'VII  ', '  ',     1.36  ,  1.22   )
   SHANNON(385) = SHANNON_type('SM+2',  '4f6 ',   'VIII ', '  ',     1.41  ,  1.27   )
   SHANNON(386) = SHANNON_type('SM+2',  '4f6 ',   'IX   ', '  ',     1.46  ,  1.32   )
   SHANNON(387) = SHANNON_type('SM+3',  '4f5 ',   'VI   ', '  ',     1.098 ,  0.958  )
   SHANNON(388) = SHANNON_type('SM+3',  '4f5 ',   'VII  ', '  ',     1.16  ,  1.02   )
   SHANNON(389) = SHANNON_type('SM+3',  '4f5 ',   'VIII ', '  ',     1.219 ,  1.079  )
   SHANNON(390) = SHANNON_type('SM+3',  '4f5 ',   'IX   ', '  ',     1.272 ,  1.132  )
   SHANNON(391) = SHANNON_type('SM+3',  '4f5 ',   'XII  ', '  ',     1.38  ,  1.24   )
   SHANNON(392) = SHANNON_type('SN+4',  '4d10',   'IV   ', '  ',     0.69  ,  0.55   )
   SHANNON(393) = SHANNON_type('SN+4',  '4d10',   'V    ', '  ',     0.76  ,  0.62   )
   SHANNON(394) = SHANNON_type('SN+4',  '4d10',   'VI   ', '  ',     0.830 ,  0.690  )
   SHANNON(395) = SHANNON_type('SN+4',  '4d10',   'VII  ', '  ',     0.89  ,  0.75   )
   SHANNON(396) = SHANNON_type('SN+4',  '4d10',   'VIII ', '  ',     0.95  ,  0.81   )
   SHANNON(397) = SHANNON_type('SR+2',  '4p6 ',   'VI   ', '  ',     1.32  ,  1.18   )
   SHANNON(398) = SHANNON_type('SR+2',  '4p6 ',   'VII  ', '  ',     1.35  ,  1.21   )
   SHANNON(399) = SHANNON_type('SR+2',  '4p6 ',   'VIII ', '  ',     1.40  ,  1.26   )
   SHANNON(400) = SHANNON_type('SR+2',  '4p6 ',   'IX   ', '  ',     1.45  ,  1.31   )
   SHANNON(401) = SHANNON_type('SR+2',  '4p6 ',   'X    ', '  ',     1.50  ,  1.36   )
   SHANNON(402) = SHANNON_type('SR+2',  '4p6 ',   'XII  ', '  ',     1.58  ,  1.44   )
   SHANNON(403) = SHANNON_type('TA+3',  '5d2 ',   'VI   ', '  ',     0.86  ,  0.72   )
   SHANNON(404) = SHANNON_type('TA+4',  '5d1 ',   'VI   ', '  ',     0.82  ,  0.68   )
   SHANNON(405) = SHANNON_type('TA+5',  '5p6 ',   'VI   ', '  ',     0.78  ,  0.64   )
   SHANNON(406) = SHANNON_type('TA+5',  '5p6 ',   'VII  ', '  ',     0.83  ,  0.69   )
   SHANNON(407) = SHANNON_type('TA+5',  '5p6 ',   'VIII ', '  ',     0.88  ,  0.74   )
   SHANNON(408) = SHANNON_type('TB+3',  '4f8 ',   'VI   ', '  ',     1.063 ,  0.923  )
   SHANNON(409) = SHANNON_type('TB+3',  '4f8 ',   'VII  ', '  ',     1.12  ,  0.98   )
   SHANNON(410) = SHANNON_type('TB+3',  '4f8 ',   'VIII ', '  ',     1.180 ,  1.040  )
   SHANNON(411) = SHANNON_type('TB+3',  '4f8 ',   'IX   ', '  ',     1.235 ,  1.095  )
   SHANNON(412) = SHANNON_type('TB+4',  '4f7 ',   'VI   ', '  ',     0.90  ,  0.76   )
   SHANNON(413) = SHANNON_type('TB+4',  '4f7 ',   'VIII ', '  ',     1.02  ,  0.88   )
   SHANNON(414) = SHANNON_type('TC+4',  '4d3 ',   'VI   ', '  ',     0.785 ,  0.645  )
   SHANNON(415) = SHANNON_type('TC+5',  '4d2 ',   'VI   ', '  ',     0.74  ,  060    )
   SHANNON(416) = SHANNON_type('TC+7',  '4p6 ',   'IV   ', '  ',     0.51  ,  0.37   )
   SHANNON(417) = SHANNON_type('TC+7',  '4p6 ',   'VI   ', '  ',     0.70  ,  0.56   )
   SHANNON(418) = SHANNON_type('TE-2',  '5p6 ',   'VI   ', '  ',     2.07  ,  2.21   )
   SHANNON(419) = SHANNON_type('TE+4',  '5s2 ',   'III  ', '  ',     0.66  ,  0.52   )
   SHANNON(420) = SHANNON_type('TE+4',  '5s2 ',   'IV   ', '  ',     0.80  ,  0.66   )
   SHANNON(421) = SHANNON_type('TE+4',  '5s2 ',   'VI   ', '  ',     1.11  ,  0.97   )
   SHANNON(422) = SHANNON_type('TE+4',  '5s2 ',   'IV   ', '  ',     0.57  ,  0.43   )
   SHANNON(423) = SHANNON_type('TE+4',  '5s2 ',   'VI   ', '  ',     0.70  ,  0.56   )
   SHANNON(424) = SHANNON_type('TH+4',  '6p6 ',   'VI   ', '  ',     1.08  ,  0.94   )
   SHANNON(425) = SHANNON_type('TH+4',  '6p6 ',   'VIII ', '  ',     1.19  ,  1.05   )
   SHANNON(426) = SHANNON_type('TH+4',  '6p6 ',   'IX   ', '  ',     1.23  ,  1.09   )
   SHANNON(427) = SHANNON_type('TH+4',  '6p6 ',   'X    ', '  ',     1.27  ,  1.13   )
   SHANNON(428) = SHANNON_type('TH+4',  '6p6 ',   'XI   ', '  ',     1.32  ,  1.18   )
   SHANNON(429) = SHANNON_type('TH+4',  '6p6 ',   'XII  ', '  ',     1.35  ,  1.21   )
   SHANNON(430) = SHANNON_type('TI+2',  '3d2 ',   'VI   ', '  ',     1.00  ,  0.86   )
   SHANNON(431) = SHANNON_type('TI+3',  '3d1 ',   'VI   ', '  ',     0.810 ,  0.670  )
   SHANNON(432) = SHANNON_type('TI+4',  '3p6 ',   'IV   ', '  ',     0.56  ,  0.42   )
   SHANNON(433) = SHANNON_type('TI+4',  '3p6 ',   'V    ', '  ',     0.65  ,  0.51   )
   SHANNON(434) = SHANNON_type('TI+4',  '3p6 ',   'VI   ', '  ',     0.745 ,  0.605  )
   SHANNON(435) = SHANNON_type('TI+4',  '3p6 ',   'VIII ', '  ',     0.88  ,  0.74   )
   SHANNON(436) = SHANNON_type('TL+1',  '6s2 ',   'VI   ', '  ',     1.64  ,  1.50   )
   SHANNON(437) = SHANNON_type('TL+1',  '6s2 ',   'VIII ', '  ',     1.73  ,  0.59   )
   SHANNON(438) = SHANNON_type('TL+1',  '6s2 ',   'XII  ', '  ',     1.84  ,  1.70   )
   SHANNON(439) = SHANNON_type('TL+3',  '5d10',   'IV   ', '  ',     0.89  ,  0.75   )
   SHANNON(440) = SHANNON_type('TL+3',  '5d10',   'VI   ', '  ',     1.025 ,  0.885  )
   SHANNON(441) = SHANNON_type('TL+3',  '5d10',   'VIII ', '  ',     1.12  ,  0.98   )
   SHANNON(442) = SHANNON_type('TM+2',  '4f13',   'VI   ', '  ',     1.17  ,  1.03   )
   SHANNON(443) = SHANNON_type('TM+2',  '4f13',   'VII  ', '  ',     1.23  ,  1.09   )
   SHANNON(444) = SHANNON_type('TM+3',  '4f12',   'VI   ', '  ',     1.020 ,  0.880  )
   SHANNON(445) = SHANNON_type('TM+3',  '4f12',   'VIII ', '  ',     1.134 ,  0.994  )
   SHANNON(446) = SHANNON_type('TM+3',  '4f12',   'IX   ', '  ',     1.192 ,  1.052  )
   SHANNON(447) = SHANNON_type('U+3 ',  '5f3 ',   'VI   ', '  ',     1.165 ,  1.025  )
   SHANNON(448) = SHANNON_type('U+4 ',  '5f2 ',   'VI   ', '  ',     1.03  ,  0.89   )
   SHANNON(449) = SHANNON_type('U+4 ',  '5f2 ',   'VII  ', '  ',     1.09  ,  0.95   )
   SHANNON(450) = SHANNON_type('U+4 ',  '5f2 ',   'VIII ', '  ',     1.14  ,  1.00   )
   SHANNON(451) = SHANNON_type('U+4 ',  '5f2 ',   'IX   ', '  ',     1.19  ,  1.05   )
   SHANNON(452) = SHANNON_type('U+4 ',  '5f2 ',   'XII  ', '  ',     1.31  ,  1.17   )
   SHANNON(453) = SHANNON_type('U+5 ',  '5f1 ',   'VI   ', '  ',     0.90  ,  0.76   )
   SHANNON(454) = SHANNON_type('U+5 ',  '5f1 ',   'VII  ', '  ',     0.98  ,  0.84   )
   SHANNON(455) = SHANNON_type('U+6 ',  '6p6 ',   'II   ', '  ',     0.59  ,  0.45   )
   SHANNON(456) = SHANNON_type('U+6 ',  '6p6 ',   'IV   ', '  ',     0.66  ,  0.52   )
   SHANNON(457) = SHANNON_type('U+6 ',  '6p6 ',   'VI   ', '  ',     0.87  ,  0.73   )
   SHANNON(458) = SHANNON_type('U+6 ',  '6p6 ',   'VII  ', '  ',     0.85  ,  0.81   )
   SHANNON(459) = SHANNON_type('U+6 ',  '6p6 ',   'VIII ', '  ',     1.00  ,  0.86   )
   SHANNON(460) = SHANNON_type('V+2 ',  '3D3 ',   'VI   ', '  ',     0.93  ,  0.79   )
   SHANNON(461) = SHANNON_type('V+3 ',  '3d2 ',   'VI   ', '  ',     0.780 ,  0.640  )
   SHANNON(462) = SHANNON_type('V+4 ',  '3d1 ',   'V    ', '  ',     0.67  ,  0.53   )
   SHANNON(463) = SHANNON_type('V+4 ',  '3d1 ',   'VI   ', '  ',     0.72  ,  0.58   )
   SHANNON(464) = SHANNON_type('V+4 ',  '3d1 ',   'VIII ', '  ',     0.86  ,  0.72   )
   SHANNON(465) = SHANNON_type('V+5 ',  '3p6 ',   'IV   ', '  ',     0.495 ,  0.355  )
   SHANNON(466) = SHANNON_type('V+5 ',  '3p6 ',   'V    ', '  ',     0.60  ,  0.46   )
   SHANNON(467) = SHANNON_type('V+5 ',  '3p6 ',   'VI   ', '  ',     0.68  ,  0.54   )
   SHANNON(468) = SHANNON_type('W+4 ',  '5d2 ',   'VI   ', '  ',     0.80  ,  0.66   )
   SHANNON(469) = SHANNON_type('W+5 ',  '5d1 ',   'VI   ', '  ',     0.76  ,  0.62   )
   SHANNON(470) = SHANNON_type('W+6 ',  '5p6 ',   'IV   ', '  ',     0.56  ,  0.42   )
   SHANNON(471) = SHANNON_type('W+6 ',  '5p6 ',   'V    ', '  ',     0.65  ,  0.51   )
   SHANNON(472) = SHANNON_type('W+6 ',  '5p6 ',   'VI   ', '  ',     0.74  ,  0.60   )
   SHANNON(473) = SHANNON_type('XE+8',  '4d10',   'IV   ', '  ',     0.54  ,  0.40   )
   SHANNON(474) = SHANNON_type('XE+8',  '4d10',   'VI   ', '  ',     0.62  ,  0.48   )
   SHANNON(475) = SHANNON_type('Y+3 ',  '4p6 ',   'VI   ', '  ',     1.040 ,  0.900  )
   SHANNON(476) = SHANNON_type('Y+3 ',  '4p6 ',   'VII  ', '  ',     1.10  ,  0.96   )
   SHANNON(477) = SHANNON_type('Y+3 ',  '4p6 ',   'VIII ', '  ',     1.159 ,  1.019  )
   SHANNON(478) = SHANNON_type('Y+3 ',  '4p6 ',   'IX   ', '  ',     1.215 ,  1.075  )
   SHANNON(479) = SHANNON_type('YB+2',  '4f14',   'VI   ', '  ',     1.16  ,  1.02   )
   SHANNON(480) = SHANNON_type('YB+2',  '4f14',   'VII  ', '  ',     1.22  ,  1.08   )
   SHANNON(481) = SHANNON_type('YB+2',  '4f14',   'VIII ', '  ',     1.28  ,  1.14   )
   SHANNON(482) = SHANNON_type('YB+3',  '4f13',   'VI   ', '  ',     1.008 ,  0.868  )
   SHANNON(483) = SHANNON_type('YB+3',  '4f13',   'VII  ', '  ',     1.065 ,  0.925  )
   SHANNON(484) = SHANNON_type('YB+3',  '4f13',   'VIII ', '  ',     1.125 ,  0.985  )
   SHANNON(485) = SHANNON_type('YB+3',  '4f13',   'IX   ', '  ',     1.182 ,  1.042  )
   SHANNON(486) = SHANNON_type('ZN+2',  '3d10',   'IV   ', '  ',     0.74  ,  0.60   )
   SHANNON(487) = SHANNON_type('ZN+2',  '3d10',   'V    ', '  ',     0.82  ,  0.68   )
   SHANNON(488) = SHANNON_type('ZN+2',  '3d10',   'VI   ', '  ',     0.880 ,  0.740  )
   SHANNON(489) = SHANNON_type('ZN+2',  '3d10',   'VIII ', '  ',     1.04  ,  0.90   )
   SHANNON(490) = SHANNON_type('ZR+4',  '4p6 ',   'IV   ', '  ',     0.73  ,  0.59   )
   SHANNON(491) = SHANNON_type('ZR+4',  '4p6 ',   'V    ', '  ',     0.80  ,  0.66   )
   SHANNON(492) = SHANNON_type('ZR+4',  '4p6 ',   'VI   ', '  ',     0.86  ,  0.72   )
   SHANNON(493) = SHANNON_type('ZR+4',  '4p6 ',   'VII  ', '  ',     0.92  ,  0.78   )
   SHANNON(494) = SHANNON_type('ZR+4',  '4p6 ',   'VIII ', '  ',     0.98  ,  0.84   )
   SHANNON(495) = SHANNON_type('ZR+4',  '4p6 ',   'IX   ', '  ',     1.03  ,  0.89   )

end subroutine set_shannon

end module shannon_module



!!------------------------------------------------------------------------------------------------------------------------------!!
module magnetic_module

!# Extrait de :
!#        C. Kittel
!#        Physique de l'etat solide
!#        Dunod universite
!#        5me edition, 1983
!#
!# ec:     electronic configuration
!# level:  niveau
!# L:      moment magnetique orbital
!# S:      moment magnetique de spin
!# J:      moment magnetique effectif
!# g:      coef. de Lande
!# p:
!# gJ:
!#
!# ion   ec      CN      SP      CR      IR
!#

 implicit none

  private

  PUBLIC ::  set_mag_3d, set_mag_4f


  ! definitions
  TYPE, PUBLIC :: mag_4f_type
   CHARACTER (LEN=4)               :: ion
   CHARACTER (LEN=4)               :: ec        ! config.
   CHARACTER (LEN=6)               :: level     ! niveau
   CHARACTER (LEN=1)               :: L         ! L
   CHARACTER (LEN=3)               :: S         ! S
   CHARACTER (LEN=4)               :: J
   character (LEN=5)               :: g         !
   character (LEN=6)               :: p         !
   character (LEN=5)               :: gJ         !
  END TYPE mag_4f_type
  TYPE (mag_4f_type), ALLOCATABLE, DIMENSION(:), PUBLIC  :: mag_4f

  INTEGER, PARAMETER, PUBLIC :: nb_mag_4f_lines = 13

  TYPE, PUBLIC :: mag_3d_type
   CHARACTER (LEN=4)               :: ion
   CHARACTER (LEN=3)               :: ec        ! config.
   CHARACTER (LEN=5)               :: level     ! niveau
   CHARACTER (LEN=5)               :: g         !
   CHARACTER (LEN=4)               :: p_calc_g  ! p_calc = g(J(J+1)**0.5
   CHARACTER (LEN=3)               :: S
   character (LEN=4)               :: p_calc_S  ! p_calc = 2[S(S+1)]**0.5
   character (LEN=3)               :: p_exp     !
  END TYPE mag_3d_type
  TYPE (mag_3d_type), ALLOCATABLE, DIMENSION(:), PUBLIC  :: mag_3d
  INTEGER, PARAMETER, PUBLIC :: nb_mag_3d_lines = 13


contains

 subroutine set_mag_3d
  if (.not. ALLOCATED( mag_3d)) ALLOCATE( mag_3d(nb_mag_3d_lines))
   !                          ion     ec       level       g    p_calc_g    S        p_calc_S        p_exp
   mag_3d(  1) = mag_3d_type('Ti3+', '3d1', '2D3/2',  '4/5  ', '1.55', '1/2', '1.73', '1.8')
   mag_3d(  2) = mag_3d_type('V4+ ', '3d1', '2D3/2',  '4/5  ', '1.55', '1/2', '1.73', '1.8')
   mag_3d(  3) = mag_3d_type('V3+ ', '3d2', '3F2  ',  '2/3  ', '1.63', '2/2', '2.83', '2.8')
   mag_3d(  4) = mag_3d_type('Cr3+', '3d3', '4F3/2',  '2/5  ', '0.77', '3/2', '3.87', '3.8')
   mag_3d(  5) = mag_3d_type('V2+ ', '3d3', '4F3/2',  '2/5  ', '0.77', '3/2', '3.87', '3.8')
   mag_3d(  6) = mag_3d_type('Mn3+', '3d4', '5D0  ',  '     ', '0   ', '4/2', '4.90', '4.9')
   mag_3d(  7) = mag_3d_type('Cr2+', '3d4', '5D0  ',  '     ', '0   ', '4/2', '4.90', '4.9')
   mag_3d(  8) = mag_3d_type('Fe3+', '3d5', '6S5/2',  '2    ', '5.92', '5/2', '5.92', '4.9')
   mag_3d(  9) = mag_3d_type('Mn2+', '3d5', '6S5/2',  '2    ', '5.92', '5/2', '5.92', '4.9')
   mag_3d( 10) = mag_3d_type('Fe2+', '3d6', '5D4  ',  '3/2  ', '6.70', '4/2', '4.90', '5.4')
   mag_3d( 11) = mag_3d_type('Co2+', '3d7', '4F9/2',  '4/3  ', '6.63', '3/2', '3.87', '4.8')
   mag_3d( 12) = mag_3d_type('Ni2+', '3d8', '3F4  ',  '5/4  ', '5.59', '2/2', '2.83', '3.2')
   mag_3d( 13) = mag_3d_type('Cu2+', '3d9', '2D5/2',  '13/10', '3.55', '1/2', '1.73', '1.9')

 end subroutine set_mag_3d


 subroutine set_mag_4f

  if (.not. ALLOCATED( mag_4f)) ALLOCATE( mag_4f(nb_mag_4f_lines))
   !                          ion     ec       level       L    S      J          g           p        gJ

   mag_4f(  1) = mag_4f_type('Ce3+', '4f1 ',  '2F5/2 ',   '3', '1/2',  '5/2 ',   '0.857' ,  ' 2.535', ' 2.14' )
   mag_4f(  2) = mag_4f_type('Pr3+', '4f2 ',  '3H4   ',   '5', '1  ',  '4   ',   '0.8  ' ,  ' 3.577', ' 3.2 ' )
   mag_4f(  3) = mag_4f_type('Nd3+', '4f3 ',  '4I9/2 ',   '6', '3/2',  '9/2 ',   '0.727' ,  ' 3.617', ' 3.27' )
   mag_4f(  4) = mag_4f_type('Pm3+', '4f4 ',  '5I4   ',   '6', '2  ',  '4   ',   '0.6  ' ,  ' 2.68 ', ' 2.4 ' )
   mag_4f(  5) = mag_4f_type('Sm3+', '4f5 ',  '6H5/2 ',   '5', '5/2',  '5/2 ',   '0.286' ,  ' 0.845', ' 0.71' )
   mag_4f(  6) = mag_4f_type('Eu3+', '4f6 ',  '7F0   ',   '3', '3  ',  '0   ',   '     ' ,  ' 0    ', '     ' )
   mag_4f(  7) = mag_4f_type('Gd3+', '4f7 ',  '8S7.2 ',   '0', '7/2',  '7/2 ',   '2    ' ,  ' 7.937', ' 7   ' )
   mag_4f(  8) = mag_4f_type('Tb3+', '4f8 ',  '7F6   ',   '3', '3  ',  '6   ',   '1.5  ' ,  ' 9.72 ', ' 9   ' )
   mag_4f(  9) = mag_4f_type('Dy3+', '4f9 ',  '6H15/2',   '5', '5/2',  '15/2',   '1.33 ' ,  '10.64 ', '10   ' )
   mag_4f( 10) = mag_4f_type('Ho3+', '4f10',  '5I8   ',   '6', '3  ',  '8   ',   '1.25 ' ,  '10.607', '10   ' )
   mag_4f( 11) = mag_4f_type('Er3+', '4f11',  '4I15/2',   '6', '3/2',  '15/2',   '1.20 ' ,  ' 9.58 ', ' 9   ' )
   mag_4f( 12) = mag_4f_type('Tm3+', '4f12',  '3H6   ',   '5', '1  ',  '6   ',   '1.166' ,  ' 7.56 ', ' 7   ' )
   mag_4f( 13) = mag_4f_type('Yb3+', '4f13',  '2F7/2 ',   '3', '1/2',  '7/2 ',   '1.143' ,  ' 4.536', ' 4   ' )

end subroutine set_mag_4f




end module magnetic_module

