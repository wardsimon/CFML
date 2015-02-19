
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

!-------------------------------------------------------
