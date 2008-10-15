module WCFGL_atomic_table
!------------------------------------------------
! Written by Laurent C.Chapon 
! June 2004. 
! Updated : 
!------------------------------------------------ 
  implicit none 
  
  real,parameter, private :: c_default(4)=(/0.5,0.5,0.5,1.0/)  
    
  contains 
!----------------------------------------------------------------------------- 
  function get_color_from_symbol(Xsymb) result(color)
    character(len=2), intent(in) :: Xsymb
    real                         :: color(4)
    
    select case (Xsymb)
      case("H ")
        color=(/0.8,0.8,0.8,1.0/)
      case("He", "HE")           
        color=(/0.8,0.8,0.8,1.0/)
!-------------------------------------    
      case("Li","LI")
        color=(/0.8,0.2,0.45,1.0/)
      case("Be", "BE")
        color=(/0.1,0.5,0.65,1.0/)
      case("B ")
        color=(/0.4,0.9,0.25,1.0/)
      case("C ")
        color=(/0.6,0.6,0.6,1.0/)
      case("N ")
        color=(/0.2,0.7,0.2,1.0/)
      case("O ")
        color=(/0.8,0.1,0.1,1.0/)
      case("F ")
        color=(/0.23,0.23,0.89,1.0/)
      case("Ne","NE")
        color=c_default
!-------------------------------------
      case("Na","NA")
        color=(/0.7,0.7,0.15,1.0/)
      case("Mg","MG")
        color=(/0.7,0.0,0.74,1.0/)
      case("Al","AL")
        color=(/0.0,0.4,0.94,1.0/)
      case("Si","SI")
        color=(/0.5,0.5,0.24,1.0/)
      case("P ")
        color=(/0.0,0.56,0.14,1.0/)
      case("S ")
        color=(/0.9,0.9,0.05,1.0/)
      case("Cl","CL")
        color=(/0.05,0.9,0.1,1.0/)
      case("Ar","AR")
        color=c_default
!-------------------------------------
      case("K ")
        color=(/0.1,0.7,0.23,1.0/)
      case("Ca","CA")
        color=(/0.75,0.1,0.74,1.0/)
      case("Sc","SC")
        color=(/0.78,0.98,0.24,1.0/)
      case("Ti","TI")
        color=(/0.28,0.38,0.48,1.0/)
      case("V ")
        color=(/0.68,0.22,0.11,1.0/)
      case("Cr","CR")
        color=(/0.18,0.45,0.87,1.0/)
      case("Mn","MN")
        color=(/0.68,0.07,0.68,1.0/)
      case("Fe","FE")
        color=(/0.62,0.31,0.01,1.0/)
      case("Co","CO")
        color=(/0.78,0.18,0.28,1.0/)
      case("Ni","NI")
        color=(/0.28,0.68,0.81,1.0/)
      case("Cu","CU")
        color=(/0.05,0.05,0.59,1.0/)
      case("Zn","ZN")
        color=(/0.85,0.15,0.59,1.0/)
      case("Ga","GA")
        color=(/0.05,0.15,0.49,1.0/)
      case("Ge","GE")
        color=(/0.45,0.45,0.49,1.0/)
      case("As","AS")
        color=(/0.75,0.75,0.79,1.0/)
      case("Se","SE")
        color=(/0.95,0.15,0.59,1.0/)
      case("Br","BR")
        color=(/0.45,0.75,0.19,1.0/)
      case("Kr","KR")
        color=(/0.95,0.65,0.89,1.0/)
!-------------------------------------
      case("Rb","RR")
        color=(/0.55,0.05,0.59,1.0/)
      case("Sr","SR")
        color=(/0.75,0.65,0.29,1.0/)
      case("Y ")
        color=(/0.85,0.65,0.59,1.0/)
      case("Zr","ZR")
        color=(/0.95,0.15,0.89,1.0/)
      case("Nb","NB")
        color=(/0.15,0.35,0.59,1.0/)
      case("Mo","MO")
        color=(/0.45,0.55,0.39,1.0/)
      case("Tc","TC")
        color=(/0.45,0.15,0.69,1.0/)
      case("Ru","RU")
        color=(/0.35,0.95,0.19,1.0/)
      case("Rh","RH")
        color=(/0.75,0.15,0.19,1.0/)
      case("Pd","PD")
        color=(/0.85,0.15,0.39,1.0/)
      case("Ag","AG")
        color=(/0.95,0.15,0.79,1.0/)
      case("Cd","CD")
        color=(/0.05,0.15,0.89,1.0/)
      case("In","IN")
        color=(/0.05,0.65,0.79,1.0/)
      case("Sn","SN")
        color=(/0.05,0.85,0.49,1.0/)
      case("Sb","SB")
        color=(/0.25,0.85,0.79,1.0/)
      case("Te","TE")
        color=(/0.15,0.45,0.39,1.0/)
      case("I ")
        color=(/0.75,0.85,0.19,1.0/)
      case("Xe","XE")
        color=(/0.05,0.45,0.59,1.0/)
!------------------------------------- 
      case("Cs","CS")
        color=(/0.15,0.35,0.59,1.0/)
      case("Ba","BA")
        color=(/0.65,0.65,0.59,1.0/)
      case("La","LA")
        color=(/0.85,0.05,0.89,1.0/)
      case("Ce","CE")
        color=(/0.95,0.05,0.29,1.0/)
      case("Pr","PR")
        color=(/0.05,0.95,0.59,1.0/)
      case("Nd","ND")
        color=(/0.15,0.15,0.89,1.0/)
      case("Pm","PM")
        color=(/0.45,0.45,0.59,1.0/)
      case("Sm","SM")
        color=(/0.25,0.85,0.09,1.0/)
      case("Eu","EU")
        color=(/0.85,0.15,0.59,1.0/)
      case("Gd","GD")
        color=(/0.15,0.75,0.19,1.0/)
      case("Tb","TB")
        color=(/0.55,0.55,0.59,1.0/)
      case("Dy","DY")
        color=(/0.85,0.85,0.06,1.0/)
      case("Ho","HO")
        color=(/0.35,0.45,0.59,1.0/)
      case("Er","ER")
        color=(/0.45,0.75,0.59,1.0/)
      case("Tm","TM")
        color=(/0.75,0.15,0.09,1.0/)
      case("Yb","YB")
        color=(/0.12,0.65,0.49,1.0/)
      case("Lu","LU")
        color=(/0.85,0.85,0.59,1.0/)
      case("Hf","HF")
        color=(/0.95,0.35,0.69,1.0/)
      case("Ta","TA")
        color=(/0.55,0.15,0.59,1.0/)
      case("W ")
        color=(/0.75,0.95,0.79,1.0/)
      case("Re", "RE")
        color=(/0.55,0.45,0.59,1.0/)
      case("Os", "OS")
        color=(/0.85,0.95,0.59,1.0/)
      case("Ir", "IR")
        color=(/0.05,0.05,0.59,1.0/)
      case("Pt", "PT")
        color=(/0.55,0.05,0.59,1.0/)
      case("Au", "AU")
        color=(/0.55,0.45,0.59,1.0/)
      case("Hg", "HG")
        color=(/0.45,0.45,0.59,1.0/)
      case("Tl", "TL")
        color=(/0.15,0.05,0.59,1.0/)
      case("Pb", "PB")
        color=(/0.05,0.75,0.59,1.0/)
      case("Bi", "BI")
        color=(/0.85,0.65,0.59,1.0/)
      case("Po", "PO")
        color=(/0.35,0.35,0.59,1.0/)
      case("At", "AT")
        color=(/0.15,0.65,0.89,1.0/)
      case("Rn", "RN")
        color=(/0.35,0.35,0.59,1.0/)
      case default
        color=c_default
    end select  
 
  end function get_color_from_symbol
!----------------------------------------------------------------------------- 
  function get_radius_from_symbol(Xsymb) result(radius)
    character(len=2), intent(in) :: Xsymb
    real                         :: radius
    
    select case (Xsymb)
      case("H ")
        radius=0.37
      case("He", "HE")           
        radius=0.32
!-------------------------------------    
      case("Li","LI")
        radius=0.90
      case("Be", "BE")
        radius=0.50
      case("B ")
        radius=0.32
      case("C ")
        radius=0.77
      case("N ")
        radius=0.75
      case("O ")
        radius=1.26
      case("F ")
        radius=1.18
      case("Ne","NE")
        radius=0.38
!-------------------------------------
      case("Na","NA")
        radius=1.21
      case("Mg","MG")
        radius=0.86
      case("Al","AL")
        radius=0.60
      case("Si","SI")
        radius=1.11
      case("P ")
        radius=1.06
      case("S ")
        radius=1.84
      case("Cl","CL")
        radius=1.67
      case("Ar","AR")
        radius=0.71
!-------------------------------------
      case("K ")
        radius=1.52
      case("Ca","CA")
        radius=1.20
      case("Sc","SC")
        radius=0.95
      case("Ti","TI")
        radius=0.74
      case("V ")
        radius=0.81
      case("Cr","CR")
        radius=0.75
      case("Mn","MN")
        radius=0.81
      case("Fe","FE")
        radius=0.75
      case("Co","CO")
        radius=0.78
      case("Ni","NI")
        radius=0.70
      case("Cu","CU")
        radius=0.70
      case("Zn","ZN")
        radius=0.88
      case("Ga","GA")
        radius=0.69
      case("Ge","GE")
        radius=1.22
      case("As","AS")
        radius=1.19
      case("Se","SE")
        radius=1.98
      case("Br","BR")
        radius=1.95
      case("Kr","KR")
        radius=0.88
!-------------------------------------
      case("Rb","RR")
        radius=1.70
      case("Sr","SR")
        radius=1.35
      case("Y ")
        radius=1.09
      case("Zr","ZR")
        radius=0.86
      case("Nb","NB")
        radius=0.82
      case("Mo","MO")
        radius=0.78
      case("Tc","TC")
        radius=0.72
      case("Ru","RU")
        radius=0.76
      case("Rh","RH")
        radius=0.74
      case("Pd","PD")
        radius=0.90
      case("Ag","AG")
        radius=1.08
      case("Cd","CD")
        radius=1.09
      case("In","IN")
        radius=0.94
      case("Sn","SN")
        radius=1.41
      case("Sb","SB")
        radius=1.38
      case("Te","TE")
        radius=2.21
      case("I ")
        radius=2.06
      case("Xe","XE")
        radius=1.08
!------------------------------------- 
      case("Cs","CS")
        radius=1.85
      case("Ba","BA")
        radius=1.52
      case("La","LA")
        radius=1.40
      case("Ce","CE")
        radius=1.30
      case("Pr","PR")
        radius=1.20
      case("Nd","ND")
        radius=1.18
      case("Pm","PM")
        radius=1.17
      case("Sm","SM")
        radius=1.15
      case("Eu","EU")
        radius=1.25
      case("Gd","GD")
        radius=1.13
      case("Tb","TB")
        radius=1.12
      case("Dy","DY")
        radius=1.10
      case("Ho","HO")
        radius=1.10
      case("Er","ER")
        radius=1.08
      case("Tm","TM")
        radius=1.07
      case("Yb","YB")
        radius=1.15
      case("Lu","LU")
        radius=1.05
      case("Hf","HF")
        radius=0.85
      case("Ta","TA")
        radius=0.83
      case("W ")
        radius=0.76
      case("Re", "RE")
        radius=0.72
      case("Os", "OS")
        radius=0.70
      case("Ir", "IR")
        radius=0.76
      case("Pt", "PT")
        radius=0.82
      case("Au", "AU")
        radius=1.00
      case("Hg", "HG")
        radius=1.16
      case("Tl", "TL")
        radius=1.48
      case("Pb", "PB")
        radius=1.33
      case("Bi", "BI")
        radius=1.46
      case("Po", "PO")
        radius=1.08
      case("At", "AT")
        radius=0.76
      case("Rn", "RN")
        radius=1.20
      case default
        radius=1.00
    end select  
 
  end function get_radius_from_symbol
!----------------------------------------------------------------------------- 
end module WCFGL_atomic_table


