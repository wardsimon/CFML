!!----
!!----
!!----
!!
SubModule (CFML_Symmetry_Tables) SpgGener
   Contains
   !!----
   !!---- GET_IT_GENERATORS
   !!----
   !!----    Provides the string "gener" containing the list of the generators
   !!----    (as given in the IT Crystallography) corresponding to the space group
   !!----    of symbol "spg". In "spg" the Hermann-Mauguin symbol or the number of the
   !!----    space group should be given. The calling program is responsible of decoding
   !!----    the string "gener". Generator are given in the Jone's Faithful notation and
   !!----    the separator is the symbol ";". 
   !!----
   !!----    An example: R-3 c
   !!----    gener = " x+1/3,y+2/3,z+2/3; -y,x-y,z; -y,-x,z+1/2"
   !!----    The variable is the string contained between the quotes.
   !!----
   !!---- 20/04/2019 
   !!
   Module Function Get_IT_Generators(Spg) Result(StrGen)
      !---- Arguments ----!
      character(len=*), intent(in)  :: spg      ! Hermann_Mauguin symbol or number of S.Group
      character(len=:), allocatable :: StrGen

      !----  Local variables ----!
      integer                 :: i, ier, numg
      character(len=len(spg)) :: symb,sp

      !> Init
      Strgen=" "
      if (len_trim(it_spg_gen(1)) <=0) call set_IT_gen()


      !> IT Number?
      read(unit=spg,fmt=*,iostat=ier) numg
      if (ier == 0) then
         if (numg > 0 .and. numg <= 230) then
            Strgen=it_spg_gen(numg)(12:)
         end if
         return
      end if   
      
      !> H-M Symbol?
      symb=u_case(spg)
      do i=1,230
         sp=u_case(it_spg_gen(i)(1:10))
         if (trim(symb) == trim(sp)) then
            strgen=it_spg_gen(i)(12:)
            return
         end if
      end do

      return
   End Function Get_IT_Generators
   
   !!--++
   !!--++ SUBROUTINE SET_IT_GEN
   !!--++
   !!--++    Fills the components of the Spg_Gen character variable
   !!--++    Called once by the public subroutine Get_Generators
   !!--++
   !!--++ 20/04/2019 
   !!
   Module Subroutine Set_It_Gen()

      it_spg_gen(  1) =  "P 1       : x,y,z "
      it_spg_gen(  2) =  "P -1      : -x,-y,-z "
      it_spg_gen(  3) =  "P 2       : -x,y,-z "
      it_spg_gen(  4) =  "P 21      : -x,y+1/2,-z "
      it_spg_gen(  5) =  "C 2       : x+1/2,y+1/2,z; -x,y,-z "
      it_spg_gen(  6) =  "P m       : x,-y,z "
      it_spg_gen(  7) =  "P c       : x,-y,z+1/2 "
      it_spg_gen(  8) =  "C m       : x+1/2,y+1/2,z; x,-y,z "
      it_spg_gen(  9) =  "C c       : x+1/2,y+1/2,z; x,-y,z+1/2 "
      it_spg_gen( 10) =  "P 2/m     : -x,y,-z; -x,-y,-z "
      it_spg_gen( 11) =  "P 21/m    : -x,y+1/2,-z; -x,-y,-z "
      it_spg_gen( 12) =  "C 2/m     : x+1/2,y+1/2,z; -x,y,-z; -x,-y,-z "
      it_spg_gen( 13) =  "P 2/c     : -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 14) =  "P 21/c    : -x,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen( 15) =  "C 2/c     : x+1/2,y+1/2,z; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 16) =  "P 2 2 2   : -x,-y,z; -x,y,-z "
      it_spg_gen( 17) =  "P 2 2 21  : -x,-y,z+1/2; -x,y,-z+1/2 "
      it_spg_gen( 18) =  "P 21 21 2 : -x,-y,z; -x+1/2,y+1/2,-z "
      it_spg_gen( 19) =  "P 21 21 21: -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2 "
      it_spg_gen( 20) =  "C 2 2 21  : x+1/2,y+1/2,z; -x,-y,z+1/2; -x,y,-z+1/2 "
      it_spg_gen( 21) =  "C 2 2 2   : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z "
      it_spg_gen( 22) =  "F 2 2 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z "
      it_spg_gen( 23) =  "I 2 2 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z "
      it_spg_gen( 24) =  "I 21 21 21: x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2 "
      it_spg_gen( 25) =  "P m m 2   : -x,-y,z; x,-y,z "
      it_spg_gen( 26) =  "P m c 21  : -x,-y,z+1/2; x,-y,z+1/2 "
      it_spg_gen( 27) =  "P c c 2   : -x,-y,z; x,-y,z+1/2 "
      it_spg_gen( 28) =  "P m a 2   : -x,-y,z; x+1/2,-y,z "
      it_spg_gen( 29) =  "P c a 21  : -x,-y,z+1/2; x+1/2,-y,z "
      it_spg_gen( 30) =  "P n c 2   : -x,-y,z; x,-y+1/2,z+1/2 "
      it_spg_gen( 31) =  "P m n 21  : -x+1/2,-y,z+1/2; x+1/2,-y,z+1/2 "
      it_spg_gen( 32) =  "P b a 2   : -x,-y,z; x+1/2,-y+1/2,z "
      it_spg_gen( 33) =  "P n a 21  : -x,-y,z+1/2; x+1/2,-y+1/2,z "
      it_spg_gen( 34) =  "P n n 2   : -x,-y,z; x+1/2,-y+1/2,z+1/2 "
      it_spg_gen( 35) =  "C m m 2   : x+1/2,y+1/2,z; -x,-y,z; x,-y,z "
      it_spg_gen( 36) =  "C m c 21  : x+1/2,y+1/2,z; -x,-y,z+1/2; x,-y,z+1/2 "
      it_spg_gen( 37) =  "C c c 2   : x+1/2,y+1/2,z; -x,-y,z; x,-y,z+1/2 "
      it_spg_gen( 38) =  "A m m 2   : x,y+1/2,z+1/2; -x,-y,z; x,-y,z "
      it_spg_gen( 39) =  "A b m 2   : x,y+1/2,z+1/2; -x,-y,z; x,-y+1/2,z "
      it_spg_gen( 40) =  "A m a 2   : x,y+1/2,z+1/2; -x,-y,z; x+1/2,-y,z "
      it_spg_gen( 41) =  "A b a 2   : x,y+1/2,z+1/2; -x,-y,z; x+1/2,-y+1/2,z "
      it_spg_gen( 42) =  "F m m 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; x,-y,z "
      it_spg_gen( 43) =  "F d d 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; x+1/4,-y+1/4,z+1/4 "
      it_spg_gen( 44) =  "I m m 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; x,-y,z "
      it_spg_gen( 45) =  "I b a 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; x+1/2,-y+1/2,z "
      it_spg_gen( 46) =  "I m a 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; x+1/2,-y,z "
      it_spg_gen( 47) =  "P m m m   : -x,-y,z; -x,y,-z; -x,-y,-z "
      it_spg_gen( 48) =  "P n n n   : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 49) =  "P c c m   : -x,-y,z; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 50) =  "P b a n   : -x+1/2,-y+1/2,z; -x+1/2,y,-z; -x,-y,-z "
      it_spg_gen( 51) =  "P m m a   : -x+1/2,-y,z; -x,y,-z; -x,-y,-z "
      it_spg_gen( 52) =  "P n n a   : -x+1/2,-y,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen( 53) =  "P m n a   : -x+1/2,-y,z+1/2; -x+1/2,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 54) =  "P c c a   : -x+1/2,-y,z; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 55) =  "P b a m   : -x,-y,z; -x+1/2,y+1/2,-z; -x,-y,-z "
      it_spg_gen( 56) =  "P c c n   : -x+1/2,-y+1/2,z; -x,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen( 57) =  "P b c m   : -x,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen( 58) =  "P n n m   : -x,-y,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen( 59) =  "P m m n   : -x+1/2,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z "
      it_spg_gen( 60) =  "P b c n   : -x+1/2,-y+1/2,z+1/2; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 61) =  "P b c a   : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen( 62) =  "P n m a   : -x+1/2,-y,z+1/2; -x,y+1/2,-z; -x,-y,-z "
      it_spg_gen( 63) =  "C m c m   : x+1/2,y+1/2,z; -x,-y,z+1/2; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 64) =  "C m c a   : x+1/2,y+1/2,z; -x,-y+1/2,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen( 65) =  "C m m m   : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z; -x,-y,-z "
      it_spg_gen( 66) =  "C c c m   : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 67) =  "C m m a   : x+1/2,y+1/2,z; -x,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z "
      it_spg_gen( 68) =  "C c c a   : x+1/2,y+1/2,z; -x+1/2,-y,z; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen( 69) =  "F m m m   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; -x,-y,-z "
      it_spg_gen( 70) =  "F d d d   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+3/4,-y+3/4,z; -x+3/4,y,-z+3/4; -x,-y,-z "
      it_spg_gen( 71) =  "I m m m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; -x,-y,-z "
      it_spg_gen( 72) =  "I b a m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x+1/2,y+1/2,-z; -x,-y,-z "
      it_spg_gen( 73) =  "I b c a   : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen( 74) =  "I m m a   : x+1/2,y+1/2,z+1/2; -x,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z "
      it_spg_gen( 75) =  "P 4       : -x,-y,z; -y,x,z "
      it_spg_gen( 76) =  "P 41      : -x,-y,z+1/2; -y,x,z+1/4 "
      it_spg_gen( 77) =  "P 42      : -x,-y,z; -y,x,z+1/2 "
      it_spg_gen( 78) =  "P 43      : -x,-y,z+1/2; -y,x,z+3/4 "
      it_spg_gen( 79) =  "I 4       : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z "
      it_spg_gen( 80) =  "I 41      : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4 "
      it_spg_gen( 81) =  "P -4      : -x,-y,z; y,-x,-z "
      it_spg_gen( 82) =  "I -4      : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z "
      it_spg_gen( 83) =  "P 4/m     : -x,-y,z; -y,x,z; -x,-y,-z "
      it_spg_gen( 84) =  "P 42/m    : -x,-y,z; -y,x,z+1/2; -x,-y,-z "
      it_spg_gen( 85) =  "P 4/n     : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,-y,-z "
      it_spg_gen( 86) =  "P 42/n    : -x+1/2,-y+1/2,z; -y,x+1/2,z+1/2; -x,-y,-z "
      it_spg_gen( 87) =  "I 4/m     : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,-y,-z "
      it_spg_gen( 88) =  "I 41/a    : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+3/4,x+1/4,z+1/4; -x,-y,-z "
      it_spg_gen( 89) =  "P 4 2 2   : -x,-y,z; -y,x,z; -x,y,-z "
      it_spg_gen( 90) =  "P 4 21 2  : -x,-y,z; -y+1/2,x+1/2,z; -x+1/2,y+1/2,-z "
      it_spg_gen( 91) =  "P 41 2 2  : -x,-y,z+1/2; -y,x,z+1/4; -x,y,-z "
      it_spg_gen( 92) =  "P 41 21 2 : -x,-y,z+1/2; -y+1/2,x+1/2,z+1/4; -x+1/2,y+1/2,-z+1/4 "
      it_spg_gen( 93) =  "P 42 2 2  : -x,-y,z; -y,x,z+1/2; -x,y,-z "
      it_spg_gen( 94) =  "P 42 21 2 : -x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z+1/2 "
      it_spg_gen( 95) =  "P 43 2 2  : -x,-y,z+1/2; -y,x,z+3/4; -x,y,-z "
      it_spg_gen( 96) =  "P 43 21 2 : -x,-y,z+1/2; -y+1/2,x+1/2,z+3/4; -x+1/2,y+1/2,-z+3/4 "
      it_spg_gen( 97) =  "I 4 2 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z "
      it_spg_gen( 98) =  "I 41 2 2  : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; -x+1/2,y,-z+3/4 "
      it_spg_gen( 99) =  "P 4 m m   : -x,-y,z; -y,x,z; x,-y,z "
      it_spg_gen(100) =  "P 4 b m   : -x,-y,z; -y,x,z; x+1/2,-y+1/2,z "
      it_spg_gen(101) =  "P 42 c m  : -x,-y,z; -y,x,z+1/2; x,-y,z+1/2 "
      it_spg_gen(102) =  "P 42 n m  : -x,-y,z; -y+1/2,x+1/2,z+1/2; x+1/2,-y+1/2,z+1/2 "
      it_spg_gen(103) =  "P 4 c c   : -x,-y,z; -y,x,z; x,-y,z+1/2 "
      it_spg_gen(104) =  "P 4 n c   : -x,-y,z; -y,x,z; x+1/2,-y+1/2,z+1/2 "
      it_spg_gen(105) =  "P 42 m c  : -x,-y,z; -y,x,z+1/2; x,-y,z "
      it_spg_gen(106) =  "P 42 b c  : -x,-y,z; -y,x,z+1/2; x+1/2,-y+1/2,z "
      it_spg_gen(107) =  "I 4 m m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; x,-y,z "
      it_spg_gen(108) =  "I 4 c m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; x,-y,z+1/2 "
      it_spg_gen(109) =  "I 41 m d  : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; x,-y,z "
      it_spg_gen(110) =  "I 41 c d  : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; x,-y,z+1/2 "
      it_spg_gen(111) =  "P -4 2 m  : -x,-y,z; y,-x,-z; -x,y,-z "
      it_spg_gen(112) =  "P -4 2 c  : -x,-y,z; y,-x,-z; -x,y,-z+1/2 "
      it_spg_gen(113) =  "P -4 21 m : -x,-y,z; y,-x,-z; -x+1/2,y+1/2,-z "
      it_spg_gen(114) =  "P -4 21 c : -x,-y,z; y,-x,-z; -x+1/2,y+1/2,-z+1/2 "
      it_spg_gen(115) =  "P -4 m 2  : -x,-y,z; y,-x,-z; x,-y,z "
      it_spg_gen(116) =  "P -4 c 2  : -x,-y,z; y,-x,-z; x,-y,z+1/2 "
      it_spg_gen(117) =  "P -4 b 2  : -x,-y,z; y,-x,-z; x+1/2,-y+1/2,z "
      it_spg_gen(118) =  "P -4 n 2  : -x,-y,z; y,-x,-z; x+1/2,-y+1/2,z+1/2 "
      it_spg_gen(119) =  "I -4 m 2  : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; x,-y,z "
      it_spg_gen(120) =  "I -4 c 2  : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; x,-y,z+1/2 "
      it_spg_gen(121) =  "I -4 2 m  : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; -x,y,-z "
      it_spg_gen(122) =  "I -4 2 d  : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; -x+1/2,y,-z+3/4 "
      it_spg_gen(123) =  "P 4/m m m : -x,-y,z; -y,x,z; -x,y,-z; -x,-y,-z "
      it_spg_gen(124) =  "P 4/m c c : -x,-y,z; -y,x,z; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen(125) =  "P 4/n b m : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x+1/2,y,-z; -x,-y,-z "
      it_spg_gen(126) =  "P 4/n n c : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x+1/2,y,-z+1/2; -x,-y,-z "
      it_spg_gen(127) =  "P 4/m b m : -x,-y,z; -y,x,z; -x+1/2,y+1/2,-z; -x,-y,-z "
      it_spg_gen(128) =  "P 4/m n c : -x,-y,z; -y,x,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen(129) =  "P 4/n m m : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,y+1/2,-z; -x,-y,-z "
      it_spg_gen(130) =  "P 4/n c c : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen(131) =  "P 42/m m c: -x,-y,z; -y,x,z+1/2; -x,y,-z; -x,-y,-z "
      it_spg_gen(132) =  "P 42/m c m: -x,-y,z; -y,x,z+1/2; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen(133) =  "P 42/n b c: -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x+1/2,y,-z; -x,-y,-z "
      it_spg_gen(134) =  "P 42/n n m: -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x+1/2,y,-z+1/2; -x,-y,-z "
      it_spg_gen(135) =  "P 42/m b c: -x,-y,z; -y,x,z+1/2; -x+1/2,y+1/2,-z; -x,-y,-z "
      it_spg_gen(136) =  "P 42/m n m: -x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen(137) =  "P 42/n m c: -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x,y+1/2,-z; -x,-y,-z "
      it_spg_gen(138) =  "P 42/n c m: -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen(139) =  "I 4/m m m : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z; -x,-y,-z "
      it_spg_gen(140) =  "I 4/m c m : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z+1/2; -x,-y,-z "
      it_spg_gen(141) =  "I 41/a m d: x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+1/4,x+3/4,z+1/4; -x+1/2,y,-z+1/2; -x,-y,-z "
      it_spg_gen(142) =  "I 41/a c d: x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+1/4,x+3/4,z+1/4; -x+1/2,y,-z; -x,-y,-z "
      it_spg_gen(143) =  "P 3       : -y,x-y,z "
      it_spg_gen(144) =  "P 31      : -y,x-y,z+1/3 "
      it_spg_gen(145) =  "P 32      : -y,x-y,z+2/3 "
      it_spg_gen(146) =  "R 3       : x+1/3,y+2/3,z+2/3; -y,x-y,z "
      it_spg_gen(147) =  "P -3      : -y,x-y,z; -x,-y,-z "
      it_spg_gen(148) =  "R -3      : x+1/3,y+2/3,z+2/3; -y,x-y,z; -x,-y,-z "
      it_spg_gen(149) =  "P 3 1 2   : -y,x-y,z; -y,-x,-z "
      it_spg_gen(150) =  "P 3 2 1   : -y,x-y,z; y,x,-z "
      it_spg_gen(151) =  "P 31 1 2  : -y,x-y,z+1/3; -y,-x,-z+2/3 "
      it_spg_gen(152) =  "P 31 2 1  : -y,x-y,z+1/3; y,x,-z "
      it_spg_gen(153) =  "P 32 1 2  : -y,x-y,z+2/3; -y,-x,-z+1/3 "
      it_spg_gen(154) =  "P 32 2 1  : -y,x-y,z+2/3; y,x,-z "
      it_spg_gen(155) =  "R 3 2     : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z "
      it_spg_gen(156) =  "P 3 m 1   : -y,x-y,z; -y,-x,z "
      it_spg_gen(157) =  "P 3 1 m   : -y,x-y,z; y,x,z "
      it_spg_gen(158) =  "P 3 c 1   : -y,x-y,z; -y,-x,z+1/2 "
      it_spg_gen(159) =  "P 3 1 c   : -y,x-y,z; y,x,z+1/2 "
      it_spg_gen(160) =  "R 3 m     : x+1/3,y+2/3,z+2/3; -y,x-y,z; -y,-x,z "
      it_spg_gen(161) =  "R 3 c     : x+1/3,y+2/3,z+2/3; -y,x-y,z; -y,-x,z+1/2 "
      it_spg_gen(162) =  "P -3 1 m  : -y,x-y,z; -y,-x,-z; -x,-y,-z "
      it_spg_gen(163) =  "P -3 1 c  : -y,x-y,z; -y,-x,-z+1/2; -x,-y,-z "
      it_spg_gen(164) =  "P -3 m 1  : -y,x-y,z; y,x,-z; -x,-y,-z "
      it_spg_gen(165) =  "P -3 c 1  : -y,x-y,z; y,x,-z+1/2; -x,-y,-z "
      it_spg_gen(166) =  "R -3 m    : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z; -x,-y,-z "
      it_spg_gen(167) =  "R -3 c    : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z+1/2; -x,-y,-z "
      it_spg_gen(168) =  "P 6       : -y,x-y,z; -x,-y,z "
      it_spg_gen(169) =  "P 61      : -y,x-y,z+1/3; -x,-y,z+1/2 "
      it_spg_gen(170) =  "P 65      : -y,x-y,z+2/3; -x,-y,z+1/2 "
      it_spg_gen(171) =  "P 62      : -y,x-y,z+2/3; -x,-y,z "
      it_spg_gen(172) =  "P 64      : -y,x-y,z+1/3; -x,-y,z "
      it_spg_gen(173) =  "P 63      : -y,x-y,z; -x,-y,z+1/2 "
      it_spg_gen(174) =  "P -6      : -y,x-y,z; x,y,-z "
      it_spg_gen(175) =  "P 6/m     : -y,x-y,z; -x,-y,z; -x,-y,-z "
      it_spg_gen(176) =  "P 63/m    : -y,x-y,z; -x,-y,z+1/2; -x,-y,-z "
      it_spg_gen(177) =  "P 6 2 2   : -y,x-y,z; -x,-y,z; y,x,-z "
      it_spg_gen(178) =  "P 61 2 2  : -y,x-y,z+1/3; -x,-y,z+1/2; y,x,-z+1/3 "
      it_spg_gen(179) =  "P 65 2 2  : -y,x-y,z+2/3; -x,-y,z+1/2; y,x,-z+2/3 "
      it_spg_gen(180) =  "P 62 2 2  : -y,x-y,z+2/3; -x,-y,z; y,x,-z+2/3 "
      it_spg_gen(181) =  "P 64 2 2  : -y,x-y,z+1/3; -x,-y,z; y,x,-z+1/3 "
      it_spg_gen(182) =  "P 63 2 2  : -y,x-y,z; -x,-y,z+1/2; y,x,-z "
      it_spg_gen(183) =  "P 6 m m   : -y,x-y,z; -x,-y,z; -y,-x,z "
      it_spg_gen(184) =  "P 6 c c   : -y,x-y,z; -x,-y,z; -y,-x,z+1/2 "
      it_spg_gen(185) =  "P 63 c m  : -y,x-y,z; -x,-y,z+1/2; -y,-x,z+1/2 "
      it_spg_gen(186) =  "P 63 m c  : -y,x-y,z; -x,-y,z+1/2; -y,-x,z "
      it_spg_gen(187) =  "P -6 m 2  : -y,x-y,z; x,y,-z; -y,-x,z "
      it_spg_gen(188) =  "P -6 c 2  : -y,x-y,z; x,y,-z+1/2; -y,-x,z+1/2 "
      it_spg_gen(189) =  "P -6 2 m  : -y,x-y,z; x,y,-z; y,x,-z "
      it_spg_gen(190) =  "P -6 2 c  : -y,x-y,z; x,y,-z+1/2; y,x,-z "
      it_spg_gen(191) =  "P 6/m m m : -y,x-y,z; -x,-y,z; y,x,-z; -x,-y,-z "
      it_spg_gen(192) =  "P 6/m c c : -y,x-y,z; -x,-y,z; y,x,-z+1/2; -x,-y,-z "
      it_spg_gen(193) =  "P 63/m c m: -y,x-y,z; -x,-y,z+1/2; y,x,-z+1/2; -x,-y,-z "
      it_spg_gen(194) =  "P 63/m m c: -y,x-y,z; -x,-y,z+1/2; y,x,-z; -x,-y,-z "
      it_spg_gen(195) =  "P 2 3     : -x,-y,z; -x,y,-z; z,x,y "
      it_spg_gen(196) =  "F 2 3     : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y "
      it_spg_gen(197) =  "I 2 3     : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y "
      it_spg_gen(198) =  "P 21 3    : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y "
      it_spg_gen(199) =  "I 21 3    : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y "
      it_spg_gen(200) =  "P m -3    : -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z "
      it_spg_gen(201) =  "P n -3    : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; -x,-y,-z "
      it_spg_gen(202) =  "F m -3    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z "
      it_spg_gen(203) =  "F d -3    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+1/4,-y+1/4,z; -x+1/4,y,-z+1/4; z,x,y; -x,-y,-z "
      it_spg_gen(204) =  "I m -3    : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z "
      it_spg_gen(205) =  "P a -3    : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; -x,-y,-z "
      it_spg_gen(206) =  "I a -3    : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; -x,-y,-z "
      it_spg_gen(207) =  "P 4 3 2   : -x,-y,z; -x,y,-z; z,x,y; y,x,-z "
      it_spg_gen(208) =  "P 42 3 2  : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2 "
      it_spg_gen(209) =  "F 4 3 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z "
      it_spg_gen(210) =  "F 41 3 2  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y+1/2,z+1/2; -x+1/2,y+1/2,-z; z,x,y; y+3/4,x+1/4,-z+3/4 "
      it_spg_gen(211) =  "I 4 3 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z "
      it_spg_gen(212) =  "P 43 3 2  : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+1/4,x+3/4,-z+3/4 "
      it_spg_gen(213) =  "P 41 3 2  : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4 "
      it_spg_gen(214) =  "I 41 3 2  : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4 "
      it_spg_gen(215) =  "P -4 3 m  : -x,-y,z; -x,y,-z; z,x,y; y,x,z "
      it_spg_gen(216) =  "F -4 3 m  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,z "
      it_spg_gen(217) =  "I -4 3 m  : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,z "
      it_spg_gen(218) =  "P -4 3 n  : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,z+1/2 "
      it_spg_gen(219) =  "F -4 3 c  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,z+1/2 "
      it_spg_gen(220) =  "I -4 3 d  : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+1/4,x+1/4,z+1/4 "
      it_spg_gen(221) =  "P m -3 m  : -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z "
      it_spg_gen(222) =  "P n -3 n  : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; y,x,-z+1/2; -x,-y,-z "
      it_spg_gen(223) =  "P m -3 n  : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; -x,-y,-z "
      it_spg_gen(224) =  "P n -3 m  : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; y+1/2,x+1/2,-z; -x,-y,-z "
      it_spg_gen(225) =  "F m -3 m  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z "
      it_spg_gen(226) =  "F m -3 c  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; -x,-y,-z "

      it_spg_gen(227) =  &
      "F d -3 m  : x+1/2,y+1/2,z; x+1/2,y,z+1/2;-x+3/4,-y+1/4,z+1/2;-x+1/4,y+1/2,-z+3/4;z,x,y;y+3/4,x+1/4,-z+1/2;-x,-y,-z"

      it_spg_gen(228) =  &
      "F d -3 c  : x+1/2,y+1/2,z; x+1/2,y,z+1/2;-x+1/4,-y+3/4,z+1/2;-x+3/4,y+1/2,-z+1/4;z,x,y;y+3/4,x+1/4,-z;-x,-y,-z"

      it_spg_gen(229) =  "I m -3 m  : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z "
      it_spg_gen(230) =  "I a -3 d  : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4; -x,-y,-z "

      return
   End Subroutine Set_It_Gen
   
End SubModule SpgGener   