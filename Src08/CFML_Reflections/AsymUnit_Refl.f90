!!----
!!----
!!----
SubModule (CFML_Reflections) RFL_008
   Contains
   !!----
   !!---- Get_Asymm_Unit_H
   !!----    Provides a reflection equivalent to the input one but
   !!----    within the asymmetric unit
   !!----
   !!---- 21/06/2019 
   !!
   Module Function Get_Asymm_Unit_H(H,SpG) Result(k)
      !---- Arguments ----!
      integer, dimension (3),  intent(in) :: h
      class(SpG_Type),         intent(in) :: SpG
      integer, dimension(3)               :: k

      !---- Local Variables ----!
      integer, dimension(3), parameter    :: nul=[0,0,0]
      integer                             :: i
      integer, dimension(3)               :: kk
      integer, dimension(3,3)             :: Op

      k=h
      do i=1,SpG%NumOps
         Op=SpG%Op(i)%Mat(1:3,1:3)
         k=matmul(h,Op)
         kk=asu_h(k,SpG)
         if (h_equal(kk,NUL)) cycle
         k=kk
         exit
      end do
   End Function Get_Asymm_Unit_H
    
   !!----
   !!---- ASU_H
   !!----    Obtain an equivalent reflection in asymmetric unit using
   !!----    simple transformation rules for each crystal system.
   !!----    When these rules are not satisfied the output is the
   !!----    (0,0,0) reflection. 
   !!----
   !!----    For obtaining a reflection within the asymmetric unit
   !!----    given an input reflection the best is to use the function: 
   !!----    Get_Asymmetric_Unit_H
   !!----
   !!--<<
   !!----    We assumed that F(hkl)=F(-h -k -l).
   !!-->>
   !!----    If and error occurs, the function returns also (0,0,0).
   !!----
   !!---- 21/06/2019 
   !!
   Module Function Asu_H(H, SpG) Result(K)
      !---- Arguments ----!
      integer, dimension (3),  intent(in) :: h
      class(SpG_Type),         intent(in) :: SpG
      integer, dimension(3)               :: k

      !---- Local  variables ----!
      character(len=2)  :: inf
      character(len=4)  :: car

      !> Init
      k=0
      
      car=l_case(adjustl(SpG%Crystalsys))
      select case (car)
         case ('tric')
            k=asu_h_triclinic(h)
            
         case ('mono')
            k=asu_h_monoclinic(h,"b")
            !> k=asu_h_monoclinic(h,"c")
            !> k=asu_h_monoclinic(h,"a")
            
         case ('orth')
            k=asu_h_orthorhombic(h)
            
         case ('trig') 
            k=asu_h_trigonal(h,trim(SpG%Laue)) 
                  
         case ('tetr') 
            k=asu_h_tetragonal(h,trim(SpG%Laue))  
                
         case ('hexa')  
            k=asu_h_hexagonal(h,trim(SpG%Laue))  
                
         case ('cubi') 
            k=asu_h_cubic(h,trim(SpG%Laue))        
      end select
      
   End Function Asu_H
   
   !!--++
   !!--++ ASU_H_CUBIC
   !!--++    Obtain a reflection in asymmetric unit for Cubic
   !!--++
   !!--++ 21/06/2019 
   !!
   Module Function Asu_H_Cubic(H, Laue) Result(K)
      !---- Argument ----!
      integer, dimension(3), intent(in) :: H
      character(len=*),      intent(in) :: Laue
      integer, dimension(3)             :: K

      !---- Local Variable ----!
      character(len=4)      :: mod_laue
      integer, dimension(3) :: hh

      !> Init
      k=0
      mod_laue=l_case(adjustl(Laue))
      if (len_trim(mod_laue) == 0) then
         return
      end if

      select case(mod_laue)
         case("m-3  ")
            !> Laue: m-3 
            !> hkl: h>l, k>l, l>=0 ; hkk: k>=0 h>=k 
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select
            if (hh(3) >=0 .and. hh(1) >= hh(3) .and. hh(2) == hh(3)) k=hh
            if (hh(3) >=0 .and. hh(1) >  hh(3) .and. hh(2) >  hh(3)) k=hh

         case("m-3m ")
            !> Laue: m-3m
            !> hkl: h >=0, k >=0, l >=0, h >=k, k >=l 
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select
            if (hh(3) >= 0 .and. hh(2) >= hh(3) .and. hh(1) >= hh(2)) k=hh
      end select

   End Function Asu_H_Cubic

   !!--++
   !!--++ ASU_H_HEXAGONAL
   !!--++    Obtain a reflection in asymmetric unit for Hexagonal
   !!--++
   !!--++ 21/06/2019 
   !!
   Module Function Asu_H_Hexagonal(H, Laue) Result(K)
      !---- Argument ----!
      integer, dimension(3), intent(in) :: h
      character(len=*),      intent(in) :: Laue
      integer, dimension(3)             :: k

      !---- Local Variable ----!
      character(len=5)      :: mod_laue
      integer, dimension(3) :: hh

      k=0
      mod_laue=l_case(adjustl(Laue))
      if (len_trim(mod_laue) == 0) then
         return
      end if

      select case(mod_laue)
         case("6/m  ")
            !> Laue: 6/m 
            !> hkl: h>0,k>0,l>=0;  0kl k>=0,l>=0 
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select
            if (hh(1) > 0 .and. hh(2) > 0 .and. hh(3) >= 0) k=hh
            if (hh(1) == 0 .and. hh(2) >= 0 .and. hh(3) >= 0) k=hh

         case("6/mmm")
            !> Laue: 6/mmm 
            !> hkl: h >=0, k >=0, l >=0, h >=k 
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select
            if (hh(2) >=0 .and. hh(1) >= hh(2) .and. hh(3) >= 0) k=hh
      end select

   End Function Asu_H_Hexagonal
   
   !!--++
   !!--++ ASU_H_TRIGONAL
   !!--++    Obtain a reflection in asymmetric unit for Trigonal
   !!--++
   !!--++ 21/06/2019 
   !!
   Module Function Asu_H_Trigonal(H, Laue) Result(K)
      !---- Argument ----!
      integer, dimension(3), intent(in) :: h
      character(len=*),      intent(in) :: Laue
      integer, dimension(3)             :: k

      !---- Local Variable ----!
      character(len=4)      :: mod_laue
      integer, dimension(3) :: hh

      k=0
      mod_laue=l_case(adjustl(Laue))
      if (len_trim(mod_laue) == 0) then
         return
      end if

      select case(mod_laue)
         case("-3  ")
            !>Laue: -3 
            !> hkl: h+k>0, l>0 ; hk0:h>0, k>=0
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h

                  end select
               case (1:)
                  hh=h
            end select
            if (hh(1) == 0 .and. hh(2) == 0 .and. hh(3) > 0) k=hh
            if (hh(1)+hh(2) > 0 .and. hh(3) > 0 ) k=hh
            if (hh(1) > 0  .and. hh(2) >= 0  .and. hh(3) == 0) k=hh

         case("-3m ")
            !> Laue: -3m 
            !> hkl: h>=0, h>=k ; hhl: h>=0,l>=0 
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select
            if (hh(1) >= hh(2) .and.  hh(2) >= 0 ) k=hh
            if (hh(1) >= 0 .and. hh(2) > 0 .and. hh(3) > 0 ) k=hh
            if (hh(1) >= 0 .and. hh(2) == hh(1) .and. hh(3) >=0) k=hh

         case("-31m")
            !> Laue: -31m 
            !>hkl: h>=0,h>=k>0 ; h0l: h>=0,l>=0 
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                          hh=h
                  end select
               case (1:)
                  hh=h
            end select
            if (hh(1) >= hh(2) .and. hh(1) >=0 .and. hh(2) > 0) k=hh
            if (hh(1) >= 0 .and. hh(2) ==0 .and. hh(3) >= 0) k=hh
      end select

   End Function Asu_H_Trigonal
   
   !!--++
   !!--++ ASU_H_TETRAGONAL
   !!--++    Obtain a reflection in asymmetric unit for Tetragonal
   !!--++
   !!--++ 21/06/2019 
   !!
   Module Function Asu_H_Tetragonal(H,Laue) Result(K)
      !---- Argument ----!
      integer, dimension(3), intent(in) :: h
      character(len=*),      intent(in) :: Laue
      integer, dimension(3)             :: k

      !---- Local Variable ----!
      character(len=5)     :: mod_laue
      integer,dimension(3) :: hh

      k=0
      mod_laue=l_case(adjustl(Laue))
      if (len_trim(mod_laue) == 0) then
         return
      end if

      select case(mod_laue)
         case("4/m  ")
            !> Laue: 4/m 
            !> hkl: h >=0, l >=0, k >=0 if h = 0 
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select
            if (hh(1) == 0 .and. hh(2) >= 0 .and. hh(3) >=0) k=hh
            if (hh(1)  > 0 .and. hh(2) >  0 .and. hh(3) >=0) k=hh

         case("4/mmm")
            !> Laue: 4/mmm 
            !> hkl: h >=0, l >=0, h >=k   (k >=0) 
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(3) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select
            if (hh(1) >=0 .and. hh(2) >=0 .and. hh(3) >=0 .and. hh(1) >= hh(2)) k=hh
      end select

   End Function Asu_H_Tetragonal
   
   !!--++
   !!--++ Asu_H_Orthorhombic
   !!--++    Obtain a reflection in asymmetric unit for Orthorhombic
   !!--++    hkl: h >=0, k >=0, l >=0
   !!--++
   !!--++ 21/06/2019 
   !!
   Module Function Asu_H_Orthorhombic(H) Result(K)
      !---- Argument ----!
      integer, dimension(3), intent(in) :: h
      integer, dimension(3)             :: k

      !---- Local Variable ----!
      integer, dimension(3) :: hh

      k=0
      !>Laue: mmm 
      !> hkl: h >=0, k >=0, l >=0 
      select case (h(1))
         case (:-1)
            hh=-h
         case (0)
            select case (h(2))
               case (:-1)
                  hh=-h
               case (0)
                  if (h(3) >= 0) then
                     hh=h
                  else
                     hh=-h
                  end if
               case (1:)
                  hh=h
            end select
         case (1:)
            hh=h
      end select

      if (hh(1) >= 0 .and. hh(2) >= 0 .and. hh(3) >= 0) k=hh

   End Function Asu_H_Orthorhombic
   
   !!--++
   !!--++ Asu_H_Monoclinic
   !!--++    (PRIVATE)
   !!--++    Obtain a reflection in asymmetric unit for Monoclinic
   !!--++    Unique axis b: hkl: k >=0, l >=0    hk0: h >=0
   !!--++    Unique axis c: hkl: k >=0, l >=0    h0l: h >=0
   !!--++    Unique axis a: hkl: h >=0, l >=0    0kl: l >=0
   !!--++
   !!--++ Update: February - 2005
   !!
   Module Function Asu_H_Monoclinic(H, Axis) Result(K)
      !---- Argument ----!
      integer, dimension(3),      intent(in) :: h
      character(len=*), optional, intent(in) :: Axis
      integer, dimension(3)                  :: k

      !---- Local Variable ----!
      character(len=1)     :: ax
      integer,dimension(3) :: hh

      k=0
      if (present(axis)) then
         ax=l_case(adjustl(axis))
         if (ax ==" ") ax="b"
      else
         ax="b"
      end if

      select case (ax)
         !> Laue: 2/m     Unique Axis: b
         !> hkl: k >=0, l >=0    hk0: h >=0
         case ("b")
            select case (h(3))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(1) >=0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select

            if (hh(3) == 0) then
               if (hh(1) >=0 ) k=hh
            else
               if (hh(2) >=0 .and. hh(3) >=0) k=hh
            end if

         !> Laue: 2/m     Unique Axis: c 
         !> hkl: k >=0, l >=0    h0l: h >=0
         case ("c")
            select case (h(3))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(2))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(1) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select

            if (hh(2) == 0) then
               if (hh(1) >= 0) k=hh
            else
               if (hh(2) >=0 .and. hh(3) >=0) k=hh
            end if

         !> Laue: 2/m     Unique Axis: c 
         !> hkl: h >=0, l >=0    0kl: l >=0
         case ("a")
            select case (h(1))
               case (:-1)
                  hh=-h
               case (0)
                  select case (h(3))
                     case (:-1)
                        hh=-h
                     case (0)
                        if (h(2) >= 0) then
                           hh=h
                        else
                           hh=-h
                        end if
                     case (1:)
                        hh=h
                  end select
               case (1:)
                  hh=h
            end select

            if (hh(1) == 0) then
               if (hh(2) >= 0) k=hh
            else
               if (hh(1) >=0 .and. hh(3) >=0) k=hh
            end if

      end select
   End Function Asu_H_Monoclinic
   
   !!--++
   !!--++ ASU_H_TRICLINIC
   !!--++    Obtain a reflection in asymmetric unit for Triclinic
   !!--++    hkl: l >=0    hk0: h >=0    0k0: k >=0
   !!--++
   !!--++ 21/06/2019 
   !!
   Module Function Asu_H_Triclinic(H) Result(K)
      !---- Argument ----!
      integer, dimension(3), intent(in) :: h
      integer, dimension(3)             :: k

      k=0
      !>Laue: -1 
      !> hkl: l >=0    hk0: h >=0    0k0: k >=0
      select case (h(3))
         case (:-1)
            k=-h
         case (0)
            select case (h(1))
               case (:-1)
                  k=-h
               case (0)
                  if (h(2) < 0) then
                     k=-h
                  else
                     k=h
                  end if
               case (1:)
                  k=h
            end select
         case (1:)
            k=h
      end select
   End Function Asu_H_Triclinic
   
End SubModule RFL_008   