!!----
!!----
!!----
!!
SubModule (CFML_gSpaceGroups) SPG_006
   Contains

   !!----
   !!---- Write_SpaceGroup_Info
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Write_SpaceGroup_Info(Grp,Lun)
      !---- Arguments ----!
      class(Spg_Type),    intent(in)   :: Grp
      integer, optional,  intent(in)   :: lun

      !---- Local Variables ----!
      integer :: iout,i,j

      !> Init
      iout=6 !To be replaced by Fortran environment value
      if(present(lun)) iout=lun

      if(Grp%d-1 > 3) then !Superspace group

            write(unit=iout,fmt="(a)")                                        "     General SuperSpace Group"
            write(unit=iout,fmt="(a)")                                        "    -------------------------"
            write(unit=iout,fmt="(a,i4)")                                     "                  Op-Dimension: ",Grp%d
            write(unit=iout,fmt="(a,i4)")                                     "               Space-Dimension: ",Grp%d-1
            write(unit=iout,fmt="(a,i4)")                                     "                  Multiplicity: ",Grp%multip
            write(unit=iout,fmt="(a,i4)")                                     "                       MagType: ",Grp%mag_type
            write(unit=iout,fmt="(a,i4)")                                     "                        NumOps: ",Grp%numops
            write(unit=iout,fmt="(a,i4)")                                     "                       Centred: ",Grp%centred
            write(unit=iout,fmt="(a,i4)")                                     "     Num. Centring translation: ",Grp%num_lat
            write(unit=iout,fmt="(a,i4)")                                     "        Num. Anti-translations: ",Grp%num_alat
            if(len_trim(Grp%Crystalsys) /= 0)   write(unit=iout,fmt="(a, a)") "                Crystal system: ",trim(Grp%Crystalsys)
            if(len_trim(Grp%centre) /= 0)       write(unit=iout,fmt="(a, a)") "            Centre of symmetry: ",trim(Grp%centre)
            if(len_trim(Grp%pg) /=0 )           write(unit=iout,fmt="(a, a)") "  Crystallographic Point group: ",trim(Grp%pg)
            if(len_trim(Grp%laue) /= 0)         write(unit=iout,fmt="(a, a)") "                    Laue class: ",trim(Grp%laue)
            if(len_trim(Grp%setting) /= 0)      write(unit=iout,fmt="(a, a)") "     Setting w.r.t. or. gener.: ",trim(Grp%setting)
            if(len_trim(Grp%mat2std) /= 0)      write(unit=iout,fmt="(a, a)") "  To Standard SuperSpace Group: ",trim(Grp%mat2std)
            if(len_trim(Grp%mag_pg) /= 0)       write(unit=iout,fmt="(a, a)") "          Magnetic Point group: ",trim(Grp%mag_pg)
                                                write(unit=iout,fmt="(a,i4)") "     Parent Space Group number: ",Grp%Parent_num
            if(len_trim(Grp%Parent_spg) /= 0)   write(unit=iout,fmt="(a, a)") "            Parent Space Group: ",trim(Grp%Parent_spg)
            if(len_trim(Grp%tfrom_parent) /= 0) write(unit=iout,fmt="(a, a)") "         Transform from Parent: ",trim(Grp%tfrom_parent)
            Select Type (myGrp => Grp)
              type is(SuperSpaceGroup_Type)
                                                      write(unit=iout,fmt="(a,i4)") "          Bravais Class number: ",myGrp%Bravais_num
                 if(len_trim(Grp%mat2std) /= 0)       write(unit=iout,fmt="(a, a)") "                 Bravais Class: ",trim(myGrp%SSG_Bravais)
                                                      write(unit=iout,fmt="(a,i4)") "       SuperSpace Group number: ",myGrp%numspg
                 if(len_trim(myGrp%SSG_nlabel) /= 0)  write(unit=iout,fmt="(a, a)") "       SuperSpace Group  Label: ",trim(myGrp%SSG_nlabel)
                 if(len_trim(myGrp%SSG_symb) /= 0)    write(unit=iout,fmt="(a, a)") "       SuperSpace Group symbol: ",trim(myGrp%SSG_symb)
              Class default
                 write(unit=iout,fmt="(a, a)")     "       SuperSpace Group symbol: ","Unknown"

            End Select
      else
                                               write(unit=iout,fmt="(a)")        "    General Space Group"
                                               write(unit=iout,fmt="(a)")        "    -------------------"
                                               write(unit=iout,fmt="(a,i4)")     "                  Op-Dimension: ",Grp%d
                                               write(unit=iout,fmt="(a,i4)")     "               Space-Dimension: ",Grp%d-1
                                               write(unit=iout,fmt="(a,i4)")     "                  Multiplicity: ",Grp%multip
                                               write(unit=iout,fmt="(a,i4)")     "                       MagType: ",Grp%mag_type
                                               write(unit=iout,fmt="(a,i4)")     "                        NumOps: ",Grp%numops
                                               write(unit=iout,fmt="(a,i4)")     "                       Centred: ",Grp%centred
                                               write(unit=iout,fmt="(a,i4)")     "     Num. Centring translation: ",Grp%num_lat
                                               write(unit=iout,fmt="(a,i4)")     "        Num. Anti-translations: ",Grp%num_alat
            if(len_trim(Grp%centre) /= 0)      write(unit=iout,fmt="(a, a)")     "            Centre of symmetry: ",trim(Grp%centre)
            if(len_trim(Grp%Crystalsys) /= 0)  write(unit=iout,fmt="(a, a)")     "                Crystal system: ",trim(Grp%Crystalsys)
            if(len_trim(Grp%pg) /= 0)          write(unit=iout,fmt="(a, a)")     "  Crystallographic Point group: ",trim(Grp%pg)
            if(len_trim(Grp%laue) /= 0)        write(unit=iout,fmt="(a, a)")     "                    Laue class: ",trim(Grp%laue)
            if(Grp%numspg /= 0)                write(unit=iout,fmt="(a,i4)")     "            Space Group number: ",Grp%numspg
            if(Grp%numshu /= 0)                write(unit=iout,fmt="(a,i4)")     "        Shubnikov Group number: ",Grp%numshu
            if(len_trim(Grp%spg_symb) /= 0)    write(unit=iout,fmt="(a, a)")     "            Space Group symbol: ",trim(Grp%spg_symb)
            if(len_trim(Grp%bns_symb) /= 0)    write(unit=iout,fmt="(a, a)")     "    Shubnikov Group BNS-symbol: ",trim(Grp%bns_symb)
            if(len_trim(Grp%bns_num) /= 0)     write(unit=iout,fmt="(a, a)")     "    Shubnikov Group BNS- label: ",trim(Grp%bns_num)
            if(len_trim(Grp%og_symb) /= 0)     write(unit=iout,fmt="(a, a)")     "    Shubnikov Group  OG-symbol: ",trim(Grp%og_symb)
            if(len_trim(Grp%mat2std) /= 0)     write(unit=iout,fmt="(a, a)")     "       To Standard Space Group: ",trim(Grp%mat2std)
            if(len_trim(Grp%setting) /= 0)     write(unit=iout,fmt="(a, a)")     "     Setting w.r.t. or. gener.: ",trim(Grp%setting)
            if(len_trim(Grp%Mat2Std_Shu) /= 0) write(unit=iout,fmt="(a, a)")     "   To Standard Shubnikov Group: ",trim(Grp%Mat2Std_Shu)
            if(len_trim(Grp%Parent_spg) /= 0)  write(unit=iout,fmt="(a, a)")     "            Parent Space Group: ",trim(Grp%Parent_spg)
            if(len_trim(Grp%tfrom_parent) /= 0)write(unit=iout,fmt="(a, a)")     "         Transform from Parent: ",trim(Grp%tfrom_parent)
      end if

      if (Len_Trim(Grp%Generators_List) /= 0) Then
         write(unit=iout,fmt="(a, a)")   "               Generators List: ",Grp%generators_list
      end if

      if (Grp%centred == 1) then
         write(unit=iout,fmt="(a)")     "                  Centre_coord: none!"
      else
         write(unit=iout,fmt="(a,10a)") "                  Centre_coord: [ ",(trim(Rational_String(Grp%centre_coord(i)))//" ",i=1,Grp%d-1),"]"
      end if

      if (Grp%anticentred == 1) then
         write(unit=iout,fmt="(a)")     "             Anti-Centre_coord: none!"
      else
         write(unit=iout,fmt="(a,10a)") "             Anti-Centre_coord: [ ",(trim(Rational_String(Grp%anticentre_coord(i)))//" ",i=1,Grp%d-1),"]"
      end if

      if (Grp%num_lat > 0) then
         write(unit=iout,fmt="(/a)")      "         Centring translations:"
         do i=1,Grp%num_lat
            write(unit=iout,fmt="(a,10a)") "                               [ ",(trim(Rational_String(Grp%Lat_tr(j,i)))//" ",j=1,Grp%d-1),"]"
         end do
      end if

      if (Grp%num_alat > 0) then
         write(unit=iout,fmt="(/a)")      "             Anti-translations:"
         do i=1,Grp%num_alat
            write(unit=iout,fmt="(a,10a)") "                               [ ",(trim(Rational_String(Grp%aLat_tr(j,i)))//" ",j=1,Grp%d-1),"]"
         end do
      end if

      write(unit=iout,fmt="(/a)")      "  Complete list of symmetry operators:"
      do i=1,Grp%Multip
         write(unit=iout,fmt="(i5,a)") i,"  ->  "//trim(Grp%Symb_Op(i))
      end do
   End Subroutine Write_SpaceGroup_Info

End SubModule SPG_006

