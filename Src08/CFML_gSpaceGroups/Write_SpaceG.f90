!!----
!!----
!!----
!!
SubModule (CFML_gSpaceGroups) SPG_006
   Contains
   
   !!----
   !!---- WRITE_SPACEG_INFO
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Write_SpaceG_Info(Grp,Lun)
      !---- Arguments ----!
      class(Spg_Type),    intent(in)   :: Grp
      integer, optional,  intent(in)   :: lun
      
      !---- Local Variables ----!
      integer :: iout,i,j

      !> Init
      iout=6 !To be replaced by Fortran environment value
      if(present(lun)) iout=lun
      
      write(unit=iout,fmt="(a)")        "    General Space Group"
      write(unit=iout,fmt="(a)")        "    -------------------"
      write(unit=iout,fmt="(a,i4)")     "                  Op-Dimension: ",Grp%d
      write(unit=iout,fmt="(a,i4)")     "               Space-Dimension: ",Grp%d-1
      write(unit=iout,fmt="(a,i4)")     "                  Multiplicity: ",Grp%multip
      write(unit=iout,fmt="(a,i4)")     "                       MagType: ",Grp%mag_type
      write(unit=iout,fmt="(a,i4)")     "                        NumOps: ",Grp%numops
      write(unit=iout,fmt="(a,i4)")     "                       Centred: ",Grp%centred
      write(unit=iout,fmt="(a,i4)")     "                       Num_Lat: ",Grp%num_lat
      write(unit=iout,fmt="(a,i4)")     "                      Num_aLat: ",Grp%num_alat
      write(unit=iout,fmt="(a, a)")     "                Crystal system: ",Grp%Crystalsys
      write(unit=iout,fmt="(a, a)")     "  Crystallographic Point group: ",Grp%pg
      write(unit=iout,fmt="(a, a)")     "                    Laue class: ",Grp%laue
      write(unit=iout,fmt="(a,i4)")     "            Space Group number: ",Grp%numspg
      write(unit=iout,fmt="(a,i4)")     "        Shubnikov Group number: ",Grp%numshu
      write(unit=iout,fmt="(a, a)")     "            Space Group symbol: ",Grp%spg_symb
      write(unit=iout,fmt="(a, a)")     "        Shubnikov Group symbol: ",Grp%shu_symb
      write(unit=iout,fmt="(a, a)")     "       To Standard Space Group: ",Grp%mat2std
      write(unit=iout,fmt="(a, a)")     "   To Standard Shubnikov Group: ",Grp%Mat2Std_Shu
      
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
         write(unit=iout,fmt="(/a)")      "        Centring translations:"
         do i=1,Grp%num_lat
            write(unit=iout,fmt="(a,10a)") "             [ ",(trim(Rational_String(Grp%Lat_tr(j,i)))//" ",j=1,Grp%d-1),"]"
         end do
      end if
      
      if (Grp%num_alat > 0) then
         write(unit=iout,fmt="(/a)")      "            Anti-translations:"
         do i=1,Grp%num_alat
            write(unit=iout,fmt="(a,10a)") "             [ ",(trim(Rational_String(Grp%aLat_tr(j,i)))//" ",j=1,Grp%d-1),"]"
         end do
      end if
      
      write(unit=iout,fmt="(/a)")      "  Complete list of symmetry operators:"
      do i=1,Grp%Multip
         write(unit=iout,fmt="(i5,a)") i,"  ->  "//trim(Grp%Symb_Op(i))
      end do
   End Subroutine Write_SpaceG_Info

End SubModule SPG_006   
   
