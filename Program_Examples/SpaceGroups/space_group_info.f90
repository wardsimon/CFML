   Program SPG_Info
     use Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                        Write_SpaceGroup
     implicit none
     character(len=20)   :: spgr
     type(Space_Group_type) :: grp_espacial

     do
       write(unit=*,fmt="(a)") " => Enter a space group: "
       write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
       read(unit=*,fmt="(a)") spgr
       if (len_trim(spgr)==0) exit
       call set_spacegroup(spgr,grp_espacial)
       call Write_SpaceGroup(grp_espacial, full=.true.)
     end do
   End Program SPG_Info
