!!----
!!---- Program: SPG_INFO
!!----          Example of simple program using CFML
!!----
!!---- Author: Juan Rodriguez-Carvajal
!!---- Revision: Nov-2008
!!
Program SPG_Info
   !---- Use Modules ----!
   use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                             Write_SpaceGroup
                                          
   !---- Variables ----!                                       
   implicit none
   
   character(len=20)      :: spgr
   type(Space_Group_type) :: grp_espacial

   !---- Procedure ----!
   do
      write(unit=*,fmt="(a)") " => Enter a space group: "
      write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
      read(unit=*,fmt="(a)") spgr
      if (len_trim(spgr)==0) exit
   
      !> Set the Space Group Information   
      call set_spacegroup(spgr,grp_espacial)
      
      !> Wrtting the SpaceGroup Information
      call Write_SpaceGroup(grp_espacial, full=.true.)
   end do
   
End Program SPG_Info
