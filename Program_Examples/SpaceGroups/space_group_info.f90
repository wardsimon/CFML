!!----
!!---- Program: SPG_INFO
!!----          Example of simple program using CFML
!!----
!!---- Author: Juan Rodriguez-Carvajal
!!---- Revision: February-2017
!!
Program SPG_Info
   !---- Use Modules ----!
   use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup,  &
                                             Write_SpaceGroup, Symmetry_Symbol, &
                                             Get_Generators_From_SpGSymbol,     &
                                             err_symm, err_symm_mess

   !---- Variables ----!
   implicit none

   character(len=20)      :: spgr
   type(Space_Group_type) :: grp_espacial
   integer,           dimension(10) :: point_op
   character(len=50), dimension(10) :: gen
   character(len=70)                :: op_symb
   character(len=6),  dimension(10) :: sgen
   integer :: ngen,i,j

   !---- Procedure ----!
   do
      write(unit=*,fmt="(a)") " => Enter a space group: "
      write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
      read(unit=*,fmt="(a)") spgr
      if (len_trim(spgr)==0) exit

      !> Setting the Space Group Information
      call set_spacegroup(spgr,grp_espacial)
      if(err_symm) then
        write(unit=*,fmt="(a)") "   ERROR: "//trim(err_symm_mess)
      end if
      !> Writing the SpaceGroup Information
      call Write_SpaceGroup(grp_espacial, full=.true.)
      ! Pointing out the operators in the HM-symbol
      call Get_Generators_From_SpGSymbol(grp_espacial,gen,point_op,ngen)
      write(unit=*,fmt="(/,a,i3,a)") " => Generators of the Space Group: # ",grp_espacial%NumSpg,"  "//trim(grp_espacial%SPG_Symb)
      do i=1,ngen
        j=point_op(i)
        call Symmetry_Symbol(grp_espacial%SymopSymb(j),op_symb)
        write(unit=*,fmt="(a,i3,a)") "    Generator #",i, "  "//gen(i)(1:30)//"Symbol: "//trim(op_symb)
      end do
      write(unit=*,fmt="(a)") "   "
   end do

End Program SPG_Info
