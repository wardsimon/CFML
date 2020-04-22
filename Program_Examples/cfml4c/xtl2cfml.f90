!!----  Tutorial for CrysFML
!!---- ./make_xtl; ./xtl test
!!----
!!----  The most primitive way to call CFML functions from C:
!!----
!!----    Pass anything as quickly as possible to a Fortran environment and leave the World of C
!!----    Come back to C only at the very end and do not try to use types from CFML or more direct
!!----    calls to CFML in C.
Subroutine Xtl2Cfml(length,line,abc,ang)
   !---- CrysFML Module ----!
   use CFML_IO_Formats,                only: Readn_set_Xtal_Structure
   use CFML_Crystal_Metrics,           only: Crystal_Cell_Type !, Write_Crystal_Cell
   use CFML_Atom_TypeDef,              only: atom_list_type
   use CFML_Crystallographic_Symmetry, only: Space_Group_type

   !---- Arguments ----!
   integer                           :: length
   character(len=length)             :: line
   real(kind=8),dimension(3)         :: abc
   real(kind=8),dimension(3)         :: ang

   !---- Local Variables ----!
   implicit none

   character(len=20)                 :: fmt,filnam
   type(Crystal_Cell_type)           :: cell
   type(Atom_List_type)              :: atoms
   type(Space_Group_type)            :: spgr


   write(fmt,'(a,i2,a)') "(a",length,",a)"
   write(filnam,fmt) line,".cfl"

   !> CrysFML routine
   call Readn_Set_Xtal_Structure(filnam,cell,spgr,atoms,mode="CFL")

   abc=cell%cell
   ang=cell%ang

   print *,abc,ang
   !call Write_Crystal_cell(cell)

End Subroutine Xtl2Cfml