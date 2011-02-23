! Tutorial for CrysFML
!./make_xtl; ./xtl test 
!
! The most primitive way to call CFML functions from C:
!   Pass anything as quickly as possible to a Fortran environment and leave the World of C 
!   Come back to C only at the very end and do not try to use types from CFML or more direct
!   calls to CFML in C.

subroutine xtl2cfml(length,line,abc,ang)
  use CFML_IO_Formats,                only: Readn_set_Xtal_Structure
  use CFML_Crystal_Metrics,           only: Crystal_Cell_Type !, Write_Crystal_Cell
  use CFML_Atom_TypeDef,              only: atom_list_type
  use CFML_Crystallographic_Symmetry, only: Space_Group_type

  implicit none
  
  integer                       :: length 
  character(length)             :: line
  real(kind=8),dimension(3)     :: abc,ang
   
  character(20)                 :: fmt,filnam
  type(Crystal_Cell_type)       :: cell
  type(Atom_List_type)          :: atoms
  type(Space_Group_type)        :: spgr
  
  write(fmt,'(a,i2,a)') "(a",length,",a)"  
  write(filnam,fmt) line,".cfl"   
  call Readn_Set_Xtal_Structure(filnam,cell,spgr,atoms,mode="CFL")
  abc=cell%cell
  ang=cell%ang
  !print *,abc,ang
  !call Write_Crystal_cell(cell)
end subroutine xtl2cfml