!!----
!!---- Menu: 0
!!---- Global items
 Module Menu_0
   !---- Use File ----!
   use CFML_GlobalDeps,  only: OPS,cp, Write_Date_Time
   character(len=5), public :: clear_string="clear" !Must be changed to "cls" for Windows
   integer, public :: i_out=11
   character(len=*),public,parameter :: fileout="CrysCalCon.log"
 End Module Menu_0