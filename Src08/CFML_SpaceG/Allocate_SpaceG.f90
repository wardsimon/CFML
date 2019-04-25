!!----
!!----
!!----
!!
SubModule (CFML_SpaceG) SPG_008
   Contains
   
   !!----
   !!---- ALLOCATE_SPACEG
   !!----
   !!---- 19/04/19
   !!
   Module Subroutine Allocate_SpaceG(D, Multip, Grp)
      !---- Arguments ----!
      integer,         intent(in)     :: d,multip
      class(Spg_Type), intent(in out) :: Grp
       
      !---- Local Variables ----! 
      integer :: i
      
      !> Init
      if (allocated(Grp%Op)) deallocate(Grp%Op)
      allocate(Grp%Op(multip))
      
      Grp%d=d
      Grp%multip=multip
      do i=1,multip
         call Allocate_Symm_Op(d,Grp%Op(i))
      end do
      
      if(allocated(Grp%Symb_Op)) deallocate(Grp%Symb_Op)
      allocate(Grp%Symb_Op(multip))
   End Subroutine Allocate_SpaceG

End SubModule SPG_008
   
