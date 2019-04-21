!!----
!!----
!!----
!!
SubModule (CFML_Groups) CFML_GRP_021
   Contains
   
   !!----
   !!---- ALLOCATE_GROUP
   !!----
   !!---- 19/04/19
   !!
   Module Subroutine Allocate_Group(D, Multip, Grp)
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
         call Allocate_Operator(d,Grp%Op(i))
      end do
      
      if(allocated(Grp%Symb_Op)) deallocate(Grp%Symb_Op)
      allocate(Grp%Symb_Op(multip))
      
      return
   End Subroutine Allocate_Group

End SubModule CFML_GRP_021   
   
