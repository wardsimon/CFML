!!----
!!----
!!----
SubModule (CFML_Atoms) Atm_001
   Contains
   
   !!----
   !!---- ALLOCATE_ATOM_LIST
   !!----    Allocation of objet A of type atom_list. 
   !!----    This procedure subroutine should be called before using an object of type atm_list
   !!----
   !!---- 12/06/2019
   !!
   Module Subroutine Allocate_Atom_List(N, A)
      !---- Arguments ----!
      integer,             intent(in)       :: n    ! Atoms in the List
      class(Atlist_type),  intent(in out)   :: A    ! Objet to be allocated

      !---- Local Variables ----!
      integer :: i

      !> Init
      if (n <= 0) then
         A%natoms=0

         !> Deallocating atom list
         if (allocated(A%Active)) deallocate(A%Active)
         if (allocated(A%Atom)) deallocate(A%Atom) 
         return
      end if
         
      !> Allocating variables       
      allocate (A%Active(n))
      A%active=.true.

      allocate (A%atom(n)) 
      do i=1,n
         call Init_Atom_Type(A%Atom(i))
      end do
            
      A%natoms=n
   End Subroutine Allocate_Atom_list
   
End SubModule Atm_001   