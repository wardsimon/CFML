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
   Module Subroutine Allocate_Atom_List(N, A,Type_Atm)
      !---- Arguments ----!
      integer,             intent(in)       :: n    ! Atoms in the List
      type(Atlist_type),   intent(in out)   :: A    ! Objet to be allocated
      character(len=*),    intent(in)       :: Type_Atm

      !---- Local Variables ----!
      integer :: i
      ! Types :: Atm_Type, Atm_Std_Type, MAtm_Std_Type, Atm_Ref_Type, MAtm_Ref_Type
      type(Atm_Type)     , dimension(n)  :: Atm
      type(Atm_Std_Type) , dimension(n)  :: Atm_Std
      type(MAtm_Std_Type), dimension(n)  :: MAtm_Std
      type(Atm_Ref_Type) , dimension(n)  :: Atm_Ref
      type(MAtm_Ref_Type), dimension(n)  :: MAtm_Ref
      !> Init
      if (n <= 0) then
         A%natoms=0

         !> Deallocating atom list
         if (allocated(A%Active)) deallocate(A%Active)
         if (allocated(A%Atom)) deallocate(A%Atom)
         return
      end if

      !> Allocating variables
      Select Case(trim(l_case(Type_Atm)))
        Case("atm")
           allocate (A%atom(n),source=Atm)
        Case("atm_std")
           allocate (A%atom(n),source=Atm_Std)
        Case("matm_std")
           allocate (A%atom(n),source=MAtm_Std)
        Case("atm_ref")
            allocate (A%atom(n),source=Atm_Ref)
       Case("matm_ref")
           allocate (A%atom(n),source=MAtm_Ref)
      End Select

      allocate (A%active(n))
      A%active=.true.

      do i=1,n
         call Init_Atom_Type(A%Atom(i))
      end do

      A%natoms=n
   End Subroutine Allocate_Atom_list

End SubModule Atm_001