!!----
!!----
!!----
SubModule (CFML_Atoms) Init_Allocating_Atoms
   Contains

   !!----
   !!---- INIT_ATOM_TYPE
   !!----    Initialize Atom_Type
   !!----
   !!---- 12/06/2019
   !!
   Module Subroutine Init_Atom_Type(Atm,d)
      !---- Arguments ----!
      class(Atm_Type), intent(in out)   :: Atm
      integer,         intent(in)       :: d !Number of k-vectors

      Atm%Lab      =" "
      Atm%ChemSymb =" "
      Atm%Z        =0
      Atm%Mult     =0
      Atm%Charge   =0
      Atm%X        =0.0_cp
      Atm%Occ      =1.0_cp

      !> Thermal model
      Atm%UType    ="B"
      Atm%ThType   ="iso"
      Atm%U_iso    = 0.0_cp
      Atm%U        = 0.0_cp

      !> Magnetic
      Atm%Magnetic =.false.
      Atm%Moment   = 0.0_cp

      !> SFac
      Atm%SfacSymb =" "
      Atm%Ind_ff   = 0

      !> Info
      Atm%AtmInfo  = " "
      Atm%wyck      = " "

      select type (Atm)
         type is (atm_std_type)
            Atm%X_Std     = 0.0_cp
            Atm%Occ_Std   = 0.0_cp
            Atm%U_iso_Std = 0.0_cp
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp

         type is (matm_std_type)
            Atm%X_Std     = 0.0_cp
            Atm%Occ_Std   = 0.0_cp
            Atm%U_iso_Std = 0.0_cp
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp

            Atm%n_oc      = 0
            Atm%n_mc      = 0
            Atm%n_dc      = 0
            Atm%n_uc      = 0
            Atm%poc_q     = 0
            Atm%pmc_q     = 0
            Atm%pdc_q     = 0
            Atm%puc_q     = 0
            Atm%Ocs       = 0.0_cp
            Atm%Ocs_std   = 0.0_cp
            Atm%Mcs       = 0.0_cp
            Atm%Mcs_std   = 0.0_cp
            Atm%Dcs       = 0.0_cp
            Atm%Dcs_std   = 0.0_cp
            Atm%Ucs       = 0.0_cp
            Atm%Ucs_std   = 0.0_cp
            if(allocated (Atm%Xs)) deallocate(Atm%Xs)
            if(allocated (Atm%Us)) deallocate(Atm%Us)
            if(allocated (Atm%Moms)) deallocate(Atm%Moms)
            allocate(Atm%Xs(3+d),Atm%Us(3+d,3+d),Atm%Moms(3+d))
            Atm%Xs=0.0; Atm%Us=0.0; Atm%Moms=0.0

         type is (atm_ref_type)
            Atm%X_Std     = 0.0_cp
            Atm%Occ_Std   = 0.0_cp
            Atm%U_iso_Std = 0.0_cp
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp

            Atm%L_X       =0
            Atm%L_Occ     =0
            Atm%L_U_iso   =0
            Atm%L_U       =0
            Atm%L_moment  =0
            Atm%M_X       =0.0_cp
            Atm%M_Occ     =0.0_cp
            Atm%M_U_iso   =0.0_cp
            Atm%M_U       =0.0_cp
            Atm%M_moment  =0

         type is (Matm_ref_type)
            Atm%X_Std     = 0.0_cp
            Atm%Occ_Std   = 0.0_cp
            Atm%U_iso_Std = 0.0_cp
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp

            Atm%n_oc      = 0
            Atm%n_mc      = 0
            Atm%n_dc      = 0
            Atm%n_uc      = 0
            Atm%poc_q     = 0
            Atm%pmc_q     = 0
            Atm%pdc_q     = 0
            Atm%puc_q     = 0
            Atm%Ocs       = 0.0_cp
            Atm%Ocs_std   = 0.0_cp
            Atm%Mcs       = 0.0_cp
            Atm%Mcs_std   = 0.0_cp
            Atm%Dcs       = 0.0_cp
            Atm%Dcs_std   = 0.0_cp
            Atm%Ucs       = 0.0_cp
            Atm%Ucs_std   = 0.0_cp

            Atm%L_X       =0
            Atm%L_Occ     =0
            Atm%L_U_iso   =0
            Atm%L_U       =0
            Atm%M_X       =0.0_cp
            Atm%M_Occ     =0.0_cp
            Atm%M_U_iso   =0.0_cp
            Atm%M_U       =0.0_cp

            Atm%L_Ocs       = 0.0_cp
            Atm%L_Mcs       = 0.0_cp
            Atm%L_Dcs       = 0.0_cp
            Atm%L_Ucs       = 0.0_cp
            Atm%M_Ocs       = 0.0_cp
            Atm%M_Mcs       = 0.0_cp
            Atm%M_Dcs       = 0.0_cp
            Atm%M_Ucs       = 0.0_cp
            if(allocated (Atm%Xs)) deallocate(Atm%Xs)
            if(allocated (Atm%Us)) deallocate(Atm%Us)
            if(allocated (Atm%Moms)) deallocate(Atm%Moms)
            allocate(Atm%Xs(3+d),Atm%Us(3+d,3+d),Atm%Moms(3+d))
            Atm%Xs=0.0; Atm%Us=0.0; Atm%Moms=0.0
      end select

   End Subroutine Init_Atom_Type
   !!----
   !!---- ALLOCATE_ATOM_LIST
   !!----    Allocation of objet A of type atom_list.
   !!----    This procedure subroutine should be called before using an object of type atm_list
   !!----
   !!---- 12/06/2019
   !!
   Module Subroutine Allocate_Atom_List(N, A,Type_Atm,d)
      !---- Arguments ----!
      integer,             intent(in)       :: n    ! Atoms in the List
      type(Atlist_type),   intent(in out)   :: A    ! Objet to be allocated
      character(len=*),    intent(in)       :: Type_Atm
      integer,             intent(in)       :: d    !Number of k-vectors

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
      A%mcomp="crystal"

      do i=1,n
         call Init_Atom_Type(A%Atom(i),d)
      end do

      A%natoms=n
   End Subroutine Allocate_Atom_list

End SubModule Init_Allocating_Atoms