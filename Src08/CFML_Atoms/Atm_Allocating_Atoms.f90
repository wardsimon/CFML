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
      Atm%wyck     = " "

      !> Free variables
      Atm%VarF     = 0.0
      Atm%active   = .true.

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
            Atm%n_bc      = 0
            Atm%n_mc      = 0
            Atm%n_dc      = 0
            Atm%n_uc      = 0
            Atm%poc_q     = 0
            Atm%pbc_q     = 0
            Atm%pmc_q     = 0
            Atm%pdc_q     = 0
            Atm%puc_q     = 0
            Atm%Ocs       = 0.0_cp
            Atm%Ocs_std   = 0.0_cp
            Atm%Bcs       = 0.0_cp
            Atm%Bcs_std   = 0.0_cp
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
            Atm%n_bc      = 0
            Atm%n_mc      = 0
            Atm%n_dc      = 0
            Atm%n_uc      = 0
            Atm%poc_q     = 0
            Atm%pbc_q     = 0
            Atm%pmc_q     = 0
            Atm%pdc_q     = 0
            Atm%puc_q     = 0
            Atm%Ocs       = 0.0_cp
            Atm%Ocs_std   = 0.0_cp
            Atm%Bcs       = 0.0_cp
            Atm%Bcs_std   = 0.0_cp
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
            Atm%L_Bcs       = 0.0_cp
            Atm%L_Mcs       = 0.0_cp
            Atm%L_Dcs       = 0.0_cp
            Atm%L_Ucs       = 0.0_cp
            Atm%M_Ocs       = 0.0_cp
            Atm%M_Bcs       = 0.0_cp
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
      integer :: i,ier
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
         case("atm")
            allocate(A%atom(n),source=Atm,stat=ier)

         case("atm_std")
            allocate(A%atom(n),source=Atm_Std,stat=ier)

         case("matm_std")
            allocate(A%atom(n),source=MAtm_Std,stat=ier)

         case("atm_ref")
            allocate(A%atom(n),source=Atm_Ref,stat=ier)

         case("matm_ref")
            allocate(A%atom(n),source=MAtm_Ref,stat=ier)
      end select

      allocate (A%active(n),stat=ier)

      if (ier /= 0) then
         Err_CFML%Ierr=1
         write(unit=Err_CFML%Msg,fmt="(a,i6,a)") "Error allocating atom List for N =",N," atoms"
         return
      end if

      A%active=.true.
      A%mcomp="crystal"

      do i=1,n
         call Init_Atom_Type(A%Atom(i),d)
      end do

      A%natoms=n
   End Subroutine Allocate_Atom_list

   !!----
   !!---- Module Subroutine Allocate_Atoms_Cell(Nasu,Mul,Dmax,Ac)
   !!----    integer, intent(in)                      :: nasu    !  In -> Number of atoms in asymmetric unit
   !!----    integer, intent(in)                      :: mul     !  In -> General multiplicity of the Space Group
   !!----    real(kind=cp),    intent(in)             :: dmax    !  In -> Maximun distance to be calculated
   !!----    type (atoms_cell_type), intent(in out)   :: Ac      !  In -> Object of type atoms_cell_type
   !!----                                                          Out -> Allocated and initialized object Ac
   !!----
   !!----    Allocation of objet "Ac" of type Atm_Cell. "Ac" contains
   !!----    components with ALLOCATABLE attribute with dimension depending
   !!----    on the input arguments "Nasu", "Mul" and "Dmax". The variables used for calculating the
   !!----    de dimensions are:
   !!--<<
   !!----          natcel=nasu*mul       and      id=idp=nint(0.74048*(dmax/r_atom)**3)
   !!-->>
   !!----    This subroutine should be called before using the subroutines of this module with
   !!----    dummy arguments of type Atm_Cell.
   !!----
   !!---- Update: February - 2005
   !!
   Module Subroutine Allocate_Atoms_Cell(Nasu,Mul,Dmax,Ac)
      !---- Arguments ----!
      integer,              intent(in)     :: nasu
      integer,              intent(in)     :: mul
      real(kind=cp),        intent(in)     :: dmax
      type (Atm_cell_type), intent(in out) :: Ac

      !---- local variables ----!
      integer :: natcel,id

      if(nasu <= 0) then !Deallocate the component of the object Ac
         if (allocated(Ac%Lab         ))   deallocate (Ac%Lab       )
         if (allocated(Ac%xyz         ))   deallocate (Ac%xyz        )
         if (allocated(Ac%charge      ))   deallocate (Ac%charge     )
         if (allocated(Ac%moment      ))   deallocate (Ac%moment     )
         if (allocated(Ac%var_free    ))   deallocate (Ac%var_free   )
         if (allocated(Ac%neighb      ))   deallocate (Ac%neighb     )
         if (allocated(Ac%neighb_atom ))   deallocate (Ac%neighb_atom)
         if (allocated(Ac%distance    ))   deallocate (Ac%distance   )
         if (allocated(Ac%trans       ))   deallocate (Ac%trans      )
         if (allocated(Ac%ddist       ))   deallocate (Ac%ddist      )
         if (allocated(Ac%ddlab       ))   deallocate (Ac%ddlab      )
      end if

      natcel=nasu*mul
      id=nint(0.74048*(dmax/r_atom)**3)
      id=max(id,natcel)

      if (.not. allocated(Ac%Lab         ))   allocate (Ac%Lab          (natcel))
      if (.not. allocated(Ac%xyz         ))   allocate (Ac%xyz         (3,natcel))
      if (.not. allocated(Ac%charge      ))   allocate (Ac%charge        (natcel))
      if (.not. allocated(Ac%moment      ))   allocate (Ac%moment        (natcel))
      if (.not. allocated(Ac%var_free    ))   allocate (Ac%var_free   (10,natcel))
      if (.not. allocated(Ac%neighb      ))   allocate (Ac%neighb        (natcel))
      if (.not. allocated(Ac%neighb_atom ))   allocate (Ac%neighb_atom(id,natcel))
      if (.not. allocated(Ac%distance    ))   allocate (Ac%distance   (id,natcel))
      if (.not. allocated(Ac%trans       ))   allocate (Ac%trans    (3,id,natcel))
      if (.not. allocated(Ac%ddist       ))   allocate (Ac%ddist      (natcel*id))
      if (.not. allocated(Ac%ddlab       ))   allocate (Ac%ddlab      (natcel*id))

      Ac%nat         = natcel
      Ac%Lab         = " "
      Ac%xyz         = 0.0
      Ac%charge      = 0.0
      Ac%moment      = 0.0
      Ac%var_free    = 0.0
      Ac%neighb      = 0
      Ac%neighb_atom = 0
      Ac%distance    = 0.0
      Ac%trans       = 0
      Ac%ddist       = 0.0
      Ac%ddlab       = " "
   End Subroutine Allocate_Atoms_Cell

End SubModule Init_Allocating_Atoms