!!----
!!---- The painfull (and better) way to call CFML functions from C:
!!----   Have access to CFML types and functions directly from C, using a binding wrapper
!!----
Module BindCFML
   !---- Use Modules ----!
   use ISO_C_BINDING

   !---- Type Definition ----!
   Type, public, bind(C) :: F_Crystal_Cell_Type
      real(kind=C_FLOAT),dimension(3)   :: cell, ang
      real(kind=C_FLOAT),dimension(3)   :: cell_std, ang_std
      real(kind=C_FLOAT),dimension(3)   :: rcell, rang
      real(kind=C_FLOAT),dimension(3,3) :: GD,GR
      real(kind=C_FLOAT),dimension(3,3) :: Cr_Orth_cel
      real(kind=C_FLOAT),dimension(3,3) :: Orth_Cr_cel
      real(kind=C_FLOAT),dimension(3,3) :: BL_M
      real(kind=C_FLOAT),dimension(3,3) :: BL_Minv
      real(kind=C_FLOAT)                :: CellVol
      real(kind=C_FLOAT)                :: RCellVol
      character(kind=C_CHAR,len=1)      :: CartType
   End Type F_Crystal_Cell_Type

   Type, public, bind(C) :: F_Atom_Type
      character(kind=C_CHAR,len=10)      :: Lab
      character(kind=C_CHAR,len=2)       :: ChemSymb
      character(kind=C_CHAR,len=4)       :: SfacSymb
      logical(kind=C_BOOL)               :: Active
      integer(kind=C_INT)                :: Z
      integer(kind=C_INT)                :: Mult
      real(kind=C_FLOAT),dimension(3)    :: X
      real(kind=C_FLOAT),dimension(3)    :: X_Std
      real(kind=C_FLOAT),dimension(3)    :: MX
      integer(kind=C_INT),dimension(3)   :: LX
      real(kind=C_FLOAT)                 :: Occ
      real(kind=C_FLOAT)                 :: Occ_Std
      real(kind=C_FLOAT)                 :: MOcc
      integer(kind=C_INT)                :: LOcc
      real(kind=C_FLOAT)                 :: Biso
      real(kind=C_FLOAT)                 :: Biso_std
      real(kind=C_FLOAT)                 :: MBiso
      integer(kind=C_INT)                :: LBiso
      character(kind=C_CHAR,len=4)       :: Utype
      character(kind=C_CHAR,len=5)       :: ThType
      real(kind=C_FLOAT),dimension(6)    :: U
      real(kind=C_FLOAT),dimension(6)    :: U_std
      real(kind=C_FLOAT)                 :: Ueq
      real(kind=C_FLOAT),dimension(6)    :: MU
      integer(kind=C_INT),dimension(6)   :: LU
      real(kind=C_FLOAT)                 :: Charge
      real(kind=C_FLOAT)                 :: Moment
      integer(kind=C_INT),dimension(5)   :: Ind
      integer(kind=C_INT)                :: NVar
      real(kind=C_FLOAT),dimension(10)   :: VarF
      character(kind=C_CHAR,len=40)      :: AtmInfo
   End Type F_Atom_Type

   Type, public, bind(C) :: F_Atom_List_Type
      integer(kind=C_INT)                :: natoms
      type(C_PTR)                        :: atom
   End type F_Atom_List_Type

   Type, public, bind(C)  :: F_Wyck_Pos_Type
      integer(kind=C_INT)                          :: multp
      character(kind=C_CHAR,len= 6)                :: site
      integer(kind=C_INT)                          :: norb
      character(kind=C_CHAR,len=40)                :: str_orig
      character(kind=C_CHAR,len=40),dimension(48)  :: str_orbit
   End Type F_Wyck_Pos_Type

   Type, public, bind(C)  :: F_Wyckoff_Type
      integer(kind=C_INT)                            :: num_orbit
      type(F_Wyck_Pos_Type), dimension(26)           :: orbit
   End Type F_Wyckoff_Type

   Type, public, bind(C)  :: F_Sym_Oper_Type
      integer(kind=C_INT),dimension(3,3) :: Rot
      real(kind=C_FLOAT), dimension(3)   :: Tr
   End Type F_Sym_Oper_Type

   Type, public, bind(C)  :: F_Space_Group_Type
      integer(kind=C_INT)                           :: NumSpg           ! Number of the Space Group
      character(kind=C_CHAR,len=20)                 :: SPG_Symb         ! Hermann-Mauguin Symbol
      character(kind=C_CHAR,len=16)                 :: Hall             ! Hall symbol
      character(kind=C_CHAR,len=12)                 :: CrystalSys       ! Crystal system
      character(kind=C_CHAR,len= 5)                 :: Laue             ! Laue Class
      character(kind=C_CHAR,len= 5)                 :: PG               ! Point group
      character(kind=C_CHAR,len= 5)                 :: Info             ! Extra information
      character(kind=C_CHAR,len=80)                 :: SG_setting       ! Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
      logical(kind=C_BOOL)                          :: Hexa             !
      character(kind=C_CHAR,len= 1)                 :: SPG_lat          ! Lattice type
      character(kind=C_CHAR,len= 2)                 :: SPG_latsy        ! Lattice type Symbol
      integer(kind=C_INT)                           :: NumLat           ! Number of lattice points in a cell
      real(kind=C_FLOAT), dimension(3,12)           :: Latt_trans       ! Lattice translations
      character(kind=C_CHAR,len=51)                 :: Bravais          ! String with Bravais symbol + translations
      character(kind=C_CHAR,len=80)                 :: Centre           ! Alphanumeric information about the center of symmetry
      integer(kind=C_INT)                           :: Centred          ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
      real(kind=C_FLOAT), dimension(3)              :: Centre_coord     ! Fractional coordinates of the inversion centre
      integer(kind=C_INT)                           :: NumOps           ! Number of reduced set of S.O.
      integer(kind=C_INT)                           :: Multip           ! Multiplicity of the general position
      integer(kind=C_INT)                           :: Num_gen          ! Minimum number of operators to generate the Group
      type(F_Sym_Oper_Type), dimension(192)         :: SymOp            ! Symmetry operators
      character(kind=C_CHAR,len=40), dimension(192) :: SymopSymb        ! Strings form of symmetry operators
      type(F_Wyckoff_Type)                          :: Wyckoff          ! Wyckoff Information
      real(kind=C_FLOAT),dimension(3,2)             :: R_Asym_Unit      ! Asymmetric unit in real(kind=cp) space
   End Type F_Space_Group_Type

   private :: F2C_Crystal_Cell_Type
   public  :: C_Readn_Set_Xtal_Structure,BindC_COMPLAINS_THAT

 Contains

   !----- Private Subroutines -----!

   !!----
   !!---- Subroutine C2F_String(cstr,fstr,length)
   !!----
   !!
   Subroutine C2F_String(cstr,fstr,length)
      !---- Arguments ----!
      implicit none
      integer(kind=C_INT), value, intent(in)            :: length
      character(kind=C_CHAR), intent(in)                :: cstr(length)
      character(len=length), intent(out)                :: fstr

      !---- Local Variables ----!
      character(kind=C_CHAR,len=9)                      :: fmt1

      write(fmt1,'(A,I6,2A)') '(', length, 'A)'
      write(fstr,fmt1) cstr

      return
   End Subroutine C2F_String

   !!----
   !!---- Subroutine F2C_String(cstr,fstr,length)
   !!----
   !!----
   !!
   Subroutine F2C_String(cstr,fstr,length)
      !---- Arguments ----!
      implicit none
      integer(kind=C_INT), value, intent(in)            :: length
      character(kind=C_CHAR), intent(out)               :: cstr(length)
      character(len=length), intent(in)                 :: fstr

      !---- Local Variables ----!
      integer(kind=4)                                   :: i

      do i=1, length
         cstr(i) = fstr(i:i)
      end do

      return
   End Subroutine F2C_String

   !!----
   !!---- Subroutine C2F_Crystal_Cell_Type(f_in, f_out)
   !!----
   !!----
   !!
   Subroutine C2F_Crystal_Cell_Type(f_in, f_out)
      !----Modules ----!
      use CFML_Crystal_Metrics,           only: Crystal_Cell_Type !, Write_Crystal_Cell

      !---- Arguments ----!
      implicit none
      type(F_Crystal_Cell_Type), intent(in)                      :: f_in
      type(Crystal_Cell_Type), intent(out)                       :: f_out

      f_out%cell = f_in%cell
      f_out%ang = f_in%ang
      f_out%cell_std = f_in%cell_std
      f_out%ang_std = f_in%ang_std
      f_out%rcell = f_in%rcell
      f_out%rang = f_in%rang
      f_out%GD = f_in%GD
      f_out%GR = f_in%GR
      f_out%Cr_Orth_cel = f_in%Cr_Orth_cel
      f_out%Orth_Cr_cel = f_in%Orth_Cr_cel
      f_out%BL_M = f_in%BL_M
      f_out%BL_Minv = f_in%BL_Minv

      return
   End Subroutine C2F_Crystal_Cell_Type

   !!----
   !!---- Subroutine F2C_Crystal_Cell_Type(f_in, f_out)
   !!----
   !!----
   !!
   Subroutine F2C_Crystal_Cell_Type(f_in, f_out)
      !---- Modules ----!
      use CFML_Crystal_Metrics,           only: Crystal_Cell_Type !, Write_Crystal_Cell

      !---- Arguments ----!
      implicit none
      type(Crystal_Cell_Type),   intent(in)                        :: f_in
      type(F_Crystal_Cell_Type), intent(out)                       :: f_out

      f_out%cell = f_in%cell
      f_out%ang = f_in%ang
      f_out%cell_std = f_in%cell_std
      f_out%ang_std = f_in%ang_std
      f_out%rcell = f_in%rcell
      f_out%rang = f_in%rang
      f_out%GD = f_in%GD
      f_out%GR = f_in%GR
      f_out%Cr_Orth_cel = f_in%Cr_Orth_cel
      f_out%Orth_Cr_cel = f_in%Orth_Cr_cel
      f_out%BL_M = f_in%BL_M
      f_out%BL_Minv = f_in%BL_Minv

      return
   End Subroutine F2C_Crystal_Cell_Type

   !---- Public Subroutines ----!

   !!----
   !!---- Subroutine C_Readn_Set_Xtal_Structure(C_filnam,C_length,C_cell,C_spgr,C_atoms)
   !!----
   !!----
   !!
   Subroutine C_Readn_Set_Xtal_Structure(C_filnam,C_length,C_cell,C_spgr,C_atoms) bind (C,name='C_Readn_Set_Xtal_Structure')
      !---- Modules ----!
      use CFML_IO_Formats,                only: Readn_set_Xtal_Structure
      use CFML_Crystal_Metrics,           only: Crystal_Cell_Type
      use CFML_Atom_TypeDef,              only: atom_list_type
      use CFML_Crystallographic_Symmetry, only: Space_Group_type

      !---- Arguments ----!
      implicit none
      character(kind=C_CHAR), intent(inout)             :: C_filnam(*)
      integer(kind=C_INT), value, intent(in)            :: C_length
      type(F_Crystal_Cell_Type),intent(out)             :: C_cell
      type(F_Space_Group_Type),intent(out)              :: C_spgr
      type(F_Atom_List_Type),intent(out)                :: C_atoms

      !---- Local Vatriables ----!
      type(Crystal_Cell_Type)                         :: F_cell
      type(Space_Group_Type)                          :: F_spgr
      type(Atom_List_Type)                            :: F_atoms
      character(len=C_length)                         :: F_filnam

      call C2F_string(C_filnam, f_filnam, C_length)
      call Readn_Set_Xtal_Structure(F_filnam,F_cell,F_spgr,F_atoms,mode="CFL")
      call F2C_Crystal_Cell_Type(F_Cell,C_Cell)

      return
   End Subroutine C_Readn_Set_Xtal_Structure

End Module BindCFML
