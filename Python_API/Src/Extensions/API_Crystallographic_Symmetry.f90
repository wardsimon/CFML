! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/API_Crystallographic_Symmetry.f90
! @brief     Symmetry groups Fortran binding
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_Crystallographic_Symmetry
  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, Write_SpaceGroup

  implicit none

  type Space_Group_Type_p
     type(Space_Group_Type), pointer :: p
  end type Space_Group_Type_p

contains

  subroutine get_object_from_arg(args, spgr_p)
    type(tuple) :: args
    type(object) :: arg_obj
    type(list) :: arg_list
    integer :: spgr_p12(12)
    type(Space_Group_type_p), intent(out) :: spgr_p

    integer :: ierror
    integer :: ii
    type(object) :: t

    ierror = args%getitem(arg_obj, 0)
    ierror = cast(arg_list, arg_obj)
    do ii=1,12
       ierror = arg_list%getitem(t, ii-1)
       ierror = cast(spgr_p12(ii), t)
    enddo
    spgr_p = transfer(spgr_p12, spgr_p)

  end subroutine get_object_from_arg

  !-------------------------------------------------------------------------
  ! Implementation of our Python methods
  !-------------------------------------------------------------------------
  ! @brief
  function crystallographic_symmetry_set_spacegroup(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(c_ptr) :: c_spgr
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ii
    integer :: ierror
    integer :: ci_spgr
    !
    type(object)       :: spgr_obj
    type(list)         :: index_obj
    character(len=:), allocatable :: spgr
    type(Space_Group_type_p) :: spgr_p
    integer :: spgr_p12(12)

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "create_space_group expects exactly 1 argument")
       call args%destroy
       return
    endif
    ierror = args%getitem(spgr_obj, 0)
    ierror = cast_nonstrict(spgr, spgr_obj)

    !
    allocate(spgr_p%p)
    call set_spacegroup(spgr,spgr_p%p)
    spgr_p12 = transfer(spgr_p, spgr_p12)
    ierror = list_create(index_obj)
    do ii=1,12
       ierror = index_obj%append(spgr_p12(ii))
    end do
    !deallocate(spgr_p%p)
    !
    ierror = dict_create(retval)
    ierror = retval%setitem("address", index_obj)
    r = retval%get_c_ptr()

  end function crystallographic_symmetry_set_spacegroup

  function crystallographic_symmetry_write_spacegroup(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_type_p) :: spgr_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_description expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_object_from_arg(args, spgr_p)

    !
    call Write_SpaceGroup(spgr_p%p, full=.true.)

    !
    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function crystallographic_symmetry_write_spacegroup

  function crystallographic_symmetry_get_latt_trans(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_type_p) :: spgr_p

    type(ndarray) :: latt_trans

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_latt_trans expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_object_from_arg(args, spgr_p)

    !
    ierror = ndarray_create(latt_trans, spgr_p%p%latt_trans)

    !
    ierror = dict_create(retval)

    !
    ierror = retval%setitem("latt_trans", latt_trans)

    r = retval%get_c_ptr()

  end function crystallographic_symmetry_get_latt_trans

  function crystallographic_symmetry_get_hexa(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_type_p) :: spgr_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_hexa expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_object_from_arg(args, spgr_p)

    !
    ierror = dict_create(retval)

    !
    ierror = retval%setitem("hexa", spgr_p%p%hexa)

    r = retval%get_c_ptr()

  end function crystallographic_symmetry_get_hexa

  function crystallographic_symmetry_get_numspg(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_type_p) :: spgr_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_numspg expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_object_from_arg(args, spgr_p)

    !
    ierror = dict_create(retval)

    !
    ierror = retval%setitem("numspg", spgr_p%p%numspg)

    r = retval%get_c_ptr()

  end function crystallographic_symmetry_get_numspg

  function crystallographic_symmetry_get_spg_symb(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_type_p) :: spgr_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_spg_symb expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_object_from_arg(args, spgr_p)

    !
    ierror = dict_create(retval)

    !
    ierror = retval%setitem("spg_symb", trim(spgr_p%p%spg_symb))

    r = retval%get_c_ptr()

  end function crystallographic_symmetry_get_spg_symb





     ! integer                                       :: NumSpg=0         ! Number of the Space Group
     ! character(len=20)                             :: SPG_Symb=" "     ! Hermann-Mauguin Symbol
     ! character(len=16)                             :: Hall=" "         ! Hall symbol
     ! character(len=90)                             :: gHall=" "        ! Generalised Hall symbol
     ! character(len=12)                             :: CrystalSys=" "   ! Crystal system
     ! character(len= 5)                             :: Laue=" "         ! Laue Class
     ! character(len= 5)                             :: PG=" "           ! Point group
     ! character(len= 5)                             :: Info=" "         ! Extra information
     ! character(len=90)                             :: SG_setting=" "   ! Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
     ! logical                                       :: Hexa=.false.     !
     ! character(len= 1)                             :: SPG_lat=" "      ! Lattice type
     ! character(len= 2)                             :: SPG_latsy=" "    ! Lattice type Symbol
     ! integer                                       :: NumLat=0         ! Number of lattice points in a cell
     ! real(kind=cp), allocatable,dimension(:,:)     :: Latt_trans       ! Lattice translations (3,12)
     ! character(len=51)                             :: Bravais=" "      ! String with Bravais symbol + translations
     ! character(len=80)                             :: Centre=" "       ! Alphanumeric information about the center of symmetry
     ! integer                                       :: Centred=0        ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
     ! real(kind=cp), dimension(3)                   :: Centre_coord=0.0 ! Fractional coordinates of the inversion centre
     ! integer                                       :: NumOps=0         ! Number of reduced set of S.O.
     ! integer                                       :: Multip=0         ! Multiplicity of the general position
     ! integer                                       :: Num_gen          ! Minimum number of operators to generate the Group
     ! type(Sym_Oper_Type), allocatable,dimension(:) :: SymOp            ! Symmetry operators (192)
     ! character(len=50),   allocatable,dimension(:) :: SymopSymb        ! Strings form of symmetry operators
     ! type(Wyckoff_Type)                            :: Wyckoff          ! Wyckoff Information
     ! real(kind=cp),dimension(3,2)                  :: R_Asym_Unit=0.0  ! Asymmetric unit in real space   A


end module API_Crystallographic_Symmetry
