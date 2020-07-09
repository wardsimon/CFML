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

    type(ndarray) :: latt_object

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
    write(*,*) spgr_p%p%Latt_trans
    ierror = ndarray_create(latt_object, spgr_p%p%Latt_trans)

    !
    ierror = dict_create(retval)
    ierror = retval%setitem("Array", latt_object)
    r = retval%get_c_ptr()


  end function crystallographic_symmetry_get_latt_trans


end module API_Crystallographic_Symmetry
