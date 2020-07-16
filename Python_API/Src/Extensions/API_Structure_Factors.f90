! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/API_Structure_Factors.f90
! @brief     CFML Structure Factors Fortran binding
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_Structure_Factors

  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use CFML_GlobalDeps,                 only: Cp
  use CFML_Structure_Factors,      only: &
       Structure_Factors, &
       Write_Structure_Factors

  use API_Crystallographic_Symmetry, only: &
       Space_Group_Type_p, &
       get_space_group_type_from_arg

  use API_Atom_TypeDef, only: &
       Atom_list_type_p, &
       get_atoms_list_from_arg

  use API_Reflections_Utilities, only: &
      Reflection_List_type_p, &
      get_reflection_list_from_arg

  implicit none

contains

  function structure_factors_structure_factors(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr)        :: r
    type(tuple)        :: args
    type(dict)         :: retval

    integer            :: num_args
    integer            :: ierror
    integer            :: ii

    type(list)   :: index_obj
    type(object) :: arg_obj

    type(Atom_list_type_p)       :: atom_list_p
    type(Space_Group_type_p)        :: spg_p
    type(Reflection_List_type_p)    :: reflection_list_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)

    if (num_args /= 3) then
       call raise_exception(TypeError, "structure_factors_structure_factors expects exactly 3 arguments")
       !@atom_list @space_group.as_fortran_object(), @reflection_list
       call args%destroy
       return
    endif

    call get_atoms_list_from_arg(args, atom_list_p, 0)

    call get_space_group_type_from_arg(args, spg_p, 1)

    call get_reflection_list_from_arg(args, reflection_list_p, 2)

    call Structure_Factors(atom_list_p%p, spg_p%p, reflection_list_p%p)

    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function structure_factors_structure_factors

  function structure_factors_write_structure_factors(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Reflection_List_type_p)    :: reflection_list_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "write_structure_factors expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_reflection_list_from_arg(args, reflection_list_p)

    !
    call Write_Structure_Factors(6, reflection_list_p%p)

    !
    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function structure_factors_write_structure_factors

end module API_Structure_Factors
