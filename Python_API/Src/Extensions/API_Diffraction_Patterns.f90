! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/API_Diffraction_Patterns.f90
! @brief     Diffraction Patterns utilities based on CrysFML
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_Diffraction_Patterns

  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use CFML_Diffraction_Patterns, only: &
       Diffraction_Pattern_Type

  use CFML_API_Calc_Powder_Patterm, only: &
       Calc_Powder_Pattern

  use API_reflection_utilities, only: &
       Reflection_List_type_p, &
       get_reflection_list_from_arg

  implicit none

  type Diffraction_Pattern_type_p
     type(Diffraction_Pattern_type), pointer :: p
  end type Diffraction_Pattern_type_p

contains

  function diffraction_patterns_compute_powder_pattern(self_ptr, args_ptr) result(r) bind(c)

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

    type(Reflection_List_type_p)     :: hkl_p
    type(Diffraction_Pattern_type_p) :: pattern_p
    integer                          :: pattern_p12(12)

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)

    if (num_args /= 2) then
       call raise_exception(TypeError, "Calc_Powder_Patter expects exactly 2 arguments")
       !Powpat_conditions, Reflectionlist
       call args%destroy
       return
    endif

    !call get_powpat_from_args(args, powpat_p, 0)
    
    call get_reflection_list_from_arg(args, hkl_p, 1)

    allocate(pattern_p%p)
    call Calc_Powder_Pattern(powpat_p%p,hkl_p%p,pattern_p%p)

    pattern_p12 = transfer(pattern_p,pattern_p12)
    ierror = list_create(index_obj)
    do ii=1,12
       ierror = index_obj%append(pattern_p12(ii))
    end do

    ierror = dict_create(retval)
    ierror = retval%setitem("address", index_obj)
    r = retval%get_c_ptr()

  end function diffraction_patterns_compute_powder_pattern


end module API_Diffraction_Patterns
 
