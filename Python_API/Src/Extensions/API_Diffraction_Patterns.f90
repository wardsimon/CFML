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

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)

    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function diffraction_patterns_compute_powder_pattern


end module API_Diffraction_Patterns
