! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/API_IO_Formats.f90
! @brief     CFML IO Formats Fortran binding
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_Reflections_Utilities
  
  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env 
  
  use CFML_Reflections_Utilities,      only: Reflection_List_type
  
  !use API_Atom_TypeDef,               only: Atom_List_Type_p  
  use API_Crystallographic_Symmetry,  only: Space_Group_Type_p
  use API_Crystal_Metrics,            only: Crystal_Cell_Type_p

  type Reflection_List_type_p
     type(Reflection_List_type), pointer :: p
  end type Reflection_List_type_p
  
  implicit none

contains

  function reflections_utilities_hkl_uni_reflist(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr)        :: r
    type(tuple)        :: args 
    type(dict)         :: retval
    integer            :: num_args
    integer            :: ierror
    integer            :: ii

    type(Crystal_Cell_type_p)       :: cell_p
    type(Space_Group_type_p)        :: spg_p

    integer                         :: cell_p12(12), spg_p12(12)
    type(list)                      :: cell_obj, spg_obj

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)

    if (num_args /= 5) then
       call raise_exception(TypeError, "hkl_uni_reflist expects exactly 5 arguments")
       !Cell, SpG, lFriedel, min, max)
       call args%destroy
       return
    endif

  end function reflections_utilities_hkl_uni_reflist

end module API_Reflections_Utilities
