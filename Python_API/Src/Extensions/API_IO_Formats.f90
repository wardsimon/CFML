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

module API_IO_Formats
  
  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env 
  
  use CFML_IO_Formats,                only: Readn_set_Xtal_structure
  use CFML_Atom_TypeDef,              only: Atom_List_Type
  
  use API_Crystallographic_Symmetry,  only: Space_Group_Type_p
  use API_Crystal_Metrics,            only: Crystal_Cell_Type_p

  implicit none


  !type definitions This should be in the API_Atom_TypeDef
  type Atom_list_type_p
     type(Atom_list_type), pointer :: p
  end type Atom_list_type_p

contains
  
  !-------------------------------------------------------------------------
  ! Implementation of our Python methods
  !-------------------------------------------------------------------------
  ! @brief 
  function IO_formats_readn_set_xtal_structure(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr)        :: r
    type(tuple)        :: args 
    type(dict)         :: retval
    integer            :: num_args
    integer            :: ierror
    integer            :: ii

    type(object)                    :: filename_obj
    character(len=:), allocatable   :: filename

    type(Crystal_Cell_type_p)       :: cell_p
    type(Space_Group_type_p)        :: spg_p
    type(Atom_list_type_p)          :: a_p

    integer                         :: cell_p12(12), spg_p12(12), a_p12(12)
    type(list)                      :: cell_obj, spg_obj, a_obj

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)

    if (num_args /= 1) then
       call raise_exception(TypeError, "readn_set_xtal_structure expects exactly 1 argument: filename")
       call args%destroy
       return
    endif

    ierror = args%getitem(filename_obj, 0)
    ierror = cast(filename, filename_obj)

    allocate(cell_p%p, spg_p%p, a_p%p)
    call Readn_set_Xtal_structure(trim(filename),cell_p%p,spg_p%p,a_p%p,Mode="CIF")

    cell_p12 = transfer(cell_p, cell_p12)
    deallocate(cell_p%p)
    ierror = list_create(cell_obj)
    do ii=1,12
       ierror = cell_obj%append(cell_p12(ii))
    end do

    spg_p12   = transfer(spg_p, spg_p12)
    deallocate(spg_p%p)
    ierror = list_create(spg_obj)
    do ii=1,12
       ierror = spg_obj%append(spg_p12(ii))
    end do

    a_p12    = transfer(a_p, a_p12)
    deallocate(a_p%p)
    ierror = list_create(a_obj)
    do ii=1,12
       ierror = a_obj%append(a_p12(ii))
    end do

    ierror = dict_create(retval)
    ierror = retval%setitem("Cell", cell_obj)
    ierror = retval%setitem("SpG",  spg_obj)
    ierror = retval%setitem("A", A_obj)

    r = retval%get_c_ptr()

  end function IO_formats_readn_set_xtal_structure
end module API_IO_Formats
