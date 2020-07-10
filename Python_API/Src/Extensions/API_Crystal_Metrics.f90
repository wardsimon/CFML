! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/API_Crystal_Metrics.f90
! @brief     Crystall Cell Fortran binding
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_Crystal_Metrics
  
  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use CFML_GlobalDeps,       only : Cp
  use CFML_Crystal_Metrics,  only: &
       Crystal_Cell_Type, &
       Write_Crystal_Cell
  
  implicit none

  !type definitions 
  type Crystal_Cell_type_p
     type(Crystal_Cell_type), pointer :: p
  end type Crystal_Cell_type_p
  
contains 

  subroutine get_cell_from_arg(args, cell_p)
    type(tuple), intent(in)                :: args
    type(Crystal_Cell_type_p), intent(out) :: cell_p
    
    type(object) :: arg_obj
    type(list) :: arg_list
    integer :: cell_p12(12)
  
    integer :: ierror
    integer :: ii
    type(object) :: t

    ierror = args%getitem(arg_obj, 0)
    ierror = cast(arg_list, arg_obj)
    do ii=1,12
       ierror = arg_list%getitem(t, ii-1)
       ierror = cast(cell_p12(ii), t)
    enddo
    cell_p = transfer(cell_p12, cell_p)

  end subroutine get_object_from_arg

  
  ! @brief Create a crystal cell from a set of lattice parameters
  function crystal_metrics_set_crystal_cell(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    
    type(tuple) :: args
    type(dict)  :: retval
    integer     :: num_args
    integer     :: ierror
    integer     :: ii
    type(list)  :: index_obj

    type(object)                :: cellv_obj, angl_obj
    real(kind=cp), dimension(3) :: cellv, angl
    type(Crystal_Cell_type_p)   :: cell_p
    integer                     :: cell_p12(12)
       
    r = C_NULL_PTR      
    call unsafe_cast_from_c_ptr(args, args_ptr)
    
    ierror = args%len(num_args)
    
    if (num_args < 2) then
       call raise_exception(TypeError, "set_crystal_Cell expects at least 2 arguments")
       call args%destroy
       return
    endif

    !> @todo optional arguments Cartype,Scell,Sangl
    if (num_args > 2) then
       call raise_exception(TypeError, "optional arguments not yet available for set_crystal_Cell")
       call args%destroy
       return
    endif
     
    ierror = args%getitem(cellv_obj, 0)        !cellv  In -> a,b,c
    ierror = cast_nonstrict(cellv, cellv_obj) 
    ierror = args%getitem(angl_obj, 1)         !angl   In -> angles of cell parameters
    ierror = cast_nonstrict(angl, angl_obj)

    call Set_Crystal_Cell(cellv, angl, cell_p%p)

    cell_p12 = transfer(cell_p,cell_p12)
    ierror = list_create(index_obj)
    do ii=1,12
       ierror = index_obj%append(cell_p12(ii))
    end do
    deallocate(cell_p%p)

    ierror = dict_create(retval)
    ierror = retval%setitem("address", index_obj)
    r = retval%get_c_ptr()

  end function crystal_metrics_set_crystal_cell

  function crystal_metrics_write_crystal_cell(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    
    type(tuple) :: args
    type(dict)  :: retval
    integer     :: num_args
    integer     :: ierror

    type(Crystal_Cell_type_p)   :: cell_p

    r = C_NULL_PTR      
    call unsafe_cast_from_c_ptr(args, args_ptr)
    
    ierror = args%len(num_args)
    
    if (num_args /= 1) then
       call raise_exception(TypeError, "xrite_crystal_Cell expects exactly 1 argument")
       call args%destroy
       return
    endif

    call get_cell_from_arg(args, cell_p)

    call Write_Crystal_Cell(Cell_p%p)

    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function crystal_metrics_write_crystal_cell
  
end module API_Crystal_Metrics
