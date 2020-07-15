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
       Write_Crystal_Cell, &
       Set_Crystal_Cell 
  
  implicit none

  !type definitions 
  type Crystal_Cell_type_p
     type(Crystal_Cell_type), pointer :: p
  end type Crystal_Cell_type_p
  
contains 

  subroutine get_cell_from_arg(args, cell_p)
    type(tuple)                            :: args
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

  end subroutine get_cell_from_arg

  
  ! @brief Create a crystal cell from a set of lattice parameters
  function crystal_metrics_set_crystal_cell(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    
    type(tuple) :: args
    type(dict)  :: retval
    integer     :: num_args
    integer     :: ierror
    integer     :: i,ii
    type(list)  :: index_obj

    type(object)                         :: cellv_obj, angl_obj
    type(ndarray)                        :: cellv_nd, angl_nd
    real(kind=cp), dimension(:), pointer :: cellv_p, angl_p
    type(Crystal_Cell_type_p)            :: cell_p
    integer                              :: cell_p12(12)
       
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
    ierror = cast(cellv_nd, cellv_obj)
    ierror = cellv_nd%get_data(cellv_p)
        
    ierror = args%getitem(angl_obj, 1)         !angl   In -> angles of cell parameters
    ierror = cast(angl_nd, angl_obj)
    ierror = angl_nd%get_data(angl_p)
       
    allocate(cell_p%p)
    call Set_Crystal_Cell(cellv_p, angl_p, cell_p%p)
    
    cell_p12 = transfer(cell_p,cell_p12)
    ierror = list_create(index_obj)
    do ii=1,12
       ierror = index_obj%append(cell_p12(ii))
    end do
    !deallocate(cell_p%p)

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
       call raise_exception(TypeError, "write_crystal_Cell expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    call get_cell_from_arg(args, cell_p)

    call Write_Crystal_Cell(Cell_p%p)

    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function crystal_metrics_write_crystal_cell

  
  ! @brief Get the parameter 'cell' from a Crystal_Cell_Type object
  function crystallographic_symmetry_get_cell(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Crystal_Cell_type_p)   :: cell_p
  
    type(ndarray) :: cell
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_cell expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    call get_cell_from_arg(args, cell_p)
    ierror = ndarray_create(cell, cell_p%p%cell)
    ierror = dict_create(retval)
    ierror = retval%setitem("cell", cell)
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_cell


  ! @brief Get the parameter 'cell' from a Crystal_Cell_Type object
  function crystallographic_symmetry_get_angl(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Crystal_Cell_type_p)   :: cell_p
    
    type(ndarray) :: angl
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_angl expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    call get_cell_from_arg(args, cell_p)
    ierror = ndarray_create(angl, cell_p%p%angl)
    ierror = dict_create(retval)
    ierror = retval%setitem("angl", angl)
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_angl
  
end module API_Crystal_Metrics
