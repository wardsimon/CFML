module API_Atom_TypeDef
  
  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use CFML_Atom_TypeDef, only: &
       Atom_List_type, &
       Atoms_Cell_type, &
       Allocate_Atom_List, &
       Atoms_Cell_To_List, &
       Write_Atom_list
  

  implicit none

  type Atom_list_type_p
     type(Atom_list_type), pointer :: p
  end type Atom_list_type_p

  type Atoms_cell_type_p
     type(Atoms_cell_type), pointer :: p
  end type Atoms_cell_type_p

contains

  subroutine get_atoms_list_from_arg(args, alist_p, indx)
    type(tuple)                            :: args
    type(Atom_list_type_p), intent(out)    :: alist_p
    integer, optional                      :: indx
    
    type(object) :: arg_obj
    type(list) :: arg_list
    integer :: alist_p12(12)
  
    integer :: ierror
    integer :: ii
    type(object) :: t

    if (present(indx)) then
       ierror = args%getitem(arg_obj, indx)
    else
       ierror = args%getitem(arg_obj, 0)
    endif
    ierror = cast(arg_list, arg_obj)
    do ii=1,12
       ierror = arg_list%getitem(t, ii-1)
       ierror = cast(alist_p12(ii), t)
    enddo
    alist_p = transfer(alist_p12, alist_p)

  end subroutine get_atoms_list_from_arg

  ! @brief Allocate an atom list for N atoms
  ! @todo Do we need to bind that function ??
  function atom_typedef_set_atomlist(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    
    type(tuple) :: args
    type(dict)  :: retval
    integer     :: num_args
    integer     :: ierror
    integer     :: i,ii
    type(list)  :: index_obj

    type(object)                         :: natom_obj
    integer                              :: natom
    type(Atom_list_type_p)               :: alist_p
    integer                              :: alist_p12(12)
       
    r = C_NULL_PTR      
    call unsafe_cast_from_c_ptr(args, args_ptr)
    
    ierror = args%len(num_args)
    
    if (num_args /=1 ) then
       call raise_exception(TypeError, "set_atomlist expects exactly 1 argument")
       call args%destroy
       return
    endif
     
    ierror = args%getitem(natom_obj, 0)        
    ierror = cast_nonstrict(natom, natom_obj)
    
       
    allocate(alist_p%p)
    call Allocate_Atom_List(natom, alist_p%p) 
    
    alist_p12 = transfer(alist_p,alist_p12)
    ierror = list_create(index_obj)
    do ii=1,12
       ierror = index_obj%append(alist_p12(ii))
    end do

    ierror = dict_create(retval)
    ierror = retval%setitem("address", index_obj)
    r = retval%get_c_ptr()

  end function atom_typedef_set_atomlist

  ! @brif Print the description of the atom list to standard output
  ! @todo Optional arguments
  function atom_typedef_write_atom_list(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    
    type(tuple) :: args
    type(dict)  :: retval
    integer     :: num_args
    integer     :: ierror

    type(Atom_list_type_p)               :: alist_p

    r = C_NULL_PTR      
    call unsafe_cast_from_c_ptr(args, args_ptr)
    
    ierror = args%len(num_args)
    
    if (num_args /= 1) then
       call raise_exception(TypeError, "write_atom_list expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    call get_atoms_list_from_arg(args, alist_p)

    call Write_Atom_List(alist_p%p)

    ierror = dict_create(retval)
    r = retval%get_c_ptr()

  end function atom_typedef_write_atom_list

end module API_Atom_TypeDef
