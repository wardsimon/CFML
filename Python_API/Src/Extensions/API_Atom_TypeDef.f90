module API_Atom_TypeDef
  
  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use CFML_Atom_TypeDef, only: &
       Atom_List_type, &
       Atoms_Cell_type, &
       Allocate_Atom_List, &
       Atoms_Cell_To_List, &
       Write_Atom_list, &
       Atom_type
  

  implicit none

  type Atom_type_p
     type(Atom_type), pointer :: p
  end type Atom_type_p
  
  type Atom_list_type_p
     type(Atom_list_type), pointer :: p
  end type Atom_list_type_p

  type Atoms_cell_type_p
     type(Atoms_cell_type), pointer :: p
  end type Atoms_cell_type_p

contains

  subroutine get_atom_list_type_from_arg(args, alist_p, indx)
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

  end subroutine get_atom_list_type_from_arg

  subroutine get_atom_type_from_arg(args, atom_p, indx)
    type(tuple)                            :: args
    type(Atom_type_p), intent(out)         :: atom_p
    integer, optional                      :: indx
    
    type(object) :: arg_obj
    type(list) :: arg_list
    integer :: atom_p12(12)
  
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
       ierror = cast(atom_p12(ii), t)
    enddo
    atom_p = transfer(atom_p12, atom_p)

  end subroutine get_atom_type_from_arg

!!$  ! @brief Allocate an atom list for N atoms
!!$  ! @todo Do we need to bind that function ??
!!$  function atom_typedef_set_atomlist(self_ptr, args_ptr) result(r) bind(c)
!!$
!!$    type(c_ptr), value :: self_ptr
!!$    type(c_ptr), value :: args_ptr
!!$    type(c_ptr) :: r
!!$    
!!$    type(tuple) :: args
!!$    type(dict)  :: retval
!!$    integer     :: num_args
!!$    integer     :: ierror
!!$    integer     :: i,ii
!!$    type(list)  :: index_obj
!!$
!!$    type(object)                         :: natom_obj
!!$    integer                              :: natom
!!$    type(Atom_list_type_p)               :: alist_p
!!$    integer                              :: alist_p12(12)
!!$       
!!$    r = C_NULL_PTR      
!!$    call unsafe_cast_from_c_ptr(args, args_ptr)
!!$    
!!$    ierror = args%len(num_args)
!!$    
!!$    if (num_args /=1 ) then
!!$       call raise_exception(TypeError, "set_atomlist expects exactly 1 argument")
!!$       call args%destroy
!!$       return
!!$    endif
!!$     
!!$    ierror = args%getitem(natom_obj, 0)        
!!$    ierror = cast_nonstrict(natom, natom_obj)
!!$    
!!$       
!!$    allocate(alist_p%p)
!!$    call Allocate_Atom_List(natom, alist_p%p) 
!!$    
!!$    alist_p12 = transfer(alist_p,alist_p12)
!!$    ierror = list_create(index_obj)
!!$    do ii=1,12
!!$       ierror = index_obj%append(alist_p12(ii))
!!$    end do
!!$
!!$    ierror = dict_create(retval)
!!$    ierror = retval%setitem("address", index_obj)
!!$    r = retval%get_c_ptr()
!!$
!!$  end function atom_typedef_set_atomlist

  function atom_typedef_del_atom_list(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror

    type(Atom_list_type_p) :: a_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "del_atom_list expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_atom_list_type_from_arg(args, a_p)

    if(allocated(a_p%p%atom)) deallocate(a_p%p%atom)
    deallocate(a_p%p)
    !
    ierror = dict_create(retval)
    r = retval%get_c_ptr()

    call args%destroy

  end function atom_typedef_del_atom_list

  function atom_typedef_get_item(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    
    type(tuple) :: args
    type(dict)  :: retval
    integer     :: num_args
    integer     :: ierror
    integer     :: ii

    type(Atom_list_type_p)               :: alist_p
    type(Atom_type_p)                    :: a_p

    type(object) :: item_obj
    integer      :: item, a_p12(12)
    type(list)   :: a_obj

    r = C_NULL_PTR      
    call unsafe_cast_from_c_ptr(args, args_ptr)
    
    ierror = args%len(num_args)
    
    if (num_args /= 2) then
       call raise_exception(TypeError, "get_item expects exactly 2 arguments")
       call args%destroy
       return
    endif
    
    call get_atom_list_type_from_arg(args, alist_p, 0)

    ierror = args%getitem(item_obj, 1)
    ierror = cast_nonstrict(item, item_obj)

    allocate(a_p%p)

    a_p%p = alist_p%p%atom(item)

    a_p12    = transfer(a_p, a_p12)
    ierror = list_create(a_obj)
    do ii=1,12
       ierror = a_obj%append(a_p12(ii))
    end do

    ierror = dict_create(retval)
    ierror = retval%setitem("Atom", a_obj)

    r = retval%get_c_ptr()
    call args%destroy

  end function atom_typedef_get_item

  function atom_typedef_get_natoms(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_list_type_p) :: atom_list_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_natoms expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_list_type_from_arg(args, atom_list_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("natoms", atom_list_type_pointer%p%natoms)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_natoms
  
  

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
    
    call get_atom_list_type_from_arg(args, alist_p)

    call Write_Atom_List(alist_p%p)

    ierror = dict_create(retval)
    r = retval%get_c_ptr()
    call args%destroy

  end function atom_typedef_write_atom_list



  !----------------------------------------------------------------------------------
  ! type(Atom_type)
  !----------------------------------------------------------------------------------
  function atom_typedef_del_atom(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror

    type(Atom_type_p) :: a_p

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "del_atom expects exactly 1 argument")
       call args%destroy
       return
    endif

    !
    call get_atom_type_from_arg(args, a_p)

    deallocate(a_p%p)
    !
    ierror = dict_create(retval)
    r = retval%get_c_ptr()

    call args%destroy

  end function atom_typedef_del_atom
    
  function atom_typedef_get_Lab(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Lab expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Lab", trim(atom_type_pointer%p%Lab))
    
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Lab
  
  
  function atom_typedef_get_ChemSymb(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_ChemSymb expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("ChemSymb", trim(atom_type_pointer%p%ChemSymb))
    
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_ChemSymb
  
  
  function atom_typedef_get_SfacSymb(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_SfacSymb expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("SfacSymb", trim(atom_type_pointer%p%SfacSymb))
    
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_SfacSymb
  
  
  function atom_typedef_get_wyck(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_wyck expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("wyck", trim(atom_type_pointer%p%wyck))
    
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_wyck
  
  
  function atom_typedef_get_Active(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Active expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Active", atom_type_pointer%p%Active)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Active
  
  
  function atom_typedef_get_Z(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Z expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Z", atom_type_pointer%p%Z)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Z
  
  
  function atom_typedef_get_Mult(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Mult expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Mult", atom_type_pointer%p%Mult)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Mult
  
  
  function atom_typedef_get_X(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: X
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_X expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(X, atom_type_pointer%p%X)
    ierror = dict_create(retval)
    ierror = retval%setitem("X", X)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_X
  
  
  function atom_typedef_get_X_Std(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: X_Std
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_X_Std expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(X_Std, atom_type_pointer%p%X_Std)
    ierror = dict_create(retval)
    ierror = retval%setitem("X_Std", X_Std)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_X_Std
  
  
  function atom_typedef_get_MX(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: MX
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_MX expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(MX, atom_type_pointer%p%MX)
    ierror = dict_create(retval)
    ierror = retval%setitem("MX", MX)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_MX
  
  
  function atom_typedef_get_LX(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: LX
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_LX expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(LX, atom_type_pointer%p%LX)
    ierror = dict_create(retval)
    ierror = retval%setitem("LX", LX)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_LX
  
  
  function atom_typedef_get_Occ(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Occ expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Occ", atom_type_pointer%p%Occ)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Occ
  
  
  function atom_typedef_get_Occ_Std(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Occ_Std expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Occ_Std", atom_type_pointer%p%Occ_Std)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Occ_Std
  
  
  function atom_typedef_get_MOcc(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_MOcc expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("MOcc", atom_type_pointer%p%MOcc)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_MOcc
  
  
  function atom_typedef_get_LOcc(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_LOcc expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("LOcc", atom_type_pointer%p%LOcc)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_LOcc
  
  
  function atom_typedef_get_Biso(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Biso expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Biso", atom_type_pointer%p%Biso)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Biso
  
  
  function atom_typedef_get_Biso_std(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Biso_std expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Biso_std", atom_type_pointer%p%Biso_std)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Biso_std
  
  
  function atom_typedef_get_MBiso(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_MBiso expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("MBiso", atom_type_pointer%p%MBiso)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_MBiso
  
  
  function atom_typedef_get_LBiso(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_LBiso expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("LBiso", atom_type_pointer%p%LBiso)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_LBiso
  
  
  function atom_typedef_get_Utype(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Utype expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Utype", trim(atom_type_pointer%p%Utype))
    
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Utype
  
  
  function atom_typedef_get_ThType(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_ThType expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("ThType", trim(atom_type_pointer%p%ThType))
    
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_ThType
  
  
  function atom_typedef_get_U(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: U
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_U expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(U, atom_type_pointer%p%U)
    ierror = dict_create(retval)
    ierror = retval%setitem("U", U)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_U
  
  
  function atom_typedef_get_U_std(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: U_std
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_U_std expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(U_std, atom_type_pointer%p%U_std)
    ierror = dict_create(retval)
    ierror = retval%setitem("U_std", U_std)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_U_std
  
  
  function atom_typedef_get_MU(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: MU
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_MU expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(MU, atom_type_pointer%p%MU)
    ierror = dict_create(retval)
    ierror = retval%setitem("MU", MU)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_MU
  
  
  function atom_typedef_get_LU(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: LU
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_LU expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(LU, atom_type_pointer%p%LU)
    ierror = dict_create(retval)
    ierror = retval%setitem("LU", LU)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_LU
  
  
  function atom_typedef_get_Ueq(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Ueq expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Ueq", atom_type_pointer%p%Ueq)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Ueq
  
  
  function atom_typedef_get_Charge(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Charge expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Charge", atom_type_pointer%p%Charge)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Charge
  
  
  function atom_typedef_get_Moment(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Moment expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Moment", atom_type_pointer%p%Moment)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Moment
  
  
  function atom_typedef_get_Ind(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: Ind
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Ind expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(Ind, atom_type_pointer%p%Ind)
    ierror = dict_create(retval)
    ierror = retval%setitem("Ind", Ind)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Ind
  
  
  function atom_typedef_get_NVar(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_NVar expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("NVar", atom_type_pointer%p%NVar)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_NVar
  
  
  function atom_typedef_get_VarF(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: VarF
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_VarF expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(VarF, atom_type_pointer%p%VarF)
    ierror = dict_create(retval)
    ierror = retval%setitem("VarF", VarF)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_VarF
  
  
  function atom_typedef_get_MVarF(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: MVarF
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_MVarF expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(MVarF, atom_type_pointer%p%MVarF)
    ierror = dict_create(retval)
    ierror = retval%setitem("MVarF", MVarF)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_MVarF
  
  
  function atom_typedef_get_LVarF(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: LVarF
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_LVarF expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(LVarF, atom_type_pointer%p%LVarF)
    ierror = dict_create(retval)
    ierror = retval%setitem("LVarF", LVarF)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_LVarF
  
  
  function atom_typedef_get_AtmInfo(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_AtmInfo expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("AtmInfo", trim(atom_type_pointer%p%AtmInfo))
    
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_AtmInfo
  
  
  function atom_typedef_get_m_xyz(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: m_xyz
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_m_xyz expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(m_xyz, atom_type_pointer%p%m_xyz)
    ierror = dict_create(retval)
    ierror = retval%setitem("m_xyz", m_xyz)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_m_xyz
  
  
  function atom_typedef_get_sm_xyz(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: sm_xyz
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_sm_xyz expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(sm_xyz, atom_type_pointer%p%sm_xyz)
    ierror = dict_create(retval)
    ierror = retval%setitem("sm_xyz", sm_xyz)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_sm_xyz
  
  
  function atom_typedef_get_Mm_xyz(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: Mm_xyz
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Mm_xyz expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(Mm_xyz, atom_type_pointer%p%Mm_xyz)
    ierror = dict_create(retval)
    ierror = retval%setitem("Mm_xyz", Mm_xyz)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Mm_xyz
  
  
  function atom_typedef_get_Lm_xyz(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Atom_Type_p) :: atom_type_pointer
  
    type(ndarray) :: Lm_xyz
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Lm_xyz expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_atom_type_from_arg(args, atom_type_pointer)
    ierror = ndarray_create(Lm_xyz, atom_type_pointer%p%Lm_xyz)
    ierror = dict_create(retval)
    ierror = retval%setitem("Lm_xyz", Lm_xyz)
    r = retval%get_c_ptr()
    call args%destroy
    
  end function atom_typedef_get_Lm_xyz
end module API_Atom_TypeDef
