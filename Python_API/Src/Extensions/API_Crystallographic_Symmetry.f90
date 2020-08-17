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

  subroutine get_space_group_type_from_arg(args, spgr_p, indx)
    type(tuple)                            :: args
    type(Space_Group_type_p), intent(out)  :: spgr_p
    integer, optional                      :: indx
    
    type(object) :: arg_obj
    type(list) :: arg_list
    integer :: spgr_p12(12)
    
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
       ierror = cast(spgr_p12(ii), t)
       call t%destroy
    enddo
    spgr_p = transfer(spgr_p12, spgr_p)

    call arg_obj%destroy
    call arg_list%destroy

  end subroutine get_space_group_type_from_arg

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

    !
    ierror = list_create(index_obj)
    do ii=1,12
       ierror = index_obj%append(spgr_p12(ii))
    end do

    ierror = dict_create(retval)
    ierror = retval%setitem("address", index_obj)
    r = retval%get_c_ptr()

    deallocate(spgr)
    call args%destroy
    call index_obj%destroy


  end function crystallographic_symmetry_set_spacegroup

  function crystallographic_symmetry_del_spacegroup(self_ptr, args_ptr) result(r) bind(c)

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
    call get_space_group_type_from_arg(args, spgr_p)

    if(allocated(spgr_p%p%Latt_trans)) deallocate(spgr_p%p%Latt_trans)
    if(allocated(spgr_p%p%SymOp)) deallocate(spgr_p%p%SymOp)
    if(allocated(spgr_p%p%SymOpSymb)) deallocate(spgr_p%p%SymOpSymb)
    deallocate(spgr_p%p)
    !
    ierror = dict_create(retval)
    r = retval%get_c_ptr()

    call args%destroy

  end function crystallographic_symmetry_del_spacegroup

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
    call get_space_group_type_from_arg(args, spgr_p)

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
    type(Space_Group_Type_p) :: space_group_type_pointer

    type(ndarray) :: latt_trans

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

    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = ndarray_create(latt_trans, space_group_type_pointer%p%latt_trans)
    ierror = dict_create(retval)
    ierror = retval%setitem("latt_trans", latt_trans)
    r = retval%get_c_ptr()

  end function crystallographic_symmetry_get_latt_trans


  function crystallographic_symmetry_get_hexa(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_hexa expects exactly 1 argument")
       call args%destroy
       return
    endif

    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("hexa", space_group_type_pointer%p%hexa)
    r = retval%get_c_ptr()

  end function crystallographic_symmetry_get_hexa


  function crystallographic_symmetry_get_numspg(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_numspg expects exactly 1 argument")
       call args%destroy
       return
    endif

    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("numspg", space_group_type_pointer%p%numspg)
    r = retval%get_c_ptr()

  end function crystallographic_symmetry_get_numspg


  function crystallographic_symmetry_get_spg_symb(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_spg_symb expects exactly 1 argument")
       call args%destroy
       return
    endif

    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("spg_symb", trim(space_group_type_pointer%p%spg_symb))

    r = retval%get_c_ptr()

  end function crystallographic_symmetry_get_spg_symb

   function crystallographic_symmetry_get_hall(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_hall expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("hall", trim(space_group_type_pointer%p%hall))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_hall
  
  
  function crystallographic_symmetry_get_ghall(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_ghall expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("ghall", trim(space_group_type_pointer%p%ghall))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_ghall
  
  
  function crystallographic_symmetry_get_CrystalSys(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_CrystalSys expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("CrystalSys", trim(space_group_type_pointer%p%CrystalSys))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_CrystalSys
  
  
  function crystallographic_symmetry_get_Laue(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Laue expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Laue", trim(space_group_type_pointer%p%Laue))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_Laue
  
  
  function crystallographic_symmetry_get_PG(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_PG expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("PG", trim(space_group_type_pointer%p%PG))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_PG
  
  
  function crystallographic_symmetry_get_Info(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_Info expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("Info", trim(space_group_type_pointer%p%Info))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_Info
  
  
  function crystallographic_symmetry_get_SG_setting(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_SG_setting expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("SG_setting", trim(space_group_type_pointer%p%SG_setting))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_SG_setting
  
  
  function crystallographic_symmetry_get_SPG_lat(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_SPG_lat expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("SPG_lat", trim(space_group_type_pointer%p%SPG_lat))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_SPG_lat
  
  
  function crystallographic_symmetry_get_SPG_latsy(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_SPG_latsy expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("SPG_latsy", trim(space_group_type_pointer%p%SPG_latsy))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_SPG_latsy
  
  
  function crystallographic_symmetry_get_NumLat(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_NumLat expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("NumLat", space_group_type_pointer%p%NumLat)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_NumLat
  
  
  function crystallographic_symmetry_get_bravais(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_bravais expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("bravais", trim(space_group_type_pointer%p%bravais))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_bravais
  
  
  function crystallographic_symmetry_get_centre(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_centre expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("centre", trim(space_group_type_pointer%p%centre))
    
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_centre
  
  
  function crystallographic_symmetry_get_centred(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_centred expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("centred", space_group_type_pointer%p%centred)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_centred
  
  
  function crystallographic_symmetry_get_centre_coord(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    type(ndarray) :: centre_coord
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_centre_coord expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = ndarray_create(centre_coord, space_group_type_pointer%p%centre_coord)
    ierror = dict_create(retval)
    ierror = retval%setitem("centre_coord", centre_coord)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_centre_coord
  
  
  function crystallographic_symmetry_get_numops(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_numops expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("numops", space_group_type_pointer%p%numops)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_numops
  
  
  function crystallographic_symmetry_get_multip(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_multip expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("multip", space_group_type_pointer%p%multip)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_multip
  
  
  function crystallographic_symmetry_get_num_gen(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_num_gen expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("num_gen", space_group_type_pointer%p%num_gen)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_num_gen
  
  
  function crystallographic_symmetry_get_SymOp(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_SymOp expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    !> @todo type(Sym_Oper_Type), allocatable,dimension(:)
    !ierror = retval%setitem("SymOp", space_group_type_pointer%p%SymOp)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_SymOp
  
  
  function crystallographic_symmetry_get_SymopSymb(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_SymopSymb expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    !> @todo character(len=50),   allocatable,dimension(:)
    !ierror = retval%setitem("SymopSymb", space_group_type_pointer%p%SymopSymb)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_SymopSymb
  
  
  function crystallographic_symmetry_get_wyckoff(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_wyckoff expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = dict_create(retval)
    !> @todo type(Wyckoff_Type)
    !ierror = retval%setitem("wyckoff", space_group_type_pointer%p%wyckoff)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_wyckoff
  
  
  function crystallographic_symmetry_get_R_asym_unit(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Space_Group_Type_p) :: space_group_type_pointer
  
    type(ndarray) :: R_asym_unit
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_R_asym_unit expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_space_group_type_from_arg(args, space_group_type_pointer)
    ierror = ndarray_create(R_asym_unit, space_group_type_pointer%p%R_asym_unit)
    ierror = dict_create(retval)
    ierror = retval%setitem("R_asym_unit", R_asym_unit)
    r = retval%get_c_ptr()
    
  end function crystallographic_symmetry_get_R_asym_unit
  




     ! integer                                       :: NumSpg=0         ! Number of the Space Group
     ! character(len=20)                             :: SPG_Symb=" "     ! Hermann-Mauguin Symbol
     ! character(len=16)                             :: Hall=" "         ! Hall symbol
     ! character(len=90)                             :: gHall=" "        ! Generalised Hall symbol
     ! character(len=12)                             :: CrystalSys=" "   ! Crystal system
     ! character(len= 5)                             :: Laue=" "         ! Laue Class
     ! character(len= 5)                             :: PG=" "           ! Point group
     ! character(len= 5)                             :: Info=" "         ! Extra information
     ! character(len=90)                             :: SG_setting=" "   ! Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
     ! logical                                       :: Hexa=.false.     !
     ! character(len= 1)                             :: SPG_lat=" "      ! Lattice type
     ! character(len= 2)                             :: SPG_latsy=" "    ! Lattice type Symbol
     ! integer                                       :: NumLat=0         ! Number of lattice points in a cell
     ! real(kind=cp), allocatable,dimension(:,:)     :: Latt_trans       ! Lattice translations (3,12)
     ! character(len=51)                             :: Bravais=" "      ! String with Bravais symbol + translations
     ! character(len=80)                             :: Centre=" "       ! Alphanumeric information about the center of symmetry
     ! integer                                       :: Centred=0        ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
     ! real(kind=cp), dimension(3)                   :: Centre_coord=0.0 ! Fractional coordinates of the inversion centre
     ! integer                                       :: NumOps=0         ! Number of reduced set of S.O.
     ! integer                                       :: Multip=0         ! Multiplicity of the general position
     ! integer                                       :: Num_gen          ! Minimum number of operators to generate the Group
     ! type(Sym_Oper_Type), allocatable,dimension(:) :: SymOp            ! Symmetry operators (192)
     ! character(len=50),   allocatable,dimension(:) :: SymopSymb        ! Strings form of symmetry operators
     ! type(Wyckoff_Type)                            :: Wyckoff          ! Wyckoff Information
     ! real(kind=cp),dimension(3,2)                  :: R_Asym_Unit=0.0  ! Asymmetric unit in real space   A


end module API_Crystallographic_Symmetry
