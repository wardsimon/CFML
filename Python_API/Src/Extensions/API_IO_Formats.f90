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
  
  use CFML_IO_Formats,                only: &
       Job_Info_type, &
       Get_Job_Info, &
       Readn_set_Xtal_structure
  
  use API_Atom_TypeDef,               only: Atom_List_Type_p  
  use API_Crystallographic_Symmetry,  only: Space_Group_Type_p
  use API_Crystal_Metrics,            only: Crystal_Cell_Type_p

  implicit none

  type Job_Info_type_p
     type(Job_Info_type), pointer :: p
  end type Job_Info_type_p

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

    type(object)                    :: filename_obj, mode_obj
    character(len=:), allocatable   :: filename, mode

    type(Crystal_Cell_type_p)       :: cell_p
    type(Space_Group_type_p)        :: spg_p
    type(Atom_list_type_p)          :: a_p
    type(Job_info_type_p)           :: job_p

    integer                         :: cell_p12(12), spg_p12(12), a_p12(12), job_p12(12)
    type(list)                      :: cell_obj, spg_obj, a_obj, job_obj

    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)

    if (num_args /= 2) then
       call raise_exception(TypeError, "readn_set_xtal_structure expects exactly 2 arguments: filename, mode")
       call args%destroy
       return
    endif

    ierror = args%getitem(filename_obj, 0)
    ierror = cast(filename, filename_obj)

    ierror = args%getitem(mode_obj, 1)
    ierror = cast_nonstrict(mode, mode_obj)

    allocate(cell_p%p, spg_p%p, a_p%p, job_p%p)
    call Readn_set_Xtal_structure(trim(filename),cell_p%p,spg_p%p,a_p%p,Mode=mode,Job_info=job_p%p)

    cell_p12 = transfer(cell_p, cell_p12)
    ierror = list_create(cell_obj)
    do ii=1,12
       ierror = cell_obj%append(cell_p12(ii))
    end do

    spg_p12   = transfer(spg_p, spg_p12)
    ierror = list_create(spg_obj)
    do ii=1,12
       ierror = spg_obj%append(spg_p12(ii))
    end do

    a_p12    = transfer(a_p, a_p12)
    ierror = list_create(a_obj)
    do ii=1,12
       ierror = a_obj%append(a_p12(ii))
    end do

    job_p12    = transfer(job_p, job_p12)
    ierror = list_create(job_obj)
    do ii=1,12
       ierror = job_obj%append(job_p12(ii))
    end do

    ierror = dict_create(retval)
    ierror = retval%setitem("Cell", cell_obj)
    ierror = retval%setitem("SpG",  spg_obj)
    ierror = retval%setitem("A",    a_obj)
    ierror = retval%setitem("JobInfo", job_obj)

    r = retval%get_c_ptr()
    call args%destroy

  end function IO_formats_readn_set_xtal_structure

  ! @brief Create an atom list from an array of CIF lines
  function IO_Formats_jobinfo_from_CIF_string_array(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r

    type(tuple) :: args
    type(dict)  :: retval
    integer     :: num_args
    integer     :: ierror

    integer     :: max_line_length
    integer     :: ii
    type(list)  :: job_obj
    type(object):: item_obj
    character(len=:), allocatable :: item

    integer                              :: nline, n_ini
    type(Job_Info_type_p)                :: job_p
    integer                              :: job_p12(12)

    character(len=:), dimension(:), allocatable :: stringarray
    type(object)                            :: stringarray_obj
    type(list)                           :: stringarray_list
    
    r = C_NULL_PTR      
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ierror = args%len(num_args)
    if (num_args /=1 ) then
       call raise_exception(TypeError, "jobinfo_from_CIF_string_array expects exactly 1 argument")
       call args%destroy
       return
    endif

    ierror = args%getitem(stringarray_obj, 0)
    ierror = cast(stringarray_list, stringarray_obj)
    ierror = stringarray_list%len(nline)

    ! Compute max line length
    max_line_length = 0
    do ii=1,nline
        ierror = stringarray_list%getitem(item_obj, ii-1)
        ierror = cast(item, item_obj)
        if (len(item) > max_line_length) then
            max_line_length = len(item)
        endif
        call item_obj%destroy
    enddo
    allocate(character(max_line_length) :: stringarray(nline))

    ! Fill stringarray
    do ii=1,nline
        ierror = stringarray_list%getitem(item_obj, ii-1)
        ierror = cast(item, item_obj)
        stringarray(ii) = item
        call item_obj%destroy
    enddo

    ! Create atom list from stringarray
    allocate(job_p%p)
    n_ini = 1
    call Get_Job_Info(stringarray, n_ini, nline+1, job_p%p)

    !
    job_p12 = transfer(job_p,job_p12)
    ierror = list_create(job_obj)
    do ii=1,12
       ierror = job_obj%append(job_p12(ii))
    end do

    !
    ierror = dict_create(retval)
    ierror = retval%setitem("JobInfo", job_obj)
    r = retval%get_c_ptr()

    !
    call args%destroy
    call job_obj%destroy
    deallocate(stringarray)
    call stringarray_obj%destroy
    call stringarray_list%destroy

  end function IO_Formats_jobinfo_from_CIF_string_array

  function IO_Formats_del_jobinfo(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    
    
  end function IO_Formats_del_jobinfo

  function IO_Formats_get_title(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_title expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("title", trim(job_info_type_pointer%p%title))
    
    r = retval%get_c_ptr()
    call args%destroy
     
  end function IO_Formats_get_title
  
  
  function IO_Formats_get_num_phases(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_num_phases expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("num_phases", job_info_type_pointer%p%num_phases)
    r = retval%get_c_ptr()
    call args%destroy
     
  end function IO_Formats_get_num_phases
  
  
  function IO_Formats_get_num_patterns(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_num_patterns expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("num_patterns", job_info_type_pointer%p%num_patterns)
    r = retval%get_c_ptr()
    call args%destroy
     
  end function IO_Formats_get_num_patterns
  
  
  function IO_Formats_get_num_cmd(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_num_cmd expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("num_cmd", job_info_type_pointer%p%num_cmd)
    r = retval%get_c_ptr()
    call args%destroy
     
  end function IO_Formats_get_num_cmd
  
  
  function IO_Formats_get_pat_typ(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: pat_typ
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_pat_typ expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(pat_typ, job_info_type_pointer%p%patt_typ)
    ierror = dict_create(retval)
    ierror = retval%setitem("pat_typ", pat_typ)
    r = retval%get_c_ptr()
    call args%destroy
    call pat_typ%destroy
     
  end function IO_Formats_get_pat_typ
  
  
  function IO_Formats_get_phas_nam(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: phas_nam
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_phas_nam expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(phas_nam, job_info_type_pointer%p%phas_nam)
    ierror = dict_create(retval)
    ierror = retval%setitem("phas_nam", phas_nam)
    r = retval%get_c_ptr()
    call args%destroy
    call phas_nam%destroy
     
  end function IO_Formats_get_phas_nam
  
  
  function IO_Formats_get_cmd(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: cmd
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_cmd expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(cmd, job_info_type_pointer%p%cmd)
    ierror = dict_create(retval)
    ierror = retval%setitem("cmd", cmd)
    r = retval%get_c_ptr()
    call args%destroy
    call cmd%destroy
     
  end function IO_Formats_get_cmd
  
  
  function IO_Formats_get_range_stl(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: range_stl
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_range_stl expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(range_stl, job_info_type_pointer%p%range_stl)
    ierror = dict_create(retval)
    ierror = retval%setitem("range_stl", range_stl)
    r = retval%get_c_ptr()
    call args%destroy
    call range_stl%destroy
     
  end function IO_Formats_get_range_stl
  
  
  function IO_Formats_get_range_q(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: range_q
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_range_q expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(range_q, job_info_type_pointer%p%range_q)
    ierror = dict_create(retval)
    ierror = retval%setitem("range_q", range_q)
    r = retval%get_c_ptr()
    call args%destroy
    call range_q%destroy
     
  end function IO_Formats_get_range_q
  
  
  function IO_Formats_get_range_d(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: range_d
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_range_d expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(range_d, job_info_type_pointer%p%range_d)
    ierror = dict_create(retval)
    ierror = retval%setitem("range_d", range_d)
    r = retval%get_c_ptr()
    call args%destroy
    call range_d%destroy
     
  end function IO_Formats_get_range_d
  
  
  function IO_Formats_get_range_2theta(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: range_2theta
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_range_2theta expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(range_2theta, job_info_type_pointer%p%range_2theta)
    ierror = dict_create(retval)
    ierror = retval%setitem("range_2theta", range_2theta)
    r = retval%get_c_ptr()
    call args%destroy
    call range_2theta%destroy
     
  end function IO_Formats_get_range_2theta
  
  
  function IO_Formats_get_range_energy(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: range_energy
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_range_energy expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(range_energy, job_info_type_pointer%p%range_energy)
    ierror = dict_create(retval)
    ierror = retval%setitem("range_energy", range_energy)
    r = retval%get_c_ptr()
    call args%destroy
    call range_energy%destroy
     
  end function IO_Formats_get_range_energy
  
  
  function IO_Formats_get_lambda(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    type(ndarray) :: lambda
    
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_lambda expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = ndarray_create(lambda, job_info_type_pointer%p%lambda)
    ierror = dict_create(retval)
    ierror = retval%setitem("lambda", lambda)
    r = retval%get_c_ptr()
    call args%destroy
    call lambda%destroy
     
  end function IO_Formats_get_lambda
  
  
  function IO_Formats_get_ratio(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_ratio expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("ratio", job_info_type_pointer%p%ratio)
    r = retval%get_c_ptr()
    call args%destroy
     
  end function IO_Formats_get_ratio
  
  
  function IO_Formats_get_dtt1(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_dtt1 expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("dtt1", job_info_type_pointer%p%dtt1)
    r = retval%get_c_ptr()
    call args%destroy
     
  end function IO_Formats_get_dtt1
  
  
  function IO_Formats_get_dtt2(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(Job_info_Type_p) :: job_info_type_pointer
  
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_dtt2 expects exactly 1 argument")
       call args%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_job_info_type_from_arg(args, job_info_type_pointer)
    ierror = dict_create(retval)
    ierror = retval%setitem("dtt2", job_info_type_pointer%p%dtt2)
    r = retval%get_c_ptr()
    call args%destroy
     
  end function IO_Formats_get_dtt2
end module API_IO_Formats
