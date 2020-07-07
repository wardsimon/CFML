! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/form_CIF.f90
! @brief     CFML IO Formats Fortran binding
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_IO_Formats !Proposed naming, to be defined
  
  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env 

  use CFML_Crystal_Metrics,           only: Crystal_Cell_Type
  use CFML_IO_Formats,                only: Readn_set_Xtal_structure
  use CFML_Atom_TypeDef,              only: Atom_List_Type

  implicit none

  !type definitions This should be in the API_Crystal_metrics
  type Crystal_Cell_type_p
     type(Crystal_Cell_type), pointer :: p
  end type Crystal_Cell_type_p

  !type definitions This should be in the API_Atom_TypeDef
  type Atom_list_type_p
     type(Atom_list_type), pointer :: p
  end type Atom_list_type_p

  !I also need type Space_Group_Type_p with is defined in symmetry.f90
  !use mymath ?? need proper module name

  private
  type(PythonModule), save      :: mod_def
  type(PythonMethodTable), save :: method_table

contains
  
  ! @brief Initialisation function for Python 3
  function PyInit_crysfml_IO() bind(c, name="PyInit_crysfml_IO") result(m)
    type(c_ptr) :: m
    m = init()
  end function PyInit_crysfml_IO

  ! @brief Initialisation function for Python 2
  subroutine init_crysfml_IO() bind(c, name="init_crysfml_IO")
    type(c_ptr) :: m
    m = init()
  end subroutine init_crysfml_IO

  ! @brief Initialisation function
  function init() result(m)
    
    type(c_ptr) :: m
    integer :: ierror

    ierror = forpy_initialize()

    call method_table%init(1)

    call method_table%add_method("Readn_set_Xtal_Structure", &                  ! method name
         "read an input file and construct the crystal structure in terms of Cell, SpG and A", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crysfml_IO_readn_set_xtal_structure))  ! address of Fortran function to add
    
    m = mod_def%init("crysfml_IO", "A Python extension for crysFML IO formats", method_table)
  end function init

  !-------------------------------------------------------------------------
  ! Implementation of our Python methods
  !-------------------------------------------------------------------------
  ! @brief 
  function crysfml_IO_readn_set_xtal_structure(self_ptr, args_ptr) result(r) bind(c)

    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr)        :: r
    type(tuple)        :: args 
    type(dict)         :: retval
    integer            :: num_args
    integer            :: ierror
    integer            :: ii

    type(object)       :: filename_obj
    character(len=*)   :: filename

    type(Crystal_Cell_type_p) :: cell_p
    type(Space_Group_type_p)  :: spg_p
    type(Atom_list_type_p)    :: a_p

    integer, dimension(12)    :: cell_p12(:), spg_p12(:), a_p(:)
    type(list)                :: cell_obj, spg_obj, a_obj

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
    ierror = cast_to_chars(filename, filename_obj)

    allocate(cell_p%p, spg_p%p, a_p%p)
    call Readn_set_Xtal_structure(trim(filename),cell_p%p,spg_p%p,a_p%p,Mode="CIF")

    cell_p12 = transfer(cell_p, cell_p12)
    ierror = list_create(cell_obj)
    do ii=1,12
       ierror = cell_obj%append(cell_p12(ii))
    end do
    
    spg_12   = transfer(spg_p, spg_p12)
    ierror = list_create(spg_obj)
    do ii=1,12
       ierror = spg_obj%append(spg_p12(ii))
    end do
    
    a_p12    = transfer(a_p, a_p12)
    ierror = list_create(a_obj)
    do ii=1,12
       ierror = a_obj%append(a_p12(ii))
    end do

    ierror = dict_create(retval)
    ierror = retval%setitem("Cell", cell_obj)
    ierror = retval%setitem("SpG",  spg_obj)
    ierror = retval%setitem("A", A_obj)

    r = retval%get_c_ptr()

  end function crysfml_IO_readn_set_xtal_structure
end module API_IO_Formats
