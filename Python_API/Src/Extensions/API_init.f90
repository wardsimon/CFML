! **************************************************************************
!
! CrysFML API
!
! @file      Src/Extensions/API_init.f90
! @brief     Initialisation functions for Fortran binding
!
! @homepage  https://code.ill.fr/scientific-software/crysfml
! @license   GNU LGPL (see LICENSE)
! @copyright Institut Laue Langevin 2020-now
! @authors   Scientific Computing Group at ILL (see AUTHORS), based on Elias Rabel work for Forpy
!
! **************************************************************************

module API_init
  use forpy_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env 
  use API_Crystallographic_Symmetry, only: crysfml_symmetry_create_space_group, crysfml_symmetry_get_description, crysfml_symmetry_get_latt_trans
  use API_IO_Formats, only: crysfml_IO_readn_set_xtal_structure
  implicit none

  PRIVATE
  type(PythonModule), save      :: mod_def
  type(PythonMethodTable), save :: method_table

CONTAINS
  ! Initialisation function for Python 3
  function PyInit_crysfml_API() bind(c, name="PyInit_crysfml_API") result(m)
  !DEC$ ATTRIBUTES DLLEXPORT :: PyInit_crysfml_API
    type(c_ptr) :: m
    m = init()
  end function PyInit_crysfml_API

  ! Initialisation function for Python 2
  subroutine init_crysfml_API() bind(c, name="init_crysfml_API")
    type(c_ptr) :: m
    m = init()
  end subroutine init_crysfml_API

  ! Initialisation function
  function init() result(m)
    type(c_ptr) :: m
    integer :: ierror
    ierror = forpy_initialize()
    call method_table%init(4)
    call method_table%add_method("create_space_group", &                  ! method name
         "Creates the space group", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crysfml_symmetry_create_space_group))  ! address of Fortran function to add
    call method_table%add_method("get_description", &                  ! method name
         "Return space group description", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crysfml_symmetry_get_description))  ! address of Fortran function to add
    call method_table%add_method("get_latt_trans", &                  ! method name
         "Get lattice transformation from space group", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crysfml_symmetry_get_latt_trans))  ! address of Fortran function to add
    call method_table%add_method("Readn_set_Xtal_Structure", &                  ! method name
         "read an input file and construct the crystal structure in terms of Cell, SpG and A", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crysfml_IO_readn_set_xtal_structure))  ! address of Fortran function to add
    m = mod_def%init("crysfml_symmetry", "A Python extension for crysFML symmetry", method_table)
  end function init


end module API_init
