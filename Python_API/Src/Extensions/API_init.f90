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
  
  use API_Crystallographic_Symmetry,only: &
       crystallographic_symmetry_set_spacegroup, &
       crystallographic_symmetry_write_spacegroup, &
       crystallographic_symmetry_get_latt_trans
  
  use API_IO_Formats, only: &
       IO_formats_readn_set_xtal_structure

  use API_Crystal_Metrics, only: &
       crystal_metrics_set_crystal_cell
  
  implicit none

  PRIVATE
  type(PythonModule), save      :: mod_def
  type(PythonMethodTable), save :: method_table

CONTAINS
  ! Initialisation function for Python 3
  function PyInit_crysfml_api() bind(c, name="PyInit_crysfml_api") result(m)
  !DEC$ ATTRIBUTES DLLEXPORT :: PyInit_crysfml_API
    type(c_ptr) :: m
    m = init()
  end function PyInit_crysfml_api

  ! Initialisation function for Python 2
  subroutine init_crysfml_api() bind(c, name="init_crysfml_api")
    type(c_ptr) :: m
    m = init()
  end subroutine init_crysfml_api

  ! Initialisation function
  function init() result(m)
    type(c_ptr) :: m
    integer :: ierror
    ierror = forpy_initialize()
    
    call method_table%init(4)
    
    call method_table%add_method("crystallographic_symmetry_set_spacegroup", &                  ! method name
         "Creates the space group", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_set_spacegroup))  ! address of Fortran function to add
    call method_table%add_method("crystallographic_symmetry_write_spacegroup", &                  ! method name
         "Return space group description", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_write_spacegroup))  ! address of Fortran function to add
    call method_table%add_method("crystallographic_symmetry_get_latt_trans", &                  ! method name
         "Get lattice transformation from space group", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_latt_trans))  ! address of Fortran function to add
    call method_table%add_method("IO_formats_readn_set_xtal_structure", &                  ! method name
         "read an input file and construct the crystal structure in terms of Cell, SpG and A", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(IO_formats_readn_set_xtal_structure))  ! address of Fortran function to add

    call method_table%add_method("crystal_metrics_set_crystal_cell", &                  ! method name
         "Creates the crystal cell", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_set_crystal_cell))  ! address of Fortran function to add
    
    m = mod_def%init("crysfml_symmetry", "A Python extension for crysFML symmetry", method_table)
  end function init


end module API_init
