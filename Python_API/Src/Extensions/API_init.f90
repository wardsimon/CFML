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
       crystallographic_symmetry_get_latt_trans, &
       crystallographic_symmetry_get_hexa, &
       crystallographic_symmetry_get_numspg, &
       crystallographic_symmetry_get_spg_symb
  
  use API_IO_Formats, only: &
       IO_formats_readn_set_xtal_structure

  use API_Crystal_Metrics, only: &
       crystal_metrics_set_crystal_cell, &
       crystal_metrics_write_crystal_cell, &
       crystal_metrics_get_cell, &
       crystal_metrics_get_ang, &
       crystal_metrics_get_lcell, &
       crystal_metrics_get_lang, &
       crystal_metrics_get_cell_std, &
       crystal_metrics_get_ang_std, &
       crystal_metrics_get_rcell, &
       crystal_metrics_get_rang, &
       crystal_metrics_get_GD, &
       crystal_metrics_get_GR, &
       crystal_metrics_get_Cr_Orth_cel, &
       crystal_metrics_get_Orth_Cr_Cel, &
       crystal_metrics_get_BL_M, &
       crystal_metrics_get_BL_Minv, &
       crystal_metrics_get_cellvol, &
       crystal_metrics_get_rcellvol, &
       crystal_metrics_get_stdvol, &
       crystal_metrics_get_CartType

  use API_Atom_TypeDef, only: &
       atom_typedef_write_atom_list
  
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
    
    call method_table%init(28)
    !--------------------------
    ! Crystallographic Symmetry (6)
    !--------------------------
    call method_table%add_method("crystallographic_symmetry_set_spacegroup", &                  ! method name
         "Creates the space group", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_set_spacegroup))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_write_spacegroup", &                  ! method name
         "Return space group description", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_write_spacegroup))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_latt_trans", &                  ! method name
         "latt_trans getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_latt_trans))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_hexa", &                  ! method name
         "hexa getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_hexa))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_numspg", &                  ! method name
         "numspg getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_numspg))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_spg_symb", &                  ! method name
         "spg_symb getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_spg_symb))  ! address of Fortran function to add
    
    !--------------------------
    ! IO formats (1)
    !--------------------------
    call method_table%add_method("IO_formats_readn_set_xtal_structure", &                  ! method name
         "read an input file and construct the crystal structure in terms of Cell, SpG and A", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(IO_formats_readn_set_xtal_structure))  ! address of Fortran function to add

    !--------------------------
    ! Atom Typedef (1)
    !--------------------------
    call method_table%add_method("atom_typedef_write_atom_list", &                  ! method name
         "Return the atom list description", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_write_atom_list))  ! address of Fortran function to add
    
    !--------------------------
    ! Crystal Metrics (20)
    !--------------------------
    call method_table%add_method("crystal_metrics_set_crystal_cell", &                  ! method name
         "Creates the crystal cell", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_set_crystal_cell))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_write_crystal_cell", &                  ! method name
         "Return the crystal cell description", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_write_crystal_cell))  ! address of Fortran function to add

    call method_table%add_method("crystal_metrics_get_cell", &                  ! method name
         "cell getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_cell))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_ang", &                  ! method name
         "angl getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_ang))  ! address of Fortran function to add

        call method_table%add_method("crystal_metrics_get_lcell", &                  ! method name
         "lcell getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_lcell))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_lang", &                  ! method name
         "lang getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_lang))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_cell_std", &                  ! method name
         "cell_std getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_cell_std))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_ang_std", &                  ! method name
         "ang_std getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_ang_std))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_rcell", &                  ! method name
         "rcell getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_rcell))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_rang", &                  ! method name
         "rang getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_rang))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_GD", &                  ! method name
         "GD getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_GD))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_GR", &                  ! method name
         "GR getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_GR))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_Cr_Orth_Cel", &                  ! method name
         "Cr_Orth_Cel getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_Cr_Orth_Cel))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_Orth_Cr_Cel", &                  ! method name
         "Orth_Cr_Cel getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_Orth_Cr_Cel))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_BL_M", &                  ! method name
         "BL_M getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_BL_M))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_BL_Minv", &                  ! method name
         "BL_Minv getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_BL_Minv))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_cellvol", &                  ! method name
         "cellvol getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_cellvol))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_rcellvol", &                  ! method name
         "rcellvol getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_rcellvol))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_stdvol", &                  ! method name
         "stdvol getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_stdvol))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_get_CartType", &                  ! method name
         "CartType getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_get_CartType))  ! address of Fortran function to add
    
    m = mod_def%init("crysfml_symmetry", "A Python extension for crysFML symmetry", method_table)
  end function init


end module API_init
