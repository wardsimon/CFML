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
       crystallographic_symmetry_del_spacegroup, &
       crystallographic_symmetry_write_spacegroup, &
       crystallographic_symmetry_get_latt_trans, &
       crystallographic_symmetry_get_hexa, &
       crystallographic_symmetry_get_numspg, &
       crystallographic_symmetry_get_spg_symb, &
       crystallographic_symmetry_get_hall, &
       crystallographic_symmetry_get_ghall, &
       crystallographic_symmetry_get_CrystalSys, &
       crystallographic_symmetry_get_Laue, &
       crystallographic_symmetry_get_PG, &
       crystallographic_symmetry_get_Info, &
       crystallographic_symmetry_get_SG_setting, &
       crystallographic_symmetry_get_SPG_lat, &
       crystallographic_symmetry_get_SPG_latsy, &
       crystallographic_symmetry_get_NumLat, &
       crystallographic_symmetry_get_bravais, &
       crystallographic_symmetry_get_centre, &
       crystallographic_symmetry_get_centred, &
       crystallographic_symmetry_get_centre_coord, &
       crystallographic_symmetry_get_numops, &
       crystallographic_symmetry_get_multip, &
       crystallographic_symmetry_get_num_gen, &
       crystallographic_symmetry_get_SymOp, &
       crystallographic_symmetry_get_SymopSymb, &
       crystallographic_symmetry_get_wyckoff, &
       crystallographic_symmetry_get_R_asym_unit, &
       crystallographic_symmetry_get_symmetry_operator_rotation_matrix, &
       crystallographic_symmetry_get_symmetry_operator_trans_matrix, &
       crystallographic_symmetry_get_wyckoff_num_orbit, &
       crystallographic_symmetry_get_wyckoff_orbits, &
       crystallographic_symmetry_get_wyckoff_multp, &
       crystallographic_symmetry_get_wyckoff_site, &
       crystallographic_symmetry_get_wyckoff_norb, &
       crystallographic_symmetry_get_wyckoff_str_orig, &
       crystallographic_symmetry_get_wyckoff_str_orbit
       
  use API_IO_Formats, only: &
       IO_formats_readn_set_xtal_structure

  use API_Structure_Factors, only: &
       structure_factors_structure_factors, &
       structure_factors_write_structure_factors

  use API_Crystal_Metrics, only: &
       crystal_metrics_set_crystal_cell, &
       crystal_metrics_del_crystal_cell, &
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
       atom_typedef_del_atom_list, &
       atom_typedef_get_item, &
       atom_typedef_get_natoms, &
       atom_typedef_atom_from_string, &
       atom_typedef_del_atom, &
       atom_typedef_get_Lab, &
       atom_typedef_get_ChemSymb, &
       atom_typedef_get_SfacSymb, &
       atom_typedef_get_wyck, &
       atom_typedef_get_Active, &
       atom_typedef_get_Z, &
       atom_typedef_get_Mult, &
       atom_typedef_get_X, &
       atom_typedef_get_X_Std, &
       atom_typedef_get_MX, &
       atom_typedef_get_LX, &
       atom_typedef_get_Occ, &
       atom_typedef_get_Occ_Std, &
       atom_typedef_get_MOcc, &
       atom_typedef_get_LOcc, &
       atom_typedef_get_Biso, &
       atom_typedef_get_Biso_std, &
       atom_typedef_get_MBiso, &
       atom_typedef_get_LBiso, &
       atom_typedef_get_Utype, &
       atom_typedef_get_ThType, &
       atom_typedef_get_U, &
       atom_typedef_get_U_std, &
       atom_typedef_get_MU, &
       atom_typedef_get_LU, &
       atom_typedef_get_Ueq, &
       atom_typedef_get_Charge, &
       atom_typedef_get_Moment, &
       atom_typedef_get_Ind, &
       atom_typedef_get_NVar, &
       atom_typedef_get_VarF, &
       atom_typedef_get_MVarF, &
       atom_typedef_get_LVarF, &
       atom_typedef_get_AtmInfo, &
       atom_typedef_get_m_xyz, &
       atom_typedef_get_sm_xyz, &
       atom_typedef_get_Mm_xyz, &
       atom_typedef_get_Lm_xyz, &
       atom_typedef_set_item, &
       atom_typedef_write_atom_list

  use API_Reflections_Utilities, only: &
       reflections_utilities_hkl_uni_reflist

  use API_Diffraction_Patterns, only: &
       diffraction_patterns_compute_powder_pattern, &
       diffraction_patterns_get_x, &
       diffraction_patterns_get_y
  
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
    
    call method_table%init(110)
    !--------------------------
    ! Diffraction Patterns (3)
    !--------------------------
    call method_table%add_method("diffraction_patterns_compute_powder_pattern", &                  ! method name
         "compute the powder diffraction pattern from some experimental conditions and a set of reflections", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(diffraction_patterns_compute_powder_pattern))  ! address of Fortran function to add

    call method_table%add_method("diffraction_patterns_get_x", &                  ! method name
         "Get x array", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(diffraction_patterns_get_x))  ! address of Fortran function to add

    call method_table%add_method("diffraction_patterns_get_y", &                  ! method name
         "Get y array", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(diffraction_patterns_get_y))  ! address of Fortran function to add

    !--------------------------
    ! Crystallographic Symmetry (37)
    !--------------------------
    call method_table%add_method("crystallographic_symmetry_set_spacegroup", &                  ! method name
         "Creates the space group", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_set_spacegroup))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_del_spacegroup", &                  ! method name
         "Delete the space group", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_del_spacegroup))  ! address of Fortran function to add

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

    call method_table%add_method("crystallographic_symmetry_get_hall", &                  ! method name
         "hall getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_hall))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_ghall", &                  ! method name
         "ghall getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_ghall))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_CrystalSys", &                  ! method name
         "CrystalSys getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_CrystalSys))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_Laue", &                  ! method name
         "Laue getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_Laue))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_PG", &                  ! method name
         "PG getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_PG))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_Info", &                  ! method name
         "Info getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_Info))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_SG_setting", &                  ! method name
         "SG_setting getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_SG_setting))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_SPG_lat", &                  ! method name
         "SPG_lat getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_SPG_lat))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_SPG_latsy", &                  ! method name
         "SPG_latsy getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_SPG_latsy))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_NumLat", &                  ! method name
         "NumLat getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_NumLat))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_bravais", &                  ! method name
         "bravais getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_bravais))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_centre", &                  ! method name
         "centre getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_centre))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_centred", &                  ! method name
         "centred getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_centred))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_centre_coord", &                  ! method name
         "centre_coord getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_centre_coord))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_numops", &                  ! method name
         "numops getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_numops))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_multip", &                  ! method name
         "multip getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_multip))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_num_gen", &                  ! method name
         "num_gen getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_num_gen))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_SymOp", &                  ! method name
         "SymOp getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_SymOp))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_SymopSymb", &                  ! method name
         "SymopSymb getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_SymopSymb))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_wyckoff", &                  ! method name
         "wyckoff getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_wyckoff))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_R_asym_unit", &                  ! method name
         "R_asym_unit getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_R_asym_unit))  ! address of Fortran function to add
    
    call method_table%add_method("crystallographic_symmetry_get_symmetry_operator_rotation_matrix", &                  ! method name
         "Symmetry operator rotation matrix getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_symmetry_operator_rotation_matrix))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_symmetry_operator_trans_matrix", &                  ! method name
         "Symmetry operator translation matrix getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_symmetry_operator_trans_matrix))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_wyckoff_num_orbit", &                  ! method name
         "Wyckoff Type num orbit getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_wyckoff_num_orbit))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_wyckoff_orbits", &                  ! method name
         "Wyckoff Type orbit getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_wyckoff_orbits))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_wyckoff_multp", &                  ! method name
         "Wyckoff Orbit multiplicity getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_wyckoff_multp))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_wyckoff_site", &                  ! method name
         "Wyckoff Orbit site getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_wyckoff_site))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_wyckoff_norb", &                  ! method name
         "Wyckoff Orbit norb getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_wyckoff_norb))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_wyckoff_str_orig", &                  ! method name
         "Wyckoff Orbit str_orig getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_wyckoff_str_orig))  ! address of Fortran function to add

    call method_table%add_method("crystallographic_symmetry_get_wyckoff_str_orbit", &                  ! method name
         "Wyckoff Orbit str_orbit getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystallographic_symmetry_get_wyckoff_str_orbit))  ! address of Fortran function to add


    !--------------------------
    ! IO formats (1)
    !--------------------------
    call method_table%add_method("IO_formats_readn_set_xtal_structure", &                  ! method name
         "read an input file and construct the crystal structure in terms of Cell, SpG and A", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(IO_formats_readn_set_xtal_structure))  ! address of Fortran function to add

    !--------------------------
    ! Atom Typedef (44)
    !--------------------------
    call method_table%add_method("atom_typedef_del_atom_list", &                  ! method name
         "Delete an atom list", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_del_atom_list))  ! address of Fortran function to add

    call method_table%add_method("atom_typedef_set_item", &                  ! method name
         "Set an item in the atom list", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_set_item))  ! address of Fortran function to add

    call method_table%add_method("atom_typedef_get_item", &                  ! method name
         "Get a specific atom", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_item))  ! address of Fortran function to add

    call method_table%add_method("atom_typedef_get_natoms", &                  ! method name
         "natoms getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_natoms))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_write_atom_list", &                  ! method name
         "Return the atom list description", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_write_atom_list))  ! address of Fortran function to add

    call method_table%add_method("atom_typedef_atom_from_string", &                  ! method name
         "Create an atom from a string", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_atom_from_string))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_del_atom", &                  ! method name
         "Delete an atom", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_del_atom))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Lab", &                  ! method name
         "Lab getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Lab))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_ChemSymb", &                  ! method name
         "ChemSymb getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_ChemSymb))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_SfacSymb", &                  ! method name
         "SfacSymb getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_SfacSymb))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_wyck", &                  ! method name
         "wyck getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_wyck))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Active", &                  ! method name
         "Active getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Active))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Z", &                  ! method name
         "Z getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Z))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Mult", &                  ! method name
         "Mult getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Mult))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_X", &                  ! method name
         "X getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_X))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_X_Std", &                  ! method name
         "X_Std getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_X_Std))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_MX", &                  ! method name
         "MX getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_MX))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_LX", &                  ! method name
         "LX getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_LX))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Occ", &                  ! method name
         "Occ getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Occ))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Occ_Std", &                  ! method name
         "Occ_Std getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Occ_Std))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_MOcc", &                  ! method name
         "MOcc getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_MOcc))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_LOcc", &                  ! method name
         "LOcc getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_LOcc))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Biso", &                  ! method name
         "Biso getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Biso))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Biso_std", &                  ! method name
         "Biso_std getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Biso_std))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_MBiso", &                  ! method name
         "MBiso getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_MBiso))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_LBiso", &                  ! method name
         "LBiso getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_LBiso))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Utype", &                  ! method name
         "Utype getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Utype))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_ThType", &                  ! method name
         "ThType getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_ThType))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_U", &                  ! method name
         "U getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_U))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_U_std", &                  ! method name
         "U_std getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_U_std))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_MU", &                  ! method name
         "MU getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_MU))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_LU", &                  ! method name
         "LU getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_LU))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Ueq", &                  ! method name
         "Ueq getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Ueq))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Charge", &                  ! method name
         "Charge getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Charge))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Moment", &                  ! method name
         "Moment getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Moment))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Ind", &                  ! method name
         "Ind getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Ind))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_NVar", &                  ! method name
         "NVar getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_NVar))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_VarF", &                  ! method name
         "VarF getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_VarF))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_MVarF", &                  ! method name
         "MVarF getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_MVarF))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_LVarF", &                  ! method name
         "LVarF getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_LVarF))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_AtmInfo", &                  ! method name
         "AtmInfo getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_AtmInfo))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_m_xyz", &                  ! method name
         "m_xyz getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_m_xyz))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_sm_xyz", &                  ! method name
         "sm_xyz getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_sm_xyz))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Mm_xyz", &                  ! method name
         "Mm_xyz getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Mm_xyz))  ! address of Fortran function to add
    
    call method_table%add_method("atom_typedef_get_Lm_xyz", &                  ! method name
         "Lm_xyz getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(atom_typedef_get_Lm_xyz))  ! address of Fortran function to add
    !--------------------------
    ! Reflection Utilities (1)
    !--------------------------
    call method_table%add_method("reflections_utilities_hkl_uni_reflist", &                  ! method name
         "Return the list of reflections", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(reflections_utilities_hkl_uni_reflist))  ! address of Fortran function to add

    !--------------------------
    ! Structure Factors (2)
    !--------------------------
    call method_table%add_method("structure_factors_structure_factors", &                  ! method name
         "Computes structure factors", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(structure_factors_structure_factors))  ! address of Fortran function to add

    call method_table%add_method("structure_factors_write_structure_factors", &                  ! method name
         "Print structure factors", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(structure_factors_write_structure_factors))  ! address of Fortran function to add
    !--------------------------
    ! Crystal Metrics (21)
    !--------------------------
    call method_table%add_method("crystal_metrics_set_crystal_cell", &                  ! method name
         "Creates the crystal cell", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_set_crystal_cell))  ! address of Fortran function to add
    
    call method_table%add_method("crystal_metrics_del_crystal_cell", &                  ! method name
         "Creates the crystal cell", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(crystal_metrics_del_crystal_cell))  ! address of Fortran function to add

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
