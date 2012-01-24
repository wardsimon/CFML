# The minimum version number insuring a correct build.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

set(LIB_NAME crysfml_common)

# The crysfml_common Fortran project.
project(${LIB_NAME} Fortran)

#################################
# Sources section
#################################

# The CFML_GlobalDeps module is OS and compiler dependent.
# WINDOWS
if(WIN32)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
        set(GLOBAL_DEPS CFML_GlobalDeps_Windows_Intel.f90)
    else()
        set(GLOBAL_DEPS CFML_GlobalDeps_Windows.f90)
    endif()
# MacOS
elseif(APPLE)
    set(GLOBAL_DEPS CFML_GlobalDeps_MacOS.f90)
# Unix
elseif(UNIX)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
        set(GLOBAL_DEPS CFML_GlobalDeps_Linux_Intel.f90)
    else()
        set(GLOBAL_DEPS CFML_GlobalDeps_Linux.f90)
    endif()
endif()

# The sources files for crysfml_common library.
set(CRYSFML_COMMON_SOURCES
    ${GLOBAL_DEPS}
    CFML_Atom_Mod.f90
    CFML_Bonds_Table.f90
    CFML_Chem_Scatt.f90
    CFML_Conf_Calc.f90
    CFML_Cryst_Types.f90
    CFML_Diffpatt.f90
    CFML_Export_Vtk.f90
    CFML_FFTs.f90
    CFML_Form_CIF.f90
    CFML_Geom_Calc.f90
    CFML_ILL_Instrm_Data.f90
    CFML_LSQ_TypeDef.f90
    CFML_MagSymm.f90
    CFML_Maps.f90
    CFML_Math_3D.f90
    CFML_Math_Gen.f90
    CFML_Molecules.f90
    CFML_Msfac.f90
    CFML_Optimization.f90
    CFML_Optimization_LSQ.f90
    CFML_Polar.f90
    CFML_Profile_Finger.f90
    CFML_Profile_Functs.f90
    CFML_Profile_TOF.f90
    CFML_Propagk.f90
    CFML_Random.f90
    CFML_Refcodes.f90
    CFML_Reflct_Util.f90
    CFML_Sfac.f90
    CFML_Spher_Harm.f90
    CFML_String_Util.f90
    CFML_SXTAL_Geom.f90
    CFML_Symmetry.f90
    CFML_Sym_Table.f90)

# Set the optimization flags.
set_source_files_properties(${CRYSFML_COMMON_SOURCES} PROPERTIES COMPILE_FLAGS ${OPT_FLAGS})

# Those files need specific optimization flags.
set_source_files_properties(CFML_Bonds_Table.f90 CFML_Chem_Scatt.f90 CFML_Sym_Table.f90 PROPERTIES COMPILE_FLAGS ${OPT_FLAGS1})        

# The crysfml_common library is common to Console and GUI mode.
add_library(${LIB_NAME} STATIC ${CRYSFML_COMMON_SOURCES})

# The directory where the CrySFML modules files will be stored.
set(CRYSFML_COMMON_MODULE_DIRECTORY ${crysfml_common_BINARY_DIR}/crysfml_common_modules)

# Sets the path where to place the mod files for the crysfml_common library.
set_target_properties(${LIB_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CRYSFML_COMMON_MODULE_DIRECTORY})