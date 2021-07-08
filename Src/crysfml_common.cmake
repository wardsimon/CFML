#################################
# Sources section
#################################

# The CFML_GlobalDeps module is OS and compiler dependent.
# WINDOWS
if(WIN32)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
       if(ARCH32)
           set(GLOBAL_DEPS CFML_GlobalDeps_Windows_Intel.f90)
       else()
           set(GLOBAL_DEPS CFML_GlobalDeps_Windows_Intel64.f90)
       endif()
       set(STRING_UTILS CFML_String_Utilities.f90)
    else()
        set(GLOBAL_DEPS CFML_GlobalDeps_Windows.f90)
        set(STRING_UTILS CFML_String_Utilities_gf.f90)
    endif()
# MacOS
elseif(APPLE)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
        set(GLOBAL_DEPS CFML_GlobalDeps_MacOS_Intel.f90)
        set(STRING_UTILS CFML_String_Utilities.f90)
    else()
        set(GLOBAL_DEPS CFML_GlobalDeps_MacOS.f90)
        set(STRING_UTILS CFML_String_Utilities_gf.f90)
    endif()
# Unix
elseif(UNIX)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
        set(GLOBAL_DEPS CFML_GlobalDeps_Linux_Intel.f90)
        set(STRING_UTILS CFML_String_Utilities.f90)
    else()
        set(GLOBAL_DEPS CFML_GlobalDeps_Linux.f90)
        set(STRING_UTILS CFML_String_Utilities_gf.f90)
    endif()
endif()

# The sources files for crysfml_common library.
set(CRYSFML_COMMON_SOURCES
    ${GLOBAL_DEPS}
    CFML_EisPack.f90
    CFML_Atom_TypeDef.f90
    CFML_Bond_Tables.f90
    CFML_BVSpar.f90
    CFML_Scattering_Chemical_Tables.f90
    CFML_BVS_Energy_Calc.f90
    CFML_Crystal_Metrics.f90
    CFML_Diffraction_Patterns.f90
    CFML_Export_Vtk.f90
    CFML_Extinction_Correction.f90
    CFML_FFT.f90
    CFML_IO_Formats.f90
    CFML_Geometry_Calc.f90
    CFML_Geometry_SXTAL.f90
    CFML_ILL_Instrm_Data.f90
    CFML_LSQ_TypeDef.f90
    CFML_Magnetic_Groups.f90
    CFML_Magnetic_Symmetry.f90
    CFML_Maps_Calculations.f90
    CFML_Math_3D.f90
    CFML_Math_General.f90
    CFML_Molecular_Crystals.f90
    CFML_Magnetic_Structure_Factors.f90
    CFML_Optimization_General.f90
    CFML_Optimization_LSQ.f90
    CFML_Percolation.f90
    CFML_Polarimetry.f90
    CFML_PowderProfiles_Finger.f90
    CFML_PowderProfiles_CW.f90
    CFML_PowderProfiles_TOF.f90
    CFML_Propagation_Vectors.f90
    CFML_Random_Generators.f90
    CFML_Keywords_Code_Parser.f90
    CFML_Reflections_Utilities.f90
    CFML_Structure_Factors.f90
    CFML_Spherical_Harmonics.f90
    ${STRING_UTILS}
    CFML_EoS.f90
    CFML_Spherical_Harmonics.f90
    CFML_Rational_Arithmetic.f90
    CFML_Crystallographic_Symmetry.f90
    CFML_Symmetry_Tables.f90)

# Set the optimization flags.
set_source_files_properties(${CRYSFML_COMMON_SOURCES} PROPERTIES COMPILE_FLAGS ${OPT_FLAGS})

# Those files need specific optimization flags.
set_source_files_properties(CFML_BVSpar.f90 CFML_Bond_Table.f90 CFML_Scattering_Chemical_Tables.f90 CFML_Symmetry_Tables.f90 PROPERTIES COMPILE_FLAGS ${OPT_FLAGS0})
set_source_files_properties(CFML_Profile_TOF.f90 PROPERTIES COMPILE_FLAGS ${OPT_FLAGS1})

#################################
# Build section
#################################

set(LIBRARY_NAME crysfml_common)

# The crysfml_common library is common to Console and GUI mode.
add_library(${LIBRARY_NAME} STATIC ${CRYSFML_COMMON_SOURCES})

# The directory where the CrySFML modules files will be stored.
set(CRYSFML_COMMON_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/crysfml_common_modules CACHE INTERNAL "")

# Sets the path where to place the mod files for the crysfml_common library.
set_target_properties(${LIBRARY_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CRYSFML_COMMON_MODULE_DIRECTORY})
