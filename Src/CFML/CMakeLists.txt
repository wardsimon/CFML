# The minimum version number insuring a correct build.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

option(GUI  "Build crysfml AND wcrysfml libraries" ON)

# The crysfml Fortran project.
project(cfml Fortran)

# Get the "resources" directory name. 
get_filename_component(RESOURCES_DIR ${PROJECT_SOURCE_DIR} PATH)
get_filename_component(RESOURCES_DIR ${RESOURCES_DIR} PATH)
set(RESOURCES_DIR ${RESOURCES_DIR}/resources)

# Append it to the CMake modules path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RESOURCES_DIR})

# Include a module defining global variables.
include(general)

# Include a module defining macros.
include(macros)

if(${GUI})
    # Find the winteracter library. If not found, stops here.
    find_package(WINTERACTER REQUIRED)
endif()

# WINDOWS
if(WIN32)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
        set(COMMON_SOURCES CFML_GlobalDeps_Windows_Intel.f90)
    else()
        set(COMMON_SOURCES CFML_GlobalDeps_Windows.f90)
    endif()
# MacOS
elseif(APPLE)
    set(COMMON_SOURCES CFML_GlobalDeps_MacOS.f90)
# Unix
elseif(UNIX)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
        set(COMMON_SOURCES CFML_GlobalDeps_Linux_Intel.f90)
    else()
        set(COMMON_SOURCES CFML_GlobalDeps_Linux.f90)
    endif()
endif()

# The sources files for crysfml_common library.
set(COMMON_SOURCES
    ${COMMON_SOURCES}
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

# The sources files for crysfml library.
set(SPECIFIC_SOURCES CFML_IO_Mess.f90 CFML_Optimization_SAn.f90)

set(COMMON_OBJS_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/crysfml_common.dir)
set(COMMON_OBJS_FILES ${COMMON_SOURCES})
add_suffix(COMMON_OBJS_FILES ${CMAKE_Fortran_OUTPUT_EXTENSION})
add_prefix(${COMMON_OBJS_DIR}/ COMMON_OBJS_FILES)

# Those file are generated during the build.
set_source_files_properties(${COMMON_OBJS_FILES} PROPERTIES GENERATED true)

# Some compiler specific flags are defined for the library to build.
if(${COMPILER_NAME} STREQUAL gfortran)

    if(WIN32)
        # By default all the source file are compiled with "-std=f2003" flag.
        set_source_files_properties(${COMMON_SOURCES} ${SPECIFIC_SOURCES} 
                                    PROPERTIES COMPILE_FLAGS "-std=f2003")

        # Excepted the following files for which the compilation flags is slightly different for various reasons.
        set_source_files_properties(CFML_Profile_TOF.f90 PROPERTIES COMPILE_FLAGS "")
        set_source_files_properties(CFML_Conf_Calc.f90 PROPERTIES COMPILE_FLAGS "")
        set_source_files_properties(CFML_Optimization_SAn.f90 PROPERTIES COMPILE_FLAGS "")
        set_source_files_properties(CFML_ILL_Instrm_Data.f90 PROPERTIES COMPILE_FLAGS "-std=gnu")        
    endif()

# Case of g95 compiler.
elseif(${COMPILER_NAME} STREQUAL g95)

    if(WIN32)
        # By default all the source file are compiled with "-std=f2003" flag.
        set_source_files_properties(${COMMON_SOURCES} ${SPECIFIC_SOURCES}
                                    PROPERTIES COMPILE_FLAGS "-std=f2003")

        # Excepted the following files for which the compilation flags is slightly different.
        set_source_files_properties(CFML_ILL_Instrm_Data.f90 PROPERTIES COMPILE_FLAGS "")
    endif()
        
endif()

# The crysfml_common library is common to Console and GUI mode.
add_library(crysfml_common STATIC ${COMMON_SOURCES})

# Sets the path where to place the mod files for the crysfml_common library.
set_target_properties(crysfml_common PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/crysfml_common)

# This directory contains the crysfml_common library mod files.
include_directories(${PROJECT_BINARY_DIR}/crysfml_common)

# The crysfml library is the CONSOLE version of the library.
add_library(crysfml STATIC ${SPECIFIC_SOURCES} ${COMMON_OBJS_FILES})

target_link_libraries(crysfml crysfml_common)

# Sets the path where to place the mod files for the crysfml library.
set_target_properties(crysfml PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/crysfml)

get_filename_component(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} PATH)

# The rules to run at install time for the library.
install(TARGETS crysfml ARCHIVE DESTINATION ${INSTALL_PREFIX}/crysfml/LibC)

# The rules to run at install time for the mod files. Take care the "/" is on purpose.
install(DIRECTORY ${PROJECT_BINARY_DIR}/crysfml_common/ ${PROJECT_BINARY_DIR}/crysfml/LibC/
        DESTINATION ${INSTALL_PREFIX}/crysfml/LibC
        FILES_MATCHING
        PATTERN "*.mod"
        PATTERN CMakeFiles EXCLUDE)
        
if(${GUI})

    # This directory contains the winteracter library mod files.
    include_directories(${WINTERACTER_INCLUDE_DIR})
            
    # The wcrysfml library is the GUI version of the library. It is linked using Winteracter.
    add_library(wcrysfml STATIC ${SPECIFIC_SOURCES} ${COMMON_OBJS_FILES})

    # The library is linked to winteracter.
    target_link_libraries(wcrysfml ${WINTERACTER_LIBRARY})
        
    # Sets the path where to place the mod files for the wcrysfml library.
    set_target_properties(wcrysfml PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/wcrysfml)

    # The rules to run at install time for the library.
    install(TARGETS wcrysfml ARCHIVE DESTINATION ${INSTALL_PREFIX}/crysfml/LibW)

    # The rules to run at install time for the mod files. Take care the "/" is on purpose.
    install(DIRECTORY ${PROJECT_BINARY_DIR}/crysfml_common/ ${PROJECT_BINARY_DIR}/crysfml/LibW/ 
            DESTINATION ${INSTALL_PREFIX}/crysfml/LibW 
            FILES_MATCHING
            PATTERN "*.mod"
            PATTERN CMakeFiles EXCLUDE)            
endif()
