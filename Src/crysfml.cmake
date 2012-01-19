# The minimum version number insuring a correct build.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# The crysfml Fortran project.
project(crysfml Fortran)

#################################
# Sources section
#################################

# The sources files for crysfml library.
set(SOURCES CFML_IO_Mess.f90
            CFML_Optimization_SAn.f90)

# Set the optimization flags.
set_source_files_properties(${SOURCES} PROPERTIES COMPILE_FLAGS ${OPT_FLAGS})

# Those file are generated during the build.
set(CRYSFML_COMMON_OBJECTS_DIR ${PROJECT_BINARY_DIR}/CMakeFiles/crysfml_common.dir)
set(CRYSFML_COMMON_OBJECTS ${CRYSFML_COMMON_SOURCES})
add_suffix(CRYSFML_COMMON_OBJECTS ${CMAKE_Fortran_OUTPUT_EXTENSION})
add_prefix(${CRYSFML_COMMON_OBJECTS_DIR}/ CRYSFML_COMMON_OBJECTS)

#################################
# Build section
#################################

# This directory contains the crysfml_common library mod files.
include_directories(${PROJECT_BINARY_DIR}/crysfml_common)

# The crysfml library is the CONSOLE version of the library.
add_library(crysfml STATIC ${SOURCES} ${CRYSFML_COMMON_OBJECTS})

# Add a dependencie to crysfml_common to keep sure that the crysfml_common library will be built first.
add_dependencies(crysfml crysfml_common)

# Sets the path where to place the mod files for the crysfml library.
set_target_properties(crysfml PROPERTIES Fortran_MODULE_DIRECTORY ${crysfml_BINARY_DIR})

#################################
# Install section
#################################

# The rules for installing the library.
install(TARGETS crysfml ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/crysfml)

# The rules for installing the mod files. Take care the "/" is on purpose.
install(DIRECTORY ${crysfml_common_BINARY_DIR}/ ${crysfml_BINARY_DIR}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/crysfml
        FILES_MATCHING
        PATTERN "*.mod"
        PATTERN CMakeFiles EXCLUDE)