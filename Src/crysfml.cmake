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
set(CRYSFML_COMMON_OBJECTS_DIR ${crysfml_common_BINARY_DIR}/CMakeFiles/crysfml_common.dir)
set(CRYSFML_COMMON_OBJECTS ${CRYSFML_COMMON_SOURCES})
add_suffix(CRYSFML_COMMON_OBJECTS ${CMAKE_Fortran_OUTPUT_EXTENSION})
add_prefix(${CRYSFML_COMMON_OBJECTS_DIR}/ CRYSFML_COMMON_OBJECTS)
set_source_files_properties(${CRYSFML_COMMON_OBJECTS} PROPERTIES GENERATED true)

#################################
# Build section
#################################

# This directory contains the crysfml_common library mod files.
include_directories(${CRYSFML_COMMON_MODULE_DIRECTORY})

# The crysfml library is the CONSOLE version of the library.
add_library(crysfml STATIC ${SOURCES} ${CRYSFML_COMMON_OBJECTS})

# Add a dependencie to crysfml_common to keep sure that the crysfml_common library will be built first.
add_dependencies(crysfml crysfml_common)

# The directory where the crysfml specific module files will be stored.
set(CRYSFML_MODULE_DIRECTORY ${crysfml_BINARY_DIR}/crysfml_modules)

# Sets the path where to place the mod files for the crysfml library.
set_target_properties(crysfml PROPERTIES Fortran_MODULE_DIRECTORY ${CRYSFML_MODULE_DIRECTORY})

#################################
# Install section
#################################

# The rules for installing the library.
install(TARGETS crysfml ARCHIVE DESTINATION ${CRYSFML_PREFIX})

# The rules for installing the mod files. Take care the "/" is on purpose.
install(DIRECTORY ${CRYSFML_COMMON_MODULE_DIRECTORY}/ ${CRYSFML_MODULE_DIRECTORY}/
        DESTINATION ${CRYSFML_PREFIX}
        FILES_MATCHING
        PATTERN "*.mod"
        PATTERN CMakeFiles EXCLUDE)
