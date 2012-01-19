# The minimum version number insuring a correct build.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# The wcrysfml Fortran project.
project(wcrysfml Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CRYSFML_SOURCE}/cmake)

#################################
# Include section
#################################

# Find the winteracter library. If not found, stops here.
find_package(WINTERACTER REQUIRED)

#################################
# Sources section
#################################

# The sources files for crysfml library.
set(SOURCES CFML_IO_Mess.f90 CFML_Optimization_SAn.f90)

# Set the optimization flags.
set_source_files_properties(${SOURCES} PROPERTIES COMPILE_FLAGS ${OPT_FLAGS})

# Those file are generated during the build.
set(CRYSFML_COMMON_OBJECTS_DIR ${PROJECT_BINARY_DIR}/CMakeFiles/crysfml_common.dir)
set(CRYSFML_COMMON_OBJECTS ${CRYSFML_COMMON_SOURCES})
add_suffix(CRYSFML_COMMON_OBJECTS ${CMAKE_Fortran_OUTPUT_EXTENSION})
add_prefix(${CRYSFML_COMMON_OBJECTS_DIR}/ CRYSFML_COMMON_OBJECTS)
set_source_files_properties(${CRYSFML_COMMON_OBJECTS} PROPERTIES GENERATED true)

#################################
# Build section
#################################

# This directory contains the crysfml_common library mod files.
include_directories(${PROJECT_BINARY_DIR}/crysfml_common)

# This directory contains the winteracter library mod files.
include_directories(${WINTERACTER_MOD_DIR})
            
# The wcrysfml library is the GUI version of the library. It is linked using Winteracter.
add_library(wcrysfml STATIC ${SOURCES} ${CRYSFML_COMMON_OBJECTS})

# Add a dependencie to crysfml_common to keep sure that the crysfml_common library will be built first.
add_dependencies(wcrysfml crysfml_common)

# The library is linked to winteracter.
target_link_libraries(wcrysfml ${WINTERACTER_LIBRARY})
        
# Sets the path where to place the mod files for the wcrysfml library.
set_target_properties(wcrysfml PROPERTIES Fortran_MODULE_DIRECTORY ${wcrysfml_BINARY_DIR})

#################################
# Install section
#################################

# The rules for installing the library.
install(TARGETS wcrysfml ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/wcrysfml)

# The rules for installing the mod files. Take care the "/" is on purpose.
install(DIRECTORY ${crysfml_common_BINARY_DIR}/ ${wcrysfml_BINARY_DIR}/ 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/wcrysfml
        FILES_MATCHING
        PATTERN "*.mod"
        PATTERN CMakeFiles EXCLUDE)            
