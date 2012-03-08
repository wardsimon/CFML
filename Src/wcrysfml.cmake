#################################
# Include section
#################################

# Find the winteracter library. If not found, stops here.
find_package(WINTERACTER REQUIRED)

#################################
# Sources section
#################################

# The sources files for crysfml library.
set(SOURCES CFML_IO_MessWin.f90
            CFML_Optimization_SAn.f90)

# Set the optimization flags.
set_source_files_properties(${SOURCES} PROPERTIES COMPILE_FLAGS ${OPT_FLAGS})

# Those file are generated during the build.
set(CRYSFML_COMMON_OBJECTS_DIR ${PROJECT_BINARY_DIR}/Src/CMakeFiles/crysfml_common.dir)
set(CRYSFML_COMMON_OBJECTS ${CRYSFML_COMMON_SOURCES})
add_suffix(CRYSFML_COMMON_OBJECTS ${CMAKE_Fortran_OUTPUT_EXTENSION})
add_prefix(${CRYSFML_COMMON_OBJECTS_DIR}/ CRYSFML_COMMON_OBJECTS)
set_source_files_properties(${CRYSFML_COMMON_OBJECTS} PROPERTIES GENERATED true)

#################################
# Build section
#################################

set(LIBRARY_NAME wcrysfml)

# This directory contains the crysfml_common library mod files.
include_directories(${CRYSFML_COMMON_MODULE_DIRECTORY})

# This directory contains the winteracter library mod files.
include_directories(${WINTERACTER_MOD_DIR})

# The wcrysfml library is the GUI version of the library.
add_library(${LIBRARY_NAME} STATIC ${SOURCES} ${CRYSFML_COMMON_OBJECTS})

# Add a dependencie to crysfml_common to keep sure that the crysfml_common library will be built first.
add_dependencies(${LIBRARY_NAME} crysfml_common)

# The library is linked to winteracter.
target_link_libraries(${LIBRARY_NAME} ${WINTERACTER_LIBRARY})

# The directory where the wcrysfml specific module files will be stored.
set(WCRYSFML_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/wcrysfml_modules)

# Sets the path where to place the mod files for the wcrysfml library.
set_target_properties(${LIBRARY_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${WCRYSFML_MODULE_DIRECTORY})

#################################
# Install section
#################################

# The rules for installing the library.
install(TARGETS ${LIBRARY_NAME} ARCHIVE DESTINATION ${WCRYSFML_PREFIX})

# The rules for installing the mod files. Take care the "/" is on purpose.
install(DIRECTORY ${CRYSFML_COMMON_MODULE_DIRECTORY}/ ${WCRYSFML_MODULE_DIRECTORY}/
        DESTINATION ${WCRYSFML_PREFIX}
        FILES_MATCHING
        PATTERN "*.mod"
        PATTERN CMakeFiles EXCLUDE)
