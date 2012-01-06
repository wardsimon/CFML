# The minimum version number insuring a correct build.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# The crysfgl Fortran project.
project(cfgl Fortran)

# The name of the library.
set(LIBRARY_NAME wcrysfgl)

# Get the "resources" directory name. 
get_filename_component(RESOURCES_DIR ${PROJECT_SOURCE_DIR} PATH)
get_filename_component(RESOURCES_DIR ${RESOURCES_DIR} PATH)
set(RESOURCES_DIR ${RESOURCES_DIR}/resources)

# Append it to the CMake modules path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RESOURCES_DIR})

# Include a module defining global variables.
include(general)

# Append it to the CMake modules path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RESOURCES_DIR})

# Find the winteracter library. If not found, stops here.
find_package(WINTERACTER REQUIRED)

# Find the wcrysfml library. If not found, stops here.
find_package(WCRYSFML REQUIRED)

# The source files.
set(SOURCES
    WCFGL_Atomic_Table.f90
    WCFGL_Atom_Tree.f90
    WCFGL_Bond_Tree.f90
    WCFGL_Chull3D.f90
    WCFGL_Constant.f90
    WCFGL_Display.f90
    WCFGL_Geometry.f90
    WCFGL_Glatom.f90
    WCFGL_Glbond.f90
    WCFGL_Glpoly.f90
    WCFGL_IO.f90
    WCFGL_Matom_Tree.f90
    WCFGL_Metrix.f90
    WCFGL_Objects_Definition.f90
    WCFGL_Poly_Tree.f90
    WCFGL_Quaternion.f90
    WCFGL_Trackball.f90
)

# The path for winteracter library mod files.
include_directories(${WINTERACTER_INCLUDE_DIR})

# The path for wcrysfml library mod files.
include_directories(${WCRYSFML_INCLUDE_DIR})

# Add the crysfgl library to the project.
add_library(${LIBRARY_NAME} STATIC ${SOURCES})

# The library is linked to winteracter and wcrysfml.
target_link_libraries(${LIBRARY_NAME} ${WINTERACTER_LIBRARY} ${WCRYSFML_LIBRARY})

get_filename_component(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} PATH)

# The rules to run at install time for the library.
install(TARGETS ${LIBRARY_NAME} ARCHIVE DESTINATION ${INSTALL_PREFIX}/crysfml/LibGL)

# The rules to run at install time for the mod files. Take care the "/" is on purpose.
install(DIRECTORY ${PROJECT_BINARY_DIR}/ 
        DESTINATION ${INSTALL_PREFIX}/crysfml/LibGL 
        FILES_MATCHING
        PATTERN "*.mod"
        PATTERN CMakeFiles EXCLUDE)