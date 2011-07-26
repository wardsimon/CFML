# CMake file for crysfgl project
# Written by Eric Pellegrini
# Institut Laue Langevin
# Grenoble, France

# The source files.
set(crysfgl_SRCS
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
include_directories(${WINTER_MOD_DIR})

# The path for crysfml_common library mod files.
include_directories(${PROJECT_BINARY_DIR}/Lib/Common)

# The path for wcrysfml library mod files.
include_directories(${PROJECT_BINARY_DIR}/Lib/LibW)

# The path for winteracter library.
link_directories(${WINTER_LIB_DIR})

# The path for crysfml_common library.
link_directories(${PROJECT_BINARY_DIR}/Lib/Common)

# The path for wcrysfml library mod files.
link_directories(${PROJECT_BINARY_DIR}/Lib/LibW)

# Add the crysfgl library to the project.
add_library(wcrysfgl STATIC ${crysfgl_SRCS})

# The crysfml_common will be linked to crysfgl library.
target_link_libraries(wcrysfgl crysfml_common)

target_link_libraries(wcrysfgl winter)

# Sets the path where to place the crysfgl library.
set_target_properties(wcrysfgl PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Lib/LibGL)

# Sets the path where to place the mod files for the crysfgl library.
set_target_properties(wcrysfgl PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Lib/LibGL)

add_custom_command(TARGET wcrysfgl
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${LIB_PATH}/LibGL
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/Lib/LibGL ${CRYSFGL_LIB_DIR})
