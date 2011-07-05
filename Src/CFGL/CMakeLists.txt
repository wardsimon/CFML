#####################################################
# Infos
#####################################################

# CMake file for crysfml project
# Written by Eric Pellegrini
# Institut Laue Langevin
# Grenoble, France

#####################################################



#####################################################
# Define the source files for crysfml library       #
#####################################################

# This will store the source files for the crysfgl library.
SET(crysfgl_LIB_SRCS "")

# The source files of the crysfgl library.
SET(crysfgl_LIB_SRCS
    ${crysfgl_LIB_SRCS}
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

#####################################################



#####################################################
# Define targets properties                         #
#####################################################

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/Lib/LibC ${CMAKE_BINARY_DIR}/Lib/LibW)

# Add the crysfgl library to the project.
ADD_LIBRARY(wcrysfgl STATIC ${crysfgl_LIB_SRCS})

# Sets the path where to place the crysfgl library.
SET_TARGET_PROPERTIES(wcrysfgl PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Lib/LibGL)

# Sets the path where to place the mod files for the crysfgl library.
SET_TARGET_PROPERTIES(wcrysfgl PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/Lib/LibGL)

IF(WIN32)
    SET_TARGET_PROPERTIES(wcrysfgl PROPERTIES COMPILE_FLAGS "/I$ENV{WINTER}/lib.if8")
ELSEIF(APPLE)
    SET_TARGET_PROPERTIES(wcrysfgl PROPERTIES COMPILE_FLAGS "-I$ENV{WINTER}/lib.ifi")        
ELSEIF(UNIX)
    SET_TARGET_PROPERTIES(wcrysfgl PROPERTIES COMPILE_FLAGS "-I$ENV{WINTER}/lib.if8")    
ENDIF()

# The crysfml_common will be linked to crysfgl library.
TARGET_LINK_LIBRARIES(wcrysfgl crysfml_common)

ADD_CUSTOM_COMMAND(
    TARGET wcrysfgl
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${LIB_PATH}/LibGL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/Lib/LibGL ${LIB_PATH}/LibGL)


# Some compiler specific flags are defined for the library to build.
# Case of Intel fortran compiler.
#IF(${COMPILER_NAME} STREQUAL lf95)
    
#    # By default, an include directive for Winteracter and wcrysfml are added in case of a build in GUI mode.
#    SET_SOURCE_FILES_PROPERTIES(${crysfgl_LIB_SRCS} PROPERTIES COMPILE_FLAGS "-mod .;$ENV{WINTER}\lib.l95;${LIB_PATH}/LibW")

#    # Excepted the following files for which the compilation flags is slightly different.
#    SET_SOURCE_FILES_PROPERTIES(${crysfgl_LIB_SRCS} PROPERTIES COMPILE_FLAGS "-nstchk")

#ENDIF()
