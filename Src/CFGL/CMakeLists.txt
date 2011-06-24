# The variable that will store the name of the source file for the built in CONSOLE mode.
SET(crysfgl_LIB_SRCS "")

# Append new sources files for the build of CrysFGL.
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

INCLUDE_DIRECTORIES(${LIB_PATH}/LibC)

# A static version of crysfgl will be built.
ADD_LIBRARY(wcrysfgl STATIC ${crysfgl_LIB_SRCS})
TARGET_LINK_LIBRARIES(wcrysfgl crysfml_common)

# Set the path where the crysfgl library and their associated mod files will be stored.
SET_TARGET_PROPERTIES(wcrysfgl PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH}/LibGL)
SET_TARGET_PROPERTIES(wcrysfgl PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_PATH}/LibGL)

# Some compiler specific flags are defined for the library to build.
# Case of Intel fortran compiler.
IF(${COMPILER_NAME} STREQUAL ifort)

    # The base output path for the libraries.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/Intel)

    # By default, an include directive for Winteracter and wcrysfml are added in case of a build in WINDOW mode.
    SET_SOURCE_FILES_PROPERTIES(${crysfgl_LIB_SRCS} PROPERTIES COMPILE_FLAGS "/I. /I$ENV{WINTER}/lib.if8 /I${LIB_PATH}/LibW")

# Case of gfortran compiler.
ELSEIF(${COMPILER_NAME} STREQUAL gfortran)

    # The base output path for the libraries.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/gfortran)

# Case of g95 compiler.
ELSEIF(${COMPILER_NAME} STREQUAL g95)

    # The base output path for the libraries.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/G95)

# Case of Lahey compiler.
ELSEIF(${COMPILER_NAME} STREQUAL lf95)

    # The base output path for the libraries.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/Lahey)

    # By default, an include directive for Winteracter and wcrysfml are added in case of a build in WINDOW mode.
    SET_SOURCE_FILES_PROPERTIES(${crysfgl_LIB_SRCS} PROPERTIES COMPILE_FLAGS "-mod .;$ENV{WINTER}\lib.l95;${LIB_PATH}/LibW")

    # Excepted the following files for which the compilation flags is slightly different.
    SET_SOURCE_FILES_PROPERTIES(${crysfgl_LIB_SRCS} PROPERTIES COMPILE_FLAGS "-nstchk")

ENDIF()
