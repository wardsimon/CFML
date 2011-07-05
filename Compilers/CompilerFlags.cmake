#####################################################
# Infos
#####################################################

# CMake file for crysfml project
# Written by Eric Pellegrini
# Institut Laue Langevin
# Grenoble, France

#####################################################



# Extract the compiler name (without extension) from its path and puts its in COMPILER_NAME variable.
GET_FILENAME_COMPONENT(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

# By default the build is in release mode (no debug).
SET(CMAKE_BUILD_TYPE_INIT Release)

# Initialize the compiler flag to a null string.
SET(CMAKE_Fortran_FLAGS_INIT "")

# Default flags for ifort compiler.
IF(${COMPILER_NAME} STREQUAL ifort)

    # The directory that will store the libraries built with ifort.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/Intel)

    IF(WIN32)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-debug:full /check /traceback /nologo")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "/O2 /nologo /Qvec-report0")
    ELSEIF(APPLE)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-c -g -warn")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-c -O -warn -vec-report0")
    ELSEIF(UNIX)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-c -g -warn")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-c -O -warn -vec-report0")
    ENDIF()
    
# Default flags for g95 compiler.
ELSEIF(${COMPILER_NAME} STREQUAL g95)

    # The directory that will store the libraries built with g95.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/G95)

    IF(WIN32)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT  "-O0 -ftrace=full")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O3 -funroll-loops -msse2")
    ELSEIF(APPLE)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT  "-c -g -Wall")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-c -O")
    ELSEIF(UNIX)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT  "-c -g -Wall")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-c -O")
    ENDIF()
    
# Default flags for gfortran compiler.
ELSEIF(${COMPILER_NAME} STREQUAL gfortran)

    # The directory that will store the libraries built with gfortran.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/gfortran)

    IF(WIN32)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-O0 -ftrace=full")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O3 -funroll-loops -msse2")
    ELSEIF(APPLE)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-c -g -Wall -m32")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-c -O -m32")
    ELSEIF(UNIX)
        SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-c -g -Wall")
        SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-c -O")    
    ENDIF()

# Default flags for lf95 compiler.
ELSEIF(${COMPILER_NAME} STREQUAL lf95)
    
    # The directory that will store the libraries built with lf95.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/Lahey)

    SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-g -chk")
    SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-o1 -nchk")
        
# Case of an unknown compiler. Throw an error.    
ELSE()

    MESSAGE(FATAL_ERROR "Unknown compiler")
    
ENDIF()

