GET_FILENAME_COMPONENT(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

SET(CMAKE_BUILD_TYPE_INIT Release)
SET(CMAKE_Fortran_FLAGS_INIT "")

IF(${COMPILER_NAME} STREQUAL ifort)

    SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "/debug:full /check /traceback /nologo")
    SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "/O2 /nologo /Qvec-report0")
    
ELSEIF(${COMPILER_NAME} STREQUAL g95)

    SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-O0 -std=f2003 -ftrace=full")
    SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O3 -std=f2003 -funroll-loops -msse2")

ELSEIF(${COMPILER_NAME} STREQUAL gfortran)

    SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-O0 -ftrace=full")
    SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O3 -funroll-loops -msse2")

ELSEIF(${COMPILER_NAME} STREQUAL lf95)

    SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-g -chk")
    SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-o1 -nchk")
        
ELSE()

    MESSAGE(FATAL_ERROR "Unknown compiler")
    
ENDIF()

