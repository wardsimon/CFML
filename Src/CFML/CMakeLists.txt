MACRO(ADD_PREFIX prefix rootlist)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${prefix}${root})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_PREFIX)

MACRO(ADD_SUFFIX rootlist suffix)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${root}${suffix})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_SUFFIX)

# This will store the source files common to a build in CONSOLE and WINDOW mode.
SET(crysfml_COMMON_LIB_SRCS "")

#####################################################
# Case of the OS and compiler dependant source files.
#####################################################
# Windows
IF(WIN32)

    # Intel Fortran compiler
    IF(${COMPILER_NAME} STREQUAL ifort)

        SET(crysfml_COMMON_LIB_SRCS ${crysfml_COMMON_LIB_SRCS} CFML_GlobalDeps_Windows_intel.f90)

    ELSE()

        SET(crysfml_COMMON_LIB_SRCS ${crysfml_COMMON_LIB_SRCS} CFML_GlobalDeps_Windows.f90)

        # Lahey Fortran compiler.
        IF(${COMPILER_NAME} STREQUAL lf95)

            SET(crysfml_COMMON_LIB_SRCS ${crysfml_COMMON_LIB_SRCS} f2kcli.f90)

        ENDIF()

    ENDIF()

# MacOS
ELSEIF(APPLE)

    SET(crysfml_COMMON_LIB_SRCS ${crysfml_COMMON_LIB_SRCS} CFML_GlobalDeps_MacOS.f90)

# Unix
ELSEIF(UNIX)

    # Intel Fortran compiler
    IF(${COMPILER_NAME} STREQUAL ifort)

        SET(crysfml_COMMON_LIB_SRCS ${crysfml_COMMON_LIB_SRCS} CFML_GlobalDeps_Linux_Intel.f90)

    ELSE()

        SET(crysfml_COMMON_LIB_SRCS ${crysfml_COMMON_LIB_SRCS} CFML_GlobalDeps_Linux.f90)

    ENDIF()

ENDIF()

####################################################################
# The other common files between a build in CONSOLE and WINDOW mode.
####################################################################
SET(crysfml_COMMON_LIB_SRCS
    ${crysfml_COMMON_LIB_SRCS}
    CFML_Atom_Mod.f90
    CFML_Bonds_Table.f90
    CFML_Chem_Scatt.f90
    CFML_Conf_Calc.f90
    CFML_Cryst_Types.f90
    CFML_Diffpatt.f90
    CFML_FFTs.f90
    CFML_Form_CIF.f90
    CFML_Geom_Calc.f90
    CFML_ILL_Instrm_Data.f90
    CFML_LSQ_TypeDef.f90
    CFML_MagSymm.f90
    CFML_Maps.f90
    CFML_Math_3D.f90
    CFML_Math_Gen.f90
    CFML_Molecules.f90
    CFML_Msfac.f90
    CFML_Optimization.f90
    CFML_Optimization_LSQ.f90
    CFML_Polar.f90
    CFML_Profile_Finger.f90
    CFML_Profile_Functs.f90
    CFML_Profile_TOF.f90
    CFML_Propagk.f90
    CFML_Random.f90
    CFML_Refcodes.f90
    CFML_Reflct_Util.f90
    CFML_Sfac.f90
    CFML_Spher_Harm.f90
    CFML_String_Util.f90
    CFML_SXTAL_Geom.f90
    CFML_Symmetry.f90
    CFML_Sym_Table.f90
)

#######################################################################################
# Case of the extra files that have to be built differently in CONSOLE and WINDOW mode.
#######################################################################################

# This will store the extra files for a build in CONSOLE mode.
SET(crysfml_CONSOLE_LIB_SRCS
    CFML_IO_Mess.f90
    CFML_Optimization_SAn.f90
)

# This will store the extra files for a build in WINDOW mode.
SET(crysfml_WINDOW_LIB_SRCS
    CFML_Optimization_SAn.f90
)
# Case of a built in WINDOW mode: some other specific source files are added.
IF(${WINDOW})

    # Case of The Lahey compiler with a RWIN build (-DRWIN=True).
    IF(${COMPILER_NAME} STREQUAL lf95 AND ${RWIN})

        SET(crysfml_WINDOW_LIB_SRCS ${crysfml_WINDOW_LIB_SRCS} CFML_IO_MessRW.f90)

    ELSE()

        SET(crysfml_WINDOW_LIB_SRCS ${crysfml_WINDOW_LIB_SRCS} CFML_IO_MessWin.f90)

    ENDIF()

    # A Winteracter-supported version of the library will be built.
    # ADD_LIBRARY(wcrysfml STATIC ${wcrysfml_LIB_SRCS})

ENDIF()

# A CONSOLE version of crysfml library will be built.
ADD_LIBRARY(crysfml_common STATIC ${crysfml_COMMON_LIB_SRCS})
SET_TARGET_PROPERTIES(crysfml_common PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH}/Common)
SET_TARGET_PROPERTIES(crysfml_common PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_PATH}/LibC)

SET(COMMON_OBJS_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/crysfml_common.dir)
SET(COMMON_OBJS_FILES ${crysfml_COMMON_LIB_SRCS})
ADD_SUFFIX(COMMON_OBJS_FILES ${CMAKE_Fortran_OUTPUT_EXTENSION})
ADD_PREFIX(${COMMON_OBJS_DIR}/ COMMON_OBJS_FILES)

SET_SOURCE_FILES_PROPERTIES(
  ${COMMON_OBJS_FILES}
  PROPERTIES GENERATED true
)

ADD_LIBRARY(crysfml STATIC ${crysfml_CONSOLE_LIB_SRCS} ${COMMON_OBJS_FILES})
TARGET_LINK_LIBRARIES(crysfml crysfml_common)
SET_TARGET_PROPERTIES(crysfml PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH}/LibC)
SET_TARGET_PROPERTIES(crysfml PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_PATH}/LibC)

# Set the path where the crysfml and wcrysfml libraries and their associated mod files will be stored.
IF(${WINDOW})

    INCLUDE_DIRECTORIES(${LIB_PATH}/LibC)
    
    ADD_LIBRARY(wcrysfml STATIC ${crysfml_WINDOW_LIB_SRCS} ${COMMON_OBJS_FILES})
    TARGET_LINK_LIBRARIES(wcrysfml crysfml_common)

    SET_TARGET_PROPERTIES(wcrysfml PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH}/LibW)
    SET_TARGET_PROPERTIES(wcrysfml PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_PATH}/LibW)

ENDIF()

# Some compiler specific flags are defined for the library to build.
# Case of Intel fortran compiler.
IF(${COMPILER_NAME} STREQUAL ifort)

    # The base output path for the libraries.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/Intel)

    # An include directive for Winteracter is added in case of a build in WINDOW mode.
    IF(${WINDOW})

        SET_TARGET_PROPERTIES(wcrysfml PROPERTIES COMPILE_FLAGS "/I$ENV{WINTER}/lib.if8")

    ENDIF()

# Case of gfortran compiler.
ELSEIF(${COMPILER_NAME} STREQUAL gfortran)

    # The base output path for the libraries.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/gfortran)

    # By default all the source file are compiled with "-std=f2003" flag.
    SET_SOURCE_FILES_PROPERTIES(${crysfml_LIB_SRCS} PROPERTIES COMPILE_FLAGS "-std=f2003")

    # Excepted the following files for which the compilation flags is slightly different.
    SET_SOURCE_FILES_PROPERTIES(CFML_Profile_TOF.f90 PROPERTIES COMPILE_FLAGS "")
    SET_SOURCE_FILES_PROPERTIES(CFML_Conf_Calc.f90 PROPERTIES COMPILE_FLAGS "")
    SET_SOURCE_FILES_PROPERTIES(CFML_Optimization_SAn.f90 PROPERTIES COMPILE_FLAGS "")
    SET_SOURCE_FILES_PROPERTIES(CFML_ILL_Instrm_Data.f90 PROPERTIES COMPILE_FLAGS "-std=gnu")

# Case of g95 compiler.
ELSEIF(${COMPILER_NAME} STREQUAL g95)

    # The base output path for the libraries.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/G95)

    # By default all the source file are compiled with "-std=f2003" flag.
    SET_SOURCE_FILES_PROPERTIES(${crysfml_LIB_SRCS} PROPERTIES COMPILE_FLAGS "-std=f2003")

    # Excepted the following files for which the compilation flags is slightly different.
    SET_SOURCE_FILES_PROPERTIES(CFML_ILL_Instrm_Data.f90 PROPERTIES COMPILE_FLAGS "")

# Case of Lahey compiler.
ELSEIF(${COMPILER_NAME} STREQUAL lf95)

    # The base output path for the libraries.
    SET(LIB_PATH ${PROJECT_SOURCE_DIR}/Lahey)

    IF(${WINDOW})

        # An include directive for Realwin is added in case of a build in WINDOW mode.
        IF(${RWIN})

            SET_TARGET_PROPERTIES(wcrysfml PROPERTIES COMPILE_FLAGS "-mod .;$ENV{RWIN}")

        # An include directive for Winteracter is added in case of a build in WINDOW mode.
        ELSE()

            SET_TARGET_PROPERTIES(wcrysfml PROPERTIES COMPILE_FLAGS "-mod .;$ENV{WINTER}/lib.l95")

        ENDIF()

    ENDIF()

ENDIF()
