#####################################################
# Define the source files for crysfml library       #
#####################################################

# This will store the source files common to a build in CONSOLE and GUI mode.
SET(crysfml_COMMON_LIB_SRCS "")


# WINDOWS
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

# Those files are in common for a build in Console or in GUI mode.
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

# Files specific to a build in CONSOLE mode.
SET(crysfml_CONSOLE_LIB_SRCS CFML_IO_Mess.f90 CFML_Optimization_SAn.f90)

# Files specific to a build in GUI mode.
IF(${GUI})

    SET(crysfml_GUI_LIB_SRCS CFML_Optimization_SAn.f90)
    
    # Case of The Lahey compiler with a RWIN build (-DRWIN=True).
    IF(${COMPILER_NAME} STREQUAL lf95 AND ${RWIN})

        SET(crysfml_WINDOW_LIB_SRCS ${crysfml_WINDOW_LIB_SRCS} CFML_IO_MessRW.f90)

    ELSE()

        SET(crysfml_WINDOW_LIB_SRCS ${crysfml_WINDOW_LIB_SRCS} CFML_IO_MessWin.f90)

    ENDIF()

ENDIF()

#####################################################



#####################################################
# Define sources properties                    #
#####################################################

SET(COMMON_OBJS_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/crysfml_common.dir)
SET(COMMON_OBJS_FILES ${crysfml_COMMON_LIB_SRCS})
ADD_SUFFIX(COMMON_OBJS_FILES ${CMAKE_Fortran_OUTPUT_EXTENSION})
ADD_PREFIX(${COMMON_OBJS_DIR}/ COMMON_OBJS_FILES)

# Those file are generated during the build.
SET_SOURCE_FILES_PROPERTIES(${COMMON_OBJS_FILES} PROPERTIES GENERATED true)

#####################################################



#####################################################
# Define targets properties                         #
#####################################################

# The crysfml_common library is common to Console and GUI mode.
ADD_LIBRARY(crysfml_common STATIC ${crysfml_COMMON_LIB_SRCS})

# The crysfml library is the CONSOLE version of the library. It is also built in GUI mode.
ADD_LIBRARY(crysfml STATIC ${crysfml_CONSOLE_LIB_SRCS} ${COMMON_OBJS_FILES})

# Sets the path where to place the crysfml_common library.
SET_TARGET_PROPERTIES(crysfml_common PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH}/Common)

# Sets the path where to place the mod files for the crysfml_common library.
SET_TARGET_PROPERTIES(crysfml_common PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_PATH}/LibC)

# Sets the path where to place the crysfml library.
SET_TARGET_PROPERTIES(crysfml PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH}/LibC)

# Sets the path where to place the mod files for the crysfml library.
SET_TARGET_PROPERTIES(crysfml PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_PATH}/LibC)

# The crysfml_common will be linked to crysfml library.
TARGET_LINK_LIBRARIES(crysfml crysfml_common)

if(${GUI})

    INCLUDE_DIRECTORIES(${LIB_PATH}/LibC)
            
    # The wcrysfml library is the GUI version of the library. It is linked using Winteracter.
    ADD_LIBRARY(wcrysfml STATIC ${crysfml_WINDOW_LIB_SRCS} ${COMMON_OBJS_FILES})

    IF(WIN32)
        SET_TARGET_PROPERTIES(wcrysfml PROPERTIES COMPILE_FLAGS "/I$ENV{WINTER}/lib.if8")        
    ELSEIF(APPLE)
        SET_TARGET_PROPERTIES(wcrysfml PROPERTIES COMPILE_FLAGS "-I$ENV{WINTER}/lib.ifi")        
    ELSEIF(UNIX)
        SET_TARGET_PROPERTIES(wcrysfml PROPERTIES COMPILE_FLAGS "-I$ENV{WINTER}/lib.if8")    
    ENDIF()


    # Sets the path where to place the mod files for the wcrysfml library.
    SET_TARGET_PROPERTIES(wcrysfml PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH}/LibW)
    
    # Sets the path where to place the mod files for the wcrysfml library.
    SET_TARGET_PROPERTIES(wcrysfml PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_PATH}/LibW)
        
    # The crysfml_common will be linked to wcrysfml library.
    TARGET_LINK_LIBRARIES(wcrysfml crysfml_common)
    
ENDIF()

#####################################################




# Some compiler specific flags are defined for the library to build.
# Case of Intel fortran compiler.
#IF(${COMPILER_NAME} STREQUAL ifort)

# Case of gfortran compiler.
#ELSEIF(${COMPILER_NAME} STREQUAL gfortran)


    # By default all the source file are compiled with "-std=f2003" flag.
#    SET_SOURCE_FILES_PROPERTIES(${crysfml_LIB_SRCS} PROPERTIES COMPILE_FLAGS "-std=f2003")

    # Excepted the following files for which the compilation flags is slightly different.
#    SET_SOURCE_FILES_PROPERTIES(CFML_Profile_TOF.f90 PROPERTIES COMPILE_FLAGS "")
#    SET_SOURCE_FILES_PROPERTIES(CFML_Conf_Calc.f90 PROPERTIES COMPILE_FLAGS "")
#    SET_SOURCE_FILES_PROPERTIES(CFML_Optimization_SAn.f90 PROPERTIES COMPILE_FLAGS "")
#    SET_SOURCE_FILES_PROPERTIES(CFML_ILL_Instrm_Data.f90 PROPERTIES COMPILE_FLAGS "-std=gnu")

# Case of g95 compiler.
#ELSEIF(${COMPILER_NAME} STREQUAL g95)

    # By default all the source file are compiled with "-std=f2003" flag.
#    SET_SOURCE_FILES_PROPERTIES(${crysfml_LIB_SRCS} PROPERTIES COMPILE_FLAGS "-std=f2003")

    # Excepted the following files for which the compilation flags is slightly different.
#    SET_SOURCE_FILES_PROPERTIES(CFML_ILL_Instrm_Data.f90 PROPERTIES COMPILE_FLAGS "")

# Case of Lahey compiler.
#ELSEIF(${COMPILER_NAME} STREQUAL lf95)


#    IF(${GUI})

#        # An include directive for Realwin is added in case of a build in GUI mode.
#        IF(${RWIN})

#            SET_TARGET_PROPERTIES(wcrysfml PROPERTIES COMPILE_FLAGS "-mod .;$ENV{RWIN}")

#        # An include directive for Winteracter is added in case of a build in GUI mode.
#        ELSE()

#            SET_TARGET_PROPERTIES(wcrysfml PROPERTIES COMPILE_FLAGS "-mod .;$ENV{WINTER}/lib.l95")

#        ENDIF()

#    ENDIF()

#ENDIF()
