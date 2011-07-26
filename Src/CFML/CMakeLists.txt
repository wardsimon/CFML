# CMake file for crysfml project
# Written by Eric Pellegrini
# Institut Laue Langevin
# Grenoble, France

# WINDOWS
if(WIN32)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
        set(crysfml_COMMON_SRCS CFML_GlobalDeps_Windows_intel.f90)
    else()
        set(crysfml_COMMON_SRCS CFML_GlobalDeps_Windows.f90)
        # Lahey Fortran compiler.
        if(${COMPILER_NAME} STREQUAL lf95)
            set(crysfml_COMMON_SRCS ${crysfml_COMMON_SRCS} f2kcli.f90)
        endif()
    endif()
# MacOS
elseif(APPLE)
    set(crysfml_COMMON_SRCS CFML_GlobalDeps_MacOS.f90)
# Unix
elseif(UNIX)
    # Intel Fortran compiler
    if(${COMPILER_NAME} STREQUAL ifort)
        set(crysfml_COMMON_SRCS CFML_GlobalDeps_Linux_Intel.f90)
    else()
        set(crysfml_COMMON_SRCS CFML_GlobalDeps_Linux.f90)
    endif()
endif()

# The sources files for crysfml_common library.
set(crysfml_COMMON_SRCS
    ${crysfml_COMMON_SRCS}
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
    CFML_Sym_Table.f90)

# The sources files for crysfml library.
set(crysfml_CONSOLE_SRCS CFML_IO_Mess.f90 CFML_Optimization_SAn.f90)

# The sources files for wcrysfml library.
if(${GUI})
    set(crysfml_GUI_SRCS CFML_Optimization_SAn.f90)    
    # Case of The Lahey compiler with a RWIN build (-DRWIN=True).
    if(${COMPILER_NAME} STREQUAL lf95 AND ${RWIN})
        set(crysfml_GUI_SRCS ${crysfml_GUI_SRCS} CFML_IO_MessRW.f90)
    else()
        set(crysfml_GUI_SRCS ${crysfml_GUI_SRCS} CFML_IO_MessWin.f90)
    endif()
endif()

set(COMMON_OBJS_DIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/crysfml_common.dir)
set(COMMON_OBJS_FILES ${crysfml_COMMON_SRCS})
add_suffix(COMMON_OBJS_FILES ${CMAKE_Fortran_OUTPUT_EXTENSION})
add_prefix(${COMMON_OBJS_DIR}/ COMMON_OBJS_FILES)

# Those file are generated during the build.
set_source_files_properties(${COMMON_OBJS_FILES} PROPERTIES GENERATED true)

# Some compiler specific flags are defined for the library to build.
if(${COMPILER_NAME} STREQUAL gfortran)

    if(WIN32)
        # By default all the source file are compiled with "-std=f2003" flag.
        set_source_files_properties(${crysfml_COMMON_SRCS} ${crysfml_CONSOLE_SRCS} 
                                    PROPERTIES COMPILE_FLAGS "-std=f2003")

        # Excepted the following files for which the compilation flags is slightly different for various reasons.
        set_source_files_properties(CFML_Profile_TOF.f90 PROPERTIES COMPILE_FLAGS "")
        set_source_files_properties(CFML_Conf_Calc.f90 PROPERTIES COMPILE_FLAGS "")
        set_source_files_properties(CFML_Optimization_SAn.f90 PROPERTIES COMPILE_FLAGS "")
        set_source_files_properties(CFML_ILL_Instrm_Data.f90 PROPERTIES COMPILE_FLAGS "-std=gnu")        
    endif()

# Case of g95 compiler.
elseif(${COMPILER_NAME} STREQUAL g95)

    if(WIN32)
        # By default all the source file are compiled with "-std=f2003" flag.
        set_source_files_properties(${crysfml_COMMON_SRCS} ${crysfml_CONSOLE_SRCS} ${crysfml_GUI_SRCS}
                                    PROPERTIES COMPILE_FLAGS "-std=f2003")

        # Excepted the following files for which the compilation flags is slightly different.
        set_source_files_properties(CFML_ILL_Instrm_Data.f90 PROPERTIES COMPILE_FLAGS "")
    endif()
        
endif()

# The crysfml_common library is common to Console and GUI mode.
add_library(crysfml_common STATIC ${crysfml_COMMON_SRCS})

# Sets the path where to place the crysfml_common library.
set_target_properties(crysfml_common PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Lib/Common)

# Sets the path where to place the mod files for the crysfml_common library.
set_target_properties(crysfml_common PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Lib/Common)

add_custom_command(TARGET crysfml_common
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${CRYSFML_COMMON_LIB_DIR}
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/Lib/Common ${CRYSFML_COMMON_LIB_DIR}
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${CRYSFML_LIB_DIR}
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/Lib/Common ${CRYSFML_LIB_DIR}
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${WCRYSFML_LIB_DIR}
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/Lib/Common ${WCRYSFML_LIB_DIR})

# This directory contains the crysfml_common library mod files.
include_directories(${PROJECT_BINARY_DIR}/Lib/Common)

# The linker path for crysfml_common library.
link_directories(${PROJECT_BINARY_DIR}/Lib/Common)

# The crysfml library is the CONSOLE version of the library. It is also built in GUI mode.
add_library(crysfml STATIC ${crysfml_CONSOLE_SRCS} ${COMMON_OBJS_FILES})

# Sets the path where to place the crysfml library.
set_target_properties(crysfml PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Lib/LibC)

# Sets the path where to place the mod files for the crysfml library.
set_target_properties(crysfml PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Lib/LibC)

# The crysfml_common will be linked to crysfml library.
target_link_libraries(crysfml crysfml_common)

add_custom_command(TARGET crysfml
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/Lib/LibC ${CRYSFML_LIB_DIR})

if(${GUI})

    # This directory contains the crysfml_common library mod files.
    include_directories(${PROJECT_BINARY_DIR}/Lib/Common)

    # This directory contains the winteracter library mod files.
    include_directories(${WINTER_MOD_DIR})

    # The linker path for crysfml_common library.
    link_directories(${PROJECT_BINARY_DIR}/Lib/Common)

    # The linker path for winteracter library.
    link_directories(${WINTER_LIB_DIR})
            
    # The wcrysfml library is the GUI version of the library. It is linked using Winteracter.
    add_library(wcrysfml STATIC ${crysfml_GUI_SRCS} ${COMMON_OBJS_FILES})

    target_link_libraries(wcrysfml crysfml_common)

    target_link_libraries(wcrysfml winter)
    
    # Sets the path where to place the wcrysfml library.
    set_target_properties(wcrysfml PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Lib/LibW)
    
    # Sets the path where to place the mod files for the wcrysfml library.
    set_target_properties(wcrysfml PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Lib/LibW)
        
    add_custom_command(TARGET wcrysfml
                       POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/Lib/LibW ${WCRYSFML_LIB_DIR})
    
endif()

#####################################################






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
