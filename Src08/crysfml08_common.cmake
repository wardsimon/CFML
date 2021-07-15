# Set source files and compiler flags for each module

# CFML_GlobalDeps
# This module is OS and compiler dependent.
if(WIN32 OR MSYS)
   # Windows
    if(${COMPILER_NAME} STREQUAL ifort)
        # Intel Fortran compiler
        set(GLOBAL_DEPS_SRC CFML_GlobalDeps_Windows_IFOR.f90)
    else()
        # GFortran compiler
        set(GLOBAL_DEPS_SRC CFML_GlobalDeps_Windows_GFOR.f90)
    endif()
elseif(APPLE)
    # MacOS
    if(${COMPILER_NAME} STREQUAL ifort)
        # Intel Fortran compiler
        set(GLOBAL_DEPS_SRC CFML_GlobalDeps_MacOS_IFOR.f90)
    else()
        # GFortran compiler
        set(GLOBAL_DEPS_SRC CFML_GlobalDeps_MacOS_GFOR.f90)
    endif()
elseif(UNIX)
    # Unix
    if(${COMPILER_NAME} STREQUAL ifort)
        # Intel Fortran compiler
        set(GLOBAL_DEPS_SRC CFML_GlobalDeps_Linux_IFOR.f90)
    else()
        # GFortran compiler
        set(GLOBAL_DEPS_SRC CFML_GlobalDeps_Linux_GFOR.f90)
    endif()
endif()
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${GLOBAL_DEPS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${GLOBAL_DEPS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Messages
file(GLOB SUBMOD_MESSAGES_SRC CFML_Messages/*.f90)
set(MESSAGES_SRC CFML_Messages.f90
                 CFML_Messages/Con_Err_Message.f90
                 CFML_Messages/Con_Info_Message.f90
                 CFML_Messages/Con_Print_Message.f90
                 CFML_Messages/Con_Wait_Message.f90
                 CFML_Messages/Con_Write_ScrollMsg.f90)
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${MESSAGES_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${MESSAGES_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# Mathematics: CFML_Maths,CFML_FFT,CFML_Random
file(GLOB SUBMOD_MATHS_SRC CFML_Maths/*.f90)
file(GLOB SUBMOD_FFT_SRC CFML_FFT/*.f90)
file(GLOB SUBMOD_RANDOM_SRC CFML_Random/*.f90)
set(MATHS_SRC CFML_Maths.f90
              CFML_FFT.f90
              CFML_Random.f90
              # CFML_Trigonometry.f90
              ${SUBMOD_MATHS_SRC}
              ${SUBMOD_FFT_SRC}
              ${SUBMOD_RANDOM_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${MATHS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${MATHS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Strings
file(GLOB SUBMOD_STRINGS_SRC CFML_Strings/*.f90)
set(STRINGS_SRC CFML_Strings.f90
                ${SUBMOD_STRINGS_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${STRINGS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${STRINGS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Rational
file(GLOB SUBMOD_RATIONAL_SRC CFML_Rational/*.f90)
set(RATIONAL_SRC CFML_Rational.f90
                ${SUBMOD_RATIONAL_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${RATIONAL_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${RATIONAL_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Metrics
file(GLOB SUBMOD_METRICS_SRC CFML_Metrics/*.f90)
set(METRICS_SRC CFML_Metrics.f90
                ${SUBMOD_METRICS_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${METRICS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${METRICS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Tables
set(SUBMOD_TABLES_1_SRC CFML_Tables/Tab_Del_ScatterT.f90
                        CFML_Tables/Tab_Get_ScatterT.f90
                        CFML_Tables/Tab_Del_BondsT.f90
                        CFML_Tables/Tab_Get_BondsT.f90
                        CFML_Tables/Tab_Del_SpgT.f90
                        CFML_Tables/Tab_Get_SpgT.f90
                        CFML_Tables/Tab_Del_BVST.f90
                        CFML_Tables/Tab_Allocating_MagneticDBase.f90
                        CFML_Tables/Tab_Read_MagneticDBase.f90
                        CFML_Tables/Tab_Allocating_SuperSpaceDBase.f90
                        CFML_Tables/Tab_Read_SSG_DBase.f90)
set(SUBMOD_TABLES_2_SRC CFML_Tables/Tab_Set_ScatterT.f90
                        CFML_Tables/Tab_Set_BondsT.f90
                        CFML_Tables/Tab_Get_SpgSymbols.f90
                        CFML_Tables/Tab_Set_SpgT.f90
                        CFML_Tables/Tab_Set_BVST.f90)
set(TABLES_1_SRC CFML_Tables_Scattering.f90
                 CFML_Tables_Bonds.f90
                 ${SUBMOD_TABLES_1_SRC})
set(TABLES_2_SRC CFML_Tables_Symmetry.f90
                 CFML_Tables_BVS.f90
                 CFML_Tables_MagneticDB.f90
                 CFML_Tables_SuperSpaceDB.f90)
set(TABLES_3_SRC ${SUBMOD_TABLES_2_SRC})

if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${TABLES_1_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${TABLES_1_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${TABLES_2_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS0} ${OPT_FLAGS2}")

    set_source_files_properties(${TABLES_3_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS0}")
else()
    set_source_files_properties(${TABLES_2_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS0}")

    set_source_files_properties(${TABLES_3_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS0}")
endif()

# CFML_gSpaceGroups
file(GLOB SUBMOD_GROUPS_SRC CFML_gSpaceGroups/*.f90)
set(GROUPS_SRC CFML_gSpaceGroups.f90
               ${SUBMOD_GROUPS_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${GROUPS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${GROUPS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Profiles
file(GLOB SUBMOD_PROFILES_SRC CFML_Profiles/Prof*.f90)
set(PROFILES_1_SRC CFML_Profiles.f90
                   ${SUBMOD_PROFILES_SRC})
set(PROFILES_2_SRC CFML_Profiles/Profile_Init_ProfVal.f90)
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${PROFILES_1_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${PROFILES_1_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${PROFILES_2_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS0} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${PROFILES_2_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS0}")
endif()

# CFML_DiffPatt
file(GLOB SUBMOD_DIFFPATT_SRC CFML_DiffPatt/*.f90)
set(DIFFPATT_SRC CFML_Diffpatt.f90
                 ${SUBMOD_DIFFPATT_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${DIFFPATT_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${DIFFPATT_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_ExtintCorr
file(GLOB SUBMOD_EXTINCORR_SRC CFML_ExtinCorr/*.f90)
set(EXTINCORR_SRC CFML_ExtinCorr.f90
                  ${SUBMOD_EXTINCORR_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${EXTINCORR_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${EXTINCORR_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_EOS
file(GLOB SUBMOD_EOS_SRC CFML_EoS/*.f90)
set(EOS_SRC CFML_EoS.f90
            ${SUBMOD_EOS_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${EOS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${EOS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Atoms
file(GLOB SUBMOD_ATOMS_SRC CFML_Atoms/*.f90)
set(ATOMS_SRC CFML_Atoms.f90
              ${SUBMOD_ATOMS_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${ATOMS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${ATOMS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Reflections
file(GLOB SUBMOD_REFLECTIONS_SRC CFML_Reflections/*.f90)
set(REFLECTIONS_SRC CFML_Reflections.f90
                 ${SUBMOD_REFLECTIONS_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${REFLECTIONS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${REFLECTIONS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Propagk
set(PROPAGK_SRC CFML_Propagk.f90)
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${PROPAGK_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${PROPAGK_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_IOForm
file(GLOB SUBMOD_IOFORM_SRC CFML_IOForm/*.f90)
set(IOFORM_SRC CFML_IOForm.f90 ${SUBMOD_IOFORM_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${IOFORM_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${IOFORM_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Geom
file(GLOB SUBMOD_GEOM_SRC CFML_Geom/*.f90)
set(GEOM_SRC CFML_Geom.f90
                 ${SUBMOD_GEOM_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${GEOM_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${GEOM_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Maps
file(GLOB SUBMOD_MAPS_SRC CFML_Maps/*.f90)
set(MAPS_SRC CFML_Maps.f90
                 ${SUBMOD_MAPS_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${MAPS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${MAPS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Optimization
file(GLOB SUBMOD_OPT_SRC CFML_Optimization/*.f90)
set(OPT_SRC CFML_Optimization.f90
                 ${SUBMOD_OPT_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${OPT_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${OPT_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Optimization_LSQ
file(GLOB SUBMOD_OPT_LSQ_SRC CFML_Optimization_LSQ/*.f90)
set(OPT_LSQ_SRC CFML_Optimization_LSQ.f90
                 ${SUBMOD_OPT_LSQ_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${OPT_LSQ_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${OPT_LSQ_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Optimization_SAnn
file(GLOB SUBMOD_OPT_SAN_SRC CFML_Optimization_SAnn/*.f90)
set(OPT_SAN_SRC CFML_Optimization_SAnn.f90
                 ${SUBMOD_OPT_SAN_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${OPT_SAN_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${OPT_SAN_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_Export_VTK
set(Export_VTK_SRC CFML_Export_VTK.f90)
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${Export_VTK_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${Export_VTK_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()
# CFML_EnBVS
file(GLOB SUBMOD_EnBVS_SRC CFML_EnBVS/*.f90)
set(EnBVS_SRC CFML_EnBVS.f90
                 ${SUBMOD_EnBVS_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${EnBVS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${EnBVS_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# Instrm_ILL
set(InstrmILL_SRC CFML_ILL_Instrm_Data.f90)
 
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${InstrmILL_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${InstrmILL_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

# CFML_SXTAL_Geom
file(GLOB SUBMOD_SXTALgeom_SRC CFML_SXTAL_Geom/*.f90)
set(SXTALgeom_SRC CFML_SXTAL_Geom.f90
                 ${SUBMOD_SXTALgeom_SRC})
if(${COMPILER_NAME} STREQUAL ifort)
    set_source_files_properties(${SXTALgeom_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGS} ${OPT_FLAGS1} ${OPT_FLAGS2}")
else()
    set_source_files_properties(${SXTALgeom_SRC}
        PROPERTIES COMPILE_FLAGS "${OPT_FLAGSC} ${OPT_FLAGS1}")
endif()

#  List of all the source files 
set(CRYSFML_COMMON_SRC
    ${GLOBAL_DEPS_SRC}
    ${MESSAGES_SRC}
    ${MATHS_SRC}
    ${STRINGS_SRC}
    ${RATIONAL_SRC}
    ${METRICS_SRC}
    ${TABLES_1_SRC}
    ${TABLES_2_SRC}
    ${TABLES_3_SRC}
    ${GROUPS_SRC}
    ${PROFILES_1_SRC}
    ${PROFILES_2_SRC}
    ${DIFFPATT_SRC}
    ${EXTINCORR_SRC}
    ${EOS_SRC}
    ${ATOMS_SRC}
    ${REFLECTIONS_SRC}
    ${PROPAGK_SRC}
    ${IOFORM_SRC}
    ${GEOM_SRC}
    ${MAPS_SRC}
    ${OPT_SRC}
    ${OPT_LSQ_SRC}
    ${OPT_SAN_SRC} 
    ${Export_VTK_SRC} 
    ${EnBVS_SRC}
    ${InstrmILL_SRC}
    ${SXTALgeom_SRC})

# Build the library
set(LIBRARY_NAME crysfml)

add_library(${LIBRARY_NAME} STATIC ${CRYSFML_COMMON_SRC})

# The directory where the CrysFML modules files will be stored.
set(CRYSFML_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Src08/crysfml_modules)

# Sets the path where to place the mod files for the crysfml_common library.
set_target_properties(${LIBRARY_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CRYSFML_MODULE_DIRECTORY})

#################################
# Install section
#################################

# The rules for installing the library.
install(TARGETS ${LIBRARY_NAME} ARCHIVE DESTINATION ${CRYSFML_PREFIX})

# The rules for installing the mod files. Take care the "/" is on purpose.
install(DIRECTORY ${CRYSFML_MODULE_DIRECTORY}/
        DESTINATION ${CRYSFML_PREFIX}
        FILES_MATCHING
        PATTERN "*.*mod"
        PATTERN CMakeFiles EXCLUDE)