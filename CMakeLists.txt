#################################
# Commad line arguments section
#################################

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# By default the built is in CONSOLE mode and 32 bits. If 64 bits have to be used one should do -D ARCH32=OFF
option(MINGW_HDF_TRICK "Help building with MinGW" OFF)
option(USE_HDF "Build crysfml and wcrysfml library" OFF)
option(GUI "Build crysfml and wcrysfml library" OFF)
option(ARCH32 "32 or 64 bit architecture" ON)
option(CRYSFML08 "Build CrysFML08" OFF)

if(DEFINED CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} CACHE STRING "The Fortran compiler. One of 'ifort', 'g95' or 'gfortran'.")
else()
    set(CMAKE_Fortran_COMPILER ifort CACHE STRING "The Fortran compiler. One of 'ifort', 'g95' or 'gfortran'.")
endif()
message(STATUS "Setting compiler to ${CMAKE_Fortran_COMPILER}")

if(DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "The build type. One of 'Release' or 'Debug'.")
else()
    set(CMAKE_BUILD_TYPE Release CACHE STRING "The build type. One of 'Release' or 'Debug'.")
endif()
message(STATUS "Setting build type to ${CMAKE_BUILD_TYPE}")

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

if (DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE PATH "The installation path.")
else()
    set(CMAKE_INSTALL_PREFIX ${CMAKE_HOME_DIRECTORY}/${COMPILER_NAME} CACHE PATH "The installation path.")
endif()
message(STATUS "Setting CMAKE_INSTALL_PREFIX to ${CMAKE_INSTALL_PREFIX}")

if (DEFINED CRYSFML_PREFIX)
    set(CRYSFML_PREFIX ${CRYSFML_PREFIX})
else()
    if(CRYSFML08)
        set(CRYSFML_PREFIX LibC08)
    else()
        set(CRYSFML_PREFIX LibC)
    endif()
endif()
set (CRYSFML_PREFIX ${CMAKE_INSTALL_PREFIX}/${CRYSFML_PREFIX})
message(STATUS "CrysFML installation path set to ${CRYSFML_PREFIX}")

# Include and library paths for HDF5. 
# They must be set by the user.
if(DEFINED HDF5_INCLUDE_PATH)   
    set(HDF5_INCLUDE_PATH ${HDF5_INCLUDE_PATH} CACHE PATH "HDF5 include path")
else()
    set(HDF5_INCLUDE_PATH "" CACHE PATH "HDF5 include path")
endif()

if(NOT DEFINED HDF5_LIBRARY_PATH)    
    set(HDF5_LIBRARY_PATH ${HDF5_LIBRARY_PATH} CACHE PATH "HDF5 library path")
else()
    set(HDF5_LIBRARY_PATH "" CACHE PATH "HDF5 library path")
endif()

#################################
# Project section
#################################

# The project is declared. It is a Fortran project.
project(libcrysfml Fortran)

# Set the a verbose make.
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

if (GUI)
    if (DEFINED WCRYSFML_PREFIX)
        set(WCRYSFML_PREFIX ${WCRYSFML_PREFIX})
    else()
        set(WCRYSFML_PREFIX LibW)
    endif()
    set (WCRYSFML_PREFIX ${CMAKE_INSTALL_PREFIX}/${WCRYSFML_PREFIX})
    message(STATUS "WCrysFML installation path set to ${WCRYSFML_PREFIX}")
endif()
if (NOT WCRYSFML_PREFIX)
    # Set the default name for wcrysfml library within the install prefix.
    set(WCRYSFML_PREFIX LibW CACHE STRING "The name of wcrysfml library within library directory")
endif()

# Includes some modules and build the project.
include(add_prefix)
include(add_suffix)    
if (NOT CRYSFML08)        
    include(set_compiler_flags)
    set_compiler_flags()
    # The subdirectories of the project.
    add_subdirectory(Src)
    #if (USE_HDF)
    #    add_subdirectory(SrcCpp)
    #endif()
else()
    include(set_compiler_flags_08)
    set_compiler_flags_08()    
    # The subdirectories of the project.
    add_subdirectory(Src08)
endif()
