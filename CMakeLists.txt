#################################
# Project section
#################################

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# Set the a verbose make.
set(CMAKE_VERBOSE_MAKEFILE ON)

IF(DEFINED CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type as ${CMAKE_BUILD_TYPE}")
  SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
  message(STATUS "Setting build type as 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
ENDIF()

# By default the built is in CONSOLE mode.
option(GUI "Build crysfml and wcrysfml library" OFF)

# The project is declared. It is a Fortran project.
project(libcrysfml Fortran)

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

if (NOT CRYSFML_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/${COMPILER_NAME})
endif()

if (NOT CRYSFML_PREFIX)
    # Set the default name for crysfml library within the install prefix.
    set(CRYSFML_PREFIX LibC CACHE STRING "The name of crysfml library within library directory")
endif()

if (NOT WCRYSFML_PREFIX)
    # Set the default name for wcrysfml library within the install prefix.
    set(WCRYSFML_PREFIX LibW CACHE STRING "The name of wcrysfml library within library directory")
endif()

#################################
# Includes section
#################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Includes some modules.
include(add_prefix)
include(add_suffix)
include(set_compiler_flags)
set_compiler_flags()

#################################
# Project building
#################################

# The subdirectories of the project.
add_subdirectory(Src)
