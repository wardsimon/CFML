#################################
# Project section
#################################

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# By default the built is in CONSOLE mode.
option(GUI "Build crysfml and wcrysfml library" OFF)

# The project is declared. It is a Fortran project.
project(libcrysfml Fortran)

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

# Set the a verbose make.
set(CMAKE_VERBOSE_MAKEFILE ON)

if (NOT CRYSFML_PREFIX)
    # Set the default name for crysfml library within the install prefix.
    set(CRYSFML_PREFIX LibC CACHE STRING "The name of crysfml library within library directory")
endif()

if (NOT WCRYSFML_PREFIX)
    # Set the default name for wcrysfml library within the install prefix.
    set(WCRYSFML_PREFIX LibW CACHE STRING "The name of wcrysfml library within library directory")
endif()

# By default the build type is set to RELEASE.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
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
