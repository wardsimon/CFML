cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# By default the built is in CONSOLE mode.
option(GUI "Build also crysfgl and wcrysfml library" OFF)

# Set the default name for crysfml library within the install prefix.
set(CRYSFML_PREFIX crysfml CACHE STRING "The name of crysfml library within library directory")

# Set the default name for wcrysfml library within the install prefix.
set(WCRYSFML_PREFIX wcrysfml CACHE STRING "The name of wcrysfml library within library directory")

# The project is declared. It is a Fortran project.
project(libcrysfml Fortran)

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

# Set the a verbose make.
set(CMAKE_VERBOSE_MAKEFILE ON)

# By default the build type is set to RELEASE.
set(CMAKE_BUILD_TYPE Release)

# Nullify all the Fortran flags.
set(CMAKE_Fortran_FLAGS "")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CRYSFML}/cmake)

# Includes some modules.
include(add_prefix)
include(add_suffix)
include(set_compiler_flags)

set_compiler_flags()

# The subdirectories of the project.
add_subdirectory(Src)
