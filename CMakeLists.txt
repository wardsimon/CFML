cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# By default the built is in CONSOLE mode.
option(GUI "Build also crysfgl and wcrysfml library" OFF)
set(CRYSFML_PREFIX lib_crysfml CACHE STRING "The name of crysfml library within library directory")
set(WCRYSFML_PREFIX lib_wcrysfml CACHE STRING "The name of wcrysfml library within library directory")

# The project is declared. It is a Fortran project.
project(libcrysfml Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CRYSFML}/cmake)

# Includes some modules.
include(general)
include(add_prefix)
include(add_suffix)
include(set_compiler_flags)

set_compiler_flags()

# The subdirectories of the project.
add_subdirectory(Src)
