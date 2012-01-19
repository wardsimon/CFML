cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# By default the built is in CONSOLE mode.
option(GUI "Build also crysfgl and wcrysfml library" OFF)

# The project is declared. It is a Fortran project.
project(libcrysfml Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CRYSFML_SOURCE}/cmake)

# Includes some modules. 
include(general)
include(add_prefix)
include(add_suffix)
include(set_compiler_flags)

set_compiler_flags()

# The subdirectories of the project.
add_subdirectory(Src)
