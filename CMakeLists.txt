# The first thing cmake does is to read, if one is present, the CMakeCache.txt file that will 
# contain the variables used for the Make process.

# If you modify one entry of the CMakeCache.txt this will result in a full rebuild of the project 
# at the next make.

SET(CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran ${CMAKE_SOURCE_DIR}/Compilers/CompilerFlags.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Set the project name and its associated language.
PROJECT(crysfml Fortran)

SET(CMAKE_VERBOSE_MAKEFILE OFF)

ADD_SUBDIRECTORY(Src)


