#####################################################
# Infos
#####################################################

# CMake file for crysfml project
# Written by Eric Pellegrini
# Institut Laue Langevin
# Grenoble, France

# When running cmake from the command line, it is possible to specify command line options to cmake 
# that will set values in the cache. This is done with a -DVARIABLE:TYPE=VALUE syntax on the command 
# line. This is useful for non-interactive nightly test builds.

# What is the CMake cache?
# The cache is best thought of as a configuration file. Indeed Unix users could consider the cache as 
# equivalent to the set of flags passed to the configure command. The first time CMake is run, it 
# produces a CMakeCache.txt file. This file contains things like the existence and location of native 
# JPEG library. The entries are added in response to certain CMake commands (e.g. FIND_LIBRARY) as they 
# are processed anywhere in CMakeLists files anywhere in the source tree. After CMake has been run, and 
# created a CMakeCache.txt file - you may edit it. The CMake GUI, will allow you to edit the options 
# easily, or you can edit the file directly. The main reason for editing the cache would be to give 
# CMake the location of a native library such as JPEG, or to stop it from using a native library and use 
# a version of the library in your source tree.

# CMake will not alter an existing entry in the cache file itself. If your CMakeLists.txt files change 
# significantly, you will need to remove the relevant entries from the cache file. If you have not already 
# hand-edited the cache file, you could just delete it before re-running CMake. 

# From http://www.cmake.org/cmake/help/runningcmake.html

#####################################################



#####################################################
# Some macros
#####################################################

MACRO(ADD_PREFIX prefix rootlist)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${prefix}${root})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_PREFIX)

MACRO(ADD_SUFFIX rootlist suffix)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${root}${suffix})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_SUFFIX)

#####################################################



# The general compiling options for all compilers will be loaded the following file.
SET(CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran ${CMAKE_SOURCE_DIR}/Compilers/CompilerFlags.cmake)

# You must use a CMake >= 2.8.4 to build this project.
# Indeed, with prior versions, the resource files were skipped from the compilation silently.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4 FATAL_ERROR)

# By default rwin library is not used.
SET(RWIN FALSE CACHE BOOL "Use rwin library (only with lf95 compiler)")

# By default the built is in CONSOLE mode.
SET(GUI FALSE CACHE BOOL "do the build in GUI mode")

# The project is declared. It is a Fortran project.
PROJECT(crysfml Fortran)

# Switch off the verbosity.
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# Some checkings.
IF (${GUI})

    IF(WIN32)        
        IF(NOT (${COMPILER_NAME} STREQUAL "ifort" OR ${COMPILER_NAME} STREQUAL "lf95"))        
            MESSAGE(FATAL_ERROR "The crysfml library can not be built in GUI mode using ${COMPILER_NAME} 
compiler on ${CMAKE_HOST_SYSTEM} platform.")            
        ENDIF()
    
    ENDIF()

ENDIF()

# Add the subdirectories making the project.
ADD_SUBDIRECTORY(Src)


