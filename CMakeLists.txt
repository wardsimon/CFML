# The first thing cmake does is to read, if one is present, the CMakeCache.txt file that will 
# contain the variables used for the Make process.

# If you modify one entry of the CMakeCache.txt this will result in a full rebuild of the project 
# at the next make.

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
#####################################################


# The general compiling options for all compilers will be loaded the following file.
SET(CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran ${CMAKE_SOURCE_DIR}/Compilers/CompilerFlags.cmake)

# You must use a CMake >= 2.8 to build this project.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# By default rwin library is not used.
SET(RWIN FALSE CACHE BOOL "Use rwin library (only with lf95 compiler)")

# By default the built is in CONSOLE mode.
SET(GUI FALSE CACHE BOOL "do the build in GUI mode")

# The project name. It is a ForTran project.
PROJECT(crysfml Fortran)

# Switch off the verbosity.
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# Add the subdirectories making the project.
ADD_SUBDIRECTORY(Src)


