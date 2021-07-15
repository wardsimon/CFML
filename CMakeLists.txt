#################################
# Commad line arguments section
#################################

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# By default the built is in CONSOLE mode and 32 bits. If 64 bits have to be used one should do -D ARCH32=OFF
option(MINGW_HDF_TRICK "Help building with MinGW" OFF)
option(USE_HDF "Build crysfml with HDF5 library" OFF)
option(GUI "Build crysfml and wcrysfml library" OFF)
option(ARCH32 "32 or 64 bit architecture" ON)
option(CRYSFML08 "Build CrysFML08" OFF)
option(HEAP_ARRAYS "Put arrays in heap instead of stack (only for Windows ifort)" OFF)
option(PROG_EX "Build the Program Examples" ON)
option(PYTHON_API "Build Python API" OFF)
option(QPARALLEL "qparallel option for Ifort on Windows" ON)
option(PYSETUP "Building from setup.py" OFF)

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

# Include and library paths for Python. 
# They must be set by the user.
# if(DEFINED PYTHON_LIBRARY_PATH)
#     set(PYTHON_LIBRARY_PATH ${PYTHON_LIBRARY_PATH} CACHE PATH "Python library path")
# else()
#     set(PYTHON_LIBRARY_PATH "" CACHE PATH "Python library path")
# endif()
# if(DEFINED PYTHON_INTERPRETER_PATH)
#     set(PYTHON_INTERPRETER_PATH ${PYTHON_INTERPRETER_PATH} CACHE PATH "Python interpreter path")
# else()
#     set(PYTHON_INTERPRETER_PATH "" CACHE PATH "Python interpreter path")
# endif()
# Use CMake v3.12+ FindPython3 to find the system-compatible interpreter *and* libraries at once
# This solves a potential issue on convoluted Windows build systems which have multiple
# python interpreters, some of which do not have build-system compatibility (I'm looking at you, msys2 python)
# In such a case FindPythonInterp might point to the wrong python interpreter and then further configuration
# steps will fail to find any (or any compatible) libraries which match the interpreter.
# find_package(Python3 COMPONENTS Interpreter Development)
# # Since FindPython3 and FindPythonInterp set different variables, fake the FindPythonInterp result
# # to prevent it running if called by, e.g., FindPybind11
# if(Python3_FOUND)
#   set(PYTHONINTERP_FOUND ON)
#   set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
#   message(STATUS "USING PYTHON: " ${PYTHON_EXECUTABLE})
# #   set(PYTHON_INTERPRETER_PATH ${Python3_EXECUTABLE})
#   set(PYTHON_VERSION_STRING "${Python3_VERSION}")
#   set(PYTHON_VERSION_MAJOR ${Python3_VERSION_MAJOR})
#   set(PYTHON_VERSION_MINOR ${Python3_VERSION_MINOR})
#   set(PYTHON_VERSION_PATCH ${Python3_VERSION_PATCH})
# endif()

# Include and library paths for HDF5.
# They must be set by the user.
if(DEFINED HDF5_INCLUDE_PATH)
    set(HDF5_INCLUDE_PATH ${HDF5_INCLUDE_PATH} CACHE PATH "HDF5 include path")
else()
    set(HDF5_INCLUDE_PATH "" CACHE PATH "HDF5 include path")
endif()

if(DEFINED HDF5_LIBRARY_PATH)
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
else()
    include(set_compiler_flags_08)
    set_compiler_flags_08()
    # The subdirectories of the project.
    add_subdirectory(Src08)
endif()

# Build the program example executables and use them as functional tests
if (PROG_EX)
  enable_testing()

  macro(add_runtest _name _bindir)
    add_test(
      NAME ${_name}
      COMMAND python3 ${CMAKE_SOURCE_DIR}/Tests/${_name}/test --binary-dir=${_bindir} --work-dir=${CMAKE_BINARY_DIR}/Tests/${_name} --verbose --log=${CMAKE_BINARY_DIR}/Tests/${_name}/runtest.stderr.log)
     #if(NOT "${_labels}" STREQUAL "")
     set_tests_properties(${_name} PROPERTIES
       RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Tests/${_name})
     #endif()
   endmacro()
  
  add_subdirectory(Program_Examples/PowderPattern)
  if(USE_HDF)
    add_subdirectory(Program_Examples/ReadNeXuS)
  endif()
endif()

if (PYTHON_API)
    set(PYTHON_API_PREFIX ${CMAKE_INSTALL_PREFIX}/Python_API)
    add_subdirectory(Python_API/Src)
endif()
