#################################
# Commad line arguments section
#################################

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

#################################
# Project section
#################################

# The project is declared. It is a Fortran project.
project(Python_API Fortran)
set(LIBRARY_NAME crysfml_api)
add_library(${LIBRARY_NAME} SHARED
  Extensions/API_init.f90
  Extensions/API_Atom_TypeDef.f90
  Extensions/API_Crystallographic_Symmetry.f90
  Extensions/API_Crystal_Metrics.f90
  Extensions/API_IO_Formats.f90
  Extensions/API_Reflections_Utilities.f90
  Extensions/API_Structure_Factors.f90
  Extensions/forpy_mod.f90)

# Forpy must be preprocessed
set_source_files_properties(forpy_mod.f90 PROPERTIES COMPILE_FLAG "-cpp")

#################################
# Dependencies section
#################################
# Add a dependency to crysfml_common to keep sure that the crysfml_common library will be built first.
add_dependencies(${LIBRARY_NAME} crysfml_common)

#find_package(Python3 REQUIRED)
#include_directories(${Python3_INCLUDE_DIRS})
if(NOT PYTHON_LIBRARY_PATH)
  #For a strange reason, find_package must be launched twice !
  find_package(Python3 COMPONENTS Interpreter Development)
  find_package(Python3 COMPONENTS Interpreter Development)
  set(PYTHON_LIBRARY_PATH ${Python3_LIBRARIES})
endif()

if (UNIX)
  add_library(Python3_LIB SHARED IMPORTED)
  set_property(TARGET Python3_LIB PROPERTY IMPORTED_LOCATION ${PYTHON_LIBRARY_PATH})
elseif(WIN32)
  add_library(Python3_LIB STATIC IMPORTED)
  set_property(TARGET Python3_LIB PROPERTY IMPORTED_LOCATION ${PYTHON_LIBRARY_PATH})
endif()

# Set include and linking
include_directories(${CRYSFML_COMMON_MODULE_DIRECTORY})
target_link_libraries(${LIBRARY_NAME} crysfml_common)
target_link_libraries(${LIBRARY_NAME} Python3_LIB)

#################################
# Install section
#################################
set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "" OUTPUT_NAME ${LIBRARY_NAME})
if(APPLE)
    set_target_properties(${LIBRARY_NAME} PROPERTIES SUFFIX ".so")
elseif(WIN32)
    set_target_properties(${LIBRARY_NAME} PROPERTIES SUFFIX ".pyd")
endif()

# The rules for installing the library.
install(DIRECTORY "" DESTINATION ${PYTHON_API_PREFIX}/CFML_api FILES_MATCHING PATTERN "*.py" PATTERN "Extensions" EXCLUDE)
if(WIN32)
install(TARGETS ${LIBRARY_NAME} RUNTIME DESTINATION ${PYTHON_API_PREFIX}/CFML_api)
else()
install(TARGETS ${LIBRARY_NAME} LIBRARY DESTINATION ${PYTHON_API_PREFIX}/CFML_api)
endif()
# Installing example files
install(DIRECTORY "../Examples" DESTINATION ${PYTHON_API_PREFIX} FILES_MATCHING PATTERN "*.py")
install(DIRECTORY "../Examples/Data" DESTINATION ${PYTHON_API_PREFIX}/Examples FILES_MATCHING PATTERN "*.cfl" PATTERN "*.cif")
