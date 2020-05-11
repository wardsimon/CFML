#################################
# Commad line arguments section
#################################

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

#################################
# Project section
#################################

# The project is declared. It is a Fortran project.
project(calc_powder Fortran)

include_directories(${CRYSFML_COMMON_MODULE_DIRECTORY})
include_directories(${CRYSFML_MODULE_DIRECTORY})

add_executable(calc_powder calc_powder.f90)

target_link_libraries(calc_powder
  crysfml
  crysfml_common)


