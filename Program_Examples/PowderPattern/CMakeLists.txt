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
target_link_libraries(calc_powder crysfml)

add_executable(simple_calc_powder Simple_calc_powder.f90)
target_link_libraries(simple_calc_powder crysfml)

add_executable(simple_calc_mag_powder Simple_calc_Mag_powder.f90)
target_link_libraries(simple_calc_mag_powder crysfml)


add_executable(laue_powder Laue_powder.f90)
target_link_libraries(laue_powder crysfml)

add_executable(corr_cells corr_cells.f90)
target_link_libraries(corr_cells crysfml)

### Unitests
add_runtest(calc_powder_STO ${CMAKE_CURRENT_BINARY_DIR})

# add_test(
#   NAME calc_powder
#   COMMAND calc_powder ${CMAKE_CURRENT_SOURCE_DIR}/SrTiO3.cfl)
# set_tests_properties(calc_powder PROPERTIES
#   RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Tests)



