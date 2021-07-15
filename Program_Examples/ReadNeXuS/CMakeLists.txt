#################################
# Commad line arguments section
#################################

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

#################################
# Project section
#################################

# The project is declared. It is a Fortran project.
project(read_nexus Fortran)
find_package(hdf5_crysfml)

include_directories(${CRYSFML_COMMON_MODULE_DIRECTORY})
include_directories(${CRYSFML_MODULE_DIRECTORY})

add_executable(read_nexus read_nexus.f90)
target_link_libraries(read_nexus crysfml)
link_directories(${HDF5_LIBRARY_DIR})
target_link_libraries(read_nexus ${HDF5_LIBRARIES})

if(WIN32 OR MSYS)
    get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)
    if(COMPILER_NAME STREQUAL ifort)
# Copy hdf5.dll hdf5_f90cstub.dll hdf5_fortran.dll
    else()
# Copy libgcc_s_seh-1.dll libgfortran-5.dll libhdf5_f90cstub-0.dll libhdf5_fortran-0.dll libhdf5-0.dll libquadmath-0.dll libszip.dl libwinpthread-1.dll zlib1.dll    
    endif()
endif()

add_test(
  NAME Test_NeXuS_reading
  COMMAND "${CMAKE_BINARY_DIR}/Program_Examples/ReadNeXuS/read_nexus" "${CMAKE_BINARY_DIR}/../Program_Examples/ReadNeXuS/400006" "D19")