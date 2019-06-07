#################################
# Command line arguments section
#################################

#cmake -DCMAKE_GENERATOR_PLATFORM=x64 ..
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

if(MINGW_HDF_TRICK)
  add_compile_options(-DMINGW_HDF_TRICK)
endif()

#################################
# Sources section
#################################

project(crysfml_cpp)

set(crysfml_cpp_sources ./NexusToCFML.cpp)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/")

add_library(crysfml_cpp STATIC ./NexusToCFML.cpp ./DataNexusLib.cpp ./blosc_filter.c ./HeaderBloc.cpp  ./IntegerBloc.cpp  ./FloatBloc.cpp ./DataBloc.cpp ./ArgumentHandler.cpp )

install(TARGETS crysfml_cpp 
   DESTINATION ${PROJECT_SOURCE_DIR}/../ifort64/LibC)
   
install(TARGETS crysfml_cpp 
   DESTINATION ${PROJECT_SOURCE_DIR}/../ifort64HDF/LibC)

#################################
# Dependencies section
#################################

# Boost
find_package(boost_crysfml)
include_directories(${Boost_INCLUDE_DIR})

# HDF5
find_package(hdf5_crysfml)
include_directories(${HDF5_INCLUDE_DIR})

# Nexus and blosc
find_package(nexus)
include_directories(${NEXUS_INCLUDE_DIR})

# Including blosc, useless on macOS and Linux
find_package(blosc)
include_directories(${BLOSC_INCLUDE_DIR})