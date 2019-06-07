#################################
# Command line arguments section
#################################

#cmake -DCMAKE_GENERATOR_PLATFORM=x64 ..
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# A "trick" has been set on the build server to build with MinGW
# See See http://hdf-forum.184993.n3.nabble.com/HDF5-and-MinGW-td3393676.html for more details about this "trick"
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
# Use proper "Findboost_cryfml.cmake"
find_package(boost_crysfml)
include_directories(${Boost_INCLUDE_DIR})

# HDF5
# Use proper "Findhdf5_cryfml.cmake"
find_package(hdf5_crysfml)
include_directories(${HDF5_INCLUDE_DIR})

# Nexus and blosc
find_package(nexus)
include_directories(${NEXUS_INCLUDE_DIR})

# Including blosc, useless on macOS and Linux
find_package(blosc)
include_directories(${BLOSC_INCLUDE_DIR})