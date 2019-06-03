#################################
# Sources section
#################################

project(crysfml_cpp)

set(crysfml_cpp_sources ./NexusToCFML.cpp)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/")
message(${CMAKE_MODULE_PATH})

add_library(crysfml_cpp STATIC ./NexusToCFML.cpp ./DataNexusLib.cpp ./blosc_filter.c ./HeaderBloc.cpp  ./IntegerBloc.cpp  ./FloatBloc.cpp ./DataBloc.cpp ./ArgumentHandler.cpp )

install(TARGETS crysfml_cpp 
   DESTINATION ${PROJECT_SOURCE_DIR}/../ifort64/LibC)
   
install(TARGETS crysfml_cpp 
   DESTINATION ${PROJECT_SOURCE_DIR}/../ifort64HDF/LibC)   


#################################
# Command line arguments section
#################################

#cmake -DCMAKE_GENERATOR_PLATFORM=x64 ..
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

if(MINGW_HDF_TRICK)
  add_compile_options(-D MINGW_HDF_TRICK)
endif()

# Boost
#set(Boost_USE_STATIC_LIBS OFF) 
#set(Boost_USE_MULTITHREADED ON)  
#set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0)
if(WIN32)
set(Boost_INCLUDE_DIRS "C:\\Projects\\fullprof\\boost\\include\\boost-1_69")
endif()
include_directories(${Boost_INCLUDE_DIRS}) 

# HDF5
find_package(HDF5)
if(WIN32)
set(HDF5_INCLUDE_DIRS "C:\\Projects\\fullprof\\hdf5\\include\\")
endif()
include_directories(${HDF5_INCLUDE_DIRS})

# Nexus and blosc
find_package(nexus)
include_directories(${NEXUS_INCLUDE_DIR})
message(${NEXUS_INCLUDE_DIR})
# Including blosc, useless on macOS and Linux
find_package(blosc)
include_directories(${BLOSC_INCLUDE_DIR})


##########################################

#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "C:\\opt\\boost\\include\\boost-1_69")
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "c:\\opt\\boost\\lib")#"C:\\opt\\boost\\lib")
#file(GLOB HDF5_LIBRARIES ${CMAKE_SOURCE_DIR}/../opt/local_x64/lib/*hdf*.lib ${CMAKE_SOURCE_DIR}/../opt/local_x64/lib/*szip*.lib  ${CMAKE_SOURCE_DIR}/../opt/local_x64/lib/*zlib*.lib)
#include_directories(externals/nexus/include)
#add_subdirectory(../externals/c-blosc)
#add_subdirectory(../externals/nexus)

##file(GLOB Boost_LIBRARIES ${CMAKE_SOURCE_DIR}/../opt/local_x64/boost/lib64-msvc-14.1/lib*gd*.lib)
##file(GLOB Boost_LIBRARIES ${CMAKE_SOURCE_DIR}/../opt/boost/lib/*gd*x64*.lib)
##file(GLOB Boost_LIBRARIES "c:\\opt\\boost\\lib\\*gd*x64*.lib")

#include_directories(SrcCpp)
#add_subdirectory(SrcCpp)
#target_link_libraries(main ${CMAKE_SOURCE_DIR}/SrcCpp/crysfml_cpp.lib ${CMAKE_SOURCE_DIR}/externals/c-blosc/blosc/Debug.blosc.lib ${CMAKE_SOURCE_DIR}/externals/nexus/src/libNeXus.a ${HDF5_LIBRARIES} C:\\local\\boost_1_61_0\\lib64-msvc-14.0\\libboost_date_time-vc140-mt-gd-1_61.lib)
#add_executable(main main.cpp)
#include_directories(SrcCpp)
#include_directories("C:\\opt\\local_x64\\boost")
#include_directories("C:\\opt\\boost\\include\\boost-1_69")
#include_directories("C:\\opt\\local_x64\\include")

#target_link_libraries(main ${CMAKE_SOURCE_DIR}/build/SrcCpp/Debug/crysfml_cpp.lib)
#target_link_libraries(main ${CMAKE_SOURCE_DIR}/externals/c-blosc/blosc/Debug/blosc.lib)
#target_link_libraries(main ${CMAKE_SOURCE_DIR}/externals/nexus/Src/Debug/NeXus.lib)
#target_link_libraries(main ${HDF5_LIBRARIES})
#target_link_libraries(main ${Boost_LIBRARIES})
#target_link_libraries(main "${CMAKE_SOURCE_DIR}/../Program Files (x86)/Microsoft Visual Studio/2017/Community/VC/Tools/MSVC/14.16.27023/lib/x64/libcpmt.lib")

#add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MTd>)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MT")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MT")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE") #/NODEFAULTLIBS")