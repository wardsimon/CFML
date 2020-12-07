# Build crysfml library.
include(crysfml08_common.cmake)

if(USE_HDF)
   add_definitions(-DUSE_HDF)
   find_package(hdf5_crysfml)
endif()

# If GUI option is True, build also wcrysfml library.
if(GUI)
    include(wcrysfml08.cmake)
endif()