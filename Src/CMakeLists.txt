#################################
# Include section
#################################

if(USE_HDF)
   add_definitions(-DUSE_HDF)
   find_package(hdf5_crysfml)
endif()

# Build crysfml_common library.
include(crysfml_common.cmake)

# Build crysfml library.
include(crysfml.cmake)

# If GUI option is True, build also wcrysfml library.
if(GUI)
    include(wcrysfml.cmake)
endif()
