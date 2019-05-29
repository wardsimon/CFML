#################################
# Include section
#################################

add_definitions(-DUSE_HDF=ON)

# Build crysfml_common library.
include(crysfml_common.cmake)

# Build crysfml library.
include(crysfml.cmake)

if(USE_HDF)
  add_compile_options(-D USE_HDF)
endif()

# If GUI option is True, build also wcrysfml library.
if(GUI)
    include(wcrysfml.cmake)
endif()

