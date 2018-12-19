#################################
# Include section
#################################


# Build crysfml_common library.
include(crysfml_common.cmake)

# Build crysfml library.
include(crysfml.cmake)

#add_compile_options(-D USE_HDF -cpp)

# If GUI option is True, build also wcrysfml library.
if(GUI)
    include(wcrysfml.cmake)
endif()

