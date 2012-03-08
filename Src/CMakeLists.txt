#################################
# Include section
#################################

# Build crysfml_common library.
include(crysfml_common.cmake)

# Build crysfml library.
include(crysfml.cmake)

# If GUI option is True, build also wcrysfml library.
if(GUI)
    include(wcrysfml.cmake)
endif()

