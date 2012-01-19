# The minimum version number insuring a correct build.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

option(GUI "Build crysfml AND wcrysfml libraries" ON)

# The cfml Fortran project.
project(cfml Fortran)

#################################
# Include section
#################################

# Add some paths to the cmake module search path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

# Build crysfml_common library.
include(crysfml_common)

# Build crysfml library.
include(crysfml)

# If GUI option is True, build also wcrysfml library.
if(GUI)
    include(wcrysfml)
endif()

