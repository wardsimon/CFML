# CMake file for crysfml project
# Written by Eric Pellegrini
# Institut Laue Langevin
# Grenoble, France

add_subdirectory(CFML)

if(${GUI})
    add_subdirectory(CFGL)
endif()
