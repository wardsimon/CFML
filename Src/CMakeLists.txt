#####################################################
# Infos
#####################################################

# CMake file for crysfml project
# Written by Eric Pellegrini
# Institut Laue Langevin
# Grenoble, France

ADD_SUBDIRECTORY(CFML)

IF(${GUI})
    ADD_SUBDIRECTORY(CFGL)
ENDIF()
