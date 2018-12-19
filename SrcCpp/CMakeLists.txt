#################################
# Sources section
#################################

project(crysfml_cpp)

set(crysfml_cpp_sources ./NexusToCFML.cpp)

add_library(crysfml_cpp STATIC ./NexusToCFML.cpp ./DataNexusLib.cpp ./blosc_filter.c ./HeaderBloc.cpp  ./IntegerBloc.cpp  ./FloatBloc.cpp ./DataBloc.cpp ./ArgumentHandler.cpp )

install(TARGETS crysfml_cpp 
   DESTINATION ${PROJECT_SOURCE_DIR}/../ifort64/LibC)
   
install(TARGETS crysfml_cpp 
   DESTINATION ${PROJECT_SOURCE_DIR}/../ifort64HDF/LibC)   

# The rules for installing the mod files. Take care the "/" is on purpose.
#install(DIRECTORY ${CRYSFML_COMMON_MODULE_DIRECTORY}/ ${CRYSFML_MODULE_DIRECTORY}/
#        DESTINATION ${CRYSFML_PREFIX}
#        FILES_MATCHING
#        PATTERN "*.mod"
#        PATTERN CMakeFiles EXCLUDE)