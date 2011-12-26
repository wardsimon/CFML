get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_Fortran_FLAGS "")

set(CMAKE_VERBOSE_MAKEFILE ON)

set(LIB_PATH $ENV{CRYSFML}/${COMPILER_NAME})

set(CRYSFML_COMMON_MOD_DIR ${LIB_PATH}/Common)
set(CRYSFML_COMMON_LIB_DIR ${LIB_PATH}/Common)

set(CRYSFML_MOD_DIR ${LIB_PATH}/LibC)
set(CRYSFML_LIB_DIR ${LIB_PATH}/LibC)

set(CRYSFGL_MOD_DIR ${LIB_PATH}/LibGL)
set(CRYSFGL_LIB_DIR ${LIB_PATH}/LibGL)

set(WCRYSFML_MOD_DIR ${LIB_PATH}/LibW)
set(WCRYSFML_LIB_DIR ${LIB_PATH}/LibW)

set(WINTER_INCLUDE_DIR $ENV{WINTER}/include)

if (IS_DIRECTORY $ENV{WINTER}/lib.i64)
	set(WINTER_MOD_DIR $ENV{WINTER}/lib.i64/)
	set(WINTER_LIB_DIR $ENV{WINTER}/lib.i64/)
elseif (IS_DIRECTORY $ENV{WINTER}/lib.if8)
	set(WINTER_MOD_DIR $ENV{WINTER}/lib.if8/)
	set(WINTER_LIB_DIR $ENV{WINTER}/lib.if8/)
else()
	message(FATAL_ERROR "Winteracter directory not found")
endif()

message("Winteracter directory found : ${WINTER_MOD_DIR}")
