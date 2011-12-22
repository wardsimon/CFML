set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn")

set(CMAKE_Fortran_FLAGS_RELEASE "-O -warn -vec-report0")

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(WINTER_MOD_DIR $ENV{WINTER}/lib.i64)
    set(WINTER_LIB_DIR $ENV{WINTER}/lib.i64)
else()
    set(WINTER_MOD_DIR $ENV{WINTER}/lib.if8)
    set(WINTER_LIB_DIR $ENV{WINTER}/lib.if8)
endif()
