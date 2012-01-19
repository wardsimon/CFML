include(LibFindMacros)

include(set_winteracter_paths)
set_winteracter_paths(WINTERACTER_PATHS)

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

if(CMAKE_SIZEOF_VOID_P MATCHES 4)
    set(ARCH64 False)
else()
    set(ARCH64 True)
endif()

set(WINTERACTER_LIBRARY_PATHS)
set(WINTERACTER_INCLUDE_PATHS)
        
foreach(f ${WINTERACTER_PATHS})
    set(WINTERACTER_INCLUDE_PATHS ${WINTERACTER_INCLUDE_PATHS} ${f}/include)
    if(COMPILER_NAME STREQUAL ifort)
        if(${ARCH64})
            set(WINTERACTER_LIBRARY_PATHS ${WINTERACTER_LIBRARY_PATHS} ${f}/lib.i64)
        else()
            set(WINTERACTER_LIBRARY_PATHS ${WINTERACTER_LIBRARY_PATHS} ${f}/lib.if8)
        endif()
    elseif(COMPILER_NAME STREQUAL g95)        
        if(${ARCH64})
            set(WINTERACTER_LIBRARY_PATHS ${WINTERACTER_LIBRARY_PATHS} ${f}/lib.g95)
        else()
            set(WINTERACTER_LIBRARY_PATHS ${WINTERACTER_LIBRARY_PATHS} ${f}/lib.g95)
        endif()        
    endif()
endforeach()
                
find_path(WINTERACTER_INCLUDE_DIR
          NAMES winparam.h
          PATHS ${WINTERACTER_INCLUDE_PATHS})

find_path(WINTERACTER_MOD_DIR
          NAMES winteracter.mod
          PATHS ${WINTERACTER_LIBRARY_PATHS})

find_library(WINTERACTER_LIBRARY
             NAMES winter wint
             PATHS ${WINTERACTER_LIBRARY_PATHS})
                          
set(WINTERACTER_PROCESS_INCLUDES WINTERACTER_INCLUDE_DIR)

set(WINTERACTER_PROCESS_MODS WINTERACTER_MOD_DIR)

set(WINTERACTER_PROCESS_LIBS WINTERACTER_LIBRARY)

libfind_process(WINTERACTER)