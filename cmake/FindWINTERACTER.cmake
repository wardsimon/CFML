include(LibFindMacros)

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

if (WIN32)

    find_path(WINTERACTER_INCLUDE_DIR
              NAMES winparam.h
              PATHS $ENV{WINTERACTER}/include
                    $ENV{WINTER}/include
                    $ENV{WINT}/include
                    $ENV{USERPROFILE}/wint/include
                    $ENV{HOMEDRIVE}/wint/include
                    $ENV{SYSTEMDRIVE}/wint/include)

    if(COMPILER_NAME STREQUAL ifort)

        find_path(WINTERACTER_MOD_DIR
                  NAMES winteracter.mod
                  PATHS $ENV{USERPROFILE}/wint/lib.i64
                        $ENV{USERPROFILE}/wint/lib.if8
                        $ENV{HOMEDRIVE}/wint/lib.i64
                        $ENV{HOMEDRIVE}/wint/lib.if8
                        $ENV{SYSTEMDRIVE}/wint/lib.i64
                        $ENV{SYSTEMDRIVE}/wint/lib.if8
                        $ENV{WINTERACTER}/lib.i64
                        $ENV{WINTERACTER}/lib.if8
                        $ENV{WINTER}/lib.i64
                        $ENV{WINTER}/lib.if8
                        $ENV{WINT}/lib.i64
                        $ENV{WINT}/lib.if8)

        find_library(WINTERACTER_LIBRARY
                     NAMES winter wint
                     PATHS $ENV{WINTERACTER}/lib.i64
                           $ENV{WINTERACTER}/lib.if8
                           $ENV{WINTER}/lib.i64
                           $ENV{WINTER}/lib.if8
                           $ENV{WINT}/lib.i64
                           $ENV{WINT}/lib.if8
                           $ENV{USERPROFILE}/wint/lib.i64
                           $ENV{USERPROFILE}/wint/lib.if8
                           $ENV{HOMEDRIVE}/wint/lib.i64
                           $ENV{HOMEDRIVE}/wint/lib.if8
                           $ENV{SYSTEMDRIVE}/wint/lib.i64
                           $ENV{SYSTEMDRIVE}/wint/lib.if8)
    
    elseif(COMPILER_NAME STREQUAL g95)

        find_path(WINTERACTER_MOD_DIR
                  NAMES winteracter.mod
                  PATHS $ENV{WINTERACTER}/lib.g95
                        $ENV{WINTER}/lib.g95
                        $ENV{WINT}/lib.g95
                        $ENV{USERPROFILE}/wint/lib.g95
                        $ENV{HOMEDRIVE}/wint/lib.g95
                        $ENV{SYSTEMDRIVE}/wint/lib.g95)

        find_library(WINTERACTER_LIBRARY
                     NAMES winter wint
                     PATHS $ENV{WINTERACTER}/lib.g95
                           $ENV{WINTER}/lib.g95
                           $ENV{WINT}/lib.g95
                           $ENV{USERPROFILE}/wint/lib.g95
                           $ENV{HOMEDRIVE}/wint/lib.g95
                           $ENV{SYSTEMDRIVE}/wint/lib.g95)
    
    endif()        
    
else()

    find_path(WINTERACTER_INCLUDE_DIR
              NAMES winparam.h
              PATHS $ENV{WINTERACTER}/include
                    $ENV{WINTER}/include
                    $ENV{WINT}/include
                    $ENV{HOME}/wint/include
                    $ENV{HOME}/include
                    /usr/local/lib/wint/include
                    /usr/lib/wint/include
                    /opt/lib/wint/include)

    if(COMPILER_NAME STREQUAL ifort)

        find_path(WINTERACTER_MOD_DIR
                  NAMES winteracter.mod
                  PATHS $ENV{WINTERACTER}/lib.i64
                        $ENV{WINTERACTER}/lib.if8                  
                        $ENV{WINTER}/lib.i64
                        $ENV{WINTER}/lib.if8
                        $ENV{WINT}/lib.i64
                        $ENV{WINT}/lib.if8
                        $ENV{HOME}/wint/lib.i64
                        $ENV{HOME}/wint/lib.if8
                        $ENV{HOME}/lib.i64
                        $ENV{HOME}/lib.if8
                        /usr/local/lib/wint/lib.i64
                        /usr/local/lib/wint/lib.if8
                        /usr/lib/wint/lib.i64
                        /usr/lib/wint/lib.if8
                        /opt/lib/wint/lib.i64
                        /opt/lib/wint/lib.if8)

        find_library(WINTERACTER_LIBRARY
                     NAMES winter wint
                     PATHS $ENV{WINTERACTER}/lib.i64
                           $ENV{WINTERACTER}/lib.if8                  
                           $ENV{WINTER}/lib.i64
                           $ENV{WINTER}/lib.if8
                           $ENV{WINT}/lib.i64
                           $ENV{WINT}/lib.if8
                           $ENV{HOME}/wint/lib.i64
                           $ENV{HOME}/wint/lib.if8
                           $ENV{HOME}/lib.i64
                           $ENV{HOME}/lib.if8
                           /usr/local/lib/wint/lib.i64
                           /usr/local/lib/wint/lib.if8
                           /usr/lib/wint/lib.i64
                           /usr/lib/wint/lib.if8
                           /opt/lib/wint/lib.i64
                           /opt/lib/wint/lib.if8)

    elseif(COMPILER_NAME STREQUAL g95)

        find_path(WINTERACTER_MOD_DIR
                  NAMES winteracter.mod
                  PATHS $ENV{WINTERACTER}/lib.g95
                        $ENV{WINTER}/lib.g95
                        $ENV{WINT}/lib.g95
                        $ENV{HOME}/wint/lib.g95
                        $ENV{HOME}/lib.g95
                        /usr/local/lib/wint/lib.g95
                        /usr/lib/wint/lib.g95
                        /opt/lib/wint/lib.g95)

        find_library(WINTERACTER_LIBRARY
                     NAMES winter wint
                     PATHS $ENV{WINTERACTER}/lib.g95
                           $ENV{WINTER}/lib.g95
                           $ENV{WINT}/lib.g95
                           $ENV{HOME}/wint/lib.g95
                           $ENV{HOME}/lib.g95
                           /usr/local/lib/wint/lib.g95
                           /usr/lib/wint/lib.g95
                           /opt/lib/wint/lib.g95)
    
    endif()

endif()
                                                  
set(WINTERACTER_PROCESS_INCLUDES WINTERACTER_INCLUDE_DIR)

set(WINTERACTER_PROCESS_MODS WINTERACTER_MOD_DIR)

set(WINTERACTER_PROCESS_LIBS WINTERACTER_LIBRARY)

libfind_process(WINTERACTER)