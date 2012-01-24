macro(set_compiler_flags)
    
    get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)
    
    if(COMPILER_NAME STREQUAL ifort)
    
        if(WIN32)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-debug:full /check /traceback /nologo")
                set(OPT_FLAGS "/Od")
                set(OPT_FLAGS1 "/Od")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "/Qvec-report0 /nologo")
                set(OPT_FLAGS "/O2")
                set(OPT_FLAGS1 "/Od")
            endif()
        elseif(APPLE)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-warn")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-warn -vec-report0")
                set(OPT_FLAGS "-O")
                set(OPT_FLAGS1 "-O0")
            endif()
        else()
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-O -warn -vec-report0")
                set(OPT_FLAGS "-O")
                set(OPT_FLAGS1 "-O0")
            endif()    
        endif()

    elseif(COMPILER_NAME STREQUAL g95)

        if(WIN32)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-ftrace=full -std=f2003")
                set(OPT_FLAGS "-O0")
                set(OPT_FLAGS1 "-O0")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-std=f2003")
                set(OPT_FLAGS "-O3 -funroll-loops -msse2")
                set(OPT_FLAGS1 "-O0")
            endif()
        elseif(APPLE)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-Wall")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-Wall")
                set(OPT_FLAGS "-O")
                set(OPT_FLAGS1 "-O0")
            endif()    
        else()
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-Wall")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-Wall")
                set(OPT_FLAGS "-O")
                set(OPT_FLAGS1 "-O0")
            endif()    
        endif()

    elseif(COMPILER_NAME STREQUAL gfortran)

        if(WIN32)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-ftrace=full")
                set(OPT_FLAGS "-O0")
                set(OPT_FLAGS1 "-O0")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "")
                set(OPT_FLAGS "-O3 -funroll-loops -msse2")
                set(OPT_FLAGS1 "-O0")
            endif()
        elseif(APPLE)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-Wall -m32")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-m32")
                set(OPT_FLAGS "-O0")
                set(OPT_FLAGS1 "-O0")
            endif()    
        else()
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-Wall")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-Wall")
                set(OPT_FLAGS "-O0")
                set(OPT_FLAGS1 "-O0")
            endif()    
        endif()
    
    endif()
        
endmacro()