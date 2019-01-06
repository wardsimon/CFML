macro(set_compiler_flags)

    # Nullify all the Fortran flags.
    set(CMAKE_Fortran_FLAGS "")

    get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

    if(COMPILER_NAME STREQUAL ifort)

        if(WIN32)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "/debug:full /traceback /nologo")
                set(OPT_FLAGS  "/Od /check:noarg_temp_created ")
                set(OPT_FLAGS1 "/Od /check:noarg_temp_created ")
                set(OPT_FLAGS2 "/Od /check:noarg_temp_created ")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "/Qopt-report=0 /nologo")
                set(OPT_FLAGS "/O2 /Qparallel")
                set(OPT_FLAGS1 "/Od")
                set(OPT_FLAGS2 "/O1 /Qparallel")
                set(OPT_FLAGS3 "/O3 /Qparallel")
            endif()
        elseif(APPLE)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-warn")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-warn -qopt-report=0")
                set(OPT_FLAGS "-O2")
                set(OPT_FLAGS1 "-O0")
                set(OPT_FLAGS2 "-O1")
            endif()
        else()
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn -traceback")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                if(USE_HDF)
                   set(CMAKE_Fortran_FLAGS_RELEASE "-warn -cpp -DUSE_HDF -qopt-report=0")
                else()
		           set(CMAKE_Fortran_FLAGS_RELEASE "-warn -cpp  -qopt-report=0")
                endif() 
                set(OPT_FLAGS "-O2")
                set(OPT_FLAGS1 "-O0")
                set(OPT_FLAGS2 "-O1")
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
                set(OPT_FLAGS2 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-Wall")
                set(OPT_FLAGS "-O3")
                set(OPT_FLAGS1 "-O0")
                set(OPT_FLAGS2 "-O1")
            endif()
        endif()

    elseif(COMPILER_NAME STREQUAL gfortran)

        if(WIN32)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-fbacktrace -ffree-line-length-none")
                set(OPT_FLAGS "-O0")
                set(OPT_FLAGS1 "-O0")
                set(OPT_FLAGS2 "-O0")
                if(USE_HDF)
                   set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -DUSE_HDF -ffree-line-length-none")
                else()
		           set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -ffree-line-length-none")
                endif() 
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                if(USE_HDF)
                   set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -DUSE_HDF -ffree-line-length-none")
                else()
		           set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -ffree-line-length-none")
                endif() 
                set(OPT_FLAGS "-O3 -funroll-loops -msse2")
                set(OPT_FLAGS1 "-O0")
                set(OPT_FLAGS2 "-O1")
            endif()
        elseif(APPLE)
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-Wall -m32 -ffree-line-length-none")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                set(CMAKE_Fortran_FLAGS_RELEASE "-m32 -ffree-line-length-none")
                set(OPT_FLAGS "-O3")
                set(OPT_FLAGS1 "-O0")
            endif()
        else()
            if(CMAKE_BUILD_TYPE STREQUAL Debug)
                set(CMAKE_Fortran_FLAGS_DEBUG "-Wall -ffree-line-length-none")
                set(OPT_FLAGS "-g")
                set(OPT_FLAGS1 "-g")
                set(OPT_FLAGS2 "-g")
            elseif(CMAKE_BUILD_TYPE STREQUAL Release)
                if(USE_HDF)
                   set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -DUSE_HDF -ffree-line-length-none")
                else()
		           set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -ffree-line-length-none")
                endif() 
                set(OPT_FLAGS "-O2")
                set(OPT_FLAGS1 "-O0")
                set(OPT_FLAGS2 "-O1")
            endif()
        endif()

    endif()

endmacro()
