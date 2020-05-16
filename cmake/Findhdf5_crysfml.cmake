# Try to find by casual method if include path or library path (or both) is not defined
if(NOT HDF5_INCLUDE_PATH OR NOT HDF5_LIBRARY_PATH)
    find_package(HDF5)
endif()

# Override include path if user defined it
if(HDF5_INCLUDE_PATH)
    set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_PATH})
    message(STATUS "Using ${HDF5_INCLUDE_PATH} as HDF5 headers (provided by user)")
endif()

# Override library path if user defined it
if(HDF5_LIBRARY_PATH)
    message(STATUS "Try to find HDF5 libraries in ${HDF5_LIBRARY_PATH} (provided by user)")
	if(APPLE)
	    set(LIBS hdf5 m dl sz z)
	elseif(UNIX)
	    set(LIBS hdf5 m dl sz z pthread)
	elseif(WIN32)
	    set(LIBS libhdf5 hdf5_fortran hdf5_f90cstub libszip libzlib)
	else()
	    message(FATAL ERROR "OS unknown")
	endif()
	
	set(HDF5_LIBRARIES)
    foreach(libname ${LIBS})
        set(temp_lib temp_lib-NOTFOUND)
        find_library(temp_lib NAMES ${libname} PATHS ${HDF5_LIBRARY_PATH})# NO_DEFAULT_PATH)
        if(${temp_lib} STREQUAL temp_lib-NOTFOUND)
	        message(FATAL_ERROR "${libname} not found")
	    else()
	        list(APPEND HDF5_LIBRARIES ${temp_lib})
	    endif()
    endforeach()
endif()

message(Found HDF5_LIBRARIES:${HDF5_LIBRARIES})
