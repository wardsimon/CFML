# Try to find by casual method if include path or library path (or both) is not defined
if(NOT HDF5_INCLUDE_PATH AND NOT HDF5_LIBRARY_PATH)
    find_package(HDF5)
endif()

# Override include path if user defined it
if(HDF5_INCLUDE_PATH)
    set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_PATH})
    message(STATUS "Using ${HDF5_INCLUDE_PATH} as HDF5 headers (provided by user)")
endif()

# Override library path if user defined it
if(HDF5_LIBRARY_PATH)
	find_library(hdf5_fortran NAME libhdf5_fortran PATHS ${HDF5_LIBRARY_PATH})
	if(hdf5_fortran STREQUAL hdf5_fortran-NOTFOUND)
	    message(FATAL_ERROR "libhdf5_fortran not found")
	endif()
	
	find_library(hdf5_f90cstub NAME libhdf5_f90cstub PATHS ${HDF5_LIBRARY_PATH})
	if(hdf5_f90cstub STREQUAL hdf5_f90cstub-NOTFOUND)
	    message(FATAL_ERROR "libhdf5_f90cstub not found")
	endif()
	
	find_library(hdf5 NAME libhdf5 PATHS ${HDF5_LIBRARY_PATH})
	if(hdf5 STREQUAL hdf5-NOTFOUND)
	    message(FATAL_ERROR "libhdf5 not found")
	endif()
	
	find_library(szip NAME libszip PATHS ${HDF5_LIBRARY_PATH})
	if(szip STREQUAL szip-NOTFOUND)
	    message(FATAL_ERROR "libszip not found")
	endif()
	
	find_library(zlib NAME libzlib PATHS ${HDF5_LIBRARY_PATH})
	if(zlib STREQUAL zlib-NOTFOUND)
	    message(FATAL_ERROR "libzlib not found")
	endif()
 
endif()

# Set HDF5_libs
set(HDF5_LIBS libhdf5_fortran libhdf5_f90cstub libhdf5 libszip libzlib)