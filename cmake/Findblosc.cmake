if(WIN32)
set(BLOSC_HINTS_DIR "$ENV{CRYSFML_BLOSC_DIR}" "$ENV{CRYSFML_BLOSC_DIR}\\include")
else()
set(BLOSC_HINTS_DIR "$ENV{CRYSFML_BLOSC_DIR}" "$ENV{CRYSFML_BLOSC_DIR}/include" /usr/local/include/blosc /usr/include/blosc /usr/local/include)
endif()

find_path(BLOSC_INCLUDE_DIR blosc.h PATHS ${BLOSC_HINTS_DIR})

#find_library(BLOSC_LIBRARY NAMES blosc)

if (BLOSC_INCLUDE_DIR)
    set(BLOSC_FOUND TRUE)
else ()
    message(STATUS "No BLOSC found.")
endif ()