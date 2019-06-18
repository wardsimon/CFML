if(WIN32)
set(NEXUS_HINTS_DIR "$ENV{CRYSFML_NEXUS_DIR}" "$ENV{CRYSFML_NEXUS_DIR}\\include")
else()
set(NEXUS_HINTS_DIR "$ENV{CRYSFML_NEXUS_DIR}" "$ENV{CRYSFML_NEXUS_DIR}/include" /usr/local/include/nexus /usr/include/nexus /usr/local/include)
endif()
find_path(NEXUS_INCLUDE_DIR napi.h HINTS ${NEXUS_HINTS_DIR})

if (NEXUS_INCLUDE_DIR)
    set(NEXUS_FOUND TRUE)
else ()
    message(STATUS "No NEXUS found.")
endif ()