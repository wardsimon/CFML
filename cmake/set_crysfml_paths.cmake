macro(set_crysfml_paths CRYSFML_PATHS)
    
    if(WIN32)
        set(PATHS $ENV{USERPROFILE}/libcrysfml/crysfml
                  $ENV{HOMEDRIVE}/libcrysfml/crysfml
                  $ENV{SYSTEMDRIVE}/libcrysfml/crysfml
                  $ENV{LIB})
    else()
        set(PATHS /usr/local/lib/libcrysfml/crysfml
                  /usr/lib/libcrysfml/crysfml
                  /opt/lib/libcrysfml/crysfml
                  /opt/libcrysfml/crysfml
                  $ENV{HOME}/lib/libcrysfml/crysfml
                  $ENV{HOME}/libcrysfml/crysfml)
    endif()

            
    set(PATHS ${PATHS} $ENV{CRYSFML}/crysfml)
    
    set(${CRYSFML_PATHS} ${PATHS})
    
    unset(PATHS)
        
endmacro()