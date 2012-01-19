macro(set_wcrysfml_paths WCRYSFML_PATHS)
    
    if(WIN32)
        set(PATHS $ENV{USERPROFILE}/libcrysfml/wcrysfml
                  $ENV{HOMEDRIVE}/libcrysfml/wcrysfml
                  $ENV{SYSTEMDRIVE}/libcrysfml/wcrysfml
                  $ENV{LIB})
    else()
        set(PATHS /usr/local/lib/libcrysfml/wcrysfml
                  /usr/lib/libcrysfml/wcrysfml
                  /opt/lib/libcrysfml/wcrysfml
                  /opt/libcrysfml/wcrysfml
                  $ENV{HOME}/lib/libcrysfml/wcrysfml
                  $ENV{HOME}/libcrysfml/wcrysfml)
    endif()
            
    set(PATHS ${PATHS} $ENV{CRYSFML}/wcrysfml)
    
    set(${WCRYSFML_PATHS} ${PATHS})
    
    unset(PATHS)
        
endmacro()