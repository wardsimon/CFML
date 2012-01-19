macro(set_crysfgl_paths CRYSFGL_PATHS)
    
    if(WIN32)
        set(PATHS $ENV{USERPROFILE}/crysfgl
                  $ENV{HOMEDRIVE}/crysfgl
                  $ENV{SYSTEMDRIVE}/crysfgl
                  $ENV{LIB})
    else()
        set(PATHS /usr/local/lib/crysfgl
                  /usr/lib/crysfgl
                  /opt/lib/crysfgl
                  /opt/crysfgl
                  $ENV{HOME}/lib/crysfgl
                  $ENV{HOME}/crysfgl)
    endif()
            
    set(PATHS ${PATHS} $ENV{CRYSFGL})
    
    set(${CRYSFGL_PATHS} ${PATHS})
    
    unset(PATHS)
        
endmacro()