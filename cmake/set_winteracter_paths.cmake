macro(set_winteracter_paths WINTERACTER_PATHS)
    
    if(WIN32)
        set(PATHS $ENV{USERPROFILE}/wint
                  $ENV{HOMEDRIVE}/wint
                  $ENV{SYSTEMDRIVE}/wint
                  $ENV{LIB})
    else()
        set(PATHS /usr/local/lib/wint
                  /usr/lib/wint
                  /opt/lib/wint
                  /opt/wint
                  $ENV{HOME}/lib/wint
                  $ENV{HOME}/wint)    
    endif()

    set(PATHS ${PATHS} $ENV{WINTERACTER} $ENV{WINTER} $ENV{WINT})

    set(${WINTERACTER_PATHS} ${PATHS})
            
endmacro()