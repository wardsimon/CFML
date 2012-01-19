macro(ADD_PREFIX prefix rootlist)

  set(outlist )
  foreach(root ${${rootlist}})
    list(APPEND outlist ${prefix}${root})
  endforeach()
  set(${rootlist} ${outlist})
  
endmacro(ADD_PREFIX)

macro(ADD_SUFFIX rootlist suffix)

  set(outlist )
  foreach(root ${${rootlist}})
    list(APPEND outlist ${root}${suffix})
  endforeach()
  set(${rootlist} ${outlist})
  
endmacro(ADD_SUFFIX)

macro(set_winteracter_paths WINTERACTER_PATHS)

    get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

    if(CMAKE_SIZEOF_VOID_P MATCHES 4)
        set(ARCH64 False)
    else()
        set(ARCH64 True)
    endif()
    
    if(WIN32)
        set(PATHS $ENV{USERPROFILE}/wint
                  $ENV{HOMEDRIVE}/wint
                  $ENV{SYSTEMDRIVE}/wint)
    else()
        set(PATHS /usr/local/lib/wint
                  /usr/lib/wint
                  /opt/lib/wint
                  /opt/wint
                  $ENV{HOME}/lib/wint
                  $ENV{HOME}/wint)    
    endif()
    
    set(TEMP)
        
    foreach(f ${PATHS})
        if(COMPILER_NAME STREQUAL ifort)
            if(${ARCH64})
                set(TEMP ${TEMP} ${f}/lib.i64)
            else()
                set(TEMP ${TEMP} ${f}/lib.if8)
            endif()
        elseif(COMPILER_NAME STREQUAL g95)        
            if(${ARCH64})
                set(TEMP ${TEMP} ${f}/lib.g95)
            else()
                set(TEMP ${TEMP} ${f}/lib.g95)
            endif()        
        endif()
        
    endforeach()
        
    set(TEMP ${TEMP} $ENV{WINTERACTER})
    
    set(${WINTERACTER_PATHS} ${TEMP})
    
    unset(TEMP)
        
endmacro()