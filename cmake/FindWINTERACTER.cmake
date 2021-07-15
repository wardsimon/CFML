include(LibFindMacros)

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

set(WINTERACTER "$ENV{WINTERACTER}")
set(WINTER "$ENV{WINTER}")
set(WINT "$ENV{WINT}")


if (WIN32 OR MSYS)
   # Windows 

    string(REGEX REPLACE "\\\\" "/" WINTERACTER "${WINTERACTER}")
    string(REGEX REPLACE "\\\\" "/" WINTER "${WINTER}")
    string(REGEX REPLACE "\\\\" "/" WINT "${WINT}")

    set(USERPROFILE "$ENV{USERPROFILE}")
    set(HOMEDRIVE "$ENV{HOMEDRIVE}")
    set(SYSTEMDRIVE "$ENV{SYSTEMDRIVE}")

    if(COMPILER_NAME STREQUAL ifort)

        if(${ARCH32})
            find_path(WINTERACTER_MOD_DIR
                      NAMES winteracter.mod
                      PATHS ${WINTERACTER}/lib.if8
                            ${WINTER}/lib.if8
                            ${WINT}/lib.if8
                            ${USERPROFILE}/wint/lib.if8
                            ${HOMEDRIVE}/wint/lib.if8
                            ${SYSTEMDRIVE}/wint/lib.if8)
        else(${ARCH32})
            find_path(WINTERACTER_MOD_DIR
                      NAMES winteracter.mod
                      PATHS ${WINTERACTER}/lib.i64
                            ${WINTER}/lib.i64
                            ${WINT}/lib.i64
                            ${USERPROFILE}/wint/lib.i64
                            ${HOMEDRIVE}/wint/lib.i64
                            ${SYSTEMDRIVE}/wint/lib.i64)
        endif(${ARCH32})

    elseif(COMPILER_NAME STREQUAL gfortran)
    
         if(${ARCH32})
            find_path(WINTERACTER_MOD_DIR
                      NAMES winteracter.mod
                      PATHS ${WINTERACTER}/lib.gnu32
                            ${WINTER}/lib.gnu32
                            ${WINT}/lib.gnu32
                            ${USERPROFILE}/wint/lib.gnu32
                            ${HOMEDRIVE}/wint/lib.gnu32
                            ${SYSTEMDRIVE}/wint/lib.gnu32)
        else(${ARCH32})
            find_path(WINTERACTER_MOD_DIR
                      NAMES winteracter.mod
                      PATHS ${WINTERACTER}/lib.gnu64
                            ${WINTER}/lib.gnu64
                            ${WINT}/lib.gnu64
                            ${USERPROFILE}/wint/lib.gnu64
                            ${HOMEDRIVE}/wint/lib.gnu64
                            ${SYSTEMDRIVE}/wint/lib.gnu64)
        endif(${ARCH32})

    endif()

    libfind_library(USER32 user32)
    libfind_library(GDI32 gdi32)
    libfind_library(COMDLG32 comdlg32)
    libfind_library(WINSPOOL winspool)
    libfind_library(WINMM winmm)
    libfind_library(SHELL32 shell32)
    libfind_library(ADVAPI32 advapi32)
    libfind_library(HTMLHELP htmlhelp)

    find_library(WINTERACTER_LIBRARY
                 NAMES winter wint
                 PATHS ${WINTERACTER_MOD_DIR})

    set(WINTERACTER_PROCESS_LIBS WINTERACTER_LIBRARY
                                 USER32_LIBRARY
                                 GDI32_LIBRARY
                                 COMDLG32_LIBRARY
                                 WINSPOOL_LIBRARY
                                 WINMM_LIBRARY
                                 SHELL32_LIBRARY
                                 ADVAPI32_LIBRARY
                                 HTMLHELP_LIBRARY)

elseif(APPLE)
   # MacOS
   if(COMPILER_NAME STREQUAL ifort)

        
	    find_path(WINTERACTER_MOD_DIR
                  NAMES winteracter.mod
                  PATHS ${WINTERACTER}/lib.ifi64                             
                        ${WINTER}/lib.ifi64
                        ${WINT}/lib.ifi64
                        ${HOME}/wint/lib.ifi64
                        ${HOME}/lib.ifi64
                        /usr/local/lib/wint/lib.ifi64
                        /usr/lib/wint/lib.ifi64
                        /opt/lib/wint/lib.ifi64)
        
		
    elseif(COMPILER_NAME STREQUAL gfortran)

        find_path(WINTERACTER_MOD_DIR
                  NAMES winteracter.mod
                  PATHS ${WINTERACTER}/lib.gni64/8.2
                        ${WINTER}/lib.gni64/8.2
                        ${WINT}/lib.gni64/8.2
                        ${HOME}/wint/lib.gni64/8.2
                        ${HOME}/lib.gni64/8.2
                        /usr/local/lib/wint/lib.gni64/8.2
                        /usr/lib/wint/lib.gni64/8.2
                        /opt/lib/wint/lib.gni64/8.2)

    endif()

    #libfind_library(XM Xm)
    #libfind_library(XT Xt)
    #libfind_library(XMU Xmu)
    #libfind_library(X11 X11)
    #libfind_library(XEXT Xext)
    #libfind_library(SM SM)
    #libfind_library(ICE ICE)
    #libfind_library(XFT Xft)
    #libfind_library(FONTCONFIG fontconfig)
    #libfind_library(XINERAMA Xinerama)
    #libfind_library(ICONV iconv)
    
    find_library(XM_LIBRARY NAMES Xm PATHS ENV LD_LIBRARY_PATH)
    find_library(XT_LIBRARY NAMES Xt PATHS ENV LD_LIBRARY_PATH)
    find_library(XMU_LIBRARY NAMES Xmu PATHS ENV LD_LIBRARY_PATH)
    find_library(X11_LIBRARY NAMES X11 PATHS ENV LD_LIBRARY_PATH)
    find_library(XEXT_LIBRARY NAMES Xext PATHS ENV LD_LIBRARY_PATH)
    #find_library(SM_LIBRARY NAMES SM PATHS ENV LD_LIBRARY_PATH)
    #find_library(ICE_LIBRARY NAMES ICE PATHS ENV LD_LIBRARY_PATH)
    find_library(XFT_LIBRARY NAMES Xft PATHS ENV LD_LIBRARY_PATH)
    find_library(FONTCONFIG_LIBRARY NAMES fontconfig PATHS ENV LD_LIBRARY_PATH)
    find_library(XINERAMA_LIBRARY NAMES Xinerama PATHS ENV LD_LIBRARY_PATH)
    #find_library(ICONV_LIBRARY NAMES iconv PATHS ENV LD_LIBRARY_PATH)
    

    find_library(WINTERACTER_LIBRARY NAMES winter wint PATHS ${WINTERACTER_MOD_DIR})

    find_library(WINTGL_LIBRARY NAMES winterGL wintGL PATHS ${WINTERACTER_MOD_DIR})

    set(WINTERACTER_PROCESS_MODS WINTERACTER_MOD_DIR)

    set(WINTERACTER_PROCESS_LIBS WINTERACTER_LIBRARY
                                 WINTGL_LIBRARY
                                 XM_LIBRARY
                                 XT_LIBRARY
                                 XMU_LIBRARY
                                 X11_LIBRARY
                                 XEXT_LIBRARY
                                 XFT_LIBRARY
                                 FONTCONFIG_LIBRARY
                                 XINERAMA_LIBRARY)
   
else()
    # Linux

    if(COMPILER_NAME STREQUAL ifort)

	    find_path(WINTERACTER_MOD_DIR
                   NAMES winteracter.mod
                   PATHS ${WINTERACTER}/lib.i64
                         ${WINTER}/lib.i64
                         ${WINT}/lib.i64
                         ${HOME}/wint/lib.i64
                         ${HOME}/lib.i64
                         /usr/local/lib/wint/lib.i64
                         /usr/lib/wint/lib.i64
                         /opt/lib/wint/lib.i64)
		
    elseif(COMPILER_NAME STREQUAL gfortran)

        find_path(WINTERACTER_MOD_DIR
                  NAMES winteracter.mod
                  PATHS ${WINTERACTER}/lib.gnu64/8.3
                        ${WINTER}/lib.gnu64/8.3
                        ${WINT}/lib.gnu64/8.3
                        ${HOME}/wint/lib.gnu64/8.3
                        ${HOME}/lib.gnu64/8.3
                        /usr/local/lib/wint/lib.gnu64/8.3
                        /usr/lib/wint/lib.gnu64/8.3
                        /opt/lib/wint/lib.gnu64/8.3)

    endif()

    #libfind_library(XM Xm)
    #libfind_library(XT Xt)
    #libfind_library(XMU Xmu)
    #libfind_library(X11 X11)
    #libfind_library(XEXT Xext)
    #libfind_library(SM SM)
    #libfind_library(ICE ICE)
    #libfind_library(XFT Xft)
    #libfind_library(FONTCONFIG fontconfig)
    #libfind_library(XINERAMA Xinerama)
    
    find_library(XM_LIBRARY NAMES Xm PATHS ENV LD_LIBRARY_PATH)
    find_library(XT_LIBRARY NAMES Xt PATHS ENV LD_LIBRARY_PATH)
    find_library(XMU_LIBRARY NAMES Xmu PATHS ENV LD_LIBRARY_PATH)
    find_library(X11_LIBRARY NAMES X11 PATHS ENV LD_LIBRARY_PATH)
    find_library(XEXT_LIBRARY NAMES Xext PATHS ENV LD_LIBRARY_PATH)
    #find_library(SM_LIBRARY NAMES SM PATHS ENV LD_LIBRARY_PATH)
    #find_library(ICE_LIBRARY NAMES ICE PATHS ENV LD_LIBRARY_PATH)
    find_library(XFT_LIBRARY NAMES Xft PATHS ENV LD_LIBRARY_PATH)
    find_library(FONTCONFIG_LIBRARY NAMES fontconfig PATHS ENV LD_LIBRARY_PATH)
    find_library(XINERAMA_LIBRARY NAMES Xinerama PATHS ENV LD_LIBRARY_PATH)

    find_library(WINTERACTER_LIBRARY NAMES winter wint PATHS ${WINTERACTER_MOD_DIR})

    find_library(WINTGL_LIBRARY NAMES winterGL wintGL PATHS ${WINTERACTER_MOD_DIR})

    set(WINTERACTER_PROCESS_MODS WINTERACTER_MOD_DIR)

    set(WINTERACTER_PROCESS_LIBS WINTERACTER_LIBRARY
                                 WINTGL_LIBRARY
                                 XM_LIBRARY
                                 XT_LIBRARY
                                 XMU_LIBRARY
                                 X11_LIBRARY
                                 XEXT_LIBRARY
                                 XFT_LIBRARY
                                 FONTCONFIG_LIBRARY
                                 XINERAMA_LIBRARY)

endif()

get_filename_component(WINTERACTER_INCLUDE_DIR ${WINTERACTER_MOD_DIR} PATH)

set(WINTERACTER_INCLUDE_DIR ${WINTERACTER_INCLUDE_DIR}/include)

set(WINTERACTER_PROCESS_INCLUDES WINTERACTER_INCLUDE_DIR)

set(WINTERACTER_PROCESS_MODS WINTERACTER_MOD_DIR)

libfind_process(WINTERACTER)
