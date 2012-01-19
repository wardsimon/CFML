include(LibFindMacros)

include(set_crysfgl_paths)
set_crysfgl_paths(CRYSFGL_PATHS)

find_path(CRYSFGL_MOD_DIR NAMES wcfgl_constant.mod PATHS ${CRYSFGL_PATHS})

find_library(CRYSFGL_LIBRARY NAMES crysfgl PATHS ${CRYSFGL_PATHS})

set(CRYSFGL_PROCESS_MODS CRYSFGL_MOD_DIR)

set(CRYSFGL_PROCESS_LIBS CRYSFGL_LIBRARY)

libfind_process(CRYSFGL)