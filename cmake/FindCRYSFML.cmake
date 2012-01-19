include(LibFindMacros)

include(set_crysfml_paths)
set_crysfml_paths(CRYSFML_PATHS)

find_path(CRYSFML_MOD_DIR NAMES cfml_globaldeps.mod PATHS ${CRYSFML_PATHS})

find_library(CRYSFML_LIBRARY NAMES crysfml PATHS ${CRYSFML_PATHS})

set(CRYSFML_PROCESS_MODS CRYSFML_MOD_DIR)

set(CRYSFML_PROCESS_LIBS CRYSFML_LIBRARY)

libfind_process(CRYSFML)