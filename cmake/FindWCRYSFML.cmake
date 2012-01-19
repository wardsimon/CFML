include(LibFindMacros)

include(set_wcrysfml_paths)
set_wcrysfml_paths(WCRYSFML_PATHS)

find_path(WCRYSFML_MOD_DIR NAMES cfml_globaldeps.mod PATHS ${WCRYSFML_PATHS})

find_library(WCRYSFML_LIBRARY NAMES wcrysfml PATHS ${WCRYSFML_PATHS})

set(WCRYSFML_PROCESS_MODS WCRYSFML_MOD_DIR)

set(WCRYSFML_PROCESS_LIBS WCRYSFML_LIBRARY)

libfind_process(WCRYSFML)