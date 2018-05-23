ifort /c m_submodules.f90  /debug:full /check /check:noarg_temp_created /traceback /nologo /I. /I"%CRYSFML%"\ifort\LibC
ifort /exe:m_submodules  *.obj C:\CrysFML\ifort\LibC\crysfml.lib


rem ifort /c m_submodules.f90 /O2 /nologo /IC:\CrysFML\ifort\LibC
rem ifort /exe:m_submodules  *.obj C:\CrysFML\ifort\LibC\crysfml.lib
del *.obj *.bak