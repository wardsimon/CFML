

ifort /c ext_submod.f90  /debug:full /check /check:noarg_temp_created /traceback /nologo /I. /I"%CRYSFML%"\ifort_debug\LibC
ifort /exe:ext_submod  *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib


rem ifort /c ext_submod.f90 /O2 /nologo /IC:\CrysFML\ifort\LibC
rem ifort /exe:ext_submod  *.obj C:\CrysFML\ifort\LibC\crysfml.lib
del *.obj *.bak