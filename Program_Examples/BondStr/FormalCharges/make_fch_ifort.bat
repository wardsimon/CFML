@echo off
rem
rem Intel Compilation
rem
   ifort Formal_Charges.f90 /warn  /c /O3 /nologo /I"%CRYSFML%"\ifort\libc
   ifort /exe:Formal_Charges *.obj "%CRYSFML%"\ifort\libc\crysfml.lib
rem   ifort Formal_Charges.f90 /warn  /c /debug=full /traceback /nologo /I"%CRYSFML%"\ifort_debug\libc
rem   ifort /exe:Formal_Charges *.obj "%CRYSFML%"\ifort_debug\libc\crysfml.lib
rem
rem Compress executable
rem
   upx Formal_Charges.exe
rem
   if exist %FULLPROF% copy Formal_Charges.exe %FULLPROF%
rem
rem Clean several files
rem
   del *.obj *.o *.mod *.map *.bak
