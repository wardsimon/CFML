@echo off
rem
rem Intel Compilation
rem
   ifort CFML_Atomic_Data.f90 /c /O1 /nologo /I"%CRYSFML%"\ifort\libc
   ifort Formal_Charges.f90   /c /O3 /nologo /I"%CRYSFML%"\ifort\libc
   ifort /exe:Formal_Charges *.obj "%CRYSFML%"\ifort\libc\crysfml.lib
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
