@echo off
rem
rem Intel Compilation
rem
   ifort bond_str.f90 /c /O3 /nologo /I"%CRYSFML%"\ifort\libc
   ifort /exe:bond_str *.obj "%CRYSFML%"\ifort\libc\crysfml.lib
rem
rem Compress executable
rem
   upx bond_str.exe
rem
rem Update FullProf Distribution
rem
   if exist %FULLPROF% copy bond_str.exe %FULLPROF%
rem
rem Clean several files
rem
   del *.obj *.o *.mod *.exe *.map *.bak
