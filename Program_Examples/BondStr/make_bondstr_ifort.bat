@echo off
rem
rem Intel Compilation
rem
   ifort BVEL_Percolation\CFML_Percolation.f90 /c /O3 /nologo /I"%CRYSFML%"\ifort\libc
   ifort Bond_Str.f90 /c /O3 /nologo /I"%CRYSFML%"\ifort\libc
   ifort /exe:Bond_Str *.obj "%CRYSFML%"\ifort\libc\crysfml.lib
rem
rem Compress executable
rem
   upx bond_str.exe
rem
rem Update FullProf Distribution
rem
rem   if exist %PROGCFML% copy Bond_Str.exe %PROGCFML%\DistFPS
rem   if exist %PROGCFML% copy Bond_Str.f90 %PROGCFML%\BondStr\Src
   if exist %FULLPROF% copy Bond_Str.exe %FULLPROF%
rem
rem Clean several files
rem
   del *.obj *.o *.mod *.exe *.map *.bak
