@echo off
rem
rem Intel Compilation
cd %CRYSFML%\Program_Examples\BondStr\Src
rem
   ifort Bond_Str.f90         /c /O3 /nologo /I"%CRYSFML%"\ifort64\libc
   ifort /exe:Bond_Str *.obj "%CRYSFML%"\ifort64\libc\crysfml.lib
rem
rem Compress executable
rem
   upx bond_str.exe
rem
rem Update FullProf Distribution
rem
   if exist %PROGCFML% copy Bond_Str.exe %PROGCFML%\DistFPS
   if exist %PROGCFML% copy Bond_Str.f90 %PROGCFML%\BondStr\Src
   if exist %FULLPROF% copy Bond_Str.exe %FULLPROF%
rem
rem Clean several files
rem
   del *.obj *.o *.mod *.exe *.map *.bak
cd %CRYSFML%\Program_Examples\BondStr\Scripts\Windows
