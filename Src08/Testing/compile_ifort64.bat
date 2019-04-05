@echo off
rem
rem Intel Compilation
rem
   ifort %1.f90         /c /O3 /nologo /I"%CRYSFML%"\ifort64\LibC08
   ifort /exe:Testing *.obj "%CRYSFML%"\ifort64\LibC08\crysfml.lib
   
   del *.obj *.mod > nul

