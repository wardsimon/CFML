@echo off
rem
rem Intel Compilation
rem
   ifort %1.f90         /c /O3 /nologo /I"%CRYSFML%"\ifort\libc08
   ifort /exe:Testing *.obj "%CRYSFML%"\ifort\libc08\crysfml.lib

