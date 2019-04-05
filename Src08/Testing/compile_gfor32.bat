@echo off
rem
rem GFortran Compilation
rem
   gfortran -c %1.f90   -I"%CRYSFML%"\gfortran\LibC08
   gfortran *.o -o Testing -L"%CRYSFML%"\gfortran\LibC08 -lcrysfml
   del *.o > nul

