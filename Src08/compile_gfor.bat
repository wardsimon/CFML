@echo off
rem
rem GFortran Compilation
rem
   gfortran -c CFML_GlobalDeps_Windows_GFOR.f90 -ffree-line-length-0
   ar cr libcrysfml.a *.o   
   del *.o > nul
   move *.mod %CRYSFML%/gfortran/LibC08
   move *.a %CRYSFML%/gfortran/LibC08
   
   

