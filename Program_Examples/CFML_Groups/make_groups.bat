@echo off
rem ****
rem ****---- Compilation for SIMILAR Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_GROUPS: Testing CFML_Groups
   echo    Syntax: make_groups [gfortran/ifort] [deb]
   goto END
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xgfortrand  goto GFORD
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem ****---- Intel Compiler ----****
:IFORT
   if x%2 == xdeb goto IFORTD
   ifort /c CFML_Rational_Groups.f90  /O2 /Qparallel /nologo /IC:\CrysFML\ifort\LibC
   ifort /c groups.f90                /O2 /Qparallel /nologo /IC:\CrysFML\ifort\LibC
   ifort /exe:groups *.obj C:\CrysFML\ifort\LibC\crysfml.lib /link /stack:64000000
rem   link /subsystem:console /out:get_standard.exe *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   goto END
:IFORTD
   ifort /c CFML_Rational_Groups.f90  /debug:full /check:noarg_temp_created /traceback  /nologo  /heap-arrays:100 /IC:\CrysFML\ifort_debug\LibC
   ifort /c groups.f90       /debug:full /check:noarg_temp_created /traceback  /nologo  /heap-arrays:100 /IC:\CrysFML\ifort_debug\LibC
   ifort /exe:groups *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
rem   link /subsystem:console /out:get_standard.exe *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c CFML_Rational_Groups.f90  -ffree-line-length-0  -I../../GFortran/LibC
   gfortran -c groups.f90       -ffree-line-length-0  -I../../GFortran/LibC
   gfortran *.o -o groups    -L../../GFortran/LibC   -lcrysfml
   goto END
:GFORD
   gfortran -c CFML_Rational_Groups.f90 -g -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC
   gfortran -c groups.f90    -g -fbacktrace  -ffree-line-length-0  -I../../GFortran/LibC
   gfortran *.o -o groups    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak *.pdb > nul
