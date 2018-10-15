@echo off
rem ****
rem ****---- Compilation for MagGroupk Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_MAGGROUP: Make MagGroupk Compilation
   echo    Syntax: make_MagGroup [f95/lf95/g95/gfortran/ifort] [deb]
   goto END
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   if x%2 == xdeb goto IFORTD
   ifort /c MagGroups_K.f90 /O2 /nologo /IC:\CrysFML\ifort\LibC
   rem ifort /exe:MagGroups_K *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   link /subsystem:console /out:MagGroups_K.exe *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   goto END
:IFORTD
   ifort /c MagGroups_K.f90 /debug:full /check /traceback  /nologo  /heap-arrays:100 /IC:\CrysFML\ifort_debug\LibC
   rem ifort /exe:MagGroups_K *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   link /subsystem:console /out:MagGroups_K.exe *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c MagGroups_K.f90    -I../../GFortran/LibC
   gfortran *.o -o MagGroups_K    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
