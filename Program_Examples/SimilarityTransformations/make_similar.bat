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
   echo    MAKE_SIMILAR: Make Similar Compilation
   echo    Syntax: make_Similar [f95/lf95/g95/gfortran/ifort] [deb]
   goto END
rem
:CONT
   if x%1 == xf95       goto F95
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem
rem ****---- Absoft Compiler ----****
:F95
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c Similar.f90    -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 *.obj -out Similar  -o1 -lib C:\crysFML\lahey\libC\crysFML  -chk
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   if x%2 == xdeb goto IFORTD
   ifort /c CFML_SPG_StdRep.f90 /O2 /nologo /IC:\CrysFML\ifort\LibC
   ifort /c Similar.f90 /O2 /Qparallel /nologo /IC:\CrysFML\ifort\LibC
   rem ifort /exe:Similar *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   link /subsystem:console /out:Similar.exe *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   goto END
:IFORTD
   ifort /c CFML_SPG_StdRep.f90 /debug:full /check:noarg_temp_created /traceback  /nologo  /heap-arrays:100 /IC:\CrysFML\ifort_debug\LibC
   ifort /c Similar.f90 /debug:full /check:noarg_temp_created /traceback  /nologo  /heap-arrays:100 /IC:\CrysFML\ifort_debug\LibC
   ifort /exe:Similar *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
rem   link /subsystem:console /out:Similar.exe *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c Similar.f90    -I../../G95/LibC
   g95 *.o -o Similar    -L../../G95/LibC   -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Similar.f90    -I../../GFortran/LibC
   gfortran *.o -o Similar    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
