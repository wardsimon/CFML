@echo off
rem ****
rem ****---- Compilation for SUBGROUPS Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_SUBGROUPS: Make SUBGROUPS Compilation
   echo    Syntax: make_subgroups [f95/lf95/g95/gfortran/ifort]
   goto END
rem
:CONT
   if x%1 == xf95       goto F95
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   goto END
rem
rem ****---- Absoft Compiler ----****
:F95
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c subgroups.f90    -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 *.obj -out subgroups  -o1 -lib C:\crysFML\lahey\libC\crysFML  -chk
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c subgroups.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
   rem ifort /exe:subgroups *.obj C:\CrysFML\Intel\LibC\crysfml.lib
   link /subsystem:console /out:subgroups.exe *.obj C:\CrysFML\Intel\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c subgroups.f90    -I../../G95/LibC
   g95 *.o -o subgroups    -L../../G95/LibC   -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c subgroups.f90    -I../../GFortran/LibC
   gfortran *.o -o subgroups    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
