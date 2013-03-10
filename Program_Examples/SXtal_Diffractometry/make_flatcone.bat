@echo off
rem ****
rem ****---- Compilation for D10_FLAT_CONE program--****
rem ****
rem **** Author: JRC
rem **** Revision: March-2013
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_flatcone: Makes D10_FLAT_CONE Compilation
   echo    Syntax: make_convcell [f95/lf95/g95/gfortran/ifort]
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
   lf95 -c d10_flat_cone.f90  -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95  *.obj -out d10_flat_cone -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC" -lib c:\CrysFML\lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c d10_flat_cone.f90   /O2 /nologo /IC:\CrysFML\ifort\LibC
   rem ifort /exe:d10_flat_cone *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   link /subsystem:console /out:d10_flat_cone.exe *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c d10_flat_cone.f90    -I../../G95/LibC
   g95 *.o -o d10_flat_cone    -L../../G95/LibC   -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c d10_flat_cone.f90    -I../../GFortran/LibC
   gfortran *.o -o d10_flat_cone    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
