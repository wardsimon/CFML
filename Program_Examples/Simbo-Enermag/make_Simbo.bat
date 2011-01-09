@echo off
rem ****
rem ****---- Compilation for SIMBO Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: Jan-2011
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_Simbo: Make SIMBO Compilation
   echo    Syntax: make_simbo [f95/lf95/g95/gfortran/ifort]
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
   lf95 -c Sup_Exc.f90    -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95 -c Simbo.f90      -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95  *.obj -out Simbo -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC" -lib c:\CrysFML\lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Sup_Exc.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
   ifort /c Simbo.f90   /O2 /nologo /IC:\CrysFML\Intel\LibC
   ifort /exe:Simbo *.obj C:\CrysFML\Intel\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c Sup_Exc.f90  -I../../G95/LibC
   g95 -c Simbo.f90    -I../../G95/LibC
   g95 *.o -o Simbo    -L../../G95/LibC   -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Sup_Exc.f90  -I../../GFortran/LibC
   gfortran -c Simbo.f90    -I../../GFortran/LibC
   gfortran *.o -o Simbo    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
