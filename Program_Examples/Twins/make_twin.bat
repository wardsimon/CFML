@echo off
rem ****
rem ****---- Compilation for TWIN_LAWS Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: Jan-2011
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_Twin: Make TWIN_LAWS Compilation
   echo    Syntax: make_Twin [f95/lf95/g95/gfortran/ifort]
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
   lf95 -c Search_TwinLaws.f90  -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95  *.obj -out Search_TwinLaws -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC" -lib c:\CrysFML\lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Search_TwinLaws.f90   /O2 /nologo /IC:\CrysFML\Intel\LibC
   ifort /exe:Search_TwinLaws *.obj C:\CrysFML\Intel\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c Search_TwinLaws.f90    -I../../G95/LibC
   g95 *.o -o Search_TwinLaws    -L../../G95/LibC   -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Search_TwinLaws.f90    -I../../GFortran/LibC
   gfortran *.o -o Search_TwinLaws    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
