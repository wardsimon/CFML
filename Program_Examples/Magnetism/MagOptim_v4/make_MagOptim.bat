@echo off
rem ****
rem ****---- Compilation for MagOptim Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: June-2012
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_MagOptim: Make MagOptim Compilation
   echo    Syntax: make_MagOptim [f95/lf95/g95/gfortran/ifort]
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
   lf95 -c Prep_Input.f90          -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95 -c Cost_MagFunctions.f90   -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95 -c MagOptim.f90            -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95  *.obj -out MagOptim_lf    -tp  -nstchk -nchk -o3 -lib c:\CrysFML\Lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Prep_Input.f90              /O2 /nologo /IC:\CrysFML\ifort\LibC
   ifort /c Cost_MagFunctions.f90       /O2 /nologo /IC:\CrysFML\ifort\LibC
   ifort /c MagOptim.f90                /O2 /nologo /IC:\CrysFML\ifort\LibC
   link /subsystem:console /stack:64000000 /out:MagOptim.exe *.obj  C:\CrysFML\ifort\LibC\CrysFML.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c Prep_Input.f90           -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95 -c Cost_MagFunctions.f90    -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95 -c MagOptim.f90             -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95  *.o -o  MagOptim_g95  -LC:\CrysFML\G95\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Prep_Input.f90           -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran -c Cost_MagFunctions.f90    -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran -c MagOptim.f90             -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran  *.o -o  MagOptim_gf  -LC:\CrysFML\GFortran\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
