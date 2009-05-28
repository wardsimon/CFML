@echo off
rem ****
rem ****---- Compilation for OPTIM_GEN Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_OPTIM_GEN: Make OPTIM_GEN Compilation
   echo    Syntax: make_optim_gen [f95/lf95/g95/gfortran/ifort]
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
   lf95 -c observ.f90             -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95 -c cost_functions.f90     -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95 -c Optim_General.f90      -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95  *.obj -out Optim_General_lf -tp  -nstchk -nchk -o3 -lib c:\CrysFML\Lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c observ.f90             /O3 /nologo /IC:\CrysFML\Intel\LibC
   ifort /c cost_functions.f90     /O3 /nologo /IC:\CrysFML\Intel\LibC
   ifort /c Optim_General.f90      /O3 /nologo /IC:\CrysFML\Intel\LibC
   ifort /exe:Optim_General_if *.obj  C:\CrysFML\Intel\LibC\CrysFML.lib /link /stack:64000000
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c observ.f90          -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95 -c cost_functions.f90  -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95 -c Optim_General.f90   -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95  *.o -o  Optim_General_g95  -LC:\CrysFML\G95\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c observ.f90          -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran -c cost_functions.f90  -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran -c Optim_General.f90   -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran  *.o -o  Optim_General_gf  -LC:\CrysFML\GFortran\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
