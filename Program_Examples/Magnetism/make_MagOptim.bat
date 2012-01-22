@echo off
rem ****
rem ****---- Compilation for Optimizing_MagStructure Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: Janurary-2012
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_MagOptim: Make Optimizing_MagStructure Program Compilation
   echo    Syntax: Make_MagOptim [lf95/g95/gfortran/ifort]
   goto END
rem
:CONT
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c Prep_input.f90          -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 -c Cost_MagFunctions.f90   -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 -c MagOptim.f90            -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 *.obj -out MagStr_Optim_lf  -o1 -lib C:\crysFML\lahey\libC\crysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Prep_input.f90        /Ox /nologo /I. /IC:\CrysFML\ifort\LibC
   ifort /c Cost_MagFunctions.f90 /Ox /nologo /I. /IC:\CrysFML\ifort\LibC
   ifort /c MagOptim.f90          /Ox /nologo /I. /IC:\CrysFML\ifort\LibC
   rem ifort /exe:MagPolar3D_if *.obj C:\CrysFML\ifort\LibC\CrysFML.lib /link /stack:102400000
   link /subsystem:console /stack:102400000 /out:MagStr_Optim.exe *.obj C:\CrysFML\ifort\LibC\CrysFML.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3  -std=f2003  -funroll-loops  -msse2   Prep_input.f90          -IC:\CrysFML\G95\LibC
   g95 -c -O3  -std=f2003  -funroll-loops  -msse2   Cost_MagFunctions.f90   -IC:\CrysFML\G95\LibC
   g95 -c -O3  -std=f2003  -funroll-loops  -msse2   MagOptim.f90            -IC:\CrysFML\G95\LibC
   g95  *.o -o MagStr_Optim_g95 -O3  -funroll-loops  -msse2  -LC:\CrysFML\G95\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3  -std=f2003  -funroll-loops  -msse2   Prep_input.f90          -IC:\CrysFML\GFortran\LibC
   gfortran -c -O3  -std=f2003  -funroll-loops  -msse2   Cost_MagFunctions.f90   -IC:\CrysFML\GFortran\LibC
   gfortran -c -O3  -std=f2003  -funroll-loops  -msse2   MagOptim.f90            -IC:\CrysFML\GFortran\LibC
   gfortran  *.o -o MagStr_Optim_gf -O3  -funroll-loops  -msse2  -LC:\CrysFML\GFortran\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
