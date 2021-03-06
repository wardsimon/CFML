@echo off
rem ****
rem ****---- Compilation for Simple_HKL_GEN Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_sHKL_GEN: Make simple_hkl_gen Compilation
   echo    Syntax: make_shkl_gen [f95/lf95/g95/gfortran/ifort]
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
   lf95 -c simple_hkl_gen.f90    -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 *.obj -out simple_hkl_gen  -o1 -lib C:\crysFML\lahey\libC\crysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c simple_hkl_gen.f90 /O2 /nologo /I%CRYSFML%\ifort\LibC
   rem ifort /exe:simple_hkl_gen *.obj %CRYSFML%\ifort\LibC\crysfml.lib
   link /subsystem:console /out:simple_hkl_gen.exe *.obj %CRYSFML%\ifort\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3 -funroll-loops  -msse2   simple_hkl_gen.f90     -IC:\CrysFML\G95\LibC
   g95  *.o -o simple_hkl_gen -O3  -funroll-loops  -msse2  -LC:\CrysFML\G95\LibC -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3 -funroll-loops  -msse2   simple_hkl_gen.f90     -IC:\CrysFML\GFortran\LibC
   gfortran  *.o -o simple_hkl_gen -O3  -funroll-loops  -msse2  -LC:\CrysFML\GFortran\LibC -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
