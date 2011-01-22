@echo off
rem ****
rem ****---- Compilation for CRYSCALC Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_CrysCalc: Make CrysCalc Compilation
   echo    Syntax: make_CrysCalc [f95/lf95/g95/gfrotran/ifort]
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
   lf95 -c menu_1.f90  -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95 -c menu_2.f90  -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95 -c menu_3.f90  -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95 -c menu_4.f90  -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95 -c calsym.f90  -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95  *.obj -out CrysCalc -tp -nomap -stchk -nchk -o1 -lib c:\CrysFML\lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c menu_1.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
   ifort /c menu_2.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
   ifort /c menu_3.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
   ifort /c menu_4.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
   ifort /c calsym.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
   rem ifort /exe:CrysCalc *.obj C:\CrysFML\Intel\LibC\crysfml.lib
   link /subsystem:console /out:CrysCalc.exe *.obj C:\CrysFML\Intel\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3 -funroll-loops  -msse2   menu_1.f90     -IC:\CrysFML\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   menu_2.f90     -IC:\CrysFML\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   menu_3.f90     -IC:\CrysFML\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   menu_4.f90     -IC:\CrysFML\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   calsym.f90     -IC:\CrysFML\G95\LibC
   g95  *.o -o cryscalc -O3  -funroll-loops  -msse2  -LC:\CrysFML\G95\LibC -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3 -funroll-loops  -msse2   menu_1.f90     -IC:\CrysFML\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   menu_2.f90     -IC:\CrysFML\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   menu_3.f90     -IC:\CrysFML\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   menu_4.f90     -IC:\CrysFML\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   calsym.f90     -IC:\CrysFML\GFortran\LibC
   gfortran *.o -o cryscalc -O3  -funroll-loops  -msse2  -LC:\CrysFML\GFortran\LibC -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
