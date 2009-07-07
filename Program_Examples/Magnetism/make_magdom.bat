@echo off
rem ****
rem ****---- Compilation for MagDom Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: June-2009
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_CALC_POWDER: Make mag_testdom Compilation
   echo    Syntax: make_magdom [lf95/g95/gfortran/ifort]
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
   lf95 -c mag_testdom.f90    -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 *.obj -out mag_testdom_lf  -o1 -lib C:\crysFML\lahey\libC\crysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c mag_testdom.f90 /Ox /nologo /I. /IC:\CrysFML\Intel\LibC
   ifort /exe:mag_testdom_if *.obj C:\CrysFML\Intel\LibC\CrysFML.lib /link /stack:1024000000
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3  -std=f2003  -funroll-loops  -msse2   mag_testdom.f90   -IC:\CrysFML\G95\LibC
   g95  *.o -o mag_testdom_g95 -O3  -funroll-loops  -msse2  -LC:\CrysFML\G95\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   echo gfortran
   gfortran -c -O3  -std=f2003  -funroll-loops  -msse2   mag_testdom.f90   -IC:\CrysFML\GFortran\LibC
   gfortran  *.o -o mag_testdom_gf -O3  -funroll-loops  -msse2  -LC:\CFML_ILL\GFortran\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
