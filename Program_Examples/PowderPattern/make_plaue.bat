@echo off
rem ****
rem ****---- Compilation for Laue_Powder Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: February-2014
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_PLaue: Make Laue_Powder  Compilation
   echo    Syntax: make_plaue [lf95/g95/gfortran/ifort]
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
   lf95 -c Laue_Powder.f90    -info  -o1 -chk -mod ".;..\..\lahey\libC"
   lf95 *.obj -out Laue_Powder_lf  -o1 -lib ..\..\lahey\libC\crysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Laue_Powder.f90 /Ox /nologo /I. /I..\..\ifort\LibC
   rem ifort /exe:Simple_Calc_Powder *.obj ..\..\ifort\LibC\CrysFML.lib /link /stack:102400000
   link /subsystem:console /stack:102400000 /out:Laue_Powder.exe *.obj ..\..\ifort\LibC\CrysFML.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3  -std=f2003  -funroll-loops  -msse2   Laue_Powder.f90   -I..\..\G95\LibC
   g95  *.o -o Laue_Powder_g95 -O3  -funroll-loops  -msse2  -L..\..\G95\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3  -std=f2003  -funroll-loops  -msse2   Laue_Powder.f90   -I..\..\GFortran\LibC
   gfortran  *.o -o Laue_Powder_gf -O3  -funroll-loops  -msse2  -L..\..\GFortran\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
