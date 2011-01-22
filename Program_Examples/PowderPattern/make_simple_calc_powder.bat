@echo off
rem ****
rem ****---- Compilation for Simple_Calculation_of_Powder_Patterns Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: June-2009
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_SIMPLE_CALC_POWDER: Make Simple_Calc_Powder Compilation
   echo    Syntax: make_Simple_Calc_Powder [lf95/g95/gfortran/ifort]
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
   lf95 -c Simple_Calc_Powder.f90    -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 *.obj -out Simple_Calc_Powder_lf  -o1 -lib C:\crysFML\lahey\libC\crysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Simple_Calc_Powder.f90 /Ox /nologo /I. /IC:\CrysFML\Intel\LibC
   rem ifort /exe:Simple_Calc_Powder *.obj C:\CrysFML\Intel\LibC\CrysFML.lib /link /stack:102400000
   link /subsystem:console /stack:102400000 /out:Simple_Calc_Powder.exe *.obj C:\CrysFML\Intel\LibC\CrysFML.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3  -std=f2003  -funroll-loops  -msse2   Simple_calc_powder.f90   -IC:\CrysFML\G95\LibC
   g95  *.o -o Simple_calc_powder_g95 -O3  -funroll-loops  -msse2  -LC:\CrysFML\G95\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3  -std=f2003  -funroll-loops  -msse2   Simple_calc_powder.f90   -IC:\CrysFML\GFortran\LibC
   gfortran  *.o -o Simple_calc_powder_gf -O3  -funroll-loops  -msse2  -LC:\CrysFML\GFortran\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
