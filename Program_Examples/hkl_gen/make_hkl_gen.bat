@echo off
rem ****
rem ****---- Compilation for HKL_GEN Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_HKL_GEN: Make HKL_GEN Compilation
   echo    Syntax: make_hkl_gen [f95/lf95/g95/gfrotran/ifort]
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
   lf95 -c hkl_gen.f90    -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 *.obj -out hkl_gen  -o1 -mod ".;C:\crysFML\lahey\LibC" -lib C:\crysFML\lahey\libC\crysFML 
   goto END
rem      
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c hkl_gen.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
   ifort /exe:hkl_gen *.obj C:\CrysFML\Intel\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   goto END
rem
   del *.obj *.o *.mod *.exe *.map
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
