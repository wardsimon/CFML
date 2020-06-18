@echo off
rem ****
rem ****---- Compilation for XRFIT Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: Jan-2011
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_XRFit: Make XRFIT Compilation
   echo    Syntax: make_xrfit [f95/lf95/g95/gfortran/ifort]
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
   lf95 -c Cw_Diffraction_Pv.f90    -tp -nomap -stchk -nchk -o1 -mod ".;%CRYSFML%\lahey\LibC"
   lf95 -c XRFit.f90      -tp -nomap -stchk -nchk -o1 -mod ".;%CRYSFML%\lahey\LibC"
   lf95  *.obj -out XRFit -tp -nomap -stchk -nchk -o1 -mod ".;%CRYSFML%\lahey\LibC" -lib c:\CrysFML\lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Cw_Diffraction_Pv.f90 /O2 /nologo /I%CRYSFML%\ifort64\LibC
   ifort /c XRFit.f90   /O2 /nologo /I%CRYSFML%\ifort64\LibC
   link /subsystem:console /stack:102400000 /out:XRFit.exe *.obj %CRYSFML%\ifort64\LibC\CrysFML.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c Cw_Diffraction_Pv.f90  -I../../G95/LibC
   g95 -c XRFit.f90    -I../../G95/LibC
   g95 *.o -o XRFit    -L../../G95/LibC   -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Cw_Diffraction_Pv.f90  -I../../GFortran/LibC
   gfortran -c XRFit.f90    -I../../GFortran/LibC
   gfortran *.o -o XRFit    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
