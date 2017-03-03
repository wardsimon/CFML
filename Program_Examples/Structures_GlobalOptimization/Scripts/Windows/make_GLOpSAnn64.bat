@echo off
rem ****
rem ****---- Compilation of GlOpSAnn Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_GLOpSAnn64: Make GLOpSAnn64 Compilation
   echo    Syntax: make_GLOpSAnn [f95/lf95/g95/gfortran/ifort]
   goto END
rem
:CONT
cd ..\..\Src
   if x%1 == xf95       goto F95
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem
rem ****---- Absoft Compiler ----****
:F95
   echo "Absoft compiler not available"
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   echo "Lahey compiler not available for 64 bits"
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c observ.f90            /O2 /nologo /heap-arrays /IC:\CrysFML\ifort64\LibC
   ifort /c cost_functions.f90    /O2 /nologo /heap-arrays /IC:\CrysFML\ifort64\LibC
   ifort /c GLOpSAnn.f90          /O2 /nologo /heap-arrays /IC:\CrysFML\ifort64\LibC
   ifort /exe:GLOpSAnn *.obj  C:\CrysFML\ifort64\LibC\CrysFML.lib /link /stack:64000000
rem   ifort  /exe:GLOpSAnn.exe *.obj  C:\CrysFML\ifort64\LibC\CrysFML.lib
   goto END
:IFORTD
   ifort /c observ.f90           /heap-arrays   /debug=full /traceback /nologo /IC:\CrysFML\ifort64_debug\LibC
   ifort /c cost_functions.f90   /heap-arrays   /debug=full /traceback /nologo /IC:\CrysFML\ifort64_debug\LibC
   ifort /c GLOpSAnn.f90         /heap-arrays   /debug=full /traceback /nologo /IC:\CrysFML\ifort64_debug\LibC
   ifort   /exe:GLOpSAnn.exe *.obj  C:\CrysFML\ifort64_debug\LibC\CrysFML.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   echo "G95 compiler not available for 64 bits"
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c observ.f90          -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran64\LibC
   gfortran -c cost_functions.f90  -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran64\LibC
   gfortran -c GLOpSAnn.f90        -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran64\LibC
   gfortran  *.o -o  GLOpSAnn_gf  -LC:\CrysFML\GFortran64\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
:END
   upx GLOpSAnn.exe
if exist %FULLPROF%  copy GLOpSAnn.exe %FULLPROF%\.
if exist %PROGCFML%  move GLOpSAnn.exe %PROGCFML%\DistFPS_64b\.
if exist %PROGCFML%  copy ..\Docs\GLOpSAnn.pdf %PROGCFML%\DistFPS\Docs\.
   del *.obj *.mod *.o *.map *.bak *.exe > nul
cd ..\Scripts\Windows
