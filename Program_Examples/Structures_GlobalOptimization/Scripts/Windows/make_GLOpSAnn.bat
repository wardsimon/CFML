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
   echo    MAKE_GLOpSAnn: Make GLOpSAnn Compilation
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
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c observ.f90             -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95 -c cost_functions.f90     -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95 -c GLOpSAnn.f90           -tp  -nstchk -nchk  -o3 -mod ".;c:\CrysFML\Lahey\LibC"
   lf95  *.obj -out GLOpSAnn_lf -tp  -nstchk -nchk -o3 -lib c:\CrysFML\Lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c observ.f90             /O2 /nologo /IC:\CrysFML\ifort\LibC
   ifort /c cost_functions.f90     /O2 /nologo /IC:\CrysFML\ifort\LibC
   ifort /c GLOpSAnn.f90           /O2 /nologo /IC:\CrysFML\ifort\LibC
   rem ifort /exe:GLOpSAnn *.obj  C:\CrysFML\ifort\LibC\CrysFML.lib /link /stack:64000000
   link /subsystem:console /stack:64000000 /out:GLOpSAnn.exe *.obj  C:\CrysFML\ifort\LibC\CrysFML.lib
   goto END
:IFORTD
   ifort /c observ.f90           /heap-arrays   /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c cost_functions.f90   /heap-arrays   /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c GLOpSAnn.f90         /heap-arrays   /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /subsystem:console  /out:GLOpSAnn.exe *.obj  C:\CrysFML\ifort_debug\LibC\CrysFML.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c observ.f90          -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95 -c cost_functions.f90  -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95 -c GLOpSAnn.f90        -O3  -std=f2003  -funroll-loops  -msse2  -IC:\CrysFML\G95\LibC
   g95  *.o -o  GLOpSAnn_g95  -LC:\CrysFML\G95\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c observ.f90          -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran -c cost_functions.f90  -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran -c GLOpSAnn.f90        -O3 -funroll-loops  -msse2  -IC:\CrysFML\GFortran\LibC
   gfortran  *.o -o  GLOpSAnn_gf  -LC:\CrysFML\GFortran\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem
:END
   upx GLOpSAnn.exe
if exist %FULLPROF%  copy GLOpSAnn.exe %FULLPROF%\. > nul
if exist %PROGCFML%  move GLOpSAnn.exe %PROGCFML%\DistFPS\. > nul
if exist %PROGCFML%  copy ..\Docs\GLOpSAnn.pdf %PROGCFML%\DistFPS\Docs\. > nul
   del *.obj *.mod *.o *.map *.bak *.exe > nul
   cd ..\Scripts\Windows
