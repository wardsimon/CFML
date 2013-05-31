@echo off
rem ****
rem ****---- Compilation for PHAS_DIAG Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: Jan-2011
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_Phas: Make PHAS_DIAG Compilation
   echo    Syntax: make_phas [f95/lf95/g95/gfortran/ifort]
   goto FIN
rem
:CONT
   if x%1 == xf95       goto F95
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   echo    Unknown compiler!
   goto FIN
rem
rem ****---- Absoft Compiler ----****
:F95
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c Sup_Exc.f90    -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95 -c Phase_Diagram.f90    -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
   lf95  *.obj -out Phase_Diagram -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC" -lib c:\CrysFML\lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Sup_Exc.f90 /O2 /nologo /IC:\CrysFML\ifort\LibC
   ifort /c Phase_Diagram.f90   /O2 /nologo /IC:\CrysFML\ifort\LibC
   rem ifort /exe:Phase_Diagram *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   link /subsystem:console /stack:64000000 /out:Phase_Diagram.exe *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c Sup_Exc.f90  -I../../G95/LibC
   g95 -c Phase_Diagram.f90    -I../../G95/LibC
   g95 *.o -o Phase_Diagram    -L../../G95/LibC   -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Sup_Exc.f90  -I../../GFortran/LibC
   gfortran -c Phase_Diagram.f90    -I../../GFortran/LibC
   gfortran *.o -o Phase_Diagram_gf    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
rem
rem  Comment the following lines if upx or %FULLPROF% are not available
rem  or if you want to conserve the object files
rem  Compression of executable
rem        upx Phase_Diagram.exe
rem  Move the excutable to a directory in the Path
rem        if exist %FULLPROF% move Phase_Diagram.exe %FULLPROF% > nul
        if exist %FULLPROF% move Phase_Diagram_gf.exe %FULLPROF% > nul
rem  Remove unnecessary files
        del *.obj *.mod *.o *.map *.bak > nul
:FIN