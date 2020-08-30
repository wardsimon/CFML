@echo off
rem ****
rem ****---- Compilation for SIMBO Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: Jan-2011
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_Simbo: Make SIMBO Compilation
   echo    Syntax: make_simbo [f95/lf95/g95/gfortran/ifort]
   goto FIN
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd     goto IFORTD
   echo    Unknown compiler!
   goto FIN
rem
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Sup_Exc.f90 /O2 /nologo /heap-arrays:100 /I%CRYSFML%\ifort64\LibC
   ifort /c Simbo.f90   /O2 /nologo /heap-arrays:100 /I%CRYSFML%\ifort64\LibC
   link /subsystem:console  /out:Simbo.exe *.obj %CRYSFML%\ifort64\LibC\CrysFML.lib
   goto END
rem
:IFORTD
   ifort /c Sup_Exc.f90 /debug:full /check /traceback  /nologo  /heap-arrays:100 /I%CRYSFML%\ifort64_debug\LibC 
   ifort /c Simbo.f90   /debug:full /check /traceback  /nologo  /heap-arrays:100 /I%CRYSFML%\ifort64_debug\LibC
   link /subsystem:console /out:Simbo.exe *.obj %CRYSFML%\ifort64_debug\LibC\CrysFML.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Sup_Exc.f90  -I../../GFortran/LibC
   gfortran -c Simbo.f90    -I../../GFortran/LibC
   gfortran *.o -o Simbo    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
rem  Comment the following lines if upx or %FULLPROF% are not available
rem  or if you want to conserve the object files
rem  Compression of executable
rem        upx Simbo.exe
rem  Move the excutable to a directory in the Path
        if exist %FULLPROF% move Simbo.exe %FULLPROF% > nul
rem  Remove unnecessary files
        del *.obj *.mod *.o *.map *.bak > nul
:FIN
