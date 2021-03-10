@echo off
rem ****
rem ****---- Compilation for ENERMAG Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: Jan-2011
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_EnerMag: Make ENERMAG Compilation
   echo    Syntax: make_EnerMag [f95/lf95/g95/gfortran/ifort]
   goto FIN
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   echo    Unknown compiler!
   goto FIN
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Sup_Exc.f90 /O2 /nologo /I%CRYSFML%\ifort64\LibC
   ifort /c EnerMag.f90   /O2 /nologo /I%CRYSFML%\ifort64\LibC
   link /subsystem:console /stack:64000000 /out:EnerMag.exe *.obj %CRYSFML%\ifort64\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Sup_Exc.f90    -I../../GFortran/LibC
   gfortran -c EnerMag.f90    -I../../GFortran/LibC
   gfortran *.o -o EnerMag    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
rem
rem  Comment the following lines if upx or %FULLPROF% are not available
rem  or if you want to conserve the object files
rem  Compression of executable
        upx EnerMag.exe
rem  Move the excutable to a directory in the Path
        if exist %FULLPROF% move EnerMag.exe %FULLPROF% > nul
rem  Remove unnecessary files
        del *.obj *.mod *.o *.map *.bak > nul
:FIN