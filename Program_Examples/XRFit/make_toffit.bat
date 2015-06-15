@echo off
rem ****
rem ****---- Compilation for TOF-FIT Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: June-2015
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_toffit: Make TOF-FIT Compilation
   echo    Syntax: make_toffit [gfortran/ifort]
   goto END
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c TOF_Diffraction.f90 /O2 /nologo /IC:\CrysFML\ifort\LibC
   ifort /c tof-fit.f90   /O2 /nologo /IC:\CrysFML\ifort\LibC
   link /subsystem:console /stack:102400000 /out:tof-fit.exe *.obj C:\CrysFML\ifort\LibC\CrysFML.lib
   goto END
rem
:IFORTD
   ifort /c TOF_Diffraction.f90 /debug:full /check /traceback  /nologo  /heap-arrays:100    /I"%CRYSFML%"\ifort_debug\LibC
   ifort /c tof-fit.f90  /debug:full /check /traceback  /nologo  /heap-arrays:100    /I"%CRYSFML%"\ifort_debug\LibC
   link /subsystem:console /stack:102400000 /out:tof-fit.exe *.obj C:\CrysFML\ifort_debug\LibC\CrysFML.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c TOF_Diffraction.f90  -I../../GFortran/LibC
   gfortran -c tof-fit.f90    -I../../GFortran/LibC
   gfortran *.o -o tof-fit    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
