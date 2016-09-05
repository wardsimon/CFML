@echo off
rem ****
rem ****---- Compilation for DIFFaX2FAULTS Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: January 2016
rem ****
rem
rem      This script is valid for both 32 and 64 bits because there is no
rem      dependence on libraries other that those related to the compiler
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_diffax2faults: Make DIFFaX2FAULTS Compilation
   echo    Syntax: make_diffax2faults   [lf95/g95/gfortran/ifort]
   goto FIN
rem
:CONT
   if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
   cd "%CRYSFML%"\Program_Examples\Faults\DIFFaX2FAULTS
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto FIN
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c DIFFaX2FAULTS.f90  -g   -tp -nomap -stchk -chk -o1 -mod .
   lf95  *.obj -out diffax2faults_lf   -tp -nomap -stchk -chk -o1
   upx diffax2faults_lf.exe
   if exist %FULLPROF% copy diffax2faults_lf.exe %FULLPROF% > nul
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
echo off
rem  ifort /c DIFFaX2FAULTS.f90  /O2 /nologo  /heap-arrays:100
rem   link  *.obj /subsystem:console /out:diffax2faults.exe
   ifort /c DIFFaX2FAULTS.f90  /O3 /nologo
   ifort /exe:diffax2faults.exe *.obj  /link /stack:64000000
   upx diffax2faults.exe
   if exist %FULLPROF% copy diffax2faults.exe %FULLPROF% > nul
   if [%TARGET_ARCH%]==[ia32] (
       if exist %PROGCFML% copy diffax2faults.exe %PROGCFML%\DistFPS\. > nul
   ) else (
       if exist %PROGCFML% copy diffax2faults.exe %PROGCFML%\DistFPS_64b\. > nul
   )
   goto END
:IFORTD
echo on
   ifort /c DIFFaX2FAULTS.f90  /debug:full /check /traceback  /nologo  /heap-arrays:100    /I"%CRYSFML%"\ifort_debug\LibC
   link  *.obj /subsystem:console /out:diffax2faultsd.exe
   upx diffax2faultsd.exe
   if exist %FULLPROF% copy diffax2faultsd.exe %FULLPROF% > nul
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3 -funroll-loops  -msse2   DIFFaX2FAULTS.f90
   g95  *.o -o diffax2faults -O3  -funroll-loops  -msse2
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O0 -g -fbounds-check  -fbacktrace    DIFFaX2FAULTS.f90    -I"%CRYSFML%"\GFortran\LibC
   gfortran *.o -o diffax2faults_gf    -static-libgfortran
   upx diffax2faults_gf.exe
   if exist %FULLPROF% copy diffax2faults_gf.exe %FULLPROF% > nul
   goto END
rem
:END
rem
rem
rem
   del *.obj *.mod *.o *.map *.bak > nul
   cd "%CRYSFML%"\Program_Examples\Faults\Scripts\Windows
:FIN
