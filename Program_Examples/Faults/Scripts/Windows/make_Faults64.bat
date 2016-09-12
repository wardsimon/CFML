@echo off
rem ****
rem ****---- Compilation for FAULTS Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: October 2011
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_Faults64: Make FAULTS Compilation
   echo    Syntax: make_Faults64  [lf95/g95/gfortran/ifort]
   goto FIN
rem
:CONT
   cd "%CRYSFML%"\Program_Examples\Faults\Src
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd     goto IFORTD
   goto FIN
rem
rem ****---- Lahey Compiler ----****
:LF95
   echo "Lahey compiler not available for 64 bits"
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
echo off
rem  ifort /c Diffax_glb.f90  /O2 /nologo  /heap-arrays:100       /I"%CRYSFML%"\ifort64\LibC
rem  ifort /c Faults_Read.f90 /O2 /nologo  /heap-arrays:100       /I"%CRYSFML%"\ifort64\LibC
rem  ifort /c Diffax_calc.f90 /O2 /nologo  /heap-arrays:100       /I"%CRYSFML%"\ifort64\LibC
rem  ifort /c Faults.f90      /O2 /nologo  /heap-arrays:100       /I"%CRYSFML%"\ifort64\LibC
rem   link  *.obj /subsystem:console /out:Faults.exe  "%CRYSFML%"\ifort64\LibC\crysfml.lib
  ifort /c Diffax_glb.f90  /O3 /nologo        /I"%CRYSFML%"\ifort64\LibC
  ifort /c Faults_Read.f90 /O3 /nologo        /I"%CRYSFML%"\ifort64\LibC
  ifort /c Diffax_calc.f90 /O3 /nologo        /I"%CRYSFML%"\ifort64\LibC
  ifort /c Faults.f90      /O3 /nologo        /I"%CRYSFML%"\ifort64\LibC
  ifort /exe:Faults.exe *.obj "%CRYSFML%"\ifort64\LibC\CrysFML.lib /link /stack:64000000

   upx Faults.exe
   if exist %FULLPROF% copy Faults.exe %FULLPROF% > nul
   if exist %PROGCFML% copy Faults.exe %PROGCFML%\DistFPS\. > nul
   if exist %CRYSFML%\Program_Examples\Faults\DistFAULTS copy Faults.exe %CRYSFML%\Program_Examples\Faults\DistFAULTS\Windows\. > nul
   goto END
:IFORTD
echo on
   ifort /c Diffax_glb.f90  /debug:full /check /traceback  /nologo  /heap-arrays:100    /I"%CRYSFML%"\ifort_debug\LibC
   ifort /c Faults_Read.f90 /debug:full /check /traceback  /nologo  /heap-arrays:100    /I"%CRYSFML%"\ifort_debug\LibC
   ifort /c Diffax_calc.f90 /debug:full /check /traceback  /nologo  /heap-arrays:100    /I"%CRYSFML%"\ifort_debug\LibC
   ifort /c Faults.f90      /debug:full /check /traceback  /nologo  /heap-arrays:100    /I"%CRYSFML%"\ifort_debug\LibC
   link  *.obj /subsystem:console /out:Faultsd.exe  "%CRYSFML%"\ifort_debug\LibC\crysfml.lib
   upx Faultsd.exe
   if exist %FULLPROF% copy Faultsd.exe %FULLPROF% > nul
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   echo "Lahey compiler not available for 64 bits"
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
rem     gfortran -c -O3 -funroll-loops  -msse2   Diffax_glb.f90    -I"%CRYSFML%"\GFortran64\LibC
rem     gfortran -c -O3 -funroll-loops  -msse2   Faults_Read.f90   -I"%CRYSFML%"\GFortran64\LibC
rem     gfortran -c -O3 -funroll-loops  -msse2   Diffax_calc.f90   -I"%CRYSFML%"\GFortran64\LibC
rem     gfortran -c -O3 -funroll-loops  -msse2   Faults.f90        -I"%CRYSFML%"\GFortran64\LibC
rem     gfortran *.o -o Faults_gf -O3  -funroll-loops  -msse2      -L"%CRYSFML%"\GFortran64\LibC -lcrysfml

   gfortran -c -O0 -g -fbounds-check  -fbacktrace    Diffax_glb.f90    -I"%CRYSFML%"\GFortran64\LibC
   gfortran -c -O0 -g -fbounds-check  -fbacktrace    Faults_Read.f90   -I"%CRYSFML%"\GFortran64\LibC
   gfortran -c -O0 -g           -fbacktrace -legacy  Diffax_calc.f90   -I"%CRYSFML%"\GFortran64\LibC
   gfortran -c -O0 -g -fbounds-check  -fbacktrace    Faults.f90        -I"%CRYSFML%"\GFortran64\LibC
   gfortran *.o -o Faults_gf     -L"%CRYSFML%"\GFortran64\LibC -lcrysfml -static-libgfortran
   upx Faults_gf.exe
   if exist %FULLPROF% copy Faults_gf.exe %FULLPROF% > nul
   goto END
rem
:END
   if exist %PROGCFML% copy Faults.exe %PROGCFML%\DistFPS_64b\. > nul
   del *.obj *.mod *.o *.map *.bak > nul
   cd "%CRYSFML%"\Program_Examples\Faults\Help
   if exist %PROGCFML% copy Faults_Manual.pdf %PROGCFML%\DistFPS\Docs\ > nul
   cd "%CRYSFML%"\Program_Examples\Faults\Scripts\Windows
:FIN
