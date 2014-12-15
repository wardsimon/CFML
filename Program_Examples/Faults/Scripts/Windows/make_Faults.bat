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
   echo    MAKE_Faults: Make FAULTS Compilation
   echo    Syntax: make_Faults  [lf95/g95/gfortran/ifort]
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
   lf95 -c Diffax_glb.f90  -g   -tp -nomap -stchk -chk -o1 -mod .;"%CRYSFML%"\Lahey\LibC
   lf95 -c Faults_Read.f90 -g   -tp -nomap -stchk -chk -o1 -mod .;"%CRYSFML%"\lahey\LibC
   lf95 -c Diffax_calc.f90 -g   -tp -nomap -stchk -chk -o1 -mod .;"%CRYSFML%"\lahey\LibC
   lf95 -c Faults.f90      -g   -tp -nomap -stchk -chk -o1 -mod .;"%CRYSFML%"\lahey\LibC
   lf95  *.obj -out Faults_lf   -tp -nomap -stchk -chk -o1 -lib   "%CRYSFML%"\lahey\LibC\CrysFML
   upx Faults_lf.exe
   if exist %FULLPROF% copy Faults_lf.exe %FULLPROF% > nul
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
echo on
  ifort /c Diffax_glb.f90  /O2 /nologo  /heap-arrays:100       /I"%CRYSFML%"\ifort\LibC
  ifort /c Faults_Read.f90 /O2 /nologo  /heap-arrays:100       /I"%CRYSFML%"\ifort\LibC
  ifort /c Diffax_calc.f90 /O2 /nologo  /heap-arrays:100       /I"%CRYSFML%"\ifort\LibC
  ifort /c Faults.f90      /O2 /nologo  /heap-arrays:100       /I"%CRYSFML%"\ifort\LibC
  link  *.obj /subsystem:console /out:Faults.exe  "%CRYSFML%"\ifort\LibC\crysfml.lib
   upx Faults.exe
   if exist %FULLPROF% copy Faults.exe %FULLPROF% > nul
   if exist %PROGCFML% copy Faults.exe %PROGCFML%\DistFPS\. > nul
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
   g95 -c -O3 -funroll-loops  -msse2   Diffax_glb.f90   -I"%CRYSFML%"\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   Faults_Read.f90  -I"%CRYSFML%"\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   Diffax_calc.f90  -I"%CRYSFML%"\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   Faults.f90       -I"%CRYSFML%"\G95\LibC
   g95  *.o -o Faults -O3  -funroll-loops  -msse2       -L"%CRYSFML%"\G95\LibC -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
rem     gfortran -c -O3 -funroll-loops  -msse2   Diffax_glb.f90    -I"%CRYSFML%"\GFortran\LibC
rem     gfortran -c -O3 -funroll-loops  -msse2   Faults_Read.f90   -I"%CRYSFML%"\GFortran\LibC
rem     gfortran -c -O3 -funroll-loops  -msse2   Diffax_calc.f90   -I"%CRYSFML%"\GFortran\LibC
rem     gfortran -c -O3 -funroll-loops  -msse2   Faults.f90        -I"%CRYSFML%"\GFortran\LibC
rem     gfortran *.o -o Faults_gf -O3  -funroll-loops  -msse2      -L"%CRYSFML%"\GFortran\LibC -lcrysfml

   gfortran -c -O0 -g -fbounds-check  -fbacktrace    Diffax_glb.f90    -I"%CRYSFML%"\GFortran\LibC
   gfortran -c -O0 -g -fbounds-check  -fbacktrace    Faults_Read.f90   -I"%CRYSFML%"\GFortran\LibC
   gfortran -c -O0 -g  -fbacktrace  -legacy  Diffax_calc.f90   -I"%CRYSFML%"\GFortran\LibC
   gfortran -c -O0 -g -fbounds-check  -fbacktrace    Faults.f90        -I"%CRYSFML%"\GFortran\LibC
   gfortran *.o -o Faults_gf     -L"%CRYSFML%"\GFortran\LibC -lcrysfml
   upx Faults_gf.exe
   if exist %FULLPROF% copy Faults_gf.exe %FULLPROF% > nul
   goto END
rem
:END
rem
rem Compress executable
rem
   del *.obj *.mod *.o *.map *.bak > nul
   cd "%CRYSFML%"\Program_Examples\Faults\Scripts\Windows
:FIN
