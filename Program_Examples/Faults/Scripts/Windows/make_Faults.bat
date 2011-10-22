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
   goto END
rem
:CONT
   cd "%CRYSFML%"\Program_Examples\Faults\Src
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c Diffax_glb.f90   -tp -nomap -stchk -nchk -o1 -mod .;"%CRYSFML%"\Lahey\LibC
   lf95 -c Faults_Read.f90  -tp -nomap -stchk -nchk -o1 -mod .;"%CRYSFML%"\lahey\LibC
   lf95 -c Diffax_calc.f90  -tp -nomap -stchk -nchk -o1 -mod .;"%CRYSFML%"\lahey\LibC
   lf95 -c Faults.f90       -tp -nomap -stchk -nchk -o1 -mod .;"%CRYSFML%"\lahey\LibC
   lf95  *.obj -out Faults  -tp -nomap -stchk -nchk -o1 -lib   "%CRYSFML%"\lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Diffax_glb.f90  /O2 /nologo         /I"%CRYSFML%"\Intel\LibC
   ifort /c Faults_Read.f90 /O2 /nologo         /I"%CRYSFML%"\Intel\LibC
   ifort /c Diffax_calc.f90 /O2 /nologo         /I"%CRYSFML%"\Intel\LibC
   ifort /c Faults.f90      /O2 /nologo         /I"%CRYSFML%"\Intel\LibC
   link  *.obj /subsystem:console /out:Faults.exe  "%CRYSFML%"\Intel\LibC\crysfml.lib
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
   gfortran -c -O3 -funroll-loops  -msse2   Diffax_glb.f90    -I"%CRYSFML%"\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   Faults_Read.f90   -I"%CRYSFML%"\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   Diffax_calc.f90   -I"%CRYSFML%"\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   Faults.f90        -I"%CRYSFML%"\GFortran\LibC
   gfortran *.o -o Faults -O3  -funroll-loops  -msse2         -L"%CRYSFML%"\GFortran\LibC -lcrysfml
   goto END
rem
:END
rem
rem Compress executable
rem
   upx Faults.exe
   if exist %FULLPROF% copy Faults.exe %FULLPROF% > nul
   del *.obj *.mod *.o *.map *.bak > nul
   cd "%CRYSFML%"\Program_Examples\Faults\Scripts\Windows
