@echo off
rem ****
rem ****---- Compilation for CRYSCALCON Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: July-2014
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_CrysCalCon: Make CrysCalcon Compilation
   echo    Syntax: make_CrysCalCon [f95/lf95/g95/gfrotran/ifort]
   goto END
rem
:CONT
   if x%1 == xf95       goto F95
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd     goto IFORTd
   if x%1 == xgfortrand  goto GFORd
   goto END
rem
rem ****---- Absoft Compiler ----****
:F95
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c menu_0.f90  -tp -nomap -stchk -nchk -o1 -mod ".;..\..\Lahey\LibC"
   lf95 -c menu_1.f90  -tp -nomap -stchk -nchk -o1 -mod ".;..\..\Lahey\LibC"
   lf95 -c menu_2.f90  -tp -nomap -stchk -nchk -o1 -mod ".;..\..\lahey\LibC"
   lf95 -c menu_3.f90  -tp -nomap -stchk -nchk -o1 -mod ".;..\..\lahey\LibC"
   lf95 -c menu_4.f90  -tp -nomap -stchk -nchk -o1 -mod ".;..\..\lahey\LibC"
   lf95 -c menu_5.f90  -tp -nomap -stchk -nchk -o1 -mod ".;..\..\lahey\LibC"
   lf95 -c calsym.f90  -tp -nomap -stchk -nchk -o1 -mod ".;..\..\lahey\LibC"
   lf95  *.obj -out CrysCalCon -tp -nomap -stchk -nchk -o1 -lib ..\..\lahey\LibC\CrysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
rem   Directory holding the source code of the FullProf library for writing PCRs
    if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
    if [%TARGET_ARCH%]==[ia32]  (
         set INC=/I"%CRYSFML%"\ifort\libC
         set CRYSLIB="%CRYSFML%"\ifort\libC\crysfml.lib
         set DistFPS=%PROGCFML%\DistFPS
      ) else (
         set INC=/I"%CRYSFML%"\ifort64\libC
         set CRYSLIB="%CRYSFML%"\ifort64\libC\crysfml.lib
         set DistFPS=%PROGCFML%\DistFPS_64b
      )
   ifort /c menu_0.f90 /O3 /nologo %INC%
   ifort /c menu_1.f90 /O3 /nologo %INC%
   ifort /c menu_2.f90 /O3 /nologo %INC%
   ifort /c menu_3.f90 /O3 /nologo %INC%
   ifort /c menu_4.f90 /O3 /nologo %INC%
   ifort /c menu_5.f90 /O3 /nologo %INC%
   ifort /c calsym.f90 /O3 /nologo %INC%
   rem ifort /exe:CrysCalcon *.obj %CRYSLIB%
   link /subsystem:console /out:CrysCalCon.exe *.obj %CRYSLIB%
   goto END
:IFORTD
    if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
    if [%TARGET_ARCH%]==[ia32]  (
         set INC=/I"%CRYSFML%"\ifort_debug\libC
         set CRYSLIB="%CRYSFML%"\ifort_debug\libC\crysfml.lib
         set DistFPS=%PROGCFML%\DistFPS
      ) else (
         set INC=/I"%CRYSFML%"\ifort64_debug\libC
         set CRYSLIB="%CRYSFML%"\ifort64_debug\libC\crysfml.lib
         set DistFPS=%PROGCFML%\DistFPS_64b
      )
   ifort /c menu_0.f90 /debug:full /check /traceback /nologo %INC%
   ifort /c menu_1.f90 /debug:full /check /traceback /nologo %INC%
   ifort /c menu_2.f90 /debug:full /check /traceback /nologo %INC%
   ifort /c menu_3.f90 /debug:full /check /traceback /nologo %INC%
   ifort /c menu_4.f90 /debug:full /check /traceback /nologo %INC%
   ifort /c menu_5.f90 /debug:full /check /traceback /nologo %INC%
   ifort /c calsym.f90 /debug:full /check /traceback /nologo %INC%
   rem ifort /exe:CrysCalcon *.obj %CRYSLIB%
   link /subsystem:console /out:CrysCalCon.exe *.obj %CRYSLIB%
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3 -funroll-loops  -msse2   menu_0.f90     -I..\..\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   menu_1.f90     -I..\..\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   menu_2.f90     -I..\..\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   menu_3.f90     -I..\..\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   menu_4.f90     -I..\..\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   menu_5.f90     -I..\..\G95\LibC
   g95 -c -O3 -funroll-loops  -msse2   calsym.f90     -I..\..\G95\LibC
   g95  *.o -o cryscalcon -O3  -funroll-loops  -msse2  -L..\..\G95\LibC -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3 -funroll-loops  -msse2   menu_0.f90     -I..\..\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   menu_1.f90     -I..\..\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   menu_2.f90     -I..\..\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   menu_3.f90     -I..\..\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   menu_4.f90     -I..\..\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   menu_5.f90     -I..\..\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   calsym.f90     -I..\..\GFortran\LibC
   gfortran *.o -o cryscalcon -O3  -funroll-loops  -msse2  -L..\..\GFortran\LibC -lcrysfml
   goto END
:GFORD
   gfortran -c -O0 -fbacktrace -funroll-loops  -msse2   menu_0.f90     -I..\..\GFortran_debug\LibC
   gfortran -c -O0 -fbacktrace -funroll-loops  -msse2   menu_1.f90     -I..\..\GFortran_debug\LibC
   gfortran -c -O0 -fbacktrace -funroll-loops  -msse2   menu_2.f90     -I..\..\GFortran_debug\LibC
   gfortran -c -O0 -fbacktrace -funroll-loops  -msse2   menu_3.f90     -I..\..\GFortran_debug\LibC
   gfortran -c -O0 -fbacktrace -funroll-loops  -msse2   menu_4.f90     -I..\..\GFortran_debug\LibC
   gfortran -c -O0 -fbacktrace -funroll-loops  -msse2   menu_5.f90     -I..\..\GFortran_debug\LibC
   gfortran -c -O0 -fbacktrace -funroll-loops  -msse2   calsym.f90     -I..\..\GFortran_debug\LibC
   gfortran *.o -o cryscalcon -fbacktrace -funroll-loops  -msse2  -L..\..\GFortran_debug\LibC -lcrysfml
   goto END
rem
:END
   upx CrysCalCon.exe
   if exist %FULLPROF% copy CrysCalCon.exe %FULLPROF% > nul
   if exist %PROGCFML% copy CrysCalCon.exe %DistFPS%\. > nul
   del *.obj *.mod *.o *.map *.bak > nul
