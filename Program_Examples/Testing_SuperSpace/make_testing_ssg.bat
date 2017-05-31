@echo off
rem ****
rem ****---- Compilation for Testin_SSG Program ----****
rem ****
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   goto END
rem
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c CFML_Rational_Arithmetic_test.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c Matrix_Mod.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c CFML_ssg_datafile.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c CFML_SuperSpaceGroups.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c testing_ssg.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /exe:testing_ssg *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c CFML_Rational_Arithmetic_test.f90   -I../../GFortran/LibC
   gfortran -c Matrix_Mod.f90   -I../../GFortran/LibC
   gfortran -c CFML_ssg_datafile.f90   -I../../GFortran/LibC
   gfortran -c CFML_SuperSpaceGroups.f90   -I../../GFortran/LibC
   gfortran -c testing_ssg.f90   -I../../GFortran/LibC
   gfortran *.o -o testing_ssg    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
