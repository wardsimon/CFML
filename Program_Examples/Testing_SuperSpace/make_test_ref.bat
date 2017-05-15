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
   ifort /c CFML_ssg_datafile.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c CFML_SuperSpaceGroups.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c test_ref.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /exe:test_ref *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c CFML_ssg_datafile.f90   -I../../GFortran/LibC
   gfortran -c CFML_SuperSpaceGroups.f90   -I../../GFortran/LibC
   gfortran -c test_ref.f90   -I../../GFortran/LibC
   gfortran *.o -o test_ref    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
