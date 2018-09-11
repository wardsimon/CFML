@echo off
rem ****
rem ****---- Compilation for Testin_SSG Program ----****
rem ****
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Matrix_Mod.f90                    /O2 /heap-arrays  /nologo /IC:\CrysFML\ifort\LibC
   ifort /c CFML_ssg_datafile.f90             /O2 /heap-arrays  /nologo /IC:\CrysFML\ifort\LibC
   ifort /c CFML_SuperSpaceGroups.f90         /O2 /heap-arrays  /nologo /IC:\CrysFML\ifort\LibC
   ifort /c testing_ssg.f90                   /O2 /heap-arrays  /nologo /IC:\CrysFML\ifort\LibC
   ifort /exe:testing_ssg *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   goto END
:IFORTD
   ifort /c Matrix_Mod.f90                    /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c CFML_ssg_datafile.f90             /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c CFML_SuperSpaceGroups.f90         /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c testing_ssg.f90                   /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /exe:testing_ssg *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Matrix_Mod.f90   -I../../GFortran/LibC
   gfortran -c CFML_ssg_datafile.f90   -I../../GFortran/LibC
   gfortran -c CFML_SuperSpaceGroups.f90   -I../../GFortran/LibC
   gfortran -c testing_ssg.f90   -I../../GFortran/LibC
   gfortran *.o -o testing_ssg    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
