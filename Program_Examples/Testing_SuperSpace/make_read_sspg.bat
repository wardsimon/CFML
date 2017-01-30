@echo off
rem ****
rem ****---- Compilation for Read_SSPG Program ----****
rem ****
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   goto END
rem
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c read_ssg_datafile.f90 /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /exe:read_ssg_datafile *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c read_ssg_datafile.f90   -I../../GFortran/LibC
   gfortran *.o -o read_ssg_datafile    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
