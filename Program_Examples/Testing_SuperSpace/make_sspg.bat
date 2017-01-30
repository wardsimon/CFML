@echo off
rem ****
rem ****---- Compilation for Super_Space_Group_Info Program ----****
rem ****
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   goto END
rem
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Super_Space_Group_Info.f90 /O2 /nologo /IC:\CrysFML\ifort\LibC
   ifort /exe:space_group_info *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Super_Space_Group_Info.f90   -I../../GFortran/LibC
   gfortran *.o -o Super_Space_Group_Info    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
