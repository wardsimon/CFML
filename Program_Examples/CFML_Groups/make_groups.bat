@echo off
rem ****
rem ****---- Compilation for SIMILAR Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
    if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
    if [%TARGET_ARCH%]==[ia32]  (
         set  INC=/I"%CRYSFML%"\ifort\LibC
         set  INCD=/I"%CRYSFML%"\ifort_debug\LibC
         set  CRYSLIB="%CRYSFML%"\ifort\LibC\crysfml.lib
         set CRYSLIBD="%CRYSFML%"\ifort_debug\libC\crysfml.lib
      ) else (
         set  INC=/I"%CRYSFML%"\ifort64\LibC
         set  INCD=/I"%CRYSFML%"\ifort64_debug\LibC
         set  CRYSLIB="%CRYSFML%"\ifort64\LibC\crysfml.lib
         set CRYSLIBD="%CRYSFML%"\ifort64_debug\libC\crysfml.lib
      )
   if not x%1 == x goto CONT
   cls
   echo    MAKE_GROUPS: Testing CFML_Groups
   echo    Syntax: make_groups [gfortran/ifort] [deb]
   goto END
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xgfortrand  goto GFORD
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem ****---- Intel Compiler ----****
:IFORT

   ifort /c CFML_Rational_Groups.f90  /O2 /Qparallel /nologo %INC%
   ifort /c CFML_Standard_Sett.f90    /O2 /Qparallel /nologo %INC%
   ifort /c groups.f90                /O2 /Qparallel /nologo %INC%
   ifort /exe:groups *.obj %CRYSLIB% /link /stack:64000000
   goto END
:IFORTD
   ifort /c CFML_Rational_Groups.f90  /debug:full /check:noarg_temp_created /traceback  /nologo  /heap-arrays:100 %INCD%
   ifort /c CFML_Standard_Sett.f90  /debug:full /check:noarg_temp_created /traceback  /nologo  /heap-arrays:100 %INCD%
   ifort /c groups.f90       /debug:full /check:noarg_temp_created /traceback  /nologo  /heap-arrays:100 %INCD%
   ifort /exe:groups *.obj %CRYSLIBD%
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c CFML_Rational_Groups.f90  -ffree-line-length-0  -I../../GFortran/LibC
   gfortran -c CFML_Standard_Sett.f90  -ffree-line-length-0  -I../../GFortran/LibC
   gfortran -c groups.f90       -ffree-line-length-0  -I../../GFortran/LibC
   gfortran *.o -o groups    -L../../GFortran/LibC   -lcrysfml
   goto END
:GFORD
   gfortran -c CFML_Rational_Groups.f90 -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC
   gfortran -c CFML_Standard_Sett.f90 -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC
   gfortran -c groups.f90    -g -fbounds-check -fbacktrace  -ffree-line-length-0  -I../../GFortran/LibC
   gfortran *.o -o groups    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak *.pdb > nul
