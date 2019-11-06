@echo off
rem ****
rem ****---- Compilation for get_UB_from2ref program--****
rem ****
rem **** Author: JRC
rem **** Revision: March-2019
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
   echo    MAKE_GetUB: Makes get_UB_from2ref Compilation
   echo    Syntax: make_GetUB [gfortran/ifort]
   goto END
rem
:CONT
   if x%1 == xgfortran   goto GFOR
   if x%1 == xgfortrand  goto GFORD
   if x%1 == xifort      goto IFORT
   if x%1 == xifortd     goto IFORTD
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c get_UB_from2ref.f90   /O2 /nologo %INC%
   ifort /exe:get_UB_from2ref *.obj %CRYSLIB%
   goto END

:IFORTD
   ifort /c get_UB_from2ref.f90   /check:all /debug:full /check:noarg_temp_created /traceback  /nologo  /heap-arrays:100 %INCD%   /warn
   ifort /exe:get_UB_from2ref *.obj %CRYSLIBD%
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c get_UB_from2ref.f90    -fbounds-check -ffree-line-length-0  -I../../GFortran/LibC
   gfortran *.o -o get_UB_from2ref    -L../../GFortran/LibC   -lcrysfml
   goto END
:GFORD
   gfortran -c get_UB_from2ref.f90    -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC
   gfortran *.o -o get_UB_from2ref    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
