@echo off
rem ****
rem ****---- Compilation for IO_FILES Program ----****
rem ****
   (set _DEBUG=Y)
rem ****
rem **** COMP=ifort|gfortran
   (set _COMP=ifort)
rem **** VER=m32|m64  for gfortran
   (set _VER=)
rem > Compilers
   if [%_COMP%]==[ifort] (
      if [%_DEBUG%]==[Y] (
         (set OPT0=/debug:full /check /check:noarg_temp_created /traceback /nologo /CB)
         (set OPT1=/debug:full /check /check:noarg_temp_created /traceback /nologo /CB)
      ) else (
         (set OPT0=/Od)
         (set OPT1=/Ox)
      )
      (set OPT2=/fpp /Qopt-report:0)
   )
rem
   if [%_COMP%]==[gfortran] (
      if [%_DEBUG%]==[Y] (
         (set OPT0=-g -O0 -std=f2008 -Wall -fdec-math -fbacktrace  -ffree-line-length-0 -fall-intrinsics)
         (set OPT1=-g -O0 -std=f2008 -Wall -fdec-math -fbacktrace  -ffree-line-length-0 -fall-intrinsics)
      ) else (
         (set OPT0=-O0 -std=f2008 -ffree-line-length-0 -fdec-math -fall-intrinsics)
         (set OPT1=-O3 -std=f2008 -ffree-line-length-0 -fdec-math -fall-intrinsics)
      )
      (set OPT2=)
   )
rem > Compilation
   if [%_COMP%]==[ifort] (
      ifort /c io_files.f90           /heap-arrays  /nologo %OPT1% /I%CRYSFML%\ifort\LibC08
      ifort /exe:io_files *.obj  %CRYSFML%\ifort\LibC08\crysfml.lib /link /stack:300000000
   )
rem
   if [%_COMP%]==[gfortran] (
      gfortran -c io_files.f90             %OPT1% -I%CRYSFML%\gfortran\LibC08
      gfortran -o io_files.exe *.o -L%CRYSFML%\gfortran\LibC08 -lcrysfml
   )
rem
   del *.obj *.mod *.o *.map *.bak > nul
