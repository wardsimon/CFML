@echo off
rem ****
rem ****---- Compilation for Simple_Calc_Mag_Powder Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2008
rem ****
rem
@echo off
rem ****
rem ****---- Compilation for Simple_Calc_Mag_Powder Program ----****
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
   echo    MAKE_Calc_mag_powder: Make Calc_mag_powder Compilation
   echo    Syntax: make_calc_mag_powder [gfortran/ifort]
   goto END
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Simple_Calc_Mag_Powder.f90 /O2 /nologo %INC%
   link /subsystem:console /out:Simple_Calc_Mag_Powder.exe *.obj %CRYSLIB%
   goto END
:IFORTD
   ifort /c Simple_Calc_Mag_Powder.f90 /O2 /nologo %INCD%
   link /subsystem:console /out:Simple_Calc_Mag_Powder.exe *.obj %CRYSLIBD%
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c Simple_Calc_Mag_Powder.f90    -I../../GFortran/LibC
   gfortran *.o -o Simple_Calc_Mag_Powder    -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
