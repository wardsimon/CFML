@echo off
rem ****
rem ****---- Compilation for Test230_groups Program ----****
rem ****
rem **** Author: JRC + JGP
rem **** Revision: Nov-2018
rem ****
rem
    if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
    if [%TARGET_ARCH%]==[ia32]  (
         set      INC=/I"%CRYSFML%"\ifort\LibC08
         set     INCD=/I"%CRYSFML%"\ifort_debug\LibC08
         set  CRYSLIB="%CRYSFML%"\ifort\LibC08\crysfml.lib
         set CRYSLIBD="%CRYSFML%"\ifort_debug\libC08\crysfml.lib
      ) else (
         set      INC=/I"%CRYSFML%"\ifort64\LibC08
         set     INCD=/I"%CRYSFML%"\ifort64_debug\LibC08
         set  CRYSLIB="%CRYSFML%"\ifort64\LibC08\crysfml.lib
         set CRYSLIBD="%CRYSFML%"\ifort64_debug\libC08\crysfml.lib
      )
   if not x%1 == x goto CONT
   cls
   echo    Testing Standard Space groups
   echo    Syntax: make_grp230 [gfortran[d]|ifort[d]] 
   goto END
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xgfortrand goto GFORD
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem ****---- Intel Compiler ----****
:IFORT

   ifort /c groups_230.f90          /O3  /nologo %INC%  /warn
   ifort /exe:grp_230 *.obj %CRYSLIB% /link /stack:256000000
rem   ifort /exe:grp_230 *.obj %CRYSLIB% 
   goto END
:IFORTD
   ifort /c groups_230.f90         /check:all /debug:full /check:noarg_temp_created /traceback  /nologo  %INCD% /warn
   ifort /exe:grp_230 *.obj %CRYSLIBD% /link /stack:256000000
rem   ifort /exe:grp_230 *.obj %CRYSLIBD% 
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3 CFML_Symmetry_Table.f90    -fbounds-check -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c -O3 CFML_Symmetry08.f90        -fbounds-check -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c -O3 CFML_Magnetic_Groups08.f90 -fbounds-check -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c -O3 CFML_Rational_Groups08.f90 -fbounds-check -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c -O3 CFML_Standard_Settings.f90 -fbounds-check -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c -O3 test230_groups.f90         -fbounds-check -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran *.o -o test230_groups_gf    -L../../GFortran/LibC08   -lcrysfml
   goto END
:GFORD
   gfortran -c CFML_Symmetry_Table.f90    -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c CFML_Symmetry08.f90        -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c CFML_Magnetic_Groups08.f90 -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c CFML_Rational_Groups08.f90 -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c CFML_Standard_Settings.f90 -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran -c test230_groups.f90         -g -fbounds-check -fbacktrace -ffree-line-length-0  -I../../GFortran/LibC08
   gfortran *.o -o test230_groups_gf    -L../../GFortran/LibC08   -lcrysfml
   goto END
rem
:END
rem   del *.obj *.mod *.o *.map *.bak *.pdb > nul
