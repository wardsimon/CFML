@echo off
rem ---------------------------------------
rem ---- CrysFML for GFortran Compiler ----
rem ---- JGP                      2019 ----
rem ---------------------------------------
rem
rem ---- INIT ----
   (set _DEBUG=N)
   (set _VER=8.1)
rem
rem GFortran for 32bits   
   (set OPTC=-m32)
rem
rem GFortran for 64bits   
rem   (set OPTC=-m64)
rem
rem ----
rem ---- Arguments ----
rem ----
:LOOP
    if [%1]==[debug] (set _DEBUG=Y)
    if [%1]==[winter] (set _WINTER=Y)
    if [%1]==[32] (set OPTC=-m32)
    if [%1]==[64] (set OPTC=-m64)
    shift
    if not [%1]==[] goto LOOP
rem
rem ---- Options
rem
   if [%_DEBUG%]==[Y] (
      if [%OPTC%]==[-m32] (set DIRECTORY=gfortran_debug) else (set DIRECTORY=gfortran64_debug)
      (set OPT0=-O0 -std=f2008 -Wall -fdec-math -fbacktrace  -ffree-line-length-0 -fall-intrinsics)
      (set OPT1=-O0 -std=f2008 -Wall -fdec-math -fbacktrace  -ffree-line-length-0 -fall-intrinsics)
   ) else (
      if [%OPTC%]==[-m32] (set DIRECTORY=gfortran) else (set DIRECTORY=gfortran64)
      (set OPT0=-O0 -std=f2008 -ffree-line-length-0 -fdec-math -fall-intrinsics)
      (set OPT1=-O3 -std=f2008 -ffree-line-length-0 -fdec-math -fall-intrinsics)
   )
   (set OPT3=)
   if [%_WINTER%]==[Y] (
      if [%OPTC%]==[-m32] (set LIBFOR=lib.gnu32/%_VER%) else (set LIBFOR=lib.gnu64/%_VER%)
      (set OPT3=-I%WINTER%\%LIBFOR%)
   )
rem
   cd %CRYSFML%\Src08
rem
   echo.
   echo **-------------------------------------**
   echo **---- CrysFML: Start Compilation  ----**
   echo **-------------------------------------**
   echo Compiler Options 
   echo OPT0:%OPT0%
   echo OPT1:%OPT1%
   echo OPT2:%OPT2%
   echo OPT3:%OPT3%
   echo.
rem
   echo .... Global Dependencies for CFML
rem
   gfortran -c %OPTC% -J.\mod CFML_GlobalDeps_Windows_GFOR.f90       %OPT1% 
rem
   echo .... Mathematics Procedures
   gfortran -c %OPTC%  -J.\mod CFML_Maths.f90                        %OPT1% 
rem 
rem   Submodules CFML_Maths   
      cd .\CFML_Maths
      gfortran -c %OPTC%  -J..\mod Co_Prime.f90                      %OPT1% 
      gfortran -c %OPTC%  -J..\mod Co_Linear.f90                     %OPT1% 
      gfortran -c %OPTC%  -J..\mod Cross_Product.f90                 %OPT1% 
      gfortran -c %OPTC%  -J..\mod Debye.f90                         %OPT1% 
      gfortran -c %OPTC%  -J..\mod Determinant.f90                   %OPT1% 
      gfortran -c %OPTC%  -J..\mod Diagonalize_SH.f90                %OPT1% 
      gfortran -c %OPTC%  -J..\mod Equal_Matrix.f90                  %OPT1% 
      gfortran -c %OPTC%  -J..\mod Equal_Vector.f90                  %OPT1% 
      gfortran -c %OPTC%  -J..\mod Factorial.f90                     %OPT1% 
      gfortran -c %OPTC%  -J..\mod Invert_Matrix.f90                 %OPT1% 
      gfortran -c %OPTC%  -J..\mod In_Limits.f90                     %OPT1% 
      gfortran -c %OPTC%  -J..\mod Is_Diagonal_Matrix.f90            %OPT1% 
      gfortran -c %OPTC%  -J..\mod Is_Null_Vector.f90                %OPT1% 
      gfortran -c %OPTC%  -J..\mod Linear_Dependent.f90              %OPT1% 
      gfortran -c %OPTC%  -J..\mod Locate.f90                        %OPT1% 
      gfortran -c %OPTC%  -J..\mod Lower_Triangular.f90              %OPT1% 
      gfortran -c %OPTC%  -J..\mod ManipVec.f90                      %OPT1% 
      gfortran -c %OPTC%  -J..\mod Mat_Cross.f90                     %OPT1% 
      gfortran -c %OPTC%  -J..\mod Modulo_Lat.f90                    %OPT1% 
      gfortran -c %OPTC%  -J..\mod Negligible.f90                    %OPT1% 
      gfortran -c %OPTC%  -J..\mod Norm.f90                          %OPT1% 
      gfortran -c %OPTC%  -J..\mod Outerprod.f90                     %OPT1% 
      gfortran -c %OPTC%  -J..\mod Pgcd.f90                          %OPT1% 
      gfortran -c %OPTC%  -J..\mod Points_In_Line2D.f90              %OPT1% 
      gfortran -c %OPTC%  -J..\mod Poly_Legendre.f90                 %OPT1% 
      gfortran -c %OPTC%  -J..\mod Rank.f90                          %OPT1% 
      gfortran -c %OPTC%  -J..\mod Resolv_System.f90                 %OPT1% 
      gfortran -c %OPTC%  -J..\mod Rotation_Axes.f90                 %OPT1% 
      gfortran -c %OPTC%  -J..\mod RowEchelon.f90                    %OPT1% 
      gfortran -c %OPTC%  -J..\mod Scalar.f90                        %OPT1% 
      gfortran -c %OPTC%  -J..\mod SistCoord_Changes.f90             %OPT1% 
      gfortran -c %OPTC%  -J..\mod Sort.f90                          %OPT1% 
      gfortran -c %OPTC%  -J..\mod Swap.f90                          %OPT1% 
      gfortran -c %OPTC%  -J..\mod Tensor_Product.f90                %OPT1% 
      gfortran -c %OPTC%  -J..\mod Trace.f90                         %OPT1% 
      gfortran -c %OPTC%  -J..\mod Upper_Triangular.f90              %OPT1% 
      gfortran -c %OPTC%  -J..\mod Vec_Length.f90                    %OPT1% 
      gfortran -c %OPTC%  -J..\mod Zbelong.f90                       %OPT1% 
      move /y *.o .. > nul
      cd ..
rem      
   echo .... Strings Procedures
   gfortran -c %OPTC%  -J.\mod  CFML_StringsUtil.f90                 %OPT1%
rem 
rem   Submodules CFML_Strings   
      cd .\CFML_Strings
      
      gfortran -c %OPTC%  -J..\mod StringTools.f90                   %OPT1% 
      gfortran -c %OPTC%  -J..\mod StringNum.f90                     %OPT1% 
      gfortran -c %OPTC%  -J..\mod StringReadKey.f90                 %OPT1% 
      gfortran -c %OPTC%  -J..\mod StringFullp.f90                   %OPT1% 
      move /y *.o .. > nul
      cd ..      
      goto END
rem
:END
rem
   echo.
   echo Creating CrysFML Library 
rem
   if [%_WINTER%]==[Y] (
     ar cr libwcrysfml.a *.o
   ) else (
     ar cr libcrysfml.a *.o
   )
rem
   echo Creating GFORT directory 
rem
   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
   if [%_WINTER%]==[Y] (
     if exist ..\%DIRECTORY%\LibW08 rmdir ..\%DIRECTORY%\LibW08 /S /Q
     mkdir ..\%DIRECTORY%\LibW08
     copy .\mod\*.mod ..\%DIRECTORY%\LibW08 > nul
     copy .\mod\*.smod ..\%DIRECTORY%\LibW08 > nul
     move *.a ..\%DIRECTORY%\LibW08 > nul
   ) else (
     if exist ..\%DIRECTORY%\LibC08 rmdir ..\%DIRECTORY%\LibC08 /S /Q
     mkdir ..\%DIRECTORY%\LibC08
     copy .\mod\*.mod ..\%DIRECTORY%\LibC08 > nul
     copy .\mod\*.smod ..\%DIRECTORY%\LibC08 > nul
     move *.a ..\%DIRECTORY%\LibC08 > nul
   )
   del *.o  *.lst *.bak > nul
rem
   echo.
   echo **---- End Compilation for CrysFML ----**
   echo.
rem

   
   

