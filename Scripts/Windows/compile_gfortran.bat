@echo off
rem ---------------------------------------
rem ---- CrysFML for GFortran Compiler ----
rem ---- JGP                      2019 ----
rem ---------------------------------------
rem
rem ---- INIT ----
   (set _DEBUG=N)
   (set _WINTER=N)
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
   cd %CRYSFML%\Src
rem
   echo.
   echo **-------------------------------------**
   echo **---- CrysFML: Start Compilation  ----**
   echo **-------------------------------------**
   echo Compiler Options 
   if [%OPTC%]==[-m32] (
      echo GFortran 32 Bits - version %_VER%
   ) else (
      echo GFortran 64 Bits - version %_VER%
   )
   echo OPT0:%OPT0%
   echo OPT1:%OPT1%
   echo OPT2:%OPT2%
   echo OPT3:%OPT3%
   echo.

rem
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   gfortran -c %OPTC% CFML_GlobalDeps_Windows.f90               %OPT1% 
rem
   gfortran -c %OPTC% CFML_math_gen.f90                         %OPT1% 
   gfortran -c %OPTC% CFML_LSQ_TypeDef.f90                      %OPT1% 
   gfortran -c %OPTC% CFML_spher_harm.f90                       %OPT1% 
   gfortran -c %OPTC% CFML_random.f90                           %OPT1% 
   gfortran -c %OPTC% CFML_ffts.f90                             %OPT1% 
   gfortran -c %OPTC% CFML_string_util.f90                      %OPT1% 
rem    gfortran -c %OPTC% CFML_Rational_Arithmetic.f90              %OPT1% %OPT2%
   if [%_WINTER%]==[Y] (
     gfortran -c %OPTC% CFML_io_messwin.f90                     %OPT1% %OPT3%
   ) else (
     gfortran -c %OPTC% CFML_io_mess.f90                        %OPT1% 
   )
   gfortran -c %OPTC% CFML_Profile_TOF.f90                      %OPT1% 
   gfortran -c %OPTC% CFML_Profile_Finger.f90                   %OPT1% 
   gfortran -c %OPTC% CFML_Profile_Functs.f90                   %OPT1% 
rem
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   gfortran -c %OPTC% CFML_math_3D.f90                          %OPT1% 
   gfortran -c %OPTC% CFML_optimization.f90                     %OPT1% 
   gfortran -c %OPTC% CFML_optimization_lsq.f90                 %OPT1% 
   gfortran -c %OPTC% CFML_sym_table.f90                        %OPT0% 
   gfortran -c %OPTC% CFML_chem_scatt.f90                       %OPT0% 
   gfortran -c %OPTC% CFML_BVSpar.f90                           %OPT0% 
   gfortran -c %OPTC% CFML_diffpatt.f90                         %OPT1% 
rem
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   gfortran -c %OPTC% CFML_bonds_table.f90                      %OPT0% 
   gfortran -c %OPTC% CFML_cryst_types.f90                      %OPT1% 
   gfortran -c %OPTC% CFML_symmetry.f90                         %OPT1% 
   gfortran -c %OPTC% -cpp CFML_ILL_Instrm_data.f90             %OPT1% 
rem
   echo .... EoS, Reflections, Atoms
rem
   gfortran -c %OPTC% CFML_Eos_Mod.f90                          %OPT1% 
   gfortran -c %OPTC% CFML_reflct_util.f90                      %OPT1% 
   gfortran -c %OPTC% CFML_atom_mod.f90                         %OPT1% 
rem
   echo .... Formats, Geometry Calculations, Molecules
rem
   gfortran -c %OPTC% CFML_geom_calc.f90                       %OPT1% 
   gfortran -c %OPTC% CFML_molecules.f90                       %OPT1% 
   gfortran -c %OPTC% CFML_form_cif.f90                        %OPT1% 
rem
   echo .... Extinction, Structure Factors, SXTAL geometry, Propag Vectors
rem
   gfortran -c %OPTC% CFML_Extinction_Correction.f90           %OPT1% 
   gfortran -c %OPTC% CFML_sfac.f90                            %OPT1% 
   gfortran -c %OPTC% CFML_sxtal_Geom.f90                      %OPT1% 
   gfortran -c %OPTC% CFML_propagk.f90                         %OPT1% 
rem
   echo .... Maps, BVS, Energy Configurations
rem
   gfortran -c %OPTC% CFML_Export_Vtk.f90                      %OPT1%
   gfortran -c %OPTC% CFML_maps.f90                            %OPT1% 
   gfortran -c %OPTC% CFML_conf_calc.f90                       %OPT1% 
rem
   echo .... Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem
   gfortran -c %OPTC% CFML_Magnetic_Groups.f90                 %OPT1% 
   gfortran -c %OPTC% CFML_magsymm.f90                         %OPT1% 
   gfortran -c %OPTC% CFML_optimization_san.f90                %OPT1%  %OPT3%
   gfortran -c %OPTC% CFML_refcodes.f90                        %OPT1% 
rem
   echo .... Magnetic Structure Factors, Polarimetry
rem
   gfortran -c %OPTC% CFML_msfac.f90                           %OPT1% 
   gfortran -c %OPTC% CFML_polar.f90                           %OPT1% 
rem
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
   echo.
   echo Creating GFORTRAN directory 
rem
   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
   if [%_WINTER%]==[Y] (
     if exist ..\%DIRECTORY%\LibW rmdir ..\%DIRECTORY%\LibW /S /Q
     mkdir ..\%DIRECTORY%\LibW
     copy *.mod ..\%DIRECTORY%\LibW\. > nul
     move *.a ..\%DIRECTORY%\LibW\. > nul
   ) else (
     if exist ..\%DIRECTORY%\LibC rmdir ..\%DIRECTORY%\LibC /S /Q
     mkdir ..\%DIRECTORY%\LibC
     copy *.mod ..\%DIRECTORY%\LibC\. > nul
     move *.a ..\%DIRECTORY%\LibC\. > nul
   )
   del *.o *.mod *.lst *.bak > nul
rem
   echo.
   echo **---- End Compilation for CrysFML ----**
   echo.
   cd %CRYSFML%\Scripts\Windows
rem
:END   
