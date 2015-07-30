@echo off
rem -------------------------------------
rem ---- CrysFML for Absoft Compiler ----
rem ---- JGP-2015                    ----
rem -------------------------------------
rem
rem ---- INIT ----
   (set _DEBUG=N)
   (set _WINTER=N)
rem
rem ---- Arguments ----
rem
:LOOP
    if [%1]==[debug] (set _DEBUG=Y)
    if [%1]==[winter] (set _WINTER=Y)
    shift
    if not [%1]==[] goto LOOP
rem
rem ---- Options
rem
   if [%_DEBUG%]==[Y] (
      (set DIRECTORY=Absoft_debug)
      (set OPT0=-g)
      (set OPT1=-g)
   ) else (
      (set DIRECTORY=Absoft)
      (set OPT0=-O0 -w)
      (set OPT1=-O3 -w)
   )
   if [%_WINTER%]==[Y] (
      (set LIBFOR=)
      (set OPT3=/I%WINTER%\%LIBFOR%)
   )
rem
   cd %CRYSFML%\Src
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   f95 -c CFML_GlobalDeps_Windows.f90    %OPT1%
rem
   f95 -c CFML_math_gen.f90              %OPT1%
   f95 -c CFML_LSQ_TypeDef.f90           %OPT1%
   f95 -c CFML_spher_harm.f90            %OPT1%
   f95 -c CFML_random.f90                %OPT1%
   f95 -c CFML_ffts.f90                  %OPT1%
   f95 -c CFML_string_util.f90           %OPT1%
   if [%_WINTER%]==[Y] (
     f95 -c CFML_io_messwin.f90          %OPT1% %OPT3%
   ) else (
     f95 -c CFML_io_mess.f90             %OPT1%
   )
   f95 -c CFML_io_mess.f90               %OPT1%
   f95 -c CFML_Profile_TOF.f90           %OPT1%
   f95 -c CFML_Profile_Finger.f90        %OPT1%
   f95 -c CFML_Profile_Functs.f90        %OPT1%
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   f95 -c CFML_math_3D.f90            %OPT1%
   f95 -c CFML_optimization.f90       %OPT1%
   f95 -c CFML_optimization_lsq.f90   %OPT1%
   f95 -c CFML_sym_table.f90          %OPT0%
   f95 -c CFML_chem_scatt.f90         %OPT0%
   f95 -c CFML_diffpatt.f90           %OPT1%
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   f95 -c CFML_bonds_table.f90        %OPT0%
   f95 -c CFML_cryst_types.f90        %OPT1%
   f95 -c CFML_symmetry.f90           %OPT1%
   f95 -c CFML_ILL_Instrm_data.f90    %OPT1%
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms
rem
   f95 -c CFML_Eos_Mod.f90            %OPT1%
   f95 -c CFML_reflct_util.f90        %OPT1%
   f95 -c CFML_atom_mod.f90           %OPT1%
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, SXTAL geometry, Propag Vectors
rem
   f95 -c CFML_sfac.f90               %OPT1%
   f95 -c CFML_geom_calc.f90          %OPT1%
   f95 -c CFML_SXTAL_geom.f90         %OPT1%
   f95 -c CFML_propagk.f90            %OPT1%
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   f95 -c CFML_Export_Vtk.f90         %OPT1%
   f95 -c CFML_maps.f90               %OPT1%
   f95 -c CFML_molecules.f90          %OPT1%
   f95 -c CFML_conf_calc.f90          %OPT1%
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   f95 -c CFML_form_cif.f90           %OPT1%
rem
   echo **---- Level 7 ----**
   echo .... Keywords Parser, Simulated Annealing, Magnetic Symmetry
rem
   f95 -c CFML_magsymm.f90            %OPT1%
   f95 -c CFML_optimization_san.f90   %OPT1%  %OPT3%
   f95 -c CFML_refcodes.f90           %OPT1%
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
   f95 -c CFML_msfac.f90        %OPT1%
   f95 -c CFML_polar.f90        %OPT1%
rem
rem
   echo **---- Crysfml Library ----**
rem
   if [%_WINTER%]==[Y] (
     lib /out:wcrysfml.lib *.obj
   ) else (
     lib /out:crysfml.lib *.obj
   )
rem
   echo **---- ifort Directory ----**
rem
   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
   if [%_WINTER%]==[Y] (
     if exist ..\%DIRECTORY%\LibW rmdir ..\%DIRECTORY%\LibW /S /Q
     mkdir ..\%DIRECTORY%\LibW
     copy *.mod ..\%DIRECTORY%\LibW > nul
     move *.lib ..\%DIRECTORY%\LibW > nul
   ) else (
     if exist ..\%DIRECTORY%\LibC rmdir ..\%DIRECTORY%\LibC /S /Q
     mkdir ..\%DIRECTORY%\LibC
     copy *.mod ..\%DIRECTORY%\LibC > nul
     move *.lib ..\%DIRECTORY%\LibC > nul
   )
   del *.obj *.mod *.lst *.bak > nul
rem
   cd %CRYSFML%\Scripts\Windows


