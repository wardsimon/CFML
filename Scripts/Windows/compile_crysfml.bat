@echo off
rem --------------------------------------------------
rem ---- CRYSTALLOGRAPHIC FORTRAN MODULES LIBRARY ----
rem --------------------------------------------------
rem
rem ---- Define Compiler ----
    if [%1]==[] goto FINAL
    (set _compiler=%1)
rem
rem ---- Options ----
    if [%_compiler%]==[ifort] (
       (set _cop=/c)
       (set _opt0=/Od /nologo /Qvec-report0)
       (set _opt1=/O2 /nologo /Qvec-report0)
rem       (set _opt0=/Od /nologo /Qopt-report:0)
rem       (set _opt1=/O2 /nologo /Qopt-report:0)
      )
rem
rem ---- Change Directory
rem
    cd %CRYSFML%\Src
rem
rem ---- Header
rem
    echo **---------------------------------------------------**
    echo **----                                           ----**
    echo **---- CRYSTALLOGRAPHIC FORTRAN MODULES LIBRARY  ----**
    echo **---- CrysFML Team                   1999-2015  ----**
    echo **---------------------------------------------------**
rem
rem ---- Compilation Zone
rem
    echo (
    echo **---- Level 0 ----**
    echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
    %_compiler% %_cop% CFML_GlobalDeps_Windows_intel.f90   %_opt1%
rem
    %_compiler% %_cop% CFML_math_gen.f90                   %_opt1%
    %_compiler% %_cop% CFML_LSQ_TypeDef.f90                %_opt1%
    %_compiler% %_cop% CFML_spher_harm.f90                 %_opt1%
    %_compiler% %_cop% CFML_random.f90                     %_opt1%
    %_compiler% %_cop% CFML_ffts.f90                       %_opt1%
    %_compiler% %_cop% CFML_string_util.f90                %_opt1%
    %_compiler% %_cop% CFML_io_mess.f90                    %_opt1%
    %_compiler% %_cop% CFML_Profile_TOF.f90                %_opt1%
    %_compiler% %_cop% CFML_Profile_Finger.f90             %_opt1%
    %_compiler% %_cop% CFML_Profile_Functs.f90             %_opt1%
rem
    echo (
    echo **---- Level 1 ----**
    echo .... Mathematical(II), Optimization, Tables, Patterns
rem
    %_compiler% %_cop% CFML_math_3D.f90                    %_opt1%
    %_compiler% %_cop% CFML_optimization.f90               %_opt1%
    %_compiler% %_cop% CFML_optimization_lsq.f90           %_opt1%
    %_compiler% %_cop% CFML_sym_table.f90                  %_opt0%
    %_compiler% %_cop% CFML_chem_scatt.f90                 %_opt0%
    %_compiler% %_cop% CFML_BVSpar.f90                     %_opt0%
    %_compiler% %_cop% CFML_diffpatt.f90                   %_opt1%
rem
    echo (
    echo **---- Level 2 ----**
    echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
    %_compiler% %_cop% CFML_bonds_table.f90                %_opt0%
    %_compiler% %_cop% CFML_cryst_types.f90                %_opt1%
    %_compiler% %_cop% CFML_symmetry.f90                   %_opt1%
    %_compiler% %_cop% CFML_ILL_Instrm_data.f90            %_opt1%
rem
    echo (
    echo **---- Level 3 ----**
    echo .... Reflections, Atoms
rem
    %_compiler% %_cop% CFML_Eos_Mod.f90                    %_opt1%
    %_compiler% %_cop% CFML_reflct_util.f90                %_opt1%
    %_compiler% %_cop% CFML_atom_mod.f90                   %_opt1%
rem
    echo (
    echo **---- Level 4 ----**
    echo .... Structure Factors, Geometry Calculations, SXTAL geometry, Propag Vectors
rem
    %_compiler% %_cop% CFML_sfac.f90                       %_opt1%
    %_compiler% %_cop% CFML_geom_calc.f90                  %_opt1%
    %_compiler% %_cop% CFML_sxtal_Geom.f90                 %_opt1%
    %_compiler% %_cop% CFML_propagk.f90                    %_opt1%
rem
    echo (
    echo **---- Level 5 ----**
    echo .... Molecules, Maps, BVS, Energy Configurations
rem
    %_compiler% %_cop% CFML_Export_Vtk.f90                 %_opt1%
    %_compiler% %_cop% CFML_maps.f90                       %_opt1%
    %_compiler% %_cop% CFML_molecules.f90                  %_opt1%
    %_compiler% %_cop% CFML_conf_calc.f90                  %_opt1%
rem
    echo (
    echo **---- Level 6 ----**
    echo .... Formats
rem
    %_compiler% %_cop% CFML_form_cif.f90                   %_opt1%
rem
    echo (
    echo **---- Level 7 ----**
    echo .... Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem
    %_compiler% %_cop% CFML_magsymm.f90                    %_opt1%
    %_compiler% %_cop% CFML_optimization_san.f90           %_opt1%
    %_compiler% %_cop% CFML_refcodes.f90                   %_opt1%
rem
    echo (
    echo **---- Level 8 ----**
    echo .... Magnetic Structure Factors, Polarimetry
rem
    %_compiler% %_cop% CFML_msfac.f90                      %_opt1%
    %_compiler% %_cop% CFML_polar.f90                      %_opt1%
rem
rem ---- CrysFML Library ----
rem
    echo (
    lib /out:crysfml.lib *.obj
rem
rem ---- Installation Directory
rem
   if [%_compiler%]==[ifort] (
      if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
      if [%TARGET_ARCH%]==[ia32] (set DIRECTORY=ifort) else (set DIRECTORY=ifort64)
   )
rem
   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
   if exist ..\%DIRECTORY%\LibC rmdir ..\%DIRECTORY%\LibC /S /Q
   mkdir ..\%DIRECTORY%\LibC
rem
   copy *.mod ..\%DIRECTORY%\LibC > nul
   move *.lib ..\%DIRECTORY%\LibC > nul
   del *.obj *.mod *.lst *.bak > nul
rem
   cd %CRYSFML%\Scripts\Windows
:FINAL