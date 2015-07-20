@echo off
rem ---------------------------------------
rem ---- CrysFML for GFortran Compiler ----
rem ---- JGP-2015                      ----
rem ---------------------------------------
rem
rem ---- INIT ----
   (set _DEBUG=N)
   if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
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
      if [%TARGET_ARCH%]==[ia32] (set DIRECTORY=GFortran_03_debug) else (set DIRECTORY=GFortran64_03_debug)
      if [%TARGET_ARCH%]==[ia32] (
         (set OPT0=-O0 -m32)
         (set OPT1=-O0 -m32)
      ) else (
         (set OPT0=-O0 -m64)
         (set OPT1=-O0 -m64)
      )
      (set OPT2=-fbacktrace -ffree-line-length-none)
   ) else (
      if [%TARGET_ARCH%]==[ia32] (set DIRECTORY=GFortran_03) else (set DIRECTORY=GFortran64_03)
      if [%TARGET_ARCH%]==[ia32] (
         (set OPT0=-O0 -m32)
         (set OPT1=-O3 -m32)
      ) else (
         (set OPT0=-O0 -m64)
         (set OPT1=-O3 -m64)
      )
      (set OPT2=-ffree-line-length-none -funroll-loops -msse2)
      (set OPT4=-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans)
   )
   (set OPT3=)
   if [%_WINTER%]==[Y] (
      if [%TARGET_ARCH%]==[ia32] (set LIBFOR=lib.gnu32) else (set LIBFOR=lib.gnu64)
      (set OPT3=/I%WINTER%\%LIBFOR%)
   )
rem
   cd %CRYSFML%\Src03
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   gfortran -c CFML_GlobalDeps_Windows.f90               %OPT1% %OPT2%
rem
   gfortran -c CFML_math_gen.f90                         %OPT1% %OPT2%
   gfortran -c CFML_LSQ_TypeDef.f90                      %OPT1% %OPT2%
   gfortran -c CFML_spher_harm.f90                       %OPT1% %OPT2%
   gfortran -c CFML_random.f90                           %OPT1% %OPT2%
   gfortran -c CFML_ffts.f90                             %OPT1% %OPT2%
   gfortran -c CFML_string_util.f90                      %OPT1% %OPT2% %OPT4%
   if [%_WINTER%]==[Y] (
     gfortran -c CFML_io_messwin.f90                     %OPT1% %OPT2% %OPT3%
   ) else (
     gfortran -c CFML_io_mess.f90                        %OPT1% %OPT2%
   )
   gfortran -c CFML_Profile_TOF.f90                      %OPT1% %OPT2%
   gfortran -c CFML_Profile_Finger.f90                   %OPT1% %OPT2%
   gfortran -c CFML_Profile_Functs.f90                   %OPT1% %OPT2%
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   gfortran -c CFML_math_3D.f90                          %OPT1% %OPT2%
   gfortran -c CFML_optimization.f90                     %OPT1% %OPT2%
   gfortran -c CFML_optimization_lsq.f90                 %OPT1% %OPT2%
   gfortran -c CFML_sym_table.f90                        %OPT0% %OPT2%
   gfortran -c CFML_chem_scatt.f90                       %OPT0% %OPT2%
   gfortran -c CFML_BVSpar.f90                           %OPT0% %OPT2%
   gfortran -c CFML_diffpatt.f90                         %OPT1% %OPT2%
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   gfortran -c CFML_bonds_table.f90                      %OPT0% %OPT2%
   gfortran -c CFML_cryst_types.f90                      %OPT1% %OPT2%
   gfortran -c CFML_symmetry.f90                         %OPT1% %OPT2%
   gfortran -c CFML_ILL_Instrm_data.f90                  %OPT1% %OPT2%
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms
rem
   gfortran -c CFML_Eos_Mod.f90                          %OPT1% %OPT2%
   gfortran -c CFML_reflct_util.f90                      %OPT1% %OPT2%
   goto FIN
   gfortran -c CFML_atom_mod.f90                         %OPT1% %OPT2%
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, SXTAL geometry, Propag Vectors
rem
   gfortran -c CFML_sfac.f90                            %OPT1% %OPT2%
   gfortran -c CFML_geom_calc.f90                       %OPT1% %OPT2%
   gfortran -c CFML_sxtal_Geom.f90                      %OPT1% %OPT2%
   gfortran -c CFML_propagk.f90                         %OPT1% %OPT2%
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   gfortran -c CFML_Export_Vtk.f90                      %OPT1% %OPT2%
   gfortran -c CFML_maps.f90                            %OPT1% %OPT2%
   gfortran -c CFML_molecules.f90                       %OPT1% %OPT2%
   gfortran -c CFML_conf_calc.f90                       %OPT1% %OPT2%
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   gfortran -c CFML_form_cif.f90                        %OPT1% %OPT2%
rem
   echo **---- Level 7 ----**
   echo .... Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem
   gfortran -c CFML_magsymm.f90                         %OPT1% %OPT2%
   gfortran -c CFML_optimization_san.f90                %OPT1% %OPT2% %OPT3%
   gfortran -c CFML_refcodes.f90                        %OPT1% %OPT2%
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
   gfortran -c CFML_msfac.f90                           %OPT1% %OPT2%
   gfortran -c CFML_polar.f90                           %OPT1% %OPT2%
rem
rem
   echo **---- Crysfml Library ----**
rem
   if [%_WINTER%]==[Y] (
     ar cr libwcrysfml.a *.o
   ) else (
     ar cr libcrysfml.a *.o
   )
rem
   echo **---- GFortran Directory ----**
rem
   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
   if [%_WINTER%]==[Y] (
     if exist ..\%DIRECTORY%\LibW rmdir ..\%DIRECTORY%\LibW /S /Q
     mkdir ..\%DIRECTORY%\LibW
     copy *.mod ..\%DIRECTORY%\LibW > nul
     move *.a ..\%DIRECTORY%\LibW > nul
   ) else (
     if exist ..\%DIRECTORY%\LibC rmdir ..\%DIRECTORY%\LibC /S /Q
     mkdir ..\%DIRECTORY%\LibC
     copy *.mod ..\%DIRECTORY%\LibC > nul
     move *.a ..\%DIRECTORY%\LibC > nul
   )
   del *.o *.mod *.lst *.bak > nul
rem
   cd %CRYSFML%\Scripts\Windows
:FIN
