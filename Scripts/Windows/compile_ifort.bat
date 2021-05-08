@echo off
rem ------------------------------------
rem ---- CrysFML for Intel Compiler ----
rem ---- JGP                   2019 ----
rem ------------------------------------
rem
rem ---- INIT ----
   (set _DEBUG=N)
   (set _WINTER=N)
   if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
rem
rem ----
rem ---- Arguments ----
rem ----
:LOOP
    if [%1]==[debug]  (set _DEBUG=Y)
    if [%1]==[winter] (set _WINTER=Y)
    shift
    if not [%1]==[] goto LOOP
rem
rem ---- Options
rem
   if [%_DEBUG%]==[Y] (
      if [%TARGET_ARCH%]==[ia32] (set DIRECTORY=ifort_debug) else (set DIRECTORY=ifort64_debug)
      (set OPT0=/debug:full /check /check:noarg_temp_created /traceback /nologo /CB)
      (set OPT1=/debug:full /check /check:noarg_temp_created /traceback /nologo /CB)
   ) else (
      if [%TARGET_ARCH%]==[ia32] (set DIRECTORY=ifort) else (set DIRECTORY=ifort64)
      (set OPT0=/Od)
      (set OPT1=/O2)
   )
rem
   (set OPT2=/fpp /Qopt-report:0)
   (set OPT3=)
   if [%_WINTER%]==[Y] (
      if [%TARGET_ARCH%]==[ia32] (
         (set OPT3=/I%WINTER%\lib.if8)
      ) else (
         (set OPT3=/I%WINTER%\lib.i64)
      )
   )
rem
   cd %CRYSFML%\Src
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
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   ifort /c CFML_GlobalDeps_Windows_intel.f90         /nologo %OPT1% %OPT2%
rem
   ifort /c CFML_Math_General.f90                     /nologo %OPT1% %OPT2%
rem NEW
   ifort /c CFML_EisPack.f90                          /nologo %OPT1% %OPT2%
rem NEW
   ifort /c CFML_LSQ_TypeDef.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_Spherical_Harmonics.f90              /nologo %OPT1% %OPT2%
   ifort /c CFML_Random_Generators.f90                /nologo %OPT1% %OPT2%
   ifort /c CFML_FFT.f90                              /nologo %OPT1% %OPT2%
   ifort /c CFML_String_Utilities.f90                 /nologo %OPT1% %OPT2%
   ifort /c CFML_Rational_Arithmetic.f90              /nologo %OPT1% %OPT2%
   if [%_WINTER%]==[Y] (
     ifort /c CFML_IO_MessagesWin.f90                 /nologo %OPT1% %OPT2% %OPT3%
   ) else (
     ifort /c CFML_IO_Messages.f90                    /nologo %OPT1% %OPT2%
   )
   ifort /c CFML_PowderProfiles_TOF.f90               /nologo %OPT1% %OPT2%
   ifort /c CFML_PowderProfiles_Finger.f90            /nologo %OPT1% %OPT2%
   ifort /c CFML_PowderProfiles_CW.f90                /nologo %OPT1% %OPT2%
rem
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   ifort /c CFML_Math_3D.f90                          /nologo %OPT1% %OPT2%
   ifort /c CFML_Optimization_General.f90             /nologo %OPT1% %OPT2%
   ifort /c CFML_Optimization_LSQ.f90                 /nologo %OPT1% %OPT2%
   ifort /c CFML_Symmetry_Tables.f90                  /nologo %OPT0% %OPT2%
   ifort /c CFML_Scattering_Chemical_Tables.f90       /nologo %OPT0% %OPT2%
   ifort /c CFML_BVSpar.f90                           /nologo %OPT0% %OPT2%
   ifort /c CFML_Diffraction_Patterns.f90             /nologo %OPT1% %OPT2%
rem
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   ifort /c CFML_Bond_Tables.f90                      /nologo %OPT0% %OPT2%
   ifort /c CFML_Crystal_Metrics.f90                  /nologo %OPT1% %OPT2%
   ifort /c CFML_Crystallographic_Symmetry.f90        /nologo %OPT1% %OPT2%
   ifort /c CFML_ILL_Instrm_Data.f90                  /nologo %OPT1% %OPT2%
   ifort /c CFML_Magnetic_Groups.f90                  /nologo %OPT1% %OPT2%
rem
   echo .... EoS, Reflections, Atoms
rem
   ifort /c CFML_Eos.f90                              /nologo %OPT1% %OPT2%
   ifort /c CFML_Reflections_Utilities.f90            /nologo %OPT1% %OPT2%
   ifort /c CFML_Atom_TypeDef.f90                     /nologo %OPT1% %OPT2%
rem
   echo .... Formats, Geometry, Molecules
rem
   ifort /c CFML_Geometry_Calc.f90                    /nologo %OPT1% %OPT2%
   ifort /c CFML_Molecular_Crystals.f90               /nologo %OPT1% %OPT2%
   ifort /c CFML_IO_Formats.f90                       /nologo %OPT1% %OPT2%
rem
   echo .... Extinction, Structure Factors, SXTAL geometry, Propag Vectors
rem
   ifort /c CFML_Extinction_Correction.f90           /nologo %OPT1% %OPT2%
   ifort /c CFML_Structure_Factors.f90               /nologo %OPT1% %OPT2%
   ifort /c CFML_Geometry_SXTAL.f90                  /nologo %OPT1% %OPT2%
   ifort /c CFML_Propagation_Vectors.f90             /nologo %OPT1% %OPT2%
rem
   echo .... Maps, BVS, Energy Configurations
rem
   ifort /c CFML_Export_Vtk.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_Maps_Calculations.f90               /nologo %OPT1% %OPT2%
   ifort /c CFML_BVS_Energy_Calc.f90                 /nologo %OPT1% %OPT2%
   ifort /c CFML_Percolation.f90                     /nologo %OPT1% %OPT2%
rem
   echo .... Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem
rem   ifort /c CFML_Magnetic_Groups.f90                 /nologo %OPT1% %OPT2%
   ifort /c CFML_Magnetic_Symmetry.f90               /nologo %OPT1% %OPT2%
   ifort /c CFML_Simulated_Annealing.f90             /nologo %OPT1% %OPT2% %OPT3%
   ifort /c CFML_Keywords_Code_Parser.f90            /nologo %OPT1% %OPT2%
rem
   echo .... Magnetic Structure Factors, Polarimetry
rem
   ifort /c CFML_Magnetic_Structure_Factors.f90      /nologo %OPT1% %OPT2%
   ifort /c CFML_Polarimetry.f90                     /nologo %OPT1% %OPT2%
rem
rem
   echo.
   echo Creating CrysFML Library
rem
   if [%_WINTER%]==[Y] (
     lib /out:wcrysfml.lib *.obj
   ) else (
     lib /out:crysfml.lib *.obj
   )
rem
   echo.
   echo Creating IFORT directory
rem
   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
rem
   if [%_WINTER%]==[Y] (
      if exist ..\%DIRECTORY%\LibW\ rmdir ..\%DIRECTORY%\LibW\ /S /Q
      mkdir ..\%DIRECTORY%\LibW\
      copy *.mod ..\%DIRECTORY%\LibW\. > nul
      move *.lib ..\%DIRECTORY%\LibW\. > nul
   ) else (
      if exist ..\%DIRECTORY%\LibC\ rmdir ..\%DIRECTORY%\LibC\ /S /Q
      mkdir ..\%DIRECTORY%\LibC\
      copy *.mod ..\%DIRECTORY%\LibC\. > nul
      move *.lib ..\%DIRECTORY%\LibC\. > nul
   )
   del *.obj *.mod *.lst *.bak > nul
rem
   echo.
   echo **---- End Compilation for CrysFML ----**
   echo.
   cd %CRYSFML%\Scripts\Windows
rem
