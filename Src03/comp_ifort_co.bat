rem
rem CrysFML for Intel Compiler (Optimization)
rem
   @echo off
rem
rem Options Definitions
rem
   set opt1=/O3 /nologo
   set opt2=%opt1% /stand:f03
rem
rem
rem
rem
   cd %CRYSFML%\Src03
rem
   echo **-----------------**
   echo **---- Level 0 ----**
   echo **-----------------**
rem
rem   echo ....  Messages, Powder Profiles
rem
   ifort /c CFML_GlobalDeps_Windows_intel.f90         %opt1%
rem
   echo .... Mathematical(I)....
   ifort /c CFML_math_gen.f90                         %opt2%
   ifort /c CFML_random.f90                           %opt2%
   ifort /c CFML_spher_harm.f90                       %opt2%
   ifort /c CFML_ffts.f90                             %opt2%

   echo .... String_Utilities ....

rem   ifort /c CFML_LSQ_TypeDef.f90      /O2 /nologo /Qvec-report0
rem   ifort /c CFML_string_util.f90      /O2 /nologo /Qvec-report0
rem   ifort /c CFML_io_mess.f90          /O2 /nologo /Qvec-report0
rem   ifort /c CFML_Profile_TOF.f90      /O2 /nologo /Qvec-report0
rem   ifort /c CFML_Profile_Finger.f90   /O2 /nologo /Qvec-report0
rem   ifort /c CFML_Profile_Functs.f90   /O2 /nologo /Qvec-report0
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
rem   ifort /c CFML_math_3D.f90          /O2 /nologo /Qvec-report0
rem   ifort /c CFML_optimization.f90     /O2 /nologo /Qvec-report0
rem   ifort /c CFML_optimization_lsq.f90 /O2 /nologo /Qvec-report0
rem   ifort /c CFML_sym_table.f90        /Od /nologo /Qvec-report0
rem   ifort /c CFML_chem_scatt.f90       /Od /nologo /Qvec-report0
rem   ifort /c CFML_BVSpar.f90           /Od /nologo /Qvec-report0
rem   ifort /c CFML_diffpatt.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
rem   ifort /c CFML_bonds_table.f90      /Od /nologo /Qvec-report0
rem   ifort /c CFML_cryst_types.f90      /O2 /nologo /Qvec-report0
rem   ifort /c CFML_symmetry.f90         /O2 /nologo /Qvec-report0
rem   ifort /c CFML_ILL_Instrm_data.f90  /O2 /nologo /Qvec-report0
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms
rem
rem   ifort /c CFML_Eos_Mod.f90          /O2 /nologo /Qvec-report0
rem   ifort /c CFML_reflct_util.f90      /O2 /nologo /Qvec-report0
rem   ifort /c CFML_atom_mod.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, SXTAL geometry, Propag Vectors
rem
rem   ifort /c CFML_sfac.f90             /O2 /nologo /Qvec-report0
rem   ifort /c CFML_geom_calc.f90        /O2 /nologo /Qvec-report0
rem   ifort /c CFML_sxtal_Geom.f90       /O2 /nologo /Qvec-report0
rem   ifort /c CFML_propagk.f90          /O2 /nologo /Qvec-report0
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
rem   ifort /c CFML_Export_Vtk.f90       /O2 /nologo /Qvec-report0
rem   ifort /c CFML_maps.f90             /O2 /nologo /Qvec-report0
rem   ifort /c CFML_molecules.f90        /O2 /nologo /Qvec-report0
rem   ifort /c CFML_conf_calc.f90        /O2 /nologo /Qvec-report0
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
rem   ifort /c CFML_form_cif.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 7 ----**
   echo .... Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem
rem   ifort /c CFML_magsymm.f90          /O2 /nologo /Qvec-report0
rem   ifort /c CFML_optimization_san.f90 /O2 /nologo /Qvec-report0
rem   ifort /c CFML_refcodes.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
rem   ifort /c CFML_msfac.f90            /O2 /nologo /Qvec-report0
rem   ifort /c CFML_polar.f90            /O2 /nologo /Qvec-report0
rem
rem
   echo **---- Crysfml Library: Console version ----**
rem
rem   lib /out:crysfml.lib *.obj
rem
   echo **---- ifort Directory ----**
rem
rem   if "%TARGET_ARCH%"==""   set TARGET_ARCH=ia32
rem   if %TARGET_ARCH% == ia32 (set DIRECTORY=ifort) else (set DIRECTORY=ifort64)
rem
rem   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
rem   if exist ..\%DIRECTORY%\LibC rmdir ..\%DIRECTORY%\LibC /S /Q
rem   mkdir ..\%DIRECTORY%\LibC
rem
rem   copy *.mod ..\%DIRECTORY%\LibC > nul
rem   move *.lib ..\%DIRECTORY%\LibC > nul
   del *.mod *.lst *.bak > nul
rem
rem   cd %CRYSFML%\Scripts\Windows
