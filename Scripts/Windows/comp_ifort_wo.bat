rem
rem CrysFML for Intel Compiler (Optimization) + WINTERACTER
rem
   @echo off
   cd %CRYSFML%\Src
rem
   if %TARGET_ARCH% == ia32 (set LIBFOR=lib.if8) else (set LIBFOR=lib.i64)
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Profile Functions
rem
   ifort /c CFML_GlobalDeps_Windows_intel.f90         /O2 /nologo /Qopt-report:0
rem
   ifort /c CFML_math_gen.f90         /O2 /nologo /Qopt-report:0
   ifort /c CFML_LSQ_TypeDef.f90      /O2 /nologo /Qopt-report:0
   ifort /c CFML_spher_harm.f90       /O2 /nologo /Qopt-report:0
   ifort /c CFML_random.f90           /O2 /nologo /Qopt-report:0
   ifort /c CFML_ffts.f90             /O2 /nologo /Qopt-report:0
   ifort /c CFML_string_util.f90      /O2 /nologo /Qopt-report:0
   ifort /c CFML_io_messwin.f90       /O2 /nologo /Qopt-report:0 /I%WINTER%\%LIBFOR%
   ifort /c CFML_Profile_TOF.f90      /O2 /nologo /Qopt-report:0
   ifort /c CFML_Profile_Finger.f90   /O2 /nologo /Qopt-report:0
   ifort /c CFML_Profile_Functs.f90   /O2 /nologo /Qopt-report:0
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   ifort /c CFML_math_3D.f90          /O2 /nologo /Qopt-report:0
   ifort /c CFML_optimization.f90     /O2 /nologo /Qopt-report:0
   ifort /c CFML_optimization_lsq.f90 /O2 /nologo /Qopt-report:0
   ifort /c CFML_sym_table.f90        /Od /nologo /Qopt-report:0
   ifort /c CFML_chem_scatt.f90       /Od /nologo /Qopt-report:0
   ifort /c CFML_BVSpar.f90           /Od /nologo /Qopt-report:0
   ifort /c CFML_diffpatt.f90         /O2 /nologo /Qopt-report:0
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   ifort /c CFML_bonds_table.f90      /Od /nologo /Qopt-report:0
   ifort /c CFML_cryst_types.f90      /O2 /nologo /Qopt-report:0
   ifort /c CFML_symmetry.f90         /O2 /nologo /Qopt-report:0
   ifort /c CFML_ILL_Instrm_data.f90  /O2 /nologo /Qopt-report:0
rem
   echo **---- Level 3 ----**
   echo .... EoS, Reflections, Atoms
rem
   ifort /c CFML_Eos_Mod.f90          /O2 /nologo /Qopt-report:0
   ifort /c CFML_reflct_util.f90      /O2 /nologo /Qopt-report:0
   ifort /c CFML_atom_mod.f90         /O2 /nologo /Qopt-report:0
rem
   echo **---- Level 4 ----**
   echo .... Formats, Geometry Calculations, Molecules
rem
   ifort /c CFML_geom_calc.f90        /O2 /nologo /Qopt-report:0
   ifort /c CFML_molecules.f90        /O2 /nologo /Qopt-report:0
   ifort /c CFML_form_cif.f90         /O2 /nologo /Qopt-report:0
rem
   echo **---- Level 5 ----**
   echo .... Extinction, Structure Factors, SXTAL geometry, Propag Vectors
rem
   ifort /c CFML_Extinction_Correction.f90 /O2 /nologo /Qopt-report:0
   ifort /c CFML_sfac.f90                  /O2 /nologo /Qopt-report:0
   ifort /c CFML_SXTAL_Geom.f90            /O2 /nologo /Qopt-report:0
   ifort /c CFML_propagk.f90               /O2 /nologo /Qopt-report:0
rem
   echo **---- Level 6 ----**
   echo .... Maps, BVS, Energy Configurations
rem
   ifort /c CFML_Export_Vtk.f90       /O2 /nologo /Qopt-report:0
   ifort /c CFML_maps.f90             /O2 /nologo /Qopt-report:0
   ifort /c CFML_conf_calc.f90        /O2 /nologo /Qopt-report:0
rem
   echo **---- Level 7 ----**
   echo .... Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem

   ifort /c CFML_Magnetic_Groups.f90  /O2 /nologo /Qopt-report:0
   ifort /c CFML_magsymm.f90          /O2 /nologo /Qopt-report:0
   ifort /c CFML_optimization_san.f90 /O2 /nologo /Qopt-report:0 /I%WINTER%\%LIBFOR%
   ifort /c CFML_refcodes.f90         /O2 /nologo /Qopt-report:0
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
   ifort /c CFML_msfac.f90            /O2 /nologo /Qopt-report:0
   ifort /c CFML_polar.f90            /O2 /nologo /Qopt-report:0
rem
rem
   echo **---- Crysfml Library: Winteracter version ----**
rem
   lib /out:wcrysfml.lib *.obj
rem
   echo **---- ifort Directory ----**
rem
   if "%TARGET_ARCH%"==""   set TARGET_ARCH=ia32
   if %TARGET_ARCH% == ia32 (set DIRECTORY=ifort) else (set DIRECTORY=ifort64)
rem
   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
   if exist ..\%DIRECTORY%\LibW rmdir ..\%DIRECTORY%\LibW /S /Q
   mkdir ..\%DIRECTORY%\LibW
rem
   copy *.mod ..\%DIRECTORY%\LibW > nul
   move *.lib ..\%DIRECTORY%\LibW > nul
   del *.obj *.mod *.lst *.bak > nul
rem
   cd %CRYSFML%\Scripts\Windows
