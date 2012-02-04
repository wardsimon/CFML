rem
rem CrysFML for Intel Compiler (Debug)
rem
   @echo off
   cd %CRYSFML%\Src
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   ifort /c CFML_GlobalDeps_Windows_intel.f90         /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   ifort /c CFML_math_gen.f90         /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_LSQ_TypeDef.f90      /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_spher_harm.f90       /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_random.f90           /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_ffts.f90             /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_string_util.f90      /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_io_mess.f90          /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_Profile_TOF.f90      /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_Profile_Finger.f90   /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_Profile_Functs.f90   /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   ifort /c CFML_math_3D.f90          /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_optimization.f90     /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_optimization_lsq.f90 /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_sym_table.f90        /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_chem_scatt.f90       /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_diffpatt.f90         /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   ifort /c CFML_bonds_table.f90      /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_cryst_types.f90      /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_symmetry.f90         /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_ILL_Instrm_data.f90  /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms, SXTAL geometry
rem
   ifort /c CFML_reflct_util.f90      /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_atom_mod.f90         /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_SXTAL_Geom.f90       /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, Propag Vectors
rem
   ifort /c CFML_sfac.f90             /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_geom_calc.f90        /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_propagk.f90          /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   ifort /c CFML_maps.f90             /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_molecules.f90        /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_conf_calc.f90        /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   ifort /c CFML_form_cif.f90         /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_Export_Vtk.f90       /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   echo **---- Level 7 ----**
   echo ....  Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem
   ifort /c CFML_magsymm.f90          /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_optimization_san.f90 /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_refcodes.f90         /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
   ifort /c CFML_msfac.f90            /debug:full /check /traceback /CB /nologo /Qvec-report0
   ifort /c CFML_polar.f90            /debug:full /check /traceback /CB /nologo /Qvec-report0
rem
rem
   echo **---- Crysfml Library: Console (DEBUG) version ----**
rem
   lib /out:crysfml.lib *.obj
rem
rem
   echo **---- Intel Directory ----**
rem
   if not exist ..\ifort_debug mkdir ..\ifort_debug
   if exist ..\ifort_debug\LibC rmdir ..\ifort_debug\LibC /S /Q
   mkdir ..\ifort_debug\LibC
rem
   copy *.mod ..\ifort_debug\LibC > nul
   move *.lib ..\ifort_debug\LibC > nul
   del *.obj *.mod *.lst *.bak > nul
rem
   cd %CRYSFML%\Scripts\Windows
