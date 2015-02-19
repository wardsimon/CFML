rem
rem CrysFML for Lahey Compiler (Debug)
rem
   @echo off
   cd %CRYSFML%\Src
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   lf95 -c f2kcli.f90                -info -g  -chk   >  out
   lf95 -c CFML_GlobalDeps_Windows.f90  -info -g  -chk   >> out
rem
   lf95 -c CFML_math_gen.f90         -info -g  -chk   >> out
   lf95 -c CFML_LSQ_TypeDef.f90      -info -g  -chk   >> out
   lf95 -c CFML_spher_harm.f90       -info -g  -chk   >> out
   lf95 -c CFML_random.f90           -info -g  -chk   >> out
   lf95 -c CFML_ffts.f90             -info -g  -chk   >> out
   lf95 -c CFML_string_util_LF.f90   -info -g  -chk   >> out
   lf95 -c CFML_io_mess.f90          -info -g  -chk   >> out
   lf95 -c CFML_Profile_TOF.f90      -info -g  -chk   >> out
   lf95 -c CFML_Profile_Finger.f90   -info -g  -chk   >> out
   lf95 -c CFML_Profile_Functs.f90   -info -g  -chk   >> out
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   lf95 -c CFML_math_3D.f90          -info -g  -chk   >> out
   lf95 -c CFML_optimization.f90     -info -g  -chk   >> out
   lf95 -c CFML_optimization_lsq.f90 -info -g  -chk   >> out
   lf95 -c CFML_sym_table.f90        -info -g  -chk   >> out
   lf95 -c CFML_chem_scatt.f90       -info -g  -chk   >> out
   lf95 -c CFML_BVSpar.f90           -info -g  -chk   >> out
   lf95 -c CFML_diffpatt.f90         -info -g  -chk   >> out
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   lf95 -c CFML_bonds_table.f90      -info -g  -chk   >> out
   lf95 -c CFML_cryst_types.f90      -info -g  -chk   >> out
   lf95 -c CFML_symmetry.f90         -info -g  -chk   >> out
   lf95 -c CFML_ILL_Instrm_data_LF.f90  -info -g  -chk   >> out
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms
rem
rem   lf95 -c CFML_Eos_Mod.f90          -info -g  -chk   >> out
   lf95 -c CFML_reflct_util.f90      -info -g  -chk   >> out
   lf95 -c CFML_atom_mod.f90         -info -g  -chk   >> out
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, SXTAL geometry, Propag Vectors
rem
   lf95 -c CFML_sfac.f90             -info -g  -chk   >> out
   lf95 -c CFML_geom_calc.f90        -info -g  -chk   >> out
   lf95 -c CFML_SXTAL_geom.f90       -info -g  -chk   >> out
   lf95 -c CFML_propagk.f90          -info -g  -chk   >> out
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   lf95 -c CFML_Export_Vtk_LF95.f90  -info -g  -chk   >> out
   lf95 -c CFML_maps.f90             -info -g  -chk   >> out
   lf95 -c CFML_molecules.f90        -info -g  -chk   >> out
   lf95 -c CFML_conf_calc.f90        -info -g  -chk   >> out
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   lf95 -c CFML_form_cif.f90         -info -g  -chk   >> out
rem
   echo **---- Level 7 ----**
   echo .... Keywords Parser, Simulated Annealing, Magnetic Symmetry
rem
   lf95 -c CFML_optimization_san.f90 -info -g  -chk   >> out
   lf95 -c CFML_magsymm.f90          -info -g  -chk   >> out
   lf95 -c CFML_refcodes.f90         -info -g  -chk   >> out
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
   lf95 -c CFML_msfac.f90            -info -g -chk    >> out
   lf95 -c CFML_polar.f90            -info -g -chk   >> out
rem
rem
   echo **---- Crysfml Library: Console version (DEBUG) ----**
rem
   lm @..\Scripts\Windows\lib_mod.lnk
rem
   echo **---- Lahey Directory ----**
rem
   if not exist ..\Lahey_debug mkdir ..\Lahey_debug
   if exist ..\Lahey_debug\LibC rmdir ..\Lahey_debug\LibC /S /Q
   mkdir ..\Lahey_debug\LibC
rem
   copy *.mod ..\Lahey_debug\LibC > nul
   move *.lib ..\Lahey_debug\LibC > nul
   del *.obj *.mod *.lst *.bak *.fwd > nul
rem
   cd %CRYSFML%\Scripts\Windows
