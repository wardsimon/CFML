rem
rem CrysFML for Lahey Compiler (Debug)
rem
   @echo off
   cd ../../cfml
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   lf95 -c f2kcli.f90                -g  -chk   >  out
   lf95 -c CFML_GlobalDeps.f90         -g  -chk   >> out
rem 
   lf95 -c CFML_math_gen.f90         -g  -chk   >> out
   lf95 -c CFML_spher_harm.f90       -g  -chk   >> out
   lf95 -c CFML_random.f90           -g  -chk   >> out
   lf95 -c CFML_ffts.f90             -g  -chk   >> out
   lf95 -c CFML_string_util.f90      -g  -chk   >> out
   lf95 -c CFML_io_mess.f90          -g  -chk   >> out
   lf95 -c CFML_Profile_TOF.f90      -g  -chk   >> out
   lf95 -c CFML_Profile_Finger.f90   -g  -chk   >> out
   lf95 -c CFML_Profile_Functs.f90   -g  -chk   >> out
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   lf95 -c CFML_math_3D.f90          -g  -chk   >> out
   lf95 -c CFML_optimization.f90     -g  -chk   >> out
   lf95 -c CFML_optimization_lsq.f90 -g  -chk   >> out
   lf95 -c CFML_sym_table.f90        -g  -chk   >> out
   lf95 -c CFML_chem_scatt.f90       -g  -chk   >> out
   lf95 -c CFML_diffpatt.f90         -g  -chk   >> out
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   lf95 -c CFML_bonds_table.f90      -g  -chk   >> out
   lf95 -c CFML_cryst_types.f90      -g  -chk   >> out
   lf95 -c CFML_symmetry.f90         -g  -chk   >> out
   lf95 -c CFML_ILL_Instrm_data.f90  -g  -chk   >> out
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms, Polarimetry
rem
   lf95 -c CFML_reflct_util.f90      -g  -chk   >> out
   lf95 -c CFML_atom_mod.f90         -g  -chk   >> out
   lf95 -c CFML_polar.f90            -g -chk   >> out
   lf95 -c CFML_SXTAL_geom.f90       -g  -chk   >> out
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, Propag Vectors
rem
   lf95 -c CFML_sfac.f90             -g  -chk   >> out
   lf95 -c CFML_geom_calc.f90        -g  -chk   >> out
   lf95 -c CFML_propagk.f90          -g  -chk   >> out
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   lf95 -c CFML_maps.f90             -g  -chk   >> out
   lf95 -c CFML_molecules.f90        -g  -chk   >> out
   lf95 -c CFML_conf_calc.f90        -g  -chk   >> out
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   lf95 -c CFML_form_cif.f90         -g  -chk   >> out
rem
   echo **---- Level 7 ----**
   echo .... Keywords Parser, Simulated Annealing, Magnetic Symmetry
rem   
   lf95 -c CFML_refcodes.f90         -g  -chk   >> out
   lf95 -c CFML_optimization_san.f90 -g  -chk   >> out
   lf95 -c CFML_magsymm.f90          -g -chk    >> out
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors
rem
   lf95 -c CFML_msfac.f90            -g -chk    >> out
rem
rem
   echo **---- Crysfml Library: Console version (DEBUG) ----**
rem
   lm @..\Scripts\Windows\lib_mod.lnk
rem
   echo **---- Lahey Directory ----**
rem
   if not exist ..\..\Lahey mkdir ..\..\Lahey
   if exist ..\..\Lahey\LibC rmdir ..\..\Lahey\LibC /S /Q
   mkdir ..\..\Lahey\LibC
rem
   copy *.mod ..\..\Lahey\LibC > nul
   move *.lib ..\..\Lahey\LibC > nul
   del *.obj *.mod *.lst *.bak > nul
rem         
   cd ..\Scripts\Windows