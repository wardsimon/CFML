rem
rem CrysFML for G95 Compiler (Debug)
rem
   @echo off
   cd ../../cfml
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   g95 -c CFML_math_gen.f90         -O0 -ftrace=full  -std=f2003    >  out
   g95 -c CFML_string_util.f90      -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_io_mess.f90          -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_random.f90           -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_ffts.f90             -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_Profile_TOF.f90      -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_Profile_Finger.f90   -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_Profile_Functs.f90   -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   g95 -c CFML_math_3D.f90          -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_spher_harm.f90       -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_optimization.f90     -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_optimization_lsq.f90 -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_sym_table.f90        -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_chem_scatt.f90       -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_diffpatt.f90         -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   g95 -c CFML_bonds_table.f90      -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_cryst_types.f90      -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_symmetry.f90         -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_ILL_Instrm_data.f90  -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms, Polarimetry
rem
   g95 -c CFML_reflct_util.f90      -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_atom_mod.f90         -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_polar.f90            -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_SXTAL_geom.f90       -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, Propag Vectors
rem
   g95 -c CFML_sfac.f90             -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_geom_calc.f90        -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_propagk.f90          -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   g95 -c CFML_maps.f90             -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_molecules.f90        -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_conf_calc.f90        -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   g95 -c CFML_form_cif.f90         -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Level 7 ----**
   echo .... Keywords Parser, Simulated Annealing, Magnetic Symmetry
rem   
   g95 -c CFML_refcodes.f90         -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_optimization_san.f90 -O0 -ftrace=full  -std=f2003    >> out
   g95 -c CFML_magsymm.f90          -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors
rem
   g95 -c CFML_msfac.f90            -O0 -ftrace=full  -std=f2003    >> out
rem
   echo **---- Crysfml Library: Console version ----**
rem
   ar cr libcrysfml.a *.o
rem
   echo **---- G95 Directory ----**
rem
   if not exist ..\..\G95 mkdir ..\..\G95
   if exist ..\..\G95\LibC rmdir ..\..\G95\LibC /S /Q
   mkdir ..\..\G95\LibC
rem
   copy *.mod ..\..\G95\LibC > nul
   move *.a   ..\..\G95\LibC > nul
   del *.o *.mod *.lst *.bak > nul
rem   
   cd ..\Scripts\Windows   
