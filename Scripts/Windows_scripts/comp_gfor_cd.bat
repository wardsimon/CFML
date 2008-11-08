rem
rem CrysFML for GFortran Compiler (Debug)
rem
   @echo off
   cd ../../cfml
rem
   echo **---- Level 0 ----**
   echo .... Mathematical,String_Utilities, Profile Functions,
   echo .... Graphical, Optimization Modules
rem
   gfortran -c CFML_math_gen.f90         -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_string_util.f90      -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_random.f90           -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_fft.f90              -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_fft_nF.f90           -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_Profile_TOF.f90      -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_Profile_Finger.f90   -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_Profile_Functs.f90   -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_optimization.f90     -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_optimization_lsq.f90 -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 1 ----**
   echo .... Math_3D, Spher_Harm
rem
   gfortran -c CFML_math_3D.f90          -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_spher_harm.f90       -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 2 ----**
   echo .... Sym_Table, Chem_Scatt, Cryst_Types, DiffPatt Modules
rem
   gfortran -c CFML_sym_table.f90        -O0                  -std=f2003
   gfortran -c CFML_chem_scatt.f90       -O0                  -std=f2003
   gfortran -c CFML_cryst_types.f90      -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_diffpatt.f90         -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 3 ----**
   echo .... Bonds_Table, Symmetry Modules
rem
   gfortran -c CFML_bonds_table.f90      -O0                  -std=f2003
   gfortran -c CFML_symmetry.f90         -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 4 ----**
   echo .... Reflct_Util, Atom_Mod Modules
rem
   gfortran -c CFML_reflct_util.f90      -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_atom_mod.f90         -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 5 ----**
   echo .... Sfac, Propagk, Geom_Calc, Maps, Molecules Modules
rem
   gfortran -c CFML_sfac.f90             -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_propagk.f90          -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_geom_calc.f90        -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_maps.f90             -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_molecules.f90        -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 6 ----**
   echo .... Form_CIF, Conf_Calc, RefCodes Modules
rem
   gfortran -c CFML_form_cif.f90         -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_conf_calc.f90        -O0                  -std=f2003
   gfortran -c CFML_refcodes.f90         -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 7 ----**
   echo .... Polar, Symm, SFac for Magnetic  Modules
rem   
   gfortran -c CFML_polar.f90            -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_magsymm.f90          -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_msfac.f90            -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 8 ----**
   echo .... ILL_Instrm_data, CFML_SXTAL_geom Modules
rem
   gfortran -c CFML_ILL_Instrm_data.f90  -O0 -ftrace=full     -std=f2003
   gfortran -c CFML_SXTAL_geom.f90       -O0 -ftrace=full     -std=f2003
rem
   echo **---- Level 9 ----**
   echo .... IO_Mess Module
rem
   gfortran -c CFML_io_mess.f90          -O0 -ftrace=full  -std=f2003
rem
   echo **---- Level 10 ----**
   echo .... Optimization Simulated Annealing Module
rem
   gfortran -c CFML_optimization_san.f90 -O0 -ftrace=full  -std=f2003
rem
   echo **---- Crysfml Library: Console (DEBUG) version ----**
rem
   ar cr libcrysfml.a *.o
rem
   echo **---- GFortran Directory ----**
rem
   if not exist ..\..\GFortran mkdir ..\..\GFortran
   if exist ..\..\GFortran\LibC rmdir ..\..\GFortran\LibC /S /Q
   mkdir ..\..\GFortran\LibC
rem
   copy *.mod ..\..\GFortran\LibC > nul
   move *.a   ..\..\GFortran\LibC > nul
   del *.o *.mod *.lst *.bak > nul
rem 
   cd ..\Scripts\Windows