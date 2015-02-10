rem
rem CrysFML for GFortran Compiler (Debug)
rem
   @echo on
   cd %CRYSFML%\Src
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   gfortran -c CFML_GlobalDeps_Windows.f90         -O0 -fbacktrace  -std=f2003 -ffree-line-length-none
rem
   gfortran -c CFML_math_gen.f90         -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_LSQ_TypeDef.f90      -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_spher_harm.f90       -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_random.f90           -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_ffts.f90             -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_string_util.f90      -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_io_mess.f90          -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
rem -std=f2003 removed in TOF because probably there is a conflict of names
rem  is erfc an intrinsic function in  F2003?
   gfortran -c CFML_Profile_TOF.f90      -O0 -fbacktrace  -std=gnu     -ffree-line-length-none 
   gfortran -c CFML_Profile_Finger.f90   -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_Profile_Functs.f90   -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   gfortran -c CFML_math_3D.f90          -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_optimization.f90     -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_optimization_lsq.f90 -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_sym_table.f90        -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_chem_scatt.f90       -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_BVSpar.f90           -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_diffpatt.f90         -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   gfortran -c CFML_bonds_table.f90      -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
   gfortran -c CFML_cryst_types.f90      -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
   gfortran -c CFML_symmetry.f90         -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
   gfortran -c CFML_ILL_Instrm_data.f90  -O0 -fbacktrace  -std=gnu    -ffree-line-length-none
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms
rem
   gfortran -c CFML_Eos_Mod.f90          -O0 -fbacktrace  -std=f2003  -ffree-line-length-none 
   gfortran -c CFML_reflct_util.f90      -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
   gfortran -c CFML_atom_mod.f90         -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, SXTAL geometry, Propag Vectors
rem
   gfortran -c CFML_sfac.f90             -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_geom_calc.f90        -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_SXTAL_geom.f90       -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
   gfortran -c CFML_propagk.f90          -O0 -fbacktrace  -std=f2003   -ffree-line-length-none
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   gfortran -c CFML_Export_Vtk.f90       -O0 -fbacktrace  -std=gnu    -ffree-line-length-none
   gfortran -c CFML_maps.f90             -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
   gfortran -c CFML_molecules.f90        -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
rem -std=f2003 removed because calls to flush subroutine
   gfortran -c CFML_conf_calc.f90        -O0 -fbacktrace              -ffree-line-length-none
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   gfortran -c CFML_form_cif.f90         -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
rem -std=f2003 removed because calls to sizeof function
rem
   echo **---- Level 7 ----**
   echo .... Keywords Parser, Simulated Annealing, Magnetic Symmetry
rem
rem -std=f2003 removed because calls to flush subroutine
   gfortran -c CFML_optimization_san.f90 -O0 -fbacktrace               -ffree-line-length-none
   gfortran -c CFML_magsymm.f90          -O0 -fbacktrace -std=f2003    -ffree-line-length-none
   gfortran -c CFML_refcodes.f90         -O0 -fbacktrace -std=f2003    -ffree-line-length-none
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
   gfortran -c CFML_msfac.f90            -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
   gfortran -c CFML_polar.f90            -O0 -fbacktrace  -std=f2003  -ffree-line-length-none
rem
   echo **---- Crysfml Library: Console (DEBUG) version ----**
rem
   ar cr libcrysfml.a *.o
rem
   echo **---- GFortran Directory ----**
rem
   if not exist ..\GFortran_debug mkdir ..\GFortran_debug
   if exist ..\GFortran_debug\LibC rmdir ..\GFortran_debug\LibC /S /Q
   mkdir ..\GFortran_debug\LibC
rem
   copy *.mod ..\GFortran_debug\LibC > nul
   move *.a   ..\GFortran_debug\LibC > nul
   del *.o *.mod *.lst *.bak > nul
rem
   cd %CRYSFML%\Scripts\Windows
