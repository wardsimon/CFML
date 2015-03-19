rem
rem CrysFML for GFortran Compiler (Optimization)
rem
   @echo off
   cd %CRYSFML%\Src03
rem
   echo **-----------------**
   echo **---- Level 0 ----**
   echo **-----------------**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   gfortran -c CFML_GlobalDeps_Windows.f90         -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >  out
rem
   echo .... Mathematical(I)....
   gfortran -c CFML_math_gen.f90                   -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
   gfortran -c CFML_random.f90                     -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
   gfortran -c CFML_spher_harm.f90                 -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
   gfortran -c CFML_ffts.f90                       -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_LSQ_TypeDef.f90      -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_spher_harm.f90       -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_string_util.f90   -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_io_mess.f90          -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem -std=f2003 removed in TOF because probably there is a conflict of names
rem  is erfc an intrinsic function in  F2003?
rem   gfortran -c CFML_Profile_TOF.f90      -O3 -ffree-line-length-none             -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_Profile_Finger.f90   -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_Profile_Functs.f90   -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
rem   gfortran -c CFML_math_3D.f90          -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_optimization.f90     -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_optimization_lsq.f90 -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_sym_table.f90        -O0 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_chem_scatt.f90       -O0 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_BVSpar.f90           -O0 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_diffpatt.f90         -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
rem   gfortran -c CFML_bonds_table.f90      -O0 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_cryst_types.f90      -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_symmetry.f90         -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_ILL_Instrm_data.f90  -O3 -ffree-line-length-none -std=gnu    -funroll-loops  -msse2  >> out
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms
rem
rem   gfortran -c CFML_Eos_Mod.f90          -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_reflct_util.f90      -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_atom_mod.f90         -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, SXTAL geometry, Propag Vectors
rem
rem   gfortran -c CFML_sfac.f90             -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_geom_calc.f90        -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_SXTAL_geom.f90       -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_propagk.f90          -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
rem   gfortran -c CFML_Export_Vtk.f90       -O3 -ffree-line-length-none -std=gnu    -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_maps.f90             -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_molecules.f90        -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem -std=f2003 removed because calls to flush subroutine
rem   gfortran -c CFML_conf_calc.f90        -O3 -ffree-line-length-none             -funroll-loops  -msse2  >> out
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
rem   gfortran -c CFML_form_cif.f90         -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem
   echo **---- Level 7 ----**
   echo .... Keywords Parser, Simulated Annealing, Magnetic Symmetry
rem
rem -std=f2003 removed because calls to flush subroutine
rem   gfortran -c CFML_optimization_san.f90 -O3 -ffree-line-length-none             -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_magsymm.f90          -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_refcodes.f90         -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
rem   gfortran -c CFML_msfac.f90            -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem   gfortran -c CFML_polar.f90            -O3 -ffree-line-length-none -std=f2003  -funroll-loops  -msse2  >> out
rem
   echo **---- Crysfml Library: Console version ----**
rem
   ar cr libcrysfml.a *.o
rem
   echo **---- GFortran Directory ----**
rem
rem   if not exist ..\GFortran  mkdir ..\GFortran
rem   if exist ..\GFortran\LibC rmdir ..\GFortran\LibC /S /Q
rem   mkdir ..\GFortran\LibC
rem
rem   copy *.mod ..\GFortran\LibC > nul
rem   move *.a   ..\GFortran\LibC > nul
rem   del *.o *.mod *.lst *.bak > nul
rem
rem   cd %CRYSFML%\Scripts\Windows
