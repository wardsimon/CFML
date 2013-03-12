rem
rem CrysFML for Absoft Compiler (Optimization)
rem
   @echo off
   cd %CRYSFML%\Src
rem
   echo **---- Level 0 ----**
   echo .... Mathematical(I), String_Utilities, Messages, Powder Profiles
rem
   f95 -c -O3 -w CFML_GlobalDeps_Windows.f90           >  out
rem
   f95 -c -O3 -w CFML_math_gen.f90           >> out
   f95 -c -O3 -w CFML_LSQ_TypeDef.f90        >> out
   f95 -c -O3 -w CFML_spher_harm.f90         >> out
   f95 -c -O3 -w CFML_random.f90             >> out
   f95 -c -O3 -w CFML_ffts.f90               >> out
   f95 -c -O3 -w CFML_string_util.f90        >> out
   f95 -c -O3 -w CFML_io_mess.f90            >> out
   f95 -c -O3 -w CFML_Profile_TOF.f90        >> out
   f95 -c -O3 -w CFML_Profile_Finger.f90     >> out
   f95 -c -O3 -w CFML_Profile_Functs.f90     >> out
rem
   echo **---- Level 1 ----**
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   f95 -c -O3 -w CFML_math_3D.f90            >> out
   f95 -c -O3 -w CFML_optimization.f90       >> out
   f95 -c -O3 -w CFML_optimization_lsq.f90   >> out
   f95 -c -O0 -w CFML_sym_table.f90          >> out
   f95 -c -O0 -w CFML_chem_scatt.f90         >> out
   f95 -c -O3 -w CFML_diffpatt.f90           >> out
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   f95 -c -O0 -w CFML_bonds_table.f90        >> out
   f95 -c -O3 -w CFML_cryst_types.f90        >> out
   f95 -c -O3 -w CFML_symmetry.f90           >> out
   rem f95 -c -O3 -w CFML_ILL_Instrm_data.f90    >> out
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms
rem
   f95 -c -O3 -w CFML_reflct_util.f90        >> out
   f95 -c -O3 -w CFML_atom_mod.f90           >> out
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, SXTAL geometry, Propag Vectors
rem
   f95 -c -O3 -w CFML_sfac.f90               >> out
   f95 -c -O3 -w CFML_geom_calc.f90          >> out
   rem f95 -c -O3 -w CFML_SXTAL_geom.f90         >> out
   f95 -c -O3 -w CFML_propagk.f90            >> out
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   f95 -c -O3 -w CFML_maps.f90               >> out
   f95 -c -O3 -w CFML_molecules.f90          >> out
   f95 -c -O3 -w CFML_conf_calc.f90          >> out
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   f95 -c -O3 -w CFML_form_cif.f90           >> out
rem
   echo **---- Level 7 ----**
   echo .... Keywords Parser, Simulated Annealing, Magnetic Symmetry
rem
   f95 -c -O3 -w CFML_refcodes.f90           >> out
   f95 -c -O3 -w CFML_optimization_san.f90   >> out
   f95 -c -O3 -w CFML_magsymm.f90            >> out
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
   f95 -c -O3 -w CFML_msfac.f90              >> out
   f95 -c -O3 -w CFML_polar.f90              >> out
rem
rem
   echo **---- Crysfml Library: Console version ----**
rem
   lib -out:crysfml.lib *.obj
rem
   echo **---- Absoft Directory ----**
rem
   if not exist ..\Absoft mkdir ..\Absoft
   if exist ..\Absoft\LibC rmdir ..\Absoft\LibC /S /Q
   mkdir ..\Absoft\LibC
rem
   copy *.mod ..\Absoft\LibC > nul
   move *.lib   ..\Absoft\LibC > nul
   del *.obj *.mod *.lst *.bak > nul
rem
   cd %CRYSFML%\Scripts\Windows
