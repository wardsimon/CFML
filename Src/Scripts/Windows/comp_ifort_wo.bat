rem
rem CrysFML for Intel Compiler (Optimization) + WINTERACTER
rem
   @echo off
   cd ../../cfml
rem
   echo **---- Level 0 ----**
   echo .... Mathematical,String_Utilities, Messages, Profile Functions
rem
   ifort /c CFML_Constant.f90         /O2 /nologo /Qvec-report0
rem
   ifort /c CFML_math_gen.f90         /O2 /nologo /Qvec-report0
   ifort /c CFML_spher_harm.f90       /O2 /nologo /Qvec-report0
   ifort /c CFML_random.f90           /O2 /nologo /Qvec-report0
   ifort /c CFML_ffts.f90             /O2 /nologo /Qvec-report0
   ifort /c CFML_string_util.f90      /O2 /nologo /Qvec-report0
   ifort /c CFML_io_messwin.f90       /O2 /nologo /Qvec-report0 /Ic:\wint\lib.if8
   ifort /c CFML_Profile_TOF.f90      /O2 /nologo /Qvec-report0
   ifort /c CFML_Profile_Finger.f90   /O2 /nologo /Qvec-report0
   ifort /c CFML_Profile_Functs.f90   /O2 /nologo /Qvec-report0
rem
   echo **---- Level 1 ----**
   echo .... Mathematical, Optimization, Tables, Patterns
rem
   ifort /c CFML_math_3D.f90          /O2 /nologo /Qvec-report0
   ifort /c CFML_optimization.f90     /O2 /nologo /Qvec-report0
   ifort /c CFML_optimization_lsq.f90 /O2 /nologo /Qvec-report0
   ifort /c CFML_sym_table.f90        /Od /nologo /Qvec-report0
   ifort /c CFML_chem_scatt.f90       /Od /nologo /Qvec-report0
   ifort /c CFML_diffpatt.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 2 ----**
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   ifort /c CFML_bonds_table.f90      /Od /nologo /Qvec-report0
   ifort /c CFML_cryst_types.f90      /O2 /nologo /Qvec-report0
   ifort /c CFML_symmetry.f90         /O2 /nologo /Qvec-report0
   ifort /c CFML_ILL_Instrm_data.f90  /O2 /nologo /Qvec-report0
rem
   echo **---- Level 3 ----**
   echo .... Reflections, Atoms, Polarimetry
rem
   ifort /c CFML_reflct_util.f90      /O2 /nologo /Qvec-report0
   ifort /c CFML_atom_mod.f90         /O2 /nologo /Qvec-report0
   ifort /c CFML_polar.f90            /O2 /nologo /Qvec-report0
   ifort /c CFML_SXTAL_Geom.f90       /O2 /nologo /Qvec-report0
rem
   echo **---- Level 4 ----**
   echo .... Structure Factors, Geometry Calculations, Propag Vectors
rem
   ifort /c CFML_sfac.f90             /O2 /nologo /Qvec-report0
   ifort /c CFML_geom_calc.f90        /O2 /nologo /Qvec-report0
   ifort /c CFML_propagk.f90          /O2 /nologo /Qvec-report0
rem
   echo **---- Level 5 ----**
   echo .... Molecules, Maps, BVS, Energy Configurations
rem
   ifort /c CFML_maps.f90             /O2 /nologo /Qvec-report0
   ifort /c CFML_molecules.f90        /O2 /nologo /Qvec-report0
   ifort /c CFML_conf_calc.f90        /O2 /nologo /Qvec-report0
rem
   echo **---- Level 6 ----**
   echo .... Formats
rem
   ifort /c CFML_form_cif.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 7 ----**
   echo .... Keywords Parser, Simulated Annealing, Magnetic Symmetry
rem   
   ifort /c CFML_refcodes.f90         /O2 /nologo /Qvec-report0
   ifort /c CFML_optimization_san.f90 /O2 /nologo /Qvec-report0 /Ic:\wint\lib.if8
   ifort /c CFML_magsymm.f90          /O2 /nologo /Qvec-report0
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors
rem
   ifort /c CFML_msfac.f90            /O2 /nologo /Qvec-report0
rem
rem
   echo **---- Crysfml Library: Winteracter version ----**
rem
   lib /out:wcrysfml.lib *.obj
rem
   echo **---- Intel Directory ----**
rem
   if not exist ..\..\Intel mkdir ..\..\Intel
   if exist ..\..\Intel\LibW rmdir ..\..\Intel\LibW /S /Q
   mkdir ..\..\Intel\LibW
rem
   copy *.mod ..\..\Intel\LibW > nul
   move *.lib ..\..\Intel\LibW > nul
   del *.obj *.mod *.lst *.bak > nul
rem  
   cd ..\Scripts\Windows  