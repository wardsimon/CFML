rem
rem CrysFML for Intel Compiler (Optimization) + WINTERACTER
rem
   @echo off
   cd ../../cfml
rem
   echo **---- Level 0 ----**
   echo .... Mathematical,String_Utilities, Profile Functions,
   echo .... Graphical, Optimization Modules
rem
   ifort /c CFML_math_gen.f90         /O2 /nologo /Qvec-report0
   ifort /c CFML_string_util.f90      /O2 /nologo /Qvec-report0
   ifort /c CFML_random.f90           /O2 /nologo /Qvec-report0
   ifort /c CFML_fft.f90              /O2 /nologo /Qvec-report0
   ifort /c CFML_fft_nF.f90           /O2 /nologo /Qvec-report0
   ifort /c CFML_Profile_TOF.f90      /O2 /nologo /Qvec-report0
   ifort /c CFML_Profile_Finger.f90   /O2 /nologo /Qvec-report0
   ifort /c CFML_Profile_Functs.f90   /O2 /nologo /Qvec-report0
   ifort /c CFML_optimization.f90     /O2 /nologo /Qvec-report0
   ifort /c CFML_optimization_lsq.f90 /O2 /nologo /Qvec-report0
rem
   echo **---- Level 1 ----**
   echo .... Math_3D, Spher_Harm
rem
   ifort /c CFML_math_3D.f90          /O2 /nologo /Qvec-report0
   ifort /c CFML_spher_harm.f90       /O2 /nologo /Qvec-report0
rem
   echo **---- Level 2 ----**
   echo .... Sym_Table, Chem_Scatt, Cryst_Types, DiffPatt Modules
rem
   ifort /c CFML_sym_table.f90        /Od /nologo /Qvec-report0
   ifort /c CFML_chem_scatt.f90       /Od /nologo /Qvec-report0
   ifort /c CFML_cryst_types.f90      /O2 /nologo /Qvec-report0
   ifort /c CFML_diffpatt.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 3 ----**
   echo .... Bonds_Table, Symmetry Modules
rem
   ifort /c CFML_bonds_table.f90      /Od /nologo /Qvec-report0
   ifort /c CFML_symmetry.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 4 ----**
   echo .... Reflct_Util, Atom_Mod Modules
rem
   ifort /c CFML_reflct_util.f90      /O2 /nologo /Qvec-report0
   ifort /c CFML_atom_mod.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 5 ----**
   echo .... Sfac, Propagk, Geom_Calc, Maps, Molecules Modules
rem
   ifort /c CFML_sfac.f90             /O2 /nologo /Qvec-report0
   ifort /c CFML_propagk.f90          /O2 /nologo /Qvec-report0
   ifort /c CFML_geom_calc.f90        /O2 /nologo /Qvec-report0
   ifort /c CFML_maps.f90             /O2 /nologo /Qvec-report0
   ifort /c CFML_molecules.f90        /O2 /nologo /Qvec-report0
rem
   echo **---- Level 6 ----**
   echo .... Form_CIF, Conf_Calc, RefCodes Modules
rem
   ifort /c CFML_form_cif.f90         /O2 /nologo /Qvec-report0
   ifort /c CFML_conf_calc.f90        /O2 /nologo /Qvec-report0
   ifort /c CFML_refcodes.f90         /O2 /nologo /Qvec-report0
rem
   echo **---- Level 7 ----**
   echo .... Polar, Symm, SFac for Magnetic  Modules
rem   
   ifort /c CFML_polar.f90            /O2 /nologo /Qvec-report0
   ifort /c CFML_magsymm.f90          /O2 /nologo /Qvec-report0
   ifort /c CFML_msfac.f90            /O2 /nologo /Qvec-report0
rem
   echo **---- Level 8 ----**
   echo .... ILL_Instrm_data, CFML_SXTAL_geom Modules
rem
   ifort /c CFML_ILL_Instrm_data.f90  /O2 /nologo /Qvec-report0
   ifort /c CFML_SXTAL_geom.f90       /O2 /nologo /Qvec-report0
rem
   echo **---- Level 9 ----**
   echo .... IO_Mess Module
rem
   ifort /c CFML_io_messwin.f90      /O2 /nologo /Qvec-report0 -I. -Ic:\wint\lib.if8
rem
   echo **---- Level 10 ----**
   echo .... Optimization Simulated Annealing Module
rem
   ifort /c CFML_optimization_san.f90 /O2 /nologo /Qvec-report0 -I. -Ic:\wint\lib.if8
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