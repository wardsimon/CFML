rem
rem CrysFML for Absoft Compiler (Optimization)
rem
   @echo off
   cd ../../cfml
rem
   echo **---- Level 0 ----**
   echo .... Mathematical,String_Utilities, Profile Functions,
   echo .... Graphical, Optimization Modules
rem
   f95 -c -O3 -w CFML_math_gen.f90           >  out
   f95 -c -O3 -w CFML_string_util.f90        >> out
   f95 -c -O3 -w CFML_random.f90             >> out
   f95 -c -O3 -w CFML_fft.f90                >> out
   f95 -c -O3 -w CFML_fft_nF.f90             >> out
   f95 -c -O3 -w CFML_Profile_TOF.f90        >> out
   f95 -c -O3 -w CFML_Profile_Finger.f90     >> out
   f95 -c -O3 -w CFML_Profile_Functs.f90     >> out
   f95 -c -O3 -w CFML_optimization.f90       >> out
   f95 -c -O3 -w CFML_optimization_lsq.f90   >> out
rem
   echo **---- Level 1 ----**
   echo .... Math_3D, Spher_Harm
rem
   f95 -c -O3 -w CFML_math_3D.f90            >> out
   f95 -c -O3 -w CFML_spher_harm.f90         >> out
rem
   echo **---- Level 2 ----**
   echo .... Sym_Table, Chem_Scatt, Cryst_Types, DiffPatt Modules
rem
   f95 -c -O0 -w CFML_sym_table.f90          >> out
   f95 -c -O0 -w CFML_chem_scatt.f90         >> out
   f95 -c -O3 -w CFML_cryst_types.f90        >> out
   f95 -c -O3 -w CFML_diffpatt.f90           >> out
rem
   echo **---- Level 3 ----**
   echo .... Bonds_Table, Symmetry Modules
rem
   f95 -c -O0 -w CFML_bonds_table.f90        >> out
   f95 -c -O3 -w CFML_symmetry.f90           >> out
rem
   echo **---- Level 4 ----**
   echo .... Reflct_Util, Atom_Mod Modules
rem
   f95 -c -O3 -w CFML_reflct_util.f90        >> out
   f95 -c -O3 -w CFML_atom_mod.f90           >> out
rem
   echo **---- Level 5 ----**
   echo .... Sfac, Propagk, Geom_Calc, Maps, Molecules Modules
rem
   f95 -c -O3 -w CFML_sfac.f90               >> out
   f95 -c -O3 -w CFML_propagk.f90            >> out
   f95 -c -O3 -w CFML_geom_calc.f90          >> out
   f95 -c -O3 -w CFML_maps.f90               >> out
   f95 -c -O3 -w CFML_molecules.f90          >> out
rem
   echo **---- Level 6 ----**
   echo .... Form_CIF, Conf_Calc, RefCodes Modules
rem
   f95 -c -O3 -w CFML_form_cif.f90           >> out
   f95 -c -O3 -w CFML_conf_calc.f90          >> out
   f95 -c -O3 -w CFML_refcodes.f90           >> out
rem
   echo **---- Level 7 ----**
   echo .... Polar, Symm, SFac for Magnetic  Modules
rem   
   f95 -c -O3 -w CFML_polar.f90              >> out
   f95 -c -O3 -w CFML_magsymm.f90            >> out
   f95 -c -O3 -w CFML_msfac.f90              >> out
rem
   echo **---- Level 8 ----**
   echo .... ILL_Instrm_data, CFML_SXTAL_geom Modules
rem
rem f95 -c -O3 -w CFML_ILL_Instrm_data.f90    >> out
rem f95 -c -O3 -w CFML_SXTAL_geom.f90         >> out
rem
   echo **---- Level 9 ----**
   echo .... IO_Mess Module
rem
   f95 -c -O3 -w CFML_io_mess.f90            >> out
rem
   echo **---- Level 10 ----**
   echo .... Optimization Simulated Annealing Module
rem
   f95 -c -O3 -w CFML_optimization_san.f90   >> out
rem
   echo **---- Crysfml Library: Console version ----**
rem
   lib -out:crysfml.lib *.obj
rem
   echo **---- Absoft Directory ----**
rem
   if not exist ..\..\Absoft mkdir ..\..\Absoft
   if exist ..\..\Absoft\LibC rmdir ..\..\Absoft\LibC /S /Q
   mkdir ..\..\Absoft\LibC
rem
   copy *.mod ..\..\Absoft\LibC > nul
   move *.lib   ..\..\Absoft\LibC > nul
   del *.obj *.mod *.lst *.bak > nul
rem
   cd ..\Scripts\Windows
   