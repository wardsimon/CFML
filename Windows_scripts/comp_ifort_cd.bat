@echo off
rem
echo **---- Level 0 ----**
echo .... Mathematical,String_Utilities, Profile Functions,
echo .... Graphical, Optimization Modules
rem
ifort /c CFML_math_gen.f90         /debug:full /check /traceback /nologo
ifort /c CFML_string_util.f90      /debug:full /check /traceback /nologo
ifort /c CFML_random.f90           /debug:full /check /traceback /nologo
ifort /c CFML_fft.f90              /debug:full /check /traceback /nologo
ifort /c CFML_fft_nF.f90           /debug:full /check /traceback /nologo
ifort /c CFML_Profile_TOF.f90      /debug:full /check /traceback /nologo
ifort /c CFML_Profile_Finger.f90   /debug:full /check /traceback /nologo
ifort /c CFML_Profile_Functs.f90   /debug:full /check /traceback /nologo
ifort /c CFML_optimization.f90     /debug:full /check /traceback /nologo
ifort /c CFML_optimization_lsq.f90 /debug:full /check /traceback /nologo
rem
echo **---- Level 1 ----**
echo .... Math_3D, Spher_Harm
rem
ifort /c CFML_math_3D.f90          /debug:full /check /traceback /nologo
ifort /c CFML_spher_harm.f90       /debug:full /check /traceback /nologo
rem
echo **---- Level 2 ----**
echo .... Sym_Table, Chem_Scatt, Cryst_Types, DiffPatt Modules
rem
ifort /c CFML_sym_table.f90        /debug:full /check /traceback /nologo
ifort /c CFML_chem_scatt.f90       /debug:full /check /traceback /nologo
ifort /c CFML_cryst_types.f90      /debug:full /check /traceback /nologo
ifort /c CFML_diffpatt.f90         /debug:full /check /traceback /nologo
rem
echo **---- Level 3 ----**
echo .... Bonds_Table, Symmetry Modules
rem
ifort /c CFML_bonds_table.f90      /debug:full /check /traceback /nologo
ifort /c CFML_symmetry.f90         /debug:full /check /traceback /nologo
rem
echo **---- Level 4 ----**
echo .... Reflct_Util, Atom_Mod Modules
rem
ifort /c CFML_reflct_util.f90      /debug:full /check /traceback /nologo
ifort /c CFML_atom_mod.f90         /debug:full /check /traceback /nologo
rem
echo **---- Level 5 ----**
echo .... Sfac, Propagk, Geom_Calc, Maps, Molecules Modules
rem
ifort /c CFML_sfac.f90             /debug:full /check /traceback /nologo
ifort /c CFML_propagk.f90          /debug:full /check /traceback /nologo
ifort /c CFML_geom_calc.f90        /debug:full /check /traceback /nologo
ifort /c CFML_maps.f90             /debug:full /check /traceback /nologo
ifort /c CFML_molecules.f90        /debug:full /check /traceback /nologo
rem
echo **---- Level 6 ----**
echo .... Form_CIF, Conf_Calc, RefCodes Modules
rem
ifort /c CFML_form_cif.f90         /debug:full /check /traceback /nologo
ifort /c CFML_conf_calc.f90        /debug:full /check /traceback /nologo
ifort /c CFML_refcodes.f90         /debug:full /check /traceback /nologo
rem
echo **---- Level 7 ----**
echo .... Polar, Symm, SFac for Magnetic  Modules
ifort /c CFML_polar.f90            /debug:full /check /traceback /nologo
ifort /c CFML_magsymm.f90          /debug:full /check /traceback /nologo
ifort /c CFML_msfac.f90            /debug:full /check /traceback /nologo
rem
echo **---- Level 7 ----**
echo .... ILL_Instrm_data, CFML_SXTAL_geom Modules
rem
ifort /c CFML_ILL_Instrm_data.f90  /debug:full /check /traceback /nologo
ifort /c CFML_SXTAL_geom.f90       /debug:full /check /traceback /nologo
rem
echo **---- Level 8 ----**
echo .... IO_Mess Module
rem
ifort /c CFML_io_mess.f90          /debug:full /check /traceback /nologo
rem
echo **---- Level 9 ----**
echo .... Optimization Simulated Annealing Module
ifort /c CFML_optimization_san.f90 /debug:full /check /traceback /nologo
rem
echo **---- Crysfml Library: Console (DEBUG) version ----**
rem
lib /out:crysfml.lib *.obj
rem
rem
echo **---- Intel Directory ----**
rem
   if not exist ..\..\Intel mkdir ..\..\Intel
   if exist ..\..\Intel\LibC rmdir ..\..\Intel\LibC /S /Q
   mkdir ..\..\Intel\LibC
rem
   copy *.mod ..\..\Intel\LibC > nul
   move *.lib ..\..\Intel\LibC > nul
   del *.obj *.mod *.lst *.bak > nul
rem 