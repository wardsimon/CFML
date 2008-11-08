@echo off
rem
echo **---- Level 0 ----**
echo .... Mathematical,String_Utilities, Profile Functions,
echo .... Graphical, Optimization Modules
rem
lf95 -c f2kcli.f90                -g  -chk   >  out
lf95 -c CFML_math_gen.f90         -g  -chk   >> out
lf95 -c CFML_string_util.f90      -g  -chk   >> out
lf95 -c CFML_random.f90           -g  -chk   >> out
lf95 -c CFML_fft.f90              -g  -chk   >> out
lf95 -c CFML_fft_nF.f90           -g  -chk   >> out
lf95 -c CFML_Profile_TOF.f90      -g  -chk   >> out
lf95 -c CFML_Profile_Finger.f90   -g  -chk   >> out
lf95 -c CFML_Profile_Functs.f90   -g  -chk   >> out
lf95 -c CFML_optimization.f90     -g  -chk   >> out
lf95 -c CFML_optimization_lsq.f90 -g  -chk   >> out
rem
echo **---- Level 1 ----**
echo .... Math_3D, Spher_Harm
rem
lf95 -c CFML_math_3D.f90          -g  -chk   >> out
lf95 -c CFML_spher_harm.f90       -g  -chk   >> out
rem
echo **---- Level 2 ----**
echo .... Sym_Table, Chem_Scatt, Cryst_Types, DiffPatt Modules
rem
lf95 -c CFML_sym_table.f90        -g  -chk   >> out
lf95 -c CFML_chem_scatt.f90       -g  -chk   >> out
lf95 -c CFML_cryst_types.f90      -g  -chk   >> out
lf95 -c CFML_diffpatt.f90         -g  -chk   >> out
rem
echo **---- Level 3 ----**
echo .... Bonds_Table, Symmetry Modules
rem
lf95 -c CFML_bonds_table.f90      -g  -chk   >> out
lf95 -c CFML_symmetry.f90         -g  -chk   >> out
rem
echo **---- Level 4 ----**
echo .... Reflct_Util, Atom_Mod Modules
rem
lf95 -c CFML_reflct_util.f90      -g  -chk   >> out
lf95 -c CFML_atom_mod.f90         -g  -chk   >> out
rem
echo **---- Level 5 ----**
echo .... Sfac, Propagk, Geom_Calc, Maps, Molecules Modules
rem
lf95 -c CFML_sfac.f90             -g  -chk   >> out
lf95 -c CFML_propagk.f90          -g  -chk   >> out
lf95 -c CFML_geom_calc.f90        -g  -chk   >> out
lf95 -c CFML_maps.f90             -g  -chk   >> out
lf95 -c CFML_molecules.f90        -g  -chk   >> out
rem
echo **---- Level 6 ----**
echo .... Form_CIF, Conf_Calc, RefCodes Modules
rem
lf95 -c CFML_form_cif.f90         -g  -chk   >> out
lf95 -c CFML_conf_calc.f90        -g  -chk   >> out
lf95 -c CFML_refcodes.f90         -g  -chk   >> out
rem
echo **---- Level 7 ----**
echo .... Polar, Symm, SFac for Magnetic  Modules
lf95 -c CFML_polar.f90            -g -chk   >> out
lf95 -c CFML_magsymm.f90          -g -chk   >> out
lf95 -c CFML_msfac.f90            -g -chk   >> out
rem
echo **---- Level 7 ----**
echo .... ILL_Instrm_data, CFML_SXTAL_geom Modules
rem
lf95 -c CFML_ILL_Instrm_data.f90  -g  -chk   >> out
lf95 -c CFML_SXTAL_geom.f90       -g  -chk   >> out
rem
echo **---- Level 8 ----**
echo .... IO_Mess Module
rem
lf95 -c CFML_io_messwin.f90       -g -chk -mod .;c:\wint\lib.l95 >> out
rem
echo **---- Level 9 ----**
echo .... Optimization Simulated Annealing Module
lf95 -c CFML_optimization_san.f90 -g -chk -mod .;c:\wint\lib.l95 >> out
rem
echo **---- Crysfml Library: Winteracter (DEBUG) version ----**
rem
lm @lib_modw.lnk
rem
echo **---- Lahey Directory ----**
rem
   if not exist ..\..\Lahey mkdir ..\..\Lahey
   if exist ..\..\Lahey\LibW rmdir ..\..\Lahey\LibW /S /Q
   mkdir ..\..\Lahey\LibW
rem
   copy *.mod ..\..\Lahey\LibW > nul
   move *.lib ..\..\Lahey\LibW > nul
   del *.obj *.mod *.lst *.bak > nul
rem 