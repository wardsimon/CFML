@echo off
rem
echo **---- Level 0 ----**
echo .... Mathematical,String_Utilities, Profile Functions,
echo .... Graphical, Optimization Modules
rem
f95 -c -g CFML_Math_Gen.f90
f95 -c -g CFML_String_Util.f90
f95 -c -g CFML_Random.f90
f95 -c -g CFML_FFT.f90
f95 -c -g CFML_FFT_nF.f90
f95 -c -g CFML_Profile_TOF.f90
f95 -c -g CFML_Profile_Finger.f90
f95 -c -g CFML_Profile_Functs.f90
rem f95 -c -g CFML_optimization.f90
f95 -c -g CFML_optimization_lsq.f90
rem
echo **---- Level 1 ----**
echo .... Math_3D, Spher_Harm
rem
f95 -c -g CFML_Math_3D.f90
f95 -c -g CFML_Spher_Harm.f90
rem
echo **---- Level 2 ----**
echo .... Sym_Table, Chem_Scatt, Cryst_Types, DiffPatt Modules
rem
f95 -c -g CFML_Sym_Table.f90
f95 -c -g CFML_Chem_Scatt.f90
f95 -c -g CFML_Cryst_Types.f90
f95 -c -g CFML_Diffpatt.f90
rem
echo **---- Level 3 ----**
echo .... Bonds_Table, Symmetry Modules
rem
f95 -c -g CFML_Bonds_Table.f90
f95 -c -g CFML_Symmetry.f90
rem
echo **---- Level 4 ----**
echo .... Reflct_Util, Atom_Mod Modules
rem
f95 -c -g CFML_Reflct_Util.f90
f95 -c -g CFML_Atom_Mod.f90
rem
echo **---- Level 5 ----**
echo .... Sfac, Propagk, Geom_Calc, Maps, Molecules Modules
rem
f95 -c -g CFML_Sfac.f90
f95 -c -g CFML_Propagk.f90
f95 -c -g CFML_Geom_Calc.f90
f95 -c -g CFML_Maps.f90
f95 -c -g CFML_Molecules.f90
rem
echo **---- Level 6 ----**
echo .... Form_CIF, Conf_Calc, RefCodes Modules
rem
f95 -c -g CFML_Form_CIF.f90
f95 -c -g CFML_Conf_Calc.f90
f95 -c -g CFML_Refcodes.f90
rem
echo **---- Level 7 ----**
echo .... Polar, Symm, SFac for Magnetic  Modules
f95 -c -g CFML_Polar.f90
f95 -c -g CFML_MagSymm.f90
f95 -c -g CFML_Msfac.f90
rem
echo **---- Level 7 ----**
echo .... ILL_Instrm_data, CFML_SXTAL_geom Modules
rem
rem f95 -c -g CFML_ILL_Instrm_data.f90
rem f95 -c -g CFML_SXTAL_geom.f90
rem
echo **---- Level 8 ----**
echo .... IO_Mess Module
rem
f95 -c -g CFML_IO_Mess.f90
rem
echo **---- Level 9 ----**
echo .... Optimization Simulated Annealing Module
f95 -c -g CFML_Optimization_SAn.f90
rem
echo **---- Crysfml Library: Console (DEBUG) version ----**
rem
lib -out:crysfml.lib *.obj
rem
rem
echo **---- Intel Directory ----**
rem
   if not exist ..\..\Absoft mkdir ..\..\Absoft
   if exist ..\..\Absoft\LibC rmdir ..\..\Absoft\LibC /S /Q
   mkdir ..\..\Absoft\LibC
rem
   copy *.mod ..\..\Absoft\LibC > nul
   move *.lib ..\..\Absoft\LibC > nul
   del *.obj *.mod *.lst *.bak > nul
rem
