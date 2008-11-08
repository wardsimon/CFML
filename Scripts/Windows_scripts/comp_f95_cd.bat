rem
rem CrysFML for Absoft Compiler (Debug)
rem
   @echo off
   cd ../../cfml
rem
   echo **---- Level 0 ----**
   echo .... Mathematical,String_Utilities, Profile Functions,
   echo .... Graphical, Optimization Modules
rem
   f95 -c -g CFML_math_gen.f90         
   f95 -c -g CFML_string_util.f90      
   f95 -c -g CFML_random.f90           
   f95 -c -g CFML_fft.f90              
   f95 -c -g CFML_fft_nF.f90           
   f95 -c -g CFML_Profile_TOF.f90      
   f95 -c -g CFML_Profile_Finger.f90   
   f95 -c -g CFML_Profile_Functs.f90   
   rem f95 -c -g CFML_optimization.f90     
   f95 -c -g CFML_optimization_lsq.f90 
rem
   echo **---- Level 1 ----**
   echo .... Math_3D, Spher_Harm
rem
   f95 -c -g CFML_math_3D.f90          
   f95 -c -g CFML_spher_harm.f90       
rem
   echo **---- Level 2 ----**
   echo .... Sym_Table, Chem_Scatt, Cryst_Types, DiffPatt Modules
rem
   f95 -c -g CFML_sym_table.f90        
   f95 -c -g CFML_chem_scatt.f90       
   f95 -c -g CFML_cryst_types.f90      
   f95 -c -g CFML_diffpatt.f90         
rem
   echo **---- Level 3 ----**
   echo .... Bonds_Table, Symmetry Modules
rem
   f95 -c -g CFML_bonds_table.f90      
   f95 -c -g CFML_symmetry.f90         
rem
   echo **---- Level 4 ----**
   echo .... Reflct_Util, Atom_Mod Modules
rem
   f95 -c -g CFML_reflct_util.f90      
   f95 -c -g CFML_atom_mod.f90         
rem
   echo **---- Level 5 ----**
   echo .... Sfac, Propagk, Geom_Calc, Maps, Molecules Modules
rem
   f95 -c -g CFML_sfac.f90             
   f95 -c -g CFML_propagk.f90          
   f95 -c -g CFML_geom_calc.f90        
   f95 -c -g CFML_maps.f90             
   f95 -c -g CFML_molecules.f90        
rem
   echo **---- Level 6 ----**
   echo .... Form_CIF, Conf_Calc, RefCodes Modules
rem
   f95 -c -g CFML_form_cif.f90         
   f95 -c -g CFML_conf_calc.f90        
   f95 -c -g CFML_refcodes.f90         
rem
   echo **---- Level 7 ----**
   echo .... Polar, Symm, SFac for Magnetic  Modules
rem   
   f95 -c -g CFML_polar.f90      
   f95 -c -g CFML_magsymm.f90    
   f95 -c -g CFML_msfac.f90      
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
   f95 -c -g CFML_io_mess.f90          
rem
   echo **---- Level 9 ----**
   echo .... Optimization Simulated Annealing Module
rem   
   f95 -c -g CFML_optimization_san.f90 
rem
   echo **---- Crysfml Library: Console (DEBUG) version ----**
rem
   lib -out:crysfml.lib *.obj
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
   cd ..\Scripts\Windows