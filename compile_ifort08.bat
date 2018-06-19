@echo off
echo.
echo -------------------------------------------------------
echo ---- Crystallographic Fortran Modules Library 2018 ----
echo ---- Compiler:  Intel Fortran Compiler XE_2017     ---- 
echo ---- CrysFML Team                                  ----
echo -------------------------------------------------------
rem
rem ---- INIT ----
   (set _DEBUG=N)
   (set _WINTER=N)
   if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
rem
rem ----
rem ---- Arguments 
rem ----
:LOOP
    if [%1]==[debug]  (set _DEBUG=Y)
    if [%1]==[winter] (set _WINTER=Y)
    shift
    if not [%1]==[] goto LOOP
rem
rem ---- Options
rem
   if [%_DEBUG%]==[Y] (
      if [%TARGET_ARCH%]==[ia32] (set DIRECTORY=ifort_debug) else (set DIRECTORY=ifort64_debug)
      (set OPT0=/debug:full /check /traceback /CB)
      (set OPT1=/debug:full /check /traceback /CB)
   ) else (
      if [%TARGET_ARCH%]==[ia32] (set DIRECTORY=ifort) else (set DIRECTORY=ifort64)
      (set OPT0=/Od)
      (set OPT1=/O2)
rem --- Testing JGP
      (set OPT1=/O3)
   )
   (set OPT2=/Qopt-report:0)
   (set OPT3=)
rem   
   if [%_WINTER%]==[Y] (
      if [%TARGET_ARCH%]==[ia32] (
         (set OPT3=/I%WINTER%\lib.if8)
      ) else (
         (set OPT3=/I%WINTER%\lib.i64)
      )
   )
rem
   cd %CRYSFML%\Src08
rem
   echo.
   echo ----
   echo ----  Set compiler options  
   echo ----
   echo OPT0:%OPT0%
   echo OPT1:%OPT1%
   echo OPT2:%OPT2%
   echo OPT3:%OPT3%
   echo ----
   echo.
   goto T1
rem
   echo .... Global Dependencies for CFML
rem
   ifort /c CFML_GlobalDeps_Windows_IFOR.f90     /nologo %OPT0% %OPT2% /module:.\mod  
rem
   echo .... Mathematics Procedures
rem
   ifort /c CFML_Math_General.f90                /nologo %OPT1% %OPT2% /module:.\mod
rem 
rem   Submodulos CFML_Math_General   
      cd .\CFML_Math_General
      ifort /c AM_Median.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Co_Linear.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Co_Prime.f90                      /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Debye.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Determinant.f90                   /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Diagonalize_SH.f90                /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Equal_Matrix.f90                  /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Equal_Vector.f90                  /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Erfc_Deriv.f90                    /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Factorial.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Invert_Matrix.f90                 /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c In_Limits.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Linear_Dependent.f90              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Locate.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Lower_Triangular.f90              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Modulo_Lat.f90                    /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Negligible.f90                    /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Norm.f90                          /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Outerprod.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Pgcd.f90                          /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Points_In_Line2D.f90              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Poly_Legendre.f90                 /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Ppcm.f90                          /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Rank.f90                          /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Scalar.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c SmoothingVec.f90                  /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Sort.f90                          /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Sph_Jn.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Spline.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Svdcmp.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Swap.f90                          /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Trace.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Upper_Triangular.f90              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Zbelong.f90                       /nologo %OPT1% %OPT2%  /module:..\mod  
      move /y *.obj .. > nul
      cd ..
rem      
   ifort /c CFML_math_3D.f90                     /nologo %OPT1% %OPT2% /module:.\mod
rem 
rem   Submodulos CFML_Math_3D
      cd CFML_math_3D
      ifort /c Cross_Product.f90                 /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Determ_3x3.f90                    /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Determ_Vec.f90                    /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Get_Cart_from_Cylin.f90           /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Get_Cart_from_Spher.f90           /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Get_Cylin_from_Cart.f90           /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Get_Spher_from_Cart.f90           /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Invert_Array3x3.f90               /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Matrix_DiagEigen.f90              /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Mat_Cross.f90                     /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Resolv_Sistem.f90                 /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Rotation_Axes.f90                 /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Tensor_Product.f90                /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Vec_Length.f90                    /nologo %OPT1% %OPT2% /module:..\mod
      move /y *.obj .. > nul
      cd ..
rem 
   ifort /c CFML_Spher_Harm.f90                  /nologo %OPT1% %OPT2% /module:.\mod
   ifort /c CFML_Random.f90                      /nologo %OPT1% %OPT2% /module:.\mod
   ifort /c CFML_Ffts.f90                        /nologo %OPT1% %OPT2% /module:.\mod
rem
   echo .... Profiles Functions
rem
   ifort /c CFML_Profile_Functs.f90              /nologo %OPT1% %OPT2% /module:.\mod
   ifort /c CFML_Profile_Finger.f90              /nologo %OPT1% %OPT2% /module:.\mod
   ifort /c CFML_Profile_TOF.f90                 /nologo %OPT1% %OPT2% /module:.\mod
rem   
   ifort /c CFML_Extinction_Correction.f90       /nologo %OPT1% %OPT2% /module:.\mod
rem
   echo .... IO Messages /String Utilities
rem
   if [%_WINTER%]==[Y] (
     ifort /c CFML_io_messwin.f90         /nologo %OPT1% %OPT2% %OPT3% /module:.\mod
   ) else (
     ifort /c CFML_io_mess.f90                   /nologo %OPT1% %OPT2% /module:.\mod
   )
   ifort /c CFML_string_util.f90                 /nologo %OPT1% %OPT2% /module:.\mod
rem  
rem   Submodulos CFML_String_Utilities
      cd  CFML_String_Utilities
      ifort /c StringFullp.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c StringTools.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c StringNum.f90                     /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c StringReadKey.f90                 /nologo %OPT1% %OPT2% /module:..\mod
      move /y *.obj .. > nul
      cd .. 
rem
   echo .... Tables definitions
rem
   ifort /c CFML_BVSpar.f90                      /nologo %OPT0% %OPT2% /module:.\mod
   ifort /c CFML_chem_scatt.f90                  /nologo %OPT0% %OPT2% /module:.\mod
   ifort /c CFML_bonds_table.f90                 /nologo %OPT0% %OPT2% /module:.\mod
   ifort /c CFML_sym_table.f90                   /nologo %OPT0% %OPT2% /module:.\mod
rem 
   echo .... Rational Arithmetic
rem
   ifort /c CFML_Rational_Arithmetic.f90         /nologo %OPT1% %OPT2% /module:.\mod
rem  
rem   Submodulos CFML_Rational_Arithmetic
      cd CFML_Rational_Arithmetic   
      ifort /c assignment.f90                    /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c constructor.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c matmul.f90                        /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c maxloc.f90                        /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c maxval.f90                        /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c minloc.f90                        /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c minval.f90                        /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c mod.f90                           /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c modulo.f90                        /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_add.f90                  /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_division.f90             /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_eq.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_ge.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_gt.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_le.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_lt.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_minus.f90                /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_multiply.f90             /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c operator_ne.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c overloads.f90                     /nologo %OPT1% %OPT2% /module:..\mod
      move /y *.obj .. > nul
      cd ..
rem
   echo .... Crystal Metrics
rem
   ifort /c CFML_crystal_metrics.f90             /nologo %OPT1% %OPT2% /module:.\mod
rem  
rem   Submodulos CFML_Crystal_Metrics
      cd CFML_Crystal_Metrics
      ifort /c genmetrics.f90                    /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ioroutines.f90                    /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ThConver.f90                      /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c Nigglicell.f90                    /nologo %OPT1% %OPT2% /module:..\mod
      move /y *.obj .. > nul
      cd ..  
rem
   echo .... Patterns Information
rem
:T1
   ifort /c CFML_diffpatt.f90                    /nologo %OPT1% %OPT2% /module:.\mod
rem  
rem   Submodulos CFML_Diffraction_Patterns

      cd CFML_Diffraction_Patterns
      ifort /c ReadPatt_ILL.f90                  /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ReadPatt_LLB.f90                  /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ReadPatt_NLS.f90                  /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ReadPatt_PSI.f90                  /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ReadPatt_ISIS.f90                 /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ReadPatt_Socabim.f90              /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ReadPatt_GSAS.f90                 /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c ReadPatt_PAN.f90                  /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c WritePatterns.f90                 /nologo %OPT1% %OPT2% /module:..\mod
      ifort /c GenPatterns.f90                   /nologo %OPT1% %OPT2% /module:..\mod
      move /y *.obj .. > nul
      cd ..
      
   goto TTT
   
   ifort /c CFML_ILL_Instrm_data.f90                  /nologo %OPT1% %OPT2%
   goto FIN
   echo .... Optimization procedures
rem
rem   ifort /c CFML_optimization.f90                     /nologo %OPT1% %OPT2%
rem   ifort /c CFML_optimization_lsq.f90                 /nologo %OPT1% %OPT2%
rem
   echo .... Patterns Information
rem
   ifort /c CFML_diffpatt.f90                         /nologo %OPT1% %OPT2%
   goto TTT
   ifort /c CFML_ILL_Instrm_data.f90                  /nologo %OPT1% %OPT2%
   goto FIN
rem
   echo .... Crystal Metrics
rem
   ifort /c CFML_symmetry.f90                         /nologo %OPT1% %OPT2%
rem
   echo **---- Level 3 ----**
   echo .... EoS, Reflections, Atoms
rem
   ifort /c CFML_Eos_Mod.f90                          /nologo %OPT1% %OPT2%
   ifort /c CFML_reflct_util.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_atom_mod.f90                         /nologo %OPT1% %OPT2%
rem
   echo **---- Level 4 ----**
   echo .... Formats, Geometry, Molecules
rem
   ifort /c CFML_geom_calc.f90                        /nologo %OPT1% %OPT2%
   ifort /c CFML_molecules.f90                        /nologo %OPT1% %OPT2%
   ifort /c CFML_form_cif.f90                        /nologo %OPT1% %OPT2%
rem
   echo **---- Level 5 ----**
   echo .... Extinction, Structure Factors, SXTAL geometry, Propag Vectors
rem
   ifort /c CFML_sfac.f90                            /nologo %OPT1% %OPT2%
   ifort /c CFML_sxtal_Geom.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_propagk.f90                         /nologo %OPT1% %OPT2%
rem
   echo **---- Level 6 ----**
   echo .... Maps, BVS, Energy Configurations
rem
   ifort /c CFML_Export_Vtk.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_maps.f90                            /nologo %OPT1% %OPT2%
   ifort /c CFML_conf_calc.f90                       /nologo %OPT1% %OPT2%
rem
   echo **---- Level 7 ----**
   echo .... Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem
   ifort /c CFML_Magnetic_Groups.f90                 /nologo %OPT1% %OPT2%
   ifort /c CFML_magsymm.f90                         /nologo %OPT1% %OPT2%
   ifort /c CFML_optimization_san.f90                /nologo %OPT1% %OPT2% %OPT3%
   ifort /c CFML_refcodes.f90                        /nologo %OPT1% %OPT2%
rem
   echo **---- Level 8 ----**
   echo .... Magnetic Structure Factors, Polarimetry
rem
   ifort /c CFML_msfac.f90                           /nologo %OPT1% %OPT2%
   ifort /c CFML_polar.f90                           /nologo %OPT1% %OPT2%
rem
rem
:TTT
   echo **---- Crysfml Library ----**
rem
   if [%_WINTER%]==[Y] (
     lib /out:wcrysfml.lib *.obj
   ) else (
     lib /out:crysfml.lib *.obj
   )
rem
   echo **---- ifort Directory ----**
rem
   if not exist ..\%DIRECTORY% mkdir ..\%DIRECTORY%
   if [%_WINTER%]==[Y] (
     if exist ..\%DIRECTORY%\LibW08 rmdir ..\%DIRECTORY%\LibW08 /S /Q
     mkdir ..\%DIRECTORY%\LibW08
     copy .\mod\*.mod ..\%DIRECTORY%\LibW08 > nul
     copy .\mod\*.smod ..\%DIRECTORY%\LibW08 > nul
     move *.lib ..\%DIRECTORY%\LibW08 > nul
   ) else (
     if exist ..\%DIRECTORY%\LibC08 rmdir ..\%DIRECTORY%\LibC08 /S /Q
     mkdir ..\%DIRECTORY%\LibC08
     copy .\mod\*.mod ..\%DIRECTORY%\LibC08 > nul
     copy .\mod\*.smod ..\%DIRECTORY%\LibC08 > nul
     move *.lib ..\%DIRECTORY%\LibC08 > nul
   )
   del *.obj *.lst *.bak > nul
rem
rem   cd %CRYSFML%\Scripts\Windows
rem
:FIN
   cd ..
