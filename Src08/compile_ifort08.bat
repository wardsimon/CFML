@echo off
rem ------------------------------------
rem ---- CrysFML for Intel Compiler ----
rem ---- JGP                   2019 ----
rem ------------------------------------
rem
rem ---- INIT ----
   (set _DEBUG=N)
   (set _WINTER=N)
   if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
rem
rem ----
rem ---- Arguments ----
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
   )
rem   
   (set OPT2=/fpp /Qopt-report:0)
   (set OPT3=)
   if [%_WINTER%]==[Y] (
      if [%TARGET_ARCH%]==[ia32] (
         (set OPT3=/I%WINTER%\lib.if8)
      ) else (
         (set OPT3=/I%WINTER%\lib.i64)
      )
   )
rem
   cd %CRYSFML%\Src08N
rem
   echo.
   echo **-------------------------------------**
   echo **---- CrysFML: Start Compilation  ----**
   echo **-------------------------------------**
   echo Compiler Options 
   echo OPT0:%OPT0%
   echo OPT1:%OPT1%
   echo OPT2:%OPT2%
   echo OPT3:%OPT3%
   echo.
rem
   echo .... Global Dependencies for CFML
rem
   ifort /c CFML_GlobalDeps_Windows_IFOR.f90         /nologo %OPT1% %OPT2% /module:.\mod
rem
   echo .... Mathematics Procedures
   ifort /c CFML_Maths.f90                           /nologo %OPT1% %OPT2% /module:.\mod
rem 
rem   Submodules CFML_Maths   
      cd .\CFML_Maths
      ifort /c Co_Prime.f90                          /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Co_Linear.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Cross_Product.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Debye.f90                             /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Determinant.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Diagonalize_SH.f90                    /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Equal_Matrix.f90                      /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Equal_Vector.f90                      /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Factorial.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Invert_Matrix.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c In_Limits.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Is_Diagonal_Matrix.f90                /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Is_Null_Vector.f90                    /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Linear_Dependent.f90                  /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Locate.f90                            /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Lower_Triangular.f90                  /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c ManipVec.f90                          /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Mat_Cross.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Modulo_Lat.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Negligible.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Norm.f90                              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Outerprod.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Pgcd.f90                              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Points_In_Line2D.f90                  /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Poly_Legendre.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Rank.f90                              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Resolv_System.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Rotation_Axes.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c RowEchelon.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Scalar.f90                            /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c SistCoord_Changes.f90                 /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Sort.f90                              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Swap.f90                              /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Tensor_Product.f90                    /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Trace.f90                             /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Upper_Triangular.f90                  /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Vec_Length.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Zbelong.f90                           /nologo %OPT1% %OPT2%  /module:..\mod
      move /y *.obj .. > nul
      cd ..
rem      
   echo .... Strings Procedures
   ifort /c CFML_StringsUtil.f90                     /nologo %OPT1% %OPT2% /module:.\mod
rem 
rem   Submodules CFML_Strings   
      cd .\CFML_Strings
      ifort /c StringFullp.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c StringNum.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c StringReadKey.f90                     /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c StringTools.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      move /y *.obj .. > nul
      cd ..
rem         
   echo .... Rational arithmetics
   ifort /c CFML_RationalArith.f90                  /nologo %OPT1% %OPT2% /module:.\mod
rem 
rem   Submodules CFML_Rational   
      cd .\CFML_Rational 
      ifort /c Assignment.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Constructor.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Equal_rational.f90                    /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Generic.f90                           /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Is_Integer.f90                        /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_add.f90                      /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_eq.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_divisor.f90                  /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_ge.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_gt.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_le.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_lt.f90                       /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_minus.f90                    /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_multiply.f90                 /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Operator_neq.f90                      /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Overloads.f90                         /nologo %OPT1% %OPT2%  /module:..\mod
      ifort /c Rational_RowEchelon.f90               /nologo %OPT1% %OPT2%  /module:..\mod
      move /y *.obj .. > nul
      cd ..        
      goto END
   
   ifort /c CFML_LSQ_TypeDef.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_spher_harm.f90                       /nologo %OPT1% %OPT2%
   ifort /c CFML_random.f90                           /nologo %OPT1% %OPT2%
   ifort /c CFML_ffts.f90                             /nologo %OPT1% %OPT2%
   ifort /c CFML_string_util.f90                      /nologo %OPT1% %OPT2%
rem   ifort /c CFML_Rational_Arithmetic.f90              /nologo %OPT1% %OPT2%
   if [%_WINTER%]==[Y] (
     ifort /c CFML_io_messwin.f90                     /nologo %OPT1% %OPT2% %OPT3%
   ) else (
     ifort /c CFML_io_mess.f90                        /nologo %OPT1% %OPT2%
   )
   ifort /c CFML_Profile_TOF.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_Profile_Finger.f90                   /nologo %OPT1% %OPT2%
   ifort /c CFML_Profile_Functs.f90                   /nologo %OPT1% %OPT2%
rem
   echo .... Mathematical(II), Optimization, Tables, Patterns
rem
   ifort /c CFML_math_3D.f90                          /nologo %OPT1% %OPT2%
   ifort /c CFML_optimization.f90                     /nologo %OPT1% %OPT2%
   ifort /c CFML_optimization_lsq.f90                 /nologo %OPT1% %OPT2%
   ifort /c CFML_sym_table.f90                        /nologo %OPT0% %OPT2%
   ifort /c CFML_chem_scatt.f90                       /nologo %OPT0% %OPT2%
   ifort /c CFML_BVSpar.f90                           /nologo %OPT0% %OPT2%
   ifort /c CFML_diffpatt.f90                         /nologo %OPT1% %OPT2%
rem
   echo .... Bonds, Crystal Metrics, Symmetry, ILL_Instr
rem
   ifort /c CFML_bonds_table.f90                      /nologo %OPT0% %OPT2%
   ifort /c CFML_cryst_types.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_symmetry.f90                         /nologo %OPT1% %OPT2%
   ifort /c CFML_ILL_Instrm_data.f90                  /nologo %OPT1% %OPT2%
rem
   echo .... EoS, Reflections, Atoms
rem
   ifort /c CFML_Eos_Mod.f90                          /nologo %OPT1% %OPT2%
   ifort /c CFML_reflct_util.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_atom_mod.f90                         /nologo %OPT1% %OPT2%
rem
   echo .... Formats, Geometry, Molecules
rem
   ifort /c CFML_geom_calc.f90                        /nologo %OPT1% %OPT2%
   ifort /c CFML_molecules.f90                        /nologo %OPT1% %OPT2%
   ifort /c CFML_form_cif.f90                        /nologo %OPT1% %OPT2%
rem
   echo .... Extinction, Structure Factors, SXTAL geometry, Propag Vectors
rem
   ifort /c CFML_Extinction_Correction.f90           /nologo %OPT1% %OPT2%
   ifort /c CFML_sfac.f90                            /nologo %OPT1% %OPT2%
   ifort /c CFML_sxtal_Geom.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_propagk.f90                         /nologo %OPT1% %OPT2%
rem
   echo .... Maps, BVS, Energy Configurations
rem
   ifort /c CFML_Export_Vtk.f90                      /nologo %OPT1% %OPT2%
   ifort /c CFML_maps.f90                            /nologo %OPT1% %OPT2%
   ifort /c CFML_conf_calc.f90                       /nologo %OPT1% %OPT2%
rem
   echo .... Magnetic Symmetry, Simulated Annealing, Keywords Parser
rem
   ifort /c CFML_Magnetic_Groups.f90                 /nologo %OPT1% %OPT2%
   ifort /c CFML_magsymm.f90                         /nologo %OPT1% %OPT2%
   ifort /c CFML_optimization_san.f90                /nologo %OPT1% %OPT2% %OPT3%
   ifort /c CFML_refcodes.f90                        /nologo %OPT1% %OPT2%
rem
   echo .... Magnetic Structure Factors, Polarimetry
rem
   ifort /c CFML_msfac.f90                           /nologo %OPT1% %OPT2%
   ifort /c CFML_polar.f90                           /nologo %OPT1% %OPT2%
rem
:END
rem
   echo.
   echo Creating CrysFML Library 
rem
   if [%_WINTER%]==[Y] (
     lib /out:wcrysfml.lib *.obj
   ) else (
     lib /out:crysfml.lib *.obj
   )
rem
   echo Creating IFORT directory 
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
   del *.obj  *.lst *.bak > nul
rem
:FIN
   echo.
   echo **---- End Compilation for CrysFML ----**
   echo.
rem
