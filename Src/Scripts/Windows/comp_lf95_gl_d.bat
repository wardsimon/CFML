rem
rem CrysFGL for Lahey Compiler (Debug) + WINTERACTER
rem
   @echo off
   cd ../../cfgl
rem
   lf95 -c  WCFGL_constant.f90           -g -chk -mod .;..\..\Lahey\LibW > out
   lf95 -c  WCFGL_atomic_table.f90       -g -chk -mod .;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_objects_definition.f90 -g -chk -mod .;c:\wint\lib.l95  >> out
   lf95 -c  WCFGL_geometry.f90    -g -chk -mod .;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_3Dchull.f90     -g -chk -mod .;c:\wint\lib.l95  >> out
   lf95 -c  WCFGL_quaternion.f90  -g -chk -mod .;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_trackball.f90   -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_metrix.f90      -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_glatom.f90      -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_glbond.f90      -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_atom_tree.f90   -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_matom_tree.f90  -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_bond_tree.f90   -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_glpoly.f90      -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_poly_tree.f90   -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_display.f90     -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
   lf95 -c  WCFGL_IO.f90   -g -chk -mod .;c:\wint\lib.l95;..\..\Lahey\LibW >> out
rem
   echo **---- CrysFGL Library: Winteracter version ----**
rem
   lm @..\Scripts\Windows\crysfgl.lnk
rem
   if exist ..\..\Lahey\LibGL rmdir ..\..\Lahey\LibGL /S /Q
   mkdir ..\..\Lahey\LibGL
rem
   copy *.mod ..\..\Lahey\LibGL > nul
   move *.lib ..\..\Lahey\LibGL > nul
   del *.obj *.mod *.bak *.lst > nul
rem
   cd ..\Scripts\Windows
