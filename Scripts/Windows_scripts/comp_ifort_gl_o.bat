rem
rem CrysFGL for Intel Compiler (Optimization) + WINTERACTER
rem
   @echo off
   cd ../../cfgl
rem
   ifort /c WCFGL_constant.f90           /Ox /nologo /I. /I..\..\Intel\LibW
   ifort /c WCFGL_atomic_table.f90       /Ox /nologo /I. /I..\..\Intel\LibW
   ifort /c WCFGL_objects_definition.f90 /Ox /nologo /I. /Ic:\wint\lib.if8
   ifort /c WCFGL_geometry.f90           /Ox /nologo /I. /I..\..\Intel\LibW
   ifort /c WCFGL_Chull3D.f90            /Ox /nologo /I. /Ic:\wint\lib.if8
   ifort /c WCFGL_quaternion.f90         /Ox /nologo /I. /I..\..\Intel\LibW
   ifort /c WCFGL_trackball.f90          /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_metrix.f90             /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_glatom.f90             /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_glbond.f90             /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_atom_tree.f90          /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_matom_tree.f90         /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_bond_tree.f90          /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_glpoly.f90             /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_poly_tree.f90          /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_display.f90            /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
   ifort /c WCFGL_IO.f90                 /Ox /nologo /I. /Ic:\wint\lib.if8 /I..\..\Intel\LibW >> out
rem
   echo **---- CrysFGL Library: Winteracter version ----**
rem
   lib /out:wcrysfgl.lib *.obj 
rem
   if exist ..\..\Intel\LibGL rmdir ..\..\Intel\LibGL /S /Q
   mkdir ..\..\Intel\LibGL
rem
   copy *.mod ..\..\Intel\LibGL > nul
   move *.lib ..\..\Intel\LibGL > nul
   del *.obj *.mod *.bak *.lst > nul
rem
   cd ..\Scripts\Windows 