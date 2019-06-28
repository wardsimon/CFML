@echo on
rem ****
rem ****---- Compilation for Testin_SSG Program ----****
rem ****
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   if x%1 == xifort64    goto IFORT64
   if x%1 == xifort64d    goto IFORT64D
   goto END
rem
rem
rem ****---- Intel Compiler ----****
:IFORTD
   ifort /c CFML_ssg_datafile.f90             /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC /warn
   ifort /c CFML_Rational_Groups.f90          /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC /warn
   ifort /c CFML_Standard_Sett.f90            /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC /warn
   ifort /c CFML_SuperSpaceGroups.f90         /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC /warn
   ifort /c test_ref.f90                      /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC /warn
   ifort /exe:test_ref *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
:IFORT
   ifort /c CFML_ssg_datafile.f90        /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC /warn
   ifort /c CFML_Rational_Groups.f90     /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC /warn
   ifort /c CFML_Standard_Sett.f90       /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC /warn
   ifort /c CFML_SuperSpaceGroups.f90    /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC /warn
   ifort /c test_ref.f90                 /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC /warn
   ifort /exe:test_ref *.obj C:\CrysFML\ifort\LibC\crysfml.lib
   goto END
:IFORT64
   ifort /c CFML_ssg_datafile.f90        /O3      /IC:\CrysFML\ifort64\LibC /warn
   ifort /c CFML_Rational_Groups.f90     /O3      /IC:\CrysFML\ifort64\LibC /warn
   ifort /c CFML_Standard_Sett.f90       /O3      /IC:\CrysFML\ifort64\LibC /warn
   ifort /c CFML_SuperSpaceGroups.f90    /O3      /IC:\CrysFML\ifort64\LibC /warn
   ifort /c test_ref.f90                 /O3      /IC:\CrysFML\ifort64\LibC /warn
   ifort /exe:test_ref *.obj C:\CrysFML\ifort64\LibC\crysfml.lib /link /stack:128000000
   goto END
:IFORT64D
   ifort /c CFML_ssg_datafile.f90        /heap-arrays /debug=full /traceback /nologo   /IC:\CrysFML\ifort64_debug\LibC /warn
   ifort /c CFML_Rational_Groups.f90     /heap-arrays /debug=full /traceback /nologo   /IC:\CrysFML\ifort64_debug\LibC /warn
   ifort /c CFML_Standard_Sett.f90       /heap-arrays /debug=full /traceback /nologo   /IC:\CrysFML\ifort64_debug\LibC /warn
   ifort /c CFML_SuperSpaceGroups.f90    /heap-arrays /debug=full /traceback /nologo   /IC:\CrysFML\ifort64_debug\LibC /warn
   ifort /c test_ref.f90                 /heap-arrays /debug=full /traceback /nologo   /IC:\CrysFML\ifort64_debug\LibC /warn
   ifort /exe:test_ref *.obj C:\CrysFML\ifort64_debug\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c CFML_ssg_datafile.f90     -ffree-line-length-none  -I../../GFortran/LibC
   gfortran -c CFML_Rational_Groups.f90  -ffree-line-length-none  -I../../GFortran/LibC
   gfortran -c CFML_Standard_Sett.f90    -ffree-line-length-none  -I../../GFortran/LibC
   gfortran -c CFML_SuperSpaceGroups.f90 -ffree-line-length-none  -I../../GFortran/LibC
   gfortran -c test_ref.f90              -ffree-line-length-none  -I../../GFortran/LibC
   gfortran *.o -o test_ref_gf -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
