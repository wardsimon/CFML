@echo on
rem ****
rem ****---- Compilation for Testin_SSG Program ----****
rem ****
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   if x%1 == xifortd    goto IFORTD
   goto END
rem
rem
rem ****---- Intel Compiler ----****
:IFORTD
   ifort /c CFML_ssg_datafile.f90             /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c CFML_Rational_Groups.f90          /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c CFML_SuperSpaceGroups.f90         /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /c test_ref.f90                      /heap-arrays /debug=full /traceback /nologo /IC:\CrysFML\ifort_debug\LibC
   ifort /exe:test_ref *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
:IFORT
   ifort /c CFML_ssg_datafile.f90        /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC
   ifort /c CFML_Rational_Groups.f90     /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC
   ifort /c CFML_SuperSpaceGroups.f90    /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC
   ifort /c test_ref.f90                 /O2     /heap-arrays  /IC:\CrysFML\ifort\LibC
   ifort /exe:test_ref *.obj C:\CrysFML\ifort_debug\LibC\crysfml.lib
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c CFML_ssg_datafile.f90     -ffree-line-length-none  -I../../GFortran/LibC
   gfortran -c CFML_Rational_Groups.f90  -ffree-line-length-none  -I../../GFortran/LibC
   gfortran -c CFML_SuperSpaceGroups.f90 -ffree-line-length-none  -I../../GFortran/LibC
   gfortran -c test_ref.f90              -ffree-line-length-none  -I../../GFortran/LibC
   gfortran *.o -o test_ref_gf -L../../GFortran/LibC   -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
