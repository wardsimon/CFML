@echo off
rem ****
rem ****---- Compilation for ILL_Pubs_WoS Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: May-2014
rem ****
rem
   if not x%1 == x goto CONT
   cls
   echo    MAKE_PUBS_WoS: Make ILL_Pubs_WoS Compilation
   echo    Syntax: make_Pubs_WoS [gfortran/ifort]
   goto END
rem
:CONT
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   goto END
rem
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c CFML.f90                /nologo
   ifort /c Data_Articles_Mod.f90   /nologo
   ifort /c ILL_Pubs_WoS.f90        /nologo
   ifort /exe:ILL_Pubs_WoS *.obj  /link /stack:64000000
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3 -funroll-loops  -msse2   CFML.f90
   gfortran -c -O3 -funroll-loops  -msse2   Data_Articles_Mod.f90
   gfortran -c -O3 -funroll-loops  -msse2   ILL_Pubs_WoS.f90
   gfortran  *.o -o ILL_Pubs_WoS -O3  -funroll-loops  -msse2
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
