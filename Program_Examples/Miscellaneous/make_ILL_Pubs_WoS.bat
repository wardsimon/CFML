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
   echo    MAKE_ILL_PUBS_WoS: Make ILL_Pubs_WoS Compilation
   echo    Syntax: ILL_Pubs_WoS [f95/lf95/g95/gfrotran/ifort]
   goto END
rem
:CONT
   if x%1 == xf95       goto F95
   if x%1 == xlf95      goto LF95
   if x%1 == xg95       goto G95
   if x%1 == xgfortran  goto GFOR
   if x%1 == xifort     goto IFORT
   goto END
rem
rem ****---- Absoft Compiler ----****
:F95
   goto END
rem
rem ****---- Lahey Compiler ----****
:LF95
   lf95 -c Data_Articles_Mod.f90  -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 -c ILL_Pubs_WoS.f90       -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
   lf95 *.obj -out ILL_Pubs_WoS -o1 -lib C:\crysFML\lahey\libC\crysFML
   goto END
rem
rem ****---- Intel Compiler ----****
:IFORT
   ifort /c Data_Articles_Mod.f90 /debug /nologo /IC:\CrysFML\Ifort\LibC
   ifort /c ILL_Pubs_WoS.f90      /debug /nologo /IC:\CrysFML\Ifort\LibC
   ifort /exe:ILL_Pubs_WoS *.obj C:\CrysFML\Ifort\LibC\crysfml.lib  /link /stack:64000000
   goto END
rem
rem **---- G95 Compiler ----**
:G95
   g95 -c -O3 -funroll-loops -ftrace=full -msse2  Data_Articles_Mod.f90  -IC:\CrysFML\G95\LibC
   g95 -c -O3 -funroll-loops -ftrace=full -msse2   ILL_Pubs_WoS.f90      -IC:\CrysFML\G95\LibC
   g95  *.o -o ILL_Pubs_WoS -O3 -ftrace=full -funroll-loops  -msse2      -LC:\CrysFML\G95\LibC -lcrysfml
   goto END
rem
rem **---- GFORTRAN Compiler ----**
:GFOR
   gfortran -c -O3 -funroll-loops  -msse2  Data_Articles_Mod.f90  -IC:\CrysFML\GFortran\LibC
   gfortran -c -O3 -funroll-loops  -msse2   ILL_Pubs_WoS.f90      -IC:\CrysFML\GFortran\LibC
   gfortran  *.o -o ILL_Pubs_WoS -O3  -funroll-loops  -msse2      -LC:\CrysFML\GFortran\LibC -lcrysfml
   goto END
rem
:END
   del *.obj *.mod *.o *.map *.bak > nul
