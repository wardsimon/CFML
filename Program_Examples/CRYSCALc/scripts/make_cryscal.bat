@echo off
rem
rem ---- CRYSCALC compilation ----
rem
   if not x%1 == x goto CONT
rem
rem Info Header
rem
   cls
   echo MAKE_CRYSCALC: Make the CRYSCALC
   echo Syntax: make_cryscalc [lf95/g95/ifort/all]  [debug]
   echo         --- compiler ---        ---  options  ---
   goto FIN
rem
rem Compilation
rem
:CONT
   echo **--------------------------------------------------------**
   echo **----                                                ----**
   echo **----                   CRYSCALC                     ----**
   echo **----                                                ----**
   echo **----          TR [CDIFX-UMR6226 Rennes]             ----**
   echo **----                (2007-2014)                     ----**
   echo **----                                                ----**
   echo **----  Compilation:                                  ----**
   echo **         make_cryscalc [lf95/g95/ifort/all]         ----**
   echo **--------------------------------------------------------**
rem
rem Compiler Version
rem
   if x%1 == xall      goto LF95_ZONE
   if x%1 == xlf95     goto LF95_ZONE
   if x%1 == xg95      goto G95_ZONE
   if x%1 == xifort    goto IFORT_ZONE
   goto FIN
rem
rem ------------------------
rem ---- LAHEY COMPILER ----
rem ------------------------
:LF95_ZONE
rem
rem CONSOLE + OPTIMIZATION
rem
   call make_comp_lf95
   call make_comp_lf95_w
   if x%1 == xlf95 goto FIN
rem
rem CONSOLE + DEBUG
rem
rem
rem ----------------------
rem ---- G95 COMPILER ----
rem ----------------------
rem
:G95_ZONE
rem
rem CONSOLE
rem
   call make_comp_g95
   if x%1 == xg95 goto FIN
rem ----------------------------------------------
rem ---- Intel Fortran COMPILER ----
rem ----------------------------------------------
rem
:IFORT_ZONE
rem
rem CONSOLE
rem
   call make_comp_ifort
   goto FIN
rem
rem >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem
:FIN
