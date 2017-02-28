@echo off
rem
rem ---- CRYSCALC compilation ----
rem
   if not x%1 == x goto CONT
rem
rem Info Header
rem
   cls
   echo MAKE_CRYSCALC64: Make the CRYSCALC for 64 bits
   echo Syntax: make_cryscalc64 [lf95/g95/ifort] [deb]
   echo                          --- compiler ---
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
   echo **         make_cryscalc64 [lf95/g95/ifort/all]       ----**
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
   if x%2 == xdeb got ifort_deb
   call make_comp_ifort64
   goto FIN
:ifort_deb
   call make_comp_ifort64d
   goto FIN
rem
rem >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem
:FIN
rem
rem UPX Compression  (all that below is already done in make_comp_ifort)
rem
rem   upx --compress-icons=0 cryscalc.exe
rem   copy cryscalc.exe ..\..\DistFPS_64b
rem   if exist %FULLPROF% copy cryscalc.exe %FULLPROF%
rem
rem   del *.obj *.mod *.bak *.map *.exe
rem
rem   cd ..\Scripts\Windows
rem

