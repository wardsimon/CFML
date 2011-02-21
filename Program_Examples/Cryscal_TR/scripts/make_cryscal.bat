@echo off
rem
rem ---- CRYSCAL compilation ----
rem
   if not x%1 == x goto CONT
rem
rem Info Header
rem
   cls
   echo MAKE_CRYSCAL: Make the CRYSCAL
   echo Syntax: make_cryscal [lf95/g95/all]  [debug]
   echo         --- compiler ---        ---  options  ---
   goto FIN
rem
rem Compilation
rem
:CONT
   echo **-------------------------------------------------------**
   echo **----                                               ----**
   echo **---- CRYSCAL                                       ----**
   echo **----                                               ----**
   echo **---- TR                                 (2007-2011)----**
   echo **-------------------------------------------------------**
rem
rem Compiler Version
rem
   if x%1 == xall      goto LF95_ZONE
   if x%1 == xlf95     goto LF95_ZONE
   if x%1 == xg95      goto G95_ZONE
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
   goto FIN
rem
rem >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem
:FIN
