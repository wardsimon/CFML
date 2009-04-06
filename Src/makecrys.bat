@echo off
rem
rem Crystallographic Fortran Modules Library
rem Author: JGP
rem Revision: Dic-2009
rem
   if not x%1 == x goto CONT
rem
   cls
   echo MAKECRYS: Making the CrysFML and CrysFGL Libraries
   echo Syntax: makecrys [f95/lf95/g95/ifort/gfortran] [con,all]
   goto FINAL
rem
:CONT
   cd .\scripts\windows
   call make_crysfml %1
   if x%2 == xall goto ALL
   goto END
:ALL
   if x%1 == xlf95 call make_crysfml %1 rwin
   call make_crysfml %1 win
   call make_crysfgl %1
:END
   cd ..\..
rem
:FINAL
