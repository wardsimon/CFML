@echo off
rem
rem Crystallographic Fortran Modules Library
rem Author: JGP
rem Revision: Nov-2008
rem
   if not x%1 == x goto CONT
rem
   cls
   echo MAKE_ALL: Making the CrysFML and CrysFGL Libraries
   echo Syntax: make_crysfml [f95/lf95/g95/ifort/gfortran] 
   goto FIN
rem
:CONT
   cd .\scripts\windows
   call make_crysfml %1 %2
   call make_crysfml %1 win %2
   call make_crysfgl %1 %2
   cd ..\..\..
rem
:FIN
