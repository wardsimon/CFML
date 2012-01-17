@echo off
rem
rem Crystallographic Fortran Modules Library
rem Author: JGP/JRC
rem Revision: Jan-2012
rem
   if not x%1 == x goto CONT
rem
   cls
   echo MAKECRYS: Making the CrysFML Library
   echo Syntax: makecrys [f95/lf95/g95/ifort/gfortran] [all,win,rw] [deb]
   echo Debug mode cannot be invoked with the "all" option
   echo The "deb" option has to be provided as second argument for debugging console programs
   goto FINAL
rem
:CONT
   cd .\Windows
   if x%2 == xwin goto WIN
   if x%2 == xrw goto RW
rem %2 may be "deb", for console debug in which case no other option is possible
   call make_crysfml %1 %2
   if x%2 == xall goto WIN
   goto END
:WIN
   if x%3 == xdeb goto WINDEB
   if x%1 == xlf95 call make_crysfml %1 rwin
   call make_crysfml %1 win
   goto END
:RW
   if x%3 == xdeb goto RWDEB
   if x%1 == xlf95 call make_crysfml %1 rwin
   goto END
:WINDEB
   if x%1 == xlf95 call make_crysfml %1 rwin deb
   call make_crysfml %1 win deb
   goto END
:RWDEB
   if x%1 == xlf95 call make_crysfml %1 rwin deb
   goto END
:END
   cd ..
rem
:FINAL
