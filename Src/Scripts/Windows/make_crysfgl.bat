@echo off
rem
rem ---- CRYSFGL LIBRARY 4.00 ----
rem
   if not x%1 == x goto CONT
rem
rem Info Header
rem
   cls
   echo MAKE_CRYSFGL: Make the CrysFGL Library
   echo Syntax: make_crysfml [lf95/ifort] [debug]
   goto FIN
rem
rem Compilation
rem
:CONT
   echo **--------------------------------------------------------**
   echo **----                                                ----**
   echo **---- CRYSTALLOGRAPHIC FORTRAN GRAPHICS LIBRARY 5.00 ----**
   echo **----                                                ----**
   echo **---- Authors: JRC-LC                     (1999-2011)----**
   echo **----                                                ----**
   echo **--------------------------------------------------------**
rem
rem Compiler Version
rem
   if x%1 == xlf95  goto LF95_ZONE
   if x%1 == xifort goto IFORT_ZONE
   goto FIN
rem
rem ------------------------
rem ---- LAHEY COMPILER ----
rem ------------------------
:LF95_ZONE
   if x%2 == xdeb goto LF95_D
   call comp_lf95_gl_o
   goto FIN
rem
rem
rem
:LF95_D
   call comp_lf95_gl_d
   goto FIN
rem
rem ------------------------
rem ---- INTEL COMPILER ----
rem ------------------------
:IFORT_ZONE
   if x%2 == xdeb goto IFORT_D
   call comp_ifort_gl_o
   goto FIN
rem
rem
rem
:IFORT_D
   call comp_ifort_gl_d
   goto FIN
rem
:FIN

