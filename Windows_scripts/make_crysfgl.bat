@echo off
if not x%1 == x goto CONT
rem 
cls
cd .\src\cfgl
echo
echo MAKE_CRYSFGL: Make the CrysFGL Library
echo Syntax: make_crysfml [lf95/ifort] [debug]
echo                      --- compiler ---        -- option --
goto FIN
rem
:CONT
cd .\src\cfgl
rem
echo **--------------------------------------------------------**
echo **----                                                ----**
echo **---- CRYSTALLOGRAPHIC FORTRAN GRAPHICS LIBRARY 3.00 ----**
echo **----                                                ----**
echo **---- Authors: JRC-LC                     (1999-2008)----**
echo **----                                                ----**
echo **--------------------------------------------------------**
rem
if x%1 == xlf95  goto LF95_ZONE
if x%1 == xifort goto IFORT_ZONE
goto FIN
rem 
rem ********************
rem ** LAHEY COMPILER **
rem ********************
rem
:LF95_ZONE
   if x%2 == xdeb goto LF95_D
   call comp_lf95_gl_o
   goto FIN
rem
:LF95_D
   call comp_lf95_gl_d
   goto FIN
rem   
rem ********************
rem ** INTEL COMPILER **
rem ********************
rem
:IFORT_ZONE
   if x%2 == xdeb goto IFORT_D
   call comp_ifort_gl_o
   goto FIN
rem
:IFORT_D
   call comp_ifort_gl_d
   goto FIN
rem
:FIN
cd ..\..

