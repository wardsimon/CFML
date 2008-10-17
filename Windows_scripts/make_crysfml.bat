@echo off
if not x%1 == x goto CONT
rem
cls
cd .\src\cfml
echo
echo MAKE_CRYSFML: Make the CrysFML Library
echo Syntax: make_crysfml [f95/lf95/g95/ifort/gfortran] [win/rwin] [debug]
echo                            --- compiler ---        ---  options  ---
goto FIN
rem
:CONT
cd .\src\cfml
rem
echo **-------------------------------------------------------**
echo **----                                               ----**
echo **---- CRYSTALLOGRAPHIC FORTRAN MODULES LIBRARY 3.00 ----**
echo **----                                               ----**
echo **---- Authors: JRC-JGP                   (1999-2008)----**
echo **-------------------------------------------------------**
rem
if x%1 == xf95      goto F95_ZONE
if x%1 == xlf95     goto LF95_ZONE
if x%1 == xifort    goto IFORT_ZONE
if x%1 == xg95      goto G95_ZONE
if x%1 == xgfortran goto GFOR_ZONE
goto FIN
rem
rem ********************
rem ** ABSOFT COMPILER **
rem ********************
rem
:F95_ZONE
if x%2 == xdeb goto F95_CD
if x%2 == xwin goto F95_WINTER
rem
rem CONSOLE + OPTIMIZATION
rem
   call comp_f95_co
   goto FIN
rem
rem CONSOLE + DEBUG
rem
:F95_CD
   call comp_f95_cd
   goto FIN
rem
:F95_WINTER
if x%3 == xdeb goto F95_WD
rem
rem WINTERACTER + OPTIMIZATION
rem
   call comp_f95_wo
   goto FIN
rem
rem WINTERACTER + DEBUG
rem
:F95_WD
   call comp_f95_wd
   goto FIN
rem
rem ********************
rem ** LAHEY COMPILER **
rem ********************
rem
:LF95_ZONE
if x%2 == xdeb goto LF95_CD
if x%2 == xwin goto LF95_WINTER
if x%2 == xrwin goto LF95_REALWIN
rem
rem CONSOLE + OPTIMIZATION
rem
   call comp_lf95_co
   goto FIN
rem
rem CONSOLE + DEBUG
rem
:LF95_CD
   call comp_lf95_cd
   goto FIN
rem
:LF95_WINTER
if x%3 == xdeb goto LF95_WD
rem
rem WINTERACTER + OPTIMIZATION
rem
   call comp_lf95_wo
   goto FIN
rem
rem WINTERACTER + DEBUG
rem
:LF95_WD
   call comp_lf95_wd
   goto FIN
rem
:LF95_REALWIN
if x%3 == xdeb goto LF95_RD
rem
rem REALWIN + OPTIMIZATION
rem
   call comp_lf95_ro
   goto FIN
rem
rem REALWIN + DEBUG
rem
:LF95_RD
   call comp_lf95_rd
   goto FIN
rem
rem
rem ********************
rem ** INTEL COMPILER **
rem ********************
rem
:IFORT_ZONE
if x%2 == xdeb goto IFORT_CD
if x%2 == xwin goto IFORT_WINTER
rem
rem CONSOLE + OPTIMIZATION
rem
   call comp_ifort_co
   goto FIN
rem
rem CONSOLE + DEBUG
rem
:IFORT_CD
   call comp_ifort_cd
   goto FIN
rem
:IFORT_WINTER
if x%3 == xdeb goto IFORT_WD
rem
rem WINTERACTER + OPTIMIZATION
rem
   call comp_ifort_wo
   goto FIN
rem
rem WINTERACTER + DEBUG
rem
:IFORT_WD
   rem call comp_ifort_wd
   goto FIN
rem
rem
rem ******************
rem ** G95 COMPILER **
rem ******************
rem
:G95_ZONE
if x%2 == xdeb goto G95_CD
rem
rem CONSOLE + OPTIMIZATION
rem
   call comp_g95_co
   goto FIN
rem
rem CONSOLE + DEBUG
rem
:G95_CD
   call comp_g95_cd
   goto FIN
rem
rem
rem ***********************
rem ** GFORTRAN COMPILER **
rem ***********************
rem
:GFOR_ZONE
if x%2 == xdeb goto GFOR_CD
rem
rem CONSOLE + OPTIMIZATION
rem
   call comp_gfor_co
   goto FIN
rem
rem CONSOLE + DEBUG
rem
:GFOR_CD
   call comp_gfor_cd
   goto FIN
rem
rem >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem
:FIN
cd ..\..
