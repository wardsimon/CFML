@echo off
rem --------------------------------------------------
rem ---- Crystallographic Fortran Modules Library ----
rem ---- JGP/JRC                             2019 ----
rem --------------------------------------------------
    if not [%1]==[] (goto CONT)
rem ---- Default Message ----
:INIT
    cls
    echo MAKECRYS: Making the CrysFML Library
    echo Syntax: "makecrys [gfortran|ifort[32|64]] [all, winter] [debug]"
    goto FINAL
rem ---- Variables ----
:CONT
    (set _COMPILER=)
    (set _DEBUG=N)
    (set _WINTER=N)
    (set _CONSOLE=Y)
    (set _MODE=32)
rem
:LOOP
rem ---- Compiler ----
    if [%1]==[gfortran] (set _COMPILER=gfortran)
    if [%1]==[gfortran32] (
       (set _MODE=32)
       (set _COMPILER=gfortran)
    )
    if [%1]==[gfortran64] (
       (set _MODE=64)
       (set _COMPILER=gfortran)
    )
    if [%1]==[ifort]    (set _COMPILER=ifort)
    if [%1]==[ifort32] (
       (set _MODE=32)
       (set _COMPILER=ifort)
    )
    if [%1]==[ifort64] (
       (set _MODE=64)
       (set _COMPILER=ifort)
    )
rem    
rem ---- Debug ----    
    if [%1]==[debug]    (set _DEBUG=Y)
rem
rem ---- Console / GUI ----    
    if [%1]==[winter] (
       (set _WINTER=Y)
       (set _CONSOLE=N)
    )
    if [%1]==[all] (
       (set _CONSOLE=Y)
       (set _WINTER=Y)
    )
    shift
    if not [%1]==[] goto LOOP
rem
rem
    echo.
    echo **---------------------------------------------**
    echo **---- CrysFML: Final compilation options  ----**
    echo **---------------------------------------------**
    echo COMPILER:%_COMPILER%
    echo    DEBUG:%_DEBUG%
    echo  CONSOLE:%_CONSOLE%
    echo   WINTER:%_WINTER%
    echo.
rem ---- Windows Subdirectory
    cd .\windows
rem
rem ----
rem ---- GFortran
rem ----
   if [%_COMPILER%]==[gfortran] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_gfortran %_MODE% ) else (call compile_gfortran %_MODE% debug)
      )
      if [%_WINTER%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_gfortran %_MODE% winter) else (call compile_gfortran %_MODE% winter debug)
      )
      goto END
   )
rem
rem ----
rem ---- Intel
rem ----
   if [%_COMPILER%]==[ifort] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_ifort) else (call compile_ifort debug)
      )
      if [%_WINTER%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_ifort winter) else (call compile_ifort winter debug)
      )
      goto END
   )
rem
:END
   cd ..
rem
:FINAL
