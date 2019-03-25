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
   echo Syntax: "makecrys [gfortran|ifort] [all, winter] [debug]"
   goto FINAL
rem ---- Variables ----
:CONT
   (set _COMPILER=)
   (set _DEBUG=N)
   (set _WINTER=N)
   (set _CONSOLE=Y)
rem
rem ---- Arguments ----
:LOOP
    if [%1]==[gfortran] (set _COMPILER=gfortran)
    if [%1]==[ifort]    (set _COMPILER=ifort)
    if [%1]==[debug]    (set _DEBUG=Y)
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
rem ---- Windows Subdirectory
   cd .\windows
rem
rem -------------------
rem ---- Compilers ----
rem -------------------
rem
rem ----
rem ---- GFortran
rem ----
   if [%_COMPILER%]==[gfortran] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_gfortran) else (call compile_gfortran debug)
      )
      if [%_WINTER%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_gfortran winter) else (call compile_gfortran winter debug)
      )
      if [%_REALWIN%]==[Y] (
         echo Sorry, option not compatible!!!
         goto END
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
      if [%_REALWIN%]==[Y] (
         echo Sorry, option not compatible!!!
         goto END
      )
      goto END
   )
rem
:END
   cd ..
rem
:FINAL
