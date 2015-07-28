@echo off
rem --------------------------------------------------
rem ---- Crystallographic Fortran Modules Library ----
rem ----               JGP/JRC                    ----
rem ----             June - 2015                  ----
rem --------------------------------------------------
   if not [%1]==[] goto CONT
rem ---- Default Message ----
:INIT
   cls
   echo MAKECRYS: Making the CrysFML Library
   echo Syntax: makecrys [f95/lf95/g95/ifort/ifort15/gfortran] [con,win,rw] [debug]
   goto FINAL
rem ---- Variables ----
:CONT
   (set _COMPILER=)
   (set _DEBUG=N)
   (set _WINTER=N)
   (set _CONSOLE=N)
   (set _REALWIN=N)
rem
rem ---- Arguments ----
:LOOP
    if [%1]==[debug] (set _DEBUG=Y)
    if [%1]==[win] (set _WINTER=Y)
    if [%1]==[con] (set _CONSOLE=Y)
    if [%1]==[rw] (set _REALWIN=Y)
    if [%1]==[f95] (set _COMPILER=f95)
    if [%1]==[lf95] (set _COMPILER=lf95)
    if [%1]==[g95] (set _COMPILER=g95)
    if [%1]==[ifort] (set _COMPILER=ifort)
    if [%1]==[gfortran] (set _COMPILER=gfortran)
    if [%1]==[ifort15] (set _COMPILER=ifort15)
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
rem ---- Absoft
rem ----
   if [%_COMPILER%]==[f95] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[Y] (call compile_ifort debug) else (call compile_ifort)
      )
      if [%_WINTER%]==[Y] (
         if [%_DEBUG%]==[Y] (call compile_ifort winter debug) else (call compile_ifort winter)
      )
      if [%_REALWIN%]==[Y] (
         echo Sorry, option not compatible!!!
         goto END
      )
      goto END
   )
rem ----
rem ---- Lahey
rem ----
   if [%_COMPILER%]==[lf95] (
      echo Compilador Lahey
      goto END
   )
rem ----
rem ---- G95
rem ----
   if [%_COMPILER%]==[g95] (
      echo Compilador Lahey
      goto END
   )
rem ---- GFortran
   if [%_COMPILER%]==[gfortran] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[Y] (call compile_gfortran debug) else (call compile_gfortran)
      )
      if [%_WINTER%]==[Y] (
         if [%_DEBUG%]==[Y] (call compile_gfortran winter debug) else (call compile_gfortran winter)
      )
      if [%_REALWIN%]==[Y] (
         echo Sorry, option not compatible!!!
         goto END
      )
      goto END
   )
rem ----
rem ---- Intel
rem ----
   if [%_COMPILER%]==[ifort] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[Y] (call compile_ifort debug) else (call compile_ifort)
      )
      if [%_WINTER%]==[Y] (
         if [%_DEBUG%]==[Y] (call compile_ifort winter debug) else (call compile_ifort winter)
      )
      if [%_REALWIN%]==[Y] (
         echo Sorry, option not compatible!!!
         goto END
      )
   )

   if [%_COMPILER%]==[ifort15] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[Y] (call compile_ifort debug 15) else (call compile_ifort 15)
      )
      if [%_WINTER%]==[Y] (
         if [%_DEBUG%]==[Y] (call compile_ifort winter debug 15) else (call compile_ifort winter 15)
      )
      if [%_REALWIN%]==[Y] (
         echo Sorry, option not compatible!!!
         goto END
      )
   )
rem
:END
   cd ..
rem
:FINAL
