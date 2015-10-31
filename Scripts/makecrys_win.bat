@echo off
rem --------------------------------------------------
rem ---- Crystallographic Fortran Modules Library ----
rem ----               JGP/JRC                    ----
rem ----        September - 2015                  ----
rem --------------------------------------------------
   if not [%1]==[] (goto CONT)
rem ---- Default Message ----
:INIT
   cls
   echo MAKECRYS: Making the CrysFML Library
   echo Syntax: makecrys [gfortran|ifort[15]|lf95] [all, winter|realwin] [debug]
   goto FINAL
rem ---- Variables ----
:CONT
   (set _COMPILER=)
   (set _DEBUG=N)
   (set _WINTER=N)
   (set _CONSOLE=Y)
   (set _REALWIN=N)
rem
rem ---- Arguments ----
:LOOP
    if [%1]==[gfortran] (set _COMPILER=gfortran)
    if [%1]==[ifort]    (set _COMPILER=ifort)
    if [%1]==[ifort15]  (set _COMPILER=ifort15)
    if [%1]==[lf95]     (set _COMPILER=lf95)
    if [%1]==[debug]    (set _DEBUG=Y)
    if [%1]==[winter] (
       (set _WINTER=Y)
       (set _CONSOLE=N)
    )
    if [%1]==[realwin] (
       (set _REALWIN=Y)
       (set _CONSOLE=N)
    )
    if [%1]==[all] (
       (set _CONSOLE=Y)
       (set _WINTER=Y)
       if [%_COMPILER%]==[lf95] (
          (set _WINTER=N)
          (set _REALWIN=Y)
       )
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

   if [%_COMPILER%]==[ifort15] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_ifort 15) else (call compile_ifort 15 debug)
      )
      if [%_WINTER%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_ifort winter 15) else (call compile_ifort winter 15 debug)
      )
      if [%_REALWIN%]==[Y] (
         echo Sorry, option not compatible!!!
         goto END
      )
      goto END
   )
rem
rem ----
rem ---- Lahey
rem ----
   if [%_COMPILER%]==[lf95] (
      if [%_CONSOLE%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_lahey) else (call compile_lahey debug)
      )
      if [%_WINTER%]==[Y] (
         echo Sorry, option not compatible!!!
         goto END
      )
      if [%_REALWIN%]==[Y] (
         if [%_DEBUG%]==[N] (call compile_lahey "realwin") else (call compile_lahey "realwin" "debug")
      )
      goto END
   )
rem
:END
   cd ..
rem
:FINAL
