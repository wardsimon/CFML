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
   echo Syntax: makecrys [f95/lf95/g95/ifort/gfortran] [con,win,rw] [deb]
   echo Debug mode cannot be invoked with the "all" option
   echo The "deb" option has to be provided as second argument for debugging console programs
   goto FINAL
rem ----
rem
:CONT
rem ---- Windows Subdirectory
   cd .\windows
rem ---- Options ----
   (set _debug=)
   if [%3]==[deb] (set _debug=y)
rem
   (set _opt1=)
   if [%2]==[deb] (
      set _debug=y
   ) else (
      if not [%2]==[] (set _opt1=%2)
   )
rem
rem ---- Compilers ----
rem
rem ---- Absoft
   if [%1]==[f95] (
      echo Compilador Lahey
      goto END
   )
rem ---- Lahey
   if [%1]==[lf95] (
      echo Compilador Lahey
      goto END
   )
rem ---- G95
   if [%1]==[g95] (
      echo Compilador Lahey
      goto END
   )
rem ---- Intel
   if [%1]==[ifort] (
      if [%_opt1%]==[con] (
         echo Intel + Console
         if not [%_debug%]==[Y] (call comp_ifort_co) else (call comp_ifort_cd)
      )
      if [%_opt1%]==[win] (
         echo Intel + Wintercater
         if not [%_debug%]==[Y] (call comp_ifort_wo) else (call comp_ifort_wd)
      )
      if [%_opt1%]==[rw] (
         echo Sorry, option not compatible!!!
         goto END
      )
      if [%_opt1%]==[] (
         echo Intel + Console + Wintercater
         if not [%_debug%]==[Y] (
            call comp_ifort_co
            call comp_ifort_wo
         ) else (
            call comp_ifort_cd
            call comp_ifort_wd
         )
      )
      goto END
   )
rem ---- GFortran
   if [%1]==[gfortran] (
      echo Compilador GFortran
      goto END
   )
rem
:END
   cd ..
rem
:FINAL
