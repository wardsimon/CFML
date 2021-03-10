@echo off
rem ****
rem **** Windows Makefile for hkl_gen Program ****
rem ****
rem **** Author: JRC (based in make_racer by EP)
rem **** Revision: July-2010
rem ****
rem The first argument is the compiler name
set COMPILER=%~1

rem No arguments - print the help
if x%COMPILER% == x goto EMPTYARGS

rem Check that the CRYSFML env variables has been previously set.
if x"%CRYSFML%" == x goto CRYSFML_PATH_UNSET

rem Check that the corresponding path really exists.
if not exist "%CRYSFML%" goto WRONG_CRYSFML_PATH

rem List of the file to compile. The last one is the program.
set SRC=(hkl_gen.f90)
set OBJ=hkl_gen.obj

if x%COMPILER% == xg95      goto G95_ZONE
if x%COMPILER% == xgfortran goto GFORTRAN_ZONE
if x%COMPILER% == xifort    goto IFORT_ZONE

rem **** None of the supported compiler was provided. ****
goto WRONGCOMPILER

rem **----------------------**
rem **---- G95 Compiler ----**
rem **----------------------**
:G95_ZONE
   set CRYSFMLOBJECTS="%CRYSFML%"\G95\LibC
   rem **** Check that CrySFML has been compiled with g95 compiler ****
   if not EXIST "%CRYSFMLOBJECTS%" goto NOTCOMPILED



   if x%2 == xdebug goto G95_CD

   goto G95_CO

:G95_CO
   set COMP_ARGS=-c -Wall -pedantic -fbounds-check -funroll-loops -msse2 -I"%CRYSFMLOBJECTS%"
   set EXEC=hkl_gen_%COMPILER%.exe
   set LINK_ARGS=%OBJ% -o %EXEC% -O -L"%CRYSFMLOBJECTS%" -lcrysfml

   goto COMPILE

:G95_CD
   set COMP_ARGS=-c -Wall -pedantic -fbounds-check -funroll-loops -msse2 -ftrace=full -I"%CRYSFMLOBJECTS%"
   set EXEC=hkl_gen_%COMPILER%_debug.exe
   set LINKERARGS=%OBJ% -o %EXEC% -O -L"%CRYSFMLOBJECTS%" -lcrysfml

   goto COMPILE

rem **---------------------------**
rem **---- GFortran Compiler ----**
rem **---------------------------**
:GFORTRAN_ZONE
   set CRYSFMLOBJECTS="%CRYSFML%"\GFortran\LibC
   rem **** Check that CrySFML has been compiled with g95 compiler ****
   if not EXIST "%CRYSFMLOBJECTS%" goto NOTCOMPILED
   if x%2 == xdebug goto GFORTRAN_CD
   goto GFORTRAN_CO

:GFORTRAN_CO
   set COMP_ARGS=-c -W -Wall -pedantic-errors -funroll-loops -msse2 -I"%CRYSFMLOBJECTS%"
   set EXEC=hkl_gen_%COMPILER%.exe
   set LINK_ARGS=*.o -o %EXEC% -O -L"%CRYSFMLOBJECTS%" -lcrysfml

   goto COMPILE


:GFORTRAN_CD
   set COMP_ARGS=-c -W -Wall -pedantic-errors -funroll-loops -msse2 -fbacktrace -I"%CRYSFMLOBJECTS%"
   set EXEC=hkl_gen_%COMPILER%_debug.exe
   set LINK_ARGS=%OBJ% -o %EXEC% -O -L"%CRYSFMLOBJECTS%" -lcrysfml

   goto COMPILE

rem **------------------------**
rem **---- Intel Compiler ----**
rem **------------------------**
:IFORT_ZONE

   if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
   if [%TARGET_ARCH%]==[ia32]  (
        if x%2 == xdebug (
          set  CRYSFMLOBJECTS="%CRYSFML%"\ifort_debug\libC
        ) else (
          set  CRYSFMLOBJECTS="%CRYSFML%"\ifort\libC
        )
     ) else (
        if x%2 == xdebug (
          set  CRYSFMLOBJECTS="%CRYSFML%"\ifort64_debug\libC
        ) else (
          set  CRYSFMLOBJECTS="%CRYSFML%"\ifort64\libC
        )
     )
   rem **** Check that CrySFML has been compiled with intel compiler ****
   if not EXIST "%CRYSFMLOBJECTS%" goto NOTCOMPILED
   if x%2 == xdebug goto IFORT_CD
   goto IFORT_CO

:IFORT_CO
   set COMP_ARGS=/c /O2 /nologo /I"%CRYSFMLOBJECTS%"
   set EXEC=hkl_gen_%COMPILER%.exe
   set LINK_ARGS=/exe:%EXEC% %OBJ% "%CRYSFMLOBJECTS%"\crysfml.lib /link /stack:64000000

   goto COMPILE

:IFORT_CD
   set COMP_ARGS=/c /debug:full /check /check:noarg_temp_created /traceback /nologo /I"%CRYSFMLOBJECTS%"
   set EXEC=hkl_gen_%COMPILER%_debug.exe
   set LINK_ARGS=/exe:%EXEC% %OBJ% "%CRYSFMLOBJECTS%"\crysfml.lib /link /stack:64000000

   goto COMPILE

:NOTCOMPILED
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo ERROR: CrysFML not compiled
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo.
        goto SYNTAX

:CRYSFML_PATH_UNSET
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo ERROR: The path to CRYSFML must be set as
        echo ERROR: an environment variable
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo.
        goto SYNTAX

:WRONG_CRYSFML_PATH
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo ERROR: The path to CRYSFML
        echo ERROR:     %CRYSFML%
        echo ERROR: does not exist
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo.
        goto SYNTAX

:WRONGCOMPILER
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo ERROR: %COMPILER% - unknown compiler
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo.
        goto SYNTAX

:EMPTYARGS
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo ERROR: No arguments provided
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo.
        goto SYNTAX

:WRONGARGS
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo ERROR: You must provide two arguments
        echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        echo.
        goto SYNTAX

:SYNTAX
        echo MAKE_pfind: Make Pfind Compilation
        echo Syntax: make_pfind [g95/gfortran/ifort] [debug]
        goto END

:COMPILE
   for %%F in %SRC% do (
      echo Compiling module %%F ...
      echo %COMPILER% %COMP_ARGS% %%F
      %COMPILER% %COMP_ARGS% %%F
      echo Done
      echo.
   )

   echo Linking program %EXEC% ...
   echo %COMPILER% %LINK_ARGS%
   %COMPILER% %LINK_ARGS%
   echo Done
   echo.

   rem Compress executable
   upx %EXEC%

   rem Update SXT Distribution
   copy %EXEC% %SXTALSOFT%\DistSXT\Windows\bin\hkl_gen.exe


   if x"%FULLPROF%" == x goto END

   if exist "%FULLPROF%" copy %EXEC% %FULLPROF%\hkl_gen.exe

rem Remove object, exe, map, bak and mod files
rem   set TODELETE=%OBJ% %MOD%
rem   del %TODELETE%
   del *.obj *.mod *.map *.bak *.exe

   goto END

:END
