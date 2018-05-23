@echo off
rem
rem GFortran Compilation (v7.3 for Windows)
rem
    echo .... compiling main module
    gfortran -c -m64 Module_Test1.f90 
    echo.
    echo.
    echo .... compiling submodules
    gfortran -c -m64 Submodule01.f90     
    gfortran -c -m64 Submodule02.f90 
    gfortran -c -m64 Submodule03.f90      
rem
    echo.
    echo.
    echo .... Making the executable!!!!
    gfortran -c -m64 main.f90
    gfortran -m64 -o Example.exe  *.o 
