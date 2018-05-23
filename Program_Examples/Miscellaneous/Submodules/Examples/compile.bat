@echo off
rem
rem Intel Compilation
rem
    echo .... compilando modulo principal
    ifort Module_Test1.f90     /c /nologo
    echo.
    echo.
    echo .... compilando submodulos
    ifort Submodule01.f90      /c /nologo
    ifort Submodule02.f90      /c /nologo
    ifort Submodule03.f90      /c /nologo
rem
    echo.
    echo.
    echo .... creando el ejecutable!!!!
    ifort main.f90      /c /nologo
    ifort /exe:Ejemplo *.obj 
