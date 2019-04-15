@echo off
rem call "C:\Program Files (x86)\IntelSWTools\parallel_studio_xe_2016.0.041\bin\"psxevars intel64
rem call "C:\Program Files (x86)\IntelSWTools\parallel_studio_xe_2017.0.036\bin\"psxevars intel64
call "C:\Program Files (x86)\IntelSWTools\parallel_studio_xe_2019.3.056\bin\"psxevars intel64
SET FULLPROF=c:\FullProf_Suite
SET SXTALSOFT=c:\SXtalSoft
SET CRYSFML=c:\CrysFML
SET PROGCFML=c:\ProgCFML
SET LIB=%LIB%;c:\CrysFML\ifort64\LibC08;c:\CrysFML\ifort64\LibW08;c:\wint\lib.i64
SET LAUESUITE=c:\Esmeralda_Laue_Suite
SET PATH=%PATH%;c:\FullProf_Suite;c:\Esmeralda_Laue_Suite;C:\Program Files\CMake\bin
SET PATH=%PATH%;C:\Program Files (x86)\Notepad++
echo "This terminal is adequate for running FullProf Suite and Cmake from the command line"
