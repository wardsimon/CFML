@echo off
rem
rem
echo **---- CrysFGL Library: Winteracter version ----**
rem
rem
if exist ..\..\Intel\LibGL rmdir ..\..\Intel\LibGL /S /Q
mkdir ..\..\Intel\LibGL
rem
rem
copy *.mod ..\..\Intel\LibGL > nul
move *.lib ..\..\Intel\LibGL > nul
del *.obj *.mod *.bak *.lst > nul
rem