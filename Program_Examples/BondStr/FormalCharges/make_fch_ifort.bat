@echo off
rem
rem Intel Compilation
rem
    if [%TARGET_ARCH%]==[] (set TARGET_ARCH=ia32)
    if [%TARGET_ARCH%]==[ia32]  (
       set DISTRIB=%PROGCFML%\DistFPS
       if [x%1 == xdeb] (
         set  INC=/I"%CRYSFML%"\ifort_debug\libC
         set CRYSLIB="%CRYSFML%"\ifort_debug\libC\crysfml.lib
       ) else (
         set  INC=/I"%CRYSFML%"\ifort\libC
         set CRYSLIB="%CRYSFML%"\ifort\libC\crysfml.lib
       )
    ) else (
       set DISTRIB=%PROGCFML%\DistFPS_64b
       if [x%1 == xdeb] (
         set INC=/I"%CRYSFML%"\ifort64_debug\LibC
         set CRYSLIB="%CRYSFML%"\ifort64_debug\LibC\crysfml.lib
       ) else (
         set INC=/I"%CRYSFML%"\ifort64\LibC
         set CRYSLIB="%CRYSFML%"\ifort64\LibC\crysfml.lib
       )
    )
if x%1 == xdeb goto DEB
  ifort Formal_Charges.f90 /warn  /c /O3 /nologo %INC%
  ifort /exe:Formal_Charges *.obj %CRYSLIB%
  goto END
:DEB
   ifort Formal_Charges.f90 /warn  /c /debug=full /traceback /nologo %INC%
   ifort /exe:Formal_Charges *.obj %CRYSLIB%
:END
rem
rem Compress executable
rem
   upx Formal_Charges.exe
rem
   if exist %FULLPROF% copy Formal_Charges.exe %FULLPROF%
   move Formal_Charges.exe %DISTRIB%
rem
rem Clean several files
rem
   del *.obj *.o *.mod *.map *.bak *.exe
