echo off
rem
   cd ..\..\Src
rem
rem  Compilation of fortran files
rem
   ifort cryscalc_mod.f90        /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort IO_console.f90          /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort text_mod.f90            /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort math_mod.f90            /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort matrix_mod.f90          /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort macros.f90              /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc_ext.f90        /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort obv_rev.f90             /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort neutrons_mod.f90        /c /debug:full /check /traceback    /nologo   /F100000000
   ifort xrays_mod.f90           /c /debug:full /check /traceback    /nologo   /F100000000
   ifort atome_mod.f90           /c /debug:full /check /traceback    /nologo   /F100000000
   ifort shannon_mod.f90         /c /debug:full /check /traceback    /nologo   /F100000000
   ifort mag_table.f90           /c /debug:full /check /traceback    /nologo   /F100000000
   ifort calculs.f90             /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort therm.f90               /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc_sfac.f90       /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc_dist.f90       /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort X_space.f90             /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort mendel.f90              /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort mu_calc.f90             /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc_symm.f90       /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort space_group.f90         /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort transf.f90              /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc_main_ifort.f90 /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc.f90            /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
REM Get CFML version and date of compilation
      ifort cryscalc_cfml_ver.f90 /c  /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
      ifort /exe:cryscalc_cfml_ver cryscalc_cfml_ver.obj  "%CRYSFML%"\ifort64_debug\LibC\CrysFML.lib
      del  cryscalc_cfml_ver.obj
      cryscalc_cfml_ver
      ifort cc_cfml_ver.f90         /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
REM
   ifort cryscalc_init.f90       /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort inter_cons.f90          /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_CFL.f90            /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_INS.f90            /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_PCR.f90            /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_KEYW.f90           /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_CIF_file.f90       /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort SIR_to_INS.f90          /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_final_y.f90        /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_CELL.f90           /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort niggli.f90              /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_SHELX_HKL.f90      /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort create_CFL.f90          /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort search_hkl.f90          /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort search_spgr.f90         /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort sort.f90                /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort stat_data.f90           /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort pgf_file.f90            /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort HELP.f90                /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort create_HTML.f90         /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort create_report.f90       /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort create_CIF.f90          /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc_lsg_cfml.f90   /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc_write.f90      /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort read_nreport.f90        /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
   ifort cryscalc_news.f90       /c /debug:full /check /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort64_debug\libc
rem
rem  Linking everything
rem
   ifort /exe:cryscalc *.obj  "%CRYSFML%"\ifort64_debug\LibC\CrysFML.lib /link /stack:64000000
   ifort /exe:cryscalc_ifort *.obj /F100000000 "%CRYSFML%"\ifort64_debug\libc\crysfml.lib
rem
rem UPX Compression
rem
   upx --compress-icons=0 cryscalc.exe
   copy cryscalc.exe ..\..\DistFPS_64b
   if exist %FULLPROF% copy cryscalc.exe %FULLPROF%
rem
   del *.obj *.mod *.bak *.map *.exe
rem
rem End of Compilation
rem
   cd ..\Scripts\Windows
