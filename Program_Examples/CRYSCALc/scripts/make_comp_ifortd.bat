

cd ..\src

del *.obj
del *.mod
del *.map


      ifort      cryscalc_mod.f90        /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      IO_console.f90          /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      text_mod.f90            /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      math_mod.f90            /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      matrix_mod.f90          /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      macros.f90              /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc_mod2.f90       /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
rem   ifort      cryscalc_ext.f90        /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc  << renamed as cryscalc_mod2.F90
      ifort      obv_rev.f90             /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      neutrons_mod.f90        /c /debug:full /check  /traceback    /nologo   /F100000000
      ifort      xrays_mod.f90           /c /debug:full /check  /traceback    /nologo   /F100000000
      ifort      atome_mod.f90           /c /debug:full /check  /traceback    /nologo   /F100000000
rem   ifort      shannon_mod.f90         /c /debug:full /check  /traceback    /nologo   /F100000000         << now included in atom_mod.F90
rem   ifort      mag_table.f90           /c /debug:full /check  /traceback    /nologo   /F100000000         << now included in atom_mod.F90
      ifort      calculs.f90             /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      therm.f90               /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc_sfac.f90       /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc_dist.f90       /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      X_space.f90             /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      mendel.f90              /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      mu_calc.f90             /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc_symm.f90       /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      space_group.f90         /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      transf.f90              /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc_main_ifort.f90 /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc.f90            /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc

REM Get CFML version and date of compilation
      ifort cryscalc_cfml_ver.f90
      del  cryscalc_cfml_ver.obj
      cryscalc_cfml_ver
      ifort      cc_cfml_ver.f90         /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
REM


      ifort      cryscalc_init.f90       /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      inter_cons.f90          /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_CFL.f90            /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_INS.f90            /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_PCR.f90            /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_KEYW.f90           /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_CIF_file.f90       /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      create_archive_CIF.f90  /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      SIR_to_INS.f90          /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_final_y.f90        /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_CELL.f90           /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      niggli.f90              /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_SHELX_HKL.f90      /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      create_CFL.f90          /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      search_hkl.f90          /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      search_spgr.f90         /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      sort.f90                /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      stat_data.f90           /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      pgf_file.f90            /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      HELP.f90                /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      create_HTML.f90         /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      create_report.f90       /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc

      ifort      create_CIF.f90          /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc_lsg_cfml.f90   /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc_write.f90      /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      read_nreport.f90        /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc
      ifort      cryscalc_news.f90       /c /debug:full /check  /traceback    /nologo   /F100000000    /I"%CRYSFML%"\ifort\libc



   ifort /exe:cryscalc       *.obj             "%CRYSFML%"\ifort_debug\LibC\CrysFML.lib /link /stack:64000000
   ifort /exe:cryscalc_ifort *.obj /F100000000 "%CRYSFML%"\ifort_debug\libc\crysfml.lib
 upx cryscalc_ifort.exe

rem del *.obj
rem del *.mod
rem del *.map

rem cd ..\scripts
