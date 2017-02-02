
cd ..\src

del *.obj
del *.mod
del *.map

rem copy c:\rw_lf95\*.lib
rem copy c:\rw_lf95\*.mod




      lf95 -c cryscalc_mod         -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c IO_rw                -o1 -info  -mod  .;c:\rw_lf95                           -nchk  -nstchk
      lf95 -c text_mod             -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c math_mod             -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c matrix_mod           -o1 -info                                               -nchk  -nstchk
      lf95 -c macros               -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c cryscalc_mod2        -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk

      lf95 -c cryscalc_ext         -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c obv_rev              -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c neutrons_mod         -o1 -info                                               -nchk  -nstchk
      lf95 -c xrays_mod            -o1 -info                                               -nchk  -nstchk
      lf95 -c atome_mod            -o1 -info                                               -nchk  -nstchk
rem   lf95 -c shannon_mod          -o1 -info                                               -nchk  -nstchk  <<  included in atome_mod.F90
rem   lf95 -c mag_table            -o1 -info                                               -nchk  -nstchk  <<  included in atome_mod.F90
      lf95 -c calculs              -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c therm                -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c cryscalc_sfac        -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c cryscalc_dist        -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c X_space              -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c mendel               -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c mu_calc              -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c cryscalc_symm        -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c space_group          -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c transf               -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c cryscalc_main_lf95_w -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c cryscalc             -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk

REM Get CFML version and date of compilation
      lf95 cryscalc_cfml_ver
      del  cryscalc_cfml_ver.obj cryscalc_cfml_ver.map
      cryscalc_cfml_ver
      lf95 -c cc_cfml_ver          -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
REM

      lf95 -c cryscalc_init        -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c inter_w              -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c read_CFL             -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c read_INS             -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c read_PCR             -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c read_KEYW            -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c read_CIF_file        -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c create_archive_CIF   -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk

      lf95 -c SIR_to_INS           -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c read_final_y         -o1 -info  -mod  .;c:\rw_lf95                           -nchk  -nstchk
      lf95 -c read_CELL            -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c niggli               -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c read_SHELX_HKL       -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c create_CFL           -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c search_hkl           -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c search_spgr          -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c sort                 -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c stat_data            -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c pgf_file             -o1 -info                                               -nchk  -nstchk
      lf95 -c HELP                 -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c create_HTML          -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c create_report        -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c create_CIF           -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c cryscalc_lsg_cfml    -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c cryscalc_write       -o1 -info  -mod  .;%CRYSFML%\Lahey\libC;c:\rw_lf95      -nchk  -nstchk
      lf95 -c read_nreport         -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk
      lf95 -c cryscalc_news        -o1 -info  -mod  .;%CRYSFML%\Lahey\libC                 -nchk  -nstchk

 lf95 -out wcryscalc *.obj   -mod %CRYSFML%\lahey\libC -lib %CRYSFML%\lahey\libC\crysfml;c:\rw_lf95\realwin.lib
 upx wcryscalc.exe

 del *.obj
 del *.mod
 del *.map

 cd ..\scripts
