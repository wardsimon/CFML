
cd ..\src

del *.obj
del *.mod
del *.map

      lf95 -c cryscal_mod       -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c IO_console        -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c text_mod          -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c math_mod          -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c matrix_mod        -o1 -info                                            -nchk  -nstchk
      lf95 -c macros            -o1 -info                                            -nchk  -nstchk

      lf95 -c cryscal_ext       -o1 -info  -mod  .;c:\crysfml\Lahey\libC;c:\rw_lf95  -nchk  -nstchk
      lf95 -c obv_rev           -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c neutrons_mod      -o1 -info                                            -nchk  -nstchk
      lf95 -c xrays_mod         -o1 -info                                            -nchk  -nstchk
      lf95 -c atome_mod         -o1 -info                                            -nchk  -nstchk
      lf95 -c shannon_mod       -o1 -info                                            -nchk  -nstchk
      lf95 -c mag_table         -o1 -info                                            -nchk  -nstchk
      lf95 -c calculs           -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal_sfac      -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal_dist      -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c X_space           -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c mendel            -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c mu_calc           -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal_symm      -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c space_group       -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c transf            -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal_main      -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal           -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal_init      -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c inter_cons        -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c read_CFL          -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c read_INS          -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c read_PCR          -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c read_KEYW         -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c read_CIF_file     -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c SIR_to_INS        -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c read_final_y      -o1 -info                                            -nchk  -nstchk
      lf95 -c read_CELL         -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c niggli            -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c read_SHELX_HKL    -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c create_CFL        -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c search_hkl        -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c search_spgr       -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c sort              -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c stat_data         -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c pgf_file          -o1 -info                                            -nchk  -nstchk
      lf95 -c HELP              -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c create_HTML       -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c create_CIF        -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal_lsg_cfml  -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal_write     -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c read_nreport      -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk
      lf95 -c cryscal_news      -o1 -info  -mod  .;c:\crysfml\Lahey\libC             -nchk  -nstchk

 lf95 -out cryscal *.obj   -mod c:\crysfml\lahey\libC -lib c:\crysfml\lahey\libC\crysfml;c:\rw_lf95\realwin.lib
 upx cryscal.exe

del *.obj
del *.mod
del *.map

cd ..\scripts
