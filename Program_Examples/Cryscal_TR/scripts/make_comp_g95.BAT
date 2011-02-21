

cd ../src

del *.o
del *.mod
del *.map



      g95 -c -O3 -funroll-loops  -msse2  cryscal_mod.f90      -Ic:\crysfml\G95\libC

      g95 -c -O3 -funroll-loops  -msse2  IO_console.f90       -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  text_mod.f90         -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  math_mod.f90         -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  matrix_mod.f90
      g95 -c -O3 -funroll-loops  -msse2  macros.f90

      g95 -c -O3 -funroll-loops  -msse2  cryscal_ext.f90      -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  obv_rev.f90          -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  neutrons_mod.f90
      g95 -c -O3 -funroll-loops  -msse2  xrays_mod.f90
      g95 -c -O3 -funroll-loops  -msse2  atome_mod.f90
      g95 -c -O3 -funroll-loops  -msse2  shannon_mod.f90
      g95 -c -O3 -funroll-loops  -msse2  mag_table.f90
      g95 -c -O3 -funroll-loops  -msse2  calculs.f90          -Ic:\crysfml\G95\libC

      g95 -c -O3 -funroll-loops  -msse2  cryscal_sfac.f90     -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  cryscal_dist.f90     -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  X_space.f90          -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  mendel.f90           -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  mu_calc.f90          -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  cryscal_symm.f90     -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  space_group.f90      -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  transf.f90           -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  cryscal_main_g95.f90 -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  cryscal.f90          -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  cryscal_init.f90     -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  inter_cons.f90       -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  read_CFL.f90         -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  read_INS.f90         -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  read_PCR.f90         -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  read_KEYW.f90        -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  read_CIF_file.f90    -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  SIR_to_INS.f90       -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  read_final_y.f90
      g95 -c -O3 -funroll-loops  -msse2  read_CELL.f90        -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  niggli.f90           -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  read_SHELX_HKL.f90   -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  create_CFL.f90       -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  search_hkl.f90       -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  search_spgr.f90      -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  sort.f90             -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  stat_data.f90        -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  pgf_file.f90
      g95 -c -O3 -funroll-loops  -msse2  HELP.f90             -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  create_HTML.f90      -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  create_CIF.f90       -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  cryscal_lsg_cfml.f90 -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  cryscal_write.f90    -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  read_nreport.f90     -Ic:\crysfml\G95\libC
      g95 -c -O3 -funroll-loops  -msse2  cryscal_news.f90     -Ic:\crysfml\G95\libC


rem creation du .EXE

   g95  *.o -o cryscal_g95   -O3  -funroll-loops  -msse2  -LC:\CrysFML\G95\LibC -lcrysfml

   upx cryscal_g95.exe

cd ../scripts

