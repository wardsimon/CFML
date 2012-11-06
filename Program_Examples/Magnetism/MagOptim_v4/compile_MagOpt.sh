g95  -c CFML_Polar.f90 -std=f2003 -I../../LibC
g95  -c CFML_Refcodes.f90 -std=f2003 -I../../LibC
g95  -c CFML_Optimization_SAn.f90 -std=f2003 -I../../LibC
g95  -c Prep_Input.f90 -std=f2003 -I../../LibC
g95  -c Cost_Magfunctions.f90 -std=f2003 -I../../LibC
g95  -c MagOptim.f90 -std=f2003 -I../../LibC

g95   *.o -o magopt -L../../LibC -lcrysfml -Wl
