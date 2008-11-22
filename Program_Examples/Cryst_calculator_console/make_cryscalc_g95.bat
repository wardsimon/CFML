g95 -c -O3    -funroll-loops  -msse2   menu_1.f90     -IC:\CrysFML\G95\LibC
g95 -c -O3    -funroll-loops  -msse2   menu_2.f90     -IC:\CrysFML\G95\LibC
g95 -c -O3    -funroll-loops  -msse2   menu_3.f90     -IC:\CrysFML\G95\LibC
g95 -c -O3    -funroll-loops  -msse2   menu_4.f90     -IC:\CrysFML\G95\LibC
g95 -c -O3    -funroll-loops  -msse2   calsym.f90     -IC:\CrysFML\G95\LibC
g95  *.o -o cryscalc -O3  -funroll-loops  -msse2  -LC:\CrysFML\G95\LibC -lcrysfml
del *.o *.mod *.bak *.map
