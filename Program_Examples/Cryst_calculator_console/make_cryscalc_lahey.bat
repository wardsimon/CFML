lf95 menu_1.f90 -c  -o1 -mod ".;c:\CrysFML\Lahey\LibC"
lf95 menu_2.f90 -c  -o1 -mod ".;c:\CrysFML\Lahey\LibC"
lf95 menu_3.f90 -c  -o1 -mod ".;c:\CrysFML\Lahey\LibC"
lf95 menu_4.f90 -c  -o1 -mod ".;c:\CrysFML\Lahey\LibC"
lf95 calsym.f90 -c  -o1 -mod ".;c:\CrysFML\Lahey\LibC"
lf95  *.obj -out cryscalc -lib c:\CrysFML\Lahey\LibC\CrysFML
del *.obj *.map *.bak
