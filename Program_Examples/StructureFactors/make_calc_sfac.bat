del *.obj
lf95 -c calc_sfac.f90  -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC"
lf95  *.obj -out calc_sfac -tp -nomap -stchk -nchk -o1 -mod ".;c:\CrysFML\lahey\LibC" -lib c:\CrysFML\lahey\LibC\CrysFML
