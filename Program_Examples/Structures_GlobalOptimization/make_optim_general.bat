del *.obj
lf95 -c observ.f90             -tp  -stchk -chk -o0 -mod ".;c:\CrysFML\Lahey\LibC" -lib c:\CrysFML\Lahey\LibC\CrysFML
lf95 -c cost_functions.f90     -tp  -stchk -chk -o0 -mod ".;c:\CrysFML\Lahey\LibC" -lib c:\CrysFML\Lahey\LibC\CrysFML
lf95 -c Optim_General.f90      -tp  -stchk -chk -o0 -mod ".;c:\CrysFML\Lahey\LibC" -lib c:\CrysFML\Lahey\LibC\CrysFML
lf95  *.obj -out Optim_General -tp  -stchk -chk -o0 -mod ".;c:\CrysFML\Lahey\LibC" -lib c:\CrysFML\Lahey\LibC\CrysFML
