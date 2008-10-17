del *.obj
del *.mod
lf95 -c subgroups.f90    -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
lf95 *.obj -out subgroups  -o1 -mod ".;C:\crysFML\lahey\LibC" -lib C:\crysFML\lahey\libC\crysFML  -chk
