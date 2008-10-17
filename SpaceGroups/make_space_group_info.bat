del *.obj
del *.mod
lf95 -c space_group_info.f90    -info  -o1 -chk -mod ".;C:\crysFML\lahey\libC"
lf95 *.obj -out space_group_info  -o1 -mod ".;C:\crysFML\lahey\LibC" -lib C:\crysFML\lahey\libC\crysFML  -chk
