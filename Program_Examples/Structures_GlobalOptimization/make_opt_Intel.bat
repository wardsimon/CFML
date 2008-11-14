del *.obj
ifort observ.f90              /c /O3 /nologo /IC:\CrysFML\Intel\LibC
ifort cost_functions.f90      /c /O3 /nologo /IC:\CrysFML\Intel\LibC
ifort Optim_General.f90       /c /O3 /nologo /IC:\CrysFML\Intel\LibC
ifort  /exe:Optim_General *.obj C:\CrysFML\Intel\LibC\CrysFML.lib /link /stack:64000000
