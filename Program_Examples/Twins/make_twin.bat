ifort /c Search_TwinLaws.f90 /O2 /nologo /IC:\CrysFML\Intel\LibC
ifort /exe:Search_TwinLaws *.obj C:\CrysFML\Intel\LibC\crysfml.lib
del *.obj
