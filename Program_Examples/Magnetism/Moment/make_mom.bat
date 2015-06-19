   ifort /c moment.f90 /Ox /nologo /I. /IC:\CrysFML\ifort\LibC
   link /subsystem:console /stack:102400000 /out:moment.exe *.obj C:\CrysFML\ifort\LibC\CrysFML.lib
   del *.obj *.mod
