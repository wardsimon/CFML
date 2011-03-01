#!/bin/sh  
# Tutorial for CrysFML     
# ./make_xtl; ./xtl test   

#gfortran -c -fPIC xtl2cfml.f95 -I../../GFortran/LibC
rm *.o *.so xtl
gfortran -c -fPIC -m32 CFML4C.f95   -I../../GFortran/LibC
gfortran -shared -o /usr/local/lib/libF4C.so -m32 CFML4C.o  -L../../GFortran/LibC -lcrysfml    
gcc -m32 xtl.c -lF4C  -o xtl
#rm *.o 
./xtl 


            