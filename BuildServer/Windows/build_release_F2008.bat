set compiler=%1
echo "Building CrysFML08 release"

rmdir build_release_F2008 /s /q
mkdir build_release_F2008
cd build_release_F2008

cmake -G "NMake Makefiles" -D CRYSFML08=ON -D ARCH32=OFF -D PROG_EX=OFF -D CMAKE_BUILD_TYPE=Release -D CMAKE_Fortran_COMPILER=ifort ..
cmake --build . --target install