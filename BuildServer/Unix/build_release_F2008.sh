compiler=$1
echo "Building CrysFML08 release"

rm -Rf build_release_F2008
mkdir build_release_F2008
cd build_release_F2008

cmake -D CRYSFML08=ON -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Release -D CMAKE_Fortran_COMPILER=ifort ..
cmake --build .