compiler=$1
echo "Building CrysFML08 with ${compiler##/*}"

rm -Rf build_${compiler##/*}_F2008
mkdir build_${compiler##/*}_F2008
cd build_${compiler##/*}_F2008

cmake -D CRYSFML08=ON -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Debug -D CMAKE_Fortran_COMPILER=${compiler} ..
cmake --build .