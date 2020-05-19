compiler=$1
echo "Building CrysFML08 with ${compiler##/*}"

rm -Rf build_${compiler##/*}_08
mkdir build_${compiler##/*}_08
cd build_${compiler##/*}_08

cmake -D CRYSFML08=ON -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Release -D CMAKE_Fortran_COMPILER=${compiler} ..
cmake --build .