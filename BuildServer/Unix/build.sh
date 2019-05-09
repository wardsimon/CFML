compiler=$1
echo "Building CrysFML with ${compiler##/*}"

rm -Rf build_${compiler##/*}
mkdir build_${compiler##/*}
cd build_${compiler##/*}

cmake -D CMAKE_Fortran_COMPILER=${compiler} ..
cmake --build .