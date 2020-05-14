compiler=$1
echo "Building CrysFML with ${compiler##/*}"

rm -Rf build_${compiler##/*}
mkdir build_${compiler##/*}
cd build_${compiler##/*}

cmake -D ARCH32=OFF -D USE_HDF=ON -D CMAKE_Fortran_COMPILER=${compiler} -D HDF5_INCLUDE_PATH=${HDF5_INCLUDE_PATH} ..
cmake --build .