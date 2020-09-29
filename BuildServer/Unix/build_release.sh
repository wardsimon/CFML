compiler=$1
echo "Building CrysFML with ${compiler*##/}"

rm -Rf build_release
mkdir build_release
cd build_release

cmake -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Release -D CMAKE_Fortran_COMPILER=ifort -DPYTHON_API=ON -DPYTHON_INTERPRETER_PATH=${PYTHON_INTERPRETER_PATH} -DPYTHON_LIBRARY_PATH=${PYTHON_LIBRARY_PATH} -D USE_HDF=ON -D HDF5_INCLUDE_PATH=${HDF5_INCLUDE_PATH} -D HDF5_LIBRARY_PATH=${HDF5_LIBRARY_PATH} ..
cmake --build . 
cmake --build . --target doxygen
eval "${TEST_TRICK_COMMAND}"
ctest
status=$?
if [ $status -ne 0 ]; then
	echo "Failure/Error during ctest"
	exit $status
fi
cmake --build . --target install
