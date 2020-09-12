set compiler=%1
echo "Building CrysFML release"

rmdir build_release /s /q
mkdir build_release
cd build_release

cmake -G "NMake Makefiles" -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Release -D CMAKE_Fortran_COMPILER=ifort -D HEAP_ARRAYS=ON -D PYTHON_API=ON -D PYTHON_LIBRARY_PATH=%PYTHON_LIBRARY_PATH% -D USE_HDF=ON -D HDF5_INCLUDE_PATH=%HDF5_INCLUDE_PATH% -D HDF5_LIBRARY_PATH=%HDF5_LIBRARY_PATH% ..
cmake --build . --target install
ctest
