set compiler=%1
echo "Building CrysFML with %compiler%"

rmdir build_%compiler% /s /q
mkdir build_%compiler%
cd build_%compiler%

if %compiler% neq gfortran (
cmake -G "NMake Makefiles" -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Debug -D CMAKE_Fortran_COMPILER=%compiler% -D USE_HDF=ON -D HDF5_INCLUDE_PATH=%HDF5_INCLUDE_PATH% -D HDF5_LIBRARY_PATH=%HDF5_LIBRARY_PATH% ..
) else (
cmake -G "MinGW Makefiles" -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Debug -D CMAKE_Fortran_COMPILER=%compiler% -D USE_HDF=ON -D HDF5_INCLUDE_PATH=%HDF5_INCLUDE_PATH% -D HDF5_LIBRARY_PATH=%HDF5_LIBRARY_PATH% ..
)

cmake --build .