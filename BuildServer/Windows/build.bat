set compiler=%1
echo "Building CrysFML with %compiler%"

rmdir build_%compiler% /s /q
mkdir build_%compiler%
cd build_%compiler%

if %compiler% neq gfortran (
cmake -G "NMake Makefiles" -D ARCH32=OFF -D USE_HDF=ON -D CMAKE_Fortran_COMPILER=%compiler% -D HDF5_INCLUDE_PATH=%HDF5_INCLUDE_PATH% ..
) else (
cmake -G "MinGW Makefiles" -D ARCH32=OFF -D USE_HDF=ON -D CMAKE_Fortran_COMPILER=%compiler% -D HDF5_INCLUDE_PATH=%HDF5_INCLUDE_PATH% ..
)

cmake --build .