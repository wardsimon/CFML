set compiler=%1
echo "Building CrysFML with %compiler%"

rmdir build_%compiler% /s /q
mkdir build_%compiler%
cd build_%compiler%

cmake -G "NMake Makefiles" -D USE_HDF=ON -D CMAKE_Fortran_COMPILER=%compiler% ..
cmake --build .