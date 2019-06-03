set compiler=%1
echo "Building CrysFML with %compiler%"

rmdir build_%compiler% /s /q
mkdir build_%compiler%
cd build_%compiler%

if %compiler% neq gfortran (

cmake -G "NMake Makefiles" -D USE_HDF=ON -D CMAKE_Fortran_COMPILER=%compiler% ..
) else (
cmake -G "NMake Makefiles" -D USE_HDF=ON -D CMAKE_Fortran_COMPILER=%compiler% -D CMAKE_CXX_COMPILER=g++.exe -D CMAKE_C_COMPILER=gcc.exe  -D MINGW_HDF_TRICK=1 ..
)

cmake --build .