set compiler=%1
echo "Building CrysFML with %compiler%"

rmdir build_%compiler% /s /q
mkdir build_%compiler%
cd build_%compiler%

if %compiler% neq gfortran (
cmake -G "NMake Makefiles" -D USE_HDF=ON -D CMAKE_Fortran_COMPILER=%compiler% ..
) else (
rem C++ and C compiles must be set. Also, a "trick" has been set on the build server to build with MinGW
rem See See http://hdf-forum.184993.n3.nabble.com/HDF5-and-MinGW-td3393676.html for more details about this "trick"
cmake -G "MinGW Makefiles" -D USE_HDF=ON -D CMAKE_Fortran_COMPILER=%compiler% -D CMAKE_CXX_COMPILER=g++.exe -D CMAKE_C_COMPILER=gcc.exe  -D MINGW_HDF_TRICK=1 ..
)

cmake --build .