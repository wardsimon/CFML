set compiler=%1
echo "Building CrysFML08 with %compiler%"

rmdir build_%compiler%_08 /s /q
mkdir build_%compiler%_08
cd build_%compiler%_08

if %compiler% neq gfortran (
cmake -G "NMake Makefiles" -D CRYSFML08=ON -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Release -D CMAKE_Fortran_COMPILER=%compiler% ..
) else (
cmake -G "MinGW Makefiles" -D CRYSFML08=ON -D ARCH32=OFF -D CMAKE_BUILD_TYPE=Release -D CMAKE_Fortran_COMPILER=%compiler% ..
)

cmake --build .