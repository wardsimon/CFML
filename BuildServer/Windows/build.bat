set compiler=%1
echo "Building CrysFML with %compiler%"

set OPTIONS=""
if %compiler% neq ifort (
    set OPTIONS=-G "NMake Makefiles"
)
rmdir build_%compiler% /s /q
mkdir build_%compiler%
cd build_%compiler%

cmake %OPTIONS% -D CMAKE_Fortran_COMPILER=%compiler% ..
cmake --build .