# -*- coding: utf-8 -*-
import os
import sys
from subprocess import CalledProcessError, check_output, check_call
import pkgutil
import re
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion
from distutils.sysconfig import get_python_inc
import distutils.sysconfig as sysconfig

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}

# We can use cmake provided from pip which (normally) gets installed at /bin
# Except that in the manylinux builds it's placed at /opt/python/[version]/bin/
# (as a symlink at least) which is *not* on the path.
# If cmake is a known module, import it and use it tell us its binary directory
if pkgutil.find_loader('cmake') is not None:
    import cmake
    CMAKE_BIN = cmake.CMAKE_BIN_DIR + os.path.sep + 'cmake'
else:
    CMAKE_BIN = 'cmake'

def get_cmake():
    return CMAKE_BIN

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class PythonExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.join('Python_API', 'Src')

class CMakeBuild(build_ext):

    def run(self):
        try:
            out = check_output([get_cmake(), '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build" +
                               " the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        rex = r'version\s*([\d.]+)'
        cmake_version = LooseVersion(re.search(rex, out.decode()).group(1))
        if cmake_version < '3.13.0':
            raise RuntimeError("CMake >= 3.13.0 is required")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cfg = "Debug" if self.debug else "Release"

        cmake_args = [
            "-DPYTHON_EXECUTABLE:FILEPATH={}".format(sys.executable),
            "-DARCH32=OFF",
            "-DCMAKE_Fortran_COMPILER=gfortran",
            "-DPYTHON_API=ON",
            "-DUSE_HDF=OFF",
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DEXAMPLE_VERSION_INFO={}".format(self.distribution.get_version()),
            "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm

        ]
        build_args = [] #["--clean-first"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        check_call(
            [get_cmake(), ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        check_call(
            [get_cmake(), "--build", ".", '--target', 'install'] + build_args, cwd=self.build_temp
        )


# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="CFML",
    version="0.0.1",
    author="Simon Ward",
    author_email="simon.ward@ess.eu",
    description="ManyLinux test of CrysFML",
    long_description="",
    ext_modules=[CMakeExtension("cfml_cmake")],
    packages=find_packages(os.path.join(os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name))), 'Python_API')),
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
