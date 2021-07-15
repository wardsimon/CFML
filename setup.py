import os
import pathlib
import pkgutil
import shutil
import struct
import sys
import setuptools
import distutils.sysconfig as sysconfig
from distutils.core import setup
from distutils.command.install_data import install_data
from subprocess import CalledProcessError, check_output, check_call

from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install_lib import install_lib
from setuptools.command.install_scripts import install_scripts

BITS = struct.calcsize("P") * 8
PACKAGE_NAME = "crysfml_api"
SOURCE_DIR = '.'
COMPILER = 'gfortran'
if os.environ.get('FC', False):
    COMPILER = os.environ.get('FC')
print(f'Compiler set to: {COMPILER}')

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


class CMakeExtension(Extension):
    """
    An extension to run the cmake build

    This simply overrides the base extension class so that setuptools
    doesn't try to build your sources for you
    """

    def __init__(self, name, sources=[]):
        super().__init__(name=name, sources=sources)


class InstallCMakeLibsData(install_data):
    """
    Just a wrapper to get the install data into the egg-info

    Listing the installed files in the egg-info guarantees that
    all of the package files will be uninstalled when the user
    uninstalls your package through pip
    """

    def run(self):
        """
        Out files are the libraries that were built using cmake
        """
        # There seems to be no other way to do this; I tried listing the
        # libraries during the execution of the InstallCMakeLibs.run() but
        # setuptools never tracked them, seems like setuptools wants to
        # track the libraries through package data more than anything...

        self.outfiles = self.distribution.data_files


class InstallCMakeLibs(install_lib):
    """
    Get the libraries from the parent distribution, use those as the outfiles

    Skip building anything; everything is already built, forward libraries to
    the installation step
    """

    def run(self):
        """
        Copy libraries from the bin directory and place them as appropriate
        """

        self.announce("Moving library files", level=3)
        # We have already built the libraries in the previous build_ext step
        self.skip_build = True

        # Depending on the files that are generated from your cmake
        # build chain, you may need to change the below code, such that
        # your files are moved to the appropriate location when the installation
        # is run

        # bin_dir = self.distribution.bin_dir
        bin_dir = os.path.abspath(os.path.join(self.distribution.bin_dir, '..'))

        libs = [os.path.join(bin_dir, _dir) for _dir in
                os.listdir(bin_dir) if
                os.path.isdir(os.path.join(bin_dir, _dir))]

        for lib in libs:
            shutil.move(lib, os.path.join(self.build_dir,
                                          os.path.basename(lib)))

        # Move the lib to the correct location.
        bin_dir = self.build_dir
        pyd_path = [os.path.join(bin_dir, _pyd) for _pyd in
                    os.listdir(bin_dir) if
                    os.path.isfile(os.path.join(bin_dir, _pyd)) and
                    os.path.splitext(_pyd)[0].startswith(PACKAGE_NAME) and
                    os.path.splitext(_pyd)[1] in [".pyd", ".so"]][0]
        shutil.move(pyd_path, os.path.join(os.path.split(pyd_path)[0], 'CFML_api', os.path.split(pyd_path)[1]))

        # Mark the libs for installation, adding them to
        # distribution.data_files seems to ensure that setuptools' record
        # writer appends them to installed-files.txt in the package's egg-info
        #
        # Also tried adding the libraries to the distribution.libraries list,
        # but that never seemed to add them to the installed-files.txt in the
        # egg-info, and the online recommendation seems to be adding libraries
        # into eager_resources in the call to setup(), which I think puts them
        # in data_files anyways.
        #
        # What is the best way?

        # These are the additional installation files that should be
        # included in the package, but are resultant of the cmake build
        # step; depending on the files that are generated from your cmake
        # build chain, you may need to modify the below code

        self.distribution.data_files = [os.path.join(self.install_dir,
                                                     os.path.basename(lib))
                                        for lib in libs]

        # Must be forced to run after adding the libs to data_files

        self.distribution.run_command("install_data")

        super().run()


class InstallCMakeScripts(install_scripts):
    """
    Install the scripts in the build dir
    """

    def run(self):
        """
        Copy the required directory to the build directory and super().run()
        """

        self.announce("Moving scripts files", level=3)

        # Scripts were already built in a previous step

        self.skip_build = True

        bin_dir = self.distribution.bin_dir
        scripts_dirs = []

        # scripts_dirs = [os.path.join(bin_dir, _dir) for _dir in
        #                 os.listdir(bin_dir) if
        #                 os.path.isdir(os.path.join(bin_dir, _dir))]
        #
        # for scripts_dir in scripts_dirs:
        #
        #     shutil.move(scripts_dir,
        #                 os.path.join(self.build_dir,
        #                              os.path.basename(scripts_dir)))

        # Mark the scripts for installation, adding them to
        # distribution.scripts seems to ensure that the setuptools' record
        # writer appends them to installed-files.txt in the package's egg-info

        self.distribution.scripts = scripts_dirs

        super().run()


class BuildCMakeExt(build_ext):
    """
    Builds using cmake instead of the python setuptools implicit build
    """

    def run(self):
        """
        Perform build_cmake before doing the 'normal' stuff
        """

        for extension in self.extensions:
            self.build_cmake(extension)
        super().run()

    def build_cmake(self, extension: Extension):
        """
        The steps required to build the extension
        """
        cfg = "Debug" if self.debug else "Release"
        self.announce("Preparing the build environment", level=3)

        build_dir = pathlib.Path(self.build_temp)

        extension_path = pathlib.Path(self.get_ext_fullpath(extension.name))

        os.makedirs(build_dir, exist_ok=True)
        os.makedirs(extension_path.parent.absolute(), exist_ok=True)

        # Now that the necessary directories are created, build

        self.announce("Configuring cmake project", level=3)

        # Change your cmake arguments below as necessary
        # Below is just an example set of arguments for building Blender as a Python module

        cmake_args = [
            '-H' + SOURCE_DIR,
            '-B' + self.build_temp,
            "-DPYTHON_EXECUTABLE:FILEPATH={}".format(sys.executable),
            "-DARCH32=OFF",
            "-DCMAKE_Fortran_COMPILER={}".format(COMPILER),
            "-DPYTHON_API=ON",
            "-DUSE_HDF=OFF",
            "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(os.path.join(build_dir), 'Release'),
            "-DPYSETUP=ON"
        ]
        build_args = []

        check_call(
            [get_cmake()] + cmake_args
        )

        self.announce("Building binaries", level=3)

        check_call(
            [get_cmake(), "--build", ".", '--target', 'install'] + build_args, cwd=self.build_temp
        )

        # Build finished, now copy the files into the copy directory
        # The copy directory is the parent directory of the extension (.pyd)
        self.announce("Moving built python module", level=3)

        bin_dir = os.path.join(os.getcwd(), COMPILER, 'Python_API', 'CFML_api')
        self.distribution.bin_dir = bin_dir

        pyd_path = [os.path.join(bin_dir, _pyd) for _pyd in
                    os.listdir(bin_dir) if
                    os.path.isfile(os.path.join(bin_dir, _pyd)) and
                    os.path.splitext(_pyd)[0].startswith(PACKAGE_NAME) and
                    os.path.splitext(_pyd)[1] in [".pyd", ".so"]][0]

        shutil.move(pyd_path, extension_path)
        # After build_ext is run, the following commands will run:
        #
        # install_lib
        # install_scripts
        #
        # These commands are subclassed above to avoid pitfalls that
        # setuptools tries to impose when installing these, as it usually
        # wants to build those libs and scripts as well or move them to a
        # different place. See comments above for additional information


setup(name="CFML",
      version="0.0.1",
      author="Simon Ward",
      author_email="simon.ward@ess.eu",
      description="The Crystallographic Fortran Modules Library (CrysFML) is a set of modules containing "
                  "procedures of interest in Crystallographic applications.",
      ext_modules=[CMakeExtension(name=PACKAGE_NAME)],
      long_description=open("./README", 'r').read(),
      long_description_content_type="text/markdown",
      keywords="crystallography, physics, neutron, diffraction",
      classifiers=["Intended Audience :: Developers",
                   "License :: OSI Approved :: "
                   "GNU Lesser General Public License v3 (LGPLv3)",
                   "Natural Language :: English",
                   "Programming Language :: Fortran",
                   "Programming Language :: Python"],
      license='LGPL',
      cmdclass={
          'build_ext':       BuildCMakeExt,
          'install_data':    InstallCMakeLibsData,
          'install_lib':     InstallCMakeLibs,
          'install_scripts': InstallCMakeScripts
      },
      setup_requires=['wheel']
      )
