#!/bin/sh
# -------------------------------------------------------------------------------------
# Compilation Script for the CrysFML library
# Author: Javier Gonzalez-Platas
# Date: Jan-2017
# OS X version
# Update: Feb 2016 - A.Filhol (more tests and error messages, added a progression bar)
# Update: May 2016 - A.Filhol ($ARCH instead of $arg in output message)
# -------------------------------------------------------------------------------------
#
cat << !
 #################################################################################
 #### make_crysfml.sh: Make the CrysFML Library using Fortran Mac Compilers   ####
 #### Jan 2017 - OS X El Capitan / Yosemite                                   ####
 #### Syntax : make_crysfml ifort:f95:g95:gfortran [m32|m64] [winter] [debug] ####
 #################################################################################

!
# Checking CrySFML Environment Variable
if [ -z "$CRYSFML" ]; then
   echo "?????"
   echo "????? Please, set the environment variable CRYSFML in your .bash_profile"
   echo "?????"
   exit 1
fi
#
# Calling Program
#
if [ -z "$1" ]; then
cat << !
 ?????
 ????? Syntax : make_crysfml f95:g95:gfortran:ifort [m32|m64] [winter] [debug]
 ????? Default: m64 debug:no wint:no
 ?????
 
!
exit 1
fi
#
# Default values for Arguments
#
COMP=""
ARCH="m64"
WINT="N"
DEBUG="N"
#
# Arguments
#
for arg in "$@"
do
   case "$arg" in
      "f95")
         COMP=$arg
         ;;
      "g95")
         COMP=$arg
         ;;
      "gfortran")
         COMP=$arg
         ;;
      "ifort")
         COMP=$arg
         ;;
      "m32")
         ARCH=$arg
         ;;
      "m64")
         ARCH=$arg
         ;;
      "deb"*)
         DEBUG="Y"
         ;;
      "win"*)
         WINT="Y"
         ;;
   esac
done
#
# Check Compiler name
#
if [ -z $COMP ]; then
   echo "?????"
   echo "????? Compiler name is wrong. Please check it!"
   echo "?????"
   exit 1
fi
#
# Detect if Directories exists
#
case "$COMP" in
   "f95")
        IDIR="Absoft"
        ;;
   "g95")
        IDIR="G95"
        ;;
   "gfortran")
        IDIR="GFortran"
        if [ "$ARCH" == "m64" ]; then
           IDIR="GFortran64"
        fi
        ;;
   "ifort")
        IDIR="ifort"
        if [ "$ARCH" == "m64" ]; then
           IDIR="ifort64"
        fi
        ;;
esac
if [ -d "$CRYSFML/$IDIR" ]; then
   if [ -d "$CRYSFML/$IDIR/LibC" ]; then
      echo " "
   else
      mkdir "$CRYSFML/$IDIR/LibC"
      mkdir "$CRYSFML/$IDIR/LibW"
   fi
else
   mkdir "$CRYSFML/$IDIR"
   mkdir "$CRYSFML/$IDIR/LibC"
   mkdir "$CRYSFML/$IDIR/LibW"
fi
#
# Compilation Options
#
if [ $DEBUG == "Y" ]; then
   case "$COMP" in
      "f95")
          OPT1="-c -g -O0"
          OPT2="-c -g -O0"
          ;;
      "g95")
          OPT1="-c -g -Wall"
          OPT2="-c -g -Wall"
          ;;
      "gfortran")
          OPT1="-c -g -fbacktrace -ffree-line-length-none"
          OPT2="-c -g -fbacktrace -ffree-line-length-none"
          ;;
      "ifort")
          OPT1="-c -g -warn -$ARCH"
          OPT2="-c -g -warn -$ARCH"
          ;;
   esac
else
   case "$COMP" in
      "f95")
          OPT1="-c -O2 -N11"
          OPT2="-c -O0 -N11"
          ;;
      "g95")
          OPT1="-c -O"
          OPT2="-c -O0"
          ;;
      "gfortran")
          OPT1="-c -O -ffree-line-length-none -funroll-loops"
          OPT2="-c -O0 -ffree-line-length-none -funroll-loops"
          ;;
      "ifort")
          OPT1="-c -O -warn -$ARCH -Qopt-report=0"
          OPT2="-c -O0 -warn -$ARCH -Qopt-report=0"
          ;;
   esac
fi
#
# External Libraries Options
#
INC=""
if [ "$WINT" == "Y" ]; then
   case "$COMP" in
      "f95")
          INC="-p$WINTER/lib.abm"
          ;;
      "g95")
          INC="-I$WINTER/lib.g95"
          ;;
      "gfortran")
          INC=""
          ;;
      "ifort")
          if [ "$ARCH" == "m32" ]; then
            INC="-I$WINTER/lib.ifi"
          else
            INC="-I$WINTER/lib.ifi64"
          fi
          ;;
   esac
fi
#
# Changing to the Src Directory
#
cd $CRYSFML/Src
#
# Compilation Process
#
echo " ----"
echo " ---- Selected architecture   : $ARCH"
echo " ---- Compiled files to folder: $CRYSFML/$IDIR/$DIRLIB"
echo " ----"
echo " "
progressionBar()
{
	n=`expr $n + 1`
	bar=${bar:0:$n}#
  printf "Compiling [$n/$ntotal]: $bar\r"
}
ntotal=41      # Nber of files to be compiled
n=0            # Progression bar intitialisation
bar=#          # Progression bar item
#------
progressionBar
if [ "$COMP" == "ifort" ]; then
   $COMP $OPT1  CFML_GlobalDeps_MacOs_Intel.f90
else
   $COMP $OPT1  CFML_GlobalDeps_MacOs.f90
fi
progressionBar
$COMP $OPT1  CFML_Math_Gen.f90
progressionBar
$COMP $OPT1  CFML_LSQ_TypeDef.f90
progressionBar
$COMP $OPT1  CFML_Spher_Harm.f90
progressionBar
$COMP $OPT1  CFML_Random.f90
progressionBar
$COMP $OPT1  CFML_FFTs.f90
progressionBar
$COMP $OPT1  CFML_String_Util.f90
progressionBar
if [ "$WINT" == "N" ]; then
   $COMP $OPT1  CFML_IO_Mess.f90       $INC
else
   $COMP $OPT1  CFML_IO_MessWin.f90    $INC
fi
progressionBar
$COMP $OPT1  CFML_Profile_TOF.f90
progressionBar
$COMP $OPT1  CFML_Profile_Finger.f90
progressionBar
$COMP $OPT1  CFML_Profile_Functs.f90
progressionBar
$COMP $OPT1  CFML_Math_3D.f90
progressionBar
$COMP $OPT1  CFML_Optimization.f90
progressionBar
$COMP $OPT1  CFML_Optimization_LSQ.f90
progressionBar
$COMP $OPT2  CFML_Sym_Table.f90
progressionBar
$COMP $OPT2  CFML_Chem_Scatt.f90
progressionBar
$COMP $OPT2  CFML_BVSpar.f90
progressionBar
$COMP $OPT1  CFML_Diffpatt.f90
progressionBar
$COMP $OPT2  CFML_Bonds_Table.f90
progressionBar
$COMP $OPT1  CFML_Cryst_Types.f90
progressionBar
$COMP $OPT1  CFML_Symmetry.f90
progressionBar
$COMP $OPT1  CFML_ILL_Instrm_Data.f90
progressionBar
$COMP $OPT1  CFML_EoS_Mod.f90
progressionBar
$COMP $OPT1  CFML_Reflct_Util.f90
progressionBar
$COMP $OPT1  CFML_Atom_Mod.f90
progressionBar
$COMP $OPT1  CFML_Geom_Calc.f90
progressionBar
$COMP $OPT1  CFML_Molecules.f90
progressionBar
$COMP $OPT1  CFML_Form_CIF.f90
progressionBar
$COMP $OPT1  CFML_Extinction_Correction.f90
progressionBar
$COMP $OPT1  CFML_Sfac.f90
progressionBar
$COMP $OPT1  CFML_SXTAL_Geom.f90
progressionBar
$COMP $OPT1  CFML_Propagk.f90
progressionBar
$COMP $OPT1  CFML_Export_Vtk.f90
progressionBar
$COMP $OPT1  CFML_Maps.f90
progressionBar
$COMP $OPT1  CFML_Conf_Calc.f90
progressionBar
$COMP $OPT1  CFML_Magnetic_Groups.f90
progressionBar
$COMP $OPT1  CFML_MagSymm.f90
progressionBar
$COMP $OPT1  CFML_Optimization_SAn.f90 $INC
progressionBar
$COMP $OPT1  CFML_Refcodes.f90
progressionBar
$COMP $OPT1  CFML_Msfac.f90
progressionBar
$COMP $OPT1  CFML_Polar.f90
#
# Making CrysFML Library
#
DIRLIB="LibC"
LIBNAME="libcrysfml.a"
if [ "$WINT" == "Y" ]; then
   DIRLIB="LibW"
   LIBNAME="libwcrysfml.a"
fi
#echo "DIRLIB: $DIRLIB  LIBNAME: $LIBNAME WINT: $WINT"
ar cr $LIBNAME *.o
mv *.mod $CRYSFML/$IDIR/$DIRLIB
mv *.a $CRYSFML/$IDIR/$DIRLIB
rm *.o
#
# Changing to default directory
#
cd $CRYSFML/Scripts
#-----
printf "\n\n"
echo "*** End of script make_crysfml ***"
