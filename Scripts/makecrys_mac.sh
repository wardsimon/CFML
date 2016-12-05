#! /bin/sh
# -------------------------------------------------------------------------------------
# Compilation for CrysFML Library
# Author: Javier Gonzalez-Platas
# Date: July 2015
# A.Filhol - 26 february 2016 (more tests and more error messages)
# A.Filhol - 7 March 2016 (support for Winteracter 64-bits) 
# A.Filhol - 3 May 2016 (added chmod +x since this bit is often lost)
# A.Filhol - 13 May 2016 (all parameters have a default value)
# A.Filhol - 30 May 2016 (added the .sh extension to allow for file content search)
#
#    make_crysfml [?:help] [ifort:f95:g95:gfortran] [m64|m32] [winter|all] [debug]
#    Defaults: ifort m64 all nodebug 
# -------------------------------------------------------------------------------------

cat << !
 ########################################################
 makecrys_MAC  : Make the CrysFML Library for Apple's OS X
 feb 2016 - OS X El Capitan
 ########################################################

!
#
# Checking CrySFML Environment Variable
if [ -z "$CRYSFML" ]; then
   echo "?????"
   echo "????? Please, set the environment variable CRYSFML in your .bash_profile"
   echo "?????"
   exit 1
fi

# -- no longer neccasry since all parameters have a default value
#
# Calling Program
#
# if [ -z "$1" ]
# then
# cat << !
# ????? Missing parameters
# ????? Syntax : make_crysfml ifort:f95:g95:gfortran [m64|m32] [winter|all] [debug]
# ????? Default: m64  winter  nodebug     (winter->libW all-> libW and libC)
# ????? CAUTION: FullProf requirements:  ifort m64 all
# !
# exit 1
# fi

#
# Default values for Arguments
#
COMP="ifort"  # compiler
ARCH="m64"    # architecture
CONS="Y"      # console version
WINT="Y"      # winteracter
DEBUG="N"     # for debugging

#
# Arguments
#
for arg in "$@"
do
   case "$arg" in
      "?"|"h"|"help")
         echo "Syntax of the command:"
         echo "  >make_crysfml [?:help] [ifort:f95:g95:gfortran] [m64|m32] [winter|all] [debug]"
         echo "Default: ifort m64 all nodebug     (winter-> libW all-> libW and libC)"
         echo "Caution: FullProf requirements = ifort m64 all"
         echo " "
         exit 1
         ;;
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
         CONS="N"
         WINT="Y"
         ;;
      "all")
         CONS="Y"
         WINT="Y"
         ;;
       *)
         echo "?????"
         echo "????? Unknown parameter: $arg"
         echo "?????"
         echo " "
         exit 1
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
# Check if Compiler is installed
#
compinst=`which $COMP`
#echo $compinst
if [ -z $compinst ]; then
   echo "?????"
   echo "????? Compiler $COMP is not found! Install it or source it."
   echo "?????*"
   exit 1
fi

echo "---- Selected options -----------------"
echo "---- Comp: $COMP, Arch: $ARCH, Wint: $WINT, Cons: $CONS, Debug: $DEBUG"
echo "---------------------------------------"
echo ""

#
chmod +x $CRYSFML/Scripts/MacOS/make_crysfml.sh
#

#
# Console version
#
if [ "$CONS" == "Y" ]; then
   OPT1=""
   if [ "$DEBUG" == "Y" ]; then
      OPT1="debug"
   fi      
   $CRYSFML/Scripts/MacOS/make_crysfml.sh $COMP $ARCH $OPT1
fi  


#
# Winteracter version
#
if [ "$WINT" == "Y" ]; then
   OPT2=""
   if [ "$DEBUG" == "Y" ]; then
      OPT2="debug"
   fi     
   OPT1="winter" 

   $CRYSFML/Scripts/MacOS/make_crysfml.sh $COMP $ARCH $OPT1 $OPT2
fi

echo "*** End of script makecrys_mac ***"