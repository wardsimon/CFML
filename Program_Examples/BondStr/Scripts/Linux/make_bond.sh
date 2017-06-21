#!/bin/sh
# Script to compile the Program: Bond_Str
#
if [ a$1 = a ]
then
cat << !
make_bond : Make the Bond_Str Program for Absoft/G95/ifort/Lahey Compilers (Linux)
Syntax    : make_bond af95:g95:ifort:lf95
!
exit
fi
#
# Define Main Directories (defined as environment variables in .bashrc)
# For instance:
# CRYSFML=$HOME/CrysFML_svn/crysfml_int
#
# BondStr Directory
#
cd $CRYSFML/Program_Examples/BondStr/Src
#
# Compiler Name
#
COMP=$1
if [ $COMP != "lf95" ]
then
   if [ $COMP != "g95" ]
   then
      if [ $COMP != "ifort" ]
      then
         if [ $COMP != "af95" ]
         then
            echo "Compiler Name was wrong!!!"
            cd ..
            exit
         fi
      fi
   fi
fi
DEBUG=nodeb
if [ $# -eq 2 ]
then
   DEBUG=$2
fi
#
# Compilation Options
#
if [ $DEBUG = "deb" ]
then
   case $COMP
   in
      'af95')
          OPT1="-c -g -O0"
            ;;
      'g95')
          OPT1="-c -g"
            ;;
      'ifort')
          OPT1="-c -g"
            ;;
      'lf95')
          OPT1="-c -g"
             ;;
   esac
else
   case $COMP
   in
      'af95')
          OPT1="-c -O2 -N11"
            ;;
      'g95')
          OPT1="-c -O"
            ;;
      'ifort')
          OPT1="-c -O -w -vec-report0"
            ;;
      'lf95')
          OPT1="-c -O --nchk --tpp"
             ;;
   esac
fi
#
# External Libraries Options
#
case $COMP
in
   'af95')
      INC="-I$CRYSFML/Absoft/LibC"
      LIB="-L$CRYSFML/Absoft/LibC"
      LIBSTATIC="-lcrysfml -lpthread"
      LIBDYNAMIC="-lcrysfml"
     ;;
   'g95')
      INC="-I$CRYSFML/G95/LibC"
      LIB="-L$CRYSFML/G95/LibC"
      LIBSTATIC="-lcrysfml -lpthread"
      LIBDYNAMIC="-lcrysfml"
     ;;
   'ifort')
      INC="-I$CRYSFML/ifort/LibC"
      LIB="-L$CRYSFML/ifort/LibC"
      LIBSTATIC="-lcrysfml"
      LIBDYNAMIC="-lcrysfml"
     ;;
   'lf95')
      INC="--mod .:$CRYSFML/Lahey/LibC"
      LIB="-L$CRYSFML/Lahey/LibC"
      LIBSTATIC="-lcrysfml -lpthread"
      LIBDYNAMIC="-lcrysfml"
     ;;
esac
#
# Compilation Process
#
echo " ########################################################"
echo " #### Bond_Str Program                         (1.0) ####"
echo " #### JRC - JGP                        CopyLeft-2009 ####"
echo " ########################################################"
$COMP $OPT1  CFML_Percolation.f90 $INC
$COMP $OPT1  Bond_Str.f90  $INC
case $COMP
in
  'af95')
     $COMP -o bond_str *.o -static $LIB $LIBSTATIC
     ;;
  'g95')
     $COMP -o bond_str *.o -static $LIB $LIBSTATIC
     ;;
  'ifort')
     $COMP -o bond_str *.o -static-intel $LIB $LIBSTATIC
     #$COMP -o bond_str *.o $LIB $LIBDYNAMIC
     ;;
  'lf95')
     $COMP *.o --out bond_str --staticlink $LIB $LIBDYNAMIC
     #$COMP *.o  --out bond_str -static $LIB $LIBSTATIC
     ;;
esac
#
# Compress BondStr Program
#
upx --compress-icons=0 bond_str
#
# Purge files
#
rm -rf *.o *.mod
cd $CRYSFML/Program_Examples/BondStr/Scripts/Linux
#
