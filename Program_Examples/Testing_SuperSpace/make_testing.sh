#!/bin/sh
# Script to compile the Program: space_group_info
#
if [ a$1 = a ]
then
cat << !
make_read_ssg_datafile : Make the read_ssg_datafile Program for Absoft/G95/GFortran/Intel/Lahey Compilers (Linux)
Syntax                : make_read_ssg_datafile af95:g95:ifort:lf95:gfortran
!
exit
fi
#
# Compiler Name
#
COMP=$1
if [ $COMP != "lf95" ]
then
   if [ $COMP != "g95" ]
   then
     if [ $COMP != "gfortran" ]
     then
        if [ $COMP != "ifort" ]
        then
           if [ $COMP != "af95" ]
           then
              echo "Compiler Name was wrong!!!"
              exit
           fi
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
      'gfortran')
          OPT1="-c -g -ffree-line-length-none"
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
      'gfortran')
	  OPT1="-c -O -ffree-line-length-none"
            ;;
      'ifort')
          OPT1="-c -O -w"
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
      INC="-I../../Absoft/LibC"
      LIB="-L../../Absoft/LibC"
      LIBSTATIC="-lcrysfml"
      LIBDYNAMIC="-lcrysfml"
     ;;
   'g95')
      INC="-I../../G95/LibC"
      LIB="-L../../G95/LibC"
      LIBSTATIC="-lcrysfml"
      LIBDYNAMIC="-lcrysfml"
     ;;
   'gfortran')

      if [ $DEBUG = "deb" ]
      then
        INC="-I$CRYSFML/GFortran64_debug/LibC"
        LIB="-L$CRYSFML/GFortran64_debug/LibC"
        LIBSTATIC="-lcrysfml"
        LIBDYNAMIC="-lcrysfml"
      else
        INC="-I$CRYSFML/GFortran64/LibC"
        LIB="-L$CRYSFML/GFortran64/LibC"
        LIBSTATIC="-lcrysfml"
        LIBDYNAMIC="-lcrysfml"
     fi
     ;;
   'ifort')
      if [ $DEBUG = "deb" ]
      then
        INC="-I$CRYSFML/ifort64_debug/LibC"
        LIB="-L$CRYSFML/ifort64_debug/LibC"
        LIBDYNAMIC="-lcrysfml"
        LIBSTATIC="-lcrysfml"
      else
        INC="-I$CRYSFML/ifort64/LibC"
        LIB="-L$CRYSFML/ifort64/LibC"
        LIBDYNAMIC="-lcrysfml"
        LIBSTATIC="-lcrysfml"
      fi
     ;;
   'lf95')
      INC="--mod .:../../Lahey/LibC"
      LIB="-L../../Lahey/LibC"
      LIBDYNAMIC="-lcrysfml"
      LIBSTATIC="-lcrysfml -lpthread"
     ;;
esac
#
# Compilation Process
#
$COMP $OPT1 CFML_Rational_Arithmetic_test.f90  $INC
$COMP $OPT1 Matrix_Mod.f90                     $INC
$COMP $OPT1 CFML_ssg_datafile.f90              $INC
$COMP $OPT1 CFML_SuperSpaceGroups.f90          $INC
$COMP $OPT1 testing_ssg.f90                    $INC

case $COMP
in
  'af95')
     $COMP *.o  -o testing_ssg -static $LIB $LIBSTATIC
     ;;
  'g95')
     $COMP *.o  -o testing_ssg  $LIB $LIBSTATIC
     ;;
  'gfortran')
     $COMP *.o  -o testing_ssg  $LIB $LIBSTATIC
     ;;
  'ifort')
     $COMP *.o -o testing_ssg -static-intel $LIB $LIBSTATIC
     ;;
  'lf95')
     $COMP *.o --out testing_ssg --staticlink $LIB $LIBDYNAMIC
     ;;
esac
rm -rf *.o *.mod
