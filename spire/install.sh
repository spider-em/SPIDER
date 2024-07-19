#!/bin/sh
#
# SOURCE:   spider/spire/install.sh
#
# PURPOSE:  This script must be run before first use of Spire or 
#           SPIDER's python tools.  It ensures that the spire tar archive
#           has been extracted, then it compiles all the '*.py' scripts.
#           It then runs '$SPIRE_BIN_DIR/setup.py' to set some Spire variables. 
#           It also intializes  XML configuration files, 
#           Spire LocalVars fire , and user's '.spire' files
#         
# NOTE:     The old PIL and Blt libraries which are used in Spire and 
#           some of SPIDER's tools are no longer compatible with recent 
#           Linux distributions. I have updated to Python2.7 and Tcl8.6 
#           to overcome the PIL incompatibility but have been unable to 
#           get the Blt forks to work with these updates.  
#           I have  been forced (a kludge) to distribute a old Python2.5 
#           with old Tcl and Blt libraries to support Spire and the five 
#           tools which need Blt widgets but do not use the PIL library.
#
# NOTE:     No spaces in bash env. variables!! (al)

# Ensure that user has set environment variable for "SPIDER_DIR"
if [ -z ${SPIDER_DIR} ]; 
   then echo "set ENVIRONMENT VARIABLE: SPIDER_DIR to location of your SPIDER distribution"; exit ;
else echo " SPIDER_DIR is:  $SPIDER_DIR";  fi

#SPIRE_DIR=$SPIDER_DIR/spire-dist

SPIRE_DIR=`pwd`
SPIRE_BIN_DIR=$SPIRE_DIR/bin
SPIRE_LIB_DIR=$SPIRE_DIR/lib
#echo " SPIRE_BIN_DIR is:  $SPIRE_BIN_DIR";  

# Ensure that user has extracted the spire.tar archive"
if [ -a ${SPIRE_DIR}/spire.tar ] 
   then  echo " Already unpacked: spire.tar " 
else
   echo " Unpacking: spire.tar.gz  please wait" 
   gunzip spire.tar.gz 
   tar xf spire.tar 
   echo " Unpacked: spire.tar " 
fi
 
# \rm *.pyc */*.pyc */*/*.pyc */*/*/*.pyc */*/*/*/*.pyc */*/*/*/*/*.pyc */*/*/*/*/*/*.pyc */*/*/*/*/*/*/*.pyc */*/*/*/*/*/*/*/*.pyc

# Compile Python library & site packages using SPIDER's Python 2.5 --------

PYTHON_VERSION=python2.5 
PYTHON=$SPIRE_BIN_DIR/$PYTHON_VERSION
PYTHON_LIB=$SPIRE_LIB_DIR/$PYTHON_VERSION 

echo ' Compiling:      '$PYTHON_LIB library
$PYTHON $PYTHON_LIB/compileall.py -f $PYTHON_LIB

# Compile Python library & site packages using SPIDER's Python 2.7 --------

PYTHON_VERSION=python2.7 
PYTHON=$SPIRE_BIN_DIR/$PYTHON_VERSION
PYTHON_LIB=$SPIRE_LIB_DIR/$PYTHON_VERSION
 
echo ' Compiling:      '$PYTHON_LIB library
$PYTHON $PYTHON_LIB/compileall.py -f $PYTHON_LIB

echo ' ' ; echo ' '

# Initialize some Spire variables 
# Set spire  and tools executables directory in users PATH --------------

$PYTHON $SPIRE_DIR/setup.py

echo ' ' 

exit
