#!/bin/sh
#
# SOURCE:   spider/tools/install.sh
#
# PURPOSE:  Create an executable script, python, to run the 
#           newly installed version of python and scripts to
#           run SPIDER's python tools using this version of
#           python
#
# NOTE: The old PIL and Blt libraries which are used in some of SPIDER's tools
#       are no longer compatible with recent Linux distributions.  I have updated
#       to Python2.7 and Tcl8.6 to overcome PIL incompatibility but have been
#       unable to get the Blt forks to work with these updates.  I have thus 
#       been forced (a kludge) to distribute a old Python2.5 with old Tcl and
#       Blt libraries to support the five tools which need Blt widgets but do not
#       use PIL library.
#
#
# NOTE:     No spaces in bash env. variables!! (al)


if [ -z ${SPIDER_DIR} ]; 
   then echo "set ENVIRONMENT VARIABLE: SPIDER_DIR to location of your SPIDER distribution"; exit ;
else echo " SPIDER_DIR is:  $SPIDER_DIR";  fi

TOOLS_DIR=$SPIDER_DIR/tools
TOOLS_LIB_DIR=$TOOLS_DIR/lib
TOOLS_BIN_DIR=$TOOLS_DIR/bin
PYTHON_BIN_DIR=$TOOLS_DIR/bin-python

# Compile Python library & site packages using SPIDER's Python2.5 --------

PYTHON_VERSION=python2.5 
PYTHON=$PYTHON_BIN_DIR/$PYTHON_VERSION
PYTHON_LIB=$TOOLS_LIB_DIR/$PYTHON_VERSION 

echo ' Compiling:      '$PYTHON_LIB library
$PYTHON $PYTHON_LIB/compileall.py $PYTHON_LIB

# Compile Python library & site packages using SPIDER's Python2.7 --------

PYTHON_VERSION=python2.7 
PYTHON=$PYTHON_BIN_DIR/$PYTHON_VERSION
PYTHON_LIB=$TOOLS_LIB_DIR/$PYTHON_VERSION
 
echo ' Compiling:      '$PYTHON_LIB library
$PYTHON $PYTHON_LIB/compileall.py $PYTHON_LIB

# Set tools executables directory in users PATH ---------
echo ' ' ; echo ' '
$PYTHON $TOOLS_DIR/install.py

# Add alternative name for montage to avoid ImageMagick command -------
ln -sf "$TOOLS_BIN_DIR"/montage "$TOOLS_BIN_DIR"/montage-spi

exit
