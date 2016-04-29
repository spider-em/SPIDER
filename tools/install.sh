#!/bin/sh
#
# SOURCE:   /usr8/spider/tools/install.sh
#
# PURPOSE:  Create an executable script, python, to run the 
#           newly installed version of python and scripts to
#           run SPIDER's python tools using this version of
#           python
#
# NOTE:     No spaces in bash env. variables!! (al)

TOOLS_DIR=$PWD

PY_BIN_DIR=$TOOLS_DIR/bin-python

TOOLS_BIN_DIR=$TOOLS_DIR/bin

# Set up the python script --------------------------------------------
PYTHON=$PY_BIN_DIR/python
echo ' 'Python in use: $PYTHON


# Put the tools installation directory into the python script ---------
replace=__TOOLS_INSTALLATION_DIRECTORY__
rm -rf $PYTHON
sed "s:$replace:${TOOLS_DIR}:" python.sh > $PYTHON

chmod 775 $PYTHON

# Compile python library & site packages using the installed python ---
cd lib/python2.5

$PYTHON compileall.py

$PYTHON compileall.py site-packages

# Make scripts for the application programs in bin --------------------
cd ../..
"$PYTHON" "install.py"

# Add alternative name for montage to avoid Imagemagick command -------
ln -sf "$TOOLS_BIN_DIR"/montage "$TOOLS_BIN_DIR"/montage-spi

echo " " ; echo " Installation finished"
