#!/bin/sh
#
# SOURCE:   /usr8/spider/spire/spire_linux-1.5.5/install.sh
#
# PURPOSE:  Create an executable script, python, that runs our 
#           newly installed version of python
#
# REQUIRES: The stub files, python.sh, spire.sh to be in the 
#           scripts directory

SPIRE_DIR=$PWD

BIN_DIR=$SPIRE_DIR/bin
LIB_DIR=$SPIRE_DIR/lib
SCRIPT_DIR=$SPIRE_DIR/scripts

# Set up the python script -----------------------
PYTHON=$BIN_DIR/python

# Put the Spire installation directory into the python script
replace=__SPIRE_INSTALLATION_DIRECTORY__
rm -rf $PYTHON
sed "s:$replace:${SPIRE_DIR}:" $SCRIPT_DIR/python.sh > $PYTHON

chmod 775 $PYTHON

cd lib/python2.5
$PYTHON compileall.py
$PYTHON compileall.py site-packages
cd ../..

# Install Spire and the spire script --------------
"$PYTHON" "setup.py"

# setup.py renames bin/spire.py to bin/spire
# But we want 'spire' to be a shell script that calls spire.py
mv $BIN_DIR/spire $BIN_DIR/spire.py 

# Put the Spire installation directory into the spire script
replace=__SPIRE_INSTALLATION_DIRECTORY__
sed "s:$replace:${SPIRE_DIR}:" $SCRIPT_DIR/spire.sh > $BIN_DIR/spire

chmod 775 $BIN_DIR/spire

#echo " " ; echo " Making: $mkapps.py" "$BIN_DIR"

# Make scripts for the application programs in bin/
"$PYTHON" "mkapps.py" "$BIN_DIR"

# Add alternative name for montage to avoid Imagemagick command
ln -sf "$BIN_DIR"/montage "$BIN_DIR"/montage-spi

echo " " ; echo " Installation finished"
