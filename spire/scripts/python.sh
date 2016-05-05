#!/bin/sh

# an executable script that runs Spire's installed version of python
#
# Copyright (C) 2006-2008  Health Research Inc.
#
# HEALTH RESEARCH INCORPORATED (HRI),
# ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455
# Email:  spider@wadsworth.org
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# the following line gets replaced by install.sh
SPIRE_DIR=__SPIRE_INSTALLATION_DIRECTORY__

BIN_DIR=$SPIRE_DIR/bin
LIB_DIR=$SPIRE_DIR/lib

tmp=$(basename $(ls -d $LIB_DIR/tcl8*))
TCL_VERSION=$(echo $tmp | cut -c4-6)

tmp=$(basename $(ls -d $LIB_DIR/python*))
PYTHON_VERSION=$(echo $tmp | cut -c7-9)

LD_LIBRARY_PATH=$LIB_DIR:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

TCL_LIBRARY="$LIB_DIR/tcl$TCL_VERSION"
export TCL_LIBRARY
TCLLIBPATH="$LIB_PATH"
export TCLLIBPATH
unset TK_LIBRARY

PYTHON="$BIN_DIR/python$PYTHON_VERSION"
unset PYTHONHOME
unset PYTHONPATH

exec "$PYTHON" "$@"
