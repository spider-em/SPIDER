#!/bin/sh

# an executable script that runs Spire with its own installed version of python
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

PYTHON="$BIN_DIR/python"

SPIRE=$BIN_DIR/spire.py

exec "$PYTHON" "$SPIRE" "$@"
