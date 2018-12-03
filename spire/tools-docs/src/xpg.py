#!/usr/bin/env python

libdir = "/usr/local/lib"   # edit this line to search a particular directory

"""
Runs command

nm libfile | grep xpg5

on all lib*** files in libdir
"""

import os, string
from commands import getoutput

os.chdir(libdir)

L = os.listdir(libdir)

for lib in L:
    if lib[:3] == 'lib':  # and not os.path.islink(lib):
        #cmd = "nm %s | grep xpg5" % lib
        cmd = "nm %s | grep xpg5_vsnprint" % lib
        out = getoutput(cmd)
        if out != None and out != "":
            print "%s =================" % lib
            print out

