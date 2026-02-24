#!/usr/bin/env python
# plotview.py: Plots number of views vs associated particles.
#
# Usage: python plotview.py data_extension

from Spider.Spiderutils import *
from Spider.Spiderscripts import extFromCommand
import os, sys

# ----------- Input files ---------
howmany = 'how_many'   # List of defocus groups

# ----------- Output files --------------
outfile = 'plot_view.gp'    # output file of gnuplot commands

# ---------- END BATCH HEADER -------------

ext = extFromCommand(sys.argv)

howmany += ext

cmds = """#!/bin/sh
set ylabel 'Particles'
set xlabel 'View'
set title 'Views for All Groups'
plot '%s' using 1:3 notitle with boxes
""" % howmany

fileWriteLines(outfile, cmds)

os.system('gnuplot -persist %s' % outfile)
