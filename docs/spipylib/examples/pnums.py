#!/usr/bin/env python

# pnums.py - gets the range of particle numbers associated with
#            each micrograph.
#
#  Output file lists total particles and 1st & last particle numbers.

from Spider.Spiderutils import *
from Spider.Spiderscripts import * 
import os, sys

# ----------- Input files --------------
FILENUMS = "../filenums"    # file numbers
coords = "coords/sndc****"  #  coordinate doc files with particle numbers

# ----------- Output files --------------
outfile = "order_picked"    # output doc file

# -------------- END BATCH HEADER -------------

ext = extFromCommand(sys.argv)
putExtension(ext)

headers = ['micrograph', '#particles', 'cum.total', 'first', 'last']

total = 0

filelist = useFilenums(readdoc(FILENUMS), coords)
lines = []
for file in filelist:
    P = readdoc(file, column=3)
    mic = filenumber(file)
    np = len(P)
    total += np
    first = P[0]
    last = P[-1]
    lines.append([mic, np, total, first, last])

writedoc(outfile, lines=lines, headers=headers)
