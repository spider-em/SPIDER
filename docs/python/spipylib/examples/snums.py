#!/usr/bin/env python

# snums.py : lists selected particle numbers from each micrograph 
#
#  Output file lists total particles and 1st & last particle numbers.
           
from Spider.Spiderutils import *
from Spider.Spiderscripts import *
import os, sys

# ----------- Input files --------------
FILENUMS = "../filenums"     # file numbers
picked   = "order_picked"    #  doc file of all picked particles
good     = "good/ngood****"  # doc files of selected particle numbers

# ----------- Output files --------------
selected = "order_selected"    # output doc file
percent  = "percent_selected"  # doc file of picked vs selected

# -------------- END BATCH HEADER -------------

ext = extFromCommand(sys.argv)
putExtension(ext) 

# count up particles, write out the order_selected file
headers = ['micrograph', '#particles', 'cum.total', 'first', 'last']

total = 0
filelist = useFilenums(FILENUMS, good)
lines = []
for file in filelist:
    P = readdoc(file)   # get the list of particle numbers
    mic = filenumber(file)
    np = len(P)
    total += np
    first = P[0]
    last = P[-1]
    lines.append([mic, np, total, first, last])

writedoc(selected, lines=lines, headers=headers)

# now compute percent selected/picked
hdrs = ['micrograph', 'picked', 'selected', 'per cent']

Mic, Pick = readdoc(picked, column=[1,2])
Sel = readdoc(selected, column=2)

lines = []
n = len(Mic)
for i in range(n):
    lines.append( [Mic[i], Pick[i], Sel[i], Sel[i]/float(Pick[i])] )
writedoc(percent, lines=lines, headers=hdrs)

# make the last line with totals
ptot = readdoc(picked, column=3)[-1]
stot = readdoc(selected, column=3)[-1]   # just want the last value
line = [ 0, ptot, stot, stot/float(ptot) ]
hdrs[0] = 'totals'
writedoc(percent, lines=[line], headers=hdrs, mode='a')  # append
