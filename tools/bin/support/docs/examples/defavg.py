#!/usr/bin/env python

# defavg.py : compute average values for the defocus groups
#
# Usage: python defavg.py data_extension
#
# Input: defsort : doc file with defocus groupings
#
# Output
#   defgroup_avg : doc file with 4 columns:
#   micrograph, defocus, defocus_group, defocus_group_average
#
from Spider.Spiderutils import *
from Spider.Spiderscripts import *
import os, sys

# ----------- Input files ---------
defsort = 'def_sort'   # doc file with defocus groupings

# ----------- Output files --------------
def_avg = 'def_avg'    # output doc file with defocus averages

order = 'order_defgrps' # list of defocus groups and their averages

# ---------- END BATCH HEADER -------------

ext = extFromCommand(sys.argv)

putExtension(ext)   # this has same effect as next three lines, adding extension to filenames
#defsort += ext
#def_avg += ext
#order   += ext

M, D, G = readdoc(defsort, column=[1,2,3])  # micrographs,defocus values,defocus groups

# make a dictionary of defocus groups and their defoci
n = len(M)
DG = {}
for i in range(n):
    g = int(G[i])
    if g in DG:
        DG[g].append(D[i])  # each defgroup points to a list of defoci
    else:
        DG[g] = [D[i]]
        
# find avg of each list of defoci
for g in DG:
    dlist = DG[g]
    sum = 0
    for d in dlist:
        sum = sum + d
    avg = sum / float(len(dlist))
    DG[g] =  avg
    
# make a list of averages and write out to order file
defgrps = DG.keys()
defgrps.sort()  # list of defocus groups
avgs = []
for defgrp in defgrps:
    avgs.append(DG[defgrp])    

writedoc(order, [defgrps, avgs], headers=['defocus_grp', 'average'])

# make lists of [micrograph, defocus, defocus group, group average]
outlists = []
for i in range(n):
    outlists.append( [ M[i], D[i], G[i], DG[int(G[i])] ] )

outlists.sort()  # sort outlists (it's sorted by micrographs, the first element)
hdrs=['micrograph', 'defocus', 'def_group', 'def_grp_avg']
# write outlists to doc file as lines 
writedoc(def_avg, lines=outlists, headers=hdrs)


