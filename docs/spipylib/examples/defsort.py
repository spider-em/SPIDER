#!/usr/bin/env python

# defsort: assign defocus values to groups
#
# Usage: python defsort.py data_extension
#
# Input : defocus doc file with 2 columns:
#   (1) micrograph#, (2) defocus value
#
# Output : defocus group doc file with 3 columns
#   (1) micrograph#, (2) defocus value, (3) defocus group
#   NB: OUTPUT FILE IS SORTED BY DEFOCUS VALUE!
from Spider.Spiderutils import *
from Spider.Spiderscripts import extFromCommand
import sys

# ----------- Input parameters ---------
maxdiff = 1000  # maximum difference for defocus values to be in the same group

# ----------- Input files --------------
defocus = "defocus"  # input defocus file

# ----------- Output files --------------
defsort = "def_sort"  # output defocus group file

# ---------- END BATCH HEADER -------------

ext = extFromCommand(sys.argv)

# now append the data extension to all the filenames  (but see also putExtension)
defocus += ext
defsort += ext  # this is equivalent to "defsort = defsort + ext"

M, D = readdoc(defocus, column=[1,2])  # get lists of micrograph numbers, and defocus values
n = len(M)

N = []              # create a new list, of (defocus, micrograph) pairs
for i in range(n):
    N.append( (D[i],M[i]) )
N.sort()           # sorts by defocus, since that's the first element

# create the output, [micrograph, defocus, def_group]
prevmic = N[0][1]
prevdef = N[0][0]
group = 1    # current def_group
G = [ (prevmic, prevdef, group) ]

for i in range(1,n):
    thismic = N[i][1]
    thisdef = N[i][0]
    if (thisdef-prevdef) > maxdiff:
        group = group + 1  # increment group number if difference exceeds maxdiff
        prevdef = thisdef
    G.append( (thismic, thisdef, group) )

#write output
labels = ["micrograph", "defocus", "def.group"]
writedoc(defsort, lines=G, headers=labels)
