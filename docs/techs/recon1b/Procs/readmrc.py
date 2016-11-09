#!/usr/bin/env python

# Updated 2015-10-18
# 
# Read output report from MRC CTFFIND, print out defocus and other values.
#
# Usage: readmrc.py mrc_report [key_number]  > outdoc.dat
#
# The micrograph number is obtained from the report filename (e.g., report003).
# The doc file key may be given as an optional 2nd argument. If none
# is given, the micrograph number is used as the key.
# 
# Output: 7 columns
# ; KEY       MIC_NUM AVG_DEFOCUS  MAJORAXIS   MINORAXIS      ANGAST    ASTIGMAT     RESOLUTION
#    1 7    2481.000   20432.565   20318.110   20547.020     -12.969    -228.910       4.012

import sys
import string
import re

re_filenum = re.compile('[0-9]+$')

def get_mic_number(s):
    a = string.rfind(s, '.')   # remove extension
    if a > -1:
        s = s[:a]
    m = re_filenum.search(s)
    b = m.span()
    i = string.atoi(s[b[0]:b[1]])
    return i

nargs = len(sys.argv[1:])
#print "nargs: %s, type: %s" % (nargs, type(sys.argv))

if nargs == 0:
    print "Usage: readmrc.py mrc_report [key]"
    sys.exit(0)

# '-v' flag will write header columns
if sys.argv[-1] == '-v' : 
    print " ; KEY       MIC_NUM AVG_DEFOCUS  MAJORAXIS   MINORAXIS      ANGAST    ASTIGMAT"

    # remove from list so that get_mic_number doesn't get confused
    del sys.argv[-1]
    nargs -= 1

if nargs > 0 :
    report = sys.argv[1]
    try: 
        mic = get_mic_number(report)
    except AttributeError:
        print "Doesn't end in number -- using '1'"
        mic = 1

if nargs > 1:
    key = string.atoi(sys.argv[2])
else:
    key = mic
#print "key: ", key

fp = open(report)
B = fp.readlines()
fp.close

for line in B:
    if string.find(line, "Final Values") > -1: break

s = string.split(line)
d1 = string.atof(s[1])
d2 = string.atof(s[2])
davg = (d1 + d2) / 2.0
angle = string.atof(s[3]) - 45
if angle < 0 : angle += 90
astig = d1 - d2 # string.atof(s[5])
res = string.atof(s[6])
nitems = 7

print "%5d%2d%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f" % (key, nitems, key, davg, d1, d2, angle, astig, res)

