#!/usr/bin/env python
#
# SOURCE:  mkfilenums.py 
#
# PURPOSE: make file numbers from a list of filenames. 
#          First file is the output file!
#
# Spider Python Library
# Copyright (C) 2006-2018  Health Research Inc., Menands, NY
# Email:    spider@health.ny.gov


import sys
import os

from Spider.Spiderutils   import writedoc, filenumber

def unique(L):
    " Remove duplicates from a list "
    D = {}
    for n in L:
        D[n] = 1
    return D.keys()

def backup(filename):
    if os.path.exists(filename):
        found_vacancy = 0
        tiebreaker = 0
        shortdir = os.path.basename(os.path.dirname(filename))
        while not found_vacancy:
            test_filename = filename + '_' + str(tiebreaker)
            if os.path.exists(test_filename):
                tiebreaker += 1
            else:
                found_vacancy = 1
                short_old = os.path.join(shortdir, os.path.basename(filename))
                short_new = os.path.join(shortdir, os.path.basename(test_filename))
                print 'Renamed', short_old, 'to', short_new
                os.rename(filename, test_filename)

nargs = len(sys.argv[1:])

if nargs < 2:
    print 'Usage: mkfilenums.py [-f] outdocfile file**pattern '
    print "\nMake file numbers from a list of filenames. First argument is the output file!"
    sys.exit()

if sys.argv[1] == '-f' :
    outfile    = sys.argv[2]
    filelist   = sys.argv[3:]

    backup(outfile)

else:
    outfile = sys.argv[1]
    filelist = sys.argv[2:]

    if os.path.exists(outfile):
        print
        print 'WARNING!', outfile, 'exists.'
        print 'Output filename should be the first argument.'
        print 'Delete pre-existing file, use a different filename, or use -f flag to back up file.'
        print
        sys.exit()

numdict  = {}
key      = 0
classnum = 1

for file in filelist:
    key          = key + 1
    numdict[key] = [filenumber(file),classnum]

#numdict = unique(numdict)    # Remove duplicates
numdict.keys().sort()

writedoc(outfile, numdict)
#writedoc(outfile, columns=[numdict])

length = len(numdict)
print length, "Keys written to", outfile
