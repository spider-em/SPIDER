#!/usr/bin/env python

import os, string, sys
from   Spider import Spiderutils

if sys.argv[1:]:
    file = sys.argv[1]
    filename = sys.argv[2]

    F = {}  # initialize dictionary
    key = 0  # initialize key

    if os.path.exists(file) : 
        input = open(file,'r')
        L = input.readlines()  # read line-by-line
        input.close()
    
        # read contents
        for line in L:
            filenum = Spiderutils.getfilenumber(line)
    	    key += 1
            F[key] = [filenum]
#           print filenum
	
        headers = ['file_number']
        if Spiderutils.writeSpiderDocFile(filename,F, headers=headers, append=0):
            print 'Wrote', key, 'keys to %s' % os.path.basename(filename)
        else:
            print "Error!", "Unable to write to %s" % os.path.basename(filename)
    else:
        print "Error!", "Unable to read %s" % file
else:
    print "Syntax: makefilenums.py inputtextfile outputspiderdoc"
    
