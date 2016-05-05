#!/usr/bin/env python

import os, sys
from commands import getoutput

if __name__ == '__main__':

    if len(sys.argv[1:]) == 0:
        print "Usage: callers.py pattern"
        sys.exit()

    pattern = sys.argv[1]

    cmd = "grep %s *.py" % pattern
    res = getoutput(cmd)
    lines = res.split("\n")

    commentchars = ['#', '"', "'"]
    called = []

    for line in lines:
        a = line.find(":")
        if a > -1:
            module = line[:a]
            prog = line[a+1:].strip()
            if prog[0] in commentchars:
                continue
            if 'import' in prog:
                continue
            called.append( (module,prog) )

    for c in called:
        print "%s: %s" % (c[0], c[1])

        
