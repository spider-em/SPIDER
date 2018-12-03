#!/usr/bin/env python
#
# MODIFICATIONS:
#    2015-03-06 -- Adapted from emancoords2spiderdoc
#
# Spider Python Library
# Copyright (C) 2006-2018  Health Research Inc., Menands, NY
# Email:  spider@health.ny.gov

import os, string, sys

from Spider            import Spiderutils
from Spider            import SpiderImagePlugin

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

if sys.argv[1:]:
    file = sys.argv[1]
    outputdoc = sys.argv[2]
    #micrograph = Image.open(sys.argv[2])
    #outputdoc = sys.argv[3]

    dictF = {}  # initialize dictionary
    key = 0  # initialize key

    if os.path.exists(file) : 
        input = open(file,'r')
        listL = input.readlines()  # read line-by-line
        input.close()
	dummy = 1

#	# EMAN x-coordinates are for the left side of the box
#	# EMAN y-coordinates are for the bottom of the box FROM THE BOTTOM OF THE MICROGRAPH
#	mic_y = micrograph.size[1]

        # read contents
        for line in listL:
    	    key += 1
	    split = line.split()

	    xdim = float(split[2])
	    ydim = float(split[3])
	    xcoord = float(split[0])
	    ycoord = float(split[1])

	    xcenter = xcoord + xdim/2
	    ycenter = ycoord + ydim/2
#	    ycenter = mic_y - (ycoord + ydim/2)

            dictF[key] = [key,xcenter,ycenter,xcenter,ycenter,dummy]
#           print filenum
	
#        headers = ['XCOORD','YCOORD','PARTICLE','PEAK_HT','XDIM','YDIM']
        headers = ['XCOORD','YCOORD','PARTICLE','PEAK_HT']
	backup(outputdoc)
        if Spiderutils.writeSpiderDocFile(outputdoc,dictF, headers=headers, append=0):
            print 'Wrote', key, 'keys to %s' % os.path.basename(outputdoc)
        else:
            print "Error!", "Unable to write to %s" % os.path.basename(outputdoc)
    else:
        print "Error!", "Unable to read %s" % file
else:
    print "syntax: emanrctcoords2spiderdoc.py input_eman_coords output_spider_doc"
#    print "syntax: emancoords2spiderdoc.py input_eman_coords input_micrograph output_spider_doc"

