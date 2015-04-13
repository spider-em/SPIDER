#!/usr/bin/env python

# convert2spi.py : convert micrograph images to spider format.
# Works on 2 tiff formats:
#    tiled tif (from ZI scanner): calls zi2spi
#    single data blocks: uses SPIDER's CP FROM RAW
# To do: convert some 'strip' format tiffs created by some scanners.
#
# Usage:
# 1) call from the command line:
#   % convert2spi.py infile outfile
#
# 2) or import this module in your script and call the function:
#
# import convert2spi
# convert2spi.convert2spider(infile, outfile)

import os, sys
from commands import getoutput
import Image

from Spider.Spiderutils import findSpider

spider = findSpider()

def isGzipped(file):
    fn, ext = os.path.splitext(file)
    if ext == '.gz': return 1
    else: return 0

def isTiff(file):
    fn, ext = os.path.splitext(file)
    if ext.lower() in ['.tif', '.tiff']: return 1
    else: return 0

def tiffSingle(infile, outfile, tile, verbose=0):
    "convert tif to Spider; for tif files with a single data block"
    a = tile[0]
    # e.g., tile = [('raw', (0, 0, 3200, 2880), 342, ('I;16', 0, 1))]
    if len(a) < 4:
        print "Unable to obtain tif header information."
        return 0
    format, ssize, offset, decoder = a[0], a[1], a[2], a[3]
    if len(ssize) == 4:
        width, height = ssize[2], ssize[3]
    else:
        print "Unable to obtain image dimensions from tif header."
        #print tile
        return 0
    x = decoder[0]
    d = x.split(';')
    itype = d[0]
    mode = d[1]
    nbits = 0
    endian = 'little'
    for n in ['8', '16', '32', '64']:
        if mode.find(n) > -1:
            nbits = int(n)
            if mode.find('B') > -1:
                endian = 'big'
            break

    #print "format: %s" % format
    #print "size: %d %d" % (width, height)
    #print "type: %s" % itype
    #print "nbits %d" % nbits
    
    btxt = "CP FROM RAW\n"
    btxt += infile + "\n"
    btxt += "%d      ; bits/pixel\n" % nbits
    btxt += "%d, %d      ; columns, rows\n" % (width, height)
    btxt += "%d       ; header bytes to skip\n" % offset
    if nbits == 16:
        if endian == 'little': btxt += '1      ; most significant byte\n'
        else: btxt +=  '2     ;  most significant byte\n'
        btxt += 'N       ; fold negatives?\n'
    btxt += outfile + "     ; output file\n"
    btxt +=  'EN D\n'
    if verbose:
        print btxt

    batchfile = 'tmp121688.spi'
    fp = open(batchfile, 'w')
    fp.write(btxt)
    fp.close()

    fn, ext = os.path.splitext(outfile)
    dataext = ext[1:]
    bat, bext = os.path.splitext(batchfile)
    batext = bext[1:]

    spicmd = "%s %s/%s @%s" % (spider, batext, dataext, bat)
    if verbose: print spicmd
    res = getoutput(spicmd)
    if verbose: print res

    os.remove(batchfile)
    
def tiffTiled(infile, outfile, verbose=0):
    "convert tif to Spider; for tiled tif files from ZI scanner"
    cmd = "zi2spi %s %s 1" % (infile, outfile)
    if verbose:
        print cmd
    res = getoutput(cmd)
    if res != "":
        print res

def convert2spider(infile, outfile, verbose=0):
    ZIPFLAG = 0   # for use later if file needs to be rezipped
    
    if isGzipped(infile):
        fn, ext = os.path.splitext(infile)
        cmd = "gunzip %s" % infile
        if verbose:
            print "gunzipping %s" % os.path.basename(infile)
        res = getoutput(cmd)
        if res != "":
            print "res"
            sys.exit()
        infile = fn  # filename without '.gz'
        ZIPFLAG = 1

    if isTiff(infile):
        im = Image.open(infile)
        if len(im.tile) == 1:
            if verbose:
                print "%s is a tif file, with a single data block" % os.path.basename(infile)
            res = tiffSingle(infile, outfile, im.tile, verbose)
        else:
            if verbose:
                print "%s is a tif file, with tiled data" % os.path.basename(infile)
            res = tiffTiled(infile, outfile, verbose)

if __name__ == '__main__':

    if spider == "":
        print "Unable to find path to Spider"
        sys.exit()

    nargs = len(sys.argv[1:])
    if nargs < 2:
        print "Usage: convert2spi.py infile outfile"
        sys.exit(0)
    else:
        infile = sys.argv[1]
        outfile = sys.argv[2]

    convert2spider(infile, outfile, verbose=1)
