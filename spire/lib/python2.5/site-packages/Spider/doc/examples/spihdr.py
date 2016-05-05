#!/usr/bin/env python

# spihdr.py : check the header of a Spider file
#
# The Spider image format is used by SPIDER software, in processing
# image data from electron microscopy and tomography.
# 
# Spider home page:
# www.wadsworth.org/spider_doc/spider/docs/master.html
#
# details about the Spider image format:
# www.wadsworth.org/spider_doc/spider/docs/image_doc.html
#
# History:
# 2004-02-05    Created BB
#
# NB needs Python2.3 or higher

import os, string, struct, sys

SpiderHeader = {'1': 'nslice',
                '2': 'nrow',
                '3': 'irec',
                '4': 'nhistrec',
                '5': 'iform',
                '6': 'imami',
                '7': 'fmax',
                '8': 'fmin',
                '9': 'av',
                '10': 'sig',
                '11': 'ihist',
                '12': 'nsam',
                '13': 'labrec',
                '14': 'iangle',
                '15': 'phi',
                '16': 'theta',
                '17': 'gamma',
                '18': 'xoff',
                '19': 'yoff',
                '20': 'zoff',
                '21': 'scale',
                '22': 'labbyt',
                '23': 'lenbyt',
                '24': 'istack',
                '25': 'NOTUSED',
                '26': 'maxim',
                '27': 'imgnum',
                '28': 'lastindx',
                '29': 'unused',
                '30': 'unused',
                '31': 'Kangle',
                '32': 'phi1',
                '33': 'theta1',
                '34': 'psi1',
                '35': 'phi2',
                '36': 'theta2',
                '37': 'psi2'
                }

def isInt(f):
    try:
        i = int(f)
        if f-i == 0: return 1
        else:        return 0
    except:
        return 0

iforms = [1,3,-11,-12,-21,-22]

# There is no magic number to identify Spider files, so just check
# a series of header locations to see if they have reasonable values.
# Returns no.of bytes in the header, if it is a valid Spider header,
# otherwise returns 0
def isSpiderHeader(t):
    h = (99,) + t   # add 1 value so can use spider header index start=1
    # header values 1,2,5,12,13,22,23 should be integers
    if not isInt(h[1]): return 0
    if not isInt(h[2]): return 0
    if not isInt(h[5]): return 0
    if not isInt(h[12]): return 0
    if not isInt(h[13]): return 0
    if not isInt(h[22]): return 0
    if not isInt(h[23]): return 0
    # check iform
    iform = int(h[5])
    if not iform in iforms: return 0
    # check other header values
    labrec = int(h[13])   # no. records in file header
    labbyt = int(h[22])   # total no. of bytes in header
    lenbyt = int(h[23])   # record length in bytes
    #print "labrec = %d, labbyt = %d, lenbyt = %d" % (labrec,labbyt,lenbyt)
    if labbyt != (labrec * lenbyt): return 0
    # looks like a valid header
    return labbyt

def isSpiderImage(filename):
    fp = open(filename,'rb')
    nbytes = 39
    f = fp.read(nbytes*4)   # read 39 * 4 bytes
    fp.close()
    bigendian = 1
    packstr = '>%df' % nbytes
    t = struct.unpack(packstr,f)    # try big-endian first
    hdrlen = isSpiderHeader(t)
    if hdrlen == 0:
        bigendian = 0
        packstr = '<%df' % nbytes
        t = struct.unpack(packstr,f)  # little-endian
        hdrlen = isSpiderHeader(t)

    printSpiderHeader(t)
    return hdrlen

def printSpiderHeader(t):
    n = len(t)
    if n > 38: n = 38
    n = 28   
    h = (99,) + t   # add 1 value so can use spider header index start=1
    if int(h[24]) == 0:
        n = 24  ###### only print up to istack

    nrow = 0
    nsam = 0

    # get width of strings
    maxwidth = 0
    for i in range(1,n+1):
        key = str(i)
        text = SpiderHeader[key]
        w = len(text)
        if w > maxwidth:
            maxwidth = w

    w = maxwidth

    for i in range(1,n+1):
        key = str(i)
        text = SpiderHeader[key].ljust(w)
        val = h[i]
        print "%3d %s %f " % (i, text, val)
        if i == 12:
            nsam = int(val)
        elif i == 2:
            nrow = int(val)
            
    print "Image size: %d, %d" % (nsam, nrow)
    istack = h[24]
    if istack > 0:
        ns = int(h[26])
        print "Stack file with %d images" % ns
    elif istack < 0:
        ns = -istack
        print "Indexed stack for %d images" % ns

#####################################################################

if __name__ == "__main__":

    if not sys.argv[1:]:
        print "Syntax: python spihdr.py imagefile"
        sys.exit(1)

    hdrlen = isSpiderImage(sys.argv[1]) / 4
    #print "The header is %d floating point numbers" % hdrlen

