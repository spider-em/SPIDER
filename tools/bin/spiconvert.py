#!/usr/bin/env python
#
# A simple SPIDER image converter.
# Output format determined by extension.

from  PIL import Image
import           os, sys

if len(sys.argv[1:]) < 2:
    print "Syntax: spiconvert.py infile outfile"
    print "To save as a SPIDER file, must specify extension: .spi"
    sys.exit()

filename = sys.argv[1]
outfile  = sys.argv[2]

im = Image.open(filename)

if im.format == 'SPIDER':
    # Convert a Spider file, output determined by extension
    im = im.convert2byte()
    im.save(outfile)
else:
    # Convert a non-Spider file
    f, ext = os.path.splitext(outfile)
    if ext == '.spi':   
        im.save(outfile, 'SPIDER') # Convert to Spider
    else:
        im.save(outfile)  # Converts according to extension

