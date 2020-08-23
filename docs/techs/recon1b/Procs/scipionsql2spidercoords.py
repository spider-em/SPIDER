#!/usr/bin/env python

# Spider Python Library
# Copyright (C) 2006  Health Research Inc.
#
# HEALTH RESEARCH INCORPORATED (HRI),
# ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455
#
# Email:  spider@wadsworth.org
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

print
print "scipionsql2spidercoords.py, Modified 2016 Nov 09"  
print

import sqlite3
import sys
import os
import time

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


def getfilenumber(filename):
    "returns file number as a string with leading zeroes "
    # copied from Spiderutils.py
    
    filename = os.path.basename(filename)
    fname,ext = os.path.splitext(filename)

    numstr = ""
    f = list(fname)
    f.reverse()
    done = 0
    for ch in f:
        if not done:
            try:
                int(ch)
                numstr = ch + numstr
            except:
                if numstr != "":
                    done = 1
    return numstr
            
            
def name2template(filename, all=0):
    " given 'mic021.dat' --> returns mic***.dat "
    " by default, only replaces number nearest extension. all !=0 replaces all"
    # copied from Spiderutils.py
    
    if len(filename) == 0: return ""
    path, basename = os.path.split(filename)
    fname,ext = os.path.splitext(basename)
    
    newfn = ""
    f = list(fname)
    f.reverse()
    if all:
        for ch in f:
            try:
                int(ch)
                newch = '*'
            except:
                newch = ch
            newfn = newch + newfn
    else:
        found = 0
        for ch in f:
            if not found:
                try:
                    int(ch)
                    newch = '*'
                except:
                    newch = ch
                    if newfn and newfn[0] == '*':
                        found = 1
            else:
                newch = ch
            newfn = newch + newfn
    fname = os.path.join(path,newfn) + ext
    return fname


def template2filename(template, n=0):  #numfile=None, n=None):
    "replaces asterisks with number: (pic***.dat, doc003.dat) returns pic003.dat"
    # template should have asterisks, num can be int or a numbered filename.
    # Like makeSpiderFilename, but it can accept a filename instead of a number.
    # copied from Spiderutils.py
    
    if type(n) == type(1):
        pass
    elif type(n) == type("string"):
        n = filenumber(n)
    else:
        print "template2filename: unable to parse input"
        return ""
    nstars = template.count("*")
    if nstars == 0:
        return template
    if len(str(n)) > nstars:
        print "template2filename: **** Warning number larger than template"
    numstr = str(n).zfill(nstars)
    sts = "*" * nstars
    filename = template.replace(sts,numstr)
    return filename


def makeDocfileHeader(filename, batext=None):
    "create the comment line used at the top of SPIDER document files"
    # copied from Spiderutils.py

    filename = os.path.basename(filename)
    fn, ext = os.path.splitext(filename)
    ext = ext[1:]
    if batext == None:
        batext = 'spl'   # Spider Python Library
    date,time,idstr = nowisthetime()
    h = " ;%s/%s   %s AT %s   %s\n" % (batext,ext,date,time,filename)
    return h


def nowisthetime():
    "return current time as tuple of 3 strings: (date, time, ID)"
    # returns e.g., ('16-OCT-03', '13:08:16', '031016130816')
    # copied from Spiderutils.py

    tt = time.localtime(time.time())
    # localtime return format: (2003, 10, 16, 12, 48, 30, 3, 289, 1)
    #t = string.split(time.asctime(tt))
    t = time.asctime(tt).split()
    # asctime return format: 'Thu Oct 16 12:50:17 2003'
    mo = t[1].upper()
    day = t[2]
    if len(day) < 2: day = '0' + day
    timestr = t[3]
    yr = t[4]
    datestr = "%s-%s-%s" % (day, mo, yr)

    yr = yr[-2:]
    # this is just to get the month as a number
    d = map(str,tt)   # stringify all numbers in the tuple
    mon = d[1]
    if len(mon) < 2: mon = '0' + mon
    #(h,m,s) = string.split(timestr,':')
    (h,m,s) = timestr.split(':')
    idstr = "%s%s%s%s%s%s" % (yr,mon,day,h,m,s)

    return (datestr, timestr, idstr)


def filenumber(file):
    "returns file number (integer nearest the file extension)"
    # copied from Spiderutils.py

    if len(file) == 0: return None
    n = getfilenumber(file)
    if n:
        return int(n)
    else:
        return None


##############################################################################################

nargs = len(sys.argv[1:])

if nargs < 4:
    print 'Usage: scipionsql2spidercoords.py sqlfile spider2micrograph_lut outputFilenums outputCoordsDocExample'
    print "\nThe file spider2micrograph_lut has the format"
    print "mic0001.dat -> /data/Micrographs/FoilHole_1108402_Data_27980690_27980691_20160102_1740_frames.spi"
    print "Assuming the local files are soft links, it can be generated using the syntax: stat -c%N mic*.dat"
    print "\nThe output file outputCoordsDocExample doesn't need to exist"
    print
    sys.exit()

#sqlfile    = 'particles-leo.sqlite'
sqlfile    = sys.argv[1]

#fileMicLUT = 'list-micrographs.txt'
fileMicLUT = sys.argv[2]
# looks like: 'mic0001.dat -> /run/media/haney.ne/SAMSUNG/EM/hTH1/january/Micrographs/Aligned/FoilHole_1108402_Data_27980690_27980691_20160102_1740_frames.spi'
# can be generated using the syntax: 'stat -c%N mic*.dat'

micList    = sys.argv[3]

#outexample = 'pkcoor_2330.dat'
outexample = sys.argv[4]

# read micrograph lookup table, of the form: 
#   mic0001.dat -> /run/media/haney.ne/SAMSUNG/EM/hTH1/january/Micrographs/Aligned/FoilHole_1108402_Data_27980690_27980691_20160102_1740_frames.spi
openMicLUT = open(fileMicLUT,'r')
listMicLUT = openMicLUT.readlines()
openMicLUT.close()

dictMicLUT = {}  # initialize

for micLine in listMicLUT:
    #strip whitespace
    micLine.strip(' \t\n\r')
    spiName = os.path.splitext(micLine.split('->')[0])[0]
    spiNum  = getfilenumber(spiName)
    #print spiNum
    epuName = os.path.basename(micLine.split('->')[1]).strip('\n')
    dictMicLUT[epuName] = spiNum

#sys.exit()

# open SQL file
conn = sqlite3.connect(sqlfile)
cur = conn.cursor()

# read column information from SQL file
cur.execute("select * from Classes")
classes = cur.fetchall()

#for classname in classes:
    #print classname

offset = 4  # first 5 columns are not listed in SQL classes

# get column numbers
micNameColumn = offset + [classinfo[1] for classinfo in classes].index('_coordinate._micName')
xcoordColumn  = offset + [classinfo[1] for classinfo in classes].index('_coordinate._x')
ycoordColumn  = offset + [classinfo[1] for classinfo in classes].index('_coordinate._y')


# read particle data from SQL file
cur.execute("select * from Objects")
table = cur.fetchall()

# initialize
prevMicNum = ''
coordsDict = {}
gloCounter = 0
listMics   = []

for particleInfo in table:
    #print particleInfo
    currMicName = particleInfo[micNameColumn]
    currMicNum  = dictMicLUT[currMicName]

    if currMicNum != prevMicNum:
        # check if there is already a dictionary entry.
        # ideally, particles are sorted already.
        if coordsDict.get(currMicNum):
            print "WARNING: Found existing particles in micrograph %s -- database probably not sorted" % currMicNum
        
        # if new micrograph
        else:
            # initialize
            coordsDict[currMicNum] = []
            listMics.append(currMicNum)
            
        #print "Now working on micrograph %s" % currMicNum
        prevMicNum = currMicNum
    
    gloCounter += 1
    coordsDict[currMicNum].append([particleInfo[xcoordColumn],particleInfo[ycoordColumn]])
    
    #if gloCounter == 53: break

#print coordsDict[currMicNum], type(coordsDict[currMicNum])

gloCounter = 0

# get micrographs
#print coordsDict.keys()
#listMics = coordsDict.keys()

for micNum in listMics:
    outtemplate = name2template(outexample)
    outname = template2filename(outtemplate,micNum)
    micPartNum = 0  # initialize
    
    with open(outname, 'w') as f:
        hdr = makeDocfileHeader(outname)
        f.write(hdr)
        
        chdr = " ;         X           Y\n"
        f.write(chdr)
        
        for coordinates in coordsDict[micNum]:
            micPartNum += 1
            line = str(micPartNum) + '  2'
            for item in coordinates:
                line += " %11g" % float(item)
            line += '\n'
            f.write(line)

    print "Wrote %s particles to micrograph %s" % (micPartNum, outname)
    gloCounter += micPartNum

print
backup(micList)
micCounter = 0

with open(micList, 'w') as filenums:
    hdr = makeDocfileHeader(micList)
    for micNum in listMics:
        micCounter += 1
        line = str(micCounter) + '  1'
        line += " %11g" % filenumber(micNum) + '\n'
        filenums.write(line)
    
print "\nFound %s particles in %s micrographs and wrote to %s" % (gloCounter, micCounter, micList)
