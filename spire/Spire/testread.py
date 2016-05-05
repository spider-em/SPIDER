import sys
import string
import time
import os
import re
from Spider.Spiderutils import fileReadLines

TypeDict = { 'R': 'Image',
             'R3': 'Volume',
             'E2': '2D Fourier',
             'O2': '2D Fourier',
             'E3': '3D Fourier',
             'O3': '3D Fourier'}

StackTypes = {'E2': 'Fourier Stack',
              'O2': 'Fourier Stack',
              'E3': 'Fourier Stack',
              'O3': 'Fourier Stack',
              'S2': '2D Stack',
              'S3': '3D Stack',
              'I2': '2D Stack',
              'I3': '3D Stack'}

def getDataTuple(s, isStack=0):
    " s is a line for binary files from spireout, of the form:  "
    " (R ) 75 75 CREATED 13-FEB-2007 AT 13:04:16  N HEADER BYTES: 1200 "
    " returns [type, size, date, time] "
    xtype = s[1:3].strip()
    if isStack and xtype in StackTypes:
        type = StackTypes[xtype]
    elif xtype in TypeDict:
        type = TypeDict[xtype]
    else:
        type = 'Unknown'
        print "spiderResult.getDataTuple: type %s not recognized" % xtype

    dim = 2
    if len(xtype) > 1:
        dim = int(xtype[1])
        
    t = s[4:].split()

    if dim == 2:
        size = (t[0],t[1])
    elif dim == 3:
        size = (t[0],t[1],t[2])
    else:
        size = ('size not found')
        print "getDataTuple: size not found"
        print s

    # assumes end of string has following format :
    # " .... CREATED 09-FEB-2007 AT 09:57:13  N HEADER BYTES:   1184"
    date = t[-7]
    time = t[-5]

    return [type, size, date, time]

def printDD(d):
    keys = d.keys()
    keys.sort()
    for k in keys:
        print k
        #print "%s: %s" % (k, str(d[k]))

def processResult(B,datext):
    if datext[0] != '.':
        datext = '.' + datext
    nlines = len(B)
    BinFileList = []
    DocFileList = []
    bflist = []  # ordered list of filenames, contains names of deleted files
    dclist = []
    ptr = 0
    StackDict = {}
    BinDict = {}   # holds file tuple info, DELETE files are removed from it

    while ptr < nlines:
	line = B[ptr].strip()

	# ------- look for new binary files (with 'N HEADER' ) -------
	# spire-1.5.1 : no longer checks if the binary file is already
	# in the dictionary. ANY time a file is rewritten (by CP or any
	# command) the N HEADER line appears. spire now keeps updating
	# dictionary every time a file appears, so that the project will
	# contain the last version. A separate dictionary is kept for
	# stacks, because the vast majority of stack calls are just
	# inserting stack elements.
	if line.find('N HEADER') > -1:
            fname = ""
	    previous_line = B[ptr-1].strip()

	    istk = previous_line.find('@')    # check if it's a stack
	    if istk > -1:
                isStack = 1
                fname = previous_line[:istk]
            else:
                isStack = 0
                fname = previous_line
                
            if fname[-4:] != datext:
                fname += datext
	    
            if not isStack: 
                # getDataTuple returns [type, size, date, time]
                newdata = getDataTuple(line, isStack)
                BinDict[fname] = [fname] + newdata
                bflist.append(fname)
                
	    elif isStack:
                if fname not in StackDict:
                    newdata = getDataTuple(line, isStack)
                    StackDict[fname] = "" #######[fname] + newdata
                    BinDict[fname] = [fname] + newdata
                    bflist.append(fname)
                else:
                    fname = ""  # skip if the stack is already in StackDict
                    
            else:
                fname = ""
                    
            if fname == "":
                ptr = ptr + 1
                continue

	# ---------- look for new doc files -----------
	elif line.find('NEW DOC FILE') > -1:
	    t = line.split()
            # "09-FEB-2007 AT 09:57:13    OPENED NEW DOC FILE: tmp/doc001.dat"
	    docfile, docdate, doctime = t[-1],t[0],t[2]
	    if docfile[-4:] != datext: docfile += datext
	    DocFileList.append([docfile, docdate, doctime])
	    dclist.append(docfile)

	# ---------- check for deleted files ------------
	elif line.find('.DELETE FILE:') > -1:
	    delfile = line.split(':')[-1].strip()
	    if delfile[-4:] != datext: delfile += datext

	    if delfile in dclist:
		idx = dclist.index(delfile)
		dclist.remove(delfile)
		del DocFileList[idx]
		
	    elif delfile in BinDict:
                del(BinDict[delfile])
		# not necessary? => bflist.remove(delfile)
		if delfile in StackDict:
                    del(StackDict[delfile])
            

        ptr = ptr + 1
        # ---------- end: while ptr < nlines --------------------

    # now use bflist to create BinFileList
    # (files only deleted from BinDict, but only bflist retains file order
    for bf in bflist:
        if bf in BinDict:
            BinFileList.append(BinDict[bf])
            del(BinDict[bf])

    for item in BinFileList:
        print item
    print "-------------------------------------------------"
    for item in DocFileList:
        print item


# -------------------------------------------------------------------

if __name__ == '__main__':

    spireout = "spireout.spi.0"
    B = fileReadLines(spireout)
    processResult(B,'dat')
