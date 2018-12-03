# Spider Python Library: Spiderscript.py
# Copyright (C) 2006, 2010 Health Research Inc.
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
import os, sys
from inspect import currentframe, getouterframes, getframeinfo
from linecache import getline
import re

from Spiderutils import fileReadLines, makeSpiderFilename, readdoc

###################################################################
#
# Batch file utilities
		
# get data extension (3 letters with dot) from command line or user.
# Exits program if a proper extension is not obtained.
def extFromCommand(args):
    " input must be sys.argv array with python program as first element "
    if len(args) < 2:
        ext = raw_input("data extension: ")
        if not ext:
            return ""
    else:
        ext = args[1]

    if len(ext) == 3 and ext[0] != '.':
        next =  '.' + ext
    elif len(ext) == 4 and ext[0] == '.':
        next = ext
    else:
        next = ""
        
    if next == "":
        print "extFromCommand: improper file extension. EXITING"
        sys.exit()
    else:
        return next
    
re_section = re.compile('# *--+ *[ \w]+ *--+')   # patterns: "# --- word ---"

_end_header = "END BATCH HEADER"

# hdrdict[section_heading] = [ [var, value, comment], [var, val, cmt], ... ]
def mkHeaderDictionary(B):
    "read python batch file header, returns a dictionary"
    hdrdict = {}
    title = None
    for line in B:
        line = line.strip()
        if line.upper().find(_end_header) > -1:
            return hdrdict
        if re_section.match(line):
            s = line.replace('#','')
            title = s.replace('-','').strip()
            hdrdict[title] = []
        elif line.find('=') > -1 and line.find('#') > -1 and line[0] != '#':
            # line has an assignment and a comment
            d = line.split('=')
            var = d[0].strip()
            vc = d[1].split('#')
            val = vc[0].strip()
            comment = vc[1].strip()
            asgn = [var, val, comment]
            if title in hdrdict:
                hdrdict[title].append(asgn)
            else:
                print "Error: line not within a section heading: %s" % line            
    return hdrdict

# Adds extension to string variables in the batch file header.
# Makes changes in the namespace of the calling module. Weird, huh?
def putExtension(ext, force=0):
    f = getouterframes(currentframe())
    caller = f[1][0]
    frameinfo = getframeinfo(caller)
    filename = frameinfo[0]
    nlines = frameinfo[1]
    vardict = caller.f_locals
    # use linecache to retrieve lines up to the call to this function
    hdr = []
    for i in range(nlines):
        hdr.append(getline(filename, i))
    vkeys = vardict.keys()
    hdrdict = mkHeaderDictionary(hdr)

    # Construct a list of variables found in the header.
    # Skip if 'directory' found in section heading or in comment.
    # Skip if already has an extension (unless force=1).
    hdrvars = []
    sections = hdrdict.keys()
    for sect in sections:
        s = sect.lower()
        if s.find('directory') > -1 or s.find('directories') > -1:
            #print "skipping %s" % sect
            continue
        lines = hdrdict[sect]
        for a in lines:
            if len(a) < 3:  # empty section?
                continue
            var = a[0]
            val = a[1]
            com = a[2]
            if type(val) != type("string"):
                #print "skipping %s" % val
                continue
            elif not force:
                # if it is a string, see if it has an extension
                name, ext2 = os.path.splitext(val)
                if ext2:
                    #print  "skipping %s" % val
                    continue
            if com.find('directory') > -1 or com.find('directories') > -1:
                #print "skipping %s" % var
                continue
            hdrvars.append(var)

    #print hdrvars
    if len(hdrvars) == 0:
        return
    # change the values of the variables in vardict from calling frame
    for hv in hdrvars:
        if hv in vkeys:
            val = vardict[hv]
            #print "%s: %s" % (str(hv), str(val))
            if type(val) == type("string"):
                newval = val + ext
                vardict[hv] = newval
       
def useFilenums(filenums, template):
    "returns filename list, with filenums substituted into template asterisks"
    " filenums: list of filenumbers   "
    " template: e.g., img****.dat, SPIDER filename template w/ asterisks "
    files = []
    if type(filenums) == type("string"):
        if os.path.exists(filenums):
            filenums = readdoc(filenums)
    for i in filenums:
        files.append(makeSpiderFilename(template, i))
    return files


# ######################################################################
# Parameter doc file utilities:
#      getParameters(parmfile)  a function that returns a dictionary.
#      P = Parameters(parmfile)  a class that returns an object.

# ---------------------------------------------------------------------
#
# getParameters: reads a parameter doc file (1 data column) and returns a
#                parameter dictionary, indexed by keys (integers).
# Each entry is a list with 2,3, or 4 elements, depending on what it finds:
# P[1] = [1, 0.0]                   # no comment was found
# P[1] = [1, 0.0, 'zip flag']       # comment found
# P[1] = [1, 0.0, 'zip flag', ['do not unzip', 'unzip']]  # comment has list
def getParameters(parmfile):
    " returns parameter dictionary, or None "
    B = fileReadLines(parmfile)
    nlines = len(B)
    PD = {}
    ptr = 0

    while ptr < nlines:
        line = B[ptr]
        if line.strip() == "":
            ptr += 1
            continue

        if isParmDataLine(line):
            d = line.split()
            key = int(float(d[0]))
            val = float(d[2])
            # is comment on the same line?
            if len(d) > 3:
                comment = ""
                for item in d[3:]:
                    comment += item + " "
                if comment[0] == ';':
                    comment = comment[1:].strip()
                comment, clist = findBrackets(comment)
                if clist == None:
                    PD[key] = [key, val, comment]
                else:
                    PD[key] = [key, val, comment, clist]
            else:
                # is comment on previous line?
                if ptr != 0:
                    pline = B[ptr-1]
                    if pline[:4] == ' ; /':
                        comment = pline[4:].strip()
                        comment, clist = findBrackets(comment)
                        if clist == None:
                            PD[key] = [key, val, comment]
                        else:
                            PD[key] = [key, val, comment, clist]
                    else:
                        # comment not found
                        print "WARNING: getParameters: no comment found. " \
                              "Is %s a valid parameter file?" % parmfile
                        PD[key] = [key, val]
        ptr += 1

    if len(PD) == 0: return None
    else: return PD
    
def isParmDataLine(line):
    "line must be of the form: 'key 1 value' (i.e., int 1 float)"
    d = line.split()
    if len(d) < 3:
        return 0
    try:
        key = int(float(d[0]))
        num = int(d[1])
        if num != 1:
            return 0
        val = float(d[2])
        return 1
    except:
        return 0

def findBrackets(s):
    " given a comment with brackets : zip flag ['do not unzip', ' unzip'] "
    " returns (comment, list)   or (comment, None) if no brackets "
    s = s.strip()
    # look for end bracket
    if s[-1] != ']': return (s, None)
    lb = s.rfind('[')
    if lb < 0: return (s, None)
    
    comment = s[:lb].strip()
    dlist = s[lb:]

    # dlist must have at least 2 items separated by comma
    if dlist.find(',') < -1: return (s, None)
    choices = None
    
    # string "['thing1', 'thing2']" must be converted to python list.
    # If they're already quoted, just eval the string
    try:
        choices = eval(dlist)
    except:
        # no quotes? convert to list
        dlist = dlist[1:-1]  # remove brackets
        d = dlist.split(',')
        choices = []
        for item in d:
            choices.append(item)

    if comment[-1] == ',':
        comment = comment[:-1].strip()

    return (comment, choices)

# ##################################################################
# A class for the Spider parameters doc file
# NB, use 'wavelength' instead of 'lambda' since the
#     latter is a reserved word in Python.
#
# Usage: P = Parameters('../params.hcc')
# then to access, use P.key
# >>> P.kv
# 300.0

class Parameters:
    def __init__(self, docfile, usecomments=0):
        PD = getParameters(docfile)
        if PD == None: return

        self.D = {}
        self.attribs = []  # list of attributes

        if not usecomments:
            # standard Spider single particle reconstruction parameter file
            self.K = {1:'zip',
                      2:'scan_format',
                      3:'width',
                      4:'height',
                      5: "pixel_size",
                      6: 'kv',
                      7: 'Cs',
                      8: 'source_size',
                      9 : 'defocus_spread',
                      10: 'astigmatism',
                      11: 'azimuth',
                      12: 'acr',
                      13: 'Gaussian_env',
                      14: 'wavelength',  #'Lambda',
                      15: 'max_spat_freq',
                      16: 'decimation',
                      17: 'window_size',
                      18: 'particle_size',
                      19: 'magnification',
                      20: 'scanning_res'}

            keys = self.K.keys()
            for key in keys:
                if key in PD:
                    value = PD[key][1]
                    self.D[key] = value
                    estr = "self.%s = %f" % (self.K[key], value)
                    exec(estr)
                    self.attribs.append(self.K[key])
                    
        elif usecomments:
            # use the comment as the attribute. USES FIRST WORD IN COMMENT
            self.K = {}
            keys = PD.keys()
            for key in keys:
                d = PD[key]
                if len(d) < 3:
                    print "ERROR: key %d has no associated comment"
                    continue
                value = d[1]
                self.D[key] = value
                comment = d[2].split()[0]  # get the first word
                if comment == "":
                    print "ERROR: key %d has no associated comment"
                    continue
                if hasattr(self, comment):
                    print "ERROR: comment for key %d already in use" % key
                    continue
                if comment == 'lambda':
                    print "RESERVED WORD WARNING: using 'wavelength' instead of 'lambda'"
                    comment = 'wavelength'
                if comment.find('.') > -1:
                    comment = comment.replace('.','_')

                estr = "self.%s = %f" % (comment, value)
                try:
                    exec(estr)
                    self.attribs.append(comment)
                except:
                    print "ERROR: unable to evaluate %s from %s" % (comment, d[2])

    
