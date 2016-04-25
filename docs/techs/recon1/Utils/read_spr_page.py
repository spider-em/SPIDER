#!/usr/bin/env python

#altered Dec 09 al for migration to /usr8

"""
Read the Single Particle Reconstruction web page:
 /usr8/spider/docs/techs/recon1/mr1.html

Generate an object containing the subdirectories and their batch files.
Expects each section to begin with:
"These procedures should be run in the <i>Directory</i> directory."
Subdirectories need to have the format:
"Subdirectory(ies): <i>dir1, dir2</i>"

Loads anything in a link in the 'Procs' directory.
---
Usage:
D = read_spr_page.getproject()
returns a dictionary D where
   D.keys() = [project directories]
   D[dir] = ( [procs], [subdirs] )
   A special key 'dirlist' has an ordered list of directories.
"""
import os,string,sys
import re

re_subdirs = re.compile("[sS]ubdirector(y|ies) *: *<i>.+</i>")

re_spire = re.compile("<!-- *SPIRE")

# source_dir : Location of the procedure files
source_dir = "/usr8/spider/docs/techs/recon1"

default_page = os.path.join(source_dir, "../Docs/mr1.html")

topLevelDir = "_top_level_project"

batdir = 'Procs'
batext = ['.spi', '.pam', '.py']

def stripTags(s):
    new = ""
    n = len(s)
    intag = 0
    for i in range(n):
        char = s[i]
        if char == '>':
            if not intag: # Then must be end tag from previous line
                new = ""
            else:
                intag = 0
        elif char == '<':
            intag = 1
        elif not intag:
            new += char
    return string.strip(new)

def get_ahref(line):
    ' Return contents of <a href = "contents"> '
    tag1 = "<a href"
    tag2 = ">"
    a = string.find(string.lower(line),tag1)
    if a > -1:
        line = line[a+len(tag1):]
        b = string.find(line,">")
        if b > 0:
            line = line[:b]
    # remove quotes
    a = string.find(line,'"')
    b = string.rfind(line,'"')
    if a > -1 and b > 0:
        line = line[a+1:b]
    return line

def getSubdirectories(line):
    " gets x,y,z from '<!-- SPIRE directory=a subdirectories=x,y,z -->' "
    subs = "subdirectories"
    if string.find(line,subs) < 0:
        return []
    
    a = string.find(line,subs)
    line = line[a+len(subs):]
    line = string.replace(line,"=","")
    line = string.replace(line,"-->","")
    s = string.split(line,",")
    ss = []
    for sub in s:
        ss.append(string.strip(sub))
    return ss

def getDirectory(line):
    " gets xxx from '<!-- SPIRE directory=xxx ...-->' "
    a = string.find(line,"=")
    if a > -1:
        line = line[a+1:]
        d =  string.split(line)
        return d[0]
    else:
        return ""
            

def get_batfile(line):
    ' get x.bat from a line of text <a href="./Procs/x.bat">'
    batfile = get_ahref(line)
    return os.path.split(batfile)[-1]

def readSPRpage(filename=None, dirobj=None):
    " returns dictionary w/ D[dirname] = list of batch files"
    if filename == None:
        filename = default_page
    fp = open(filename,'r')
    B = fp.readlines()
    fp.close()

    # Get batch files for the top level project directory
    dirstring = "<!-- SPIRE "
    endstring = '<!-- SPIRE end -->'  # don't bother after this point
    
    if dirobj == None:
        D = {}
        dirlist = [] # for retaining the order of the directories
    else:
        D = dirobj
        if D.has_key('dirlist'):
            dirlist = D['dirlist']
        else:
            dirlist = []
            D['dirlist'] = dirlist

    lendstr = len(dirstring)

    i = 0
    n = len(B)
    s = B[i]
    procs = []; oldprocs = []
    subs=[]; oldsubs = []
    dirname = ""

    while string.find(s,endstring) < 0 and i < n-1:
        a = re_spire.search(s)

        # Found a new directory heading
        if a:
            while string.find(s, "-->") < 0: # Keep reading lines til end comment
                i += 1
                s = s + B[i]
            newdirname = getDirectory(s)
            if newdirname == "":
                print " read_spr_page.py: unable to get directory from %s" % str(d)
                break

            newsubs = getSubdirectories(s)
            
            if dirname != "":
                if D.has_key(dirname):
                    oldprocs, oldsubs = D[dirname]
                else:
                    oldprocs, oldsubs = [],[]

                    procs = oldprocs + procs
                    subs = oldsubs + subs
                D[dirname] = (procs, subs)
                if dirname not in dirlist:
                    dirlist.append(dirname)
            dirname = newdirname
            subs = newsubs
            procs = []

        # Look for procedures 
        elif string.find(s,batdir) > -1:

	    if string.find(s, 'class="project"') > -1:  # Added al Oct 2015
            	proc = get_batfile(s)
            	fname,ext = os.path.splitext(proc)
            	if ext in batext and proc not in procs:
			procs.append(proc)
                
        i = i + 1
        s = B[i]
        # end while

    if dirname != "":
        if D.has_key(dirname):
            oldprocs, oldsubs = D[dirname]
            for p in oldprocs:
                if p not in procs:
                    procs.append(p)
            for s in oldsubs:
                if s not in subs:
                    subs.append(s)

        D[dirname] = (procs, subs)   # finish up the last set
        if dirname not in dirlist:
            dirlist.append(dirname)
 
    D['dirlist'] = dirlist
    return D

def getproject(webpages=None):
    if webpages == None:
        webpages = [os.path.join(source_dir,"Docs/mr1.html")]
    elif type(webpages) == type("string"):
        webpages = [webpages]

    D = {}

    for webpage in webpages:
        webdir, webfile = os.path.split(webpage)
        os.chdir(webdir)

        D = readSPRpage(webpage, D)

    return D

if __name__ == '__main__':

    D = getproject()

    dirs = D['dirlist']
    print dirs
    for dir in dirs:
        print dir
        procs, subdirs = D[dir]
        for proc in procs:
            print "     " + proc
        for sub in subdirs:
            print "     " + sub + "/"
