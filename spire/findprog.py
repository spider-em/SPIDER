#!/usr/bin/env python

import sys, os, re
from commands import getoutput

def testSpider(spiderpath):
    file = 'test6637'
    ext = ".bat"
    filename = file + ext
    cmd = 'echo "en d" > %s' % filename
    os.system(cmd)
    spicmd = "%s bat/dat @%s" % (spiderpath, file)
    #print spicmd

    success = 0
    output = getoutput(spicmd)
    if output.find('Results file') > 0:
        success = 1
        #print output
        log = "LOG" + ext
        if os.path.exists(log):
            os.remove(log)
    os.remove(filename)
    return success


#def substituteVars(pathstring):
    #" replaces $PATTERN with os.environ[PATTERN] "
    #s = pathstring.split('/')
    #if pathstring[0] == "/":
        #new = "/"
    #else:
        #new = ""
        
    #for d in s:       
        #a = d.find("$")
        #if a > -1:
            #var = d[a+1:]
            #if var in os.environ:
                #value = os.environ[var]
                #d = d.replace("$"+var,value)
        #new += d + "/"
        
    #return new[:-1]



def unquote(s):
    " remove quotes, if any"
    p = s.strip()
    if p[0] == "'" or p[0] == '"':
        p = p[1:]
    if p[-1] == "'" or p[-1] == '"':
        p = p[:-1]
    return p

def findwhich(program, verbose=0):
    cmd = 'which %s' % program
    res = getoutput(cmd)
    if verbose:
        print res
    wlist = res.split('\n')
    progpath = []
    for which in wlist:
        if program in which and os.path.exists(which):
            return [which]
        if which.find('aliased') > -1:
            dw = which.split('aliased to')
            for d in dw:
                if d.find(program) > -1:
                    progpath.append(unquote(d))
        notin = program + ' not in'
        noprog = "no %s in" % program
        notfound = "not found"
        if which.find(notin) > -1 or which.find(notfound) > -1 or which.find(noprog) > -1 :
            return []
        else:
            return progpath

def findpath(program, verbose=0):
    " returns a path, i.e., no command-line arguments"
    found = 0
    progpath = []
    # see if SPIDER environmental variables are set
    if program == 'spider' and 'SPIDER_LOC' in os.environ:
        spider_loc = os.path.join(os.environ['SPIDER_LOC'], program)
        if os.path.exists(spider_loc):
            if spider_loc not in progpath:
                progpath.append(spider_loc)
                
    # try the which command
    whichlist = findwhich(program)
    if len(whichlist) != 0:
        for w in whichlist:
            if not os.path.exists(w):
                if verbose: print "which: %s : not a valid path" % w
            else:
                found = 1
                if w not in progpath:
                    if w not in progpath:
                        progpath.append(w)
                    if verbose:
                        print "The 'which' command located spider executable %s" % w
            
    if len(progpath) > 0:
        return progpath

    # look for aliases in user's .cshrc, .bashrc files
    user = os.environ['USER']
    shell = os.environ['SHELL']

    if shell[-3:] == 'csh':
        sh = 'csh'
        rcfile = "~%s/.cshrc" % user
    else:
        sh = 'sh'
        rcfile = "~%s/.bashrc" % user

    rcfile = os.path.expanduser(rcfile)

    if os.path.exists(rcfile):
        if verbose: print "Looking for %s in %s" % (program, rcfile)
        fp = open(rcfile,'r')
        B = fp.readlines()
        fp.close()
        if sh == 'csh':
            for line in B:
                # look for "alias program ...
                s = line.split()
                if len(s) > 2 and s[0] == 'alias' and s[1] == program:
                    a = line.replace('alias','',1)
                    b = a.replace(program,'',1).strip()
                    p = unquote(b)
                    if verbose: print "found %s in %s" % (p, rcfile)
                    # if len == 1, then see if it's a path
                    if len(p) == 1:
                        if os.path.exists(p):
                            if p not in progpath:
                                progpath.append(p)
                                found = 1
                        else:
                            if verbose: print "...but it's not a valid path."
                    # if len > 1 then it may be a command with args
                    else:
                        progpath.append(p)
                        
                # can also try looking for .cshrc-spider
                if line.find('.cshrc-spider') > -1:
                    s = line.split()
                    for item in s:
                        if os.path.exists(item):
                            if verbose: print "checking %s" % item
                            fp = open(item,'r')
                            C = fp.readlines()
                            fp.close()
                            for cline in C:
                                cs = cline.split()
                                if len(cs) > 2 and cs[0] == 'alias' and cs[1] == program:
                                    a = cline.replace('alias','',1)
                                    b = a.replace(program,'',1).strip()
                                    p = unquote(b)
                                    if len(p) == 1:
                                        if p not in progpath:
                                            if verbose: print "found %s in %s" % (p, item)
                                            if os.path.exists(p):
                                                if p not in progpath:
                                                    progpath.append(p)
                                                found = 1
                                            else:
                                                if verbose: print "but %s not a valid path" % p
                                    else:
                                        if p not in progpath:
                                            progpath.append(p)
                    
        elif sh == 'sh':
            re_alias = re.compile('alias +%s=' % program)
            for line in B:
                if re_alias.search(line):
                    s = line.split('=')
                    prog = s[0].split()[1]
                    if prog == program:
                        p = unquote(s[1])
                        if verbose: print "found %s in %s" % (p, rcfile)
                        if os.path.exists(p):
                            if p not in progpath:
                                progpath = p
                            found = 1
                        else:
                            if verbose: print "but %s not a valid path" % p
    return progpath
        

if __name__ == '__main__':
    
    pathlist = findpath('spider', verbose=1)
    print pathlist
    if len(pathlist) == 0:
        print "spider executable not found"
        sys.exit()
    for path in pathlist:
        if testSpider(path):
            print "%s can successfuly run Spider" % (path)
    print
    # jweb ---------------------------
    
    JWEB_DIR = ""
    if 'JWEB_DIR' in os.environ:
        JWEB_DIR = os.environ['JWEB_DIR']
        
    jwebcmd = ""
    
    jw  = findpath('jweb', verbose=1)
    if len(jw) > 0:
        print "The following are possible aliases for jweb:"
        for item in jw:
            print item
            if item.find(JWEB_DIR) > -1:
                jwebcmd = item
    else:
        print "unable to find jweb."

    if jwebcmd != "":
        print "try running jweb with:"
        print jwebcmd

