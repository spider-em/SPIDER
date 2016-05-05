#!/usr/bin/env python

# setup.py
#
# Installs Spire and Spider packages in Python's site-packages directory.
#
# spire and other executables are put in the local directory: spire-?.?.?/bin

import os,sys, glob
from commands import getoutput
from shutil import copy, copytree, rmtree
import findprog, testinstall

SPIRE_DIR = os.getcwd()
SPIRE_INSTALL_DIR = SPIRE_DIR

BIN_DIR = os.path.join(SPIRE_INSTALL_DIR, "bin")
if not os.path.exists(BIN_DIR):
    os.mkdir(BIN_DIR,0775)

GUITOOLS_DIR = os.path.join(SPIRE_DIR, "guitools")
CONFIG_DIR   = os.path.join(SPIRE_DIR, "config")
BATCH_DIR    = os.path.join(SPIRE_DIR, "batchfiles")

def checkFileAccess(f):
    if os.access(f,os.R_OK) and os.access(f,os.W_OK):
        return 1
    else:
        return 0

def mkdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname,0775)

def getPythonPath():
    for p in sys.path:
        base = os.path.basename(p)
        if base == 'site-packages':
            return p
    return ""

def getPILpath():
    for p in sys.path:
        base = os.path.basename(p)
        if base == 'PIL':
            return p
    return ""

def newsrc2config(oldconfig, newconfig, srclist):
    fp = open(oldconfig, 'r')
    B = fp.readlines()
    fp.close()

    fp = open(newconfig,'w')
    for line in B:
        if line.find("<Locations>") > -1:
            fp.write("<Locations>\n")
            for item in srclist:
                fp.write(item)
        elif line.find("<location>") > -1: # Removes pre-existing locations
            continue
        else:
            fp.write(line)
    fp.close()

def copy2bin(dirname, dest='bin'):
    " copy all *.py files to spire/bin "
    if not os.path.exists(dirname):
        return
    progs = os.listdir(dirname)
    for file in progs:
        path,ext = os.path.splitext(file)
        if ext == '.py':
            pp = os.path.join(dirname, file)
            copy(pp, dest)

# ##################################################################
# test if installation has Spire prerequisites

res = testinstall.testinstall()
if res == 'ok':
    print " All required components are present"
else:
    print " Some components are missing from this version of Python:"
    if 'tk' in res:
        print " Tcl/Tk 8.4 or higher"
    if 'blt' in res:
        print " Blt"
    if 'python' in res:
        print " Python version must be 2.3 or higher"
    if 'pmw' in res:
        print " Pmw (Python Megawidgets)"
    if 'image' in res:
        print " PIL (Python Imaging Library)"
    if 'numpy' in res:
        print "numpy (Numerical Python)"
    if 'gnuplot' in res:
        print " Gnuplot.py"
    if 'db' in res:
        print " Berkeley Database"
        
    print " Try the full install of Spire with its own copy of Python and all required packages."
    #sys.exit()

#################################################################
# Get Spider, jweb commands
print "\n ----------------------------------------------"
spidercmd = ""
pathlist = findprog.findpath('spider')

for cmd in pathlist:
    if findprog.testSpider(cmd):
        print " %s can successfuly run Spider" % (cmd)
        spidercmd = cmd
        break
    
if spidercmd != "":
    print "\n Is this the correct command to run spider?: %s" % spidercmd
print " Please enter the full path of the command to run SPIDER on this system,"
#print "e.g., the results of the Unix commands 'which spider' or 'alias spider'"
if spidercmd != "":
    print " or just hit Enter if the above command is acceptable"
cmd = raw_input(" SPIDER command: ")
if cmd != "":
    spidercmd = cmd
elif cmd == "" and spidercmd != "":
    cmd = spidercmd

if findprog.testSpider(cmd):
    print " %s can successfuly run Spider" % (cmd)
        
# jweb ----------------------------
jwebcmd = "$JWEB_DIR/j2re1.4.2_06/bin/java -Xmx512m -Djava.util.prefs.syncInterval=2000000 web/StartWeb"
"""
os.system(cmd)             ; system(cmd) can have variables like $JWEB_DIR, but
os.path.exists($JWEB_DIR)  ; os.path.exists() cannot resolve them.

"""
j = findprog.substituteVars(jwebcmd.split()[0])
if not os.path.exists(j):
    jlist = findprog.findpath('jweb')
    if len(jlist) > 0:
        print " The following commands all seem to be able to run jweb:"
        for item in jlist:
            print item
            if item.find('JWEB_DIR') > -1:
                jwebcmd = item
if jwebcmd != "":
    print "\n Is this the correct command to run jweb?: %s" % jwebcmd
print " Please enter the full path of the command to run jweb on this system,"
#print "e.g., the results of the Unix commands 'which jweb' or 'alias jweb'"
if jwebcmd != "":
    print " or just hit Enter if the above command is acceptable"
cmd = raw_input(" JWeb command: ")
if cmd != "":
    jwebcmd = cmd
elif cmd == "" and jwebcmd != "":
    cmd = jwebcmd
if jwebcmd != "" and jwebcmd.find('SPIDUI') < 0:
    jwebcmd += " SPIDUI &"

#print jwebcmd

###################################################################
# put the <location> tag in the configuration files
print "\n Writing <location> tags in configuration files"
os.chdir(CONFIG_DIR)

# for each config file XYZ.xml, if it finds a corresponding directory 
# named ../batchfiles/XYZ, it will write that location in the config file.

configfiles = glob.glob("*.xml")
for config in configfiles:
    basename, ext = os.path.splitext(config)
    batchfiles = os.path.join(BATCH_DIR, basename)
    if basename == "simple2":    # a special case
        batchfiles = batchfiles[:-1]
    if os.path.exists(batchfiles):
        print config
        batchsrc = ["<location>" + batchfiles + "</location>\n"]
        newsrc2config(config, config, batchsrc)
os.chdir(SPIRE_DIR)
                                      
###################################################################
# write values to Spire/LocalVars.py
localvars = 'Spire/LocalVars.py'
if os.path.exists(localvars):
    os.remove(localvars)
    
date = getoutput("date")
localtxt = '"""\n'
localtxt += " LocalVars.py : written %s\n" % date
localtxt += '"""\n\n'
localtxt += "spidercmd = '%s'\n" % spidercmd
localtxt += "jwebcmd = '%s'\n" % jwebcmd
configdir = os.path.join(SPIRE_INSTALL_DIR, 'config')
localtxt += "defaultconfigdir = '%s'\n" % configdir
host = os.uname()[1]   #os.environ['HOST']
localtxt += "hostlist = ['%s']\n" % host
localtxt += "ssh_hosts = []" + "\n"

fp = open(localvars,'w')
fp.write(localtxt)
fp.close()
print localtxt
print " Written to: %s" % localvars

###################################################################
# check if there is a .spire preference file
helpurl = 'helpURL: http://www.wadsworth.org/spider_doc/spider/docs/'

prefs_file = os.path.join(os.environ['HOME'],".spire")

if os.path.exists(prefs_file):
    fp = open(prefs_file,'r')
    B = fp.readlines()
    fp.close()
    P = ""
    for line in B:
        if line.find('helpURL:') > -1:
            P += helpurl + "\n"
        elif line.find('spider:') > -1:
            P += "spider: " + spidercmd + "\n"
        elif line.find('configdir:') > -1:
            P += 'configdir: ' + configdir + "\n"
        else:
            P += line

    fp = open(prefs_file,'w')
    fp.writelines(P)
    fp.close()

# ##################################################################
# Copy Spire files to Python site-packages

print "\n === Installing Spire ==="
os.chmod('Spire/spire.py', 0775) # make spire.py executable
sitepackages = getPythonPath()
print " site-packages location: %s" % sitepackages

if sitepackages == "":
    print " Please copy the Spire directory to the python/site-packages directory"
else:
    if not checkFileAccess(sitepackages):
        print " Installation requires write permission for %s" % sitepackages
        print " Please copy the Spire directory to the python/site-packages directory"
    else: 
        # copy the Spire package
        spirepath = os.path.join(sitepackages,"Spire")
        if os.path.exists(spirepath):
            print " Spire package already installed."
            r = raw_input(" Delete the old Spire and install new version? (y/n): ")
            res = r.strip()
            if len(res) == 0 or res[0] == 'y' or res[0] == 'Y':
                rmtree(spirepath, ignore_errors=1)
            else:
                sys.exit()

        src = 'Spire'
        print " Copying %s to %s" % (src, sitepackages)
        copytree(src, spirepath)
        
        # ------- Copy the Spider package -------
        spirepath = os.path.join(sitepackages,"Spider")
        if os.path.exists(spirepath):
            #print "Spider package already installed."
            #r = raw_input("delete the old Spider and install new version? (y/n): ")
            #res = r.strip()
            if len(res) == 0 or res[0] == 'y' or res[0] == 'Y':
                rmtree(spirepath, ignore_errors=1)
            else:
                sys.exit()

        src = 'Spider'
        print " Copying: %s to: %s" % (src, sitepackages)
        copytree(src, spirepath)
     
        # This is done by mkapps.py also
        # pyplot.py and montage.py are also put in the Spider directory
        #for prog in ["montage.py", "pyplot.py"]:
        #    src = os.path.join(GUITOOLS_DIR,prog)
        #    if os.path.exists(src):
        #        print " Copying: %s to: %s" % (src, spirepath)
        #        copy(src, spirepath)
        #    else:
        #        print " Unable to find: %s For: %s" % (src, spirepath)

        # ------- copy SpiderImagePlugin.py to PIL -------
        pil = os.path.join(sitepackages, "PIL")
        if os.path.exists("SpiderImagePlugin.py") and os.path.exists(pil):
            copy("SpiderImagePlugin.py", pil)
        
###############################################################
# Put binaries in bin directory, add to path in .cshrc

# Copy executables to bin
copy2bin(GUITOOLS_DIR, dest=BIN_DIR)

copy('Spire/spire.py', os.path.join(BIN_DIR,'spire'))

path = os.environ['PATH']
if BIN_DIR in path:
    print " Setup complete"    # you're done
    sys.exit()

print " Executables in %s must be in your path." % (BIN_DIR)    
print " You can either"
print " 1) add %s to your path, or" % BIN_DIR
print " 2) copy the files in %s to another directory located on your path" % BIN_DIR

shell = os.environ['SHELL']
if shell[-3:] == 'csh':
    cshrc = os.path.expanduser("~" + os.environ['USER'] + "/.cshrc")
    print #cshrc
    if os.path.exists(cshrc):
        print " The following line can be added to your .cshrc file:"
        newpath = " set path =(%s $path)" % BIN_DIR
        print newpath
        yn = raw_input(" Shall I add it for you? (y/n): ")
        yn = yn.lower()
        if yn == 'y' or yn == 'yes' or yn == 'Y':
            old_cshrc = cshrc + ".old"
            copy(cshrc, old_cshrc)
            fp = open(cshrc,"a")
            fp.write("\n# added by Spire installation %s\n" % date)
            fp.write(newpath + "\n")
            fp.close()
            print " Installation complete. Type:"
            print " source %s" % cshrc
            print " to be able to run spire"
            
elif shell[-4:] == 'bash':
    bashrc = os.path.expanduser("~" + os.environ['USER'] + "/.bashrc")
    print #bashrc
    if os.path.exists(bashrc):
        print " The following line can be added to your .bashrc file:"
        newpath = "export PATH=%s:$PATH" % BIN_DIR
        print newpath
        yn = raw_input(" Shall I add it for you? (y/n): ")
        yn = yn.lower()
        if yn == 'y' or yn == 'yes':
            old_bashrc = bashrc + ".old"
            copy(bashrc, old_bashrc)
            fp = open(bashrc,"a")
            fp.write("\n# added by Spire installation %s \n" % date)
            fp.write(newpath + "\n")
            fp.close()
            print " Installation complete. Type:"
            print " source %s" % bashrc
            print " to be able to run spire"
