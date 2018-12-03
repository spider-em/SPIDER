#!/usr/bin/env python
#
# SOURCE:   spider/spire/setup.py
#
# PURPOSE:  Sets: some Spire variables for spider command, XML configuration
#           files, Spire LocalVars, and users home dir.  '.spire' files
#           Some of these may be obsolete??

import               findprog
import               os,sys, glob
from commands import getoutput
from shutil   import copy, copytree, rmtree

#################################################################

print "\n ----------------------------------------------"

SPIRE_DIR    = os.getcwd()
BIN_DIR      = os.path.join(SPIRE_DIR, "bin")
CONFIG_DIR   = os.path.join(SPIRE_DIR, "config")
BATCH_DIR    = os.path.join(SPIRE_DIR, "batchfiles")

print " SPIRE_DIR:    %s" % (SPIRE_DIR)
print " BIN_DIR:      %s" % (BIN_DIR)
print " CONFIG_DIR:   %s" % (CONFIG_DIR)
print " BATCH_DIR:    %s" % (BATCH_DIR)

#################################################################

print "\n ----------------------------------------------"

# Get python executable 

Pythonprog        = os.path.join(BIN_DIR, "python")
#print " Pythonprog: %s" % Pythonprog
if not os.path.exists(Pythonprog):
    print " Unable to find executable python: %s" % Pythonprog
    sys.exit(1)

#################################################################

# Get SPIDER command

spidercmd = ""
pathlist  = findprog.findpath('spider')

for cmd in pathlist:
    if findprog.testSpider(cmd):
        print " %s can successfully run SPIDER" % (cmd)
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
    print " %s can successfully run SPIDER" % (cmd)
   

###################################################################

# Put the <location> tag in the configuration files

print "\n Writing <location> tags in configuration files"
os.chdir(CONFIG_DIR)

# For each config file XYZ.xml, if it finds a corresponding directory 
# named ../batchfiles/XYZ, it will write that location in the config file.
# As distributed in 2018 this does NOTHING! al

configfiles = glob.glob("*.xml")
for config in configfiles:
    basename, ext = os.path.splitext(config)
    batchfiles = os.path.join(BATCH_DIR, basename)
    
    #print " batchfilesR:    %s" % (batchfiles)
    
    if basename == "simple2":    # A special case
        batchfiles = batchfiles[:-1]
    if os.path.exists(batchfiles):
        print config
        batchsrc = ["<location>" + batchfiles + "</location>\n"]
        newsrc2config(config, config, batchsrc)
	
os.chdir(SPIRE_DIR)
                                      
###################################################################

# Write some values to Spire/LocalVars.py.
# Am uncertain where these are used currently. al

localvars = 'lib/python2.5/site-packages/Spire/LocalVars.py'
if os.path.exists(localvars):
    os.remove(localvars)
    
date      = getoutput("date")
localtxt  = '"""\n'
localtxt += " LocalVars.py : written %s\n" % date
localtxt += '"""\n\n'
localtxt += "spidercmd = '%s'\n" % spidercmd
localtxt += "jwebcmd = ''\n"
configdir = os.path.join(SPIRE_DIR, 'config')
localtxt += "defaultconfigdir = '%s'\n" % configdir
host      = os.uname()[1]     #os.environ['HOST']
localtxt += "hostlist = ['%s']\n" % host
localtxt += "ssh_hosts = []" + "\n"

fp = open(localvars,'w')
fp.write(localtxt)
fp.close()
print localtxt
print " Written to: %s" % localvars

###################################################################

# Check if there is a .spire preference file in the users home directory
# If so add a SPIDER web page link.   Am uncertain how this is used
# currently. al

helpurl    = 'helpURL: http://spider.wadsworth.org/spider_doc/spider/docs/'

prefs_file = os.path.join(os.environ['HOME'],".spire")

if os.path.exists(prefs_file):
    fp = open(prefs_file,'r')
    B  = fp.readlines()
    fp.close()
    
    P  = ""
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

print " " 
print " Prefs_file: %s" % prefs_file


###################################################################

# Add to path in users shell startup script 

path = os.environ['PATH']
if BIN_DIR in path:
    print " Setup complete"    # Done
    sys.exit()

print " Executables in: %s must be in your path." % (BIN_DIR)    

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
            fp.write("\n# added by SPIDER python tools installation %s\n" % date)
            fp.write(newpath + "\n")
            fp.close()
            print " Installation complete. Type:"
            print " source %s" % cshrc
            print " to be able to run SPIDER python tools"
            
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
            fp.write("\n# added by SPIDER python tools installation %s \n" % date)
            fp.write(newpath + "\n")
            fp.close()
            print " Installation complete. Type:"
            print " source %s" % bashrc
            print " to be able to run SPIDER python tools"

print " " 
print " Installation complete" 
print " " 






 
  
  
  
  
