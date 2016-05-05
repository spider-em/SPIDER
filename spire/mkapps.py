
#  PURPOSE: Find the Python executable for Spire, in spire_linux-1.5.3/bin/python.
#           Create scripts for the *.py programs in the spire bin directory.
#           Move the *.py programs into the scripts directory.
#
# USAGE:    mkapps.py BIN_DIR_FOR_EXECS TARGET_DIR_FOR_PY_PROGS

import os, sys
from shutil import copy


apps = ['scatter.py', 'pyplot.py', 'ctfdemo.py', 'xplor.py', 'montagefromdoc.py',
        'ctfgroup.py', 'ctfcircle.py', 'ctfmatch.py', 'classavg.py', 'montage.py',
        'qview.py', 'verifybyview.py', 'backup.py', 'emancoords2spiderdoc.py', 
        'mkfilenums.py', 'spiconvert.py', 'binarytree.py', 'viewstack.py', 'xmippsel2spiderdoc.py']

_TARGETDIR = 'scripts'

SPIRE_INSTALL_DIR = os.getcwd()
GUITOOLS_DIR      = os.path.join(SPIRE_INSTALL_DIR, "guitools")
#print " GUITOOLS_DIR: %s" % GUITOOLS_DIR

BIN_DIR = os.path.join(SPIRE_INSTALL_DIR, "bin")
if not os.path.exists(BIN_DIR):
    os.mkdir(BIN_DIR,0775)

# Find the python executable for spire

numargs = len(sys.argv[1:])

if numargs == 0:
    d = os.getcwd()    # This script is either in SPIREHOME or SPIREHOME/scripts
    path,basename = os.path.split(d)
    if basename.find("spire") > -1:
        spiredir  = basename
        targetdir = os.path.join(spiredir,_TARGETDIR)
    else:
        path, basename = os.path.split(path)
        if basename.find("spire") > -1:
            spiredir  = basename
            targetdir = os.path.join(spiredir,_TARGETDIR)
        else:
            print " Usage: mkapps.py spire_dir"
            sys.exit(1)
    Pythonprog = os.path.join(spiredir,"bin/python")

elif numargs == 1:   # 1st arg is the bin directory
    bindir = sys.argv[1]
    spiredir,basename = os.path.split(bindir)
    targetdir         = os.path.join(spiredir,_TARGETDIR)
    Pythonprog        = os.path.join(bindir, "python")

else:
    bindir = sys.argv[1]
    spiredir,basename = os.path.split(bindir)
    targetdir  = sys.argv[2]
    Pythonprog = os.path.join(bindir, "python")
            
if not os.path.exists(Pythonprog):
    print " Unable to find executable python in: %s" % Pythonprog
    sys.exit(1)
    
os.chdir(GUITOOLS_DIR)
targetdir = os.path.join(spiredir,BIN_DIR)

for app in apps:
    appscript, ext = os.path.splitext(app)
    if ext != ".py" or not os.path.exists(app):
        continue
    app_path = os.path.join(targetdir,app)
    txt      = "#!/bin/sh\n\n"
    txt     += 'exec "%s" "%s" "$@"' % (Pythonprog, app_path)
    appscript = os.path.join(targetdir,appscript) 
    
    print " Creating script: %s" % appscript
    
    try:
        fp = open(appscript, 'w')
        fp.write(txt)
        fp.close()
        os.chmod(appscript, 0775)
	
        #print " Try to copy: %s to: %s" % (app, targetdir)
        copy(app,targetdir)
        #print " Copying: %s to: %s" % (app, targetdir)
        #if os.path.exists(os.path.join(targetdir,app)):
        #   os.remove(app)
    except:
        print " Unable to create: %s script" % appscript

    
