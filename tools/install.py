#!/usr/bin/env python
#
# SOURCE:  spider/tools/install.py
#
# PURPOSE:  Puts executables in directory: spider/tools/bin

import            os,sys, glob
from    commands  import getoutput
from    shutil    import copy, copytree, rmtree

# 15 apps
apps = ['binarytree.py',               'montagefromdoc.py',
        'classavg.py',                 'montage.py',
        'ctfcircle.py',                'pyplot.py',
        'qview.py',
        'spiconvert.py',
        'verifybyview.py',
        'emanrctcoords2spiderdoc.py',  'viewstack.py',
        'mkapps.py',                   'xmippsel2spiderdoc.py',
        'mkfilenums.py',               'xplor.py']

# 5 apps
appsold = ['ctfdemo.py',     'ctfgroup.py',
           'scatter.py',     'ctfmatch.py', 
           'emancoords2spiderdoc.py']

tools_install_dir = os.getcwd()
bindir            = os.path.join(tools_install_dir, "bin")
pybindir          = os.path.join(tools_install_dir, "bin-python")

# Python executable for tools
Pythonprog        = os.path.join(pybindir, "python")
#print " Pythonprog: %s" % Pythonprog

if not os.path.exists(Pythonprog):
    print " Unable to find executable python: %s" % Pythonprog
    sys.exit(1)
    
# Add to path in shell startup script -------------------------------

path = os.environ['PATH']
if bindir in path:
    print " Setup complete"    # Done
    sys.exit()

print " Executables in: %s must be in your path." % (bindir)    
print " You can either"
print " 1) Add: %s to your path, or" % bindir
print " 2) Copy the files in: %s to another directory" % bindir
print "    located on your path" 

shell = os.environ['SHELL']
if shell[-3:] == 'csh':
    cshrc = os.path.expanduser("~" + os.environ['USER'] + "/.cshrc")
    print #cshrc
    if os.path.exists(cshrc):
        print " The following line can be added to your .cshrc file:"
        newpath = " set path =(%s $path)" % bindir
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
        newpath = "export PATH=%s:$PATH" % bindir
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







 
  
  
  
  
