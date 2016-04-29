#!/usr/bin/env python

# SOURCE: testinstall.py
# PURPOSE: Test if packages required for tools are present

import os,sys,shelve,whichdb

RequiredPythonVersion = 2.3
RequiredTkVersion     = 8.4

def testinstall(verbose=1):
    needed = []

    print " \n Checking for required python packages. -----------" 
    
    # Python version --------------------------------
    pyversion = sys.version.split()[0]
    try:
        pv = float(pyversion[0:3])
        if verbose: print " Found Python %s" % pyversion 
        if pv < RequiredPythonVersion:
            if verbose:
                print " These tools require Python version %s or higher." % str(RequiredPythonVersion)
            needed.append('python')
    except:
        print " Unable to get Python version"
        return ['python']

    # Check for PIL, Numeric, and Gnuplot packages
    try:
        import Image
        if verbose: print " Found Image.%s" % Image.VERSION
    except:
        needed.append('image')
    try:
        import numpy
        if verbose: print " Found numpy.%s" % numpy.__version__
    except:
        needed.append('numpy')
    try:
        import Gnuplot
        if verbose: print " Found Gnuplot.%s" % Gnuplot.__version__
    except:
        needed.append('gnuplot')
        
    # Pmw ---------------------------------------------
    pmwflag = 0
    try:
        import Pmw
        if verbose: print " Found Pmw.%s" % Pmw._version
        pmwflag = 1
    except:
        needed.append('pmw')

    # Tkinter ------------------------------------------
    tkflag = 0
    try:
        import Tkinter
        if verbose: print " Found Tkinter with Tk %s" % str(Tkinter.TkVersion)
        if Tkinter.TkVersion < RequiredTkVersion:
            print "Tools require Tcl/Tk version %s or higher." % str(RequiredTkVersion)
            needed.append('tk')
        else:
            tkflag = 1
    except:
        needed.append('tk')

    print " Need packages:" % needed
    print "  "  
    
