#!/usr/bin/env python

# Test if packages required for Spire are present

import os,sys,shelve,whichdb

RequiredPythonVersion = 2.3
RequiredTkVersion = 8.4

# a class to test the database
class AC:
    var1 = 42
    def acprint(self):
        print "hello"

def has_dbhash():
    " returns 1 if has a dbhash database "
    dbflag = 0
    filename = '_mytmpdata63636'

    db = shelve.open(filename)

    # put in some test data (incl. a class object)
    k = 'key1'
    db[k] = [1,0,2,0,3,0]  
    ac = AC()
    k = 'key2'
    db[k] = ac
    db.close()

    dbtype = whichdb.whichdb(filename)
    if dbtype == 'dbhash':
        dbflag = 1 
    os.system("\\rm %s*" % filename)
    return dbflag

###################################################################
def testinstall(verbose=1):
    needed = []

    # Python version --------------------------------
    pyversion = sys.version.split()[0]
    try:
        pv = float(pyversion[0:3])
        if verbose: print "...found Python %s" % pyversion 
        if pv < RequiredPythonVersion:
            if verbose:
                print "Spire requires Python version %s or higher." % str(RequiredPythonVersion)
            needed.append('python')
    except:
        print "unable to get Python version"
        return ['python']

    # Database/Shelve test  --------------------------
    """
    dbflag = has_dbhash()
    if not dbflag:
        if verbose:
            print "shelve module does not have dbhash database."
            print "Ask your Unix adminstrator to install db4x-devel-4.3.xx.rpm"
        needed.append('db')
    else:
        if verbose:
            print "...found dbhash database for Python shelve module."
    """

    # check for PIL, Numeric, and Gnuplot packages
    try:
        import Image
        if verbose: print "...found Image.%s" % Image.VERSION
    except:
        needed.append('image')
    try:
        import numpy
        if verbose: print "...found numpy.%s" % numpy.__version__
    except:
        needed.append('numpy')
    try:
        import Gnuplot
        if verbose: print "...found Gnuplot.%s" % Gnuplot.__version__
    except:
        needed.append('gnuplot')
        
    # Pmw ---------------------------------------------
    pmwflag = 0
    try:
        import Pmw
        if verbose: print "...found Pmw.%s" % Pmw._version
        pmwflag = 1
    except:
        needed.append('pmw')

    # Tkinter ------------------------------------------
    tkflag = 0
    try:
        import Tkinter
        if verbose: print "...found Tkinter with Tk %s" % str(Tkinter.TkVersion)
        if Tkinter.TkVersion < RequiredTkVersion:
            print "Spire requires Tcl/Tk version %s or higher." % str(RequiredTkVersion)
            needed.append('tk')
        else:
            tkflag = 1
    except:
        needed.append('tk')

    # Blt ---------------------------------------------
    if tkflag:
        root = Tkinter.Tk()            
        if pmwflag:
            if Pmw.Blt.haveblt(root):
                if verbose: print "...found BLT"
            else:
                needed.append('blt')
        else:
            needed.append('blt')
        root.destroy()
    else:
        needed += ['blt']

    if len(needed) > 0:
        return needed
    else:
        return "ok"

#########################################################################
if __name__ == '__main__':

    res = testinstall(verbose=1)
    if res == 'ok':
        print "All components required by Spire can be located."
    else:
        print "The following components were not found:"
        for item in res:
            print item

