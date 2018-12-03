#!/usr/bin/env python
import                os, sys
from   Tkinter import *

root = Tk()  # There are problems if this is not in the __main__ module
import Spire.GB as GB
import Spire.GG as GG
import Spire.spiderMain as spiderMain

######################################################################

GG.topwindow = Toplevel(root)

root.withdraw()

args    = sys.argv[1:]
project = 0
small   = 0    # Small font (not yet implemented)
if len(args) > 0:
    if '-s' in args:
        small = 1
        i     = args.index('-s')
        del args[i]
    if len(args) > 0: project = args[0]
    
# Create instance of GUI 
GG.spidui = spiderMain.SpiderGUI(GG.topwindow, project, usesmallfont=small)

if project:
    GG.spidui.openProject(project)
    #GG.spidui.readUserPrefs()

root.mainloop()
